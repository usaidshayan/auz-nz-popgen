#!/bin/bash -e
#SBATCH --account=uow04177
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%j.out

cd "${SLURM_SUBMIT_DIR}"
: "${REF_NAME:?REF_NAME not set}"
: "${SAMPLE_SET:?SAMPLE_SET not set}"
: "${N_CHUNKS:?N_CHUNKS not set}"
module purge
module load Miniforge3
source $(conda info --base)/etc/profile.d/conda.sh
export PYTHONNOUSERSITE=1
conda activate /nesi/project/uow04177/${USER}/conda_envs/ch2_zarr
set -o pipefail

REFTAG=$(echo ${REF_NAME} | sed -E 's/^[0-9]+_//; s/_masked$//')
TAG=${SAMPLE_SET}_on_${REFTAG}
OUT=vcf/${TAG}
FINAL=${OUT}/${TAG}.allsites.vcf.gz
T=${SLURM_CPUS_PER_TASK}
if [ -f "${FINAL}.done" ]; then echo "already done, skipping"; exit 0; fi

rm -f ${OUT}/chunk_list.txt
SUM=0
for i in $(seq 1 ${N_CHUNKS}); do
  n=$(printf "%04d" ${i})
  [ -f ${OUT}/chunks/chunk_${n}.done ] || { echo "chunk ${n} is not finished"; exit 1; }
  echo ${OUT}/chunks/chunk_${n}.vcf.gz >> ${OUT}/chunk_list.txt
  SUM=$(( SUM + $(bcftools index -n ${OUT}/chunks/chunk_${n}.vcf.gz) ))
done

bcftools concat -f ${OUT}/chunk_list.txt -Oz --threads ${T} -o ${FINAL}.tmp
mv ${FINAL}.tmp ${FINAL}
bcftools index -t ${FINAL}

echo "== checks"
NREC=$(bcftools index -n ${FINAL})
echo "records in joined VCF: ${NREC} | sum over chunks: ${SUM}"
[ "${NREC}" = "${SUM}" ] || { echo "RECORD COUNT MISMATCH"; exit 1; }
bcftools query -l ${FINAL} | sort > ${OUT}/samples_in_vcf.txt
sort lists/${SAMPLE_SET}.txt > ${OUT}/samples_expected.txt
if ! diff -q ${OUT}/samples_in_vcf.txt ${OUT}/samples_expected.txt > /dev/null; then echo "SAMPLE MISMATCH"; exit 1; fi
echo "samples: $(wc -l < ${OUT}/samples_in_vcf.txt), same names as lists/${SAMPLE_SET}.txt"
CHR=$(head -1 regions/${REF_NAME}/chunk_0001.tsv | cut -f1)
bcftools view -H -r ${CHR}:1-1000000 ${FINAL} \
  | awk -v c=${CHR} '{if ($5==".") inv++; else var++} END{printf "first 1 Mb of %s: %d records without an ALT allele, %d with one\n", c, inv, var}'
bcftools view -h ${FINAL} | grep '^##FORMAT' | cut -d, -f1
touch ${FINAL}.done
echo "finished ${FINAL}"
