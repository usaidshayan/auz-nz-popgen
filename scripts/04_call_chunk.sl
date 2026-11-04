#!/bin/bash -e
#SBATCH --account=uow04177
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%A_%a.out

cd "${SLURM_SUBMIT_DIR}"
: "${REF_NAME:?REF_NAME not set}"
: "${SAMPLE_SET:?SAMPLE_SET not set}"
: "${MODE:?MODE not set (test or full)}"
module purge
module load Miniforge3
source $(conda info --base)/etc/profile.d/conda.sh
export PYTHONNOUSERSITE=1
conda activate /nesi/project/uow04177/${USER}/conda_envs/ch2_zarr
set -o pipefail

REFTAG=$(echo ${REF_NAME} | sed -E 's/^[0-9]+_//; s/_masked$//')
TAG=${SAMPLE_SET}_on_${REFTAG}
OUT=vcf/${TAG}
if [ "${MODE}" = "test" ]; then
  REGIONS=regions/${REF_NAME}/test.tsv
  VCF=${OUT}/test/test.vcf.gz
else
  NNN=$(printf "%04d" ${SLURM_ARRAY_TASK_ID})
  REGIONS=regions/${REF_NAME}/chunk_${NNN}.tsv
  VCF=${OUT}/chunks/chunk_${NNN}.vcf.gz
fi
DONE=${VCF%.vcf.gz}.done
echo "set: ${TAG} | regions: ${REGIONS} | output: ${VCF}"
bcftools --version | head -1 || true

if [ -f "${DONE}" ]; then echo "already done, skipping"; exit 0; fi

START=$(date +%s)
# all sites with coverage are written (no variants-only flag); indel calling is off so each position has one record
# reads are skipped if unmapped, secondary, QC-fail, duplicate or supplementary (3844 = 4+256+512+1024+2048)
bcftools mpileup -f ref/${REF_NAME}.fa -b ${OUT}/bams.txt -R ${REGIONS} \
    -q 20 -Q 20 -I --ns 3844 -a FORMAT/AD,FORMAT/DP -Ou \
  | bcftools call -m -f GQ -Ou \
  | bcftools annotate -x INFO,FORMAT/PL -Oz --threads 2 -o ${VCF}.tmp
mv ${VCF}.tmp ${VCF}
bcftools index -t ${VCF}
touch ${DONE}
echo "finished ${VCF}: $(bcftools index -n ${VCF}) records in $(( $(date +%s) - START )) s"
