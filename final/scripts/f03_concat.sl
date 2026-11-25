#!/bin/bash -e
#SBATCH --job-name=f03_concat
#SBATCH --account=uow04177
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%A_%a.out
# join the chunk files of one set: beagle + mafs (one header each) and bcf. array index = line in lists/sets.tsv
cd "${SLURM_SUBMIT_DIR}"
set -o pipefail
E=/nesi/project/uow04177/${USER}/conda_envs/ch2_zarr/bin
read TAG REF LIST <<< "$(sed -n "${SLURM_ARRAY_TASK_ID}p" lists/sets.tsv)"
OUT=gl/${TAG}/${TAG}
if [ -f ${OUT}.done ]; then echo "already done, skipping"; exit 0; fi
EXPECTED=$(awk -v t=${TAG} '$1 == t' gl/tasks.tsv | wc -l)
DONE=$(ls gl/${TAG}/chunks/*.done | wc -l)
echo "${TAG}: ${DONE} of ${EXPECTED} chunks finished"
[ ${DONE} = ${EXPECTED} ] || { echo "chunks missing"; exit 1; }
for EXT in beagle mafs; do
  FIRST=1
  for f in $(ls gl/${TAG}/chunks/chunk_*.${EXT}.gz | sort); do
    if [ ${FIRST} = 1 ]; then zcat ${f}; FIRST=0; else zcat ${f} | tail -n +2; fi
  done | ${E}/bgzip -@ 8 > ${OUT}.${EXT}.gz.tmp
  mv ${OUT}.${EXT}.gz.tmp ${OUT}.${EXT}.gz
done
ls gl/${TAG}/chunks/chunk_*.bcf | sort > gl/${TAG}/bcfs.txt
${E}/bcftools concat -f gl/${TAG}/bcfs.txt --threads 8 -Ob -o ${OUT}.bcf
${E}/bcftools index ${OUT}.bcf
NB=$(zcat ${OUT}.beagle.gz | tail -n +2 | wc -l)
NM=$(zcat ${OUT}.mafs.gz | tail -n +2 | wc -l)
NV=$(${E}/bcftools index -n ${OUT}.bcf)
echo "${TAG}: ${NB} snps in beagle, ${NM} in mafs, ${NV} in bcf"
[ ${NB} = ${NM} ] && [ ${NB} = ${NV} ] || { echo "line count mismatch"; exit 1; }
touch ${OUT}.done
