#!/bin/bash -e
#SBATCH --job-name=s6_concat
#SBATCH --account=uow04177
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%A_%a.out
# join chunk files of one set into <TAG>.beagle.gz and <TAG>.mafs.gz, one header each
cd "${SLURM_SUBMIT_DIR}"
set -o pipefail
BGZIP=/nesi/project/uow04177/${USER}/conda_envs/ch2_zarr/bin/bgzip
read TAG REF LIST ZARR <<< "$(sed -n "${SLURM_ARRAY_TASK_ID}p" sets.tsv)"
OUT=gl/${TAG}/${TAG}
if [ -f ${OUT}.done ]; then echo "already done, skipping"; exit 0; fi
EXPECTED=$(awk -v t=${TAG} '$1==t' gl/tasks.tsv | wc -l)
DONE=$(ls gl/${TAG}/chunks/*.done | wc -l)
echo "${TAG}: ${DONE} of ${EXPECTED} chunks finished"
[ ${DONE} = ${EXPECTED} ] || { echo "CHUNKS MISSING"; exit 1; }
for EXT in beagle mafs; do
  FIRST=1
  for f in $(ls gl/${TAG}/chunks/chunk_*.${EXT}.gz | sort); do
    if [ ${FIRST} = 1 ]; then zcat ${f}; FIRST=0; else zcat ${f} | tail -n +2; fi
  done | ${BGZIP} -@ 8 > ${OUT}.${EXT}.gz.tmp
  mv ${OUT}.${EXT}.gz.tmp ${OUT}.${EXT}.gz
done
NB=$(zcat ${OUT}.beagle.gz | tail -n +2 | wc -l)
NM=$(zcat ${OUT}.mafs.gz | tail -n +2 | wc -l)
echo "${TAG}: ${NB} SNPs in beagle, ${NM} in mafs"
[ ${NB} = ${NM} ] || { echo "LINE COUNT MISMATCH"; exit 1; }
touch ${OUT}.done
