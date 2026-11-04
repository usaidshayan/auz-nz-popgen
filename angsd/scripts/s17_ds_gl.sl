#!/bin/bash -e
#SBATCH --job-name=s17_ds_gl
#SBATCH --account=uow04177
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%A_%a.out
# gls for hilli_on_hilli with one sample swapped for its downsampled bam, at the original snp sites only
# same read and depth filters as s6, major/minor taken from the original sites
cd "${SLURM_SUBMIT_DIR}"
module purge
module load angsd/0.935-GCC-9.2.0
read SET CHUNK <<< "$(sed -n "${SLURM_ARRAY_TASK_ID}p" robust/gl/tasks.tsv)"
source gl/hilli_on_hilli/depth.env
mkdir -p robust/gl/${SET}
OUT=robust/gl/${SET}/${CHUNK}
if [ -f ${OUT}.done ]; then echo "already done, skipping"; exit 0; fi
angsd -b robust/gl/${SET}.bams.txt -ref ref/01_hilli_masked.fa -rf gl/hilli_on_hilli/chunks/${CHUNK}.rf \
  -sites robust/gl/hilli_on_hilli.sites -out ${OUT} \
  -GL 2 -doGlf 2 -doMajorMinor 3 -doMaf 1 \
  -minMapQ 20 -minQ 20 -remove_bads 1 -uniqueOnly 1 -baq 1 -C 50 \
  -doCounts 1 -setMinDepth ${MIN_DEPTH} -setMaxDepth ${MAX_DEPTH} -minInd ${MIN_IND} -nThreads 4
touch ${OUT}.done
