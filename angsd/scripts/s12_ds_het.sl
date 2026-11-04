#!/bin/bash -e
#SBATCH --job-name=s12_ds_het
#SBATCH --account=uow04177
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --output=logs/%x_%A_%a.out
# downsample one bam (pairs kept together) to a target depth, rerun the stage 10 het estimate
# downsampled bam deleted after
cd "${SLURM_SUBMIT_DIR}"
SAMTOOLS=/nesi/project/uow04177/${USER}/conda_envs/ch2_zarr/bin/samtools
module purge
module load angsd/0.935-GCC-9.2.0
read S REF TARGET <<< "$(sed -n "${SLURM_ARRAY_TASK_ID}p" het_ds/tasks.tsv)"
OUT=het_ds/${S}__${TARGET}x
if [ -f ${OUT}.done ]; then echo "already done, skipping"; exit 0; fi
MEAN=$(awk -F'\t' -v s=${S} '$1==s {print $6}' ../qc/bam_summary_${REF}.tsv)
FRAC=$(awk -v t=${TARGET} -v m=${MEAN} 'BEGIN {printf "%.4f", t / m}')
MAXD=$(awk -v t=${TARGET} 'BEGIN {d = int(3 * t + 1); if (d < 5) d = 5; print d}')
echo "sample ${S} | ${REF} | depth ${MEAN} -> ${TARGET} (keep ${FRAC} of pairs) | site depth 2-${MAXD}"
${SAMTOOLS} view -b -@ 4 --subsample ${FRAC} --subsample-seed 11 -o ${OUT}.bam bam_pvp/${REF}/${S}.bam
${SAMTOOLS} index ${OUT}.bam
angsd -i ${OUT}.bam -ref ref/${REF}.fa -anc ref/${REF}.fa -out ${OUT} \
  -dosaf 1 -GL 2 -minMapQ 20 -minQ 20 -remove_bads 1 -uniqueOnly 1 -baq 1 -C 50 \
  -doCounts 1 -setMinDepth 2 -setMaxDepth ${MAXD} -nThreads 4
realSFS ${OUT}.saf.idx -fold 1 -P 4 -nSites 20000000 > ${OUT}.ml
awk -v s=${S} -v t=${TARGET} '{h += $2; for (i = 1; i <= NF; i++) x += $i}
  END {printf "%s at %sx: %.0f sites, heterozygosity %.6f\n", s, t, x, h / x}' ${OUT}.ml
rm -f ${OUT}.bam ${OUT}.bam.bai
touch ${OUT}.done
