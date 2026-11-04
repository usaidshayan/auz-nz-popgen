#!/bin/bash -e
#SBATCH --job-name=s10_het
#SBATCH --account=uow04177
#SBATCH --time=16:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --output=logs/%x_%A_%a.out
# genome-wide het for one sample on one assembly, no genotype calls
# angsd saf -> realSFS single-sample sfs, het = middle class / all sites
# sites need >= 2 reads and <= 3x the sample's mean depth (drops repeats)
# realSFS in 20 mb blocks, summed for the genome value
cd "${SLURM_SUBMIT_DIR}"
module purge
module load angsd/0.935-GCC-9.2.0
read S REF <<< "$(sed -n "${SLURM_ARRAY_TASK_ID}p" het/tasks.tsv)"
OUT=het/${REF}/${S}
if [ -f ${OUT}.done ]; then echo "already done, skipping"; exit 0; fi
MEAN=$(awk -F'\t' -v s=${S} '$1==s {print $6}' ../qc/bam_summary_${REF}.tsv)
MAXD=$(awk -v m=${MEAN} 'BEGIN {d = int(3 * m + 1); if (d < 5) d = 5; print d}')
echo "sample ${S} | ${REF} | mean depth ${MEAN} | site depth 2-${MAXD}"
angsd -i bam_pvp/${REF}/${S}.bam -ref ref/${REF}.fa -anc ref/${REF}.fa -out ${OUT} \
  -dosaf 1 -GL 2 -minMapQ 20 -minQ 20 -remove_bads 1 -uniqueOnly 1 -baq 1 -C 50 \
  -doCounts 1 -setMinDepth 2 -setMaxDepth ${MAXD} -nThreads 4
realSFS ${OUT}.saf.idx -fold 1 -P 4 -nSites 20000000 > ${OUT}.ml
awk -v s=${S} '{h += $2; for (i = 1; i <= NF; i++) t += $i; n++}
  END {printf "%s: %d blocks, %.0f sites, heterozygosity %.6f\n", s, n, t, h / t}' ${OUT}.ml
touch ${OUT}.done
