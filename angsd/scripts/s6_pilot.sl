#!/bin/bash -e
#SBATCH --job-name=s6_pilot
#SBATCH --account=uow04177
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%A_%a.out
# depth pilot on the 2 mb test region, same read filters as the main run
# depth window = 1st-99th percentile of total depth (drops collapsed repeats / paralogs)
# site needs reads in at least half the samples (-minInd)
cd "${SLURM_SUBMIT_DIR}"
module purge
module load angsd/0.935-GCC-9.2.0
read TAG REF LIST ZARR <<< "$(sed -n "${SLURM_ARRAY_TASK_ID}p" sets.tsv)"
mkdir -p gl/${TAG}/chunks
sed "s#^#bam_pvp/${REF}/#; s#\$#.bam#" ${LIST} > gl/${TAG}/bams.txt
NIND=$(wc -l < ${LIST})
awk '{print $1":"$2"-"$3}' regions_pvp/${REF}/test.tsv > gl/${TAG}/pilot.rf
angsd -b gl/${TAG}/bams.txt -ref ref/${REF}.fa -rf gl/${TAG}/pilot.rf -out gl/${TAG}/pilot \
  -minMapQ 20 -minQ 20 -remove_bads 1 -uniqueOnly 1 -baq 1 -C 50 \
  -doCounts 1 -doDepth 1 -maxDepth 1000 -nThreads 4
awk -v n=${NIND} '{for (i = 2; i <= NF; i++) {c[i-1] = $i; tot += $i}
    for (d = 1; d < NF; d++) {cum += c[d]; if (lo == "" && cum >= 0.01 * tot) lo = d; if (hi == "" && cum >= 0.99 * tot) hi = d}
    printf "MIN_DEPTH=%d\nMAX_DEPTH=%d\nMIN_IND=%d\n", lo, hi, int((n + 1) / 2)}' gl/${TAG}/pilot.depthGlobal > gl/${TAG}/depth.env
echo "${TAG}: ${NIND} samples"; cat gl/${TAG}/depth.env
