#!/bin/bash -e
#SBATCH --job-name=f01_pilot
#SBATCH --account=uow04177
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%A_%a.out
# depth pilot on 5 random 10 mb chunks (seed 1), cleaned bams, same read filters as the gl run
# max depth = 99th percentile of total depth, light touch, only cuts collapsed repeats
# site needs reads in at least half the samples (-minInd). array index = line in lists/sets.tsv
cd "${SLURM_SUBMIT_DIR}"
module purge
module load angsd/0.935-GCC-9.2.0
read TAG REF LIST <<< "$(sed -n "${SLURM_ARRAY_TASK_ID}p" lists/sets.tsv)"
mkdir -p gl/${TAG}/chunks
sed "s#^#bam/${REF}/#; s#\$#.bam#" ${LIST} > gl/${TAG}/bams.txt
NIND=$(wc -l < ${LIST})
ls regions/${REF}/chunk_*.tsv | awk 'BEGIN {srand(1)} {print rand() "\t" $0}' | sort -k1,1g | head -5 | cut -f2 \
  | xargs cat | awk '{print $1":"$2"-"$3}' > gl/${TAG}/pilot.rf
echo "${TAG}: ${NIND} samples, pilot on $(awk -F'[:-]' '{s += $3 - $2 + 1} END {print s}' gl/${TAG}/pilot.rf) bp"
angsd -b gl/${TAG}/bams.txt -ref ref/${REF}.fa -rf gl/${TAG}/pilot.rf -out gl/${TAG}/pilot \
  -minMapQ 20 -minQ 20 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -baq 1 -C 50 \
  -doCounts 1 -doDepth 1 -maxDepth 2000 -nThreads 4
awk -v n=${NIND} '{for (i = 2; i <= NF; i++) {c[i-1] = $i; tot += $i}
    for (d = 1; d < NF; d++) {cum += c[d]; if (hi == "" && cum >= 0.99 * tot) hi = d}
    printf "MAX_DEPTH=%d\nMIN_IND=%d\n", hi, int((n + 1) / 2)}' gl/${TAG}/pilot.depthGlobal > gl/${TAG}/depth.env
cat gl/${TAG}/depth.env
