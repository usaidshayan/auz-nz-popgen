#!/bin/bash -e
#SBATCH --job-name=s6_gl
#SBATCH --account=uow04177
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%A_%a.out
# genotype likelihoods for one 10 mb chunk of one set
# gatk gl model (-GL 2), major/minor from likelihoods, snp if lrt p < 1e-6
# reads: mapq >= 20, bq >= 20, -remove_bads, unique only, baq + -C 50
# out: beagle gls (for pcangsd / ngsadmix / ngsrelate) + mafs
cd "${SLURM_SUBMIT_DIR}"
module purge
module load angsd/0.935-GCC-9.2.0
read TAG REF CHUNK <<< "$(sed -n "${SLURM_ARRAY_TASK_ID}p" gl/tasks.tsv)"
source gl/${TAG}/depth.env
OUT=gl/${TAG}/chunks/${CHUNK}
echo "set ${TAG} | ${CHUNK} | depth ${MIN_DEPTH}-${MAX_DEPTH} | minInd ${MIN_IND}"
if [ -f ${OUT}.done ]; then echo "already done, skipping"; exit 0; fi
awk '{print $1":"$2"-"$3}' regions_pvp/${REF}/${CHUNK}.tsv > ${OUT}.rf
START=$(date +%s)
angsd -b gl/${TAG}/bams.txt -ref ref/${REF}.fa -rf ${OUT}.rf -out ${OUT} \
  -GL 2 -doGlf 2 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 \
  -minMapQ 20 -minQ 20 -remove_bads 1 -uniqueOnly 1 -baq 1 -C 50 \
  -doCounts 1 -setMinDepth ${MIN_DEPTH} -setMaxDepth ${MAX_DEPTH} -minInd ${MIN_IND} -nThreads 4
touch ${OUT}.done
echo "finished ${CHUNK}: $(zcat ${OUT}.mafs.gz | tail -n +2 | wc -l) SNPs in $(( $(date +%s) - START )) s"
