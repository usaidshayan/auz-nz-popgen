#!/bin/bash -e
#SBATCH --job-name=f02_gl
#SBATCH --account=uow04177
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%A_%a.out
# genotype likelihoods for one 10 mb chunk of one set (gl/tasks.tsv)
# gatk gl model, reference allele as major (-doMajorMinor 4), snp if lrt p < 1e-6, no maf filter here
# reads: mapq >= 20, bq >= 20, proper pairs, baq + -C 50. sites: <= 99th pct depth, reads in half the flies
# out: beagle (pcangsd / ngsadmix / ngsrelate / ngsld), mafs, bcf with gl + posterior calls (kestrel, g123)
cd "${SLURM_SUBMIT_DIR}"
module purge
module load angsd/0.935-GCC-9.2.0
read TAG REF CHUNK <<< "$(sed -n "${SLURM_ARRAY_TASK_ID}p" gl/tasks.tsv)"
source gl/${TAG}/depth.env
OUT=gl/${TAG}/chunks/${CHUNK}
echo "set ${TAG} | ${CHUNK} | max depth ${MAX_DEPTH} | minInd ${MIN_IND}"
if [ -f ${OUT}.done ]; then echo "already done, skipping"; exit 0; fi
awk '{print $1":"$2"-"$3}' regions/${REF}/${CHUNK}.tsv > ${OUT}.rf
angsd -b gl/${TAG}/bams.txt -ref ref/${REF}.fa -rf ${OUT}.rf -out ${OUT} \
  -GL 2 -doGlf 2 -doMajorMinor 4 -doMaf 1 -SNP_pval 1e-6 -doPost 1 -doGeno 8 -dobcf 1 --ignore-RG 0 \
  -minMapQ 20 -minQ 20 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -baq 1 -C 50 \
  -doCounts 1 -setMaxDepth ${MAX_DEPTH} -minInd ${MIN_IND} -nThreads 4
rm -f ${OUT}.geno.gz
touch ${OUT}.done
echo "finished ${CHUNK}: $(zcat ${OUT}.mafs.gz | tail -n +2 | wc -l) snps"
