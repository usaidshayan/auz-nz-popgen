#!/bin/bash -e
#SBATCH --job-name=f04b_pca_lowhet
#SBATCH --account=uow04177
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/%x_%A_%a.out
# pca without the 7 hilli flies with excess het-like sites (lists/hilli_highhet.txt), same snps and pcangsd command as f04
# array index = line in lists/sets.tsv (1, 2)
cd "${SLURM_SUBMIT_DIR}"
set -o pipefail
GL=/nesi/project/uow04177/${USER}/conda_envs/ch2_gl
BGZIP=/nesi/project/uow04177/${USER}/conda_envs/ch2_zarr/bin/bgzip
read TAG REF LIST <<< "$(sed -n "${SLURM_ARRAY_TASK_ID}p" lists/sets.tsv)"
O=pca/${TAG}/lowhet
COLS=$(awk 'NR == FNR {k[$1]; next} ($1 in k) {printf ",%d-%d", 3 * FNR + 1, 3 * FNR + 3}' lists/${TAG}_lowhet.txt ${LIST})
zcat gl/${TAG}/${TAG}.beagle.gz | cut -f1-3${COLS} | ${BGZIP} -@ 8 > ${O}.beagle.gz
echo "${TAG}: $(( ( $(zcat ${O}.beagle.gz | head -1 | wc -w) - 3 ) / 3 )) flies kept"
${GL}/bin/pcangsd -b ${O}.beagle.gz -o ${O} -t ${SLURM_CPUS_PER_TASK} --maf 0.05 --admix --inbreed-samples > ${O}.pcangsd.out 2>&1
${GL}/bin/python scripts/f04_cov2pcs.py ${O} lists/${TAG}_lowhet.txt | head -1
rm -f ${O}.beagle.gz
