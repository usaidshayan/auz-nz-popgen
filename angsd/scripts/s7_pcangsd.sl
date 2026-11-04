#!/bin/bash -e
#SBATCH --job-name=s7_pcangsd
#SBATCH --account=uow04177
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=logs/%x_%A_%a.out
# pca from genotype likelihoods (pcangsd), low-depth samples don't get pulled to the origin
# no. of pcs picked by pcangsd
# --admix: k = pcs + 1, a hybrid would sit between the species clusters
# --inbreed-samples: per-sample inbreeding. maf < 0.05 dropped
cd "${SLURM_SUBMIT_DIR}"
GL=/nesi/project/uow04177/${USER}/conda_envs/ch2_gl
read TAG REF LIST ZARR <<< "$(sed -n "${SLURM_ARRAY_TASK_ID}p" sets.tsv)"
echo "set ${TAG}: $(zcat gl/${TAG}/${TAG}.mafs.gz | tail -n +2 | wc -l) SNPs in the Beagle file"
${GL}/bin/pcangsd -b gl/${TAG}/${TAG}.beagle.gz -o pcangsd/${TAG} -t ${SLURM_CPUS_PER_TASK} \
  --maf 0.05 --admix --inbreed-samples
${GL}/bin/python scripts/s7_cov2pcs.py pcangsd/${TAG} ${LIST}
