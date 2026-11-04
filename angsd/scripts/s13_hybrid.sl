#!/bin/bash -e
#SBATCH --job-name=s13_hybrid
#SBATCH --account=uow04177
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=logs/%x_%A_%a.out
# hybrid index + interspecific het from the stage 6 beagle, both species together
# array index = line in sets.tsv (5 = all_on_hilli, 6 = all_on_stygia)
cd "${SLURM_SUBMIT_DIR}"
GL=/nesi/project/uow04177/${USER}/conda_envs/ch2_gl
read TAG REF LIST ZARR <<< "$(sed -n "${SLURM_ARRAY_TASK_ID}p" sets.tsv)"
mkdir -p hybrid
${GL}/bin/python scripts/s13_hybrid_index.py gl/${TAG}/${TAG}.beagle.gz ${LIST} hybrid/${TAG}
