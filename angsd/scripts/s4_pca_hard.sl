#!/bin/bash -e
#SBATCH --job-name=s4_pca_hard
#SBATCH --account=uow04177
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --output=logs/%x_%A_%a.out
cd "${SLURM_SUBMIT_DIR}"
read TAG REF LIST ZARR <<< "$(sed -n "${SLURM_ARRAY_TASK_ID}p" sets.tsv)"
echo "set ${TAG}: samples ${LIST}, hard calls from zarr_pvp/${ZARR}.zarr"
/nesi/project/uow04177/${USER}/conda_envs/ch2_zarr/bin/python scripts/s4_pca_hard.py zarr_pvp/${ZARR}.zarr ${LIST} pca_hard/${TAG}
