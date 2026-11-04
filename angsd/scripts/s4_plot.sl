#!/bin/bash -e
#SBATCH --job-name=s4_plot
#SBATCH --account=uow04177
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/%x_%j.out
cd "${SLURM_SUBMIT_DIR}"
/nesi/project/uow04177/${USER}/conda_envs/ch2_gl/bin/python scripts/pca_plot.py pca_hard .pcs.tsv tables/S4_pca_hardcalls.png
