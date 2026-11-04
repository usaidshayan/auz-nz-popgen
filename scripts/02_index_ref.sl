#!/bin/bash -e
#SBATCH --account=uow04177
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%j.out

cd "${SLURM_SUBMIT_DIR}"
: "${REF_NAME:?REF_NAME not set, submit with --export=ALL,REF_NAME=...}"
module purge
module load Miniforge3
source $(conda info --base)/etc/profile.d/conda.sh
export PYTHONNOUSERSITE=1
conda activate /nesi/project/uow04177/${USER}/conda_envs/ch2_zarr

bwa index ref/${REF_NAME}.fa
samtools faidx ref/${REF_NAME}.fa
ls -l ref/${REF_NAME}.fa*
