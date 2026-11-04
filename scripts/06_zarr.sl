#!/bin/bash -e
#SBATCH --account=uow04177
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%j.out

cd "${SLURM_SUBMIT_DIR}"
: "${TAG:?TAG not set, for example hilli_on_hilli}"
module purge
module load Miniforge3
source $(conda info --base)/etc/profile.d/conda.sh
export PYTHONNOUSERSITE=1
conda activate /nesi/project/uow04177/${USER}/conda_envs/ch2_zarr
set -o pipefail

VCF=vcf/${TAG}/${TAG}.allsites.vcf.gz
SET=${TAG%%_on_*}
/nesi/project/uow04177/${USER}/conda_envs/ch2_zarr/bin/python scripts/to_zarr.py ${VCF} zarr/${TAG}.zarr
/nesi/project/uow04177/${USER}/conda_envs/ch2_zarr/bin/python scripts/zarr_check.py zarr/${TAG}.zarr ${VCF} lists/${SET}.txt
