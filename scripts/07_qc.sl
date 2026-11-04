#!/bin/bash -e
#SBATCH --account=uow04177
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%j.out

cd "${SLURM_SUBMIT_DIR}"
: "${TAG:?TAG not set, for example all_on_hilli}"
PY=${PY:-/nesi/project/uow04177/${USER}/conda_envs/ch2_zarr/bin/python}
mkdir -p qc
${PY} scripts/sample_qc.py zarr/${TAG}.zarr > qc/sample_qc_${TAG}.txt
${PY} scripts/site_qc.py zarr/${TAG}.zarr > qc/site_qc_${TAG}.txt
cat qc/sample_qc_${TAG}.txt qc/site_qc_${TAG}.txt
