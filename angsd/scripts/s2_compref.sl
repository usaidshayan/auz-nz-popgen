#!/bin/bash -e
#SBATCH --job-name=s2_compref
#SBATCH --account=uow04177
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --output=logs/%x_%j.out
# join the four masked calliphora assemblies into one ref, contigs prefixed with species (hilli__, stygia__ ...)
cd "${SLURM_SUBMIT_DIR}"
module purge
module load Miniforge3
source $(conda info --base)/etc/profile.d/conda.sh
export PYTHONNOUSERSITE=1
conda activate /nesi/project/uow04177/${USER}/conda_envs/ch2_zarr
set -o pipefail
if [ -f ref/compref.done ]; then echo "already done, skipping"; exit 0; fi
rm -f ref/compref.fa
n=0
for FA in ../../genome_assemblies/0[1-4]_*_masked.fa; do
  SP=$(basename ${FA} .fa | sed -E 's/^[0-9]+_//; s/_masked$//')
  echo "adding ${FA} as ${SP}__"
  sed "s/^>\([^[:space:]]*\).*/>${SP}__\1/" ${FA} >> ref/compref.fa
  n=$((n+1))
done
[ ${n} = 4 ] || { echo "expected 4 assemblies, found ${n}"; exit 1; }
samtools faidx ref/compref.fa
awk '{split($1,a,"__"); bp[a[1]]+=$2; nc[a[1]]++} END{for (k in bp) printf "%s: %d contigs, %.1f Mb\n", k, nc[k], bp[k]/1e6}' ref/compref.fa.fai
bwa index ref/compref.fa
touch ref/compref.done
echo "finished compref index"
