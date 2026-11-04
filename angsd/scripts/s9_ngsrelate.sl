#!/bin/bash -e
#SBATCH --job-name=s9_ngsrelate
#SBATCH --account=uow04177
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --output=logs/%x_%A_%a.out
# pairwise relatedness from gls (ngsrelate v2) for every pair in one set, stage 6 beagle + angsd freqs
# checking if the mt keira samples that split off in the pcas are close kin, and for duplicates
cd "${SLURM_SUBMIT_DIR}"
set -o pipefail
GL=/nesi/project/uow04177/${USER}/conda_envs/ch2_gl
read TAG REF LIST ZARR <<< "$(sed -n "${SLURM_ARRAY_TASK_ID}p" sets.tsv)"
N=$(wc -l < ${LIST})
zcat gl/${TAG}/${TAG}.mafs.gz | awk 'NR==1 {for (i = 1; i <= NF; i++) if ($i == "knownEM") c = i; next} {print $c}' > relate/${TAG}.freq
echo "${TAG}: ${N} samples, $(wc -l < relate/${TAG}.freq) sites with frequencies"
${GL}/bin/ngsRelate -G gl/${TAG}/${TAG}.beagle.gz -n ${N} -f relate/${TAG}.freq -z ${LIST} \
  -O relate/${TAG}.res -p ${SLURM_CPUS_PER_TASK}
${GL}/bin/python scripts/s9_top.py relate/${TAG}.res
