#!/bin/bash -e
#SBATCH --job-name=s15_evaladmix
#SBATCH --account=uow04177
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --output=logs/%x_%A_%a.out
# evalAdmix on the best ngsadmix run per k (from tables/S8_ngsadmix_runs.tsv), same thinned beagle
# minMaf and misTol match ngsadmix. residual correlations near 0 = the k fits, + within a group = structure missed
# array index = line in sets.tsv (5 = all_on_hilli, 6 = all_on_stygia)
cd "${SLURM_SUBMIT_DIR}"
EA=/nesi/project/uow04177/${USER}/tools/evalAdmix/evalAdmix
read TAG REF LIST ZARR <<< "$(sed -n "${SLURM_ARRAY_TASK_ID}p" sets.tsv)"
mkdir -p evaladmix
for K in 2 3 4 5; do
  SEED=$(awk -v t=${TAG} -v k=${K} '$1 == t && $2 == k {print $6}' tables/S8_ngsadmix_runs.tsv)
  O=evaladmix/${TAG}.K${K}.corres.txt
  [ -s ${O} ] && continue
  ${EA} -beagle admix/${TAG}.thin5.beagle.gz -fname admix/${TAG}.K${K}.s${SEED}.fopt.gz \
    -qname admix/${TAG}.K${K}.s${SEED}.qopt -minMaf 0.05 -misTol 0.05 -P ${SLURM_CPUS_PER_TASK} -o ${O} \
    > evaladmix/${TAG}.K${K}.log 2>&1
  echo "K${K} (seed ${SEED}) done"
done
