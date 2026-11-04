#!/bin/bash -e
#SBATCH --job-name=s8_ngsadmix
#SBATCH --account=uow04177
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --output=logs/%x_%A_%a.out
# admixture from gls (ngsadmix), k 1-5 (croft et al. 2024 used k 2-5 with snmf), 10 seeds per k
# both species together, on each assembly. array index = line in sets.tsv (5 = all_on_hilli, 6 = all_on_stygia)
# every 5th snp of the stage 6 beagle to cut ld and run time, maf >= 0.05
# a hybrid would carry both species' clusters at k = 2
cd "${SLURM_SUBMIT_DIR}"
module purge
module load angsd/0.935-GCC-9.2.0
read TAG REF LIST ZARR <<< "$(sed -n "${SLURM_ARRAY_TASK_ID}p" sets.tsv)"
mkdir -p admix
B=admix/${TAG}.thin5.beagle.gz
if [ ! -s ${B} ]; then
  zcat gl/${TAG}/${TAG}.beagle.gz | awk 'NR == 1 || (NR - 2) % 5 == 0' | gzip > ${B}.tmp
  mv ${B}.tmp ${B}
fi
echo "${TAG}: $(( $(zcat ${B} | wc -l) - 1 )) snps after thinning"
for K in 1 2 3 4 5; do
  for SEED in $(seq 1 10); do
    O=admix/${TAG}.K${K}.s${SEED}
    [ -s ${O}.qopt ] && continue
    NGSadmix -likes ${B} -K ${K} -seed ${SEED} -minMaf 0.05 -P ${SLURM_CPUS_PER_TASK} -outfiles ${O} > /dev/null 2>&1
    echo "K${K} seed ${SEED}: $(grep -o 'like=[^ ]*' ${O}.log | tail -1)"
  done
done
