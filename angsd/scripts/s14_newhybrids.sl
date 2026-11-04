#!/bin/bash -e
#SBATCH --job-name=s14_newhybrids
#SBATCH --account=uow04177
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%A_%a.out
# newhybrids on 10 replicate sets of 200 random loci (croft et al. 2024), both species together
# NLOC=500 as a second run, the lowest-depth samples have calls at only a few % of sites
# (2000 loci underflows, newhybrids returns nan for the best covered samples)
# no prior species assignment, jeffreys priors on pi and theta, 50k burn-in + 100k sweeps
# array index = line in sets.tsv (5 = all_on_hilli, 6 = all_on_stygia)
cd "${SLURM_SUBMIT_DIR}"
GL=/nesi/project/uow04177/${USER}/conda_envs/ch2_gl
NH=/nesi/project/uow04177/${USER}/tools/newhybrids/newhybrids-no-gui-linux.exe
read TAG REF LIST ZARR <<< "$(sed -n "${SLURM_ARRAY_TASK_ID}p" sets.tsv)"
NLOC=${NLOC:-200}
D=newhybrids/${TAG}_${NLOC}
[ -s ${D}/r10.txt ] || ${GL}/bin/python scripts/s14_nh_input.py admix/${TAG}.thin5.beagle.gz ${LIST} ${D} ${NLOC}
run() {
  R=$1
  mkdir -p ${D}/r${R}_run
  cd ${D}/r${R}_run
  [ -s aa-PofZ.txt ] && return 0
  ${NH} -d ../r${R}.txt --no-gui --burn-in 50000 --num-sweeps 100000 --pi-prior Jeffreys --theta-prior Jeffreys \
    -s $((10#${R} * 7)) $((10#${R} * 13)) > run.log 2>&1 || true
  [ -s aa-PofZ.txt ] && echo "replicate ${R} done" || echo "replicate ${R} FAILED"
}
export -f run
export D NH
seq -w 1 10 | xargs -P ${SLURM_CPUS_PER_TASK} -I{} bash -c 'run {}'
