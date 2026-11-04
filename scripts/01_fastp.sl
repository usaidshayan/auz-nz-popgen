#!/bin/bash -e
#SBATCH --job-name=fastp
#SBATCH --account=uow04177
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --output=logs/fastp_%A_%a.out

cd "${SLURM_SUBMIT_DIR}"
module purge
module load fastp/0.23.4-GCC-11.3.0

S=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)
SUFFIX=DKDL230000819-1A_HVYFMDSX5_L3
echo "sample: ${S}"

if [ -f "trimmed/${S}.done" ]; then echo "already done, skipping"; exit 0; fi

fastp \
  --in1 raw_data/${S}_${SUFFIX}_1.fq.gz --in2 raw_data/${S}_${SUFFIX}_2.fq.gz \
  --out1 trimmed/${S}_R1.fq.gz --out2 trimmed/${S}_R2.fq.gz \
  --detect_adapter_for_pe --trim_poly_g --length_required 50 \
  --thread ${SLURM_CPUS_PER_TASK} \
  --json qc/fastp/${S}.fastp.json --html qc/fastp/${S}.fastp.html \
  --report_title "${S}"

touch "trimmed/${S}.done"
