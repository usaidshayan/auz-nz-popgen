#!/bin/bash -e
#SBATCH --account=uow04177
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=logs/%x_%A_%a.out

cd "${SLURM_SUBMIT_DIR}"
: "${REF_NAME:?REF_NAME not set, submit with --export=ALL,REF_NAME=...}"
module purge
module load Miniforge3
source $(conda info --base)/etc/profile.d/conda.sh
export PYTHONNOUSERSITE=1
conda activate /nesi/project/uow04177/${USER}/conda_envs/ch2_zarr
set -o pipefail

S=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)
REF=ref/${REF_NAME}.fa
OUT=bam/${REF_NAME}
QC=qc/bam/${REF_NAME}
T=${SLURM_CPUS_PER_TASK}
echo "sample: ${S} | reference: ${REF_NAME}"
bwa 2>&1 | grep Version || true
samtools --version | head -1 || true

if [ -f "${OUT}/${S}.done" ]; then echo "already done, skipping"; exit 0; fi

mkdir -p ${OUT} ${QC} tmp/${REF_NAME}/${S}
bwa mem -t ${T} -R "@RG\tID:${S}\tSM:${S}\tLB:${S}\tPL:ILLUMINA\tPU:HVYFMDSX5.3" \
    ${REF} trimmed/${S}_R1.fq.gz trimmed/${S}_R2.fq.gz \
  | samtools fixmate -m -u -@ ${T} - - \
  | samtools sort -u -@ ${T} -T tmp/${REF_NAME}/${S}/sort - \
  | samtools markdup -@ ${T} - ${OUT}/${S}.bam

samtools index -@ ${T} ${OUT}/${S}.bam
samtools quickcheck -v ${OUT}/${S}.bam
samtools flagstat -@ ${T} ${OUT}/${S}.bam > ${QC}/${S}.flagstat.txt
samtools coverage ${OUT}/${S}.bam > ${QC}/${S}.coverage.tsv
rm -rf tmp/${REF_NAME}/${S}
touch ${OUT}/${S}.done
echo "finished ${S}"
