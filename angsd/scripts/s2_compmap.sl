#!/bin/bash -e
#SBATCH --job-name=s2_compmap
#SBATCH --account=uow04177
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=36G
#SBATCH --output=logs/%x_%A_%a.out
# competitive mapping, each pair aligned once against all four assemblies combined
# mapq >= 20 = placed on one species, < 20 = conserved / ambiguous
# split of placed reads checks for mislabelled species, hybrids, contamination
# unmapped pairs kept for the stage 3 screen
cd "${SLURM_SUBMIT_DIR}"
module purge
module load Miniforge3
source $(conda info --base)/etc/profile.d/conda.sh
export PYTHONNOUSERSITE=1
conda activate /nesi/project/uow04177/${USER}/conda_envs/ch2_zarr
set -o pipefail

S=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)
BAM=compmap/bam/${S}.bam
QC=compmap/qc
T=${SLURM_CPUS_PER_TASK}
echo "sample: ${S} | reference: ref/compref.fa (4 assemblies)"
bwa 2>&1 | grep Version || true
if [ -f "compmap/bam/${S}.done" ]; then echo "already done, skipping"; exit 0; fi
mkdir -p tmp/compmap/${S}

bwa mem -t ${T} -R "@RG\tID:${S}\tSM:${S}\tLB:${S}\tPL:ILLUMINA\tPU:HVYFMDSX5.3" \
    ref/compref.fa trimmed/${S}_R1.fq.gz trimmed/${S}_R2.fq.gz \
  | samtools fixmate -m -u -@ ${T} - - \
  | samtools sort -u -@ ${T} -T tmp/compmap/${S}/sort - \
  | samtools markdup -@ ${T} - ${BAM}
samtools index -@ ${T} ${BAM}
samtools flagstat -@ ${T} ${BAM} > ${QC}/${S}.flagstat.txt
samtools coverage ${BAM} > ${QC}/${S}.coverage.tsv

# primary, non-dup reads: placed (mapq >= 20) per assembly, ambiguous (mapq < 20), unmapped
samtools view -@ 4 -F 3844 ${BAM} \
  | awk -v s=${S} 'BEGIN{OFS="\t"} {split($3,a,"__"); if ($5>=20) q[a[1]]++; else low++}
      END{for (k in q) print s,"mapq20",k,q[k]; print s,"mapq_lt20","any",low+0}' > ${QC}/${S}.assign.tsv
printf "%s\tunmapped\tnone\t%s\n" ${S} $(samtools view -c -f 4 -F 0x900 ${BAM}) >> ${QC}/${S}.assign.tsv
cat ${QC}/${S}.assign.tsv

# pairs where neither read maps anywhere, for the contam screen
samtools view -u -f 12 -F 0x900 ${BAM} \
  | samtools collate -u -O - tmp/compmap/${S}/coll \
  | samtools fastq -n -0 /dev/null -s /dev/null \
      -1 compmap/unmapped/${S}_R1.fq.gz -2 compmap/unmapped/${S}_R2.fq.gz -
rm -rf tmp/compmap/${S}
touch compmap/bam/${S}.done
echo "finished ${S}"
