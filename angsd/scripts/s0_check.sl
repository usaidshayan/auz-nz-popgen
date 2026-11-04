#!/bin/bash -e
#SBATCH --job-name=s0_check
#SBATCH --account=uow04177
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --output=logs/%x_%j.out
# integrity check from raw reads down to the stage 13 outputs before writing anything up
# part 1: md5 + line count of every raw and trimmed fastq (a broken gzip fails zcat)
# part 2: s0_check.py for refs, bams, gls and downstream tables
cd "${SLURM_SUBMIT_DIR}"
set -o pipefail
module purge
module load SAMtools/1.23.1-GCC-12.3.0
GL=/nesi/project/uow04177/${USER}/conda_envs/ch2_gl
mkdir -p check

fq() {
  f=$1
  m=$(md5sum ${f} | cut -d' ' -f1)
  if n=$(zcat ${f} | wc -l); then ok=ok; else ok=BROKEN; fi
  echo -e "$(basename ${f})\t$(dirname ${f})\t${m}\t${n}\t${ok}"
}
export -f fq
if [ ! -s check/fastq.tsv ]; then
  ls ../raw_data/*.fq.gz ../trimmed/*.fq.gz | xargs -P ${SLURM_CPUS_PER_TASK} -I{} bash -c 'fq {}' > check/fastq.tsv.tmp
  mv check/fastq.tsv.tmp check/fastq.tsv
fi
echo "fastq done: $(wc -l < check/fastq.tsv) files"

if [ ! -s check/bam_headers.tsv ]; then
  for R in 01_hilli_masked 03_stygia_masked; do
    for S in $(cat samples.txt); do
      B=bam_pvp/${R}/${S}.bam
      samtools quickcheck ${B} && q=ok || q=FAIL
      [ ${B}.bai -nt ${B} ] && i=ok || i=STALE
      echo -e "${R}\t${S}\t${q}\t${i}"
      samtools view -H ${B} > check/hdr_${R}_${S}.txt
    done
  done > check/bam_headers.tsv
fi

${GL}/bin/python scripts/s0_check.py | tee check/report.txt
