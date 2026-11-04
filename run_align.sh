#!/bin/bash
# usage: bash run_align.sh test   (index references if needed, then align WE_CH_141 to each reference)
#        bash run_align.sh all    (align every sample in samples.txt to each reference)
cd "$(dirname "$(readlink -f "$0")")" || exit 1
MODE=$1
REFS="01_hilli_masked 03_stygia_masked"

if [ "$MODE" = "test" ]; then
  n=$(grep -n -x WE_CH_141 samples.txt | cut -d: -f1)
  for R in $REFS; do
    mkdir -p ref bam/$R qc/bam/$R tmp/$R logs
    ln -sf "$(readlink -f ../genome_assemblies/$R.fa)" ref/$R.fa
    DEP=""
    if [ ! -s ref/$R.fa.sa ] || [ ! -s ref/$R.fa.fai ]; then
      IDX=$(sbatch --parsable --job-name=index_$R --export=ALL,REF_NAME=$R scripts/02_index_ref.sl)
      DEP="--dependency=afterok:$IDX"
      echo "index job for $R: $IDX"
    else
      echo "index for $R already present"
    fi
    sbatch --job-name=align_$R $DEP --array=$n --export=ALL,REF_NAME=$R scripts/03_align.sl
  done
elif [ "$MODE" = "all" ]; then
  for R in $REFS; do
    sbatch --job-name=align_$R --array=1-$(wc -l < samples.txt)%4 --export=ALL,REF_NAME=$R scripts/03_align.sl
  done
else
  echo "usage: bash run_align.sh test|all"
fi
