#!/bin/bash
# live progress view, ctrl+c to exit. usage: bash monitor.sh
cd "$(dirname "$(readlink -f "$0")")" || exit 1
while true; do
  clear
  date
  echo "== queue"
  squeue --me
  echo
  echo "== reference indexes finished"
  ls ref/*.sa 2>/dev/null || echo "(none yet)"
  echo
  echo "== alignments"
  for R in 01_hilli_masked 03_stygia_masked; do
    n=$(ls bam/$R/*.done 2>/dev/null | wc -l)
    echo "$R: $n samples finished"
    for f in $(find logs -name "align_${R}_*.out" -mmin -10 2>/dev/null | sort); do
      grep -q '^finished' "$f" && continue
      S=$(grep -m1 '^sample:' "$f" | awk '{print $2}')
      tot=$(awk -F'\t' -v s="$S" '$1==s {print $3}' qc/fastp_summary.tsv)
      grep Processed "$f" | awk -v s="$S" -v tot="$tot" '{r+=$3; t+=$9} END { if (r>0 && tot>0) printf "   %-16s %5.1f%% of reads aligned (%d of %d), rough ETA for this phase: %.0f min\n", s, 100*r/tot, r, tot, (tot-r)/(r/t)/60; else printf "   %-16s starting\n", s }'
    done
  done
  sleep 30
done
