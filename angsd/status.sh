#!/bin/bash
# one-shot progress report for the ch2 angsd work. usage: bash angsd/status.sh
cd ~/00_nesi_projects/uow04177_nobackup/ch2_pre-vs-post_wgr/angsd
date
echo "== queue"
squeue --me -h -o "%T" | sort | uniq -c
squeue --me -o "%.20i %.14j %.9T %.10M %R" | head -12
echo "== stage progress"
printf "S2 index:             %s\n" "$([ -f ref/compref.done ] && echo done || echo running)"
printf "S2 competitive maps:  %s of 36\n" "$(ls compmap/bam/*.done 2>/dev/null | wc -l)"
printf "S2 summary table:     %s\n" "$([ -s tables/S2_competitive_mapping.tsv ] && echo done || echo not yet)"
printf "S4 hard-call PCAs:    %s of 6\n" "$(ls pca_hard/*.pcs.tsv 2>/dev/null | wc -l)"
printf "S4 figure:            %s\n" "$([ -s tables/S4_pca_hardcalls.png ] && echo done || echo not yet)"
printf "S6 depth pilots:      %s of 6\n" "$(ls gl/*/depth.env 2>/dev/null | wc -l)"
printf "S6 GL chunks:         %s of %s\n" "$(ls gl/*/chunks/*.done 2>/dev/null | wc -l)" "$(cat gl/tasks.tsv 2>/dev/null | wc -l)"
printf "S6 joined sets:       %s of 6\n" "$(ls gl/*/*.done 2>/dev/null | wc -l)"
printf "S10 heterozygosity:   %s of 72\n" "$(ls het/*/*.done 2>/dev/null | wc -l)"
printf "S10 summary table:    %s\n" "$([ -s tables/S10_heterozygosity.tsv ] && echo done || echo not yet)"
echo "== jobs that failed, timed out or ran out of memory (since 30 Oct)"
sacct -X -n -S 2026-10-30 --format=JobID%20,JobName%14,State%14,Elapsed,ExitCode \
  | grep -E "FAILED|TIMEOUT|OUT_OF_ME|NODE_FAIL" || echo "none"
