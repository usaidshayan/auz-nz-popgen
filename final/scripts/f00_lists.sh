#!/bin/bash
# sample lists, sets, saf groups and fst pairs for the final run
# hilli is run with and without mt keira (mka), stygia once (sample 35 kept)
cd "$(dirname "$(readlink -f "$0")")/.." || exit 1
M=samples_meta.tsv
code() { echo "$1" | cut -d_ -f1; }
for sp in hilli stygia; do
  awk -v s=${sp} '$2 == s {print $1}' ${M} | sort > lists/${sp}.txt
  for c in AUS NZ; do awk -v s=${sp} -v c=${c} '$2 == s && $4 == c {print $1}' ${M} | sort > lists/${sp}_${c}.txt; done
  for p in $(cut -d_ -f1 lists/${sp}.txt | sort -u); do grep "^${p}_" lists/${sp}.txt > lists/${sp}_${p}.txt; done
done
grep -v '^MKA_' lists/hilli.txt > lists/hilli_noMKA.txt
grep -v '^MKA_' lists/hilli_AUS.txt > lists/hilli_AUS_noMKA.txt

printf "hilli\t01_hilli_masked\tlists/hilli.txt\nhilli_noMKA\t01_hilli_masked\tlists/hilli_noMKA.txt\nstygia\t03_stygia_masked\tlists/stygia.txt\n" > lists/sets.tsv

# group, set it takes depth filters from, list
{
  for p in BL MKA SE DU PI WE; do echo -e "hilli_${p}\thilli\tlists/hilli_${p}.txt"; done
  for p in BL CA HO DU PI WE; do echo -e "stygia_${p}\tstygia\tlists/stygia_${p}.txt"; done
  echo -e "hilli_AUS\thilli\tlists/hilli_AUS.txt"
  echo -e "hilli_AUS_noMKA\thilli_noMKA\tlists/hilli_AUS_noMKA.txt"
  echo -e "hilli_NZ\thilli\tlists/hilli_NZ.txt"
  echo -e "stygia_AUS\tstygia\tlists/stygia_AUS.txt"
  echo -e "stygia_NZ\tstygia\tlists/stygia_NZ.txt"
  echo -e "hilli\thilli\tlists/hilli.txt"
  echo -e "hilli_noMKA\thilli_noMKA\tlists/hilli_noMKA.txt"
  echo -e "stygia\tstygia\tlists/stygia.txt"
} > lists/groups.tsv

# every pair of populations within a species, then aus vs nz
{
  for sp in hilli stygia; do
    [ ${sp} = hilli ] && P="BL MKA SE DU PI WE" || P="BL CA HO DU PI WE"
    set -- ${P}
    while [ $# -gt 1 ]; do a=$1; shift; for b in "$@"; do echo -e "${sp}_${a}\t${sp}_${b}"; done; done
  done
  echo -e "hilli_AUS\thilli_NZ"
  echo -e "hilli_AUS_noMKA\thilli_NZ"
  echo -e "stygia_AUS\tstygia_NZ"
} > lists/fst_pairs.tsv
wc -l lists/*.txt lists/*.tsv | tail -40

# sensitivity: hilli country pools without the 7 flies with excess het-like sites at matched depth (f19)
printf "MKA_AUS_CH_254\nPI_K_CH_204\nPI_K_CH_205\nSE_AUS_CH_209\nSE_AUS_CH_210\nWE_CH_131\nWE_CH_141\n" > lists/hilli_highhet.txt
for g in hilli_AUS hilli_AUS_noMKA hilli_NZ; do grep -vxf lists/hilli_highhet.txt lists/${g}.txt > lists/${g}_lowhet.txt; done
printf "hilli_AUS_lowhet\thilli\tlists/hilli_AUS_lowhet.txt\nhilli_AUS_noMKA_lowhet\thilli_noMKA\tlists/hilli_AUS_noMKA_lowhet.txt\nhilli_NZ_lowhet\thilli\tlists/hilli_NZ_lowhet.txt\n" >> lists/groups.tsv
printf "hilli_AUS_lowhet\thilli_NZ_lowhet\nhilli_AUS_noMKA_lowhet\thilli_NZ_lowhet\n" >> lists/fst_pairs.tsv
