#!/bin/bash -e
#SBATCH --job-name=s17_prep
#SBATCH --account=uow04177
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%j.out
# how much to trust the gl pca
# 1. share of snp sites where each sample has any reads (gl triplet not flat)
# 2. split each beagle into alternating 1 mb blocks, pcangsd on each half separately
# 3. two deep hilli (nz WE_CH_131, aus BL_AUS_CH_214) downsampled to WE_CH_143's depth, gls redone at the same sites
cd "${SLURM_SUBMIT_DIR}"
SAMTOOLS=/nesi/project/uow04177/${USER}/conda_envs/ch2_zarr/bin/samtools
module purge
module load angsd/0.935-GCC-9.2.0
mkdir -p robust/half robust/gl robust/ds

for TAG in hilli_on_hilli stygia_on_stygia; do
  (
  zcat gl/${TAG}/${TAG}.beagle.gz | awk -v tag=${TAG} -v a="gzip -c > robust/half/${TAG}.A.beagle.gz" \
    -v b="gzip -c > robust/half/${TAG}.B.beagle.gz" -v inf=robust/${TAG}.informative.tsv '
    NR == 1 {print | a; print | b; n = (NF - 3) / 3; next}
    {for (j = 0; j < n; j++) {i = 4 + 3 * j; if ($i != $(i + 1) || $i != $(i + 2)) c[j]++}
     t++
     p = $1; sub(/.*_/, "", p); key = substr($1, 1, length($1) - length(p) - 1) ":" int(p / 1e6)
     if (key != last) {k++; last = key}
     if (k % 2) print | a; else print | b}
    END {for (j = 0; j < n; j++) printf "%s\t%d\t%d\t%d\n", tag, j, c[j], t > inf}'
  ) &
done

TARGET=$(awk -F'\t' '$1=="WE_CH_143" {print $6}' ../qc/bam_summary_01_hilli_masked.tsv)
for S in WE_CH_131 BL_AUS_CH_214; do
  MEAN=$(awk -F'\t' -v s=${S} '$1==s {print $6}' ../qc/bam_summary_01_hilli_masked.tsv)
  FRAC=$(awk -v t=${TARGET} -v m=${MEAN} 'BEGIN {printf "%.4f", t / m}')
  echo "${S}: depth ${MEAN} -> ${TARGET} (keep ${FRAC} of pairs)"
  if [ ! -f robust/ds/${S}.bam.bai ]; then
    ${SAMTOOLS} view -b -@ 6 --subsample ${FRAC} --subsample-seed 11 -o robust/ds/${S}.bam bam_pvp/01_hilli_masked/${S}.bam
    ${SAMTOOLS} index robust/ds/${S}.bam
  fi
  sed "s#bam_pvp/01_hilli_masked/${S}.bam#robust/ds/${S}.bam#" gl/hilli_on_hilli/bams.txt > robust/gl/ds_${S}.bams.txt
done

zcat gl/hilli_on_hilli/hilli_on_hilli.mafs.gz | awk 'NR > 1 {print $1"\t"$2"\t"$3"\t"$4}' > robust/gl/hilli_on_hilli.sites
angsd sites index robust/gl/hilli_on_hilli.sites
cut -f1 robust/gl/hilli_on_hilli.sites | uniq > robust/gl/hilli_on_hilli.sites.chr
for S in WE_CH_131 BL_AUS_CH_214; do
  for C in $(ls gl/hilli_on_hilli/chunks/*.rf | xargs -n1 basename | sed 's/.rf$//'); do echo -e "ds_${S}\t${C}"; done
done > robust/gl/tasks.tsv
wait
echo "tasks: $(wc -l < robust/gl/tasks.tsv)"
cat robust/*.informative.tsv
