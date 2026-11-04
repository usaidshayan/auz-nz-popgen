#!/bin/bash -e
#SBATCH --job-name=s17_pcangsd
#SBATCH --account=uow04177
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --output=logs/%x_%j.out
# pcangsd for the s17 checks, same settings as s7
# MODE=half: each half of hilli_on_hilli and stygia_on_stygia
# MODE=ds: the two downsampled sets, plus the original gls cut to the same sites as the baseline
cd "${SLURM_SUBMIT_DIR}"
GL=/nesi/project/uow04177/${USER}/conda_envs/ch2_gl
run() {
  ${GL}/bin/pcangsd -b $1 -o $2 -t ${SLURM_CPUS_PER_TASK} --maf 0.05
  ${GL}/bin/python scripts/s7_cov2pcs.py $2 $3
}
if [ "${MODE}" = half ]; then
  for TAG in hilli_on_hilli stygia_on_stygia; do
    L=$(awk -v t=${TAG} '$1==t {print $3}' sets.tsv)
    for H in A B; do run robust/half/${TAG}.${H}.beagle.gz robust/half/${TAG}.${H} ${L}; done
  done
else
  for S in WE_CH_131 BL_AUS_CH_214; do
    D=robust/gl/ds_${S}
    ( zcat $(ls ${D}/chunk_*.beagle.gz | head -1) | head -1
      for f in ${D}/chunk_*.beagle.gz; do zcat $f | tail -n +2; done ) | gzip -c > ${D}.beagle.gz
    echo "ds_${S}: $(zcat ${D}.beagle.gz | tail -n +2 | wc -l) sites"
  done
  # sites present in both downsampled sets, original gls at those sites
  zcat robust/gl/ds_WE_CH_131.beagle.gz | cut -f1 | tail -n +2 | sort > robust/gl/a.txt
  zcat robust/gl/ds_BL_AUS_CH_214.beagle.gz | cut -f1 | tail -n +2 | sort > robust/gl/b.txt
  comm -12 robust/gl/a.txt robust/gl/b.txt > robust/gl/common.txt
  rm robust/gl/a.txt robust/gl/b.txt
  for X in base:gl/hilli_on_hilli/hilli_on_hilli.beagle.gz ds_WE_CH_131:robust/gl/ds_WE_CH_131.beagle.gz \
           ds_BL_AUS_CH_214:robust/gl/ds_BL_AUS_CH_214.beagle.gz; do
    N=${X%%:*}; B=${X#*:}
    zcat ${B} | awk 'NR == FNR {k[$1]; next} FNR == 1 || ($1 in k)' robust/gl/common.txt - \
      | gzip -c > robust/gl/${N}.common.beagle.gz
    echo "${N}: $(zcat robust/gl/${N}.common.beagle.gz | tail -n +2 | wc -l) common sites"
    run robust/gl/${N}.common.beagle.gz robust/gl/${N} lists/hilli.txt
  done
fi
