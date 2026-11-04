#!/bin/bash
# usage: bash run_call.sh <REF_NAME> <SAMPLE_SET> test|full
#   REF_NAME     01_hilli_masked or 03_stygia_masked
#   SAMPLE_SET   all, hilli or stygia (reads lists/<SAMPLE_SET>.txt)
#   test         calls the first 2 Mb of the longest scaffold only
#   full         calls the whole genome in chunks, then joins them (re-run it to resubmit only missing chunks)
cd "$(dirname "$(readlink -f "$0")")" || exit 1
REF=$1; SET=$2; MODE=$3
CHUNK_BP=${CHUNK_BP:-10000000}
if [ -z "$REF" ] || [ -z "$SET" ] || { [ "$MODE" != test ] && [ "$MODE" != full ]; }; then
  echo "usage: bash run_call.sh <REF_NAME> <SAMPLE_SET> test|full"; exit 1
fi
[ -s ref/$REF.fa.fai ] || { echo "missing ref/$REF.fa.fai"; exit 1; }
[ -s lists/$SET.txt ] || { echo "missing lists/$SET.txt"; exit 1; }

REFTAG=$(echo $REF | sed -E 's/^[0-9]+_//; s/_masked$//')
TAG=${SET}_on_${REFTAG}
OUT=vcf/$TAG
mkdir -p $OUT/chunks $OUT/test logs regions/$REF

# every sample in the list needs a finished bam on this reference
BAD=0
while read -r S; do
  [ -f bam/$REF/$S.done ] || { echo "no finished BAM for $S on $REF"; BAD=1; }
done < lists/$SET.txt
[ $BAD = 0 ] || exit 1
sed "s#^#bam/$REF/#; s#\$#.bam#" lists/$SET.txt > $OUT/bams.txt
echo "$TAG: $(wc -l < lists/$SET.txt) samples on $REF"

# test region: first 2 Mb of the longest scaffold
sort -k2,2nr ref/$REF.fa.fai | head -1 | awk '{e=($2<2000000)?$2:2000000; printf "%s\t1\t%d\n", $1, e}' > regions/$REF/test.tsv

if [ "$MODE" = test ]; then
  rm -f $OUT/test/test.vcf.gz $OUT/test/test.vcf.gz.tbi $OUT/test/test.done
  sbatch --job-name=calltest_$TAG --array=1 --export=ALL,REF_NAME=$REF,SAMPLE_SET=$SET,MODE=test scripts/04_call_chunk.sl
  exit 0
fi

# chunks of CHUNK_BP bases in reference order (small scaffolds are grouped, long ones are split)
rm -f regions/$REF/chunk_*.tsv
awk -v cs=$CHUNK_BP -v d=regions/$REF 'BEGIN{n=1; cap=cs}
  {pos=1; while (pos<=$2) {take=$2-pos+1; if (take>cap) take=cap;
     f=sprintf("%s/chunk_%04d.tsv", d, n); printf "%s\t%d\t%d\n", $1, pos, pos+take-1 > f;
     pos+=take; cap-=take; if (cap==0) {close(f); n++; cap=cs}}}' ref/$REF.fa.fai
N=$(ls regions/$REF/chunk_*.tsv | wc -l)
echo "$N chunks of up to $CHUNK_BP bases"

TODO=""
for i in $(seq 1 $N); do
  n=$(printf "%04d" $i)
  [ -f $OUT/chunks/chunk_$n.done ] || TODO="$TODO,$i"
done
TODO=${TODO#,}

DEP=""
if [ -n "$TODO" ]; then
  JOB=$(sbatch --parsable --job-name=call_$TAG --array=${TODO}%10 --export=ALL,REF_NAME=$REF,SAMPLE_SET=$SET,MODE=full scripts/04_call_chunk.sl)
  echo "calling array $JOB, chunks: $TODO"
  DEP="--dependency=afterok:$JOB"
else
  echo "all $N chunks are already done"
fi
CJOB=$(sbatch --parsable --job-name=concat_$TAG $DEP --export=ALL,REF_NAME=$REF,SAMPLE_SET=$SET,N_CHUNKS=$N scripts/05_concat.sl)
echo "join job $CJOB"
if [ -n "$ZARR" ]; then
  ZJOB=$(sbatch --parsable --job-name=zarr_$TAG --dependency=afterok:$CJOB --export=ALL,TAG=$TAG scripts/06_zarr.sl)
  QJOB=$(sbatch --parsable --job-name=qc_$TAG --dependency=afterok:$ZJOB --export=ALL,TAG=$TAG scripts/07_qc.sl)
  echo "zarr job $ZJOB, qc job $QJOB"
fi
