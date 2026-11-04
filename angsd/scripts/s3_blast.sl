#!/bin/bash -e
#SBATCH --job-name=s3_blast
#SBATCH --account=uow04177
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=logs/%x_%j.out
# kraken left most unmapped reads unclassified in some samples (we_cs_139 99%, mka_211 68%, ho_204 69%)
# its nt db is from 2020, so blast 500 random unmapped r1 reads per sample against the current nt
# se_aus_ch_209 as a normal-mapping comparison
cd "${SLURM_SUBMIT_DIR}"
module purge
module load BLAST/2.17.0-GCC-15.2.0
DB=/opt/nesi/db/blast/2026-08
export BLASTDB=${DB}
GL=/nesi/project/uow04177/${USER}/conda_envs/ch2_gl
mkdir -p compmap/blast
for S in WE_CS_139 MKA_AUS_CH_211 HO_AUS_CS_204 WE_CS_133 SE_AUS_CH_209; do
  Q=compmap/blast/${S}.sub.fa
  zcat compmap/unmapped/${S}_R1.fq.gz | ${GL}/bin/python -c "
import sys, random
random.seed(1)
recs = []
for i, line in enumerate(sys.stdin):
    if i % 4 == 0: name = line[1:].split()[0]
    elif i % 4 == 1:
        r = (name, line.strip())
        if len(recs) < 500: recs.append(r)
        else:
            j = random.randint(0, i // 4)
            if j < 500: recs[j] = r
for n, s in recs: print('>%s\n%s' % (n, s))
" > ${Q}
  blastn -task megablast -db ${DB}/nt -query ${Q} -evalue 1e-10 -max_target_seqs 5 -max_hsps 1 \
    -num_threads ${SLURM_CPUS_PER_TASK} \
    -outfmt "6 qseqid sseqid pident length evalue bitscore staxids sscinames sskingdom" \
    > compmap/blast/${S}.blast.tsv
  echo "${S}: $(cut -f1 compmap/blast/${S}.blast.tsv | sort -u | wc -l) of 500 reads hit"
done
