#!/bin/bash -e
#SBATCH --job-name=s3_kraken
#SBATCH --account=uow04177
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=230G
#SBATCH --output=logs/%x_%j.out
# what the reads that map to none of the four calliphora assemblies are (stage 2 unmapped pairs)
# kraken2 against nesi nt (has eukaryotes, so fly reads not in the assemblies still get placed)
# db loaded once, low-mapping eight first, the rest as the baseline
# confidence 0.1 to cut spurious hits
cd "${SLURM_SUBMIT_DIR}"
module purge
module load Kraken2/2.1.3-GCC-11.3.0
DB=/opt/nesi/db/Kraken2/nt
mkdir -p compmap/kraken
kraken2 --version | head -1
LOW="WE_CS_133 DU_F_CH_26 WE_CS_139 MKA_AUS_CH_211 WE_CH_141 DU_F_CH_28 DU_F_CH_24 HO_AUS_CS_204"
for S in ${LOW} $(grep -vxF -f <(tr ' ' '\n' <<< "${LOW}") samples.txt); do
  if [ -s compmap/kraken/${S}.report ]; then echo "${S} done, skipping"; continue; fi
  kraken2 --db ${DB} --threads ${SLURM_CPUS_PER_TASK} --paired --gzip-compressed --confidence 0.1 \
    --report compmap/kraken/${S}.report.tmp --output /dev/null \
    compmap/unmapped/${S}_R1.fq.gz compmap/unmapped/${S}_R2.fq.gz 2> compmap/kraken/${S}.log
  mv compmap/kraken/${S}.report.tmp compmap/kraken/${S}.report
  echo "${S}: $(grep -m1 'sequences classified' compmap/kraken/${S}.log)"
done
/nesi/project/uow04177/${USER}/conda_envs/ch2_gl/bin/python scripts/s3_summary.py
