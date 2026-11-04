# auz-nz-popgen

Population genomics of two blowflies, *Calliphora hilli* and *Calliphora stygia*, sampled in Australia and New Zealand. Whole-genome resequencing of 36 flies at low coverage, analysed mainly with genotype likelihoods. PhD chapter 2.

## Samples

18 flies per species, 3 from each of 6 locations.

| species | Australia | New Zealand |
|---|---|---|
| *C. hilli* | Blackheath NSW, Mt Keira NSW, Seaford VIC | Dunedin, Pirongia, Wellington |
| *C. stygia* | Blackheath NSW, Canberra ACT, Hobart TAS | Dunedin, Pirongia, Wellington |

Sample names are `LOCATION_[AUS_]C{H,S}_ID` (CH = *C. hilli*, CS = *C. stygia*). The sample lists are in `samples.txt`, `lists/` and `angsd/lists/`.

Paired-end Illumina reads (NovaSeq, one lane). Every sample is mapped to both the *C. hilli* and the *C. stygia* masked assemblies.

## Layout

```
run_align.sh        trim + align, submits scripts/01-03
run_call.sh         all-sites calling in chunks, submits scripts/04-07
monitor.sh          live view of the alignment jobs
scripts/            read processing, hard-call vcfs, zarr stores, qc
angsd/scripts/      genotype likelihood analyses, stages s0-s17
angsd/sets.tsv      the six sample set x assembly combinations
angsd/status.sh     progress across the angsd stages
envs/               conda environments
```

Raw reads, bams, vcfs, zarr stores, genotype likelihoods and results are not in the repo.

## Read processing and hard calls (`scripts/`)

| step | script | what |
|---|---|---|
| 1 | `01_fastp.sl` | adapter and poly-g trimming, reads < 50 bp dropped |
| 2 | `02_index_ref.sl` | bwa + samtools index for each assembly |
| 3 | `03_align.sl` | bwa mem, fixmate, sort, markdup |
| 4 | `04_call_chunk.sl` | bcftools mpileup / call over all sites, 10 mb chunks |
| 5 | `05_concat.sl` | join chunks into one all-sites vcf per set |
| 6 | `06_zarr.sl` | vcf to zarr (scikit-allel), checked against the vcf |
| 7 | `07_qc.sl` | per-sample and per-site depth, missingness, heterozygosity |

```
bash run_align.sh test     # one sample on both assemblies
bash run_align.sh all
bash run_call.sh 01_hilli_masked hilli full
ZARR=1 bash run_call.sh 03_stygia_masked stygia full
```

## Genotype likelihood analyses (`angsd/scripts/`)

Most samples are at ~1-4x effective depth, too low for reliable heterozygous calls, so the population analyses use genotype likelihoods. Each stage has a `.sl` job and a summary script. Array jobs index into `angsd/sets.tsv`.

| stage | what | tools |
|---|---|---|
| s0 | integrity checks, raw reads through to the final tables | md5sum, samtools |
| s2 | competitive mapping against four *Calliphora* assemblies (species id, contamination) | bwa, samtools |
| s3 | what the unmapped reads are | kraken2, blastn |
| s4 | pca on hard calls, for comparison | scikit-allel |
| s6 | depth pilot, then genotype likelihoods in 10 mb chunks | angsd |
| s7 | pca from genotype likelihoods | pcangsd |
| s8 | admixture, k 1-5, 10 seeds per k | ngsadmix |
| s9 | pairwise relatedness | ngsrelate |
| s10 | per-sample genome-wide heterozygosity | angsd, realsfs |
| s12 | heterozygosity after downsampling, checks depth effects | samtools, angsd |
| s13 | hybrid index and interspecific heterozygosity (triangle plot) | python |
| s14 | hybrid class assignment, 10 replicates of 200 / 500 random loci | newhybrids |
| s15 | fit of the admixture models | evaladmix |
| s16 | evidence for the low-mapping samples in one table | python |
| s17 | robustness of the gl pca (split halves, downsampling) | angsd, pcangsd |

The hybrid analysis follows Croft et al. 2024 (*Ecology and Evolution*): combined pca, admixture at k 2-5 and NewHybrids on 200 loci x 10 replicates. Here each step runs on genotype likelihoods instead of hard calls.

## Software

Run on NeSI (Slurm). The environments are in `envs/`:

- `ch2_zarr`: bwa 0.7.19, samtools / bcftools 1.24, scikit-allel, zarr
- `ch2_gl`: pcangsd 1.36.4, ngsrelate 2.3.1, numpy, pandas, matplotlib

NeSI modules: angsd 0.935 (with ngsadmix and realsfs), fastp 0.23.4, Kraken2 2.1.3, BLAST 2.17.0, SAMtools 1.23.1. evalAdmix and NewHybrids are compiled from source.

```
conda env create -f envs/ch2_zarr.yaml
conda env create -f envs/ch2_gl.yaml
```

Paths in the scripts point to the NeSI project folders and need changing to run elsewhere.
