# integrity checks, run from s0_check.sl after the fastq md5s and bam headers are written to check/
# prints pass / warn / fail per check, a summary count at the end
import glob
import gzip
import os
import re

import numpy as np
import pandas as pd

REFS = ["01_hilli_masked", "03_stygia_masked"]
res = []


def rep(status, name, detail=""):
    res.append(status)
    print("%-4s %s%s" % (status, name, (": " + detail) if detail else ""), flush=True)


def chk(cond, name, detail=""):
    rep("PASS" if cond else "FAIL", name, detail)


# samples and metadata
samples = [l.strip() for l in open("samples.txt") if l.strip()]
meta = pd.read_csv("samples_meta.tsv", sep="\t")
chk(len(samples) == 36 and len(set(samples)) == 36, "36 unique samples in samples.txt")
chk(sorted(meta["sample"]) == sorted(samples), "samples_meta.tsv matches samples.txt")
code = meta["sample"].str.extract(r"_(C[HS])_")[0].map({"CH": "hilli", "CS": "stygia"})
chk((code == meta.species).all(), "species in metadata matches CH / CS code in sample name",
    ", ".join(meta["sample"][code != meta.species]))
print(meta.groupby(["species", "country"]).size().to_string())
m = meta.set_index("sample")
for f in sorted(glob.glob("lists/*.txt")):
    lst = [l.strip() for l in open(f) if l.strip()]
    parts = os.path.basename(f)[:-4].split("_")
    want = m.index
    if parts[0] != "all":
        want = want[m.species == parts[0]]
    if len(parts) > 1:
        want = want[m.loc[want, "country"] == parts[1]]
    chk(sorted(lst) == sorted(want) and len(set(lst)) == len(lst), "list %s matches metadata" % f, "%d samples" % len(lst))

# fastq
fq = pd.read_csv("check/fastq.tsv", sep="\t", header=None, names=["file", "dir", "md5", "lines", "gz"])
fq["kind"] = np.where(fq.dir.str.contains("raw_data"), "raw", "trimmed")
fq["sample"] = fq.file.str.replace(r"(_DKDL.*|_R[12]\.fq\.gz)$", "", regex=True)
fq["read"] = np.where(fq.file.str.contains(r"(_1\.fq|_R1\.fq)"), "R1", "R2")
chk(len(fq) == 144, "144 fastq files (36 x 2 raw + 36 x 2 trimmed)", "%d found" % len(fq))
chk((fq.gz == "ok").all(), "all fastq gzip streams decompress", ", ".join(fq.file[fq.gz != "ok"]))
chk((fq.lines % 4 == 0).all(), "all fastq line counts divisible by 4", ", ".join(fq.file[fq.lines % 4 != 0]))
dup = fq[fq.md5.duplicated(keep=False)]
chk(dup.empty, "no two fastq files identical (md5)", ", ".join(dup.file))
chk(set(fq["sample"]) == set(samples), "fastq sample names match samples.txt")
w = fq.pivot_table(index="sample", columns=["kind", "read"], values="lines", aggfunc="sum") // 4
for k in ["raw", "trimmed"]:
    bad = w.index[w[(k, "R1")] != w[(k, "R2")]]
    chk(len(bad) == 0, "%s R1 and R2 read counts equal" % k, ", ".join(bad))
fp = pd.read_csv("../qc/fastp_summary.tsv", sep="\t").set_index("sample")
raw_tot = w[("raw", "R1")] + w[("raw", "R2")]
trim_tot = w[("trimmed", "R1")] + w[("trimmed", "R2")]
chk((raw_tot == fp.reads_before.loc[raw_tot.index]).all(), "raw reads match fastp reads_before")
chk((trim_tot == fp.reads_after.loc[trim_tot.index]).all(), "trimmed reads match fastp reads_after")
fq[["kind", "sample", "read", "file", "md5", "lines"]].to_csv("check/fastq_md5.tsv", sep="\t", index=False)


# references
def fasta_lengths(path):
    L, name, n = {}, None, 0
    with open(path) as f:
        for line in f:
            if line[0] == ">":
                if name:
                    L[name] = n
                name, n = line[1:].split()[0], 0
            else:
                n += len(line.strip())
    L[name] = n
    return L


fai = {}
for r in REFS:
    fa = fasta_lengths("ref/%s.fa" % r)
    fi = pd.read_csv("ref/%s.fa.fai" % r, sep="\t", header=None, usecols=[0, 1], names=["c", "l"])
    fai[r] = fi
    chk(list(fa) == list(fi.c) and list(fa.values()) == list(fi.l), "%s fai matches fasta" % r,
        "%d contigs, %.1f mb" % (len(fa), sum(fa.values()) / 1e6))
    chk(os.path.getmtime("../ref/%s.fa.bwt" % r) > os.path.getmtime("ref/%s.fa" % r), "%s bwa index newer than fasta" % r)
cr = [l[1:].split()[0] for l in open("ref/compref.fa") if l[0] == ">"]
for r, p in [("01_hilli_masked", "hilli__"), ("03_stygia_masked", "stygia__")]:
    sub = [c[len(p):] for c in cr if c.startswith(p)]
    chk(sub == list(fai[r].c), "competitive ref holds every %s contig" % r, "%d" % len(sub))

# bams
bh = pd.read_csv("check/bam_headers.tsv", sep="\t", header=None, names=["ref", "sample", "quick", "index"])
chk(len(bh) == 72, "72 bams (36 x 2 refs)")
chk((bh.quick == "ok").all(), "all bams pass samtools quickcheck", ", ".join(bh["sample"][bh.quick != "ok"]))
chk((bh["index"] == "ok").all(), "all bam indexes newer than bam", ", ".join(bh["sample"][bh["index"] != "ok"]))
bad = []
for _, b in bh.iterrows():
    h = open("check/hdr_%s_%s.txt" % (b.ref, b["sample"])).read().splitlines()
    sq = [(re.search(r"SN:(\S+)", l).group(1), int(re.search(r"LN:(\d+)", l).group(1))) for l in h if l.startswith("@SQ")]
    rg = [l for l in h if l.startswith("@RG")]
    cl = [l for l in h if l.startswith("@PG\tID:bwa")]
    s = b["sample"]
    probs = []
    if not h[0].startswith("@HD") or "SO:coordinate" not in h[0]:
        probs.append("not coordinate sorted")
    if sq != list(zip(fai[b.ref].c, fai[b.ref].l)):
        probs.append("@SQ differs from fai")
    if len(rg) != 1 or "SM:%s\t" % s not in rg[0] + "\t":
        probs.append("read group sample wrong")
    if len(cl) != 1 or ("ref/%s.fa" % b.ref) not in cl[0] or ("trimmed/%s_R1.fq.gz" % s) not in cl[0] \
            or ("trimmed/%s_R2.fq.gz" % s) not in cl[0]:
        probs.append("bwa command line names wrong ref or reads")
    if not any("markdup" in l for l in h):
        probs.append("no markdup")
    fs = open("../qc/bam/%s/%s.flagstat.txt" % (b.ref, s)).read()
    prim = int(re.search(r"(\d+) \+ \d+ primary\n", fs).group(1))
    if prim != trim_tot[s]:
        probs.append("primary reads %d != trimmed reads %d" % (prim, trim_tot[s]))
    if probs:
        bad.append("%s/%s (%s)" % (b.ref, s, "; ".join(probs)))
chk(not bad, "bam headers: sorted, @SQ = fai, read group, bwa ref and reads, markdup, all trimmed reads present",
    " | ".join(bad))
fs_time = max(os.path.getmtime("bam_pvp/%s/%s.bam" % (r, s)) for r in REFS for s in samples)
for r in REFS:
    chk(os.path.getmtime("ref/%s.fa" % r) < fs_time, "%s fasta older than the bams" % r)

# regions used for the gl chunks
for r in REFS:
    rg = pd.concat([pd.read_csv(f, sep="\t", header=None, names=["c", "s", "e"])
                    for f in sorted(glob.glob("regions_pvp/%s/chunk_*.tsv" % r))])
    L = fai[r].set_index("c").l
    per = (rg.e - rg.s + 1).groupby(rg.c).sum()
    full = (per == L.reindex(per.index)).all()
    cov = (rg.e - rg.s + 1).sum() / L.sum()
    chk(full and rg.c.nunique() == len(L) and not rg.duplicated(["c", "s"]).any(),
        "%s gl chunks tile the whole assembly once" % r, "%.4f of bases, %d contigs" % (cov, rg.c.nunique()))

# genotype likelihoods
sets = pd.read_csv("sets.tsv", sep="\t", header=None, names=["tag", "ref", "list", "zarr"])
tasks = pd.read_csv("gl/tasks.tsv", sep="\t", header=None, names=["tag", "ref", "chunk"])
gl_sites = {}
for _, s in sets.iterrows():
    lst = [l.strip() for l in open(s.list) if l.strip()]
    bams = [l.strip() for l in open("gl/%s/bams.txt" % s.tag) if l.strip()]
    chk(bams == ["bam_pvp/%s/%s.bam" % (s.ref, x) for x in lst], "%s bams.txt = %s order on %s" % (s.tag, s.list, s.ref))
    nt = (tasks.tag == s.tag).sum()
    nd = len(glob.glob("gl/%s/chunks/*.done" % s.tag))
    chk(nt == nd, "%s all gl chunks finished" % s.tag, "%d of %d" % (nd, nt))
    pos_b, n, gmax, nan, bad_hdr = [], 0, 0.0, 0, False
    for i, ch in enumerate(pd.read_csv("gl/%s/%s.beagle.gz" % (s.tag, s.tag), sep="\t", chunksize=500000)):
        if i == 0:
            bad_hdr = list(ch.columns[:3]) != ["marker", "allele1", "allele2"] or len(ch.columns) != 3 + 3 * len(lst)
        g = ch.iloc[:, 3:].to_numpy(np.float64).reshape(len(ch), len(lst), 3)
        nan += int(np.isnan(g).sum())
        gmax = max(gmax, float(np.nanmax(np.abs(g.sum(2) - 1))))
        pos_b.append(ch.marker.to_numpy())
        n += len(ch)
    chk(not bad_hdr, "%s beagle header has %d individuals" % (s.tag, len(lst)))
    chk(nan == 0 and gmax < 1e-3, "%s beagle gls finite and sum to 1" % s.tag, "max deviation %.2g" % gmax)
    mk = np.concatenate(pos_b)
    mf = pd.read_csv("gl/%s/%s.mafs.gz" % (s.tag, s.tag), sep="\t", usecols=[0, 1])
    chk(len(mf) == n and (mk == (mf.chromo + "_" + mf.position.astype(str)).to_numpy()).all(),
        "%s beagle and mafs hold the same sites" % s.tag, "%d snps" % n)
    order = {c: k for k, c in enumerate(fai[s.ref].c)}
    ci = mf.chromo.map(order)
    L = fai[s.ref].set_index("c").l
    sorted_ok = ci.notna().all() and ((ci.diff().fillna(0) > 0) | ((ci.diff().fillna(0) == 0) & (mf.position.diff().fillna(1) > 0))).all()
    chk(sorted_ok and (mf.position <= mf.chromo.map(L)).all(), "%s sites sorted, unique and inside contigs" % s.tag)
    gl_sites[s.tag] = n
    del mk, mf, pos_b

# downstream
for _, s in sets.iterrows():
    lst = [l.strip() for l in open(s.list) if l.strip()]
    lg = open("pcangsd/%s.log" % s.tag).read()
    mm = re.search(r"Loaded (\d+) sites and (\d+) samples", lg)
    chk(int(mm.group(1)) == gl_sites[s.tag] and int(mm.group(2)) == len(lst), "%s pcangsd read the full beagle" % s.tag)
    pcs = pd.read_csv("pcangsd/%s.pcs.tsv" % s.tag, sep="\t")
    cov = np.loadtxt("pcangsd/%s.cov" % s.tag)
    chk(list(pcs["sample"]) == lst and cov.shape == (len(lst), len(lst)), "%s pca sample labels follow beagle order" % s.tag)
    rel = pd.read_csv("relate/%s.res" % s.tag, sep="\t")
    chk(len(rel) == len(lst) * (len(lst) - 1) // 2 and set(rel.ida) | set(rel.idb) == set(lst),
        "%s ngsrelate has every pair" % s.tag)
    top = rel.sort_values("KING", ascending=False).iloc[0]
    rep("PASS" if top.KING < 0.354 else "FAIL", "%s no duplicate samples (KING < 0.354)" % s.tag,
        "max KING %.3f %s-%s, max rab %.3f" % (top.KING, top.ida, top.idb, rel.rab.max()))
for tag in ["all_on_hilli", "all_on_stygia"]:
    logs = sorted(glob.glob("admix/%s.K*.s*.log" % tag))
    dims = {re.search(r"AFTER filtering\): nsites=(\d+) nind=(\d+)", open(f).read()).groups() for f in logs}
    thin = sum(1 for _ in gzip.open("admix/%s.thin5.beagle.gz" % tag, "rt")) - 1
    chk(len(logs) == 50 and len(dims) == 1 and list(dims)[0][1] == "36",
        "%s ngsadmix 50 runs, same sites and 36 individuals" % tag, "thinned %d, after maf %s" % (thin, list(dims)[0][0]))
    chk(abs(thin - gl_sites[tag] / 5) <= 1, "%s thinned beagle is every 5th site of the full beagle" % tag)
    h = pd.read_csv("hybrid/%s.hybrid_index.tsv" % tag, sep="\t")
    chk(list(h["sample"]) == [l.strip() for l in open("lists/all.txt") if l.strip()] and
        (h.species.values == m.loc[h["sample"], "species"].values).all(), "%s hybrid index rows match samples and species" % tag)
for f in ["tables/S2_competitive_mapping.tsv", "tables/T1_mapping_depth.tsv", "tables/S3_unmapped_kraken.tsv",
          "tables/S10_heterozygosity.tsv", "tables/S8_ngsadmix_K2_q.tsv"]:
    t = pd.read_csv(f, sep="\t")
    chk(sorted(t["sample"]) == sorted(samples) and (t.species.values == m.loc[t["sample"], "species"].values).all(),
        "%s has all 36 samples with the right species" % f)
for r in ["all_on_hilli", "all_on_stygia", "hilli_on_hilli", "stygia_on_stygia"]:
    a = [l.strip() for l in open("vcf_pvp/%s/samples_in_vcf.txt" % r)]
    b = [l.strip() for l in open("vcf_pvp/%s/samples_expected.txt" % r)]
    v = "vcf_pvp/%s/%s.allsites.vcf.gz" % (r, r)
    chk(a == b and os.path.getmtime(v + ".tbi") >= os.path.getmtime(v), "%s vcf samples and index ok" % r)

print("\n%d pass, %d warn, %d fail" % (res.count("PASS"), res.count("WARN"), res.count("FAIL")))
