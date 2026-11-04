import glob, os, re
import pandas as pd

pd.set_option("display.width", 250)
REFS = {"hilli": "01_hilli_masked", "stygia": "03_stygia_masked"}


def count(pattern, text):
    m = re.search(pattern, text, re.M)
    return int(m.group(1)) if m else 0


def summarise(ref):
    rows = []
    for flag_path in sorted(glob.glob("qc/bam/%s/*.flagstat.txt" % ref)):
        sample = os.path.basename(flag_path).replace(".flagstat.txt", "")
        cov_path = "qc/bam/%s/%s.coverage.tsv" % (ref, sample)
        if not os.path.exists("bam/%s/%s.done" % (ref, sample)) or not os.path.exists(cov_path):
            continue
        txt = open(flag_path).read()
        primary = count(r"^(\d+) \+ \d+ primary$", txt)
        mapped = count(r"^(\d+) \+ \d+ primary mapped", txt)
        proper = count(r"^(\d+) \+ \d+ properly paired", txt)
        dups = count(r"^(\d+) \+ \d+ primary duplicates", txt)
        total_len = depth_sum = covered = 0
        with open(cov_path) as fh:
            next(fh)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                length = int(f[2]) - int(f[1]) + 1
                total_len += length
                depth_sum += float(f[6]) * length
                covered += int(f[4])
        rows.append([sample, primary, round(100.0 * mapped / primary, 2), round(100.0 * proper / primary, 2),
                     round(100.0 * dups / primary, 2), round(depth_sum / total_len, 2),
                     round(100.0 * covered / total_len, 2)])
    cols = ["sample", "primary_reads", "pct_mapped", "pct_properly_paired", "pct_duplicates", "mean_depth", "pct_genome_covered"]
    return pd.DataFrame(rows, columns=cols)


meta = pd.read_csv("samples_meta.tsv", sep="\t")
tab = meta
for name, ref in REFS.items():
    d = summarise(ref)
    d.to_csv("qc/bam_summary_%s.tsv" % ref, sep="\t", index=False)
    print("%s: %d of %d samples finished" % (ref, len(d), len(meta)))
    d = d[["sample", "pct_mapped", "pct_properly_paired", "pct_duplicates", "mean_depth"]]
    d.columns = ["sample"] + ["%s_%s" % (c, name) for c in d.columns[1:]]
    tab = tab.merge(d, on="sample", how="left")
tab["mapped_diff_hilli_minus_stygia"] = (tab["pct_mapped_hilli"] - tab["pct_mapped_stygia"]).round(2)
tab.to_csv("qc/ref_compare.tsv", sep="\t", index=False)

cols = ["pct_mapped_hilli", "pct_mapped_stygia", "pct_properly_paired_hilli", "pct_properly_paired_stygia",
        "mean_depth_hilli", "mean_depth_stygia"]
print("\nMean per species (columns) for each measure (rows):")
print(tab.groupby("species")[cols].mean().round(1).T.to_string())

tab["best_mapped"] = tab[["pct_mapped_hilli", "pct_mapped_stygia"]].max(axis=1)
show = ["sample", "species", "location", "pct_mapped_hilli", "pct_mapped_stygia", "mean_depth_hilli", "mean_depth_stygia"]
print("\nSix samples with the lowest best-case mapping rate:")
print(tab.sort_values("best_mapped").head(6)[show].to_string(index=False))
print("\nPer-sample table written to qc/ref_compare.tsv")
