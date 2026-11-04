# stage 10 summary: angsd het on own vs other assembly, next to hard-call het (dp >= 4)
# from the main sample qc (../qc), ordered by depth
import glob
import os

import numpy as np
import pandas as pd

pd.set_option("display.width", 250)
REF = {"hilli": "01_hilli_masked", "stygia": "03_stygia_masked"}
meta = pd.read_csv("samples_meta.tsv", sep="\t")
rows = []
for f in glob.glob("het/*/*.ml"):
    if not os.path.exists(f.replace(".ml", ".done")):
        continue
    a = np.loadtxt(f, ndmin=2)
    blocks = a[:, 1] / a.sum(axis=1)
    rows.append({"sample": os.path.basename(f)[:-3], "ref": os.path.basename(os.path.dirname(f)), "sites_M": a.sum() / 1e6,
                 "het": a[:, 1].sum() / a.sum(), "het_block_sd": blocks.std()})
h = pd.DataFrame(rows).merge(meta[["sample", "species", "country"]], on="sample")
h["assembly"] = np.where(h["ref"] == h["species"].map(REF), "own", "other")
w = h.pivot_table(index=["sample", "species", "country"], columns="assembly", values=["het", "sites_M"]).reset_index()
w.columns = ["_".join(c).strip("_") for c in w.columns]
t1 = pd.read_csv("tables/T1_mapping_depth.tsv", sep="\t")[["sample", "effDP_own"]]
hard = pd.concat([pd.read_csv("../qc/sample_qc_%s_on_%s.tsv" % (s, s), sep="\t") for s in REF])
w = w.merge(t1, on="sample").merge(hard[["sample", "het_per_site_DP_ge4"]].drop_duplicates("sample"), on="sample")
w = w.sort_values("effDP_own")
w.round(6).to_csv("tables/S10_heterozygosity.tsv", sep="\t", index=False)
print(w.round(5).to_string(index=False))
for col in [c for c in w.columns if c.startswith("het")]:
    print("correlation of %s with effective depth: r = %.2f" % (col, w[col].corr(w["effDP_own"])))
print("written to tables/S10_heterozygosity.tsv")
