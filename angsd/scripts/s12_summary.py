# stage 12 summary: het at each downsampled depth next to the full-depth stage 10 value
import glob
import os

import numpy as np
import pandas as pd

pd.set_option("display.width", 250)
rows = []
for f in glob.glob("het_ds/*.ml"):
    if not os.path.exists(f.replace(".ml", ".done")):
        continue
    s, t = os.path.basename(f)[:-3].split("__")
    a = np.loadtxt(f, ndmin=2)
    rows.append({"sample": s, "target": "het_at_" + t, "het": a[:, 1].sum() / a.sum()})
d = pd.DataFrame(rows).pivot(index="sample", columns="target", values="het").reset_index()
full = pd.read_csv("tables/S10_heterozygosity.tsv", sep="\t")[["sample", "species", "effDP_own", "het_own"]]
d = full.merge(d, on="sample").sort_values(["species", "het_own"])
d.round(5).to_csv("tables/S12_downsampled_het.tsv", sep="\t", index=False)
print(d.round(5).to_string(index=False))
print("written to tables/S12_downsampled_het.tsv")
