# stage 7: pcangsd covariance -> pc coordinates (same layout as stage 4 so one plot script does both)
# also prints admixture and per-sample inbreeding next to names
# usage: python s7_cov2pcs.py <out_prefix> <sample_list>
import glob
import sys

import numpy as np
import pandas as pd

out, lst = sys.argv[1:3]
names = [l.strip() for l in open(lst) if l.strip()]
meta = pd.read_csv("samples_meta.tsv", sep="\t")
C = np.loadtxt(out + ".cov")
val, vec = np.linalg.eigh(C)
order = val.argsort()[::-1]
val, vec = np.clip(val[order], 0, None), vec[:, order]
k = min(10, len(names))
pcs = pd.DataFrame(vec[:, :k] * np.sqrt(val[:k]), columns=["PC%d" % (i + 1) for i in range(k)])
pcs.insert(0, "sample", names)
pcs = meta.merge(pcs, on="sample")
pcs.round(5).to_csv(out + ".pcs.tsv", sep="\t", index=False)
var = 100 * val / val.sum()
pd.DataFrame({"PC": range(1, len(var) + 1), "pct_var": var.round(3)}).to_csv(out + ".var.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print("%s: variance explained PC1-PC5 (%%): %s" % (out, ", ".join("%.1f" % v for v in var[:5])))
tab = pcs[["sample", "species", "country", "PC1", "PC2"]].copy()
for q in sorted(glob.glob(out + ".admix.*.Q")):
    Q = np.loadtxt(q, ndmin=2)
    for j in range(Q.shape[1]):
        tab["Q%d_of_%d" % (j + 1, Q.shape[1])] = Q[:, j]
try:
    tab["F_inbreed"] = np.loadtxt(out + ".inbreed.samples")
except OSError:
    pass
print(tab.sort_values("PC1").round(3).to_string(index=False))
