# pc1 vs pc2, six panels for one pca type. colour = location, marker = species
# labels the five <1.7x samples and anything > 2.5 sd from its species centre
# usage: python pca_plot.py <dir> <suffix e.g. .pcs.tsv> <out.png>
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

d, suffix, out = sys.argv[1:4]
SETS = ["hilli_on_hilli", "hilli_on_stygia", "stygia_on_stygia", "stygia_on_hilli", "all_on_hilli", "all_on_stygia"]
LOW = {"WE_CH_143", "MKA_AUS_CH_211", "MKA_AUS_CH_212", "DU_F_CS_29", "WE_CS_133"}
locs = sorted(pd.read_csv("samples_meta.tsv", sep="\t")["location"].unique())
col = dict(zip(locs, plt.cm.tab10(np.arange(len(locs)))))
mark = {"hilli": "o", "stygia": "^"}
fig, axes = plt.subplots(3, 2, figsize=(12, 15))
for ax, tag in zip(axes.ravel(), SETS):
    try:
        t = pd.read_csv("%s/%s%s" % (d, tag, suffix), sep="\t")
        var = pd.read_csv("%s/%s.var.tsv" % (d, tag), sep="\t")["pct_var"].tolist()
    except FileNotFoundError:
        ax.set_title(tag + " (missing)")
        continue
    for _, r in t.iterrows():
        ax.scatter(r.PC1, r.PC2, color=col[r.location], marker=mark[r.species], s=45,
                   edgecolor="black" if r.country == "NZ" else "none")
    for sp, g in t.groupby("species"):
        z = np.hypot((g.PC1 - g.PC1.mean()) / (g.PC1.std() + 1e-9), (g.PC2 - g.PC2.mean()) / (g.PC2.std() + 1e-9))
        for _, r in g[(z > 2.5) | g["sample"].isin(LOW)].iterrows():
            ax.annotate(r["sample"], (r.PC1, r.PC2), fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax.set_title(tag)
    ax.set_xlabel("PC1 (%.1f%%)" % var[0])
    ax.set_ylabel("PC2 (%.1f%%)" % var[1])
handles = [plt.Line2D([], [], color=col[l], marker="s", ls="") for l in locs]
handles += [plt.Line2D([], [], color="grey", marker=m, ls="") for m in mark.values()]
fig.legend(handles, locs + list(mark), loc="lower center", ncol=5, fontsize=8)
fig.suptitle("%s (black edge = NZ)" % out)
fig.tight_layout(rect=(0, 0.05, 1, 0.97))
fig.savefig(out, dpi=150)
print("wrote", out)
