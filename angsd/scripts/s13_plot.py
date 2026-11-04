# triangle plot, hybrid index (share of c. stygia alleles) vs interspecific het, one panel per assembly
# expected classes drawn as open squares, the eight low-mapping samples ringed
# usage: python s13_plot.py <out prefix>
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

out = sys.argv[1]
for f in ["Regular", "Bold", "Italic", "BoldItalic"]:
    fm.fontManager.addfont("/usr/share/fonts/urw-base35/NimbusSans-%s.otf" % f)
plt.rcParams.update({
    "font.family": "Nimbus Sans", "font.size": 7, "axes.labelsize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6,
    "legend.fontsize": 6, "axes.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
    "mathtext.fontset": "custom", "mathtext.rm": "Nimbus Sans", "mathtext.it": "Nimbus Sans:italic",
    "mathtext.bf": "Nimbus Sans:bold", "mathtext.cal": "Nimbus Sans", "mathtext.sf": "Nimbus Sans",
    "pdf.fonttype": 42, "savefig.dpi": 600,
})
CM = 1 / 2.54
LOW = {"WE_CS_133", "DU_F_CH_26", "WE_CS_139", "MKA_AUS_CH_211", "WE_CH_141", "DU_F_CH_28", "DU_F_CH_24", "HO_AUS_CS_204"}
COL = {"hilli": "#DDAA33", "stygia": "#004488"}
MARK = {"hilli": "o", "stygia": "^"}
SP = {"hilli": r"$\it{C.}$ $\it{hilli}$", "stygia": r"$\it{C.}$ $\it{stygia}$"}
CLASSES = [(0, 0, ""), (1, 0, ""), (0.5, 1, "F1"), (0.5, 0.5, "F2"), (0.25, 0.5, "BC"), (0.75, 0.5, "BC")]

fig, axs = plt.subplots(1, 2, figsize=(17.8 * CM, 8.2 * CM), layout="constrained")
for ax, tag, letter in zip(axs, ["all_on_hilli", "all_on_stygia"], "AB"):
    t = pd.read_csv("hybrid/%s.hybrid_index.tsv" % tag, sep="\t")
    ax.plot([0, 0.5, 1], [0, 1, 0], color="0.6", lw=0.6, zorder=1)
    for x, y, lab in CLASSES:
        ax.scatter(x, y, marker="s", s=30, facecolors="none", edgecolors="0.45", linewidths=0.6, zorder=2)
        if lab:
            ax.annotate(lab, (x, y), xytext=(5, 3), textcoords="offset points", fontsize=6, color="0.35")
    for _, r in t.iterrows():
        ax.scatter(r.stygia_fraction, r.q12, marker=MARK[r.species], s=22, c=COL[r.species],
                   edgecolors="black" if r["sample"] in LOW else "0.15", linewidths=1.0 if r["sample"] in LOW else 0.4,
                   alpha=0.9, zorder=3)
    ax.set_xlim(-0.04, 1.04)
    ax.set_ylim(-0.04, 1.06)
    ax.set_xlabel(r"Hybrid index (proportion of $\it{C.}$ $\it{stygia}$ alleles)")
    ax.set_ylabel("Interspecific heterozygosity")
    ax.set_title("%s reference" % SP[tag.split("_on_")[1]], fontsize=7)
    ax.text(-0.16, 1.03, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")
h = [Line2D([], [], ls="", marker=MARK[s], ms=4.5, mfc=COL[s], mec="0.15", mew=0.4) for s in COL]
h += [Line2D([], [], ls="", marker="o", ms=4.5, mfc="white", mec="black", mew=1.0),
      Line2D([], [], ls="", marker="s", ms=4.5, mfc="none", mec="0.45", mew=0.6)]
fig.legend(h, [SP["hilli"], SP["stygia"], "<90% mapped to own assembly", "expected class"],
           loc="outside lower center", ncol=4, frameon=False, handletextpad=0.3, columnspacing=1.2)
for ext in ["pdf", "png"]:
    fig.savefig("%s.%s" % (out, ext), bbox_inches="tight", pad_inches=0.02)
print("wrote %s.{pdf,png}" % out)
