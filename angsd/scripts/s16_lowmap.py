# the eight samples mapping < 90% to their own assembly, pulled together for the reply to ang
# figure: pcangsd pc1 vs pc2, both species on each assembly (top) and each species on its own assembly (bottom),
# the eight ringed and labelled
# table: one row per sample with every line of evidence (mapping, pca, ngsadmix, hybrid index, newhybrids, unmapped reads)
import glob

import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

LOW = ["WE_CS_133", "DU_F_CH_26", "WE_CS_139", "MKA_AUS_CH_211", "WE_CH_141", "DU_F_CH_28", "DU_F_CH_24", "HO_AUS_CS_204"]
meta = pd.read_csv("samples_meta.tsv", sep="\t").set_index("sample")

# table
t1 = pd.read_csv("tables/T1_mapping_depth.tsv", sep="\t").set_index("sample")
s2 = pd.read_csv("tables/S2_competitive_mapping.tsv", sep="\t").set_index("sample")
s3 = pd.read_csv("tables/S3_unmapped_kraken.tsv", sep="\t").set_index("sample")
q = pd.read_csv("tables/S8_ngsadmix_K2_q.tsv", sep="\t").set_index("sample")
nh = pd.read_csv("tables/S14_newhybrids.tsv", sep="\t")
d = meta[["species", "country", "location"]].copy()
d["depth_own"] = t1.effDP_own
d["mapped_own_pct"] = t1.mapped_own
d["mapped_other_pct"] = t1.mapped_other
d["competitive_own_share_pct"] = s2.own_share
d["competitive_other_species_pct"] = np.where(d.species == "hilli", s2.placed_stygia.loc[d.index], s2.placed_hilli.loc[d.index])
d["unmapped_pct"] = s3.pct_unmapped
d["unmapped_bacteria_pct_of_all_reads"] = (s3.bacteria * s3.pct_unmapped / 100).round(2)
d["unmapped_top_genera"] = s3.top_genera
for tag in ["all_on_hilli", "all_on_stygia"]:
    r = tag.split("_on_")[1]
    pcs = pd.read_csv("pcangsd/%s.pcs.tsv" % tag, sep="\t").set_index("sample")
    # position along pc1 between the two species means, 0 = own species mean, 1 = other species mean
    m = pcs.groupby("species").PC1.mean()
    own, oth = m.loc[d.species].values, m.loc[np.where(d.species == "hilli", "stygia", "hilli")].values
    d["pc1_shift_to_other_%s_ref" % r] = ((pcs.PC1.loc[d.index].values - own) / (oth - own)).round(3)
    d["ngsadmix_K2_own_%s_ref" % r] = [q.loc[s, "%s_K2_%s" % (tag, d.species[s])] for s in d.index]
    h = pd.read_csv("hybrid/%s.hybrid_index.tsv" % tag, sep="\t").set_index("sample")
    d["other_species_ancestry_%s_ref" % r] = h.other_ancestry
    d["interspecific_het_%s_ref" % r] = h.q12
    for n in [200, 500]:
        x = nh[nh.set == "%s_%d" % (tag, n)].set_index("sample")
        d["newhybrids_%d_min_p_pure_own_%s_ref" % (n, r)] = x.p_pure_own_min
d.insert(3, "low_mapping", d.index.isin(LOW))
d = d.sort_values(["low_mapping", "species", "mapped_own_pct"], ascending=[False, True, True])
d.to_csv("tables/S16_per_sample_evidence.tsv", sep="\t", float_format="%.4f")
print(d[d.low_mapping].T.to_string())

# figure
for f in ["Regular", "Bold", "Italic", "BoldItalic"]:
    fm.fontManager.addfont("/usr/share/fonts/urw-base35/NimbusSans-%s.otf" % f)
plt.rcParams.update({
    "font.family": "Nimbus Sans", "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6, "axes.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "mathtext.fontset": "custom",
    "mathtext.rm": "Nimbus Sans", "mathtext.it": "Nimbus Sans:italic", "mathtext.bf": "Nimbus Sans:bold",
    "mathtext.cal": "Nimbus Sans", "mathtext.sf": "Nimbus Sans", "pdf.fonttype": 42, "savefig.dpi": 600,
})
CM = 1 / 2.54
COL = {"hilli": "#DDAA33", "stygia": "#004488"}
MARK = {"hilli": "o", "stygia": "^"}
SP = {"hilli": r"$\it{C.}$ $\it{hilli}$", "stygia": r"$\it{C.}$ $\it{stygia}$"}
fig, axs = plt.subplots(2, 2, figsize=(17.8 * CM, 15.5 * CM), layout="constrained")
for ax, tag, letter in zip(axs.flat, ["all_on_hilli", "all_on_stygia", "hilli_on_hilli", "stygia_on_stygia"], "ABCD"):
    t = pd.read_csv("pcangsd/%s.pcs.tsv" % tag, sep="\t")
    var = pd.read_csv("pcangsd/%s.var.tsv" % tag, sep="\t")["pct_var"].tolist()
    for _, r in t.iterrows():
        low = r["sample"] in LOW
        ax.scatter(r.PC1, r.PC2, marker=MARK[r.species], s=22, c=COL[r.species], alpha=0.9, zorder=3,
                   edgecolors="black" if low else "0.15", linewidths=1.0 if low else 0.4)
    # numbers instead of names, samples sitting on top of each other share one label
    lo = t[t["sample"].isin(LOW)].copy()
    lo["n"] = [LOW.index(x) + 1 for x in lo["sample"]]
    span = np.array([np.ptp(t.PC1), np.ptp(t.PC2)])
    done = set()
    for k, r in lo.iterrows():
        if k in done:
            continue
        near = lo[(np.abs(lo.PC1 - r.PC1) < 0.03 * span[0]) & (np.abs(lo.PC2 - r.PC2) < 0.03 * span[1])]
        done |= set(near.index)
        ax.annotate(",".join(str(x) for x in sorted(near.n)), (near.PC1.mean(), near.PC2.mean()), xytext=(5, 3),
                    textcoords="offset points", fontsize=6, zorder=4)
    ax.set_xlabel("PC1 (%.1f%%)" % var[0])
    ax.set_ylabel("PC2 (%.1f%%)" % var[1])
    ax.margins(0.12)
    who, ref = tag.split("_on_")
    ax.set_title("%s, %s reference" % ("Both species" if who == "all" else SP[who], SP[ref]))
    ax.text(-0.14, 1.03, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")
h = [Line2D([], [], ls="", marker=MARK[s], ms=4.5, mfc=COL[s], mec="0.15", mew=0.4) for s in COL]
h.append(Line2D([], [], ls="", marker="o", ms=4.5, mfc="white", mec="black", mew=1.0))
key = ["%d  %s" % (i + 1, x) for i, x in enumerate(LOW)]
blank = Line2D([], [], ls="")
first = list(zip(h, [SP["hilli"], SP["stygia"], "<90% of reads mapped to own assembly"])) + [(blank, "")]
items = [x for c in range(4) for x in [first[c], (blank, key[c]), (blank, key[c + 4])]]
fig.legend(*zip(*items), loc="outside lower center", ncol=4, frameon=False, handletextpad=0.3, columnspacing=1.5,
           handlelength=1.0)
for ext in ["pdf", "png"]:
    fig.savefig("tables/S16_low_mapping_pca.%s" % ext, bbox_inches="tight", pad_inches=0.02)
print("wrote tables/S16_per_sample_evidence.tsv, tables/S16_low_mapping_pca.{pdf,png}")
