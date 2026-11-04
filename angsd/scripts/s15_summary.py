# evaladmix residual correlation heatmaps, k 2-5 x both assemblies, samples ordered as in the ngsadmix plot
# also the mean residual within and between groups (species x country) per k
import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SETS, KS = ["all_on_hilli", "all_on_stygia"], [2, 3, 4, 5]
names = [l.strip() for l in open("lists/all.txt") if l.strip()]
meta = pd.read_csv("samples_meta.tsv", sep="\t").set_index("sample").loc[names].reset_index()
order = meta.sort_values(["species", "country", "location", "sample"]).index.values
grp = (meta.species + "_" + meta.country).values

for f in ["Regular", "Bold", "Italic", "BoldItalic"]:
    fm.fontManager.addfont("/usr/share/fonts/urw-base35/NimbusSans-%s.otf" % f)
plt.rcParams.update({"font.family": "Nimbus Sans", "font.size": 7, "xtick.labelsize": 4, "ytick.labelsize": 4,
                     "axes.linewidth": 0.6, "mathtext.fontset": "custom", "mathtext.rm": "Nimbus Sans",
                     "mathtext.it": "Nimbus Sans:italic", "mathtext.bf": "Nimbus Sans:bold",
                     "mathtext.cal": "Nimbus Sans", "pdf.fonttype": 42, "savefig.dpi": 600})
SP = {"hilli": r"$\it{C.}$ $\it{hilli}$", "stygia": r"$\it{C.}$ $\it{stygia}$"}
CM = 1 / 2.54
fig, axs = plt.subplots(2, 4, figsize=(17.8 * CM, 9.0 * CM), layout="constrained")
GL = [r"$\it{C.\ hilli}$ AUS", r"$\it{C.\ hilli}$ NZ", r"$\it{C.\ stygia}$ AUS", r"$\it{C.\ stygia}$ NZ"]
rows = []
for j, tag in enumerate(SETS):
    for i, K in enumerate(KS):
        c = np.genfromtxt("evaladmix/%s.K%d.corres.txt" % (tag, K), missing_values="NA", filling_values=np.nan)
        off = ~np.eye(len(c), dtype=bool)
        for g in sorted(set(grp)):
            w = np.outer(grp == g, grp == g) & off
            rows.append([tag, K, g, np.nanmean(c[w])])
        rows.append([tag, K, "all pairs", np.nanmean(c[off])])
        rows.append([tag, K, "max |r|", np.nanmax(np.abs(c[off]))])
        ax = axs[j, i]
        cc = c[np.ix_(order, order)].copy()
        cc[np.tril_indices(len(cc))] = np.nan
        im = ax.imshow(cc, cmap="RdBu_r", vmin=-0.1, vmax=0.1, interpolation="nearest")
        for b in [9, 18, 27]:
            ax.axhline(b - 0.5, color="0.3", lw=0.4)
            ax.axvline(b - 0.5, color="0.3", lw=0.4)
        ax.set_xticks([])
        ax.set_yticks([4, 13, 22, 31])
        ax.set_yticklabels(GL if i == 0 else [], fontsize=6)
        ax.tick_params(length=0)
        ax.set_title("K = %d" % K, fontsize=7)
        if i == 0:
            ax.text(-0.62, 0.5, "%s reference" % SP[tag.split("_on_")[1]], transform=ax.transAxes, rotation=90,
                    va="center", ha="center", fontsize=7)
fig.colorbar(im, ax=axs, shrink=0.5, label="Correlation of residuals")
fig.savefig("tables/S15_evaladmix.pdf", bbox_inches="tight", pad_inches=0.02)
fig.savefig("tables/S15_evaladmix.png", bbox_inches="tight", pad_inches=0.02)
t = pd.DataFrame(rows, columns=["set", "K", "group", "mean_residual_corr"])
t = t.pivot_table(index=["set", "K"], columns="group", values="mean_residual_corr").round(4)
t.to_csv("tables/S15_evaladmix_summary.tsv", sep="\t")
print(t.to_string())
