# scree for the single-species pcangsd runs: % variance per pc against the no-structure level (1 / (n - 1))
# each bar named by the fly carrying most of that axis when one fly carries >= 50% of it
import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np

for f in ["Regular", "Bold", "Italic", "BoldItalic"]:
    fm.fontManager.addfont("/usr/share/fonts/urw-base35/NimbusSans-%s.otf" % f)
plt.rcParams.update({"font.family": "Nimbus Sans", "font.size": 7, "axes.linewidth": 0.6, "axes.spines.top": False,
                     "axes.spines.right": False, "mathtext.fontset": "custom", "mathtext.it": "Nimbus Sans:italic",
                     "mathtext.rm": "Nimbus Sans", "pdf.fonttype": 42, "savefig.dpi": 600})
CM = 1 / 2.54
fig, axs = plt.subplots(1, 2, figsize=(17.8 * CM, 6.5 * CM), layout="constrained")
for ax, (tag, sp, col) in zip(axs, [("hilli_on_hilli", "hilli", "#DDAA33"), ("stygia_on_stygia", "stygia", "#004488")]):
    names = [l.strip() for l in open("lists/%s.txt" % sp) if l.strip()]
    v, w = np.linalg.eigh(np.loadtxt("pcangsd/%s.cov" % tag))
    o = v.argsort()[::-1]
    v, w = np.clip(v[o], 0, None), w[:, o]
    var = 100 * v / v.sum()
    k = 10
    ax.bar(range(1, k + 1), var[:k], color=col, width=0.7)
    ax.axhline(100 / (len(names) - 1), color="0.3", lw=0.6, ls="--", label="no structure, 1 / (n - 1)")
    ax.legend(loc="upper right", frameon=False, fontsize=6)
    for i in range(k):
        s = w[:, i] ** 2
        if s.max() >= 0.5:
            ax.text(i + 1, var[i] + 0.4, "%s\n%.0f%%" % (names[s.argmax()], 100 * s.max()), rotation=90,
                    ha="center", va="bottom", fontsize=4.5)
    ax.set_xticks(range(1, k + 1))
    ax.set_xlabel("PC")
    ax.set_ylabel("Variance explained (%)")
    ax.set_title(r"$\it{C.\ %s}$, own reference" % sp)
fig.savefig("figures/pca/pcangsd/scree_single_species.pdf", bbox_inches="tight", pad_inches=0.02)
fig.savefig("figures/pca/pcangsd/scree_single_species.png", bbox_inches="tight", pad_inches=0.02)
