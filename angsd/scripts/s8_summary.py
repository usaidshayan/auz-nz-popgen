# stage 8 summary: best ngsadmix run per k and set, how many of the 10 seeds reach it,
# per-sample q at the best run, and a bar plot k 2-5 for both assemblies
# clusters matched across k so colours carry over (hungarian on correlation with the k - 1 run)
import glob
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

SETS, KS = ["all_on_hilli", "all_on_stygia"], [1, 2, 3, 4, 5]
names = [l.strip() for l in open("lists/all.txt") if l.strip()]
meta = pd.read_csv("samples_meta.tsv", sep="\t").set_index("sample").loc[names].reset_index()
order = meta.sort_values(["species", "country", "location", "sample"]).index.values

runs, best_q = [], {}
for tag in SETS:
    for K in KS:
        ll = {}
        for f in glob.glob("admix/%s.K%d.s*.log" % (tag, K)):
            m = re.findall(r"best like=(\S+)", open(f).read())
            if m:
                ll[int(re.search(r"\.s(\d+)\.log$", f).group(1))] = float(m[-1])
        if not ll:
            continue
        seed = max(ll, key=ll.get)
        top = max(ll.values())
        runs.append([tag, K, len(ll), top, sum(top - v < 2 for v in ll.values()), seed])
        best_q[(tag, K)] = np.loadtxt("admix/%s.K%d.s%d.qopt" % (tag, K, seed), ndmin=2)
r = pd.DataFrame(runs, columns=["set", "K", "n_runs", "best_loglik", "runs_within_2ll", "best_seed"])
r["delta_ll_from_k_minus_1"] = r.groupby("set").best_loglik.diff()
r.to_csv("tables/S8_ngsadmix_runs.tsv", sep="\t", index=False, float_format="%.2f")
print(r.to_string(index=False))

# align cluster labels across k, and put the hilli cluster first at k = 2
for tag in SETS:
    prev = None
    for K in KS[1:]:
        if (tag, K) not in best_q:
            continue
        q = best_q[(tag, K)]
        if prev is None:
            if q[meta.species.values == "hilli", 0].mean() < 0.5:
                q = q[:, ::-1]
        else:
            c = np.corrcoef(prev.T, q.T)[:prev.shape[1], prev.shape[1]:]
            _, cols = linear_sum_assignment(-np.nan_to_num(c))
            q = q[:, list(cols) + [j for j in range(K) if j not in cols]]
        best_q[(tag, K)] = prev = q

qt = meta[["sample", "species", "country", "location"]].copy()
for (tag, K), q in sorted(best_q.items()):
    if K == 2:
        qt["%s_K2_hilli" % tag] = q[:, 0].round(4)
        qt["%s_K2_stygia" % tag] = q[:, 1].round(4)
qt.to_csv("tables/S8_ngsadmix_K2_q.tsv", sep="\t", index=False)

for f in ["Regular", "Bold", "Italic", "BoldItalic"]:
    fm.fontManager.addfont("/usr/share/fonts/urw-base35/NimbusSans-%s.otf" % f)
plt.rcParams.update({"font.family": "Nimbus Sans", "font.size": 7, "xtick.labelsize": 5, "ytick.labelsize": 6,
                     "axes.linewidth": 0.6, "mathtext.fontset": "custom", "mathtext.rm": "Nimbus Sans",
                     "mathtext.it": "Nimbus Sans:italic", "mathtext.bf": "Nimbus Sans:bold",
                     "mathtext.cal": "Nimbus Sans", "pdf.fonttype": 42, "savefig.dpi": 600})
PAL = ["#DDAA33", "#004488", "#BB5566", "#228833", "#AA4499"]
SP = {"hilli": r"$\it{C.}$ $\it{hilli}$", "stygia": r"$\it{C.}$ $\it{stygia}$"}
CM = 1 / 2.54
ks = [k for k in KS[1:] if all((t, k) in best_q for t in SETS)]
fig, axs = plt.subplots(len(ks), 2, figsize=(17.8 * CM, (2.0 * len(ks) + 2.6) * CM), sharex=True,
                        layout="constrained", squeeze=False)
labels = meta["sample"].values[order]
for j, tag in enumerate(SETS):
    axs[0, j].set_title("Both species, %s reference" % SP[tag.split("_on_")[1]], fontsize=7)
    for i, K in enumerate(ks):
        ax = axs[i, j]
        q = best_q[(tag, K)][order]
        bottom = np.zeros(len(order))
        for c in range(K):
            ax.bar(np.arange(len(order)), q[:, c], bottom=bottom, width=1.0, color=PAL[c], linewidth=0)
            bottom += q[:, c]
        ax.set_xlim(-0.5, len(order) - 0.5)
        ax.set_ylim(0, 1)
        ax.set_yticks([0, 1])
        ax.set_ylabel("K = %d" % K, rotation=0, ha="right", va="center")
        ax.axvline(18 - 0.5, color="white", lw=1.2)
        for s in ["top", "right"]:
            ax.spines[s].set_visible(False)
    axs[-1, j].set_xticks(np.arange(len(order)))
    axs[-1, j].set_xticklabels(labels, rotation=90)
fig.savefig("tables/S8_ngsadmix.pdf", bbox_inches="tight", pad_inches=0.02)
fig.savefig("tables/S8_ngsadmix.png", bbox_inches="tight", pad_inches=0.02)
print("wrote tables/S8_ngsadmix.{pdf,png}, tables/S8_ngsadmix_runs.tsv, tables/S8_ngsadmix_K2_q.tsv")
