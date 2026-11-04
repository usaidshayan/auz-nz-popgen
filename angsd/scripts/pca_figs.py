# pca figures for the six sets, pnas sizing (8.7 cm single, 17.8 cm double), nimbus sans
# per set: pc1 vs pc2 on its own, and pc1 vs pc2 + pc2 vs pc3 side by side
# colour = location (warm = aus, cool = nz), marker = species, every sample numbered with a key
# usage: python pca_figs.py <pc dir> <out dir> [pc pair for panel b, default 2,3]
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.text import OffsetFrom

d, out = sys.argv[1:3]
pb = tuple(int(x) for x in (sys.argv[3] if len(sys.argv) > 3 else "2,3").split(","))
os.makedirs(out, exist_ok=True)

for f in ["Regular", "Bold", "Italic", "BoldItalic"]:
    fm.fontManager.addfont("/usr/share/fonts/urw-base35/NimbusSans-%s.otf" % f)
plt.rcParams.update({
    "font.family": "Nimbus Sans", "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "mathtext.fontset": "custom", "mathtext.rm": "Nimbus Sans", "mathtext.it": "Nimbus Sans:italic",
    "mathtext.bf": "Nimbus Sans:bold", "mathtext.cal": "Nimbus Sans", "mathtext.sf": "Nimbus Sans",
    "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 600,
})
CM = 1 / 2.54

LOCS = [  # north to south within country
    ("BlackheathNSW_AUS", "Blackheath, NSW", "#882255"),
    ("MtKeiraNSW_AUS", "Mt Keira, NSW", "#D55E00"),
    ("Canberra_AUS", "Canberra, ACT", "#E69F00"),
    ("SeafordVIC_AUS", "Seaford, VIC", "#CC79A7"),
    ("HobartTAS_AUS", "Hobart, TAS", "#DDCC77"),
    ("PirongiaKaniwhaniwha", "Pirongia", "#009E73"),
    ("Wellington", "Wellington", "#56B4E9"),
    ("DunedinFairfield", "Dunedin", "#0072B2"),
]
NAME = {k: n for k, n, _ in LOCS}
COL = {k: c for k, _, c in LOCS}
MARK = {"hilli": "o", "stygia": "^"}
SP = {"hilli": r"$\it{C.}$ $\it{hilli}$", "stygia": r"$\it{C.}$ $\it{stygia}$"}
SETS = ["hilli_on_hilli", "stygia_on_stygia", "hilli_on_stygia", "stygia_on_hilli", "all_on_hilli", "all_on_stygia"]


def title(tag, n):
    who, ref = tag.split("_on_")
    who = "Both species" if who == "all" else SP[who]
    return "%s (n = %d), %s reference" % (who, n, SP[ref])


def scatter(ax, t, var, i, j):
    for _, r in t.iterrows():
        ax.scatter(r["PC%d" % i], r["PC%d" % j], s=22 if r.species == "hilli" else 26,
                   marker=MARK[r.species], c=COL[r.location], edgecolors="0.15",
                   linewidths=0.4, alpha=0.9, zorder=3)
    ax.set_xlabel("PC%d (%.1f%%)" % (i, var[i - 1]))
    ax.set_ylabel("PC%d (%.1f%%)" % (j, var[j - 1]))
    ax.margins(0.12)
    ax.locator_params(nbins=5)


def label(ax, x, y, nums, gap=0.04, fs=5, rows=6):
    # each point gets its sample number (key under the plot). points close together (single linkage on
    # axis-scaled distance) share a block of numbers beside the group, up to six per column, with leader lines.
    # each block goes on whichever side (then nudged up or down) clears the points and blocks already placed.
    # call after the layout is drawn
    tr, px = ax.transData, ax.figure.dpi / 72
    p = np.c_[(x - x.min()) / np.ptp(x), (y - y.min()) / np.ptp(y)]
    grp = np.arange(len(p))
    for a in range(len(p)):
        for b in range(len(p)):
            if np.hypot(*(p[a] - p[b])) < gap:
                grp[grp == grp[b]] = grp[a]
    d = tr.transform(np.c_[x, y])
    lh, cw, pad, r = 1.35 * fs * px, 1.6 * fs * px, 7 * px, 4 * px
    box = ax.get_window_extent()
    placed = [(u - r, v - r, u + r, v + r) for u, v in d]
    hit = lambda q: sum(max(0, min(q[2], o[2]) - max(q[0], o[0])) * max(0, min(q[3], o[3]) - max(q[1], o[1]))
                        for o in placed)
    out = lambda q: max(0, box.x0 - q[0]) + max(0, q[2] - box.x1) + max(0, box.y0 - q[1]) + max(0, q[3] - box.y1)
    for g in sorted(np.unique(grp), key=lambda g: -(grp == g).sum()):
        k = np.where(grp == g)[0]
        k = k[np.argsort(-d[k, 1])]
        nc = -(-len(k) // rows)
        nr = -(-len(k) // nc)
        w, h = cw * nc, lh * nr
        cy = d[k, 1].mean()
        best = None
        for shift in [0, 1, -1, 2, -2, 3, -3, 4, -4]:
            for side, far in [(1, 1), (-1, 1), (1, 3), (-1, 3), (1, 6), (-1, 6)]:
                x0 = (d[k, 0].max() + far * pad) if side > 0 else (d[k, 0].min() - far * pad - w)
                y0 = cy - h / 2 + shift * lh * 1.5
                q = (x0, y0, x0 + w, y0 + h)
                cost = hit(q) + 50 * out(q) * h + abs(shift) * 0.5 + far * 0.5
                if best is None or cost < best[0]:
                    best = (cost, side, q)
        _, side, q = best
        placed.append(q)
        for n, s in enumerate(k):
            c = n // nr if side > 0 else nc - 1 - n // nr
            tx, ty = q[0] + (c + 0.5) * cw, q[3] - (n % nr + 0.5) * lh
            ax.annotate(str(nums[s]), (x[s], y[s]), xytext=tuple(tr.inverted().transform((tx, ty))),
                        textcoords="data", ha="center", va="center", fontsize=fs, color="0.1", zorder=4,
                        arrowprops=dict(arrowstyle="-", color="0.55", lw=0.3, shrinkA=0, shrinkB=2))


def legend(fig, t, ncol):
    # key: one block per species x location, bold location name then its samples by number,
    # marker = species, colour = location. blocks dealt across the columns in turn
    blank = lambda: Line2D([], [], ls="")
    blocks = []
    for sp in sorted(t.species.unique()):
        for k, _, _ in LOCS:
            u = t[(t.species == sp) & (t.location == k)].sort_values("n")
            if len(u) == 0:
                continue
            g = [(blank(), r"$\bf{%s}$" % NAME[k].replace(" ", r"\ "))]
            for _, r in u.iterrows():
                g.append((Line2D([], [], ls="", marker=MARK[sp], ms=4, mfc=COL[k], mec="0.15", mew=0.4),
                          "%d  %s" % (r.n, r["sample"])))
            blocks.append(g)
    cols = [[x for b in blocks[c::ncol] for x in b] for c in range(ncol)]
    rows = max(len(c) for c in cols)
    items = [x for c in cols for x in c + [(blank(), "")] * (rows - len(c))]
    fig.legend(*zip(*items), loc="outside lower center", ncol=ncol, frameon=False, fontsize=5,
               handletextpad=0.3, columnspacing=1.0, labelspacing=0.25, handlelength=0.8)


def save(fig, name):
    for ext in ["pdf", "png"]:
        fig.savefig("%s/%s.%s" % (out, name, ext), bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print("wrote", "%s/%s.{pdf,png}" % (out, name))


# sample numbers fixed across every figure: species, country, location (as in LOCS), sample
meta = pd.read_csv("samples_meta.tsv", sep="\t")
meta["o"] = meta.location.map({k: i for i, (k, _, _) in enumerate(LOCS)})
NUM = dict(zip(meta.sort_values(["species", "o", "sample"])["sample"], range(1, len(meta) + 1)))

for tag in SETS:
    t = pd.read_csv("%s/%s.pcs.tsv" % (d, tag), sep="\t")
    t["n"] = t["sample"].map(NUM)
    var = pd.read_csv("%s/%s.var.tsv" % (d, tag), sep="\t")["pct_var"].tolist()
    both = t.species.nunique() > 1

    fig = plt.figure(figsize=(8.7 * CM, 9.6 * CM if both else 9.0 * CM), layout="constrained")
    ax = fig.add_subplot()
    scatter(ax, t, var, 1, 2)
    ax.set_title(title(tag, len(t)), pad=4)
    legend(fig, t, 3)
    fig.canvas.draw()
    label(ax, t.PC1.values, t.PC2.values, t.n.values)
    save(fig, "%s_pc1_pc2" % tag)

    fig, axs = plt.subplots(1, 2, figsize=(17.8 * CM, 8.4 * CM), layout="constrained")
    for ax, (i, j), letter in zip(axs, [(1, 2), pb], "AB"):
        scatter(ax, t, var, i, j)
        ax.text(-0.2, 1.04, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")
    fig.suptitle(title(tag, len(t)), fontsize=7.5)
    legend(fig, t, 6)
    fig.canvas.draw()
    for ax, (i, j) in zip(axs, [(1, 2), pb]):
        label(ax, t["PC%d" % i].values, t["PC%d" % j].values, t.n.values)
    save(fig, "%s_pc1_pc2__pc%d_pc%d" % (tag, *pb))
