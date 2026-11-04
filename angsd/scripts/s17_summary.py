# s17 checks on the gl pca
# 1. share of snp sites with reads, per sample, next to depth
# 2. split halves (alternating 1 mb blocks) against each other and the full run
# 3. downsampled WE_CH_131 / BL_AUS_CH_214 against themselves at full depth, same sites
# position = distance to own location's centroid (without the sample) on pc1-pc2,
# scaled by the distance between the aus and nz centroids
import numpy as np
import pandas as pd

dp = pd.read_csv("tables/T1_mapping_depth.tsv", sep="\t").set_index("sample")
rows = []


def pcs(f):
    return pd.read_csv(f, sep="\t").set_index("sample")


def offset(t, s):
    loc = t.location[s]
    own = t[(t.location == loc) & (t.index != s)][["PC1", "PC2"]].mean()
    au, nz = t[t.country == "AUS"][["PC1", "PC2"]].mean(), t[t.country == "NZ"][["PC1", "PC2"]].mean()
    if t.species.iloc[0] == "hilli":
        au = t[(t.country == "AUS") & (t.location != "MtKeiraNSW_AUS")][["PC1", "PC2"]].mean()
    return np.hypot(*(t.loc[s, ["PC1", "PC2"]] - own)) / np.hypot(*(au - nz))


print("1. share of snp sites where each sample has reads")
for tag, lst in [("hilli_on_hilli", "lists/hilli.txt"), ("stygia_on_stygia", "lists/stygia.txt")]:
    names = [l.strip() for l in open(lst) if l.strip()]
    x = pd.read_csv("robust/%s.informative.tsv" % tag, sep="\t", header=None, names=["set", "i", "with_reads", "sites"])
    x["sample"] = [names[i] for i in x.i]
    x["pct_sites_with_reads"] = 100 * x.with_reads / x.sites
    x["effDP"] = dp.effDP_own.loc[x["sample"]].values
    x = x.sort_values("pct_sites_with_reads")
    print(x[["set", "sample", "effDP", "with_reads", "sites", "pct_sites_with_reads"]].round(2).to_string(index=False))
    x[["set", "sample", "effDP", "with_reads", "sites", "pct_sites_with_reads"]].round(3).to_csv(
        "tables/S17_sites_with_reads_%s.tsv" % tag, sep="\t", index=False)

print("\n2. split halves")
for tag in ["hilli_on_hilli", "stygia_on_stygia"]:
    full = pcs("pcangsd/%s.pcs.tsv" % tag)
    a, b = pcs("robust/half/%s.A.pcs.tsv" % tag), pcs("robust/half/%s.B.pcs.tsv" % tag)
    for pc in ["PC1", "PC2", "PC3"]:
        r_ab = abs(np.corrcoef(a[pc], b.loc[a.index, pc])[0, 1])
        r_af = abs(np.corrcoef(a[pc], full.loc[a.index, pc])[0, 1])
        r_bf = abs(np.corrcoef(b[pc], full.loc[b.index, pc])[0, 1])
        rows.append([tag, pc, r_ab, r_af, r_bf])
        print("%s %s: |r| half A vs B %.3f, A vs full %.3f, B vs full %.3f" % (tag, pc, r_ab, r_af, r_bf))
    for s in (["WE_CH_143"] if tag.startswith("hilli") else ["DU_F_CS_29", "WE_CS_133"]):
        print("  %s offset from own location: full %.3f, half A %.3f, half B %.3f"
              % (s, offset(full, s), offset(a, s), offset(b, s)))

pd.DataFrame(rows, columns=["set", "pc", "r_half_A_vs_B", "r_half_A_vs_full", "r_half_B_vs_full"]).round(4).to_csv(
    "tables/S17_split_half.tsv", sep="\t", index=False)
rows = []

print("\n3. downsampled to WE_CH_143's depth, same sites")
base = pcs("robust/gl/base.pcs.tsv")
print("baseline offsets (all hilli, same sites):")
off = pd.Series({s: offset(base, s) for s in base.index}).sort_values()
print(off.round(3).to_string())
for s in ["WE_CH_131", "BL_AUS_CH_214"]:
    ds = pcs("robust/gl/ds_%s.pcs.tsv" % s)
    sgn = np.sign(np.corrcoef(base.PC1, ds.loc[base.index, "PC1"])[0, 1]), \
        np.sign(np.corrcoef(base.PC2, ds.loc[base.index, "PC2"])[0, 1])
    print("%s: full depth offset %.3f -> downsampled %.3f (WE_CH_143 %.3f); PC1/PC2 %.2f,%.2f -> %.2f,%.2f"
          % (s, offset(base, s), offset(ds, s), offset(base, "WE_CH_143"), base.PC1[s], base.PC2[s],
             sgn[0] * ds.PC1[s], sgn[1] * ds.PC2[s]))
    rows.append([s, base.PC1[s], sgn[0] * ds.PC1[s], base.PC2[s], sgn[1] * ds.PC2[s], sgn[1] * ds.PC2[s] / base.PC2[s],
                 offset(base, s), offset(ds, s), offset(base, "WE_CH_143")])
nz = base[(base.country == "NZ") & (base.index != "WE_CH_143")].PC2.mean()
print("WE_CH_143 PC2 %.3f vs other nz hilli mean %.3f (ratio %.2f)" % (base.PC2["WE_CH_143"], nz, base.PC2["WE_CH_143"] / nz))
pd.DataFrame(rows, columns=["sample", "PC1_full", "PC1_downsampled", "PC2_full", "PC2_downsampled", "PC2_ratio",
                            "offset_full_depth", "offset_downsampled", "offset_WE_CH_143"]).round(4).to_csv(
    "tables/S17_downsampled_pca.tsv", sep="\t", index=False)
