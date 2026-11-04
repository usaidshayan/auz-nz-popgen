# stage 4 (ang t2): pca on hard-called snps from the all-sites zarr, one sample subset
# site filters within the subset: biallelic snp, called in >= 72% (mask v1, 13 of 18),
# depth <= 3x median (collapsed repeats), mac >= 2 (singletons mostly errors at this depth)
# thinned at random to ~300k snps (~1 per 1.5 kb) instead of ld pruning
# patterson scaled, missing set to 0. low call rate pulls a sample to the origin so call rate
# is written next to the coordinates
# usage: python s4_pca_hard.py <store.zarr> <sample_list> <out_prefix>
import sys

import numpy as np
import pandas as pd
import zarr

store, lst, out = sys.argv[1:4]
TARGET = 300000
rng = np.random.default_rng(1)
z = zarr.open_group(store, mode="r")
all_names = [str(s) for s in z["samples"][:]]
names = [l.strip() for l in open(lst) if l.strip()]
idx = np.array([all_names.index(s) for s in names])
n = len(idx)
need = int(np.ceil(0.72 * n))
gt_arr, dp_arr = z["calldata/GT"], z["calldata/DP"]
N, step = gt_arr.shape[0], gt_arr.chunks[0]
starts = list(range(0, N, step))


def load(s):
    gt = gt_arr[s:s + step][:, idx, :]
    dp = np.clip(dp_arr[s:s + step][:, idx].astype(np.int32), 0, None).sum(axis=1)
    called = (gt[:, :, 0] >= 0) & (gt[:, :, 1] >= 0)
    a0, a1, ax = (gt == 0).sum((1, 2)), (gt == 1).sum((1, 2)), (gt >= 2).sum((1, 2))
    snp = (ax == 0) & (np.minimum(a0, a1) >= 2) & (called.sum(axis=1) >= need)
    return gt, called, dp, snp


# pass 1 (every 20th chunk): depth cap and expected no. of passing snps
dps, snps = [], []
for s in starts[::20]:
    _, _, dp, snp = load(s)
    dps.append(dp)
    snps.append(snp)
dps, snps = np.concatenate(dps), np.concatenate(snps)
cap = 3 * np.median(dps[dps > 0])
expected = (snps & (dps <= cap)).sum() * len(starts) / len(starts[::20])
p_keep = min(1.0, TARGET / max(expected, 1))
print("%s | %s | %d samples | %d sites in store" % (out, store, n, N))
print("filters: called in >= %d of %d, total depth <= %.0f (3x median), MAC >= 2, biallelic" % (need, n, cap))
print("expected passing SNPs: %.2f million; keeping each with probability %.4f" % (expected / 1e6, p_keep))

# pass 2 (all chunks): alt allele counts for kept snps
X, n_pass = [], 0
for s in starts:
    gt, called, dp, snp = load(s)
    ok = snp & (dp <= cap)
    n_pass += ok.sum()
    ok &= rng.random(len(ok)) < p_keep
    g = (gt[ok] == 1).sum(axis=2).astype(np.float32)
    g[~called[ok]] = np.nan
    X.append(g)
X = np.concatenate(X)
print("SNPs passing: %d | used for PCA: %d" % (n_pass, X.shape[0]))

call_rate = 100 * (~np.isnan(X)).mean(axis=0)
p = np.nanmean(X, axis=1) / 2
Xs = (X - 2 * p[:, None]) / np.sqrt(2 * p * (1 - p))[:, None]
Xs[~np.isfinite(Xs)] = 0
U, S, _ = np.linalg.svd(Xs.T, full_matrices=False)
k = min(10, n)
pcs = pd.DataFrame(U[:, :k] * S[:k], columns=["PC%d" % (i + 1) for i in range(k)])
pcs.insert(0, "call_rate", call_rate.round(1))
pcs.insert(0, "sample", names)
meta = pd.read_csv("samples_meta.tsv", sep="\t")
pcs = meta.merge(pcs, on="sample")
var = 100 * S ** 2 / (S ** 2).sum()
pcs.round(4).to_csv(out + ".pcs.tsv", sep="\t", index=False)
pd.DataFrame({"PC": range(1, len(var) + 1), "pct_var": var.round(3)}).to_csv(out + ".var.tsv", sep="\t", index=False)
print("variance explained PC1-PC5 (%%): %s" % ", ".join("%.1f" % v for v in var[:5]))
pd.set_option("display.width", 250)
print(pcs[["sample", "species", "country", "location", "call_rate", "PC1", "PC2", "PC3"]].sort_values("PC1").round(3).to_string(index=False))
