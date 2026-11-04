# usage: python scripts/site_qc.py <store.zarr>
# site-level depth and coverage summary, and how many sites survive some candidate masks
# (every 20th block of sites along the genome, scaled up to the whole store)
import sys

import allel
import numpy as np
import zarr

z = zarr.open_group(sys.argv[1], mode="r")
gt_arr, dp_arr = z["calldata/GT"], z["calldata/DP"]
n_sites, k = gt_arr.shape[0], gt_arr.shape[1]
step = gt_arr.chunks[0]
tot_dp, n_called, seg = [], [], []
for s in list(range(0, n_sites, step))[::20]:
    gt = gt_arr[s:s + step]
    dp = np.clip(dp_arr[s:s + step].astype(np.int32), 0, None)
    tot_dp.append(dp.sum(axis=1))
    n_called.append(((gt[:, :, 0] >= 0) & (gt[:, :, 1] >= 0)).sum(axis=1))
    seg.append(allel.GenotypeArray(gt).count_alleles().is_segregating())
tot_dp, n_called, seg = np.concatenate(tot_dp), np.concatenate(n_called), np.concatenate(seg)
n = len(tot_dp)
scale = n_sites / n
med = np.median(tot_dp[tot_dp > 0])
print("%d samples, %d sites in store, %d sampled (x%.1f to scale up)" % (k, n_sites, n, scale))
print("total depth per site across samples, percentiles 1/5/25/50/75/95/99:",
      np.percentile(tot_dp, [1, 5, 25, 50, 75, 95, 99]).round(1).tolist())
print("median total depth over sites with any reads: %.1f (about %.1f per sample)" % (med, med / k))
print("samples with a called genotype per site, percentiles 5/25/50/75/95:",
      np.percentile(n_called, [5, 25, 50, 75, 95]).astype(int).tolist())
print()
print("%-34s %12s %14s %18s" % ("candidate mask", "sites kept", "kept (millions)", "segregating (millions)"))
for frac in (0.5, 0.67, 0.83, 1.0):
    need = int(np.ceil(frac * k))
    for mult in (2, 3):
        keep = (n_called >= need) & (tot_dp <= mult * med)
        print("%-34s %11.1f%% %14.1f %18.2f" % (
            "called in >=%d samples, depth<=%dx median" % (need, mult), 100 * keep.mean(),
            keep.sum() * scale / 1e6, (keep & seg).sum() * scale / 1e6))
