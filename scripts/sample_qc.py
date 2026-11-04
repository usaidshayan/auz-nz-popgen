# usage: python scripts/sample_qc.py <store.zarr>
# per-sample coverage, missingness and heterozygosity from every 20th block of sites along the genome
import os
import sys

import numpy as np
import pandas as pd
import zarr

store = sys.argv[1]
tag = os.path.basename(store.rstrip("/")).replace(".zarr", "")
reftag = tag.split("_on_")[1]
z = zarr.open_group(store, mode="r")
samples = [str(s) for s in z["samples"][:]]
gt_arr, dp_arr = z["calldata/GT"], z["calldata/DP"]
n_sites = gt_arr.shape[0]
step = gt_arr.chunks[0]
starts = list(range(0, n_sites, step))[::20]
k = len(samples)
tot = 0
dp_sum, miss, ge4, het, called4 = (np.zeros(k) for _ in range(5))
for s in starts:
    gt = gt_arr[s:s + step]
    dp = dp_arr[s:s + step].astype(np.int32)
    called = (gt[:, :, 0] >= 0) & (gt[:, :, 1] >= 0)
    ok = dp >= 4
    c4 = called & ok
    tot += gt.shape[0]
    dp_sum += np.clip(dp, 0, None).sum(axis=0)
    miss += (~called).sum(axis=0)
    ge4 += ok.sum(axis=0)
    called4 += c4.sum(axis=0)
    het += (c4 & (gt[:, :, 0] != gt[:, :, 1])).sum(axis=0)
res = pd.DataFrame({"sample": samples, "mean_DP": dp_sum / tot, "pct_missing": 100 * miss / tot,
                    "pct_DP_ge4": 100 * ge4 / tot, "het_per_site_DP_ge4": het / np.where(called4 == 0, np.nan, called4)})
if os.path.exists("samples_meta.tsv"):
    res = res.merge(pd.read_csv("samples_meta.tsv", sep="\t")[["sample", "species", "country"]], on="sample", how="left")
col = "pct_mapped_" + reftag
if os.path.exists("qc/ref_compare.tsv"):
    res = res.merge(pd.read_csv("qc/ref_compare.tsv", sep="\t")[["sample", col]], on="sample", how="left")
res = res.round(4).sort_values("mean_DP")
os.makedirs("qc", exist_ok=True)
res.to_csv("qc/sample_qc_%s.tsv" % tag, sep="\t", index=False)
pd.set_option("display.width", 250)
print("%s: %d sites sampled of %d" % (tag, tot, n_sites))
print(res.to_string(index=False))
print("written to qc/sample_qc_%s.tsv" % tag)
