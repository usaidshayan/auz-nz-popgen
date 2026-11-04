# usage: python scripts/zarr_check.py <store.zarr> <source.vcf.gz> <lists/SET.txt>
import subprocess
import sys

import allel
import zarr

store, vcf, lst = sys.argv[1:4]
z = zarr.open_group(store, mode="r")
samples = [str(s) for s in z["samples"][:]]
gt = z["calldata/GT"]
n_vcf = int(subprocess.check_output(["bcftools", "index", "-n", vcf]).decode().strip())
n_zarr = z["variants/POS"].shape[0]
expected = sorted(l.strip() for l in open(lst) if l.strip())
print("GT shape (sites, samples, ploidy):", gt.shape)
print("records in VCF: %d | sites in zarr: %d | equal: %s" % (n_vcf, n_zarr, n_vcf == n_zarr))
print("sample names match %s: %s" % (lst, sorted(samples) == expected))
n = min(1000000, gt.shape[0])
g = allel.GenotypeArray(gt[:n])
ac = g.count_alleles()
called = ac.sum(axis=1) > 0
seg = ac.is_segregating()
print("first %d sites: %d segregating, %d invariant, %d with no called genotype in any sample"
      % (n, seg.sum(), (called & ~seg).sum(), (~called).sum()))
print("fraction of missing genotypes in that slice: %.3f" % g.is_missing().mean())
dp = z["calldata/DP"][:n]
print("mean per-sample DP in that slice (all sites): %.2f" % dp[dp >= 0].mean())
