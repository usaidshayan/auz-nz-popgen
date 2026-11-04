# usage: python scripts/to_zarr.py <input.vcf.gz> <output.zarr>
import os
import sys

import allel
import numcodecs

vcf, out = sys.argv[1], sys.argv[2]
if os.path.exists(out):
    sys.exit("output exists, remove it first: " + out)
fields = ["samples", "variants/CHROM", "variants/POS", "variants/REF", "variants/ALT", "variants/QUAL",
          "calldata/GT", "calldata/GQ", "calldata/DP", "calldata/AD"]
allel.vcf_to_zarr(vcf, out, fields=fields, alt_number=3, chunk_length=2**18,
                  compressor=numcodecs.Blosc(cname="zstd", clevel=1, shuffle=-1), log=sys.stdout)
print("wrote", out)
