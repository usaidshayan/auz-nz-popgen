# newhybrids input, following croft et al. 2024 (200 random loci, 10 replicates)
# genotypes called from the thinned beagle: most likely genotype if its normalised gl >= 0.9, else missing
# at low depth a het with one read of each allele passes easily while a hom needs ~4 reads,
# so calls lean towards hets, which if anything pushes samples towards hybrid classes
# loci: >= 50% of samples called overall and >= 9 of 18 called in each species, maf >= 0.05 among calls
# (asking for near-complete calls left ~800 loci, biased to deep sites and fixed paralog hets)
# usage: python s14_nh_input.py <thin beagle.gz> <sample list> <out dir> [loci per replicate, default 200]
import os
import sys

import numpy as np
import pandas as pd

beagle, lst, out = sys.argv[1:4]
NLOC = int(sys.argv[4]) if len(sys.argv) > 4 else 200
NREP, PMIN = 10, 0.9
names = [l.strip() for l in open(lst) if l.strip()]
sp = pd.read_csv("samples_meta.tsv", sep="\t").set_index("sample").loc[names].species.values
n = len(names)

keep_m, keep_g = [], []
for ch in pd.read_csv(beagle, sep="\t", chunksize=200000):
    g = ch.iloc[:, 3:].to_numpy(np.float32).reshape(len(ch), n, 3)
    call = g.argmax(2)
    call[g.max(2) < PMIN] = -1
    ok = call >= 0
    ac = np.where(ok, call, 0).sum(1)
    af = ac / np.maximum(2 * ok.sum(1), 1)
    sel = (ok.sum(1) >= 0.5 * n) & (ok[:, sp == "hilli"].sum(1) >= 9) & (ok[:, sp == "stygia"].sum(1) >= 9) \
        & (np.minimum(af, 1 - af) >= 0.05)
    keep_m.append(ch.marker.to_numpy()[sel])
    keep_g.append(call[sel])
mk = np.concatenate(keep_m)
gt = np.concatenate(keep_g)
print("%d loci pass the call filters" % len(mk))

code = {-1: "0", 0: "11", 1: "12", 2: "22"}
os.makedirs(out, exist_ok=True)
for r in range(1, NREP + 1):
    idx = np.sort(np.random.default_rng(r).choice(len(mk), NLOC, replace=False))
    with open("%s/r%02d.txt" % (out, r), "w") as f:
        f.write("NumIndivs %d\nNumLoci %d\nDigits 1\nFormat Lumped\n\n" % (n, NLOC))
        f.write("LocusNames %s\n\n" % " ".join("L%d" % (k + 1) for k in range(NLOC)))
        for i, s in enumerate(names):
            f.write("%d n %s %s\n" % (i + 1, s, " ".join(code[c] for c in gt[idx, i])))
    pd.Series(mk[idx]).to_csv("%s/r%02d.loci.txt" % (out, r), index=False, header=False)
    miss = (gt[idx] < 0).mean(0)
    print("replicate %d: missing per sample %.2f-%.2f" % (r, miss.min(), miss.max()))
