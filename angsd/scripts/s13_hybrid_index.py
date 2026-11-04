# hybrid index and interspecific heterozygosity from gls, for the triangle plot
# species allele freqs by em on the beagle gls. a site is diagnostic for a sample if the two species
# differ by >= 0.9 in freq, with that sample left out of its own species' freq (no circularity)
# and >= 10 of 18 per species with data. the sample's genotype class freqs at its diagnostic sites
# (hom own, het, hom other) are fitted by em on its gls, so low depth is handled by the likelihoods
# hybrid index = share of other-species alleles, q12 = interspecific het
# pure ~ (0, 0), f1 (0.5, 1), f2 (0.5, 0.5), backcross (0.25, 0.5). se by delete-block jackknife
# usage: python s13_hybrid_index.py <beagle.gz> <sample list> <out prefix>
import sys

import numpy as np
import pandas as pd

beagle, lst, out = sys.argv[1:4]
DIFF, MIN_N, BLOCK = 0.9, 10, 200000
names = [l.strip() for l in open(lst) if l.strip()]
meta = pd.read_csv("samples_meta.tsv", sep="\t").set_index("sample").loc[names]
sp = meta.species.values
ih, ist = np.where(sp == "hilli")[0], np.where(sp == "stygia")[0]
n = len(names)


def em_freq(g, iters=40):
    # g: sites x inds x 3 gls. returns minor allele freq and each ind's posterior dosage
    p = np.clip((g[:, :, 1] + 2 * g[:, :, 2]).mean(1) / 2, 0.01, 0.99)
    for _ in range(iters):
        pr = np.stack([(1 - p) ** 2, 2 * p * (1 - p), p ** 2], 1)[:, None, :]
        w = g * pr
        w /= w.sum(2, keepdims=True)
        dos = w[:, :, 1] + 2 * w[:, :, 2]
        p = np.clip(dos.mean(1) / 2, 1e-6, 1 - 1e-6)
    return p, dos


keep = [[] for _ in range(n)]
blocks = [[] for _ in range(n)]
start = 0
for ch in pd.read_csv(beagle, sep="\t", chunksize=250000):
    g = ch.iloc[:, 3:].to_numpy(np.float32).reshape(len(ch), n, 3)
    has = (g.max(2) - g.min(2)) > 1e-3
    ok = (has[:, ih].sum(1) >= MIN_N) & (has[:, ist].sum(1) >= MIN_N)
    g, has = g[ok], has[ok]
    blk = (start + np.where(ok)[0]) // BLOCK
    start += len(ch)
    p_h, d_h = em_freq(g[:, ih])
    p_s, d_s = em_freq(g[:, ist])
    for idx, other_p, d in [(ih, p_s, d_h), (ist, p_h, d_s)]:
        m = len(idx)
        tot = d.sum(1)
        for j, i in enumerate(idx):
            p_own = (tot - d[:, j]) / (2 * (m - 1))
            delta = other_p - p_own
            sel = (np.abs(delta) >= DIFF) & has[:, i]
            gi = g[sel, i, :]
            flip = delta[sel] < 0  # other species carries the major allele, so swap hom classes
            gi[flip] = gi[flip][:, ::-1]
            keep[i].append(gi)
            blocks[i].append(blk[sel])
    print("%d sites read" % start, flush=True)


def em_pi(L, pi=None, tol=1e-8, iters=2000):
    pi = np.full(3, 1 / 3) if pi is None else pi.copy()
    for _ in range(iters):
        w = L * pi
        w /= w.sum(1, keepdims=True)
        new = w.mean(0)
        if np.abs(new - pi).max() < tol:
            return new
        pi = new
    return pi


def stats(pi):
    return pi[1] / 2 + pi[2], pi[1]


rows = []
for i, s in enumerate(names):
    L = np.concatenate(keep[i]).astype(np.float64)
    b = np.concatenate(blocks[i])
    pi = em_pi(L)
    h, q = stats(pi)
    ub = np.unique(b)
    jk = np.array([stats(em_pi(L[b != k], pi)) for k in ub])
    B = len(ub)
    se = np.sqrt((B - 1) / B * ((jk - jk.mean(0)) ** 2).sum(0))
    h_sty = h if sp[i] == "hilli" else 1 - h
    rows.append([s, sp[i], meta.country.iloc[i], meta.location.iloc[i], len(L), h, se[0], q, se[1], h_sty,
                 pi[0], pi[1], pi[2]])
    print("%s %s sites %d other-species ancestry %.4f (se %.4f) q12 %.4f (se %.4f)"
          % (s, sp[i], len(L), h, se[0], q, se[1]), flush=True)
t = pd.DataFrame(rows, columns=["sample", "species", "country", "location", "n_diag_sites", "other_ancestry",
                                "other_ancestry_se", "q12", "q12_se", "stygia_fraction", "p_hom_own", "p_het",
                                "p_hom_other"])
t.to_csv(out + ".hybrid_index.tsv", sep="\t", index=False, float_format="%.5f")
print("wrote", out + ".hybrid_index.tsv")
