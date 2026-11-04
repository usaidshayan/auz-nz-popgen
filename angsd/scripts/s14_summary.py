# newhybrids summary: posterior class probabilities per sample, mean (min-max) over the 10 replicates
# for both assemblies and both locus counts (200 as in croft et al. 2024, 500)
# runs are unsupervised so pure_0 / pure_1 can swap between replicates; each replicate is turned so
# that c. hilli is pure_0 (and 0_bx is the hilli backcross) before averaging
import glob

import numpy as np
import pandas as pd

CL = ["pure_hilli", "pure_stygia", "F1", "F2", "bx_hilli", "bx_stygia"]
rows, rep = [], []
meta = pd.read_csv("samples_meta.tsv", sep="\t").set_index("sample")
for tag in ["all_on_hilli_200", "all_on_stygia_200", "all_on_hilli_500", "all_on_stygia_500"]:
    reps = []
    for f in sorted(glob.glob("newhybrids/%s/r*_run/aa-PofZ.txt" % tag)):
        p = pd.read_csv(f, sep="\t", skiprows=1, header=None).iloc[:, 1:8]
        p.columns = ["sample"] + CL
        p = p.set_index("sample")
        assert p.notna().all().all(), "nan in " + f
        if p.loc[meta.loc[p.index, "species"] == "hilli", "pure_hilli"].mean() < 0.5:
            p = p[["pure_stygia", "pure_hilli", "F1", "F2", "bx_stygia", "bx_hilli"]]
            p.columns = CL
        reps.append(p)
        rep.append([tag, f.split("/")[2], (p.idxmax(axis=1) == "pure_" + meta.loc[p.index, "species"]).sum(),
                    p[["F1", "F2", "bx_hilli", "bx_stygia"]].sum(axis=1).max()])
    a = np.stack([r.values for r in reps])
    for i, s in enumerate(reps[0].index):
        own = CL.index("pure_" + meta.loc[s, "species"])
        hyb = a[:, i, 2:].sum(1)
        rows.append([tag, s, meta.loc[s, "species"], len(reps), a[:, i, own].mean(), a[:, i, own].min(),
                     hyb.mean(), hyb.max(), CL[int(a[:, i].mean(0).argmax())]])
t = pd.DataFrame(rows, columns=["set", "sample", "species", "n_reps", "p_pure_own_mean", "p_pure_own_min",
                                "p_any_hybrid_mean", "p_any_hybrid_max", "assigned"])
t.to_csv("tables/S14_newhybrids.tsv", sep="\t", index=False, float_format="%.4f")
r = pd.DataFrame(rep, columns=["set", "replicate", "n_assigned_pure_own", "max_p_any_hybrid"])
r.to_csv("tables/S14_newhybrids_replicates.tsv", sep="\t", index=False, float_format="%.4f")
print(r.to_string(index=False))
print(t.groupby("set").agg(n=("sample", "size"), min_pure_own=("p_pure_own_min", "min"),
                           max_hybrid=("p_any_hybrid_max", "max")).to_string())
print(t[t.assigned != "pure_" + t.species].to_string(index=False) or "every sample assigned pure to its own species")
