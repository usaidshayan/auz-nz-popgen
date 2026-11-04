# task 1 (ang, 23 oct): mapping to own vs other assembly, unmapped % for both, by own-assembly coverage
# alignment depth = samtools coverage mean (all reads), effective depth = what bcftools used
# (mapq/bq >= 20, no dups), from the main sample qc tables (../qc)
import pandas as pd

pd.set_option("display.width", 250)
P = "../qc/"
t = pd.read_csv(P + "ref_compare.tsv", sep="\t")
other = {"hilli": "stygia", "stygia": "hilli"}
eff = {}
for sp in ("hilli", "stygia"):
    eff[("own", sp)] = pd.read_csv(P + "sample_qc_%s_on_%s.tsv" % (sp, sp), sep="\t").set_index("sample")["mean_DP"]
    eff[("other", sp)] = pd.read_csv(P + "sample_qc_all_on_%s.tsv" % other[sp], sep="\t").set_index("sample")["mean_DP"]
rows = []
for _, r in t.iterrows():
    own, oth = r["species"], other[r["species"]]
    rows.append({
        "sample": r["sample"], "species": own, "origin": r["country"],
        "mapped_own": r["pct_mapped_" + own], "mapped_other": r["pct_mapped_" + oth],
        "unmapped_own": 100 - r["pct_mapped_" + own], "unmapped_other": 100 - r["pct_mapped_" + oth],
        "depth_own": r["mean_depth_" + own], "depth_other": r["mean_depth_" + oth],
        "depth_kept_pct": 100 * r["mean_depth_" + oth] / r["mean_depth_" + own],
        "effDP_own": eff[("own", own)].get(r["sample"]), "effDP_other": eff[("other", own)].get(r["sample"]),
    })
tab = pd.DataFrame(rows).sort_values("depth_own").round(2)
tab.to_csv("tables/T1_mapping_depth.tsv", sep="\t", index=False)
print(tab.to_string(index=False))
print("\nMeans by species:")
print(tab.groupby("species")[["mapped_own", "mapped_other", "unmapped_own", "unmapped_other",
                              "depth_own", "depth_other", "effDP_own", "effDP_other"]].mean().round(2).to_string())
print("\nown mapping higher than other in %d of %d samples; own depth higher in %d of %d"
      % ((tab.mapped_own > tab.mapped_other).sum(), len(tab), (tab.depth_own > tab.depth_other).sum(), len(tab)))
top = tab.sort_values("depth_own").iloc[-1]
print("highest coverage: %s, %.2fx alignment depth, %.2fx effective depth" % (top["sample"], top["depth_own"], top["effDP_own"]))
print("written to tables/T1_mapping_depth.tsv")
