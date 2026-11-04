# stage 2 summary: where each sample's reads land when all four assemblies compete
import glob
import pandas as pd

pd.set_option("display.width", 250)
meta = pd.read_csv("samples_meta.tsv", sep="\t").set_index("sample")
rows = []
for f in sorted(glob.glob("compmap/qc/*.assign.tsv")):
    d = pd.read_csv(f, sep="\t", header=None, names=["sample", "kind", "asm", "n"])
    s = d["sample"].iloc[0]
    q = d[d.kind == "mapq20"].set_index("asm")["n"]
    low = int(d.loc[d.kind == "mapq_lt20", "n"].sum())
    un = int(d.loc[d.kind == "unmapped", "n"].sum())
    tot = q.sum() + low + un
    r = {"sample": s, "species": meta.loc[s, "species"], "origin": meta.loc[s, "country"]}
    for a, n in q.items():
        r["placed_" + a] = 100 * n / q.sum()          # % of placed reads
    r["pct_placed"] = 100 * q.sum() / tot
    r["pct_ambiguous"] = 100 * low / tot
    r["pct_unmapped"] = 100 * un / tot
    cov = pd.read_csv("compmap/qc/%s.coverage.tsv" % s, sep="\t")
    cov["asm"] = cov["#rname"].str.split("__").str[0]
    cov["len"] = cov["endpos"] - cov["startpos"] + 1
    for a, g in cov.groupby("asm"):
        r["depth_" + a] = (g["meandepth"] * g["len"]).sum() / g["len"].sum()
    rows.append(r)
tab = pd.DataFrame(rows).fillna(0)
place = sorted(c for c in tab.columns if c.startswith("placed_"))
tab["best_asm"] = tab[place].idxmax(axis=1).str.replace("placed_", "")
tab["own_share"] = [row["placed_" + row["species"]] for _, row in tab.iterrows()]
tab["flag"] = ""
tab.loc[tab.best_asm != tab.species, "flag"] += "MOST READS ON ANOTHER SPECIES; "
tab.loc[(tab.best_asm == tab.species) & (tab.own_share < 90), "flag"] += "own share <90% (hybrid or mix?); "
tab.loc[tab.pct_unmapped > 10, "flag"] += "unmapped >10%; "
tab = tab.sort_values(["species", "own_share"]).round(2)
tab.to_csv("tables/S2_competitive_mapping.tsv", sep="\t", index=False)
cols = ["sample", "species", "origin"] + place + ["pct_placed", "pct_ambiguous", "pct_unmapped", "flag"]
print("Share of confidently placed reads (MAPQ>=20) per assembly, then % of all reads placed / ambiguous / unmapped")
print(tab[cols].to_string(index=False))
print("\nMean depth on each assembly (all MAPQ):")
print(tab[["sample", "species"] + sorted(c for c in tab.columns if c.startswith("depth_"))].to_string(index=False))
print("written to tables/S2_competitive_mapping.tsv")
