# stage 3 summary: kraken2 groups for the unmapped pairs, as % of unmapped and % of all pairs
# insecta minus calliphoridae / diptera kept separate so fly reads missing from the assemblies show up
import glob
import os
import re

import pandas as pd

GROUPS = [  # (column, taxid)
    ("bacteria", 2), ("archaea", 2157), ("viruses", 10239), ("eukaryota", 2759),
    ("fungi", 4751), ("plants", 33090), ("metazoa", 33208), ("insecta", 50557),
    ("diptera", 7147), ("calliphoridae", 7371), ("mammalia", 40674), ("human", 9606),
]
s2 = pd.read_csv("tables/S2_competitive_mapping.tsv", sep="\t")
rows = []
for f in sorted(glob.glob("compmap/kraken/*.report")):
    s = os.path.basename(f).replace(".report", "")
    r = pd.read_csv(f, sep="\t", header=None, names=["pct", "clade", "direct", "rank", "taxid", "name"])
    r["name"] = r["name"].str.strip()
    clade = dict(zip(r.taxid, r.clade))
    unc, root = clade.get(0, 0), clade.get(1, 0)
    n = unc + root
    d = {"sample": s, "unmapped_pairs": n, "classified": 100 * root / n}
    for col, t in GROUPS:
        d[col] = 100 * clade.get(t, 0) / n
    d["other_euk"] = d["eukaryota"] - d["fungi"] - d["plants"] - d["metazoa"]
    d["insecta_not_diptera"] = d["insecta"] - d["diptera"]
    d["diptera_not_calliphoridae"] = d["diptera"] - d["calliphoridae"]
    g = r[r["rank"] == "G"].sort_values("clade", ascending=False).head(3)
    d["top_genera"] = "; ".join("%s %.1f%%" % (a, 100 * b / n) for a, b in zip(g.name, g.clade))
    prim = re.search(r"(\d+) \+ \d+ primary$", open("compmap/qc/%s.flagstat.txt" % s).read(), re.M)
    d["pct_of_all_pairs"] = 100 * n / (int(prim.group(1)) / 2) if prim else float("nan")
    rows.append(d)
t = pd.DataFrame(rows)
t = s2[["sample", "species", "origin", "pct_unmapped", "own_share"]].merge(t, on="sample")
t = t.sort_values("pct_unmapped", ascending=False)
t.round(2).to_csv("tables/S3_unmapped_kraken.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
cols = ["sample", "species", "pct_unmapped", "classified", "bacteria", "fungi", "plants", "other_euk",
        "calliphoridae", "diptera_not_calliphoridae", "insecta_not_diptera", "human", "viruses", "top_genera"]
print(t[cols].round(1).to_string(index=False))
