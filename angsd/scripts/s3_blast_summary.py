# blast of 500 unmapped r1 reads per sample against nt (2026-08), best hit per read
# reads with a hit, split by kingdom, and the top genera
import pandas as pd

rows = []
for s in ["WE_CS_133", "WE_CS_139", "MKA_AUS_CH_211", "HO_AUS_CS_204", "SE_AUS_CH_209"]:
    b = pd.read_csv("compmap/blast/%s.blast.tsv" % s, sep="\t", header=None,
                    names=["q", "s", "pid", "len", "e", "bits", "taxid", "name", "kingdom"])
    b = b.sort_values("bits", ascending=False).drop_duplicates("q")
    b["genus"] = b.name.str.split().str[0]
    k = b.kingdom.value_counts()
    g = b.genus.value_counts().head(5)
    rows.append({"sample": s, "reads_blasted": 500, "reads_with_hit": len(b),
                 "bacteria": int(k.get("Bacteria", 0)), "eukaryota": int(k.get("Eukaryota", 0)),
                 "viruses": int(k.get("Viruses", 0)), "median_pident": round(b.pid.median(), 1),
                 "top_genera": "; ".join("%s %d" % (n, c) for n, c in g.items())})
t = pd.DataFrame(rows)
t.to_csv("tables/S3b_unmapped_blast.tsv", sep="\t", index=False)
print(t.to_string(index=False))
