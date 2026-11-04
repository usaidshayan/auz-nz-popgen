import os
import pandas as pd

raw = pd.read_excel("N=40SampleList.xlsx", usecols=[0, 1, 2]).dropna(subset=["Sample code"])
raw.columns = ["sample_code", "location", "species"]
meta = pd.DataFrame({
    "sample": raw["sample_code"].str.replace("-", "_"),
    "species": raw["species"].str.lower(),
    "location": raw["location"],
})
meta["country"] = meta["location"].apply(lambda x: "AUS" if x.endswith("_AUS") else "NZ")

files = [l.strip() for l in open("samples.txt") if l.strip()]
print("samples in sheet:", len(meta), "| sample pairs in samples.txt:", len(files))
print("in sheet but no reads:", sorted(set(meta["sample"]) - set(files)))
print("reads but not in sheet:", sorted(set(files) - set(meta["sample"])))
bad = meta[(meta["sample"].str.contains("_CH_") & (meta["species"] != "hilli")) |
           (meta["sample"].str.contains("_CS_") & (meta["species"] != "stygia"))]
print("code/species mismatches:", len(bad))

meta.to_csv("samples_meta.tsv", sep="\t", index=False)
os.makedirs("lists", exist_ok=True)
meta["sample"].to_csv("lists/all.txt", index=False, header=False)
for sp in ("hilli", "stygia"):
    meta.loc[meta["species"] == sp, "sample"].to_csv("lists/%s.txt" % sp, index=False, header=False)
print(meta.groupby(["country", "species"]).size().to_string())
