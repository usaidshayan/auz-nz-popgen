import glob, json, os
header = ["sample", "reads_before", "reads_after", "pct_kept", "Gb_after", "Q30_pct", "GC_pct", "dup_pct"]
rows = []
for path in sorted(glob.glob("qc/fastp/*.fastp.json")):
    d = json.load(open(path))
    sample = os.path.basename(path).replace(".fastp.json", "")
    b = d["summary"]["before_filtering"]
    a = d["summary"]["after_filtering"]
    rows.append([sample, b["total_reads"], a["total_reads"], round(100.0 * a["total_reads"] / b["total_reads"], 1),
                 round(a["total_bases"] / 1e9, 3), round(100.0 * a["q30_rate"], 1),
                 round(100.0 * a["gc_content"], 1), round(100.0 * d["duplication"]["rate"], 1)])
with open("qc/fastp_summary.tsv", "w") as out:
    out.write("\t".join(header) + "\n")
    for r in rows:
        out.write("\t".join(str(x) for x in r) + "\n")
print("".join("%-16s" % h for h in header))
for r in rows:
    print("".join("%-16s" % x for x in r))
print("samples summarised:", len(rows))
