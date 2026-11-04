# stage 9 summary: most related pairs from one ngsrelate run
# rab ~0.5 parent-offspring / full sibs, 0.25 half sibs, ~0 unrelated
# king ~0.25 first degree, 0.125 second degree, 0.5 duplicate. ~0.5 on both = duplicate
import sys

import pandas as pd

pd.set_option("display.width", 250)
r = pd.read_csv(sys.argv[1], sep="\t")
cols = [c for c in ["ida", "idb", "nSites", "rab", "KING", "theta", "Fa", "Fb"] if c in r.columns]
print("%s: %d pairs; rab median %.3f" % (sys.argv[1], len(r), r["rab"].median()))
print(r.sort_values("rab", ascending=False)[cols].head(15).round(3).to_string(index=False))
