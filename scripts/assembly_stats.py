import sys

def summarise(path):
    lens, n_count, lower_count = [], 0, 0
    lines = []
    with open(path, "rb") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(b">"):
                if lines:
                    seq = b"".join(lines)
                    lens.append(len(seq))
                    n_count += seq.count(b"N")
                    lower_count += sum(seq.count(c) for c in (b"a", b"c", b"g", b"t", b"n"))
                lines = []
            else:
                lines.append(line)
    if lines:
        seq = b"".join(lines)
        lens.append(len(seq))
        n_count += seq.count(b"N")
        lower_count += sum(seq.count(c) for c in (b"a", b"c", b"g", b"t", b"n"))
    total = sum(lens)
    acc = 0
    n50 = 0
    for length in sorted(lens, reverse=True):
        acc += length
        if acc >= total / 2:
            n50 = length
            break
    return len(lens), total, n50, max(lens), 100.0 * n_count / total, 100.0 * lower_count / total

print("file\tn_seqs\ttotal_bp\tN50\tlongest\tpct_N\tpct_lowercase")
for path in sys.argv[1:]:
    n, total, n50, longest, pct_n, pct_low = summarise(path)
    print("%s\t%d\t%d\t%d\t%d\t%.2f\t%.2f" % (path.split("/")[-1], n, total, n50, longest, pct_n, pct_low))
