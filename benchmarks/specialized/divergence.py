"""Estimates pairwise 21-mer Jaccard similarity between genomes via MinHash sketches."""
import json
import random
import sys
from pathlib import Path

MASK = (1 << 64) - 1
K = 21
SKETCH = 4000


def kmers(seq: bytes, limit: int):
    h = 0
    out = []
    mask = (1 << (2 * K)) - 1
    code = {65: 0, 67: 1, 71: 2, 84: 3}
    run = 0
    for b in seq[:limit]:
        v = code.get(b)
        if v is None:
            run = 0
            h = 0
            continue
        h = ((h << 2) | v) & mask
        run += 1
        if run >= K:
            out.append(h)
    return out


def sketch(seq: bytes, limit: int, seed: int = 1):
    hs = kmers(seq, limit)
    hs = {((x * 0x9E3779B97F4A7C15) & MASK) for x in hs}
    return sorted(hs)[:SKETCH]


def read_records(path: Path):
    recs, cur = [], []
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if cur:
                    recs.append(b"".join(cur))
                    cur = []
            else:
                cur.append(line.strip())
    if cur:
        recs.append(b"".join(cur))
    return recs


def main():
    path = Path(sys.argv[1])
    limit = int(sys.argv[2]) if len(sys.argv) > 2 else 1_000_000
    recs = read_records(path)
    print(f"{len(recs)} records, {limit} bp each")
    sks = [sketch(r, limit) for r in recs]
    sets = [set(s) for s in sks]
    for i in range(len(recs)):
        row = []
        for j in range(len(recs)):
            inter = len(sets[i] & sets[j])
            union = len(sets[i] | sets[j])
            row.append(inter / union if union else 0.0)
        if i in (0, 1, 2, 12, 23):
            print(f"g{i:02d} nearest: " + " ".join(f"{j}:{v:.3f}" for j, v in enumerate(row) if j != i and v > 0.3))
    # mean pairwise jaccard and implied divergence (1 - J)^(1/K) approximation
    vals = []
    for i in range(len(recs)):
        for j in range(i + 1, len(recs)):
            inter = len(sets[i] & sets[j])
            union = len(sets[i] | sets[j])
            vals.append(inter / union if union else 0.0)
    mean = sum(vals) / len(vals)
    d = 1 - mean ** (1.0 / K)
    print(f"mean pairwise Jaccard = {mean:.4f}; implied per-base divergence ~ {d*100:.2f}%")


if __name__ == "__main__":
    main()
