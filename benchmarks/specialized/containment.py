"""Estimates local per-base divergence between genomes with a 21-mer Bloom containment test.

Rearranged genomes break whole-window Jaccard, so instead every 20 kb window of one genome is
tested against ALL k-mers of another genome. Containment of homologous windows is ~(1-d)^k.
"""
import sys
from pathlib import Path

K = 21
BITS = 1 << 28
MASK64 = (1 << 64) - 1
CODE = {65: 0, 67: 1, 71: 2, 84: 3}


def read_records(path: Path):
    recs, cur = [], []
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if cur:
                    recs.append(b"".join(cur))
                    cur = []
            else:
                cur.append(line.strip().upper())
    if cur:
        recs.append(b"".join(cur))
    return recs


def kmers(seq):
    h = 0
    run = 0
    mask = (1 << (2 * K)) - 1
    for b in seq:
        v = CODE.get(b)
        if v is None:
            run = 0
            h = 0
            continue
        h = ((h << 2) | v) & mask
        run += 1
        if run >= K:
            yield h


def mark(table, h):
    idx = h & (BITS - 1)
    table[idx >> 3] |= 1 << (idx & 7)


def probe(table, h):
    idx = h & (BITS - 1)
    return (table[idx >> 3] >> (idx & 7)) & 1


def add_all(table, seq):
    for x in kmers(seq):
        h1 = (x * 0x9E3779B97F4A7C15) & MASK64
        h2 = (x * 0xC2B2AE3D27D4EB4F) & MASK64
        mark(table, h1)
        mark(table, h2)


def contains(table, seq):
    total = 0
    hit = 0
    for x in kmers(seq):
        h1 = (x * 0x9E3779B97F4A7C15) & MASK64
        h2 = (x * 0xC2B2AE3D27D4EB4F) & MASK64
        total += 1
        if probe(table, h1) and probe(table, h2):
            hit += 1
    return hit / total if total else 0.0


def main():
    path = Path(sys.argv[1])
    recs = read_records(path)
    present = [i for i, r in enumerate(recs) if len(r) > 100_000]
    print(f"{len(recs)} records")
    table = bytearray(BITS // 8)
    target = present[int(sys.argv[2]) if len(sys.argv) > 2 else 1]
    print(f"indexing record {target} ({len(recs[target])} bp)")
    add_all(table, recs[target])
    print(f"self containment = {contains(table, recs[target]):.4f}")
    windows = 120
    for src in present:
        if src == target:
            continue
        step = max(1, (len(recs[src]) - 60_000) // windows)
        vals = []
        for w in range(windows):
            start = w * step
            vals.append(contains(table, recs[src][start:start + 20_000]))
        vals.sort()
        p90 = vals[int(0.9 * (len(vals) - 1))]
        if p90 > 0.05:
            d = 1 - p90 ** (1.0 / K)
            print(f"record {src:2d}: median containment {vals[len(vals)//2]:.3f} p90 {p90:.3f}"
                  f" -> divergence ~{d*100:.2f}%")


if __name__ == "__main__":
    main()
