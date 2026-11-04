#!/usr/bin/env python3
"""Choose collection orders from sampled k-mer Jaccard similarity."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

K = 16
STRIDE = 251
COMP = bytes.maketrans(b"ACGT", b"TGCA")


def kmers(path: Path):
    seq = path.read_bytes()
    seen = set()
    for i in range(0, len(seq) - K + 1, STRIDE):
        kmer = seq[i : i + K]
        if any(c not in b"ACGT" for c in kmer):
            continue
        rc = kmer.translate(COMP)[::-1]
        seen.add(min(kmer, rc))
    return seen


def jaccard(a: set[bytes], b: set[bytes]) -> float:
    if not a or not b:
        return 0.0
    inter = len(a & b)
    return inter / (len(a) + len(b) - inter)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--genomes", type=Path, default=Path("prototype/rcmdl/work/genomes_masked"))
    parser.add_argument("--out", type=Path, default=Path("prototype/rcmdl/work/orders"))
    args = parser.parse_args()
    paths = sorted(args.genomes.glob("*.seq"))
    names = [p.name for p in paths]
    sets = [kmers(p) for p in paths]
    n = len(paths)
    sim = [[1.0 if i == j else jaccard(sets[i], sets[j]) for j in range(n)] for i in range(n)]
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "jaccard.json").write_text(json.dumps({"names": names, "sim": sim}, indent=2) + "\n")

    print("mean similarity to others (centrality):")
    centrality = [(sum(sim[i]) - 1.0) / (n - 1) for i in range(n)]
    for i in sorted(range(n), key=lambda i: -centrality[i]):
        print("%2d %-8s %.4f" % (i, names[i], centrality[i]))

    orders = {}
    orders["name"] = names
    orders["size-asc"] = [p.name for p in sorted(paths, key=lambda p: p.stat().st_size)]
    orders["size-desc"] = [p.name for p in sorted(paths, key=lambda p: -p.stat().st_size)]

    start = max(range(n), key=lambda i: centrality[i])
    greedy = [start]
    remaining = set(range(n)) - {start}
    while remaining:
        nxt = max(remaining, key=lambda j: max(sim[i][j] for i in greedy))
        greedy.append(nxt)
        remaining.remove(nxt)
    orders["greedy-similarity"] = [names[i] for i in greedy]

    # Prefer placing a near-duplicate immediately after its best match: pure single-linkage
    # chain can bury a cluster late, so start from the most central genome and always attach
    # the best remaining neighbour of the most recently placed genome when it is strong.
    chain = [start]
    remaining = set(range(n)) - {start}
    while remaining:
        last = chain[-1]
        best = max(remaining, key=lambda j: sim[last][j])
        if sim[last][best] >= 0.5:
            nxt = best
        else:
            nxt = max(remaining, key=lambda j: max(sim[i][j] for i in chain))
        chain.append(nxt)
        remaining.remove(nxt)
    orders["greedy-chain"] = [names[i] for i in chain]

    for name, order in orders.items():
        (args.out / (name + ".txt")).write_text("\n".join(order) + "\n")
        print(name, order)


if __name__ == "__main__":
    main()
