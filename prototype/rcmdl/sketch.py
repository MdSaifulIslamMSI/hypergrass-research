#!/usr/bin/env python3
"""Phase-insensitive k-mer sketches for panel similarity and ordering.

Hashes every canonical k-mer, samples by hash (not by position), and reports
Jaccard/containment. Sampling by hash makes the estimate insensitive to indels,
which broke the fixed-stride proxy.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

K = 16
BASE = np.uint64(0x9E3779B185EBCA87)
MASK = np.uint64(0xFF)  # ~1/256 of k-mers, ~20k per 5 Mb genome
LUT = np.zeros(256, dtype=np.uint8)
for code, base in enumerate(b"ACGT"):
    LUT[base] = code


def rolling_hashes(codes: np.ndarray, k: int) -> np.ndarray:
    n = len(codes) - k + 1
    h = np.zeros(n, dtype=np.uint64)
    for j in range(k):
        h = h * BASE + codes[j : j + n].astype(np.uint64)
    return h


def mix(h: np.ndarray) -> np.ndarray:
    h = h ^ (h >> np.uint64(30))
    h = h * np.uint64(0xBF58476D1CE4E5B9)
    h = h ^ (h >> np.uint64(27))
    h = h * np.uint64(0x94D049BB133111EB)
    return h ^ (h >> np.uint64(31))


def sketch(path: Path) -> np.ndarray:
    seq = np.frombuffer(path.read_bytes(), dtype=np.uint8)
    codes = LUT[seq]
    forward = rolling_hashes(codes, K)
    reverse = rolling_hashes((3 - codes)[::-1], K)[::-1]
    canonical = mix(np.minimum(forward, reverse))
    return np.unique(canonical[(canonical & MASK) == 0])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--genomes", type=Path, default=Path("prototype/rcmdl/work/genomes_masked"))
    parser.add_argument("--out", type=Path, default=Path("prototype/rcmdl/work/orders"))
    args = parser.parse_args()

    paths = sorted(args.genomes.glob("*.seq"))
    names = [p.name for p in paths]
    sketches = [sketch(p) for p in paths]
    print("sketch sizes:", [len(s) for s in sketches])

    n = len(paths)
    containment = [[0.0] * n for _ in range(n)]
    jaccard = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                containment[i][j] = 1.0
                jaccard[i][j] = 1.0
                continue
            inter = np.intersect1d(sketches[i], sketches[j], assume_unique=True).size
            containment[i][j] = inter / len(sketches[i])
            jaccard[i][j] = inter / (len(sketches[i]) + len(sketches[j]) - inter)

    # Prim order: start at the most central genome, then repeatedly add the genome
    # best covered by the already-placed set.
    centrality = [float(np.mean([containment[i][j] for j in range(n) if j != i])) for i in range(n)]
    start = int(np.argmax(centrality))
    order = [start]
    remaining = set(range(n)) - {start}
    while remaining:
        nxt = max(remaining, key=lambda j: max(containment[j][i] for i in order))
        order.append(nxt)
        remaining.remove(nxt)

    # Nearest-neighbour path: follow the strongest unvisited predecessor chain.
    chain = [start]
    remaining = set(range(n)) - {start}
    while remaining:
        nxt = max(remaining, key=lambda j: containment[j][chain[-1]])
        chain.append(nxt)
        remaining.remove(nxt)

    args.out.mkdir(parents=True, exist_ok=True)
    payload = {
        "k": K,
        "mask": int(MASK),
        "names": names,
        "sketch_sizes": [len(s) for s in sketches],
        "containment": containment,
        "jaccard": jaccard,
        "centrality": centrality,
        "prim_order": [names[i] for i in order],
        "chain_order": [names[i] for i in chain],
    }
    (args.out / "sketch.json").write_text(json.dumps(payload, indent=2) + "\n")
    (args.out / "prim.txt").write_text("\n".join(payload["prim_order"]) + "\n")
    (args.out / "chain.txt").write_text("\n".join(payload["chain_order"]) + "\n")

    print("\ncentrality (mean containment of this genome by others):")
    for i in sorted(range(n), key=lambda i: -centrality[i]):
        print("%2d %-8s %.4f" % (i, names[i], centrality[i]))
    print("\nprim:", payload["prim_order"])
    print("chain:", payload["chain_order"])


if __name__ == "__main__":
    main()
