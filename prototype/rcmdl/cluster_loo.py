#!/usr/bin/env python3
"""Leave-cluster-out benchmark for the RC-MDL prototype.

Standard leave-one-out leaves 24 near-identical genomes in the panel, so a target
with a twin is nearly free. This harness removes the target *and every genome in its
near-identical cluster* from the panel, which is the stricter "unseen lineage" test.

Clusters come from the hash-sampled 16-mer Jaccard matrix written by sketch.py.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import time
from pathlib import Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def invoke(cmd: list[str]) -> tuple[dict, float]:
    start = time.perf_counter()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.perf_counter() - start
    if proc.returncode != 0:
        raise RuntimeError(
            "command failed (%d): %s\n%s\n%s"
            % (proc.returncode, " ".join(cmd), proc.stdout, proc.stderr)
        )
    lines = [line for line in proc.stdout.splitlines() if line.startswith("{")]
    if not lines:
        raise RuntimeError("no stats JSON from: %s" % " ".join(cmd))
    return json.loads(lines[-1]), elapsed


def clusters_at(jaccard: list[list[float]], threshold: float) -> list[int]:
    n = len(jaccard)
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i in range(n):
        for j in range(i + 1, n):
            if jaccard[i][j] >= threshold:
                a, b = find(i), find(j)
                if a != b:
                    parent[b] = a
    return [find(i) for i in range(n)]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--genomes", type=Path, default=Path("prototype/rcmdl/work/saureus/genomes"))
    parser.add_argument("--sketch", type=Path, default=Path("prototype/rcmdl/work/orders-saureus/sketch.json"))
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--threshold", type=float, default=0.98)
    parser.add_argument("--java", default="java")
    parser.add_argument("--classpath", default="prototype/rcmdl/bin")
    parser.add_argument("--heap", default="1400m")
    parser.add_argument("--k", type=int, default=16)
    parser.add_argument("--step", type=int, default=8)
    args = parser.parse_args()

    genomes = sorted(args.genomes.glob("*.seq"))
    sketch = json.loads(args.sketch.read_text())
    names = sketch["names"]
    if names != [p.name for p in genomes]:
        raise SystemExit("sketch names do not match genome directory order")
    roots = clusters_at(sketch["jaccard"], args.threshold)
    members: dict[int, list[int]] = {}
    for i, root in enumerate(roots):
        members.setdefault(root, []).append(i)

    args.out.mkdir(parents=True, exist_ok=True)
    base = [args.java, "-Xmx" + args.heap, "-cp", args.classpath, "Rcmdl"]
    list_path = args.out / "panel.txt"
    records = []

    for index, target in enumerate(genomes):
        excluded = members[roots[index]]
        panel = [g for i, g in enumerate(genomes) if i not in excluded]
        list_path.write_text(
            "\n".join(p.name for p in panel) + "\n", encoding="utf-8"
        )
        archive = args.out / ("%02d-%s.rcmd" % (index, target.stem))
        restored = args.out / ("%02d-%s.out" % (index, target.stem))
        compress = base + [
            "compress",
            "--target",
            str(target),
            "--out",
            str(archive),
            "--panel",
            str(args.genomes),
            "--panel-list",
            str(list_path),
            "--k",
            str(args.k),
            "--step",
            str(args.step),
        ]
        stats, wall = invoke(compress)
        dstats, dwall = invoke(
            base
            + [
                "decompress",
                "--archive",
                str(archive),
                "--panel",
                str(args.genomes),
                "--out",
                str(restored),
            ]
        )
        actual = sha256(restored)
        expected = sha256(target)
        if actual != expected:
            raise RuntimeError("round-trip mismatch for %s" % target)
        size = archive.stat().st_size
        record = {
            "index": index,
            "target": target.name,
            "target_bytes": target.stat().st_size,
            "target_sha256": expected,
            "excluded": [genomes[i].name for i in excluded],
            "panel_genomes": len(panel),
            "archive_bytes": size,
            "bpb": size * 8.0 / target.stat().st_size,
            "panel_match_fraction": stats["panel_match_bases"] / max(1, target.stat().st_size),
            "literal_fraction": stats["literal_bases"] / max(1, target.stat().st_size),
            "mismatches": stats["mismatches"],
            "compress_seconds": stats["compress_seconds"],
            "compress_wall_seconds": wall,
            "decompress_seconds": dstats["decompress_seconds"],
            "decompress_wall_seconds": dwall,
            "verified": True,
        }
        records.append(record)
        restored.unlink()
        (args.out / "results.json").write_text(
            json.dumps(records, indent=2) + "\n", encoding="utf-8"
        )
        print(
            "%2d %-6s exclude=%d panel=%2d %9d  %8.4f bpb  %5.1fs"
            % (
                index,
                target.stem,
                len(excluded) - 1,
                len(panel),
                size,
                record["bpb"],
                wall,
            ),
            flush=True,
        )

    total = sum(r["archive_bytes"] for r in records)
    target_bases = sum(r["target_bytes"] for r in records)
    summary = {
        "genomes": len(records),
        "threshold": args.threshold,
        "target_bases": target_bases,
        "archive_bytes": total,
        "bpb": total * 8.0 / target_bases,
        "verified": all(r["verified"] for r in records),
        "clusters": [
            [genomes[i].name for i in idx] for idx in members.values()
        ],
    }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
