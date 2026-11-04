#!/usr/bin/env python3
"""Greedy marginal-order search for the progressive collection.

Prim order uses k-mer containment as a proxy. This harness measures the real
compressed marginal: at each step it codes every remaining candidate against the
genomes already placed and appends the cheapest one.

The resulting order can be fed to collection.py with --order-file. The sum of the
recorded marginals is the total the collection run will produce.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from pathlib import Path


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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--genomes", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--sketch", type=Path, default=None, help="optional sketch.json for shortlists")
    parser.add_argument("--candidates", type=int, default=0, help="evaluate only the N most contained candidates (0 = all)")
    parser.add_argument("--java", default="java")
    parser.add_argument("--classpath", default="prototype/rcmdl/bin")
    parser.add_argument("--heap", default="1400m")
    parser.add_argument("--k", type=int, default=16)
    parser.add_argument("--step", type=int, default=8)
    args = parser.parse_args()

    genomes = sorted(args.genomes.glob("*.seq"))
    if not genomes:
        raise SystemExit("no genomes")
    by_name = {p.name: p for p in genomes}
    containment = None
    if args.sketch:
        sketch = json.loads(args.sketch.read_text())
        names = sketch["names"]
        if names != [p.name for p in genomes]:
            raise SystemExit("sketch names do not match genome directory order")
        containment = sketch["containment"]

    args.out.mkdir(parents=True, exist_ok=True)
    base = [args.java, "-Xmx" + args.heap, "-cp", args.classpath, "Rcmdl"]
    list_path = args.out / "panel.txt"
    scratch = args.out / "candidate.rcmd"

    order: list[str] = []
    records: list[dict] = []
    total = 0

    while len(order) < len(genomes):
        remaining = [p for p in genomes if p.name not in order]
        if containment is not None and args.candidates > 0 and len(order) > 0:
            placed = [i for i, p in enumerate(genomes) if p.name in order]
            remaining.sort(
                key=lambda p: max(
                    containment[genomes.index(p)][i] for i in placed
                ),
                reverse=True,
            )
            remaining = remaining[: args.candidates]

        best = None
        for candidate in remaining:
            cmd = base + [
                "compress",
                "--target",
                str(candidate),
                "--out",
                str(scratch),
                "--k",
                str(args.k),
                "--step",
                str(args.step),
            ]
            if order:
                list_path.write_text("\n".join(order) + "\n", encoding="utf-8")
                cmd += ["--panel", str(args.genomes), "--panel-list", str(list_path)]
            stats, wall = invoke(cmd)
            size = scratch.stat().st_size
            row = {
                "candidate": candidate.name,
                "archive_bytes": size,
                "wall_seconds": wall,
            }
            if best is None or (size, candidate.name) < (best["archive_bytes"], best["candidate"]):
                best = row
        assert best is not None
        order.append(best["candidate"])
        total += best["archive_bytes"]
        records.append(
            {
                "step": len(order) - 1,
                "target": best["candidate"],
                "archive_bytes": best["archive_bytes"],
                "cumulative_bytes": total,
                "panel_genomes": len(order) - 1,
                "evaluated": len(remaining),
            }
        )
        (args.out / "order.txt").write_text("\n".join(order) + "\n", encoding="utf-8")
        (args.out / "results.json").write_text(
            json.dumps(records, indent=2) + "\n", encoding="utf-8"
        )
        print(
            "%2d %-6s %9d  cum %10d  (evaluated %d)"
            % (len(order) - 1, best["candidate"], best["archive_bytes"], total, len(remaining)),
            flush=True,
        )
        scratch.unlink(missing_ok=True)

    if scratch.exists():
        scratch.unlink()
    summary = {
        "genomes": len(order),
        "target_bases": sum(by_name[name].stat().st_size for name in order),
        "archive_bytes": total,
        "bpb": total * 8.0 / sum(by_name[name].stat().st_size for name in order),
        "order": order,
    }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
