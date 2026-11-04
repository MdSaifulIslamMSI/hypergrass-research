#!/usr/bin/env python3
"""Sweep AGC parameters on one leave-one-out target and report the marginal size."""
from __future__ import annotations

import argparse
import json
import subprocess
import time
from pathlib import Path

CONFIGS = [
    ("default", ["-a"], ["-a"]),
    ("fallback", ["-a", "-f", "0.05"], ["-a", "-f", "0.05"]),
    ("min15", ["-a", "-l", "15"], ["-a"]),
    ("min15-fallback", ["-a", "-l", "15", "-f", "0.05"], ["-a", "-f", "0.05"]),
    ("seg50k", ["-a", "-s", "50000"], ["-a"]),
    ("seg100k", ["-a", "-s", "100000"], ["-a"]),
    ("seg200k", ["-a", "-s", "200000"], ["-a"]),
    ("seg500k", ["-a", "-s", "500000"], ["-a"]),
    ("seg1m", ["-a", "-s", "1000000"], ["-a"]),
    ("seg1m-min15", ["-a", "-s", "1000000", "-l", "15"], ["-a"]),
    ("k25-seg200k", ["-a", "-k", "25", "-s", "200000"], ["-a"]),
    ("plain-seg200k", ["-s", "200000"], []),
    ("seg200k-min15-fallback", ["-a", "-s", "200000", "-l", "15", "-f", "0.05"], ["-a", "-f", "0.05"]),
]


def run(cmd: list[str]) -> float:
    t = time.perf_counter()
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise SystemExit("failed: %s\n%s\n%s" % (" ".join(map(str, cmd)), p.stdout, p.stderr))
    return time.perf_counter() - t


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta-dir", type=Path, default=Path("prototype/rcmdl/work/fasta"))
    parser.add_argument("--out", type=Path, default=Path("prototype/rcmdl/work/agc-sweep"))
    parser.add_argument("--agc", type=Path, default=Path("prototype/rcmdl/vendor/agc/agc.exe"))
    parser.add_argument("--target", type=int, default=0)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--configs", default="", help="space or comma separated config names")
    args = parser.parse_args()

    target = args.fasta_dir / ("%02d.fa" % args.target)
    others = [p for p in sorted(args.fasta_dir.glob("*.fa")) if p != target]
    if not others:
        raise SystemExit("no panel files in %s" % args.fasta_dir)
    args.out.mkdir(parents=True, exist_ok=True)
    listfile = args.out / ("panel-%02d.list" % args.target)
    listfile.write_text("\n".join(str(p.resolve()) for p in others[1:]) + "\n")
    wanted = {x for x in args.configs.replace(",", " ").split() if x} or None

    rows = []
    for name, opts, append_opts in CONFIGS:
        if wanted is not None and name not in wanted:
            continue
        base = args.out / ("base-%s.agc" % name)
        full = args.out / ("full-%s.agc" % name)
        t_create = run(
            [str(args.agc), "create", *opts, "-t", str(args.threads), "-i", str(listfile), "-o", str(base), str(others[0])]
        )
        t_append = run(
            [str(args.agc), "append", *append_opts, "-t", str(args.threads), "-o", str(full), str(base), str(target)]
        )
        delta = full.stat().st_size - base.stat().st_size
        row = {
            "config": name,
            "opts": opts,
            "base_bytes": base.stat().st_size,
            "full_bytes": full.stat().st_size,
            "marginal_bytes": delta,
            "bpb": delta * 8.0 / target.stat().st_size,
            "create_seconds": t_create,
            "append_seconds": t_append,
        }
        rows.append(row)
        print(row, flush=True)
        base.unlink()
        full.unlink()
    (args.out / ("results-%02d.json" % args.target)).write_text(json.dumps(rows, indent=2) + "\n")


if __name__ == "__main__":
    main()
