#!/usr/bin/env python3
"""Parameter and panel-size ablations for the RC-MDL prototype, with round-trip checks."""
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
        raise RuntimeError("failed: %s\n%s\n%s" % (" ".join(cmd), proc.stdout, proc.stderr))
    lines = [line for line in proc.stdout.splitlines() if line.startswith("{")]
    if not lines:
        raise RuntimeError("no stats: %s" % " ".join(cmd))
    return json.loads(lines[-1]), elapsed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=["panel-size", "k", "step", "approx"], required=True)
    parser.add_argument("--target", type=int, default=0)
    parser.add_argument("--genomes", type=Path, default=Path("prototype/rcmdl/work/genomes"))
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--java", default=str(Path.home() / "tools/jdk-27/bin/java.exe"))
    parser.add_argument("--classpath", default="prototype/rcmdl/bin")
    parser.add_argument("--heap", default="1400m")
    args = parser.parse_args()

    target = args.genomes / ("%02d.seq" % args.target)
    args.out.mkdir(parents=True, exist_ok=True)
    base = [args.java, "-Xmx" + args.heap, "-cp", args.classpath, "Rcmdl"]
    if args.kind == "panel-size":
        configs = [("panel-%02d" % n, {"panel_limit": n}) for n in (1, 2, 4, 8, 16, 23)]
    elif args.kind == "k":
        configs = [("k-%02d" % k, {"k": k}) for k in (12, 16, 20, 24)]
    elif args.kind == "step":
        configs = [("step-%02d" % s, {"step": s}) for s in (4, 8, 16, 32)]
    else:
        configs = [("approx", {}), ("no-approx", {"no_approx": True})]

    records = []
    for name, cfg in configs:
        archive = args.out / (name + ".rcmd")
        restored = args.out / (name + ".out")
        cmd = base + [
            "compress",
            "--target",
            str(target),
            "--panel",
            str(args.genomes),
            "--exclude",
            target.name,
            "--out",
            str(archive),
            "--k",
            str(cfg.get("k", 16)),
            "--step",
            str(cfg.get("step", 8)),
        ]
        if "panel_limit" in cfg:
            cmd += ["--panel-limit", str(cfg["panel_limit"])]
        if cfg.get("no_approx"):
            cmd.append("--no-approx")
        stats, wall = invoke(cmd)
        dec, dwall = invoke(
            base
            + ["decompress", "--archive", str(archive), "--panel", str(args.genomes), "--out", str(restored)]
        )
        if sha256(restored) != sha256(target):
            raise RuntimeError("round-trip mismatch: " + name)
        record = {
            "config": name,
            "archive_bytes": archive.stat().st_size,
            "bpb": archive.stat().st_size * 8.0 / target.stat().st_size,
            "panel_genomes": stats["panel_genomes"],
            "panel_bases": stats["panel_bases"],
            "panel_match_bases": stats["panel_match_bases"],
            "literal_bases": stats["literal_bases"],
            "panel_tokens": stats["panel_tokens"],
            "self_tokens": stats["self_tokens"],
            "mismatches": stats["mismatches"],
            "compress_seconds": stats["compress_seconds"],
            "wall_seconds": wall,
            "decompress_seconds": dec["decompress_seconds"],
            "verified": True,
        }
        records.append(record)
        print(record, flush=True)
        (args.out / "results.json").write_text(json.dumps(records, indent=2) + "\n")
        restored.unlink()


if __name__ == "__main__":
    main()
