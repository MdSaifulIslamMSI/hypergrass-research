#!/usr/bin/env python3
"""Run the RC-MDL prototype: self-contained (S) and leave-one-out (R) modes.

Every archive is decompressed and compared byte-for-byte against the target via
SHA-256 before its numbers are accepted.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def invoke(cmd: list[str]) -> tuple[dict, float, str]:
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
    return json.loads(lines[-1]), elapsed, proc.stderr


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--genomes", type=Path, default=Path("prototype/rcmdl/work/genomes"))
    parser.add_argument("--out", type=Path, required=True)
    default_java = Path.home() / "tools/jdk-27/bin/java.exe"
    parser.add_argument(
        "--java",
        default=str(default_java) if default_java.exists() else "java",
        help="Java executable (defaults to the pinned JDK 27 on Windows, else PATH)",
    )
    parser.add_argument("--heap", default="1500m")
    parser.add_argument("--mode", choices=["loo", "self", "both"], default="both")
    parser.add_argument("--k", type=int, default=16)
    parser.add_argument("--step", type=int, default=8)
    parser.add_argument("--panel-limit", type=int, default=0)
    parser.add_argument("--limit", type=int, default=0, help="only the first N targets")
    parser.add_argument("--no-approx", action="store_true")
    args = parser.parse_args()

    genomes = sorted(args.genomes.glob("*.seq"))
    if args.limit > 0:
        genomes = genomes[: args.limit]
    if not genomes:
        raise SystemExit("no genomes found in %s" % args.genomes)
    args.out.mkdir(parents=True, exist_ok=True)
    classpath = Path("prototype/rcmdl/bin")
    base = [args.java, "-Xmx" + args.heap, "-cp", str(classpath), "Rcmdl"]

    records = []
    for index, target in enumerate(genomes):
        hashes = {"target": sha256(target)}
        result = {
            "target": target.name,
            "target_bytes": target.stat().st_size,
            "target_sha256": hashes["target"],
        }
        modes = []
        if args.mode in ("self", "both"):
            modes.append("self")
        if args.mode in ("loo", "both") and len(genomes) > 1:
            modes.append("loo")

        for mode in modes:
            archive = args.out / ("%02d-%s-%s.rcmd" % (index, target.stem, mode))
            restored = args.out / ("%02d-%s-%s.out" % (index, target.stem, mode))
            compress = base + [
                "compress",
                "--target",
                str(target),
                "--out",
                str(archive),
                "--k",
                str(args.k),
                "--step",
                str(args.step),
            ]
            if args.no_approx:
                compress.append("--no-approx")
            if mode == "loo":
                compress += ["--panel", str(args.genomes), "--exclude", target.name]
                if args.panel_limit > 0:
                    compress += ["--panel-limit", str(args.panel_limit)]
            cstats, cwall, _ = invoke(compress)
            decompress = base + [
                "decompress",
                "--archive",
                str(archive),
                "--out",
                str(restored),
            ]
            if mode == "loo":
                decompress += ["--panel", str(args.genomes)]
            dstats, dwall, _ = invoke(decompress)
            restored_hash = sha256(restored)
            if restored_hash != hashes["target"]:
                raise RuntimeError("round-trip mismatch for %s (%s)" % (target, mode))
            if restored.stat().st_size != target.stat().st_size:
                raise RuntimeError("length mismatch for %s (%s)" % (target, mode))
            result[mode] = {
                "archive_bytes": archive.stat().st_size,
                "bpb": archive.stat().st_size * 8.0 / target.stat().st_size,
                "panel_genomes": cstats["panel_genomes"],
                "panel_bases": cstats["panel_bases"],
                "panel_match_bases": cstats["panel_match_bases"],
                "self_match_bases": cstats["self_match_bases"],
                "literal_bases": cstats["literal_bases"],
                "panel_tokens": cstats["panel_tokens"],
                "self_tokens": cstats["self_tokens"],
                "literal_runs": cstats["literal_runs"],
                "mismatches": cstats["mismatches"],
                "compress_seconds": cstats["compress_seconds"],
                "compress_wall_seconds": cwall,
                "decompress_seconds": dstats["decompress_seconds"],
                "decompress_wall_seconds": dwall,
                "verified": True,
            }
            restored.unlink()
            print(
                "%-3d %-8s %-4s %9d  %7.4f bpb  %6.2fs  verified"
                % (
                    index,
                    target.stem,
                    mode,
                    archive.stat().st_size,
                    result[mode]["bpb"],
                    cwall,
                ),
                flush=True,
            )
        records.append(result)
        (args.out / "results.json").write_text(
            json.dumps(records, indent=2) + "\n", encoding="utf-8"
        )

    summary = {"genomes": len(records), "k": args.k, "step": args.step, "modes": {}}
    for mode in ("self", "loo"):
        rows = [r[mode] for r in records if mode in r]
        if not rows:
            continue
        target_bases = sum(r["target_bytes"] for r in records if mode in r)
        archive_bytes = sum(r["archive_bytes"] for r in rows)
        summary["modes"][mode] = {
            "targets": len(rows),
            "target_bases": target_bases,
            "archive_bytes": archive_bytes,
            "bpb": archive_bytes * 8.0 / target_bases if target_bases else 0.0,
            "mean_bpb": sum(r["bpb"] for r in rows) / len(rows),
            "panel_match_fraction": sum(r["panel_match_bases"] for r in rows) / target_bases
            if target_bases
            else 0.0,
            "literal_fraction": sum(r["literal_bases"] for r in rows) / target_bases
            if target_bases
            else 0.0,
            "mismatch_bases": sum(r["mismatches"] for r in rows),
            "compress_seconds": sum(r["compress_seconds"] for r in rows),
            "decompress_seconds": sum(r["decompress_seconds"] for r in rows),
        }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    if any(
        mode.get("verified") is not True
        for record in records
        for mode in (record.get("self"), record.get("loo"))
        if mode is not None
    ):
        sys.exit("unverified records present")


if __name__ == "__main__":
    main()
