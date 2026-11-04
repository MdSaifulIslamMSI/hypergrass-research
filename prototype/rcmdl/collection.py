#!/usr/bin/env python3
"""Progressive collection archive for the RC-MDL prototype.

Genome 0 is coded self-contained. Genome i is coded against genomes 0..i-1 only.
The sum of the per-genome archives is a self-contained collection archive that can
be compared directly with a collection codec that stores every genome.
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--genomes", type=Path, default=Path("prototype/rcmdl/work/genomes"))
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--java", default=str(Path.home() / "tools/jdk-27/bin/java.exe"))
    parser.add_argument("--classpath", default="prototype/rcmdl/bin")
    parser.add_argument("--heap", default="1400m")
    parser.add_argument("--order", choices=["name", "size-asc", "size-desc"], default="name")
    parser.add_argument("--order-file", type=Path)
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument(
        "--fixed-anchors",
        type=int,
        default=0,
        help="code the first N genomes progressively, then every other genome "
        "against those N anchors only (0 = full progressive collection)",
    )
    parser.add_argument("--k", type=int, default=16)
    parser.add_argument("--step", type=int, default=8)
    parser.add_argument("--no-approx", action="store_true")
    args = parser.parse_args()

    genomes = sorted(args.genomes.glob("*.seq"))
    if args.order == "size-asc":
        genomes.sort(key=lambda p: p.stat().st_size)
    elif args.order == "size-desc":
        genomes.sort(key=lambda p: p.stat().st_size, reverse=True)
    if args.order_file:
        wanted = [line.strip() for line in args.order_file.read_text().splitlines() if line.strip()]
        by_name = {p.name: p for p in genomes}
        missing = [name for name in wanted if name not in by_name]
        if missing:
            raise SystemExit("order file names not found: %s" % missing)
        genomes = [by_name[name] for name in wanted]
    if args.limit > 0:
        genomes = genomes[: args.limit]
    if not genomes:
        raise SystemExit("no genomes")
    args.out.mkdir(parents=True, exist_ok=True)
    base = [args.java, "-Xmx" + args.heap, "-cp", args.classpath, "Rcmdl"]
    list_path = args.out / "panel.txt"
    records = []
    total = 0
    for i, target in enumerate(genomes):
        archive = args.out / ("%02d-%s.rcmd" % (i, target.stem))
        restored = args.out / ("%02d-%s.out" % (i, target.stem))
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
        if i > 0:
            if args.fixed_anchors > 0 and i >= args.fixed_anchors:
                panel = genomes[: args.fixed_anchors]
            else:
                panel = genomes[:i]
            list_path.write_text(
                "\n".join(p.name for p in panel) + "\n", encoding="utf-8"
            )
            compress += ["--panel", str(args.genomes), "--panel-list", str(list_path)]
        stats, wall = invoke(compress)
        decompress = base + ["decompress", "--archive", str(archive), "--out", str(restored)]
        if i > 0:
            decompress += ["--panel", str(args.genomes)]
        dstats, dwall = invoke(decompress)
        if sha256(restored) != sha256(target):
            raise RuntimeError("round-trip mismatch for %s" % target)
        size = archive.stat().st_size
        total += size
        record = {
            "index": i,
            "target": target.name,
            "target_bytes": target.stat().st_size,
            "archive_bytes": size,
            "bpb": size * 8.0 / target.stat().st_size,
            "cumulative_bytes": total,
            "panel_genomes": stats["panel_genomes"],
            "panel_bases": stats["panel_bases"],
            "panel_match_bases": stats["panel_match_bases"],
            "self_match_bases": stats["self_match_bases"],
            "literal_bases": stats["literal_bases"],
            "panel_tokens": stats["panel_tokens"],
            "self_tokens": stats["self_tokens"],
            "mismatches": stats["mismatches"],
            "compress_seconds": stats["compress_seconds"],
            "compress_wall_seconds": wall,
            "decompress_seconds": dstats["decompress_seconds"],
            "decompress_wall_seconds": dwall,
            "verified": True,
        }
        records.append(record)
        bases = sum(r["target_bytes"] for r in records)
        record["cumulative_bpb"] = total * 8.0 / bases
        restored.unlink()
        (args.out / "results.json").write_text(
            json.dumps(records, indent=2) + "\n", encoding="utf-8"
        )
        print(
            "%2d %-6s %-8s %9d  %8.4f bpb  cum %10d  %7.4f bpb  %5.1fs"
            % (
                i,
                target.stem,
                "self" if i == 0 else "panel",
                size,
                record["bpb"],
                total,
                record["cumulative_bpb"],
                wall,
            ),
            flush=True,
        )
    summary = {
        "genomes": len(records),
        "order": args.order_file.name if args.order_file else args.order,
        "fixed_anchors": args.fixed_anchors,
        "k": args.k,
        "step": args.step,
        "target_bases": sum(r["target_bytes"] for r in records),
        "archive_bytes": total,
        "bpb": total * 8.0 / sum(r["target_bytes"] for r in records),
        "first_genome_bytes": records[0]["archive_bytes"],
        "marginal_bytes": total - records[0]["archive_bytes"],
        "verified": all(r["verified"] for r in records),
    }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
