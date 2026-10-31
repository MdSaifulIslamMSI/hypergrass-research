"""Anchor-bank and block-size ablation on the 121.8 MB pangenome.

Tests the hypothesis recorded in docs/BENCHMARKS.md: the default 4 MiB anchor
bank cannot represent a full ~5 Mbp record, and 2 MiB blocks cannot exploit the
same long-range matches as xz's 64 MiB dictionary. Each configuration runs in a
fresh locked-container process; the harness verifies byte-exact restoration.

Usage:
    python -B benchmarks/ablation/run.py --run-id anchor-ablation-1
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "pangenome"))
from run import (  # noqa: E402
    image_id,
    parse_fasta,
    run_operation,
    sha256_file,
    to_container,
    write_json,
)

PROJECT = Path(__file__).resolve().parent.parent.parent
PANGENOME = PROJECT / "benchmarks/pangenome"
BASELINE_REPORT = PANGENOME / "runs/pangenome-scale-2/report.json"


def baseline_bytes() -> tuple[int | None, str | None]:
    if not BASELINE_REPORT.exists():
        return None, None
    report = json.loads(BASELINE_REPORT.read_text(encoding="utf-8"))
    for entry in report.get("results", []):
        if entry.get("profile") == "xz-9" and entry.get("status") == "success":
            return entry["artifact_bytes"], str(BASELINE_REPORT.relative_to(PROJECT))
    return None, None


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-id", default=None)
    parser.add_argument("--banks", default="4,8,16")
    parser.add_argument("--blocks", default="2,4,8")
    parser.add_argument("--timeout", type=float, default=1800.0)
    args = parser.parse_args(argv)

    manifest = json.loads((PANGENOME / "manifest.json").read_text(encoding="utf-8"))
    fasta = PANGENOME / manifest["fasta"]["path"]
    if not fasta.exists():
        raise SystemExit(f"missing pangenome FASTA: {fasta}")
    if fasta.stat().st_size != manifest["fasta"]["size_bytes"]:
        raise SystemExit("pangenome size does not match manifest")
    digest = sha256_file(fasta)
    if digest != manifest["fasta"]["sha256"]:
        raise SystemExit("pangenome SHA-256 does not match manifest")
    if fasta.stat().st_size < 100_000_000:
        raise SystemExit("ablation requires the 100 MB+ pangenome")

    if not (PROJECT / "bin/hypergrass/HyperGraSS.class").exists():
        raise SystemExit("compiled classes not found; run build.ps1 or build.sh")
    if not (PROJECT / "bin/ablation/AnchorBankAblation.class").exists():
        raise SystemExit(
            "ablation harness not compiled; run:\n"
            "  javac -cp bin -d bin benchmarks/ablation/AnchorBankAblation.java"
        )

    image = image_id(PROJECT)
    inspect = subprocess.run(
        ["docker", "image", "inspect", image, "--format", "{{.Os}}/{{.Architecture}}"],
        capture_output=True, text=True, timeout=60,
    )
    if inspect.returncode != 0 or inspect.stdout.strip() != "linux/amd64":
        raise SystemExit("locked image is not available or not linux/amd64")

    banks = [int(value) for value in args.banks.split(",")]
    blocks = [int(value) for value in args.blocks.split(",")]
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_id = args.run_id or f"anchor-ablation-{stamp}"
    run_dir = PROJECT / "benchmarks/ablation/runs" / run_id
    if run_dir.exists():
        raise SystemExit(f"run directory already exists: {run_dir}")
    run_dir.mkdir(parents=True)

    headers, _, lengths = parse_fasta(fasta)
    baseline, baseline_source = baseline_bytes()
    report = {
        "schema_version": 1,
        "kind": "anchor-bank-ablation",
        "run_id": run_id,
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "image_id": image,
        "input": {
            "path": str(fasta.relative_to(PROJECT)),
            "size_bytes": fasta.stat().st_size,
            "sha256": digest,
            "records": len(headers),
            "bases": sum(lengths),
        },
        "baseline": {"profile": "xz-9", "archive_bytes": baseline, "source": baseline_source},
        "configs": {"banks_mb": banks, "blocks_mb": blocks},
        "results": [],
        "status": "running",
    }
    write_json(run_dir / "report.json", report)
    print(f"run {run_id}: {len(banks) * len(blocks)} configurations")

    for bank_mb in banks:
        for block_mb in blocks:
            name = f"bank{bank_mb}-block{block_mb}"
            work = run_dir / "work" / name
            work.mkdir(parents=True, exist_ok=True)
            archive = work / "archive.hgs"
            restored = work / "restored.fasta"
            command = [
                "java", "-Xmx1536m", "-XX:ActiveProcessorCount=1", "-cp", "/work/bin",
                "ablation.AnchorBankAblation", to_container(PROJECT, fasta),
                str(bank_mb), str(block_mb), to_container(PROJECT, archive),
                to_container(PROJECT, restored),
            ]
            print(f"  {name} ...", flush=True)
            operation = run_operation(
                PROJECT, image, run_dir, name, command, work, args.timeout,
                memory="2g",
            )
            entry = {
                "bank_mb": bank_mb,
                "block_mb": block_mb,
                "status": operation["status"],
                "wall_seconds": operation.get("wall_seconds"),
                "cpu_seconds": operation.get("cpu_seconds"),
                "peak_rss_bytes": operation.get("peak_rss_bytes"),
                "error": operation.get("error"),
            }
            if operation["status"] == "success":
                metrics_path = run_dir / "ops" / name / "stdout.log"
                lines = [
                    line for line in metrics_path.read_text(encoding="utf-8").splitlines()
                    if line.strip()
                ]
                if not lines:
                    entry["status"] = "error"
                    entry["error"] = "harness produced no metrics"
                else:
                    try:
                        entry.update(json.loads(lines[-1]))
                    except json.JSONDecodeError as error:
                        entry["status"] = "error"
                        entry["error"] = f"invalid harness output: {error}"
                    if entry.get("exact") is not True:
                        entry["status"] = "error"
                        entry["error"] = "round trip was not byte-exact"
            if entry.get("archive_bytes") and baseline:
                entry["delta_vs_xz_bytes"] = entry["archive_bytes"] - baseline
                entry["delta_vs_xz_percent"] = (
                    (entry["archive_bytes"] - baseline) / baseline * 100.0
                )
                entry["beats_xz"] = entry["archive_bytes"] < baseline
            report["results"].append(entry)
            write_json(run_dir / "report.json", report)
            if entry["status"] == "success":
                print(
                    f"    {entry['archive_bytes']:,} B, "
                    f"{entry['delta_vs_xz_percent']:+.2f}% vs xz, "
                    f"{entry['wall_seconds']:.1f}s wall, "
                    f"{entry['peak_rss_bytes'] / 1048576:.0f} MiB RSS"
                )
            else:
                print(f"    {entry['status']}: {entry.get('error')}")

    succeeded = [entry for entry in report["results"] if entry["status"] == "success"]
    report["status"] = "success" if len(succeeded) == len(report["results"]) else "partial"
    report["finished_utc"] = datetime.now(timezone.utc).isoformat()
    write_json(run_dir / "report.json", report)

    if succeeded:
        winner = min(succeeded, key=lambda entry: entry["archive_bytes"])
        print(
            f"best: bank {winner['bank_mb']} MiB, block {winner['block_mb']} MiB, "
            f"{winner['archive_bytes']:,} B "
            f"({winner['delta_vs_xz_percent']:+.2f}% vs xz)"
        )
    print(f"report: {run_dir / 'report.json'}")
    return 0 if report["status"] == "success" else 1


if __name__ == "__main__":
    raise SystemExit(main())
