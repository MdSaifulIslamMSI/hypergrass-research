"""Run the scale profile on an independent pangenome.

Reuses the measurement and verification code from
``benchmarks/pangenome/run.py`` so the independent datasets are tested with the
same container, archive checks, and query probes. The only difference is the
manifest and output location.

Usage:
    python -B benchmarks/independent/run.py --dataset salmonella --run-id salmonella-scale-1
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import importlib.util
import json
from pathlib import Path
import subprocess
import sys


PROJECT = Path(__file__).resolve().parent.parent.parent
DATASETS = Path(__file__).resolve().parent / "data"
RUNS = Path(__file__).resolve().parent / "runs"

SPEC = importlib.util.spec_from_file_location(
    "pangenome_run", PROJECT / "benchmarks/pangenome/run.py"
)
pangenome_run = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(pangenome_run)

PROFILE_IDS = ("hypergrass", "hypergrass-pangenome", "gzip-9", "xz-9")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--run-id", default=None)
    parser.add_argument("--profiles", default=",".join(PROFILE_IDS))
    parser.add_argument("--timeout", type=float, default=3600.0)
    args = parser.parse_args(argv)

    dataset_dir = DATASETS / args.dataset
    manifest_path = dataset_dir / "manifest.json"
    if not manifest_path.exists():
        raise SystemExit(f"missing dataset manifest: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    fasta = dataset_dir / manifest["fasta"]["path"]
    if not fasta.exists():
        raise SystemExit(f"missing pangenome FASTA: {fasta}")
    if fasta.stat().st_size != manifest["fasta"]["size_bytes"]:
        raise SystemExit("pangenome size does not match manifest")
    digest = pangenome_run.sha256_file(fasta)
    if digest != manifest["fasta"]["sha256"]:
        raise SystemExit("pangenome SHA-256 does not match manifest")
    if fasta.stat().st_size < 100_000_000:
        raise SystemExit("independent scale test requires a 100 MB+ pangenome")

    image = pangenome_run.image_id(PROJECT)
    inspect = subprocess.run(
        ["docker", "image", "inspect", image, "--format", "{{.Os}}/{{.Architecture}}"],
        capture_output=True, text=True, timeout=60,
    )
    if inspect.returncode != 0 or inspect.stdout.strip() != "linux/amd64":
        raise SystemExit("locked image is not available or not linux/amd64")
    if not (PROJECT / "bin/hypergrass/HyperGraSS.class").exists():
        raise SystemExit("compiled classes not found; run build.ps1 or build.sh")

    wanted = set(args.profiles.split(","))
    profiles = [p for p in pangenome_run.PROFILES if p["id"] in wanted]
    missing = wanted - {p["id"] for p in profiles}
    if missing:
        raise SystemExit(f"unknown profiles: {', '.join(sorted(missing))}")
    if not profiles:
        raise SystemExit("no profiles selected")

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_id = args.run_id or f"{args.dataset}-{stamp}"
    run_dir = RUNS / args.dataset / run_id
    if run_dir.exists():
        raise SystemExit(f"run directory already exists: {run_dir}")
    run_dir.mkdir(parents=True)

    report = {
        "schema_version": 1,
        "kind": "independent-scale-test",
        "run_id": run_id,
        "dataset": args.dataset,
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "image_id": image,
        "manifest": manifest,
        "fasta": str(fasta.relative_to(PROJECT)),
        "profiles": [p["id"] for p in profiles],
        "results": [],
        "status": "running",
    }
    pangenome_run.write_json(run_dir / "report.json", report)
    print(
        f"run {run_id}: {len(profiles)} profiles on "
        f"{manifest['fasta']['size_bytes']:,} bytes"
    )

    failures = 0
    for profile in profiles:
        print(f"  {profile['id']} ...", flush=True)
        try:
            entry = pangenome_run.run_profile(
                PROJECT, image, run_dir, fasta, digest,
                manifest["fasta"]["size_bytes"], profile, args.timeout,
            )
        except Exception as error:
            entry = {
                "profile": profile["id"],
                "status": "error",
                "error": f"{type(error).__name__}: {error}",
            }
        report["results"].append(entry)
        if entry["status"] != "success":
            failures += 1
            print(f"    FAILED: {entry.get('error')}", flush=True)
        else:
            print(
                f"    {entry['artifact_bytes']:,} B, {entry['ratio']:.2f}x, "
                f"{entry['wall_seconds']:.1f}s, "
                f"{entry['peak_rss_bytes'] / 1048576:.0f} MiB RSS",
                flush=True,
            )
        pangenome_run.write_json(run_dir / "report.json", report)

    report["status"] = "success" if failures == 0 else "partial"
    report["failures"] = failures
    report["finished_utc"] = datetime.now(timezone.utc).isoformat()
    pangenome_run.write_json(run_dir / "report.json", report)
    print(f"report: {run_dir / 'report.json'}")
    return 0 if failures == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
