#!/usr/bin/env python3
"""Permanent Kaggle-Only End-to-End Orchestrator for HyperGraSS / EIMDL.

Replaces all local unit/integration/benchmark runs:
  1. Syntax-compiles Java locally only for preflight typo detection (--release 21 -Xlint:all)
     NEVER executing any Java tests or compression benchmarks on the local machine.
  2. Stages and pushes `mdsaifulislammsi/hypergrass-eimdl-src` to Kaggle Datasets.
  3. Waits for Kaggle dataset extraction (`ready`).
  4. Pushes selected Kaggle kernel(s):
       --target shootout  (default, CPU-only, uses 0 GPU/TPU quota: all 7 Verify* suites + 15+ world compressors incl. GraSS)
       --target all       (shootout + gpu_largescale + tpu_scale)
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

REPO_ROOT = Path(__file__).resolve().parents[2]
DATASET_SLUG = "mdsaifulislammsi/hypergrass-eimdl-src"
KERNELS = {
    "shootout": REPO_ROOT / "benchmarks" / "kaggle" / "world_shootout",
    "gpu": REPO_ROOT / "benchmarks" / "kaggle" / "gpu_largescale",
    "tpu": REPO_ROOT / "benchmarks" / "kaggle" / "tpu_scale",
}


def run(cmd: list[str], *, check: bool = True) -> subprocess.CompletedProcess[str]:
    print("+", " ".join(cmd), flush=True)
    return subprocess.run(cmd, check=check, text=True)


def stage_dataset(stage_dir: Path) -> None:
    if stage_dir.exists():
        shutil.rmtree(stage_dir)
    stage_dir.mkdir(parents=True)
    meta = {
        "id": DATASET_SLUG,
        "title": "HyperGraSS EIMDL Src",
        "licenses": [{"name": "apache-2.0"}],
        "isPrivate": True,
    }
    (stage_dir / "dataset-metadata.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")

    for helper in ("kaggle_common.py", "world_compressors.py"):
        shutil.copy2(REPO_ROOT / "benchmarks" / "kaggle" / helper, stage_dir / helper)

    sub = stage_dir / "hypergrass_eimdl_src"
    (sub / "src" / "hypergrass").mkdir(parents=True)
    (sub / "benchmarks" / "corpus").mkdir(parents=True)
    (sub / "benchmarks" / "pangenome").mkdir(parents=True)
    (sub / "benchmarks" / "independent" / "data" / "salmonella").mkdir(parents=True)
    (sub / "benchmarks" / "independent" / "data" / "klebsiella").mkdir(parents=True)
    (sub / "grass_fixtures").mkdir(parents=True)

    for rel in (
        "build.sh",
        "pom.xml",
        "benchmarks/corpus/sources.lock.json",
        "benchmarks/pangenome/manifest.json",
        "benchmarks/independent/data/salmonella/manifest.json",
        "benchmarks/independent/data/klebsiella/manifest.json",
    ):
        shutil.copy2(REPO_ROOT / rel, sub / rel)

    for jf in sorted((REPO_ROOT / "src" / "hypergrass").glob("*.java")):
        shutil.copy2(jf, sub / "src" / "hypergrass" / jf.name)

    grass_dir = Path.home() / ".gemini/antigravity/brain/92f6b277-dea3-42a9-b123-90c65fe25a1c/scratch/GraSS/src"
    if grass_dir.is_dir():
        for fn in ("In.fa", "BuEb", "GCF_001884535.1.fna"):
            src = grass_dir / fn
            if src.exists():
                shutil.copy2(src, sub / "grass_fixtures" / fn)


def main() -> None:
    ap = argparse.ArgumentParser(description="Push and verify HyperGraSS/EIMDL end-to-end strictly on Kaggle.")
    ap.add_argument("--target", choices=["shootout", "gpu", "tpu", "all"], default="shootout")
    ap.add_argument("--message", "-m", default="E2E world shootout + 7-suite cloud verification")
    ap.add_argument("--skip-dataset", action="store_true")
    args = ap.parse_args()

    if not args.skip_dataset:
        with tempfile.TemporaryDirectory(prefix="hgs_kstage_") as tmp:
            st = Path(tmp) / "stage"
            stage_dataset(st)
            run(["kaggle", "datasets", "version", "-p", str(st), "-r", "zip", "-m", args.message])
        for _ in range(24):
            res = subprocess.run(
                ["kaggle", "datasets", "status", DATASET_SLUG], capture_output=True, text=True, check=False
            )
            if "ready" in (res.stdout + res.stderr).lower():
                print(f"Dataset {DATASET_SLUG} ready.", flush=True)
                break
            time.sleep(5)

    targets = ["shootout", "gpu", "tpu"] if args.target == "all" else [args.target]
    for t in targets:
        run(["kaggle", "kernels", "push", "-p", str(KERNELS[t])])

    print("All requested Kaggle E2E runs pushed. No local tests were executed.")


if __name__ == "__main__":
    main()
