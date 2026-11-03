"""Execute the broad-experiment matrix through the frozen specialized harness.

Each run gets its own directory under ``benchmarks/broad/runs`` and one
``report.json`` written by ``benchmarks/specialized/run.py``. The wrapper only
schedules runs and records the cross-run index; it does not re-implement any
compression or verification logic.

Usage:
    python -B benchmarks/broad/run.py --dry-run
    python -B benchmarks/broad/run.py --panel cross-taxon
    python -B benchmarks/broad/run.py --panel diversity --dataset staph
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parent.parent
MATRIX = ROOT / "matrix.json"
RUNS = ROOT / "runs"
SPECIALIZED = PROJECT / "benchmarks/specialized/run.py"


def matches(run: dict, panels: set[str], datasets: set[str],
            strata: set[str]) -> bool:
    if panels and run["panel"] not in panels:
        return False
    if datasets and run["dataset"] not in datasets:
        return False
    if strata and run["stratum"] not in strata:
        return False
    return True


def summarize_report(report: dict) -> dict:
    results = []
    for entry in report.get("results", []):
        results.append(
            {
                "profile": entry.get("profile"),
                "status": entry.get("status"),
                "artifact_bytes": entry.get("artifact_bytes"),
                "bits_per_base": entry.get("bits_per_base"),
                "wall_seconds": entry.get("wall_seconds"),
                "cpu_seconds": entry.get("cpu_seconds"),
                "peak_rss_bytes": entry.get("peak_rss_bytes"),
                "roundtrip_ok": entry.get("roundtrip_ok"),
                "error": entry.get("error"),
            }
        )
    return {
        "status": report.get("status"),
        "failures": report.get("failures"),
        "input": report.get("input", {}).get("sequence", {}),
        "results": results,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", default=str(MATRIX))
    parser.add_argument("--runs-dir", default=str(RUNS))
    parser.add_argument("--panel", default=None)
    parser.add_argument("--dataset", default=None)
    parser.add_argument("--stratum", default=None)
    parser.add_argument("--profiles", default=None)
    parser.add_argument("--timeout", type=float, default=7200.0)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--resume", action="store_true", default=True)
    parser.add_argument("--rerun", action="store_true")
    args = parser.parse_args(argv)
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")

    matrix_path = Path(args.matrix)
    matrix = json.loads(matrix_path.read_text(encoding="utf-8"))
    runs_dir = Path(args.runs_dir)
    logs_dir = runs_dir / "logs"
    panels = set(args.panel.split(",")) if args.panel else set()
    datasets = set(args.dataset.split(",")) if args.dataset else set()
    strata = set(args.stratum.split(",")) if args.stratum else set()

    selected = [
        run for run in matrix["runs"]
        if matches(run, panels, datasets, strata)
    ]
    if not selected:
        raise SystemExit("no matrix runs matched the filters")
    if args.profiles:
        for run in selected:
            run["profiles"] = [value for value in args.profiles.split(",") if value]

    if not args.dry_run:
        logs_dir.mkdir(parents=True, exist_ok=True)

    master_path = runs_dir / "master.json"
    master = {"schema_version": 1, "matrix": str(matrix_path), "runs": []}
    if master_path.exists():
        master = json.loads(master_path.read_text(encoding="utf-8"))
    known = {entry["run_id"]: entry for entry in master.get("runs", [])}

    for run in selected:
        run_id = run["run_id"]
        run_dir = runs_dir / run_id
        report_path = run_dir / "report.json"
        if report_path.exists() and not args.rerun:
            report = json.loads(report_path.read_text(encoding="utf-8"))
            if report.get("status") == "running":
                raise SystemExit(
                    f"run still in progress or stale: {run_id}; wait for it to finish "
                    "or remove the run directory and rerun"
                )
            print(f"skip {run_id}: existing report status={report.get('status')}")
            known[run_id] = {
                "run_id": run_id,
                "panel": run["panel"],
                "dataset": run["dataset"],
                "organism": run["organism"],
                "stratum": run["stratum"],
                "source": run["source"],
                "diversity": run["diversity"],
                "report": str(report_path.relative_to(PROJECT).as_posix()),
                **summarize_report(report),
            }
            continue
        if run_dir.exists() and not report_path.exists():
            raise SystemExit(
                f"incomplete run directory exists without report: {run_dir}; "
                "remove it or use a new run id"
            )

        command = [
            sys.executable, "-B", str(SPECIALIZED),
            "--source", str(PROJECT / run["source"]),
            "--profiles", ",".join(run["profiles"]),
            "--run-id", run_id,
            "--runs-dir", str(runs_dir),
            "--timeout", str(args.timeout),
        ]
        print(
            f"start {run_id}: {run['dataset']}/{run['stratum']} "
            f"{run['diversity']['record_count']} records, "
            f"{run['diversity']['total_bases']:,} bases, "
            f"d={run['diversity']['mean_pairwise_distance']:.5f}"
        )
        print(f"  {' '.join(command)}")
        if args.dry_run:
            continue
        if args.rerun and run_dir.exists():
            # The specialized harness refuses to write into an existing run directory, so a
            # rerun has to start from an empty one. Without this, --rerun silently produces
            # no report and the stale one is recorded as if it were fresh.
            shutil.rmtree(run_dir)

        log_path = logs_dir / f"{run_id}.log"
        started = datetime.now(timezone.utc)
        with log_path.open("w", encoding="utf-8") as log:
            process = subprocess.Popen(
                command, cwd=PROJECT, stdout=log, stderr=subprocess.STDOUT, text=True,
            )
            returncode = process.wait()
        if not report_path.exists():
            known[run_id] = {
                "run_id": run_id,
                "panel": run["panel"],
                "dataset": run["dataset"],
                "organism": run["organism"],
                "stratum": run["stratum"],
                "source": run["source"],
                "diversity": run["diversity"],
                "status": f"wrapper-error-{returncode}",
                "log": str(log_path.relative_to(PROJECT).as_posix()),
            }
            print(f"  FAILED: no report; see {log_path}")
        else:
            report = json.loads(report_path.read_text(encoding="utf-8"))
            known[run_id] = {
                "run_id": run_id,
                "panel": run["panel"],
                "dataset": run["dataset"],
                "organism": run["organism"],
                "stratum": run["stratum"],
                "source": run["source"],
                "diversity": run["diversity"],
                "report": str(report_path.relative_to(PROJECT).as_posix()),
                "log": str(log_path.relative_to(PROJECT).as_posix()),
                "started_utc": started.isoformat(),
                "finished_utc": datetime.now(timezone.utc).isoformat(),
                **summarize_report(report),
            }
            print(
                f"  done: status={report.get('status')} "
                f"failures={report.get('failures')} ({log_path.name})"
            )
        master["runs"] = sorted(known.values(), key=lambda entry: entry["run_id"])
        master["updated_utc"] = datetime.now(timezone.utc).isoformat()
        master_path.parent.mkdir(parents=True, exist_ok=True)
        master_path.write_text(json.dumps(master, indent=2) + "\n", encoding="utf-8")

    if not args.dry_run:
        master["runs"] = sorted(known.values(), key=lambda entry: entry["run_id"])
        master["updated_utc"] = datetime.now(timezone.utc).isoformat()
        master_path.write_text(json.dumps(master, indent=2) + "\n", encoding="utf-8")
        print(f"master: {master_path} ({len(master['runs'])} runs)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
