"""End-to-end query benchmark on the 121.8 MB pangenome.

The archive-size benchmark answers one question. This one asks whether
HyperGraSS is useful for a concrete workload: serve many small record slices
from a stored pangenome, with one full restore.

Every query runs in a fresh locked-container process. HyperGraSS answers from
its page index. gzip and xz must decompress the stream up to the requested
record, which is the honest cost of a non-indexed format.

Usage:
    python -B benchmarks/endtoend/run.py --run-id endtoend-1
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import statistics
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "pangenome"))
from run import (  # noqa: E402
    image_id,
    parse_fasta,
    run_operation,
    sequence_slice,
    sha256_file,
    to_container,
    write_json,
)

PROJECT = Path(__file__).resolve().parent.parent.parent
PANGENOME = PROJECT / "benchmarks/pangenome"
SCALE_RUN = PANGENOME / "runs/pangenome-scale-2"
PANGENOME_RUN = PANGENOME / "runs/pangenome-scale-3"
SCALE_REPORT = SCALE_RUN / "report.json"
QUERY_LENGTH = 2_000

PROFILES = (
    {
        "id": "hypergrass",
        "kind": "hgs",
        "heap": "1g",
        "archive": "archive.hgs",
        "memory": "2g",
        "scale_run": SCALE_RUN,
    },
    {
        "id": "hypergrass-pangenome",
        "kind": "hgs",
        "heap": "1g",
        "archive": "archive.hgs",
        "memory": "2g",
        "scale_run": PANGENOME_RUN,
    },
    {
        "id": "hypergrass-lean",
        "kind": "hgs",
        "heap": "8m",
        "archive": "archive.hgs",
        "memory": "256m",
        "scale_run": SCALE_RUN,
    },
    {
        "id": "gzip-9",
        "kind": "gzip",
        "archive": "archive.gz",
        "memory": "2g",
        "scale_run": SCALE_RUN,
    },
    {
        "id": "xz-9",
        "kind": "xz",
        "archive": "archive.xz",
        "memory": "2g",
        "scale_run": SCALE_RUN,
    },
)


def percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    index = max(0, min(len(ordered) - 1, math.ceil(fraction * len(ordered)) - 1))
    return ordered[index]


def query_plan(fasta: Path) -> tuple[list[dict], int]:
    headers, starts, lengths = parse_fasta(fasta)
    records = list(range(0, len(headers), 2))
    if records[-1] != len(headers) - 1:
        records.append(len(headers) - 1)
    queries = []
    for record in records:
        available = lengths[record]
        length = min(QUERY_LENGTH, available)
        offset = min(max(0, int(available * 0.37)), available - length)
        queries.append(
            {
                "record": record,
                "offset": offset,
                "length": length,
                "expected_stream_start": starts[record] + offset,
            }
        )
    return queries, len(headers)


def dominates(left: dict, right: dict) -> bool:
    axes = ("archive_bytes", "query_wall_median_seconds", "restore_wall_seconds",
            "query_peak_rss_bytes")
    no_worse = all(left[axis] <= right[axis] for axis in axes)
    strictly_better = any(left[axis] < right[axis] for axis in axes)
    return no_worse and strictly_better


def breakeven_queries(hgs: dict, baseline: dict) -> int | None:
    """Queries needed for one restore + Q queries to favour HGS on wall time."""
    hgs_rate = hgs["query_wall_median_seconds"]
    base_rate = baseline["query_wall_median_seconds"]
    if base_rate <= hgs_rate:
        return None
    extra_restore = hgs["restore_wall_seconds"] - baseline["restore_wall_seconds"]
    return max(0, math.ceil(extra_restore / (base_rate - hgs_rate)))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-id", default=None)
    parser.add_argument("--timeout", type=float, default=900.0)
    parser.add_argument("--profiles", default=None)
    parser.add_argument("--max-queries", type=int, default=None)
    args = parser.parse_args(argv)

    manifest = json.loads((PANGENOME / "manifest.json").read_text(encoding="utf-8"))
    fasta = PANGENOME / manifest["fasta"]["path"]
    if not fasta.exists():
        raise SystemExit(f"missing pangenome FASTA: {fasta}")
    digest = sha256_file(fasta)
    if digest != manifest["fasta"]["sha256"]:
        raise SystemExit("pangenome SHA-256 does not match manifest")
    if not SCALE_REPORT.exists():
        raise SystemExit(f"missing scale report: {SCALE_REPORT}")
    scale_cache = {}
    scale_report_paths = []

    def scale_entry_for(profile: dict) -> dict:
        report_path = profile["scale_run"] / "report.json"
        if report_path not in scale_cache:
            if not report_path.exists():
                raise SystemExit(f"missing scale report: {report_path}")
            scale_cache[report_path] = {
                entry["profile"]: entry
                for entry in json.loads(report_path.read_text(encoding="utf-8"))["results"]
            }
            scale_report_paths.append(str(report_path.relative_to(PROJECT)))
        return scale_cache[report_path].get(profile["id"])

    image = image_id(PROJECT)
    inspect = subprocess.run(
        ["docker", "image", "inspect", image, "--format", "{{.Os}}/{{.Architecture}}"],
        capture_output=True, text=True, timeout=60,
    )
    if inspect.returncode != 0 or inspect.stdout.strip() != "linux/amd64":
        raise SystemExit("locked image is not available or not linux/amd64")
    if not (PROJECT / "bin/hypergrass/HyperGraSS.class").exists():
        raise SystemExit("compiled classes not found; run build.ps1 or build.sh")

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_id = args.run_id or f"endtoend-{stamp}"
    run_dir = PROJECT / "benchmarks/endtoend/runs" / run_id
    if run_dir.exists():
        raise SystemExit(f"run directory already exists: {run_dir}")
    run_dir.mkdir(parents=True)

    queries, record_count = query_plan(fasta)
    if args.max_queries is not None:
        queries = queries[: args.max_queries]
    profiles = [
        profile for profile in PROFILES
        if args.profiles is None or profile["id"] in args.profiles.split(",")
    ]
    if not profiles:
        raise SystemExit("no profiles selected")
    report = {
        "schema_version": 1,
        "kind": "end-to-end-query",
        "run_id": run_id,
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "image_id": image,
        "input": {
            "path": str(fasta.relative_to(PROJECT)),
            "size_bytes": fasta.stat().st_size,
            "sha256": digest,
            "records": record_count,
        },
        "workload": {
            "query_length": QUERY_LENGTH,
            "query_count": len(queries),
            "queries": queries,
            "restore": "one full decompress per profile",
        },
        "scale_reports": [],
        "profiles": [],
        "status": "running",
    }
    write_json(run_dir / "report.json", report)
    print(f"run {run_id}: {len(queries)} queries x {len(profiles)} profiles")

    for profile in profiles:
        name = profile["id"]
        scale_entry = scale_entry_for(profile)
        if not scale_entry or scale_entry.get("status") != "success":
            raise SystemExit(f"scale report has no successful profile {name}")
        archive = profile["scale_run"] / "work" / name / profile["archive"]
        if not archive.exists():
            raise SystemExit(f"missing archive: {archive}")
        if archive.stat().st_size != scale_entry["artifact_bytes"]:
            raise SystemExit(f"archive size changed: {archive}")
        if sha256_file(archive) != scale_entry["artifact_sha256"]:
            raise SystemExit(f"archive SHA-256 changed: {archive}")

        work = run_dir / "work" / name
        work.mkdir(parents=True, exist_ok=True)
        archive_c = to_container(PROJECT, archive)
        query_results = []
        print(f"  {name} ...", flush=True)
        for query in queries:
            record = query["record"]
            start = query["offset"]
            end = start + query["length"]
            expected = sequence_slice(
                fasta, query["expected_stream_start"], query["length"]
            )
            if profile["kind"] == "hgs":
                command = [
                    "java", f"-Xmx{profile['heap']}", "-XX:ActiveProcessorCount=1",
                    "-cp", "/work/bin", "hypergrass.HyperGraSS", "query",
                    archive_c, str(record), str(start), str(end),
                ]
            else:
                command = [
                    "python3", "-B", "/work/benchmarks/endtoend/extract.py",
                    profile["kind"], archive_c, str(record), str(start), str(end),
                ]
            query_file = work / f"query-{record}.seq"
            operation = run_operation(
                PROJECT, image, run_dir, f"{name}/query-{record}",
                command, work, args.timeout, memory=profile["memory"],
                stdout_artifact=query_file,
            )
            observed = query_file.read_bytes() if query_file.exists() else b""
            observed = observed.rstrip(b"\r\n")
            query_results.append(
                {
                    "record": record,
                    "offset": start,
                    "length": query["length"],
                    "status": operation["status"],
                    "wall_seconds": operation.get("wall_seconds"),
                    "cpu_seconds": operation.get("cpu_seconds"),
                    "peak_rss_bytes": operation.get("peak_rss_bytes"),
                    "exact": observed == expected,
                    "error": operation.get("error"),
                }
            )
        successful = [
            query for query in query_results if query["status"] == "success"
        ]
        exact = sum(1 for query in successful if query["exact"])
        wall_values = [query["wall_seconds"] for query in successful]
        rss_values = [
            query["peak_rss_bytes"] for query in successful
            if query["peak_rss_bytes"] is not None
        ]
        restore = scale_entry["decompress"]
        profile_result = {
            "id": name,
            "kind": profile["kind"],
            "archive_bytes": scale_entry["artifact_bytes"],
            "archive_sha256": scale_entry["artifact_sha256"],
            "dependencies": 0 if profile["kind"] == "hgs" else 1,
            "queries": query_results,
            "query_successful": len(successful),
            "query_exact": exact,
            "query_wall_median_seconds": statistics.median(wall_values) if wall_values else None,
            "query_wall_p95_seconds": percentile(wall_values, 0.95) if wall_values else None,
            "query_wall_max_seconds": max(wall_values) if wall_values else None,
            "query_total_wall_seconds": sum(wall_values) if wall_values else None,
            "query_peak_rss_bytes": max(rss_values) if rss_values else None,
            "restore_wall_seconds": restore["wall_seconds"],
            "restore_cpu_seconds": restore["cpu_seconds"],
            "restore_peak_rss_bytes": restore["peak_rss_bytes"],
        }
        report["profiles"].append(profile_result)
        write_json(run_dir / "report.json", report)
        print(
            f"    {profile_result['query_exact']}/{len(queries)} exact, "
            f"query median {profile_result['query_wall_median_seconds'] * 1000:.0f} ms, "
            f"restore {profile_result['restore_wall_seconds']:.1f} s, "
            f"{profile_result['archive_bytes']:,} B"
        )

    by_id = {profile["id"]: profile for profile in report["profiles"]}
    report["scale_reports"] = scale_report_paths
    hgs_profiles = [
        by_id[name] for name in
        ("hypergrass", "hypergrass-pangenome", "hypergrass-lean")
        if name in by_id
    ]
    baselines = [by_id[name] for name in ("gzip-9", "xz-9") if name in by_id]
    speedups = {}
    breakeven = {}
    for hgs in hgs_profiles:
        for baseline in baselines:
            key = f"{hgs['id']} vs {baseline['id']}"
            speedups[key] = (
                baseline["query_wall_median_seconds"]
                / hgs["query_wall_median_seconds"]
            )
            breakeven[key] = breakeven_queries(hgs, baseline)

    frontier = []
    for candidate in report["profiles"]:
        if not any(
            dominates(other, candidate)
            for other in report["profiles"]
            if other is not candidate
        ):
            frontier.append(candidate["id"])

    report["analysis"] = {
        "query_speedup": speedups,
        "breakeven_queries_one_restore": breakeven,
        "pareto_frontier": frontier,
        "pareto_axes": [
            "archive_bytes", "query_wall_median_seconds",
            "restore_wall_seconds", "query_peak_rss_bytes",
        ],
    }
    report["status"] = "success"
    report["finished_utc"] = datetime.now(timezone.utc).isoformat()
    write_json(run_dir / "report.json", report)

    print("speedups (baseline query median / HGS query median):")
    for key, value in speedups.items():
        print(f"  {key}: {value:.1f}x")
    print("break-even queries for one restore + Q queries:")
    for key, value in breakeven.items():
        print(f"  {key}: {'never' if value is None else value}")
    print(f"Pareto frontier: {', '.join(frontier)}")
    print(f"report: {run_dir / 'report.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
