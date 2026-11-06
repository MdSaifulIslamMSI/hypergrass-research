"""JARVIS3 level control on the 25% subset and the full flat panels.

This is the missing tuned-baseline control for the EIMDL comparison.  It runs
JARVIS3 at an explicit ``-l`` level (default 31) on the exact same inputs as the
published runs, with the same 2 GiB container cap, a fresh-container decode,
streaming byte/SHA-256 verification, and cgroup/rusage metrics.

Usage:
    python -B benchmarks/independent/run_jarvis_level.py \
        --scope subset --level 31 --run-id jarvis-l31-subset-1

    python -B benchmarks/independent/run_jarvis_level.py \
        --scope full --level 31 --run-id jarvis-l31-full-1 --timeout 1800
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil

import run_modes as modes
import run_trio as trio


PROJECT = modes.PROJECT
DATASETS = modes.DATASETS
RUNS = modes.RUNS
JARVIS3 = modes.JARVIS3
pangenome_run = modes.pangenome_run


def opt(value: float | None, spec: str = ".5f") -> str:
    """Format an optional summary metric without crashing aborted runs."""
    return "n/a" if value is None else format(value, spec)


def load_panel(dataset: str) -> tuple[dict, Path, int, str]:
    data_dir = DATASETS / dataset
    manifest = json.loads((data_dir / "manifest.json").read_text(encoding="utf-8"))
    fasta = data_dir / manifest["fasta"]["path"]
    fasta_bytes = fasta.stat().st_size
    fasta_sha = modes.sha256_file(fasta)
    if (
        fasta_bytes != manifest["fasta"]["size_bytes"]
        or fasta_sha != manifest["fasta"]["sha256"]
    ):
        raise SystemExit(f"{dataset}: FASTA does not match its locked manifest")
    return manifest, fasta, fasta_bytes, fasta_sha


def run_jarvis(
    image: str,
    run_dir: Path,
    record: dict,
    level: int,
    timeout: float,
    key: str,
    mode: str,
) -> dict:
    archive = run_dir / "ops" / key / "archive.jc"
    restored = run_dir / "ops" / key / "restored.seq"
    compress = pangenome_run.run_operation(
        PROJECT,
        image,
        run_dir,
        f"{key}/compress",
        [
            JARVIS3,
            "-l",
            str(level),
            "-o",
            pangenome_run.to_container(PROJECT, archive),
            pangenome_run.to_container(PROJECT, PROJECT / record["path"]),
        ],
        run_dir,
        timeout,
        "2g",
    )
    decompress = None
    if compress["status"] == "success" and archive.exists():
        decompress = pangenome_run.run_operation(
            PROJECT,
            image,
            run_dir,
            f"{key}/decompress",
            [
                JARVIS3,
                "-d",
                "-o",
                pangenome_run.to_container(PROJECT, restored),
                pangenome_run.to_container(PROJECT, archive),
            ],
            run_dir,
            timeout,
            "2g",
        )
    return modes.finish(
        modes.make_entry(mode, f"jarvis3-l{level}", record),
        compress,
        decompress,
        archive,
        restored,
        record,
        run_dir,
        key,
    )


def base_report(
    dataset: str,
    args: argparse.Namespace,
    image: str,
    fasta: Path,
    fasta_bytes: int,
    fasta_sha: str,
    scope: str,
) -> dict:
    return {
        "schema_version": 1,
        "benchmark": "jarvis3-level-control",
        "dataset": dataset,
        "run_id": args.run_id,
        "scope": scope,
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "provenance": pangenome_run.provenance(PROJECT),
        "image": {
            "tag": "hypergrass-corpus:jdk27",
            "id": image,
            "memory": "2g",
            "cpus": 2,
        },
        "tools": {
            "jarvis3": {
                "path": JARVIS3,
                "version": "3.7",
                "level": args.level,
            }
        },
        "inputs": {
            "fasta": {
                "path": str(fasta.relative_to(PROJECT)),
                "size_bytes": fasta_bytes,
                "sha256": fasta_sha,
            }
        },
        "protocol": (
            "JARVIS3 is run at the explicit -l level on the identical canonical "
            "flat bytes used by the published comparison. Every archive is "
            "decompressed in a fresh container and checked byte-for-byte and by "
            "SHA-256 against the input."
        ),
        "results": [],
        "summaries": {},
        "status": "running",
        "failures": [],
        "finished_utc": None,
    }


def run_subset(args: argparse.Namespace, image: str) -> int:
    failures = 0
    for dataset in args.datasets:
        dataset_failures = 0
        manifest, fasta, fasta_bytes, fasta_sha = load_panel(dataset)
        run_dir = RUNS / dataset / args.run_id
        if run_dir.exists() and not args.force:
            raise SystemExit(f"{run_dir} exists; use --force to rerun")
        if run_dir.exists():
            shutil.rmtree(run_dir)
        run_dir.mkdir(parents=True)

        records = modes.split_to_flat(fasta, run_dir / "genomes")
        total_bases = sum(record["bases"] for record in records)
        if total_bases != manifest["total_bases"]:
            raise SystemExit(
                f"{dataset}: split bases {total_bases:,} != "
                f"manifest {manifest['total_bases']:,}"
            )
        indices = modes.sample_indices(len(records), args.genomes)
        targets = [records[index] for index in indices]

        report = base_report(
            dataset, args, image, fasta, fasta_bytes, fasta_sha, "subset"
        )
        report["inputs"]["sample"] = {
            "policy": f"{args.genomes} evenly spaced records per panel",
            "indices": indices,
            "targets": [record["name"] for record in targets],
        }
        report["inputs"]["genomes"] = records
        report_path = run_dir / "report.json"
        modes.write_json(report_path, report)

        print(
            f"{dataset}: {len(records)} genomes, {total_bases:,} bases; "
            f"targets {indices}",
            flush=True,
        )

        entries: list[dict] = []

        def record_entry(entry: dict) -> None:
            entries.append(entry)
            report["results"] = entries
            modes.write_json(report_path, report)
            rss = (entry.get("decompress") or {}).get("cgroup_memory_peak_bytes")
            print(
                f"  S {entry['target']:26s} "
                f"{entry['archive_bytes'] if entry['archive_bytes'] else '-':>10} B "
                f"{entry['bits_per_base'] if entry['bits_per_base'] else 0:.5f} bpb "
                f"rss={((rss or 0) / 1048576):.0f}MiB {entry['status']}",
                flush=True,
            )

        for index, record in enumerate(targets, 1):
            print(f"  S {index:02d}/{len(targets)} {record['name']}", flush=True)
            entry = run_jarvis(
                image,
                run_dir,
                record,
                args.level,
                args.timeout,
                f"mode-s/jarvis3-l{args.level}/{Path(record['name']).stem}",
                "S",
            )
            record_entry(entry)
            if entry["status"] != "success":
                failures += 1
                dataset_failures += 1

        key = f"S/jarvis3-l{args.level}"
        summary = modes.summarize(entries)
        report["summaries"] = {key: summary}
        report["status"] = "success" if not dataset_failures else "failed"
        report["failures"] = [
            entry["target"] for entry in entries if entry["status"] != "success"
        ]
        report["finished_utc"] = datetime.now(timezone.utc).isoformat()
        modes.write_json(report_path, report)

        cg_c = summary["compress_cgroup_peak"] or {}
        cg_d = summary["decompress_cgroup_peak"] or {}
        print(
            f"{dataset} {key}: {summary['verified_targets']}/{summary['targets']} "
            f"verified, {summary['total_archive_bytes']:,} B total, "
            f"{opt(summary['total_bpb'])} bpb total, "
            f"wall C/D {summary['compress_wall_seconds_total']:.1f}/"
            f"{summary['decompress_wall_seconds_total']:.1f} s, "
            f"CPU C/D {summary['compress_cpu_seconds_total']:.1f}/"
            f"{summary['decompress_cpu_seconds_total']:.1f} s, "
            f"speed C/D {opt(summary['compress_mib_per_second'], '.2f')}/"
            f"{opt(summary['decompress_mib_per_second'], '.2f')} MiB/s, "
            f"cgroup C/D {cg_c.get('median_bytes', 0) / 1048576:.0f}-"
            f"{cg_c.get('max_bytes', 0) / 1048576:.0f}/"
            f"{cg_d.get('median_bytes', 0) / 1048576:.0f}-"
            f"{cg_d.get('max_bytes', 0) / 1048576:.0f} MiB, "
            f"oom={summary['cgroup_oom_kills']}",
            flush=True,
        )
        print(f"report: {report_path}", flush=True)

    return 2 if failures else 0


def run_full(args: argparse.Namespace, image: str) -> int:
    failures = 0
    for dataset in args.datasets:
        manifest, fasta, fasta_bytes, fasta_sha = load_panel(dataset)
        run_dir = RUNS / dataset / args.run_id
        if run_dir.exists() and not args.force:
            raise SystemExit(f"{run_dir} exists; use --force to rerun")
        if run_dir.exists():
            shutil.rmtree(run_dir)
        run_dir.mkdir(parents=True)

        flat = trio.flat_sequence(fasta, run_dir / "input.seq")
        if flat["bases"] != manifest["total_bases"]:
            raise SystemExit(
                f"{dataset}: flat bases {flat['bases']:,} != "
                f"manifest {manifest['total_bases']:,}"
            )
        record = {
            "name": "full-panel.seq",
            "path": flat["path"],
            "size_bytes": flat["size_bytes"],
            "bases": flat["bases"],
            "sha256": flat["sha256"],
        }

        report = base_report(
            dataset, args, image, fasta, fasta_bytes, fasta_sha, "full"
        )
        report["inputs"]["flat"] = flat
        report_path = run_dir / "report.json"
        modes.write_json(report_path, report)

        print(
            f"{dataset}: full flat panel {flat['bases']:,} bases "
            f"({flat['non_acgt_bases_masked']} masked)",
            flush=True,
        )
        entry = run_jarvis(
            image,
            run_dir,
            record,
            args.level,
            args.timeout,
            f"mode-c/jarvis3-l{args.level}/full-panel",
            "C",
        )
        if entry["status"] != "success":
            failures += 1

        key = f"C/jarvis3-l{args.level}"
        summary = modes.summarize([entry])
        report["results"] = [entry]
        report["summaries"] = {key: summary}
        report["status"] = "success" if entry["status"] == "success" else "failed"
        report["failures"] = [] if entry["status"] == "success" else ["full-panel"]
        report["finished_utc"] = datetime.now(timezone.utc).isoformat()
        modes.write_json(report_path, report)

        cg_c = summary["compress_cgroup_peak"] or {}
        cg_d = summary["decompress_cgroup_peak"] or {}
        print(
            f"{dataset} {key}: {summary['verified_targets']}/{summary['targets']} "
            f"verified, {summary['total_archive_bytes']:,} B total, "
            f"{opt(summary['total_bpb'])} bpb total, "
            f"wall C/D {summary['compress_wall_seconds_total']:.1f}/"
            f"{summary['decompress_wall_seconds_total']:.1f} s, "
            f"CPU C/D {summary['compress_cpu_seconds_total']:.1f}/"
            f"{summary['decompress_cpu_seconds_total']:.1f} s, "
            f"speed C/D {opt(summary['compress_mib_per_second'], '.2f')}/"
            f"{opt(summary['decompress_mib_per_second'], '.2f')} MiB/s, "
            f"cgroup C/D {cg_c.get('median_bytes', 0) / 1048576:.0f}-"
            f"{cg_c.get('max_bytes', 0) / 1048576:.0f}/"
            f"{cg_d.get('median_bytes', 0) / 1048576:.0f}-"
            f"{cg_d.get('max_bytes', 0) / 1048576:.0f} MiB, "
            f"oom={summary['cgroup_oom_kills']}",
            flush=True,
        )
        print(f"report: {report_path}", flush=True)

    return 2 if failures else 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--datasets", nargs="+", default=["salmonella", "klebsiella"])
    parser.add_argument("--scope", choices=("subset", "full"), required=True)
    parser.add_argument("--level", type=int, default=31)
    parser.add_argument("--genomes", type=int, default=6, help="targets per panel")
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--timeout", type=float, default=1800.0)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)

    image = pangenome_run.image_id(PROJECT)
    if args.scope == "subset":
        return run_subset(args, image)
    return run_full(args, image)


if __name__ == "__main__":
    raise SystemExit(main())
