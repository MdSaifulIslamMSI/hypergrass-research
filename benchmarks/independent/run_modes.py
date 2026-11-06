"""Mode S / Mode R subset comparison with JARVIS3 and original GraSS.

The full campaign (96 EIMDL runs) is cut to 25%: six evenly spaced genomes per
locked panel.  On that subset:

* **Mode S** is compared like-for-like against JARVIS3 ``-l 7`` and original
  GraSS ``-e2``.  All three receive the identical per-genome flat sequence
  (headers removed, uppercased, non-ACGT masked to ``A``), because JARVIS3
  accepts raw DNA only and GraSS's raw path is headerless.
* **Mode R** is EIMDL-only leave-one-out against the other 24/22 genomes.
  JARVIS3 and GraSS have no reference mode, so their Mode S columns are
  context only, not a like-for-like ranking against Mode R.

Every operation runs in a fresh locked container; every archive is decompressed
in a different fresh container and checked byte-for-byte by the standalone
streaming verifier.  Container cgroup peak and rusage RSS are recorded for
every phase.

Usage:
    python -B benchmarks/independent/run_modes.py --run-id modes-subset-1
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import shutil
import statistics


PROJECT = Path(__file__).resolve().parent.parent.parent
DATASETS = Path(__file__).resolve().parent / "data"
RUNS = Path(__file__).resolve().parent / "runs"
BSC = PROJECT / "benchmarks/specialized/vendor/libbsc/build/bsc"
GRASS_CLASSES = "/work/benchmarks/specialized/vendor/GraSS/build"
JAVA = "/opt/java/openjdk/bin/java"
JARVIS3 = "/usr/local/bin/JARVIS3"

GRASS_COMMIT = "6d93700f5a2de8841f7cf2c20064413c80df2c4a"
LIBBSC_COMMIT = "1dc676540cb40a41562af6f2db04bd1adce97f90"
GRASS_HEAP = "1500m"

MODE_S_HEAP = "1200m"
MODE_S_MEMORY = "1800m"
MODE_R_HEAP = "1600m"
MODE_R_MEMORY = "2200m"
READ_CHUNK = 8 << 20

# Same canonical normalization as run_trio.py / specialized.run.sequence_only.
NORMALIZE = bytes((value if value in b"ACGT" else ord("A")) for value in range(256))
ACGT = (b"A", b"C", b"G", b"T")


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


pangenome_run = load_module("pangenome_run", PROJECT / "benchmarks/pangenome/run.py")
specialized_run = load_module(
    "specialized_run", PROJECT / "benchmarks/specialized/run.py"
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(READ_CHUNK), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, value: dict) -> None:
    pangenome_run.write_json(path, value)


def eimdl_cmd(heap: str) -> list[str]:
    return [
        JAVA,
        f"-Xmx{heap}",
        "-XX:ActiveProcessorCount=1",
        "-cp",
        "/work/bin",
        "hypergrass.EimdlCli",
    ]


def split_to_flat(fasta: Path, out_dir: Path) -> list[dict]:
    """Write every record as its canonical flat sequence, streaming.

    Byte-for-byte equivalent to ``run_trio.flat_sequence`` applied per record.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    records: list[dict] = []
    out = None
    digest = None
    header = ""
    name = ""
    path = None
    bases = 0
    masked = 0
    source_bytes = 0

    def finalize() -> None:
        assert out is not None and path is not None and digest is not None
        out.close()
        records.append(
            {
                "name": name,
                "path": str(path.relative_to(PROJECT)),
                "size_bytes": path.stat().st_size,
                "bases": bases,
                "non_acgt_masked": masked,
                "sha256": digest.hexdigest(),
                "header": header,
                "source_record_bytes": source_bytes,
            }
        )

    with fasta.open("rb") as stream:
        for line in stream:
            if line.startswith(b">"):
                if out is not None:
                    finalize()
                header = line[1:].rstrip(b"\r\n").decode("ascii", "replace")
                accession = header.split()[0] if header.split() else f"record{len(records)}"
                safe = re.sub(r"[^A-Za-z0-9_.-]", "_", accession)
                name = f"{len(records):02d}_{safe}.seq"
                path = out_dir / name
                out = path.open("wb")
                digest = hashlib.sha256()
                bases = 0
                masked = 0
                source_bytes = len(line)
                continue
            if out is None:
                raise SystemExit("FASTA has sequence data before its first header")
            stripped = line.strip()
            if not stripped:
                continue
            upper = stripped.upper()
            sequence = upper.translate(NORMALIZE)
            out.write(sequence)
            digest.update(sequence)
            bases += len(sequence)
            source_bytes += len(line)
            masked += len(sequence) - sum(upper.count(base) for base in ACGT)

    if out is not None:
        finalize()
    if not records:
        raise SystemExit("FASTA contained no records")
    if len({record["name"] for record in records}) != len(records):
        raise SystemExit("duplicate genome file names after splitting")
    return records


def sample_indices(count: int, wanted: int) -> list[int]:
    """Evenly spaced indices, always including the first and last."""
    if wanted >= count:
        return list(range(count))
    if wanted <= 1:
        return [0]
    indices = sorted({round(i * (count - 1) / (wanted - 1)) for i in range(wanted)})
    if len(indices) != wanted:
        raise SystemExit(
            f"could not select {wanted} distinct indices from {count} records"
        )
    return indices


def make_entry(mode: str, tool: str, record: dict) -> dict:
    return {
        "mode": mode,
        "tool": tool,
        "target": record["name"],
        "input": {
            "path": record["path"],
            "size_bytes": record["size_bytes"],
            "bases": record["bases"],
            "sha256": record["sha256"],
        },
        "status": "error",
        "archive_bytes": None,
        "bits_per_base": None,
        "compression_ratio": None,
        "roundtrip_ok": False,
        "verified_sha256": False,
        "compress": None,
        "decompress": None,
        "verification": None,
        "error": None,
    }


def finish(
    entry: dict,
    compress: dict,
    decompress: dict | None,
    archive: Path,
    restored: Path,
    record: dict,
    run_dir: Path,
    key: str,
) -> dict:
    entry["compress"] = compress
    entry["decompress"] = decompress
    if compress["status"] != "success" or not archive.exists():
        entry["error"] = compress.get("error") or "compression failed"
        return entry
    entry["archive_bytes"] = archive.stat().st_size
    entry["bits_per_base"] = entry["archive_bytes"] * 8 / record["bases"]
    entry["compression_ratio"] = record["size_bytes"] / entry["archive_bytes"]

    verification = None
    if decompress is not None and decompress["status"] == "success" and restored.exists():
        verification = specialized_run.independent_verification(
            PROJECT / record["path"],
            restored,
            archive,
            record["sha256"],
            run_dir / "ops" / key / "verification.json",
        )
    entry["verification"] = verification
    entry["roundtrip_ok"] = bool(verification and verification.get("verdict") == "pass")
    entry["verified_sha256"] = bool(
        verification
        and verification.get("verdict") == "pass"
        and verification.get("original_sha256") == record["sha256"]
    )
    if entry["verified_sha256"]:
        entry["status"] = "success"
    else:
        entry["error"] = "round-trip verification failed"
    return entry


def run_eimdl(
    image: str,
    run_dir: Path,
    mode: str,
    record: dict,
    panel: Path | None,
    timeout: float,
) -> dict:
    stem = Path(record["name"]).stem
    key = f"mode-{mode.lower()}/eimdl/{stem}"
    heap, memory = (MODE_S_HEAP, MODE_S_MEMORY) if mode == "S" else (MODE_R_HEAP, MODE_R_MEMORY)
    archive = run_dir / "ops" / key / "archive.hgs2"
    restored = run_dir / "ops" / key / "restored.seq"

    compress_args = [
        "compress",
        "--target",
        pangenome_run.to_container(PROJECT, PROJECT / record["path"]),
        "--out",
        pangenome_run.to_container(PROJECT, archive),
    ]
    if mode == "S":
        compress_args += ["--single"]
    else:
        assert panel is not None
        compress_args += [
            "--panel",
            pangenome_run.to_container(PROJECT, panel),
            "--exclude",
            record["name"],
        ]
    compress_args += ["--threads", "1"]

    compress = pangenome_run.run_operation(
        PROJECT,
        image,
        run_dir,
        f"{key}/compress",
        eimdl_cmd(heap) + compress_args,
        run_dir,
        timeout,
        memory,
    )
    decompress = None
    if compress["status"] == "success" and archive.exists():
        decompress_args = [
            "decompress",
            "--archive",
            pangenome_run.to_container(PROJECT, archive),
            "--out",
            pangenome_run.to_container(PROJECT, restored),
        ]
        if mode == "R":
            assert panel is not None
            decompress_args += [
                "--panel",
                pangenome_run.to_container(PROJECT, panel),
            ]
        decompress = pangenome_run.run_operation(
            PROJECT,
            image,
            run_dir,
            f"{key}/decompress",
            eimdl_cmd(heap) + decompress_args,
            run_dir,
            timeout,
            memory,
        )
    return finish(
        make_entry(mode, "eimdl", record),
        compress,
        decompress,
        archive,
        restored,
        record,
        run_dir,
        key,
    )


def run_jarvis3(image: str, run_dir: Path, record: dict, timeout: float) -> dict:
    stem = Path(record["name"]).stem
    key = f"mode-s/jarvis3/{stem}"
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
            "7",
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
    return finish(
        make_entry("S", "jarvis3", record),
        compress,
        decompress,
        archive,
        restored,
        record,
        run_dir,
        key,
    )


def run_grass(image: str, run_dir: Path, record: dict, timeout: float) -> dict:
    stem = Path(record["name"]).stem
    key = f"mode-s/grass/{stem}"
    grass_op = run_dir / "ops" / key / "compress"
    grass_op.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(BSC, grass_op / "bsc")
    archive = grass_op / "FinalBsc.bsc"
    grass_dec = run_dir / "ops" / key / "decompress"
    restored = grass_dec / "Final.SSG"
    invocation = (
        f"{JAVA} -Xmx{GRASS_HEAP} -XX:ActiveProcessorCount=1 "
        f"-cp {GRASS_CLASSES} SSG.SSG"
    )
    compress = pangenome_run.run_operation(
        PROJECT,
        image,
        run_dir,
        f"{key}/compress",
        [
            "/bin/sh",
            "-c",
            "chmod +x bsc && exec "
            + invocation
            + " compress "
            + pangenome_run.to_container(PROJECT, PROJECT / record["path"]),
        ],
        grass_op,
        timeout,
        "2g",
    )
    decompress = None
    if compress["status"] == "success" and archive.exists():
        grass_dec.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(archive, grass_dec / "FinalBsc.bsc")
        shutil.copyfile(BSC, grass_dec / "bsc")
        decompress = pangenome_run.run_operation(
            PROJECT,
            image,
            run_dir,
            f"{key}/decompress",
            [
                "/bin/sh",
                "-c",
                "chmod +x bsc && exec " + invocation + " decompress",
            ],
            grass_dec,
            timeout,
            "2g",
        )
    return finish(
        make_entry("S", "grass", record),
        compress,
        decompress,
        archive,
        restored,
        record,
        run_dir,
        key,
    )


def rss_stats(entries: list[dict], phase: str, field: str) -> dict | None:
    values = [
        (entry.get(phase) or {}).get(field) for entry in entries
    ]
    values = [value for value in values if value]
    if not values:
        return None
    return {
        "median_bytes": statistics.median(values),
        "max_bytes": max(values),
    }


def phase_stats(entries: list[dict], phase: str, field: str) -> dict | None:
    values = [(entry.get(phase) or {}).get(field) for entry in entries]
    values = [value for value in values if value is not None]
    if not values:
        return None
    return {
        "median": statistics.median(values),
        "max": max(values),
    }


def summarize(entries: list[dict]) -> dict:
    ok = [entry for entry in entries if entry["status"] == "success"]
    total_bases = sum(entry["input"]["bases"] for entry in entries)
    total_input = sum(entry["input"]["size_bytes"] for entry in entries)
    total_archive = sum(entry["archive_bytes"] or 0 for entry in ok)
    bpbs = [entry["bits_per_base"] for entry in ok]
    compress_wall = sum(
        (entry.get("compress") or {}).get("wall_seconds") or 0 for entry in ok
    )
    decompress_wall = sum(
        (entry.get("decompress") or {}).get("wall_seconds") or 0 for entry in ok
    )
    compress_cpu = sum(
        (entry.get("compress") or {}).get("cpu_seconds") or 0 for entry in ok
    )
    decompress_cpu = sum(
        (entry.get("decompress") or {}).get("cpu_seconds") or 0 for entry in ok
    )
    return {
        "targets": len(entries),
        "verified_targets": len(ok),
        "all_verified_sha256": len(ok) == len(entries),
        "total_bases": total_bases,
        "total_input_bytes": total_input,
        "total_archive_bytes": total_archive,
        "total_bpb": (total_archive * 8 / total_bases) if total_bases and ok else None,
        "median_bpb": statistics.median(bpbs) if bpbs else None,
        "min_bpb": min(bpbs) if bpbs else None,
        "max_bpb": max(bpbs) if bpbs else None,
        "compress_wall_seconds_total": compress_wall,
        "decompress_wall_seconds_total": decompress_wall,
        "compress_cpu_seconds_total": compress_cpu,
        "decompress_cpu_seconds_total": decompress_cpu,
        "compress_wall_seconds": phase_stats(ok, "compress", "wall_seconds"),
        "decompress_wall_seconds": phase_stats(ok, "decompress", "wall_seconds"),
        "compress_cpu_seconds": phase_stats(ok, "compress", "cpu_seconds"),
        "decompress_cpu_seconds": phase_stats(ok, "decompress", "cpu_seconds"),
        "compress_mib_per_second": (
            (total_input / 1048576 / compress_wall) if compress_wall else None
        ),
        "decompress_mib_per_second": (
            (total_input / 1048576 / decompress_wall) if decompress_wall else None
        ),
        "compress_rusage_rss": rss_stats(ok, "compress", "peak_rss_bytes"),
        "decompress_rusage_rss": rss_stats(ok, "decompress", "peak_rss_bytes"),
        "compress_cgroup_peak": rss_stats(ok, "compress", "cgroup_memory_peak_bytes"),
        "decompress_cgroup_peak": rss_stats(ok, "decompress", "cgroup_memory_peak_bytes"),
        "cgroup_oom_kills": sum(
            (entry.get(phase) or {}).get("cgroup_oom_kills") or 0
            for entry in ok
            for phase in ("compress", "decompress")
        ),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--datasets", nargs="+", default=["salmonella", "klebsiella"])
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--genomes", type=int, default=6, help="targets per panel")
    parser.add_argument("--timeout", type=float, default=900.0)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)

    image = pangenome_run.image_id(PROJECT)
    failures = 0
    for dataset in args.datasets:
        data_dir = DATASETS / dataset
        manifest = json.loads((data_dir / "manifest.json").read_text(encoding="utf-8"))
        fasta = data_dir / manifest["fasta"]["path"]
        fasta_bytes = fasta.stat().st_size
        fasta_sha = sha256_file(fasta)
        if (
            fasta_bytes != manifest["fasta"]["size_bytes"]
            or fasta_sha != manifest["fasta"]["sha256"]
        ):
            raise SystemExit(f"{dataset}: FASTA does not match its locked manifest")

        run_dir = RUNS / dataset / args.run_id
        report_path = run_dir / "report.json"
        if run_dir.exists() and not args.force:
            raise SystemExit(f"{run_dir} exists; use --force to rerun")
        if run_dir.exists():
            shutil.rmtree(run_dir)
        run_dir.mkdir(parents=True)

        records = split_to_flat(fasta, run_dir / "genomes")
        total_bases = sum(record["bases"] for record in records)
        if total_bases != manifest["total_bases"]:
            raise SystemExit(
                f"{dataset}: split bases {total_bases:,} != "
                f"manifest {manifest['total_bases']:,}"
            )
        indices = sample_indices(len(records), args.genomes)
        targets = [records[index] for index in indices]

        report = {
            "schema_version": 1,
            "benchmark": "independent-pangenome-mode-subset-rss",
            "dataset": dataset,
            "run_id": args.run_id,
            "started_utc": datetime.now(timezone.utc).isoformat(),
            "provenance": pangenome_run.provenance(PROJECT),
            "image": {
                "tag": "hypergrass-corpus:jdk27",
                "id": image,
                "mode_s": {"heap": MODE_S_HEAP, "memory": MODE_S_MEMORY},
                "mode_r": {"heap": MODE_R_HEAP, "memory": MODE_R_MEMORY},
                "baselines": {"memory": "2g"},
                "cpus": 2,
            },
            "tools": {
                "eimdl": {"entry": "hypergrass.EimdlCli", "threads": 1},
                "jarvis3": {"path": JARVIS3, "version": "3.7", "level": 7},
                "grass": {
                    "repo": "https://github.com/AnshuOishik/GraSS",
                    "commit": GRASS_COMMIT,
                    "heap": GRASS_HEAP,
                    "mode": "default-e2",
                },
                "libbsc": {
                    "repo": "https://github.com/IlyaGrebnov/libbsc",
                    "commit": LIBBSC_COMMIT,
                    "bsc_sha256": sha256_file(BSC),
                },
            },
            "inputs": {
                "fasta": {
                    "path": str(fasta.relative_to(PROJECT)),
                    "size_bytes": fasta_bytes,
                    "sha256": fasta_sha,
                    "records": len(records),
                    "bases": total_bases,
                },
                "sample": {
                    "policy": f"{args.genomes} evenly spaced records per panel",
                    "indices": indices,
                    "targets": [record["name"] for record in targets],
                },
                "genomes": records,
            },
            "protocol": (
                "The full 96-run mode campaign is cut to 25%: six evenly spaced "
                "genomes per panel. Every record is written as its canonical flat "
                "sequence (headers removed, uppercased, non-ACGT masked to A). Mode S "
                "compresses each target alone with EIMDL --single, JARVIS3 -l 7, and "
                "original GraSS -e2 on the identical bytes. Mode R compresses each "
                "target with EIMDL --panel against the other N-1 genomes and "
                "--exclude <target>. JARVIS3 and GraSS expose no reference mode, so "
                "their columns are not a like-for-like baseline for Mode R. Every "
                "archive is decompressed in a fresh container and checked against the "
                "target's exact bytes and SHA-256."
            ),
            "results": [],
            "summaries": {},
            "comparison": [],
            "status": "running",
        }
        write_json(report_path, report)
        print(
            f"{dataset}: {len(records)} genomes, {total_bases:,} bases; "
            f"targets {indices}",
            flush=True,
        )

        entries: list[dict] = []

        def record_entry(entry: dict) -> None:
            entries.append(entry)
            report["results"] = entries
            write_json(report_path, report)
            rss = (entry.get("decompress") or {}).get("cgroup_memory_peak_bytes")
            print(
                f"  {entry['mode']} {entry['tool']:8s} {entry['target']:26s} "
                f"{entry['archive_bytes'] if entry['archive_bytes'] else '-':>10} B "
                f"{entry['bits_per_base'] if entry['bits_per_base'] else 0:.5f} bpb "
                f"rss={((rss or 0) / 1048576):.0f}MiB {entry['status']}",
                flush=True,
            )

        # --- Mode S: EIMDL vs JARVIS3 vs GraSS on identical flat targets ----
        for index, record in enumerate(targets, 1):
            print(f"  S {index:02d}/{len(targets)} {record['name']}", flush=True)
            record_entry(run_eimdl(image, run_dir, "S", record, None, args.timeout))
            record_entry(run_jarvis3(image, run_dir, record, args.timeout))
            record_entry(run_grass(image, run_dir, record, args.timeout))

        # --- Mode R: EIMDL leave-one-out against the full panel -------------
        panel = run_dir / "genomes"
        for index, record in enumerate(targets, 1):
            print(f"  R {index:02d}/{len(targets)} {record['name']}", flush=True)
            record_entry(run_eimdl(image, run_dir, "R", record, panel, args.timeout))

        summaries = {}
        for mode, tools in (("S", ("eimdl", "jarvis3", "grass")), ("R", ("eimdl",))):
            for tool in tools:
                subset = [
                    entry
                    for entry in entries
                    if entry["mode"] == mode and entry["tool"] == tool
                ]
                summaries[f"{mode}/{tool}"] = summarize(subset)

        by_target: dict[str, dict] = {}
        for entry in entries:
            by_target.setdefault(entry["target"], {})[
                f"{entry['mode'].lower()}_{entry['tool']}"
            ] = entry["archive_bytes"]
        comparison = []
        for target in [record["name"] for record in targets]:
            row = {"target": target, **by_target.get(target, {})}
            if row.get("s_eimdl") and row.get("s_jarvis3"):
                row["s_eimdl_vs_jarvis3_pct"] = (
                    100 * (row["s_eimdl"] - row["s_jarvis3"]) / row["s_jarvis3"]
                )
            if row.get("s_eimdl") and row.get("s_grass"):
                row["s_eimdl_vs_grass_pct"] = (
                    100 * (row["s_eimdl"] - row["s_grass"]) / row["s_grass"]
                )
            if row.get("r_eimdl") and row.get("s_eimdl"):
                row["r_vs_s_pct"] = (
                    100 * (row["r_eimdl"] - row["s_eimdl"]) / row["s_eimdl"]
                )
            comparison.append(row)

        failed = [entry["target"] for entry in entries if entry["status"] != "success"]
        failures += len(failed)
        report = {
            **report,
            "results": entries,
            "summaries": summaries,
            "comparison": comparison,
            "status": "success" if not failed else "failed",
            "failures": failed,
            "finished_utc": datetime.now(timezone.utc).isoformat(),
        }
        write_json(report_path, report)

        for key, summary in summaries.items():
            def mib(value: float | None) -> str:
                return f"{value / 1048576:.0f}" if value else "-"

            cg_c = summary["compress_cgroup_peak"] or {}
            cg_d = summary["decompress_cgroup_peak"] or {}
            print(
                f"{dataset} {key}: {summary['verified_targets']}/{summary['targets']} "
                f"verified, {summary['total_archive_bytes']:,} B total, "
                f"{summary['total_bpb']:.5f} bpb total, "
                f"{summary['median_bpb']:.5f} bpb median, "
                f"wall C/D {summary['compress_wall_seconds_total']:.1f}/"
                f"{summary['decompress_wall_seconds_total']:.1f} s, "
                f"CPU C/D {summary['compress_cpu_seconds_total']:.1f}/"
                f"{summary['decompress_cpu_seconds_total']:.1f} s, "
                f"speed C/D {summary['compress_mib_per_second']:.2f}/"
                f"{summary['decompress_mib_per_second']:.2f} MiB/s, "
                f"RSS compress {mib(cg_c.get('median_bytes'))}/"
                f"{mib(cg_c.get('max_bytes'))} MiB (med/max), "
                f"restore {mib(cg_d.get('median_bytes'))}/"
                f"{mib(cg_d.get('max_bytes'))} MiB, "
                f"oom={summary['cgroup_oom_kills']}",
                flush=True,
            )
        print(f"report: {report_path}", flush=True)

    return 2 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
