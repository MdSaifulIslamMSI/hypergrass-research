"""EIMDL vs JARVIS3 vs original GraSS on the two locked independent pangenomes.

JARVIS3 v3.7 accepts raw DNA only (it rejects FASTA with "currently JARVIS3 only
supports DNA format"), and the original GraSS raw path consumes a headerless
sequence.  So the comparison feeds all three codecs the *identical* canonical
flat sequence derived from the locked multi-FASTA: headers removed, letters
uppercased, non-ACGT replaced with A.  That is the same normalization used by
``benchmarks/specialized/run.py``.

An optional ``--with-fasta-row`` adds an EIMDL-only row over the untouched
multi-FASTA.  It is off by default: EIMDL reconstructs FASTA line/header bytes
through an unbuffered stream, which is extremely slow on a Windows bind mount
and would dominate the run without changing any archive size.

Every compression is measured in a fresh locked container; every archive is
decompressed in a different fresh container; and the restored bytes are checked
against the exact input by the streaming standalone verifier
(``benchmarks/verify/verify.py``).

Usage:
    python -B benchmarks/independent/run_trio.py --run-id trio-1
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import sys


PROJECT = Path(__file__).resolve().parent.parent.parent
DATASETS = Path(__file__).resolve().parent / "data"
RUNS = Path(__file__).resolve().parent / "runs"
BSC = PROJECT / "benchmarks/specialized/vendor/libbsc/build/bsc"
GRASS_CLASSES = "/work/benchmarks/specialized/vendor/GraSS/build"
JAVA = "/opt/java/openjdk/bin/java"
JARVIS3 = "/usr/local/bin/JARVIS3"

# Pinned by benchmarks/specialized/run.py and the local vendor checkout.
GRASS_COMMIT = "6d93700f5a2de8841f7cf2c20064413c80df2c4a"
LIBBSC_COMMIT = "1dc676540cb40a41562af6f2db04bd1adce97f90"

# EIMDL Mode C needs the heap proven by benchmarks/fusion (2200 MiB on a
# comparable 122 MB panel).  The container cap leaves room for JVM overhead.
EIMDL_HEAP = "2200m"
EIMDL_MEMORY = "2800m"

READ_CHUNK = 8 << 20


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


def flat_sequence(fasta: Path, destination: Path) -> dict:
    """Stream the canonical flat sequence so host RAM stays flat.

    Byte-for-byte equivalent to ``specialized.run.sequence_only`` but written in
    chunks instead of joining the whole 122 MB file in memory.
    """
    table = bytes(
        (value if value in b"ACGT" else ord("A")) for value in range(256)
    )
    counts = {base: 0 for base in b"ACGT"}
    total = 0
    destination.parent.mkdir(parents=True, exist_ok=True)
    with fasta.open("rb") as stream, destination.open("wb") as out:
        for line in stream:
            if line.startswith(b">"):
                continue
            stripped = line.strip()
            if not stripped:
                continue
            upper = stripped.upper()
            for base in counts:
                counts[base] += upper.count(bytes([base]))
            masked = upper.translate(table)
            out.write(masked)
            total += len(masked)
    digest = sha256_file(destination)
    return {
        "path": str(destination.relative_to(PROJECT)),
        "size_bytes": total,
        "bases": total,
        "sha256": digest,
        "non_acgt_bases_masked": total - sum(counts.values()),
        "base_counts": {chr(key): value for key, value in counts.items()},
    }


def last_json_line(path: Path) -> dict | None:
    if not path.exists():
        return None
    for line in reversed(path.read_text(encoding="utf-8", errors="replace").splitlines()):
        stripped = line.strip()
        if not stripped.startswith("{"):
            continue
        try:
            return json.loads(stripped)
        except json.JSONDecodeError:
            continue
    return None


def verify(
    original: Path, restored: Path, archive: Path, expected: str, report: Path
) -> dict:
    return specialized_run.independent_verification(
        original, restored, archive, expected, report
    )


def row(
    label: str,
    kind: str,
    input_path: Path,
    input_sha: str,
    input_bytes: int,
    compress: dict,
    decompress: dict | None,
    archive: Path,
    verification: dict | None,
    tool: dict,
) -> dict:
    archive_bytes = archive.stat().st_size if archive.exists() else None
    return {
        "profile": label,
        "kind": kind,
        "tool": tool,
        "input": {
            "path": str(input_path.relative_to(PROJECT)),
            "bytes": input_bytes,
            "sha256": input_sha,
        },
        "status": "success"
        if compress["status"] == "success" and verification
        and verification.get("verdict") == "pass"
        else "failure",
        "archive_bytes": archive_bytes,
        "bits_per_base": (archive_bytes * 8 / input_bytes) if archive_bytes else None,
        "compression_ratio": (input_bytes / archive_bytes) if archive_bytes else None,
        "roundtrip_ok": bool(
            verification and verification.get("verdict") == "pass"
        ),
        "verified_sha256": bool(
            verification
            and verification.get("verdict") == "pass"
            and verification.get("original_sha256") == input_sha
        ),
        "compress": compress,
        "decompress": decompress,
        "verification": verification,
    }


def eimdl_fasta_row(
    image: str,
    run_dir: Path,
    timeout: float,
    eimdl_cmd: list[str],
    fasta: Path,
    fasta_sha: str,
    fasta_bytes: int,
    bases: int,
) -> dict:
    """Optional EIMDL-only row over the untouched multi-FASTA."""
    archive = run_dir / "ops" / "eimdl-fasta" / "compress" / "archive.hgs2"
    restored = run_dir / "ops" / "eimdl-fasta" / "decompress" / "restored.fasta"
    compress = pangenome_run.run_operation(
        PROJECT,
        image,
        run_dir,
        "eimdl-fasta/compress",
        eimdl_cmd
        + [
            "compress",
            "--target",
            pangenome_run.to_container(PROJECT, fasta),
            "--out",
            pangenome_run.to_container(PROJECT, archive),
            "--collection",
            "--threads",
            "1",
        ],
        run_dir,
        timeout,
        EIMDL_MEMORY,
    )
    decompress = None
    verification = None
    if compress["status"] == "success" and archive.exists():
        decompress = pangenome_run.run_operation(
            PROJECT,
            image,
            run_dir,
            "eimdl-fasta/decompress",
            eimdl_cmd
            + [
                "decompress",
                "--archive",
                pangenome_run.to_container(PROJECT, archive),
                "--out",
                pangenome_run.to_container(PROJECT, restored),
            ],
            run_dir,
            timeout,
            EIMDL_MEMORY,
        )
        if decompress["status"] == "success" and restored.exists():
            verification = verify(
                fasta,
                restored,
                archive,
                fasta_sha,
                run_dir / "ops" / "eimdl-fasta" / "verification.json",
            )
    entry = row(
        "eimdl-fasta",
        "eimdl",
        fasta,
        fasta_sha,
        fasta_bytes,
        compress,
        decompress,
        archive,
        verification,
        {
            "artifact": "HGS2 Mode C, FASTA-preserving",
            "stats": last_json_line(
                run_dir / "ops" / "eimdl-fasta" / "compress" / "stdout.log"
            ),
        },
    )
    entry["bases"] = bases
    entry["bits_per_base"] = (
        entry["archive_bytes"] * 8 / bases if entry["archive_bytes"] else None
    )
    return entry


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--datasets", nargs="+", default=["salmonella", "klebsiella"]
    )
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--timeout", type=float, default=3600.0)
    parser.add_argument("--with-fasta-row", action="store_true")
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
        if fasta_bytes != manifest["fasta"]["size_bytes"] or fasta_sha != manifest["fasta"]["sha256"]:
            raise SystemExit(
                f"{dataset}: FASTA does not match its locked manifest "
                f"({fasta_bytes:,} B, {fasta_sha})"
            )

        run_dir = RUNS / dataset / args.run_id
        report_path = run_dir / "report.json"
        if run_dir.exists() and not args.force:
            raise SystemExit(f"{run_dir} exists; use --force to rerun")
        if run_dir.exists():
            shutil.rmtree(run_dir)
        run_dir.mkdir(parents=True)

        records, bases = specialized_run.parse_input(fasta)
        if bases != manifest["total_bases"]:
            raise SystemExit(
                f"{dataset}: {bases:,} bases != manifest {manifest['total_bases']:,}"
            )
        flat = flat_sequence(fasta, run_dir / "input.seq")
        if flat["bases"] != bases:
            raise SystemExit(f"{dataset}: flat sequence lost bases")

        protocol = (
            "All codecs receive the identical canonical flat sequence (headers "
            "removed, uppercased, non-ACGT masked to A), because JARVIS3 v3.7 "
            "accepts raw DNA only and original GraSS's raw path is headerless. "
            "Every archive is decompressed in a fresh container and compared by "
            "an independent streaming verifier."
        )
        if args.with_fasta_row:
            protocol += (
                " An extra EIMDL-only row archives the untouched multi-FASTA and "
                "is verified against the FASTA hash; it is not part of the "
                "like-for-like ranking."
            )
        report = {
            "schema_version": 1,
            "benchmark": "independent-pangenome-trio",
            "dataset": dataset,
            "run_id": args.run_id,
            "started_utc": datetime.now(timezone.utc).isoformat(),
            "provenance": pangenome_run.provenance(PROJECT),
            "image": {
                "tag": "hypergrass-corpus:jdk27",
                "id": image,
                "memory": EIMDL_MEMORY,
                "cpus": 2,
            },
            "tools": {
                "eimdl": {"entry": "hypergrass.EimdlCli", "heap": EIMDL_HEAP},
                "jarvis3": {
                    "path": JARVIS3,
                    "version": "3.7",
                    "level": 7,
                    "input": "raw DNA only",
                },
                "grass": {
                    "repo": "https://github.com/AnshuOishik/GraSS",
                    "commit": GRASS_COMMIT,
                    "mode": "default-e2",
                    "heap": "1500m",
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
                    "records": records,
                    "bases": bases,
                    "source": "locked manifest, re-fetched by accessions",
                },
                "flat": flat,
            },
            "protocol": protocol,
            "results": [],
            "status": "running",
        }
        write_json(report_path, report)
        print(
            f"{dataset}: {records} records, {bases:,} bases, "
            f"FASTA {fasta_bytes:,} B, flat {flat['size_bytes']:,} B "
            f"({flat['non_acgt_bases_masked']} masked)",
            flush=True,
        )

        flat_path = run_dir / "input.seq"
        entries = []

        eimdl_cmd = [
            JAVA,
            f"-Xmx{EIMDL_HEAP}",
            "-XX:ActiveProcessorCount=1",
            "-cp",
            "/work/bin",
            "hypergrass.EimdlCli",
        ]

        # --- EIMDL, Mode C, identical flat input -----------------------------
        archive = run_dir / "ops" / "eimdl-flat" / "compress" / "archive.hgs2"
        restored = run_dir / "ops" / "eimdl-flat" / "decompress" / "restored.seq"
        compress = pangenome_run.run_operation(
            PROJECT,
            image,
            run_dir,
            "eimdl-flat/compress",
            eimdl_cmd
            + [
                "compress",
                "--target",
                pangenome_run.to_container(PROJECT, flat_path),
                "--out",
                pangenome_run.to_container(PROJECT, archive),
                "--collection",
                "--threads",
                "1",
            ],
            run_dir,
            args.timeout,
            EIMDL_MEMORY,
        )
        decompress = None
        verification = None
        if compress["status"] == "success" and archive.exists():
            decompress = pangenome_run.run_operation(
                PROJECT,
                image,
                run_dir,
                "eimdl-flat/decompress",
                eimdl_cmd
                + [
                    "decompress",
                    "--archive",
                    pangenome_run.to_container(PROJECT, archive),
                    "--out",
                    pangenome_run.to_container(PROJECT, restored),
                ],
                run_dir,
                args.timeout,
                EIMDL_MEMORY,
            )
            if decompress["status"] == "success" and restored.exists():
                verification = verify(
                    flat_path,
                    restored,
                    archive,
                    flat["sha256"],
                    run_dir / "ops" / "eimdl-flat" / "verification.json",
                )
        entries.append(
            row(
                "eimdl-flat",
                "eimdl",
                flat_path,
                flat["sha256"],
                flat["size_bytes"],
                compress,
                decompress,
                archive,
                verification,
                {
                    "artifact": "HGS2 Mode C",
                    "stats": last_json_line(
                        run_dir / "ops" / "eimdl-flat" / "compress" / "stdout.log"
                    ),
                },
            )
        )
        write_json(report_path, {**report, "results": entries})

        # --- JARVIS3, level 7, identical flat input --------------------------
        archive = run_dir / "ops" / "jarvis3/compress" / "archive.jc"
        restored = run_dir / "ops" / "jarvis3/decompress" / "restored.seq"
        compress = pangenome_run.run_operation(
            PROJECT,
            image,
            run_dir,
            "jarvis3/compress",
            [
                JARVIS3,
                "-l",
                "7",
                "-o",
                pangenome_run.to_container(PROJECT, archive),
                pangenome_run.to_container(PROJECT, flat_path),
            ],
            run_dir,
            args.timeout,
            "2g",
        )
        decompress = None
        verification = None
        if compress["status"] == "success" and archive.exists():
            decompress = pangenome_run.run_operation(
                PROJECT,
                image,
                run_dir,
                "jarvis3/decompress",
                [
                    JARVIS3,
                    "-d",
                    "-o",
                    pangenome_run.to_container(PROJECT, restored),
                    pangenome_run.to_container(PROJECT, archive),
                ],
                run_dir,
                args.timeout,
                "2g",
            )
            if decompress["status"] == "success" and restored.exists():
                verification = verify(
                    flat_path,
                    restored,
                    archive,
                    flat["sha256"],
                    run_dir / "ops" / "jarvis3" / "verification.json",
                )
        entries.append(
            row(
                "jarvis3-l7",
                "jarvis3",
                flat_path,
                flat["sha256"],
                flat["size_bytes"],
                compress,
                decompress,
                archive,
                verification,
                {"artifact": "JARVIS3 .jc", "level": 7},
            )
        )
        write_json(report_path, {**report, "results": entries})

        # --- original GraSS, default -e2, identical flat input --------------
        grass_op = run_dir / "ops" / "grass/compress"
        grass_op.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(BSC, grass_op / "bsc")
        archive = grass_op / "FinalBsc.bsc"
        grass_dec = run_dir / "ops" / "grass/decompress"
        restored = grass_dec / "Final.SSG"
        grass_invocation = (
            f"{JAVA} -Xmx1500m -XX:ActiveProcessorCount=1 "
            f"-cp {GRASS_CLASSES} SSG.SSG"
        )
        compress = pangenome_run.run_operation(
            PROJECT,
            image,
            run_dir,
            "grass/compress",
            [
                "/bin/sh",
                "-c",
                "chmod +x bsc && exec "
                + grass_invocation
                + " compress "
                + pangenome_run.to_container(PROJECT, flat_path),
            ],
            grass_op,
            args.timeout,
            "2g",
        )
        decompress = None
        verification = None
        if compress["status"] == "success" and archive.exists():
            grass_dec.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(archive, grass_dec / "FinalBsc.bsc")
            shutil.copyfile(BSC, grass_dec / "bsc")
            decompress = pangenome_run.run_operation(
                PROJECT,
                image,
                run_dir,
                "grass/decompress",
                [
                    "/bin/sh",
                    "-c",
                    "chmod +x bsc && exec "
                    + grass_invocation
                    + " decompress",
                ],
                grass_dec,
                args.timeout,
                "2g",
            )
            if decompress["status"] == "success" and restored.exists():
                verification = verify(
                    flat_path,
                    restored,
                    archive,
                    flat["sha256"],
                    run_dir / "ops" / "grass" / "verification.json",
                )
        entries.append(
            row(
                "grass-e2",
                "grass",
                flat_path,
                flat["sha256"],
                flat["size_bytes"],
                compress,
                decompress,
                archive,
                verification,
                {
                    "artifact": "FinalBsc.bsc",
                    "grass_commit": GRASS_COMMIT,
                    "libbsc_commit": LIBBSC_COMMIT,
                },
            )
        )
        write_json(report_path, {**report, "results": entries})

        if args.with_fasta_row:
            entries.append(
                eimdl_fasta_row(
                    image,
                    run_dir,
                    args.timeout,
                    eimdl_cmd,
                    fasta,
                    fasta_sha,
                    fasta_bytes,
                    bases,
                )
            )
            write_json(report_path, {**report, "results": entries})

        failed = [e["profile"] for e in entries if e["status"] != "success"]
        failures += len(failed)
        report = {
            **report,
            "results": entries,
            "ranking_flat_input": [
                {
                    "rank": i + 1,
                    "profile": e["profile"],
                    "archive_bytes": e["archive_bytes"],
                    "bits_per_base": e["bits_per_base"],
                }
                for i, e in enumerate(
                    sorted(
                        (e for e in entries if e["archive_bytes"]),
                        key=lambda e: e["archive_bytes"],
                    )
                )
            ],
            "status": "success" if not failed else "failed",
            "failures": failed,
            "finished_utc": datetime.now(timezone.utc).isoformat(),
        }
        write_json(report_path, report)

        print(f"{dataset} results:", flush=True)
        for e in entries:
            print(
                f"  {e['profile']:12s} {e['status']:8s} "
                f"{e['archive_bytes'] if e['archive_bytes'] is not None else '-':>12} B  "
                f"{e['bits_per_base']:.4f} bpb  "
                f"wall={(e['compress'].get('wall_seconds') or 0):.1f}s  "
                f"sha256={e['verified_sha256']}",
                flush=True,
            )
        print(f"report: {report_path}", flush=True)

    return 2 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
