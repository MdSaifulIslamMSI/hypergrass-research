"""Run specialized DNA compressors on the >100 MB pangenome.

JARVIS3 and GeCo3 both document FASTA input and use only the A/C/G/T sequence,
so this harness feeds them the exact same multi-FASTA as the HyperGraSS scale
test. Their archives exclude metadata; HyperGraSS and xz archives include it.
That asymmetry slightly favors the specialized tools, which is the conservative
direction for a HyperGraSS comparison.

A small prefix slice can be used as a pilot to estimate full-scale runtime
before committing hours to a 121.8 MB run.

Usage:
    python -B benchmarks/specialized/run.py \
        --run-id specialized-pilot-5m --slice-bases 5000000 --timeout 1800
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys


PROJECT = Path(__file__).resolve().parent.parent.parent
DATA = Path(__file__).resolve().parent / "data"
RUNS = Path(__file__).resolve().parent / "runs"
SOURCE = PROJECT / "benchmarks/pangenome/data/pangenome.fasta"

SPEC = importlib.util.spec_from_file_location(
    "pangenome_run", PROJECT / "benchmarks/pangenome/run.py"
)
pangenome_run = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(pangenome_run)

PROFILES = (
    {"id": "jarvis3-l7", "kind": "jarvis3", "level": "7"},
    {"id": "jarvis3-l26", "kind": "jarvis3", "level": "26"},
    {"id": "geco3-l3", "kind": "geco3", "level": "3"},
    {"id": "geco3-l6", "kind": "geco3", "level": "6"},
    {"id": "geco3-l9", "kind": "geco3", "level": "9"},
    {"id": "xz-9", "kind": "xz", "level": "9"},
    {
        "id": "hypergrass-pangenome",
        "kind": "hypergrass",
        "level": "consensus-anchor",
        "options": [
            "--pangenome", "--single-thread", "--anchor-single-page",
            "--block-mib", "32", "--anchor-mib", "9",
        ],
        "jvm": ["-Dhypergrass.matchMin=104"],
        "heap": "1536m",
        "memory": "2g",
    },
    {
        "id": "hypergrass-bounded-256m",
        "kind": "hypergrass",
        "level": "bounded-256m",
        "options": [
            "--pangenome", "--bounded", "--single-thread", "--anchor-single-page",
            "--block-mib", "2", "--anchor-mib", "9",
        ],
        "jvm": ["-Dhypergrass.matchMin=104"],
        "heap": "192m",
        "memory": "256m",
    },
)

VERIFIER = Path(__file__).resolve().parent.parent / "verify" / "verify.py"


def independent_verification(
    original: Path,
    restored: Path,
    archive: Path,
    expected_sha256: str,
    report_path: Path,
) -> dict:
    """Run the standalone verifier as its own process (never imports the codec)."""
    command = [
        sys.executable, "-B", str(VERIFIER),
        "--original", str(original),
        "--restored", str(restored),
        "--archive", str(archive),
        "--expect-sha256", expected_sha256,
        "--report", str(report_path),
    ]
    run = subprocess.run(
        command, cwd=PROJECT, capture_output=True, text=True, timeout=7200,
    )
    if report_path.exists():
        return json.loads(report_path.read_text(encoding="utf-8"))
    return {
        "verdict": "fail",
        "problems": [f"verifier produced no report (exit {run.returncode})"],
        "stderr": run.stderr[-2000:],
    }


def slice_fasta(source: Path, target_bases: int, destination: Path) -> tuple[int, int, int]:
    """Copy sequence until target_bases have been written.

    The last record may be truncated, which keeps the slice close to the
    requested size even when individual genomes are several Mbp.
    """
    destination.parent.mkdir(parents=True, exist_ok=True)
    records = 0
    bases = 0
    with source.open("rb") as stream, destination.open("wb") as out:
        for line in stream:
            if line.startswith(b">"):
                if bases >= target_bases:
                    break
                out.write(line)
                records += 1
                continue
            sequence = line.strip()
            remaining = target_bases - bases
            if len(sequence) > remaining:
                out.write(sequence[:remaining] + b"\n")
                bases += remaining
                break
            out.write(line)
            bases += len(sequence)
    if bases < target_bases:
        raise SystemExit(f"slice has only {bases:,} bases; target was {target_bases:,}")
    return records, bases, destination.stat().st_size


def sequence_only(fasta: Path, destination: Path) -> dict:
    """Extract one wrapped A/C/G/T sequence, masking rare ambiguity codes."""
    chunks = []
    with fasta.open("rb") as stream:
        for line in stream:
            if line.startswith(b">"):
                continue
            stripped = line.strip()
            if stripped:
                chunks.append(stripped)
    upper = b"".join(chunks).upper()
    counts = {base: upper.count(bytes([base])) for base in b"ACGT"}
    non_acgt = len(upper) - sum(counts.values())
    table = bytes((value if value in b"ACGT" else ord("A")) for value in range(256))
    masked = upper.translate(table)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(masked)
    sequence_digest = hashlib.sha256(masked).hexdigest()
    return {
        "path": str(destination.relative_to(PROJECT)),
        "size_bytes": len(masked),
        "bases": len(masked),
        "sha256": sequence_digest,
        "non_acgt_bases_masked": non_acgt,
        "base_counts": {chr(key): value for key, value in counts.items()},
    }


def sequence_sha256(path: Path) -> str:
    payload = path.read_bytes()
    return hashlib.sha256(payload.translate(None, b"\r\n\t ")).hexdigest()


def parse_input(path: Path) -> tuple[int, int]:
    """Return (records, bases) for a FASTA file without loading it into memory."""
    records = 0
    bases = 0
    with path.open("rb") as stream:
        for line in stream:
            if line.startswith(b">"):
                records += 1
            else:
                bases += len(line.strip())
    return records, bases


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-id", default=None)
    parser.add_argument("--source", default=str(SOURCE))
    parser.add_argument("--slice-bases", type=int, default=0)
    parser.add_argument("--profiles", default=",".join(p["id"] for p in PROFILES))
    parser.add_argument("--timeout", type=float, default=7200.0)
    parser.add_argument("--memory", default="2g")
    args = parser.parse_args(argv)

    source = Path(args.source)
    if not source.exists():
        raise SystemExit(f"missing source FASTA: {source}")

    image = pangenome_run.image_id(PROJECT)
    wanted = set(args.profiles.split(","))
    profiles = [p for p in PROFILES if p["id"] in wanted]
    missing = wanted - {p["id"] for p in profiles}
    if missing:
        raise SystemExit(f"unknown profiles: {', '.join(sorted(missing))}")
    if not profiles:
        raise SystemExit("no profiles selected")

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_id = args.run_id or f"specialized-{stamp}"
    run_dir = RUNS / run_id
    if run_dir.exists():
        raise SystemExit(f"run directory already exists: {run_dir}")
    run_dir.mkdir(parents=True)

    if args.slice_bases:
        fasta_path = run_dir / "input.fasta"
        records, bases, size = slice_fasta(source, args.slice_bases, fasta_path)
    else:
        fasta_path = source
        records, bases = parse_input(source)
        size = source.stat().st_size

    fasta_digest = pangenome_run.sha256_file(fasta_path)
    sequence = sequence_only(fasta_path, run_dir / "input.seq")
    report = {
        "schema_version": 1,
        "kind": "specialized-scale-test",
        "run_id": run_id,
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "image_id": image,
        "source": str(source.relative_to(PROJECT)),
        "sliced": bool(args.slice_bases),
        "input": {
            "fasta": {
                "path": str(fasta_path.relative_to(PROJECT)),
                "size_bytes": size,
                "sha256": fasta_digest,
                "records": records,
                "bases": bases,
            },
            "sequence": sequence,
            "normalization": (
                "headers removed, letters uppercased, non-ACGT replaced with A; "
                "both tools receive this sequence and are verified against it"
            ),
        },
        "profiles": [p["id"] for p in profiles],
        "results": [],
        "status": "running",
    }
    if not args.slice_bases:
        baseline_path = PROJECT / "benchmarks/pangenome/runs/pangenome-scale-2/report.json"
        if baseline_path.exists():
            baseline = json.loads(baseline_path.read_text(encoding="utf-8"))
            for entry in baseline.get("results", []):
                if entry.get("profile") == "xz-9" and entry.get("status") == "success":
                    report["baseline"] = {
                        "profile": "xz-9",
                        "archive_bytes": entry["artifact_bytes"],
                        "source": str(baseline_path.relative_to(PROJECT)),
                    }
    pangenome_run.write_json(run_dir / "report.json", report)
    print(
        f"run {run_id}: {len(profiles)} profiles on {sequence['size_bytes']:,} bytes, "
        f"{sequence['bases']:,} bases, {sequence['non_acgt_bases_masked']} masked"
    )

    failures = 0
    for profile in profiles:
        print(f"  {profile['id']} ...", flush=True)
        entry = {
            "profile": profile["id"],
            "kind": profile["kind"],
            "level": profile["level"],
            "input_bytes": sequence["size_bytes"],
            "input_sha256": sequence["sha256"],
            "input_bases": sequence["bases"],
            "status": "error",
            "roundtrip_ok": False,
            "artifact_bytes": None,
            "bits_per_base": None,
            "wall_seconds": None,
            "cpu_seconds": None,
            "peak_rss_bytes": None,
            "error": None,
        }
        artifact: Path | None = None
        restored: Path | None = None
        memory = profile.get("memory", args.memory)
        try:
            if profile["kind"] == "jarvis3":
                op_dir = run_dir / "ops" / profile["id"]
                artifact = op_dir / "archive.jc"
                restored = op_dir / "restored.seq"
                compress = pangenome_run.run_operation(
                    PROJECT, image, run_dir, f"{profile['id']}/compress",
                    ["JARVIS3", "-l", profile["level"], "-o",
                     pangenome_run.to_container(PROJECT, artifact),
                     pangenome_run.to_container(PROJECT, run_dir / "input.seq")],
                    run_dir, args.timeout, memory,
                )
                entry.update(compress)
                if compress["status"] == "success" and artifact.exists():
                    decompress = pangenome_run.run_operation(
                        PROJECT, image, run_dir, f"{profile['id']}/decompress",
                        ["JARVIS3", "-d", "-o",
                         pangenome_run.to_container(PROJECT, restored),
                         pangenome_run.to_container(PROJECT, artifact)],
                        run_dir, args.timeout, memory,
                    )
                    entry["decompress"] = decompress
                    entry["roundtrip_ok"] = (
                        decompress["status"] == "success"
                        and restored.exists()
                        and sequence_sha256(restored) == sequence["sha256"]
                    )
            elif profile["kind"] == "geco3":
                op_dir = run_dir / "ops" / profile["id"]
                op_dir.mkdir(parents=True, exist_ok=True)
                local_input = op_dir / "input.seq"
                shutil.copyfile(run_dir / "input.seq", local_input)
                artifact = op_dir / "input.seq.co"
                restored = op_dir / "input.seq.de"
                compress = pangenome_run.run_operation(
                    PROJECT, image, run_dir, f"{profile['id']}/compress",
                    ["GeCo3", "-l", profile["level"],
                     pangenome_run.to_container(PROJECT, local_input)],
                    op_dir, args.timeout, args.memory,
                )
                entry.update(compress)
                if compress["status"] == "success" and artifact.exists():
                    decompress = pangenome_run.run_operation(
                        PROJECT, image, run_dir, f"{profile['id']}/decompress",
                        ["GeDe3", pangenome_run.to_container(PROJECT, artifact)],
                        op_dir, args.timeout, memory,
                    )
                    entry["decompress"] = decompress
                    entry["roundtrip_ok"] = (
                        decompress["status"] == "success"
                        and restored.exists()
                        and sequence_sha256(restored) == sequence["sha256"]
                    )
            elif profile["kind"] == "xz":
                op_dir = run_dir / "ops" / profile["id"]
                artifact = op_dir / "archive.xz"
                restored = op_dir / "restored.seq"
                compress = pangenome_run.run_operation(
                    PROJECT, image, run_dir, f"{profile['id']}/compress",
                    ["xz", "-9", "-c",
                     pangenome_run.to_container(PROJECT, run_dir / "input.seq")],
                    run_dir, args.timeout, memory,
                    stdout_artifact=artifact,
                )
                entry.update(compress)
                if compress["status"] == "success" and artifact.exists():
                    decompress = pangenome_run.run_operation(
                        PROJECT, image, run_dir, f"{profile['id']}/decompress",
                        ["xz", "-dc", pangenome_run.to_container(PROJECT, artifact)],
                        run_dir, args.timeout, memory,
                        stdout_artifact=restored,
                    )
                    entry["decompress"] = decompress
                    entry["roundtrip_ok"] = (
                        decompress["status"] == "success"
                        and restored.exists()
                        and sequence_sha256(restored) == sequence["sha256"]
                    )
            elif profile["kind"] == "hypergrass":
                op_dir = run_dir / "ops" / profile["id"]
                artifact = op_dir / "archive.hgs"
                restored = op_dir / "restored.seq"
                java = ["java", f"-Xmx{profile.get('heap', '4g')}",
                        "-XX:ActiveProcessorCount=1", *profile.get("jvm", []),
                        "-cp", "/work/bin", "hypergrass.HyperGraSS"]
                compress = pangenome_run.run_operation(
                    PROJECT, image, run_dir, f"{profile['id']}/compress",
                    java + ["compress", *profile["options"],
                            pangenome_run.to_container(PROJECT, run_dir / "input.seq"),
                            pangenome_run.to_container(PROJECT, artifact)],
                    run_dir, args.timeout, memory,
                )
                entry.update(compress)
                if compress["status"] == "success" and artifact.exists():
                    decompress = pangenome_run.run_operation(
                        PROJECT, image, run_dir, f"{profile['id']}/decompress",
                        java + ["decompress",
                                pangenome_run.to_container(PROJECT, artifact),
                                pangenome_run.to_container(PROJECT, restored)],
                        run_dir, args.timeout, memory,
                    )
                    entry["decompress"] = decompress
                    entry["roundtrip_ok"] = (
                        decompress["status"] == "success"
                        and restored.exists()
                        and sequence_sha256(restored) == sequence["sha256"]
                    )
            else:
                raise ValueError(f"unknown kind: {profile['kind']}")
        except Exception as error:
            entry["error"] = f"{type(error).__name__}: {error}"

        if entry.get("artifact_bytes") is None:
            if profile["kind"] == "jarvis3":
                artifact = run_dir / "ops" / profile["id"] / "archive.jc"
            elif profile["kind"] == "geco3":
                artifact = run_dir / "ops" / profile["id"] / "input.seq.co"
            elif profile["kind"] == "hypergrass":
                artifact = run_dir / "ops" / profile["id"] / "archive.hgs"
            else:
                artifact = run_dir / "ops" / profile["id"] / "archive.xz"
            if artifact.exists():
                entry["artifact_bytes"] = artifact.stat().st_size
        if entry.get("artifact_bytes"):
            entry["bits_per_base"] = entry["artifact_bytes"] * 8 / sequence["bases"]
        if artifact is not None and restored is not None:
            if artifact.is_file() and restored.is_file() and restored.stat().st_size > 0:
                verification = independent_verification(
                    run_dir / "input.seq", restored, artifact, sequence["sha256"],
                    run_dir / "ops" / profile["id"] / "verification.json",
                )
                entry["verification"] = verification
                entry["restored_sha256"] = verification.get("restored_sha256")
                entry["roundtrip_ok"] = bool(entry["roundtrip_ok"]) and (
                    verification.get("verdict") == "pass"
                )
        if entry["status"] != "success" or not entry["roundtrip_ok"]:
            failures += 1
            print(f"    FAILED: status={entry['status']} roundtrip={entry['roundtrip_ok']} "
                  f"error={entry.get('error')}", flush=True)
        else:
            print(
                f"    {entry['artifact_bytes']:,} B, "
                f"{entry['bits_per_base']:.3f} bits/base, "
                f"{entry['wall_seconds']:.1f}s, "
                f"{entry['peak_rss_bytes'] / 1048576:.0f} MiB RSS",
                flush=True,
            )
        report["results"].append(entry)
        pangenome_run.write_json(run_dir / "report.json", report)

    report["status"] = "success" if failures == 0 else "partial"
    report["failures"] = failures
    report["finished_utc"] = datetime.now(timezone.utc).isoformat()
    pangenome_run.write_json(run_dir / "report.json", report)
    print(f"report: {run_dir / 'report.json'}")
    return 0 if failures == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
