#!/usr/bin/env python3
"""Optional stdlib-only archive-size check; compile HyperGraSS externally first.

From the project root: python benchmarks/panman/run.py
--classpath selects an existing build (relative entries resolve from the project
root). Each run replaces report.json, data/sars_20.hgs,
data/sars_20.restored.fa, and the two PanMAN outputs; command logs are retained
in a unique logs/ subdirectory. Only missing, small pinned fixtures are fetched.
No tools, images, or dependencies are installed or built.
"""

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import subprocess
import sys
import time
import urllib.request
import uuid


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
SOURCE_COMMIT = "1c3c6204134af94c10a447fe7c6182636c19e31a"
SOURCE_URL = "https://raw.githubusercontent.com/TurakhiaLab/panman/" + SOURCE_COMMIT + "/test/"
IMAGE = "quay.io/biocontainers/panman@sha256:ae2e1ab9f4ef57d3bef248e5580f95a615a5aff9bb53be05cac03a0c56ea0a27"
FIXTURES = {
    "sars_20.fa": (609255, "7a5cf9dffcd22d03ab6f50e02e2750ae9f58bacb328440cedf4b068c013cbefe"),
    "sars_20.json": (114274, "6d1390e22d73b81d007d0368acd8a812bcf7c05c35470595385c93be0c605e5d"),
    "sars_20.nwk": (1075, "f6fcf515751f0018988fc25673f12c86ac96b9c02b620eca055184974126d660"),
}
SEQUENCE_BYTES = frozenset(b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz-*.")
LIMITATIONS = [
    "Sequence-storage archive sizes only, not biological analysis or feature parity.",
    "PanMAN retains tree/alignment/evolutionary structure that HGS does not.",
    "PanMAN uses the supplied graph JSON and tree; generating these inputs is excluded.",
    "FASTA bytes are the size-ratio denominator for both archives, despite different stored information.",
    "Host Java versus Linux Docker wall times cannot support speed claims; Docker startup is included.",
    "Single-run wall times also include process/JVM startup; there are no warmup or throughput claims.",
    "Docker CPU/memory limits are not peak-memory measurements; peak memory is not measured.",
    "Slice verification exercises only the HyperGraSS CLI, not PanMAN random access.",
    "Results apply only to this pinned 20-record demo, not larger collections or other inputs.",
]


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def file_info(path):
    data = path.read_bytes()
    return {"path": str(path), "size_bytes": len(data), "sha256": sha256(data)}


def parse_fasta(path):
    """Strict ASCII FASTA: unique first-token IDs and nonempty records.

    Allow blank separators between records or at EOF, as emitted by PanMAN.
    Preserve complete header bytes and sequence case/symbols; remove only line
    endings and concatenate sequence lines. Letters and conventional -*. symbols
    are accepted, without interpreting or normalizing them.
    """
    records = {}
    identifier = None
    header = None
    sequence = bytearray()
    separator = False

    def finish():
        if identifier is not None:
            if not sequence:
                raise ValueError(f"{path}: empty sequence for {identifier!r}")
            records[identifier] = (header, bytes(sequence))

    with path.open("rb") as stream:
        for number, raw in enumerate(stream, 1):
            line = raw.removesuffix(b"\n").removesuffix(b"\r")
            if not line:
                if identifier is None or not sequence:
                    raise ValueError(f"{path}:{number}: blank line before sequence")
                separator = True
                continue
            if line.startswith(b">"):
                finish()
                header = line[1:]
                if (not header or header[:1] in (b" ", b"\t")
                        or any(c != 9 and not 32 <= c <= 126 for c in header)):
                    raise ValueError(f"{path}:{number}: invalid ASCII header")
                identifier = header.split()[0]
                if identifier in records:
                    raise ValueError(f"{path}:{number}: duplicate ID {identifier!r}")
                sequence = bytearray()
                separator = False
            else:
                if separator or identifier is None or not set(line) <= SEQUENCE_BYTES:
                    raise ValueError(f"{path}:{number}: invalid FASTA sequence line")
                sequence.extend(line)
    finish()
    if not records:
        raise ValueError(f"{path}: no FASTA records")
    return records


def check_records(expected, actual):
    if expected.keys() != actual.keys():
        raise ValueError("PanMAN extraction record IDs differ from the pinned FASTA")
    for identifier, (header, sequence) in expected.items():
        other_header, other_sequence = actual[identifier]
        if header != other_header:
            raise ValueError(f"PanMAN full header mismatch for {identifier!r}")
        if sequence != other_sequence:
            raise ValueError(f"PanMAN case-sensitive sequence mismatch for {identifier!r}")


def fixtures(report, timeout):
    for name, (size, digest) in FIXTURES.items():
        path = HERE / "data" / name
        url = SOURCE_URL + name
        downloaded = not path.exists()
        if downloaded:
            request = urllib.request.Request(url, headers={"User-Agent": "HyperGraSS-size-benchmark"})
            with urllib.request.urlopen(request, timeout=timeout) as response:
                data = response.read(size + 1)  # Bound downloads to the pinned size, plus one byte.
            if len(data) != size or sha256(data) != digest:
                raise ValueError(f"Downloaded fixture failed pinned size/SHA-256 check: {name}")
            path.write_bytes(data)
        info = file_info(path)
        report["fixtures"][name] = dict(info, source_url=url, downloaded=downloaded)
        if info["size_bytes"] != size or info["sha256"] != digest:
            raise ValueError(f"Fixture failed pinned size/SHA-256 check: {name}; not replacing it")


def build_identity(classpath):
    """Hash actual compiled inputs, not mutable source files or an invented version."""
    entries = []
    for entry in classpath.split(os.pathsep):
        path = (ROOT / (entry or ".")).resolve()
        if path.is_dir():
            files = sorted(path.rglob("*.class"))
            if not files:
                raise ValueError(f"No compiled classes in {path}; compile externally first")
            manifest = [dict(file_info(f), path=f.relative_to(path).as_posix()) for f in files]
            entries.append({"path": str(path), "class_files": manifest})
        elif path.is_file():
            entries.append(file_info(path))
        else:
            raise ValueError(f"Missing classpath entry: {path} (wildcards are not supported)")
    return entries


class Commands:
    def __init__(self, report, logs, timeout):
        self.report = report
        self.logs = logs
        self.timeout = timeout

    def run(self, name, argv, timeout=None):
        stem = f"{len(self.report['commands']) + 1:03d}-{name}"
        stdout = self.logs / (stem + ".stdout.log")
        stderr = self.logs / (stem + ".stderr.log")
        entry = {
            "name": name, "argv": argv, "cwd": str(ROOT),
            "timeout_seconds": self.timeout if timeout is None else timeout,
            "stdout_log": str(stdout.relative_to(HERE)),
            "stderr_log": str(stderr.relative_to(HERE)),
            "returncode": None,
        }
        self.report["commands"].append(entry)
        started = time.perf_counter()
        try:
            with stdout.open("wb") as out, stderr.open("wb") as err:
                result = subprocess.run(
                    argv, cwd=ROOT, stdin=subprocess.DEVNULL, stdout=out, stderr=err,
                    timeout=entry["timeout_seconds"], check=False,
                )
            entry["returncode"] = result.returncode
            if result.returncode:
                raise RuntimeError(f"{name} exited {result.returncode}; see {stderr}")
        except (OSError, subprocess.TimeoutExpired, RuntimeError, KeyboardInterrupt) as exc:
            entry["error"] = str(exc) or type(exc).__name__
            raise
        finally:
            entry["wall_seconds"] = time.perf_counter() - started
            (self.logs / (stem + ".command.json")).write_text(
                json.dumps(entry, indent=2) + "\n", encoding="utf-8",
            )
        return stdout.read_bytes()

    def panman(self, name, arguments):
        container = "hypergrass-panman-" + self.logs.name.lower() + "-" + name
        argv = [
            "docker", "run", "--rm", "--pull=never", "--network", "none",
            "--cpus", "2", "--memory", "2g", "--name", container,
            "--mount", f"type=bind,source={HERE.as_posix()},target=/work",
            "--workdir", "/work", IMAGE, "panmanUtils", *arguments,
        ]
        try:
            return self.run(name, argv)
        except (Exception, KeyboardInterrupt):
            # Killing a timed-out Docker client alone need not stop its container.
            try:
                self.run(name + "-cleanup", ["docker", "rm", "-f", container], timeout=30)
            except (Exception, KeyboardInterrupt):
                pass  # Preserve the original failure; cleanup failure is logged separately.
            raise


def run_benchmark(args, report, commands):
    fixtures(report, min(args.timeout, 30))
    source = HERE / "data" / "sars_20.fa"
    records = parse_fasta(source)
    if len(records) != 20:
        raise ValueError(f"Expected 20 FASTA records, found {len(records)}")
    report["records"] = [
        {"id": identifier.decode("ascii"), "header": header.decode("ascii"),
         "sequence_bytes": len(sequence), "sequence_sha256": sha256(sequence)}
        for identifier, (header, sequence) in records.items()
    ]
    build = build_identity(args.classpath)
    report["versions"]["hypergrass"] = {
        "version": None, "note": "CLI has no version option; actual compiled inputs are hashed",
        "classpath": args.classpath, "relative_entries_base": str(ROOT), "build_inputs": build,
    }
    java = commands.run("java-version", ["java", "-version"])
    java_stderr = HERE / report["commands"][-1]["stderr_log"]
    report["versions"]["java"] = (java + java_stderr.read_bytes()).decode("utf-8", "replace").strip()
    report["versions"]["docker"] = commands.run("docker-version", ["docker", "version"]).decode("utf-8", "replace").strip()
    image = json.loads(commands.run("image-inspect", ["docker", "image", "inspect", IMAGE]))[0]
    report["versions"]["panman_image"] = {
        "reference": IMAGE, "id": image["Id"], "repo_digests": image.get("RepoDigests", []),
        "os": image["Os"], "architecture": image["Architecture"],
    }
    if image["Os"] != "linux":
        raise ValueError("The pinned PanMAN image must be a Linux image")
    version = commands.panman("panman-version", ["--version"]).decode("utf-8", "replace").strip()
    report["versions"]["panman"] = version
    if not re.search(r"\b0\.1\.4\b", version):
        raise ValueError(f"Expected PanMAN 0.1.4, got {version!r}")

    archive = HERE / "panman" / "sars_20.panman"
    extracted = HERE / "info" / "sars_20_0.fasta"
    archive.unlink(missing_ok=True)
    commands.panman("panman-create", [
        "-P", "/work/data/sars_20.json", "-N", "/work/data/sars_20.nwk",
        "-o", "sars_20", "--threads", "2",
    ])
    report["artifacts"]["panman"] = file_info(archive)
    if not report["artifacts"]["panman"]["size_bytes"]:
        raise ValueError("PanMAN produced an empty archive")
    extracted.unlink(missing_ok=True)
    commands.panman("panman-extract", [
        "-I", "/work/panman/sars_20.panman", "--fasta", "-o", "sars_20", "--threads", "2",
    ])
    report["artifacts"]["panman_extracted"] = file_info(extracted)
    check_records(records, parse_fasta(extracted))
    report["verification"]["panman_full_header_sequence_map"] = {"passed": True, "records": 20}

    hgs = HERE / "data" / "sars_20.hgs"
    restored = HERE / "data" / "sars_20.restored.fa"
    java_cli = ["java", "-cp", args.classpath, "hypergrass.HyperGraSS"]
    hgs.unlink(missing_ok=True)
    compression_options = ["--optimize"] if args.optimize else []
    report["versions"]["hypergrass"]["compression_options"] = compression_options
    commands.run("hgs-compress", [*java_cli, "compress", *compression_options, str(source), str(hgs)])
    report["artifacts"]["hgs"] = file_info(hgs)
    if not report["artifacts"]["hgs"]["size_bytes"]:
        raise ValueError("HyperGraSS produced an empty archive")
    restored.unlink(missing_ok=True)
    commands.run("hgs-decompress", [*java_cli, "decompress", str(hgs), str(restored)])
    restored_info = file_info(restored)
    report["artifacts"]["hgs_restored"] = restored_info
    original_info = report["fixtures"]["sars_20.fa"]
    if (restored_info["sha256"] != original_info["sha256"]
            or restored.read_bytes() != source.read_bytes()):
        raise ValueError("HyperGraSS decompression failed byte-exact/SHA-256 verification")
    report["verification"]["hgs_byte_exact_round_trip"] = {"passed": True, "sha256": restored_info["sha256"]}

    queries = report["verification"]["hgs_slices"] = {"passed": False, "checks": []}
    for index, (identifier, (_, sequence)) in enumerate(records.items()):
        for start in (0, len(sequence) // 2, max(0, len(sequence) - 128)):
            end = min(len(sequence), start + 128)
            expected = sequence[start:end]
            got = commands.run(f"hgs-query-{index}-{start}", [
                *java_cli, "query", str(hgs), str(index), str(start), str(end),
            ])
            # The CLI appends exactly one host-platform newline, not FASTA formatting.
            if got not in (expected + b"\n", expected + b"\r\n"):
                raise ValueError(f"HyperGraSS query mismatch: record {index}, [{start}, {end})")
            queries["checks"].append({
                "record_index": index, "id": identifier.decode("ascii"),
                "start_inclusive": start, "end_exclusive": end, "sha256": sha256(expected),
            })
    queries.update(passed=True, records=20, count=len(queries["checks"]))
    if build_identity(args.classpath) != build:
        raise ValueError("Compiled classpath changed during the run; results are not attributable to one build")
    # Recheck inputs after both tools ran, rather than assuming the shared mount was unchanged.
    for name, (size, digest) in FIXTURES.items():
        info = file_info(HERE / "data" / name)
        if info["size_bytes"] != size or info["sha256"] != digest:
            raise ValueError(f"Fixture changed during the run: {name}")
    raw_size = original_info["size_bytes"]
    hgs_size = report["artifacts"]["hgs"]["size_bytes"]
    panman_size = report["artifacts"]["panman"]["size_bytes"]
    report["archive_size_comparison"] = {
        "fasta_bytes": raw_size, "hgs_bytes": hgs_size, "panman_bytes": panman_size,
        "fasta_over_hgs": raw_size / hgs_size, "fasta_over_panman": raw_size / panman_size,
        "panman_over_hgs": panman_size / hgs_size,
        "hgs_bytes_saved_vs_panman": panman_size - hgs_size,
        "hgs_percent_smaller_than_panman": 100 * (panman_size - hgs_size) / panman_size,
    }


def argument_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--classpath", default=str(ROOT / "bin"), help="Existing Java classpath; relative entries use the project root")
    parser.add_argument("--timeout", type=float, default=300, help="Per-command timeout in seconds (default: 300); fixture requests use at most 30")
    parser.add_argument("--optimize", action="store_true", help="Enable HyperGraSS's bounded extra archive-size search (more CPU)")
    return parser


def main():
    parser = argument_parser()
    args = parser.parse_args()
    if not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error("--timeout must be finite and positive")
    run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "-" + uuid.uuid4().hex[:8]
    logs = HERE / "logs" / run_id
    report_path = HERE / "report.json"
    report = {
        "schema_version": 1, "status": "running", "run_id": run_id,
        "started_utc": datetime.now(timezone.utc).isoformat(), "source_commit": SOURCE_COMMIT,
        "versions": {"python": sys.version, "host": platform.platform()},
        "sequence_contract": {
            "panman": "Exact full-header and case-sensitive sequence-byte map for all 20 records; record order and line wrapping may differ",
            "hypergrass": "Exact original FASTA bytes, including headers, order, wrapping, and line endings",
            "transformations": "None: no ambiguity normalization, case folding, gap removal, trimming, reverse complement, or circular rotation",
        },
        "limits": {"panman_cpus": 2, "panman_memory": "2g", "panman_threads": 2,
                   "network": "none in Docker; HTTPS only for missing pinned fixtures", "peak_memory_measured": False},
        "limitations": LIMITATIONS, "fixtures": {}, "artifacts": {}, "verification": {}, "commands": [],
    }
    try:
        for directory in (logs, HERE / "data", HERE / "panman", HERE / "info"):
            directory.mkdir(parents=True, exist_ok=True)
        # Invalidate an older success report before attempting any external command.
        report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        run_benchmark(args, report, Commands(report, logs, args.timeout))
        report["status"] = "success"
    except (Exception, KeyboardInterrupt) as exc:
        report["status"] = "failed"
        report["error"] = f"{type(exc).__name__}: {exc}"
        print(report["error"], file=sys.stderr)
    report["finished_utc"] = datetime.now(timezone.utc).isoformat()
    try:
        report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    except OSError as exc:
        print(f"Cannot write report {report_path}: {exc}", file=sys.stderr)
        return 1
    print(f"{report['status']}: {report_path}")
    if report["status"] == "success":
        sizes = report["archive_size_comparison"]
        print(f"FASTA {sizes['fasta_bytes']} B; HGS {sizes['hgs_bytes']} B; PanMAN {sizes['panman_bytes']} B")
    return 0 if report["status"] == "success" else 1


if __name__ == "__main__":
    sys.exit(main())
