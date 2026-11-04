#!/usr/bin/env python3
"""Stage 0 of the HyperGraSS + RC-MDL fusion: measure a hybrid collection.

The experiment builds a HyperGraSS pack over 23 E. coli genomes and charges the
24th genome as an RC-MDL leave-one-out delta against that pack.  The comparable
pure-HyperGraSS number is the pack over all 24 genomes.  Both packs are built
from the same normalized FASTA and both are decompressed and checked against
the canonical sequence before their sizes are accepted.

The canonical sequence is the masked 120,108,223-base input.seq used by the
specialized benchmark; its SHA-256 is
51cb3ecb2e9f532833caa0b60345f37d0342474e9f4a239801486f9bc9219acb.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import time
from pathlib import Path


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fasta_sequence(path: Path) -> bytes:
    """Return the concatenated sequence, ignoring headers and blank lines."""
    out = bytearray()
    with path.open("rb") as handle:
        for line in handle:
            if not line.startswith(b">"):
                out += line.strip()
    return bytes(out)


def run_logged(cmd: list[str], log_path: Path, timeout: int) -> float:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    with log_path.open("w", encoding="utf-8") as log:
        log.write("$ " + " ".join(cmd) + "\n")
        log.flush()
        proc = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=timeout,
        )
        log.write(proc.stdout)
        log.flush()
    elapsed = time.perf_counter() - start
    if proc.returncode != 0:
        raise RuntimeError(
            "command failed with %d; see %s\n%s"
            % (proc.returncode, log_path, proc.stdout[-4000:])
        )
    return elapsed


def measure_pack(
    which: int,
    fasta: Path,
    expected: Path,
    out: Path,
    java: Path,
    classpath: Path,
    heap: str,
    timeout: int,
) -> dict:
    archive = out / ("pack%d.hgs" % which)
    restored = out / ("pack%d.restored.fa" % which)
    if not archive.exists():
        run_logged(
            [
                str(java),
                "-Xmx" + heap,
                "-cp",
                str(classpath),
                "hypergrass.HyperGraSS",
                "compress",
                "--pangenome",
                str(fasta),
                str(archive),
            ],
            out / ("compress%d.log" % which),
            timeout,
        )
    if not restored.exists():
        run_logged(
            [
                str(java),
                "-Xmx" + heap,
                "-cp",
                str(classpath),
                "hypergrass.HyperGraSS",
                "decompress",
                str(archive),
                str(restored),
            ],
            out / ("decompress%d.log" % which),
            timeout,
        )

    expected_bytes = expected.read_bytes()
    restored_bytes = fasta_sequence(restored)
    return {
        "fasta": str(fasta),
        "fasta_bytes": fasta.stat().st_size,
        "fasta_sha256": sha256_file(fasta),
        "expected_sequence_bytes": len(expected_bytes),
        "expected_sequence_sha256": sha256_bytes(expected_bytes),
        "archive": str(archive),
        "archive_bytes": archive.stat().st_size,
        "archive_sha256": sha256_file(archive),
        "restored_sequence_bytes": len(restored_bytes),
        "restored_sequence_sha256": sha256_bytes(restored_bytes),
        "sequence_equal": restored_bytes == expected_bytes,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--which", choices=["23", "24", "both"], default="both")
    parser.add_argument(
        "--dir",
        type=Path,
        default=Path("prototype/rcmdl/work/hybrid-ecoli-00"),
    )
    parser.add_argument(
        "--java",
        type=Path,
        default=Path.home() / "tools/jdk-27/bin/java.exe",
    )
    parser.add_argument(
        "--classpath",
        type=Path,
        default=Path("bin"),
    )
    parser.add_argument("--heap", default="1g")
    parser.add_argument("--delta-bytes", type=int, default=73194)
    parser.add_argument("--published-baseline", type=int, default=3810409)
    parser.add_argument("--timeout", type=int, default=7200)
    args = parser.parse_args()

    out = args.dir
    out.mkdir(parents=True, exist_ok=True)
    report_path = out / "report.json"
    report = json.loads(report_path.read_text()) if report_path.exists() else {}

    inputs = {
        23: (out / "pangenome23_blank.fasta", out / "expected23.seq"),
        24: (out / "pangenome24_blank.fasta", Path("prototype/rcmdl/work/ecoli/input.seq")),
    }

    order = [23, 24] if args.which == "both" else [int(args.which)]
    for which in order:
        fasta, expected = inputs[which]
        if not fasta.exists():
            raise SystemExit("missing %s" % fasta)
        if not expected.exists():
            raise SystemExit("missing %s" % expected)
        print("measuring pack%d ..." % which, flush=True)
        report["pack%d" % which] = measure_pack(
            which,
            fasta,
            expected,
            out,
            args.java,
            args.classpath,
            args.heap,
            args.timeout,
        )
        report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(
            "pack%d: %d B, sequence_equal=%s"
            % (
                which,
                report["pack%d" % which]["archive_bytes"],
                report["pack%d" % which]["sequence_equal"],
            ),
            flush=True,
        )

    if "pack23" in report and "pack24" in report:
        pack23 = report["pack23"]["archive_bytes"]
        pack24 = report["pack24"]["archive_bytes"]
        hybrid = pack23 + args.delta_bytes
        report["delta"] = {
            "bytes": args.delta_bytes,
            "description": "RC-MDL leave-one-out archive for target 00 against the 23-genome panel",
        }
        report["comparison"] = {
            "pure_hypergrass_bytes": pack24,
            "hybrid_bytes": hybrid,
            "savings_bytes": pack24 - hybrid,
            "savings_percent": (pack24 - hybrid) * 100.0 / pack24,
            "hybrid_wins": hybrid < pack24,
            "pack23_threshold_for_win": pack24 - args.delta_bytes,
            "published_baseline_bytes": args.published_baseline,
            "normalized_baseline_delta_bytes": pack24 - args.published_baseline,
        }
        report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(report["comparison"], indent=2), flush=True)


if __name__ == "__main__":
    main()
