"""Independent reconstruction verifier for HyperGraSS benchmark evidence.

Runs as its own process and never imports the codec. It streams both files in fixed-size chunks,
compares every byte, and hashes the full contents with SHA-256. The verdict is written as JSON so
benchmark reports can carry a machine-readable proof that reconstruction was byte-exact.

Usage:
  python benchmarks/verify/verify.py --original IN --restored OUT --archive ARCHIVE
      [--expect-sha256 HEX] [--report REPORT.json] [--chunk-bytes N]

Exit code 0 means: same length, identical bytes, archive present, and (when supplied) the expected
SHA-256 matched. Any other outcome exits 1 with a JSON verdict describing the first difference.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path


def sha256_file(path: Path, chunk_bytes: int) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            block = handle.read(chunk_bytes)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def compare_files(left: Path, right: Path, chunk_bytes: int) -> tuple[bool, int | None, int]:
    """Returns (identical, first_mismatch_offset, bytes_compared)."""
    compared = 0
    with open(left, "rb") as a, open(right, "rb") as b:
        while True:
            block_a = a.read(chunk_bytes)
            block_b = b.read(chunk_bytes)
            if block_a != block_b:
                if len(block_a) != len(block_b):
                    common = min(len(block_a), len(block_b))
                    for i in range(common):
                        if block_a[i] != block_b[i]:
                            return False, compared + i, compared + common
                    return False, compared + common, compared + common
                for i in range(len(block_a)):
                    if block_a[i] != block_b[i]:
                        return False, compared + i, compared
                return False, compared, compared
            if not block_a:
                break
            compared += len(block_a)
    return True, None, compared


def verify(original: Path, restored: Path, archive: Path, expect_sha256: str | None,
           chunk_bytes: int) -> dict:
    started = time.time()
    problems: list[str] = []

    for label, path in (("original", original), ("restored", restored), ("archive", archive)):
        if not path.exists():
            problems.append(f"{label} missing: {path}")
    if problems:
        return {
            "verdict": "fail",
            "problems": problems,
            "seconds": round(time.time() - started, 3),
        }

    original_bytes = original.stat().st_size
    restored_bytes = restored.stat().st_size
    archive_bytes = archive.stat().st_size

    identical, mismatch_at, compared = compare_files(original, restored, chunk_bytes)
    original_sha256 = sha256_file(original, chunk_bytes)
    restored_sha256 = sha256_file(restored, chunk_bytes)
    archive_sha256 = sha256_file(archive, chunk_bytes)

    if original_bytes != restored_bytes:
        problems.append(f"length differs: original {original_bytes} != restored {restored_bytes}")
    if not identical:
        problems.append(f"bytes differ at offset {mismatch_at}")
    if original_sha256 != restored_sha256:
        problems.append("SHA-256 differs")
    if expect_sha256 is not None and original_sha256.lower() != expect_sha256.lower():
        problems.append("expected input SHA-256 does not match the original file")
    if archive_bytes == 0:
        problems.append("archive is empty")

    verdict = {
        "verdict": "pass" if not problems else "fail",
        "problems": problems,
        "original": str(original),
        "restored": str(restored),
        "archive": str(archive),
        "original_bytes": original_bytes,
        "restored_bytes": restored_bytes,
        "archive_bytes": archive_bytes,
        "bytes_compared": compared,
        "byte_identical": identical,
        "first_mismatch_offset": mismatch_at,
        "original_sha256": original_sha256,
        "restored_sha256": restored_sha256,
        "archive_sha256": archive_sha256,
        "expected_sha256": expect_sha256,
        "chunk_bytes": chunk_bytes,
        "seconds": round(time.time() - started, 3),
        "comparison": "streaming byte comparison independent of the codec",
    }
    return verdict


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original", type=Path, required=True)
    parser.add_argument("--restored", type=Path, required=True)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--expect-sha256", default=None)
    parser.add_argument("--report", type=Path, default=None)
    parser.add_argument("--chunk-bytes", type=int, default=8 << 20)
    args = parser.parse_args(argv)

    verdict = verify(args.original, args.restored, args.archive, args.expect_sha256,
                     args.chunk_bytes)
    payload = json.dumps(verdict, indent=2) + "\n"
    if args.report is not None:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(payload, encoding="utf-8")
    sys.stdout.write(payload)
    return 0 if verdict["verdict"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
