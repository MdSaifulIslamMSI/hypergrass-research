"""Extract one FASTA record slice from a gzip or xz stream.

This is the baseline side of the end-to-end query benchmark. It answers the
same question HyperGraSS answers directly: return bases [start, end) of record
N as sequence characters, with headers and line endings removed.

The stream is decompressed sequentially and stopped once the requested slice is
complete. A query for a late record therefore pays for most of the archive,
which is the honest cost of a non-indexed format.
"""

from __future__ import annotations

import subprocess
import sys

READ_SIZE = 1024 * 1024


def main(argv: list[str]) -> int:
    if len(argv) != 6:
        print(
            "usage: extract.py {gzip|xz} ARCHIVE RECORD START END",
            file=sys.stderr,
        )
        return 2
    tool, archive, record_text, start_text, end_text = argv[1:]
    if tool not in ("gzip", "xz"):
        print(f"unsupported tool: {tool}", file=sys.stderr)
        return 2
    record = int(record_text)
    start = int(start_text)
    end = int(end_text)
    if record < 0 or start < 0 or end <= start:
        print("invalid coordinates", file=sys.stderr)
        return 2

    process = subprocess.Popen(
        [tool, "-dc", archive],
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
    )
    assert process.stdout is not None
    current = -1
    position = 0
    out = bytearray()
    found = False
    try:
        for line in process.stdout:
            if line.startswith(b">"):
                if found:
                    break
                current += 1
                if current == record:
                    found = True
                    position = 0
                continue
            if not found:
                continue
            clean = line.rstrip(b"\r\n")
            line_start = position
            line_end = position + len(clean)
            if line_end > start and line_start < end:
                left = max(0, start - line_start)
                right = min(len(clean), end - line_start)
                out.extend(clean[left:right])
            position = line_end
            if position >= end:
                break
    finally:
        process.stdout.close()
        process.terminate()
        process.wait()

    if not found:
        print(f"record {record} not found", file=sys.stderr)
        return 3
    if len(out) != end - start:
        print(
            f"short slice: expected {end - start} bases, got {len(out)}",
            file=sys.stderr,
        )
        return 4
    sys.stdout.buffer.write(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
