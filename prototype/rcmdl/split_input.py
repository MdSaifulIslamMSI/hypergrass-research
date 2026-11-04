#!/usr/bin/env python3
"""Split the canonical concatenated input.seq into per-genome .seq files."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    records = json.loads(args.manifest.read_text(encoding="utf-8"))["records"]
    payload = args.input.read_bytes()
    args.out.mkdir(parents=True, exist_ok=True)
    offset = 0
    for index, record in enumerate(records):
        length = int(record["length"])
        chunk = payload[offset : offset + length]
        if len(chunk) != length:
            raise SystemExit(f"record {index} truncated: {len(chunk)} != {length}")
        (args.out / f"{index:02d}.seq").write_bytes(chunk)
        offset += length
    if offset != len(payload):
        raise SystemExit(f"unconsumed input bytes: {len(payload) - offset}")
    print(f"wrote {len(records)} sequences to {args.out} ({offset} bytes)")


if __name__ == "__main__":
    main()
