#!/usr/bin/env python3
"""Split a multi-record FASTA panel into one flat .seq file per record."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def parse_fasta(path: Path):
    name = None
    chunks: list[bytes] = []
    with path.open("rb") as handle:
        for line in handle:
            if line.startswith(b">"):
                if name is not None:
                    yield name, b"".join(chunks)
                name = line[1:].strip().split(None, 1)[0].decode("ascii", "replace")
                chunks = []
            else:
                stripped = line.strip()
                if stripped:
                    chunks.append(stripped)
    if name is not None:
        yield name, b"".join(chunks)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    rows = []
    for index, (name, sequence) in enumerate(parse_fasta(args.fasta)):
        path = args.out / f"{index:02d}.seq"
        path.write_bytes(sequence)
        rows.append({"index": index, "name": name, "bases": len(sequence), "path": path.name})
    (args.out / "index.json").write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {len(rows)} sequences to {args.out}")


if __name__ == "__main__":
    main()
