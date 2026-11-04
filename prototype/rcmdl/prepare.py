#!/usr/bin/env python3
"""Extract the 24 E. coli genomes from the pangenome FASTA into flat .seq files."""
from __future__ import annotations

import argparse
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
                chunks.append(line.strip())
    if name is not None:
        yield name, b"".join(chunks)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fasta",
        default="benchmarks/pangenome/data/pangenome.fasta",
        type=Path,
    )
    parser.add_argument(
        "--out",
        default="prototype/rcmdl/work/genomes",
        type=Path,
    )
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    index = []
    for i, (name, seq) in enumerate(parse_fasta(args.fasta)):
        path = args.out / f"{i:02d}.seq"
        path.write_bytes(seq)
        index.append((i, name, len(seq), path.name))
    with (args.out / "index.tsv").open("w", encoding="utf-8") as handle:
        for row in index:
            handle.write("\t".join(map(str, row)) + "\n")
    print(f"wrote {len(index)} genomes to {args.out}")


if __name__ == "__main__":
    main()
