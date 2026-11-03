"""Fetch a pool of complete genomes for the broad pangenome experiment.

Unlike ``benchmarks/independent/fetch.py``, this keeps a larger, randomly
sampled pool per species so that the same pool can supply both a natural
~120 MB collection and low/medium/high diversity strata. Records are stored
verbatim with byte offsets, which lets ``select.py`` write subsets without
re-parsing or modifying the source.

Usage:
    python -B benchmarks/broad/fetch.py \
        --name staph --organism "Staphylococcus aureus" --pool-size 72
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
import time
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request, urlopen

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
BATCH = 5
SUMMARY_BATCH = 50
INTERVAL = 0.45
TIMEOUT = 180
UA = "HyperGraSS-broad-benchmark/1.0"


def get(url: str) -> bytes:
    request = Request(url, headers={"User-Agent": UA})
    with urlopen(request, timeout=TIMEOUT) as response:
        return response.read()


def esearch(query: str, max_ids: int) -> list[str]:
    url = (
        f"{EUTILS}/esearch.fcgi?db=nuccore&retmode=json"
        f"&retmax={max_ids}&term={quote(query)}"
    )
    return json.loads(get(url))["esearchresult"]["idlist"]


def esummary(ids: list[str]) -> list[dict]:
    rows: list[dict] = []
    for start in range(0, len(ids), SUMMARY_BATCH):
        batch = ids[start : start + SUMMARY_BATCH]
        url = (
            f"{EUTILS}/esummary.fcgi?db=nuccore&retmode=json"
            f"&id={','.join(batch)}"
        )
        payload = json.loads(get(url))
        result = payload["result"]
        for uid in result["uids"]:
            item = result[uid]
            rows.append(
                {
                    "uid": uid,
                    "accession": item.get("accessionversion", ""),
                    "length": int(item.get("slen", 0)),
                    "title": item.get("title", ""),
                }
            )
        time.sleep(INTERVAL)
    return rows


def efetch(ids: list[str]) -> bytes:
    url = (
        f"{EUTILS}/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text"
        f"&id={','.join(ids)}"
    )
    return get(url)


def split_fasta(payload: bytes) -> list[bytes]:
    """Split a FASTA payload into verbatim records beginning with '>'."""
    payload = payload.replace(b"\r\n", b"\n")
    records: list[bytes] = []
    current = bytearray()
    for line in payload.splitlines(keepends=True):
        if line.startswith(b">"):
            if current:
                records.append(bytes(current))
            current = bytearray()
        if current or line.startswith(b">"):
            current.extend(line if line.endswith(b"\n") else line + b"\n")
    if current:
        records.append(bytes(current))
    return records


def record_facts(record: bytes) -> dict:
    header = b""
    sequence = bytearray()
    for line in record.splitlines():
        if line.startswith(b">"):
            header = line[1:].strip()
        else:
            sequence.extend(line.strip())
    sequence_bytes = bytes(sequence).upper()
    non_acgt = sum(1 for value in sequence_bytes if value not in b"ACGT")
    gc = sequence_bytes.count(b"G") + sequence_bytes.count(b"C")
    accession = header.split(b" ")[0].decode("ascii", "replace")
    title = header[len(header.split(b" ")[0]) :].strip().decode("ascii", "replace")
    return {
        "accession": accession,
        "header": header.decode("ascii", "replace"),
        "title": title,
        "sequence_length": len(sequence_bytes),
        "non_acgt_bases": non_acgt,
        "gc_fraction": (gc / len(sequence_bytes)) if sequence_bytes else 0.0,
        "sequence_sha256": hashlib.sha256(sequence_bytes).hexdigest(),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--name", required=True)
    parser.add_argument("--organism", required=True)
    parser.add_argument("--query", default=None)
    parser.add_argument("--pool-size", type=int, default=72)
    parser.add_argument("--min-length", type=int, default=800_000)
    parser.add_argument("--max-length", type=int, default=9_000_000)
    parser.add_argument("--max-ids", type=int, default=300)
    parser.add_argument("--seed", type=int, default=20261102)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)

    query = args.query or (
        f'"{args.organism}"[Organism] AND "complete genome"[Title] '
        f"AND refseq[Filter]"
    )
    root = Path(__file__).resolve().parent
    dataset = root / "data" / args.name
    fasta = dataset / "pangenome.fasta"
    manifest_path = dataset / "manifest.json"
    if fasta.exists() and manifest_path.exists() and not args.force:
        print(f"exists: {fasta} ({fasta.stat().st_size:,} bytes); use --force to replace")
        return 0

    ids = esearch(query, args.max_ids)
    print(f"search: {len(ids)} candidate records")
    time.sleep(INTERVAL)
    rows = esummary(ids)

    eligible = []
    seen = set()
    for row in rows:
        accession = row["accession"]
        if not accession or accession in seen:
            continue
        if not args.min_length <= row["length"] <= args.max_length:
            continue
        seen.add(accession)
        eligible.append(row)
    if not eligible:
        raise SystemExit("no eligible complete genomes after length filtering")

    rng = random.Random(args.seed)
    if len(eligible) <= args.pool_size:
        selected = eligible
        selection = "all eligible records"
    else:
        selected = rng.sample(eligible, args.pool_size)
        selection = (
            f"uniform random sample of {args.pool_size} eligible records "
            f"with seed {args.seed}"
        )
    selected.sort(key=lambda row: row["accession"])
    print(
        f"eligible: {len(eligible)} records; selected: {len(selected)} "
        f"({sum(row['length'] for row in selected):,} bp)"
    )

    dataset.mkdir(parents=True, exist_ok=True)
    temporary = fasta.with_suffix(".fasta.tmp")
    records: list[dict] = []
    written = 0
    with temporary.open("wb") as out:
        for start in range(0, len(selected), BATCH):
            batch = selected[start : start + BATCH]
            payload = efetch([row["uid"] for row in batch])
            if not payload.startswith(b">"):
                raise SystemExit(f"unexpected efetch payload for {batch[0]['accession']}")
            for record in split_fasta(payload):
                facts = record_facts(record)
                if not args.min_length <= facts["sequence_length"] <= args.max_length:
                    continue
                offset = out.tell()
                out.write(record)
                written += len(record)
                records.append(
                    {
                        "uid": next(
                            (row["uid"] for row in batch
                             if row["accession"] == facts["accession"]),
                            None,
                        ),
                        "accession": facts["accession"],
                        "title": facts["title"],
                        "sequence_length": facts["sequence_length"],
                        "non_acgt_bases": facts["non_acgt_bases"],
                        "gc_fraction": facts["gc_fraction"],
                        "sequence_sha256": facts["sequence_sha256"],
                        "file_offset": offset,
                        "file_length": len(record),
                    }
                )
            print(f"  fetched {min(start + len(batch), len(selected))}/{len(selected)} "
                  f"({written:,} bytes)", flush=True)
            time.sleep(INTERVAL)

    temporary.replace(fasta)
    data = fasta.read_bytes()
    observed_records = sum(1 for line in data.splitlines() if line.startswith(b">"))
    observed_bases = sum(
        len(line.strip()) for line in data.splitlines() if not line.startswith(b">")
    )
    if observed_records != len(records):
        raise SystemExit(f"record count mismatch: {observed_records} != {len(records)}")
    if observed_bases != sum(row["sequence_length"] for row in records):
        raise SystemExit("base count mismatch")
    if not records:
        raise SystemExit("no records written")

    digest = hashlib.sha256(data).hexdigest()
    manifest = {
        "schema_version": 1,
        "kind": "broad-pool",
        "name": args.name,
        "organism": args.organism,
        "source": f"NCBI Nucleotide (RefSeq complete genomes: {args.organism})",
        "query": query,
        "selection": (
            f"{selection}; distinct accessions; "
            f"{args.min_length:,}-{args.max_length:,} bp; sorted by accession"
        ),
        "seed": args.seed,
        "candidate_count": len(ids),
        "eligible_count": len(eligible),
        "pool_size": len(records),
        "total_bases": sum(row["sequence_length"] for row in records),
        "fasta": {
            "path": "pangenome.fasta",
            "size_bytes": fasta.stat().st_size,
            "sha256": digest,
        },
        "records": records,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {fasta} ({fasta.stat().st_size:,} bytes), {len(records)} records")
    print(f"wrote {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
