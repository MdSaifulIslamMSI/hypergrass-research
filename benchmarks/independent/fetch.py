"""Fetch an independent multi-strain bacterial pangenome from NCBI.

This is the independent-replication counterpart to
``benchmarks/pangenome/fetch.py``. It takes an arbitrary NCBI organism query
and writes a self-describing dataset under ``benchmarks/independent/data``.

Usage:
    python -B benchmarks/independent/fetch.py \
        --name salmonella \
        --query '"Salmonella enterica"[Organism] AND "complete genome"[Title] AND refseq[Filter]'
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request, urlopen

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
BATCH = 5
INTERVAL = 0.45
TIMEOUT = 120
UA = "HyperGraSS-independent-benchmark/1.0"


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
    if not ids:
        return []
    url = f"{EUTILS}/esummary.fcgi?db=nuccore&retmode=json&id={','.join(ids)}"
    payload = json.loads(get(url))
    result = payload["result"]
    rows = []
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
    return rows


def efetch(ids: list[str]) -> bytes:
    url = (
        f"{EUTILS}/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text"
        f"&id={','.join(ids)}"
    )
    return get(url)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--name", required=True)
    parser.add_argument("--query")
    parser.add_argument("--organism")
    parser.add_argument("--target-bytes", type=int, default=120_000_000)
    parser.add_argument("--min-length", type=int, default=3_000_000)
    parser.add_argument("--max-length", type=int, default=7_500_000)
    parser.add_argument("--max-ids", type=int, default=140)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)
    if not args.query:
        if not args.organism:
            raise SystemExit("provide --query or --organism")
        args.query = (
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

    ids = esearch(args.query, args.max_ids)
    print(f"search: {len(ids)} candidate records")
    time.sleep(INTERVAL)
    rows = esummary(ids)
    time.sleep(INTERVAL)

    selected = []
    total = 0
    seen = set()
    for row in rows:
        accession = row["accession"]
        if not accession or accession in seen:
            continue
        if not args.min_length <= row["length"] <= args.max_length:
            continue
        seen.add(accession)
        selected.append(row)
        total += row["length"]
        if total >= args.target_bytes:
            break
    if total < args.target_bytes:
        raise SystemExit(f"only {total:,} bases available across {len(selected)} records")
    print(f"selected: {len(selected)} records, {total:,} bases")

    dataset.mkdir(parents=True, exist_ok=True)
    temporary = fasta.with_suffix(".fasta.tmp")
    written = 0
    records = []
    with temporary.open("wb") as out:
        for start in range(0, len(selected), BATCH):
            batch = selected[start : start + BATCH]
            payload = efetch([row["uid"] for row in batch])
            if not payload.startswith(b">"):
                raise SystemExit(f"unexpected efetch payload for {batch[0]['accession']}")
            out.write(payload if payload.endswith(b"\n") else payload + b"\n")
            written += len(payload)
            for row in batch:
                records.append(
                    {
                        "accession": row["accession"],
                        "length": row["length"],
                        "title": row["title"],
                    }
                )
            print(f"  fetched {start + len(batch)}/{len(selected)} ({written:,} bytes)")
            time.sleep(INTERVAL)

    temporary.replace(fasta)
    data = fasta.read_bytes()
    observed_records = sum(1 for line in data.splitlines() if line.startswith(b">"))
    observed_bases = sum(
        len(line.strip()) for line in data.splitlines() if not line.startswith(b">")
    )
    if observed_records != len(records):
        raise SystemExit(f"record count mismatch: {observed_records} != {len(records)}")
    if observed_bases != total:
        raise SystemExit(f"base count mismatch: {observed_bases:,} != {total:,}")

    digest = hashlib.sha256(data).hexdigest()
    size = fasta.stat().st_size
    if size < 100_000_000:
        raise SystemExit(f"pangenome is only {size:,} bytes; target was 100 MB+")
    manifest = {
        "schema_version": 1,
        "source": f"NCBI Nucleotide (RefSeq complete genomes: {args.name})",
        "query": args.query,
        "selection": (
            f"NCBI relevance order; distinct accessions; "
            f"{args.min_length}-{args.max_length} bp; stop at target bases"
        ),
        "target_bases": args.target_bytes,
        "total_bases": total,
        "fasta": {
            "path": "pangenome.fasta",
            "size_bytes": size,
            "sha256": digest,
        },
        "records": records,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {fasta} ({size:,} bytes), {len(records)} records")
    print(f"wrote {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
