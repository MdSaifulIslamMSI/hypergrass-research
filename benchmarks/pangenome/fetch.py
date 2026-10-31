"""Build a reproducible >100 MB multi-strain E. coli pangenome from NCBI.

The selection is deterministic for a fixed NCBI snapshot: complete RefSeq
E. coli records are returned in relevance order, then selected in order until
the summed sequence length reaches TARGET_BASES. Accession versions, lengths,
and SHA-256 hashes are written to manifest.json; the FASTA itself is ignored by
git. Re-running against a changed NCBI snapshot legitimately changes the lock.

Usage:
    python -B benchmarks/pangenome/fetch.py [--target-bytes N] [--force]
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
QUERY = '"Escherichia coli"[Organism] AND "complete genome"[Title] AND refseq[Filter]'
TARGET_BASES = 120_000_000
MAX_IDS = 140
BATCH = 5
INTERVAL = 0.45
TIMEOUT = 90
UA = "HyperGraSS-scale-benchmark/1.0"


def get(url: str) -> bytes:
    request = Request(url, headers={"User-Agent": UA})
    with urlopen(request, timeout=TIMEOUT) as response:
        return response.read()


def esearch() -> list[str]:
    url = f"{EUTILS}/esearch.fcgi?db=nuccore&retmode=json&retmax={MAX_IDS}&term={quote(QUERY)}"
    payload = json.loads(get(url))
    return payload["esearchresult"]["idlist"]


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
    url = f"{EUTILS}/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id={','.join(ids)}"
    return get(url)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-bytes", type=int, default=TARGET_BASES)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)

    root = Path(__file__).resolve().parent
    data = root / "data"
    fasta = data / "pangenome.fasta"
    manifest_path = root / "manifest.json"
    if fasta.exists() and not args.force:
        print(f"exists: {fasta} ({fasta.stat().st_size:,} bytes); use --force to replace")
        return 0

    ids = esearch()
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
        if row["length"] < 3_000_000 or row["length"] > 7_500_000:
            continue
        seen.add(accession)
        selected.append(row)
        total += row["length"]
        if total >= args.target_bytes:
            break
    if total < args.target_bytes:
        raise SystemExit(f"only {total:,} bases available across {len(selected)} records")
    print(f"selected: {len(selected)} records, {total:,} bases")

    data.mkdir(parents=True, exist_ok=True)
    temporary = fasta.with_suffix(".fasta.tmp")
    written = 0
    records = []
    with temporary.open("wb") as out:
        for start in range(0, len(selected), BATCH):
            batch = selected[start : start + BATCH]
            payload = efetch([row["uid"] for row in batch])
            text = payload.decode("ascii")
            if not text.startswith(">"):
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
    text = fasta.read_text(encoding="ascii")
    observed_records = text.count(">")
    observed_bases = sum(len(line.strip()) for line in text.splitlines() if not line.startswith(">"))
    if observed_records != len(records):
        raise SystemExit(f"record count mismatch: {observed_records} != {len(records)}")
    if observed_bases != total:
        raise SystemExit(f"base count mismatch: {observed_bases:,} != {total:,}")
    digest = hashlib.sha256(fasta.read_bytes()).hexdigest()
    size = fasta.stat().st_size
    if size < 100_000_000:
        raise SystemExit(f"pangenome is only {size:,} bytes; target was 100 MB+")
    manifest = {
        "schema_version": 1,
        "source": "NCBI Nucleotide (RefSeq complete E. coli genomes)",
        "query": QUERY,
        "selection": "NCBI relevance order; distinct accessions; 3-7.5 Mbp; stop at target bases",
        "target_bases": args.target_bytes,
        "total_bases": total,
        "fasta": {
            "path": "data/pangenome.fasta",
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
