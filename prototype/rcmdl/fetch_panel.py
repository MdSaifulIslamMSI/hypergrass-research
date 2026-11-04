#!/usr/bin/env python3
"""Fetch a deterministic multi-genome panel from NCBI RefSeq complete genomes."""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request, urlopen

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UA = "RC-MDL-panel-fetch/1.0"
TIMEOUT = 120


def get(url: str) -> bytes:
    last = None
    for attempt in range(4):
        try:
            request = Request(url, headers={"User-Agent": UA})
            with urlopen(request, timeout=TIMEOUT) as response:
                return response.read()
        except Exception as exc:  # network flakiness is expected
            last = exc
            time.sleep(2.0 * (attempt + 1))
    raise RuntimeError("GET failed: %s (%s)" % (url, last))


def esearch(query: str, retmax: int) -> list[str]:
    url = f"{EUTILS}/esearch.fcgi?db=nuccore&retmode=json&retmax={retmax}&term={quote(query)}"
    return json.loads(get(url))["esearchresult"]["idlist"]


def esummary(ids: list[str]) -> list[dict]:
    if not ids:
        return []
    url = f"{EUTILS}/esummary.fcgi?db=nuccore&retmode=json&id={','.join(ids)}"
    result = json.loads(get(url))["result"]
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--query", required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--target-bases", type=int, default=100_000_000)
    parser.add_argument("--min-length", type=int, default=1_000_000)
    parser.add_argument("--max-length", type=int, default=10_000_000)
    parser.add_argument("--max-ids", type=int, default=200)
    parser.add_argument("--batch", type=int, default=5)
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    fasta = args.out / "panel.fasta"
    manifest = args.out / "manifest.json"

    ids = esearch(args.query, args.max_ids)
    print("search: %d candidates" % len(ids), flush=True)
    time.sleep(0.5)
    rows = esummary(ids)
    time.sleep(0.5)

    selected = []
    total = 0
    seen = set()
    for row in rows:
        accession = row["accession"]
        if not accession or accession in seen:
            continue
        if not (args.min_length <= row["length"] <= args.max_length):
            continue
        seen.add(accession)
        selected.append(row)
        total += row["length"]
        if total >= args.target_bases:
            break
    if total < args.target_bases:
        raise SystemExit("only %d bases across %d records" % (total, len(selected)))
    print("selected: %d records, %d bases" % (len(selected), total), flush=True)

    temporary = fasta.with_suffix(".fasta.tmp")
    written = 0
    with temporary.open("wb") as out:
        for start in range(0, len(selected), args.batch):
            batch = selected[start : start + args.batch]
            payload = efetch([row["uid"] for row in batch])
            if not payload.startswith(b">"):
                raise SystemExit("unexpected efetch payload")
            out.write(payload if payload.endswith(b"\n") else payload + b"\n")
            written += len(payload)
            print(
                "  fetched %d/%d (%d bytes)" % (min(start + args.batch, len(selected)), len(selected), written),
                flush=True,
            )
            time.sleep(0.5)
    temporary.replace(fasta)

    text = fasta.read_text(encoding="ascii")
    observed_records = text.count(">")
    observed_bases = sum(len(line.strip()) for line in text.splitlines() if not line.startswith(">"))
    if observed_records != len(selected):
        raise SystemExit("record count mismatch: %d != %d" % (observed_records, len(selected)))
    if observed_bases != total:
        raise SystemExit("base count mismatch: %d != %d" % (observed_bases, total))

    digest = hashlib.sha256(fasta.read_bytes()).hexdigest()
    manifest.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "query": args.query,
                "selection": "NCBI relevance order; distinct accessions; length %d-%d; stop at target bases"
                % (args.min_length, args.max_length),
                "target_bases": args.target_bases,
                "total_bases": total,
                "fasta": {
                    "path": fasta.name,
                    "size_bytes": fasta.stat().st_size,
                    "sha256": digest,
                },
                "records": selected,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print("wrote %s (%d bytes), %s" % (fasta, fasta.stat().st_size, manifest), flush=True)


if __name__ == "__main__":
    main()
