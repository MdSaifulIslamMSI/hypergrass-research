"""Re-fetch the exact accessions locked in an independent-pangenome manifest.

The original fetch.py selects records by NCBI relevance order, which is not a
permanent snapshot.  This script fetches only the accession versions recorded in
manifest.json, in manifest order, and fails unless the resulting FASTA matches
the locked byte size and SHA-256.
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


def uid_for(accession: str) -> str:
    term = quote(f"{accession}[accn]")
    url = f"{EUTILS}/esearch.fcgi?db=nuccore&retmode=json&retmax=5&term={term}"
    ids = json.loads(get(url))["esearchresult"]["idlist"]
    if not ids:
        raise SystemExit(f"no NCBI UID for {accession}")
    return ids[0]


def efetch(uids: list[str]) -> bytes:
    url = (
        f"{EUTILS}/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text"
        f"&id={','.join(uids)}"
    )
    return get(url)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    root = Path(__file__).resolve().parent
    dataset = root / "data" / args.dataset
    manifest_path = dataset / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    fasta = dataset / manifest["fasta"]["path"]
    if fasta.exists() and not args.force:
        print(f"exists: {fasta} ({fasta.stat().st_size:,} bytes)")
        return

    accessions = [row["accession"] for row in manifest["records"]]
    uids = []
    for i, accession in enumerate(accessions):
        uids.append(uid_for(accession))
        print(f"uid {i + 1}/{len(accessions)}: {accession} -> {uids[-1]}", flush=True)
        time.sleep(INTERVAL)

    temporary = fasta.with_suffix(".fasta.tmp")
    written = 0
    with temporary.open("wb") as out:
        for start in range(0, len(uids), BATCH):
            payload = efetch(uids[start : start + BATCH])
            if not payload.startswith(b">"):
                raise SystemExit(f"unexpected efetch payload at batch {start}")
            out.write(payload if payload.endswith(b"\n") else payload + b"\n")
            written += len(payload)
            print(
                f"  fetched {min(start + BATCH, len(uids))}/{len(uids)} ({written:,} bytes)",
                flush=True,
            )
            time.sleep(INTERVAL)
    temporary.replace(fasta)

    data = fasta.read_bytes()
    digest = hashlib.sha256(data).hexdigest()
    expected = manifest["fasta"]
    print(f"{fasta}: {len(data):,} bytes, sha256 {digest}")
    if len(data) != expected["size_bytes"]:
        raise SystemExit(f"size mismatch: {len(data):,} != {expected['size_bytes']:,}")
    if digest != expected["sha256"]:
        raise SystemExit(f"sha256 mismatch: {digest} != {expected['sha256']}")
    print("locked manifest verified")


if __name__ == "__main__":
    main()
