"""Freeze/verify the held-out corpus using only the Python standard library.

``root`` is the directory containing protocol.json. Run this file without
arguments to prepare data and source archives, or use --data-only,
--sources-only, or --verify (both locks, offline). Existing locks and artifacts
are never replaced; delete nothing to "repair" a failed verification.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import time
from urllib.parse import urlencode, urlsplit
from urllib.request import HTTPRedirectHandler, Request, build_opener


SCHEMA_VERSION = 1
REQUEST_TIMEOUT = 30
FASTA_LIMIT = 32 * 1024 * 1024
SOURCE_LIMIT = 128 * 1024 * 1024
READ_CHUNK = 64 * 1024
REQUEST_INTERVAL = 0.35  # Slightly slower than NCBI's unauthenticated 3/s limit.
SOURCES = (
    ("jarvis3", "4dcb37b05b1a1482c2bf864d4968244817b0ad0a"),
    ("geco3", "69642a0a913e543599f35aee372e300601e4cbbe"),
)
_last_request: float | None = None


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(READ_CHUNK), b""):
            digest.update(chunk)
    return digest.hexdigest()


class _HTTPSRedirectHandler(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        if urlsplit(newurl).scheme != "https":
            raise ValueError(f"Refusing non-HTTPS redirect: {newurl}")
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def _download(url: str, max_bytes: int) -> bytes:
    """Bounded HTTPS reads: 30s socket timeout, elapsed guard, and byte cap."""
    global _last_request
    if urlsplit(url).scheme != "https":
        raise ValueError("Downloads require HTTPS")
    if _last_request is not None:
        delay = REQUEST_INTERVAL - (time.monotonic() - _last_request)
        if delay > 0:
            time.sleep(delay)
    _last_request = time.monotonic()
    deadline = _last_request + REQUEST_TIMEOUT
    request = Request(url, headers={
        "User-Agent": "hypergrass-held-out-corpus/1 (Python stdlib)",
        "Accept-Encoding": "identity",
    })
    opener = build_opener(_HTTPSRedirectHandler())
    with opener.open(request, timeout=REQUEST_TIMEOUT) as response:
        if urlsplit(response.geturl()).scheme != "https":
            raise ValueError("Download did not remain on HTTPS")
        if response.headers.get("Content-Encoding", "identity") != "identity":
            raise ValueError("Unexpected Content-Encoding; refusing transformed bytes")
        length = response.headers.get("Content-Length")
        if length is not None and not 0 <= int(length) <= max_bytes:
            raise ValueError(f"Download exceeds {max_bytes} bytes: {url}")
        chunks = []
        size = 0
        while True:
            if time.monotonic() >= deadline:
                raise TimeoutError(f"Download exceeded {REQUEST_TIMEOUT}s: {url}")
            chunk = response.read1(min(READ_CHUNK, max_bytes + 1 - size))
            if time.monotonic() >= deadline:
                raise TimeoutError(f"Download exceeded {REQUEST_TIMEOUT}s: {url}")
            if not chunk:
                break
            chunks.append(chunk)
            size += len(chunk)
            if size > max_bytes:
                raise ValueError(f"Download exceeds {max_bytes} bytes: {url}")
        if length is not None and size != int(length):
            raise ValueError(f"Truncated download: {url}")
    return b"".join(chunks)


def _read_bounded(path: Path, limit: int) -> bytes:
    with path.open("rb") as stream:
        data = stream.read(limit + 1)
    if len(data) > limit:
        raise ValueError(f"File exceeds {limit} bytes: {path}")
    return data


def _load_json(path: Path) -> dict:
    value = json.loads(_read_bounded(path, 1024 * 1024))
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def _protocol(root: Path) -> tuple[dict, str]:
    data = _read_bounded(root / "protocol.json", 1024 * 1024)
    protocol = json.loads(data)
    if not isinstance(protocol, dict) or protocol.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("Unsupported protocol schema_version")
    datasets = protocol.get("datasets")
    if not isinstance(datasets, list) or not datasets:
        raise ValueError("Protocol must contain datasets")
    seen = set()
    for spec in datasets:
        if not isinstance(spec, dict):
            raise ValueError("Invalid protocol dataset")
        dataset_id = spec.get("id")
        if not isinstance(dataset_id, str) or not re.fullmatch(r"[a-z0-9][a-z0-9_-]*", dataset_id):
            raise ValueError("Invalid protocol dataset id")
        if dataset_id in seen:
            raise ValueError(f"Duplicate dataset id: {dataset_id}")
        seen.add(dataset_id)
        if not isinstance(spec.get("label"), str):
            raise ValueError(f"Missing dataset label: {dataset_id}")
        accessions = spec.get("accessions")
        if not isinstance(accessions, list) or not accessions:
            raise ValueError(f"Missing accessions: {dataset_id}")
        for accession in accessions:
            if not isinstance(accession, str) or not re.fullmatch(r"[A-Z][A-Z0-9_]*\.[1-9][0-9]*", accession):
                raise ValueError(f"Accession must be version-pinned: {accession!r}")
        if len(set(accessions)) != len(accessions):
            raise ValueError(f"Duplicate requested accession: {dataset_id}")
    return protocol, sha256_bytes(data)


def _parse_fasta(data: bytes, accessions: list[str]) -> tuple[bytes, list[dict]]:
    """Remove only header lines and CR/LF endings, never normalize sequence."""
    records = []
    sequence = []
    seen_ids = set()
    seen_headers = set()
    for line in re.split(rb"\r\n|\n|\r", data):
        if line.startswith(b">"):
            header_bytes = line[1:]
            if not header_bytes or any(byte < 32 or byte > 126 for byte in header_bytes):
                raise ValueError("FASTA header must be nonempty printable ASCII")
            header = header_bytes.decode("ascii")
            record_id = header.partition(" ")[0]
            if not record_id:
                raise ValueError("FASTA header must start with its accession-version")
            if record_id in seen_ids or header in seen_headers:
                raise ValueError(f"Duplicate FASTA record id/header: {record_id}")
            seen_ids.add(record_id)
            seen_headers.add(header)
            records.append({"id": record_id, "header": header, "sequence_bytes": 0})
        elif line:
            if not records:
                raise ValueError("Sequence before the first FASTA header")
            if any(byte < 32 or byte > 126 for byte in line):
                raise ValueError("Sequence contains non-printable/non-ASCII bytes")
            sequence.append(line)
            records[-1]["sequence_bytes"] += len(line)
    if [record["id"] for record in records] != accessions:
        raise ValueError(f"FASTA record IDs/order do not match requested {accessions!r}")
    if any(record["sequence_bytes"] == 0 for record in records):
        raise ValueError("Empty FASTA record")
    return b"".join(sequence), records


def _file_entry(path: str, data: bytes) -> dict:
    return {"path": path, "size_bytes": len(data), "sha256": sha256_bytes(data)}


def _dataset_entry(spec: dict, fasta: bytes) -> tuple[dict, bytes]:
    sequence, records = _parse_fasta(fasta, spec["accessions"])
    base = f"data/{spec['id']}"
    entry = {
        "id": spec["id"],
        "label": spec["label"],
        "accessions": spec["accessions"],
        "fasta": _file_entry(base + ".fasta", fasta),
        "sequence": {
            **_file_entry(base + ".seq", sequence),
            "raw_acgt": all(byte in b"ACGT" for byte in sequence),
        },
        "records": records,
    }
    return entry, sequence


def _require_manifest(locked: dict, actual: dict, path: Path) -> None:
    # JSON comparison also distinguishes booleans from integers (True == 1).
    if json.dumps(locked, sort_keys=True) != json.dumps(actual, sort_keys=True):
        raise ValueError(f"Immutable {path.name} mismatch: size/hash/metadata changed")


def _write_once(path: Path, data: bytes) -> None:
    if path.exists():
        if path.stat().st_size != len(data) or sha256_file(path) != sha256_bytes(data):
            raise ValueError(f"Refusing to replace changed cached file: {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as stream:
        stream.write(data)


def _write_lock(path: Path, manifest: dict) -> None:
    if path.exists():
        _require_manifest(_load_json(path), manifest, path)
        return
    with path.open("xb") as stream:
        stream.write((json.dumps(manifest, indent=2, ensure_ascii=True) + "\n").encode("ascii"))


def verify_data(root: Path) -> dict:
    """Verify protocol, byte counts, hashes, derivation, and all record metadata."""
    root = Path(root)
    protocol, protocol_hash = _protocol(root)
    lock_path = root / "data.lock.json"
    locked = _load_json(lock_path)
    if locked.get("protocol_sha256") != protocol_hash:
        raise ValueError("Immutable data.lock.json protocol_sha256 mismatch")
    entries = []
    for spec in protocol["datasets"]:
        base = root / "data" / spec["id"]
        fasta = _read_bounded(base.with_suffix(".fasta"), FASTA_LIMIT * len(spec["accessions"]))
        entry, sequence = _dataset_entry(spec, fasta)
        cached_sequence = _read_bounded(base.with_suffix(".seq"), len(sequence))
        if cached_sequence != sequence:
            raise ValueError(f"Cached sequence differs from exact FASTA derivation: {spec['id']}")
        entries.append(entry)
    actual = {"schema_version": SCHEMA_VERSION, "protocol_sha256": protocol_hash, "datasets": entries}
    _require_manifest(locked, actual, lock_path)
    return locked


def prepare_data(root: Path) -> dict:
    """Fetch each unique accession once, or strictly verify an existing lock."""
    root = Path(root)
    if (root / "data.lock.json").exists():
        return verify_data(root)
    protocol, protocol_hash = _protocol(root)
    downloaded = {}
    for spec in protocol["datasets"]:
        for accession in spec["accessions"]:
            if accession not in downloaded:
                url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urlencode({
                    "db": "nuccore", "id": accession, "rettype": "fasta", "retmode": "text",
                })
                fasta = _download(url, FASTA_LIMIT)
                _parse_fasta(fasta, [accession])
                downloaded[accession] = fasta
    entries = []
    artifacts = []
    for spec in protocol["datasets"]:
        fasta = b"".join(downloaded[accession] for accession in spec["accessions"])
        entry, sequence = _dataset_entry(spec, fasta)
        entries.append(entry)
        artifacts.extend(((entry["fasta"]["path"], fasta), (entry["sequence"]["path"], sequence)))
    for path, data in artifacts:
        _write_once(root / path, data)
    manifest = {"schema_version": SCHEMA_VERSION, "protocol_sha256": protocol_hash, "datasets": entries}
    _write_lock(root / "data.lock.json", manifest)
    return verify_data(root)


def _source_entry(source_id: str, commit: str, archive: bytes) -> dict:
    return {
        "id": source_id,
        "commit": commit,
        "url": f"https://codeload.github.com/cobilab/{source_id}/tar.gz/{commit}",
        **_file_entry(f"vendor/{source_id}.tar.gz", archive),
    }


def _verify_sources(root: Path) -> dict:
    lock_path = root / "sources.lock.json"
    locked = _load_json(lock_path)
    entries = []
    for source_id, commit in SOURCES:
        archive = _read_bounded(root / "vendor" / f"{source_id}.tar.gz", SOURCE_LIMIT)
        entries.append(_source_entry(source_id, commit, archive))
    _require_manifest(locked, {"schema_version": SCHEMA_VERSION, "sources": entries}, lock_path)
    return locked


def prepare_sources(root: Path) -> dict:
    """Freeze commit-pinned archives; an existing lock permits no network repair."""
    root = Path(root)
    if (root / "sources.lock.json").exists():
        return _verify_sources(root)
    entries = []
    artifacts = []
    for source_id, commit in SOURCES:
        url = f"https://codeload.github.com/cobilab/{source_id}/tar.gz/{commit}"
        archive = _download(url, SOURCE_LIMIT)
        entries.append(_source_entry(source_id, commit, archive))
        artifacts.append((f"vendor/{source_id}.tar.gz", archive))
    for path, archive in artifacts:
        _write_once(root / path, archive)
    manifest = {"schema_version": SCHEMA_VERSION, "sources": entries}
    _write_lock(root / "sources.lock.json", manifest)
    return _verify_sources(root)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--sources-only", action="store_true")
    group.add_argument("--data-only", action="store_true")
    group.add_argument("--verify", action="store_true")
    args = parser.parse_args(argv)
    root = Path(__file__).resolve().parent
    try:
        if args.verify:
            data = verify_data(root)
            sources = _verify_sources(root)
            print(f"Verified {len(data['datasets'])} datasets and {len(sources['sources'])} source archives (offline)")
        else:
            if not args.sources_only:
                data = prepare_data(root)
                print(f"Prepared {len(data['datasets'])} datasets: data.lock.json", flush=True)
            if not args.data_only:
                sources = prepare_sources(root)
                print(f"Prepared {len(sources['sources'])} source archives: sources.lock.json", flush=True)
    except (OSError, ValueError) as error:
        parser.exit(1, f"error: {error}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
