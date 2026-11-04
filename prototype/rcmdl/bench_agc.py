#!/usr/bin/env python3
"""AGC marginal-cost benchmark on the same 24-genome E. coli panel.

For every target i:
  base_i   = AGC archive of the other 23 genomes (adaptive mode)
  full_i   = base_i with target i appended
  marginal = |full_i| - |base_i|

That marginal is the closest AGC analogue of an RC-MDL leave-one-out archive:
the same 23-genome panel, only the target's bytes charged to the target.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import time
from pathlib import Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run(cmd: list[str]) -> float:
    start = time.perf_counter()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.perf_counter() - start
    if proc.returncode != 0:
        raise RuntimeError(
            "command failed (%d): %s\nSTDOUT:\n%s\nSTDERR:\n%s"
            % (proc.returncode, " ".join(str(c) for c in cmd), proc.stdout, proc.stderr)
        )
    return elapsed


def read_fasta_seq(path: Path) -> bytes:
    chunks = []
    with path.open("rb") as handle:
        for line in handle:
            if not line.startswith(b">"):
                chunks.append(line.strip())
    return b"".join(chunks)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta-dir", type=Path, default=Path("prototype/rcmdl/work/fasta"))
    parser.add_argument("--seq-dir", type=Path, default=Path("prototype/rcmdl/work/genomes"))
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--agc", type=Path, default=Path("prototype/rcmdl/vendor/agc/agc.exe"))
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--create-opts", default="-a", help="extra create options, space separated")
    parser.add_argument("--append-opts", default="-a", help="extra append options, space separated")
    parser.add_argument("--targets", default="", help="comma list of indices; default all")
    parser.add_argument("--full", action="store_true", help="also build one 24-genome archive")
    args = parser.parse_args()

    fastas = sorted(args.fasta_dir.glob("*.fa"))
    if not fastas:
        raise SystemExit("no fasta files")
    if args.targets:
        wanted = {int(x) for x in args.targets.split(",") if x.strip()}
        fastas = [p for i, p in enumerate(fastas) if i in wanted]
    args.out.mkdir(parents=True, exist_ok=True)
    create_opts = args.create_opts.split()
    append_opts = args.append_opts.split()
    records = []

    if args.full:
        archive = args.out / "full24.agc"
        ref, others = fastas[0], fastas[1:]
        listfile = args.out / "full24.list"
        listfile.write_text("\n".join(str(p.resolve()) for p in others) + "\n", encoding="utf-8")
        t = run(
            [
                str(args.agc),
                "create",
                *create_opts,
                "-t",
                str(args.threads),
                "-i",
                str(listfile),
                "-o",
                str(archive),
                str(ref.resolve()),
            ]
        )
        record = {
            "mode": "full",
            "genomes": len(fastas),
            "archive_bytes": archive.stat().st_size,
            "bases": sum(read_fasta_seq(p).__len__() for p in fastas),
            "seconds": t,
        }
        records.append(record)
        print(record, flush=True)
        (args.out / "results.json").write_text(
            json.dumps(records, indent=2) + "\n", encoding="utf-8"
        )

    for target in fastas:
        index = int(target.stem)
        others = [p for p in sorted(args.fasta_dir.glob("*.fa")) if p != target]
        ref, rest = others[0], others[1:]
        base = args.out / ("%02d-base.agc" % index)
        full = args.out / ("%02d-full.agc" % index)
        listfile = args.out / ("%02d-panel.list" % index)
        listfile.write_text("\n".join(str(p.resolve()) for p in rest) + "\n", encoding="utf-8")
        t_create = run(
            [
                str(args.agc),
                "create",
                *create_opts,
                "-t",
                str(args.threads),
                "-i",
                str(listfile),
                "-o",
                str(base),
                str(ref.resolve()),
            ]
        )
        t_append = run(
            [
                str(args.agc),
                "append",
                *append_opts,
                "-t",
                str(args.threads),
                "-o",
                str(full),
                str(base.resolve()),
                str(target.resolve()),
            ]
        )
        listset = subprocess.run(
            [str(args.agc), "listset", str(full.resolve())],
            capture_output=True,
            text=True,
            check=True,
        )
        names = [line.strip() for line in listset.stdout.splitlines() if line.strip()]
        sample = target.stem if target.stem in names else None
        if sample is None:
            for name in names:
                if name.split()[0] == target.stem:
                    sample = name
                    break
        verified = False
        restored_len = 0
        if sample is not None:
            restored = args.out / ("%02d-restored.fa" % index)
            with restored.open("wb") as handle:
                subprocess.run(
                    [str(args.agc), "getset", str(full.resolve()), sample],
                    stdout=handle,
                    check=True,
                )
            seq = read_fasta_seq(restored)
            wanted = read_fasta_seq(target)
            verified = seq == wanted
            restored_len = len(seq)
            restored.unlink()
        delta = full.stat().st_size - base.stat().st_size
        record = {
            "mode": "loo",
            "index": index,
            "target": target.name,
            "target_bytes": target.stat().st_size,
            "base_bytes": base.stat().st_size,
            "full_bytes": full.stat().st_size,
            "marginal_bytes": delta,
            "bpb": delta * 8.0 / target.stat().st_size,
            "create_seconds": t_create,
            "append_seconds": t_append,
            "sample": sample,
            "verified": verified,
            "restored_bytes": restored_len,
        }
        records.append(record)
        print(
            "%2d %-4s base %9d  full %9d  delta %8d  %7.4f bpb  verified=%s"
            % (
                index,
                target.stem,
                base.stat().st_size,
                full.stat().st_size,
                delta,
                record["bpb"],
                verified,
            ),
            flush=True,
        )
        (args.out / "results.json").write_text(
            json.dumps(records, indent=2) + "\n", encoding="utf-8"
        )
        base.unlink()
        full.unlink()


if __name__ == "__main__":
    main()
