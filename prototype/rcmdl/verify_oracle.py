#!/usr/bin/env python3
"""Verify the AGC per-target oracle: re-run each target's best swept configuration
and decode the appended sample byte-for-byte.

The sweep records marginal sizes but deletes its archives. This harness rebuilds
only the winning configuration for each target, checks that the marginal reproduces,
and verifies the decoded sequence.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from pathlib import Path

from sweep_agc import CONFIGS


def run(cmd: list[str]) -> float:
    start = time.perf_counter()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.perf_counter() - start
    if proc.returncode != 0:
        raise RuntimeError(
            "command failed (%d): %s\n%s\n%s"
            % (proc.returncode, " ".join(map(str, cmd)), proc.stdout, proc.stderr)
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
    parser.add_argument("--fasta-dir", type=Path, required=True)
    parser.add_argument("--sweep-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--agc", type=Path, default=Path("prototype/rcmdl/vendor/agc/agc.exe"))
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument(
        "--config",
        default="",
        help="verify one fixed sweep config for every target instead of each target's best",
    )
    args = parser.parse_args()

    configs = {name: (opts, append_opts) for name, opts, append_opts in CONFIGS}
    fastas = sorted(args.fasta_dir.glob("*.fa"))
    args.out.mkdir(parents=True, exist_ok=True)
    records = []

    for target in fastas:
        index = int(target.stem)
        sweep = json.loads((args.sweep_dir / ("results-%02d.json" % index)).read_text())
        if args.config:
            matches = [row for row in sweep if row["config"] == args.config]
            if not matches:
                raise RuntimeError("config %s missing for %s" % (args.config, target.name))
            best = matches[0]
        else:
            best = min(sweep, key=lambda row: row["marginal_bytes"])
        opts, append_opts = configs[best["config"]]
        others = [p for p in fastas if p != target]
        ref, rest = others[0], others[1:]
        base = args.out / ("%02d-base.agc" % index)
        full = args.out / ("%02d-full.agc" % index)
        listfile = args.out / ("%02d-panel.list" % index)
        listfile.write_text("\n".join(str(p.resolve()) for p in rest) + "\n", encoding="utf-8")
        t_create = run(
            [
                str(args.agc),
                "create",
                *opts,
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
        delta = full.stat().st_size - base.stat().st_size
        if delta != best["marginal_bytes"]:
            raise RuntimeError(
                "marginal did not reproduce for %s: %d != %d"
                % (target.name, delta, best["marginal_bytes"])
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
        if sample is None:
            raise RuntimeError("target sample not found in %s" % full)
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
        restored.unlink()
        base.unlink()
        full.unlink()
        listfile.unlink()
        record = {
            "index": index,
            "target": target.name,
            "target_bytes": target.stat().st_size,
            "config": best["config"],
            "marginal_bytes": delta,
            "bpb": delta * 8.0 / target.stat().st_size,
            "create_seconds": t_create,
            "append_seconds": t_append,
            "verified": verified,
        }
        records.append(record)
        (args.out / "results.json").write_text(
            json.dumps(records, indent=2) + "\n", encoding="utf-8"
        )
        print(
            "%2d %-4s %-24s %8d  %7.4f bpb  verified=%s"
            % (index, target.stem, best["config"], delta, record["bpb"], verified),
            flush=True,
        )

    target_bases = sum(r["target_bytes"] for r in records)
    marginal_bytes = sum(r["marginal_bytes"] for r in records)
    summary = {
        "targets": len(records),
        "target_bases": target_bases,
        "marginal_bytes": marginal_bytes,
        "bpb": marginal_bytes * 8.0 / target_bases,
        "verified": sum(1 for r in records if r["verified"]),
        "all_verified": all(r["verified"] for r in records),
    }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
