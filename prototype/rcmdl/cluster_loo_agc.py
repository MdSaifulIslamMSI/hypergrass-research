#!/usr/bin/env python3
"""AGC leave-cluster-out benchmark, matching cluster_loo.py's panel exclusion."""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from pathlib import Path

from bench_agc import read_fasta_seq, run


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta-dir", type=Path, required=True)
    parser.add_argument("--cluster-summary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--agc", type=Path, default=Path("prototype/rcmdl/vendor/agc/agc.exe"))
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--create-opts", default="-a")
    parser.add_argument("--append-opts", default="-a")
    args = parser.parse_args()

    fastas = sorted(args.fasta_dir.glob("*.fa"))
    clusters = json.loads(args.cluster_summary.read_text())["clusters"]
    cluster_of: dict[str, set[str]] = {}
    for members in clusters:
        stems = {Path(name).stem for name in members}
        for stem in stems:
            cluster_of[stem] = stems

    args.out.mkdir(parents=True, exist_ok=True)
    create_opts = args.create_opts.split()
    append_opts = args.append_opts.split()
    records = []

    for target in fastas:
        index = int(target.stem)
        excluded = cluster_of.get(target.stem, {target.stem})
        panel = [p for p in fastas if p.stem not in excluded]
        if not panel:
            raise SystemExit("empty panel for %s" % target)
        ref, rest = panel[0], panel[1:]
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
        delta = full.stat().st_size - base.stat().st_size
        base_bytes = base.stat().st_size
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
            "excluded": sorted(excluded - {target.stem}),
            "panel_genomes": len(panel),
            "base_bytes": base_bytes,
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
            "%2d %-4s exclude=%d panel=%2d %9d  %7.4f bpb  verified=%s"
            % (index, target.stem, len(record["excluded"]), len(panel), delta, record["bpb"], verified),
            flush=True,
        )

    summary = {
        "targets": len(records),
        "target_bases": sum(r["target_bytes"] for r in records),
        "marginal_bytes": sum(r["marginal_bytes"] for r in records),
        "verified": sum(1 for r in records if r["verified"]),
        "all_verified": all(r["verified"] for r in records),
    }
    summary["bpb"] = summary["marginal_bytes"] * 8.0 / summary["target_bases"]
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
