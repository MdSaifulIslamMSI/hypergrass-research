#!/usr/bin/env python3
"""Offline end-to-end smoke test for benchmarks/hprc_giab/run_hprc_giab_bench.py using native eimdl."""

from __future__ import annotations

import json
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile


def main() -> None:
    repo = Path(__file__).resolve().parents[2]
    bench_script = repo / "benchmarks" / "hprc_giab" / "run_hprc_giab_bench.py"

    env_bin = os.environ.get("EIMDL_BIN")
    candidates = ([Path(env_bin)] if env_bin else []) + [
        repo / "native" / "eimdl-rs" / "target" / "release" / "eimdl",
        repo / "native" / "eimdl-rs" / "target" / "release" / "eimdl.exe",
        repo / "native" / "eimdl-rs" / "target" / "debug" / "eimdl",
        repo / "native" / "eimdl-rs" / "target" / "debug" / "eimdl.exe",
    ]
    eimdl_bin = next((c for c in candidates if c.is_file()), None)
    assert eimdl_bin is not None, f"Native eimdl binary not built under {repo / 'native' / 'eimdl-rs' / 'target'}"

    dry = subprocess.run([sys.executable, str(bench_script), "--eimdl", str(eimdl_bin), "--dry-run"], check=True, capture_output=True, text=True)
    dry_json = json.loads(dry.stdout)
    assert "tools_detected" in dry_json

    rng = random.Random(2026)
    bases = "ACGT"
    ref_list = [bases[rng.randrange(4)] for _ in range(40_000)]
    ref_seq = "".join(ref_list)
    alt_list = list(ref_list)
    truth_rows = []
    for pos in range(500, 39_500, 2000):
        old_b = alt_list[pos]
        new_b = bases[(bases.index(old_b) + 1) & 3]
        alt_list[pos] = new_b
        truth_rows.append((pos + 1, old_b, new_b))
    alt_seq = "".join(alt_list)

    with tempfile.TemporaryDirectory(prefix="hprc_giab_smoke_") as td:
        td_p = Path(td)
        pdir = td_p / "panel"
        pdir.mkdir()
        (pdir / "chr21.fa").write_text(f">chr21\n{ref_seq}\n", encoding="utf-8")

        hprc_fa = td_p / "hprc_chr21_4haps.fa"
        hprc_fa.write_text(
            f">hap1\n{ref_seq}\n>hap2\n{alt_seq}\n>hap3\n{ref_seq}\n>hap4\n{alt_seq}\n",
            encoding="utf-8",
        )

        fq_path = td_p / "hg002_sim.fastq"
        qual = "I" * 150
        with fq_path.open("w", encoding="ascii", newline="\n") as f:
            for idx in range(5_000):  # ~18.75x read coverage
                s = rng.randrange(0, len(alt_seq) - 150)
                rd = list(alt_seq[s : s + 150])
                if idx % 23 == 0:
                    rd[70] = "T" if rd[70] != "T" else "A"
                f.write(f"@HG002.{idx}\n{''.join(rd)}\n+\n{qual}\n")

        truth_vcf = td_p / "HG002.truth.vcf"
        lines = ["##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
        for p1, rb, ab in truth_rows:
            lines.append(f"chr21\t{p1}\t.\t{rb}\t{ab}\t.\tPASS\t.")
        truth_vcf.write_text("\n".join(lines) + "\n", encoding="utf-8")

        truth_bed = td_p / "HG002.highconf.bed"
        truth_bed.write_text("21\t0\t40000\n", encoding="utf-8")

        out_json = td_p / "report.json"
        subprocess.run(
            [
                sys.executable,
                str(bench_script),
                "--eimdl",
                str(eimdl_bin),
                "--hprc-fasta",
                str(hprc_fa),
                "--panel",
                str(pdir),
                "--giab-fastq",
                str(fq_path),
                "--giab-truth-vcf",
                str(truth_vcf),
                "--giab-truth-bed",
                str(truth_bed),
                "--min-support",
                "2",
                "--threads",
                "2",
                "--out",
                str(out_json),
            ],
            check=True,
        )

        rep = json.loads(out_json.read_text(encoding="utf-8"))
        assert rep["hprc"]["runs"][0]["sha256_verified"] is True
        assert rep["giab_hg002"]["sha256_verified"] is True
        conc = rep["giab_hg002"]["giab_exact_concordance"]
        assert conc["precision"] >= 0.95, f"Precision too low: {conc}"
        assert conc["recall"] >= 0.85, f"Recall too low: {conc}"
        print(f"PASS: HPRC/GIAB smoke suite verified (F1={conc['f1']}, P={conc['precision']}, R={conc['recall']})")


if __name__ == "__main__":
    main()
