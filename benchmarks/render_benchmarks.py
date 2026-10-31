"""Render docs/BENCHMARKS.md tables from audited reports.

The interpretation section is intentionally left as a marker: prose claims are
reviewed by a human, not generated from the same script that computes them.

Usage:
    python -B benchmarks/render_benchmarks.py \
        --corpus benchmarks/corpus/runs/corpus-full-1/report.json \
        --scale benchmarks/pangenome/runs/pangenome-scale-2/report.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from summarize import CORPUS_ACCEPTED, corpus_coverage, corpus_tables, scale_tables  # noqa: E402


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--corpus", required=True)
    parser.add_argument("--scale")
    parser.add_argument("--doc", default="docs/BENCHMARKS.md")
    args = parser.parse_args(argv)

    doc = Path(args.doc)
    text = doc.read_text(encoding="utf-8")
    corpus = json.loads(Path(args.corpus).read_text(encoding="utf-8"))
    if corpus.get("status") not in CORPUS_ACCEPTED:
        raise SystemExit(f"corpus report is not finished: {corpus.get('status')}")

    infeasible = sum(
        1 for entry in corpus["results"]
        if entry["profile"] in ("geco3-l9", "jarvis3-l26")
        and entry["status"] not in ("success", "not_applicable")
    )
    present, planned, complete, total = corpus_coverage(corpus)
    if corpus.get("status") in ("success", "partial"):
        status = (
            f"`{corpus['status']}`; {infeasible} high-effort `geco3-l9` or "
            f"`jarvis3-l26` cells exceeded the cap or memory limit and are "
            f"marked infeasible. All {present:,} results cover {complete} of "
            f"{total} datasets; every other applicable cell passed round-trip "
            f"verification."
        )
    else:
        status = (
            f"`{corpus['status']}` — the run was stopped before the runner wrote "
            f"a final status. It contains {present:,} of {planned:,} planned "
            f"results and {complete} of {total} complete datasets, with "
            f"{infeasible} infeasible high-effort cells. This is a partial "
            f"matrix, not a finished benchmark."
        )
    text = text.replace("<!-- RUN-STATUS -->", status)
    text = text.replace("<!-- CORPUS-TABLES -->", "\n".join(corpus_tables(corpus)))
    if args.scale:
        scale = json.loads(Path(args.scale).read_text(encoding="utf-8"))
        if scale.get("status") not in ("success", "partial"):
            raise SystemExit(f"scale report is not finished: {scale.get('status')}")
        text = text.replace("<!-- SCALE-TABLES -->", "\n".join(scale_tables(scale)))
    doc.write_text(text, encoding="utf-8")
    print(f"rendered {doc}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
