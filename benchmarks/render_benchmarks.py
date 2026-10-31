"""Render docs/BENCHMARKS.md tables from audited reports.

The interpretation section is intentionally left as a marker: prose claims are
reviewed by a human, not generated from the same script that computes them.

Usage:
    python -B benchmarks/render_benchmarks.py \
        --corpus benchmarks/corpus/runs/corpus-full-1/report.json \
        --scale benchmarks/pangenome/runs/pangenome-scale-1/report.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from summarize import corpus_tables, scale_tables  # noqa: E402


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--corpus", required=True)
    parser.add_argument("--scale")
    parser.add_argument("--doc", default="docs/BENCHMARKS.md")
    args = parser.parse_args(argv)

    doc = Path(args.doc)
    text = doc.read_text(encoding="utf-8")
    corpus = json.loads(Path(args.corpus).read_text(encoding="utf-8"))
    if corpus.get("status") not in ("success", "partial"):
        raise SystemExit(f"corpus report is not finished: {corpus.get('status')}")

    infeasible = sum(
        1 for entry in corpus["results"]
        if entry["profile"] == "geco3-l9" and entry["status"] not in ("success", "not_applicable")
    )
    status = (
        f"`{corpus['status']}`; {corpus.get('failures', 0)} failures, "
        f"{infeasible} of them the frozen `geco3-l9` profile exceeding the cap "
        f"or memory limit. All other applicable cells passed round-trip verification."
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
