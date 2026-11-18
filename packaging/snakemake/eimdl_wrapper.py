"""Snakemake wrapper for HyperGraSS / EIMDL v3 (`eimdl` Native CLI)."""

from __future__ import annotations

import subprocess
from pathlib import Path


def main() -> None:
    subcmd = str(snakemake.params.get("subcommand", "compress-omics"))
    threads = int(getattr(snakemake, "threads", 1) or 1)
    panel = snakemake.input.get("panel", None)
    extra = str(snakemake.params.get("extra", "")).split()

    cmd: list[str] = ["eimdl", subcmd]

    if subcmd in ("compress", "compress-omics"):
        cmd += ["--target", str(snakemake.input.target), "--out", str(snakemake.output.archive)]
        if panel:
            cmd += ["--panel", str(panel)]
        elif subcmd == "compress":
            cmd.append("--collection")
        if subcmd == "compress-omics" and "--cds" not in extra:
            cmd.append("--cds")
        cmd += ["--threads", str(threads)]
    elif subcmd in ("decompress", "decompress-omics"):
        cmd += ["--archive", str(snakemake.input.archive), "--out", str(snakemake.output.sequence)]
        if panel:
            cmd += ["--panel", str(panel)]
        cmd += ["--threads", str(threads)]
    elif subcmd == "variants":
        cmd += ["--archive", str(snakemake.input.archive), "--panel", str(panel), "--out", str(snakemake.output.vcf)]
    elif subcmd == "analyze":
        cmd += ["--archive", str(snakemake.input.archive), "--json"]
        if panel:
            cmd += ["--panel", str(panel)]
        cmd += ["--threads", str(threads)]
    elif subcmd == "kmer":
        query = str(snakemake.params["query"])
        cmd += ["--archive", str(snakemake.input.archive), "--query", query, "--both-strands", "--json"]
        if panel:
            cmd += ["--panel", str(panel)]
    else:
        raise ValueError(f"Unsupported EIMDL subcommand: {subcmd}")

    cmd += extra
    if subcmd in ("analyze", "kmer") and hasattr(snakemake.output, "json"):
        out_txt = subprocess.check_output(cmd, text=True)
        Path(str(snakemake.output.json)).write_text(out_txt, encoding="utf-8")
    else:
        subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
