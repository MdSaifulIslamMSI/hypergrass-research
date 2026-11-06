#!/usr/bin/env python3
"""Deterministic generator for the 8-corpus Non-Existing / Synthetic / Xenobiological
& Adversarial Genome Suite used to evaluate HyperGraSS + RC-MDL (EIMDL / HGS2)
beyond natural RefSeq organisms.

Every generated file is 100% deterministic (seeded via SHAKE-256 / Random(20261106))
and spans biological, mathematical, synthetic-biology, and cryptographic regimes:
  1. 01_euk_centromere.fasta     : Eukaryotic HOR satellite + soft-masked transposons + exons
  2. 02_viral_rna_acgu.fasta     : Positive-sense ACGU RNA viral quasispecies + stem-loop hairpins
  3. 03_chimeric_iupac_crlf.fa   : Synthetic construct mixing T+U, IUPAC RYSWKMBDHVN, CRLF, mixed widths
  4. 04_xeno_at88.fasta          : Xenobiological 88% AT endosymbiont chromosome with 3-periodic ORFs
  5. 05_xeno_gc86.fasta          : Xenobiological 86% GC extremophile chromosome with CpG + G4 motifs
  6. 06_denovo_grammar.fasta     : De novo synthetic regulatory/CRISPR/operon/inverted-repeat chromosome
  7. 07_hypermutator_panel/      : 6-strain hypermutator pangenome (2.5% Kimura+Wobble SNPs + micro-indels + inversions)
  8. 08_crypto_white_noise.seq   : Cryptographic SHAKE-256 white-noise DNA (exact 2.0000 bpb entropy ceiling test)
"""

from __future__ import annotations

import hashlib
import json
import random
from pathlib import Path

BASES = "ACGT"
RNA_BASES = "ACGU"
RC_MAP = str.maketrans("ACGTacgt", "TGCAtgca")
TRANSITION = {"A": "G", "G": "A", "C": "T", "T": "C"}
TRANSVERSIONS = {
    "A": ("C", "T"),
    "G": ("C", "T"),
    "C": ("A", "G"),
    "T": ("A", "G"),
}


def revcomp(seq: str) -> str:
    return seq.translate(RC_MAP)[::-1]


def wrap_fasta(header: str, seq: str, width: int = 80, line_sep: str = "\n") -> str:
    lines = [f">{header}"]
    for i in range(0, len(seq), width):
        lines.append(seq[i : i + width])
    return line_sep.join(lines) + line_sep


def mutate_biological(seq: str, snp_rate: float, indel_rate: float, rng: random.Random) -> str:
    out: list[str] = []
    i = 0
    n = len(seq)
    while i < n:
        if indel_rate > 0.0 and rng.random() < indel_rate:
            if rng.random() < 0.55 and i + 60 < n:
                # Micro-deletion (1..45 bp, enriched for 3n in-frame deletions)
                dlen = rng.choice([1, 2, 3, 3, 6, 6, 9, 12, 15, 21, 33, 45])
                i += dlen
                continue
            else:
                # Micro-insertion (1..24 bp)
                ilen = rng.choice([1, 2, 3, 3, 6, 9, 12, 18, 24])
                for _ in range(ilen):
                    out.append(BASES[rng.randrange(4)])
        b = seq[i]
        b_up = b.upper()
        # 3-periodic codon wobble enrichment (position i % 3 == 2 mutates 3.2x more often)
        local_rate = snp_rate * (2.1 if (i % 3 == 2) else 0.45)
        if b_up in BASES and rng.random() < local_rate:
            if rng.random() < 0.76:
                nb = TRANSITION[b_up]
            else:
                nb = TRANSVERSIONS[b_up][rng.randrange(2)]
            out.append(nb.lower() if b.islower() else nb)
        else:
            out.append(b)
        i += 1
    return "".join(out)


def generate_all(out_dir: Path) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}

    # 1. Eukaryotic centromeric HOR satellite + soft-masked transposons + coding exons (~600 kbp)
    rng1 = random.Random(2026110601)
    monomer27 = "".join(BASES[rng1.randrange(4)] for _ in range(27))
    hor162 = "".join(mutate_biological(monomer27, 0.04, 0.0, rng1) for _ in range(6))
    transposon = "".join(BASES[rng1.randrange(4)] for _ in range(4200))
    exon_codons = ["ATG", "GCT", "GGT", "CTG", "TTC", "AAG", "GAC", "GAG", "CAG", "ATC", "GTG", "TAA"]
    contigs1: list[str] = []
    for c in range(3):
        parts: list[str] = []
        # Soft-masked pericentromeric HOR array
        for _ in range(450):
            parts.append(mutate_biological(hor162, 0.012, 0.0005, rng1).lower())
        # Interspersed autonomous + RC transposons + coding gene islands
        for block in range(18):
            if block % 2 == 0:
                tp = transposon if block % 4 == 0 else revcomp(transposon)
                parts.append(mutate_biological(tp, 0.015, 0.001, rng1).lower())
            gene = "".join(exon_codons[rng1.randrange(len(exon_codons))] for _ in range(950))
            parts.append(mutate_biological(gene, 0.08, 0.0, rng1))
        contigs1.append(wrap_fasta(f"syn_euk_chr{c+1} pericentromere_and_exons", "".join(parts), 80))
    p1 = out_dir / "01_euk_centromere.fasta"
    p1.write_bytes("".join(contigs1).encode("ascii"))
    paths["01_euk_centromere"] = p1

    # 2. Viral RNA (ACGU) positive-sense quasispecies with stem-loop hairpins (~300 kbp)
    rng2 = random.Random(2026110602)
    master_dna_parts: list[str] = []
    for _ in range(120):
        stem = "".join(BASES[rng2.randrange(4)] for _ in range(28))
        loop = "AAAC"
        hairpin = stem + loop + revcomp(stem)
        orf_chunk = "".join(BASES[rng2.randrange(4)] for _ in range(190))
        master_dna_parts.append(hairpin + orf_chunk)
    master_dna = "".join(master_dna_parts)  # ~30 kbp viral genome
    rna_records: list[str] = []
    for s in range(10):
        mut = master_dna if s == 0 else mutate_biological(master_dna, 0.008, 0.0004, rng2)
        rna = mut.replace("T", "U")
        rna_records.append(wrap_fasta(f"syn_ribo_virus_isolate_{s:02d} segment=complete_rna", rna, 70))
    p2 = out_dir / "02_viral_rna_acgu.fasta"
    p2.write_bytes("".join(rna_records).encode("ascii"))
    paths["02_viral_rna_acgu"] = p2

    # 3. Synthetic chimeric T+U + full IUPAC + variable width + CRLF construct (~250 kbp)
    rng3 = random.Random(2026110603)
    iupac = "RYSWKMBDHVNryswkmbdhvn"
    cassette = "".join(BASES[rng3.randrange(4)] for _ in range(5000))
    rec3_bytes = bytearray()
    for r in range(6):
        seq_chars: list[str] = []
        base_seq = cassette + revcomp(cassette) + "".join(BASES[rng3.randrange(4)] for _ in range(31000))
        for idx, ch in enumerate(base_seq):
            if idx % 499 == 17:
                seq_chars.append("N" * 25)
            elif idx % 311 == 5:
                seq_chars.append(iupac[rng3.randrange(len(iupac))])
            elif idx % 97 == 0 and ch == "T":
                seq_chars.append("U")  # synthetic dU/rU incorporation!
            elif (idx // 120) % 7 == 0:
                seq_chars.append(ch.lower())
            else:
                seq_chars.append(ch)
        s_full = "".join(seq_chars)
        width = [60, 70, 80][r % 3]
        sep = "\r\n" if (r % 2 == 0) else "\n"
        rec3_bytes.extend(wrap_fasta(f"syn_chimera_construct_{r} iupac_mixed=true", s_full, width, sep).encode("ascii"))
    p3 = out_dir / "03_chimeric_iupac_crlf.fa"
    p3.write_bytes(bytes(rec3_bytes))
    paths["03_chimeric_iupac_crlf"] = p3

    # 4. Xenobiological Extreme 88% AT Genome (~500 kbp)
    rng4 = random.Random(2026110604)
    at_weights_by_phase = [
        ("A", 0.43, "T", 0.43, "G", 0.08, "C", 0.06),
        ("A", 0.45, "T", 0.41, "G", 0.07, "C", 0.07),
        ("A", 0.48, "T", 0.46, "G", 0.03, "C", 0.03),  # 94% AT at wobble 3rd codon site!
    ]
    at_seq: list[str] = []
    for i in range(500_000):
        ph = at_weights_by_phase[i % 3]
        u = rng4.random()
        if u < ph[1]:
            at_seq.append(ph[0])
        elif u < ph[1] + ph[3]:
            at_seq.append(ph[2])
        elif u < ph[1] + ph[3] + ph[5]:
            at_seq.append(ph[4])
        else:
            at_seq.append(ph[6])
    p4 = out_dir / "04_xeno_at88.fasta"
    p4.write_bytes(wrap_fasta("xeno_endosymbiont_AT88 complete_chromosome", "".join(at_seq), 80).encode("ascii"))
    paths["04_xeno_at88"] = p4

    # 5. Xenobiological Extreme 86% GC Genome (~500 kbp)
    rng5 = random.Random(2026110605)
    gc_seq: list[str] = []
    g4_motif = "GGGAGGGTGGGAGGG"
    while len(gc_seq) < 500_000:
        if rng5.random() < 0.004:
            gc_seq.extend(g4_motif)
            continue
        ph = len(gc_seq) % 3
        u = rng5.random()
        gc_prob = 0.93 if ph == 2 else 0.825
        if u < gc_prob * 0.51:
            gc_seq.append("G")
        elif u < gc_prob:
            gc_seq.append("C")
        elif u < gc_prob + (1.0 - gc_prob) * 0.5:
            gc_seq.append("A")
        else:
            gc_seq.append("T")
    p5 = out_dir / "05_xeno_gc86.fasta"
    p5.write_bytes(
        wrap_fasta("xeno_halophile_GC86 complete_chromosome", "".join(gc_seq[:500_000]), 80).encode("ascii")
    )
    paths["05_xeno_gc86"] = p5

    # 6. De Novo Regulatory / Operon / Palindromic / CRISPR Grammar Chromosome (~500 kbp)
    rng6 = random.Random(2026110606)
    promoters = [
        "TTGACANNNNNNNNNNNNNNNNNTATAAT".replace("N", "A"),
        "TTGACAGCGCGCGCGCGCGCGCGTATAAT",
        "AGGAGGTAAACCCGGGTTTAAACCC",
    ]
    crispr_repeat = "GTTTTAGAGCTATGCTGTTTTGAATGGTCCCAAAAC"
    gene_bank = ["".join(BASES[rng6.randrange(4)] for _ in range(2400)) for _ in range(28)]
    denovo_parts: list[str] = []
    while sum(len(x) for x in denovo_parts) < 500_000:
        # Operon with promoter + 3 paralogous genes (forward or inverted) + rho-independent terminator hairpin
        prom = promoters[rng6.randrange(len(promoters))]
        denovo_parts.append(prom)
        for _ in range(3):
            g = gene_bank[rng6.randrange(len(gene_bank))]
            if rng6.random() < 0.4:
                g = revcomp(g)
            denovo_parts.append(mutate_biological(g, 0.025, 0.0008, rng6))
        term_stem = "".join(BASES[rng6.randrange(4)] for _ in range(18))
        denovo_parts.append(term_stem + "TTAA" + revcomp(term_stem) + "TTTTTTTT")
        # Occasional CRISPR array (direct repeat + variable 32-bp spacer)
        if rng6.random() < 0.25:
            for _ in range(14):
                spacer = "".join(BASES[rng6.randrange(4)] for _ in range(32))
                denovo_parts.append(crispr_repeat + spacer)
    denovo_str = "".join(denovo_parts)[:500_000]
    p6 = out_dir / "06_denovo_grammar.fasta"
    p6.write_bytes(wrap_fasta("denovo_synthetic_operon_crispr_chr1", denovo_str, 80).encode("ascii"))
    paths["06_denovo_grammar"] = p6

    # 7. Hypermutator Indel + Inversion Pangenome (6 strains x ~250 kbp)
    rng7 = random.Random(2026110607)
    anc7 = "".join(BASES[rng7.randrange(4)] for _ in range(250_000))
    panel7_dir = out_dir / "07_hypermutator_panel"
    panel7_dir.mkdir(parents=True, exist_ok=True)
    p7_concat: list[str] = []
    for s in range(6):
        cur = anc7
        if s > 0:
            # Large 35-kbp chromosomal inversion on odd strains + 2.5% Kimura+Wobble SNPs + micro-indels
            if s % 2 == 1:
                cur = cur[:80_000] + revcomp(cur[80_000:115_000]) + cur[115_000:]
            cur = mutate_biological(cur, 0.025, 0.0012, rng7)
        (panel7_dir / f"strain_{s:02d}.seq").write_bytes(cur.encode("ascii"))
        p7_concat.append(wrap_fasta(f"hypermutator_strain_{s:02d}", cur, 80))
    p7_fa = out_dir / "07_hypermutator_pangenome.fasta"
    p7_fa.write_bytes("".join(p7_concat).encode("ascii"))
    paths["07_hypermutator_panel_dir"] = panel7_dir
    paths["07_hypermutator_pangenome"] = p7_fa

    # 8. Cryptographic White-Noise DNA Horizon (500,000 bases via SHAKE-256 -> exact 2.000000 bpb)
    raw_bytes = hashlib.shake_256(b"HyperGraSS-EIMDL-Law-IV-Horizon-2026").digest(125_000)
    crypto_chars: list[str] = []
    for byte_val in raw_bytes:
        crypto_chars.append(BASES[byte_val & 3])
        crypto_chars.append(BASES[(byte_val >> 2) & 3])
        crypto_chars.append(BASES[(byte_val >> 4) & 3])
        crypto_chars.append(BASES[(byte_val >> 6) & 3])
    p8 = out_dir / "08_crypto_white_noise.seq"
    p8.write_bytes("".join(crypto_chars).encode("ascii"))
    paths["08_crypto_white_noise"] = p8

    manifest = {
        k: {
            "path": str(v),
            "bytes": v.stat().st_size if v.is_file() else sum(f.stat().st_size for f in v.glob("*.seq")),
            "sha256": hashlib.sha256(v.read_bytes()).hexdigest() if v.is_file() else "directory",
        }
        for k, v in paths.items()
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return paths


if __name__ == "__main__":
    out = Path("benchmarks/fusion/synthetic_data")
    generated = generate_all(out)
    print(f"Generated {len(generated)} deterministic synthetic/non-existing corpora in {out}")
