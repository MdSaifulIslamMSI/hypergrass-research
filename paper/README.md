# Research Monograph Build & Automated Numerical Audit (`paper/`)

This directory contains the reproducible publication pipeline for the **22-page two-column IEEE/ACM research monograph** on the three-generation **`HyperGraSS` / `RC-MDL` / `HGS2` (`EIMDL`)** genomic compression architecture.

## Contents

| File | Description |
|---|---|
| [`HGS2_EIMDL_Research_Paper.pdf`](HGS2_EIMDL_Research_Paper.pdf) | Compiled 22-page two-column PDF monograph (`0` `Overfull \hbox` warnings) |
| [`HGS2_EIMDL_Research_Paper.tex`](HGS2_EIMDL_Research_Paper.tex) | Generated LaTeX source file (built deterministically from committed JSON evidence) |
| [`references.bib`](references.bib) | BibTeX bibliography covering genomic compression, syncmers, Kimura/codon models, and MDL |
| [`generate_figures.py`](generate_figures.py) | Generates all 7 publication-grade 300 DPI figures (`figures/fig1_regime_overview.png` .. `fig7_gen1_gen2_foundations.png`) directly from committed JSON results |
| [`table_builders.py`](table_builders.py) | Builds the per-genome Mode R LOO, Cluster-LOO, Mode C, Mode S, and independent benchmark LaTeX tables directly from JSON |
| [`build_latex_paper.py`](build_latex_paper.py) | Generates `HGS2_EIMDL_Research_Paper.tex` and `references.bib` |
| [`build_markdown_paper.py`](build_markdown_paper.py) | Generates the full Markdown monograph [`../docs/PAPER-DRAFT.md`](../docs/PAPER-DRAFT.md) |
| [`audit_paper_numbers.py`](audit_paper_numbers.py) | Executes **1,274 automated numerical assertions** verifying that every byte count, base count, and bit-rate across `.tex`, `docs/PAPER-DRAFT.md`, `README.md`, and `HANDOFF.md` matches the committed JSON benchmark reports |
| [`figures/`](figures/) | `fig1_regime_overview.png` through `fig7_gen1_gen2_foundations.png` (PNG + PDF vector outputs) |

## Rebuilding Everything End-to-End

From the repository root (`C:\Users\mdsai\.gemini\antigravity\scratch\hypergrass`):

```powershell
# 1. Regenerate all 7 publication figures (PNG + PDF) from the committed JSON reports
python paper/generate_figures.py

# 2. Rebuild the LaTeX monograph (paper/HGS2_EIMDL_Research_Paper.tex) and Markdown monograph (docs/PAPER-DRAFT.md)
python paper/build_latex_paper.py
python paper/build_markdown_paper.py

# 3. Compile the 22-page PDF via pdflatex + bibtex
Push-Location paper
pdflatex -interaction=nonstopmode HGS2_EIMDL_Research_Paper.tex
bibtex HGS2_EIMDL_Research_Paper
pdflatex -interaction=nonstopmode HGS2_EIMDL_Research_Paper.tex
pdflatex -interaction=nonstopmode HGS2_EIMDL_Research_Paper.tex
Pop-Location

# 4. Run the automated numerical audit (verifies 100% of numbers against JSON benchmark reports)
python paper/audit_paper_numbers.py
```
