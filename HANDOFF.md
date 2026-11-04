# HyperGraSS / RC-MDL research handoff

Snapshot: **2026-11-04**
Repository: `https://github.com/MdSaifulIslamMSI/hypergrass-research` (private)
Release asset: `hypergrass-handoff-v1.tar.gz`

This package is the compact, byte-exact research state for the lossless single-genome
compressor project. It contains all source, docs, result JSON, manifests, and the canonical
raw inputs needed to reproduce the measurements. The ~29 GB of benchmark corpus and run
output is **not** included because it is regenerable from the manifests and scripts; the
evidence (result JSON, report JSON, summaries) is included.

## Canonical hashes

| File | Bytes | SHA-256 |
|---|---:|---|
| `benchmarks/specialized/runs/full-all-methods-v1/input.seq` (masked *E. coli*) | 120,108,223 | `51cb3ecb2e9f532833caa0b60345f37d0342474e9f4a239801486f9bc9219acb` |
| `benchmarks/pangenome/data/pangenome.fasta` (raw *E. coli* panel) | 121,825,884 | `f03f0dd10c541b99ddf63d18053e3a9c6009910a81700ba18f15e24dcb7d2285` |
| `prototype/rcmdl/work/mtb/panel.fasta` | 107,257,786 | `4a62d9cf8563b01fe926d5a779dc5c44cfb3c04fd6e5ea4d07939056f905618f` |
| `prototype/rcmdl/work/saureus/panel.fasta` | 72,060,360 | `9413289ee51e253b6adf6c389d04f20765ad1d5ac815b33d2552d9b7732adbd5` |
| `prototype/rcmdl/vendor/agc.zip` (AGC 3.2.5) | 802,120 | `4467b4183b37affe60b0444dee255761f085c060935c2b2e29d1b7c4cd683acd` |

## Results in one screen

| Experiment | Result | Verification |
|---|---|---|
| RC-MDL LOO, *E. coli* (24 genomes, 120,108,223 bases) | **1,780,761 B / 0.1186 bits/base** | 24/24 byte-exact |
| AGC 3.2.5 LOO, same panel (`-s 200000`) | 8,384,117 B / 0.5584 bits/base | 24/24 byte-exact |
| RC-MDL vs AGC, *E. coli* | **4.71x smaller** | median 4.09x; AGC wins only targets 13, 14 |
| RC-MDL LOO, *M. tuberculosis* (24 genomes, 105,745,015 bases) | **122,998 B / 0.009305 bits/base** | 24/24 byte-exact |
| AGC 3.2.5 LOO, *M. tuberculosis* | 1,390,893 B | 11.3x larger than RC-MDL |
| Self-contained single genome (*E. coli* 00) | HyperGraSS 1,136,083 B vs RC-MDL 1,192,494 B | HyperGraSS wins by 4.7% |
| Progressive collection (24 *E. coli*) | HyperGraSS 3,810,409 B vs RC-MDL Prim **3,972,975 B** | HyperGraSS wins by 4.3% |
| Panel-size ablation (target 00) | 1/2/4/8/16/23 genomes → 249,491 / 185,439 / 146,758 / 102,931 / 87,254 / 73,194 B | all byte-exact |
| *S. aureus* LOO | incomplete (cancelled at ~5/25) | partial evidence preserved |

The central claim is the reference-conditioned one: when a 23-genome panel is already
present, RC-MDL costs 0.12 bits/base on *E. coli* and 0.0093 bits/base on *M. tuberculosis*,
4.7x and 11.3x smaller than AGC. HyperGraSS still wins the self-contained single-genome and
collection settings. See `prototype/rcmdl/RESULTS.md` for the full tables, caveats, and the
lower-bound analysis of the collection gap.

## Restore from the release asset

```powershell
git clone https://github.com/MdSaifulIslamMSI/hypergrass-research.git
cd hypergrass-research
gh release download handoff-v1 --pattern 'hypergrass-handoff-v1.tar.gz'
tar -xzf hypergrass-handoff-v1.tar.gz
```

Verify every file against `MANIFEST.sha256` (the archive contains it):

```powershell
Get-Content MANIFEST.sha256 | ForEach-Object {
  $hash, $path = $_ -split '  ', 2
  $actual = (Get-FileHash $path -Algorithm SHA256).Hash.ToLower()
  if ($actual -ne $hash) { throw "mismatch: $path" }
}
```

Then rebuild the per-genome files the codec expects:

```powershell
python prototype/rcmdl/split_panel.py --fasta benchmarks/pangenome/data/pangenome.fasta --out prototype/rcmdl/work/genomes
python prototype/rcmdl/split_input.py --input benchmarks/specialized/runs/full-all-methods-v1/input.seq --manifest benchmarks/pangenome/manifest.json --out prototype/rcmdl/work/genomes_masked
python prototype/rcmdl/split_panel.py --fasta prototype/rcmdl/work/mtb/panel.fasta --out prototype/rcmdl/work/mtb/genomes
python prototype/rcmdl/split_panel.py --fasta prototype/rcmdl/work/saureus/panel.fasta --out prototype/rcmdl/work/saureus/genomes
```

Compile and run the canonical *E. coli* leave-one-out experiment:

```powershell
javac -d prototype/rcmdl/bin2 prototype/rcmdl/Rcmdl.java
python prototype/rcmdl/run.py --genomes prototype/rcmdl/work/genomes_masked --out prototype/rcmdl/work/runs/loo24-masked --heap 1200m
```

The harness writes `results.json` and `summary.json`; it verifies each decoded target with an
independent SHA-256 before reporting `verified`.

## Regenerating the bulk data

The ~29 GB under `benchmarks/**/data` and `benchmarks/**/runs` is reproducible:

- **E. coli panel:** `python benchmarks/pangenome/fetch.py`. `benchmarks/pangenome/manifest.json`
  locks the NCBI query, accession versions, lengths, and the FASTA SHA-256.
- **M. tuberculosis panel:** `python prototype/rcmdl/fetch_panel.py --query '"Mycobacterium tuberculosis"[Organism] AND "complete genome"[Title] AND refseq[Filter]' --out prototype/rcmdl/work/mtb --target-bases 100000000`
- **S. aureus panel:** same command with the *Staphylococcus aureus* query and
  `--out prototype/rcmdl/work/saureus --target-bases 70000000`.
- **AGC 3.2.5:** the exact Windows build is in `prototype/rcmdl/vendor/agc.zip`
  (SHA-256 above); extract `agc.exe` to `prototype/rcmdl/vendor/agc/agc.exe`.
- **HyperGraSS runs:** `python benchmarks/specialized/run.py` with the provenance recorded in
  each `report.json`. The canonical HyperGraSS collection baseline is
  `benchmarks/specialized/runs/hypergrass-lz-1/report.json`.

The result/report JSONs are the evidence. The bulk files are inputs and intermediates.

## Recommended way to continue (free first)

1. **Private GitHub repository** — done. Code, docs, result JSON, manifests, and restore
   scripts are versioned and available on any machine.
2. **GitHub Release asset** — done. The canonical raw inputs travel as one verified
   `tar.gz`, so no laptop copy is the only copy.
3. **Free, no credit card: GitHub Codespaces.** Open the repository, choose
   **Code → Codespaces → Create codespace on master**. The included
   `.devcontainer/` downloads and verifies the 122 MiB handoff, rebuilds the per-genome
   files, and compiles `Rcmdl.java` automatically. GitHub Free includes **120 core-hours per
   month**; a 2-core/8 GB machine is roughly 60 hours/month, a 4-core/8 GB machine roughly
   30 hours/month. The 32 GB disk is enough for the compact handoff and one experiment at a
   time. It is not intended to hold the full 29 GB benchmark corpus.
4. **Free persistent VM, if you can pass card verification: Oracle Cloud Always Free.**
   The Ampere A1 allowance is 4 ARM OCPUs, 24 GB RAM, 200 GB storage, and 10 TB/month
   egress at no charge within the always-free limits. A card is used for identity
   verification, not a charge. Build AGC from source on ARM (or run the x86 binary under
   emulation); RC-MDL and HyperGraSS are Java and portable. Capacity varies by region.
5. **Student?** The GitHub Student Developer Pack may add Codespaces hours and cloud
   credits if you are eligible.
6. **Paid only if the free quota is insufficient: Google Compute Engine.** For reference,
   `e2-standard-4` (4 vCPU, 16 GB) or `e2-custom-4-8192` (4 vCPU, 8 GB), 100 GB
   `pd-balanced`, `us-central1`, Spot. Roughly $0.04–0.06/hr Spot, $0.13–0.16/hr on-demand,
   plus about $10/month for the disk.
7. **Bulk storage only if exact run artifacts must be kept: Google Cloud Storage.** Standard
   is roughly $0.02/GB/month (29 GB ≈ $0.58/month), Nearline ≈ $0.01/GB/month, Coldline
   ≈ $0.004/GB/month. Do not pay to store regenerable panels; re-fetch them on the compute
   machine instead.
8. **Local fallback:** the laptop now has ~190 GB free and the remaining checkout is only
   ~11 MB. Fetch one panel at a time, keep the result JSON, and delete the bulk outputs
   afterwards.

**GCP note:** project `billy-b674c` has billing disabled and the Compute Engine API is not
enabled. That only blocks the paid GCE/GCS path; GitHub Codespaces needs no billing and is
already available.

## What was removed from the laptop

After the release asset was uploaded and its SHA-256 verified, the regenerable benchmark
data/runs (~29 GB) and the prototype work data (~1 GB) were deleted locally. The private
repository plus the release asset are the off-laptop source of truth. The exact deletion
list and freed space are in the session report.
