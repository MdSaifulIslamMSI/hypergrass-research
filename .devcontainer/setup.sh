#!/usr/bin/env bash
# Download, verify, and restore the compact handoff into a fresh Codespace.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

ASSET="/tmp/hypergrass-handoff-v1.tar.gz"
EXPECTED="a3842a892560a261df7d62704979479a4f9d2f24a043a8ad0f48370c8f50046f"

if [ ! -f "$ASSET" ]; then
  echo "Downloading the verified handoff release..."
  gh release download handoff-v1 \
    --pattern 'hypergrass-handoff-v1.tar.gz' \
    --dir /tmp \
    --clobber
fi

ACTUAL="$(sha256sum "$ASSET" | awk '{print $1}')"
if [ "$ACTUAL" != "$EXPECTED" ]; then
  echo "SHA-256 mismatch: expected $EXPECTED, got $ACTUAL" >&2
  exit 1
fi

EXTRACT="/tmp/hypergrass-handoff-v1"
rm -rf "$EXTRACT"
mkdir -p "$EXTRACT"
tar -xzf "$ASSET" -C "$EXTRACT" --strip-components=1

mkdir -p benchmarks/pangenome/data
mkdir -p benchmarks/specialized/runs/full-all-methods-v1
mkdir -p prototype/rcmdl/work/mtb
mkdir -p prototype/rcmdl/work/saureus

cp "$EXTRACT/benchmarks/pangenome/data/pangenome.fasta" benchmarks/pangenome/data/pangenome.fasta
cp "$EXTRACT/benchmarks/specialized/runs/full-all-methods-v1/input.seq" benchmarks/specialized/runs/full-all-methods-v1/input.seq
cp "$EXTRACT/prototype/rcmdl/work/mtb/panel.fasta" prototype/rcmdl/work/mtb/panel.fasta
cp "$EXTRACT/prototype/rcmdl/work/saureus/panel.fasta" prototype/rcmdl/work/saureus/panel.fasta

python prototype/rcmdl/split_panel.py \
  --fasta benchmarks/pangenome/data/pangenome.fasta \
  --out prototype/rcmdl/work/genomes
python prototype/rcmdl/split_input.py \
  --input benchmarks/specialized/runs/full-all-methods-v1/input.seq \
  --manifest benchmarks/pangenome/manifest.json \
  --out prototype/rcmdl/work/genomes_masked
python prototype/rcmdl/split_panel.py \
  --fasta prototype/rcmdl/work/mtb/panel.fasta \
  --out prototype/rcmdl/work/mtb/genomes
python prototype/rcmdl/split_panel.py \
  --fasta prototype/rcmdl/work/saureus/panel.fasta \
  --out prototype/rcmdl/work/saureus/genomes

mkdir -p prototype/rcmdl/bin
javac -d prototype/rcmdl/bin prototype/rcmdl/Rcmdl.java

cat <<'EOF'
Ready.

Canonical E. coli leave-one-out run:
  python prototype/rcmdl/run.py \
    --genomes prototype/rcmdl/work/genomes_masked \
    --out prototype/rcmdl/work/runs/loo24-masked \
    --heap 1200m

Results: prototype/rcmdl/work/runs/loo24-masked/results.json
Guide:    HANDOFF.md and prototype/rcmdl/RESULTS.md
EOF
