#!/usr/bin/env sh
set -eu
if [ "${CI:-}" != "true" ] || [ "${GITHUB_ACTIONS:-}" != "true" ]; then
  echo "MANDATORY POLICY: Never run tests locally." >&2
  echo "Run all verification suites and world-compressor tests end-to-end on Kaggle:" >&2
  echo "  python benchmarks/kaggle/run_kaggle_e2e.py" >&2
  exit 1
fi
cd "$(dirname "$0")"
if [ "${SKIP_BUILD:-0}" != "1" ]; then
  ./build.sh
fi

echo "== VerifyAllEdgeCases =="
java -cp bin hypergrass.VerifyAllEdgeCases

echo "== VerifyDonorIndex =="
java -cp bin hypergrass.VerifyDonorIndex

echo "== VerifyAnchorCodecs =="
java -Xmx512m -cp bin hypergrass.VerifyAnchorCodecs

echo "== VerifyArchiveAccounting =="
java -Xmx512m -cp bin hypergrass.VerifyArchiveAccounting

echo "== VerifyDonorCosts =="
java -Xmx512m -cp bin hypergrass.VerifyDonorCosts

echo "== VerifyCompressionLaws =="
java -Xmx512m -cp bin hypergrass.VerifyCompressionLaws

echo "== VerifyFusionEngine =="
java -Xmx512m -cp bin hypergrass.VerifyFusionEngine

echo "All verification suites passed."
