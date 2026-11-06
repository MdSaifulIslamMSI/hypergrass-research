#!/usr/bin/env sh
set -eu
cd "$(dirname "$0")"
./build.sh

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
