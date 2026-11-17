#!/usr/bin/env sh
set -eu
cd "$(dirname "$0")"
JAVA_RELEASE="${JAVA_RELEASE:-25}"
rm -rf bin
mkdir -p bin
javac --release "$JAVA_RELEASE" -Xlint:all -d bin src/hypergrass/*.java
printf 'Built %s classes (release %s) into bin/\n' "$(find bin/hypergrass -name '*.class' | wc -l)" "$JAVA_RELEASE"
