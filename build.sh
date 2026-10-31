#!/usr/bin/env sh
set -eu
cd "$(dirname "$0")"
mkdir -p bin
javac --release 21 -Xlint:all -d bin src/hypergrass/*.java
printf 'Built %s classes into bin/\n' "$(find bin/hypergrass -name '*.class' | wc -l)"
