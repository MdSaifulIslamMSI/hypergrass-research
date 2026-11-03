#!/usr/bin/env sh
# Build the static-musl native-image pilot.
#
# This is an opt-in experiment, not the default build. The default Java build
# remains `build.sh`/`build.ps1` with --release 27. GraalVM 25 is the newest
# native-image release available, so the pilot compiles the same sources with
# --release 25; the resulting archive bytes are verified against the JVM.
set -eu
cd "$(dirname "$0")"

IMAGE="ghcr.io/graalvm/native-image-community:25-muslib@sha256:bd08fee5f6effb210e829f434e5b5341d9b5b321ddc47363a241853723c810eb"
BUILD_MEMORY="${BUILD_MEMORY:-3500m}"

docker run --rm --platform linux/amd64 --memory="$BUILD_MEMORY" \
    -v "$PWD:/work" -w /work --entrypoint /bin/sh "$IMAGE" -lc '
set -eu
rm -rf bin/native-classes bin/native
mkdir -p bin/native-classes bin/native
javac --release 25 -Xlint:all -d bin/native-classes src/hypergrass/*.java
native-image \
    --static --libc=musl --no-fallback --gc=serial \
    -O3 \
    -march=compatibility \
    -J-Xmx2500m \
    -cp bin/native-classes \
    -o bin/native/hypergrass \
    hypergrass.HyperGraSS
sha256sum bin/native/hypergrass
'
