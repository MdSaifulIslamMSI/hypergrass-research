$ErrorActionPreference = 'Stop'
# Opt-in static-musl native-image pilot (GraalVM 25 LTS / --release 25).
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $root
$image = 'ghcr.io/graalvm/native-image-community:25-muslib@sha256:bd08fee5f6effb210e829f434e5b5341d9b5b321ddc47363a241853723c810eb'
$buildMemory = if ($env:BUILD_MEMORY) { $env:BUILD_MEMORY } else { '3500m' }
$script = @'
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
'@
docker run --rm --platform linux/amd64 --memory=$buildMemory `
    -v "${root}:/work" -w /work --entrypoint /bin/sh $image -lc $script
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
