# Corpus image build

Run from the HyperGraSS project root with Python 3.11+ and a Docker daemon
supporting Linux/amd64 images and BuildKit. Do not start the build until the
source-preparation step has finished and a build has been authorized.

```text
python benchmarks/corpus/build_image.py
```

No source downloads are performed. The script requires `sources.lock.json`
with a `sources` list containing exactly these two entries:

| ID | Commit | Local archive |
| --- | --- | --- |
| `jarvis3` | `4dcb37b05b1a1482c2bf864d4968244817b0ad0a` | `vendor/jarvis3.tar.gz` |
| `geco3` | `69642a0a913e543599f35aee372e300601e4cbbe` | `vendor/geco3.tar.gz` |

Each entry must contain `id`, `commit`, `url` (HTTPS provenance), `path`,
`size_bytes`, and `sha256`. Paths are relative to this directory. Both size
and SHA-256 are checked before invoking Docker and again inside the image.
The archives are the upstream snapshots with one top-level directory;
both pinned Makefiles were inspected and place their binaries in `src/`.

## Image and artifacts

- Base: `linux/amd64`,
  `eclipse-temurin@sha256:56f10ddc224d7b900698422ffe82432ed0e178f606bd66859ed1d179825dfbe3`
  (Alpine 3.24).
- Packages: Python 3, `build-base` (including GCC, make, and musl development
  files), GNU gzip, xz, and bash. Package installation is the only networked
  `RUN` step; compilation and environment recording have networking disabled.
- Both tools use `make -C src CC=gcc` with command-line `CFLAGS` overrides,
  including recursive make invocations. JARVIS3 uses `-Wall -Ofast -DPROGRESS`;
  GeCo3/GeDe3 use `-Wall -O3`. Upstream native CPU flags are not used.
- Installed tools: `/usr/local/bin/JARVIS3`, `/usr/local/bin/GeCo3`, and
  `/usr/local/bin/GeDe3`.
- Complete upstream source trees and licenses: `/opt/sources/jarvis3/` and
  `/opt/sources/geco3/`. The source lock, original archives, and build helper
  remain under `/opt/benchmark/`.
- Image environment record: `/opt/benchmark-environment.json`. It contains
  source provenance, portable flags, executable SHA-256 hashes, GCC/Java/
  Python/gzip/xz version output, Alpine release, and installed APK versions.
- Host output: `benchmarks/corpus/build.log` (appended per build) and
  `benchmarks/corpus/image.lock.json` (immutable image ID, build-input hashes,
  and the complete environment record).
- Accepted image tag: `hypergrass-corpus:locked`.

The corpus directory is the Docker build context, not the project root.
`.dockerignore` permits only the Dockerfile, ignore file, build helper,
source lock, and the two archives. Corpus inputs, results, Java code, and
other workspace files are not sent to the builder.

The script builds by immutable image ID, checks the Linux/amd64 image, and
extracts metadata with a short-lived `--rm --network=none --pull=never`
container with **no mounts**. It tags only an accepted image. An existing
lock must match the build inputs, image ID, and environment exactly; an
incompatible lock is an error, not an implicit update. To intentionally
replace one, first archive the previous lock and build log, then explicitly
remove the lock before rebuilding. No containers or images are pruned.

The base and source revisions are pinned, but the Alpine package repository
is not a package snapshot. A fresh build can resolve different package
versions; those versions are recorded, and differing rebuilds are rejected
rather than replacing a lock. Retain/export the locked image if exact
long-term availability is required.

## Runtime boundary

The image uses `/work` as its working directory. Benchmark runtime containers
must have networking disabled and bind **only this corpus directory** to
`/work`, never the whole HyperGraSS workspace. For example, replace
`ABSOLUTE_CORPUS_PATH` with the absolute path to `benchmarks/corpus`:

```text
docker run --rm --pull=never --network=none --platform linux/amd64 --mount "type=bind,source=ABSOLUTE_CORPUS_PATH,target=/work" hypergrass-corpus:locked python3 --version
```

Prefer the `image_id` from `image.lock.json` over the mutable tag for measured
runs. This tooling does not stage Java classes or change benchmark runners.

## Host validation

```text
python -B -m unittest discover -s benchmarks/corpus -p test_build_image.py -v
python -B benchmarks/corpus/build_image.py --verify-sources benchmarks/corpus
```

The tests use temporary directories and mock every Docker invocation; they
never build images, start containers, or download sources. The second command
only validates local source locks. Actual image build and compressor runtime
smoke tests require separate authorization; host tests do not establish
successful compilation or compression round trips.
