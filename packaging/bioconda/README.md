# Bioconda Packaging Recipe (`packaging/bioconda/eimdl/`)

This directory contains the reproducible **Bioconda** recipe for installing **`eimdl`** (zero-crate Native `C23 / C++23 / Rust 2024` engine), **`hypergrass`**, and **`eimdl-java`** (`Java 25 LTS` `.jar` wrappers) completely offline from tagged release tarballs.

## Files

| File | Purpose |
|---|---|
| [`meta.yaml`](eimdl/meta.yaml) | Conda/Bioconda recipe pinning `{{ compiler('c') }}`, `{{ compiler('cxx') }}`, `{{ compiler('rust') }} >=1.85`, `openjdk >=25`, and immutable release tarball + `SHA-256` |
| [`build.sh`](eimdl/build.sh) | Offline build installing `$PREFIX/bin/eimdl`, `$PREFIX/share/hypergrass/hypergrass.jar`, `$PREFIX/bin/hypergrass`, and `$PREFIX/bin/eimdl-java` |
| [`run_test.sh`](eimdl/run_test.sh) | End-to-end Bioconda post-install verification testing CLI versions, bidirectional Native $\leftrightarrow$ Java `.hgs3` round-trip SHA-256 parity, compressed-domain `analyze`, and lazy zero-target `kmer` |

## Why Offline Bioconda Builds Work Without Vendoring
- **Native Engine ([`Cargo.toml`](../../../native/eimdl-rs/Cargo.toml)):** `0` external Rust crates — `cargo build --release --offline` builds hermetically inside offline Bioconda containers with zero `cargo vendor` bundle bloat.
- **Java Engine ([`src/hypergrass/`](../../../src/hypergrass/)):** `0` external JAR dependencies — `javac --release 25` + `jar` builds directly against `conda-forge`'s `openjdk >=25` (`25.0.2+`).
- **Compiler Env Support ([`build.rs`](../../../native/eimdl-rs/build.rs)):** Explicitly honors Bioconda's `$CC`, `$CXX`, and `$AR` compiler wrappers (with a `gcc`/`g++` shim in `build.sh` for prior source tags).

## Submitting to `bioconda/bioconda-recipes`
1. Ensure [`meta.yaml`](eimdl/meta.yaml) `version` and `sha256` point to the desired published GitHub Release tarball (currently `v3.0.0`, ready to bump cleanly to `v3.1.0` after Kaggle v7 completes).
2. Fork `bioconda/bioconda-recipes`, copy `packaging/bioconda/eimdl/` into `recipes/eimdl/`, lint with `bioconda-utils lint recipes/eimdl`, and open a Pull Request against `bioconda:master`.
