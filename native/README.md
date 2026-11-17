# EIMDL / HyperGraSS Native Tri-Language Engine (`C23 + C++23 + Rust 2024`)

Zero-external-crate native `.hgs2`-wire-compatible implementation of the Generation 3 **EIMDL / `HGS2`** engine, interoperating bidirectionally (SHA-256 identical) with Java [`EimdlCli`](../src/hypergrass/EimdlCli.java).

## Architecture

| Layer | Files | Responsibility |
|---|---|---|
| **C23 SIMD & Entropy Kernels** | [`c_kernels/hgs_simd.h`](c_kernels/hgs_simd.h), [`c_kernels/hgs_simd.c`](c_kernels/hgs_simd.c), [`c_kernels/hgs_deflate.h`](c_kernels/hgs_deflate.h), [`c_kernels/hgs_deflate.c`](c_kernels/hgs_deflate.c), [`c_kernels/hgs_ans.h`](c_kernels/hgs_ans.h), [`c_kernels/hgs_ans.c`](c_kernels/hgs_ans.c) | Runtime-CPUID AVX2 & SWAR 2-bit DNA pack/unpack, C23 `<stdbit.h>` / `constexpr` LUTs / `[[nodiscard]]`, 64-bit SWAR reverse-complement, open-syncmer scan, Slice-by-8 IEEE CRC-32, SHA-256, RFC 1951/1950 inflate/deflate, and **4-stripe interleaved Byte-rANS** for DNA (Order-0/2/4/6 + Law IV raw 2-bit bypass) & pangenome-conditioned Phred quality scores |
| **C++23 64-Bit Mosaic, Novelty & Inspect Engine** | [`cpp_engine/hgs_cpp_engine.hpp`](cpp_engine/hgs_cpp_engine.hpp), [`cpp_engine/hgs_cpp_engine.cpp`](cpp_engine/hgs_cpp_engine.cpp) | `HgsSyncmerIndex64` (>2.14 Gbp), `<bit>` intrinsics, parallel Two-Gear Kimura-$V_4$ + syntenic splice + self-repeat `BitRange` Mosaic codec, `FastNoveltyCoder`, parallel multi-segment `4x`-rANS novelty wrapper, and **zero-target-allocation compressed-domain Mosaic token inspector** (`hgs_cpp_mosaic_inspect`) |
| **Rust 2024 Streaming CLI & Engines (`eimdl-rs`)** | [`eimdl-rs/Cargo.toml`](eimdl-rs/Cargo.toml), [`eimdl-rs/build.rs`](eimdl-rs/build.rs), [`eimdl-rs/src/ffi.rs`](eimdl-rs/src/ffi.rs), [`eimdl-rs/src/fasta.rs`](eimdl-rs/src/fasta.rs), [`eimdl-rs/src/hgs2.rs`](eimdl-rs/src/hgs2.rs), [`eimdl-rs/src/hgs3.rs`](eimdl-rs/src/hgs3.rs), [`eimdl-rs/src/compute.rs`](eimdl-rs/src/compute.rs), [`eimdl-rs/src/main.rs`](eimdl-rs/src/main.rs) | Rust 2024 Edition (`rust-version = "1.85"`, zero external crates) CLI & `.hgs2`/`.hgs3` container muxers, lossless FASTQ + Phred quality coupling + CDS ORF digest (`hgs3.rs`), compressed-domain `analyze` / VCF v4.2 `variants` / exact `kmer` / slice `query` (`compute.rs`), and SHA-256 reference panel pinning |

## Offline Build & Unit Tests

- Linux / macOS / MSYS2 (`gcc`/`g++` or `clang`/`clang++` supporting C23 + C++23 + `rustc 1.85+`), or Windows (`zig` C23/C++23 fallback + `cargo`):

```bash
cd native/eimdl-rs
cargo clippy --all-targets --offline -- -D warnings
cargo build --release --offline
cargo test --release --offline
```

## CLI Usage

```bash
# Single genome (Mode S)
./target/release/eimdl compress --single --target genome.fa --out genome.hgs2 --threads 4
./target/release/eimdl decompress --archive genome.hgs2 --out restored.fa --threads 4

# Reference-conditioned panel (Mode R)
./target/release/eimdl compress --target target.seq --panel panel_dir/ --exclude target.seq --out target.hgs2
./target/release/eimdl decompress --archive target.hgs2 --panel panel_dir/ --out restored.seq

# Pangenome collection / bounded memory (Mode C)
./target/release/eimdl compress --collection --bounded --target pangenome.fa --out pangenome.hgs2 --threads 4
./target/release/eimdl decompress --archive pangenome.hgs2 --out restored.fa --threads 4
./target/release/eimdl inspect --archive pangenome.hgs2

# Multi-omic HGS3 (lossless FASTQ with coupled Phred + 4x-rANS, or FASTA + CDS ORF digest)
./target/release/eimdl compress-omics --target reads.fastq --panel panel_dir/ --cds --out reads.hgs3 --threads 4
./target/release/eimdl decompress-omics --archive reads.hgs3 --panel panel_dir/ --out restored.fastq --threads 4

# Compressed-domain analytics, VCF v4.2 variant extraction, exact k-mer search & random slicing
./target/release/eimdl analyze --archive target.hgs2 --panel panel_dir/ --json
./target/release/eimdl variants --archive target.hgs2 --panel panel_dir/ --out calls.vcf
./target/release/eimdl kmer --archive target.hgs2 --panel panel_dir/ --query ACGTACGTACGT --both-strands
./target/release/eimdl query --archive target.hgs2 --panel panel_dir/ --unit 0 --start 1000 --end 1200
```
