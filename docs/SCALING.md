# Scaling HyperGraSS Beyond 2 GiB

Status: active. This note records the scaling limits of the current Java codec, the
order in which they should be removed, and the evidence required before each step is
accepted. It is deliberately language-neutral: the first three blockers are format and
memory-architecture problems, not compiler problems.

## Target scale

| Tier | Input | Deployment | Acceptance |
|---|---|---|---|
| A | 10 GB | one node, 2-8 threads, 256 MiB-1 GiB RSS | byte-exact round trip; bounded RSS; no archive >2 GiB limit |
| B | 1 TB | sharded archive, external-memory anchor/index | bounded RSS per shard; shard-parallel encode; random access by record |
| C | 100 TB | distributed storage + index service | out of scope for a single codec process; requires a partitioned format |

Tier A is the next milestone. Tier B requires a format change (sharding), not a language
change. Tier C is a storage-system design, not a compression-kernel design.

## Current hard limits

| Limit | Where | Consequence |
|---|---|---|
| Unbounded path reads the whole input into one `byte[]` | `HyperGraSS.compressFile` | cannot exceed ~2 GiB |
| Bounded path used `int` block offsets | `FastaStreamFactorizer.BlockSlice` | fixed: fields are now `long` |
| `partitionFile` rejected inputs >2 GiB | `FastaStreamFactorizer.partitionFile` | fixed: guard removed |
| Archive is assembled in one `ByteOut` (`byte[]`) | `HyperGraSS.serializeMultiBlockArchive` | archive must fit ~2 GiB; all encoded blocks are retained |
| In-memory decoder rejects >2 GiB | `HyperGraSS.decompressBytes` | `decompressFile` streams and is the supported path |
| Random-access query uses `int` record positions | `HyperGraSS.query` | a single record >2 GiB cannot be queried |
| Global anchor bank + seed index | `CanonicalSyncmerEngine.AnchorBankBuilder` | memory grows with pangenome diversity, not input size |

The archive footer already stores `totalRawBytes` as a `varLong` and block file offsets as
`long` deltas, so the on-disk format does not need a version bump for Tier A.

## Phase 0 - measure before changing anything

Instrument one bounded run at 1 GB, 2.2 GB and 10 GB with:

- wall time, CPU time, peak RSS;
- bytes read/written (including the extra anchor and vote passes);
- time in partition, factorize, anchor build, encode, serialize, verify;
- archive bytes by channel (`ArchiveByteAccounting` once it supports the current flags).

Gate: a written baseline table. No optimization is accepted without a before/after on the
same locked image and the same input.

## Phase 1 - 64-bit bounded ingestion (this change)

- `BlockSlice.startByte/endByte` are `long`.
- `partitionFile` accepts files larger than 2 GiB.
- `factorizeSlice` rejects only a single block larger than the in-memory `int` limit.
- The in-memory path (`partitionIntoBlocks`, `compressBytes`) intentionally stays 32-bit;
  it is not the scalable path.

Gate:

1. all existing verification suites pass;
2. a 40 MB recompression is byte-identical to the pre-change archive (format compatibility);
3. a >2 GiB synthetic input compresses and restores with SHA-256 equality;
4. peak RSS stays at or below the 256 MiB bounded budget.

### Phase 1 measurements (2026-11-03, host JVM, bounded profile)

Input is 55 concatenated copies of the campy-high panel sequence. It is a **scale
probe, not a ratio probe**: the 73.9x comes from repeated content and must not be
quoted as a compression result.

| Metric | Value |
|---|---|
| Input | 2,229,710,670 B (2.077 GiB) |
| Archive | 30,171,539 B |
| Blocks / records | 132 / 1 |
| Anchor | 2,800,811 bp, 1 page |
| Compress wall | 864.52 s (2.58 MB/s) |
| Decompress wall | 238.19 s (9.36 MB/s) |
| Peak sampled working set (compress) | 235 MB |
| Peak heap cap | `-Xmx160m` |
| CRC-32 | `0xEE111322`, verified on restore |
| SHA-256 in / out | `3078E87C…4B54F7` / identical |

Byte compatibility on the 40.5 MB campy-high input: `1,125,629 B`, SHA-256
`14DD7E94…E1C5D6`, identical to the pre-change archive, wall 42.16 s.

Gates 2, 3 and 4 pass. Gate 1 now passes as well: `VerifyArchiveAccounting` was extended
to parse DNA modes 4-6 and auxiliary modes 4-7, so the suite reports
`ALL 34 ACCOUNTING SUITES PASSED (3326 checks, 1721 malformed archives rejected)`.

**The memory ceiling is gone; throughput is now the binding limit.** Extrapolating the
measured single-threaded rate:

| Tier | Input | Single-thread wall |
|---|---|---|
| A | 10 GB | ~65 min |
| B | 1 TB | ~4.5 days |
| C | 100 TB | ~15 months |

That is why parallel blocks and shards matter more than the archive writer. Phase 2
removes the archive-size ceiling; it does not make 1 TB finish in a useful time.

Phase 3's first measurement also sets the expectation honestly: the bounded path makes
four passes over the input (anchor build, two consensus passes, encode) and only the
last is parallel. Threads alone therefore cap out well below linear, which is why
sharding - not more cores - is the Tier B answer.

Phase 1 left one risk open: the archive writer materialized every encoded block. A
2.2 GB input with a compressible payload stays under 2 GiB of archive, but a 10 GB
high-entropy input would not. That is closed by Phase 2 below.

## Phase 2 - streaming archive writer (implemented)

`compressFileBounded` now writes header, blocks and footer straight to a
`FileChannel`. `writeArchiveHeader` is shared with the in-memory serializer, so
the header bytes are produced by one code path. Peak resident payload is one
block instead of the whole archive.

- a `FileChannel` writer emits the header and each block as it is encoded;
- the footer is built after the last block from the running offset;
- the block directory still carries `long` file offsets;
- inputs below 2 GiB produce byte-identical archives (verified on campy-high).

Gate: 10 GB input at bounded RSS; archive >2 GiB round-trips; byte-identical
output to the pre-change writer on inputs below 2 GiB.

Status: byte-identity and the >2 GiB *input* round trip pass. The campy-high archive
rebuilt through the streaming writer is `1,125,629 B` with SHA-256
`14DD7E94…E1C5D6`, identical to the pre-change archive. The >2 GiB *archive* leg
still needs a >=8 GB high-entropy input, which is the Tier A acceptance run.

## Phase 3 - parallel block pipeline (implemented, `--threads N`)

`--threads N` encodes a window of `2N` blocks in parallel and appends them in
block order. The default is `1`, which is the published bounded profile and
keeps the cross-block gap channel (`FLAG_SHARED_GAP`). With `N > 1` the flag is
cleared: each block is encoded independently and the decoder builds no gap
history. That is a format-level decision recorded in the archive, not a hidden
JVM setting.

Gate: 2-4x wall-time improvement on 2-8 cores; identical round trip; RSS =
base + N x per-worker budget; sequential output unchanged.

Status: implemented and covered by `VerifyAllEdgeCases` (both flag states verified, plus
determinism of the sequential writer). First measurement on the 40.5 MB campy-high input,
which has only 2 blocks:

| Mode | Archive | Wall | Speedup |
|---|---|---|---|
| `--threads 1` | 1,125,629 B | 30.6 s | 1.00x |
| `--threads 2` | 1,134,556 B | 26.6 s | 1.15x |

A 1.15x end-to-end win from 2 threads on 2 blocks is the expected shape, not a defect
to tune away: the anchor build and the two consensus passes are still serial, so
Amdahl caps this path. The encode pass is the only parallel stage. That is the argument
for Phase 4 sharding rather than for more threads: with four passes over the input and
one of them parallel, no thread count rescues Tier B.

The +0.79% archive growth is the shared gap channel being switched off, recorded in
`FLAG_SHARED_GAP` rather than hidden in a JVM flag.

## Phase 4 - sharded anchors

At Tier B a single global anchor is the wrong data structure. Partition by clade or by
record range, build one anchor per shard, and store a manifest with:

- shard boundaries (record index and raw offset);
- per-shard anchor and block directory;
- a global record index for random access.

Gate: 100 GB encode at bounded RSS; per-shard round trip; a query resolves through the
manifest without loading other shards.

## Phase 5 - resumability

Add a checkpoint manifest written after every completed block. On restart, validate the
manifest and continue from the first incomplete block.

Gate: kill the process at 50%; resume reprocesses less than 5% of the input; final
archive is byte-identical to an uninterrupted run.

## Phase 6 - native kernels (optional, evidence-gated)

Only after Phase 0 identifies a hot loop with a measured gap of at least 1.5x versus a
tuned Java implementation. Extract that kernel, not the codec:

- **Rust** if parser safety, no-GC memory control and a standalone CLI matter most;
- **C++** if the deciding factor is an existing library (htslib, libdeflate, zstd, ISA-L);
- **never both** in the same release path.

Gate: same locked container, byte-identical archives, 1.5x on the isolated kernel and a
measurable end-to-end win. If the kernel does not clear the gate, keep Java.

## Decision log

- A Java + C++ + Rust split is rejected. Three build systems, three ABIs and three
  conformance surfaces for one format is a cost multiplier, not a scaling strategy.
- The 2 GiB ceiling is a format/memory-architecture bug. Fixing it is worth more than a
  native rewrite of the hot loops.
- Language is a 1.5-3x lever. Removing the 32-bit ceiling, streaming the archive and
  sharding the anchor are 10-1000x levers.
- Byte-exact verification remains the non-negotiable gate for every phase.
