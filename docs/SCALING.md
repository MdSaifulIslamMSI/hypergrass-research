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

## What 1000x actually costs

100 TB is 1000x the 100 GB class. At the measured 2.58 MB/s single-thread rate that is
**10,800 core-hours** (449 core-days), and the compute is the cheap part: order
$10^2-10^3 at spot prices. The binding resources are elsewhere.

| Shard | Shards for 100 TB | Wall per shard | Cores for a 1 h wall |
|---|---:|---:|---:|
| 10 GB | 10,000 | 1.1 h | ~10,800 |
| 100 GB | 1,000 | 10.8 h | unreachable; the shard wall is 10.8 h |

Three costs dominate, none of them FLOPs:

1. **Read amplification.** The bounded path reads every byte four times:
   `partitionFile`, anchor build, one consensus pass, encode. 100 TB becomes
   400 TB of reads, so a one-hour wall needs ~110 GB/s of aggregate storage
   bandwidth. Verification is a second 100 TB workload: each shard is restored
   and compared, which writes 100 TB and reads it back unless verification is
   fused with decompression.
2. **Shard-local anchors.** Each shard builds its own bank, so the 10 GB ratio is
   not the 120 MB ratio. A fixed anchor budget covers a falling fraction of a
   growing shard. This is the one number that decides whether 1000x costs
   archive size, and it has not been measured.
3. **Resumability and verification.** A 10.8 h shard that dies at hour ten must
   resume, not restart, and every shard must be independently byte-exact.
   Neither exists yet.

More threads inside a shard do not fix this: the measured 2-thread speedup is 1.15x
because only the encode pass is parallel, while the partition scan, anchor build and
consensus pass are serial. Massive core count helps only after the work is
partitioned, and partitioning is a format change (Phase 4), not a scheduler setting.

A 1000x *single-stream* speedup is a different claim: 2.58 MB/s to 2.6 GB/s is
~1.2 cycles per base at 3 GHz, below the cost of the range coder alone, and the serial
passes cap thread scaling near 4x. That needs a different compression model, not more
cores.

## Current hard limits

| Limit | Where | Consequence |
|---|---|---|
| Unbounded path reads the whole input into one `byte[]` | `HyperGraSS.compressFile` | cannot exceed ~2 GiB |
| Bounded path used `int` block offsets | `FastaStreamFactorizer.BlockSlice` | fixed: fields are now `long` |
| `partitionFile` rejected inputs >2 GiB | `FastaStreamFactorizer.partitionFile` | fixed: guard removed |
| Archive is assembled in one `ByteOut` (`byte[]`) | `HyperGraSS.serializeMultiBlockArchive` | archive must fit ~2 GiB; all encoded blocks are retained |
| In-memory decoder rejects >2 GiB | `HyperGraSS.decompressBytes` | `decompressFile` streams and is the supported path |
| Random-access query uses `int` record positions | `HyperGraSS.query` | a single record >2 GiB cannot be queried |
| Global anchor bank + seed index | `CanonicalSyncmerEngine.AnchorBankBuilder` | memory grows with pangenome diversity, not input size: a saturated 9 MiB anchor carries a ~67 MB index and exhausts a 160 MB heap on a 123 MB diverse pangenome |

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

### The bound is set by anchor saturation, not input size

Gate 4 passes on the probe above because that input is 55 copies of one sequence: the
anchor bank never saturates. It does not generalize to diverse input.

The bank costs about 1 byte per base. Its seed index costs about **7.1 bytes per base**:
`SEED_STEP = 4` indexes both strands, and the open-addressed table runs at 70% load
with 8-byte slots (9,437,184 bases -> 4,718,592 entries -> 8,388,608 slots x 8 B =
67.1 MB). At the 9 MiB cap that is 9 MB of bank plus 67 MB of index. With a 16 MiB raw
block and its 16 MiB DNA array, the builder holds ~115 MB live, and a 160 MB heap has
no room for the next block.

Measured on the 122,948,162 B *E. coli* medium stratum (24 genomes, 9 MiB anchor,
16 MiB blocks):

- `OutOfMemoryError` in `FastaStreamFactorizer.factorize` at 18.7 s in the container
  and 6.0 s on the host JVM;
- heap 154 MB with 112 MB live at the failure, per `-Xlog:gc`;
- the same input on the same 160 MB heap with `--anchor-mib 4` runs past that point
  (still alive at 60 s, archive streaming), so the index is the cause.

A diverse 123 MB pangenome fails where a repetitive 2.08 GB one passes. The bound
tracks novelty, not bytes. For the gate this means the 9 MiB anchor and the 160 MB heap
cannot both hold on diverse input; either cap the anchor (about 34 MB of index at 4 MiB,
at a size cost) or raise the heap past the 256 MiB RSS budget. That trade is measured
before the published profile changes.

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
three factorization passes over the input (anchor build, one consensus pass, encode)
plus the partitioning scan, and only the encode pass is parallel. Threads alone
therefore cap out well below linear, which is why sharding - not more cores - is the
Tier B answer.

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
to tune away: the partition scan, anchor build and consensus pass are still serial, so
Amdahl caps this path. The encode pass is the only parallel stage. That is the argument
for Phase 4 sharding rather than for more threads: with four passes over the input and
one of them parallel, no thread count rescues Tier B.

The +0.79% archive growth is the shared gap channel being switched off, recorded in
`FLAG_SHARED_GAP` rather than hidden in a JVM flag.

### Where the wall time goes

Temporary stage instrumentation on the 40.5 MB campy-high input (2 blocks, host JVM,
same 9 MiB anchor), removed before commit. Shares are more useful than absolute times
because host and container wall times differ by about 2.5x.

| Stage | Time | Share |
|---|---:|---:|
| `partitionFile` | 0.61-0.63 s | 5% |
| anchor build | 1.06-1.10 s | 9% |
| consensus (1 pass) | 1.02-1.05 s | 9% |
| encode | 8.68-8.71 s | 76% |
| total | 11.4-12.1 s | |

Inside the encode stage, block factorization plus encoding is only 3.65 s (1.50 s for
the 16 MiB block, 2.15 s for the 22.7 MiB block). The remaining ~5 s is the anchor
header. It is not a trial: the header must emit one encoded anchor. With
`--anchor-single-page` that payload is one 2,800,646 bp context-mixing stream, which
costs 4.94 s; the fast 2-bit encode it beats costs 0.25 s. The loser is the only
overhead, so sampling the decision would save 0.25 s, not 5 s.

The real trade is anchor page size and codec:

| Anchor mode | Archive | Header | Total | Δ size |
|---|---:|---:|---:|---:|
| `--anchor-single-page` (one CM page) | 1,125,629 B | 5.19 s | 11.70 s | — |
| `--anchor-cm-pages` (64 KiB CM pages) | 1,130,790 B | 2.38 s | 8.85 s | +0.46% |
| default 64 KiB fast pages | 1,151,824 B | 0.20 s | 6.85 s | +2.33% |

Host JVM, 1 GiB heap, 40.5 MB campy-high. At the 160 MiB bounded heap the same three
modes total 33.1 s, 25.1 s and 17.4 s: the ordering holds, GC makes each slower.

In the locked 256 MiB / 2-CPU container on the 120.1 Mbp *E. coli* high input, the
same option costs +0.489% (4,636,907 -> 4,659,584 B) and produced the same archive
bytes on three runs. Tuned compress wall was 144.2 and 158.4 s; paged was 112.3,
116.5 and 122.8 s. Decompress fell from 66.3-71.4 s to 40.3-53.3 s. Peak RSS was
272.5-273.0 MiB for tuned and 272.0-272.6 MiB for paged, so the option does not
change memory: this input sits at the 256 MiB cap in both modes. Absolute wall time
swings by 3x between container sessions (the same tuned profile took 49.7 s in the
broad panel), so only same-session comparisons are used.

`--anchor-cm-pages` is the middle trade: it keeps the context-mixing codec, resets it
every 64 KiB page, and costs 0.46% in size for 1.32x in wall. Parallel CM pages were
slower (9.71 s vs 8.85 s serial at 1 GiB) because every page allocates its own model
and the common pool contends, so CM pages are encoded serially. Fast pages are the
fastest option and already the default when `--anchor-single-page` is omitted. The
published bounded profile keeps `--anchor-single-page` so the 256 MiB gate and the
broad-panel archives stay comparable; `--anchor-cm-pages` is opt-in.

Ordered by measured impact, with the corrected accounting:

1. **Anchor header** (44% of wall). Three measured modes: `--anchor-cm-pages` gives
   1.32x for +0.46%; fast pages give 1.71x for +2.33%. Sampling the codec decision is
   not worth it: it saves 0.25 s, not 5 s.
2. **Block encode** (32% of wall). Already parallel; scales with the number of blocks
   and the per-worker memory budget, not with the anchor.
3. **Anchor build + consensus** (19% of wall). Two serial reads. Fuse them, or build
   the bank per shard so the work partitions with the data.
4. **Partition scan** (5% of wall). Fuse it with the anchor build so the input is read
   once instead of twice.

The header cannot be removed, only made cheaper. Even the fast-page floor is 6.85 s,
and the remaining serial reads cap threads; sharding is still the structural fix, and
no single-node mode reaches 2x on this input.

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
