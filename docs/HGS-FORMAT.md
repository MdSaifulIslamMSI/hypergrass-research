# HyperGraSS Archive Container Specifications (`.hgs2` v2.0 & `.hgs` v1.0)

Status: Production implementation-defined binary container specifications documented from the unified **Generation 3 `HGS2` / `EIMDL` engine** ([`HyperGraSSFusionEngine.java`](../src/hypergrass/HyperGraSSFusionEngine.java), [`EimdlMosaicCodec.java`](../src/hypergrass/EimdlMosaicCodec.java), [`FastNoveltyCoder.java`](../src/hypergrass/FastNoveltyCoder.java)) and the **Generation 1 `HyperGraSS` 1.0.0 engine** ([`HyperGraSS.java`](../src/hypergrass/HyperGraSS.java)). Both decoders enforce strict structural validation, CRC-32 verification, and (for Mode R in `.hgs2`) cryptographic 32-byte `SHA-256` donor-panel manifest verification before emitting output.

---

# Part I — Unified `.hgs2` Container Format (Version 2, Magic `"HGS2"`)

## 1. Overview & Regime Separation

An `.hgs2` archive begins with the 4-byte ASCII magic `48 47 53 32` (`"HGS2"`) and losslessly encodes any raw byte stream (FASTA, multi-FASTA, RNA `A/C/G/U`, mixed IUPAC/masked sequences, or headerless `.seq` streams) under one of three explicit operational regimes:

| Mode ID (`flags & 0x03`) | Regime Name | Reference Conditioning | Description |
|---:|---|---|---|
| `0` (`MODE_SINGLE`) | **Mode S** | $\mathcal{R} = \emptyset$ (self-contained) | Single genome, chromosome, or synthetic sequence encoded via bidirectional ($\pm$) self-mosaic, tandem periodicity, and `FastNoveltyCoder` (`6F-ORF-HMM` + `FastCRM`). |
| `1` (`MODE_REF`) | **Mode R** | $\mathcal{R} \neq \emptyset$ (external panel) | Reference-conditioned delta archive relative to an external donor panel pinned in the header by a compressed donor-name directory and a mandatory **32-byte `SHA-256` Merkle manifest digest**. |
| `2` (`MODE_COLLECTION`) | **Mode C** | $\mathcal{R} = \emptyset$ (self-contained) | Multi-genome pangenome collection ordered along a causal syncmer-Jaccard Prim minimum spanning tree (`rankToUnit[r]`, `maxAncestorRank[r]`) with a unified pan-genome novelty superstream. |

## 2. Top-Level `.hgs2` Physical Layout

```text
offset  type                   field                        condition
------  ---------------------  ---------------------------  ----------------------------------------
0       bytes[4]               magic = "HGS2" (48 47 53 32) always
4       byte                   headerFlags                  bits 0..1 = mode (0,1,2);
                                                            bit 2 (0x04) = FLAG_RAW_FLAT_ACGT;
                                                            bit 3 (0x08) = FLAG_RAW_RNA_U;
                                                            bits 4..7 must be 0 (else rejected)
5       longLE (8 B)           rawByteLen                   exact uncompressed file byte length (>= 0)
13      intLE  (4 B)           rawCrc32                     IEEE 802.3 CRC-32 of uncompressed bytes
17      uvarint                totalRecords                 number of FASTA records (or 1 if flat)
...     uvarint                numUnits                     number of coding units U (>= 0)
...     [PanelPinningBlock]    modeRPinning                 present iff mode == 1 (MODE_REF)
...     [CollectionDagBlock]   modeCDagPermutation          present iff mode == 2 && numUnits > 1
...     [UnitFramingBlock]     unitMetadataOrLengths        present iff numUnits > 0
...     [RankMosaicStreams]    mosaicPayloads[0..U-1]       present iff numUnits > 0 (in rank order)
...     [NoveltySuperstream]   unifiedNoveltyPayload        present iff numUnits > 0
```

A valid `.hgs2` decoder (`HyperGraSSFusionEngine.parseHeaderAndDirectory`) rejects any archive shorter than 18 bytes, any non-zero bit in `headerFlags & ~0x0F`, any `mode > 2`, any negative length, or any trailing byte after `unifiedNoveltyPayload`.

## 3. Mode R Cryptographic Reference-Panel Pinning (`mode == 1`)

Present immediately after `numUnits` when and only when `mode == 1`:

```text
panelCount          uvarint      number of donor sequences M in the external reference panel
compNamesLen        uvarint      byte length of Deflate-compressed donor name table
compNames           bytes[n]     Deflate-compressed UTF-8 donor filenames (newline-delimited)
panelManifestSha256 bytes[32]    256-bit SHA-256 Merkle digest over ordered donor names,
                                 base lengths, and 2-bit packed DNA word banks
```

Prior to decoding any mosaic token (`verifyAndLoadPanelForDecode`), the decoder loads the external `--panel` directory or file, resolves the exact ordered donor list `panelNames[0..M-1]`, recomputes the 32-byte `SHA-256` manifest digest, and aborts immediately with an `IllegalStateException` if any donor is missing, reordered, or modified by even a single base.

## 4. Mode C Causal Phylogenetic DAG Permutation (`mode == 2 && numUnits > 1`)

When compressing a multi-unit collection in Mode C, `buildPhylogeneticPrimOrder` computes bottom-2048 canonical open-syncmer sketches for all $U$ units, selects the medoid root $r = 0$, and orders $r = 0, \dots, U-1$ by Prim's maximum-similarity spanning tree:

```text
rankToUnit[0..U-1]       uvarint[U]   bijective permutation mapping phylogenetic rank r -> input unit u
maxAncestorRank[0..U-1]  uvarint[U]   causal ancestor cap for rank r (must satisfy 0 <= anc <= r)
```

The decoder verifies that `rankToUnit` is a valid permutation of $\{0, \dots, U-1\}$ and that `0 <= maxAncestorRank[r] <= r` for every rank $r$, guaranteeing strictly causal DAG reconstruction without cycles.

## 5. Per-Unit Framing & 6-Channel Factorized Metadata (`numUnits > 0`)

- **Case A — Fast Flat Sequence Path (`FLAG_RAW_FLAT_ACGT = 0x04` set):**
  Used when the input is a headerless, newline-free uppercase `ACGT` (or `ACGU` when `FLAG_RAW_RNA_U = 0x08` is set) stream streamed directly via `tryStreamFlatFile`:
  ```text
  unitDnaLens[0..U-1]    uvarint[U]   exact base length of each input unit u
  ```
- **Case B — General 6-Channel Factorized Path (`FLAG_RAW_FLAT_ACGT` clear):**
  Used for all FASTA, multi-FASTA, soft-masked, `N`-scaffolded, or IUPAC-containing inputs factorized by `FastaStreamFactorizer`:
  ```text
  for each input unit u in 0 .. U-1:
    firstRecordIdx       uvarint      0-based index of the first FASTA record in unit u
    recordCount          uvarint      number of FASTA records in unit u
    unitMetaLen          uvarint      byte length of serialized 6-channel metadata
    unitMetaBytes        bytes[mLen]  6-channel metadata stream (identical to Section 11.1 below)
  ```

## 6. Rank-Ordered `EimdlMosaicCodec` Streams (`mosaicPayloads[0..U-1]`)

For each phylogenetic rank $r \in \{0, \dots, U-1\}$:

```text
mosaicLen                uvarint      byte length of rank r's encoded mosaic stream
mosaicBytes              bytes[mLen]  self-contained EimdlMosaicCodec stream for rank r
```

Inside each `mosaicBytes` payload:
1. **Mosaic Header:**
   - `configByte` (`1 byte`): encodes active EIMDL law flags (`allowApprox`, `allowSplice`, `allowV4`, `allowWobble`, `allowSelf`, `allowTandem`, `hasPanel`).
   - `targetBases` (`uvarint`): exact DNA base length $n_r$ reconstructed by rank $r$.
   - `noveltyBases` (`uvarint`, present iff `targetBases > 0`): exact number of novel/literal + insertion-splice bases consumed by rank $r$ from the unified novelty superstream.
2. **Bitwise Order-1 Adaptive Range-Coded Phrase & Evolutionary Residual Stream (`BitRangeEnc`):**
   - Phrase kind decision trees (`LITERAL`, `PANEL` copy, `SELF` copy, `TANDEM` period $1 \le P \le 64$).
   - Panel donor ID (conditioned on `preferredTwinDonor` and similarity order `mtf[0..3]`), strand bit ($\pm$), reference coordinate delta or global position, and token span.
   - **Law IIIc Syntenic Micro-Indel Splices:** `hasSplice` flag conditioned on 8-state `(kind, spanBin)` prior, splice sub-block length, and $(\Delta_{\text{tgt}}, \Delta_{\text{ref}})$ insertion/deletion offsets.
   - **Law IIa Bijective Codon-Wobble Substitution Gaps:** quotient $q = \lfloor \text{gap}/3 \rfloor$ and codon-phase remainder $r = (\text{gap}+1)\bmod 3$ conditioned on local divergence regime.
   - **Law I Klein Four-Group ($V_4$) Substitution Operators:** non-zero XOR operator $\delta = y \oplus x \in \{1, 2, 3\}$ conditioned on reference trinucleotide context $c = (x_{-1}, x_0, x_{+1})$ (including $\text{CpG}$ hypermutability) and divergence regime:
     - First binary decision: **Kimura Transition** ($\delta = 2$) vs. **Transversion** ($\delta \in \{1, 3\}$).
     - Second binary decision (if transversion): **Keto-Amino** ($\delta = 1$) vs. **Watson-Crick Complement** ($\delta = 3$).

## 7. Unified Pan-Genome Novelty Superstream (`FastNoveltyCoder`)

All novel literal runs and insertion-splice bases from ranks $r = 0, \dots, U-1$ are concatenated in rank order into a single contiguous buffer of length `totalNoveltyBases` ($\sum_{r=0}^{U-1} \text{noveltyLen}_r$) and encoded at the end of the `.hgs2` archive:

```text
totalNoveltyBases        uvarint      total novelty bases N_nov across all ranks
noveltyPayloadLen        uvarint      byte length of the FastNoveltyCoder payload
noveltyPayload           bytes[nLen]  segmented FastNoveltyCoder stream
```

Inside `noveltyPayload` (when `totalNoveltyBases > 0`):

```text
totalSegs                uvarint      1 (serial) or 1 + workerCount (warm-primed parallel stripes)
modelFlags               byte         bit 0 (0x01) = boundedMem; bit 1 (0x02) = panSuper
segLens[0..totalSegs-1]  uvarint[]    base count of each segment (sum must equal totalNoveltyBases)
compLens[0..totalSegs-1] uvarint[]    compressed byte length of each segment
segPayloads[0..totalSegs-1]:
  segMode                byte         0x00 = SEG_CRM (6F-ORF-HMM + multi-order Markov/repeat mixer)
                                      0x01 = SEG_RAW_2BIT (Law IV incompressible raw 2-bit bypass)
  segBody                bytes[...]   if SEG_CRM: carryless range-coded stream;
                                      if SEG_RAW_2BIT: ceil(segLen / 4) bytes of 2-bit packed DNA
```

When `totalSegs > 1`, segment `0` (`warmBases`) is decoded first to populate the shared `WarmSnapshot` (`6F-ORF-HMM` posteriors, Markov tables, and prefix repeat rings), after which segments `1 .. totalSegs-1` are decoded concurrently across worker threads.

---

# Part II — Legacy `.hgs` Container Format (Version 1, Magic `"HGS\x01"`)

Status: implementation-defined format, documented from the HyperGraSS 1.0.0
encoder and decoder (`src/hypergrass/HyperGraSS.java`). A decoder must reject
malformed input; the encoder emits only the canonical choices described here.

## 8. Scope

An `.hgs` file stores an arbitrary byte stream losslessly. When the input is
FASTA, the encoder additionally factorizes headers, line layout, soft-masking,
scaffold `N` runs, and non-ACGT symbols into separate channels, then compresses
the DNA stream with a stride-sampled 16-mer anchor model. The decoder reconstructs the
original bytes exactly, including headers, wrapping, case, and line endings.

The format has two physical layouts:

| Layout | Used for | Header |
| --- | --- | --- |
| Compact single-block | Small inputs or `--single` superblock | 13-byte header + one payload |
| Multi-block | General inputs | Header, optional anchor bank, blocks, footer, 12-byte trailer |

Both layouts start with the same 4-byte magic `48 47 53 01` (`"HGS"` + version
byte `0x01`). Readers must check the magic and the layout flag before parsing
anything else.

## 9. Primitive Encodings (Shared by `.hgs` and `.hgs2`)

- **`uvarint`** — unsigned LEB128, 7 bits per byte, least-significant group
  first, high bit set on every byte except the last. Values are limited to the
  range of the field's declared type.
- **`svarint`** — signed ZigZag followed by `uvarint`:
  `zz = (value << 1) ^ (value >> 31)` for 32-bit values.
- **`intLE` / `longLE`** — fixed-width little-endian 32-bit / 64-bit integers.
- **`bytes[n]`** — `n` raw bytes with no length prefix unless stated.
- **`aux`** — a length-delimited auxiliary channel (`VarIntCodec`):

  | First byte | Meaning |
  | --- | --- |
  | `0` | Empty. No further bytes. |
  | `1` | Raw. Followed by `uvarint` length and that many bytes. |
  | `2` | Deflate. Followed by `uvarint` uncompressed length, `uvarint` compressed length, and that many bytes of a raw DEFLATE stream (`nowrap=true`, level 9). |
  | `3..6` | Specialized entropy-floor codecs (RLE, Delta-ZigZag, Nibble-pack, Order-0/Order-1 Range Coding). |
  | `7` | Cross-block LZ77 + range-coded stream decoded against the preceding 256 KiB of gap plaintext (when `SHARED_GAP` is set). |

## 10. `.hgs` Common Header, Compact Layout & Multi-Block Layout

```text
offset  size  field
0       4     magic = 48 47 53 01
4       1     mode flags (bit field)
```

Mode flags:

| Bit | Name | Meaning |
| --- | --- | --- |
| `0x01` | `HAS_ANCHOR` | Multi-block layout with an anchor bank. |
| `0x02` | `LEAN_MODE` | Streaming/low-memory encoder settings. |
| `0x04` | `COMPACT_SINGLE` | Compact single-block layout. |
| `0x08` | `SHARED_GAP` | Block gap channels share a 256 KiB LZ history across blocks. |
| `0x10` | `COMPACT_TUNED` | Compact single archive using the tuned single-genome `QCRM` profile. |
| `0x20` | `QCRM_SEGMENTED` | Tuned `QCRM` payload split into independently coded segments. |

### 10.1 Compact Single-Block Layout

```text
offset  size  field
0       4     magic = 48 47 53 01
4       1     flags (includes 0x04 COMPACT_SINGLE)
5       4     rawLen, intLE
9       4     crc32, intLE
13      ...   payload
```

### 10.2 Multi-Block Layout

```text
offset  size  field
0       4     magic = 48 47 53 01
4       1     flags
5       ...   optional anchor bank (present iff HAS_ANCHOR)
...     ...   block payloads, concatenated
...     ...   footer body (aux-encoded)
...     8     footerOffset, longLE
...     4     FOOTER_MAGIC = 0x31464748 ("HGF1" in little-endian bytes)
```

### 10.3 Anchor Bank (`HAS_ANCHOR`)

```text
fast2Bit       byte        0 = 2-bit packed pages, 1 = context-mixer-coded pages
anchorLen      uvarint     total anchor bases
numPages       uvarint     number of pages (ANCHOR_PAGE_SIZE = 65,536 bases)
primarySeedLen uvarint
numDonors      uvarint
donorGaps      aux         per donor: uvarint length, then uvarint(pos - prevPos - 1)
donorDiffs     aux         one byte per element, low 2 bits used
pageBaseLen[]  uvarint     one per page, in page order
pageEncLen[]   uvarint     one per page, in page order
pagePayload[]  bytes       one per page, in page order
```

### 10.4 Footer Body

```text
numBlocks              uvarint
totalRawBytes          uvarint (64-bit)
totalRecords           uvarint
crc32                  intLE
entries[numBlocks]:
  offsetDelta          uvarint (64-bit)
  compressedLen        uvarint
  firstRecordDelta     uvarint
  recordCount          uvarint
  firstRecordBaseOffset uvarint (64-bit)
  blockSeqChars        uvarint
```

## 11. Block Payload & 6-Channel `FastaStreamFactorizer` Metadata

### 11.1 Auxiliary Metadata Channels (`FastaStreamFactorizer`)

Used identically by `.hgs` block payloads and `.hgs2` non-flat unit metadata:

```text
flags          byte        0x01 RNA mode, 0x02 headers, 0x04 lowercase, 0x08 N-runs, 0x10 exceptions
dnaLen         uvarint
totalSeqChars  uvarint
numSeqLines    uvarint
```

- **Channel 1 — Headers (iff `0x02`):** `uvarint numHeaders`, then `aux headerBlock` containing per-header `uvarint seqLineDelta` followed by front-coded `(sharedPrefixLen, suffixLen, suffixBytes)` (excluding leading `>` and EOL).
- **Channel 2 — Line Layout (always present):** `aux lineLengths` (RLE of sequence line widths) and `aux lineEols` (RLE of EOL codes: `0` = none, `1` = `LF`, `2` = `CRLF`).
- **Channel 3 — Lowercase Soft-Mask Intervals (iff `0x04`):** `uvarint numIntervals`, `uvarint minLength`, `aux gaps`, `aux lengths` (`length - minLength`).
- **Channel 4 — `N` Runs (iff `0x08`):** `uvarint numRuns`, `aux runs` of `(uvarint acgtGap, uvarint runLength)`.
- **Channel 5 — IUPAC Exception Symbols (iff `0x10`):** `uvarint numExceptions`, `aux gaps`, `aux symbols` (1 raw byte per exception).

### 11.2 `.hgs` DNA Payload Modes

| Mode | Layout |
| --- | --- |
| `0` | Empty DNA stream. No further bytes. |
| `1` | `uvarint compressedLen`, then that many bytes of a `QuantumContextRepeatMixer` stream that decodes to `dnaLen` bases. |
| `2` | Seed matches with per-match literals and SNP residuals (`CanonicalSyncmerEngine`). |
| `3` | Anchor-coordinate haplotype matrix with donor reuse (`CanonicalSyncmerEngine`). |

## 12. Reference Implementation Files

- [`src/hypergrass/HyperGraSSFusionEngine.java`](../src/hypergrass/HyperGraSSFusionEngine.java) — `.hgs2` container, Mode R SHA-256 manifest pinning, Mode C Prim DAG permutation, random-access query, byte accounting.
- [`src/hypergrass/EimdlMosaicCodec.java`](../src/hypergrass/EimdlMosaicCodec.java) — `.hgs2` Laws I–IV mosaic phrase, micro-indel splice, wobble bijection, and Klein $V_4$ substitution coder.
- [`src/hypergrass/FastNoveltyCoder.java`](../src/hypergrass/FastNoveltyCoder.java) — `.hgs2` `6F-ORF-HMM` + `FastCRM` parallel novelty superstream and `SEG_RAW_2BIT` safety bypass.
- [`src/hypergrass/HyperGraSS.java`](../src/hypergrass/HyperGraSS.java) — `.hgs` v1 container, footer, CLI, `--lean` streaming, and automatic `.hgs2` magic dispatch.
- [`src/hypergrass/FastaStreamFactorizer.java`](../src/hypergrass/FastaStreamFactorizer.java) — 6-channel lossless FASTA/RNA/DNA stream factorization.
- [`src/hypergrass/VarIntCodec.java`](../src/hypergrass/VarIntCodec.java) — LEB128 varints, 2-bit packing, and `aux` modes `0..7`.

---

## 13. `.hgs3` v3.0 Multi-Omic Container, 4-Stripe SIMD rANS, & Compressed-Domain Genomics

The `.hgs3` container (magic `"HGS3"`) generalizes `.hgs2` across **multi-omic biological streams** while maintaining 100% deterministic bidirectional interoperability between Native (`C23` + `C++23` + `Rust 2024`) and Java 25 LTS (`Hgs3OmicCodec`, `HgsAnsCoder`, `CompressedDomainEngine`):

1. **Header (`HGS3`):**
   - Bytes `0..4`: `"HGS3"`
   - Byte `4`: `flags` — bits `0..1` = mode (`0` Single, `1` Ref-Panel, `2` Collection MST), bit `2` (`0x04`) = `FLAG_RAW_FLAT_ACGT`, bit `3` (`0x08`) = `FLAG_RAW_RNA_U`, bit `4` (`0x10`) = `FLAG3_FASTQ`, bit `5` (`0x20`) = `FLAG3_HAS_CDS`
   - `rawBytes` (`uint64 LE`), `crc32` (`uint32 LE`), `totalRecords` (`uvarint`), `numUnits` (`uvarint`), plus Mode R SHA-256 panel manifest pin block (identical to Section 3).
2. **Parallel 4-Stripe Byte-rANS Novelty Coder (`FLAG_ANS4_NOVELTY = 0x10`, `hgs_ans.c` / `HgsAnsCoder.java`):**
   - Overflow-free 32-bit Byte-rANS (`PROB_BITS = 12`, `ANS_L = 2^23`, `XMAX_MUL = 524288`), partitioned into 4 independent stripes mixing adaptive Order-0/2/4/6 contexts across `1,000,000`-base deterministic segments, falling back strictly to `ANS_MODE_RAW_2BIT` (`0x00`) under Law IV noise ceilings.
3. **Lossless Pangenome-Conditioned FASTQ (`FLAG3_FASTQ = 0x10`):**
   - **Metadata Stream:** Common read-ID prefix stripping + incremental LCP deltas, uniform/mixed `CRLF`/`LF` & `+` comment flags, fixed/variable read lengths, and sparse/XOR ambiguity & soft-mask restoration (`deflate`/`zlib`).
   - **DNA Stream:** Reference/panel Mosaic (`EimdlMosaicCodec`) + 4-stripe rANS novelty (`FLAG_ANS4_NOVELTY`).
   - **Coupled Phred Quality Stream (`ANS_MODE_RANS4_QUAL = 3`, `ANS_MODE_RANS4_QUAL_BWT = 4`):**
     - `ANS_MODE_RANS4_QUAL` (`3`): hierarchical `(same, coarse[9], fine[8], nibble[16])` residual coding conditioned on Mosaic variant/literal `mismatch_mask`, DNA `homo_bins`, `prev_qbin`, and read cycle phase — regenerated deterministically from the DNA Mosaic tokens at zero mask overhead.
     - `ANS_MODE_RANS4_QUAL_BWT` (`4`, v3.1): pangenome-stratified (`mismatch_mask == 0` matched partition vs. `mismatch_mask == 1` variant/novelty partition at zero mask overhead) linear-time `SA-IS` Burrows-Wheeler Transform per `2,000,000`-byte block -> active-alphabet `MTF` -> bijective `RUNA`/`RUNB` zero-run collapsing -> 3-context 4-stripe interleaved Byte-rANS, strictly gated by Law IV to pick `min(raw, mode 3, mode 4)`.
4. **Optional 6-Frame CDS Proteome Summary (`FLAG3_HAS_CDS = 0x20`):**
   - Appends `(orfCount, longestOrfAa, totalOrfAa, gcPermille, proteomeSha256[32])` for genomic/transcriptomic validation without disturbing byte-exact reconstruction.
5. **Compressed-Domain Algebra (`HGS2` & `HGS3`):**
   - `analyze` computes Kimura `Ti/Tv` / $\kappa$ proxy, codon-wobble mod-3 periodicity, indel spectrum, donor recombination switches, and donor base shares directly from compressed tokens.
   - `variants` emits strand-normalized, 1-base anchored **VCF v4.2** calls directly without target FASTA reconstruction.
   - `kmer` executes exact forward + reverse-complement k-mer queries across conserved panel spans (zero target buffer materialization), SNP-patched local spans, and span boundary stitches.
