# HyperGraSS Archive Format (`.hgs`), version 1

Status: implementation-defined format, documented from the HyperGraSS 1.0.0
encoder and decoder. It is not yet an external standard. A decoder must reject
malformed input; the encoder emits only the canonical choices described here.

## 1. Scope

An `.hgs` file stores an arbitrary byte stream losslessly. When the input is
FASTA, the encoder additionally factorizes headers, line layout, soft-masking,
scaffold `N` runs, and non-ACGT symbols into separate channels, then compresses
the DNA stream with a syncmer-anchored model. The decoder reconstructs the
original bytes exactly, including headers, wrapping, case, and line endings.

The format has two physical layouts:

| Layout | Used for | Header |
| --- | --- | --- |
| Compact single-block | Small inputs chosen by the encoder | 13-byte header + one payload |
| Multi-block | General inputs | Header, optional anchor bank, blocks, footer, 12-byte trailer |

Both layouts start with the same 4-byte magic `48 47 53 01` (`"HGS"` + version
byte `0x01`). Readers must check the magic and the layout flag before parsing
anything else.

## 2. Primitive encodings

- **`uvarint`** — unsigned LEB128, 7 bits per byte, least-significant group
  first, high bit set on every byte except the last. Values are limited to the
  range of the field's declared type.
- **`svarint`** — signed ZigZag followed by `uvarint`:
  `zz = (value << 1) ^ (value >> 31)` for 32-bit values.
- **`intLE` / `longLE`** — fixed-width little-endian 32-bit / 64-bit integers.
- **`bytes[n]`** — `n` raw bytes with no length prefix unless stated.
- **`aux`** — a length-delimited auxiliary channel:

  | First byte | Meaning |
  | --- | --- |
  | `0` | Empty. No further bytes. |
  | `1` | Raw. Followed by `uvarint` length and that many bytes. |
  | `2` | Deflate. Followed by `uvarint` uncompressed length and a raw DEFLATE stream (`nowrap=true`, level 9). |

  The encoder chooses raw when the payload is at most 12 bytes, otherwise it
  picks whichever of raw or Deflate is strictly smaller. A decoder must accept
  raw even when Deflate would also have been valid.

## 3. Common header

```
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

`COMPACT_SINGLE` and `HAS_ANCHOR` are mutually exclusive. `COMPACT_SINGLE`
implies no anchor bank. `LEAN_MODE` affects the payload encoding choices, not
the framing.

## 4. Compact single-block layout

```
offset  size  field
0       4     magic
4       1     flags = 0x04
5       4     rawLen, intLE
9       4     crc32, intLE
13      ...   payload
```

`rawLen` is the exact decompressed byte length. `crc32` is the standard
CRC-32 (IEEE 802.3, as in `java.util.zip.CRC32`) of the decompressed bytes.
`payload` is a block payload as defined in section 6. The archive ends after
the payload; there is no footer or trailer.

## 5. Multi-block layout

```
offset  size  field
0       4     magic
4       1     flags
5       ...   optional anchor bank (present iff HAS_ANCHOR)
...     ...   block payloads, concatenated
...     ...   footer body (aux-encoded)
...     8     footerOffset, longLE
...     4     FOOTER_MAGIC = 0x31464748 ("HGF1" in little-endian bytes)
```

The last 12 bytes are always the trailer. `footerOffset` is the absolute file
offset of the footer body. The footer body is `aux`-encoded; its decoded bytes
are the directory. Block offsets in the directory are absolute.

### 5.1 Anchor bank

Present only when `HAS_ANCHOR` is set:

```
fast2Bit      byte        0 = context-mixer-coded pages, 1 = 2-bit packed pages
anchorLen     uvarint     total anchor bases
numPages      uvarint     number of pages
primarySeedLen uvarint
numDonors     uvarint
donorGaps     aux         per donor: uvarint length, then per element
                          uvarint(pos - prevPos - 1)
donorDiffs    aux         one byte per element, low 2 bits used
pageBaseLen[] uvarint     one per page, in page order
pageEncLen[]  uvarint     one per page, in page order
pagePayload[] bytes       one per page, in page order
```

`ANCHOR_PAGE_SIZE` is 65,536 bases. The final page may be shorter. When
`fast2Bit == 1`, each page payload is an `aux` stream containing 2-bit packed
bases (A=0, C=1, G=2, T/U=3), most-significant base first within each byte.
When `fast2Bit == 0`, each page payload is a `QuantumContextRepeatMixer` stream
as described in section 6.3. Global donor entries encode a position and a
2-bit difference; they are an internal reuse catalog, not part of the public
API.

### 5.2 Footer body

The decoded footer body is:

```
numBlocks      uvarint
totalRawBytes  uvarint (64-bit)
totalRecords   uvarint
crc32          intLE
entries[numBlocks]:
  offsetDelta  uvarint (64-bit)
  compressedLen uvarint
  firstRecordDelta uvarint
  recordCount  uvarint
  firstRecordBaseOffset uvarint (64-bit)
  blockSeqChars uvarint
```

For entry 0, `offsetDelta` is the absolute block offset. For later entries it
is `blockOffset - previousBlockEnd`, where `previousBlockEnd` is the previous
entry's `blockOffset + compressedLen`. `firstRecordDelta` is relative to the
previous entry's `firstRecordIdx`; entry 0 is relative to 0. This keeps the
directory small and makes block order explicit.

`firstRecordBaseOffset` is the byte offset of the block's first record within
the decompressed sequence stream. It is what makes direct record and range
queries possible without decompressing earlier blocks.

`blockSeqChars` is the number of sequence characters contributed by the block;
it is used to map a record-relative query range to the block's DNA stream.

## 6. Block payload

A block payload is:

```
auxMetadata    section 6.1
dnaPayload     section 6.2
```

### 6.1 Auxiliary metadata channels

```
flags          byte
dnaLen         uvarint
totalSeqChars  uvarint
numSeqLines    uvarint
```

Flag bits: `0x01` RNA mode, `0x02` headers present, `0x04` lowercase intervals
present, `0x08` `N` runs present, `0x10` exception symbols present.

Channel 1 — headers, present iff `0x02`:

```
numHeaders     uvarint
headerBlock    aux
```

`headerBlock` contains, for each header: `uvarint seqLineDelta`; then, for each
header in order, `uvarint sharedPrefixLen`, `uvarint suffixLen`, and the suffix
bytes. Headers are front-coded against the previous header. Header bytes
include the leading `>` and exclude the line ending.

Channel 2 — line layout (always present):

```
lineLengths    aux    RLE of sequence line widths
lineEols       aux    RLE of per-line EOL codes
```

EOL codes: `0` = none, `1` = LF, `2` = CRLF. The RLE streams contain
`(value, runLength)` pairs of `uvarint`s; see `writeRleInts` / `writeRleBytes`
in `FastaStreamFactorizer`.

Channel 3 — lowercase intervals, present iff `0x04`:

```
numIntervals   uvarint
minLength      uvarint
gaps           aux    uvarint gap to each interval
lengths        aux    uvarint (length - minLength)
```

Gaps are measured in uppercase bases between lowercase intervals. Case is
restored from these intervals, not stored per base.

Channel 4 — `N` runs, present iff `0x08`:

```
numRuns        uvarint
runs           aux    per run: uvarint ACGT gap, uvarint run length
```

Channel 5 — exception symbols, present iff `0x10`:

```
numExceptions  uvarint
gaps           aux    uvarint gap before each exception
symbols        aux    one raw byte per exception
```

Exception symbols are bytes that are neither A/C/G/T/U/N nor lowercase
variants. They are stored verbatim and reinserted at the recorded gaps.

The DNA stream contains only canonical base codes. `U` is folded to `T` in RNA
mode and restored on decode; lowercase is reconstructed from channel 3.

### 6.2 DNA payload

The first byte selects the DNA encoding:

| Mode | Layout |
| --- | --- |
| `0` | Empty DNA stream. No further bytes. |
| `1` | `uvarint rcLen`, then a `QuantumContextRepeatMixer` stream of `rcLen` bytes. |
| `2` | Syncmer matches with per-match literals and SNP residuals. |
| `3` | Anchor-coordinate haplotype matrix with donor reuse. |

Modes 2 and 3 are internal to `CanonicalSyncmerEngine`. Their subfields use
`aux` streams for literal lengths, strand bits, match geometry, SNP counts,
SNP gaps, and SNP differences. A conforming decoder should use the reference
implementation rather than reimplement these modes from this summary.

### 6.3 `QuantumContextRepeatMixer` stream

Mode 1 and non-2-bit anchor pages contain a self-delimiting context-mixer
stream. The stream is stored as `uvarint length` followed by that many bytes
in mode 1; anchor pages store the bytes inside an `aux` wrapper. The mixer is
a byte-oriented adaptive model; its internal state is reset per stream. The
decoder must consume exactly the declared number of output bytes and reject
truncated or overlong input.

## 7. Compression and verification rules

- The encoder computes CRC-32 over the exact decompressed bytes and stores it
  in the header (compact layout) or footer (multi-block layout).
- The decoder must verify the CRC before returning data.
- The decoder must reject unknown mode bits, impossible lengths, truncated
  varints, and directory entries that point outside the file.
- A round trip is successful only when the restored byte sequence is identical
  to the input, not merely when it decodes without error.

## 8. Query semantics

`HyperGraSS query <archive.hgs> <recordIndex> <start> <end>` returns a slice of
one record. `start` and `end` are zero-based, half-open, and relative to the
record's sequence characters; `end = -1` means the record end. The query reads
the footer, locates the block containing the record via
`firstRecordBaseOffset`, decodes only that block (plus any anchor pages it
needs), and writes the requested bases to stdout followed by a newline. Query
output is therefore not a FASTA record; it is the raw base slice.

## 9. Determinism and compatibility

- For a fixed input and fixed options, the encoder is deterministic: the same
  bytes and the same archive size are produced.
- The `--optimize` flag performs a bounded extra search and may choose a
  different, smaller archive; it does not change the format.
- `--lean` changes encoder memory behavior and payload choices but not the
  framing.
- Version byte `0x01` identifies this layout. A future incompatible layout
  must use a different magic/version byte; readers must not guess.

## 10. Reference implementation

- `src/hypergrass/HyperGraSS.java` — container, footer, CLI, query.
- `src/hypergrass/FastaStreamFactorizer.java` — six-channel metadata layout.
- `src/hypergrass/CanonicalSyncmerEngine.java` — DNA modes 0–3, anchor bank.
- `src/hypergrass/VarIntCodec.java` — varints, 2-bit packing, `aux` streams.
- `src/hypergrass/QuantumContextRepeatMixer.java` — context-mixer stream.
