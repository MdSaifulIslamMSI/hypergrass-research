# Specialized DNA compressors at pangenome scale

## Question

The frozen corpus compares HyperGraSS with JARVIS3 and GeCo3 on individual
genomes, where the margins are narrow. Does either specialized DNA compressor
close the gap on a 121.8 MB multi-strain pangenome, where most of the redundant
information is between genomes rather than within one genome?

## Protocol

- Input: the same 24-record *E. coli* pangenome as the scale test. For the
  specialized tools the FASTA is converted to one flat sequence: headers
  removed, line breaks removed, letters uppercased, and the 2,357 non-A/C/G/T
  symbols (0.002%) replaced with `A`. This is the raw-sequence contract JARVIS3
  requires; GeCo3 also strips headers and compresses only A/C/G/T.
- The normalization is deliberately conservative for HyperGraSS: the
  specialized tools receive slightly less information (no headers and no rare
  ambiguity symbols) than the 4,885,851-byte multi-FASTA archive.
- Same locked Linux image, 2 CPUs, 2 GiB container limit, fresh process per
  operation.
- Profiles: JARVIS3 level 7, GeCo3 levels 3 and 6, and `xz -9` on the identical
  sequence-only bytes. GeCo3 level 9 and JARVIS3 level 26 are excluded at this
  scale: in the frozen corpus their 4.6 MB *E. coli* cells reached about
  2.14 GB RSS and the 180 s cap, so they are not feasible on 120 MB under the
  same 2 GiB container. Those are infeasibility findings, not size results.
- Each decompressed sequence is hashed with SHA-256 and must equal the input
  sequence exactly.

<!-- SPECIALIZED-TABLES:START -->
<!-- SPECIALIZED-TABLES:END -->

## Findings

The hypothesis was wrong. JARVIS3 level 7 — a general DNA compressor with a
context-mixing model — compressed the identical 120.1 Mbp sequence to
**3,896,361 bytes (0.260 bits/base)**, **20.3% smaller than the tuned
HyperGraSS archive**, and restored the sequence exactly. The table gives the
sizes and times for both tools and for GeCo3 levels 3 and 6 on the same bytes.

The trade is speed, not size. HyperGraSS compressed in 30.6 s versus 95.4 s
(3.1x faster) and restored in 16.4 s versus 84.2 s (5.1x faster), at 419 MiB
versus 286 MiB peak restore RSS. The tuned HyperGraSS archive is 25.4% larger
than JARVIS3's. So the 41.8–46.2% margin over xz is real, but xz is a
general-purpose baseline. Against a specialized DNA compressor the tuned
profile loses on archive size at this scale.

This does **not** overturn the frozen-corpus result. On single genomes,
JARVIS3 and GeCo3 are close to HyperGraSS, and on some cells they are better.
What it overturns is the stronger interpretation of the pangenome result: the
16 MiB anchor bank does not make HyperGraSS the smallest archive when many
near-identical genomes are stored together. A specialized context-mixing coder
with no explicit pangenome model finds at least as much redundancy.

JARVIS3 also has no indexed query interface. Its restore is a full 84.2 s
decompression, so it is not a competitor for the slice-query workload in
[`END-TO-END.md`](END-TO-END.md); it is a better choice for store-once,
restore-rarely archives, which is the workload where HyperGraSS's size claim
was strongest.

## Limits

- One pangenome, one run, one host.
- The sequence-only normalization removes headers, wrapping, case, and rare
  ambiguity symbols. The HyperGraSS archive includes those bytes. The
  asymmetry favors the specialized tools.
- JARVIS3 level 26 and GeCo3 level 9 are not measured at this scale, so this is
  not a comparison against either tool's maximum-compression setting.
- The specialized tools have no equivalent of the record-aware anchor bank or
  the page index; this measures their out-of-the-box sequence mode.
- Times are indicative; the original scale test already documents large
  wall-time variance on this host.

## Reproduce

```sh
python -B benchmarks/specialized/run.py \
  --run-id specialized-scale-1 \
  --profiles jarvis3-l7,geco3-l3,geco3-l6,xz-9 \
  --timeout 7200
python -B benchmarks/specialized/summarize.py --run specialized-scale-1
```

The sequence extract and archives are gitignored; `report.json` is the
published evidence.
