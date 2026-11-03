# Bounded Memory Sweep (2026-11-03)

This note records the first measurement that separates the codec's live set
from the JVM's native overhead. It replaces the earlier single observation that
the 32 MiB profile needed "a 256 MiB heap."

## Input and profile

| Item | Value |
|---|---|
| Input | 120,108,223 bases, SHA-256 `51cb3ecb2e9f532833caa0b60345f37d0342474e9f4a239801486f9bc9219acb` |
| Image | `hypergrass-corpus:jdk27`, `sha256:580e58b6ee16d43f4dd99a2871b1c920548b7c6e9b0c695a01f521c17701a4ae` |
| Profile | `--pangenome --bounded --single-thread --anchor-single-page --seed-step 8 --block-mib 32 --anchor-mib 9` |
| JVM | `SerialGC`, 32 MiB code cache, 64 MiB metaspace cap, 2 compilers, 512 KiB stacks |
| Archive | 3,824,166 B, SHA-256 `c7400109056a52a610f78960eb19f36f8ab04a0ceb01e8b6384fa02125acd617` |

The profile is the only bounded candidate that passes the 2% size gate: it is
0.36% above the unbounded 3,810,409 B and 1.85% below JARVIS3 l7.

## Natural heap sweep

Each row is one observation under a permissive 512 MiB cgroup cap. A row is a
success only when both compress and restore complete and the independent
verifier reports byte equality.

| Heap | Archive | Compress wall | Compress peak RSS | Restore wall | Restore peak RSS | Round trip |
|---:|---:|---:|---:|---:|---:|:---:|
| 208m | 3,824,166 B | 143.9 s | 362.1 MiB | OOM at 71.6 s | 255.4 MiB | no |
| 224m | 3,824,166 B | 47.2 s | 374.8 MiB | OOM at 21.0 s | 262.7 MiB | no |
| 232m | 3,824,166 B | 52.3 s | 375.4 MiB | 36.6 s | 286.6 MiB | yes |
| 240m | 3,824,166 B | 53.5 s | 391.0 MiB | 29.9 s | 298.5 MiB | yes |

The encoder live set fits a 208 MiB heap, though that heap pays 143.9 s.
The decoder needs more than 224 MiB and no more than 232 MiB; the exact boundary
is between those values. The smallest completing configuration still peaks
375.4 MiB compressing and 286.6 MiB restoring.

The earlier 256 MiB-heap observation recorded 339.1 / 302.7 MiB. That sample is
not a ceiling: the same archive under the 232 and 240 MiB heaps peaked at 375.4
and 391.0 MiB compressing. Peak RSS moves with GC scheduling even when the
archive bytes are identical, so only the worst observation is a gate value.

## Native-memory attribution

`-XX:NativeMemoryTracking=summary -XX:+PrintNMTStatistics` at a 240 MiB heap:

| Process | Heap | Committed at exit | Peak RSS | Dominant transient |
|---|---:|---:|---:|---|
| Compress | 240 MiB | 262.6 MiB total, 22.6 MiB outside heap | 389.3 MiB | `malloc` peak 109.0 MiB; compiler arena 56.0 MiB; `Other` 50.5 MiB in 2 allocations; shared class space 13.9 MiB |
| Restore | 240 MiB | 271.8 MiB total, 31.8 MiB outside heap | 298.5 MiB | `malloc` peak 21.0 MiB; compiler arena 16.9 MiB; arena chunk 17.4 MiB; `Other` 2.2 MiB; shared class space 13.9 MiB |

The encoder's exit-time committed memory is only 262.6 MiB, but transient
compiler and malloc arenas push peak RSS to 389.3 MiB. Those arenas are JVM
overhead, not codec data. The decoder's overhead is smaller but still about
32 MiB at exit.

## Strict 256 MiB cap

With `--memory-swap` equal to `--memory`, no tested JVM heap completes:

| Heap | Result | Wall | Wrapper RSS | cgroup peak | `oom_kill` |
|---:|---|---:|---:|---:|---:|
| 224m | OOM-killed | 6.4 s | 247.1 MiB | 256.0 MiB | 1 |
| 240m | OOM-killed | 5.9 s | 247.8 MiB | 256.0 MiB | 1 |
| 248m | OOM-killed | 5.8 s | 247.9 MiB | 256.0 MiB | 1 |

The cap kills the process around 6 s, before the natural peak. A larger heap
does not help: the JVM commits native memory on top of the heap, and the
transient compressor arenas alone exceed the remaining budget. No JVM heap
setting in this range is a 256 MiB mode.

## Native runtime (GraalVM 25 musl)

The native build keeps the codec and archive format unchanged and removes the
JVM runtime tax. `build-native.ps1` compiles the sources with `javac --release
25` and links them with `native-image --static --libc=musl --no-fallback
--gc=serial -O3 -march=compatibility` from the pinned
`ghcr.io/graalvm/native-image-community:25-muslib` image. The binary is
13,698,160 B, SHA-256
`508f4d53d5d81563f408dd59f8ffdbb79ae60fd86ab4ff1eb14dfdbc3530b657`.

Two opt-in memory reductions make the profile fit:

- `FastaStreamFactorizer.reconstructTo` writes headers, sequence, and line
  endings straight to the output stream, so the decoder no longer allocates a
  third block-sized output array.
- `--split-tail-block` cuts a huge final record at the target block size. On
  this input it changes the layout from 32 + 32 + 56 MiB to four blocks:
  32 + 32 + 32 + 18.5 MiB. The default partition is unchanged.

The natural observation used a 192 MiB heap under a permissive 512 MiB cap:

| Direction | Archive | Wall | Peak RSS | cgroup peak | Exact |
|---|---:|---:|---:|---:|:---:|
| Compress | 3,864,464 B | 143.3 s | 240.1 MiB | 250.8 MiB | yes |
| Restore | same | 29.0 s | 225.4 MiB | 235.3 MiB | yes |

The archive is `8f351911f318da52c368196a47d63dc4cf0ced3791b749e0ca1bd9d7e13c3705`,
1.42% above the unbounded 3,810,409 B and inside the 2% allowance. The
independent verifier compared all 120,108,223 bytes. That pilot used a
512 MiB cap; the strict 256 MiB gate below is the binding measurement.

## Strict 256 MiB gate (native)

`native-gate-3` interleaves three JARVIS3 l7 runs with three native bounded
runs in the same `hypergrass-corpus:jdk27` image on the identical input. The
native candidate runs under `--memory=256m --memory-swap=256m`. JARVIS3 keeps
the usual 2 GiB limit because it cannot finish under 256 MiB; that asymmetry is
what the bounded requirement tests, and it is recorded rather than hidden.

| Profile | Archive | Compress median | Worst compress RSS | Worst restore RSS | Worst cgroup peak | Exact |
|---|---:|---:|---:|---:|---:|:---:|
| JARVIS3 l7 (2 GiB limit) | 3,896,361 B | 58.3 s | 286.0 MiB | 285.9 MiB | 294.8 MiB | 3/3 |
| Native bounded (256 MiB limit) | 3,864,464 B | 54.8 s | 240.1 MiB | 225.5 MiB | 251.1 MiB | 3/3 |

The native archive is 1.42% above the unbounded 3,810,409 B (the 2% limit is
3,886,617 B) and 0.82% below JARVIS3. Every gate check passes on the worst of
the three interleaved repeats:

| Check | Value | Limit | Verdict |
|---|---:|---:|:---:|
| status | success / success | success | pass |
| compress peak RSS | 240.1 MiB | 256.0 MiB | pass |
| restore peak RSS | 225.5 MiB | 256.0 MiB | pass |
| size ratio vs unbounded | 1.014 | 1.020 | pass |
| wall ratio vs JARVIS3 median | 0.940 | 2.000 | pass |

All three archives are byte-identical at
`8f351911f318da52c368196a47d63dc4cf0ced3791b749e0ca1bd9d7e13c3705`; each of
the three restores was compared over all 120,108,223 bytes by the separate
verifier process. The cgroup compress headroom is only **4.9 MiB** (251.1 of
256.0 MiB) even though the `ru_maxrss` headroom is 15.9 MiB, so the margin is
real but thin. The native binary is 13,698,160 B, SHA-256
`508f4d53d5d81563f408dd59f8ffdbb79ae60fd86ab4ff1eb14dfdbc3530b657`.

Archive bytes on the default path are unchanged: recompressing the 40.5 MB
campy-high collection with the current `bin/` at a 160 MiB heap and the original
bounded flags produced `1,125,629 B` with SHA-256
`14DD7E9413C10B8414E90EE831666B420D81A9A3BA2A9F8DBBC1BED149E1C5D6`, an exact
match to the published pre-change archive. `--split-tail-block` is opt-in, so
existing archives keep their bytes.

## Consequence

The live data fits 256 MiB: the encoder needs at most 208 MiB of heap and the
decoder 225-232 MiB. The JVM does not: its runtime adds 143-154 MiB on the
encoder and about 32-55 MiB on the decoder, dominated by JIT compiler arenas,
malloc transients, and shared class metadata. The native image removes that
tax without changing the archive format, and the native profile now completes
3/3 under a strict no-swap 256 MiB cap.

The remaining costs are explicit: +1.42% archive versus the unbounded profile,
4.9 MiB of cgroup headroom, one host, and no cross-taxon native run. The JVM
bounded mode remains an *E. coli*-calibrated prototype; the native gate does
not retroactively make the JVM path a 256 MiB mode.

## Evidence

- `benchmarks/specialized/runs/heap-natural-208/report.json`
- `benchmarks/specialized/runs/heap-natural-224/report.json`
- `benchmarks/specialized/runs/heap-natural-232/report.json`
- `benchmarks/specialized/runs/heap-natural-240/report.json`
- `benchmarks/specialized/runs/heap-nmt-240b/report.json`
- `benchmarks/specialized/runs/heap-nmt-240b/ops/hypergrass-bounded-wide-r0/{compress,decompress}/stdout.log`
- `benchmarks/specialized/runs/heap-sweep-224/report.json`
- `benchmarks/specialized/runs/heap-sweep-240/report.json`
- `benchmarks/specialized/runs/heap-sweep-248/report.json`
- `benchmarks/specialized/runs/native-smoke-5m-split/report.json`
- `benchmarks/specialized/runs/native-natural-192-split/report.json`
- `benchmarks/specialized/runs/native-gate-3/report.json`
- `benchmarks/specialized/runs/native-gate-3/ops/hypergrass-bounded-native-r{0,1,2}/{compress,decompress}/metrics.json`
