# Third-party codecs

GraSS and libbsc are fetched at fixed commits and built inside the locked
image. The checkouts, the `bsc` binary, and the class files are gitignored;
`fetch_vendor.ps1` recreates them.

| Tool | Commit | Build |
| --- | --- | --- |
| GraSS | `6d93700f5a2de8841f7cf2c20064413c80df2c4a` | `javac -d build src/*.java` |
| libbsc | `1dc676540cb40a41562af6f2db04bd1adce97f90` | g++ 15.2, `-O2`, manual object list from the GraSS README |

`bsc` is version 3.4.0 (21 October 2026), the current upstream at the time of
this measurement. GraSS hard-codes `./bsc e FinalTar.tar FinalBsc.bsc -e2`,
so it uses a 25 MB block, LZP preprocessing (`-H15`, `-M72`), and the adaptive
QLFC entropy coder. The binary is a musl build; `ldd` is not applicable.

## Why the flat input works

`input.seq` has no `>` header and no line breaks. GraSS treats it as one raw
sequence: its `seqExtraction` writes `F1.SSG` beginning with a literal `0`,
and `Decompress.finalDecomp` takes the raw path that writes one line to
`Final.SSG`. On a 2,000,000-base pilot this reproduced the input byte for
byte. In other words, the codec consumes the same 120,108,223-byte
`input.seq` as every other profile; no FASTA wrapper is inserted.

`Compress.bscCompress` shells out to `tar` and `./bsc`, so the Java process
must run with a working directory that contains `bsc`. The harness copies
`vendor/bin/bsc` into each operation directory for that reason.

## Measured result

On the full 120,108,223-byte pangenome (`runs/grass-1`) the unmodified
pipeline produced 24,777,592 B (1.650 bits/base, 4.85x) and restored it
byte-exactly. Compression took 224.5 s at 1,254.5 MiB peak RSS; restore took
165.5 s at 1,019.3 MiB. The archive is 6.50x the HyperGraSS archive and 6.36x
JARVIS3 level 7. See [`docs/SPECIALIZED.md`](../../../docs/SPECIALIZED.md)
for the full comparison.

GraSS's `tar` call embeds the F1.SSG/F2.SSG modification times, so the archive
bytes and size are not fully deterministic. Two 2 MB reruns produced identical
F1/F2 streams (same sizes and SHA-256) but archives of 493,330 and 493,326 B.
The 24,777,592 B figure is the size of the committed `grass-1` run.

## Rebuild

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File benchmarks/specialized/vendor/fetch_vendor.ps1
```
