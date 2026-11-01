# Independent pangenome replication

## Question

The 41.8% win over xz on the *E. coli* pangenome was measured after the
16 MiB / 8 MiB profile was selected on that same pangenome. That is a tuning
result on one input. Does it survive on pangenomes that played no part in
choosing the profile?

## Protocol

- Two additional NCBI RefSeq complete-genome pangenomes:
  *Salmonella enterica* (25 records, 120,527,038 bp, 122,251,418-byte FASTA)
  and *Klebsiella pneumoniae* (23 records, 122,700,601 bp, 124,455,339-byte
  FASTA).
- Both use the same selection rule as the *E. coli* set: NCBI relevance order,
  distinct accessions, records between 3 and 7.5 Mbp, stop once the target is
  reached. Accessions, lengths, and SHA-256 hashes are in
  `benchmarks/independent/data/<name>/manifest.json`.
- The locked Linux image, 2 CPUs, 2 GiB container limit, 1 GiB Java heap,
  fresh process per operation, and SHA-256 round-trip check are the same as the
  published scale test.
- Profiles: the published default (`hypergrass`), the scale-tuned
  `--pangenome` profile (16 MiB bank / 8 MiB blocks), `gzip -9`, and `xz -9`.
  No parameter was changed after seeing these results.
- Each HyperGraSS archive is additionally queried for records 0, the middle
  record, and the last record; every slice must match the source bytes.
- One full run per dataset. Times are indicative, as in the original scale
  test; archive sizes and round trips are the stable measurements.

<!-- INDEPENDENT-TABLES:START -->
| Dataset | Profile | Archive | Ratio | Bits/base | Wall | Peak RSS | Exact | Δ vs xz |
| --- | --- | ---: | ---: | ---: | ---: | ---: | :---: | ---: |
| klebsiella | hypergrass | 8,571,498 B | 14.52x | 0.559 | 65.0 s | 1147 MiB | yes | +33.12% |
| klebsiella | hypergrass-pangenome | 3,464,423 B | 35.92x | 0.226 | 10.6 s | 912 MiB | yes | -46.19% |
| klebsiella | gzip-9 | 36,388,919 B | 3.42x | 2.373 | 201.1 s | 10 MiB | yes | +465.16% |
| klebsiella | xz-9 | 6,438,744 B | 19.33x | 0.420 | 521.3 s | 736 MiB | yes | — |
| salmonella | hypergrass | 7,317,099 B | 16.71x | 0.486 | 60.9 s | 922 MiB | yes | +19.56% |
| salmonella | hypergrass-pangenome | 3,443,600 B | 35.50x | 0.229 | 28.6 s | 732 MiB | yes | -43.73% |
| salmonella | gzip-9 | 35,969,250 B | 3.40x | 2.387 | 97.6 s | 10 MiB | yes | +487.71% |
| salmonella | xz-9 | 6,120,260 B | 19.97x | 0.406 | 184.6 s | 734 MiB | yes | — |

Tuned profile vs xz: median -44.96% across 2 independent pangenomes.
<!-- INDEPENDENT-TABLES:END -->

## Findings

The tuned profile is smaller than xz on both independent pangenomes:

- *Salmonella enterica*: 3,443,600 B vs 6,120,260 B, **43.7% smaller**.
- *Klebsiella pneumoniae*: 3,464,423 B vs 6,438,744 B, **46.2% smaller**.

The original *E. coli* pangenome was 41.8% smaller (4,885,851 B vs
8,394,900 B). That measurement is not part of this independent count because
the profile was selected on it. The two independent margins have a median of
44.96% and a range of 43.73–46.19%; the direction and rough magnitude of the
original result therefore replicate.

The default profile loses to xz on both new pangenomes, as it did on *E. coli*:
+19.56% on Salmonella and +33.12% on Klebsiella, versus +9.78% on *E. coli*.
The failure of the 4 MiB bank is also not an *E. coli* accident.

Compression time moved in the same direction as size on these runs: the tuned
profile compressed Salmonella in 28.6 s versus 184.6 s for xz, and Klebsiella
in 10.6 s versus 521.3 s. Restore time was mixed: 9.5 s versus 6.9 s on
Salmonella, and 9.9 s versus 11.4 s on Klebsiella. Do not read these times too
precisely. The gzip baseline took 97.6 s on Salmonella and 201.1 s on
Klebsiella for almost the same input size, which shows how much host variance
is present. Archive sizes and round trips are the stable measurements.

Two details matter for interpretation:

1. **This is still a small, convenience sample.** Three pangenomes from one
   family of Gram-negative bacteria, each built from an NCBI relevance-ordered
   selection, do not represent all pangenomes. The result supports
   "replicates on two additional *Enterobacteriaceae* pangenomes", not
   "wins on pangenomes in general".
2. **The margin depends on the input, not just the profile.** The tuned
   profile's archive size, compression time, and restore time differ between
   the three inputs. The comparison that is stable is the direction: tuned
   < xz < default on every pangenome measured here, with the default losing to
   xz.

## Limits

- Two independent datasets, one run each, one host.
- The datasets were fetched on 2026-11-01. NCBI relevance order can change;
  the manifests lock the exact accessions and hashes that were measured, not a
  permanent snapshot of the query.
- The FASTA inputs include headers and wrapping; the corpus comparison does
  not. These numbers are comparable to the published 121.8 MB scale test, not
  to the sequence-only corpus.
- Query probes are correctness checks, not a query-latency workload. The
  end-to-end workload is measured separately on the *E. coli* pangenome.

## Reproduce

```sh
python -B benchmarks/independent/fetch.py --name salmonella --organism "Salmonella enterica"
python -B benchmarks/independent/fetch.py --name klebsiella --organism "Klebsiella pneumoniae"
python -B benchmarks/independent/run.py --dataset salmonella --run-id salmonella-scale-1
python -B benchmarks/independent/run.py --dataset klebsiella --run-id klebsiella-scale-1
python -B benchmarks/independent/summarize.py
```

FASTAs are gitignored; the manifests and `report.json` files are the published
evidence.
