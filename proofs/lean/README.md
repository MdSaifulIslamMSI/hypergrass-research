# EIMDL Formal Proofs in Lean 4 (`proofs/lean/`)

Zero-dependency, zero-`sorry` **Lean 4** (`leanprover/lean4:v4.36.0`) formal verification of the core algebraic, combinatorial, DAG, safety-ceiling, and adaptive arithmetic invariants underlying **EIMDL** ([`docs/COMPRESSION-LAW.md`](../../docs/COMPRESSION-LAW.md)).

## Verified Modules & Key Theorems

| Module | EIMDL Law | Core Formal Theorems Verified (`0 sorry`) |
|---|---|---|
| [`EimdlProofs/Law1KleinV4.lean`](EimdlProofs/Law1KleinV4.lean) | **Law I** | `diff_strand_invariant` ($(r \oplus 3) \oplus (t \oplus 3) = r \oplus t$), `transition_iff_delta2`, `revComp_involutive` ($\forall s : \text{List Base}$), `canonicalHash_revComp_invariant`, plus SAT-certified **`BitVec 64` SWAR** kernel theorems `swarReverseComplementWord32_involutive` and `swar_v4_diff_strand_invariant` mirroring `Packed2BitSequence.reverseComplementWord32` across all $2^{64}$ register words and all 32 lanes (`swar_lane_spec`). |
| [`EimdlProofs/Law2CodonWobble.lean`](EimdlProofs/Law2CodonWobble.lean) | **Law II** | Unbounded bijection theorems `wobble_left_inverse`, `wobble_right_inverse` ($\mathbb{N} \leftrightarrow \mathbb{N} \times \text{Fin } 3$, strengthening Java's bounded `0..2_000_000` loop check to $\forall g \in \mathbb{N}$), `bioPhase_involutive`, and `antisenseStop_iff_revComp_senseStop` across all 64 triplets. |
| [`EimdlProofs/Law3ModeCDag.lean`](EimdlProofs/Law3ModeCDag.lean) | **Law III** | Strict topological descent (`ancestor_strictly_smaller`), guaranteed Mode C acyclic pangenome reference graph (`causal_dag_acyclic`), bounded chain depth (`chainDepth_le_rank`), and MDL panel extension bound. |
| [`EimdlProofs/Law4UniversalBound.lean`](EimdlProofs/Law4UniversalBound.lean) | **Laws IV & VII** | Universal non-expansion ceiling `eimdl_universal_ceiling` ($L \le \text{hdr} + \lceil n/4 \rceil$), exact `500,000 bp -> 125,075 B` cryptographic noise bound (`eimdl_500k_crypto_bound`), and disjoint Stage-1/Stage-2 memory reduction (`disjointPeak_strict_saving`). |
| [`EimdlProofs/Law6TwoAdicTaper.lean`](EimdlProofs/Law6TwoAdicTaper.lean) | **Law VI** | Fixed-point separation (`classicStep0_fixed_point_iff`: $31$ vs `taper20_fixed_point_iff`: $11$), closed invariant band `[11, 4085]` across arbitrary binary streams (`taper20_stream_invariant`, proving zero-width division collapse impossibility), strict monotone Lyapunov convergence, and reciprocal table soundness (`recip_mul_le_exact`). |

## Build & Check Locally / in CI

```bash
cd proofs/lean
lake build
lake exe eimdl_proofs_check
```

Verified automatically on every push by the **`Formal Lean 4 Proofs (Laws I-VII)`** job in [`.github/workflows/ci.yml`](../../.github/workflows/ci.yml) and audited for zero `sorry` / zero unproved axioms by [`paper/audit_paper_numbers.py`](../../paper/audit_paper_numbers.py).
