/-
  EimdlProofs.Law6TwoAdicTaper
  Formal machine-checked proofs for Law VIa & Law VIc of EIMDL (`docs/COMPRESSION-LAW.md`):
    1. Exact integer fixed-point comparison of classic shift-5 updates (`31/4096`)
       vs two-adic tapered adaptation (`tau = 20` -> `11/4096`)
    2. General closed-form fixed-point characterization for any `0 <= tau <= 31`:
       `taperStep0 tau p = p <-> p <= 31 - tau`
    3. Strict invariant band `[11, 4085]` under BOTH `b = 0` and `b = 1` updates
       over ANY finite bitstream (proving zero division / zero interval collapse impossibility)
    4. Strict monotone descent (`p >= 12`) and ascent (`p <= 4084`) convergence to attractors
    5. Law VIc fixed-point reciprocal table ceiling soundness (`(2^30) / d`)
-/
namespace EimdlProofs.Law6

/-- Probability scale `M = 2^12 = 4096` (`BitRangeEnc` / `BitRangeDec`). -/
def probScale : Nat := 4096

/-- Classic shift-5 update on bit `b = 0`: `p - (p >>> 5)`. -/
def classicStep0 (p : Nat) : Nat :=
  p - p / 32

/-- Two-adic tapered shift-5 update on bit `b = 0`: `p - ((p + tau) >>> 5)`. -/
def taperStep0 (tau p : Nat) : Nat :=
  p - (p + tau) / 32

/-- Symmetric two-adic tapered shift-5 update on bit `b = 1`: `p + ((4096 - p + tau) >>> 5)`. -/
def taperStep1 (tau p : Nat) : Nat :=
  p + (probScale - p + tau) / 32

/-- Classic update stalls prematurely at all `p <= 31` (`31 / 4096 = 0.757%`). -/
theorem classicStep0_fixed_point_iff (p : Nat) :
    classicStep0 p = p ↔ p ≤ 31 := by
  unfold classicStep0
  omega

/-- General theorem for any two-adic taper `tau <= 31`: exact fixed point set is `[0, 31 - tau]`. -/
theorem taperStep0_fixed_point_iff {tau : Nat} (ht : tau ≤ 31) (p : Nat) :
    taperStep0 tau p = p ↔ p ≤ 31 - tau := by
  unfold taperStep0
  omega

/-- EIMDL specialization (`tau = 20`): stalling threshold drops cleanly from `31` to `11` (`0.269%`). -/
theorem taper20_fixed_point_iff (p : Nat) :
    taperStep0 20 p = p ↔ p ≤ 11 :=
  taperStep0_fixed_point_iff (tau := 20) (by decide) p

/-- Symmetric upper fixed point on `b = 1` with `tau = 20`: threshold rises cleanly to `4085`. -/
theorem taper20_upper_fixed_point_iff {p : Nat} (hle : p ≤ 4096) :
    taperStep1 20 p = p ↔ 4085 ≤ p := by
  unfold taperStep1 probScale
  omega

/-- Core Arithmetic Range-Coder Safety Invariant:
    Once `p` lies in `[11, 4085]`, BOTH `b = 0` and `b = 1` tapered updates stay inside `[11, 4085]`,
    so `0 < p < 4096` holds unconditionally forever (no zero-width degeneracy). -/
theorem taper20_band_closed (p : Nat) (hlo : 11 ≤ p) (hhi : p ≤ 4085) :
    (11 ≤ taperStep0 20 p ∧ taperStep0 20 p ≤ 4085) ∧
    (11 ≤ taperStep1 20 p ∧ taperStep1 20 p ≤ 4085) := by
  unfold taperStep0 taperStep1 probScale
  omega

/-- Execute an arbitrary sequence of bits `bits : List Bool` through the adaptive two-adic coder. -/
def runTaperStream (tau : Nat) (p0 : Nat) (bits : List Bool) : Nat :=
  bits.foldl (fun p b => if b then taperStep1 tau p else taperStep0 tau p) p0

/-- Multi-step Range-Coder Safety across ANY finite binary message of arbitrary length! -/
theorem taper20_stream_invariant (p0 : Nat) (bits : List Bool)
    (hlo : 11 ≤ p0) (hhi : p0 ≤ 4085) :
    11 ≤ runTaperStream 20 p0 bits ∧ runTaperStream 20 p0 bits ≤ 4085 := by
  induction bits generalizing p0 with
  | nil => exact ⟨hlo, hhi⟩
  | cons b rest ih =>
    simp only [runTaperStream, List.foldl_cons]
    have hstep := taper20_band_closed p0 hlo hhi
    cases b
    · exact ih (taperStep0 20 p0) hstep.1.1 hstep.1.2
    · exact ih (taperStep1 20 p0) hstep.2.1 hstep.2.2

/-- Strictly monotone Lyapunov descent toward `11` on `b = 0` whenever `p >= 12`. -/
theorem taper20_strict_descent (p : Nat) (h : 12 ≤ p) : taperStep0 20 p < p := by
  unfold taperStep0
  omega

/-- Strictly monotone Lyapunov ascent toward `4085` on `b = 1` whenever `p <= 4084`. -/
theorem taper20_strict_ascent (p : Nat) (h : p ≤ 4084) : p < taperStep1 20 p := by
  unfold taperStep1 probScale
  omega

/-- Helper iterating `taperStep0 20` on an all-zeros stream. -/
def iterTaper20Zeros : Nat → Nat → Nat
  | 0, p => p
  | n + 1, p => iterTaper20Zeros n (taperStep0 20 p)

/-- Helper iterating `classicStep0` on an all-zeros stream. -/
def iterClassicZeros : Nat → Nat → Nat
  | 0, p => p
  | n + 1, p => iterClassicZeros n (classicStep0 p)

theorem iterTaper20Zeros_add (a b p : Nat) :
    iterTaper20Zeros (a + b) p = iterTaper20Zeros b (iterTaper20Zeros a p) := by
  induction a generalizing p with
  | zero => simp [iterTaper20Zeros]
  | succ a ih =>
    rw [Nat.succ_add, iterTaper20Zeros, iterTaper20Zeros, ih]

theorem iterClassicZeros_add (a b p : Nat) :
    iterClassicZeros (a + b) p = iterClassicZeros b (iterClassicZeros a p) := by
  induction a generalizing p with
  | zero => simp [iterClassicZeros]
  | succ a ih =>
    rw [Nat.succ_add, iterClassicZeros, iterClassicZeros, ih]

theorem iterTaper20Zeros_fixed_11 (n : Nat) : iterTaper20Zeros n 11 = 11 := by
  induction n with
  | zero => rfl
  | succ n ih => simp [iterTaper20Zeros, taperStep0, ih]

theorem iterClassicZeros_fixed_31 (n : Nat) : iterClassicZeros n 31 = 31 := by
  induction n with
  | zero => rfl
  | succ n ih => simp [iterClassicZeros, classicStep0, ih]

/-- Asymptotic theorem matching `VerifyCompressionLaws.java`: for ANY horizon `n >= 150`
    (in particular `n = 4096`), `taper20` hits `11/4096` whereas classic stalls at `31/4096`. -/
theorem java_loop_all_zeros_reaches_exact_fixed_point (n : Nat) (hn : 150 ≤ n) :
    iterTaper20Zeros n 2048 = 11 ∧ iterClassicZeros n 2048 = 31 := by
  have ht1 : iterTaper20Zeros 25 2048 = 924 := by decide
  have ht2 : iterTaper20Zeros 25 924 = 416 := by decide
  have ht3 : iterTaper20Zeros 25 416 = 185 := by decide
  have ht4 : iterTaper20Zeros 25 185 = 81 := by decide
  have ht5 : iterTaper20Zeros 25 81 = 36 := by decide
  have ht6 : iterTaper20Zeros 25 36 = 11 := by decide
  have hbase1 : iterTaper20Zeros 150 2048 = 11 := by
    have hs : 150 = 25 + 25 + 25 + 25 + 25 + 25 := rfl
    rw [hs, iterTaper20Zeros_add, iterTaper20Zeros_add, iterTaper20Zeros_add,
        iterTaper20Zeros_add, iterTaper20Zeros_add, ht1, ht2, ht3, ht4, ht5, ht6]
  have hc1 : iterClassicZeros 25 2048 = 932 := by decide
  have hc2 : iterClassicZeros 25 932 = 429 := by decide
  have hc3 : iterClassicZeros 25 429 = 202 := by decide
  have hc4 : iterClassicZeros 25 202 = 100 := by decide
  have hc5 : iterClassicZeros 25 100 = 55 := by decide
  have hc6 : iterClassicZeros 25 55 = 31 := by decide
  have hbase2 : iterClassicZeros 150 2048 = 31 := by
    have hs : 150 = 25 + 25 + 25 + 25 + 25 + 25 := rfl
    rw [hs, iterClassicZeros_add, iterClassicZeros_add, iterClassicZeros_add,
        iterClassicZeros_add, iterClassicZeros_add, hc1, hc2, hc3, hc4, hc5, hc6]
  have hsplit : n = 150 + (n - 150) := by omega
  constructor
  · rw [hsplit, iterTaper20Zeros_add, hbase1, iterTaper20Zeros_fixed_11]
  · rw [hsplit, iterClassicZeros_add, hbase2, iterClassicZeros_fixed_31]

theorem java_loop_4096_zeros_reaches_exact_fixed_point :
    iterTaper20Zeros 4096 2048 = 11 ∧ iterClassicZeros 4096 2048 = 31 :=
  java_loop_all_zeros_reaches_exact_fixed_point 4096 (by decide)

/-! ### Law VIc: Fixed-Point Reciprocal Table Soundness (`FastNoveltyCoder.RECIP`) -/

/-- Exact integer quotient lower bound for `RECIP[d] = 2^30 / d`. -/
def recipEntry (d : Nat) : Nat := (2 ^ 30) / d

/-- Precomputed reciprocal never overestimates integer division: `((num * (2^30 / d)) / 2^30) <= num / d`. -/
theorem recip_mul_le_exact (num d : Nat) (hd : 0 < d) :
    (num * recipEntry d) / (2 ^ 30) ≤ num / d := by
  unfold recipEntry
  have hK : 0 < (2 ^ 30 : Nat) := by decide
  by_cases hkd : (2 ^ 30 : Nat) / d = 0
  · simp [hkd]
  · have hpos : 0 < (2 ^ 30 : Nat) / d := Nat.pos_of_ne_zero hkd
    have hmul : d * ((2 ^ 30 : Nat) / d) ≤ 2 ^ 30 := Nat.mul_div_le (2 ^ 30) d
    have hmod : num % d < d := Nat.mod_lt num hd
    have hlt1 : (num % d) * ((2 ^ 30 : Nat) / d) < d * ((2 ^ 30 : Nat) / d) :=
      Nat.mul_lt_mul_of_pos_right hmod hpos
    have hlt2 : (num % d) * ((2 ^ 30 : Nat) / d) < 2 ^ 30 := Nat.lt_of_lt_of_le hlt1 hmul
    have hsplit : num = d * (num / d) + num % d := (Nat.div_add_mod num d).symm
    have hbound : num * ((2 ^ 30 : Nat) / d) < (num / d + 1) * (2 ^ 30) := by
      calc
        num * ((2 ^ 30 : Nat) / d)
            = (d * (num / d) + num % d) * ((2 ^ 30 : Nat) / d) :=
              congrArg (· * ((2 ^ 30 : Nat) / d)) hsplit
        _   = (num / d) * (d * ((2 ^ 30 : Nat) / d)) + (num % d) * ((2 ^ 30 : Nat) / d) := by
              rw [Nat.add_mul]
              congr 1
              ac_rfl
        _   ≤ (num / d) * (2 ^ 30) + (num % d) * ((2 ^ 30 : Nat) / d) := by
              exact Nat.add_le_add_right (Nat.mul_le_mul_left (num / d) hmul) _
        _   < (num / d) * (2 ^ 30) + 2 ^ 30 := Nat.add_lt_add_left hlt2 _
        _   = (num / d + 1) * (2 ^ 30) := by rw [Nat.succ_mul]
    have hdiv : (num * ((2 ^ 30 : Nat) / d)) / (2 ^ 30) < num / d + 1 :=
      (Nat.div_lt_iff_lt_mul hK).2 hbound
    exact Nat.lt_succ_iff.mp hdiv

end EimdlProofs.Law6
