/-
  EimdlProofs.Law2CodonWobble
  Formal machine-checked proofs for Law II of EIMDL (`docs/COMPRESSION-LAW.md`):
    1. Exact unbounded bijection psi : Nat <-> Nat x Fin 3 between inter-SNP gap
       and (quotient, wobble phase) (`VerifyCompressionLaws.java` line 164)
    2. Synonymous 3rd-codon step equivalence: (gap + 1) % 3 = 0 iff wobblePhase = 0
    3. 6F-ORF-HMM strand-mirrored biological wobble phase involution (`bioPhase`)
    4. Exact reverse-complement mapping of sense {TAA, TAG, TGA} and antisense
       {TTA, CTA, TCA} stop codons across all 64 trinucleotides
-/
import EimdlProofs.Law1KleinV4

namespace EimdlProofs.Law2

open EimdlProofs.Law1

/-- Homologous substitution gap decomposition `psi(g) = (floor(g/3), (g+1) mod 3)`. -/
def wobbleEncode (g : Nat) : Nat × Fin 3 :=
  (g / 3, ⟨(g + 1) % 3, Nat.mod_lt _ (by decide)⟩)

/-- Exact inverse `psi^{-1}(q, r) = 3*q + ((r + 2) mod 3)`. -/
def wobbleDecode (p : Nat × Fin 3) : Nat :=
  3 * p.1 + (p.2.val + 2) % 3

/-- Left inverse across ALL natural numbers `g : Nat` (strengthens Java bounded check). -/
theorem wobble_left_inverse (g : Nat) : wobbleDecode (wobbleEncode g) = g := by
  unfold wobbleDecode wobbleEncode
  dsimp only
  omega

/-- Right inverse across ALL `(q, r) : Nat x Fin 3`. -/
theorem wobble_right_inverse (p : Nat × Fin 3) : wobbleEncode (wobbleDecode p) = p := by
  rcases p with ⟨q, ⟨r, hr⟩⟩
  unfold wobbleEncode wobbleDecode
  have h1 : (3 * q + (r + 2) % 3) / 3 = q := by omega
  have h2 : (3 * q + (r + 2) % 3 + 1) % 3 = r := by omega
  exact Prod.ext h1 (Fin.ext h2)

/-- Wobble quotient preservation under inverse. -/
theorem wobble_bijective :
    (∀ g, wobbleDecode (wobbleEncode g) = g) ∧
    (∀ p, wobbleEncode (wobbleDecode p) = p) :=
  ⟨wobble_left_inverse, wobble_right_inverse⟩

/-- Third-codon synonymous wobble predicate `(gap + 1) % 3 = 0` matches phase 0 iff. -/
theorem wobble_phase_zero_iff (g : Nat) :
    (wobbleEncode g).2.val = 0 ↔ (g + 1) % 3 = 0 := by
  unfold wobbleEncode
  dsimp only
  omega

/-- Injectivity of `wobbleEncode` over `Nat`. -/
theorem wobble_injective {g1 g2 : Nat} (h : wobbleEncode g1 = wobbleEncode g2) : g1 = g2 := by
  calc
    g1 = wobbleDecode (wobbleEncode g1) := (wobble_left_inverse g1).symm
    _  = wobbleDecode (wobbleEncode g2) := by rw [h]
    _  = g2                             := wobble_left_inverse g2

/-! ### 6F-ORF-HMM Phase Mirror & Sense/Antisense Stop-Triplet Involution -/

/-- `bioPhase` maps coordinate phase `c in {0,1,2}` to antisense wobble `2 - c` on `-` strand. -/
def bioPhase (isMinus : Bool) (c : Fin 3) : Fin 3 :=
  if isMinus then ⟨2 - c.val, by omega⟩ else c

/-- Double strand flip restores exact biological reading-frame phase. -/
theorem bioPhase_involutive (isMinus : Bool) (c : Fin 3) :
    bioPhase isMinus (bioPhase isMinus c) = c := by
  cases isMinus <;> revert c <;> decide

/-- Sense stop triplets `{TAA, TAG, TGA}` under `A=0, C=1, G=2, T=3`. -/
def isSenseStop (t : List Base) : Bool :=
  t == [3#2, 0#2, 0#2] || t == [3#2, 0#2, 2#2] || t == [3#2, 2#2, 0#2]

/-- Antisense reverse-complement stop triplets `{TTA, CTA, TCA}` terminating in forward phase 2. -/
def isAntisenseStop (t : List Base) : Bool :=
  t == [3#2, 3#2, 0#2] || t == [1#2, 3#2, 0#2] || t == [3#2, 1#2, 0#2]

/-- For every 3-base codon `[a, b, c]`, antisense stop iff reverse complement is a sense stop. -/
theorem antisenseStop_iff_revComp_senseStop (a b c : Base) :
    isAntisenseStop [a, b, c] = isSenseStop (revComp [a, b, c]) := by
  revert a b c; decide

end EimdlProofs.Law2
