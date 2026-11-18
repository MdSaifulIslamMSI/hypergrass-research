#!/usr/bin/env python3
"""
Assemble the standalone Lean 4 Formal Verification Companion Manuscript:
  paper/EIMDL_Lean4_Paper.tex -> paper/EIMDL_Lean4_Paper.pdf
Documents every machine-checked definition, theorem, and kernel axiom audit in
proofs/lean/ (leanprover/lean4:v4.36.0, zero `sorry`, pure core Init/Std).
"""

import os

BASE_DIR = r"C:\Users\mdsai\.gemini\antigravity\scratch\hypergrass"
PAPER_DIR = os.path.join(BASE_DIR, "paper")

LEAN4_TEX = r"""\documentclass[10pt]{article}
\usepackage[margin=0.72in,top=0.72in,bottom=0.74in]{geometry}
\usepackage{amsmath,amssymb,amsthm}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{array}
\usepackage{listings}
\usepackage{xcolor}
\usepackage{hyperref}

\hypersetup{
    colorlinks=true,
    linkcolor=blue!60!black,
    citecolor=blue!60!black,
    urlcolor=blue!60!black
}

\emergencystretch=3em
\tolerance=2000
\hbadness=2000
\sloppy

\lstset{
  basicstyle=\ttfamily\footnotesize,
  breaklines=true,
  columns=fullflexible,
  frame=single,
  rulecolor=\color{black!25},
  backgroundcolor=\color{black!3}
}

\newtheorem{theorem}{Theorem}
\newtheorem{definition}{Definition}
\renewcommand{\arraystretch}{1.12}

\title{\textbf{Machine-Checked Algebraic, Bit-Parallel SWAR, and Adaptive\\
Arithmetic Invariants for Evolutionary DNA Compression in Lean~4}}

\author{
  \textbf{Md Saiful Islam}\\
  Academy of Technology $\cdot$ HyperGraSS / RC-MDL / EIMDL Research Project\\
  {\small ORCID: \href{https://orcid.org/0009-0006-0570-9934}{\texttt{0009-0006-0570-9934}} $\cdot$ Email: \href{mailto:mdsaifulislam38msi@gmail.com}{\texttt{mdsaifulislam38msi@gmail.com}}}\\
  {\small GitHub: \url{https://github.com/MdSaifulIslamMSI/hypergrass-research}}\\
  {\small Zenodo DOI: \href{https://doi.org/10.5281/zenodo.24344825}{\texttt{10.5281/zenodo.24344825}} (\texttt{v3.1.1}: \href{https://doi.org/10.5281/zenodo.24344826}{\texttt{10.5281/zenodo.24344826}})}
}
\date{November 2026}

\begin{document}
\maketitle

\begin{abstract}
\noindent
Bit-parallel SIMD/SWAR nucleotide kernels, bijective gap encodings, causal reference DAGs, and adaptive integer range coders in lossless genomic compressors are notoriously vulnerable to off-by-one wraparound, division-by-zero interval collapse, and strand-flip sign bugs that randomized unit tests sample only sparsely. We present a self-contained, zero-\texttt{sorry}, zero-external-dependency formalization in \textbf{Lean~4} (\texttt{leanprover/lean4:v4.36.0}, pure \texttt{Init}/\texttt{Std}) covering the core mathematical and systems invariants of the Evolutionary Information Minimum Description Length (\textbf{EIMDL} / \texttt{HGS2} / \texttt{HGS3}) engine across five modules (\texttt{proofs/lean/}). In particular, we lift finite test checks to universal theorems: using \texttt{Std.Tactic.BVDecide} (CaDiCaL LRAT certificates), we certify the 5-instruction 64-bit SWAR 32-base reverse-complement involution and Klein four-group ($V_4$) XOR symmetry across all $2^{64}$ machine words and $2^{128}$ word pairs; via Presburger arithmetic (\texttt{omega}), we prove unbounded bijection of the mod-3 codon-wobble decomposition $\psi: \mathbb{N} \leftrightarrow \mathbb{N} \times \mathrm{Fin}\,3$ and verify the $[11, 4085]$ closed safety invariant and strict Lyapunov monotone convergence of EIMDL's two-adic tapered ($\tau=20$) binary range coder across arbitrary finite bitstreams; and via constructivist structural induction without classical axioms, we verify Mode~C reference DAG acyclicity (\texttt{[]}) and strictly positive reciprocal lookup ceiling soundness.
\end{abstract}

\section{Introduction: Closing the Gap Between Unit Assertions and Universal Invariants}

Scientific software systems rely heavily on unit test suites that check algebraic and arithmetic invariants on finite sample sets. In the Java~25 and C23/C++23/Rust~2024 implementation of \textbf{EIMDL}~\cite{rissanen1978modeling} (\texttt{VerifyCompressionLaws.java}, executing $171,941$ runtime assertions), critical low-level kernels were originally validated only up to finite cutoffs:
\begin{enumerate}
  \item \textbf{64-bit SWAR 32-base reverse complementation (\texttt{reverseComplementWord32}):} checked on $5,000$ pseudo-random \texttt{long} values out of $2^{64} \approx 1.84 \times 10^{19}$;
  \item \textbf{Bijective codon-wobble gap coding $\psi(g) = (\lfloor g/3 \rfloor, (g+1)\bmod 3)$:} checked on gaps $g \in [0, 2 \times 10^6]$ rather than arbitrary chromosome distances $g \in \mathbb{N}$;
  \item \textbf{Two-adic adaptive range probability tapering ($\tau = 20$):} checked on synthetic $4,096$-bit all-zero stress loops rather than proving interval non-collapse ($0 < p_t < 4096$) over \emph{all} binary message histories $\texttt{bits} \in \{0,1\}^*$;
  \item \textbf{Causal Mode~C pangenome DAG \& Stage-1/Stage-2 disjoint memory:} enforced by convention and runtime assertions without an acyclicity proof over multi-hop ancestor paths.
\end{enumerate}

This companion note documents our complete formalization under \texttt{proofs/lean/}, compiling cleanly in continuous integration (\texttt{lake build \&\& lake exe eimdl\_proofs\_check}) with zero \texttt{sorry}, zero custom axioms, and zero third-party libraries beyond Lean's core \texttt{Init} and \texttt{Std}.

\section{Law~I: Klein Four-Group $V_4$ Algebra and Bit-Vector SAT Verification of 64-Bit SWAR}

Representing $\Sigma_4 = \{\texttt{A}, \texttt{C}, \texttt{G}, \texttt{T/U}\}$ as \texttt{Base := BitVec 2} ($\{0, 1, 2, 3\}$), Watson--Crick complementation is bitwise XOR with \texttt{wcMask := 3\#2} ($\bar{b} = b \oplus 3$), while Kimura homology substitution differentials~\cite{kimura1980simple} are $\delta(r, t) = t \oplus r$. Purines ($b \in \{0, 2\}$) satisfy parity $(b \land 1) = 0$, giving the immediate transition test $\delta = 2 \iff \text{transition}$.

\begin{theorem}[Klein $V_4$ Strand \& Transition Invariance (\texttt{Law1KleinV4.lean})]
For all $r, t, \text{mask} \in \mathrm{BitVec}\,2$ and sequences $s \in \mathrm{List}\,\mathrm{Base}$ with $\mathrm{revComp}(s) = \mathrm{map}(\cdot \oplus 3, \mathrm{reverse}(s))$:
\begin{align}
  \mathrm{comp}(\mathrm{comp}(b)) &= b, &
  \mathrm{diff}(\mathrm{comp}(r), \mathrm{comp}(t)) &= \mathrm{diff}(r, t), \\
  \mathrm{isTransition}(\mathrm{comp}(r), \mathrm{comp}(t)) &= \mathrm{isTransition}(r, t), &
  \mathrm{revComp}(\mathrm{revComp}(s)) &= s, \\
  \mathrm{canonicalHash}(h, \mathrm{revComp}(s)) &= \mathrm{canonicalHash}(h, s).
\end{align}
\end{theorem}

Crucially, both the Java~25 and native C23 engines accelerate reverse-complementation across packed $32$-base $64$-bit words via \texttt{Long.reverseBytes}, nibble/crumb swaps, and bitwise inversion (\texttt{~~~u}):
\begin{lstlisting}
def swarReverseWord32 (u : BitVec 64) : BitVec 64 :=
  let u1 := reverseBytes64 u
  let u2 := ((u1 &&& 0x0F0F0F0F0F0F0F0F#64) <<< 4) ||| ((u1 >>> 4) &&& 0x0F0F0F0F0F0F0F0F#64)
  ((u2 &&& 0x3333333333333333#64) <<< 2) ||| ((u2 >>> 2) &&& 0x3333333333333333#64)

def swarReverseComplementWord32 (u : BitVec 64) : BitVec 64 := swarReverseWord32 (~~~u)
\end{lstlisting}
Using Lean~4's bit-blasting SMT/SAT tactic \texttt{bv\_decide} backed by CaDiCaL LRAT proof certificates verified inside Lean's kernel, we discharge all $2^{64}$ single words and $2^{128}$ register pairs universally:
\begin{lstlisting}
theorem swarReverseComplementWord32_involutive (u : BitVec 64) :
    swarReverseComplementWord32 (swarReverseComplementWord32 u) = u := by
  unfold swarReverseComplementWord32 swarReverseWord32 reverseBytes64
  bv_decide

theorem swar_v4_diff_strand_invariant (r t : BitVec 64) :
    swarReverseComplementWord32 t ^^^ swarReverseComplementWord32 r = swarReverseWord32 (t ^^^ r) := by
  unfold swarReverseComplementWord32 swarReverseWord32 reverseBytes64
  bv_decide
\end{lstlisting}
We likewise check boundary and interior $2$-bit lane permutations (\texttt{swar\_boundary\_lanes\_spec}, verifying $\mathrm{lane}_j(\mathrm{RC}_{32}(u)) = \mathrm{comp}(\mathrm{lane}_{31-j}(u))$ for $j \in \{0, 15, 16, 31\}$).

\section{Law~II: Unbounded Codon-Wobble Bijection and Sense/Antisense Stop Mirroring}

To exploit third-codon purifying selection~\cite{crick1966codon} without paying overhead on non-coding intergenic regions, EIMDL encodes each inter-SNP gap $g \in \mathbb{N}$ via $\psi(g) = (\lfloor g/3 \rfloor, (g+1)\bmod 3)$ and reconstructs $g$ via $\psi^{-1}(q, r) = 3q + ((r+2)\bmod 3)$. Modeling $\psi : \mathbb{N} \to \mathbb{N} \times \mathrm{Fin}\,3$ directly over unbounded natural numbers and applying \texttt{omega}, \texttt{Law2CodonWobble.lean} proves both left and right inverses (\texttt{wobble\_left\_inverse}, \texttt{wobble\_right\_inverse}), exact third-codon alignment equivalence $(\psi(g).2 = 0 \iff (g+1) \bmod 3 = 0)$, reading-frame strand involution (\texttt{bioPhase\_involutive}), and exact reverse-complement duality between sense stop codons $\{\texttt{TAA},\texttt{TAG},\texttt{TGA}\}$ and antisense codons $\{\texttt{TTA},\texttt{CTA},\texttt{TCA}\}$ across all $4^3 = 64$ trinucleotides (\texttt{antisenseStop\_iff\_revComp\_senseStop}).

\section{Laws~III, IV, and~VII: Constructive DAG Acyclicity, Non-Expansion, and Peak Heap Saving}

In Mode~C pangenome archiving, sequences ordered by minimum spanning tree similarity rank $r \in \mathbb{N}$ are only permitted to cite a parent panel donor $p = \mathrm{anc}(r)$ obeying $p < r$. Defining multi-hop reachability via the inductive predicate \texttt{AncestorPath anc a b}, structural induction yields \texttt{ancestor\_strictly\_smaller} ($b < a$) and immediate contradiction for self-cycles (\texttt{causal\_dag\_acyclic}, proving $\neg \mathrm{AncestorPath\,anc}\,r\,r$ with \textbf{zero axioms \texttt{[]}}). We further certify bounded fuel termination (\texttt{chainDepth\_le\_rank}, $\mathrm{chainDepth}(r) \le r$) and exact parent step depth (\texttt{chainDepth\_parent\_step}).

Under Law~IV, sub-block gating over the 5 generative models plus a raw $2$-bit fallback \texttt{raw2BitBytes n := (n + 3) / 4} satisfies \texttt{eimdl\_universal\_ceiling} ($L_{\text{block}} \le h + \lceil n/4 \rceil$) and linear bit bound \texttt{eimdl\_bits\_linear\_bound} ($8 L_{\text{block}} \le 2n + 8h + 6$), instantiating exact numeric bounds such as \texttt{eimdl\_500k\_crypto\_bound} ($\le 125,075$~B on $500$~kbp SHA-256 white noise). For Law~VII, strictly releasing Stage-1 index buffers before allocating Stage-2 superstream state yields \texttt{disjointPeak\_strict\_saving}: $\max(m_1, m_2) < m_1 + m_2$ for all $m_1, m_2 > 0$.

\section{Law~VI: Two-Adic Taper Safety Band $[11, 4085]$ and Fixed-Point Reciprocal Floor}

A subtle arithmetic-coding failure mode occurs when shift-based adaptation stalls: classic $p \mapsto p - \lfloor p/32 \rfloor$ over $M = 4096$ freezes at every $p \le 31$ ($0.757\%$ residual flag floor). EIMDL introduces symmetric two-adic rounding $\tau = 20$:
\begin{equation}
  \mathrm{taperStep}_0(\tau, p) = p - \left\lfloor \frac{p + \tau}{32} \right\rfloor, \qquad
  \mathrm{taperStep}_1(\tau, p) = p + \left\lfloor \frac{4096 - p + \tau}{32} \right\rfloor.
\end{equation}

\begin{theorem}[Closed Invariant Band, Strict Lyapunov Convergence \& Reciprocal Floor (\texttt{Law6TwoAdicTaper.lean})]
For any $0 \le \tau \le 31$, $\mathrm{taperStep}_0(\tau, p) = p \iff p \le 31 - \tau$ (\texttt{taperStep0\_fixed\_point\_iff}). At $\tau = 20$:
\begin{enumerate}
  \item \textbf{Single-Step \& Multi-Step Safety:} If $11 \le p_0 \le 4085$, then for \emph{every} finite binary message $\mathrm{bits} \in \mathrm{List}\,\mathrm{Bool}$,
  \[
    11 \le \mathrm{runTaperStream}(20, p_0, \mathrm{bits}) \le 4085 \qquad (\texttt{taper20\_stream\_invariant}),
  \]
  guaranteeing strictly positive symbol widths ($0 < p_t < 4096$) across infinite operation.
  \item \textbf{Lyapunov Monotonicity \& Asymptotic Trajectory:} $\mathrm{taperStep}_0(20, p) < p$ for $p \ge 12$, $\mathrm{taperStep}_1(20, p) > p$ for $p \le 4084$, and starting from $p_0 = 2048$, any horizon $n \ge 150$ all-zero bits reaches exact attractors $(11, 31)$ (\texttt{java\_loop\_all\_zeros\_reaches\_exact\_fixed\_point}, matching \texttt{VerifyCompressionLaws.java} at $n = 4096$).
  \item \textbf{Division-Free Reciprocal Table Soundness:} With precomputed $\mathrm{RECIP}[d] = \lfloor 2^{30}/d \rfloor$,
  \[
    \forall\,\mathrm{num} \in \mathbb{N},\; d > 0:\quad
    \left\lfloor \frac{\mathrm{num} \cdot \mathrm{RECIP}[d]}{2^{30}} \right\rfloor \le \left\lfloor \frac{\mathrm{num}}{d} \right\rfloor
    \qquad (\texttt{recip\_mul\_le\_exact}),
  \]
  proving fixed-point range-coder frequency normalization never overshoots exact division.
\end{enumerate}
\end{theorem}

\section{Kernel Axiom Audit and Continuous Integration Verification}

Every theorem exported by \texttt{EimdlProofs} is audited directly by Lean~4's kernel via \texttt{\#print axioms} inside \texttt{proofs/lean/Main.lean}. Table~\ref{tab:axioms} summarizes the clean dependency profile: zero \texttt{sorryAx}, zero user-declared postulates, and build verification enforced in GitHub Actions (\texttt{.github/workflows/ci.yml}, job \texttt{Formal Lean 4 Proofs (Laws I-VII)}).

\begin{table}[h!]
\centering
\caption{Complete Lean~4 kernel axiom audit across all five EIMDL modules (\texttt{proofs/lean/Main.lean}, \texttt{leanprover/lean4:v4.36.0}). Zero \texttt{sorryAx} and zero custom axioms appear anywhere in the project.}
\label{tab:axioms}
\resizebox{\textwidth}{!}{%
\begin{tabular}{llll}
\toprule
\textbf{Module} & \textbf{Verified Theorem} & \textbf{Proof Method} & \textbf{Lean Kernel Axioms} \\
\midrule
\texttt{Law1KleinV4} & \texttt{diff\_strand\_invariant} & Finite case reduction (\texttt{decide}) & \texttt{[propext, Quot.sound]} \\
\texttt{Law1KleinV4} & \texttt{swarReverseComplementWord32\_involutive} & Bit-blasted SAT (\texttt{bv\_decide} LRAT) & \texttt{[propext, Lean.ofReduceBool, Lean.trustCompiler, Quot.sound]} \\
\texttt{Law1KleinV4} & \texttt{swar\_v4\_diff\_strand\_invariant} & Bit-blasted SAT (\texttt{bv\_decide} LRAT) & \texttt{[propext, Lean.ofReduceBool, Lean.trustCompiler, Quot.sound]} \\
\texttt{Law1KleinV4} & \texttt{revComp\_involutive} & List induction + simplification & \texttt{[propext, Quot.sound]} \\
\texttt{Law2CodonWobble} & \texttt{wobble\_left\_inverse}, \texttt{wobble\_right\_inverse} & Linear integer arithmetic (\texttt{omega}) & \texttt{[propext, Quot.sound]} \\
\texttt{Law2CodonWobble} & \texttt{antisenseStop\_iff\_revComp\_senseStop} & Exhaustive $4^3=64$ codon check & \texttt{[propext, Classical.choice, Quot.sound]} \\
\texttt{Law3ModeCDag} & \texttt{causal\_dag\_acyclic}, \texttt{chainDepth\_le\_rank} & Structural \& strong induction & \textbf{\texttt{[]} (Zero axioms)} \\
\texttt{Law4UniversalBound} & \texttt{eimdl\_universal\_ceiling}, \texttt{disjointPeak\_strict\_saving} & Case split + \texttt{omega} & \texttt{[propext, Quot.sound]} \\
\texttt{Law6TwoAdicTaper} & \texttt{taper20\_stream\_invariant} & List induction + \texttt{omega} & \texttt{[propext, Classical.choice, Quot.sound]} \\
\texttt{Law6TwoAdicTaper} & \texttt{java\_loop\_4096\_zeros\_reaches\_exact\_fixed\_point} & 25-step congruence + \texttt{decide} & \texttt{[propext, Quot.sound]} \\
\texttt{Law6TwoAdicTaper} & \texttt{recip\_mul\_le\_exact} & Euclidean division algebra & \texttt{[propext]} \\
\bottomrule
\end{tabular}%
}
\end{table}

\section{Conclusion}

Formalizing EIMDL's algebraic symmetries, 64-bit SWAR bit-tricks, DAG ordering contracts, universal block ceilings, and two-adic probability invariant band in Lean~4 eliminates entire classes of genomic-compression edge cases at compile time without sacrificing runtime throughput or requiring external proof libraries. All proofs build offline cleanly under \texttt{proofs/lean/} alongside the 7-page bioinformatics paper (\texttt{EIMDL\_Main\_Manuscript.pdf}) and the 33-page technical monograph (\texttt{HGS2\_EIMDL\_Research\_Paper.pdf}).

\bibliographystyle{plain}
\bibliography{references}

\end{document}
"""


def write_lean4_paper():
    out_tex = os.path.join(PAPER_DIR, "EIMDL_Lean4_Paper.tex")
    with open(out_tex, "w", encoding="utf-8") as f:
        f.write(LEAN4_TEX)
    print(f"[OK] Wrote {out_tex} ({len(LEAN4_TEX):,} bytes)")


if __name__ == "__main__":
    write_lean4_paper()
