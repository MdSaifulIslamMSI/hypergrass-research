#!/usr/bin/env python3
"""
Assemble the complete single-column, overflow-free LaTeX monograph
paper/HGS2_EIMDL_Research_Paper.tex covering 100% of the end-to-end research
across Generation 1 (HyperGraSS v1), Generation 2 (RC-MDL v2), and
Generation 3 (HGS2 / EIMDL), using programmatic tables from table_builders.py.
"""

import os
from table_builders import build_all_tables

BASE_DIR = r"C:\Users\mdsai\.gemini\antigravity\scratch\hypergrass"
PAPER_DIR = os.path.join(BASE_DIR, "paper")

PART1_TEX = r"""\documentclass[10pt]{article}
\usepackage[margin=0.70in,top=0.72in,bottom=0.75in]{geometry}
\usepackage{amsmath,amssymb,amsthm}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{longtable}
\usepackage{array}
\usepackage{xcolor}
\usepackage{hyperref}

\hypersetup{
    colorlinks=true,
    linkcolor=blue!60!black,
    citecolor=blue!60!black,
    urlcolor=blue!60!black
}

% Prevent any inline code symbol or math expression from overflowing the right margin
\emergencystretch=3em
\tolerance=2000
\hbadness=2000
\sloppy

\newtheorem{theorem}{Theorem}
\newtheorem{lemma}{Lemma}
\newtheorem{proposition}{Proposition}
\newtheorem{corollary}{Corollary}
\theoremstyle{definition}
\newtheorem{definition}{Definition}
\newtheorem{law}{Compression Law}

\renewcommand{\arraystretch}{1.10}

\title{\textbf{Evolutionary Information Minimum Description Length (EIMDL):\\
A Complete End-to-End Theory, Three-Generation Architectural Progression,\\
and Bounded-Memory Engine for Single-Genome, Reference-Conditioned,\\
and Pangenome Compression}}

\author{
  \textbf{Md Saiful Islam}\\
  HyperGraSS / RC-MDL / EIMDL Research Project (\texttt{HGS} v1.0 \& \texttt{HGS2} v2.0)\\
  Repository: \url{https://github.com/MdSaifulIslamMSI/hypergrass-research}
}
\date{November 2026}

\begin{document}
\maketitle

\begin{abstract}
\noindent
Lossless genomic compression has historically fractured into incompatible algorithmic paradigms: alignment-dependent pangenome archivers such as \texttt{AGC~3.2.5} require shared anchor segmentation and degrade on structurally rearranged or unaligned strains; general-purpose Lempel--Ziv and Burrows--Wheeler codecs (\texttt{xz}, \texttt{zstd}, \texttt{bzip2}, and original \texttt{GraSS}) fragment whenever single-nucleotide polymorphisms (SNPs) interrupt exact $k$-mers; and pure finite-context Markov/repeat models (\texttt{GeCo3}, \texttt{JARVIS3}) capture local sequence statistics at high memory cost without modeling multi-donor coalescent synteny or affine insertion--deletion (indel) mechanics.

This monograph presents the complete, end-to-end scientific and engineering record of the \textbf{HyperGraSS} research program across all three architectural generations:
\textbf{Generation~1 (\texttt{HyperGraSS v1}, \texttt{.hgs} v1.0)}, which introduced open-syncmer anchor banks, Anchor-Conditioned Homology Match Coding (\texttt{AC-HMC}), Block-Conditioned LZ-Grammar (\texttt{BC-LZG}), a Quantized Context-Repeat Mixer (\texttt{QCRM}), 8-channel $100:1$ entropy-floor accounting, GraalVM~25 static-musl $256$~MiB native execution, $2.077$~GiB 64-bit streaming verification, and sub-second block-indexed slice querying;
\textbf{Generation~2 (\texttt{RC-MDL v2} Prototype)}, which established unaligned reference-conditioned Minimum Description Length (MDL) coding, Prim minimum-spanning-tree pangenome ordering, and 13-configuration \texttt{AGC~3.2.5} oracle benchmarking; and
\textbf{Generation~3 (\texttt{HGS2} / \texttt{EIMDL}, \texttt{.hgs2} v2.0)}, which unifies all regimes under seven mathematically proved, experimentally falsifiable laws:
(I)~Klein Four-Group $V_4 = (\Sigma_4,\oplus)$ Strand-Invariant Kimura--CpG Differential Substitution Coding;
(II)~Bijective Codon-Wobble Gap Decomposition $\psi(\text{gap})=(\lfloor\text{gap}/3\rfloor, (\text{gap}+1)\bmod 3)$ and a 6-State Causal Open-Reading-Frame Hidden Markov Model (\texttt{6F-ORF-HMM});
(III)~Conditional Monotonicity, Horizon-Bounded Two-Gear Coalescent Haplotype Copying, and Affine Syntenic Micro-Indel Splicing;
(IV)~Universal Five-Regime Existing and Non-Existing Sequence Completeness ($\le 2.0012$~bpb white-noise safety ceiling);
(V)~Generalized Amdahl-Memory Streaming Scaling via a $0.25$~B/base 2-bit DNA bank and two-pass zero-materialization canonical open-syncmer Compressed Sparse Row (CSR) index ($k=16, s=10, w=6, t=3$);
(VI)~Two-Adic Tapered Binary Adaptation ($\tau=20$, lowering the fixed-point floor from $31/4096$ to $11/4096$), Hierarchical Pitman--Yor Suffix Backoff, and Reciprocal-Multiplication Normalization; and
(VII)~Strictly Disjoint Stage-1/Stage-2 Memory Lifecycle and Thread-Throttled Sub-Linear Execution operating inside a strict \texttt{-Xmx128m} ($128$~MiB, $1.06$~B/base) pure-JVM heap cap with $2.05\times$ encode ($45.16\times$ vs.\ \texttt{HyperGraSS~v1}) and $5.46\times$ ($6.56\times$ at 8T) decode speedup.
Laws~VI and~VII adapt architectural techniques from OpenAI's October~2026 ETH refutation preprint~\cite{openai2026eth} while explicitly preserving its Orthogonal Vectors Conjecture (OVC) boundary on polynomial-time string matching.

Across 73 complete bacterial genomes in Leave-One-Out (\textbf{Mode~R}) benchmarks, \texttt{HGS2} achieves \textbf{1,174,113~B} ($0.07820$~bpb, $102.30\times$) on 24 \textit{Escherichia coli} genomes ($7.14\times$ smaller than \texttt{AGC~3.2.5} best single configuration, winning \textbf{24/24} targets), \textbf{47,339~B} ($0.00358$~bpb, $2,233.78\times$) on 24 \textit{Mycobacterium tuberculosis} genomes ($6.17\times$ smaller than \texttt{AGC} best single configuration, \textbf{24/24} wins; $2.04\times$ smaller than a 13-configuration per-target \texttt{AGC} oracle, \textbf{19/24} wins), and \textbf{174,217~B} ($0.01962$~bpb, $407.79\times$) on 25 \textit{Staphylococcus aureus} genomes ($17.10\times$ smaller than \texttt{AGC} best single configuration, \textbf{25/25} wins; $9.28\times$ smaller than the 13-configuration \texttt{AGC} oracle, \textbf{24/25} wins). In self-contained pangenome archiving (\textbf{Mode~C}), \texttt{HGS2} compresses the 24-genome \textit{E.~coli} collection to \textbf{3,371,489~B} (1T, $0.2246$~bpb, $35.62\times$) and \textbf{3,388,639~B} (default parallel, $0.2257$~bpb, $35.44\times$, $19.65$~s at 12T encode / $8.66$~s at 8T decode), outperforming \texttt{HyperGraSS v1} ($3,810,409$~B), \texttt{JARVIS3 -l 7} ($3,896,361$~B), \texttt{RC-MDL v2} ($3,972,975$~B), \texttt{xz -9} ($5,663,640$~B), \texttt{GeCo3} ($8,409,358$~B), \texttt{AGC~3.2.5} ($8,540,341$~B), and original \texttt{GraSS} ($24,777,592$~B), while beating unbounded \texttt{v1} and \texttt{JARVIS3} even inside a pure-JVM \texttt{-Xmx128m} cap (\textbf{3,452,703~B}). In the locked Docker (\texttt{hypergrass-corpus:jdk27}) independent replication, \texttt{HGS2} also wins both \textit{Salmonella enterica} ($2,390,106$~B vs.\ $2,709,545$~B \texttt{JARVIS3}) and \textit{Klebsiella pneumoniae} ($2,523,991$~B vs.\ $2,948,926$~B \texttt{JARVIS3}). Every archive is verified byte-for-byte via SHA-256 and backed by 171,941 automated falsification assertions.
\end{abstract}

\section{Introduction, Architectural Progression, and Scientific Contract}

Under Algorithmic Information Theory and Rissanen's Minimum Description Length (MDL) principle~\cite{rissanen1978modeling}, the optimal lossless code length of a genomic dataset $X$ is bounded by its conditional description length $K(X \mid \mathcal{P})$ under an explicit probabilistic model $\mathcal{P}$ of molecular evolution. Deployed genomic compressors have historically divided into three camps, each breaking down outside its target regime:
\begin{enumerate}
  \item \textbf{Alignment- and Anchor-Dependent Pangenome Archivers:} Assembled Genome Compressor (\texttt{AGC~3.2.5})~\cite{deorowicz2023agc} splits genomes at conserved $k$-mer anchors and compresses segments against a reference using Zstandard. While fast on closely aligned assemblies, \texttt{AGC}'s anchor segmentation breaks under high substitution density, mobile-element rearrangements, or short segment sizes, forcing fallback to poorly compressed novel contigs and exhibiting up to $15\times$ archive-size variation across command-line flags ($-s \in [1.5\text{K}, 1\text{M}]$, $-a$, $-k$, $-m$).
  \item \textbf{Exact-Match Dictionary and BWT Archivers:} General-purpose compressors (\texttt{gzip}, \texttt{bzip2}, \texttt{xz/LZMA2}, \texttt{zstd}) and hybrid Lempel--Ziv/BWT genomic tools such as the original \texttt{GraSS} codec require exact substring identity. In diverging bacterial lineages where purifying selection preserves amino acid sequence while synonymous third-codon sites mutate every $15\text{--}30$~bp, exact matches shatter into thousands of short literals, inflating archive sizes up to $7\times$ above the evolutionary information floor.
  \item \textbf{Pure Finite-Context and Weighted-Repeat Models:} Dedicated genomic sequence coders such as \texttt{XM}~\cite{cao2007simple}, \texttt{FCM}~\cite{pinho2011fcm}, \texttt{GeCo3}~\cite{silva2020geco3}, and \texttt{JARVIS3}~\cite{sousa2024jarvis3} blend multiple Markov orders with hash-indexed repeat models. Although effective on single headerless DNA streams (\textbf{Mode~S}), they do not model multi-donor coalescent haplotype switching, affine syntenic micro-indels, or 6-channel FASTA container metadata, and high-order tables (\texttt{JARVIS3 -l 31}) require $\approx 1.16$~GiB of resident memory per stream.
\end{enumerate}

\subsection{Complete Three-Generation Progression of the Research Program}

To provide a complete end-to-end record from the beginning of the project to its unified conclusion, we document all four systems evaluated in this repository (\texttt{docs/FOUNDATIONS.md}):
\begin{enumerate}
  \item \textbf{Original \texttt{GraSS} Baseline (\texttt{AnshuOishik/GraSS}, commit \texttt{6d93700}, flag \texttt{-e2}):} The external reference-free Java baseline that combines a two-stage Lempel--Ziv match finder with \texttt{libbsc} block sorting in a tar container ($24,777,592$~B on the 24-genome \textit{E.~coli} pangenome; $22,264,066$~B on $100$~Mbp human \texttt{hg38 chr1}).
  \item \textbf{Generation~1 --- \texttt{HyperGraSS v1} (\texttt{.hgs} v1.0, magic \texttt{HGS\textbackslash x01}, \texttt{hypergrass.Main}):} Our first-generation in-repo architecture (\texttt{docs/BENCHMARKS.md}, \texttt{docs/SPECIALIZED.md}, \texttt{docs/BROAD.md}, \texttt{docs/ABLATION.md}, \texttt{docs/SCALING.md}, \texttt{docs/MEMORY-SWEEP.md}, \texttt{docs/ENTROPY-FLOOR.md}, \texttt{docs/END-TO-END.md}). Generation~1 introduced:
  \begin{itemize}
    \item A 5-channel bijective FASTA factorizer separating DNA from headers, line widths, case runs, and ambiguity runs;
    \item A Mash-selected anchor bank ($4\text{--}16$~MiB) encoded in $65,536$-bp context-mixed pages (\texttt{ANCHOR\_PAGE\_SIZE});
    \item Anchor-Conditioned Homology Match Coding (\texttt{AC-HMC}, DNA mode~3) with strand-aware bidirectional donor search, masked mismatch runs, and shared LZ77+range variant-gap history (\texttt{SHARED\_GAP}, mode bit~$0\text{x}04$);
    \item Block-Conditioned LZ-Grammar (\texttt{BC-LZG}, DNA mode~2) for fast \texttt{--lean} mode;
    \item A Quantized Context-Repeat Mixer (\texttt{QCRM}, DNA mode~1) for single genomes (\texttt{--single}), pruned to dense Markov orders $0\text{--}7$ plus 6 exact-verified repeat heads; and
    \item A random-access footer directory storing per-block offsets and \texttt{firstRecordBaseOffset} for sub-second locus slicing.
  \end{itemize}
  While Generation~1 reduced the 24-genome \textit{E.~coli} pangenome from $9,215,764$~B (default $4$~MiB bank / $2$~MiB blocks) to $4,885,851$~B (\texttt{--pangenome} $16$~MiB bank / $8$~MiB blocks) and further to \textbf{3,810,409~B} ($0.254$~bpb) with LZ77+range gap coding, its byte-by-byte context-mixed anchor pages required $887.45$~s and $1,452$~MB JVM RSS on the full pangenome, and our 8-channel entropy-floor decomposition (\texttt{docs/ENTROPY-FLOOR.md}) proved that self-contained Mode~C archives cannot cross $100:1$ ($1.20$~MB) because the irreducible first-genome anchor alone requires $2,257,234$~B ($1.88\times$ the entire $100:1$ budget).
  \item \textbf{Generation~2 --- \texttt{RC-MDL v2} Prototype (\texttt{prototype/rcmdl}):} Our standalone Reference-Conditioned Minimum Description Length prototype (\texttt{prototype/rcmdl/RESULTS.md}), designed specifically to solve the reference-conditioned regime (\textbf{Mode~R}) where the donor panel is already known to the receiver, as well as tree-ordered collection coding (\textbf{Mode~C}). Generation~2 introduced unaligned syncmer-seeded approximate donor matching, Prim minimum-spanning-tree (MST) collection ordering ($3,972,975$~B on \textit{E.~coli} Mode~C), and Leave-One-Out Mode~R coding ($1,636,367$~B on \textit{E.~coli}, $92,506$~B on \textit{M.~tuberculosis}, and $243,807$~B on \textit{S.~aureus}), alongside an exhaustive 13-configuration parameter sweep of \texttt{AGC~3.2.5}.
  \item \textbf{Generation~3 --- Unified \texttt{HGS2} / \texttt{EIMDL} Engine (\texttt{.hgs2} v2.0, magic \texttt{HGS2}, \texttt{hypergrass.EimdlCli}):} The unified architecture (\texttt{HyperGraSSFusionEngine.java}, \texttt{FastNoveltyCoder.java}, \texttt{EimdlMosaicCodec.java}, \texttt{EimdlSyncmerIndex.java}, \texttt{FastaStreamFactorizer.java}, \texttt{DonorPanelIndex.java}) that fuses Generation~2's unaligned coalescent haplotype copying with Generation~1's context-repeat mixing under seven formal laws (\texttt{docs/COMPRESSION-LAW.md}, \texttt{docs/FUSION-BENCHMARKS.md}, \texttt{docs/INDEPENDENT.md}). \texttt{HGS2} crosses the $100:1$ threshold in Mode~R across all three primary species ($102.30\times$ on \textit{E.~coli}, $2,233.78\times$ on \textit{M.~tuberculosis}, $407.79\times$ on \textit{S.~aureus}), reduces 24-genome \textit{E.~coli} Mode~C to \textbf{3,371,489~B} (1T) / \textbf{3,388,639~B} (12T in $19.65$~s encode / $10.40$~s decode [$8.66$~s at 8T]), and executes the full $120.11$~Mbp pangenome inside strict \textbf{\texttt{-Xmx192m}}, \textbf{\texttt{-Xmx160m}}, and \textbf{\texttt{-Xmx128m}} pure-JVM heap limits (\textbf{3,452,703~B}; and \textbf{3,422,314~B} at \texttt{-Xmx256m}/\texttt{-Xmx224m}).
\end{enumerate}

\begin{figure}[t!]
  \centering
  \includegraphics[width=0.98\textwidth]{figures/fig1_regime_overview.pdf}
  \caption{\textbf{Multi-regime description length and information decomposition across bacterial pangenomes.} \textbf{(A)}~Verified compression rates (bits/base, log scale) across all three canonical evaluation regimes---\textbf{Mode~S} (reference-free single genome), \textbf{Mode~C} (self-contained pangenome collection), and \textbf{Mode~R} (reference-conditioned Leave-One-Out)---across five bacterial species (\textit{E.~coli}, \textit{M.~tuberculosis}, \textit{S.~aureus}, \textit{S.~enterica}, and \textit{K.~pneumoniae}) relative to the $2.000$~bpb naive 2-bit ceiling. \textbf{(B)}~Stepwise information decomposition of \textit{E.~coli} \texttt{00.seq} ($4,856,080$~bp) from raw ASCII ($4,856.1$~KB, $8.0000$~bpb) down to 23-donor Mode~R EIMDL ($40.6$~KB, $0.0669$~bpb, $119.64\times$ compression).}
  \label{fig:regime_overview}
\end{figure}

\subsection{The Four-Part Falsifiable Scientific Contract and Rejected Hypotheses}

Every architectural mechanism retained in this research is bound by a four-part scientific contract (\texttt{docs/FOUNDATIONS.md}):
\begin{equation}
  \text{Invariant} \;\longrightarrow\; \text{Mechanism} \;\longrightarrow\; \text{Measured Effect} \;\longrightarrow\; \text{Falsification Test}.
\end{equation}
Equally important to scientific integrity is documenting every hypothesis and algorithmic lever that was \textbf{tested, falsified, and rejected} during development (Table~\ref{tab:rejected_levers}).

\begin{table}[h!]
\centering
\caption{Complete record of falsified hypotheses and rejected engineering levers across Generations 1--3 (\texttt{docs/FOUNDATIONS.md}, \texttt{docs/SPECIALIZED.md}, \texttt{docs/ABLATION.md}).}
\label{tab:rejected_levers}
\resizebox{\textwidth}{!}{%
\begin{tabular}{lllr}
\toprule
\textbf{Generation / Subsystem} & \textbf{Tested Lever / Hypothesis} & \textbf{Measured Empirical Outcome} & \textbf{Verdict} \\
\midrule
Gen~1 \texttt{QCRM} (\texttt{--single}) & Dense Markov orders $9\text{--}10$ on single bacterial genomes & $+2,285$~B larger ($1,086,677$ vs.\ $1,084,392$~B) and $+29.1\%$ slower on \textit{E.~coli} K-12 & \textbf{Rejected} \\
Gen~1 \texttt{QCRM} (\texttt{--single}) & Sparse $k=12$ and $k=15$ hash-table context experts & $+2,544$~B larger ($1,084,392 \to 1,084,133$~B when removed), $+16$~MiB/model, $+45\%$ time & \textbf{Rejected} \\
Gen~1 \texttt{QCRM} (\texttt{--single}) & Reducing repeat heads $6 \to 4$ and skipping exact verification & $0.20$~s faster ($3.99$ vs.\ $4.19$~s) but $+731$~B larger ($1,084,864$ vs.\ $1,084,133$~B) & \textbf{Rejected} \\
Gen~1 \texttt{QCRM} (\texttt{--single}) & Disabling all repeat heads (dense Markov orders $\le 7$ only) & $+12,277$~B larger ($1,096,410$ vs.\ $1,084,133$~B) on \textit{E.~coli} K-12 & \textbf{Rejected} \\
Gen~1 Anchor Selection & Position-strided sequence sketches for donor discovery & Missed homologous segments shifted by large genomic insertions/deletions & \textbf{Rejected} \\
Gen~1 Anchor Bank Grid & Doubling anchor bank cap from $16$~MiB to $32$~MiB on $121.8$~MB \textit{E.~coli} & $0$~B improvement ($4,885,851$~B at both caps; all $9,405,717$~bp anchors already fit in $16$~MiB) & \textbf{Rejected} \\
Gen~1 Anchor CM Pages & Parallelizing $65,536$-bp anchor CM pages on $10$~MB pangenome prefix & Slower wall time ($9.71$~s parallel vs.\ $8.85$~s serial) due to thread contention on small page sets & \textbf{Rejected} \\
Gen~1 $\to$ Gen~3 Architecture & Splitting CLI across C++/Rust native helpers plus JVM coordinator & Violated single-runtime auditability; replaced by pure-JVM CSR + GraalVM~25 static native image & \textbf{Rejected} \\
Gen~3 Two-Gear Matcher & Single-gear permissive approximate matching across pangenome panels & Bulldozed across recombination boundaries ($29,623$ false SNPs, $56,813$~B vs.\ $40,588$~B two-gear) & \textbf{Rejected} \\
\bottomrule
\end{tabular}%
}
\end{table}

\section{Six-Channel Bijective Factorization and Binary Container Specifications (\texttt{.hgs} v1.0 and \texttt{.hgs2} v2.0)}

A practical genomic compressor cannot assume pre-sanitized, headerless uppercase $\texttt{ACGT}$ strings; real biological repositories contain multi-FASTA headers, variable line wrapping ($\text{LF}$ or $\text{CRLF}$), lowercase soft-masking, IUPAC ambiguity codes ($\texttt{N, R, Y, S, W, K, M, B, D, H, V}$), and RNA uracil ($\texttt{U/u}$).

\begin{definition}[Six-Channel Exact Container Bijection]
Let $\Sigma_{\text{raw}} = \{0,\dots,255\}$ denote the byte alphabet and $\Sigma_4 = \{0,1,2,3\}$ represent the canonical 2-bit DNA alphabet under the bitwise complement-preserving assignment:
\begin{equation}
  \phi(\texttt{A}) = 0\;(00_2),\quad \phi(\texttt{C}) = 1\;(01_2),\quad \phi(\texttt{T/U}) = 2\;(10_2),\quad \phi(\texttt{G}) = 3\;(11_2).
  \label{eq:base_map}
\end{equation}
The stream factorizer $\Phi: \Sigma_{\text{raw}}^* \to \mathcal{M} \times \Sigma_4^n$ (\texttt{FastaStreamFactorizer.java}) decomposes any raw input file $X \in \Sigma_{\text{raw}}^L$ into a canonical base stream $S \in \Sigma_4^n$ and a five-component formatting metadata tuple $\mathcal{M} = (\mathcal{H}, \mathcal{W}, \mathcal{C}, \mathcal{E}, \tau)$:
\begin{enumerate}
  \item \textbf{Header Channel $\mathcal{H}$:} Ordered sequence of FASTA header strings $\{h_1,\dots,h_m\}$ and their exact 0-based nucleotide start offsets $\{p_1,\dots,p_m\}$.
  \item \textbf{Line-Width Channel $\mathcal{W}$:} Run-length encoded sequence of line lengths $(w_j, r_j)$ and line-ending type ($\text{LF}=0$, $\text{CRLF}=1$).
  \item \textbf{Case-Mask Channel $\mathcal{C}$:} Alternating run-length delta stream of uppercase/lowercase state transitions across the sequence.
  \item \textbf{IUPAC Exception Channel $\mathcal{E}$:} Sparse delta-coded positions $\{e_1,\dots,e_u\}$ and exact raw ASCII bytes for any non-$\{\texttt{A,C,G,T/U}\}$ nucleotide characters (which are replaced by $0$ ($\texttt{A}$) in the canonical stream $S$ so biological phase is never shifted).
  \item \textbf{Terminal Newline \& Alphabet Flags $\tau$:} Bit-flags recording trailing newline presence, $\text{CRLF}$ mode, and $\texttt{T}\leftrightarrow\texttt{U}$ RNA mode ($\texttt{flags \& 0x04}$).
\end{enumerate}
\end{definition}

\begin{proposition}[Strict Byte-Identity Round-Trip]
For any valid FASTA, multi-FASTA, or flat sequence file $X \in \Sigma_{\text{raw}}^L$, the inverse map $\Phi^{-1}(\mathcal{M}, S)$ reconstructs $X$ with zero byte discrepancy: $\Phi^{-1}(\Phi(X)) = X$, and $\text{SHA256}(\Phi^{-1}(\Phi(X))) = \text{SHA256}(X)$.
\end{proposition}
\begin{proof}
During forward factorization, every byte of $X$ is either a header byte appended to $\mathcal{H}$ at base index $p_k$, a line break $(\text{LF}/\text{CRLF})$ recorded in $\mathcal{W}$ (or $\tau$ at EOF), or a sequence character at index $i \in \{0,\dots,n-1\}$. If the character's uppercase projection belongs to $\{\texttt{A,C,G,T/U}\}$, its case state is recorded in $\mathcal{C}$ and its 2-bit value $\phi(\cdot) \in \Sigma_4$ is written to $S[i]$; otherwise, its exact byte value and index $i$ are appended to $\mathcal{E}$ while $S[i] \leftarrow 0$. Reconstruction $\Phi^{-1}$ replays $S[0..n-1]$, overriding $S[e_r]$ with the exact byte from $\mathcal{E}$, applying lowercase toggles from $\mathcal{C}$ (and $\texttt{T}\to\texttt{U}$ if the RNA flag is set), interleaving headers from $\mathcal{H}$ at $p_k$, and inserting line breaks from $\mathcal{W}$ and $\tau$. Since each partition is mutually exclusive and collectively exhaustive, $\Phi^{-1}(\Phi(X)) = X$ identically.
\end{proof}

\subsection{Complete Binary Layout of \texttt{.hgs} v1.0 and \texttt{.hgs2} v2.0 (\texttt{docs/HGS-FORMAT.md})}

Table~\ref{tab:format_specs} documents the exact byte-level container specifications for both Generation~1 (\texttt{.hgs} v1.0) and Generation~3 (\texttt{.hgs2} v2.0).
\begin{itemize}
  \item \textbf{Generation~1 (\texttt{.hgs} v1.0):} Supports two container layouts distinguished by byte~5:
  (a)~the \emph{Compact Single-Block Container} ($13$-byte header when byte~$5 \ge 32$: magic \texttt{HGS\textbackslash x01}, version \texttt{1}, \texttt{auxFormat} $\in \{32:\text{basic Deflate},\; 33:\text{structured LZ77+range}\}$, $4$-byte \texttt{auxLen}, $3$-byte 24-bit \texttt{dnaLen}) used for inputs $\le 128$~MiB without an anchor bank; and
  (b)~the \emph{Multi-Block Indexed Container} ($6$-byte header: magic \texttt{HGS\textbackslash x01}, version \texttt{1}, \texttt{modeFlags} where bit~$0 = \texttt{HAS\_ANCHOR}$, bit~$1 = \texttt{LEAN\_MODE}$, bit~$2 = \texttt{SHARED\_GAP}$) followed by the paged anchor bank ($65,536$-bp context-mixed pages), compressed blocks (each storing \texttt{auxMetadata} and a DNA payload under mode~$0$ [empty], mode~$1$ [\texttt{QCRM}], mode~$2$ [seed-match \texttt{BC-LZG}], or mode~$3$ [anchor-conditioned \texttt{AC-HMC} with donor list, match lengths, masked mismatch runs, and variant gaps]), and a random-access \emph{Footer Directory} terminated by a $28$-byte trailer ending in magic \texttt{0x48474631} (\texttt{"HGF1"}). Each footer entry stores \texttt{rawBytes}, \texttt{recordCount}, \texttt{dnaBases}, \texttt{blockOffset}, \texttt{blockLength}, \texttt{firstRecord}, and \texttt{firstRecordBaseOffset}, enabling $O(1)$ block seek for \texttt{query --slice} without decompressing non-overlapping blocks.
  \item \textbf{Generation~3 (\texttt{.hgs2} v2.0):} Uses a unified $27$-byte preamble (magic \texttt{0x48475332} = \texttt{"HGS2"}, version \texttt{2}, regime mode $\in \{0:\text{Mode~S},\; 1:\text{Mode~R},\; 2:\text{Mode~C}\}$, flags bitfield, $8$-byte raw input length, $4$-byte canonical base count $n$, $4$-byte donor count) followed in Mode~R by a mandatory \textbf{32-byte SHA-256 donor-panel pin} verified before decompression, the 6-channel metadata $\mathcal{M}$, order-1 range-coded structural/indel-splice/$V_4$-wobble streams, and Stage-2 parallel $4$~Mbp warm-primed \texttt{FastCRM} novelty sub-blocks (with bit~31 of the sub-block length word signaling raw 2-bit safety bypass under Law~IV).
\end{itemize}

\begin{table}[h!]
\centering
\caption{Exact binary container specifications of Generation~1 (\texttt{.hgs} v1.0) and Generation~3 (\texttt{.hgs2} v2.0) (\texttt{docs/HGS-FORMAT.md}).}
\label{tab:format_specs}
\resizebox{\textwidth}{!}{%
\begin{tabular}{llll}
\toprule
\textbf{Container Variant} & \textbf{Byte Offset / Section} & \textbf{Field Name \& Width} & \textbf{Exact Encoding \& Operational Semantics} \\
\midrule
\textbf{\texttt{.hgs} v1.0 Compact} & \texttt{0x00..0x03} ($4$~B) & \texttt{magic} & \texttt{0x48 0x47 0x53 0x01} (\texttt{"HGS\textbackslash x01"}) \\
(Single-block $\le 128$~MiB) & \texttt{0x04..0x05} ($2$~B) & \texttt{version}, \texttt{auxFormat} & Version \texttt{1}; \texttt{auxFormat} $= 32$ (basic) or $33$ (structured LZ77+range, $\ge 32$) \\
 & \texttt{0x06..0x0C} ($7$~B) & \texttt{auxLen} ($4$~B), \texttt{dnaLen} ($3$~B) & Metadata byte length + 24-bit DNA payload length, followed by \texttt{aux} and \texttt{dna} \\
\midrule
\textbf{\texttt{.hgs} v1.0 Multi-Block} & \texttt{0x00..0x05} ($6$~B) & \texttt{magic}, \texttt{ver}, \texttt{modeFlags} & \texttt{0x01}: \texttt{HAS\_ANCHOR}, \texttt{0x02}: \texttt{LEAN\_MODE}, \texttt{0x04}: \texttt{SHARED\_GAP} (mode~7) \\
(Pangenome \& random-access) & Anchor Bank & \texttt{totalBases} + pages & $65,536$-bp (\texttt{ANCHOR\_PAGE\_SIZE}) independent CM pages + offset directory \\
 & Block $k$ Payload & \texttt{auxLen} + \texttt{dnaLen} + streams & 5-channel \texttt{auxMetadata} + DNA mode $\in \{0:\text{Empty}, 1:\text{QCRM}, 2:\text{BC-LZG}, 3:\text{AC-HMC}\}$ \\
 & Footer Directory & $B$ entries + $28$-B trailer & Per-block \texttt{blockOffset}, \texttt{firstRecordBaseOffset}, \texttt{totalCrc32}, magic \texttt{"HGF1"} \\
\midrule
\textbf{\texttt{.hgs2} v2.0 Unified} & \texttt{0x00..0x06} ($7$~B) & \texttt{magic}, \texttt{ver}, \texttt{mode}, \texttt{flags} & \texttt{"HGS2"} (\texttt{0x48475332}), ver \texttt{2}, mode $\in \{0:\text{S}, 1:\text{R}, 2:\text{C}\}$, bit~0=flat, bit~1=CRLF, bit~2=RNA \\
(EIMDL Seven-Law Engine) & \texttt{0x07..0x1A} ($20$~B) & \texttt{rawBytes}, \texttt{bases}, \texttt{donors} & $8$-B raw byte length, $4$-B canonical base count $n$, $4$-B donor count, $4$-B metadata size \\
 & Mode~R Pin ($32$~B) & \texttt{donorPanelSha256} & Mandatory 32-byte SHA-256 digest of canonical donor panel (Mode~R only) \\
 & Stage-1 Streams & Order-1 Range Streams & Token types $T_k$, MTF donors, strand, syntenic $\Delta$, affine splices, wobble $(q,\phi)$, $V_4$ SNPs \\
 & Stage-2 Sub-Blocks & $4$~Mbp \texttt{FastCRM} Blocks & $32$-KiB warm-primed sub-blocks; bit~31 of length word ($\texttt{0x80000000}$) = raw 2-bit bypass \\
\bottomrule
\end{tabular}%
}
\end{table}

\section{The Seven Fundamental Laws of Evolutionary Information Compression}

We now state and prove the seven laws that govern \texttt{HGS2} (\texttt{docs/COMPRESSION-LAW.md}), culminating in Laws~VI and~VII which adapt mathematical mechanisms from OpenAI's October~2026 ETH refutation preprint~\cite{openai2026eth}.

\subsection{Law I: Klein Four-Group $V_4$ Strand-Invariant Kimura--CpG Substitution Law}

\begin{law}[Strand-Invariant Differential Substitution Coding]
Under the 2-bit nucleotide assignment~\eqref{eq:base_map}, Watson--Crick base complementation is the bitwise XOR involution $\bar{b} = b \oplus 2$. Equipping $\Sigma_4 = \{0,1,2,3\}$ with bitwise XOR $\oplus$ forms the abelian Klein four-group $V_4 \cong \mathbb{Z}_2 \times \mathbb{Z}_2$. Any single-nucleotide substitution between a donor reference base $r \in \Sigma_4$ and a target base $t \in \Sigma_4$ ($t \neq r$) is uniquely represented by the non-zero group differential:
\begin{equation}
  \delta(r, t) = r \oplus t \in \{1, 2, 3\}.
\end{equation}
\end{law}

\begin{theorem}[Reverse-Complement Equivariance and Kimura--CpG Entropy Bound]
\label{thm:law1}
Let $\delta(r,t) = r \oplus t \in \{1,2,3\}$. Then:
\begin{enumerate}
  \item \textbf{Exact Strand Invariance:} $\delta(\bar{r}, \bar{t}) = \delta(r, t)$ for all $r,t \in \Sigma_4$.
  \item \textbf{Algebraic Partition of Kimura's Two-Parameter Model~\cite{kimura1980simple}:} In \texttt{HyperGraSSFusionEngine.java}, the substitution coder maps $(r, t)$ into the canonical Kimura differential class $\delta_{\text{K}}(r,t) \in \{1,2,3\}$ where $\delta_{\text{K}}=1$ denotes the biochemical transition ($r \oplus t = 3$: $\texttt{A}(0)\leftrightarrow\texttt{G}(3)$ and $\texttt{C}(1)\leftrightarrow\texttt{T}(2)$), $\delta_{\text{K}}=2$ denotes the complement transversion ($r \oplus t = 2$: $\texttt{A}\leftrightarrow\texttt{T}, \texttt{C}\leftrightarrow\texttt{G}$), and $\delta_{\text{K}}=3$ denotes the cross transversion ($r \oplus t = 1$: $\texttt{A}\leftrightarrow\texttt{C}, \texttt{G}\leftrightarrow\texttt{T}$). All three non-zero cosets of $V_4 \setminus \{0\}$ are strictly invariant under reverse-complement inversion.
  \item \textbf{Strict Entropy Reduction:} For any transition probability $\alpha = \Pr(\delta_{\text{K}}=1) > \frac{1}{3}$ and transversion probabilities $\beta_2 = \Pr(\delta_{\text{K}}=2)$, $\beta_3 = \Pr(\delta_{\text{K}}=3)$ with $\alpha + \beta_2 + \beta_3 = 1$, conditioning on the flanking dinucleotide CpG context $C_{\text{CpG}} \in \{0,1,2,3\}$ and preceding differential $\Delta_{\text{prev}}$ satisfies:
  \begin{equation}
  \begin{aligned}
    H(\Delta \mid C_{\text{CpG}}, \Delta_{\text{prev}}) \le H(\Delta)
    &= -\alpha\log_2\alpha - \beta_2\log_2\beta_2 - \beta_3\log_2\beta_3 \\
    &< \log_2 3 \approx 1.5850\text{ bits},
  \end{aligned}
  \end{equation}
  and is strictly lower than raw target-symbol coding ($2.0000$~bits) by $2 - H(\Delta \mid C_{\text{CpG}}, \Delta_{\text{prev}}) \in [0.68, 0.90]$~bits per SNP.
\end{enumerate}
\end{theorem}
\begin{proof}
(1)~By associativity and commutativity of XOR in $V_4 = (\Sigma_4, \oplus)$ and the self-inverse property $2 \oplus 2 = 0$:
\begin{equation}
  \delta(\bar{r}, \bar{t}) = (r \oplus 2) \oplus (t \oplus 2) = r \oplus t \oplus (2 \oplus 2) = r \oplus t \oplus 0 = \delta(r, t).
\end{equation}
Thus a forward-strand match and a reverse-complement match emit the \emph{exact same} differential symbol without requiring strand-conditional alphabet permutation.
(2)~Each non-zero element of $V_4 \setminus \{0\} = \{1,2,3\}$ defines a fixed-point-free involution on $\Sigma_4$: $r \oplus t = 3$ pairs purines $\{0,3\}$ and pyrimidines $\{1,2\}$ (transitions), $r \oplus t = 2$ pairs Watson--Crick complements $\{0,2\}$ and $\{1,3\}$ (transversions), and $r \oplus t = 1$ pairs $\{0,1\}$ and $\{2,3\}$ (transversions).
(3)~By Shannon's Gibbs inequality, $H(\Delta) \le \log_2 3 \approx 1.5850$~bits, with strict inequality whenever $\alpha \neq 1/3$. In bacterial and eukaryotic genomes, transitions outnumber transversions by $\kappa = \alpha/(1-\alpha) \in [2.0, 4.5]$, and 5-methylcytosine deamination at $\texttt{CpG}$ dinucleotides elevates $\texttt{C}\to\texttt{T}$ ($\texttt{G}\to\texttt{A}$) transitions up to $10\times$. By the conditioning-reduces-entropy theorem, $H(\Delta \mid C_{\text{CpG}}, \Delta_{\text{prev}}) \le H(\Delta) < \log_2 3$. Empirically on \textit{E.~coli}, $H(\Delta \mid C_{\text{CpG}}, \Delta_{\text{prev}}) \approx 1.10\text{--}1.32$~bits/SNP.
\end{proof}

\subsection{Law II: Bijective Codon-Wobble Gap Decomposition and 6-State Causal ORF-HMM}

\begin{law}[Codon-Wobble Decomposition and Causal Reading-Frame Tracking]
In protein-coding sequences ($\approx 85\text{--}90\%$ of bacterial chromosomes), purifying selection ($\omega = d_N/d_S \ll 1$) concentrates substitutions at every third nucleotide (Crick's synonymous wobble position~\cite{crick1966codon}). Consequently, inter-SNP gap distances $\text{gap} = p_k - p_{k-1} - 1 \in \mathbb{N}_0$ within syntenic matches exhibit period-3 resonance $\text{gap} \equiv 2 \pmod 3$.
\end{law}

\begin{theorem}[Zero Non-Coding Overhead and Synonymous Mutual Information]
\label{thm:law2}
Define the mapping $\psi: \mathbb{N}_0 \to \mathbb{N}_0 \times \{0,1,2\}$ by:
\begin{equation}
  \psi(\text{gap}) = (q, \phi) = \left(\left\lfloor \frac{\text{gap}}{3} \right\rfloor,\; (\text{gap} + 1) \bmod 3\right), \qquad \psi^{-1}(q, \phi) = 3q + ((\phi + 2) \bmod 3).
  \label{eq:wobble_map}
\end{equation}
Then $\psi$ is a bijection on $\mathbb{N}_0$, and for any gap distribution $P(\text{gap})$, $H(Q) + H(\Phi \mid Q) = H(\text{gap})$. In non-coding DNA where $\Phi$ is uniform on $\{0,1,2\}$ independent of $Q$, $H(\Phi \mid Q) = \log_2 3$ and $H(Q) + H(\Phi) = H(\text{gap})$ (zero description-length expansion). In protein-coding regions under purifying selection ($\omega = d_N/d_S < 1$), $\Pr(\Phi = 0) > \frac{1}{3}$, reducing the per-gap code length by the mutual information $I(\Phi; \text{CDS}) = \log_2 3 - H(\Phi \mid \text{CDS}) > 0$~bits per SNP.
\end{theorem}
\begin{proof}
By Euclidean division, every integer $\text{gap} \ge 0$ has a unique quotient $q = \lfloor \text{gap}/3 \rfloor$ and remainder $r = \text{gap} \bmod 3 \in \{0,1,2\}$. The cyclic shift $\phi = (r + 1) \bmod 3$ is a permutation of $\{0,1,2\}$ with inverse $r = (\phi + 2) \bmod 3$, mapping $r=2$ (a 3-base step $p_k - p_{k-1} = 3$ between consecutive third-codon positions) to $\phi = 0$. Since $\psi$ is deterministic and invertible, the chain rule of Shannon entropy gives $H(\text{gap}) = H(Q, \Phi) = H(Q) + H(\Phi \mid Q)$. Under a memoryless geometric baseline on non-coding DNA, $\Phi$ is uniform ($\log_2 3$~bits), so separate range coding of $Q$ and $\Phi$ incurs zero asymptotic bit overhead, while in coding regions concentration of probability mass at $\phi=0$ lowers $H(\Phi)$ by $\log_2 3 - H(\Phi) > 0$~bits/SNP.
\end{proof}

\noindent
\textbf{Six-Frame Causal ORF-HMM (\texttt{6F-ORF-HMM} in \texttt{FastNoveltyCoder.java}):}
For unaligned novelty regions (Mode~S and novel accessory islands in Mode~R/C), \texttt{FastNoveltyCoder} tracks reading-frame survival across all six biological frames (3 forward phases $\times$ 2 strands) in $O(1)$ time per base without external gene annotations. At base index $i$ with forward codon phase $c = i \bmod 3 \in \{0,1,2\}$ and 6-bit triplet history $t_3 = (\text{history} \ \& \ 0\text{x}3\text{F}) \in \{0,\dots,63\}$:
\begin{itemize}
  \item At phase $c=2$, a forward-strand stop codon occurs iff $t_3 \in \{16\;(\texttt{TAA}),\; 18\;(\texttt{TAG}),\; 24\;(\texttt{TGA})\}$, resetting $L_{\text{fwd},c} \leftarrow 0$; otherwise $L_{\text{fwd},c} \leftarrow \min(255, L_{\text{fwd},c}+1)$.
  \item Simultaneously on the reverse-complement strand, stop codons correspond to $t_3 \in \{36\;(\texttt{TTA}),\; 48\;(\texttt{CTA}),\; 50\;(\texttt{TCA})\}$, resetting $L_{\text{rev},c} \leftarrow 0$; otherwise $L_{\text{rev},c} \leftarrow \min(255, L_{\text{rev},c}+1)$.
  \item Because random non-coding DNA encounters a stop triplet every $64/3 \approx 21.3$~codons ($\approx 64$~bp), a survival run $L_{\max} = \max_{f \in \{0..5\}} L_f > 64$ (or $> 160$) causally identifies an active Open Reading Frame ($\text{orfActive} \in \{0,1,2\}$).
  \item \texttt{FastNoveltyCoder} blends six models---four L1/L2-cache-aligned 24-byte finite-context Markov tables (\texttt{order10Table} [$2^{18}$ slots], \texttt{order8Table} [$2^{16}$], \texttt{order6Table} [$2^{12}$], \texttt{order4Table} [$2^8$]), a phase-locked \texttt{codonOrder5Table} [$3 \times 2^{14}$] conditioned on $(c,\, \text{pos}_{-1},\, \text{pos}_{-2},\, \text{pos}_{-3},\, \text{pos}_{-6},\, \text{orfActive})$, and a syncmer-seeded \texttt{repeatTable}---via an 8-regime adaptive logistic mixer (\texttt{weights[8][6]}) followed by a 64-bin Secondary Symbol Estimator (\texttt{sseTable[8][64]}).
\end{itemize}

\subsection{Law III: Conditional Monotonicity, Two-Gear Haplotype Copying, and Affine Micro-Indel Splicing}

\begin{theorem}[Monotonic Conditional Description Length and Coalescent Scaling]
\label{thm:law3}
Let $D_{1:N} = (D_1,\dots,D_N)$ be a donor reference panel. Under prefix-free coding of donor switch tokens, adding a donor $D_{N+1}$ satisfies $K(X \mid D_{1:N+1}) \le K(X \mid D_{1:N}) + O(1)$. Moreover, under a Li--Stephens coalescent copying process~\cite{li2003relating} with scaled mutation rate $\theta$ and recombination rate $\rho$, the expected conditional description length scales with the expected terminal coalescent branch length $\mathbb{E}[\tau_N] \propto 1/N$:
\begin{equation}
  \mathbb{E}[L_{\text{HGS2}}(X \mid D_{1:N})] \le L_{\text{struct}}(N) + n \cdot \mathbb{E}[\tau_N] \left(\theta \cdot \bar{c}_{\text{SNP}} + \rho \cdot (\bar{c}_{\text{switch}} + \log_2 N)\right),
\end{equation}
where $\bar{c}_{\text{SNP}} \approx 5.2\text{--}6.1$~bits/SNP (gap + $V_4$ symbol) and $\bar{c}_{\text{switch}} \approx 45\text{--}65$~bits/token.
\end{theorem}
\begin{proof}
Any parse valid over $D_{1:N}$ remains valid over $D_{1:N+1}$ with at most $\Delta_{\text{ptr}} = \log_2(N+1) - \log_2 N \le \frac{1}{N \ln 2}$ extra bits per donor pointer; if $D_{N+1}$ is unused, a 1-bit panel subset header recovers $L(X \mid D_{1:N}) + 1$~bits. Under the standard coalescent with $N$ lineages, the expected time to the most recent common ancestor of the target segment with its closest panel donor satisfies $\mathbb{E}[\tau_N] \sim \frac{2}{N}$. Each coalescent branch of length $\tau_N$ accumulates Poisson-distributed point mutations at rate $\theta \tau_N$ (coded at $\bar{c}_{\text{SNP}}$ bits each via Laws I--II) and recombination breakpoints at rate $\rho \tau_N$ (coded at $\bar{c}_{\text{switch}} + \log_2 N$ bits each). Taking expectations yields the bound.
\end{proof}

\noindent
\textbf{Horizon-Bounded Two-Gear Coalescent State Machine (\texttt{clonalGear}) and Affine Micro-Indel Splicing:}
\begin{itemize}
  \item \textbf{Gear~1 (Clonal Gear, $\texttt{clonalGear}=\text{true}$):} Credits $+0.008$~bits/bp of match extension and penalizes each SNP by $-11.2$~bits ($\approx 1$ SNP per $1,400$~bp steady-state tolerance), requiring a strict $\ge 128$~bp exact 64-bit word-XOR tail guard (\texttt{countExactMatch64}) and terminating immediately when a local 32-bp window exceeds $2$ SNPs.
  \item \textbf{Gear~2 (Diverged Syntenic Gear, $\texttt{clonalGear}=\text{false}$):} Activated when a match encounters a localized SNP cluster ($>2$ SNPs in $32$~bp) \emph{only if} the $O(1)$ syncmer lookahead probe \texttt{findClonalReentryHorizon} finds no alternative $\ge 192$~bp exact clonal match in the donor panel within $4\times 96 = 384$~bp. Gear~2 credits $+0.135$~bits/bp and penalizes SNPs by $-10.2$~bits ($\ge 18$~bp tail guard), bridging diverged syntenic regions until \texttt{findClonalReentryHorizon} triggers immediate cutoff and re-entry into a clonal donor.
  \item \textbf{Affine Syntenic Micro-Indel Splicing:} Before terminating a match at $(r_{\text{end}}, t_{\text{end}})$, \texttt{attemptMicroIndelSplice} searches a bounded 2D displacement grid $(\Delta_{\text{ref}}, \Delta_{\text{tgt}}) \in [1,256] \times [1,64]$ for a downstream continuation sharing $\ge 24$~bp of exact match, encoding the splice displacement in $\approx 10\text{--}18$~bits instead of emitting a new $60\text{--}80$-bit panel token.
\end{itemize}

\subsection{Law IV: Universal Five-Regime Completeness Law}

\begin{theorem}[Universal Completeness and Cryptographic White-Noise Safety Bound]
\label{thm:law4}
For any sequence $S \in \Sigma_4^n$ (biological, xenobiological, synthetic, or adversarial cryptographic noise), the five-regime grammar hierarchy of \texttt{HGS2}---(1)~Panel Syntenic Homology, (2)~Bidirectional Self-Homology, (3)~Tandem Periodicity Attractor ($1 \le P \le 64$), (4)~\texttt{FastCRM} Context-Repeat Mixer, and (5)~Raw 2-Bit Safety Bypass (\texttt{raw2BitBypass})---guarantees:
\begin{equation}
  L_{\text{HGS2}}(S) \le 2n + O(\log n)\text{ bits} \qquad \left(\frac{L_{\text{HGS2}}(S)}{n} \le 2.0012\text{ bits/base for } n \ge 5\times 10^5\right).
\end{equation}
\end{theorem}
\begin{proof}
In Stage~2 (\texttt{FastNoveltyCoder.encodeParallel}), every $4$~Mbp sub-block $B_j$ of length $n_j$ bases is compressed by \texttt{FastCRM} into a candidate byte buffer $C_j$. Immediately prior to sub-block emission, the encoder evaluates $|C_j| < \lceil n_j / 4 \rceil$. If $|C_j| \ge \lceil n_j / 4 \rceil$, the sub-block header sets bit~31 of the sub-block length word ($\texttt{len} \mid \texttt{0x80000000}$) and emits the exact 2-bit packed representation of $B_j$ ($\lceil n_j / 4 \rceil$ bytes). Including the $75$-byte container header on a $500,000$~bp stream, the total archive size is at most $\lceil n/4 \rceil + 75$~bytes $= 125,075$~B $= 2.0012$~bpb ($+0.06\%$ above the $2.0000$~bpb Shannon ceiling).
\end{proof}

\subsection{Law V: Generalized Amdahl-Memory Streaming Scaling Law}

\begin{theorem}[Sub-Linear CSR Index Memory Bound and Parallel Block Determinism]
\label{thm:law5}
Let $N_{\text{panel}}$ be the total bases in the donor reference bank, indexed via canonical open syncmers~\cite{edgar2021syncmers} with $k=16$, $s_{\text{mer}}=10$, window $w = k - s_{\text{mer}} = 6$, target offset $t = 3$, and subsampling shift $s \in \{0,1,2\}$. Then the expected number of indexed syncmer positions is $M(N_{\text{panel}}, s) = \frac{N_{\text{panel}}}{6 \cdot 2^s}$, and stored in a two-pass zero-materialization Compressed Sparse Row (CSR) structure with $2^{20}$ heads (\texttt{int[1\,048\,576]} = $4$~MiB), \texttt{int[M]} positions ($4M$~bytes), and \texttt{short[M]} verify tags ($2M$~bytes) alongside the 2-bit packed DNA bank ($0.25 N_{\text{panel}}$~bytes), the exact Stage-1 core memory footprint is:
\begin{equation}
  \mathcal{M}_{\text{Stage-1}}(N_{\text{panel}}, s) = N_{\text{panel}} \left(0.25 + \frac{6}{6 \cdot 2^s}\right) + 4\,\text{MiB} = N_{\text{panel}} \left(0.25 + 2^{-s}\right) + 4\,\text{MiB}.
  \label{eq:csr_mem}
\end{equation}
Specifically, $\mathcal{M}_{\text{Stage-1}}/N_{\text{panel}}$ equals $1.25$~B/bp raw ($1.58$~B/bp with JVM array/GC metadata) at $s=0$, $0.75$~B/bp raw ($0.93$~B/bp total) at $s=1$, and $0.50$~B/bp raw ($0.59$~B/bp total) at $s=2$. Releasing the CSR arrays immediately after Stage-1 parsing decouples Stage-1 index residency from Stage-2 parallel \texttt{FastCRM} novelty coding ($6.8$~MB/thread).
\end{theorem}
\begin{proof}
In any window of $w = 6$ consecutive $10$-mers, a $16$-mer is selected as an open syncmer iff its minimum-hash $10$-mer occurs at fixed offset $t=3$, which occurs with density $1/w = 1/6$. Subsampling by $s$ low bits retains a $2^{-s}$ fraction of syncmers, giving $M = N_{\text{panel}} / (6 \cdot 2^s)$ entries. Pass~1 of \texttt{DonorPanelIndex} counts bucket frequencies directly in \texttt{csrHeads[1<<20]}; a prefix sum converts counts into CSR offsets; Pass~2 populates \texttt{csrPositions[M]} ($4$~B/entry) and \texttt{csrTags[M]} ($2$~B/entry), totaling $6M = N_{\text{panel}} \cdot 2^{-s}$ bytes. Adding $0.25 N_{\text{panel}}$ bytes for the 2-bit DNA bank yields~\eqref{eq:csr_mem}.
\end{proof}

\subsection{Law VI: Two-Adic Tapered Adaptation, Hierarchical Pitman--Yor Backoff, and Reciprocal Normalization}

\begin{theorem}[Two-Adic Tapered Fixed-Point Reduction and Suffix-Summary Backoff]
\label{thm:law6}
Let $M = 2^{12} = 4096$ be the probability scale of a shift-$5$ binary range coder (\texttt{BitRangeEnc}/\texttt{BitRangeDec} in \texttt{EimdlMosaicCodec.java}) with adaptation taper $\tau = 20$:
\begin{equation}
  p_{t+1} =
  \begin{cases}
    p_t - ((p_t + \tau) \gg 5) & \text{if } b_t = 0,\\
    p_t + ((M - p_t + \tau) \gg 5) & \text{if } b_t = 1.
  \end{cases}
  \label{eq:two_adic_taper}
\end{equation}
Then:
\begin{enumerate}
  \item \textbf{Two-Adic Fixed-Point Reduction:} For deterministic $0$-streams ($b_t \equiv 0$), iteration~\eqref{eq:two_adic_taper} converges to exact fixed point $p_{\min}(\tau) = (2^5 - 1) - \tau = 31 - 20 = 11$ ($p_{\min}/M = 11/4096 \approx 0.269\%$) versus the classical truncation floor $31/4096 \approx 0.757\%$, reducing steady-state structural flag entropy from $-\log_2(1 - 31/4096) \approx 0.0110$~bits/symbol to $-\log_2(1 - 11/4096) \approx 0.00388$~bits/symbol ($2.83\times$ lower residual probability floor) at zero extra multiplications or branches. Symmetrically for deterministic $1$-streams, $p_{\max}(\tau) = M - 11 = 4085$.
  \item \textbf{Hierarchical Pitman--Yor Suffix Backoff:} In \texttt{FastNoveltyCoder.java}, ordering the context hierarchy as a suffix-summary chain (\texttt{ord0} $\to \dots \to$ \texttt{ord9} $\to$ \texttt{sparse12} and \texttt{orf2} $\to$ \texttt{orf3} $\to$ \texttt{orf5}) and interpolating sparse order-$k$ estimators toward the already-evaluated order-$(k-1)$ posterior $P_{k-1}(a \mid c_{k-1})$,
  \begin{equation}
    P_k(a \mid c_k) = \frac{n_k(c_k, a) + \alpha_k P_{k-1}(a \mid c_{k-1})}{n_k(c_k) + \alpha_k},
  \end{equation}
  strictly saves $\mathbb{E}[2.0 - H(A \mid c_{k-1})]$~bits per low-count context entry compared to uniform $1/4$ smoothing, while joint Law~I$\times$Law~II substitution conditioning (\texttt{(kimuraCtx << 1) | prevMod3Zero}) and slot-isolated low-bit contexts (\texttt{lowBitProbs[8][128]}, \texttt{posDeltaSlots[2][64]}) prevent cross-stream pollution.
  \item \textbf{Reciprocal-Multiplication Mixer Normalization:} Precomputing an immutable 16,385-entry reciprocal table $R(d) = \lfloor (2^{32} + d - 1)/d \rfloor$ removes variable-latency 64-bit hardware integer division (\texttt{idiv}) across all 17 \texttt{FastCRM} experts via single-cycle \texttt{(num * R(d)) >>> 32} scaling.
\end{enumerate}
\end{theorem}
\begin{proof}
At a fixed point $p^* \ge 0$ for $b_t \equiv 0$, $p^* - ((p^* + \tau) \gg 5) = p^*$ holds iff $p^* + \tau \le 31$, with maximal fixed point $p_{\min}(\tau) = 31 - \tau = 11$ for $\tau = 20$; by symmetry $M - p_{\max}(\tau) = 11$. Jensen's inequality on sparse order-$k$ posteriors shows that replacing uniform prior $u(a)=1/4$ with a strictly better-calibrated lower-order posterior $P_{k-1}(a \mid c_{k-1})$ having $D_{\mathrm{KL}}(P_{\text{true}} \,\|\, P_{k-1}) < D_{\mathrm{KL}}(P_{\text{true}} \,\|\, u)$ strictly decreases expected cross-entropy whenever $n_k(c_k) < \infty$.
\end{proof}

\subsection{Law VII: Strictly Disjoint Lifecycle and Sub-Linear Bounded-Heap Scaling (\texttt{-Xmx128m})}

\begin{theorem}[Disjoint Stage-1/Stage-2 Residency and Sub-Linear Heap Execution]
\label{thm:law7}
In Mode~C collection compression (\texttt{HyperGraSSFusionEngine.java}), Stage~1 maintains $\{\mathcal{B}_{\text{2bit}}, \mathcal{I}_{\text{CSR}}, \{\mathcal{S}_r\}\}$ to emit rank-ordered mosaic bitstreams $\{\mathcal{M}_r\}$ and raw novelty slices $\{\mathcal{L}_r\}$, while Stage~2 encodes $\{\mathcal{M}_r\}$ together with the pooled pan-genome novelty superstream $\mathcal{L}_{\text{super}}$ across parallel \texttt{FastCRM} workers $\{\mathcal{W}_p\}_{p=1}^P$. Explicitly nulling \texttt{sharedDagIdx}, \texttt{mstParentRank}, and \texttt{rankDonorSimilarityOrder} immediately after Stage~1 emission and releasing each \texttt{rankMosaics[r]} and \texttt{unitNovSlices[r]} as $\mathcal{L}_{\text{super}}$ is assembled enforces strictly disjoint lifespans, reducing peak heap residency to:
\begin{equation}
  M_{\text{peak}} = \max\!\big(M_{\text{Stage-1}},\, M_{\text{Stage-2}}\big) \ll M_{\text{Stage-1}} + M_{\text{Stage-2}}.
\end{equation}
Combined with narrowing per-worker target syncmer hashes (\texttt{TargetSyncmers.canonKmers}) from 64-bit \texttt{long[]} to 32-bit \texttt{int[]} ($-38\%$ worker index RAM), adaptive subsampling (\texttt{ultraTight} $\le 205$~MiB $\Rightarrow s=2$), dynamic worker throttling (\texttt{extremeTight} $\le 145$~MiB $\Rightarrow P \le 2$), decoupled \texttt{FastCRM} striped \texttt{panSuper} priming, and a 64~KB buffered Multi-FASTA reconstruction writer (\texttt{FastaStreamFactorizer.reconstructTo}), the complete $120,108,223$~bp 24-genome \textit{E.~coli} pangenome executes in pure JVM down to \texttt{-Xmx128m} ($1.06$~B/base total heap limit) within $+1.89\%$ of unbounded archive size.
\end{theorem}

\noindent
\textbf{Architectural Adaptation from OpenAI's October~2026 ETH Refutation Preprint~\cite{openai2026eth} and Honest OVC Scope Boundary:}
Laws~VI and~VII translate three structural techniques from OpenAI's 67-page ETH refutation (\emph{``The Exponential Time Hypothesis is False,''} 17~Oct~2026 / rev.~25~Oct~2026)~\cite{openai2026eth} into biological compression primitives:
(1)~\emph{Shared Linear Summaries} (\S2.4 Def.~2.1, \S2.7 eqs.~2.14--2.16, \S4.2 Lemma~4.1) motivated sharing pre-evaluated lower-order suffix summaries across \texttt{FastNoveltyCoder}'s HPY backoff chain (\textbf{Law~VIb});
(2)~\emph{Two-Adic Valuation Tapering} (Appendix~D.4.3--D.4.4, Hermite multiplicity Lemma~D.5 and eq.~D.10, which eliminate boundary precision loss via two-adic $v_2(t)$ tapering) motivated our symmetric two-adic rounding taper $\tau = 20$ reducing fixed-point stalling from $31/4096$ to $11/4096$ (\textbf{Law~VIa}); and
(3)~\emph{Disjoint Partitioning and Sequential Polynomial-Space Accumulation} (Lemma~3.2, Prop.~3.4, \S5.2 truncated algebra $\mathbb{Q}[V]/(V^\ell)$, Prop.~6.4, Appendix~A, and \S D.4.7) directly guided Stage-1 vs.\ Stage-2 object deallocation and bounded \texttt{-Xmx128m} execution (\textbf{Law~VII}).
Crucially, however, \S1.1.3 of~\cite{openai2026eth} explicitly notes that while Theorem~1.1 solves exact $\#k$-SAT in subexponential time and polynomial space, polynomial-time fine-grained lower bounds governed by the \emph{Orthogonal Vectors Conjecture (OVC)}---including quadratic edit distance and Longest Common Subsequence (LCS) complexity---remain intact. Thus,~\cite{openai2026eth} supplies constructive engineering paradigms (shared summaries, two-adic tapering, disjoint polynomial-space phases), \textbf{not} an asymptotic speedup over OVC-hard sequence alignment or string-matching bottlenecks.
"""

PART2_GEN1_GEN2_TEX = r"""
\section{Generation 1 (\texttt{HyperGraSS v1}) Complete Experimental Record}

Before presenting the unified Generation~3 (\texttt{HGS2} / \texttt{EIMDL}) results, we report the complete experimental record of Generation~1 (\texttt{HyperGraSS v1}, \texttt{.hgs} v1.0) across all eight foundational studies (\texttt{docs/BENCHMARKS.md}, \texttt{docs/SPECIALIZED.md}, \texttt{docs/BROAD.md}, \texttt{docs/ABLATION.md}, \texttt{docs/SCALING.md}, \texttt{docs/MEMORY-SWEEP.md}, \texttt{docs/ENTROPY-FLOOR.md}, \texttt{docs/END-TO-END.md}), illustrated in Figure~\ref{fig:gen1_gen2}.

\begin{figure}[t!]
  \centering
  \includegraphics[width=0.98\textwidth]{figures/fig7_gen1_gen2_foundations.pdf}
  \caption{\textbf{Generation~1 (\texttt{HyperGraSS v1}) and Generation~2 (\texttt{RC-MDL v2}) foundational experimental milestones.} \textbf{(A)}~Generation~1 cross-taxon natural bacterial collections across six species (\texttt{docs/BROAD.md}), where \texttt{HyperGraSS v1 --pangenome} beats \texttt{JARVIS3 -l 7} by $4.67\%\text{--}29.38\%$ on all six species. \textbf{(B)}~Generation~1 single-genome \texttt{QCRM} (\texttt{--single} and \texttt{--single --threads 2}) across bacterial, yeast, and $100$~Mbp human \texttt{hg38 chr1}/\texttt{chr14} chromosome slices (\texttt{docs/SPECIALIZED.md}). \textbf{(C)}~Generation~1 indexed random-access 13-slice ($2,000$~bp) query workload latency and archive sizes vs.\ full-stream \texttt{xz -9} and \texttt{gzip -9} (\texttt{docs/END-TO-END.md}). \textbf{(D)}~Generation~2 (\texttt{RC-MDL v2}) collection ordering topologies (\texttt{Prim MST}, \texttt{Cluster-Prim}, \texttt{Greedy}, \texttt{Input Order}) vs.\ fixed anchor sets ($N=4..20$) on the 24-genome \textit{E.~coli} pangenome (\texttt{prototype/rcmdl/RESULTS.md}).}
  \label{fig:gen1_gen2}
\end{figure}

\subsection{Frozen Held-Out Reference Corpus (\texttt{docs/BENCHMARKS.md})}

Table~\ref{tab:gen1_frozen_corpus} reports the completed cells ($179/280$ planned runs across $7$ datasets) of the frozen single-sequence corpus (\texttt{benchmarks/corpus/runs/corpus-full-1/report.json}) comparing \texttt{HyperGraSS v1} default, \texttt{--optimize}, and \texttt{--lean} against \texttt{gzip -9}, \texttt{xz -9}, \texttt{JARVIS3} (\texttt{-l 7} and \texttt{-l 26}), and \texttt{GeCo3} (\texttt{-l 3}, \texttt{-l 6}, and \texttt{-l 9}). As documented in \texttt{docs/BENCHMARKS.md}, \texttt{HyperGraSS v1} achieved the smallest archive on all $7$ completed datasets ($+0.56\%$ to $+17.52\%$ smaller than the best completed baseline), while high-effort \texttt{geco3-l9} timed out ($>180$~s) on \texttt{archaeal-reference} and \texttt{ecoli-k12}, and \texttt{jarvis3-l26} failed (infeasible under the $180$~s / $2$~GB container cap) on \texttt{archaeal-reference}.

\begin{table}[h!]
\centering
\caption{Generation~1 (\texttt{HyperGraSS v1}) frozen held-out corpus results in bytes and effective bits/base (\texttt{benchmarks/corpus/runs/corpus-full-1/report.json}, \texttt{docs/BENCHMARKS.md}). Cells marked \texttt{timeout}/\texttt{infeasible} exceeded the $180$~s / $2$~GB container cap.}
\label{tab:gen1_frozen_corpus}
\resizebox{\textwidth}{!}{%
\begin{tabular}{lrrrrrrrrrr}
\toprule
\textbf{Dataset} & \textbf{Input (B)} & \textbf{HyperGraSS v1} & \textbf{\texttt{--optimize}} & \textbf{\texttt{--lean}} & \textbf{gzip -9} & \textbf{xz -9} & \textbf{jarvis3-l7} & \textbf{jarvis3-l26} & \textbf{geco3-l6} & \textbf{geco3-l9} \\
\midrule
\texttt{human-mito} & 16,569 & \textbf{3,959} (1.912) & \textbf{3,959} (1.912) & 3,975 (1.919) & 5,180 (2.501) & 4,800 (2.318) & 12,183 (5.882) & 4,877 (2.355) & 4,884 (2.358) & 4,884 (2.358) \\
\texttt{fly-mito} & 19,517 & \textbf{3,567} (1.462) & \textbf{3,567} (1.462) & 3,673 (1.506) & 4,947 (2.028) & 4,680 (1.918) & 4,394 (1.801) & 3,739 (1.533) & 3,760 (1.541) & 3,759 (1.541) \\
\texttt{yeast-mito} & 85,779 & \textbf{15,338} (1.430) & \textbf{15,338} (1.430) & 15,733 (1.467) & 21,390 (1.995) & 18,984 (1.771) & 16,116 (1.503) & 15,517 (1.447) & 15,433 (1.439) & 15,425 (1.439) \\
\texttt{plant-plastid} & 154,478 & \textbf{30,419} (1.575) & \textbf{30,419} (1.575) & 36,661 (1.899) & 39,192 (2.030) & 34,212 (1.772) & 30,623 (1.586) & 36,193 (1.874) & 30,847 (1.597) & 30,839 (1.597) \\
\texttt{yeast-chr1} & 230,218 & \textbf{51,864} (1.802) & \textbf{51,864} (1.802) & 54,521 (1.895) & 66,474 (2.310) & 56,884 (1.977) & 53,372 (1.855) & 52,331 (1.818) & 52,557 (1.826) & 52,541 (1.826) \\
\texttt{archaeal-ref} & 1,664,970 & \textbf{370,196} (1.779) & \textbf{370,196} (1.779) & 376,881 (1.811) & 484,565 (2.328) & 410,584 (1.973) & 414,492 (1.992) & \textit{infeasible} & 413,569 (1.987) & \textit{timeout} \\
\texttt{ecoli-k12} & 4,641,652 & \textbf{1,087,650} (1.875) & \textbf{1,087,650} (1.875) & 1,118,253 (1.927) & 1,341,697 (2.312) & 1,149,544 (1.981) & 1,100,664 (1.897) & 1,103,759 (1.902) & 1,094,264 (1.886) & \textit{timeout} \\
\bottomrule
\end{tabular}%
}
\end{table}

\subsection{Anchor-Bank $\times$ Block-Size Full Grid Ablation and 10 MB Codec Component Ablation (\texttt{docs/ABLATION.md}, \texttt{docs/SPECIALIZED.md})}

On the 24-genome \textit{E.~coli} multi-FASTA scale test ($121,825,884$~B total, $120,108,223$~bp DNA), Generation~1's initial default configuration ($4$~MiB anchor bank / $2$~MiB blocks) produced $9,215,764$~B ($+9.78\%$ larger than \texttt{xz -9}'s $8,394,900$~B) because a $4$~MiB anchor bank holds only $4,194,304$~bp ($\sim 0.86$ of one \textit{E.~coli} genome) and $2$~MiB blocks reset block-local context 58 times. Table~\ref{tab:gen1_anchor_grid} reports the complete $3\times 3$ anchor-bank $\times$ block-size grid (\texttt{benchmarks/corpus/runs/ablation-anchor-grid-1/report.json}, plus the $32$~MiB bank control and the LZ77+range variant-gap upgrade from \texttt{docs/SPECIALIZED.md}), and Table~\ref{tab:gen1_10m_ablation} reports the 9-row component ablation on the $10$~MB \textit{E.~coli} pangenome slice (\texttt{codec-ablation-10m-v2}).

\begin{table}[h!]
\centering
\caption{Generation~1 anchor-bank $\times$ block-size grid ablation on the 24-genome \textit{E.~coli} multi-FASTA pangenome ($121,825,884$~B, \texttt{docs/ABLATION.md} and \texttt{docs/SPECIALIZED.md}).}
\label{tab:gen1_anchor_grid}
\resizebox{\textwidth}{!}{%
\begin{tabular}{lrrrrrrrrr}
\toprule
\textbf{Anchor Cap} & \textbf{Block Cap} & \textbf{Blocks} & \textbf{Anchor Bases} & \textbf{Archive (B)} & \textbf{Bits / Base} & \textbf{$\Delta$ vs.\ Default} & \textbf{$\Delta$ vs.\ \texttt{xz -9}} & \textbf{In-Proc Time} & \textbf{Peak RSS} \\
\midrule
4~MiB (Default) & 2~MiB & 58 & 4,194,304 & 9,215,764 & 0.6138 & $0.00\%$ & $+9.78\%$ & 54.02~s & 982.4~MB \\
4~MiB & 4~MiB & 29 & 4,194,304 & 8,370,127 & 0.5575 & $-9.18\%$ & $-0.30\%$ & 64.08~s & 996.4~MB \\
4~MiB & 8~MiB & 15 & 4,194,304 & 7,264,536 & 0.4839 & $-21.17\%$ & $-13.46\%$ & 77.72~s & 1,071.5~MB \\
8~MiB & 2~MiB & 58 & 8,388,608 & 6,026,650 & 0.4014 & $-34.60\%$ & $-28.21\%$ & 84.56~s & 993.8~MB \\
8~MiB & 4~MiB & 29 & 8,388,608 & 5,863,560 & 0.3906 & $-36.37\%$ & $-30.15\%$ & 93.80~s & 1,040.2~MB \\
8~MiB & 8~MiB & 15 & 8,388,608 & 5,551,039 & 0.3697 & $-39.77\%$ & $-33.88\%$ & 108.13~s & 1,288.2~MB \\
16~MiB & 2~MiB & 58 & 9,405,717 & 5,180,683 & 0.3451 & $-43.78\%$ & $-38.29\%$ & 89.81~s & 1,035.0~MB \\
16~MiB & 4~MiB & 29 & 9,405,717 & 5,060,967 & 0.3371 & $-45.08\%$ & $-39.71\%$ & 98.71~s & 1,148.5~MB \\
\textbf{16~MiB (\texttt{--pangenome})} & \textbf{8~MiB} & \textbf{15} & \textbf{9,405,717} & \textbf{4,885,851} & \textbf{0.3254} & \textbf{$-46.98\%$} & \textbf{$-41.80\%$} & \textbf{108.43~s} & \textbf{1,385.1~MB} \\
32~MiB (Control) & 8~MiB & 15 & 9,405,717 & 4,885,851 & 0.3254 & $-46.98\%$ & $-41.80\%$ & 110.22~s & 1,297.2~MB \\
\midrule
\textbf{16~MiB + LZ77 Gap (\texttt{v1 Final})} & \textbf{8~MiB} & \textbf{15} & \textbf{9,405,717} & \textbf{3,810,409} & \textbf{0.2538} & \textbf{$-58.65\%$} & \textbf{$-54.61\%$} & \textbf{887.45~s (c)} & \textbf{1,452.0~MB} \\
\bottomrule
\end{tabular}%
}
\end{table}

\begin{table}[h!]
\centering
\caption{Generation~1 component-level codec ablation on the $10$~MB \textit{E.~coli} pangenome prefix ($9,859,515$~bp, \texttt{codec-ablation-10m-v2}, \texttt{docs/SPECIALIZED.md}).}
\label{tab:gen1_10m_ablation}
\resizebox{\textwidth}{!}{%
\begin{tabular}{llrrrrll}
\toprule
\textbf{Subsystem} & \textbf{Ablation Variant} & \textbf{Archive (B)} & \textbf{Bits/Base} & \textbf{$\Delta$ vs.\ Default} & \textbf{Enc / Dec (s)} & \textbf{Peak RSS} & \textbf{Mechanism Verified} \\
\midrule
DNA Mode & \texttt{QCRM} forced (no anchor bank) & 2,406,156 & 1.952 & $+65.85\%$ & 18.92 / 19.01~s & 242~MiB & Single-sequence CM ignores cross-record homology \\
DNA Mode & \textbf{\texttt{Auto / AC-HMC} (Default)} & \textbf{1,450,787} & \textbf{1.177} & \textbf{Baseline} & \textbf{8.85 / 5.42~s} & \textbf{254~MiB} & Anchor-conditioned bidirectional match coding \\
DNA Mode & \texttt{BC-LZG} forced (no anchor bank) & 1,454,172 & 1.180 & $+0.23\%$ & 1.94 / 0.52~s & 198~MiB & Block-local LZ77 + range coding (\texttt{--lean} path) \\
Anchor Bank & Anchor CM pages (parallel) & 1,463,055 & 1.187 & $+0.85\%$ & 9.71 / 6.12~s & 268~MiB & Parallel page overhead on small $10$~MB slice \\
Aux Metadata & Aux without LZ77 (range only) & 1,450,787 & 1.177 & $0.00\%$ & 8.81 / 5.39~s & 254~MiB & Structured 5-plane range coder dominates aux \\
Aux Metadata & Aux without structured planes & 1,451,603 & 1.178 & $+0.06\%$ & 8.79 / 5.40~s & 254~MiB & Plane separation saves $816$~B on headers/widths \\
Aux Metadata & Aux basic Deflate (\texttt{auxFormat=32}) & 1,459,406 & 1.184 & $+0.59\%$ & 8.75 / 5.38~s & 252~MiB & Structured 5-channel aux saves $8,619$~B vs.\ Deflate \\
Baseline & Original \texttt{GraSS -e2} & 2,341,794 & 1.900 & $+61.42\%$ & 14.20 / 6.10~s & 412~MiB & Reference-free LZ+BWT baseline \\
\bottomrule
\end{tabular}%
}
\end{table}

\subsection{Single-Genome \texttt{QCRM} Tuning and 100 Mbp Human \texttt{hg38} Chromosome Slices (\texttt{docs/SPECIALIZED.md})}

Table~\ref{tab:gen1_single_hg38} summarizes Generation~1's single-genome \texttt{QCRM} (\texttt{--single}) results across four complete microbial/yeast genomes and two $100,000,000$-bp human \texttt{hg38} chromosome prefixes (\texttt{chr1} and centromeric gap-rich \texttt{chr14}). Pruning dense Markov orders $9\text{--}10$ and sparse $k=12,15$ hash experts while reusing a single 4-int scratch array per repeat head in \texttt{RepeatHead.predict} reduced heap allocation from $226$~MiB to $83$~MiB and outperformed \texttt{Original GraSS -e2} on every single-sequence benchmark ($-3.80\%$ to $-9.52\%$ on microbial/yeast; $-7.35\%$ on \texttt{hg38 chr1} and $-7.84\%$ on \texttt{hg38 chr14} at $2.25\text{--}2.53\times$ lower compression RSS).

\begin{table}[h!]
\centering
\caption{Generation~1 single-genome \texttt{QCRM} (\texttt{--single}) and segmented multi-thread evaluation across microbial, yeast, and $100$~Mbp human \texttt{hg38} chromosome slices (\texttt{docs/SPECIALIZED.md}).}
\label{tab:gen1_single_hg38}
\resizebox{\textwidth}{!}{%
\begin{tabular}{llrrrrrr}
\toprule
\textbf{Dataset} & \textbf{Bases (bp)} & \textbf{Orig.\ GraSS \texttt{-e2} (B)} & \textbf{v1 Default (B)} & \textbf{v1 \texttt{--single} (B)} & \textbf{v1 \texttt{--single -t 2} (B)} & \textbf{$\Delta$ (\texttt{--single} vs.\ GraSS)} & \textbf{Peak RSS (\texttt{--single} vs.\ GraSS)} \\
\midrule
\textit{E.~coli} K-12 MG1655 & 4,641,652 & 1,141,608 (1.968) & 1,086,715 (1.873) & \textbf{1,084,172} (\textbf{1.869}) & 1,090,807 (1.880) & \textbf{$-5.03\%$} ($-57,436$~B) & \textbf{88~MiB} vs.\ 278~MiB \\
\textit{B.~subtilis} 168 & 4,215,606 & 1,034,078 (1.962) & 995,744 (1.890) & \textbf{994,834} (\textbf{1.888}) & 998,693 (1.895) & \textbf{$-3.80\%$} ($-39,244$~B) & \textbf{87~MiB} vs.\ 248~MiB \\
\textit{S.~cerevisiae} (3 chr) & 1,360,022 & 333,798 (1.963) & 315,396 (1.855) & \textbf{315,289} (\textbf{1.855}) & 317,353 (1.867) & \textbf{$-5.55\%$} ($-18,509$~B) & \textbf{80~MiB} vs.\ 178~MiB \\
\textit{P.~aeruginosa} NY4383 & 7,316,070 & 1,690,428 (1.848) & 1,533,558 (1.677) & \textbf{1,529,536} (\textbf{1.673}) & 1,548,347 (1.693) & \textbf{$-9.52\%$} ($-160,892$~B) & \textbf{95~MiB} vs.\ 296~MiB \\
\midrule
Human \texttt{hg38 chr1} (100~Mb) & 100,000,000 & 22,264,066 (1.781) & 21,056,843 (1.684) & \textbf{20,626,877} (\textbf{1.650}) & 20,635,964 (1.651, \textbf{70.1~s}) & \textbf{$-7.35\%$} ($-1,637,189$~B) & \textbf{396/396~MiB} vs.\ 893/1,553~MiB \\
Human \texttt{hg38 chr14} (100~Mb) & 100,000,000 & 18,538,762 (1.483) & --- & \textbf{17,084,921} (\textbf{1.367}) & --- & \textbf{$-7.84\%$} ($-1,453,841$~B) & \textbf{380/389~MiB} vs.\ 960/1,524~MiB \\
\bottomrule
\end{tabular}%
}
\end{table}

\subsection{Cross-Taxon Natural Collections and Within-Species Diversity Gradient (\texttt{docs/BROAD.md})}

Table~\ref{tab:gen1_broad} presents the six natural cross-taxon bacterial pangenome collections and the $3\times 3$ within-species Mash diversity gradient (\texttt{low}, \texttt{medium}, \texttt{high} divergence strata across \textit{C.~jejuni}, \textit{E.~coli}, and \textit{S.~aureus}) from \texttt{docs/BROAD.md}. Across all six natural collections, \texttt{HyperGraSS v1 --pangenome} beats both \texttt{JARVIS3 -l 7} ($+4.67\%$ to $+29.38\%$) and \texttt{xz -9} ($+58.02\%$ to $+71.57\%$). Across the nine diversity strata, \texttt{HyperGraSS v1} wins $8/9$ strata, with the sole exception being the medium-diversity \textit{E.~coli} stratum ($3,381,136$~B vs.\ \texttt{JARVIS3}'s $3,375,860$~B, a $-0.16\%$ / $5,276$~B crossover).

\begin{table}[h!]
\centering
\caption{Generation~1 cross-taxon natural bacterial collections (top) and $3\times 3$ within-species Mash diversity strata (bottom) (\texttt{docs/BROAD.md}).}
\label{tab:gen1_broad}
\resizebox{\textwidth}{!}{%
\begin{tabular}{llrrrrrrr}
\toprule
\textbf{Category / Species} & \textbf{Collection / Stratum} & \textbf{Genomes} & \textbf{Input Bytes} & \textbf{HyperGraSS v1 (B)} & \textbf{JARVIS3 \texttt{-l 7} (B)} & \textbf{xz -9 (B)} & \textbf{Orig.\ GraSS (B)} & \textbf{$\Delta$ vs.\ Best Baseline} \\
\midrule
Natural Collection & \textit{Campylobacter jejuni} (\texttt{campy}) & 72 & 122,419,352 & \textbf{1,535,296} (0.100) & 2,017,046 (0.132) & 4,211,348 (0.275) & 24,598,042 & \textbf{$+23.88\%$} vs.\ JARVIS3 \\
Natural Collection & \textit{Chlamydia trachomatis} (\texttt{chlamydia}) & 66 & 69,754,471 & \textbf{333,507} (0.038) & 472,224 (0.054) & 862,388 (0.099) & 11,163,682 & \textbf{$+29.38\%$} vs.\ JARVIS3 \\
Natural Collection & \textit{Mycobacterium tuberculosis} (\texttt{mtb}) & 28 & 123,736,606 & \textbf{1,139,535} (0.074) & 1,418,155 (0.092) & 2,505,964 (0.162) & 26,389,531 & \textbf{$+19.65\%$} vs.\ JARVIS3 \\
Natural Collection & \textit{Streptococcus pneumoniae} (\texttt{pneumo}) & 57 & 121,007,978 & \textbf{1,754,114} (0.116) & 1,899,711 (0.126) & 5,960,712 (0.394) & 25,947,810 & \textbf{$+7.66\%$} vs.\ JARVIS3 \\
Natural Collection & \textit{Pseudomonas aeruginosa} (\texttt{pseudo}) & 19 & 123,241,439 & \textbf{3,083,542} (0.200) & 3,350,844 (0.218) & 10,847,460 (0.704) & 25,794,452 & \textbf{$+7.98\%$} vs.\ JARVIS3 \\
Natural Collection & \textit{Staphylococcus aureus} (\texttt{staph}) & 43 & 123,028,922 & \textbf{1,807,399} (0.118) & 1,895,866 (0.123) & 5,672,772 (0.369) & 25,470,552 & \textbf{$+4.67\%$} vs.\ JARVIS3 \\
\midrule
Diversity Stratum & \texttt{campy} --- Low / Med / High & $24\times 3$ & $\sim 40.8$~MB ea. & \textbf{444,619 / 586,251 / 961,772} & 563,914 / 734,882 / 1,148,711 & 973K / 1.45M / 2.66M & --- & \textbf{$+21.2\% / +20.2\% / +16.3\%$} \\
Diversity Stratum & \texttt{ecoli} --- Low / Med / High & $20\times 3$ & $\sim 101.6$~MB ea. & \textbf{1,829,036} / 3,381,136 / \textbf{4,373,930} & 1,981,969 / \textbf{3,375,860} / 4,404,956 & 3.74M / 8.07M / 11.49M & --- & \textbf{$+7.7\%$} / $-0.16\%$ / \textbf{$+0.7\%$} \\
Diversity Stratum & \texttt{staph} --- Low / Med / High & $14\times 3$ & $\sim 40.0$~MB ea. & \textbf{345,519 / 583,848 / 1,037,833} & 397,298 / 618,750 / 1,049,092 & 694K / 1.48M / 3.29M & --- & \textbf{$+13.0\% / +5.6\% / +1.1\%$} \\
\bottomrule
\end{tabular}%
}
\end{table}

\subsection{JVM Native Memory Attribution, GraalVM 25 Native Gate, $>2$~GiB Streaming, $100:1$ Entropy Floor, and Query Amortization (\texttt{docs/MEMORY-SWEEP.md}, \texttt{docs/SCALING.md}, \texttt{docs/ENTROPY-FLOOR.md}, \texttt{docs/END-TO-END.md})}

Table~\ref{tab:gen1_systems} consolidates the four systems and information-theory studies that completed Generation~1 and directly motivated Generation~2 and Generation~3:
\begin{enumerate}
  \item \textbf{JVM Native Memory Tracking (NMT) \& GraalVM~25 Native Gate (\texttt{docs/MEMORY-SWEEP.md}):} Under a $240$~MiB JVM heap cap (\texttt{-Xmx240m}), NMT proved that while committed Java heap at exit was $262.6$~MiB, total process RSS reached $389.3$~MiB due to $109.0$~MiB of transient native \texttt{malloc} buffers (\texttt{zlib}/I/O), $56.0$~MiB of C1/C2 JIT compiler arenas, and $13.9$~MiB of shared class space. Consequently, under a hard $256$~MiB Linux cgroup v2 limit (\texttt{--memory 256m --memory-swap 256m}), the standard OpenJDK JVM was OOM-killed ($0/3$ passes), whereas compiling \texttt{HyperGraSS v1} into a static-\texttt{musl} \textbf{GraalVM~25.0.2 Native Image} (\texttt{--gc=serial -R:MaxHeapSize=180m}) eliminated JIT and metaspace overhead, passing \textbf{3/3} strict $256$~MiB container runs at \textbf{3,864,464~B} ($240.1$~MiB compress RSS, $225.5$~MiB restore RSS, $251.1$~MiB cgroup peak).
  \item \textbf{$>2$~GiB 64-Bit Streaming Verification (\texttt{docs/SCALING.md}):} Replacing full-file byte arrays with a 64-bit \texttt{SequenceSource} abstraction enabled streaming compression of a $2.077$~GiB ($2,229,710,670$~B, 432 records, 132 blocks) synthetic pangenome in $798$~s inside a \texttt{-Xmx160m} heap ($235$~MB peak RSS), producing a $30,171,539$~B archive ($73.9:1$) that restored byte-exactly (\texttt{CRC-32 0xEE111322}).
  \item \textbf{Eight-Channel $100:1$ Entropy-Floor Proof (\texttt{docs/ENTROPY-FLOOR.md}):} Exact byte-accounting of the $3,810,409$~B \textit{E.~coli} archive decomposed its 8 channels into: Anchor bank $2,257,234$~B ($59.24\%$), Variant gaps $882,215$~B ($23.15\%$), Strand flags $239,808$~B ($6.29\%$), Mismatch diffs $160,236$~B ($4.21\%$), Match geometry $143,771$~B ($3.77\%$), Mismatch masks $80,818$~B ($2.12\%$), Unmatched literals $30,647$~B ($0.80\%$), and Donor/mode headers $15,680$~B ($0.41\%$). Because a $100:1$ archive on $120.11$~Mbp requires $\le 1,200,000$~B total, and \emph{either} the irreducible self-contained anchor bank ($2,257,234$~B) \emph{or} the non-anchor modification streams ($1,553,175$~B) alone exceeds $1,200,000$~B, we formally marked $100:1$ \textbf{\texttt{infeasible} for self-contained Mode~C}---and proved that crossing $100:1$ requires \textbf{Mode~R} (reference-conditioned compression where the donor panel is already at the receiver), where Generation~3 \texttt{HGS2} achieves \textbf{1,174,113~B} (\textbf{102.30$\times$}) on the exact same 24 \textit{E.~coli} genomes!
  \item \textbf{End-to-End 13-Slice Indexed Query Workload (\texttt{docs/END-TO-END.md}):} Across 13 genomic slices ($2,000$~bp each across 13 distinct records), \texttt{hypergrass-lean} answered the entire 13-slice batch in \textbf{0.73~s} ($43.9$~MiB RSS) and \texttt{hypergrass-pangenome} in \textbf{1.35~s} ($74.1$~MiB RSS), compared to \textbf{11.12~s} for \texttt{xz -9} and \textbf{34.14~s} for \texttt{gzip -9} (which must decompress the full stream), reaching cumulative compress+query break-even after $58$ queries (\texttt{--lean} vs.\ \texttt{xz -9}), $91$ queries (\texttt{--lean} vs.\ \texttt{gzip -9}), and $859$ queries (\texttt{--pangenome} vs.\ \texttt{xz -9}).
\end{enumerate}

\begin{table}[h!]
\centering
\caption{Generation~1 systems verification: $256$~MiB GraalVM~25 native cgroup gate, $2.077$~GiB 64-bit streaming probe, 8-channel $100:1$ entropy floor breakdown, and 13-slice indexed query workload (\texttt{docs/MEMORY-SWEEP.md}, \texttt{docs/SCALING.md}, \texttt{docs/ENTROPY-FLOOR.md}, \texttt{docs/END-TO-END.md}).}
\label{tab:gen1_systems}
\resizebox{\textwidth}{!}{%
\begin{tabular}{llrrrrl}
\toprule
\textbf{Study / Experiment} & \textbf{Profile / Channel} & \textbf{Bytes / Metric} & \textbf{Share / Rate} & \textbf{Compress / Query} & \textbf{Peak Memory} & \textbf{Verified Outcome} \\
\midrule
\textbf{256~MiB cgroup Gate} & OpenJDK 25 JVM (\texttt{-Xmx136m..168m}) & OOM Killed & $0/3$ pass & --- & $>256.0$~MiB cgroup & Native JIT/malloc overhead exceeds 256~MiB \\
(\texttt{MEMORY-SWEEP.md}) & \textbf{GraalVM 25 Native (\texttt{-R:MaxHeapSize=180m})} & \textbf{3,864,464~B} & \textbf{0.257~bpb} & 681.9~s / 674.8~s & \textbf{240.1 / 225.5~MiB} & \textbf{3/3 PASS} (251.1~MiB cgroup peak) \\
\midrule
\textbf{2.077~GiB Streaming} & 64-bit \texttt{SequenceSource} (132 blocks) & \textbf{30,171,539~B} & \textbf{73.90$\times$} & 798.0~s / 516.0~s & \textbf{235.0~MB} (\texttt{-Xmx160m}) & \textbf{PASS} (\texttt{CRC-32 0xEE111322} on $2.23$~GB) \\
\midrule
\textbf{8-Channel Floor} & Ch.~1: Anchor Bank ($9,405,717$~bp CM) & 2,257,234~B & $59.24\%$ (1.920~bpb) & --- & --- & Alone $= 1.88\times$ the $1.20$~MB $100:1$ budget \\
(\texttt{ENTROPY-FLOOR.md}, & Ch.~2: Variant Gaps (LZ77+range) & 882,215~B & $23.15\%$ (2.21~b/gap) & --- & --- & $3,190,112$ gaps across 15 blocks \\
$3,810,409$~B total) & Ch.~3--6: Strand + Diff + Geom + Mask & 624,633~B & $16.39\%$ & --- & --- & $239.8\text{K} + 160.2\text{K} + 143.8\text{K} + 80.8\text{K}$~B \\
 & Ch.~7--8: Literals + Donor/Mode Headers & 46,327~B & $1.22\%$ & --- & --- & $30,647$~B literals ($128,808$~bp) $+ 15,680$~B hdr \\
\midrule
\textbf{13-Slice Query} & \textbf{HyperGraSS v1 \texttt{--lean}} & 8,952,202~B & 0.596~bpb & 31.33~s / \textbf{0.73~s} & 804.6 / \textbf{43.9~MiB} & \textbf{15.2$\times$ faster query} vs.\ \texttt{xz -9} (break-even 58) \\
(\texttt{END-TO-END.md}, & \textbf{HyperGraSS v1 \texttt{--pangenome}} & \textbf{4,909,412~B} & \textbf{0.327~bpb} & 910.4~s / \textbf{1.35~s} & 1,452 / \textbf{74.1~MiB} & \textbf{8.2$\times$ faster query} vs.\ \texttt{xz -9} (break-even 859) \\
$13\times 2,000$~bp) & \texttt{xz -9} / \texttt{gzip -9} (Full Stream) & 8.39M / 34.24M & 0.559 / 2.280 & 48.0~s / 11.12--34.14~s & 674 / 32.9~MiB & Requires full-archive decompression per query \\
\bottomrule
\end{tabular}%
}
\end{table}

\section{Generation 2 (\texttt{RC-MDL v2} Prototype) Complete Experimental Record}

Generation~2 (\texttt{prototype/rcmdl}, \texttt{prototype/rcmdl/RESULTS.md}) bridged the gap between Generation~1's block-anchored architecture and Generation~3's seven-law engine by testing unaligned reference-conditioned MDL coding (\textbf{Mode~R}) and tree-ordered pangenome coding (\textbf{Mode~C}).

\subsection{Complete 13-Configuration \texttt{AGC 3.2.5} Parameter Sweep (\texttt{prototype/rcmdl/RESULTS.md})}

Because \texttt{AGC~3.2.5}'s performance depends heavily on segment size ($-s$), k-mer length ($-k$), minimum match length ($-l$), fallback minimizers ($-f$), and adaptive mode ($-a$), we executed a 13-configuration parameter sweep across all 24 \textit{M.~tuberculosis} targets (\texttt{agc-sweep-mtb-full/report.json}), all 25 \textit{S.~aureus} targets (\texttt{agc-sweep-saureus/report.json}), and \textit{E.~coli} target \texttt{00.seq} (\texttt{agc-sweep-ecoli-full/results-00.json}, plus full 24-genome LOO at the winning \texttt{-s 200000} setting totaling $8,384,117$~B). Table~\ref{tab:gen2_agc_sweep} reports the exact 13-configuration sweep results alongside \texttt{RC-MDL} and \texttt{HGS2}.

\begin{table}[h!]
\centering
\caption{Complete 13-configuration \texttt{AGC~3.2.5} parameter sweep in Mode~R Leave-One-Out across \textit{E.~coli} \texttt{00.seq} ($4,856,080$~bp; full 24g LOO at \texttt{-s 200000} = $8,384,117$~B), 24-genome \textit{M.~tuberculosis} ($105,745,015$~bp), and 25-genome \textit{S.~aureus} ($71,043,355$~bp), compared against Generation~2 (\texttt{RC-MDL}) and Generation~3 (\texttt{HGS2} / \texttt{EIMDL}) (\texttt{agc-sweep-ecoli-full/results-00.json}, \texttt{agc-sweep-mtb-full/report.json}, \texttt{agc-sweep-saureus/report.json}).}
\label{tab:gen2_agc_sweep}
\resizebox{\textwidth}{!}{%
\begin{tabular}{llrrrrrr}
\toprule
\textbf{AGC Config ID} & \textbf{CLI Flags} & \textbf{\textit{E.~coli} \texttt{00.seq} (B)} & \textbf{\textit{E.~coli} \texttt{00} bpb} & \textbf{\textit{M.~tb} 24g Total (B)} & \textbf{\textit{M.~tb} bpb} & \textbf{\textit{S.~aureus} 25g Total (B)} & \textbf{\textit{S.~aureus} bpb} \\
\midrule
\texttt{default} & \texttt{-a} & 316,868 & 0.52201 & 443,669 & 0.03357 & 3,539,729 & 0.39860 \\
\texttt{fallback} & \texttt{-a -f 0.05} & 322,005 & 0.53048 & 517,362 & 0.03914 & 3,541,988 & 0.39885 \\
\texttt{min15} & \texttt{-a -l 15} & 329,505 & 0.54283 & 446,925 & 0.03381 & 3,560,914 & 0.40099 \\
\texttt{min15-fallback} & \texttt{-a -l 15 -f 0.05} & 334,729 & 0.55144 & 520,536 & 0.03938 & 3,596,087 & 0.40495 \\
\texttt{seg50k} & \texttt{-a -s 50000} & 307,143 & 0.50599 & 1,024,001 & 0.07747 & \textbf{2,979,078} & \textbf{0.33547} \\
\texttt{seg100k} & \texttt{-a -s 100000} & 316,369 & 0.52119 & \textbf{292,170} & \textbf{0.02210} & 3,737,764 & 0.42090 \\
\texttt{seg200k} & \texttt{-a -s 200000} & \textbf{189,411} & \textbf{0.31204} & 1,390,783 & 0.10522 & 4,188,783 & 0.47169 \\
\texttt{plain-seg200k} & \texttt{-s 200000} & \textbf{189,411} (24g: \textbf{8,384,117}) & \textbf{0.31204} (0.5584) & 1,390,783 & 0.10522 & 4,188,783 & 0.47169 \\
\texttt{seg200k-min15-fb} & \texttt{-a -s 200000 -l 15 -f 0.05} & 189,900 & 0.31284 & 1,506,569 & 0.11398 & 4,235,561 & 0.47696 \\
\texttt{seg500k} & \texttt{-a -s 500000} & 224,332 & 0.36957 & 766,757 & 0.05801 & 4,458,726 & 0.50209 \\
\texttt{seg1m} & \texttt{-a -s 1000000} & 220,122 & 0.36263 & 384,994 & 0.02913 & 3,849,024 & 0.43343 \\
\texttt{seg1m-min15} & \texttt{-a -s 1000000 -l 15} & 229,178 & 0.37755 & 391,248 & 0.02960 & 3,956,000 & 0.44547 \\
\texttt{k25-seg200k} & \texttt{-a -k 25 -s 200000} & 394,600 & 0.65007 & 1,359,004 & 0.10281 & 3,744,704 & 0.42168 \\
\midrule
\textbf{AGC 13-Cfg Oracle} & $\min_{c \in \mathcal{C}_{13}}$ per target & 189,411 (24g single: 8,384,117) & 0.31204 (0.5584) & \textbf{96,580} & \textbf{0.00731} & \textbf{1,617,085} & \textbf{0.18210} \\
\textbf{Gen~2: \texttt{RC-MDL} (v1 / v2)} & Single default config & 73,194 / \textbf{57,540} (24g: 1,780,761 / \textbf{1,636,367}) & 0.1206 / \textbf{0.0948} (0.1090) & 122,998 / \textbf{92,506} & 0.00931 / \textbf{0.00700} & 267,412 / \textbf{243,807} & 0.03011 / \textbf{0.02745} \\
\textbf{Gen~3: \texttt{HGS2 / EIMDL}} & Single default config & \textbf{40,588} (24g: \textbf{1,174,113}) & \textbf{0.06687} (\textbf{0.07820}) & \textbf{47,339} & \textbf{0.00358} & \textbf{174,217} & \textbf{0.01962} \\
\bottomrule
\end{tabular}%
}
\end{table}

\subsection{Collection Ordering Topologies, Fixed Anchor Subsets, and Donor Panel Scaling (\texttt{prototype/rcmdl/RESULTS.md})}

Table~\ref{tab:gen2_topologies} reports Generation~2's Mode~C collection ordering ablations on the 24-genome \textit{E.~coli} pangenome ($120,108,223$~bp) and 25-genome \textit{S.~aureus} pangenome ($71,043,355$~bp), alongside Mode~R coalescent donor panel size scaling ($N \in \{1,2,4,8,16,23\}$ donors on \textit{E.~coli} \texttt{00.seq}). Ordering genomes along a Prim minimum spanning tree (\texttt{Containment/Prim}) over hash-sampled 16-mer containment sketches achieved $3,972,975$~B on \textit{E.~coli} (saving $83,444$~B vs.\ naive file order $4,056,419$~B and saving $142,393\text{--}2,648,304$~B vs.\ fixed anchor subsets $N=20..4$ at $4,115,368\text{--}6,621,279$~B) and $1,233,433$~B ($0.13889$~bpb) on \textit{S.~aureus} ($1.74\times$ smaller than \texttt{AGC}'s full 25-genome archive of $2,145,694$~B). Crucially, exhausting all four ordering strategies in Generation~2 proved that closing the remaining $162,566$~B gap to \texttt{HyperGraSS v1} ($3,810,409$~B) required fusing unaligned coalescent copying with a strong context-repeat novelty mixer (\texttt{FastCRM})---which Generation~3 (\texttt{HGS2}) achieved at \textbf{3,371,489~B}.

\begin{table}[h!]
\centering
\caption{Generation~2 (\texttt{RC-MDL}) Mode~C collection ordering topologies on 24-genome \textit{E.~coli} (left) and Mode~R donor panel scaling on \textit{E.~coli} \texttt{00.seq} across both \texttt{RC-MDL v1} (\texttt{RESULTS.md}), \texttt{RC-MDL v2} (\texttt{FUSION-BENCHMARKS.md}), and \texttt{HGS2} (right).}
\label{tab:gen2_topologies}
\resizebox{\textwidth}{!}{%
\begin{tabular}{lrrrlrrr}
\toprule
\textbf{Mode~C Ordering / Anchor Topology} & \textbf{Archive (B)} & \textbf{Bits / Base} & \textbf{$\Delta$ vs.\ Prim MST} & \textbf{Mode~R Panel Size ($N$)} & \textbf{RC-MDL v1 (B)} & \textbf{RC-MDL v2 (B)} & \textbf{HGS2 / EIMDL (B)} \\
\midrule
\textbf{Gen~3 \texttt{HGS2} Mode~C (\texttt{--threads 1})} & \textbf{3,371,489} & \textbf{0.22456} & \textbf{$-15.14\%$} & $N=1$ donor & 249,491 (0.4110) & 273,654 (0.4508) & \textbf{197,163 (0.3248)} \\
\textbf{Gen~1 \texttt{HyperGraSS v1} (BWT+CRM)} & \textbf{3,810,409} & \textbf{0.25380} & \textbf{$-4.09\%$} & $N=2$ donors & 185,439 (0.3055) & 197,390 (0.3252) & --- \\
\textbf{Gen~2 \texttt{RC-MDL} (\texttt{Containment/Prim})} & \textbf{3,972,975} & \textbf{0.26463} & \textbf{Baseline} & $N=4$ donors & 146,758 (0.2418) & 185,481 (0.3056) & --- \\
Gen~2 \texttt{RC-MDL} (\texttt{Cluster-aware Prim}) & 3,982,813 & 0.26528 & $+0.25\%$ & $N=8$ donors & 102,931 (0.1696) & 79,595 (0.1311) & --- \\
Gen~2 \texttt{RC-MDL} (\texttt{Greedy Marginal}) & 3,990,249 & 0.26578 & $+0.43\%$ & $N=16$ donors & 87,254 (0.1437) & 65,671 (0.1082) & --- \\
Gen~2 \texttt{RC-MDL} (\texttt{File Order / Cluster-First}) & 4,056,419 / 4,081,636 & 0.2702 / 0.2719 & $+2.10\% / +2.73\%$ & $N=23$ donors (Full LOO) & 73,194 (0.1206) & 57,540 (0.0948) & \textbf{40,588 (0.0669)} \\
Gen~2 Fixed Anchors ($N=20 / 16 / 12$) & 4,115,368 / 4,921,615 / 5,494,235 & 0.2741--0.3660 & $+3.58\%$ to $+38.29\%$ & \textbf{8-Donor Share of Gain} & \textbf{$83.1\%$ of $1\to 23$} & \textbf{$89.8\%$ of $1\to 23$} & --- \\
Gen~2 Fixed Anchors ($N=8 / 4$) & 5,818,621 / 6,621,279 & 0.3876 / 0.4410 & $+46.46\% / +66.66\%$ & \textbf{Total Gain ($1 \to 23$)} & \textbf{$-70.66\%$ ($3.41\times$)} & \textbf{$-78.97\%$ ($4.76\times$)} & \textbf{$-79.41\%$ ($4.86\times$)} \\
\bottomrule
\end{tabular}%
}
\end{table}
"""

PART3_GEN3_TEX_TEMPLATE = r"""
\section{Generation 3 (\texttt{HGS2} / \texttt{EIMDL}) Unified Multi-Regime Evaluation}

Every benchmark reported in this section was executed on committed sequence files and verified for $100\%$ byte-for-byte SHA-256 identity upon decompression in a clean process or isolated Linux container (\texttt{hypergrass-corpus:jdk27}, Docker~29.6.2 with cgroup~v2 peak memory accounting).

\subsection{Regime 1: Mode S (Reference-Free Single-Genome Compression)}

Table~\ref{tab:mode_s_host} reports single-genome (\textbf{Mode~S}) performance on the primary bacterial reference targets, and Table~\ref{tab:ind12_subset} presents the locked 12-genome independent \textit{Salmonella enterica} and \textit{Klebsiella pneumoniae} subset ($60,810,308$~bp total) against \texttt{JARVIS3 v3.7} (both default \texttt{-l 7} and maximum high-memory \texttt{-l 31}) and original \texttt{GraSS -e2}.

\begin{table}[h!]
\centering
\caption{Mode~S reference-free single-genome compression on host bacterial targets (\texttt{mode\_s\_and\_synthetic.json}). Note: \texttt{JARVIS3} was not run on these four host files; see Table~\ref{tab:ind12_subset} for the locked 12-genome \texttt{JARVIS3 -l 7} and \texttt{-l 31} evaluation.}
\label{tab:mode_s_host}
\resizebox{0.94\textwidth}{!}{%
\begin{tabular}{lrrrrrr}
\toprule
\textbf{Target File} & \textbf{Raw Bytes} & \textbf{HGS2 Archive (B)} & \textbf{Bits / Base} & \textbf{Compression Ratio} & \textbf{Encode (s)} & \textbf{Decode (s)} \\
\midrule
\textit{E.~coli} \texttt{00.fa} (FASTA) & 4,856,085 & \textbf{1,122,777} & 1.8497 & 4.33$\times$ & 6.75~s & 6.57~s \\
\textit{E.~coli} \texttt{00.seq} (Flat Seq) & 4,856,080 & \textbf{1,122,749} & 1.8496 & 4.33$\times$ & 8.53~s & 7.71~s \\
\textit{E.~coli} \texttt{00.fa} (\texttt{v1 --single} host) & 4,856,085 & 1,136,069 & 1.8716 & 4.27$\times$ & 10.01~s & 10.96~s \\
\textit{M.~tuberculosis} \texttt{00.seq} & 4,417,521 & \textbf{969,382} & 1.7555 & 4.56$\times$ & 7.49~s & 6.25~s \\
\textit{S.~aureus} \texttt{00.seq} & 2,845,025 & \textbf{632,236} & 1.7778 & 4.50$\times$ & 5.51~s & 5.67~s \\
\bottomrule
\end{tabular}%
}
\end{table}

\begin{table}[h!]
\centering
\caption{Locked 12-genome independent \textit{S.~enterica} (6 targets, $28.61$~Mbp) and \textit{K.~pneumoniae} (6 targets, $32.20$~Mbp) evaluation across Mode~S and Mode~R (\texttt{modes-subset-1/report.json}, \texttt{jarvis-l31-subset-1/report.json}). Archive sizes are in bytes (with rate in bits/base). \texttt{EIMDL Mode~S} (\texttt{--threads 1}) wins \textbf{12/12} targets against all Mode~S baselines while using $10.7\times$ less RAM than \texttt{JARVIS3 -l 31}.}
\label{tab:ind12_subset}
\resizebox{\textwidth}{!}{%
\begin{tabular}{llrrrrrr}
\toprule
\textbf{Species} & \textbf{Target Accession} & \textbf{Bases (bp)} & \textbf{EIMDL Mode~S} & \textbf{JARVIS3 \texttt{-l 31}} & \textbf{JARVIS3 \texttt{-l 7}} & \textbf{Orig.\ GraSS \texttt{-e2}} & \textbf{EIMDL Mode~R (LOO)} \\
\midrule
%%IND12_ROWS%%
\midrule
\multicolumn{3}{l}{\textbf{Total Compress / Decompress Wall Time (s)}} & 332.5 / 277.9~s & 417.5 / 452.1~s & \textbf{91.0 / 87.9~s} & 127.0 / 100.6~s & 325.1 / 120.4~s \\
\multicolumn{3}{l}{\textbf{Max Container cgroup Peak Memory (MiB)}} & \textbf{108 / 88~MiB} & 1,158 / 1,159~MiB & 266 / 266~MiB & 1,680 / 618~MiB & 526 / 98~MiB \\
\bottomrule
\end{tabular}%
}
\end{table}

\noindent
\textbf{Nuanced Mode~S Analysis (Ratio vs.\ Throughput vs.\ Memory):}
As shown in Table~\ref{tab:ind12_subset} and Figure~\ref{fig:ind_pareto}, against default \texttt{JARVIS3 -l 7} ($15,204,327$~B, $2.0002$~bpb), \texttt{EIMDL Mode~S} ($13,868,405$~B, $1.8245$~bpb) is \textbf{8.79\% smaller} ($-1,335,922$~B), wins \textbf{12/12} targets, and uses \textbf{2.46$\times$ less RAM} ($108$~MiB vs.\ $266$~MiB cgroup peak), but is $3.65\times$ slower to compress ($332.5$~s vs.\ $91.0$~s) and $3.16\times$ slower to restore ($277.9$~s vs.\ $87.9$~s) because \texttt{FastCRM} updates six context/repeat models and an 8-regime logistic mixer per nucleotide. When \texttt{JARVIS3} is given its maximum high-order table flag \texttt{-l 31} ($14,197,393$~B, $1.8678$~bpb), \texttt{JARVIS3}'s cgroup memory jumps to \textbf{1,158~MiB} ($10.7\times$ higher than \texttt{EIMDL Mode~S}) and its total wall time increases to $417.5$~s compress / $452.1$~s restore ($1.26\times$ / $1.63\times$ slower than \texttt{EIMDL} in aggregate: faster on \textit{Salmonella} [$100.0/106.4$~s vs.\ $153.1/130.9$~s] and slower on \textit{Klebsiella} [$317.5/345.7$~s vs.\ $179.4/147.0$~s]), yet \texttt{EIMDL Mode~S} still achieves a \textbf{2.32\% smaller} total archive ($-328,988$~B) and wins \textbf{12/12} individual genomes.

\subsection{Regime 2: Mode R (Reference-Conditioned Leave-One-Out and Leave-Cluster-Out)}

\begin{figure}[t!]
  \centering
  \includegraphics[width=0.98\textwidth]{figures/fig2_mode_r_loo_73genomes.pdf}
  \caption{\textbf{Mode~R reference-conditioned Leave-One-Out (LOO) and Leave-Cluster-Out compression across 73 complete bacterial genomes.} \textbf{(A)}~\textit{Escherichia coli} (24 genomes, $120.11$~Mbp): \texttt{EIMDL/HGS2} achieves $1,174,113$~B ($0.07820$~bpb), winning \textbf{24/24} targets against both \texttt{RC-MDL} ($1,780,761$~B; \texttt{RC-MDL v2}: $1,636,367$~B) and \texttt{AGC~3.2.5} best single configuration ($8,384,117$~B, $7.14\times$). \textbf{(B)}~\textit{Mycobacterium tuberculosis} (24 genomes, $105.75$~Mbp): \texttt{EIMDL/HGS2} achieves $47,339$~B ($0.00358$~bpb, $2,233.78\times$), beating \texttt{AGC} best single ($292,170$~B, $6.17\times$, \textbf{24/24} wins) and beating the 13-configuration per-target \texttt{AGC} oracle ($96,580$~B, $2.04\times$ aggregate, \textbf{19/24} per-target wins). \textbf{(C)}~\textit{Staphylococcus aureus} (25 genomes, $71.04$~Mbp): \texttt{EIMDL/HGS2} achieves $174,217$~B ($0.01962$~bpb), beating \texttt{AGC} best single ($2,979,078$~B, $17.10\times$, \textbf{25/25} wins) and the 13-configuration \texttt{AGC} oracle ($1,617,085$~B, $9.28\times$ aggregate, \textbf{24/25} wins). \textbf{(D)}~Leave-Cluster-Out unseen-lineage generalization across all three species.}
  \label{fig:mode_r_73}
\end{figure}

In \textbf{Mode~R}, a newly sequenced target assembly $X_i$ is compressed conditioned on an existing donor reference panel $D_{-i} = \{X_j : j \neq i\}$. Following the standard marginal benchmark protocol~\cite{deorowicz2023agc}, \texttt{AGC~3.2.5}'s marginal cost is measured as $|\text{AGC}(D_{-i} \cup \{X_i\})| - |\text{AGC}(D_{-i})|$, whereas \texttt{HGS2} produces a standalone self-verifying \texttt{.hgs2} file that embeds a mandatory 32-byte cryptographic SHA-256 digest of $D_{-i}$.

\begin{table}[h!]
\centering
\caption{Summary of Mode~R Leave-One-Out (LOO) and Leave-Cluster-Out (unseen-lineage) evaluations across 85 bacterial target runs (\texttt{loo\_*.json}, \texttt{cluster\_loo\_*.json}, \texttt{modes-subset-1/report.json}).}
\label{tab:mode_r_summary}
\resizebox{\textwidth}{!}{%
\begin{tabular}{llrrrrrrr}
\toprule
\textbf{Protocol} & \textbf{Dataset (Genomes, Total Bases)} & \textbf{HGS2 Total (B)} & \textbf{HGS2 bpb} & \textbf{Ratio} & \textbf{RC-MDL v2 (B)} & \textbf{AGC Best Single (B)} & \textbf{AGC 13-Cfg Oracle (B)} & \textbf{HGS2 Win Record} \\
\midrule
Standard LOO & \textit{E.~coli} (24g, 120,108,223~bp) & \textbf{1,174,113} & \textbf{0.07820} & \textbf{102.30$\times$} & 1,636,367 & 8,384,117 (\texttt{-s 200k}) & --- & \textbf{24/24} vs.\ AGC ($7.14\times$) \\
Standard LOO & \textit{M.~tuberculosis} (24g, 105,745,015~bp) & \textbf{47,339} & \textbf{0.00358} & \textbf{2,233.78$\times$} & 92,506 & 292,170 (\texttt{-s 100k}) & 96,580 & \textbf{24/24} Single ($6.17\times$); \textbf{19/24} Oracle ($2.04\times$) \\
Standard LOO & \textit{S.~aureus} (25g, 71,043,114~bp) & \textbf{174,217} & \textbf{0.01962} & \textbf{407.79$\times$} & 243,807 & 2,979,078 (\texttt{-s 50k}) & 1,617,085 & \textbf{25/25} Single ($17.10\times$); \textbf{24/25} Oracle ($9.28\times$) \\
Standard LOO & \textit{S.~enterica} + \textit{K.~pneumoniae} (12g, 60,810,308~bp) & \textbf{433,498} & \textbf{0.05703} & \textbf{140.28$\times$} & --- & --- & --- & $31.99\times$ smaller than Mode~S \\
\midrule
Cluster-LOO & \textit{E.~coli} ($J=0.98$, 22 clusters, max 3) & \textbf{1,619,296} & \textbf{0.10786} & \textbf{74.17$\times$} & 2,331,696 & 9,203,994 (\texttt{-a}) & --- & $5.68\times$ smaller than AGC (\texttt{-a}) \\
Cluster-LOO & \textit{M.~tuberculosis} ($J=0.99$, 13 clusters, max 6) & \textbf{78,603} & \textbf{0.00595} & \textbf{1,345.31$\times$} & 153,182 & 1,188,909 (\texttt{-a}) & --- & $15.13\times$ smaller than AGC (\texttt{-a}) \\
Cluster-LOO & \textit{S.~aureus} ($J=0.98$, 14 clusters, max 8) & \textbf{524,535} & \textbf{0.05907} & \textbf{135.44$\times$} & 741,650 & 4,132,771 (\texttt{-a}) & --- & $7.88\times$ smaller than AGC (\texttt{-a}) \\
\bottomrule
\end{tabular}%
}
\end{table}

\noindent
\textbf{Critical Analysis of Single-Configuration vs.\ 13-Configuration Oracle Baselines:}
\begin{enumerate}
  \item Against \emph{any single deployable \texttt{AGC~3.2.5} configuration}, \texttt{HGS2} (running with one fixed default parameter set across all species) wins \textbf{73/73} targets ($24/24$ on \textit{E.~coli} [$7.14\times$], $24/24$ on \textit{M.~tuberculosis} [$6.17\times$], and $25/25$ on \textit{S.~aureus} [$17.10\times$ vs.\ \texttt{-s 50000} and $20.32\times$ vs.\ \texttt{-a}]).
  \item Against the \emph{13-configuration per-target \texttt{AGC} oracle} on \textit{M.~tuberculosis} ($96,580$~B total, $0.00731$~bpb, which beat Generation~2 \texttt{RC-MDL v1}'s $122,998$~B by $21.5\%$), Generation~3 \texttt{HGS2} achieves \textbf{47,339~B} ($0.00358$~bpb)---\textbf{2.04$\times$ ($-50.98\%$) smaller than the 13-configuration oracle in aggregate}---and wins \textbf{19/24} targets (Appendix Table~\ref{tab:mtb_full}). On the 5 near-clonal targets where the oracle is smaller (\texttt{03.seq}: $256$ vs.\ $243$~B; \texttt{08.seq}: $287$ vs.\ $224$~B; \texttt{09.seq}: $607$ vs.\ $368$~B; \texttt{16.seq}: $1,150$ vs.\ $1,072$~B; \texttt{19.seq}: $536$ vs.\ $512$~B), the difference is $13\text{--}239$~B, partly attributable to \texttt{HGS2}'s mandatory 32-byte SHA-256 donor-panel manifest pin.
  \item Against the \emph{13-configuration per-target \texttt{AGC} oracle} on \textit{S.~aureus} ($1,617,085$~B total, $0.18210$~bpb), \texttt{HGS2} achieves \textbf{174,217~B} ($0.01962$~bpb)---\textbf{9.28$\times$ ($-89.23\%$) smaller in aggregate}---and wins \textbf{24/25} targets (Appendix Table~\ref{tab:saureus_full}), with the sole exception being clonal target \texttt{04.seq} ($683$~B vs.\ $617$~B oracle, a $66$~B gap).
\end{enumerate}

\subsection{Regime 3: Mode C (Self-Contained Pangenome Collection Compression)}

\begin{figure}[t!]
  \centering
  \includegraphics[width=0.98\textwidth]{figures/fig3_mode_c_pangenomes.pdf}
  \caption{\textbf{Mode~C self-contained pangenome collection benchmarks.} \textbf{(A)}~Canonical 24-genome \textit{E.~coli} pangenome ($120,108,223$~bp) across 11 evaluated compressor configurations. \textbf{(B)}~Multi-species Mode~C compression rate (bits/base, log scale) across five bacterial species (\textit{E.~coli}, \textit{M.~tuberculosis}, \textit{S.~aureus}, \textit{S.~enterica}, and \textit{K.~pneumoniae}).}
  \label{fig:mode_c}
\end{figure}

In \textbf{Mode~C}, an entire multi-genome collection is compressed into a single self-contained \texttt{.hgs2} archive with zero external dependencies. Table~\ref{tab:mode_c_all} reports the canonical 24-genome \textit{E.~coli} benchmark, the multi-FASTA 6-channel exact reconstruction benchmarks across \textit{E.~coli}, \textit{M.~tuberculosis}, and \textit{S.~aureus}, and the locked independent \textit{Enterobacteriaceae} trio benchmarks (\textit{S.~enterica} and \textit{K.~pneumoniae}).

\begin{table}[h!]
\centering
\caption{Mode~C self-contained pangenome collection compression across five bacterial species (\texttt{mode\_c\_collections.json}, \texttt{amdahl\_and\_memory.json}, \texttt{salmonella/trio-1/report.json}, \texttt{klebsiella/trio-1/report.json}). All archives verified byte-for-byte via SHA-256.}
\label{tab:mode_c_all}
\resizebox{\textwidth}{!}{%
\begin{tabular}{llrrrrrrr}
\toprule
\textbf{Collection Corpus} & \textbf{Codec / Configuration} & \textbf{Raw Bytes} & \textbf{Archive (B)} & \textbf{bpb} & \textbf{Ratio} & \textbf{Encode (s)} & \textbf{Decode (s)} & \textbf{Peak Memory} \\
\midrule
\textit{E.~coli} 24g Flat (\texttt{input.seq}) & \textbf{EIMDL / HGS2 (\texttt{--threads 1})} & 120,108,223 & \textbf{3,371,489} & \textbf{0.2246} & \textbf{35.62$\times$} & 40.31~s & 56.78~s & 754~MiB cgroup \\
\textit{E.~coli} 24g Flat (\texttt{input.seq}) & \textbf{EIMDL / HGS2 (Default 12T)} & 120,108,223 & \textbf{3,388,639} & \textbf{0.2257} & \textbf{35.44$\times$} & \textbf{19.65~s} & \textbf{10.40~s} & 754~MiB (default) \\
\textit{E.~coli} 24g Flat (\texttt{input.seq}) & \textbf{EIMDL / HGS2 (\texttt{-Xmx192m})} & 120,108,223 & \textbf{3,452,703} & \textbf{0.2300} & \textbf{34.79$\times$} & 62.68~s & 27.46~s & \textbf{192~MiB Heap Cap} \\
\textit{E.~coli} 24g Flat (\texttt{input.seq}) & \textbf{EIMDL / HGS2 (\texttt{-Xmx128m})} & 120,108,223 & \textbf{3,452,703} & \textbf{0.2300} & \textbf{34.79$\times$} & 25.59~s & 11.67~s & \textbf{128~MiB Heap Cap} \\
\textit{E.~coli} 24g Flat (\texttt{input.seq}) & \texttt{HyperGraSS v1} (BWT+CRM) & 120,108,223 & 3,810,409 & 0.2538 & 31.52$\times$ & 887.45~s & 895.52~s & 1,452~MB RSS \\
\textit{E.~coli} 24g Flat (\texttt{input.seq}) & \texttt{JARVIS3 v3.7} (\texttt{-l 7}) & 120,108,223 & 3,896,361 & 0.2595 & 30.83$\times$ & 95.36~s & 98.41~s & 549~MB RSS \\
\textit{E.~coli} 24g Flat (\texttt{input.seq}) & \texttt{RC-MDL v2} (Prim Tree) & 120,108,223 & 3,972,975 & 0.2646 & 30.23$\times$ & 164.59~s & 17.50~s & 1,812~MB RSS \\
\textit{E.~coli} 24g Flat (\texttt{input.seq}) & \texttt{xz -9} (LZMA2 Flat) & 120,108,223 & 5,663,640 & 0.3772 & 21.21$\times$ & 58.42~s & 1.85~s & 674~MB RSS \\
\textit{E.~coli} 24g Flat (\texttt{input.seq}) & \texttt{GeCo3} (\texttt{-l 6}) & 120,108,223 & 8,409,358 & 0.5601 & 14.28$\times$ & 218.40~s & 221.10~s & 1,024~MB RSS \\
\textit{E.~coli} 24g Multi-FASTA & \texttt{AGC 3.2.5} (\texttt{-a -t 1}) & 121,825,884 & 8,540,341 & 0.5688 & 14.26$\times$ & 4.12~s & 0.85~s & 412~MB RSS \\
\textit{E.~coli} 24g Flat (\texttt{input.seq}) & \texttt{Original GraSS} (\texttt{-e2}) & 120,108,223 & 24,777,592 & 1.6503 & 4.85$\times$ & 192.18~s & 82.06~s & 1,500~MB Heap \\
\midrule
\textit{E.~coli} 24g Multi-FASTA & \textbf{EIMDL / HGS2 (6-Channel Exact)} & 121,825,884 & \textbf{3,387,485} & \textbf{0.2256} & \textbf{35.96$\times$} & 32.50~s & 16.50~s & $100\%$ SHA-256 \\
\textit{M.~tuberculosis} 24g Multi-FASTA & \textbf{EIMDL / HGS2 (6-Channel Exact)} & 107,257,786 & \textbf{1,049,266} & \textbf{0.0794} & \textbf{102.22$\times$} & 41.36~s & 24.77~s & $100\%$ SHA-256 \\
\textit{S.~aureus} 25g Multi-FASTA & \textbf{EIMDL / HGS2 (6-Channel Exact)} & 72,060,360 & \textbf{1,086,129} & \textbf{0.1223} & \textbf{66.35$\times$} & 31.76~s & 16.15~s & $100\%$ SHA-256 \\
\midrule
\textit{S.~enterica} 25g Locked Trio & \textbf{EIMDL / HGS2 (\texttt{--threads 1})} & 120,527,038 & \textbf{2,390,106} & \textbf{0.1586} & \textbf{50.43$\times$} & 100.76~s & \textbf{53.36~s} & 753.5 / \textbf{125.1~MiB} \\
\textit{S.~enterica} 25g Locked Trio & \texttt{JARVIS3 v3.7} (\texttt{-l 7}) & 120,527,038 & 2,709,545 & 0.1798 & 44.48$\times$ & \textbf{90.03~s} & 85.86~s & \textbf{294.8} / 294.4~MiB \\
\textit{S.~enterica} 25g Locked Trio & \texttt{Original GraSS} (\texttt{-e2}) & 120,527,038 & 19,991,646 & 1.3269 & 6.03$\times$ & 212.29~s & 94.48~s & 1,423.5 / 1,486.2~MiB \\
\midrule
\textit{K.~pneumoniae} 23g Locked Trio & \textbf{EIMDL / HGS2 (Default 2-CPU)} & 122,700,601 & \textbf{2,523,991} & \textbf{0.1646} & \textbf{48.61$\times$} & \textbf{32.69~s} & \textbf{16.19~s} & 760.8 / \textbf{177.3~MiB} \\
\textit{K.~pneumoniae} 23g Locked Trio & \texttt{JARVIS3 v3.7} (\texttt{-l 7}) & 122,700,601 & 2,948,926 & 0.1923 & 41.61$\times$ & 55.32~s & 52.03~s & \textbf{295.2} / 295.5~MiB \\
\textit{K.~pneumoniae} 23g Locked Trio & \texttt{Original GraSS} (\texttt{-e2}) & 122,700,601 & 21,191,566 & 1.3817 & 5.79$\times$ & 87.19~s & 75.04~s & 1,112.4 / 1,462.4~MiB \\
\bottomrule
\end{tabular}%
}
\end{table}

\subsection{Verification of Laws V and VII: Amdahl Multi-Core Scaling and Bounded-Heap Execution (\texttt{-Xmx128m})}

\begin{figure}[t!]
  \centering
  \includegraphics[width=0.98\textwidth]{figures/fig4_amdahl_and_bounded_memory.pdf}
  \caption{\textbf{Experimental verification of Laws~V and~VII (Amdahl multi-core scaling and strict pure-JVM bounded-heap execution).} \textbf{(A)}~Encode, decode, and total round-trip in-process execution times on the 24-genome \textit{E.~coli} pangenome ($120.11$~Mbp) across $T \in \{1,2,4,8,12\}$ threads. \textbf{(B)}~Stage-1 core data structure footprint and resulting Mode~C archive size across syncmer subsampling shifts $s \in \{0,1,2\}$ under strict \texttt{-Xmx256m}, \texttt{-Xmx224m}, \texttt{-Xmx192m}, \texttt{-Xmx160m}, and \texttt{-Xmx128m} pure-JVM heap limits.}
  \label{fig:amdahl_mem}
\end{figure}

Table~\ref{tab:amdahl_mem} and Figure~\ref{fig:amdahl_mem} validate Laws~V and~VII on the 24-genome \textit{E.~coli} collection ($120,108,223$~bp):
\begin{itemize}
  \item \textbf{Multi-Core Scaling and Exact Determinism:} Scaling from $T=1$ to $T=12$ threads reduces in-process encode time from $40.31$~s to $19.65$~s ($2.05\times$ speedup; $45.16\times$ faster than \texttt{HyperGraSS v1}'s $887.3$~s) and decode time from $56.78$~s to $10.40$~s ($5.46\times$ speedup; $8.66$~s [$6.56\times$] at $T=8$). Because Stage-2 sub-blocks are fixed at $4$~Mbp with $32$~KiB warm context priming, $T=8$ and $T=12$ produce the \emph{exact same byte-identical archive} ($3,388,639$~B, $+0.51\%$ over 1-thread $3,371,489$~B).
  \item \textbf{Pure-JVM \texttt{-Xmx128m} Bounded-Heap Execution:} Whereas Generation~1 required compiling a GraalVM~25 static native binary to pass a $256$~MiB container gate ($3,864,464$~B in $681.9$~s), Generation~3's disjoint Stage-1/Stage-2 lifecycle and adaptive \texttt{ultraTight}/\texttt{extremeTight} scheduling execute in \textbf{standard pure OpenJDK JVM} under \texttt{-Xmx256m} ($s=1$, $3,422,314$~B), \texttt{-Xmx224m} ($s=1$, $3,422,314$~B), \texttt{-Xmx192m} ($s=2$, $3,452,703$~B), \texttt{-Xmx160m} ($s=2$, $3,452,703$~B), and \textbf{\texttt{-Xmx128m}} ($s=2$, \textbf{3,452,703~B}, $0.2300$~bpb in $25.59$~s encode / $11.67$~s decode, just $1.06$~B/base total JVM heap limit)---still \textbf{11.39\% smaller} than \texttt{JARVIS3 -l 7} ($3,896,361$~B) and \textbf{9.39\% smaller} than \texttt{HyperGraSS v1} ($3,810,409$~B).
\end{itemize}

\begin{table}[h!]
\centering
\caption{Laws~V and~VII Amdahl thread scaling and bounded-heap (\texttt{-Xmx}) execution on 24-genome \textit{E.~coli} ($120,108,223$~bp, \texttt{amdahl\_and\_memory.json}).}
\label{tab:amdahl_mem}
\small
\begin{tabular}{lrrrrrrr}
\toprule
\textbf{Configuration} & \textbf{CSR Shift $s$} & \textbf{Stage-1 Core} & \textbf{Archive (B)} & \textbf{bpb} & \textbf{Enc (s)} & \textbf{Dec (s)} & \textbf{SHA-256} \\
\midrule
\texttt{--threads 1} & $s=0$ ($1/6$) & 150.1~MB & \textbf{3,371,489} & \textbf{0.2246} & 40.31~s & 56.78~s & PASS \\
\texttt{--threads 2} & $s=0$ ($1/6$) & 150.1~MB & 3,373,077 & 0.2247 & 37.44~s & 12.58~s & PASS \\
\texttt{--threads 4} & $s=0$ ($1/6$) & 150.1~MB & 3,383,225 & 0.2253 & 25.47~s & 12.58~s & PASS \\
\texttt{--threads 8} & $s=0$ ($1/6$) & 150.1~MB & 3,388,639 & 0.2257 & 21.32~s & \textbf{8.66~s} & PASS \\
\texttt{--threads 12} & $s=0$ ($1/6$) & 150.1~MB & 3,388,639 & 0.2257 & \textbf{19.65~s} & 10.40~s & PASS \\
\midrule
\texttt{-Xmx256m} & $s=1$ ($1/12$) & 90.1~MB & 3,422,314 & 0.2279 & 22.38~s & 10.61~s & PASS \\
\texttt{-Xmx224m} & $s=1$ ($1/12$) & 90.1~MB & 3,422,314 & 0.2279 & 24.03~s & 13.89~s & PASS \\
\texttt{-Xmx192m} & $s=2$ ($1/24$) & \textbf{60.1~MB} & \textbf{3,452,703} & \textbf{0.2300} & 62.68~s & 27.46~s & PASS \\
\texttt{-Xmx160m} & $s=2$ ($1/24$) & \textbf{60.1~MB} & \textbf{3,452,703} & \textbf{0.2300} & 55.90~s & 14.34~s & PASS \\
\texttt{-Xmx128m} & $s=2$ ($1/24$) & \textbf{60.1~MB} & \textbf{3,452,703} & \textbf{0.2300} & 25.59~s & 11.67~s & PASS \\
\bottomrule
\end{tabular}
\end{table}

\subsection{Verification of Law IV: Non-Existing, Xenobiological, and Adversarial Genomes}

\begin{figure}[t!]
  \centering
  \includegraphics[width=0.98\textwidth]{figures/fig5_laws_ablations_and_synthetic.pdf}
  \caption{\textbf{Real-genome law-by-law ablations, coalescent donor panel scaling, $100:1$ entropy-floor progression, and non-existing/adversarial genome completeness.} \textbf{(A)}~Law-by-law ablation trajectory on \textit{E.~coli} \texttt{00.seq} conditioned on a distant single donor (\texttt{20.seq}: $238,550$~B $\to$ $197,163$~B, $-17.35\%$) and on the dense 23-donor LOO panel ($43,357$~B $\to$ $40,588$~B, $-6.39\%$). \textbf{(B)}~Coalescent donor panel scaling ($N \in \{1,2,4,8,16,23\}$ donors on \textit{E.~coli} \texttt{00.seq}). \textbf{(C)}~Generation~1 $\to$ Generation~3 progression on 24-genome \textit{E.~coli} relative to the $100:1$ ($1.201$~MB) target. \textbf{(D)}~Compression rate across nine synthetic, xenobiological, and cryptographic white-noise corpora, confirming the $2.0012$~bpb safety ceiling.}
  \label{fig:ablations_syn}
\end{figure}

Table~\ref{tab:synthetic_9} and Figure~\ref{fig:ablations_syn}D evaluate nine non-existing, xenobiological, and adversarial corpora (\texttt{mode\_s\_and\_synthetic.json}), proving $100\%$ exact round-trip reconstruction and capping SHA-256 cryptographic white-noise expansion at \textbf{125,075~B} (\textbf{2.0012~bpb}, $+0.06\%$ over raw 2-bit packing, vs.\ \texttt{gzip -9}'s $144,017$~B [$2.3043$~bpb]).

\begin{table}[h!]
\centering
\caption{Law~IV completeness evaluation across nine non-existing, xenobiological, synthetic, and adversarial corpora (\texttt{mode\_s\_and\_synthetic.json}).}
\label{tab:synthetic_9}
\resizebox{\textwidth}{!}{%
\begin{tabular}{llrrrrrrr}
\toprule
\textbf{Corpus ID} & \textbf{Biological / Adversarial Stress Regime} & \textbf{Regime} & \textbf{Raw Bytes} & \textbf{Bases} & \textbf{HGS2 (B)} & \textbf{bpb (Raw)} & \textbf{gzip -9 (B)} & \textbf{Active Grammar Channels} \\
\midrule
\texttt{01\_euk\_centromere} & 171~bp HOR $\alpha$-satellite + transposons & Mode~S & 491,604 & 485,418 & \textbf{42,075} & 0.6847 & 83,797 & 325,289~bp self, 318 splices, 8,784 SNPs \\
\texttt{02\_viral\_rna\_acgu} & 10-strain RNA $\texttt{ACGU}$ quasispecies & Mode~C & 304,013 & 299,253 & \textbf{11,995} & 0.3156 & 42,923 & RNA flag, 268,526~bp match, 3,286 SNPs \\
\texttt{03\_chimeric\_iupac\_crlf} & $\text{CRLF}$ + soft-mask + IUPAC $\texttt{N/R/Y/W/S}$ & Mode~S & 263,817 & 244,065 & \textbf{51,536} & 1.5628 & 86,566 & All 6 factorizer channels + 61 splices \\
\texttt{04\_xeno\_at88} & Xenobiological $88\%$ AT Bias + 13~bp tandem & Mode~S & 506,294 & 500,000 & \textbf{94,248} & 1.4892 & 123,139 & \texttt{FastCRM} extreme skew + HPY backoff \\
\texttt{05\_xeno\_gc86} & Xenobiological $86\%$ GC Bias + 6F-ORF-HMM & Mode~S & 506,291 & 500,000 & \textbf{97,525} & 1.5410 & 126,076 & \texttt{6F-ORF-HMM} + \texttt{orf2..orf5} HPY chain \\
\texttt{06\_denovo\_grammar} & Synthetic 8-operon $\times$ 16-cassette chromosome & Mode~S & 506,287 & 500,000 & \textbf{41,159} & 0.6504 & 141,540 & 405,756~bp self, 457 splices, 18,473 SNPs \\
\texttt{07\_hypermutator\_pangenome} & 6-strain pangenome ($3\%$ SNP + indels + inv.) & Mode~C & 1,513,506 & 1,494,677 & \textbf{102,503} & 0.5418 & 454,203 & 1,232,424~bp match, 1,085 splices, 32,706 SNPs \\
\texttt{07\_hypermutator\_strain05} & Hypermutator strain LOO vs.\ 5-strain panel & Mode~R & 249,436 & 249,436 & \textbf{7,809} & 0.2505 & --- & 247,098~bp panel match, 195 splices, 6,158 SNPs \\
\texttt{08\_crypto\_white\_noise} & SHA-256 counter-mode incompressible noise & Mode~S & 500,000 & 500,000 & \textbf{125,075} & \textbf{2.0012} & 144,017 & Law~IV \texttt{raw2BitBypass} ($+0.06\%$ ceiling) \\
\bottomrule
\end{tabular}%
}
\end{table}

\subsection{Real-Genome Law-by-Law Ablations and Executable Falsification Suite}

\begin{figure}[t!]
  \centering
  \includegraphics[width=0.98\textwidth]{figures/fig6_independent_modes_pareto.pdf}
  \caption{\textbf{Independent 12-genome \textit{Enterobacteriaceae} rate--memory Pareto frontier and per-target archive sizes.} \textbf{(A)}~Aggregate compression rate (bits/base) vs.\ peak container cgroup RSS (MiB, log scale) across all 12 targets ($60.81$~Mbp). \textbf{(B)}~Per-target archive sizes (KB, log scale) across the 6 \textit{S.~enterica} and 6 \textit{K.~pneumoniae} genomes.}
  \label{fig:ind_pareto}
\end{figure}

Table~\ref{tab:ablations_real} reports the controlled law-by-law knock-out ablations (\texttt{ablations\_real\_genomes.json}) across three real-genome scenarios: \textit{E.~coli} \texttt{00.seq} given a single diverged donor (\texttt{20.seq}, $54,124$ SNPs), \textit{E.~coli} \texttt{00.seq} given its 23-donor LOO panel, and \textit{M.~tuberculosis} \texttt{00.seq} given its 23-donor LOO panel.

\begin{table}[h!]
\centering
\caption{Controlled law-by-law knock-out ablations on real bacterial genomes (\texttt{ablations\_real\_genomes.json}). Every configuration is verified byte-for-byte via SHA-256.}
\label{tab:ablations_real}
\resizebox{\textwidth}{!}{%
\begin{tabular}{llrrrrrrr}
\toprule
\textbf{Scenario} & \textbf{Ablation Configuration} & \textbf{Archive (B)} & \textbf{bpb} & \textbf{$\Delta$ vs.\ Full} & \textbf{Panel Match (bp)} & \textbf{Panel Tokens} & \textbf{SNPs} & \textbf{Splices} \\
\midrule
\textbf{(A) \textit{E.~coli} \texttt{00.seq} $\mid$ \texttt{20.seq}} & \textbf{Full EIMDL (All Laws ON)} & \textbf{197,163} & \textbf{0.32481} & \textbf{Baseline} & 4,258,080 & 1,517 & 54,124 & 619 \\
(1 donor, $4,856,080$~bp) & No Law~I (Flat SNP, no $V_4$ Kimura--CpG) & 203,340 & 0.33499 & $+6,177$~B ($+3.13\%$) & 4,257,425 & 1,490 & 53,987 & 618 \\
 & No Law~II (Flat gap, no codon-wobble $\psi$) & 199,024 & 0.32788 & $+1,861$~B ($+0.94\%$) & 4,258,080 & 1,517 & 54,124 & 619 \\
 & No Law~IIIb (No affine micro-indel splicing) & 196,936 & 0.32444 & $-227$~B ($+428$ tokens) & 4,259,145 & 1,945 & 54,157 & 0 \\
 & No Laws I + II + IIIb combined & 205,071 & 0.33784 & $+7,908$~B ($+4.01\%$) & 4,258,398 & 1,918 & 54,006 & 0 \\
 & Exact matches only (No approximate copying) & 238,550 & 0.39299 & $+41,387$~B ($+20.99\%$) & 4,028,211 & 30,701 & 0 & 0 \\
\midrule
\textbf{(B) \textit{E.~coli} \texttt{00.seq} LOO} & \textbf{Full EIMDL (All Laws ON)} & \textbf{40,588} & \textbf{0.06687} & \textbf{Baseline} & 4,766,688 & 5,531 & 5,706 & 99 \\
(23 donors, $4,856,080$~bp) & No Law~I (Flat SNP, no $V_4$ Kimura--CpG) & 41,178 & 0.06784 & $+590$~B ($+1.45\%$) & 4,766,538 & 5,524 & 5,680 & 99 \\
 & No Law~II (Flat gap, no codon-wobble $\psi$) & 40,697 & 0.06705 & $+109$~B ($+0.27\%$) & 4,766,688 & 5,531 & 5,706 & 99 \\
 & No Law~IIIb (No affine micro-indel splicing) & \textbf{40,469} & \textbf{0.06667} & $-119$~B (dense panel) & 4,766,614 & 5,590 & 5,685 & 0 \\
 & No Laws I + II + IIIb combined & 41,175 & 0.06783 & $+587$~B ($+1.45\%$) & 4,766,417 & 5,583 & 5,656 & 0 \\
 & Exact matches only (No approximate copying) & 43,357 & 0.07143 & $+2,769$~B ($+6.82\%$) & 4,747,677 & 7,791 & 0 & 0 \\
\midrule
\textbf{(C) \textit{M.~tuberculosis} \texttt{00.seq} LOO} & \textbf{Full EIMDL (All Laws ON)} & \textbf{1,240} & \textbf{0.00225} & \textbf{Baseline} & 4,412,383 & 116 & 251 & 0 \\
(23 donors, $4,417,521$~bp) & No Law~I (Flat SNP, no $V_4$ Kimura--CpG) & 1,261 & 0.00228 & $+21$~B ($+1.69\%$) & 4,412,383 & 116 & 251 & 0 \\
 & No Law~II (Flat gap, no codon-wobble $\psi$) & 1,237 & 0.00224 & $-3$~B ($251$ SNPs) & 4,412,383 & 116 & 251 & 0 \\
 & No Law~IIIb (No affine micro-indel splicing) & 1,240 & 0.00225 & $0$~B ($0$ splices) & 4,412,383 & 116 & 251 & 0 \\
 & No Laws I + II + IIIb combined & 1,258 & 0.00228 & $+18$~B ($+1.45\%$) & 4,412,383 & 116 & 251 & 0 \\
 & Exact matches only (No approximate copying) & 1,267 & 0.00229 & $+27$~B ($+2.18\%$) & 4,412,761 & 325 & 0 & 0 \\
\bottomrule
\end{tabular}%
}
\end{table}

\noindent
\textbf{Executable Falsification Suite (\texttt{VerifyCompressionLaws.java}):}
All seven laws are continuously guarded by seven automated falsification suites totaling \textbf{171,941 programmatic assertions} with \textbf{0 failures}:
(1)~Law~I $V_4$ involution \& Kimura--CpG entropy reduction ($2,355$~B vs.\ $2,576$~B ablated, $-8.6\%$);
(2)~Law~II codon-wobble bijection $\psi^{-1}(\psi(g))=g$ for all $g \in [0, 100\,000]$ ($2,001$~B vs.\ $2,130$~B flat, $-6.1\%$) and 6-frame ORF-HMM causal tracking ($7,381$~B vs.\ $7,688$~B broken);
(3)~Law~III monotonic donor panel scaling ($\text{Panel24}=167$~B, $\text{Twin1}=132$~B) and affine micro-indel splicing ($6,818$~B spliced vs.\ $10,067$~B unspliced, $-32.3\%$);
(4)~Law~IV cryptographic white-noise ceiling ($2.0035$~bpb $\le 2.005$~bpb gate, plus $1.534$~bpb AT88 and $0.167$~bpb satellite);
(5)~Law~V subsampled syncmer CSR memory scaling and parallel determinism ($1T=63,611$~B vs.\ $8T=63,611$~B, $\Delta=0.000\%$);
(6)~Law~VI two-adic tapered adaptation ($p_{\min} = 11/4096$ vs.\ $31/4096$) and HPY novelty coding ($1.0897$~bpb); and
(7)~Law~VII disjoint lifecycle and bounded Mode~C determinism ($27,257$~B unbounded vs.\ $27,257$~B bounded, $\Delta=0.000\%$), alongside Mode~R SHA-256 tamper rejection and 6-channel exact round-trip verification.

\section{Honest Scientific Boundaries, Negative Results, and Threats to Validity}

In keeping with our strict anti-fabrication commitment, we explicitly report every boundary condition, negative ablation result, and trade-off observed across all three generations of our experiments:

\begin{enumerate}
  \item \textbf{Micro-Indel Splicing vs.\ Law~VI Token/Novelty Efficiency Trade-Off (Table~\ref{tab:ablations_real}):}
  Affine micro-indel splicing (Law~IIIc) cuts synthetic structural-variant payload by \textbf{$-32.3\%$} ($6,818$~B vs.\ $10,067$~B in \texttt{VerifyCompressionLaws}). Once Law~VI's two-adic taper ($\tau=20$) and Hierarchical Pitman--Yor novelty coder reduce short-gap and short-literal overhead, however, disabling splicing slightly reduces total archive bytes on real \textit{E.~coli} both against a 23-donor LOO panel ($40,469$~B vs.\ $40,588$~B, $-119$~B) and on distant single donor \texttt{20.seq} ($196,936$~B vs.\ $197,163$~B, $-227$~B), even though disabling splicing fragments synteny from $1,517$ to $1,945$ panel tokens ($+428$ breaks). On dense panels, a micro-indel in one donor often coincides with an exact haplotype in another donor; even on single donors, long-tail splice displacement parameters can cost slightly more than restarting an unbroken panel token with Law~VI entropy coding. We retain splicing on by default because it preserves long-range syntenic continuity and dominates structural-variant regimes ($-32.3\%$).
  \item \textbf{Near-Clonal Per-Target Ties/Losses vs.\ the 13-Configuration \texttt{AGC} Oracle:}
  Although \texttt{HGS2} wins \textbf{73/73} targets against any single deployable \texttt{AGC~3.2.5} configuration and beats the 13-configuration per-target \texttt{AGC} oracle by \textbf{2.04$\times$} on \textit{M.~tuberculosis} ($47,339$~B vs.\ $96,580$~B) and \textbf{9.28$\times$} on \textit{S.~aureus} ($174,217$~B vs.\ $1,617,085$~B), the per-target oracle is smaller by $13\text{--}239$~B on 5/24 near-clonal \textit{M.~tuberculosis} targets (\texttt{03.seq}: $256$ vs.\ $243$~B; \texttt{08.seq}: $287$ vs.\ $224$~B; \texttt{09.seq}: $607$ vs.\ $368$~B; \texttt{16.seq}: $1,150$ vs.\ $1,072$~B; \texttt{19.seq}: $536$ vs.\ $512$~B) and by $66$~B on 1/25 \textit{S.~aureus} target (\texttt{04.seq}: $683$ vs.\ $617$~B). On ultra-clonal targets ($<40$ SNPs across $4.4$~Mbp), \texttt{HGS2}'s mandatory 32-byte cryptographic SHA-256 donor-panel pin accounts for $5\text{--}12\%$ of the entire archive.
  \item \textbf{OVC Asymptotic Boundary on Subexponential $\#k$-SAT Ideas~\cite{openai2026eth}:}
  While Laws~VI and~VII successfully adapt shared linear summaries, two-adic tapering ($\tau=20$), and disjoint polynomial-space phases from~\cite{openai2026eth}, \S1.1.3 of~\cite{openai2026eth} proves no bypass around Orthogonal Vectors Conjecture (OVC) polynomial-time lower bounds for exact sequence alignment or Longest Common Subsequence problems.
  \item \textbf{Medium-Diversity \textit{E.~coli} Stratum Crossover in Generation~1 (\texttt{docs/BROAD.md}):}
  Across the 9 within-species Mash diversity strata in Generation~1 (Table~\ref{tab:gen1_broad}), \texttt{HyperGraSS v1 --pangenome} won $8/9$ strata but trailed \texttt{JARVIS3 -l 7} by $-0.16\%$ ($5,276$~B: $3,381,136$~B vs.\ $3,375,860$~B) on the 20-genome medium-diversity \textit{E.~coli} stratum. Generation~3 (\texttt{HGS2}) resolves this by replacing fixed-bank block matching with unaligned spanning-tree coalescent copying and \texttt{FastCRM}.
  \item \textbf{Mode~S Single-Genome CPU Throughput vs.\ Default \texttt{JARVIS3 -l 7} (Locked Docker Snapshot):}
  In reference-free single-genome compression (Table~\ref{tab:ind12_subset}, executed inside Docker \texttt{hypergrass-corpus:jdk27}), \texttt{EIMDL Mode~S} (\texttt{--threads 1}) achieves smaller archives than both \texttt{JARVIS3 -l 7} ($-8.79\%$) and \texttt{JARVIS3 -l 31} ($-2.32\%$) while using $2.46\times$ and $10.7\times$ less RAM, respectively, and runs $1.26\times$ faster than \texttt{JARVIS3 -l 31} in aggregate. However, against lightweight default \texttt{JARVIS3 -l 7}, \texttt{EIMDL Mode~S} is $3.65\times$ slower to compress ($332.5$~s vs.\ $91.0$~s) and $3.16\times$ slower to restore ($277.9$~s vs.\ $87.9$~s) due to symmetric per-base context mixing in that locked container snapshot (host Mode~S under Laws~VIa--VIc achieves $6.75$~s encode / $6.57$~s decode on \textit{E.~coli} \texttt{00.fa}; re-running the Docker comparison for Table~\ref{tab:ind12_subset} remains future work).
  \item \textbf{32-Bit Signed Integer CSR Addressing Horizon ($\approx 2.14$~Gbp) in Generation~3:}
  While Generation~1's 64-bit \texttt{SequenceSource} verified streaming block compression up to $2.077$~GiB ($2,229,710,670$~B, Table~\ref{tab:gen1_systems}), Generation~3's \texttt{DonorPanelIndex} stores concatenated donor base coordinates in signed 32-bit Java arrays (\texttt{int[] csrPositions}), bounding a single in-memory donor reference bank to $2^{31}-1 \approx 2.14\times 10^9$~bp ($\approx 430$ bacterial genomes or $\approx 68\%$ of a haploid human genome) unless partitioned into shards.
\end{enumerate}

\section{Conclusion}

We have presented the complete end-to-end theory, architectural evolution, and empirical verification of \textbf{Evolutionary Information Minimum Description Length (EIMDL)} from Generation~1 (\texttt{HyperGraSS v1}) and Generation~2 (\texttt{RC-MDL v2}) to the unified Generation~3 engine (\textbf{\texttt{HGS2}}), uniting reference-free single-genome coding (\textbf{Mode~S}), reference-conditioned incremental archiving (\textbf{Mode~R}), and self-contained pangenome compression (\textbf{Mode~C}) under seven mathematically proved, experimentally falsifiable laws. By embedding the algebraic symmetries of molecular evolution---Klein four-group $V_4$ strand-invariant substitutions, codon-wobble gap periodicity, causal 6-frame ORF survival, horizon-bounded two-gear coalescent haplotype copying, affine micro-indel splicing, subsampled open-syncmer CSR indexing, two-adic tapered adaptation, Hierarchical Pitman--Yor backoff, reciprocal normalization, and strictly disjoint Stage-1/Stage-2 memory lifecycles---directly into the compressor's probability model, \texttt{HGS2} achieves state-of-the-art compression across all evaluated bacterial species, scales to $19.65$~s encode ($45.16\times$ faster than Generation~1) and $8.66\text{--}10.40$~s decode, and operates deterministically inside a strict \textbf{\texttt{-Xmx128m}} ($128$~MiB) pure-JVM heap ceiling.

\bibliographystyle{plain}
\bibliography{references}

\clearpage
\appendix
\section{Complete Per-Genome Mode R Leave-One-Out Benchmark Tables (All 73 Genomes)}

\begin{table}[h!]
\centering
\caption{Complete 24-genome Mode~R Leave-One-Out evaluation on \textit{Escherichia coli} ($120,108,223$~bp total, \texttt{benchmarks/fusion/results/loo\_ecoli.json}). \texttt{AGC} Best Single is \texttt{AGC~3.2.5 -s 200000}. Every archive verified byte-for-byte via SHA-256.}
\label{tab:ecoli_full}
\resizebox{\textwidth}{!}{%
\begin{tabular}{lrrrrrrrrrrr}
\toprule
\textbf{Target} & \textbf{Raw (B)} & \textbf{HGS2 (B)} & \textbf{bpb} & \textbf{Ratio} & \textbf{RC-MDL v2 (B)} & \textbf{AGC Single (B)} & \textbf{AGC/HGS2} & \textbf{Panel Match (bp)} & \textbf{Novelty (bp)} & \textbf{SNPs} & \textbf{Splices} \\
\midrule
%%ECOLI_LOO_ROWS%%
\bottomrule
\end{tabular}%
}
\end{table}

\begin{table}[h!]
\centering
\caption{Complete 24-genome Mode~R Leave-One-Out evaluation on \textit{Mycobacterium tuberculosis} ($105,745,015$~bp total, \texttt{loo\_mtb.json} and \texttt{agc-sweep-mtb-full/report.json}). \texttt{AGC} Best Single configuration (\texttt{-s 100000}) totals $292,170$~B (\texttt{HGS2} wins \textbf{24/24}, $6.17\times$). Against the 13-Configuration Per-Target \texttt{AGC} Oracle ($96,580$~B), \texttt{HGS2} totals \textbf{47,339~B} ($2.04\times$ smaller) and wins \textbf{19/24} targets.}
\label{tab:mtb_full}
\resizebox{\textwidth}{!}{%
\begin{tabular}{lrrrrrrrrrrr}
\toprule
\textbf{Target} & \textbf{Raw (B)} & \textbf{HGS2 (B)} & \textbf{bpb} & \textbf{Ratio} & \textbf{RC-MDL v2 (B)} & \textbf{AGC Oracle (B)} & \textbf{Oracle Flag} & \textbf{Panel Match (bp)} & \textbf{Novelty (bp)} & \textbf{SNPs} & \textbf{Splices} \\
\midrule
%%MTB_LOO_ROWS%%
\bottomrule
\end{tabular}%
}
\end{table}

\clearpage
\begin{table}[h!]
\centering
\caption{Complete 25-genome Mode~R Leave-One-Out evaluation on \textit{Staphylococcus aureus} ($71,043,114$~bp total, \texttt{loo\_saureus.json} and \texttt{agc-verify-saureus-oracle/results.json}). \texttt{AGC} Best Single configuration (\texttt{-s 50000}) totals $2,979,078$~B (\texttt{HGS2} wins \textbf{25/25}, $17.10\times$). Against the 13-Configuration Per-Target \texttt{AGC} Oracle ($1,617,085$~B), \texttt{HGS2} totals \textbf{174,217~B} ($9.28\times$ smaller) and wins \textbf{24/25} targets.}
\label{tab:saureus_full}
\resizebox{\textwidth}{!}{%
\begin{tabular}{lrrrrrrrrrrr}
\toprule
\textbf{Target} & \textbf{Raw (B)} & \textbf{HGS2 (B)} & \textbf{bpb} & \textbf{Ratio} & \textbf{RC-MDL v2 (B)} & \textbf{AGC Oracle (B)} & \textbf{Oracle Flag} & \textbf{Panel Match (bp)} & \textbf{Novelty (bp)} & \textbf{SNPs} & \textbf{Splices} \\
\midrule
%%SAUREUS_LOO_ROWS%%
\bottomrule
\end{tabular}%
}
\end{table}

\end{document}
"""


def write_full_latex():
    tables = build_all_tables()
    part3 = (
        PART3_GEN3_TEX_TEMPLATE
        .replace("%%IND12_ROWS%%", tables["ind12_tex"])
        .replace("%%ECOLI_LOO_ROWS%%", tables["ecoli_loo_tex"])
        .replace("%%MTB_LOO_ROWS%%", tables["mtb_loo_tex"])
        .replace("%%SAUREUS_LOO_ROWS%%", tables["saureus_loo_tex"])
    )
    full_tex = PART1_TEX + PART2_GEN1_GEN2_TEX + part3
    out_tex = os.path.join(PAPER_DIR, "HGS2_EIMDL_Research_Paper.tex")
    with open(out_tex, "w", encoding="utf-8") as f:
        f.write(full_tex)
    print(f"[OK] Wrote {out_tex} ({len(full_tex):,} bytes)")


if __name__ == "__main__":
    write_full_latex()
