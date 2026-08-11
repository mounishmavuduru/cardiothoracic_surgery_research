"""Build Supplement 2 (full technical report) from the audited master manuscript.

The condensed Original Investigation cites this for everything that could not fit
in 3000 words. Body is byte-identical to the master; only the title block is
replaced with a supplement banner, and the table of contents is kept because a
long supplement benefits from one.
"""
import sys

SRC, DST = sys.argv[1], sys.argv[2]
t = open(SRC, encoding="utf-8").read()

BANNER = r"""\hypersetup{
  pdftitle={Supplement 2. Full technical report --- AtrialSpectralBench},
  pdfauthor={Mounish Mavuduru},
  pdfsubject={Online-only supplement to the JAMA Cardiology Original Investigation}
}

\begin{document}
\thispagestyle{empty}

\begin{center}
{\Large\bfseries Supplement 2. Full technical report\par}
\vspace{6pt}
{\large Incremental predictive value and perturbation-theoretic validity of a
spectral fragility index for atrial reentry\par}
\vspace{4pt}
Mounish Mavuduru
\end{center}

\vspace{1em}
\noindent This supplement is the complete technical report behind the Original
Investigation. It contains the full methods, the substrate provenance audit and
its four controls, the independent-solver cross-check, the resolution-convergence
study, the cross-medium localization results, the power analyses, the formal
(proof-assistant) verification of the algebraic core, and the dated record of
every correction made during this work. Every number in the main article appears
here with its derivation and its stored source file.

\noindent Sections are numbered independently of the main article. References to
\code{results/*.json} are files in the public project repository.

\vspace{1em}
\hrule
\vspace{1em}
"""

start = t.index(r"\title{")
end = t.index(r"\maketitle") + len(r"\maketitle")
t = t[:start] + BANNER.strip() + t[end:]

open(DST, "w", encoding="utf-8", newline="").write(t)
print("wrote", DST)
print("banner in place:", "Supplement 2. Full technical report" in t)
print("body intact:", r"\section{Introduction}" in t)
print("toc kept:", r"\tableofcontents" in t)
