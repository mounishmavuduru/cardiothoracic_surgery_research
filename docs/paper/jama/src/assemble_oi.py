"""Assemble the condensed JAMA Cardiology Original Investigation.

Sections come verbatim from the drafted fragments; citations are injected by
exact string match and the script aborts without writing if any pattern does not
match exactly once, so a stale fragment cannot silently drop a reference.
"""
import os
import re
import sys

SP = os.path.dirname(os.path.abspath(__file__))
DST = sys.argv[1]


def read(name):
    return open(os.path.join(SP, name), encoding="utf-8-sig").read().replace("﻿", "")


parts = read("jama_parts.tex")
preamble_title = parts.split("%% ===================== FLOATS =====================")[0]
floats = parts.split("%% ===================== FLOATS =====================")[1]
floats, artinfo = floats.split("%% ===================== ARTICLE INFORMATION =====================")

keypoints = read("keypoints.tex")
abstract = read("abstract.tex")
intro = read("intro.tex")
methods = read("methods.tex")
results = read("results.tex")
discussion = read("discussion.tex")
bib = read("bib.tex")

# (target-name, old, new, expected_count, why)
CITES = [
    ("methods",
     "is reported in accordance with TRIPOD+AI, and a completed checklist",
     "is reported in accordance with TRIPOD+AI~\\cite{tripodai2024}, and a completed checklist",
     1, "reporting guideline"),
    ("methods",
     "The Roney virtual cohort (Zenodo 5801337)",
     "The Roney virtual cohort~\\cite{roney2022} (Zenodo 5801337)",
     1, "cohort 1"),
    ("methods",
     "The UW/Boyle late-gadolinium-enhancement cohort (Dryad",
     "The UW/Boyle late-gadolinium-enhancement cohort~\\cite{uwboyle2025} (Dryad",
     1, "cohort 2"),
    ("methods",
     "mapped from the shipped\nimage-intensity ratio",
     "mapped from the shipped\nimage-intensity ratio~\\cite{khurram2014}",
     1, "IIR definition"),
    ("methods",
     "$\\partial\\lam_2/\\partial w_{ij}=(\\phii_{2,i}-\\phii_{2,j})^{2},",
     "$\\partial\\lam_2/\\partial w_{ij}=(\\phii_{2,i}-\\phii_{2,j})^{2},",
     0, "no-op placeholder"),
    ("methods",
     "was a monodomain Mitchell--Schaeffer inducibility verdict",
     "was a monodomain Mitchell--Schaeffer~\\cite{mitchell2003} inducibility verdict",
     1, "solver"),
    ("methods",
     "compared by the paired DeLong test in its fast midrank form",
     "compared by the paired DeLong test~\\cite{delong1988} in its fast midrank"
     " form~\\cite{sunxu2014}",
     1, "test + implementation"),
    ("results",
     "and we flag it as anomalous",
     "and we flag it as anomalous~\\cite{kumar2012,oral2008,marquardt2018,kawai2019,liu2020}",
     1, "inducibility literature"),
    ("results",
     "$\\rhoval=\\lVert\\Delta L\\rVert/(\\lam_3-\\lam_2)$ is $O(1)$",
     "$\\rhoval=\\lVert\\Delta L\\rVert/(\\lam_3-\\lam_2)$, the Davis--Kahan"
     " ratio~\\cite{daviskahan1970}, is $O(1)$",
     1, "validity radius provenance"),
    ("results",
     "within $0.5\\%$ of the Weyl\n2-manifold exponent",
     "within $0.5\\%$ of the Weyl\n2-manifold exponent~\\cite{weyl1912}",
     1, "Weyl law"),
    ("intro",
     "(Fiedler,\n1973)~\\cite{fiedler1973}",
     "(Fiedler,\n1973)~\\cite{fiedler1973}",
     0, "already cited"),
]

named = {"intro": intro, "methods": methods, "results": results,
         "discussion": discussion, "abstract": abstract}
failed = []
for tgt, old, new, want, why in CITES:
    if want == 0:
        continue
    got = named[tgt].count(old)
    if got != want:
        failed.append((tgt, why, want, got))
        continue
    named[tgt] = named[tgt].replace(old, new)

if failed:
    print("ABORTED -- nothing written. Citation pattern mismatches:")
    for tgt, why, want, got in failed:
        print(f"  {tgt}: want {want} got {got} | {why}")
    sys.exit(1)

doc = "".join([
    preamble_title.rstrip(), "\n\n",
    "%% ===================== KEY POINTS =====================\n",
    keypoints.strip(), "\n\n\\newpage\n\n",
    "%% ===================== ABSTRACT =====================\n",
    "\\section*{Abstract}\n",
    named["abstract"].strip(), "\n\n\\newpage\n\n",
    "%% ===================== MAIN TEXT =====================\n",
    named["intro"].strip(), "\n\n",
    named["methods"].strip(), "\n\n",
    named["results"].strip(), "\n\n",
    named["discussion"].strip(), "\n\n",
    "%% ===================== ARTICLE INFORMATION =====================",
    artinfo.rstrip(), "\n\n",
    "%% ===================== REFERENCES =====================\n",
    "\\begingroup\\small\n", bib.strip(), "\n\\endgroup\n\n",
    "%% ===================== FLOATS =====================\n",
    "\\newpage\n", floats.strip(), "\n\n",
    "\\end{document}\n",
])

# ---- word counts (JAMA counts main text only: Introduction..Conclusions) ----


def wc(tex):
    t = re.sub(r"(?<!\\)%.*", "", tex)
    t = re.sub(r"\\(begin|end)\{[^}]*\}", " ", t)
    t = re.sub(r"\\cite\{[^}]*\}", " ", t)
    t = re.sub(r"\\code\{([^}]*)\}", r"\1", t)
    t = re.sub(r"\$[^$]*\$", " X ", t)          # each math atom counts as one word
    t = re.sub(r"\\[a-zA-Z]+\*?(\[[^\]]*\])?", " ", t)
    t = re.sub(r"[{}~&\\]", " ", t)
    return len([w for w in t.split() if re.search(r"[A-Za-z0-9]", w)])


main = "\n".join(named[k] for k in ("intro", "methods", "results", "discussion"))
counts = {
    "main text": wc(main),
    "abstract": wc(named["abstract"]),
    "key points": wc(keypoints.replace(r"\noindent\textbf{Key Points}", "")),
}
nrefs = bib.count("\\bibitem")
subs = {"WCMAIN": str(counts["main text"]), "WCABS": str(counts["abstract"]),
        "WCKP": str(counts["key points"]), "NREFS": str(nrefs)}
for k, v in subs.items():
    if doc.count(k) != 1:
        print(f"ABORTED -- placeholder {k} appears {doc.count(k)} times")
        sys.exit(1)
    doc = doc.replace(k, v)

open(DST, "w", encoding="utf-8", newline="").write(doc)

for k, v in counts.items():
    print(f"{k:12s} {v}")
print("references  ", nrefs)
print("floats      ", floats.count("\\begin{table}") + floats.count("\\begin{figure}"))
print("wrote", DST)
