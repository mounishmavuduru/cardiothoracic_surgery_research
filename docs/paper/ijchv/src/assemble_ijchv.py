"""Assemble the IJC Heart & Vasculature version from the verified JAMA build.

The main text, tables, figures and bibliography are not re-drafted: the JAMA
assembler is run first and its output is converted. Only the title page, the
abstract (250 words, Background/Methods/Results/Conclusions), the back matter
and the citation style change. Every substitution must match the expected
number of times or the script aborts without writing.

    python docs/paper/ijchv/src/assemble_ijchv.py docs/paper/ijchv/IJCHV_manuscript.tex
"""
import os
import re
import subprocess
import sys
import tempfile

SP = os.path.dirname(os.path.abspath(__file__))
JAMA_ASSEMBLER = os.path.join(SP, "..", "..", "jama", "src", "assemble_oi.py")
DST = sys.argv[1]


def read(name):
    return open(os.path.join(SP, name), encoding="utf-8-sig").read()


def abort(msg):
    print("ABORTED -- nothing written.", msg)
    sys.exit(1)


tmp = os.path.join(tempfile.mkdtemp(), "oi.tex")
r = subprocess.run([sys.executable, JAMA_ASSEMBLER, tmp], capture_output=True, text=True)
if r.returncode != 0:
    abort("JAMA assembler failed:\n" + r.stdout + r.stderr)
doc = open(tmp, encoding="utf-8").read()


def between(text, start, end):
    i, j = text.find(start), text.find(end)
    if i < 0 or j < 0 or j < i:
        abort(f"markers not found in order: {start!r} .. {end!r}")
    return i, j


def swap(text, old, new, want=1):
    got = text.count(old)
    if got != want:
        abort(f"want {want} got {got}: {old[:60]!r}")
    return text.replace(old, new)


M = "%% ===================== {} ====================="

# title page (drops the JAMA Key Points box with it)
i, j = between(doc, M.format("TITLE PAGE"), M.format("ABSTRACT"))
doc = doc[:i] + read("front.tex").strip() + "\n\n" + doc[j:]

# abstract
i, j = between(doc, M.format("ABSTRACT"), M.format("MAIN TEXT"))
abstract = read("abstract.tex").strip()
doc = (doc[:i] + M.format("ABSTRACT") + "\n\\section*{Abstract}\n" + abstract
       + "\n\n\\newpage\n\n" + doc[j:])

# back matter
i, j = between(doc, M.format("ARTICLE INFORMATION"), M.format("REFERENCES"))
doc = doc[:i] + M.format("BACK MATTER") + "\n" + read("back.tex").strip() + "\n\n" + doc[j:]

# preamble: metadata, numbered square-bracket citations, line numbers for review
doc = swap(doc, "pdfsubject={Original Investigation submitted to JAMA Cardiology}",
           "pdfsubject={Original Article submitted to IJC Heart \\& Vasculature}")
doc = swap(doc, """  pdfkeywords={atrial fibrillation; reentry; biomarker; prediction model;
               spectral graph theory; algebraic connectivity; negative result;
               pre-registration; TRIPOD+AI}""",
           """  pdfkeywords={atrial fibrillation; reentry; spectral graph theory;
               prediction model; negative result; pre-registration}""")
doc = swap(doc, """% AMA citation style: superscript numerals, no brackets; bibliography numbered "1."
\\makeatletter
\\renewcommand\\@cite[2]{\\textsuperscript{#1\\if@tempswa , #2\\fi}}
\\renewcommand\\@biblabel[1]{#1.}
\\makeatother
""", """% numbered citations in square brackets, placed before punctuation
\\usepackage[numbers,sort&compress]{natbib}
\\usepackage{lineno}
\\linenumbers
\\onehalfspacing
""")

# JAMA supplement naming; the TRIPOD+AI checklist is not part of this package
doc = swap(doc, "Full methods appear in eMethods in Supplement 2.",
           "Full methods appear in the Supplementary material.")
doc = swap(doc, ", and a completed checklist accompanies\nthis article.", ".")

# bracketed integer ranges read as citation numbers in a bracket-citation journal
doc = swap(doc, "(interquartile range $[1082,7384]$; full range\n$[817,16554]$)",
           "(interquartile range $1082$--$7384$; full range\n$817$--$16554$)")

# AMA puts the superscript after punctuation; bracketed numerals go before it
doc, n1 = re.subn(r"([.,;:])\\cite\{([^}]*)\}", r"~\\cite{\2}\1", doc)
doc, n2 = re.subn(r"(?<=[^\s~(])\\cite\{", r"~\\cite{", doc)
if re.search(r"[.,;:]~?\\cite\{", doc):
    abort("a citation still follows punctuation")


def wc(tex):
    t = re.sub(r"(?<!\\)%.*", "", tex)
    t = re.sub(r"\\(begin|end)\{[^}]*\}", " ", t)
    t = re.sub(r"\\cite\{[^}]*\}", " ", t)
    t = re.sub(r"\\code\{([^}]*)\}", r"\1", t)
    t = re.sub(r"\$[^$]*\$", " X ", t)
    t = re.sub(r"\\[a-zA-Z]+\*?(\[[^\]]*\])?", " ", t)
    t = re.sub(r"[{}~&\\]", " ", t)
    return len([w for w in t.split() if re.search(r"[A-Za-z0-9]", w)])


i, j = between(doc, M.format("MAIN TEXT"), M.format("BACK MATTER"))
counts = {"main text": wc(doc[i:j]), "abstract": wc(abstract)}
if counts["abstract"] > 250:
    abort(f"abstract is {counts['abstract']} words; the journal allows 250")
nrefs = doc.count("\\bibitem{")
for k, v in {"WCMAIN": counts["main text"], "WCABS": counts["abstract"], "NREFS": nrefs}.items():
    doc = swap(doc, k, str(v))
for leftover in ("Key Points", "JAMA", "Article Information"):
    if leftover in doc.split("\\begin{document}")[1]:
        abort(f"JAMA-specific text left in the body: {leftover!r}")

open(DST, "w", encoding="utf-8", newline="").write(doc)
for k, v in counts.items():
    print(f"{k:12s} {v}")
print("references  ", nrefs)
print(f"citations moved before punctuation: {n1}; spacing added: {n2}")
print("wrote", DST)
