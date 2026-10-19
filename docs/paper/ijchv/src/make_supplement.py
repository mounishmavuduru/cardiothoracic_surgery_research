import sys
src, dst = sys.argv[1], sys.argv[2]
t = open(src, encoding="utf-8").read()
subs = [
    ("\\graphicspath{{figures/}}", "\\graphicspath{{../jama/figures/}}"),
    ("pdftitle={Supplement 2. Full technical report", "pdftitle={Supplementary material. Full technical report"),
    ("pdfsubject={Online-only supplement to the JAMA Cardiology Original Investigation}",
     "pdfsubject={Supplementary material to the IJC Heart \\& Vasculature Original Article}"),
    ("{\\Large\\bfseries Supplement 2. Full technical report\\par}",
     "{\\Large\\bfseries Supplementary material. Full technical report\\par}"),
    ("technical report behind the Original\nInvestigation.", "technical report behind the main\narticle."),
]
for old, new in subs:
    n = t.count(old)
    if n != 1:
        sys.exit(f"ABORTED -- {n} matches for {old[:50]!r}")
    t = t.replace(old, new)
open(dst, "w", encoding="utf-8", newline="").write(t)
print("ok")
