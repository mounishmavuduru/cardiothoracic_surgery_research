# IJC Heart & Vasculature submission

Built from the verified JAMA Cardiology article (`../jama/`), which stays the
single source of the main text, tables, figures and bibliography.

## Build

```
python docs/paper/ijchv/src/assemble_ijchv.py docs/paper/ijchv/IJCHV_manuscript.tex
python docs/paper/ijchv/src/make_supplement.py docs/paper/jama/Supplement2_full_technical_report.tex docs/paper/ijchv/IJCHV_Supplementary_material.tex
cd docs/paper/ijchv && tectonic IJCHV_manuscript.tex && tectonic IJCHV_Supplementary_material.tex
python docs/paper/jama/src/verify_oi.py docs/paper/ijchv/IJCHV_manuscript.tex docs/paper/manuscript.tex
```

`assemble_ijchv.py` runs the JAMA assembler and then changes only the title page
(`src/front.tex`), the abstract (`src/abstract.tex`, at most 250 words under
Background/Methods/Results/Conclusions), the back matter (`src/back.tex`) and the
citation style (numbered square brackets placed before punctuation). Each
substitution must match the expected number of times or nothing is written.

## Differences from the JAMA version

* No Key Points box. Six keywords instead of ten.
* The title page carries the IJC author-responsibility footnote.
* The back matter has CRediT, competing-interest, funding, ethics, data-availability
  and generative-AI declarations.
* Line numbers and 1.5 spacing, for review.
* The rho ranges `[1082,7384]` and `[817,16554]` are written with dashes, so they
  cannot be read as citation numbers.
* The TRIPOD+AI checklist (JAMA Supplement 1) is left out. That file contains
  garbled location stubs and audit-style comments, and it needs rewriting before
  it goes to any journal.

## Upload set

| file | role |
|---|---|
| `IJCHV_manuscript.pdf` | manuscript |
| `IJCHV_Supplementary_material.pdf` | supplementary material |
| `Highlights.txt` | highlights (3--5 bullets, at most 85 characters each) |
| `IJCHV_Declaration_of_interest.pdf` | declaration of interest |
| `IJCHV_Author_agreement.pdf` | author agreement; sign it only once no other journal is considering the paper |
| figures in `figures/` | if the system asks for separate figure files |
