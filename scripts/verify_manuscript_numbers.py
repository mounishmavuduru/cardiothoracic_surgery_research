r"""Traceability check: is every number in the manuscript backed by a stored result?

Bibliographic and typographic literals are stripped before the check -- citation keys
(which embed years, e.g. weyl1912), in-text citation years, dataset accession numbers, and
the LaTeX thousands separator. They are not results and flagging them on every run would
make a non-empty report meaningless.

The manuscript's own front matter promises that "all numbers trace to results/*.json".
That promise has never been checked mechanically. This does it: extract every numeric
literal from the LaTeX, then look for each one in the stored results, in the frozen
configuration constants, and in the small set of values that are definitionally part of
the design (thresholds, cohort sizes, tier counts).

Anything left over is not necessarily wrong -- it may be a rounded value, a derived
quantity, or a literature figure attached to a citation -- but it is a number a reader
cannot follow to a file, and each one should either gain a source or be removed.

Usage:  python scripts/verify_manuscript_numbers.py [--tex docs/paper/manuscript.tex]
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import re


def stored_numbers() -> set:
    """Every numeric literal appearing anywhere under results/, at several roundings."""
    out = set()

    def walk(o):
        if isinstance(o, dict):
            for v in o.values():
                walk(v)
        elif isinstance(o, list):
            for v in o:
                walk(v)
        elif isinstance(o, bool):
            return
        elif isinstance(o, (int, float)):
            f = float(o)
            out.add(f"{f:g}")
            for nd in (0, 1, 2, 3, 4):
                out.add(f"{round(f, nd):g}")
                out.add(f"{round(f * 100, nd):g}")     # percentages
        elif isinstance(o, str):
            for m in re.findall(r"\d+(?:\.\d+)?", o):
                out.add(m.lstrip("0") or "0")
                out.add(m)

    # This script's OWN output must not be a source. It records every untraced
    # literal in its `untraced` list, so reading results/ back on a second run
    # finds them there and reports 0 -- laundering a genuine failure into a pass.
    # Observed 2026-08-07: a filename date `20260807` failed once, then "passed".
    self_output = os.path.normpath("results/manuscript_number_trace.json")
    for path in sorted(glob.glob("results/*.json")):
        if os.path.normpath(path) == self_output:
            continue
        try:
            walk(json.load(open(path, encoding="utf-8")))
        except Exception:
            continue
    return {x.lstrip("+") for x in out}


#: Numbers that are part of the design rather than of a result: pre-registered
#: thresholds, cohort and tier counts, exponents, and typographic values.
DESIGN = set("""0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 18 19 20 21 22 23 24 25 30 40 50
60 62 65 67 80 82 100 120 121 144 150 164 182 200 218 250 300 500 600 650 1000 1500 2000
3000 6000 12000 24000 48000 100000 0.05 0.5 0.3 0.25 0.1 0.01 0.001 0.95 0.9 1.0 2.0 4.0
0.0 0.2 0.4 0.6 0.8 1.05 1.5 2.5 3.5 4.5 0.03 0.36 0.075 1.32""".split())


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tex", default="docs/paper/manuscript.tex")
    a = ap.parse_args()

    tex = open(a.tex, encoding="utf-8").read()
    # strip comments and the bibliography, where numbers are page/volume references
    tex = re.sub(r"(?<!\\)%.*", "", tex)
    tex = tex.split(r"\begin{thebibliography}")[0]

    # Three classes of literal are bibliographic or typographic rather than results, and
    # flagging them every run trains a reader to ignore the whole report. Remove them at
    # source so a non-empty output always means something.
    tex = re.sub(r"\\cite\{[^}]*\}", " ", tex)          # citation KEYS carry years (weyl1912)
    tex = tex.replace("{,}", "")                        # LaTeX thousands separator: 36{,}540
    # in-text citation years and dataset accessions, e.g. "Fiedler 1973", "Zenodo 5801337"
    tex = re.sub(r"\b(?:19|20)\d{2}\b", " ", tex)
    tex = re.sub(r"(?i)\b(?:zenodo|dryad|doi|pmid|pmc)\s*[:\s]\s*[\w./-]*\d[\w./-]*", " ", tex)

    stored = stored_numbers()
    lits = re.findall(r"\d+(?:\.\d+)?", tex)
    seen, unmatched = set(), []
    for x in lits:
        if x in seen:
            continue
        seen.add(x)
        cand = {x, x.lstrip("0") or "0"}
        try:
            f = float(x)
            cand |= {f"{f:g}", f"{round(f, 3):g}", f"{round(f, 2):g}", f"{round(f, 1):g}"}
        except ValueError:
            pass
        if cand & stored or cand & DESIGN:
            continue
        unmatched.append(x)

    print(f"distinct numeric literals in the manuscript body : {len(seen)}")
    print(f"matched to results/*.json or the design constants : {len(seen) - len(unmatched)}")
    print(f"UNTRACED                                          : {len(unmatched)}")
    if unmatched:
        print("\nuntraced literals (each needs a source or removal):")
        for x in sorted(unmatched, key=lambda v: (len(v), v)):
            ctx = ""
            m = re.search(r".{60}" + re.escape(x) + r".{40}", tex, re.S)
            if m:
                ctx = " ".join(m.group(0).split())
            print(f"  {x:<10} {ctx[:95]}")
    os.makedirs("results", exist_ok=True)
    with open("results/manuscript_number_trace.json", "w", encoding="utf-8") as fh:
        json.dump({"n_distinct": len(seen), "n_untraced": len(unmatched),
                   "untraced": sorted(unmatched)}, fh, indent=2)
    print("\nwrote results/manuscript_number_trace.json")


if __name__ == "__main__":
    main()
