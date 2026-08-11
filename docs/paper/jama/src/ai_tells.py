"""Scan the article body against the Wikipedia 'Signs of AI writing' checklist,
restricted to the items that can apply to a journal manuscript. Reports counts
and context; changes nothing."""
import re
import sys
from collections import Counter

t = open(sys.argv[1], encoding="utf-8").read()
body = t.split("%% ===================== KEY POINTS =====================")[1]
body = body.split("%% ===================== REFERENCES =====================")[0]
body += t.split("%% ===================== FLOATS =====================")[-1]
body = re.sub(r"(?<!\\)%.*", "", body)
flat = " ".join(re.sub(r"\$[^$]*\$", " X ", body).split())

VOCAB = ["additionally", "align with", "boasts", "bolstered", "crucial",
         "deep dive", "delve", "emphasizing", "emphasising", "enduring",
         "enhance", "fostering", "garner", "highlight", "interplay",
         "intricate", "intricacies", "landscape", "meticulous", "pivotal",
         "robust", "showcase", "tapestry", "testament", "underscore",
         "valuable", "vibrant", "groundbreaking", "renowned", "diverse array",
         "profound", "rich ", "commitment to", "notably", "importantly"]
COPULA = ["serves as", "stands as", "functions as", "operates as",
          "represents a", "marks a", "refers to", "boasts", "features a"]
SIGNIF = ["testament", "underscores", "highlights the", "reflects broader",
          "plays a key role", "plays a crucial role", "setting the stage",
          "turning point", "indelible", "deeply rooted", "ongoing legacy",
          "paradigm", "sheds light"]
SUPERFICIAL = ["highlighting", "underscoring", "emphasizing", "emphasising",
               "ensuring", "reflecting", "symbolizing", "symbolising",
               "contributing to", "cultivating", "encompassing", "enhancing",
               "resonate with", "valuable insight"]
VAGUE = ["experts argue", "some critics", "observers have", "industry reports",
         "it is widely", "it is generally", "many researchers", "studies show"]
NEGPAR = [(r"\bnot just\b[^.;]{0,50}\bbut\b", "not just X but Y"),
          (r"\bnot only\b[^.;]{0,50}\bbut\b", "not only X but Y"),
          (r"\brather than\b", "X rather than Y"),
          (r"\bis not\b[^.;]{0,30};\s*it is\b", "it is not X; it is Y"),
          (r"\bis not\b[^.;]{0,25},\s*(it|they)\s+(is|are)\b", "is not X, it is Y"),
          (r"\band not\b", "and not"),
          (r"\bnot\b[^.;]{0,30}\bbut\b", "not X but Y")]

print("=== overused AI vocabulary ===")
for w in VOCAB:
    n = len(re.findall(rf"\b{re.escape(w)}", flat, re.I))
    if n:
        print(f"  {w}: {n}")

print("\n=== copula avoidance ===")
for w in COPULA:
    n = len(re.findall(rf"\b{re.escape(w)}\b", flat, re.I))
    if n:
        print(f"  {w}: {n}")

print("\n=== significance / legacy inflation ===")
for w in SIGNIF + SUPERFICIAL + VAGUE:
    n = len(re.findall(rf"\b{re.escape(w)}", flat, re.I))
    if n:
        print(f"  {w}: {n}")

print("\n=== negative parallelism ===")
for pat, name in NEGPAR:
    hits = re.findall(pat, flat, re.I)
    if hits:
        print(f"  {name}: {len(hits)}")

print("\n=== em dashes / boldface / curly quotes ===")
print(f"  em dash (---): {body.count('---')}")
print(f"  \\textbf in body: {len(re.findall(r'\\textbf', body))}")
print(f"  curly quotes: {len(re.findall(chr(8216)+'|'+chr(8217)+'|'+chr(8220)+'|'+chr(8221), body))}")

print("\n=== elegant variation: near-synonyms for the same quantity ===")
for group in [["increment", "gain", "improvement", "uplift", "boost"],
              ["null", "negative result", "absence of effect"],
              ["cohort", "sample", "dataset", "set of patients"],
              ["mesh", "anatomy", "geometry"]]:
    counts = {w: len(re.findall(rf"\b{w}", flat, re.I)) for w in group}
    used = {w: c for w, c in counts.items() if c}
    if len(used) > 1:
        print(f"  {used}")
