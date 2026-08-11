"""Scan the condensed article for the cadences that mark machine-drafted prose:
recycled contrast frames, adverb openers, tricolons, cleft closers, and em-dash
density. Reports counts and locations; changes nothing."""
import re
import sys
from collections import Counter

src = open(sys.argv[1], encoding="utf-8").read()
body = src.split("%% ===================== KEY POINTS =====================")[1]
body = body.split("%% ===================== REFERENCES =====================")[0]
body = re.sub(r"(?<!\\)%.*", "", body)
prose = re.sub(r"\\[a-zA-Z]+\*?(\[[^\]]*\])?", " ", body)
prose = re.sub(r"[{}$]", " ", prose)
flat = " ".join(prose.split())
sentences = re.split(r"(?<=[.!?])\s+(?=[A-Z(])", flat)

print(f"sentences: {len(sentences)}   words: {len(flat.split())}")

FRAMES = [
    (r"\brather than\b", "rather than"),
    (r"\bnot\b[^.;]{0,40}\bbut\b", "not X but Y"),
    (r"\binstead of\b", "instead of"),
    (r"\bIt is (worth|important|notable)\b", "it is worth/important"),
    (r"\b(Importantly|Notably|Crucially|Decisively|Critically|Ultimately|Moreover|Furthermore)\b",
     "stock adverb"),
    (r"---", "em dash"),
    (r"\bWhat is [a-z]+ is\b|\bWhat [a-z]+ is\b", "cleft closer"),
    (r"\bdoes not by itself\b", "does-not-by-itself"),
    (r"\bthe (contribution|result|point|thesis) here is\b", "meta framing"),
    (r"\bboth\b[^.;]{0,60}\band\b", "both...and"),
    (r"\bis not a\b", "is not a"),
    (r"\bnot (only|merely|simply)\b", "not only/merely"),
]
print("\n=== recycled frames ===")
for pat, name in FRAMES:
    hits = re.findall(pat, flat)
    if hits:
        print(f"  {name:22s} {len(hits)}")

print("\n=== tricolons (A, B, and C) ===")
tri = re.findall(r"\b[\w-]+(?:\s+[\w-]+){0,3},\s+[\w-]+(?:\s+[\w-]+){0,3},\s+and\s+[\w-]+", flat)
print(f"  {len(tri)}")
for t in tri[:12]:
    print(f"    {t}")

print("\n=== sentence openers used more than once ===")
op = Counter(" ".join(s.split()[:2]) for s in sentences)
for k, v in op.most_common():
    if v > 1:
        print(f"  {v}x  {k}")

print("\n=== long sentences (>45 words) ===")
for s in sentences:
    n = len(s.split())
    if n > 45:
        print(f"  [{n}] {s[:110]}...")

print("\n=== semicolon and colon density ===")
print(f"  semicolons {flat.count(';')}   colons {flat.count(':')}")
