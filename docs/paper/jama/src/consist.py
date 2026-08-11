"""House-style consistency for the article body: numerals set in math mode,
thousands separators from five digits up, no stray bare numbers in prose."""
import re
import sys

t = open(sys.argv[1], encoding="utf-8").read()
body = t.split("%% ===================== KEY POINTS =====================")[1]
body = body.split("%% ===================== REFERENCES =====================")[0]
body = re.sub(r"(?<!\\)%.*", "", body)

# blank every math span so what remains is prose
prose = re.sub(r"\$[^$]*\$", lambda m: " " * len(m.group(0)), body)
prose = re.sub(r"\\begin\{equation\}.*?\\end\{equation\}", " ", prose, flags=re.S)

ALLOWED = re.compile(
    r"(Zenodo|dryad|10\.5061|Table~?|Figure~?|Supplement|TRIPOD\+AI|eMethods"
    r"|includegraphics|linewidth|\\S)")

print("=== bare numerals in prose (should be captions/labels only) ===")
n = 0
for m in re.finditer(r"(?<![\w.$])\d+(?:\.\d+)?", prose):
    ctx = " ".join(prose[max(0, m.start() - 60):m.end() + 30].split())
    if ALLOWED.search(ctx):
        continue
    print(f"  [{m.group(0)}] ...{ctx}...")
    n += 1
print(f"  ({n} found)")

print("\n=== thousands separators ===")
sep = sorted(set(re.findall(r"\d+\{,\}\d{3}", body)))
unsep = sorted({x for x in re.findall(r"(?<![\d.{])\d{4,}(?![\d}])", body)})
print(f"  separated:   {sep}")
print(f"  unseparated: {unsep}")
print("  rule: separate at five digits and above")
bad = [x for x in unsep if len(x) >= 5 and not x.startswith("5801337")]
bad += [x for x in sep if len(x.replace("{,}", "")) <= 4]
print(f"  violations:  {bad if bad else 'none'}")

print("\n=== spelling variants ===")
for a, b in [("in-silico", "in silico"), ("pre-registered", "preregistered"),
             ("shape-family", "shape family"), ("out-of-fold", "out of fold"),
             ("labeller", "labeler"), ("artefact", "artifact")]:
    ca, cb = len(re.findall(a, body)), len(re.findall(rf"\b{b}\b", body))
    if ca and cb:
        print(f"  MIXED: {a}={ca}  {b}={cb}")
print("  (only mixed pairs shown)")
