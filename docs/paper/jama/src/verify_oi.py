"""Every numeral in the condensed article must already appear in the audited
master manuscript. The master is itself traced to results/*.json, so matching
against it is the whole chain. Anything unmatched is printed with context and
must be resolved by hand -- there is no whitelist that hides a miss."""
import re
import sys

OI, MASTER = sys.argv[1], sys.argv[2]
oi = open(OI, encoding="utf-8").read()
master = open(MASTER, encoding="utf-8").read()

# main article only: drop preamble/title page (word counts, phone, ORCID, address)
# and the bibliography (page/volume numbers, DOIs, PMIDs) -- both are checked
# separately by eye.
body = oi.split("%% ===================== KEY POINTS =====================")[1]
body = body.split("%% ===================== REFERENCES =====================")[0] + \
    oi.split("%% ===================== FLOATS =====================")[-1]
body = re.sub(r"(?<!\\)%.*", "", body)


def normalise(text):
    # LaTeX thousands separators vanish; every other comma becomes a boundary so
    # that "[1082,7384]" cannot fuse into one bogus numeral.
    return text.replace("{,}", "").replace("\\,", "").replace(",", " ")


def numerals(text):
    out = []
    t = normalise(text)
    # blank out control sequences so "\to0.345" still yields 0.345
    t = re.sub(r"\\[a-zA-Z]+", lambda m: " " * len(m.group(0)), t)
    for m in re.finditer(r"(?<![\w.])(\d+(?:\.\d+)?)", t):
        out.append((m.group(1), m.start()))
    return out


mset = {n for n, _ in numerals(master)}

# The condensed article may quote a stored result the master states only in
# prose (e.g. the gradient-boosting clinical interval). Those are traced to the
# result files directly, at every rounding the article could legitimately use.
import glob
import json


def walk(o, sink):
    if isinstance(o, dict):
        for v in o.values():
            walk(v, sink)
    elif isinstance(o, (list, tuple)):
        for v in o:
            walk(v, sink)
    elif isinstance(o, bool):
        pass
    elif isinstance(o, (int, float)):
        sink.append(float(o))


vals = []
for p in glob.glob("results/*.json"):
    try:
        walk(json.load(open(p, encoding="utf-8")), vals)
    except Exception:
        pass
for v in vals:
    if v != v or v in (float("inf"), float("-inf")):
        continue
    for dp in range(0, 5):
        mset.add(f"{round(v, dp):.{dp}f}".rstrip(".") if dp else f"{round(v):d}")
        mset.add(f"{abs(round(v, dp)):.{dp}f}")
        mset.add(str(abs(round(v, dp))))
# the master writes some values only inside \code{} paths or with different
# thousands separators; normalise both sides the same way, above.

missing = []
seen = set()
for n, pos in numerals(body):
    if n in mset or n in seen:
        continue
    seen.add(n)
    nb = normalise(body)
    ctx = " ".join(nb[max(0, pos - 90):pos + 60].split())
    missing.append((n, ctx))

print(f"distinct numerals in article body: {len({n for n, _ in numerals(body)})}")
print(f"unmatched against master: {len(missing)}")
for n, ctx in missing:
    print(f"\n  [{n}] ...{ctx}...")
sys.exit(1 if missing else 0)
