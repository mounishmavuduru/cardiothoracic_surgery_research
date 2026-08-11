"""Compare the numeral multiset of two builds of the article. A prose rewrite
must not add, drop, or alter a single number; anything reported here is either a
mistake or a change I have to justify in the notebook."""
import re
import sys
from collections import Counter


def nums(path):
    t = open(path, encoding="utf-8").read()
    t = t.split("%% ===================== KEY POINTS =====================")[1]
    t = re.sub(r"(?<!\\)%.*", "", t)
    t = t.replace("{,}", "").replace("\\,", "").replace(",", " ")
    t = re.sub(r"\\[a-zA-Z]+", lambda m: " " * len(m.group(0)), t)
    return Counter(re.findall(r"(?<![\w.])(\d+(?:\.\d+)?)", t))


a, b = nums(sys.argv[1]), nums(sys.argv[2])
added, dropped = b - a, a - b
print(f"baseline distinct {len(a)}  new distinct {len(b)}")
print(f"added:   {dict(added) if added else 'none'}")
print(f"dropped: {dict(dropped) if dropped else 'none'}")
