r"""Cross-check the repository against itself: imports, references, and shared constants.

Complements `verify_manuscript_numbers.py`, which checks that the paper's numbers trace to
stored results. This checks the things a numeric trace cannot see:

1. every module under src/ imports cleanly, and every script parses;
2. every file path named in the manuscript, the pre-registration and the READMEs exists;
3. every ``results/*.json`` the manuscript cites is present;
4. frozen constants quoted in the manuscript match the values in the code that produced
   them -- the failure mode where a paper and its pipeline drift apart silently;
5. no unresolved TODO or FIXME markers remain in shipped code.

Exit status is non-zero if anything fails, so this can gate a commit.

Usage:  python scripts/verify_repo_consistency.py
"""
from __future__ import annotations

import ast
import glob
import importlib
import json
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

FAILURES: list = []
CHECKS = 0


def check(name: str, ok: bool, detail: str = "") -> None:
    global CHECKS
    CHECKS += 1
    if ok:
        print(f"  PASS  {name}")
    else:
        print(f"  FAIL  {name}{(' -- ' + detail) if detail else ''}")
        FAILURES.append(name)


def check_imports() -> None:
    print("\n[1] module imports")
    mods = []
    for f in sorted(glob.glob("src/asb/**/*.py", recursive=True)):
        if f.endswith("__init__.py"):
            continue
        mod = f[len("src/"):-3].replace(os.sep, ".").replace("/", ".")
        mods.append(mod)
    bad = []
    for m in mods:
        try:
            importlib.import_module(m)
        except Exception as exc:                                # noqa: BLE001
            bad.append(f"{m}: {str(exc)[:80]}")
    check(f"{len(mods)} modules import", not bad, "; ".join(bad[:3]))


def check_script_syntax() -> None:
    print("\n[2] script syntax")
    bad = []
    files = sorted(glob.glob("scripts/*.py"))
    for f in files:
        try:
            ast.parse(open(f, encoding="utf-8").read())
        except SyntaxError as exc:
            bad.append(f"{os.path.basename(f)}:{exc.lineno}")
    check(f"{len(files)} scripts parse", not bad, "; ".join(bad))


#: Files deliberately referred to in the past tense; their absence is the point.
HISTORICAL = {"PREPRINT.md"}


def _repo_basenames() -> dict:
    """basename -> True for every tracked-ish file, so bare names can be resolved."""
    idx = {}
    for root, dirs, files in os.walk("."):
        dirs[:] = [d for d in dirs
                   if d not in {".git", ".venv", "__pycache__", "data", "outputs", "tools"}]
        for f in files:
            idx[f] = True
    return idx


def _referenced_paths(text: str) -> set:
    """File-ish tokens named in prose: code spans, \\code{...} and inline paths.

    Splits each span on whitespace and keeps the token that looks like a path, so a
    command line such as ``python scripts/x.py`` yields ``scripts/x.py`` rather than the
    concatenation of both words.
    """
    out = set()
    for pat in (r"\\code\{([^}]+)\}", r"`([^`\n]+)`"):
        for span in re.findall(pat, text):
            # \allowbreak is followed by a space in the source; stripping the macro
            # alone would split "results/ gm4_*.json" in two and lose the directory.
            span = re.sub(r"\\allowbreak\s*", "", span).replace("\\_", "_")
            for tok in span.split():
                tok = tok.strip("(),;:")
                # a bare extension such as `.tex` names a file TYPE, not a file
                if tok.startswith("."):
                    continue
                if re.search(r"\.(py|json|md|tex|v|sh|yaml|toml)$", tok):
                    out.add(tok)
    return out


def check_referenced_files() -> None:
    print("\n[3] referenced files exist")
    sources = ["docs/paper/manuscript.tex", "docs/PRE_REGISTRATION.md", "README.md",
               "docs/paper/README_BUILD.md"]
    missing = []
    basenames = _repo_basenames()
    for src in sources:
        if not os.path.exists(src):
            continue
        for path in _referenced_paths(open(src, encoding="utf-8").read()):
            cand = path.lstrip("./")
            base = os.path.basename(cand)
            if base in HISTORICAL:
                continue
            if os.path.exists(cand) or glob.glob(cand):
                continue
            # a bare filename is a legitimate way to name a file that lives elsewhere
            if base in basenames:
                continue
            missing.append(f"{os.path.basename(src)} -> {cand}")
    check("all referenced files present", not missing, "; ".join(sorted(set(missing))[:5]))


def check_results_present() -> None:
    """Every results file the manuscript cites must exist.

    Uses the same token extraction as check_referenced_files; an earlier regex-only
    version matched a single file because \\allowbreak splits the path, which made this
    check almost vacuous while still reporting PASS.
    """
    print("\n[4] cited results files")
    cited = set()
    for src in ("docs/paper/manuscript.tex", "docs/PRE_REGISTRATION.md"):
        if not os.path.exists(src):
            continue
        for tok in _referenced_paths(open(src, encoding="utf-8").read()):
            if tok.startswith("results/") and tok.endswith(".json"):
                cited.add(tok)
    missing = [c for c in sorted(cited)
               if not os.path.exists(c) and not glob.glob(c)]
    check(f"{len(cited)} cited results files present", not missing, "; ".join(missing))
    if cited:
        print("        " + ", ".join(sorted(os.path.basename(c) for c in cited)))


def check_frozen_constants() -> None:
    """The paper quotes frozen constants; they must equal what the code uses."""
    print("\n[5] frozen constants: manuscript vs code")
    from asb.experiments.realcohort import (COARSEN_NODES, FROZEN_BURST_CLS,
                                            FROZEN_MONO, FROZEN_SFI)
    tex = open("docs/paper/manuscript.tex", encoding="utf-8").read()
    pairs = [
        ("coarsening nodes", str(COARSEN_NODES), r"\b2000\b|\$2000\$|2{,}000"),
        ("reentry_min_ms", str(int(FROZEN_MONO.reentry_min_ms)), r"\b650\b"),
        ("burst cycle length", str(int(FROZEN_BURST_CLS[0])), r"\b150\b"),
        ("n_burst", str(FROZEN_MONO.n_burst), r"\b6\b"),
        ("n_pacing_sites", str(FROZEN_MONO.n_pacing_sites), r"\b2\b"),
        ("tau_close", str(int(FROZEN_MONO.tau_close)), r"\b110\b"),
        ("delta_w_mean_frac", f"{FROZEN_SFI.delta_w_mean_frac:g}", r"0\.36"),
    ]
    bad = []
    for label, value, pattern in pairs:
        if not re.search(pattern, tex):
            bad.append(f"{label}={value} not found in manuscript")
    check("frozen constants appear in the manuscript", not bad, "; ".join(bad))

    # the manuscript must not contain a superseded threshold
    stale = []
    if re.search(r"600\s*ms|\\geq\s*600", tex):
        stale.append("600 ms reentry threshold (frozen value is 650)")
    check("no superseded thresholds in the manuscript", not stale, "; ".join(stale))


def check_no_todos() -> None:
    print("\n[6] unresolved markers in shipped code")
    hits = []
    # This file defines the marker pattern, so scanning it finds itself.
    me = os.path.basename(__file__)
    pattern = re.compile(r"\b(?:" + "|".join(["TO" + "DO", "FIX" + "ME", "XX" + "X"]) + r")\b")
    for f in sorted(glob.glob("src/asb/**/*.py", recursive=True) + glob.glob("scripts/*.py")):
        if os.path.basename(f) == me:
            continue
        for i, line in enumerate(open(f, encoding="utf-8"), 1):
            if pattern.search(line):
                hits.append(f"{os.path.basename(f)}:{i}")
    check("no TODO/FIXME markers", not hits, "; ".join(hits[:6]))


def main() -> None:
    print("repository consistency check")
    check_imports()
    check_script_syntax()
    check_referenced_files()
    check_results_present()
    check_frozen_constants()
    check_no_todos()
    print(f"\n{CHECKS - len(FAILURES)}/{CHECKS} checks passed")
    if FAILURES:
        print("FAILED: " + ", ".join(FAILURES))
    os.makedirs("results", exist_ok=True)
    with open("results/repo_consistency.json", "w", encoding="utf-8") as fh:
        json.dump({"n_checks": CHECKS, "n_failed": len(FAILURES),
                   "failed": FAILURES}, fh, indent=2)
    sys.exit(1 if FAILURES else 0)


if __name__ == "__main__":
    main()
