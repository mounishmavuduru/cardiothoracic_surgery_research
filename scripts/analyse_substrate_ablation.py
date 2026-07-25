"""Paired analysis of the UW 2x2 substrate ablation.

Arms share the same 82 subjects, so the correct test is McNemar's exact test on the
discordant pairs, not a two-sample proportion test. Also reports the factorial
decomposition and, importantly, the per-subject CHURN between arms -- if the labeller's
verdicts are unstable, a change in the aggregate rate can hide wholesale relabelling.

Usage:  python scripts/analyse_substrate_ablation.py [--in results/uw_substrate_ablation.json]
"""
from __future__ import annotations

import argparse
import json
from itertools import combinations

from scipy.stats import binomtest

RONEY_RATE = 20 / 62  # fibre-carrying cohort, identical frozen protocol


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="src", default="results/uw_substrate_ablation.json")
    ap.add_argument("--write", action="store_true", help="write the stats back into the JSON")
    a = ap.parse_args()

    d = json.load(open(a.src, encoding="utf-8"))
    rows = d["rows"]
    arms = list(dict.fromkeys(r["arm"] for r in rows))
    by = {arm: {r["subject"]: bool(r["inducible"]) for r in rows if r["arm"] == arm}
          for arm in arms}
    subs = sorted(by[arms[0]])
    assert all(sorted(by[arm]) == subs for arm in arms), "arms cover different subjects"
    n = len(subs)

    print(f"n = {n} subjects, {len(arms)} arms\n")
    counts = {arm: sum(by[arm].values()) for arm in arms}
    for arm in arms:
        print(f"  {arm:<34} {counts[arm]:>3}/{n}  {counts[arm]/n:6.1%}")

    print("\nMcNemar exact test on discordant pairs")
    tests = {}
    for x, y in combinations(arms, 2):
        b = sum(by[x][s] and not by[y][s] for s in subs)
        c = sum(by[y][s] and not by[x][s] for s in subs)
        p = binomtest(min(b, c), b + c, 0.5).pvalue if (b + c) else 1.0
        tests[f"{x} vs {y}"] = {"only_first": b, "only_second": c, "p": p}
        print(f"  {x[0]} vs {y[0]}:  discordant {b}/{c}   p = {p:.4f}"
              f"{'  *' if p < 0.05 else ''}")

    A, B, C, D = (counts[arm] for arm in arms)
    print("\nfactorial decomposition (subjects turning inducible)")
    print(f"  opening the atrial orifices : {((B - A) + (D - C)) / 2:+.1f}")
    print(f"  varying the fibre field     : {((C - A) + (D - B)) / 2:+.1f}")
    print(f"  interaction                 : {(D - C) - (B - A):+d}")

    print("\nper-subject churn (verdict stability)")
    for arm in arms[1:]:
        kept = sum(by[arms[0]][s] and by[arm][s] for s in subs)
        gained = sum(by[arm][s] and not by[arms[0]][s] for s in subs)
        print(f"  {arms[0][0]} -> {arm[0]}: {kept}/{A} baseline positives retained, "
              f"{gained} new")

    closed = (D / n - A / n) / (RONEY_RATE - A / n)
    print(f"\ngap to the fibre-carrying cohort (Roney {RONEY_RATE:.1%})")
    print(f"  {A}/{n} = {A/n:.1%}  ->  {D}/{n} = {D/n:.1%}   closes {closed:.0%} of the gap")

    if a.write:
        d["paired_stats"] = {
            "counts": counts, "mcnemar": tests,
            "effect_open_orifices": ((B - A) + (D - C)) / 2,
            "effect_vary_fibres": ((C - A) + (D - B)) / 2,
            "interaction": (D - C) - (B - A),
            "fraction_of_roney_gap_closed": closed,
            "baseline_positives_retained_in_full_fix": sum(
                by[arms[0]][s] and by[arms[-1]][s] for s in subs),
        }
        json.dump(d, open(a.src, "w", encoding="utf-8"), indent=2)
        print(f"\nwrote paired_stats into {a.src}")


if __name__ == "__main__":
    main()
