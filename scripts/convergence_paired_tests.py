"""Paired tests between adjacent mesh-resolution tiers.

The manuscript compares the six convergence tiers using marginal Clopper-Pearson
intervals and says so explicitly: those intervals discard the pairing, and
overlap between two of them is not a test of their difference.  It then declines
to say whether the large 1500 -> 3000 drop is real, because the paired
comparison "that would settle it" was not reported.

It is recoverable.  The same 24 subjects are simulated at every tier, so for any
adjacent pair the 2x2 table is uniquely determined by three stored quantities:

    a + b = n_A          (subjects inducible at the coarser tier)
    a + c = n_B          (subjects inducible at the finer tier)
    a + d = n_agree      (subjects with the same verdict at both)

with a + b + c + d = 24, giving

    b = ((n_A - n_B) + (24 - n_agree)) / 2
    c = ((n_B - n_A) + (24 - n_agree)) / 2

No per-subject data is needed.  Every reconstructed cell must come out a
non-negative integer; the script asserts that, which is a real check on the
stored summaries rather than a formality.

Cochran's Q across all six tiers is NOT computed: it needs the full 24x6 verdict
matrix, which is not stored.  Only the adjacent-pair tests are derivable, and
the script reports exactly those.

Writes results/convergence_paired_tests.json.
"""
from __future__ import annotations

import json
import math
import os
from typing import Dict, List

SRC = "results/convergence_summary_roney.json"
OUT = "results/convergence_paired_tests.json"


def mcnemar_exact(b: int, c: int) -> float:
    """Two-sided exact McNemar: binomial(b+c, 0.5) tail on the discordant pairs."""
    n = b + c
    if n == 0:
        return 1.0
    k = min(b, c)
    tail = sum(math.comb(n, i) for i in range(0, k + 1)) / (2 ** n)
    return float(min(1.0, 2.0 * tail))


def holm(pvals: List[float]) -> List[float]:
    """Holm-Bonferroni adjusted p-values, preserving input order."""
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    adj = [0.0] * m
    running = 0.0
    for rank, idx in enumerate(order):
        val = (m - rank) * pvals[idx]
        running = max(running, val)
        adj[idx] = min(1.0, running)
    return adj


def main() -> None:
    with open(SRC, encoding="utf-8") as fh:
        d = json.load(fh)

    n = int(d["inducibility_rate"][str(d["resolutions"][0])]["n"])
    counts: Dict[int, int] = {}
    for res in d["resolutions"]:
        rate = d["inducibility_rate"][str(res)]["rate"]
        k = round(rate * n)
        # the stored rate is rounded to 3 dp; k must reproduce it
        assert abs(k / n - rate) < 0.5 / n, f"tier {res}: rate {rate} inconsistent with n={n}"
        counts[int(res)] = int(k)

    rows = []
    for entry in d["consecutive_agreement"]:
        a_res, b_res = int(entry["a"]), int(entry["b"])
        assert int(entry["n"]) == n, "tiers do not share the same subject count"
        n_agree = round(entry["agreement"] * n)
        n_a, n_b = counts[a_res], counts[b_res]

        discordant = n - n_agree
        two_b = (n_a - n_b) + discordant
        two_c = (n_b - n_a) + discordant
        assert two_b % 2 == 0 and two_c % 2 == 0, f"{a_res}->{b_res}: 2x2 not integral"
        b, c = two_b // 2, two_c // 2
        a = n_a - b
        dd = n_agree - a
        assert min(a, b, c, dd) >= 0, f"{a_res}->{b_res}: negative cell in reconstructed 2x2"
        assert a + b + c + dd == n, f"{a_res}->{b_res}: cells do not sum to n"

        rows.append({
            "coarse": a_res,
            "fine": b_res,
            "n": n,
            "n_inducible_coarse": n_a,
            "n_inducible_fine": n_b,
            "n_agree": n_agree,
            "lost_on_refinement": b,      # inducible coarse, not fine
            "gained_on_refinement": c,    # not inducible coarse, inducible fine
            "mcnemar_exact_p": mcnemar_exact(b, c),
        })

    adj = holm([r["mcnemar_exact_p"] for r in rows])
    for r, p in zip(rows, adj):
        r["holm_adjusted_p"] = p
        r["resolvable_at_05"] = bool(p < 0.05)

    out = {
        "description": (
            "Exact McNemar tests between adjacent mesh-resolution tiers on the same 24 "
            "subjects. The 2x2 for each pair is reconstructed from the stored marginals "
            "and agreement count, which determine it uniquely; every cell is verified to "
            "be a non-negative integer. Cochran's Q across all six tiers is not computed "
            "because the full 24x6 verdict matrix is not stored."
        ),
        "source": SRC,
        "n_subjects": n,
        "n_comparisons": len(rows),
        "multiplicity": "Holm-Bonferroni across the adjacent-pair family",
        "pairs": rows,
        "n_resolvable": sum(1 for r in rows if r["resolvable_at_05"]),
    }

    os.makedirs("results", exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)

    print("paired tests between adjacent resolution tiers (n=24, same subjects)\n")
    print(f"{'pair':>16} {'coarse':>7} {'fine':>5} {'lost':>5} {'gained':>7} "
          f"{'McNemar p':>10} {'Holm p':>8}  verdict")
    for r in rows:
        print(f"{str(r['coarse'])+'->'+str(r['fine']):>16} "
              f"{r['n_inducible_coarse']:>7} {r['n_inducible_fine']:>5} "
              f"{r['lost_on_refinement']:>5} {r['gained_on_refinement']:>7} "
              f"{r['mcnemar_exact_p']:>10.4f} {r['holm_adjusted_p']:>8.4f}  "
              f"{'RESOLVABLE' if r['resolvable_at_05'] else 'not resolvable'}")
    print(f"\n{out['n_resolvable']}/{len(rows)} adjacent transitions are resolvable at alpha=0.05")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
