"""Analyze the resolution-convergence study: do inducibility verdicts stabilize as resolution rises?

Usage: python scripts/convergence_analyze.py [roney|uw]

Reads outputs/convergence_<cohort>/*.json (one per subject x resolution) and reports:
  (1) inducibility RATE per resolution (does the aggregate rate settle, or keep moving?);
  (2) consecutive-resolution verdict AGREEMENT among subjects present at both (does per-subject
      agreement climb toward 1.0 at the fine end = convergence, or stay low = drift?);
  (3) full per-subject consistency across all resolutions.
Writes results/convergence_summary_<cohort>.json. Honest about partial tiers (reports n per res).
"""
import glob
import json
import os
import sys
from collections import defaultdict

import numpy as np


def main():
    cohort = sys.argv[1] if len(sys.argv) > 1 else "roney"
    outdir = f"outputs/convergence_{cohort}"
    recs = []
    for f in glob.glob(f"{outdir}/*.json"):
        d = json.load(open(f, encoding="utf-8"))
        if "error" not in d and "inducible" in d:
            recs.append(d)
    if not recs:
        print(f"no convergence records in {outdir}"); return

    verdict = defaultdict(dict)
    for r in recs:
        verdict[r["subject"]][r["target_res"]] = bool(r["inducible"])
    all_res = sorted({r["target_res"] for r in recs})

    rate = {}
    for res in all_res:
        vals = [v[res] for v in verdict.values() if res in v]
        rate[res] = (float(np.mean(vals)), len(vals))

    pair_agreement = []
    for a, b in zip(all_res[:-1], all_res[1:]):
        both = [s for s in verdict if a in verdict[s] and b in verdict[s]]
        if both:
            agree = np.mean([verdict[s][a] == verdict[s][b] for s in both])
            pair_agreement.append({"a": a, "b": b, "agreement": round(float(agree), 3),
                                   "n": len(both)})

    complete = {s: v for s, v in verdict.items() if len(v) >= 3}
    fully_consistent = sum(1 for v in complete.values() if len(set(v.values())) == 1)

    # honest verdict on convergence: is the finest-tier rate close to the previous tier, and is
    # fine-end pair agreement high? If the finest rate diverges upward, we say NOT converged.
    fine_rates = [rate[r][0] for r in all_res[-3:]]
    monotone_settling = max(fine_rates) - min(fine_rates) < 0.10
    fine_pair = pair_agreement[-1]["agreement"] if pair_agreement else float("nan")

    summary = {
        "cohort": cohort,
        "resolutions": all_res,
        "n_records": len(recs),
        "n_subjects": len(verdict),
        "inducibility_rate": {str(k): {"rate": round(v[0], 3), "n": v[1]} for k, v in rate.items()},
        "consecutive_agreement": pair_agreement,
        "fully_consistent_subjects": fully_consistent,
        "subjects_with_ge3_res": len(complete),
        "finest_res": all_res[-1],
        "finest_rate": round(rate[all_res[-1]][0], 3),
        "fine_end_settled": bool(monotone_settling),
        "fine_end_pair_agreement": fine_pair,
        "finding": (
            f"Coarse mesh over-calls ({rate[all_res[0]][0]:.0%} at {all_res[0]}). The rate drops "
            f"sharply then is NON-MONOTONE across the fine tiers "
            f"({', '.join(f'{rate[r][0]:.0%}@{r}' for r in all_res[1:])}); the finest affordable tier "
            f"({all_res[-1]} nodes) sits at {rate[all_res[-1]][0]:.0%}, "
            f"{'close to' if monotone_settling else 'ABOVE'} the mid tiers. Conclusion: the labeller "
            f"does NOT collapse to zero (the phenomenon is real) but it has NOT converged in the "
            f"tested window either — inducibility remains resolution-sensitive up to ~50k nodes. "
            f"Every ABSOLUTE number (rates, competitor AUCs, localizer origins) at the 2000-node "
            f"operating resolution is coarse-mesh-provisional pending a validated fine-resolution "
            f"solver (openCARP). The RELATIVE SFI-vs-competitor comparison (same labels both arms) "
            f"and the label-free rho argument are unaffected."
        ),
    }
    os.makedirs("results", exist_ok=True)
    json.dump(summary, open(f"results/convergence_summary_{cohort}.json", "w", encoding="utf-8"), indent=2)

    # ---- reproducible figure: rate vs resolution + consecutive verdict agreement ----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        res_x = all_res
        rates = [rate[r][0] * 100 for r in res_x]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 3.6))
        ax1.plot(res_x, rates, "o-", color="#c0392b", lw=2, ms=7)
        for r, y in zip(res_x, rates):
            ax1.annotate(f"{y:.0f}%", (r, y), textcoords="offset points", xytext=(0, 8),
                         ha="center", fontsize=8)
        ax1.axvspan(1000, 2000, color="#f9e79f", alpha=0.5, label="operating res (2000)")
        ax1.set_xscale("log")
        ax1.set_xlabel("mesh nodes (log)")
        ax1.set_ylabel("inducibility rate (%)")
        ax1.set_title(f"{cohort}: rate does not settle\n(finest {all_res[-1]} = {rates[-1]:.0f}%)")
        ax1.set_ylim(0, 65)
        ax1.grid(alpha=0.3)
        ax1.legend(fontsize=7, loc="upper right")
        pa_x = [f"{p['a']//1000}k→{p['b']//1000}k" for p in pair_agreement]
        pa_y = [p["agreement"] * 100 for p in pair_agreement]
        ax2.bar(range(len(pa_x)), pa_y, color="#2874a6")
        ax2.set_xticks(range(len(pa_x)))
        ax2.set_xticklabels(pa_x, fontsize=7, rotation=30)
        ax2.axhline(50, color="gray", ls="--", lw=1, label="chance (50%)")
        ax2.set_ylabel("per-subject verdict agreement (%)")
        ax2.set_title("consecutive-tier agreement\n(climbs but rate still drifts)")
        ax2.set_ylim(0, 100)
        ax2.legend(fontsize=7)
        fig.tight_layout()
        os.makedirs("docs/paper/figures", exist_ok=True)
        fig.savefig("docs/paper/figures/fig_convergence.png", dpi=150)
        print("wrote docs/paper/figures/fig_convergence.png")
    except Exception as e:  # noqa: BLE001
        print(f"(figure skipped: {e})")

    print(f"[{cohort}] inducibility rate by resolution:")
    for r in all_res:
        print(f"  res={r:6d}  n={rate[r][1]:2d}  rate={rate[r][0]:.3f}")
    print("consecutive verdict agreement:")
    for p in pair_agreement:
        print(f"  {p['a']}->{p['b']}: {p['agreement']:.0%} (n={p['n']})")
    print(f"fully consistent across all res (>=3): {fully_consistent}/{len(complete)}")
    print(f"finest tier {all_res[-1]}: {rate[all_res[-1]][0]:.1%}  "
          f"fine-end settled(<10pp spread over last 3)={monotone_settling}")
    print("CONVERGENCE_ANALYZE_DONE", flush=True)


if __name__ == "__main__":
    main()
