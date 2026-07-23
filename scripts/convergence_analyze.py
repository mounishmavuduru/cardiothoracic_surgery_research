"""Analyze the resolution-convergence study: do inducibility verdicts stabilize as resolution rises?

Reads outputs/convergence/*.json (one per subject×resolution) and reports, per consecutive
resolution pair, the fraction of subjects whose inducible/not verdict AGREES — the key signal is
whether that agreement climbs toward 1.0 at the fine end (convergence) or stays low (drift). Runs on
partial data too (only uses resolutions with enough subjects). Writes results/convergence_summary.json.
"""
import glob
import json
import os
from collections import defaultdict

import numpy as np


def main():
    recs = []
    for f in glob.glob("outputs/convergence/*.json"):
        d = json.load(open(f))
        if "error" not in d and "inducible" in d:
            recs.append(d)
    if not recs:
        print("no convergence records yet"); return

    # verdict[subject][res] = bool
    verdict = defaultdict(dict)
    sustained = defaultdict(dict)
    for r in recs:
        verdict[r["subject"]][r["target_res"]] = bool(r["inducible"])
        sustained[r["subject"]][r["target_res"]] = float(r.get("sustained_ms", 0.0))
    all_res = sorted({r["target_res"] for r in recs})

    # inducibility rate per resolution
    rate = {}
    for res in all_res:
        vals = [v[res] for v in verdict.values() if res in v]
        rate[res] = (float(np.mean(vals)), len(vals))

    # consecutive-resolution verdict agreement (subjects present at both)
    pair_agreement = []
    for a, b in zip(all_res[:-1], all_res[1:]):
        both = [s for s in verdict if a in verdict[s] and b in verdict[s]]
        if both:
            agree = np.mean([verdict[s][a] == verdict[s][b] for s in both])
            pair_agreement.append({"res_a": a, "res_b": b, "n_subjects": len(both),
                                   "verdict_agreement": float(agree)})

    # full-consistency across ALL resolutions a subject has (need >=3 res)
    complete = {s: v for s, v in verdict.items() if len(v) >= 3}
    fully_consistent = sum(1 for v in complete.values() if len(set(v.values())) == 1)

    summary = {
        "n_records": len(recs), "resolutions": all_res, "n_subjects": len(verdict),
        "inducibility_rate_by_res": {str(k): {"rate": round(v[0], 3), "n": v[1]} for k, v in rate.items()},
        "consecutive_pair_agreement": pair_agreement,
        "fully_consistent_subjects": fully_consistent,
        "subjects_with_ge3_res": len(complete),
        "interpretation": ("verdict_agreement climbing toward 1.0 at the fine end => the labeller "
                           "converges (coarse-mesh artifacts, not phenomenon instability); staying low "
                           "=> the inducibility phenomenon itself is resolution-sensitive and needs a "
                           "validated fine-resolution solver (openCARP)."),
    }
    os.makedirs("results", exist_ok=True)
    json.dump(summary, open("results/convergence_summary.json", "w"), indent=2)
    print("inducibility rate by resolution:",
          {k: (round(v[0], 2), v[1]) for k, v in rate.items()}, flush=True)
    print("consecutive verdict agreement:")
    for p in pair_agreement:
        print(f"  {p['res_a']}->{p['res_b']}: {p['verdict_agreement']:.0%} (n={p['n_subjects']})", flush=True)
    print(f"fully consistent across all res (>=3): {fully_consistent}/{len(complete)}", flush=True)
    print("CONVERGENCE_ANALYZE_DONE", flush=True)


if __name__ == "__main__":
    main()
