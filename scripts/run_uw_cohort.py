"""Label the UW/Boyle cohort and run the expanded (cross-cohort) GM1."""
import time

from asb.experiments.uw_cohort import run_gm1_expanded

if __name__ == "__main__":
    t0 = time.time()
    m = run_gm1_expanded(n_jobs=2, n_boot=10000)
    print(f"=== UW EXPANDED GM1 DONE in {time.time() - t0:.0f}s ===", flush=True)
    for name, c in m["cohorts"].items():
        print(f"[{name:9s}] n={c['n_subjects']:>3} inducible={c.get('n_inducible', 0):>3} "
              f"({c.get('inducible_fraction', 0):.2f})", flush=True)
        for cmp in ("competitors_vs_+SFI", "fibrosis_vs_+SFI"):
            cv = c.get(cmp)
            if not cv:
                continue
            for clf in ("lr", "gbt"):
                d = cv[clf]
                print(f"    {cmp:22s} {clf}: base={d['grouped_auc_base']:.3f} "
                      f"+SFI={d['grouped_auc_sfi']:.3f} dAUC={d['grouped_delta_auc']:+.3f} "
                      f"p={d['delong_p']:.3f} {'MET' if d['endpoint_met'] else 'null'}", flush=True)
    print("=== UW_DONE ===", flush=True)
