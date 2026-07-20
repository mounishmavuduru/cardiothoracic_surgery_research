# GM1 — SFI vs competitors on real Roney anatomy (monodomain-MS labels)

> Labels are a monodomain Mitchell-Schaeffer **simulator verdict**, not clinical POAF and not openCARP (deferred).

- Subjects: **62** (patient-held-out GroupKFold), inducible **20** (0.3226)
- SFI features added: **9**; competitors: 6

## competitors_vs_+SFI

| clf | grouped base | grouped +SFI | grouped ΔAUC | DeLong p | naive base | naive +SFI | boot ΔAUC [95% CI] | endpoint |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| lr | 0.8321 | 0.8024 | -0.0298 | 0.6434 | 0.8488 | 0.7738 | -0.0313 [-0.1667, 0.0857] | not met |
| gbt | 0.8202 | 0.8560 | 0.0357 | 0.2570 | 0.8036 | 0.8357 | 0.0352 [-0.0286, 0.0964] | not met |

## fibrosis_vs_+SFI

| clf | grouped base | grouped +SFI | grouped ΔAUC | DeLong p | naive base | naive +SFI | boot ΔAUC [95% CI] | endpoint |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| lr | 0.7810 | 0.7976 | 0.0167 | 0.5944 | 0.8179 | 0.7774 | 0.0169 [-0.0452, 0.0821] | not met |
| gbt | 0.7917 | 0.8607 | 0.0690 | 0.0635 | 0.7756 | 0.8369 | 0.0688 [-0.0036, 0.1440] | not met |
