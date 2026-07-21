# GM3 — validity radius + subspace/exact SFI recovery (real Roney, monodomain-MS labels)

> Labels are a monodomain Mitchell-Schaeffer **simulator verdict**, not clinical POAF, not openCARP.

- Subjects **62**, inducible **20** (0.3226), patient-held-out GroupKFold.

## Part A — does a correctly-computed SFI recover the signal?

| add-on to competitors | clf | grouped base | grouped +SFI | grouped ΔAUC | DeLong p | boot ΔAUC [95% CI] | endpoint |
| --- | --- | --- | --- | --- | --- | --- | --- |
| single_vector | lr | 0.8321 | 0.8036 | -0.0286 | 0.6558 | -0.0295 [-0.1643, 0.0881] | not met |
| single_vector | gbt | 0.8202 | 0.8202 | -0.0000 | 1.0000 | -0.0000 [-0.0298, 0.0286] | not met |
| subspace | lr | 0.8321 | 0.8107 | -0.0214 | 0.4467 | -0.0212 [-0.0821, 0.0310] | not met |
| subspace | gbt | 0.8202 | 0.8310 | 0.0107 | 0.8122 | 0.0108 [-0.0857, 0.0905] | not met |
| exact | lr | 0.8321 | 0.8048 | -0.0274 | 0.3302 | -0.0270 [-0.0857, 0.0262] | not met |
| exact | gbt | 0.8202 | 0.8095 | -0.0107 | 0.6402 | -0.0103 [-0.0595, 0.0345] | not met |

## Part B — validity radius (synthetic, exact vs first-order vs subspace)

- λ2=0.0384, λ3=0.1522, gap=0.1138; **ρ\* (first-order 10% error) = 2.9998846240069104**

| ρ=‖ΔL‖/gap | exact Δλ2 | first-order err | subspace err | Weyl ok |
| --- | --- | --- | --- | --- |
| 0.3000 | 0.0000 | 0.0106 | 11.9115 | True |
| 0.6000 | 0.0000 | 0.0211 | 11.7737 | True |
| 1.4999 | 0.0001 | 0.0528 | 11.3602 | True |
| 2.9999 | 0.0003 | 0.1056 | 10.6717 | True |
| 5.9998 | 0.0006 | 0.2109 | 9.2971 | True |
| 10.4996 | 0.0013 | 0.3682 | 7.2447 | True |
| 16.4994 | 0.0031 | 0.5751 | 4.5443 | True |
| 22.4991 | 0.0078 | 0.7719 | 1.9762 | True |
| 28.4989 | 0.0277 | 0.9189 | 0.0586 | True |

> Real-cohort ρ (from GM1 validity scan) has median ≈ 2422 — far right of this table, i.e. deep in the regime where only the exact recompute is trustworthy.

## Part B(ii) — near-degeneracy: subspace SFI vs single Fiedler vector

As the bridge shrinks the λ2–λ3 gap closes; error is vs the exact drop in the degeneracy-invariant sum λ2+λ3.

| bridge | gap(λ3−λ2) | single-vector err | subspace err |
| --- | --- | --- | --- |
| 1.0000 | 0.1673 | 0.0846 | 0.0681 |
| 0.3000 | 0.0660 | 0.0740 | 0.0713 |
| 0.1000 | 0.0239 | 0.1421 | 0.0028 |
| 0.0300 | 0.0074 | 0.0760 | 0.0754 |
| 0.0100 | 0.0025 | 0.1417 | 0.0003 |
| 0.0030 | 0.0007 | 0.0069 | 0.1004 |
