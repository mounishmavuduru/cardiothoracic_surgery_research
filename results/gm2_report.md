# GM2 — hotspot ↔ reentry-origin localization (real Roney, monodomain-MS labels)

> Localizer argmax → geodesic distance (mm) to the monodomain reentry origin, vs a UAC rotational/shift spatial null. Simulator verdict, not clinical POAF.

- Inducible subjects localized: **20**

### Primary endpoint — argmax geodesic error vs UAC spatial null

| localizer | median error (mm) | norm (×diam) | null median | null 5th-pct | p(median) | **endpoint** |
| --- | --- | --- | --- | --- | --- | --- |
| grad_phi2 | 48.322 | 0.472 [0.350,0.656] | 0.523 | 0.406 | 0.229 | not met |
| perron | 59.155 | 0.605 [0.523,0.687] | 0.520 | 0.405 | 0.903 | not met |
| combined | 54.138 | 0.565 [0.528,0.683] | 0.523 | 0.411 | 0.738 | not met |
| wdegree | 59.337 | 0.644 [0.531,0.683] | 0.517 | 0.411 | 0.978 | not met |
| fibrosis | 52.256 | 0.574 [0.426,0.655] | 0.526 | 0.412 | 0.768 | not met |
| fibrosis_grad | 69.133 | 0.727 [0.604,0.813] | 0.526 | 0.417 | 1.000 | not met |

### Keep/delete each spectral claim by measured correlation (dossier p.39)

Origin score-rank = fraction of nodes scoring ≥ the reentry origin (small ⇒ origin is a hotspot). One-sided Wilcoxon vs the no-association median 0.5.

| claim (localizer) | mean rank | perm-null mean | **perm p** | Wilcoxon p | verdict |
| --- | --- | --- | --- | --- | --- |
| grad_phi2 | 0.284 | 0.510 | 0.000 | 0.002 | **KEEP** |
| perron | 0.639 | 0.515 | 0.963 | 0.989 | **DELETE** |
| combined | 0.626 | 0.507 | 0.967 | 0.988 | **DELETE** |
| wdegree | 0.725 | 0.505 | 1.000 | 0.999 | **DELETE** |
| fibrosis | 0.326 | 0.531 | 0.000 | 0.011 | **KEEP** |
| fibrosis_grad | 0.716 | 0.515 | 0.997 | 0.997 | **DELETE** |

Interpretation: the strict localization endpoint (argmax below the spatial-null 5th pct) is not met by any field. But the keep/delete test shows the **Fiedler gradient |∇φ₂|** has a genuine spatial association with reentry origins, while **Perron localization is deleted** — so the hypothesized |∇φ₂|∩Perron hotspot is *worse* than |∇φ₂| alone (the Perron factor is uncorrelated). Fibrosis localizes comparably; spectral-radius (a scalar) is not a localizer.
