# GM4 — TRANSFER: spectral-fragility calculus on a second excitable medium

> Second medium: a 2-D FitzHugh-Nagumo neural network with an epileptic-focus (hyperexcitable-lesion) instability. FHN labels are simulator verdicts, not clinical. The spectral/SFI/baseline code is reused **verbatim** from the cardiac pipeline.

- Networks **80**, unstable **42** (0.5250), each network its own group.

## Part A — predict FHN instability (GM1-analog, grouped)

| add-on | clf | grouped base | grouped +SFI | ΔAUC | DeLong p | endpoint |
| --- | --- | --- | --- | --- | --- | --- |
| competitors vs +single_vector | lr | 0.8878 | 0.8941 | 0.0063 | 0.7420 | not met |
| competitors vs +single_vector | gbt | 0.8753 | 0.8997 | 0.0244 | 0.3202 | not met |
| lesion vs +single_vector | lr | 0.8878 | 0.9010 | 0.0132 | 0.5286 | not met |
| lesion vs +single_vector | gbt | 0.8634 | 0.8985 | 0.0351 | 0.2105 | not met |
| competitors vs +subspace | lr | 0.8878 | 0.8628 | -0.0251 | 0.3882 | not met |
| competitors vs +subspace | gbt | 0.8753 | 0.8590 | -0.0163 | 0.4587 | not met |
| lesion vs +subspace | lr | 0.8878 | 0.8766 | -0.0113 | 0.6011 | not met |
| lesion vs +subspace | gbt | 0.8634 | 0.8640 | 0.0006 | 0.9827 | not met |
| competitors vs +exact | lr | 0.8878 | 0.8847 | -0.0031 | 0.8775 | not met |
| competitors vs +exact | gbt | 0.8753 | 0.8628 | -0.0125 | 0.2429 | not met |
| lesion vs +exact | lr | 0.8878 | 0.8678 | -0.0201 | 0.3802 | not met |
| lesion vs +exact | gbt | 0.8634 | 0.8459 | -0.0175 | 0.3421 | not met |

## Part B — localize the instability origin (GM2-analog, N=42)

| localizer | mean origin-rank | perm-null | perm p | verdict |
| --- | --- | --- | --- | --- |
| grad_phi2 | 0.5456 | 0.5037 | 0.8333 | **DELETE** |
| perron | 0.3622 | 0.4974 | 0.0033 | **KEEP** |
| combined | 0.3633 | 0.5012 | 0.0000 | **KEEP** |
| fibrosis | 0.5991 | 0.5045 | 0.9933 | **DELETE** |
| fibrosis_grad | 0.5212 | 0.5038 | 0.6500 | **DELETE** |

## Part C — validity radius on the network cohort (GM3-analog)

- ρ = ‖ΔL‖/(λ3−λ2): median **1852.4387** (min 908.6008, max 7245.5859) over 24 networks.
- Synthetic validity boundary ρ* (first-order 10% error) ≈ 2.9998846240069104.
