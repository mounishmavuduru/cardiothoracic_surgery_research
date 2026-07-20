# Pre-Registration — AtrialSpectralBench

**Status:** frozen before any inducibility labels are generated or inspected.
**Scope:** 100% in silico. Zero patients, zero IRB — only a documented human-data
*exemption* determination for public de-identified geometries, filed before any experiment.
Ground-truth labels are an EP **simulator's inducibility verdict** (openCARP when
available; the `mock_ep` stand-in during development), never clinical POAF.

This document is written *before* results so the hypothesis is falsifiable and the
analysis is not tuned to the outcome. Deviations must be logged, dated, and justified in
`notebooks/lab_notebook.md`.

---

## 1. The pre-registered falsifiable hypothesis

> **The per-region Spectral Fragility Index (SFI) adds incremental inducibility-classification
> AUC and reentry-origin localization *beyond* the published fibrosis-heterogeneity feature
> set** (fibrosis burden, spatial entropy, patch-size, and the min-cut / percolation
> substrate-isolation family).

The mechanism under test is the closed-form Fiedler-value sensitivity
`dλ₂/dw_ij = (φ_i − φ_j)²`, aggregated under a stochastic fibrosis-weighted diffuse
perioperative-uncoupling field Δw into a per-region SFI, with hotspots at
`|∇φ₂| ∩ Perron`.

**Accepted reportable null.** A clean null — **SFI ≡ re-encoded fibrosis**, i.e. SFI adds
no incremental signal beyond the fibrosis-heterogeneity baseline — is an **explicitly
accepted, reportable outcome**. It is stated out loud here so that reporting the null is a
success of the protocol, not a failure to be buried. The circularity concern (edge weights
are reduced in fibrotic nodes *and* the simulator label is driven by that same fibrosis
field, so SFI could merely re-encode fibrosis burden) is precisely what the incremental
test and the held-out DeLong comparison are designed to expose.

---

## 2. Primary endpoint (headline number to beat)

On **shape-family-held-out (GroupKFold)** folds, adding the per-region SFI to the published
fibrosis-heterogeneity feature set:

- **ΔAUC ≥ 0.05** in openCARP-inducibility classification, **AND**
- **DeLong paired test p < 0.05** for that ΔAUC (paired on the held-out predictions of the
  baseline vs baseline+SFI models).

Both conditions must hold for the primary endpoint to be declared met. AUC is reported
**both grouped and naive** so that PCA-resample train↔test leakage is visible rather than
hidden (naive AUC is expected to be optimistic; grouped AUC is the endpoint).

## 3. Secondary endpoint — mechanism localization (falsifiable spatial test)

The SFI hotspot map (`|∇φ₂| ∩ Perron`) localizes the simulator's reentry-initiation site
with **median geodesic/UAC error below the rotational/shift spatial-null's 5th percentile**.
Reported with a bootstrap confidence interval. Each linear-spectral claim is kept or
**deleted** by its measured correlation; "spectral radius = reentry" is dropped unless it
empirically correlates.

---

## 4. Leakage control — shape-family GroupKFold

- Cross-validation is **GroupKFold with the group = base shape family**. Every
  statistical-shape variant (including PCA-resampled near-duplicate anatomies) is assigned
  to its **parent base family's group**, so no shape family appears in both train and test.
- Both **grouped AND naive** AUC are reported; the gap between them quantifies the leakage
  that naive resampling would otherwise inflate.
- The classifier hyperparameters are chosen in a **nested** inner GroupKFold loop so tuning
  never sees the outer test fold.
- The λ₂−λ₃ spectral-gap distribution across the cohort is reported; when the gap is small
  the single Fiedler vector is ill-defined and analysis auto-switches to the
  **spectral-projector (subspace) SFI**.

## 5. Real competitor baselines (kill the strawman)

The baseline feature set the SFI must beat is composed of published, non-trivial competitors —
not a strawman:

1. **Fibrosis burden** — global fibrotic fraction.
2. **Fibrosis spatial entropy** — heterogeneity of the fibrosis field.
3. **Fibrosis patch size** — mean connected fibrotic-patch size.
4. **Deterministic global min-cut** — Stoer–Wagner substrate-isolation value
   (the spectral cousin of the Cheeger/λ₂ story; guards against "you renamed min-cut").
5. **Percolation threshold** — edge-removal fraction at which the giant component breaks.
6. **λ₂-alone** — the Fiedler value on its own, without the per-region SFI aggregation
   (isolates what the *SFI construction* adds over the bare algebraic connectivity).

SFI is evaluated strictly as an **add-on / ablation** to this set, on the **same** simulator
labels, with the incremental ΔAUC (Section 2) as the test statistic.

---

## 6. Analysis integrity commitments

- All stochastic steps are seeded (`rng` / `seed`); the exact config
  (`configs/default.yaml`) and commit hash accompany every reported number.
- The primary and secondary endpoints, the CV scheme, the baseline set, and the accepted
  null above are **fixed before labels are inspected**.
- Claims are scoped to "a cheap spectral surrogate for a specific EP-simulator's
  inducibility verdict." No FDA / VICTRE / ISCT language; nothing here is validated against
  real post-operative AF. The 2018 LA Segmentation Challenge geometries are used **only** as
  an anatomy-realism check, never as clinical validation.
