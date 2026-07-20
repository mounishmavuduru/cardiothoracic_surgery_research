# EXPERIMENTAL_PLAN — earning "generational" for AtrialSpectralBench

The software engine is the *instrument*. This file is the *science*: the ordered
experiments that decide whether the Spectral Fragility Index (SFI) is a
generational contribution or an honest null. Each of the four generational-makers
has a **pre-registered endpoint** and an explicit **decision gate** (what a
negative result means and where the project pivots).

## The thesis (four claims, four makers)

The closed-form Fiedler fragility `dλ₂/dw_ij = (φ_i − φ_j)²`, aggregated as the SFI,
is a surrogate for surgery-triggered reentry that is:

- **GM1 — predictive:** beats the published fibrosis-heterogeneity features head-to-head.
- **GM2 — mechanistically localizing:** its hotspots find where reentry actually starts.
- **GM3 — validity-bounded:** we know exactly when the linear derivative is trustworthy and when it breaks (**the honest core**).
- **GM4 — transferable:** the same fragility calculus works on a second excitable medium (**the generational headline**).

## Ground rules (pre-registration; anti-p-hacking)

1. Endpoints, feature sets, CV scheme, and the Δw stress-field calibration are **frozen in `docs/PRE_REGISTRATION.md` and git-tagged BEFORE looking at any openCARP label.**
2. Always shape-family `GroupKFold`; always report **grouped AND naive** AUC.
3. A **null result is a valid, reported outcome.** Never tune to beat an endpoint.
4. `mock_ep` labels are for plumbing only. **Every scientific claim uses openCARP labels.**
5. No silent caps: `log()` exactly what was and wasn't scaled (cohort size, samples, labels).

## Dependency graph

```
E0 verify engine ─▶ E1 real labels ─▶ ┌─ GM1 predict ─┐
                                      ├─ GM2 localize ─┤─▶ GM4 transfer
                                      └─ GM3 validity ─┘
                    E6 UQ  and  E7 credibility  wrap  GM1–GM3
```

---

## Prerequisites

### E0 — Verify the ground-truth engine  *(gate for everything downstream)*
- Reproduce the **Niederer N-version** cardiac-EP benchmark in openCARP (PMID 21969679): activation times within the published multi-simulator consensus tolerance.
- Analytic spectral gates (path/ring/grid `λ_k = 2−2cos(kπ/N)`; already CI) stay green.
- **Mesh/step convergence:** λ₂ and the inducibility label are stable under mesh refinement and Δt reduction (ASME V&V-40 calculation verification).
- **GATE:** if openCARP fails Niederer, or labels are not mesh-converged → **STOP**. The ground truth is untrustworthy; nothing downstream is interpretable.

### E1 — Generate real ground-truth labels  *(replaces `mock_ep`)*
- **Anatomy:** Rodero cohort (Zenodo 4506930) via atrialmtk → LA graphs with fibres + Universal Atrial Coordinates; fibrosis fields from Roney (Zenodo 5801337).
- **Stress:** calibrate the diffuse Δw uncoupling field to **literature-measured perioperative CV/ERP changes** (cite; freeze the magnitude before labelling).
- **Labels:** openCARP monodomain + Mitchell-Schaeffer, S1–S2 / burst pacing from a few sites, pre- and post-stress → binary `inducible` + reentry-origin coordinate. Start with dozens of anatomies; scale as compute allows; **reuse published Roney-cohort labels where possible.**
- **GATE:** inducible fraction lands in a physiologically plausible band (~10–40%, near the pooled ~26% post-CABG POAF incidence, PMC10823463) and is reproducible. If 0% or 100% → recalibrate the **stress field**, never the endpoint.

---

## The four generational-makers

### GM1 (E2) — PREDICT: SFI beats the real competitors head-to-head  ★ PRIMARY
- **Competitors (the real strawman-killers):** fibrosis burden, fibrosis spatial entropy, patch-size, deterministic min-cut, percolation threshold, and **λ₂-alone**.
- **Design:** nested `GroupKFold` (shape-family held out); classifiers = logistic regression + gradient-boosted trees; compare feature sets `{fibrosis-heterogeneity}` vs `{fibrosis-heterogeneity + SFI}` on the **same openCARP labels**.
- **Pre-registered endpoint:** grouped **ΔAUC ≥ 0.05 with DeLong p < 0.05**; report grouped AND naive AUC (naive exposes PCA-resample leakage).
- **Scale:** as many labeled subjects as compute allows; DeLong bootstrap / null resamples **≥ 10,000**.
- **DECISION GATE:**
  - **WIN** → GM1 achieved; the head-to-head is the spine of the paper. Proceed to GM2/GM4 with confidence.
  - **NULL** (SFI ≡ re-encoded fibrosis) → **report it honestly** and pivot the narrative to *"when and why does spectral fragility collapse to fibrosis burden?"* — a solid, publishable result, **not** generational. GM3 still stands and becomes the lead contribution.

### GM2 (E3) — LOCALIZE: the mechanism is measured, not asserted
- Colocalize the `|∇φ₂| ∩ Perron` hotspot map with openCARP reentry-origin coordinates.
- **Null:** UAC rotational/shift spatial-null; statistic = median geodesic error hotspot→origin, reported against the **null 5th percentile** with a CI.
- Correlate **each** linear-spectral claim (Fiedler gradient, Perron localization, spectral radius) with the nonlinear label and **keep or delete each by its measured correlation.** Drop "spectral radius = reentry" unless it empirically correlates.
- **Endpoint:** median geodesic error **< null 5th percentile**.
- **GATE:** pass → mechanistic localization is generational-maker #2. Fail → drop the localization claim; keep the predictive claim (GM1) only. Do not over-claim a mechanism the data won't support.

### GM3 (E4) — BOUND: the validity radius  *(the honest core — a result either way)*
- Sweep the stress magnitude `‖ΔL‖` and the `λ₂–λ₃` gap; measure agreement between first-order SFI and the **exact** Δλ₂ as a function of `‖ΔL‖ / (λ₃−λ₂)`.
- Show the empirical validity boundary matches the theory: **Weyl** bounds the eigenvalue move, **Davis–Kahan / Bauer–Fike** bound the eigenvector/subspace term. Show the **subspace-SFI auto-switch** restores accuracy near degeneracy, and recompute λ₂ exactly for discrete (Maze) cuts to show where the linear form breaks.
- **Endpoint:** a quantified validity radius (a `‖ΔL‖` threshold) with a stated error bound inside it.
- **NOTE:** this holds **whether or not GM1/GM2 win.** "Exactly when this seconds-on-a-laptop derivative is a trustworthy surrogate for an overnight PDE" is itself generational-grade honesty and is the fallback lead if GM1 nulls.

### GM4 (E5) — TRANSFER: a general spectral-fragility calculus  *(the headline)*
- Port `dλ₂/dw` fragility to a **second excitable medium**: a simplified reaction-diffusion / epilepsy-like network graph with its own induced-instability ground truth.
- Show the same fragility derivative predicts induced instability there, above that medium's own spatial-null.
- **Endpoint:** fragility predicts instability in medium #2 above its null.
- **DEPENDENCY:** atrial core (GM1 or GM3) must hold first. **GATE:** transfers → *"not one biomarker but a fragility calculus for any excitable network whose stability lives in a spectral gap"* (the generational claim). Doesn't → report as heart-specific; still a strong result.

---

## Supporting experiments

### E6 — Uncertainty quantification (descoped)
Morris elementary-effects screening + a Gaussian-process/PCE surrogate over `{conduction velocity, fibrosis threshold, fibre field (rule-based vs DT-MRI atlas), Δw stress magnitude}`. Report the robustness of the GM1/GM2 conclusions to each. **Not** full Sobol.

### E7 — Credibility case (ASME V&V-40 / FDA CM&S framing, honestly scoped)
Literature range-checks (POAF ~26%; LA CV 0.3–1.2 m/s; AUC vs Bieging 0.78 / Varela 0.68); anatomy-realism vs the public 2018 LA Segmentation Challenge; independent-geometry generalization on the Vascular Model Repository; dated notebook; preprint + two-audience poster. **Delete all FDA/VICTRE/"in-silico clinical trial" claims** — this is a surrogate for a simulator's verdict, not clinical POAF.

---

## Scale — the honest answer to "why not 10,000 / 50,000"

- **Code tests:** raw count is not a quality metric. Convert the core invariants to **property-based tests** (`hypothesis`, added as a dev dependency) that generate **thousands of random graphs per invariant** — Laplacian PSD + zero-row-sum, the derivative identity `dλ₂/dw=(φ_i−φ_j)²` vs central finite-difference, the subspace fallback under forced degeneracy, Perron positivity. This is the legitimate "10,000+ checks," and it is cheap. Keep the 89 curated example tests as fast regression anchors alongside it.
- **Scientific N:** scale where it buys statistical power — virtual-cohort size, Monte-Carlo Δw draws for the SFI expectation, spatial-null resamples, bootstrap (**≥ 10,000** each). openCARP labels are the bottleneck, so scale the **cheap** spectral/SFI/feature computations to thousands and use as many openCARP labels as compute + reused Roney labels allow. Log exactly what was capped.

---

## Division of labor

- **Claude Code (networked session):** E0 verification, E1 wiring + **small** real-label batches, build every GM1–GM4 experiment as a **one-command config-driven harness**, add property-based tests, scale the cheap parts (spectral/SFI/MC/null to ≥10⁴), reproduce Niederer.
- **Claude Science (deferred):** the **full** openCARP inducibility sweep (large cohort × pacing × pre/post × stochastic realizations), full Morris/GP UQ, Bayesian hierarchical validation, publication-scale figures.

## The universal sequence (every idea)
**Verify → Reproduce → Innovate → Validate → Quantify uncertainty → Write → Rehearse the interview.** File the SRC/human-data exemption before any experiment. Keep the dated notebook. Budget the last month for the interview — 25 of 100 points.
