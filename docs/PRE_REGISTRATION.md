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

---

## 7. Real-data addendum — FROZEN 2026-07-20, before any SFI-vs-label analysis (GM1)

This section pins every real-data choice **before** the SFI↔inducibility relationship is
inspected. It is git-tagged (`prereg-frozen-YYYYMMDD`). Any later change is a logged,
dated deviation in `notebooks/lab_notebook.md`. The frozen synthetic `mock_ep` endpoints of
Sections 1–6 are unchanged; this addendum instantiates them on real anatomy.

### 7.1 Deviation logged: ground-truth simulator
openCARP is **not installable** on this CPU-only ephemeral container (no apt/conda/pip
package; a PETSc source build is out of scope) and the user defers the full openCARP sweep
to Claude Science. The Roney archive (Zenodo 5801337) ships **no** inducibility label files,
so published labels cannot be reused directly. Ground-truth labels are therefore generated
by a **genuine monodomain Mitchell–Schaeffer reaction–diffusion solver**
(`asb.labels.monodomain`, `source='monodomain_ms'`) — exactly plan §8.4's "monodomain +
phenomenological Mitchell–Schaeffer on coarsened meshes". This is a *simulator verdict*,
in the same category as openCARP, and is **never** clinical POAF. openCARP stays the wired,
one-command deferred target (`asb.labels.opencarp`). Reentry here is governed by
wavelength (CV×ERP), source–sink mismatch and unidirectional block — **never by λ₂** — so a
predictive SFI→label link is not circular.

### 7.2 Anatomy + substrate (real, frozen)
- **Cohort:** Roney LA virtual cohort, Zenodo **5801337** (patient-derived LA surfaces).
  Each `Mesh_<id>.vtk` carries `UAC1/UAC2` (∈[0,1]), `IIR` (LGE fibrosis proxy) and
  `fiber_endo`/`fiber_epi`. Coordinates microns → mm. Coarsened by grid vertex-clustering
  to **~2000 nodes**, carrying all fields.
- **Fibrosis map (frozen):** `fibrosis = clip((IIR − 1.0)/(1.32 − 1.0), 0, 1)` — the
  Khurram 2014 IIR healthy(1.0)→dense-scar(1.32) ramp.
- **Grouping (CV):** each patient mesh is its own `shape_family`, so **GroupKFold =
  patient-held-out** (the natural leakage-free grouping). Grouped **and** naive AUC both
  reported.

### 7.3 Δw diffuse uncoupling field (frozen, literature-calibrated)
Perioperative inflammation slows atrial CV ~**20 %** (best estimate; range 10–30 %; Heida
2021 PMID 34198544 regional; Clayton 2018 remodelled). With monodomain CV∝√D, a fractional
CV drop `f_CV` implies a fractional conductance drop `f_D = 1 − (1 − f_CV)²`. Frozen:
- `delta_w_mean_frac = 0.36`  (from `f_CV = 0.20` → `f_D = 1 − 0.8² = 0.36`; band 0.19–0.51),
- `delta_w_cov = 0.5` (heterogeneity — the perioperative effect is dispersion, not uniform),
- `fibrosis_coupling = 1.0` (larger drops in fibrotic tissue),
- `validity_safety = 0.25` (first-order guard `‖ΔL‖ ≤ 0.25·(λ₃−λ₂)`).
Δw = 0.36 is a **large** perturbation, likely beyond the first-order validity radius — the
SFI feature is used as a *ranking* feature and the exact validity boundary is mapped in GM3;
sensitivity to Δw over 0.19–0.51 is reported (E6).

### 7.4 Monodomain-MS ground-truth protocol (frozen)
Membrane (Mitchell & Schaeffer 2003, atrial-tuned): `tau_in=0.3, tau_out=6, tau_open=120,
tau_close=110, v_gate=0.13` ms; fibrosis shortens `tau_close` by up to 50 %
(`fibrosis_erp_shortening=0.5`). Diffusion: anisotropic conductance-weighted cotangent
Laplacian, `d0=0.20` (healthy along-fibre CV ≈ 0.87 m/s), `w_ref=0.3`, implicit
backward-Euler, `dt=0.05` ms. Induction battery: **burst pacing**, cycle length **150 ms**,
6 beats, from **2** random sites; inducible iff self-sustained supra-threshold activity
persists **≥ 650 ms** after the last stimulus (observation window 1000 ms). Verified: APD90
218/121 ms (healthy/fibrotic), planar CV 0.4–1.2 m/s, monotone reentry vs fibrosis.
**Calibration gate PASSED:** inducible fraction **0.33** (10–40 % band, near 26 % POAF);
Spearman(fibrosis, sustained-reentry) = 0.75. The battery deliberately excludes the
half-field cross-field inducer (anchors substrate-independent geometric reentry).

### 7.5 GM1 head-to-head (frozen)
- **Competitor baseline set** (must be beaten): `fibrosis_burden`, `fibrosis_spatial_entropy`,
  `fibrosis_patch_size`, `min_cut_value`, `percolation_threshold`, `lambda2_alone`.
- **Feature sets compared:** `{competitors}` vs `{competitors + per-region SFI}` (SFI as a
  strict add-on/ablation), and additionally `{fibrosis-heterogeneity}` vs `{+SFI}`.
- **Classifiers:** logistic regression **and** gradient-boosted trees, nested (inner
  GroupKFold model selection never sees the outer test fold).
- **Primary endpoint (unchanged):** grouped **ΔAUC ≥ 0.05 AND DeLong p < 0.05**; report
  grouped **and** naive AUC; bootstrap/null resamples **≥ 10 000**.
- **Accepted null (unchanged):** SFI ≡ re-encoded fibrosis is an explicitly reportable
  outcome; we will not tune to beat the endpoint.

---

## 8. Real *clinical outcome* addendum — FROZEN 2026-07-24, before any label-linked analysis

**What changed.** On 2026-07-24 Prof. Patrick M. Boyle (UW) shared per-patient 2-year
post-ablation arrhythmia-recurrence outcomes for the 82-patient UW cohort whose LA meshes
are public on Dryad (`10.5061/dryad.kkwh70sg0`). `docs/DATASETS.md:49` had recorded these
labels as **withheld**; that is no longer true. Sections 1–7 govern *simulator* endpoints
and are unchanged. This section governs the **clinical** endpoint and is written before any
feature↔outcome association has been computed, inspected, or plotted.

**Integrity statement.** At the time of freezing, the only properties of the outcome column
that have been examined are its **marginal counts** and their split across the publisher's
predefined cohorts (below). Marginals are required to compute power and are not a
feature–outcome association. No SFI value, no fibrosis value, and no mesh-derived quantity
has been placed alongside the outcome column.

### 8.1 Data provenance and custody
- **Labels:** `Patient_ID, Recurrence_Rhythm_2yr` ∈ {`NR`, `AF`, `AFL`}, n = 82,
  IDs 1–15 and 21–87. `NR` = no recurrence within the 2-year follow-up; `AF`/`AFL` = first
  documented recurrence rhythm was fibrillation / flutter respectively.
- **Anatomy:** the matching public Dryad meshes (`ID001–ID015` holdout, `ID021–ID087`
  training; pre- and post-ablation per patient). The ID sets match the label file exactly.
- **Custody:** the label file is **human-subject outcome data shared under restricted
  terms**, not covered by the CC0 that applies to the meshes. It lives under the gitignored
  `data/` tree and is additionally protected by explicit `.gitignore` rules. It **will not be
  committed, redistributed, or included in any release artifact** unless Prof. Boyle
  explicitly permits it in writing. Derived aggregate statistics may be published; the
  per-patient column may not.

### 8.2 Observed marginals (the only thing inspected pre-freeze)

| Cohort | n | NR | AF | AFL | recurrence |
|---|---|---|---|---|---|
| All | 82 | 34 | 35 | 13 | 48 (58.5 %) |
| Dryad holdout `ID001–015` | 15 | 7 | 6 | 2 | 8 (53.3 %) |
| Dryad training `ID021–087` | 67 | 27 | 29 | 11 | 40 (59.7 %) |

### 8.3 Power — computed from these marginals *before* the analysis
Paired binormal ROC simulation (4 000 replicates, DeLong two-sided, α = 0.05, n₁ = 48,
n₀ = 34), sweeping the shared-noise correlation `r` between the nested baseline and
baseline+SFI scores. Power to detect an incremental ΔAUC:

| base AUC | r | Δ=+0.03 | Δ=+0.05 | Δ=+0.07 | Δ=+0.10 | Δ=+0.15 |
|---|---|---|---|---|---|---|
| 0.70 | 0.80 | 0.11 | 0.24 | 0.43 | 0.75 | 0.99 |
| 0.75 | 0.80 | 0.12 | 0.27 | 0.49 | 0.82 | 0.99 |
| 0.70 | 0.50 | 0.09 | 0.13 | 0.22 | 0.41 | 0.79 |
| 0.75 | 0.50 | 0.08 | 0.14 | 0.26 | 0.48 | 0.87 |

**Consequence, stated in advance.** At n = 82 this study has roughly **25 % power** at the
ΔAUC ≥ 0.05 threshold used for the simulator endpoint in §2/§7.5. That threshold is
therefore **not reused here** — at this sample size it would be effectively unfalsifiable,
and a null against it would carry almost no information. The **minimum detectable effect is
ΔAUC ≈ 0.10** (optimistic, strongly correlated nested models) to **≈ 0.15** (conservative).
A null result below that band will be reported as **inconclusive — underpowered**, never as
evidence that SFI adds nothing. This is fixed now precisely so it cannot be renegotiated
after seeing the answer.

### 8.4 Primary endpoint (clinical)
Binary outcome `recurrence = (rhythm ≠ NR)`, i.e. AF and AFL both count as events — matching
the source study's definition of recurrent atrial arrhythmia.

Features are computed from the **pre-ablation mesh only** (the clinically useful question is
prediction from baseline substrate; the post-ablation mesh encodes the delivered treatment).

- **Baseline set:** the §7.5 competitors recomputed on UW anatomy — `fibrosis_burden`,
  `fibrosis_spatial_entropy`, `fibrosis_patch_size`, `min_cut_value`,
  `percolation_threshold`, `lambda2_alone`.
- **Test set:** baseline **+ per-region SFI** as a strict add-on.
- **Estimation:** patient-held-out cross-validation over all 82 (each patient its own group),
  nested model selection, logistic regression and gradient-boosted trees, ≥ 10 000 bootstrap
  resamples. All 82 are used for the primary because §8.3 shows the study is already at the
  edge of detectability; splitting further would forfeit power we cannot spare.
- **Declared met iff:** ΔAUC ≥ 0.10 **and** DeLong paired p < 0.05 (Holm-adjusted, §8.7).
- SFI parameters remain the §7.3 frozen values. SFI is closed-form and never fitted to the
  outcome, so it contributes no tuning degrees of freedom.

### 8.5 Pre-specified secondary analyses (fixed list, no additions later)
1. **Publisher's split, honoured.** Fit on `ID021–087` (n = 67), evaluate **once** on the
   untouched `ID001–015` (n = 15). With 8 events this is **descriptive only** — reported with
   a confidence interval and explicitly not powered. Its value is that the split was chosen
   by someone else, before we existed.
2. **Ablation-induced spectral change.** Δ(SFI) = post-ablation − pre-ablation, testing
   whether the lesion set's effect on spectral fragility predicts recurrence. This is the
   SFI-native analogue of the source study's finding that post-ablation substrate carries
   signal.
3. **Mechanism split (exploratory).** Among the 48 recurrers, does SFI separate `AF` (35)
   from `AFL` (13)? SFI indexes distributed reentry vulnerability, whereas flutter is
   typically macro-reentrant and often right-atrial — so a *negative* result here is
   mechanistically informative. Labelled exploratory; n = 13 in the smaller arm.

### 8.6 Known substrate limitations of the UW meshes — disclosed before analysis
- **Fibrosis is categorical, not continuous.** UW files carry a per-*cell* `elemTag`
  (`{111:0.0, 115:1.0, 164:0.5, 199:1.0}`, `uw_boyle.py:68`) averaged onto vertices, whereas
  the Roney cohort carries continuous LGE `IIR`. The fibrosis competitors are therefore
  **measured more coarsely here than in §7.5**.
- **This biases the comparison toward SFI**, because a degraded baseline is easier to beat.
  Mitigation, fixed now: we will report the absolute AUC of the fibrosis-only baseline and
  compare it against the source study's published fibrosis-derived performance. If our
  baseline is materially weaker than theirs, any ΔAUC we observe will be reported as
  **confounded by baseline degradation**, not as a win for SFI.
- **UAC is a PCA surrogate**, not anatomical (`uw_boyle.py:180-192`), so any region-based or
  localization claim on this cohort is weaker than on Roney and will be labelled as such.
- **Amendment 2026-07-24 (same day, still before any label join): the tag map is probably
  wrong, and this blocks the analysis.** Having downloaded the meshes, they contain a per-cell
  `elemTag` and a `VECTORS fiber` array, and no continuous LGE field.
  *(Correction, same day: an earlier draft of this amendment said there was no fibre array.
  That was wrong — the array is present; it is its **contents** that are degenerate, see the
  final bullet below.)* Measured over all 164 meshes: tag 111 59.0 % → 47.3 % (pre → post ablation), 115 20.9 % →
  14.9 %, 199 absent → 17.9 %, and **164 20.0 % → 19.9 %**. Tag 164 is untouched by ablation,
  which is not fibrosis behaviour; it is far more consistent with a non-myocardial structure
  (mitral annulus / PV sleeves). The loader nevertheless assigns it `fibrosis = 0.5`
  (`uw_boyle.py:68`), injecting a near-constant ~20 % pseudo-fibrosis into every patient —
  which both inflates `fibrosis_burden` and *shrinks its between-patient variance*, weakening
  the very baseline SFI must beat, and corrupts the Δw field SFI is computed from.
  **The confirmatory analysis will not be run until the tag semantics are confirmed with the
  data provider** (`docs/outreach/boyle_reply_followup.md`, Q2). This correction is being made
  blind to the outcome column, which has still never been joined to any feature.
- **Amendment 2026-07-24 (b): the released fibre field is degenerate.** All 164 meshes carry a `VECTORS fiber`
  array whose value is a constant `(1, 0, 0)` for every element — measured mean directional
  spread is exactly `0.0`. `edge_weights_from_fibres` makes conduction anisotropic
  (along 1.0 / cross 0.3) *relative to the local fibre direction*, so a single global
  direction reduces the anisotropy to a fixed coordinate bias with no anatomical content and
  removes all fibre heterogeneity — a principal substrate for unidirectional block and hence
  for reentry initiation. `uw_boyle.py` now measures the spread, records
  `fibres_are_degenerate` in mesh metadata, and raises a `RuntimeWarning`, so the condition can
  no longer pass unnoticed. Note this is **upstream data**, not a local bug: the defect is in
  the public Dryad deposit. It does not affect the §8 clinical endpoint, which uses real
  outcomes and needs no simulator label — but it does affect every existing in-silico result
  that includes UW subjects.
- **Amendment 2026-07-24 (c): measured attribution — both defects matter, neither explains the
  anomaly.** `scripts/uw_substrate_ablation.py` ran the 2×2 (orifices sealed/open × fibres
  constant/varying) over all 82 subjects under the otherwise-frozen protocol;
  `scripts/analyse_substrate_ablation.py` does the paired statistics. Results:

  | | constant fibres | varying fibres |
  |---|---|---|
  | orifices sealed (as released) | **6/82 (7.3 %)** | 8/82 (9.8 %) |
  | orifices opened | 8/82 (9.8 %) | **14/82 (17.1 %)** |

  The sealed/constant arm **reproduces the recorded 6/82 exactly**, which validates the setup.
  Each defect alone converts ~4 subjects; together they convert 8, so the two are
  **super-additive** (interaction +4) — consistent with reentry needing both an anatomical
  obstacle and anisotropic heterogeneity to anchor.

  Three cautions are recorded so that no stronger claim is made later than the data supports:
  (i) **the change is not significant** — the strongest contrast, sealed/constant vs
  open/varying, gives McNemar exact **p = 0.077**; (ii) it closes only **39 %** of the gap to
  Roney's 32.3 %, so most of the discrepancy is still unexplained, with the categorical
  three-level fibrosis (against Roney's continuous 276–1071-level IIR) the obvious remaining
  suspect; (iii) **only 2 of the original 6 positives survive the repair** while 12 new ones
  appear, i.e. per-subject verdicts churn heavily — consistent with the labeller instability
  already documented in the mesh-convergence study (only 7/24 subjects consistent across six
  resolution tiers). An earlier draft of amendment (b) called the fibre defect "the leading
  explanation" for the 7 %. The measurement does not support that phrasing and it has been
  withdrawn here and in the manuscript.
- **Amendment 2026-07-24 (d): census confirming the tag-164 identification, and a second
  self-correction.** The identification was initially made on one mesh. Repeated over all 82:
  tag 164 resolves into **four to six large components** (five in 75 meshes, four in five,
  six in two — the usual pulmonary-vein variation, a left common trunk giving four and a
  right middle vein six), each a topological disc; it borders fibrotic elements on **~0.08 %
  of its incident edges** against **14–20 % by chance**; and removing it opens the surface by
  a median of **7** boundary loops. Two meshes, **ID040 and ID049**, are topologically
  pathological (χ before removal −16 and −18).

  Earlier same-day drafts of this section and of `notebooks/lab_notebook.md` claimed "exactly
  five components", "**zero** adjacency to fibrosis" and "χ = −3, exactly five openings".
  Those were generalisations from mesh ID001 and from a mean fraction rounded to three
  decimals; all three are withdrawn. The identification itself stands and is better supported
  than before, since the component spread matches known anatomy. Guards are pinned in
  `tests/test_uw_boyle.py`, which is what caught the error.
- **Amendment 2026-07-25 (f): the third and last candidate is refuted too, in the wrong
  direction.** `scripts/roney_fibrosis_quantization.py` destroyed the fibrosis *gradation*
  on Roney while holding anatomy, fibres and every solver parameter fixed:

  | arm | inducible | rate | mean burden | distinct levels |
  |---|---|---|---|---|
  | continuous, as shipped (control) | 20/62 | 32.3 % | 0.333 | 1318 |
  | binary at f > 0.5 | 31/62 | 50.0 % | 0.290 | 2 |
  | binary, burden held fixed | **40/62** | **64.5 %** | 0.333 | 2 |

  The control reproduces 20/62. The prediction was that coarsening fibrosis would *suppress*
  reentry by removing graded border zones, and so explain the UW cohort's low rate.
  It does the opposite: coarsening roughly **doubles** inducibility, and the burden-matched
  arm — which changes nothing but the gradation — is the strongest of the three. A sharp
  binary boundary is evidently a stronger source–sink discontinuity than a smooth ramp.

  **Consequence.** All three candidate explanations for the UW anomaly are now dead:
  sealed orifices (+2/82, n.s.), degenerate fibres (no effect on the rate, p = 1.0), and
  coarse fibrosis (wrong direction). Two of the three, imposed deliberately on a clean
  cohort, would *raise* the rate. The anomaly is unexplained and is reported as an open
  problem. No further explanation will be offered in any document without a control
  behind it.

- **Amendment 2026-07-25 (g): openCARP is now installable, and the substitution is no longer
  forced.** §7.1 recorded openCARP as uninstallable on the then-current container, which is
  why the monodomain Mitchell–Schaeffer solver was substituted. On the current machine a WSL
  Ubuntu 26.04 environment (8 cores, 15 GB) is available, and the openCARP v19.0 AppImage
  runs there without root once extracted (`--appimage-extract`; FUSE is present but
  `libfuse.so.2` is not). `openCARP -buildinfo` returns cleanly. This is recorded as
  *feasibility only*. Switching the ground-truth labeller is a major protocol change that
  would invalidate every in-silico number in the manuscript, and it is not being made
  mid-analysis. It is the single highest-value change available, given amendment (e)'s
  finding that per-subject labels are noisy, and it should be scheduled deliberately with
  its own re-validation.

  **Validation completed 2026-07-25** (`scripts/opencarp_validate.sh`, reproducible):
  the binary runs; **MitchellSchaeffer is among the 58 available ionic models**, i.e. the
  *same membrane model* as `asb.labels.monodomain`, so a switch is like-for-like rather
  than a change of physics; single-cell gives **APD90 = 246 ms** against the in-house
  solver's 218 ms; and a tissue run propagates a planar wave across a 10,201-node slab
  (0 → 544 → 1526 → … → 10,201 activated nodes). Two gotchas are recorded so the next
  run does not rediscover them: the AppImage must be `--appimage-extract`ed because
  `libfuse.so.2` is absent even though `fusermount` is present, and MitchellSchaeffer is
  a **normalised** model (Vm ≈ 0–1), so a physiological 250 µA/cm² stimulus makes the
  parabolic solve diverge with NaN — 60 works.

  **One design risk is flagged now, before it is designed in.** openCARP assigns
  conductivity through discrete `gregion` element tags, so a continuous fibrosis field
  must be binned to be used. Amendment (f) established that binarising fibrosis roughly
  *doubles* inducibility on this pipeline. A coarse binning would therefore import
  precisely the artefact we just characterised. Any openCARP labeller must use enough
  bins to approximate the continuous field, and must demonstrate insensitivity to the bin
  count before its labels are used for anything.
- **Amendment 2026-07-24 (e): the fibre hypothesis is REFUTED by direct control, and the
  label is noisy per subject.** `scripts/roney_fibre_control.py` destroyed only the fibre
  field on the Roney cohort — which ships a real one — holding anatomy, fibrosis and every
  solver parameter fixed:

  | arm | inducible | rate | mean fibre spread |
  |---|---|---|---|
  | real fibres (control) | **20/62** | 32.3 % | 0.954 |
  | fibres set to a constant | **20/62** | 32.3 % | 0.000 |

  The control arm reproduces the recorded 20/62 exactly. The treatment arm is *identical in
  count* (McNemar b=5, c=5, **p = 1.0**). Fibre degeneracy therefore **does not drive the
  inducibility rate**, and amendment (b)'s framing of it as the explanation for the UW
  cohort's 7 % is withdrawn in full. It remains a genuine data-quality defect in the public
  deposit, worth reporting to the provider, but it is not the cause of the anomaly. The
  +4-subject fibre effect seen in the UW 2×2 (amendment (c)) is best read as the same churn
  described below rather than as a fibre effect; it was already non-significant (p = 0.69).

  **The more consequential finding is what the control exposes about the labeller.** Removing
  the fibre field changes the *dynamics* substantially — sustained-reentry duration moves on
  50/62 subjects, by up to 661 ms — and flips 10/62 individual verdicts, yet leaves the
  aggregate rate untouched. Combined with the UW repair retaining only 2 of 6 positives, and
  with the mesh-convergence study's 7/24 subjects consistent across six resolutions, the
  monodomain-MS inducibility label is **stable in aggregate and unstable per subject**. Since
  the §2/§7.5 predictive endpoint is scored on per-subject labels, non-differential label
  noise biases ΔAUC toward zero. The in-silico null must therefore be reported as partly a
  statement about the labeller, not only about SFI. This does **not** affect the §8 clinical
  endpoint, which uses real outcomes and no simulator label.
- **Benchmark context.** The source study reports ROC AUC **0.80 ± 0.04** using 89 features
  including EHR/clinical risk factors that we do not hold. We are **not** claiming to beat
  that model, and will not present our mesh-only AUC as if it were comparable.

### 8.7 Multiplicity, stopping rule, and accepted null
- **Three** confirmatory tests are pre-specified (§8.4 primary, §8.5.1, §8.5.2). Holm
  correction across those three. §8.5.3 is exploratory and reported unadjusted, labelled.
- **The confirmatory analysis is run once.** The pipeline will be built, unit-tested, and
  dry-run end-to-end on *shuffled* labels before the real column is ever joined. Any change
  after the real join is a logged, dated deviation in `notebooks/lab_notebook.md` and
  demotes the affected test to exploratory.
- **Accepted reportable null (unchanged in spirit).** "SFI adds no incremental clinical
  predictive value beyond fibrosis and geometry" is an explicitly accepted, publishable
  outcome. Given §8.3, the most likely honest verdict is *inconclusive at this sample size*,
  and we commit in advance to reporting exactly that rather than reaching for a subgroup
  that reaches significance.
