# AtrialSpectralBench: a pre-registered falsification of a spectral fragility biomarker for atrial reentry, and the validity-radius law that explains why it fails

*Working preprint draft. All numbers trace to `results/*.json` and `notebooks/lab_notebook.md`.
This is an in-silico computational study; labels are electrophysiology-simulator verdicts, never
clinical outcomes.*

## Abstract

We ask whether a closed-form spectral graph biomarker — the **Spectral Fragility Index (SFI)**,
the Fiedler-value sensitivity `∂λ₂/∂w_ij = (φ_i − φ_j)²` aggregated under a fibrosis-weighted
conduction-uncoupling field — predicts reentry inducibility in atrial tissue beyond established
fibrosis and connectivity features. Pre-registered before any labels were seen, the primary
endpoint was ΔAUC ≥ 0.05 with DeLong p < 0.05 on shape-family-held-out folds. **The endpoint is not
met.** On two independent real left-atrial cohorts (Roney virtual cohort, n=62; UW/Boyle LGE-MRI
cohort, n=82) and on up to 100,000 synthetic excitable networks, SFI adds no practically meaningful
predictive value (ΔAUC ≈ +0.003; a subspace variant reaches p<10⁻⁴ only at N≥10⁴, an effect so small
it would require ~4,150 patients to detect at 80% power — clinically undetectable). We explain the
null mechanistically: the first-order SFI is valid only while `ρ = ‖ΔL‖/(λ₃−λ₂) = O(1)`, and real
excitable media violate this by ~10³ (atrial median ρ ≈ 2422). We show the spectral gap collapses by
a Weyl scaling law (real atrial surface: gap ∝ N⁻¹·⁰⁵, matching the 2-manifold prediction), so ρ
grows lawfully with resolution — the biomarker fails *predictably*, and *more* as the medium is
meshed finer. A single component survives: the Fiedler-gradient localizer `|∇φ₂|` locates the
instability origin above a spatial null (p<10⁻³) and **transfers across three dynamical media**
(cardiac monodomain, FitzHugh–Nagumo excitable networks, Kuramoto oscillators) — though it is not
uniquely better than substrate imaging, for a reason we make precise. Finally, we demonstrate that
naive methodology (random cross-validation, p-value thresholds, no effect-size gate) would have
reported SFI as a working biomarker; our pre-registered protocol correctly rejects it. The
contribution is a rigorously falsified biomarker, a validity criterion and scaling law that explain
the failure and generalize to any spectral-graph fragility marker, and a worked demonstration of the
false positives standard methodology produces.

## 1. Introduction

Algebraic connectivity `λ₂` (Fiedler 1973) governs the stability of coupled excitable and oscillator
networks (the master-stability-function result, Pecora & Carroll 1998; synchronizability eigenratio,
Barahona & Pecora 2002) and has been used as a real-time seizure detector (Bomela et al. 2020). The
per-edge sensitivity `∂λ₂/∂w_ij = (φ_i − φ_j)²` is a classical identity underlying Fiedler-vector
edge-editing heuristics (Ghosh & Boyd 2006). It is therefore natural to hypothesize that aggregating
this sensitivity over a fibrosis-weighted perioperative uncoupling field yields a per-region
biomarker of atrial-fibrillation reentry fragility. We test this hypothesis and, crucially,
pre-register the null — "SFI merely re-encodes fibrosis" — as an accepted, reportable outcome, so
that a negative result is a success of the protocol rather than a file-drawer casualty.

The mathematics we use is classical; we claim no new theorem. What is ours is (i) a leakage-controlled
incremental test of SFI against real competitors and a spatial null on real anatomy; (ii) a validity
criterion and Weyl scaling law that explain the failure and predict it for the whole class of spectral
fragility biomarkers; (iii) a cross-medium transfer analysis isolating the one component that
generalizes; and (iv) a falsification-protocol demonstration.

## 2. Methods (summary; full detail in `docs/SFI_THEORY.md`, `docs/PRE_REGISTRATION.md`)

**Substrate.** Real left-atrial surfaces are loaded into weighted conduction graphs (nodes = mesh
vertices, edges = triangle edges, weights = fibre-anisotropic conductance attenuated by local
fibrosis) from two public cohorts: the Roney virtual cohort (Zenodo 5801337) and the UW/Boyle
LGE-MRI cohort (Dryad 10.5061/dryad.kkwh70sg0, 82 distinct AF patients). Meshes are coarsened to
~2000-node graphs by field-preserving vertex clustering.

**Labels.** Ground truth is a monodomain Mitchell–Schaeffer reaction–diffusion inducibility verdict
(self-sustained reentry ≥ 600 ms after pacing), run on each real anatomy. These are simulator
verdicts, not clinical POAF, and not the deferred openCARP gold standard.

**SFI.** Per-region `SFI(R) ≈ Σ_{(i,j)∈R} E[Δw_ij](φ₂,ᵢ − φ₂,ⱼ)²` (first order); a basis-independent
subspace/projector variant is used near eigenvalue degeneracy; discrete edits are recomputed exactly.

**Evaluation protocol (the guards).** Pre-registered: shape-family GroupKFold (grouped *and* naive
AUC reported), DeLong paired test (Sun & Xu 2014 fast midrank), ≥10⁴ stratified bootstrap, a
rotational/torus-shift spatial null for localization, tie-robust permutation nulls, and a keep/delete
rule deciding each sub-claim by *measured* correlation. Effect sizes are gated, not just p-values.

**Transfer media.** The identical calculus is applied to FitzHugh–Nagumo excitable networks and
Kuramoto phase-oscillator networks (a different dynamical class), where `λ₂` governs synchronizability.

## 3. Results

### 3.1 The predictive null (primary endpoint not met), replicated across cohorts and scale

| cohort / medium | N | SFI ΔAUC vs full competitors | verdict |
|---|---|---|---|
| Roney (real) | 62 | −0.030 / +0.036 (lr/gbt), p≥0.26 | null |
| UW/Boyle (real, external) | 82 | −0.057 / +0.104, p≥0.29 | null (underpowered, 6 inducible) |
| Combined real | 144 | −0.004 / +0.051, p≥0.16 | null |
| FHN networks | 2,000 | ≤ +0.003, p≥0.12 | null |
| FHN networks | up to 100,000 | subspace ΔAUC +0.0028 (p<10⁻⁴ at N≥10⁴) | null (effect negligible) |
| Kuramoto networks | 500 | ≤ +0.0036, p≥0.38 | null |

The pre-registered ΔAUC ≥ 0.05 endpoint is not met anywhere. A secondary comparison shows SFI *does*
beat a fibrosis-heterogeneity-only baseline (combined GBT ΔAUC +0.103, p=0.009), so **SFI is not
merely re-encoded fibrosis** — it carries connectivity information beyond the substrate — **but that
information is redundant with existing connectivity features** (λ₂-alone, min-cut, percolation), which
is why it adds nothing over the full competitor set. The precise statement: *SFI ≡ connectivity
information already available from standard graph features, not ≡ fibrosis.*

**Power analysis (`results/power_analysis.json`).** The subspace effect (+0.0028 AUC) requires
~4,150 cases for 80% power; the single-vector effect (+0.0007) requires ~23,000. Real AF imaging
cohorts are ~200–1,000. The effect is statistically real at N=10⁵ and clinically undetectable — the
distinction between statistical and practical significance, quantified.

### 3.2 The validity radius (why it fails) — `docs/SFI_THEORY.md` §4, §9

The first-order SFI's relative error scales as `Θ(ρ)`, `ρ = ‖ΔL‖/(λ₃−λ₂)` (second-order
Rayleigh–Schrödinger; the Fiedler vector rotates by `≤ √2·ρ`, Davis–Kahan). Empirically the error
crosses 10% at ρ* ≈ 3. On real excitable media the gap is tiny, so **ρ is enormous — atrial median
ρ ≈ 2422, neural ρ ≈ 1852**, ~10³× past ρ*. The linear SFI is therefore *provably* outside its
validity radius on real anatomy: the null is mechanistically necessary, not incidental. (Weyl's
inequality, Davis–Kahan, and second-order perturbation theory are classical; applying ρ as an
*a priori* biomarker-validity screen and measuring its violation on real tissue is the contribution.)

### 3.3 The ρ scaling law — `results/rho_scaling.json`, `docs/SFI_THEORY.md` §10

The gap collapses with system size by Weyl's law:

| medium | fitted gap exponent | reference |
|---|---|---|
| real atrial surface (2-manifold) | N⁻¹·⁰⁵ | Weyl 2-manifold −1.0 (within 5%) |
| 2-D random-geometric | N⁻¹·⁸¹ | (faster than naive 2/d) |
| 3-D random-geometric | N⁻⁰·⁸⁴ | (faster than naive 2/d) |

So `ρ` grows lawfully with resolution — the biomarker fails *more*, predictably, as the medium is
meshed finer (the opposite of "more resolution helps"). Only the real manifold matches Weyl exactly;
we claim no more.

### 3.4 Localization transfers across three media — `results/gm4_*_metrics.json`

The Fiedler-gradient localizer `|∇φ₂|` locates the instability origin above a rotational spatial null
in all three media; centrality baselines (Perron, weighted degree, the `∩`-hotspot) are deleted by
measured correlation.

| medium | grad_φ₂ origin-rank | perm-p | verdict | fibrosis rank |
|---|---|---|---|---|
| cardiac / FHN (100k, converged) | 0.174 | <10⁻⁵ | KEEP | ~0.10 |
| Kuramoto oscillators | 0.091 | <10⁻⁴ | KEEP | 0.147 |

The pattern (predict-null + `|∇φ₂|`-KEEP) replicates across excitable **and** oscillator dynamics.
**Honest caveat:** the substrate field also localizes, and a controlled dissociation experiment
(`asb.transfer.bottleneck`) shows why they coincide — in a diffusively-coupled medium the Fiedler
bottleneck *is* the low-conduction lesion, which fibrosis imaging detects. `|∇φ₂|` is a purely
connectivity-derived localizer that beats graph-centrality baselines and, in the oscillator medium,
localizes more sharply than the raw substrate — but it does not find a site the substrate misses.
(Whether wavefront-curvature source–sink block can dissociate them on a healthy isthmus is under test.)

### 3.5 The falsification protocol — `results/falsification_protocol.json`, §11

Each guard, removed, manufactures a specific false positive: shape-family leakage inflates AUC by
**+0.117** (naive random-CV 0.756 vs grouped 0.640); the effect-size gate catches the **p-value trap**
(subspace-SFI p=2×10⁻²⁴ but ΔAUC +0.003, clinically useless). A naive analyst would have reported SFI
as a working biomarker; the pre-registered protocol correctly rejects it.

### 3.6 Uncertainty quantification — `results/e6_metrics.json`

Morris elementary-effects screening: the label is most sensitive to diffusion `d0`, then fibrosis
density. A GP surrogate has LOO R²=0.35 (genuine anatomical variance beyond the screened knobs). The
predictive null is robust across uncoupling strengths Δw ∈ {0.18, 0.36, 0.54}.

## 4. Discussion

**What is genuinely contributed.** (1) A pre-registered, leakage-controlled, effect-size-gated
falsification of a plausible spectral biomarker on two real cohorts and 10⁵ networks. (2) A validity
criterion `ρ` and a Weyl scaling law that explain the failure and predict it for the entire class of
spectral fragility biomarkers, testable a priori in any domain (seizure-focus localization, brain
connectomics, power-grid fragility): compute ρ; if ρ ≫ 1, the linear spectral marker cannot be
trusted regardless of in-sample correlation. (3) A connectivity-derived localizer that transfers
across three dynamical media. (4) A demonstration that standard methodology produces exactly the
false positives this protocol catches.

**What is not claimed.** The mathematics (Fiedler sensitivity, Weyl, Davis–Kahan, master stability
function) is classical; no new theorem is asserted. The localizer is not uniquely better than
substrate imaging. Labels are a monodomain simulator, not clinical outcomes and not openCARP.

**Limitations.** Simulator (not clinical) labels; monodomain, not the deferred openCARP gold standard;
the UW cohort's universal atrial coordinates are a PCA surrogate; UW inducibility is low, so the
UW-only predictive test is underpowered; ρ* and the random-geometric gap exponents are empirical.

## 5. Conclusion

A closed-form spectral fragility index does not predict atrial reentry inducibility beyond existing
features, on real anatomy and at scale — and we show, via a validity radius and a Weyl scaling law,
why no first-order spectral perturbation biomarker can on real excitable media. The honest core is
not a discovery but a boundary: where this class of biomarker is valid, and why real tissue lies far
outside it.

## References (all classical / established; verify formatting before submission)

Fiedler (1973) *Czech. Math. J.*; Ghosh & Boyd (2006) *IEEE CDC*; Pecora & Carroll (1998) *PRL*;
Barahona & Pecora (2002) *PRL*; Bomela, Wang, et al. (2020) *Sci. Rep.*; Weyl (1912); Davis & Kahan
(1970) *SIAM J. Numer. Anal.*; DeLong, DeLong & Clarke-Pearson (1988) *Biometrics* (fast form: Sun &
Xu 2014 *IEEE SPL*); Mitchell & Schaeffer (2003) *Bull. Math. Biol.*; Roney et al. LA virtual cohort
(Zenodo 5801337); Bifulco, Boyle et al. (2025) UW atrial meshes (Dryad 10.5061/dryad.kkwh70sg0).
