# AtrialSpectralBench: a pre-registered test of a spectral fragility biomarker for atrial reentry — a feature-redundancy null, an estimator-validity boundary, and the false positives naive methodology would have produced

> ## ⚠ SUPERSEDED — read `manuscript.tex` instead
>
> **`docs/paper/manuscript.tex` is the canonical paper.** This Markdown draft is kept for
> history only. It predates the substrate audit of 2026-07-24/25 and several of its headline
> numbers were computed on a defective substrate. Specifically, in the sections below:
>
> - the real-cohort rows use **Roney n=62** (it is now n=100, all released meshes) and a
>   **UW substrate whose atrial orifices were sealed** by a mis-read element tag;
> - the **combined ΔAUC = +0.051 (p = 0.155)** that "exceeds the 0.05 threshold" has been
>   recomputed to **+0.012 (p = 0.449)** on n=182 and no longer approaches the gate;
> - the secondary endpoint **ΔAUC = +0.103, p = 0.009** ("not merely re-encoded fibrosis")
>   has been **withdrawn**: it is +0.017 (p = 0.452) with the full cohort;
> - the framing "underpowered and inconclusive, not a clean null" is therefore obsolete;
>   the real-anatomy result is now a measured null.
>
> Corrected values, the audit that found the defect, and the three controls that failed to
> explain the residual are in `manuscript.tex` and `docs/PRE_REGISTRATION.md` §8.6. Individual
> superseded numbers below are flagged inline; nothing here should be cited.

*Working preprint draft (SUPERSEDED, see above). This is an in-silico computational study;
labels are electrophysiology-simulator verdicts, never clinical outcomes.*

## Abstract

We ask whether a closed-form spectral graph biomarker — the **Spectral Fragility Index (SFI)**,
the Fiedler-value sensitivity `∂λ₂/∂w_ij = (φ_i − φ_j)²` aggregated under a fibrosis-weighted
conduction-uncoupling field — predicts reentry inducibility in atrial tissue beyond established
fibrosis and connectivity features. Pre-registered before any labels were seen, the primary
endpoint was ΔAUC ≥ 0.05 with DeLong p < 0.05 on shape-family-held-out folds. **The endpoint is not
met — but the strength of the negative differs by regime, and we are careful to separate them.** On
up to 100,000 synthetic excitable networks SFI adds no practically meaningful predictive value
(subspace ΔAUC ≈ +0.003; ~4,150 cases needed at 80% power — statistically real only because N is
huge, so *clinically undetectable* in that cohort). On two real cohorts (182 patients: Roney
n=100, UW/Boyle n=82) the combined gradient-boosting estimate is **+0.012** with a bootstrap
interval [−0.017, +0.044] that excludes the threshold, so the real-anatomy result is a measured
null rather than a shortfall of power. (An earlier draft reported +0.051 here and called the
real-cohort evidence inconclusive; that figure was an artefact of a defect in one cohort's
released substrate — see the banner at the top of this file.) Where the null *is* clean (the synthetic cohort), the mechanism is **feature redundancy**:
the exact per-region Δλ₂ SFI — which has no first-order approximation error — also adds nothing over
the competitor set, so SFI's information is already carried by standard connectivity features. A
separate, label-free result explains why the *single-vector* SFI is additionally ill-conditioned on
real tissue: it is valid only while `ρ = ‖ΔL‖/(λ₃−λ₂) = O(1)`, and real atria have ρ ≈ 2422 (~10³×
past the empirical ρ*≈3), so by Davis–Kahan the Fiedler vector is unstable. We do **not** claim ρ
proves non-predictiveness — ρ bounds the accuracy of the Δλ₂ *magnitude*, not the classification
content of a ranking feature. The spectral gap collapses ∝ N⁻¹·⁰⁵ on the real atrial surface
(matching the 2-manifold Weyl exponent), so ρ grows lawfully with resolution. One component partially
survives: the Fiedler-gradient localizer `|∇φ₂|` — while the *strict* pre-registered localization
endpoint is null for all fields — is a weak but real rank localizer (origin in the top ~18%) that
beats graph-centrality baselines and transfers across three dynamical media, though it coincides with
(does not beat) substrate imaging, by a mechanism we make precise. Finally, naive methodology (random
CV, p-value thresholds, no effect-size gate) would have reported SFI as working; our pre-registered
protocol rejects it. The contribution is a rigorously-scoped negative with its mechanism (feature
redundancy), a label-free validity criterion for the single-vector estimator, a scaling
characterization, and a worked demonstration of the false positives standard methodology produces.
**A key limitation surfaced by our own convergence study: the labeller over-calls inducibility at the
operating mesh resolution and, on the completed 6-tier study (1500–48000 nodes, n=24), does not
converge within the tested window — the inducibility rate falls sharply from the coarse mesh (58% at
1500) but stays resolution-sensitive and even rises again at the finest affordable tier (37.5% at
~50k nodes). So all absolute rates and AUCs are coarse-mesh-provisional pending the deferred openCARP
validation; the relative SFI-vs-competitor comparison (same labels both arms) and the label-free
validity argument are unaffected.**

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
(self-sustained reentry ≥ 650 ms after pacing), run on each real anatomy. These are simulator
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

| cohort / medium | N | SFI ΔAUC vs full competitors (lr / gbt) | verdict |
|---|---|---|---|
| Roney (real) | 100 | −0.021 / +0.009 (p 0.444 / 0.625) | endpoint not met (34/100 inducible) |
| UW/Boyle (real) | 82 | +0.000 / +0.035 (p 1.000 / 0.669) | uninformative (8 positives, 9.8%) |
| Combined real | 182 | −0.007 / **+0.012** (p 0.704 / 0.449) | endpoint not met; no longer near the gate |

> *Superseded rows, kept for the record:* Roney 62, −0.030/+0.036 (p 0.64/0.26);
> UW/Boyle 82, −0.057/+0.104 (p 0.29/0.66) on 6 positives; Combined 144, −0.004/**+0.051**
> (p 0.92/0.155). These were computed with the UW atrial orifices sealed by a mis-read
> element tag and with only 62 of the 100 Roney meshes on disk.
| FHN networks | 2,000 | ≤ +0.003 (p ≥ 0.12) | clean null |
| FHN networks | up to 100,000 | subspace +0.0028 (p<10⁻⁴ at N≥10⁴) | clean null (effect negligible) |
| Kuramoto networks | 500 | ≤ +0.0036 (p ≥ 0.38) | clean null |

**The negative is not uniform, and we separate the regimes.** On the *synthetic* cohorts
(up to 100,000 networks) the null is clean: the effect converges to ΔAUC ≈ +0.003, statistically
resolvable only because N is enormous. On the *real* cohorts the combined gradient-boosting
estimate is **+0.012** (DeLong p=0.449, bootstrap CI [−0.017, +0.044]), roughly a quarter of the
pre-registered 0.05 threshold and with an interval that excludes it. With 42 inducible of 182 this
is the largest real-anatomy arm in the study, and the null here is a measurement rather than a
shortfall of power. The UW cohort remains uninformative rather than null (8 positives), so it
functions as an underpowered external check of the *labeller*, not of SFI.

> *Superseded framing, kept for the record:* this paragraph previously read "underpowered and
> inconclusive, not a demonstrated null", on the basis of a combined estimate of **+0.051**
> (p=0.155) with 26 inducible of 144. That estimate was an artefact of the sealed-orifice UW
> substrate; correcting it gives +0.036 at the original cohort size and +0.012 with all 100
> Roney meshes. The claim that real anatomy "cannot exclude a moderate effect" no longer holds.
> Note also that the phrase "competitor base AUC ≈ chance" for UW was simply wrong: the measured
> competitor base AUC is 0.834 (lr) / 0.748 (gbt). Only the *fibrosis-only* baseline is 0.562.

**Where the clean (synthetic) null comes from — feature redundancy, not the validity radius.** The
**exact** per-region Δλ₂ SFI — which carries *no* first-order approximation error and is immune to the
ρ argument of §3.2 — *also* adds nothing over the competitor set. So the mechanism of the null is that
SFI's information is already present in standard connectivity features (λ₂-alone, min-cut, percolation),
a label-dependent finding, not a consequence of ρ.

**Secondary endpoint (pre-specified) — WITHDRAWN.** This section previously reported that SFI
beats a *fibrosis-heterogeneity-only* baseline (combined GBT ΔAUC **+0.103, p=0.009**) and read
that as showing SFI is **not merely re-encoded fibrosis**. That result does not survive and the
claim is withdrawn. Correcting the UW substrate at the original cohort size leaves it largely
intact (+0.093, p=0.013), so it was not principally a substrate artefact; extending the Roney arm
from 62 to all 100 released meshes collapses it to **+0.017 (p=0.452)**. An effect that disappears
when 38 further real subjects are added was a small-sample fluctuation, and it was one of ~30
DeLong comparisons besides. The discarded estimate was carried by the UW arm, which under the
sealed substrate returned +0.340 (p=0.026) on six positives — a value that should not have been
believed at the time. On the full cohort no comparison, primary or secondary, shows SFI adding
information beyond either baseline.

**Power analysis (`results/power_analysis.json`).** In the *synthetic* cohort the subspace effect
(+0.0028) needs ~4,150 cases for 80% power and the single-vector effect (+0.0007) ~23,000 — clinically
undetectable *there*. We do **not** transport "clinically undetectable" to the real cohorts, which
are a separate regime. On the corrected substrate their point estimates are +0.009 (Roney, n=100)
and +0.012 (combined, n=182), both with intervals that contain zero and exclude the 0.05 gate.
(An earlier draft cited +0.036–0.051 here as "compatible with a detectable effect that we are
simply underpowered to confirm"; those estimates are superseded and that reading no longer holds.)

### 3.2 The validity radius (why it fails) — `docs/SFI_THEORY.md` §4, §9

The first-order SFI's relative error scales as `Θ(ρ)`, `ρ = ‖ΔL‖/(λ₃−λ₂)` (second-order
Rayleigh–Schrödinger; the Fiedler vector rotates by `≤ √2·ρ`, Davis–Kahan). Empirically the error
crosses 10% at ρ* ≈ 3. On real excitable media the gap is tiny, so **ρ is enormous — atrial median
ρ ≈ 2422**, ~10³× past ρ*. **What this proves, and what it does not.** It proves the *single-vector*
first-order SFI is a numerically **invalid estimator of Δλ₂** on real tissue and that `φ₂` is
ill-conditioned (Davis–Kahan) — a label-free statement about the estimator. It does **not** prove
SFI is non-predictive: ρ bounds the accuracy of the Δλ₂ *magnitude*, not the classification content
of a monotone ranking feature, and indeed the **exact** Δλ₂ SFI (no ρ limitation) is *also*
non-predictive (§3.1) — so the predictive null is separately due to feature redundancy. Two further
honest caveats: (i) ρ's numerator `‖ΔL‖` carries the fibrosis-weighted Δw field, so ρ is *not* fully
label-free — only its denominator `λ₃−λ₂` is; (ii) ρ*≈3 is calibrated on a 16-node path graph with a
problem-dependent constant, so "10³× past ρ*" compares an atrial value to a path-graph threshold and
should be read as order-of-magnitude, not exact. (Weyl, Davis–Kahan, second-order perturbation theory
are classical; the contribution is using ρ as an a-priori *estimator-validity* screen and measuring
its violation on real tissue.)

### 3.3 The ρ scaling law — `results/rho_scaling.json`, `docs/SFI_THEORY.md` §10

The gap collapses with system size by Weyl's law:

| medium | fitted gap exponent | reference |
|---|---|---|
| real atrial surface (2-manifold) | N⁻¹·⁰⁵ | Weyl 2-manifold −1.0 (within 5%) |
| 2-D random-geometric | N⁻¹·⁸¹ | (faster than naive 2/d) |
| 3-D random-geometric | N⁻⁰·⁸⁴ | (faster than naive 2/d) |

So `ρ` grows lawfully with resolution — the estimator degrades *more*, predictably, as the medium is
meshed finer (the opposite of "more resolution helps"). **Honest scope:** each exponent is a 5-point
fit; the 3-D random-geometric series is non-monotone (the fit is a loose upper bound, not a law);
these gaps are computed on the *unweighted* combinatorial Laplacian, whereas ρ uses the weighted
conduction Laplacian. **Only the real atrial 2-manifold cleanly matches Weyl (−1.05 vs −1.0)**; the
random-geometric exponents are reported as "at least as fast as Weyl," not law confirmations.

### 3.4 Localization transfers across three media — `results/gm4_*_metrics.json`

**The strict pre-registered localization endpoint is null for every field, including `|∇φ₂|`**
(median argmax-geodesic error is *not* below the spatial-null 5th percentile: 0.472 vs 0.406). What
survives is a weaker statement: `|∇φ₂|` is a **rank** localizer — the true origin sits in its top
~18% — that beats graph-centrality baselines (Perron, weighted degree, `∩`-hotspot, deleted by
measured correlation) and holds across media.

| medium | grad_φ₂ origin-rank | perm-p | fibrosis rank | n origins |
|---|---|---|---|---|
| FHN networks (100k) | 0.176 | <10⁻⁵ | ~0.10 | 36,540 |
| Kuramoto oscillators (500) | 0.091 | <10⁻⁴ | 0.147 | 243 |
| cardiac monodomain | ~0.28 | <10⁻³ | (co-localizes) | **only ~20** |

The weak-rank pattern replicates across excitable **and** oscillator dynamics, but three honesties
apply: (i) the strict endpoint is null, so `|∇φ₂|` is **not a precise localizer**; (ii) the substrate
field co-localizes and a controlled dissociation experiment (`asb.transfer.bottleneck`) shows *why* —
in a diffusively-coupled medium the Fiedler bottleneck *is* the low-conduction lesion — so `|∇φ₂|`
matches, and does **not** beat, substrate imaging (the transfer is partly a near-tautology); (iii) the
cardiac arm rests on only ~20 origins. The honest claim is a purely connectivity-derived *rank*
localizer that equals substrate localization and beats centrality baselines — nothing stronger.

### 3.5 The falsification protocol — `results/falsification_protocol.json`, §11

Each guard, removed, manufactures a specific false positive: shape-family leakage inflates AUC by
**+0.117** (naive random-CV 0.756 vs grouped 0.640); the effect-size gate catches the **p-value trap**
(subspace-SFI p=2×10⁻²⁴ but ΔAUC +0.003, clinically useless). A naive analyst would have reported SFI
as a working biomarker; the pre-registered protocol correctly rejects it.

### 3.6 Uncertainty quantification — `results/e6_metrics.json`

Morris elementary-effects screening: the label is most sensitive to diffusion `d0`, then fibrosis
density. A GP surrogate has LOO R²=0.35 (large residual anatomical/label variance beyond the screened
knobs — a label-noise caveat, not a positive). The synthetic null is robust across uncoupling
strengths Δw ∈ {0.18, 0.36, 0.54}.

### 3.7 External inducibility-rate anchor (base-rate plausibility, not clinical validation)

Clinically observed AF inducibility is not a single number — it spans ~5% to ~85% depending on
population, protocol aggressiveness, and the sustained-AF definition (Kumar et al. 2012, *Circ
Arrhythm Electrophysiol*, PMID 22528047, show burst-sustained 14.8% vs decremental-sustained 41.2%
*in the same patients*). Against the **verified** literature, our labeller's **Roney ~20–32%** sits
squarely in the moderate-burst band: Marquardt et al. 2018 (*J Atr Fibrillation*, PMID 30455834)
report 30.6% in patients without prior AF; Kumar 29.5% (sustained); between the no-AF-control floor
(Oral et al. 2008, *JCE*, PMID 18266669: 5% in controls) and the AF/aggressive ceiling (Kawai et al.
2019, *J Arrhythm*, PMID 31007786: 51% persistent AF; Oral 84% with isoproterenol), consistent with
~42% post-ablation (Liu et al. 2020, *JAHA*, PMID 32654581). The **UW ~7%** is at/below the control
floor for a cohort that *is* AF patients — flagged as anomalous and attributable to the substrate
rather than to biology. ("Rule-based fibres" in an earlier draft was a euphemism: the released
meshes carry a *constant* fibre vector, not a rule-based field. Repairing that and the sealed
orifices lifts the rate to 17.1%, still short of the 32.3% seen on fibre-carrying anatomy, so the
anomaly is only partly accounted for and the in-silico inducer may also be under-aggressive.) These are in-silico verdicts with a simulator-scaled ≥650 ms cutoff, so this is a base-rate
**plausibility check, not clinical validation.** (One candidate source, Darma et al. 2020, is *not*
cited: its adversarial citation check did not verify.)

## 4. Discussion

**What is genuinely contributed.** (1) A pre-registered, leakage-controlled, effect-size-gated
evaluation of a plausible spectral biomarker: a **clean null at synthetic scale** (10⁵ networks) with
its mechanism (**feature redundancy** — even the exact Δλ₂ SFI adds nothing), and a **measured
null** on 182 patients of real anatomy (combined +0.012, CI [−0.017, +0.044]).
(2) A **label-free validity criterion** `ρ` for the single-vector estimator — it proves `φ₂` is
ill-conditioned on real tissue (ρ ≈ 2422 ≫ ρ*≈3, Davis–Kahan), a portable a-priori check for any
Fiedler-based marker — with the honest scope that ρ bounds estimator accuracy, *not* predictive
content, and that ρ's numerator is not itself label-free. (3) A weak-but-real, purely
connectivity-derived **rank** localizer that equals substrate localization and beats centrality
baselines across three media. (4) A demonstration that standard methodology manufactures the exact
false positives this protocol catches (leakage +0.117 AUC; the p-value trap).

**What is NOT claimed.** The mathematics is classical; no new theorem. ρ does **not** prove SFI is
non-predictive — only that the single-vector estimator is ill-conditioned; the null is separately due
to feature redundancy. The real-cohort result is **not** a demonstrated null — it is underpowered.
The localizer is **not** precise (strict endpoint null) and does **not** beat substrate imaging. The
ρ scaling "law" is a clean match to Weyl only on the real atrial 2-manifold. Labels are a monodomain
simulator, not clinical outcomes and not openCARP.

**Limitations.** (a) **The labeller over-calls inducibility at the operating mesh resolution, and does
not converge in the tested window** — our completed 6-tier convergence study (Roney, n=24 at each of
1500/3000/6000/12000/24000/48000 nodes) shows the inducibility rate drops sharply from the coarse mesh
(58% at 1500) to 17% at 3000, but is then non-monotone across the fine tiers (25%, 21%, 25%) and rises
to **37.5% at the finest affordable tier (~50k nodes)** — a 17-percentage-point spread over the last
three tiers, so it has *not* settled. Per-subject verdict agreement between adjacent tiers does climb
(50%→79% by 24k→48k) and only 7/24 subjects are fully consistent across all six tiers; the rate does
*not* collapse to zero (the phenomenon is real) but remains resolution-sensitive up to ~50k nodes. The
main results ran at 2000 nodes, so all *absolute* rates and AUCs are coarse-mesh-provisional pending
the deferred openCARP validation (the *relative* SFI-vs-competitor comparison scores the same labels
in both arms, so this bias cancels there). (b) Simulator, not clinical, labels
(the pre-registered openCARP ground truth was substituted with monodomain Mitchell–Schaeffer — a
logged deviation). (c) The UW cohort uses PCA-surrogate UAC; its released meshes carry a *degenerate*
fibre field (a constant vector in all 164) and an element tag we mis-read as tissue when it is the
caps over the atrial orifices. Repairing both raises inducibility only from 6/82 to 14/82, and all
three candidate explanations for its low rate are now excluded by direct control, so the cohort's
predictive test is uninformative and the anomaly is an open problem. (d) ~30 DeLong comparisons were
run; the one apparent positive (fibrosis_vs_+SFI, +0.103 at p=0.009) has since been **withdrawn** —
it is +0.017 (p=0.452) on the full cohort. (e) ρ* and the random-geometric gap exponents are
empirical. (f) The inducibility label is reproducible in aggregate but unstable per subject, which
attenuates every ΔAUC reported here.

## 5. Conclusion

A closed-form spectral fragility index does not add practically meaningful reentry-inducibility
prediction beyond existing features, either at synthetic scale or on 182 patients of real anatomy,
and the mechanism is feature redundancy. Independently, we show the single-vector estimator is
provably ill-conditioned on real tissue (ρ ≫ 1). The core is not a discovery but a carefully-scoped
boundary — where a Fiedler-based estimator is numerically valid, why real tissue lies far outside
it, and how easily naive methodology would have mistaken redundant signal for a working biomarker.
This project produced such a false positive itself (+0.103 at p=0.009) and its own protocol
withdrew it. The absolute rates await the openCARP validation, which is now feasible.

## References (all classical / established; verify formatting before submission)

Fiedler (1973) *Czech. Math. J.*; Ghosh & Boyd (2006) *IEEE CDC*; Pecora & Carroll (1998) *PRL*;
Barahona & Pecora (2002) *PRL*; Bomela, Wang, et al. (2020) *Sci. Rep.*; Weyl (1912); Davis & Kahan
(1970) *SIAM J. Numer. Anal.*; DeLong, DeLong & Clarke-Pearson (1988) *Biometrics* (fast form: Sun &
Xu 2014 *IEEE SPL*); Mitchell & Schaeffer (2003) *Bull. Math. Biol.*; Roney et al. LA virtual cohort
(Zenodo 5801337); Bifulco, Boyle et al. (2025) UW atrial meshes (Dryad 10.5061/dryad.kkwh70sg0).

*Clinical inducibility-rate anchors (adversarially verified, PMIDs confirmed):* Kumar et al. (2012)
*Circ Arrhythm Electrophysiol* 5(3):531 (PMID 22528047); Marquardt et al. (2018) *J Atr Fibrillation*
11(1):1837 (PMID 30455834); Oral et al. (2008) *J Cardiovasc Electrophysiol* 19(5):466 (PMID 18266669);
Kawai et al. (2019) *J Arrhythm* 35(2):223 (PMID 31007786); Liu et al. (2020) *JAHA* 9(14):e015260
(PMID 32654581). *Δw perioperative CV-slowing calibration:* Heida et al. (2021) *J Clin Med* (PMID
34198544; value as-cited, not independently re-confirmed).
