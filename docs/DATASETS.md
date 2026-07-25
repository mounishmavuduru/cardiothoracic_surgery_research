# Public atrial-anatomy datasets for cohort expansion (verified)

Candidate open datasets of **distinct real (or realistic) left-atrial anatomies** to grow
the cohort beyond the current Roney 100. Verified via a literature/repository sweep.
**Distinctness:** Real (independent patients) > SSM (statistical-shape-model draws) >
Variant (perturbations of one shape). Grouped CV weights *distinct* anatomies, so prefer
Real over SSM/Variant (100 000 variants of one shape ≈ 1 effective group — see ROADMAP §5).

| # | Dataset | Content (N distinct / format / fibrosis / fibres / UAC) | Access | License | Distinctness |
|---|---|---|---|---|---|
| — | **Roney LA virtual cohort** (in use) | 100 patient LA; `.vtk`; IIR-fibrosis ✓; fibres ✓; UAC ✓ | Zenodo 5801337 | CC-BY-4.0 | Real (100) |
| 1★ | **UW/Boyle LGE-MRI LA meshes** | **82 distinct patients** ×2 (pre/post-ablation) = 164 `.vtk`; LGE fibrosis/scar ✓; fibres ✗; UAC ✗ | Dryad `10.5061/dryad.kkwh70sg0` | CC0 (default) | **Real (82)** |
| 2 | **Nagel/Karlsruhe bi-atrial SSM** | ~195 bi-atrial volumetric `.vtk`; fibres ✓; wall thickness/bridges ✓; fibrosis ✗ (add rule-based) | Zenodo 4309957 → 5571925 | CC-BY-4.0 | SSM |
| 3 | **LAScarQS 2022** (MICCAI) | 194 LGE-MRIs (LA + scar masks, NIfTI; **not meshes**); no fibres/UAC | zmiclab.github.io/projects/lascarqs22 (registration) | Challenge terms | Real (needs meshing) |
| 4 | **2018 LA Segmentation Challenge** | 154 LGE-MRIs (100 train + 54 test); `.nrrd` images + LA masks (**not meshes**) | cardiacatlas.org/atriaseg2018 (register + cite) | Registration | Real (needs meshing) — **anatomy-realism check (E7)** |
| 5 | **Roney Human Atrial Fibre Atlas** | 7 ex-vivo; `.vtk` + CARP; DT-MRI fibres ✓ (gold std); UAC ✓; fibrosis ✗ | Zenodo 3764917 | CC-BY-4.0 | Real (7; reference fibres) |
| 6 | **CEMRG 4-chamber (sex/disease)** | 50 patient 4-chamber (26 healthy/24 HF); tet + fibres; LA extractable | Zenodo 19351401 | CC-BY-4.0 | Real (50) |
| 7 | **Rodero 24 whole-heart** | 24 patient 4-chamber; tet + rule-based fibres; LA extractable | Zenodo 3890033 | CC-BY-4.0 | Real (24) |
| 8 | **MMWHS** | 120 whole-heart (~40 with LA labels); images+masks (**not meshes**) | zmiclab.github.io/zxh/0/mmwhs (register) | Registration | Real (needs meshing) |
| 9 | **Rodero 1000 synthetic 4-chamber** | 1000 `.vtk` (~27.5 GB); LA extractable; no fibres/UVC in batch | Zenodo 4506930 | CC-BY-4.0 | **Variant** (PCA ±2σ; least independent) |

## Recommended additions (in order)

1. **UW/Boyle Dryad (`10.5061/dryad.kkwh70sg0`) — add first.** ~82 additional *distinct real
   patient* LA anatomies, already `.vtk` with LGE fibrosis → drops into the Roney loader with
   minimal work. Missing UAC + fibres, so either (a) run through atrialmtk/UAC to add them, or
   (b) use rule-based fibres + a coarse UAC surrogate for the graph. **Nearly triples distinct
   real N (100 → ~182).** Biggest, cheapest credibility win for generalization.
2. **Nagel bi-atrial SSM (Zenodo 5571925)** for simulation-ready volumetric meshes with fibres
   (SSM, so grouped as its own shape families).
3. **2018 LA Segmentation Challenge** — use as the **anatomy-realism check (E7)**, not as extra
   training cohort (it's masks, needs meshing).

## Integration notes
- The Roney loader (`asb.substrate.roney`) already parses `.vtk` POLYDATA + point arrays; the
  Dryad meshes are `.vtk` but lack the `UAC1/UAC2/IIR/fiber_*` arrays, so a small adapter is
  needed (LGE→fibrosis map; rule-based fibres via a Laplace–Dirichlet field; UAC surrogate).
- Each new *dataset* is a natural extra GroupKFold stratum; keep patient = group.
- Deferred to Claude Science where the full atrialmtk→UAC→fibre→openCARP pipeline runs at scale.

## Real clinical-outcome labels (recurrence / POAF) — availability audit (2026-07-23)

> **UPDATE 2026-07-24 — superseded for the UW/Boyle row.** Prof. Patrick M. Boyle responded to
> the outreach request and shared the per-patient 2-year recurrence outcomes for all 82 UW
> patients (`NR` / `AF` / `AFL`, keyed to the public Dryad mesh IDs). The verdict below —
> "no public dataset ties usable atrial geometry to a real clinical outcome" — remains true
> of *public downloads*, but is no longer true of this project: the author-request path (1)
> worked. The labels are restricted-use and are **not** redistributable; see
> `docs/PRE_REGISTRATION.md` §8 for the custody terms and the frozen analysis plan.

**Verdict: no public dataset ties usable atrial geometry to a real clinical outcome label in a
single open deposit.** Mesh deposits are released *because* outcome labels were stripped for
confidentiality. Implications for a "predict real recurrence" pivot: not feasible off-the-shelf.

| dataset | geometry | real outcome in public download? | note |
|---|---|---|---|
| UW/Boyle Dryad `10.5061/dryad.kkwh70sg0` | 82 pts, LGE meshes ✓ | **NO in public download** — recurrence withheld ("protecting confidential patient information"); code repo ships an empty template. **Obtained 2026-07-24 by author request** (restricted use, not redistributable) | paper Comms Med 2025 `10.1038/s43856-025-01058-4`, model AUROC 0.80 |
| Roney/CEMRG 100 LA models (cemrg.com) | 100 pts, meshes+fibres+UAC ✓ | **NO in download**, but outcome EXISTS in study (34/99 recurred, 1 yr) | best target for an author **request**; paper `10.1161/CIRCEP.121.010253` (PMC8845531) |
| UK Biobank cardiac MRI | LA from cine CMR (no LGE/fibrosis) | **YES** — incident AF/stroke via linked records | application-gated + fee; different question (new-onset, not recurrence/POAF) |
| AtriaSeg 2018, LAScarQS 2022 | LGE + masks/scar ✓ | **NO** — segmentation benchmarks only | — |
| Zenodo SSM/atlas meshes (4309957, 5801337, 3890033, …) | meshes ✓ | **NO** — shape-model/atlas instances | — |

**Only realistic paths to real outcomes + atrial geometry:** (1) email the Roney/CEMRG or Boyle
groups requesting recurrence labels keyed to their already-public meshes (highest value; geometry is
ML-ready and they demonstrably hold the labels); (2) UK Biobank application (real outcomes at scale,
but gated, no fibrosis, incident-AF not recurrence); (3) a formal data-use agreement / clinical
collaboration. None is a same-day download.
