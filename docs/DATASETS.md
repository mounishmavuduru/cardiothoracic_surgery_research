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
