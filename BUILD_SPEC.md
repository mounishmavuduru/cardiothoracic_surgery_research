# BUILD_SPEC — AtrialSpectralBench engine (frozen interface contract)

This is the single source of truth for the AtrialSpectralBench Python package. Every
build agent implements exactly these signatures against the frozen types in
`src/asb/types.py` and config in `src/asb/config.py` (already written — do **not**
modify them). Import package name is `asb` (src layout, editable-installed).

## 0. Environment & rules (read first)

- Python 3.11, CPU-only, no GPU. Deps are installed in a venv at `.venv` via
  `pip install -e ".[dashboard,dev]"`. Core deps: numpy, scipy, networkx,
  scikit-learn, pandas, matplotlib, plotly, pyyaml, tqdm; dev: pytest.
- **No network at build time** (zenodo.org / opencarp.org / github are blocked by
  policy). So: no downloads, no openCARP execution. Real-data loaders and the
  openCARP runner are written as clean, well-documented interfaces that **raise an
  informative error** if the data/binary is absent — they are NOT executed here.
- Determinism: every stochastic function takes an explicit `rng: np.random.Generator`
  or `seed: int`. No global `np.random` state, no `Date.now`, no wall-clock seeds.
- `mock_ep` labels are a **development stand-in**, never clinical POAF. Every public
  surface that emits or consumes them must say so in its docstring.
- Style: type hints, numpy-style docstrings, small pure functions, no I/O in compute
  functions (I/O lives in loaders / pipeline / figures). Match the backbone files.
- Each module ships a `tests/test_<module>.py`. Tests must pass with `pytest -q`.
- **Smoke-test isolation (parallel phase):** when self-testing your module, depend
  only on `asb.types`, `asb.config`, numpy/scipy, and your own module. Do NOT import
  sibling `asb.*` modules other agents are writing concurrently; construct small
  local inputs instead. Cross-module integration is verified in the final phase.

## 1. Canonical types (already implemented — import, do not redefine)

- `asb.types.AtrialGraph` — coords(n,3), edges(m,2 int, i<j), weights(m,), fibrosis(n,),
  uac(n,2), region(n,), shape_family(str); `.n_nodes`, `.n_edges`, `.adjacency()->csr`,
  `.degree()->(n,)`.
- `asb.types.AtrialMesh` — points, faces, fibres, uac, fibrosis, region, shape_family.
- `asb.types.InducibilityLabel` — inducible(bool), reentry_origin(int|None),
  protocol(str), source('mock_ep'|'opencarp').
- `asb.config.Config` with `.load(path)`, `.to_yaml(path)`; sub-configs
  `SFIConfig, CohortConfig, LabelConfig, EvalConfig`.

## 2. Modules to implement

### 2.1 `src/asb/graph.py` — Laplacians & weighted-graph construction
```
def laplacian(G: AtrialGraph, kind: str = "combinatorial") -> scipy.sparse.csr_matrix
    # kind in {"combinatorial" (L=D-W), "sym" (I - D^-1/2 W D^-1/2), "rw" (I - D^-1 W)}.
def laplacian_from_adjacency(W, kind="combinatorial") -> csr
def edge_weights_from_fibres(mesh: AtrialMesh, edges, *, along=1.0, cross=0.3,
                             fibrosis_floor=0.05) -> np.ndarray
    # conductance high along fibre, low cross-fibre, multiplied by (1-fibrosis) with a floor.
def cut_edges(G: AtrialGraph, edge_ids) -> AtrialGraph
    # return a copy with those edge weights set ~0 (discrete surgical cut, e.g. Maze).
```
Contract: `laplacian` is symmetric PSD with zero row sums for the combinatorial kind.

### 2.2 `src/asb/spectral.py` — eigenpairs & spectral features
```
def smallest_eigpairs(L, k=6) -> (vals: (k,), vecs: (n,k))   # ascending, via eigsh(sigma=0) or shift-invert
def fiedler(L) -> (lambda2: float, phi2: (n,))               # 2nd-smallest eigenpair
def spectral_gap(L) -> float                                  # lambda3 - lambda2
def perron(W) -> (rho: float, v: (n,))                        # dominant eigenpair of nonnegative W; v >= 0
def inverse_participation_ratio(v) -> float                   # IPR = sum v^4 / (sum v^2)^2
def spectral_features(G: AtrialGraph) -> dict                 # gap, spectral_entropy, cheeger_estimate,
                                                              # n_near_zero, perron_ipr, heat_kernel_trace(t)
def cheeger_estimate(L) -> float                              # lambda2/2 <= h <= sqrt(2 lambda2) bounds midpoint
```
Contract: on the path graph P_N (N nodes), `smallest_eigpairs` reproduces
`lambda_k = 2 - 2 cos(k*pi/N)`, k=0..N-1. Perron vector is entrywise nonnegative.

### 2.3 `src/asb/sfi.py` — **THE PROTECTED NOVEL SEED**
```
def edge_fragility(phi2, edges) -> np.ndarray                 # (m,) = (phi2_i - phi2_j)^2  == dλ2/dw_ij
def perturbation_field(G, cfg: SFIConfig, rng) -> dict        # per-edge E[Δw_ij] and cov from fibrosis-weighted
                                                              # diffuse stress; Δw grows with local fibrosis.
def sfi_region(G, phi2, edges, region_ids, expected_dw) -> dict
    # per-region SFI(R) = sum_{(i,j) in R} E[Δw_ij] (phi2_i - phi2_j)^2  (first-order analytic expectation).
def sfi_monte_carlo(G, L, cfg, rng) -> dict                   # sample Δw fields, recompute Δλ2 exactly per draw,
                                                              # return per-region mean/std of Δλ2 (MC ground truth for the analytic expectation).
def validity_radius_ok(dL_norm: float, lambda2: float, lambda3: float, safety: float) -> bool
    # True iff dL_norm <= safety*(lambda3 - lambda2)  (first-order single-vector regime valid).
def subspace_sfi(L, edges, expected_dw, k_dim=2) -> np.ndarray
    # projector/Davis-Kahan generalization using the invariant subspace spanning {phi2,...}; well-defined when λ2 near-degenerate.
def exact_delta_lambda2(L, edge_id, dw) -> float              # recompute λ2 exactly after editing one edge (honest boundary for large edits).
def hotspot_map(G, phi2, perron_v) -> np.ndarray             # per-node score = |grad phi2| (edge-avg) * perron localization; predicts reentry initiation.
```
Contracts (hard tests):
- **Derivative identity:** for a small edge perturbation ε on a simple λ2,
  `(exact_delta_lambda2(L, e, ε) - 0)/ε ≈ edge_fragility(phi2, edges)[e]` to O(ε).
  Verify via finite difference on a random weighted graph with a clear spectral gap.
- Analytic sfi_region expectation agrees with sfi_monte_carlo mean within MC error
  in the small-Δw regime.
- `validity_radius_ok` flips to False and code path switches to `subspace_sfi` when
  the λ2–λ3 gap is made small (near-degenerate construction).

### 2.4 `src/asb/substrate/mesh.py` — mesh <-> graph
```
def mesh_to_graph(mesh: AtrialMesh, *, along=1.0, cross=0.3) -> AtrialGraph
    # edges = unique triangle edges; weights from edge_weights_from_fibres; carries
    # fibrosis/uac/region/shape_family through.
def mesh_edges(faces) -> np.ndarray                          # (m,2) unique undirected edges, i<j.
def geodesic_edge_length(points, edges) -> np.ndarray
```

### 2.5 `src/asb/substrate/synthetic.py` — synthetic atrial anatomy (no patient data)
```
def make_base_atrium(seed: int, n_nodes=1200, family: str = None) -> AtrialMesh
    # closed LA-like surface (ellipsoidal body + a few pulmonary-vein stubs), a smooth
    # rule-based fibre field, UAC-like (alpha,beta) coords, region labels. Deterministic in seed.
def shape_variant(base: AtrialMesh, seed: int) -> AtrialMesh  # perturb shape modes; SAME shape_family as base.
def paint_fibrosis(mesh: AtrialMesh, seed: int, burden=0.2, n_patches=5, patch_scale=0.15) -> AtrialMesh
    # patchy fibrosis field in [0,1]; returns a copy with .fibrosis set.
def make_cohort(cfg: CohortConfig) -> list[AtrialMesh]        # n_base families x n_variants, each fibrosis-painted, family-tagged.
```
Contract: every mesh is a valid closed-ish surface (all vertices used by >=1 face);
`shape_variant` preserves `shape_family`; fibrosis in [0,1].

### 2.6 `src/asb/substrate/loaders.py` — real open-data loaders (interface only)
```
def load_rodero_cohort(root: str) -> list[AtrialMesh]         # Zenodo 4506930 (1000 SSM meshes)
def load_roney_meshes(root: str) -> list[AtrialMesh]          # Zenodo 5801337 (LA + fibrosis + UAC + fibres)
def load_fibre_atlas(root: str) -> list[AtrialMesh]           # Zenodo 3764917 (DT-MRI fibre atlas, CARP format)
```
Contract: if `root` is missing/empty, raise `DataUnavailableError` (define it here)
with the exact Zenodo record + download instructions in the message. Parse VTK/CARP
formats when present. These are NOT run in this environment (network blocked).

### 2.7 `src/asb/labels/mock_ep.py` — EP inducibility stand-in (DEVELOPMENT ONLY)
```
def induce(G: AtrialGraph, cfg: LabelConfig, rng) -> InducibilityLabel
    # Lightweight excitable-media surrogate on the graph: a monodomain-style
    # reaction-diffusion (Mitchell-Schaeffer / FitzHugh-Nagumo) OR eikonal + wavelength
    # (CV x ERP) reentry criterion, with S1-S2/burst pacing from a few sites. Returns a
    # binary inducible flag + reentry-origin node. source='mock_ep'.
def paced_activation(G, site, cfg) -> np.ndarray             # per-node activation times (helper).
```
Docstring MUST state: stand-in for openCARP; not clinical POAF; used only to exercise
the downstream pipeline until real labels exist. Should be deterministic in rng and
tend to label low-conductance / high-fibrosis, high-SFI atria as more inducible (so the
pipeline has learnable signal), WITHOUT trivially copying SFI (avoid circularity — base
it on wavelength/source-sink, not on λ2 directly).

### 2.8 `src/asb/labels/opencarp.py` — real solver runner (interface only)
```
def write_carp_mesh(mesh: AtrialMesh, out_dir: str) -> dict   # .pts/.elem/.lon
def make_pacing_protocol(cfg: LabelConfig) -> dict
def run_opencarp(mesh_dir: str, protocol: dict, *, opencarp_bin=None) -> InducibilityLabel
    # If opencarp not found -> raise OpenCARPUnavailableError with install pointer. Not run here.
def reproduce_niederer_benchmark(...) -> dict                 # code-verification target (deferred run).
```

### 2.9 `src/asb/baselines.py` — competitor features (kill the strawman)
```
def fibrosis_burden(G) -> float
def fibrosis_spatial_entropy(G, n_bins=16) -> float
def fibrosis_patch_size(G) -> float                          # mean connected fibrotic-patch size.
def min_cut_value(G) -> float                                # deterministic global min-cut (networkx / Stoer-Wagner).
def percolation_threshold(G, n_steps=25, rng=None) -> float  # edge-removal fraction at which giant component breaks.
def lambda2_alone(G) -> float
def baseline_feature_vector(G) -> dict                       # all of the above, named.
```

### 2.10 `src/asb/features.py` — assemble the design matrix
```
def spectral_feature_vector(G) -> dict                       # from spectral.spectral_features + fiedler stats.
def sfi_feature_vector(G, cfg: SFIConfig, rng) -> dict        # per-region SFI summary stats (max/mean/top-k).
def subject_features(G, cfg, rng) -> dict                     # merge baseline + spectral + SFI; namespaced keys.
def build_design_matrix(graphs, labels, cfg, rng) -> (pandas.DataFrame X, np.ndarray y, groups: list[str])
    # groups = shape_family per subject (for GroupKFold). y = inducible flags.
```

### 2.11 `src/asb/evaluation.py` — leakage-controlled stats
```
def nested_group_kfold_auc(X, y, groups, feature_sets: dict, cfg: EvalConfig, seed) -> dict
    # feature_sets maps name -> column list; e.g. {"fibrosis": [...], "fibrosis+SFI": [...]}.
    # returns grouped AND naive AUC per feature set (expose PCA-resample leakage).
def delong_test(y_true, prob_a, prob_b) -> dict              # {"auc_a","auc_b","delta","z","p"} paired DeLong.
def roc_auc_ci(y_true, prob, n_boot, seed) -> dict
def colocalization_vs_null(hotspot_scores, origin_node, uac, n_null, seed) -> dict
    # geodesic/UAC distance from predicted hotspot to true reentry origin vs a rotational/
    # shift UAC spatial-null; returns observed distance, null 5th pct, and pass flag.
def calibration_curve_data(y_true, prob, n_bins=10) -> dict
```
Contract: `delong_test` matches a known reference case; grouped AUC <= naive AUC on a
constructed leaky dataset; colocalization returns a valid p/percentile in [0,1].

### 2.12 `src/asb/figures.py` — matplotlib/plotly figures
```
def fig_fiedler_atrium(mesh, phi2, out) ; fig_sfi_vs_origin(mesh, sfi, origin, out)
def fig_eigen_spectrum(vals, out) ; fig_roc_panel(results, out) ; fig_sensitivity_bars(sobol, out)
```
Each writes a PNG (matplotlib) and returns the path; pure rendering, no compute.

### 2.13 `src/asb/pipeline.py` + `src/asb/cli.py` — end-to-end (config-driven)
```
def run(cfg: Config) -> dict
    # cohort -> mesh_to_graph -> spectral + SFI -> mock_ep labels -> build_design_matrix ->
    # nested_group_kfold_auc (fibrosis vs fibrosis+SFI) + DeLong + colocalization -> figures.
    # writes outputs/results_report.md + PNGs + metrics.json. Fully runnable on synthetic data.
# cli.main(): `asb run --config configs/default.yaml`  (also `python -m asb.pipeline`).
```

### 2.14 `src/asb/dashboard/app.py` — Streamlit
Interactive: pick a synthetic atrium, show φ2 coloring, SFI hotspot map, per-subject
inducibility + risk readout. Import-guarded so the core package works without streamlit.

## 3. Tests (tests/) — analytic gates are HARD CI gates
- `test_spectral.py`: path/ring/grid Laplacian eigenvalues vs closed form
  (path/ring: `2-2cos(k*pi/N)` / `2-2cos(2*pi*k/N)`; grid = pairwise sums); Perron nonneg.
- `test_graph.py`: Laplacian PSD, zero row-sum (combinatorial), sym-normalized spectrum in [0,2].
- `test_sfi.py`: derivative identity via finite difference; analytic vs MC agreement;
  validity-radius switch to subspace SFI on near-degenerate λ2; exact recompute for a cut.
- `test_baselines.py`, `test_substrate.py`, `test_labels.py`, `test_evaluation.py`: shape/þrange/determinism.
- `test_pipeline.py`: end-to-end on a tiny synthetic cohort produces metrics.json + a report.
- Property tests seeded and fast (<~60s total).

## 4. Definition of done
`pip install -e ".[dashboard,dev]"` then `pytest -q` all green; `asb run --config
configs/default.yaml` produces `outputs/results_report.md`, `outputs/metrics.json`, and
the five figures on synthetic data; dashboard imports without error.
