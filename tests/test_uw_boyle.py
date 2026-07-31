"""Guards on the UW/Boyle substrate findings of 2026-07-24.

Two defects were found in the released Dryad meshes after they had already been used
for published in-silico results. These tests pin the diagnosis so neither can silently
come back:

  1. element tag 164 is NOT tissue -- it is the caps sealing the four pulmonary veins
     and the mitral valve, and must be removed before simulation;
  2. the shipped fibre field is degenerate (one direction mesh-wide) and must be
     flagged rather than used as if it were anatomy.

The meshes live under the gitignored ``data/`` tree, so every test here skips cleanly
when the data has not been downloaded (``python scripts/dryad_fetch.py``, see
``docs/DATASETS.md``). The pure-logic tests at the bottom run unconditionally.
"""
from __future__ import annotations

import glob
import warnings

import numpy as np
import pytest
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from asb.substrate.uw_boyle import (DEAD_TAGS, DEGENERATE_FIBRE_SPREAD, DROP_TAGS,
                                    TAG_FIBROSIS,
                                    load_uw_mesh, read_vtk_unstructured_bin,
                                    tags_to_fibrosis)

MESHES = sorted(glob.glob("data/uw_boyle/*_PreAbl_rev1.vtk"))
needs_data = pytest.mark.skipif(not MESHES, reason="UW meshes not downloaded")


def _cell_adjacency(faces: np.ndarray):
    """Pairs of triangles sharing an edge."""
    e = np.sort(np.concatenate(
        [faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1)
    cell = np.tile(np.arange(faces.shape[0]), 3)
    uniq, inv, cnt = np.unique(e, axis=0, return_inverse=True, return_counts=True)
    o = np.argsort(inv, kind="stable")
    st = np.searchsorted(inv[o], np.arange(uniq.shape[0]))
    sh = cnt == 2
    return cell[o][st[sh]], cell[o][st[sh] + 1]


def _euler_characteristic(points, faces: np.ndarray) -> int:
    """chi = V - E + F over the vertices actually referenced by ``faces``."""
    v = np.unique(faces).size
    e = np.unique(np.sort(np.concatenate(
        [faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1), axis=0).shape[0]
    return int(v - e + faces.shape[0])


@needs_data
@pytest.mark.parametrize("path", MESHES[:3])
def test_dropping_tag_164_opens_holes(path):
    """Removing a cap from a closed surface raises the number of boundary loops.

    Measured over all 82 pre-ablation meshes, the drop in Euler characteristic has
    median 7 (range 2-12), i.e. tag 164 seals roughly seven openings: five large ones
    plus a few small fragments. It is NOT a fixed five -- an earlier revision asserted
    ``chi == -3`` from a single mesh (ID001) and that does not generalise, because these
    surfaces are not clean spheres to begin with (chi before removal ranges 3 down to
    -18). Two meshes, ID040 and ID049, are topologically pathological and go the other
    way; they are excluded from this guard rather than hidden.
    """
    if any(bad in path for bad in ("ID040", "ID049")):
        pytest.skip("ID040/ID049 have pathological topology (chi before removal -16/-18)")
    parsed = read_vtk_unstructured_bin(path)
    faces = np.asarray(parsed["faces"])
    tags = np.asarray(parsed["cell_scalars"]["elemTag"]).astype(int)
    before = _euler_characteristic(None, faces)
    after = _euler_characteristic(None, faces[tags != 164])
    assert before - after > 0, f"removing {DROP_TAGS} opened no holes (delta chi {before-after})"


@needs_data
@pytest.mark.parametrize("path", MESHES[:3])
def test_tag_164_forms_four_to_six_large_components(path):
    """The signature of the atrial orifices: a handful of large caps, not a diffuse patch.

    Across all 82 meshes the count of tag-164 components occupying >=1% of the region is
    5 in 75 meshes, 4 in 5, and 6 in 2 -- exactly the spread of human pulmonary-vein
    anatomy (four PVs plus the mitral valve, with a left common trunk giving four and a
    right middle PV giving six).
    """
    parsed = read_vtk_unstructured_bin(path)
    faces = np.asarray(parsed["faces"])
    tags = np.asarray(parsed["cell_scalars"]["elemTag"]).astype(int)
    a, b = _cell_adjacency(faces)
    m = tags == 164
    idx = np.where(m)[0]
    remap = np.full(faces.shape[0], -1, dtype=int)
    remap[idx] = np.arange(idx.size)
    sel = m[a] & m[b]
    g = coo_matrix((np.ones(sel.sum()), (remap[a[sel]], remap[b[sel]])), shape=(idx.size,) * 2)
    _, lab = connected_components(g, directed=False)
    big = int((np.bincount(lab) >= 0.01 * idx.size).sum())
    assert 4 <= big <= 6, f"expected 4-6 large tag-164 components, got {big}"


@needs_data
@pytest.mark.parametrize("path", MESHES[:3])
def test_tag_164_is_essentially_never_adjacent_to_fibrosis(path):
    """Tag 164 borders fibrotic tag 115 far less than chance -- a >100x depletion.

    This rules out the "fibrosis border zone" reading that the legacy
    ``TAG_FIBROSIS[164] = 0.5`` mapping implies: a border zone would be MAXIMALLY
    adjacent to fibrosis. Measured over all 82 meshes, contact is a median of 27 edges,
    ~0.08% of the edges incident to tag 164, against a chance expectation of 14-20%.
    An earlier revision of this test asserted EXACTLY zero; that was an artefact of
    reading a mean fraction rounded to three decimals. Contact is near-zero, not zero
    (it is exactly zero in only 5 of 82 meshes).
    """
    parsed = read_vtk_unstructured_bin(path)
    faces = np.asarray(parsed["faces"])
    tags = np.asarray(parsed["cell_scalars"]["elemTag"]).astype(int)

    edges = np.sort(np.concatenate(
        [faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1)
    cell_of = np.tile(np.arange(faces.shape[0]), 3)
    uniq, inv, cnt = np.unique(edges, axis=0, return_inverse=True, return_counts=True)
    order = np.argsort(inv, kind="stable")
    starts = np.searchsorted(inv[order], np.arange(uniq.shape[0]))
    shared = cnt == 2
    a, b = cell_of[order][starts[shared]], cell_of[order][starts[shared] + 1]

    touching = np.sum(((tags[a] == 164) & (tags[b] == 115))
                      | ((tags[a] == 115) & (tags[b] == 164)))
    incident = np.sum((tags[a] == 164) | (tags[b] == 164))
    frac = touching / max(incident, 1)
    # Chance level is the fibrotic fraction, ~0.14-0.20. Observed is ~0.0008.
    assert frac < 0.01, f"tag 164 touches fibrosis on {frac:.4%} of its edges"


@needs_data
@pytest.mark.parametrize("path", MESHES[:3])
def test_released_fibre_field_is_flagged_degenerate(path):
    """The shipped fibre array is a single direction; the loader must say so out loud."""
    with pytest.warns(RuntimeWarning, match="degenerate"):
        mesh = load_uw_mesh(path)
    assert mesh.meta["fibres_are_degenerate"] is True
    assert mesh.meta["fibre_spread"] < DEGENERATE_FIBRE_SPREAD
    assert np.unique(mesh.fibres.round(6), axis=0).shape[0] == 1


@needs_data
def test_drop_tags_can_be_disabled_for_reproducibility(path=None):
    """``drop_tags=()`` must still reproduce the pre-fix substrate exactly."""
    p = MESHES[0]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        kept = load_uw_mesh(p, drop_tags=())
        dropped = load_uw_mesh(p)
    assert kept.meta["n_cells_dropped"] == 0
    assert dropped.meta["n_cells_dropped"] > 0
    assert dropped.faces.shape[0] < kept.faces.shape[0]
    # The legacy mapping inflated fibrosis by scoring the caps as half-fibrotic.
    assert kept.fibrosis.mean() > dropped.fibrosis.mean()


def test_tag_164_is_excluded_from_simulation_by_default():
    """Pure logic: whatever TAG_FIBROSIS says about 164, DROP_TAGS must remove it."""
    assert 164 in DROP_TAGS
    # The discredited mapping is retained only so it can be reproduced deliberately.
    assert TAG_FIBROSIS[164] == 0.5


def test_tags_to_fibrosis_averages_cell_tags_onto_vertices():
    """Two triangles sharing an edge: shared vertices average, private ones do not."""
    faces = np.array([[0, 1, 2], [1, 2, 3]])
    tags = np.array([111, 115])          # healthy, fibrotic
    fib = tags_to_fibrosis(tags, faces, 4, mapping=TAG_FIBROSIS)
    assert fib[0] == pytest.approx(0.0)   # only in the healthy triangle
    assert fib[3] == pytest.approx(1.0)   # only in the fibrotic triangle
    assert fib[1] == pytest.approx(0.5)   # shared
    assert fib[2] == pytest.approx(0.5)   # shared


def test_dead_tags_match_the_cohort_authors_legend():
    """Pure logic: 164 and 199 are both electrically dead per the legend of 2026-07-30.

    P. M. Boyle supplied the elemTag legend in correspondence: 111 atrial non-fibrotic,
    115 atrial fibrotic, 164 veins/valves "electrically non-conductive/dead", 199 ablation
    scar "also electrically non-conductive/dead". ``DEAD_TAGS`` records that; ``DROP_TAGS``
    records what this study removes, and the two differ by 199 alone because 199 never
    occurs in the pre-ablation meshes (see the next test).
    """
    assert DEAD_TAGS == (164, 199)
    assert 164 in DROP_TAGS
    assert 199 not in DROP_TAGS
    assert set(DROP_TAGS).issubset(set(DEAD_TAGS))


@needs_data
def test_no_pre_ablation_mesh_contains_ablation_scar():
    """The assumption that lets DROP_TAGS omit 199 without changing any result.

    If this ever fails, ``DROP_TAGS`` must gain 199 and every label cache must be
    rebuilt -- ``uw_drop_tags`` is part of the frozen-config hash, so the change is not
    free. Checked over the whole cohort rather than a sample, because the cost of the
    assumption being wrong is that ablation scar gets simulated as living tissue.
    """
    for path in MESHES:
        tags = np.asarray(read_vtk_unstructured_bin(path)["cell_scalars"]["elemTag"]).astype(int)
        present = set(np.unique(tags).tolist())
        assert 199 not in present, f"{path} contains ablation scar (tag 199)"
        assert present <= {111, 115, 164}, f"{path} has unexpected tags {present - {111, 115, 164}}"


@needs_data
@pytest.mark.parametrize("path", MESHES[:2])
def test_loading_a_mesh_with_unremoved_dead_tags_raises(path):
    """The guard that stops post-ablation meshes being simulated as if scar were tissue."""
    with pytest.raises(ValueError, match="electrically dead"):
        load_uw_mesh(path, drop_tags=(199,))     # 164 present but not removed
