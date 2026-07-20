"""End-to-end pipeline test on a TINY synthetic cohort.

Runs :func:`asb.pipeline.run` with an overridden, minimal config (small cohort,
low-resolution meshes) into a temporary output directory and asserts that the
artifacts (``metrics.json``, ``results_report.md``, the five PNG figures) are
produced and that the reported AUC values are valid floats in ``[0, 1]``.
"""
from __future__ import annotations

import json
import os

import numpy as np

from asb.config import Config
from asb.pipeline import run


def _tiny_cfg(outputs_dir: str) -> Config:
    cfg = Config()
    # Small but NON-degenerate: this cohort/resolution reliably yields BOTH
    # inducible and non-inducible mock_ep labels, so the AUC/DeLong path is
    # actually exercised (a single-class cohort silently produces NaN AUCs).
    cfg.cohort.n_base = 5
    cfg.cohort.n_variants = 3
    cfg.cohort.n_nodes = 162  # icosphere level -> small, fast, but enough spread
    cfg.eval.n_null = 50      # keep the spatial-null resampling small and fast
    cfg.eval.n_bootstrap = 50
    cfg.outputs_dir = outputs_dir
    return cfg


def test_pipeline_writes_outputs_and_valid_aucs(tmp_path):
    out = str(tmp_path / "outputs")
    cfg = _tiny_cfg(out)

    metrics = run(cfg)

    # --- Artifacts on disk ------------------------------------------------- #
    metrics_path = os.path.join(out, "metrics.json")
    report_path = os.path.join(out, "results_report.md")
    assert os.path.isfile(metrics_path)
    assert os.path.isfile(report_path)
    for fig in (
        "fig_fiedler_atrium.png",
        "fig_sfi_vs_origin.png",
        "fig_eigen_spectrum.png",
        "fig_roc_panel.png",
        "fig_sensitivity_bars.png",
    ):
        assert os.path.isfile(os.path.join(out, fig)), fig

    # --- metrics.json is valid JSON and matches the return value ----------- #
    with open(metrics_path) as fh:
        on_disk = json.load(fh)
    assert on_disk["cohort"]["n_subjects"] == 15
    # Guard the headline path: a single-class cohort makes AUC/DeLong NaN and
    # would let a mis-calibrated label model pass silently.
    assert on_disk["cohort"]["both_classes_present"] is True

    # --- AUC values are valid floats in [0, 1] ----------------------------- #
    assert set(metrics["auc"].keys()) == {"fibrosis", "fibrosis+SFI"}
    for name, res in metrics["auc"].items():
        for kind in ("grouped_auc", "naive_auc"):
            val = res[kind]
            assert isinstance(val, float)
            assert np.isfinite(val), f"{name}/{kind} not finite"
            assert 0.0 <= val <= 1.0, f"{name}/{kind}={val} out of range"

    # --- The report is non-empty and carries the mock_ep disclaimer -------- #
    with open(report_path) as fh:
        report = fh.read()
    assert "AtrialSpectralBench" in report
    assert "mock_ep" in report


def test_pipeline_deterministic_auc(tmp_path):
    """Two runs of the same tiny config give identical, finite AUCs."""
    m1 = run(_tiny_cfg(str(tmp_path / "a")))
    m2 = run(_tiny_cfg(str(tmp_path / "b")))
    # Both classes present => real (non-NaN) AUCs, so equality is meaningful.
    assert m1["cohort"]["both_classes_present"] is True
    for name in ("fibrosis", "fibrosis+SFI"):
        for kind in ("grouped_auc", "naive_auc"):
            assert np.isfinite(m1["auc"][name][kind])
            assert m1["auc"][name][kind] == m2["auc"][name][kind]
