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
    cfg.cohort.n_base = 3
    cfg.cohort.n_variants = 2
    cfg.cohort.n_nodes = 42  # nearest icosphere resolution -> tiny, fast meshes
    cfg.eval.n_null = 50     # keep the spatial-null resampling small and fast
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
    assert on_disk["cohort"]["n_subjects"] == 6

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
    """Two runs of the same tiny config give identical AUCs."""
    m1 = run(_tiny_cfg(str(tmp_path / "a")))
    m2 = run(_tiny_cfg(str(tmp_path / "b")))
    assert m1["auc"] == m2["auc"]
