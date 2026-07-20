"""Configuration schema for AtrialSpectralBench.

All experiment knobs live here so runs are reproducible from a single YAML file.
FROZEN CONTRACT: field names are referenced across modules; keep in sync with
BUILD_SPEC.md.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

import yaml


@dataclass
class SFIConfig:
    """Spectral Fragility Index parameters."""

    n_monte_carlo: int = 256
    # Mean fractional conductance drop under diffuse perioperative stress.
    delta_w_mean_frac: float = 0.15
    # Coefficient of variation of the stochastic Delta-w field.
    delta_w_cov: float = 0.5
    # How strongly local fibrosis scales the expected uncoupling Delta-w.
    fibrosis_coupling: float = 1.0
    # Validity radius: require ||Delta L|| <= validity_safety * (lambda3 - lambda2)
    # before trusting the first-order (single-vector) SFI; else use subspace SFI.
    validity_safety: float = 0.25


@dataclass
class CohortConfig:
    """Synthetic cohort generation parameters."""

    n_base: int = 8         # number of base atrial anatomies (shape families)
    n_variants: int = 6     # statistical-shape variants per base
    n_nodes: int = 1200     # approximate mesh resolution per atrium
    seed: int = 0


@dataclass
class LabelConfig:
    """Mock-EP label generation parameters (development stand-in for openCARP)."""

    protocol: str = "S1S2"
    n_pacing_sites: int = 4
    source: str = "mock_ep"  # switch to 'opencarp' when the real solver is wired in


@dataclass
class EvalConfig:
    """Evaluation and statistics parameters."""

    n_splits: int = 5
    n_bootstrap: int = 1000
    n_null: int = 500       # rotational/shift spatial-null resamples for colocalization


@dataclass
class Config:
    seed: int = 0
    outputs_dir: str = "outputs"
    sfi: SFIConfig = field(default_factory=SFIConfig)
    cohort: CohortConfig = field(default_factory=CohortConfig)
    labels: LabelConfig = field(default_factory=LabelConfig)
    eval: EvalConfig = field(default_factory=EvalConfig)

    @classmethod
    def load(cls, path: str) -> "Config":
        with open(path) as fh:
            raw: dict[str, Any] = yaml.safe_load(fh) or {}
        return cls(
            seed=raw.get("seed", 0),
            outputs_dir=raw.get("outputs_dir", "outputs"),
            sfi=SFIConfig(**(raw.get("sfi") or {})),
            cohort=CohortConfig(**(raw.get("cohort") or {})),
            labels=LabelConfig(**(raw.get("labels") or {})),
            eval=EvalConfig(**(raw.get("eval") or {})),
        )

    def to_yaml(self, path: str) -> None:
        with open(path, "w") as fh:
            yaml.safe_dump(asdict(self), fh, sort_keys=False)
