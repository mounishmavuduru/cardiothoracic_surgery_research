"""AtrialSpectralBench (asb) — a spectral-fragility benchmark for surgery-triggered
atrial fibrillation.

The protected novel seed is the closed-form Fiedler-value sensitivity
dλ₂/dw_ij = (φ_i − φ_j)², aggregated under a stochastic fibrosis-weighted diffuse
perioperative-uncoupling field into a per-region Spectral Fragility Index (SFI).

Everything is in-silico: zero patients, zero IRB. Ground-truth labels come from an
EP simulator (openCARP when available; an in-repo `mock_ep` stand-in otherwise).
`mock_ep` labels are a development convenience and are NEVER clinical POAF.
"""
from __future__ import annotations

__version__ = "0.1.0"

from asb.types import AtrialGraph, AtrialMesh, InducibilityLabel  # noqa: E402
from asb.config import Config  # noqa: E402

__all__ = ["AtrialGraph", "AtrialMesh", "InducibilityLabel", "Config", "__version__"]
