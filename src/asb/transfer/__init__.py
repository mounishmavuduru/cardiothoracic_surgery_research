"""GM4 second excitable medium — porting the spectral-fragility calculus off the heart.

Contains a 2-D excitable **neural-network** generator (:mod:`asb.transfer.network`)
and a FitzHugh--Nagumo instability simulator (:mod:`asb.transfer.fhn`). The point of
this package is that the fragility calculus (``asb.sfi``), the spectral engine
(``asb.spectral``) and the substrate baselines (``asb.baselines``) are **graph-
universal** and are reused here *verbatim* on a completely different excitable medium.
"""
