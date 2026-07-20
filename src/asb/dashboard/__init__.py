"""AtrialSpectralBench interactive dashboard package.

The dashboard depends on Streamlit, which is an *optional* extra. Importing this
package (or :mod:`asb.dashboard.app`) never imports Streamlit at module load time, so
the core ``asb`` package works with or without the dashboard extra installed. Launch
the app with::

    streamlit run src/asb/dashboard/app.py
"""
from __future__ import annotations

__all__: list[str] = []
