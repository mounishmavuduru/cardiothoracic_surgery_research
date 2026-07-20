"""Command-line interface for AtrialSpectralBench.

Exposes the ``asb`` console script (declared in ``pyproject.toml`` as
``asb = "asb.cli:main"``). The single subcommand today is ``run``::

    asb run --config configs/default.yaml

which executes the full synthetic-data benchmark in :func:`asb.pipeline.run` and
writes ``metrics.json``, ``results_report.md`` and the figures under the config's
``outputs_dir``. The same benchmark is reachable via ``python -m asb.pipeline``.
"""
from __future__ import annotations

import argparse
from typing import List, Optional

from asb.config import Config
from asb.pipeline import run

__all__ = ["main", "build_parser"]


def build_parser() -> argparse.ArgumentParser:
    """Construct the ``asb`` argument parser."""
    parser = argparse.ArgumentParser(
        prog="asb", description="AtrialSpectralBench command-line interface.")
    sub = parser.add_subparsers(dest="command", required=True)

    run_p = sub.add_parser(
        "run", help="Run the synthetic-data benchmark end to end.")
    run_p.add_argument(
        "--config", type=str, default=None,
        help="Path to a YAML config (defaults to the built-in Config()).")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    """Entry point for the ``asb`` console script.

    Parameters
    ----------
    argv : list of str, optional
        Argument vector (defaults to ``sys.argv[1:]``).

    Returns
    -------
    int
        Process exit code (0 on success).
    """
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "run":
        cfg = Config.load(args.config) if args.config else Config()
        metrics = run(cfg)
        print(f"Wrote {metrics['metrics_path']}")
        print(f"Wrote {metrics['report_path']}")
        for name, p in metrics["figures"].items():
            print(f"Wrote {p}")
        return 0

    parser.error(f"unknown command {args.command!r}")
    return 2  # unreachable; argparse.error exits.


if __name__ == "__main__":
    raise SystemExit(main())
