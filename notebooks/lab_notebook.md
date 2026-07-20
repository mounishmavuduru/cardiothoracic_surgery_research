# AtrialSpectralBench — Lab Notebook

> **Why this file exists.** This is a dated research notebook. Science-fair and
> peer judges routinely inspect the lab notebook to verify that the work was done
> incrementally, honestly, and reproducibly — that hypotheses were written down
> *before* results, that dead ends and negative results were recorded, and that
> every claimed figure traces back to a seeded, config-driven run. Keep it append-only:
> add new dated entries, never silently rewrite old ones. If an earlier entry was wrong,
> add a later entry that corrects it and say so.

## How to keep this notebook

- **One entry per working session**, newest at the bottom (or top — pick one and be
  consistent). Always start with a real date (`YYYY-MM-DD`).
- Record: what you set out to do, what you actually did, the exact command / config /
  seed, what you observed (numbers, not vibes), what surprised you, and the next step.
- Log **negative and null results** with equal weight — a clean null
  ("SFI ≡ re-encoded fibrosis") is an accepted, reportable outcome for this project.
- Reference the commit hash and the config file for anything reproducible.
- Never invent data or backfill a date. If you did nothing on a day, there is no entry.

---

## Entry template (copy for each new session)

### YYYY-MM-DD — <short title>

- **Goal:** what question or task this session targets.
- **Setup:** branch / commit hash, config file, seed(s), environment notes.
- **Did:** the concrete steps taken (commands run, files changed, experiments launched).
- **Observed:** measured results — metrics, AUC/ΔAUC, plots produced, test pass/fail,
  error messages. Paste the numbers, not a summary.
- **Interpretation:** what the numbers mean; whether they support or falsify the
  pre-registered hypothesis; caveats and threats to validity.
- **Surprises / dead ends:** anything unexpected, including things that did not work.
- **Next:** the single most important next action.

---

## Entries

### YYYY-MM-DD — Phase 0: scaffold & environment

- **Goal:** stand up the repository skeleton so every later phase is one command away.
- **Setup:** branch `claude/third-idea-project-plan-*`; Python 3.11 CPU-only venv at
  `.venv`; config `configs/default.yaml` emitted from `asb.config.Config`.
- **Did:** created the venv, editable-installed the package with `[dashboard,dev]`
  extras, verified `import asb, asb.types, asb.config` and the default-config emit,
  wrote LICENSE / README / Makefile / CI / Dockerfile / this notebook / the
  pre-registration document.
- **Observed:** _fill in — package version, numpy/scipy/sklearn versions, whether the
  backbone import and config emit succeeded._
- **Interpretation:** _fill in._
- **Surprises / dead ends:** _fill in._
- **Next:** implement the spectral engine (Phase 1) and its analytic CI gates
  (`λ_k = 2 − 2cos(kπ/N)` for path/ring/grid).
