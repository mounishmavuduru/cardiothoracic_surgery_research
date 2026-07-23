# AtrialSpectralBench — developer convenience targets.
# Uses a local .venv (Python 3.11, CPU-only). All stochastic code is seeded.

PYTHON ?= python3
VENV   := .venv
BIN    := $(VENV)/bin
CONFIG ?= configs/default.yaml

.PHONY: venv install test run dashboard paper overleaf-zip clean

## venv: create the local virtual environment
venv:
	$(PYTHON) -m venv $(VENV)
	$(BIN)/python -m pip install --upgrade pip

## install: editable install with dashboard + dev extras
install: venv
	$(BIN)/python -m pip install -e ".[dashboard,dev]"

## test: run the analytic + property test suite (hard CI gates)
test:
	$(BIN)/python -m pytest -q

## run: end-to-end pipeline on the synthetic cohort
run:
	$(BIN)/python -m asb.cli run --config $(CONFIG)

## dashboard: launch the Streamlit dashboard (needs the [dashboard] extra)
dashboard:
	$(BIN)/streamlit run src/asb/dashboard/app.py

## paper: compile the LaTeX manuscript to docs/paper/manuscript.pdf (needs a TeX install)
paper:
	cd docs/paper && latexmk -pdf -interaction=nonstopmode manuscript.tex

## overleaf-zip: bundle manuscript.tex + figures/ (top-level) for Overleaf "Upload Project"
overleaf-zip:
	cd docs/paper && rm -f AtrialSpectralBench_overleaf.zip && \
	  zip -r AtrialSpectralBench_overleaf.zip manuscript.tex figures/ -x "figures/.*"

## clean: remove build/venv/cache artifacts
clean:
	rm -rf $(VENV) *.egg-info src/*.egg-info .pytest_cache
	find . -type d -name __pycache__ -prune -exec rm -rf {} +
