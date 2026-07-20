# AtrialSpectralBench — reproducible CPU-only image.
# Build: docker build -t atrialspectralbench .
# Test:  docker run --rm atrialspectralbench          # runs pytest -q (default CMD)
# Run:   docker run --rm atrialspectralbench asb run --config configs/default.yaml
FROM python:3.11-slim

# Avoid interactive prompts and keep Python output unbuffered.
ENV PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

WORKDIR /app

# Copy the project and install it (core + dev extras) editable.
COPY . /app
RUN python -m pip install --upgrade pip && \
    python -m pip install -e ".[dev]"

# Default: run the analytic + property test suite (hard CI gates).
CMD ["pytest", "-q"]
