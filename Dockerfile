# Vivum (tmesim) Studio -- self-hosting image.
#
# Ships the simulator + Studio web UI. Does NOT ship any data cache: real
# tissue data (data_cache/) is gitignored and not baked into this image, in
# keeping with this project's "real data, explicitly acquired, never
# fabricated or silently bundled" discipline -- mount it at runtime instead
# (see docs/DEPLOYMENT.md).
#
# Build:
#   docker build -t vivum-studio .
# Run (local-only, no auth needed):
#   docker run -p 8000:8000 -v $(pwd)/data_cache:/app/data_cache vivum-studio
# Run (network-exposed, auth required -- see docs/DEPLOYMENT.md):
#   docker run -p 8000:8000 -v $(pwd)/data_cache:/app/data_cache \
#     -e TMESIM_STUDIO_API_KEY=<your-key> \
#     vivum-studio --host 0.0.0.0

FROM python:3.11-slim

# libgomp1: scipy/scikit-image wheels on slim images sometimes need it at
# import time even though pip resolves clean; cheap to include, avoids a
# confusing runtime ImportError on first request instead of at build time.
RUN apt-get update && apt-get install -y --no-install-recommends \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Dependencies first (their own layer, so `docker build` doesn't reinstall
# scipy/rdkit on every source-code change). Editable installs via
# pyproject.toml alone need PEP 660 support (pip >= 21.3) -- verified locally
# that an old pip fails here with a confusing "setup.py not found" error even
# though this project has no setup.py by design; upgrading pip first avoids
# depending on whatever pip version happens to ship with the base image.
COPY pyproject.toml ./
COPY tmesim/__init__.py tmesim/__init__.py
RUN pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir -e '.[app]'

COPY . .
# data_cache/ and app/runs/ are gitignored (real tissue data, run history) --
# create empty mount points so the app starts cleanly even with nothing
# mounted yet (it will just report those sections as unavailable, honestly,
# rather than fail to boot).
RUN mkdir -p data_cache app/runs

EXPOSE 8000

ENTRYPOINT ["python", "-m", "app.server"]
CMD ["--host", "127.0.0.1", "--port", "8000"]
