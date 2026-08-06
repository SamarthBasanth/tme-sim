# Vivum (tmesim) Studio -- self-hosting image.
#
# Ships the simulator + Studio web UI. Does NOT ship any data cache: real
# tissue data (data_cache/) is gitignored and not baked into this image, in
# keeping with this project's "real data, explicitly acquired, never
# fabricated or silently bundled" discipline -- mount it at runtime instead
# (see docs/DEPLOYMENT.md).
#
# (That claim depends on .dockerignore, which excludes data_cache/ from the
# build context. Docker does not read .gitignore -- without .dockerignore this
# `COPY . .` would bake the whole 5.3 GB cache in and the paragraph above
# would be false. Do not delete it.)
#
# Build:
#   docker build -t vivum-studio .
# Run (local-only): bind the PUBLISHED port to the host loopback. Inside the
# container the server must still listen on 0.0.0.0 -- a socket bound to the
# container's own 127.0.0.1 is unreachable from a published port, since -p
# forwards into the container's network namespace.
#   docker run -p 127.0.0.1:8000:8000 -v $(pwd)/data_cache:/app/data_cache vivum-studio
# Run (network-exposed, auth required -- see docs/DEPLOYMENT.md):
#   docker run -p 8000:8000 -v $(pwd)/data_cache:/app/data_cache \
#     -e TMESIM_STUDIO_API_KEY=<your-key> \
#     vivum-studio
#
# NOTE the default CMD listens on 0.0.0.0, which app/server/__main__.py treats
# as non-loopback: with no TMESIM_STUDIO_API_KEY set it AUTO-GENERATES a key
# and prints it once rather than serving unauthenticated. Set the env var
# yourself unless you intend to read that key out of `docker logs`.

FROM python:3.11-slim

# CPython block-buffers stdout when it is a pipe, which is every `docker run`
# without a TTY. The startup banner -- including an auto-generated API key,
# which is printed exactly once and never persisted -- would otherwise sit in
# an unflushed buffer while uvicorn's stderr logs appear normally, making the
# container look healthy while every request 401s with a key nobody can read.
ENV PYTHONUNBUFFERED=1

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
#
# -c constraints.txt is what makes this image REPRODUCIBLE rather than merely
# buildable. pyproject.toml declares bounds, not versions, so without it two
# builds a month apart resolve different numpy/scipy/pandas and are not
# numerically identical -- which matters here because this repo asserts
# bit-exact reference results (see constraints.txt's header, and
# docs/REPRODUCE.md 1a). constraints.txt must therefore be COPYed into this
# same dependency layer, alongside pyproject.toml and BEFORE the install; it
# changes about as rarely, so the layer-cache behaviour is unchanged.
#
# Every version it pins has a cp311 manylinux x86_64 wheel (checked against
# PyPI, not assumed), so this resolves from wheels on python:3.11-slim with no
# compiler in the image.
COPY pyproject.toml constraints.txt ./
COPY tmesim/__init__.py tmesim/__init__.py
RUN pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir -c constraints.txt -e '.[app]'

COPY . .
# data_cache/ and app/runs/ are gitignored (real tissue data, run history) --
# create empty mount points so the app starts cleanly even with nothing
# mounted yet (it will just report those sections as unavailable, honestly,
# rather than fail to boot).
RUN mkdir -p data_cache app/runs

EXPOSE 8000

ENTRYPOINT ["python", "-m", "app.server"]
# 0.0.0.0 = "this container's interface", NOT "exposed to the world". The real
# boundary is the publish spec on `docker run` (-p 127.0.0.1:8000:8000 keeps it
# host-local). Binding 127.0.0.1 here would make the container unreachable
# through any published port. See the header for the API-key consequence.
CMD ["--host", "0.0.0.0", "--port", "8000"]
