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

# SECURITY: this image ran as root (uid=0) until it was measured, which is the
# exact opposite of the control docs/DEPLOYMENT.md 2 names for the one route
# this application deliberately does not confine. `POST /api/sections/scan`
# takes an arbitrary filesystem path and reports on any directory the SERVER
# PROCESS can read; there is no allowlisted root, and the API key makes it
# authenticated, not safe. The documented mitigation is therefore "run the
# server as a user that cannot read anything you would mind a tester listing"
# -- and as root, `scan {"path":"/etc"}` and `scan {"path":"/root"}` both
# answered 200 with a valid key. A non-root uid is the Docker half of that
# advice; it does not confine the path, it shrinks what the path can reach.
#
# The uid is pinned (10001, not "whatever useradd picks next") because it is
# part of the runtime contract: a HOST directory bind-mounted over a writable
# path keeps the HOST's ownership, so an operator who wants to persist
# app/runs/ that way needs a number to chown to. See docs/DEPLOYMENT.md 1c.
#
# Deliberately placed AFTER the dependency layer rather than next to apt-get:
# above `RUN pip install` it would invalidate the one layer worth caching (a
# ~4-minute reinstall of the whole scipy/rdkit stack) for a line that never
# changes. Here it costs a fraction of a second per source-code rebuild.
#
# HOME is pinned because Docker resolves it from /etc/passwd for USER but not
# every runtime that consumes an OCI image does, and matplotlib writes a font
# cache into HOME on first import.
RUN useradd --create-home --uid 10001 --shell /usr/sbin/nologin vivum
ENV HOME=/home/vivum

# data_cache/ and app/runs/ are gitignored (real tissue data, run history) --
# create empty mount points so the app starts cleanly even with nothing
# mounted yet (it will just report those sections as unavailable, honestly,
# rather than fail to boot).
#
# The chown is what makes USER below actually work, and it is deliberately
# narrow. Everything this application writes resolves to <repo>/app/<something>
# (app/server/cache.py runs_dir + private_makedirs, report.py, projects.py,
# compounds.py, tmesim/protocol/resolve.py uploads_dir) -- app/runs/,
# app/projects/, app/reports/, app/compounds/ and app/uploads/, all created on
# demand at 0700. So /app/app must be writable BY THE DIRECTORY, not just by
# its contents: four of those five do not exist in the image at all and are
# mkdir'd on first use. They are not pre-created here on purpose -- the CI
# workflow asserts they are ABSENT, because their presence is how a
# .dockerignore regression that bakes the operator's unpublished SMILES into a
# layer would show up.
#
# data_cache/ is chowned for `docker exec <container> vivum acquire`, which
# downloads into it. A read-only bind mount over it ignores this ownership
# (the host's wins), which is the documented and recommended shape.
#
# /app itself stays root-owned and read-only to the app: the source tree, the
# static front end and site-packages are not things a running server should be
# able to rewrite.
RUN mkdir -p data_cache app/runs \
    && chown -R vivum:vivum /app/data_cache /app/app

# --- the one tissue section, FETCHED FROM SOURCE AT BUILD TIME -------------- #
# Without this the image starts fine and can run NOTHING: .dockerignore
# excludes data_cache/ (5.3 GB locally), so a hosted container has no tissue and
# every section is unavailable. That is correct for a local deploy that bind-
# mounts the host's cache, and useless for a hosted one, which has no host.
#
# FETCHED, NOT COPIED, and the distinction is deliberate. `vivum acquire`
# downloads xenium_319 from the ORIGINAL source, so the image carries
# data obtained under the source's own terms rather than a redistribution of one
# developer's local copy. It also cannot silently ship a section that has been
# modified on someone's laptop.
#
# WHAT THIS ACTUALLY COSTS: 87 MB, MEASURED IN THE IMAGE. `--sections
# xenium_319` selects the SPATIAL section (16 MB) and `acquire` fetches the four
# TABULAR sources alongside it regardless -- Ivy GAP 37 MB, TCGA expression
# 33 MB, GDSC dose-response 1.9 MB, TCGA clinical 64 KB. Those are not optional
# padding: the GDSC file is where a compound's MEASURED IC50 comes from, and
# without it the drug panel falls back to flagged placeholders.
#
# The commit that added this step said "~16.6 MB", which was the section's
# download size mistaken for the step's cost. .github/workflows/docker.yml
# asserts an allowlist of exactly these five entries plus a 300 MB ceiling, so
# the real figure is now pinned by CI rather than by a comment.
#
# ARG so a local build can skip it: --build-arg FETCH_SECTION=0 gives the old
# no-tissue image for the bind-mount workflow, which is still what
# docs/DEPLOYMENT.md describes.
ARG FETCH_SECTION=1
RUN if [ "$FETCH_SECTION" = "1" ]; then \
      python -m tmesim.cli acquire --sections xenium_319 --timeout 180 \
      && chown -R vivum:vivum /app/data_cache ; \
    else \
      echo "[build] FETCH_SECTION=0 -- image ships with NO tissue; bind-mount one" ; \
    fi

EXPOSE 8000

# Everything above this line runs as root (apt, pip, COPY). Everything the
# SERVER does runs as vivum. Verified in a real container, not assumed: the
# five 0700 runtime trees are created and written, a read-only data_cache
# bind mount is readable, and a named volume on /app/app/runs inherits this
# ownership so it is writable with no host-side chown.
USER vivum

ENTRYPOINT ["python", "-m", "app.server"]
# 0.0.0.0 = "this container's interface", NOT "exposed to the world". The real
# boundary is the publish spec on `docker run` (-p 127.0.0.1:8000:8000 keeps it
# host-local). Binding 127.0.0.1 here would make the container unreachable
# through any published port. See the header for the API-key consequence.
CMD ["--host", "0.0.0.0", "--port", "8000"]
