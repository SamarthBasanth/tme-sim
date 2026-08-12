# Deployment — running Studio for a team, not just yourself

Studio (`app/server` + `app/web`) was built as a local diagnostic workstation
tool: `python -m app.server` binds to `127.0.0.1` by default and, until this
document, had no authentication at all. This page covers what changes when a
pharma/biotech or research team wants to reach it from somewhere other than
the machine it's running on.

**What this does NOT cover.** Any use involving real, identifiable patient
data or a clinical decision. Every dataset this project ships with is open,
public research data (10x Genomics/GEO Xenium and Visium sections) with no
patient-identifiable information — see `docs/DATA_INVENTORY.md`. If your team
plans to reconstruct anything else, that's a different risk posture than what
is set up here, and needs its own review before you follow this guide.

---

## 1. The three ways to run this

### 1a. Local, single-user (the default, unchanged)

```bash
vivum studio
```

Binds to `127.0.0.1`. No API key is required or generated. Nothing here
changes this path — every test and every piece of documentation that predates
this page still applies exactly as written.

### 1b. On your own network, one shared key

```bash
vivum studio --host 0.0.0.0 --api-key <a-real-secret>
```

or generate one automatically:

```bash
vivum studio --host 0.0.0.0
```

```
  api key     Xk3n... [generated (not persisted -- save it now)]
              every request needs header: X-API-Key: Xk3n...
```

Passing any `--host` other than `127.0.0.1`/`localhost`/`::1` with no
`--api-key` and no `TMESIM_STUDIO_API_KEY` environment variable **auto-
generates a random key and prints it once** rather than serving
unauthenticated — this project's own principle ("a failure a user caused
comes back as a diagnostic, never a silent default") applied to the
deployment surface, not just the protocol language. The generated key is
**not written to disk anywhere** — copy it from the terminal output before
you lose the process.

Open `http://<host>:<port>/` in a browser. The frontend shows an **"Access
key required"** gate on first load if the server needs one; paste the key in
once and it's remembered in that browser's `localStorage` for next time. A
wrong key shows an inline error and lets you retry, rather than failing
silently or falling back to fake demo data (the earlier version of this app
would degrade to synthetic stub data on ANY unreachable/rejected request,
which would have silently shown plausible-looking fake numbers instead of a
real reconstruction — fixed as part of adding this gate).

If you'd rather set the key yourself (e.g. to something a whole team already
shares) instead of taking the generated one:

```bash
TMESIM_STUDIO_API_KEY=<your-key> vivum studio --host 0.0.0.0
```

An already-exported `TMESIM_STUDIO_API_KEY` is always respected over the
auto-generate-or-refuse logic above, so this composes correctly with Docker's
`-e` / systemd's `Environment=` / your shell profile.

**What this is, and isn't.** One shared secret, checked on every `/api/*`
request. It stops a stranger on the network from using the tool; it does
**not** give you separate accounts, audit logs, or per-user permissions. If
your team needs that, this is the wrong layer to build it in — put a real
auth proxy in front (see §3) rather than extending this shared-key check.

### 1c. Docker

```bash
docker build -t vivum-studio .
docker run -p 8000:8000 \
  -v $(pwd)/data_cache:/app/data_cache \
  -e TMESIM_STUDIO_API_KEY=<your-key> \
  vivum-studio --host 0.0.0.0
```

**Honesty note on this section specifically.** Docker is still not installed
in the environment this was written in, so everything below was first
established *without* building an image — by reproducing the Dockerfile's
steps outside a container in throwaway Python 3.11 virtualenvs. That is no
longer the only evidence: `.github/workflows/docker.yml` has since built and
exercised the image on a real linux/amd64 runner, and the measurements are in
"What the build settled" below. The static audit is recorded first because it
is what found the defects. What it verified:

* The dependency layer resolves and installs. `pip install -e '.[app]'` from
  `pyproject.toml` plus a stub `tmesim/__init__.py` succeeds on Python 3.11
  (it needs pip >= 21.3 / PEP 660 — confirmed an old pip fails with a
  confusing "setup.py not found", which is why the Dockerfile upgrades pip
  first).
* The layer-caching trick is sound: installing editable against a stub
  package and copying the rest of `tmesim/` in afterwards still imports
  correctly, because setuptools' editable finder maps the top-level package
  directory.
* The `data_cache/` mount target is right. The loaders resolve
  `data_cache/...` relative to the working directory, so `WORKDIR /app` plus
  `-v $(pwd)/data_cache:/app/data_cache` lands where the code looks. With the
  cache absent, `/api/sections` reports each section unavailable *with a
  reason* rather than failing to boot. (Corollary: do not override the
  container's working directory.)
* `app/runs/`, the static front end at `/`, and the `TMESIM_STUDIO_API_KEY`
  gate (401 without the header, 200 with it) all behave correctly at the
  paths the image uses.
* Python 3.11 — the base image's version, newer than the 3.9 used for
  development — runs the suite green when the full dependency set is present.

**Four real defects were found and fixed this way, without ever building the
image**, which is worth stating because it is the argument for doing this
audit at all:

1. **There was no `.dockerignore`.** Docker does not read `.gitignore`, so
   `COPY . .` was baking in `data_cache/` (5.3 GB of real tissue), `app/runs/`
   (647 MB of run history) and a 618 MB macOS virtualenv — making this
   section's own "does NOT ship any data cache" claim false, and the build
   context large enough to be impractical. Fixed by adding `.dockerignore`.
2. **The default `CMD` bound the container's own loopback** (`--host
   127.0.0.1`), which is unreachable through a published port, so the
   documented `docker run -p 8000:8000 …` produced a container nothing could
   talk to. Now binds `0.0.0.0`; use `-p 127.0.0.1:8000:8000` to keep it
   host-local.
3. **`rdkit` and `h5py` were undeclared runtime dependencies.** Measured on an
   otherwise-identical environment: without `rdkit`, 61 tests fail, `ribociclib`
   and `palbociclib` compile to `mw=0.0` and are permanently unusable, and
   bit-exact reference results shift; without `h5py`, GBM-Space ingestion
   cannot even import. Both now declared in `pyproject.toml`.
4. **The auto-generated API key never reached `docker logs`.** CPython
   block-buffers stdout when it is a pipe, so the startup banner (and the key,
   which is printed once and never persisted) stayed in an unflushed buffer
   while uvicorn's stderr logs appeared normally — a container that looks
   healthy while 401ing every request with a key nobody can read. Fixed with
   `ENV PYTHONUNBUFFERED=1`.

**What the build settled.** The `image` job of `.github/workflows/docker.yml`
passed on `ubuntu-latest` in run `31147942003` (`9ae3e10`, 2026-08-07) and
again in run `31152915639` (`70b2ed6`, 2026-08-07). Both are `linux/amd64`
only — that is the single platform the job builds. Measured:

* **`docker build` completes.** Run `31147942003` built every layer from
  cold in 100 s wall (`4.50 kB` Dockerfile, `7.35 MB` build context — the
  `.dockerignore` working, against the ~6.5 GB context it would otherwise be).
  `apt-get` and pip-as-root both succeed; pip warns about running as root and
  proceeds.
* **The cp311 `manylinux` wheels exist and install.** In that same cold build
  the dependency layer took 28.5 s and resolved to wheels, not source:
  `numpy-2.0.2`, `scipy-1.13.1`, `scikit_image-0.24.0`, `scikit_learn-1.6.1`,
  `imagecodecs-2024.12.30` and `h5py-3.14.0` as
  `cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64`,
  `rdkit-2025.9.2` as `cp311-cp311-manylinux_2_28_x86_64`, `pandas-2.3.3` as
  `manylinux_2_24_x86_64.manylinux_2_28_x86_64`. The only line reading
  `Building wheels for collected packages:` is `tmesim` itself, so nothing
  compiles — which matters, because the image ships no compiler. The versions
  installed are exactly `constraints.txt`'s.
* **The slim image's system libraries are sufficient.** Inside the built
  image — which installs `libgomp1` and nothing else — all of `numpy`,
  `scipy`, `scipy.ndimage`, `skimage`, `skimage.measure`, `tifffile`,
  `imagecodecs`, `h5py`, `rdkit`, `rdkit.Chem`, `rdkit.Chem.Descriptors`,
  `sklearn`, `networkx`, `pandas`, `matplotlib`, `fastapi`, `uvicorn` and
  `httpx` import. rdkit in particular imports **without** `libXrender` or
  `libXext`, and `Descriptors.MolWt` on ribociclib's SMILES returns
  `434.5480000000002` in-image. Note the shape of this result: the installed
  set is *sufficient*. Nothing removes `libgomp1` to test whether it is
  *necessary*, so that line is still precautionary.
* **The published port reaches the app.** `docker run -d -p
  127.0.0.1:8000:8000` then polling: `/api/health` returned 200 after 2 s.
  This is the assertion that would fail if the `CMD` regressed to binding the
  container's own loopback (defect 2 above), because a container-loopback
  socket is unreachable through any published port.
* **BuildKit cache reuse works.** Run `31152915639` imported the `type=gha`
  cache manifest and reported `CACHED` for the `apt-get`, `WORKDIR`, both
  `COPY`s and the whole pip layer, restoring the 266 MB site-packages layer
  rather than re-resolving it. Total build fell from 100 s to 53 s.
* **PID 1 handles `SIGTERM`.** `docker stop -t 30` returned in 1 s with exit
  code 0 (0 s in run `31147942003`). Not 137, so the exec-form `ENTRYPOINT`
  is delivering the signal to uvicorn rather than letting the 30 s timeout
  expire into a `SIGKILL`.
* **Final image size: 980 MB** — 980,609,100 bytes in run `31152915639`,
  980,558,715 in run `31147942003`.

The same job also re-confirms the four defects above stay fixed, in the image
rather than by reasoning about it: the `.dockerignore` exclusions hold and
`data_cache/`/`app/runs/` exist but are **empty** — emptiness, not absence,
being the real test, since a deleted `.dockerignore` line would let the
`mkdir -p` succeed over copied tissue data (defect 1); `/api/health` answers
through the published port (defect 2); ribociclib and palbociclib compile to
`mw=434.548` and `mw=447.543` with non-null diffusivity and no `reason`
(defect 3); and a 32-character generated key is readable from `docker logs`
of a still-*running* container — the only state in which the buffering bug is
detectable — and then authenticates (defect 4). Alongside them, the §1b key
gate answers 401 / 401 / 200 for absent / wrong / correct `X-API-Key` while
`GET /` stays ungated at 200, and `/api/sections` returns its 7 built-ins
with the cache absent.

**Still unverified, because the workflow has no step that exercises it:**

* **Bind mounts.** Every container the job runs has `data_cache/` absent —
  nothing is ever mounted. So bind-mount ownership, SELinux relabeling (the
  `:z`/`:Z` question on RHEL/Fedora hosts), and whether a *populated*
  `-v $(pwd)/data_cache:/app/data_cache` actually reads back correctly inside
  a container are all still only verified outside one.
* **Named-volume persistence.** The `-v vivum_runs:/app/app/runs` form below
  is never run.
* **`docker exec <container> vivum acquire`.** Never run.
* **Whether `libgomp1` is required.** Sufficient, not shown necessary.
* **Any architecture other than `linux/amd64`.** The job pins
  `platforms: linux/amd64`; no arm64 image has ever been built.

**One earlier run of this workflow failed, and it was not the image.** Run
`31122338583` (`99bb366`, 2026-08-06) is recorded as a failure, but both of
its jobs are `cancelled` with `runner_id: 0`, no runner name, no steps
recorded, and a log download that 404s — **not one step of it ever
executed**. Nothing in this repository ran, so nothing in this repository can
be blamed for it, and nothing was changed to fix it: the only relevant commit
between that run and the first green one (`0fef386`) added `constraints.txt`
and the `-c` flag, which cannot affect whether a job is assigned a runner.
Both runs since were picked up within ~3 s of being created and completed, so
the condition has not recurred — but it was never diagnosed either, and a
repeat is not ruled out. A red on this workflow with no step logs means this,
not a broken image; re-run it before investigating the Dockerfile.

**The cross-architecture question is now measured, and the answer is that the
references do not hold.** This repo's bit-exact reference values were
generated on macOS/arm64. The informational `reference-tests-linux` job ran in
both green runs and was **red in both, with byte-identical numbers** — so the
linux/amd64 values are themselves deterministic, they simply are not the
stored ones. `6 failed, 203 passed, 81 skipped, 12 deselected` both times.
The job's fingerprint printed numpy 2.0.2, scipy 1.13.1, scikit-learn 1.6.1,
scikit-image 0.24.0 and pandas 2.3.3 — exactly `constraints.txt`'s versions,
**so this is not version drift**. What differs is the machine: CPython
3.11.15 on Linux x86_64 / glibc 2.39, numpy built against `scipy-openblas`
(`OpenBLAS 0.3.27 USE64BITINT DYNAMIC_ARCH NO_AFFINITY Zen MAX_THREADS=64`).
Five of the six are tolerance-free comparisons in
`tests/test_time_resolved_pk.py`, and the divergence is narrower than
"the references fail":

* Every one of them is an **AUC accumulator** — `res.extras["mean_auc"]` in
  four parametrised cases, `res.auc.sum()` in the fifth. A summation order,
  which is what a different BLAS/SIMD reduction changes.
* `kill_fraction` is asserted on the line immediately above `mean_auc` in all
  four parametrised cases and **passed bit-for-bit** every time.
* The magnitudes are 10, 11, 14, 80 and 81 ULP — 1.2e-15 to 1.4e-14
  relative. E.g. `20642672.153042287` against the stored
  `20642672.15304199`.

The sixth red is
`test_visium_pseudobulk_mgmt_is_recovered_and_labelled_pseudobulk` failing
with "no cached Visium GBM sections found" — absent data, not floats, and a
gap in that job's `--deselect` list rather than a result.

This narrows the caveat; it does not remove it. What is measured is the
*synthetic* references: the failing cases run on `_build_selftest_micro()` and
`make_flat_micro()`, while the real-tissue bitwise test is deselected because
`data_cache/` is not available in CI. So the real-tissue references remain
unmeasured across architectures. Treat a bitwise mismatch inside the image as
this known ~1e-14 architecture difference **only** if it is of this magnitude
in one of these accumulators; anything larger, or in `kill_fraction`, is
unexplained and is not accounted for by anything above.

**Also note, given the above:** this used to say that no dependency version was
pinned anywhere, so two builds a month apart resolved different numpy/scipy
versions and were not numerically identical. That was measured and true — a
fresh resolve gave numpy 2.4.6 / scipy 1.17.1 / pandas 3.0.5 against the
2.0.2 / 1.13.1 / 2.3.3 every published number was produced with, and 5 of 6
bitwise-reference failures went away once the versions were pinned back.

There is now a **`constraints.txt`** at the repo root pinning the numerically
load-bearing packages exactly, and the `Dockerfile` installs with
`-c constraints.txt`. Two builds a month apart now install the same numeric
stack. What that does *not* settle is the cross-architecture question above —
and the CI result is the demonstration, not a worry about one: the linux job
installs with the same `-c`, its fingerprint printed exactly the reference
versions, and the AUC references diverged anyway. Pinning the versions does
not pin the CPU, the BLAS vendor or its thread count. See
[`REPRODUCE.md` §1a](REPRODUCE.md) for what is pinned, what deliberately is
not, and how to regenerate it.

Build it on your own host and open an issue if any of this doesn't hold —
the numbers above are one runner image on `linux/amd64`, not a guarantee
about yours.

`data_cache/` is gitignored and not baked into the image (this project's
"real data, explicitly acquired, never silently bundled" rule extends to the
container) — mount your own populated cache, or fetch it inside the running
container. Note the `ENTRYPOINT` is `python -m app.server`, so
`docker run vivum-studio vivum acquire` does **not** work (it becomes
`python -m app.server vivum acquire` and exits 2). Use:

```bash
docker exec <container> vivum acquire
```

**Writable state does not survive container replacement.** Run history,
`app/runs/_ingested_sections.json` (the file §1d tells you to populate),
`app/projects/`, `app/reports/` and `app/uploads/` all live inside the
image's filesystem. Mount them if you want them to persist:

```bash
docker run -p 127.0.0.1:8000:8000 \
  -v $(pwd)/data_cache:/app/data_cache \
  -v vivum_runs:/app/app/runs \
  vivum-studio
```

---

## 1d. Restart to pick up sections ingested by another process

`app/server/ingest_registry.py` reads the ingested-section table (real
GBM-Space / breast / etc. sections, on top of the 7 built-in ones) **once per
process** and caches it in memory. If a long-running Studio server is already
up and a separate process registers a new section — a data-ingestion script,
or `POST /api/sections/scan` hit against a short-lived second instance — the
running server keeps serving its stale, older list until restarted. There is
no live-reload. If `/api/sections` looks thinner than you expect (only the 7
built-ins, none of the sections you know you've ingested), restart the
server before assuming the data is missing.

To register a section that already has real data in `data_cache/` but isn't
in the dropdown yet (this repo's breast-cancer Xenium section, for example,
was ingested via direct Python calls in earlier work rather than through this
route, so it needed registering separately):

```bash
curl -X POST http://127.0.0.1:8000/api/sections/scan \
  -H "Content-Type: application/json" \
  -d '{"path": "data_cache/<your-section-dir>", "register": true, "section_id": "<id>", "label": "<label>"}'
```

This persists to `app/runs/_ingested_sections.json` (gitignored, local state)
and is picked up immediately by the process that ran the scan — other
already-running processes still need the restart above.

## 2. What is NOT handled for you

* **TLS.** This serves plain HTTP. Put a reverse proxy in front for anything
  crossing an untrusted network (Caddy and nginx both terminate TLS for a
  single backend in ~10 lines of config; this project deliberately doesn't
  bundle one, to avoid maintaining certificate/domain logic orthogonal to the
  simulator itself).
* **Multi-user accounts.** One shared key, not per-user identity. See §1b.
* **Rate limiting / abuse protection.** None. A reverse proxy is also the
  right place for this if you need it.
* **Horizontal scaling.** `app/server/cache.py`'s reconstruction cache and
  `app/server/jobs.py`'s job queue are both single-process, in-memory
  structures — running two replicas behind a load balancer would give each
  one its own cache and queue, silently duplicating reconstruction work and
  splitting run history. Run one instance.

---

## 3. A reasonable reference setup

For a small pharma/research team (a handful of people, trusted network or VPN,
no need for per-user accounts):

```
[team's machines] --HTTPS--> [Caddy/nginx, TLS termination] --HTTP, localhost--> [vivum studio --host 127.0.0.1]
```

i.e. keep Studio itself bound to loopback (no API key needed at all — the
proxy is the boundary), and let the reverse proxy own TLS + whatever
authentication your organization already has for internal tools (basic auth,
an SSO gateway, a VPN requirement). This is usually less work than exposing
Studio directly with `--host 0.0.0.0` and its own key, and centralizes auth
in infrastructure you're probably already running rather than in this
application.

Use §1b's own API key instead when there's no existing reverse proxy /
internal-tools auth to sit behind, and a shared secret is an acceptable bar
for who can reach it.
