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
in the environment this was written in, so **no image has ever been built
from this `Dockerfile` and no container has ever been run.** What *has* been
verified, by reproducing the Dockerfile's steps outside a container in
throwaway Python 3.11 virtualenvs:

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

**Still unverified, and needing one real build to settle:** that
`docker build` completes at all (build context, `apt-get`, pip as root); that
linux `manylinux` wheels exist and install for `scipy`/`scikit-image`/
`imagecodecs`/`h5py`/`rdkit` on cp311 (every wheel tested here was
macOS/arm64); whether those wheels need shared libraries absent from
`python:3.11-slim` (the `libgomp1` line is precautionary, not measured, and
rdkit has historically wanted `libXrender`/`libXext` on slim images); that
`-p 8000:8000` actually reaches the app; bind-mount ownership and SELinux
relabeling; BuildKit cache reuse; PID-1 signal handling on `docker stop`; and
final image size.

One thing a build will **not** settle: this repo's bit-exact reference tests
were generated on macOS/arm64. On linux/amd64 they may differ for BLAS/SIMD
reasons unrelated to the container. Treat a bitwise mismatch inside the image
as unexplained until measured, not as a container defect.

**Also note, given the above:** this used to say that no dependency version was
pinned anywhere, so two builds a month apart resolved different numpy/scipy
versions and were not numerically identical. That was measured and true — a
fresh resolve gave numpy 2.4.6 / scipy 1.17.1 / pandas 3.0.5 against the
2.0.2 / 1.13.1 / 2.3.3 every published number was produced with, and 5 of 6
bitwise-reference failures went away once the versions were pinned back.

There is now a **`constraints.txt`** at the repo root pinning the numerically
load-bearing packages exactly, and the `Dockerfile` installs with
`-c constraints.txt`. Two builds a month apart now install the same numeric
stack. What that does *not* settle is the cross-architecture question in the
paragraph above: the references were generated on macOS/arm64, and pinning the
versions does not pin the CPU, the BLAS vendor or its thread count. See
[`REPRODUCE.md` §1a](REPRODUCE.md) for what is pinned, what deliberately is
not, and how to regenerate it.

Build it once yourself and open an issue if any of this doesn't hold.

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
