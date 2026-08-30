# Deployment — running Studio for a team, not just yourself

Studio (`app/server` + `app/web`) was built as a local research workstation
tool: `python -m app.server` binds to `127.0.0.1` by default and, until this
document, had no authentication at all. This page covers what changes when a
pharma/biotech or research team wants to reach it from somewhere other than
the machine it's running on.

**What this does NOT cover: identifiable patient data, or a clinical
decision.** Vivum is not built for either and must not be used with either. It
has one shared secret and no accounts, no audit log, no encryption at rest, and
no per-record access control — §2 is the full list of what it does not do for
you. Software with those gaps has no business holding data that identifies a
person, and declining that scope is more honest than claiming a bar this
application does not clear. Vivum reconstructs a tissue section and ranks
compounds against it; it performs no diagnosis and produces nothing that is
intended to select therapy for a patient.

Every dataset this project ships with is open, public research data (10x
Genomics/GEO Xenium and Visium sections) with no patient-identifiable
information — see `docs/DATA_INVENTORY.md`.

**If you want the commands rather than the model, read
[`RUNBOOK_DEPLOY.md`](RUNBOOK_DEPLOY.md) instead.** This page explains what each
mode *means* and what it does not do for you; the runbook is the shortest path
from a fresh clone to a finished simulation, with every command executed and
every output copied from the session that ran it. It also lists what that
exercise found broken, which on the last pass was seven things — including the
command this page led with in §1b.

## 0. Field testing on your own tissue

A field tester may point Vivum at **de-identified tissue only**. That is a real
constraint rather than a formality, because a spatial-transcriptomics export
carries identifiers in more places than the obvious columns. Before you run
`vivum ingest <dir>` or use the `/ingest` page, check all four of these:

1. **The directory path itself.** `GET /api/sections` returns `cache_dir` — the
   absolute path of the directory you registered — to every client that can
   reach the API, and `--register` writes it to
   `app/runs/_ingested_sections.json`. `/data/2026-03-14_Smith_J/` is an
   identifier, and registering it hands that identifier to everyone who can
   reach the server. Copy or rename to a study code first.
2. **The section id and the label.** You choose both (`--id` / `--label`, or the
   corresponding fields on the `/ingest` form), and both default to something
   derived from the directory name. Both are shown in the Studio section picker,
   and the label is the `<h1>` of every exported report.
3. **Sample identifiers inside the export.** A Space Ranger metrics CSV has a
   `SampleID` column holding whatever the sequencing operator typed, and 10x
   output filenames are routinely prefixed with that same string — this repo's
   own cache has `zh1007inf_metrics.csv` and
   `Xenium_FFPE_Human_Breast_Cancer_Rep1_cells.csv.gz`. Vivum does not read
   those fields, but they travel with the directory.
4. **Everything else in the directory.** Vivum reads the count matrix and the
   coordinates table and ignores the rest — but "ignores" is not "removes". A
   Visium export also ships tissue and fiducial images, QC summaries and run
   logs; whatever you put in the directory stays in the directory.

De-identified here means what the custodian of the sample would recognise as
de-identified: no name, no medical record number, no dates tied to an
individual, no accession number from your institution's systems, in any of the
four places above. If you are not the person who can make that call for a
sample, don't ingest it. This section scopes the tool; it is not a statement
about your obligations, and it does not replace whatever review your institution
requires before tissue-derived data leaves the system it was generated in.

### What the deployment model does for you, and what it leaves to you

The deployment model genuinely helps here, and it is worth stating plainly
rather than defensively: **the sample never leaves the machine you run Vivum on
unless you deliberately expose it.** Two things that were checked rather than
assumed:

* **Sections are registered by path, not by upload.** `vivum ingest` and the
  `/ingest` page both take a directory on the *server's own* filesystem; the
  bytes are read where they lie. There is no HTTP upload route for tissue at
  all, so there is no copy of the sample to lose track of.
* **The server makes no outbound network requests.** No telemetry, no crash
  reporting, no analytics, no license check anywhere in `app/server/`. The
  Studio front end loads no third-party assets either — fonts are vendored under
  `app/web/fonts/` — so a browser rendering it contacts nothing but your own
  server. (Note that `vivum acquire`, a separate command, *does* download public
  datasets; that is the one deliberate exception and it fetches, never sends.)

Which means the exposure surface is exactly the machine and the network you put
it on, and both of those are yours:

* **The default is loopback, and widening it puts tissue-derived data on the
  wire.** `vivum studio` binds `127.0.0.1`. Everything in §1b and §3 below is
  about deliberately widening that — and once you do, the browser is a remote
  client, so cell coordinates, phenotype calls and field arrays for the selected
  section travel to it over plain HTTP unless you have put TLS in front (§2).
  The raw export stays put; what it reconstructs into does not. If a field test
  doesn't need remote access, don't widen it.
* **The disk is the boundary, and the application is not guarding it.** Nothing
  is encrypted at rest. `data_cache/`, `app/runs/`, `app/reports/`,
  `app/projects/` and the ingested-section table are ordinary files with your
  umask on them. Host-level disk encryption and filesystem permissions are the
  control here, because Vivum supplies none.
* **Exports are the moment data crosses off the machine.** A downloaded report
  or a shared run carries the section label and its provenance with it. That
  happens because a person chose to; the server never sends anything on its own.
* **Removing a section is not deleting data.** `vivum ingest --remove <id>` and
  `DELETE /api/sections/<id>` drop the entry from the registry and nothing else.
  Your data directory is untouched, and so are reports already written and
  reconstructions already cached.

---

## 1. The three ways to run this

### 1a. Local, single-user (the default, unchanged)

```bash
vivum studio
```

Binds to `127.0.0.1`. No API key is required or generated. Nothing here
changes this path — every test and every piece of documentation that predates
this page still applies exactly as written.

This path is unauthenticated on purpose, and since it is easy to stop reading
here, the server now says so on every start with no key:

```
  api key     NONE -- no authentication on any /api/ route

  !!  WARNING: NO API KEY -- every /api/ route is UNAUTHENTICATED.
  ...
```

It is the right default *as long as nothing forwards traffic to that port*.
The moment something does — a tunnel, a proxy, a port forward — the loopback
bind stops meaning anything and no key is generated to compensate. Read **§1e**
before you demo this to anyone who is not sitting at the machine.

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
  tissue roots /path/to/repo/data_cache  [defaulted for this non-loopback bind]
              tissue outside that directory will be REFUSED; set TMESIM_STUDIO_DATA_ROOTS to include yours
```

Both of those lines are real output, copied from a run. Note the **second**
one: a non-loopback bind now confines the scanner as well as authenticating it
(§1c, §4). If your samples live outside this repo, that line is telling you the
server will refuse them until you say where they are.

**That command did not run until 2026-08-22, and this page led with it anyway.**
`--api-key` was an argument of `python -m app.server` only; `vivum studio` is a
thin delegate and its subparser never forwarded the flag, so the first thing
this section told a reader to type exited 2 with `unrecognized arguments:
--api-key`. It was found by running the documentation rather than reading it,
and fixed in the CLI rather than in the prose — `vivum studio` is the branded
entry point and the flag matters most on exactly the non-loopback bind this
section is about. Verified after the fix: the banner reads `[set via
--api-key]`, an unkeyed request is `401` and a keyed one is `200`.

All three spellings now work and are equivalent, except that `--api-key` beats
an exported variable (§4):

```bash
vivum studio --host 0.0.0.0 --api-key <a-real-secret>
TMESIM_STUDIO_API_KEY=<a-real-secret> vivum studio --host 0.0.0.0
python -m app.server --host 0.0.0.0 --api-key <a-real-secret>
```

The environment variable is the one that also composes with Docker's `-e` and
systemd's `Environment=`, which is why the rest of this page uses it.

**One thing to know before you rely on reading that key out of a log file.**
The startup banner goes to **stdout**, and CPython block-buffers stdout when it
is a pipe — which is what systemd, `nohup`, and any supervisor give you.
Measured on this repo, `vivum studio` with stdout redirected to a file:

| | first lines of the log |
|---|---|
| default | `INFO: Uvicorn running on ...` at line 4, `Vivum Studio backend` at line **5** |
| `PYTHONUNBUFFERED=1` | `Vivum Studio backend` at line 1, `INFO: Uvicorn running` at line **30** |

So by default the banner — including an **auto-generated API key, which is
printed exactly once and never persisted** — lands *after* the server is
already accepting requests, and until the buffer fills it is not in the file at
all. The Docker image sets `ENV PYTHONUNBUFFERED=1` for precisely this reason
and the comment in the `Dockerfile` describes the failure; nothing sets it on
the `vivum studio` path. **If you are running this under systemd or behind any
redirect, set `PYTHONUNBUFFERED=1` too**, or pass `--api-key` yourself so there
is no generated key you need to go looking for.

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

**Honesty note on this section specifically.** This section was written before
Docker was installed in the environment it was written in, and said for a long
time that no image had ever been built from this `Dockerfile`. **One has now
been built and run.** On `linux/arm64`, measured rather than reasoned about:

* `docker build` **completes**: 4 m 40 s cold, 3.9 s when the dependency layer
  is cached. (A rebuild that also re-copies the source tree costs more than the
  cached-layer figure — the build that added the non-root `USER` below took
  2 m 39 s with the dependency layer still cached.)
* The final image is **1.38 GB** — it was 1.67 GB before the `.dockerignore`
  fix below landed — of which `/app` is **10 MB**. Both figures are from
  `docker images` / `docker system df`; see the note on `docker image inspect`
  at the end of this list before quoting a size from anywhere else.
* **All 18 native dependencies import inside `python:3.11-slim` with only
  `libgomp1`.** In particular rdkit needs **no** `libXrender`/`libXext` there,
  which the Dockerfile's own comment had flagged as a historical risk.
* `/api/health` answers **200 about 2 s after `docker run`**, through a
  published port — so the `CMD` really does bind an address `-p` can reach.
  `docker stop -t 30` **exits 0 immediately**: PID 1 handles SIGTERM, there is
  no 30 s wait and no SIGKILL.
* A real end-to-end run — reconstruction plus transport plus metrics, submitted
  over the published port and polled to `done` — completes in-container on
  **40,887 cells**. On the stock `GET /api/defaults` spec (`xenium_319`,
  temozolomide, `dx 12 um`, 6 h window, `plasma_peak 50 uM`) that is
  **`runtime_s` 47.13 cold**, on a 4-CPU VM-backed daemon. The reported
  `kill_fraction` was **0.11221232991768859**, bit-identical across two
  submissions in the same container.

  This line used to read "in 4 s on 40,887 cells" and that could not be
  reproduced. The defect is not really the number: it is that **a runtime was
  quoted with no RunSpec attached**, and on this product a runtime without a
  spec is not a measurement — the window, the mesh and the compound each move
  it. Quote the spec or don't quote the time. Do not read the figure above as a
  benchmark either; it is one machine, and the point of measuring it was to
  establish that a real run finishes in a container at all.
* With **no cache mounted**, `/api/sections` reports **7 sections, 0
  available**, each with a reason, and the server still boots. Mounting
  `data_cache` read-only makes **all 7 available**.
* The server runs as **`vivum`, uid 10001, not root** (see the `USER` note
  further down). All five 0700 runtime trees are created and written by that
  user, a read-only `data_cache` bind mount stays readable, and a named volume
  on `/app/app/runs` inherits the ownership and is writable with no host-side
  `chown`.

One measurement tool to distrust: `docker image inspect --format '{{.Size}}'`
reported **300,307,441 bytes for the image `docker images` and `docker system
df` both call 1.38 GB** — an under-report of roughly 4.6×, on a
containerd-snapshotter daemon where that field is the compressed content size
rather than the on-disk size. CI reads that field. Do not quote a size from it
without checking it against `docker images` on the same daemon.

What *had* been verified earlier, by reproducing the Dockerfile's steps outside
a container in throwaway Python 3.11 virtualenvs, and which the build above is
consistent with:

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

**Still genuinely unverified — and the architecture is the reason.** Every
number above was measured on **`linux/arm64`**. The open question is whether
`manylinux` **cp311 x86_64** wheels exist and install for `scipy` /
`scikit-image` / `imagecodecs` / `h5py` / `rdkit`, and whether *those* wheels
want shared libraries this image does not carry. An arm64 build cannot answer
it: it resolves a different set of wheels from the same version pins. **Only a
CI run on `linux/amd64` settles it** — that is what
`.github/workflows/docker.yml` is for, and until it has gone green on that
platform this remains open, not "probably fine".

Three smaller things the local build also did not exercise: **SELinux
relabeling** of bind mounts (there is no SELinux on the machine this was built
on, so a `:z`/`:Z` suffix may be needed on RHEL/Fedora hosts), **BuildKit cache
reuse** (the local build used the legacy builder; CI uses BuildKit with a
`type=gha` cache), and any **host bind mount over a writable path**, which
interacts with the non-root `USER` — see the ownership note below.

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

**The container runs as `vivum`, uid 10001 — not root.** It used to run as
root, which was the exact opposite of the control §2 names for
`POST /api/sections/scan`: that route takes an arbitrary filesystem path, is
not confined to a root, and the API key makes it *authenticated, not safe*, so
§2's mitigation is "run the server as a user that cannot read anything you
would mind a tester listing". As root it could read all of it. Measured
before/after, with a valid key, against a `0700` root-owned directory holding
a valid Xenium bundle:

| | `scan {"path": "/root/probe"}` |
|---|---|
| root image | `platform: xenium`, `confidence: 1.0`, `ingestable: true`, and the absolute path of every file inside it |
| uid 10001 | `platform: unknown`, `confidence: 0.0`, `ingestable: false`, `"/root/probe is not a directory and does not exist"` |

Read what that table does and does not say. It is about the **uid**, and at the
time it was measured the uid was the only control in the image: the route was
unconfined, so `scan {"path": "/etc"}` answered `200` under both rows, because
`/etc` is world-readable and a non-root user can list it perfectly well. What
the `USER` line changed is the *set of files the process can reach*, not the
route's behaviour.

**That is no longer the whole story, and this paragraph used to say the
opposite of what the image now does.** It read: "the `path`-root restriction §2
named now exists, but it is opt-in and nothing here sets it." Something does
now. The default `CMD` binds `0.0.0.0`, and on a non-loopback bind
`app/server/__main__.py` fills in `TMESIM_STUDIO_DATA_ROOTS` with
`<repo>/data_cache` — which inside this image is `/app/data_cache`, the mount
point — exactly as it fills in an API key. So a stock `docker run` gets **both**
halves of §2's advice without being asked, and says so on startup:

```
  api key     verify-key-2026  [set via TMESIM_STUDIO_API_KEY]
  tissue roots /app/data_cache  [defaulted for this non-loopback bind]
              tissue outside that directory will be REFUSED; set TMESIM_STUDIO_DATA_ROOTS to include yours
```

Measured in that container, with a valid key:

| `scan {"path": ...}` | |
|---|---|
| `/etc` | **403** `path_outside_data_roots` — refused before the filesystem is touched |
| `/app` | **403** — the repo root is outside the allowlist, and being the *parent* of a root does not admit it |
| `/app/data_cache/../etc` | **403** — the path is resolved first, so `..` cannot walk out |
| `data_cache/<a real sample>` | **200** — resolves inside the root, and is scanned |

The cost is the one you would expect and it is the reason this is worth
knowing before a field test: **a sample that is not under the mounted
`data_cache` is refused**, including one bind-mounted somewhere else in the
container. Add it explicitly, remembering that this is `os.pathsep`-separated
and that §1d's re-register command needs `data_cache` in the list too:

```bash
docker run -p 127.0.0.1:8000:8000 \
  -v $(pwd)/data_cache:/app/data_cache:ro \
  -v /srv/tissue:/tissue:ro \
  -e TMESIM_STUDIO_DATA_ROOTS=/app/data_cache:/tissue \
  -e TMESIM_STUDIO_API_KEY=<your-key> \
  vivum-studio
```

Two things that are **not** true of this default, both checked rather than
assumed. It does not fire on a loopback bind: `docker run vivum-studio --host
127.0.0.1` still prints `api key NONE` and still answers `scan {"path":
"/etc"}` with `200`, so the §1a workflow is untouched. And it does not
overwrite a value you set — including an empty one. `-e
TMESIM_STUDIO_DATA_ROOTS=` (no value) is read as "you asked for confinement and
named no directory", and the server refuses **every** scan rather than falling
back to the default root or to unconfined:

```
  tissue roots NONE -- every scan is refused  [set via TMESIM_STUDIO_DATA_ROOTS]
```

A scan of `/app/data_cache` — the directory the automatic default *would* have
allowed — then returns `403`. That is the intended behaviour (§4) and it is
what makes the auto-default safe to have: it can only ever add confinement to a
deployment that specified none.

**What this costs you, and it is only the one thing:** a **host directory**
bind-mounted over a writable path arrives with the *host's* ownership, and uid
10001 almost certainly does not own it.

* **Use a named volume, and there is nothing to do** — `-v
  vivum_runs:/app/app/runs`, as in the example above. Docker seeds a new named
  volume from the image's directory, ownership included, so it comes up owned
  by uid 10001 and writable with no host-side step. **Verified in a real
  container**: a full run wrote its history into exactly that mount. This is
  the recommended shape, and it is the reason the `chown` in the `Dockerfile`
  is worth having.
* **A host bind mount over `app/runs/`** (or `app/projects/`, `app/reports/`,
  `app/compounds/`, `app/uploads/`) **will fail without a host-side step.**
  Measured, mounting an ordinary host directory:

  ```
  /app/app/runs is owned by uid 0 inside the container; the server runs as 10001
  touch: cannot touch '/app/app/runs/probe': Permission denied
  ```

  On a Linux daemon, where host and container uids are the same numbers, the
  fix is `sudo chown -R 10001:10001 <that directory>` on the host before
  mounting it. **That fix was not verified here** — it could not be: this was
  measured on a VM-backed daemon (colima/Lima on macOS), where the mount is
  translated and the host directory showed up as `uid 0` regardless of who owns
  it on the host, so no host-side `chown` can produce uid 10001 through it.
  **On that kind of daemon, use a named volume; a host bind mount over a
  writable path is not a workable shape.** If you need a host bind mount on a
  Linux daemon, do the `chown` and check it before you rely on it:

  ```bash
  docker run --rm -v "$PWD/runs:/app/app/runs" --entrypoint sh vivum-studio \
    -c 'touch /app/app/runs/.probe && echo WRITABLE || echo NOT-WRITABLE'
  ```

**`docker run --user "$(id -u):$(id -g)"` is not a drop-in substitute for that
`chown`** — this was tried and measured, and it fails:

```
PermissionError: [Errno 13] Permission denied: '/app/app/runs/.probe'
PermissionError: [Errno 13] Permission denied: '/app/app/reports'
```

`/app/app` is owned by uid 10001, and four of the five runtime trees do not
exist in the image at all — they are created on first use, *inside* that
directory. An arbitrary overriding uid can neither write the existing ones nor
create the missing ones. Override the user only if you are prepared to mount
every one of the five writable paths from the host with ownership to match —
which is strictly more work than the named volume.

`data_cache/` is unaffected: it is mounted **read-only** in every example here,
and a read-only mount only needs to be *readable*, which an ordinary
`0644`-inside-`0755` host cache is. Verified: with `-v
$(pwd)/data_cache:/app/data_cache:ro` the non-root container reports all 7
sections available.

**Check the host port is free before you publish onto it, because Docker may
not tell you.** Measured: with a plain `python -m http.server` already holding
`127.0.0.1:8305`, `docker run -d -p 127.0.0.1:8305:8000 vivum-studio` **exited
0, the container stayed `running`, and every request to that port reached the
*host* process instead.** `/api/health` answered — with the host process's
JSON. Nothing in `docker ps`, `docker logs` or the response says the container
is not the thing you are talking to; you can spend a long time reading a stale
server's data and calling it a Vivum bug.

That was measured on a VM-backed daemon (colima/Lima on macOS), where the
published port is forwarded by a host-side proxy rather than bound by the
daemon; a daemon that binds the port itself will usually refuse with `bind:
address already in use` instead. Do not assume which kind you have. The
pre-check costs one command and is correct on both:

```bash
lsof -nP -iTCP:8000 -sTCP:LISTEN     # must print nothing before you publish :8000
```

Empty output means the port is free. Anything else — an old `vivum studio`, a
previous container, an unrelated dev server — take the other port, or stop that
process first.

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

## 1e. Tunnels and reverse proxies — where the loopback default stops working

**Read this before you demo Vivum to anyone remote.** It is the most likely way
a field test goes wrong, and nothing in the software can detect it for you.

The auto-generate rule in §1b keys off the **bind address**. The thing that
actually matters is **reachability**. Those are the same right up until you do
the single most ordinary thing anyone does when showing software to a colleague:

```bash
vivum studio                 # binds 127.0.0.1, so no key is generated
ngrok http 8000              # ...and now it is on the internet
```

`ngrok`, `cloudflared`, `ssh -R`, VS Code / JetBrains port forwarding, and any
reverse proxy that does not authenticate all do this. (Docker is *not* on this
list: the image's default `CMD` binds `0.0.0.0`, so §1b's rule fires and a key is
generated — see §1c.) In every case on the list the server still sees
`--host 127.0.0.1`, `is_loopback` is still true, **no key is generated**, and a
public URL now reaches an unauthenticated API. The startup banner warns on every
keyless start (§1a) precisely because that is the only place the risk can be
surfaced — the process cannot see the tunnel.

What is on the other end of that URL is worse than "your tissue data":

* `POST /api/sections/scan` takes an **arbitrary filesystem path in the request
  body** and reports on any directory the server process can read. Unauthenticated
  and public, that is a filesystem reconnaissance primitive on the host, not just
  access to the sections you meant to share. Verified by request, not inferred.
  `TMESIM_STUDIO_DATA_ROOTS` (§4) bounds which directories it will read, and is
  worth setting here as well as a key. **The automatic confinement does not
  help you here**: like the automatic key, it keys off the bind address, and
  the whole point of this section is that you are on a loopback bind. Both
  defaults are off, and a tunnel makes neither of them fire. Set both by hand.
* The destructive routes are equally open: `DELETE /api/sections/{id}` unregisters
  an ingested section, `DELETE /api/compounds/{id}` forgets a stored compound, and
  `DELETE /api/projects/{id}` removes a project. `DELETE /api/runs/{id}` and
  `DELETE /api/reports/{id}` cancel work in flight. Anyone who finds the URL can do
  all of it. (Per §0, unregistering a section does not touch the data directory —
  but it does mean someone else can empty your section list and your compound
  library while you are demoing.)

**So: if anything forwards traffic to Vivum, set a key, regardless of what you
bound to.**

```bash
TMESIM_STUDIO_API_KEY=<a-real-secret> vivum studio   # still 127.0.0.1, now gated
```

This is also why §3's reverse-proxy setup says the proxy is the boundary: that
is only true if the proxy **actually authenticates**. A proxy that terminates TLS
and forwards everything is not a boundary — it is a tunnel with a certificate.
If the thing in front of Vivum does not check who the caller is, Vivum must.

## 2. What is NOT handled for you

* **TLS.** This serves plain HTTP. Put a reverse proxy in front for anything
  crossing an untrusted network (Caddy and nginx both terminate TLS for a
  single backend in ~10 lines of config; this project deliberately doesn't
  bundle one, to avoid maintaining certificate/domain logic orthogonal to the
  simulator itself).
* **Multi-user accounts.** One shared key, not per-user identity. See §1b.
* **Confining `POST /api/sections/scan` to a directory root — on a loopback
  bind, unless you turn it on.** The key makes `scan` **authenticated, not
  safe**. It takes an arbitrary `path` in the request body, so with no
  allowlisted root any tester holding the key can enumerate and read any
  directory the *server process* can read, anywhere on the host, not just under
  `data_cache/`. The allowlist that stops that is `TMESIM_STUDIO_DATA_ROOTS`,
  §4.

  **Whether it is on by default depends entirely on the bind address, so work
  out which case you are in.** On a **non-loopback** bind — `--host 0.0.0.0`,
  which is every stock `docker run` — it is filled in with `<repo>/data_cache`
  automatically, alongside the auto-generated key, and the startup banner names
  the directory it chose. On a **loopback** bind it is **off**, deliberately:
  bring-your-own tissue is the product's main input, §0 tells you to copy a
  sample to a study-code directory that is almost never inside this repo, and
  confining the single-user workflow by default would break it silently. The
  first unconfined scan in each process prints a warning naming this variable.

  **So on a loopback bind, until you set it, treat the key as equivalent to
  shell read access on the host and hand it out on that basis** — and note from
  §1e that a loopback bind is exactly the case where a tunnel can publish that
  API with no key generated either. Set it before you give the key to anyone
  who is not you. The other control, which applies either way, is §0's: run the
  server as a user that cannot read anything you would mind a tester listing.

  **The Docker image now does that for you** — it runs as `vivum`, uid 10001,
  after a spell of running as root, i.e. of doing the precise opposite of the
  sentence above. §1c has the measured before/after and the one thing it costs
  (a `chown` on host bind mounts over writable paths). Running `vivum studio`
  directly still runs as **whoever you are**, so on that path the sentence is
  still advice you have to take yourself: don't run it as root, and don't run
  it as the account that can read the rest of the department's data.
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

**The load-bearing word above is *authentication*, not *proxy*.** "No API key
needed at all" holds only while the proxy refuses unauthenticated callers. Drop
the basic-auth block, point the VPN-only hostname at the public internet, or
swap Caddy for a plain tunnel, and this diagram becomes §1e: a public URL onto an
unauthenticated API, with the server still reporting a loopback bind and still
generating no key. If you are not certain the front end authenticates, set
`TMESIM_STUDIO_API_KEY` as well — two checks cost nothing.

Three more things this setup needs that the loopback default does *not* hand
you, all covered in §4: set `TMESIM_STUDIO_ALLOWED_HOSTS` to the hostname the
proxy forwards (the automatic loopback allowlist is `localhost,127.0.0.1`, which
will **400 every proxied request** arriving with your real hostname unless the
proxy rewrites `Host`), set `TMESIM_STUDIO_DATA_ROOTS` to the directories your
samples live in (everyone the proxy admits can otherwise make this server read
any directory it can read), and leave `TMESIM_STUDIO_EXPOSE_TRACEBACKS` unset.

---

## 4. Environment variables

Ten. There is no config file; these are the whole surface, and only five are
security-relevant.

This section said "Six" while the code read ten — `TMESIM_MAX_GRID_CELLS` had
been load-bearing for a while and was never listed here, which is the failure
this repository keeps having: the document that was true when it was written.
`tests/test_app_server.py::TestEnvVarsAreDocumented` now enumerates every
`TMESIM_*` read anywhere in `app/server/` and fails if one is missing from this
file, so the count above cannot drift again without a red test.

**Four of the six are read once, at process start, and need a restart to
change.** `TMESIM_STUDIO_API_KEY`, `TMESIM_STUDIO_ALLOWED_HOSTS` and
`TMESIM_STUDIO_ALLOWED_ORIGINS` are read while the app is being built — each
decides whether a middleware is installed at all, so there is nothing left to
re-read afterwards — and `TMESIM_RECON_CACHE_SIZE` is read at import, when the
cache object is constructed.

**`TMESIM_STUDIO_EXPOSE_TRACEBACKS` and `TMESIM_STUDIO_DATA_ROOTS` are the
exceptions: they are read per request**, inside the
`GET /api/runs/{id}/traceback` and `POST /api/sections/scan` handlers
(`app/server/api.py`). Practically that means the value in force is whatever is
in the *live process's* environment at the moment the route is called, not
whatever was set when it started — so a debugging session that mutates the
process environment takes effect immediately, and, less comfortably, so would
anything else that mutates it. It also means there is no restart-shaped moment
at which an operator gets to notice that traces are being served, or that a
scan root was widened.

`vivum studio` / `python -m app.server` writes **three** of them into its own
environment before the app is built — `TMESIM_STUDIO_API_KEY`,
`TMESIM_STUDIO_ALLOWED_HOSTS` and, since the confinement default landed,
`TMESIM_STUDIO_DATA_ROOTS` — and in all three cases a value you exported
yourself **beats the automatic behaviour and is never clobbered**, so `-e` on
`docker run`, `Environment=` in a systemd unit and a plain shell export all
compose correctly. "Never clobbered" includes an **empty** value: exporting
`TMESIM_STUDIO_DATA_ROOTS=` is read as a confinement you asked for and left
empty, and every scan is refused rather than the default root being filled in
behind you (verified; §1c). An explicit `--api-key` on the command line still
overrides the exported key — that is the one deliberate exception, it is the
ordering you'd expect, and it is available only on `python -m app.server`
(§1b).

### `TMESIM_MAX_GRID_CELLS`

**The one that keeps a small host alive**, and the reason this section was
audited at all.

Measured on a 512 MB instance (Render's free tier), flagship Xenium GBM, real
container: `dx = 30 um` is 46,070 voxels and peaks at 220 MiB — comfortable.
`dx = 12 um` is 286,200 voxels and the container is **OOM-KILLED**. Not the run
failing: the *process* dies, `/healthz` stops answering, and one ordinary
request takes the service down for everyone. The default ceiling of 40,000,000
does not catch that, because 286,200 voxels is a perfectly reasonable ask on a
machine with room for it.

* **Unset:** 40,000,000. That is the "no request should ever ask for this"
  number, not the "this machine can afford it" one. Those are different numbers
  and only the operator knows the second.
* **An integer:** used verbatim. `100000` permits `dx >= ~20 um` on the
  flagship section and refuses `dx = 12`.
* **`auto`:** adopts 46,070 — the largest grid *measured* to survive a host this
  size — but **only when a limit at or below 512 MB is actually detected**. On a
  larger or undeclared host it falls back to the default and says which branch
  ran, so `auto` never hides its own inaction. It applies a measurement rather
  than fitting one: the only two points that exist are "46,070 survived" and
  "286,200 died", so it is deliberately conservative about everything between
  them. An operator who wants that range back sets a number.
* **Anything else:** the process **refuses to start**. A typo must not fall back
  to 40,000,000 on a small box.

Detection alone never moves the ceiling — that is opt-in, by name. What
detection always does is *warn*: when the configured ceiling still permits the
286,200-voxel spec on a host at or below 512 MB, the startup banner says so and
so does the Account panel's Limits tab. On a correctly configured host both are
silent.

On Render's free tier, set `TMESIM_MAX_GRID_CELLS=auto` or `100000`.

### `TMESIM_HOST_MEMORY_MB`

An override for the memory detection above, in megabytes. Read only to answer
"how big is this host", never to set a ceiling by itself.

Detection reads cgroup v2 `memory.max`, then cgroup v1
`memory.limit_in_bytes`, and treats `max` and v1's huge sentinel as *no limit
declared*. Set this when a platform enforces a limit it does not expose through
cgroups. Undetected memory is reported as **not declared**, never as
"unlimited" — a host can be capped by something this process cannot see.

### `TMESIM_BILLING_PROVIDER`

Names a payment provider for `GET /api/account`'s billing block. **Setting it
does not enable billing and cannot**: with no implementation mounted the block
reports `configured-but-not-implemented` and `enabled` stays `false`. It exists
so that wiring a provider later is a config change plus a real implementation,
and so a half-finished integration is a state an operator can *see* rather than
a silent one. No charge is made, no payment method is stored, no usage is
metered.

### `TMESIM_STUDIO_API_KEY_SOURCE`

**Written by the process, not for you to set.** `python -m app.server` records
which of its three key sources won — `generated`, `flag` or `environment` — so
`GET /api/account` can report it weeks later, when the startup banner that said
the same thing is scrollback nobody has. An operator who has *lost* the key
needs to know whether they chose it or the process invented one at boot and
printed it exactly once. An unrecognised value is reported as `unknown` rather
than mapped to a plausible guess.

### `TMESIM_STUDIO_API_KEY`

The shared secret checked against the `X-API-Key` header on every `/api/*` route,
and on every non-`GET`/`HEAD`/`OPTIONS` request whatever its path. Static front-end
files are deliberately *not* gated, so the page can load and prompt for the key
rather than failing to render.

* **Unset or empty:** the gate is not installed at all. Every route is open. On a
  loopback bind this is the intended §1a workflow; the startup banner now prints
  an unmissable warning either way.
* **Set to a non-empty value:** required on every API request; a wrong or missing
  header is a `401`. Compared with `secrets.compare_digest`, so it does not leak
  its prefix by timing. Leading/trailing whitespace is stripped before the compare.
* **Must set it when:** anything other than you on that machine can reach the port
  — including through a tunnel or proxy (§1e). A non-loopback `--host` with no key
  generates one for you and prints it once; a loopback bind behind a tunnel does
  not, which is exactly the gap this variable closes by hand.
* `--api-key` sets the same thing and takes precedence over the exported value.
  `--api-key ''` explicitly disables the gate even on a non-loopback bind; the
  warning still prints, which is the point of it. `vivum studio` forwards the
  flag to `python -m app.server` and preserves all three states, including the
  `--api-key ''` one — it did not accept the flag at all until 2026-08-22, when
  running §1b's own headline command found that out (§1b).

### `TMESIM_STUDIO_ALLOWED_HOSTS`

Comma-separated hostnames for Starlette's `TrustedHostMiddleware`. This is
defence in depth against DNS rebinding — a hostile domain that resolves to
`127.0.0.1` becomes same-origin with your server, which CORS cannot stop. It is
**not** the fix for a Host-header auth bypass; that is fixed at the gate itself.

* **Unset:** the middleware is not installed and there is no host validation at
  all — *except* that on an **IPv4 loopback** bind `vivum studio` fills it in with
  `localhost,127.0.0.1` before the app is built, so the default local workflow is
  covered without you doing anything. **`--host ::1` is deliberately excluded
  from that fill-in and therefore gets no validation** — see the sharp edge at
  the end of this entry. On a **non-loopback** bind it is left unset
  on purpose: the code does not know the hostname your deployment answers to, and
  guessing wrong would refuse every request. **A proxied or `--host 0.0.0.0`
  deployment therefore gets no host validation until you set this.**
* **Exported as an empty string:** that counts as "you set it", so the automatic
  loopback fill-in is skipped and you get no validation. Unset it rather than
  blanking it.
* **Set:** only listed names are accepted; anything else is a `400 Invalid host
  header`. A literal `*` in the list disables the middleware entirely rather than
  matching everything — same effect, but worth knowing.
* **Must set it when:** you bind non-loopback or sit behind a proxy. Use the
  hostname the *browser* sends, which is the one the proxy forwards.
* **Sharp edge, and it is a real one: an IPv6 loopback bind gets no host
  validation at all.** `vivum studio --host ::1` is *excluded* from the
  automatic loopback fill-in above, deliberately, so `TrustedHostMiddleware` is
  never installed on it. The server works — `http://[::1]:<port>/` returns
  **200** — but it has **no DNS-rebinding protection**, unlike the IPv4
  loopback default. Measured, same build, same request:

  | bind | `Host: evil.example.com` |
  |---|---|
  | `--host 127.0.0.1` | **400** `Invalid host header` |
  | `--host ::1` | **200** |

  **Why it is excluded rather than fixed here.** Starlette matches on
  `host.split(":")[0]` to strip the port. That is right for a name or an IPv4
  literal and wrong for an IPv6 one: a browser sends `Host: [::1]:8000`, which
  splits to `[`. No allowlist entry can ever equal `[`, so installing the
  middleware on an `::1` bind would reject **100 % of its own traffic** with
  `400 Invalid host header` — a supported flag broken by its own protection.
  (That is not hypothetical: it is what this paragraph used to describe as
  verified behaviour, before the exclusion landed.) A middleware that cannot
  express the bind buys no safety and costs the whole server, so IPv6 loopback
  gets the same treatment as a non-loopback bind: unset, and yours to set.

  **So on `::1`, set it yourself, and set it to a NAME:** a literal still cannot
  match, for the same `split(":")` reason.

  ```bash
  TMESIM_STUDIO_ALLOWED_HOSTS=localhost vivum studio --host ::1   # reach it at http://localhost:8000/
  ```

  Measured on that command: `http://localhost:<port>/` → **200**,
  `Host: evil.example.com` → **400** (the protection is now on), and
  `http://[::1]:<port>/` → **400**. That last one is the price, not a bug —
  turning the check on is exactly what stops the bracketed literal working.
  Reach it by the name from then on.

  If you don't need IPv6, **`--host 127.0.0.1` is the safer default** — not
  because `::1` is broken, but because it is the one that comes with the
  rebinding check already on.

### `TMESIM_STUDIO_ALLOWED_ORIGINS`

Comma-separated origins for CORS. **CORS is off by default now.** It used to be
`allow_origins=["*"]`, which meant any web page you happened to visit could read
*and delete* your whole compound library through `127.0.0.1` — binding loopback is
no defence there, because the browser is a local process.

* **Unset:** no `CORSMiddleware` is installed and no cross-origin grant exists.
  This is correct for every normal deployment, because the front end is served
  same-origin from this same app.
* **Set:** exactly the origins listed are allowed, for `GET`/`POST`/`DELETE`/
  `OPTIONS` with the `content-type` and `x-api-key` headers. **`*` is filtered out
  and ignored, not honoured** — an entry of `*` contributes nothing, and a value
  of just `*` leaves CORS off entirely.
* **Must set it when:** you run the split dev-port workflow, with a dev front end
  on its own port talking to this backend. That is broken until you set it:

  ```bash
  TMESIM_STUDIO_ALLOWED_ORIGINS=http://localhost:5173 vivum studio
  ```

  Full origins including scheme and port, not bare hostnames.

### `TMESIM_STUDIO_DATA_ROOTS`

Directories `POST /api/sections/scan` is allowed to read, separated by `:`
(`os.pathsep`) — **not** commas like the two list-valued variables above,
because a directory name may legitimately contain a comma. This is the
confinement §2 says the API key does not give you: the key makes `scan`
*authenticated*, this makes it *bounded*.

It applies to the **HTTP route** — which is the trust boundary, and also what
the `/ingest` page calls. `vivum ingest <dir>` on the command line is not
confined by it and is not meant to be: that is the operator running a program as
themselves on their own machine, with no boundary to enforce. So the two
surfaces can disagree about a directory, and when they do it is because one of
them is being driven by a remote caller.

The check is on the **resolved** path. A leading `~`, every `..` component and
every symlink are followed first, and it is the directory that path really
designates which must be inside a root — so `data_cache/../..`, a symlink
planted inside an allowed root, and a `..` that traverses *through* a symlink
are all refused. Containment is by path component, so a root of `/data` allows
`/data` itself and `/data/anything`, and does **not** allow `/data_other`.

* **Unset on a NON-LOOPBACK bind:** `vivum studio` / `python -m app.server`
  fills it in with `<repo>/data_cache` before the app is built and prints the
  directory it chose, on the same principle as the auto-generated API key: the
  moment the server is reachable from another machine, an unconfined
  filesystem-scan route is not a default anyone chose. **If your samples live
  outside this repo you must set it yourself** — the fill-in is a floor for a
  deployment that specified nothing, not a guess at where your tissue is, and
  it will `403` a directory you legitimately meant to scan. Inside the Docker
  image this resolves to `/app/data_cache`, the documented mount point (§1c).
* **Unset on a LOOPBACK bind:** no confinement. `scan` reads any directory the
  server process can read, exactly as it always has, and **the first scan in
  each process prints a warning on stderr** naming this variable and the path
  it just read. This is
  the same trade as §1a's unauthenticated loopback default: bring-your-own
  tissue is the product's main input and the documented workflow points the
  scanner at wherever the operator's export actually lives (§0 tells you to copy
  a sample to a study-code directory, which is almost never inside this repo), so
  confining by default would silently break the single-user workflow this tool
  was built for. The default therefore stays permissive and says so out loud.
* **Set to one or more directories:** only paths resolving inside one of them
  are scanned. Anything else is a `403` whose body names *your* path, says it
  resolved outside the allowlist, names this variable and says how to add a
  root — and never discloses the roots themselves. The server's stderr gets the
  resolved path; the client does not (it would otherwise leak the operator's
  home directory to anyone who sent `~/x`). One thing a caller *can* learn, and
  it is unavoidable in a diagnostic that tells an operator what to do: one bit
  per request about whether a given directory is inside the allowlist.
* **Set to a value that names no usable directory** — empty, whitespace, a bare
  `:`, as an `-e TMESIM_STUDIO_DATA_ROOTS=` with no value produces — **refuses
  every scan** rather than falling back to unconfined. You asked for
  confinement; a typo must not quietly hand back the behaviour you were turning
  off. Unset it to go back to unconfined.
* **Must set it when:** anyone who is not you holds the API key — an external
  field tester, a shared team deployment, anything behind §3's proxy. Without
  it, §2's sentence stands: the key is equivalent to read access on that host as
  the server's uid. Point it at the directory tree your samples live in:

  ```bash
  TMESIM_STUDIO_DATA_ROOTS=/srv/tissue:/data/study_a vivum studio --host 0.0.0.0
  ```

  **Include this repo's own `data_cache/` if you use it.** §1d's re-register
  command passes `data_cache/<your-section-dir>`, and a roots list that names
  only your tissue tree will `403` it. Two entries, not one:

  ```bash
  TMESIM_STUDIO_DATA_ROOTS=/srv/tissue:$PWD/data_cache vivum studio --host 0.0.0.0
  ```

* **Three things it does not do**, stated so nobody reads more into it. It
  bounds the *directory* that is scanned, not the files inside it — a symlink
  *within* an allowed sample directory still points wherever it points, and the
  ingester follows it. It is a check at one instant, so a symlink re-pointed
  between the check and the read is not caught. And it is not retroactive:
  sections already in `app/runs/_ingested_sections.json` keep working and keep
  reporting their absolute `cache_dir` through `GET /api/sections` even if that
  directory is outside every root — setting this variable stops new directories
  being reached, it does not un-register old ones (`vivum ingest --remove <id>`
  does). All three are far smaller than the "any path at all" they replace; none
  is closable at this layer. The control that covers them is still §0's: run the
  server as a user that cannot read anything you would mind a tester listing.
* **One behavioural difference when it is set:** the scanner is handed the
  *resolved* path, so that the directory checked and the directory opened cannot
  be two different things. A section registered through a symlink therefore
  records its real directory in `cache_dir` (and, if you do not pass
  `section_id`, derives its id from the real directory's name). Unset, the path
  you sent is passed through untouched, exactly as before.
* Read per request, not at startup — but the whole allowlist is process
  environment, so in practice it changes when you restart.

### `TMESIM_STUDIO_EXPOSE_TRACEBACKS`

Whether `GET /api/runs/{id}/traceback` returns the Python stack for a failed run.

* **Unset (or any value other than `1`):** the stack is **not served**. The
  response carries the typed diagnosis in `error` — which is what the product
  actually promises and is always returned — plus a `traceback_withheld`
  explanation. The unscrubbed stack always goes to the server's stderr for whoever
  is at the console.
* **Set to exactly `1`:** the stack is returned, with the repo root replaced by
  `<repo>`. Only that literal value enables it; `true`/`yes` do not.
* **Must set it when:** you are debugging a failure locally and the typed error is
  not enough. **Leave it unset on any shared or proxied deployment** —
  `traceback.format_exc()` carries absolute host paths, which disclose the repo
  root, the operator's OS username and the internal module layout to any client
  that can reach the route.

### `TMESIM_RECON_CACHE_SIZE`

How many reconstructions the in-memory LRU holds. **Default `4`** when unset. Each
entry holds a full cell population and several float64 grids, so this is the main
memory dial on the process — raise it on a machine with headroom and a team
switching between sections, lower it to `1` on a small container. Equivalent to
`--cache-size`, which sets this variable before the cache is imported. Read once,
at import; changing it needs a restart. Not a security control.

* **Valid values are integers `>= 1`.** Whitespace is stripped, so `" 8 "` is 8.
* **Unset:** `4`, silently. That is the documented default, not a mistake.
* **Anything else — empty, non-numeric, `0`, negative —** is refused with a
  diagnostic on **stderr** naming the variable, the value it was given and the
  valid range, and the process **continues on the default `4`**. The server's
  own `recon cache maxsize=` banner line then shows the value that really took
  effect, so the two never disagree.

  ```
  [vivum] TMESIM_RECON_CACHE_SIZE='' is not usable: it is empty. Valid values are
  integers >= 1 (each entry holds a full cell population and several float64
  grids; 1 is the minimum, 4 is the default). Continuing with 4. ...
  ```

  Both halves of that were defects until recently, in opposite directions.
  **Empty or non-numeric killed the process at import** with a bare
  `ValueError: invalid literal for int() with base 10: ''` and a traceback
  pointing at `app/server/cache.py` — no mention of the variable or of what
  would have been acceptable. `-e TMESIM_RECON_CACHE_SIZE=` on a `docker run`
  is enough to produce it, and in Docker that is a **restart loop** whose only
  symptom is a stack trace from a cache module. Meanwhile **`0` and `-1` were
  accepted silently** and are not capacities at all: the eviction loop drops
  the entry it has just inserted, so every run rebuilds a reconstruction that
  costs 1 s to ~90 s while `/api/health` reports a perfectly healthy
  `maxsize=0`.
* **Why it diagnoses and continues rather than exiting.** This is the one
  variable here where continuing is clearly right: it changes no result, no
  protocol number and no access decision — only how much work the process
  re-uses. Exiting would also convert a cosmetic misconfiguration into an
  unimportable package for the CLI and the test suite, which import the same
  module. A variable whose bad value *could* change a result or open access
  should exit instead; this one cannot.

*(One further variable, `TMESIM_PROTOCOL_MAX_WINDOW_S`, is read by the engine
rather than the server — it caps how long a transport solve may cover, default
24 h — and the value in force is reported in every result. It is documented with
the protocol language, not here.)*
