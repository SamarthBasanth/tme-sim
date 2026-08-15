# Deployment — running Studio for a team, not just yourself

Studio (`app/server` + `app/web`) was built as a local research workstation
tool: `python -m app.server` binds to `127.0.0.1` by default and, until this
document, had no authentication at all. This page covers what changes when a
pharma/biotech or research team wants to reach it from somewhere other than
the machine it's running on.

**What this does NOT cover: identifiable patient data, or a clinical
decision.** Codon is not built for either and must not be used with either. It
has one shared secret and no accounts, no audit log, no encryption at rest, and
no per-record access control — §2 is the full list of what it does not do for
you. Software with those gaps has no business holding data that identifies a
person, and declining that scope is more honest than claiming a bar this
application does not clear. Codon reconstructs a tissue section and ranks
compounds against it; it performs no diagnosis and produces nothing that is
intended to select therapy for a patient.

Every dataset this project ships with is open, public research data (10x
Genomics/GEO Xenium and Visium sections) with no patient-identifiable
information — see `docs/DATA_INVENTORY.md`.

## 0. Field testing on your own tissue

A field tester may point Codon at **de-identified tissue only**. That is a real
constraint rather than a formality, because a spatial-transcriptomics export
carries identifiers in more places than the obvious columns. Before you run
`codon ingest <dir>` or use the `/ingest` page, check all four of these:

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
   `Xenium_FFPE_Human_Breast_Cancer_Rep1_cells.csv.gz`. Codon does not read
   those fields, but they travel with the directory.
4. **Everything else in the directory.** Codon reads the count matrix and the
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
rather than defensively: **the sample never leaves the machine you run Codon on
unless you deliberately expose it.** Two things that were checked rather than
assumed:

* **Sections are registered by path, not by upload.** `codon ingest` and the
  `/ingest` page both take a directory on the *server's own* filesystem; the
  bytes are read where they lie. There is no HTTP upload route for tissue at
  all, so there is no copy of the sample to lose track of.
* **The server makes no outbound network requests.** No telemetry, no crash
  reporting, no analytics, no license check anywhere in `app/server/`. The
  Studio front end loads no third-party assets either — fonts are vendored under
  `app/web/fonts/` — so a browser rendering it contacts nothing but your own
  server. (Note that `codon acquire`, a separate command, *does* download public
  datasets; that is the one deliberate exception and it fetches, never sends.)

Which means the exposure surface is exactly the machine and the network you put
it on, and both of those are yours:

* **The default is loopback, and widening it puts tissue-derived data on the
  wire.** `codon studio` binds `127.0.0.1`. Everything in §1b and §3 below is
  about deliberately widening that — and once you do, the browser is a remote
  client, so cell coordinates, phenotype calls and field arrays for the selected
  section travel to it over plain HTTP unless you have put TLS in front (§2).
  The raw export stays put; what it reconstructs into does not. If a field test
  doesn't need remote access, don't widen it.
* **The disk is the boundary, and the application is not guarding it.** Nothing
  is encrypted at rest. `data_cache/`, `app/runs/`, `app/reports/`,
  `app/projects/` and the ingested-section table are ordinary files with your
  umask on them. Host-level disk encryption and filesystem permissions are the
  control here, because Codon supplies none.
* **Exports are the moment data crosses off the machine.** A downloaded report
  or a shared run carries the section label and its provenance with it. That
  happens because a person chose to; the server never sends anything on its own.
* **Removing a section is not deleting data.** `codon ingest --remove <id>` and
  `DELETE /api/sections/<id>` drop the entry from the registry and nothing else.
  Your data directory is untouched, and so are reports already written and
  reconstructions already cached.

---

## 1. The three ways to run this

### 1a. Local, single-user (the default, unchanged)

```bash
codon studio
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
codon studio --host 0.0.0.0 --api-key <a-real-secret>
```

or generate one automatically:

```bash
codon studio --host 0.0.0.0
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
TMESIM_STUDIO_API_KEY=<your-key> codon studio --host 0.0.0.0
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
docker build -t codon-studio .
docker run -p 8000:8000 \
  -v $(pwd)/data_cache:/app/data_cache \
  -e TMESIM_STUDIO_API_KEY=<your-key> \
  codon-studio --host 0.0.0.0
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
* A real end-to-end run completed **in-container in 4 s on 40,887 cells**.
* With **no cache mounted**, `/api/sections` reports **7 sections, 0
  available**, each with a reason, and the server still boots. Mounting
  `data_cache` read-only makes **all 7 available**.
* The server runs as **`codon`, uid 10001, not root** (see the `USER` note
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
`docker run codon-studio codon acquire` does **not** work (it becomes
`python -m app.server codon acquire` and exits 2). Use:

```bash
docker exec <container> codon acquire
```

**Writable state does not survive container replacement.** Run history,
`app/runs/_ingested_sections.json` (the file §1d tells you to populate),
`app/projects/`, `app/reports/` and `app/uploads/` all live inside the
image's filesystem. Mount them if you want them to persist:

```bash
docker run -p 127.0.0.1:8000:8000 \
  -v $(pwd)/data_cache:/app/data_cache \
  -v codon_runs:/app/app/runs \
  codon-studio
```

**The container runs as `codon`, uid 10001 — not root.** It used to run as
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

Read what that table does and does not say. **The route is still unconfined and
still answers `200`** — `scan {"path": "/etc"}` returns `200` under both,
because `/etc` is world-readable and a non-root user can list it perfectly
well. What changed is the *set of files the process can reach*, not the
route's behaviour. §2's deferred `path`-root restriction is still deferred, and
the key is still equivalent to read access as uid 10001 on that host.

**What this costs you, and it is only the one thing:** a **host directory**
bind-mounted over a writable path arrives with the *host's* ownership, and uid
10001 almost certainly does not own it.

* **Use a named volume, and there is nothing to do** — `-v
  codon_runs:/app/app/runs`, as in the example above. Docker seeds a new named
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
  docker run --rm -v "$PWD/runs:/app/app/runs" --entrypoint sh codon-studio \
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
`127.0.0.1:8305`, `docker run -d -p 127.0.0.1:8305:8000 codon-studio` **exited
0, the container stayed `running`, and every request to that port reached the
*host* process instead.** `/api/health` answered — with the host process's
JSON. Nothing in `docker ps`, `docker logs` or the response says the container
is not the thing you are talking to; you can spend a long time reading a stale
server's data and calling it a Codon bug.

That was measured on a VM-backed daemon (colima/Lima on macOS), where the
published port is forwarded by a host-side proxy rather than bound by the
daemon; a daemon that binds the port itself will usually refuse with `bind:
address already in use` instead. Do not assume which kind you have. The
pre-check costs one command and is correct on both:

```bash
lsof -nP -iTCP:8000 -sTCP:LISTEN     # must print nothing before you publish :8000
```

Empty output means the port is free. Anything else — an old `codon studio`, a
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

**Read this before you demo Codon to anyone remote.** It is the most likely way
a field test goes wrong, and nothing in the software can detect it for you.

The auto-generate rule in §1b keys off the **bind address**. The thing that
actually matters is **reachability**. Those are the same right up until you do
the single most ordinary thing anyone does when showing software to a colleague:

```bash
codon studio                 # binds 127.0.0.1, so no key is generated
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
* The destructive routes are equally open: `DELETE /api/sections/{id}` unregisters
  an ingested section, `DELETE /api/compounds/{id}` forgets a stored compound, and
  `DELETE /api/projects/{id}` removes a project. `DELETE /api/runs/{id}` and
  `DELETE /api/reports/{id}` cancel work in flight. Anyone who finds the URL can do
  all of it. (Per §0, unregistering a section does not touch the data directory —
  but it does mean someone else can empty your section list and your compound
  library while you are demoing.)

**So: if anything forwards traffic to Codon, set a key, regardless of what you
bound to.**

```bash
TMESIM_STUDIO_API_KEY=<a-real-secret> codon studio   # still 127.0.0.1, now gated
```

This is also why §3's reverse-proxy setup says the proxy is the boundary: that
is only true if the proxy **actually authenticates**. A proxy that terminates TLS
and forwards everything is not a boundary — it is a tunnel with a certificate.
If the thing in front of Codon does not check who the caller is, Codon must.

## 2. What is NOT handled for you

* **TLS.** This serves plain HTTP. Put a reverse proxy in front for anything
  crossing an untrusted network (Caddy and nginx both terminate TLS for a
  single backend in ~10 lines of config; this project deliberately doesn't
  bundle one, to avoid maintaining certificate/domain logic orthogonal to the
  simulator itself).
* **Multi-user accounts.** One shared key, not per-user identity. See §1b.
* **Confining `POST /api/sections/scan` to a directory root.** The key makes
  `scan` **authenticated, not safe**. It takes an arbitrary `path` in the request
  body, and there is no allowlisted root — so any tester holding the key can
  enumerate and read any directory the *server process* can read, anywhere on the
  host, not just under `data_cache/`. Restricting `path` to a configured root was
  identified and deliberately **deferred as a MEDIUM**; it is not implemented, and
  the shared key does not substitute for it. Treat the key as equivalent to shell
  read access on the host, and hand it out on that basis. The controls that do
  exist are the ones in §0: run the server as a user that cannot read anything you
  would mind a tester listing.

  **The Docker image now does that for you** — it runs as `codon`, uid 10001,
  after a spell of running as root, i.e. of doing the precise opposite of the
  sentence above. §1c has the measured before/after and the one thing it costs
  (a `chown` on host bind mounts over writable paths). Running `codon studio`
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
[team's machines] --HTTPS--> [Caddy/nginx, TLS termination] --HTTP, localhost--> [codon studio --host 127.0.0.1]
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

Two more things this setup needs that the loopback default does *not* hand you,
both covered in §4: set `TMESIM_STUDIO_ALLOWED_HOSTS` to the hostname the proxy
forwards (the automatic loopback allowlist is `localhost,127.0.0.1`, which will
**400 every proxied request** arriving with your real hostname unless the proxy
rewrites `Host`), and leave `TMESIM_STUDIO_EXPOSE_TRACEBACKS` unset.

---

## 4. Environment variables

Five. There is no config file; these are the whole surface, and only the first
four are security-relevant.

**Four of the five are read once, at process start, and need a restart to
change.** `TMESIM_STUDIO_API_KEY`, `TMESIM_STUDIO_ALLOWED_HOSTS` and
`TMESIM_STUDIO_ALLOWED_ORIGINS` are read while the app is being built — each
decides whether a middleware is installed at all, so there is nothing left to
re-read afterwards — and `TMESIM_RECON_CACHE_SIZE` is read at import, when the
cache object is constructed.

**`TMESIM_STUDIO_EXPOSE_TRACEBACKS` is the exception: it is read per request**,
inside the `GET /api/runs/{id}/traceback` handler
(`app/server/api.py`). Practically that means the value in force is whatever is
in the *live process's* environment at the moment the route is called, not
whatever was set when it started — so a debugging session that mutates the
process environment takes effect immediately, and, less comfortably, so would
anything else that mutates it. It also means there is no restart-shaped moment
at which an operator gets to notice that traces are being served.

`codon studio` / `python -m app.server` writes two of them into its own
environment (`TMESIM_STUDIO_API_KEY`, `TMESIM_STUDIO_ALLOWED_HOSTS`), and in both
cases a value you exported yourself **beats the automatic behaviour and is never
clobbered** — so `-e` on `docker run`, `Environment=` in a systemd unit and a
plain shell export all compose correctly. An explicit `--api-key` on the command
line still overrides the exported key; that is the one deliberate exception, and
it is the ordering you'd expect.

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
  warning still prints, which is the point of it.

### `TMESIM_STUDIO_ALLOWED_HOSTS`

Comma-separated hostnames for Starlette's `TrustedHostMiddleware`. This is
defence in depth against DNS rebinding — a hostile domain that resolves to
`127.0.0.1` becomes same-origin with your server, which CORS cannot stop. It is
**not** the fix for a Host-header auth bypass; that is fixed at the gate itself.

* **Unset:** the middleware is not installed and there is no host validation at
  all — *except* that on an **IPv4 loopback** bind `codon studio` fills it in with
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
  validation at all.** `codon studio --host ::1` is *excluded* from the
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
  TMESIM_STUDIO_ALLOWED_HOSTS=localhost codon studio --host ::1   # reach it at http://localhost:8000/
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
  TMESIM_STUDIO_ALLOWED_ORIGINS=http://localhost:5173 codon studio
  ```

  Full origins including scheme and port, not bare hostnames.

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
  [codon] TMESIM_RECON_CACHE_SIZE='' is not usable: it is empty. Valid values are
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
