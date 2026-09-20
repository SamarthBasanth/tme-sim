# Runbook — fresh clone to a completed run

`docs/DEPLOYMENT.md` explains the deployment *model*: what each mode means,
what is not handled for you, what every environment variable does. This page is
the other thing — **the shortest sequence of commands that ends in a finished
simulation, with nothing in it that has not been executed.**

Every command below was run, in this order, on **macOS 14 / arm64 with a
colima-backed Docker daemon (linux/arm64, 4 CPUs, 5.7 GB)**, and every output
quoted is copied from a real session. Where a number is machine-specific it
says so.

**Two sessions, and the second one moved numbers the first one recorded.**
Paths A and B were executed on **2026-08-22**. §H (hosted) and every build
figure in §A1 were executed on **2026-09-20**, on the same machine, against the
same Dockerfile plus the one line §H4 exists for. The 2026-08-22 image size and
build times were measured *before* the Dockerfile grew its build-time
`vivum acquire` step and are no longer what the default build produces; §A1 now
carries the 2026-09-20 figures and says what changed. Nothing else was
re-measured, so treat a §A/§B number as of 2026-08-22 unless it says otherwise.

**One substitution in those transcripts, and only one.** That session ran under
the name **Codon**; the product was renamed back to **Vivum** later on
2026-08-22 (it was Vivum until 2026-08-14). The product name has been
substituted throughout the commands and quoted output below — the image tag, the
container name, the container's unix account, the startup banner. Nothing else
was touched: no timing, no line number, no verdict, no uid.

**What this runbook does not establish.** It was executed on **one
architecture**. `linux/amd64` resolves a different set of wheels from the same
version pins and is settled only by `.github/workflows/docker.yml`; until that
has gone green, treat amd64 as open. See `docs/DEPLOYMENT.md` §1c.

---

## Which path to take

| | use it when | ends at |
|---|---|---|
| **H. Hosted (Render)** | someone who is not you needs a URL | §H4, a live instance you have verified |
| **A. Docker** | you want the shipped artifact, isolation, and a non-root uid | §A6, a `done` run |
| **B. Local venv** | you are developing, or Docker is not available | §B3, a `done` run |

A and B end in the same place. Path A is what a field tester should get on
their own machine; **path H is the ten-minute one**, and it is the only one
that puts the product somewhere a person can click.

Read §H0 before §H1. It is one paragraph and it is the difference between a
demo and a place to keep work.

---

## 0. Before either path: is the port free?

Do this first. It costs one command and it is the failure that wastes the most
time, because on a VM-backed daemon Docker will **not** tell you:

```bash
lsof -nP -iTCP:8000 -sTCP:LISTEN
```

Empty output means free. This is not hypothetical — running it while writing
this runbook printed:

```
COMMAND   PID           USER   FD   TYPE             DEVICE SIZE/OFF NODE NAME
Python  32865 samarthbasanth    8u  IPv4 0x...      0t0  TCP 127.0.0.1:8000 (LISTEN)
```

so everything below uses **8401** instead. If you skip this check and publish
onto an occupied port, `docker run` still exits 0, the container still says
`running`, and every request you make reaches the *other* process — see
`docs/DEPLOYMENT.md` §1c for the measured version of that trap.

§0 is for paths A and B. **Path H publishes no local port and you can skip it.**

---

## H. Hosted, on Render — ten minutes

`render.yaml` in the repository root is this deployment as code. Read it: every
line of it is a thing that used to be a setting someone had to remember, and
the two that matter most fail *silently* when forgotten.

### H0. What you are deploying, and what it now keeps

**Corrected 2026-09-20.** This section previously opened "a hosted instance is
a place to DEMONSTRATE this product, not a place to keep work", and explained
at length why a disk could not fix it. That was true when it was written and is
no longer true. It is corrected here rather than quietly rewritten, because the
old paragraph is exactly the kind of text people act on — someone reading it
would not have added the disk that now works.

**What changed.** `TMESIM_STATE_DIR` makes `runs/`, `projects/` and `reports/`
children of ONE directory (`app/server/cache.py::state_dir`). The old obstacle
was real: Render allows one disk per service, the stores were siblings under
`/app/app`, and that parent also holds `server/` and `web/`, so mounting it
would hide the application. One configurable root removes the obstacle
entirely. `render.yaml` now declares the disk and sets the variable.

**The part that was not obvious, and cost a dead container to find.** A volume
mounted at a path the image does not already contain is created **root-owned**.
This image runs as `vivum` (uid 10001). The first thing the server does with
`TMESIM_STATE_DIR` is `mkdir` inside it, so the result was:

```
PermissionError: [Errno 13] Permission denied: '/var/lib/vivum/runs'
```

and an exited container — **the whole service dead at boot, caused by the one
setting whose purpose is to stop a redeploy destroying people's projects.**
Whoever turned persistence on would have taken the site down.

Two changes answer it. The image pre-creates `/var/lib/vivum` owned by uid
10001, so the mounted volume inherits that ownership. And
`app/server/cache.py::state_dir_problem` tests the directory with a real write
at boot and, if it cannot, exits **78** (`EX_CONFIG`) with one sentence naming
the running uid, the directory's owner, and the fix — instead of a traceback
out of `os.makedirs`. It deliberately does **not** fall back to the ephemeral
default: a fallback would look like working persistence and lose the work
anyway.

**What is measured, and what is still open.** Measured 2026-09-20 in the built
image under `--memory=512m`, on a volume mounted at `/var/lib/vivum`:

| Checked | Result |
|---|---|
| Mount ownership under `USER vivum` | inherited uid 10001; `runs/`, `projects/`, `reports/` created `0700` |
| Boot | `/healthz` 200 in ~2s, `OOMKilled=false` |
| Route table | `/` workspace, `/studio` Studio, `/index.html` Studio |
| Banner | `run history /var/lib/vivum/runs` |
| Survives a redeploy | project created via API, container `rm -f`'d and replaced on the same volume, project present with `created_at` byte-identical |
| Unwritable state dir | exit 78, one sentence, no traceback |

**Still open, and do not assume otherwise: that was Docker's volume semantics,
not Render's.** No Render account is attached to this repository, so Render's
own disk ownership is unverified from here. If Render mounts disks root-owned
regardless of the image, this service refuses to start and says so in one line
— see §H7 for the check and the two fixes.

**One more thing a hosted instance still does not do:** adding *your own*
tissue, because a section registers **by filesystem path on the server's own
machine** — see §H6.

Two things this deploy *does* do that a bare container would not:

* **It can actually run a simulation.** The image fetches one tissue section
  from its original source at build time (`ARG FETCH_SECTION=1`), so a hosted
  instance has tissue with no host to bind-mount from. Adding **your own**
  tissue to a hosted instance is a different problem, because a section
  registers **by filesystem path on the server's own machine** — see §H6.
* **It will not die on the first run.** See §H4, which is the check for that.

### H1. What is manual and what is automatic

Everything in the left column happens once, by a human, in a browser. Nothing
in it is triggered by pushing code.

| Manual, once | Automatic, from then on |
|---|---|
| Connect the repo and create the Blueprint (§H2) | Render reads `render.yaml` and creates the service |
| — | The API key is generated by Render and kept in the service's environment |
| — | `healthCheckPath: /healthz` is set from the file, not the dashboard |
| — | `docker build` runs on Render's builder; `.dockerignore` applies |
| — | One tissue section is fetched from source during the build |
| — | `$PORT` is assigned by Render and honoured by the image |
| — | The grid ceiling defaults to `auto` from the image's `ENV` |
| Read the key out of the dashboard (§H3) | — |
| Run the four checks in §H4 and read them | — |
| — | Every push to the linked branch rebuilds and redeploys |
| Change `autoDeployTrigger` to `checksPass` once anyone depends on this (§H5) | — |

The asymmetry is deliberate. Every automatic row used to be a manual one, and
two of them — the health-check path and the grid ceiling — have each taken a
service down in this project's history when a human forgot.

### H2. Create the Blueprint (manual)

Render dashboard → **New** → **Blueprint** → pick this repository → apply.
Render reads `render.yaml` from the repository root and creates one web
service, `vivum-studio`.

Nothing about this is triggered by a push. A Blueprint is created once by
hand; after that, pushes deploy.

**Choose the instance size here.** The image protects itself on a small one
(§H4) rather than assuming a large one, so a 512 MB instance is a legitimate
choice — it will refuse a grid it cannot afford instead of dying. It is the
resolution that costs you, not the service.

The first build fetches one tissue section over the network and installs the
whole scipy/rdkit stack. For scale, **measured 2026-09-20 on this repository's
own machine** (macOS/arm64, colima, 4 CPUs, builder cache empty, base image
already pulled):

```
docker build -t vivum-deploy:audit .          exit 0, 83 s
docker images                                 1.54 GB
docker build --build-arg FETCH_SECTION=0 .    exit 0, 19 s   (no tissue)
```

Render builds `linux/amd64` on its own hardware, so treat those as an order of
magnitude, not a prediction. What they *do* establish is that the build is
minutes, not hours, and that the tissue step is most of the difference.

### H3. Read the key (manual)

Dashboard → the service → **Environment** → `TMESIM_STUDIO_API_KEY`.

`render.yaml` sets `generateValue: true` on it, and that is not a convenience.
The image's `CMD` binds `0.0.0.0`, which `app/server/__main__.py` treats as
non-loopback: with no key set it **auto-generates one and prints it exactly
once to stdout, persisting it nowhere**. On a platform whose logs roll, the
service is then unusable by everyone including the person who deployed it.
Render minting the value instead means it is sitting in the dashboard,
readable and rotatable, and the startup banner says `[set via
TMESIM_STUDIO_API_KEY]` rather than `[generated (not persisted -- save it
now)]`.

**Before you give that key to anyone, read `docs/DEPLOYMENT.md` §2.** It
authenticates `POST /api/sections/scan`; it does not *confine* it. Treat it as
closer to read access as the server's uid than to a password.

### H4. Verify, in four curls

Do all four. Each one corresponds to a failure this project has actually had,
and every one of them is silent — the dashboard says "live" for all four.

```bash
B=https://<your-service>.onrender.com; K="X-API-Key: <the key from H3>"

curl -s -o /dev/null -w "%{http_code}\n" $B/healthz                 # 200
curl -s -o /dev/null -w "%{http_code}\n" -I $B/healthz              # 200  (HEAD)
curl -s -o /dev/null -w "%{http_code}\n" $B/                        # 200, no key
curl -s -H "$K" $B/api/account | python3 -m json.tool | grep grid   # see below
```

1. **`/healthz`, not `/api/health`.** Every `/api/` route is behind the key and
   a platform health checker cannot send a header. Pointed at `/api/health` a
   check got 401 forever, concluded the service never came up, and killed a
   container that was serving correctly. `render.yaml` sets the right path; this
   confirms it arrived.
2. **HEAD, not just GET.** `HEAD /healthz` returned **404** until the route was
   registered for both methods — while the key middleware had already been
   written to exempt HEAD on that path, so the gate was expecting a probe the
   router refused. A health check configured for HEAD would have marked a
   perfectly healthy service down. Nothing else in this repository probes with
   HEAD.
3. **The page loads without a key.** Static assets are deliberately ungated so
   the page can render and *prompt* for the key. A 401 here means first load is
   broken for everyone.
4. **The grid ceiling actually bound.** This is the one that decides whether the
   first run a visitor clicks takes the service down.

For (4), the field to read is
`deployment.limits.max_grid_cells_source`. **Measured 2026-09-20 in a real
container, reading `app.server.api._MAX_GRID_CELLS` directly:**

| the container | ceiling | `max_grid_cells_source` |
|---|---|---|
| `--memory=512m` | **46,070** | `auto (cgroup v2 memory.max)` |
| no memory limit | 40,000,000 | `auto (no small-host limit detected)` |
| `-e TMESIM_MAX_GRID_CELLS=100000` | 100,000 | `TMESIM_MAX_GRID_CELLS` |

The image sets `ENV TMESIM_MAX_GRID_CELLS=auto`, so a small host is protected
without anyone remembering.

**What that is worth was measured the same day, by turning it back off.**
`docker run --memory=512m -e TMESIM_MAX_GRID_CELLS=` — an empty value, which
restores the untouched default and is exactly how an image without this `ENV`
behaves. With the ceiling back at 40,000,000, a POST of the stock
`/api/defaults` spec (`dx_um = 12`, the value the product itself proposes) was
accepted with `202`, and the container was
**OOM-killed within a second**: `State.OOMKilled=true`, exit **137**, `/healthz`
unreachable. Not the run failing; the process. With the `ENV` in place the same
POST is a **400** naming the 46,070 ceiling and the container is still serving.
The first click a first user makes was the one that killed the instance.

If line (4) reads `auto (no small-host limit detected)` on an instance you know
to be small, the host is capping memory without declaring it through cgroup v2
`memory.max` or cgroup v1 `memory.limit_in_bytes`. Fix it by declaring the size
with `TMESIM_HOST_MEMORY_MB`, or by setting an integer ceiling outright — both
in `docs/DEPLOYMENT.md` §4. Note that setting the variable at all silences the
startup line that used to say "no memory limit is declared … if this host is
small, set `TMESIM_MAX_GRID_CELLS`", because that warning fires only while the
ceiling is the untouched default. That is why this check is a step in a runbook
rather than a line in a log.

### H5. Redeploys, and when to make them stricter

`autoDeployTrigger: commit` means every push to the linked branch rebuilds and
redeploys. That is the right setting while the product is being built, because
a deploy is how you find out.

**Change it to `checksPass` the moment anyone else depends on the URL.** The
two workflows in `.github/workflows/` are what that setting would then gate:

* `docker.yml` — the image on `linux/amd64`: that it builds at all, that it
  does not run as root, that `.dockerignore` kept the operator's tissue and
  unpublished compounds out of a layer.
* `product.yml` — the product a visitor meets: both HTML pages served **from
  inside a built container**, every `app/web/js/**` module parsing as an ES
  module, every file under `app/web/` reachable at its own URL, `$PORT`
  honoured, `/healthz` answering GET and HEAD unauthenticated, and the claim
  surface (`README.md`, `docs/**`, the four root contracts).

`product.yml` is a separate workflow on purpose: `docker.yml` skips pushes that
only touch `docs/**`, which is exactly where the claim surface lives.

### H6. Your own tissue on a hosted instance

A section registers **by filesystem path on the server's own machine**. Nothing
is transferred to the server; the bytes are read where they lie. A hosted
container has no filesystem of yours to read, so §A7's bind mount has no
equivalent here — that path exists on a machine you control, and a hosted
instance is not one.

The honest position: use a hosted instance for the section the image carries,
and use path A on your own machine for your own tissue. Read
`docs/DEPLOYMENT.md` §0 first either way — a spatial-transcriptomics export
carries identifiers in the **directory name**, the section id and label, and in
files this product never reads but which travel with the directory.

### H7. When it is wrong, it is usually one of these

| Symptom | Almost certainly |
|---|---|
| Every request hangs; the dashboard says live | `$PORT` is being ignored. The `CMD` must carry **no** `--port`; `app/server/__main__.py` defaults `--port` from the environment. A flag there overrides the platform and nothing listens where it routes. |
| The service is killed and restarted in a loop, logs full of 401 | The health check is pointed at `/api/health` instead of `/healthz`. `render.yaml` sets it; a dashboard edit can override it. |
| Health check fails, but every request you make works | The probe uses HEAD. See §H4 (2). |
| The instance dies mid-run, exit 137 | OOM. §H4 (4). |
| `/studio` is a 404 while `/` works | The static mount does not append `.html` to a path and does not know `/studio`; it needs its own route registered before the mount at `/`. `product.yml` fails on exactly this. |
| Projects you created yesterday are gone | The disk, or the variable that points at it. See the two checks below. |
| The service will not start and the log says `REFUSING TO START` | The disk is mounted somewhere this process cannot write. This is the good failure — the sentence names the running uid and the mount's owner. See below. |

**The disk, in two checks.** Both are one line each and answer different
questions.

*Is persistence actually on?* **Ask the running service, from anywhere.**
`GET /api/health` reports the two paths it is actually using, so this needs no
dashboard and no log access:

```bash
curl -s -H "X-API-Key: $KEY" https://<service>/api/health \
  | python3 -c 'import json,sys; d=json.load(sys.stdin); print(d["runs_dir"]); print(d["projects_dir"])'
```

Persistence is on when those read `/var/lib/vivum/runs` and
`/var/lib/vivum/projects`. If they read `/app/app/runs` and
`/app/app/projects`, `TMESIM_STATE_DIR` is not set on this service and
**nothing is being kept** — every deploy destroys it. The startup banner prints
the same fact as `run history …` if you do have log access.

Measured on the live service 2026-09-20, immediately after a push that added
the `disk:` block to `render.yaml`: `/app/app/runs`. The file had no effect,
which is the §H1 asymmetry doing exactly what it says — **a service created by
hand in the dashboard auto-deploys on push and ignores this file.** Do not
infer from a green deploy that a `render.yaml` change was applied; ask the
service. The most likely cause is
the one in the header of `render.yaml`: a service created by hand in the
dashboard auto-deploys on push but **ignores `render.yaml` entirely**, so the
`disk:` block and the variable below it never applied. Check
dashboard → the service → Settings; a Blueprint-managed service says so. Fix it
by adding the disk and the variable in the dashboard, or by recreating the
service from the Blueprint (§H2).

*Did the mount come up root-owned?* If Render mounts disks root-owned
regardless of the image, the service refuses to start and the log carries one
sentence naming `uid 10001`, the directory's actual owner, and this fix. Two
options, in order:

1. Mount the disk at **`/var/lib/vivum`**, which the image pre-creates and
   chowns to uid 10001 (`Dockerfile`). A volume mounted on a path the image
   already owns inherits that ownership — measured. This is what `render.yaml`
   does and why the path is not arbitrary.
2. If the platform overrides that anyway, the mount has to be chowned to 10001
   before the app starts. That needs a root step ahead of `USER vivum`, which
   this image deliberately does not have — treat it as a change to make
   consciously, not a quick patch.

Do **not** "fix" it by unsetting `TMESIM_STATE_DIR` to get the service back up.
That starts cleanly and silently returns you to losing every project on every
deploy, which is the failure this whole section exists to prevent.

---

## A. Docker

### A1. Build

```bash
docker build -t vivum-studio .
```

Exits 0. **Re-measured 2026-09-20** on the machine above, because the figures
this section used to carry were taken before the Dockerfile grew its build-time
`vivum acquire` step and were no longer what the default build produces:

| | exit | wall clock | `docker images` |
|---|---|---|---|
| builder cache empty, base image already pulled | 0 | **83 s** | **1.54 GB** |
| rebuild after a one-line Dockerfile edit (pip layer cached, tissue re-fetched) | 0 | 48 s | 1.54 GB |
| `--build-arg FETCH_SECTION=0` (no tissue) | 0 | 19 s | — |

Those are wall-clock on a shared laptop, not benchmarks: a repeat of the third
row later the same day, while several other builds and test runs were using the
same machine, took **187 s** for the identical command and the identical result.
Read the column as an order of magnitude and the `exit` column as the fact.

This section previously read **58 s / 51 s / 1.38 GB**. Those were honest
measurements of an image that did not yet fetch a tissue section at build time;
the ~160 MB and the extra ~30 s are that step, and `docker images` is where they
show up. The default build is `FETCH_SECTION=1`, so 1.54 GB is what you get
unless you ask for the other one.

> **Check the exit code. Do not grep the output for a success string.** A
> previous session in this repository read a build log for "Compiled
> successfully" while the build was exiting 1 the whole time. In `zsh`,
> `${PIPESTATUS[0]}` is empty — redirect to a file and test `$?`.

> Do not read a size out of `docker image inspect --format '{{.Size}}'`. On this
> containerd-snapshotter daemon it reported **365,591,455 bytes** for the
> 1.54 GB image measured above — a 4.2x under-report, because that field is
> compressed content size (it was 301,296,504 against 1.38 GB in the earlier
> session; the ratio is a property of the graph driver, not of the image).
> `.github/workflows/docker.yml` reads that field and prints both. Cross-check
> against `docker images`.

> `docker build --progress=plain` is a **BuildKit** flag. This daemon ships the
> legacy builder with no `buildx` component, where it exits **125** with
> `unknown flag: --progress` — in under a second, which reads exactly like a
> fast build if you are only watching the clock.

### A2. Start it

`data_cache/` is mounted **read-only** — the container never needs to write to
your tissue — and run history goes to a **named volume**, which is the only
shape that works without a host-side `chown` (`docs/DEPLOYMENT.md` §1c).

```bash
docker run -d --name vivum \
  -p 127.0.0.1:8401:8000 \
  -v "$(pwd)/data_cache:/app/data_cache:ro" \
  -v vivum_runs:/app/app/runs \
  -e TMESIM_STUDIO_API_KEY=<a-real-secret> \
  vivum-studio
```

### A3. Read the banner — it tells you the security posture it chose

```bash
docker logs vivum
```

```
Vivum Studio backend
  repo        /app
  run history /app/app/runs
  recon cache maxsize=4
  serving     http://0.0.0.0:8000
  your tissue http://0.0.0.0:8000/ingest
  api key     <a-real-secret>  [set via TMESIM_STUDIO_API_KEY]
              every request needs header: X-API-Key: <a-real-secret>
  tissue roots /app/data_cache  [defaulted for this non-loopback bind]
              tissue outside that directory will be REFUSED; set TMESIM_STUDIO_DATA_ROOTS to include yours
INFO:     Uvicorn running on http://0.0.0.0:8000 (Press CTRL+C to quit)
```

Two lines are load-bearing.

* **`api key`** — if you had not passed one, this would read `[generated (not
  persisted -- save it now)]` and the key would appear exactly once, here.
  Copy it before you lose the process; nothing writes it to disk.
* **`tissue roots`** — the default `CMD` binds `0.0.0.0`, so the scanner is
  confined to the mount point automatically. **If your samples are somewhere
  else, this line is telling you they will be refused.** §A7 covers that.

### A4. Confirm the gate is real

```bash
curl -s -o /dev/null -w "%{http_code}\n" http://127.0.0.1:8401/api/health                       # 401
curl -s -o /dev/null -w "%{http_code}\n" -H "X-API-Key: wrong" http://127.0.0.1:8401/api/health # 401
curl -s -H "X-API-Key: <a-real-secret>" http://127.0.0.1:8401/api/health                        # 200
curl -s -o /dev/null -w "%{http_code}\n" http://127.0.0.1:8401/                                 # 200
```

The last one is not a bug: static front-end files are deliberately ungated so
the page can load and *prompt* for the key rather than failing to render.

### A5. Confirm the tissue actually mounted

```bash
curl -s -H "X-API-Key: <a-real-secret>" http://127.0.0.1:8401/api/sections
```

With `data_cache/` mounted read-only, all **7** built-in sections report
`available: true`. With nothing mounted the server still boots and reports each
one unavailable *with a reason* — that is the honest-failure path, not a
broken container, so check this before concluding the image is at fault.

### A6. Run a real simulation

Not a health check. This is a reconstruction plus a transport solve plus
metrics, on real tissue.

```bash
B=http://127.0.0.1:8401; K="X-API-Key: <a-real-secret>"

SPEC=$(curl -s -H "$K" $B/api/defaults)
RID=$(curl -s -X POST -H "$K" -H "Content-Type: application/json" \
        -d "$SPEC" $B/api/runs \
      | python3 -c "import json,sys; print(json.load(sys.stdin)['run_id'])")

while :; do
  ST=$(curl -s -H "$K" $B/api/runs/$RID/status \
       | python3 -c "import json,sys; d=json.load(sys.stdin); print(d['status'], d['stage'], d['progress'])")
  echo "$ST"; case "$ST" in done*|error*) break;; esac; sleep 3
done

curl -s -H "$K" $B/api/runs/$RID | python3 -m json.tool | head -40
```

Observed progression: `queued` -> `reconstruct` -> `simulate` -> `metrics` ->
`done`. On the stock `/api/defaults` spec (`xenium_319`, temozolomide,
`dx 12 um`, 6 h, `plasma_peak 50 uM`):

```
n_cells        40887      n_tumor 23812
grid           [424, 675] @ dx 12.0 um
runtime_s      47.13
kill_fraction  0.11221232991768859
```

Two submissions in the same container returned a **bit-identical**
`kill_fraction`. The binary payload routes answer too — `/geometry` 368,003 B,
`/field/oxygen` and `/field/drug_auc` 1,144,828 B each, `/mechanism` 45,772 B —
which is what tells you the front end will actually render this run and not
just list it.

**Read the run's caveats before you quote any of it.** The submit response
returns `warnings` immediately, and this section's own panel measures no
hypoxia, so absolute kill must not be quoted off it under any state route. That
is a property of the tissue, not a defect in the deployment.

Confirm the history landed in the volume rather than in the container's
writable layer:

```bash
docker exec vivum sh -c 'id; ls /app/app/runs/'
```

```
uid=10001(vivum) gid=10001(vivum) groups=10001(vivum)
_cfl_speeds.json  _section_meta.json  r_00001
```

The server runs as **uid 10001, not root**, and the named volume inherited that
ownership with no host-side step.

### A7. Only if your own tissue lives outside `data_cache/`

Sections register **by filesystem path on the server's own machine**. Nothing
is uploaded; the bytes are read where they lie. So the directory has to be
visible *inside the container*, and it has to be inside the scan allowlist —
which, per §A3, defaulted to `/app/data_cache` alone.

```bash
docker run -d --name vivum \
  -p 127.0.0.1:8401:8000 \
  -v "$(pwd)/data_cache:/app/data_cache:ro" \
  -v /srv/tissue:/tissue:ro \
  -v vivum_runs:/app/app/runs \
  -e TMESIM_STUDIO_DATA_ROOTS=/app/data_cache:/tissue \
  -e TMESIM_STUDIO_API_KEY=<a-real-secret> \
  vivum-studio
```

`:`-separated (`os.pathsep`), **not** commas — a directory name may contain a
comma. Keep `/app/data_cache` in the list even if your samples are elsewhere,
or the re-register command in `docs/DEPLOYMENT.md` §1d will `403`.

Before pointing this at anything, read `docs/DEPLOYMENT.md` §0: a
spatial-transcriptomics export carries identifiers in the **directory name**,
the section id and label, and files Vivum never reads but which travel with the
directory.

### A8. Stop

```bash
docker stop vivum      # exits 0 immediately; PID 1 handles SIGTERM, no 30 s wait
```

---

## B. Local virtualenv (no Docker)

```bash
python3 -m venv .venv
./.venv/bin/python -m pip install --upgrade pip          # PEP 660 needs pip >= 21.3
./.venv/bin/python -m pip install -c constraints.txt -e '.[app]'
```

`requires-python` is `>=3.9`. This was executed from scratch in a throwaway
venv on **Python 3.9.6** (the container is 3.11, so both ends of the supported
range are now exercised — 3.9 here, 3.11 in §A1) and installed clean. The
numerically load-bearing versions it resolved:

```
numpy 2.0.2 | scipy 1.13.1 | pandas 2.3.3 | rdkit 2025.09.2
```

Those are the pinned ones, and `import tmesim, fastapi, uvicorn, rdkit, h5py,
skimage` and `import app.server.api` all succeed afterwards.

`-c constraints.txt` is not optional if you care about reproducing published
numbers. `pyproject.toml` declares bounds, not versions; without the
constraints file two installs a month apart resolve different
numpy/scipy/pandas and are not numerically identical. See `docs/REPRODUCE.md`
§1a. Note that pip resolves *unpinned* packages freely even with the file —
fastapi/starlette/uvicorn came out at different versions here than in the
image, which is expected and does not touch a number.

### B2. Start it

```bash
./.venv/bin/vivum studio --port 8401
```

Binds `127.0.0.1`. **No key is generated and none is required** — that is the
intended single-user workflow, and the banner says so in as many words on every
keyless start. Read `docs/DEPLOYMENT.md` §1e before you put *anything* in front
of this: ngrok, `ssh -R`, VS Code port forwarding and any non-authenticating
reverse proxy all make it public while the process still reports a loopback
bind and still generates no key.

Verified on this path: `/api/health` 200 without any key, `/` 200, and
`/api/sections` reporting **38 sections, 38 available** — more than the
container's 7, because the host has ingested sections registered in
`app/runs/_ingested_sections.json` that the image does not carry.

To gate it anyway — all three are equivalent:

```bash
./.venv/bin/vivum studio --port 8401 --api-key <a-real-secret>
TMESIM_STUDIO_API_KEY=<a-real-secret> ./.venv/bin/vivum studio --port 8401
./.venv/bin/python -m app.server --port 8401 --api-key <a-real-secret>
```

Verified: banner reads `api key SECRET123 [set via --api-key]`, an unkeyed
`/api/health` is **401** and a keyed one is **200**.

> `vivum studio --api-key` did not exist until 2026-08-22 — it exited 2 with
> `unrecognized arguments` while being the headline command in
> `docs/DEPLOYMENT.md` §1b. It was found by running the documentation.

### B2a. If you redirect stdout, set `PYTHONUNBUFFERED=1`

**This is the one thing on the local path that can cost you an API key.** The
startup banner goes to stdout, and CPython block-buffers stdout when it is a
pipe — systemd, `nohup`, `| tee`, any supervisor. Measured, same command, only
the environment differing:

| | log ordering |
|---|---|
| default | `INFO: Uvicorn running` line 4, `Vivum Studio backend` line **5** |
| `PYTHONUNBUFFERED=1` | `Vivum Studio backend` line 1, `INFO: Uvicorn running` line **30** |

By default the banner lands *after* the server is already serving, and until
the buffer fills it is not in the file at all. On a non-loopback bind that
banner carries an auto-generated key printed exactly once and persisted
nowhere. The Docker image sets `ENV PYTHONUNBUFFERED=1` for exactly this reason;
nothing sets it here.

```bash
PYTHONUNBUFFERED=1 ./.venv/bin/vivum studio --host 0.0.0.0
```

### B3. The rest is identical

`/api/health`, `/api/sections` and the run submission in **§A4-A6** are
identical here — same routes, same payloads, same spec. Point `B` at
`http://127.0.0.1:8401` and drop the `X-API-Key` header if you started without
a key.

---

## The static export

`results/studio_export.html` is the whole Studio in one file: no server, no
network request, real precomputed runs. Rebuild it with

```bash
./.venv/bin/python scripts/build_studio_export.py
```

and verify it by driving it, not by opening it:

```bash
node scripts/studio_export/drive.mjs
```

That launches headless Chrome, loads the file from `file://`, switches runs,
decodes every field overlay and pushes it to the GPU, picks a cell, and asserts
that a withdrawn claim is bannered. It also asserts **0 non-`file:`/`data:` network
requests**, which is the property that makes the file safe to hand to someone.

The build refuses rather than shipping a quiet defect. It fails if the baked
front end grows an API route the static adapter does not answer, if a run is
baked without its claim-set verdict, or if a retired sentence appears anywhere
in the page's own copy. All three were confirmed by planting them.

**Two things to know before you rebuild and commit it.**

* **It bakes `app/web/` as it is on disk, including uncommitted work.** The
  export is a build artifact of the live front end, so rebuilding it in the
  middle of someone's UI change bakes that change into a git-tracked 10 MB
  file. Rebuild when the front end is settled, not opportunistically.
* **The run set moved, and it was previously wrong.** The build reads
  `results/_export_run_ids.json`, which points at `results/export_runs/`
  (`r_00001`-`r_00003`: temozolomide, cisplatin, doxorubicin — kill 1.84 %,
  43.62 %, 25.95 %). The exports that had been committed were built from
  something else: `r_00027`-`r_00029`, **all three temozolomide with a
  byte-identical `kill_fraction` of 0.025323366369897515**, i.e. one run baked
  three times, which left the A/B comparison panel comparing a run against
  itself. Building from the manifest fixes that and is the reproducible path.

  The one thing it costs: those manifest runs carry no `timelapse.bin`, so the
  time-lapse scrubber has no series in the current export. That degrades
  cleanly — `main.js::loadTimelapse` clears the scrubber and the reason rides
  in `RunResult.timelapse`, and the driver confirms no uncaught exception — but
  it is a real feature the older file demonstrated. Restoring it means
  re-running the three export runs with time-lapse capture enabled, which is a
  product decision about the flagship artifact, not a build fix.

---

## Everything this runbook found broken

Recorded because a runbook whose commands all worked first time is a runbook
that was not actually run.

| # | What | Where | Status |
|---|---|---|---|
| 1 | The export baked the **disk** view of a run, not the client view, so three runs shipped a claim withdrawn on 2026-08-21 with no banner | `scripts/build_studio_export.py::bake_runs` | fixed — goes through `Job.result_view()` |
| 2 | The export **could not be rebuilt at all**: the static API adapter was missing 8 exports the front end had grown | `scripts/studio_export/api_baked.js` | fixed — honest refusals, plus a documented no-op for `setApiKey` |
| 3 | The static reports panel still carried the withdrawn ranking claim as **current page copy** | `scripts/studio_export/reports_static.js` | fixed — matches the live panel |
| 4 | `vivum studio --host 0.0.0.0 --api-key <secret>` — §1b's headline command — exited 2 | `tmesim/cli.py::cmd_studio` subparser | fixed in the CLI (`ef46459`), re-verified here |
| 5 | §1c and §2 asserted scan confinement was "opt-in and nothing here sets it"; the image now sets it on every default run | `docs/DEPLOYMENT.md` | fixed |
| 6 | §1c quoted a 4 s in-container run with no RunSpec attached; not reproducible | `docs/DEPLOYMENT.md` | fixed — 47.13 s with the spec stated |
| 7 | Startup banner (and any auto-generated key) is reordered behind uvicorn when stdout is a pipe; the `PYTHONUNBUFFERED` fix exists only in the `Dockerfile` | `app/server/__main__.py::main` | documented here and in §1b; **not fixed in code** |

Nothing in the container failed to build, start, or run.

### And what the 2026-09-20 session found

| # | What | Where | Status |
|---|---|---|---|
| 8 | On a 512 MB container the stock `/api/defaults` spec was **accepted** and OOM-killed the process in about a second — exit 137, `/healthz` unreachable. The first run a first user clicks. | image had no ceiling default | fixed — `ENV TMESIM_MAX_GRID_CELLS=auto` in the `Dockerfile`; the same POST is now a 400 naming the ceiling |
| 9 | §A1's **1.38 GB / 58 s / 51 s** were measured before the build-time tissue fetch existed and no longer described the default build | `docs/RUNBOOK_DEPLOY.md` | fixed — re-measured, and what changed is stated |
| 10 | `GET /` served the **Studio** even with `app/web/home.html` present in the image, and `GET /studio` was a **404** | `app/server/api.py::_mount_web` | fixed the same day, elsewhere — `StaticFiles(html=True)` resolves `/` to the directory index and never appends `.html`, so both URLs needed explicit routes ahead of the mount. Re-measured after they landed: `/` serves `home.html` and `/studio` serves `index.html`, byte for byte, and `/studio?project=…` is 200. `product.yml` asserts exactly that, and was red on this before the fix. |
| 11 | A hosted instance keeps nothing: no disk, and one disk cannot cover `app/runs/` and `app/projects/` because neither path is configurable | `cache.py::runs_dir`, `projects.py::projects_dir` | fixed 2026-09-20 — `TMESIM_STATE_DIR` makes `runs/`, `projects/` and `reports/` children of one root, so one disk covers all three. Enabling it then **killed the container at boot**: a volume mounted where the image has no directory is root-owned, the image runs as uid 10001, `PermissionError` out of `os.makedirs`. Answered by pre-creating and chowning `/var/lib/vivum` in the image, and by `state_dir_problem` exiting 78 with a sentence instead of a traceback. Proven by creating a project, destroying the container and replacing it on the same volume. Render's own disk ownership remains unverified — §H7. |
| 12 | `node --check` — the obvious spelling of an ES-module syntax check — parsed a duplicate `export const` and an `await` outside an async function **without error** | CI that did not exist yet | fixed by construction — `product.yml` uses `node --input-type=module --check`, confirmed by planting the duplicate export in a real module and watching the wrong flag stay green |

There is one thing this session could not settle, and it is the same one §A
could not: **nothing here has been deployed to Render.** `render.yaml` was
checked against the Blueprint specification and its keys are current, and the
image it deploys was built and exercised — but no account is attached to this
repository, so §H2 and §H3 are the two steps in this runbook that have been
written and not run. Everything in §H4 was executed against a local container
on the port a platform would assign.

Finding 7 is the one still open in code. It is documented rather than worked
around because the fix is a decision about the server's own output discipline —
flushing the banner, or forcing line buffering — and that is not a thing to
change quietly in a runbook.
