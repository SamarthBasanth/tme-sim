# Runbook — fresh clone to a completed run

`docs/DEPLOYMENT.md` explains the deployment *model*: what each mode means,
what is not handled for you, what every environment variable does. This page is
the other thing — **the shortest sequence of commands that ends in a finished
simulation, with nothing in it that has not been executed.**

Every command below was run, in this order, on **macOS 14 / arm64 with a
colima-backed Docker daemon (linux/arm64, 4 CPUs, 5.7 GB)** on **2026-08-22**,
and every output quoted is copied from that session. Where a number is
machine-specific it says so.

**What this runbook does not establish.** It was executed on **one
architecture**. `linux/amd64` resolves a different set of wheels from the same
version pins and is settled only by `.github/workflows/docker.yml`; until that
has gone green, treat amd64 as open. See `docs/DEPLOYMENT.md` §1c.

---

## Which path to take

| | use it when | ends at |
|---|---|---|
| **A. Docker** | you want the shipped artifact, isolation, and a non-root uid | §A6, a `done` run |
| **B. Local venv** | you are developing, or Docker is not available | §B3, a `done` run |

Both end in the same place. Path A is what a field tester should get.

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

---

## A. Docker

### A1. Build

```bash
docker build -t codon-studio .
```

Exits 0. Measured on the machine above: **58 s** for `docker build --no-cache`
(every layer re-run, wheels re-downloaded) and **51 s** for an ordinary build
reusing the apt layer. Both produce an image `docker images` reports as
**1.38 GB**.

> Do not read a size out of `docker image inspect --format '{{.Size}}'`. On this
> containerd-snapshotter daemon it reported **301,296,504 bytes** for the same
> image — a 4.6x under-report, because that field is compressed content size.
> CI reads that field. Cross-check against `docker images`.

### A2. Start it

`data_cache/` is mounted **read-only** — the container never needs to write to
your tissue — and run history goes to a **named volume**, which is the only
shape that works without a host-side `chown` (`docs/DEPLOYMENT.md` §1c).

```bash
docker run -d --name codon \
  -p 127.0.0.1:8401:8000 \
  -v "$(pwd)/data_cache:/app/data_cache:ro" \
  -v codon_runs:/app/app/runs \
  -e TMESIM_STUDIO_API_KEY=<a-real-secret> \
  codon-studio
```

### A3. Read the banner — it tells you the security posture it chose

```bash
docker logs codon
```

```
Codon Studio backend
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
docker exec codon sh -c 'id; ls /app/app/runs/'
```

```
uid=10001(codon) gid=10001(codon) groups=10001(codon)
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
docker run -d --name codon \
  -p 127.0.0.1:8401:8000 \
  -v "$(pwd)/data_cache:/app/data_cache:ro" \
  -v /srv/tissue:/tissue:ro \
  -v codon_runs:/app/app/runs \
  -e TMESIM_STUDIO_DATA_ROOTS=/app/data_cache:/tissue \
  -e TMESIM_STUDIO_API_KEY=<a-real-secret> \
  codon-studio
```

`:`-separated (`os.pathsep`), **not** commas — a directory name may contain a
comma. Keep `/app/data_cache` in the list even if your samples are elsewhere,
or the re-register command in `docs/DEPLOYMENT.md` §1d will `403`.

Before pointing this at anything, read `docs/DEPLOYMENT.md` §0: a
spatial-transcriptomics export carries identifiers in the **directory name**,
the section id and label, and files Codon never reads but which travel with the
directory.

### A8. Stop

```bash
docker stop codon      # exits 0 immediately; PID 1 handles SIGTERM, no 30 s wait
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
./.venv/bin/codon studio --port 8401
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
./.venv/bin/codon studio --port 8401 --api-key <a-real-secret>
TMESIM_STUDIO_API_KEY=<a-real-secret> ./.venv/bin/codon studio --port 8401
./.venv/bin/python -m app.server --port 8401 --api-key <a-real-secret>
```

Verified: banner reads `api key SECRET123 [set via --api-key]`, an unkeyed
`/api/health` is **401** and a keyed one is **200**.

> `codon studio --api-key` did not exist until 2026-08-22 — it exited 2 with
> `unrecognized arguments` while being the headline command in
> `docs/DEPLOYMENT.md` §1b. It was found by running the documentation.

### B2a. If you redirect stdout, set `PYTHONUNBUFFERED=1`

**This is the one thing on the local path that can cost you an API key.** The
startup banner goes to stdout, and CPython block-buffers stdout when it is a
pipe — systemd, `nohup`, `| tee`, any supervisor. Measured, same command, only
the environment differing:

| | log ordering |
|---|---|
| default | `INFO: Uvicorn running` line 4, `Codon Studio backend` line **5** |
| `PYTHONUNBUFFERED=1` | `Codon Studio backend` line 1, `INFO: Uvicorn running` line **30** |

By default the banner lands *after* the server is already serving, and until
the buffer fills it is not in the file at all. On a non-loopback bind that
banner carries an auto-generated key printed exactly once and persisted
nowhere. The Docker image sets `ENV PYTHONUNBUFFERED=1` for exactly this reason;
nothing sets it here.

```bash
PYTHONUNBUFFERED=1 ./.venv/bin/codon studio --host 0.0.0.0
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
| 4 | `codon studio --host 0.0.0.0 --api-key <secret>` — §1b's headline command — exited 2 | `tmesim/cli.py::cmd_studio` subparser | fixed in the CLI (`ef46459`), re-verified here |
| 5 | §1c and §2 asserted scan confinement was "opt-in and nothing here sets it"; the image now sets it on every default run | `docs/DEPLOYMENT.md` | fixed |
| 6 | §1c quoted a 4 s in-container run with no RunSpec attached; not reproducible | `docs/DEPLOYMENT.md` | fixed — 47.13 s with the spec stated |
| 7 | Startup banner (and any auto-generated key) is reordered behind uvicorn when stdout is a pipe; the `PYTHONUNBUFFERED` fix exists only in the `Dockerfile` | `app/server/__main__.py::main` | documented here and in §1b; **not fixed in code** |

Nothing in the container failed to build, start, or run.

Finding 7 is the one still open in code. It is documented rather than worked
around because the fix is a decision about the server's own output discipline —
flushing the banner, or forcing line buffering — and that is not a thing to
change quietly in a runbook.
