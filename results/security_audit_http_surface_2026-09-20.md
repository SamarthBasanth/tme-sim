# Security audit — the HTTP surface that shipped 2026-09-20

Authorised review of the operator's own deployment, run 2026-09-20 against the live
service at `https://vivum.onrender.com` (read-only probes) and against a local instance
for everything intrusive. The live demo was never subjected to a denial-of-service; the
one resource finding below was demonstrated locally and is reproduced here as a
measurement, not as an attack on the live host.

**Scope — the surface with no prior security pass.** `results/security_audit.md`
(2026-08-14) predates all of it.

- `app/server/workspace.py` — `GET /api/workspace`, the home page's single boot call
- `app/server/projects.py` — the project store
- `app/server/cross_section.py`, `app/server/uncertainty_api.py`
- `app/server/account.py` — `GET /api/limits`, `GET /api/account`
- the static routes in `app/server/api.py::_mount_web` (`/`, `/studio`, `/index.html`)

**Counts: 1 HIGH (fixed) · 2 MEDIUM (reported) · 1 LOW (reported) · 6 CLEAN.**

The headline is that the four defect classes this codebase has actually suffered —
path traversal, exception text in the browser, an unconfined scanner, an ungated
route — are all **closed and verified closed**. The new finding is a different shape:
a store that validated its text for *emptiness* and never for *size*.

---

## HIGH — H1. One request permanently bloats the home page's boot call

**Status: FIXED in `app/server/projects.py`. Test: `tests/test_project_text_limits.py`.**

`ProjectStore.create` / `update` / `duplicate` / `add_specimen` / `update_specimen`
checked their free-text fields for emptiness and nothing else. A project name was
refused at zero characters and accepted at twenty million.

That is not cosmetic, because of where the text goes. `ProjectStore._summary` puts
`notes` into **every row** of `store.list()`, and `store.list()` is what
`GET /api/workspace` returns — the single boot call for the new home page.

### The triggering request

```
POST /api/projects
X-API-Key: <key>
Content-Type: application/json
Transfer-Encoding: chunked

{"name": "pwn", "notes": "<20 MB of text>"}
```

### What was measured (local instance, 2026-09-20)

| | before | after one request |
|---|---|---|
| `GET /api/workspace` | 163 B, 0.006 s | **20,972,002 B**, 0.092 s |
| `project.json` on disk | — | 20,971,682 B |

### What an attacker gets

A **persistent** denial of the product's front door. The payload is on disk, so it
survives a restart and a redeploy; it is served to *every* visitor on *every* boot;
and it is removed only by finding the project directory and deleting it by hand.
Repeating the request multiplies it. The cost to the attacker is one HTTP call.

This needs the API key — but the key is a single shared secret handed to field
testers, and `docs/DEPLOYMENT.md` §2 characterises it as *read* access as the
server's uid. This is a write that degrades the service for everyone else, which is
strictly more than that document leads an operator to expect.

### Why the 4 MB body cap did not stop it

`api.py::_cap_request_body` refuses on the **`Content-Length` header**. A request sent
`Transfer-Encoding: chunked` does not have one, so the cap never fires. Measured on
the same endpoint, same body:

```
Content-Length: 8388620            -> 413   (cap fires)
Transfer-Encoding: chunked, 8 MB   -> 403   (cap never fires; body fully buffered)
```

That is a separate defect, in a file this pass does not own — see M1.

### The fix

Two named constants and one gate, `ProjectStore._text`, applied at all five entry
points:

- `MAX_NAME_CHARS = 200` — project name and specimen label
- `MAX_NOTES_CHARS = 20_000` — project notes and specimen notes

Three properties worth stating:

1. **The store decides what it will persist.** The caps are deliberately *not* a
   restatement of the body cap. A request-size limit answers "how much may you send
   me"; this answers "how much am I willing to store and then hand to every future
   reader". The second question is the store's, and its answer must hold even when
   the first one is wrong — which, per M1, it currently is.
2. **The refusal keeps the `"<field>: <what is wrong>"` shape.** That prefix is
   load-bearing: `app/web/js/home/createProject.js::attachFieldError` parses it to
   decide which input to mark red. A message that loses it is silently demoted to a
   form-level error. The two pre-existing required-field wordings are preserved
   byte for byte.
3. **`duplicate` truncates instead of refusing.** It appends `" (copy)"`, so a source
   name *at* the cap would otherwise exceed it — and the caller supplied no name, so a
   400 would be a refusal they cannot act on.

### Plant evidence

- Test written first. Against the unfixed store (constants added, enforcement not):
  **9 failed, 2 passed** — failing on *behaviour*, not on import. The end-to-end
  assertion read `30060550 B for 3 projects, budget 68792 B`.
- Fix applied: **11 passed**.
- Defect planted back (the `notes` cap in `create` alone, replaced by
  `str(notes or "")`): **2 failed** — `test_create_refuses_oversize_notes` and
  `test_one_project_cannot_make_the_workspace_boot_call_unbounded`. Exactly the unit
  refusal and the end-to-end property, which is the pair that should catch it.
- Plant reverted; `grep PLANTED` clean; **11 passed**.
- Regression: `test_projects.py`, `test_workspace.py`, `test_home_page.py`,
  `test_api_hostile_input.py`, `test_project_text_limits.py` — **492 passed**.
- Verified against a running server: the same chunked 20 MB request now returns 400,
  nothing is written to the projects directory, and `GET /api/workspace` stays at
  419 B.

---

## MEDIUM — M1. The request body cap is bypassed by `Transfer-Encoding: chunked`

**Status: REPORTED, not fixed — `app/server/api.py` is out of scope for this pass.**

`api.py::_cap_request_body` (≈line 1830) reads `Content-Length` and refuses above
`MAX_REQUEST_BODY_BYTES` (4 MB). A chunked request has no `Content-Length`, so the
body reaches `await request.json()` in full. The docstring is honest that the header
"may be absent or a lie" and points at `_capped_body` — but **`_capped_body` has
exactly one user, `app/server/protocol_api.py`**. Every other JSON body on the server
is uncapped in practice: `api.py:468`, `:2010` (`/api/sections/scan`), `:2178`,
`:2188`, `:2314`.

There is a second-order amplifier on `/api/sections/scan`. `cache.PathNotAllowed`
carries the caller's path verbatim and `as_detail()` emits it **twice** (once inside
`message`, once as `path`), so the refusal is larger than the request:

```
request  8,388,620 B  (chunked)   ->   response 16,777,881 B, status 403
```

An 8 MB request costs the server 8 MB of buffered body plus a 16 MB response it
builds in memory. Suggested fix: count the stream for every route, not just the
protocol router — i.e. make `_capped_body` the default path — and truncate the echoed
path in `PathNotAllowed.as_detail()`.

## MEDIUM — M2. The redaction sweep is RED on `main`, so it is not guarding anything

**Status: REPORTED, not fixed — `app/server/cache.py` is out of scope and is being
edited by another agent.**

`tests/test_response_redaction.py` currently **fails**:

```
FAILED test_no_undeclared_exception_text_on_the_http_surface
  app/server/cache.py:295  state_dir_problem  [f-string]  f"...{type(exc).__name__}: {exc}..."
```

This is **not** an artefact of another agent's in-flight edit to `cache.py`:
`git show HEAD:app/server/cache.py` carries the same f-string at the same line 295,
and the line does not appear in the working-tree diff. The check is red on committed
`main`. The rest of that run was green — `1 failed, 387 passed` across
`test_api_hostile_input.py`, `test_scan_confinement.py` and `test_response_redaction.py`.

`state_dir_problem` is called only from `__main__.py:202` and prints to stderr before
`return 78`, so **this is not itself a live HTTP leak**. The finding is the state, not
the string: the guard that exists to stop exception text reaching the browser is red,
and a red guard cannot tell anyone that a *new* leak has appeared. It needs either an
`EXEMPTIONS` entry with a reason, or the call moved out of the swept tree.

This is `results/codon_vacuous_green_checks.md`'s pattern in its other form — not a
check that passes without testing its subject, but a check whose failure has stopped
carrying information.

### The sweep can also be evaded, and the evasions are the dangerous ones

`leak_kind` matches f-strings, `str`/`repr`/`format`, `.format`/`.join`, and `%`/`+`.
It is blind to bare attribute access. Verified by calling the sweep's own
`scan_source` on synthetic modules:

| shape | swept? |
|---|---|
| `detail={"error": exc.filename}` on an `OSError` | **EVADES** |
| `detail={"error": exc.args[0]}` | **EVADES** |
| `{"error": exc.strerror, "where": exc.filename}` | **EVADES** |
| `{"error": traceback.format_exc()}` | **EVADES** |
| `{"error": f"{exc}"}` (control) | CAUGHT |

The first is the one that matters: **`OSError.filename` *is* the absolute server path**
that `redact_exception` exists to hide, and it reaches a response without tripping the
guard. Four current sites rely on bare `.message`/`.detail` and pass undeclared rather
than declared — `cross_section.py:1085`, `api.py:2702`, `api.py:2131`,
`projects.py:722`. All four are safe today because those attributes carry authored
refusal text; none is protected by the sweep if that changes.

## LOW — L1. `SPECIMEN_ID_RE` is defined and never applied

`app/server/projects.py:126` defines `SPECIMEN_ID_RE` and nothing uses it;
`_find_specimen` relies on an in-document equality scan. Harmless today because
`specimen_id` never becomes a path segment. Recorded because it is precisely the
shape of the traversal this file already suffered — the file's own comment at
`_dir` says the guard "existed and was not wired to the boundary that needed it".

---

## CLEAN — what was checked and found sound

A clean result is a finding. Each of these was tested, not read.

**C1. Path traversal is closed, across every encoding tried.** `GET` and `DELETE` on
`/api/projects/{id}` with `..`, `%2e%2e`, `%252e%252e`, `%2e%2e%2f%2e%2e`, `.%2e`,
`%c0%ae%c0%ae`, fullwidth `．．`, `%00`, `..%2f..%2fetc`, `%2Fetc%2Fpasswd`, a
fullwidth digit, RTL-override `%e2%80%ae`, and `p_00001/../../etc` — **every one 404
or 405. No 200, no 500.** The NUL byte that used to raise `ValueError` and 500 now
returns 404.

**C2. No newer module builds a path from an unvalidated request value.** An
exhaustive map of `os.path.join` / `Path(` / `open(` / `os.remove` / `shutil.rmtree` /
`os.listdir` / `glob` across all five files:

| file | path-building sites | request-derived | validated first |
|---|---|---|---|
| `cross_section.py` | **0** — does not even import `os` | — | — |
| `workspace.py` | **0** — delegates to `ProjectStore` | — | — |
| `account.py` | 1 (`open`) | no — two cgroup constants | — |
| `uncertainty_api.py` | 4 | `run_id` | yes, `_RUN_ID_RE` then `jobs.get()` |
| `projects.py` | 13 | `project_id`, `run_id` | yes, `PROJECT_ID_RE` / `RUN_ID_RE` |

`PROJECT_ID_RE` is enforced *inside* `_dir`, so no call site can skip it. The two
modules the brief flagged hardest, `cross_section.py` and `uncertainty_api.py`, are
the two with the least filesystem contact.

**C3. Scanner confinement holds, including against symlinks and relative paths.**
`confine_scan_path` resolves with `realpath` (never `realpath(abspath(...))`) and
compares with `commonpath`. Tested directly: `..` escape, a symlink pointing outside
the root, symlink-plus-`..`, an absolute outside path, `~`, `/`, `/etc/passwd` and a
bare relative path — **all refused**. A literal `%2e%2e` directory *inside* the root
is allowed, which is correct: it resolves inside the root and is not an escape.

**C4. The confinement is actually in force on the live deployment.** The Docker
`ENTRYPOINT` is `python -m app.server`, so `__main__.py:189` runs and auto-sets
`TMESIM_STUDIO_DATA_ROOTS` on the non-loopback bind. Confirmed on the live service:
`GET /api/account` → `deployment.data_roots_confined: true`.

**C5. No new route is accidentally ungated, and HEAD is handled.** All of
`/api/workspace`, `/api/account`, `/api/limits`, `/api/cross-section`,
`/api/uncertainty/spans`, `/api/projects`, `/api/sections/scan` return 401 without the
key, on the live service and locally. `HEAD /api/workspace` → 401. `POST /healthz` →
401. Bypass shapes `//api/…`, `/./api/…`, `/api/../api/…`, `/%2e/api/…`,
`/api%2fworkspace`, `/API/…` yield 401, 404 or 405 — never data. The Starlette
Host-header auth bypass (GHSA-86qp-5c8j-p5mr, C1 in the 2026-08-14 audit) is fixed:
`Host: 127.0.0.1:8000/../`, `x/../` and `testserver/../..` all return **401**.

**C6. The ungated routes leak nothing.** `/healthz`, `/api/health`, `/openapi.json`,
`/docs` and `/ingest` answer 200 anonymously **by design**, and
`tests/test_api_hostile_input.py::test_the_ungated_documentation_routes_disclose_no_operational_detail`
already pins that they carry no key, no repo root and no `runs_dir` / `projects_dir` /
`reconstruction_cache`. Verified live: anonymous `GET /api/health` returns exactly
`{"ok":true}`; the operational detail requires the key. `/ingest` contains one
*example* path (`/data/xenium_run_2026_08/…`), not a real one.
`account.key_fingerprint` is a truncated HMAC-SHA256 under a fixed label, not a
digest of the key — non-reversing, and nothing authenticates against it.

---

## Residual risk this pass did not close

- **M1 and M2 are open** and both live in files outside this pass's ownership.
- The reference-panel allocation whose size is a function of the panel rather than
  the grid (so the grid ceiling structurally cannot catch it) was **not** examined
  here — another agent owns `tmesim/data/real/scref_gbm.py`.
- H1's caps bound *one* document. They do not bound the *number* of projects; a key
  holder can still create many. That is a quota question, not a validation one, and
  this deployment has no identity model to hang a quota on — which
  `GET /api/account` already says in its own payload
  (`"attribution": "none — … no run can be attributed to a person"`).
