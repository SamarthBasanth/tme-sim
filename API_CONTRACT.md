# tmesim Studio — API contract (v1)

Frozen interface between the backend (`app/server/`) and the front end (`app/web/`).
Both sides code to this. Do not change a route or payload shape without updating this file.

Stack: **FastAPI + uvicorn** backend, **vanilla ES modules + WebGL2** front end (no build
step, no npm — it must run from a clean checkout with only the existing venv plus fastapi).

Run: `./.venv/bin/python -m app.server` → serves the UI at `http://127.0.0.1:8000`.

---

## Design constraints (both sides)

The **visual** reference is diagnostic-imaging software — HALO, Visiopharm, Paige, Sectra —
**not** consumer web. That is a statement about layout and density and nothing else: Vivum is
a research tool, not a diagnostic product. Concretely:

- **Square corners** (≤2px radius). No pills, no large rounded cards.
- **Dense information.** Small type, tight leading, high data-ink ratio. A clinician's
  workstation shows everything at once; it does not hide data behind whitespace.
- **Monospace for every number**, `font-variant-numeric: tabular-nums`, always unit-labelled.
- **Restrained palette.** One accent. Semantic colour (ok / warn / critical) is separate from
  the accent and is used only for state, never decoration.
- **No** cream/terracotta, no gradient hero, no emoji as UI, no centred marketing layout,
  no Inter/Space Grotesk as the identity face.
- Every displayed number carries provenance: measured / computed / assumed.

---

## Data model

```jsonc
// Section — a reconstructable specimen
{ "id": "xenium_319", "label": "Xenium GBM · 319-panel",
  "n_cells": 40887, "extent_um": [8120.0, 5100.0], "panel": "hBrain_254g", "available": true }
```

```jsonc
// RunSpec — what the user submits
{ "section": "xenium_319",
  "dx_um": 12.0,
  "drug": "temozolomide",          // or null for reconstruction-only
  "sensitivity_mode": "uniform",   // "uniform" | "determinant"
  "vascular_mode": "discrete",     // "discrete" | "kde"
  "oxygen_mode": "pde",            // "pde" | "vessel_distance" | "markers"
  "state_mode": "oxygen",          // "oxygen" | "vessel_distance" | "markers" — cell-state route.
                                   // Mirrors oxygen_mode. Default "oxygen" is bit-identical to
                                   // the pre-option pipeline. Part of the reconstruction cache
                                   // key (it rewrites cells.state). sim.engine RE-BANDS states
                                   // from its own oxygen solve, so this moves the REPORTED band
                                   // occupancies, not kill_fraction. On the 319 panel NO route
                                   // clears the bar (results/cell_state_impact.md §5), so
                                   // absolute kill must not be quoted there — /api/validate
                                   // warns, and the run returns that caveat.
  "mgmt": null,                    // optional float, log2 CPM
  "pk": { "plasma_peak": 1.0, "half_life_s": 14400, "n_doses": 1, "dose_interval_s": 86400 },
  "total_time_s": 600.0,
  "resistance_off": false }         // optional, default false. The COUNTERFACTUAL arm of the
                                    // rank attribution: PDConfig(hypoxia_resistance=1,
                                    // quiescent_resistance=1), i.e. this tumour's oxygen-driven
                                    // resistance switched OFF. Everything else is bit-identical
                                    // to the same spec with it false, so
                                    // rho = kill / kill(resistance off) is a MEASURED ratio, not
                                    // an assertion. A true arm's kill_fraction is NOT a
                                    // prediction for the tissue and must not be quoted as one;
                                    // /api/validate warns and the run carries that caveat.
                                    // POST /api/reports generates both arms per drug itself.

  "attribute_mechanism": true,      // optional, default TRUE. Compute the per-cell mechanism
                                    // attribution (tmesim.sim.mechanism) and serve it on
                                    // GET /api/runs/{id}/mechanism. It is a read-out of arrays
                                    // cells.pkpd.apply_drug_kill already computed, takes no extra
                                    // RNG draw, and leaves the run BIT-IDENTICAL — which is why
                                    // it is on by default and why switching it off only removes
                                    // the answer to "why", never changes kill_fraction.
  "timelapse_frames": 24,           // optional, default 24, 0 disables, hard ceiling 40. Number of
                                    // drug-CONCENTRATION snapshots captured through the window and
                                    // served on GET /api/runs/{id}/timelapse. The payload is
                                    // n_frames full f32 grids, so it is reduced to whatever fits
                                    // the 48 MB budget and the run reports by how much — never
                                    // silently truncated, never silently skipped.
  "time_resolved": false }          // optional, default false (the historical Stage-0 path).
                                    // PDEConfig.time_resolved: drive the transport source with the
                                    // real plasma(t) curve instead of its window average, so
                                    // tissue concentration genuinely rises AND falls. Both paths
                                    // deliver the same TOTAL drug, so a difference is a SCHEDULE
                                    // effect, not a dose effect (results/dose_schedule.md). A
                                    // time-lapse on the static path shows the drug ENTER and
                                    // SPREAD toward a quasi-steady state but never CLEAR; only
                                    // this flag makes clearance visible. Every stored calibration
                                    // in the repo was produced on the static path.
```

```jsonc
// RunResult
{ "run_id": "r_00007", "spec": {...}, "status": "done",
  "report": { ... reconstruction_report ... },
  "metrics": { "kill_fraction": 0.115, "mean_tumor_auc": 1.2e7,
               "therapeutic_coverage_1x": 0.31, "penetration_length_um": 206.0 },
  "provenance": { "alpha": 148.6, "ic50_uM": 606.0, "ic50_source": "GDSC measured (n=51)",
                  "sensitivity_note": "...", "caveats": ["..."],
                  "state_mode": "oxygen",        // the cell-state route that was used
                  "cell_state": { ... },         // reconstruction_report()['cell_state'] verbatim:
                                                 // mode, source, is_measured, caveat, marker
                                                 // coverage detail, tumour band fractions
                  "state_unmeasured": [],        // states this panel cannot measure — assigned to
                                                 // NO cell; NOT "computed and found to be zero"
                  "state_drives_kill": false },  // sim.engine re-bands internally; see above
  "mechanism": { ... },            // the aggregate half of the attribution — the SAME JSON block
                                   // embedded in the TMSM payload, minus the per-cell labels.
                                   // null when there is none; `mechanism_reason` then says why.
  "mechanism_reason": null,
  "timelapse": { "available": true, "n_frames": 24, "interval_s": 939.1,
                 "realised_interval_s": 939.1, "window_s": 21600.0,
                 "frame_times_s": [...], "grid": [422, 674], "dx_um": 12.0,
                 "quantity": "drug concentration C (uM)", "payload_bytes": 27302880,
                 "plateau_note": "..." },   // frame metadata; the frames are binary
  "runtime_s": 41.2 }
```

---

## Routes

| Method | Path | Body → Returns |
|---|---|---|
| `GET`  | `/api/sections` | → `Section[]` |
| `GET`  | `/api/drugs` | → `[{name, mw, logp, ic50_uM, ic50_source, bbb_penetrance}]` |
| `GET`  | `/api/defaults` | → a valid `RunSpec` with current production defaults |
| `POST` | `/api/runs` | `RunSpec` → `{run_id}` — enqueues, returns immediately |
| `GET`  | `/api/runs` | → `RunResult[]` (summary; newest first) |
| `GET`  | `/api/runs/{id}` | → `RunResult` |
| `GET`  | `/api/runs/{id}/status` | → `{status, progress, stage, elapsed_s, log_tail}` — poll target |
| `DELETE` | `/api/runs/{id}` | cancel a queued/running job |
| `GET`  | `/api/runs/{id}/geometry` | → binary cell payload (see below) |
| `GET`  | `/api/runs/{id}/field/{name}` | `name` ∈ vascular_density, ecm_density, oxygen, drug_auc → field payload |
| `GET`  | `/api/runs/{id}/mechanism` | → binary TMSM: per-cell mechanism labels + the aggregate breakdown |
| `GET`  | `/api/runs/{id}/timelapse` | → binary TMSL: drug-concentration frames + their metadata |
| `POST` | `/api/validate` | `RunSpec` → `{ok, errors[], warnings[]}` — for the code editor |
| `POST` | `/api/sections/scan` | `{path, measure_detection_rates?, register?, section_id?, label?}` → `ScanResult` — ingest ANY directory |
| `DELETE` | `/api/sections/{id}` | remove an **ingested** section; **400** for a built-in one |
| `GET`  | `/api/sections/{id}/capability` | → `Capability` — what this panel can support, before a project exists |
| `GET`  | `/api/browse` | `?path=&q=&show_hidden=` → `{confined, places, path, crumbs, entries[], n_entries, truncated, diagnosis}` — the NAMES in one server-side directory; registers nothing |
| `GET`  | `/api/projects` | → `Project[]` (summary; newest first) |
| `POST` | `/api/projects` | `{name, notes}` → `Project` — **201** |
| `GET`  | `/api/projects/{id}` | → `Project` with its `Specimen[]` (capability refreshed on read) |
| `PATCH` | `/api/projects/{id}` | `{name?, notes?}` → `Project` |
| `DELETE` | `/api/projects/{id}` | → `{deleted, runs_retained[]}` — runs are **never** deleted |
| `POST` | `/api/projects/{id}/specimens` | `{section, label?, notes?}` → `Specimen` — **201**, **409** if that section is already a specimen |
| `GET`  | `/api/projects/{id}/specimens/{sid}` | → `Specimen` + its resolved `runs[]` |
| `PATCH` | `/api/projects/{id}/specimens/{sid}` | `{label?, notes?}` → `Specimen` |
| `DELETE` | `/api/projects/{id}/specimens/{sid}` | → `{removed, runs_retained[]}` |
| `POST` | `/api/projects/{id}/specimens/{sid}/runs` | `{run_id}` → `{attached}` — **201**; **409** if the run's `spec.section` ≠ the specimen's section |
| `DELETE` | `/api/projects/{id}/specimens/{sid}/runs/{run_id}` | → `{detached}` — the run itself is untouched |
| `GET`  | `/api/projects/{id}/compare` | `?metric=kill_fraction` → `Comparison` (below) |
| `POST` | `/api/reports` | `{section, drugs[], spec_overrides}` → `{report_id, status, run_ids[], warnings[]}` — **202**, enqueues, returns immediately |
| `GET`  | `/api/reports` | → `ReportSummary[]` (newest first) |
| `GET`  | `/api/reports/{id}` | → `Report` (the full JSON document; also the poll target) |
| `GET`  | `/api/reports/{id}/html` | → `text/html` — the rendered document, self-contained; **409** until every arm is `done` |
| `DELETE` | `/api/reports/{id}` | cancel a queued/running report (cancels every arm) |
| `POST` | `/api/protocol/validate` | `{source}` → `{ok, diagnostics[], n_runs, max_runs}` — parse only, no compute |
| `POST` | `/api/protocol/plan` | `{source}` → the expanded run list, run count and estimated cost; **executes nothing** |
| `POST` | `/api/protocol/run` | `{source}` → `{run_id, n_runs, warnings[], estimated_cost, truncation}` — **202**, enqueues the whole sweep as ONE job |
| `GET`  | `/api/protocol/runs/{id}` | → results + the sweep grid; safe to poll while running |
| `GET`  | `/api/protocol/sweeps/{id}` | → the same job in the editor's projection (alias) |
| `DELETE` | `/api/protocol/sweeps/{id}` | cancel a queued/running sweep |
| `POST` | `/api/protocol/upload` | `{filename, content \| content_base64}` or multipart `file` → `{path, preview}` |
| `GET`  | `/api/limits` | → what this server refuses; see below |
| `GET`  | `/api/account` | → this deployment: key, limits, usage inventory, billing |
| `POST` | `/api/preview` | `RunSpec` (+ optional `readout`) → ONE coarse readout with its measured agreement — **synchronous, and not a run** |
| `GET`  | `/api/preview/agreement` | → the error bar and the refusals, having solved nothing |
| `POST` | `/api/memo` | `RunSpec` or `{spec, allow_scaling?}` → `{mode, is_solve: false, …}` — **always 200**; a miss is an answer |
| `GET`  | `/api/auth/me` | → `{authenticated, via, user, …}` — **never 401**; the page asks this to decide what to draw |
| `POST` | `/api/auth/login` | `{email, password}` → sets the `vivum_session` cookie |
| `POST` | `/api/auth/logout` | → clears it; **200 even with no session** |
| `POST` | `/api/auth/signup` | `{email, password}` + the operator's bootstrap token → `201`, and does NOT sign you in |
| `POST` | `/api/auth/password` | `{old_password, new_password}` → session only; the key cannot change an account's password |

`status` ∈ `queued | running | done | failed | cancelled`. On `failed`, `RunResult.error`
carries the message; the UI must surface it, never swallow it.

### Pages — the routes that are not `/api/`

Registered on the app BEFORE the `StaticFiles` mount at `/`, because starlette
matches in registration order and a mount at `/` swallows every path under it.
The route table for the first three is `docs/WORKSPACE_CONTRACT.md` §1.

| Path | Serves | Notes |
|---|---|---|
| `/` | `app/web/home.html` | the workspace home. Falls back to the Studio when `home.html` is absent, rather than 500ing on the front door |
| `/studio` | `app/web/index.html` | the Studio. `/index.html` keeps working byte for byte, so every old bookmark survives |
| `/ingest` | `app/web/ingest.html` | bring-your-own-tissue. **The same page `/tissue` serves** |
| `/tissue` | `app/web/ingest.html` | the same page under a name that says what it is |
| `/login` | `app/web/login.html` | sign-in (`docs/AUTH_CONTRACT.md` §1) |
| `/healthz` | `{"ok": true}` | liveness, GET **and HEAD**. Point a platform health check HERE, not at `/api/health` |
| `/docs`, `/redoc`, `/openapi.json` | FastAPI's own | the map of the API. Ungated — see below |

**`/ingest` and `/tissue` are one page at two URLs, and the duplication is the
URL and not the page.** Until 2026-10-10 `/ingest` served a different,
self-contained form built as a Python string in `app/server/api.py`, whose only
input was a text box you typed an absolute server path into from memory;
`/tissue` served the directory browser that replaced it. Two pages did one job.
The form is deleted and both routes now read the one file
`app/web/ingest.html`, with the filename imported from `app/server/browse.py`
(`BROWSE_PAGE`) so the two URLs cannot drift apart. `/ingest` keeps answering
**200 with the page** rather than redirecting, because it is the path named in
the startup banner, `README.md`, `docs/QUICKSTART.md`, `docs/DEPLOYMENT.md` §0,
`docs/WORKSPACE_CONTRACT.md` §1 and the home page's own links — a 200 leaves all
of them literally correct, while a 3xx would only work for a client that
follows redirects. Each page route falls back to a served "not built" notice if
its file is missing, so a backend deployed ahead of the front end says what is
absent instead of 404ing a path the banner printed.

### Authentication — two credentials, one gate

`docs/AUTH_CONTRACT.md` IS THE SPECIFICATION; this is the part of it the front
end codes to. The gate is middleware in `create_app` and is installed **only
when `TMESIM_STUDIO_API_KEY` is set** — unset means no gate at all, which is the
local single-user default.

When it is installed it accepts **EITHER** credential:

* **`X-API-Key`** — the pre-existing shared secret. Every documented `curl`
  flow, `render.yaml` and CI depend on it, and nothing about that path changed.
* **the `vivum_session` cookie** — issued by `POST /api/auth/login`, HMAC-signed,
  14 days, no refresh-on-use. It identifies the person at the keyboard and
  **partitions nothing**: projects, runs and specimens stay shared, which is why
  `/api/auth/me` ships `shared_workspace: true` and a sentence saying so.

The key is tested first because it is the cheap test. A session fault of any
kind — no cookie, a tampered or expired one, a deleted or disabled account, an
unreadable store — reads as "no session" and falls through to the key, so the
gate fails closed exactly as it did before sessions existed.

Refused requests get **401** with
`{"detail": "missing or incorrect credentials: send an X-API-Key header, or
sign in at /login for a browser session"}` — the same body for a wrong key, a
truncated key and no key at all, so the key cannot be binary-searched a
character at a time. The front end keys off the status, not the sentence, and
sends the user to `/login`.

**The ungated paths, and nothing else:**

| Path | Why it must answer without a credential |
|---|---|
| `/healthz`, `/api/health` (GET/HEAD) | a platform health checker cannot send a header. `/healthz` returns a literal `{"ok": true}`; `/api/health` returns that much to an anonymous caller and the runs directory, run counts and cache statistics only to one holding the key — the ROUTE is exempt, the DATA is not |
| `/login`, `/api/auth/*` (any method) | a credential gate in front of the page that obtains the credential is a loop |
| every non-`/api/` GET/HEAD | the pages and static assets above. A page must be able to load in order to prompt for a key |

Two consequences stated so they are decisions rather than oversights:

* **Exempting the `/api/auth/*` PREFIX moves the burden into those handlers.**
  `POST /api/auth/logout` and `POST /api/auth/password` are not public
  operations; the middleware has simply stopped protecting them, and each
  re-establishes its own credential — the password route's 401 is its own.
* **Everything else fails closed.** Any `/api/` path, and any request whose
  method is not GET/HEAD/OPTIONS at any path, needs a credential. So a future
  non-`/api/` route that mutates state is authenticated by default rather than
  exposed by default. OPTIONS stays exempt because a CORS preflight cannot carry
  the header.

The gate decides on `get_route_path(request.scope)` — the same string the router
dispatches on — and never on `request.url.path`, which is rebuilt from the
client-supplied `Host` header and can be made to not start with `/api/` while
the router still dispatches the endpoint (GHSA-86qp-5c8j-p5mr, which was live
here).

`docs/DEPLOYMENT.md` §2 carries the operator's view of the same facts.

### The deployment's own surface — `GET /api/limits`, `GET /api/account`

Every other route in this contract answers a question about TISSUE. These two
answer questions about the INSTANCE, and they are the backing for Studio's
Account panel (`app/web/js/ui/account.js`, server side `app/server/account.py`).

Both are gated like any `/api/` route — either credential, see
**Authentication** above. That is deliberate: a key fingerprint, a
usage count and a deployment's refusal thresholds are operational detail about
a machine. The "how do I get a key" copy that an UNAUTHENTICATED visitor needs
therefore lives in the client, because by definition that visitor cannot call
these.

Both are `noFallback` on the front end: they never answer from `js/mock.js`.
A fabricated ceiling would make the spec editor confidently green on a spec the
real server refuses, and a fabricated billing block is a claim about money.

`GET /api/limits` →

```
{max_grid_cells, max_grid_cells_default, max_grid_cells_source,
 large_grid_cells, min_dx_um, dx_kill_ladder_coarse_um,
 dx_penetration_ladder_coarse_um, pde_dt_s,
 grid_cells_formula, field_bytes_per_grid_cell, grid_cells_requires}
```

`max_grid_cells_source` is `"default"` or `"TMESIM_MAX_GRID_CELLS"`.
`large_grid_cells` is the ALLOWED-BUT-SLOW threshold and is a different number
from the refusal ceiling; a client showing only the ceiling lets a user walk
into a ten-minute reconstruction with no warning.

`grid_cells_formula` is published as text because the client has to reimplement
it — it cannot call `/api/validate` on every keystroke — so there are two
copies of one formula by construction. `tests/test_app_server.py::
TestAccountAndLimits::test_limits_formula_matches_validate` checks the
published formula against the real refusal boundary on a real section, which is
the only thing that keeps the two honest.

It is a route of its own rather than a field on `/api/defaults` because the
`/api/defaults` body is adopted verbatim as a `RunSpec` and posted straight back
to `/api/validate`, which warns on unrecognised fields.

`GET /api/account` → `{contract, auth, deployment, usage, billing}`

* `auth` — `{required, header, scheme, source, key_fingerprint,
  caller_fingerprint, matches, self_service_signup, rotatable_via_api,
  how_to_rotate}`. The fingerprints are HMAC-SHA256 truncated to 12 hex, and
  are IDENTIFIERS: nothing authenticates against them, and the key itself never
  appears in the payload. `source` is `environment`, `--api-key`,
  `generated-at-startup`, or `unknown` — the distinction an operator who has
  LOST the key needs. Both `self_service_signup` and `rotatable_via_api` are
  `false`: there is no user table and no key-issuing endpoint.
* `deployment` — `{single_shared_key, data_roots_confined, limits,
  started_at_unix, uptime_s}`.
* `usage` — a **current inventory, not a meter**, and it says so in its own
  `kind` field. Deleting a run decrements it. `attribution` states the
  limitation outright: one shared key means no run can be attributed to a
  person.
* `billing` — scaffolding, labelled as scaffolding. `enabled` is always
  `false`; naming a provider in `TMESIM_BILLING_PROVIDER` moves `status` to
  `configured-but-not-implemented` and does NOT flip `enabled`. `blocked_on`
  names the preconditions (per-caller identity; a durable usage record) rather
  than leaving them for the UI to remember.

### Pruned payloads — `410 Gone`

A binary payload can be deleted to reclaim disk by
`scripts/prune_run_payloads.py`, which records the deletion in the run's own
`result.json` under `pruned_payloads`. `_binary` reads that record FIRST and
answers **410** with the truth — this run did produce it, it was removed on this
date, re-run the spec — instead of falling through to a reason that says the
payload was never produced. That fall-through sentence is false for a pruned
run, and a confident false reason is worse than the missing file it explains.

### TPL — `POST /api/protocol/*`

The Tumour Protocol Language (`PROTOCOL_LANGUAGE.md`, frozen). Backend:
`tmesim/protocol/`; editor: `app/web/js/ui/protocol.js`. Both code to that document.

**Every route takes the same body**, and a bare `text/plain` script is accepted too:

```jsonc
{ "source": "drug TMZ from smiles \"...\"\n...",
  "max_runs": 64,          // optional; the sweep budget for this request
  "max_window_s": 86400 }  // optional; longest transport window per run
```

**Diagnostic** — the only error shape these routes ever return. `line` and `col` are
1-based, `end_col` exclusive; `line: 0` means the message is not attributable to a
position (a backend fault). A user mistake is **never** a 500 and never a traceback.

```jsonc
{ "severity": "error",          // "error" | "warning"
  "line": 7, "col": 8, "end_col": 11,
  "code": "missing-unit",       // stable, keyable by the UI
  "message": "protocol stupp_adjuvant: dose: a unit is REQUIRED ...",
  "hint": "write `dose 150 mg/m2`" }
```

`POST /api/protocol/validate` → `{ok, diagnostics[], n_runs, max_runs}`. Parse, units,
name resolution, specimen capability and the run budget. It does **not** compile a
molecule or touch the solver, so it is fast enough for keystroke-level use. `ok: false`
still returns 200 — a rejected script is a normal answer, not a transport error.

`POST /api/protocol/plan` → everything `/validate` returns **plus**:

```jsonc
{ "runs": [ { "index": 0, "label": "TMZ/stupp_adjuvant@xenium_319 [dose=150 mg/m2]",
              "axes": {"dose": "150 mg/m2"}, "dx_um": 30.0,
              "spec": { ...RunSpec... },
              "dose_conversion": { "pk": {"plasma_peak", "half_life_s",
                                          "n_doses", "dose_interval_s"},
                                   "amount_mg", "bioavailability",
                                   "volume_of_distribution_L", "formula",
                                   "conversions": ["150 mg/m2 x 1.73 m2 = 259.5 mg", ...],
                                   "assumptions": [ {quantity, value, unit,
                                                     provenance, source} ],
                                   "infusion": null,
                                   "window_s", "protocol_duration_s",
                                   "n_administrations_outside_window",
                                   "time_resolved": true, "time_resolved_reason" },
              "reconstruction_key": ["xenium_319", 30.0, "discrete", "pde", "oxygen"] } ],
  "drugs": { "TMZ": { "identity": "temozolomide", "identity_route": "...",
                      "drug_params": {...},
                      "provenance": { "ic50": {"value", "unit",
                                               "provenance": "measured",
                                               "source": "GDSC glioma screen, ..."},
                                      "hill": {"provenance": "flagged-default", ...} } } },
  "capability": { "xenium_319:oxygen:uniform": { ...panel support report... } },
  "sweep_axes": [ {"param": "dose", "n_points": 4, "values": [...]} ],
  "expects": [ {"metric": "kill", "comparator": ">", "value": 0.05} ],
  "estimated_cost": { "n_runs", "n_reconstructions", "reconstructions_cached",
                      "total_pde_timesteps", "approx_grid_cells",
                      "wall_time_s": null, "wall_time_reason": "..." } }
```

`provenance` is per **field**, and the tag is one of `measured` / `computed` /
`user-supplied` / `flagged-default`. `estimated_cost.wall_time_s` is deliberately
`null` with a reason — the backend does not benchmark in advance and will not guess.

`POST /api/protocol/run` — **202**, or **400** with
`detail: {message, ok: false, diagnostics[]}` for a rejected script. A whole sweep is
**one** job on the same `JobManager` and the same `RECON_CACHE` as `/api/runs`, so the
runs are serialised on one worker and a dose sweep on one specimen **reconstructs once**.

`GET /api/protocol/runs/{id}` → `{run_id, status, stage, progress, elapsed_s, error,
source, n_runs, grid[], rows[], sweep_axes[], expects[], diagnostics[],
reconstructions, plan, truncation, log_tail[]}`. **409** if `{id}` is a `RunSpec` run
rather than a protocol run.

### `truncation` — the window the solve actually covered

`{truncated, n_runs, n_runs_truncated, window_s, worst_run_index, message, runs[]}`,
carried on `/plan`, `/run`, and on `/runs/{id}` and `/sweeps/{id}` **while the job is
still queued**. Verified against a live plan rather than transcribed:

```
run TMZ with stupp on xenium_319   (75 mg/m2 every 24 h, repeat 30, dx 30 um)
  -> truncated: True   n_runs_truncated: 1/1   window_s: 86400.0
     "This protocol spans 30 d but the transport solve covers the first 24 h.
      29 of 30 administrations fall OUTSIDE the observed window and deliver no
      drug in this result. The window is a compute limit ... not a claim that
      the later doses do not matter."
```

A transport solve covers at most `TMESIM_PROTOCOL_MAX_WINDOW_S` (default 86,400 s). A
longer protocol is truncated to it — a Stupp course is ~168 days, and at `dt = 10 s` that
is 1.45e6 steps, which this backend will not start.

`tmesim/protocol/dosing.py` has always said the truncation "is reported (with how many
doses fell outside it) rather than hidden". Traced end to end on 2026-09-03 it reached
`/plan`'s `assumptions[]` and **nothing else**: `/run` — the commit point, where a user
spends the compute — said nothing, because the sentence lived on a `DoseConversion` and
`/run` plans with `compile_drugs=False`. An operator going validate → run was never told
that **29 of a Stupp course's 30 administrations fall outside the observed window**. The
key above is that fix; the arithmetic now lives in one place (`window_plan()`) so the
number and the sentence cannot drift apart.

It is REPORTED, not refused. The contract's posture here is report-not-hide, and a
mandatory acknowledgement would break existing clients.

> **2026-08-21 — `ranking` removed from this payload.** It carried the `rank by` order of a
> multi-compound `compare`. Both are withdrawn from TPL (`PROTOCOL_LANGUAGE.md` §5): a
> `compare` now names ONE compound and its runs come from a `sweep`, so there is no
> cross-compound order to return. `rows[]` still carries every run's metrics, and the
> defensibility disclosure that used to ride on `ranking.note` is now in `diagnostics[]`,
> emitted per run whether or not a `compare` was written. Each entry of `plan.runs[]` also
> loses its `rank` field for the same reason.

`grid` is the progress grid and fills in live:

```jsonc
{ "index": 0, "label": "...", "axes": {"dose": "150 mg/m2"},
  "status": "queued|running|done|failed|cancelled",
  "metrics": {"kill_fraction", "therapeutic_coverage_1x", "mean_tumor_auc",
              "penetration_length_um"},
  "error": null, "runtime_s": 41.2, "reconstruction_cache_hit": true }
```

`rows[]` carries the full result per run: `metrics`, `reported` (just the metrics the
script asked for), `extra_metrics`, `expectations` (each `expect` scored pass/fail —
**reported, never enforced**), `report`, and `provenance` with `drug`, `dose_to_pk`,
`sensitivity`, `capability`, `cell_state`, `solver` and `caveats`.
`reconstructions` reports `{built[], n_built, n_cache_hits}` so the caching claim is
checkable from the response.

`POST /api/protocol/upload` → `{ok, path, bytes, use_in_script, preview}`. `path` is
what a `drug X from file "<path>"` statement should use, verbatim. Structures
(`.smi`, `.mol`, `.sdf`) preview as `{kind: "structure", smiles, descriptors}`;
parameter files (`.json`, `.yaml`) as `{kind: "parameters", parameters}`. The body is
JSON rather than multipart so no new dependency enters this checkout. Filenames are
basenamed and sanitised; anything that will not parse is **400** and is not kept.

#### Two projections in one payload

`app/web/js/protocol/api.js` froze its own wire shape for these routes before this
backend landed, and the two do not agree: it wants `errors[]` / `warnings[]` with
**0-based, half-open** columns, `run_count`, `estimated_cost_s`, `axes`, `refusal`, and a
sweep under `/api/protocol/sweeps/{id}`. Rather than break one client, **every response
carries both projections under disjoint keys**, and `app/server/protocol_wire.py` owns the
second one:

| this contract | the editor's | note |
|---|---|---|
| `diagnostics[]` — `col`/`end_col` **1-based** | `errors[]` + `warnings[]` — `col`/`end_col` **0-based, half-open** | the two conventions never share a key |
| `n_runs` | `run_count` | same number |
| `estimated_cost.wall_time_s` = `null` | `estimated_cost_s` + `cost_note` | the editor's is an ESTIMATE from throughput constants measured once on this repo; the contract's stays null on purpose |
| `sweep_axes[]` (`values` as strings) | `axes[]` (`values` as numbers) + `runs[].coords[]` | `coords` are indices into `axes[i].values` |
| the first budget error in `diagnostics[]` | `refusal: {code, line, message, remedy, remedy_text, remedy_line}` | `remedy_text` is TPL to insert verbatim, or `null` |
| `GET /api/protocol/runs/{id}` | `GET /api/protocol/sweeps/{id}` | the same job; `/runs/{id}` returns both |

`POST /api/protocol/run` also accepts `confirm_run_count` and answers **409** if the
program now expands to a different number — plan-before-run means the operator confirmed
a *specific* count. `POST /api/protocol/upload` accepts `multipart/form-data` with a
`file` part as well as JSON; multipart is parsed in-process rather than by
`python-multipart`, which is not in this checkout and must not be added.

`runs[].run_id` in a sweep payload is always `null`, with `run_id_reason` alongside: a TPL
sweep is **one** job whose points share a cached reconstruction and must be serialised on
the single worker, so there is no per-point run to open in the viewer. The point's
metrics, provenance and caveats travel in the sweep payload instead.

**TPL runs always set `PDEConfig.time_resolved = True`** and say so in
`provenance.dose_to_pk.time_resolved_reason`. With the engine default (`False`) the
solver replaces the plasma curve by its window mean, which makes `every`, `repeat`,
`cycles` and `over` bit-identical inputs — the artefact `results/dose_schedule.md` was
regenerated to remove. This is the one solver setting the protocol path changes
relative to `/api/runs`.

### The fast preview — `POST /api/preview`, `GET /api/preview/agreement`

**A coarser solve shown fast, and not a model.** `POST /api/preview` runs the
SAME solver on `dx = 30 um` instead of the shipped `dx = 12 um`, synchronously,
and returns ONE readout: the delivery spread `E_p90/E_p10` over the tumour cell
set. Median wall clock 8.6678 s against 58.7711 s for the real run. It is
synchronous on purpose — a preview with its own queue, poll and cancel would be
a second run-shaped object a client could mistake for a run.

Body is a `RunSpec`, plus an optional `readout`. `dx_um` is **required and is
the dx the REAL run will use**, because the payload has to name which measured
pair it is quoting.

```jsonc
// 200
{ "is_preview": true, "quotable": false, "preview_label": "PREVIEW - A COARSER SOLVE, …",
  "readout": "E_p90_over_p10", "readout_name": "delivery spread  P90 / P10",
  "value": 3.41, "unit": "x", "null_reason": null, "n_tumor": 12043,
  "derivation": "extra_metrics.E_p90_uM / extra_metrics.E_p10_uM from a dx = 30 um solve …",
  "agreement": { /* the measured dx = 30 vs dx = 12 figures */ },
  "preview_dx_um": 30.0, "target_dx_um": 12.0, "preview_grid": [170, 271],
  "preview_wall_clock_s": 8.67, "section": "xenium_319", "drug": "temozolomide",
  "refused": [ {"readout": "…", "why": "…"} ],
  "caveats": [ /* the engine's own, verbatim */ ], "claim_set": { /* … */ },
  "superseded_by": "the real run. …" }
```

**Nothing is enqueued, nothing is persisted, no `run_id` is minted** — the
module calls `execute_run` directly and drops `_payloads`. So a preview number
cannot be opened in the Read rail, attached to a Specimen, exported, put in a
report or counted on a project card, which are the four routes by which a fast
number becomes a quoted one. `quotable: false` and the `agreement` block are
both checked on the FINISHED body before it is returned; a body that has lost
either is not served.

Every refusal is a measurement, not a preference, and each carries its own
figure:

| Refused | The measurement |
|---|---|
| spot (Visium) sections | `dx_um` is not a lever there at all — a preview at 30 and a run at 12 are both clamped to 100 um and are the IDENTICAL solve at the IDENTICAL cost |
| `penetration_length_um` | disagrees by up to 16.0405 % (17.3158 um) over the 12 pairs |
| `penetration_index` | has NO settled `dx = 30` / `dx = 12` agreement: this repository's own three measurements of nominally the same section are 1.04 %, 26.0 % and 102 % apart, in two directions |
| any absolute level | `median_auc_over_peak` agrees to 0.3215 % — the best of any readout — and is still refused, because agreement is not defensibility |
| anything whose name contains `kill` | never a defensible output here, at any resolution |
| a pair nobody measured | the envelope is `dx = 30` vs `dx = 12`, two Xenium sections, six compounds — 12 pairs. Outside it this refuses rather than widening the tolerance |

Refusals come back as `{"error": "<the sentence>", …}` with the fields naming
what was asked (`readout`, `field`, `dx_um`, `previewable`); the default status
is **422**, with **409** for a run already at `dx = 30` and **503** when the
measured agreement is not on this server at all. That last one is the design:
the figures are read at request time from
`experiments/exact_memo/coarse_baseline.json` and are deliberately NOT copied
into the source, so **if the measurement is gone the feature is gone** rather
than the error bar being gone.

`GET /api/preview/agreement` serves the same error bar and the same refusal
list having solved nothing, so a surface can show what a preview would and
would not tell you before spending the seconds — and so the ABSENCE of a
measurement is readable rather than only discoverable by being refused:

```
{available, agreement, unavailable_reason, preview_dx_um, target_dx_um,
 previewable, preview_label, measured_sections, measured_compounds,
 spot_refusal, refused[]}
```

### The exact memo — `POST /api/memo`

**Answer a `RunSpec` without solving it, in the two cases where that is exact.**
There is no RunSpec-keyed run cache in this backend, so re-submitting a solved
spec re-solves it. This route answers from the store instead, and it is not a
model:

* **`exact_memo`** — the stored record of a run whose RunSpec, server-side
  solver constants, engine fingerprint and claim-set version all match. The
  numbers are bit-identical to a re-solve.
* **`algebraic_scale`** — another run's numbers rearranged by an exact algebraic
  identity for a change of dose and/or half-life. In the shipped
  (`time_resolved: false`) path the whole PK block reaches the solver as ONE
  float, `plasma_mean`, and the transport chain is linear in it with `C(0) = 0`,
  so every exposure ratio is invariant and every magnitude is proportional. The
  payload states where it is bit-identical (a dose-only change by a power of
  two) and where it agrees only to a measured number of ULPs.
* **`miss`** — nothing answers this spec exactly and **nothing was
  approximated**.

Body is a bare `RunSpec`, or `{"spec": <RunSpec>, "allow_scaling": <bool>}` — a
RunSpec has no `spec` key, so the two forms are unambiguous. **`allow_scaling`
defaults to `false`.** A caller that asked whether the answer is already known
and was handed an algebraic rearrangement of a DIFFERENT run would have been
told something it did not ask about, and that is the one distinction this
surface must never blur.

**Always 200.** A miss is an answer ("nothing here answers this; submit it"),
not an error, and every shape of non-answer carries the same
`mode` / `reason` / `disclosure` fields so a client has one branch to write. The
spec goes through `/api/validate`'s own validator first, so this surface can
never answer a spec the run surface would reject; a rejected spec comes back as
a miss carrying the full `validate` verdict.

```jsonc
{ "mode": "exact_memo",            // exact_memo | algebraic_scale | miss
  "is_solve": false,               // always, in every mode
  "scaling_requested": false,
  "spec": { /* … */ }, "key": "…", "shape_key": "…", "key_components": { /* … */ },
  "disclosure": "NOT A SOLVE: this is the stored record of run r_00007, …",
  "source_run_id": "r_00007",      // absent on a miss
  "validate_warnings": [] }
```

`is_solve: false` and `disclosure` exist because a memo hit and a fresh solve
agree on every NUMBER and do **not** agree on provenance: a hit carries the
original execution's `run_id`, its timestamps and its payload endpoints. A
client that recorded a hit as a new run would mint a second citation to one
execution. An `algebraic_scale` additionally carries `derivation` (the two
plasma means, their ratio, `bit_exact`, the identity in words, and the
`measured_by` script), `readouts`, `copied_verbatim` and **`refused`** — the
readouts the identity does not cover, which were NOT estimated: anything
compared against a threshold (`therapeutic_coverage_*`), anything through the
Hill kill model, and the exponential-plus-floor FIT (`penetration_length_um`).

This route deliberately does **not** hook `POST /api/runs`. A memo on the
submit path would make a second click return the first click's `run_id`.

### Generic sample ingestion — `POST /api/sections/scan`

The seven built-in sections are hardcoded rows. This route is how a **new specimen** —
one no table anywhere mentions — becomes a section, without editing code.
`tmesim/data/ingest.py` does the work: it decides what a directory is, loads it through
the **existing** loaders wherever the layout matches (so a built-in section ingests to
bit-identical arrays), validates it, and reports what it can and cannot support.

**Layouts recognised.** Xenium output bundles (`cells.csv.gz` + `cell_feature_matrix/`
or `*cell_feature_matrix.tar.gz` + `gene_panel.json`); Space Ranger / Visium outputs
(`spatial/` with `tissue_positions[_list].csv` + `scalefactors_json.json`, plus
`filtered_feature_bc_matrix/`); and plain MEX triples (`matrix.mtx` + `features.tsv` +
`barcodes.tsv`, gzipped or not) beside a coordinates table. A directory that only
*contains* samples is reported as `platform: "multi_sample"` with its children listed —
never guessed at.

```jsonc
// ScanResult (register: false — the default, a read-only scan)
{ "detection": { "path": "/data/study7_s042", "platform": "xenium",
                 "layout": "xenium_bundle_extracted", "confidence": 1.0,
                 "issues": [], "notes": ["…"], "contained_samples": [],
                 "ingestable": true },
  "files": { "cells": "…/cells.csv.gz", "matrix_dir": "…", "panel": "…" },
  "capability": { /* CapabilityReport, below */ },
  "capability_text": "==== CAPABILITY REPORT …",   // the same thing, printable
  "sample": null,                                   // populated on a deep scan
  "error": null }
```

`measure_detection_rates` (default `false`) is the cheap/expensive switch, and it changes
what the report is allowed to *claim*:

- **`false`** — reads filenames and `features.tsv.gz` only; milliseconds even on the
  851k-cell section. Reports **panel PRESENCE**. Every rate-dependent verdict comes back
  `available: null` with `"UNVERIFIED"` in its reason. It does **not** guess.
- **`true`** — performs a full ingest and measures the **detection RATE** of every marker
  set from the count matrix. This is what turns *"MKI67 is on the panel"* into
  *"MKI67+TOP2A fire in 1.34% of cells, below the 5% floor"*.

`register: true` implies a deep scan and adds the sample to the section registry, returning
`{ok, section_id, section, capability, capability_text, error}`. Registration is **refused**
(`ok: false`, nothing stored) when the directory is not ingestable, fails validation, or the
id collides — and a **built-in id can never be shadowed**, so nothing an ingested section
does can change what a published number was measured from. Ingested sections appear in
`GET /api/sections` *after* the seven built-ins with the extra fields `ingested: true`,
`platform`, `layout`, `detect_confidence`, `detect_issues[]`, `runnable`, `capability`,
`validation`; they are otherwise ordinary sections and go through the same `/api/validate`
and run path. They persist in `app/runs/_ingested_sections.json` **by path** — if the bytes
go away the section reports `available: false` with the reason, exactly like a built-in one
with a missing cache.

`platform: "mex"` samples register as `runnable: false`: they are ingestable and their
capability report is real, but `reconstruct_real` accepts `xenium` or `visium` only.

> `path` is a **server-side filesystem path** — nothing is uploaded, the bytes are read
> where they lie — and this route reads **any path the server process can read**. The
> server binds to `127.0.0.1` by default; do not expose this route without the credential gate
> and a reverse proxy in front (`docs/DEPLOYMENT.md` §1b, §3). Note also that a registered
> path is echoed back to every client as `cache_dir` on `GET /api/sections`, so the path,
> the `section_id` and the `label` must all be study codes rather than anything that
> identifies a person (`docs/DEPLOYMENT.md` §0) — hence the naming in the examples here.

```jsonc
// CapabilityReport — what this specimen can and cannot support, BEFORE a run
{ "sample_id": "study7_s042", "platform": "xenium", "n_genes": 319,
  "detection_rates_measured": true,
  "panel": { "selected_map": "brain_319", "selected_coverage": 1.0, "fits": true,
             "coverage": { "brain_319": 1.0, "io_480": 0.5 }, "verdict": "…" },
  "cell_types":   { "TUMOR": { "markers": [...], "present": [...], "missing": [],
                               "detection_rate": 0.7204 }, /* … */ },
  "proliferation":{ "present": ["MKI67","TOP2A","PCNA"], "call_set": ["MKI67","TOP2A"],
                    "detection_rate": 0.0134, "floor": 0.05,
                    "scoring_population": "all cells" },
  "hypoxia":      { "present": [], "missing": ["VEGFA","CA9","ADM","NDRG1","SLC2A1"] },
  "mgmt": {...}, "ecm": {...}, "endothelial": {...},
  "supports": {
    "hypoxia_measured":  { "available": false, "reason": "hypoxia UNMEASURED: …",
                           "evidence": "panel presence" },
    "state_mode:markers":{ "available": false,
                           "reason": "proliferation detection 1.34% is below the 5.0% …",
                           "evidence": "measured detection rate 1.3354%" },
    /* oxygen_mode:markers, sensitivity_mode:determinant, cell_typing,
       vascular_field, ecm_field, mgmt_from_panel, absolute_kill */ },
  "headline": [ { "capability": "…", "status": "YES|NO|CAVEAT", "reason": "…",
                  "evidence": "…" } ] }
```

`available` is deliberately **three-valued**: `true` supported, `false` refused with the
reason, `null` *not knowable from what was read* (a caveat, never silently a yes). The
flagship 319-gene panel comes back `false` on `hypoxia_measured`, `state_mode:markers`,
`oxygen_mode:markers` and `absolute_kill` — which is `results/cell_state_impact.md` §1.1
and §5, stated before a run instead of discovered after quoting a number.

One honest asymmetry: rates here are over **all cells**, because phenotyping has not run
yet at ingest time, while `reconstruct.cell_state` re-measures over **TUMOR cells only**.
On the 319 panel that is 1.34% here vs 1.42% at run time — same verdict, different number,
and both are labelled.

### Finding the tissue — `GET /api/browse`

The other half of the same capability. A section registers **by filesystem
path**, so the input to `POST /api/sections/scan` is an absolute path on the
machine running this server — and until this route existed the interaction was
typing it from memory into a text box. This lets the operator NAVIGATE the
server's directories instead. **It adds no capability:** registration is still
`POST /api/sections/scan`, and this adds a way to FIND a directory and not one
new thing you can do to one.

It cannot be a file picker, and that is what the product does rather than a
limitation of this route. An `<input type="file">` cannot express "that
directory over there on the host" at all: it would either be a lie about what
happens or a multi-gigabyte transfer of data that is already on the right
machine.

```jsonc
// GET /api/browse?path=/data/study7&q=xen&show_hidden=0
{ "confined": true,                      // is TMESIM_STUDIO_DATA_ROOTS set?
  "places": [ {"name": "data_cache", "path": "/app/data_cache",
               "why": "allowed by TMESIM_STUDIO_DATA_ROOTS"} ],
  "path": "/data/study7",                // the ONE absolute path in the body, resolved
  "crumbs": [ {"name": "data", "path": "/data"} ],
  "entries": [ {"name": "s042", "kind": "dir"}, {"name": "notes.txt", "kind": "file"} ],
  "n_entries": 2, "truncated": false,
  "diagnosis": { /* the registry's own `detect` verdict, paths stripped */ } }
```

* **An entry is `{name, kind}` and nothing else** — no size, no timestamp, no
  owner, no path. A size column would make this a read primitive over
  directories it is only meant to let you walk. `kind` is `dir`, `file`, `link`
  or `unknown`; **a symlink is never followed** to compute a type, because
  `stat`ing the target of a link pointing out of the allowlist is precisely the
  disclosure the allowlist exists to prevent. Clicking one issues an ordinary
  request and is refused if it leaves.
* **Entries carry no paths.** The caller gets `path` (the one directory it asked
  about, resolved) and child NAMES, and composes the next request itself.
* **No `path` returns the navigation anchors and no listing** — the allowlisted
  roots when `TMESIM_STUDIO_DATA_ROOTS` is set, otherwise a handful of plausible
  starting points. `confined: true` with an EMPTY `places` is the
  set-but-usable-for-nothing state and refuses everything by design: an operator
  who set the variable asked for confinement, and a typo must not hand back the
  behaviour they were switching off.
* `q` is a case-insensitive substring filter applied **before** the 500-entry
  cap, so a big directory stays navigable rather than merely truncated. Over the
  cap the body says `truncated: true` and `diagnosis` is **omitted** — the
  registry's `detect` walks every child, and a directory too big to list is too
  big to walk.
* `show_hidden` is off by default (`1`/`true`/`yes`/`on` to enable).

**Confinement is `TMESIM_STUDIO_DATA_ROOTS` through the same
`cache.confine_scan_path` the scan route uses, and there is no second
mechanism.** The check is on the `realpath`, so `..`, a leading `~` and every
symlink are resolved BEFORE the decision and the directory opened is the
directory that was allowed. The allowed roots ARE disclosed, as navigation
anchors, to a caller that already holds the credential — you cannot navigate to
a place you are not told exists — while the 403 refusal still never names the
allowlist.

Refusals, all in the `{code, message}` shape the scan route already uses:
**403** `path_outside_data_roots` (naming the caller's own text, the environment
variable and the way out) and `directory_unreadable`; **400** `path_not_a_path`;
**404** `no_such_directory` and `not_a_directory`. A 404 echoes the caller's own
text and never the resolved path. When `TMESIM_STUDIO_DATA_ROOTS` is unset this
route warns **once per process** on stderr, separately from the scan route's
warning, because an operator told "scan is unconfined" has not been told that a
directory browser is unconfined too.

The page that drives this route is `/ingest` and `/tissue` (see **Pages**).

### Projects, specimens, and the cross-specimen comparison

A **Project** is a study; a **Specimen** is a patient, bound to one `section` id. Specimens hold
`run_ids` **into `app/runs/`** — runs are *associated, never duplicated*, and neither deleting a
project nor detaching a run touches the run itself. Projects persist as one JSON document per
project under `app/projects/<id>/project.json` (gitignored, like `app/runs/`).

```jsonc
// Project
{ "id": "p_00001", "name": "Two-patient GBM therapy panel", "notes": "…",
  "created_at": 1.7e9, "updated_at": 1.7e9,
  "specimens": [ /* Specimen */ ] }
```
```jsonc
// Specimen — a section id, its capability report, and the runs that belong to it
{ "id": "s_00001", "label": "Patient A · 319 panel", "section": "xenium_319",
  "notes": "", "created_at": 1.7e9, "run_ids": ["r_00007", "r_00008"],
  "capability": { /* Capability */ } }
```
```jsonc
// Capability — what this specimen's PANEL can support. Derived from the panel gene
// list (cheap: the memoised features.tsv.gz symbols, never the count matrix).
{ "section": "xenium_319", "panel": "hBrain_254g + hBrain_65g", "source": "xenium",
  "n_cells": 40887, "n_genes": 319, "available": true,
  "genes_known": true,
  "proliferation_present": ["MKI67","TOP2A"], "cycle_specific_present": ["MKI67","TOP2A"],
  "hypoxia_present": [], "hypoxia_missing": ["VEGFA","HIF1A","CA9","ADM","NDRG1","SLC2A1"],
  "state_routes": {
    "oxygen":          { "supported": true,  "is_measured": false, "reason": "…" },
    "vessel_distance": { "supported": true,  "is_measured": false, "reason": "…" },
    "markers":         { "supported": true,  "is_measured": true,  "reason": "…" } },
  "supports_cell_state": false,          // panel measures hypoxia AND the marker route runs
  "absolute_kill_supported": false,      // ← the load-bearing verdict
  "absolute_kill_reason": "…" }
```

`absolute_kill_supported` is `false` unless the panel carries at least one hypoxia-response
transcript. That is `results/cell_state_impact.md` §5: with no hypoxia gene, the `oxygen`
default's HYPOXIC band is empty for a *bookkeeping* reason, `markers` declares hypoxia
UNMEASURED, and `vessel_distance` is a geometric proxy — so **no** route clears the bar, and
kill moves by up to **0.3252 absolute (1.77×)** purely with the choice of route (§3.2). The
319-gene flagship panel fails this test; the 480-gene IO panel passes it (VEGFA).

```jsonc
// Comparison — GET /api/projects/{id}/compare
{ "project_id": "p_00001", "metric": "kill_fraction",
  "specimens": [ { "id": "s_00001", "label": "…", "section": "…",
                   "capability": {…}, "absolute_kill_supported": false,
                   "unsupported_state_routes": [],
                   "labels": ["ABSOLUTE KILL NOT SUPPORTED BY THIS PANEL"],
                   "n_runs_attached": 3, "n_runs_done": 3,
                   "runs_not_done": [], "runs_missing": [] } ],
  "matched_sets": [ { "signature": "dx_um=30.0; sensitivity_mode=uniform; …",
                      "is_primary": true, "specimens": ["s_00001","s_00002"],
                      "drugs": ["cisplatin","doxorubicin","temozolomide"],
                      "n_matched_drugs": 3 } ],

  // 1. THE PRIMARY RESULT.
  "ranking": { "metric": "kill_fraction", "drugs": ["cisplatin","doxorubicin","temozolomide"],
               "by_specimen": { "s_00001": { "order": ["cisplatin","doxorubicin","temozolomide"],
                                             "values": { "cisplatin": 0.37, … } } },
               "agrees": true,
               "pairwise": [ { "a": "s_00001", "b": "s_00002", "agrees": true,
                               "kendall_tau_b": 1.0, "swapped": [] } ],
               "verdict": "RANKING AGREES across every specimen pair. …" },

  // 2. Absolute values — readable ONLY where `comparable` is true.
  "kill": [ { "drug": "temozolomide", "reference": "s_00001",
              "runs": { "s_00001": "r_00007", "s_00002": "r_00010" },
              "metrics": { "kill_fraction": { "values": {…}, "delta": {…}, "ratio": {…} }, … },
              "comparable": false,
              "blockers": ["Patient A (xenium_319): …carries none of the hypoxia-response transcripts…"],
              "caveats": ["…different panels…", "…uniform sensitivity is an upper bound…"],
              "state_measured": { "s_00001": false, "s_00002": false } } ],

  // 3. Composition — cell-type fractions, vascular statistics, measured state fractions.
  "composition": { "from_run": { "s_00001": "r_00007" },
                   "rows": [ { "key": "frac_endothelial", "label": "ENDOTHELIAL fraction",
                               "unit": "frac", "fmt": "pct",
                               "values": { "s_00001": 0.0974, "s_00002": 0.0418 },
                               "provenance": { "s_00001": { "kind": "computed", "why": "…" } },
                               "comparable": false, "reason": "…",
                               "reference": "s_00001", "delta": {…}, "ratio": {…} } ] },
  "headline": "RANKING AGREES … ABSOLUTE KILLS ARE FLAGGED NOT-COMPARABLE …",
  "problems": [] }
```

Rules the comparison enforces, and why:

- **Matching.** Two runs are compared only when their RunSpecs are identical apart from
  `section` and `drug` (`spec_signature` — mesh, all four modes, `mgmt`, the whole `pk` block
  and `total_time_s`). Different mesh or different state route ⇒ different matched set.
- **Rank agreement is the primary result.** `results/heldout_patient.md` §3.3 found the kill
  ordering IDENTICAL on two real GBM sections; `results/cell_state_impact.md` §3.3 found it
  INVARIANT across all three cell-state routes. It is restricted to the drugs *every* compared
  specimen has a matched run for. Disagreements list the swapped pairs.
- **`comparable: false` for absolute kill** unless *every* specimen's panel has
  `absolute_kill_supported` **and** every matched run actually used a MEASURED state route. A
  measured kill is never silently differenced against an assumed one.
- **Specimens whose panel cannot support a route are labelled**, from the capability report, in
  `specimens[].labels` and `unsupported_state_routes` — the comparison never compares a measured
  quantity against an assumed one without saying so.
- **Composition rows carry their own verdict.** `mean_ecm_density` is *always*
  `comparable: false` (a max-normalised KDE over whichever matrix transcripts the panel carries
  — 11 genes vs 3 — has no common absolute scale, `heldout_patient.md` §3.5). `mean_oxygen` is
  *always* `comparable: false` (uncalibrated on every section). Cell-type fractions are
  `comparable: false` across different panels (§1.2). The ENDOTHELIAL *fraction* and *areal
  density* rows sit adjacent because on those two sections they move in opposite directions
  (0.43× and 2.67×, §1.1) and quoting one alone inverts the conclusion.

### Binary payloads (performance-critical)

`/geometry` returns `application/octet-stream`:
```
magic "TMSG" (4B) | version u32 | n_cells u32 | extent_x f32 | extent_y f32
| x f32[n] | y f32[n] | type u8[n]
```
`/field/{name}` returns:
```
magic "TMSF" (4B) | version u32 | rows u32 | cols u32 | dx_um f32
| min f32 | max f32 | values f32[rows*cols]
```
`/mechanism` returns:
```
magic "TMSM" (4B) | version u32 | n_cells u32 | n_labels u32 | json_len u32
| aggregate json[json_len] (UTF-8) | label u8[n_cells]
```
One byte per cell, in **exactly the geometry payload's cell order**, so the front end uploads
it as a third vertex attribute against the position buffer it already holds. The aggregate
breakdown travels in the same response as a short JSON block, because a legend without its
counts is decoration and a second round trip for four integers is not worth a request.

`/timelapse` returns:
```
magic "TMSL" (4B) | version u32 | n_frames u32 | interval_s f32 | window_s f32
| n_frames x [ frame_index u32 | time_s f32 | <a complete TMSF block> ]
```
The per-frame payload is **byte-identical to what `/field/{name}` serves** — the existing TMSF
shape prefixed with a frame index and a time stamp, not a new grid format. The front end decodes
a frame with the TMSF reader it already ships, and each frame carries its own finite min/max so a
scrubber can window per-frame or across the series without a second pass. Frames are the drug
**concentration** field (uM), which is a different quantity from `/field/drug_auc`: AUC only ever
accumulates, so it can never show drug leaving.

All little-endian. The front end reads these into typed arrays and uploads straight to WebGL.
JSON for 850k cells is not acceptable.

---

## Mechanism attribution

`kill_fraction` says how much happened. It cannot say **why**, and the kill model already knows:
`cells.pkpd.apply_drug_kill` computes a local exposure, a Hill term, a resistance multiplier, a
kill probability and a Bernoulli draw for every cell, and discards all five one line later.
`tmesim.sim.mechanism` picks them back up. It is a **read-out of the model's own arithmetic, not
a second model** — `apply_drug_kill` fills a `diagnostics` dict with the exact arrays it decided
on (the same objects, not recomputed copies) and every label is a comparison between two of them.

Per-cell labels, and the `u8` code each occupies on the wire:

| code | label | limiter | meaning |
|---|---|---|---|
| 0 | `KILLED` | — | the engine's own `killed` mask, not re-derived |
| 1 | `SURVIVED_UNREACHED` | **delivery** | `E/IC50 < unreached_ratio` — the drug never arrived |
| 2 | `SURVIVED_UNDERDOSED` | **dose** | drug arrived, `E` stayed below `IC50` |
| 3 | `SURVIVED_RESISTANT` | **resistance** | exposure adequate, `sensitivity x state <= resistant_modulation` |
| 4 | `SURVIVED_STOCHASTIC` | **chance** | both adequate; the cell drew above its `p_kill` |
| 5 | `UNATTRIBUTED` | — | not explained by the bands. **Counted, never reassigned.** |
| 6 | `ALREADY_DEAD` | — | tumour cell dead before the drug; counts in `kill_fraction`, not killed by the drug |
| 7 | `NOT_TUMOR` | — | non-tumour cell; `kill_fraction` is a tumour-only quantity |

**Thresholds are parameters, and the two that would otherwise be arbitrary are derived from the
run's own configuration rather than picked:**

* `adequate_ratio = 1.0` is not a choice — `E == IC50` is the definition of the half-effect
  point and exactly where the Hill term reaches half its ceiling.
* `unreached_ratio = 0.1` is the one conventional constant; at `hill = 1.2` it is 5.9% of the
  drug's achievable kill probability. It is exposed so a caller can move it and watch the
  breakdown move.
* `resistant_modulation` **derives** as `max(pd.hypoxia_resistance, pd.quiescent_resistance)` —
  the mildest resistance penalty the kill model itself applies (0.5 on the production `PDConfig`).
* `stochastic_min_p_kill` **derives** as `f(adequate_ratio) * resistant_modulation`, the lowest
  `p_kill` the CHANCE branch can produce at `max_kill_rate == 1`. Consequence: `UNATTRIBUTED` is
  **empty by construction** under the defaults, which is what makes it a live check rather than
  decoration — it fills up precisely when something outside the bands suppressed the kill.

**Precedence.** A cell can be short of drug *and* resistant. One label fits in a byte, so
exposure bands bind first: below `IC50` no restoration of sensitivity lifts the Hill term past
half its ceiling, so the exposure is the binding constraint. Those cells are not hidden —
`joint_limited.count` reports how many exposure-limited survivals are *also* below the resistance
threshold, i.e. how many would still not have died on dose alone.

**Reconciliation is mandatory and bit-exact.** `mechanism.reconcile` rebuilds `kill_fraction`
from the labels with the *same* integer arithmetic `CellPopulation.viable_tumor_fraction` uses —
`1.0 - (n_tumor - (KILLED + ALREADY_DEAD)) / n_tumor` — and compares with `==`, not `allclose`.
An attribution that cannot reproduce the number it claims to explain is describing a run that did
not happen, so `simulate_response` raises `MechanismReconciliationError` rather than returning a
plausible-looking breakdown, and no payload is written.

---

## Therapy reports

A **Report** is the deliverable: one ranked, mechanistically attributed drug panel for ONE
specimen, with the provenance of every input and the limits of the output attached. It is the
in-product form of `experiments/patient_report.py`, and it runs the *same* attribution —
`tmesim.report.attribution.attribute`, imported by both. The science is shared, never forked.

**A report is made of runs.** `POST /api/reports` enqueues **2 × len(drugs)** ordinary
`RunSpec`s on the ordinary job queue: for each drug the real arm, and the
`resistance_off: true` counterfactual whose kill is the denominator of
`rho = kill / kill(resistance off)`. Every arm shares one
`(section, dx, vascular_mode, oxygen_mode, state_mode)` key, so the expensive reconstruction is
built once and every later arm is a cache hit. The arms are visible in `/api/runs`, cancellable,
and persisted like any other run. Nothing bypasses the queue or the cache.

`sensitivity_mode` is **forced to `uniform`** and cannot be overridden: it is the only mode in
which a cross-drug ranking is fair, because an identical per-cell factor cancels out of the
comparison. `drug` and `resistance_off` are likewise owned by the report (they are swept). Any
other `RunSpec` field may be supplied in `spec_overrides`; overriding an owned field is a
warning, not an error, and the override is reported back as ignored.

```jsonc
// POST /api/reports body
{ "section": "xenium_319",
  "drugs": ["temozolomide","lomustine","cisplatin","methotrexate","doxorubicin"],  // default: this panel
  "spec_overrides": { "dx_um": 12.0, "total_time_s": 21600.0 },                    // any RunSpec field
  "label": "Patient A · five-drug panel" }                                         // optional
```

```jsonc
// Report — GET /api/reports/{id}
{ "report_id": "rep_00001", "status": "running", "section": "xenium_319",
  "drugs": [...], "n_runs": 10, "run_ids": ["r_1", ...],
  "arms": { "cisplatin": { "real": "r_5", "counterfactual": "r_6" }, ... },
  "progress": { "runs_total": 10, "runs_done": 4, "fraction": 0.4,
                "by_status": {"done": 4, "queued": 6}, "stage": "simulate (r_00007)" },
  "spec_base": { ...RunSpec with drug=null... }, "spec_overrides": {...},
  "sensitivity_mode": "uniform",
  "not_for_clinical_use": true, "not_for_clinical_use_text": "...",

  // THE LIMITS. First section of every rendering, never a footnote. `applies` is a property
  // of THIS section and THIS spec; an item that does not apply is still returned, saying why.
  "limits": [ { "id": "ranking_is_the_output", "severity": "critical", "applies": true,
                "title": "The RANKING is the defensible output. The absolute kill is not.",
                "body": "...", "source": "results/cell_state_impact.md §3.2, §3.3, §5" }, ... ],

  // THE RANKING, sorted by predicted kill. Each row carries its attribution AND its provenance.
  "ranking": [ { "rank": 1, "name": "cisplatin",
                 "kill": 0.3699, "kill_no_resistance": 0.7377, "rho": 0.5014,
                 "tumor_auc": 1.69e6, "c_tumor_uM": 78.4, "ic50": 24.6,
                 "L_um": 160.0, "d_eff": 225.0, "bbb": 0.858, "pen": 0.420,
                 "delivery_decades": 0.98, "potency_decades": -0.06,
                 "hypoxia_decades": -0.01, "total_decades": 0.91,
                 "dominant": "DELIVERY (advantage)", "gap_reason": "Versus **doxorubicin** ...",
                 "provenance": { "run_id": "r_5", "counterfactual_run_id": "r_6",
                                 "alpha": 148.6, "alpha_source": "...",
                                 "ic50": { "ic50_uM": 24.6, "confidence": "measured",
                                           "source": "GDSC measured (n=34 glioma cell lines)",
                                           "n_cell_lines": 34, "note": "...",
                                           "provenance_kind": "measured" },
                                 "sensitivity_mode": "uniform",
                                 "sensitivity_upper_bound": true,
                                 // "unmodelled - MGMT not on panel" where that is the truth
                                 "mgmt": { "state": "not applicable", "detail": "...", "value": null },
                                 "seed": 0, "state_mode": "oxygen", "caveats": ["..."] } } ],

  "attribution": { "terms": ["delivery_decades","potency_decades","hypoxia_decades"],
                   "term_labels": {...}, "method": "C_tumor_i = ...",
                   "implementation": "tmesim.report.attribution.attribute",
                   "check": { "spearman_total_vs_kill": 1.0, "faithful": "faithful", ... } },
  "microenvironment": { ...reconstruction_report verbatim... },
  "caveats": ["...the de-duplicated union of every arm's own caveats..."],
  "top_ranked": { "drug": "cisplatin", "kill": 0.3699, "dominant": "DELIVERY (advantage)" },
  "error": null }
```

`GET /api/reports/{id}/html` renders that document as **one self-contained HTML file**: the
Studio's own `app/web/css/tokens.css` is inlined, and there is no external stylesheet, webfont,
image, script or fetch of any kind. It must render identically off a USB stick. It carries the
**NOT FOR CLINICAL USE** banner at the top *and* the bottom, and the limits block as section 1.
It returns **409** with the current progress until every arm is `done` — a partial panel is not
a ranking, so no partial document is ever produced.

Reports persist as `app/reports/<report_id>/report.json` (+ `report.html` once done), gitignored
like `app/runs/`.

---

## Front-end responsibilities

- WebGL2 point rendering, **must hold 850k cells at ≥30 fps** with zoom/pan and per-cell picking.
- Channel compositing (cell types as toggleable layers), field overlays with a colour ramp
  and a real legend carrying units.
- A **code/JSON editor** for `RunSpec` validated against `POST /api/validate`, kept in sync
  with the form controls in both directions.
- Run queue with live status polling, cancel, and side-by-side comparison of two runs.
- Every metric shown with its provenance chip.
- A **project / specimen browser** and a **cross-specimen comparison view**. Rank agreement is
  presented first and largest; an absolute-kill delta whose `comparable` is `false` must be
  rendered with its blockers, never as a bare number.
- A **Reports panel**: trigger a report for the current section, watch its arms build, view the
  rendered document inline, and download it as one file. The stub in `js/mock.js` deliberately
  does **not** synthesise reports — a fabricated clinical-looking document is the one thing this
  product must never produce — so the panel says so and refuses rather than inventing one.
- A **mechanism colour mode** (`colour by: mechanism`) painting every cell from the TMSM label
  payload through the existing colour-by-attribute path in `gl/shaders.js` — a third vertex
  attribute, not a second renderer — with a legend naming each category, its count and its
  share, the delivery/dose/resistance/chance roll-up as one stacked bar, and the bit-exact
  reconciliation stated on screen. Legend rows toggle their category, and the cull happens in
  the vertex stage like every other channel. Every category carries a visible **label and a
  number**: the palette clears its CVD target on the dark surround but three slots fall under
  3:1 on the light one, so identity is never colour-alone.
- A **time-lapse scrubber** (play / pause / seek) over the TMSL frames, windowed across the
  whole series rather than per frame so scrubbing shows the drug moving instead of the ramp
  re-normalising. The frames are drug **concentration**, and the legend says so — it is a
  different quantity from the `drug_auc` overlay and must never share its label.
- `js/mock.js` synthesises **neither** payload. A fabricated mechanism map would be a picture of
  a reason that did not exist; a run without one says why instead.

## Backend responsibilities

- Wrap the EXISTING engine. Do **not** reimplement science. Import
  `tmesim.reconstruct.real_pipeline`, `tmesim.sim.engine`, `tmesim.validation.*`.
- Run jobs in a worker thread/process so the API stays responsive; stream stage + progress.
- Cache reconstructions by `(section, dx, vascular_mode, oxygen_mode, state_mode)` — they are
  expensive and frequently reused across drug runs. `state_mode` is in the key because it
  rewrites `cells.state`; serving one state route from another's entry would report the wrong
  band occupancies, which is the largest single term in an absolute kill.
- Persist projects under `app/projects/` (gitignored). A specimen stores `run_id`s; runs live
  only in `app/runs/` and are never copied into a project.
- Never invent a number. If a quantity is unavailable, return `null` and a reason.
