# tmesim Studio — API contract (v1)

Frozen interface between the backend (`app/server/`) and the front end (`app/web/`).
Both sides code to this. Do not change a route or payload shape without updating this file.

Stack: **FastAPI + uvicorn** backend, **vanilla ES modules + WebGL2** front end (no build
step, no npm — it must run from a clean checkout with only the existing venv plus fastapi).

Run: `./.venv/bin/python -m app.server` → serves the UI at `http://127.0.0.1:8000`.

---

## Design constraints (both sides)

This is diagnostic imaging software. Reference points are HALO, Visiopharm, Paige, Sectra —
**not** consumer web. Concretely:

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
  "mgmt": null,                    // optional float, log2 CPM
  "pk": { "plasma_peak": 1.0, "half_life_s": 14400, "n_doses": 1, "dose_interval_s": 86400 },
  "total_time_s": 600.0 }
```

```jsonc
// RunResult
{ "run_id": "r_00007", "spec": {...}, "status": "done",
  "report": { ... reconstruction_report ... },
  "metrics": { "kill_fraction": 0.115, "mean_tumor_auc": 1.2e7,
               "therapeutic_coverage_1x": 0.31, "penetration_length_um": 206.0 },
  "provenance": { "alpha": 148.6, "ic50_uM": 606.0, "ic50_source": "GDSC measured (n=51)",
                  "sensitivity_note": "...", "caveats": ["..."] },
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
| `POST` | `/api/validate` | `RunSpec` → `{ok, errors[], warnings[]}` — for the code editor |

`status` ∈ `queued | running | done | failed | cancelled`. On `failed`, `RunResult.error`
carries the message; the UI must surface it, never swallow it.

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
Both little-endian. The front end reads these into typed arrays and uploads straight to WebGL.
JSON for 850k cells is not acceptable.

---

## Front-end responsibilities

- WebGL2 point rendering, **must hold 850k cells at ≥30 fps** with zoom/pan and per-cell picking.
- Channel compositing (cell types as toggleable layers), field overlays with a colour ramp
  and a real legend carrying units.
- A **code/JSON editor** for `RunSpec` validated against `POST /api/validate`, kept in sync
  with the form controls in both directions.
- Run queue with live status polling, cancel, and side-by-side comparison of two runs.
- Every metric shown with its provenance chip.

## Backend responsibilities

- Wrap the EXISTING engine. Do **not** reimplement science. Import
  `tmesim.reconstruct.real_pipeline`, `tmesim.sim.engine`, `tmesim.validation.*`.
- Run jobs in a worker thread/process so the API stays responsive; stream stage + progress.
- Cache reconstructions by `(section, dx, vascular_mode, oxygen_mode)` — they are expensive
  and frequently reused across drug runs.
- Never invent a number. If a quantity is unavailable, return `null` and a reason.
