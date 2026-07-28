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

`status` ∈ `queued | running | done | failed | cancelled`. On `failed`, `RunResult.error`
carries the message; the UI must surface it, never swallow it.

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
{ "detection": { "path": "/data/Patient_042", "platform": "xenium",
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

> `path` is a **server-side filesystem path**. The server binds to `127.0.0.1` and is a
> local workstation tool; do not expose this route on a shared host.

```jsonc
// CapabilityReport — what this specimen can and cannot support, BEFORE a run
{ "sample_id": "patient_042", "platform": "xenium", "n_genes": 319,
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
