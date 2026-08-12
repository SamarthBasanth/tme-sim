# tmesim Studio

A local diagnostic-imaging application over the tmesim engine: reconstruct one
real glioblastoma section, simulate a drug on it, and inspect the result at
single-cell resolution. `app/server/` is the backend; `app/web/` is the front
end. Both sides implement the frozen [`API_CONTRACT.md`](../API_CONTRACT.md) v1.

```
./.venv/bin/python -m pip install -e '.[app]'      # fastapi + uvicorn + httpx
./.venv/bin/python -m app.server                   # http://127.0.0.1:8000
```

Binds to `127.0.0.1` only. This is a workstation tool that serves patient tissue
reconstructions; pass `--host` explicitly if you really mean to expose it.

**Not for clinical use.** Every number the API returns is a prediction of a
physical model, on tissue that carries no treatment-outcome label. Nothing here
has been scored against a measured drug response, because no such measurement
exists for these sections. See the `provenance.caveats` list on every result.

---

## What the backend is, and is not

It is a **wrapper**. `tmesim/core`, `tmesim/physics`, `tmesim/sim`,
`tmesim/molecular` and `tmesim/reconstruct` are untouched. Every quantity the API
reports comes out of a production entry point:

| Quantity | Produced by |
|---|---|
| geometry, cell types, tissue fields | `reconstruct.real_pipeline.reconstruct_real` |
| the reconstruction report | `reconstruct.real_pipeline.reconstruction_report` |
| per-cell sensitivity (+ MGMT scalar) | `reconstruct.real_pipeline.apply_drug_sensitivity` |
| perivascular-sink constant α | `validation.calibration.calibrate_alpha` |
| compiled `DrugParams` | `validation.calibration.calibrated_drugparams` → `molecular.binding_site_barrier.apply_barrier_to_drugparams` |
| measured IC50 + its source | `molecular.potency.resolve_ic50` |
| kill fraction, AUC field | `sim.engine.simulate_response` |
| delivery metrics | `validation.delivery_metrics.compute_all` |
| analytic penetration length | `validation.pk_reality.penetration_length_um` |

`tests/test_app_server.py::test_api_run_is_numerically_identical_to_the_cli`
drives that chain by hand, the way `experiments/patient_report.py` does, and
asserts **exact** equality (`==`, not a tolerance) against the same run submitted
through the API. If the API ever forks into a second, subtly different pipeline,
that test fails.

### Files

| File | Responsibility |
|---|---|
| `payloads.py` | TMSG / TMSF binary wire formats |
| `cache.py` | section registry + the reconstruction LRU |
| `jobs.py` | run execution, the queue, the worker, on-disk history |
| `projects.py` | projects / specimens (multi-patient) + the cross-specimen comparison |
| `report.py` | therapy reports: panel assembly, the limits block, the self-contained HTML |
| `api.py` | routes, drug/defaults registries, `POST /api/validate` |
| `__main__.py` | `python -m app.server` |

---

## Routes

All of `API_CONTRACT.md`:

| Method | Path | Notes |
|---|---|---|
| `GET` | `/api/sections` | `Section[]`; unavailable sections carry `available:false` + `reason` |
| `GET` | `/api/drugs` | real IC50 + `ic50_source`, MW, logP, compiled `bbb_penetrance` |
| `GET` | `/api/defaults` | a valid `RunSpec` at the production defaults |
| `POST` | `/api/runs` | `202` + `{run_id, status, warnings}`; rejects an invalid spec with `400` |
| `GET` | `/api/runs` | newest first; `?limit=` (default 200) |
| `GET` | `/api/runs/{id}` | full `RunResult` |
| `GET` | `/api/runs/{id}/status` | `{status, progress, stage, elapsed_s, log_tail}`; `?log_lines=` |
| `DELETE` | `/api/runs/{id}` | cancel |
| `GET` | `/api/runs/{id}/geometry` | TMSG binary |
| `GET` | `/api/runs/{id}/field/{name}` | TMSF binary; `name ∈ vascular_density, ecm_density, oxygen, drug_auc` |
| `GET` | `/api/runs/{id}/mechanism` | TMSM binary; per-cell mechanism labels + the aggregate breakdown |
| `GET` | `/api/runs/{id}/timelapse` | TMSL binary; drug-concentration frames + their metadata |
| `POST` | `/api/validate` | `{ok, errors[], warnings[]}` |
| `GET` | `/api/sections/{id}/capability` | what this panel can support, before a project exists |
| `GET`/`POST` | `/api/projects` | list / create a study |
| `GET`/`PATCH`/`DELETE` | `/api/projects/{id}` | one study; delete detaches runs, never deletes them |
| `POST` | `/api/projects/{id}/specimens` | bind a section as a patient (`409` if already bound) |
| `GET`/`PATCH`/`DELETE` | `/api/projects/{id}/specimens/{sid}` | one specimen + its resolved runs |
| `POST`/`DELETE` | `/api/projects/{id}/specimens/{sid}/runs[/{run_id}]` | associate / detach a run (`409` on a section mismatch) |
| `GET` | `/api/projects/{id}/compare` | cross-specimen comparison; see below |
| `POST` | `/api/reports` | `202` + `{report_id, status, run_ids, warnings}`; enqueues 2 runs per drug |
| `GET` | `/api/reports` | newest first; `?limit=` (default 200) |
| `GET` | `/api/reports/{id}` | the full report document; also the poll target |
| `GET` | `/api/reports/{id}/html` | the rendered, self-contained document; `409` until every arm is done |
| `DELETE` | `/api/reports/{id}` | cancel a queued/running report (cancels every arm) |

Two additive diagnostics the UI need not use: `GET /api/health` (cache and queue
counters) and `GET /api/runs/{id}/traceback` (the stack for a `failed` run).
`app/web/` is mounted at `/` **last**, so it can never shadow `/api`.

Responses are strict JSON: non-finite floats are serialised as `null`, never as
the invalid `NaN` token, and a `null` always travels with a reason
(`metric_reasons`, `*_reason`, or `provenance.caveats`).

---

## Sections

Enumerated from the bytes actually present in `data_cache/`. A missing cache
gives `available:false` with a specific reason instead of an exception.

| id | source | cache | cells |
|---|---|---|---|
| `xenium_319` | Xenium | `data_cache/xenium_gbm` | 40,887 |
| `xenium_io` | Xenium | `data_cache/xenium_gbm_io` | 816,769 |
| `gse324268` | Xenium | `data_cache/gse324268_xenium` | 851,509 |
| `visium_10x_parent_gbm` | Visium | `data_cache/visium_gbm` | placed by deconvolution |
| `visium_gbm1` / `gbm2` / `gbm3` | Visium | `data_cache/visium_gbm/GBM{1,2,3}` | placed by deconvolution |

`n_cells` and `extent_um` are read from `cells.csv.gz` centroids (Xenium) or
`tissue_positions` (Visium), so they are the raw measured extent and do **not**
move with `dx` the way `Domain.extent_um` does.

For a Visium section `n_cells` is `null` with a reason on a cold checkout: spots
are not cells, and the cell population is *placed* by deconvolution at
reconstruction time. Once a section has been reconstructed once, the observed
count is remembered in `app/runs/_section_meta.json` and reported.

---

## Jobs

`POST /api/runs` validates, enqueues, and returns. A single background **worker
thread** drains the queue.

* **Thread, not process.** The expensive object is the cached reconstruction — a
  graph of large numpy arrays. A process pool would have to rebuild or pickle it
  per run, which is the entire cost the cache exists to avoid. The heavy inner
  loops (numpy, `scipy.sparse` factorise/solve) release the GIL, so the API stays
  responsive during a run.
* **One worker.** Runs on the same key share the cached `Microenvironment`, and
  `apply_drug_sensitivity` writes `cells.sensitivity` on it. Serialising the
  worker makes that safe by construction, and stops N runs fighting over the same
  cores. Every run additionally restores the pristine sensitivity **and cell-state**
  array snapshots taken at build time, so one run can never inherit another's
  assumption — asserted by
  `test_a_run_never_inherits_the_previous_run_sensitivity`.

**Stages** are `reconstruct → sensitivity → compile_drug → simulate → metrics →
payloads → done`, with `progress` fixed per stage. Progress is deliberately
**stage-granular**: the engine entry points are single blocking calls with no
progress callback, and this backend does not modify them, so a within-stage
percentage would be a fabricated number. `log_tail` carries the real call log —
which entry point ran with which arguments, and what it returned.

**Cancel.** A queued job is cancelled immediately. A running job is cancelled
*cooperatively*: the flag is checked at every stage boundary, so a job already
inside `simulate_response` finishes that solve first. The engine is not patched
to accept an interrupt. `DELETE` says exactly this in its `reason`.

**Persistence.** On completion a run writes `app/runs/<run_id>/`:

```
result.json            the full RunResult + the whole log
geometry.bin           TMSG
field_vascular_density.bin
field_ecm_density.bin
field_oxygen.bin
field_drug_auc.bin     (only when the spec had a drug)
```

The manager reloads these at startup, so history — including the payloads the UI
renders — survives a restart, and run ids keep counting rather than colliding.

---

## Reconstruction cache

Keyed by exactly `(section, dx_um, vascular_mode, oxygen_mode, state_mode)`, as the
contract specifies. `dx_um` is rounded to 6 dp so `12.0` and `12.000000001` share an
entry.

`state_mode` is in the key because it **rewrites `cells.state`** — `'vessel_distance'`
and `'markers'` both overwrite it, `'oxygen'` leaves it exactly as the reconstructor
produced it. Two specs that differ only in `state_mode` are therefore two different
reconstructions, and serving one from the other's entry would report the wrong band
occupancies. That is not a staleness nicety: the band occupancy is the largest single
term in every absolute kill number here, and it moves the kill by up to 0.3252 absolute
(1.77x relative) across the three routes with tissue, drug, transport and RNG held
identical (`results/cell_state_impact.md` §3.2). A cached reconstruction under a
different `state_mode` is **wrong to reuse**, not merely out of date.

Measured build times on this machine:

| section | cells | dx | build |
|---|---|---|---|
| `xenium_319` | 40,887 | 30 µm | ~1–2 s |
| `xenium_io` | 816,769 | 30 µm | ~20 s |
| `gse324268` | 851,509 | 30 µm | ~32 s |
| `visium_*` | placed | 100 µm | ~9–14 s |

A drug run at `dx = 30 µm` on `xenium_319` is ~6 s cold and ~4 s warm; at the
production `dx = 12 µm` the reconstruction dominates, which is why this cache is
the single biggest responsiveness win in the backend. It is in-process (a disk
round-trip of these arrays would cost most of the saving) and LRU-bounded at 4
entries by default — each holds a full cell population and several float64 grids.
Override with `--cache-size` or `TMESIM_RECON_CACHE_SIZE`.

Two jobs wanting the same key do not both reconstruct: the second blocks on a
per-key lock and then finds the first one's entry.

### The one place the wrapper reaches past `reconstruct_real`

`reconstruct_real` builds the **discrete** segmented-vessel field
(`real_fields._DEFAULT_VASCULAR_MODE`) and exposes no `vascular_mode` parameter,
so `vascular_mode='kde'` cannot be expressed through it. Rather than reimplement
the assembler, `build_reconstruction` runs the production reconstruction and then
swaps in `real_fields.vascular_density_from_cells(..., mode='kde')` — the same
public function on the same measured endothelial cells — and re-installs oxygen
via `real_pipeline._install_oxygen`, because the `pde` and `vessel_distance`
routes are driven *by* the vascular field. The cost (a discarded discrete field
and one discarded oxygen solve) is logged, and a note lands in
`provenance.caveats`. `vascular_mode='discrete'` — the default — is the
unmodified production call, which is what CLI parity is asserted on.

`vascular_mode` has no meaning for a Visium section (`reconstruct_from_visium`
builds vascular density from the measured endothelial *spot fraction* and has no
mode switch); a `'kde'` request there is reported as **not applied** rather than
silently ignored.

---

## Binary payloads

Little-endian, no padding, laid out so the front end slices one `ArrayBuffer`
into typed arrays and uploads straight to WebGL2.

```
geometry:  "TMSG" | version u32 | n_cells u32 | extent_x f32 | extent_y f32
                  | x f32[n] | y f32[n] | type u8[n]
field:     "TMSF" | version u32 | rows u32 | cols u32 | dx_um f32
                  | min f32 | max f32 | values f32[rows*cols]
mechanism: "TMSM" | version u32 | n_cells u32 | n_labels u32 | json_len u32
                  | aggregate json[json_len] | label u8[n_cells]
timelapse: "TMSL" | version u32 | n_frames u32 | interval_s f32 | window_s f32
                  | n_frames x [ frame u32 | time_s f32 | <a whole TMSF block> ]
```

Header sizes are 20 B, 28 B, 20 B and 20 B (+8 B per time-lapse frame).

**TMSM** is one byte per cell in the geometry payload's order, so the front end
binds it as a third vertex attribute against the position buffer it already
holds — 850 kB against a 3.4 MB position buffer, and the measured presented
frame rate at 850k cells is unchanged by it. The aggregate breakdown rides in
the same response as a short JSON block; a legend without its counts is
decoration, and a second round trip for four integers is not worth a request.

**TMSL** does not invent a grid format: each frame is a **byte-identical
`encode_field` block** prefixed with its index and time stamp, so the client
decodes a frame with the TMSF reader it already ships. `test_mechanism.py::
test_tmsl_frames_are_byte_identical_to_the_tmsf_field_payload` asserts that
literally, by searching for the standalone TMSF bytes inside the container.
Frames are drug **concentration** (uM), not AUC — AUC only accumulates, so it
can never show drug leaving. Field values are **row-major**
(`values[r*cols + c]`), matching the engine's `(ny, nx)` grid. `min` / `max` are
the **finite** range, so the UI can build a colour ramp without a second pass; a
field with no finite entry reports `nan` for both rather than a fabricated `0..1`.

Why binary: `gse324268` has 851,509 cells. Its geometry payload is
`20 + 9 × 851,509 = 7.66 MB`, versus roughly 40 MB of JSON text that then has to
be parsed into JS numbers.

The layout is verified by round-trip, not trusted:
`TestPayloadRoundTrip::test_geometry_header_is_little_endian_and_exact` unpacks
every header field at its literal byte offset, checks the magic bytes, checks
that a big-endian read of the cell count gives a *different* answer, and pins the
total length; `test_field_round_trip_and_row_major` re-reads `values[5][7]` from
the flat body at `5*cols + 7`.

---

## Provenance

Every `RunResult.provenance` carries, mandatorily:

* `alpha` and `alpha_source` — the single perivascular-sink constant and what
  measurement it was fit to;
* `ic50_uM`, `ic50_source` (`"GDSC measured (n=51)"` / `"curated literature"` /
  `"flagged default (LOW confidence)"`), plus `ic50_confidence`, `ic50_note`
  and `ic50_n_cell_lines`;
* `sensitivity_mode`, `sensitivity_upper_bound`, `sensitivity_note`,
  `sensitivity_determinant_gene`, and the `mgmt_*` block when an MGMT value was
  supplied — including `mgmt_applied:false` **with a reason** when the drug is
  not MGMT-calibrated;
* `caveats[]` — the oxygen route's own falsification note, the upper-bound
  warning, the grid-convergence status at this `dx`, the Stage-0 PK
  simplification, and the absence of any outcome label;
* the solver settings the frozen `RunSpec` cannot carry — `pde_dt_s`,
  `pde_save_every`, `pde_convection`, `pde_time_resolved`, `seed` — read off the
  engine's own dataclass defaults rather than retyped, plus the list of engine
  entry points that produced the result.

`RunSpec` has no solver fields, so the backend pins them at the values every
tmesim experiment uses (`dt = 10 s`, `save_every = 500`, `convection = True`,
`seed = 0`) and leaves `PDEConfig.time_resolved` at the engine default. Those are
the same constants `experiments/patient_report.py` and
`validation/delivery_metrics.py` run, which is the other half of why CLI parity
holds.

Unavailable quantities are `null` **with a reason**, never a substituted default.
Concretely: a `drug:null` spec returns all four contract metrics as `null` with
`metric_reasons` explaining that they need a drug, and `/field/drug_auc` returns
`404` saying the same; a `penetration_length_um` whose exp+floor fit did not
converge is `null` with that stated, not replaced by the analytic length (which
is reported separately in `extra_metrics`).

One honesty note the API surfaces itself: `sim.engine.simulate_response`
re-solves the oxygen PDE internally and overwrites `fields.oxygen`. So with
`oxygen_mode` set to `vessel_distance` or `markers`, the field served by
`/field/oxygen` is the **reconstruction's** oxygen and did *not* drive the kill.
That is stated in `provenance.caveats` and as a `/api/validate` warning.

The same is true of `state_mode`, and for the same reason: `simulate_response`
re-bands every living cell through `cells.agents.update_cell_states` from its own
oxygen solve before applying the kill. So `state_mode` moves the **reported** band
occupancies and their provenance, not this run's `kill_fraction`. It is returned as
`provenance.cell_state` (the pipeline's own block, verbatim: mode, source,
`is_measured`, caveat, marker-coverage detail, tumour band fractions), plus
`provenance.state_unmeasured` and `provenance.state_drives_kill = false`, and it is
stated as a `/api/validate` warning. Propagating a state route into a kill means
re-running `cells.pkpd.apply_drug_kill` on the engine's own transport result, which is
what `validation.cell_state_impact` does — the API deliberately does not fake it.

---

## `POST /api/validate`

Messages are written to be shown verbatim in the UI's code editor. Errors mean
the run would fail or be meaningless; warnings mean it would run but needs care.

Caught as **errors**: unknown or unavailable `section`; unknown `drug`; any
invalid mode; non-numeric or non-positive `dx_um`; `total_time_s` below the
solver step `dt` (which would take zero steps and report every exposure as 0);
`oxygen_mode='markers'` on Visium, or on a panel carrying none of
`HYPOXIA_MARKER_GENES` (the 319-gene brain panel carries none, so hypoxia there
is *unmeasured*);
`state_mode='markers'` on Visium (no per-cell transcripts) or on a panel carrying
neither a proliferation nor a hypoxia-response transcript (nothing would be measured
at all); `state_mode='vessel_distance'` on a panel with no endothelial marker in any
registered panel map — no cell can then be typed ENDOTHELIAL, the distance transform
has no point to measure from, and **every** cell would be labelled HYPOXIC by
construction, so `reconstruct_real` refuses rather than returning a 100%-hypoxic
section that looks like a result; `pk.plasma_peak <= 0`, `pk.half_life_s <= 0`, `pk.n_doses < 1`
or fractional, `pk.dose_interval_s <= 0` with multiple doses; a non-numeric,
negative, or non-`methylated`/`unmethylated` `mgmt`.

Caught as **warnings**: `dx` outside the convergence regimes from
`results/numerical_rigor.md` §2 (20 µm for `kill_fraction`, ~15 µm for
`penetration_index`); a `dx` that over-resolves the 100 µm Visium spot pitch and
will therefore be overridden by the loader; grid sizes that make a run take
minutes; `sensitivity_mode='determinant'` (not comparable across drugs, and inert
on Visium); the surrogate oxygen routes; unrecognised spec keys; `mgmt` on a
non-MGMT-calibrated drug; clinically unattainable PK; a dose interval longer than
the window; every non-default `state_mode` (it moves reported band occupancies, not
the kill, and forces a cache rebuild); `state_mode='markers'` on a panel missing one
axis — naming the axis, the genes, and the fact that the state is then assigned to
**no cell**, which is not the same statement as "computed and found to be zero"; and,
on any panel carrying no hypoxia-response transcript, that **absolute kill must not be
quoted from that section under any state route**, because no route clears the bar
there and the kill moves up to 0.3252 absolute across the three
(`results/cell_state_impact.md` §5). The cross-drug **kill** ranking is invariant
across all three **state** routes and is the defensible output of that
comparison. It is **not** invariant across the vascular delivery law: correcting
the one-way source `S = PS*Cp` to `S = PS*(Cp - C)` re-weights every drug by its
own `bbb_penetrance` (47x spread on the GBM panel) and changes the GBM top-1,
cisplatin -> vincristine
([`results/transvascular_exchange_revalidation.md`](../results/transvascular_exchange_revalidation.md)
§4.1). The warning strings say so.

### The advection CFL

`physics.diffusion.simulate_drug_transport` integrates advection with an explicit
upwind step and raises when `dt·max|v|/dx > 1` — tens of seconds into a run,
after the reconstruction has already been paid for. `max|v|` is the peak
interstitial Darcy speed, a property of a *built* reconstruction, so it is
measured with the production `physics.pressure.solve_ifp` whenever one is built
and recorded in `app/runs/_cfl_speeds.json` (so the evidence survives a restart).

Validation then uses that evidence in three labelled tiers:

1. **this exact key** has been built → the speed is this spec's own measurement;
   a violation is an **error** quoting the CFL number and the minimum safe `dx`;
2. **same section, same `vascular_mode`, other mesh** → `solve_ifp` reads the
   vascular field and the domain and never reads oxygen, so those keys share this
   spec's vascular source; a violation is still an **error**, labelled with the
   key it was measured on;
3. **only a different `vascular_mode`, or nothing** → a **warning** that the CFL
   is unverified, and why.

No speed is ever assumed. On `xenium_319` the measured peak is ≈0.60 µm/s, so at
`dt = 10 s` the limit is `dx ≳ 6 µm` — comfortably below the production 12 µm.

---

## Projects and specimens — the multi-patient layer

`projects.py`. A **project** is a study; a **specimen** is a patient, bound to
one `section` id and to the *capability report* of that section's panel. A
specimen holds `run_id`s **into `app/runs/`** — runs are associated, never
copied. Deleting a project or detaching a run therefore never deletes a run, and
attaching a run whose `spec.section` differs from the specimen's section is a
`409`, because filing another patient's run under a specimen would silently
corrupt every comparison built on the study. Projects persist as one JSON
document per project under `app/projects/<id>/project.json` (gitignored).

`GET /api/projects/{id}/compare` is the point of the whole layer, and its shape
is dictated by what the evidence in `results/` supports:

* **Rank agreement is the primary result.** `results/heldout_patient.md` §3.3
  found the kill ordering IDENTICAL on two real GBM sections;
  `results/cell_state_impact.md` §3.3 found it INVARIANT across all three
  cell-state routes. The comparison reports the per-specimen ordering, a
  Kendall tau-b per pair, and the swapped drug pairs when it disagrees.
* **Absolute kill is flagged, not differenced.** `comparable` is `false` unless
  *every* specimen's panel carries a hypoxia-response transcript **and** every
  matched run actually used a MEASURED state route. That is
  `results/cell_state_impact.md` §5: with no hypoxia gene no route clears the
  bar, and kill moves by up to **0.3252 absolute (1.77×)** with the choice of
  route alone. The 319-gene flagship panel fails this; the 480-gene IO panel
  passes it (VEGFA). A measured kill is never differenced against an assumed one.
* **Only identical RunSpecs are matched** — `spec_signature` folds in the mesh,
  every mode, `mgmt`, `resistance_off`, the whole `pk` block and the window, and
  fails closed on any RunSpec field it has never heard of.
* **Composition rows carry their own verdict.** `mean_ecm_density` and
  `mean_oxygen` are *always* `comparable: false` (a max-normalised KDE over
  whichever matrix transcripts the panel carries has no common absolute scale;
  the oxygen level is uncalibrated on every section — §3.5). The ENDOTHELIAL
  *fraction* and *areal density* rows sit adjacent because on the real pair they
  moved in opposite directions, 0.43× and 2.67× (§1.1).

Front end: `app/web/js/ui/projects.js` + `css/projects.css`, a modal reached
from the `projects` button in the top bar.

---

## Therapy reports

A **report** is the deliverable: one ranked, mechanistically attributed drug
panel for ONE specimen, with the provenance of every input and the limits of the
output attached. It is the in-product form of `experiments/patient_report.py`.

**`report.py` simulates nothing.** A report is assembled out of runs:

* every arm is an ordinary `RunSpec` on the ordinary queue, so it shares the one
  worker, the reconstruction cache, cancellation and `app/runs/` persistence;
* each drug contributes **two** arms — the real run, and the
  `resistance_off: true` counterfactual (`PDConfig(hypoxia_resistance=1,
  quiescent_resistance=1)`) whose kill is the denominator of
  `rho = kill / kill(resistance off)`;
* all 2N arms share one cache key, so the reconstruction is built once and every
  later arm is a hit — that is the whole reason a report goes *through* the queue
  rather than around it;
* when every arm is `done` the panel is assembled by
  `tmesim.report.attribution.attribute` — the **same function object**
  `experiments/patient_report.py` calls, asserted in the tests. Two
  implementations of an attribution would be two sciences.

`sensitivity_mode` is forced to `uniform` and cannot be overridden: it is the
only mode in which a cross-drug ranking is fair. `drug` and `resistance_off` are
swept per arm. Any other `RunSpec` field may be passed in `spec_overrides`;
overriding an owned field is a warning and is reported back as ignored.

**The document carries its own limits, prominently.** `build_limits()` returns
the block that opens every rendering — the NOT FOR CLINICAL USE banner, the
ranking-is-the-output rule (`results/cell_state_impact.md` §3.2/§5, kill moving
0.3252 absolute across state routes), whether *this* panel can support an
absolute kill at all (the 319-gene panel carries no hypoxia gene and calls
proliferation in 1.42% of tumour cells against a published 19-40%), the interval
absolute kill is defensible as (TMZ 0.0048-0.1442, 30x —
`results/kill_calibration.md` §5.1), the uniform-mode upper bound, the oxygen
PDE's falsification (`results/oxygen_demotion.md`), and the absence of any
outcome label anywhere in this repository. An item that does not apply to a
section is still returned, saying why — a caveat that disappears cannot be
audited.

> **Outstanding on the ranking-is-the-output limit.** Its body in
> `app/server/report.py::build_limits` (`id: ranking_is_the_output`) still reads
> "INVARIANT across all three cell-state routes", unscoped. The claim is true of
> the state routes and is not true of the vascular delivery law, which changes
> the GBM top-1 — see
> [`results/transvascular_exchange_revalidation.md`](../results/transvascular_exchange_revalidation.md)
> §4.1. That string was outside the file set of the pass that scoped the rest of
> these surfaces and still needs the same treatment.

`render_html()` inlines `app/web/css/tokens.css` and emits **one self-contained
file**: no external stylesheet, webfont, image, script or fetch of any kind, so
it reads the same off a USB stick as it does in the product. The banner is at
the top *and* the bottom; the limits are section 1, above the ranking.

Reports persist to `app/reports/<id>/report.json` (+ `report.html` once done),
gitignored like `app/runs/`.

Front end: `app/web/js/ui/reports.js` + `css/reports.css`, a third dock column
with a modal document viewer. `js/mock.js` deliberately does **not** synthesise
reports — a fabricated clinical-looking document is the one thing this product
must never produce — so the panel says so and refuses.

---

## Tests

```
./.venv/bin/python -m pytest tests/test_app_server.py -q      # 30 passed, ~90 s
./.venv/bin/python -m pytest tests/test_projects.py -q        # 41 passed, ~2 s
./.venv/bin/python -m pytest tests/test_app_reports.py -q     # 34 passed, ~60 s
```

`tests/test_app_server.py` covers, in rising order of what a break would cost:

1. **payload round-trip** — literal byte offsets, magic, endianness, row-major
   ordering, the finite-range rule, and malformed-input rejection;
2. **job lifecycle** — submit → queued → running → done with live stage/progress/
   log polling, queued-run cancel, failure surfacing with a traceback, and
   history reloaded by a fresh manager over the same directory;
3. **cache** — one reconstruction reused across drugs (miss then hit, asserted on
   the counters *and* on `report.reconstruction_cache_hit`), a different `dx`
   being a different key, key normalisation, and sensitivity isolation between
   runs;
4. **CLI parity** — exact equality against the hand-driven production chain.

The file skips itself cleanly if `fastapi`/`httpx` are absent, and each
end-to-end test skips if `data_cache/xenium_gbm` is not present, so it never
breaks a checkout without the extra or without the data.

`tests/test_app_reports.py` covers three separable things, because they fail
separately: (1) the science is *imported, not re-implemented* — the Studio and
`experiments/patient_report.py` must hold the identical function object, and the
decomposition's own `delivery + potency == log10 drive - mean` identity is
checked on numbers it did not choose; (2) every mandatory limit is present with
its measured number and survives rendering, the banner is at top and bottom, and
the limits are positioned before the first kill figure in the document; (3) the
rendered file references nothing external at all. The end-to-end arm asserts the
two arms per drug are real, spec-identical apart from `resistance_off`, that the
reconstruction was built once, that the counterfactual never kills *less* than
the real arm (if it did, `rho` would not be a resistance ratio), and that every
row's numbers came out of its own arm rather than being recomputed.

`tests/test_projects.py` covers the multi-patient layer on the same principle:

1. **persistence** — a project survives a fresh `ProjectStore` over the same
   directory, with its specimens, capability blocks and run associations intact,
   and ids do not collide after a restart;
2. **association** — the project document holds ids and no result body, attach is
   idempotent, a run from another section is a `409`, and detach / remove-specimen
   / delete-project all report the runs they retained;
3. **the comparison's capability-aware flagging** — rank agreement and the swapped
   pairs, the 319-vs-IO panel asymmetry, measured-against-assumed refusal, the
   composition verdicts, and an UNMEASURED state band coming back `null` rather
   than a measured-looking `0.0`.

The flagging tests build their `RunResult` views by hand, so the whole decision
table runs in milliseconds; the end-to-end path is exercised by the two real
sections wherever they are cached.
