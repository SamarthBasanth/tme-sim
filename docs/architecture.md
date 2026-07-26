# tmesim — Architecture

This document describes how data flows through `tmesim`, what each module is
responsible for, the core data model, the physics, and — importantly — the honest
limitations. It is a companion to the frozen public API in
[`../ARCHITECTURE_CONTRACT.md`](../ARCHITECTURE_CONTRACT.md); where the two ever
disagree, the contract wins.

Conventions throughout: length in **micrometers (µm)**, time in **seconds**,
concentration in **µM**, diffusivity in **µm²/s**. All grid fields are 2D numpy
arrays of shape `(ny, nx)`, indexed `[row = y, col = x]` (image/numpy convention).

---

## 1. Pipeline data flow

```
 data/  ──►  reconstruct/  ──►  Microenvironment  ──►  physics/ + cells/  ──►  sim/  ──►  validation/
 (sample)     (build TME)        (core object)          (transport+kill)      (engine)   (beat baseline?)
                                                                                             │
                                                                                         viz/ + report
```

The `Microenvironment` (in `tmesim.core`) is the hinge of the whole system. It is
the **output** of reconstruction and the **input** to simulation. Every module
either produces one, consumes one, or transforms the fields inside one. Because
the object and all result containers are frozen in `core`, the physics, cell,
simulation, and validation modules can be developed and tested against a fixed
contract.

A single patient makes this trip once per `simulate_response` call; a cohort makes
it `n` times inside `run_stage0_experiment`, after which the predictions are scored
against baselines to produce the verdict.

---

## 2. Module responsibilities and public API

Signatures below are summarized from the contract. Every function takes/returns
`core` types; none invent their own return types.

### `tmesim.data` — get a cohort of patients

- **`synthetic.generate_synthetic_patient(patient_id, seed=0, grid=(120,120), dx=15.0, subtype=None) -> PatientRecord`**
  Builds a biologically-plausible synthetic GBM microenvironment with a *known*
  ground-truth response: a tumor mass + infiltrative margin, stroma, immune cells,
  a vascular network (sparser in the core → hypoxia), spatially varying ECM, and
  per-cell `sensitivity` correlated with metadata (MGMT methylation, subtype). The
  `observed_response` is a **noisy function of the true biology** (transport ×
  intrinsic sensitivity) — deliberately not trivially recoverable, so the test is
  fair. `features` carries non-spatial scalars for the baselines.
- **`synthetic.generate_synthetic_cohort(n=40, seed=0, drug_names=None) -> list[PatientRecord]`**
  A cohort spanning subtypes / vascularization / MGMT status for the Stage 0
  experiment.
- **`acquire.acquire_all(cache_dir="data_cache", timeout=60) -> dict`**
  Best-effort *live* download of open GBM data (Ivy GAP, TCGA-GBM via the GDC API,
  DepMap/GDSC glioma drug sensitivity) into `./data_cache/`. Degrades gracefully:
  returns `{source: {status: ok|failed|skipped, path, note}}` and never crashes.
  This is the only module permitted network access.
- **`loaders.load_real_cohort(cache_dir="data_cache") -> list[PatientRecord]`**
  Parses whatever was cached into `PatientRecord`s. Returns `[]` if nothing usable
  is present (the caller then falls back to synthetic).

### `tmesim.reconstruct` — sample → digital microenvironment

- **`segmentation.segment_and_phenotype(image, dx=15.0, markers=None) -> CellPopulation`**
  Segments nuclei from a (multi-channel) tissue image with scikit-image and
  assigns `CellType` by marker intensity. Robust to plain grayscale.
- **`spatial.build_tissue_fields(domain, cells, vascular_density=None) -> TissueFields`**
  Kernel-density-estimates cells into an ECM-density field (stroma-weighted) and,
  if not supplied, a vascular-density field. Leaves `oxygen`/`ifp` as `None` (the
  physics modules fill them).
- **`vasculature.build_vascular_network(domain, seed=0, tumor_mask=None) -> np.ndarray`**
  A `(ny, nx)` vascular-density field in `[0,1]`: tortuous and realistically
  **sparser inside the tumor core**, which is what drives both hypoxia and poor
  delivery downstream.
- **`assemble.assemble_microenvironment(domain, cells, fields, patient_id="reconstructed", metadata=None) -> Microenvironment`**
  Shape-checks and bundles the validated pieces into a `Microenvironment`.

### `tmesim.physics` — continuum fields

- **`diffusion.simulate_drug_transport(micro, drug, source, pde, velocity=None) -> TransportResult`**
  The core numerics; see §4.
- **`oxygen.solve_oxygen(micro, cfg) -> np.ndarray`** — steady-state
  reaction-diffusion O₂ field in `[0,1]`. **`oxygen.label_hypoxia(micro, oxygen, cfg) -> None`**
  sets `CellState.HYPOXIC` in place on cells below threshold.
- **`pressure.solve_ifp(micro) -> (ifp, velocity)`** — interstitial fluid pressure
  (mmHg) and the derived convective velocity field `(2, ny, nx)`; see §4.3.

### `tmesim.cells` — agent-based response

- **`agents.update_cell_states(micro, oxygen, cfg) -> None`** — relabel
  VIABLE/HYPOXIC/QUIESCENT from local oxygen, in place, on non-dead cells.
- **`pkpd.apply_drug_kill(micro, transport, pd, drug) -> np.ndarray`** — the
  pharmacodynamic kill model; see §5. Marks killed cells DEAD in place and returns
  the per-cell kill-probability array for diagnostics. Deterministic given the
  engine's RNG.

### `tmesim.treatment` — dosing and delivery

- **`drugs.get_drug(name) -> DrugParams`** — `KeyError`-safe lookup from
  `DRUG_LIBRARY`.
- **`drugs.plasma_profile(pk, times) -> np.ndarray`** — vascular (blood)
  concentration in µM at each time; multi-dose exponential decay.
- **`drugs.delivery_source(micro, drug, plasma_conc) -> np.ndarray`** — the
  `(ny, nx)` spatial delivery-rate source `= vascular_density × bbb_penetrance ×
  (plasma_conc − tissue)` proxy, scaled to `plasma_conc = 1`; the engine multiplies
  it by `plasma(t)`.

### `tmesim.sim` — the integrator

- **`engine.simulate_response(micro, cfg, rng=None) -> SimulationResult`** — the
  full per-patient pipeline (§3). Deterministic given `cfg.seed`.

### `tmesim.validation` — the kill-test

- **`baselines.baseline_mean(records) -> np.ndarray`** — constant prior.
- **`baselines.baseline_intrinsic_sensitivity(records) -> np.ndarray`** — mean
  per-cell sensitivity only, *no transport*.
- **`baselines.baseline_expression_logreg(train, test) -> np.ndarray`** — logistic
  regression on the `features` dict.
- **`metrics`** — `roc_auc`, `brier`, `calibration_curve(bins=10)`, `spearman`.
- **`experiment.run_stage0_experiment(records, cfg, out_dir="results") -> ExperimentResult`**
  Runs the sim per patient, scores it against every baseline (with cross-validation
  where a baseline needs training), and sets `beat_baseline`/`verdict` **honestly**:
  the sim must exceed the *best* baseline on `primary_metric = roc_auc` by a margin
  to pass. Writes a JSON.

### `tmesim.viz` — figures and report

- **`plots.plot_microenvironment(micro, out_path)`**, **`plots.plot_simulation(result, out_path)`**,
  **`plots.plot_validation(exp, out_path)`**.
- **`report.build_report(exp, sim_example, micro_example, out_dir="results") -> str`**
  — writes `results/report.md` embedding the figures and the honest verdict.

### `tmesim.cli`

`python -m tmesim.cli demo` (one synthetic patient + plots) · `run --n 40 --data
synthetic|real|auto` · `acquire` (attempt real data).

---

## 3. Core data model

All defined in `tmesim.core`, frozen, dependency-light.

| Type | Role |
|---|---|
| `Domain(nx, ny, dx=10.0, tissue_mask=None)` | The 2D finite-difference grid. `.shape=(ny,nx)`, `.extent_um`, `.zeros()`, `.ones()`, `.coords_um()`, `.to_index(x,y)`. `tissue_mask` (bool `(ny,nx)`) marks tissue for no-flux boundaries. |
| `CellType(IntEnum)` | `TUMOR, STROMAL, IMMUNE, ENDOTHELIAL, NECROTIC`. |
| `CellState(IntEnum)` | `VIABLE, HYPOXIC, QUIESCENT, DEAD`. |
| `CellPopulation(x, y, cell_type, state, sensitivity, exposure=None)` | Vectorized agents as parallel numpy arrays (fast for tens of thousands of cells). `sensitivity ∈ [0,1]` encodes biology transport cannot (MGMT, resistance clones). `exposure` accumulates local drug AUC. Helpers: `.n`, `.tumor_mask()`, `.alive_mask()`, `.viable_tumor_fraction()`, `.copy()`. |
| `TissueFields(ecm_density, vascular_density, oxygen=None, ifp=None, velocity=None)` | Continuum fields on the grid. ECM raises tortuosity (lowers effective diffusivity); vasculature is both drug and oxygen source. |
| `Microenvironment(domain, fields, cells, patient_id, metadata)` | The central object. `metadata` carries non-spatial features (subtype, `mgmt_methylated`, expression) used by baselines and the sensitivity model. `.summary()`. |
| `DrugParams` | Transport + PD parameters (diffusivity, MW, `ic50`, `hill`, `binding`, `clearance`, `bbb_penetrance`). See `DRUG_LIBRARY` for temozolomide-/doxorubicin-/cisplatin-like defaults. |
| `PKConfig, PDEConfig, OxygenConfig, PDConfig, SimConfig` | Stage configs, bundled by `SimConfig` (which also holds `seed`). |
| `TransportResult` | `auc (ny,nx)`, `final_concentration`, `times`, `mean_concentration_history`, `penetration_index ∈ [0,1]`. |
| `SimulationResult` | `predicted_response ∈ [0,1]` (the scored number), `viable_tumor_fraction`, `kill_fraction`, `viable_history`, `times`, `transport`, `fields` (named arrays for plotting), `extras`. |
| `PatientRecord` | `patient_id`, `micro`, `observed_response` (binary label or continuous −log IC50; the harness handles both), `features` (flat scalars for baselines). |
| `ExperimentResult` | `model_scores`, `baseline_scores`, `per_patient`, `primary_metric="roc_auc"`, `verdict`, `beat_baseline`. |

**The primary simulation readout** is `viable_tumor_fraction` — the fraction of the
*original* tumor cells still alive (1.0 = no kill, 0.0 = complete kill). The engine
maps `kill_fraction = 1 − viable_tumor_fraction` to `predicted_response`, and that
prediction is what the validation harness scores.

---

## 4. The physics

The physical model is intentionally the *minimal* mechanistic set that captures why
delivery fails in GBM: drug transport limited by tortuosity and clearance,
oxygen-driven resistance, and pressure-opposed inward flux.

### 4.1 Drug transport — convection–diffusion–reaction

`physics.diffusion.simulate_drug_transport` solves for the drug concentration
field `C(x, y, t)`:

```
∂C/∂t = ∇·(D_eff ∇C) − ∇·(v C) + S − (clearance + binding) C
```

- **Effective diffusivity with tortuosity.** `D_eff = D·(1 − α·ecm_density)`. Dense
  extracellular matrix raises tortuosity and lowers the effective diffusion
  coefficient — drug moves more slowly through matrix-rich tissue.
- **Convection.** The advective term `−∇·(vC)` carries drug with the interstitial
  velocity field `v` from `pressure.solve_ifp` (optional, gated by
  `PDEConfig.convection`).
- **Source.** `S` is the `(ny, nx)` delivery-rate field from `treatment.delivery_source`
  (`vascular_density × bbb_penetrance × plasma proxy`), scaled by the plasma time
  profile `plasma(t)`.
- **Reaction sink.** `(clearance + binding)·C` captures metabolism/clearance and
  tissue binding/uptake as first-order loss.
- **Boundaries.** No-flux at the domain edges and at the `tissue_mask` boundary.
- **Stability.** An unconditionally stable implicit/ADI scheme is preferred;
  otherwise a CFL-safe `dt` is enforced with a clear error if violated.

The integral `∫ C dt` is accumulated into `TransportResult.auc` — the exposure
field the kill model consumes — along with a `penetration_index` summarizing how
far the drug reached relative to a well-mixed ideal.

### 4.2 Oxygen — steady-state reaction–diffusion

`physics.oxygen.solve_oxygen` solves for a steady-state O₂ field in `[0,1]`:

```
0 = D_O2 ∇²O − consumption·(tumor-cell density) + supply·(vascular_density)
```

Vasculature supplies oxygen; tumor cells consume it. Because the vascular network
is sparse in the tumor core, the core goes hypoxic. `label_hypoxia` then flags
cells sitting in `oxygen < hypoxia_threshold` as `HYPOXIC` — a state the kill model
treats as drug-resistant. This is the mechanism by which the *same* poorly-perfused
region is simultaneously under-dosed and intrinsically more resistant.

### 4.3 Interstitial fluid pressure — the Jain effect

`physics.pressure.solve_ifp` solves a Darcy / poroelastic-lite Poisson problem for
interstitial fluid pressure `p` (mmHg), elevated in the tumor core, and returns the
convective velocity `v = −K ∇p`. Elevated core pressure produces an *outward*
velocity that **reduces net inward drug delivery** — the classic Jain effect, a
first-order reason systemic agents fail to accumulate in solid tumor cores. This
velocity feeds back into the transport advection term (§4.1).

---

## 5. The pharmacodynamic kill model

`cells.pkpd.apply_drug_kill` converts exposure into cell death. For each cell it
samples the local drug AUC from `transport.auc` at the cell's grid index, then maps
exposure to a **kill probability** through a Hill function:

```
p_kill = max_kill_rate · sensitivity · resistance · [ AUC^h / (IC50^h + AUC^h) ]
```

- **Hill dose–response** on accumulated exposure, using the drug's `ic50` and `hill`
  coefficient.
- **Per-cell `sensitivity`** scales the effect — this is where intrinsic biology
  (MGMT status, resistance clones) enters, independent of delivery.
- **Resistance modulation** from cell state: hypoxic cells are multiplied by
  `PDConfig.hypoxia_resistance` (< 1) and quiescent cells by
  `PDConfig.quiescent_resistance`, so poorly-oxygenated regions resist kill even
  where drug reaches them.

Killed cells are marked `DEAD` in place; the kill is applied over the treatment
course and is deterministic given the engine's seeded `numpy.default_rng`.

### The per-patient engine loop

`sim.engine.simulate_response` executes, in order:

1. Oxygen field + hypoxia labeling (`physics.oxygen`).
2. IFP + velocity (`physics.pressure`).
3. Delivery source from the plasma profile (`treatment.delivery_source`).
4. Drug-transport PDE → `TransportResult` (`physics.diffusion`).
5. Pharmacodynamic kill over the course (`cells.pkpd`).
6. `predicted_response ∈ [0,1]` computed from `kill_fraction` via a documented map.

`SimulationResult.fields` is populated with oxygen / IFP / AUC / vascular arrays for
plotting. The whole loop is deterministic given `cfg.seed`.

---

## 6. Honest limitations

These are load-bearing caveats, stated plainly because the value of a Stage 0
kill-test depends on not overclaiming.

- **Parameterization is uncalibrated.** Diffusivities, BBB penetrance, IC50s, and
  oxygen/consumption/supply rates are order-of-magnitude biologically plausible,
  not fit to patient data. Which of these the verdict actually depends on is a
  Stage 1 sensitivity-analysis question, not something this scaffold answers.
- **The model is 2D.** Reconstructing a faithful 3D microenvironment from sparse
  samples is an unsolved inference problem, so 2D is a deliberate dodge (a
  whole-slide section or organoid cross-section is naturally 2D). Real 3D geometry
  — vascular connectivity, out-of-plane transport — is not represented.
- **Synthetic ground truth is not independent evidence.** On synthetic cohorts, the
  `observed_response` is generated from the same broad biological assumptions the
  simulator encodes (transport × sensitivity, plus noise). A synthetic pass shows
  the machinery is correct and the signal is *recoverable when present* — it is a
  test of the code and the harness, **not** of real-world predictive power. The
  noise is calibrated so the answer is not trivially recoverable, which keeps the
  test fair but does not make it clinical.
- **The validation gap is the whole point.** Beating baselines on synthetic data is
  necessary but not sufficient. A meaningful result requires *controlled-access,
  paired sample→outcome data* (retrospective at Stage 3, prospective at Stage 5).
  Until then, `beat_baseline = True` on a synthetic cohort means "the machinery
  works," not "the model predicts patient response."
- **Reduced biology.** Immune dynamics, adaptive resistance, drug combinations,
  and treatment scheduling beyond a simple multi-dose PK profile are out of scope
  for Stage 0 and would be added only if a real-data residual demands them.

The harness is built to *fail loudly*. A null or negative verdict on this pipeline
is a legitimate, publishable outcome — and reporting it honestly is the point of
running Stage 0 before spending on controlled-access data.
