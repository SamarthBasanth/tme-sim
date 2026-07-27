# tmesim — Architecture

How data flows through `tmesim` in the **real build**: what each module does, the
core data model, the physics, and the honest limitations. It is a companion to the
frozen public API in [`../ARCHITECTURE_CONTRACT.md`](../ARCHITECTURE_CONTRACT.md)
and the molecular-layer contract in
[`../MOLECULAR_CONTRACT.md`](../MOLECULAR_CONTRACT.md); where any of them disagree,
the contracts win. Every equation and its provenance is in
[`MATHEMATICAL_MODEL.md`](MATHEMATICAL_MODEL.md).

**Real data only.** `tmesim/data/synthetic.py` (the synthetic cohort generator) and
the fabricated-geometry loader are deleted. Nothing in the pipeline can produce an
invented patient. The synthetic-era scripts are retired to
[`../archive/synthetic_era/`](../archive/synthetic_era/README.md).

Conventions throughout: length in **micrometers (µm)**, time in **seconds**,
concentration in **µM**, diffusivity in **µm²/s**. All grid fields are 2D numpy
arrays of shape `(ny, nx)`, indexed `[row = y, col = x]` (image/numpy convention).

**Production resolution is `dx = 12 µm`.** The calibrated perivascular boundary
layer is ~45 µm thick; coarser grids do not resolve it. See §7.

---

## 1. Pipeline data flow

```
 data/real/  ──►  reconstruct/real_*  ──►  Microenvironment  ──►  physics/ + cells/  ──►  sim/  ──►  validation/
 (real bytes)      (real whole niche)       (core object)          (transport + kill)      (engine)   (vs real PK)
                                                  ▲                                                        │
                                        molecular/ (structure → DrugParams)                            viz/ + results/
```

The `Microenvironment` (in `tmesim.core`) is the hinge. It is the **output** of
reconstruction and the **input** to simulation; every other module produces one,
consumes one, or transforms fields inside one. Because the object and all result
containers are frozen in `core`, the physics, cell, simulation and validation
modules develop against a fixed contract.

Two orthogonal inputs meet at the engine:

* the **tissue** — real cells at real coordinates with real vascular / ECM / oxygen
  fields, from `reconstruct.real_pipeline.reconstruct_real`;
* the **drug** — a `DrugParams` compiled from a verified molecular structure and a
  *measured* GDSC/literature IC50, via `molecular.compile_drug.molecule_to_drugparams`
  plus `validation.calibration.calibrated_drugparams` and the binding-site barrier.

Nothing in the engine knows about genes, drugs by name, or datasets. That
separation is what lets the reconstruction and the pharmacology be validated
independently of each other.

---

## 2. Module responsibilities

### `tmesim.data.real` — the real caches

Each loader reads a genuinely downloaded open dataset and returns plain
numpy/pandas. Provenance for every byte is in [`real_data/`](real_data/).

- **`singlecell_spatial.load(path=None, ...) -> dict`** — the 10x **Xenium** GBM
  section: 40,887 segmented cells, measured `(x, y)` centroids in µm, a 319-gene
  panel count matrix, and the panel spec. Cache: `data_cache/xenium_gbm/`.
- **`visium.list_samples() / sample_info(id) / load_sample(id) -> dict`** — the
  four real **Visium** GBM sections (`10x_parent_gbm`, `GBM1`, `GBM2`, `GBM3`):
  MEX matrix + `spatial/` positions + scalefactors, converted to µm. Cache:
  `data_cache/visium_gbm/`.
- **`scref_gbm.load_reference(...)`** — the **Darmanis GSE84465** scRNA-seq atlas
  (23,460 genes × 3,589 annotated cells, 7 cell types): the deconvolution
  reference. Cache: `data_cache/scref_gbm/`.
- **`ivygap.load_structures(...) -> IvyGapStructures`** and
  **`ivygap.structure_signatures(...)`** — **Ivy GAP** anatomic-structure RNA-seq
  (270 samples, 37 tumors) mapped to coarse niche regions (margin / core /
  perivascular / necrotic). Cache: `data_cache/ivygap/`.
- **`drug_response.load_hgcc(...) / load_gdsc(...) / summarize()`** — the measured
  potency layer: **HGCC** (262 drugs × 85 patient-derived cultures) and **GDSC2**
  glioma dose-response. `SOURCES` also records the ranked *gated* datasets that
  would close the response-validation loop.

### `tmesim.data` — acquisition and the base spatial reader

- **`acquire.acquire_all(cache_dir="data_cache", timeout=60) -> dict`** — the only
  module permitted network access. Best-effort live download of TCGA-GBM clinical
  and expression (GDC REST), the Ivy GAP expression zip, and GDSC2 glioma
  dose-response. Never crashes; returns `{source: {status, path, note}}`.
- **`spatial_ingest`** — the original single-sample Visium reader that
  `data.real.visium` generalises.

### `tmesim.molecular` — structure → `DrugParams`

Verified structures and SMILES (`library`), RDKit-or-fallback descriptors
(`descriptors`), Stokes–Einstein hydrodynamic radius (`stokes_einstein`), Brownian
/ Ogston–Renkin hindered diffusion through an ECM mesh (`brownian`, `ecm_mesh`),
BBB permeability rules (`permeability`), measured GDSC potency (`potency`), tissue
binding (`binding`) and the saturable antibody binding-site barrier
(`binding_site_barrier`).

- **`compile_drug.build_profile(molecule, ecm_reference=0.4, ic50=None, ...) -> MolecularProfile`**
- **`compile_drug.molecule_to_drugparams(molecule_or_profile, ...) -> DrugParams`**
  — maps `effective_diffusivity → DrugParams.diffusivity`, the measured IC50, the
  structure-derived `bbb_penetrance`, and tissue binding.

### `tmesim.reconstruct` — real tissue → real whole niche

The entry point:

- **`real_pipeline.reconstruct_real(source='xenium'|'visium', sample=None, dx_um=30.0) -> Microenvironment`**
  Dispatches to the Xenium or Visium path. *Note:* the Visium path snaps `dx` to
  the ~100 µm spot pitch when `dx < 50`; experiments that need matched-resolution
  Visium (e.g. `experiments/cross_tissue.py`) drive the deconvolution and field
  construction directly so they can simulate on the dx = 12 µm grid.
- **`real_pipeline.apply_drug_sensitivity(micro, drug_name, sample=None) -> dict`**
  Overwrites `micro.cells.sensitivity` in place with the DRUG-SPECIFIC per-cell
  value (§6). Lives here, not in the engine, on purpose: `core`, `physics` and
  `sim.engine` stay drug- and gene-agnostic.
- **`real_pipeline.reconstruction_report(micro) -> dict`** — cell-type fractions,
  mean vascular/ECM/oxygen over the tissue mask, the grid, the ECM marker
  transcripts used, and the **sensitivity model in force** (every absolute kill
  fraction is conditional on it, and the report says so).

Xenium path (single-cell resolution):

- **`real_xenium.phenotype_cells(...)`** — **threshold-gated** marker typing into
  `TUMOR / STROMAL / IMMUNE / ENDOTHELIAL` (TUMOR: EGFR/SOX2/OLIG2/MKI67;
  ENDOTHELIAL: PECAM1 …). The gate replaced a marker-argmax classifier that
  inflated ENDOTHELIAL to 20.6% via single-marker PECAM1 bias; the gated result,
  9.7%, is corroborated independently by the Visium deconvolution's 9.77%.
- **`marker_panels.select_marker_map(genes)`** — panel-aware marker maps. The
  hardcoded map was written for the 319-gene Human Brain panel; a second cached
  section (`data_cache/xenium_gbm_io`, the 480-gene Immuno-Oncology panel) carries
  none of PECAM1/SOX2/OLIG2/AQP4/GJA1/MOG, on which the brain map would collapse
  the reconstruction to TUMOR + IMMUNE with no vasculature and therefore no oxygen
  source. Selection is by panel content, and the 319-gene behaviour is unchanged
  bit-for-bit.
- **`real_xenium.tissue_mask_from_cells(...)`** — `Domain.tissue_mask` as the union
  of fixed-radius physical discs around measured centroids (`geometric` mode), so
  the simulated domain is **independent of `dx`** — which is what makes the grid
  ladder in `numerical_rigor` a clean convergence test.
- **`real_xenium.per_cell_sensitivity(expr, genes, drug) -> (array, info)`** — the
  measured determinant-gene model of §6.
- **`real_xenium.reconstruct_from_xenium(sample_path=None, dx_um=30.0, ...) -> Microenvironment`**

Visium path (spot resolution):

- **`deconvolve.build_signatures(...)` / `deconvolve_spot(...)` / `deconvolve_sample(sample_id) -> (coords_um, fractions, type_names)`**
  — NNLS deconvolution of each spot against the Darmanis atlas signatures.
- **`real_visium.reconstruct_from_visium(sample_id, dx_um=100.0) -> Microenvironment`**
  — places cells at continuous `(x, y)` inside each spot footprint in the
  deconvolved proportions.

Fields (both paths):

- **`real_fields.vascular_density_from_cells(domain, cells, bandwidth_um=40)`** —
  Gaussian KDE of the **measured ENDOTHELIAL cell positions**, masked to tissue and
  divided by its own maximum. It therefore carries *geometry*, not an absolute
  perfused fraction.
- **`real_fields.ecm_density_from_positions(domain, x, y, weights=None, bandwidth_um=40)`**
  — the production ECM field: an **expression-weighted** KDE of 11 real matrix
  transcripts (COL1A2, COL12A1, COL25A1, POSTN, VCAN, DCN, FBLN1, LAMA2, THBS1,
  SPON1, TGFBI). `ecm_density_from_cells` is the stromal-cell proxy used where
  per-cell transcripts are unavailable (Visium).
- **`real_fields.build_real_fields(...) -> TissueFields`** — assembles ECM +
  vasculature and leaves `oxygen = None` deliberately, so
  **`real_fields.add_oxygen(micro, cfg)`** can solve it with the *unmodified*
  `physics.oxygen.solve_oxygen`. Physics stays in the physics module.

Niche regions:

- **`region_label.region_signatures(...)` / `label_regions(...)` / `region_field(...)`**
  — assign each location to its best-correlating **Ivy GAP** anatomic region
  (margin / core / perivascular / necrotic) from real per-structure signatures.
  Visium shares 283/311 signature genes with the panel; Xenium only 41/311, so
  region calls on Xenium rest on a much thinner overlap and are treated as weaker.

**Legacy, not on the real path:** `segmentation.py` (image → phenotyped cells),
`spatial.py` and `vasculature.py` (the synthetic field/vessel generators, now only
reachable through `data/spatial_ingest.py`'s fallback), and `assemble.py` (still
used by `real_visium` for shape-checked bundling).

### `tmesim.physics` — continuum fields

- **`diffusion.simulate_drug_transport(micro, drug, source, pde, velocity=None) -> TransportResult`**
  — the core numerics; see §4.1.
- **`oxygen.solve_oxygen(micro, cfg) -> np.ndarray`** — steady-state
  reaction–diffusion O₂ field in `[0,1]`.
  **`oxygen.label_hypoxia(micro, oxygen, cfg)`** — flags `CellState.HYPOXIC` in place.
- **`pressure.solve_ifp(micro) -> (ifp, velocity)`** — interstitial fluid pressure
  (mmHg) and the derived convective velocity `(2, ny, nx)`; see §4.3.

### `tmesim.cells` — agent-based response

- **`agents.update_cell_states(micro, oxygen, cfg)`** — relabel VIABLE / HYPOXIC /
  QUIESCENT from local oxygen, in place. **`agents.state_fractions(micro)`** reports
  the banding.
- **`pkpd.apply_drug_kill(micro, transport, pd, drug) -> np.ndarray`** — the Hill
  kill model (§5); marks cells DEAD in place, returns per-cell kill probabilities.

### `tmesim.treatment` — dosing and delivery

- **`drugs.get_drug(name) -> DrugParams`**, **`drugs.plasma_profile(pk, times)`**,
  **`drugs.delivery_source(micro, drug, plasma_conc) -> np.ndarray`** — the
  `(ny, nx)` delivery-rate field `vascular_density × bbb_penetrance × plasma proxy`.

### `tmesim.sim` — the integrator

- **`engine.simulate_response(micro, cfg, rng=None) -> SimulationResult`** — the
  full per-tissue pipeline (§5), run on a private copy of `micro`. Deterministic
  given `cfg.seed`.

### `tmesim.validation` — the honesty layer

This is where the real build spends most of its code. Nothing here modifies
`core` / `physics` / `molecular` / `reconstruct` / `sim.engine`; every module only
*drives* them and scores the output against a real measurement.

| Module | What it checks | Report |
|---|---|---|
| `calibration` | fits the single α on doxorubicin's real anchor; re-validates out-of-sample | `results/calibration.md` |
| `pk_anchors` | the 7 PMID-verified tissue-PK anchors and the `d_half ↔ L` conversion | — |
| `pk_multidrug` | predicted vs real cited penetration bands for 6 drugs | `results/pk_multidrug.md` |
| `pk_reality` | predicted penetration vs published tissue PK, end to end | `results/pk_reality.md` |
| `glioma_regional_pk` | predicted tissue:plasma partition vs measured human glioma regional PK | — |
| `reconstruction_fidelity` | cross-platform composition concordance, permutation tests on niche labels, deconvolution stability | `results/reconstruction_fidelity.md` |
| `vascular_validation` | reconstructed vessel geometry vs published human-GBM microvascular morphometry | — |
| `cellstate_validation` | falsifies the oxygen-banded cell-state model against measured per-cell MKI67/TOP2A | — |
| `delivery_metrics` | the five **defensible** delivery metrics that replace `penetration_index` | `results/delivery_metrics.md` |
| `vessel_profile` | AUC vs distance-from-vessel fits, with a SOURCE NULL control | `results/vessel_profile.md` |
| `numerical_rigor` | mass conservation, dx ladder, dt ladder, parameter elasticities | `results/numerical_rigor.md` |
| `uncertainty` | Latin-hypercube propagation of the reconstruction knobs + α | `results/uncertainty.md` |
| `pharmacology_real` | measured GDSC vs measured HGCC (no model in the loop) | `results/pharmacology_real.md` |

**Retired:** `experiment.py` (the Stage-0 cohort harness), `baselines.py` and
`metrics.py` scored a synthetic cohort against cheap non-spatial baselines. They
still import cleanly but nothing on the real path invokes them — there is no real
cohort with paired outcomes to run them on. Their conclusions are preserved in
`results/STAGE0_FINDINGS.md`.

### `tmesim.viz` and `tmesim.cli`

`viz.plots` (microenvironment / simulation figures), `viz.molecular_plots`,
`viz.explainer` (the multiscale explainer figures), `viz.patient_report`.

`python -m tmesim.cli {reconstruct, simulate, molecules, molecular-run, acquire}`
— see `--help`, and [`REPRODUCE.md`](REPRODUCE.md) for tested invocations.

---

## 3. Core data model

All defined in `tmesim.core`, frozen, dependency-light.

| Type | Role |
|---|---|
| `Domain(nx, ny, dx=10.0, tissue_mask=None)` | The 2D finite-volume grid. `.shape=(ny,nx)`, `.extent_um`, `.coords_um()`, `.to_index(x,y)`. `tissue_mask` marks tissue for no-flux boundaries; on the real path it is built geometrically so it does not move with `dx`. |
| `CellType(IntEnum)` | `TUMOR, STROMAL, IMMUNE, ENDOTHELIAL, NECROTIC`. |
| `CellState(IntEnum)` | `VIABLE, HYPOXIC, QUIESCENT, DEAD`. |
| `CellPopulation(x, y, cell_type, state, sensitivity, exposure=None)` | Vectorized agents as parallel numpy arrays — on the Xenium path these are **40,887 measured cells at measured coordinates**. `sensitivity ∈ [0,1]` encodes biology transport cannot; `exposure` accumulates local drug AUC. |
| `TissueFields(ecm_density, vascular_density, oxygen=None, ifp=None, velocity=None)` | Continuum fields. ECM raises tortuosity; vasculature is both the drug and the oxygen source. |
| `Microenvironment(domain, fields, cells, patient_id, metadata)` | The central object. `metadata` carries `real_pipeline` provenance (source, sample, markers used) and `sensitivity_model_applied`. |
| `DrugParams` | Transport + PD parameters (`diffusivity`, `molecular_weight`, `ic50`, `hill`, `binding`, `clearance`, `bbb_penetrance`), produced by the molecular layer, not hand-written. |
| `PKConfig, PDEConfig, OxygenConfig, PDConfig, SimConfig` | Stage configs, bundled by `SimConfig` (which holds `seed`). |
| `TransportResult` | `auc (ny,nx)`, `final_concentration`, `times`, `mean_concentration_history`, `penetration_index` (**see §7 — do not rank drugs by it**). |
| `SimulationResult` | `predicted_response ∈ [0,1]`, `viable_tumor_fraction`, `kill_fraction`, `viable_history`, `times`, `transport`, `fields`, `extras`. |

`PatientRecord` and `ExperimentResult` remain in `core` for the retired cohort
harness; nothing on the real path emits them.

---

## 4. The physics

The minimal mechanistic set that captures why delivery fails in GBM: transport
limited by tortuosity and consumption, oxygen-driven resistance, and
pressure-opposed inward flux. Full derivations in
[`MATHEMATICAL_MODEL.md`](MATHEMATICAL_MODEL.md) §3–5.

### 4.1 Drug transport — convection–diffusion–reaction

```
∂C/∂t = ∇·(D_eff ∇C) − ∇·(v C) + S − (clearance + binding) C
```

- **Effective diffusivity with tortuosity.** `D_eff = D·(1 − 0.6·ecm_density)`, with
  `D` itself compiled from the molecule (Stokes–Einstein + Ogston–Renkin mesh
  hindrance), not assumed.
- **Convection.** `−∇·(vC)` with `v` from `pressure.solve_ifp`, gated by
  `PDEConfig.convection`. On this tissue interstitial speeds are < 0.7 µm/s and
  change the penetration index by < 3e-4, so the audits run with advection off to
  isolate the diffusion–reaction scheme.
- **Source.** `S = vascular_density × bbb_penetrance × plasma(t)` — i.e. the drug
  enters through the **real endothelial-derived vascular field**.
- **Reaction sink.** `(clearance + binding)·C`, where the binding term carries the
  calibrated perivascular consumption (§7).
- **Boundaries.** No-flux at the domain edge and at the `tissue_mask` interface.
- **Scheme.** Implicit (backward-Euler) finite volume — unconditionally stable, and
  **exactly conservative by construction** (zero row sums). Verified on the real
  irregular mask: relative mass error ~1e-13, concentrations never negative.

`∫C dt` accumulates into `TransportResult.auc` — the exposure field the kill model
consumes.

### 4.2 Oxygen — steady-state reaction–diffusion

```
0 = D_O2 ∇²O − consumption·(tumor-cell density) + supply·(vascular_density)
```

Solved on the **real** vascular field, so hypoxia appears where the measured
endothelium is sparse. `label_hypoxia` bands cells by `o < t` (HYPOXIC),
`t ≤ o < 2t` (QUIESCENT), `o ≥ 2t` (VIABLE) with `t = OxygenConfig.hypoxia_threshold`.
The resulting hypoxic fraction is **grid-sensitive** — see §7.

### 4.3 Interstitial fluid pressure — the Jain effect

`pressure.solve_ifp` solves a Darcy/poroelastic-lite Poisson problem for `p`
(mmHg), elevated in the tumor core, and returns `v = −K ∇p`. Elevated core pressure
produces an outward velocity that reduces net inward drug delivery.

---

## 5. Pharmacodynamics and the engine loop

`cells.pkpd.apply_drug_kill` samples `transport.auc` at each cell's grid index and
maps exposure to a kill probability:

```
p_kill = max_kill_rate · sensitivity · resistance · [ AUC^h / (IC50^h + AUC^h) ]
```

with `ic50` and `hill` from the **measured** potency layer, `sensitivity` per §6,
and `resistance` the oxygen-state multiplier (`PDConfig.hypoxia_resistance`,
`quiescent_resistance`). Those multipliers are load-bearing: the patient report's
own counterfactual shows they cost every panel drug ~50% of its predicted kill —
which is why `validation.cellstate_validation` exists to falsify them against
measured per-cell MKI67/TOP2A.

`sim.engine.simulate_response` executes, in order:

1. Oxygen field + hypoxia labeling (`physics.oxygen`).
2. IFP + velocity (`physics.pressure`).
3. Delivery source from the plasma profile (`treatment.delivery_source`).
4. Drug-transport PDE → `TransportResult` (`physics.diffusion`).
5. Pharmacodynamic kill over the course (`cells.pkpd`).
6. `predicted_response ∈ [0,1]` from `kill_fraction`.

Deterministic given `cfg.seed`; only step 5 draws from the RNG, and the numerical
audits average it over seeds 0–3.

---

## 6. Per-cell sensitivity — measured where measurable, uniform where not

`reconstruct.real_pipeline.apply_drug_sensitivity` sets `cells.sensitivity` from a
**measured determinant transcript on the 319-gene Xenium panel**, and only where
the panel carries one:

| drug(s) | determinant | why |
|---|---|---|
| doxorubicin, etoposide, mitoxantrone | **TOP2A** | TOP2 poisons trap the enzyme–DNA cleavage complex, so cytotoxicity scales with TOP2A level (PMID 18574145); TOP2A-low cells are resistant |
| cisplatin, lomustine, carmustine, methotrexate, vincristine | **MKI67** | proxy for cell-cycle dependence |
| **temozolomide**, procarbazine | **none — uniform 1.0** | the canonical determinant is **MGMT**, which is not on the panel, and methylation is not measurable from these transcripts |
| bevacizumab | **none — uniform 1.0** | anti-VEGF antibody; no direct tumor-cell cytotoxicity to gate |

Panel reality check: PRESENT — MKI67, TOP2A, EGFR, TP53, AQP4, PCNA, CDK1.
ABSENT — MGMT, ABCB1, ABCG2, CDKN2A, ATM, MSH2, MLH1, PARP1.

The **default reconstruction is drug-agnostic** (uniform 1.0 for every drug), so
drug *rankings* are unaffected by the sensitivity model and absolute kills are an
equal no-resistance upper bound. Applying the drug-specific model makes absolute
kill fractions **conditional** on it — never a measured response. Visium
reconstructions have no per-cell transcripts, so sensitivity stays uniform 1.0 and
that is reported.

---

## 7. Resolution, and what `penetration_index` actually measures

Both facts must be carried together; either alone is misleading.

**(a) `penetration_index` is not converged on coarse grids.** The calibrated
doxorubicin penetration length is ~45 µm, so a grid with `dx > ~45 µm` cannot
resolve the perivascular boundary layer and reports an artificially flat (high)
index. On the real Xenium section it falls monotonically
0.285 (dx 60) → 0.2666 (30) → 0.1991 (15) → **0.1934 (12)** → 0.1869 (10). Values
quoted at dx = 30 µm sit ~15–25% high; values at dx = 100 µm are far outside the
converged regime entirely.

**(b) Even converged, it is a scale-invariant shape statistic.**
`penetration_index = mean(AUC) / p99(AUC)` is a ratio: multiply the entire exposure
field by any constant and it does not move. At dx = 12 µm it therefore reports the
**vascular geometry**, identical for every drug on one section — it collapses to
**0.194–0.196 across the whole five-drug panel** while their delivered exposure
spans **~122×**. Correlating it against vascular density is close to circular.

**Consequence.** Compare drugs by **`mean_tissue_auc` / `mean_tumor_auc`** (µM·s,
absolute) or **`therapeutic_coverage` at 1× IC50** (the fraction of tumor cells
whose time-averaged exposure `E = AUC/T` clears the drug's own IC50 — dimensionally
sound, grid-converged, and genuinely drug-discriminating). Both are in
`validation.delivery_metrics`. The earlier ranking "temozolomide penetrates best,
doxorubicin barely reaches the tumor bulk" is **retracted**, as is the cross-tissue
"kill spans 2% → 24%" headline that mixed a dx = 12 µm Xenium run with dx = 100 µm
Visium runs.

`kill_fraction` **is** grid-converged (plateaus by dx ≲ 20 µm; 0.0035 relative
movement between the two finest grids).

**The hypoxic fraction is grid-sensitive and only order-of-magnitude.** The same
Visium tissue reports 100% hypoxic at dx = 100 µm and 0% at dx = 12 µm; the Xenium
section reported ≈94% coarse and ≈13% at dx = 12 µm. Earlier attribution of small
predicted kills to hypoxic resistance was wrong — they are driven mainly by limited
delivered exposure versus these drugs' IC50s.

---

## 8. Calibration — one constant, the e-fold convention

The tissue loss rate is modelled as `k_eff = k_passive + α · uptake_proxy(drug)`,
with **α the single free constant in the whole pipeline**. It is fit **once** to
doxorubicin's real perivascular anchor (Primeau 2005, PMID 16361566) and then held
fixed for every other drug and every tissue.

The anchor is published as a **half-distance** `d_half`; the simulator's
penetration length `L = √(D_eff / k_eff)` is an **e-folding length**:

```
d_half = L · ln 2        =>        L = d_half / ln 2 ≈ 1.443 · d_half
```

Fitting α directly to the raw half-distance skips that conversion and inflates α by
exactly `1/ln²2 = 2.081×`. **That is the origin of the old α = 309.6; the correct
value is α ≈ 148.6** (`calibrate_alpha()` → 148.608 against a target e-folding
length of 64.92 µm). `calibrate_alpha` now takes an e-folding target so the mistake
cannot recur silently. Any report still quoting 309.6 is stale.

With α frozen, the model **predicts** out-of-sample: mitoxantrone 48.2 µm and
topotecan 83.6 µm (a chemically distinct class), both inside their real cited
bands.

---

## 9. Honest limitations

- **The response axis is unvalidated.** Reconstruction fidelity and transport
  physics are validated against real measurements; **predicted drug response
  against real patient outcomes is not**. Every open spatial GBM dataset is an
  untreated resection with no paired response, and every open drug-response
  dataset is dissociated cells with no microenvironment. This is a data-access
  gap (HTAN + dbGaP/DUA, PDMR/organoid MTA), not an engineering gap.
- **2D.** A section is naturally 2D, but out-of-plane vascular connectivity is
  invisible: a vessel just above the plane makes a point look "far" when it is not.
  Distances are Euclidean, not geodesic around the irregular mask.
- **Fields are inferred, not measured.** Vascular density, ECM and oxygen are
  reconstructed from real markers and real perfusion signal. `vascular_density` is
  max-normalised, so it encodes geometry, not an absolute perfused fraction, and is
  not comparable across tissues as an absolute.
- **One anchor sets the delivery scale.** α's own measurement uncertainty is
  plausibly ±50–100%. It has the *smallest* local elasticity of the tuned knobs
  (|E_pen| ≈ 0.11, |E_kill| ≈ 0.20 over ±20%) but the widest prior range, so it is a
  leading *global* uncertainty.
- **Reconstruction gates dominate local sensitivity.** Endothelial typing threshold
  (|E_pen| ≈ 0.52) and KDE bandwidth (|E_kill| ≈ 0.44) move the prediction more than
  the solver does. The mesh is not the bottleneck once dx ≲ 12 µm.
- **`bbb_penetrance` is a structure-only ordering**, not a partition ratio — only
  partially concordant with measured brain:plasma ratios (Spearman 0.50), because
  there is no efflux or carrier transport in the model.
- **Simplified systemic PK.** A single time-averaged sustained plasma level drives
  the transport solve; it preserves total delivered drug and relative ranking but
  is not a time-resolved infusion.
- **Reduced biology.** No adaptive resistance, no immune dynamics, no explicit
  drug-combination synergy beyond spatial complementarity.
- **One section per tumor, no replicates.** Where a claim rests on a single
  section, the write-ups say so.

The design principle behind all of the above: **the modules that could flatter the
result and the modules that check it are separate, and the checking modules are
allowed to fail loudly.** Retractions in this repo (the ENDOTHELIAL argmax bug, the
penetration-index ranking, the cross-tissue spread, the α unit error) were all
found that way.
