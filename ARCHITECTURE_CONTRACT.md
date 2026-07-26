# tmesim — Architecture Contract (v1)

**Read this fully before writing code.** It is the frozen interface every module
builds against. The `tmesim.core` package is already implemented and MUST NOT be
modified. Everything else consumes/produces `core` types with the exact
signatures below. Do not invent alternative return types — use the `core.results`
containers.

## The pipeline (data flow)

```
 data/  ──►  reconstruct/  ──►  Microenvironment  ──►  physics/ + cells/  ──►  sim/  ──►  validation/
 (sample)     (build TME)        (core object)          (transport+kill)      (engine)   (beat baseline?)
                                                                                            │
                                                                                        viz/ + report
```

The **Stage 0 kill-test** we are proving: *does a mechanistic microenvironment +
drug-transport simulation predict observed drug response better than a cheap
non-spatial baseline?* If it does not beat the baseline, the idea is not ready —
so the validation harness must be honest and adversarial, not flattering.

Cancer focus: **glioblastoma** (drug delivery / diffusion is the crux; BBB and
infiltrative margins make the microenvironment matter most).

## Core types (already implemented in `tmesim.core`)

- `Domain(nx, ny, dx=10.0, tissue_mask=None)` — 2D grid, fields are `(ny, nx)`.
  `.shape`, `.extent_um`, `.zeros()`, `.ones()`, `.coords_um()`, `.to_index(x,y)`.
- `CellType(IntEnum)`: TUMOR, STROMAL, IMMUNE, ENDOTHELIAL, NECROTIC.
- `CellState(IntEnum)`: VIABLE, HYPOXIC, QUIESCENT, DEAD.
- `CellPopulation(x, y, cell_type, state, sensitivity, exposure=None)` — parallel
  numpy arrays. `.n`, `.tumor_mask()`, `.alive_mask()`, `.viable_tumor_fraction()`, `.copy()`.
- `TissueFields(ecm_density, vascular_density, oxygen=None, ifp=None, velocity=None)`.
- `Microenvironment(domain, fields, cells, patient_id="synthetic", metadata={})` — `.summary()`.
- `DrugParams(name, diffusivity, molecular_weight, ic50, hill, binding, clearance, bbb_penetrance)`.
- `PKConfig, PDEConfig, OxygenConfig, PDConfig, SimConfig`, and `DRUG_LIBRARY` dict.
- Result containers in `core.results`: `TransportResult`, `SimulationResult`,
  `PatientRecord`, `ExperimentResult`. **Return these, don't reinvent.**

Units: length **um**, time **seconds**, concentration **uM**, diffusivity **um^2/s**.

## Module APIs (implement exactly these signatures)

### `tmesim/data/synthetic.py`
```python
def generate_synthetic_patient(patient_id: str, seed: int = 0,
                               grid=(120, 120), dx: float = 15.0,
                               subtype: str | None = None) -> PatientRecord:
    """Build a biologically-plausible synthetic GBM Microenvironment with a KNOWN
    ground-truth response, so the pipeline can be validated end to end.
    Must create: a tumor mass + infiltrative margin, stroma, immune cells, a
    vascular network (sparser in tumor core -> hypoxia), spatially varying ECM,
    and per-cell `sensitivity` correlated with metadata (e.g. mgmt_methylated,
    subtype). observed_response must be a *noisy function of true biology*
    (transport + intrinsic sensitivity), NOT trivially recoverable, so the test
    is fair. Populate `features` with non-spatial scalars for baselines."""

def generate_synthetic_cohort(n: int = 40, seed: int = 0,
                              drug_names: list | None = None) -> list[PatientRecord]:
    """A cohort spanning subtypes/vascularization/MGMT for the Stage 0 experiment."""
```

### `tmesim/data/acquire.py` and `tmesim/data/loaders.py`
```python
# acquire.py — best-effort LIVE download of OPEN glioblastoma data into ./data_cache/.
# Attempt, in order, and DEGRADE GRACEFULLY (log what worked/failed, never crash):
#   - Ivy GAP (Allen): open RNA-seq + anatomic structure data
#     https://glioblastoma.alleninstitute.org/  (gene_expression zip)
#   - TCGA-GBM open expression + clinical via GDC API https://api.gdc.cancer.gov/
#   - DepMap / GDSC glioma cell-line drug sensitivity (open CSVs) for response ground truth
def acquire_all(cache_dir: str = "data_cache", timeout: int = 60) -> dict:
    """Return {source_name: {'status': 'ok'|'failed'|'skipped', 'path': ..., 'note': ...}}."""

# loaders.py — parse whatever acquire.py cached into a list[PatientRecord].
def load_real_cohort(cache_dir: str = "data_cache") -> list[PatientRecord]:
    """Best-effort. If no usable real data, return [] (caller falls back to synthetic)."""
```
Also write `docs/data_acquisition_plan.md` covering the controlled-access datasets
(HTAN, spatial Visium/Xenium GBM on GEO, PDO drug screens) that a funded team
would use for real validation, with exact accession IDs / URLs where known.

### `tmesim/reconstruct/segmentation.py`
```python
def segment_and_phenotype(image: np.ndarray, dx: float = 15.0,
                          markers: dict | None = None) -> CellPopulation:
    """From a (multi-channel) tissue image, segment nuclei and assign CellType by
    marker intensity. Use skimage. Must run on synthetic images from
    data.synthetic (which can emit an image) and be robust to a plain grayscale."""
```

### `tmesim/reconstruct/spatial.py`, `vasculature.py`, `assemble.py`
```python
# spatial.py
def build_tissue_fields(domain: Domain, cells: CellPopulation,
                        vascular_density: np.ndarray | None = None) -> TissueFields:
    """Kernel-density the cells into ECM density (stroma-weighted) and, if not
    supplied, a vascular_density field. Returns TissueFields (oxygen/ifp left None)."""

# vasculature.py
def build_vascular_network(domain: Domain, seed: int = 0,
                           tumor_mask: np.ndarray | None = None) -> np.ndarray:
    """Return a (ny,nx) vascular_density field in [0,1]; realistic tortuous, and
    SPARSER inside the tumor core (drives hypoxia + poor delivery). networkx ok."""

# assemble.py
def assemble_microenvironment(domain, cells, fields, patient_id="reconstructed",
                              metadata=None) -> Microenvironment:
    """Bundle validated pieces into a Microenvironment (with shape checks)."""
```

### `tmesim/physics/diffusion.py`  (CORE NUMERICS — must be stable & tested)
```python
def simulate_drug_transport(micro: Microenvironment, drug: DrugParams,
                            source: np.ndarray, pde: PDEConfig,
                            velocity: np.ndarray | None = None) -> TransportResult:
    """Solve convection-diffusion-reaction for drug concentration C(x,y,t):
        dC/dt = div(D_eff grad C) - div(v C) + S - (clearance+binding) C
    - D_eff reduced by ecm_density (tortuosity): D_eff = D*(1 - alpha*ecm).
    - `source` is a (ny,nx) delivery rate from vasculature*bbb_penetrance*plasma(t)
      (the treatment module supplies the time profile; here accept a static source
      scaled internally by an optional plasma(t) — keep it simple & documented).
    - No-flux at domain/tissue-mask boundaries. Unconditionally stable preferred
      (implicit/ADI) or enforce a CFL-safe dt with clear error if violated.
    Accumulate and return TransportResult.auc = integral C dt."""
```

### `tmesim/physics/oxygen.py`, `physics/pressure.py`
```python
# oxygen.py
def solve_oxygen(micro: Microenvironment, cfg: OxygenConfig) -> np.ndarray:
    """Steady-state reaction-diffusion O2 field (ny,nx) in [0,1]; source from
    vascular_density, sink from local tumor-cell density. Return the field AND it
    is the caller's job to store into micro.fields.oxygen."""
def label_hypoxia(micro: Microenvironment, oxygen: np.ndarray, cfg: OxygenConfig) -> None:
    """Set CellState.HYPOXIC on cells sitting in oxygen<threshold (in place)."""

# pressure.py
def solve_ifp(micro: Microenvironment) -> tuple[np.ndarray, np.ndarray]:
    """Interstitial fluid pressure (mmHg, elevated in tumor core) via a Darcy/
    poroelastic-lite Poisson solve; return (ifp (ny,nx), velocity (2,ny,nx)=-K grad p).
    Elevated core IFP must REDUCE net inward delivery (the classic Jain effect)."""
```

### `tmesim/cells/agents.py`, `cells/pkpd.py`
```python
# agents.py
def update_cell_states(micro, oxygen, cfg: OxygenConfig) -> None:
    """Relabel VIABLE/HYPOXIC/QUIESCENT from local oxygen (in place, non-dead cells)."""

# pkpd.py
def apply_drug_kill(micro: Microenvironment, transport: TransportResult,
                    pd: PDConfig, drug: DrugParams) -> np.ndarray:
    """Map each cell's local drug AUC (sampled from transport.auc at the cell's
    grid index) to a kill PROBABILITY via a Hill function on exposure, modulated
    by per-cell `sensitivity`, hypoxia/quiescence resistance. Mark killed cells
    DEAD (in place). Return per-cell kill-probability array for diagnostics.
    Deterministic given a seed passed via numpy default_rng in the engine."""
```

### `tmesim/treatment/drugs.py`
```python
def get_drug(name: str) -> DrugParams:  # from DRUG_LIBRARY, KeyError-safe with message
def plasma_profile(pk: PKConfig, times: np.ndarray) -> np.ndarray:
    """Vascular (blood) concentration uM at each time — multi-dose exponential decay."""
def delivery_source(micro: Microenvironment, drug: DrugParams,
                    plasma_conc: float) -> np.ndarray:
    """(ny,nx) instantaneous delivery-rate source =
       vascular_density * bbb_penetrance * (plasma_conc - tissue) proxy; return the
       spatial source scaled to plasma_conc=1, engine multiplies by plasma(t)."""
```

### `tmesim/sim/engine.py`  (INTEGRATOR — wave 2)
```python
def simulate_response(micro: Microenvironment, cfg: SimConfig,
                      rng=None) -> SimulationResult:
    """Full pipeline for ONE patient+drug:
      1. oxygen field + hypoxia labeling
      2. IFP + velocity
      3. build delivery source (treatment.delivery_source) with plasma profile
      4. drug transport PDE (physics.diffusion) -> TransportResult
      5. pharmacodynamic kill (cells.pkpd) over the treatment course
      6. compute predicted_response in [0,1] from kill_fraction (documented map)
    Populate SimulationResult.fields with oxygen/ifp/auc/vascular for plotting.
    Must be deterministic given cfg.seed."""
```

### `tmesim/validation/baselines.py`, `metrics.py`, `experiment.py`  (wave 2)
```python
# baselines.py  — the cheap models the sim MUST beat
def baseline_mean(records: list[PatientRecord]) -> np.ndarray  # constant prior
def baseline_intrinsic_sensitivity(records) -> np.ndarray      # mean per-cell sensitivity only (no transport)
def baseline_expression_logreg(train, test) -> np.ndarray      # logistic reg on features dict
# each returns predicted responses aligned to the given records.

# metrics.py
def roc_auc(y_true, y_score) -> float
def brier(y_true, y_score) -> float
def calibration_curve(y_true, y_score, bins=10) -> tuple
def spearman(y_true, y_score) -> float   # for continuous outcomes

# experiment.py
def run_stage0_experiment(records: list[PatientRecord], cfg: SimConfig,
                          out_dir: str = "results") -> ExperimentResult:
    """For each patient: run simulate_response -> predicted_response. Score the sim
    vs every baseline with cross-validation where a baseline needs training.
    Fill ExperimentResult, set beat_baseline / verdict HONESTLY (sim must exceed the
    BEST baseline on primary_metric=roc_auc by a margin to 'pass'). Save a json."""
```

### `tmesim/viz/plots.py`, `viz/report.py`
```python
def plot_microenvironment(micro, out_path) -> None      # cells + vascular + ecm panels
def plot_simulation(result: SimulationResult, out_path) -> None  # oxygen/ifp/auc/kill + viable curve
def plot_validation(exp: ExperimentResult, out_path) -> None     # ROC + calibration + bar of metrics
# report.py
def build_report(exp, sim_example, micro_example, out_dir="results") -> str:
    """Write results/report.md embedding the figures + the honest verdict. Return path."""
```

### `tmesim/cli.py`  (wave 2)
```python
# argparse CLI: `python -m tmesim.cli demo`   (one synthetic patient, plots)
#               `python -m tmesim.cli run --n 40 --data synthetic|real|auto`
#               `python -m tmesim.cli acquire` (attempt real data)
```

## Rules for every agent
- Only create/modify the files assigned to you. Never touch `tmesim/core/*`.
- Code to the signatures above even if a dependency isn't built yet.
- Every public function gets a docstring and type hints. Add a `if __name__ ==
  '__main__':` self-test where feasible.
- Prefer vectorized numpy/scipy. Keep runtime for one synthetic patient (120x120
  grid) under a few seconds.
- Determinism: take a seed / np.random.default_rng, never global randomness.
- No network except in `data/acquire.py`.
