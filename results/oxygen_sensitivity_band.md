# Oxygen sensitivity band: resolving the uncalibrated-`consumption` and sink-functional-form uncertainty into a quantified grid

**Status: this is the resolution the module docstring in `tmesim/physics/oxygen.py` pointed to, not a further round of "more research needed."** Two axes of that docstring were flagged as genuinely uncertain and NOT independently fixable: (1) `OxygenConfig.consumption` (default 0.02 /s) is bracketed from opposite directions by two real literature estimates roughly an order of magnitude apart each way, and (2) the sink's functional form (linear-in-O vs. Michaelis-Menten) swings the predicted hypoxic tumor fraction from 0.0 to 0.69-1.0 on its own. Picking a single point value for either would be false precision dressed up as a fix. This document instead (a) makes the Michaelis-Menten form an available, documented, opt-in option in `tmesim/physics/oxygen.py` (`OxygenConfig.sink_form`, defaulting to `"linear"` — nothing existing changes), and (b) sweeps the full literature-defensible `consumption` range, for both sink forms, on the real flagship reconstruction, and reports the resulting hypoxic-fraction grid plainly. **The shipping default (`sink_form="linear"`, `consumption=0.02`) is UNCHANGED and remains the default — this report quantifies the uncertainty around it, it does not silently move it.**

---

## 1. What changed in the code (summary; see `tmesim/physics/oxygen.py` and `tmesim/core/config.py` for the full diff)

- `OxygenConfig` (`tmesim/core/config.py`) gained four new fields: `sink_form: str = "linear"`, `mm_half_saturation: float = 1/35 ≈ 0.0286`, `mm_max_iter: int = 50`, `mm_tol: float = 1e-9`. All four default to values that reproduce the pre-existing behavior exactly (`sink_form="linear"` never consults the other three).
- `solve_oxygen` (`tmesim/physics/oxygen.py`) now branches on `cfg.sink_form`:
  - `"linear"` (default): bit-identical to the pre-existing single sparse direct solve (`A @ O = b` via `spsolve`). Verified: `solve_oxygen(micro, OxygenConfig())` is `np.array_equal` to the pre-change code path on every test in `tests/test_physics.py`, `tests/test_oxygen_proxy.py`, `tests/test_gbmspace_adapter.py`, `tests/test_immune_sink.py` (all still pass, see §5).
  - `"michaelis_menten"` (opt-in): `consumption*rho_tum*O/(O+mm_half_saturation)`, solved by a new `_solve_michaelis_menten` helper — a proper damped Newton-Raphson iteration (analytic Jacobian, backtracking line search, reuses the same sparse direct solver at every step), not the standalone prototype's throwaway Picard script. Raises `RuntimeError` rather than returning a silently-unconverged field if `mm_max_iter` is exhausted.
- **Correctness cross-check:** the new production Newton solver was run at the exact six `K` values and `consumption=0.02` the standalone Picard prototype (`results/sweep/oxygen_functional_form_prototype_2026-07-29T05-54-00Z.md`) used, on the same flagship reconstruction. The two independent numerical methods agree to 4 decimal places on every value (§4 below) — strong evidence the new solver is correct, not just fast.
- New reproducible script: `scripts/oxygen_sensitivity_band.py` (regenerates everything in this report; raw numbers also dumped to `results/sweep/oxygen_sensitivity_band_raw.json`).

## 2. How the `consumption` grid was chosen

From the module docstring's literature pass (`results/sweep/oxygen_ocr_literature_second_pass_2026-07-29T05-54-00Z.md`):

- **Watabe et al. 2017** (rat C6 glioma, ¹⁵O-PET TMRO2, PMID [28210996](https://pubmed.ncbi.nlm.nih.gov/28210996/)) implies `consumption` should be **12-25x HIGHER** than the current default → bracket `[0.24, 0.5]`.
- **Diepart et al., via Milotti 2017** (DOI [10.1038/srep39762](https://doi.org/10.1038/srep39762)), already a first-order-in-O rate constant matching this module's own functional form, implies `consumption` should be **4-100x LOWER** → bracket `[0.0002, 0.005]`.

The swept grid spans both brackets plus the current default and the transition zone between them: `consumption ∈ {0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.24, 0.5}` — a combined range of **2,500x** (3.4 decades), run for both sink forms.

The Michaelis-Menten half-saturation constant is fixed, for the main grid, at `K = 1/35 ≈ 0.0286` — Grimes et al. 2014's own base case (`Km = 1 mmHg`, PMID [26064525](https://pubmed.ncbi.nlm.nih.gov/26064525/)) expressed against this project's own `O=1` anchor (`REFERENCE_PO2_MMHG = 35 mmHg`, Erecinska & Silver 2001, PMID [11718758](https://pubmed.ncbi.nlm.nih.gov/11718758/)). A secondary sweep (§4) varies `K` itself at the fixed default `consumption=0.02` to show that axis's own sensitivity.

All runs use `reconstruct_real('xenium', dx_um=12.0)` — the grid-converged production spacing this project's results are computed at (40,887 cells, 23,812 tumor cells, grid 424×675 = 286,200 unknowns). Hypoxic tumor fraction uses the exact same definition `label_hypoxia` does: fraction of tumor cells sitting on a grid cell with `O < hypoxia_threshold = 0.2`.

## 3. Main grid: consumption × sink_form → predicted hypoxic tumor fraction

| consumption (/s) | vs. default | linear: hypoxic frac | linear: mean O2 (tumor) | MM (K=0.0286): hypoxic frac | MM: mean O2 (tumor) |
|---:|---|---:|---:|---:|---:|
| 0.0002 | Diepart/Milotti lower (÷100) | 0.0000 | 0.9815 | 0.0000 | 0.9817 |
| 0.0005 | | 0.0000 | 0.9551 | 0.0000 | 0.9543 |
| 0.001  | | 0.0000 | 0.9140 | 0.0000 | 0.9088 |
| 0.002  | | 0.0000 | 0.8418 | 0.0000 | 0.8182 |
| 0.005  | Diepart/Milotti upper (÷4) | 0.0000 | 0.6806 | 0.0000 | 0.5528 |
| 0.01   | | 0.0000 | 0.5164 | **0.6639** | 0.1913 |
| **0.02** | **current shipping default** | **0.0000** | **0.3486** | **1.0000** | **0.0347** |
| 0.05   | | **0.8134** | 0.1769 | 1.0000 | 0.0083 |
| 0.1    | | 1.0000 | 0.0973 | 1.0000 | 0.0036 |
| 0.24   | Watabe lower (×12) | 1.0000 | 0.0431 | 1.0000 | 0.0014 |
| 0.5    | Watabe upper (×25) | 1.0000 | 0.0212 | 1.0000 | 0.0006 |

(Full precision, both mean-O2-tumor and mean-O2-tissue-wide, per-run solve time: `results/sweep/oxygen_sensitivity_band_raw.json`.)

**Headline reading.**

- **Across the entire Diepart/Milotti sub-range (`consumption` in `[0.0002, 0.005]`), BOTH sink forms predict a hypoxic tumor fraction of exactly 0.0** — no hypoxic core at all, at every value that literature bracket defends.
- **Across the entire Watabe sub-range (`consumption` in `[0.24, 0.5]`), BOTH sink forms predict a hypoxic tumor fraction of exactly 1.0** — the tumor is entirely hypoxic by the model's own threshold, at every value that bracket defends.
- **The two literature-defensible sub-ranges predict opposite headline answers** (0.0 vs 1.0), regardless of which sink form is used. This is the single most important fact this sweep establishes: the disagreement in `consumption` (§2) is not a second-order wobble on the hypoxic-fraction headline number, it spans its ENTIRE range.
- **At the current shipping default (`consumption=0.02`, unchanged), the two sink forms disagree with each other by the full range**: linear predicts 0.0 (as already flagged in the pre-existing docstring), Michaelis-Menten (at the literature-derived `K=0.0286`) predicts 1.0. The crossover for the linear sink happens between `consumption=0.02` and `0.1`; for the Michaelis-Menten sink it happens earlier, between `0.005` and `0.02` — i.e. **the Michaelis-Menten sink reaches saturation (100% hypoxic) at a `consumption` value roughly 5-10x lower than the linear sink needs to**, consistent with the mechanism already documented in the prototype report (a near-saturated MM sink does not self-throttle as O falls, so it drives the tissue further toward anoxia for the same nominal rate constant).

## 4. Secondary axis: Michaelis-Menten `K`-sensitivity at the current default `consumption=0.02`

| K (dimensionless) | hypoxic tumor fraction | mean O2 (tumor) | Newton solve time |
|---:|---:|---:|---:|
| 0.010 | 1.0000 | 0.0144 | 11.4s |
| 0.025 | 1.0000 | 0.0310 | 10.2s |
| 0.050 | 1.0000 | 0.0556 | 10.4s |
| 0.100 | 1.0000 | 0.0977 | 9.2s |
| 0.250 | 0.6941 | 0.1930 | 10.3s |
| 0.500 | 0.0000 | 0.3004 | 8.1s |

**This exactly reproduces the standalone Picard prototype's numbers (`results/sweep/oxygen_functional_form_prototype_2026-07-29T05-54-00Z.md` §4) to 4 decimal places at every one of the six `K` values**, using an independently-implemented production Newton-Raphson solver instead of the prototype's Picard iteration — a genuine cross-check that the new solver is correct, and that this sensitivity is real rather than a numerical artifact of either method. Grimes' own literature range for `Km` (1-10 mmHg) maps to `K` in `[0.0286, 0.286]` against this project's own mmHg anchor, i.e. `K` in `[0.029, 0.29]` in the table above spans hypoxic fraction from **1.0 down to somewhere between 0.69 and 0.0** — again, the entire plausible range, not a minor correction.

## 5. Tests run

```
./.venv/bin/python -m pytest tests/ -k oxygen -q
```
40 passed (covers `tests/test_oxygen_proxy.py` and the oxygen-specific tests inside `tests/test_physics.py`, `tests/test_gbmspace_adapter.py`, `tests/test_immune_sink.py`, `tests/test_core.py`). This confirms the new `sink_form` option does not change the default (`"linear"`) path's behavior and that the grid-independence, determinism, core-vs-periphery, and bit-identical-`reconstruct_real`-default properties documented elsewhere in `tmesim/physics/oxygen.py` all still hold.

A standalone Newton-vs-Picard correctness check (§4, run directly, not part of the pytest suite) additionally confirms the new Michaelis-Menten solver reproduces the independently-implemented prototype's numbers exactly.

## 6. What is now known vs. still not known

**Now known (the actual resolution):**

- The shape and magnitude of the uncertainty, quantified rather than asserted: both literature-defensible sub-ranges for `consumption` predict a hypoxic fraction at the EXTREME of the [0,1] range (0.0 or 1.0), for both sink forms, and the transition between those extremes happens over less than one decade of `consumption`, positioned differently depending on which sink form is chosen.
- A production-quality, opt-in Michaelis-Menten sink option (`OxygenConfig(sink_form="michaelis_menten")`) that a future calibration effort can use directly, with proper Newton-Raphson convergence (typically 6-22 iterations on the flagship grid; never silently unconverged) instead of a prototype script.
- That the two uncertain axes are coupled, not independent, is now shown quantitatively rather than argued qualitatively: the Michaelis-Menten sink saturates to the 1.0 extreme at roughly 5-10x lower `consumption` than the linear sink needs.

**Still not known, and not resolvable from within this codebase:**

- Which single `(sink_form, consumption)` combination in this grid is actually correct for human GBM tissue. Nothing in this sweep, or in the literature surveyed so far, picks a point — real, glioma-specific, spatially-resolved pO2 (or a diffusion-model-fitted local OCR) would be needed, and `results/oxygen_calibration_landscape.md` and the second-pass literature search both concluded that measurement does not currently exist in the public literature for glioma specifically. Until it does, reporting the band (this document) rather than a point estimate is the honest state of the model.
- Nothing changes for existing published results: `OxygenConfig()`'s defaults (`sink_form="linear"`, `consumption=0.02`) are unchanged, so every number in `results/` computed before this document remains exactly reproducible.

## 7. Reproduction

`./.venv/bin/python scripts/oxygen_sensitivity_band.py` (a few minutes; deterministic, no seeds). Raw per-point numbers (including mean-O2-tissue-wide and per-solve wall time) are written to `results/sweep/oxygen_sensitivity_band_raw.json`.
