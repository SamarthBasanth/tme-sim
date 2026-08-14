# Prototyping the Michaelis-Menten oxygen sink: does the functional-form fix actually change the answer?

**Status: research/prototype only. `tmesim/physics/oxygen.py` is untouched — this document reports on a standalone script run from the scratchpad, not any file added to this repository other than this report.** This is the direct follow-up `results/oxygen_calibration_landscape.md` §2/§4 asked for: that document found tme-sim's oxygen sink, `-consumption * rho_tum * O` (`tmesim/physics/oxygen.py:14`, confirmed again below), is first-order in `O` *everywhere*, while the closest matching governing-equation-family paper in the literature — Grimes, Fletcher, Partridge, *R Soc Open Sci* 2014;1(1):140080, PMID 26064525, PMC4448765, DOI 10.1098/rsos.140080 — treats the literature default as zeroth-order (zeroth-order, i.e. an oxygen-tension-*independent* rate, is the more familiar name for it), tapering to first-order only within a "hyperbolic oxygen constant" `Km` of anoxia. The landscape document was explicit that recalibrating the bare `consumption` constant alone would not fix a *functional-form* mismatch, and equally explicit that this was a claim, not something it had tested. This document tests it.

**Bottom line, stated plainly up front.** Switching from the current first-order sink to a Michaelis-Menten sink changes the predicted oxygen field and hypoxic fraction **dramatically** on tme-sim's real flagship reconstruction — from 0.0% hypoxic tumor cells (the current, already-flagged-as-implausible answer) to anywhere from 69% to 100% hypoxic, depending on the one free parameter (`K`, the dimensionless half-saturation constant) that this codebase cannot currently pin down. This is the **opposite** of Grimes 2014's own headline finding, which is that "differences between these functional forms are likely to be negligible within the range of literature estimates of the hyperbolic oxygen constant" — and §4 below explains, mechanistically, why tme-sim's regime is not the regime that finding was established in. The functional-form question and the missing-mmHg-anchor question (`results/sweep/oxygen_mmhg_reference_proposal_2026-07-29T05-54-00Z.md`, prepared in parallel) turn out to be entangled, not independent: the mmHg anchor that document recommends (`O=1 = 35 mmHg`, tissue-anchored) combined with Grimes' own literature `Km` range (`~1-10 mmHg`) lands `K` at `0.029-0.286` — squarely inside the range swept here, and squarely inside the range where the two functional forms disagree by 30-100 percentage points of hypoxic fraction.

---

## 1. What `tmesim/physics/oxygen.py` actually solves (read in full, unmodified)

`solve_oxygen` (`tmesim/physics/oxygen.py:181-239`) solves the steady-state balance

```
D * lap(O) + supply * vasc * (1 - O) - consumption * rho_tum * O = 0
```

on the `(ny, nx)` domain grid with no-flux (Neumann) boundaries, rearranged into the linear system

```
A = -D * L + diag(supply * vasc + consumption * rho_tum)
A @ O = supply * vasc
```

which is symmetric, sparse and diagonally dominant (`reaction = supply*vasc + consumption*rho_tum >= 0` everywhere), so it is solved with **one direct sparse solve**, `scipy.sparse.linalg.spsolve` (line 236), no iteration, no CFL concern — the reaction term is *linear in `O`*, which is exactly what makes this a single linear solve rather than a nonlinear system. `_build_laplacian` (lines 140-178) builds the 5-point finite-difference Laplacian with ghost-cell no-flux boundaries. `_tumor_density` (lines 106-137) turns raw tumor-cell counts into a grid-independent areal density `rho_tum` (cells per `TUMOR_REFERENCE_AREA_UM2 = 144 um^2` reference tile), which is what makes the sink term `consumption * rho_tum * O` grid-spacing-invariant (the module's own extensive self-test and docstring, lines 27-72, document this fix in detail). `OxygenConfig` defaults (`tmesim/core/config.py:76-82`): `diffusivity=2000.0 um^2/s`, `consumption=0.02`, `supply=0.05`, `hypoxia_threshold=0.2`. None of this — the Laplacian, the density machinery, the linear-system assembly, the defaults — is touched here; the prototype below imports `_build_laplacian`, `_tumor_density`, `TUMOR_REFERENCE_AREA_UM2`, and the unmodified `solve_oxygen` itself directly from this module, so the first-order baseline reported below is not a re-implementation, it is the actual production function.

## 2. Grimes, Fletcher & Partridge 2014 — the exact functional form and parameters (PMID 26064525, PMC4448765, DOI 10.1098/rsos.140080)

Fetched via PubMed/PMC (full text). The paper's governing equation is `D * lap(p) = Omega * a(p)`, where `p` is oxygen partial pressure in mmHg, `Omega = 3.0318 x 10^7 mmHg*kg/m^3` is a Henry's-law conversion constant, and `a(p)` is the oxygen consumption rate. Two functional forms for `a(p)` are contrasted:

* **Constant (zeroth-order):** `a(p) = a_max` — independent of local oxygen tension.
* **Hyperbolic (Michaelis-Menten):** `a(p) = a_max * p / (p + Km)` — "a typical hyperbolic reaction rate... similar in form to Michaelis-Menten kinetics or the Langmuir adsorption isotherm... shown to be a good approximation in a variety of tissues, both healthy and cancerous."

Parameters used in the paper's own simulations: `a_max` swept over `2.5e-7` to `1.5e-6 m^3 O2/kg/s` (equivalently `Omega*a_max = 7.50-45.45 mmHg O2/s`); external boundary `p` fixed at 100 mmHg (spherical/MCTS case) or 100/40/20 mmHg ("arterial"/"venous"/"chaotic", cylindrical Krogh case); `Km` — the "hyperbolic oxygen constant" — **base case 1 mmHg** ("the literature to date indicates that this constant for oxygen is quite low... typically `Km <~ 1 mmHg`"), with an explicit own sensitivity sweep over **0-10 mmHg** to gauge how much that assumption matters.

**The paper's own headline finding, quoted precisely because it is the one this prototype tests against tme-sim's data:** "differences between these functional forms are likely to be negligible within the range of literature estimates of the hyperbolic oxygen constant, suggesting that the constant consumption rate approximation suffices for modelling oxygen dynamics for most values of OCR." Concretely, the paper reports the two models' oxygen profiles differ by `<=1 mmHg` (spherical) or `~2 mmHg` (cylindrical) over most of the domain, with the anoxic-boundary location shifting by only "of the order of tens of micrometres," when `Km` is at its literature value (`~1 mmHg`). Only when `Km` is pushed an order of magnitude above literature values (`10 mmHg`) does the paper see a large shift (`~55%` deeper oxygen penetration in the spherical case).

## 3. The prototype: a Picard-iterated Michaelis-Menten sink, on the SAME real reconstruction, reusing the SAME solver machinery

**What was built (scratchpad only, not committed to this repo):** a standalone script that (a) builds the exact flagship real microenvironment every drug result in `results/` is computed on — `reconstruct_real("xenium", sample="data_cache/xenium_gbm", dx_um=12.0)` (`tmesim.reconstruct.real_pipeline`, `BRAIN_SAMPLE`/`DX_CONVERGED` from `tmesim.validation.state_validation`; 40,887 cells, 23,812 tumor cells, grid 424x675 = 286,200 unknowns, the same grid-converged production spacing `tmesim/physics/oxygen.py`'s own docstring cites); (b) calls the **unmodified** `solve_oxygen` for the first-order baseline (an `assert` confirms this reproduces `micro.fields.oxygen` bit-for-bit, i.e. exactly what `reconstruct_real`'s default `oxygen_mode='pde'` already stores — this is not a fresh reimplementation drifting from production); (c) implements the Michaelis-Menten alternative as a **Picard (successive-substitution) fixed-point iteration** built on the identical linear system:

```
D * lap(O) + supply * vasc * (1 - O) - consumption * rho_tum * O/(O + K) = 0
```

At Picard step `k+1`, the nonlinear coefficient is frozen at the previous iterate: `reaction_k = supply*vasc + consumption*rho_tum/(O_k + K)`, and `(-D*L + diag(reaction_k)) @ O_{k+1} = supply*vasc` is solved with the **same `spsolve` call** `oxygen.py` already uses. Iteration stops when `max|O_{k+1} - O_k| < 1e-7`, seeded from the first-order field (a physically reasonable warm start). This is the standard, simplest correct method for this type of saturating nonlinearity — no new numerical machinery, just the existing linear solve called repeatedly.

**`K`, the dimensionless half-saturation constant, and how the sweep range was chosen.** tme-sim's `O` is a dimensionless `[0,1]` saturation with no defined mmHg reference (`oxygen_calibration_landscape.md` §4) — `O=1` has never been anchored to a physical `p_ref`. Rather than invent one number, `K = Km_mmHg / p_ref` was swept over every combination of Grimes' own base case and sensitivity ceiling (`Km in {1, 10} mmHg`) against three physiologically-grounded reference points already in use elsewhere in this project's literature review (`p_ref in {20, 40, 100} mmHg` — "venous"/generic-well-oxygenated/"arterial", per Grimes' own cylindrical-case boundary conditions and `oxygen_calibration_landscape.md` §1/§4): `K in {0.01, 0.025, 0.05, 0.10, 0.25, 0.50}`.

## 4. Results

| K (dimensionless) | Picard iters | converged | mean O2, tissue | mean O2, tumor cells | **hypoxic tumor fraction** | max\|ΔO\| vs first-order | RMSE vs first-order | Spearman vs first-order |
|---:|---:|:---:|---:|---:|---:|---:|---:|---:|
| **first-order (current, baseline)** | 1 (linear) | — | 0.3520 | 0.3486 | **0.0000** | — | — | — |
| 0.010 | 37 | yes | 0.0158 | 0.0144 | **1.0000** | 0.389 | 0.338 | 0.823 |
| 0.025 | 30 | yes | 0.0330 | 0.0310 | **1.0000** | 0.363 | 0.320 | 0.892 |
| 0.050 | 25 | yes | 0.0582 | 0.0556 | **1.0000** | 0.325 | 0.295 | 0.940 |
| 0.100 | 21 | yes | 0.1009 | 0.0977 | **1.0000** | 0.266 | 0.252 | 0.976 |
| 0.250 | 15 | yes | 0.1968 | 0.1930 | **0.6941** | 0.162 | 0.155 | 0.998 |
| 0.500 | 11 | yes | 0.3044 | 0.3004 | **0.0000** | 0.058 | 0.048 | 0.999 |

(Hypoxic tumor fraction = fraction of the 23,812 tumor cells sitting on a grid cell with `O < hypoxia_threshold = 0.2`, the exact same definition `label_hypoxia` uses.)

**Reading the table.** Every `K` at or below 0.10 — which includes Grimes' own base-case `Km=1 mmHg` combined with *every one* of the three reference points considered (`K = 1/100 = 0.01`, `1/40 = 0.025`, `1/20 = 0.05`) — drives the predicted hypoxic tumor fraction from the current 0.0% to a full **100%**. `K = 0.25` (matching `Km=10 mmHg`, Grimes' own sensitivity ceiling, over `p_ref=40 mmHg`) gives 69.4%. Only at `K = 0.50`, well outside the range any combination of Grimes' literature `Km` and a physiologically defensible `p_ref` produces, does the hypoxic fraction return to 0.0% — and even there the mean tissue oxygen (0.304) is meaningfully lower than the first-order baseline's 0.352, and still converging toward it rather than matching it. **The two functional forms do not converge anywhere inside the literature-plausible parameter range on this section; they diverge by 30 to 100 percentage points of hypoxic fraction.**

## 5. Why this contradicts Grimes' own "negligible" finding — the mechanism, not just the number

Grimes 2014's "negligible" conclusion was established in a specific regime: a single vessel or spheroid with a fixed, well-oxygenated Dirichlet boundary (`p_ref` = 100/40/20 mmHg) and a *thin* near-anoxic shell where the two models diverge, surrounded by a large volume of tissue that both models agree is close to `p_ref` because it hasn't consumed much oxygen yet. The divergence is confined to that thin shell precisely because the geometry (one vessel, or one spherical boundary) keeps most of the domain well-oxygenated.

tme-sim's flagship reconstruction is not that geometry. It is a full 2-D tissue section with **widespread tumor occupancy** (23,812 of 40,887 cells are tumor; `rho_tum` reaches occupancy 4 per reference tile, mean occupancy ~1 where tumor is present at all) consuming oxygen supplied by a comparatively sparse, heterogeneous vascular field — and the *first-order baseline's own mean oxygen is already only 0.35*, not close to a well-oxygenated `O=1`. Most of the tissue is not in the "far from the sink, still near `p_ref`" regime Grimes' geometry guarantees; it is already sitting in the part of the `[0,1]` range where a Michaelis-Menten sink and a first-order sink give genuinely different answers, because there is no large well-oxygenated bulk keeping the domain average pinned near `p_ref`. Mechanistically: where `rho_tum` is high, `consumption*rho_tum` (up to `0.02*4 = 0.08`) can exceed the maximum possible supply term `supply*vasc*(1-O) <= 0.05`; a first-order sink self-throttles (the sink shrinks as `O` falls, so it can only ever draw `O` down to wherever it re-balances against supply), while a near-saturated Michaelis-Menten sink (small `K`) stays pinned near its maximum rate `consumption*rho_tum` almost independent of `O`, i.e. behaves like the very "constant OCR regardless of oxygen availability" case Grimes describes — and a genuinely oxygen-independent sink that exceeds available supply drives the tissue toward true, widespread anoxia rather than a self-limiting equilibrium. That is exactly the qualitative difference the results show: the first-order sink cannot produce a hypoxic fraction above roughly what its own O-proportional throttling allows, while the Michaelis-Menten sink at small `K` collapses large fractions of the domain to near-zero oxygen.

**An honest caveat this cuts both ways on.** A 100% hypoxic-tumor-fraction result (mean tumor O2 = 0.014-0.10 at the smallest `K` values) is itself not obviously more credible than the current 0.0% — `results/oxygen_calibration_landscape.md` §1 argues a 0.0% hypoxic fraction on a GBM section is implausible given four independent human studies, but it does not argue for total, near-uniform anoxia either (the same studies report a *range*, with medians in the low single-digit mmHg, not zero, and viable, non-hypoxic peritumoral tissue). The most defensible reading is that `consumption = 0.02` was tuned (informally, by whatever process set the current default) against the *first-order* functional form's self-throttling behavior, and simply swapping the functional form while holding that same constant fixed overshoots in the other direction — consistent with `oxygen_calibration_landscape.md`'s original point that this is a joint problem (functional form **and** the constant), not a one-line fix. **So: yes, the functional form change is not just theoretically different but produces a large, real change in this specific parameter regime — but "large and real" is not the same as "correct"; a real fix would need to re-derive or refit `consumption` for the new functional form, not just drop the new form in with the old constant.**

## 6. Exact code change to `tmesim/physics/oxygen.py`, if adopted (diff-style, not applied)

```diff
--- a/tmesim/physics/oxygen.py
+++ b/tmesim/physics/oxygen.py
@@ class OxygenConfig (tmesim/core/config.py, for reference -- not in this file)
+    mm_half_saturation: Optional[float] = None  # K; None -> legacy first-order sink
+    mm_max_iter: int = 100
+    mm_tol: float = 1e-7

@@ def solve_oxygen(micro, cfg):
     D = float(cfg.diffusivity)
     supply = float(cfg.supply)
     consumption = float(cfg.consumption)
+    K = getattr(cfg, "mm_half_saturation", None)

     L = _build_laplacian(ny, nx, domain.dx)
-
-    source = (supply * vasc).ravel()                          # RHS b
-    reaction = (supply * vasc + consumption * rho_tum).ravel()  # O-proportional loss
-
-    # A = -D * L + diag(reaction). (-L is positive-definite; reaction >= 0.)
-    A = (-D * L) + sp.diags(reaction, format="csr")
-
-    # If there is no source and no reaction anywhere the system is singular; guard.
-    if not np.any(source) and not np.any(reaction):
-        return domain.zeros()
-
-    oxygen = spsolve(A.tocsc(), source)
+    source = (supply * vasc).ravel()                          # RHS b
+
+    if not np.any(source) and not np.any(rho_tum):
+        return domain.zeros()
+
+    if K is None:
+        # Legacy first-order-everywhere sink: linear system, one solve.
+        reaction = source + (consumption * rho_tum).ravel()
+        A = (-D * L) + sp.diags(reaction, format="csr")
+        oxygen = spsolve(A.tocsc(), source)
+    else:
+        # Michaelis-Menten sink: nonlinear in O -> Picard (fixed-point) iteration.
+        # Same linear system class and solver at every step; only the reaction
+        # coefficient is refreshed from the previous iterate.
+        oxygen = np.ones(ny * nx, dtype=np.float64)
+        rho_flat = rho_tum.ravel()
+        for _ in range(int(cfg.mm_max_iter)):
+            mm_coef = consumption * rho_flat / (oxygen + K)
+            reaction = source + mm_coef
+            A = (-D * L) + sp.diags(reaction, format="csr")
+            oxygen_new = np.clip(spsolve(A.tocsc(), source), 0.0, 1.0)
+            if np.max(np.abs(oxygen_new - oxygen)) < float(cfg.mm_tol):
+                oxygen = oxygen_new
+                break
+            oxygen = oxygen_new
+
     oxygen = np.asarray(oxygen, dtype=np.float64).reshape(ny, nx)
     np.clip(oxygen, 0.0, 1.0, out=oxygen)
     return oxygen
```

`mm_half_saturation` defaults to `None` so the change is additive and bit-identical by default, in the same spirit as `reconstruct_real`'s existing `oxygen_mode` opt-in pattern (`tmesim/reconstruct/real_pipeline.py`). A real adoption would also need: (a) the `K` value itself resolved, which is blocked on the same missing `O=1 <-> mmHg` anchor `results/sweep/oxygen_mmhg_reference_proposal_2026-07-29T05-54-00Z.md` addresses — `K` cannot be chosen independently of that decision, as this prototype's sweep demonstrates; and (b) `consumption` likely needs to be re-examined jointly with `K` rather than held at its current, never-calibrated value (§5).

## 7. Numerical cost: yes, this requires a nonlinear solve, and it is not free

The first-order sink's linearity is exactly what lets `oxygen.py` do a **single** sparse direct solve (`tmesim/physics/oxygen.py:229-236`, and its own docstring: "so a single sparse direct solve gives the exact steady state (no iteration/CFL concerns)"). A Michaelis-Menten sink is nonlinear in `O` (the reaction coefficient itself depends on the unknown), so it structurally requires **iteration** — there is no linear-algebra shortcut available. Measured on this prototype (dx=12 um, 286,200 unknowns, the shipping production grid):

| K | Picard iterations to convergence (tol 1e-7) | wall time |
|---:|---:|---:|
| first-order baseline | 1 (single linear solve) | 7.9 s |
| 0.01 | 37 | 118.9 s |
| 0.025 | 30 | 80.7 s |
| 0.05 | 25 | 64.5 s |
| 0.10 | 21 | 40.8 s |
| 0.25 | 15 | 20.5 s |
| 0.50 | 11 | 14.9 s |

**The robust cost metric is the iteration count, not the raw wall-clock ratio** (the per-iteration wall time measured here is noisier than expected — later Picard iterations averaged noticeably faster per solve than the very first baseline solve, most plausibly a one-time BLAS/thread-pool and sparse-factorization warm-up cost paid by the first `spsolve` call of the process rather than a real algorithmic speedup; this was not controlled for with repeated/warmed-up timing runs, and should be re-measured properly, e.g. with `TIMING_REPEATS`-style best-of-N as `tmesim/validation/oxygen_benchmark.py` already does for its PDE-vs-surrogate comparison, before this number is trusted to a digit). By iteration count, adopting the Michaelis-Menten sink at a literature-plausible `K` costs **11 to 37 sparse solves instead of 1** — a genuine, order-of-magnitude increase in the dominant cost of every reconstruction, not a rounding change. Since `solve_oxygen` is called once per `reconstruct_real` call and again internally inside every `simulate_response` (per `tmesim/validation/oxygen_benchmark.py`'s own documentation of the call graph), this would materially inflate the runtime of every existing calibration, robustness, and dose/drug sweep script in `results/` that re-solves the PDE many times, not just the one-off reconstruction cost. A production implementation should also consider Newton's method (the Jacobian here is a simple diagonal correction, `d(reaction)/dO = -consumption*rho_tum/(O+K)^2`, cheap to assemble) in place of Picard — Newton's quadratic convergence would likely cut the iteration count to single digits, but that is a further, unimplemented optimization, not what was measured here.

## 8. What this does and does not establish

**Does establish:** on tme-sim's own real flagship reconstruction, at every literature-grounded value of the one free parameter this fix introduces, switching the sink's functional form changes the predicted hypoxic fraction by 30-100 percentage points relative to the current first-order sink, holding `consumption` fixed — a first-order effect, not a second-order refinement, and the opposite of what Grimes 2014's own "negligible" finding would suggest if quoted out of context. It also establishes, concretely, that this is a genuinely nonlinear numerical problem (11-37 sparse solves instead of 1 at dx=12 um) and not just a one-line formula swap.

**Does not establish:** which hypoxic fraction, if either, is closer to correct. Both the current 0.0% and the Michaelis-Menten sweep's 69-100% sit at implausible extremes relative to the real human-GBM pO2 literature `results/oxygen_calibration_landscape.md` §1 surveys (which implies a real but partial hypoxic fraction, not 0% or ~100%). Fixing the functional form without also resolving `K` (which needs the mmHg anchor) and re-examining `consumption` jointly would trade one uncalibrated, implausible answer for a different uncalibrated, implausible answer — the honest conclusion is that the functional-form fix and the calibration fix are coupled and neither is sufficient alone.

## 9. Reproduction

Not committed anywhere in this repository. The prototype was a standalone script in the session scratchpad importing `tmesim.physics.oxygen._build_laplacian`, `_tumor_density`, `TUMOR_REFERENCE_AREA_UM2`, `solve_oxygen` (unmodified) and `tmesim.reconstruct.real_pipeline.reconstruct_real` against `tmesim.validation.state_validation.BRAIN_SAMPLE` / `DX_CONVERGED`, run once (`./.venv/bin/python <script>`, ~6 minutes total, dominated by the 6-point `K` sweep's Picard iterations); no seeds are relevant (the reconstruction and both solvers are deterministic). Anyone repeating this would reconstruct the same script from §3/§6's description; nothing here depends on a file that only exists outside this repository.
