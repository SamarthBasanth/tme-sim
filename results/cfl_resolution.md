# The CFL discrepancy, resolved: two correct measurements of two different tissues

**Two records state the maximum interstitial speed of the "same" configuration and straddle the advection stability limit at `dx = 8 um`. Re-measured on the shipped configuration, one of them reproduces bit-for-bit and the other does not reproduce at all — and the reason is not the vascular field, the tissue mask, the amplitude or the link distance. It is that the cell-typing rule changed 1 h 38 min after the record that disagrees was generated.**

> **This is a reconciliation, not a calibration.** Nothing in `tmesim/` was modified. `DX_UM`, `dt`, `alpha` and every frozen snapshot are untouched. `results/grid_convergence_gbm.md` and `results/vessel_promotion_provenance.md` are not edited by this file; §6 states what each needs.

## 0. Bottom line

| question | answer |
|---|---|
| Are the GBM ladder's fine rungs (`dx = 8`, `dx = 7` at `dt = 10 s`) admissible? | **YES.** Re-measured CFL **0.939** and **0.955** on the shipped configuration — bit-identical to what that ladder recorded. |
| Which measurement describes the shipped configuration? | **`results/grid_convergence_gbm.md`.** Its speeds reproduce to the last bit of a float64 (§2). |
| Then what is `results/vessel_promotion_provenance.md` measuring? | **The reconstruction as it stood at commit `ceb96c0`**, before the cell-classifier arity fix. Its four discrete speeds AND the four KDE speeds in `results/vessel_promotion.md` §1 all reproduce exactly when the legacy typing rule is restored (§3). |
| Does the GBM ladder's conclusion stand? | **YES — it stands, unqualified.** No re-run was needed, but `dx = 7` was re-run at `dt = 5 s` anyway: the scored ratio moves by at most **0.00007 %** and every sign of every refinement step is unchanged (§4). The non-monotonicity is a property of the mesh, not of the timestep. |
| Is anything stale? | **YES.** `tmesim/reconstruct/real_fields.py` lines 146-149 assert, in the present tense and about the shipped field, that `dt = 10 s` violates the CFL at `dx <= 8`. On the shipped field it does not. That claim is now a fact about a superseded reconstruction (§6). |

**Neither record contains an error.** Both measured correctly; they measured different tissues, four weeks apart, and neither said which. That is the whole of the discrepancy.

## 1. The exact configuration re-measured

Every number in §2 comes from `reconstruct_real(source='xenium', dx_um=dx)` with **no keyword overridden**, followed by `validation.kill_calibration.prepare_tissue` — the drug-independent prologue (`solve_oxygen` -> `update_cell_states` -> `solve_ifp`) that the published pipeline itself runs. Stated in full, because the whole disagreement turns on which of these differed:

| knob | shipped value used here | where the default lives |
|---|---|---|
| vascular mode | **`discrete`** (segmented vessels, exactly zero away from a vessel) | `real_fields.VASCULAR_MODES` / `_DEFAULT_VASCULAR_MODE` |
| vascular amplitude | **`unit`** (1.0 inside a vessel, 0.0 outside) | `vessel_network.discrete_vascular_field` |
| vessel radius | **7.5 um** | `vessel_network.VESSEL_RADIUS_UM` |
| link distance | **20.0 um** | `vessel_network.LINK_DISTANCE_UM` |
| geometry / radius mode | **`binary` / `data`** | `vessel_network.discrete_vascular_field` |
| tissue mask | **`geometric`**, 30 um discs (grid-independent) | `real_pipeline.TISSUE_MASK_MODE` |
| oxygen / state route | **`pde` / `oxygen`** | `real_pipeline.OXYGEN_MODE_DEFAULT`, `STATE_MODE_DEFAULT` |
| **cell typing rule** | **`matched_null_z`** | `real_xenium.TYPING_AGGREGATION` |
| convection, seed | **on, 0** | `PDEConfig`, `SEED` |

Reconstructed population, identical at every `dx` (the mask is grid-independent): **40,887 cells, 23,812 tumour, 4,295 endothelial, 3,091 segmented vessel objects**. Tissue pixels 124,074 / 193,752 / 435,924 / 569,385 / 775,035 at `dx = 15 / 12 / 8 / 7 / 6`.

**Control on the arm construction.** Where a variant arm below rebuilds the vascular field through `vascular_density_from_cells`, the shipped arm rebuilds it the same way with the shipped defaults and asserts the result is `np.array_equal` to the field `reconstruct_real` built itself. That check **passed at all five `dx`**, so the variant arms differ from the shipped arm by their named keyword and by nothing else.

**The statistic.** `max(|vx| + |vy|)` over the whole grid — byte-for-byte the quantity `physics/diffusion.py` line 361 computes and compares against `dx/dt` before it raises. The ladder's script computes exactly this (`experiments/grid_convergence_gbm.py` line 600); the promotion record has no surviving script, but its published values are reproduced by this L1 statistic and not by the Euclidean one (§3.1), so it used the same form. The Euclidean maximum `max(sqrt(vx^2+vy^2))` is reported alongside it because it runs **24-29 % lower** on this section, and quoting that as "max speed" would make an inadmissible grid look safe.

## 2. Re-measured speeds and CFL on the shipped configuration

| `dx` (um) | grid | **max speed, L1 (um/s)** | Euclidean max (um/s) | `dt_max = dx/speed` (s) | **CFL at `dt = 10 s`** | **CFL at `dt = 5 s`** | solver at `dt = 10 s` |
|---:|---|---:|---:|---:|---:|---:|---|
| 15 | 339x541 | 0.7077 | 0.5019 | 21.20 | **0.472** | 0.236 | OK |
| 12 **(published)** | 424x675 | 0.6336 | 0.4798 | 18.94 | **0.528** | 0.264 | OK |
| 8 | 635x1012 | 0.7513 | 0.5321 | 10.65 | **0.939** | 0.470 | OK |
| 7 | 725x1157 | 0.6684 | 0.4900 | 10.47 | **0.955** | 0.477 | OK |
| 6 | 846x1349 | 0.6773 | 0.4803 | 8.86 | **1.129** | 0.564 | **REFUSES** |

**Against `results/grid_convergence_gbm.md` §2 — bit-identical.** Compared to the persisted probe that file was generated from (`results/sweep/gbmgrid_cfl_probe.json`), at every `dx` the two records share:

| `dx` | this re-measurement | the ladder's stored probe | equal as float64? |
|---:|---|---|---|
| 12 | 0.6336488624409844 | 0.6336488624409844 | **yes** |
| 10 | 0.6629421712103596 | 0.6629421712103596 | **yes** |
| 8 | 0.7513136752848846 | 0.7513136752848846 | **yes** |
| 7 | 0.668377690832811 | 0.668377690832811 | **yes** |
| 6 | 0.6772910174041108 | 0.6772910174041108 | **yes** |

Not "agrees to four figures" — the same float64, at every `dx` the two runs share. (`dx = 15` is not on that ladder, so there is nothing to compare it against; `dx = 10` was added here for §3.3.)

**Against `results/vessel_promotion_provenance.md` — does not reproduce.** That record states 0.730 / 0.779 / 0.843 / 0.868 at `dx = 15 / 12 / 8 / 6`; the shipped configuration gives 0.7077 / 0.6336 / 0.7513 / 0.6773. The record is high by **+3.2 % / +23.0 % / +12.2 % / +28.1 %** respectively, and the `dx = 8` excess is what flips that rung from CFL 0.939 to 1.054.

**So the fine rungs are admissible, with the margin stated plainly.** `dx = 8` runs at 94 % of the limit and `dx = 7` at 96 %. There is 4.5 % of headroom at the finest rung: a section with slightly faster interstitial flow would refuse it. The floor is a property of this section, not of the solver.

## 3. The cause, found by measurement rather than by assumption

Five candidate causes were run as explicit arms at the same five `dx`. Four are eliminated and the fifth reproduces the disputed record exactly.

| arm | `dx=15` | `dx=12` | `dx=8` | `dx=7` | `dx=6` | verdict |
|---|---:|---:|---:|---:|---:|---|
| **shipped** | 0.7077 | 0.6336 | 0.7513 | 0.6684 | 0.6773 | the reference |
| KDE vascular field | 0.2739 | 0.2429 | 0.2842 | 0.2718 | 0.2698 | **not it** |
| legacy `occupancy` tissue mask | 0.8498 | 0.8867 | 0.7241 | 0.6500 | 0.3871 | **not it** |
| `amplitude='per_vessel_cells'` | 0.5032 | 0.5863 | 0.6370 | 0.5192 | 0.5593 | **not it** |
| `link_distance_um=40` | 0.7077 | 0.6336 | 0.7513 | 0.6684 | 0.6773 | **not it** |
| **legacy `mean_z` cell typing** | **0.7301** | **0.7794** | **0.8429** | 0.8402 | **0.8678** | **this is it** |

**The KDE explanation, checked first because it was the cheapest, is wrong.** Today's KDE field runs at 0.2429-0.2842 um/s — *below* both disputed measurements, not between them. It cannot be a mislabelled arm of the discrete record, and it is not what the ladder measured either (the ladder reproduces bit-identically as the discrete field, §2).

**The tissue mask is not it, and the disputed record itself proves it.** `results/vessel_promotion_provenance.md` states tissue-pixel counts of 124,074 / 193,752 / 435,924 / 775,035 at `dx = 15 / 12 / 8 / 6`, and the shipped `geometric` mask reproduces **all four exactly** today (§1). The mask has not moved between the two runs, so it cannot be what separates them. (The `occupancy` arm is a genuinely different domain and its speeds match neither record. It is also a clean demonstration of why the mask was switched: `geometric` holds the simulated tissue at **27.90 mm² at every `dx`** — 27.917 / 27.900 / 27.897 / 27.899 / 27.900 / 27.901 mm² at `dx = 15 / 12 / 10 / 8 / 7 / 6` — while `occupancy` shrinks it from 8.58 mm² to 1.47 mm² over the same refinement.)

**The link distance cannot be it, structurally — and the arm shows why.** Doubling `link_distance_um` from 20 to 40 um genuinely changes the segmentation, merging **3,091 vessel objects into 2,060**, and the resulting speed is nonetheless **bit-identical to `shipped` at all five `dx`**. With `amplitude='unit'` the painted field is the union of 7.5 um discs around every endothelial cell that belongs to some vessel object, and single-linkage clustering labels every endothelial cell (singletons included); the linking rule re-partitions those cells into objects but cannot change their union. Link distance is a no-op for this field and therefore cannot move the CFL number at all.

### 3.1 The classifier arity fix, and the timeline that settles it

| commit | date | what |
|---|---|---|
| `ceb96c0` | **2026-07-28 01:03:46** | vessel promotion — the run `results/vessel_promotion_provenance.md` records |
| `4ad221f` | **2026-07-28 02:42:15** | *"Classifier arity fix: endothelial spread 2.91x -> 1.49x"* — `TYPING_AGGREGATION` changed from `mean_z` to `matched_null_z` |
| `88c7882` / `823fcbe` | 2026-08-09 | the GBM grid-convergence ladder |

**One hour and thirty-eight minutes** separate the promotion run from the change that invalidates its speed table. Restoring the legacy rule through the public `reconstruct_from_xenium(typing_aggregation='mean_z')` keyword — no source modified — reproduces both rows of `results/vessel_promotion.md` §1 exactly:

| | `dx = 15` | `dx = 12` | `dx = 8` | `dx = 6` |
|---|---:|---:|---:|---:|
| discrete, **recorded** in `vessel_promotion.md` | 0.730 | 0.779 | 0.843 | 0.868 |
| discrete, **re-measured under `mean_z`** | **0.7301** | **0.7794** | **0.8429** | **0.8678** |
| KDE, **recorded** in `vessel_promotion.md` | 0.361 | 0.407 | 0.509 | 0.551 |
| KDE, **re-measured under `mean_z`** | **0.3610** | **0.4073** | **0.5093** | **0.5512** |
| recorded `dt_max`, discrete (s) | 20.5 | 15.4 | 9.5 | 6.9 |
| re-measured `dt_max` under `mean_z` (s) | **20.54** | **15.40** | **9.49** | **6.91** |

**All eight cells of that table reproduce to every digit published** — both the discrete arm and the KDE arm, at all four `dx`. Under the shipped typing, none of the eight reproduces. **Measurement A is correct — for the tissue it was computed on.**

### 3.2 What actually changed, and the honest limit of the attribution

The typing rule decides which measured cells are endothelium and which are tumour. Under `mean_z` this section reconstructs as **3,984 endothelial and 30,174 tumour** cells; under `matched_null_z` as **4,295 endothelial and 23,812 tumour**. Those two numbers are already printed on the face of the two documents — `vessel_promotion_provenance.md` says 30,174 tumour and 3,984 endothelial, `grid_convergence_gbm.md` says 23,812 tumour — so the discrepancy was legible in the record before anything was re-run.

Both feed `solve_ifp`: the endothelial cells set the vascular source, and the tumour cells set `_tumor_field`, which raises vessel leakiness and suppresses lymphatic drainage inside the tumour. A 26.7 % larger tumour field is a materially different pressure problem.

**This arm does not separate the two.** It moves the endothelial and tumour populations together, because that is what the classifier change did. The measurement establishes *that the typing rule is the cause* and reproduces the disputed numbers from it; it does not apportion the effect between the vascular source and the drainage field. That decomposition was not needed to answer the question and was not run.

### 3.3 The other stale claim this touches

`dea609e` corrected `results/delivery_metrics.md` and the `DX_UM` comment in `kill_calibration.py`, which had said *"`dx = 10 um` violates the advection CFL condition at `dt = 10 s`"*. That correction's **direction and floor are confirmed here** — on the shipped configuration `dx = 7` is the finest admissible grid at `dt = 10 s` and `dx = 6` refuses. Its **attributed cause is not**. It guesses the claim "most likely predates the tissue-mask switch from `occupancy` to `geometric`"; `dx = 10` was therefore run under all three candidate reconstructions:

| arm at `dx = 10`, `dt = 10 s` | max speed (um/s) | CFL | solver |
|---|---:|---:|---|
| shipped | 0.6629 | 0.663 | OK |
| legacy `occupancy` tissue mask | 0.7004 | 0.700 | OK |
| legacy `mean_z` cell typing | 0.8343 | 0.834 | OK |

**None of them makes `dx = 10` infeasible**, so neither the mask nor the classifier explains the original claim. Note the third row: even the reconstruction that *does* reproduce the disputed 0.843 at `dx = 8` runs `dx = 10` at CFL 0.834. The origin of the `dx = 10` claim remains unexplained; what is now measured is that no arm tested reproduces it.

## 4. The finest rung, re-run at half the timestep

§2 says no re-run was *required*: `dx = 8` and `dx = 7` are admissible at `dt = 10 s`, so neither was measuring instability. But `dx = 7` is the rung the ladder never gave a `dt` control of its own — §7.3 halved `dt` at `dx = 12` and `dx = 8` and stopped there — and it is the rung with the least CFL headroom on the whole ladder (0.955). So it was re-run at `dt = 5 s`, which halves the Courant number to 0.477, through `grid_convergence_gbm.score()` — the published scoring path, imported unmodified, so the numbers below are directly comparable with that file's §3.

**3 transport solves, 4,320 steps each, `dx = 7 um`, everything else at its published value** (`alpha = 148.6081576545`, 6 h window, `save_every = 10`, seed 0, convection on, p99 bath matching, same reporting concentrations):

| drug | ratio at `dt = 10 s` (CFL 0.955) | ratio at `dt = 5 s` (CFL 0.477) | delta | delta % | `F_pen` | `F_state` | exposure-field change | verdict |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| temozolomide | 0.25850806 | 0.25850808 | +1.92e-08 | **+0.00001 %** | 0.470219 -> 0.470219 | 0.540673 | 1.18e-04 | PASS -> PASS |
| cisplatin | 0.27260976 | 0.27260988 | +1.20e-07 | **+0.00004 %** | 0.496181 -> 0.496181 | 0.540673 | 1.16e-04 | PASS -> PASS |
| doxorubicin | 0.20561309 | 0.20561323 | +1.36e-07 | **+0.00007 %** | 0.372324 -> 0.372324 | 0.540673 | 1.16e-04 | PASS -> PASS |

**`dt` reached the solver**: step count doubles (2,160 -> 4,320), the CFL number halves, the solves take 421 / 691 / 561 s, and the per-cell exposure vector changes (relative change in its sum, 1.16-1.18e-04). The field moved; the score did not.

**Worst movement from halving `dt` at the finest rung: 0.00007 %.** That is the same order as the 0.0003 % §7.3 found at `dx = 12` and `dx = 8`, on a rung with half the CFL headroom — so the ladder's *finest* rung behaves exactly like its middle ones under a timestep change. It is not a special case, and it is not near an edge in the way the disputed record implied.

### 4.1 The ladder's shape is unchanged

Substituting the `dt = 5 s` value for `dx = 7` into the published ladder and re-reading the sequence fine -> coarse:

| drug | published sequence of step directions | with `dx = 7` taken at `dt = 5 s` | fine-rung span (`dx <= 18`) |
|---|---|---|---|
| temozolomide | down down up down up down down | **identical** | 6.60 % -> **6.60 %** |
| cisplatin | down down up down up down down | **identical** | 7.38 % -> **7.38 %** |
| doxorubicin | down down up down up down down | **identical** | 6.74 % -> **6.74 %** |

**Every sign of every refinement step is unchanged, and the fine-rung spans are unchanged to three significant figures.** The non-monotonicity is a property of the mesh, not of the timestep.

### 4.2 The verdict on the GBM ladder, plainly

**The conclusion of `results/grid_convergence_gbm.md` stands. It is not an artefact and it is not provisional.**

Specifically: the panel is not in an asymptotic convergence regime; the curve is non-monotone in `dx` for all three drugs; the power-law fits are degenerate and their `q0` values are not converged limits; the ratio oscillates 6.6-7.4 % below `dx = 18`; and the 3/3 PASS survives because the literature bands are wide, not because the numerics are clean. **None of that changes.** The two fine rungs it rests on were admissible when it ran, they are admissible now, and the finest of them moves by 7e-5 % when its Courant number is halved.

The provenance record's closing instruction — *"until then the two finest rungs of that ladder should be read as provisional"* — was the right call on the evidence available when it was written. That provisionality is now discharged: **the rungs are sound, and the three documents that cite the ladder need no correction on this account.**

## 5. Reconciling the dt controls — why they were so quiet

`results/grid_convergence_gbm.md` §7.3 halved `dt` at `dx = 12` and `dx = 8` and found the scored ratio moved by at most **0.0003 %**. The brief asks how that can be, if the fine rungs were unstable.

**It is not a puzzle once §2 lands: the fine rungs were not unstable, and a quiet dt control is exactly what a stable rung produces.** Three things line up and none of them needed explaining away:

1. **The solver is a hard gate, not a warning.** `physics/diffusion.py` line 363 raises on `CFL > 1`. A rung that returned a number had `CFL <= 1` in the run that produced it — that was never in doubt. The only live question was whether the ladder ran the shipped configuration, and §2 answers it to the bit.
2. **The `dx = 8` control halves a CFL of 0.939 to 0.470** — a factor-2 change in the advection Courant number, on a rung that the disputed record would have made inadmissible. The exposure field moved 1.2e-4 relative; the score moved 1e-7. A field that was ringing at CFL 0.939 could not have done that.
3. **The scheme is IMEX, not fully explicit.** The upwind flux divergence enters the right-hand side (`diffusion.py` line 452) of an implicit diffusion-plus-reaction solve (line 453), so the implicit operator damps the explicit term rather than amplifying it. `CFL <= 1` is the standalone-upwind bound and is therefore conservative here.

**The one real mechanism in the neighbourhood, and its measured size.** Upwind advection carries numerical diffusion proportional to `(1 - CFL)`, so holding `dt` fixed while refining `dx` *does* vary the artificial diffusion across the ladder — from 0.53 of the limit at `dx = 12` to 0.955 at `dx = 7`. That is a legitimate candidate mechanism for a non-monotone curve and it is not a stability question at all. §7.3's control bounds it at `dx = 12` and `dx = 8`; §4 extends the bound to `dx = 7`, the rung where it would be largest, and finds 0.00007 %. **Measured, the mechanism is real but far too small to explain a 6.6-7.4 % oscillation** — three to four orders of magnitude too small.

**One caveat worth recording.** The solver clips negative undershoot every step (`np.maximum(c_new, 0.0, out=c_new)`, line 455). A nonlinear clip can mask mild oscillation rather than let it grow, so "the run completed" is weaker evidence of cleanliness than it looks. The `dt` controls, which compare a near-limit CFL against a comfortable one on the same mesh, are the evidence that carries weight here — not the absence of a crash.

## 6. What the two source documents need

Neither is wrong and neither needs a retraction on this point. **This file does not edit them.**

* **`results/vessel_promotion_provenance.md`** — its §"A CFL discrepancy this record exposes" poses the question this file answers; the answer is §3. Its speed table and its cell counts are correct **for commit `ceb96c0`** and do not describe the shipped reconstruction. It already says "anything below that matters should be re-measured rather than cited", which was the right instruction and is now carried out.
* **`results/vessel_promotion.md`** §1 — the same applies to its speed table, which reproduces exactly under `mean_z` and not at all under the shipped typing. **The wider issue is not the CFL:** every number in that file was computed on a section with 30,174 tumour and 3,984 endothelial cells, and today's reconstruction has 23,812 and 4,295. The promotion decision may well survive re-measurement — nothing here tests it — but its numbers are not reproducible on `HEAD`, and that is worth knowing independently of this file's question.
* **`results/grid_convergence_gbm.md`** — needs **nothing changed**. Its §2 CFL table reproduces bit-identically, its fine rungs are admissible, and §4 tests the one rung it never gave a `dt` control of its own. If anything, its §2.1 could gain a pointer here, since that section attributes the stale `dx = 10` claim to the tissue mask and §3.3 measures that the mask does not do it.
* **`tmesim/reconstruct/real_fields.py` lines 146-149** — this is the item that should move, and it was deliberately left alone. It reads: *"the max interstitial speed is ~1.6-2x the KDE's and dt_max falls to 20.5 / 15.4 / 9.5 / 6.9 s ... The historical dt = 10 s VIOLATES the CFL for this field at dx <= 8; use dt <= 5 s."* Measured on the field that docstring describes: `dt_max` is 21.20 / 18.94 / 10.65 / 8.86 s and `dt = 10 s` is admissible down to `dx = 7`. The guidance is conservative rather than dangerous, but it is stated as a present-tense property of the shipped field and it is a property of a superseded one. Its owner should replace the numbers or date them. (Its **qualitative** claim survives and its quantitative one does not: the discrete field really is faster than the KDE, but by **2.5-2.6x** on the shipped reconstruction against the **1.6-2.0x** the docstring quotes — which is itself exactly the ratio the legacy arm reproduces, 1.57-2.02x.)

## 7. Limitations of this file

1. **One section, one patient.** Every speed here is the 10x Xenium GBM section. The 4.5 % headroom at `dx = 7` is that section's, not the solver's.
2. **The classifier arm moves two populations at once** (§3.2) and does not decompose the effect between the vascular source and the tumour drainage field.
3. **`max(|vx| + |vy|)` is a stability statistic, not an accuracy statistic.** Nothing here says the velocity field is right — only that the solver's own admissibility test is satisfied and that two records of it now have a single explanation.
4. **This file does not re-audit the promotion decision**, only the speed table it published. Whether discrete-over-KDE survives the current typing is a separate re-run that nobody has done.
5. **The re-run in §4 covers `dx = 7` only.** `dx = 8` already had a `dt` control in the ladder's own §7.3 and was not repeated. `dt = 5 s` was chosen because it halves the CFL and is the value the promotion work used; no third `dt` was run, so this bounds the temporal error rather than extrapolating it.

---

*Generated by `experiments/cfl_resolution.py` — **38 prologue probes** across 7 configuration arms and 6 grid spacings, plus **3 transport solves** at `dx = 7 um`, `dt = 5 s` (421 / 691 / 561 s), seed 0, deterministic. Raw records: `results/sweep/cfl_speed_probe.jsonl` and `.json` (one line per arm x dx), `results/sweep/cfl_rerun_records.jsonl` (the `dt`-halved fine rung). `tmesim/` unmodified; `results/grid_convergence_gbm.md`, `results/vessel_promotion_provenance.md` and `results/vessel_promotion.md` unmodified.*
