# Stale blockers — tests declined on capabilities that now exist

**Date: 2026-09-20. Status: SWEEP. Nothing here is a new scientific result.** No
anchor was registered, no verdict moved, no existing file was edited. Every
"this now exists" line below was verified by **importing or executing the code**
in this pass, and the output is pasted. Where a blocker still holds, it is
recorded as still holding — a sweep that recovers everything is as useless as one
that recovers nothing.

**Calibration guard, start of pass:**

```
$ .venv/bin/python -c "from tmesim.validation.calibration import calibrate_alpha; print(repr(calibrate_alpha()))"
START calibrate_alpha() = 22.223477711172563
ASSERT OK
```

*(The closing assertion is at the end of this file, §6.)*

---

## 0. What this sweep is, and why it is a different defect from the one already tracked

Claim-surface drift — the pattern this repository has been tracking since
`results/claim_surface_audit_2026-08-30.md`, and through
`results/verdict_inversion_census.md`, `results/grid_dependence_audit.md` and
`results/withdrawal_banner_log_2026-09-03.md` — is **the document overstating
what the code does**. Every guard, every audit and every retraction so far has
been pointed in that direction.

**This sweep looks the other way.** It asks: where does a `results/` document
*decline a whole class of test* on the grounds that the engine cannot do
something — and the engine can do it today? That failure costs nothing visible.
Nothing goes red. No claim inflates. What is lost is **verification that was
available and never taken**, which is exactly the "overly cynical" error
`results/verification_routes_2026-09-20.md` was commissioned to avoid and which
that document itself names in its §1.3.

`results/verification_routes_2026-09-20.md` §1.3 found three. This sweep grepped
all 149 files in `results/` for the shapes a declined test takes — `cannot`,
`not possible`, `no facility`, `would require`, `out of scope because`,
`blocked on`, `not testable`, `NOT TESTABLE`, `deferred`, `not actionable`,
`has no`, `not implemented`, `structurally cannot`, `cannot be expressed` — and
triaged every hit against the current tree.

**Counts: 27 numbered rows, 28 entries (row 1 splits). 9 NOW POSSIBLE.
7 PARTIALLY. 10 STILL BLOCKED. 2 already closed by the documents themselves.**

### 0.1 The finding that surprised me most

**The median gap between a blocker being written and the capability landing is
ONE DAY.** In three cases the capability **already existed when the blocker was
written**. This is not a story about documents rotting over months. It is a
story about a document and a module being produced in the same burst of work by
different threads, with nothing reconciling them — the exact failure mode
`results/` already knows as territory drift, arriving on the understatement side
where no guard is watching.

And the pattern is live, not historical: **`results/verification_routes_2026-09-20.md`,
written today, carries a stale blocker of its own** (row 9 below).

---

## 1. The table

Dates are the date the file was added to the repository (`git log --diff-filter=A`),
which for these files is the date they were written. "Capability" is the file and
line that removes the blocker, with the date it landed. Every NOW POSSIBLE and
PARTIALLY row was checked by running the code, not by reading a name.

| # | What was declined | Where, and when | Blocker as written | Capability, verified | Gap | Verdict |
|---|---|---|---|---|---|---|
| 1a | Tunggal 1999 / Cowan 2001 multicellular-layer flux tests, on four of the engine's own panel drugs | `results/external_anchor_survey.md:132-133,307-310` (2026-08-09); `results/anchor_registration.md:490` (2026-08-10); `results/axis_inventory.md:274` (2026-08-21) | *"the engine … has **no 1-D slab geometry**. Computing this would require a new solve configuration"*; *"NOT ACTIONABLE without new code"* | `experiments/isomer_penetration_map.py:775 _slab` + `:795 mcc_ratio` (2026-08-22) — a 140 µm slab solved through the frozen `simulate_response` | 13 / 12 / 1 d | **NOW POSSIBLE** |
| 1b | …specifically the *donor/receiver steady-state flux ratio* ("% of the rate through the blank membrane") | same | same | `bath_diffusion` takes a **scalar** `c_bath` (`moncayo_explant.py:528`); a two-valued Dirichlet raises | — | **PARTIALLY** |
| 2 | A true fixed-concentration face for Kyle 2004 arm A; the construction bias was estimated analytically instead of solved | `results/kyle2004_penetration.md:432` (2026-08-22) | *"Kyle's near face is a stirred reservoir at fixed concentration (a Dirichlet condition); **the engine has no Dirichlet inflow**"* | `tmesim/validation/moncayo_explant.py:528 bath_diffusion` (2026-07-27), Crank-verified at `:608`; Kyle's one-sided geometry is reachable **exactly** by mirror symmetry | capability **preceded** the doc by 26 d | **NOW POSSIBLE** |
| 3 | Comparing the solved oxygen field to measured human-glioma pO₂ in mmHg | `results/oxygen_calibration_landscape.md` §4 (2026-07-28) | *"there is no `p_ref` … anywhere in `OxygenConfig`. **Before any mmHg number … can be plugged in, someone has to decide this mapping**"* | `tmesim/physics/oxygen.py:215 REFERENCE_PO2_MMHG = 35.0`, `:225 to_mmhg` (2026-07-29) | **1 day** | **NOW POSSIBLE** |
| 4 | Testing the zeroth-order / Michaelis–Menten oxygen sink against the shipped first-order one | `results/oxygen_calibration_landscape.md` §4 item 2 (2026-07-28) | *"moving to a Michaelis–Menten or zeroth-order sink would make the system nonlinear and require iteration, **which is a real implementation cost**"* | `tmesim/core/config.py:170 sink_form`, `:179 mm_half_saturation`; `tmesim/physics/oxygen.py:222 SINK_FORMS = ("linear", "michaelis_menten")` with a Newton–Raphson solve (2026-07-29) | **1 day** | **NOW POSSIBLE** |
| 5 | Any test of regional barrier integrity (contrast-enhancing vs non-enhancing) | `results/glioma_regional_pk.md:116` (2026-07-27) | *"**The model has NO region-specific BBB mechanism** … There is no term in the molecular layer, the transport PDE, or `DrugParams` that represents regional barrier integrity … **it cannot even express the answer**"* | `tmesim/reconstruct/barrier.py:117 L_PARACELLULAR_LEAK = 0.238`, `:179 effective_penetrance`, `:264 bbb_disruption_field`, `:358 barrier_delivery_source` (2026-07-27); wired at `tmesim/treatment/drugs.py:260` | **same day** | **NOW POSSIBLE** — but see §3.1, the obvious test is spent |
| 6 | Any temozolomide-sensitivity term other than uniform 1.0 | `results/regimen.md:428` (2026-07-28) | *"**This model has no MGMT term at all** (MGMT is not on the Xenium panel; sensitivity is uniform 1.0)"* | `tmesim/molecular/mgmt.py:250 mgmt_sensitivity`, `:403 is_mgmt_calibrated` (2026-07-27); wired at `tmesim/reconstruct/real_pipeline.py:59` | capability **preceded** the doc by 1 d | **PARTIALLY** — scalar yes, depletion/regeneration dynamic no |
| 7 | Sequential-vs-concurrent scheduling, and any fractionation contrast | `results/combination.md:212` (2026-07-27) | *"schedule (sequential vs concurrent …) **cannot be represented at all**"* | `tmesim/core/config.py:87 time_resolved`; `PKConfig(n_doses, dose_interval)`; `tmesim/protocol/dosing.py` (2026-07-28) expands `dose / every / repeat / rest / cycles` and infusions | **1 day** | **NOW POSSIBLE** on the transport axis |
| 8 | Multi-cycle / combination pharmacodynamics | `results/kill_calibration_breast.md:336` (2026-07-29) | *"this pipeline has no combination or multi-cycle pharmacodynamic model"* | PK/transport half exists (row 7). PD map is still a static Hill on window-averaged exposure — `results/regimen.md` §7 says so itself | — | **PARTIALLY** |
| 9 | Ingesting the 30 GBM-Space Visium sections | `results/gbmspace_ingestion.md:146,235` (2026-07-28) — **and restated as live in `results/verification_routes_2026-09-20.md` R12, written today** | *"Reading it correctly would require either the … a processed AnnData format that **requires new reader code**"*; R12: *"ingestion currently blocked on an `.h5ad` reader … **Add ~3 [days] for an `.h5ad` reader**, which unlocks 30 more"* | `tmesim/data/adapters/gbmspace.py` (2026-07-28) reads `.h5ad` with `h5py` alone. **All 30 sections are already materialised on disk** and `ingest.detect()` returns `visium / spaceranger_outs` | same day; **54 days** for the R12 restatement | **NOW POSSIBLE** |
| 10 | Replicating the n = 1 claims on a second section | `results/CLAIM_INVENTORY.md` §6 item 1 (2026-08-12); `results/verification_routes_2026-09-20.md` R12 (2026-09-20) | *"A second section, same tumour type. Almost every defensible claim above is n = 1 in tissue. This is the cheapest, largest improvement available"*; R12: *"currently costing four sections of replication for free"* | `results/multisection_generalisation.json` (written 2026-08-23, `alpha = 22.223477711172563`) holds **seven completed sections**, `status = ok`, including all four Visium GBM. **No `results/*.md` references it** | 11 d / **28 d** | **NOW POSSIBLE — and largely already run** |
| 11 | A genuinely 3-D solve | `results/dimensionality.md:315,326` (2026-07-28) | *"which would require a 3-D reconstruction we do not have"*; *"a genuinely 3-D solve would require serial sections, **which we do not have**"* | `tmesim/reconstruct/registration_3d.py` (2026-07-28); `tmesim/validation/dimensionality.py:530 solve_helmholtz_dct` — *"Works in any number of dimensions"* — plus `:580 build_network` / `:643 control_transverse_is_exact` | same day | **PARTIALLY** — `core.Domain` is still strictly 2-D |
| 12 | Ingesting the Greenwald CODEX layer for serial-section registration | `results/landscape_gbm_serial_sections.md:179` (2026-07-28) | *"The CODEX layer is **not** ingestible without new code: it is OME-TIFF multiplex imagery plus a QuPath project"* | `tmesim/data/adapters/codex.py` (2026-07-28): `read_ome_metadata`, `extract_nucleus_positions`, `load_section`. The file is on disk at `data_cache/greenwald_gbm_codex/ZH1007_NEC_B.ome.tif` | same day | **NOW POSSIBLE — already done** (`results/gbm_3d_registration.md`) |
| 13 | Composing an efflux term into the transport chain | `results/two_route_transport.md:524` (2026-08-13) | *"**This model has no efflux term**, and `efflux.py` exists in the package but is not in this chain"* | `tmesim/treatment/drugs.py:252` imports `efflux_factor` under `efflux=True`. But `efflux_ratio_for('paclitaxel')` → `None` — no record for the very compound the passage is about | capability preceded by 17 d | **PARTIALLY** |
| 14 | `tail90` at any `dx ≠ 12 µm` | `results/tail90_falsification.md:62` (2026-08-10) | *"`tail90` therefore **cannot be computed** at `dx ≠ 12 µm` … without re-running the transport solves … the fix is one line in each ladder script — persist `e_unit_by_key`"* | The one-line fix **landed the same day** (`results/small_fixes.md` §3): `experiments/grid_convergence_gbm.py:210 E_UNIT_STORE`, `:421 store_e_unit`, `:645`. But `results/sweep/gbmgrid_e_unit.npz` **does not exist** — no ladder has been re-run | same day | **PARTIALLY** — code blocker gone, data blocker intact |
| 15 | Anything requiring the `PS_SCALE_S_INV` conversion | `tmesim/validation/scorecard.py:735` | *"THE MODEL CANNOT EXPRESS THE MEASURED QUANTITY … an implicit, underived conversion of exactly 1 1/s"* | `tmesim/validation/ps_scale.calibrate_ps_scale()` runs and the scorecard row **already calls it** (`scorecard.py:730-742`). Shipped default is still `drugs.py:134 PS_SCALE_S_INV = 1.0`; the two anchors disagree **8.70×** | — | **STILL BLOCKED** (correctly; the calibration is already surfaced *in* the row) |
| 16 | Ionisation / pH-partitioning in the transport path | `results/axis_inventory.md:294` (2026-08-21); `results/transport_anchor_prereg.md:564,733`; `results/brain_slice_anchor.md:277` | *"**The engine has no pH field at all**"*; *"the engine has no ionisation term in this path"* | `tmesim/molecular/ion_trapping.py` exists (2026-07-28) and was exercised by `results/ion_trapping_impact.md` — but `ion_trapping.py:60` and `immune_sink.py:40` both state *"no shipped code path imports this module"*, and the grep confirms it. There is no pH **field** anywhere in `core/`, `physics/` or `reconstruct/` | — | **PARTIALLY** — scalar mechanism exists off-path; the field does not |
| 17 | CosMx / Stereo-seq / Visium HD ingestion | `results/landscape_serial_sections.md:7,55,68` (2026-07-28) | *"they would need to be reshaped into the generic MEX-plus-coordinates form, or a new loader written"* | `ingest.detect` (`ingest.py:318`) still dispatches to exactly `_detect_xenium` / `_detect_visium` / `_detect_mex`, and its own error message still lists only those three | — | **STILL BLOCKED**, softened: `adapters/gbmspace.py` and `adapters/codex.py` establish the materialise-into-a-recognised-layout pattern, so the cost is an adapter, not a loader rewrite |
| 18 | A TAM drug **re-release** term | `results/immune_sink_impact.md:415` (2026-07-28) | *"A re-release term would partly or wholly cancel the sink and could flip its sign. **That term is not implemented**, and it is the single largest structural omission here"* | Confirmed absent from `tmesim/molecular/immune_sink.py` — the module carries the literature account of the depot in prose only | — | **STILL BLOCKED** |
| 19 | A per-compound tissue-capacity term | `results/two_route_promotion.md:453`; `results/two_route_efflux.md:372,840` (2026-08-13/14) | *"The model has no capacity term"* | Confirmed absent. `results/brain_slice_anchor.md` (2026-08-21) records the pre-registered FAIL against 33 measurements, `GMFE_floor = 4.77×` | — | **STILL BLOCKED** (already re-verified by the repo) |
| 20 | A measured vessel radius for the PS decomposition | `results/permeability_split.md:370` | *"**this repository does not carry a measured vessel radius**, which is why `surface_to_volume_per_um` has no default radius"* | Confirmed: no default radius | — | **STILL BLOCKED** (data) |
| 21 | A quantitative GBM ECM/collagen area-fraction anchor | `results/reconstruction_validation_roadmap.md:127` item 7 (2026-07-28) | *"blocked on the same KDE-normalisation problem **already solved once for vasculature**"* | The vascular half **was** solved (`results/area_conservation_promotion.md`, 2026-08-15: `discrete_vascular_field(geometry='area')` default from 2026-08-14, area conserved to 6.3e-14 %). ECM is still a KDE: `reconstruct/real_fields.py:245 ecm_density_from_cells` → `:72 _kde_on_grid` | — | **STILL BLOCKED** — but there is now a worked template |
| 22 | The within-section DCE-MRI regional experiment | `results/landscape_imaging_derived.md:287,358` (2026-07-28) | *"none of our three sections has any paired MRI … **cannot be run on the data we currently hold, at any voxel count**"* | No paired MRI has been acquired | — | **STILL BLOCKED** (data; requires prospective collection) |
| 23 | A spatial lysosomal / acidic-phospholipid field | `results/moncayo_explant.md:163` (2026-07-27) | *"**no measured mg/g or mol% acidic-phospholipid content exists for tumour tissue** and no published PBPK model carries tumour-specific lysosomal parameters"* | Unchanged | — | **STILL BLOCKED** (literature) |
| 24 | Scoring the Ivy GAP zone labelling | `results/niche_response.md:184` (2026-07-27) | *"A pathologist's zone map of this section does not exist in the open data, so **the labelling cannot be scored**"* | Unchanged | — | **STILL BLOCKED** (data) |
| 25 | Regrowth-from-sanctuary | `results/vascular_validation.md:232` (2026-07-27) | *"Regrowth-from-sanctuary … **cannot be represented at all**"* | Unchanged — no proliferation is simulated during the exposure window (`results/regimen.md` §7) | — | **STILL BLOCKED** (structural) |
| 26 | Pruijn 2008, "blocked on retrieval" | `results/anchor_registration.md:488` (2026-08-10) | *"blocked on retrieval and not on judgement"* | **Already closed by the file's own 2026-08-14 banner** | — | closed, recorded |
| 27 | Reconciling `vivum validate` with `kill_calibration_breast.md` over zero-width bands | `results/CLAIM_INVENTORY.md` §6 item 0 (2026-08-12) | *"the cheapest item on this list"* | **Already closed** (`results/zero_width_bands.md`), and already noted in `verification_routes_2026-09-20.md` §1.3 | — | closed, recorded |

**Class counts, one per entry, 28 entries over 27 numbered rows:**

| verdict | n | rows |
|---|---:|---|
| **NOW POSSIBLE** | **9** | 1a, 2, 3, 4, 5, 7, 9, 10, 12 |
| **PARTIALLY** | **7** | 1b, 6, 8, 11, 13, 14, 16 |
| **STILL BLOCKED** | **10** | 15, 17, 18, 19, 20, 21, 22, 23, 24, 25 |
| closed by the document itself | 2 | 26, 27 |

**Ten of twenty-six live rows are still genuinely blocked** — five of them on
data or literature that does not exist, which no amount of re-reading fixes.
That ratio is the answer to "is this just cynicism in reverse?" It is not: the
sweep recovers about a third of what it examined and leaves the rest where it
found it.

---

## 2. The verifications, pasted

Nothing above is asserted from a function name. This is the evidence.

### 2.1 Row 2 — the Dirichlet bath, and Kyle's geometry solved with one

```
$ .venv/bin/python -c "from tmesim.validation.moncayo_explant import verify_against_analytic_slab; verify_against_analytic_slab()"
  [slab] t=     900s  numeric=0.33807  analytic=0.33931  rel_err=0.37%
  [slab] t=    3600s  numeric=0.66440  analytic=0.66791  rel_err=0.53%
  [slab] t=   14400s  numeric=0.97258  analytic=0.97717  rel_err=0.47%
  [slab] t=   86400s  numeric=1.00000  analytic=1.00000  rel_err=0.00%
worst_rel_err = 0.005261854604708477   passed = True
```

**Kyle's one-sided geometry is reachable exactly.** A slab with a Dirichlet face
at `x = 0` and a sealed face at `x = T` is, by mirror symmetry, identical to a
two-sided bath slab of thickness `2T` — which is what `bath_diffusion` already
builds. Run on daunorubicin at the arm-A configuration (`ecm = 0.0875`,
`T = 140 µm`, `dx = 2 µm`, 30 µm windows, 1 h AUC, `first_order_loss=True`, the
same `k = clearance + binding` the engine uses at `physics/diffusion.py:315`):

```
daunorubicin: D=28.2728  clearance=7.430e-05  binding=6.384e-03
engine source-column near:far = 3.6457        <- reproduces results/accuracy_gap_decomposition.md:227 exactly
Dirichlet-face near:far       = 3.5902
symmetry check (mirror residual) = 1.41e-15   <- the reduction is exact, not approximate
```

So the construction term, **measured on a solved field at the matching 1 h AUC
readout, is ×0.9848**. `results/accuracy_gap_decomposition.md:229` records
×0.9890 for the same term — but that number came from the *analytic* steady-state
transfer function `R_of_lambda` (`experiments/accuracy_gap_decomposition.py:105`),
not from a Dirichlet solve. They are the same sign and the same order; they are
not the same object, and the solved one did not need to be estimated.

**The residual for row 1b.** A donor/receiver pair needs two different Dirichlet
values. It is refused today:

```
$ bath_diffusion(mask, 5.0, drug, c_bath=np.array([1.0, 0.0]))
TWO-VALUED DIRICHLET REFUSED: TypeError only length-1 arrays can be converted to Python scalars
```

`moncayo_explant.py:576` builds the boundary vector as
`L_TR @ np.full(len(r_idx), float(c_bath))` — a single scalar broadcast over the
whole ring. That is a genuine, if small, obstacle and it is why row 1b is
PARTIALLY and not NOW POSSIBLE.

### 2.2 Row 1a — a 1-D slab solved through the frozen public solver

```
alpha 22.223477711172563
SLAB SOLVE THROUGH simulate_response: nx= 71  thickness 140.0  dx 2.0
doxorubicin near:far ratio = 2.9976871659491016
```

This is `simulate_response` — the shipped path — on `experiments/isomer_penetration_map.py:775 _slab`.
The geometry the 2026-08-09/08-10 declines said did not exist runs today, and its
value reproduces `results/accuracy_gap_decomposition.md:92`'s 2.9976 to four
decimals.

**One caution before anyone runs the Cowan arm.** Cowan 2001's
concentration-invariance from 10 µM to 3 mM tests *linearity of the sink*. The
engine's sink is `k = clearance + binding`, a constant, at
`physics/diffusion.py:315` — **linear by construction**. Running it would produce
a PASS that could not have failed: the same vacuous-green-check shape
`results/anchor_registration.md` already flags for the in-sample doxorubicin row
(`scorecard.py:511`, *"the model reproduces it by construction and could not have
failed"*). Register Cowan as a consistency check, never as a validation.

### 2.3 Rows 3 and 4 — mmHg, and both sink forms

```
REFERENCE_PO2_MMHG = 35.0
to_mmhg([0.0, 0.2, 0.4, 1.0]) = [ 0.  7. 14. 35.]
```

The `O < 0.2` hypoxia band is therefore `< 7 mmHg`, directly commensurable with
the Eppendorf literature. And both sink forms solve, on the same probe geometry,
through the same entry point:

```
SINK_FORMS = ('linear', 'michaelis_menten')
OxygenConfig fields: [diffusivity, consumption, supply, hypoxia_threshold,
                      sink_form, mm_half_saturation, mm_max_iter, mm_tol]

linear            median pO2=  0.453 mmHg   frac<2.5=1.0000   frac<10=1.0000   min=0.3386
michaelis_menten  median pO2=  0.012 mmHg   frac<2.5=1.0000   frac<10=1.0000   min=0.0011
```

> **Read that block as a capability demonstration and nothing else.** It was run
> on a synthetic 64×64 lattice with a regular vessel grid, not on any
> reconstructed section. **It is not a statement about GBM.** What it establishes
> is that the chain *solve → `to_mmhg` → sub-2.5 / sub-10 mmHg fraction* executes
> end to end under both sink forms, which is precisely what
> `oxygen_calibration_landscape.md` §4 said could not be done without a modelling
> decision and an implementation cost. Both were paid on 2026-07-29.

Note also `config.py:179 mm_half_saturation = 1.0/35.0` — the Michaelis constant
is already expressed against `REFERENCE_PO2_MMHG`, i.e. `Km = 1 mmHg`, which is
the value `oxygen_calibration_landscape.md` §2 itself quotes from the literature.

### 2.4 Row 5 — the regional barrier term

`glioma_regional_pk.md` §5 listed three things that "would have to be added".
Items 1 and 2 are in the tree, in the form it specified:

```
L_PARACELLULAR_LEAK = 0.238
effective_penetrance(0.065, 1.0) = 0.303     (cefazolin, non-enhancing 0.065 -> CE)
effective_penetrance(0.602, 1.0) = 0.840     (levetiracetam)
field form: effective_penetrance(0.065, [0.0, 0.5, 1.0]) = [0.065 0.184 0.303]
```

`results/revalidation.md:113` already refers to the model's "**pre-fix** regional
prediction", so the fix was known to at least one document. `glioma_regional_pk.md`
was never annotated.

### 2.5 Rows 6 and 13 — MGMT and efflux, both wired, both partial

```
temozolomide MGMT-calibrated: True
mgmt_sensitivity(status='methylated')   = 1.0
mgmt_sensitivity(status='unmethylated') = 0.39232336956521796

efflux_factor('temozolomide') = 0.6666666666666666
efflux_factor('paclitaxel')   = 1.0          <- efflux_ratio_for('paclitaxel') is None
efflux_factor('doxorubicin')  = 1.0          <- efflux_ratio_for('doxorubicin') is None
```

Both are PARTIALLY, for opposite reasons. MGMT: the **static scalar** exists and
is wired, but `regimen.md`'s actual argument is about *depletion and regeneration
with a time constant*, which is a state variable and is not implemented — so the
half of the claim that reads "sensitivity is uniform 1.0" is stale and the half
that reads "no schedule mechanism" is not. Efflux: the **term** exists and is
wired opt-in on the delivery leg, but it is not composed into the two-route
membrane leg, and there is **no measured efflux ratio for paclitaxel**, the one
compound `two_route_transport.md:524` was worried about. That second obstacle is
a data gap, not a code gap, and no amount of wiring removes it.

### 2.6 Row 7 — the schedule

```
simulate_drug_transport params: [micro, drug, source, pde, velocity, source_profile, capture_frames, exchange_k]
PKConfig fields: [plasma_peak, half_life, n_doses, dose_interval]

time_resolved=False   mean AUC=2452.1004   near:far=2.9744
time_resolved=True    mean AUC=2483.1282   near:far=2.9447
```

Two doses, 1 h apart, 1800 s half-life, on the arm-A slab. The two paths give
different fields, so the schedule is genuinely resolved rather than nominally
available. `tmesim/protocol/dosing.py` additionally expands
`dose / every / repeat / rest / cycles` and infusions into `PKConfig`, with the
four hidden conversion quantities written into provenance.

### 2.7 Row 9 — the 30 GBM-Space sections

```
materialised sections: 30      with .gbmspace_materialized marker: 30      raw .h5ad remaining: 2

detect('data_cache/gbmspace_visium/AT10-BRA-5-FO-1_1')
  -> {"platform": "visium", "layout": "spaceranger_outs"}

ingest(...) -> SampleRecord
  sample_id: AT10-BRA-5-FO-1_1 | platform: visium | n_cells(spots): 2956 | n_genes: 36601
  extent_um: (5373.4, 5300.6) | coord units: micrometer
```

The ingest capability report on that section, verbatim from the record:

| capability | status |
|---|---|
| `vascular_field` | **YES** — measurable from PECAM1/VWF/FLT1/PLVAP/CD34/CLDN5, endothelial detection 35.69 % |
| `cell_typing` | **YES** — marker map `brain_319` fits exactly (coverage 1.00) |
| `hypoxia_measured` | **YES** — detection 94.93 %, band size uncalibrated |
| `ecm_field` | **YES** — real expression-weighted KDE over 11 ECM genes |
| `mgmt_from_panel` | **YES** — *"MGMT IS on this panel, so temozolomide sensitivity can be read spatially"*, detection 16.75 % |
| `spatial_spread` | **NO** — per-cell statistics invalid on spot data, and the record says so at length |

That last row is the point `verification_routes` R12 makes about not discarding a
whole section because one statistic is refused — and the ingest layer **already
draws the line correctly, per capability**. The sections are not blocked. They
are loaded, typed, and waiting.

### 2.8 Row 10 — seven sections already solved, and never written up

```
results/multisection_generalisation.json   (mtime 2026-08-23, 129 KB)
  alpha: 22.223477711172563     drug: temozolomide

section                        status    n_cells    platform   spatial_spread_valid   penetration_length_um
xenium_319@dx30                ok         40,887    xenium     True                   277.69
xenium_io@dx30                 ok        816,769    xenium     True                   133.16
gse324268@dx30                 ok        851,509    xenium     True                   103.18
visium_10x_parent_gbm@dx30     ok         41,560    visium     False                   66.27
visium_gbm1@dx30               ok         57,050    visium     False                   71.87
visium_gbm2@dx30               ok         34,707    visium     False                  412.91
visium_gbm3@dx30               ok         21,279    visium     False                  301.20
xenium_319@dx12                ok         40,887    xenium     True                   278.49
```

Every run carries an `oxygen_block`, a `cell_state_block`, `celltype_fractions`,
`tissue_coverage_fraction`, `penetration_fit_error` and `kill_fraction`. The
`alpha` in the file is the pinned 22.223477711172563. `spatial_spread_valid` is
`False` on exactly the four spot sections and `True` on the three Xenium ones —
the refusal applied per-statistic, as it should be.

**No file in `results/`, `docs/`, `README.md`, `app/` or `landing-3d/` mentions
`multisection_generalisation`.** A grep for the string across all of them returns
nothing but this file. The run happened, the numbers are on disk, and the claim
surface does not know.

**And the report generator already exists.** `experiments/multisection_report.py`
opens with: *"Render `results/multisection_generalisation.md` from the sweep's
JSON … Every number in the document is interpolated from that file — nothing is
typed by hand, so the prose cannot drift away from the run that produced it."*
Its `_OUT` is `results/multisection_generalisation.md` (line 24). The data, the
solver run and the renderer are all present; the document was never produced.
I did not run it — `results/` is not this pass's territory to write into beyond
this one file — but the command is
`./.venv/bin/python experiments/multisection_report.py`.

I am **not** reading a finding out of the 66–413 µm spread here. Interpreting it
needs a pre-registration (what counts as replication, what spread is tolerable,
and whether a spot section's synthesised vessel positions make its `L` comparable
to a Xenium section's at all — the ingest record says they do not for the
per-cell statistics, and that question has to be settled before, not after).
What I am reporting is that **the solves are spent and the interpretation is not**.

### 2.9 Row 11 — a real 3-D solve

```
3-D DCT solve ran. shape (64, 64, 64)  ndim 3
   r_um      numeric      analytic exp(-r/L)/(4 pi D r)    ratio
   18.0   1.394548e-05     1.294136e-05                   1.0776
   36.0   6.183978e-06     5.682420e-06                   1.0883
   60.0   3.335888e-06     2.867230e-06                   1.1635

control_transverse_is_exact() -> grid (96,96,96), dx 6.0 um,
                                 max_rel_abs_error = 4.727583579881759e-16, exact = True
```

A 3-D steady transport solve on a 64³ grid, against the closed-form 3-D point
source. The 3-D→2-D reduction control is exact to 4.7e-16 on 96³. This matters
for `verification_routes` R10: that route scores the Secomb rat network as
geometry only ("scores a solved field? **no**"). It does not have to be. A
solver and a 3-D vessel-painting routine (`dimensionality.py:580 build_network`,
`:623 paint_network`) are both present.

**The residual is real and is why this is PARTIALLY.** `core/domain.py:22` —
*"A **2D** finite-difference grid"* — and `reconstruct/registration_3d.py`'s own
docstring: *"there is no z-axis anywhere in `tmesim.core`, and this module does
not add one there."* The shipped engine is 2-D. The analysis layer is not.

### 2.10 Row 12 — the CODEX OME-TIFF

```
read_ome_metadata('data_cache/greenwald_gbm_codex/ZH1007_NEC_B.ome.tif') ->
  shape_cyx (43, 7330, 8205)   axes CYX   43 channels   3 pyramid levels
  um_per_px 0.5066976883697085
  filename_metadata_mismatch: "OME-XML TiffData/UUID blocks self-reference ['ZH1017_INF_B.ome.tif'], which does NOT inclu..."
```

920 MB on disk, read. The adapter even carries its own provenance warning about
a filename/metadata mismatch in the source file, which is the kind of thing a
"not ingestible" verdict would never have surfaced.

### 2.11 Row 14 — the honest half of the tail90 story

The fix landed:

```
experiments/grid_convergence_gbm.py:210   E_UNIT_STORE = .../gbmgrid_e_unit.npz
experiments/grid_convergence_gbm.py:421   def store_e_unit(...)
experiments/grid_convergence_gbm.py:645   store_e_unit(k, e_unit)   # persist BEFORE scoring
```

The store did not:

```
$ ls results/sweep/gbmgrid_e_unit.npz results/sweep/gridconv_breast_e_unit.npz
ls: results/sweep/gbmgrid_e_unit.npz: No such file or directory
ls: results/sweep/gridconv_breast_e_unit.npz: No such file or directory
```

So `tail90_falsification.md`'s cost estimate — 1.5–2.5 h of ladder re-runs — is
**still the correct cost of the first run**, and free only thereafter. Anyone
quoting "the fix landed, so it's free now" would be wrong. This is the row that
most sharply distinguishes *a capability existing in code* from *a capability
existing in data*, and it is why every row above was checked against the
filesystem and not only against the source.

### 2.12 Rows 15–25 — the blockers that hold

Checked and confirmed still real, with the check that confirmed it:

* **Row 15, PS scale.** `drugs.py:134 PS_SCALE_S_INV = 1.0` unchanged.
  `calibrate_ps_scale()` returns `k0 = 0.16121537317661533` with
  `disagreement_fold = 8.703923900118907` and `irreducible_residual = 2.95×`
  over `n_anchors = 2`. The scorecard row **already calls this function and
  prints the disagreement**, so the NOT TESTABLE verdict is current, not stale.
* **Row 16, pH.** No `ph_field` / `extracellular_ph` / `pHe` anywhere in
  `physics/`, `reconstruct/` or `core/`. `molecular/ion_trapping.py:60` and
  `molecular/immune_sink.py:40` both declare *"no shipped code path imports this
  module"*; the importer grep returns only `validation/ion_trapping_impact.py`.
  `basicity('doxorubicin')` returns `pka_source='not_in_table'`,
  `trapping_ratio=1.0` — the module is present and the table does not cover the
  flagship compound.
* **Row 17, platforms.** `ingest.py:318 detect` dispatches to exactly three
  detectors and its own failure message enumerates exactly three layouts.
* **Row 18, TAM re-release.** Absent from `molecular/immune_sink.py`.
* **Row 19, capacity.** Absent. `brain_slice_anchor.md` is the standing FAIL.
* **Row 21, ECM.** `reconstruct/real_fields.py:245 ecm_density_from_cells`
  delegates to `:72 _kde_on_grid`. The vascular analogue was promoted to an
  area-conserving discrete field on 2026-08-14; ECM was not.
* **Rows 22–25.** Data and structural absences, unchanged. Do not re-propose.

---

## 3. Caveats that must travel with the recoveries

### 3.1 Row 5's obvious test is spent

`L_PARACELLULAR_LEAK = 0.238` is the **mean of the three fitted values**
(0.211 / 0.230 / 0.272) from levetiracetam, cefazolin and methotrexate —
exactly the three drug/assay pairs `glioma_regional_pk.md` §4 measured. Scoring
the barrier term against Oh 2026 or Blakeley 2008 is therefore **circular**,
under the same rule `transport_anchor_prereg.md` §0.1 applies to Pruijn: a
dataset that selected a mechanism is spent as a test of that mechanism. The
capability is real; the free test is not. Row 5 is NOW POSSIBLE **only against
out-of-sample regional data**, and `barrier.py` itself flags the disruption field
as *"a documented proxy, not a measurement"*.

### 3.2 Row 4 is not a licence to flip a default

`sink_form="michaelis_menten"` being available is not evidence that it is right.
The value of row 4 is that the two forms can now be **A/B'd with everything else
held identical**, which makes any difference attributable to the functional form
and to nothing else. That is the attribution `verification_routes` R9 wants and
currently prices at 2–3 days of new work.

### 3.3 Nothing here touches the withdrawn claims

Cross-compound ranking remains withdrawn (2026-08-21), absolutely. Absolute kill
fraction remains not a defensible output. None of the nine NOW POSSIBLE rows
creates a path back to either, and running all nine would not. No row above is
clinical, diagnostic, "digital twin", personalised-medicine, or an AI/ML claim.
Tissue registers by filesystem path throughout — `data_cache/gbmspace_visium/...`,
`data_cache/greenwald_gbm_codex/...` — never by upload.

### 3.4 Concurrent workstreams are already moving on some of this — check before starting

This pass ran alongside others. Files that did not exist when
`results/verification_routes_2026-09-20.md` was written (mtime **03:19:22**) and
do exist now, by mtime alone:

| file | mtime | bears on |
|---|---|---|
| `results/mms_transport_verification.md` | 03:21:48 | R1 — whose text says *"**No `results/` document records an MMS outcome**"*. That sentence was true for two and a half minutes |
| `experiments/visium_section_level.py` | 03:27:51 | R12 / rows 9–10 below — **this is the same territory as my ★2 and ★3. Check it before starting either** |
| `experiments/packing_density_anchor.py` | 03:36:03 | R4, `verification_routes`' own #1-ranked route |
| `results/packing_density_anchor.md` | 03:39:24 | as above |

**I have not read any of them and claim nothing about their contents** — they
belong to other workstreams and this file did not touch them. They are recorded
here only because (a) their existence is a verifiable fact that a reader of
`verification_routes` needs, and (b) `experiments/visium_section_level.py`
overlaps two of my three recommendations, so acting on them without checking
would duplicate someone's work.

It also sharpens §0.1's point rather than blunting it. The gap between a blocker
being written and its capability existing has now been observed at **two and a
half minutes**. The problem was never that people forget over time. It is that
parallel threads do not reconcile, and the claim surface has no mechanism that
makes them.

### 3.5 Two of these are already-run work, not new work

Rows 10 and 12 are not "tests that could be run". They are **tests that were
run and never written down**. That is a different and cheaper defect than the
other seven, and it should be fixed by writing, not by computing.

---

## 4. The three highest-value NOW POSSIBLE items

Ranked by what a day buys, not by how interesting the physics is.

### ★ 1 — Rows 3 + 4 together: the oxygen field in mmHg, under both sink forms

**What it buys.** `verification_routes` R9 is the repository's best shot at
converting the oxygen story from a *coefficient* failure into a *field* failure,
and it currently prices at 2–3 engineer-days on the assumption that the mmHg
mapping is the hard part. **It is not — the mapping shipped on 2026-07-29 and
the Michaelis–Menten alternative shipped in the same commit.** Both were verified
executing in §2.3.

That collapses R9 from "build the conversion, then score" to "score, twice". And
the second run is the one that matters: `oxygen_calibration_landscape.md` §4
item 2 *derived on paper* that a first-order-everywhere sink "will systematically
under-predict hypoxia", and the engine reports **HYPOXIC = 0.0000 at every grid
spacing** (`CLAIM_INVENTORY` C5) — the exact predicted signature. Running the
same section under `sink_form="linear"` and `sink_form="michaelis_menten"`, with
every other input byte-identical, is a **controlled test of that derivation**.
A difference is the functional form. There is no other candidate.

**Cost.** ~1 engineer-day, against R9's 2–3. The solver, the conversion, the
alternative sink and the sections are all in place; the work is the
pre-registered distributional statistic (median, sub-2.5 mmHg fraction,
sub-10 mmHg fraction) and the band, with Nozue 1997's 1.5–5.6 mmHg
inter-laboratory spread as the floor the bar must clear.

**What a FAIL says.** Exactly what R9 says it would, plus attribution: if the
linear arm fails and the MM arm does not, the repository has a measured
confirmation of a mechanism it derived analytically and has been carrying as a
suspicion since 2026-07-28.

### ★ 2 — Row 10: register the seven sections that are already solved

**What it buys.** `CLAIM_INVENTORY.md` §6 ranks "a second section, same tumour
type" **first** among everything that would change the honest assessment, and
calls it "the cheapest, largest improvement available". Seven sections are
already computed, at the pinned `alpha`, with `status = ok`, and **no `results/`
file mentions them**. The n = 1-in-tissue caveat currently attached to nearly
every surviving claim comes off the ones that earn it, and stays on the ones that
do not.

**Cost.** Zero solves, and the renderer that turns the JSON into
`results/multisection_generalisation.md` is already written
(`experiments/multisection_report.py`, §2.8). The remaining work is a
pre-registration — what counts as replication, decided before looking at the
spread — and the spot-vs-Xenium commensurability question. Call it 1–2 days
against `verification_routes` R12's 4–6 engineer-days, because R12 budgets for
running what has already been run.

**What a FAIL says.** R12 already says it, and says it well: a claim that holds
on section 1 and not on sections 2–7 would be this project's first retraction
driven by *tissue* variability rather than configuration, and it would be the
honest bound on what one section licenses. That failure is available today for
the price of writing it down.

### ★ 3 — Row 9: the 30 GBM-Space sections, which are not blocked

**What it buys.** It turns item 2's n = 7 into n = 37, from 30 additional
patients' sections, at **whole-transcriptome depth (36,601 genes)** rather than a
319-gene panel — which means `mgmt_from_panel` is available on every one of them
(§2.7), and so is a real ECM field over 11 genes and a vascular field over six
endothelial markers. `verification_routes` R12 prices this at "+~3 [days] for an
`.h5ad` reader". The reader shipped 2026-07-28, the conversion has already been
run, all 30 carry the `.gbmspace_materialized` marker, and `ingest.detect()`
returns `spaceranger_outs`. **The three days are already spent.**

**Cost.** ~0 days of loader work. The real cost is compute (30 reconstructions
and solves) and the same pre-registration item 2 needs. Note the honest limit up
front: these are spot sections, so `spatial_spread_valid` is `False` on every
one of them and the per-cell statistics stay refused — which is correct, and is
not a reason to decline the section-level quantities.

---

## 5. What should be done about the mechanism, not the instances

The 27 rows are symptoms. The mechanism is that **`results/` has no guard on the
understatement side**, and the date evidence says the failure is fast rather than
slow: median gap one day, three cases where the capability preceded the
blocker, and one case (row 9) where a document written **today** restated a
blocker that has been false for 54 days.

The overstatement side is covered — `claim_surface_audit_2026-08-30.md`,
`verdict_inversion_census.md`, `grid_dependence_audit.md`, the withdrawal banner
log. Nothing scans for a document saying *the engine cannot X* while `X` is
importable. That check is mechanical in a way most claim-surface checks are not:
a sentence containing "has no" / "cannot" / "not implemented" plus a module or
symbol name is a testable assertion, and `import` is the test.

I am **not** proposing to build it here — this pass wrote one file and changed no
code. But it is worth recording that the cheapest version of this guard is a
grep for capability-shaped negations in `results/` that names a `tmesim.` symbol,
checked against `importlib`, and that it would have caught rows 3, 4, 5, 6, 7,
9, 12 and 13 — eight of the sixteen recoverable-or-partial rows — on the day
each was written.

And one thing this sweep confirms about the existing guards: `results/` and
`archive/` being immutable dated records is **correct and should not change**.
Rows 26 and 27 are the proof — both were closed by the repository's own banner
convention, which annotates without rewriting. The fix for a stale blocker is a
new dated file, exactly like this one, not an edit to the old one.

---

## 6. Closing calibration guard

```
$ .venv/bin/python -c "from tmesim.validation.calibration import calibrate_alpha; print(repr(calibrate_alpha()))"
END calibrate_alpha() = 22.223477711172563
ASSERT OK
```

Identical to the opening value in §0. No code under `tmesim/`, `experiments/`,
`app/`, `docs/` or `data_cache/` was modified in this pass, and no file in
`results/` or `archive/` other than this new one was created or touched.
