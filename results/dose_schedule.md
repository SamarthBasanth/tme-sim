# DOSE & SCHEDULE on ONE REAL TUMOR

> ## STATUS: this is a TWO-ARM document. Read the arm label on every number.
>
> **The transvascular delivery law is an OPEN decision that has not been made**,
> so this report carries both candidate laws side by side and labels every
> affected number with the arm it came from. Neither arm is "history".
>
> | arm | delivery law | steady state | status as of this writing |
> |---|---|---|---|
> | **exchange OFF** | one-way leak, `S = PS·Cp` | `PS·Cp/(clearance+binding)` — **unbounded above plasma** | **the SHIPPED default** — `tmesim/core/config.py`, `transvascular_exchange: bool = False`; verified at runtime for this run, `PDEConfig().transvascular_exchange` → `False`. Every unlabelled number elsewhere in this repo is this arm. **Known to be physically wrong**: with no back-flux, tissue concentration can exceed the plasma level, which passive delivery cannot do. |
> | **exchange ON** | gradient-driven, `S = PS·(Cp − C)` | `PS·Cp/(PS+clearance+binding) ≤ Cp` **by construction** | `PDEConfig(transvascular_exchange=True)`. Physically correct *form*, same `PS`, no new parameter. **Not the default**, and its absolute scale is conditional on an unjustified implicit `k₀ = 1 s⁻¹` (`results/transvascular_exchange_revalidation.md` §6.2). |
>
> **The prior version of this file carried a conditional-retraction banner.
> That banner is superseded by this rewrite:** the body no longer states the
> one-way result as if it were the only result. Both arms were measured for
> this rewrite from one code path, in one process per stage, with **nothing but
> the delivery law changed** (`experiments/dose_schedule_arms.py`).
>
> On the runs below the one-way arm puts temozolomide at **1.10× plasma
> exposure on average and 3.07× at maximum**, and cisplatin at **2.61× / 11.55×**.
> The exchange arm puts them at **0.49× / 0.81×** and **0.61× / 0.96×** — inside
> the bound, as its algebra requires. *Every one-way figure in that sentence is
> a physical impossibility for a passively delivered drug; no exchange figure
> is.* That is measured at **four** grid resolutions in §8 and holds at all of
> them.
>
> ### ⚠ A HEADLINE CLAIM OF THIS DOCUMENT DIES UNDER THE CORRECTED LAW
>
> §6 of every earlier version led with:
>
> > "Delivery is not the binding constraint here … **at high enough dose >95 %
> > of tumour cells clear IC50 for both drugs**"
>
> **That is FALSE for temozolomide under the corrected law, and it is not a
> near miss.** At the top of the swept range, 30× nominal:
>
> | temozolomide, IC50 coverage at 30× | exchange OFF (shipped) | exchange ON (corrected) |
> |---|---:|---:|
> | dx = 12 µm (primary) | **97.80 %** | **3.43 %** |
> | across dx = 16 / 12 / 10 / 8 µm | 97.66 – 97.98 % | 2.73 – 3.72 % |
>
> The claim survives under the shipped one-way law with ~3 percentage points of
> margin, and fails under exchange by **~92 percentage points**. The grid
> oscillation is ±0.5 pp (§8) — **two orders of magnitude smaller than the
> gap**, so this is a physics result and not a mesh artefact. For **cisplatin
> the claim survives under BOTH laws** (99.87–99.88 % at 30× at every dx).
>
> The honest headline is therefore: *at high enough dose >95 % of tumour cells
> clear IC50 for **cisplatin under either law and for temozolomide only under
> the one-way law**.* It is not a statement about "both drugs", and it never
> should have been quoted without the law attached.
>
> ### What else changes, and what does not
>
> * **SURVIVES, with its number roughly halved.** Every schedule conclusion
>   (§4–5): fractionation stays small, non-zero, monotone in *n*, still exactly
>   the terminal residual, and the bolus-reaches-deeper refutation survives with
>   a *larger* relative margin. Schedule-is-third-order-vs-dose is reinforced.
> * **DEGRADES, and this is stated rather than smoothed over.** The
>   bound-vs-measured agreement that is this module's strongest evidence gets
>   materially *worse* under exchange: measured/bound goes **1.0073 → 1.1194**
>   (§4.4). The one-way law has one loss rate; exchange has `k_loss + PS(x,y)`,
>   which varies in space, so the single-τ derivation the bound rests on loses
>   its premise. Agreement to 0.7 % becomes disagreement by 12 %.
> * **RESCALED, not invalidated — the previous banner overstated this.** §3's
>   delivery-geometry table is regenerated below **in full, both arms**. Its
>   qualitative claims all survive; its numbers move, in the direction the
>   previous banner named (exchange removes drug fastest exactly where the
>   source is strongest, so it flattens the field: temozolomide penetration
>   index **0.5112 → 0.7140, +39.7 %**).
> * **Cisplatin's §6 verdict INVERTS** — sensitivity-limited under the one-way
>   law, dose-limited under exchange. That is now presented as two arms rather
>   than as a retraction (§6.2).
>
> ### The silent-no-op trap, and how it was defeated
>
> `experiments/dose_schedule.py` — the module that produced every earlier
> version of this file — is **the module on which the whole silent-no-op
> failure mode was first measured**, and is declared `INERT` in
> `tmesim/validation/propagation.py::DECLARATIONS`. Its `unit_transport` helper
> calls `simulate_drug_transport` **directly**, so
> `PDEConfig(transvascular_exchange=True)` never reaches the solver there:
>
> | run | `auc_sum` |
> |---|---:|
> | exchange OFF | 1.554447793e+09 |
> | **flag set, `exchange_k` NOT passed** | **1.554447793e+09 — bit-identical** |
> | honest ON (`exchange_k` passed) | 6.86608671e+08 |
>
> **Setting the flag on this path and comparing would have been a comparison of
> a run to itself**, and it would have looked clean. So this rewrite does not
> set the flag; `experiments/dose_schedule_arms.py` injects
> `exchange_k = exchange_coefficient(work, drug)` explicitly — byte-for-byte the
> field `sim.engine.simulate_response` builds at its step 4b — at a **single**
> call site that always names `exchange_k=`. It is registered as a
> **non-production** `WiringPoint` in `tmesim/validation/propagation.py`
> (`production=False`), which `tests/test_propagation.py` requires.
>
> Four independent checks, all of them run rather than asserted:
>
> 1. **Counted at the boundary.** Every solve ran inside
>    `propagation.instrument()`. `exchange_k` crossed the function boundary on
>    **1 of 1** solver calls in every ON solve and **0 of 1** in every OFF
>    solve. A solve that failed this raised rather than reporting.
> 2. **The arms are not bit-identical.** Required *before* any comparison was
>    reported, using `propagation.digest` (a **bit-level** structural hash, not
>    `allclose`). Temozolomide `b69a5772…` → `1a8d22d6…`; cisplatin
>    `0c3e9ac3…` → `1a9ebc1e…`. Asserted in code, not eyeballed.
> 3. **Reconciled against the PRODUCTION wiring.** The one declared
>    `PROPAGATED` path — the full `sim.engine.simulate_response` with
>    `PDEConfig(transvascular_exchange=True)`, which builds its own `PS` field —
>    was re-run end to end at **8 points, 4 per arm**, and reproduced this
>    script's fast-path kill with **|Δ| = 0.0e+00 everywhere** (§1.1). A
>    hand-injected field that differed from the engine's in any way could not
>    do that.
> 4. **Reconciled against an independent module.** The ON-arm 1× kills and
>    tissue/plasma ratios reproduce `results/mechanism_attribution.md` §3, which
>    measured the exchange arm through a *different* script and a *different*
>    code path: temozolomide `0.0072232488` here vs `0.007223248782126657`
>    there; mean tumour AUC `3.3553e5` vs `3.355e5`; ratios 0.49× / 0.81× vs
>    0.49× / 0.81×.
>
> ### Telling a REAL zero from a SOLVER-ARTEFACT zero — by construction
>
> **This file has been burned by a `0.00e+00` once already** (see the
> static-source retraction below), and under the corrected law temozolomide's
> leg-A kill spread prints `0.00e+00` again at the reference dose. A real zero
> and an artefact zero are *cosmetically identical in a table*, so they are
> separated mechanically here, not by eye:
>
> | test | OFF arm | ON arm | what it decides |
> |---|---|---|---|
> | the three schedules' AUC fields are **bit-level distinct** | ✅ 3 distinct digests | ✅ 3 distinct digests | if they were identical the zero would BE the artefact. They are not, so the inputs genuinely differ |
> | **continuous** observable: mean tumour exposure spread | 0.7724 % | **0.3792 %** | not quantised by counting cells — non-zero and monotone in *n* in both arms |
> | kill spread vs the **counting quantum** `1/23,812 = 4.20e-05` | 8.40e-05 = 2 quanta | **0.00e+00 = 0 quanta** | the ON zero is *below the quantum*, i.e. arithmetic of counting, not of the solver |
> | **positive control**: the superseded STATIC source, run deliberately | 0.00e+00 | 0.00e+00 | the artefact, reproduced on purpose next to the real measurement |
> | up-dose ladder recovers the effect | 8.4e-05 → 1.26e-03 | 0.00e+00 → 2.52e-04 | the effect reappears once it clears the quantum |
>
> **The ON zero is a counting-quantum zero, not a solver zero**, and the
> evidence is the digest row, which is a property of the fields and not of the
> table. The up-dose ladder is retained for exactly this reason and must not be
> dropped.
>
> **One thing the earlier hazard note got optimistic about:** in the ON arm the
> ladder is *itself* heavily quantised — 1× gives 0 quanta, 3× gives 1, 10× and
> 30× both give 6 — so it is no longer strictly monotone in dose. Under the
> corrected law the kill-fraction readout of the schedule effect is at the edge
> of what cell counting can express at all, and **the continuous exposure column
> is the only trustworthy readout of §4–5 in the ON arm.**

> ## THIS VERSION ALSO SUPERSEDES THE STATIC-SOURCE VERSION OF THIS FILE
>
> **What changed.** Every version before the time-resolved rewrite was produced
> with a transport solver driven by a *static* source: the engine collapsed the
> whole multi-dose plasma profile to a single time-averaged level before the
> tissue ever saw it. Two schedules with the same cumulative plasma exposure
> were therefore **bit-identical inputs**, and the reported kill spread across
> 1 / 3 / 10 fractions, and across a 192-fold plasma-half-life range, was
> `0.00e+00` — *exactly* zero, to machine precision.
>
> **That zero was an arithmetic artefact of the solver, not a result about this
> tumor.** The old file said so in its caveats but still drew the conclusions
> "fractionation changes nothing" and "sustained infusion does NOT reach further
> here". Neither was earned: the simulator was structurally incapable of
> answering either question. Both are retracted and replaced in §4–5.
>
> **What is done instead.** The source is rescaled by `plasma(t)` at every
> timestep (`PDEConfig.time_resolved=True`; the solver's `source_profile`
> argument), so tissue concentration genuinely rises and falls and each schedule
> gets its own transport solve. Cumulative delivered drug is held *exactly*
> equal between the compared paths. The static path remains the default
> elsewhere in the repo and is bit-for-bit unchanged.

> **Vascular source / units / grid.** `vascular_mode = discrete`
> (`reconstruct.real_fields._DEFAULT_VASCULAR_MODE`, promoted in
> [`results/vessel_promotion.md`](vessel_promotion.md)) · **alpha = 148.6082**,
> the e-folding-length convention `L = d_half / ln 2` · per-cell sensitivity
> mode = **`determinant`** (MKI67 proxy for cisplatin; temozolomide has no
> determinant on this panel and stays at 1.0 — flagged UPPER BOUND in the body)
> · **dx = 12 µm primary**, with dx = 16 / 10 / 8 µm measured in §8.
>
> Every absolute exposure / AUC / coverage / kill below is on the **discrete-source
> unit** and is NOT comparable with the KDE-era numbers this file carried before
> that regeneration (the tissue-integrated source total is 3.36× smaller; that is
> a change of *unit*, not of physics).

**Predicted dose-response, fractionation and infusion-duration behaviour for
temozolomide and cisplatin on the real 10x Xenium glioblastoma section.** This
is the question the simulator exists to inform — not *which* drug, but *how
much* and *how often* — asked of one real patient's tissue.

![dose & schedule](dose_schedule.png)

> **The figure is the exchange-OFF arm only.** It is generated by
> `experiments/dose_schedule.py`, which is the INERT path and therefore cannot
> draw the ON arm. Read the figure as the shipped-default column of the tables
> below, not as a summary of this document.

> ## READ THIS FIRST
>
> * **Research prototype. NOT clinical dosing guidance.** Nothing here should
>   influence the care of any patient, and no number here is a dose recommendation.
> * **Zero outcome validation.** This Xenium section is an *untreated* resection.
>   There is no measured drug response for either drug on this tissue, so every
>   kill fraction below is a model prediction with nothing to score it against —
>   **under either delivery law.**
> * **A kill quoted without its delivery law is ambiguous**, and for cisplatin
>   the two laws give opposite operational verdicts (§6.2). Quote the arm.
> * **Absolute kills are conditional on the sensitivity model.** For temozolomide
>   the per-cell sensitivity determinant (MGMT) is *not on this panel*, so
>   sensitivity is uniform 1.0 — i.e. temozolomide resistance is **unmodelled**,
>   and its absolute kill is correspondingly optimistic. For cisplatin
>   sensitivity is an MKI67-derived *mechanistic proxy*, not a measured
>   chemosensitivity.
> * **Absolute kills and every "% of ceiling" additionally inherit a
>   non-monotone grid dependence through the cell-state label** (§8.3). The kill
>   ceiling itself oscillates 0.5305–0.5502 across dx = 16/12/10/8, and the
>   published dx = 12 sits at the **maximum** of that oscillation. **IC50
>   coverage does not** — it is a pure exposure statistic and does not touch the
>   cell-state label at all, which is why the headline claim above is stated in
>   coverage rather than in kill.
> * **Pharmacodynamics is still a single static Hill map on cumulative exposure.**
>   Time resolution was added to *transport*, not to the kill model. There is no
>   inter-fraction DNA repair, no cell-cycle redistribution, no regrowth and no
>   normal-tissue recovery — and those, not transport, are what give
>   fractionation its clinical benefit. The §4–5 effect is a *lower bound on how
>   much schedule could matter*.
> * The perivascular gradient the depth numbers rest on partly **inherits the
>   shape of the vascular-density source field** — the confound documented in
>   `tmesim/validation/vessel_profile.py`. Depth figures are order-of-magnitude,
>   and §3.1 shows the source-null column is resolution-limited in both arms.
> * **~90 % of tumor cells are labelled QUIESCENT by the oxygen-banding model**,
>   carrying a 0.5× kill multiplier. That single labelling assumption dominates
>   every absolute kill and the kill ceiling below, in both arms.

---

## 1. Methods

* **Tissue — one real tumor.** The open 10x Xenium GBM section: 40,887 measured
  cells (23,812 tumor), reconstructed **once per grid** and reused for every
  sweep and both arms, so all differences are attributable to dose/schedule/law
  and not to the tissue. At the primary `dx = 12 µm` (grid 424×675): mean
  vascular density 0.0268, mean ECM 0.0759, mean oxygen 0.3520. After the oxygen
  solve **0.0 % of tumor cells are hypoxic and 89.96 % quiescent**.
* **Vascular geometry.** Distance from the nearest endothelial cell over tissue
  (PRIMARY direct-cell metric): median **44.1 µm**, p90 **85.3 µm**, max
  **789.3 µm**. The legacy p95 density mask resolves to 2.68 % of tissue area.
* **Drugs.** `calibrated_drugparams(name, alpha = 148.6082)` then
  `apply_barrier_to_drugparams` — the identical compilation pipeline
  `real_drug_panel.py` uses. **The drug compilation is delivery-law-independent:**
  every parameter in §1.2 is bit-identical between the arms, which is the point —
  the arms differ in the transport operator and in nothing else.
* **Numerics.** `dt = 10 s`, convection on. Dose-response window 6 h;
  schedule / half-life window 12 h, held EXACTLY constant across every schedule
  and half-life compared. Seed 0, fully deterministic (verified: two OFF runs of
  the schedule stage on different days returned bit-identical digests).
* **Grid, and a CORRECTION.** Earlier versions of this file stated that
  "`dx = 10 µm` is deliberately not used — it violates the upwind advection CFL
  condition at this `dt`". **That is false and is retracted here.** Measured
  directly on this section at `dt = 10 s`:

  | dx (µm) | grid | max(\|vx\|+\|vy\|) (µm/s) | CFL = max·dt/dx | solver |
  |---:|---|---:|---:|---|
  | 16 | 318×507 | 0.662716 | 0.4142 | OK |
  | 12 | 424×675 | 0.633649 | 0.5280 | OK |
  | **10** | 508×810 | 0.662942 | **0.6629** | **OK — the retracted claim is false** |
  | 8 | 635×1012 | 0.751314 | 0.9391 | OK |
  | 7 | 725×1157 | 0.668378 | 0.9548 | OK |
  | 6 | 846×1349 | 0.677291 | **1.1288** | **refused by the solver** |

  The true floor on this section is **dx = 7 µm**, not 10. All four of
  16/12/10/8 µm were actually run for §8; none was refused. The stale sentence
  most likely predates the tissue-mask switch from `occupancy` to `geometric`,
  which changed the velocity field.
* **Plasma driver.** `PDEConfig.time_resolved = True`; the transport source is
  rescaled by `plasma(t)` at EVERY timestep. The plasma curve is sampled on the
  solver's own step grid, whose mean is exactly the `plasma_mean` the static path
  used, so **both paths deliver exactly the same total drug**.
* **Dose axis.** Nominal plasma peak = 50 µM (the repo-wide `PKConfig` default),
  swept 0.05×–30× with a 4 h plasma half-life. The *same* nominal level is used
  for both drugs so the axis means one thing; it is roughly temozolomide's
  clinical order of magnitude but **supratherapeutic for cisplatin**.

### 1.1 One PDE solve per plasma-curve SHAPE — exact under BOTH laws

The transport PDE is **linear** in the source and starts from `C = 0`. Adding
`exchange_k` adds `PS(x,y)·C` to the reaction term, which is **still linear and
still independent of the source amplitude**, so the operator `A = I − dt·L + dt·K`
is unchanged by dose and the shape × amplitude split survives the corrected law
untouched:

* **Amplitude is free.** All 7 dose points and the matched-exposure /
  matched-administered-dose / matched-peak legs of a given schedule share ONE
  solve and are recovered by exact rescaling.
* **Shape is not free.** 1 vs 3 vs 10 fractions, and a 0.25 h vs a 48 h
  half-life, are DIFFERENT shapes and each needs its own solve.

That linearity argument is *verified rather than asserted*, in both arms. The
**full** `sim.engine.simulate_response` — the one declared PROPAGATED path, which
builds its own `PS` field from `PDEConfig.transvascular_exchange` — was re-run
end to end at 8 points:

| verification point | arm | full-engine kill | fast-path kill | abs delta | `exchange_k` reached solver |
|---|---|---|---|---|---|
| temozolomide dose 1× (6 h window) | OFF | 0.0183940870 | 0.0183940870 | 0.0e+00 | 0 / 1 |
| cisplatin dose 1× (6 h window) | OFF | 0.2550394759 | 0.2550394759 | 0.0e+00 | 0 / 1 |
| temozolomide schedule n=10 matched-AUC (12 h) | OFF | 0.0117167815 | 0.0117167815 | 0.0e+00 | 0 / 1 |
| temozolomide half-life 0.25 h matched-exposure (12 h) | OFF | 0.0118847640 | 0.0118847640 | 0.0e+00 | 0 / 1 |
| temozolomide dose 1× (6 h window) | **ON** | 0.0072232488 | 0.0072232488 | 0.0e+00 | **1 / 1** |
| cisplatin dose 1× (6 h window) | **ON** | 0.1401394255 | 0.1401394255 | 0.0e+00 | **1 / 1** |
| temozolomide schedule n=10 matched-AUC (12 h) | **ON** | 0.0044935327 | 0.0044935327 | 0.0e+00 | **1 / 1** |
| temozolomide half-life 0.25 h matched-exposure (12 h) | **ON** | 0.0044935327 | 0.0044935327 | 0.0e+00 | **1 / 1** |

Largest discrepancy anywhere: **0.0e+00**. The four OFF rows reproduce this
file's previously published verification table exactly, so the baseline is the
same baseline; the four ON rows are what proves the hand-injected `exchange_k`
is the production field.

### 1.2 Drug parameters as compiled onto this tissue — identical in both arms

| drug | MW (Da) | D (µm²/s) | IC50 (µM) | hill | BBB | `k_loss` (1/s) | 1/k | sensitivity determinant | median sensitivity (tumor) | kill ceiling |
|---|---|---|---|---|---|---|---|---|---|---|
| temozolomide | 194.2 | 278 | 606.05 | 1.20 | 0.1254 | 3.019e-03 | 331 s | none (uniform 1.0) | 1.000 | 0.5502 |
| cisplatin | 298.0 | 236 | 24.55 | 1.20 | 0.8577 | 8.786e-03 | 114 s | MKI67 | 0.589 | 0.3254 |

**What the exchange law adds — and it is not a constant.** The ON arm adds
`PS(x,y) = vascular_density · bbb_penetrance` to the loss rate. On this section
`PS` is **zero over most of the tissue** (the source is the discrete vessel
field, so its median and p95 are both exactly 0) and large only at vessels:

| drug | `PS` mean / p95 / max (1/s) | `k_eff = k_loss + PS`: min / mean / max | `τ_eff = 1/k_eff`: max / mean / min |
|---|---|---|---|
| temozolomide | 3.364e-03 / 0.000 / 0.1254 | 3.019e-03 / 6.383e-03 / 0.1284 | 331 s / **157 s** / 7.8 s |
| cisplatin | 2.301e-02 / 0.000 / 0.8577 | 8.786e-03 / 3.180e-02 / 0.8665 | 114 s / **31 s** / 1.2 s |

This is the single most important structural consequence of the correction:
**the one-way law has ONE loss rate and the exchange law has a spatially varying
one, spanning a factor of 43 (temozolomide) and 99 (cisplatin) across the
tissue.** Every closed form in this document that assumed a single exponential —
the analytic `√(D_eff/k)`, the `t95 = ln(20)/k` equilibration time, and the
`τ/T` schedule bound — has no single-valued counterpart in the ON arm. Where a
closed form is quoted for the ON arm below it is evaluated at the tissue-mean
`k_eff` and labelled as such; where a measurement exists it is preferred.

### 1.3 What actually sets the ceiling — the decomposition (law-independent)

The kill ceiling is the kill fraction reached as exposure → ∞. It is a property
of the sensitivity model and this tissue's oxygenation, **not of delivery**, so
it is identical in both arms — which is exactly why it does not rescue the
temozolomide coverage failure.

**temozolomide** — ceiling 0.5502

| tumor-cell band | fraction of tumor cells | resistance multiplier | mean asymptotic kill probability | contribution to ceiling |
|---|---|---|---|---|
| VIABLE | 10.0 % | ×1.00 | 1.000 | 0.1004 |
| HYPOXIC | 0.0 % | ×0.80 | 0.000 | 0.0000 |
| QUIESCENT | 90.0 % | ×0.50 | 0.500 | 0.4498 |

**cisplatin** — ceiling 0.3254

| tumor-cell band | fraction of tumor cells | resistance multiplier | mean asymptotic kill probability | contribution to ceiling |
|---|---|---|---|---|
| VIABLE | 10.0 % | ×1.00 | 0.593 | 0.0595 |
| HYPOXIC | 0.0 % | ×0.35 | 0.000 | 0.0000 |
| QUIESCENT | 90.0 % | ×0.50 | 0.296 | 0.2659 |

*The temozolomide HYPOXIC multiplier reads ×0.80 rather than the ×0.35 earlier
versions printed, because `cells.pkpd.effective_hypoxia_resistance` honours a
per-drug override for temozolomide. It is immaterial here — the band holds 0.0 %
of tumor cells and contributes exactly 0.0000 to the ceiling in every arm and at
every grid — but the earlier number was wrong and is corrected rather than left.*

---

## 2. Dose-response on this patient's tissue — both arms

`IC50 coverage` = fraction of tumor cells whose *time-averaged tissue exposure*
(`AUC / window`, the quantity the engine's Hill mapping consumes) is at or above
the drug's own IC50. `IC50 depth` = the deepest distance from a vessel at which
the binned exposure profile still clears IC50 (`+` = still above IC50 at the far
end of the profiled range, i.e. a lower bound; `—` = below IC50 even at the
vessel wall). The profiled range ends at **148.6 µm** (p99.5 of
distance-to-vessel), so `138+` means "still above IC50 in the last usable bin".

### 2.1 temozolomide

**exchange OFF — the shipped default**

| dose (× nominal) | plasma peak (µM) | window-mean plasma (µM) | kill fraction | % of kill ceiling | tumor cells ≥ IC50 | tumor cells ≥ 0.1·IC50 | median exposure / IC50 | IC50 depth (µm) |
|---|---|---|---|---|---|---|---|---|
| 0.05× | 2.5 | 1.55 | 0.00055 | 0.1 % | 0.00 % | 0.00 % | 0.003 | — |
| 0.1× | 5.0 | 3.11 | 0.00139 | 0.3 % | 0.00 % | 0.00 % | 0.006 | — |
| 0.3× | 15.0 | 9.33 | 0.00466 | 0.8 % | 0.00 % | 0.00 % | 0.017 | — |
| 1× | 50.0 | 31.09 | 0.01839 | 3.3 % | 0.00 % | 2.95 % | 0.056 | — |
| 3× | 150.0 | 93.28 | 0.06161 | 11.2 % | 0.00 % | 97.80 % | 0.168 | — |
| 10× | 500.0 | 310.95 | 0.18986 | 34.5 % | 2.95 % | 99.87 % | 0.559 | — |
| 30× | 1500.0 | 932.85 | 0.35486 | 64.5 % | **97.80 %** | 99.87 % | 1.678 | 138+ |

**exchange ON — the corrected law**

| dose (× nominal) | plasma peak (µM) | window-mean plasma (µM) | kill fraction | % of kill ceiling | tumor cells ≥ IC50 | tumor cells ≥ 0.1·IC50 | median exposure / IC50 | IC50 depth (µm) |
|---|---|---|---|---|---|---|---|---|
| 0.05× | 2.5 | 1.55 | 0.00008 | 0.0 % | 0.00 % | 0.00 % | 0.001 | — |
| 0.1× | 5.0 | 3.11 | 0.00055 | 0.1 % | 0.00 % | 0.00 % | 0.003 | — |
| 0.3× | 15.0 | 9.33 | 0.00172 | 0.3 % | 0.00 % | 0.00 % | 0.008 | — |
| 1× | 50.0 | 31.09 | 0.00722 | 1.3 % | 0.00 % | 0.00 % | 0.025 | — |
| 3× | 150.0 | 93.28 | 0.02394 | 4.4 % | 0.00 % | 3.43 % | 0.076 | — |
| 10× | 500.0 | 310.95 | 0.09067 | 16.5 % | 0.00 % | 99.86 % | 0.255 | — |
| 30× | 1500.0 | 932.85 | 0.22976 | 41.8 % | **3.43 %** | 99.87 % | 0.764 | **—** |

**The coverage failure, stated plainly.** Under the corrected law temozolomide
**never clears IC50 anywhere in this tumour at any dose in the swept range** —
not at 30× nominal, not at the vessel wall. The `IC50 depth` column is `—` in
every row of the ON table, where the OFF table reaches 138+ µm at 30×. The
mechanism is visible in the median column: the median cell sits at **0.764× IC50
at 30×** under exchange versus 1.678× under the one-way law, i.e. the corrected
law moves the median tumour cell from comfortably above its IC50 to below it,
and the 3.43 % that do clear are the perivascular tail.

This is a **headline claim of this document dying**, not a rescaling: 97.80 % →
3.43 % is not "the same conclusion with a smaller number", it is the opposite
conclusion. Under the one-way law the operational reading was *"the drug can
reach essentially the whole tumour once enough is given"*. Under exchange the
reading is *"no achievable dose in this sweep gets temozolomide to a
pharmacologically meaningful concentration anywhere on this section"*.

**Saturation.** OFF: kill has NOT saturated within the swept range; coverage
passes 90 % of tumor cells by 30×. ON: kill has not saturated either (0.00722 →
0.22976, still climbing), but coverage never passes 4 %. **Under exchange the
sweep is dose-limited all the way up and never becomes delivery-sufficient.**

**What 30× buys over 1×.** OFF: kill 0.01839 → 0.35486 (19.29× relative). ON:
kill 0.00722 → 0.22976 (31.81× relative) — the corrected law makes the drug
*more* dose-responsive over this range, because it pushes the whole population
further down the steep limb of the Hill curve.

**At the approximate clinical plasma band.** Published Cmax for temozolomide is
of order 30–70 µM, i.e. **0.6×–1.4×** nominal. That band lies between the
measured 0.3× and 3× rows: OFF kill 0.00466–0.06161, ON 0.00172–0.02394, and
**IC50 coverage is 0.00 % at every measured dose in that bracket in both arms.**
*Order-of-magnitude literature reference for scale only — not fitted, not
patient-specific, not a dose recommendation. No number in this paragraph is
interpolated; the bracketing measured rows are quoted instead.*

### 2.2 cisplatin

**exchange OFF — the shipped default**

| dose (× nominal) | plasma peak (µM) | window-mean plasma (µM) | kill fraction | % of kill ceiling | tumor cells ≥ IC50 | tumor cells ≥ 0.1·IC50 | median exposure / IC50 | IC50 depth (µm) |
|---|---|---|---|---|---|---|---|---|
| 0.05× | 2.5 | 1.55 | 0.03624 | 11.1 % | 0.00 % | 91.49 % | 0.159 | — |
| 0.1× | 5.0 | 3.11 | 0.06980 | 21.5 % | 0.50 % | 99.55 % | 0.319 | — |
| 0.3× | 15.0 | 9.33 | 0.16110 | 49.5 % | 44.53 % | 99.87 % | 0.956 | 45 |
| 1× | 50.0 | 31.09 | 0.25504 | 78.4 % | **99.55 %** | 99.87 % | 3.185 | 138+ |
| 3× | 150.0 | 93.28 | 0.29783 | 91.5 % | 99.87 % | 99.87 % | 9.556 | 138+ |
| 10× | 500.0 | 310.95 | 0.31375 | 96.4 % | 99.87 % | 99.87 % | 31.852 | 138+ |
| 30× | 1500.0 | 932.85 | 0.31786 | 97.7 % | 99.87 % | 99.87 % | 95.555 | 138+ |

**exchange ON — the corrected law**

| dose (× nominal) | plasma peak (µM) | window-mean plasma (µM) | kill fraction | % of kill ceiling | tumor cells ≥ IC50 | tumor cells ≥ 0.1·IC50 | median exposure / IC50 | IC50 depth (µm) |
|---|---|---|---|---|---|---|---|---|
| 0.05× | 2.5 | 1.55 | 0.00727 | 2.2 % | 0.00 % | 0.00 % | 0.039 | — |
| 0.1× | 5.0 | 3.11 | 0.01491 | 4.6 % | 0.00 % | 4.67 % | 0.079 | — |
| 0.3× | 15.0 | 9.33 | 0.04918 | 15.1 % | 0.00 % | 99.65 % | 0.236 | — |
| 1× | 50.0 | 31.09 | 0.14014 | 43.1 % | **4.67 %** | 99.87 % | 0.787 | — |
| 3× | 150.0 | 93.28 | 0.23476 | 72.2 % | 99.65 % | 99.87 % | 2.361 | 138+ |
| 10× | 500.0 | 310.95 | 0.29284 | 90.0 % | 99.87 % | 99.87 % | 7.869 | 138+ |
| 30× | 1500.0 | 932.85 | 0.31178 | 95.8 % | **99.87 %** | 99.87 % | 23.608 | 138+ |

**Cisplatin is the drug for which the headline claim SURVIVES.** 99.87 % of
tumor cells clear IC50 at 30× under both laws, and at every grid tested (§8).
What the corrected law changes is *where on the dose axis* that happens: IC50
coverage at the nominal dose falls **99.55 % → 4.67 %**, saturation moves from
~3× to ~10× (90 % of the 30× kill), and the 30×/1× kill ratio goes **1.25 →
2.22**. That is a shift of the same curve along the dose axis, not the loss of a
result — and it is the inversion §6.2 is about.

---

## 3. Delivery geometry: what the exposure field actually looks like — both arms

Because the field is linear in dose, **dose moves the amplitude of the
perivascular exposure profile and nothing else** — the shape is fixed by
`D_eff`, the loss rate, and the vascular geometry. This is exact in **both**
arms: adding `exchange_k` changes the loss rate but not its dependence on the
source amplitude.

**The profile is NOT a single exponential, and saying otherwise would be the
main way to mislead with this figure.** It is a steep perivascular drop sitting
on a flat, non-zero FLOOR that persists to the far edge of the tissue. The floor
is *not* a diffusion front: the engine's delivery term is the **continuous**
vascular-density field, so a pixel far from the thresholded vessel mask is still
fed by diffuse low-density source. All three characterisations are therefore
reported together.

**This is the table §4.3 of `results/transvascular_exchange_revalidation.md`
listed as INVALIDATED. It is regenerated here in full rather than deleted — and
the finding is that it was RESCALED, not invalidated.** Every qualitative claim
below survives; the numbers move.

| drug | arm | empirical half-distance | **source-null half-distance** | near-field exponential L (d ≤ p90) | near-field fit R² | floor-model L | floor level (% of perivascular) | far-field plateau / perivascular | analytic √(D_eff/k) | equilibration t95 (**measured**) | penetration index |
|---|---|---|---|---|---|---|---|---|---|---|---|
| temozolomide | **OFF** | n/a (never halves) | 12 µm (= 1 bin) | 305.5 µm | 0.9923 | 283.8 µm | 0.0 % | 0.611 | 296 µm | 1010 s (16.8 min) | **0.5112** |
| temozolomide | **ON** | n/a (never halves) | 12 µm (= 1 bin) | 438.7 µm | 0.9934 | 396.7 µm | 0.0 % | 0.695 | 204 µm ‡ | 520 s (8.7 min) | **0.7140** |
| cisplatin | **OFF** | 106.7 µm | 12 µm (= 1 bin) | 160.7 µm | 0.9892 | 149.6 µm | 2.6 % | 0.427 | 160 µm | 360 s (6.0 min) | **0.4284** |
| cisplatin | **ON** | 136.1 µm | 12 µm (= 1 bin) | 243.1 µm | 0.9959 | 229.4 µm | 0.0 % | 0.547 | 84 µm ‡ | 150 s (2.5 min) | **0.7210** |

‡ *`√(D_eff/k)` has no single value in the ON arm — `k_eff` spans a factor of 43
(temozolomide) and 99 (cisplatin) across the tissue (§1.2). The value printed is
evaluated at the tissue-mean `k_eff` and is a summary, not a length that exists
anywhere in particular. The OFF values are exact.*

**Perivascular amplitude, per µM of window-mean plasma** (the quantity the
shape statistics are normalised by): temozolomide **1.283 → 0.544** (−57.6 %),
cisplatin **3.487 → 0.732** (−79.0 %). Exchange takes most of its bite out of
the vessel-adjacent amplitude, which is exactly why it flattens the profile.

### 3.1 Six things to take from that table

1. **The corrected law FLATTENS the field, by ~40 %, and this is the single
   most robust number in the section.** Penetration index
   **0.5112 → 0.7140 (+39.7 %)** for temozolomide and **0.4284 → 0.7210
   (+68.3 %)** for cisplatin. The temozolomide figure reproduces at every grid
   tested to within 0.1 pp on three of four rungs (+39.6 / +39.7 / +37.1 /
   +39.6 %, §8). Mechanism: exchange removes drug fastest exactly where the
   source is strongest, so the perivascular peak is cut far harder than the
   avascular floor.
2. **Every fitted decay length gets LONGER, not shorter, and the OFF→ON RATIO is
   mesh-stable even where the lengths are not.** Near-field L
   305.5 → 438.7 µm (temozolomide, ×1.436) and 160.7 → 243.1 µm (cisplatin,
   ×1.513). Across dx = 16/12/8 the temozolomide ratio is ×1.431 / ×1.436 /
   ×1.430 — stable to 0.4 % over a 2× mesh range, while the absolute lengths
   drift ~5 %. **Quote the ratio, not the length.**
3. **The far-field plateau RISES**: 0.611 → 0.695 (temozolomide), 0.427 → 0.547
   (cisplatin). More of the field is floor under the corrected law, so the
   warning that "drug reaches X µm" is mostly a statement about where the floor
   sits relative to IC50 is *more* true after the correction, not less.
4. **The source-null column is resolution-limited and was over-read.** It comes
   out at exactly **one bin width at every resolution tested** — 16 µm at
   dx = 16, 12 µm at dx = 12, 10 µm at dx = 10, 8 µm at dx = 8. It is therefore
   not a measured length at all; it is the statement "the raw vascular-density
   field halves inside the first bin", i.e. `≤ 1 bin`. The earlier table's
   "12 µm" was the bin width. The *qualitative* claim it was cited for — that
   the measured gradient is far shallower than the source and so the steep
   near-vessel drop is inherited from where the endothelial cells are —
   survives and is if anything strengthened, because the drug profiles are
   9–37× broader than one bin in every arm.
5. **The empirical half-distance is the fragile statistic in this table and
   should not be quoted as a verdict.** It exists only when the far-field
   plateau sits below 0.5, and the plateau is not converged: temozolomide OFF
   runs 0.652 / 0.611 / 0.599 / 0.569 across dx = 16 / 12 / 10 / 8. So
   temozolomide "never halves" at dx = 16, 12 and 10 but **does** halve at
   dx = 8 (144.9 µm), and cisplatin ON "never halves" at dx = 16 but does at
   dx = 12 (136.1 µm) and dx = 8 (134.2 µm). **Existence of this number is
   decided by the mesh, in both arms.** This is the same last-bin fragility
   `results/vessel_profile.md` documents for the same statistic; it is a
   property of the readout, not of the delivery law, and it is why the two
   temozolomide `n/a` entries above must not be read as "exchange destroyed the
   penetration depth".
6. **Tissue equilibration roughly halves, and stops being a single
   exponential.** Measured t95 under a sustained source: 1010 → 520 s
   (temozolomide) and 360 → 150 s (cisplatin). The OFF measurements confirm the
   closed form to within 2 % (1010 vs `ln(20)/k` = 992 s; 360 vs 341 s), which
   validates the measurement method; **the ON values have no closed form to
   check against**, because `k_eff` is spatially varying. That is the same fact
   as §1.2 and it is what breaks the §4 bound.

For scale: median distance-to-vessel on this section is 44.1 µm, p90 85.3 µm,
max 789.3 µm. Bins beyond the p90 hold few pixels each, so the far tail is noisy
as well as flat — the near-field exponential is deliberately restricted to
d ≤ p90 (7 usable bins).

---

## 4. Schedule: fractionation at MATCHED TOTAL DOSE — both arms

One 12 h course, 4 h plasma half-life, split into 1 / 3 / 10 fractions. Three
legs, because "matched dose" can mean two different things and the distinction
turns out to be the whole story:

* **A — matched cumulative exposure** (the key comparison): the per-dose peak is
  rescaled so the *window-averaged plasma level*, i.e. total plasma AUC, is
  identical across schedules.
* **B — matched administered dose**: the sum of the per-dose peaks is held fixed
  (`peak = nominal / n`).
* **C — matched peak** (a control): each fraction reaches the full nominal peak,
  so total exposure *rises* with the number of doses.

> **DROPPED for compute budget:** the schedule and half-life legs were NOT run
> for cisplatin, in either arm — each plasma-curve SHAPE needs its own transport
> solve at the schedule window. Only temozolomide is reported. What was dropped
> is named rather than extrapolated.

### 4.1 temozolomide, exchange OFF (the shipped default)

| schedule | per-dose peak, leg A (µM) | window-mean plasma, A | **kill, A (matched exposure)** | tumor ≥ IC50, A | window-mean plasma, B | kill, B | window-mean plasma, C | kill, C |
|---|---|---|---|---|---|---|---|---|
| 1 dose (bolus) | 50.00 | 21.0444 | **0.011801** | 0.00 % | 21.0444 | 0.011801 | 21.0444 | 0.011801 |
| 3 doses q4.0h | 20.59 | 21.0444 | **0.011759** | 0.00 % | 17.0359 | 0.009575 | 51.1077 | 0.032505 |
| 10 doses q1.2h | 7.04 | 21.0444 | **0.011717** | 0.00 % | 14.9463 | 0.008105 | 149.4628 | 0.100790 |

Leg B falls 0.011801 → 0.008105 (31.32 % relative). Leg C rises 0.011801 →
0.100790. **Leg A kill spread 8.40e-05 (= 2 counting quanta), tissue-exposure
spread 0.7724 %.**

### 4.2 temozolomide, exchange ON (the corrected law)

| schedule | per-dose peak, leg A (µM) | window-mean plasma, A | **kill, A (matched exposure)** | tumor ≥ IC50, A | window-mean plasma, B | kill, B | window-mean plasma, C | kill, C |
|---|---|---|---|---|---|---|---|---|
| 1 dose (bolus) | 50.00 | 21.0444 | **0.004494** | 0.00 % | 21.0444 | 0.004494 | 21.0444 | 0.004494 |
| 3 doses q4.0h | 20.59 | 21.0444 | **0.004494** | 0.00 % | 17.0359 | 0.003612 | 51.1077 | 0.012431 |
| 10 doses q1.2h | 7.04 | 21.0444 | **0.004494** | 0.00 % | 14.9463 | 0.002940 | 149.4628 | 0.041282 |

Leg B falls 0.004494 → 0.002940 (34.58 % relative). Leg C rises 0.004494 →
0.041282. **Leg A kill spread `0.00e+00`, tissue-exposure spread 0.3792 %.**

> ### ⚠ THAT `0.00e+00` IS NOT THE RETRACTED ARTEFACT — and here is the proof,
> not the reassurance
>
> The kill level fell 2.6× under the corrected law (0.011801 → 0.004494), and
> the schedule effect is a fixed *relative* change in exposure, so the absolute
> kill spread fell with it — from 2 counting quanta to **under one**. The
> quantum is `1/23,812 = 4.20e-05`; a spread of 4.4e-05 in exposure-equivalent
> terms simply cannot be expressed by counting whole cells.
>
> | evidence | value | verdict |
> |---|---|---|
> | the three schedules' AUC fields, bit-level digests | `7251a51f…`, `c3bb1f30…`, `a62f6065…` | **three DISTINCT fields** — the inputs are genuinely different, so this cannot be the "same run twice" artefact |
> | mean tumour exposure (continuous, not quantised) | 4.556461e5 → 4.546535e5 → 4.539250e5 µM·s | **non-zero and strictly monotone decreasing in n** |
> | exposure spread | **0.3792 %** | non-zero |
> | kill spread ÷ quantum | 0.00e+00 ÷ 4.20e-05 = 0 quanta | below what cell counting can express |
> | **positive control** — the superseded STATIC source, same arm | all three schedules share ONE field, digest `9e1055b2…`; spread `0.00e+00` | **this is what the artefact looks like: ONE digest, not three** |
>
> The distinguishing test is the digest row and it is structural: the artefact
> produces **one** field for three schedules, the real result produces **three**.
> That is checked in code (`schedule_fields_all_distinct`, which raises if it
> fails), not read off a table.

### 4.3 The same comparison against the STATIC source — the bug, visible, in both arms

`mean tumor exposure` is the continuous quantity the kill map consumes; it is
not quantised by cell counting, so it is the honest place to read the size of
the schedule effect. **The static column is bit-identical down the whole table
BY CONSTRUCTION — that is the bug, deliberately reproduced.**

| arm | schedule | mean tumor exposure, TIME-RESOLVED (µM·s) | mean tumor exposure, STATIC (µM·s) | Δ vs static | kill, time-resolved | kill, static | terminal residual C(T)/delivered | penetration index |
|---|---|---|---|---|---|---|---|---|
| OFF | 1 dose (bolus) | 1.050561e+06 | 1.044857e+06 | +0.5459 % | 0.011801 | 0.011759 | 0.2315 % | 0.511733 |
| OFF | 3 doses q4.0h | 1.045917e+06 | 1.044857e+06 | +0.1014 % | 0.011759 | 0.011759 | 0.6674 % | 0.511408 |
| OFF | 10 doses q1.2h | 1.042508e+06 | 1.044857e+06 | −0.2248 % | 0.011717 | 0.011759 | 0.9873 % | 0.511139 |
| **ON** | 1 dose (bolus) | 4.556461e+05 | 4.544130e+05 | +0.2714 % | 0.004494 | 0.004494 | **0.1014 %** | 0.714967 |
| **ON** | 3 doses q4.0h | 4.546535e+05 | 4.544130e+05 | +0.0529 % | 0.004494 | 0.004494 | **0.2924 %** | 0.714342 |
| **ON** | 10 doses q1.2h | 4.539250e+05 | 4.544130e+05 | −0.1074 % | 0.004494 | 0.004494 | **0.4326 %** | 0.713892 |

Relative exposure spread across the fractionation range: **0.7724 % (OFF)** and
**0.3792 % (ON)**, monotonically decreasing with the number of fractions in both
arms. The static column's spread is exactly 0 in both.

> **The terminal-residual column is computed with the CORRECT per-voxel
> elimination rate in each arm.** Under the one-way law the whole loss is the
> scalar `k_loss = clearance + binding`; under exchange it is
> `k_loss + PS(x,y)`. Using the scalar in the ON arm — which is what the
> obvious code reuse would have done — understates eliminated drug and
> **overstates the residual by ~2.3×** (it would have printed 0.2294 / 0.6598 /
> 0.9745 % instead of 0.1014 / 0.2924 / 0.4326 %). Both are recorded in the run
> payload (`residual_fraction` and `residual_fraction_kloss_only`); the OFF arm
> is bit-identical between them, which is the check that the corrected form
> reduces to the published one.

### 4.4 The bound — which SURVIVES in the OFF arm and DEGRADES in the ON arm

Integrating the transport equation over `[0, T]` with `C(0) = 0`, for `A = AUC`
and `M = k − div(D grad ·) + div(v ·)`:

```
    A = M⁻¹ ( S_total − C(T) ),      S_total = ∫₀ᵀ S dt
```

`S_total` is **identical** for any two schedules at matched cumulative plasma
exposure. Therefore:

1. **The entire schedule effect is the terminal residual `C(T)`.** Nothing else
   about the timing can reach the answer. **True in both arms** — the derivation
   uses only linearity, which exchange preserves.
2. **Earlier delivery wins**, so exposure decreases monotonically with the
   number of fractions. **Confirmed in both arms** (the tables above).
3. **The effect is bounded by the tissue time constant over the window.** This
   is where the arms part company. The bound is measured, not assumed: the
   constant-infusion (static-source) residual **saturates** it by construction,
   so it is quoted as measured rather than as `τ/T`.

| | exchange OFF | exchange ON |
|---|---:|---:|
| analytic `τ/T = (1/k_loss)/12 h` | 0.7667 % | *no single value — `k_eff` spans 43×* |
| **measured constant-infusion residual (the bound)** | **0.7668 %** | **0.3387 %** |
| **measured exposure spread, 1 → 10 fractions** | **0.7724 %** | **0.3792 %** |
| measured ÷ bound | **1.0073** | **1.1194** |

**In the OFF arm the analytic prediction, the measured bound and the measured
spread all agree to 0.7 %** — three numbers, one derived before it was measured.
That agreement is the strongest single piece of evidence in this module that the
time-resolved source is doing real physics.

**In the ON arm that agreement degrades to 12 %, and the measured spread
EXCEEDS the bound.** This is stated rather than smoothed over, because the
previous banner described it as "the same relation". It is not the same
relation. The reason is structural and is the same fact as §1.2 and §3.1(6):
the `τ/T` argument needs *a* τ, and under exchange the loss rate varies by a
factor of 43 across the tissue, so the tissue-integrated constant-infusion
residual is an average over modes rather than an envelope over them. **Under the
corrected law the schedule effect is still small, still monotone and still
exactly the terminal residual — but it is no longer cleanly bounded by a
one-line closed form, and any future claim of "bounded by τ/T" must be
restricted to the one-way arm.**

### 4.5 Absolute kill spread up the dose axis — the ladder, kept deliberately

The schedule effect is a fixed *relative* change in exposure, so the *absolute*
kill spread grows with dose. Rescaling leg A up the dose axis is free (exposure
is exactly linear in amplitude, §1.1). **Every static-source row is exactly zero
in both arms — that is the control, not an omission.**

| leg-A dose (× the matched-exposure level) | arm | 1 dose | 3 doses | 10 doses | spread, TIME-RESOLVED | spread, STATIC | spread ÷ quantum |
|---|---|---|---|---|---|---|---|
| 1× | OFF | 0.011801 | 0.011759 | 0.011717 | **8.40e-05** | 0.00e+00 | 2 |
| 3× | OFF | 0.041492 | 0.041408 | 0.041156 | **3.36e-04** | 0.00e+00 | 8 |
| 10× | OFF | 0.138166 | 0.137872 | 0.137452 | **7.14e-04** | 0.00e+00 | 17 |
| 30× | OFF | 0.295187 | 0.294431 | 0.293927 | **1.26e-03** | 0.00e+00 | 30 |
| 1× | **ON** | 0.004494 | 0.004494 | 0.004494 | **0.00e+00** | 0.00e+00 | **0** |
| 3× | **ON** | 0.015958 | 0.015958 | 0.015916 | **4.20e-05** | 0.00e+00 | 1 |
| 10× | **ON** | 0.060726 | 0.060642 | 0.060474 | **2.52e-04** | 0.00e+00 | 6 |
| 30× | **ON** | 0.173316 | 0.173190 | 0.173064 | **2.52e-04** | 0.00e+00 | 6 |

**Read the ON ladder honestly.** In the OFF arm the ladder does what it was
built to do: the spread lifts cleanly off the quantum and grows monotonically,
2 → 8 → 17 → 30 quanta. In the ON arm it is quantisation almost all the way up —
0 → 1 → 6 → 6 quanta, and **not monotone between 10× and 30×**. The ladder still
demonstrates that the effect is real and non-zero, but under the corrected law
**the kill-fraction readout of the schedule effect is at the edge of what
counting 23,812 cells can express, and §4.3's continuous exposure column is the
only quantitatively trustworthy readout of this section in the ON arm.**

### 4.6 Interpretation

**(a) RETRACTED, and it stays retracted: "fractionation changes nothing".** The
static-source version reported a leg-A spread of exactly `0.00e+00` and
concluded that this model predicts fractionation changes nothing. That was a
model artefact — the solver was fed the time-average of the plasma curve, so 1
dose and 10 doses were literally the same input array. It is now measured in
both arms, and in both arms the underlying exposure effect is non-zero and
monotone.

**(b) Why the effect is small HERE.** Tissue equilibration is minutes against a
12 h window in both arms — and it is *shorter* under the corrected law (measured
t95 1010 → 520 s, §3), so `τ/T` is smaller and **the corrected physics makes
schedule matter LESS, not more.** This is the opposite of the intuition that
motivated the leg, and it is the direction both arms agree on.

**(c) Leg B is a window artifact, in both arms.** Kill falls with fractionation
there only because later fractions have less of the 12 h observation window left
in which to be delivered and counted. Extend the window past full washout and
leg B collapses onto leg A. It is a measurement-window effect, not a schedule
effect — and note it is *slightly larger* under exchange (31.32 % → 34.58 %),
because faster tissue clearance makes the truncation bite harder.

**(d) What is still missing, and it is the big one.** Time resolution was added
to *transport*. The kill map remains a **static Hill function of cumulative
exposure**: no inter-fraction DNA repair, no cell-cycle redistribution, no
regrowth, no resistance induction, no normal-tissue recovery. Those are the
mechanisms that give fractionation its actual clinical benefit, and they are
pharmacodynamics, not transport. Read §4–5 as **the transport-only floor on how
much schedule can matter**, in either arm.

---

## 5. Infusion duration proxy: plasma half-life, bolus vs sustained — both arms

Half-life is swept 0.25–48 h (a 192-fold range) with the peak rescaled to hold
total window exposure fixed, plus a matched-peak control. The hypothesis under
test is the standard transport-limited one: *a sustained low concentration
should reach further than a high, brief one.*

> **DROPPED for compute budget:** not run for cisplatin, in either arm.

### 5.1 temozolomide, both arms

| arm | plasma half-life (h) | peak, matched exposure (µM) | **kill (matched total exposure)** | tumor ≥ IC50 | mean tumor exposure (µM·s) | terminal residual C(T)/delivered | penetration index | far-field plateau | kill (matched peak) |
|---|---|---|---|---|---|---|---|---|---|
| OFF | 0.25 | 697.48 | **0.011885** | 0.00 % | 1.053028e+06 | 0.0000 % | 0.511937 | 0.6120 | 0.000504 |
| OFF | 4 | 50.00 | **0.011801** | 0.00 % | 1.050561e+06 | 0.2315 % | 0.511733 | 0.6117 | 0.011801 |
| OFF | 48 | 22.92 | **0.011759** | 0.00 % | 1.045534e+06 | 0.7032 % | 0.511382 | 0.6112 | 0.028893 |
| **ON** | 0.25 | 697.48 | **0.004494** | 0.00 % | 4.561734e+05 | 0.0000 % | 0.715299 | 0.6963 | 0.000084 |
| **ON** | 4 | 50.00 | **0.004494** | 0.00 % | 4.556461e+05 | 0.1014 % | 0.714967 | 0.6959 | 0.004494 |
| **ON** | 48 | 22.92 | **0.004494** | 0.00 % | 4.545600e+05 | 0.3104 % | 0.714288 | 0.6952 | 0.010919 |

Kill spread across the 192-fold half-life range at matched exposure: **1.26e-04
(OFF, = 3 quanta)** and **0.00e+00 (ON, = 0 quanta)**. Exposure spread:
**0.7167 % (OFF)** and **0.3549 % (ON)**. The ON zero is the same
counting-quantum zero as §4.2, and it is separated from the artefact by the same
digest test: the three half-life shapes produce three **distinct** AUC fields
(`6c1f7155…`, `7251a51f…`, `cf4c1611…`) in the ON arm, and the static control
produces one.

**Every half-life now has its OWN exposure field, so depth is measurable — in
both arms.** Under the static source there was a single AUC field for the whole
sweep, so "depth is unchanged across half-life" was true by construction and
carried no information.

### 5.2 The naive transport-limited hypothesis is REFUTED, in both arms, and MORE strongly under the corrected law

The penetration index moves from the bolus-like 0.25 h half-life to the
sustained 48 h one:

| arm | penetration index, t½ = 0.25 h | t½ = 48 h | shift | relative |
|---|---:|---:|---:|---:|
| OFF | 0.511937 | 0.511382 | **−0.000554** | −0.108 % |
| **ON** | 0.715299 | 0.714288 | **−0.001010** | **−0.141 %** |

**Both arms say the BOLUS reaches relatively deeper, not the sustained
infusion**, and the corrected law says it **1.8× more strongly in absolute terms
and 1.3× more strongly in relative terms.** The same ordering holds on the
schedule leg (penetration index falls monotonically with the number of
fractions: −0.000593 OFF, −0.001076 ON).

**What the old version said, and why it was not an answer.** It reported the
depth as identical across the whole 192-fold range and concluded "sustained
infusion does NOT reach further here … it is exposure that moves depth, never
the shape of the plasma curve". But there was only ONE exposure field in that
run — the plasma shape had been averaged away before the solve. The invariance
was the input being invariant. **Retracted, and it stays retracted under both
laws.**

**Why the bolus can never lose here, under either law.** From
`A = M⁻¹(S_total − C(T))`: at matched cumulative exposure `S_total` is
identical, so the schedules differ by `M⁻¹` applied to their terminal residuals.
The bolus has washed out (`C(T) ≈ 0`, measured 0.0000 % in both arms); the
sustained infusion has not. The exposure the sustained infusion loses is
`M⁻¹C_ss = M⁻²S_rate` — the source field smoothed **twice**, i.e. a *flatter*
field than the exposure field itself. Subtracting a flatter field from a peakier
one leaves something steeper. **`M` is still a positive operator with
`exchange_k` added**, so the argument carries over intact; what changes is that
`M` is now stiffer, which is why the effect is larger.

**Do not over-read the size of this.** The shifts are in the fourth decimal
place of a normalised index, against a perivascular gradient that §3 shows is
partly inherited from the source geometry, and the penetration index itself
carries a ~0.5 %-per-rung grid drift (§8). The value of this leg is **that the
question is answerable and the answer has a sign and a magnitude that agree
across two different delivery laws** — not that infusion duration is an
actionable lever on this tumor.

---

## 6. THE VERDICT for this tumor: dose-limited or delivery-limited?

> ### Headline — and it is now TWO headlines, because the arms disagree
>
> | | exchange OFF (the shipped default) | exchange ON (the corrected law) |
> |---|---|---|
> | **temozolomide** | **DOSE-limited** over the swept range, then resistance-limited. Coverage 0.00 % → **97.80 %** by 30× | **DELIVERY-limited, categorically.** Coverage never exceeds **3.43 %** at any dose swept, and IC50 is not cleared even at the vessel wall |
> | **cisplatin** | **SENSITIVITY-limited.** 99.55 % clear IC50 already at 1×; 30× buys only 1.25× more kill | **DOSE-limited at the nominal dose.** Coverage at 1× falls to 4.67 %; 30× buys 2.22× more kill |
> | ">95 % clear IC50 at high enough dose" | **TRUE for both drugs** | **TRUE for cisplatin (99.87 %), FALSE for temozolomide (3.43 %)** |
> | "response is DOSE-limited, never delivery-limited" | holds | **FALSE — temozolomide is delivery-limited under the corrected law** |
>
> **The single sentence this document used to lead with — *"Delivery is not the
> binding constraint here … at high enough dose >95 % of tumour cells clear IC50
> for both drugs"* — is arm-dependent in its first clause and false in its
> second.** Under bounded physics delivery IS the binding constraint for
> temozolomide on this section, at every dose tested and at every grid tested.
>
> **What survives both laws, unchanged:**
> * **Schedule is a THIRD-ORDER lever.** Fractionation moves tissue exposure by
>   0.7724 % (OFF) / 0.3792 % (ON) against a dose lever that moves kill by 19×
>   (OFF) / 32× (ON) over the swept range. The corrected law makes schedule
>   matter *less*, not more (§4.6b).
> * **The bolus reaches relatively deeper than the sustained infusion**, with a
>   larger margin under the corrected law (§5.2).
> * **The kill ceiling is not a transport number** and no dose touches it. For
>   temozolomide it is 0.5502 in both arms, of which ~82 % is the quiescent band
>   alone. **This is the one part of the analysis the delivery law cannot move**
>   — and it is also the part that carries the largest grid dependence (§8.3).
> * **Cisplatin is never delivery-limited** in either arm — its coverage at 30×
>   is 99.87 % under both laws.
>
> **Load-bearing caveat, unchanged and now quantified:** the ~90 % quiescent
> fraction is a *model* label from the simulated oxygen field, not a
> measurement. §8.3 shows it is also **non-monotone in grid resolution**,
> swinging 89.96–93.90 % across dx = 16/12/10/8, and that the published
> dx = 12 sits at the favourable end. Every absolute kill and every "% of
> ceiling" in this document inherits that, in both arms.

### 6.1 temozolomide — the verdict INVERTS between the laws

**exchange OFF: DOSE-limited over the swept range, then sensitivity-limited.**
At the nominal dose 0.00 % of tumor cells clear IC50, rising to 97.80 % by 30× —
the drug *can* reach essentially the whole tumor once enough of it is given.
Kill climbs 0.01839 → 0.35486 (19.29×) and at 30× has reached 64.5 % of the
0.5502 ceiling without flattening.

**exchange ON: DELIVERY-limited, categorically.** Coverage is 0.00 % at every
dose up to 10× and **3.43 % at 30×**. The `IC50 depth` column is `—` in every
row: exposure does not clear IC50 *even at the vessel wall* at any swept dose.
The median tumour cell sits at 0.764× IC50 at 30×. Kill still climbs (0.00722 →
0.22976, 31.81×) because the Hill term is steep down there, but it climbs from
a population that is entirely sub-therapeutic.

**This reproduces, through a completely different code path, what
`results/mechanism_attribution.md` §3.2 found per-cell:** temozolomide is
97.10 % delivery-limited under the one-way law and **100.00 %** under exchange,
with the DOSE bucket emptying entirely. Two modules, two scripts, one
conclusion.

Geometry check: near-field decay length **305.5 µm (OFF) / 438.7 µm (ON)**
against median distance-to-vessel 44.1 µm and p90 85.3 µm. But see §3: with a
far-field plateau at 0.611 (OFF) / 0.695 (ON) of the perivascular value, large
depths are the **floor** clearing IC50, not a diffusion front reaching further —
and under exchange the floor never clears IC50 at all.

### 6.2 cisplatin — the verdict INVERTS, and this is the fourth independent reproduction

**exchange OFF: NEITHER dose- nor delivery-limited — SENSITIVITY-limited.**
Already at the nominal dose 99.55 % of tumor cells see an exposure above IC50;
raising the dose 30-fold moves kill only 0.25504 → 0.31786 (1.25×). What caps
the response is the ceiling (0.3254).

**exchange ON: DOSE-limited at the nominal dose.** Coverage at 1× falls
**99.55 % → 4.67 %**, kill at 1× goes **78.4 % → 43.1 % of ceiling**, saturation
moves from **~3× to ~10×**, and the 30×/1× kill ratio goes **1.25 → 2.22**.

**Which arm to believe reduces to which delivery law is right, and on that the
evidence is one-sided.** The one-way arm puts cisplatin at **2.61× plasma
exposure on average and 11.55× at maximum** on this section. Passive
transvascular delivery cannot concentrate a drug above its plasma level. The
exchange arm is bounded by construction (0.61× / 0.96×) and does not breach it —
**at all four grid resolutions tested** (§8.2). That is an argument about the
physics, not about this module's arithmetic, but it is the argument that decides
which column describes a real tumour.

The same inversion was reproduced independently in `results/mechanism_attribution.md`
§3.1 (RESISTANCE 96.42 % → 3.64 %, DOSE 0.44 % → 95.65 %) and in the GBM 3-D:2-D
and cross-drug ranking re-runs. **Do not quote cisplatin as sensitivity-limited
or resistance-limited in external material without naming the one-way law.**

Note that the one-way arm's own "Important qualification" already anticipated
this direction on clinical-dose grounds: cisplatin's approximate clinical Cmax
(2–6 µM) is only 0.04×–0.12× nominal, where the OFF sweep also gives 0.00 %
coverage. **Corrected physics says the dose-limited reading applies at the
nominal dose too, not only down at clinical Cmax.**

### 6.3 What this means operationally (within the model's limits)

* **Would more drug help?** Read the `% of kill ceiling` column in §2 — *and the
  arm label*. For cisplatin the two arms give opposite answers at the nominal
  dose (78.4 % vs 43.1 % of ceiling), which is the difference between "escalation
  is pure added toxicity" and "the tumour is under-exposed". That distinction is
  the one thing this simulator is positioned to make, and it is currently
  undecided pending the delivery-law decision.
* **Is the tissue transport-saturated?** DOSE moves amplitude only — exactly, in
  both arms, because the operator is linear in both. SCHEDULE moves the shape by
  under a percent in both arms. So the geometry of delivery is set by `D_eff`,
  the loss rate, and the vascular geometry — and under the corrected law the
  loss rate is no longer a constant, which is the substantive change.
* **The lever that is NOT in this analysis** is the normal-tissue side of the
  therapeutic index. Every statement above is about tumor kill alone.

---

## 7. Limitations (the full list)

1. **No outcome validation, at all.** Untreated resection; no measured response
   for either drug, under either delivery law.
2. **Not clinical dosing guidance.** Research prototype.
3. **THE DELIVERY LAW IS UNDECIDED.** The shipped default is known to be
   physically wrong (it puts tissue above plasma); the corrected law's absolute
   scale rests on an unjustified implicit `k₀ = 1 s⁻¹`
   (`results/transvascular_exchange_revalidation.md` §6.2). Neither arm is
   "the" answer and this document does not pretend otherwise.
4. **`PS` and `bbb_penetrance` are perfectly confounded** in the ON arm: scaling
   `vascular_density` by *f* gives a bit-identical `PS` field to scaling
   `bbb_penetrance` by *f*. The model cannot distinguish a wrong vessel map from
   a wrong permeability. `PS` is also still a blood-**brain**-barrier term.
5. **The `τ/T` schedule bound is an OFF-arm result.** Under exchange the loss
   rate varies by 43× across the tissue, the measured spread exceeds the measured
   bound by 12 %, and no single-τ closed form applies (§4.4).
6. **Schedule resolution is window-conditional.** The `τ/T` ratio means a longer
   window shrinks the measured schedule effect and a shorter one inflates it.
   The 12 h window was chosen for compute, not physics. It is held exactly
   constant across everything compared, so the comparisons are sound, but the
   absolute percentages must be quoted with the window attached.
7. **Pharmacodynamics is a single static Hill map** — no inter-fraction repair,
   cell-cycle redistribution, regrowth, resistance evolution, immune
   contribution or normal-tissue toxicity. Time-resolving *transport* removed a
   structural blind spot; it did not make this a fractionation model.
8. **Sensitivity model.** Temozolomide's determinant (MGMT) is not on this panel,
   so TMZ resistance is unmodelled (uniform 1.0) and its absolute kill is
   optimistic — which makes the ON-arm coverage failure a statement about an
   *optimistic* model. Cisplatin uses an MKI67 proliferation proxy.
9. **The perivascular gradient partly inherits the source-field shape**, and the
   source-null column that was cited as evidence for this is resolution-limited
   (one bin at every dx, §3.1(4)). The qualitative confound stands; the
   quantitative comparison was over-read.
10. **The empirical half-distance is mesh-fragile to the point of
    non-existence** — it exists or not depending on dx, in both arms (§3.1(5)).
    Do not quote it as a verdict.
11. **Cell-state banding dominates the ceiling AND is non-monotone in grid**
    (§8.3). 90 % of tumor cells are labelled QUIESCENT from the *simulated*
    oxygen field, not a measurement, and the label swings 89.96–93.90 % across
    the ladder.
12. **Absolute quantities do NOT inherit ratio-cancellation protection** (§8.4).
    Coverage percentages are absolute, not ratios.
13. **Scope was cut for the compute budget**: 7 dose points; schedule and
    half-life legs for temozolomide only, in **both** arms; 3 schedules; 3
    half-lives; 12 h schedule window. Nothing was extrapolated to fill the gaps —
    the gaps are named.
14. **One calibration constant** (`alpha = 148.6082`, fit once on doxorubicin
    under the ONE-WAY law) sets the perivascular sink for every drug in **both**
    arms. It has not been refitted under exchange, so the ON arm inherits an
    OFF-arm calibration. `results/transvascular_exchange_revalidation.md` §5
    documents that this is a live limitation, not a neutral choice.
15. **Single section, single patient, 2D.** No inter-patient variability, no
    out-of-plane transport, and the nominal plasma level is a shared reference
    rather than a per-drug clinical PK model.

---

## 8. Grid sensitivity — because the headline claim is a THRESHOLD crossing

A threshold claim is the most grid-fragile kind of claim there is: a few percent
of movement in a scored quantity can carry it across the line. This document's
headline is exactly that shape (">95 % of cells clear IC50"), so it is checked
at **four** grid resolutions rather than asserted at one.

**Everything in §1–7 above is dx = 12 µm.** All four rungs were run through the
same script, both arms, both drugs, with only `dx` changed
(`experiments/dose_schedule_arms.py --stage geom --dx {16,12,10,8} --skip-eq`).
`dx = 10 µm` is included because the claim that it violates CFL is false (§1);
`dx = 6 µm` is genuinely refused by the solver at `dt = 10 s` (CFL 1.1288), and
`dt` was **not** cut, so nothing below is confounded between spatial and
temporal error.

> **This is a targeted robustness check on one claim, not a convergence study.**
> The full GBM and breast convergence ladders are owned by other agents and are
> not duplicated here. What is claimed below is only: *does the verdict move?*

### 8.1 The headline claim, at four resolutions

**temozolomide, IC50 coverage at 30× nominal — the number the headline rests on**

| dx (µm) | grid | exchange OFF | exchange ON |
|---:|---|---:|---:|
| 16 | 318×507 | 97.98 % | 3.72 % |
| **12 (published)** | 424×675 | **97.80 %** | **3.43 %** |
| 10 | 508×810 | 97.95 % | 2.73 % |
| 8 | 635×1012 | 97.66 % | 3.30 % |
| **spread** | | **0.32 pp** | **0.99 pp** |

**cisplatin, IC50 coverage**

| dx (µm) | at 30×, OFF | at 30×, ON | at 1×, OFF | at 1×, ON |
|---:|---:|---:|---:|---:|
| 16 | 99.88 % | 99.88 % | 99.60 % | 4.02 % |
| **12 (published)** | 99.87 % | 99.87 % | 99.55 % | 4.67 % |
| 10 | 99.87 % | 99.87 % | 99.72 % | 4.20 % |
| 8 | 99.87 % | 99.87 % | 99.65 % | 4.93 % |

**The verdict does not move at any rung, in either direction:**

* temozolomide **OFF** clears 95 % at every rung, with **2.66–2.98 pp of
  margin**. The claim survives, but the margin is only ~9× the grid spread, so
  it should be quoted to two significant figures (**"~98 %"**) and not as
  "97.80 %".
* temozolomide **ON** fails at every rung by **91.3–92.3 pp**. The grid spread
  is 0.99 pp — **93× smaller than the gap**. *This failure cannot be a mesh
  artefact.*
* cisplatin clears 95 % at 30× at every rung under **both** laws (99.87–99.88 %,
  spread 0.01 pp).
* cisplatin's 1× inversion (99.55 % → 4.67 %) reproduces at every rung
  (99.55–99.72 % → 4.02–4.93 %).

**The movement is NON-MONOTONE**, in both arms and both drugs — 16 → 12 → 10 → 8
does not walk in one direction. A three-rung sample (16/12/8) would have looked
monotone in the ON arm and been misleading; the dx = 10 rung is what exposes it.
**No order of accuracy is extracted from these four points and none should be**:
a Richardson extrapolation on a non-monotone sequence is meaningless, and the
GBM ladder that other agents own reports the same pathology (power-law fits go
degenerate). What is claimed here is a *bound on the wobble*, not a limit.

### 8.2 The physics argument is mesh-robust — which matters more than the numbers

The argument that decides *which arm describes a real tumour* is that the
one-way law puts tissue concentration above plasma, which passive delivery
cannot do. That is checked at every rung:

| dx (µm) | temozolomide mean / max tissue AUC ÷ plasma AUC | cisplatin mean / max |
|---:|---|---|
| | **OFF** → **ON** | **OFF** → **ON** |
| 16 | 1.1001 / 2.9201 → 0.4849 / 0.7951 | 2.6012 / 10.6616 → 0.5870 / 0.9529 |
| **12** | **1.1043 / 3.0713 → 0.4900 / 0.8137** | **2.6111 / 11.5496 → 0.6097 / 0.9605** |
| 10 | 1.0905 / 3.0090 → 0.4892 / 0.8103 | 2.5784 / 11.1716 → 0.6138 / 0.9584 |
| 8 | 1.1016 / 3.0534 → 0.4915 / 0.8102 | 2.6048 / 11.1741 → 0.6181 / 0.9591 |

**The one-way arm exceeds 1.0 at every resolution and the exchange arm never
does.** No rung comes close to flipping either statement: the closest ON
approach to the bound is cisplatin's max at 0.9605, and the closest OFF approach
from above is temozolomide's mean at 1.0905. The `C ≤ Cp` bound is not a
numerical accident of the published grid.

Likewise the §3 flattening result: the temozolomide penetration-index gain is
**+39.6 / +39.7 / +37.1 / +39.6 %** across the four rungs — three of four agree
to 0.1 pp. **"Exchange flattens the field by ~40 %" is a two-significant-figure
result**, and the OFF→ON *ratio* of near-field decay lengths is stabler still
(×1.431 / ×1.436 / ×1.416 / ×1.430) than either length alone (~5 % drift).

### 8.3 Where the grid dependence really bites: the cell-state label

The largest grid effect in this document is not in transport at all. The
oxygen-banded cell-state label is **non-monotone in dx and swings by 4
percentage points**, and it multiplies straight into the kill ceiling and
therefore into every absolute kill and every "% of ceiling":

| dx (µm) | quiescent fraction of tumor cells | kill ceiling, temozolomide | kill ceiling, cisplatin |
|---:|---:|---:|---:|
| 16 | 93.34 % | 0.5333 | 0.3154 |
| **12 (published)** | **89.96 %** | **0.5502** ← ladder MAXIMUM | **0.3254** ← ladder MAXIMUM |
| 10 | 93.90 % | 0.5305 ← ladder minimum | 0.3137 ← ladder minimum |
| 8 | 90.24 % | 0.5488 | 0.3245 |

**The published dx = 12 sits at the favourable end of the oscillation for both
drugs**, exactly as the GBM convergence ladder reports for its own scored
ratios. The ceiling is anti-correlated with the quiescent fraction rung for rung,
which identifies the mechanism unambiguously: it is `F_state`, not transport.
Ceiling spread is 3.7 % (temozolomide) and 3.7 % (cisplatin) of the published
value, one-directional in the sense that the published rung is the top.

**Consequence, stated plainly: every absolute kill and every "% of kill ceiling"
in §2, §4, §5 and §6 carries a ~4 % downward revision risk from `F_state`
alone**, in both arms, before any transport error is counted. Temozolomide's
30× kill across the ladder is 0.34323 / 0.35486 / 0.33899 / 0.35327 (OFF) and
0.22132 / 0.22976 / 0.22203 / 0.22972 (ON) — a 4.7 % / 3.8 % spread with the
published rung at the top.

**This is exactly why the headline claim is stated in COVERAGE and not in
kill.** IC50 coverage is computed from the exposure field against a per-drug
constant, over a fixed set of 23,812 tumor cells; **it never touches the
cell-state label.** It is immune to `F_state` by construction, and the ladder
confirms it: coverage does not track the quiescent fraction rung for rung
(dx = 10 has the ladder's *highest* quiescent fraction and its *lowest* ON
coverage). Coverage has its own, smaller grid wobble (§8.1) from the
discretisation of the vascular source field, and that is the only grid error the
headline carries.

### 8.4 A protection this document does NOT get — and a `0.00e+00` bounded rather than asserted

Three independent studies in this repo have now found the same cancellation
structure: the GBM `dt`-control (halving `dt` moves the scored ratio by ≤0.0003 %
while provably reaching the solver and moving the exposure field by 1.2e-04),
the boundary-layer study (a domain **mean** carries essentially zero grid error,
2.1e-11, while the p99-normalised **ratio** carries +8.3 %), and the diffusivity
sweep (99.6 % cancellation against a 3.07× exposure move). The mechanism is that
p99 bath-matching removes multiplicative **scale** error exactly, leaving only
**shape** error.

**None of that protection applies here.** IC50 coverage, kill fraction, mean
tumour AUC and the tissue/plasma ratio are **absolute** quantities, not
bath-matched ratios. A stable ratio elsewhere in the repo licenses nothing about
the absolute numbers in this document, and §8.1/§8.3 are the direct measurements
that replace that borrowed reassurance.

**Bounding the numerical contribution to the `0.00e+00`, rather than asserting
it is small.** The boundary-layer study measured that under-resolving the
perivascular layer biases tissue concentration **DOWN**, with
`|rel err| = 1.427e-3 · dx^1.939 · δ^-1.145`. That is the **same direction as a
spurious zero**, so the two are actively confusable and the bound is worth
computing rather than waving at. The exchange law introduces exactly such a
layer, `δ = √(D_eff/k_eff)`, and this run measures its own `k_eff` (§1.2), so δ
is not borrowed:

| drug, ON arm | `D_eff` (µm²/s) | `k_eff` max (1/s) | δ (µm) | `\|rel err\|` at dx = 12 | at dx = 8 |
|---|---:|---:|---:|---:|---:|
| temozolomide | 265.3 | 0.1284 | 45.5 | **0.223 %** | 0.102 % |
| cisplatin | 225.2 | 0.8665 | 16.1 | **0.732 %** | 0.333 % |

*(δ evaluated at the tissue **maximum** `k_eff`, i.e. the thinnest layer on the
section and therefore the worst case; everywhere else the layer is thicker and
the error smaller.)*

So the numerical contribution to a concentration in the ON arm is bounded at
**≲0.25 % (temozolomide) and ≲0.75 % (cisplatin)** at the published grid, biased
downward. Against that:

* the temozolomide leg-A **exposure** spread is **0.3792 %** — the same order as
  the temozolomide bound, so **the exposure spread in the ON arm is NOT resolved
  with margin, and §4.4's 12 % bound-vs-measured disagreement must not be
  interpreted more finely than that.** It is a real, monotone, digest-verified
  effect; it is not a number to quote to four significant figures.
* the temozolomide coverage gap from 95 % is **~92 pp**, i.e. **~400× the
  numerical bound.** The headline failure is not numerical.

**And the `0.00e+00` itself is decided by neither of these**, which is the point
of doing it structurally: it is decided by the three schedule AUC fields having
**three distinct bit-level digests** where the artefact produces **one** (§4.2).
A downward numerical bias could shrink a spread; it cannot merge three distinct
fields into one, and it cannot manufacture the monotone ordering that the
continuous exposure column shows.

---

## 9. Reproduce, and provenance

Both arms of §2–§5, the propagation guard, the arms-differ assertion and the
real-vs-artefact zero test, in one process per stage:

```bash
# sections 2-3, both arms, both drugs (dx = 12 primary)
./.venv/bin/python -u experiments/dose_schedule_arms.py --out geom.json  --stage geom  --dx 12
# sections 4-5, both arms, temozolomide
./.venv/bin/python -u experiments/dose_schedule_arms.py --out sched.json --stage sched --dx 12
# section 1.1 engine cross-checks, both arms
./.venv/bin/python -u experiments/dose_schedule_arms.py --out eng.json --stage engine --dx 12 \
    --resume geom.json,sched.json
# section 8 grid ladder
for DX in 16 10 8; do
  ./.venv/bin/python -u experiments/dose_schedule_arms.py --out geom_dx$DX.json \
      --stage geom --dx $DX --skip-eq
done
# the flag-propagation guard itself
./.venv/bin/python -m pytest tests/test_propagation.py -q
```

The single-arm (shipped-default) path through the original generator, which
produces the figure and the OFF column only:

```bash
./.venv/bin/python experiments/dose_schedule.py
```

**`experiments/dose_schedule.py` cannot produce the ON column.** It is declared
`INERT` in `tmesim/validation/propagation.py::DECLARATIONS` — it calls
`simulate_drug_transport` directly, so setting
`PDEConfig(transvascular_exchange=True)` there is a silent no-op that returns a
bit-identical field. `experiments/dose_schedule_arms.py` is registered as a
**non-production** `WiringPoint` (`production=False`) precisely because it has
to inject `exchange_k` by hand to measure that path at all.

### Provenance

Every number in this file comes from runs executed for this rewrite, on `main`,
at `alpha = 148.6082`, seed 0, `dt = 10 s`. Nothing is interpolated,
extrapolated or carried over. AUC-field blake2b-128 digests — the evidence that
the arms are distinct runs rather than one run reported twice:

| field (dx = 12, 6 h dose window) | exchange OFF | exchange ON |
|---|---|---|
| temozolomide | `b69a57725d8c…` | `1a8d22d68996…` |
| cisplatin | `0c3e9ac35d9f…` | `1a9ebc1e5dae…` |
| `auc_sum`, temozolomide | 4.621386757e+09 | 2.050549950e+09 |
| `auc_sum`, cisplatin | 1.092751945e+10 | 2.551702175e+09 |

Schedule-shape digests (12 h window, temozolomide), the three-vs-one test of
§4.2: OFF `23caf577f65a…` / `a70163121216…` / `318fe18c7e14…` against static
`a9036609b1c5…`; ON `7251a51fdbc5…` / `c3bb1f309193…` / `a62f6065e07a…` against
static `9e1055b2f540…`.

**Four independent reproductions, recorded rather than assumed:**

1. **The OFF arm reproduces the superseded version of this file exactly** where
   the quantity is physical: every kill fraction, every IC50 coverage, every
   median E/IC50, every leg-A/B/C kill, every ladder spread, every penetration
   index, and all three previously published engine cross-checks
   (`0.0183940870`, `0.0117167815`, `0.0118847640`) match to all published
   digits.
2. **The ON arm reproduces `results/mechanism_attribution.md` §3**, measured
   through a different script and a different code path: temozolomide 1× kill
   `0.0072232488` here vs `0.007223248782126657` there; mean tumour AUC 3.3553e5
   vs 3.355e5; tissue/plasma 0.49× / 0.81× vs 0.49× / 0.81×; cisplatin 2.61× /
   11.55× (OFF) and 0.61× / 0.96× (ON) in both.
3. **The corrected terminal-residual bound reproduces the figure recorded in the
   superseded retraction banner**: constant-infusion residual **0.3387 %** ON,
   measured spread **0.3792 %**, against **0.7668 %** / **0.7724 %** OFF.
4. **Determinism**: the schedule stage was run twice, on different days, and
   returned bit-identical AUC digests for every shape in both arms.

**Two stated discrepancies against the superseded file, not reconciled away:**

* **The distance-binned profile columns of §3 have drifted**, and it is upstream
  of the delivery law. This run measures median distance-to-vessel **44.09 µm**,
  p90 **85.26 µm**, max **789.32 µm** against the superseded file's 43 / 85 /
  786 µm. The vessel-distance map itself changed slightly, so the binned profile
  statistics changed with it — in **both** arms identically. This is why
  temozolomide's OFF empirical half-distance reads `n/a` here where the
  superseded file printed 145 µm, and why the OFF plateau reads 0.611 against
  0.5628. **The OFF↔ON comparison is unaffected** — both arms are binned on the
  same map — but a cross-file comparison of those specific columns is not apples
  to apples. §3.1(5) and §8.1 show this statistic is mesh-fragile anyway.
* **The `IC50 depth` column tops out at 138+ µm here** against the superseded
  file's 146 / 150+, for the same reason: one fewer usable distance bin under
  the current binning (profiled range 148.6 µm).

**The temozolomide HYPOXIC resistance multiplier** in §1.3 reads ×0.80 here
against the superseded ×0.35, because `cells.pkpd.effective_hypoxia_resistance`
honours a per-drug override. It changes nothing — the band holds 0.0 % of tumor
cells at every grid and contributes exactly 0.0000 — but the old value was wrong.

*Generated by `experiments/dose_schedule_arms.py`: 4 grid resolutions × 2
delivery laws, 8 dose-window solves + 4 sustained-source solves at dx = 12,
14 schedule-window solves, 8 full-engine verification runs, seed 0,
deterministic. Figure `dose_schedule.png` is the exchange-OFF arm only, from
`experiments/dose_schedule.py`.*

**Not for clinical use.** This tissue is an untreated resection with no outcome
label. Every number here is a prediction of the physical model and has never
been scored against a measured treatment response — under either delivery law.
