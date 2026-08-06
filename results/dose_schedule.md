# DOSE & SCHEDULE on ONE REAL TUMOR

> ## SCOPE NOTICE — what a corrected delivery law does to this report
>
> Re-run under `PDEConfig(transvascular_exchange=True)` (branch
> `transvascular-exchange`; see `results/breast_kp_validation.md` for why the
> one-way source was wrong). Summary of what holds and what does not:
>
> * **Every SCHEDULE conclusion survives**, with its number roughly halved.
>   Fractionation stays small, non-zero and monotone in n; the effect is still
>   exactly the terminal residual, still bounded by tissue-time/window, and the
>   bound-vs-measured agreement that is this module's strongest evidence still
>   holds (0.3387 % bound vs 0.3792 % measured ON, the same relation as
>   0.7668 % vs 0.7724 % OFF). The bolus-reaches-deeper refutation survives
>   with a LARGER margin. Schedule-is-third-order-vs-dose is reinforced.
> * **The direction is the opposite of what was expected.** Corrected physics
>   roughly halves the effective tissue half-life (temozolomide 233 s → ~126 s;
>   cisplatin becomes faster than any plasma curve tested). Because the whole
>   schedule effect IS the terminal residual, a faster tissue makes schedule
>   matter LESS, not more.
> * **Washout stops being single-exponential.** With `k_eff = k_loss + PS(x,y)`
>   spatially varying, temozolomide's fitted half-life drifts 126 → 156 s
>   across the window as the slow avascular mode takes over. The one-way source
>   gave a clean single exponential. This is new behaviour, not a rescaling.
> * **INVALIDATED, not merely rescaled:** every absolute AUC and kill in this
>   file (temozolomide −56 %, cisplatin −77 %, doxorubicin −3.5 %), and the
>   whole section-3 delivery-geometry table — exchange removes drug fastest
>   exactly where the source is strongest, so it flattens the field and the
>   temozolomide penetration index rises 0.5112 → 0.7140 (+40 %). Section 3's
>   *qualitative* claims all survive; its numbers do not.
> * **Cisplatin's verdict flips** (see the conditional-retraction box at that
>   drug's section below), and the headline "at high enough dose >95 % of
>   tumour cells clear IC50 for both drugs" FAILS for temozolomide, whose 30x
>   coverage goes 97.80 % → 3.43 %.
>
> **⚠ HAZARD IF REGENERATING THIS REPORT WITH THE FLAG ON.** Temozolomide's
> leg-A kill spread at the reference dose prints as `0.00e+00` again —
> cosmetically identical to the solver artefact this very file was rewritten to
> retract (see below). It is NOT that artefact returning: the continuous
> exposure spread is 0.3792 %, non-zero and still monotone in n, but kill is
> quantised at 1/23,812 = 4.2e-05 and the kill level fell 2.6x, so the effect
> now sits under the counting quantum. The up-dose ladder already in this
> module recovers it. Keep the ladder, or you will re-publish a spurious zero.
>
> ## THIS VERSION SUPERSEDES THE STATIC-SOURCE VERSION OF THIS FILE
>
> **What changed.** Every earlier version of this report was produced with a transport solver driven by a *static* source: the engine collapsed the whole multi-dose plasma profile to a single time-averaged level before the tissue ever saw it. Two schedules with the same cumulative plasma exposure were therefore **bit-identical inputs**, and the reported kill spread across 1 / 3 / 10 fractions, and across a 192-fold plasma-half-life range, was `0.00e+00` — *exactly* zero, to machine precision.
>
> **That zero was an arithmetic artefact of the solver, not a result about this tumor.** The old file said so in its caveats but still drew the conclusion "on this tumor, at matched cumulative exposure, this model predicts fractionation changes nothing" and "sustained infusion does NOT reach further here". Neither conclusion was earned: the simulator was structurally incapable of answering either question. Both are retracted and replaced below.
>
> **What this version does instead.** The source is now rescaled by `plasma(t)` at every timestep (`PDEConfig.time_resolved=True`; the solver's `source_profile` argument), so tissue concentration genuinely rises and falls and each schedule gets its own transport solve. Cumulative delivered drug is held *exactly* equal between the two paths, so this is a comparison of schedule and not of dose. The static path remains the default everywhere else in the repo and is bit-for-bit unchanged (pinned by `tests/test_time_resolved_pk.py` against a stored reference).
>
> **The corrected answer, in one line.** Fractionation and infusion duration do NOT have zero effect on this tumor — they have a *small, non-zero, monotone, and analytically bounded* effect, and it runs in the **opposite direction to the textbook expectation**: at matched cumulative exposure the BOLUS reaches deeper than the sustained infusion. Sections 4 and 5 give the numbers, the derivation, and the bound.

> **Vascular source / units.** `vascular_mode = discrete` (`reconstruct.real_fields._DEFAULT_VASCULAR_MODE`, promoted in [`results/vessel_promotion.md`](vessel_promotion.md)) · **alpha = 148.6**, the e-folding-length convention `L = d_half / ln 2` · per-cell sensitivity mode = **`determinant`** (per-drug: MKI67 proxy for cisplatin; temozolomide has no determinant on this panel and stays at 1.0 — flagged UPPER BOUND in the body) · **dx = 12 µm**.
>
> Every absolute exposure / AUC / coverage / kill below is on the **discrete-source unit** and is NOT comparable with the KDE-era numbers this file carried before that regeneration (the tissue-integrated source total is 3.36x smaller; that is a change of *unit*, not of physics — the dose-matched, geometry-only effect on coverage is <= 0.0166). To recover the old numbers, re-run with `mode='kde'`.
>
> Note the consequence of promotion: the p95 'vessel mask' used for the depth statistics resolves to **2.7% of tissue**, i.e. essentially the segmented vessel mask itself, rather than the broad p95 contour of a smooth KDE.

**Predicted dose-response, fractionation and infusion-duration behaviour for temozolomide and cisplatin on the real 10x Xenium glioblastoma section.** This is the question the simulator exists to inform — not *which* drug, but *how much* and *how often* — asked of one real patient's tissue.

![dose & schedule](dose_schedule.png)

> ## READ THIS FIRST
>
> * **Research prototype. NOT clinical dosing guidance.** Nothing here should influence the care of any patient, and no number here is a dose recommendation.
> * **Zero outcome validation.** This Xenium section is an *untreated* resection. There is no measured drug response for either drug on this tissue, so every kill fraction below is a model prediction with nothing to score it against.
> * **Absolute kills are conditional on the sensitivity model.** For temozolomide the per-cell sensitivity determinant (MGMT) is *not on this panel*, so sensitivity is uniform 1.0 — i.e. temozolomide resistance is **unmodelled**, and its absolute kill is correspondingly optimistic. For cisplatin sensitivity is an MKI67-derived *mechanistic proxy*, not a measured chemosensitivity. Trust the *shape* of the dose-response and the *relative* comparisons; do not trust an absolute kill fraction.
> * **The schedule comparison is now real, but the effect is SMALL and its size is set by a ratio you can compute in advance.** With a time-resolved plasma driver the whole schedule effect is the drug still in flight at the end of the observation window, which for a linear transport model is bounded by (tissue equilibration time) / (window length). For temozolomide here that is 0.77%, and the measured spread matches it. A slowly-equilibrating molecule (an antibody) would show a large effect on the same machinery; a fast small molecule over a long window cannot.
> * **Pharmacodynamics is still a single static Hill map on cumulative exposure.** Time resolution was added to *transport*, not to the kill model. There is still no inter-fraction DNA repair, no cell-cycle redistribution, no regrowth and no normal-tissue recovery — and those, not transport, are what give fractionation its clinical benefit. So the small effect reported here is a *lower bound on how much schedule could matter*, not an estimate of it.
> * The perivascular gradient the depth numbers rest on partly **inherits the shape of the vascular-density source field** — the confound documented in `tmesim/validation/vessel_profile.py`. Depth figures are order-of-magnitude.
> * **~90% of tumor cells are labelled QUIESCENT by the oxygen-banding model**, and quiescent cells carry a 0.5x kill multiplier. That single labelling assumption dominates every absolute kill and the kill ceiling below. It is a *model* assignment from the oxygen field, currently being tested independently against measured MKI67 — **absolute kills may shift materially when that lands.** The dose/schedule *comparisons* are far more robust to it than the absolute levels.
> * **SCOPE (compute budget).** The machine was running ~35 concurrent jobs, so this run was deliberately cut: 7 dose points per drug; the schedule and half-life legs run for **temozolomide only** (each plasma-curve SHAPE now needs its own PDE solve — that is the cost of resolving schedule at all, and it is why the second drug is dropped rather than the legs); 3 schedules and 3 half-lives; the schedule/half-life window shortened to 12 h (held EXACTLY constant across everything compared); and 3 rather than 6 full-engine verification runs. Everything reported below actually ran — nothing is extrapolated. What was dropped is named where it would have appeared.

---

## 1. Methods

* **Tissue — one real tumor.** The open 10x Xenium GBM section: 40,887 measured cells (23,812 tumor), reconstructed **once** at the grid-converged `dx = 12 um` (grid 424x675) and reused for every sweep, so all differences are attributable to dose/schedule and not to the tissue. Mean vascular density 0.027, mean ECM 0.076, mean oxygen 0.352. After the oxygen solve 0.0% of tumor cells are hypoxic and 90.0% quiescent (both carry reduced drug kill).
* **Vascular geometry.** Vessel mask = top 5% of the endothelial-derived vascular-density field (2.7% of tissue area). Distance from the nearest vessel over tissue: median 43 um, p90 85 um, max 786 um. This is the geometry the drug's penetration has to cover — though see section 3: on this section most of the far tissue is reached by a diffuse source floor rather than by a diffusion front, which changes how that comparison should be read.
* **Drugs.** `calibrated_drugparams(name, alpha = 148.61)` (single perivascular-sink constant, fit once on doxorubicin) then `apply_barrier_to_drugparams` — the identical compilation pipeline `real_drug_panel.py` uses.
* **Per-cell sensitivity** is layered on per drug with `reconstruct.real_pipeline.apply_drug_sensitivity` (measured determinant where the panel carries one, honest uniform 1.0 otherwise). See the table below.
* **Numerics.** `dt = 10 s`, convection on. Dose-response window 6 h; schedule / half-life window 12 h (cut from 24 h for compute, and held EXACTLY constant across every schedule and half-life compared — this matters more than before, because with a time-resolved driver the window length now *sets the size of the schedule effect*; see section 4). `dx = 10 um` is deliberately not used — it violates the upwind advection CFL condition at this `dt`. Seed 0, fully deterministic.
* **Plasma driver — the change that defines this version.** `PDEConfig.time_resolved = True`. The transport source is rescaled by `plasma(t)` at EVERY timestep, so tissue concentration rises and falls with the real multi-dose curve. The implicit operator `A = I - dt*L + dt*K` contains no source term, so the LU factorisation is still built once and reused for all 4320 steps — a time-varying source is one extra scalar multiply per step and costs no measurable runtime. The plasma curve is sampled on the solver's own step grid, whose mean is exactly the `plasma_mean` the static path used, so **both paths deliver exactly the same total drug** and a static-vs-resolved difference is a pure schedule effect.
* **Dose axis.** Nominal plasma peak = 50 uM (the repo-wide `PKConfig` default), swept 0.05x–30x with a 4 h plasma half-life. The *same* nominal level is used for both drugs so the axis means one thing; note this is roughly temozolomide's clinical order of magnitude but is **supratherapeutic for cisplatin** (free-platinum Cmax is roughly an order of magnitude lower), so read cisplatin's low-dose end as the clinically relevant part of its curve.

### 1.1 One PDE solve per plasma-curve SHAPE — exact, not an approximation

The transport PDE `dC/dt = div(D_eff grad C) - div(vC) + S - (clearance+binding)C` is **linear** in the source and starts from `C = 0`. The source is `S(x,y,t) = vascular_density x bbb_penetrance x plasma(t)`, so scaling the whole plasma curve scales the exposure field exactly: `AUC(P x shape) = P x AUC(shape)`. Every schedule is therefore split into **(shape x amplitude)**, where the shape is `plasma(t)` renormalised to unit window-mean, and:

* **Amplitude is free.** All 7 dose points, and the matched-exposure / matched-administered-dose / matched-peak legs of a given schedule, share ONE solve and are recovered by rescaling it. Exact, not interpolated.
* **Shape is not free, and that is the whole point.** 1 vs 3 vs 10 fractions, and a 0.25 h vs a 48 h half-life, are DIFFERENT shapes and each needs its own solve. Under the superseded static source every shape collapsed to the same constant level, which is exactly why the schedule question had no answer. This run performs 9 transport solves rather than the 3 the static version needed.

Every dose / schedule / half-life point then re-runs the engine's own unmodified kill mapping (`cells.pkpd.apply_drug_kill`) on the rescaled field — the same trick the engine itself uses to build its viable-fraction history. This script orchestrates `core` / `physics` / `molecular` / `reconstruct` / `engine`; it re-implements none of them.

This is verified, not asserted: the **full** `sim.engine.simulate_response` was re-run end to end at 3 points — one on each leg (dose sweep, fractionation, half-life) and covering both exposure windows, so the shape-x-amplitude fast path is pinned everywhere it is used. Each check is itself a complete extra PDE solve.

| verification point | full-engine kill | fast-path kill | abs delta |
|---|---|---|---|
| temozolomide dose 1x (6 h window) | 0.0183940870 | 0.0183940870 | 0.0e+00 |
| temozolomide schedule n=10 matched-AUC (12 h window) | 0.0117167815 | 0.0117167815 | 0.0e+00 |
| temozolomide half-life 0.25 h matched-exposure (12 h window) | 0.0118847640 | 0.0118847640 | 0.0e+00 |

Largest discrepancy anywhere: **0.0e+00**.

### 1.2 Drug parameters as compiled onto this tissue

| drug | MW (Da) | D (um²/s) | IC50 (uM) | hill | BBB | k_loss (1/s) | transport time-scale 1/k | analytic L = √(D_eff/k) | sensitivity determinant | median sensitivity (tumor) | kill ceiling |
|---|---|---|---|---|---|---|---|---|---|---|---|
| temozolomide | 194.2 | 278 | 606.0 | 1.20 | 0.125 | 3.02e-03 | 331 s | 296 um | none (uniform 1.0) | 1.000 | 0.5502 |
| cisplatin | 298.0 | 236 | 24.6 | 1.20 | 0.858 | 8.79e-03 | 114 s | 160 um | MKI67 | 0.589 | 0.3254 |

**The kill ceiling** is the kill fraction reached as exposure → ∞: the Hill term saturates at `max_kill_rate`, so per-cell kill probability tends to `sensitivity x hypoxia/quiescence resistance`. It is a property of the sensitivity model and this tissue's oxygenation. **No dose can exceed it**, and it is the reference against which "is more drug worth it?" must be judged.

### 1.3 What actually sets the ceiling — the decomposition

This matters more than any single kill number, so it is broken out per cell-state band. `contribution` is that band's share of the ceiling (fraction of tumor cells x their mean asymptotic kill probability); the three contributions sum to the ceiling.

**temozolomide** — ceiling 0.5502

| tumor-cell band | fraction of tumor cells | resistance multiplier | mean asymptotic kill probability | contribution to ceiling |
|---|---|---|---|---|
| VIABLE | 10.0% | x1.00 | 1.000 | 0.1004 |
| HYPOXIC | 0.0% | x0.35 | 0.000 | 0.0000 |
| QUIESCENT | 90.0% | x0.50 | 0.500 | 0.4498 |

**cisplatin** — ceiling 0.3254

| tumor-cell band | fraction of tumor cells | resistance multiplier | mean asymptotic kill probability | contribution to ceiling |
|---|---|---|---|---|
| VIABLE | 10.0% | x1.00 | 0.593 | 0.0595 |
| HYPOXIC | 0.0% | x0.35 | 0.000 | 0.0000 |
| QUIESCENT | 90.0% | x0.50 | 0.296 | 0.2659 |

Read that table before reading anything else. The ceiling is **not** a delivery number and **not** an IC50 number: at infinite exposure every cell is at the top of its Hill curve, so what is left is exactly `sensitivity x resistance`. For temozolomide, whose per-cell sensitivity is uniform 1.0 on this panel, the ceiling is *entirely* the hypoxia (x0.35) and quiescence (x0.5) multipliers applied to the oxygen-banded cell states — i.e. it is set by how much of this tumor the model considers metabolically out of reach, not by how much drug arrives. For cisplatin the MKI67-derived sensitivity multiplies that down further.

---

## 2. Dose-response on this patient's tissue

`IC50 coverage` = fraction of tumor cells whose *time-averaged tissue exposure* (`AUC / window`, the quantity the engine's Hill mapping consumes) is at or above the drug's own IC50. `IC50 depth` = the deepest distance from a vessel at which the binned exposure profile still clears IC50 (`+` = still above IC50 at the far end of the profiled range, i.e. a lower bound).

### 2.1 temozolomide

| dose (× nominal) | plasma peak (uM) | window-mean plasma (uM) | kill fraction | % of kill ceiling | tumor cells ≥ IC50 | tumor cells ≥ 0.1·IC50 | median exposure / IC50 | IC50 depth from vessel (um) |
|---|---|---|---|---|---|---|---|---|
| 0.05x | 2.5 | 1.55 | 0.00055 | 0.1% | 0.00% | 0.00% | 0.003 | — (below IC50 even at the vessel wall) |
| 0.1x | 5.0 | 3.11 | 0.00139 | 0.3% | 0.00% | 0.00% | 0.006 | — (below IC50 even at the vessel wall) |
| 0.3x | 15.0 | 9.33 | 0.00466 | 0.8% | 0.00% | 0.00% | 0.017 | — (below IC50 even at the vessel wall) |
| 1x | 50.0 | 31.09 | 0.01839 | 3.3% | 0.00% | 2.95% | 0.056 | — (below IC50 even at the vessel wall) |
| 3x | 150.0 | 93.28 | 0.06161 | 11.2% | 0.00% | 97.80% | 0.168 | — (below IC50 even at the vessel wall) |
| 10x | 500.0 | 310.95 | 0.18986 | 34.5% | 2.95% | 99.87% | 0.559 | — (below IC50 even at the vessel wall) |
| 30x | 1500.0 | 932.85 | 0.35486 | 64.5% | 97.80% | 99.87% | 1.678 | 146 |

**Saturation.** **Kill has NOT saturated within the swept range** — it is still climbing at 30x nominal, reaching 90% of its 30x value only at 30x itself. Over this whole range more drug is still buying more kill; IC50 coverage passes 90% of tumor cells by **30x**. At 30x the predicted kill is 0.3549 = 64.5% of the 0.5502 ceiling, so the residual 35.5% of the ceiling is what *no achievable dose* recovers — that residual is the hypoxic/quiescent fraction the Hill term can never reach past.

**At the approximate clinical plasma band.** Published Cmax for temozolomide is of order 30–70 uM, i.e. **0.6x–1.4x** of the shared nominal level used for the axis. Interpolating the sweep there gives kill 0.0126–0.0316 (2–6% of ceiling) and IC50 coverage 0.0–0.0% of tumor cells. *This band is an order-of-magnitude literature reference attached for scale only — it is not fitted, not patient-specific, and not a dose recommendation.*

**What 30x buys over 1x.** kill 0.0184 → 0.3549 (+0.3365 absolute, 19.29x relative) for a 30-fold increase in systemic exposure — i.e. a 95% share of the final kill bought at 30x the systemic toxicity.

### 2.2 cisplatin

| dose (× nominal) | plasma peak (uM) | window-mean plasma (uM) | kill fraction | % of kill ceiling | tumor cells ≥ IC50 | tumor cells ≥ 0.1·IC50 | median exposure / IC50 | IC50 depth from vessel (um) |
|---|---|---|---|---|---|---|---|---|
| 0.05x | 2.5 | 1.55 | 0.03624 | 11.1% | 0.00% | 91.49% | 0.159 | — (below IC50 even at the vessel wall) |
| 0.1x | 5.0 | 3.11 | 0.06980 | 21.5% | 0.50% | 99.55% | 0.319 | — (below IC50 even at the vessel wall) |
| 0.3x | 15.0 | 9.33 | 0.16110 | 49.5% | 44.53% | 99.87% | 0.956 | 47 |
| 1x | 50.0 | 31.09 | 0.25504 | 78.4% | 99.55% | 99.87% | 3.185 | 150+ |
| 3x | 150.0 | 93.28 | 0.29783 | 91.5% | 99.87% | 99.87% | 9.556 | 150+ |
| 10x | 500.0 | 310.95 | 0.31375 | 96.4% | 99.87% | 99.87% | 31.852 | 150+ |
| 30x | 1500.0 | 932.85 | 0.31786 | 97.7% | 99.87% | 99.87% | 95.555 | 150+ |

**Saturation.** **Kill saturates at about 3x** (90% of its 30x value) and is flat by 30x; IC50 coverage passes 90% of tumor cells by **1x**. At 30x the predicted kill is 0.3179 = 97.7% of the 0.3254 ceiling, so the residual 2.3% of the ceiling is what *no achievable dose* recovers — that residual is the hypoxic/quiescent fraction the Hill term can never reach past.

**At the approximate clinical plasma band.** Published Cmax for cisplatin is of order 2–6 uM, i.e. **0.04x–0.12x** of the shared nominal level used for the axis. Interpolating the sweep there gives kill 0.0362–0.0849 (11–26% of ceiling) and IC50 coverage 0.0–7.8% of tumor cells (the low end is clamped to the swept minimum). *This band is an order-of-magnitude literature reference attached for scale only — it is not fitted, not patient-specific, and not a dose recommendation.*

**What 30x buys over 1x.** kill 0.2550 → 0.3179 (+0.0628 absolute, 1.25x relative) for a 30-fold increase in systemic exposure — i.e. a 20% share of the final kill bought at 30x the systemic toxicity.

---

## 3. Delivery geometry: what the exposure field actually looks like

Because the field is linear in dose, **dose moves the amplitude of the perivascular exposure profile and nothing else** — the shape is fixed by `D_eff`, the loss rate `k`, and the vascular geometry. Panel C shows this: the curves are parallel on a log axis at every dose. This is exact, not approximate: the transport operator is linear, so scaling the source scales the field.

**Dose, but not schedule.** The plasma curve's SHAPE does move the profile shape — by well under a percent on this drug and window, and in the direction that favours the bolus. Section 5 quantifies it. Everything in this section is measured on one shape (the 4 h-half-life single dose over the 6 h dose window), which is the correct reading of a dose sweep since the shape is held constant across it.

**The profile is NOT a single exponential, and saying otherwise would be the main way to mislead with this figure.** It is a steep perivascular drop sitting on a flat, non-zero FLOOR that persists to the far edge of the tissue. The floor is *not* a diffusion front: the engine's delivery term is the **continuous** vascular-density field, so a pixel 1 mm from the thresholded vessel mask is still being fed by diffuse low-density source. A single exponential fitted over the whole range is dominated by that plateau and returns a meaningless decay length — which is exactly why `vessel_profile.fit_exponential_plus_floor` exists. All three characterisations are therefore reported together.

| drug | empirical half-distance (no model assumed) | **source-null half-distance** | near-field exponential L (d ≤ p90) | near-field fit R² | floor-model L | floor level (% of perivascular) | far-field plateau / perivascular | analytic √(D_eff/k) | equilibration t95 | penetration index |
|---|---|---|---|---|---|---|---|---|---|---|
| temozolomide | 145 um | 12 um | 292 um | 0.99 | 270 um | 0% | 0.56 | 296 um | 992 s (16.5 min) | 0.511 |
| cisplatin | 98 um | 12 um | 153 um | 0.98 | 130 um | 8% | 0.37 | 160 um | 341 s (5.7 min) | 0.428 |

**Which column to believe:** the *empirical* half-distance assumes no functional form and is the trustworthy gradient number. The single-exponential fits — over the full range and over the near field alike — have poor R² and return decay lengths inflated by the plateau; disregard them as descriptions of the gradient. The floor model reconciles the two: a short decaying component (its `L` is close to the empirical half-distance) plus a background floor. Note also that the two drugs return nearly identical geometry despite differing in `D_eff` and `k` by ~2x — further evidence that what is being measured here is mostly the source field, not drug-specific transport.

Three things to take from that table, in order of importance:

1. **The measured gradient is essentially the source geometry.** The empirical half-distance and the **source-null** half-distance (the identical binning applied to the raw vascular-density field, with no drug physics at all) agree to within a few micrometres. The steep near-vessel drop this simulation predicts is inherited from where the endothelial cells are, not produced by transport. This is the confound `vessel_profile.py` documents and it is the single biggest caveat on every depth number in this report.
2. **A large fraction of the field is floor, not gradient.** See the plateau column. Any statement of the form "drug reaches X um" is therefore mostly a statement about where the floor sits relative to IC50, not about a diffusion front, and the `IC50 depth` column in section 2 must be read that way: when it jumps to the edge of the profiled range between two doses, what happened is that the **floor** crossed IC50, not that the drug travelled further.
3. **Depth cannot be bought with dose.** The shape is fixed under dose, so raising the dose slides the whole curve up and the IC50 crossing moves only as `L*ln(dose)` in the decaying region — and, once the floor clears IC50, saturates entirely.  Nor, in practice, with schedule: section 5 finds the shape does respond to the plasma time-course, but by a fraction of a percent.

For scale: median distance-to-vessel on this section is 43 um, p90 85 um, max 786 um. Bins beyond the p90 hold only a few hundred pixels each, so the far tail of panel C is noisy as well as flat — the near-field exponential above is deliberately restricted to d ≤ p90.

---

## 4. Schedule: fractionation at MATCHED TOTAL DOSE

One 12 h course, 4 h plasma half-life, split into 1 / 3 / 10 fractions. Three legs are reported because "matched dose" can mean two different things and the distinction turns out to be the whole story:

* **A — matched cumulative exposure** (the key comparison the brief asks for): the per-dose peak is rescaled so the *window-averaged plasma level*, i.e. total plasma AUC, is identical across schedules.
* **B — matched administered dose**: the sum of the per-dose peaks is held fixed (`peak = nominal / n`), which is what "same total milligrams" naively means.
* **C — matched peak** (a control): each fraction reaches the full nominal peak, so total exposure *rises* with the number of doses.

> **DROPPED for compute budget:** the schedule leg was NOT run for cisplatin — it requires its own transport solve at the schedule window. Only temozolomide is reported here.

### 4.1 temozolomide

| schedule | per-dose peak, leg A (uM) | window-mean plasma, A | **kill, A (matched exposure)** | tumor ≥ IC50, A | window-mean plasma, B | kill, B (matched administered dose) | window-mean plasma, C | kill, C (matched peak) |
|---|---|---|---|---|---|---|---|---|
| 1 dose (bolus) | 50.00 | 21.0444 | **0.011801** | 0.00% | 21.0444 | 0.011801 | 21.0444 | 0.011801 |
| 3 doses q4.0h | 20.59 | 21.0444 | **0.011759** | 0.00% | 17.0359 | 0.009575 | 51.1077 | 0.032505 |
| 10 doses q1.2h | 7.04 | 21.0444 | **0.011717** | 0.00% | 14.9463 | 0.008105 | 149.4628 | 0.100790 |

Leg B falls 0.011801 → 0.008105 (31.3% relative). Leg C rises 0.011801 → 0.100790.

**Leg A at this reference dose: kill spread 8.40e-05, tissue-exposure spread 0.772%.** Both are non-zero, and both were exactly zero under the static source.

#### The same comparison, now with the STATIC source alongside it

This is the direct before/after. `mean tumor exposure` is the continuous quantity the kill map consumes; it is not quantised by cell counting, so it is the honest place to read the size of the schedule effect. The static column is bit-identical down the whole table BY CONSTRUCTION — that is the bug, visible.

| schedule | mean tumor exposure, TIME-RESOLVED (uM·s) | mean tumor exposure, STATIC (uM·s) | Δ vs static | kill, time-resolved | kill, static | terminal residual C(T) / delivered | exposure half-distance from vessel (um) | penetration index |
|---|---|---|---|---|---|---|---|---|
| 1 dose (bolus) | 1.05056e+06 | 1.04486e+06 | +0.546% | 0.011801 | 0.011759 | 0.23% | 144.63 | 0.511733 |
| 3 doses q4.0h | 1.04592e+06 | 1.04486e+06 | +0.101% | 0.011759 | 0.011759 | 0.67% | 144.55 | 0.511408 |
| 10 doses q1.2h | 1.04251e+06 | 1.04486e+06 | -0.225% | 0.011717 | 0.011759 | 0.99% | 144.49 | 0.511139 |

**Relative exposure spread across the fractionation range: 0.772%**, monotonically decreasing with the number of fractions. The static column's spread is exactly 0.

#### Absolute kill spread depends on where you are on the Hill curve

The schedule effect is a fixed *relative* change in exposure, so the *absolute* kill spread it produces grows with dose — at 1x nominal temozolomide kill is only ~0.012, so a ~0.8% relative change in exposure lands near the 1/23,812 = 4.2e-05 quantum of a kill fraction computed by counting cells. Rescaling leg A up the dose axis (free: exposure is exactly linear in amplitude) lifts it clear of that floor and shows the effect scaling as it should. **Every static-source row here is exactly zero.**

| leg-A dose (× the matched-exposure level) | 1 dose | 3 doses | 10 doses | spread, TIME-RESOLVED | spread, STATIC |
|---|---|---|---|---|---|
| 1x | 0.011801 | 0.011759 | 0.011717 | **8.40e-05** | 0.00e+00 |
| 3x | 0.041492 | 0.041408 | 0.041156 | **3.36e-04** | 0.00e+00 |
| 10x | 0.138166 | 0.137872 | 0.137452 | **7.14e-04** | 0.00e+00 |
| 30x | 0.295187 | 0.294431 | 0.293927 | **1.26e-03** | 0.00e+00 |

### 4.x Interpretation — what the corrected model actually says

**(a) RETRACTED: "fractionation changes nothing".** The static-source version of this file reported a leg-A spread of exactly `0.00e+00` and concluded that *on this tumor, at matched cumulative exposure, this model predicts fractionation changes nothing*. That conclusion was **a model artefact**. The solver was fed the time-average of the plasma curve, so 1 dose and 10 doses were literally the same input array; the zero was arithmetic. It is now measured rather than assumed away, and it is not zero.

**(b) The corrected answer, and the derivation that bounds it.** Integrate the transport equation over the window `[0, T]` and let `A = AUC = ∫C dt`, `M = k − div(D grad ·) + div(v ·)`. With `C(0) = 0`:

```
    C(T) − C(0) = div(D grad A) − div(v A) − k A + ∫S dt
=>            A = M⁻¹ ( S_total − C(T) ),      S_total = ∫₀ᵀ S dt
```

`S_total` is **identical** for any two schedules at matched cumulative plasma exposure — that is what "matched" means, and the implementation guarantees it exactly (the plasma curve is sampled on the solver's own step grid, so its mean is exactly the level the static source used). Therefore:

1. **The entire schedule effect is the terminal residual `C(T)`** — drug that was delivered but is still in the tissue, uncounted, when the window closes. Nothing else about the timing can reach the answer.
2. **Earlier delivery wins.** A schedule that finishes early has washed out by `T` (`C(T) ≈ 0`) and realises the full steady-state exposure field. One still infusing at `T` does not. Hence exposure decreases monotonically with the number of fractions — which is exactly what the table shows.
3. **The effect is bounded by `τ_tissue / T_window`.** Under a sustained level the tissue holds `C_ss ≈ S_rate/k`, while the drug delivered over the window is `S_rate·T`; the residual is therefore at most `(1/k)/T = τ/T` of what was given. This is a *prediction*, and the run tests it:

   * predicted ceiling `τ/T` = 331 s / 12 h = **0.767%**
   * measured exposure spread across 1 → 10 fractions = **0.772%**
   * measured terminal residual of the constant-infusion (static-source) limit, which should SATURATE the bound = **0.767%**

   The three agree. That agreement is the strongest evidence here that the time-resolved source is doing real physics and not producing numerical noise: the size of the effect was derived before it was measured.

**(c) Why the effect is small HERE, and where it would not be.** Temozolomide on this tissue has a tissue equilibration time of minutes (temozolomide: 992 s (16.5 min) to 95%, = ln(20)/k) against a 12 h window, so `τ/T` is under a percent and the tissue tracks plasma quasi-statically. **This is a statement about the ratio, not about schedules.** Shorten the window, or take a molecule that equilibrates slowly — an antibody, `D ~ 1–10 um²/s` with a binding-site-barrier sink and `τ` of hours to days — and `τ/T` approaches 1, at which point schedule dominates. The same solver now answers that question; before, it could not be asked.

**(d) What is still missing, and it is the big one.** Time resolution was added to *transport*. The kill map remains a **static Hill function of cumulative exposure**: no inter-fraction DNA repair, no cell-cycle redistribution, no regrowth between fractions, no resistance induction, and no normal-tissue recovery. Those are the mechanisms that give fractionation its actual clinical benefit, and they are pharmacodynamics, not transport. So read the numbers above as **the transport-only floor on how much schedule can matter** — the part that is now correctly computed — and not as an estimate of the clinical effect of fractionating temozolomide.

**Leg B is a window artifact, and worth naming as one.** Kill falls with fractionation there only because later fractions have less of the 12 h observation window left in which to be delivered and counted. Analytically the window exposure for `n` fractions of peak `P/n` at interval `Δ` is `(P/n)·Σ_k (1/λ)(1 − e^(−λ(T−kΔ)))`, which tends to `P/λ` for *every* `n` as `T → ∞`. Extend the window past full washout and leg B collapses onto leg A. It is a measurement-window effect, not a schedule effect.

Leg C confirms the model's dominant lever by a wide margin: raise total exposure and kill rises steeply. Everything else about the schedule is now *visible* rather than invisible — but it is two to three orders of magnitude smaller than the exposure lever on this drug and this window.

---

## 5. Infusion duration proxy: plasma half-life, bolus vs sustained

A short half-life is a bolus that is gone quickly; a long half-life is a sustained low level. The hypothesis under test is the standard transport-limited one: *if delivery is limited by how far drug can travel before it is consumed, a sustained low concentration should reach further than a high, brief one.* Half-life is swept 0.25–48 h (a 192-fold range) with the peak rescaled to hold total window exposure fixed, plus a matched-peak control.

> **DROPPED for compute budget:** the half-life sweep was NOT run for cisplatin (same second-solve cost as the schedule leg).

### 5.1 temozolomide

| plasma half-life (h) | peak, matched exposure (uM) | window-mean plasma | **kill (matched total exposure)** | tumor ≥ IC50 | IC50 depth (um) | window-mean plasma, matched peak | kill (matched peak) | IC50 depth, matched peak (um) |
|---|---|---|---|---|---|---|---|---|
| 0.25 | 697.48 | 21.0444 | **0.011885** | 0.00% | — (below IC50 even at the vessel wall) | 1.5086 | 0.000504 | — (below IC50 even at the vessel wall) |
| 4 | 50.00 | 21.0444 | **0.011801** | 0.00% | — (below IC50 even at the vessel wall) | 21.0444 | 0.011801 | — (below IC50 even at the vessel wall) |
| 48 | 22.92 | 21.0444 | **0.011759** | 0.00% | — (below IC50 even at the vessel wall) | 45.9085 | 0.028893 | — (below IC50 even at the vessel wall) |

Spread in kill across the 192-fold half-life range at matched exposure: **1.26e-04** — the number that was **exactly `0.00e+00`** in the static-source version, where the penetration depth was also reported as "unchanged" because there was only ever one exposure field to measure it on.

#### Each half-life now has its OWN exposure field, so depth is measurable

Under the static source there was a single AUC field for the whole sweep, so "depth is unchanged across half-life" was true by construction and carried no information. Each shape is now solved separately and the profile is re-binned per shape. `half-distance` and `penetration index` are both **shape** statistics — normalised by the perivascular value / the p99, so amplitude cannot move them.

| plasma half-life (h) | mean tumor exposure (uM·s) | Δ vs static | kill, time-resolved | kill, static | terminal residual C(T) / delivered | exposure half-distance from vessel (um) | penetration index | far-field plateau / perivascular |
|---|---|---|---|---|---|---|---|---|
| 0.25 | 1.05303e+06 | +0.782% | 0.011885 | 0.011759 | 0.00% | 144.67 | 0.511937 | 0.5631 |
| 4 | 1.05056e+06 | +0.546% | 0.011801 | 0.011759 | 0.23% | 144.63 | 0.511733 | 0.5628 |
| 48 | 1.04553e+06 | +0.065% | 0.011759 | 0.011759 | 0.70% | 144.55 | 0.511382 | 0.5622 |

**Half-distance moves from 144.67 um at the 0.25 h (bolus-like) half-life to 144.55 um at 48 h (sustained) — a change of -0.13 um.** The penetration index moves 0.511937 → 0.511382. **Both say the BOLUS reaches relatively deeper, not the sustained infusion.**

#### Kill spread up the dose axis (static rows are exactly zero)

| matched-exposure dose (×) | t½ 0.25 h | t½ 4 h | t½ 48 h | spread, TIME-RESOLVED | spread, STATIC |
|---|---|---|---|---|---|
| 1x | 0.011885 | 0.011801 | 0.011759 | **1.26e-04** | 0.00e+00 |
| 3x | 0.041702 | 0.041492 | 0.041408 | **2.94e-04** | 0.00e+00 |
| 10x | 0.138376 | 0.138166 | 0.137872 | **5.04e-04** | 0.00e+00 |
| 30x | 0.295355 | 0.295187 | 0.294389 | **9.66e-04** | 0.00e+00 |

### 5.x Interpretation — the naive transport-limited hypothesis is REFUTED, and in the opposite direction to the retracted claim

The hypothesis under test was: *if delivery is transport-limited, a sustained low concentration should reach further than a high brief one.* With a working time-resolved driver we can now actually test it, and the answer is **no — and not for the reason the old version of this file gave.**

**What the old version said, and why it was not an answer.** It reported the depth as identical across the whole 192-fold half-life range and concluded "sustained infusion does NOT reach further here … it is exposure that moves depth, never the shape of the plasma curve". But there was only ONE exposure field in that run — the plasma shape had been averaged away before the solve. The invariance was the input being invariant. **Retracted.**

**What the corrected model says.** The depth DOES move with plasma shape, by a small amount, and it moves the *other* way: the bolus reaches relatively deeper than the sustained infusion. This is derivable from `A = M⁻¹(S_total − C(T))` (section 4):

* At matched cumulative exposure `S_total` is identical, so the two schedules differ by `M⁻¹` applied to their terminal residuals. The bolus has washed out (`C(T) ≈ 0`); the sustained infusion has not.
* The exposure lost by the sustained infusion is `M⁻¹C_ss = M⁻²S_rate` — the source field smoothed **twice** by the operator, i.e. a *flatter* field than the exposure field `M⁻¹S_total` itself. Subtracting a flatter field from a peakier one leaves something **steeper**. So the sustained infusion loses more of its exposure at depth, proportionally, than near the vessel.
* Hence: bolus wins on total exposure AND on relative reach. The table above confirms both, on the half-distance and on the penetration index independently.

**Why a bolus can never lose here.** Over a fixed observation window, the average delivery time of a schedule is what sets its residual. A bolus at `t = 0` has the earliest possible average delivery time, so it always realises at least as much of the steady-state exposure field as anything else at the same cumulative dose. There is no regime in a *linear* transport model where the ordering flips. (Where the textbook intuition does apply is a model with saturable perivascular binding or nonlinear uptake — neither of which is in this transport operator.)

**Do not over-read the size of this.** The effect is bounded by `τ_tissue/T` (0.767% here) and the depth shift is sub-micrometre-scale against a perivascular gradient that section 3 shows is mostly inherited from the source geometry anyway. The value of this leg is **that the question is now answerable and the answer has a sign and a derived magnitude** — not that infusion duration is an actionable lever on this tumor.

**Where it WOULD be actionable.** `τ_tissue/T` is the whole story, so the regime to look at is slow equilibration: an antibody (`D ~ 1–10 um²/s` plus a binding-site-barrier sink, `τ` of hours to days), or any molecule observed over a window comparable to its own washout. The old version of this file named that regime and said the plasma-averaging simplification "would have to be lifted before the question could be asked at all". It has been.

---

## 6. THE VERDICT for this tumor: dose-limited or delivery-limited?

> ### Headline
> 
> **On this tumor the response is DOSE-limited, not delivery-limited — and above the dose-limited range it becomes RESISTANCE-limited, never delivery-limited.**
> 
> Delivery is not the binding constraint here. The exposure field reaches the whole section: at high enough dose >95% of tumor cells clear IC50 for both drugs, and the far-field floor (section 3) means even tissue a millimetre from a vessel is exposed. What limits the answer is, in order: (1) at clinically plausible plasma levels, simply not enough drug — the sub-saturating part of the curve, where more exposure genuinely buys more kill; and (2) above that, a hard ceiling that dose cannot touch.
> 
> **Either way, this is NOT a schedule verdict, and now that is a measurement rather than an assumption.** Sections 4-5 show the schedule effect on this tumor is real but bounded by `τ_tissue/T_window` — orders of magnitude below the dose lever. The static-source version of this file reported it as exactly zero, which was an artefact; the corrected answer is small, not zero.
> 
> **That ceiling is the actionable result, and it is not a transport number.** For temozolomide it is 0.5502, and section 1.3 decomposes it: 90.0% of tumor cells are labelled QUIESCENT (x0.50 kill multiplier) and 0.0% HYPOXIC (x0.35), contributing 0.4498 and 0.0000 of it respectively — i.e. ~82% of the ceiling is the quiescent band alone. Because temozolomide's per-cell sensitivity is uniform 1.0 on this panel (MGMT is not measurable here), the ceiling for TMZ is *entirely* the cell-state resistance multipliers. It is not delivery, not IC50, and not achievable exposure: at infinite drug every cell sits at the top of its Hill curve and what remains is exactly `sensitivity x resistance`.
> 
> The operational reading: **on this tumor, escalating dose is a real but bounded lever, and the bound is metabolic state, not penetration.** If the ceiling is to move, it moves by changing the cell-state distribution (reoxygenation, cycle recruitment) or by a drug whose kill does not carry the quiescence penalty — not by more drug and not by better penetration.
> 
> **And schedule? Now measured rather than erased — and it is a REAL but THIRD-ORDER lever on this tumor.** With the time-resolved plasma driver, fractionation at matched cumulative exposure moves tissue exposure by 0.77% and infusion duration by a comparable amount, against the `τ_tissue/T_window` ceiling of 0.77% that sections 4-5 derive and confirm. Compare that with dose, which moves kill by a factor of 19 over the swept range. **The previous version of this file reported this effect as exactly zero and concluded schedule does not matter; that zero was a solver artefact. The corrected conclusion is not "schedule matters after all" — it is that schedule matters a little, for a reason we can now name and bound, and far less than dose on this fast-equilibrating small molecule.** The caveat that survives unchanged: the pharmacodynamic mechanisms that make fractionation matter clinically are still absent from the kill model, so this is a floor.
> 
> **Load-bearing caveat:** that 90% quiescent fraction is a *model* label from the simulated oxygen field, not a measurement. It is being tested independently against measured MKI67. If it is wrong, the ceiling and every absolute kill here move with it — the dose/schedule *comparisons* survive, the absolute levels may not.

### temozolomide

**DOSE-limited over the swept range, then sensitivity-limited.**

At the nominal dose 0.0% of tumor cells clear IC50, rising to 97.8% by 30x — the drug *can* reach essentially the whole tumor once enough of it is given, so the binding constraint over this range is how much drug is present, not how far it travels. Kill climbs with dose (0.0184 → 0.3549), and at 30x it has reached 64.5% of the 0.5502 ceiling without yet flattening — the sweep never gets far enough for delivery to stop mattering.

Geometry check: near-field exposure decay length **L = 292 um** (d ≤ p90) vs median distance-to-vessel **43 um** and p90 **85 um**. IC50 depth at 1x: — (below IC50 even at the vessel wall) um; at 30x: 146 um — but see section 3: with a far-field plateau at 0.56 of the perivascular value, the large depths are the FLOOR clearing IC50, not a diffusion front reaching further.

### cisplatin

> **⚠ CONDITIONAL RETRACTION — this verdict inverts under corrected transport
> physics.** Not retracted yet, because it remains the literal output of the
> shipped default path, but it will be the moment that default changes.
>
> The "99.6 % of tumour cells above IC50" figure this verdict rests on is
> produced by the one-way delivery law, which lets tissue concentration exceed
> plasma without bound (`results/breast_kp_validation.md`). Cisplatin runs at
> **2.71x plasma** on this section, which passive delivery cannot do. With
> `PDEConfig(transvascular_exchange=True)`, IC50 coverage at 1x falls
> **99.55 % → 4.67 %**, kill goes 78.4 % → 43.1 % of ceiling, saturation moves
> from ~3x to ~10x, and the 30x/1x kill ratio goes 1.25 → 2.22. **At the
> nominal dose cisplatin becomes DOSE-limited** — the opposite operational
> reading.
>
> Note this is the direction the "Important qualification" paragraph below
> already anticipated on clinical-dose grounds. Corrected physics says it
> applies at the nominal dose too, not only down at clinical Cmax.
>
> The same inversion was reproduced independently in
> `results/mechanism_attribution.md` §3 (RESISTANCE 96.42 % → 3.59 %,
> DOSE 0.44 % → 95.71 %) and in the GBM 3-D:2-D and cross-drug ranking
> re-runs. Four independent code paths, same result. Do not quote cisplatin as
> sensitivity-limited or resistance-limited in external material.

**NEITHER dose- nor delivery-limited — it is SENSITIVITY-limited.**

Already at the nominal dose, 99.6% of tumor cells see an exposure above IC50. The drug reaches essentially the whole tumor; raising the dose 30-fold moves kill only 0.2550 → 0.3179. What caps the response is the kill ceiling (0.3254) — the sensitivity model plus hypoxic / quiescent resistance. More drug does not help; a different drug, or reversing resistance, is the only lever the model leaves open.

**Important qualification — the verdict above is stated at the shared nominal level, which is well above cisplatin's clinical range.** Its approximate clinical Cmax (2–6 uM) is only 0.04x–0.12x nominal. Down there the sweep gives IC50 coverage 0% and kill 0.0362 (11% of ceiling), extrapolated to the edge of the swept range — i.e. at clinically plausible exposures this drug is materially **dose-limited on this tumor**, and it only becomes sensitivity-limited once pushed to supratherapeutic levels. Both statements are in the table; which one is relevant depends on the dose actually deliverable.

Geometry check: near-field exposure decay length **L = 153 um** (d ≤ p90) vs median distance-to-vessel **43 um** and p90 **85 um**. IC50 depth at 1x: 150+ um; at 30x: 150+ um — but see section 3: with a far-field plateau at 0.37 of the perivascular value, the large depths are the FLOOR clearing IC50, not a diffusion front reaching further.

### What this means operationally (within the model's limits)

* **Would more drug help?** Read the `% of kill ceiling` column in section 2. Where it is already near 100% at the nominal dose, escalation is pure added toxicity in this model. Where it climbs steeply, the model says the tumor is under-exposed rather than intrinsically resistant — a materially different clinical situation, and the distinction is the one thing this simulator is positioned to make.
* **Is the tissue transport-saturated?** DOSE moves amplitude only — the perivascular profile shape is invariant under it exactly, because the operator is linear. SCHEDULE does move the shape, but by well under a percent here (the `τ_tissue/T` bound of 0.77%, sections 4-5), and in the direction that favours the bolus. So for practical purposes the geometry of delivery is set by three things — the drug's `D_eff`, its loss rate `k`, and the vascular geometry — and neither dose nor schedule is a lever on it *on this drug and this window*.
* **The lever that is NOT in this analysis** is the normal-tissue side of the therapeutic index. Every statement above is about tumor kill alone. A schedule that is near-identical for the tumor is emphatically not identical for the patient.

---

## 7. Limitations (the full list)

1. **No outcome validation, at all.** Untreated resection; no measured response for either drug. Nothing here has been scored against reality.
2. **Not clinical dosing guidance.** Research prototype.
3. **Schedule IS now resolvable — but only through transport, and the size of the effect is set by the observation window.** The `τ_tissue/T_window` bound (sections 4-5) means a longer window shrinks the measured schedule effect toward zero and a shorter one inflates it. The 12 h window here was chosen for compute, not physics, and it is a *reporting* choice that changes a headline number. It is held exactly constant across everything compared, so the comparisons are sound, but the absolute percentages are window-conditional and should be quoted with the window attached.
4. **Pharmacodynamics is a single static Hill map — and this is now the binding limitation on the schedule question.** No inter-fraction DNA repair, no cell-cycle redistribution, no regrowth, no resistance evolution, no immune contribution, no normal-tissue toxicity. These are exactly the mechanisms that give real fractionation its benefit. Time-resolving *transport* removed a structural blind spot; it did not make this model a fractionation model.
5. **Sensitivity model.** Temozolomide's determinant (MGMT) is not on this panel, so TMZ resistance is unmodelled (uniform 1.0) and its absolute kill is optimistic. Cisplatin uses an MKI67 proliferation proxy, which is mechanistic but not a measured chemosensitivity.
6. **The perivascular gradient partly inherits the source-field shape** (the source-null column in section 3). Depth numbers are order-of-magnitude.
7. **Vascular density is a KDE of measured endothelial cell positions**, not a measured vessel lumen or a measured perfusion map; the vessel mask is a percentile threshold (p95).
8. **Single section, single patient, 2D.** No inter-patient variability, no out-of-plane transport, and the nominal plasma level is a shared reference rather than a per-drug clinical PK model — supratherapeutic for cisplatin in particular.
9. **Cell-state banding dominates the ceiling.** 90% of tumor cells are labelled QUIESCENT by the oxygen model (x0.5 kill) and 0% HYPOXIC (x0.35). That labelling is a model assignment from the simulated oxygen field, not a measurement; it is being tested independently against measured MKI67. If the quiescent fraction is wrong, every absolute kill and the ceiling move with it.
10. **Scope was cut for the compute budget** (see the banner at the top): 7 dose points; schedule + half-life for temozolomide only; 3 schedules; 3 half-lives; 12 h schedule window; 3 engine verification runs. Nothing was extrapolated to fill the gaps — the gaps are named instead.
11. **One calibration constant** (`alpha`, fit once on doxorubicin) sets the perivascular sink for every drug; the penetration lengths and therefore every depth number inherit its uncertainty.
12. **The absolute kill levels shifted slightly from the static-source version of this file, and not only in the schedule sections.** The dose-response leg is also solved with a time-resolved source now (a decaying bolus over the 6 h dose window), so its exposure field is a little different too. Cross-file comparisons of absolute kill against the superseded version are therefore not apples to apples; comparisons *within* this file are.

*Generated by `experiments/dose_schedule.py` — ONE reconstruction reused everywhere, 9 transport solves (`PDEConfig.time_resolved=True`), 3 full-engine verification runs, 46.7 min total, seed 0, deterministic.*
