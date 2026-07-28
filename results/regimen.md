# REAL TEMOZOLOMIDE REGIMENS, AND COMBINATION SCHEDULING, ON ONE REAL TUMOR

> ## WHY THIS FILE EXISTS, AND WHY IT COULD NOT HAVE EXISTED YESTERDAY
>
> Until the time-resolved plasma driver landed, `sim.engine` collapsed the whole multi-dose plasma profile to a single **time-averaged** level before the transport solver ever saw it. Two regimens with the same cumulative dose were therefore **bit-identical inputs**: the measured fractionation spread was exactly `0.00e+00`, an arithmetic artefact and not a result about any tumor. Every schedule question, and every combination-scheduling question, was **structurally unanswerable** — and the combination and regimen work done before today was comparing identical inputs.
>
> **Every solve in this report uses `PDEConfig(time_resolved=True)`.** There is no static-source path anywhere in `experiments/regimen.py`. The source is rescaled by `plasma(t)` at every timestep, so tissue concentration genuinely rises and falls and each schedule shape gets its own transport solve.
>
> **The headline, in one line.** At matched cumulative dose over a full 28-day cycle, this model predicts that the real temozolomide regimens are indistinguishable: the total spread in mean tumour exposure across all 10 arms is **0.000359%**, against the analytic bound `tau_tissue / T_window` = **0.013692%**. That is not a claim that schedule does not matter clinically. It is a claim that **transport cannot be the reason it does** — and the mechanisms that make fractionation matter (MGMT depletion, inter-fraction repair, cell-cycle redistribution, immune effects, normal-tissue recovery) are all ABSENT from this model. Section 7 says so at length.

> **Configuration.** `vascular_mode = discrete` (promoted in [`results/vessel_promotion.md`](vessel_promotion.md)) · **alpha = 148.6** (single perivascular-sink constant, fit ONCE on doxorubicin) · per-cell sensitivity mode = **`uniform`** (the pipeline default and the mode `results/combination.md` used — it makes the cross-drug comparison fair by making every agent equally an upper bound; the engine still applies the hypoxia/quiescence resistance factor) · **dx = 12 um** · `dt = 10 s` · convection on · seed 0, fully deterministic.
>
> Absolute exposures are on the **discrete-source unit** and are not comparable with KDE-era numbers elsewhere in the repo. Absolute kills over a 28-day window are also not comparable with the 6 h-window numbers in `results/dose_schedule.md` or `results/combination.md`, because the engine's Hill map consumes `AUC / T_window` — see the limitation on this in section 8. Comparisons *within* this file are apples to apples.

> ## READ THIS FIRST
>
> * **Research prototype. NOT clinical dosing guidance.** Nothing here should influence the care of any patient. No number below is a dose recommendation, and the real regimens are modelled to make a *model* comparison interpretable, not to evaluate them.
> * **Zero outcome validation.** This 10x Xenium section is an **untreated** resection. There is no measured drug response for any agent on this tissue, so every kill fraction is a prediction with nothing to score it against.
> * **Absolute kill on this section is NOT defensible, and that is not my hedge — it is this repo's own finding.** `results/cell_state_impact.md` concludes verbatim: *"on the 319-gene flagship section, stop quoting an absolute kill fraction or a kill ceiling as a result at all: NO route clears the bar there, and the honest output of that section is a ranking, a coverage map and a delivery metric"*, and records that absolute kill moves by up to 0.3252 absolute (1.77x relative) purely from the choice of cell-state route, with tissue, drug, transport and RNG all held identical. **This is that section.** The state route here is the shipping default, `state_mode = 'oxygen'`. Every absolute kill and the kill ceiling below therefore carry an unresolved uncertainty of that size and are printed ONLY because the comparisons are computed from them. The comparisons — which are what this file is for — survive; the levels do not.
> * **~90% of tumor cells are labelled QUIESCENT by the oxygen-banding model** (x0.5 kill multiplier). That single *model* labelling dominates the kill ceiling of 0.5502 and every absolute kill here — and per `results/cell_state_impact.md` the oxygen route's VIABLE band on this section moves from 0.0000% to 10.0412%, i.e. empties entirely, purely with bookkeeping.
> * **Time resolution was added to TRANSPORT, not to pharmacodynamics.** The kill map is still a static Hill function of window-averaged exposure. The schedule effect reported here is a **transport-only floor**.
> * **Bliss independence is an ASSUMPTION, not a result.** The combinations model SPATIAL complementarity only. No biochemical synergy or antagonism is or can be represented.

---

## 1. Method

* **Tissue — one real tumor, reconstructed ONCE.** The open 10x Xenium GBM section: 40,887 measured cells (23,812 tumor), grid (424, 675) at the grid-converged `dx = 12 um`, mean vascular density 0.027, mean ECM 0.076, mean oxygen 0.352. Reconstruction took 11.1 s and the drug-independent prologue (oxygen -> cell-state banding -> IFP -> Darcy velocity) 12.7 s. **Both were run once and reused by every arm**, so all differences below are drug/schedule and not tissue. That reuse is the whole reason a 28-day cycle is affordable at all; the cost it avoids is one full reconstruct + prologue (~24 s) per arm across 20 regimen arms plus every combination arm.
* **Drug compilation.** `calibrated_drugparams(name, alpha = 148.61)` then `apply_barrier_to_drugparams` — the identical pipeline `real_drug_panel.py` and `dose_schedule.py` use.
* **Numerics.** `dt = 10 s`, convection on. `dx = 10 um` is deliberately not used — it violates the upwind advection CFL condition at this `dt`.
* **Total cost.** 21 transport solves, 136,080 solver steps, 109.6 min wall clock.

### 1.1 One solve per single dose, then EXACT superposition

A 28-day cycle at `dt = 10 s` is 241,920 steps per shape, and there are 10 regimen shapes. Solving each directly is not affordable. It is also not necessary, and the reason is exact rather than approximate.

The discrete solve is **linear** in the plasma profile (the operator `A = I - dt*L + dt*K` does not involve the source, the scheme starts from `C = 0`, and the source enters as `dt * plasma(t_n) * u(x,y)`) and it is **time-invariant** (`D_eff`, `k` and `v` are static). `plasma_profile` builds a multi-dose curve as a superposition of shifted single-dose exponentials. Therefore, for doses of peak `P_k` at times `t_k` observed over `[0, T]`:

```
    AUC_field  =  sum_k  P_k * A1(T - t_k)
    C(T)       =  sum_k  P_k * Cend1(T - t_k)
```

where `A1(s)` is the exposure field from ONE unit-peak dose observed over a window of length `s`. Every dose time here is an exact multiple of `dt`, so this identity holds **exactly on the discrete grid** — it is not interpolated and not extrapolated. Amplitude is likewise exactly free, which is the same trick `dose_schedule.py` uses for its dose axis.

`A1(s)` saturates once the plasma curve has washed out AND the tissue has cleared, so a 28-dose regimen is built from a handful of short solves. **The saturation is certified, not assumed**, with an exact bound: a unit dose eventually deposits all its drug, so

```
    (A1(inf) - A1(S)) / A1(inf)  =  2^(-S/t_half)  +  C(S) / (k * A1(S))
                                    \_ plasma  _/     \_ tissue residual _/
```

and BOTH terms are computed from the solved fields. The plasma term alone would NOT be sufficient — a molecule whose plasma has washed out may still be equilibrating in tissue, and that second term is exactly the terminal residual through which schedule acts in the first place.

| single-dose kernel | plasma t1/2 (h) | windows solved (h) | saturation window S (h) | un-delivered plasma at S | tissue still holding at S | certified error | substitution USED? |
|---|---:|---|---:|---:|---:|---:|---|
| temozolomide | 1.80 | 6, 12, 18, 24, 48 | 48 | 9.39e-09 | 3.45e-10 | **9.73e-09** | **YES** (28-day regimen leg) |
| temozolomide | 4.00 | 6, 12, 24 | 24 | 1.56e-02 | 2.57e-04 | 1.59e-02 | no — exact grid only |
| cisplatin | 4.00 | 6, 12, 24 | 24 | 1.56e-02 | 8.75e-05 | 1.57e-02 | no — exact grid only |
| doxorubicin | 4.00 | 6, 12, 24 | 24 | 1.56e-02 | 2.30e-05 | 1.56e-02 | no — exact grid only |
| cisplatin | 0.50 | 12, 24 | 24 | 3.55e-15 | 1.63e-16 | 3.72e-15 | no — exact grid only |
| doxorubicin | 30.00 | 12, 24 | 24 | 5.74e-01 | 2.60e-04 | 5.75e-01 | no — exact grid only |

Read the last column before the second-to-last. Only the temozolomide kernel that drives the 28-day regimen leg ever *uses* the saturated field, and its certified error is the one that matters. Every other kernel is read at exact grid points, so its 'certified error' column is informational — the doxorubicin 30 h row in particular shows a large number precisely because a 30 h half-life has NOT washed out by 24 h, which is exactly why the code refuses to substitute there and solves each window it needs.

The substitution is only ever *used* for the temozolomide kernel driving the 28-day regimen leg, where the certified error is **9.7e-09**. Every other leg reads exact grid entries; the kernel refuses to interpolate and raises instead.

### 1.2 Verification — three independent checks, each a full extra solve

The superposition is not taken on faith.

| # | check | statistic | result |
|---|---|---|---|
| V1 | direct 2-dose q12h solve vs A1(24h)+A1(12h), 24 h window | relative L1 over tissue | **1.74e-15** |
| V2 | FULL sim.engine.simulate_response (time_resolved=True), temozolomide 2 doses q12h, 24 h window | kill fraction, engine vs fast path | engine 0.0051234672 vs 0.0051234672, delta **0.0e+00** |
| V2b | (same run) | exposure field, relative L1 | **2.11e-15** |
| V3 | doxorubicin dose DELAYED 12 h in a 24 h solve vs A1(12h) | relative L1 over tissue | **0.00e+00** |

V1 tests LINEARITY (a genuine two-dose profile against the sum of two single-dose kernels). V3 tests TIME-INVARIANCE on its own, by driving the solver with a manually shifted profile — zeros for 12 h, then the dose — and comparing against the un-shifted kernel over the remaining window. V2 puts the **unmodified** `sim.engine.simulate_response` with `time_resolved=True` against the fast path end to end, so nothing about the orchestration here is unpinned.

### 1.3 Drug parameters as compiled onto this tissue

| drug | MW (Da) | D (um^2/s) | IC50 (uM) | hill | BBB | k_loss (1/s) | tissue time-scale tau = 1/k |
|---|---:|---:|---:|---:|---:|---:|---:|
| temozolomide | 194.2 | 278 | 606.0 | 1.2 | 0.125 | 3.02e-03 | 331 s |
| cisplatin | 298.0 | 236 | 24.6 | 1.2 | 0.858 | 8.79e-03 | 114 s |
| doxorubicin | 543.5 | 185 | 0.5 | 1.2 | 0.021 | 3.33e-02 | 30 s |

**Kill ceiling on this tumor = 0.5502** — the kill fraction reached as exposure -> infinity, where the Hill term saturates and per-cell kill probability tends to `sensitivity x resistance`. **No dose and no schedule can exceed it.** It is quoted here only because the `% of ceiling` column in section 3 is computed from it; per `results/cell_state_impact.md` a kill ceiling on THIS section is explicitly not a defensible result, and the decomposition below should be read as "what the model's bookkeeping currently implies", not as a property of the tumor.

| tumor-cell band | fraction of tumor cells | mean asymptotic kill probability | contribution to ceiling |
|---|---:|---:|---:|
| VIABLE | 10.0% | 1.000 | 0.1004 |
| HYPOXIC | 0.0% | 0.000 | 0.0000 |
| QUIESCENT | 90.0% | 0.500 | 0.4498 |

Temozolomide's per-cell sensitivity determinant (MGMT) is **not on this panel**, so under `uniform` sensitivity the ceiling is *entirely* the hypoxia/quiescence multipliers. Temozolomide resistance is therefore **unmodelled** and its absolute kill is correspondingly optimistic.

---

## 2. The real regimens, and the PK conversion

### 2.1 What is actually given to patients

Every regimen below is a real, published schedule. The 28-day comparison window is one adjuvant cycle; the Stupp concomitant phase runs 42 days in practice and only its first 28 days are modelled, so its cumulative figure here is per 28 days and not per course.

| regimen | dose | dosing days of the 28-day window | doses | cumulative mg/m2 per 28 d | vs Stupp adjuvant C1 | source |
|---|---|---|---:|---:|---:|---|
| Stupp adjuvant, cycle 1 | 150 mg/m2/d | 1-5 | 5 | **750** | 1.00x | Stupp 2005 (PMID 15758009) |
| Stupp adjuvant, cycles 2-6 | 200 mg/m2/d | 1-5 | 5 | **1000** | 1.33x | Stupp 2005 (PMID 15758009) |
| Stupp concomitant (with RT) | 75 mg/m2/d | 1-28 (continuous) | 28 | **2100** | 2.80x | Stupp 2005 (PMID 15758009) |
| Dose-dense 21/28, cycle 1 (RTOG 0525) | 75 mg/m2/d | 1-21 | 21 | **1575** | 2.10x | Gilbert 2013, RTOG 0525 (PMID 24101040) |
| Dose-dense 21/28, cycles 2+ (RTOG 0525) | 100 mg/m2/d | 1-21 | 21 | **2100** | 2.80x | Gilbert 2013, RTOG 0525 (PMID 24101040) |
| Dose-intense 7 on / 7 off (Wick) | 150 mg/m2/d | 1-7, 15-21 | 14 | **2100** | 2.80x | Wick 2004 (PMID 15184628); Wick 2007 (PMID 17664483) |
| DIRECTOR arm A (7/7) | 120 mg/m2/d | 1-7, 15-21 | 14 | **1680** | 2.24x | Weller 2015, DIRECTOR (PMID 25655102) |
| DIRECTOR arm B (21/28) | 80 mg/m2/d | 1-21 | 21 | **1680** | 2.24x | Weller 2015, DIRECTOR (PMID 25655102) |
| Continuous low-dose (RESCUE) | 50 mg/m2/d | 1-28 (continuous) | 28 | **1400** | 1.87x | Perry 2010, RESCUE (PMID 20308655) |
| Single large dose (CONTROL) | 1000 mg/m2, once | 1 | 1 | **1000** | 1.33x | not a real regimen — CONTROL |

Two facts from that table are worth stating before any simulation output:

1. **The dose-dense and metronomic regimens are not schedule variants of the standard one — they deliver 2-3x the cumulative dose.** RTOG 0525's dose-dense arm delivered ~2.1x the standard arm's cumulative temozolomide. Any comparison that does not hold cumulative dose fixed is measuring dose, not schedule. Section 3 shows exactly that, and section 4 then holds dose fixed.
2. **The two DIRECTOR arms ARE cumulative-dose-matched by design** (120 mg/m2 x 14 d = 80 mg/m2 x 21 d = 1,680 mg/m2 per cycle). That trial is the closest real-world analogue of the matched-dose comparison in section 4, and it was **negative** (median TTF 1.8 vs 2.0 months, OS 9.8 vs 10.6 months) — though underpowered, having closed early at n=105.

### 2.2 The PK conversion, and what it buries

mg/m2 has to become a plasma curve the transport solver can consume. The conversion is a single arithmetic chain, stated in full so it can be checked or rejected:

| step | value | source |
|---|---|---|
| Oral bioavailability | treated as **1.0** | FDA label: TEMODAR IV over 90 min is **bioequivalent** to 150 mg/m2 oral for both Cmax and AUC (7.5 vs 7.3 ug/mL; 23.4 vs 24.6 ug*h/mL) -> F ~ 0.95-1.0. Baker 1999 (PMID 10037179): *"rapid oral absorption and high systemic availability"* |
| Plasma elimination half-life | **1.8 h** | FDA label (*"mean elimination half-life of 1.8 hours"*); Beale 1999 (PMID 10501912) 1.8 h; Baker 1999 1.9 h; Ostermann 2004 (PMID 15173079) 2.1 h population estimate. Verified range 1.8-2.1 h |
| Molar mass | **194.15 g/mol** | PubChem CID 5394, C6H6N6O2 |
| Anchor: plasma AUC(0-inf) after 150 mg/m2 | **23.4 ug*h/mL** = 120.5 uM*h | FDA label, geometric mean |
| Dose linearity | assumed linear over 75-250 mg/m2/d | FDA label: *"linear kinetics over the therapeutic dosing range of 75 mg/m2/day to 250 mg/m2/day"* |
| Model curve | 1-compartment, instantaneous absorption, mono-exponential elimination -> `AUC = peak / lambda` | `tmesim.treatment.drugs.plasma_profile` |
| **Resulting peak at 150 mg/m2** | **46.4 uM** | `peak = AUC x lambda` |

**The conversion is anchored on AUC, not on Cmax, and the discrepancy is reported rather than hidden.** Tissue exposure in this model is driven by the plasma integral, so AUC is the right invariant to preserve. The price is that the model's implied peak, **46.4 uM**, exceeds the measured geometric-mean Cmax of **38.6 uM** (7.5 ug/mL) by a factor of **1.20**. That is the cost of instantaneous absorption: real oral temozolomide has a median tmax of ~1 h, so the model front-loads a dose that really arrives over about an hour.

**Direction of that error, explicitly.** It makes every modelled dose *more* bolus-like than reality. `results/dose_schedule.md` establishes that in this transport operator the bolus always wins, so the bias favours the more concentrated regimens. It is identical across regimens at matched cumulative dose, so it **cancels from the comparison** and survives only in the absolute levels.

**Other assumptions this conversion buries, named:**

* **No active-species split.** Temozolomide's cytotoxic species is MTIC, formed by pH-dependent hydrolysis; MTIC AUC is 2-4% of parent AUC (Beale 1999). The model carries the parent curve only.
* **No brain-compartment step.** Ostermann 2004 measured `AUC_CSF / AUC_plasma = 20%` in n=35 glioma patients. That partitioning is *not* applied on top of the plasma curve here — the model's own `bbb_penetrance` term (0.125) plays that role, and the two are not the same object. Do not read the plasma peaks below as tumour concentrations.
* **No plasma protein binding correction** (~15% bound per the label).
* **Population-average mg/m2 -> plasma mapping**, not patient-specific.
* **50 mg/m2 (RESCUE) and the rescaled matched-dose arms fall outside the label's stated linear range of 75-250 mg/m2/day.** Rows where that happens are flagged in the tables below.

---

## 3. AS PRESCRIBED — this is a DOSE comparison, not a schedule comparison

Each regimen at its own real dose, over one 28-day window. Cumulative dose differs by 2.8x across the table, so this is a dose axis with a schedule axis confounded into it. It is reported first precisely so that section 4 can be read as the schedule question with the dose confound removed.

| regimen | cumulative mg/m2 / 28 d | peak per dose (uM) | mean tumour AUC (uM*s) | E[kill] | kill (1 draw) | % of ceiling | tumor cells >= IC50 | median E / IC50 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Stupp adjuvant, cycle 1 | 750 | 46.4 | 2.514e+06 | 0.000272 | 0.000168 | 0.05% | 0.00% | 0.00163 |
| Stupp adjuvant, cycles 2-6 | 1000 | 61.9 | 3.352e+06 | 0.000384 | 0.000378 | 0.07% | 0.00% | 0.00217 |
| Stupp concomitant (with RT) | 2100 | 23.2 | 7.04e+06 | 0.000935 | 0.001134 | 0.17% | 0.00% | 0.00456 |
| Dose-dense 21/28, cycle 1 (RTOG 0525) | 1575 | 23.2 | 5.28e+06 | 0.000662 | 0.000756 | 0.12% | 0.00% | 0.00342 |
| Dose-dense 21/28, cycles 2+ (RTOG 0525) | 2100 | 30.9 | 7.04e+06 | 0.000935 | 0.001134 | 0.17% | 0.00% | 0.00456 |
| Dose-intense 7 on / 7 off (Wick) | 2100 | 46.4 | 7.04e+06 | 0.000935 | 0.001134 | 0.17% | 0.00% | 0.00456 |
| DIRECTOR arm A (7/7) | 1680 | 37.1 | 5.632e+06 | 0.000716 | 0.000756 | 0.13% | 0.00% | 0.00365 |
| DIRECTOR arm B (21/28) | 1680 | 24.8 | 5.632e+06 | 0.000716 | 0.000756 | 0.13% | 0.00% | 0.00365 |
| Continuous low-dose (RESCUE) * | 1400 | 15.5 | 4.693e+06 | 0.000575 | 0.000630 | 0.10% | 0.00% | 0.00304 |
| Single large dose (CONTROL) * | 1000 | 309.4 | 3.352e+06 | 0.000384 | 0.000378 | 0.07% | 0.00% | 0.00217 |

`*` = per-dose mg/m2 outside the label's stated linear-PK range (75-250 mg/m2/day); the linear conversion is an extrapolation for that row.

**Ranking is monotone in cumulative dose, and in nothing else** (checked, not asserted: sorting the table by cumulative mg/m2 gives a non-decreasing E[kill] at every step). The top arm (Dose-dense 21/28, cycles 2+ (RTOG 0525), 2100 mg/m2) reaches E[kill] 0.000935; the bottom (Stupp adjuvant, cycle 1, 750 mg/m2) reaches 0.000272, a 3.44x range across a 2.8x range in cumulative dose. That is the entire content of this table: in a transport model with a static Hill kill map, cumulative dose is the only lever visible at this scale.

> **This is where the model and the clinic diverge, and it is worth being blunt about it.** RTOG 0525's dose-dense arm delivered ~2.1x the cumulative dose of the standard arm and produced **no overall-survival benefit** (median OS 14.9 vs 16.6 months, HR 1.03, P = .63) while doubling grade >= 3 toxicity (53% vs 34%). A model in which outcome is a monotone function of delivered exposure — which is exactly what this model is — **cannot reproduce that trial**. Section 7 lists what would have to be added before it could.

Note also that IC50 coverage never exceeds 0.00% of tumor cells on any arm. Temozolomide's cached GDSC glioma IC50 is 606 uM against a median window-averaged tumour exposure of 2.764 uM even on the highest-dose arm. This is consistent with `results/dose_schedule.md`, which finds temozolomide dose-limited across its whole swept range and never reaching IC50 at clinical plasma levels.

---

## 4. AT MATCHED CUMULATIVE DOSE — schedule as the ONLY variable

Every regimen is now rescaled to deliver exactly **1000 mg/m2 over the 28-day window** (the standard Stupp adjuvant cycle at 200 mg/m2 x 5 days). Only the *timing* of that fixed total differs. This is the comparison the question actually asks for.

| regimen | doses | mg/m2 per dose | peak per dose (uM) | plasma AUC COUNTED in window (uM*s) | vs best | **mean tumour AUC (uM*s)** | vs best | tissue AUC per counted plasma AUC | vs best | E[kill] | tumor >= IC50 | terminal residual C(T) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Stupp adjuvant, cycle 1 | 5 | 200.0 | 61.88 | 2894156.154 | +0.000000% | **3352293.869** | -0.000000% | 1.158297511 | +0.000000% | 0.00038424 | 0.00% | 3.450e-08% |
| Stupp adjuvant, cycles 2-6 | 5 | 200.0 | 61.88 | 2894156.154 | +0.000000% | **3352293.869** | -0.000000% | 1.158297511 | +0.000000% | 0.00038424 | 0.00% | 3.450e-08% |
| Stupp concomitant (with RT) * | 28 | 35.7 | 11.05 | 2894146.139 | -0.000346% | **3352281.839** | -0.000359% | 1.158297363 | -0.000013% | 0.00038424 | 0.00% | 1.275e-05% |
| Dose-dense 21/28, cycle 1 (RTOG 0525) * | 21 | 47.6 | 14.73 | 2894156.154 | -0.000000% | **3352293.869** | -0.000000% | 1.158297511 | +0.000000% | 0.00038424 | 0.00% | 3.450e-08% |
| Dose-dense 21/28, cycles 2+ (RTOG 0525) * | 21 | 47.6 | 14.73 | 2894156.154 | -0.000000% | **3352293.869** | -0.000000% | 1.158297511 | +0.000000% | 0.00038424 | 0.00% | 3.450e-08% |
| Dose-intense 7 on / 7 off (Wick) * | 14 | 71.4 | 22.10 | 2894156.154 | +0.000000% | **3352293.869** | +0.000000% | 1.158297511 | +0.000000% | 0.00038424 | 0.00% | 3.450e-08% |
| DIRECTOR arm A (7/7) * | 14 | 71.4 | 22.10 | 2894156.154 | +0.000000% | **3352293.869** | +0.000000% | 1.158297511 | +0.000000% | 0.00038424 | 0.00% | 3.450e-08% |
| DIRECTOR arm B (21/28) * | 21 | 47.6 | 14.73 | 2894156.154 | -0.000000% | **3352293.869** | -0.000000% | 1.158297511 | +0.000000% | 0.00038424 | 0.00% | 3.450e-08% |
| Continuous low-dose (RESCUE) * | 28 | 35.7 | 11.05 | 2894146.139 | -0.000346% | **3352281.839** | -0.000359% | 1.158297363 | -0.000013% | 0.00038424 | 0.00% | 1.275e-05% |
| Single large dose (CONTROL) * | 1 | 1000.0 | 309.41 | 2894156.154 | -0.000000% | **3352293.869** | -0.000000% | 1.158297511 | +0.000000% | 0.00038424 | 0.00% | 3.450e-08% |

**Total spread in mean tumour exposure: 0.000359%.** The analytic bound is `tau_tissue / T_window` = 331 s / 28 d = 0.013692%, and the measured spread sits under it, as it must.

**Reading the three 'vs best' columns — the difference decomposes cleanly.** Matching *administered* cumulative mg/m2 is the clinically meaningful match, but it is not quite the same thing as matching the drug the solver actually integrates. Two distinct terminal effects separate them, and both columns are here so neither is hidden:

1. **Drug not yet delivered when the window closes** — a dose given on day 28 has not finished being absorbed and eliminated by day 28, so part of its plasma curve falls outside the integration. That is the *plasma AUC counted in window* column, and it is a **measurement-window artefact**, not a schedule effect.
2. **Drug delivered but still resident in tissue at T** — the terminal residual `C(T)`, the last column. This one IS the schedule effect: it is the `C(T)` in `A = M^-1 (S_total - C(T))`.

The **tissue AUC per counted plasma AUC** column divides the first effect out and leaves only the second, so it is the cleanest possible statement of the pure schedule effect at matched dose. Its spread across all arms is **0.000013%**.

**The Bernoulli-drawn kill fraction is IDENTICAL (0.000378) for every arm**, because the exposure differences above are far below the 1/23,812 = 4.20e-05 quantum of a kill fraction obtained by counting cells. Expected kill (the mean per-cell probability) is reported instead throughout, because it is continuous and carries no counting noise. This is worth stating plainly: **at matched cumulative dose this model does not distinguish the regimens even in principle at the resolution of a cell count.**

### 4.1 Why the spread is where it is, exactly

`results/dose_schedule.md` section 4b derives it. Integrating the transport PDE over `[0, T]` with `C(0) = 0`, for `A = AUC` and `M = k - div(D grad .) + div(v .)`:

```
    A = M^-1 ( S_total - C(T) )
```

`S_total` — the drug the solver actually injects over `[0, T]` — is **identical** for any two schedules at matched *counted* exposure, and when it is, the entire schedule effect is the terminal residual `C(T)`. That is the pure-transport term, and here it is **0.000013%**.

Matching *administered* mg/m2 rather than counted exposure adds the second term, and on a 28-day window that second term **dominates**: 0.000346% of the total 0.000359% is drug whose plasma curve had not finished running when the window closed, versus 0.000013% of genuine terminal residual — a ratio of about 27 to 1. **The larger of the two effects is therefore a bookkeeping artefact of where the window edge falls, not physics.** Both are reported so that neither can be quoted as the other, and section 4.3 re-measures the pure term alone at three window lengths.

Concretely, the mechanism is visible in the table: the 5/28, 21/28, 7/7 and single-bolus arms all finish dosing at least seven days before the window closes, so their drug has washed out completely, they realise the **full** steady-state exposure field, and they are identical to one another **to every stored digit** — four schedule families ranging from 1 to 21 fractions, indistinguishable. Only the two continuous-daily arms, which dose on day 28 itself, differ at all, and only in the last dose's tail.

### 4.2 The same comparison up the dose axis

The schedule effect is a fixed **relative** change in exposure, so the **absolute** kill spread grows with where you sit on the Hill curve. At clinical amplitude temozolomide's predicted kill on this section is tiny and a sub-0.02% relative change lands far below the 1/23,812 = 4.20e-05 quantum of a cell-counted kill fraction. Rescaling up the dose axis is exactly free (exposure is linear in amplitude) and lifts it clear. **The multipliers below are a diagnostic amplitude, NOT a dose — 1000x clinical is not a regimen anyone could give.** Expected kill (mean per-cell probability) is reported because it carries no cell-counting quantisation.

| amplitude (x clinical) | Stupp adjuvant, cycles 2-6 | Dose-dense 21/28, cycles 2+ (RTOG 0525) | Dose-intense 7 on / 7 off (Wick) | Stupp concomitant (with RT) | Single large dose (CONTROL) | spread |
|---:|---:|---:|---:|---:|---:|---:|
| 1x | 0.00038424 | 0.00038424 | 0.00038424 | 0.00038424 | 0.00038424 | **1.65e-09** |
| 10x | 0.00601987 | 0.00601987 | 0.00601987 | 0.00601984 | 0.00601987 | **2.56e-08** |
| 100x | 0.08095723 | 0.08095723 | 0.08095723 | 0.08095694 | 0.08095723 | **2.93e-07** |
| 1000x | 0.39577448 | 0.39577448 | 0.39577448 | 0.39577401 | 0.39577448 | **4.67e-07** |

Every one of those spreads would have been **exactly `0.00e+00`** under the static source, for every amplitude, because the inputs were bit-identical. They are now small, non-zero and bounded — which is a different statement from the one the old engine forced.

### 4.3 The bound is a ratio, so it moves with the observation window

`tau_tissue / T_window` is the whole story, which means the size of the measured schedule effect is a **reporting choice** as much as a physical one. The same burst-vs-continuous contrast — now matched on **counted plasma exposure**, so this is the PURE transport term of section 4.1 with the window-truncation artefact removed — re-measured over three windows:

| window | burst 5/28 mean tumour AUC | continuous daily mean tumour AUC | single bolus | measured spread | bound tau/T |
|---|---:|---:|---:|---:|---:|
| 7 d | 3352293.9 | 3352292.2 | 3352293.9 | **0.000051%** | 0.054766% |
| 14 d | 3352293.9 | 3352293 | 3352293.9 | **0.000026%** | 0.027383% |
| 28 d | 3352293.9 | 3352293.4 | 3352293.9 | **0.000013%** | 0.013692% |

The measured spread halves when the window doubles, exactly as `1/T` predicts. **Quote any schedule percentage from this model with its window attached, or it means nothing.**

### 4.4 Where the effect IS resolvable: intra-day fractionation

To show the machinery is measuring physics rather than round-off, the same matched-cumulative-exposure comparison is run at a 24 h window, where `tau/T` = 0.3834% and the effect is comfortably above any numerical floor. One day's dose given as 1, 2 or 4 fractions:

| schedule | peak per dose (uM) | plasma AUC in window (uM*s) | mean tumour AUC (uM*s) | E[kill] | tumor >= IC50 | terminal residual |
|---|---:|---:|---:|---:|---:|---:|
| 1 dose (bolus) | 61.88 | 578775 | 670391.41 | 0.003020 | 0.00% | 0.000% |
| 2 doses q12h | 31.09 | 578775 | 670269.29 | 0.003020 | 0.00% | 0.018% |
| 4 doses q6h | 15.91 | 578775 | 669687.97 | 0.003016 | 0.00% | 0.104% |

**Exposure spread 0.1049%** against the bound 0.3834%, and it is **monotone: the bolus wins.** That is the same sign `results/dose_schedule.md` found in every regime it tested (its section 5 records the bolus beating the sustained infusion across a 192-fold plasma-half-life range, which is the **opposite** of the textbook transport-limited intuition). The reason is structural — the earlier a schedule finishes, the more of the steady-state exposure field it realises inside the window, and there is no saturable perivascular sink in this transport operator to reverse it. The terminal-residual column shows the mechanism directly: it rises with the number of fractions.

| amplitude | 1 dose (bolus) | 2 doses q12h | 4 doses q6h | spread (expected kill) |
|---:|---:|---:|---:|---:|
| 1x | 0.003020 | 0.003020 | 0.003016 | **3.76e-06** |
| 10x | 0.043898 | 0.043889 | 0.043848 | **5.03e-05** |
| 100x | 0.311115 | 0.311086 | 0.310947 | **1.69e-04** |
| 1000x | 0.523118 | 0.523112 | 0.523085 | **3.26e-05** |

Note that the absolute kill spread **peaks and then falls** — it is largest at 100x (1.69e-04) and smaller again at 1000x (3.26e-05). A fixed *relative* change in exposure produces the largest *absolute* change in kill where the Hill curve is steepest, and the top amplitude here is already pressed against the kill ceiling (0.5502), where every schedule converges on the same saturated answer. **Escalating dose does not amplify a schedule difference indefinitely; it eventually erases it.**

---

## 5. COMBINATIONS, re-run with time resolution

> ### THE COMBINATION ASSUMPTION: BLISS INDEPENDENCE — read this before any number
>
> Each agent is simulated **separately** on the same real tissue, giving every tumor cell its own exposure AUC and its own kill probability from the engine's own Hill pharmacodynamics. The combination is formed at the exposure level by
>
> ```
> p_combo(cell) = 1 - (1 - p_A(cell)) * (1 - p_B(cell))
> ```
>
> **NO pharmacological synergy or antagonism is modelled.** There is no shared target, no DNA-repair crosstalk, no efflux competition, no drug-drug PK interaction, and drug B does not alter drug A's transport. What is modelled is **SPATIAL complementarity only** — whether the two delivery fields reach *different cells*. Every gain below is geographic or arithmetic, never biochemical. Real combinations may synergise or antagonise and this model is blind to both.
>
> Because Bliss independence *is* the statement that the per-cell draws are independent, all set statistics are exact expectations under that coupling (`E|A n B| = sum p_A p_B`, `E|A u B| = sum [1-(1-p_A)(1-p_B)]`) and carry **no Monte-Carlo noise**.

### 5.1 Single agents

Two PK conventions are run, because **plasma half-life is the parameter a stagger acts on** and answering the staggering question only on a shared placeholder would be evasive:

| convention | temozolomide | cisplatin | doxorubicin | rationale |
|---|---|---|---|---|
| `shared` | 50.0 uM, 4.0 h | 50.0 uM, 4.0 h | 50.0 uM, 4.0 h | identical to `results/combination.md`, so the ONLY differences from that report are time resolution and staggering |
| `clinical` | 46.4 uM, 1.8 h | 8.0 uM, 0.5 h | 2.8 uM, 30.0 h | each agent on its own measured plasma half-life — a **60x span**, which is exactly the axis staggering acts on |

**Provenance of the clinical PK, with its weaknesses stated:**

* **temozolomide** — 1.8 h half-life (FDA label; Beale 1999, PMID 10501912), peak from the section 2.2 conversion at 150 mg/m2. Solid.
* **cisplatin** — 0.5 h for **free / ultrafilterable** platinum (Himmelstein 1981, PMID 7194166: *"t1/2 of 0.3 to 0.5 hr"*; total platinum is triphasic with a terminal t1/2 >= 24 h and is the wrong species here). The peak is **DERIVED, not measured**: unchanged cisplatin Cmax 1.5-2 ug/mL after 80 mg/m2 (Nagai 1996, PMID 8995510) divided by Himmelstein's cisplatin/filterable-Pt ratio of 0.6-0.8 gives ~6-11 uM; 8 uM is the midpoint. Flag this. Also note that free-Pt elimination is dominated by irreversible covalent protein binding, which a single first-order clearance term mis-specifies.
* **doxorubicin** — 30.0 h terminal half-life (Benjamin 1977, PMID 856462, *"triphasic ... mean half-life of 30 hr"*; Camaggi 1988 gives t1/2-gamma 48.4 h; interpatient CV ~65%). Peak 2.8 uM is the measured 1509 ng/mL at 35 mg/m2 IV bolus (Muller 1993, PMID 8339389). **The peak is dominated by sampling time** — doxorubicin's distribution half-life is ~4 min (Camaggi 1988 t1/2-alpha 4.8 min), so any reported Cmax is as much a statement about when blood was drawn as about the drug. This is the least trustworthy number in the report. Speth 1988 (PMID 3042244) adds that intracellular doxorubicin concentrations reached in vivo sit **below** the in-vitro clonogenic-kill threshold and that the optimal anthracycline schedule *"has never been investigated in a prospective, randomised manner"* — i.e. for this agent, plasma exposure is known NOT to determine effect, which is precisely the assumption this model rests on.

Single agents over the 24 h combination window, dosed at t = 0:

**`shared` convention**

| drug | peak (uM) | plasma t1/2 (h) | plasma AUC in window (uM*s) | mean tumour AUC (uM*s) | E[kill] | tumor >= IC50 | median E / IC50 | terminal residual |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| temozolomide | 50.00 | 4.0 | 1.023e+06 | 1.184e+06 | 0.005943 | 0.00% | 0.0215 | 0.026% |
| cisplatin | 50.00 | 4.0 | 1.023e+06 | 2.789e+06 | 0.310409 | 71.23% | 1.2164 | 0.009% |
| doxorubicin | 50.00 | 4.0 | 1.023e+06 | 1.755e+04 | 0.133878 | 2.73% | 0.3511 | 0.002% |

**`clinical` convention**

| drug | peak (uM) | plasma t1/2 (h) | plasma AUC in window (uM*s) | mean tumour AUC (uM*s) | E[kill] | tumor >= IC50 | median E / IC50 | terminal residual |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| temozolomide | 46.41 | 1.8 | 4.341e+05 | 5.028e+05 | 0.002142 | 0.00% | 0.0091 | 0.000% |
| cisplatin | 8.00 | 0.5 | 2.081e+04 | 5.677e+04 | 0.007319 | 0.00% | 0.0248 | 0.000% |
| doxorubicin | 2.80 | 30.0 | 1.857e+05 | 3185 | 0.024017 | 0.00% | 0.0637 | 0.026% |

### 5.2 Pairs, with the dose-matched controls that decide the question

Two controls, because "the combination kills more" is not a result on its own:

* **2x dose of the better single agent** — the same total drug, in one agent. If the pair cannot beat this, it is not worth two drugs.
* **Bliss(better + ITSELF)**, i.e. `1 - (1 - p_best)^2` — a hypothetical arm whose two components have *perfectly redundant* delivery. This is the control that separates spatial complementarity from arithmetic: under a saturating Hill curve, **any** second independent hit beats **any** second dose, including the same drug twice. If `combo(A,B) <= selfcombo(B,B)`, the apparent benefit is the independence assumption meeting a saturating curve, not geography.

**`shared` convention, simultaneous dosing, 24 h window**

| pair | A alone | B alone | A+B (Bliss) | better | 2x better | Bliss(better+ITSELF) | combo / better | combo / 2x-dose | combo / self-combo |
|---|---:|---:|---:|---|---:|---:|---:|---:|---:|
| temozolomide + cisplatin | 0.0059 | 0.3104 | 0.3142 | cisplatin | 0.4081 | 0.5103 | 1.01x | 0.77x | **0.62x** |
| temozolomide + doxorubicin | 0.0059 | 0.1339 | 0.1388 | doxorubicin | 0.2246 | 0.2428 | 1.04x | 0.62x | **0.57x** |
| cisplatin + doxorubicin | 0.3104 | 0.1339 | 0.3939 | cisplatin | 0.4081 | 0.5103 | 1.27x | 0.97x | **0.77x** |

**`clinical` convention, simultaneous dosing, 24 h window**

| pair | A alone | B alone | A+B (Bliss) | better | 2x better | Bliss(better+ITSELF) | combo / better | combo / 2x-dose | combo / self-combo |
|---|---:|---:|---:|---|---:|---:|---:|---:|---:|
| temozolomide + cisplatin | 0.0021 | 0.0073 | 0.0094 | cisplatin | 0.0165 | 0.0146 | 1.29x | 0.57x | **0.65x** |
| temozolomide + doxorubicin | 0.0021 | 0.0240 | 0.0261 | doxorubicin | 0.0510 | 0.0470 | 1.09x | 0.51x | **0.55x** |
| cisplatin + doxorubicin | 0.0073 | 0.0240 | 0.0311 | doxorubicin | 0.0510 | 0.0470 | 1.29x | 0.61x | **0.66x** |

The last column is the one that matters. A ratio at or below 1.00 means combining two *different* drugs did no better than the model's own independence arithmetic applied to a single drug against itself.

### 5.3 STAGGERED vs SIMULTANEOUS — the genuinely new question

The old engine could not pose this at all: with the plasma curve averaged away, "drug B at t = 0" and "drug B at t = 12 h" were the same input array. Four arms over a 24 h window:

* **simultaneous** — both agents at `t = 0`.
* **B delayed, matched administered dose** — same peak, B given at `t = 12 h`. B then has less window left in which to be delivered *and counted*, so it loses exposure. That is a **measurement-window artefact**, the same one `results/dose_schedule.md` names as its 'leg B'. It is reported because it is what a naive "same milligrams" comparison would produce, and it should not be mistaken for a schedule effect.
* **B delayed, matched COUNTED exposure** — B's peak is rescaled so its plasma AUC *inside the window* matches the simultaneous arm. `S_total` is then identical and the only remaining difference is the terminal residual. **This is the leg that isolates schedule.**
* **A delayed, matched counted exposure** — the reverse order, because which agent goes first is a real clinical choice.

Jaccard is the overlap of the killed populations, `E|A n B| / E|A u B|`. Because independent-Bernoulli overlap shrinks automatically as kill probabilities shrink, the raw Jaccard proves nothing on its own; the **delivery-shuffle null** (200 permutations, 20 strata) permutes only drug B's delivery term within strata of the shared cell-intrinsic killability factor, so each cell keeps its own intrinsic killability and only B's delivery *geography* is scrambled. `obs/null` < 1 = complementary territory, ~1 = no spatial structure, > 1 = redundant.

**`shared` convention**

| pair | arm | A kill | B kill | combo kill | **Jaccard** | delivery-shuffle null | obs/null | both (% of killed) | A only | B only |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| temozolomide + cisplatin | simultaneous | 0.0059 | 0.3104 | 0.3142 | **0.006979** | 0.006664 | 1.047x | 0.70% | 1.19% | 98.11% |
| temozolomide + cisplatin | B delayed, matched administered | 0.0059 | 0.2919 | 0.2958 | **0.007009** | 0.006667 | 1.051x | 0.70% | 1.31% | 97.99% |
| temozolomide + cisplatin | B delayed, matched counted exposure | 0.0059 | 0.3103 | 0.3140 | **0.006979** | 0.006664 | 1.047x | 0.70% | 1.19% | 98.11% |
| temozolomide + cisplatin | A delayed, matched counted exposure | 0.0059 | 0.3104 | 0.3142 | **0.006962** | 0.006648 | 1.047x | 0.70% | 1.19% | 98.11% |
| temozolomide + doxorubicin | simultaneous | 0.0059 | 0.1339 | 0.1388 | **0.007364** | 0.006645 | 1.108x | 0.74% | 3.55% | 95.72% |
| temozolomide + doxorubicin | B delayed, matched administered | 0.0059 | 0.1210 | 0.1260 | **0.007372** | 0.006624 | 1.113x | 0.74% | 3.98% | 95.28% |
| temozolomide + doxorubicin | B delayed, matched counted exposure | 0.0059 | 0.1339 | 0.1388 | **0.007365** | 0.006645 | 1.108x | 0.74% | 3.55% | 95.72% |
| temozolomide + doxorubicin | A delayed, matched counted exposure | 0.0059 | 0.1339 | 0.1388 | **0.007348** | 0.006630 | 1.108x | 0.73% | 3.54% | 95.73% |
| cisplatin + doxorubicin | simultaneous | 0.3104 | 0.1339 | 0.3939 | **0.127867** | 0.118985 | 1.075x | 12.79% | 66.01% | 21.20% |
| cisplatin + doxorubicin | B delayed, matched administered | 0.3104 | 0.1210 | 0.3857 | **0.118428** | 0.109992 | 1.077x | 11.84% | 68.63% | 19.53% |
| cisplatin + doxorubicin | B delayed, matched counted exposure | 0.3104 | 0.1339 | 0.3939 | **0.127849** | 0.118968 | 1.075x | 12.78% | 66.02% | 21.20% |
| cisplatin + doxorubicin | A delayed, matched counted exposure | 0.3103 | 0.1339 | 0.3938 | **0.127859** | 0.118971 | 1.075x | 12.79% | 66.00% | 21.21% |

**`clinical` convention**

| pair | arm | A kill | B kill | combo kill | **Jaccard** | delivery-shuffle null | obs/null | both (% of killed) | A only | B only |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| temozolomide + cisplatin | simultaneous | 0.0021 | 0.0073 | 0.0094 | **0.002186** | 0.001936 | 1.129x | 0.22% | 22.47% | 77.31% |
| temozolomide + cisplatin | B delayed, matched administered | 0.0021 | 0.0073 | 0.0094 | **0.002186** | 0.001936 | 1.129x | 0.22% | 22.47% | 77.31% |
| temozolomide + cisplatin | B delayed, matched counted exposure | 0.0021 | 0.0073 | 0.0094 | **0.002186** | 0.001936 | 1.129x | 0.22% | 22.47% | 77.31% |
| temozolomide + cisplatin | A delayed, matched counted exposure | 0.0021 | 0.0073 | 0.0094 | **0.002186** | 0.001935 | 1.129x | 0.22% | 22.47% | 77.32% |
| temozolomide + doxorubicin | simultaneous | 0.0021 | 0.0240 | 0.0261 | **0.002686** | 0.002315 | 1.160x | 0.27% | 7.94% | 91.79% |
| temozolomide + doxorubicin | B delayed, matched administered | 0.0021 | 0.0126 | 0.0147 | **0.002521** | 0.002158 | 1.168x | 0.25% | 14.30% | 85.45% |
| temozolomide + doxorubicin | B delayed, matched counted exposure | 0.0021 | 0.0240 | 0.0261 | **0.002686** | 0.002315 | 1.160x | 0.27% | 7.95% | 91.79% |
| temozolomide + doxorubicin | A delayed, matched counted exposure | 0.0021 | 0.0240 | 0.0261 | **0.002685** | 0.002314 | 1.160x | 0.27% | 7.94% | 91.79% |
| cisplatin + doxorubicin | simultaneous | 0.0073 | 0.0240 | 0.0311 | **0.008499** | 0.006675 | 1.273x | 0.85% | 22.71% | 76.44% |
| cisplatin + doxorubicin | B delayed, matched administered | 0.0073 | 0.0126 | 0.0198 | **0.007101** | 0.005513 | 1.288x | 0.71% | 36.27% | 63.02% |
| cisplatin + doxorubicin | B delayed, matched counted exposure | 0.0073 | 0.0240 | 0.0311 | **0.008498** | 0.006674 | 1.273x | 0.85% | 22.71% | 76.44% |
| cisplatin + doxorubicin | A delayed, matched counted exposure | 0.0073 | 0.0240 | 0.0311 | **0.008499** | 0.006675 | 1.273x | 0.85% | 22.71% | 76.44% |

**The answer, stated plainly.** Staggering changes the overlap of the killed populations by:

| pair | convention | Jaccard, simultaneous | Jaccard, B delayed (matched counted) | absolute change | relative change |
|---|---|---:|---:|---:|---:|
| temozolomide + cisplatin | `shared` | 0.006979 | 0.006979 | +2.35e-07 | **+0.003%** |
| temozolomide + doxorubicin | `shared` | 0.007364 | 0.007365 | +3.53e-08 | **+0.000%** |
| cisplatin + doxorubicin | `shared` | 0.127867 | 0.127849 | -1.77e-05 | **-0.014%** |
| temozolomide + cisplatin | `clinical` | 0.002186 | 0.002186 | -1.45e-12 | **-0.000%** |
| temozolomide + doxorubicin | `clinical` | 0.002686 | 0.002686 | -5.81e-08 | **-0.002%** |
| cisplatin + doxorubicin | `clinical` | 0.008499 | 0.008498 | -6.54e-07 | **-0.008%** |

Largest relative change in Jaccard from staggering, across all six pair-by-convention cells: **0.014%**.

> **And this is a STRUCTURAL result, not an empirical one — which is why it should be believed and why it is also nearly useless clinically.** In this model the transport operator contains **no drug-drug coupling of any kind**: each agent's exposure field is solved independently, and Bliss independence then combines the two per-cell kill probabilities without reference to timing. Delaying agent B therefore cannot change *where* B goes — only, very slightly, *how much* of B's steady-state field is realised inside the observation window, and that is the terminal residual again, bounded by `tau/T`. **Staggering can only move combination overlap in this model through an amplitude change, never through a geometry change.** The clinical reasons to stagger — resistance-reversal windows, cell-cycle synchronisation, sequential DNA-damage-and-repair-block, non-overlapping toxicity — are all absent. So take the ~0% as a statement about the model's structure and NOT as evidence that sequencing does not matter.

The one place staggering does move something materially is the **matched-administered-dose** row, and that movement is a window artefact in exactly the sense named above — the delayed agent simply has less of the window left in which to be counted. The largest such shift is **cisplatin + doxorubicin (`clinical`)**, Jaccard 0.008499 -> 0.007101 (16.44%), against a largest matched-COUNTED shift of 0.014%. The gap between those two numbers is the whole reason both legs are reported.

---

## 6. Verdict

> **1. At matched cumulative dose over a real 28-day cycle, this model cannot tell the real temozolomide regimens apart.** The spread in mean tumour exposure is 0.000359% across 10 arms spanning 1 to 28 fractions, against the derived bound 0.013692% — and of that spread, only 0.000013% is transport at all; the rest is the window edge cutting off the last dose's plasma tail. Every regimen that finishes dosing more than a day before the window closes — 1 fraction, 5, 14, 21 — is identical to every other one **to the last stored digit**.
>
> **2. That indistinguishability is a structural fact about transport, not a finding about temozolomide.** The entire schedule channel in a linear transport model is the terminal residual, bounded by `tau_tissue / T_window`. Temozolomide equilibrates in this tissue in 331 s against a 28-day window, so the ratio is ~1e-4. The same machinery would show a large effect for a slowly-equilibrating molecule (an antibody, `tau` of hours to days) or over a short window — section 4.4 demonstrates that by shrinking the window to 24 h and recovering a 0.105% effect from the identical code.
>
> **3. Where the model DOES separate regimens, it is separating them by dose.** As prescribed, the dose-dense and metronomic schedules deliver 2-3x the cumulative temozolomide of the standard adjuvant cycle, and the predicted ranking is monotone in that. **RTOG 0525 showed that this is the wrong prediction**: 2.1x the dose, no survival benefit, double the grade >=3 toxicity. Do not read section 3's ranking as a recommendation; read it as a demonstration of what an exposure-driven model does when the biology is not exposure-driven.
>
> **4. Staggering does not change which cells a combination kills.** Largest relative Jaccard change across six pair-by-convention cells: 0.014%. This is structural — there is no drug-drug coupling in the transport operator and Bliss independence is timing-blind — so it is a strong statement about the model and a weak one about sequencing therapy.
>
> **5. Against the dose-matched controls, the combinations do not survive.** Of the 6 pair-by-convention cells, **0 beat 2x the better single agent** and **0 beat Bliss(better + ITSELF)** (combo/self-combo range 0.55x-0.77x). This reproduces `results/combination.md`'s finding under a time-resolved solver: the apparent combination benefit over 2x dosing is the independence arithmetic meeting a saturating Hill curve, not spatial complementarity. The overlap ratios say the same thing from the other side: `obs/null` ranges 1.05x-1.29x across every pair, convention and arm — i.e. the two agents' killed populations co-localise **more** than a delivery-scrambled null in every single arm, not less. There is no spatial complementarity to exploit on this section.

---

## 7. THE HONEST CEILING: what makes fractionation work is not in this model

This section is not a caveat appended to a result. It is the result's load-bearing context, and the numbers above should not be quoted without it.

**The schedule effect measured here is a TRANSPORT-ONLY FLOOR.** Time resolution was added to the transport solve. The pharmacodynamic map is unchanged: a static Hill function of window-averaged cumulative exposure. Everything that gives fractionation and dose-density their clinical rationale is therefore missing:

| missing mechanism | why it is a schedule mechanism | evidence it is the real story |
|---|---|---|
| **MGMT / AGAT depletion and regeneration** | MGMT is a *consumable, stoichiometric, suicide* repair protein. Protracted dosing depletes it and it regenerates between doses — so the effective sensitivity is a state variable with its own time constant, which is precisely what a schedule acts on. **This model has no MGMT term at all** (MGMT is not on the Xenium panel; sensitivity is uniform 1.0). | Tolcher 2003 (PMID 12671695): mean AGAT activity fell 72% after 7 days and 73% after 21 days of TMZ in PBMCs, and stayed depressed >= 7 days after dosing. Note the plateau: 7 -> 21 days bought essentially nothing, which itself undercuts the dose-dense premise. |
| **Inter-fraction DNA repair** | The entire radiobiological basis of fractionation. Absent. | — |
| **Cell-cycle redistribution and S-phase gating** | TMZ cytotoxicity is MMR-dependent and requires replication cycles, so *when* a cell is hit relative to its cycle matters. This model's cell states are oxygen-banded and static through the whole window. | mechanism widely described; not independently sourced in this report |
| **Regrowth / repopulation between fractions** | No proliferation is simulated during the window at all. | — |
| **Immune contribution, and TMZ-induced lymphopenia** | Both absent. The second is a *schedule-dependent harm* that dose-dense regimens are known to cause. | Gilbert 2013: grade >= 3 toxicity 34% -> 53% on dose-dense, largely lymphopenia and fatigue |
| **Normal-tissue recovery — the entire denominator of the therapeutic index** | Every number in this report is tumour kill alone. A schedule that is near-identical for the tumor is emphatically not identical for the patient. | Gilbert 2013, as above |

**And the clinical record says exposure is not the variable.** Ostermann 2004 (PMID 15173079), which measured both plasma and CSF temozolomide in n=35 glioma patients, concluded that *"systemic or cerebral exposures are not better predictors than the cumulative dose alone for both efficacy and safety"* — and RTOG 0525 then showed cumulative dose itself does not predict survival. A model whose output is a monotone function of tumour exposure is, by construction, modelling the wrong variable for this question. **What it can legitimately do is bound the transport contribution, which is what it does above.**

> **Do not present any number in this file as a clinical recommendation.** A 0.0001% exposure difference between two regimens is not a reason to choose one. It is a demonstration that transport is not where the difference lives.

---

## 8. Limitations (the full list)

1. **No outcome validation, at all.** Untreated resection; no measured response for any agent on this tissue. Nothing here has been scored against reality.
2. **Not clinical dosing guidance.** Research prototype. No regulatory, diagnostic or treatment-selection use.
3. **Absolute kill on this section is not defensible, per this repo's own `results/cell_state_impact.md`.** That report finds that on the 319-gene flagship section — this one — no cell-state route clears the bar, absolute kill moves by up to 0.3252 absolute (1.77x) purely with the state route, and its explicit recommendation is to **stop quoting an absolute kill fraction or a kill ceiling as a result** here. Rankings and relative comparisons are what this file is for; every absolute kill and the ceiling are unit-carrying intermediates printed so the comparisons can be checked.
4. **The schedule effect is bounded by `tau_tissue / T_window` and is therefore window-conditional.** 331 s / 28 d = 0.013692%. Section 4.3 shows the measured value tracking `1/T` across three windows. Any schedule percentage from this model is meaningless without its window.
5. **Pharmacodynamics is a single static Hill map on window-averaged exposure — and this is the binding limitation on the whole question.** No MGMT depletion, no inter-fraction repair, no cell-cycle redistribution, no regrowth, no resistance evolution, no immune contribution, no normal-tissue toxicity. Section 7.
6. **Kill depends on `AUC / T_window`, so lengthening the observation window mechanically lowers predicted kill at fixed total dose.** This is why the 28-day-window kills here are far below the 6 h-window kills in `results/dose_schedule.md` and `results/combination.md`. The window is held EXACTLY constant across everything compared, so the comparisons are sound; cross-file comparisons of absolute kill are not.
7. **Bliss independence is an assumption.** No synergy, no antagonism, no drug-drug PK interaction, no shared target, no efflux competition. The combination results model SPATIAL complementarity only and are a delivery-only floor, not a forecast of combination efficacy.
8. **Staggering can only act through amplitude in this model.** There is no drug-drug coupling in the transport operator and Bliss independence is timing-blind, so the ~0% overlap change is a structural consequence and not an empirical discovery about sequencing.
9. **The PK conversion is anchored on AUC and overstates the peak by 1.20x** (instantaneous absorption vs a real tmax of ~1 h). It also carries no MTIC species, no CSF partitioning, no protein binding, and no patient-specific mg/m2 mapping. Section 2.2.
10. **Two of the three clinical peaks in the combination leg are weak.** Cisplatin's free-platinum peak is DERIVED (from unchanged-cisplatin Cmax divided by a species ratio), not measured. Doxorubicin's peak is sampling-time-dominated (distribution t1/2 ~4 min) with interpatient CV ~65%. Treat both as order-of-magnitude.
11. **Per-cell sensitivity is uniform 1.0.** Temozolomide's determinant (MGMT) is not on this panel, so TMZ resistance is unmodelled and its absolute kill is optimistic. The `uniform` mode makes the cross-drug combination comparison fair by making every agent equally an upper bound; it does not make any of them right.
12. **Cell-state banding dominates the ceiling.** 90.0% of tumor cells are labelled QUIESCENT (x0.5 kill) and 0.0% HYPOXIC (x0.35) by the oxygen model. That labelling is a *model* assignment from a simulated oxygen field, not a measurement.
13. **The perivascular gradient partly inherits the vascular source-field shape** (`tmesim/validation/vessel_profile.py`), and the exposure geometry every spatial statistic rests on inherits that confound.
14. **One calibration constant** (alpha = 148.6, fit once on doxorubicin) sets the perivascular sink for every drug.
15. **Single section, single patient, 2D.** No inter-patient variability, no out-of-plane transport, no cohort, no replicate sections.
16. **Regimens are modelled over ONE 28-day cycle, not a full course.** Six adjuvant cycles, the 42-day concomitant phase, and any cycle-to-cycle carry-over are not represented. With no regrowth and no resistance evolution in the model, repeating the cycle would add nothing but exposure.
17. **The single-dose superposition is exact for this solver, but it is exact BECAUSE the solver is linear and time-invariant.** Any future change that makes transport nonlinear (saturable binding, concentration-dependent efflux, drug-induced vascular change) invalidates it. V1/V2/V3 in section 1.2 would catch that, and they should be re-run before the method is reused.

---

## References

1. Stupp R, Mason WP, van den Bent MJ, et al. Radiotherapy plus concomitant and adjuvant temozolomide for glioblastoma. *N Engl J Med* 2005;352(10):987-96. PMID 15758009. doi:10.1056/NEJMoa043330
2. Gilbert MR, Wang M, Aldape KD, et al. Dose-dense temozolomide for newly diagnosed glioblastoma: a randomized phase III clinical trial. *J Clin Oncol* 2013;31(32):4085-91. PMID 24101040. doi:10.1200/JCO.2013.49.6968
3. Wick W, Steinbach JP, Kuker WM, et al. One week on/one week off: a novel active regimen of temozolomide for recurrent glioblastoma. *Neurology* 2004;62(11):2113-5. PMID 15184628. doi:10.1212/01.wnl.0000127617.89363.84
4. Wick A, Felsberg J, Steinbach JP, et al. Efficacy and tolerability of temozolomide in an alternating weekly regimen in patients with recurrent glioma. *J Clin Oncol* 2007;25(22):3357-61. PMID 17664483. doi:10.1200/JCO.2007.10.7722
5. Weller M, Tabatabai G, Kastner B, et al. MGMT promoter methylation is a strong prognostic biomarker for benefit from dose-intensified temozolomide rechallenge in progressive glioblastoma: the DIRECTOR trial. *Clin Cancer Res* 2015;21(9):2057-64. PMID 25655102. doi:10.1158/1078-0432.CCR-14-2737
6. Perry JR, Belanger K, Mason WP, et al. Phase II trial of continuous dose-intense temozolomide in recurrent malignant glioma: RESCUE study. *J Clin Oncol* 2010;28(12):2051-7. PMID 20308655. doi:10.1200/JCO.2009.26.5520
7. Baker SD, Wirth M, Statkevich P, et al. Absorption, metabolism, and excretion of 14C-temozolomide following oral administration to patients with advanced cancer. *Clin Cancer Res* 1999;5(2):309-17. PMID 10037179
8. Beale P, Judson I, Moore S, et al. Effect of gastric pH on the relative oral bioavailability and pharmacokinetics of temozolomide. *Cancer Chemother Pharmacol* 1999;44(5):389-94. PMID 10501912. doi:10.1007/s002800050994
9. Ostermann S, Csajka C, Buclin T, et al. Plasma and cerebrospinal fluid population pharmacokinetics of temozolomide in malignant glioma patients. *Clin Cancer Res* 2004;10(11):3728-36. PMID 15173079. doi:10.1158/1078-0432.CCR-03-0807
10. Tolcher AW, Gerson SL, Denis L, et al. Marked inactivation of O6-alkylguanine-DNA alkyltransferase activity with protracted temozolomide schedules. *Br J Cancer* 2003;88(7):1004-11. PMID 12671695. doi:10.1038/sj.bjc.6600827
11. TEMODAR (temozolomide) US Prescribing Information, sections 2.1 and 12.3 (dose escalation rule; Cmax 7.5 ug/mL and AUC 23.4 ug*h/mL after 150 mg/m2; elimination half-life 1.8 h; linear kinetics 75-250 mg/m2/day; IV/oral bioequivalence).
12. Himmelstein KJ, Patton TF, Belt RJ, et al. Clinical kinetics of intact cisplatin and some related species. *Clin Pharmacol Ther* 1981;29(5):658-64. PMID 7194166. doi:10.1038/clpt.1981.91
13. Nagai N, Kinoshita M, Ogata H, et al. Relationship between pharmacokinetics of unchanged cisplatin and nephrotoxicity after intravenous infusions of cisplatin to cancer patients. *Cancer Chemother Pharmacol* 1996;39(1-2):131-7. PMID 8995510. doi:10.1007/s002800050548
14. Benjamin RS, Riggs CE, Bachur NR. Plasma pharmacokinetics of adriamycin and its metabolites in humans with normal hepatic and renal function. *Cancer Res* 1977;37(5):1416-20. PMID 856462
15. Camaggi CM, Comparsi R, Strocchi E, et al. Epirubicin and doxorubicin comparative metabolism and pharmacokinetics. *Cancer Chemother Pharmacol* 1988;21(3):221-8. PMID 3162849. doi:10.1007/BF00262774
16. Muller C, Chatelut E, Gualano V, et al. Cellular pharmacokinetics of doxorubicin in patients with chronic lymphocytic leukemia. *Cancer Chemother Pharmacol* 1993;32(5):379-84. PMID 8339389. doi:10.1007/BF00735923
17. Speth PA, van Hoesel QG, Haanen C. Clinical pharmacokinetics of doxorubicin. *Clin Pharmacokinet* 1988;15(1):15-31. PMID 3042244. doi:10.2165/00003088-198815010-00002
18. PubChem CID 5394 (temozolomide, C6H6N6O2, MW 194.15); CID 31703 (doxorubicin, MW 543.5).

---

*Generated by `experiments/regimen.py` — ONE reconstruction reused everywhere, 21 transport solves (`PDEConfig(time_resolved=True)` on every one), 136,080 solver steps, 109.6 min total, dx = 12 um, dt = 10 s, seed 0, deterministic. No file in `tmesim/` was modified.*
