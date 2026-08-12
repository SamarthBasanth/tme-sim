# The jointly calibrated oxygen configuration is grid-labile — **more** so than the shipped one, not less

> **One-line answer.** `results/oxygen_supply_calibration.md` §6.3 left one question open: it
> achieved `F_state` span exactly 0.0000 by correcting `supply` alone, admitted that span was
> degenerate (100 % of cells in one band), and stated that **"the jointly calibrated
> configuration was NOT grid-tested and that is the biggest open question left."** It has now
> been grid-tested. **It is worse.** `F_state` spans **0.1575** across dx = 30/18/12/8 at the
> jointly calibrated low-CMRO₂ rung and **0.0542** at the high one, against the shipped
> **0.0461** and the derived-supply-only **0.0000**. **61.5 %** of tumour cells change
> oxygen band somewhere on that ladder, against 9.9 % under the shipped configuration.
>
> **The band-edge question, answered directly:** jointly calibrated mean tumour O2 does *not*
> park next to an edge in the way the shipped configuration parks next to 0.4. At dx = 12 it
> sits at **0.2615** (low) and **0.1614** (high), i.e. 0.062 and 0.039 from the HYPOXIC/QUIESCENT
> edge at **0.2** — comparable to the shipped configuration's 0.051 from the edge at 0.4. And
> the *absolute* crowding near the edges barely differs: **16.8 %** of cells lie within ±0.02 of
> an edge under joint-lo, against **15.6 %** shipped. **The mean did not get closer to an edge
> and the population did not get more crowded. What changed by an order of magnitude is the
> grid noise on each cell** — median per-cell grid σ of **0.0793** (2.77 mmHg) under joint-lo
> against **0.0085** (0.30 mmHg) shipped, **9.4×**. Lability is `population near edge × grid
> noise`, and the joint calibration multiplies the second factor, not the first.
>
> **The mechanism is structural and was predictable from the constants alone (§4).** Correcting
> the sink by 500× shortens the oxygen diffusion length `L = sqrt(D / (consumption·ρ))` from
> **316 µm to 14 µm** at a dense tile, and correcting the source shortens the source
> equilibration length `sqrt(D/supply)` from **200 µm to 3.2 µm**. Under the shipped constants
> the sink length is **10.5 to 39.5 mesh cells** across this ladder and the source length 6.7 to
> 25 — the field is smooth and every structure in it is resolved. Under the jointly calibrated
> constants the sink length is **0.47 to 1.77 mesh cells** and the source length is **under
> half a cell at every rung**. **The jointly calibrated oxygen field is not under-resolved by
> accident; it is under-resolved because a physically correct vessel is a boundary-layer object
> and this project's meshes cannot represent one.**
>
> **So the answer to the question the task poses — "does the lability move or disappear?" — is
> neither. It grows, and it grows for a different reason than the one that caused it.**
> `results/oxygen_level_shift.md`'s mechanism (a mean parked a hair inside a band edge) is a
> *distribution-position* explanation. The joint configuration's is a *discretisation-error*
> explanation. **That distinguishes the two hypotheses the task set up: the band-edge problem is
> NOT purely a consequence of one wrong constant, and it is NOT purely structural to the banded
> design either. It is structural to banding a field that is not grid-converged** — and fixing
> the constants makes the field less grid-converged, not more.
>
> **The consequence, deliberately last and reported separately in §5:** the 3/3 PASS that
> `oxygen_supply_calibration.md` §6.2 found at dx = 12 does **not** survive the ladder. Cisplatin
> — the drug that report flagged as sitting 3.7 % above its band floor — **FAILs at dx = 30 on
> both joint arms and at dx = 18 on the high arm**, giving **3/4 and 2/4 rungs**. The decisive
> statistic `grid_convergence_gbm.md` §4 defines, grid movement over distance to the nearest band
> edge, is **0.65× shipped** against **1.25× and 4.52×** jointly calibrated.

**Status: measurement of a hypothetical configuration. Nothing was changed.** `OxygenConfig.supply`
is still `0.05`, `consumption` still `0.02`, `hypoxia_threshold` still `0.2`, `REFERENCE_PO2_MMHG`
still `35.0`. `tmesim/` is untouched and `pytest` was therefore not re-run. The only intervention
is `oxygen_supply_calibration._force_cfg`, imported unmodified from the work being extended: a
harness-level wrapper that hands `physics.oxygen.solve_oxygen` a
`dataclasses.replace(cfg, supply=…, consumption=…)`. The solve, the banding
(`cells.agents.update_cell_states`) and the whole scoring path are unmodified library code.

* Generator: `experiments/oxygen_joint_grid_test.py`
* Raw: `results/sweep/oxjoint_repro.json`, `oxjoint_ladder.json`, `oxjoint_bandedge.json`,
  `oxjoint_mechanism.json`, `oxjoint_fineprobe.json`, `oxjoint_verdict.json`,
  `oxjoint_o2_percell.npz` (per-tumour-cell O2 at every rung of every arm),
  `oxjoint_verdict_records.jsonl`
* Section: `data_cache/xenium_gbm`, 40,887 cells / 23,812 tumour, identical at every rung.

---

## 0. Reproduction — the precondition

Nothing below is worth reading unless the jointly calibrated configuration can be re-created
from `results/oxygen_supply_calibration.md` and its generator. It can, exactly. The four
configurations were re-solved at dx = 12 µm through the harness and checked against that
report's published numbers; `experiments/oxygen_joint_grid_test.py` **raises** if any check
fails, so a mis-specified arm cannot propagate.

| configuration | `supply` | `consumption` | mean tumour O2 | median pO₂ (mmHg) | HYPOXIC frac | `F_state` | vs `oxygen_supply_calibration.md` |
|---|---:|---:|---:|---:|---:|---:|---|
| shipped | 0.05 | 0.02 | 0.348628 | 12.107 | 0.0000 | 0.550206 | §6.1: 0.3486 / 0.5502 ✓ |
| derived `supply` only | 194 | 0.02 | 0.992519 | 34.762 | 0.0000 | 1.000000 | §6.1: 0.9925 / 1.0000 ✓ |
| **jointly calibrated, CMRO₂ low** | **194** | **10** | **0.261476** | **8.193** | **0.4022** | **0.524853** | §5.1 raw: exact to 1e-6 ✓ |
| **jointly calibrated, CMRO₂ high** | **194** | **20** | **0.161351** | **4.526** | **0.7076** | **0.416744** | §5.1 raw: exact to 1e-6 ✓ |

**12 of 12 checks pass**, the two joint arms to a relative tolerance of **1e-6** against
`results/sweep/oxsupply_sinkmatch.json`. The scored-ratio harness reproduces that report's
`killcal` arm too — §5 below returns 0.3322 / 0.2769 / 0.2211 at dx = 12 for the low arm
against its published 0.33218 / 0.277 / 0.221.

**Neither `supply` nor `consumption` was re-derived here.** They are taken as established by
`results/oxygen_supply_calibration.md`: `supply ≈ 194 s⁻¹` from Secomb 2020's measured cerebral
`K_TO₂,cap` = 1.213 (PMID [31423930](https://pubmed.ncbi.nlm.nih.gov/31423930/)) via
`supply = K_TO₂,cap/(α·f_v)`; `consumption = 10` and `20` as the two rungs of that report's
sink sweep whose implied CMRO₂ (2.227 and 2.748 cm³O₂·(100 cm³·min)⁻¹) fall inside Watabe 2017's
PET-measured C6-glioma TMRO₂ band of 1.9–2.8 (PMID
[28210996](https://pubmed.ncbi.nlm.nih.gov/28210996/)). The measured human-GBM median those two
rungs bracket is Lally 2006's 5.1 mmHg (PMID
[17207315](https://pubmed.ncbi.nlm.nih.gov/17207315/)). **All of that is inherited, including
its UNVERIFIED flags; none of it is re-litigated here.**

---

## 1. The ladder — numerical behaviour, reported before any verdict

dx = 30 / 18 / 12 / 8 µm, matching `results/grid_convergence_gbm.md`'s rungs (18/12/8 is its
constant-ratio r = 1.5 Richardson triple) so the two ladders are directly comparable. **The
oxygen solve is steady state — it carries no CFL condition and no `dt`** — so §1–§4 are
dt-independent by construction. `dt` enters only in §5, where it is stated per rung.

### 1.1 Jointly calibrated, CMRO₂ low (`supply` 194, `consumption` 10)

| dx (µm) | mean tumour O2 | mmHg | median mmHg | p10 / p90 mmHg | HYPOXIC | QUIESCENT | VIABLE | **`F_state`** |
|---:|---:|---:|---:|---|---:|---:|---:|---:|
| 30 | 0.11500 | 4.03 | 2.19 | 0.32 / 9.81 | 0.8267 | 0.1381 | 0.0352 | **0.39357** |
| 18 | 0.20027 | 7.01 | 5.70 | 1.79 / 13.82 | 0.6015 | 0.3026 | 0.0959 | **0.45771** |
| **12** | **0.26148** | **9.15** | **8.19** | 3.22 / 16.23 | 0.4022 | 0.4274 | 0.1704 | **0.52485** |
| 8 | 0.28055 | 9.82 | 9.00 | 3.71 / 17.08 | 0.3469 | 0.4470 | 0.2061 | **0.55102** |

### 1.2 Jointly calibrated, CMRO₂ high (`supply` 194, `consumption` 20)

| dx (µm) | mean tumour O2 | mmHg | median mmHg | p10 / p90 mmHg | HYPOXIC | QUIESCENT | VIABLE | **`F_state`** |
|---:|---:|---:|---:|---|---:|---:|---:|---:|
| 30 | 0.07238 | 2.53 | 0.89 | 0.07 / 6.45 | 0.9131 | 0.0660 | 0.0209 | **0.37350** |
| 18 | 0.12210 | 4.27 | 2.90 | 0.61 / 9.65 | 0.8156 | 0.1577 | 0.0268 | **0.39106** |
| **12** | **0.16135** | **5.65** | **4.53** | 1.30 / 11.34 | 0.7076 | 0.2466 | 0.0458 | **0.41674** |
| 8 | 0.17385 | 6.09 | 5.08 | 1.57 / 11.99 | 0.6609 | 0.2854 | 0.0536 | **0.42768** |

### 1.3 The two control arms, same ladder

| dx (µm) | shipped: mean O2 | shipped `F_state` | derived-only: mean O2 | derived-only `F_state` |
|---:|---:|---:|---:|---:|
| 30 | 0.33636 | 0.50412 | 0.97280 | 1.00000 |
| 18 | 0.34817 | 0.54015 | 0.98884 | 1.00000 |
| 12 | 0.34863 | 0.55021 | 0.99252 | 1.00000 |
| 8 | 0.34824 | 0.54882 | 0.99329 | 1.00000 |

The shipped column reproduces `results/grid_convergence_gbm.md` §5 (0.5041 / 0.5401 / 0.5502 /
0.5488) to four decimals — that agreement is the control that this harness computes the same
`F_state` the published ladder computes, and it is why the joint arms can be read against it.

### 1.4 The headline that matters most on this ladder

**The measured pO₂ agreement that motivated the joint calibration is itself a function of the
mesh.** `oxygen_supply_calibration.md` §5.1's central claim is that the two Watabe-matched sink
rungs bracket Lally's measured 5.1 mmHg with nothing fitted to it. Across this ladder the
bracket moves:

| dx (µm) | median pO₂, CMRO₂ low | median pO₂, CMRO₂ high | **does the pair bracket 5.1 mmHg?** |
|---:|---:|---:|---|
| 30 | 2.19 | 0.89 | **NO — both below the measurement** |
| 18 | 5.70 | 2.90 | yes, but narrowly from beneath |
| **12 (the published rung)** | **8.19** | **4.53** | **yes** |
| 8 | 9.00 | 5.08 | yes, but the high rung is now **0.4 % below** 5.1 |

At dx = 30 the jointly calibrated GBM is a median 0.9–2.2 mmHg tumour — hypoxic well past
anything Lally measured. At dx = 8 it is 5.1–9.0 mmHg. **The three-way agreement between
Secomb, Watabe and Lally that `oxygen_supply_calibration.md` §5.1 reports holds at dx = 18, 12
and 8 and fails outright at dx = 30, and the bracket's lower rail has drifted onto the
measurement itself by dx = 8.** That is not a criticism of the derivation — it is the same
discretisation error §4 explains — but it means the agreement is quoted at one resolution, the
whole bracket slides upward under refinement, and how much further it slides below dx = 8 was
not measured for the median (only for `F_state`, §2.1, where the movement continues to dx ≈ 5).

---

## 2. `F_state` span and monotonicity — against all three known configurations

| configuration | `F_state`, dx = 30 → 8 | **span** | span, relative | monotone in dx? |
|---|---|---:|---:|---|
| shipped (0.05, 0.02) | 0.5041, 0.5401, 0.5502, 0.5488 | **0.0461** | 9.1 % | **NO** (turns between 12 and 8) |
| derived `supply` only (194, 0.02) | 1.0000, 1.0000, 1.0000, 1.0000 | **0.0000** | 0.0 % | trivially yes — **degenerate** |
| **jointly calibrated, CMRO₂ low (194, 10)** | 0.3936, 0.4577, 0.5249, 0.5510 | **0.1575** | **40.0 %** | **yes, monotone increasing** |
| **jointly calibrated, CMRO₂ high (194, 20)** | 0.3735, 0.3911, 0.4167, 0.4277 | **0.0542** | 14.5 % | **yes, monotone increasing** |

**Two things are true at once and they have to be said together.**

1. **The joint arms are monotone where the shipped arm is not.** That is a genuine improvement
   in *shape*: a monotone ladder is at least consistent with an asymptotic convergence regime,
   which `grid_convergence_gbm.md` §3 says the shipped ladder is not (it calls the non-monotone
   ladder the reason its Richardson and power-law extrapolations are outside their assumptions).
2. **The joint arms move 1.2× to 3.4× further doing it.** Monotone but three times as large is
   not better for any use this project makes of `F_state`. The shipped 0.0461 is oscillation
   inside a plateau; the joint-lo 0.1575 is a curve that has not arrived.

### 2.1 Has it arrived by dx = 8? — the fine probe

The ladder leaves `F_state` still climbing at its finest rung, and "still moving at dx = 8" and
"converging to a limit just below dx = 8" are different findings. Extra oxygen-only rungs (no
transport, no `dt`, no CFL) settle it:

| dx (µm) | shipped `F_state` | shipped Δ | joint-lo mean O2 | **joint-lo `F_state`** | **joint-lo Δ** | joint-hi `F_state` |
|---:|---:|---:|---:|---:|---:|---:|
| 30 | 0.50412 | — | 0.11500 | 0.39357 | — | 0.37350 |
| 18 | 0.54015 | +0.03603 | 0.20027 | 0.45771 | +0.06414 | 0.39106 |
| 12 | 0.55021 | +0.01006 | 0.26148 | 0.52485 | +0.06714 | 0.41674 |
| 8 | 0.54882 | −0.00139 | 0.28055 | 0.55102 | +0.02617 | 0.42768 |
| 6 | 0.53901 | −0.00981 | 0.28647 | 0.55862 | +0.00760 | 0.42959 |
| 5 | 0.54351 | +0.00450 | 0.28785 | 0.56132 | +0.00270 | 0.42938 |
| 4 | 0.54267 | −0.00084 | 0.28795 | **0.56169** | **+0.00037** | **0.42945** |

**The joint-lo sequence converges — the step size runs 0.064 → 0.067 → 0.026 → 0.008 → 0.003 →
0.0004, falling monotonically from dx = 18 downwards, and it converges to `F_state` ≈ 0.5617 at
dx = 4–5 µm.** So the
jointly calibrated `F_state` at the grid this project publishes on (dx = 12, 0.5249) is
**6.6 % below its own converged value**, and its dx = 30 value (0.3936) is **30 % below** it.
The high-CMRO₂ arm converges the same way, to 0.4294, leaving dx = 12 **3.0 % low**.
The whole 0.1575 ladder span is real discretisation error, not noise. The shipped arm, by
contrast, oscillates inside a ±0.005 band from dx = 18 all the way down to dx = 4 and never
settles — consistent with `grid_convergence_gbm.md`'s "plateau with oscillation" reading, and
the reason its span cannot be interpreted as an error estimate.

**This is the cleanest single contrast in the file.** The shipped configuration has small,
non-convergent grid noise. The jointly calibrated configuration has large, genuinely convergent
grid error. Neither is grid-converged at dx = 12; they fail to be converged in opposite ways.

---

## 3. The band-edge question, answered directly

### 3.1 Where the mean sits

| configuration | mean tumour O2 (dx = 12) | distance to edge 0.2 | distance to edge 0.4 | nearer edge |
|---|---:|---:|---:|---|
| shipped | 0.3486 | +0.1486 | **−0.0514** | **0.4 (QUIESCENT/VIABLE)** |
| derived `supply` only | 0.9925 | +0.7925 | +0.5925 | 0.4, but 0.59 away |
| **jointly calibrated, low** | **0.2615** | **+0.0615** | −0.1385 | **0.2 (HYPOXIC/QUIESCENT)** |
| **jointly calibrated, high** | **0.1614** | **−0.0386** | −0.2386 | **0.2 (HYPOXIC/QUIESCENT)** |

**So the task's hypothesis is half right.** The jointly calibrated mean does move off the 0.4
edge and onto the 0.2 edge — but it does not park *appreciably closer* to it. The low arm sits
**0.0615** from its nearest edge against the shipped configuration's **0.0514**, i.e. slightly
*further* away; the high arm sits **0.0386**, slightly nearer. On the mean alone the three
configurations are indistinguishable in how exposed they are, and the shipped one is in the
middle of them. Yet their band-change counts differ by a factor of six (§3.4). **The mean is
the wrong statistic and this is the demonstration of why.**

### 3.2 How much of the population is within grid-noise distance of an edge

The grid-noise scale is **measured, not assumed**. This section carries measured single-cell
coordinates and `reconstruct_real`'s tissue mask is the dx-independent `'geometric'` one, so the
cell population is identical at every rung (40,887 / 23,812 — asserted in
`grid_convergence_gbm.md` §2.1). Cell *i* at dx = 30 is therefore the same cell as cell *i* at
dx = 8, and the spread of its own four O2 values is grid error and nothing else. `σ_i` is that
per-cell sample standard deviation across the four rungs.

| configuration | median per-cell σ | in mmHg | **within 1σ of an edge** | within 2σ | within 3σ |
|---|---:|---:|---:|---:|---:|
| shipped | 0.00846 | 0.30 | **12.18 %** (all at 0.4) | 17.13 % | 20.77 % |
| derived `supply` only | 0.00860 | 0.30 | **0.00 %** | 0.00 % | 0.00 % |
| **jointly calibrated, low** | **0.07926** | **2.77** | **66.45 %** (45.8 % at 0.2, 32.0 % at 0.4) | 81.74 % | 87.98 % |
| **jointly calibrated, high** | **0.04964** | **1.74** | **41.97 %** (35.9 % at 0.2, 13.9 % at 0.4) | 56.55 % | 65.01 % |

**Two-thirds of the tumour population is within one grid-σ of a band edge under the jointly
calibrated low-CMRO₂ configuration**, against one-eighth under the shipped one. That, and not
the position of the mean, is the answer to the task's question.

### 3.3 Why — it is the noise, not the crowding

The σ-scaled measure conflates two things, so here is the same population measured against a
**fixed** window that has no noise model in it:

| configuration | within ±0.005 of an edge | ±0.01 | ±0.02 | ±0.05 | **median per-cell grid σ** |
|---|---:|---:|---:|---:|---:|
| shipped | 4.86 % | 9.87 % | **15.63 %** | 41.79 % | **0.00846** |
| **jointly calibrated, low** | 4.38 % | 8.34 % | **16.78 %** | 42.29 % | **0.07926** |
| jointly calibrated, high | 2.72 % | 5.31 % | 11.00 % | 28.42 % | 0.04964 |

> **The jointly calibrated distribution is no more crowded against the band edges than the
> shipped one — 16.8 % vs 15.6 % within ±0.02, 42.3 % vs 41.8 % within ±0.05. The entire
> difference in lability is the 9.4× increase in per-cell grid noise.**

This is the finding that separates the two competing explanations the task named. It is **not**
that the joint configuration parks next to a different edge (§3.1 shows it parks slightly
further away). It is **not** that banding is inherently fragile for any field (the shipped
field is banded by the same rule with 12 % near-edge population and moves 3× less). It is that
**the jointly calibrated oxygen field is an order of magnitude less grid-converged**, and
banding is a step-function operator, so it converts a field's discretisation error into a
population's relabelling at a rate set by the local density of cells at the edge.

Measured, that density (cells per unit O2, in a ±0.01 window, as a fraction of the tumour
population):

| configuration | density at edge 0.2 | density at edge 0.4 |
|---|---:|---:|
| shipped | 0.000 | **4.937** |
| jointly calibrated, low | **2.772** | 1.398 |
| jointly calibrated, high | **2.203** | 0.454 |

Under the shipped configuration one edge is loaded and the other is empty — exactly the picture
`results/oxygen_level_shift.md` §6 describes. Under the joint configuration **both edges are
loaded**, and the total is comparable. So the joint calibration does not remove the mechanism;
it spreads it over two edges instead of one, and then multiplies the driver.

### 3.4 The direct measurement — cells that actually change band

No noise model at all: count the cells whose band index differs between two rungs.

| configuration | dx 30↔18 | dx 18↔12 | **dx 12↔8** | **changed band somewhere on the ladder** |
|---|---:|---:|---:|---:|
| shipped | 1,878 (7.9 %) | 565 (2.4 %) | **122 (0.5 %)** | **2,349 (9.86 %)** |
| derived `supply` only | 0 | 0 | **0** | **0 (0.00 %)** |
| **jointly calibrated, low** | 8,431 (35.4 %) | 7,740 (32.5 %) | **3,365 (14.1 %)** | **14,650 (61.52 %)** |
| **jointly calibrated, high** | 4,436 (18.6 %) | 4,559 (19.2 %) | **2,673 (11.2 %)** | **8,424 (35.38 %)** |

**The single refinement from the published dx = 12 µm to dx = 8 µm relabels 14.1 % of the
tumour population under the jointly calibrated configuration, against 0.5 % under the shipped
one — a factor of 27.6.** Under the shipped configuration the fine end of the ladder has
essentially stopped relabelling cells; under the joint configuration it has not.

---

## 4. Why — the mechanism, from the constants alone

This needed no ladder to predict, and it is the same structural argument
`grid_convergence_gbm.md` §2 makes for the *drug* field, applied to the *oxygen* field. Two
lengths govern `D·∇²O + supply·vasc·(1−O) − consumption·ρ·O = 0`:

* **source equilibration length** `sqrt(D/supply)` — how far from a vessel pixel the source acts;
* **sink diffusion length** `L = sqrt(D/(consumption·ρ))` — how far oxygen travels before it is consumed.

With `D = 2000 µm²/s` (shipped `OxygenConfig.diffusivity`) and the section's own measured tumour
density (tissue-mean ρ = 0.1229, p90 = 1.0, max = 4.0 cells per 144 µm² tile):

| configuration | source length | sink `L` at mean ρ | **sink `L` at p90 ρ** | **mesh cells per `L`, dx = 30 / 18 / 12 / 8** |
|---|---:|---:|---:|---|
| shipped (0.05, 0.02) | **200.0 µm** | 902 µm | **316.2 µm** | **10.5 / 17.6 / 26.4 / 39.5** |
| derived `supply` only (194, 0.02) | **3.2 µm** | 902 µm | 316.2 µm | 10.5 / 17.6 / 26.4 / 39.5 |
| **jointly calibrated, low (194, 10)** | **3.2 µm** | 40.3 µm | **14.1 µm** | **0.47 / 0.79 / 1.18 / 1.77** |
| **jointly calibrated, high (194, 20)** | **3.2 µm** | 28.5 µm | **10.0 µm** | **0.33 / 0.56 / 0.83 / 1.25** |

> **Under the shipped constants the oxygen field's shortest structural length is 10–40 mesh
> cells and the field is smooth. Under the jointly calibrated constants it is 0.3–1.8 mesh
> cells — below one cell at three of the four rungs — and the field is a boundary layer the
> mesh cannot represent.**

And this is not a fixable oversight in the calibration. It is what the calibration *means*.
`oxygen_supply_calibration.md` §1 established that the shipped `supply` gives a source
equilibration length of 200 µm against a 7.5 µm vessel radius, and called that the signature of
the error: *"the model's vessels are not vessels, they are weak volumetric sources sitting in a
bath."* Correcting `supply` to 194 s⁻¹ makes the length **3.2 µm — smaller than the vessel
radius, which is exactly the point**: the vessel now pins its own oxygen, as a real vessel does.
But a source that equilibrates within 3.2 µm and a sink that consumes within 10–14 µm together
describe a field whose entire structure lives inside one 12 µm grid cell.

**Correcting the physics and resolving the physics are the same problem, and this project has
only done the first half.** The honest statement of §1–§3 is therefore not "the joint
calibration is grid-labile" but **"a physically correct oxygen field on this section requires a
mesh several times finer than any this project runs, and the jointly calibrated numbers at
dx = 12 µm are a coarse-grid rendering of it."** The fine probe in §2.1 supports that: at
dx = 4–5 µm — about 3 cells per sink length — `F_state` finally settles.

---

# 5. Consequence — the GBM 3-D:2-D verdict, reported separately and last

Everything above is a statement about an oxygen field and a discretisation. Nothing in it was
chosen with reference to anything in this section, and §1–§4 should be readable and disagreeable
without this section existing.

24 transport solves through `kill_calibration`'s own unmodified scoring path (imported via
`experiments/grid_convergence_gbm.py`: `prepare_tissue` → `unit_transport` → `decompose` →
`ratio_test`), `alpha = 148.6082` unchanged, the literature bands frozen and **identical at
every rung** — they depend only on the frozen `LIT_2D3D` table and on (IC50, hill, c_bath),
none of which dx or the oxygen configuration touches. **A verdict can therefore only move
because the ratio moved.** Exposure fingerprints are asserted distinct across the ladder, so dx
is not a no-op.

**`dt`**: 10 s (the published value) at dx = 30 / 18 / 12; **`dt` = 5 s at dx = 8**, stated
because whether `dt` = 10 s is admissible for the discrete oxygen field at dx ≤ 8 is an open
question owned by another workstream right now and this ladder must not inherit it.
`grid_convergence_gbm.md` §7.3 measured that halving `dt` at dx = 8 moves the scored ratio by
0.0003 %, so the rung remains comparable with its published counterpart. Measured CFL at
dx = 8 / `dt` = 5 was **0.470**.

**Cross-check before anything is read off this section**: at dx = 12 this harness returns
temozolomide **0.3322 / 0.3424**, cisplatin **0.2769 / 0.2164**, doxorubicin **0.2211 / 0.1699**
for the low / high arms — against `results/sweep/oxsupply_killcal.json`'s 0.33218 / 0.34240,
0.277 / 0.21644, 0.221 / 0.170. **It is the same pipeline producing the same numbers.**

### 5.1 The ladder — CMRO₂ low (`supply` 194, `consumption` 10)

| drug | band | dx = 30 | dx = 18 | dx = 12 | dx = 8 (dt 5) | span | monotone | **PASS** |
|---|---|---|---|---|---|---:|---|---|
| temozolomide | 0.1590–0.4622 | 0.3476 PASS | 0.3589 PASS | 0.3322 PASS | 0.3278 PASS | 0.0311 | no | **4/4** |
| **cisplatin** | 0.2086–1.5004 | **0.1951 FAIL** | 0.2485 PASS | 0.2769 PASS | 0.2805 PASS | 0.0855 | yes | **3/4** |
| doxorubicin | 0.0575–1.0000 | 0.1328 PASS | 0.1931 PASS | 0.2211 PASS | 0.2290 PASS | 0.0962 | yes | **4/4** |

Distance above the band floor: cisplatin **−6.5 % / +19.1 % / +32.7 % / +34.5 %**.

### 5.2 The ladder — CMRO₂ high (`supply` 194, `consumption` 20)

| drug | band | dx = 30 | dx = 18 | dx = 12 | dx = 8 (dt 5) | span | monotone | **PASS** |
|---|---|---|---|---|---|---:|---|---|
| temozolomide | 0.1590–0.4622 | 0.3575 PASS | 0.3722 PASS | 0.3424 PASS | 0.3320 PASS | 0.0402 | no | **4/4** |
| **cisplatin** | 0.2086–1.5004 | **0.1810 FAIL** | **0.2078 FAIL** | 0.2164 PASS | 0.2150 PASS | 0.0355 | no | **2/4** |
| doxorubicin | 0.0575–1.0000 | 0.1199 PASS | 0.1569 PASS | 0.1699 PASS | 0.1731 PASS | 0.0532 | yes | **4/4** |

Distance above the band floor: cisplatin **−13.2 % / −0.4 % / +3.8 % / +3.1 %**.

### 5.3 What that says

**1. Cisplatin flips, exactly as `oxygen_supply_calibration.md` §6.2 predicted it would.** That
report wrote of the high-CMRO₂ arm: *"it gets closer to its band edge on cisplatin: at
`consumption = 20` the predicted 0.21644 sits 3.7 % above the band floor of 0.20864. A modestly
different `consumption` would fail cisplatin instead of temozolomide."* **It did not take a
different `consumption`. It took a different mesh.** The high arm's cisplatin is **FAIL at
dx = 30 and FAIL at dx = 18** — and the dx = 18 miss is by **0.4 %**, 0.2078 against a floor of
0.2086, which is smaller than the 0.65 % the ratio moves between dx = 12 and dx = 8. **The
published 3/3 at dx = 12 is 2 of 4 rungs, and the rung immediately coarser than the published
one is a hair's breadth FAIL.**

**2. The low arm is more robust but not clean either** — cisplatin FAILs at dx = 30, so 3/4.

**3. Neither joint arm reproduces `grid_convergence_gbm.md`'s clean result, and the reason is
`F_state`, not `F_penetration`.** `F_penetration` is essentially unchanged between the shipped
configuration and both joint arms at every rung (dx = 12: 0.4678 / 0.4980 / 0.3688, identical to
the published shipped values to four decimals — the oxygen configuration does not touch drug
transport). **Every difference in this section comes through `F_state`**, which is precisely the
quantity §2 showed spans 0.1575:

| | shipped | joint lo | joint hi |
|---|---|---|---|
| `F_state` at dx = 30 → 8 (cisplatin/doxorubicin) | 0.5041 → 0.5488 | **0.3936 → 0.5510** | **0.3735 → 0.4277** |
| cisplatin ratio, dx = 30 → 8 | 0.2358 → 0.2689 | 0.1951 → 0.2805 | 0.1810 → 0.2150 |
| **worst-drug grid movement / distance to nearest band edge** | **0.65×** (`grid_convergence_gbm.md` §4) | **1.25×** | **4.52×** |

> **`grid_convergence_gbm.md` §4 identified the single number that decides whether a verdict is
> grid-fragile: the grid-induced movement of the ratio divided by the distance to the nearest
> band edge. It measured a worst case of 0.65× on the shipped configuration and concluded "for
> every drug the entire grid-induced movement is smaller than the distance to the nearest band
> edge… that, and not the smoothness of the curve, is why the GBM 3/3 PASS survives." Under the
> jointly calibrated configuration that number is 1.25× and 4.52× — both above 1, i.e. the mesh
> alone can move the ratio across a band edge, and on the high arm it can do so four times over.
> The condition that protected the published verdict is gone, and cisplatin's two FAILs are that
> condition failing in practice.**

**4. A structural side-effect worth recording: temozolomide's per-drug hypoxia override becomes
load-bearing for the first time.** `PDConfig.hypoxia_resistance` is 0.35 globally, but
temozolomide carries a per-drug override of **0.80**, derived in `core/config.py` from
Musah-Eroje & Watson 2019 (PMID [30694423](https://pubmed.ncbi.nlm.nih.gov/30694423/)) and
Wang et al. 2017 (PMID [28801626](https://pubmed.ncbi.nlm.nih.gov/28801626/)). On the shipped
GBM field the HYPOXIC band is **empty**, so that override is numerically inert and all three
drugs share `F_state` = 0.5502 (visible in `grid_convergence_gbm.md` §3, where the column is
identical across drugs). Under the jointly calibrated field 35–91 % of tumour cells are HYPOXIC
and the override splits the drugs apart: at dx = 12, temozolomide runs at `F_state` **0.7059 /
0.7352** against **0.5249 / 0.4167** for cisplatin and doxorubicin — a **+34 % / +76 %** relative
difference. **Temozolomide's ratios in §5.1–5.2 are held up by a per-drug constant that has
never before affected a published GBM number.**

How much does it matter? The decomposition is `ratio ≈ F_penetration × F_state × interaction`
with `interaction` measured at 0.992–1.009 on every temozolomide rung, so the counterfactual is a
clean linear rescale. Stripping the override (i.e. running temozolomide at the global 0.35 like
the other two drugs) would give:

| | dx = 30 | dx = 18 | dx = 12 | dx = 8 |
|---|---:|---:|---:|---:|
| temozolomide as scored, low arm | 0.3476 (+118.6 % over floor) | 0.3589 | 0.3322 | 0.3278 (+106.1 %) |
| **without the override, low arm** | **0.1787 (+12.4 %)** | 0.2255 (+41.8 %) | 0.2470 (+55.3 %) | 0.2554 (+60.6 %) |
| **without the override, high arm** | **0.1702 (+7.1 %)** | 0.1920 (+20.8 %) | 0.1941 (+22.1 %) | 0.1958 (+23.2 %) |

**Temozolomide would still be 4/4 without the override — but its worst rung's margin over the
band floor collapses from +106 % to +7 %.** So the override is not what makes temozolomide pass; it is
what makes temozolomide pass *comfortably*, and it is doing that work for the first time in this
project's history.

**5. The direction of the whole effect.** Both joint arms *raise* temozolomide's ratio well
above the shipped 0.2626 (to 0.33–0.37) and *lower* cisplatin's and doxorubicin's. So the
jointly calibrated configuration does not shift the panel uniformly — it re-orders it, because
temozolomide's hypoxia multiplier is 2.3× the other two drugs' and the hypoxic fraction went
from 0 to 0.35–0.91.

**6. Read against §1–§4, none of this is a verdict about the drugs.** The ratio moves because
`F_state` moves; `F_state` moves because the oxygen field is under-resolved by a factor of
several in mesh spacing (§4); and that under-resolution is a direct consequence of the constants
being physically correct. **The honest summary of §5 is not "the joint calibration fails
cisplatin" but "the joint calibration's verdicts are not measurable on the meshes this project
can currently run."**

---

## 6. Limitations

1. **One section, one patient, one tumour type.** Every number is the `xenium_gbm` section.
2. **`consumption = 10–20 s⁻¹` is `oxygen_supply_calibration.md`'s extrapolation, not an adopted
   value**, and is 500–1000× the shipped default and 20–80× above the higher of the two brackets
   the previous OCR pass derived. Its limitation 4 carries over in full.
3. **The sink is first-order in O throughout.** `results/oxygen_sensitivity_band.md` and Grimes
   et al. 2014 (PMID [26064525](https://pubmed.ncbi.nlm.nih.gov/26064525/)) both establish that
   is the wrong functional form. A Michaelis–Menten sink saturates at low O2 and would give a
   *longer* effective diffusion length in exactly the hypoxic regions where §4's `L` is
   shortest, so it would plausibly reduce — not remove — the grid sensitivity measured here.
   **This is the single most likely way §4's conclusion could be softened, and it was not
   tested.**
4. **The 7.5 µm vessel-painting radius is known-wrong** (`results/vascular_segmentation_audit.md`)
   and `supply = 194` inherits that error linearly through `f_v`. `oxygen_supply_calibration.md`
   §7 says fix the geometry first; this file did not, because its job was to grid-test the
   configuration that report actually specified.
5. **The 35 mmHg anchor is used throughout**, and `oxygen_supply_calibration.md` §2 argues the
   compartment-correct anchor for `O = 1` is ~55–60 mmHg. Every mmHg figure here would scale by
   ~1.7× under that re-anchoring, and the band edges would move to 11.75 / 23.5 mmHg — which
   would change §3 substantially, since it is a statement about where a distribution sits
   relative to fixed edges. **§3 is conditional on the current anchor.**
6. **The per-cell grid σ of §3.2 is estimated from four rungs.** It is a spread across a
   deliberately wide dx range (30 → 8), not a local discretisation-error estimate, and it is
   inflated by the coarse rung. The fixed-window measure of §3.3 and the direct band-change
   counts of §3.4 carry no such dependence and reach the same conclusion, which is why all
   three are reported.
7. **The ladder is 4 spatial rungs, and `dt` is not converged.** `dt` = 10 s at dx = 30/18/12
   and 5 s at dx = 8 (§5), so the finest verdict rung differs from the others in two ways, not
   one. `grid_convergence_gbm.md` §7.3 bounds that at 0.0003 % of the scored ratio **on the
   shipped configuration**; it was not re-measured under the jointly calibrated one, where the
   oxygen field (and hence `F_state`, and hence the score) is different. The oxygen results of
   §1–§4 are unaffected — that solve is steady state and carries no `dt` at all.
8. **No richer convergence analysis was attempted on the verdict ladder.** Four rungs with a
   non-monotone drug (temozolomide on both arms, cisplatin on the high arm) do not support
   Richardson extrapolation, and none was performed. The `F_state` fine probe (§2.1) does reach
   a converged value, but the *scored ratio* at dx = 4–6 µm was not measured — that would be
   another 18 transport solves.
9. **This file grid-tests a configuration; it does not endorse one.** Whether the jointly
   calibrated constants are right is `oxygen_supply_calibration.md`'s question, and its own
   §5 caveats stand unmodified. **Nothing here argues for keeping `supply = 0.05`** — §4 of that
   report showed the shipped value is wrong by ~3.7 orders of magnitude with an internal
   signature that needs no external data. What this file adds is that the correction cannot be
   evaluated on a 12 µm mesh, so **`results/vessel_geometry_fix.md`'s sub-grid vessel path and a
   mesh-convergence plan are prerequisites for the joint re-derivation, not follow-ups to it** —
   which sharpens step 6 of that report's own recommended ordering.

---

## 7. Reproduce

```
./.venv/bin/python experiments/oxygen_joint_grid_test.py repro ladder bandedge \
                   mechanism fineprobe verdict
```

`repro` **raises** if the jointly calibrated configuration does not reproduce
`results/oxygen_supply_calibration.md`'s published numbers, so a mis-specified arm cannot reach
any table. Each stage checkpoints to `results/sweep/oxjoint_<stage>.json`; `verdict` additionally
appends one JSON line per transport solve to `results/sweep/oxjoint_verdict_records.jsonl` and is
restart-safe.
