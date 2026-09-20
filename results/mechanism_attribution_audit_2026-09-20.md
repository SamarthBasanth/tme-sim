# Mechanism attribution: the decision procedure, the provenance of every threshold, and how far a cell sits from flipping

**2026-09-20.** An audit of `tmesim/sim/mechanism.py` — the module that produces the product's
headline sentence ("of the tumour cells that survived, N were DELIVERY-limited / DOSE-limited /
RESISTANCE-limited / limited by CHANCE") and the labels a click on a finding isolates.

**Nothing was changed.** No parameter was moved, no threshold retuned, no engine file edited.
Three guards were added to `tests/test_mechanism.py` and each was plant-tested (§7). The
measurement ran the shipped path (`simulate_response`, `transvascular_exchange = False`) and
read the kill model's own arrays by wrapping `apply_drug_kill` at the engine's import site
inside the measurement process; no repository file was touched to obtain them, and the
re-derivation of `ratio` and `p_kill` from the captured `E` was verified **bit-identical** to
the engine's own arrays before anything was computed on it (§8).

`calibrate_alpha()` asserted live at `22.223477711172563` at the START and the END of both
measurement runs (§8).

---

## Headline

| finding | number |
|---|---|
| Tumour cells on the flagship section at distance **exactly zero** from a resistance-axis boundary | **23,812 of 23,812 — all of them** |
| Their modulation spectrum | two values only: **0.5** (22,074 cells) and **1.0** (1,738 cells) |
| `resistant_modulation`, derived | **0.5** — bit-identical to the value 92.7 % of the cells carry |
| Relative perturbation of modulation that **inverts the dominant limiter** for doxorubicin and cisplatin | **1 part in 10¹²** (resistance → chance, 97.5 % / 98.4 % of survivors re-labelled) |
| Survivors at exactly zero distance — doxorubicin / cisplatin / temozolomide | **98.55 % / 99.64 % / 0.00 %** |
| `unreached_ratio` at which temozolomide's headline inverts DOSE → DELIVERY | **0.2918** (shipped value 0.1, a convention) |
| Cells the DELIVERY band holds for cisplatin (and all but 2 for doxorubicin) | the **140** cells that get essentially no drug — flat across every `unreached_ratio` from 0.01 to 0.9 |
| Temozolomide "DOSE-limited" survivors that are **also** below the resistance threshold | **19,921 of 21,268 = 93.67 %** |
| Does a **bit-exact reconciliation** constrain the split at all? | **No.** Randomising all four survival labels reconciles bit-exactly on all three real runs |
| Cisplatin's margin from being the 8th of 12 inverting compounds | **0.22 %** of exposure (would make the census 66.7 %, not 58.3 %) |

**The verdict in one sentence: the exposure axis of this classifier is robust and the biology
axis is not a measurement at all.** On the flagship section every single tumour cell sits on an
exact tie of the resistance test, resolved not by data but by the direction of an inclusive
inequality, and the threshold it is compared against is derived from the same constant that
generated the value being tested. Where the verdict is RESISTANCE — which is 2 of the 3
compounds measured here and 9 of 13 in `results/verdict_inversion_census.md` — it is not a
finding about the tissue.

> **Not a compound comparison.** Three compounds appear below because the classifier behaves
> differently in different exposure regimes and one compound cannot show that. Nothing here
> ranks them, orders them, or says anything about which would work. Cross-compound ranking was
> withdrawn 2026-08-21 and stays withdrawn.

> **`kill_fraction` appears below only as the quantity the reconciliation is checked against.**
> It is not a result here and absolute kill is not a defensible output of this project.

---

## 1. The decision procedure, exactly as coded

`classify()` (`tmesim/sim/mechanism.py`) reads five arrays the kill model filled
(`alive`, `killed`, `ratio`, `modulation`, `p_kill`) and writes one `uint8` per cell. Assignment
is by **successive overwrite in source order**, so the order is the precedence. `ratio = E/IC50`
and `modulation = sensitivity x resistance`, both taken from `apply_drug_kill`'s own objects.

| # | test, in the order it is applied | quantity | threshold | label if true |
|---|---|---|---|---|
| 0 | `cell_type != TUMOR` | — | — | `NOT_TUMOR` (7) — the initial fill |
| 1 | `tumor & ~alive` | pre-drug state | — | `ALREADY_DEAD` (6) |
| 2 | `tumor & killed` | the engine's own mask | — | `KILLED` (0) |
| 3 | survivor **and** any of `ratio`, `modulation`, `p_kill` non-finite | — | — | `UNATTRIBUTED` (5) |
| 4 | `ratio < unreached_ratio` | `E/IC50` | `0.1` | `SURVIVED_UNREACHED` (1) → **delivery** |
| 5 | `ratio < adequate_ratio` | `E/IC50` | `1.0` | `SURVIVED_UNDERDOSED` (2) → **dose** |
| 6 | `modulation <= resistant_modulation` **and** `modulation < 1.0` | `sensitivity x resistance` | `0.5` derived | `SURVIVED_RESISTANT` (3) → **resistance** |
| 7 | `p_kill >= stochastic_min_p_kill` | `p_kill` | `0.25` derived | `SURVIVED_STOCHASTIC` (4) → **chance** |
| 8 | anything left | — | — | `UNATTRIBUTED` (5) |

`survived = tumor & alive & ~killed`; a cell can only reach tests 4–8 through it. Tests 4 and 5
are **two bands on one axis** (`E/IC50`), which the module states in every payload
(`exposure_axis_note`). Tests 6 and 7 are reached **only** by cells that passed test 5, i.e.
only by cells with `E >= IC50`.

### 1.1 What happens at a tie, measured at 1-ULP resolution

Every boundary was probed at the exact threshold and one unit-in-the-last-place either side,
through the shipped `classify()` at the production `PDConfig` and `hill = 1.2`:

| cell | label |
|---|---|
| `ratio` exactly `0.1` | `SURVIVED_UNDERDOSED` |
| `ratio` = `0.1 - 1 ulp` | `SURVIVED_UNREACHED` |
| `ratio` exactly `1.0` | `SURVIVED_STOCHASTIC` |
| `ratio` = `1.0 - 1 ulp` | `SURVIVED_UNDERDOSED` |
| `modulation` exactly `0.5` (a QUIESCENT cell) | **`SURVIVED_RESISTANT`** |
| `modulation` = `0.5 + 1 ulp` | **`SURVIVED_STOCHASTIC`** |
| `modulation` exactly `1.0` (a VIABLE cell) | `SURVIVED_STOCHASTIC` |
| `modulation` = `1.0 - 1 ulp` | `SURVIVED_STOCHASTIC` |
| `p_kill` exactly `0.25` | `SURVIVED_STOCHASTIC` |
| `p_kill` = `0.25 - 1 ulp` | `UNATTRIBUTED` |

Three of the four boundaries are half-open upward (`<`), so the boundary value belongs to the
*higher* band. The resistance test is the exception: `m <= resistant_modulation` is inclusive
*downward*, so a cell **exactly on** the threshold is called RESISTANCE-limited. That single
asymmetry decides the label of 92.7 % of the flagship's tumour cells (§3.1).

### 1.2 Over-determination is reported, not hidden

A cell can be short of drug *and* below the resistance threshold. Exposure binds first; the
module's stated reason is counterfactual leverage (below IC50 no restoration of sensitivity
lifts the Hill term past half its ceiling). The count of such cells is returned as
`joint_limited` — for temozolomide on the flagship it is **19,921 of 21,268 exposure-limited
survivors, 93.67 %**. That number is in the payload and `js/ui/finding.js` renders its note; it
is not in the headline sentence, which says DOSE.

---

## 2. Provenance of every threshold, honestly

| threshold | value | where the number comes from | verdict |
|---|---|---|---|
| `adequate_ratio` | `1.0` | `E == IC50` is the definition of the half-effect point and is exactly where the Hill term equals half its ceiling, for any `h`. | **A definition.** Not a choice. Solid. |
| `unreached_ratio` | `0.1` | Nothing. It was born at `0.1` in the module's first commit (`ac5cb4c`, 2026-07-28) with no derivation, no citation and no measurement. The stated justification — "at `h = 1.2` that is 5.9 % of the drug's achievable kill probability" — is arithmetic *downstream* of the choice (`0.1^1.2/(1+0.1^1.2) = 0.0593`), not a reason for it. Any constant produces some percentage. | **An unexamined convention, correctly labelled as one** by the module itself ("the one genuinely conventional constant here"). Its consequences are measured in §3.3 and they are not small. |
| `resistant_modulation` | `0.5`, derived | `max(pd.hypoxia_resistance, pd.quiescent_resistance)` = `max(0.35, 0.5)`. Described as "the mildest resistance penalty the kill model itself applies". | **Derived, and the derivation is circular on this section** (§3.1) — and the description is **false for temozolomide** (§2.1). |
| `stochastic_min_p_kill` | `0.25`, derived | `f(adequate_ratio) * resistant_modulation` = `0.5 * 0.5`. `f(1.0) = 0.5` for every `h`, so this is always `0.5 * resistant_modulation`. | **Derived and internally consistent.** But see §2.2: it makes `UNATTRIBUTED` unreachable. |

Both "derived" values are derived from `PDConfig` — and `PDConfig.hypoxia_resistance = 0.35`
and `quiescent_resistance = 0.5` are themselves, in `core/config.py`'s own words, **"UNVALIDATED
PRIORS with no citation"**. Deriving a threshold from an uncited prior is better than picking a
second uncited number, because the two then move together; it is not a measurement.

### 2.1 The per-drug hypoxia override is invisible to the derivation

`DrugParams.hypoxia_resistance` is a per-drug override and `pkpd.effective_hypoxia_resistance`
makes it win over the global default. Temozolomide sets it
(`core.config.TMZ_HYPOXIA_RESISTANCE = 0.80`, literature-derived), on both the synthetic library
entry and the real compile path. So on a temozolomide run the mildest resistance penalty the
kill model actually applies is **0.80**, while `resolve()` derives the band at **0.5** from a
`PDConfig` that run did not use. `engine.simulate_response` passes `pd=cfg.pd` and not the drug,
so the override cannot reach the derivation even in principle — and `diagnostics` carries the
applied value (`hypoxia_resistance: 0.8`, confirmed live in the run) which the module never
reads.

Consequence, and it moves a label: a hypoxic temozolomide cell with sensitivity 1.0 has
`modulation = 0.80`, above the band, so an adequately-exposed one is reported to the user as
**CHANCE** — "its neighbour with identical inputs died" — when the kill model demonstrably
suppressed its kill probability by 20 % for being hypoxic. Asserted by
`test_per_drug_hypoxia_override_is_invisible_to_the_derived_threshold`, plant-tested in §7.

**This defect does not bite on the flagship section**, because that section has **zero** hypoxic
cells (§3.1) — hypoxia is UNMEASURED on the 319-gene panel (`results/cell_state_impact.md`) and
the oxygen route banded none. It is latent, not active, and it is reported as latent.

### 2.2 `UNATTRIBUTED` is the canary, and it is wired to a condition the product cannot reach

The module's argument for `UNATTRIBUTED` is that it is "empty by construction" at
`max_kill_rate == 1`, which is what makes a non-zero count informative. The arithmetic is
correct. But **no product path sets `max_kill_rate` to anything but 1.0** (the only other write
in the repo is a read-back in `app/server/uncertainty_api.py`), and the solver does not produce
non-finite exposures. So its count is structurally zero in every shipped run — it was 0 in all
three runs here and in all 26 census arms. A canary that cannot sing is not evidence that the
air is clean; it is decoration until someone exposes `max_kill_rate`.

---

## 3. Distance to boundary, on real runs

Flagship arm, taken unchanged from `results/permeability_split.md` §4.2: 10x Xenium human GBM
section, `dx = 12 um`, `repair = 1`, `dt = 10 s`, 6 h window, convection ON, `save_every = 500`,
seed 0, `PKConfig(plasma_peak=50 uM, half_life=14400 s, n_doses=1)`, **`transvascular_exchange =
False` (shipped)**, compounds compiled by the product's own path
(`calibrated_drugparams` → `apply_barrier_to_drugparams`), nothing scaled. 40,887 cells, 23,812
tumour.

**Metric.** For each surviving tumour cell, the `|log10|` of the smallest multiplicative change
in the quantity that governs its label which would move it into a different category — i.e. the
fold change needed to flip it. The Bernoulli outcome (`alive`, `killed`) is held exactly as the
run produced it, so this isolates *classifier* fragility from *simulation* fragility. `0.30103`
is a factor of 2; `0.0` is an exact tie.

The denominator is survivors with a **finite** distance. The 140 cells whose exposure is exactly
`0.0` (§3.4) have infinite multiplicative distance — no scaling of zero is anything but zero —
so they are excluded rather than counted as robust: 21,128 of temozolomide's 21,268 survivors,
12,532 of doxorubicin's 12,672, 11,366 of cisplatin's 11,506.

### 3.1 The resistance axis: every cell is on the line

The flagship section's per-cell `sensitivity` is uniform `1.0` (real Xenium reconstruction sets
a neutral sensitivity; `reconstruct/real_pipeline.py`). So `modulation = resistance`, and
`resistance` takes one of three discrete values. Measured over all 23,812 tumour cells:

| `modulation` | cells | % | what it is | `== resistant_modulation` (0.5)? |
|---|---|---|---|---|
| `0.5` | **22,074** | 92.70 % | the QUIESCENT multiplier | **yes, exactly** |
| `1.0` | **1,738** | 7.30 % | VIABLE, nothing suppressed | no — disqualified by the `m < 1.0` clause, also an exact tie |
| the hypoxia multiplier | **0** | 0.00 % | no hypoxic cells on this section | — |

There is **not one cell in the interior of the resistance classifier.** 92.7 % sit exactly on
`resistant_modulation` and are called RESISTANT by the inclusive `<=`; the other 7.3 % sit
exactly on `1.0` and are excluded by the strict `<`. The threshold is
`max(0.35, quiescent_resistance)` and `quiescent_resistance` is `0.5` — so the band is derived
from the very constant that generates the values being tested. It is an identity, not a
comparison.

What that costs, measured by scaling `modulation` by `1 + 1e-12` and re-classifying:

| compound | survivors | re-labelled by a 1e-12 nudge | dominant limiter before → after |
|---|---|---|---|
| doxorubicin | 12,672 | **12,350 (97.46 %)** | **resistance (0.9746) → chance (0.9861)** |
| cisplatin | 11,506 | **11,325 (98.43 %)** | **resistance (0.9843) → chance (0.9878)** |
| temozolomide | 21,268 | 0 (0.00 %) | dose → dose (no survivor reaches the resistance test) |

A nudge of 1 part in 10¹² is below any physical or numerical meaning. **For doxorubicin and
cisplatin the headline limiter is decided by the direction of an inequality at exact equality,
and nothing else.** Scaling by `0.99` or `0.9` changes nothing at all — the flip is one-sided,
which is the signature of a tie rather than of a distribution near a threshold.

### 3.2 The full distance distribution (not a summary statistic)

Fraction of surviving tumour cells whose nearest label-changing boundary is the stated fold
change away:

| distance to nearest boundary | temozolomide | doxorubicin | cisplatin |
|---|---|---|---|
| **exactly 0 (tie)** | 0 (0.00 %) | **12,350 (98.55 %)** | **11,325 (99.64 %)** |
| < 1.0233x | 4 (0.02 %) | 1 (0.01 %) | 0 |
| 1.0233–1.05x | 4 (0.02 %) | 0 | 0 |
| 1.05–1.10x | 7 (0.03 %) | 2 (0.02 %) | 0 |
| 1.10–1.25x | 41 (0.19 %) | 8 (0.06 %) | 0 |
| 1.25–1.5x | 279 (1.32 %) | 10 (0.08 %) | 0 |
| 1.5–2x | 1,991 (9.42 %) | 6 (0.05 %) | 0 |
| 2–5x | **18,802 (88.99 %)** | 155 (1.24 %) | 41 (0.36 %) |
| 5–10x | 0 | 0 | 0 |
| > 10x | 0 | 0 | 0 |

Percentiles of the same distribution, in `|log10|` fold:

| | p1 | p5 | p10 | p25 | **p50** | p75 | p90 | p99 |
|---|---|---|---|---|---|---|---|---|
| temozolomide | 0.1493 | 0.2533 | 0.2955 | 0.3532 | **0.4197 (2.63x)** | 0.4611 | 0.4848 | 0.4986 |
| doxorubicin | 0.0 | 0.0 | 0.0 | 0.0 | **0.0 (exact tie)** | 0.0 | 0.0 | 0.3010 |
| cisplatin | 0.0 | 0.0 | 0.0 | 0.0 | **0.0 (exact tie)** | 0.0 | 0.0 | 0.0 |

**Read plainly: for two of the three compounds the modal cell — and the 99th-percentile cell for
cisplatin — is at distance exactly zero.** For temozolomide the modal cell is a comfortable
2.6-fold from any boundary and the classification is genuinely robust *on the exposure axis*;
that robustness is real and it is the one place this module earns the word.

The 146 doxorubicin and 41 cisplatin CHANCE cells all sit at exactly `log10 2` — they are the
`modulation == 1.0` cells, a factor of two from the resistance band, which is the same tie seen
from the other side.

### 3.3 The exposure axis, and what the conventional 0.1 buys

Scaling exposure `E` by `k` and re-classifying (Bernoulli outcome held fixed):

| k | temozolomide re-labelled | doxorubicin | cisplatin |
|---|---|---|---|
| 0.99 | 1 (0.00 %) | 1 (0.01 %) | 0 |
| 1.01 | 2 (0.01 %) | 0 | 0 |
| 1.10 | 8 (0.04 %) | 3 (0.02 %) | 0 |
| 1.25 | 35 (0.16 %) | 9 (0.07 %) | 0 |
| 2.00 | 880 (4.14 %) | 22 (0.17 %) | 0 |
| 0.50 | 1,446 (6.80 %) | 347 (2.74 %) | 1 (0.01 %) |

The exposure bands are **not** on a knife edge at ±1–25 %. That is a genuine robustness result.

The convention itself is a different matter. Sweeping `unreached_ratio` and re-classifying:

| `unreached_ratio` | temozolomide delivery : dose (of exposure-limited) | dominant limiter |
|---|---|---|
| 0.01 – 0.05 | 140 : 21,128 (0.66 %) | dose |
| **0.1 (shipped)** | **150 : 21,118 (0.71 %)** | **dose** |
| 0.2 | 1,596 : 19,672 (7.50 %) | dose |
| **0.291808 — the measured crossing** | **delivery reaches half of all survivors** | **inverts here** |
| 0.3 | 11,617 : 9,651 (54.6 %) | **delivery** |
| 0.5 | 20,398 : 870 (95.9 %) | delivery |
| 0.9 | 21,268 : 0 (100 %) | delivery |

**Temozolomide's headline inverts from DOSE to DELIVERY at `unreached_ratio = 0.2918082`** —
which is, necessarily, just the median survivor `E/IC50` (`p50 = 0.2918052`), because the
DELIVERY band is "below `u`". The census
(`results/verdict_inversion_census.md`) reports the same DOSE → DELIVERY inversion for
temozolomide caused by a missing physics term. **There is a second, independent route to exactly
the same reversal that costs nothing but moving an admittedly conventional constant from 0.1 to
0.29, and it is not disclosed anywhere.** The panel prints the threshold value
(`js/ui/finding.js`, `js/ui/inspect.js`) but says nothing about how close the run is to it.

For doxorubicin and cisplatin the sweep does almost nothing (doxorubicin 140 → 172 delivery
cells across the whole range 0.01 → 0.9; cisplatin a flat **140 at every single value**). The
reason is §3.4.

### 3.4 The DELIVERY band is, for most compounds, a fixed 140-cell geometric fact

On the flagship section **140 tumour cells receive exposure exactly `0.0`**. Measured directly
on the temozolomide run: `n(E/IC50 == 0) = 140`, `n(< 0.001) = 140`, `n(< 0.01) = 140` — one set
of cells, identical at every cut. The `unreached_ratio` sweep shows the same floor for the other
two: cisplatin's DELIVERY band is **exactly 140 cells at every value from 0.01 to 0.9**, and
doxorubicin's is 140 at 0.01 rising only to 172 at 0.9. The census's OFF arm shows an
`SURVIVED_UNREACHED` count of 140–150 for **ten of the thirteen** shipped compounds (the
exceptions are procarbazine at 174, methotrexate at 286, and the degenerate bevacizumab at
23,812, whose source term is identically zero everywhere).

For any compound whose remaining exposure sits above `unreached_ratio`, "DELIVERY-limited" is
therefore **not a threshold call at all** — it is the count of cells the source term never
reaches, which is a property of the tissue mask and the vessel geometry and is the same number
whatever drug is run and whatever the convention is set to. It is reported as a *finding* about
the compound. For cisplatin, "1.2 % of survivors were DELIVERY-limited" is arithmetically
`140 / 11,506` and would be that for any drug.

---

## 4. The reconciliation: verified independently, and what it does not constrain

### 4.1 It is exact, and independently so

Recomputed outside the module, from the engine's own post-kill population, and compared as IEEE
doubles in hex:

| compound | `n_tumor` from cells | from labels | engine `kill_fraction` (hex) | rebuilt from labels (hex) | `==` |
|---|---|---|---|---|---|
| temozolomide | 23,812 | 23,812 | `0x1.b59a9901ad590p-4` | `0x1.b59a9901ad590p-4` | **True** |
| doxorubicin | 23,812 | 23,812 | `0x1.df0f2e248d914p-2` | `0x1.df0f2e248d914p-2` | **True** |
| cisplatin | 23,812 | 23,812 | `0x1.0899c7d65f61cp-1` | `0x1.0899c7d65f61cp-1` | **True** |

The claim in the module docstring is true as stated, and `attribute()` really does raise rather
than return on a mismatch (`test_a_lying_attribution_raises_instead_of_returning`).

### 4.2 And it says nothing whatsoever about the headline

`reconcile()` computes `1.0 - (n_tumor - (KILLED + ALREADY_DEAD)) / n_tumor` over
`labels != NOT_TUMOR`. **The four survival categories never enter that arithmetic.** So an
attribution whose entire delivery/dose/resistance/chance split is random noise reconciles
exactly as bit-exactly as the real one. Demonstrated on all three real runs by replacing every
survivor's label with a uniform draw over the four codes:

| compound | scrambled split (delivery : dose : resistance : chance) | still bit-exact? |
|---|---|---|
| temozolomide | 5,267 : 5,280 : 5,476 : 5,245 | **yes** |
| doxorubicin | 3,087 : 3,146 : 3,309 : 3,130 | **yes** |
| cisplatin | 2,802 : 2,843 : 3,009 : 2,852 | **yes** |

Putting *all* survivors in one bucket reconciles too, for any of the four.

This matters because `js/ui/finding.js` — the first screen — tells a reader: *"It invents no
number. Every figure is read off `run.mechanism`, which reconciles bit-exactly with the engine's
own labels."* Both halves of that sentence are true and the second is **not** evidence for the
first as applied to the split. Bit-exactness constrains the KILL TOTAL; it constrains the WHY by
exactly zero. The property is now asserted by
`test_reconciliation_is_blind_to_the_survival_split` so the two cannot be quietly conflated in
code (§7).

---

## 5. Is the panel's caveat strong enough? A number.

`js/ui/finding.js` carries a compound-aware critical caveat: *"THIS VERDICT REVERSES FOR
&lt;compound&gt;"* / *"it reverses for 7 of the 12 shipped compounds"*, with the 58.3 % census
figure, the 57 % churn, the mechanical stable/unstable rule and both `results/` sources. **As a
disclosure of the transvascular-exchange inversion it is strong, specific, quantified and
correctly placed.** Nothing below asks for it to be softened.

It is not strong enough in three measurable respects, and all three are about fragilities it
does not mention.

**(a) It discloses one route to reversal. There are at least three, and the other two are
cheaper.**

| route to the same inversion | disclosed? | cost |
|---|---|---|
| switch on the missing transvascular exchange term | **yes**, prominently | a physics term the repo does not have the parameters for |
| move `unreached_ratio` from 0.1 to 0.2918 | **no** | changing one admittedly conventional constant |
| nudge `modulation` by 1 part in 10¹² | **no** | nothing; it is a tie |

**(b) The 58.3 % figure has a compound 0.22 % away from making it 66.7 %.** Bisecting a uniform
scaling of exposure until each compound's dominant limiter flips:

| compound | uniform exposure scale at which the headline flips | the exchange arm's measured mean AUC ratio | margin | census verdict |
|---|---|---|---|---|
| temozolomide | 0.3427 | 0.1401 | −59.1 % (well past) | inverts |
| doxorubicin | 0.1787 | 0.8373 | +368.6 % (well clear) | holds |
| **cisplatin** | **0.05431** | **0.05443** | **+0.22 %** | **holds** |

Cisplatin holds its dominant limiter by **two parts in a thousand** of exposure. Had the
exchange arm cost it 0.2 % more, the census headline would read 8 of 12 (66.7 %), not 7 of 12
(58.3 %). *The two perturbations are not identical — the exchange law reshapes the field rather
than scaling it uniformly — so this is an order-of-magnitude statement about the margin, not a
recomputation of the census.* It is still the right order: the census's own number is one
compound away from moving, and the panel presents 58.3 % as a measured rate without a margin.

**(c) The panel's own mechanical rule is contradicted by that compound.** The caveat offers a
reader a way to know which verdicts to trust: *"every compound that inverts loses most of its
exposure under the split (0.06 to 0.63 of the original), and every compound that holds keeps
most of it (0.69 to 0.93)."* Checked against `results/verdict_inversion_census_data.json`: the
inverting set spans 0.0623–0.6331 ✓ and the named stable four span 0.6931–0.9273 ✓ — but
**cisplatin holds its dominant limiter while retaining 0.0544, the lowest retention on the whole
panel.** It is counted among the five non-inverters that make the rate 7 of 12, and it is
correctly left out of `VERDICT_STABLE_COMPOUNDS` (it *does* invert on the narrower delivery-vs-
dose axis, so `verdict_stability("cisplatin")` returns `stable: False`) — so the product acts on
8 unstable compounds while its sentence says 7 reverse. A cisplatin user is shown both
statements at once. The rule as worded is false for it; the flag as implemented is right.

*(b) and (c) are observations about `app/server/jobs.py` and `app/web/js/ui/finding.js`, which
are not this audit's territory. Nothing there was changed. They are reported for whoever owns
them.*

**What the caveat should carry, stated as a claim someone can check:** on the flagship section,
where the verdict is RESISTANCE it is carried by cells sitting exactly on the threshold and
inverts under a 1e-12 perturbation, and where the verdict is DOSE-vs-DELIVERY for temozolomide
it inverts at `unreached_ratio = 0.29` against a shipped convention of 0.10. Neither is in the
panel.

---

## 6. What this audit does **not** say

* It does not say the attribution is wrong. Every label follows from the kill model's own
  arithmetic, the re-derivation was bit-exact, and the reconciliation holds.
* It does not say the exposure axis is fragile. Measured, it is not: scaling exposure by 0.8 or
  by 1.25 re-labels 0.10 % and 0.16 % of temozolomide's survivors respectively.
* It does not propose a fix. Moving `<=` to `<`, or deriving `resistant_modulation` from the
  effective per-drug resistance, would change the shipped verdict for two of three compounds
  measured here — that is a decision with a claim surface attached, not a tidy-up, and it is not
  made here.
* It says nothing about any patient, any response, or any compound's merit.
* One section, one arm, three compounds. The tie in §3.1 is a property of *this* reconstruction
  (uniform `sensitivity = 1.0`, discrete state multipliers); a section with per-cell sensitivity
  would spread the modulation distribution and the tie would dissolve. That is exactly why it is
  reported as a property of the flagship and not of the classifier in general — and the flagship
  is the section the Studio opens on.

---

## 7. Guards added, and the plant test for each

Three guards added to `tests/test_mechanism.py`. Each pins a fact measured above; each was
planted against, confirmed REACHED, watched go red, and reverted (`tmesim/sim/mechanism.py`
sha256 `013068ef…` before and after, verbatim below).

| guard | pins | plant that must make it red |
|---|---|---|
| `test_reconciliation_is_blind_to_the_survival_split` | §4.2 — bit-exactness constrains the total only | make `reconcile()` sensitive to the CHANCE count |
| `test_per_drug_hypoxia_override_is_invisible_to_the_derived_threshold` | §2.1 — the known defect, and the label it moves | make `resolve()` honour the 0.80 override |
| `test_below_ic50_no_biology_can_reach_the_resistance_or_chance_bands` | §1 / §3.4 — the precedence's consequence | test biology before the exposure bands |

```
013068efc871326e7c83c70e354643ea2648ac272f80a0a9951451dc673da411  tmesim/sim/mechanism.py

### PLANT P1 - reconcile() made SENSITIVE to the survival split
planted P1
--- pytest -k reconciliation_is_blind ---
E       assert False is True

tests/test_mechanism.py:586: AssertionError
=========================== short test summary info ============================
FAILED tests/test_mechanism.py::test_reconciliation_is_blind_to_the_survival_split
1 failed, 25 deselected in 0.81s
--- reverted, sha: 013068efc871326e7c83c70e354643ea2648ac272f80a0a9951451dc673da411 ---
.                                                                        [100%]
1 passed, 25 deselected in 0.51s

### PLANT P2 - resolve() silently honours a per-drug 0.80 override
planted P2
--- pytest -k per_drug_hypoxia ---
E        +  where 0.8 = MechanismThresholds(unreached_ratio=0.1, adequate_ratio=1.0, resistant_modulation=0.8, stochastic_min_p_kill=0.4).resistant_modulation

tests/test_mechanism.py:630: AssertionError
=========================== short test summary info ============================
FAILED tests/test_mechanism.py::test_per_drug_hypoxia_override_is_invisible_to_the_derived_threshold
1 failed, 25 deselected in 0.69s
--- reverted, sha: 013068efc871326e7c83c70e354643ea2648ac272f80a0a9951451dc673da411 ---
.                                                                        [100%]
1 passed, 25 deselected in 0.55s

### PLANT P3 - precedence inverted, biology tested before exposure
planted P3
--- pytest -k below_ic50 ---
E         Use -v to get more diff

tests/test_mechanism.py:668: AssertionError
=========================== short test summary info ============================
FAILED tests/test_mechanism.py::test_below_ic50_no_biology_can_reach_the_resistance_or_chance_bands
1 failed, 25 deselected in 0.67s
--- reverted, sha: 013068efc871326e7c83c70e354643ea2648ac272f80a0a9951451dc673da411 ---
.                                                                        [100%]
1 passed, 25 deselected in 0.57s

### FINAL STATE
013068efc871326e7c83c70e354643ea2648ac272f80a0a9951451dc673da411  tmesim/sim/mechanism.py
.............                                                            [100%]
13 passed, 13 deselected in 5.71s
```

Line 586 is `assert noise["bit_exact"] is True`, line 630 is
`assert th.resistant_modulation == 0.5`, line 668 is the label-order assertion — each is the
assertion that carries its guard's claim, so each guard is REACHED and not merely present.

---

## 8. Reproduction, and the frozen-state assertions

Both measurement runs asserted `calibrate_alpha()` live at the top and the bottom, and
`PDEConfig.transvascular_exchange` at the shipped `False`:

```
[START] calibrate_alpha() = 22.223477711172563
...
[END]   calibrate_alpha() = 22.223477711172563
```

Self-checks performed before any number above was computed, all three compounds:

| check | temozolomide | doxorubicin | cisplatin |
|---|---|---|---|
| `ratio` re-derived from captured `E` is bit-identical to the engine's | True | True | True |
| `p_kill` re-derived is bit-identical to the engine's | True | True | True |
| labels re-derived through the shipped `classify()` are bit-identical to `res.mechanism["labels"]` | True | True | True |

Every counterfactual in §3 is therefore a re-classification of the *same* run through the *same*
`classify()`, with one input perturbed and the Bernoulli outcome held fixed — not a second model
and not a second solve.

`kill_fraction` for the record, as the reconciliation target only and not as a result:
temozolomide `0.10683688896354782`, doxorubicin `0.4678313455400638`, cisplatin
`0.5167982529816899`. The temozolomide value is bit-identical to the one stored in
`results/verdict_inversion_census_data.json`, which is the independent check that this audit ran
the same arm the census did.
