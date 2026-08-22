# Mechanistic result visualisation — per-cell attribution

> ## NUMBERS SUPERSEDED 2026-08-22 — the METHOD stands, the MEASUREMENTS do not
>
> **Last regenerated at `3f42f03` (2026-08-08). Two physics changes landed after
> it and neither triggered a regeneration:** `d77d169` (2026-08-13, area-conserving
> vessel source) and `fd2a9e7` (2026-08-14, two-route transport ON). Both move
> `e_unit` for every drug, so both move every exposure and therefore every band
> assignment in this file.
>
> **The disagreement this causes is not cosmetic — it flips the headline verdict.**
> §3.2 here reports temozolomide on the flagship section as **97.10 % DELIVERY-limited**
> with a median exposure shortfall of **17.9×**. The current engine reports the same
> compound on the same section as **99.8 % DOSE-limited**
> (`results/synthesised_positions.md` §5.2). A reader taking §3.2's verdict as
> current would reach the opposite conclusion about what limits this drug.
>
> **What is unaffected:** the attribution METHOD — the six labels, the band
> definitions, the reconciliation against `kill_fraction`, and the argument that a
> per-cell mechanism is recoverable at all. That is what this file is for and none
> of it depends on the numbers. The per-cell attribution is now the strongest
> surviving output in the product, and it is strong because of the method
> described here.
>
> **Where the current numbers live:** `results/synthesised_positions.md` §5.2 for
> the flagship verdict, and the Studio's own "Delivery or dose" panel, which
> recomputes live rather than quoting.
>
> Found on 2026-08-22 while sourcing a marketing claim. This is the tenth
> instance of a results file outliving the physics it describes, and the first
> that the 2026-08-21 banner sweep missed — that sweep looked for documents making
> *reachability* claims, and this one makes a *measurement* claim, which is a
> different shape and was not searched for.

**What this adds.** A run used to return scalars: `kill_fraction`, therapeutic coverage,
mean tumour AUC. Those say *how much* happened. They cannot say **why**, and the kill model
already knows: `tmesim/cells/pkpd.py` computes a local exposure, a Hill term, a resistance
multiplier, a kill probability and a Bernoulli draw for every single cell, and discards all
five one line later. `tmesim/sim/mechanism.py` picks them back up and labels every tumour cell
with the reason it died, or the reason it did not.

It is a **read-out of the model's own arithmetic, not a second model**. `apply_drug_kill` now
takes an optional `diagnostics=` dict and fills it with the exact arrays it decided on — the
same objects, not recomputed copies — and every label is a comparison between two of them.
The run is bit-identical with the attribution on or off, which is asserted rather than claimed
(`tests/test_mechanism.py::test_attribution_is_bit_identical_to_a_run_without_it` compares
`kill_fraction`, `viable_tumor_fraction` and the whole AUC field).

---

> ## STATUS — this report is measured under BOTH transvascular delivery laws
>
> **There are two delivery laws in this engine and the choice between them is
> an OPEN DECISION.** Every quantitative result below is therefore reported in
> two arms, and neither is presented as "the" answer:
>
> | arm | law | steady state | status as of this writing |
> |---|---|---|---|
> | **one-way** | `S = PS*Cp` | `PS*Cp / (clearance + binding)` — **unbounded above plasma** | **the SHIPPED default** (`PDEConfig.transvascular_exchange = False`, `tmesim/core/config.py`) |
> | **exchange** | `S = PS*(Cp - C)` | `PS*Cp / (PS + clearance + binding)` — **bounded by `Cp`** | opt-in (`PDEConfig(transvascular_exchange=True)`) |
>
> **The one-way law is known to be physically wrong** — it has no path back from
> the interstitium to the blood, so tissue concentration can exceed the plasma
> level, which passive delivery cannot do. `results/breast_kp_validation.md`
> found it against a measured human microdialysis anchor;
> `results/transvascular_exchange_revalidation.md` is the consolidated
> re-validation. On the runs in §3 below the one-way arm puts **cisplatin at
> 2.61x plasma exposure on average and 11.55x at maximum**, and even
> temozolomide at 1.10x mean / 3.07x max. The exchange arm puts them at 0.61x /
> 0.96x and 0.49x / 0.81x — inside the bound, as its algebra requires.
>
> **The prior version of this file carried a conditional-retraction notice.
> That notice is superseded by this rewrite:** the body no longer states the
> one-way result as if it were the only result. Both arms were re-measured for
> this rewrite, in one process, from one code path, with nothing but the flag
> changed (`experiments/mechanism_attribution_arms.py`).
>
> **What to quote.** If you quote a mechanism verdict from this file you must
> name the law it came from, because for cisplatin **the two laws give opposite
> verdicts** (§3.1). "Cisplatin is resistance-limited" is not a quotable
> statement about this tissue; "under the shipped one-way law cisplatin is
> resistance-limited, and under the corrected exchange law it is dose-limited"
> is.
>
> **What is NOT arm-dependent:** the category definitions (§1), the bit-exact
> reconciliation property (§2), and the two guards in §4. Those are properties
> of the attribution itself and hold under either law — verified in all four
> runs below.

---

## 1. Categories and thresholds

| code | label | limiter | condition |
|---|---|---|---|
| 0 | `KILLED` | — | the engine's own `killed` mask |
| 1 | `SURVIVED_UNREACHED` | delivery | `E/IC50 < unreached_ratio` |
| 2 | `SURVIVED_UNDERDOSED` | dose | `unreached_ratio <= E/IC50 < adequate_ratio` |
| 3 | `SURVIVED_RESISTANT` | resistance | exposure adequate, `m <= resistant_modulation` **and** `m < 1` |
| 4 | `SURVIVED_STOCHASTIC` | chance | both adequate, `p_kill >= stochastic_min_p_kill` |
| 5 | `UNATTRIBUTED` | — | none of the above. **Counted, never reassigned.** |
| 6 | `ALREADY_DEAD` | — | dead before the drug; counts in `kill_fraction`, not killed by it |
| 7 | `NOT_TUMOR` | — | non-tumour cell |

`m = sensitivity * resistance`, the exact multiplicative modulation the kill model applied.

**Two of the four thresholds are not choices, and two more are derived from the run's own
configuration rather than picked:**

| threshold | value | where it comes from |
|---|---|---|
| `adequate_ratio` | 1.0 | the *definition* of IC50 — `E == IC50` is exactly where the Hill term equals half its ceiling |
| `unreached_ratio` | 0.1 | the one conventional constant. At `hill = 1.2` it is 5.9% of the drug's achievable kill probability. Exposed as a parameter so it can be moved and the breakdown watched to move |
| `resistant_modulation` | **derived** 0.5 | `max(pd.hypoxia_resistance, pd.quiescent_resistance)` — the *mildest* resistance penalty the kill model itself applies |
| `stochastic_min_p_kill` | **derived** 0.25 | `f(adequate_ratio) * resistant_modulation` — the lowest `p_kill` the CHANCE branch can produce at `max_kill_rate == 1` |

These four resolved identically in all four runs of §3 — the thresholds are read off `PDConfig`
and the Hill exponent, neither of which the delivery law touches.

The last derivation is what makes `UNATTRIBUTED` a live check rather than decoration: under the
production `PDConfig` it is **empty by construction**, and it fills up precisely when something
outside the bands suppressed the kill (a `max_kill_rate` ceiling below 1, a clip, a non-finite
exposure). It was worth having — see §4.

**Precedence.** Exposure bands bind before biology, because below IC50 no restoration of
sensitivity lifts the Hill term past half its ceiling. Cells that are *both* short of drug and
resistant are not hidden: `joint_limited.count` reports how many exposure-limited survivals
would still not have died on dose alone. **This precedence rule is what makes the §3.1
inversion possible, and §3.5 says so explicitly.**

---

## 2. Bit-exact reconciliation — holds in both arms

`mechanism.reconcile` rebuilds `kill_fraction` from the labels with the **same integer
arithmetic** `CellPopulation.viable_tumor_fraction` uses —
`1.0 - (n_tumor - (KILLED + ALREADY_DEAD)) / n_tumor` — and compares with `==`, not `allclose`.
`simulate_response` **raises** `MechanismReconciliationError` rather than returning a breakdown
that does not add up, and no payload is written for such a run.

Measured on the real Xenium 319 section, `dx = 12 um`, 6 h window, `dt = 10 s`, seed 0,
40,887 measured cells of which 23,812 are TUMOR:

| drug | delivery law | `kill_fraction` (engine) | rebuilt from labels | `\|delta\|` | UNATTRIBUTED |
|---|---|---|---|---|---|
| temozolomide | one-way | `0.01839408701495049` | `0.01839408701495049` | `0.0` | 0 |
| temozolomide | exchange | `0.007223248782126657` | `0.007223248782126657` | `0.0` | 0 |
| cisplatin | one-way | `0.43658659499412056` | `0.43658659499412056` | `0.0` | 0 |
| cisplatin | exchange | `0.23118595666050734` | `0.23118595666050734` | `0.0` | 0 |

Identical to the last bit of the mantissa, not to a tolerance, in **every** arm.

**This is the point worth separating from everything else in this file.** The reconciliation is
a correctness property of the *attribution*, not of the *physics*: it asserts that the label
array describes the run that actually happened. It is therefore **independent of which delivery
law is in force**, and the table above is the evidence rather than the assertion. The delivery
law changes *what* the labels say; it cannot make them stop adding up. The same holds for
`UNATTRIBUTED = 0`: no cell in any arm fell outside the bands.

> **These are the time-resolved-PK numbers, and that matters when comparing files.** The runs
> in this report use `time_resolved = true` with `plasma_peak = 50 µM`, `half_life = 4 h`
> (stated in §3 below). Other files in `results/` report the **static-source** config for the
> same drug and section — cisplatin `0.4361666386695784`, temozolomide `0.01835209138249627`,
> both **one-way-law** numbers. Both pairs were re-run directly from
> `reconstruct_real('xenium', dx_um=12)` and each reproduces to the last bit, so the ~0.1% gap
> (≈10 of 23,812 tumour cells) is a **config difference, not a reproducibility failure**.
>
> Worth naming, because "reconciles bit-exactly" appears throughout this report and could be
> misread: it is a claim that the *label array* rebuilds *this run's* `kill_fraction`. It is
> **not** a claim that two runs of different configs — or of different delivery laws — agree.
> A kill quoted without its PK mode **and its delivery law** is ambiguous.

---

## 3. The finding: temozolomide vs cisplatin on the flagship section, under both laws

All four runs: `xenium_319`, `dx = 12 um`, 6 h, `dt = 10 s`, convection on,
`plasma_peak = 50 uM`, `half_life = 4 h`, `sensitivity_mode = uniform`, `time_resolved = true`,
seed 0. The **only** difference between the arms of a drug is
`PDEConfig.transvascular_exchange`.

| | TMZ · one-way | TMZ · exchange | cisplatin · one-way | cisplatin · exchange |
|---|---|---|---|---|
| compiled IC50 | 606.05 uM | 606.05 uM | 24.55 uM | 24.55 uM |
| compiled `bbb_penetrance` | 0.1254 | 0.1254 | 0.8577 | 0.8577 |
| mean tumour AUC (uM·s) | 7.710e5 | 3.355e5 | 1.826e6 | 4.157e5 |
| mean tissue AUC / plasma AUC | **1.10x** | 0.49x | **2.61x** | 0.61x |
| max tissue AUC / plasma AUC | **3.07x** | 0.81x | **11.55x** | 0.96x |
| therapeutic coverage 1× | 0.000 | 0.000 | **0.9955** | **0.0467** |
| `kill_fraction` | 0.0184 | 0.0072 | 0.4366 | 0.2312 |
| tumour cells killed / 23,812 | 438 | 172 | 10,396 | 5,505 |
| surviving tumour cells | 23,374 | 23,640 | 13,416 | 18,307 |
| **survival: DELIVERY** | **97.10 %** | **100.00 %** | 0.23 % | 0.17 % |
| **survival: DOSE** | 2.90 % | 0.00 % | 0.44 % | **95.65 %** |
| **survival: RESISTANCE** | 0.00 % | 0.00 % | **96.42 %** | 3.64 % |
| **survival: CHANCE** | 0.00 % | 0.00 % | 2.91 % | 0.54 % |
| UNATTRIBUTED | 0 | 0 | 0 | 0 |
| tumour E/IC50 — p05 / p50 / p95 / max | 0.037 / 0.056 / 0.090 / 0.157 | 0.020 / 0.025 / 0.032 / 0.042 | 1.79 / 3.19 / 5.84 / 14.37 | 0.56 / 0.79 / 1.00 / 1.21 |
| exposure shortfall at the median cell | 17.9× | 39.3× | 0.31× (past IC50) | 1.27× |
| penetration index | 0.5112 | 0.7140 | 0.4285 | 0.7210 |

`plasma AUC` is the integral of the driving plasma curve over the same 6 h window on the same
10 s step grid, `6.7165e5 uM·s`, identical in all four runs — the delivery law changes what
reaches the tissue, not what is in the blood. The tissue ratios are taken over tissue-mask
voxels of the AUC field; per-tumour-cell they are 1.15× / 0.50× (temozolomide) and 2.72× /
0.62× (cisplatin), with maxima 3.05× / 0.81× and 11.35× / 0.95×. **Every one-way figure in
those two rows is a physical impossibility for a passively delivered drug; no exchange figure
is.**

Raw counts behind the survival shares:

| label | TMZ · one-way | TMZ · exchange | cisplatin · one-way | cisplatin · exchange |
|---|---:|---:|---:|---:|
| `KILLED` | 438 | 172 | 10,396 | 5,505 |
| `SURVIVED_UNREACHED` | 22,696 | 23,640 | 31 | 31 |
| `SURVIVED_UNDERDOSED` | 678 | 0 | 59 | 17,511 |
| `SURVIVED_RESISTANT` | 0 | 0 | 12,936 | 666 |
| `SURVIVED_STOCHASTIC` | 0 | 0 | 390 | 99 |
| `UNATTRIBUTED` | 0 | 0 | 0 | 0 |

**Both arms are real runs, and that was checked rather than assumed.** Setting
`PDEConfig(transvascular_exchange=True)` is a **silent no-op** on any path that calls
`simulate_drug_transport` directly (`tmesim/validation/propagation.py`, and
`results/transvascular_exchange_revalidation.md` §6, where it fired for real). This module goes
through `sim.engine.simulate_response`, the one declared `PROPAGATED` path, and every run above
was executed inside `propagation.instrument()`: the ON arms were handed `exchange_k` on 1 of 1
solver calls, the OFF arms on 0 of 1, and the two arms' AUC fields were required to differ
bit-wise before any comparison was reported. The OFF arms additionally reproduce this file's
previously published `kill_fraction` values to the last bit (§2), so the baseline is the same
baseline.

### 3.1 The striking result: cisplatin's verdict INVERTS between the two laws

| cisplatin | one-way (shipped) | exchange |
|---|---:|---:|
| survival: RESISTANCE | **96.42 %** | 3.64 % |
| survival: DOSE | 0.44 % | **95.65 %** |
| IC50 coverage at 1× | 99.55 % | 4.67 % |
| median tumour E/IC50 | 3.19 | 0.79 |

**Under the shipped one-way law cisplatin is a biology problem. Under the corrected exchange
law it is a dose problem. These are opposite readings of the same tissue, the same drug, the
same dose and the same RNG stream.** 12,270 survivors move out of the RESISTANCE bucket
(12,936 → 666) and the UNDERDOSED bucket goes 59 → 17,511.

**The mechanism of the inversion is the documented precedence rule, and it is worth being
precise about because it is not a threshold artefact.** Exposure bands bind before biology
(§1). Under the one-way law essentially every tumour cell is driven past IC50 —
`coverage_1x = 0.9955`, `p05 E/IC50 = 1.79`, so the *fifth percentile* cell is already at 1.8×
IC50 — which leaves the UNREACHED and UNDERDOSED bands holding 90 cells between them out of
13,416 survivors, and by precedence essentially every survivor falls through to RESISTANCE.
Under the exchange law the same population sits at
`p95 E/IC50 = 1.00`: essentially the whole distribution is below IC50, so the exposure bands
capture it and the survivors are dose-limited. The label followed the exposure distribution,
which followed the delivery law. **Neither arm's label is "wrong given its physics"; the
physics is what is in question.**

Which arm to believe therefore reduces to which delivery law is right, and on that the evidence
is one-sided: the one-way arm puts cisplatin at **2.61× plasma exposure on average and 11.55×
at maximum**, and its peak tissue concentration at **555.0 uM against a 50 uM plasma peak**
(§5). Passive transvascular delivery cannot concentrate a drug above its plasma level. The
exchange arm is bounded by construction (0.61× mean, 0.96× max, 45.9 uM peak) and does not
breach it. **That is an argument about the physics, not about the attribution** — the
attribution reports faithfully in both arms — but it is the argument that decides which column
of this table describes a real tumour.

### 3.2 Temozolomide is DELIVERY-limited under both laws, and more so under exchange

The brief expected temozolomide to be "overwhelmingly dose- or resistance-limited". **It is
not, under either delivery law.** By the documented bands it is 97.10 % DELIVERY-limited under
the one-way law and **100.00 %** under exchange — the DOSE bucket empties entirely (678 → 0
cells), because the whole population drops further below the UNREACHED threshold.

The reason is not a marginal threshold call in either arm. **No tumour cell anywhere on the
section reaches even 16 % of IC50 under the one-way law (`max E/IC50 = 0.157`), or even 5 %
under exchange (`max E/IC50 = 0.042`)**, and the median cell is short by a factor of 17.9× and
39.3× respectively. There is no `unreached_ratio` in `[0.157, 1.0]` that would move a single
cell out of the exposure-limited regime and into RESISTANCE or CHANCE under either law; only
the split *between* UNREACHED and UNDERDOSED moves with that constant, and under exchange even
that split is gone because every cell is below 0.1.

`therapeutic_coverage_1x = 0.0` is the same statement in the pre-existing scalar, in both arms.
The attribution is what makes it actionable: it says the 0.0 is not a rounding of "nearly
enough", it is 23,374 (one-way) / 23,640 (exchange) individual cells that never saw a
pharmacologically meaningful concentration.

**This is the half of the report that the physics correction strengthens rather than
overturns.** The conclusion was law-independent to begin with — a drug two orders of magnitude
short of its own IC50 does not become resistance-limited when you cut its exposure a further
2.3× (mean tumour AUC 7.710e5 → 3.355e5 uM·s) — and the
corrected law moves it from 97.10 % to a categorical 100.00 %.

The cross-drug contrast survives in both arms but **changes shape**: under the one-way law it
is *delivery-limited vs resistance-limited*, and under exchange it is *delivery-limited vs
dose-limited*. What is common to both is that **temozolomide's problem is upstream of the
tumour's biology and cisplatin's is not the same problem as temozolomide's.**

### 3.3 What "DELIVERY-limited" does and does not mean here

This must be read precisely, because the categories describe *the kill model*, not biology in
general. Nothing in this subsection is arm-dependent — it is a property of the kill model's
functional form:

* `UNREACHED` and `UNDERDOSED` are **two bands on one axis**, `E/IC50`. The kill model reaches
  `p_hill` through that ratio alone, so it **has no separate potency term** and cannot
  distinguish "too little drug arrived" from "the drug is too weak". Neither can this
  attribution, and it says so in every payload (`exposure_axis_note`).
* What *is* separable within a run: IC50 is a single per-drug constant, so the entire
  cell-to-cell **spread** in `E/IC50` is delivery, while the population's **position** relative
  to 1.0 is set jointly by delivery and potency. Temozolomide's shortfall decomposes upstream
  into a 6.8× `bbb_penetrance` disadvantage (0.125 vs 0.858) and a 24.7× IC50 disadvantage
  (606 vs 24.6 uM) against cisplatin — but that decomposition comes from the drug compilation,
  not from the per-cell kill, and this module does not assert it.

So the honest headline is: **temozolomide fails on this tissue because exposure never comes
within an order of magnitude of its potency**, and the delivery/dose split is a statement about
*how far short* it fell, not about which upstream term caused the shortfall.

### 3.4 The resistance state is nearly universal, and it is the section's known weak point

`joint_limited` — exposure-limited survivors that are *also* at or below the resistance
threshold, i.e. cells a dose increase alone would not have killed:

| | one-way | exchange |
|---|---|---|
| temozolomide | 21,073 / 23,374 (**90.2 %**) | 21,283 / 23,640 (**90.0 %**) |
| cisplatin | 89 / 90 (98.9 %, of a tiny set) | 16,303 / 17,542 (**92.9 %**) |

That is the engine's internal oxygen re-banding calling almost every tumour cell HYPOXIC or
QUIESCENT, and the number barely moves with the delivery law — as expected, since the oxygen
solve does not depend on it. **Note what this does to the §3.1 inversion:** under the exchange
law cisplatin is 95.65 % DOSE-limited, but 92.9 % of those exposure-limited survivors are
*also* below the resistance threshold. "Dose-limited" is the label the precedence rule assigns;
it is **not** a claim that raising the dose would kill them.

It also lands directly on this section's documented limitation: the 319-gene panel **measures
no hypoxia**, so no state route clears the bar (`results/cell_state_impact.md` §5) and
`kill_fraction` on this section carries an unresolved uncertainty of up to 0.3252 absolute
purely from the choice of state route. The one-way arm's "96.4 % RESISTANCE-limited" inherits
that uncertainty in full — and so does the exchange arm's 3.64 %, in the other direction.
**The cross-drug contrast is the defensible output under either law; the absolute size of the
resistance bucket is not, under either law.**

### 3.5 A structural limit of this decomposition, made visible by the two arms

The precedence rule means the four categories are **not four independent measurements**. Where
a drug's whole exposure distribution sits relative to IC50 decides how many categories can be
occupied at all:

* **Population entirely below IC50** → RESISTANCE and CHANCE are **exactly zero cells**, not
  small. Temozolomide in both arms is this case (0 and 0).
* **Population almost entirely above IC50** → DELIVERY and DOSE are near-empty and nearly
  every survivor falls through to RESISTANCE or CHANCE. Cisplatin under the one-way law is this
  case (31 + 59 exposure-limited cells out of 13,416 survivors).
* **Population straddling IC50** → all four can be occupied. Cisplatin under exchange is this
  case: `max E/IC50 = 1.21` means a minority of cells *do* clear IC50, which is why RESISTANCE
  (666) and CHANCE (99) are small but **non-zero** rather than zero.

Under bounded physics the first two cases collapse toward the first: this section's own
exchange-arm coverages are 0.000 (temozolomide) and 0.0467 (cisplatin), and
`results/transvascular_exchange_revalidation.md` §4.3 reports temozolomide's coverage failing
even at a 30× dose (97.80 % → 3.43 %). So "delivery vs dose vs resistance" degenerates toward a
two-category answer more often than a four-row table suggests. Read the breakdown together with the `E/IC50` percentile row above it;
the percentiles are what tell you how many of the four categories were even reachable.

---

## 4. What the attribution refuses to do

`UNATTRIBUTED` was 0 on all four real runs — both drugs, both delivery laws — as the derived
thresholds guarantee under the production `PDConfig`. Two guards exist so that it is a real
check and not a formality. Both are properties of the labelling and are unaffected by the
delivery law:

* **A ceiling the bands do not model is not called "bad luck".** Drop `max_kill_rate` to 0.4
  and a cell at exactly the adequate boundary with no resistance has `p_kill = 0.2`, below the
  derived CHANCE floor of 0.25. It is labelled `UNATTRIBUTED` and counted, not filed under
  CHANCE. (`test_unexplained_survival_is_unattributed_not_reassigned`)
* **A defect found by this work.** `resistant_modulation` derives from the run's `PDConfig`, so
  the **counterfactual resistance-off arm** (`resistance_off: true`, used by the rank
  attribution) derives a threshold of 1.0 — under which `m <= threshold` is true for *every*
  cell, and 100 % of survivors in the one run where resistance was switched OFF would have been
  reported as RESISTANCE-limited with an empty CHANCE bucket: the exact inversion of the truth.
  `is_resistance_limited` therefore requires `m < 1.0` as well — "something was actually
  suppressed" — and the degenerate case is pinned by
  `test_resistance_off_arm_does_not_file_every_survivor_as_resistant`. The complementary test
  checks the guard does not evict QUIESCENT cells, which sit at `m == 0.5` exactly on the
  threshold.

---

## 5. Time-lapse

`PDEConfig.time_resolved` (`results/dose_schedule.md`) lets tissue concentration genuinely rise
**and fall**, so a time-lapse can show the drug clear rather than plateau.
`simulate_drug_transport(..., capture_frames=N)` keeps N snapshots of the **concentration**
field `C` — not the AUC field, which only accumulates and so can never show drug leaving —
evenly spaced over the window and always including `t = 0` and `t = total_time`. Capture is a
`copy` at selected steps: no arithmetic changes, and
`test_frame_capture_does_not_perturb_the_solve` asserts the AUC and final-concentration fields
are bit-identical with capture on and off.

Measured on the four runs above (24 frames, 939.1 s nominal spacing snapped to the 10 s step
grid → 940.0 s, 6 h window, 424×675 grid, 27.5 MB payload each), peak tissue concentration by
frame, against a **50 uM plasma peak**:

```
temozolomide  one-way    0 → 144.3 → 142.4 → 136.4 → … →  57.8 →  55.2 uM
temozolomide  exchange   0 →  39.0 →  37.3 →  35.7 → … →  15.1 →  14.4 uM
cisplatin     one-way    0 → 555.0 → 530.5 → 507.0 → … → 214.8 → 205.3 uM
cisplatin     exchange   0 →  45.9 →  43.9 →  42.0 → … →  17.8 →  17.0 uM
```

Drug enters, spreads, and then clears in every arm — a 2.6× fall for one-way temozolomide and
2.7× for the other three across the window. **The clearance shape is not what the delivery law
changes; the level is.** Under the one-way law both drugs peak *above* the 50 uM plasma peak
(2.9× for temozolomide, 11.1× for cisplatin); under exchange neither does (0.78× and 0.92×),
which is the bound `C <= Cp` showing up directly in a filmable quantity.

On the **default static** plasma source the same series rises monotonically toward a
quasi-steady state instead (temozolomide peaks at 96.1 uM at `t = T` under the one-way law —
carried over from this file's earlier static-source measurement; the static path was **not**
re-run in the exchange arm for this rewrite, so no exchange counterpart to that number is
claimed), because the delivery rate is held at its window average by construction. The run says so:
`RunResult.timelapse.plateau_note` carries that sentence whenever `time_resolved` is false, and
`/api/validate` warns before the run is submitted.

---

## 6. Reproduce

Both arms of §2, §3 and §5, in one process, with the propagation guard and the
non-bit-identity assertion built in:

```bash
./.venv/bin/python experiments/mechanism_attribution_arms.py --out arms.json
```

The single-arm (shipped-default) path through the API, and the unit tests:

```bash
./.venv/bin/python -m app.server          # then POST the spec below, or:
./.venv/bin/python -m pytest tests/test_mechanism.py -q
./.venv/bin/python -m pytest tests/test_propagation.py -q   # the flag-propagation guard
```

```jsonc
{ "section": "xenium_319", "dx_um": 12.0, "drug": "temozolomide",
  "sensitivity_mode": "uniform", "vascular_mode": "discrete",
  "oxygen_mode": "pde", "state_mode": "oxygen", "mgmt": null,
  "pk": {"plasma_peak": 50.0, "half_life_s": 14400.0,
         "n_doses": 1, "dose_interval_s": 86400.0},
  "total_time_s": 21600.0,
  "attribute_mechanism": true, "timelapse_frames": 24, "time_resolved": true }
```

The RunSpec has **no field for the delivery law** — it follows the engine default, which is
`transvascular_exchange = False` as of this writing. A run submitted through the API therefore
produces the *one-way* column of §3. If that default is flipped, the API will produce the
*exchange* column and this sentence, not the tables, is what needs editing.

`GET /api/runs/{id}/mechanism` (TMSM, 45.8 kB for 40,887 cells) and
`GET /api/runs/{id}/timelapse` (TMSL, 27.5 MB for 24 frames) serve the two payloads; the wire
formats are frozen in `API_CONTRACT.md`.

### Provenance and one stated discrepancy

Every number in §2, §3 and §5 comes from the four runs executed by
`experiments/mechanism_attribution_arms.py` for this rewrite, on `main`. The one-way arms
reproduce this file's previously published `kill_fraction`, survival breakdown and time-lapse
series exactly. AUC-field blake2b-128 digests, the evidence that the arms are four distinct
runs rather than two: temozolomide `399b0da9…` (one-way) / `36c7b4f6…` (exchange), cisplatin
`3da8b796…` (one-way) / `033d2f2a…` (exchange).

**Stated gap, not filled.** An earlier, independent exchange-arm re-run of this module — done
on branch `transvascular-exchange`, recorded in `results/transvascular_exchange_revalidation.md`
§4.2 and in the conditional-retraction notice this rewrite replaces — reports cisplatin
**RESISTANCE 3.59 % / DOSE 95.71 %** and 657 surviving RESISTANT cells, against the
**3.64 % / 95.65 %** and 666 cells measured here; and mean/max tissue-to-plasma for the one-way
arm as 2.60× / 11.53× against 2.61× / 11.55× here. The `kill_fraction` (0.2311), the IC50
coverage (4.67 %) and the median `E/IC50` (0.79) agree between the two. **The residual
~0.05–0.06 pp disagreement in the survival shares is not explained here** — the two re-runs
were made from different working trees, and no single rescale relates the two tissue/plasma
pairs. (The superseded notice's own arithmetic is also internally inconsistent: it reports
12,943 survivors leaving the resistance bucket and 657 remaining, which sums to 13,600 rather
than the 13,416 survivors its own 96.42 % share implies. The counts measured here,
12,936 → 666 of 13,416, do sum.) None of this touches the inversion, which is a
92-percentage-point move, but it is recorded rather than reconciled away. The numbers in this
file are the ones this file measured.

**Not for clinical use.** This tissue is an untreated resection with no outcome label. Every
number here is a prediction of the physical model and has never been scored against a measured
treatment response — under either delivery law.
