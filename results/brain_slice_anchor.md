# TEST B executed — rat brain-slice tissue partitioning (`Vu,brain`)

**Status: RUN. Pre-registered. FAIL.**

This is the execution record for **TEST B** of `results/transport_anchor_prereg.md`.
That document was written by another agent and is binding here; this one does not
amend it and may not (prereg §9.5 reserves amendment to the coordinator). Every
threshold, exclusion code and statistic below is transcribed from it. Nothing was
tuned, nothing was substituted for a more flattering alternative, and no compound
was dropped for scoring badly.

* script: `experiments/brain_slice_anchor.py` (`--stage prereg | floor | engine`)
* frozen execution record: `results/sweep/transport_anchor_prereg.json` (prereg §7.3)
* engine state at execution: commit `544734e469fa5a22f710248cd10a73db394e59d3`,
  every §2.3 constant asserted against live code and matching
* scope: `SCOPE_ANALYTIC`, conditional `coefficient-not-field`. The transport solver
  was **not entered** — call counters on
  `physics.diffusion.simulate_drug_transport` and `sim.engine.simulate_response`
  both read **0**.

---

## 0. Summary

### 0.1 The headline finding — `GMFE_floor` (prereg §3.5a)

> ## `GMFE_floor` = **4.77×**  (90 % BCa CI **3.60 – 6.98**, n = 33)
>
> The best unsigned geometric-mean fold error **any constant predictor could
> achieve** on measured rat brain-slice `Vu,brain`. The minimising constant is
> **39 mL·g brain⁻¹**.
>
> This is the price the engine pays for having **no per-compound tissue-capacity
> term at all**, quantified against measurement for the first time. It contains no
> engine quantity, no descriptor and no code version, and it therefore **does not
> move** under any resolution of the live `HD`-convention investigation, under any
> recalibration of `alpha`, or under any bug fix.
>
> *Measured reproducibility floor: **1.25×** (§3.4a; n = 2 replicate pairs, an
> order-of-magnitude bound, not a calibrated uncertainty).*

Read plainly, and this is the answer to *"what would a per-compound capacity term be
worth?"*: **the best constant still carries a 4.77× geometric-mean fold error, and
no amount of re-centring a constant can remove any of it.** A perfect per-compound
term could at best drive that to the assay's own 1.25× floor — so the headroom a
capacity term is competing for is roughly **4.77× down to 1.25×**, and the headroom
a *better constant* is competing for is **zero**.

### 0.2 The verdict

| | |
|---|---|
| **THE TEST (B)** | `GMFE(engine) / GMFE_floor ≤ 3.0` |
| `GMFE(engine)` | **57.50×** |
| `GMFE_floor` | **4.77×** |
| **ratio** | **12.05** |
| **90 % BCa CI on the ratio** | **[7.21, 18.93]** — entirely **above** 3.0 |
| **bucket** | **B-B — MISCALIBRATED** (prereg §4.1) |
| **verdict** | **FAIL** |
| **BW1 (reachability)** | **TRIGGERED** — with an interpretation note; see §4.1 |

The interval lies wholly on one side of the bar, so §3.6's resolution rule is
satisfied and the verdict stands; this is not B-C.

**Under prereg §4.1 and §8.5 this must appear in `codon validate` and must make it
exit 1. That registration has NOT been done by this run** — it requires editing
`tmesim/validation/scorecard.py`, which was outside this run's territory. It is the
outstanding action, and until it is done the failure is recorded here but not
enforced. **The required next action on the model is to constrain the coefficient —
not to widen the bar.**

### 0.3 What failed, in one sentence

The engine's implied unbound brain tissue distribution volume is **1.0 mL·g
brain⁻¹**, a value **39× below** the best constant available and **below every one
of the 33 retained measurements**, whose range is 3.2 – 4184 mL·g brain⁻¹.

---

## 1. The scored set — the §3.1 cascade, applied in the registered order

Applied mechanically to `data_cache/brain_slice/vu_brain.csv`
(sha256 `a00d75a0151fb4d8…`, recorded in full in the frozen JSON), whose
`INCLUSION_RULE.md` was frozen before extraction.

| step | rule | n retained | dropped |
|---|---|---:|---|
| SOURCE | every row in `vu_brain.csv` | 77 | — |
| **E2** | non-rat species excluded; species never pooled | 60 | 17 mouse rows |
| **E1** | values imported from another study excluded | 59 | oxycodone (`imported:Fridén 2011`) |
| **E5** | rows whose `Vu,brain` was **derived from a pH-partition model** excluded | 35 | 24 Fridén 2014 rows |
| **E6** | cyclic peptides excluded | **33** | SFTI-1, Kalata B1 |

**THE SCORED SET: n = 33.** Range **3.2 → 4184 mL·g brain⁻¹**, spread **1307.5×** —
the figure quoted throughout, as §3.1 requires. **The 9596× figure is not quoted as
the headline**; it describes the with-peptide set and appears only as a labelled
sensitivity arm.

Every exclusion is one the pre-registration fixed in advance on stated properties.
**No compound was excluded for how it scored, and none was added.** The cascade was
also run as a cross-check against the dataset's own `scoreable` flag: the two agree
exactly, and every one of the 33 retained rows carries `scoreable = yes`.

`kp_uu_brain_DO_NOT_SCORE.csv` was **not opened for scoring** (§3.4 E3). Test C
(`Kp,uu,cell`) was **not run** — 4 scoreable values, §3.8, NOT RUNNABLE.

**§8.2 — does the dataset separate its own compounds?** 10 of 33 rows report an SD;
median relative SD **11.4 %** (fold factor 1.114×). Measured spread 1308× against
3 × median uncertainty = 3.34× → **PASS**. Against 3 × the 1.25× reproducibility
floor = 3.75× → **PASS**. Both readings of the criterion clear it by ~350×
(interpretation note in §6.3).

**§8.1 structural requirement** — the retained set spans benzodiazepines,
phenothiazines, SSRIs/SNRIs, a taxane, a cardiac glycoside, a zwitterionic
anticonvulsant, an arylacetic-acid NSAID and 23 structurally unrelated PET
radioligands. Far more than two scaffolds. Satisfied.

---

## 2. GATE 1 — does the engine express this quantity? (prereg §3.3)

**GATE 1 CLEARS.** The engine has a commensurable quantity, so the outcome is not
B-D and BW3 is not triggered.

The four candidates §3.3 lists resolve as follows:

| candidate | verdict |
|---|---|
| `binding.tissue_binding_rate` | **rejected — wrong dimension.** 1/s. A rate is not a volume (`molecular/types.py:64`, `binding.py:32`) |
| `binding.plasma_free_fraction` | **rejected — wrong compartment.** Dimensionless, plasma |
| `two_route.cell_partition` (K) | **rejected as the scored object.** §3.2 forbids reading `Vu,brain` as a constraint on `K`; it is a component of the lump, not the lump |
| `ion_trapping` retardation `R` | **not shipped.** No shipped code path imports `molecular/ion_trapping.py`; the module is default-off |

**The quantity actually scored** is the engine's **tissue capacity term**, written
down explicitly in the derivation inside `tmesim/molecular/two_route.py`
(`retardation_factor`, lines ~669–684):

```
s_eff = (1 - phi) * 1  +  phi * K * (1 + K_b)
```

— matrix (interstitial) capacity normalised to 1, inclusion (cell) capacity
`K·(1+K_b)`, combined by an exact volume average because, in that module's own
words, *"the amount stored per unit volume at equilibrium is a volume average and
nothing else."* This **is** the lumped tissue content per unit unbound interstitial
concentration, which is what `Vu,brain` measures.

At the §2.3 frozen shipping state — `phi = 0.7892`, `K = cell_partition = 1.0`, and
`K_b` absent because `compile_drug.build_profile` never passes
`bound_to_free_ratio` — this evaluates to

```
s_eff = (1 - 0.7892) + 0.7892 x 1.0 x 1.0  =  1.000000   mL interstitial-equivalent per mL tissue
```

**exactly**, for every compound.

**The four gate items, each checked:**

1. *Ratio of tissue content to unbound interstitial concentration, same direction* —
   **holds.** `s_eff` is defined as stored amount per unit tissue volume referenced
   to the matrix-equivalent potential.
2. *Dimensions reduce to mL·g⁻¹ under a stated, pre-registered assumption* —
   **holds.** The assumption, written into the frozen JSON **before** any
   comparison was computed (§3.3 item 2, §8.4): **tissue density ρ = 1.00 g·mL⁻¹**,
   the water-equivalent convention. It requires no citation and cannot be tuned.
   → **engine `Vu,brain` = 1.000 mL·g brain⁻¹**.
   *Sensitivity:* at ρ = 1.04 g·mL⁻¹ the value is 0.962 mL·g brain⁻¹ — a 4 % move,
   far inside the 1.25× reproducibility floor. **The verdict is identical under
   either density.**
3. *Both unbound-referenced on the denominator* — **holds.** The engine's matrix
   concentration is the free interstitial concentration; there is no interstitial
   binding capacity in the convention.
4. *Equilibrium partitioning, not an active-transport steady state* — **holds.**
   `s_eff` is derived at local equilibrium and efflux ships default-off
   (`efflux_ratio=None` is bit-identical to omitting it).

**The prediction is constant, and this was verified rather than assumed.**
`two_route_diffusivity` was evaluated over 210 descriptor combinations spanning
logP −4…8, HBD 0…16 and hydrodynamic radius 0.15…3.0 nm. **One** distinct capacity
state appeared across all of them (`retardation_factor = 1.0`,
`cell_partition_effective = 1.0`, `phi = 0.7892`). Spearman ρ is therefore
**NOT COMPUTABLE**, exactly as §3.5 registered — not 0.00, and not a failure.

---

## 3. THE TEST (B)

```
GMFE(p)     = exp( mean_i | ln( p / Vu_i ) | )
GMFE_floor  = min over constants c of GMFE(c)      -- attained at the median of ln(Vu)
THE TEST    : GMFE(engine) / GMFE_floor  <=  3.0
```

The closed form for the floor (the minimiser of a mean absolute deviation in log
space is the median) was **checked numerically** against a 200 001-point dense scan
of the constant; the scan minimum exceeds the closed form by 1.3e-07, i.e. the
closed form is the minimum.

| quantity | value | beside it |
|---|---:|---|
| n | 33 rows (31 distinct compounds — see §6.4) | |
| measured range | 3.2 – 4184 mL·g brain⁻¹ | spread 1307.5× |
| engine constant | **1.000 mL·g brain⁻¹** | ρ = 1.00 g·mL⁻¹, declared in advance |
| `GMFE(engine)` | **57.50×** | *reproducibility floor 1.25×* |
| `GMFE_floor` | **4.77×** (CI 3.60 – 6.98) | *reproducibility floor 1.25×* |
| **ratio** | **12.05** | bar = 3.0 |
| **90 % BCa CI on the ratio** | **[7.21, 18.93]** | BCa, 20 000 resamples, seed 20260821, z₀ = −0.0827, a = −0.0328 |

**The interval lies entirely above 3.0. Bucket B-B — MISCALIBRATED. FAIL.**

The bar was registered as *forgiving by construction* (§3.5) on a wide-spread set.
It was missed by **4×**.

---

## 4. Withdrawal conditions (prereg §4.2)

### 4.1 BW1 — reachability: **TRIGGERED**

> *BW1 — the assumed value is outside the measured range entirely. The engine's
> implied `Vu,brain` lies outside `[min, max]` of the retained measurements.*

Engine **1.0** against the retained range **[3.2, 4184]** → **outside**.

**The claim BW1 withdraws** — stated in §4.2 and quoted here so the withdrawal is
legible:

> *The engine's treatment of tissue drug capacity is a physically defensible
> approximation, rather than a placeholder whose value was never checked against
> anything.*

**That claim is withdrawn, not narrowed.** §4.2 is explicit that this is a
reachability argument, *"decisive because it is not a statistical claim at all —
one comparison, no resampling."*

**An interpretation this run had to make, flagged rather than buried.** §4.2's
illustrative paragraph — *"the rat range spans 0.436 → 4184 mL·g⁻¹, so it straddles
unity. Any engine constant inside that interval clears BW1"* — quotes the range of
the **n = 35 pre-E6** rat set. BW1's operative sentence says *"of the **retained**
measurements"*, and the retained set is the post-cascade n = 33, whose minimum is
**3.2** because the 0.436 value is SFTI-1, a cyclic peptide that **§3.4 E6 excludes**.
**Unity is inside the illustrative range and outside the retained one.** I applied
the operative sentence, because "retained" has a single fixed meaning in this
document and the following paragraph is labelled an illustration. Recorded so the
coordinator can overrule it: **under the illustrative range BW1 would NOT trip.**
It changes nothing about THE TEST, which fails independently and by 4×.

### 4.2 BW2 — proxy absorption by the membrane leg: **NOT EVALUABLE**

BW2 requires B-S9 — Spearman ρ between measured `Vu,brain` and Crippen logP.
**B-S9 is NOT COMPUTABLE:** the frozen dataset carries no SMILES column, and 23 of
the 33 retained rows are radiolabelled PET tracers. Resolving structures now would
be a curation step taken **after** the measured column was read, which §7.1
prohibits in substance if not in letter. I declined to do it.

BW2 was in any case already weakened by §6.3 to **hypothesis-generating only**,
because `descriptor_falsification.md` found the shipping descriptors degrade
materially. **BW2 is therefore neither triggered nor cleared. It remains open**, and
if the coordinator wants it closed, the honest route is to freeze a SMILES column
into the dataset under a dated amendment and re-run — not to source structures
inside this run.

### 4.3 BW3 — units/role error: **NOT TRIGGERED**

GATE 1 found a commensurable **partition**, not a rate wearing a partition's
clothes. This is not the `_run_ps_scale` signature.

### 4.4 What did NOT trigger withdrawal, and is not being counted as evidence

Per §4.3: the failure to predict per-compound variation (arithmetic, §3.5); the
large `GMFE_floor` (a measured finding, not a failure to clear a bar); the mouse
column (§3.4 E2, reported separately below); `Kp,uu,brain` (E3, not computed); the
single worst outlier (B-S6 — reported, **not removed**).

**One §4.3 escape hatch was checked and does NOT apply.** §4.3 exempts *"a miss
concentrated in the base stratum"*, because the engine has no ionisation term in
this path. The miss is **not** confined to bases: the non-base stratum's ratio is
**6.58**, still more than double the bar (§5, B-S3). The failure is general, not an
ionisation artefact.

---

## 5. Secondaries (prereg §3.7) — none replaces THE TEST

*The 1.25× reproducibility floor prints beside every fold-error figure, per §3.4a.
It is estimated from **n = 2 replicate pairs** — paroxetine 714/876 = 1.23×,
diazepam 17.8/20 = 1.12× — and is an order-of-magnitude bound, not a calibrated
uncertainty.*

| # | measure | result |
|---|---|---|
| **B-S1** | `GMFE_floor` + minimising constant | **4.77×** (90 % BCa CI 3.60 – 6.98); constant **39 mL·g brain⁻¹** |
| **B-S2** | centring ratio (engine ÷ optimal constant) | **0.0256 — the engine is 39.0× low.** Diagnostic, not an acceptance criterion |
| **B-S3** | stratified by ion class | see below |
| **B-S4** | direction of bias | **one-directional. Below 33 of 33 measurements.** Signed geometric bias 0.0174 (57.5× low) |
| **B-S5** | within-3× / within-10× coverage | **0.0 % / 21.2 %** |
| **B-S6** | worst-case fold error | **4184×**, `[11C]Sertraline` (measured 4184 mL·g brain⁻¹). Reported, **not removed** |
| **B-S7** | reproducibility floor | **1.25×**, printed beside every figure above |
| **B-S8** | Spearman ρ (engine) | **NOT COMPUTABLE** — zero-variance prediction, verified over 210 descriptor inputs. Not 0.00, not a failure |
| **B-S9** | ρ(measured `Vu,brain`, Crippen logP) | **NOT COMPUTABLE** — no SMILES in the frozen dataset (§4.2 above) |
| **B-S10** | sensitivity arms | see below |

### B-S3 — ion-class strata (diagnostic only)

Only 8 of the 33 rows carry an ion-class assignment; the 25 PET-tracer and CE-123
rows do not. **No stratum clears n ≥ 10, so none carries a verdict, and strata were
not used to select any subset for the headline.**

| stratum | n | `GMFE(engine)` | `GMFE_floor` | ratio | median `Vu` |
|---|---:|---:|---:|---:|---:|
| base | 5 | 295.9× | 4.22× | **70.15** | 714 |
| non-base | 3 | 10.4× | 1.58× | **6.58** | 14.1 |
| unassigned | 25 | 50.9× | 4.42× | **11.51** | 39 |

The miss is **largest** in bases — the pre-registered expected direction (§3.7
B-S3) — but it is **not confined** to them (§4.4 above).

### B-S10 — sensitivity arms. Reported; **never** the verdict

| arm | n | `GMFE_floor` | `GMFE(engine)` | ratio | spread |
|---|---:|---:|---:|---:|---:|
| **HEADLINE (E1, E2, E5, E6 applied)** | **33** | **4.77×** | **57.50×** | **12.05** | **1308×** |
| + oxycodone (E1 relaxed) | 34 | 4.87× | 53.24× | 10.94 | 1308× |
| + cyclic peptides (E6 relaxed) | 35 | 5.16× | 53.92× | 10.45 | 9596× |
| + model-derived rows (E5 relaxed) | 57 | 3.53× | 42.97× | 12.17 | 2461× |
| mouse arm (E2 relaxed — **separate**, never pooled) | 17 | 7.57× | 23.57× | 3.11 | 1608× |

**Every arm fails.** No relaxation of any exclusion rescues the test, and no
exclusion is doing load-bearing work for the verdict. Two notes:

* The **with-peptides** arm has a minimum of 0.436 mL·g brain⁻¹, so **BW1 would not
  trigger on it** — the only place any arm disagrees with the headline, and it is
  recorded prominently here per §3.4 E6.
* The **mouse** arm at ratio 3.11 sits just above the bar. It is a separate species
  arm, is not pooled, and carries no verdict (§3.4 E2). It is not evidence that the
  rat result is marginal.

---

## 6. Things this run had to interpret — read this section first if you read only one

The pre-registration is unusually tight, but it is not gapless. Every place I made
a call the plan did not cover is listed here rather than absorbed into the numbers.

### 6.1 Which engine quantity GATE 1 selects — **the largest call in the run**

§3.3 lists four candidates and rejects three by inspection; §8.4 then hands *"the
engine quantity chosen as commensurable"* to the executor as a slot to fill. The
quantity I selected — the two-route **capacity term `s_eff`** — **is not on §3.3's
list.** I selected it because it is the only object in the engine that is the
*lumped* tissue-content-per-unbound-interstitial-concentration that §3.2 says
`Vu,brain` constrains, and because §3.5's own reasoning (*"`cell_partition` is fixed
at 1.0"*, *"retardation was built and refuted"*) is reasoning about precisely the two
inputs of `s_eff`.

**Why the choice does not move the number.** Under the §2.3-frozen shipping state
**every** dimensionless partition/capacity quantity in the engine equals exactly
**1.0**: `cell_partition` = 1.0; `two_route.retardation_factor` = 1.0 (`K_b` never
supplied); `ion_trapping` retardation = 1.0 (module not shipped); `s_eff` = 1.0.
**The verdict is invariant to which of them a reader prefers.** I label the scored
object `s_eff` rather than `K` specifically so that §3.2's prohibition is respected:
**nothing here constrains `cell_partition`, which remains unconstrained (§3.8).**

**The alternative reading, and what it would change.** A reader could hold that a
capacity fixed at 1 by a *normalisation convention* is not a quantity the engine
"expresses" at all, and that GATE 1 therefore fails → **bucket B-D, NOT TESTABLE,
reason *the model cannot express the measured quantity***. I did not take that
reading, because `s_eff` passes all four gate items on its own terms. **Both
readings are non-passes and neither supports the claim**, so the choice does not
launder anything — but the coordinator should know it was available.

### 6.2 The density assumption

§8.4 lists *"the density/conversion assumption used"* as a slot to fill, and §3.3
item 2 requires it be written down before the comparison. I fixed **ρ = 1.00
g·mL⁻¹** and wrote it into the frozen JSON before running either scoring stage.
I chose the water-equivalent convention because it needs no citation and so cannot
be tuned. It is immaterial: at ρ = 1.04 the engine value is 0.962 and every verdict
in this document is unchanged.

### 6.3 §8.2's "3 × the median reported measurement uncertainty"

The phrase does not say whether the uncertainty is absolute or relative, or how a
spread in fold units compares to it. Only 10 of 33 rows report an SD at all. I
computed it both ways and **both clear by roughly two and a half orders of
magnitude**, so I did not need to resolve the ambiguity and did not.

### 6.4 Duplicate compounds were **not** de-duplicated

The 33 rows are **31 distinct compounds**: paroxetine appears twice (714 Loryan,
876 Fridén) and diazepam twice (17.8, 20) — the two pairs that *are* the 1.25×
reproducibility floor. The dataset's rule R4 keeps duplicates as separate rows, and
**the §3.1 cascade registers no de-duplication step**, so I did not add one. I have
deliberately **not** computed a de-duplicated arm, because inventing an unregistered
arm after seeing the result is the failure mode this exercise exists to prevent.
Flagging it so the coordinator can register one if it matters — it cannot plausibly
move a 12.05 to below 3.0.

### 6.5 A leak in the process, disclosed

§7.3 requires the frozen JSON be written *"before the scoring script may read any
measured column."* The JSON was written before either scoring stage ran — but while
**verifying the cascade counts** I had already printed the measured `Vu,brain`
column to confirm that the retained n and spread matched §3.1's stated 33 / 1308×.
**No exclusion, threshold, statistic or engine-quantity choice was decided by that
inspection** — the exclusion codes are transcribed verbatim from §3.4 and the
cascade is mechanical — but the ordering was not clean and pretending otherwise
would defeat the point. Recorded.

### 6.6 A file created outside the granted territory

`results/sweep/transport_anchor_prereg.json` did not exist and §7.3 mandates it. I
created it. It is new, not an edit, and it is the binding copy of the bars the plan
requires; flagging it because it is a second artefact where one was authorised.

---

## 7. An observation the plan did not anticipate — and it changes the shape of the fix

**Registered here as an observation, not as a result, and it did not affect the
verdict.**

§4.2 frames the claim BW1 kills as being about *"a placeholder whose value was never
checked against anything."* What GATE 1 actually found is narrower and stranger:
the engine's tissue capacity is **1.0 by construction, not by choice**. From
`two_route.py`'s own derivation:

> *"This module's published output is a CONDUCTANCE, computed against a convention
> in which the no-binding capacity is 1."*

Nobody set this value. It is the normalisation implied by publishing a conductance
and declining to publish a capacity. That has two consequences for §4.1's *"required
next action is to constrain the coefficient":*

1. **Neither of the two available knobs is a clean one.** `s_eff = (1-φ) + φ·K·(1+K_b)`
   has exactly two adjustable inputs, and both are blocked:
   * **`K` (`cell_partition`) is not a capacity-only knob.** It also sets the
     inclusion conductivity in `two_route_factor` (`sigma_c = cell_partition *
     cytoplasm_D`). Raising `K` to lift the capacity would move `D_eff` and re-open
     the Pruijn HT-29 fit — and §3.2 independently forbids using `Vu,brain`, a
     lumped quantity, to set `K`, which is one component of the lump. **This test
     does not license moving `K` and must not be cited as if it did.**
   * **`K_b` (intracellular retardation) is the per-compound capacity term the
     engine lacks** — and it was built and **refuted** on Pruijn HT-29
     (`results/intracellular_retardation.md`: R² 0.718 → 0.606; matched-MW rank
     agreement 77.4 % → 51.6 %). So the term this FAIL points at and the term the
     earlier work rejected are the **same term**, and the two datasets point in
     opposite directions.

   That is a real tension and resolving it is a model decision, not an execution
   decision. I flag it and stop.
2. **B-S4 says the offset is one-directional** — below 33 of 33 measurements,
   signed geometric bias 57.5× low. Under §4.1 that is the *"single corrected
   value"* shape rather than the *"needs a per-compound term and a full blast-radius
   re-run"* shape. **But B-S1 says a single corrected value buys at most a fall from
   57.5× to 4.77×**, and 4.77× is still nearly 4× the assay's own noise floor. **A
   re-centred constant would pass THE TEST by construction and would still be wrong
   by 4.8× on average.** That is the whole content of §3.5a's insistence that
   `GMFE_floor` be a headline, and it is why re-centring must not be reported as a
   fix.

---

## 8. What may and may not be claimed

**Claimable.** The engine's assumed tissue drug capacity is miscalibrated: it lies
below every one of 33 measured rat brain-slice `Vu,brain` values and 39× below the
best constant available, and it misses the pre-registered baseline-relative bar by
4×, with the 90 % interval entirely on the failing side. The best any constant could
do on this measurement is `GMFE_floor` = 4.77×.

**Not claimable.** §3.9's PASS sentence does not apply and must not be quoted.
Nothing here is evidence about: the transport solver (never entered), any solved
concentration field, vessel spacing, transvascular delivery, interstitial
convection, any absolute kill fraction, anything clinical or diagnostic, `D_eff`
accuracy, whether the two-route law generalises, or human tumour (§3.10 — rat, brain,
wrong pH environment, and *the anchor is weakest exactly where the quantity matters
most*).

**`cell_partition` (K) remains unconstrained by any measurement** (§3.2, §3.8). This
FAIL does not constrain it either. Test C is still NOT RUNNABLE on 4 values.

---

## 9. Reproduction

```
./.venv/bin/python experiments/brain_slice_anchor.py --stage prereg   # freeze
./.venv/bin/python experiments/brain_slice_anchor.py --stage floor    # GMFE_floor
./.venv/bin/python experiments/brain_slice_anchor.py --stage engine   # THE TEST
```

Bootstrap: BCa, 20 000 resamples over compounds, seed **20260821**, 90 % interval,
identical resamples for `GMFE_floor` and for THE TEST ratio (§3.6). Both arms of the
ratio were computed in the **same engine state**, commit
`544734e469fa5a22f710248cd10a73db394e59d3`, per §9.5 rider 1. No parameter,
exclusion, threshold or descriptor choice was changed by looking at a residual.

**Rider 3 disclosure — the working tree was not clean.** Sibling agents were live in
this repository during the run and the tree carries uncommitted edits, including one
to `tmesim/validation/propagation.py`. That file is **not** in §2.3's frozen set, is
not on the `molecular.two_route` / `compile_drug` path this test uses, and every
§2.3 constant was asserted against live code at both the `prereg` and `engine`
stages and matched. `GMFE_floor` is in any case engine-independent, and both arms of
THE TEST ratio were computed inside a single process in a single engine state. I do
not believe rider 3 was breached, but the tree was dirty and saying so is cheaper
than asserting it was not.

**Nothing in `tmesim/`, `data_cache/`, `docs/` or
`results/transport_anchor_prereg.md` was modified by this run.** No assertion was
deleted or weakened: `./.venv/bin/python -m pytest -q` → **1076 passed, 0 failed**
(555.85 s), the unchanged baseline.

Three files were added and nothing else:
`experiments/brain_slice_anchor.py`, `results/sweep/transport_anchor_prereg.json`,
and this document.
