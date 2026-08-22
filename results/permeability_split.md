# Separating ENTRY from SPREAD — the permeability split, and the thing it found instead

**Status: RUN AND MEASURED.** Driver: `experiments/permeability_split.py`.
Data: `results/permeability_split_data.json`. Code:
`tmesim/treatment/transvascular.py`. Tests: `tests/test_permeability_split.py`
(16 pass; four defects planted and watched to bite, §6).

**Nothing shipped was changed.** `PDEConfig.transvascular_exchange` is still
`False`, `transvascular.SPLIT_DEFAULT` is `False`,
`drugs.PS_SCALE_S_INV` is still `1.0`, and `calibrate_alpha()` is asserted live
at `22.223477711172563` at the top and the bottom of every stage — the run
raises if any of them has moved. **No parameter was fitted to any anchor.** The
entry-scale sweeps below construct hypothetical compounds and adopt none of
them.

---

## 0. Verdict up front

**The split does not close the accuracy gap, and the reason is structural
rather than empirical. It does invert a per-cell verdict on real tissue, at the
shipped parameters, with no scaling required.**

| | question | answer |
|---|---|:--|
| **A** | is entry conflated with spread? | **the brief had the polarity backwards** — entry contributes *exactly nothing* to shape; `D_eff`/`k_loss` does both jobs |
| **B** | does the split move the Kyle gradient? | **NO — and it cannot.** Best available move −0.15 %, *away* from the measurement, against gaps of 3.67× and 13.7× |
| **C** | does the split invert a per-cell verdict? | **YES.** Daunorubicin's dominant limiter flips DOSE → DELIVERY on the flagship section at the shipped `PS` |
| **D** | does the split need an unmeasured parameter? | **it trades an unbounded one for a saturating one.** It does not introduce a measured `P`; nobody has one |
| **E** | *(not the assignment)* what IS the accuracy gap? | **an error of ≤ 2.43× in ONE length**, `λ = sqrt(D_eff/k_loss)` |

**E is the most valuable thing in this file and it is a by-product.** It is
reported first, in §1, because a 14× failure and a 2.4× failure are different
problems and the project has been carrying the former.

---

## 1. The by-product that matters most: the 14× is a ≤ 2.43× error in one length

`results/kyle2004_penetration.md` §6 records the engine as **~14× too shallow**
against Kyle 2004 (daunorubicin near:far 3.65 against a measured 50 ± 15) and
**~3.7× under-separating** the doxorubicin/daunorubicin pair. Those are ratio
comparisons, and a ratio across a 140 µm slab is *exponential* in the decay
length. So a modest error in one length presents as a large error in the ratio.

**Reproduction first, before any argument.** On the arm-A slab exactly as
`experiments/isomer_penetration_map.py` builds it — 140 µm thick, `dx = 2 µm`,
uniform `ecm_density = 0.0875`, 1 h exposure window, convection OFF, exchange
OFF, routed through `sim.engine.simulate_response` — doxorubicin's near:far
exposure ratio is **2.9977**, against the **3.00** published in
`results/kyle2004_penetration.md` §6.1. The construction is reproduced before
it is reasoned on.

**The map.** On that same slab and configuration, sweeping the ratio against
`λ = sqrt(D_eff / (clearance + binding))` alone, by two *independent* routes to
each λ — vary `k` at fixed `D`, and vary `D` at fixed `k`:

| λ (µm) | ratio via `k` | ratio via `D` | disagreement |
|---:|---:|---:|---:|
| 20.000 | 205.60 | 243.60 | 1.185× |
| 25.000 | 64.65 | 73.55 | 1.138× |
| 27.000 | 45.76 | 51.42 | 1.124× |
| 30.000 | 29.65 | 32.80 | 1.106× |
| 35.000 | 16.92 | 18.33 | 1.084× |
| 40.000 | 11.09 | 11.82 | 1.066× |
| 45.000 | 7.98 | 8.39 | 1.051× |
| 50.000 | 6.14 | 6.38 | 1.039× |
| 60.000 | 4.16 | 4.24 | 1.019× |
| **72.487** | **2.9978** | **2.9978** | 1.000× |
| 90.000 | 2.24 | 2.19 | 1.021× |

**Where the measurement falls on the model's own coordinate** (same slab,
`dx = 2 µm`, ecm 0.0875, 1 h, convection OFF, exchange OFF):

| compound | Kyle Table 1 (fluorescence) | implied λ | shipped λ | fold error in λ |
|---|---:|---:|---:|---:|
| doxorubicin | 11 ± 5 | **40.15 µm** | 72.49 µm | **1.81×** |
| epirubicin | 12 ± 3 | 39.22 µm | 72.49 µm | 1.85× |
| daunorubicin | 50 ± 15 | **26.55 µm** | 64.40 µm | **2.43×** |

**The ~14× absolute-gradient failure and the ~3.7× fold failure are, in this
geometry, a single length that is between 1.8× and 2.4× too long.** That is a
findable error. "The map is 14× too shallow" is a true sentence about a ratio
and a misleading sentence about the model.

**Four things this does NOT say, and they matter:**

1. **λ cannot tell you whether the error is in `D_eff` or in `k_loss`.** It is
   their ratio. The two-route columns above are the only handle this
   construction gives and they agree to within **1.185×** at the far end and
   1.00× at the shipped point — far too weak to separate `D` from `k`. A
   sibling is decomposing the accuracy gap in `experiments/accuracy_gap_*`;
   this file hands the question over, it does not answer it.
2. **The ratio is very nearly, but not exactly, a function of λ alone.** The
   two routes differ by up to 18.5 % at λ = 20 µm. The 1 h window is ~18 loss
   time constants at the shipped `k` but fewer when `D` is what moves, so the
   transient is not identically resolved along the two routes.
3. **Nothing here is a fit.** No constant was changed and no λ was adopted. The
   table maps the model's own readout onto the model's own coordinate and reads
   off where the measurement sits. `calibrate_alpha()` is asserted live at
   `22.223477711172563` before and after.
4. **Doxorubicin is the α anchor.** Its row is weaker evidence than
   daunorubicin's for the reason `results/kyle2004_penetration.md` §7 already
   states.

---

## 2. The entanglement, measured — and the brief had the polarity backwards

The assignment stated that "one parameter is doing both jobs" — entry setting
the level, spread setting the shape, conflated. **The measurement says the
opposite, and it says it exactly.**

Under the shipped one-way law

```
dC/dt = div(D_eff grad C) - div(vC) + PS*Cp - (clearance + binding)*C
PS    = vascular_density * bbb_penetrance * PS_SCALE_S_INV
```

`PS` enters a linear PDE only as a multiplicative factor on the source, so the
solved field is **positively homogeneous of degree 1** in the entry term:
`C(x, t; c·PS) == c · C(x, t; PS)` for every `c > 0`.

Measured, same slab and configuration as §1 (140 µm, `dx = 2 µm`, ecm 0.0875,
1 h window, convection OFF, **exchange OFF**), sweeping `PS` by scaling
`bbb_penetrance` — which reaches the engine through `delivery_source` and
nowhere else, so scaling it *is* scaling `PS` and touches neither `D_eff` nor
`clearance` nor `binding`:

| `c` | near:far | AUC at the wall (µM·s) |
|---:|:---|---:|
| 1e−06 | 2.997687165949104 | 1.850775e−02 |
| 1e−03 | 2.9976871659491064 | 1.850775e+01 |
| 1e+00 | 2.9976871659491016 | 1.850775e+04 |
| 1e+03 | 2.9976871659490905 | 1.850775e+07 |
| 1e+06 | 2.9976871659491144 | 1.850775e+10 |

**Twelve decades of entry term. Twelve decades of level. Zero change in shape:**
ratio spread `7.99e−15`, worst whole-profile shape deviation `4.45e−14` — both
at floating-point round-off.

> **The parameter doing two jobs is `D_eff`/`k_loss`, not `PS`.** Entry sets the
> level and contributes *exactly nothing* to the gradient. Spread sets the
> shape *and*, once `PS` has fixed the scale, the perivascular level. The brief
> was wrong about which way round it runs, and this file records that rather
> than quietly answering a different question.

This is also why the defect was invisible for so long: `tmesim/validation/ps_scale.py`
already notes that the GBM penetration anchor is `L = sqrt(D_eff/k)`, in which
`PS` does not appear. The homogeneity result above is the general statement of
that observation — **no shape-based anchor of any kind can constrain the entry
term.**

---

## 3. The split does not move the Kyle gradient, and cannot

The gradient-driven exchange law `S = PS·(Cp − C)` — algebraically the same
source minus a first-order loss at the same rate — is the split. It already
exists (`treatment.drugs.exchange_coefficient`,
`PDEConfig.transvascular_exchange`, wired in `sim/engine.py`) and ships OFF.

Same slab, same configuration, **exchange ON**, sweeping the same entry scale:

| `c` | wall `PS` (1/s) | near:far | AUC at the wall (µM·s) |
|---:|---:|---:|---:|
| 1e−02 | 2.0568e−04 | 2.9977 | 1.848707e+02 |
| 1e−01 | 2.0568e−03 | 2.9976 | 1.830300e+03 |
| 1e+00 | 2.0568e−02 | 2.9972 | 1.664552e+04 |
| 1e+01 | 2.0568e−01 | 2.9952 | 8.732947e+04 |
| 1e+02 | 2.0568e+00 | 2.9934 | 1.517310e+05 |
| 1e+03 | 2.0568e+01 | 2.9931 | 1.638057e+05 |

**Largest move available: −0.15 %, in the SHALLOWER direction,** against a gap
of 3.67× for doxorubicin (11 / 3.00) and 13.7× for daunorubicin (50 / 3.65).
The ratio saturates at 2.9931 as the wall opens without bound.

**Why, and this is the sentence that should stop anyone trying this again:**

> `PS` is supported only where `vascular_density > 0`. Kyle's gradient is
> measured where it is not. The exchange term adds a sink *at the vessels* and
> leaves the tissue interior governed by `sqrt(D_eff/k)` exactly as before. The
> split is **structurally incapable** of reaching the interior, so it cannot be
> part of the gradient error at any magnitude.

Reason #1 in the assignment — "it is probably part of the ~14× accuracy gap" —
is **refuted**. It is not part of it, and §1 says what is instead.

---

## 4. It DOES invert a per-cell verdict, on real tissue, at the shipped `PS`

Flagship 10x Xenium human GBM section, `dx = 12 µm`, 424 × 675,
`repair = 1`, `dt = 10 s`, 6 h window, **convection ON**, through
`simulate_response(..., attribute=True)`. **No scaling of anything** — the
compiled `DrugParams` as shipped, exchange OFF versus exchange ON.

| compound | law | kill fraction | unreached / tumour | delivery : dose | dominant limiter |
|---|:--|---:|---:|---:|:--|
| doxorubicin | OFF | 0.045817 | 0.3668 | 0.384 : 0.616 | dose |
| doxorubicin | **ON** | 0.039140 | **0.4774** | **0.497 : 0.503** | dose |
| epirubicin | OFF | 0.045817 | 0.3668 | 0.384 : 0.616 | dose |
| epirubicin | **ON** | 0.039140 | 0.4774 | 0.497 : 0.503 | dose |
| **daunorubicin** | OFF | 0.044263 | 0.4098 | 0.429 : 0.571 | **dose** |
| **daunorubicin** | **ON** | 0.037502 | **0.5209** | **0.541 : 0.459** | **delivery — INVERTS** |
| mitoxantrone | OFF | 0.025113 | 0.7690 | 0.789 : 0.211 | delivery |
| mitoxantrone | **ON** | 0.021124 | 0.8255 | 0.843 : 0.157 | delivery |

**Daunorubicin's headline mechanism verdict inverts.** Under the shipped
one-way law the run reports that most survivors are **DOSE-limited** — the drug
arrived and there was not enough of it. Under the split, on the same tissue at
the same parameters, most survivors are **DELIVERY-limited** — the drug never
arrived. Those are opposite clinical readings and they lead to opposite
actions.

**Per-cell, at the shipped `PS`:**

| compound | cells changing label | of which DOSE → DELIVERY |
|---|---:|---:|
| doxorubicin | 2,757 / 40,887 (6.743 %) | 2,598 |
| epirubicin | 2,757 / 40,887 (6.743 %) | 2,598 |
| daunorubicin | 2,766 / 40,887 (6.765 %) | 2,605 |
| mitoxantrone | 1,400 / 40,887 (3.424 %) | 1,305 |

The other transitions are `KILLED → SURVIVED_UNDERDOSED` (55–124) and
`KILLED → SURVIVED_UNREACHED` (35–40). The direction is uniform: the split
moves cells *away* from "we killed it" and *toward* "the drug never got there".

**This is the whole remaining case for the split**, and it is the one the
landing page already carries as its first limit. It is not a caveat. It is a
reported conclusion changing sign under a law the engine ships in the OFF
position.

### 4.1 And the deeper reason: the shipped law makes every absolute verdict unbounded

Same section and configuration, doxorubicin, sweeping the entry scale by `c`:

| `c` | law | kill fraction | mean AUC (µM·s) | KILLED | UNREACHED | UNDERDOSED |
|---:|:--|---:|---:|---:|---:|---:|
| 1e−02 | OFF | 0.000042 | 4.8716e+02 | 1 | 23,811 | 0 |
| 1e−02 | ON | 0.000042 | 4.8617e+02 | 1 | 23,811 | 0 |
| 1e−01 | OFF | 0.003654 | 4.8716e+03 | 87 | 23,725 | 0 |
| 1e−01 | ON | 0.003486 | 4.7753e+03 | 83 | 23,729 | 0 |
| 1e+00 | OFF | 0.045817 | 4.8716e+04 | 1,091 | 8,735 | 13,986 |
| 1e+00 | ON | 0.039140 | 4.0790e+04 | 932 | 11,368 | 11,512 |
| 1e+01 | OFF | 0.287964 | 4.8716e+05 | 6,857 | 157 | 7,014 |
| 1e+01 | ON | 0.154124 | 1.7748e+05 | 3,670 | 196 | 19,865 |
| 1e+02 | OFF | 0.502688 | 4.8716e+06 | 11,970 | 141 | 12 |
| 1e+02 | ON | 0.219343 | 2.8661e+05 | 5,223 | 157 | 15,720 |
| 1e+03 | OFF | 0.531665 | 4.8716e+07 | 12,660 | 140 | 1 |
| 1e+03 | ON | 0.233874 | 3.1552e+05 | 5,569 | 154 | 13,355 |

Read the OFF column of mean AUC: `4.8716e+02, 4.8716e+03, 4.8716e+04,
4.8716e+05, 4.8716e+06, 4.8716e+07`. **Exactly proportional over five decades,
on real tissue with convection ON** — the homogeneity result of §2, surviving
the full production path with advection, a real vascular field and a real
tissue mask. The ON column saturates: `4.8617e+02 → 3.1552e+05`, a factor of
649 for a factor of 100,000 in entry.

**Why that is a correctness problem rather than a scaling nuisance.** Per-cell
DOSE / DELIVERY / RESISTANCE calls are thresholds on an *absolute* exposure
(`E / IC50`). Under the shipped law every absolute exposure is proportional to
`PS`, and `PS` carries an **implicit, underived factor of exactly 1 s⁻¹**
(`drugs.PS_SCALE_S_INV`; the dimensional audit is
`results/transvascular_exchange_revalidation.md` §6.2). So the per-cell map is
arbitrary up to an unvalidated multiplicative constant, and the sweep above
shows it: the same tissue and the same compound go from **23,811 unreached
cells (every viable tumour cell on the section) to 140** across five decades of
a number nobody has measured.

Formally, the elasticity of exposure with respect to the entry term is

* one-way: **exactly 1**, at every `PS`, without bound;
* split: `k_loss / (PS + k_loss)`, which → 0 as the wall opens.

(`transvascular.entry_elasticity`, asserted in `tests/test_permeability_split.py`.)

**What stays safe under the shipped law.** Everything that is a *ratio within
one map* — "this region got N× less than that region", the compound rank by
half-distance, the shape of the reach curve. Those are exactly invariant to the
unmeasured scale, by §2, to `4.45e−14`. That invariance is a real asset and it
is why the map claim survived while the absolute claims did not.

---

## 5. Honest scope — the split does not buy a measured entry parameter

A dimensionally honest entry term is

```
PS [1/s]  =  P [µm/s] × (S/V) [1/µm],      S/V = 2·φ / R  for cylinders
```

`tmesim/treatment/transvascular.py` writes that decomposition down
(`ps_from_permeability`, `surface_to_volume_per_um`,
`implied_permeability_um_s`). **This repository has neither factor.**

Inverting the shipped law against it, `φ` cancels exactly and the implied wall
permeability is a per-compound constant:

| compound | `bbb_penetrance` | implied `P` at R = 2 µm | R = 4 µm | R = 8 µm |
|---|---:|---:|---:|---:|
| doxorubicin | 0.02057 | 0.020568 µm/s | 0.041136 | 0.082271 |
| epirubicin | 0.02057 | 0.020568 | 0.041136 | 0.082271 |
| daunorubicin | 0.02542 | 0.025424 | 0.050848 | 0.10170 |
| mitoxantrone | 0.02252 | 0.022521 | 0.045042 | 0.090084 |

**These are not measurements and not predictions. They are what the shipped law
would have to be claiming if it were a permeability law**, and every factor in
them is contested:

* `bbb_penetrance` is a dimensionless brain:plasma **partition** coefficient,
  and `core/config.py`, `reconstruct/barrier.py` ("It is NOT a permeability
  with physical units") and `molecular/permeability.py` all say so;
* the `1 s⁻¹` is undocumented and nothing in the repo derives it;
* `R` is the caller's — **this repository does not carry a measured vessel
  radius**, which is why `surface_to_volume_per_um` has no default radius.

`tmesim/validation/ps_scale.py` does hold two independently *measured* rodent
brain PS products (methotrexate 3.3e−5 /s, PMID 524251; temozolomide 1.6e−3
/s, PMID 39690422) and a calibrated scale of 0.16122 1/s fit to them. That
calibration is not adopted here and is not adopted anywhere; and it cannot fix
the ordering, because the two anchors disagree by 8.70× and the engine's
cross-drug permeability ordering is wrong *in both directions* (that module
states both).

**So, plainly:**

> The split does **not** replace an unmeasured entry parameter with a measured
> one. It replaces an **unbounded** unmeasured parameter with a **saturating**
> one. The gain is in falsifiability — `C ≤ Cp` holds by construction for every
> drug at every dose, which the one-way law violated by 7.80× against a measured
> 0.61× for 5-fluorouracil (`results/breast_kp_validation.md`) — and not in
> accuracy. Anyone who reads §4 as "the split makes delivery more accurate" has
> read it wrong.

At the shipped wall `PS = 2.0568e−02 1/s` and doxorubicin's
`D_eff = 26.459 µm²/s`, `k_loss = 5.0355e−03 1/s`: the perivascular decay
length is 72.49 µm under the one-way law and **32.15 µm under the split**,
which is 2.7× the flagship section's `dx = 12 µm` and therefore resolved on
that grid at the shipped `PS`. It stops being resolved above `PS ≈ 0.179 1/s`.
The discretisation question is owned by
`experiments/boundary_layer_resolution.py`, not by this file.

---

## 6. The checks, and the defects planted to prove they are not blind

`tests/test_permeability_split.py`, 16 tests, all passing. The PDE checks use a
25 × 3 slab at `dx = 8 µm`, uniform `ecm_density = 0.3`, `dt = 2 s`, 1200 s
window, convection OFF, routed through `simulate_response` so the flag reaches
the solver by the one production wiring point.

**Four defects were planted and each was watched to bite.**

| # | defect | verbatim result |
|---|---|---|
| 1 | `delivery_pair` ignores `split=True` and returns `(source, None)` — the INERT defect | **3 failed**: `assert None is not None` / `AssertionError: split=True must build an exchange field` |
| 2 | `SPLIT_DEFAULT` flipped to `True` — the shipped-default defect | **2 failed**: `assert True is False`, and `split=False` stopped being bit-identical |
| 3 | source and back-flux built from different `ps_scale` — the drift `drugs.exchange_coefficient` warns about in prose | **1 failed**: `test_split_on_carries_the_same_ps`, `assert np.array_equal(ex * cp, src)` |
| 4 | `exchange_k` dropped at the **engine** boundary (monkeypatched `simulate_drug_transport`) — the real defect `tmesim/validation/propagation.py` documents | **2 failed**: `the split did not change the gradient SHAPE across 8 decades of PS (max relative shape change 1.39e-15)`, and `the split law exceeded the plasma bound: 3.61839e+10 > 58300` |

**Defect 4 is the one that mattered and the reason the checks are written the
way they are.** A test that compared "flag ON" against "flag OFF" would have
been **green** under defect 4 — it would have compared a run to itself, and
`tmesim/validation/propagation.py` records that exact failure happening for
real (`auc_sum 1.554447793e+09` bit-identical between the two arms). The check
instead asserts the *shape-invariance property directly*, which an inert flag
cannot satisfy.

**One check is blind by design and it is recorded rather than hidden.**
`test_one_way_law_is_exactly_homogeneous_in_the_entry_term` **passed** under
defect 4. It must: it is an OFF-arm test and an inert flag makes the ON arm
into the OFF arm, which is precisely what it asserts. It is not coverage of the
split and is not counted as such.

**A fifth failure was real and is recorded because it was mine, not the
physics'.** The plasma-bound check first failed at `58254.2 > 1200` — because
it hard-coded `Cp = 1 µM` while the default `PKConfig(plasma_peak=50,
half_life=14400, n_doses=1)` gives a window mean of 48.586 µM over 1200 s. The
bound holds at 58254.2 against 58302.8 (0.08 % below). The test now computes
the bound from the run's own `PKConfig` on the solver's own step grid.

`tests/test_propagation.py` passes (54 tests) with the new module and
experiment in the tree: `tmesim/treatment/transvascular.py` never calls
`simulate_drug_transport`, and `experiments/permeability_split.py` routes every
solve through `simulate_response`, so neither is a wiring point and
`tmesim/sim/engine.py` remains the only production one.

---

## 7. Recommendation — handed over, not switched

**Do not turn `transvascular_exchange` on as part of this work, and I have not.**

The evidence splits cleanly:

**For turning it on.** `C ≤ Cp` holds by construction. The entry elasticity
falls from exactly 1 to `k/(PS + k)`, so the underived 1 s⁻¹ stops propagating
without bound into every absolute per-cell call. §4.1 shows that bound is not
academic: five decades of an unmeasured scale currently take the same tissue
from 23,811 unreached cells — every viable tumour cell on the section — to 140.

**Against turning it on now.**

1. It **inverts a shipped conclusion** (§4, daunorubicin), so it is a
   re-validation exercise and not a flag flip. Every absolute number in
   `results/` is on the other law.
2. It does **not** improve accuracy against the only out-of-sample spatial
   anchor the project has (§3), and it moves it 0.15 % the wrong way.
3. It does **not** introduce a measured entry parameter (§5). Turning it on
   with `PS = φ · bbb_penetrance · 1 s⁻¹` makes the *form* of the law right
   while the *value* stays a mislabelled partition coefficient times an
   underived unit. `results/transvascular_exchange_revalidation.md` §6.2
   already established that the value is wrong in both directions across drugs
   and that no single rescale fixes it.
4. The prerequisite is a per-compound wall permeability `P [µm/s]` and a tissue
   `S/V [1/µm]`. **Neither exists in this repository and I did not find a
   source for either.** That is the early answer the brief asked for: the split
   *can* be made — the code is here and it works — but it cannot be made
   *honest* without a parameter nobody has measured.

**What I would do instead, and it is §1's finding, not this one's.** The
accuracy gap is a single length that is 1.8–2.4× too long. `λ` cannot say
whether that sits in `D_eff` or in `k_loss`, and that is now the highest-value
open question in the transport model. It is a much smaller target than "the map
is 14× too shallow", and unlike the entry term it is reachable by the anchors
the project already has.

---

## 8. Provenance

* Driver: `experiments/permeability_split.py`, stages `entanglement`,
  `lambda`, `tissue`, `scope`. Data: `results/permeability_split_data.json`.
* Slab geometry, compound compilation and the Kyle Table 1 values are taken
  **unchanged** from `experiments/isomer_penetration_map.py` so §1's 2.9977 is
  directly comparable with `results/kyle2004_penetration.md` §6.1's 3.00.
* Real-tissue arm: `validation.vessel_profile.build_arm_microenvironment(
  dx_um=12.0, repair=1)` — the `dx12_repaired` arm of
  `results/kyle2004_penetration.md` §7.
* Per-cell labels: `tmesim.sim.mechanism.attribute`, which reconciles
  bit-exactly with the engine's own kill fraction or the run fails.
* Kyle CJ *et al.*, Clin Cancer Res 2004 — quoted through
  `experiments/isomer_penetration_map.py`'s `KYLE_TABLE1`, not re-retrieved.
  Measured PS anchors quoted through `tmesim/validation/ps_scale.py`
  (PMID 524251, PMID 39690422), not re-retrieved. **No citation is introduced
  by this file.**
* `calibrate_alpha()` asserted live at `22.223477711172563` in every stage.
  `PDEConfig().transvascular_exchange`, `transvascular.SPLIT_DEFAULT` and
  `drugs.PS_SCALE_S_INV` asserted at their shipped values (`False`, `False`,
  `1.0`) at the top and bottom of every run.
