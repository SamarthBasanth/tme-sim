# The oxygen route's downstream effect grew 12x — it is a normalisation change, not a physics improvement, and it manufactured `F_state`

> **One-line answer.** The growth is real arithmetic on the shipped configuration and it is caused by commit
> `ceb96c0` (discrete vessels promoted to default), confirmed by running the benchmark's own
> kill-delta computation on both sides of the switch. But it is **not** evidence that the oxygen
> route matters more than we thought. **98.6 % of it is the tissue-integrated oxygen source
> falling 2.78x**, which drags the mean tumour O2 from 0.5945 to 0.3486 — across the
> `QUIESCENT` band ceiling at O2 = 0.4 — and flips `F_state` from **exactly 1.0000** to
> **0.5502**. The physically better vessel geometry, which is the actual argument for the
> promotion, contributes **1.4 %**. Neither 0.026 nor 0.3252 is a measurement of anything: both
> are readings of an uncalibrated level against an absolute banding threshold.
>
> **The single most consequential consequence, verified by re-running the unmodified generator
> rather than projected:** `results/kill_calibration.md`'s headline **"3-D:2-D ratio test PASSES
> for all 3 drugs"** is **2/3 under the immediately preceding configuration**, and the drug that
> fails is temozolomide — the one that report calls "the matched anchor".

**Status: diagnosis. Nothing was tuned, nothing was reverted, and no file outside the three
this report owns was edited.** `ceb96c0` is left in place; `DX_UM`, `alpha`, `OxygenConfig`,
`PDConfig` and every frozen `DRUG_LIBRARY` snapshot are untouched. The only intervention
anywhere is a harness-level rebinding of the name `real_pipeline.vascular_density_from_cells`
to a wrapper that forwards `mode='kde'` — a supported public argument of the unmodified
library function, and the exact historical default `ceb96c0` replaced.

* Generator: `experiments/oxygen_level_shift.py`
* Raw: `results/sweep/oxlevel_fields.json`, `oxlevel_downstream.json`, `oxlevel_invariants.json`,
  `oxlevel_factorial.json`, `oxlevel_cliff.json`, `oxlevel_fstate.json`, `oxlevel_killcal.json`
* Section: `data_cache/xenium_gbm`, 40,887 cells (23,812 tumour), dx = 12 um, grid 424x675,
  alpha = 148.6082, seed 0. Drugs: temozolomide, cisplatin, doxorubicin.

---

## 1. What changed, and the experiment that confirms it

`results/oxygen_demotion.md` §6 records the superseded downstream arm as kill deltas of
**-0.0039 / -0.0262 / -0.0160** (temozolomide / cisplatin / doxorubicin), max |delta| **0.0262**,
and states in the same sentence that those numbers "were computed pre-`4ad221f` and
pre-`ceb96c0`". Regenerating after `a9bb703` gives max |delta| **0.3252**. That is 12.41x.

`ceb96c0` is the cause. Measured by running `validation.oxygen_benchmark.downstream_effect()`
— unmodified — on both sides of the vascular-mode switch, same section, same dx, same drugs,
same seed:

| arm | temozolomide | cisplatin | doxorubicin | **max abs delta kill** | mean O2 (PDE) | mean O2 (surrogate) |
|---|---:|---:|---:|---:|---:|---:|
| **`mode='discrete'`** (production, post-`ceb96c0`) | 0.0184 -> 0.0315 (**+0.0131**, +71.4 %) | 0.4362 -> 0.7614 (**+0.3252**, +74.6 %) | 0.2595 -> 0.4593 (**+0.1998**, +77.0 %) | **0.3252** | 0.3520 | 0.6216 |
| **`mode='kde'`** (historical, pre-`ceb96c0`) | 0.1032 -> 0.1002 (**-0.0030**, -2.9 %) | 0.9289 -> 0.8888 (**-0.0401**, -4.3 %) | 0.7466 -> 0.7206 (**-0.0260**, -3.5 %) | **0.0401** | 0.5983 | 0.6216 |

(The mean-O2 columns are `reconstruction_report`'s `mean_oxygen`, which is **tissue-masked**.
§2's table reports whole-grid and at-tumour-cell means of the same fields, which is why 0.3520,
0.35389 and 0.3486 all appear in this document — three different averages of one field, not
three measurements.)

**Ratio 8.11x**, and the KDE arm reproduces the superseded numbers' **sign, ordering and order
of magnitude** (-0.0030 / -0.0401 / -0.0260 against -0.0039 / -0.0262 / -0.0160). The residual
0.0262 -> 0.0401 (1.53x) is not attributable to `ceb96c0` and is not isolated here; the
superseded document names `4ad221f` (the classifier-arity fix, which changes which cells are
endothelial and therefore both fields) as the other intervening change.

**The flagging note understated it in one respect that matters: the sign reversed.** Under the
historical field the surrogate killed *less* than the PDE on all three drugs; under the
production field it kills *more* on all three. This is not a magnitude change, it is a
direction change.

### 1.1 Why it can only be the PDE arm

Structural, then measured. `physics.oxygen.solve_oxygen` reads
`micro.fields.vascular_density` (`source = supply * vasc`), which is the field `ceb96c0`
replaced. `physics.oxygen_proxy.oxygen_from_vessel_distance` defaults to
`source='endothelial'` and reads **only** `micro.cells` positions — it never touches
`vascular_density`. So the surrogate arm must be invariant to the vascular mode, and it is:
the surrogate field is **bit-identical** (SHA-256 of the raw buffer) under `discrete`, `kde`
and the dose-matched variant. Every number in the table above moved because arm A moved.

---

## 2. The mechanism, end to end

| quantity (dx = 12 um, over tissue) | KDE (historical) | **discrete (production)** |
|---|---:|---:|
| tissue-integrated vascular source | 14464.0 | **5199.0** (2.782x smaller) |
| support — fraction of tissue with nonzero source | 99.780 % | **2.683 %** |
| mean vascular density | 0.07465 | **0.02683** |
| mean O2, whole grid | 0.6002 | **0.35389** |
| **mean O2, at tumour cells** | **0.5945** | **0.3486** |
| max O2 anywhere on the grid | 0.6978 | **0.4296** |
| tumour bands HYPOXIC / QUIESCENT / VIABLE | 0.0000 / **0.0000** / **1.0000** | 0.0000 / **0.8996** / 0.1004 |
| **`F_state`** | **1.0000** | **0.5502** |

`cells.agents.update_cell_states` bands on **absolute** thresholds: `HYPOXIC` below
`t = 0.2`, `QUIESCENT` on `[t, 2t) = [0.2, 0.4)`, `VIABLE` above. The historical field put the
whole tumour above 0.4; the production field puts 89.96 % of it inside `[0.2, 0.4)`. Under the
production field the **maximum** oxygen value anywhere on the grid is 0.4296 — the entire
section is at most marginally viable.

`kill_fraction` is `mean_i[p_i * sens_i * resist_i]`, linear in the per-cell resistance
multiplier, so the kill ratio between two oxygen routes is close to their `F_state` ratio.
That is what the downstream table is measuring:

| arm | `F_state` (PDE) | `F_state` (surrogate) | ratio | measured relative kill change |
|---|---:|---:|---:|---|
| discrete | 0.5502 | 0.9550 | **1.7357** | +71.4 % / +74.6 % / +77.0 % |
| kde | **1.0000** | 0.9550 | **0.9550** | -2.9 % / -4.3 % / -3.5 % |

The surrogate's own bands (H 0.002562 / Q 0.086679 / V 0.910759, `F_state` 0.9550) are
identical under both vascular modes and reproduce `kill_calibration.CELL_STATE_ROUTES`'
published 0.9550 exactly.

**The whole 12x is this one number.** Under the historical field the oxygen banding was doing
*nothing at all* — `F_state` was exactly 1.0000, every tumour cell VIABLE — so the PDE-vs-surrogate
disagreement had almost nothing to disagree about.

---

## 3. Is the new number more correct, or less?

Split the question the way the field itself splits.

### 3.1 The support: discrete is unambiguously the entitled representation

Oxygen enters tissue at vessel walls. The KDE supplies oxygen over **99.780 %** of the
tissue — it is a 40 um Gaussian smear of endothelial cell positions, max-normalised to the
densest spot *in this particular section*. As an oxygen source that is not defensible, and
`real_fields.VASCULAR_MODES`' own documentation calls the 99.7 % support "the structural
confound … and the reason the default moved away from it". The discrete field's 2.683 %
support is at least the right kind of object. **On support, `ceb96c0` moved in the right
direction.**

### 3.2 The level: that is where the 12x actually came from, and it is not a physics result

The controlled comparison — hold the total source fixed, change only *where* it sits. The
KDE field scaled **down** to the discrete field's tissue total (down, because `solve_oxygen`
clips — see §3.4):

| field | source total | support | mean O2 (tumour) | H / Q / V | `F_state` |
|---|---:|---:|---:|---|---:|
| KDE, as shipped historically | 14464.0 | 99.780 % | 0.5945 | 0.0000 / 0.0000 / 1.0000 | **1.0000** |
| **KDE, rescaled to the discrete total** | **5199.0** | **99.780 %** | **0.3495** | 0.0000 / 0.9117 / 0.0883 | **0.5441** |
| discrete, as shipped now | 5199.0 | 2.683 % | 0.3486 | 0.0000 / 0.8996 / 0.1004 | **0.5502** |

At matched total source the two fields agree on mean tumour O2 to **0.26 %** and on `F_state`
to **1.1 %**, despite a **37x** difference in support. Decomposing the full move:

* total source (2.78x): **0.4559 of the 0.4498 `F_state` change — 101 %**
* geometry (99.78 % -> 2.683 % support, at matched total): **0.0061 — 1.4 %**

(The two exceed 100 % because they push in opposite directions slightly; the honest reading is
that **the dose term is essentially all of it and the geometry term is noise beside it**.)

This is the same conclusion `ceb96c0` reached for *delivery* — "That is ~97 % dose, not
geometry" — now shown to hold for the *oxygen* channel too. `ceb96c0` flagged the
unit change loudly for every absolute exposure number and did not flag it for oxygen. It
should have: the oxygen solve's source is `supply * vascular_density` and inherits the
normalisation change linearly.

### 3.3 So: neither number is more correct, and that is the finding

`OxygenConfig.supply = 0.05` is a relative rate multiplying a field that **both** conventions
document as relative — KDE 1.0 means "the densest endothelial-KDE spot in this section",
discrete 1.0 means "a perfused vessel wall leaking at the reference rate". `ceb96c0` changed
the referent of the unit and kept the coefficient. Nothing was re-derived, because nothing was
noticed. Meanwhile the cell-state bands are tested against **absolute** cut-offs at 0.2 and 0.4.

The consequence is that the oxygen route's downstream effect is **unidentified**, not
measured. Sweeping a scalar multiplier on the production source (equivalent to sweeping
`supply`; `results/sweep/oxlevel_cliff.json`):

| source x | mean O2 (tumour) | HYPOXIC | QUIESCENT | VIABLE | `F_state` |
|---:|---:|---:|---:|---:|---:|
| 0.50 | 0.21263 | **0.2676** | 0.7324 | 0.0000 | 0.4599 |
| 0.75 | 0.28734 | 0.0000 | **1.0000** | 0.0000 | 0.5000 |
| **1.00 (shipped)** | **0.34863** | 0.0000 | 0.8996 | 0.1004 | **0.5502** |
| (KDE convention) | 0.5945 | 0.0000 | 0.0000 | **1.0000** | **1.0000** |

A factor of 2 on an uncalibrated coefficient walks `F_state` across essentially its entire
attainable range (0.46 to 1.00; `PDConfig`'s multipliers bound it to [0.35, 1.0]). This is
the source-side twin of `results/oxygen_sensitivity_band.md`, which swept the **sink**
(`consumption`) over 2,500x and found the two literature-defensible sub-ranges predict
opposite headline answers. **That study reports only the HYPOXIC fraction, which is 0.0000 at
the shipping default and 0.0000 under the KDE — so the one band that did not move is the only
one it reports, and it is structurally blind to the band this entire 12x lives in.** The
source axis it never swept is the one `ceb96c0` moved.

### 3.4 A hazard worth stating loudly (not a code change)

`physics/oxygen.py:470` does `vasc = np.clip(vasc, 0.0, 1.0)` before forming the source.
`vascular_density` is documented as a [0,1] fraction so the clip is defensible in itself, but
it has a consequence nobody has written down: **`discrete_vascular_field`'s `match_total_to=`
control is a silent no-op for the oxygen solve whenever it scales the field up.** Measured —
the dose-matched discrete field (scale factor 2.7821, so every vessel pixel goes to 2.7821 and
clips back to 1.0) produces an oxygen field **bit-identical** to the unmatched one
(SHA `2fb1bae6cc4a7f36`, mean 0.3538863603 to ten digits, both). The same control applied to
*transport* works exactly as advertised (dose-matched discrete transport gives cisplatin kill
0.9284 against the KDE's 0.9289). Anyone who uses `match_total_to` to make a dose-matched
claim about oxygen will get a number that looks like a result and is an artefact of the clip.
This is why §3.2's control scales the KDE **down** instead.

### 3.5 The oxygen solve inherits "Wrong #2" in full

`results/vascular_segmentation_audit.md` established that the production `vascular_density`
is bit-identical across linking distances because `discrete_vascular_field` never reads
cluster labels, and that the delivery source is still the per-endothelial-cell 7.5 um disc
that `vessel_geometry.py` calls "Wrong #2". **The oxygen solve consumes that same array**, and
inherits both properties exactly. Measured on this section:

```
link= 10.0 um  n_vessel_objects=4149  vasc SHA 1a30cb3950c0c6ee  oxygen SHA 2fb1bae6cc4a7f36
link= 20.0 um  n_vessel_objects=3091  vasc SHA 1a30cb3950c0c6ee  oxygen SHA 2fb1bae6cc4a7f36
link= 50.0 um  n_vessel_objects=1765  vasc SHA 1a30cb3950c0c6ee  oxygen SHA 2fb1bae6cc4a7f36
link= 80.0 um  n_vessel_objects= 825  vasc SHA 1a30cb3950c0c6ee  oxygen SHA 2fb1bae6cc4a7f36
```

Collapsing 4,149 vessel objects into 825 leaves the solved oxygen field **bit-identical**. So
"the oxygen PDE now runs on discrete vessels" is not true in any sense the solve can detect:
it runs on a union of 7.5 um discs around endothelial *cells*, and the vessel segmentation is
invisible to it. The support area that results — **2.683 %** of tissue — sits above the
2.1 +/- 0.3 % the module's own `PUBLISHED` block cites, and it is set by an unvalidated
per-cell painting rule rather than by any measured caliber.

**Net verdict.** The new number is a better *representation* resting on the same *uncalibrated
level*, and the level is what moved. The correct thing to say about the oxygen route's effect
on absolute kill is not 0.026 and not 0.3252 but **"somewhere in [0.026, 0.38], set by a
coefficient nobody has fitted"**.

---

## 4. The two invariants, verified independently

Both claims in the flagging note hold. Re-derived from `delivery_metrics.compute_all` on the
engine's own transport result rather than read off the benchmark's summary
(`results/sweep/oxlevel_invariants.json`):

**Coverage is exactly invariant — `==`, not `isclose`.** Under both vascular modes, for all
three drugs, all four coverage statistics (`coverage_0.5x`, `coverage_1x`, `coverage_2x`,
`median_E_over_ic50`) compare **exactly equal** between the PDE and surrogate arms. Sample:
cisplatin `0.9954644716949437 == 0.9954644716949437` under discrete,
`0.998698135393919 == 0.998698135393919` under KDE. The per-drug AUC buffers hash identically
across the two oxygen arms by construction — the surrogate arm re-bands cells on the PDE arm's
own `transport` object.

**Cross-drug ranking is unchanged.** `[cisplatin, doxorubicin, temozolomide]` by kill under the
PDE, by kill under the surrogate, and by coverage — under **both** vascular modes. Four
rankings, one order.

**And the two facts together are more informative than either.** Coverage cannot move because
the oxygen field is not an input to transport at all. Ranking cannot move because the oxygen
field enters as a per-cell resistance multiplier that is **drug-independent** — the same
`F_state` multiplies all three drugs (the measured relative kill changes span only
71.4-77.0 %, a 5.6 pp spread across drugs against a 74 % common shift). So ranking invariance
here is **structural, not empirical luck**, and "ranking defensible, level not" survives — but
the exact same structure is what makes the level maximally exposed: a near-common multiplier
is invisible to a ranking and goes straight into every absolute number.

---

## 5. Affected published claims

**Not edited — listed.** The claim class is: *any absolute kill, kill-derived ratio, or
`F_state` figure computed on the GBM brain section through the oxygen-banded cell-state route.*
The size of that quantity's route uncertainty has grown from +/-0.026 to +/-0.33.

### 5.1 The one that changes a verdict — verified by re-run, not projected

`results/kill_calibration.md` headlines **"The 3-D:2-D ratio test PASSES for all 3 drugs"**,
with the temozolomide prediction "essentially on the literature median" and the statement
"Nothing in that prediction was tuned". Running the **unmodified**
`validation.kill_calibration.run()` + `analyse()` under `mode='kde'` — nothing else changed:

| drug | predicted ratio, KDE | literature band | verdict, KDE | published (discrete) | published verdict |
|---|---:|---|---|---:|---|
| **temozolomide** | **0.46784** | 0.15899 - 0.46221 | **FAIL** | 0.2626 | PASS |
| cisplatin | 0.49795 | 0.20864 - 1.50035 | PASS | 0.2791 | PASS |
| doxorubicin | 0.36883 | 0.05746 - 1.00000 | PASS | 0.2077 | PASS |

**2/3, and the drug that fails is the one with the matched anchor.** The published `3/3 PASS`
is a consequence of `ceb96c0`. This does not make the published result wrong — the discrete
source is the better representation (§3.1) and the run is honest — but the report's framing
that nothing in the prediction was tuned needs the rider that the prediction *does* depend on
an uncalibrated normalisation convention that was changed for unrelated reasons on
2026-07-28 and never re-derived.

The report is not unaware of the sensitivity: §5.1 sweeps the forced quiescent fraction and
correctly bounds `F_state` to [0.5, 1.0]. **But it propagates that sweep only into absolute
kill (2.00x, and concludes cell state is "non-dominant"), never into §3's verdict** — and the
`F_state = 1.0000` endpoint it treats as a hypothetical worst case is precisely the
configuration this project shipped until 2026-07-28.

### 5.2 Files quoting a figure this moves

| file | what in it | note |
|---|---|---|
| `results/kill_calibration.md` | the 3/3 PASS, `F_state` = 0.5502 in every decomposition row, implied IC50 folds | §5.1; verdict-changing |
| `results/cell_state_impact.md` | kills 0.0184 / 0.4362 / 0.2595, ceilings 0.5502-0.9561, the "1.74x route swing" headline | **already documents the ~3.4x supply change at §-line 83** — the most honest existing treatment; it does not record that the KDE arm gives `F_state` **exactly** 1.0000 |
| `results/grid_convergence_gbm.md` | `F_state` = 0.5502 at dx = 12, the ladder, the 3/3 PASS | §6 — the non-monotonicity it diagnoses is created by `ceb96c0` |
| `results/tail90_falsification.md` | `F_state` carrying 25-51 % of var(ln ratio) on GBM; quiescent 0.8996-0.9918 | §6 |
| `results/p99_convention.md` | `F_state` 0.5305-0.5502 within dx = 10-12 on GBM | §6 |
| `results/dimensionality.md` | propagates `F_penetration x F_state` = 0.2626 | inherits §5.1 |
| `results/dose_schedule.md`, `results/regimen.md`, `results/flagship_markers.md`, `results/delivery_source_geometry.md` | the 0.5502 temozolomide kill ceiling | ceiling is definitionally `sens x resist`, so it *is* `F_state` |
| `results/real_drug_panel.md`, `results/delivery_metrics.md`, `results/patient_report_xenium.md`, `results/immune_sink_impact.md`, `results/heldout_patient.md`, `results/heldout_patient_validation.md`, `results/mechanism_attribution.md`, `results/alpha_identifiability.md`, `results/breast_cancer_flagship.md` | absolute GBM kills 0.0184 / 0.4362 / 0.2595 | levels only; every **coverage** and **ranking** claim in these is untouched (§4) |
| `results/oxygen_demotion.md` §6 | the -0.0039 / -0.0262 / -0.0160 deltas | already marked stale by its own text; the sign, not just the digits, has flipped |
| `docs/DATA_ACCESS_PLAN.md`, `docs/outreach/02_erasmus_lamfers.md` | quote `F_state` 0.9550 | route-comparison figures |
| `app/server/report.py`, `results/viewer_data.json`, `app/runs/*/result.json` | Studio-emitted absolute kills | product surface; `app/runs/` is gitignored run artefacts |

**Not affected, and worth saying so explicitly:** every therapeutic-coverage number anywhere in
`results/`, and every cross-drug ranking. §4 verifies both directly rather than by assertion.

**Also worth recording:** `results/oxygen_benchmark.md` does not exist. `a9bb703` moved the
generator's default output there and describes running it, but the regenerated document was
never committed, so the 0.3252 has until now existed only in a commit message.

---

## 6. The `F_state` connection — same story, third direction, and it is causal

`results/grid_convergence_gbm.md` §5 identifies non-monotone `F_state` as the mechanism behind
the GBM ladder's non-monotonicity, measuring quiescent 0.8996 at dx = 12 and 0.9918 at dx = 30.
`results/tail90_falsification.md` §6 measures `F_state` carrying 25-51 % of var(ln ratio) on
GBM against 0-11 % on breast, and routes the problem to "whoever owns the oxygen/cell-state
path". Running that same ladder under both vascular fields
(`results/sweep/oxlevel_fstate.json`):

| dx (um) | discrete: mean O2 | discrete: quiescent | **discrete: `F_state`** | kde: mean O2 | kde: quiescent | **kde: `F_state`** | surrogate: `F_state` |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 30 | 0.3399 | 0.9918 | 0.5041 | 0.6056 | 0.0000 | **1.0000** | 0.9560 |
| 20 | 0.3489 | 0.9196 | 0.5402 | 0.5985 | 0.0000 | **1.0000** | 0.9560 |
| 16 | 0.3528 | 0.9334 | 0.5333 | 0.6027 | 0.0000 | **1.0000** | 0.9551 |
| **12 (published)** | 0.3539 | **0.8996** | **0.5502** | 0.6002 | 0.0000 | **1.0000** | 0.9550 |
| 10 | 0.3519 | 0.9390 | 0.5305 | 0.5984 | 0.0000 | **1.0000** | 0.9549 |
| 8 | 0.3536 | 0.9024 | 0.5488 | 0.5988 | 0.0000 | **1.0000** | 0.9554 |

The discrete column reproduces `grid_convergence_gbm.md` §5 exactly (0.9918 at 30, 0.8996 at
12, 0.9390 at 10, 0.9024 at 8), which is the cross-check that this is the same quantity.

**The KDE column is exactly 1.0000 at every rung. Span zero. The surrogate column spans 0.11 %.**

**So `F_state`'s grid dependence is not a property of the section, the solver, or the mesh — it
is created by `ceb96c0`, and it did not exist under either of the other two oxygen routes.**
The reason is now mechanical: the discrete field parks mean tumour O2 at 0.349, a hair inside
the `[0.2, 0.4)` `QUIESCENT` band with 90-99 % of tumour cells in one band, so a grid change
that nudges the field by a few percent reshuffles thousands of cells across the 0.4 edge. The
KDE parks it at 0.595, clear of the edge, so nothing crosses at any dx.

This is the same finding as §3 from a third direction: **an uncalibrated level tested against
an absolute threshold**. It gives `tail90_falsification.md`'s referral a concrete target — the
question for the oxygen workstream is not "why is `F_state` non-monotone" but "why is the mean
tumour O2 sitting on a band edge, and what fixes the level".

**One thing this does NOT establish**, and the tempting inference is wrong:
`results/p99_convention.md` reports breast `F_state` at 0.5049-0.5152, so the breast section
sits **inside the quiescent band too**. Landing in the band is therefore not by itself the
GBM/breast discriminator that `tail90_falsification.md` §6 measures at 10x; what differs is the
*lability* of the labelling within the band. Testing that on breast was not done here.

---

## 7. Limitations

1. **One section, one dx for the causal test.** All downstream and factorial numbers are the
   40,887-cell GBM brain section at dx = 12 um. The `F_state` ladder covers dx = 8-30 um.
2. **`4ad221f` is not isolated.** The 1.53x residual between the superseded 0.0262 and this
   report's KDE-arm 0.0401 is attributed by difference, not by experiment. Isolating it needs a
   checkout, which this task was not permitted to do.
3. **`kill_calibration` was re-run under the KDE only.** The discrete arm is read from
   `results/kill_calibration.md` rather than re-run; the cross-check that this is the same
   pipeline is that the KDE run reproduces temozolomide kill 0.1032252646, matching this
   report's own §1 downstream table to seven digits.
4. **The dose/geometry split is measured at one matched total.** §3.2 matches the KDE down to
   the discrete total; the reverse direction is unavailable because of the clip (§3.4).
5. **Nothing here calibrates a level.** This report says the level is unidentified and shows
   how far it moves. It does not propose a value for `supply`, and deliberately did not touch
   it, `consumption`, `hypoxia_threshold` or the `PDConfig` multipliers.
6. **No defect requiring a code change was found.** `tmesim/` is untouched by this work. The
   clip in §3.4 is a documented-behaviour hazard rather than a bug, and the "Wrong #2"
   inheritance in §3.5 belongs to the vessel-geometry workstream that is already live in
   `tmesim/reconstruct/vessel_geometry.py`. The `subgrid` geometry that agent owns is **not**
   on the default path, so none of it affects the numbers here.

---

## 8. Reproduce

```
./.venv/bin/python experiments/oxygen_level_shift.py fields downstream invariants factorial cliff fstate killcal
```

Each stage checkpoints to `results/sweep/oxlevel_<stage>.json` as it goes. Total wall clock on
this machine: about 16 minutes, dominated by `downstream` (2 x 118 s), `invariants` (235 s) and
`factorial` (9 engine runs, 395 s).

`./.venv/bin/python -m pytest -q` -> **754 passed**, 8 warnings, 297 s — the stated baseline.
Nothing under `tmesim/` was modified by this work.
