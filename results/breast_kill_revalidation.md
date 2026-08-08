# Breast 3-D:2-D kill calibration, re-validated under corrected transvascular physics

**One physics change, ten drugs, zero verdict flips — and one result worth more than the survival itself.**

This is the full record of re-running `results/kill_calibration_breast.md` under the
corrected vascular delivery law. It lifts the conditional annotation that file carried
(its §6 limitation 0), which said the ratio verdicts should survive a per-drug delivery
scale error but that *"'largely' is not 'exactly', and no one has quantified the
residual."* The residual is now quantified.

It is a companion to `results/transvascular_exchange_revalidation.md`, the consolidated
record of all eleven re-validations, and does not restate or contradict it. Where the
two touch — the ranking claim in particular — this file defers to it.

---

## 1. The change, and what was re-run

`tmesim/treatment/drugs.py::delivery_source` was a ONE-WAY leak:

    S = PS * Cp                     (no path back from interstitium to blood)

whose steady state `PS*Cp / (clearance + binding)` is unbounded above the plasma level.
The physically correct form is gradient-driven exchange:

    S = PS * (Cp - C)  =  PS*Cp  -  PS*C

the same source minus a first-order loss at the same rate, so it introduces no new
parameter. Steady state becomes `PS*Cp/(PS + clearance + binding) <= Cp`, bounded by
construction. It is behind `PDEConfig(transvascular_exchange=True)`, default off.

**Re-run:** all ten drugs of the breast panel, on the same real 10x Xenium breast
section (Janesick *et al.* 2023, PMID 38114474, Rep1, dx = 12 um, 167,780 cells /
89,785 tumour), same `alpha = 148.6082`, same cell-state labelling, both arms, and the
published 3-D:2-D decomposition recomputed from each.

### 1.1 The runs genuinely differ

`results/transvascular_exchange_revalidation.md` §6 documents a real verification trap:
several modules call `simulate_drug_transport` directly, for which the `PDEConfig` flag
is a **silent no-op**, so an unguarded re-validation reports "nothing changed" when in
fact nothing was tested. Five of the eleven re-validations would have fallen into it.

This one was guarded explicitly, and passed all three checks:

| guard | result |
|---|---|
| PS field non-zero | asserted, per drug |
| ON vs OFF exposure fields not bit-identical | asserted, per drug |
| per-cell exposure change, median across tumour cells | **2.0% - 89.2%** across the panel |
| per-cell exposure change, maximum | **5.5% - 94.1%** |

The arms differ by up to 94% cellwise. Nothing below is a no-op reported as a survival.

### 1.2 The OFF arm reproduces the published report

Every drug reproduces `results/kill_calibration_breast.md` **to 4 decimal places** with
exchange off. `F_state = 0.5063` for all ten (identical by construction — the cell-state
labelling comes from the oxygen solve, which exchange never touches), 167,780 cells /
89,785 tumour, `alpha = 148.6082`. The published verdict structure is confirmed: **6
formally testable drugs, 5 PASS + docetaxel FAIL, 4 untestable** (cyclophosphamide,
tamoxifen, vinorelbine, palbociclib).

## 2. The result — per drug, off to on

Ratios are at the reporting concentration (1xIC50), the same convention §3 of the
published report uses. Implied IC50 folds are recovered from the ratio by the report's
own Hill inversion (`r/2 = 1/(1 + f^1.20)`, i.e. `f = (2/r - 1)^(1/1.2)`); that formula
reproduces every published OFF fold exactly, so the ON folds are derived, not measured
independently.

| drug | ratio OFF | ratio ON | change | fold OFF | fold ON | band (folds) | verdict OFF | verdict ON |
|---|---:|---:|---:|---:|---:|---|---|---|
| doxorubicin | 0.2301 | 0.2383 | 1.036x | 5.47x | 5.30x | 1.0-18.8 (drug-matched) | **PASS** | **PASS** |
| paclitaxel | 0.1503 | 0.1522 | 1.013x | 8.10x | 8.01x | 0.4-25.0 (broad) | **PASS** | **PASS** |
| cyclophosphamide | 0.2838 | 0.4119 | **1.451x** | 4.48x | 3.08x | none drug-matched | NOT TESTABLE | NOT TESTABLE |
| tamoxifen | 0.1719 | 0.2154 | 1.253x | 7.17x | 5.82x | none drug-matched | NOT TESTABLE | NOT TESTABLE |
| docetaxel | 0.1484 | 0.1501 | 1.012x | 8.19x | 8.11x | 8.3-58.0 (breast-specific) | **FAIL** | **FAIL** |
| epirubicin | 0.2301 | 0.2383 | 1.036x | 5.47x | 5.30x | 1.0-76.0 (breast-specific) | **PASS** | **PASS** |
| 5-fluorouracil | 0.3379 | 0.4613 | 1.365x | 3.77x | 2.73x | 0.4-25.0 (broad) | **PASS** | **PASS** |
| vinorelbine | 0.1577 | 0.1612 | 1.022x | 7.76x | 7.60x | none drug-matched | NOT TESTABLE | NOT TESTABLE |
| cisplatin | 0.3085 | 0.4409 | **1.429x** | 4.13x | 2.87x | 0.4-6.0 (drug-matched) | **PASS** | **PASS** |
| palbociclib | 0.2193 | 0.2997 | 1.367x | 5.73x | 4.25x | none drug-matched | NOT TESTABLE | NOT TESTABLE |

**Zero verdict flips.** Every band membership is unchanged. The four drugs with no
drug-matched band remain untestable — exchange cannot create literature, and their
context-only comparison against the broad band (ratio range 0.0412-1.5004) still lands
inside it in both arms.

Panel-wide, the scored-ratio spread barely moves: **2.28x OFF -> 3.07x ON**
(0.3379/0.1484 -> 0.4613/0.1501). Compare that with the delivery spread, §4.

## 3. Docetaxel does not flip — and moves the wrong way to be rescued

Docetaxel is the only FAIL, and it is the obvious candidate for a physics change to
rescue. It is not rescued, and the reason is directional, not marginal.

Its breast-specific band is 8.3-58.0x, i.e. a ratio band of **0.0152-0.1463**. Its
failure mode is that the model predicts a **SMALLER** 3-D penalty than measured: the
predicted ratio sits **above** the band's upper edge (equivalently, its implied fold
8.19x sits below the measured minimum 8.3x).

| docetaxel | OFF | ON |
|---|---:|---:|
| ratio | 0.1484 | 0.1501 |
| upper band edge (fold 8.3x) | 0.1463 | 0.1463 |
| margin outside | −0.0022 | −0.0039 |
| % outside the edge | 1.5% | 2.6% |

**Every ratio in the panel moves UP under exchange** (§2, all ten). Since docetaxel
fails by being too high, an upward move can only push it further out. The margin roughly
doubles. There is no version of this physics change that rescues docetaxel; a rescue
would require the ratios to move down.

(Recomputing the margin from the 4-dp published ratios gives 0.0021 -> 0.0038; the
values in the table are the run's own full-precision figures. The disagreement is
rounding.)

**Caveat, stated because it cuts against the strength of this file.** Docetaxel's FAIL
is knife-edge in **both** physics. 0.1484 against a band edge of 0.1463 is 1.5% outside;
ON it is 2.6% outside. This is not a drug that fails by a wide margin and would have to
travel far to pass — it is a drug that sits just outside its band in both arms, on a
band built from a single paper's two cell lines (see the published report's limitation
13). The correct reading is "the verdict did not change and the physics pushes it
further from a pass", **not** "docetaxel robustly fails".

## 4. The most fragile PASS moves away from its edge

Cisplatin is the PASS closest to a boundary, and it is also the drug the corrected law
hits hardest elsewhere in the repo — its GBM mechanism verdict inverts entirely
(`results/transvascular_exchange_revalidation.md` §4.2). On the breast 3-D:2-D test it
moves in the safe direction.

Its band is 0.4-6.0x, i.e. a ratio band of **0.2086-1.5004**. The binding edge is the
lower one:

| cisplatin | OFF | ON |
|---|---:|---:|
| ratio | 0.3085 | 0.4409 |
| distance to `ratio_lo` (0.2086) | 0.0998 | **0.2322** |

The margin more than doubles. The single most fragile verdict in the panel becomes less
fragile, not more.

## 5. The headline: cancellation

The verdicts surviving is the smaller result. The larger one is *why* they survive, and
it is measurable.

Absolute delivery on this section changes enormously under the corrected law — by up to
**9.58x** for a single drug, and the cross-drug delivery spread collapses from **970x to
104x**. The scored ratios move by at most **1.45x**, and their spread moves 2.28x ->
3.07x. The error largely divides out.

Quantified per drug, with the surviving fraction measured in log space
(`log(ratio change) / log(delivery change)`, the natural measure for a multiplicative
error) and the cancellation factor its reciprocal:

| drug | delivery change | scored-ratio change | residual surviving | cancellation |
|---|---:|---:|---:|---:|
| 5-fluorouracil | **9.58x** | 1.37x | 13.8% | **7.3x** |
| cisplatin | **8.86x** | 1.43x | 16.4% | **6.1x** |
| cyclophosphamide | **6.35x** | 1.45x | 20.1% | **5.0x** |
| the five low-permeability drugs<br>(doxorubicin, paclitaxel, docetaxel, epirubicin, vinorelbine; `bbb_penetrance` 0.019-0.033) | 1.02-1.07x | 1.01-1.04x | ~48% | ~2.1x |

**The narrow claim this supports, stated exactly.** For precisely the drugs whose
delivery was most wrong, **84-86% of the delivery error cancels in the scored ratio**.
The 3-D:2-D ratio is a within-drug, same-tissue comparison — the same drug, the same
section, the same cell-state labelling, differing only in the geometry being compared —
so a multiplicative error in that drug's delivery appears in both halves and largely
divides out. That is direct, measured evidence that Vivum's actual output, *relative
comparison on one tissue*, is substantially insulated from absolute parameter error.

**Three things this does NOT say, all of which matter more than the claim itself.**

1. **The residual is not zero.** It ranges from 13.8% (5-fluorouracil) to roughly 48%
   for the low-permeability drugs, where the delivery change and the ratio change are
   both small. Cancellation is large where the error is large, which is the useful
   direction, but "insulated" is not "immune".

2. **This does not show ranking is invariant.** That claim was **retracted as
   measured-false** by this same physics change: correcting the delivery law changes the
   top-1 GBM drug (cisplatin -> vincristine), because it re-weights the panel by
   permeability. See `results/transvascular_exchange_revalidation.md` §4.1 — top-1
   changes, top-3 changes, and two of the four swapped pairs cross gaps of 40.5% and
   35.6% of the leader's kill against a Bernoulli sd of ~0.003. On breast, top-1
   survives but the top-3 set changes and 7 of 45 pairs are discordant.

3. **The distinction is common-mode versus drug-specific, and it is the whole point.**
   Ratio outputs are robust to **COMMON-MODE** parameter error — an error that
   multiplies both halves of the same comparison. Cross-drug **RANKING** is not robust
   to **DRUG-SPECIFIC** transport parameters — the exchange attenuation
   `(clearance+binding)/(PS + clearance+binding)` is driven by `bbb_penetrance`, which
   spans 47x across the GBM panel, so this perturbation is differential by construction.
   Both statements are true simultaneously, and this file's result is entirely the first
   one. Nothing here licenses the second.

That is also why the cancellation is *unequal* in the table above: the drugs that cancel
best are exactly the high-permeability drugs whose delivery moved most. The cancellation
is a property of the ratio construction, not of the panel.

## 6. Scope — what this re-validates, and what it does not

**This re-validates the breast 3-D:2-D kill calibration under the corrected
transvascular law, and nothing else.** Specifically:

* It covers the §3 ratio verdicts of `results/kill_calibration_breast.md`, at the 1xIC50
  reporting convention, on this one section.
* It does **NOT** re-validate the absolute 3-D kill numbers in that report's §4, or the
  `F_penetration` values in its §2 decomposition. Those are levels, not ratios; a
  common-mode error does not cancel out of a level. They still inherit the delivery
  error and are not rescued by anything here.
* It does **NOT** address the **third known defect**: `bbb_penetrance` is a
  *dimensionless* brain:plasma partition coefficient that the engine uses as a
  *per-second rate*, obtaining `s^-1` via an implicit, undocumented, uncalibrated
  conversion factor of exactly `1 s^-1`. Nothing in the repo derives or justifies that
  scale, and the exchange term added by this branch inherits it — `k_0` multiplies both
  influx and back-flux, so the term is at least self-consistent, but the steady state
  `k_0*PS'*Cp / (k_0*PS' + k)` depends on `k_0`. **A calibrated `k_0` would change every
  exchange-on number in this file.** They are conditional on `k_0 = 1 s^-1`. Full
  diagnosis: `results/transvascular_exchange_revalidation.md` §6.2.
* It does **NOT** address the permeability VALUE, only the law's FORM. PS is still
  `bbb_penetrance`, a blood-**brain**-barrier term, applied to breast tissue. That is a
  second, independent defect and it remains open (§5 of the same file).
* It does **NOT** address **tissue diffusivity**, which is under separate investigation.
* It says nothing about clinical response. The published report's §7 was explicitly not
  a validation and this changes none of that.

## 7. Provenance

Ten drugs, two arms each, on `data_cache/xenium_breast_janesick/rep1`, dx = 12 um,
seed 0, `PDEConfig(transvascular_exchange=True)` versus default. Baseline reproduction
asserted to 4 dp against the published report before any comparison; both arms proven to
differ before any conclusion was drawn. This is one of the eleven parallel
re-validations recorded in `results/transvascular_exchange_revalidation.md` (§3, row
"Breast 3-D:2-D kill calibration"), and the source of that row's summary line
*"5/6 PASS + docetaxel FAIL, zero flips. Docetaxel moves the wrong way to be rescued."*

Related: `results/kill_calibration_breast.md` (the report re-validated),
`results/breast_kp_validation.md` (the external anchor failure that forced the physics
change), `results/kill_calibration.md` (the GBM module this mirrors, re-validated 3/3
PASS).
