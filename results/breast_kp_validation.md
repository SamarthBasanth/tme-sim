# Breast transvascular delivery: the model FAILS its first external anchor

**A falsification, reported as one.** This is the breast flagship's first
external *transport* validation — as opposed to the 3-D:2-D potency-ratio
tests in `results/kill_calibration_breast.md`, which are a different kind of
claim. It does not pass. The engine's predicted interstitial:plasma AUC ratio
misses the measured human value for **both** anchored drugs, in **opposite
directions**, and the root cause is structural rather than a mis-set
parameter.

Produced by `experiments/breast_kp_anchor.py` (driver only; nothing under
`tmesim/` was modified to produce it). Data: `results/breast_kp_anchor_data.json`.

> ### Headline
>
> The model predicts interstitial:plasma AUC ratios spanning **970x** across
> the breast panel (0.0080 to 7.80). The measurement says two chemically
> unrelated drugs both sit at **~0.60**, a spread of **1.02x**.
>
> On the two anchored drugs specifically:
>
> | drug | predicted | measured | verdict |
> |---|---:|---|---|
> | 5-fluorouracil | **7.7956** | 0.61 ± 0.11 (n=10, PMID 9205062) | **FAIL, 12.8x too HIGH** |
> | methotrexate | **0.0571** | 0.60 ± 0.20 (n=9, PMID 9679959) | **FAIL, 10.5x too LOW** |
>
> The model says these two drugs differ in transvascular delivery by
> **136x**. The measurement says they differ by **1.02x**. Because the two
> errors point in opposite directions, **no single calibration constant can
> fix this.** It is not a scale error.

## 1. What was compared, and why it is a like-for-like comparison

The anchor (full provenance in `results/breast_transport_anchors.md` §1) is
in vivo microdialysis, with the probe placed directly in the primary breast
tumour, reporting interstitial drug AUC as a percentage of plasma AUC:

* Müller M *et al.*, *Cancer Res* 1997;57(13):2598-601, **PMID 9205062** —
  5-FU, n=10, **61 ± 11 %** (mean ± SE).
* Müller M *et al.*, *Cancer Res* 1998;58(14):2982-5, **PMID 9679959** —
  methotrexate, n=9, **60 ± 20 %** (mean ± SE, CV = 100 %).

The engine quantity is

    ratio = mean_tumor_auc [uM*s] / (plasma_mean [uM] * total_time [s])

**These are the same physical quantity, not merely analogous.** The transport
PDE's state variable is the *free interstitial* concentration, which is
exactly what a microdialysis probe samples — not total tissue homogenate,
which would include intracellular and protein-bound drug. This is also why
the paclitaxel tissue:plasma figure of 4-70x (PMID 24670687) is NOT used as
an anchor here: that one is total tissue concentration including
tubulin-bound intracellular drug, a different quantity that must not be
pooled with these.

Run conditions: real Janesick Xenium breast section (`rep1`, 167,780 cells),
`dx = 12 um`, 6 h window, `dt = 10 s`, convection on, seed 0,
`sensitivity_mode='uniform'`. IC50 is uncontrolled and irrelevant — AUC is a
pure transport quantity and never enters the Hill kill term. `alpha` is the
value fit once on GBM doxorubicin penetration and reused verbatim; it was
deliberately **not** refit to this anchor.

## 2. The full panel

| drug | MW | bbb_penetrance | predicted ratio | vs 0.60 |
|---|---:|---:|---:|---:|
| docetaxel | 807.9 | 0.0190 | 0.0080 | 0.01x |
| paclitaxel | 853.9 | 0.0191 | 0.0087 | 0.01x |
| vinorelbine | 778.9 | 0.0326 | 0.0168 | 0.03x |
| epirubicin | 543.5 | 0.0206 | 0.0415 | 0.07x |
| doxorubicin | 543.5 | 0.0206 | 0.0415 | 0.07x |
| **methotrexate** | 454.4 | 0.0225 | **0.0571** | **0.10x** |
| tamoxifen | 371.5 | 0.6020 | 0.3163 | 0.53x |
| palbociclib | 447.5 | 0.5788 | 0.8712 | 1.45x |
| cyclophosphamide | 261.1 | 0.9735 | 4.2514 | 7.09x |
| cisplatin | 298.0 | 0.8577 | 6.7200 | 11.20x |
| **5-fluorouracil** | 130.1 | 0.6327 | **7.7956** | **12.99x** |

The ordering is essentially monotone in `bbb_penetrance`. That is the tell.

## 3. Two distinct root causes, both structural

### 3.1 A blood-BRAIN-barrier permeability is being applied to breast tissue

`bbb_penetrance` is a structure-derived estimate of how well a molecule
crosses the **blood-brain barrier**. `tmesim/treatment/drugs.py`'s
`delivery_source` multiplies the vascular source by it:

    S(x, y) = vascular_density(x, y) * bbb_penetrance * plasma_conc

There is **no tissue-type switch anywhere in the compile or transport path** —
verified directly: neither `molecular/compile_drug.py` nor
`physics/diffusion.py` branches on tissue type. Breast tumour vasculature has
no blood-brain barrier. It is, if anything, unusually leaky.

This is why the model puts methotrexate at 0.0571: MTX's `bbb_penetrance` is
0.0225, which is a *correct and well-sourced* fact about the brain (MTX is
dosed intrathecally for CNS disease precisely because it does not cross the
BBB) and an *irrelevant* one about breast. The same mechanism drives the
whole bottom half of the table: docetaxel, paclitaxel and vinorelbine are all
predicted to essentially not reach breast tissue, on the strength of a
barrier that is not present in breast.

Directly contradicted by measurement: paclitaxel is measured at 4-70x plasma
in human breast tumours (PMID 24670687, total tissue), while this model
predicts 0.0087 of plasma.

### 3.2 The delivery source is one-way, so the ratio is unconstrained

The source term is a pure influx proportional to plasma. Nothing carries drug
back from interstitium to blood. `tmesim/treatment/drugs.py:117` says so
itself, in the source:

> *"(A fuller model would use `(plasma_conc - tissue_conc)` as the driving
> gradient; the transport PDE already removes drug via clearance/binding, so
> this linear-leak form keeps the treatment module decoupled from the tissue
> state and is what the contract specifies.)"*

With influx fixed and the only losses being clearance and binding, the
steady-state tissue concentration is set by (influx rate)/(loss rate). That
quantity has **no reason whatsoever to land near the plasma concentration**,
and for a drug with slow clearance it grows without bound relative to it.

That is how 5-FU, cisplatin and cyclophosphamide end up at ratios of 4-8,
i.e. the model asserts that free interstitial drug sits 4-8x *above* plasma.
For passively-delivered small molecules with no active uptake mechanism that
is not physically available — you cannot passively concentrate a solute
several-fold above its source.

**This defect was documented but never quantified.** The parenthetical above
has been in the source all along. What was missing was any external
measurement that made its consequence visible. That is precisely what an
anchor is for.

## 4. Robustness: it is not the plasma-averaging convention

The default `time_resolved=False` collapses the plasma curve to its window
mean. Re-running 5-FU with `time_resolved=True`, so the genuine decaying
plasma curve drives the solve step by step:

| drug | time_resolved | ratio |
|---|---|---:|
| 5-fluorouracil | False | 7.7956 |
| 5-fluorouracil | True | **7.8241** |
| methotrexate | False | 0.0571 |
| methotrexate | True | **0.0572** |

A 0.4 % change on 5-FU and 0.2 % on methotrexate, against discrepancies of
12.8x and 10.5x. The Stage-0 plasma-averaging convention is not the cause.

## 5. What this does and does not invalidate

**Does NOT invalidate:**

* **Cross-drug RANKING on the same tissue**, which is this project's stated
  defensible output. A monotone error shared across drugs preserves order.
  But note the caveat in §6 — the error is monotone in `bbb_penetrance`, not
  constant, so ranking is preserved only to the extent that `bbb_penetrance`
  happens to rank breast delivery correctly, which this result gives no
  reason to believe.
* **The GBM results.** GBM is the tissue `bbb_penetrance` is *for*, and GBM's
  own transport anchor (doxorubicin fit → mitoxantrone predicted
  out-of-sample) is a penetration-depth test that passed. Nothing here
  touches it.
* **The spatial/geometric findings** (e.g. the immune-sink exposure spread),
  which are about relative exposure *within* one tissue at fixed drug and are
  insensitive to a per-drug scale factor.

**DOES invalidate:**

* Any absolute or cross-drug **breast** delivery magnitude from this engine.
* Any reading of breast `mean_tumor_auc`, `therapeutic_coverage_1x` or
  `kill_fraction` as an absolute quantity. These were already flagged
  non-absolute; this supplies a measured number showing how far off.

## 6. What NOT to do about it

**Do not refit `alpha` to this anchor.** `alpha` is a binding/uptake
amplitude fit once to GBM doxorubicin penetration depth. Bending it to hit
0.60 on breast would (a) break the GBM chain it was fit for, and (b) hide a
structural defect behind a tuned constant — the exact failure mode
`results/SESSION_AUDIT.md` exists to prevent. The two errors point in
opposite directions, so no single constant can absorb them anyway.

**Do not simply swap `bbb_penetrance` for a breast number.** That fixes §3.1
and leaves §3.2 untouched: the ratio would still be an unconstrained
influx/loss balance, still free to exceed 1 without limit. Both defects have
to be addressed for the anchor to mean anything.

**Both of those have since been done, and the residual is the alpha-scaled
binding sink `k = clearance + binding` — specifically its 4.78x drug-to-drug
variation, against a measured permeability and a measured anchor that are both
flat to 1.02x.** The chain is `results/breast_tissue_permeability.md`
(measured permeability makes it worse, which localises the defect) →
`results/binding_sink_refuted.md` (the sink is refuted in magnitude and
inverted in direction) → **`results/breast_binding_anchors.md`**, the search
for an independent *measurement* of tissue binding or clearance — the only
thing that can replace `alpha` here without fitting it to this anchor.

The honest fix is a real transvascular exchange term driven by
`(plasma_conc - tissue_conc)` with a tissue-appropriate permeability-surface
-area product, replacing the one-way leak. That is a change to a frozen
contract and a genuine design decision, not a patch. It is written down here
rather than made unilaterally.

## 7. Status

**Breast transport physics is UNVALIDATED and now measurably wrong in
magnitude.** The potency-ratio side (`results/kill_calibration_breast.md`,
5 of 6 testable drugs passing) stands unchanged — it tests a different claim.

The anchor test itself is now a reusable, falsifiable harness. Any future fix
can be scored against it by re-running `experiments/breast_kp_anchor.py`; the
target is both anchored drugs inside their bands and a panel spread far below
970x.
