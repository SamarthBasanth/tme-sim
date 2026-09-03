# The doxorubicin anchor identifies ONE number, and it is not alpha

> # ⚠️ WITHDRAWN SUBJECT — BANNER ADDED 2026-09-03. THE CROSS-COMPOUND ORDERING THIS DOCUMENT RESTS ON NO LONGER EXISTS AS A CLAIM.
>
> **This banner is the only edit made to this file.** Nothing below it has been
> changed, reworded, re-scored or removed. It is preserved byte-identical as the
> record of what was believed and measured at the time, because a record that was
> wrong is part of the record.
>
> **What is below.** A diagnosis of what the doxorubicin anchor actually identifies — that `alpha` and the anchor compound's binding prior are not separately identified, that only their product is, and that `alpha` is a gauge coordinate rather than a measurement. **Those findings are not affected by this banner and are not withdrawn.** What is affected is the closing line of its headline block (line 19): having established that the absolutes do not survive the perturbation it studies, it names a cross-compound ordering as the one output that does hold. That output no longer exists. The sentence predates the withdrawal and is the only place in the file that reads as a live claim.
>
> **What was withdrawn, and when.** On **2026-08-21** Vivum withdrew the ordering
> of one compound against another. The withdrawal is absolute rather than
> conditional: there is no grid spacing, no section and no compound set for which
> this product orders compounds against each other. It is not a confidence
> downgrade and it is not satisfied by a hedge, a stated limit or a stated
> uncertainty — the output was removed, and `compare-drugs` was removed with it.
> `README.md` carries the withdrawal in full.
>
> **The two measurements that withdrew it.**
>
> * [`results/pk_ranking_robustness.md`](pk_ranking_robustness.md) — per-compound
>   plasma exposure across the retained compounds spans **113.1x**, while the
>   published order first breaks at a spread of **1.53x**. The variation that real
>   per-compound dosing supplies is far wider than the margin the order tolerates.
> * [`results/capacity_ranking_robustness.md`](capacity_ranking_robustness.md) —
>   the engine carries **no per-compound tissue-capacity term at all**. Under a
>   *uniform* capacity — every compound given the same value, which is exactly the
>   case the "a constant error cancels out of a ranking" defence covers — the
>   published order still breaks, at **R = 36.3** on the GBM 10-drug panel and
>   **R = 5.28** on the breast 10-drug panel, both below the best available
>   constant of **39**.
>
> **Neither reason is repairable inside this document**, and nothing below is
> re-run or re-scored by this banner.

**A diagnosis. No parameter was changed, nothing was tuned to an anchor, and no code in `tmesim/` was edited.** Everything below is either closed-form algebra read off `tmesim/validation/calibration.py` or a run executed by `experiments/alpha_identifiability.py`.

> ### Headline
> 
> **`alpha` and the anchor drug's binding prior are not separately identified — only their product is.** Perturbing doxorubicin's raw binding input by a factor `f` moves `alpha` by exactly `1/f`; across six perturbation factors spanning 100x the product `alpha·b_dox` is invariant to **2.1e-16** relative error and `L_dox` to **0.0e+00**. The anchor pins one number, `k_dox = 3.331353e-02 /s`. It pins nothing else.
> 
> **Any improvement to the anchor drug's chemistry is re-exported as a global multiplier on every OTHER drug.** That is not a side effect; it is the only channel through which `b_dox` can act, because the anchor equation forces `alpha·b_dox` to stay put.
> 
> **Worse: `alpha` is a gauge coordinate, not a measurement.** Scaling *every* drug's binding prior by `t` and refitting moves `alpha` over four orders of magnitude (1.4861 to 14860.82) while leaving every downstream number — every `k`, every `L`, the out-of-sample prediction — identical to 2.8e-14 µm. Publishing `alpha = 148.60815765445545` to 17 figures publishes a coordinate.
> 
> **And the branch would have broken the project's own falsification test.** Under the reported anchor-binding drop, `calibration.revalidate()`'s held-out mitoxantrone prediction falls 48.2 → 26.4 µm and exits the published 36–108 µm band (`PASS` → `PASS(3x)`).
> 
> **The choice of anchor DRUG moves the whole panel.** Refitting on mitoxantrone — the only other drug in this repo's cited compendium with a measured perivascular half-distance — gives `alpha = 66.29` against the published `148.61`, a **0.446x** shift, and every drug's penetration length moves by **1.480–1.496x**.
> 
> **Does it cancel? The question has two answers, because the perturbation is not the one it looks like.** A uniform ×3.3478 on every drug’s sink is common-mode and cancels almost completely in cross-drug ratios (median **73.8%** cancelled against the anchor drug; ranking fully intact). But **what a change to `b_dox` actually does is move nine drugs and leave the tenth — the anchor — exactly where it was**, and against the anchor that cancels **0.0%**.
> 
> **Absolutes do not survive either way** (kill 0.272–0.595×), and the **3-D:2-D kill ratio** — the one external test absolute kill passes — still takes 37.7% of the perturbation, because renormalising to the best-perfused cell removes a global concentration but not a change in gradient SHAPE. The cross-drug **ranking is the one thing that holds**, with one caveat: the anchor-exempt arm flips top-1 across a 0.03% margin — a near-tie, not a reversal.

---

## 1. What is being fit, against what, with how many degrees of freedom

`calibration.calibrate_alpha()` is a closed-form solve of **one scalar equation in one scalar unknown**:

```
k_d(alpha) = c_d + alpha * b_d                      (per drug d)
L_d(alpha) = sqrt(D_eff,d / k_d(alpha))

FIT:   L_dox(alpha) = L*                            ONE equation
SOLVE: alpha = (D_eff,dox / L*^2 - c_dox) / b_dox    ONE unknown
```

with, on the doxorubicin anchor as the code actually computes it:

| symbol | meaning | value |
|---|---|---:|
| `L*` | fit target, e-folding length | 64.92128 µm |
| `D_eff,dox` | `D·(1 − 0.6·ECM)` at ECM = 0.4 | 140.40895 µm²/s |
| `c_dox` | `clearance`, NOT amplified | 6.308002e-05 /s |
| `b_dox` | raw `tissue_binding_rate` prior | 2.237458e-04 /s |
| `k*` | `D_eff/L*²`, the fitted quantity | **3.331353e-02 /s** |
| `alpha` | the solved constant | **148.60815765445545** |
| `alpha·b_dox` | the anchor's amplified sink | **3.325045e-02 /s** |

**Free parameters vs independent constraints, counted exactly.** Write the panel's physics in terms of what the PDE actually sees: the per-drug amplified sink `s_d ≡ alpha·b_d`, one per drug, `N = 10` on this GBM panel. The parameter vector is `(alpha, b_1 … b_N)` — `N + 1` numbers. The map to observables is `(alpha, {b_d}) ↦ {alpha·b_d}`, which has a **one-dimensional fibre** (§1.2 verifies it): `(alpha/t, {t·b_d})` gives identical physics for every `t`. So `N + 1` parameters carry only `N` distinguishable quantities, and the calibration supplies **one** equation. The ledger:

| | count | status |
|---|---:|---|
| parameters `(alpha, b_1…b_N)` | N + 1 = 11 | — |
| gauge directions (unobservable by construction) | 1 | `alpha` alone is pure gauge |
| physically distinguishable quantities `s_d` | N = 10 | — |
| independent data constraints in `calibrate_alpha()` | **1** | the doxorubicin depth |
| → identified | **1** | `s_dox`, hence `k_dox` |
| → NOT identified | **9** | every ratio `b_d/b_dox`, `d ≠ dox` |

* **Identified — exactly one number:** `k_dox = c_dox + alpha·b_dox = 3.331353e-02 /s`. This is fixed by the anchor and is *invariant to any change in `b_dox`*. A better `b_dox` cannot move it, which is why the branch saw zero improvement at the anchor.
* **NOT identified — nine numbers:** the ratios `s_d/s_dox = b_d/b_dox`. These come entirely from `molecular/binding.py`'s logP/HBD sigmoid and have **zero** data constraint here. `results/binding_sink_refuted.md` §2 already established that this sigmoid is refuted in direction, not merely in magnitude — so the nine unconstrained numbers are not merely unmeasured, they are the ones independently shown to be wrong.
* **NOT identified — `alpha` itself:** it is a gauge coordinate (§1.2). Restricted to the two-parameter sub-problem `(ln alpha, ln b_dox)`, the Jacobian of the single fitted observable is `(1, 1)` — rank 1 in a 2-dimensional space, with a null direction `(δ ln alpha, δ ln b_dox) = (t, −t)` along which the fit is *exactly* flat.

**The consequence, algebraically.** Send `b_dox → f·b_dox` with every other drug's prior untouched. The anchor equation gives

```
alpha'      = (k* − c_dox)/(f·b_dox) = alpha / f            (exactly)
s_dox'      = alpha'·(f·b_dox) = alpha·b_dox                (invariant)
s_d'        = alpha'·b_d       = s_d / f      for d != dox  (GLOBAL RESCALE)
k_d'/k_d    = (c_d + s_d/f)/(c_d + s_d)       -> 1/f  as s_d >> c_d
```

**So `alpha` is not a physical constant that a better `b_dox` refines. It is a sink whose only job is to absorb `b_dox`, and it re-emits that absorption onto every drug it was never fit to.** `alpha`'s docstring calls it "the one factor that lifts the deliberately-conservative Stage-0 binding prior up to the true perivascular sink strength" — that reading is only available for doxorubicin. For every other drug `alpha` is a doxorubicin-specific correction applied to a different molecule.

### 1.1 Proof by perturbation (six factors, spanning 100×)

`experiments/alpha_identifiability.py` monkeypatches `calibration._compile` to scale **doxorubicin's raw binding only**, then calls the unmodified `calibrate_alpha()`. If the degeneracy is exact, `alpha_ratio × f = 1`.

| `f` (scale on `b_dox`) | refit `alpha` | `alpha`/`alpha₀` | **`alpha_ratio × f`** | `alpha·b_dox` rel. error | `L_dox` rel. error | sink multiplier on every other drug |
|---:|---:|---:|---:|---:|---:|---:|
| 0.10000 | 1486.08158 | 10.00000 | **1.000000000000000** | 0.00e+00 | 0.00e+00 | 10.00000× |
| 0.29870 | 497.51039 | 3.34780 | **1.000000000000000** | 0.00e+00 | 0.00e+00 | 3.34780× |
| 0.50000 | 297.21632 | 2.00000 | **1.000000000000000** | 0.00e+00 | 0.00e+00 | 2.00000× |
| 2.00000 | 74.30408 | 0.50000 | **1.000000000000000** | 0.00e+00 | 0.00e+00 | 0.50000× |
| 3.34780 | 44.38980 | 0.29870 | **1.000000000000000** | 0.00e+00 | 0.00e+00 | 0.29870× |
| 10.00000 | 14.86082 | 0.10000 | **1.000000000000000** | 2.09e-16 | 0.00e+00 | 0.10000× |

Every row of column 4 is 1 to within double-precision round-off. The anchor drug's penetration length does not move at all — **it cannot**, it is the fit. The entire effect of "improving doxorubicin's binding input" lands on the nine drugs that were not improved.

For the specific case reported by the prior agent — a measured-logP substitution that dropped doxorubicin's raw binding by 3.3478× — this reproduces `alpha` moving by 3.34780×, i.e. the reciprocal, to 15 decimal places. **The 0.04% agreement they observed was not a coincidence and not an approximation. It is an identity, and the residual is bounded in closed form.**

The only channel that can break the reciprocal is the anchor's *clearance*, which `molecular/binding.py::clearance_rate` also makes a function of logP, and which enters the numerator:

```
alpha'/alpha  =  (1/f) · (k* − c'_dox)/(k* − c_dox)
```

Since `c_dox/k* = 1.893525e-03`, that second factor is confined to `[0.998103, 1.001897]` — the reciprocal identity cannot be violated by more than **0.1897%** no matter what a measured logP does to doxorubicin. The observed 0.04% sits inside that bound. There is no room in this calibration for a change to the anchor drug to do anything except rescale everything else.

### 1.1.1 The closed form predicts BOTH concurrently-reported alphas, in opposite directions

Two independent workstreams hit this from opposite sides. The identity of §1 predicts both from doxorubicin's descriptors alone — no simulation, no refit — using the public `binding.tissue_binding_rate` (doxorubicin: MW 543.525, logP 0.001, HBD 6):

| reported case | source | implied `f` on `b_dox` | `alpha` predicted by the identity | `alpha` reported | rel. error |
|---|---|---:|---:|---:|---:|
| doxorubicin logP 0.001 (RDKit Crippen) -> 1.27 (HSDB/Hansch) | `results/logp_provenance.md section 5.3` | 2.16428 | **68.66404** | 68.63680 | 3.97e-04 |
| doxorubicin raw binding dropped 3.3478x | `the measured-logP branch reported to this run` | 0.29870 | **497.51039** | 497.30899 | 4.05e-04 |

Both land inside the 0.1897% clearance bound derived above. **The two reports disagree about the sign of the correction to doxorubicin's logP** — one raises it (Crippen 0.001 → Hansch 1.27, `alpha` falls to ~68.6, every other drug's `L` grows ~1.47×), the other lowers it (`alpha` rises to ~497, every other drug's `L` shrinks ~0.55×). **That disagreement is not the problem this file is about, and resolving it would not fix anything.** Whichever is right, the correction has exactly one place to go, and it is not doxorubicin.

### 1.2 `alpha` is a GAUGE quantity — its value carries no physical content

The degeneracy is not only about the anchor. Scale **every** drug's raw binding prior by the same `t` — the whole `tissue_binding_rate` output band moved bodily — and refit:

```
b_d -> t·b_d  for ALL d     =>     alpha -> alpha/t     =>     s_d = alpha·b_d  UNCHANGED
```

| `t` (scale on EVERY `b_d`) | refit `alpha` | `alpha`/`alpha₀` | `alpha_ratio × t` | worst \|Δ`L`\| over the whole panel | worst rel. Δ`k` | mitoxantrone predicted `L` |
|---:|---:|---:|---:|---:|---:|---:|
| 0.01 | 14860.81577 | 100 | 1.000000000000000 | **0.000e+00 µm** | 1.436e-16 | 48.2086416245 µm |
| 0.2 | 743.04079 | 5 | 1.000000000000000 | **2.842e-14 µm** | 2.599e-16 | 48.2086416245 µm |
| 5 | 29.72163 | 0.2 | 1.000000000000000 | **0.000e+00 µm** | 1.436e-16 | 48.2086416245 µm |
| 100 | 1.48608 | 0.01 | 1.000000000000000 | **2.842e-14 µm** | 2.599e-16 | 48.2086416245 µm |

**`alpha` moves over four orders of magnitude and every single downstream number — every `k`, every `L`, the out-of-sample mitoxantrone prediction, its verdict — is bit-identical.** `alpha` is a *coordinate*, not a measurement. The pair `(alpha, {b_d})` has an exact one-parameter gauge freedom `(alpha/t, {t·b_d})`, and only the gauge-invariant products `s_d = alpha·b_d` reach the PDE.

This is what makes the 17-significant-figure publication of `alpha = 148.60815765445545` in `README.md`, `API_CONTRACT.md`, `tmesim/core/config.py` and `tmesim/validation/revalidation.py` misleading in kind rather than merely in precision. **The number is not wrong; it is not the sort of thing that can be right.** It is only meaningful jointly with the exact `binding.py` sigmoid that produced `b_dox`, and it silently changes whenever that sigmoid does — which is precisely the event that triggered this investigation. The quantity worth freezing and publishing is the gauge-invariant one:

> `s_dox = alpha·b_dox = 3.325045e-02 /s`  (equivalently `k_dox = 3.331353e-02 /s`)

### 1.3 What the calibration's own out-of-sample test can and cannot see

`calibration.revalidate()` predicts mitoxantrone with `alpha` frozen and it PASSES *inside* the published band. What does that test actually constrain? Mitoxantrone's predicted `L` depends on `alpha·b_mitox = s_dox·(b_mitox/b_dox)`. So it constrains the **ratio** `b_mitox/b_dox` — one of the N−1 gauge-invariant quantities — and is completely blind to the gauge (§1.2 proves that empirically: the mitoxantrone column there does not move at all).

**It is NOT blind to an anchor-only change, and this is the finding with the sharpest consequence for the branch.** Changing `b_dox` alone changes the ratio `b_mitox/b_dox`, so it moves the one out-of-sample prediction this project has:

| `f` on `b_dox` only | `alpha` | mitoxantrone predicted `L` | ratio to the 72.1 µm mid | verdict | TMZ/dox fold |
|---:|---:|---:|---:|:--:|---:|
| 1.00000 | 148.60816 | 48.2086 µm | 0.6683 | **PASS** | 4.0747 |
| 0.29870 | 497.51039 | 26.3571 µm | 0.3654 | **PASS(3x)** | 2.2421 |
| 10.00000 | 14.86082 | 151.7658 µm | 2.1039 | **PASS(3x)** | 11.8986 |

**The measured-logP correction to doxorubicin knocks the module's own out-of-sample test out of the published band.** At the reported 3.3478× binding drop, mitoxantrone's predicted `L` falls 48.21 → 26.36 µm against a published e-fold-equivalent band of 36.1–108.2 µm. It is **below the band**, and survives only under the looser 3× order-of-magnitude lens (`PASS` → `PASS(3x)`). `results/calibration.md`'s headline sentence — "the predicted penetration falls **inside** the band — the one constant transfers to a second drug it was never fit to" — would no longer be true. This is an independent, code-internal reason the branch is not mergeable, and it is stronger than a magnitude argument because it is the project's own declared falsification test.

**The structural blind spot that remains.** Nothing in this fit — not the anchor, not the out-of-sample test — sees the gauge, and nothing sees `PS` either (`results/transvascular_exchange_revalidation.md` §6.2: with `transvascular_exchange = False`, "PS appears nowhere in that formula at all"). Two independent unidentified directions live in the same one-equation fit.

---

## 2. Blast radius of the ANCHOR CHOICE

`tmesim/validation/pk_anchors.py` is this repo's own cited compendium of measured tissue-PK anchors. **It contains exactly three drugs with a measured perivascular half-distance — doxorubicin, mitoxantrone, topotecan — and the last two get their number from a range the source paper reports *jointly for all three drugs*.** Temozolomide has no penetration-depth anchor anywhere in the repo (its only entry is a brain:plasma AUC partition, PMID 19861433, which is a different physical quantity and is flagged as such in that file). So the anchor cannot be swapped to temozolomide without inventing data. The swap below uses mitoxantrone, and every target is a published half-distance converted with the module's own `d_half/ln2` rule.

*(Data-hygiene note found while doing this: the repo cites PMID 23680920 two different ways — `pk_anchors.py` as "Patel KJ, Tredan O, Tannock IF … 2013;72(1):127-38" and `calibration.py` as "Rocca … Tannock … 2013;72(3):591-601". Same PMID, different first author, volume issue and pages. One of them is wrong. Not load-bearing for anything below — the numeric anchor is identical in both — but it should be reconciled.)*

| anchor | target `L*` | refit `alpha` | vs baseline | `L_dox` | `L_tmz` | `L_tmz/L_dox` | drugs over the 5e-3 /s ceiling | worst |
|---|---:|---:|---:|---:|---:|---:|:--:|---:|
| doxorubicin, Primeau 40-50 um midpoint (THE BASELINE) | 64.92 µm | **148.608** | 1.0000× | 64.92 µm | 264.53 µm | 4.075 | 9/10 | 17.8× |
| doxorubicin, Primeau band LOW end | 57.71 µm | **188.157** | 1.2661× | 57.71 µm | 235.57 µm | 4.082 | 9/10 | 22.6× |
| doxorubicin, Primeau band HIGH end | 72.13 µm | **120.319** | 0.8096× | 72.13 µm | 293.33 µm | 4.066 | 9/10 | 14.4× |
| doxorubicin, Patel/Tredan/Tannock 25-75 um midpoint | 72.13 µm | **120.319** | 0.8096× | 72.13 µm | 293.33 µm | 4.066 | 9/10 | 14.4× |
| MITOXANTRONE, 25-75 um midpoint (ANCHOR-DRUG SWAP) | 72.13 µm | **66.292** | 0.4461× | 97.09 µm | 391.43 µm | 4.032 | 7/10 | 7.9× |
| MITOXANTRONE, band LOW end | 36.07 µm | **265.617** | 1.7874× | 48.58 µm | 198.71 µm | 4.090 | 10/10 | 31.8× |
| MITOXANTRONE, band HIGH end | 108.20 µm | **29.380** | 0.1977× | 145.45 µm | 573.04 µm | 3.940 | 4/10 | 3.5× |

**The anchor-drug swap on its own moves `alpha` by 0.4461×** (148.608 → 66.292). Two independent things are mixed into that number and both matter:

1. the target length differs (64.92 vs 72.13 µm), because the two papers report different half-distances; and
2. the anchor drug's own `b`, `c` and `D` differ, and *those* enter only through the degeneracy of §1.

**Per-drug penetration lengths under the mitoxantrone anchor, relative to the published panel:**

| drug | `L` (baseline) | `L` (mitoxantrone anchor) | fold |
|---|---:|---:|---:|
| temozolomide | 264.53 µm | 391.43 µm | 1.4797× |
| lomustine | 100.46 µm | 149.89 µm | 1.4921× |
| carmustine | 139.52 µm | 207.76 µm | 1.4890× |
| procarbazine | 100.12 µm | 149.48 µm | 1.4930× |
| cisplatin | 142.93 µm | 212.86 µm | 1.4893× |
| methotrexate | 75.24 µm | 112.47 µm | 1.4948× |
| etoposide | 75.33 µm | 112.57 µm | 1.4943× |
| vincristine | 36.35 µm | 54.40 µm | 1.4963× |
| irinotecan | 47.50 µm | 71.04 µm | 1.4955× |
| doxorubicin | 64.92 µm | 97.09 µm | 1.4955× |

And the *published band itself* is wider than the anchor-drug swap. Taking the two ends of the SAME doxorubicin measurement (40 µm and 50 µm) gives `alpha` = 188.16 and 120.32 — a 1.564× spread from the reported uncertainty of the single anchor, before any question of which drug to anchor on. Taking the mitoxantrone band's two ends (25–75 µm) gives 29.38 to 265.62, a **9.0× range**.

**`alpha = 148.60815765445545` is quoted to 17 significant figures in `README.md`, `API_CONTRACT.md`, `tmesim/core/config.py` and `tmesim/validation/revalidation.py`. Its defensible precision is about one significant figure.** The digits are reproducible, not meaningful: they are a deterministic function of an RDKit logP and a 40–50 µm literature range.

### 2.1 The two measured anchors are NOT mutually consistent, and the out-of-sample PASS is an artefact of which band is wider

Fit on one anchor, then ask whether the OTHER anchor's published band is hit. The answer depends entirely on which way round you do it:

| fit on | predict | predicted `L` | published e-fold band | band width | verdict |
|---|---|---:|---:|---:|:--:|
| doxorubicin (64.92 µm) | mitoxantrone | 48.21 µm | 36.1–108.2 µm | **3.00×** | **INSIDE** |
| mitoxantrone (72.13 µm) | doxorubicin | 97.09 µm | 57.7–72.1 µm | **1.25×** | **OUTSIDE** |

**The model cannot fit both measurements at once.** Anchored on doxorubicin it hits mitoxantrone; anchored on mitoxantrone it puts doxorubicin at 97.1 µm, **1.35× above the top** of doxorubicin's own published band. That is a real, if modest, falsification of the single-sink model on the only two penetration measurements it has — and it is invisible in the current direction of the test.

**And the direction that passes is the easy one.** The held-out mitoxantrone band spans **3.00×** (36–108 µm); the doxorubicin band it was fitted to spans only **1.25×**. `pk_anchors.py`'s own note explains why: the 25–75 µm range is reported *jointly for three drugs at one timepoint*, not as a mitoxantrone-specific value, and the file grades that entry `confidence: medium` for exactly this reason. **A prediction landing inside a 3× band is a weak pass**, and `results/calibration.md` presents it as the module's headline validation without quoting the band width. The honest version of that claim is: the one constant transfers to a second, chemically-adjacent drug within a factor of 3, in the only direction where the target is 3× wide.

### 2.2 The internal contradiction gets worse under every alternative anchor

`results/binding_sink_refuted.md` §4 showed that after the `alpha` amplification most calibrated drugs breach `molecular/binding.py`'s own documented [1e-5, 5e-3] /s band. The counts in the table above are on the 10-drug `experiments/real_drug_panel.py` panel, so they are not the same denominator as that file's 19-drug or 13-drug counts — but the point is the *movement*: from 4/10 (worst 3.5×) to 10/10 (worst 31.8×) across nothing but the two ends of one published range. **The repo's own internal self-contradiction is a function of an arbitrary anchor choice.** That is a stronger statement than the contradiction itself: a quantity whose breach of its own documented bounds can be dialled from 4/10 to 10/10 by picking a different end of a literature range is not a measured quantity.

---

## 3. Does it cancel?

**Setup.** Real 10x Xenium GBM section, `dx = 12` µm, `dt = 10` s, 6 h window, convection on, seed 0 — identical to `tmesim/validation/kill_calibration.py`, whose `prepare_tissue`, `unit_transport`, `decompose` and `kill_with` are **imported, not re-implemented**. 20 transport solves at `alpha = 148.60816` and `alpha = 497.51039` (**3.3478×**). Sanity: the reduced path reproduces the full engine (`simulate_response` on temozolomide gives kill 0.018352, matched here to all printed digits) and reproduces `results/transvascular_exchange_revalidation.md` §2's published one-way numbers exactly (cisplatin 2.7122, temozolomide 1.1403, lomustine 1.3881).

### 3.1 The perturbation is TWO different perturbations, and they answer the question differently

This is the part the framing of the question gets wrong, and it matters.

| arm | what moves | is it common-mode? |
|---|---|---|
| **A — uniform** | every drug’s sink × 3.3478 | YES — this is the gauge direction of §1.2 |
| **B — anchor** | nine drugs’ sinks × 3.3478, **doxorubicin exempt** | NO — one drug in ten does not move |

**Arm B is what actually happens when `b_dox` changes**, and it is *not* a global multiplier. Under `b_dox → b_dox/g`, `alpha → g·alpha`, the anchor’s `s_dox = alpha·b_dox` and its `clearance` are both untouched, so doxorubicin’s `DrugParams` are bit-identical and its solve *is* the baseline solve; the other nine get exactly the perturbed solve. Both arms are therefore exact reconstructions from the same 20 runs, with no extra simulation and no approximation.

**A global multiplier on every drug’s sink is common-mode. A change to the ANCHOR drug’s binding is not — it is common-mode on 9 of 10 drugs with the 10th pinned.** That distinction is the whole answer to “does it cancel?”

`alpha` scales only the `binding` half of `k = clearance + binding`, so even in arm A `k'/k = (c + g·s)/(c + s) ≤ g`:

| drug | `k` baseline | `k` ×g arm | `k` fold | purity `ln(k'/k)/ln g` (arm A) | purity (arm B) |
|---|---:|---:|---:|---:|---:|
| temozolomide | 3.0191e-03 /s | 9.9713e-03 /s | 3.3027× | 0.98879 | 0.98879 |
| lomustine | 1.9518e-02 /s | 6.5087e-02 /s | 3.3347× | 0.99677 | 0.99677 |
| carmustine | 1.0473e-02 /s | 3.4843e-02 /s | 3.3270× | 0.99483 | 0.99483 |
| procarbazine | 2.0071e-02 /s | 6.6981e-02 /s | 3.3372× | 0.99736 | 0.99736 |
| cisplatin | 8.7864e-03 /s | 2.9237e-02 /s | 3.3275× | 0.99496 | 0.99496 |
| methotrexate | 2.6697e-02 /s | 8.9215e-02 /s | 3.3417× | 0.99849 | 0.99849 |
| etoposide | 2.3921e-02 /s | 7.9906e-02 /s | 3.3405× | 0.99818 | 0.99818 |
| vincristine | 8.9146e-02 /s | 2.9823e-01 /s | 3.3454× | 0.99942 | 0.99942 |
| irinotecan | 6.0239e-02 /s | 2.0140e-01 /s | 3.3434× | 0.99890 | 0.99890 |
| doxorubicin | 3.3314e-02 /s | 1.1138e-01 /s | 3.3434× | 0.99890 | 0.00000 |

Arm A purity 0.9888–0.9994 — common-mode to within 1.12%, genuinely the class of error that has been cancelling in this repo. Arm B’s doxorubicin row is **0.00000**: the anchor does not move at all.

### 3.2 Absolute outputs — the full size of the error

Arm A and arm B differ only in the doxorubicin row (arm B’s is 1.0000× by construction).

| drug | mean tumour exposure (µM per µM plasma) | fold, arm A | fold, arm B | kill @ nominal dose | fold, arm A | fold, arm B |
|---|---:|---:|---:|---:|---:|---:|
| temozolomide | 1.1403 → 0.34957 | **0.3066×** | 0.3066× | 0.018352 → 0.004997 | **0.2723×** | 0.2723× |
| lomustine | 1.3881 → 0.41388 | **0.2982×** | 0.2982× | 0.237065 → 0.089703 | **0.3784×** | 0.3784× |
| carmustine | 2.536 → 0.76215 | **0.3005×** | 0.3005× | 0.403116 → 0.220267 | **0.5464×** | 0.5464× |
| procarbazine | 1.2272 → 0.36562 | **0.2979×** | 0.2979× | 0.024945 → 0.007307 | **0.2929×** | 0.2929× |
| cisplatin | 2.7122 → 0.81564 | **0.3007×** | 0.3007× | 0.436167 → 0.259407 | **0.5947×** | 0.5947× |
| methotrexate | 0.023426 → 0.0069273 | **0.2957×** | 0.2957× | 0.019192 → 0.005459 | **0.2845×** | 0.2845× |
| etoposide | 0.045267 → 0.013394 | **0.2959×** | 0.2959× | 0.098564 → 0.028851 | **0.2927×** | 0.2927× |
| vincristine | 0.0070309 → 0.0019992 | **0.2843×** | 0.2843× | 0.417815 → 0.221485 | **0.5301×** | 0.5301× |
| irinotecan | 0.06909 → 0.020049 | **0.2902×** | 0.2902× | 0.064505 → 0.018646 | **0.2891×** | 0.2891× |
| doxorubicin | 0.017128 → 0.0050402 | **0.2943×** | 1.0000× | 0.259491 → 0.102637 | **0.3955×** | 1.0000× |

**Nothing quoted as an absolute survives either arm.** Kill moves 0.272–0.595× and mean exposure 0.284–0.307×.

### 3.3 Ratio-type outputs — what survives

`surviving fraction = |Δln(ratio-type output)| / |Δln(absolute output)|`; `cancellation = 1 − surviving`. Same construction as `results/transvascular_exchange_revalidation.md`. Two normalisers: against the anchor drug, and against the panel geometric mean (what a relative score effectively uses).

**Arm A (uniform / gauge direction) — it cancels essentially completely.**

| drug | kill/kill_dox | surviving | **cancelled** | kill/panel-geomean | surviving | **cancelled** |
|---|---:|---:|---:|---:|---:|---:|
| temozolomide | 0.0707 → 0.0487 | 28.7% | **71.3%** | 0.1673 → 0.1228 | 23.8% | **76.2%** |
| lomustine | 0.9136 → 0.8740 | 4.6% | **95.4%** | 2.1610 → 2.2036 | 2.0% | **98.0%** |
| carmustine | 1.5535 → 2.1461 | 53.5% | **46.5%** | 3.6746 → 5.4109 | 64.0% | **36.0%** |
| procarbazine | 0.0961 → 0.0712 | 24.5% | **75.5%** | 0.2274 → 0.1795 | 19.3% | **80.7%** |
| cisplatin | 1.6809 → 2.5274 | 78.5% | **21.5%** | 3.9759 → 6.3724 | 90.8% | **9.2%** |
| methotrexate | 0.0740 → 0.0532 | 26.2% | **73.8%** | 0.1749 → 0.1341 | 21.1% | **78.9%** |
| etoposide | 0.3798 → 0.2811 | 24.5% | **75.5%** | 0.8985 → 0.7087 | 19.3% | **80.7%** |
| vincristine | 1.6101 → 2.1579 | 46.1% | **53.9%** | 3.8086 → 5.4408 | 56.2% | **43.8%** |
| irinotecan | 0.2486 → 0.1817 | 25.3% | **74.7%** | 0.5880 → 0.4580 | 20.1% | **79.9%** |
| doxorubicin *(normaliser)* | 1.0000 → 1.0000 | — | — | 2.3654 → 2.5213 | 6.9% | **93.1%** |

**Arm B (what the branch actually does) — against the anchor it does not cancel at all.**

| drug | kill/kill_dox | surviving | **cancelled** | kill/panel-geomean | surviving | **cancelled** |
|---|---:|---:|---:|---:|---:|---:|
| temozolomide | 0.0707 → 0.0193 | 100.0% | **0.0%** | 0.1673 → 0.1119 | 30.9% | **69.1%** |
| lomustine | 0.9136 → 0.3457 | 100.0% | **0.0%** | 2.1610 → 2.0084 | 7.5% | **92.5%** |
| carmustine | 1.5535 → 0.8488 | 100.0% | **0.0%** | 3.6746 → 4.9316 | 48.7% | **51.3%** |
| procarbazine | 0.0961 → 0.0282 | 100.0% | **0.0%** | 0.2274 → 0.1636 | 26.8% | **73.2%** |
| cisplatin | 1.6809 → 0.9997 | 100.0% | **0.0%** | 3.9759 → 5.8079 | 72.9% | **27.1%** |
| methotrexate | 0.0740 → 0.0210 | 100.0% | **0.0%** | 0.1749 → 0.1222 | 28.5% | **71.5%** |
| etoposide | 0.3798 → 0.1112 | 100.0% | **0.0%** | 0.8985 → 0.6460 | 26.9% | **73.1%** |
| vincristine | 1.6101 → 0.8535 | 100.0% | **0.0%** | 3.8086 → 4.9589 | 41.6% | **58.4%** |
| irinotecan | 0.2486 → 0.0719 | 100.0% | **0.0%** | 0.5880 → 0.4175 | 27.6% | **72.4%** |
| doxorubicin *(exempt — does not move)* | 1.0000 → 1.0000 | — | — | 2.3654 → 5.8098 | — | — |

**The headline cancellation numbers.** Against the anchor drug: arm A cancels a median **73.8%** of the kill error; arm B cancels a median **0.0%**. Against the panel geometric mean: arm A **79.4%**, arm B **71.5%**. Compare `results/transvascular_exchange_revalidation.md`’s 84–86% for the transvascular fix.

**Cross-drug ranking.**

| arm | quantity | Spearman ρ | Kendall τ | discordant / 45 | top-1 | top-3 set | top-1 margin base → arm |
|---|---|---:|---:|---:|:--:|:--:|---:|
| A uniform | kill_fraction | 1.0000 | 1.0000 | 0 | HELD | HELD | 4.21% → 14.62% |
| A uniform | mean exposure | 1.0000 | 1.0000 | 0 | HELD | HELD | — |
| B anchor | kill_fraction | 0.9273 | 0.8667 | 3 | **CHANGED** | **CHANGED** | 4.21% → 0.03% |
| B anchor | mean exposure | 0.9636 | 0.9111 | 2 | HELD | HELD | — |

**Rank order and relative magnitude are different claims, and only the first one survives.** Arm A preserves the kill ORDER perfectly (ρ = 1.0000, 0/45 discordant) while the kill of individual drugs *relative to doxorubicin* still moves by up to 78% of the absolute error — cisplatin’s kill:kill_dox goes 1.6809 → 2.5274. Any claim of the form “drug X kills N× as much as drug Y” is NOT protected by the ranking result. This is the same distinction `results/transvascular_exchange_revalidation.md` §4.1 insists on — “any future claim must say WHICH quantity it is about” — and it applies here in the ratio direction rather than the delivery-vs-kill direction.

* kill order, baseline:  cisplatin > vincristine > carmustine > doxorubicin > lomustine > etoposide > irinotecan > procarbazine > methotrexate > temozolomide
* kill order, arm A:     cisplatin > vincristine > carmustine > doxorubicin > lomustine > etoposide > irinotecan > procarbazine > methotrexate > temozolomide
* kill order, arm B:     doxorubicin > cisplatin > vincristine > carmustine > lomustine > etoposide > irinotecan > procarbazine > methotrexate > temozolomide

**Arm B changes top-1** — but read the margin before reacting to it. The baseline leader wins by 4.21% of its own kill; under arm B the new leader wins by 0.03%. **This is a near-tie, not a reversal.** It is categorically unlike `results/transvascular_exchange_revalidation.md` §4.1, where two swapped pairs crossed gaps of 40.5% and 35.6% of the leader’s kill — roughly 50-sigma. The honest statement is that arm B leaves the ranking intact except at a boundary that was already inside the noise.

### 3.4 The 3-D:2-D kill ratio — 37.7% of the sink error reaches it

`results/kill_calibration.md`’s headline is a **ratio** (`F_penetration × F_state`), and passing it is the only external validation absolute kill has. `F_state` is transport-independent and exactly invariant here. `F_penetration` is not: it is `mean_i p_hill(E_i) / p_hill(C_bath)` with the field renormalised so the p99 tumour cell sits at the bath — a renormalisation that removes a *global concentration* factor but **cannot** remove a change in the SHAPE of the perivascular gradient, and changing the sink changes exactly that shape. Arm A and arm B are identical here except for doxorubicin, which does not move in arm B.

| drug | `F_penetration` @1×IC50 base → ×g | fold | 3-D:2-D `ratio_exact` base → ×g | fold |
|---|---:|---:|---:|---:|
| temozolomide | 0.5993 → 0.4999 | **0.8342×** | 0.3348 → 0.2802 | **0.8369×** |
| lomustine | 0.4412 → 0.3396 | **0.7697×** | 0.2478 → 0.1914 | **0.7723×** |
| carmustine | 0.4932 → 0.3961 | **0.8032×** | 0.2765 → 0.2228 | **0.8060×** |
| procarbazine | 0.4406 → 0.3387 | **0.7686×** | 0.2475 → 0.1909 | **0.7711×** |
| cisplatin | 0.4980 → 0.4028 | **0.8090×** | 0.2791 → 0.2266 | **0.8118×** |
| methotrexate | 0.3939 → 0.2851 | **0.7237×** | 0.2216 → 0.1609 | **0.7258×** |
| etoposide | 0.3945 → 0.2853 | **0.7231×** | 0.2219 → 0.1610 | **0.7252×** |
| vincristine | 0.2595 → 0.1642 | **0.6328×** | 0.1465 → 0.0928 | **0.6338×** |
| irinotecan | 0.3129 → 0.2025 | **0.6471×** | 0.1765 → 0.1144 | **0.6484×** |
| doxorubicin | 0.3688 → 0.2557 | **0.6932×** | 0.2077 → 0.1443 | **0.6950×** |

**At worst 37.7% of the sink perturbation survives into the 3-D:2-D ratio**, i.e. it cancels 62.3% — real protection, but not the 84–86% the transvascular fix got, and the residual lands directly on the number `kill_calibration.md` compares to a 2–5× literature band. Temozolomide’s predicted 3-D:2-D ratio moves 0.3348 → 0.2802. The band is wide enough (`kill_calibration.md` §1 concedes the 2-D reference alone has a ~19× IQR) that the verdict does not flip at this perturbation size — but “PASSES for all 3 drugs” should not be read as robust to the calibration, because the margin it passes by is comparable to the anchor’s own published uncertainty (§2).

---

## 4. What would identify it

The degeneracy is `rank(constraints) = 1` against `N+1` parameters. It is broken only by **adding independent constraints of the right kind**, never by adding parameters. Specifically:

**(a) The one that breaks the `alpha`/`b_dox` split — a direct binding measurement on the anchor drug.** The split is unidentifiable from *any* number of penetration depths, because `alpha` and `b_dox` enter every penetration length only through `s_d = alpha·b_d`. What is needed is a measurement of the sink *itself*, in the right compartment: the free-interstitial first-order loss rate `k` (units 1/s), or equivalently an unbound tissue-to-plasma partition `Kp,uu` plus a transvascular `PS`, for doxorubicin in tumour tissue. `results/binding_sink_refuted.md` §1 already identifies the exact measurement class: intratumoural microdialysis with a resolvable washout phase, the design Blakeley *et al.* (PMID 18787762) used for methotrexate ECF elimination. One such number for doxorubicin converts `alpha` from a fitted scale into a *testable prediction*, and the current model already fails that test by ~94× on the one drug where it exists.

**(b) The one that breaks the N−1 unconstrained ratios — a second drug with an independently measured penetration depth, in the same tissue, from a DIFFERENT chemical class.** This is the constraint the mitoxantrone test is trying and failing to be. It fails on two counts, both documented in the code: mitoxantrone is an anthracenedione DNA-intercalator chemically adjacent to doxorubicin (`calibration.py`'s own honesty caveat), and its 25–75 µm number comes from a range PMID 23680920 reports jointly for doxorubicin, mitoxantrone and topotecan — so it is not fully independent of the anchor either. The requirement is specific: **a perivascular immunofluorescence-vs-CD31-distance profile for a WEAKLY-bound, deeply-penetrating drug — temozolomide is the right one — measured in glioma tissue by the Tannock-lab method.** That single measurement would constrain `b_tmz/b_dox`, which currently rests on nothing but a logP/HBD sigmoid, and it spans the widest arm of the panel. It does not exist today; §2 above had to be run on mitoxantrone for exactly that reason.

**(c) What would NOT help, and must not be attempted.** Fitting `alpha` jointly with anything else to the same doxorubicin anchor; fitting a per-drug correction to the same anchor; or refitting `alpha` on a breast anchor (explicitly ruled out by `results/binding_sink_refuted.md` §6). Each adds parameters to a rank-1 constraint and produces a model that fits better and predicts no more. Two anchors and one free constant is already the maximum defensible ratio here — it is at rank, not below it.

**(d) The measurement that would make the whole term obsolete rather than identified.** `results/binding_sink_refuted.md` §2 establishes that the residual is **active efflux and intracellular trapping**, which a first-order interstitial loss term does not physically represent and no logD/logP heuristic can rank. Identifying `alpha` better makes a structurally wrong term more precisely wrong. The measurement that changes the model rather than its constants is a per-drug efflux-transporter substrate characterisation (P-gp/BCRP/MRP) coupled to a two-compartment interstitial/intracellular sink — at which point `alpha` has no job left to do.

---

## 5. Verdict

**`alpha` is not a fitted physical constant; it is a doxorubicin-shaped hole into which every revision of the anchor drug's chemistry disappears, and out of which it re-emerges as a global rescale of nine other drugs.** The published value 148.60815765445545 is exactly reproducible and *not* a measurement of anything: it is a gauge coordinate (§1.2), it is defensible to about one significant figure even as a coordinate (§2), and the single physical quantity behind it is `k_dox`.

**On mergeability of the measured-logP branch: the prior agent’s verdict stands, and this sharpens why.** The branch’s zero improvement at the anchor is not a disappointing result; it is a proof that the change had nowhere to go except onto the other nine drugs. And §1.3 adds an independent, code-internal reason: the change pushes `calibration.revalidate()`’s held-out mitoxantrone prediction out of its published band.

**On the product.** The cross-drug ranking — the repo’s declared primary output — holds under the common-mode arm (kill Spearman 1.0000, 0/45 discordant, top-1 held) and flips top-1 across a near-tie under the anchor-exempt arm (3/45 discordant, margin 0.03%). That is consistent with `results/transvascular_exchange_revalidation.md` §4.1’s explanation of *why* — common-mode survives, drug-specific does not — and this perturbation is common-mode on 9 of 10 drugs.

**But the honest reading is not “it cancels”.** Everything quoted as an absolute moves by 0.272–0.595×. Against the anchor drug itself, the arm that actually occurs cancels 0.0%. And the 3-D:2-D ratio test — a ratio, and the only external validation absolute kill has — still takes 37.7% of the perturbation. **The product is safe only in the one output it has always been safe in, and that safety was never in question here.**

**Nothing here licenses changing `alpha`.** The correct action is to stop describing it as calibrated-and-validated, and to stop publishing it as a number. It is a gauge coordinate attached to *one* identified physical quantity, `k_dox`, fixed by one drug, one measurement and one half-distance→e-fold convention, with a rank-1 constraint and nine unconstrained ratios supplied by a heuristic this repo has already refuted in direction.

---

*Generated by `experiments/alpha_identifiability.py`. Anchor values read from `tmesim/validation/pk_anchors.py`; PDE machinery imported from `tmesim/validation/kill_calibration.py`. No file under `tmesim/` was modified.*
