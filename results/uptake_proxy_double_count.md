# Does `uptake_proxy` already carry the capacity ordering — and is the ranking failure double-counted?

**A sensitivity study on top of two already-published results. No model change.** Nothing in `tmesim/` was modified, no default moved, and `calibrate_alpha()` is still `22.223477711172563` — asserted against live code at the start of every stage.

* script: `experiments/uptake_proxy_double_count.py` (`--stage role | overlap | residual | report`)
* frozen record: `results/sweep/uptake_proxy_double_count.json`
* the two findings this answers to: `results/brain_slice_anchor.md` (TEST B, **FAIL**) and `results/capacity_ranking_robustness.md` sec. 6.2, which raised the double-count possibility and explicitly left it open
* the `E[kill]`-versus-capacity tables are **re-used, not re-run**, from `results/sweep/capacity_ranking_robustness.json`. Sec. 5 states exactly what that inherits and what was done to license it.

---

## 0. The answer

> ## NO — and the reason is a physics distinction, not a statistic.
>
> **`uptake_proxy` is wired into a SINK, not a capacity.** It returns `drug.binding`, and `drug.binding` reaches the solver at exactly one place in the shipped transport path — `physics/diffusion.py:315`, `k = clearance + binding` — as a first-order removal term with no return path and no bound pool (the opt-in `molecular/immune_sink.py` route builds the identical scalar and adds to it; same form). Measured on the shipped solver rather than argued: multiplying `binding` by 10 moves the steady-state concentration to **0.1089×** its shipping value, while a capacity of 10 leaves it at **0.999999999998×**. A capacity stores solute and gives it back; this term destroys it. **They are different operators, and the double-counting defence is weaker than it looks before you check which one is implemented.**
>
> **The defence cannot reach the arm that does most of the damage.** The published failure includes a *uniform*-capacity arm — every compound given the SAME value, so zero per-compound ordering — and the full published order breaks there at **R = 36.3** (GBM 10-drug) and **R = 5.28** (breast 10-drug), both inside the measured range and below the best available constant of 39. **No overlap between `uptake_proxy` and capacity, however strong, can shrink a term that has no per-compound content.**
>
> **Pushed to its limit, the defence makes the result WORSE on two of three panels.** Give every compound a capacity exactly proportional to its `uptake_proxy` — the maximal double-count, where the entire per-compound shape is the engine's own and none of it can be new physics — and the published order breaks at a geometric-mean capacity of **14.7** (GBM 10-drug, against 36.3 uniform) and **10** (GBM 5-drug, against **435** uniform). Re-applying the engine's own ordering is not a null operation; it is a reordering force in its own right.
>
> **Against the measured quantity there is no test to run: n = 1.** Exactly 1 of the 35 rat brain-slice compounds appears on either shipping panel (`docetaxel`), the frozen dataset carries no SMILES column, and `uptake_proxy` is a function of Crippen logP and HBD — it cannot be evaluated without a structure. **Whether the proxy orders compounds the way measured capacity does is unknown, and this study does not resolve it in either direction** (sec. 2.1).
>
> **One thing the defence genuinely explains, reported here and not buried:** the breast panel's top-1 change from `vinorelbine` to `doxorubicin` under the structural `K_b` rule. On that panel `uptake_proxy` and `K_b` really are correlated (Spearman **+0.780**, p = 0.008) and applying the rule *residually* — stripping what the proxy already predicts — restores `vinorelbine` as the top pick; and on the GBM 10-drug panel the same operation restores the full published order. **Those two named findings are substantially double-counted and should be discounted.** The breast full ordering is not restored, and no uniform, proportional or sampled arm moves at all (sec. 4).

**What changes in what may be said.**

* The `capacity_ranking_robustness.md` sec. 6.2 sentence — *"the sensitivity measured above is partly a re-application of physics the model already has"* — is **too generous and should be narrowed.** It is true of the two named structural-`K_b` arms re-tested here, and of nothing else. It is false of the uniform arms, false of the sampled arms, and the proportional arm points the other way.
* The published effect sizes should **not** be relabelled an upper bound on this basis. Sec. 3 quantifies the correction and it is small, two-signed, and smaller than the disagreement between the study's own arms.
* A **new** defect is established here that is not in either prior document: the engine holds a per-compound quantity that its own docstring calls *"a capacity ordering"* and applies it as an irreversible sink. That is not a double-count. It is a mis-assignment, and sec. 3.4 says why it cannot simply be moved.

---

## 1. What `uptake_proxy` actually is

| | |
|---|---|
| returns | `drug.binding` — nothing else |
| which is set from | `molecular.binding.tissue_binding_rate(descriptors)` |
| computed from | `logp`, `hbd` — RDKit Crippen logP and H-bond-donor count, and nothing else |
| functional form | `score=(logP-2.5)/2.5+(HBD-2)/4 ; rate=1e-5*(5e-3/1e-5)**sigmoid(score)` |
| range | clamped to `1e-05` – `5e-03` 1/s |
| per-compound? | **yes** — a pure function of two descriptors |
| enters the physics at | `physics/diffusion.py:315  k = clearance + binding` — the only place it reaches the transport solver |
| with what functional form | `-k*C   (first-order removal, no return path)` |
| bound pool / return path | **False** |

`alpha` multiplies it into `binding` in `calibration.calibrated_drugparams` (`binding := alpha * uptake_proxy`, `clearance := k_passive`), and `binding` is summed with `clearance` into the PDE reaction rate. **It is per-compound and it dominates `k`**: across the GBM panel the proxy spans **30.1×** and the resulting `k` spans **26.8×**, with Spearman +1.000 and a log–log slope of 0.966 between them. **The engine's per-compound sink ordering IS the `uptake_proxy` ordering.** That much of the brief is exactly right.

### 1.1 Capacity or sink? Measured on the shipped solver, not asserted

This is the distinction everything else turns on. A tissue capacity `s_eff` stores solute reversibly: it retards the approach to equilibrium and leaves the equilibrium alone (`capacity_ranking_robustness.md` sec. 2b). A sink removes solute and moves the equilibrium. A homogeneous probe was run to steady state through `sim.engine.simulate_response` three ways:

| arm | steady-state mean C | ratio to shipping |
|---|---:|---:|
| shipping (`k` = 0.0011 /s) | 236.083052 | 1 |
| capacity R = 10 (D, k, S all ÷ 10) | 236.083052 | **0.999999999998** |
| `binding` × 10 | 25.712016 | **0.1089** |

**A capacity is invisible at equilibrium. `binding` is not.** The `binding`×10 arm lands at 0.1089, which is `k_base/k_new` to five figures — the exact signature of a first-order sink and not of a store.

**This is the unflattering half of the answer and it is stated first.** `uptake_proxy`'s own docstring says it is *"intentionally a capacity ordering, not an absolute rate"*. The module docstring **~190 lines above it, in the same file**, describes the same physics as *"a near-irreversible perivascular sink"*. Both sentences are in `tmesim/validation/calibration.py`; they describe different operators; **the code implements the sink.** So the premise of the double-counting defence — that the engine already carries this physics *as capacity* — is false as stated. What it carries is a capacity-flavoured *ordering* doing a sink's job.

## 2. Does it order compounds the way measured capacity does?

### 2.1 Against measured `Vu,brain` — n = 1, so there is no test

Of the 35 distinct rat, scoreable compounds in `data_cache/brain_slice/vu_brain.csv`, **1** appears on either shipping panel: **docetaxel**.

* `docetaxel` — measured `Vu,brain` = 777 mL·g⁻¹ (rat, Loryan 2013); `uptake_proxy` = 1.0032e-03 /s (the **largest** on the breast panel); Crippen logP 3.26; structural `K_b` = 90.06. The engine assumes a capacity of **1.0** for it.

**One point is not an ordering.** No rank statistic, no log-linear fit and no significance is computable at n = 1, and none is reported. The single point is consistent with the proxy ranking a high-capacity compound high — and a single consistent point is not evidence.

The frozen brain-slice dataset carries NO SMILES column (recorded in results/brain_slice_anchor.md sec. 4.2, where B-S9 was declared NOT COMPUTABLE for exactly this reason), 23 of its 33 retained rows are radiolabelled PET tracers, and `uptake_proxy` is a function of RDKit-computed Crippen logP and HBD -- it cannot be evaluated without a structure. Sourcing structures now would be a curation step taken after the measured column was read; that report declined it and this one declines it too. No proxy value for a measured compound is estimated, guessed or substituted here.

**The brief's framing has to be corrected here.** *"How well does `uptake_proxy` order compounds the way measured capacity does?"* is **not answerable from anything in this repository**. Everything below compares the proxy to the repo's own *structural* `K_b` rule — a QSAR fitted to liver microsomes, not a measurement on these compounds, and one `results/intracellular_retardation.md` has already refuted on measured Pruijn HT-29 transport (rank agreement 77.4 % → 51.6 %). **Agreement or disagreement with it bounds nothing about the real quantity.** A prior overlap of n ≈ 10 at p = 0.35 was not evidence of agreement and is not evidence of disagreement; that reading was correct and it is preserved.

### 2.2 Against the structural `K_b` rule

| panel | n | Spearman ρ | p | Pearson r on logs | p | R² | log–log slope |
|---|---:|---:|---:|---:|---:|---:|---:|
| GBM 10-drug | 10 | +0.333 | 0.347 | +0.630 | 0.051 | 0.397 | 0.760 |
| breast 10-drug | 10 | +0.780 | 0.008 | +0.741 | 0.014 | 0.549 | 1.422 |

**The figures I was given are right, and they are also incomplete.** `capacity_ranking_robustness.md` sec. 6.2's Spearman **+0.333** (p = 0.35) and Pearson-on-logs **+0.630** (p = 0.051) reproduce exactly (+0.3333, p = 0.3466; +0.6303, p = 0.0507). **They are the GBM panel only.** On the breast panel — the panel where the top-1 pick actually changed — the overlap is much stronger and is significant: Spearman **+0.780** (p = 0.008), Pearson-on-logs **+0.741** (p = 0.014), R² = 0.549. Quoting the GBM number alone understates the overlap exactly where it matters, and sec. 4.1 is the consequence.

**Correlation is the wrong statistic here; the SLOPE is the right one.** Capacity and sink enter the finite-window exposure only through the dimensionless group `k·T/R` (`capacity_ranking_robustness.md` sec. 2), so the capacity that cancels against the sink is one *proportional* to it — not one merely rank-correlated with it. The log–log slope of `K_b` on `uptake_proxy` is **0.76** on GBM and **1.42** on breast. Neither is 1. On breast the capacity grows nearly half again as fast as the sink that would absorb it, which is why a strong rank correlation still leaves a large reordering force.

**And the two terms contradict each other on a specific compound.** `uptake_proxy` is monotone **increasing** in logP; the `K_b` rule (`125·10^(0.072·logP² + 0.067·logP − 1.126)`) is a **parabola** in logP whose minimum sits near logP = −0.47, so it turns back up for hydrophilic compounds. `temozolomide` (logP −2.08) is therefore the **smallest** sink in the engine and sits mid-table on capacity (5th smallest of 10) under the structural rule. This is the second defect the brief anticipated: two terms in the same engine, both claiming to describe intracellular sequestration, disagreeing about the same molecule.

## 3. Is the ranking failure double-counted?

### 3.1 What CAN be double-counted, stated before any number

A per-compound capacity assignment decomposes in log space into a **level** (its geometric mean) and a **shape** (the per-compound departures from it). Only the shape is a candidate:

* **the level is not carried by the proxy.** `alpha` is fitted to a penetration *length*, and a uniform capacity is exactly invariant on penetration length (`capacity_ranking_robustness.md` sec. 2b, verified to 1e-9). Nothing in the calibration sees the capacity level, so nothing about it can have been applied twice;
* **the shape is the only candidate**, and only to the extent the proxy predicts it.

**That alone settles the largest part of the published failure**, because the uniform arm has no shape at all:

| panel | full order holds up to R = | top-1 holds up to R = |
|---|---:|---:|
| GBM 10-drug | **36.3** | > 4184 (never changes on the grid) |
| GBM 5-drug | **435** | 435 |
| breast 10-drug | **5.28** | > 4184 (never changes on the grid) |

Two of those three are below the best available constant of 39 and all three are inside the measured range [3.2, 4184]. **The double-count defence has no purchase on any of it.**

### 3.2 The strongest possible defence: let the proxy BE the capacity

Give every compound a capacity exactly proportional to its `uptake_proxy`, normalised to a geometric mean of `c`, and sweep `c`. This is the double-count at its limit: the per-compound shape is *entirely* the engine's own, so *none* of it can be new physics. If the ordering held here, the defence would be complete.

| panel | full order breaks at c = | vs. the UNIFORM arm | top-1 changes at c = |
|---|---:|---:|---:|
| GBM 10-drug | **14.7** | 36.3 | 53.4 |
| GBM 5-drug | **10** | 435 | never on [1, 4184] |
| breast 10-drug | **3.16** (clipped; 7.57 unclipped) | 5.28 | never on [1, 4184] |

**It does not hold — and on the two GBM panels it is markedly worse than a uniform capacity.** On the 5-drug sub-panel, the one `results/ranking_robustness.md` actually scores, a uniform capacity preserves the published order to R = 435; a capacity proportional to the engine's own uptake ordering breaks it at **c = 10**. Only on the breast panel is the proportional arm the more benign of the two (7.57 against 5.28), and by less than a factor of 1.5.

**Why re-applying the engine's own ordering is not a null operation.** Under `R_i ∝ k_i`, the group `k_i·T/R_i` becomes the *same* number for every compound, so every compound receives the *same* AUC attenuation. At the shipping state it does not: attenuation ranges over 0.9075–0.9965 across the GBM panel. Equalising it therefore *relatively boosts the slow-sink compounds*, which is itself a reordering. **A capacity that perfectly mirrors the sink does not cancel — it flattens, and flattening reorders.** This is the quantitative core of the answer and it is why "the effect is partly double-counted, so treat it as an upper bound" does not follow.

*Disclosure: on the breast panel the proportional rule assigns some compounds a capacity below 1 for small `c`, which is off the bottom of the swept grid and was clipped to 1. The table reports the break point both with clipping and restricted to the un-clipped region (c ≥ 7.33); the un-clipped figure is the conservative one and is the one quoted above.*

### 3.3 The sampled arm, applied residually

The published Monte Carlo draws each compound's capacity independently from the 33 measured values. The residual version strips from every draw the component linearly predictable from `log uptake_proxy` and keeps the geometric mean (sec. 3.1). Same seed (20260821), same 40 000 draws, same bootstrap, same decision-grade bar. The vectorised residualiser is asserted against a row-wise `polyfit` on four draws.

| panel | arm | P(top-1 = published leader) | P(full published order) |
|---|---|---:|---:|
| GBM 10-drug (`vincristine`) | as published | 0.730 [0.726, 0.734] | 0.0047 [0.0042, 0.0053] |
| GBM 10-drug (`vincristine`) | **residual** | 0.765 [0.762, 0.769] | 0.0055 [0.0049, 0.0061] |
| GBM 5-drug (`cisplatin`) | as published | 0.584 [0.580, 0.588] | 0.2002 [0.1969, 0.2034] |
| GBM 5-drug (`cisplatin`) | **residual** | 0.566 [0.562, 0.570] | 0.2397 [0.2363, 0.2431] |
| breast 10-drug (`vinorelbine`) | as published | 0.797 [0.794, 0.801] | 0.0006 [0.0004, 0.0008] |
| breast 10-drug (`vinorelbine`) | **residual** | 0.792 [0.789, 0.796] | 0.0000 [0.0000, 0.0000] |

**The published figures reproduce** — 0.005 / 0.200 / 0.001 against 0.0047 / 0.2002 / 0.0006 here, on the same seed. **The residual correction moves nothing across any bar**, and on the breast panel it moves the wrong way (to 0.0000). Residualising shrinks the geometric spread of the applied capacity by a factor of only **0.896** on the GBM panel, because a draw from the measured distribution is independent of `uptake_proxy` by construction: regressing it on the proxy removes only the sampling-noise overlap (expected R² ≈ 1/(n−2) ≈ 0.12 at n = 10).

**This arm is the weakest evidence in this document and is labelled so.** It cannot be strong: the sampled capacities were never linked to structure, so there was never a real double-count in it to remove. It is reported because its absence is itself the answer to "is the sampled instability inflated by the proxy?" — no, it is not.

*Disclosure: residualising can push a capacity outside the swept grid; 1.05 % of GBM-panel values were clipped to [1, 4184] (raw draws never leave it, because they ARE measured values).*

### 3.4 Is it cleanly separable? No — and what that implies is worse

The residual construction above is a **statistical** separation. No physical one exists, and that is a result rather than an evasion:

1. **The two terms are different operators (sec. 1.1).** Taking `alpha·uptake_proxy` out of `k` and re-expressing it as a capacity is not a change of variables: it moves the steady state, which a capacity provably cannot do. There is no algebra that converts one into the other.
2. **The sink is load-bearing for a different anchor.** `alpha` exists solely to make doxorubicin's penetration length hit ~64.9 µm. A capacity has **no** effect on penetration length. So the sink cannot be reduced in favour of a capacity without breaking the one anchor the engine is fitted to. **The two terms are not substitutes and cannot be traded off** — which is exactly why the double-count cannot be netted out.
3. **The implication runs against the engine, not for it.** If `uptake_proxy` is a capacity ordering, as its docstring says, then the engine is currently using a capacity ordering to set a quantity (penetration length) that a real capacity would leave untouched, and leaving the slot a real capacity belongs in (`s_eff`) pinned at 1.0. **That is not a double-count that cancels; it is a mis-assignment that would have to be undone before either term could be calibrated honestly.** It also means `alpha` = 22.22 is absorbing whatever the capacity slot is missing, which is a second reason not to read `alpha` as a physical constant.

## 4. Where the defence DOES do work — reported in full

### 4.1 The breast panel's top-1 change: substantially double-counted

`capacity_ranking_robustness.md`'s headline per-compound finding is that the structural `K_b` rule changes the breast top-1 pick from `vinorelbine` to `doxorubicin`. That is the arm where `uptake_proxy` and `K_b` are genuinely correlated (Spearman +0.780, p = 0.008, R² = 0.549 on logs). Applying the same rule **residually**:

| | full structural `K_b` | **residual** |
|---|---|---|
| spread of `s_eff` across the panel | 874.7× | 141.1× |
| top-1 | `doxorubicin` | `vinorelbine` |
| = published order? | **NO** | **NO** |

* full-rule order: doxorubicin > cisplatin > docetaxel > paclitaxel > epirubicin > vinorelbine > 5-fluorouracil > cyclophosphamide > palbociclib > tamoxifen
* residual order: vinorelbine > paclitaxel > doxorubicin > docetaxel > epirubicin > cisplatin > palbociclib > cyclophosphamide > 5-fluorouracil > tamoxifen

**The top-1 change is partly an artefact of double-counting.** Strip the part of the structural rule that `uptake_proxy` already predicts — which on this panel is over half the log-variance — and `vinorelbine` comes back as the top pick. The residual spread falls from 874.7× to 141.1×, and that is enough to undo the swap. **This is the one place in either document where the double-count defence is load-bearing, and it should be said plainly.**

**What it does not buy.** The residual arm still does **not** reproduce the published full ordering. It is a QSAR residualised against a QSAR on 10 points, with neither side measured on these compounds in this tissue. And it is one of four named arms: the uniform R = 39 arm changes the breast order too (sec. 3.1) and is untouched by any of this.

### 4.2 The GBM panels

* **GBM 10-drug** — full-rule top-1 `vincristine`, residual top-1 `vincristine`; residual order **matches** the published one (full-rule order does not match). `s_eff` spread 26.66× → 17.24×; log–log slope of `s_eff` on the proxy 0.736, R² 0.401.
* **GBM 5-drug** — full-rule top-1 `cisplatin`, residual top-1 `cisplatin`; residual order **matches** the published one (full-rule order matches). `s_eff` spread 3.006× → 2.983×; log–log slope of `s_eff` on the proxy -0.030, R² 0.004.

Two things to read off that, one in each direction. **For the defence:** on the GBM 10-drug panel residualising also puts the full published order back — the structural rule's single adjacent-pair swap (`irinotecan`/`procarbazine`) is inside the part the proxy predicts, so that swap is double-counted too. **Against it:** on the GBM 5-drug panel the proxy explains essentially none of the structural rule (R² = 0.004, log–log slope -0.030), so there is nothing there to remove — which is the honest reading of a Spearman of +0.333 at p = 0.35, and a reminder that the two GBM results come from the same rule on overlapping compounds.

**Net.** The published study ran four named per-compound arms that broke the order (the structural `K_b` rule and its access-limited variant, on the GBM 10-drug and breast panels). This study re-tested the two structural ones: residualising **restores the GBM 10-drug full order** and **restores the breast top-1** (though not the breast full order). The two access-limited variants were not re-tested — they produced identical orderings to the structural rule on every panel, so the same conclusion is likely, but it is not measured here. **That is the whole of what the double-count defence buys. It buys nothing on the uniform arms (sec. 3.1), nothing on the sampled arms (sec. 3.3), and it is contradicted by the proportional arm (sec. 3.2).**

## 5. Verified, assumed, inherited

**Verified against live code or a live solve.**

* `calibrate_alpha()` = `22.223477711172563`, asserted at the start of every stage. `two_route.DEFAULT_CELL_VOLUME_FRACTION` = 0.7892 and `s_eff` at the shipping state = 1.0, likewise.
* `uptake_proxy(dp) == dp.binding`, and `dp.binding == molecular.binding.tissue_binding_rate(molecule)` — asserted for the compiled molecule, not read off a docstring, and re-asserted for every GBM-panel drug in `stage_overlap`.
* `binding` reaches the PDE at one place as `−k·C` with no return path (`physics/diffusion.py:315`); the sink-versus-capacity distinction was measured through `sim.engine.simulate_response` (sec. 1.1) and the two arms are asserted to differ in the predicted direction.
* Every correlation, slope, spread and overlap count in sec. 2 is computed here from live compiles and the frozen CSV — none is transcribed. The three figures I was handed (+0.333, p = 0.35; +0.630, p = 0.051; `s_eff` = 1.0) all reproduce.
* **The published named-arm orderings are reproduced exactly by this script's own pipeline before any new arm is read off it** — asserted in `stage_residual`, and the published Monte-Carlo probabilities reproduce on the same seed (sec. 3.3). That is the check that licenses re-using the frozen tables.

**Inherited, not re-run.** The `E[kill]`-versus-capacity tables come from `results/sweep/capacity_ranking_robustness.json` — 220 transport solves per panel on the real Xenium reconstructions, at dx = 15 µm, dt = 20 s, convection **off**, 6 h window. This study does not re-solve them. Every ordering statement here is only as good as that record, including its convection-off limitation, and its own worst absolute interpolation error of 3.5e-04 is the relevant numerical bound. **If that record is wrong, sec. 3.2–4 are wrong with it**; sec. 1 and 2 are independent of it.

**Assumed, and attackable.**

* That removing a double-count means stripping the log-linearly predictable component and keeping the level. Sec. 3.1 argues the level cannot have been double-counted; a reader who disagrees would get a *smaller* residual capacity and a more benign answer, so this is the load-bearing assumption and it is stated where it can be attacked.
* That a *linear* fit in log space is the right way to project out the proxy. A monotone but non-linear relationship would leave residual structure this construction does not remove.
* That the structural `K_b` rule is an acceptable stand-in for measured capacity. It is not one, it has been refuted on other data, and sec. 2.1 says so before any number is quoted.

**Not claimable from this document.**

* That `uptake_proxy` agrees, or disagrees, with measured tissue capacity. **n = 1.** Nothing was measured about that here, and no structure was sourced to make it n > 1.
* That the published ranking failure is smaller than reported — except for the single named breast arm in sec. 4.1, under the assumptions listed there.
* That a capacity term should, or should not, be added to the engine. `results/intracellular_retardation.md` refuted the only available per-compound rule and nothing here reverses that.
* Any absolute kill number; anything about `cell_partition` (K), which remains unconstrained by any measurement; anything clinical.

## 6. Reproduction

```
./.venv/bin/python experiments/uptake_proxy_double_count.py --stage all
```

Monte-Carlo seed **20260821**, bootstrap seed **20260822**, 40 000 draws, 4 000 bootstrap replicates — the same seeds and counts as `capacity_ranking_robustness.md`, so the two are directly comparable.

**Nothing in `tmesim/`, `data_cache/`, `docs/`, `app/`, `landing-3d/` or any existing `results/*.md` was modified.** Two files were added: `experiments/uptake_proxy_double_count.py` and this document, plus the frozen record `results/sweep/uptake_proxy_double_count.json`.

**No assertion was deleted or weakened.** `./.venv/bin/python -m pytest -q` → **1171 passed, 0 failed** (407.62 s), the unchanged baseline.

**Rider disclosure — the working tree was not clean.** Sibling agents were live in this repository during the run and the tree carries other untracked files (`experiments/pk_ranking_robustness.py`, `results/axis_inventory.md`). No tracked file under `tmesim/` was modified (`git status --porcelain` shows only untracked additions), and every constant this study leans on is asserted against LIVE code at the start of every stage and matched. Saying the tree was dirty is cheaper than asserting it was not.
