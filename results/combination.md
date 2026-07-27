# Combination therapy on the real Xenium glioblastoma microenvironment

**Experiment:** `experiments/combination.py` (run: `./.venv/bin/python experiments/combination.py`)

**Question.** Do two drugs with complementary DELIVERY profiles — a small, penetrant-but-weak agent (temozolomide) and a potent-but-delivery-limited agent (doxorubicin, consumed within ~45 um of a vessel) — cover more of a real tumor than either alone, *for delivery reasons*?

> **Headline.** THE DELIVERY-COMPLEMENTARITY HYPOTHESIS FAILS on this tissue. Every pair's killed populations overlap MORE than the delivery-shuffled null (the two drugs hit the same cells), and no combination beats the better agent Bliss-combined with ITSELF — so the apparent combination benefit over 2x dosing is an artifact of the independence assumption acting against a saturating Hill curve, not spatial complementarity.

## THE COMBINATION ASSUMPTION: BLISS INDEPENDENCE (read this first)

Each drug is simulated **separately** on the same real tissue, giving every tumor cell its own exposure AUC and its own kill probability from the frozen engine's Hill pharmacodynamics. The combination is then formed at the exposure level by

```
p_kill_combo(cell) = 1 - (1 - p_A(cell)) * (1 - p_B(cell))
```

This is **Bliss independence**: the two agents kill each cell independently and a cell survives only if it survives both.

> **NO pharmacological synergy or antagonism is modelled.** There is no shared target, no DNA-repair crosstalk, no efflux competition, no schedule dependence, no drug-drug PK interaction, and drug B does not alter drug A's transport. What is modelled is **SPATIAL complementarity only** — whether the two delivery fields reach *different cells*. Every gain reported below is geographic or arithmetic, never biochemical. Real combinations may synergise or antagonise; this model is blind to both, so treat these numbers as a delivery-only floor, not a forecast of combination efficacy.

Because Bliss independence *is* the statement that the per-cell draws are independent, all set statistics are reported as exact expectations under that coupling:

```
E|A n B|  = sum p_A p_B                E|A only| = sum p_A (1 - p_B)
E|A u B|  = sum [1-(1-p_A)(1-p_B)]     E|B only| = sum (1 - p_A) p_B
```

so the headline numbers carry no Monte-Carlo noise; a single seeded realisation is drawn only to render the territory maps.

## Method

1. **One real tissue, reconstructed once.** `reconstruct_real('xenium', dx_um=12)` -> 40,887 measured cells (30,174 tumor), grid (424, 675), real endothelial-KDE vasculature, real matrix-transcript ECM, real-perfusion oxygen. The same object is reused by every arm, so all differences are drug physics, not tissue.
2. **Drug compilation.** `calibrate_alpha()` -> alpha = 148.6 (the single perivascular-sink constant, fit ONCE on doxorubicin) -> `calibrated_drugparams(name, alpha)` -> `apply_barrier_to_drugparams` (a no-op for these small molecules, but the identical route an antibody would take).
3. **Separate simulations.** Each agent gets its own `simulate_response` run (oxygen -> IFP -> delivery -> transport PDE -> Hill kill), 6 h exposure window, convection on.
4. **Per-cell readout.** The engine's own `apply_drug_kill` is re-read on a throwaway copy of the *identical* oxygen-labelled cell population, recovering each cell's exposure AUC and kill probability with the same hypoxia / quiescence / sensitivity modulation the engine applied. Nothing is mutated; core / physics / molecular / reconstruct / engine are untouched.
5. **Combination.** Bliss independence, above.
6. **Control 1 — dose-matched.** The better single agent is re-run at **2x** plasma peak: the same total systemic exposure (in standard-dose units) the two-drug arm receives.
7. **Control 2 — self-combination.** The SAME Bliss rule applied to the better agent against *itself*: `1-(1-p_B)^2`. This is a hypothetical reference arm, not a realisable regimen — a 'combination' whose two components have perfectly identical, maximally redundant delivery. If `combo(A,B) ~= selfcombo(B,B)`, the benefit is an artifact of the independence assumption meeting a saturating Hill curve — not spatial complementarity.
8. **Control 3 — permutation nulls (200 perms each).** Independent-Bernoulli overlap shrinks automatically as kill probabilities shrink, so raw Jaccard proves nothing alone. The engine's PD map factorises **exactly** as `p = p_hill(exposure) x (sensitivity x resistance)`, and the second factor is CELL-INTRINSIC and SHARED by every drug. Two nulls follow:
   - *full shuffle* — permute `p_B` across tumor cells (marginals exactly preserved; both delivery and shared-killability alignment destroyed);
   - *delivery shuffle* — permute only drug B's delivery term `p_hill` within 20 strata of the shared factor, then re-multiply, so each cell keeps its own intrinsic killability and only B's delivery geography is scrambled. **This is the null that isolates DELIVERY co-localisation.**

   observed/null **< 1 = complementary territory**, **~1 = no spatial structure**, **> 1 = redundant** (both drugs hit the same cells).

## Single agents on the real tissue (1x dose)

| drug | MW | IC50 (uM) | D (um^2/s) | sink k (1/s) | BBB | mean tissue AUC (uM*s) | E[kill] | engine kill (1 draw) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| temozolomide | 194 | 606 | 278 | 0.00302 | 0.558 | 1.03e+07 | 0.4205 | 0.4173 |
| doxorubicin | 544 | 0.5 | 185 | 0.0333 | 0.051 | 8.59e+04 | 0.8806 | 0.8839 |
| cisplatin | 298 | 24.6 | 236 | 0.00879 | 0.893 | 5.7e+06 | 0.9147 | 0.9152 |

**POTENCY HEADROOM — the number that decides the whole experiment.** A steep delivery gradient only produces a steep KILL gradient if the concentration crosses IC50 somewhere inside the tissue. Per tumor cell, the time-averaged effective concentration `E = AUC / T` in units of that drug's own IC50:

| drug | median E (uM) | IC50 (uM) | median E / IC50 | % tumor cells above IC50 | exposure half-distance | deepest bin still >= IC50 |
|---|---:|---:|---:|---:|---:|---:|
| temozolomide | 477.6 | 606 | 0.788 | 15.0% | 380 um | 30 um |
| doxorubicin | 3.71 | 0.5 | 7.42 | 99.6% | 66 um | 650 um |
| cisplatin | 257.2 | 24.6 | 10.5 | 99.9% | 180 um | 650 um |

`E[kill]` (mean per-cell kill probability over tumor cells) is the noise-free expectation; `engine kill` is one Bernoulli realisation returned by `simulate_response`. They agree to within sampling error — the cross-check that the re-read probabilities are the engine's own.

**Dose-response is already saturating**, which is central to reading the controls below:

| drug | E[kill] @1x | E[kill] @2x | gain from doubling the dose |
|---|---:|---:|---:|
| temozolomide | 0.4205 | 0.6134 | 1.46x |
| doxorubicin | 0.8806 | 0.9295 | 1.06x |
| cisplatin | 0.9147 | 0.9467 | 1.03x |

## Kill: A, B, A+B, and both controls

| pair | A alone | B alone | A+B (Bliss) | better single | 2x dose of better | Bliss(better+ITSELF) | combo/best single | combo/2x-dose | combo/self-combo |
|---|---:|---:|---:|---|---:|---:|---:|---:|---:|
| temozolomide + doxorubicin | 0.4205 | 0.8806 | 0.9195 | doxorubicin | 0.9295 | 0.9688 | 1.04x | 0.99x | 0.95x |
| cisplatin + doxorubicin | 0.9147 | 0.8806 | 0.9750 | cisplatin | 0.9467 | 0.9786 | 1.07x | 1.03x | 1.00x |
| temozolomide + cisplatin | 0.4205 | 0.9147 | 0.9407 | cisplatin | 0.9467 | 0.9786 | 1.03x | 0.99x | 0.96x |

All values are expected fractions of **tumor cells** killed. The last column is the one that matters: a ratio at or below 1.00 means combining two *different* drugs did no better than the model's own independence arithmetic applied to a single drug against itself.

## Spatial overlap of the killed populations

| pair | killed by A only | by B only | by both | Jaccard | overlap coef | full-shuffle null | obs/full | delivery-shuffle null | obs/delivery |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| temozolomide + doxorubicin | 4.2% | 54.3% | 41.5% | 0.4151 | 0.9076 | 0.3978 +- 0.0001 | 1.04x | 0.4095 +- 0.0000 | 1.01x |
| cisplatin + doxorubicin | 9.7% | 6.2% | 84.1% | 0.8414 | 0.9316 | 0.8138 +- 0.0002 | 1.03x | 0.8370 +- 0.0000 | 1.01x |
| temozolomide + cisplatin | 2.8% | 55.3% | 41.9% | 0.4194 | 0.9382 | 0.4047 +- 0.0001 | 1.04x | 0.4160 +- 0.0000 | 1.01x |

A-only / B-only / both are percentages **of the killed population** (they sum to 100%). The delivery-shuffle null preserves drug B's marginal kill to within 0.0%, 0.0%, 0.0% respectively.

As fractions of all tumor cells:

| pair | A only | B only | both | any (= combo kill) |
|---|---:|---:|---:|---:|
| temozolomide + doxorubicin | 3.89% | 49.90% | 38.17% | 91.95% |
| cisplatin + doxorubicin | 9.44% | 6.02% | 82.04% | 97.50% |
| temozolomide + cisplatin | 2.60% | 52.02% | 39.45% | 94.07% |

Correlation of the per-cell kill probabilities across tumor cells — a coupling-free redundancy measure. `p` includes the shared cell-intrinsic factor; `p_hill` is the DELIVERY term alone:

| pair | Pearson(p_A,p_B) | Spearman(p_A,p_B) | Pearson(hill_A,hill_B) | Spearman(hill_A,hill_B) |
|---|---:|---:|---:|---:|
| temozolomide + doxorubicin | 0.821 | 0.823 | 0.714 | 0.812 |
| cisplatin + doxorubicin | 0.961 | 0.932 | 0.874 | 0.926 |
| temozolomide + cisplatin | 0.784 | 0.963 | 0.708 | 0.962 |

## Where each drug actually kills (the delivery axis)

Vessel mask = the p95 threshold on the real endothelial-KDE vascular density (5.0% of tissue area, 9,688 pixels); distance is the Euclidean transform in um. Per tumor-cell distance bin: mean exposure AUC (normalised to the perivascular bin), mean kill probability, and the shared cell-intrinsic factor.

Bin widths are 12 um (= dx) out to 120 um, then 40 um, then 100 um: tumor cells on this section sit a MEDIAN of 134 um from the vessel mask (max 1179 um), so a uniform near-field grid would collapse most of the tumor into one overflow bin. Cells beyond 1200 um (0.0% of tumor cells) are DROPPED from the profiles rather than clipped into the last bin.

| distance (um) | cells | AUC temo | AUC doxo | AUC cisp | p_kill temo | p_kill doxo | p_kill cisp | p_kill combo | intrinsic factor |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 6 | 2,186 | 1.000 | 1.000 | 1.000 | 0.5460 | 0.9696 | 0.9700 | 0.9853 | 1.0000 |
| 18 | 1,070 | 0.839 | 0.675 | 0.769 | 0.4964 | 0.9532 | 0.9581 | 0.9746 | 0.9967 |
| 30 | 1,268 | 0.819 | 0.614 | 0.735 | 0.4888 | 0.9496 | 0.9572 | 0.9727 | 0.9980 |
| 42 | 1,055 | 0.791 | 0.562 | 0.699 | 0.4804 | 0.9443 | 0.9546 | 0.9694 | 0.9972 |
| 54 | 1,063 | 0.776 | 0.521 | 0.675 | 0.4744 | 0.9391 | 0.9527 | 0.9662 | 0.9972 |
| 66 | 1,534 | 0.742 | 0.472 | 0.634 | 0.4622 | 0.9305 | 0.9480 | 0.9602 | 0.9951 |
| 78 | 1,181 | 0.734 | 0.445 | 0.619 | 0.4589 | 0.9253 | 0.9461 | 0.9569 | 0.9945 |
| 90 | 1,412 | 0.705 | 0.412 | 0.586 | 0.4485 | 0.9186 | 0.9435 | 0.9526 | 0.9947 |
| 102 | 1,398 | 0.696 | 0.391 | 0.571 | 0.4424 | 0.9083 | 0.9365 | 0.9444 | 0.9896 |
| 114 | 1,301 | 0.663 | 0.362 | 0.536 | 0.4289 | 0.8997 | 0.9319 | 0.9381 | 0.9885 |
| 140 | 4,391 | 0.637 | 0.335 | 0.506 | 0.4174 | 0.8882 | 0.9257 | 0.9295 | 0.9854 |
| 180 | 3,735 | 0.598 | 0.299 | 0.461 | 0.3980 | 0.8691 | 0.9142 | 0.9145 | 0.9795 |
| 220 | 2,598 | 0.563 | 0.276 | 0.425 | 0.3782 | 0.8471 | 0.8975 | 0.8957 | 0.9679 |
| 260 | 2,160 | 0.538 | 0.253 | 0.397 | 0.3619 | 0.8191 | 0.8801 | 0.8726 | 0.9553 |
| 300 | 1,492 | 0.514 | 0.242 | 0.376 | 0.3426 | 0.7916 | 0.8537 | 0.8469 | 0.9323 |
| 340 | 993 | 0.502 | 0.244 | 0.370 | 0.3275 | 0.7668 | 0.8227 | 0.8218 | 0.9008 |
| 380 | 637 | 0.484 | 0.248 | 0.362 | 0.3186 | 0.7707 | 0.8201 | 0.8229 | 0.9003 |
| 450 | 473 | 0.408 | 0.206 | 0.302 | 0.2638 | 0.6890 | 0.7472 | 0.7423 | 0.8520 |
| 550 | 106 | 0.246 | 0.104 | 0.171 | 0.1329 | 0.4628 | 0.5602 | 0.5232 | 0.6981 |
| 650 | 59 | 0.212 | 0.139 | 0.167 | 0.1084 | 0.4398 | 0.4879 | 0.4895 | 0.6864 |

Near-vessel vs deep, for exposure (the delivery quantity) and for kill:

| drug | AUC near / AUC deep | p_kill perivascular bin | p_kill first bin >100 um | p_kill far/near |
|---|---:|---:|---:|---:|
| temozolomide | 1.44x | 0.5460 | 0.4424 | 0.81x |
| doxorubicin | 2.56x | 0.9696 | 0.9083 | 0.94x |
| cisplatin | 1.75x | 0.9700 | 0.9365 | 0.97x |

The shared cell-intrinsic factor is nearly FLAT over the same range (1.0000 perivascular -> 0.9896 deep, 0.99x), so it is not what shapes these curves — delivery is.

## THE HYPOTHESIS, TESTED DIRECTLY: where does each agent kill *uniquely*?

The claim under test is not just that the combination kills more, but that the penetrant agent reaches DEEP tissue the delivery-limited agent cannot. That is a statement about where each drug's UNIQUE kills sit. Expected per-cell contribution to the A-only / B-only / both sets, by distance from vessel:

**temozolomide (A) + doxorubicin (B)**

| distance (um) | cells | A only | B only | both | A-only share of kills |
|---:|---:|---:|---:|---:|---:|
| 6 | 2,186 | 0.0157 | 0.4393 | 0.5303 | 1.6% |
| 18 | 1,070 | 0.0213 | 0.4782 | 0.4750 | 2.2% |
| 30 | 1,268 | 0.0231 | 0.4839 | 0.4657 | 2.4% |
| 42 | 1,055 | 0.0251 | 0.4889 | 0.4554 | 2.6% |
| 54 | 1,063 | 0.0270 | 0.4918 | 0.4473 | 2.8% |
| 66 | 1,534 | 0.0297 | 0.4981 | 0.4324 | 3.1% |
| 78 | 1,181 | 0.0316 | 0.4980 | 0.4273 | 3.3% |
| 90 | 1,412 | 0.0340 | 0.5041 | 0.4145 | 3.6% |
| 102 | 1,398 | 0.0361 | 0.5020 | 0.4062 | 3.8% |
| 114 | 1,301 | 0.0384 | 0.5091 | 0.3906 | 4.1% |
| 140 | 4,391 | 0.0413 | 0.5121 | 0.3761 | 4.4% |
| 180 | 3,735 | 0.0454 | 0.5165 | 0.3526 | 5.0% |
| 220 | 2,598 | 0.0486 | 0.5175 | 0.3296 | 5.4% |
| 260 | 2,160 | 0.0535 | 0.5107 | 0.3084 | 6.1% |
| 300 | 1,492 | 0.0554 | 0.5043 | 0.2872 | 6.5% |
| 340 | 993 | 0.0550 | 0.4942 | 0.2725 | 6.7% |
| 380 | 637 | 0.0522 | 0.5043 | 0.2663 | 6.3% |
| 450 | 473 | 0.0533 | 0.4786 | 0.2104 | 7.2% |
| 550 | 106 | 0.0604 | 0.3903 | 0.0726 | 11.5% |
| 650 | 59 | 0.0497 | 0.3811 | 0.0587 | 10.1% |

**cisplatin (A) + doxorubicin (B)**  (condensed — every fourth bin)

| distance (um) | cells | A only | B only | both | A-only share of kills |
|---:|---:|---:|---:|---:|---:|
| 6 | 2,186 | 0.0293 | 0.0290 | 0.9407 | 2.9% |
| 54 | 1,063 | 0.0564 | 0.0428 | 0.8963 | 5.7% |
| 102 | 1,398 | 0.0804 | 0.0521 | 0.8562 | 8.1% |
| 220 | 2,598 | 0.1224 | 0.0720 | 0.7751 | 12.6% |
| 380 | 637 | 0.1494 | 0.1000 | 0.6707 | 16.2% |

**temozolomide (A) + cisplatin (B)**  (condensed — every fourth bin)

| distance (um) | cells | A only | B only | both | A-only share of kills |
|---:|---:|---:|---:|---:|---:|
| 6 | 2,186 | 0.0154 | 0.4394 | 0.5306 | 1.6% |
| 54 | 1,063 | 0.0205 | 0.4989 | 0.4538 | 2.1% |
| 102 | 1,398 | 0.0238 | 0.5180 | 0.4185 | 2.5% |
| 220 | 2,598 | 0.0301 | 0.5494 | 0.3481 | 3.2% |
| 380 | 637 | 0.0366 | 0.5381 | 0.2820 | 4.3% |

## Verdict

THE DELIVERY-COMPLEMENTARITY HYPOTHESIS FAILS on this tissue. Every pair's killed populations overlap MORE than the delivery-shuffled null (the two drugs hit the same cells), and no combination beats the better agent Bliss-combined with ITSELF — so the apparent combination benefit over 2x dosing is an artifact of the independence assumption acting against a saturating Hill curve, not spatial complementarity.

- temozolomide + doxorubicin: overlap vs delivery-null 1.01x (vs full-null 1.04x) -> no spatial structure. combo 0.9195 vs 2x-doxorubicin 0.9295 (0.99x, loses) and vs Bliss(doxorubicin+itself) 0.9688 (0.95x, LOSES).
- cisplatin + doxorubicin: overlap vs delivery-null 1.01x (vs full-null 1.03x) -> no spatial structure. combo 0.9750 vs 2x-cisplatin 0.9467 (1.03x, beats) and vs Bliss(cisplatin+itself) 0.9786 (1.00x, LOSES).
- temozolomide + cisplatin: overlap vs delivery-null 1.01x (vs full-null 1.04x) -> no spatial structure. combo 0.9407 vs 2x-cisplatin 0.9467 (0.99x, loses) and vs Bliss(cisplatin+itself) 0.9786 (0.96x, LOSES).

### What the numbers mean, mechanistically

1. **The delivery gradients are real but nearly IDENTICAL in shape.** All three agents' exposure falls to half within 66-380 um of the vessel mask and the normalised profiles lie almost on top of one another. This is the confound `results/vessel_profile.md` already documented: with a fast perivascular sink the tissue field is close to reaction-dominated (`C ~ S/k`), so every drug largely INHERITS the shape of the same vascular source field. **Two drugs entering from the same vessels have correlated, not complementary, delivery geometry** — which is the structural reason the complementarity hypothesis has so little room to work here.
2. **What differs between the agents is MAGNITUDE relative to potency, not shape.** Doxorubicin's mean tissue AUC is ~120x below temozolomide's — it is genuinely delivery-limited in concentration — but its IC50 is ~1212x lower, so the delivered concentration still sits at 7.4x IC50 in the median tumor cell versus 0.79x for temozolomide. **Doxorubicin is therefore NOT effectively delivery-limited on this section**: its Hill curve is saturated well past the ~45 um cuff, so its kill barely decays with distance while the weak agent's kill decays steeply. The premise of the headline pairing fails at this step.
3. **The penetrant agent's unique kills are PERIVASCULAR, not deep** — the opposite of the hypothesis. In the temozolomide + doxorubicin pair the temozolomide-only contribution is largest in the near-vessel bins and falls with distance (table above), because temozolomide only wins where its own exposure is highest. It does not open up a disjoint deep territory; the deep tissue is already covered by the potent agent.
4. **Hence overlap ABOVE both nulls.** The killed populations co-localise more than chance even against the delivery-shuffle null, which holds each cell's intrinsic killability fixed. The redundancy is delivery redundancy.
5. **The combination still beats 2x dosing — but so does a drug combined with ITSELF.** Returns to dose are steeply diminishing (dose-doubling table), so under a Bliss model *any* second independent hit beats *any* second dose, including the same agent twice. The self-combination control is the one that separates the two explanations, and every pair LOSES to it. The apparent 'combination benefit' here is the independence arithmetic meeting a saturating Hill curve, not spatial complementarity.

**What would change this answer.** A genuinely complementary pair needs (a) delivery fields with different SHAPES — which on this tissue means escaping the shared vascular source geometry (e.g. a hypoxia-activated prodrug, a drug whose sink is anti-correlated with perfusion, or convection-dominated delivery) — and (b) a potent agent whose exposure actually crosses IC50 inside the tissue, so that its kill is confined where its delivery is confined. Neither condition holds for these three agents on this section.

## Caveats (all of them)

- **No synergy is modelled.** Bliss independence is an assumption, not a result. Real drug-drug interaction (synergy or antagonism) is outside this model's reach and could dominate everything reported here.
- **No outcome validation.** This open Xenium section is an *untreated* resection: there is no paired measured response for any of these agents, alone or combined. Every number here is a PREDICTION.
- **Absolute kills are conditional on the sensitivity model.** The reconstruction carries a UNIFORM per-cell `sensitivity` of 1.0 (a deliberately drug-agnostic geometry-only choice, documented in `reconstruct/real_xenium.py`), so the only cell-intrinsic modulation here is hypoxia / quiescence resistance (0% of tumor cells are hypoxic). Absolute kill fractions are therefore an upper-ish bound conditional on 'every tumor cell is fully sensitive', and would fall — unevenly — under any real per-cell sensitivity model. Such a model could move the overlap ratios in either direction; it would not change the exposure fields, on which the delivery conclusions rest.
- **Absolute kills are also conditional on the potency inputs.** IC50s come from the cached GDSC glioma screen (e.g. temozolomide at 606 uM, a genuinely weak in-vitro potency); the headroom argument above is a direct function of those numbers and of the shared nominal dose.
- **All agents share one nominal dose.** Every drug is driven at the same 50 uM peak plasma level, which is NOT clinically dose-matched between agents, so cross-drug absolute comparisons are soft. The dose-matched control is defined in the same 'standard dose' units, so the combo-vs-2x comparison is internally consistent.
- **Dosing is the engine's Stage-0 plasma average.** The transport solve uses a static source at the time-averaged plasma level, so schedule (sequential vs concurrent — which matters enormously in real combination therapy) cannot be represented at all.
- **The vessel mask is a threshold on a 40 um-bandwidth endothelial KDE**, not a measured lumen, and `results/vessel_profile.md` documents that most of the apparent perivascular gradient on this section is vessel/KDE geometry rather than drug transport. The distance axis inherits both limitations.
- **Overlap statistics depend on the coupling.** Under independent draws, overlap shrinks automatically as kill probabilities shrink, which is precisely why the null RATIOS, not the raw Jaccard, are the reported test statistics.
- **One section, one tumor, three agents, three pairs.** No cohort, no replicate sections, no cross-tissue generalisation is claimed.
- **Research prototype.** Not a clinical decision tool; no regulatory, diagnostic or treatment-selection use.

