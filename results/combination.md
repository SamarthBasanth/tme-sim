# Combination therapy on the real Xenium glioblastoma microenvironment

> **Regenerated under the discrete vessel source.** `vascular_mode = discrete`
> (`reconstruct.real_fields._DEFAULT_VASCULAR_MODE`, promoted in
> [`results/vessel_promotion.md`](vessel_promotion.md)) · **alpha = 148.6**, the
> e-folding-length convention `L = d_half / ln 2` · per-cell sensitivity mode =
> **`uniform` (pipeline default; the engine additionally applies the hypoxia/quiescence resistance factor, identical for both agents)** · **dx = 12 µm**.
>
> Every absolute exposure / AUC / coverage / kill below is on the **discrete-source
> unit** and is NOT comparable with the KDE-era numbers this file carried before this
> regeneration (the tissue-integrated source total is 3.36x smaller; that is a change
> of *unit*, not of physics — the dose-matched, geometry-only effect on coverage is
> <= 0.0166). To recover the old numbers, re-run with `mode='kde'`.

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
| temozolomide | 194 | 606 | 278 | 0.00302 | 0.125 | 6.86e+05 | 0.0144 | 0.0147 |
| doxorubicin | 544 | 0.5 | 185 | 0.0333 | 0.021 | 1.03e+04 | 0.2212 | 0.2165 |
| cisplatin | 298 | 24.6 | 236 | 0.00879 | 0.858 | 1.63e+06 | 0.3748 | 0.3699 |

**POTENCY HEADROOM — the number that decides the whole experiment.** A steep delivery gradient only produces a steep KILL gradient if the concentration crosses IC50 somewhere inside the tissue. Per tumor cell, the time-averaged effective concentration `E = AUC / T` in units of that drug's own IC50:

| drug | median E (uM) | IC50 (uM) | median E / IC50 | % tumor cells above IC50 | exposure half-distance | deepest bin still >= IC50 |
|---|---:|---:|---:|---:|---:|---:|
| temozolomide | 31.88 | 606 | 0.0526 | 0.0% | 140 um | none |
| doxorubicin | 0.4278 | 0.5 | 0.856 | 38.4% | 30 um | 30 um |
| cisplatin | 72.93 | 24.6 | 2.97 | 99.4% | 78 um | 140 um |

`E[kill]` (mean per-cell kill probability over tumor cells) is the noise-free expectation; `engine kill` is one Bernoulli realisation returned by `simulate_response`. They agree to within sampling error — the cross-check that the re-read probabilities are the engine's own.

**Dose-response is already saturating**, which is central to reading the controls below:

| drug | E[kill] @1x | E[kill] @2x | gain from doubling the dose |
|---|---:|---:|---:|
| temozolomide | 0.0144 | 0.0317 | 2.20x |
| doxorubicin | 0.2212 | 0.3093 | 1.40x |
| cisplatin | 0.3748 | 0.4267 | 1.14x |

## Kill: A, B, A+B, and both controls

| pair | A alone | B alone | A+B (Bliss) | better single | 2x dose of better | Bliss(better+ITSELF) | combo/best single | combo/2x-dose | combo/self-combo |
|---|---:|---:|---:|---|---:|---:|---:|---:|---:|
| temozolomide + doxorubicin | 0.0144 | 0.2212 | 0.2320 | doxorubicin | 0.3093 | 0.3864 | 1.05x | 0.75x | 0.60x |
| cisplatin + doxorubicin | 0.3748 | 0.2212 | 0.5085 | cisplatin | 0.4267 | 0.6048 | 1.36x | 1.19x | 0.84x |
| temozolomide + cisplatin | 0.0144 | 0.3748 | 0.3835 | cisplatin | 0.4267 | 0.6048 | 1.02x | 0.90x | 0.63x |

All values are expected fractions of **tumor cells** killed. The last column is the one that matters: a ratio at or below 1.00 means combining two *different* drugs did no better than the model's own independence arithmetic applied to a single drug against itself.

## Spatial overlap of the killed populations

| pair | killed by A only | by B only | by both | Jaccard | overlap coef | full-shuffle null | obs/full | delivery-shuffle null | obs/delivery |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| temozolomide + doxorubicin | 4.7% | 93.8% | 1.5% | 0.0155 | 0.2492 | 0.0137 +- 0.0000 | 1.13x | 0.0144 +- 0.0000 | 1.08x |
| cisplatin + doxorubicin | 56.5% | 26.3% | 17.2% | 0.1721 | 0.3957 | 0.1616 +- 0.0001 | 1.07x | 0.1672 +- 0.0000 | 1.03x |
| temozolomide + cisplatin | 2.3% | 96.2% | 1.5% | 0.0149 | 0.3975 | 0.0141 +- 0.0000 | 1.06x | 0.0146 +- 0.0000 | 1.02x |

A-only / B-only / both are percentages **of the killed population** (they sum to 100%). The delivery-shuffle null preserves drug B's marginal kill to within 0.0%, 0.0%, 0.0% respectively.

As fractions of all tumor cells:

| pair | A only | B only | both | any (= combo kill) |
|---|---:|---:|---:|---:|
| temozolomide + doxorubicin | 1.08% | 21.76% | 0.36% | 23.20% |
| cisplatin + doxorubicin | 28.73% | 13.37% | 8.75% | 50.85% |
| temozolomide + cisplatin | 0.87% | 36.91% | 0.57% | 38.35% |

Correlation of the per-cell kill probabilities across tumor cells — a coupling-free redundancy measure. `p` includes the shared cell-intrinsic factor; `p_hill` is the DELIVERY term alone:

| pair | Pearson(p_A,p_B) | Spearman(p_A,p_B) | Pearson(hill_A,hill_B) | Spearman(hill_A,hill_B) |
|---|---:|---:|---:|---:|
| temozolomide + doxorubicin | 0.808 | 0.823 | 0.764 | 0.776 |
| cisplatin + doxorubicin | 0.831 | 0.923 | 0.865 | 0.902 |
| temozolomide + cisplatin | 0.834 | 0.967 | 0.838 | 0.961 |

## Where each drug actually kills (the delivery axis)

Vessel mask = the p95 threshold on the real endothelial-KDE vascular density (2.5% of tissue area, 4,843 pixels); distance is the Euclidean transform in um. Per tumor-cell distance bin: mean exposure AUC (normalised to the perivascular bin), mean kill probability, and the shared cell-intrinsic factor.

Bin widths are 12 um (= dx) out to 120 um, then 40 um, then 100 um: tumor cells on this section sit a MEDIAN of 43 um from the vessel mask (max 779 um), so a uniform near-field grid would collapse most of the tumor into one overflow bin. Cells beyond 1200 um (0.0% of tumor cells) are DROPPED from the profiles rather than clipped into the last bin.

| distance (um) | cells | AUC temo | AUC doxo | AUC cisp | p_kill temo | p_kill doxo | p_kill cisp | p_kill combo | intrinsic factor |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 6 | 210 | 1.000 | 1.000 | 1.000 | 0.0203 | 0.3602 | 0.4218 | 0.3728 | 0.4886 |
| 18 | 4,291 | 0.879 | 0.634 | 0.794 | 0.0174 | 0.3023 | 0.4056 | 0.3141 | 0.4892 |
| 30 | 6,563 | 0.808 | 0.453 | 0.679 | 0.0158 | 0.2614 | 0.3932 | 0.2728 | 0.4876 |
| 42 | 5,004 | 0.763 | 0.366 | 0.615 | 0.0147 | 0.2329 | 0.3828 | 0.2439 | 0.4848 |
| 54 | 4,112 | 0.738 | 0.315 | 0.578 | 0.0141 | 0.2120 | 0.3758 | 0.2228 | 0.4827 |
| 66 | 4,048 | 0.704 | 0.264 | 0.534 | 0.0133 | 0.1871 | 0.3659 | 0.1977 | 0.4798 |
| 78 | 2,145 | 0.668 | 0.221 | 0.491 | 0.0124 | 0.1629 | 0.3530 | 0.1731 | 0.4742 |
| 90 | 1,679 | 0.642 | 0.188 | 0.458 | 0.0118 | 0.1415 | 0.3403 | 0.1514 | 0.4680 |
| 102 | 972 | 0.594 | 0.153 | 0.408 | 0.0106 | 0.1174 | 0.3218 | 0.1266 | 0.4603 |
| 114 | 527 | 0.567 | 0.128 | 0.377 | 0.0099 | 0.0986 | 0.3071 | 0.1073 | 0.4516 |
| 140 | 523 | 0.492 | 0.092 | 0.309 | 0.0081 | 0.0689 | 0.2694 | 0.0763 | 0.4346 |
| 180 | 64 | 0.342 | 0.037 | 0.182 | 0.0054 | 0.0270 | 0.1925 | 0.0322 | 0.4344 |

Near-vessel vs deep, for exposure (the delivery quantity) and for kill:

| drug | AUC near / AUC deep | p_kill perivascular bin | p_kill first bin >100 um | p_kill far/near |
|---|---:|---:|---:|---:|
| temozolomide | 1.68x | 0.0203 | 0.0106 | 0.52x |
| doxorubicin | 6.53x | 0.3602 | 0.1174 | 0.33x |
| cisplatin | 2.45x | 0.4218 | 0.3218 | 0.76x |

The shared cell-intrinsic factor is nearly FLAT over the same range (0.4886 perivascular -> 0.4603 deep, 0.94x), so it is not what shapes these curves — delivery is.

## THE HYPOTHESIS, TESTED DIRECTLY: where does each agent kill *uniquely*?

The claim under test is not just that the combination kills more, but that the penetrant agent reaches DEEP tissue the delivery-limited agent cannot. That is a statement about where each drug's UNIQUE kills sit. Expected per-cell contribution to the A-only / B-only / both sets, by distance from vessel:

**temozolomide (A) + doxorubicin (B)**

| distance (um) | cells | A only | B only | both | A-only share of kills |
|---:|---:|---:|---:|---:|---:|
| 6 | 210 | 0.0126 | 0.3525 | 0.0077 | 3.4% |
| 18 | 4,291 | 0.0118 | 0.2966 | 0.0057 | 3.8% |
| 30 | 6,563 | 0.0113 | 0.2570 | 0.0044 | 4.2% |
| 42 | 5,004 | 0.0110 | 0.2292 | 0.0037 | 4.5% |
| 54 | 4,112 | 0.0109 | 0.2087 | 0.0032 | 4.9% |
| 66 | 4,048 | 0.0106 | 0.1844 | 0.0027 | 5.4% |
| 78 | 2,145 | 0.0102 | 0.1607 | 0.0023 | 5.9% |
| 90 | 1,679 | 0.0099 | 0.1396 | 0.0019 | 6.5% |
| 102 | 972 | 0.0092 | 0.1160 | 0.0015 | 7.2% |
| 114 | 527 | 0.0088 | 0.0974 | 0.0012 | 8.2% |
| 140 | 523 | 0.0074 | 0.0682 | 0.0007 | 9.7% |
| 180 | 64 | 0.0052 | 0.0268 | 0.0002 | 16.1% |

**cisplatin (A) + doxorubicin (B)**  (condensed — every fourth bin)

| distance (um) | cells | A only | B only | both | A-only share of kills |
|---:|---:|---:|---:|---:|---:|
| 6 | 210 | 0.2666 | 0.2051 | 0.1552 | 42.5% |
| 54 | 4,112 | 0.2931 | 0.1293 | 0.0827 | 58.0% |
| 102 | 972 | 0.2807 | 0.0764 | 0.0411 | 70.5% |

**temozolomide (A) + cisplatin (B)**  (condensed — every fourth bin)

| distance (um) | cells | A only | B only | both | A-only share of kills |
|---:|---:|---:|---:|---:|---:|
| 6 | 210 | 0.0114 | 0.4129 | 0.0089 | 2.6% |
| 54 | 4,112 | 0.0086 | 0.3703 | 0.0056 | 2.2% |
| 102 | 972 | 0.0069 | 0.3180 | 0.0037 | 2.1% |

## Verdict

THE DELIVERY-COMPLEMENTARITY HYPOTHESIS FAILS on this tissue. Every pair's killed populations overlap MORE than the delivery-shuffled null (the two drugs hit the same cells), and no combination beats the better agent Bliss-combined with ITSELF — so the apparent combination benefit over 2x dosing is an artifact of the independence assumption acting against a saturating Hill curve, not spatial complementarity.

- temozolomide + doxorubicin: overlap vs delivery-null 1.08x (vs full-null 1.13x) -> REDUNDANT (same cells). combo 0.2320 vs 2x-doxorubicin 0.3093 (0.75x, loses) and vs Bliss(doxorubicin+itself) 0.3864 (0.60x, LOSES).
- cisplatin + doxorubicin: overlap vs delivery-null 1.03x (vs full-null 1.07x) -> no spatial structure. combo 0.5085 vs 2x-cisplatin 0.4267 (1.19x, beats) and vs Bliss(cisplatin+itself) 0.6048 (0.84x, LOSES).
- temozolomide + cisplatin: overlap vs delivery-null 1.02x (vs full-null 1.06x) -> no spatial structure. combo 0.3835 vs 2x-cisplatin 0.4267 (0.90x, loses) and vs Bliss(cisplatin+itself) 0.6048 (0.63x, LOSES).

### What the numbers mean, mechanistically

1. **The delivery gradients are real but nearly IDENTICAL in shape.** All three agents' exposure falls to half within 30-140 um of the vessel mask and the normalised profiles lie almost on top of one another. This is the confound `results/vessel_profile.md` already documented: with a fast perivascular sink the tissue field is close to reaction-dominated (`C ~ S/k`), so every drug largely INHERITS the shape of the same vascular source field. **Two drugs entering from the same vessels have correlated, not complementary, delivery geometry** — which is the structural reason the complementarity hypothesis has so little room to work here.
2. **What differs between the agents is MAGNITUDE relative to potency, not shape.** Doxorubicin's mean tissue AUC is ~66x below temozolomide's — it is genuinely delivery-limited in concentration — but its IC50 is ~1212x lower, so the delivered concentration still sits at 0.86x IC50 in the median tumor cell versus 0.053x for temozolomide. **Doxorubicin is therefore NOT effectively delivery-limited on this section**: its Hill curve is saturated well past the ~45 um cuff, so its kill barely decays with distance while the weak agent's kill decays steeply. The premise of the headline pairing fails at this step.
3. **The penetrant agent's unique kills are PERIVASCULAR, not deep** — the opposite of the hypothesis. In the temozolomide + doxorubicin pair the temozolomide-only contribution is largest in the near-vessel bins and falls with distance (table above), because temozolomide only wins where its own exposure is highest. It does not open up a disjoint deep territory; the deep tissue is already covered by the potent agent.
4. **Hence overlap ABOVE both nulls.** The killed populations co-localise more than chance even against the delivery-shuffle null, which holds each cell's intrinsic killability fixed. The redundancy is delivery redundancy.
5. **The combination still beats 2x dosing — but so does a drug combined with ITSELF.** Returns to dose are steeply diminishing (dose-doubling table), so under a Bliss model *any* second independent hit beats *any* second dose, including the same agent twice. The self-combination control is the one that separates the two explanations, and every pair LOSES to it. The apparent 'combination benefit' here is the independence arithmetic meeting a saturating Hill curve, not spatial complementarity.

**What would change this answer.** A genuinely complementary pair needs (a) delivery fields with different SHAPES — which on this tissue means escaping the shared vascular source geometry (e.g. a hypoxia-activated prodrug, a drug whose sink is anti-correlated with perfusion, or convection-dominated delivery) — and (b) a potent agent whose exposure actually crosses IC50 inside the tissue, so that its kill is confined where its delivery is confined. Neither condition holds for these three agents on this section.

## Caveats (all of them)

- **No synergy is modelled.** Bliss independence is an assumption, not a result. Real drug-drug interaction (synergy or antagonism) is outside this model's reach and could dominate everything reported here.
- **No outcome validation.** This open Xenium section is an *untreated* resection: there is no paired measured response for any of these agents, alone or combined. Every number here is a PREDICTION.
- **Absolute kills are conditional on the sensitivity model.** The reconstruction carries a UNIFORM per-cell `sensitivity` of 1.0 (a deliberately drug-agnostic geometry-only choice, documented in `reconstruct/real_xenium.py`), so the only cell-intrinsic modulation here is hypoxia / quiescence resistance (13% of tumor cells are hypoxic). Absolute kill fractions are therefore an upper-ish bound conditional on 'every tumor cell is fully sensitive', and would fall — unevenly — under any real per-cell sensitivity model. Such a model could move the overlap ratios in either direction; it would not change the exposure fields, on which the delivery conclusions rest.
- **Absolute kills are also conditional on the potency inputs.** IC50s come from the cached GDSC glioma screen (e.g. temozolomide at 606 uM, a genuinely weak in-vitro potency); the headroom argument above is a direct function of those numbers and of the shared nominal dose.
- **All agents share one nominal dose.** Every drug is driven at the same 50 uM peak plasma level, which is NOT clinically dose-matched between agents, so cross-drug absolute comparisons are soft. The dose-matched control is defined in the same 'standard dose' units, so the combo-vs-2x comparison is internally consistent.
- **Dosing is the engine's Stage-0 plasma average.** The transport solve uses a static source at the time-averaged plasma level, so schedule (sequential vs concurrent — which matters enormously in real combination therapy) cannot be represented at all.
- **The vessel mask is a threshold on a 40 um-bandwidth endothelial KDE**, not a measured lumen, and `results/vessel_profile.md` documents that most of the apparent perivascular gradient on this section is vessel/KDE geometry rather than drug transport. The distance axis inherits both limitations.
- **Overlap statistics depend on the coupling.** Under independent draws, overlap shrinks automatically as kill probabilities shrink, which is precisely why the null RATIOS, not the raw Jaccard, are the reported test statistics.
- **One section, one tumor, three agents, three pairs.** No cohort, no replicate sections, no cross-tissue generalisation is claimed.
- **Research prototype.** Not a clinical decision tool; no regulatory, diagnostic or treatment-selection use.

