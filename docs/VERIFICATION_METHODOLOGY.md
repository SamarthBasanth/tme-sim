# Verification methodology for scientific simulation codebases

**What this document is.** An adversarial verification pass over this project found **eleven distinct
defects** in a simulator whose test suite was green throughout, and retired two of its components
outright. The defects are specific to tumour-microenvironment modelling. The *patterns* are not.

This is a write-up of the patterns, each stated so it can be applied to a computational model in any
field, and each illustrated by the real case from this project with the numbers it actually produced.
Every number below is traceable to a file in `results/` or `docs/`; the index is at the end.

**What it is not.** It is not a claim that this methodology finds all defects, that it is cheap, or
that it was executed correctly here. Section 12 documents what it cost and the analyses inside the
verification pass that were themselves wrong and had to be retracted — including the top-level
diagnosis that motivated the single most expensive experiment in the project.

---

## 0. The ledger

The pass ran against a model with a passing test suite (**55 tests green at the time**; the suite now
collects 114). Test count was not the constraint. Every defect below was found by a check that a unit
test cannot express.

| # | Defect | Class | How it was found |
|---|---|---|---|
| 1 | Endothelial fraction 2× too high (single-marker argmax bias): 20.6% → 9.7% | wrong value | external comparison to volume-EM tissue counting |
| 2 | `penetration_index` not grid-converged **and definitionally unable to discriminate** | metric below its own noise floor | signal-vs-drift ratio (pattern 5) |
| 3 | Hypoxic fraction a grid artifact: 94% at dx=30 → 13% at dx=12 → 7% at dx=8 | invariance violation | grid-refinement invariance check |
| 4 | α off by exactly 1/ln²2 — half-distance vs e-fold convention | unit/convention error | writing the derivation down |
| 5 | Per-cell sensitivity a placeholder suppressing all kills ~5× | placeholder shipped as model | provenance audit |
| 6 | Grid-dependent oxygen sink (raw count, scales as dx²) | invariance violation | dimensional analysis |
| 7 | Grid-dependent tissue mask — **was partially cancelling #6** | compensating error | invariance check *after* fixing #6 |
| 8 | Unfair cross-drug comparison: TMZ advantaged because its determinant is *unmeasured* | measurement-coverage confound | asking what differs between arms besides the thing under test |
| 9 | `bbb_penetrance` arithmetic mean cannot express a veto | aggregation-form error | external validation outlier (66×) traced to source |
| 10 | `membrane_permeability` returned 4.2e7 µm/s for imatinib (42 m/s) | missing physical ceiling | plausibility check on units |
| 11 | Oxygen PDE **falsified as value-add** — beaten by a kd-tree query | trivial-baseline failure | baseline comparison (pattern 2) |

Two of the eleven were not repairs but retirements: `penetration_index` (#2) was withdrawn as a
delivery metric, and the oxygen reaction–diffusion solve (#11) was shown to be out-predicted by a
one-line distance transform computed from its own inputs.

**The generalisable claim.** Nine of these eleven are invisible to unit tests by construction,
because a unit test asserts that a function returns what the author believes it should return. The
checks that found them assert something different: that a result is *invariant* to things it should
not depend on, that it *beats* something cheaper, and that it is *not reachable* by the pipeline's
own inputs by a shorter route.

---

## 1. Compensating errors

> **The pattern.** Two independent defects that push a result in opposite directions make the result
> look *more* stable than either error alone would allow. Fixing one of them makes the output look
> worse. The test suite is green the whole time. This is the most dangerous defect class in
> scientific code, because the usual signal that something is wrong — the answer moving when it
> should not — has been suppressed.

### The case

Two grid-dependence bugs coexisted in the oxygen model.

**Defect 6 — the sink.** The tumour-cell density feeding the oxygen consumption term was a *raw count
per grid cell* rather than a density per unit area. Counts scale as Δx², so the consumption strength
was proportional to Δx². On the real Xenium reconstruction the hypoxic tumour fraction ran
**0.937 / 0.364 / 0.134** at Δx = 30 / 20 / 12 µm; on a controlled synthetic tissue with identical
physical extent and identical cells it ran **0.998 → 0.000** over Δx = 30 → 10 µm.

**Defect 7 — the mask.** `Domain.tissue_mask` was defined as "grid cells containing at least one
segmented cell". Coverage of the bounding box was **50.6% / 32.7% / 13.9%** at the same three
spacings — a **3.64× spread**. Because the vascular field is masked to tissue, the total oxygen
*supply* shrank as the mesh was refined.

The sink error inflated hypoxia as Δx grew; the mask error deflated supply as Δx shrank. They pulled
in opposite directions on the same quantity.

### The number that makes the pattern legible

With the sink fixed but the mask still grid-dependent, mean O₂ at the fixed real cell positions ran
**0.522 / 0.426 / 0.251** across Δx = 30/20/12 µm — i.e. **after a correct fix, the grid dependence
of the output got worse, not better.** A team using "did the numbers get more stable?" as its
acceptance signal would have reverted the correct fix.

After the geometric mask (union of 30 µm discs around measured centroids), coverage runs
**67.2 / 67.4 / 67.7%** — spread **1.007×** — and mean O₂ runs 0.570 / 0.565 / 0.558, spread 1.02×.

### The downstream damage

Defect 7 was not confined to the oxygen field. It corrupted the production kill fraction by **2.25×**,
and — discovered a full experiment later — it was the sole cause of a headline finding in
`delivery_metrics.md` (pattern 4 below). One grid-dependent boolean mask propagated into three
separate published conclusions.

### How to transfer this

1. **A stable output is not evidence of a correct model.** Stability is evidence of *either* a
   correct model *or* a balanced set of errors. The two are distinguished only by checking each
   component's invariance separately, not the pipeline's end-to-end output.
2. **Check invariance component-by-component, not end-to-end.** The end-to-end grid study
   (`numerical_rigor.md` §2) was *itself contaminated* by the shrinking tissue mask — the study
   designed to detect grid dependence was measuring a domain that changed with the grid. The
   defect was found by asserting a property of the mask alone (`test_geometric_tissue_mask_
   coverage_is_dx_independent`), not of the simulation.
3. **When a correct fix makes the output look worse, that is information, not a regression.** Log it
   and look for the counterpart. In this project the sequence was: fix the sink → grid dependence
   worsens → conclude a second grid-dependent term exists → find the mask.
4. **The invariance checks worth writing are the ones about things that are not physics.** Grid
   spacing, time step, tolerance, random seed, unit system, array ordering, chunk size, worker
   count. A result that moves with any of these has a defect, regardless of whether the movement is
   in a plausible direction.

**The concrete check.** For each discretisation or bookkeeping parameter *p* that carries no
scientific meaning, and for each intermediate field *f* the pipeline computes, assert
`f(p₁) ≈ f(p₂)` under a declared tolerance — one assertion per field, not one for the pipeline. This
is what makes compensating errors decompose.

---

## 2. Trivial baselines are mandatory

> **The pattern.** Every model component must beat the cheapest thing that could plausibly work on
> the same inputs. If it does not, its complexity is not earning anything, and no amount of internal
> validation ("it converges", "it conserves mass", "it is significant") changes that.

### The case

The model solves a steady reaction–diffusion PDE for oxygen — a sparse SPD linear system, solved
directly, mass-conservative, grid-converged, deterministic, all verified. It was scored against
measured hypoxia-response transcripts (VEGFA, HIF1A) in 816,769 cells of a real GBM section, against
three baselines built from **the solver's own inputs**:

| predictor | what it is | Spearman ρ vs measured hypoxia score |
|---|---|---:|
| **simulated O₂ (the PDE)** | the reaction–diffusion solve under test | **−0.1132** |
| `vascular_kde` | the solver's own source term, *before* it is diffused | **−0.1270** |
| `dist_endothelial` | one `cKDTree` query to the nearest endothelial cell; no physics | **−0.1503** |
| `local_cell_density` | pure local crowding; ignores vessels entirely | −0.0192 |

All three non-trivial predictors are significant at the permutation floor (p = 0.00050). The PDE
loses to the kd-tree query by **33%**, and it loses to its own undiffused input. Partialling out
distance-to-vessel leaves the PDE with ρ = −0.0660; partialling out the PDE leaves distance with
ρ = −0.1192. **The PDE is the redundant variable in the pair.**

The failure is not a resolution problem. Refining the mesh over a 2.5× range of Δx moves the
correlation by **less than 0.0002** — the field has converged onto an answer that a distance
transform beats.

`local_cell_density` at −0.0192 is the control that makes the comparison meaningful: it establishes
that not every spatial variable correlates with the marker, so the other three numbers are signal
rather than spatial-autocorrelation artefact.

### The second instance, in a different subsystem

`moncayo_explant.md` §3.4: the simulator's slice-level prediction of drug accumulation reached
ΔR² = **6.7%** beyond the random effects. The one-number geometric covariate it is standing in for —
the section's mean distance to its own edge — reached **4.3%** in the same statistical frame with no
simulator at all. The simulator's prediction correlates with that covariate at Spearman −0.948
(p = 9.4e-06). The transport model was a section-thickness proxy with extra steps.

### How to transfer this

1. **The baseline must use the same inputs as the model.** A baseline the model can beat because it
   is fed worse data proves nothing. `dist_endothelial` is computed from *the PDE's own source
   cells*. That is what makes losing to it decisive.
2. **Include the model's own input as a baseline.** "Does the transformation add information, or
   destroy it?" is a distinct and often more damaging question than "does the model beat nothing?".
   Here the diffusion step took a field correlating at −0.1270 and returned one correlating at
   −0.1132.
3. **Statistical significance is not a baseline.** The PDE's p-value was at the estimator's floor.
   Significance measures distance from zero; a baseline measures distance from *the cheap
   alternative*, which is the only comparison that licenses complexity.
4. **Include a null-spatial or null-mechanism control** (here `local_cell_density`) so a passing
   baseline comparison cannot be attributed to a generic property of the data.
5. **Losing to a baseline does not falsify the mechanism** — it falsifies the *value-add of computing
   it this way*. The honest verdict recorded here was "PARTIALLY SUPPORTED: not falsified, out-
   competed", and the consequence drawn was to derive hypoxia from measured markers where they exist
   rather than to delete the solver.

**The concrete check.** For every model component, name the cheapest predictor of the same target
that can be computed from the same inputs in one line, and score both. Write the baseline *before*
running the model, so the bar cannot be set at whatever the model achieves.

---

## 3. Circularity: when your label and your prediction share an input

> **The pattern.** A prediction that appears to validate the model is a restatement of one of its
> inputs, because the *label* being predicted and the *quantity* predicting it are both functions of
> the same measured signal. This produces large effect sizes and small p-values, and it is not
> detectable from the statistics.

### The case

The claim under test: *perivascular regions of a real GBM section receive more drug.* The result was
strong and consistent across all five drugs — exposure ratio **1.84×–2.62×** perivascular vs
core+necrotic, Cohen's *d* = +1.21 to +1.40, Cliff's δ = +0.61 to +0.69, spatially-structured
permutation p ≤ 0.004 (using a rotation null, not an i.i.d. shuffle, because the i.i.d. null is
anti-conservative for spatial data).

It is a pipeline consistency check, not a prediction:

- the `perivascular` **label** is assigned by rank correlation against Ivy GAP signatures that are
  driven by endothelial transcripts — `PECAM1`, `VWF`, `CD34`, `CLDN5`, `FLT1`;
- the drug **source term** is `vascular_density × bbb_penetrance × plasma`, where `vascular_density`
  is a kernel density estimate of the *deconvolved endothelial spot fraction* — the same endothelial
  marker expression, at the same spots.

"Drug arrives where the endothelial signal is" is the model's input, restated.

### The diagnostic that makes it decidable

Strip out every contrast that involves the shared signal and re-score. The zones whose labels do
*not* feed the vascular field are `core` (from `EGFR`/`SOX2`/`OLIG2`/`MKI67`) and `necrotic` (from
HIF targets `CA9`/`VEGFA`/`ADM`/`NDRG1`):

| contrast | ratio | Cliff's δ | p (rotation null) |
|---|---:|---:|---:|
| perivascular vs core | 2.15× | +0.701 | 0.008 |
| perivascular vs necrotic | 1.72× | +0.565 | 0.013 |
| **margin vs necrotic** | **0.89×** | **+0.003** | **0.32** |
| **core vs necrotic** | **0.80×** | **−0.243** | **0.70** |

**The non-circular half fails.** The textbook expectation is that the necrotic zone is the
worst-delivered compartment; on this section it receives 1.25× *more* predicted drug than the core,
and the difference is not significant. Once the shared input is removed, only the contrast that
contains it survives.

### The independent confirmation, from a different experiment

`vessel_network_impact.md` §8.1 replaced the smooth endothelial KDE with a discrete segmented vessel
field — a step function over vessel objects rather than a blur of the marker. That construction
**eliminates the circularity** (the label and the source no longer share the smoothed signal), which
is one of only four reasons that report gives for adopting the discrete field at all. The
circularity was real and it was fixable by changing the *representation*, not the statistics.

### How to transfer this

1. **Trace every label and every predictor back to the raw measured channel.** Circularity is a
   property of the provenance graph, not of the numbers. If two nodes share an ancestor that is doing
   the work, no statistic will tell you.
2. **The tell is an effect size that is too clean.** Cliff's δ of +0.69 across five chemically
   distinct drugs, all landing within 0.08 of each other, is a geometry result wearing a
   pharmacology costume.
3. **Build the ablation into the experiment, not the discussion.** The useful artefact here is the
   table above — the pre-registered set of contrasts that do *not* share an input — reported beside
   the headline, not as a caveat after it.
4. **Report the circular result as what it is.** "The pipeline propagates real vascular structure
   end-to-end" is a genuine and useful statement. It is not "the model predicts delivery." Both
   sentences describe the same numbers.

**The concrete check.** For every claim of the form "X predicts Y", write out the derivation chains
of X and Y to the raw data. If they intersect, either re-derive one of them from a disjoint channel,
or reclassify the result as a consistency check.

---

## 4. Source nulls: what does your measurement return when the mechanism is off?

> **The pattern.** A fitted parameter that looks like a physical quantity may be measuring the
> smoothing, binning, or resolution of the pipeline rather than the physics. The only way to know is
> to run the identical fit on an input where the mechanism is switched off, and compare.

### The case

`penetration_length_L` is fitted from `AUC(d) = C₀·exp(−d/L) + floor` against distance from the
nearest vessel. It is dimensionally a physical length in µm, and drugs with 2× differences in
effective diffusivity and 10× differences in loss rate should differ in it.

They did not:

| drug | fitted L (µm) |
|---|---:|
| temozolomide | 28.2 |
| lomustine | 27.2 |
| cisplatin | 27.7 |
| methotrexate | 26.7 |
| doxorubicin | 26.5 |
| **SOURCE NULL** — identical fit on the raw `vascular_density` field, no transport solve at all | **23.7** |

Cross-drug spread **1.07×**; drug/null ratio 1.12–1.19×. The source field was a Gaussian KDE of
endothelial positions at ~40 µm bandwidth. **The "penetration length" was the smoothing bandwidth.**

The source null costs nothing to compute — it is the same fitting function applied to a field the
pipeline already has in memory — and it converts an apparently physical measurement into an obviously
degenerate one.

### The source null also solved the follow-up mystery, without running the model

Months of experiment later, the same fit on the same section returned L = 178 / 109 / 58 µm for three
drugs — a 3.07× spread, against a source null of 21.3 µm. The published 26–28 µm result no longer
reproduced. The cause was traced **without a single transport solve**, because the source null is a
fit to the source field alone: refitting it under the *retired* occupancy tissue mask returned
**23.68 µm**, matching the published **23.7 µm** to three significant figures.

The degenerate penetration length was a symptom of defect #7 (the grid-dependent tissue mask), not of
the KDE bandwidth. A cheap null was both the detector of the original problem and the instrument that
identified its true cause.

### The strongest form of the null: a source that admits no fit

`vessel_network_impact.md` §5 makes the point sharply. Under the discrete vessel field the source is
a step function — exactly zero beyond one bin — so the three-parameter exponential fit **cannot be
applied to it at all**. The report records the source null as `undefined`, and treats that as the
finding:

> Under the KDE the null is a real, competing explanation; under the discrete source there is nothing
> left to explain the drug's L except transport.

A representation that makes the null undefined is strictly better than one that makes the null small.

### How to transfer this

1. **Every fitted quantity needs a mechanism-off control run through the identical estimator.** Not a
   theoretical argument that the mechanism should dominate — the same code path, on the same data,
   with the mechanism disabled.
2. **Report the null beside the value, always.** "L = 27.7 µm" is a claim. "L = 27.7 µm against a
   source null of 23.7 µm" is a refutation of that claim, in the same line.
3. **Candidate nulls, by pipeline stage:** the raw input field before any solve; the output with the
   coefficient of interest set to zero; a shuffled or rotated version of the spatial structure; a
   constant field of matched total mass; the point-spread function or bin width of the instrument.
4. **When the null is close to the signal, the metric is measuring the pipeline.** Ratio 1.12× is not
   "close enough to be encouraging". It is the answer.

---

## 5. Scale-invariant metrics that cannot discriminate

> **The pattern.** A metric that is a *ratio* or a *normalised shape* is invariant to the thing you
> care about by construction. It will still vary a little across conditions, and that variation can
> be smaller than the metric's own numerical noise. Before trusting any metric, check that its signal
> exceeds its own noise floor.

### The case

`penetration_index = mean(AUC) / p99(AUC)` over the tissue mask. Multiply the entire exposure field by
any constant and it does not move — a property asserted arithmetically in the module's self-test:
scaling the field by 10⁻³ and 10³ leaves it identical to 1e-12.

Across a five-drug panel whose mean delivered exposure differs by **120×**:

| | `penetration_index` (retired) | `therapeutic_coverage` @1× IC50 (replacement) |
|---|---|---|
| range across the panel | 0.1946 – 0.1959 | 0.3108 – 0.9248 |
| **across-drug spread** | **0.0014** | **0.614** |
| grid drift 15→12 µm, worst drug, absolute | 0.0183 | 0.0281 |
| **separation = spread / drift** | **0.082** | **21.9** |

**The grid drift is 13× the entire across-drug spread.** Refining the mesh moves the metric an order
of magnitude more than changing the drug does. It was not measuring delivery; it was measuring the
vascular geometry of the section, which is identical for all five drugs.

The consequence was not cosmetic. Ranked by `penetration_index` the panel comes out
`cisplatin > lomustine > temozolomide > methotrexate > doxorubicin`. Ranked by the replacement metric
it comes out `cisplatin > doxorubicin > lomustine > methotrexate > temozolomide` — which is the
*identical* order the simulated kill fraction produces from an independent computation. The retired
metric's ranking was noise.

### The rule, pre-declared

`delivery_metrics.md` §3b states the bar before the numbers: `separation = spread / worst-drug grid
drift`, both on the same scale, must be **≥ 3**. Under that rule, five of eighteen candidate metrics
fail — and three of the five failures are *informative*, not bugs:

- `penetration_length_L` fails because it measures the source (pattern 4).
- `lowO₂/highO₂ exposure` fails at 0.335–0.337 across the panel because that is the *correct* answer:
  the perfusion gradient is a property of the tissue and all five drugs ride the same vascular field.
- `penetration_index` fails because it is scale-invariant by construction.

A metric that fails a discrimination bar is not necessarily wrong. It may be faithfully measuring
something that does not vary across your conditions. The error is quoting it as though it does.

### The scale rule matters too

Bounded fractions and unbounded quantities cannot be judged on the same convergence scale. Relative
error is ill-conditioned for a fraction — it diverges as the fraction approaches zero even when the
answer moved by a fraction of a point. `delivery_metrics.md` declares, before the numbers: bounded
[0,1] metrics are judged on **absolute** change with tolerance 0.03; everything else on **relative**
change with tolerance 5%. Both columns are printed for every metric so any other rule can be applied
to the same table.

The honest reading was then recorded even where it was awkward: the chosen headline metric moves
0.0281 absolute (inside tolerance) but 8.3% relative (outside the other tolerance), driven entirely
by the least-covered drug — so it is quoted as "converged to ±0.03 in coverage units, **not** to 5%
relative", with that error bar attached everywhere it appears.

### How to transfer this

1. **Ask what the metric is invariant to, and whether that is the thing you are trying to measure.**
   A dimensionless ratio is invariant to overall magnitude; a normalised shape is invariant to
   amplitude; a percentile ratio is invariant to units. If your effect lives in the invariant
   direction, the metric cannot see it.
2. **Compute signal-to-own-noise before interpreting any metric.** Across-condition spread divided by
   the metric's drift under a parameter that carries no meaning. Declare the threshold first.
3. **Prefer metrics computed over a fixed population.** `therapeutic_coverage` is stable partly
   because it is evaluated at the positions of a fixed set of 30,174 measured cells, so mesh
   refinement cannot change the denominator. Metrics averaged over a rasterised mask inherit every
   change in that mask.
4. **Prefer metrics with units that make the comparison meaningful.** Coverage compares a
   concentration to a concentration (exposure vs IC50); the retired index compared a field to itself.
5. **Retire, do not patch.** `penetration_index` was not wrong — it is a faithful measurement of
   exposure-field uniformity. The fix was to stop quoting it as a delivery number and, if a
   uniformity statistic is wanted, to rename it `exposure_uniformity_index` and say what it is.

---

## 6. Pre-declared bars

> **The pattern.** A threshold chosen after seeing the result is not a threshold. Goalpost-moving is
> invisible after the fact — including to the person doing it — so the only defence is to write the
> bar down before the run.

### The case

`ranking_robustness.md` asks whether a drug ranking is decision-grade. §1 defines the answer's
acceptance criterion **before §4 reports any probability**:

> Decision-grade is defined here, in advance, as: the recommended drug holds rank 1 with probability
> whose 90% credible interval lies entirely above **0.90**, and the runner-up's lies entirely below
> **0.10**.

The result, from 36 tissue draws × 256 IC50 draws = 9,216 complete re-rankings:

- P(cisplatin = rank 1) = **0.927**, cluster-bootstrap 90% CI **[0.924, 0.931]** — clears 0.90.
- P(doxorubicin = rank 1) = **0.044**, 90% CI **[0.042, 0.047]** — clears 0.10.
- **Verdict: the top-1 call is decision-grade.**
- P(top-2 set = {cisplatin, doxorubicin}) = **0.565**, 90% CI [0.560, 0.571] — **fails the same bar**,
  and is reported as failing.

The bar passed one claim and failed another in the same report. That is the evidence that it was a
bar and not a rationalisation.

The pre-declaration also constrained an earlier judgement in the opposite direction. `SESSION_AUDIT.md`
§6 item 7 recorded, before the Monte Carlo existed, that cisplatin (0.241) and doxorubicin (0.231)
differed by ~4% and that "the top-1 call is **probably not** decision-grade and should be reported as
a tie until shown otherwise." The measurement then contradicted that prior. Having written the prior
down is what makes the contradiction informative rather than invisible.

### Other pre-declarations in the same body of work

- `delivery_metrics.md` §2: convergence tolerances (0.03 absolute for bounded fractions, 5% relative
  otherwise) stated before the convergence table, with both scales printed for every metric.
- `delivery_metrics.md` §3b: the discrimination bar (`separation ≥ 3`) stated before the
  discrimination table.
- `moncayo_explant.md` §4: "*Passing (a) at 1/175th of the observed magnitude is not a pass*" —
  written into the description of the test, before the 175× shortfall was measured.
- `moncayo_explant.md` §1: the variance-decomposition table showing that 89.9% of the measured
  variance lives *inside* a slice, printed **before** the fit, establishing the structural ceiling on
  what a slice-resolved prediction could possibly explain. "We state this before reporting the fit
  rather than after."

### How to transfer this

1. **Write the acceptance criterion into the analysis module, not the write-up.** A constant in code
   that the report reads is harder to revise than a sentence in a discussion section.
2. **State the ceiling before the result.** If the design structurally caps what can be explained,
   compute the cap first. It reframes a "6.7% variance explained" from a modest success into a number
   that has to be compared against 5.0%.
3. **A bar is only credible if you can point at something it failed.** If every pre-declared bar in a
   document is passed, the bars were probably set afterwards.
4. **Record the prior expectation, including when it turns out wrong.** The value is in the audit
   trail, not in being right.

---

## 7. Refusing your own positives

> **The pattern.** Adversarial self-review means applying the same scepticism to results that support
> the model as to results that attack it. In practice this means declining to report significant
> findings, declining better-fitting models, and declining acceptance criteria — including ones
> imposed from above — when they are protecting an error.

Three instances from this project, each a different form of the discipline.

### 7a. Declining a significant result on power grounds

`state_validation.md` §2.3 tested whether cells the model labels VIABLE are more proliferative. The
point estimate went the right way: MKI67-positive in **2.22%** of VIABLE cells against 0.444%
section-wide — a **5.0× rate ratio**, permutation **p = 0.0078**, surviving depth stratification
(p = 0.0069).

It was reported as **not to be believed**, for reasons computed in the same table:

- it rests on **4 positive cells out of 180**; Wilson interval 0.87–5.57%;
- the minimum rate ratio the contrast could detect at α=0.05, power=0.8 was **4.14×**, and the
  measured biology it is supposed to detect is a *depletion* of at most 2× (Brat 2004, PMID 14871821)
  — a test that needs 4.1× to fire cannot see an effect of that size in either direction;
- **TOP2A**, the other proliferation marker on the same cells, shows nothing: 1.11%, 2/180, p = 1.0.

The verdict recorded was **UNDERPOWERED** — explicitly *not* "the state model is right" and *not*
"the state model is wrong". And the better-powered contrast on the same section (HYPOXIC vs
QUIESCENT, n = 4,051 vs 25,943, minimum detectable ratio 1.71×) points the *wrong* way: both
proliferation markers are higher in the model's HYPOXIC band, the opposite of the measured tissue
biology. That contradiction is reported in the same section as the favourable point estimate.

> Reporting it as support would be the exact error this module exists to prevent.

### 7b. Declining a better-fitting model

`mgmt_scalar.md` §1 fits a per-patient sensitivity scalar from measured MGMT expression. Two forms
were fitted; leave-one-out cross-validation is unambiguous:

| model | LOO MAE (percentage points of viability) |
|---|---:|
| constant (what the model did before) | 16.98 |
| **linear in viability** | **13.24** |
| exponential in kill (chosen) | 14.52 |

The linear model wins out-of-sample by 1.28 pp. It was rejected, because it predicts **zero kill at
x = 5.08 and negative kill above that** — and **three of the four real Visium sections in this project
sit at x ≥ 4.0**. An unbounded multiplier is not usable in the engine, and the region where it
breaks is the region the data occupies.

The trade is stated in the report rather than hidden: 1.28 pp of cross-validated accuracy, bought for
a bounded multiplier. Choosing the worse-fitting model was correct *because of where the data lives*,
not on general principle.

### 7c. Declining an acceptance criterion set from above

`glioma_regional_pk.md` P2 records the sharpest case. The brief for that work — set by the
orchestrator, not by the person doing the analysis — required that aggregate fold-bias not worsen
against the pre-fix values of **2.02× overall / 1.01× excluding one outlier**.

After the two fixes, the model gives 1.07× overall (better) but **0.73× ex-outlier** — farther from
1.0 than the pre-fix 1.01×. The report states plainly: **"That requirement is not met on the signed
ex-eribulin measure."**

It then refuses the requirement:

> The pre-fix 1.01× was two large errors cancelling — temozolomide +2.86× and cefazolin +1.66×
> annihilating adavosertib 0.26× — so it is the wrong thing to protect.

The unsigned error over the same five drugs was 2.06×. On every non-cancelling measure the fixed
model is better: unsigned geometric-mean fold error 3.66× → **1.99×** (2.06× → 1.55× ex-outlier),
worst case 65.6× → **6.9×**, within-10× band 5/6 → **6/6**.

This is pattern 1 recurring at the level of the *acceptance metric*: a signed, symmetric error
measure lets over- and under-prediction cancel, so it rewards compensating errors. An acceptance
criterion built on it will actively defend them.

### How to transfer this

1. **Compute the minimum detectable effect for every negative or weakly-positive test, and compare it
   to the effect size the literature says exists.** A test that cannot resolve the real effect
   produces "no evidence", not "evidence of no effect" — and its occasional significant result is
   noise.
2. **Require corroboration from a second readout on the same units.** MKI67 fired; TOP2A did not.
   Two markers of the same biology on the same 180 cells disagreeing is decisive against the one that
   fired.
3. **Prefer bounded, physically admissible functional forms over better-fitting unbounded ones when
   the data sits near the boundary.** State the cross-validation cost of doing so.
4. **Never accept a signed/symmetric error metric as an acceptance criterion.** Use unsigned error,
   worst case, and coverage bands. If a criterion is imposed on you, and it is a cancelling measure,
   say so and report both.
5. **The verdict vocabulary matters.** This project distinguished FALSIFIED / PARTIALLY SUPPORTED /
   UNDERPOWERED / NOT TESTABLE. Collapsing "underpowered" into either "supported" or "refuted" is the
   most common way an honest analysis becomes a dishonest one.

---

## 8. Structural vs tunable failure

> **The pattern.** When a model misses a measurement, the first question is not *how far off* but
> *whether any parameter choice could have reached it*. If the measured value lies outside the
> model's mathematically reachable set, calibration is wasted effort and the answer is a new
> mechanism. This distinction is usually cheap to establish analytically and it determines whether
> you tune or rebuild.

### Case A — a reachability floor (GIST imatinib)

The model predicts a tumour:normal drug partition ratio from a 1-D reaction–diffusion slab:

```
E_t / E_n  =  [(L/R_t)·tanh(R_t/L)] / [(L/R_n)·tanh(R_n/L)]
```

As `L → 0` (infinitely poor penetration — the limit most favourable to a transport explanation of low
tumour uptake), `E = (L/R)·tanh(R/L) → L/R`, so the ratio tends to `R_n/R_t`. As `L → ∞` it tends
to 1. **The transport-only prediction is therefore confined to `[R_n/R_t, 1]` whatever L is.**

Sweeping L over five orders of magnitude confirms it numerically:

| L (µm) | 1 | 5 | 10 | 20 | 36 | 100 | 400 |
|---|---|---|---|---|---|---|---|
| predicted liver-pair ratio | 0.208 | 0.211 | 0.246 | 0.374 | 0.580 | 0.900 | 0.993 |

- Floor at the nominal geometry (normal-liver half-spacing 12.5 µm, tumour 60 µm): **0.208**.
- Floor at the most model-favourable end of the literature band (10 µm / 120 µm): **0.0833**.
- Measured ratios in 5 matched human liver pairs: **0.023, 0.043, 0.044, 0.047, 0.142**.

**5 of 5 below the nominal floor; 4 of 5 below even the extreme floor.** No penetration length and no
intervessel geometry in the literature band can reproduce them. Verdict: falsified as a sufficient
explanation, and the failure names the missing mechanism — active hepatic uptake (OCT1/OATP) plus
lysosomal ion-trapping of a weak base, concentrating drug in the *normal* side. Neither is a diffusion
process; neither exists in the model.

Note what this bought: a definitive negative on **n = 6** matched pairs, of which 5 are the same
tissue pairing — a dataset the report itself says carries "approximately ONE independent
discriminating comparison" and "cannot rank two plausible models". A reachability argument extracts a
falsification from a dataset too small to fit anything.

### Case B — a symmetry argument (Moncayo explants)

The claim under test: lysosomal sequestration of weak bases drives spatial drug heterogeneity in
bath-dosed tumour explants. Implemented exactly, as a retardation factor `R = 1 + φ_lys·(ratio − 1)`
from Henderson–Hasselbalch ion trapping.

Predicted coefficient of variation of accumulation, against measurement:

| | predicted (weak base) | measured, IS-noise-corrected |
|---|---:|---:|
| CV | 2.75e-03 | 0.481 |

Off by **175×**. Two parameter sweeps say it is not a tuning problem:

| φ_lys | 0.005 | 0.014 | 0.05 |
|---|---|---|---|
| predicted CV | 6.79e-05 | 2.74e-03 | 5.62e-02 |

At φ_lys = 0.05 — at or above every value in the cell-biology literature — it is still **8.6× short**.
Sweeping the ECM tortuosity over (0.2, 0.4, 0.6) moves it across 1.53e-03 / 2.75e-03 / 5.09e-03.

But the sweep is not the argument. The argument is one line:

> With `φ_lys` **uniform in space**, R is uniform, so at equilibrium `C_total = R × C_bath`
> **everywhere** and the dispersion goes to zero however large R becomes. Uniform sequestration buys
> a longer time constant, never an equilibrium pattern.

A uniform retardation factor **cannot** produce spatial dispersion at equilibrium at any R. And the
system *is* at equilibrium: the Fourier number at the 24 h endpoint has median 7.4 with trapping
(minimum 1.13), 18/18 sections above 1, with mean tissue saturation 0.998.

The same report then makes the structural point about the *other* half of the mechanism: the model
recovers the weak-base / non-base chemistry correctly from SMILES alone (predicted heterogeneity
ratio **46,209×** between weak bases and olaparib, matching the literature pKa cross-check exactly)
but cannot distinguish niraparib from rucaparib — predicted CV ratio 0.994 — because the
Henderson–Hasselbalch trapping ratio saturates at 10^(7.2−4.7) ≈ 316 for any sufficiently strong base.
**Any pH-partition model has this ceiling.** That is a structural resolution limit, not a
parameterisation.

The missing ingredient is named precisely: a *spatial field* of acidic-compartment content. And the
report establishes that this is unexplored rather than unmodelled — no measured mg/g or mol%
acidic-phospholipid content exists for tumour tissue and no published PBPK model carries
tumour-specific lysosomal parameters.

### Case C — a ceiling artefact, correctly *not* treated as a physics error

`glioma_regional_pk.md` §6: the model predicts a tissue:plasma partition ≤ 1 by construction.
Adavosertib's measured contrast-enhancing ratio is **4.18**. The model is 3.8× low.

That is not a tuning failure and not a physics failure — it is a **quantity mismatch**. 4.18 comes
from tissue homogenate (total drug, including intracellular and lipid-bound); the model predicts a
free-drug crossing fraction. A free-drug partition of 4.18 would require active uptake into brain.
The model's 0.639 is in fact close to adavosertib's *non-enhancing* measured ratio of 0.74 (0.86×).

Every row of that scorecard is therefore labelled with its assay modality, and the modality-matched
subset is scored separately. Distinguishing "wrong number" from "different quantity" is the same
skill as distinguishing "wrong parameter" from "missing mechanism".

### How to transfer this

1. **Derive the reachable set analytically before fitting.** Take limits of your model's expression
   with respect to the free parameters. If the measurement lies outside `[min, max]`, stop tuning.
2. **Look for symmetry arguments that kill a mechanism outright.** "A spatially uniform multiplier
   cannot produce spatial structure at equilibrium" is a stronger statement than any parameter sweep,
   costs nothing to derive, and cannot be argued with.
3. **Sweep anyway, and report both.** The sweep is what makes the analytic argument checkable by a
   reader who does not trust the algebra.
4. **Check whether the measurement and the prediction are the same physical quantity** before
   concluding either is wrong.
5. **A structural failure is more valuable than a tunable one.** It names the missing mechanism, it
   survives on tiny n, and it cannot be undone by a future recalibration.

---

## 9. Dose- and effort-matched controls

> **The pattern.** When comparing two implementations of a component, they must be matched on
> everything except the property under test. If arm B differs from arm A in both *structure* and
> *total input*, the measured difference is a sum of the two, and attributing it to structure can
> overstate the effect by an order of magnitude or more.

### The case

The single highest-value structural fix identified by the audit was to replace the smooth endothelial
KDE drug source with a discrete segmented vessel network. Two arms were the obvious design:

- **A** — 40 µm Gaussian KDE of endothelial positions, divided by its own maximum. Nonzero over
  99.71% of tissue.
- **B** — amplitude 1.0 within 7.5 µm of a segmented vessel, exactly 0.0 elsewhere. Nonzero over
  2.50% of tissue.

**A delivers 3.36× more total drug into the tissue than B**, because a max-normalised KDE spread over
100% of the tissue integrates to more than a unit source on 2.5% of it. Any A→B difference is part
geometry and part dose.

A third arm was added: **C** — arm B rescaled so the tissue-integrated total matches A.

| doxorubicin, therapeutic coverage @1× IC50 | value |
|---|---:|
| A (KDE) | 0.9692 |
| B (discrete, naive comparison) | 0.3836 → **Δ = −0.5855** |
| C (discrete, dose-matched) | 0.9526 → **Δ = −0.0166** |

**The naive comparison overstates the structural effect by 35×.** The dose-matched change of −0.0166
is inside the ±0.03 grid-convergence band already established for that metric — i.e. the structural
change is, on the headline metric, not detectable.

The report puts this first, above every scientific finding:

> Every structural claim below is scored on arm C, not arm B. […] Reporting A→B as the structural
> effect would have overstated it 35× on the headline metric — this is the single most important
> methodological point in the report.

### The matched control also validated a *positive* result

The same design confirms which effects are real. Fitted decay length L is scale-invariant, and arms B
and C — which differ *only* by the 3.36× dose factor — agree on L to within **2.7%**. So the
widening of cross-drug separation from 3.12× (A) to 6.78× (C) is pure geometry, not dose. Matched
controls are not only for deflating claims.

### The cost of the control, stated honestly

Arm C is not free of side effects, and the report says so: matching total drug also scales the oxygen
and interstitial-pressure source, so arm C is a slightly hyperoxic tissue. Coverage and the exposure
profile do not depend on the oxygen field, so they are unaffected; kill fraction does, which is why
kill fraction is **not** a headline in that report. Both arms are reported; neither is hidden.

### How to transfer this

1. **Enumerate what differs between arms besides the thing under test.** Total input mass, compute
   budget, number of parameters, training data volume, random draws, resolution. Any of these can
   carry the effect.
2. **Add the matched arm even when it is a third of the compute.** This study ran 3 source fields ×
   3 drugs = 9 transport solves at 2815 s total; arm C was one third of that and it changed the
   conclusion.
3. **State which arm the claims are scored on, at the top.** Not in the methods.
4. **Report the side effects of matching.** Matching one quantity usually unmatches another. Name it
   and route the affected metrics around it.
5. **This generalises past "dose".** Effort-matched, parameter-count-matched, FLOP-matched,
   sample-size-matched, and wall-clock-matched controls are the same pattern.

---

## 10. Provenance and bookkeeping as correctness

> **The pattern.** A large fraction of defects in scientific code are not physics errors. They are
> conventions, units, aggregation choices, and stale artefacts. They produce numbers that are wrong
> by exact, derivable factors, and they are invisible unless provenance is treated as part of
> correctness.

### 10a. A 165× swing from a bookkeeping choice

`state_validation.md` §2.2 tabulates the model's tumour-cell state fractions across two tissue-mask
conventions and three grid spacings:

| tissue mask | Δx (µm) | HYPOXIC | QUIESCENT | VIABLE |
|---|---|---:|---:|---:|
| `geometric` | 30 | 0.00% | 1.81% | **98.19%** |
| `geometric` | 20 | 0.00% | 3.81% | 96.19% |
| `geometric` | 12 | 0.00% | 5.08% | 94.92% |
| `occupancy` | 30 | 0.00% | 6.71% | 93.29% |
| `occupancy` | 20 | 0.00% | 25.97% | 74.03% |
| `occupancy` | 12 | 13.43% | 85.98% | **0.60%** |

**The VIABLE fraction runs 0.60% to 98.19% — a 165× swing — and neither the mask nor the mesh is a
measurement.** Both are bookkeeping. The report states the consequence: under the current shipping
default the HYPOXIC state is *empty at every Δx*, so the hypoxia-resistance multiplier has no cells
to act on; the 86% quiescent / 13.4% hypoxic split quoted in an earlier results file comes from the
legacy convention that the reconstruction layer has since abandoned.

Downstream, ~89% of one drug's kill ceiling derives from the quiescent band via a ×0.5 multiplier
applied to a label with no measured support and a 165× convention-dependent size. The report draws
the distinction that matters: the dose/schedule **comparisons** are robust because they share the
multiplier and it largely cancels; the **levels** are not defensible as measurements.

### 10b. An error of exactly 1/ln²2

The calibration constant α had two live values in the tree: **309.6** and **148.6**. They differ by
exactly **1/ln²2 = 2.081**.

The cause: published perivascular immunofluorescence studies report `d_½`, the distance at which
intensity falls to half. The model's `L` is the 1/e distance. For `C(x) = C₀e^{−x/L}`,
`d_½ = L·ln2`, so `L = d_½/ln2 ≈ 1.443·d_½`. One caller fitted the raw 45 µm half-distance *as if it
were* an e-fold length; another passed the converted 64.9 µm target. Since α enters through
`L = √(D_eff/k)` with `k ∝ α`, the convention error squares into a factor 1/ln²2 in α.

Two properties make this instructive:

- **The error is exact and derivable.** An unexplained factor of 2.081 between two computations of
  the same constant is a convention error, not noise. Discrepancies with clean closed forms should be
  chased to their algebra, not averaged.
- **The old value passed its own out-of-sample test.** With α = 309.6, mitoxantrone's predicted
  L = 33.4 µm was compared against a 25–75 µm band — a band of *raw half-distances*. Both sides of
  the comparison used the wrong convention, which is why it passed. With α = 148.6, mitoxantrone
  (48 µm) and chemically distinct topotecan (84 µm) both land inside the correctly converted 36–108 µm
  band, 2/2. **A consistently-applied wrong convention validates itself.** The fix is now pinned by
  `test_calibrate_alpha_target_is_an_efold_length`.

### 10c. Results files going stale against a changing codebase

`vessel_network_impact.md` §3 opens by re-measuring the status quo before measuring the fix, and the
status quo does not reproduce:

| drug | L (µm), this run | L as published in `delivery_metrics.md` |
|---|---:|---|
| temozolomide | 178.3 | 26.5–28.2 |
| cisplatin | 108.6 | 26.5–28.2 |
| doxorubicin | 58.0 | 26.5–28.2 |
| cross-drug spread | **3.07×** | 1.07× |

Two published documents disagreed about the same quantity on the same section because an unrelated
fix (the tissue mask, defect #7) had changed every distance-binned fit in the tree. The report says
so directly:

> `SESSION_AUDIT.md` §3 bullet 1 is stale. […] That is a claim in a published results file that is
> now wrong in the conservative direction, and finding it is arguably the most useful thing this
> study did.

`vessel_network_impact.md` §8.4 then lists, **in dependency order**, the seven modules that must be
re-run if the vascular field ever becomes the default — starting with the calibration, because α is
fitted against the current field and every absolute exposure inherits it.

### 10d. Aggregation form as a defect

Defect #9 is not a bug in the ordinary sense. `bbb_penetrance` combined four descriptor desirabilities
with an **arithmetic mean**. Eribulin's two disqualifying properties (MW 730, desirability ~0.01;
TPSA 146, ~0.01) were averaged against its two favourable ones (logP 3.44, ~0.86; 2 HBD, ~0.78),
rescuing it to 0.334 against a measured 0.0051 — a **65.6×** miss.

> A mean cannot express a veto.

Switching to a geometric mean — no new parameter, and identical to the arithmetic mean when all four
desirabilities are equal, so it bites only on heterogeneous profiles — cut the miss to **6.9×**, took
the within-10× band from 5/6 to **6/6**, and improved temozolomide from 2.86× to 1.22× *without being
targeted at it*. Four candidate aggregations were scored against the same measured panel before one
was chosen, and the table of all four is in the report.

The same root pathology — unbounded or unshaped use of a single descriptor — produced defect #10, in
the same file: `membrane_permeability` returned **4.23e+07 µm/s** for imatinib (42 metres per second)
because it exponentiated a computed logP with a clamp far too high and no physical ceiling. The fix
put an unstirred aqueous boundary layer in series (`1/P = 1/P_mem + 1/P_aq`), which is the real
physical ceiling on passive permeability. Maximum over the drug library after the fix: **27.0 µm/s**.

The residual is reported too: three drugs now sit within 15% of each other because all three are
boundary-layer limited, whereas physically they differ by orders of magnitude. The output is "a
bounded, correctly-ordered-at-the-bottom ranking, not a calibrated rate. That is a smaller error than
42 m/s, but it is not zero."

### How to transfer this

1. **Treat unit and convention conversion as a first-class, tested operation.** A test that pins the
   *convention* (`target is an e-fold length`) is worth more than a test that pins the value.
2. **An unexplained clean factor between two computations of one constant is a convention error.**
   Chase 2.081, 1.443, 2, π/2, ln2 to their algebra.
3. **Bookkeeping parameters must be swept and their swing published beside the result they feed.**
   "The VIABLE fraction is 0.60%" is unusable; "0.60%, on a range of 0.60–98.19% across two
   conventions and three meshes, neither of which is a measurement" is honest.
4. **Every generated artefact needs the code version, the configuration and the inputs that produced
   it.** Every report in this project opens with its module path, its run command, its resolution and
   the value of the shared calibration constant, and several open with a provenance box naming what
   changed since the previous revision and which numbers are superseded.
5. **Re-measure the baseline before measuring the fix.** `vessel_network_impact.md` did this and
   discovered the defect it was built to fix had already been cured by an unrelated change.
6. **Maintain a dependency-ordered re-run list** for any change to a shared field or constant.
7. **Aggregation form is a modelling decision with a defect surface.** Mean, geometric mean, min,
   product and weighted forms encode different logics (compensatory vs vetoing). Score the candidates
   against measurement; do not inherit the default.

---

## 11. Self-application checklist

Run this against your own model. It is ordered by cost: the first six items are cheap and catch the
majority of what was found here.

**Invariance (cheap, catches compensating errors)**

- [ ] For every parameter with **no scientific meaning** — mesh spacing, time step, solver tolerance,
      seed, chunk size, worker count, array layout, unit system — is every *intermediate field*
      invariant to it within a declared tolerance? One assertion per field, not one for the pipeline.
- [ ] Do any of your densities, rates or sinks depend on a **count** where they should depend on a
      **density**? Check the dimensions of every term you wrote yourself.
- [ ] Is any mask, domain, or population **defined by the discretisation** rather than by the
      measurement? If so, its coverage is a free parameter.
- [ ] When you fix a defect and the output gets *less* stable, do you look for the compensating
      counterpart before reverting?

**Baselines and nulls (cheap, catches unearned complexity)**

- [ ] For every model component, what is the **cheapest predictor of the same target from the same
      inputs** — a nearest-neighbour query, a linear fit, a single geometric covariate — and does the
      component beat it?
- [ ] Does the component beat **its own input**, un-transformed?
- [ ] Is there a **null-mechanism control** in the comparison (a predictor that should not work) to
      establish that the effect is not generic?
- [ ] For every **fitted parameter**, what does the identical estimator return with the mechanism
      switched off? Is that number reported beside the fit?
- [ ] Is there a representation in which the null becomes **undefined** rather than merely small?

**Metrics (cheap)**

- [ ] What is each metric **invariant to** — and is the effect you care about in that direction?
- [ ] Is each metric's **across-condition spread ≥ 3× its drift** under a meaningless parameter?
      (Declare the threshold first.)
- [ ] Are bounded fractions judged on **absolute** change and unbounded quantities on **relative**
      change, with both printed?
- [ ] Are per-unit statistics computed over a **fixed population** so refinement cannot change the
      denominator?
- [ ] Do the metric's **units** make the comparison meaningful (a concentration against a
      concentration), or is it a field compared to itself?

**Circularity (cheap, needs a provenance graph)**

- [ ] For each "X predicts Y", do the derivation chains of X and Y back to raw measurement
      **intersect**?
- [ ] Is there a set of contrasts that **excludes** the shared input, and is it reported beside the
      headline?
- [ ] Are effect sizes suspiciously **uniform across conditions** that should differ?

**Structure vs tuning (cheap, high value)**

- [ ] For each miss: take the **limits** of the model's expression over its free parameters. Is the
      measurement inside the reachable set?
- [ ] Is there a **symmetry argument** that rules the mechanism out entirely (e.g. "a spatially
      uniform factor cannot produce spatial structure at equilibrium")?
- [ ] Are the prediction and the measurement **the same physical quantity** (free vs total, snapshot
      vs integral, per-cell vs per-area)?
- [ ] Does the model's own predicted feature size sit **below the resolution of the measurement**
      you are validating against? (Here: predicted L = 35.9 µm against a 50 µm raster = 0.72 pixels.)

**Controls (moderate cost)**

- [ ] What differs between your arms **besides the thing under test** — total input, compute,
      parameter count, sample size?
- [ ] Is there a **matched arm**, and are the claims scored on it, stated at the top?
- [ ] What did matching **unmatch**, and which metrics are routed around that?
- [ ] Is the statistical **unit of analysis** the exchangeable unit? (Here: slices, not ROIs — 11 not
      112; spots, not cells; patients, not samples. Getting this wrong inflated one ΔR² to a
      nominally significant 33.1% that a slice-exchangeable permutation dissolved to p = 0.091.)
- [ ] For spatial data, is the null **spatially structured** (rotation/reflection) rather than i.i.d.?

**Pre-declaration (free)**

- [ ] Is the acceptance criterion written into the **code**, before the run?
- [ ] Is the **structural ceiling** on what could be explained computed and printed before the fit?
- [ ] Can you point at something a pre-declared bar **failed**?

**Self-review (free, hardest)**

- [ ] For every weak positive: what is the **minimum detectable effect**, and how does it compare to
      the effect size the literature says exists?
- [ ] Does a **second readout on the same units** corroborate it?
- [ ] Is your verdict vocabulary rich enough to say **UNDERPOWERED** and **NOT TESTABLE**, distinctly
      from supported and refuted?
- [ ] Is any acceptance metric you are being held to a **signed or cancelling** one? If so, report
      unsigned error and worst case as well, and say why.
- [ ] Is a better-fitting model **inadmissible** (unbounded, negative, non-physical) in the region
      your data actually occupies?

**Provenance (moderate)**

- [ ] Is every unit/convention conversion **tested as a convention**, not just as a value?
- [ ] Do two computations of the same constant differ by a **clean closed-form factor**?
- [ ] Do bookkeeping parameters have their **swing published** beside the results they feed?
- [ ] Does every generated artefact carry the code version, configuration and inputs?
- [ ] Do you **re-measure the baseline** before measuring a fix?
- [ ] Is there a **dependency-ordered re-run list** for shared fields and constants?
- [ ] Was each **aggregation form** (mean / geometric mean / min / product) scored against
      measurement rather than inherited?

---

## 12. What this methodology costs

The methodology is expensive, and several of its own products were wrong. Both facts belong in the
write-up, because a methodology reported only through its successes is subject to exactly the
selection bias it exists to prevent.

### 12a. Compute

Single-run wall-clock for the verification analyses cited above, excluding failed and re-run
attempts:

| analysis | cost |
|---|---|
| `vessel_network_impact` — 9 transport solves at dx = 12 µm | **2,815 s** |
| `ranking_robustness` — 36 tissue draws × 5 drugs, 9,216 re-rankings | **42.4 min** |
| `delivery_metrics` — 3 grids × 5 drugs | **661 s** |
| `state_validation` — 6 reconstructions of an 816,769-cell section + a 480-gene scan | **~8 min** |
| the data-landscape sweep behind `DATA_HUNT.md` | **39 agents, 164 datasets** |

Not counted: the runs that produced the retracted analyses in 12c, the grid ladders re-run after each
of defects 6 and 7, and the re-runs triggered by α changing by a factor of 2.08.

**The load-bearing observation about that cost**: the *decisive* checks were nearly all cheap. The
source null is one curve fit on a field already in memory. The trivial baseline is one `cKDTree`
query. The GIST falsification is a limit of a `tanh`. The Moncayo falsification is one line of
algebra about uniform multipliers. The reachability floor extracted a definitive negative from n = 6.
The expensive experiment — the 2,815 s discrete-vessel study — was the one whose motivating premise
turned out to be wrong (12c).

**The transferable scheduling rule: run the cheap structural checks first.** Analytic reachability
bounds, symmetry arguments, source nulls, one-line baselines and invariance assertions cost minutes
and can eliminate the need for the expensive experiment entirely.

### 12b. Human cost

The pass produced 33 results documents, several exceeding 25 KB, against a codebase of 86 Python
modules in 10 subpackages. Most of that text is caveats, nulls, and negative results. The eleven defects were found
because the model was attacked, and attacking a model is slower than building one — this project's
own summary observes that a version of it that had skipped the audit "would still be reporting nine
wrong numbers with beautiful figures attached."

### 12c. Analyses that were themselves wrong

**The orchestrator's "vasculature blur" diagnosis — largely refuted.** Three independent analyses
converged on one root cause: the vascular source field is a smooth ~40 µm-bandwidth Gaussian KDE, and
that was blamed for (1) all drugs fitting the same penetration length, (2) a flat far-field exposure
floor at ~12% of peak, and (3) the perivascular circularity. Replacing the KDE with a discrete
segmented vessel network — the single highest-value structural fix in the audit — refuted **two of
the three**:

1. **The flat floor is not a KDE artefact — it is diffusion physics.** With a source that is exactly
   zero over 97.5% of the tissue, far-field exposure moves only **1.12×** (25.99% → 23.18% of
   perivascular). The reason is arithmetic: reaction–diffusion length is 73–296 µm and the diffusion
   length over the 6 h window is ~4,000 µm, against a median intervessel gap of 40 µm — one to two
   orders of magnitude larger. No tumour cell in any arm receives zero drug; the essentially-zero
   fraction is 0.11–0.17% in every arm. The floor is a property of the tissue, not of the smoothing.
2. **"Penetration length cannot discriminate drugs" was stale**, and the cause was defect #7 (the
   tissue mask), not the vasculature. Demonstrated without running a solve, by refitting the source
   null under the retired mask: 23.68 µm, matching the published 23.7 µm to three significant figures.
3. **Only the circularity diagnosis stood**, and the discrete field does eliminate it.

The diagnosis that motivated the most expensive experiment in the pass was two-thirds wrong, and the
experiment's most useful output was discovering that.

**A model that fit better and was rejected as a knob.** In the GIST validation, model M2 (transport ×
tissue composition) lowered geometric-mean fold error from 7.92× to 5.95×. It was refused, for two
reasons stated together: it **inverts the ordering** (Spearman −0.655 vs +0.655), predicting the
best-perfused pair to be the worst; and swapping the arbitrary normal-tissue surrogate used to build
its tumour composition moves the answer by **4.5×** (0.080 to 0.357) — "a free choice among three
equally-unjustified surrogates that spans the answer is a knob, not a model."

**A rank correlation that was one data point.** The same report's M1 achieves Spearman ρ = 0.655
across 6 patients. It is reported as evidence of nothing: the model predicts an *identical* value for
all five liver pairs, so the entire correlation comes from a single intestinal pair sorting above a
tied block. Within the five liver pairs the measurement spreads 6.1× and the model spreads 1.0× —
0% of the within-tissue variance.

**A calibration that was attempted and abandoned.** Absolute quantitation of the GIST MSI data was
genuinely attempted and then rejected, with seven independently sufficient blockers documented. The
report names the tempting alternative explicitly: sorting six unlabelled standard spots to force
monotonicity and fitting a line "would have produced a beautiful R² that meant nothing."

**Two nulls disagreeing, reported as the finding.** In the Moncayo analysis, the permutation test on
|t| gives p = 0.044 and the one on ΔR² gives p = 0.074, while the unweighted slice-level Spearman is
0.009. Rather than pick one, the report states: "An effect that appears or vanishes depending on
whether slices are weighted by ROI count is not a validated prediction; it is 11 points arranged so
that one weighting scheme finds a trend."

### 12d. What the methodology did not do

- It did not find these defects **before** they were published internally. Every headline number
  reported before the pass was wrong in at least one way.
- It does not address **model-form error**, and this is plausibly the largest remaining term. Every
  Monte Carlo draw in the uncertainty analysis is computed inside one fixed structure; if that
  structure is wrong, no amount of parameter sampling detects it. The KDE→discrete swap alone moves
  one drug's coverage by a structural amount that **no prior in the sampling could produce**.
- It cannot manufacture **statistical power**. Several tests came back UNDERPOWERED or NOT TESTABLE
  and stayed that way; the honest output was a named measurement that would settle them.
- It cannot manufacture **data**. The pass established that no open dataset anywhere pairs measured
  spatial drug distribution with measured treatment response in the same tumour. The response axis of
  this model is not validatable with open data, and no amount of methodology changes that.
- Verifying a component does not validate it. Mass conserved to 4e-14, unconditionally stable,
  dt-converged and grid-converged — all true of the oxygen solve, and all compatible with being
  out-predicted by a kd-tree query.

### 12e. The net position

Eleven defects, two components retired, four external validations run against measured human data —
one partial pass on magnitude with an outright fail on regional structure, one pass on potency, and
two structural negatives that name their missing mechanisms. The output is not a
better simulator; it is a simulator with a quantified error budget and an explicit list of what it
cannot do — which the project's own summary calls "a fundamentally different object from a simulator
that produces confident numbers."

Whether that trade is worth its cost depends on what the model is for. For a model that will inform a
decision, the error budget *is* the product. For a model that will only ever generate hypotheses, most
of section 12a is unaffordable and the checklist items marked *cheap* in section 11 are the ones to
run.

---

## 13. Evidence index

Every number in this document is traceable to one of these files.

| pattern | primary evidence |
|---|---|
| 1 — compensating errors | `docs/MATHEMATICAL_MODEL.md` §4.2; `results/numerical_rigor.md` provenance box; `results/SESSION_AUDIT.md` §2 |
| 2 — trivial baselines | `results/state_validation.md` §1.2; `results/moncayo_explant.md` §3.4 |
| 3 — circularity | `results/niche_response.md` §4c, §6; `results/vessel_network_impact.md` §8.1 |
| 4 — source nulls | `results/delivery_metrics.md` §3a, §5; `results/vessel_network_impact.md` §3.1, §5 |
| 5 — metric noise floors | `results/delivery_metrics.md` §2, §3b, §4 |
| 6 — pre-declared bars | `results/ranking_robustness.md` §0, §1; `results/delivery_metrics.md` §2, §3b; `results/moncayo_explant.md` §1, §4 |
| 7 — refusing your own positives | `results/state_validation.md` §2.3; `results/mgmt_scalar.md` §1; `results/glioma_regional_pk.md` P2 |
| 8 — structural vs tunable | `results/gist_imatinib.md` §3.3; `results/moncayo_explant.md` §4.3; `results/glioma_regional_pk.md` §6 |
| 9 — matched controls | `results/vessel_network_impact.md` verdict, §2, §6 |
| 10 — provenance | `results/state_validation.md` §2.2; `docs/MATHEMATICAL_MODEL.md` §7.3, §8.3; `results/vessel_network_impact.md` §3, §8.4; `results/glioma_regional_pk.md` P1, P4 |
| 12 — costs and retractions | `results/SESSION_AUDIT.md` §3; `results/vessel_network_impact.md` §3, §4; `results/gist_imatinib.md` §3.3, §4, §5.3; `results/moncayo_explant.md` §3.3 |

Supporting: `results/reconstruction_fidelity.md` (what survived), `results/gbm_tmz_exvivo.md` (an
external validation that passed, and the forfeit it quantified), `results/DATA_HUNT.md` (the data
landscape), `docs/MATHEMATICAL_MODEL.md` (the derivations whose act of writing found defects 4, 6
and 7).
