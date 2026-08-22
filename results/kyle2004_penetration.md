# Does the map reproduce a MEASURED penetration difference? — mitoxantrone, daunorubicin, epirubicin and doxorubicin against Kyle 2004

**Status: RUN AND SCORED against a pre-registration frozen before any field was
solved.** Driver: `experiments/isomer_penetration_map.py`. Registration:
`results/sweep/kyle2004_penetration_prereg.json`. Data:
`results/kyle2004_penetration_data.json`. Figure:
`results/kyle2004_penetration_map.png`.

**Nothing under `tmesim/`, `app/`, `landing-3d/`, `docs/` or `data_cache/` was
modified. No default, constant or shipped parameter was changed.
`calibrate_alpha()` is asserted live at `22.223477711172563` in every stage and
holds. Nothing was fitted to the measurement.**

---

## 0. Verdict up front

**The map reproduces the measured ORDERING. It does not reproduce the measured
MAGNITUDE, and on the reconstructed tissue the compound difference does not clear
its own noise floor.**

**3 of 6 scored criteria pass.**

| | criterion | result |
|---|---|:--:|
| **C0** | epimer control — the engine must not separate a measured-identical pair | **PASS** |
| **C1** | **direction** — daunorubicin placed as the worse penetrator, both arms | **PASS** |
| **C5** | the map arm measures transport, not vessel geometry | **PASS** |
| **C2** | **noise floor** — separation bigger than nuisance spread | **FAIL** (clean-slab arm passes; **real-tissue map arm fails**) |
| **C3** | **magnitude** — modelled fold vs measured fold | **FAIL** — 1.22× against a measured **4.55×** |
| **C4** | absolute near:far ratios | **FAIL** — daunorubicin 3.65 against a measured **50 ± 15** |
| C6 | far-side disclosure | not scoreable — reported, and it **disagrees** |
| C7 | near-side accumulation | not scoreable — **registered in advance as expected to fail; it failed** |

**The three numbers that matter.**

1. **Ordering: right.** On one real tissue section the engine gives half-distances
   of doxorubicin 45.2 = epirubicin 45.2 > daunorubicin 38.9 > mitoxantrone
   31.5 µm. That is Kyle's measured order exactly, tie included — four compounds,
   three distinct ranks. And it was **not** free: the shipping transport law, on
   its own, predicts the opposite order for the doxorubicin/daunorubicin pair
   (§2.6). The uptake sink overturns it (§5).

2. **Magnitude: ~3.7× too small.** The measurement separates doxorubicin from
   daunorubicin by **4.55×** (Table 1 fluorescence; 2.5× on the independent
   effect-based column). The engine separates them by **1.22×**. That misses the
   strict bar set by the measurement's own SDs (**[2.19, 10.83]**) and also the
   deliberately weak 3× bar.

3. **Absolute gradients: ~14× too shallow, and the error is in the OPTIMISTIC
   direction.** Daunorubicin's modelled near:far exposure ratio across a 140 µm
   slab is **3.65** against a measured **50 ± 15**. Even doxorubicin — the drug α
   is fit to — comes out at **3.00** against **11 ± 5**.

   **The model makes the drug look like it penetrates better than it does.** Say
   it in those words, because on a product whose new claim is *"here is where the
   drug does NOT reach"*, an error that flatters reach is the worst possible sign.
   A map that is too shallow under-reports the size of the untreated reservoir —
   it will show a tumour as covered when the measurement says it is not. If this
   error had gone the other way it would merely make the product pessimistic and
   therefore safe to ship with a caveat. It does not go the other way.

**And the finding that most constrains the product claim (C2).** In the clean
single-vessel slab the two compounds' nuisance intervals are **disjoint**
(doxorubicin 2.87–3.08, daunorubicin 3.50–3.76): the *physics* separates them
robustly. On the **reconstructed tissue map** they are **not** — doxorubicin spans
39.4–63.1 µm and daunorubicin 32.7–53.9 µm across the registered nuisance grid,
and those overlap heavily. **The ordering in §7 is reproducible only if every
analysis choice is held fixed. Allow the grid, the mask, the bin width or the
distance metric to move, and the compound difference disappears into the
method.**

**Arm A is clean and arm B is not — and that split is the actionable half of this
file.** The transport law separates these compounds robustly; the *reconstruction
readout* does not separate them on the section. **That localises the problem to
the spatial pipeline rather than to the physics**, which is a far more tractable
thing to fix: it is a question about grids, masks, binning and distance metrics,
not about diffusivity laws or sinks. Anyone reading this as "the transport model
is unreliable" has read it wrong. The transport model is *compressed* (§6, ~3.7×)
and *too shallow* (~14×); it is the map's spatial reduction that is *noisy*.

**One readout of this pipeline was separately tested and DOES clear its noise
floor.** `results/delivery_spread_noise.md` puts the panel's `P90/P10` delivery
spread through the same class of perturbation and finds temozolomide
(1.969–2.171) and doxorubicin (3.770–3.996) **disjoint with a gap ~7× wider than
either compound's own spread**. That metric is a percentile ratio over a fixed
cell population, so the distance-binning and distance-metric choices that break
the half-distance readout cannot touch it. **The C2 failure here is a property of
the distance-binned spatial readout specifically, not of every delivery number
this engine produces** — and it should not be quoted as though it were.

**What this means for the repositioning.** The claim *"within one map, this region
got N× less than that region"* is **not** validated by this test in its
quantitative form — the one out-of-sample gradient the engine predicts is ~14×
too shallow, so an N read off the map is not trustworthy. What survives is
weaker and should be stated as such: the map gets the *shape* qualitatively right
(exposure falls with distance from vessels, and it falls faster for compounds
measured to penetrate less), and it clears the source-null control, so it is
measuring transport rather than redrawing the vasculature. That is a real result
and it is worth having. It is not "the map is right."

**Nothing here was fitted, softened or re-scoped after the fact.** The bar was
frozen in `results/sweep/kyle2004_penetration_prereg.json` before any field was
solved, the run asserts live code against it, and the two criteria that were
registered as *expected to fail* (C7) and *not scoreable* (C6) are reported
exactly as registered.

---

## 1. What was being tested

Codon's surviving claim is a **map**: given one real tumour tissue section and one
compound, where the drug reaches and how much of the tumour it never reaches. The
claim is defensible because within a single map the uncertain per-compound
absolutes cancel — "this region got 10× less than that region" needs no absolute
calibration.

That claim has never been checked. `results/two_route_promotion.md` §8.9 states
the gap precisely:

> *"Nothing here says anything about a solved concentration field... no public
> dataset measures a spatially resolved free-interstitial concentration field in
> tumour tissue that this engine could be scored against."*

The engine's transport law has been scored **as a PDE coefficient** — R² = 0.718
against 71 measured multicellular-layer diffusion coefficients. It has never been
scored **as a solved spatial field**. This file is that test.

The instrument is a compound pair, and the logic matters: a measured
compound-to-compound difference in gradient steepness is the only public data
that constrains gradient *shape*, which is what the map claims. This is not a
ranking exercise — no compound is being ranked by predicted efficacy. It asks
whether a model of transport reproduces a measured transport difference.

---

## 2. Choosing the pair — and the two sources that were rejected

This section is longer than the result, deliberately. Choosing the pair was the
hardest part of the work and two apparently strong candidates failed on criteria
that only became visible on reading the Methods.

**Two gates, applied to every candidate:**

1. **Dose parity.** Were the compounds compared at matched exposure, or was one
   boosted for detectability? Read the Methods, not the Abstract.
2. **Is the measured quantity actually transport?** Not transport-plus-activation,
   not a summed signal from two chemical species.

### 2.1 REJECTED — Trédan 2009, mitoxantrone vs AQ4N (PMID 19176397)

The most attractive candidate on paper. Mitoxantrone and AQ4N are **exact
isomers** (both C22H28N4O6, MW 444.5), measured head-to-head in one laboratory by
quantitative immunofluorescence against distance from blood vessels *and* through
multilayered cultures, with the direction stated verbatim in the abstract.

**It fails both gates.**

* **Gate 1.** Multicellular culture: *"1 mg/mL for mitoxantrone and 10 mg/mL for
  AQ4N"*. In vivo: *"40 mg/kg and 200 mg/kg, respectively, to facilitate
  detection"*. A 5–10× dose asymmetry, chosen for detectability. No transport
  comparison survives that.
* **Gate 2.** The measured signal is *"AQ4N/AQ4 (which have similar fluorescent
  properties and cannot be separated)"* — the **sum of prodrug and active
  metabolite** — and at 48 h it concentrates in hypoxic regions with the
  discussion crediting stable DNA intercalation. Transport, bioreductive
  activation and binding are inseparable in that number. **A transport model
  cannot be scored against a quantity that is not transport.**
* It also tabulates no penetration number of any kind.

### 2.2 REJECTED — Grantab 2006, Table 2 (PMID 16424039)

Properly tabulated, four drugs, standard deviations, n = 3–5 — everything the
noise floor needs. It fails gate 1, and the paper says so itself:

> *"The initial concentration of drugs in compartment 1 was as follows: 25 μmol/L
> paclitaxel, 5 μmol/L doxorubicin, 10 μmol/L methotrexate, and 77 μmol/L 5-FU.
> These concentrations approximate those that can be achieved in serum after in
> vivo administration **and permit sensitive detection of the drug in compartment
> 2**"*

A **15.4-fold** spread across the compared compounds, explicitly driven by
detectability. Its dose-matched axis is **cell line, not drug** — the
epithelioid-vs-round contrast at fixed drug and dose is clean and is the paper's
actual claim. It is the right dataset for a different experiment.

### 2.3 REJECTED — Tunggal 1999 (PMID 10389947)

The most-cited such comparison (doxorubicin and mitoxantrone at 5–10 % of a
cell-free control, methotrexate ~30 %, 5-FU ~40 %). It fails gate 1 (≈9 / 10 / 50
/ ≈78 µM, an 8.6-fold spread with no rationale, and the methotrexate
concentration is self-contradictory between the Methods and the Fig. 4 legend),
and its headline pair fails gate 2 for a subtler reason: **doxorubicin was counted
as ¹⁴C radiolabel (parent plus any labelled metabolite) while mitoxantrone was
quantified by HPLC at 600 nm (parent only)** — two different measurands presented
as one comparison. It also prints no numeric SD or SEM anywhere, so it cannot
supply a measurement-level noise floor.

### 2.4 SELECTED — Kyle 2004, Table 1 (PMID 15342419)

> Kyle AH, Huxham LA, Chiam ASJ, Sim DH, Minchinton AI. "Direct assessment of drug
> penetration into tissue using a novel application of three-dimensional cell
> culture." *Cancer Res* 2004;64(17):6304-9.
> [DOI](https://doi.org/10.1158/0008-5472.CAN-04-1099)

Retrieved first-hand from the publisher (`curl` with a browser User-Agent;
`WebFetch` 403s on aacrjournals.org). **Table 1**, verbatim header — *"Drug
exposure to the near side relative to the far side of MCCs as estimated from the
effect-based assay or drug derived fluorescence"*:

| drug | effect-based assay | drug fluorescence |
|---|---:|---:|
| DOX | ~12 | **11 ± 5 (SD)** |
| EPI | ~10 | **12 ± 3 (SD)** |
| DAU | ~30 | **50 ± 15 (SD)** |
| MIT | >30 | *(not measurable)* |

HCT-116 multilayered cell culture, **140 ± 15 µm (SD)** thick when frozen
immediately after exposure, 1 h exposure to one side with the other sealed,
near/far read over *"the first 30 μm of tissue"* on each side.

* **Gate 1 — PASSES, and this is the only one of the four sources that does.**
  Verbatim: *"the drug was then added to the growth medium at concentrations of
  0.3, 1, 3, and 10 μmol/L"* — the identical four-point series for all four
  drugs — and both Table 1 columns are anchored at the **same 10 μmol/L** point.
* **Gate 2 — PASSES on the fluorescence column.** *"Only DOX, EPI, and DAU, which
  all possess the same fluorophore, could be examined. Emission intensities of the
  three drugs in aqueous solution were within 10% of each other."* One species,
  one fluorophore, cross-calibrated to 10 %, spatially resolved.

### 2.5 The caveat the source carries, stated before any result

Kyle's ratio is **not pure diffusion**. It is diffusion against cellular
sequestration, and the authors say the second term dominates:

> *"DOX appears to accumulate the least over the 1-hour exposure with EPI and DAU
> reaching ∼1.5 and ∼7 times higher levels, respectively, on the near sides."*

> *"Differences in the effects of the drugs on the far sides of the cultures were
> found to be statistically indistinguishable, P > 0.05, except at the
> 10 μmol/L × hour exposure, P < 0.005 (ANOVA), where the effect of EPI was
> significantly greater than DOX."*

**So the 12-vs-30 spread is a near-side phenomenon, not a far-side one.** That is
admissible here only because the engine models the same composite — a transport
coefficient `D_eff` against a loss/uptake sink `k`, giving `L = sqrt(D_eff/k)` —
and it would not be admissible against a pure-diffusion model. It is turned into
two extra registered criteria (§9) rather than left as a footnote, and it is
disclosed everywhere the result is quoted.

### 2.6 Why the pair is a risky test rather than a rigged one

Daunorubicin is doxorubicin minus the 14-hydroxyl:

| | PubChem CID | formula | MW | HBD | XLogP |
|---|---:|---|---:|---:|---:|
| doxorubicin | 31703 | C27H29NO11 | 543.5 | 6 | 1.3 |
| **daunorubicin** | **30323** | **C27H29NO10** | **527.5** | **5** | **1.8** |

Daunorubicin is **more lipophilic with one fewer hydrogen-bond donor**. The
shipping two-route law opens the transcellular route as a function of exactly
those two descriptors, with a decrement of 1.34 log units per donor. **On its own
it therefore predicts daunorubicin's `D_eff` to be the larger — i.e. that
daunorubicin penetrates further. The measurement says the opposite, by 4.5×.**
That prediction, and the fact that only the uptake sink could overturn it, were
registered in the pre-registration before the run.

And the pair is nearly matched in mass (3 % apart), so the **superseded**
molecular-weight-only law of `results/diffusivity_falsification.md` is essentially
flat between them — registered in advance as predicting **1.01×** against a
measured 4.55×.

### 2.7 Out-of-sample status, checked rather than asserted

* **None of the four compounds is among the 71 Pruijn 2008 HT-29 compounds** that
  `results/transport_anchor_prereg.md` §0.1 declares **spent**. Checked against
  the table in `experiments/diffusivity_falsification.py`: 67 SN-numbered
  benzotriazine di-N-oxides plus four miscellaneous compounds, and **zero**
  occurrences of doxorubicin, epirubicin, daunorubicin or mitoxantrone.
* **Daunorubicin is new to this repository.** Nothing has ever been fit to it.
* **Doxorubicin is IN-SAMPLE and its row is not evidence.** `calibrate_alpha()`
  is fit on doxorubicin, pinning its penetration length to 64.92 µm by
  construction. `scorecard.py::_run_perivascular_penetration` already declines to
  score it for exactly this reason and this file inherits that judgement. What is
  out-of-sample is the **ratio**, which is entirely a prediction about
  daunorubicin — the same structure the repository already accepts for
  mitoxantrone and topotecan in `pk_multidrug.run_penetration`.
* **Epirubicin** is out-of-sample as a fit but is numerically doxorubicin (§4).

---

## 3. The pre-registration

`results/sweep/kyle2004_penetration_prereg.json` was written before any
concentration field was solved. Every later run **asserts live code against it**
and raises on any drift, so a criterion cannot be edited once the outcome is
known. Every stage after `prereg` refuses to start until the file exists.

Registered before the run:

| | criterion | bar |
|---|---|---|
| **C0** | epimer control | the engine must **not** separate doxorubicin from epirubicin (rel. diff < 1e-6); the measurement does not resolve them either (11 ± 5 vs 12 ± 3) |
| **C1** | **direction — PRIMARY** | daunorubicin must be the **worse** penetrator in both arms |
| **C2** | noise floor | the two compounds' readout intervals must be **disjoint** over the full nuisance grid |
| **C3** | magnitude | strict: R_sim inside the measurement's own ±1 SD interval **[2.19, 10.83]**; loose: within 3× of the central **4.55** |
| **C4** | absolute ratios | each compound's modelled near:far ratio inside its Table 1 ±1 SD band — **only the daunorubicin row is evidence** |
| **C5** | is the map arm informative? | each compound's half-distance must beat the **source null** by ≥ 2× |
| **C6** | far-side disclosure | **NOT SCOREABLE** — the paper prints no numeric far-side values. Reported so a right answer cannot be reached by a mechanism the measurement excludes |
| **C7** | near-side accumulation | **registered as EXPECTED TO FAIL** before the run — see §9 |

The 3× loose bar was registered as *deliberately weak*: a 3× band around a 4.5×
effect admits R_sim below 1.0, i.e. the wrong direction. **C1 carries the verdict
and C3 qualifies it.**

---

## 4. The compounds, compiled through one path

Every compound goes through the identical pipeline — `from_smiles` →
`build_profile` → `molecule_to_drugparams` → `calibrated_drugparams(alpha)` →
`apply_barrier_to_drugparams` — with the same ECM reference and the same IC50, so
the only thing that differs between calls is the SMILES string. Three of the four
structures are **reused verbatim from the repository** so this file cannot
introduce a second, silently diverging definition of a molecule the codebase
already carries; only daunorubicin is new.

| drug | MW | logP | HBD | TPSA | `D_inter` | two-route hindrance | `P_mem` (µm/s) | `D_eff` | `k_loss` (1/s) | `L` (µm) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| doxorubicin | 543.52 | 0.001 | 6 | 206.1 | 184.749 | 0.1512 | 9.63e-07 | 27.926 | 5.035e-03 | **64.92** |
| epirubicin | 543.52 | 0.001 | 6 | 206.1 | 184.749 | 0.1512 | 9.63e-07 | 27.926 | 5.035e-03 | **64.92** |
| daunorubicin | 527.53 | 1.029 | 5 | 185.8 | 187.035 | 0.1512 | 2.27e-04 | 28.273 | 6.458e-03 | **57.68** |
| mitoxantrone | 444.49 | −0.139 | 8 | 163.2 | 200.639 | 0.1512 | 1.56e-09 | 30.328 | 9.868e-03 | **48.33** |

**C0 — the epimer control PASSES, and it is worth reading twice.** Doxorubicin and
epirubicin compile to *bitwise identical* transport parameters (relative
difference exactly `0.000e+00` in both `D_eff` and `k_loss`). The measurement does
not resolve them either (11 ± 5 vs 12 ± 3 overlap within 1 SD). The model and the
measurement agree that there is no difference here.

**But this is agreement by blindness, not by insight.** The engine's descriptors
cannot see stereochemistry, so it *could not* have separated them. The control
establishes that the readout does not manufacture differences; it establishes
nothing about predictive power, and it must not be quoted as if it did.

---

## 5. Result 1 — the decomposition, and the most important finding in this file

**The promoted two-route transport law is inert for this entire compound class.**

Read the `two-route hindrance` column above: **0.1512 for all four compounds**,
identical. That is the `P → 0` impermeable-cell Maxwell limit
`(1−φ)/(1+φ/2) = 0.15115` that `results/two_route_promotion.md` §3 identifies.
Membrane permeabilities span `1.6e-09` to `2.3e-04` µm/s — six orders of magnitude
— and **none of them is large enough to open the transcellular route at all**.

The consequence is exact, not approximate:

```
D_eff fold (DAU/DOX) = 28.273 / 27.926 = 1.0124
D_inter fold (DAU/DOX) = 187.035 / 184.749 = 1.0124      <- identical
superseded MW-only fold                    = 1.0124      <- identical
```

**Every bit of the engine's diffusivity separation between these compounds comes
from molecular weight — that is, from the law `diffusivity_falsification.md`
refuted.** The two-route physics that was built to create spread at matched
molecular weight, scored at R² = 0.718, and flipped on as the shipping default
contributes *nothing* here.

So where does the correct direction come from? Entirely from the uptake sink:

```
L = sqrt(D_eff / k)
L fold (DAU/DOX) = sqrt(1.0124 / 1.2826) = 0.8885
                          ^^^^^^   ^^^^^^
                          D fold   k fold  (6.458e-3 / 5.035e-3)
```

The diffusivity term pushes daunorubicin the **wrong** way by 1.2 %; the sink
pushes it the **right** way by 28 %; the sink wins. The same decomposition holds
for mitoxantrone (`D` fold 1.086 the wrong way, `k` fold 1.960 the right way,
net `L` fold 0.744).

**This matters for how the result may be described.** `k = clearance + α·uptake_proxy`,
and `results/alpha_identifiability.md` §1 proves that of the 11 parameters in that
calibration exactly **one is identified** (`k_dox`), **nine are not** (*"zero data
constraint"*), and α itself is a gauge coordinate. The direction this experiment
recovers is therefore produced by a **structure-only prior that the repository's
own analysis says is unconstrained by any data**. It came out right. It was not
*earned* right, and no claim here should imply otherwise.

---

## 6. Result 2 — arm A, the MCC analogue

Arm B tells you the ordering. It cannot tell you the *magnitude*, because a
half-distance on a reconstructed section and a near:far ratio across a cultured
slab are different quantities. Arm A removes that gap by reproducing **Kyle's
geometry** and computing **Kyle's quantity**, so the comparison needs no
conversion and no assumption of exponentiality.

A uniform slab at the real section's mean ECM (0.0875), **140 µm** thick — the
measured MCC thickness — with a one-sided planar source and a sealed far side.
The engine's boundaries are no-flux (Neumann), which *is* the sealed far side.
Exposure is integrated over 1 h, matching the measurement, and the readout is the
mean exposure over the near 30 µm divided by the mean over the far 30 µm —
Kyle's own windows. The solver is the same frozen `simulate_drug_transport`; the
construction is the one `vessel_profile.isolated_vessel_profile` already uses.

**One structural difference, disclosed:** Kyle's near face is a stirred reservoir
at fixed concentration (a Dirichlet condition); the engine has no Dirichlet
inflow, so the exposed face is represented by the vascular source occupying the
first column. For a decaying profile the two produce the same shape; they are not
identical, and the mismatch biases the modelled ratio **downward** relative to a
clean fixed-concentration boundary. That bias is common to all four compounds, so
it cancels from the fold comparison (C3) and does **not** cancel from the
absolute comparison (C4) — which is why C3 carries the magnitude verdict.

**The repository's own propagation guard shaped this arm.**
`tests/test_propagation.py::test_no_undeclared_paths[transvascular_exchange]`
fired on the first version of this script, which called `simulate_drug_transport`
directly and so would have made the `transvascular_exchange` flag a silent no-op
on an undeclared path. Rather than widen an exemption or add a declaration, arm A
was **routed through `simulate_response`** — the sanctioned fix the test itself
names, and the better construction anyway, because arm A and arm B now share one
code path. The two routes agree to 0.08 % (2.9977 vs 3.00 for doxorubicin), so
nothing was bought or lost numerically. The guard did its job on a file written
by someone who had not thought about that flag.

### 6.1 The result

| compound | **engine near:far** | **measured (fluorescence)** | in band? | measured (effect-based) |
|---|---:|---:|:--:|---:|
| doxorubicin | **3.00** | 11 ± 5 | no | ~12 |
| epirubicin | **3.00** | 12 ± 3 | no | ~10 |
| **daunorubicin** | **3.65** | **50 ± 15** | **no** | ~30 |
| mitoxantrone | **5.15** | *(not measurable)* | — | >30 |

**C1 — direction. PASS.** Daunorubicin's modelled ratio (3.65) exceeds
doxorubicin's (3.00): the engine places it as the worse penetrator, which is what
the measurement says. Mitoxantrone is placed worst of all (5.15), consistent with
its `>30`. The registered risk did not materialise — the transport law's own
prediction was overturned by the sink (§5).

**C3 — magnitude. FAIL, on both bars.**

```
R_sim      = 3.6457 / 2.9976 = 1.216
R_measured = 50 / 11         = 4.545      (+-1 SD interval 2.19 - 10.83)
                                           effect-based column: 30/12 = 2.50
```

`R_sim = 1.22` is outside the measurement's own ±1 SD interval and outside a 3×
band around the central value. **The engine under-separates this pair by about
3.7×.** That is not a surprise in kind — `results/two_route_promotion.md` §8.1
registered that the model *"under-separates a congeneric series systematically —
by roughly 4×"* against the Pruijn data. **This is an independent laboratory, an
independent assay and an independent compound class landing on the same
shortfall.** The corroboration is worth more than the individual number.

**C4 — absolute ratios. FAIL.** Only the daunorubicin row is evidence
(doxorubicin is the α fit; epirubicin is numerically doxorubicin). Daunorubicin's
**3.65 against 35–65** is **~14× too shallow**.

Part of that is construction: the source-term boundary biases the modelled ratio
down relative to Kyle's stirred reservoir, and doxorubicin's own 3.00-against-11
shows the offset is present even where the model cannot be wrong about `L`. But
the bias is common to all four compounds, so it cannot explain C3 — **the fold
comparison, which cancels it, still fails by 3.7×.**

---

## 7. Result 3 — arm B, the map on real tissue

The product's actual output. All four compounds on **one** reconstructed section —
the 10x Xenium human GBM, hBrain 319-gene panel, at the grid-converged
`dx = 12 µm`, 424 × 675, 201,674 tissue pixels, 4,295 measured endothelial cells,
mean ECM 0.0875. **The tissue is built once and reused**, so it is byte-identical
between compounds and every difference is the molecule.

The readout is the same quantity the immunofluorescence experiments report:
exposure binned by distance from the nearest *measured endothelial cell*
(`vessel_profile.direct_vessel_distance_um`, unmodified), reduced to the
**model-free** empirical half-distance — the distance at which the binned mean
first falls to half its perivascular value, no functional form assumed.

| compound | empirical d½ (µm) | fitted `L` (µm) | measured near:far (Kyle Table 1) | measured rank |
|---|---:|---:|---|---|
| doxorubicin | **45.2** | 49.7 | 11 ± 5 (~12) | best (tied) |
| epirubicin | **45.2** | 49.7 | 12 ± 3 (~10) | best (tied) |
| daunorubicin | **38.9** | 43.2 | 50 ± 15 (~30) | worse |
| mitoxantrone | **31.5** | 35.5 | — (>30, effect-based) | worst |
| *SOURCE NULL — `vascular_density` itself, no drug physics at all* | *10.4* | — | — | — |

**The engine reproduces the measured rank order exactly, including the tie.**
Kyle's ordering is doxorubicin ≈ epirubicin (best), then daunorubicin, then
mitoxantrone. The map gives 45.2 = 45.2 > 38.9 > 31.5. Four compounds, three
distinct ranks plus one tie, all correct.

**C5 — the real-section arm is informative. PASS.** Every compound clears the
source null by more than the registered 2× margin (4.35×, 4.35×, 3.74×, 3.03×).
`results/vessel_profile.md` warned that on this section the drug field can
inherit the shape of the vascular source field; here it does not — the null
halves at 10.4 µm while the drugs halve at 31–45 µm, so transport is doing real
work and the ordering is not vessel geometry in disguise.

**But the separation is small.** The whole span from the best to the worst
penetrator is 45.2 → 31.5 µm, a **1.43×** spread in half-distance, and the
doxorubicin-to-daunorubicin gap that the measurement puts at 4.5× is **1.16×** on
this map. §6 makes that comparison on the measurement's own scale, where it is
apples-to-apples.

**Doxorubicin's 45.2 µm sits inside the measured 40–50 µm Primeau band
(PMID 16361566). That is not evidence and must not be quoted as agreement** —
α was fit to that very anchor, so the value could not have come out otherwise.

---

## 8. Result 4 — the noise floor

A difference is only real if it is bigger than the spread produced by things that
differ between the two runs **for reasons that have nothing to do with the
compound**. Two independent noise grids were registered.

**Arm A** — grid resolution `dx ∈ {1, 2, 4} µm` × exposure window
`∈ {1 h, 6 h}` = 6 cells per compound.

**Arm B** — tissue grid and mask repair (`dx12_repaired`, `dx12_native`,
`dx30_native`) × exposure vs final concentration × distance metric (primary
direct-endothelial-cell vs legacy percentile-mask) × distance bin width
`∈ {5, 10, 20} µm` = 36 cells per compound.

The registered bar is that the two compounds' readout intervals must be
**disjoint** — not merely that their central values differ.

| arm | readout | doxorubicin span | daunorubicin span | disjoint? |
|---|---|---|---|:--:|
| **A** — clean slab | near:far ratio | **2.87 – 3.08** | **3.50 – 3.76** | **YES** |
| **B** — real tissue map | empirical d½ (µm) | **39.4 – 63.1** | **32.7 – 53.9** | **NO** |

**C2 — FAIL, and the split between the two arms is the single most
decision-relevant result in this file.**

* **The physics separates the compounds robustly.** In the clean slab the
  intervals do not touch: across a 4× change in grid resolution and a 6× change in
  exposure window, every daunorubicin cell sits above every doxorubicin cell.
  The separation is a property of the model, not of a lucky setting.
* **The reconstructed-tissue map does not.** The same two compounds, on the same
  section, span 39.4–63.1 µm and 32.7–53.9 µm across the registered nuisance grid
  — 36 cells each, no compound ever failing to halve — and those ranges overlap
  across most of their length. The **central** values are correctly ordered
  (45.2 vs 38.9 µm at the shipped settings, §7). The **ordering is not robust to
  the analysis choices**: the tissue grid, whether the cell-derived mask is
  morphologically repaired, exposure versus final concentration, the
  distance-bin width, and direct-endothelial-cell versus percentile-mask distance
  together move each compound's half-distance by more than the compound difference
  itself.

**So the honest statement is: the map's compound-to-compound difference is inside
its own methodological noise.** A user who re-ran the same two compounds at
`dx = 30 µm` on the native mask with 20 µm bins could get the opposite ordering.
Nothing about this is hidden by the shipped defaults — it is only visible because
the nuisance grid was registered in advance and run.

---

## 9. Result 5 — the near side and the far side

§2.5 recorded the source's own caveat before any result was seen: Kyle's spread
is a **near-side** phenomenon, with the far sides statistically indistinguishable.
Two criteria were registered to stop the engine reaching a right answer by a
mechanism the measurement excludes.

### 9.1 C6 — the far side. NOT SCOREABLE, and it disagrees.

| | doxorubicin | epirubicin | daunorubicin | mitoxantrone |
|---|---:|---:|---:|---:|
| **engine far-side exposure, relative to doxorubicin** | 1.000 | 1.000 | **0.862** | **0.401** |

Engine far-side spread across the four compounds: **2.49-fold**. The measurement
finds the far sides *"statistically indistinguishable, P > 0.05"*.

**The engine therefore gets its (already too small) near:far ordering partly by a
route the measurement excludes.** It is not scoreable — the paper prints no
numeric far-side values, and no number will be read off a figure to manufacture a
bar — but it is a discrepancy and it points the same way as C3: what the engine
is short of is *near-side* steepness, and it partly compensates by depleting the
far side instead.

### 9.2 C7 — the near side. Registered in advance as expected to fail. It failed.

| | doxorubicin | epirubicin | daunorubicin |
|---|---:|---:|---:|
| **measured near-side accumulation, relative to doxorubicin** | 1.0 | ~1.5 | **~7** |
| **engine near-side exposure, relative to doxorubicin** | 1.0 | 1.0 | **1.05** |

The measurement's central mechanism is that daunorubicin accumulates **~7× more**
in the first cell layers. The engine produces **1.05×** — essentially nothing —
and it cannot do otherwise, for the reason registered before the run: it solves
for **free interstitial** drug and represents cellular uptake as a **loss**
(`k = clearance + α·uptake_proxy`), not as a bound reservoir. The retardation term
that would have made it a reservoir was built and **refuted** in
`results/intracellular_retardation.md` and ships OFF. A stronger sink *lowers*
the engine's near-side concentration where the measurement *raises* total tissue
drug.

**This does not invalidate C1 or C3**, and the reason was registered in advance: a
near:far *ratio* is large for a high-uptake compound under both descriptions — a
store inflates the numerator, a sink steepens the gradient. But it does explain
C3's shortfall mechanistically. **The engine is missing the very mechanism the
authors identify as producing the effect**, and it recovers only 1.22× of a 4.55×
difference because a sink is a weaker instrument for this than a store.

---

## 10. What this licenses saying

### 10.1 Established, and defensible

* **The engine reproduces the measured penetration ORDERING of four compounds on
  one real tissue section**, including the correct tie on a measured-identical
  epimer pair: 45.2 = 45.2 > 38.9 > 31.5 µm against Kyle's 11 ≈ 12 < 50 < (>30).
* **That ordering was not free.** The shipping two-route transport law, taken
  alone, predicts the *opposite* order for the doxorubicin/daunorubicin pair,
  because daunorubicin is more lipophilic with one fewer donor. This was
  registered before the run. The uptake sink overturns it.
* **The map arm is measuring transport, not vasculature.** Every compound clears
  the source-null half-distance (10.4 µm) by 3.0–4.4×, against a registered 2×
  bar.
* **The engine under-separates a congeneric series by ~3.7× on this dataset** —
  1.22× modelled against 4.55× measured — which independently corroborates the
  ~4× under-separation `results/two_route_promotion.md` §8.1 measured against a
  different laboratory, assay and compound class.
* **The engine's absolute perivascular gradient is ~14× too shallow** for the one
  compound where it is a genuine prediction (daunorubicin 3.65 against 50 ± 15).
  The error is in the optimistic direction.
* **On the reconstructed section the compound difference does not clear the
  nuisance spread**, while in the clean single-vessel geometry it does. The
  physics separates these compounds; the map's spatial readout does not separate
  them above its own methodological noise.
* **The promoted two-route transcellular term is inert for this whole compound
  class** — hindrance 0.1512 for all four, the impermeable-cell limit — so the
  engine's diffusivity separation here is entirely the superseded
  molecular-weight law (§5).

### 10.2 NOT established, and must not be claimed

* **That the map's absolute penetration depths are right.** α is a gauge fit on
  doxorubicin and `results/two_route_promotion.md` §8.3 records two external
  anchors that now disagree about where the engine sits. Doxorubicin's agreement
  with anything is circular.
* **That the two-route transport law is validated by this.** §5 shows it is
  *inert* for this compound class — the transcellular route never opens, the
  diffusivity separation is pure molecular weight, and the direction comes from
  the uptake sink instead. This test says nothing either way about the law that
  was promoted on the Pruijn data.
* **That the direction was earned.** It is produced by `α·uptake_proxy`, a
  structure-only prior that `results/alpha_identifiability.md` §1 proves is
  **unidentified** for every drug except doxorubicin (*"zero data constraint"*).
  It pointed the right way. That is not the same as being constrained to.
* **That this generalises.** n = 4 compounds, one structural class
  (anthracyclines and one anthracenedione), one measurement, one laboratory, one
  tissue section, one tumour type.
* **That Kyle's ratio is pure diffusion.** It is diffusion against cellular
  sequestration and the authors attribute the ranking mainly to the latter
  (§2.5). The engine reproduces the ratio through a *loss* term where the
  measurement's mechanism is a *store* (§9).
* **Anything about clinical response.** Out of scope, as everywhere else in this
  repository.

---

### 10.3 What would most improve this

In the order that would change a conclusion fastest:

1. **A second measured pair from a different chemical class.** The whole result
   rests on one measurement of one anthracycline series. `results/axis_inventory.md`
   §7's test applies: this axis is per-compound, and every per-compound axis
   swept in this repository has moved a claim.
2. **A bound-drug compartment.** §9.2 shows the engine is missing the mechanism
   the authors identify. `results/intracellular_retardation.md` refuted the
   retardation term *against a free-drug diffusion coefficient*, which was the
   right call for that anchor; Kyle's quantity is **total tissue drug**, which is
   a different measurand and does not inherit that refutation. Whether a store
   term that is wrong for Pruijn is right for Kyle is now an open, answerable
   question.
3. **Stabilising the map's spatial readout.** C2's arm-B failure is a property of
   the reduction (grid, mask, bin, metric), not of the physics — arm A is clean.
   That is the cheapest of the three to attack and it is the one blocking the
   product claim.

---

## 11. Reproducing this

```
./.venv/bin/python experiments/isomer_penetration_map.py --stage prereg
./.venv/bin/python experiments/isomer_penetration_map.py --stage compounds
./.venv/bin/python experiments/isomer_penetration_map.py --stage all
```

`--stage all` takes ~25 min: 4 solves on the real Xenium GBM section at
`dx = 12 µm`, 28 slab solves for arm A and its grid sweep, and 6 further
real-section solves across three tissue arms for the noise floor. It writes
`results/kyle2004_penetration_data.json` and
`results/kyle2004_penetration_map.png`, and it **refuses to start** unless
`results/sweep/kyle2004_penetration_prereg.json` exists and matches live code.

### 11.1 Suite state, stated precisely

**Nothing in this work modified any file outside `experiments/` and `results/`.**
`git status` confirms the only additions are
`experiments/isomer_penetration_map.py`, this file,
`results/sweep/kyle2004_penetration_prereg.json`,
`results/kyle2004_penetration_data.json` and
`results/kyle2004_penetration_map.png`. **No assertion was deleted or weakened
anywhere.**

A concurrent agent is removing the ranking feature from `app/` and
`tmesim/cli.py` while this ran, so the suite is not quiet:

| observation | owner |
|---|---|
| `tests/test_propagation.py` — **53 passed** | **mine.** `test_no_undeclared_paths[transvascular_exchange]` fired on the first version of this script and is the reason arm A routes through `simulate_response` (§6). **Fixed properly, in my own file, by routing rather than by declaring an exemption.** |
| `tests/test_app_reports.py` — collection error, `cannot import name 'DEFAULT_PANEL'` | the concurrent agent's in-flight edit to `app/server/report.py` |
| `tests/test_app_server.py::TestApiKeyAuth::test_the_gate_uses_the_same_path_the_router_dispatches_on` — failed in the full run, **passes in isolation** | the same agent editing `app/server/` mid-run; a file-changed-under-the-run artefact, not a real regression |

Full run excluding the collection error: **1165 passed, 1 failed**, that one being
the `app/server` row above. Re-run once `app/` settles.
