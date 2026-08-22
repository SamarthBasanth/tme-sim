# Two-route tissue diffusivity — the transcellular route the transport model never had

> # ⚠️ CORRECTED 2026-08-14 — THIS LAW NOW SHIPS. THE HEADER'S "NOTHING IS WIRED IN" PARAGRAPH IS FALSE.
>
> **This banner is the only edit made to this file.** No number, table or verdict
> below it has been changed or re-run. Everything below was true on
> **2026-08-13**, when this law was default-OFF.
>
> **What changed.** Commit `fd2a9e7` (2026-08-14, *"Flip two-route transport
> ON"*) made it the default. `compile_drug.build_profile` now calls
> `two_route.two_route_diffusivity` on every compile, passing `logp` and `hbd`.
> The promotion — its gating conditions, its verdict propagation, and the values
> that moved — is `results/two_route_promotion.md`, which supersedes this file
> wherever the two disagree about the shipping configuration.
>
> **Specifically false now:** the header's *"Nothing is wired into
> `compile_drug.py`, no frozen snapshot is touched, `DrugParams.diffusivity` is
> bit-identical, and `calibrate_alpha()` still returns **148.60815765445545**."*
> All four clauses are false of the shipping engine. `alpha` is
> **22.223477711172563**. The superseded value is still exactly recoverable by
> passing `two_route=False`, and that is asserted in `tests/test_two_route.py`.
>
> **What still stands:** the physics, the citation, the limit check, and above all
> the pre-registration discipline — nothing was fitted to Pruijn, and the
> pre-registered sweep file is untouched. §4.4's disclosure of the one parameter
> with no primary source remains the least comfortable and most important
> section here; `results/two_route_efflux.md` resolved it with a **sourced** value
> that scores slightly worse.
>
> **A qualification added later, which this file could not have carried.**
> `results/descriptor_falsification.md` finds that the published blind R² was
> computed with the **paper's own measured descriptors**, not the descriptors the
> shipping pipeline computes. Every R² quoted below inherits that qualification.
> Read it before quoting a number from this file.

**What this file is.** The physics that closes the defect
`results/diffusivity_falsification.md` established, its citation, the limit check
that proves the algebra, what it scores against 71 out-of-sample measurements,
and the propagation of the correction through every published breast verdict.

Driver: `experiments/two_route_transport.py`. Code changed:
`tmesim/molecular/two_route.py` (new) and one additive function in
`tmesim/molecular/ecm_mesh.py`. Nothing is wired into `compile_drug.py`, no
frozen snapshot is touched, `DrugParams.diffusivity` is bit-identical, and
`calibrate_alpha()` still returns **148.60815765445545** (asserted in
`tests/test_two_route.py`).

> **The rule this file was written under.** The Pruijn dataset is the cleanest
> out-of-sample transport anchor this project has. **Nothing was fitted to it.**
> Every parameter was frozen in `results/sweep/tworoute_prereg.json` before the
> model was scored, and the scoring run *asserts* the live code still matches that
> file (§2.3). The paper's own Eqn-1 QSAR was deliberately not adopted: it is
> fitted to HT29 and the paper itself shows its asymptotes moving up to 100-fold
> across cell lines. §4.4 discloses the one parameter with no primary source,
> together with the value that would maximise the score — and §4.4 is the least
> comfortable section in this file, so read it.

---

## 0. Bottom line

Scored against **Pruijn 2008** (doi:10.1071/CH08240), 71 compounds, HT29
multicellular layers, one method, one cell line. 68 of the 71 carry a logP₇.₄ in
Table S1 and are scoreable; the other three are excluded and named in §4.

| | shipping engine | + impermeable cells | **+ transcellular route** | measured |
|---|---:|---:|---:|---|
| spread across the set | 2.32× | 2.32× | **7.11×** | **29.7×** |
| correlation, log₁₀ **R²** | 0.022 | 0.022 | **0.725** | — |
| Spearman | +0.083 | +0.083 | **+0.842** | — |
| geometric-mean bias | 3.40× high | 0.51× low | **1.25× high** | — |
| within 3× of measurement | 38 % | 68 % | **96 %** | — |
| tirapazamine, µm²/s | 287.4 (7.1× over) | 43.4 | **47.9 (1.19×)** | **40.4** |
| spread at **identical MW 305**, n=4 | 1.00× | 1.00× | **2.69×** | **3.8×** |
| spread at **identical MW 263**, n=4 | 1.00× | 1.00× | **1.45×** | **3.7×** |

**The one number that is fully parameter-free.** The model has exactly one
constant with no primary source (§2.4). Set it to zero — logP only, no
hydrogen-bond term, nothing unsourced anywhere in the model — and **R² is still
0.345 with Spearman +0.703 and 6.4× of spread** (§4.4). That is the floor this
work stands on if the convention is not believed; everything above it is the
convention's contribution and is disclosed.

**The matched-MW answer is YES.** At fixed molecular weight, where the shipping
model is flat *by construction* and the impermeable-cell correction is flat
*by construction* (Maxwell's factor contains no solute property at all), the
two-route model produces 1.10–2.69× against a measured 1.4–3.8×. Across every
exact-MW pair in the set it orders **24 of 31 pairs correctly (77 %)**; a
molecular-weight-only model orders **0 of 31**, because it cannot distinguish
them at all.

**Blast radius (§5), both α arms, published `dx = 12 µm` grid, 12 solved points
in each of two correction schemes: exactly 1 of 12 verdicts moves — docetaxel
FAIL → PASS, and only in `frozen`.** In `recal`, the mode the codebase actually
runs, **no verdict moves**, and docetaxel moves *further outside* its band
(1.5 % → ≥3.7 % outside). Doxorubicin has 6 hydrogen-bond donors, its membrane is
shut, and its transcellular gain is **1.000×**, so the correction at the α
calibration anchor is numerically the cell-exclusion correction already analysed
in `results/cellular_exclusion.md` §5.5. What is new is the **differential**
across the panel — 4.94× — and it is almost entirely paclitaxel, whose row is the
model's most questionable output (§5.1).

**Merge recommendation: MERGE the physics, do NOT flip the default.** §7.

---

## 1. The defect, restated in one paragraph

`results/diffusivity_falsification.md` measured that `D_eff` is a function of
molecular weight and ECM fraction only, reproduces 8 % of a measured 30-fold
spread, correlates with the measurement at R² = 0.028, and is exactly flat where
the measurement moves 3.7× at identical MW. It diagnosed the cause as structural:
the two dominant predictors of the measurement — hydrogen-bond donor count
(R² = 0.305) and logP₇.₄ (R² = 0.290) — are computed by
`descriptors.compute_descriptors` and discarded before the transport
calculation. Pruijn's own stated mechanism is that hydrophilic compounds are
confined to the **paracellular** route between cells while lipophilic ones also
access the **transcellular** route through membranes. The engine has no membrane.

`results/cellular_exclusion.md` reached the same place from the opposite
direction — a measured 12-fold MCL hindrance against a modelled 1.77-fold — and
added cells as *impermeable* inclusions. That was the right zeroth order and it
is what this work builds on. But Maxwell's impermeable factor
`(1−φ)/(1+φ/2)` contains **no solute property whatsoever**, so it moves every
drug by the same 6.61× and cannot, even in principle, produce spread. That file
says so itself, in its own §6: *"It does not make the hindrance drug-specific …
The real per-drug variation … is untouched by this fix and remains the sharpest
outstanding falsification of the diffusivity model."*

This file closes that.

---

## 2. The physics

### 2.1 The law

`tmesim/molecular/two_route.py::two_route_factor` implements **Maxwell's
effective-medium result generalised to a permeable inclusion phase with an
interfacial (membrane) resistance**.

The tissue is a two-phase composite: an interstitial **matrix** of conductivity
`σ_m` — the ECM-hindered interstitial diffusivity the engine already computes —
containing spherical **cells** of radius `a` at volume fraction `φ`, of interior
conductivity `σ_c = K·D_cyto`, separated from the matrix by a membrane of
permeability `P`.

```
σ_c*     = σ_c / (1 + σ_c/(a·P))                    [interior + interfacial
                                                     resistances in series]
β        = (σ_c* − σ_m) / (σ_c* + 2σ_m)
D_eff/σ_m = (1 + 2βφ) / (1 − βφ)                    [Maxwell]
```

> Maxwell JC. *A Treatise on Electricity and Magnetism*, Vol. 1, Art. 314
> (Clarendon Press, 1873) — the effective conductivity of a dispersion of spheres
> of one conductivity in a matrix of another. Diffusion in a two-phase medium
> solves the same Laplace problem as conduction, so the result transfers
> unchanged; Hashin & Shtrikman (1962) state the transfer to **diffusivity**
> explicitly and that clause is already quoted in `ecm_mesh.cell_exclusion`.
>
> **Hasselman DPH, Johnson LF. "Effective Thermal Conductivity of Composites with
> Interfacial Thermal Barrier Resistance." *J Composite Materials*
> 1987;21(6):508–515.** [DOI](https://doi.org/10.1177/002199838702100602) — the
> generalisation to a finite interfacial conductance `h`. Their published Eqn 8 is
> ```
> k_eff/k_m = [2(k_p/k_m − k_p/(a·h) − 1)φ + k_p/k_m + 2k_p/(a·h) + 2]
>           / [ (1 − k_p/k_m + k_p/(a·h))φ + k_p/k_m + 2k_p/(a·h) + 2]
> ```
> with their `h` this module's `P`. `two_route.hasselman_johnson_hindrance`
> writes that form out **literally**, and the two forms are asserted equal to
> machine precision — see §2.2 [D]. The series-resistance substitution above is
> therefore not an analogy; it is an algebraic identity with a published result.

### 2.2 The limits, checked in code rather than in prose

Run `./.venv/bin/python experiments/two_route_transport.py --stage limits`.
Every line below is an assertion in `tests/test_two_route.py`, not a claim.

| limit | requirement | result |
|---|---|---|
| **[A] `P → 0`** | must reproduce `ecm_mesh.cell_exclusion(0.7892)["hindrance"]` | `0.15115445288971746` vs `0.1511544528897175`, **&#124;diff&#124; = 2.8e-17** |
| **[A′] `P → 0` at every φ** | same, for φ ∈ {0, 0.1, 0.4, 0.6, 0.7892, 0.9} | all < 1e-15 |
| **[B] `φ → 0`** | no cells ⇒ factor 1 | 1.0 exactly |
| **[C] zero contrast** | `D_cell = D_interstitium`, `P → ∞` ⇒ factor 1 | 1.0 exactly |
| **[D] Hasselman–Johnson Eqn 8** | two algebraic forms of one published result | max &#124;diff&#124; **2.2e-16** over `P` ∈ [1e-4, 1e4] |
| **[E] dilute Maxwell slope** | `1 − 3φ/2` survives the generalisation | slope 1.5 to 1e-3 |
| **[F] monotone, bounded** | ↑ in `P`, ↓ in `φ`, never above the `P → ∞` ceiling | asserted |

**[A] is the load-bearing one.** It means the two-route law is a *strict
generalisation* of the law `results/cellular_exclusion.md` already shipped, not a
competitor to it. Any conclusion that file reached is recovered exactly by
setting `P = 0`, and the new physics can only ever move a compound in one
direction: faster than paracellular-only, never slower.

**[E′] The ceiling is a registered, falsifiable prediction.** As `P → ∞` the
factor saturates at the Maxwell value for permeable inclusions of conductivity
`K·D_cyto`. At the probe values used in the self-test that ceiling is **5.31×**
the impermeable value. No compound, however lipophilic, can exceed it, and it is
set by the cytoplasm/interstitium conductivity contrast — a *measured* quantity
(Kao 1993), not a knob. This is why §0 reports a bounded 7.11× against a measured
29.7× and why that shortfall was predicted in the pre-registration before the
scoring run.

### 2.3 The parameters, and where each one comes from

All seven are frozen in `results/sweep/tworoute_prereg.json`. `--stage prereg`
asserts the live code still matches it and **fails the run if any has drifted**,
so "we tuned it afterwards" is not a silent option.

| parameter | value | source |
|---|---:|---|
| `cell_volume_fraction` φ | **0.7892** | **Reused, not re-derived** — `results/cellular_exclusion.md` §3.2, Delesse's principle on the repo's own Xenium breast section (167,780 cells, Σ`cell_area` = 3.2495e7 µm²). |
| `cell_radius_um` a | **7.8517** | Equivalent-circle radius of the **mean** segmented cell area on the same section, √(193.675/π). Derived from the repo's own data. |
| `cytoplasm_D_over_water` | **0.27** | Kao HP, Abney JR, Verkman AS, *J Cell Biol* 1993;120(1):175–184, PMID 8416987, [DOI](https://doi.org/10.1083/jcb.120.1.175). |
| `cell_partition_K` | **1.0** | **Not measured.** An assumption, labelled as one. |
| `membrane_viscosity_pa_s` | **1.0** | Adopted unchanged from `permeability.py::_ETA_MEMBRANE_PA_S`. |
| `bilayer_thickness_nm` | **4.0** | Adopted unchanged from `permeability.py`. |
| `logk_membrane_max` | **3.0** | Adopted unchanged from `permeability.py::_LOGP_KP_MAX`. |
| `hbd_logp_decrement` | **1.0** | **A DECLARED CONVENTION WITH NO PRIMARY SOURCE FOR ITS MAGNITUDE.** See §2.4 and §4.4. |

The Kao value is quoted from its abstract verbatim: *"the relative diffusion
coefficients for BCECF translation and rotation in cytoplasm were 0.27 ± 0.01
(SEM, n = 24, 23 °C) and 0.78 ± 0.03"*. It is a *mechanistic* number, not a
black box — the paper decomposes it into fluid-phase viscosity (0.78), transient
binding (0.81) and collisions with cell solids (0.40), whose product 0.25
reproduces the measurement.

Two of these deserve a note on **direction of bias**:

* φ = 0.7892 is an *upper* estimate (Xenium's expansion segmentation tiles the
  section), and the literature tumour band implies 0.40–0.72. §4.5 sweeps it.
* `a` inherits the same inflation, and a larger `a` makes the membrane resistance
  group `σ_c/(a·P)` **larger**, i.e. *less* transcellular transport. The bias is
  conservative for this model.

### 2.4 The membrane leg — and why it is not `permeability.membrane_permeability`

```
log₁₀ K_mem = min( logP − 1.0·n_HBD , 3.0 )
D_mem       = kB·T / (6π·η_mem·r)          η_mem = 1.0 Pa·s
P           = K_mem · D_mem / d            d = 4 nm
```

Every constant is `permeability.py`'s own, asserted identical in
`tests/test_two_route.py`, so **no new tunable membrane number is introduced**.
Only the *topology* differs, for a stated physical reason:

`permeability.membrane_permeability` puts a 30 µm **unstirred aqueous boundary
layer** in series with the bilayer. That is right for a solute arriving from a
stirred bath. Inside tissue there is no such layer — the aqueous resistance
around a cell **is** the interstitium, and the interstitium is already the matrix
phase of the effective-medium problem this permeability feeds. Reusing that
function would double-count the aqueous leg and, worse, cap every lipophilic
compound at the same boundary-layer-limited ~22 µm/s. That function's own
docstring records the resulting compression as a known defect: *"the top of the
ranking is compressed … within 15 % of each other"*. The compression is precisely
the discrimination this module exists to produce.

`permeability.py`'s **lower** clamp at logP = −3 is also deliberately not
adopted. It is documented there as a numerical guard, not physics, and a floor at
K = 1e-3 would give urea a bilayer permeability of 0.126 µm/s — forcing the
transcellular route *open* for exactly the hydrophilic compounds whose exclusion
from it is the mechanism under test.

**The `hbd_logp_decrement` is where the honesty lives.** Its *structure* is
forced by two primary results, both verified against their own abstracts:

> Walter A, Gutknecht J. "Permeability of small nonelectrolytes through lipid
> bilayer membranes." *J Membr Biol* 1986;90(3):207–217. PMID 3735402,
> [DOI](https://doi.org/10.1007/BF01870127). 22 solutes across egg-PC/decane
> bilayers: *"High correlations were observed with hexadecane and olive oil
> (r = 0.95 and 0.93), but not octanol and ether (r = 0.75 and 0.74)."* Octanol is
> the wrong partition solvent for a bilayer, because — unlike an alkane — wet
> octanol accepts hydrogen bonds. The correction from octanol toward an alkane is
> therefore dominated by the solute's hydrogen-bond **donor** capacity.

> Conradi RA, Hilgers AR, Ho NF, Burton PS. *Pharm Res* 1991;8(12):1453–1460,
> PMID 1808606, [DOI](https://doi.org/10.1023/a:1015825912542) and *Pharm Res*
> 1992;9(3):435–439, PMID 1614980,
> [DOI](https://doi.org/10.1023/a:1015867608405). Over peptides spanning logP
> −2.2 to +2.8, *"no correlation was found between the observed flux and the
> apparent lipophilicity"* while *"a strong correlation was found for the flux …
> and the total number of hydrogen bonds the peptide could potentially make with
> water"*; sequential N-methylation raised permeability with *"only slight
> difference in the octanol-water partition coefficient"*.

Neither source states a **per-donor decrement in log units** that this work could
verify against a primary text, and none was invented. **1.0 log unit per donor is
a round number fixed before the model was scored**, recorded as such in the
pre-registration, and labelled a convention everywhere it appears in the code.
§4.4 sweeps it and discloses the optimum. **The disclosure is not flattering and
is not hidden.**

---

## 3. What the model does to the drugs the project cares about

`--stage drugs`, at the repo default `ecm_reference = 0.4`:

| drug | MW | logP | HBD | log K_mem | P (µm/s) | engine `D` | paracellular | **two-route** | transcellular gain |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| temozolomide | 194.15 | −2.08 | 1 | −3.08 | 0.124 | 278.0 | 42.0 | **42.9** | 1.02× |
| cisplatin | 298.03 | +0.20 | 2 | −1.80 | 2.02 | 236.2 | 35.7 | **49.0** | **1.37×** |
| doxorubicin | 543.52 | 0.00 | 6 | −6.00 | 1.06e-4 | 184.7 | 27.9 | **27.9** | **1.000×** |
| tirapazamine | 178.15 | −0.92 | 1 | −1.92 | 1.86 | 287.3 | 43.4 | **55.9** | 1.29× |

**Doxorubicin's membrane is shut.** Six hydrogen-bond donors at logP 0 puts
`K_mem` at 1e-6 and `P` at 1e-4 µm/s, so its transcellular gain is 1.000× to
four decimals. That is a physically sensible result — doxorubicin is a
notoriously poorly-permeating, strongly tissue-binding anthracycline — and it has
a large consequence for §5: **the correction at the α calibration anchor is
numerically the cell-exclusion correction already analysed.**

**Tirapazamine, the compound the project registered as an anchor.** Compiled from
its own SMILES rather than from the Pruijn MW, so the "before" column is
bit-identical to `results/external_anchor_survey.md` §3:

| | `D` (µm²/s) | vs measured 40.4 (HT29, Pruijn Table S1) |
|---|---:|---|
| engine, shipping | 287.3 | **7.1× over** |
| paracellular only (impermeable cells) | 43.4 | 1.08× |
| **two-route** | **55.9** | **1.38×** |

In the Pruijn run itself (where TPZ is compiled from the table's MW 178 and the
table's measured logP₇.₄ = −0.38 rather than RDKit's Crippen −0.92) the same
compound lands at **47.9 µm²/s, 1.19× the measurement**, with a transcellular gain
of 1.10×.

> **A real caveat, not a footnote.** The engine's `logp` is RDKit's Crippen
> estimate of the **neutral** octanol/water logP; Pruijn's column is a **measured
> logP₇.₄**, a distribution coefficient at pH 7.4. They differ, and for cisplatin
> the gap is severe — Crippen gives +0.20 against an experimental ≈ −2.2, which
> `permeability.py` already records. The membrane leg inherits that error
> wholesale, so **cisplatin's 1.37× gain above is built on a logP that is roughly
> 2.4 log units too lipophilic**. This model does not fix descriptors and does not
> claim to.

---

## 4. Scored against Pruijn 2008 — once, blind

`--stage pruijn`. 68 of the 71 compounds carry a logP₇.₄ in Table S1; SN28875,
SN28874 and SN28871 have `np` in both the measured and calculated columns and are
excluded from the two-route scoring (they set neither extreme, so the measured
29.7× spread is unchanged). The measured values, their provenance and the
extraction validation are all in `experiments/diffusivity_falsification.py` and
are **imported from there, not re-typed**.

> Every "shipping engine" number in this section is recomputed on those **68**
> compounds, so it differs slightly from the published 71-compound figures in
> `results/diffusivity_falsification.md` (R² 0.022 here vs 0.028 there, spread
> 2.32× vs 2.40×, bias 3.40× vs 3.49×). The three dropped compounds are not the
> reason the engine fails, and dropping them makes the shipping model look
> marginally *worse*, not better — i.e. the comparison is not being staged.

### 4.1 Spread — the statistic the falsification turned on

| model | min | max | **fold** |
|---|---:|---:|---:|
| **measured `D_mcl`** | 12.8 | 380.0 | **29.7×** |
| engine, shipping (MW-only) | 185.4 | 430.8 | 2.32× |
| + impermeable cells (paracellular only) | 28.0 | 65.1 | 2.32× |
| **+ transcellular route (two-route)** | **34.6** | **246.4** | **7.11×** |
| *(the paper's own Eqn 1 — not adopted)* | 15.0 | 339.0 | 22.6× |

The engine reproduced 7.8 % of the measured spread. The two-route model
reproduces **24 %**. That is a **3.1× improvement and still a 4.2× shortfall**,
and the shortfall was predicted before the run (§2.2 [E′]): the transcellular
route saturates, so a bounded model cannot reproduce an unbounded-looking spread.

Note the middle row. **Adding impermeable cells changes the spread by exactly
zero** — it is a uniform multiplier — which is the cleanest possible confirmation
that the previous correction, however large, could not have addressed this defect.

### 4.2 Correlation with measurement

| model | Pearson r (log₁₀) | **R²** | Spearman |
|---|---:|---:|---:|
| engine, shipping | +0.149 | **0.022** | +0.083 |
| + impermeable cells | +0.149 | **0.022** | +0.083 |
| **two-route** | **+0.852** | **0.725** | **+0.842** |
| *(paper's own Eqn 1, fitted to this data)* | +0.918 | 0.843 | +0.897 |

R² goes **0.022 → 0.725** on a dataset no parameter has seen. The paper's own
QSAR reaches 0.843 — but it was *fitted to these 71 points*, so the honest
comparison is: a physics model with no free parameter fitted to this data
recovers 86 % of the R² of a regression fitted directly to it.

### 4.3 Absolute accuracy

| model | ratio min | median | max | within 3× | geometric bias |
|---|---:|---:|---:|---:|---:|
| engine, shipping | 0.76 | 3.65 | 17.84 | 38 % | **3.40× high** |
| + impermeable cells | 0.11 | 0.55 | 2.70 | 68 % | **0.51× low** |
| **two-route** | **0.41** | **1.19** | **3.98** | **96 %** | **1.25× high** |

The impermeable-cell correction **overshoots** — it lands 0.51× low, i.e. it
crosses the measurement rather than reaching it, which is the same overshoot
`results/cellular_exclusion.md` §5.3 reported against the direct `D_tissue`
measurements and flagged as diagnostic. Re-opening the transcellular route
recovers it: the bias goes 3.40 high → 0.51 low → **1.25 high**, and 65 of 68
compounds land within 3× of measurement.

**The worst six residuals**, so the failures are visible:

| compound | measured | two-route | ratio | logP₇.₄ | HBD |
|---|---:|---:|---:|---:|---:|
| SN28860 | 39.9 | 158.9 | **3.98×** | +2.15 | 2 |
| SN28141 | 38.2 | 121.5 | 3.18× | +1.53 | 2 |
| SN28565 | 12.8 | 40.3 | 3.15× | +0.15 | 2 |
| SN27140 | 54.1 | 146.3 | 2.71× | +1.84 | 2 |
| SN28433 | 13.9 | 37.5 | 2.70× | −2.62 | 2 |
| SN26434 | 304.0 | 125.6 | **0.41×** | +1.06 | 2 |

Five of the six are over-predictions of **lipophilic, high-MW** compounds — the
model opens the transcellular route for them and the measurement says it should
not have, or not that far. Every one of the six has exactly 2 HBD, so the
residual is not the H-bond term; it is the lipophilicity term at the top end,
where §2.4's clamp at logP 3 is the only thing holding it. That is the sharpest
outstanding defect in this model.

### 4.4 The matched-MW test — the kill shot, and the disclosure

**Can the model produce any spread at fixed molecular weight?**

| MW | n | measured range (µm²/s) | measured fold | engine | paracellular | **two-route** |
|---:|---:|---|---:|---:|---:|---:|
| 192 | 3 | 54.3 – 75.0 | 1.4× | 1.00× | 1.00× | **1.10×** |
| 213 | 3 | 35.4 – 87.0 | 2.5× | 1.00× | 1.00× | **1.35×** |
| **263** | **4** | **13.9 – 50.9** | **3.7×** | 1.00× | 1.00× | **1.45×** |
| 289 | 3 | 59.0 – 83.5 | 1.4× | 1.00× | 1.00× | **1.82×** |
| **305** | **4** | **45.2 – 171.0** | **3.8×** | 1.00× | 1.00× | **2.69×** |

**Yes.** At MW 305 the model produces 2.69× against a measured 3.8×; at MW 263,
1.45× against 3.7×. Over **all 31 exact-MW pairs in the set the model orders 24
correctly (77 %)**. A MW-only model orders 0 of 31 — not badly, but *not at all*,
because the two members of every pair are the same number.

At MW 289 the model produces 1.82× where the measurement moves only 1.4× — it
over-spreads that tie. Both directions of error occur.

#### The disclosure sweep, and why it is uncomfortable

`hbd_logp_decrement` is the one parameter with no primary source for its
magnitude. If it had been tuned to this dataset, it would sit at the value that
maximises R². Here is the sweep:

| `hbd_decrement` | spread | **R²** | Spearman | geo bias | within 3× |
|---:|---:|---:|---:|---:|---:|
| 0.0 *(logP only — no convention at all)* | 6.4× | 0.345 | +0.703 | 2.20× | 66 % |
| 0.5 | 6.5× | 0.553 | +0.833 | 1.74× | 85 % |
| **1.0 — ADOPTED** | 7.1× | **0.725** | +0.842 | 1.25× | 96 % |
| 1.5 | 8.4× | 0.711 | +0.819 | 0.96× | 97 % |
| 2.0 | 8.7× | 0.638 | +0.802 | 0.81× | 96 % |
| 3.0 | 8.8× | 0.426 | +0.596 | 0.68× | 84 % |

**The adopted value is the sweep optimum.** That is the single most attackable
fact in this file and it is stated here rather than buried. Three things are
offered against the obvious suspicion, and none of them is proof:

1. The value was written into `results/sweep/tworoute_prereg.json` **before** the
   scoring run, and `--stage prereg` asserts the live code still matches it. A
   reader has to take the ordering on trust; the file is at least a fixed artefact
   that any later tuning would visibly change.
2. **1.0 is a round number.** A value chosen to maximise R² would land somewhere
   like 1.15, not exactly on a sweep grid point that was also the obvious default.
3. **The optimum is shallow and it is not the optimum on every metric.** R² is
   0.725 at 1.0 and 0.711 at 1.5 — a 2 % difference. On geometric bias 1.5 is
   *better* (0.96× vs 1.25×) and on within-3× 1.5 is better (97 % vs 96 %). A
   tuner optimising accuracy rather than correlation would have picked 1.5.

**The claim that does not depend on this parameter at all** is the top row:
with `hbd_decrement = 0` the model has **no unsourced constant anywhere**, and R²
still goes **0.022 → 0.345** with a Spearman of **+0.703** and 6.4× of spread.
That is the number to quote if the convention is not believed, and it is still a
16× improvement in explained variance.

### 4.5 Sensitivity to φ_cell — the widest-uncertainty input

| φ_cell | spread | **R²** | geo bias |
|---:|---:|---:|---:|
| 0.40 (bottom of the literature band) | 2.7× | 0.746 | 2.26× |
| 0.50 | 3.2× | 0.763 | 2.00× |
| 0.60 | 4.0× | 0.757 | 1.74× |
| 0.72 (top of the literature band) | 5.6× | 0.738 | 1.43× |
| 0.7635 (200 µm mask) | 6.5× | 0.730 | 1.32× |
| **0.7892 (adopted)** | **7.1×** | **0.725** | **1.25×** |
| 0.85 (50 µm mask) | 9.3× | 0.712 | 1.08× |

**R² is essentially flat at 0.71–0.76 across the entire plausible range of φ.**
That matters more than it looks: `results/cellular_exclusion.md` §4.3 established
φ_cell as "the dominant uncertainty in the hindrance law", and it still dominates
the *absolute* answer (geo bias moves 2.26× → 1.08×) — but it does **not** drive
the correlation. The correlation comes from the membrane leg, and the membrane
leg is what is new here. The adopted φ is also **not** the one that minimises
bias (0.85 would be) nor the one that maximises R² (0.50 would be).

---

## 5. Blast radius

### 5.1 How the correction reaches the engine — and why it is not one multiplier

Unlike the cell-exclusion term, the two-route factor is **compound-specific**.
Per-drug scales on the published breast panel, `--stage blast`:

| drug | MW | logP | HBD | log K_mem | P (µm/s) | **two-route scale** | fold | **vs doxorubicin** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| doxorubicin | 543.5 | 0.00 | 6 | −6.00 | 1.06e-4 | **0.1512** | 6.62× | **1.000×** |
| epirubicin | 543.5 | 0.00 | 6 | −6.00 | 1.06e-4 | 0.1512 | 6.62× | 1.000× |
| 5-fluorouracil | 130.1 | −0.80 | 2 | −2.80 | 0.270 | 0.1571 | 6.37× | 1.039× |
| cisplatin | 298.0 | +0.20 | 2 | −1.80 | 2.02 | 0.2077 | 4.82× | 1.374× |
| docetaxel | 807.9 | +3.26 | 5 | −1.74 | 1.68 | 0.2214 | 4.52× | 1.465× |
| **paclitaxel** | 853.9 | **+3.74** | 4 | −0.26 | **49.4** | **0.7466** | 1.34× | **4.939×** |

**spread of the correction across the panel: 4.94×**, and it is almost entirely
one drug.

> **Read the paclitaxel row as the model's most questionable output, not as a
> result.** Crippen puts paclitaxel at logP +3.74 with 4 donors, which opens its
> membrane almost completely (`P` = 49 µm/s, gain 4.94×, i.e. 93 % of the way to
> the §2.2 [E′] ceiling). Paclitaxel is in reality a P-glycoprotein substrate with
> notoriously poor tumour penetration. **This model has no efflux term**, and
> `efflux.py` exists in the package but is not in this chain. Passive
> physicochemistry says paclitaxel should cross membranes freely; measured biology
> says it does not get far into tissue. That contradiction is unresolved here, it
> is the largest single number in the blast radius, and it is the reason §7
> recommends against flipping the default.

`alpha` is fitted in closed form so that doxorubicin's `L = sqrt(D_eff/k_eff)`
reproduces a **measured** 64.92 µm penetration depth (Primeau 2005, PMID
16361566), which makes `alpha` nearly proportional to `D`.
`results/alpha_identifiability.md` establishes that `alpha` is a gauge coordinate:
it absorbs a change at the anchor and re-exports it as a global multiplier on
every other drug. **So the doxorubicin column is absorbed and the "vs
doxorubicin" column is what survives.**

And doxorubicin's transcellular gain is **1.000×** (§3). Therefore:

* the **absolute** perturbation at the anchor is numerically the cell-exclusion
  perturbation, so `results/cellular_exclusion.md` §5.2 applies verbatim —
  including the warning that the `frozen` arm breaks the measured Primeau anchor
  and that `recal` is the arm the codebase actually runs;
* what is **genuinely new** is the differential, and it is **not** small:
  paclitaxel moves 4.94× relative to the anchor and docetaxel 1.47×, while
  epirubicin and 5-fluorouracil move essentially not at all.

**One resolution hazard is *better* here than under cell exclusion, and that is
worth stating.** `cellular_exclusion.md` §5.4 recorded that in the `frozen` arm
its uniform 6.6× correction drove both taxanes below one grid cell per boundary
layer (`L/dx` 2.4 → 0.93 and 0.91 at `dx = 12 µm`), i.e. onto a grid that no
longer resolves the gradient producing the number. The two-route correction is
much gentler on precisely those two drugs, because their membranes are open:
`L` scales as √(scale), so paclitaxel goes to `L/dx` = 2.4 × √0.7466 = **2.07**
and docetaxel to 2.4 × √0.2214 = **1.13**, rather than to 0.93 and 0.91. Both are
still marginal; neither is below one.

### 5.2 Both arms, published `dx = 12 µm` grid

All 12 (mode × drug) points solved on the published `dx = 12 µm` grid, 458×628,
167,780 cells / 89,785 tumour, same `dt`, same 6 h window, same seed, same
`PDConfig`, same p99 bath-matching convention, same frozen literature bands.
`alpha`: published **148.6082** in `frozen`, re-solved in `recal`. Four of the
twelve rows (paclitaxel and 5-fluorouracil, both arms) are published **NOT
TESTABLE** — single-study bands — leaving **8 scored PASS/FAIL verdicts**.

**Table A — the doxorubicin-anchored *uniform* correction.** Directly comparable
with `results/cellular_exclusion.md` §5.5, and since doxorubicin's transcellular
gain is 1.000× it is numerically the same perturbation.

| mode | drug | ratio, shipping | ratio, corrected | verdict shipping | verdict corrected | **moved?** |
|---|---|---:|---:|---|---|:--:|
| frozen | doxorubicin | 0.2301 | 0.1403 | PASS | PASS | no |
| frozen | paclitaxel | 0.1503 | 0.0743 | NOT TESTABLE | NOT TESTABLE | no |
| **frozen** | **docetaxel** | **0.1484** | **0.0727** | **FAIL** | **PASS** | **YES** |
| frozen | epirubicin | 0.2301 | 0.1403 | PASS | PASS | no |
| frozen | 5-fluorouracil | 0.3379 | 0.2556 | NOT TESTABLE | NOT TESTABLE | no |
| frozen | cisplatin | 0.3085 | 0.2162 | PASS | PASS | no |
| recal | doxorubicin | 0.2301 | 0.2330 | PASS | PASS | no |
| recal | paclitaxel | 0.1503 | 0.1537 | NOT TESTABLE | NOT TESTABLE | no |
| **recal** | **docetaxel** | **0.1484** | **0.1517** | **FAIL** | **FAIL** | **no** |
| recal | epirubicin | 0.2301 | 0.2330 | PASS | PASS | no |
| recal | 5-fluorouracil | 0.3379 | 0.3353 | NOT TESTABLE | NOT TESTABLE | no |
| recal | cisplatin | 0.3085 | 0.3071 | PASS | PASS | no |

**Table B — each drug at *its own* two-route `D`, which is what a real wiring
would do.** The `frozen` rows are exact. The `recal` rows are **approximate for
the five non-anchor drugs**: the harness re-solves `alpha` at each swept `D`, so
it reports `alpha` = 32.68 for docetaxel where a real wiring would use
doxorubicin's 22.2241 everywhere. The `alpha` column is printed so the
approximation is visible rather than hidden.

| mode | drug | scale | ratio, shipping | ratio, corrected | verdict → verdict | α | **moved?** |
|---|---|---:|---:|---:|---|---:|:--:|
| frozen | doxorubicin | 0.1512 | 0.2301 | 0.1403 | PASS → PASS | 148.6082 | no |
| frozen | paclitaxel | 0.7466 | 0.1503 | 0.1366 | NOT TESTABLE | 148.6082 | no |
| **frozen** | **docetaxel** | 0.2214 | **0.1484** | **0.0864** | **FAIL → PASS** | 148.6082 | **YES** |
| frozen | epirubicin | 0.1512 | 0.2301 | 0.1403 | PASS → PASS | 148.6082 | no |
| frozen | 5-fluorouracil | 0.1571 | 0.3379 | 0.2574 | NOT TESTABLE | 148.6082 | no |
| frozen | cisplatin | 0.2077 | 0.3085 | 0.2321 | PASS → PASS | 148.6082 | no |
| recal | doxorubicin | 0.1512 | 0.2301 | 0.2330 | PASS → PASS | 22.2241 | no |
| recal | paclitaxel | 0.7466 | 0.1503 | 0.1503 | NOT TESTABLE | *110.8825* | no |
| **recal** | **docetaxel** | 0.2214 | **0.1484** | **0.1503** | **FAIL → FAIL** | *32.6790* | **no** |
| recal | epirubicin | 0.1512 | 0.2301 | 0.2330 | PASS → PASS | 22.2241 | no |
| recal | 5-fluorouracil | 0.1571 | 0.3379 | 0.3355 | NOT TESTABLE | *23.1069* | no |
| recal | cisplatin | 0.2077 | 0.3085 | 0.3076 | PASS → PASS | *30.6355* | no |

**Exactly 1 of 12 verdicts moves, in both tables: docetaxel FAIL → PASS, and only
in `frozen`.** In `recal` — the mode the codebase actually runs, because
`calibrate_alpha()` is called at run time — **no verdict moves in either table.**

**The docetaxel row is rigorous even where Table B's `recal` is approximate.**
The ratio is monotonically increasing in `D` at fixed `alpha` (Table A's frozen
column shows 0.1484 → 0.0727 as `D` falls), and docetaxel's own scale (0.2214) is
**higher** than the uniform one (0.1512). So at the true `alpha` = 22.2241 its
recal ratio must exceed Table A's 0.1517. Against a band top edge of 0.1463 that
means docetaxel goes from **1.5 % outside its band to at least 3.7 % outside**.
**The correction does not rescue the panel's only FAIL; it deepens it.** That is
the opposite of what a rigged correction produces, and it is the most useful line
in this section.

Two caveats carried over unchanged from `results/cellular_exclusion.md` §5.5,
because they still bind:

1. **The `recal` movements are smaller than the discretisation error underneath
   them.** `results/grid_convergence_breast.md` measures docetaxel's error at
   `dx = 12` as −3.4 % against its Richardson limit; the recal move here is
   +1.3 %. The *verdict* is safe (docetaxel is FAIL at every grid tested and this
   pushes it further out), but the *number* is not resolved to the precision
   quoted.
2. **The `frozen` movements are large and clean but their arm is not.** `frozen`
   silently breaks the measured Primeau penetration anchor (§5.1), and the one
   verdict that moves moves only there.

### 5.3 What this buys and does not buy

`results/dtissue_sensitivity.md` §8 established the structural point and it
applies verbatim: because `alpha` is fitted so that `D/k` reproduces a measured
penetration depth, the measurement constrains only the *ratio*. A uniform change
to `D` is absorbed. **What is not absorbed is the compound-to-compound
differential** — and that is precisely what this work adds and what the previous
correction had none of.

`results/cellular_exclusion.md` §6 also warns that threshold statistics such as
IC50 coverage move much harder than the p99-normalised ratio scored here, which
cancels a uniform exposure change exactly. That warning carries over unchanged and
those reports were not re-run here.

---

## 6. What this does and does not license saying

**Established, and defensible:**

* The transport model now has a transcellular route, built from an effective-medium
  result with a primary citation, whose impermeable limit reproduces the law
  already shipped to 2.8e-17.
* Against 71 same-method same-cell-line measurements it moves R² from 0.022 to
  **0.725**, geometric bias from 3.40× high to 1.25× high, and within-3×
  agreement from 38 % to 96 % — with no parameter fitted to that dataset.
* It produces spread at **identical molecular weight**, which no version of this
  model previously could, and orders 77 % of exact-MW pairs correctly.
* Tirapazamine goes from **7.1× over** its measured HT29 value to **1.19×**.
* The two properties the falsification identified as discarded — logP and
  hydrogen-bond donor count — now reach the diffusivity calculation, through a
  mechanism (membrane permeability) rather than through a regression.

**Not established, and must not be claimed:**

* **That the engine's diffusivity is now correct.** It reproduces 24 % of the
  measured spread, not 100 %; five of its six worst residuals are systematic
  over-predictions of lipophilic high-MW compounds; and the whole test is one
  cell line.
* **That the `hbd_logp_decrement` is a physical constant.** It is a convention
  with no source for its magnitude, and it happens to sit at the sweep optimum
  (§4.4). The parameter-free variant (`decrement = 0`) reaches R² = 0.345, and
  that is the defensible floor.
* **That this fixes predicted concentration fields.** This scores one PDE
  coefficient, not a solved field, and §5 shows `alpha` absorbs most of the
  absolute change at the anchor.
* **That the descriptors feeding it are right.** Crippen logP is ~2.4 log units
  wrong for cisplatin, and the model inherits that error linearly.
* Anything about clinical response. Out of scope.

**Known missing physics inside this model:**

1. **No intracellular partition or binding.** `K = 1` is an assumption. A drug
   that partitions into cells (`K > 1`) both conducts better *and* is retarded by
   the extra capacity; this model represents the first and not the second. For a
   steady-state flux measurement that is the right basis, but the Pruijn fit is to
   a *transient*, so a retardation term is a real omission.
2. **One φ_cell and one cell radius for the whole section.** A per-voxel
   cellularity field is derivable from the same `cell_area` column but lives in
   `tmesim/reconstruct/`, which this work does not own.
3. **Spherical cells, no junctions.** The paracellular route in a real epithelium
   is gated by tight junctions, which this model has no term for.
4. **The double-count `results/cellular_exclusion.md` §5.3 identified is not
   resolved.** `brownian.py`'s 300 nm obstacle field and `ALPHA_TORTUOSITY = 0.6`
   may be partly standing in for cell-scale obstruction; this work inherits that
   ambiguity, because it consumes the same interstitial `D` as its matrix
   conductivity.

---

## 7. Merge recommendation

**MERGE the physics. Do NOT flip the default in this change.**

**Merge now.** `tmesim/molecular/two_route.py` is new and is imported by nothing
in the shipping path. `ecm_mesh.tissue_hindrance_two_route` is purely additive;
no existing function's signature, default or return value changes, and
`tissue_hindrance` is untouched. Nothing calls either from `compile_drug.py`, so
`DrugParams.diffusivity` is bit-identical, every frozen `DRUG_LIBRARY` snapshot
is valid, `calibrate_alpha()` still returns 148.60815765445545, and every
published number stands. A test asserts each of those, including one that fails
if anyone imports `two_route` into `compile_drug` without re-running §5.

**Do not flip the default in this change**, for four reasons that are about
evidence:

1. **The one unsourced parameter sits at its own sweep optimum (§4.4).** Before
   this becomes a default, `hbd_logp_decrement` should be replaced by a sourced
   octanol→alkane conversion — the Abraham Δlog P analysis (*J Pharm Sci*
   1994;83:1085–1100) is the right primary source and was **not** obtainable
   within this run. Until it is, the defensible variant is `decrement = 0`
   (R² = 0.345), not the adopted one.
2. **Five of the six worst residuals are lipophilic over-predictions (§4.3).**
   The top of the logP range is held only by the clamp at 3, and that clamp was
   adopted from a different function for a different purpose.
3. **Every unresolved item in `results/cellular_exclusion.md` §7 is still
   unresolved**, and this model composes on top of all of them: the matrix
   conductivity it consumes is the possibly-double-counted interstitial `D`, and
   φ_cell is still one global constant above the literature band.
4. **The descriptor layer is wrong for at least one drug in the panel by ~2.4 log
   units** and is owned by another piece of work.

**The order of work that follows from this file:**

1. Obtain the Abraham 1994 Δlog P coefficients (or Seiler's I_H increments) and
   replace the declared convention with a sourced conversion. This is a
   literature question, not a compute question, and it is the single highest-value
   item because it converts §4.4's disclosure into a citation.
2. Add the retardation basis (`K > 1` with the accompanying capacity term) and
   report both bases, so the transient-fit commensurability question in §6.1 is
   answered rather than assumed.
3. Score the model against a **second** cell line. Pruijn's own paper reports the
   asymptote coefficients moving up to 100-fold across SiHa, HCT8-Ea and HCT8-Ra;
   a model claiming to be physics rather than a fit should be tested on that
   movement, and this run tested exactly one cell line.
4. Resolve the `brownian.py` / `ALPHA_TORTUOSITY` double-count
   (`cellular_exclusion.md` §7 item 2), which sets the matrix conductivity this
   law consumes.
5. Only then flip the default, in `recal` mode, and regenerate the coverage-type
   reports that `cellular_exclusion.md` §6 flags as unprotected by the α gauge.

---

## 8. Reproducing this

```
./.venv/bin/python experiments/two_route_transport.py --stage prereg   # <1 s
./.venv/bin/python experiments/two_route_transport.py --stage limits   # <1 s
./.venv/bin/python experiments/two_route_transport.py --stage drugs    # ~2 s
./.venv/bin/python experiments/two_route_transport.py --stage pruijn   # ~20 s
./.venv/bin/python experiments/two_route_transport.py --stage blast    # long
./.venv/bin/python -m tmesim.molecular.two_route                       # self-test
./.venv/bin/python -m pytest tests/test_two_route.py -q                # 44 tests
```
