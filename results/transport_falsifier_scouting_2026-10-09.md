# Transport-falsifier scouting, 2026-10-09 — measurements that could SCORE or BREAK the solver

**Status: SCOUTING PASS. Nothing here is a result, a test, or a pre-registration.**
No code was changed. No existing file in `results/` or `archive/` was edited. This
file is the only thing created. Nothing below is registered in
`tmesim/validation/scorecard.py`.

**Commission.** Find *"a measurement that could FALSIFY the transport solver"* —
drug or tracer concentration as a function of distance from vasculature,
extractable as numbers — by broadening off human GBM along five axes: (1) other
tumour types, (2) preclinical models, (3) non-tumour tissue, (4) engineered /
phantom systems with known diffusivity, (5) MSI repositories beyond METASPACE.
Explicitly instructed **not** to bring back input-side data (tissue sections,
single-cell references, panels), because `results/data_scouting_2026-10-09.md`
§5 established that side is saturated and a further find there has *negative*
value.

**Calibration guard, start of pass:**

```
$ .venv/bin/python -c "from tmesim.validation.calibration import calibrate_alpha; a=calibrate_alpha(); assert a==22.223477711172563, a; print('alpha OK',a)"
alpha OK 22.223477711172563
```

*(Closing guard in §8.)*

---

## 0. What I read before searching, and the one axis that was already CLOSED

Two files from earlier passes bound this one, and the second **removes a whole
class of candidate before any searching**:

* **`results/data_scouting_2026-10-09.md`** (today, earlier). Its verified
  negatives were taken as given and **not re-run**: human brain-tumour IFP in
  mmHg; a temozolomide perivascular µm profile; human CED concentration-vs-
  distance; drug PET in human glioma; and METASPACE (0 hits for temozolomide,
  glioblastoma, brain tumour, doxorubicin, paclitaxel, olaparib).
* **`results/verification_routes_2026-09-20.md`** — 14 routes, and it already
  occupies much of my assigned ground: R4/R5/R6 are multicellular-layer
  transport (axis 4), R7 is real-time iontophoresis in human glioma (axis 3),
  R8 is an in-human drug concentration-vs-distance profile digitised from a
  CC-BY figure, R10 is a rat tumour 3-D vascular network, R14 is CIVO. Nothing
  below duplicates a route there; where I touch one, I say so.

### 0.1 The correction that matters most in this file

`verification_routes_2026-09-20.md` R5 and R7 propose scoring
`ALPHA_TORTUOSITY = 0.6` against measured multicellular-layer and in-situ
tortuosity, describing the constant as *"currently constrained by **nothing**"*
and as *"carrying no citation at all"*. **Two other `results/` files contradict
that, and one of them forbids the comparison outright:**

1. **`results/uncertainty.md` line 46** already maps the constant to measurement:
   *"a ~2.5x diffusivity reduction at full ECM density … measured lambda in brain
   / tumour extracellular space spans ~1.3-2.2 (Sykova & Nicholson, Physiol Rev
   2008), i.e. a reduction of ~0.4-0.8 — which is the range sampled here, with
   the frozen 0.6 near its centre."* So the constant is **not** uncited and
   **not** unconstrained; it was sanity-checked against measured brain
   tortuosity.
2. **`results/reconstruction_validation_roadmap.md` line 98** goes further:
   *"The vascular KDE and ECM KDE fields are relative, max-normalised
   constructions with **no absolute physical scale** … **No published measurement
   can validate an amplitude that is defined only relative to the densest point
   in one section.** This is not a missing-data gap at all — it is a
   **construction** choice."*

I verified the construction in code. `tmesim/reconstruct/real_fields.py`
`_kde_on_grid` ends:

```
dmax = float(density.max())
if dmax > 0.0:
    density /= dmax
```

with the docstring *"divided by its own max so the densest region reads 1.0"*,
and `tmesim/core/microenvironment.py:112` documents `ecm_density` as
*"[0,1] **relative** extracellular-matrix density"*. So `ecm` attains exactly
1.0 at the densest pixel of **every** section regardless of that section's actual
matrix content, and `heldout_patient.md` §3.5 already shows the consequence
empirically — mean ECM 0.076 vs 0.457 across two panels of the *same* tumour
type, *"panel content, not tissue biology"*.

**Consequence for this pass, stated plainly.** Any candidate whose value is an
*absolute in-situ tortuosity* is dead on arrival against the shipped spatial
field, and `verification_routes` R7 is proposing a comparison that
`reconstruction_validation_roadmap` prohibits. That cross-document inconsistency
sits in the route list a future agent would work from, which is this project's
documented recurring defect (`CLAIM_INVENTORY.md`, the claim-surface-drift
thread). **I did not edit either file. It is recorded here and nowhere else.**

This also redirected the pass productively: what *survives* is a measurement of
how D changes **with a physically stated matrix quantity over a range** (a
slope), or a measurement of a term that is **not** relative — the membrane leg,
the clearance, the source. §2 is ranked on that basis.

### 0.2 The three legs the engine actually ships, read from code

Needed to say what any candidate scores. From `tmesim/molecular/compile_drug.py`
and the compiled provenance string for doxorubicin (read-only lookup, **not a
registered run**):

| leg | value for doxorubicin | where |
|---|---|---|
| Stokes–Einstein `free_D` | **351.31 µm²/s** at T = 310 K, **η = 1.2e-3 Pa·s** | `stokes_einstein.py` |
| 2-D Brownian micro-tortuosity at φ=0.4 | → 239.0 | `brownian.py` |
| **Ogston–Renkin fibre gate** | `mesh_h = 0.773` | `ecm_mesh.ogston_renkin_hindrance` |
| two-route cell term | φ_cell = 0.7892, **P_mem = 9.63e-07 µm/s**, factor 0.1512 | `two_route.py` |
| = `effective_D` | **27.93 µm²/s** (total hindrance 12.7×) | |
| spatial ECM factor, applied on top | `D_eff = D·(1 − 0.6·ecm)`, so **D_eff/D ∈ [0.40, 1.00]** | `physics/diffusion.py:75-88` |

**Two structural facts I did not find recorded anywhere in `results/`:**

* `ogston_renkin_hindrance` is called at `compile_drug.py:174` with **only the
  solute radius** — `r_f = 2.0 nm` and `φ_f = 0.05` are hard-wired defaults for
  **every compound and every tissue**. The fibre gate therefore cannot respond to
  a tissue's actual matrix content at all; all tissue-dependence lives in the
  max-normalised `(1 − 0.6·ecm)` factor. **The engine has no per-tissue,
  physically-scaled matrix hindrance anywhere.**
* The field named **`free_diffusivity` is not a free-solution value.** It carries
  η = 1.2e-3 Pa·s (interstitial fluid), against water at 37 °C ≈ 0.691e-3 Pa·s —
  a built-in **1.74×** penalty. Anyone comparing `free_D` to a published
  free-solution D will mis-score it by that factor, in a direction the field name
  conceals. This is the convention-not-travelling-with-the-number failure of
  `external_anchor_survey.md` §0.1, pre-loaded in a variable name. **It bit this
  pass** — see the sign flip in §2.2.

---

## 1. Ranking rule

Ranked by **capacity to make Vivum's position worse**, per the commission — not
by size, novelty or convenience. Two properties are weighted above everything:

* **Attributable** — would a disagreement name a *term*, or only say "something
  is wrong"? A non-attributable disagreement is nearly worthless here.
* **Reachable** — is the measured quantity inside the set of values the shipped
  configuration *can* emit? A measurement outside the reachable set is a
  refutation of the model's **form**, derivable without running anything, in the
  sense `docs/VERIFICATION_METHODOLOGY.md` §8 already has vocabulary for.

---

## 2. Ranked candidates

| # | candidate | axis | what it scores | extractable? | licence | refutes? |
|---:|---|---|---|---|---|---|
| **1** | **Shafieyan 2014** cartilage, D vs fluid volume fraction | 3/4 | **the LINEAR FORM of `ALPHA_TORTUOSITY`** | **yes — body text, verified verbatim** | free in PMC, no OA tag | **yes, form-level** |
| **2** | **Degerstedt 2024** DOX in defined gels + measured cell permeability | **4** | **`P_mem` (935× gap)**; corroborates the fibre gate | **yes — body text** | **CC BY 4.0** | **yes, attributable** |
| **3** | **Kyle 2014** camptothecins, 100 µm, one tissue | 1/2 | **the missing per-compound term** | partial — 3 ratios, free abstract | paywalled | **yes, attributable** |
| **4** | **Shi 2014** 6 solutes, radial gradients around single microvessels | 3 | `D_eff/D_free` vs the **0.40 floor** | **yes — free abstract** | paywalled | **yes, form-level** |
| **5** | **Dreher 2006** inert dextrans, k ≈ 0, permeability reported separately | 1/2 | **separates `D_eff` from the source term** | partial — abstract scalars | free PDF via Unpaywall | **yes, attributable** |
| **6** | **Zenodo 5776125** PEGPH20 + paclitaxel per-pixel MSI | 2/5 | relative prediction under a **matrix perturbation** | **yes — per-pixel `.xlsb`** | **CC BY 4.0** | yes, relative only |
| 7 | **Patel & Tannock 2009** DOX vs vessel distance, **Table 1** | 1/2 | λ, 8 profiles | **yes — a real table** | **CC BY 2.0** | partly — see §2.7 |
| 8 | **Smith 2017** eLife, radial point source, D pinned externally | 3 | **k alone**, and the convection term | figure-only | **CC BY** | yes, attributable to k |
| 9 | **Kutuzov 2018** single capillaries, radial profiles | 3 | **near-vessel `D_eff` and source symmetry** | figure-only | free in PMC | yes, structural |
| 10 | **Kyle 2007** taxanes, xenograft + multilayered culture | 1/2/4 | λ for two bound compounds | 2 scalars, abstract | paywalled | yes, not sharply |
| 11 | **Chary & Jain 1989** FRAP, measured interstitial velocity | 3 | **`v` magnitude AND sign** | abstract | free in PMC | yes, attributable to v |
| 12 | **Nugent & Jain 1984** size series, normal vs tumour | 1/3 | direction of the tumour-vs-normal effect | abstract, **one printed exponent unusable** | paywalled | see §5.2 |
| — | **Micromachines 2026** 3D microtissue paclitaxel | 4 | **nothing — no spatial profile** | n/a | CC BY | **rejected, §4** |

---

### 2.1 — #1 Shafieyan 2014: the one find that refutes a FORM, not a value

> Shafieyan Y, Khosravi N, Moeini M, Quinn TM. "Diffusion of MRI and CT Contrast
> Agents in Articular Cartilage under Static Compression." *Biophys J* 2014 Jul
> 15;**107(2):485–492**. **PMID 25028890**, **PMCID PMC4104049**, DOI
> 10.1016/j.bpj.2014.04.041.

**Provenance warning — read this before citing.** This paper was reported to me
as *"Pouran B, Arbabi V, Zadpoor AA, Weinans H … Biophys J 2014;107(9)"*. That
attribution is **wrong**: the authors are Shafieyan, Khosravi, Moeini & Quinn and
the pages are 107(2):485–492. I verified the correct record by PubMed `efetch`
and the PMCID, pages, DOI and open-access status (`isOpenAccess: N`,
`license: None`, `hasSuppl: N`) by the EuropePMC `core` record. **The numbers
below I then read myself, verbatim, from the PMC full text.** Had the misattribution
propagated, this repo would have shipped a citation to the wrong authors.

- **System:** bovine articular cartilage explant discs; sodium iodide, sodium
  diatrizoate, Gd-DTPA; static compression 0–60 %; radial transport into a
  well-mixed bath, modelled by the paper as 1-D radial diffusion.
- **What was measured, verbatim:** *"Diffusivity increased from 400 ± 123 to
  778 ± 111, 110 ± 19 to 155 ± 14, and 124 ± 28 to 195 ± 26 μm²/s when fluid
  volume fraction increased from 0.66 to 0.97, 0.78 to 0.94, and 0.80 to 0.95
  for sodium iodide, sodium diatrizoate, and Gd-DTPA, respectively (Fig. 6)."*
  And: *"Increasing static compression from 0% to 50% was associated with a
  decrease in the diffusivity of sodium iodide from 838 ± 111 μm²/s to
  412 ± 95 μm²/s … Diffusivity of Gd-DTPA decreased from 187 ± 31 to
  133 ± 32 μm²/s … 0% to 60%."* Sodium diatrizoate's trend was **not**
  significant (*p* = 0.061) — the paper says so and so do I.
- **Why this is the top find.** It is **D against a physically stated matrix
  quantity, over a range** — a *slope*, not a point. That is the one shape of
  measurement §0.1 leaves standing, because a slope can be compared without
  committing to an absolute scale for `ecm`.

**What the engine can and cannot reach.** Declaring `ecm ≡ 1 − φ_fluid` (solid
fraction), `D_eff/D = 1 − 0.6·ecm` gives, over each paper's own fluid-fraction
interval:

| solute | φ_fluid range | measured D ratio | engine @ 0.6 permits | too flat by | implied coefficient |
|---|---|---:|---:|---:|---:|
| sodium iodide | 0.66 → 0.97 | **1.945×** | 1.234× | **1.58×** | **1.50** |
| sodium diatrizoate | 0.78 → 0.94 | 1.409× | 1.111× | 1.27× | **1.64** |
| Gd-DTPA | 0.80 → 0.95 | 1.573× | 1.102× | 1.43× | **2.16** |

**All three implied coefficients exceed 1.0, and that is the finding.** For
`a > 1`, `D_eff = D·(1 − a·ecm)` goes **negative** above `ecm = 1/a`
(0.46–0.67), i.e. over much of the field, and the shipped floor
`max(1e-6, 1e-3·D)` would clamp it. **The measured slope is therefore not
reachable by re-tuning `ALPHA_TORTUOSITY` — it breaks the linear form.** Three
solutes, one paper, same direction.

- **Scores:** the functional form of the spatial ECM coupling, and only that.
- **Caveats, all of which must be in any pre-registration.** (i) The mapping
  `ecm ≡ 1 − φ_fluid` is **mine, not the paper's**, and §0.1 is precisely why a
  mapping onto a max-normalised field is a *declaration* that must be fixed
  before the residual is visible. (ii) Compression changes pore size and
  tortuosity *as well as* fluid fraction, so this is not a pure volume-fraction
  sweep; the paper itself reports GAG content as co-varying. (iii) Bovine
  cartilage is avascular and far denser than brain. **None of these rescues the
  form**, because the inadmissibility of `a > 1` is independent of the mapping's
  offset — it is a statement about slope.

---

### 2.2 — #2 Degerstedt 2024: the engineered system the commission asked for, CC BY, and a 935× gap

> Degerstedt O, O'Callaghan P, Clavero AL, Gråsjö J, Eriksson O, Sjögren E,
> Hansson P, Heindryckx F, Kreuger J, Lennernäs H. "Quantitative imaging of
> doxorubicin diffusion and cellular uptake in biomimetic gels with human liver
> tumor cells." *Drug Deliv Transl Res* 2024;**14(4):970–983**. **PMID 37824040**,
> **PMCID PMC10927899**, DOI 10.1007/s13346-023-01445-1.
> **Licence verified verbatim from the full-text XML: CC BY 4.0** (*"© The
> Author(s) 2023 … Creative Commons Attribution 4.0 International License"*).

This is axis 4 as the commission described it: *"a hydrogel … with a KNOWN
diffusivity and a measured concentration profile … score the solver against a
system whose ground truth is established independently."* I retrieved the full
text (153 KB) and read the numbers myself.

- **System:** a miniaturised chip; doxorubicin at 20 µM; **cell-free** gels of
  stated composition — 1 % w/v low-melting-point agarose (pH 7.41) and a
  "biomimetic cirrhotic" gel of **rat-tail type I collagen 5 mg/mL + fibrinogen**
  (pH 6.86 ± 0.07) — then the same gel with **Huh7 or HepG2 cells at 1 × 10⁶ (LD)
  or 2 × 10⁶ (HD) cells/mL**.
- **Measured, verbatim:** agarose **501 ± 77 µm²/s** (*n* = 8); cirrhotic gel in
  PBS **373 ± 108** (*n* = 7, *p* = 0.019); same gel in DMEM **256 ± 30**
  (*n* = 6, *p* = 0.028); + LD Huh7 **168 ± 16**; + LD HepG2 **310 ± 13**;
  + HD HepG2 **328 ± 74** (all *n* = 4). Apparent cell permeability
  **P_cell = 9.00 ± 0.74 × 10⁻⁴ µm/s** for HepG2. And
  *"Clinically relevant and quantifiable doxorubicin concentration gradients
  (1–20 µM) were established in the chip within one hour."*

**(a) The membrane leg: a 935× disagreement, and it points the same way as a known failure.**

| | doxorubicin cell permeability |
|---|---|
| **engine, shipped** | **9.63 × 10⁻⁷ µm/s** (`two_route`, compiled) |
| **Degerstedt, measured** | **9.00 × 10⁻⁴ µm/s** |
| ratio | **engine 935× LOWER** |

The paper independently notes a *published in-silico* value of
**220 × 10⁻⁴ µm/s** — 24× **above** their own measurement — and reports that
swapping their measured value in *"reduced the observed tumor intracellular DOX
concentrations by a factor ten (from approximately 5 to 0.5 µM)."* So the
literature spread for this one parameter covers ~4.4 orders of magnitude and the
engine sits at the extreme bottom.

**Why this is the most consequential number in the file.**
`tmesim/validation/calibration.py`'s own docstring says the structure-only priors
give a loss rate *"roughly **300x** faster"* than needed, that *"the missing
physics is … rapid cellular uptake plus high-capacity intracellular / DNA
binding"*, and that the logP permeability model *"spans ~5 orders of magnitude …
and mis-ranks doxorubicin"*. α exists to patch exactly that hole. An **external,
measured** doxorubicin permeability 935× above the engine's is independent
evidence that the missing ~300× lives in the **membrane/uptake leg** — a named,
shipped parameter — rather than in an unmeasurable "binding capacity". That
reframes α from a mechanistic constant into a patch for a measurably wrong
permeability, which is a materially **worse** position for the product and is
attributable.
*Honest discount:* their P_cell is an **apparent** permeability derived from
intracellular accumulation (`P = k_cell·V_cell/(A_cell·60)`), so it lumps active
transport, binding and partition, whereas the engine's `P_mem` is a *passive*
permeability. Doxorubicin is a P-gp substrate and a DNA intercalator, so these
are not strictly the same object. A 935× gap is far wider than that caveat
covers, but the comparison must be registered with it stated.

**(b) The fibre gate is CORROBORATED — recorded because a confirmation is also a result.**
Engine `mesh_h = 0.773` for doxorubicin's radius. Measured matrix-only contrast
(both gels in PBS, so buffer is matched): 373/501 = **0.745**. **Agreement within
3.8 %.** The honest reading is narrow: the gate is fixed at φ_f = 0.05 / r_f = 2 nm
*universally* (§0.2), so it happens to be right for a dilute collagen gel and is
applied unchanged to tumour tissue with no way to tell the difference.

**(c) The convention trap, and the sign flip it caused.** My first comparison was
engine `free_D` 351.31 vs measured 501 in agarose → "engine 1.43× **LOW**". That
is wrong, because `free_D` already carries η = 1.2e-3 (§0.2). Corrected to water
at 37 °C the engine's value is **609.8 µm²/s** → "engine 1.22× **HIGH**". The
disagreement **changes sign** under the viscosity convention. Nothing should be
claimed here without fixing the convention first; this is logged as a worked
example of the defect in §0.2.

**(d) The cell-density arm CANNOT test the packing term — stated so nobody tries.**
2 × 10⁶ cells/mL at ~15 µm diameter is a cell volume fraction of **~0.3 %**,
against the engine's **φ_cell = 0.7892** — roughly 200× apart. The engine's
Maxwell term predicts ~1.00× hindrance at 0.3 %, so the paper's own conclusion
(*"None of the observed effects of the added cells on DOX diffusion in the media
were statistically significant when compared to the same gel matrix without
cells"*) is **consistent with the engine, not adverse to it**. The HD arms only
got noisier (*"CV increasing from 9 to 41% (p = 0.014) and 4% to 23%
(p = 0.017)"*), and LD Huh7 (168) vs LD HepG2 (310) move in *opposite*
directions at identical density, which is line-specific confounding plus an
apparent-D readout biased by cellular uptake. **This arm is not a test of
`φ_cell`.** R4 (Grantab 2012) remains the only route that is.

- **Extractability:** the paper has **zero tables** (verified: `table-wrap`
  count 0); every value above is in abstract or body prose, which is extractable
  as text. The **µM concentration gradients are figure-only**, and data
  availability is *"available from the corresponding author on reasonable
  request"*. So the *coefficient* arm is runnable from the text today; the
  **field** arm — a measured C(x,t) in a matrix of known composition, which is
  the genuinely new thing — needs digitisation under the three conditions
  `verification_routes` R8 already sets, or one email.

---

### 2.3 — #3 Kyle 2014: one tissue, one distance, a 19× spread across compounds

> Kyle AH, Baker JH, Gandolfo MJ, Reinsberg SA, Minchinton AI. "Tissue
> penetration and activity of camptothecins in solid tumor xenografts." *Mol
> Cancer Ther* 2014 Nov;**13(11):2727–37**. **PMID 25143448**, DOI
> 10.1158/1535-7163.MCT-14-0475. No PMCID; **paywalled**.
> Verified by me via PubMed `efetch`.

- **System:** HCT116 and HT29 colorectal xenografts in mice; topotecan,
  irinotecan, irinophore C (liposomal irinotecan), and SN38; **native drug
  fluorescence mapped against the underlying vasculature**.
- **Verbatim:** *"in tissue 100 μm away from vessels, it reached 94% ± 5% of
  levels seen around blood vessels, whereas irinotecan and irinophore C were
  found to reach only 41% ± 10% and 5% ± 2%, respectively"*; SN38 *"penetrate
  tissue similarly to topotecan"*.
- **Why it ranks this high.** It holds **tissue, vasculature, model and method
  fixed** and varies only the compound, at one named distance. C(100)/C(0) spans
  **19×** across three agents. A single fitted α plus Stokes–Einstein size
  scaling cannot span that; the spread has to be carried by a per-compound `k`
  or a per-compound tissue-capacity term — and `results/brain_slice_anchor.md`
  established that the engine **has no per-compound tissue-capacity term**
  (implied `Vu,brain` = 1.000 for every compound, `GMFE_floor = 4.77×`). This is
  a pre-loaded, attributable, multi-point test of the exact term known to be
  missing.
- **Exclude irinophore C** from any small-molecule comparison — it is a
  nanoparticle.
- **A second sentence in the same abstract that cuts against the product's
  narrative, and should be quoted whenever this paper is:** *"all three agents
  were able to initially inhibit proliferation uniformly throughout the tumors,
  and it was their rate of washout … that correlated with activity."*
  Penetration depth did **not** determine effect. That is independent support for
  this project having withdrawn cross-compound ranking and absolute kill
  fraction, and it is a reason not to quietly re-acquire them.
- **Extractability:** three ratios at one distance, from the free abstract. Full
  profiles are figure-only behind an AACR paywall. So this is a 3-point test, not
  a profile test. Say so.

---

### 2.4 — #4 Shi 2014: six measured ratios, and four of them are below the engine's floor

> Shi L, Zeng M, Sun Y, Fu BM. "Quantification of blood-brain barrier solute
> permeability and brain transport by multiphoton microscopy." *J Biomech Eng*
> 2014 Mar;**136(3):031005**. **PMID 24193698**, DOI 10.1115/1.4025892.
> No PMCID; **paywalled**, abstract free. Verified by me via PubMed `efetch`,
> full abstract retrieved.

- **System:** laser-scanning multiphoton microscopy, **normal rat brain tissue
  100–250 µm below the pia**; sodium fluorescein (376 Da), FITC-dextran 4k/20k/
  40k/70k, Alexa-488 IgG (~160 kDa), injected via the ipsilateral carotid.
- **A true distance axis, verbatim:** *"P and Deff were determined from the rate
  of tissue solute accumulation and the **radial concentration gradient around
  individual microvessels** in the brain tissue."*
- **Verbatim numbers:** *"The mean Deff values were 33.2, 4.4, 1.3, 0.89, 0.59,
  and 0.47 × 10-7 cm2/s … The corresponding mean ratio of Deff to the free
  diffusion coefficient Dfree, Deff/Dfree, were **0.46, 0.19, 0.12, 0.12, 0.11,
  and 0.11**"*, and *"there is no significant difference in Deff/Dfree between
  solutes with molecular weights from 20,000 to 160,000 Da"*.
- **The reachability statement.** The shipped `D_eff/D` is confined to
  **[0.40, 1.00]** (§0.2). The small-solute ratio **0.46** is inside that band but
  requires `ecm ≈ 0.90`; the flagship section's **mean** ECM is 0.076
  (`heldout_patient.md`), giving `D_eff/D ≈ 0.954` — i.e. the engine is **~2.1×
  too fast in diffusivity** on its own flagship tissue against this measurement,
  which alone would make λ = √(D_eff/k) about **1.45× too long**. The repo
  independently reports λ **2.45× too long** (`kyle2004_penetration.md`), with
  `D_eff/k` ~6× too large. **This is the first external measurement that
  partitions that 6× between D_eff and k** — roughly 2.1× in D_eff, leaving ~2.9×
  in k — and attribution is the thing `verification_routes` says the current
  anchor set *"cannot make at all"*.
- **The four macromolecule ratios (0.19 → 0.11) are OUTSIDE the reachable set
  entirely.** No value of the constant and no `ecm ∈ [0,1]` reproduces 0.11.
  That is a form refutation — but an honest one: Vivum's compounds are small
  molecules, so this refutes the claim that `(1 − 0.6·ecm)` *is a tortuosity
  law*, not the shipped configuration's operating point. Both halves must be
  reported; collapsing them into one verdict would be an over-claim.
- **Cross-checks that make this safe to use.** Every value is 4.7–332 µm²/s,
  well under free-water self-diffusion (~2300 µm²/s). Its 70 kDa dextran value
  (5.9 µm²/s) sits within ~25 % of the independent Nicholson/Tao value
  (7.5 µm²/s) quoted in Smith 2017 (§2.8). And the control-cortex λ = 1.55 in the
  repo's held Vargová 2003 gives 1/λ² = **0.416**, agreeing with Shi's 0.46 for a
  small solute from a different method. Three independent readings converge.
- **An impeachment that must travel with it:** Kutuzov 2018 (§2.9) criticises
  Shi's radial-symmetry assumption directly. §2.9 states it.

---

### 2.5 — #5 Dreher 2006: the only candidate that separates D_eff from the source term

> Dreher MR, Liu W, Michelich CR, Dewhirst MW, Yuan F, Chilkoti A. "Tumor
> vascular permeability, accumulation, and penetration of macromolecular drug
> carriers." *J Natl Cancer Inst* 2006 Mar 1;**98(5):335–44**. **PMID 16507830**,
> DOI 10.1093/jnci/djj070. No PMCID; free PDF via Unpaywall per EuropePMC.
> Verified by me via PubMed `efetch`.

- **System:** dorsal skinfold window chamber, intravital laser-scanning confocal,
  fluorophore-dextrans **3.3 kDa → 2 MDa**, in solid tumour.
- **Verbatim:** permeability fell *"from 154 x 10(-7) cm/s, 95% confidence
  interval [CI] = 134 to 174 x 10(-7) cm/s, for 3.3-kDa dextran to 1.7 x 10(-7)
  cm/s, 95% CI = 0.7 to 2.6 x 10(-7) cm/s for 2-MDa dextran"*; and
  *"Dextrans of 3.3 and 10 kDa penetrated deeply (greater than 35 microm) and
  homogeneously into tumor tissue from the vessel wall. After a 30-minute period,
  a high concentration was observed only approximately 15 microm from the vessel
  wall"* for the 40–70 kDa range.
- **Why it is structurally special, and why it ranks above papers with better
  numbers:** the tracers are **inert**, so `k ≈ 0`, which removes the clearance
  term from the comparison; **and** the vascular permeability is reported
  *separately* from the penetration depth. So `D_eff` and the **vessel source
  term** can be constrained **independently of each other** — which none of the
  drug papers in this file permit, and which is exactly the attribution the
  repo's three failing solved-field rows cannot currently make.
- **Weaknesses, plainly:** transient (30 min, not steady state), macromolecular
  (so it constrains the tortuosity *form* rather than a small-molecule λ), and
  the depths are abstract scalars with the profiles figure-only.

---

### 2.6 — #6 Zenodo 5776125: a deposited, CC-BY, per-pixel drug map with a matrix ON/OFF perturbation

> *"PEGylated recombinant human hyaluronidase (PEGPH20) pre-treatment improves
> intra-tumour distribution and efficacy of paclitaxel in preclinical models"* —
> **DOI 10.5281/zenodo.5776125**. Paper PMC8434701, DOI
> 10.1186/s13046-021-02070-x. Istituto Mario Negri, Milan.
> **I verified the Zenodo record myself via the Zenodo API: title, DOI,
> `licence = cc-by-4.0`, publication date 2021-09-10, 37 files
> (e.g. `IALU182_1.xlsb` 20,657,046 B).**

- **What it is:** MALDI MSI of **paclitaxel** in mouse subcutaneous xenografts —
  SKOV3, SKOV3/HAS3 (hyaluronan-overexpressing) and BxPC3 — at a **100 × 100 µm**
  raster. Reported to me, and consistent with the file listing: the `.xlsb`
  workbooks carry per-pixel long tables plus tissue/ROI **masks** (sheet names
  `info, picco, fit, pixNr, plot, mask, maskT, maskS, ROI`). **That per-pixel
  structure was verified by the strand that found it, by downloading and binary-
  parsing one file; I verified the record and file manifest but did NOT parse a
  workbook myself.** Treat the internal structure as reported-not-verified.
- **Why it ranks here despite having no vessels:** it is a **controlled
  perturbation of the transport medium** — enzymatic hyaluronan degradation
  (PEGPH20) *plus* a genetic HA-overexpression pair — with per-pixel drug
  concentration read out in the same models. That is a *relative, within-compound*
  prediction under a measured matrix change, which is the regime
  `vivum-validation-state` concedes is the defensible one, and it is the closest
  deposited analogue to R4's packing-density logic.
- **The disqualifying caveat, and it is the structural finding of axis 5:**
  **no vessel or blood channel was co-imaged.** Reported to me as zero full-text
  hits for CD31, PECAM, meca-32, CD34, Hoechst 33342, lectin or "endotheli*",
  and zero hits for "distance". **As deposited it cannot be placed on a
  distance-from-vessel axis.** Sibling deposits **5158808** (niraparib/olaparib,
  CC BY, 11 files — verified by me) and **13375209** (pamiparib, CC BY, 38 files
  — verified by me, but the files are **PNG/TIFF + small `.xlsx`**, not per-pixel
  raw, so intensity→concentration is not recoverable) share the same absence.

---

### 2.7 — #7 Patel & Tannock 2009: the best extractable table, and the most compromised

> Patel KJ, Tannock IF. "The influence of P-glycoprotein expression and its
> inhibitors on the distribution of doxorubicin in breast tumors." *BMC Cancer*
> 2009 Oct 6;**9:356**. **PMID 19807929**, **PMCID PMC2770566**, DOI
> 10.1186/1471-2407-9-356. **CC BY 2.0.**
> Verified by me via PubMed `efetch`; **Table 1 read by me from the PMC HTML.**

**This is NOT the held Patel/Trédan/Tannock 2013 (PMID 23680920)** — different
year, journal, title and compound set. Confirmed distinct.

- **System:** EMT6 murine mammary and PgP-overexpressing AR1; MCF-7 human breast
  xenografts and PgP sublines; doxorubicin 25 mg/kg IV ± verapamil or PSC 833.
- **Table 1, as I read it** — DOX uptake at 10–20 / 50–60 / 110–120 µm from the
  nearest vessel, with vessel density and a fitted gradient:
  EMT6 wild type **29.4 ± 6.1 / 15.6 ± 1.9 / 11.3 ± 2.4**, gradient
  **−0.23 ± 0.08**, vessel density 4.7 ± 1.2; AR1 (PgP+) **17.1 ± 2.9**, vessel
  density 6.0 ± 1.4. Eight tumour × treatment groups in all.
- **Three honest discounts, which together drop it from the top tier.**
  (a) **Same laboratory, same autofluorescence method and the same drug as
  Primeau 2005**, the in-sample anchor α is solved to — so it tests
  tissue-to-tissue transfer, not the method, and the doxorubicin rows inherit the
  in-sample problem the scorecard already handles by declining doxorubicin's own
  penetration anchor. (b) The PgP rows add an **active-efflux** term, so a
  mismatch there is **not attributable**; only the wild-type/DOX-only rows are
  usable. (c) The paper's own ceiling: *"At distances greater than about 90 - 100
  μm from blood vessels estimates of doxorubicin fluorescence are subject to
  noise."* **If the engine's λ is much longer than 100 µm, this dataset
  physically cannot falsify it.** Given the engine's λ is reportedly 2.45× too
  long, that is a live risk and must be checked *before* any fit.
- One finding worth carrying regardless: *"there were steeper gradients … in
  wild-type tumors compared to PgP overexpressing tumors"* — efflux **flattens**
  the gradient, which is the opposite of the intuitive direction.

---

### 2.8 — #8 Smith 2017: a radial point source with D pinned from outside, so it constrains k alone

> Smith AJ, Yao X, Dix JA, Jin BJ, Verkman AS. "Test of the 'glymphatic'
> hypothesis demonstrates diffusive and aquaporin-4-independent solute transport
> in rodent brain parenchyma." *eLife* 2017 Aug 21;**6:e27679**. **PMID 28826498**,
> **PMCID PMC5578736**, DOI 10.7554/eLife.27679. **CC BY.**
> Identity verified by me via PubMed `efetch`; the figure-level detail below was
> reported to me and is **not independently verified**.

- Mouse/rat brain; fluorescent dextrans (10/70/2000 kDa) and Alexa-647
  ovalbumin; intracisternal and intraparenchymal injection. Three profile sets
  were reported: depth from the cortical surface; a **penetrating-arteriole**
  paravascular profile; and *"fluorescence as a function of radial distance from
  the center"* of an intraparenchymal injection with per-dextran 10 %-decay
  distances.
- **Why it is attributable to k specifically:** it uses **independently measured**
  D values from the Nicholson lineage (reported: 16, 51 and 7.5 µm²/s for
  ovalbumin, 10 kDa and 70 kDa dextran). With `D_eff` pinned externally, a
  measured decay distance constrains **`k` alone**. That is the cleanest
  λ = √(D_eff/k) test located in this pass, and it is CC BY so digitisation is
  licence-clean.
- It also bears on the **convection** term: the paper's conclusion is that
  transport is purely diffusive with no convective component over 60 min, and
  `PDEConfig.convection` ships **ON**.
- **Extractability: figure-only.** Reported, and worth repeating because it was
  checked specifically: **no data-availability statement and no source-data
  files**; the only additional file is a transparent reporting form.

---

### 2.9 — #9 Kutuzov 2018: two structural threats to the source term

> Kutuzov N, Flyvbjerg H, Lauritzen M. "Contributions of the glycocalyx,
> endothelium, and extravascular compartment to the blood-brain barrier." *PNAS*
> 2018 Oct 2;**115(40):E9429–E9438**. **PMID 30217895**, **PMCID PMC6176561**,
> DOI 10.1073/pnas.1802155115. Free in PMC; no OA licence tag reported.
> Identity verified by me via PubMed `efetch`; the two quoted findings were
> reported to me and are **not independently verified**.

Two reported findings, both aimed at the **vessel source term** rather than at a
coefficient:

1. *"In the extravascular compartment, co-localized with astrocyte endfeet,
   diffusion coefficients of the dyes were an order of magnitude lower than in
   the brain parenchyma."* The engine carries **one** `D_eff` field driven by
   `ecm`; a ~10×-lower-D shell in the first few micrometres around every vessel
   is exactly where the source sits, and a max-normalised stromal-cell KDE will
   not reproduce it. On a 12 µm grid it would also be sub-cell.
2. *"In 6 out of 22 selected vessels, the initial circular symmetry of the
   fluorescence distribution was preserved."* The engine's source assumes radially
   symmetric egress. Measured: **27 %**.

**It also impeaches §2.4.** The paper reportedly criticises Shi 2014's symmetry
assumption directly (*"their image data … showed an absence of the symmetry,
which challenges the accuracy of their estimates"*). Candidates 4 and 9 partly
impeach each other and **must be cited together or not at all.**

---

### 2.10 — #10–#12, briefly

- **#10 Kyle 2007** (PMID **17473214**, *Clin Cancer Res* 2007;13(9):2804–10,
  paywalled — verified by me). *"In xenografts, taxanes were mapped relative to
  blood vessels to obtain drug profiles as a function of distance from
  vasculature"*; *"little drug reaching further than 100 microm"*, paclitaxel
  *"up to 2-fold greater penetration than docetaxel"*. **Also has a multilayered-
  cell-culture arm** measuring depth into tissue, so it straddles axes 1/2 and 4.
  Note the near-miss: PMID 17473213 circulates for this paper and is **wrong**.
- **#11 Chary & Jain 1989** (PMID **2748592**, PMCID PMC297627, *PNAS*
  1989;86(14):5385–9 — verified by me, verbatim). *"The average fluid velocity
  was about 0.6 microns/s, and albumin diffusion coefficients were 5.8 +/- 1.3 x
  10(-7) cm2/s and 6.3 +/- 1.9 x 10(-7) cm2/s in normal and neoplastic tissues …
  The interstitial fluid flow, in general, was **directed into postcapillary
  venules**."* Three things: a measured **`v` = 0.6 µm/s**; D essentially
  **identical** in normal and neoplastic tissue (supporting this whole pass's
  premise that transport physics is not tumour-specific); and a flow field
  **converging on** vessels, whereas the engine treats segmented vessels as a
  source with no corresponding sink — so the **sign** of the advective term near
  a subset of vessels may be wrong. Attributable to `v`.
  *Dissent, recorded deliberately:* its D/D_free for albumin is ~0.7–0.95, i.e.
  nearly unhindered — the **opposite** direction from Shi and Fu, and it would
  make 0.6 look too steep. The literature spread for similar-sized solutes across
  tissues is ~6×, so **the constant cannot be validated by picking one study.**
- **#12 Nugent & Jain 1984** (PMID **6197161**, *Cancer Res* 1984;44(1):238–44,
  paywalled, abstract verified by me — see §5.2 for the defect). FITC-BSA and
  FITC-dextrans Mr 19,400–71,800 plus sodium fluorescein, in normal granulation
  tissue and VX2 carcinoma in a rabbit ear chamber. *"A one-dimensional diffusion
  model adequately described extravascular transport in both normal and tumor
  tissue"*; **"Macromolecular transport in tumor tissue was hindered to a lesser
  extent than in normal tissue"**. With Chary & Jain, that is **two measurements
  from the same lab saying tumour matrix is *more* permissive**, while the engine
  raises tortuosity where stromal/collagen density is high. Against that, Song
  2017 reports tumour tortuosity *rising* with progression. The size regimes
  differ (macromolecules vs Gd-DTPA), so this is an **unresolved sign question,
  not a refutation** — and the engine's direction is therefore not externally
  determined.

---

## 3. The structural finding of this pass

**The refutation-grade object the commission describes — a deposited drug raster
with a co-imaged vessel channel — does not exist in any public repository
checked.** Four deposited drug spatial maps were found and three verified by me
as CC-BY records; **not one co-imaged endothelium or blood.** Meanwhile the one
high-resolution tumour **vessel** dataset located (S-BIAD500 / idr0085, MF-HREM,
CC BY 4.0, ~2 µm, *"vascular networks in tumour xenograft models"* — reported,
not verified by me) has **no drug**.

The bottleneck is therefore **not** a scarcity of drug MSI deposits. It is that
**MSI drug papers do not stain for endothelium.** Of 15 papers verified against
primary records on axes 1–2, exactly **two** put distance-resolved numbers in a
table (§2.7 and Liang 2013), and **zero** deposited the underlying profiles.

That has a concrete consequence for any future Vivum claim: a
distance-from-vessel prediction tested against public data must either (a) use a
**non-vessel geometric axis** (spheroid edge, GI mucosal depth, cartilage/explant
surface — §2.1 is exactly this), or (b) test a **relative, perturbation-based**
prediction (§2.6), or (c) accept **one or two scalars from an abstract** per paper
rather than a profile. It confirms, from a third independent direction, the
conclusion of `external_anchor_survey.md` §5 that closing this gap *"may require
generating data, not finding it."*

---

## 4. Checked and REJECTED — so nobody re-chases the title

> Fariha R, Hamze J, Okoh OD, Rothkopf E, Tripathi A. "Integrated
> Experimental-Computational Framework for Drug Transport Quantification in 3D
> Microtissues™." *Micromachines (Basel)* 2026;**17(3):332**. **PMID 41900218**,
> **PMCID PMC13028606**, DOI 10.3390/mi17030332. **CC BY** (verified verbatim).

On its title — "drug transport quantification", "3D microtissues", with a
finite-element diffusion-reaction model — this is the ideal axis-4 candidate. **I
read the full text. It does not work:** it is an **LC-MS/MS bulk quantification
method paper**. Searching the full text returns **no diffusion coefficient**, no
*"distance from"*, and **no depth or spatial profile** of any kind; the compound
is paclitaxel measured as whole-microtissue content, and the FEM is used only for
*"mechanistic context"*. Data availability: *"Additional data can be made
available upon reasonable request."*
**Verdict: not a transport measurement. Do not re-chase on the title.**
(Same failure mode as Cribaro 2021 in `data_scouting_2026-10-09.md` §4.)

---

## 5. Traps and published defects caught — this is as valuable as the finds

### 5.1 A compilation table that must be discarded
**Namazi H, Kulish VV, Wong A, Nazeri S.** *BioMed Res Int* 2016;2016:8437247,
PMID 27376087, PMCID PMC4916326, CC BY. Its **Table 1** looks exactly like the
benchmark this project wants — "experimental penetration depth" for Dox-LTSL 78,
Doxil 34, cisplatin 20, docetaxel 100, paclitaxel 100, 5-FU 24, mitoxantrone 37,
Herceptin 150, topotecan 100 µm. **It is not usable as data.** Reported to me,
from the full text: no per-row citation, no per-row tumour model, species or
method; the table is cited collectively to refs 46–52, of which **two are books**
and two are nanoparticle/stent papers, and the only traceable rows are Kyle 2007,
Baker 2008 and Kyle 2014 — already listed above. **The cisplatin 20 µm, 5-FU
24 µm and mitoxantrone 37 µm rows have no retrievable primary source.** Same
failure mode as the Huber/ResearchGate numbers in `data_scouting_2026-10-09.md`
§3.2. **UNVERIFIED by me; do not use any row.**

### 5.2 A published exponent whose printed sign is physically impossible
Nugent & Jain 1984's abstract states the dextran fit as *"D = a(Mr)b, in both
normal tissue (a = 10(6) and b = -2.96) and tumor (a = 2.51 X 10(-2) and
b = 1.14)"*. **A positive exponent means D increases with molecular weight**,
which is impossible. Checking the normal-tissue parameters at Mr = 19,400 gives
~2.0 × 10⁻⁷ cm²/s ≈ 20 µm²/s — physically sane. The tumour parameters as printed
give a nonsensical magnitude; with the sign flipped to **−1.14** they give
~32.5 µm²/s, which is sane **and** reproduces the paper's own stated conclusion
that tumour is hindered *less* than normal tissue. **My inference, UNVERIFIED
against the paper body** (Cancer Res 1984, paywalled, no PMCID). Recorded because
anyone lifting `b = 1.14` into this repo would ship a model where diffusivity
rises with molecular weight. Third defect of this class after the ResearchGate
µm figures and the Geiger 0.12 cm²/s unit error.

### 5.3 A misattribution caught in transit
§2.1's paper was reported to me under the wrong authors and wrong pages. Caught
by verifying against PubMed before writing. Recorded in §2.1 rather than silently
fixed, because the *mechanism* — a citation gaining authority at each hop — is
this project's documented recurring defect.

### 5.4 A variable name that conceals a 1.74× convention
`MolecularProfile.free_diffusivity` is **not** a free-solution diffusivity
(§0.2, §2.2c). The disagreement with Degerstedt's measurement **changes sign**
depending on which viscosity convention is used. Flagged as a live hazard, not a
finding.

### 5.5 Minor citation corrections to material reported to me
- *Braz J Med Biol Res* Song 2017 is **50(7):e5403**, not 50(8) (verified).
- MassIVE **MSV000083713** exists with the reported title *"Antiretroviral GI
  Tract Distribution by IR-MALDESI"* (verified by me via MassIVE
  `QueryDatasets`), but the record returns **`"protected":"1"`**, which does not
  obviously match the "public" reading reported to me. **Its licence, access
  terms and pixel size are all UNVERIFIED.** Do not treat it as open without
  checking.

---

## 6. Verified EMPTY — record so nobody re-hunts

**Repositories** (exact queries and counts reported to me by the strand that ran
them; I independently verified only the four Zenodo records in §2.6):
MetaboLights (MSI studies exist, all endogenous; `"drug distribution" AND
imaging` = 0) · PRIDE (`imzML` 10, all endogenous; `MALDI imaging drug
distribution` = 0) · **MassIVE/GNPS** (complete 20,359-dataset index downloaded
and regex-searched: `drug distribution|drug penetration|intratumoral
distribution` = **0**, `autoradiograph` = **0**; exactly one drug map) ·
**IDR** (all 252 containers enumerated — **no tissue drug maps, no MSI**) ·
BioImage Archive (`"drug penetration"` = 0; the two MSI collections are glycans
and lipids, and the ccRCC mIF panel has **no CD31**) · **Zenodo `imzML`** (44
records, all titles enumerated, **zero drug maps**) · **Zenodo elemental
imaging** (`"LA-ICP-MS" AND imaging` = 29, all ice cores and software — **no
tissue platinum/gadolinium/ruthenium imaging**, and this was the brief's "ideal
object") · Dryad · figshare · OSF · **EMPIAR** (0 hits) · **TCIA** (complete
modality list enumerated — **no slide-microscopy and no autoradiography modality
exists**; only non-radiology holdings are H&E co-registration collections) ·
SSBD · **no dedicated MSI archive exists besides METASPACE** (the 2014 proposal
PMC4357645 is built *on* PRIDE). **Cell Image Library is the one item that could
not be closed out — UNVERIFIED.**
*Dead endpoints, recorded to save time:* `empiar/api/search/*` (404),
`api.osf.io/v2/search/` (gone), MassIVE `ProteoSAFe/result_json.jsp` (500).
**I hit two more myself: EuropePMC `/PMCID/fullTextXML` returned HTTP 500
repeatedly today** (the same failure `data_scouting_2026-10-09.md` §2.5 hit on
PMC2629387 — it is the endpoint, not the article), and **PMC PDF downloads are
bot-gated** ("Preparing to download…", a challenge not a paywall, exactly as
`diffusivity_falsification.md` diagnosed for Pruijn 2008). Use
`pmc.ncbi.nlm.nih.gov/articles/PMCxxxx/` HTML with a browser user-agent instead.

**Literature lanes** (reported to me, with searches run):
- **Imaging mass cytometry of platinum vs distance-from-vessel: EMPTY.** Chang
  2016 (PMID 27812005, pancreatic PDX, cisplatin, 1 µm spot) is the only
  IMC-platinum biodistribution paper and has **no endothelial marker at all**.
- **LA-ICP-MS platinum vs distance-from-vessel: EMPTY.** Every platinum elemental
  imaging study in tumours is compartment-based (e.g. Schaier 2025, PMID
  40575293, CC BY: 13.4 µM parenchyma vs **32.6 µM** CAF/vascular vs 87.8 µM
  necrotic — a 2.43× ratio and no distance axis) or depth-from-peritoneal-surface.
  **Report this lane closed**: it is the right community with the right
  instruments and it has **never** published a vessel-distance axis.
- **NanoSIMS / SIMS drug-vs-vessel-distance: EMPTY** (no hits at all).
- **Deposited distance-binned profiles on axes 1–2: EMPTY.** Nobody in this field
  deposits them.
- **Pancreas vessel-distance: EMPTY.** Fanchon 2020 is region-based (hypoxic vs
  oxic) and uses Meca-32, not CD31.
- Head and neck, and sarcoma: nothing with a vessel-distance axis.
- **Nicholson / Hrabětová / Tao raw diffusion curves: not deposited anywhere.**
  That literature publishes **fitted λ and α only**; the realistic route to raw
  curves is an author request. **The repo's known-open item — Vargová 2003's
  paywalled grade III/IV values — is NOT closed by this pass.**
- **Skin and vitreous: wrong geometry/medium**, and explicitly not recommended.
  Skin depth profiling reports `D/H²`, not D, and has no vessel source; vitreous
  is a ~99 %-water gel with no matrix fraction to test.

**One live lead not closed:** Strittmatter N, Sutton D, et al., *Anal Chem*
2022;**94(3):1795–1803**, PMID **35005896**, DOI 10.1021/acs.analchem.1c04579 —
DESI MSI of gemcitabine and its phosphorylated metabolites in KPC pancreatic
tumours **co-registered with immunofluorescence on the same section**. The
best-positioned pancreas candidate in the literature; **whether it carries a
vessel-distance axis is UNVERIFIED** (ACS returned 403, no PMC copy). Worth one
retry with institutional access. Likewise **NCT02567253** lists cisplatin tissue
penetration distance by LA-ICP-MS as a *primary outcome*; **no results
publication was found — UNVERIFIED whether one exists.**

---

## 7. Honest read — could any of this actually refute a Vivum claim?

**Yes. Three of them, and two do it at the level of a model's *form* rather than
a constant's value — which is the stronger kind.**

1. **§2.1 Shafieyan 2014 refutes the linear form of the spatial ECM coupling.**
   Three solutes, same direction, implied coefficients 1.50 / 1.64 / 2.16, all
   **> 1 and therefore inadmissible** in `D_eff = D·(1 − a·ecm)`. The measured
   sensitivity of D to matrix content is 1.27–1.58× steeper than the shipped law
   permits over the paper's own interval. This cannot be absorbed by recalibration.
   Its weakness is the declared `ecm ≡ 1 − φ_fluid` mapping, which §0.1 shows is a
   *declaration* about a unitless field — so **the pre-registration is the whole
   experiment**, and this project withdrew a headline claim on 2026-08-21 over
   exactly a post-hoc substitution.
2. **§2.2 Degerstedt 2024 puts an external, measured number against a shipped
   parameter and the engine is 935× off.** `P_mem = 9.63e-07 µm/s` vs a measured
   `9.00e-04 µm/s`. It is attributable to a named leg, it is CC BY, the numbers
   are in body prose, and it points in the **same direction** as the known
   2.45×-too-long λ. Most importantly it reframes α: `calibration.py` says α
   patches a ~300× missing sink attributed to "binding capacity"; this suggests
   the hole is a measurably wrong **permeability** instead. **That is a worse
   position for the product than the one it currently documents**, which is why
   it ranks second rather than being filed as a curiosity.
3. **§2.4 Shi 2014 partitions a known failure for the first time.** Four of its
   six measured `D_eff/D_free` ratios (0.19 → 0.11) lie **outside** the shipped
   reachable band [0.40, 1.00] altogether, and its small-solute 0.46 — corroborated
   independently by the repo's own Vargová λ (0.416) and by Fu 2005 (0.27 in a
   different tissue and species) — implies the engine is ~2.1× too fast in
   diffusivity on its flagship section. Set against the repo's `D_eff/k` ~6×, that
   splits the error ~2.1× in `D_eff` and ~2.9× in `k`. **Attribution is the thing
   `verification_routes` says the current anchor set cannot do at all.**

**And the honest counterweights, which matter as much:**

- **§2.2b is a CONFIRMATION** — the Ogston–Renkin gate predicts 0.773 against a
  measured 0.745, within 3.8 %. Recorded as a result, not buried. The narrow
  reading is that the gate is right *by construction-independence*: it is fixed at
  φ_f = 0.05 for every tissue (§0.2), so being right for a dilute collagen gel
  says nothing about tumour.
- **§2.2d is a non-test, and I nearly reported it as a refutation.** The
  cell-density arm sits at φ ≈ 0.3 % against the engine's φ_cell = 0.7892 — ~200×
  apart, in a regime where the engine *also* predicts no effect. Checking the
  actual seeding density is what turned a headline into a null.
- **§2.11's dissent is real.** Chary & Jain measure albumin as nearly unhindered,
  the opposite direction from Shi and Fu. With Nugent & Jain's "tumour hindered
  *less* than normal", the literature does not even agree on the **sign** of the
  tumour-vs-normal effect. **The 0.6 constant cannot be validated — or
  refuted — by picking one study**, and a pass that quoted only Shi would have
  manufactured a verdict.
- **Nothing here closes the real gap.** No candidate in this file is a
  spatially-resolved *free-interstitial drug concentration field* in tumour
  tissue. §3 confirms that conclusion from a third independent direction.
- **The input/output imbalance that `data_scouting_2026-10-09.md` §5 named is
  genuinely corrected here, and that is this pass's main structural
  contribution.** Every candidate above scores an **output** or a **constitutive
  law**: the ECM coupling's form, `P_mem`, `k`, `v` and its sign, `D_eff`, the
  source term's symmetry. **Not one is a tissue section.** Three of the twelve are
  also cheap: §2.1 and §2.4 need no download at all, only a pre-registered
  mapping.

**The wrong direction, restated because the commission asked.** Nothing above is
a patient-outcome or treatment-response dataset, and if one is ever recommended
for this product the answer is no, for the reason `data_scouting_2026-10-09.md`
§6 already gives and `DATA_HUNT.md` establishes from 164 datasets: **no open
dataset pairs a measured spatial drug distribution with a per-specimen response
in the same tumour.** A model fitted across that gap has the transport physics as
decoration and the fitted layer absorbs solver error, which makes the simulator
unfalsifiable — the exact opposite of what every item in §2 is for. Nothing here
supports a cross-compound ordering (withdrawn 2026-08-21) or an absolute kill
fraction; §2.3's own finding that **washout rate, not penetration depth,
correlated with activity** is an argument against re-acquiring either.

---

## 8. What I verified myself, and what I did not

**This pass ran four strands: I took axis 4 and all code reading myself, and
three parallel strands covered axes 1–2, axis 3, and axis 5.** Per the standard
`data_scouting_2026-10-09.md` §7 set, the split is stated rather than blurred.

**Verified by me against a primary record or against the code:**
everything in §0.1, §0.2 and §4; the engine values in §0.2 and §2.2 by read-only
lookup; the Ogston–Renkin and reachability arithmetic in §2.1, §2.2 and §2.4 by
running the engine's own pure functions; PMIDs **25028890** (and its correct
authors/pages/PMCID/licence), **37824040** (incl. CC BY verbatim and all gel
numbers from full text), **25143448**, **19807929** (incl. Table 1 values read
from PMC), **24193698** (full abstract), **16507830**, **17473214**, **15971705**,
**2748592**, **28678913**, **6197161**, **30217895**, **28826498**, **41900218**;
and Zenodo **5776125 / 5158808 / 13375209 / 4944035** by the Zenodo API (title,
DOI, licence, date, file manifest). MassIVE **MSV000083713**'s existence and
title by `QueryDatasets`.

**Reported to me and NOT independently verified** — do not quote onward via this
file: the per-pixel internal structure of the Zenodo `.xlsb` workbooks; the
absence of vessel markers in those three papers' full texts; all repository query
counts in §6; the figure-level detail of §2.8 and §2.9; §5.1's reference-list
analysis; §2.1's and §2.4's body-text detail beyond what I quote; the Song 2017
unit defect; and every literature lane declared empty in §6. §5.2's sign
inference is **mine and unverified against the paper body**.

**Not reached:** Strittmatter 2022's vessel axis; NCT02567253 results; Cell Image
Library; Vargová 2003's grade III/IV values (still open); and no attempt was made
on `verification_routes` R7's retrieval, which remains the half-day that decides
that route.

---

## 9. Closing calibration guard

```
$ .venv/bin/python -c "from tmesim.validation.calibration import calibrate_alpha; a=calibrate_alpha(); assert a==22.223477711172563, a; print('alpha OK',a)"
alpha OK 22.223477711172563
```

Unchanged from §0.

No code, data loader, scorecard entry or existing dated record was touched by
this pass. The only file this pass created or modified is this one.

**Nothing in this file is registered, pre-registered, or scored.** Every item is
a recommendation. The two items worth acting on first (§2.1 and §2.4) are
specifically the ones whose entire value depends on **declaring the mapping and
the band before the data is in hand** — because both compare a measured absolute
to a field that §0.1 proves has no absolute scale, and that is the one move this
project has already been burned by.

**Territory note.** A sweep for files modified today returned
`results/data_scouting_2026-10-09.md`, `results/surrogate_feasibility_2026-10-09.md`,
`experiments/surrogate_feasibility/*`, `docs/AUTH_CONTRACT.md` and
`docs/DEPLOYMENT.md`. This pass did not open or write any of them. Flagged rather
than investigated, since another owner may hold them and this project's fan-out
discipline is exclusive file ownership.
