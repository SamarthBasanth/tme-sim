# The held-out patient, run and scored — what it actually tests, and what it does not

*Written after running `tmesim/validation/heldout_patient.py` end to end at three grid
spacings. Every number below was observed in this run. Nothing was tuned, no parameter
was changed, and no result file other than this one was written.*

> ### Headline
>
> **The anchor is real, it runs clean, and it reproduces its published report bit-for-bit
> — but it is not the thing it was described as.** It was left in `OUT_OF_SCOPE`
> "on runtime grounds only". Runtime was never the real objection and the runtime
> figure was attached to the wrong half: the part that carries every external
> measurement costs **126 s and zero PDE solves**. It is now registered.
>
> **The 10 PDE solves are not where the anchor is anchored.** All three
> external, published comparisons in this module (PMIDs 38723085, 36611806,
> 27583442) are *reconstruction and vessel-geometry* quantities. They are
> computed before the transport solver is called and are numerically unchanged
> by it. The five-drug panel — the only part that invokes the transport PDE — is
> scored against **no measured number at all**; it compares the two sections with
> each other. So this is **not** the project's first out-of-sample PDE-invoking
> anchor. It is a genuine out-of-sample *tissue* anchor whose scored rows never
> touch the solver, exactly like the 11 quick-tier rows it was supposed to outrank.
>
> **The verdict on the genuinely held-out rows is FAIL.** Vascular area fraction
> lands at **10.98 %** against a measured **2.1 ± 0.3 %** (5.23×, outside the
> 1.8–2.4 % band). Microvessel density lands at **233 /mm²** against **65 (35–95)**
> (3.59×). Neither is a near miss and neither is grid-dependent.
>
> **One scored claim is not blind, and it is the only one that "replicates".** The
> module reports the endothelial fraction (7.07 %, 0.88× the volume-EM 8.05 %) as
> REPLICATES. The production typing rule that produces that number was selected in
> `results/classifier_fix.md` on a cross-section consistency criterion computed
> **with this section in it**, and that selection is precisely what moved this
> section from **4.18 % (0.52×, outside the module's own declared 1.5× tolerance —
> a FAIL)** to 7.07 % (a PASS). It is registered NOT TESTABLE, not as a pass.
>
> **A separate, load-bearing arithmetic error surfaced.** `results/state_validation.md`
> § 1.2 states that a distance-to-vessel baseline "beats" the solved oxygen field.
> The measured numbers are ρ = **−0.1523** (PDE) versus **−0.1204** (baseline), and
> that same document's table header says "a correct predictor gives a NEGATIVE ρ".
> **|−0.1523| > |−0.1204|: the solved field out-predicts the baseline.** The
> document contradicts itself internally (its own § 3 summary line says the solved
> field wins) and `heldout_patient.py` — which decides this arithmetically rather
> than narratively — returns **REFUTED**.

---

## 0. Correction to the brief

The file is **`tmesim/validation/heldout_patient.py`**, not `experiments/heldout_patient.py`
(that path does not exist). It is a library module with a `__main__`, run as
`./.venv/bin/python -m tmesim.validation.heldout_patient`.

It had also **already been run at least once** — `results/heldout_patient.md` and
`results/heldout_patient.png` are committed. What had never happened is registration:
it was invisible to `vivum validate`, so nothing re-ran it and nothing would have
noticed it breaking.

## 1. What the anchor tests, against what measurement

Both sections are driven through the identical frozen production path
(`reconstruct.real_pipeline.reconstruct_real` → `physics.solve_oxygen` →
`sim.engine.simulate_response`) at identical settings, so section 1's numbers are
re-derived rather than quoted.

| | Section 1 (DEVELOPMENT) | Section 2 (HELD OUT) |
|---|---|---|
| cache | `data_cache/xenium_gbm` | `data_cache/xenium_gbm_io` |
| panel | hBrain_254g + 65-gene add-on, 319 targets | Human Immuno-Oncology hImmune_v1 + Add-on B, 480 targets |
| cells | 40,887 | 816,769 |
| marker map | `brain_319` | `io_480` |

The two panels share **60 genes**. Every composition difference is jointly a
difference in biology and in what is measurable, and neither section carries a
drug-response outcome label.

### 1.1 The three external, published comparisons — all verified

| what | measured value | source | verified |
|---|---|---|---|
| endothelial % of cells | **8.05 %** (4,604 EC of 57,180 cells, 1.05 mm³ volume EM, n=1) | Shapson-Coe A *et al.*, *Science* 2024;384(6696):eadk4858, [PMID 38723085](https://pubmed.ncbi.nlm.nih.gov/38723085/), [DOI](https://doi.org/10.1126/science.adk4858) | PMID resolves to the stated paper |
| CD34+ vascular area fraction | **2.1 ± 0.3 %**, IDH-wt GBM, n=30 | *Cells* 2022;12(1):11, [PMID 36611806](https://pubmed.ncbi.nlm.nih.gov/36611806/), [DOI](https://doi.org/10.3390/cells12010011) | PMID resolves; n=30 IDH1-wt GBM and CD34 microvessels confirmed in the abstract |
| CD31 whole-slide MVD | **65 /mm²** (band 35–95), human GBM, n=20 | Marien KM *et al.*, *PLoS One* 2016;11(9):e0161496, [PMID 27583442](https://pubmed.ncbi.nlm.nih.gov/27583442/), [DOI](https://doi.org/10.1371/journal.pone.0161496) | PMID resolves; whole-slide CD31 MVD method paper incl. glioblastoma |
| doxorubicin perivascular half-distance (the α anchor, upstream) | **40–50 µm** | Primeau AJ *et al.*, *Clin Cancer Res* 2005;11(24):8782-8, [PMID 16361566](https://pubmed.ncbi.nlm.nih.gov/16361566/), [DOI](https://doi.org/10.1158/1078-0432.CCR-05-1664) | verified verbatim from the abstract |

*Article metadata retrieved from PubMed.*

**One citation defect, not a fabrication.** `heldout_patient.py`'s docstring attributes
PMID 36611806 to "Angelucci C et al."; `endothelial_density.py` attributes it to
"Annese T et al.". PubMed gives the author list as **Tamma R, Ingravallo G, Annese T,
d'Amati A, Lorusso L, Ribatti D** — so the first-author attribution is wrong in both
places. The PMID, journal, volume, article number and the measured quantity are all
correct. I did not edit either file; flagging it here.

### 1.2 The pass criteria, as declared in code before the run

* endothelial fraction: fold error within **1.5×** of 8.0 % → band **5.33 %–12.0 %**
  (`ENDO_VOLUME_EM["tolerance_fold"]`).
* vascular area fraction: inside **1.8–2.4 %** (`vessel_network.PUBLISHED`).
* MVD: inside **35–95 /mm²** (`vessel_network.PUBLISHED`).
* drug ordering: section 2's kill ordering equals section 1's — an *internal*
  consistency bar with no measured number behind it.

## 2. Is it genuinely held out? Partly — and the part that "passes" is the part that is not

This is the finding that outranks the verdict, so it comes first.

### 2.1 What IS genuinely held out

None of these was refitted, inspected against, or tuned on section 2:

* the **20 µm vessel linking distance** and the vessel radius model (section-1 and
  literature values);
* the **30 µm tissue-mask radius**;
* **α = 148.6**, fit once in closed form on doxorubicin (§ 4);
* every drug's IC50, MW and BBB penetrance (external GDSC/literature);
* `dx = 12 µm`, chosen for convergence on section 1's geometry.

So the two **vascular-morphometry** rows are genuine out-of-sample tests. That is
what I registered.

### 2.2 Breach 1 — the typing rule was selected with this section in the criterion

`results/classifier_fix.md` § 0 changed the production cell-typing rule from `mean_z`
to `matched_null_z` and reports the effect on all three cached sections:

| endothelial % of cells | s1 | **s2 (held out)** | s3 | spread |
|---|---:|---:|---:|---:|
| before (`mean_z`) | 9.74 | **4.18** | 12.19 | 2.91× |
| after (`matched_null_z`) | 10.50 | **7.07** | 8.75 | 1.49× |

Scored against this module's own declared 1.5× tolerance (band 5.33 %–12.0 %):

| | s1 | **s2** | s3 |
|---|---|---|---|
| before | PASS (1.22×) | **FAIL (0.52×)** | FAIL (1.52×) |
| after | PASS (1.31×) | **PASS (0.88×)** | PASS (1.09×) |

**The rule change is exactly what converts the held-out section's verdict from FAIL
to PASS**, and it does the same on s3. `classifier_fix.md` defends this by saying "the
target was consistency, not the 8.0 % volume-EM anchor, which was never optimised
against". That is true as stated and it does not neutralise the problem: minimising
cross-section *spread* around a section (s1) that already sits near the anchor
necessarily pulls the other sections toward the anchor. The post-fix values cluster
at 7.07–10.50 %, i.e. around ≈ 8.8 %, next to the 8.05 % reference. Agreement produced
that way is not out-of-sample evidence.

### 2.3 Breach 2 — the IO marker map was written against this section

`reconstruct/marker_panels.py`'s docstring records three candidate STROMAL marker sets
**measured on this section**, with the selected one marked `(chosen)`:

```
STROMAL markers                            TUMOR  STROM  IMMUNE  ENDO
SLC1A3,S100B,PTN            (chosen)       79.09   9.98    6.74  4.18
RGS5,PDGFRB,NOTCH3,SLC1A3,S100B,PTN        86.57   2.68    6.85  3.90
RGS5,PDGFRB,NOTCH3                         86.18   4.02    6.76  3.04
```

`results/heldout_patient.md` § 1.2 already discloses this. The choice was made on
internal comparability, not against any external anchor — but the map has seen the
section, so the composition rows are not blind.

### 2.4 Breach 3 — the hypoxia half is not out-of-sample at all, by construction

The 319-gene panel measures no hypoxia-response gene, so the VEGFA/HIF1A test can only
run on section 2 — **the section `results/state_validation.md`'s finding was made on**.
It is a re-derivation through a second code path, not a held-out test. The module says
so itself.

### 2.5 The structural point: the PDE and the anchors do not overlap

Tracing every one of the module's ten scored claims:

| claim | invokes transport PDE? | scored against a published measurement? |
|---|:--:|:--:|
| endothelial fraction vs 8.05 % | no | **yes** (PMID 38723085) |
| endothelial areal density | no | no (labelled CONTEXT) |
| vascular area fraction vs 2.1 % | no | **yes** (PMID 36611806) |
| MVD vs 65 /mm² | no | **yes** (PMID 27583442) |
| kill ordering s1 vs s2 | **yes** | no — section-vs-section |
| coverage recovers kill ordering | **yes** | no — internal |
| coverage still discriminating | **yes** | no — internal |
| ordering is not the IC50 ordering | **yes** | no — internal diagnostic |
| simulated O₂ vs measured VEGFA/HIF1A | oxygen solve | measured transcripts, but **within-section** |
| VEGFA rank 1/480 | oxygen solve | same |

**Every published external number is reached without the transport solver, and every
transport-solver number is scored against nothing published.** The brief's premise —
that this is a genuine out-of-sample *PDE-invoking* anchor worth more than the whole
quick tier — does not survive reading the module. It is worth registering for a
different and still real reason: it is the only **out-of-sample tissue** in the
registry.

## 3. The verdict

### 3.1 Reproduction

The run reproduces `results/heldout_patient.md` **exactly** — same 10 verdicts, and
section 1's kill fractions also reproduce `results/real_drug_panel.md`'s published
numbers (0.0184 / 0.2371 / 0.4362 / 0.0192 / 0.2595) digit for digit. The harness is
running the shipped pipeline, not a variant.

### 3.2 The externally scored rows — both FAIL, neither marginally

| row | section 1 (dev) | **section 2 (HELD OUT)** | measured band | verdict on the held-out section |
|---|---:|---:|---|:--|
| vascular area fraction | 1.25 % (0.59×) | **10.98 % (5.23×)** | 1.8–2.4 % (2.1 ± 0.3) | **FAIL** |
| microvessel density | 110.8 /mm² (1.70×) | **233.3 /mm² (3.59×)** | 35–95 (65) | **FAIL** |
| endothelial cell fraction | 10.50 % (1.31×) | 7.07 % (0.88×) | 5.33–12.0 (module's constant 8.0; H01's exact value is 8.05) | module says PASS; **NOT TESTABLE** (§ 2.2, § 3.4) |

**The vascular area fraction is the load-bearing failure.** It is the only anchor here
that tests the vascular *field* rather than the point pattern, and it has no threshold
to pick — `vessel_geometry.vascular_area_fraction` paints a sub-grid footprint from
data-derived radii. Section 1 misses the published band **low** (0.59×) and the
held-out section misses it **high** (5.23×). The two errors run in opposite directions,
which is why the existing follow-up
(`results/sweep/vascular_area_fraction_generalization_2026-07-30T03-18-42Z.md`) found
a genuine null: no section-invariant radius model closes both.

**MVD is not a held-out failure.** Section 1 misses the same band in the same direction
(110.8 /mm²). Both sections agree with each other and disagree with the literature, so
this is a systematic segmentation property — endothelial cells over-split into vessel
objects — that the second patient did not reveal. It is registered as FAIL because it
is a real miss against a measured value, with that reading in its note.

### 3.3 The transport-PDE claims — no external bar, so no pass or fail

Kill fractions, `uniform` sensitivity mode, α = 148.6, dt = 10 s over 6 h, convection on:

| drug | s1 @ dx12 | s1 @ dx24 | **s2 @ dx12** | **s2 @ dx24** |
|---|---:|---:|---:|---:|
| cisplatin | 0.4362 | 0.4236 | **0.4248** | **0.4240** |
| doxorubicin | 0.2595 | 0.2445 | **0.3603** | **0.3560** |
| lomustine | 0.2371 | 0.2281 | **0.3486** | **0.3468** |
| temozolomide | 0.0184 | 0.0184 | **0.0797** | **0.0793** |
| methotrexate | 0.0192 | 0.0198 | **0.0657** | **0.0656** |

The kill ordering changes on the held-out patient — `methotrexate > temozolomide`
on section 1 becomes `temozolomide > methotrexate` on section 2. Stated honestly:
**that verdict is decided by two near-zero kills separated by 4 %** on section 1
(0.0192 vs 0.0184). It is a real reordering, it is grid-stable (§ 5), and it is also
the weakest possible way for an ordering claim to break.

Two things worth naming about this half:

* **`therapeutic_coverage @1× IC50` has saturated.** On the held-out section **5 of 5**
  drugs sit at either 0.0000 or above 0.95. `results/delivery_metrics.md` flagged this
  saturation as "latent, not active"; at this exposure scale it is active, and the
  metric is behaving as a binary indicator.
* **The coverage ordering fails on *both* sections**, so it is a disagreement with a
  prior configuration, not a generalisation failure.

### 3.4 The hypoxia row — and an error it exposes in `results/state_validation.md`

Measured on all 816,769 cells, every predictor sign-aligned so that a correct predictor
gives a negative ρ:

| predictor | ρ vs hypoxia score, dx=12 | dx=24 |
|---|---:|---:|
| **simulated O₂ (the PDE)** | **−0.1523** | **−0.1538** |
| dist_endothelial (one cKDTree query, no physics) | −0.1204 | −0.1204 |
| vascular_kde (the undiffused source term) | −0.1111 | −0.0855 |

`results/state_validation.md` § 1.2 reads these as:

> "**`dist_endothelial` (−0.1204) beats the PDE (−0.1523) by −21%.**"
> "**`vascular_kde` (−0.1111) also beats the PDE.** … the diffusion step *destroys*
> predictive information rather than adding it."

**That is backwards.** With every predictor sign-aligned — as that document's own table
header states — the strongest predictor is the one with the largest |ρ|, and the PDE has
it. The "−21 %" is itself the tell: 0.1204/0.1523 = 0.79, i.e. the baseline is **79 % as
strong**, which is weaker, not better. The same document contradicts itself in its § 3
summary ("the solved field beats every trivial baseline (best: dist_endothelial at
−0.1204)"). The nested test agrees with the arithmetic: partial ρ(score, simO₂ | dist)
= **−0.1227** while partial ρ(score, −dist | simO₂) = **−0.0792**.

`heldout_patient.py` gets this right because it decides the verdict arithmetically
(`abs(dist_rho) > abs(pde_rho)`) rather than in prose, and returns **REFUTED**.

Two honest limits on how far to push this. It is a **within-section re-derivation**, so
it is not out-of-sample evidence, and ρ ≈ 0.15 is a weak association however it is
ordered — "the PDE wins" here means "wins a contest between three weak predictors".
Also, the positive tail of the gene ranking contains `FLT1` and `PLVAP`, which are
themselves ENDOTHELIAL markers that *build* the oxygen source, so only the negative
tail is informative. **I did not edit `results/state_validation.md`** — I do not own it.

**A minor numerical defect, recorded.** The permutation kernel
(`cellstate_validation.py:504`) emits `divide by zero`, `overflow` and `invalid value`
RuntimeWarnings on this section, and all three reported p-values come out at
**0.00050 = 1/2001**, the resolution floor of 2,000 draws. Those p-values therefore
carry no information beyond "below the estimator's floor". Effect sizes are computed by
a different kernel on all cells and are unaffected.

## 4. Conditionals inherited

### 4.1 α — and what it does to independence

`calibration.calibrate_alpha()` is a closed-form solve of **one scalar equation in one
scalar unknown**, fit to the doxorubicin perivascular anchor (PMID 16361566, 40–50 µm
half-distance → L\* = 64.92 µm e-folding). This run used **α = 148.6**, not refitted for
section 2.

Per `results/alpha_identifiability.md`:

* α is a **gauge coordinate**. Scaling every drug's binding prior by *t* and refitting
  moves α over four orders of magnitude while every downstream number is identical to
  2.8e-14 µm.
* The anchor identifies **exactly one number**: `k_dox = 3.331353e-02 /s`.
* **Nine per-drug binding ratios `b_d/b_dox` are unconstrained** — they come from
  `molecular/binding.py`'s logP/HBD sigmoid with zero data constraint, and
  `results/binding_sink_refuted.md` shows that sigmoid is refuted in *direction*, not
  merely magnitude.

**The split matters here, and it is unusually clean:**

* **The registered rows are NOT downstream of α.** Composition and vessel morphometry
  are computed from cell positions and transcript counts before any drug parameter is
  compiled. Their independence of the calibration is total.
* **The drug panel IS entirely downstream of α** — every drug goes through
  `calibrated_drugparams(name, alpha)`. And the *cross-drug ordering* is the one thing
  the module calls transportable, so the transportable claim rests **precisely on the
  nine ratios the anchor does not constrain**. A held-out ordering result is not
  evidence about the sink model; it is evidence conditional on a sink model that has
  been independently falsified in direction.

### 4.2 The other conditionals

* **Delivery law.** `PDEConfig.transvascular_exchange` defaults to `False`, so the drug
  panel ran the one-way leak `S = PS·Cp`. The module reaches the solver through
  `sim.engine.simulate_response`, which is the declared **production** wiring point in
  `validation/propagation.py`, so the flag is not silently inert here. **I injected no
  `exchange_k` and added no `WiringPoint`** — nothing in this work needed one.
* **Absolute oxygen is uncalibrated on both sections.** `OxygenConfig.consumption` has
  never been fitted to a measured oxygen tension, so the 33.8 %-vs-0 % hypoxic-fraction
  gap is a statement about the sink constant and cell density, not about the tumours.
  Every absolute kill inherits it.
* **ECM is not comparable across the panels.** The IO panel carries 3 of the 11 matrix
  transcripts the ECM field is built from (COL1A2 and POSTN, both fibrillar-collagen
  markers, are absent), and the field is a max-normalised KDE over whichever markers are
  measured. Mean ECM 0.457 vs 0.076 is panel content, not biology.
* **No outcome label** on either section. Nothing here validates a kill prediction.

## 5. Grid check

### 5.1 The registered rows are grid-invariant — measured, not assumed

Section 2, re-run at three spacings spanning a **9× range in cell area**:

| dx (µm) | grid | endothelial % | vascular area fraction % | MVD /mm² | vessel objects | tumour > 100 µm |
|---:|---|---:|---:|---:|---:|---:|
| 24 | 453 × 479 | 7.0685 | 10.9822 | 233.2702 | 20,875 | 0.003532 |
| 12 | 905 × 956 | 7.0685 | 10.9819 | 233.2739 | 20,875 | 0.003532 |
| 8 | 1356 × 1433 | 7.0685 | **10.9908** | **233.4677** | 20,875 | 0.003532 |

Maximum spread: endothelial fraction **0** (bit-identical), area fraction **0.08 %**,
MVD **0.08 %**, vessel count and tumour-distance fraction **0**. This is structural, not
luck: cell typing is per-cell and never touches the grid; vessel objects come from
single-linkage clustering on cell coordinates; the area fraction uses a sub-grid
analytic footprint. The only grid-dependent term is the tissue-mask area denominator.

**The verdicts cannot move.** The area fraction would need to fall **4.6×** and MVD
**2.5×** to re-enter their bands; the observed grid sensitivity is 0.08 %. Neither row
sits near a band edge, so the coarser/finer requirement is satisfied and answered.

This is the opposite of the concurrent breast-panel finding (docetaxel −2.4 %
discretisation error at dx = 12, every drug monotone in dx) — and consistently so,
because those are transport-PDE quantities and these are not.

### 5.2 The PDE half is grid-sensitive, but every verdict held

All **10 verdicts are identical at dx = 24 and dx = 12**. Kill fractions move ≤ 1.2 % on
the held-out section and ≤ 6.1 % on section 1 (doxorubicin, the largest). The near-ties
keep their direction:

| near-tie (section 2) | dx = 12 | dx = 24 | order preserved? |
|---|---|---|:--:|
| coverage: lomustine vs doxorubicin | 0.9702 vs 0.9691 (**0.11 %**) | 0.9603 vs 0.9519 (0.9 %) | yes |
| kill: doxorubicin vs lomustine | 0.3603 vs 0.3486 (3.4 %) | 0.3560 vs 0.3468 (2.7 %) | yes |
| kill: temozolomide vs methotrexate | 0.0797 vs 0.0657 (21 %) | 0.0793 vs 0.0656 (21 %) | yes |

The 0.11 % coverage margin is genuinely a coin-flip at dx = 12 and should not be cited;
it widens and keeps its sign at dx = 24, so the reported ordering is not an artefact of
that one comparison.

### 5.3 CFL — the finest admissible grid, measured

The explicit upwind advection step raises if `dt·(|vx|+|vy|)/dx > 1`. Measured maximum
speeds from `solve_ifp`:

| section | dx (µm) | grid | max ‖v‖₁ (µm/s) | CFL @ dt = 10 s |
|---|---:|---|---:|---:|
| section 1 | 12 | 424 × 675 | 0.6336 | 0.5280 |
| section 1 | 8 | 635 × 1012 | 0.7513 | 0.9391 |
| section 2 | 12 | 905 × 956 | 0.5869 | 0.4891 |
| section 2 | 8 | 1356 × 1433 | 0.7882 | **0.9853** |

**dx = 8 µm is the finest reachable grid at the frozen dt = 10 s**, with essentially no
margin on the held-out section (0.985). dx = 6 is refused. Going finer would require
cutting dt, which confounds spatial with temporal error, so it was not done.

This also corrects the module's own docstring: `DX_UM`'s comment calls dx = 12 "the
finest grid satisfying the advection CFL at dt = 10 s". It is not — dx = 8 satisfies it.
Note that max ‖v‖₁ **grows** as the mesh refines (0.587 → 0.788 on section 2) because the
IFP gradient is better resolved, so the CFL limit tightens faster than 1/dx.

**No dx = 8 drug-panel arm was run.** It would cost several hours of solver time to
refine claims that carry no external measurement, while the rows that do carry one are
already shown invariant to 0.08 % over dx 24→8.

## 6. What was registered

`heldout-patient-geometry` is now a **full-tier, geometry-scope** anchor group in
`tmesim/validation/scorecard.py`. Measured runtime **125.8 s**, **zero PDE solves**.

```
SUMMARY   PASS 0   FAIL 2   NOT TESTABLE 1   (of 3 rows)
VERDICT   2 anchor(s) FAIL against measured data. Exiting 1.
```

| row | status |
|---|---|
| vascular area fraction (held-out section) — 10.98 % vs 1.8–2.4 % | **FAIL** |
| microvessel density (held-out section) — 233 /mm² vs 35–95 | **FAIL** |
| endothelial cell fraction (held-out section) — 7.07 % vs 8.05 % | **NOT TESTABLE** |

The endothelial row is NOT TESTABLE for two independent reasons, both carried in the
row: the registry's existing `endothelial-cellularity` group already reports this exact
comparison NOT TESTABLE (a 3-D whole-cell volume census against 2-D polygons on a ~5 µm
section, reference is normal cortex, n = 1) — a second section does not repair that —
and § 2.2's non-blindness. Registering it as a PASS would have put two contradictory
verdicts on the same published number in the same scorecard.

**Why only the geometry half.** The drug panel is scored against no measured number and
is downstream of α; the hypoxia check is model-vs-marker within one section, the same
category `state_validation.py` is already excluded under. The `OUT_OF_SCOPE` entry was
rewritten to say this, replacing the "~18 min and 10 PDE solves — too slow" reason,
which was wrong twice: runtime was not the objection, and the 18 min belonged to the
half that carries no anchor.

**No double counting.** This is structurally the same move the registry already makes
with `breast-vascular-morphometry`: the same anchors, a different specimen. What is new
is the specimen, not the measurement.

## 7. Reproduction

```bash
./.venv/bin/python -m tmesim.validation.heldout_patient          # full module, ~18 min nominal
./.venv/bin/python -m tmesim.cli validate --only heldout-patient-geometry   # registered rows, ~126 s
```

Deterministic (`seed = 0` everywhere; permutation seed 0). This run was executed at
dx = 12 and dx = 24 for the full module and at dx = 24 / 12 / 8 for the geometry rows,
with output directed outside `results/` so that `results/heldout_patient.md` was not
overwritten — which is also what made the bit-for-bit reproduction check possible.

**Files changed by this work:** `results/heldout_patient_validation.md` (this file) and
`tmesim/validation/scorecard.py` (registration only). Nothing else was modified; in
particular `results/state_validation.md`, `results/classifier_fix.md`,
`tmesim/validation/heldout_patient.py` and `tmesim/reconstruct/marker_panels.py` are
untouched despite the defects reported above.

## 8. What could make this wrong

1. **n = 2.** Two sections cannot separate "the pipeline does not generalise" from
   "these two tumours differ". A *failure* to replicate is well supported by one
   discrepancy; a *success* rests on one agreement and is much weaker.
2. **Different patients AND different panels** (60 shared genes). The composition rows
   are confounded by panel content; § 2.3 and the ECM note in § 4.2 name the specific
   places.
3. **The published bands are cohort means**, which `vessel_network.segmentation_report`
   itself calls "context, not a target" for a single patient. A 5.23× miss survives that
   caveat comfortably; a 1.2× miss would not have.
4. **Under-resolution under-estimates tissue concentration** (`|rel err| = 1.427e-3 ·
   dx^1.939 · δ^-1.145`, from the boundary-layer study). This is a candidate contributor
   for any PDE quantity that comes out low — but the two registered rows are not PDE
   quantities and are measured invariant to 0.08 %, so it does not apply to them.
5. **Pericytes are unlabelled on the held-out section.** `CellType` has no MURAL member,
   so RGS5/PDGFRB/NOTCH3+ mural cells fall into the ungated TUMOR default, inflating the
   tumour count and hence the oxygen sink.
6. **2-D, in-plane, single section** on both. Distances are Euclidean and in-plane.
