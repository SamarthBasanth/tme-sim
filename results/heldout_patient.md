# The held-out patient — running the whole validated pipeline on a second real GBM section

> **Regenerated under the discrete vessel source.** `vascular_mode = discrete`
> (`reconstruct.real_fields._DEFAULT_VASCULAR_MODE`, promoted in
> [`results/vessel_promotion.md`](vessel_promotion.md)) · **alpha = 148.6**, the
> e-folding-length convention `L = d_half / ln 2` · per-cell sensitivity mode =
> **`uniform` (both sections)** · **dx = 12 µm (both sections)**.
>
> Every absolute exposure / AUC / coverage / kill below is on the **discrete-source
> unit** and is NOT comparable with the KDE-era numbers this file carried before this
> regeneration (the tissue-integrated source total is 3.36x smaller; that is a change
> of *unit*, not of physics — the dose-matched, geometry-only effect on coverage is
> <= 0.0166). To recover the old numbers, re-run with `mode='kde'`.
>
> **This file was NOT carrying a stale banner but was stale anyway** — it was produced at 23:59 on the KDE default, ~30 min before `real_fields._DEFAULT_VASCULAR_MODE` was flipped to `discrete`, and the post-promotion banner sweep missed it. It has now been re-run end to end on both sections under the discrete source. **The headline generalisation result is unchanged: the kill ordering (cisplatin > doxorubicin > lomustine > methotrexate > temozolomide) REPLICATES on the held-out 816,769-cell section**, and `therapeutic_coverage` @1x IC50 still fails to recover it on either section.

Every validated result in this project was developed on **one** tissue: the 319-gene 10x Xenium GBM section (40,887 cells). This report re-runs the entire pipeline — reconstruction, vessel segmentation, the five-drug simulation and the hypoxia check — on a **second real patient the pipeline was never developed on**: the 480-gene Immuno-Oncology Xenium GBM section (816,769 cells). Both sections are driven through the identical frozen production path at the identical settings, so **section 1's numbers here are re-derived, not quoted**.

**5 of 7 claims that could be scored as generalisation replicate.** A further 2 fail on **both** sections — those are not generalisation failures, they are claims that do not hold in this configuration at all, and they are called out separately because conflating the two would be the easiest way to make this report say something it should not. The headline, one line each:

| what | section 1 (development) | section 2 (held out) | verdict |
|---|---|---|:--|
| ENDOTHELIAL fraction matches human volume-EM 8.0% (PMID 38723085) | 9.74 % | 4.18 % | **DOES NOT REPLICATE** |
| ENDOTHELIAL areal density (the oxygen/drug SOURCE) is comparable | 142.79 cells/mm2 | 381.94 cells/mm2 | **CONTEXT** |
| Vascular area fraction matches measured 2.1 +/- 0.3% (PMID 36611806) | 2.50 % | 5.79 % | **DOES NOT REPLICATE** |
| MVD lands in the published human GBM range 35-95/mm2 (PMID 27583442) | 110.82 vessels/mm2 | 104.17 vessels/mm2 | **NOT REPRODUCED ON EITHER** |
| Kill ordering is cisplatin > doxorubicin > lomustine > methotrexate > temozolomide | cisplatin > doxorubicin > lomustine > methotrexate > temozolomide  | cisplatin > doxorubicin > lomustine > methotrexate > temozolomide  | **REPLICATES** |
| therapeutic_coverage @1x IC50 recovers the kill ordering | cisplatin > doxorubicin > lomustine > temozolomide > methotrexate  | cisplatin > lomustine > doxorubicin > temozolomide > methotrexate  | **NOT REPRODUCED ON EITHER** |
| therapeutic_coverage @1x IC50 is still a DISCRIMINATING headline metric | 3/5 drugs pinned at 0 or >0.9  | 3/5 drugs pinned at 0 or >0.9  | **HOLDS** |
| The ordering is NOT simply the IC50 ordering read back out | doxorubicin > methotrexate > cisplatin > lomustine > temozolomide  | doxorubicin > methotrexate > cisplatin > lomustine > temozolomide  | **HOLDS** |
| Simulated oxygen is BEATEN by a distance-to-nearest-vessel baseline (state_validation.md Test 1) | UNTESTABLE (no hypoxia gene on the 319 panel) Spearman rho | rho_PDE -0.1135 vs rho_dist -0.1503 Spearman rho | **CONFIRMED** |
| VEGFA is the most anti-correlated gene on the panel; HIF1A is null | n/a — neither gene is on the 319 panel  | VEGFA rank 1/480, HIF1A rank 233/480  | **CONFIRMED** |

> **THE CONFOUND, FIRST.** These are **different patients** *and* **different panels** — the two Xenium panels share only **60 genes**. The IO panel carries none of PECAM1 / SOX2 / OLIG2 / TOP2A / AQP4 / GJA1 / MOG, so the endothelial gate runs on VWF/FLT1/PLVAP instead of PECAM1 and the stromal gate on SLC1A3/S100B/PTN instead of AQP4/GJA1/MOG. Every composition difference below is therefore **jointly** a difference in biology and a difference in what is measurable. Where the two can be separated, this report separates them; where they cannot, it says so. Neither section carries a drug-response outcome label, so nothing here validates a kill prediction — only the machinery that produces one.

## 0. What was held out, and what was not

This is an out-of-sample test of the **pipeline**, not of a fitted model. Being explicit about what that does and does not cover:

| held out from section 2 | NOT held out |
|---|---|
| the reconstruction was never inspected on it while the marker gate, the typing threshold or the 30 µm tissue-mask radius were chosen | the IO-panel **marker map** in `reconstruct/marker_panels.py` was written for this section (in `results/state_validation.md`), by choosing among measured candidate marker sets — see §1.2 |
| the 20 µm vessel linking distance and the 7.5 µm vessel radius were fixed against section 1 and are used here unchanged | the oxygen solver's `consumption` constant has never been fitted to a measured oxygen tension on **either** section |
| the perivascular-sink calibration **α = 148.6** was fit once on doxorubicin against section 1 and is **not** refitted here | the drug IC50s are external (GDSC/literature) and are the same numbers on both |
| every threshold in `delivery_metrics` (the p95 vessel threshold, the 1× IC50 coverage bar) | `dx = 12 µm` was chosen for grid convergence on section 1; convergence is **not** re-established here (§6) |

## 1. Reconstruction — does the marker-gated classifier survive a different panel?

| | Section 1 (development) | Section 2 (HELD OUT) |
|---|---|---|
| dataset | `data_cache/xenium_gbm` | `data_cache/xenium_gbm_io` |
| panel | hBrain_254g + 65-gene add-on (319 targets), Xenium v1.3.0 | Human Immuno-Oncology Profiling hImmune_v1 + Add-on B (480 targets), Xenium v2.0.0 |
| marker map selected | `brain_319` | `io_480` |
| marker coverage (brain / IO map) | 1.00 / 0.50 | 0.50 / 1.00 |
| cells | 40,887 | 816,769 |
| grid @ dx = 12 µm | (424, 675) | (905, 956) |
| tissue extent (µm) | 8100 × 5088 | 11472 × 10860 |
| tissue area (geometric mask) | 27.9 mm² | 89.5 mm² |
| **cell density** | 1465 /mm² | 9127 /mm² |
| TUMOR | 73.80% | 79.09% |
| STROMAL | 7.92% | 9.98% |
| IMMUNE | 8.54% | 6.74% |
| ENDOTHELIAL | 9.74% | 4.18% |
| mean vascular density | 0.0250 | 0.0579 |
| mean ECM density | 0.0759 | 0.4571 |
| ECM markers on panel | COL1A2, COL12A1, COL25A1, POSTN, VCAN, DCN, FBLN1, LAMA2, THBS1, SPON1, TGFBI | VCAN, DCN, FBLN1 |
| mean oxygen (tissue) | 0.2955 | 0.1232 |
| tissue coverage of bbox | 67.70% | 71.83% |

### 1.1 The endothelial check — the one external anchor the reconstruction has

Human cortex volume electron microscopy counts **8.0%** of cells as endothelial (Shapson-Coe A et al., Science 2024;384(6696):eadk4858 — 'A petavoxel fragment of human cerebral cortex reconstructed at nanoscale resolution', [PMID 38723085](https://pubmed.ncbi.nlm.nih.gov/38723085/)). Section 1's **9.74%** matching that — rather than the 1–2% that dissociation-based scRNA-seq reports — is the single strongest external check the reconstruction has, and `results/SESSION_AUDIT.md` §4 leans on it.

**On the held-out patient the classifier calls 4.18% endothelial — 0.52× the volume-EM value, against section 1's 1.22×.**

| | section 1 | section 2 | ratio |
|---|---:|---:|---:|
| ENDOTHELIAL, % of cells | 9.74% | 4.18% | 0.43× |
| ENDOTHELIAL, cells per mm² | 143 | 382 | 2.67× |
| all cells per mm² | 1465 | 9127 | 6.23× |

**These two rows point in opposite directions, and that is the finding.** The endothelial *fraction* falls while the endothelial *areal density* rises, because section 2's tissue is far more cellular. The oxygen and drug source term is built from a KDE of endothelial **positions**, so it depends on the density, not the fraction — the compartment that matters physically did not shrink. But the *fraction* is the quantity that was validated against volume EM, and that quantity does not replicate.

### 1.2 Why this is a weaker test than it looks — and still a real one

The IO-panel marker map was not written blind. `results/state_validation.md` §'MURAL / PERICYTE' records that three candidate STROMAL marker sets were measured on **this section** and the astro-lineage set was chosen partly because it kept the two sections comparable. So the map has seen section 2. What it has **not** seen is any external anchor: the choice was made on internal comparability, never against the 8.0% volume-EM number or against any vascular measurement. The endothelial gate in particular (`VWF, FLT1, PLVAP`) was picked because those are the canonical endothelial genes the panel carries, not because they produced a particular fraction.

| class | section 1 markers used | section 2 markers used |
|---|---|---|
| TUMOR | EGFR, SOX2, OLIG2, MKI67 | EGFR, PDGFRA, IL13RA2, MET |
| IMMUNE | PTPRC, CD68, AIF1, CD163 | PTPRC, CD68, AIF1, CD163 |
| ENDOTHELIAL | PECAM1 | VWF, FLT1, PLVAP |
| STROMAL | MOG, AQP4, GJA1 | SLC1A3, S100B, PTN |

Three candidate explanations for the endothelial gap, and what separates them:

1. **Real biology.** Two different tumours genuinely differ in vascularity. Partly supported — but the *areal* endothelial density went **up**, not down, so this cannot be the whole story.
2. **Denominator.** Section 2 has 6.2× the cell density, so the same vasculature is divided by many more non-endothelial cells. The denominator alone is arithmetically sufficient to explain the whole gap: divide section 2's endothelial areal density by section 1's total cell density and the fraction would read 26.1% — far ABOVE 8.0%, not the 4.18% observed. The observed fraction is set as much by how many cells the section contains as by how many of them are endothelial.
3. **Panel content.** A 3-marker mean-z gate (VWF/FLT1/PLVAP) behaves differently from a 1-marker gate (PECAM1). Measured detection rates — section 1: PECAM1 24.57%; section 2: VWF 4.15%, FLT1 9.42%, PLVAP 4.36% — so on section 2 the gate is the mean z of three sparse markers, which a cell strongly positive for one but silent for the other two can fail. That pushes the endothelial call **down** for a purely technical reason.

**Explanations 2 and 3 are measurement artefacts, not biology, and neither can be excluded here.** The honest statement is that *the endothelial cell-type fraction is not a portable quantity across panels*, and the 8.0% agreement on section 1 should be read as a property of that section and that gate, not as a validated property of the classifier.

## 2. Vascular morphometry — the only vascular quantity scored against a measurement

`reconstruct/vessel_network.py`, unmodified: measured endothelial cells single-linkage clustered at **20 µm** into discrete vessel objects, then painted as a source of radius **7.5 µm**. Both constants were fixed against section 1 and are used here unchanged.

| quantity | section 1 | section 2 (HELD OUT) | published | range | s1 verdict | s2 verdict | source |
|---|---:|---:|---:|---|:--|:--|---|
| microvessel density, human GBM (n=20), CD31 whole-slide, unbiased systematic-uniform-random sampling (vessels/mm^2) | 110.82 | **104.17** | 65.00 | 35.0–95.0 | HIGH | **HIGH** | [PMID 27583442](https://pubmed.ncbi.nlm.nih.gov/27583442/) |
| 1/sqrt(MVD) spacing implied by the MVD anchor above (um) | 94.99 | **97.98** | 124.00 | 102.6–169.0 | LOW | **LOW** | [PMID 27583442](https://pubmed.ncbi.nlm.nih.gov/27583442/) |
| median tumour-cell-to-nearest-vessel distance implied by the MVD anchor under a random (Poisson) vessel arrangement, median = sqrt(ln2 / (pi q)) (um) | 43.29 | **39.62** | 58.00 | 48.2–79.4 | LOW | **LOW** | [PMID 27583442](https://pubmed.ncbi.nlm.nih.gov/27583442/) |
| CD34+ vascular area fraction, human IDH-wt GBM (n=30), whole-field positive-pixel count (% of tissue area) | 2.50 | **5.79** | 2.10 | 1.8–2.4 | HIGH | **HIGH** | [PMID 36611806](https://pubmed.ncbi.nlm.nih.gov/36611806/) |

### 2.1 The load-bearing row

`results/vessel_network_impact.md` calls the **vascular area fraction** the load-bearing row, because it is the only one that tests the *field* rather than the point pattern, and because the vessel radius was fixed from published vessel-caliber imaging **before** the field was built — so it is an untuned comparison. Section 1 scored **2.50%** against a measured whole-section CD34 value of **2.1 ± 0.3%**, i.e. 1.19×.

**Section 2 scores 5.79% — 2.76× the measured value, and 2.32× section 1's.** That is outside the published band by a wide margin, and it is the clearest single failure of generalisation in this report.

The cause is not mysterious, and it is worth stating precisely because it is *not* 'the second patient is more vascular by that factor'. The area fraction is (number of vessel objects) × (area painted per object) ÷ (tissue area). Section 2 has 34,179 endothelial cells over 89.5 mm² versus section 1's 3,984 over 27.9 mm² — 2.67× the areal density — and the fixed 7.5 µm radius converts that density directly into painted area. The radius is a **grid-constrained constant**: at dx = 12 µm one grid cell already has an equivalent disc radius of 6.8 µm, so 7.5 µm is the smallest radius that reliably paints ≥1 cell per vessel. On a denser section that floor is hit by more vessels, and the area fraction inflates. **The 1.19× agreement on section 1 was therefore partly a coincidence of that section's endothelial density matching the grid-imposed radius.**

### 2.2 Vessel objects and the clustering that makes them comparable to histology

| | section 1 | section 2 |
|---|---:|---:|
| endothelial cells | 3,984 | 34,179 |
| discrete vessel objects | 3,092 | 9,322 |
| tissue area (mm²) | 27.9 | 89.5 |
| MVD (/mm²) | 110.8 | 104.2 |
| endothelial cells per vessel | 1.29 | 3.67 |
| singleton objects | 84.8% | 42.9% |
| median object extent, multi-cell (µm) | 18.4 | 22.7 |
| median aspect ratio, multi-cell | 2.29 | 2.56 |
| tumour cells > 100 µm from a vessel | 5.32% | 6.60% |
| tumour cells > 150 µm from a vessel | 0.49% | 1.14% |
| intervessel NN median (µm) | 40.1 | 31.0 |
| tumour→vessel median (µm) | 43.3 | 39.6 |
| tumour→vessel p90 (µm) | 86.8 | 88.8 |
| clustering index vs Poisson (median) | 0.97× | 0.86× |
| clustering index vs Poisson (>100 µm tail) | 1.73× | 1.74× |

**One thing genuinely improves out of sample.** `results/vascular_validation.md` §3.2 found that on section 1 the clustering step 'barely bites' — 1.29 endothelial cells per putative vessel and 85% singletons — so the reconstruction was resolving scattered endothelial nuclei rather than vessel walls. On section 2 it merges 3.67 cells per object with 43% singletons and a median multi-cell extent of 23 µm at aspect ratio 2.6: **elongated, multi-cell objects that actually look like vessel profiles.** The denser section supports the segmentation the sparse one could not. That is a genuine, if backhanded, out-of-sample success — and it is also part of why the area fraction inflates.

### 2.3 Linking-distance sweep — the free parameter, swept on both sections

| link (µm) | s1 vessels | s1 MVD | s1 cells/vessel | s2 vessels | s2 MVD | s2 cells/vessel |
|---:|---:|---:|---:|---:|---:|---:|
| 10 | 3,910 | 140.1 | 1.02 | 24,034 | 268.6 | 1.42 |
| 15 | 3,577 | 128.2 | 1.11 | 13,265 | 148.2 | 2.58 |
| 20 | 3,092 | 110.8 | 1.29 | 9,322 | 104.2 | 3.67 |
| 25 | 2,707 | 97.0 | 1.47 | 7,604 | 85.0 | 4.49 |
| 30 | 2,456 | 88.0 | 1.62 | 6,531 | 73.0 | 5.23 |
| 40 | 2,120 | 76.0 | 1.88 | 4,868 | 54.4 | 7.02 |

On section 2 the published MVD range (35–95/mm²) is reached at a linking distance of ~25–40 µm rather than the 20 µm fixed on section 1. **No linking distance was re-tuned for this report** — reporting the sweep is the honest alternative to tuning it.

## 3. The drug panel — does the ordering hold on a second patient?

Five drugs, **uniform** per-cell sensitivity (the only fair cross-drug mode — every drug carries 1.0 at every cell, so the factor cancels and the ranking isolates delivery × potency), dx = 12 µm, α = 148.6 fit once on section 1 and **not refitted**, PDE dt = 10 s over 6 h with convection on. Identical settings on both sections.

### 3.1 Section 1 (development)

| drug | IC50 (µM) | BBB | D_eff (µm²/s) | mean tumour AUC (µM·s) | **coverage @1× IC50** | coverage @0.5×/2× | median E/IC50 | **kill** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| temozolomide | 606.047 | 0.125 | 265.3 | 7.06e+05 | **0.0000** | 0.000 / 0.000 | 0.05 | **0.0147** |
| lomustine | 50 | 0.973 | 247.4 | 8.7e+05 | **0.2175** | 0.804 / 0.021 | 0.73 | **0.1975** |
| cisplatin | 24.5503 | 0.858 | 225.4 | 1.69e+06 | **0.9936** | 0.998 / 0.853 | 2.97 | **0.3699** |
| methotrexate | 12 | 0.023 | 189.8 | 1.48e+04 | **0.0000** | 0.000 / 0.000 | 0.05 | **0.0162** |
| doxorubicin | 0.5 | 0.021 | 176.3 | 1.08e+04 | **0.3836** | 0.826 / 0.064 | 0.86 | **0.2165** |

### 3.2 Section 2 (HELD OUT)

| drug | IC50 (µM) | BBB | D_eff (µm²/s) | mean tumour AUC (µM·s) | **coverage @1× IC50** | coverage @0.5×/2× | median E/IC50 | **kill** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| temozolomide | 606.047 | 0.125 | 201.7 | 1.64e+06 | **0.0000** | 0.000 / 0.000 | 0.12 | **0.0274** |
| lomustine | 50 | 0.973 | 188.1 | 1.95e+06 | **0.8013** | 0.957 / 0.365 | 1.66 | **0.2210** |
| cisplatin | 24.5503 | 0.858 | 171.4 | 3.85e+06 | **0.9976** | 0.999 / 0.980 | 6.80 | **0.3162** |
| methotrexate | 12 | 0.023 | 144.3 | 3.25e+04 | **0.0000** | 0.001 / 0.000 | 0.11 | **0.0276** |
| doxorubicin | 0.5 | 0.021 | 134.1 | 2.36e+04 | **0.7908** | 0.933 / 0.477 | 1.92 | **0.2284** |

### 3.3 The ordering

| ranked by | section 1 | section 2 (HELD OUT) | same? |
|---|---|---|:--:|
| simulated **kill** | cisplatin > doxorubicin > lomustine > methotrexate > temozolomide | cisplatin > doxorubicin > lomustine > methotrexate > temozolomide | ✅ |
| **coverage @1× IC50** | cisplatin > doxorubicin > lomustine > temozolomide > methotrexate | cisplatin > lomustine > doxorubicin > temozolomide > methotrexate | ❌ |
| mean tumour AUC (delivery alone) | cisplatin > lomustine > temozolomide > methotrexate > doxorubicin | cisplatin > lomustine > temozolomide > methotrexate > doxorubicin | ✅ |
| 1/IC50 (potency alone) | doxorubicin > methotrexate > cisplatin > lomustine > temozolomide | doxorubicin > methotrexate > cisplatin > lomustine > temozolomide | ✅ |
| AUC/IC50 (delivery × potency) | cisplatin > doxorubicin > lomustine > methotrexate > temozolomide | cisplatin > doxorubicin > lomustine > methotrexate > temozolomide | ✅ |

*(The 'same?' column compares the two SECTIONS with each other. Whether each ranking agrees with the KILL ranking — the separate question — is read down the rows: on both sections `coverage` swaps doxorubicin and lomustine relative to kill, and `AUC/IC50` reproduces kill exactly.)*

**Published ordering** (`results/real_drug_panel.md`, uniform mode): cisplatin > doxorubicin > lomustine > methotrexate > temozolomide. Section 1's kill fractions here (temozolomide 0.0147, lomustine 0.1975, cisplatin 0.3699, methotrexate 0.0162, doxorubicin 0.2165) reproduce that report's numbers exactly, which is the check that this harness is running the same pipeline and not a variant of it.

**The kill ordering is IDENTICAL on both patients.** That is the headline, and the immediate next question is whether it is identical for a *trivial* reason.

### 3.4 Is the ordering trivially the IC50 ordering?

An ordering that just reproduces the potency ranking would be reproducible by a spreadsheet with no geometry, no transport and no patient in it. The diagnostic is whether delivery ever overturns potency.

| | section 1 | section 2 |
|---|---:|---:|
| IC50 fold range across the panel | 1212× | 1212× |
| mean-tumour-AUC fold range | 156× | 163× |
| ρ(kill, 1/IC50) — potency alone | +0.50 | +0.50 |
| ρ(kill, mean AUC) — delivery alone | +0.30 | +0.30 |
| ρ(kill, AUC/IC50) — the product | +1.00 | +1.00 |
| kill spread (max − min) | 0.355 | 0.289 |
| coverage spread (max − min) | 0.994 | 0.998 |

**It is NOT the IC50 ordering.** Potency alone would rank doxorubicin > methotrexate > cisplatin > lomustine > temozolomide; the simulation ranks cisplatin > doxorubicin > lomustine > methotrexate > temozolomide. Delivery spans 163× and potency 1212×, and the product AUC/IC50 tracks kill at ρ = +1.00 — the delivery × potency claim from `results/SESSION_AUDIT.md` §4 survives on the held-out patient.

### 3.5 Why the absolute numbers move so much, and why that is NOT a ranking result

| | section 1 | section 2 |
|---|---:|---:|
| mean oxygen over tissue | 0.2955 | 0.1232 |
| model-HYPOXIC tumour fraction | 12.68% | 96.90% |
| model-VIABLE tumour fraction | 0.00% | 0.00% |
| mean ECM density | 0.076 | 0.457 |

Two mechanisms move every absolute kill on section 2, and **neither is a drug property**, so both largely cancel out of the ranking:

1. **Oxygen.** Mean tissue O₂ is 0.123 on section 2 against 0.295 on section 1 — the denser tissue is a bigger sink — so 97% of tumour cells fall in the model's HYPOXIC band versus 13%, and `PDConfig.hypoxia_resistance` (×0.35) plus `quiescent_resistance` (×0.5) then suppress every drug's kill together. **`results/state_validation.md` §1.6 already established that the absolute oxygen level is uncalibrated** — `OxygenConfig.consumption` has never been fitted to a measured oxygen tension — so this difference is a statement about the sink constant and the cell density, not about the two tumours.
2. **ECM.** The IO panel carries only 3 of the 11 matrix transcripts the ECM field is built from (VCAN, DCN, FBLN1) against section 1's 11 (COL1A2, COL12A1, COL25A1, POSTN, VCAN, DCN, FBLN1, LAMA2, THBS1, SPON1, TGFBI) — **COL1A2 and POSTN, the two fibrillar-collagen markers, are absent from the IO panel entirely.** Mean ECM comes out at 0.457 vs 0.076, which scales D_eff to 0.76× of its section-1 value. That is **panel content, not tissue biology** — the ECM field is a max-normalised KDE over whichever matrix markers happen to be measured, so an 11-marker field and a 3-marker field have no common absolute scale even on identical tissue. This is the clearest case in the report of a difference that must NOT be read as biology.

**Consequence.** Absolute kills and absolute AUCs are *not* comparable between the two sections. The **ordering** and the **within-section spread** are the transportable quantities, and those are what §3.3 reports.

## 4. The hypoxia check, on its home turf

The 319-gene brain panel measures **no** hypoxia-response gene — no VEGFA, no HIF1A, no CA9, no ADM, no NDRG1, no SLC2A1 — so on section 1 hypoxia is **unmeasured**, and this module declines to produce a number for it rather than substituting a proxy. The IO panel measures VEGFA and HIF1A, so section 2 is the *only* place this test can run, and it is the section `results/state_validation.md` Test 1 was found on. This is therefore a **re-derivation on the same data**, not an out-of-sample test — its value is that it is computed by a different module through the same production path, so it either confirms or refutes that report arithmetically.

Measured genes: **VEGFA** (detected in 63.36% of cells), **HIF1A** (detected in 54.93% of cells). Their mutual Spearman ρ is **+0.0393** — close to zero, so they are *not* two readouts of one latent hypoxia variable, and the composite is reported alongside the per-gene results, never instead of them.

### 4.1 Every predictor, sign-aligned to oxygen (negative = correct)

| predictor | what it is | ρ vs hypoxia score | ρ vs VEGFA | ρ vs HIF1A | ρ with log-depth partialled out |
|---|---|---:|---:|---:|---:|
| **SIMULATED O₂ (the PDE)** | the reaction-diffusion solve under test | **-0.1135** | -0.1464 | +0.0067 | -0.1283 |
| vascular_kde | the solver's own input field, undiffused | **-0.1504** | -0.1386 | -0.0755 | -0.1589 |
| dist_endothelial | one cKDTree query on the PDE's own source cells; no physics | **-0.1503** | -0.1647 | -0.0290 | -0.1545 |
| local_cell_density | pure local crowding; ignores vessels entirely | **-0.0192** | -0.0031 | -0.0439 | -0.0336 |

Permutation p-values (label shuffling, 2,000 draws on a seeded 50,000-cell subsample — the kernel is O(n × draws); every effect size above is on **all** 816,769 cells, so subsampling the null can only make a p-value more conservative):

| predictor | ρ on the subsample | permutation p |
|---|---:|---:|
| simulated_o2 | -0.1088 | 0.00050 |
| dist_endothelial | -0.1549 | 0.00050 |
| vascular_kde | -0.1511 | 0.00050 |

**CONFIRMED.** The simulated field is significantly anti-correlated with measured hypoxia response (ρ = -0.1135), but a one-line `cKDTree` query on its own source cells scores -0.1503 — 1.32× the solved field — and the undiffused source term itself scores -0.1504, so the reaction–diffusion solve is not what carries the signal. `results/state_validation.md` §1.2 is reproduced.

### 4.2 Nested test — does the PDE add anything beyond the raw geometry?

| partial Spearman | ρ |
|---|---:|
| hypoxia score vs **simulated O₂**, distance-to-vessel held fixed | -0.0631 |
| hypoxia score vs **−distance-to-vessel**, simulated O₂ held fixed | -0.1174 |
| simulated O₂ vs −distance-to-vessel (how nested the two are) | +0.3689 |

**Depth control.** The simulated field correlates with total transcript count — a purely technical variable — at ρ = +0.0952; partialling log-depth out (last column of §4.1) does not change the ordering, so it is not a depth artefact.

### 4.3 The 480-gene ranking

Every one of the 480 panel genes was depth-normalised and rank-correlated against the simulated oxygen at each of the 816,769 cells.

| rank (1 = most negative) | gene | ρ vs simulated O₂ | detection |
|---:|---|---:|---:|
| 1 | VEGFA **←** | -0.1464 | 63.36% |
| 2 | PTN | -0.1241 | 94.57% |
| 3 | SOX9 | -0.1238 | 79.41% |
| 4 | NCAM1 | -0.1150 | 83.24% |
| 5 | S100B | -0.1045 | 96.95% |
| 6 | IL1RAP | -0.1022 | 45.74% |
| 7 | FEZ1 | -0.0894 | 83.21% |
| 8 | CCND1 | -0.0811 | 77.43% |
| … | … | … | … |
| 233 | HIF1A **←** | +0.0067 | 54.93% |
| … | … | … | … |
| 480 | S100A9 | +0.1849 | 7.67% |
| 479 | A2M | +0.1674 | 51.44% |
| 478 | RGS5 | +0.1541 | 12.45% |
| 477 | SPARC | +0.1511 | 78.80% |
| 476 | FLT1 | +0.1335 | 9.42% |

* **VEGFA is rank 1 of 480** at ρ = -0.14641. Of every gene the panel measures, the canonical hypoxia-response gene is (among) the ones the simulated oxygen field most strongly anti-correlates with. Reproduced.
* **HIF1A is rank 233 of 480** at ρ = +0.00667 — null. Expected on mechanism: HIF-1α is regulated by oxygen-dependent prolyl hydroxylation and VHL-mediated proteolysis of the **protein** (Ivan et al., *Science* 2001;292:464-468, [PMID 11292862](https://pubmed.ncbi.nlm.nih.gov/11292862/); Jaakkola et al., *Science* 2001;292:468-472, [PMID 11292861](https://pubmed.ncbi.nlm.nih.gov/11292861/)), so its mRNA has no obligation to track oxygen. Reproduced.
* **The circularity, restated.** The positive end of that table is the vascular and myeloid programme, and `FLT1` is one of the three ENDOTHELIAL markers that *builds* the oxygen source term — so its positive correlation is partly the model correlating with itself. Only the negative end is informative.

### 4.4 Field-vs-field agreement, and the state bands

Averaging both fields over fixed 60 µm tissue tiles (≥5 measured cells per tile, 24,254 tiles) gives Spearman ρ = **-0.2718** (Pearson -0.3155) — about 2.4× the per-cell value, the expected signature of per-cell dropout rather than of a better model. It is the most favourable framing of the test and is reported as such.

State occupancy over all cells at this configuration: **HYPOXIC** 96.01%, **QUIESCENT** 3.99%.

## 5. Everything, side by side

| quantity | Section 1 (development) | Section 2 (HELD OUT) | ratio / note |
|---|---:|---:|---|
| **PROVENANCE** | | | |
| cells | 40,887 | 816,769 | 19.98× |
| panel gene targets | 319 | 480 | 1.50× |
| marker map selected | brain_319 | io_480 | shared genes between the two panels: 60 |
| tissue area (mm²) | 27.9 | 89.5 | 3.21× |
| cell density (/mm²) | 1465 | 9127 | 6.23× |
| **COMPOSITION** | | | |
| TUMOR (%) | 73.80 | 79.09 | 1.07× |
| STROMAL (%) | 7.92 | 9.98 | 1.26× |
| IMMUNE (%) | 8.54 | 6.74 | 0.79× |
| ENDOTHELIAL (%) | 9.74 | 4.18 | 0.43× |
| ENDOTHELIAL vs 8.0% volume-EM | 1.22× | 0.52× | PASS → **LOW** |
| endothelial areal density (/mm²) | 143 | 382 | 2.67× |
| **FIELDS** | | | |
| mean vascular density | 0.0250 | 0.0579 | 2.32× |
| mean ECM density | 0.0759 | 0.4571 | different matrix genes on the two panels — NOT comparable as biology |
| mean oxygen (tissue) | 0.2955 | 0.1232 | absolute level is uncalibrated on both |
| model-HYPOXIC tumour fraction | 12.68% | 96.90% | a statement about the sink constant, not the tumours |
| **VASCULATURE** | | | |
| discrete vessel objects | 3,092 | 9,322 | 3.01× |
| MVD (/mm²) | 110.8 | 104.2 | published 65 (35–95): HIGH → **HIGH** |
| intervessel spacing 1/√MVD (µm) | 95.0 | 98.0 | published 124 (103–169): LOW → **LOW** |
| intervessel NN median (µm) | 40.1 | 31.0 | 0.77× |
| tumour→vessel median (µm) | 43.3 | 39.6 | published 58 (48–79): LOW → **LOW** |
| tumour cells > 100 µm from a vessel (%) | 5.32 | 6.60 | 1.24× |
| **vascular area fraction (%)** | 2.50 | 5.79 | measured **2.1 ± 0.3**: HIGH (1.19×) → **HIGH** (2.76×) |
| endothelial cells per vessel object | 1.29 | 3.67 | section 2's objects look far more like real vessel profiles |
| **DRUG PANEL (uniform sensitivity)** | | | |
| temozolomide — coverage @1× IC50 | 0.0000 | 0.0000 |  |
| lomustine — coverage @1× IC50 | 0.2175 | 0.8013 | 3.68× |
| cisplatin — coverage @1× IC50 | 0.9936 | 0.9976 | 1.00× |
| methotrexate — coverage @1× IC50 | 0.0000 | 0.0000 |  |
| doxorubicin — coverage @1× IC50 | 0.3836 | 0.7908 | 2.06× |
| temozolomide — mean tumour AUC (µM·s) | 7.06e+05 | 1.64e+06 | 2.32× |
| lomustine — mean tumour AUC (µM·s) | 8.7e+05 | 1.95e+06 | 2.24× |
| cisplatin — mean tumour AUC (µM·s) | 1.69e+06 | 3.85e+06 | 2.28× |
| methotrexate — mean tumour AUC (µM·s) | 1.48e+04 | 3.25e+04 | 2.20× |
| doxorubicin — mean tumour AUC (µM·s) | 1.08e+04 | 2.36e+04 | 2.18× |
| temozolomide — kill | 0.0147 | 0.0274 | 1.86× |
| lomustine — kill | 0.1975 | 0.2210 | 1.12× |
| cisplatin — kill | 0.3699 | 0.3162 | 0.85× |
| methotrexate — kill | 0.0162 | 0.0276 | 1.71× |
| doxorubicin — kill | 0.2165 | 0.2284 | 1.06× |
| kill ordering | cisplatin > doxorubicin > lomustine > methotrexate > temozolomide | cisplatin > doxorubicin > lomustine > methotrexate > temozolomide | **IDENTICAL** |
| coverage ordering | cisplatin > doxorubicin > lomustine > temozolomide > methotrexate | cisplatin > lomustine > doxorubicin > temozolomide > methotrexate | **DIFFERENT** |
| ρ(kill, AUC/IC50) | +1.00 | +1.00 | delivery × potency |
| **HYPOXIA** | | | |
| hypoxia-response genes on panel | none | VEGFA, HIF1A | the 319 panel cannot test this at all |
| ρ(simulated O₂, hypoxia score) | n/a | -0.1135 |  |
| ρ(dist-to-endothelial, hypoxia score) | n/a | -0.1503 | **beats the PDE** |
| VEGFA rank of all panel genes | n/a | 1/480 |  |
| HIF1A rank of all panel genes | n/a | 233/480 |  |

## 6. Verdict — what generalises to a second patient, and what does not

### 6.1 Replicates / confirmed

Two of these are confirmations of a *negative* result — the hypoxia rows confirm that the reaction-diffusion solve is out-predicted by a distance transform. 'Confirmed' means the prior finding reproduced, not that the model performed well.

* **Kill ordering is cisplatin > doxorubicin > lomustine > methotrexate > temozolomide** — identical on both sections
* **therapeutic_coverage @1x IC50 is still a DISCRIMINATING headline metric** — coverage still spreads the panel on both sections
* **The ordering is NOT simply the IC50 ordering read back out** — potency alone would rank doxorubicin > methotrexate > cisplatin > lomustine > temozolomide; the simulation ranks cisplatin > doxorubicin > lomustine > methotrexate > temozolomide. rho(kill, 1/IC50) = +0.50, rho(kill, AUC/IC50) = +1.00
* **Simulated oxygen is BEATEN by a distance-to-nearest-vessel baseline (state_validation.md Test 1)** — the one-line cKDTree baseline scores 1.32x the solved field
* **VEGFA is the most anti-correlated gene on the panel; HIF1A is null** — VEGFA rho -0.1464, HIF1A rho +0.0067

### 6.2 Does NOT generalise

* **ENDOTHELIAL fraction matches human volume-EM 8.0% (PMID 38723085)** — 4.18% is 1.91x from 8.0% (outside the declared 1.5x tolerance); section 1 was 9.74% (1.22x)
* **Vascular area fraction matches measured 2.1 +/- 0.3% (PMID 36611806)** — 5.79% vs 2.1% is 2.76x and lies OUTSIDE the published 1.8-2.4% band; section 1 gave 2.50% (1.19x, also outside but only just)

### 6.3 Fails on BOTH sections — not a generalisation failure

These are the rows where the held-out patient behaved exactly like the development section, and both disagree with something this project has claimed. Reading them as 'it did not generalise' would be wrong; they are disagreements with a *prior configuration*, and the second patient's role here is to show the disagreement is not patient-specific.

* **MVD lands in the published human GBM range 35-95/mm2 (PMID 27583442)** — 104/mm2 vs published 65 (35-95); section 1 gave 111/mm2. BOTH sections are high, and by a similar amount, so this is a systematic property of the segmentation (endothelial cells are over-split into vessel objects) rather than something the held-out patient revealed — and it is the one place where the two sections AGREE with each other while both disagreeing with the literature.
* **therapeutic_coverage @1x IC50 recovers the kill ordering** — coverage does NOT reproduce the kill ordering on EITHER section in this configuration (both swap doxorubicin and lomustine), so this is not a held-out failure — it is a disagreement with `results/delivery_metrics.md`, which ran at a different exposure scale where coverage was not saturated. See §3.3

### 6.4 Untestable on one section or the other

* **Hypoxia against measured markers** is impossible on the 319-gene section: it carries no VEGFA, no HIF1A, no CA9, no ADM, no NDRG1, no SLC2A1. The honest statement for section 1 remains that hypoxia is *unmeasured*, not that it has been simulated.
* **The TOP2A-based doxorubicin sensitivity determinant** cannot be evaluated on section 2 — TOP2A is absent from the IO panel. This is why the whole comparison runs in `uniform` mode; `determinant` mode is not even definable in a comparable way across these two panels.
* **MGMT** — and therefore any measured temozolomide sensitivity — is absent from **both** panels. TMZ's kill is an upper bound on both sections for the same reason.
* **Outcome.** Neither section has a drug-response label. Nothing here validates a kill prediction; the drug panel tests whether the machinery behaves consistently on a second real tissue, which is a different and much weaker claim.

### 6.5 The one-paragraph reading

**The response layer travels; the reconstruction layer does not.** The drug ordering and the delivery × potency structure behind it survive a second real patient on a different panel — but they were always going to be robust, because they are dominated by drug properties (IC50, BBB penetrance, molecular weight) that are identical on both sections, and the tissue only modulates them. The quantities that are genuinely *about the tissue* — the endothelial cell-type fraction and the vascular area fraction, i.e. the two numbers this project cites as its external validations — are the ones that fail to replicate. That is the honest asymmetry, and it points at where the next work is: the reconstruction's external anchors are section-specific, and the two constants that produce them (the typing threshold and the 7.5 µm vessel radius) are entangled with the grid and with panel sparsity in a way that a single section could never expose.

## 7. Caveats — what could make these conclusions wrong

1. **Different patients AND different panels.** Only 60 genes are shared. Every composition difference is confounded by panel content; §1.2 and §3.5 name the specific places where a difference is more plausibly measurement than biology (the endothelial gate arity, the ECM marker set).
2. **n = 2.** Two sections cannot separate 'this pipeline does not generalise' from 'these two tumours differ'. A claim that a quantity *fails* to replicate is well supported by a single discrepancy; a claim that one *does* replicate rests on a single agreement and is much weaker.
3. **The IO marker map is not fully blind** (§1.2): it was written against this section, though never against an external anchor.
4. **Grid convergence is not re-established here.** dx = 12 µm was shown converged on section 1's geometry (`results/numerical_rigor.md` §2). Section 2 is denser and its vessel objects are larger; the convergence argument is inherited, not re-run. A dx sweep on section 2 is the obvious next check and it was not affordable in this run.
5. **The absolute oxygen level is uncalibrated on both sections**, so the 97% vs 13% hypoxic-fraction gap is a statement about `OxygenConfig.consumption` and cell density, not about the two tumours. Every absolute kill inherits it.
6. **Typing circularity.** `FLT1` and `PLVAP` are IO-panel ENDOTHELIAL markers and therefore help build the oxygen source, so the positive tail of the gene ranking is partly the model correlating with its own inputs. Only the negative tail is informative.
7. **Pericytes are unlabelled on section 2.** `CellType` has no MURAL member, so RGS5/PDGFRB/NOTCH3-positive mural cells fall into the ungated TUMOR default, inflating the tumour count and hence the oxygen sink. That is part of why section 2 is so much more hypoxic.
8. **2D, in-plane, single section** on both. Distances are Euclidean and in-plane, not geodesic and not 3D.
9. **Subsampling.** Only the permutation null is subsampled (50,000 of 816,769 cells, 2,000 draws); every effect size, fraction and simulation number is computed on all cells. Plot (c) decimates the CDF for rendering only.

## 8. Reproduction

```bash
./.venv/bin/python -m tmesim.validation.heldout_patient
```

Deterministic (`seed = 0` in every `SimConfig`; permutation seed 0). Runtime 34 min, dominated by 10 PDE solves — 28 min of it on section 2's 905×956 grid.

**Files added by this validation (nothing else was modified):**

* `tmesim/validation/heldout_patient.py` — this analysis
* `results/heldout_patient.md`, `results/heldout_patient.png`

**Data provenance.** Both sections are open 10x Genomics Xenium datasets, CC BY 4.0, no login or DUA:

* Section 1 — 10x Xenium Human Brain — Glioblastoma section (v1.3.0); panel hBrain_254g base + 65-gene add-on (319 gene targets)
* Section 2 — 10x Xenium Human Brain Cancer (GBM) FFPE — IO panel (v2.0.0); panel Human Immuno-Oncology Profiling base (hImmune_v1, 380 genes) + Add-on B (100 genes) = 480 gene targets

**Figure.** `heldout_patient.png`.
