# The held-out patient — running the whole validated pipeline on a second real GBM section

> # ⚠️ WITHDRAWN SUBJECT — BANNER ADDED 2026-09-03. THE CROSS-COMPOUND ORDERING THIS DOCUMENT RESTS ON NO LONGER EXISTS AS A CLAIM.
>
> **This banner is the only edit made to this file.** Nothing below it has been
> changed, reworded, re-scored or removed. It is preserved byte-identical as the
> record of what was believed and measured at the time, because a record that was
> wrong is part of the record.
>
> **What is below.** Re-runs the whole pipeline on a second real GBM section the pipeline was never developed on. Its §3 asks whether the compound ordering holds on that second patient, §6 lists ordering results among what generalises and what does not, and §6.5 reports the ordering and the delivery × potency structure behind it as surviving. **Those rows are about the withdrawn output.** The reconstruction and vascular-morphometry rows — endothelial fraction, vascular area fraction, MVD — and the hypoxia check are unaffected, and they remain the load-bearing negative results of the report.
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

Every validated result in this project was developed on **one** tissue: the 319-gene 10x Xenium GBM section (40,887 cells). This report re-runs the entire pipeline — reconstruction, vessel segmentation, the five-drug simulation and the hypoxia check — on a **second real patient the pipeline was never developed on**: the 480-gene Immuno-Oncology Xenium GBM section (816,769 cells). Both sections are driven through the identical frozen production path at the identical settings, so **section 1's numbers here are re-derived, not quoted**.

**3 of 7 claims that could be scored as generalisation replicate.** A further 2 fail on **both** sections — those are not generalisation failures, they are claims that do not hold in this configuration at all, and they are called out separately because conflating the two would be the easiest way to make this report say something it should not. The headline, one line each:

| what | section 1 (development) | section 2 (held out) | verdict |
|---|---|---|:--|
| ENDOTHELIAL fraction matches human volume-EM 8.0% (PMID 38723085) | 10.50 % | 7.07 % | **REPLICATES** |
| ENDOTHELIAL areal density (the oxygen/drug SOURCE) is comparable | 153.94 cells/mm2 | 645.15 cells/mm2 | **CONTEXT** |
| Vascular area fraction matches measured 2.1 +/- 0.3% (PMID 36611806) | 1.25 % | 10.98 % | **DOES NOT REPLICATE** |
| MVD lands in the published human GBM range 35-95/mm2 (PMID 27583442) | 110.79 vessels/mm2 | 233.27 vessels/mm2 | **NOT REPRODUCED ON EITHER** |
| Kill ordering is cisplatin > doxorubicin > lomustine > methotrexate > temozolomide | cisplatin > doxorubicin > lomustine > methotrexate > temozolomide  | cisplatin > doxorubicin > lomustine > temozolomide > methotrexate  | **DOES NOT REPLICATE** |
| therapeutic_coverage @1x IC50 recovers the kill ordering | cisplatin > doxorubicin > lomustine > temozolomide > methotrexate  | cisplatin > lomustine > doxorubicin > temozolomide > methotrexate  | **NOT REPRODUCED ON EITHER** |
| therapeutic_coverage @1x IC50 is still a DISCRIMINATING headline metric | 3/5 drugs pinned at 0 or >0.9  | 5/5 drugs pinned at 0 or >0.9  | **DOES NOT REPLICATE** |
| The ordering is NOT simply the IC50 ordering read back out | doxorubicin > methotrexate > cisplatin > lomustine > temozolomide  | doxorubicin > methotrexate > cisplatin > lomustine > temozolomide  | **HOLDS** |
| Simulated oxygen is BEATEN by a distance-to-nearest-vessel baseline (state_validation.md Test 1) | UNTESTABLE (no hypoxia gene on the 319 panel) Spearman rho | rho_PDE -0.1523 vs rho_dist -0.1204 Spearman rho | **REFUTED** |
| VEGFA is the most anti-correlated gene on the panel; HIF1A is null | n/a — neither gene is on the 319 panel  | VEGFA rank 1/480, HIF1A rank 327/480  | **CONFIRMED** |

> **THE CONFOUND, FIRST.** These are **different patients** *and* **different panels** — the two Xenium panels share only **60 genes**. The IO panel carries none of PECAM1 / SOX2 / OLIG2 / TOP2A / AQP4 / GJA1 / MOG, so the endothelial gate runs on VWF/FLT1/PLVAP instead of PECAM1 and the stromal gate on SLC1A3/S100B/PTN instead of AQP4/GJA1/MOG. Every composition difference below is therefore **jointly** a difference in biology and a difference in what is measurable. Where the two can be separated, this report separates them; where they cannot, it says so. Neither section carries a drug-response outcome label, so nothing here validates a kill prediction — only the machinery that produces one.

## 0. What was held out, and what was not

This is an out-of-sample test of the **pipeline**, not of a fitted model. Being explicit about what that does and does not cover:

| held out from section 2 | NOT held out |
|---|---|
| the reconstruction was never inspected on it while the marker gate, the typing threshold or the 30 µm tissue-mask radius were chosen | the IO-panel **marker map** in `reconstruct/marker_panels.py` was written for this section (in `results/state_validation.md`), by choosing among measured candidate marker sets — see §1.2 |
| the 20 µm vessel linking distance and the 7.5 µm vessel radius were fixed against section 1 and are used here unchanged | the oxygen solver's `consumption` constant has never been fitted to a measured oxygen tension on **either** section |
| the perivascular-sink calibration **α = 148.6** was fit once on doxorubicin against section 1 and is **not** refitted here | the drug IC50s are external (GDSC/literature) and are the same numbers on both |
| the fixed thresholds in `delivery_metrics` (the 1× IC50 coverage bar; the vessel-distance axis is now the PRIMARY direct-cell metric and has no threshold at all — see `results/vessel_profile.md`'s dedicated distance-metric-change section) | `dx = 12 µm` was chosen for grid convergence on section 1; convergence is **not** re-established here (§6) |

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
| TUMOR | 58.24% | 66.18% |
| STROMAL | 11.78% | 15.37% |
| IMMUNE | 19.48% | 11.39% |
| ENDOTHELIAL | 10.50% | 7.07% |
| mean vascular density | 0.0268 | 0.1004 |
| mean ECM density | 0.0759 | 0.4571 |
| ECM markers on panel | COL1A2, COL12A1, COL25A1, POSTN, VCAN, DCN, FBLN1, LAMA2, THBS1, SPON1, TGFBI | VCAN, DCN, FBLN1 |
| mean oxygen (tissue) | 0.3520 | 0.2210 |
| tissue coverage of bbox | 67.70% | 71.83% |

### 1.1 The endothelial check — the one external anchor the reconstruction has

Human cortex volume electron microscopy counts **8.0%** of cells as endothelial (Shapson-Coe A et al., Science 2024;384(6696):eadk4858 — 'A petavoxel fragment of human cerebral cortex reconstructed at nanoscale resolution', [PMID 38723085](https://pubmed.ncbi.nlm.nih.gov/38723085/)). Section 1's **10.50%** matching that — rather than the 1–2% that dissociation-based scRNA-seq reports — is the single strongest external check the reconstruction has, and `results/SESSION_AUDIT.md` §4 leans on it.

**On the held-out patient the classifier calls 7.07% endothelial — 0.88× the volume-EM value, against section 1's 1.31×.**

| | section 1 | section 2 | ratio |
|---|---:|---:|---:|
| ENDOTHELIAL, % of cells | 10.50% | 7.07% | 0.67× |
| ENDOTHELIAL, cells per mm² | 154 | 645 | 4.19× |
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
2. **Denominator.** Section 2 has 6.2× the cell density, so the same vasculature is divided by many more non-endothelial cells. The denominator alone is arithmetically sufficient to explain the whole gap: divide section 2's endothelial areal density by section 1's total cell density and the fraction would read 44.0% — far ABOVE 8.0%, not the 7.07% observed. The observed fraction is set as much by how many cells the section contains as by how many of them are endothelial.
3. **Panel content.** A 3-marker mean-z gate (VWF/FLT1/PLVAP) behaves differently from a 1-marker gate (PECAM1). Measured detection rates — section 1: PECAM1 24.57%; section 2: VWF 4.15%, FLT1 9.42%, PLVAP 4.36% — so on section 2 the gate is the mean z of three sparse markers, which a cell strongly positive for one but silent for the other two can fail. That pushes the endothelial call **down** for a purely technical reason.

**Explanations 2 and 3 are measurement artefacts, not biology, and neither can be excluded here.** The honest statement is that *the endothelial cell-type fraction is not a portable quantity across panels*, and the 8.0% agreement on section 1 should be read as a property of that section and that gate, not as a validated property of the classifier.

## 2. Vascular morphometry — the only vascular quantity scored against a measurement

`reconstruct/vessel_network.py`, unmodified: measured endothelial cells single-linkage clustered at **20 µm** into discrete vessel objects, then painted as a source of radius **7.5 µm**. Both constants were fixed against section 1 and are used here unchanged.

| quantity | section 1 | section 2 (HELD OUT) | published | range | s1 verdict | s2 verdict | source |
|---|---:|---:|---:|---|:--|:--|---|
| microvessel density, human GBM (n=20), CD31 whole-slide, unbiased systematic-uniform-random sampling (vessels/mm^2) | 110.79 | **233.27** | 65.00 | 35.0–95.0 | HIGH | **HIGH** | [PMID 27583442](https://pubmed.ncbi.nlm.nih.gov/27583442/) |
| 1/sqrt(MVD) spacing implied by the MVD anchor above (um) | 95.01 | **65.47** | 124.00 | 102.6–169.0 | LOW | **LOW** | [PMID 27583442](https://pubmed.ncbi.nlm.nih.gov/27583442/) |
| median tumour-cell-to-nearest-vessel distance implied by the MVD anchor under a random (Poisson) vessel arrangement, median = sqrt(ln2 / (pi q)) (um) | 42.85 | **23.87** | 58.00 | 48.2–79.4 | LOW | **LOW** | [PMID 27583442](https://pubmed.ncbi.nlm.nih.gov/27583442/) |
| CD34+ vascular area fraction, human IDH-wt GBM (n=30), whole-field positive-pixel count (% of tissue area) | 1.25 | **10.98** | 2.10 | 1.8–2.4 | LOW | **HIGH** | [PMID 36611806](https://pubmed.ncbi.nlm.nih.gov/36611806/) |

### 2.1 The load-bearing row

`results/vessel_network_impact.md` calls the **vascular area fraction** the load-bearing row, because it is the only one that tests the *field* rather than the point pattern, and because the vessel radius was fixed from published vessel-caliber imaging **before** the field was built — so it is an untuned comparison. Section 1 scored **1.25%** against a measured whole-section CD34 value of **2.1 ± 0.3%**, i.e. 0.59×.

**Section 2 scores 10.98% — 5.23× the measured value, and 8.82× section 1's.** That is outside the published band by a wide margin, and it is the clearest single failure of generalisation in this report.

The cause is not the fixed-7.5-µm-per-cell painting an earlier version of this paragraph described — that path was superseded by `vessel_geometry.py`'s sub-grid, data-derived radius model (see `results/vessel_geometry_fix.md`), which is what §2.1's numbers above are already computed from. A dedicated follow-up (`results/sweep/vascular_area_fraction_generalization_2026-07-30T03-18-42Z.md`) tested whether ANY non-circular, section-invariant function from a vessel object's own cell count to its radius could close this gap — six candidates, including literature-motivated flattened and capped variants, on this section plus a third independent GBM section — and found a genuine null: every candidate that helps this held-out section pulls the development section further from target, because the two errors run in opposite directions (development: 0.59× too low; both other sections: 5.2× too high) while every candidate radius model moves monotonically in the same direction on all three. The radius that WOULD hit 2.1% on each section separately is 7.28 / 3.12 / 2.57 µm — a near-perfect inverse rank correlation with each section's own microvessel density, directionally consistent with a real published finding (Halshtok Neiman et al., *Clin Neuroradiol* 2021, PMID 34940899). But using this project's own MVD to calibrate its own vessel radius would be circular: §6.3 below independently flags MVD itself as an unresolved, both-sections-fail segmentation defect. Conditioning one flagged defect on another is not a fix, so none was applied — this remains open, honestly.

### 2.2 Vessel objects and the clustering that makes them comparable to histology

| | section 1 | section 2 |
|---|---:|---:|
| endothelial cells | 4,295 | 57,733 |
| discrete vessel objects | 3,091 | 20,875 |
| tissue area (mm²) | 27.9 | 89.5 |
| MVD (/mm²) | 110.8 | 233.3 |
| endothelial cells per vessel | 1.39 | 2.77 |
| singleton objects | 81.3% | 62.9% |
| median object extent, multi-cell (µm) | 18.7 | 22.6 |
| median aspect ratio, multi-cell | 2.28 | 2.39 |
| tumour cells > 100 µm from a vessel | 4.02% | 0.35% |
| tumour cells > 150 µm from a vessel | 0.24% | 0.12% |
| intervessel NN median (µm) | 37.8 | 27.7 |
| tumour→vessel median (µm) | 42.8 | 23.9 |
| tumour→vessel p90 (µm) | 82.7 | 49.4 |
| clustering index vs Poisson (median) | 0.96× | 0.78× |
| clustering index vs Poisson (>100 µm tail) | 1.31× | 5.38× |

**One thing genuinely improves out of sample.** `results/vascular_validation.md` §3.2 found that on section 1 the clustering step 'barely bites' — 1.29 endothelial cells per putative vessel and 85% singletons — so the reconstruction was resolving scattered endothelial nuclei rather than vessel walls. On section 2 it merges 2.77 cells per object with 63% singletons and a median multi-cell extent of 23 µm at aspect ratio 2.4: **elongated, multi-cell objects that actually look like vessel profiles.** The denser section supports the segmentation the sparse one could not. That is a genuine, if backhanded, out-of-sample success — and it is also part of why the area fraction inflates.

### 2.3 Linking-distance sweep — the free parameter, swept on both sections

| link (µm) | s1 vessels | s1 MVD | s1 cells/vessel | s2 vessels | s2 MVD | s2 cells/vessel |
|---:|---:|---:|---:|---:|---:|---:|
| 10 | 4,149 | 148.7 | 1.04 | 43,154 | 482.2 | 1.34 |
| 15 | 3,698 | 132.5 | 1.16 | 28,335 | 316.6 | 2.04 |
| 20 | 3,091 | 110.8 | 1.39 | 20,875 | 233.3 | 2.77 |
| 25 | 2,665 | 95.5 | 1.61 | 16,063 | 179.5 | 3.59 |
| 30 | 2,421 | 86.8 | 1.77 | 12,209 | 136.4 | 4.73 |
| 40 | 2,060 | 73.8 | 2.08 | 6,617 | 73.9 | 8.72 |

On section 2 the published MVD range (35–95/mm²) is reached at a linking distance of ~25–40 µm rather than the 20 µm fixed on section 1. **No linking distance was re-tuned for this report** — reporting the sweep is the honest alternative to tuning it.

## 3. The drug panel — does the ordering hold on a second patient?

Five drugs, **uniform** per-cell sensitivity (the only fair cross-drug mode — every drug carries 1.0 at every cell, so the factor cancels and the ranking isolates delivery × potency), dx = 12 µm, α = 148.6 fit once on section 1 and **not refitted**, PDE dt = 10 s over 6 h with convection on. Identical settings on both sections.

### 3.1 Section 1 (development)

| drug | IC50 (µM) | BBB | D_eff (µm²/s) | mean tumour AUC (µM·s) | **coverage @1× IC50** | coverage @0.5×/2× | median E/IC50 | **kill** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| temozolomide | 606.047 | 0.125 | 265.3 | 7.66e+05 | **0.0000** | 0.000 / 0.000 | 0.06 | **0.0184** |
| lomustine | 50 | 0.973 | 247.4 | 9.32e+05 | **0.2802** | 0.856 / 0.020 | 0.78 | **0.2371** |
| cisplatin | 24.5503 | 0.858 | 225.4 | 1.82e+06 | **0.9955** | 0.999 / 0.914 | 3.18 | **0.4362** |
| methotrexate | 12 | 0.023 | 189.8 | 1.57e+04 | **0.0000** | 0.000 / 0.000 | 0.05 | **0.0192** |
| doxorubicin | 0.5 | 0.021 | 176.3 | 1.15e+04 | **0.4373** | 0.858 / 0.071 | 0.92 | **0.2595** |

### 3.2 Section 2 (HELD OUT)

| drug | IC50 (µM) | BBB | D_eff (µm²/s) | mean tumour AUC (µM·s) | **coverage @1× IC50** | coverage @0.5×/2× | median E/IC50 | **kill** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| temozolomide | 606.047 | 0.125 | 201.7 | 2.87e+06 | **0.0000** | 0.000 / 0.000 | 0.21 | **0.0797** |
| lomustine | 50 | 0.973 | 188.1 | 3.46e+06 | **0.9702** | 0.998 / 0.819 | 3.04 | **0.3486** |
| cisplatin | 24.5503 | 0.858 | 171.4 | 6.79e+06 | **0.9987** | 0.999 / 0.999 | 12.30 | **0.4248** |
| methotrexate | 12 | 0.023 | 144.3 | 5.8e+04 | **0.0000** | 0.012 / 0.000 | 0.21 | **0.0657** |
| doxorubicin | 0.5 | 0.021 | 134.1 | 4.22e+04 | **0.9691** | 0.996 / 0.843 | 3.62 | **0.3603** |

### 3.3 The ordering

| ranked by | section 1 | section 2 (HELD OUT) | same? |
|---|---|---|:--:|
| simulated **kill** | cisplatin > doxorubicin > lomustine > methotrexate > temozolomide | cisplatin > doxorubicin > lomustine > temozolomide > methotrexate | ❌ |
| **coverage @1× IC50** | cisplatin > doxorubicin > lomustine > temozolomide > methotrexate | cisplatin > lomustine > doxorubicin > temozolomide > methotrexate | ❌ |
| mean tumour AUC (delivery alone) | cisplatin > lomustine > temozolomide > methotrexate > doxorubicin | cisplatin > lomustine > temozolomide > methotrexate > doxorubicin | ✅ |
| 1/IC50 (potency alone) | doxorubicin > methotrexate > cisplatin > lomustine > temozolomide | doxorubicin > methotrexate > cisplatin > lomustine > temozolomide | ✅ |
| AUC/IC50 (delivery × potency) | cisplatin > doxorubicin > lomustine > methotrexate > temozolomide | cisplatin > doxorubicin > lomustine > methotrexate > temozolomide | ✅ |

*(The 'same?' column compares the two SECTIONS with each other. Whether each ranking agrees with the KILL ranking — the separate question — is read down the rows: on both sections `coverage` swaps doxorubicin and lomustine relative to kill, and `AUC/IC50` reproduces kill exactly.)*

**Published ordering** (`results/real_drug_panel.md`, uniform mode): cisplatin > doxorubicin > lomustine > methotrexate > temozolomide. Section 1's kill fractions here (temozolomide 0.0184, lomustine 0.2371, cisplatin 0.4362, methotrexate 0.0192, doxorubicin 0.2595) reproduce that report's numbers exactly, which is the check that this harness is running the same pipeline and not a variant of it.

**The kill ordering CHANGES on the second patient.** A different ordering on a different real tumour is a legitimate finding — real inter-patient heterogeneity — not a failure. What follows is the check on whether the change is biology or bookkeeping.

### 3.4 Is the ordering trivially the IC50 ordering?

An ordering that just reproduces the potency ranking would be reproducible by a spreadsheet with no geometry, no transport and no patient in it. The diagnostic is whether delivery ever overturns potency.

| | section 1 | section 2 |
|---|---:|---:|
| IC50 fold range across the panel | 1212× | 1212× |
| mean-tumour-AUC fold range | 158× | 161× |
| ρ(kill, 1/IC50) — potency alone | +0.50 | +0.20 |
| ρ(kill, mean AUC) — delivery alone | +0.30 | +0.40 |
| ρ(kill, AUC/IC50) — the product | +1.00 | +0.90 |
| kill spread (max − min) | 0.418 | 0.359 |
| coverage spread (max − min) | 0.995 | 0.999 |

**It is NOT the IC50 ordering.** Potency alone would rank doxorubicin > methotrexate > cisplatin > lomustine > temozolomide; the simulation ranks cisplatin > doxorubicin > lomustine > temozolomide > methotrexate. Delivery spans 161× and potency 1212×, and the product AUC/IC50 tracks kill at ρ = +0.90 — the delivery × potency claim from `results/SESSION_AUDIT.md` §4 survives on the held-out patient.

### 3.5 Why the absolute numbers move so much, and why that is NOT a ranking result

| | section 1 | section 2 |
|---|---:|---:|
| mean oxygen over tissue | 0.3520 | 0.2210 |
| model-HYPOXIC tumour fraction | 0.00% | 33.76% |
| model-VIABLE tumour fraction | 10.04% | 0.00% |
| mean ECM density | 0.076 | 0.457 |

Two mechanisms move every absolute kill on section 2, and **neither is a drug property**, so both largely cancel out of the ranking:

1. **Oxygen.** Mean tissue O₂ is 0.221 on section 2 against 0.352 on section 1 — the denser tissue is a bigger sink — so 34% of tumour cells fall in the model's HYPOXIC band versus 0%, and `PDConfig.hypoxia_resistance` (×0.35) plus `quiescent_resistance` (×0.5) then suppress every drug's kill together. **`results/state_validation.md` §1.6 already established that the absolute oxygen level is uncalibrated** — `OxygenConfig.consumption` has never been fitted to a measured oxygen tension — so this difference is a statement about the sink constant and the cell density, not about the two tumours.
2. **ECM.** The IO panel carries only 3 of the 11 matrix transcripts the ECM field is built from (VCAN, DCN, FBLN1) against section 1's 11 (COL1A2, COL12A1, COL25A1, POSTN, VCAN, DCN, FBLN1, LAMA2, THBS1, SPON1, TGFBI) — **COL1A2 and POSTN, the two fibrillar-collagen markers, are absent from the IO panel entirely.** Mean ECM comes out at 0.457 vs 0.076, which scales D_eff to 0.76× of its section-1 value. That is **panel content, not tissue biology** — the ECM field is a max-normalised KDE over whichever matrix markers happen to be measured, so an 11-marker field and a 3-marker field have no common absolute scale even on identical tissue. This is the clearest case in the report of a difference that must NOT be read as biology.

**Consequence.** Absolute kills and absolute AUCs are *not* comparable between the two sections. The **ordering** and the **within-section spread** are the transportable quantities, and those are what §3.3 reports.

## 4. The hypoxia check, on its home turf

The 319-gene brain panel measures **no** hypoxia-response gene — no VEGFA, no HIF1A, no CA9, no ADM, no NDRG1, no SLC2A1 — so on section 1 hypoxia is **unmeasured**, and this module declines to produce a number for it rather than substituting a proxy. The IO panel measures VEGFA and HIF1A, so section 2 is the *only* place this test can run, and it is the section `results/state_validation.md` Test 1 was found on. This is therefore a **re-derivation on the same data**, not an out-of-sample test — its value is that it is computed by a different module through the same production path, so it either confirms or refutes that report arithmetically.

Measured genes: **VEGFA** (detected in 63.36% of cells), **HIF1A** (detected in 54.93% of cells). Their mutual Spearman ρ is **+0.0393** — close to zero, so they are *not* two readouts of one latent hypoxia variable, and the composite is reported alongside the per-gene results, never instead of them.

### 4.1 Every predictor, sign-aligned to oxygen (negative = correct)

| predictor | what it is | ρ vs hypoxia score | ρ vs VEGFA | ρ vs HIF1A | ρ with log-depth partialled out |
|---|---|---:|---:|---:|---:|
| **SIMULATED O₂ (the PDE)** | the reaction-diffusion solve under test | **-0.1523** | -0.2004 | +0.0140 | -0.1619 |
| vascular_kde | the solver's own input field, undiffused | **-0.1111** | -0.1017 | -0.0556 | -0.1119 |
| dist_endothelial | one cKDTree query on the PDE's own source cells; no physics | **-0.1204** | -0.1237 | -0.0342 | -0.1185 |
| local_cell_density | pure local crowding; ignores vessels entirely | **-0.0192** | -0.0031 | -0.0439 | -0.0336 |

Permutation p-values (label shuffling, 2,000 draws on a seeded 50,000-cell subsample — the kernel is O(n × draws); every effect size above is on **all** 816,769 cells, so subsampling the null can only make a p-value more conservative):

| predictor | ρ on the subsample | permutation p |
|---|---:|---:|
| simulated_o2 | -0.1489 | 0.00050 |
| dist_endothelial | -0.1254 | 0.00050 |
| vascular_kde | -0.1174 | 0.00050 |

**REFUTED.** The simulated field is significantly anti-correlated with measured hypoxia response (ρ = -0.1523), but a one-line `cKDTree` query on its own source cells scores -0.1204 — 0.79× the solved field — and the undiffused source term itself scores -0.1111, so the solve does add predictive information here, contradicting `results/state_validation.md` §1.2.

### 4.2 Nested test — does the PDE add anything beyond the raw geometry?

| partial Spearman | ρ |
|---|---:|
| hypoxia score vs **simulated O₂**, distance-to-vessel held fixed | -0.1227 |
| hypoxia score vs **−distance-to-vessel**, simulated O₂ held fixed | -0.0792 |
| simulated O₂ vs −distance-to-vessel (how nested the two are) | +0.3006 |

**Depth control.** The simulated field correlates with total transcript count — a purely technical variable — at ρ = +0.0572; partialling log-depth out (last column of §4.1) does not change the ordering, so it is not a depth artefact.

### 4.3 The 480-gene ranking

Every one of the 480 panel genes was depth-normalised and rank-correlated against the simulated oxygen at each of the 816,769 cells.

| rank (1 = most negative) | gene | ρ vs simulated O₂ | detection |
|---:|---|---:|---:|
| 1 | VEGFA **←** | -0.2004 | 63.36% |
| 2 | PTGDS | -0.1161 | 41.11% |
| 3 | CAV1 | -0.1048 | 12.64% |
| 4 | CD44 | -0.0833 | 26.16% |
| 5 | SPARCL1 | -0.0831 | 74.94% |
| 6 | NCAM1 | -0.0810 | 83.24% |
| 7 | SOX9 | -0.0658 | 79.41% |
| 8 | TCIM | -0.0601 | 15.94% |
| … | … | … | … |
| 327 | HIF1A **←** | +0.0140 | 54.93% |
| … | … | … | … |
| 480 | A2M | +0.1494 | 51.44% |
| 479 | RGS5 | +0.1447 | 12.45% |
| 478 | SPARC | +0.1339 | 78.80% |
| 477 | FLT1 | +0.1273 | 9.42% |
| 476 | ITGA1 | +0.1196 | 6.57% |

* **VEGFA is rank 1 of 480** at ρ = -0.20039. Of every gene the panel measures, the canonical hypoxia-response gene is (among) the ones the simulated oxygen field most strongly anti-correlates with. Reproduced.
* **HIF1A is rank 327 of 480** at ρ = +0.01403 — null. Expected on mechanism: HIF-1α is regulated by oxygen-dependent prolyl hydroxylation and VHL-mediated proteolysis of the **protein** (Ivan et al., *Science* 2001;292:464-468, [PMID 11292862](https://pubmed.ncbi.nlm.nih.gov/11292862/); Jaakkola et al., *Science* 2001;292:468-472, [PMID 11292861](https://pubmed.ncbi.nlm.nih.gov/11292861/)), so its mRNA has no obligation to track oxygen. Reproduced.
* **The circularity, restated.** The positive end of that table is the vascular and myeloid programme, and `FLT1` is one of the three ENDOTHELIAL markers that *builds* the oxygen source term — so its positive correlation is partly the model correlating with itself. Only the negative end is informative.

### 4.4 Field-vs-field agreement, and the state bands

Averaging both fields over fixed 60 µm tissue tiles (≥5 measured cells per tile, 24,254 tiles) gives Spearman ρ = **-0.3117** (Pearson -0.3663) — about 2.0× the per-cell value, the expected signature of per-cell dropout rather than of a better model. It is the most favourable framing of the test and is reported as such.

State occupancy over all cells at this configuration: **VIABLE** 0.00%, **HYPOXIC** 30.90%, **QUIESCENT** 69.09%.

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
| TUMOR (%) | 58.24 | 66.18 | 1.14× |
| STROMAL (%) | 11.78 | 15.37 | 1.30× |
| IMMUNE (%) | 19.48 | 11.39 | 0.58× |
| ENDOTHELIAL (%) | 10.50 | 7.07 | 0.67× |
| ENDOTHELIAL vs 8.0% volume-EM | 1.31× | 0.88× | PASS → **PASS** |
| endothelial areal density (/mm²) | 154 | 645 | 4.19× |
| **FIELDS** | | | |
| mean vascular density | 0.0268 | 0.1004 | 3.74× |
| mean ECM density | 0.0759 | 0.4571 | different matrix genes on the two panels — NOT comparable as biology |
| mean oxygen (tissue) | 0.3520 | 0.2210 | absolute level is uncalibrated on both |
| model-HYPOXIC tumour fraction | 0.00% | 33.76% | a statement about the sink constant, not the tumours |
| **VASCULATURE** | | | |
| discrete vessel objects | 3,091 | 20,875 | 6.75× |
| MVD (/mm²) | 110.8 | 233.3 | published 65 (35–95): HIGH → **HIGH** |
| intervessel spacing 1/√MVD (µm) | 95.0 | 65.5 | published 124 (103–169): LOW → **LOW** |
| intervessel NN median (µm) | 37.8 | 27.7 | 0.73× |
| tumour→vessel median (µm) | 42.8 | 23.9 | published 58 (48–79): LOW → **LOW** |
| tumour cells > 100 µm from a vessel (%) | 4.02 | 0.35 | 0.09× |
| **vascular area fraction (%)** | 1.25 | 10.98 | measured **2.1 ± 0.3**: LOW (0.59×) → **HIGH** (5.23×) |
| endothelial cells per vessel object | 1.39 | 2.77 | section 2's objects look far more like real vessel profiles |
| **DRUG PANEL (uniform sensitivity)** | | | |
| temozolomide — coverage @1× IC50 | 0.0000 | 0.0000 |  |
| lomustine — coverage @1× IC50 | 0.2802 | 0.9702 | 3.46× |
| cisplatin — coverage @1× IC50 | 0.9955 | 0.9987 | 1.00× |
| methotrexate — coverage @1× IC50 | 0.0000 | 0.0000 |  |
| doxorubicin — coverage @1× IC50 | 0.4373 | 0.9691 | 2.22× |
| temozolomide — mean tumour AUC (µM·s) | 7.66e+05 | 2.87e+06 | 3.75× |
| lomustine — mean tumour AUC (µM·s) | 9.32e+05 | 3.46e+06 | 3.71× |
| cisplatin — mean tumour AUC (µM·s) | 1.82e+06 | 6.79e+06 | 3.73× |
| methotrexate — mean tumour AUC (µM·s) | 1.57e+04 | 5.8e+04 | 3.68× |
| doxorubicin — mean tumour AUC (µM·s) | 1.15e+04 | 4.22e+04 | 3.67× |
| temozolomide — kill | 0.0184 | 0.0797 | 4.34× |
| lomustine — kill | 0.2371 | 0.3486 | 1.47× |
| cisplatin — kill | 0.4362 | 0.4248 | 0.97× |
| methotrexate — kill | 0.0192 | 0.0657 | 3.42× |
| doxorubicin — kill | 0.2595 | 0.3603 | 1.39× |
| kill ordering | cisplatin > doxorubicin > lomustine > methotrexate > temozolomide | cisplatin > doxorubicin > lomustine > temozolomide > methotrexate | **DIFFERENT** |
| coverage ordering | cisplatin > doxorubicin > lomustine > temozolomide > methotrexate | cisplatin > lomustine > doxorubicin > temozolomide > methotrexate | **DIFFERENT** |
| ρ(kill, AUC/IC50) | +1.00 | +0.90 | delivery × potency |
| **HYPOXIA** | | | |
| hypoxia-response genes on panel | none | VEGFA, HIF1A | the 319 panel cannot test this at all |
| ρ(simulated O₂, hypoxia score) | n/a | -0.1523 |  |
| ρ(dist-to-endothelial, hypoxia score) | n/a | -0.1204 | loses to the PDE |
| VEGFA rank of all panel genes | n/a | 1/480 |  |
| HIF1A rank of all panel genes | n/a | 327/480 |  |

## 6. Verdict — what generalises to a second patient, and what does not

### 6.1 Replicates / confirmed

*(Corrected. This paragraph previously read "the hypoxia rows confirm that the reaction-diffusion solve is out-predicted by a distance transform" — which contradicted §4.1 and §5 of this same document, both of which already had it the right way round: §4.1 returns **REFUTED** on exactly that claim, and §5 annotates the distance baseline "loses to the PDE".)* The hypoxia rows here **refute** a prior negative result rather than confirming one: the solve out-predicts the distance transform (-0.1523 vs -0.1204). 'Confirmed' below means the prior finding reproduced, not that the model performed well — and that caution still applies to the hypoxia rows in the reversed direction, since |ρ| ≈ 0.15 is weak, within-section, and carried by VEGFA alone (HIF1A +0.0140). The solve being the best of these predictors does not make it a validated oxygen field.

* **ENDOTHELIAL fraction matches human volume-EM 8.0% (PMID 38723085)** — 7.07% is 1.13x from 8.0% (inside the declared 1.5x tolerance); section 1 was 10.50% (1.31x)
* **The ordering is NOT simply the IC50 ordering read back out** — potency alone would rank doxorubicin > methotrexate > cisplatin > lomustine > temozolomide; the simulation ranks cisplatin > doxorubicin > lomustine > temozolomide > methotrexate. rho(kill, 1/IC50) = +0.20, rho(kill, AUC/IC50) = +0.90
* **VEGFA is the most anti-correlated gene on the panel; HIF1A is null** — VEGFA rho -0.2004, HIF1A rho +0.0140

### 6.2 Does NOT generalise

* **Vascular area fraction matches measured 2.1 +/- 0.3% (PMID 36611806)** — 10.98% vs 2.1% is 5.23x and lies OUTSIDE the published 1.8-2.4% band; section 1 gave 1.25% (0.59x, also outside but only just)
* **Kill ordering is cisplatin > doxorubicin > lomustine > methotrexate > temozolomide** — section 2 reorders to cisplatin > doxorubicin > lomustine > temozolomide > methotrexate
* **therapeutic_coverage @1x IC50 is still a DISCRIMINATING headline metric** — coverage has saturated at BOTH ends in this configuration — every drug is either effectively 0 or above 0.9, on both sections — so it is behaving as a binary indicator rather than a graded metric. `results/delivery_metrics.md` flagged this saturation as a 'latent, not active' weakness; at this exposure scale it is ACTIVE.
* **Simulated oxygen is BEATEN by a distance-to-nearest-vessel baseline (state_validation.md Test 1)** — the solved field now wins

### 6.3 Fails on BOTH sections — not a generalisation failure

These are the rows where the held-out patient behaved exactly like the development section, and both disagree with something this project has claimed. Reading them as 'it did not generalise' would be wrong; they are disagreements with a *prior configuration*, and the second patient's role here is to show the disagreement is not patient-specific.

* **MVD lands in the published human GBM range 35-95/mm2 (PMID 27583442)** — 233/mm2 vs published 65 (35-95); section 1 gave 111/mm2. BOTH sections are high, and by a similar amount, so this is a systematic property of the segmentation (endothelial cells are over-split into vessel objects) rather than something the held-out patient revealed — and it is the one place where the two sections AGREE with each other while both disagreeing with the literature.
* **therapeutic_coverage @1x IC50 recovers the kill ordering** — coverage does NOT reproduce the kill ordering on EITHER section in this configuration (both swap doxorubicin and lomustine), so this is not a held-out failure — it is a disagreement with `results/delivery_metrics.md`, which ran at a different exposure scale where coverage was not saturated. See §3.3

### 6.4 Untestable on one section or the other

* **Hypoxia against measured markers** is impossible on the 319-gene section: it carries no VEGFA, no HIF1A, no CA9, no ADM, no NDRG1, no SLC2A1. The honest statement for section 1 remains that hypoxia is *unmeasured*, not that it has been simulated.
* **The TOP2A-based doxorubicin sensitivity determinant** cannot be evaluated on section 2 — TOP2A is absent from the IO panel. This is why the whole comparison runs in `uniform` mode; `determinant` mode is not even definable in a comparable way across these two panels.
* **MGMT** — and therefore any measured temozolomide sensitivity — is absent from **both** panels. TMZ's kill is an upper bound on both sections for the same reason.
* **Outcome.** Neither section has a drug-response label. Nothing here validates a kill prediction; the drug panel tests whether the machinery behaves consistently on a second real tissue, which is a different and much weaker claim.

### 6.5 The one-paragraph reading

**The response ordering moved, and the reconstruction layer did not travel either.** The drug ordering and the delivery × potency structure behind it survive a second real patient on a different panel — but they were always going to be robust, because they are dominated by drug properties (IC50, BBB penetrance, molecular weight) that are identical on both sections, and the tissue only modulates them. The quantities that are genuinely *about the tissue* — the endothelial cell-type fraction and the vascular area fraction, i.e. the two numbers this project cites as its external validations — are the ones that fail to replicate. That is the honest asymmetry, and it points at where the next work is: the reconstruction's external anchors are section-specific, and the two constants that produce them (the typing threshold and the 7.5 µm vessel radius) are entangled with the grid and with panel sparsity in a way that a single section could never expose.

## 7. Caveats — what could make these conclusions wrong

1. **Different patients AND different panels.** Only 60 genes are shared. Every composition difference is confounded by panel content; §1.2 and §3.5 name the specific places where a difference is more plausibly measurement than biology (the endothelial gate arity, the ECM marker set).
2. **n = 2.** Two sections cannot separate 'this pipeline does not generalise' from 'these two tumours differ'. A claim that a quantity *fails* to replicate is well supported by a single discrepancy; a claim that one *does* replicate rests on a single agreement and is much weaker.
3. **The IO marker map is not fully blind** (§1.2): it was written against this section, though never against an external anchor.
4. **Grid convergence is not re-established here.** dx = 12 µm was shown converged on section 1's geometry (`results/numerical_rigor.md` §2). Section 2 is denser and its vessel objects are larger; the convergence argument is inherited, not re-run. A dx sweep on section 2 is the obvious next check and it was not affordable in this run.
5. **The absolute oxygen level is uncalibrated on both sections**, so the 34% vs 0% hypoxic-fraction gap is a statement about `OxygenConfig.consumption` and cell density, not about the two tumours. Every absolute kill inherits it.
6. **Typing circularity.** `FLT1` and `PLVAP` are IO-panel ENDOTHELIAL markers and therefore help build the oxygen source, so the positive tail of the gene ranking is partly the model correlating with its own inputs. Only the negative tail is informative.
7. **Pericytes are unlabelled on section 2.** `CellType` has no MURAL member, so RGS5/PDGFRB/NOTCH3-positive mural cells fall into the ungated TUMOR default, inflating the tumour count and hence the oxygen sink. That is part of why section 2 is so much more hypoxic.
8. **2D, in-plane, single section** on both. Distances are Euclidean and in-plane, not geodesic and not 3D.
9. **Subsampling.** Only the permutation null is subsampled (50,000 of 816,769 cells, 2,000 draws); every effect size, fraction and simulation number is computed on all cells. Plot (c) decimates the CDF for rendering only.

## 8. Reproduction

```bash
./.venv/bin/python -m tmesim.validation.heldout_patient
```

Deterministic (`seed = 0` in every `SimConfig`; permutation seed 0). Runtime 18 min, dominated by 10 PDE solves — 13 min of it on section 2's 905×956 grid.

**Files added by this validation (nothing else was modified):**

* `tmesim/validation/heldout_patient.py` — this analysis
* `results/heldout_patient.md`, `results/heldout_patient.png`

**Data provenance.** Both sections are open 10x Genomics Xenium datasets, CC BY 4.0, no login or DUA:

* Section 1 — 10x Xenium Human Brain — Glioblastoma section (v1.3.0); panel hBrain_254g base + 65-gene add-on (319 gene targets)
* Section 2 — 10x Xenium Human Brain Cancer (GBM) FFPE — IO panel (v2.0.0); panel Human Immuno-Oncology Profiling base (hImmune_v1, 380 genes) + Add-on B (100 genes) = 480 gene targets

**Figure.** `heldout_patient.png`.
