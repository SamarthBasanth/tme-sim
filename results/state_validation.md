# Do the simulated oxygen field and the oxygen-banded cell states carry real biology?

Two independent falsification tests against MEASURED markers in real GBM tissue, on two different open 10x Xenium sections.

| | Test 1 — simulated oxygen field | Test 2 — cell-state labelling |
|---|---|---|
| **Verdict** | **PARTIALLY SUPPORTED** | **UNDERPOWERED** |
| Section | 10x Xenium Human Brain Cancer (GBM) FFPE — IO panel (v2.0.0) | 10x Xenium Human Brain — Glioblastoma section (v1.3.0) |
| Cells | 816,769 (480-gene panel) | 40,887 (30,174 tumour; 319-gene panel) |
| Measured against | VEGFA + HIF1A (hypoxia response) | MKI67 + TOP2A (proliferation) |
| Headline statistic | Spearman ρ = -0.1132 (p = 0.00050) | MKI67+ in VIABLE 2.22% (4/180) |
| Killer detail | a **one-line distance-to-vessel baseline scores -0.1503**, i.e. 1.33× the solved field | the test could only have detected a **≥4.1× rate ratio** |

**In one sentence each.** The simulated oxygen field *is* significantly anti-correlated with measured hypoxia response, but it is out-predicted by a physics-free distance-to-nearest-vessel baseline computed from its own inputs — so the reaction–diffusion solve is not what carries the signal. And the proliferation test on the 319-gene section is dropout-limited to the point where it cannot distinguish a correct cell-state model from noise; that is a power failure, not a vindication and not a refutation.

> **The two sections are DIFFERENT PATIENTS.** Test 1 and Test 2 are two experiments, not two views of one tumour. Nothing here transfers from one to the other.

## 0. What is under test, and why it matters

The chain is: `physics.solve_oxygen` produces an oxygen field → `cells.agents.update_cell_states` bands it into HYPOXIC / QUIESCENT / VIABLE at `t = OxygenConfig.hypoxia_threshold` = 0.2 → `cells.pkpd` multiplies drug kill by `PDConfig.hypoxia_resistance` = 0.35 and `PDConfig.quiescent_resistance` = 0.5.

`results/dose_schedule.md` shows the consequence: for temozolomide — whose per-cell sensitivity is uniform 1.0 because MGMT is not on the panel — the kill CEILING is *entirely* those two multipliers, and ~89% of it comes from the quiescent band alone. If the oxygen field or the banding is wrong, every absolute kill number in this repository inherits the error.

Both tests drive the **frozen production path** (`reconstruct.real_pipeline.reconstruct_real` → `physics.solve_oxygen` → `cells.agents.update_cell_states`). Nothing is re-implemented here.

### 0.1 Regression guard — the IO panel map did not disturb the 319-gene section

Test 1 needs a marker map for the 480-gene Immuno-Oncology panel (`reconstruct.marker_panels`), because the brain-panel map collapses that section to TUMOR 93.1% / IMMUNE 6.9% with **no ENDOTHELIAL cells at all** — and the endothelial cells are the oxygen SOURCE, so that failure would have silently produced a supply-free oxygen field. Adding it is only admissible if it is inert on the other section:

| class | this run | historical reference | drift |
|---|---|---|---|
| TUMOR | 73.80% | 73.80% | 0.001 pp |
| STROMAL | 7.92% | 7.92% | 0.002 pp |
| IMMUNE | 8.54% | 8.54% | 0.004 pp |
| ENDOTHELIAL | 9.74% | 9.74% | 0.004 pp |

`select_marker_map` returns the **same dict object** as `BRAIN_PANEL_MARKER_MAP` for the 319-gene panel (`is` identity: True), so the historical typing is bit-for-bit unchanged. Max drift 0.004 pp. **PASS**

## 1. TEST 1 — simulated oxygen vs measured hypoxia markers

**Section.** 10x Xenium Human Brain Cancer (GBM) FFPE — IO panel (v2.0.0)  
**Panel.** Human Immuno-Oncology Profiling base (hImmune_v1, 380 genes) + Add-on B (100 genes) = 480 gene targets  
**Panel selection.** `io_480` (marker coverage: {'brain_319': 0.5, 'io_480': 1.0})  
**Cells.** 816,769 × 480 gene targets

Measured hypoxia-response genes on this panel:

| gene | cells with ≥1 transcript |
|---|---|
| VEGFA | 63.36% |
| HIF1A | 54.93% |

The two genes correlate with **each other** at Spearman ρ = +0.0393. That is close to zero: VEGFA and HIF1A are *not* two readouts of one latent 'hypoxia' variable in these data, so the composite score is reported alongside — never instead of — the per-gene results.

### 1.1 The interrupted run's claim, adjudicated

> *"VEGFA is the single most anti-correlated gene of all 480 with simulated oxygen — but HIF1A is exactly null."*

**Both halves CONFIRMED.** Every one of the 480 panel genes was depth-normalised and rank-correlated against the simulated oxygen at each cell (dx = 12 µm, all 816,769 cells):

| rank (1 = most negative) | gene | ρ vs simulated O₂ | detection |
|---|---|---|---|
| 1 | VEGFA **←** | -0.1479 | 63.36% |
| 2 | SOX9 | -0.1227 | 79.41% |
| 3 | PTN | -0.1197 | 94.57% |
| 4 | NCAM1 | -0.1114 | 83.24% |
| 5 | IL1RAP | -0.1040 | 45.74% |
| 6 | S100B | -0.0992 | 96.95% |
| 7 | FEZ1 | -0.0874 | 83.21% |
| 8 | CCND1 | -0.0836 | 77.43% |
| … | … | … | … |
| 249 | HIF1A **←** | +0.0085 | 54.93% |
| … | … | … | … |
| 480 | S100A9 | +0.1869 | 7.67% |
| 479 | A2M | +0.1668 | 51.44% |
| 478 | SPARC | +0.1487 | 78.80% |
| 477 | RGS5 | +0.1430 | 12.45% |
| 476 | IFITM3 | +0.1262 | 66.57% |

* **VEGFA is rank 1 of 480** at ρ = -0.14793. Of every gene the panel measures, the canonical hypoxia-response gene is the one the simulated oxygen field most strongly anti-correlates with. That is a real, non-trivial hit and it is the strongest positive evidence in this report.
* **HIF1A is rank 249 of 480** at ρ = +0.00850 — null, and if anything the *wrong sign*. This is expected on mechanism, not a surprise: HIF-1α is regulated by oxygen-dependent prolyl hydroxylation and VHL-mediated proteolysis of the **protein** (Ivan et al., Science 2001;292(5516):464-468 (PMID 11292862, doi:10.1126/science.1059817); Jaakkola et al., Science 2001;292(5516):468-472 (PMID 11292861, doi:10.1126/science.1059796)), so its mRNA has no obligation to track oxygen. HIF1A mRNA should not have been treated as a hypoxia readout in the first place, and the composite score is diluted by including it.

**The circularity, stated plainly.** The positive end of that table is the vascular and myeloid programme — RGS5, FLT1, PDGFRB, A2M, IGFBP7, CD163. `FLT1` is one of the three ENDOTHELIAL markers that *builds* the oxygen source term, so its positive correlation is partly the model correlating with itself. Only the NEGATIVE end of the ranking is informative, which is why VEGFA's rank-1 position is the claim worth making.

VEGFA is the right measured readout here: Fischer et al., Brain Pathol 2005;15(4):297-310 (PMID 16389942, doi:10.1111/j.1750-3639.2005.tb00115.x) — *in GBM, vascular involution produces necrosis and hypoxia, which induces HIF-1alpha and VEGF expression specifically in the perinecrotic pseudopalisading glioma cells*.

### 1.2 (a) Is simulated oxygen anti-correlated with measured hypoxia — and does it beat a trivial baseline?

Per-cell Spearman ρ against the measured hypoxia score (mean z of depth-normalised VEGFA and HIF1A), computed on **all 816,769 cells** at dx = 12 µm. Every predictor is sign-aligned to oxygen (larger = better perfused), so a correct predictor gives a NEGATIVE ρ.

| predictor | what it is | ρ vs hypoxia score | ρ vs VEGFA | ρ vs HIF1A | ρ with log-depth partialled out |
|---|---|---|---|---|---|
| **SIMULATED O₂ (the PDE)** | the reaction-diffusion solve under test | **-0.1132** | -0.1479 | +0.0085 | -0.1283 |
| vascular_kde | the solver's own input field, undiffused | **-0.1270** | -0.1456 | -0.0085 | -0.1320 |
| dist_endothelial | one cKDTree query on the PDE's own source cells; no physics | **-0.1503** | -0.1647 | -0.0290 | -0.1545 |
| local_cell_density | pure local crowding; ignores vessels entirely | **-0.0192** | -0.0031 | -0.0439 | -0.0336 |

Permutation p-values (label shuffling, 2,000 draws on a seeded random subsample of 50,000 cells — subsampled because the kernel is O(n × draws); every effect size above is on the full population, so subsampling can only make a p-value *more* conservative):

| predictor | ρ on the subsample | permutation p |
|---|---|---|
| simulated_o2 | -0.1090 | 0.00050 |
| dist_endothelial | -0.1549 | 0.00050 |
| vascular_kde | -0.1281 | 0.00050 |

**The answer to (a) is yes — and it does not matter.** The simulated field is genuinely, significantly anti-correlated with measured hypoxia response (ρ = -0.1132, p = 0.00050, floor of the estimator). But so is every trivial baseline built from the same inputs, and two of them do it BETTER:

* **`dist_endothelial` (-0.1503) beats the PDE (-0.1132) by 33%.** One `cKDTree` query — the distance from each cell to the nearest endothelial cell — predicts measured VEGFA/HIF1A better than the reaction–diffusion solve does.
* **`vascular_kde` (-0.1270) also beats the PDE.** That is the solver's OWN SOURCE TERM, sampled before it is diffused. The solve takes a field correlating at -0.1270 and returns one correlating at -0.1132: the diffusion step *destroys* predictive information rather than adding it.
* **`local_cell_density` (-0.0192) is near zero**, which is the sanity floor — not every spatial variable correlates with the marker, so the effects above are real signal, not autocorrelation artefacts.

**Nested test — does the PDE add anything beyond the raw geometry?**

| partial Spearman | ρ |
|---|---|
| hypoxia score vs **simulated O₂**, distance-to-vessel held fixed | -0.0660 |
| hypoxia score vs **−distance-to-vessel**, simulated O₂ held fixed | -0.1192 |
| simulated O₂ vs −distance-to-vessel (how nested the two predictors are) | +0.3459 |

Distance-to-vessel retains 0.1192 of explanatory power after the simulated field is conditioned out; the simulated field retains only 0.0660 after distance is conditioned out. **The PDE is the redundant variable in that pair.**

**Depth control.** The simulated field correlates with total transcript count — a purely technical variable — at ρ = +0.0975. Partialling log-depth out of the table above changes nothing material (simulated O₂ goes -0.1132 → -0.1283; distance goes -0.1503 → -0.1545), so the ordering is not a depth artefact.

### 1.3 (b) Are model-HYPOXIC cells enriched for measured VEGFA/HIF1A vs model-VIABLE?

State occupancy on this section at dx = 12 µm (all cells): **HYPOXIC** 73.68%, **QUIESCENT** 26.18%, **VIABLE** 0.13%.

| marker | contrast | n (A vs B) | mean A | mean B | Cohen's d | Cliff's δ | permutation p |
|---|---|---|---|---|---|---|---|
| hypoxia_score | hypoxic vs viable | 601,832 vs 1,090 | +0.0474 | -0.5724 | **+0.8621** | +0.5129 | 0.00050 |
| hypoxia_score | hypoxic vs quiescent | 601,832 vs 213,847 | +0.0474 | -0.1305 | **+0.2514** | +0.1439 | 0.00050 |
| hypoxia_score | hypoxic vs viable TUMOUR | 491,814 vs 411 | +0.0669 | -0.4707 | **+0.7280** | +0.4244 | 0.00050 |
| VEGFA | hypoxic vs viable | 601,832 vs 1,090 | +0.0906 | -0.7993 | **+0.8589** | +0.5320 | 0.00050 |
| VEGFA | hypoxic vs quiescent | 601,832 vs 213,847 | +0.0906 | -0.2508 | **+0.3452** | +0.1820 | 0.00050 |
| VEGFA | hypoxic vs viable TUMOUR | 491,814 vs 411 | +0.1380 | -0.7221 | **+0.8039** | +0.4843 | 0.00050 |
| HIF1A | hypoxic vs viable | 601,832 vs 1,090 | +0.0042 | -0.3455 | **+0.3468** | +0.1892 | 0.02549 |
| HIF1A | hypoxic vs quiescent | 601,832 vs 213,847 | +0.0042 | -0.0102 | **+0.0144** | +0.0051 | 0.01299 |
| HIF1A | hypoxic vs viable TUMOUR | 491,814 vs 411 | -0.0042 | -0.2194 | **+0.2063** | +0.1180 | 0.07146 |

Depth-stratified control (labels permuted **within** total-transcript deciles, so only comparably-deep cells are ever compared):

| marker | HYPOXIC vs VIABLE, depth-stratified difference | p |
|---|---|---|
| hypoxia_score | +0.57470 | 0.00050 |
| VEGFA | +0.84227 | 0.00050 |
| HIF1A | +0.30714 | 0.00950 |

**The answer to (b) is yes for VEGFA, no for HIF1A.** Model-HYPOXIC cells carry substantially more measured VEGFA than model-VIABLE cells (d = +0.859, Cliff's δ = +0.532, p = 0.00050), and the effect survives depth stratification. For HIF1A the same contrast is only d = +0.347, and the tumour-only version is not significant at all (p = 0.07146).

**Two things blunt this result.** First, the VIABLE group holds only 1,090 of 816,769 cells (0.133%), so "HYPOXIC vs VIABLE" is a comparison against a sliver. The better-balanced HYPOXIC-vs-QUIESCENT contrast is far weaker (VEGFA d = +0.345, HIF1A d = +0.014). Second — and decisively — §1.2 already showed the same separation is obtained, more strongly, by thresholding distance-to-vessel. A large Cohen's d here is evidence that the *bands* sit somewhere sensible in space, not that the PDE computed them.

### 1.4 (c) Field-vs-field agreement on the grid

The per-cell test is dominated by single-cell dropout. Averaging both fields over fixed 60 µm tissue tiles (≥5 measured cells per tile) suppresses that noise and asks the coarser question: does a REGION the model calls hypoxic express more hypoxia response? Tiles are a fixed physical size at every dx, so the three rows are comparable to each other rather than to the mesh.

| dx (µm) | grid | tiles used | Spearman ρ | Pearson r |
|---|---|---|---|---|
| 30 | 363×383 | 24,254 | -0.2741 | -0.3188 |
| 20 | 543×574 | 24,254 | -0.2740 | -0.3187 |
| 12 | 905×956 | 24,254 | -0.2736 | -0.3184 |

Agreement is ~2.4× stronger at the field level (ρ ≈ -0.274 vs -0.113 per cell), which is the expected signature of per-cell dropout rather than of a better model. It is the most favourable framing of the test and it is reported as such.

### 1.5 (d) Does agreement improve as the grid is refined?

**No — it is flat to four decimal places.**

| dx (µm) | grid | mean O₂ (tissue) | HYPOXIC / QUIESCENT / VIABLE (all cells) | ρ(sim O₂, hypoxia score) | ρ(dist-to-vessel, hypoxia score) |
|---|---|---|---|---|---|
| 30 | 363×383 | 0.1757 | 73.81% / 26.06% / 0.13% | -0.1132 | -0.1503 |
| 20 | 543×574 | 0.1772 | 72.71% / 27.14% / 0.15% | -0.1133 | -0.1503 |
| 12 | 905×956 | 0.1760 | 73.68% / 26.18% / 0.13% | -0.1132 | -0.1503 |

This is the *good* news about the numerics and the *bad* news about the model. With the geometric tissue mask and the areal oxygen sink both in place, the field is genuinely grid-converged: mean O₂ moves by 0.2% and the correlation by <0.0002 over a 2.5× range of dx. So the disagreement with measured biology is **not** a resolution problem — refining the mesh cannot rescue it. The model has converged onto an answer that a distance transform beats.

### 1.6 Verdict — Test 1

## **PARTIALLY SUPPORTED**

* the sign and significance are right (rho = -0.1132, p = 0.00050) but a trivial physics-free baseline (dist_endothelial, rho = -0.1503) matches or BEATS the solved field, so the reaction-diffusion solve is not what is carrying the signal
* the composite is carried by ONE gene: VEGFA rho=-0.1479, HIF1A rho=+0.0085 — HIF1A is null

**Consequences, stated plainly.**

1. **The oxygen field is not falsified — it is out-competed.** Its sign is right, its significance is overwhelming, and VEGFA being its rank-1 anti-correlate out of 480 genes is a genuine hit. But the reaction–diffusion solve is not what produces that hit: its own undiffused source term and a plain distance transform both do better. Whatever the PDE is adding, it is not predictive accuracy about hypoxia.
2. **Hypoxia-driven resistance in the kill model is therefore not supported *as physics*.** The `PDConfig.hypoxia_resistance` multiplier is applied to cells selected by a field that a one-line geometric heuristic selects better. Any kill number that leans on the *spatial pattern* of hypoxia is, at best, an expensive re-derivation of "far from a vessel".
3. **Hypoxia should be derived from measured markers where they exist.** On a section that measures VEGFA, the defensible per-cell hypoxia covariate is the measured VEGFA level, not the simulated oxygen — and where no such marker is on the panel (the 319-gene section: no VEGFA, no CA9, no ADM, no NDRG1, no SLC2A1), the honest statement is that hypoxia is **unmeasured**, not that it has been simulated.
4. **Drop HIF1A from any hypoxia score.** Its mRNA is null against the field (rank 249/480) for a known mechanistic reason, and including it dilutes the composite from ρ ≈ -0.148 (VEGFA alone) to -0.113.
5. **The absolute oxygen level remains uncalibrated.** Nothing here tests it. `OxygenConfig.consumption` has never been fitted to a measured oxygen tension, and on this section the model calls 73.7% of all cells HYPOXIC while on the 319-gene section under the shipping mask it calls 0.0% hypoxic (§2.2). Those two numbers cannot both be right.

## 2. TEST 2 — cell-state model vs measured proliferation

**Section.** 10x Xenium Human Brain — Glioblastoma section (v1.3.0)  
**Panel.** hBrain_254g base + 65-gene add-on (319 gene targets)  
**Cells.** 40,887 total, 30,174 tumour  
**Headline configuration.** tissue mask `occupancy`, dx = 12 µm — *the configuration `results/dose_schedule.md` actually ran under*, which is the one whose kill ceilings are at stake.

### 2.1 The measured proliferating fraction, under two documented rules

| marker | cells at exactly 0 | rule `nonzero` (count ≥ 1) | rule `percentile` (top 20% of depth-normalised expression) |
|---|---|---|---|
| MKI67 | 99.56% | 0.444% (134/30,174) | 0.444% — **DEGENERATE** |
| TOP2A | 98.94% | 1.057% (319/30,174) | 1.057% — **DEGENERATE** |

**The percentile rule is mathematically degenerate on this panel, and that IS the result.** The rule was written to bracket the published Ki-67 index — "call the top 20% of cells proliferating" — but with >98.9% of tumour cells at exactly zero counts, the 80th percentile of expression *is* zero. The threshold falls inside a tie mass at the origin and the rule collapses back onto `count ≥ 1`. There is no percentile of this data that corresponds to a Ki-67 index, because the data does not reach that far.

**Against the literature.** Published Ki-67 labelling indices in glioblastoma (PubMed, all verified):

| cohort | Ki-67 index | source |
|---|---|---|
| GBM, n=35 (prognostic cutoff) | 19% | Fudaba et al., World Neurosurg 2021;149:e469-e480 (PMID 33567368, doi:10.1016/j.wneu.2021.02.002) |
| GBM, n=44 (cohort split point) | 27% | Alkhaibary et al., Hematol Oncol Stem Cell Ther 2018;12(2):82-88 (PMID 30552865, doi:10.1016/j.hemonc.2018.11.001) |
| epithelioid GBM, n=16 (median, range 10-50%) | 30% | Zheng et al., Zhonghua Bing Li Xue Za Zhi 2021;50(3):229-235 (PMID 33677887, doi:10.3760/cma.j.cn112151-20200617-00479) |
| GBM, n=70 (median) | 40% | Caramanti et al., Cureus 2022;14(4):e24246 (PMID 35602791, doi:10.7759/cureus.24246) |

The panel detects MKI67 in **0.444%** of tumour cells — **43× to 90× below** the published range. This is not a claim that this tumour has a Ki-67 index of 0.44%; it is a statement about probe sensitivity and transcript dropout on a targeted in-situ panel. The reference standard is itself imprecise (Nielsen et al., Diagn Pathol 2018;13:38 (PMID 29885671, doi:10.1186/s13000-018-0711-2) — Ki-67 LI scoring in gliomas agrees poorly between observers: kappa 0.32 for whole-slide and as low as 0.04 for hot-spot scoring (n=235 gliomas)), but not by two orders of magnitude.

**Model VIABLE fraction, for comparison.** 0.60% in the headline configuration, but 0.60% – 98.19% across the mask × dx grid (§2.2). **Do not read any agreement between "model VIABLE %" and "MKI67-detected %" as corroboration** — it would be a coincidence between two unrelated small numbers, one of which is a bookkeeping artefact and the other a dropout rate.

### 2.2 The model's state vector is set by bookkeeping, not by tissue

| tissue mask | dx (µm) | mean O₂ (tissue) | HYPOXIC | QUIESCENT | VIABLE |
|---|---|---|---|---|---|
| `geometric` | 30 | 0.5868 | 0.00% | 1.81% | 98.19% |
| `geometric` | 20 | 0.5819 | 0.00% | 3.81% | 96.19% |
| `geometric` | 12 | 0.5755 | 0.00% | 5.08% | 94.92% |
| `occupancy` | 30 | 0.5259 | 0.00% | 6.71% | 93.29% |
| `occupancy` | 20 | 0.4248 | 0.00% | 25.97% | 74.03% |
| `occupancy` ← | 12 | 0.2498 | 13.43% | 85.98% | 0.60% |

Read that table before reading any statistic below it. The tumour-cell VIABLE fraction runs from 0.60% to 98.19% — a 165× swing — across two tissue-mask conventions and three grid spacings. **Neither the mask nor the mesh is a measurement.** Under the current shipping default (`geometric`) the HYPOXIC state is *empty at every dx*, so `PDConfig.hypoxia_resistance` has no cells to act on at all; the 86% QUIESCENT / 13.4% HYPOXIC split that `results/dose_schedule.md` quotes comes from the legacy `occupancy` mask, whose coverage shrinks as the mesh is refined. Those kill ceilings are conditional on a convention that the reconstruction layer has since abandoned.

### 2.3 Are model-VIABLE cells more proliferative? — and could we have told?

| marker | state | rate | k / n | Wilson 95% CI |
|---|---|---|---|---|
| MKI67 | HYPOXIC | 0.592% | 24 / 4,051 | 0.398–0.880% |
| MKI67 | QUIESCENT | 0.409% | 106 / 25,943 | 0.338–0.494% |
| MKI67 | VIABLE | 2.222% | 4 / 180 | 0.867–5.574% |
| TOP2A | HYPOXIC | 1.358% | 55 / 4,051 | 1.045–1.763% |
| TOP2A | QUIESCENT | 1.010% | 262 / 25,943 | 0.895–1.139% |
| TOP2A | VIABLE | 1.111% | 2 / 180 | 0.305–3.960% |

| marker | contrast | Cohen's d | Cliff's δ | permutation p | depth-stratified p | MINIMUM DETECTABLE rate ratio | expected positives |
|---|---|---|---|---|---|---|---|
| MKI67 | VIABLE vs rest | +0.2712 | +0.0179 | 0.0078 | 0.0069 | **4.14×** | 0.8 |
| MKI67 | HYPOXIC vs QUIESCENT | +0.0500 | +0.0018 | 0.1187 | 0.0148 | **1.71×** | 18.0 |
| TOP2A | VIABLE vs rest | +0.0687 | +0.0006 | 1.0000 | 0.2781 | **3.03×** | 1.9 |
| TOP2A | HYPOXIC vs QUIESCENT | +0.0571 | +0.0035 | 0.0478 | 0.4867 | **1.46×** | 42.8 |

**The power analysis is the result.** The "minimum detectable rate ratio" column is the smallest enrichment the contrast could have found at α = 0.05, power = 0.8 given the observed positivity rate and the group sizes.

* For **MKI67, VIABLE vs rest**, the VIABLE group holds 180 tumour cells and would contain 0.8 MKI67-positive cells under the null. The test could only have detected a **≥4.1× rate ratio**.
* The measured biology this is supposed to detect is nowhere near that size: Brat et al., Cancer Res 2004;64(3):920-927 (PMID 14871821, doi:10.1158/0008-5472.CAN-03-2073) measured that hypoxic pseudopalisading GBM cells are *5–50% LESS* proliferative than adjacent astrocytoma — a rate ratio of roughly 0.5–0.95, i.e. a **depletion of at most 2×**. A test that needs 4.1× to fire cannot see an effect of that size in either direction.
* Therefore **no conclusion about the VIABLE label can be drawn from this section.** Not "the state model is wrong"; not "the state model is right". The measurement is not capable of adjudicating it.

**A positive result is visible, and it must not be believed.** The point estimate for MKI67 goes the RIGHT way — 2.22% positive in VIABLE against 0.444% section-wide, a 5.0× rate ratio, permutation p = 0.0078, surviving depth stratification (p = 0.0069). It rests on **4 positive cells out of 180**. Its Wilson interval spans 0.87–5.57%, one Poisson fluctuation moves it from 4 to 2 or 7, and TOP2A — the other proliferation marker on the same cells — shows nothing (1.11%, 2/180, p = 1.0000). This is precisely the regime the minimum-detectable-effect calculation says is uninterpretable, and reporting it as support would be the exact error this module exists to prevent.

**What the test CAN rule out.** The better-powered HYPOXIC-vs-QUIESCENT contrast (n = 4,051 vs 25,943) needs only 1.71× for MKI67 and 1.46× for TOP2A, which is inside the plausible range. It rules out a **large** proliferation deficit in the model's HYPOXIC band — and it does so in the wrong direction:

* **MKI67+ is HIGHER in model-HYPOXIC (0.592%, 24/4,051) than in model-QUIESCENT (0.409%, 106/25,943)** — the opposite of the measured biology. Permutation p = 0.1187; depth-stratified p = 0.0148.
* **TOP2A+ is HIGHER in model-HYPOXIC (1.358%, 55/4,051) than in model-QUIESCENT (1.010%, 262/25,943)** — the opposite of the measured biology. Permutation p = 0.0478; depth-stratified p = 0.4867.

Neither reversal is large, and both sit close to the depth-stratification boundary, so this is a weak contradiction rather than a strong one. But it is a contradiction, and it points the same way as the continuous test:

| continuous test | Spearman ρ | p | ρ with log-depth removed |
|---|---|---|---|
| simulated O₂ vs MKI67 (per tumour cell) | +0.0078 | 0.1812 | +0.0075 |
| simulated O₂ vs TOP2A (per tumour cell) | +0.0037 | 0.5284 | +0.0033 |

Rank-correlating the **simulated oxygen** against measured proliferation across all 30,174 tumour cells — the best-powered version of the question, because it uses every cell and does not depend on where the band edges fall — returns essentially nothing, with a *positive* sign where a negative one was not expected either. The oxygen field carries no per-cell proliferation information on this section.

### 2.4 dx-dependence of the Test 2 contrast

| mask | dx (µm) | VIABLE fraction | MKI67+ in VIABLE | MKI67+ elsewhere | min. detectable rate ratio |
|---|---|---|---|---|---|
| `geometric` | 30 | 98.19% | 0.446% | 0.366% | 2.8× |
| `geometric` | 20 | 96.19% | 0.451% | 0.261% | 2.3× |
| `geometric` | 12 | 94.92% | 0.454% | 0.261% | 2.1× |
| `occupancy` | 30 | 93.29% | 0.455% | 0.296% | 2.0× |
| `occupancy` | 20 | 74.03% | 0.439% | 0.459% | 1.6× |
| `occupancy` | 12 | 0.60% | 2.222% | 0.433% | 4.1× |

The measured MKI67+ rate inside the model's VIABLE band is ~0.44–0.46% at five of the six configurations — i.e. **indistinguishable from the section-wide rate of 0.444%** — because at those settings VIABLE contains 74–98% of the tumour and the "contrast" is a population against itself. Only in the sixth configuration (the dose_schedule one) does VIABLE become small enough to be a real subgroup, and there it is too small to test. There is no dx at which this comparison is both meaningful and powered.

### 2.5 Verdict — Test 2

## **UNDERPOWERED**

* MKI67 VIABLE-vs-rest could only have detected a 4.1x rate ratio (0.8 positives expected in VIABLE under the null)
* TOP2A VIABLE-vs-rest could only have detected a 3.0x rate ratio (1.9 positives expected in VIABLE under the null)
* the decisive VIABLE-vs-rest comparison cannot distinguish 'the state labelling is right' from 'the state labelling is noise'
* independently of power, the model's VIABLE fraction is not a measurement: it runs 0.60% to 98.19% across two tissue-mask conventions and three grid spacings, all of which are bookkeeping choices

**Consequences, stated plainly.**

1. **The quiescent labelling is UNSUPPORTED — in the specific sense that it has not been, and on this section cannot be, tested.** The measurement that would adjudicate it is 40–90× too insensitive. Nobody should read this as clearance.
2. **The kill ceilings inherit that.** `results/dose_schedule.md` derives ~89% of the temozolomide kill ceiling from the quiescent band, via a ×0.5 multiplier applied to a label with no measured support and a 165× convention-dependent size. Every absolute kill fraction and every ceiling in that report carries an unquantified error of that magnitude. The dose/schedule **comparisons** are far more robust — they share the same multiplier and it largely cancels — but the **levels** are not defensible as measurements.
3. **The one thing the data does say points the wrong way.** In the better-powered HYPOXIC-vs-QUIESCENT contrast, both measured proliferation markers are *higher* in the model's HYPOXIC band. That is the opposite of the measured tissue biology (Brat et al., Cancer Res 2004;64(3):920-927, PMID 14871821).
4. **What would fix it.** Either (a) a proliferation readout with a real detection rate — IHC/IF Ki-67 on a serial section, or a whole-transcriptome platform — or (b) derive the state labels from measured markers (a per-cell proliferation score) instead of from the simulated oxygen, and report the states as *observed* rather than *predicted*.

## 3. Caveats — what could make these conclusions wrong

1. **The two sections are different patients.** Test 1 runs on `data_cache/xenium_gbm_io` (816,769 cells, IO panel); Test 2 on `data_cache/xenium_gbm` (40,887 cells, brain panel). They are different tumours with different vascular architectures and different panels. Nothing here is a paired observation, and the two verdicts do not reinforce each other.
2. **VEGFA and HIF1A are hypoxia-RESPONSE genes, not an oxygen tension.** They are the transcriptional output of a pathway with its own kinetics and its own non-hypoxic inducers. A simulated field could be a perfectly good oxygen tension and still correlate weakly with them. What the trivial-baseline control establishes is narrower and safer: *whatever* VEGFA is reporting, a distance transform predicts it better than the PDE does — and that comparison is internal, so it is immune to this caveat.
3. **Targeted-panel dropout.** A zero is "no probe fired", not "no transcript". This is fatal for Test 2 (MKI67 at 0.44%) and merely limiting for Test 1 (VEGFA at 63.4%, which is why Test 1 is powered and Test 2 is not).
4. **Typing circularity.** `FLT1` and `PLVAP` are IO-panel ENDOTHELIAL markers, so they help build the oxygen source; the positive tail of the 480-gene ranking is partly the model correlating with its own inputs. `MKI67` is one of the four brain-panel TUMOR markers, so an MKI67-high cell is marginally more likely to be typed TUMOR — that inflates the MKI67+ rate *uniformly across states*, so it cannot manufacture a between-state contrast, but it does inflate the marginal rate slightly.
5. **`SLC1A3`/`S100B`/`PTN` (the IO STROMAL class) are expressed by malignant astrocyte-like glioma cells too**, so the IO section's STROMAL compartment is not clean. The brain panel's `AQP4`/`GJA1` have the same ambiguity, so the two sections are consistently rather than differently wrong here.
6. **Pericytes are unlabelled on the IO section.** The panel measures `RGS5`/`PDGFRB`/`NOTCH3` but `CellType` has no MURAL member, so mural cells fall into the ungated TUMOR default. That inflates the tumour count and hence the oxygen sink.
7. **The absolute oxygen level is uncalibrated.** `OxygenConfig.consumption` has never been fitted against a measured oxygen tension. The two sections disagree by orders of magnitude in hypoxic fraction (73.7% of IO cells vs 0.0% of brain-panel cells under the shipping mask), which is a statement about the sink constant and the cell density, not about the two tumours.
8. **Subsampling.** Permutation p-values on the IO section use a seeded random subsample of 50,000 of 816,769 cells with 2,000 draws. All effect sizes and all Spearman ρ are on the full population. Subsampling the null can only make a p-value more conservative, and every p reported here is already at the estimator's resolution floor.

## 4. Reproduction

```bash
python -m tmesim.validation.state_validation
```

Deterministic: the only randomness is the seeded permutation / subsampling (`PERM_SEED = 0`). Runtime ~8 minutes, dominated by six reconstructions of the 816,769-cell IO section and the 480-gene correlation scan.

**Files added by this validation (nothing else was modified):**

* `tmesim/validation/state_validation.py` — this analysis
* `tmesim/reconstruct/marker_panels.py` — the panel-aware marker map (regression-guarded in §0.1)
* `results/state_validation.md`, `results/state_validation.png`

**Data provenance.** Both sections are open 10x Genomics Xenium datasets, CC BY 4.0, no login or DUA:

* Test 1 — 10x Xenium Human Brain Cancer (GBM) FFPE — IO panel (v2.0.0); panel: Human Immuno-Oncology Profiling base (hImmune_v1, 380 genes) + Add-on B (100 genes) = 480 gene targets
* Test 2 — 10x Xenium Human Brain — Glioblastoma section (v1.3.0); panel: hBrain_254g base + 65-gene add-on (319 gene targets)

Neither section carries a paired drug-response outcome label, so nothing here validates a kill prediction — only the internal machinery that produces one.

**Figure.** `results/state_validation.png`.
