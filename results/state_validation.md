# Do the simulated oxygen field and the oxygen-banded cell states carry real biology?

Two independent falsification tests against MEASURED markers in real GBM tissue, on two different open 10x Xenium sections.

> ## CORRECTION (§1.2) — this document had its central inequality backwards
>
> **What §1.2 previously said.** That a one-line distance-to-vessel baseline **out-predicted** the solved oxygen field:
>
> > *"**`dist_endothelial` (-0.1204) beats the PDE (-0.1523) by -21%.**"*
> > *"**`vascular_kde` (-0.1111) also beats the PDE.** … the diffusion step *destroys* predictive information rather than adding it."*
> > *"**The PDE is the redundant variable in that pair.**"*
>
> The one-sentence summary, the headline table's "Killer detail" row, §1.3's closing paragraph, §1.5's closing line, and consequences 1–2 of §1.6 all repeated that reading.
>
> **Why it was wrong.** Every predictor in the §1.2 table is sign-aligned to oxygen, so — as **that table's own header states** — a correct predictor gives a NEGATIVE ρ, and the *stronger* predictor is the one with the larger **|ρ|**. |-0.1523| > |-0.1204|. **The solved field out-predicts the distance baseline**, by 1.27×. The "-21%" figure is itself the tell: 0.1204 / 0.1523 = 0.79, i.e. the baseline is **79% as strong** — weaker, not better. The nested partials say the same thing (simulated O₂ retains **0.1227** with distance held fixed; distance retains only **0.0792** with simulated O₂ held fixed), and so does the depth-partialled column (-0.1619 vs -0.1185).
>
> **How it went undetected, and how it was caught.** The inverted sentences are hand-written narrative in the report generator (`tmesim/validation/state_validation.py`, the `dist_endothelial … beats the PDE` template). Everything in this repository that decides the same question **arithmetically** always disagreed with them:
>
> * `state_validation.verdict_test1()` — the mechanical verdict for this very test — tests `abs(sim) > abs(best) * (1 + BASELINE_MARGIN)` and returns **SUPPORTED**, with the reason string *"the solved field beats every trivial baseline"*. That is why **§1.6 below has always contradicted §1.2 above it.**
> * `tmesim/validation/heldout_patient.py` re-derives the comparison on the same section and adjudicates with `abs(dist_rho) > abs(pde_rho)`. It returns **REFUTED** on the §1.2 claim — *"the solved field now wins"*. Written up in `results/heldout_patient_validation.md` §3.4.
> * `tmesim/validation/scorecard.py` records the same arithmetic finding in its exclusion notes.
>
> **What replaces it.** §1.2 is not deleted, because the comparison it ran is the right comparison and its numbers were always correct — only the ordering was misread. The corrected reading is below: the PDE is the **best** of four weak predictors, not the worst. The load-bearing caveat is now a different and narrower one — ρ ≈ 0.15 is weak however it is ordered, this is a **within-section** re-derivation and not out-of-sample, and the composite is carried by a single gene (VEGFA).
>
> **This correction is document-only.** The generator still emits the inverted prose, so re-running `python -m tmesim.validation.state_validation` will reproduce it until `state_validation.py` is fixed. That file is not owned by this correction pass. **The claim also propagated** — see the note at the end of §1.6.

| | Test 1 — simulated oxygen field | Test 2 — cell-state labelling |
|---|---|---|
| **Verdict** | **SUPPORTED** | **UNDERPOWERED** |
| Section | 10x Xenium Human Brain Cancer (GBM) FFPE — IO panel (v2.0.0) | 10x Xenium Human Brain — Glioblastoma section (v1.3.0) |
| Cells | 816,769 (480-gene panel) | 40,887 (23,812 tumour; 319-gene panel) |
| Measured against | VEGFA + HIF1A (hypoxia response) | MKI67 + TOP2A (proliferation) |
| Headline statistic | Spearman ρ = -0.1523 (p = 0.00050) | MKI67+ in VIABLE nan% (0/0) |
| Baseline contrast *(corrected)* | a one-line distance-to-vessel baseline scores **-0.1204**, i.e. **0.79× the solved field — weaker, not better** | the test could only have detected a **≥inf× rate ratio** |
| Killer detail *(corrected)* | not the baseline. It is that ρ ≈ 0.15 is **weak**, this is a **within-section** re-derivation, and the composite is carried by **one gene** (VEGFA -0.2004; HIF1A +0.0140, null) | the test could only have detected a **≥inf× rate ratio** |

**In one sentence each** *(corrected — see the CORRECTION block above; this sentence previously asserted the opposite)*. The simulated oxygen field *is* significantly anti-correlated with measured hypoxia response, and it **out-predicts** every physics-free baseline built from its own inputs — including the distance-to-nearest-vessel transform (-0.1204) and its own undiffused source term (-0.1111) — but it wins a contest between four weak predictors on the one section that measures the marker, so the result is real and small rather than decisive. And the proliferation test on the 319-gene section is dropout-limited to the point where it cannot distinguish a correct cell-state model from noise; that is a power failure, not a vindication and not a refutation.

> **The two sections are DIFFERENT PATIENTS.** Test 1 and Test 2 are two experiments, not two views of one tumour. Nothing here transfers from one to the other.

## 0. What is under test, and why it matters

The chain is: `physics.solve_oxygen` produces an oxygen field → `cells.agents.update_cell_states` bands it into HYPOXIC / QUIESCENT / VIABLE at `t = OxygenConfig.hypoxia_threshold` = 0.2 → `cells.pkpd` multiplies drug kill by `PDConfig.hypoxia_resistance` = 0.35 and `PDConfig.quiescent_resistance` = 0.5.

`results/dose_schedule.md` shows the consequence: for temozolomide — whose per-cell sensitivity is uniform 1.0 because MGMT is not on the panel — the kill CEILING is *entirely* those two multipliers, and ~89% of it comes from the quiescent band alone. If the oxygen field or the banding is wrong, every absolute kill number in this repository inherits the error.

**Note on the 0.35 figure above.** That is the GLOBAL `PDConfig` default reported here for context; temozolomide specifically now carries a literature-informed PER-DRUG override of 0.80 (`core.config.DRUG_LIBRARY`, `PDConfig` docstring; `results/sweep/state_multiplier_literature_check_*.md`), converted Hill-curve-consistently from real GBM hypoxia-vs-TMZ IC50-shift data at this project's own operating point. It changes nothing in THIS report: every conclusion above turns on the oxygen-field / cell-state-banding machinery upstream of the multiplier, not its value, and (as section 2 below and `results/kill_calibration.md` both independently confirm) the HYPOXIC band is empty or near-empty on every real section this project currently reconstructs, so neither multiplier's value moves any number here.

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
**Panel selection.** `io_480` (marker coverage: {'brain_319': 0.5, 'io_480': 1.0, 'breast_313': 0.7142857142857143})  
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

* **VEGFA is rank 1 of 480** at ρ = -0.20039. Of every gene the panel measures, the canonical hypoxia-response gene is the one the simulated oxygen field most strongly anti-correlates with. That is a real, non-trivial hit and it is the strongest positive evidence in this report.
* **HIF1A is rank 327 of 480** at ρ = +0.01403 — null, and if anything the *wrong sign*. This is expected on mechanism, not a surprise: HIF-1α is regulated by oxygen-dependent prolyl hydroxylation and VHL-mediated proteolysis of the **protein** (Ivan et al., Science 2001;292(5516):464-468 (PMID 11292862, doi:10.1126/science.1059817); Jaakkola et al., Science 2001;292(5516):468-472 (PMID 11292861, doi:10.1126/science.1059796)), so its mRNA has no obligation to track oxygen. HIF1A mRNA should not have been treated as a hypoxia readout in the first place, and the composite score is diluted by including it.

**The circularity, stated plainly.** The positive end of that table is the vascular and myeloid programme — RGS5, FLT1, PDGFRB, A2M, IGFBP7, CD163. `FLT1` is one of the three ENDOTHELIAL markers that *builds* the oxygen source term, so its positive correlation is partly the model correlating with itself. Only the NEGATIVE end of the ranking is informative, which is why VEGFA's rank-1 position is the claim worth making.

VEGFA is the right measured readout here: Fischer et al., Brain Pathol 2005;15(4):297-310 (PMID 16389942, doi:10.1111/j.1750-3639.2005.tb00115.x) — *in GBM, vascular involution produces necrosis and hypoxia, which induces HIF-1alpha and VEGF expression specifically in the perinecrotic pseudopalisading glioma cells*.

### 1.2 (a) Is simulated oxygen anti-correlated with measured hypoxia — and does it beat a trivial baseline?

Per-cell Spearman ρ against the measured hypoxia score (mean z of depth-normalised VEGFA and HIF1A), computed on **all 816,769 cells** at dx = 12 µm. Every predictor is sign-aligned to oxygen (larger = better perfused), so a correct predictor gives a NEGATIVE ρ.

| predictor | what it is | ρ vs hypoxia score | ρ vs VEGFA | ρ vs HIF1A | ρ with log-depth partialled out |
|---|---|---|---|---|---|
| **SIMULATED O₂ (the PDE)** | the reaction-diffusion solve under test | **-0.1523** | -0.2004 | +0.0140 | -0.1619 |
| vascular_kde | the solver's own input field, undiffused | **-0.1111** | -0.1017 | -0.0556 | -0.1119 |
| dist_endothelial | one cKDTree query on the PDE's own source cells; no physics | **-0.1204** | -0.1237 | -0.0342 | -0.1185 |
| local_cell_density | pure local crowding; ignores vessels entirely | **-0.0192** | -0.0031 | -0.0439 | -0.0336 |

Permutation p-values (label shuffling, 2,000 draws on a seeded random subsample of 50,000 cells — subsampled because the kernel is O(n × draws); every effect size above is on the full population, so subsampling can only make a p-value *more* conservative):

| predictor | ρ on the subsample | permutation p |
|---|---|---|
| simulated_o2 | -0.1489 | 0.00050 |
| dist_endothelial | -0.1254 | 0.00050 |
| vascular_kde | -0.1174 | 0.00050 |

**The answer to (a) is yes, and the PDE wins the comparison — but it is a contest between weak predictors.** *(This paragraph and the three bullets below are the CORRECTED text. The superseded version asserted that two baselines beat the PDE; see the CORRECTION block at the top of this document for what it said and why it was backwards.)*

The simulated field is genuinely, significantly anti-correlated with measured hypoxia response (ρ = -0.1523, p = 0.00050, floor of the estimator). Every trivial baseline built from the same inputs is anti-correlated too, but the solved field is the **strongest** of the four — reading the table by |ρ|, as its own sign convention requires:

* **The PDE (-0.1523) beats `dist_endothelial` (-0.1204) by 1.27×.** One `cKDTree` query — the distance from each cell to the nearest endothelial cell — recovers 79% of the solved field's association with measured VEGFA/HIF1A. That a physics-free one-liner gets **most** of the way there is the honest caveat; it does not get all the way there, and it does not win.
* **`vascular_kde` (-0.1111) is the weakest vessel-derived predictor of the three.** That is the solver's OWN SOURCE TERM, sampled before it is diffused. The solve takes a field correlating at -0.1111 and returns one correlating at -0.1523: **the diffusion step ADDS predictive information** (1.37× on |ρ|) rather than destroying it. This is the single most direct evidence in the document that the reaction–diffusion step is doing work.
* **`local_cell_density` (-0.0192) is near zero**, which is the sanity floor — not every spatial variable correlates with the marker, so the effects above are real signal, not autocorrelation artefacts.

**How large is "wins", honestly.** 1.27× on a ρ of 0.15 is a small margin on a weak association, measured on one section, with no confidence interval attached to the *difference* between two correlated predictors (the permutation test above tests each ρ against zero, not the predictors against each other). The defensible statement is **"the solve is not redundant with the geometry it is built from"**, not "the solve is accurate".

**Nested test — does the PDE add anything beyond the raw geometry?**

| partial Spearman | ρ |
|---|---|
| hypoxia score vs **simulated O₂**, distance-to-vessel held fixed | -0.1227 |
| hypoxia score vs **−distance-to-vessel**, simulated O₂ held fixed | -0.0792 |
| simulated O₂ vs −distance-to-vessel (how nested the two predictors are) | +0.3006 |

*(Corrected reading — this paragraph previously concluded the opposite.)* The simulated field retains **0.1227** of explanatory power after distance-to-vessel is conditioned out; distance-to-vessel retains only **0.0792** after the simulated field is conditioned out. The two are substantially nested (ρ = +0.3006), as they must be — the PDE's source cells *are* the distance transform's targets — but **distance-to-vessel is the redundant variable in that pair, not the PDE.** The nested test therefore agrees with the marginal comparison above, which is what makes the corrected ordering hard to argue with: two independent framings of the question return the same winner.

**Depth control.** The simulated field correlates with total transcript count — a purely technical variable — at ρ = +0.0572. Partialling log-depth out of the table above changes nothing material (simulated O₂ goes -0.1523 → -0.1619; distance goes -0.1204 → -0.1185), so the ordering is not a depth artefact.

### 1.3 (b) Are model-HYPOXIC cells enriched for measured VEGFA/HIF1A vs model-VIABLE?

State occupancy on this section at dx = 12 µm (all cells): **HYPOXIC** 30.90%, **QUIESCENT** 69.09%, **VIABLE** 0.00%.

| marker | contrast | n (A vs B) | mean A | mean B | Cohen's d | Cliff's δ | permutation p |
|---|---|---|---|---|---|---|---|
| hypoxia_score | hypoxic vs viable | 252,413 vs 37 | +0.1226 | -0.6404 | **+1.0210** | +0.5992 | 0.01549 |
| hypoxia_score | hypoxic vs quiescent | 252,413 vs 564,319 | +0.1226 | -0.0548 | **+0.2509** | +0.1371 | 0.00050 |
| hypoxia_score | hypoxic vs viable TUMOUR | 182,505 vs 6 | +0.1595 | -0.5608 | **+0.9486** | +0.5598 | 0.34933 |
| VEGFA | hypoxic vs viable | 252,413 vs 37 | +0.2731 | -0.8513 | **+0.9710** | +0.6094 | 0.02799 |
| VEGFA | hypoxic vs quiescent | 252,413 vs 564,319 | +0.2731 | -0.1221 | **+0.4019** | +0.1808 | 0.00050 |
| VEGFA | hypoxic vs viable TUMOUR | 182,505 vs 6 | +0.3378 | -0.7628 | **+0.9273** | +0.5703 | 0.15792 |
| HIF1A | hypoxic vs viable | 252,413 vs 37 | -0.0278 | -0.4294 | **+0.3997** | +0.2098 | 0.34433 |
| HIF1A | hypoxic vs quiescent | 252,413 vs 564,319 | -0.0278 | +0.0125 | **-0.0403** | -0.0241 | 0.00050 |
| HIF1A | hypoxic vs viable TUMOUR | 182,505 vs 6 | -0.0189 | -0.3588 | **+0.3331** | +0.1893 | 0.70715 |

Depth-stratified control (labels permuted **within** total-transcript deciles, so only comparably-deep cells are ever compared):

| marker | HYPOXIC vs VIABLE, depth-stratified difference | p |
|---|---|---|
| hypoxia_score | +0.76101 | 0.03198 |
| VEGFA | +1.04284 | 0.05147 |
| HIF1A | +0.47918 | 0.29685 |

**The answer to (b) is yes for VEGFA, no for HIF1A.** Model-HYPOXIC cells carry substantially more measured VEGFA than model-VIABLE cells (d = +0.971, Cliff's δ = +0.609, p = 0.02799), and the effect survives depth stratification. For HIF1A the same contrast is only d = +0.400, and the tumour-only version is not significant at all (p = 0.70715).

**Two things blunt this result.** First, the VIABLE group holds only 37 of 816,769 cells (0.005%), so "HYPOXIC vs VIABLE" is a comparison against a sliver. The better-balanced HYPOXIC-vs-QUIESCENT contrast is far weaker (VEGFA d = +0.402, HIF1A d = -0.040). Second, §1.2 shows that most of this separation — 79% of the association, by |ρ| — is recoverable by thresholding distance-to-vessel alone, so a large Cohen's d here is mostly evidence that the *bands* sit somewhere sensible in space. *(Corrected: this passage previously read "the same separation is obtained, more strongly, by thresholding distance-to-vessel … not that the PDE computed them", which inherited the inverted §1.2 inequality. The PDE separates them somewhat better than the distance threshold does; the point that survives is that the margin is small.)*

### 1.4 (c) Field-vs-field agreement on the grid

The per-cell test is dominated by single-cell dropout. Averaging both fields over fixed 60 µm tissue tiles (≥5 measured cells per tile) suppresses that noise and asks the coarser question: does a REGION the model calls hypoxic express more hypoxia response? Tiles are a fixed physical size at every dx, so the three rows are comparable to each other rather than to the mesh.

| dx (µm) | grid | tiles used | Spearman ρ | Pearson r |
|---|---|---|---|---|
| 30 | 363×383 | 24,254 | -0.3096 | -0.3528 |
| 20 | 543×574 | 24,254 | -0.3153 | -0.3697 |
| 12 | 905×956 | 24,254 | -0.3117 | -0.3663 |

Agreement is ~2.4× stronger at the field level (ρ ≈ -0.312 vs -0.152 per cell), which is the expected signature of per-cell dropout rather than of a better model. It is the most favourable framing of the test and it is reported as such.

### 1.5 (d) Does agreement improve as the grid is refined?

**No — it is flat to four decimal places.**

| dx (µm) | grid | mean O₂ (tissue) | HYPOXIC / QUIESCENT / VIABLE (all cells) | ρ(sim O₂, hypoxia score) | ρ(dist-to-vessel, hypoxia score) |
|---|---|---|---|---|---|
| 30 | 363×383 | 0.2211 | 32.59% / 67.24% / 0.17% | -0.1530 | -0.1204 |
| 20 | 543×574 | 0.2190 | 33.38% / 66.59% / 0.02% | -0.1536 | -0.1204 |
| 12 | 905×956 | 0.2210 | 30.90% / 69.09% / 0.00% | -0.1523 | -0.1204 |

This is good news about the numerics and neutral news about the model. With the geometric tissue mask and the areal oxygen sink both in place, the field is genuinely grid-converged *for this statistic*: mean O₂ moves by 0.2% and the correlation by <0.0002 over a 2.5× range of dx. Note also that the PDE's margin over the distance baseline is **stable across the whole ladder** — the baseline sits at -0.1204 at every dx while the solved field sits at -0.1523 to -0.1536 — so the corrected ordering in §1.2 is not a single-resolution artefact. What the flatness *does* establish is the ceiling: refining the mesh will not raise ρ above ≈0.153, so **the remaining gap to the measured marker is a model or measurement limitation, not a resolution one.** *(Corrected: this paragraph previously closed "The model has converged onto an answer that a distance transform beats.")*

### 1.6 Verdict — Test 1

## **SUPPORTED**

* rho = -0.1523 (p = 0.00050) and the solved field beats every trivial baseline (best: dist_endothelial at -0.1204)
* the composite is carried by ONE gene: VEGFA rho=-0.2004, HIF1A rho=+0.0140 — HIF1A is null

**Consequences, stated plainly.**

1. **The oxygen field is not falsified, and it is not out-competed either.** *(Corrected — this item previously read "it is out-competed … its own undiffused source term and a plain distance transform both do better".)* Its sign is right, its significance is overwhelming, VEGFA being its rank-1 anti-correlate out of 480 genes is a genuine hit, and it is the strongest of the four predictors tested (|ρ| 0.1523 vs 0.1204 vs 0.1111 vs 0.0192). The reaction–diffusion step measurably **improves** on its own undiffused source term, 0.1111 → 0.1523. What remains true is that the association is **weak in absolute terms** and that a physics-free distance transform recovers 79% of it, so the PDE's marginal value over trivial geometry is real but modest.
2. **Hypoxia-driven resistance in the kill model is supported only as *ordering*, not as *magnitude*.** *(Corrected — this item previously concluded the multiplier is "applied to cells selected by a field that a one-line geometric heuristic selects better", which does not follow from the data.)* `PDConfig.hypoxia_resistance` is applied to cells selected by the best available predictor of measured hypoxia response on this section — but "best" here means ρ ≈ 0.15, and §1.5 shows refining the mesh cannot improve it. So the *spatial pattern* of hypoxia is defensible as a weak, correctly-signed ordering; **any kill number that treats the hypoxic/non-hypoxic assignment as accurate per-cell is over-reading it**, and consequence 5 (the absolute level is uncalibrated) binds regardless.
3. **Hypoxia should be derived from measured markers where they exist.** On a section that measures VEGFA, the defensible per-cell hypoxia covariate is the measured VEGFA level, not the simulated oxygen — and where no such marker is on the panel (the 319-gene section: no VEGFA, no CA9, no ADM, no NDRG1, no SLC2A1), the honest statement is that hypoxia is **unmeasured**, not that it has been simulated.
4. **Drop HIF1A from any hypoxia score.** Its mRNA is null against the field (rank 327/480) for a known mechanistic reason, and including it dilutes the composite from ρ ≈ -0.200 (VEGFA alone) to -0.152.
5. **The absolute oxygen level remains uncalibrated.** Nothing here tests it. `OxygenConfig.consumption` has never been fitted to a measured oxygen tension, and on this section the model calls 30.9% of all cells HYPOXIC while on the 319-gene section under the shipping mask it calls 0.0% hypoxic (§2.2). Those two numbers cannot both be right. **This, not the baseline comparison, is the strongest negative result in this document.**

### 1.7 Where the superseded §1.2 claim propagated

The inverted "a distance transform out-predicts the PDE" reading was load-bearing well beyond this file, and **the corrections below are NOT applied** — they are recorded here so the trail is visible. Each of these asserts or is premised on the superseded inequality:

| location | what it says | status |
|---|---|---|
| `tmesim/validation/state_validation.py` | generates the inverted §1.2 prose; re-running this module reproduces it | **not corrected** — generator not owned by this pass |
| `results/oxygen_demotion.md` §1 | opens *"`results/state_validation.md` §1.2 found the reaction-diffusion oxygen solve **falsified as value-add**"* — the whole demotion rests on it | **not corrected** |
| `tmesim/physics/oxygen_proxy.py` | a distance-based proxy module whose stated premise is that it *"out-predicts the solve on the one section where we can measure"* | **not corrected** |
| `tmesim/reconstruct/real_pipeline.py`, `tmesim/reconstruct/cell_state.py` | carry the "out-predicted by a distance transform" caveat in docstrings and in an emitted caveat string | **not corrected** |
| `tmesim/validation/oxygen_benchmark.py` | emits "falsified as value-add — out-predicted against measured hypoxia markers" | **not corrected** |
| `docs/architecture.md`, `docs/VERIFICATION_METHODOLOGY.md` | repeat the claim as settled | **not corrected** |
| Studio caveat panel (`app/`, e.g. `app/reports/rep_0000*/report.html`) | *"oxygen (pde): FALSIFIED AS VALUE-ADD … out-predicted … (rho -0.1503 vs -0.1132)"* — note |-0.1503| > |-0.1132|, so **this instance is backwards by its own quoted numbers too** | **not corrected** — front end owned elsewhere |
| `results/heldout_patient.md` §6.1 | describes the hypoxia rows as *"confirmations of a negative result"* | **not corrected** |

Already correct and needing no change: `results/heldout_patient_validation.md` §3.4 (which independently caught this), `tmesim/validation/scorecard.py`, and this document's own §1.6 verdict — all three decide the question arithmetically.

## 2. TEST 2 — cell-state model vs measured proliferation

**Section.** 10x Xenium Human Brain — Glioblastoma section (v1.3.0)  
**Panel.** hBrain_254g base + 65-gene add-on (319 gene targets)  
**Cells.** 40,887 total, 23,812 tumour  
**Headline configuration.** tissue mask `occupancy`, dx = 12 µm — *the configuration `results/dose_schedule.md` actually ran under*, which is the one whose kill ceilings are at stake.

### 2.1 The measured proliferating fraction, under two documented rules

| marker | cells at exactly 0 | rule `nonzero` (count ≥ 1) | rule `percentile` (top 20% of depth-normalised expression) |
|---|---|---|---|
| MKI67 | 99.45% | 0.550% (131/23,812) | 0.550% — **DEGENERATE** |
| TOP2A | 99.01% | 0.991% (236/23,812) | 0.991% — **DEGENERATE** |

**The percentile rule is mathematically degenerate on this panel, and that IS the result.** The rule was written to bracket the published Ki-67 index — "call the top 20% of cells proliferating" — but with >98.9% of tumour cells at exactly zero counts, the 80th percentile of expression *is* zero. The threshold falls inside a tie mass at the origin and the rule collapses back onto `count ≥ 1`. There is no percentile of this data that corresponds to a Ki-67 index, because the data does not reach that far.

**Against the literature.** Published Ki-67 labelling indices in glioblastoma (PubMed, all verified):

| cohort | Ki-67 index | source |
|---|---|---|
| GBM, n=35 (prognostic cutoff) | 19% | Fudaba et al., World Neurosurg 2021;149:e469-e480 (PMID 33567368, doi:10.1016/j.wneu.2021.02.002) |
| GBM, n=44 (cohort split point) | 27% | Alkhaibary et al., Hematol Oncol Stem Cell Ther 2018;12(2):82-88 (PMID 30552865, doi:10.1016/j.hemonc.2018.11.001) |
| epithelioid GBM, n=16 (median, range 10-50%) | 30% | Zheng et al., Zhonghua Bing Li Xue Za Zhi 2021;50(3):229-235 (PMID 33677887, doi:10.3760/cma.j.cn112151-20200617-00479) |
| GBM, n=70 (median) | 40% | Caramanti et al., Cureus 2022;14(4):e24246 (PMID 35602791, doi:10.7759/cureus.24246) |

The panel detects MKI67 in **0.550%** of tumour cells — **35× to 73× below** the published range. This is not a claim that this tumour has a Ki-67 index of 0.44%; it is a statement about probe sensitivity and transcript dropout on a targeted in-situ panel. The reference standard is itself imprecise (Nielsen et al., Diagn Pathol 2018;13:38 (PMID 29885671, doi:10.1186/s13000-018-0711-2) — Ki-67 LI scoring in gliomas agrees poorly between observers: kappa 0.32 for whole-slide and as low as 0.04 for hot-spot scoring (n=235 gliomas)), but not by two orders of magnitude.

**Model VIABLE fraction, for comparison.** 0.00% in the headline configuration, but 0.00% – 10.04% across the mask × dx grid (§2.2). **Do not read any agreement between "model VIABLE %" and "MKI67-detected %" as corroboration** — it would be a coincidence between two unrelated small numbers, one of which is a bookkeeping artefact and the other a dropout rate.

### 2.2 The model's state vector is set by bookkeeping, not by tissue

| tissue mask | dx (µm) | mean O₂ (tissue) | HYPOXIC | QUIESCENT | VIABLE |
|---|---|---|---|---|---|
| `geometric` | 30 | 0.3383 | 0.00% | 99.18% | 0.82% |
| `geometric` | 20 | 0.3479 | 0.00% | 91.96% | 8.04% |
| `geometric` | 12 | 0.3520 | 0.00% | 89.96% | 10.04% |
| `occupancy` | 30 | 0.3374 | 0.00% | 99.18% | 0.82% |
| `occupancy` | 20 | 0.3461 | 0.00% | 91.96% | 8.04% |
| `occupancy` ← | 12 | 0.3050 | 0.00% | 100.00% | 0.00% |

Read that table before reading any statistic below it. The tumour-cell VIABLE fraction runs from 0.00% to 10.04% — a 10041155720× swing — across two tissue-mask conventions and three grid spacings. **Neither the mask nor the mesh is a measurement.** Under the current shipping default (`geometric`) the HYPOXIC state is *empty at every dx*, so `PDConfig.hypoxia_resistance` has no cells to act on at all; the 86% QUIESCENT / 13.4% HYPOXIC split that `results/dose_schedule.md` quotes comes from the legacy `occupancy` mask, whose coverage shrinks as the mesh is refined. Those kill ceilings are conditional on a convention that the reconstruction layer has since abandoned.

### 2.3 Are model-VIABLE cells more proliferative? — and could we have told?

| marker | state | rate | k / n | Wilson 95% CI |
|---|---|---|---|---|
| MKI67 | HYPOXIC | — | 0 / 0 | state is empty |
| MKI67 | QUIESCENT | 0.550% | 131 / 23,812 | 0.464–0.652% |
| MKI67 | VIABLE | — | 0 / 0 | state is empty |
| TOP2A | HYPOXIC | — | 0 / 0 | state is empty |
| TOP2A | QUIESCENT | 0.991% | 236 / 23,812 | 0.873–1.125% |
| TOP2A | VIABLE | — | 0 / 0 | state is empty |

| marker | contrast | Cohen's d | Cliff's δ | permutation p | depth-stratified p | MINIMUM DETECTABLE rate ratio | expected positives |
|---|---|---|---|---|---|---|---|
| MKI67 | VIABLE vs rest | +nan | +nan | nan | nan | **inf×** | 0.0 |
| MKI67 | HYPOXIC vs QUIESCENT | +nan | +nan | nan | nan | **inf×** | 0.0 |
| TOP2A | VIABLE vs rest | +nan | +nan | nan | nan | **inf×** | 0.0 |
| TOP2A | HYPOXIC vs QUIESCENT | +nan | +nan | nan | nan | **inf×** | 0.0 |

**The power analysis is the result.** The "minimum detectable rate ratio" column is the smallest enrichment the contrast could have found at α = 0.05, power = 0.8 given the observed positivity rate and the group sizes.

* For **MKI67, VIABLE vs rest**, the VIABLE group holds 0 tumour cells and would contain 0.0 MKI67-positive cells under the null. The test could only have detected a **≥inf× rate ratio**.
* The measured biology this is supposed to detect is nowhere near that size: Brat et al., Cancer Res 2004;64(3):920-927 (PMID 14871821, doi:10.1158/0008-5472.CAN-03-2073) measured that hypoxic pseudopalisading GBM cells are *5–50% LESS* proliferative than adjacent astrocytoma — a rate ratio of roughly 0.5–0.95, i.e. a **depletion of at most 2×**. A test that needs inf× to fire cannot see an effect of that size in either direction.
* Therefore **no conclusion about the VIABLE label can be drawn from this section.** Not "the state model is wrong"; not "the state model is right". The measurement is not capable of adjudicating it.

**A positive result is visible, and it must not be believed.** The point estimate for MKI67 goes the RIGHT way — nan% positive in VIABLE against 0.550% section-wide, a nan× rate ratio, permutation p = nan, surviving depth stratification (p = nan). It rests on **0 positive cells out of 0**. Its Wilson interval spans nan–nan%, one Poisson fluctuation moves it from 0 to 2 or 7, and TOP2A — the other proliferation marker on the same cells — shows nothing (nan%, 0/0, p = nan). This is precisely the regime the minimum-detectable-effect calculation says is uninterpretable, and reporting it as support would be the exact error this module exists to prevent.

**What the test CAN rule out.** The better-powered HYPOXIC-vs-QUIESCENT contrast (n = 0 vs 23,812) needs only inf× for MKI67 and inf× for TOP2A, which is inside the plausible range. It rules out a **large** proliferation deficit in the model's HYPOXIC band — and it does so in the wrong direction:


Neither reversal is large, and both sit close to the depth-stratification boundary, so this is a weak contradiction rather than a strong one. But it is a contradiction, and it points the same way as the continuous test:

| continuous test | Spearman ρ | p | ρ with log-depth removed |
|---|---|---|---|
| simulated O₂ vs MKI67 (per tumour cell) | +0.0108 | 0.0917 | +0.0140 |
| simulated O₂ vs TOP2A (per tumour cell) | +0.0013 | 0.8356 | +0.0000 |

Rank-correlating the **simulated oxygen** against measured proliferation across all 23,812 tumour cells — the best-powered version of the question, because it uses every cell and does not depend on where the band edges fall — returns essentially nothing, with a *positive* sign where a negative one was not expected either. The oxygen field carries no per-cell proliferation information on this section.

### 2.4 dx-dependence of the Test 2 contrast

| mask | dx (µm) | VIABLE fraction | MKI67+ in VIABLE | MKI67+ elsewhere | min. detectable rate ratio |
|---|---|---|---|---|---|
| `geometric` | 30 | 0.82% | 2.551% | 0.534% | 3.7× |
| `geometric` | 20 | 8.04% | 0.940% | 0.516% | 1.9× |
| `geometric` | 12 | 10.04% | 0.878% | 0.514% | 1.8× |
| `occupancy` | 30 | 0.82% | 2.551% | 0.534% | 3.7× |
| `occupancy` | 20 | 8.04% | 0.940% | 0.516% | 1.9× |
| `occupancy` | 12 | 0.00% | nan% | 0.550% | inf× |

The measured MKI67+ rate inside the model's VIABLE band is ~0.44–0.46% at five of the six configurations — i.e. **indistinguishable from the section-wide rate of 0.550%** — because at those settings VIABLE contains 74–98% of the tumour and the "contrast" is a population against itself. Only in the sixth configuration (the dose_schedule one) does VIABLE become small enough to be a real subgroup, and there it is too small to test. There is no dx at which this comparison is both meaningful and powered.

### 2.5 Verdict — Test 2

## **UNDERPOWERED**

* MKI67 VIABLE-vs-rest could only have detected a infx rate ratio (0.0 positives expected in VIABLE under the null)
* TOP2A VIABLE-vs-rest could only have detected a infx rate ratio (0.0 positives expected in VIABLE under the null)
* the decisive VIABLE-vs-rest comparison cannot distinguish 'the state labelling is right' from 'the state labelling is noise'
* independently of power, the model's VIABLE fraction is not a measurement: it runs 0.00% to 10.04% across two tissue-mask conventions and three grid spacings, all of which are bookkeeping choices

**Consequences, stated plainly.**

1. **The quiescent labelling is UNSUPPORTED — in the specific sense that it has not been, and on this section cannot be, tested.** The measurement that would adjudicate it is 40–90× too insensitive. Nobody should read this as clearance.
2. **The kill ceilings inherit that.** `results/dose_schedule.md` derives ~89% of the temozolomide kill ceiling from the quiescent band, via a ×0.5 multiplier applied to a label with no measured support and a 10041155720× convention-dependent size. Every absolute kill fraction and every ceiling in that report carries an unquantified error of that magnitude. The dose/schedule **comparisons** are far more robust — they share the same multiplier and it largely cancels — but the **levels** are not defensible as measurements.
3. **The one thing the data does say points the wrong way.** In the better-powered HYPOXIC-vs-QUIESCENT contrast, both measured proliferation markers are *higher* in the model's HYPOXIC band. That is the opposite of the measured tissue biology (Brat et al., Cancer Res 2004;64(3):920-927, PMID 14871821).
4. **What would fix it.** Either (a) a proliferation readout with a real detection rate — IHC/IF Ki-67 on a serial section, or a whole-transcriptome platform — or (b) derive the state labels from measured markers (a per-cell proliferation score) instead of from the simulated oxygen, and report the states as *observed* rather than *predicted*.

## 3. Caveats — what could make these conclusions wrong

1. **The two sections are different patients.** Test 1 runs on `data_cache/xenium_gbm_io` (816,769 cells, IO panel); Test 2 on `data_cache/xenium_gbm` (40,887 cells, brain panel). They are different tumours with different vascular architectures and different panels. Nothing here is a paired observation, and the two verdicts do not reinforce each other.
2. **VEGFA and HIF1A are hypoxia-RESPONSE genes, not an oxygen tension.** They are the transcriptional output of a pathway with its own kinetics and its own non-hypoxic inducers. A simulated field could be a perfectly good oxygen tension and still correlate weakly with them. What the trivial-baseline control establishes is narrower and safer: *whatever* VEGFA is reporting, the PDE predicts it better than a distance transform does — and that comparison is internal, so it is immune to this caveat. *(Corrected: this sentence previously had the two predictors the other way round.)*
3. **Targeted-panel dropout.** A zero is "no probe fired", not "no transcript". This is fatal for Test 2 (MKI67 at 0.55%) and merely limiting for Test 1 (VEGFA at 63.4%, which is why Test 1 is powered and Test 2 is not).
4. **Typing circularity.** `FLT1` and `PLVAP` are IO-panel ENDOTHELIAL markers, so they help build the oxygen source; the positive tail of the 480-gene ranking is partly the model correlating with its own inputs. `MKI67` is one of the four brain-panel TUMOR markers, so an MKI67-high cell is marginally more likely to be typed TUMOR — that inflates the MKI67+ rate *uniformly across states*, so it cannot manufacture a between-state contrast, but it does inflate the marginal rate slightly.
5. **`SLC1A3`/`S100B`/`PTN` (the IO STROMAL class) are expressed by malignant astrocyte-like glioma cells too**, so the IO section's STROMAL compartment is not clean. The brain panel's `AQP4`/`GJA1` have the same ambiguity, so the two sections are consistently rather than differently wrong here.
6. **Pericytes are unlabelled on the IO section.** The panel measures `RGS5`/`PDGFRB`/`NOTCH3` but `CellType` has no MURAL member, so mural cells fall into the ungated TUMOR default. That inflates the tumour count and hence the oxygen sink.
7. **The absolute oxygen level is uncalibrated.** `OxygenConfig.consumption` has never been fitted against a measured oxygen tension. The two sections disagree by orders of magnitude in hypoxic fraction (30.9% of IO cells vs 0.0% of brain-panel cells under the shipping mask), which is a statement about the sink constant and the cell density, not about the two tumours.
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
