# Demoting the oxygen PDE: the benchmark, the cost, and the downstream effect

`results/state_validation.md` §1.2 found the reaction-diffusion oxygen solve **falsified as value-add** — out-predicted against measured hypoxia markers by a one-line distance transform built from its own source cells, and by its own undiffused source. This report acts on that. It benchmarks every oxygen route the codebase can build, prices them, measures whether the choice changes a drug prediction at all, and makes a recommendation.

**The PDE has not been deleted.** It is still the default (`reconstruct_real(oxygen_mode='pde')`, bit-identical to the pre-option pipeline). What has changed is that there is now an honest alternative, it is labelled a SURROGATE, and the default's own falsification travels with it in `reconstruction_report()['oxygen']['caveat']`.

## Headline

| | |
|---|---|
| **Best predictor of measured VEGFA** | `vessel_distance_exact` at rho = -0.1647 (p = 0.00050) |
| **The PDE** | rho = -0.1479 (p = 0.00050), rank 5 of 7 |
| **Cost ratio** | the PDE is 11x to 31x slower than the surrogate |
| **Downstream effect on predicted kill** | at most 0.0262 absolute across three drugs |
| **Downstream effect on therapeutic coverage** | 0 — coverage is oxygen-invariant *by construction* |
| **Recommendation** | keep `'pde'` as the default; use `oxygen_mode='vessel_distance'` for any NEW hypoxia claim. See §5. |

> **The two sections are DIFFERENT PATIENTS.** The benchmark runs on the 480-gene IO section (816,769 cells); the downstream drug arm runs on the 319-gene brain section (40,887 cells), because that is the section every flagship drug result in `results/` was computed on. Nothing transfers between them.

## 1. Benchmark — every oxygen route against measured VEGFA

**Section.** `data_cache/xenium_gbm_io` — 816,769 cells, 480-gene IO panel, dx = 12 um, grid 905x956.  
**Target.** Depth-normalised VEGFA (detected in 63.36% of cells). HIF1A (54.93%) is reported only inside the historical composite column — `state_validation.md` §1.1 showed it ranks 249/480 at rho = +0.0085 and dilutes the signal, for the mechanistic reason that oxygen sensing acts on the HIF-1alpha PROTEIN, not its mRNA (Ivan et al., Science 2001;292(5516):464-468, PMID 11292862, doi:10.1126/science.1059817; Jaakkola et al., Science 2001;292(5516):468-472, PMID 11292861, doi:10.1126/science.1059796).

Every predictor is **sign-aligned to oxygen** (larger = better perfused), so a correct predictor gives a NEGATIVE rho. All rho are on the full 816,769 cells; permutation p uses a fixed seeded subsample of 50,000 cells with 2,000 draws (the kernel is O(n x draws); subsampling the null can only make a p more conservative).

| rank | route | rho vs VEGFA | permutation p | rho vs VEGFA+HIF1A composite | rho with log-depth removed | what it is |
|---|---|---|---|---|---|---|
| 1 | **vessel_distance_exact** | **-0.1647** | 0.00050 | -0.1503 | -0.1663 | identical rho to raw distance-to-vessel: the map is monotone |
| 2 | **dist_endothelial** | **-0.1647** | 0.00050 | -0.1503 | -0.1663 | one cKDTree query on the PDE's own source cells; no physics |
| 3 | **vessel_distance** | **-0.1636** | 0.00050 | -0.1489 | -0.1651 | one cKDTree query per grid cell; no solve. This is what oxygen_mode='vessel_distance' installs |
| 4 | **vessel_network** | **-0.1521** | 0.00050 | -0.1366 | -0.1535 | more physical vessel definition; discards singleton endothelium |
| 5 | **pde** | **-0.1479** | 0.00050 | -0.1132 | -0.1548 | sparse steady-state solve; the current shipping default |
| 6 | **vascular_kde** | **-0.1456** | 0.00050 | -0.1270 | -0.1476 | the solver's own input field, undiffused |
| 7 | **local_cell_density** | **-0.0031** | 0.59770 | -0.0192 | -0.0093 | pure local crowding; ignores vessels entirely |
| — | ~~markers~~ (CIRCULAR) | -0.6346 | 0.00050 | -0.4556 | -0.6411 | CIRCULAR against this target — the only HYPOXIA_MARKER_GENES member on this panel IS VEGFA, so this correlates VEGFA with a smoothed copy of itself. Reported as a CEILING, not a competitor |

**Reading the table.**

* **The vessel-distance surrogate wins.** `vessel_distance_exact` at rho = -0.1647 against the PDE's -0.1479 — a 1.11x larger effect, from a kd-tree query instead of a sparse linear solve. Ranks 1 and 2 are the SAME predictor written two ways: `vessel_distance_exact` is `exp(-d/L)` and `dist_endothelial` is `-d`, and a monotone map cannot change a Spearman rho. They are both listed only to make that identity visible rather than claimed.
* **Meshing costs almost nothing here.** The grid-sampled surrogate (-0.1636) and the exact per-cell surrogate (-0.1647) differ by 0.0012, so binning the distance to a 12 um mesh is not what separates the routes.
* **The 'more physical' vessel definition is WORSE.** Measuring distance to a SEGMENTED vessel object rather than to any endothelial cell scores -0.1521 against -0.1636. Segmentation discards 4,003 of 34,179 endothelial cells as singletons, and on this section those discarded cells were carrying signal. It is offered as an option and is NOT the default — the plainer definition is the one that was benchmarked and the one that wins.
* **The mapping cannot be tuned to win.** `exp(-d/L)` is strictly monotone in distance, so **every** monotone map of distance gives the identical Spearman rho. The constant `L` sets where the HYPOXIC/QUIESCENT/VIABLE band edges fall and nothing else. This is why it is anchored to published morphometry rather than fitted — fitting it would be unfalsifiable decoration.
* **The marker route is a CEILING, not a competitor.** rho = -0.6346. The only `HYPOXIA_MARKER_GENES` member on this panel is VEGFA (absent: CA9, ADM, NDRG1, SLC2A1, LDHA, PGK1, ANGPTL4), so this correlates VEGFA with a kernel-smoothed copy of itself. It measures how much of the per-cell marker signal survives 40 um smoothing — i.e. an upper bound on what ANY spatial field could score here — and it must never be quoted as a validated predictor.

### 1.1 Calibration of the surrogate mapping

| quantity | value | provenance |
|---|---|---|
| hypoxia onset distance | 150 um | Thomlinson & Gray, Br J Cancer 1955;9(4):539-549 (PMID 13304213, doi:10.1038/bjc.1955.55) — tumour cords viable to ~100-180 um from stroma, necrotic beyond |
| hypoxia threshold | 0.20 | `OxygenConfig.hypoxia_threshold`, unchanged |
| derived length scale L | 93.20 um | `L = onset / ln(1/threshold)`, so a cell at the published onset distance lands exactly on the existing hypoxia threshold — no new tuning knob |
| VIABLE within | 85.4 um | `L * ln(1/(2*threshold))`, the QUIESCENT/VIABLE edge |
| the PDE's OWN length scale | 310 um | `sqrt(D / (consumption * rho_tum))` at the mean tumour areal density |

**The last row is a finding, not a footnote.** The solved field's own effective diffusion length is ~310 um, 2.1x the measured oxygen-diffusion distance. `OxygenConfig.consumption` has never been fitted to an oxygen tension, and at its current value the solve low-pass-filters its own source at a length scale the tissue does not have. That is a mechanism for §1.2's observation that the diffusion step DESTROYS predictive information: it is not adding physics, it is blurring. Independent support that vessel distance is the right axis in the first place: Helmlinger et al., Nat Med 1997;3(2):177-182 (PMID 9018236, doi:10.1038/nm0297-177) measured interstitial pO2 profiles BETWEEN ADJACENT VESSELS at <=10 um resolution in vivo — and, honestly, also found heterogeneous profile shapes and hypoxic well-perfused vessels, which is the ceiling on what geometry alone can know.

### 1.2 Nested test on the VEGFA-only target

| partial Spearman | rho |
|---|---|
| VEGFA vs **simulated O2**, distance-to-vessel held fixed | -0.0983 |
| VEGFA vs **-distance-to-vessel**, simulated O2 held fixed | -0.1224 |
| simulated O2 vs -distance-to-vessel (how nested the pair is) | 0.3459 |

Distance keeps more explanatory power after the PDE is conditioned out than the PDE keeps after distance is conditioned out. **The PDE is the redundant member of the pair**, reproducing §1.2's conclusion on the cleaner VEGFA-only target.

## 2. Runtime — a sparse solve against a kd-tree query

Best of 3 on the same reconstructed microenvironment (`data_cache/xenium_gbm_io`), so the shared reconstruction cost is excluded. The surrogate's timing INCLUDES building the kd-tree.

| dx (um) | grid | unknowns | PDE solve (s) | vessel-distance (s) | speed-up |
|---|---|---|---|---|---|
| 30 | 363x383 | 139,029 | 1.441 | 0.1275 | **11x** |
| 20 | 543x574 | 311,682 | 4.394 | 0.2571 | **17x** |
| 12 | 905x956 | 865,180 | 20.428 | 0.6556 | **31x** |

The gap widens with refinement, as it must: the sparse direct solve's cost grows superlinearly in the number of unknowns while the kd-tree query grows as `O(n_grid log n_vessels)`. At the production dx = 12 um the solve is the single most expensive step in a reconstruction and the surrogate is 31x cheaper.

## 3. Downstream — does the oxygen route change a drug prediction?

**Section.** `data_cache/xenium_gbm` — 40,887 cells, 319-gene brain panel, dx = 12 um, grid 424x675, geometric tissue mask (the shipping default). Perivascular-sink alpha = 148.61, calibrated once on doxorubicin as in every flagship run.

**How arm B is computed.** `sim.engine.simulate_response` re-solves the PDE internally and `sim/engine.py` is not modified here, so arm B re-bands the cells with the surrogate via the unmodified `cells.agents.update_cell_states` and re-runs the unmodified `cells.pkpd.apply_drug_kill` on **the engine's own transport result**, with the engine's own seeded generator. This is valid because nothing in `physics.pressure`, `physics.diffusion` or `treatment.drugs` reads `fields.oxygen` — oxygen reaches the kill through exactly one channel, the state banding. That is not assumed: arm A is recomputed through the same path and **must reproduce the engine's `kill_fraction` bit-for-bit or the run aborts**. It did, for all three drugs.

Mean oxygen over tissue: PDE 0.5755, surrogate 0.6170. The two fields agree spatially at Spearman rho = 0.2433 over the grid.

### 3.1 Tumour-cell state occupancy under each route

| route | HYPOXIC | QUIESCENT | VIABLE |
|---|---|---|---|
| `oxygen_mode='pde'` | 0.00% | 5.08% | 94.92% |
| `oxygen_mode='vessel_distance'` | 0.46% | 10.19% | 89.35% |

The banding is identical for every drug (it is set before any drug is applied), so one row each suffices. The resistance multipliers this feeds are `hypoxia_resistance = 0.35` and `quiescent_resistance = 0.5`.

**Note the PDE's HYPOXIC band is EMPTY on this section under the shipping `geometric` tissue mask** — a point `results/state_validation.md` §2.2 already made. So `PDConfig.hypoxia_resistance` currently has no cells to act on at all in the flagship configuration, and the entire pde-vs-surrogate difference below runs through the QUIESCENT band. Any statement anywhere in `results/` about hypoxia-driven resistance on this section is, under this mask, describing an empty set.

### 3.2 Predicted kill and therapeutic coverage

| drug | IC50 (uM) | kill (pde) | kill (vd) | delta kill | delta % | cov@1xIC50 (pde) | cov@1xIC50 (vd) | median E/IC50 (pde) | median E/IC50 (vd) |
|---|---|---|---|---|---|---|---|---|---|
| temozolomide | 606 | 0.1151 | 0.1113 | -0.0039 | -3.4% | 0.0000 | 0.0000 | 0.177 | 0.177 |
| cisplatin | 24.6 | 0.9125 | 0.8862 | -0.0262 | -2.9% | 0.9986 | 0.9986 | 10.06 | 10.06 |
| doxorubicin | 0.5 | 0.7531 | 0.7371 | -0.0160 | -2.1% | 0.9692 | 0.9692 | 3.003 | 3.003 |

`median E/IC50` is carried alongside coverage because coverage@1xIC50 is degenerate at one end or the other for these drugs on this section (temozolomide 0.0000 — no tumour cell reaches its IC50; cisplatin and doxorubicin near-saturated above 0.96), so an identical pair of coverage numbers could be read as a saturation artefact. The median potency margin is un-saturated, spans two orders of magnitude across the panel, and is *also* identical to every digit across the two arms.

**Therapeutic coverage does not move at all, and that is structural.** `delivery_metrics.therapeutic_coverage` is a pure function of the per-cell AUC and the drug's IC50, and the AUC field has no oxygen dependence anywhere in the pipeline. So **no oxygen model of any kind — the PDE, the surrogate, a perfect measurement — can change a delivery metric in this codebase.** Every coverage number already published is therefore immune to this whole question, which is the most useful thing in this report for anyone auditing `results/delivery_metrics.md`, `results/dose_schedule.md` or `results/vessel_network_impact.md`.

**Predicted kill moves a little, and always downwards**: -0.0262 to -0.0039 absolute, i.e. 3.4% relative at worst. The direction is not a coincidence — the surrogate bands more tumour into the resistant QUIESCENT band than the PDE does (§3.1), and every resistant cell carries a 0.5x kill multiplier.

**The cross-drug KILL ranking is unchanged BY THE OXYGEN ROUTE**: cisplatin > doxorubicin > temozolomide under the PDE, cisplatin > doxorubicin > temozolomide under the surrogate — identical. That is the expected structure: the resistance multipliers are shared across drugs, so they largely cancel out of a comparison and do not cancel out of a level.

**Read the reason, and do not over-extend the claim.** "Shared across drugs, so it cancels out of a comparison" is exactly why this test cannot certify a DRUG-SPECIFIC perturbation. The vascular delivery law is one: correcting the one-way source `S = PS*Cp` to the gradient-driven `S = PS*(Cp - C)` re-weights each drug by its own `bbb_penetrance`, which spans 47x across the GBM panel, and it does change the GBM top-1 (cisplatin -> vincristine) and the top-3 set on both panels. This is a KILL ranking; the delivery (AUC) ranking is more stable than it (breast: rho 1.000 for AUC against 0.830 for kill). See [`results/transvascular_exchange_revalidation.md`](transvascular_exchange_revalidation.md) §4.1.

**So: how load-bearing is the oxygen field?** *Weakly.* A 2-3% relative shift in absolute kill is far smaller than the acknowledged uncertainty in those levels — `results/state_validation.md` §2.2 shows the VIABLE fraction on this section swings 165x across two tissue-mask conventions and three grid spacings, all of which are bookkeeping choices. The oxygen ROUTE is a second-order knob next to the tissue-mask convention. It is not zero, though, so absolute kill fractions must be reported with the oxygen mode attached, which is why `reconstruction_report()['oxygen']['mode']` now exists.

### 3.3 What this does NOT show

Neither section carries a paired drug-response outcome label, so nothing here says which arm's kill is *right*. The comparison is internal: it measures how much the answer depends on a modelling choice, which is a sensitivity, not a validation.

## 4. What was wired, and what did not change

| | |
|---|---|
| new module | `tmesim/physics/oxygen_proxy.py` — `oxygen_from_vessel_distance`, `vessel_distance_um`, `hypoxia_from_markers`, `diffusion_length_um`, `pde_diffusion_length_um` |
| new option | `reconstruct_real(..., oxygen_mode=...)` with `'pde'` / `'vessel_distance'` / `'markers'`, plus `oxygen_kwargs=None` |
| default | **`'pde'`, unchanged and BIT-IDENTICAL** — asserted by `tests/test_oxygen_proxy.py` |
| recorded in | `reconstruction_report()['oxygen']` -> `mode`, `source`, `is_surrogate`, `caveat`, `detail` (and a flat `['oxygen_mode']`) |
| NOT touched | `tmesim/core`, `tmesim/sim/engine.py`, `tmesim/physics/oxygen.py` |

`oxygen_mode='markers'` **raises** on the 319-gene brain panel rather than degrading silently, because that panel carries no VEGFA, CA9, ADM, NDRG1 or SLC2A1: on that section hypoxia is UNMEASURED, and the honest interface says so. `hypoxia_from_markers` returns `None` with the same reason when called directly.

A caller who passes `oxygen_mode='vessel_distance'` straight into `sim.engine.simulate_response` will NOT get a surrogate-driven kill — the engine re-solves the PDE internally. That is documented in `reconstruct_real`'s docstring, and the supported route is the re-banding path §3 uses.

## 5. RECOMMENDATION

**Keep `'pde'` as the shipping default. Use `oxygen_mode='vessel_distance'` for every NEW hypoxia claim, and stop making per-cell hypoxia claims from the solved field.**

The reasoning, in order of weight:

1. **The surrogate is better and cheaper on the one section where we can measure** — rho -0.1647 vs -0.1479, at 31x less compute at dx = 12 um. There is no accuracy argument for the solve.
2. **But the drug outputs barely notice, and flipping the default would invalidate the provenance of everything already published.** §3 shows therapeutic coverage is EXACTLY invariant and kill moves by at most 0.0262 absolute (3.4% relative) with the cross-drug KILL ranking unchanged BY THIS ROUTE (a shared-multiplier perturbation; it says nothing about the drug-specific delivery law, which does change the GBM top-1 — [`results/transvascular_exchange_revalidation.md`](transvascular_exchange_revalidation.md) §4.1). The measurable gain on the drug side is therefore close to nil, while the cost is re-deriving every kill number in `results/` (§5.1) and breaking the correspondence between those documents and the code. Silent improvement is still silent change.
   The honest asymmetry: the surrogate is unambiguously better as a *hypoxia predictor* (§1) and only marginally different as a *drug-kill input* (§3). So the recommendation splits along exactly that line — flip it for hypoxia claims, leave it for the kill pipeline until that pipeline is being re-run for some other reason anyway.
3. **The measured route is the right long-run answer and is not available yet.** `hypoxia_from_markers` is the only route here that is a measurement, but the flagship drug section's panel does not carry a hypoxia-response gene, and on the IO section that does, benchmarking it is circular. A panel or platform carrying CA9/ADM/NDRG1/SLC2A1 alongside a drug-relevant readout would settle this properly.

### 5.1 If the default WERE flipped, this is the re-run list

Everything whose numbers depend on the per-cell state banding — i.e. on `PDConfig.hypoxia_resistance` / `quiescent_resistance` having cells to act on. Coverage-only and AUC-only results are immune (§3.2) and would not need re-running:

| artefact | why | immune? |
|---|---|---|
| `results/dose_schedule.md` | its kill ceilings are ~89% the quiescent multiplier | **no — re-run** |
| `results/state_validation.md` | Test 1 IS this comparison | **no — re-run** |
| `results/report.md`, `results/patient_report_*.md` | quote absolute kill | **no — re-run** |
| `results/ranking_robustness.md`, `results/uncertainty.md` | resample kill | **no — re-run** |
| `results/delivery_metrics.md` | AUC + coverage only | yes |
| `results/vessel_network_impact.md` | coverage + penetration length | mostly — only its `hypoxia` block moves |
| `results/numerical_rigor.md` | would need the surrogate's own dx sweep | **no — re-run** |

### 5.2 What would change our mind

* A section with a **measured oxygen tension** (phosphorescence-quenching pO2, or pimonidazole IHC on a serial section). Every route here is uncalibrated in absolute level; a tension would adjudicate levels, not just orderings, and could rehabilitate the PDE by fitting `consumption`.
* A panel carrying **several** hypoxia-response genes, so the marker route can be benchmarked non-circularly by leave-one-gene-out.
* Evidence that the solve matters somewhere the distance transform cannot reach — e.g. a geometry with a supply-limited vessel where the (1-O) saturation term does real work. Nothing in the two cached sections exercises that.

## 6. Caveats — what could make this wrong

1. **VEGFA is a hypoxia-RESPONSE gene, not an oxygen tension.** It is the transcriptional output of a pathway with its own kinetics and its own non-hypoxic inducers (Fischer et al., Brain Pathol 2005;15(4):297-310, PMID 16389942, doi:10.1111/j.1750-3639.2005.tb00115.x). A field could be a perfectly good oxygen tension and still score poorly here. What the comparison establishes is narrower and immune to this: *whatever* VEGFA is reporting, a distance transform predicts it better than the solve does, and that is an internal comparison on identical inputs.
2. **The two sections are different patients**, with different vascular architectures and different panels. The benchmark (§1) and the downstream effect (§3) are two experiments, not two views of one tumour.
3. **Circularity in the vessel term.** `FLT1` and `PLVAP` are IO-panel ENDOTHELIAL markers, so the cells that define 'distance to a vessel' are themselves typed from the panel. Both the PDE and the surrogate inherit this equally, so it cannot explain the gap between them — but neither route's absolute rho is free of it.
4. **Targeted-panel dropout.** A zero is 'no probe fired', not 'no transcript'. VEGFA fires in 63.4% of cells, which is enough to power this test and not enough to make any single cell's value trustworthy.
5. **No absolute calibration anywhere.** `OxygenConfig.consumption` has never been fitted to a measured oxygen tension, and the surrogate's `L` is anchored to published morphometry, not to a tension. Demoting the PDE does not calibrate anything.
6. **The surrogate has no supply saturation and no sink.** It cannot represent a vessel that is itself hypoxic, or tissue whose consumption varies. Helmlinger et al. measured exactly those effects in vivo. The claim is only that on these data the extra physics costs accuracy rather than buying it.

## 7. Reproduction

```bash
python -m tmesim.validation.oxygen_benchmark
```

Deterministic (the only randomness is the seeded permutation / subsampling, `PERM_SEED = 0`, and the engine's `seed = 0`). Runtime ~10 min, dominated by the 816,769-cell reconstructions and the three transport solves.

**Files added or changed by this work:**

* `tmesim/physics/oxygen_proxy.py` — new; the surrogate and the marker route
* `tmesim/validation/oxygen_benchmark.py` — new; this analysis
* `tmesim/reconstruct/real_pipeline.py` — additive `oxygen_mode` / `oxygen_kwargs` keyword, default `'pde'` and bit-identical; the oxygen provenance block in `reconstruction_report`
* `tests/test_oxygen_proxy.py` — new
* `results/oxygen_demotion.md` — this report

**Data provenance.** Both sections are open 10x Genomics Xenium datasets, CC BY 4.0, no login or DUA: the Human Brain Cancer (GBM) FFPE IO-panel section (v2.0.0, 480 targets) and the Human Brain Glioblastoma section (v1.3.0, 319 targets).
