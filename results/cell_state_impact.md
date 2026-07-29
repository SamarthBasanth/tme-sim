# Where the cell state comes from, and what it is worth

`results/dose_schedule.md` derives ~89% of the temozolomide kill ceiling from the QUIESCENT band alone; `results/state_validation.md` shows that band is set by a tissue-mask convention and a mesh, not by tissue; `results/oxygen_demotion.md` shows the field it is banded from is beaten by a distance transform. So the largest single term in every absolute kill number this project reports is a bookkeeping choice. This report replaces it where the data allows, refuses where it does not, and measures what the difference is worth.

**The default has not been changed.** `reconstruct_real(state_mode='oxygen')` is still the shipping route and does not touch `cells.state` at all, so it is bit-identical to the pre-option pipeline (asserted by `tests/test_cell_state.py`). What is new is that a MEASURED route exists, that it declares what its panel cannot measure, and that the default's own falsification travels with it in `reconstruction_report()['cell_state']['caveat']`.

## Headline

| | |
|---|---|
| **Panels that can call proliferation per cell** | 2 of 3 — the 319-gene flagship panel is NOT one of them |
| **Panels that measure hypoxia at all** | 2 of 3 — on the 319-gene flagship section hypoxia is UNMEASURED |
| **VIABLE-fraction swing across mask x dx (flagship section)** | oxygen **unbounded (the band empties at one setting)**, vessel-distance 1.00x, markers **1.00x** |
| **Temozolomide kill across routes** | 0.0184 to 0.0315 (1.71x) |
| **Temozolomide kill CEILING across routes** | 0.5502 to 0.9561 (1.74x) |
| **Cross-drug ranking** | **INVARIANT** across all three routes |
| **Recommendation** | see §5 |

> **The three sections are THREE DIFFERENT PATIENTS** on three different panels. Nothing transfers between them and three sections is a sample size of three. None of them carries a drug-response outcome label, so nothing here says which route is RIGHT — only how much the answer depends on the choice.

## 1. Marker coverage, panel by panel

A state can only be derived from a marker the panel actually assays, and only at a detection rate that is biology rather than probe sensitivity. Detection is over TUMOUR cells (the population every kill number is about); a zero is "no probe fired", not "no transcript".

| section | panel | cells (tumour) | MKI67 | TOP2A | PCNA | VEGFA | CA9 | ADM | NDRG1 | SLC2A1 |
|---|---|---|---|---|---|---|---|---|---|---|
| `brain_319` | 319 genes | 40,887 (23,812) | 0.55% | 0.99% | 49.53% | — | — | — | — | — |
| `io_480` | 480 genes | 816,769 (540,533) | 34.06% | — | — | 65.07% | — | — | — | — |
| `gse324268_366` | 366 genes | 851,509 (515,887) | 13.34% | 25.72% | 49.18% | 45.52% | 1.03% | 10.33% | — | 12.63% |

`—` means the gene is **not on that panel**. Reading the table:

* **`brain_319`** — 10x Xenium Human Brain — Glioblastoma (v1.3.0), 319-gene brain panel.  
  Proliferation: `neighbourhood` over ['MKI67', 'TOP2A'] (union 1.419% of tumour cells). Hypoxia: `UNMEASURED` — no gene of the set is on this panel.  **UNMEASURED on this panel: HYPOXIC.**
* **`io_480`** — 10x Xenium Human Brain Cancer (GBM) FFPE (v2.0.0), 480-gene IO panel.  
  Proliferation: `per_cell` over ['MKI67'] (union 34.058% of tumour cells). Hypoxia: `per_cell` over ['VEGFA'] (union 65.07%).
* **`gse324268_366`** — GEO GSE324268 GBMtumor07, 366-gene brain + GBM add-on panel.  
  Proliferation: `per_cell` over ['MKI67', 'TOP2A'] (union 31.145% of tumour cells). Hypoxia: `per_cell` over ['VEGFA', 'CA9', 'ADM', 'SLC2A1'] (union 52.58%).

### 1.1 Why the flagship section cannot make a per-cell proliferation call

On the 319-gene panel MKI67 fires in 0.550% of tumour cells and TOP2A in 0.991%, for a union of 1.419%. Published GBM Ki-67 labelling indices are 19-40% (Fudaba et al., World Neurosurg 2021;149:e469-e480, PMID 33567368; Alkhaibary et al., Hematol Oncol Stem Cell Ther 2018;12(2):82-88, PMID 30552865; Zheng et al., Zhonghua Bing Li Xue Za Zhi 2021;50(3):229-235, PMID 33677887; Caramanti et al., Cureus 2022;14(4):e24246, PMID 35602791) — i.e. this panel sits **13x below the lower bound**. A per-cell "this cell is cycling" call at that rate reports probe sensitivity, not biology, so `state_from_markers` refuses to make one and falls back to a neighbourhood aggregation (§2.4 prices that fallback).

The same panel carries **PCNA at 49.53%**, which is inside the published Ki-67 range and looks like a rescue. It is not, and it is deliberately excluded from the default call: PCNA protein has a long half-life and the transcript is expressed in non-cycling and DNA-repairing cells. A 90x gap between PCNA and MKI67 detection **on the same cells** is a specificity difference, not a biological one — using PCNA would buy a per-cell call by measuring something other than proliferation. It is reported here and available via `include_low_specificity=True`; it is not the default.

### 1.2 The IO panel *can* — and it lands on the published index

On the 480-gene IO section MKI67 fires in 34.06% of tumour cells. That is **inside** the published 19-40% Ki-67 range, from a completely independent measurement (in-situ transcript detection vs IHC protein scoring on other cohorts). It is the single strongest piece of external corroboration in this report, and it is worth being precise about what it does and does not establish: it establishes that **this panel's MKI67 probe is not dropout-limited**, which is exactly the failure that disqualifies the 319-gene panel. It does not establish that any individual cell's call is correct. Note also that the two sections are different patients, so the comparison is panel-vs-panel, not tumour-vs-tumour.

### 1.3 The hypoxia band is measurable on two panels and its SIZE is uncalibrated on both

The 319-gene flagship panel carries no VEGFA, CA9, ADM, NDRG1 or SLC2A1, so on the section every published kill number was computed from, **hypoxia is UNMEASURED**. It has not been simulated away, and `state_from_markers` assigns no cell the HYPOXIC state there rather than manufacturing one.

Where the markers do exist, the honest limit is different and just as important. "Positive for a hypoxia-response transcript" is traceable to a measurement, but it is **not calibrated to a hypoxic fraction**: VEGFA, ADM and SLC2A1 all have real non-hypoxic expression, and VEGFA alone fires in 65.1% of IO-section tumour cells. A 60%-hypoxic tumour is not a finding, it is a threshold artefact. Of the set only **CA9** is near-restricted to severe hypoxia, and where it is on the panel it fires in 1.03% of tumour cells — a plausible perinecrotic fraction, but too sparse for a per-cell call on its own. The marker route's HYPOXIC band is therefore an ENRICHMENT label with an uncalibrated size, and this report never quotes its magnitude as a hypoxic fraction.

## 2. Three routes to a state, on the same cells

Every configuration is ONE reconstruction (the frozen `reconstruct_real` path) whose cells, solved field and vessel geometry feed all three routes, so the comparison isolates the state rule. The oxygen and vessel-distance routes both end in the UNMODIFIED `cells.agents.update_cell_states`; the marker route never touches a grid. Fractions are over TUMOUR cells.

| route | what it is | what it is traceable to |
|---|---|---|
| `oxygen` | the current default: band the SIMULATED oxygen field at `t = 0.2` / `2t` | nothing measured — a PDE with an unfitted `consumption`, gated by a tissue-mask convention |
| `vessel_distance` | band `exp(-d/L)`, `L = 93.2` um, in distance to the nearest measured endothelial cell | the MEASURED endothelial positions + the Thomlinson & Gray 150 um oxygen-diffusion distance (PMID 13304213) |
| `markers` | proliferative/quiescent from measured proliferation transcripts, hypoxic from measured hypoxia-response transcripts | MEASURED transcripts, with an explicit UNMEASURED verdict where the panel is silent |

### 2.1 State fractions across every configuration

| section | mask | dx (um) | grid | mean O2 | oxygen V/Q/H | vessel_distance V/Q/H | markers V/Q/H |
|---|---|---|---|---|---|---|---|
| `brain_319` | `geometric` | 30 | 170x271 | 0.3383 | 0.82 / 99.18 / 0.00 | 91.29 / 8.45 / 0.26 | 29.80 / 70.20 / 0.00 |
| `brain_319` | `geometric` | 20 | 255x406 | 0.3479 | 8.04 / 91.96 / 0.00 | 91.26 / 8.50 / 0.24 | 29.80 / 70.20 / 0.00 |
| `brain_319` | `geometric` | 12 | 424x675 | 0.3520 | 10.04 / 89.96 / 0.00 | 91.08 / 8.67 / 0.26 | 29.80 / 70.20 / 0.00 |
| `brain_319` | `occupancy` | 30 | 170x271 | 0.3374 | 0.82 / 99.18 / 0.00 | 91.29 / 8.45 / 0.26 | 29.80 / 70.20 / 0.00 |
| `brain_319` | `occupancy` | 20 | 255x406 | 0.3461 | 8.04 / 91.96 / 0.00 | 91.26 / 8.50 / 0.24 | 29.80 / 70.20 / 0.00 |
| `brain_319` | `occupancy` | 12 | 424x675 | 0.3050 | 0.00 / 100.00 / 0.00 | 91.08 / 8.67 / 0.26 | 29.80 / 70.20 / 0.00 |
| `io_480` | `geometric` | 30 | 363x383 | 0.2211 | 0.12 / 64.35 / 35.53 | 99.16 / 0.72 / 0.12 | 11.36 / 23.57 / 65.07 |
| `io_480` | `geometric` | 20 | 543x574 | 0.2190 | 0.01 / 63.73 / 36.26 | 99.15 / 0.73 / 0.12 | 11.36 / 23.57 / 65.07 |
| `io_480` | `geometric` | 12 | 905x956 | 0.2210 | 0.00 / 66.23 / 33.76 | 99.14 / 0.74 / 0.12 | 11.36 / 23.57 / 65.07 |
| `gse324268_366` | `geometric` | 30 | 479x353 | 0.2886 | 9.13 / 72.79 / 18.09 | 99.88 / 0.12 / 0.01 | 13.43 / 33.99 / 52.58 |
| `gse324268_366` | `geometric` | 20 | 718x529 | 0.2862 | 8.55 / 71.91 / 19.54 | 99.88 / 0.11 / 0.01 | 13.43 / 33.99 / 52.58 |
| `gse324268_366` | `geometric` | 12 | 1196x880 | 0.2870 | 8.52 / 71.93 / 19.55 | 99.88 / 0.11 / 0.01 | 13.43 / 33.99 / 52.58 |

> **These oxygen-route numbers are NOT the ones in `results/state_validation.md` §2.2, and that is expected.** That table was computed under the KDE vascular source; the shipping default is now the DISCRETE segmented-vessel source (`reconstruct.real_fields._DEFAULT_VASCULAR_MODE`, promoted in `results/vessel_promotion.md`), whose tissue-integrated supply is ~3.4x smaller, so the whole oxygen field — and therefore every band edge — sits lower. The cross-check that the path here is the frozen one is `results/dose_schedule.md` §1.3, which ran under the discrete source: it quotes 10.0% VIABLE / 90.0% QUIESCENT and a temozolomide kill ceiling of 0.5502 at `geometric`, dx = 12 um, and both are reproduced to four decimals in §3 below. The consequence for this report is that the swing figures below are re-measured under the current source rather than quoted from the older document.

### 2.2 dx-stability — the property the bar is about

Range of the tumour VIABLE fraction across the three grid spacings (30, 20, 12 um) at a FIXED tissue mask. The mesh is a numerical convenience; a state fraction that moves with it is not a property of the tissue.

| section | mask | route | VIABLE min | VIABLE max | swing | range (pp) |
|---|---|---|---|---|---|---|
| `brain_319` | `geometric` | `oxygen` | 0.8231% | 10.0412% | 12.20x | 9.22 |
| `brain_319` | `geometric` | `vessel_distance` | 91.0759% | 91.2859% | 1.00x | 0.21 |
| `brain_319` | `geometric` | `markers` | 29.8043% | 29.8043% | 1.00x | 0.00 |
| `brain_319` | `occupancy` | `oxygen` | 0.0000% | 8.0380% | unbounded (the band empties at one setting) | 8.04 |
| `brain_319` | `occupancy` | `vessel_distance` | 91.0759% | 91.2859% | 1.00x | 0.21 |
| `brain_319` | `occupancy` | `markers` | 29.8043% | 29.8043% | 1.00x | 0.00 |
| `io_480` | `geometric` | `oxygen` | 0.0011% | 0.1225% | 110.33x | 0.12 |
| `io_480` | `geometric` | `vessel_distance` | 99.1421% | 99.1579% | 1.00x | 0.02 |
| `io_480` | `geometric` | `markers` | 11.3612% | 11.3612% | 1.00x | 0.00 |
| `gse324268_366` | `geometric` | `oxygen` | 8.5228% | 9.1253% | 1.07x | 0.60 |
| `gse324268_366` | `geometric` | `vessel_distance` | 99.8763% | 99.8829% | 1.00x | 0.01 |
| `gse324268_366` | `geometric` | `markers` | 13.4336% | 13.4336% | 1.00x | 0.00 |

And across BOTH tissue-mask conventions x three dx on the flagship section — the same six-configuration grid over which `results/state_validation.md` §2.2 measured its 165x under the KDE-era vascular source, re-measured here under the current discrete source:

| route | VIABLE min | VIABLE max | swing | configurations |
|---|---|---|---|---|
| `oxygen` | 0.0000% | 10.0412% | unbounded (the band empties at one setting) | 6 |
| `vessel_distance` | 91.0759% | 91.2859% | 1.00x | 6 |
| `markers` | 29.8043% | 29.8043% | 1.00x | 6 |

**The marker route is mesh-independent by construction, and the table measures it rather than asserting it.** It is computed from measured cell coordinates and measured counts in micrometres; no grid enters. The oxygen route is not: the tissue mask gates the vascular source, so it sets the absolute oxygen level, and the mask is grid-dependent under the `occupancy` convention.

### 2.3 Cross-section stability at the production configuration (`geometric`, dx = 12 um)

| route | `brain_319` | `io_480` | `gse324268_366` | spread |
|---|---|---|---|---|
| `oxygen` | 10.0412% | 0.0011% | 8.5228% | 9045.96x |
| `vessel_distance` | 91.0759% | 99.1421% | 99.8806% | 1.10x |
| `markers` | 29.8043% | 11.3612% | 13.4336% | 2.62x |

Three different patients on three different panels, so a spread here is not by itself a fault — real tumours differ. What it does show is that no route produces a *panel-independent* number: the marker route's fractions move with which markers the panel carries and at what detection rate, which is why the coverage block in §1 has to travel with any number derived from it.

### 2.4 Pricing the one convention the marker route carries

On the flagship panel the per-cell proliferation call is inadmissible, so the label comes from a neighbourhood enrichment rule whose one free constant is the window radius. That constant is swept rather than asserted:

Section-wide tumour positive rate: 1.419% (MKI67 or TOP2A). Windows below 10 cells fall back to the section rate.

| window radius (um) | median cells/window | windows below the floor | resulting proliferative fraction |
|---|---|---|---|
| 50 | 8.0 | 67.8% | 4.83% |
| 75 | 16.0 | 10.8% | 18.45% |
| 100 **(default)** | 27.0 | 0.9% | 29.80% |
| 150 | 58.0 | 0.2% | 42.50% |
| 200 | 101.0 | 0.1% | 38.80% |

**The fraction swings 8.79x across the radius**, and that is structural rather than fixable: as R -> 0 only the positive cells themselves are enriched, so the fraction tends to the detection rate; as R -> infinity every window tends to the section rate, so the fraction tends to 0. **The ordering is measured; the level is not.** The default 100 um is chosen by a stated rule — the smallest swept radius at which >95% of windows meet the 10-cell floor on the sparsest section — and 50 um (the radius `validation.state_validation` uses for its crowding baseline) is degenerate here, with 68% of windows below the floor.

So the marker route on THIS panel buys mesh-independence and loses a window radius. That is a smaller and better-behaved convention than the oxygen route's mask x mesh swing (§2.2, unbounded here because the VIABLE band empties; 165x under the KDE-era source) — it is bounded, it is one number, and its two limits are known analytically — but it is a convention, and this report does not pretend otherwise.

### 2.5 Do the routes agree with each other?

Per-cell agreement and Cohen's kappa over tumour cells, at `geometric`, dx = 12 um.

| section | pair | agreement | Cohen's kappa |
|---|---|---|---|
| `brain_319` | `oxygen` vs `vessel_distance` | 17.81% | 0.0104 |
| `brain_319` | `oxygen` vs `markers` | 65.50% | -0.0188 |
| `brain_319` | `vessel_distance` vs `markers` | 34.26% | 0.0154 |
| `io_480` | `oxygen` vs `vessel_distance` | 0.26% | -0.0028 |
| `io_480` | `oxygen` vs `markers` | 40.71% | 0.0502 |
| `io_480` | `vessel_distance` vs `markers` | 11.46% | -0.0007 |
| `gse324268_366` | `oxygen` vs `vessel_distance` | 8.61% | 0.0001 |
| `gse324268_366` | `oxygen` vs `markers` | 34.02% | -0.0288 |
| `gse324268_366` | `vessel_distance` vs `markers` | 13.48% | 0.0002 |

Kappa near zero means the two routes are labelling the same cells no better than their marginal rates alone would predict — i.e. they agree on *how many* cells are resistant and disagree on *which*. That distinction matters downstream: a kill fraction depends almost entirely on the marginals, while any spatial claim ("this region is drug-refractory") depends on the identities, and the identities are where the routes part company.

## 3. Downstream — what the choice is worth in predicted kill

**Section.** `data_cache/xenium_gbm` — 40,887 cells (23,812 tumour), geometric tissue mask, dx = 12 um, grid 424x675. Perivascular-sink alpha = 148.61, calibrated once on doxorubicin as in every flagship run. Per-cell sensitivity is the reconstruction default (uniform 1.0), so for all three drugs the ceiling is the state labelling and nothing else — which is precisely the quantity under test.

**How the alternative arms are computed.** `sim.engine.simulate_response` re-bands the cells from its own oxygen solve and `sim/engine.py` is not modified, so each alternative arm installs its states on a private copy and re-runs the unmodified `cells.pkpd.apply_drug_kill` on **the engine's own transport result**, with the engine's own seeded generator. Valid because the cell state reaches the kill through exactly one channel — the `PDConfig` resistance multipliers — and the AUC field has no state dependence anywhere. Not assumed: the oxygen arm is recomputed through the same path and **must reproduce the engine bit-for-bit or the run aborts**. It did, for all three drugs.

### 3.1 Tumour-cell band occupancy under each route

| route | VIABLE | QUIESCENT | HYPOXIC | resistance multiplier applied |
|---|---|---|---|---|
| `oxygen` | 10.04% | 89.96% | 0.00% | x1.00 / x0.50 / x0.35 |
| `vessel_distance` | 91.08% | 8.67% | 0.26% | x1.00 / x0.50 / x0.35 |
| `markers` | 29.80% | 70.20% | 0.00% | x1.00 / x0.50 / x0.35 |

**The default route's HYPOXIC band is EMPTY on this section under the shipping mask**, so `PDConfig.hypoxia_resistance = 0.35` currently acts on no cells at all in the flagship configuration — the point `results/state_validation.md` §2.2 and `results/oxygen_demotion.md` §3.1 both made. The marker route leaves it empty too, but for a completely different and more honest reason: this panel carries no hypoxia-response gene, so hypoxia is **UNMEASURED** rather than computed-and-found-to-be-zero. Those are not the same statement and they should never be reported as if they were.

### 3.2 Predicted kill and the kill CEILING, by route

The **ceiling** is the kill fraction as exposure -> infinity: the Hill term saturates, so the per-cell probability tends to `max_kill_rate x sensitivity x resistance`. No dose can exceed it, and with uniform sensitivity it is *entirely* the state labelling. It is the cleanest measure of what this whole question is worth.

| drug | IC50 (uM) | kill (oxygen) | kill (vessel_distance) | kill (markers) | ceiling (oxygen) | ceiling (vessel_distance) | ceiling (markers) |
|---|---|---|---|---|---|---|---|
| temozolomide | 606 | 0.0184 | 0.0315 | 0.0219 | 0.5502 | 0.9561 | 0.6490 |
| cisplatin | 24.6 | 0.4362 | 0.7614 | 0.5155 | 0.5502 | 0.9550 | 0.6490 |
| doxorubicin | 0.5 | 0.2595 | 0.4593 | 0.3073 | 0.5502 | 0.9550 | 0.6490 |

**Absolute kill moves by up to 0.3252 absolute (1.77x relative) purely from the choice of state route**, with the tissue, the drug, the transport solve and the RNG all held identical. Nothing physical changed between those columns; only the label on each cell did.

### 3.3 The cross-drug ranking

| route | ranking by predicted kill |
|---|---|
| `oxygen` | cisplatin > doxorubicin > temozolomide |
| `vessel_distance` | cisplatin > doxorubicin > temozolomide |
| `markers` | cisplatin > doxorubicin > temozolomide |

**The ranking is INVARIANT and the levels are not. That is the honest and useful result.** The resistance multipliers are shared across drugs — every drug pays the same `x0.5` in the same QUIESCENT cells — so they very largely cancel out of a comparison and do not cancel at all out of a level. Concretely: a claim of the form *"cisplatin beats doxorubicin beats temozolomide on this tissue"* survives every state route tested here, while a claim of the form *"temozolomide kills 1.8% of tumour cells"* does not survive any of them.

This is the same structure `results/oxygen_demotion.md` §3.2 found for the oxygen ROUTE — ranking identical, levels moving — but roughly 23x larger: that report's kill moved at most 3.4% relative, and changing the state RULE moves it 77%. The reason is structural: swapping the oxygen field nudges a field the bands are then cut from, while swapping the state rule moves the band OCCUPANCIES themselves (§3.1: VIABLE 10.04% -> 91.08% -> 29.80%).

### 3.4 What this does NOT show

None of the three sections carries a paired drug-response outcome, so nothing here says which route's kill is correct. This is a sensitivity analysis: it measures how much the answer depends on a modelling choice. A route can be perfectly stable, fully traceable to a measurement, and still wrong about drug response.

## 4. The bar, scored

The bar, stated in advance: **state fractions must be stable across dx and across sections, and must be traceable to a measurement or explicitly declared unmeasured.**

| criterion | `oxygen` | `vessel_distance` | `markers` |
|---|---|---|---|
| stable across dx (VIABLE swing, worst of the 4 section x mask combinations; best in brackets) | unbounded (the band empties at one setting) [1.07x] | 1.00x [1.00x] | 1.00x [1.00x] |
| stable across the mask x dx grid (flagship) | unbounded (the band empties at one setting) | 1.00x | 1.00x |
| stable across the three sections (`geometric`, dx = 12 um) | 9045.96x | 1.10x | 2.62x |
| traceable to a measurement | **no** — an unfitted PDE gated by a bookkeeping mask | **partly** — measured endothelial positions + a published diffusion distance; the *band edges* are published, the *field* is a proxy | **yes, where the panel allows** — and an explicit UNMEASURED verdict where it does not |
| declares what it cannot know | no — it labels every cell regardless | no — it labels every cell regardless | **yes** — refuses HYPOXIC on the 319-gene panel, refuses a per-cell proliferation call below the detection floor |
| free constants | `consumption` (unfitted), `hypoxia_threshold`, the tissue-mask convention, the mesh | `L` (derived from a published distance), `hypoxia_threshold`, the mesh | the detection floor, and a window radius only where the per-cell call is inadmissible |

**Read the cross-section row carefully — it is the one criterion no route passes cleanly, and for three different reasons.** Three different patients on three different panels *should* differ, so a spread there is not automatically a fault. But the reasons are not equivalent: the oxygen route's 9045.96x spread includes a section where its VIABLE band is essentially empty, which is a model failure; the vessel-distance route's 1.10x is small only because it calls ~91-100% of every section VIABLE, i.e. it is stable and nearly uninformative; and the marker route's 2.62x tracks which markers each panel carries and at what detection rate, which is a property of the ASSAY and is reported alongside every number it produces (§1). So: **stability across sections is not established for any route, and the marker route is the only one whose spread is attributable to something a reader can inspect.** A fourth section on a panel matched to one of these would test it properly; three sections on three panels cannot.

## 5. RECOMMENDATION

**Adopt `state_mode='markers'` as the recommended route on any panel that clears the detection floor — the IO 480-gene and GSE324268 366-gene sections do. Keep `'oxygen'` as the shipping default until the dependent results are re-run (§5.2). And on the 319-gene flagship section, stop quoting an absolute kill fraction or a kill ceiling as a result at all: NO route clears the bar there, and the honest output of that section is a ranking, a coverage map and a delivery metric.**

In order of weight:

1. **On the two panels that can support it, the marker route clears every part of the bar that three sections are capable of testing.** It is exactly dx-stable (it never touches a grid — measured, not asserted), it is traceable to measured transcripts, and it declares its gaps. The cross-section criterion is not established for ANY route here (§4), and that is the honest gap in this recommendation rather than a point in the marker route's favour. The vessel-distance route is a genuine improvement on the PDE — measured vessel positions, a published length scale, and a flagship-grid swing of 1.00x where the oxygen route's is unbounded (the band empties at one setting) — but it is still a geometric proxy for a state, not a measurement of one, and on these sections it is nearly uninformative: it calls 91-100% of tumour cells VIABLE everywhere, so its stability is partly the stability of a near-constant.
2. **On the 319-gene flagship section NO route clears the bar, and that is the most important sentence in this report.** Hypoxia is unmeasured there; proliferation detection is 13x below the published Ki-67 lower bound, so the marker route falls back to a neighbourhood rule whose level moves 8.79x with the window radius; and the oxygen route's VIABLE band moves from 0.0000% to 10.0412% — i.e. it empties entirely — purely with bookkeeping. Every absolute kill fraction and every kill ceiling published on that section — `results/dose_schedule.md`, `results/report.md`, `results/patient_report_*.md` — inherits an unresolved uncertainty of that size. The **comparisons** in those documents are fine (§3.3); the **levels** are not, and should carry the state route and this caveat attached.
3. **Do not flip the shipping default silently.** §3.2 shows kill moves by up to 0.3252 absolute between routes, which is far too large to change under a reader who is comparing against an older document. The option exists, its provenance is recorded in `reconstruction_report()['cell_state']`, and flipping it is a deliberate act with a re-run list (§5.2).

### 5.1 What measurement would settle it

The bar is unmet on the flagship section for a data reason, not a modelling reason, so more modelling cannot fix it. What would:

* **A proliferation readout with a real detection rate on THAT section** — Ki-67 IHC/IF on a serial section, scored per cell and registered to the Xenium coordinates, or the same tumour on a whole-transcriptome platform. That converts the neighbourhood fallback into a per-cell call and removes the window radius entirely. The IO panel already demonstrates the target is reachable: MKI67 at 34.1%, inside the published range.
* **Any hypoxia-response gene on the flagship panel** — one of CA9, ADM, NDRG1, SLC2A1 would move HYPOXIC from unmeasured to measured there. CA9 is the one worth adding: it is the most hypoxia-restricted of the set, so it constrains the band's SIZE and not merely its ordering.
* **A measured oxygen tension** (pimonidazole IHC or phosphorescence-quenching pO2 on a serial section) would calibrate the hypoxic band's absolute size for every route at once — the one thing none of them currently has — and would let `OxygenConfig.consumption` be fitted instead of assumed.
* **A paired drug-response outcome.** Everything above makes the state labelling defensible; only an outcome makes the kill defensible.

### 5.2 If the default WERE flipped, this is the re-run list

Anything whose numbers depend on the per-cell band occupancy, i.e. on the `PDConfig` resistance multipliers having cells to act on. Coverage-only and AUC-only results are immune — the AUC field has no state dependence (`results/oxygen_demotion.md` §3.2 established the same immunity for the oxygen route):

| artefact | why | immune? |
|---|---|---|
| `results/dose_schedule.md` | its kill ceilings are ~89% the quiescent multiplier | **no — re-run** |
| `results/report.md`, `results/patient_report_*.md` | quote absolute kill | **no — re-run** |
| `results/ranking_robustness.md`, `results/uncertainty.md` | resample kill | **no — re-run** (their RANKINGS should survive, §3.3) |
| `results/state_validation.md`, `results/oxygen_demotion.md` | both are ABOUT the oxygen-banded states | **no — re-run** |
| `results/delivery_metrics.md` | AUC + coverage only | yes |
| `results/vessel_network_impact.md` | coverage + penetration length | mostly |

## 6. Caveats — what could make this wrong

1. **Three sections, three patients, three panels.** Nothing transfers between them, and a cross-section "spread" mixes real biological difference with panel difference. n = 3.
2. **No outcome label anywhere.** Every kill number in §3 is a prediction with nothing to score it against. This report ranks routes on stability and traceability, which are necessary and not sufficient.
3. **Hypoxia-response transcripts are not an oxygen tension.** They are the output of a programme with its own kinetics and its own non-hypoxic inducers, so the marker HYPOXIC band's SIZE is uncalibrated (§1.3). HIF1A is excluded on mechanism (oxygen sensing acts on the protein: Ivan et al., Science 2001;292(5516):464-468, PMID 11292862; Jaakkola et al., Science 2001;292(5516):468-472, PMID 11292861).
4. **Targeted-panel dropout.** A zero is "no probe fired". This is what disqualifies the per-cell call on the 319-gene panel and it also means no individual cell's marker call is trustworthy on any panel — only rates and orderings are.
5. **The neighbourhood binarisation is a convention**, spanning 8.79x across the window radius (§2.4). Where a level is needed rather than an ordering, that axis should be treated as unmeasured.
6. **Typing circularity.** `MKI67` is one of the four brain-panel TUMOR markers, so an MKI67-high cell is marginally more likely to be typed TUMOR; that inflates the marginal positive rate inside the tumour population on that section. It cannot manufacture a between-route difference (all routes share the typing) but it does inflate one input to the marker route.
7. **The one-label collapse.** Proliferation and hypoxia are two axes and `core.CellState` has one slot; the marker route resolves the clash by letting HYPOXIC win. That precedence is a schema-forced convention, stated in the returned provenance, and it is the only choice in the marker route that is not a marker threshold.
8. **The detection floor is a judgement.** 5% is the published Ki-67 lower bound divided by ~4, and it is the one number in the marker route that is neither a marker threshold nor a published anchor. Both the raw rates and the resulting call level are always returned so a reader can apply a different tolerance. The flagship panel misses it by 3.5x, so admitting it would mean setting the floor at 1.42% — 13x below the published lower bound, i.e. accepting a rate at which the call cannot be distinguished from probe noise. The IO and GSE324268 panels clear it by 6-7x, so on those two the choice of tolerance does not matter either.

## 7. Reproduction

```bash
python -m tmesim.validation.cell_state_impact
```

Deterministic — the only randomness is the engine's seeded generator (`seed = 0`). Runtime ~15 min, dominated by the reconstructions of the two ~800k-cell sections across the dx sweep.

**Files added or changed by this work:**

* `tmesim/reconstruct/cell_state.py` — new; the marker-derived state and its coverage report
* `tmesim/validation/cell_state_impact.py` — new; this analysis
* `tmesim/reconstruct/real_pipeline.py` — additive `state_mode` / `state_kwargs` keyword, default `'oxygen'` and bit-identical (it does not touch `cells.state` at all), plus the `cell_state` provenance block in `reconstruction_report`
* `tests/test_cell_state.py` — new
* `results/cell_state_impact.md` — this report

**Data provenance.** Three open single-cell spatial GBM sections: the 10x Xenium Human Brain Glioblastoma section (v1.3.0, 319 targets, CC BY 4.0), the 10x Xenium Human Brain Cancer GBM FFPE IO-panel section (v2.0.0, 480 targets, CC BY 4.0), and GEO GSE324268 GBMtumor07 (366 targets). None carries a paired drug-response outcome label.
