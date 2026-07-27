# Predicted drug response by anatomic niche region on real glioblastoma tissue

**Experiment:** `experiments/niche_response.py` (run: `./.venv/bin/python experiments/niche_response.py`)

**Figure:** `results/niche_response.png`

Everywhere else in this repo a drug gets **one** number per tumour (`kill_fraction`). That number is a spatial average — and a spatial average is exactly what hides the microenvironment. This experiment asks the question the whole model exists to answer: **where inside a real tumour does each drug actually work?**

Every spot of a real 10x Visium glioblastoma section is assigned to a real Ivy GAP anatomic zone (`margin` / `core` / `perivascular` / `necrotic`) by rank correlation against the laser-microdissected Ivy GAP per-structure FPKM signatures (Puchalski et al. 2018, *Science*). The real drug panel is then simulated on the real reconstructed microenvironment, and every readout is resolved **per zone**.

## 1. The real tissue and its real niche zones

- **Sample:** `10x_parent_gbm` — Parent_Visium_Human_Glioblastoma, Glioblastoma (primary resection, fresh-frozen).
- **Measurement:** 3,468 Visium spots × 36,601 genes; ~100 µm centre-to-centre spot pitch, 55 µm spot diameter.
- **Zone signatures:** real Ivy GAP per-anatomic-structure FPKM, reduced to a 311-gene region-discriminative panel. **283 / 311 panel genes are present in this Visium sample** — this is why Visium is the primary platform here. The Xenium GBM section shares only 41/311 panel genes, so its zone assignment rests on a much thinner correlation; it is deliberately out of scope for this experiment.
- **Assignment confidence:** mean top-1 Spearman correlation to the winning signature = 0.395; mean top-1 − top-2 margin = 0.062. These are *soft* assignments — see the caveats.

### Zone composition and local microenvironment

| zone | spots | % of tissue | grid points | cells | tumour cells | mean vascular density | mean oxygen (rel.) |
|---|---:|---:|---:|---:|---:|---:|---:|
| **margin** | 90 | 2.6% | 4,041 | 1,079 | 1,006 | 0.0890 | 0.5291 |
| **core** | 786 | 22.7% | 34,742 | 9,400 | 8,582 | 0.0778 | 0.5060 |
| **perivascular** | 964 | 27.8% | 44,398 | 11,568 | 8,048 | 0.2839 | 0.6854 |
| **necrotic** | 1,628 | 46.9% | 72,632 | 19,513 | 16,798 | 0.1036 | 0.5654 |

Marker validation of the assignment (mean per-spot log CP10k of canonical markers, by **assigned** zone — these markers are part of the signature panel, so this is an internal consistency check, not independent validation):

| zone | `PECAM1` | `VWF` | `CA9` | `VEGFA` | `MKI67` | `MBP` |
|---|---:|---:|---:|---:|---:|---:|
| **margin** | 0.326 | 0.914 | 0.017 | 1.265 | 0.000 | 0.389 |
| **core** | 0.277 | 0.755 | 0.037 | 1.533 | 0.049 | 0.160 |
| **perivascular** | 0.829 | 2.018 | 0.018 | 1.125 | 0.076 | 0.087 |
| **necrotic** | 0.171 | 0.468 | 0.224 | 2.563 | 0.047 | 0.059 |

- `PECAM1` in **perivascular** = 0.829 vs rest 0.210 — higher, as expected
- `VWF` in **perivascular** = 2.018 vs rest 0.574 — higher, as expected
- `CA9` in **necrotic** = 0.224 vs rest 0.026 — higher, as expected
- `VEGFA` in **necrotic** = 2.563 vs rest 1.306 — higher, as expected
- `MKI67` in **core** = 0.049 vs rest 0.056 — **NOT higher — inconsistent**
- `MBP` in **margin** = 0.389 vs rest 0.090 — higher, as expected

## 2. Simulation setup

- **Reconstruction:** real cells at real spot positions (41,560 cells, tumour 34,434), real vascular field (KDE of the deconvolved ENDOTHELIAL spot fraction), real ECM proxy (STROMAL fraction), real-perfusion oxygen from the unmodified physics solver.
- **Grid:** (473, 410) at dx = 12 µm (extent 4920 × 5676 µm).
- **Drug compilation (identical to the flagship panel):** `calibrate_alpha()` → α = 148.6 (fit **once** on doxorubicin) → `calibrated_drugparams(name, α)` (real GDSC/literature IC50) → `apply_barrier_to_drugparams` (saturable binding-site barrier).
- **Engine:** `simulate_response`, PDE dt = 10 s over 6 h, convection on, seed 0.

**Resolution — read this.** The PDE runs at the grid-**converged** spacing dx = 12 µm established in [`results/numerical_rigor.md`](numerical_rigor.md) §2 (`penetration_index` only settles once the ~45 µm perivascular boundary layer is resolved). **This buys numerical convergence, not biological resolution.** Every field in this experiment ultimately derives from Visium spots on a ~100 µm pitch, so the *information content* is spot-limited to ~100 µm no matter how fine the mesh is. Zone labels are replicated over each spot's ~60 µm footprint — the sub-spot structure you see in the AUC field comes from the transport physics acting on a spot-resolution source, not from sub-spot measurement. Note also that `reconstruct_real`'s Visium branch deliberately **clamps** `dx_um < 50` back to the 100 µm spot pitch for exactly this reason; the fine-grid run calls the same two public steps (`reconstruct_from_visium` + `add_oxygen`) directly, and the 100 µm pipeline grid is run as an explicit resolution control (§5).

## 3. Drug × niche zone — where each drug actually works

| drug | zone | mean tissue AUC (µM·s) | mean per-cell exposure (µM·s) | mean tumour-cell exposure (µM·s) | predicted kill of tumour cells in zone | expected kill (mean *p*) | vascular density | oxygen |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| temozolomide | margin | 1.33e+07 | 1.33e+07 | 1.34e+07 | 0.4851 | 0.4926 | 0.0890 | 0.5291 |
| temozolomide | core | 1.2e+07 | 1.21e+07 | 1.2e+07 | 0.4067 | 0.4124 | 0.0778 | 0.5060 |
| temozolomide | perivascular | 2.57e+07 | 2.6e+07 | 2.46e+07 | 0.6348 | 0.6365 | 0.2839 | 0.6854 |
| temozolomide | necrotic | 1.49e+07 | 1.51e+07 | 1.45e+07 | 0.4736 | 0.4772 | 0.1036 | 0.5654 |
| lomustine | margin | 3.65e+06 | 3.66e+06 | 3.66e+06 | 0.7684 | 0.7563 | 0.0890 | 0.5291 |
| lomustine | core | 2.98e+06 | 2.98e+06 | 2.95e+06 | 0.6192 | 0.6238 | 0.0778 | 0.5060 |
| lomustine | perivascular | 8.28e+06 | 8.44e+06 | 7.59e+06 | 0.8697 | 0.8681 | 0.2839 | 0.6854 |
| lomustine | necrotic | 3.75e+06 | 3.79e+06 | 3.56e+06 | 0.7021 | 0.7013 | 0.1036 | 0.5654 |
| cisplatin | margin | 7.61e+06 | 7.64e+06 | 7.67e+06 | 0.9493 | 0.9453 | 0.0890 | 0.5291 |
| cisplatin | core | 6.3e+06 | 6.31e+06 | 6.26e+06 | 0.8393 | 0.8432 | 0.0778 | 0.5060 |
| cisplatin | perivascular | 1.59e+07 | 1.62e+07 | 1.48e+07 | 0.9648 | 0.9651 | 0.2839 | 0.6854 |
| cisplatin | necrotic | 7.87e+06 | 7.95e+06 | 7.54e+06 | 0.9005 | 0.9013 | 0.1036 | 0.5654 |
| methotrexate | margin | 3.04e+05 | 3.04e+05 | 3.04e+05 | 0.4483 | 0.4818 | 0.0890 | 0.5291 |
| methotrexate | core | 2.49e+05 | 2.49e+05 | 2.46e+05 | 0.3846 | 0.3855 | 0.0778 | 0.5060 |
| methotrexate | perivascular | 7.44e+05 | 7.6e+05 | 6.74e+05 | 0.6887 | 0.6919 | 0.2839 | 0.6854 |
| methotrexate | necrotic | 3.17e+05 | 3.21e+05 | 2.98e+05 | 0.4469 | 0.4493 | 0.1036 | 0.5654 |
| doxorubicin | margin | 1.05e+05 | 1.05e+05 | 1.05e+05 | 0.9095 | 0.8905 | 0.0890 | 0.5291 |
| doxorubicin | core | 8.7e+04 | 8.69e+04 | 8.57e+04 | 0.7590 | 0.7637 | 0.0778 | 0.5060 |
| doxorubicin | perivascular | 2.68e+05 | 2.74e+05 | 2.41e+05 | 0.9528 | 0.9515 | 0.2839 | 0.6854 |
| doxorubicin | necrotic | 1.11e+05 | 1.12e+05 | 1.04e+05 | 0.8422 | 0.8401 | 0.1036 | 0.5654 |

Whole-tumour reference (the number the flagship panel would report, i.e. what the per-zone resolution is hiding):

| drug | MW (Da) | real IC50 (µM) | bbb | whole-tumour kill | penetration index | whole-tissue mean AUC (µM·s) | tumour hypoxic fraction |
|---|---:|---:|---:|---:|---:|---:|---:|
| temozolomide | 194 | 606 | 0.558 | 0.4949 | 0.325 | 1.56e+07 | 0.0% |
| lomustine | 234 | 50 | 0.973 | 0.7225 | 0.234 | 4.26e+06 | 0.0% |
| cisplatin | 298 | 24.6 | 0.893 | 0.9017 | 0.259 | 8.65e+06 | 0.0% |
| methotrexate | 454 | 12 | 0.115 | 0.4879 | 0.219 | 3.68e+05 | 0.0% |
| doxorubicin | 544 | 0.5 | 0.051 | 0.8493 | 0.213 | 1.3e+05 | 0.0% |

*(The per-cell kill reconstruction reproduces the engine's own `kill_fraction` to 0.00e+00 absolute — the per-zone kills are a decomposition of the engine's result, not a re-implementation.)*

### Zone ordering by delivery and by kill

| drug | zones ranked by mean tumour-cell exposure (high → low) | zones ranked by predicted kill (high → low) | expected ordering recovered? |
|---|---|---|:--:|
| temozolomide | perivascular > necrotic > margin > core | perivascular > margin > necrotic > core | no |
| lomustine | perivascular > margin > necrotic > core | perivascular > margin > necrotic > core | no |
| cisplatin | perivascular > margin > necrotic > core | perivascular > margin > necrotic > core | no |
| methotrexate | perivascular > margin > necrotic > core | perivascular > margin > necrotic > core | no |
| doxorubicin | perivascular > margin > necrotic > core | perivascular > margin > necrotic > core | no |

Textbook expectation for reference: `perivascular > margin > core > necrotic`.

## 4. Is delivery/kill really highest perivascularly and lowest in core/necrotic?

**Contrast:** `perivascular` vs `core + necrotic`, one-sided in the expected (positive) direction.

**Unit of analysis:** the **spot**, not the cell. All ~12 cells inside one Visium spot sample essentially the same drug field, so treating cells as independent would inflate *n* by ~12× and every p-value with it. n = 964 perivascular spots vs 2,414 core+necrotic spots.

**Two nulls are reported.** `p_perm` is the standard i.i.d. label shuffle — it is **anti-conservative for spatial data**, because a smooth field beats a spatially-white null almost by construction. `p_spatial` rigidly rotates and reflects the *entire* label map about the tissue centroid, preserving patch size, class balance and spatial autocorrelation, and only randomising the label map's registration to the drug field. **`p_spatial` is the p-value to quote.** Both use the add-one correction, so the floor is 1/10,001 and 1/1,001 respectively.

### 4a. Delivery (per-spot drug exposure AUC)

| drug | mean AUC perivascular | mean AUC core+necrotic | ratio | Δ | Cohen's *d* | Cliff's δ | *p*<sub>perm</sub> | *p*<sub>spatial</sub> |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| temozolomide | 2.6e+07 | 1.41e+07 | 1.84× | 1.19e+07 | +1.211 | +0.609 | 0.0001 | 0.004 |
| lomustine | 8.43e+06 | 3.53e+06 | 2.39× | 4.9e+06 | +1.361 | +0.676 | 0.0001 | 0.003 |
| cisplatin | 1.62e+07 | 7.42e+06 | 2.18× | 8.75e+06 | +1.316 | +0.661 | 0.0001 | 0.003 |
| methotrexate | 7.59e+05 | 2.98e+05 | 2.55× | 4.61e+05 | +1.389 | +0.684 | 0.0001 | 0.003 |
| doxorubicin | 2.73e+05 | 1.04e+05 | 2.62× | 1.69e+05 | +1.401 | +0.687 | 0.0001 | 0.003 |

### 4b. Kill (per-spot mean predicted kill probability of tumour cells)

| drug | mean *p*<sub>kill</sub> perivascular | mean *p*<sub>kill</sub> core+necrotic | ratio | Δ | Cohen's *d* | Cliff's δ | *p*<sub>perm</sub> | *p*<sub>spatial</sub> |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| temozolomide | 0.6478 | 0.4621 | 1.40× | 0.1857 | +1.117 | +0.604 | 0.0001 | 0.007 |
| lomustine | 0.8784 | 0.6835 | 1.29× | 0.1949 | +0.977 | +0.667 | 0.0001 | 0.006 |
| cisplatin | 0.9678 | 0.8850 | 1.09× | 0.08276 | +0.559 | +0.653 | 0.0001 | 0.006 |
| methotrexate | 0.7138 | 0.4398 | 1.62× | 0.274 | +1.283 | +0.675 | 0.0001 | 0.006 |
| doxorubicin | 0.9565 | 0.8207 | 1.17× | 0.1358 | +0.718 | +0.676 | 0.0001 | 0.006 |

### 4c. Pairwise zone contrasts (delivery, first drug in the panel as reference)

| contrast | ratio | Cohen's *d* | Cliff's δ | *p*<sub>spatial</sub> |
|---|---:|---:|---:|---:|
| perivascular vs necrotic | 1.72× | +1.061 | +0.565 | 0.013 |
| perivascular vs core | 2.15× | +1.341 | +0.701 | 0.008 |
| perivascular vs margin | 1.95× | +1.085 | +0.671 | 0.16 |
| margin vs necrotic | 0.89× | -0.195 | +0.003 | 0.32 |
| core vs necrotic | 0.80× | -0.347 | -0.243 | 0.7 |

## 5. Resolution control (dx = 100 µm, `reconstruct_real`'s native Visium grid)

The same panel, the same zone labels, on the assembler's own spot-pitch grid:

| drug | perivasc/(core+necrotic) AUC ratio @12 µm | @100 µm | Cliff's δ @12 µm | @100 µm | whole-tumour kill @12 µm | @100 µm |
|---|---:|---:|---:|---:|---:|---:|
| temozolomide | 1.84× | 1.85× | +0.609 | +0.605 | 0.4949 | 0.5105 |
| lomustine | 2.39× | 2.38× | +0.676 | +0.672 | 0.7225 | 0.7389 |
| cisplatin | 2.18× | 2.18× | +0.661 | +0.656 | 0.9017 | 0.9176 |
| methotrexate | 2.55× | 2.53× | +0.684 | +0.679 | 0.4879 | 0.5011 |
| doxorubicin | 2.62× | 2.60× | +0.687 | +0.681 | 0.8493 | 0.8639 |

The perivascular/(core+necrotic) contrast has the **same direction at both resolutions** for every drug (12 µm ratios 1.84–2.62× vs 100 µm ratios 1.85–2.60×). Refining the mesh sharpens the near-vessel gradient the PDE can represent, but it adds no measurement information: the source field is spot-limited at both resolutions.

## 6. Honest interpretation — does microenvironment structure produce regionally different drug response?

**Yes — the expected perivascular gradient appears, and it is large.** For every drug in the panel, mean per-spot exposure is higher in the perivascular zone than in the core+necrotic zones, by 1.84–2.62×, with Cliff's δ from +0.61 to +0.69 and a spatially-structured permutation p ≤ 0.004. The microenvironment does **not** deliver drug uniformly: the same systemic dose produces materially different local exposure depending on which anatomic zone a tumour cell sits in.

**The full four-zone ordering is only partly recovered.** Observed orderings by mean tumour-cell exposure: temozolomide: perivascular > necrotic > margin > core; lomustine: perivascular > margin > necrotic > core; cisplatin: perivascular > margin > necrotic > core; methotrexate: perivascular > margin > necrotic > core; doxorubicin: perivascular > margin > necrotic > core — against the textbook `perivascular > margin > core > necrotic`. The perivascular-vs-core/necrotic contrast is the robust part; the finer ordering is not.

**But a large part of this is circular, and that has to be said plainly.** The `perivascular` zone label is driven by endothelial/angiogenic transcripts (`PECAM1`, `VWF`, `CD34`, `CLDN5`, `FLT1`, …), and the model's drug **source term is proportional to the vascular density field**, which is itself a KDE of the deconvolved ENDOTHELIAL spot fraction — i.e. of the *same* endothelial marker expression at the *same* spots. So 'drug arrives where the endothelial signal is' is close to a restatement of the model's input, not an independent discovery. The perivascular result should be read as a **consistency check that the pipeline propagates real vascular structure end-to-end**, not as evidence that the model predicts delivery. The genuinely non-circular comparisons are the ones that do not involve the perivascular label — `core` vs `necrotic` and `margin` vs `necrotic` (necrotic is called from HIF targets `CA9`/`VEGFA`/`ADM`/`NDRG1`, core from `EGFR`/`SOX2`/`OLIG2`/`MKI67`, neither of which feeds the vascular field) — see §4c.

**The non-circular half of the prediction FAILS, and this is the most informative result here.** Stripping out every contrast that involves the perivascular label leaves: **margin vs necrotic**: 0.89× (δ=+0.00, p_spatial=0.32), **core vs necrotic**: 0.80× (δ=-0.24, p_spatial=0.7). The textbook expectation is that the necrotic/perinecrotic zone is the worst-perfused, worst-delivered compartment in the tumour — it is not. On this section the `necrotic` zone receives **more** predicted drug than the `core` (1.25× more, Cliff's δ=+0.24), and the difference is not significant against the spatial null (p_spatial = 0.7). Mechanically this follows straight from the reconstructed fields: mean vascular density is 0.1036 in the necrotic zone vs 0.0778 in the core, and mean oxygen 0.5654 vs 0.5060. Biologically that is defensible — perinecrotic HIF-high tissue (`CA9`/`VEGFA`/`ADM`) is precisely the tissue that *drives* microvascular proliferation, so it sits next to vessels — but it means **this model does not reproduce the classic 'necrotic core is a drug sanctuary' picture**, and the honest reading is that once you remove the vascular-marker circularity, only the perivascular-vs-everything contrast survives, and the finer zone structure does not.

**Exposure differences convert into kill differences.** Predicted kill of the tumour cells inside the perivascular zone exceeds that in core+necrotic for every drug, by 1.09–1.62×. This is the spatially explicit version of the whole claim: a single whole-tumour `kill_fraction` averages over cells whose predicted fates differ by that factor. Note however that the **kill contrast is much smaller than the exposure contrast** (exposure 1.84–2.62× vs kill 1.09–1.62×): the Hill dose-response compresses it, because at these predicted exposures most drugs sit high on the sigmoid where extra drug buys little extra kill. A delivery gradient only matters when the tissue straddles the drug's IC50: here the largest kill contrast is **methotrexate** (1.62×) and the smallest is **cisplatin** (1.09×), which tracks where each drug's achieved exposure sits on its own dose-response curve rather than how steep its delivery gradient is.

**Hypoxia is not the mechanism here.** On this section at dx = 12 µm the predicted tumour hypoxic fraction is 0.0%, so the hypoxia-resistance multiplier in the PD model is essentially inactive and the regional kill differences are driven almost entirely by **delivered exposure**, not by hypoxic drug resistance. The Visium-derived vascular field is a spot-resolution KDE and is consequently much smoother than the Xenium single-cell field, which is why this section reads as well-oxygenated. Do not transfer the hypoxia numbers from the Xenium panel to this one.

**What would make this a real result rather than a coherent one:** paired region-resolved drug measurement (e.g. MALDI-MSI or LC-MS of a treated resection, zone-annotated), which would let the per-zone AUC column be scored instead of merely reported. Nothing in this experiment is validated against a measured outcome.

## 7. Caveats — what this does NOT show

- **No outcome labels anywhere.** This Visium section is an untreated primary resection. There is no measured drug response, and no measured per-zone drug concentration, for any of the five drugs. Every kill and exposure number is a **prediction of the physics on real geometry**, not a validated outcome.
- **Zone labels are signature-correlation assignments, not pathology.** Each spot is given the zone whose Ivy GAP signature it rank-correlates with best. The mean winning correlation is only 0.395 and the mean margin over the runner-up is 0.062 — these are **soft, low-margin assignments**, not annotated histology. A pathologist's zone map of this section does not exist in the open data, so the labelling cannot be scored. Treat zone membership as a graded score that has been hardened into a class.
- **Circularity in the perivascular arm** (see §6): the perivascular label and the model's vascular source term are both functions of endothelial marker expression at the same spots. The perivascular-vs-rest contrast is therefore partly a restatement of an input.
- **Spot-limited resolution.** Visium spots are 55 µm across on a ~100 µm pitch. Real perivascular gradients in GBM are ~50 µm wide — at or below the measurement scale. The 12 µm mesh is a numerical convergence requirement and creates no sub-spot information; zone labels are replicated over each spot's footprint. Any claim about structure finer than ~100 µm is unsupported by the data, whatever the figure looks like.
- **Deconvolution stands between the counts and the cells.** Cell types (and hence the vascular and ECM fields) come from reference-based deconvolution of spot expression, not from observed cells. Errors there propagate into every per-zone number.
- **Necrotic zones contain fewer viable cells than the model assumes.** The reconstruction places a fixed bag of cells per spot apportioned by deconvolved fractions; it does not down-weight the truly acellular fraction of a necrotic zone. Per-zone kill *fractions* in the necrotic zone are therefore computed over a cell population that is probably too viable and too numerous.
- **One section, one patient.** n = 1 tissue. Nothing here establishes that these zone contrasts generalise; the cross-tissue question is not addressed.
- **Statistical caveats.** The spot is treated as the unit of analysis, which is better than the cell but still not independent — neighbouring spots are strongly autocorrelated. That is why `p_spatial` (rotation null) is quoted rather than `p_perm` (i.i.d. null); the rotation null is itself approximate, since nearest-neighbour label transfer after rotation slightly perturbs class counts near the tissue boundary. The five drugs are not independent tests — they share one tissue, one vascular field and one α calibration, so their p-values are effectively one test repeated with different drug parameters.
- **A single calibration constant sets the delivery scale.** α is fit once on a doxorubicin penetration measurement; the absolute AUC of every zone rides on that one anchor. Zone *ratios* are far more trustworthy than zone *absolutes*.
- **Not every canonical marker separates the way it should.** `MKI67` in core = 0.049 vs rest 0.056 — NOT higher — inconsistent. A zone whose own canonical marker does not enrich in it is a zone whose assignment should be treated with extra suspicion — the correlation-based call is being carried by the other ~280 panel genes, not by the marker the zone is named for.
- **Some zones are small.** `margin` has only 90 spots. Their per-zone means are noisy and their pairwise contrasts are underpowered — visible in §4c, where contrasts involving these zones do not reach significance against the spatial null.
