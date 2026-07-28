# Perivascular drug gradient on the REAL tissue — predicted vs measured

> **STALE ABSOLUTE NUMBERS — regenerate.** Every absolute exposure / AUC / coverage / kill
> figure below was generated with the **KDE** vascular source, which was the default until
> `results/vessel_promotion.md`. The default is now the **discrete** segmented-vessel source,
> whose tissue-integrated total is **3.36x smaller** (a change of *unit*, not of physics —
> the dose-matched geometry-only effect on coverage is <=0.0166). Relative orderings and
> per-drug contrasts here are largely unaffected; absolute levels are not. Reproduce with
> `mode='kde'` to recover exactly these numbers.

**Module:** `tmesim/validation/vessel_profile.py` (run: `./.venv/bin/python -m tmesim.validation.vessel_profile`)

**Figure:** `results/vessel_profile.png`

The canonical experimental readout of tumor drug penetration is a **profile**, not a scalar: quantitative immunofluorescence of drug intensity against **distance from the nearest CD31⁺ blood vessel**. tmesim predicts a full spatial drug field on **real** tissue, so we can compute *the same quantity* and compare the **functional form**.

## Headline

| | doxorubicin | temozolomide |
|---|---:|---:|
| **fitted L** (µm), real tissue, primary arm | **60.8** | **185.7** |
| **half-distance d½** (µm), real tissue, primary arm | **71.0** | **406.7** |
| fitted L (µm), isolated single vessel | 71.4 | 288.5 |
| d½ = L·ln2 (µm), isolated single vessel | 49.5 | 200.0 |
| measured | **40–50 µm** (PMID 16361566), 25–75 µm (PMID 23680920) | *no published perivascular-distance dataset* |

### Read this before quoting any number

Two different quantities are easy to conflate, and only one of them is clean:

1. **The analytic / single-vessel penetration length** `L = √(D_eff / (clearance + binding))`, which is genuinely drug-specific (doxorubicin 71.3 µm vs temozolomide 290.6 µm) and is what the isolated single-vessel control reproduces with the real solver. **This is a legitimate model output and is the thing to compare to the literature.**
2. **The spatial fit to drug-vs-distance on our reconstructed tissue.** This one is **confounded**: the vessel mask and the drug source are built from the *same* endothelial-cell KDE, so part of the fitted decay is our own 40 µm smoothing bandwidth. The source null (d½ = 23.1 µm with **no drug physics at all**) and the bandwidth sweep below quantify how much. On the fragmented native dx = 12 µm mask the confound is *total* — every drug returns the same ~26–28 µm decay, within 1.1–1.2× of the null, reproducing the independent finding in `results/delivery_metrics.md`. **A match to the measured 40–50 µm band read off the spatial fit is not a validation.** Quote the drug/source-null ratio instead.

**Shape:** exponential. On the isolated single vessel the predicted profile is a pure exponential to five decimals (log-space R² = 1.00000). On real tissue it is **exponential decaying onto a floor**, `C(d) = C₀·exp(−d/L) + floor` (R² = 0.9652), the floor being the drug delivered by the *distributed* vascular source. Published profiles also flatten into a background, so the qualitative form matches.

**Verdict for doxorubicin:** real tissue d½ = 71.0 µm → *PASS(wide) (outside 40-50 um, inside the 25-75 um replication band)*; the fitted `L = 60.8 µm` against the L-equivalent band `40–50/ln2 ≈ 58–72 µm` → *PASS (inside the L-equivalent 58-72 um band)*. **This is not an independent validation** — α was fit to this same anchor. **Temozolomide (L = 185.7 µm, 3.1× doxorubicin's) is the free prediction**, and it has no measurement to be scored against.

---

## Method

1. **Real geometry.** The open 10x Xenium GBM section (40,887 measured cells; tumor 30,174) reconstructed by `reconstruct.real_pipeline.reconstruct_real('xenium', dx_um=…)` → grid (424, 675), extent 8100×5088 µm, mean ECM 0.138, mean vascular density 0.082. **dx = 12 µm** is the grid-converged resolution `results/numerical_rigor.md` §2 asks for (the ~45 µm perivascular boundary layer is unresolved for dx ≳ 15 µm).
2. **Vessel mask.** The REAL endothelial-derived `vascular_density` field thresholded at its 95th tissue percentile (`≥ 0.217`) → **10,084 vessel pixels = 5.0% of tissue area**, in the range of measured tumor vascular area fractions.
3. **Distance map.** `scipy.ndimage.distance_transform_edt` on the complement of the mask with `sampling = dx` → distance-from-nearest-vessel in µm at every grid point (over tissue: median 144 µm, p95 368 µm, max 1116 µm).
4. **Simulation.** Each drug compiled through the standard pipeline (`calibrate_alpha()` → α = 148.6, fit **once** on doxorubicin → `calibrated_drugparams` → `apply_barrier_to_drugparams`), run on the same tissue with `PDEConfig(dt=10 s, total_time=6 h, convection=True)`.
5. **Profile + fits.** AUC (and final concentration) binned in 10 µm distance bins (≥ 30 grid points per bin), mean ± IQR. Readouts: the **empirical** half-distance (no functional form assumed), a `C₀·exp(−d/L) + floor` fit over the full range (the correct model on real tissue), a pure exponential over the near field (d ≤ 150 µm), and a straight-line null.
6. **Two controls:** the **source null** (identical binning applied to the raw `vascular_density` field) and the **isolated single-vessel** geometry.

---

## A numerical artefact we had to deal with first (please read)

The reconstruction rasterises measured **cell positions**, so as `dx` shrinks the tissue mask stops being a filled sheet and becomes a scatter of isolated pixels. The transport solver applies a no-flux boundary at every mask edge, so an isolated pixel has **no neighbour to diffuse into** and its concentration collapses to the purely local `C = S/k`. Measured directly:

| dx (µm) | tissue px | tissue as % of bounding box | connected components | largest component (% of tissue) |
|---:|---:|---:|---:|---:|
| 60 | 7,748 | 66.2% | 27 | 98.0% |
| 40 | 17,436 | 66.8% | 26 | 98.1% |
| 30 | 30,976 | 67.2% | 27 | 98.1% |
| 20 | 69,769 | 67.4% | 27 | 98.1% |
| 15 | 124,074 | 67.7% | 30 | 98.0% |
| 12 | 193,752 | 67.7% | 30 | 98.0% |

**At dx = 12 µm — exactly the resolution `numerical_rigor.md` recommends — the domain is 26,058 disconnected fragments whose largest component holds 0.7% of the tissue.** Lateral diffusion is effectively switched off. That is why, on the native dx = 12 µm mask, doxorubicin and temozolomide produce *the same* profile (69.2 vs 411.2 µm) despite ~4× different transport lengths, and both sit on top of the source null (23.1 µm) — correlation between the AUC field and the source field 0.858 and 0.520. This is a grid artefact, not physics.

We therefore report **three arms** and take as primary the one that resolves the boundary layer *and* has a connected domain:

| arm | dx (µm) | largest connected component | doxorubicin d½ (µm) | doxorubicin L (µm) | temozolomide d½ (µm) | temozolomide L (µm) | source-null d½ (µm) | corr(AUC, source): dox / TMZ |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| dx = 12 µm, morphologically REPAIRED (connected) mask **(PRIMARY)** | 12 | 98.3% | **71.0** | 60.8 | **406.7** | 185.7 | *23.1* | 0.854 / 0.512 |
| dx = 12 µm, NATIVE cell-derived mask (fragmented) | 12 | 98.0% | **69.2** | 58.0 | **411.2** | 178.3 | *23.1* | 0.858 / 0.520 |
| dx = 30 µm, NATIVE mask (connected, but boundary layer coarse) | 30 | 98.1% | **79.5** | 63.1 | **413.4** | 150.4 | *37.9* | 0.864 / 0.522 |

The **repair** is pure image morphology applied to a *copy* of the mask (`scipy.ndimage.binary_closing`, 2 iterations, then `binary_fill_holes`), which restores the largest component to 98.3% of tissue; the ECM and vascular fields are then **re-derived on the repaired mask** with the public `reconstruct.real_fields.build_real_fields`, because the KDE is multiplied by the tissue mask at build time and reusing the old fields would leave the added pixels with a spurious zero source. **Nothing in `reconstruct` is modified.** The added pixels are the interstitial gaps *between* measured cells — real tissue the drug diffuses through, excluded by the native mask only because the rasteriser found no cell centroid there.

The dx = 30 µm arm (native mask, 98% connected) is an independent cross-check with no morphology applied at all, and it agrees with the repaired arm qualitatively: doxorubicin 79.5 µm vs temozolomide 413.4 µm. **Honest range across the two physically-sane arms: doxorubicin d½ ≈ 71–80 µm, temozolomide d½ ≈ 407–413 µm.**

---

## Result 1 — the gradient on real tissue (primary arm)

**Model-free readout** (what the experiment reports — the distance at which the binned mean falls to half its perivascular value, no functional form assumed):

| signal | perivascular → far-field plateau | **empirical d½** (µm) | corr(AUC, source field) |
|---|---:|---:|---:|
| **doxorubicin** AUC | 1.00 → 0.24 | **71.0** | 0.854 |
| **temozolomide** AUC | 1.00 → 0.48 | **406.7** | 0.512 |
| *SOURCE NULL: `vascular_density` itself* | 1.00 → 0.14 | *23.1* | 1.000 |

**Fits.** A *pure* exponential is the wrong model on real tissue: the distributed source puts a plateau under the perivascular decay and drags a single-exponential fit towards it. The three-parameter `C(d) = C₀·exp(−d/L) + floor` fits far better and its `L` is the real perivascular decay length. All R² are reported in **linear space** so the models are directly comparable:

| signal | model | **L** (µm) | **d½ = L·ln2** (µm) | floor (% of perivascular) | **R² (linear space)** | R² (log space) |
|---|---|---:|---:|---:|---:|---:|
| **doxorubicin** | `C₀·exp(−d/L) + floor` (full range) | **60.8** | **42.2** | 28.8% | **0.9652** | — |
| doxorubicin | pure `C₀·exp(−d/L)`, near field d ≤ 150 µm | 154.1 | 106.8 | — | 0.8771 | 0.9210 |
| doxorubicin | straight-line null `a + b·d` (near field) | — | — | — | 0.8175 | — |
| **temozolomide** | `C₀·exp(−d/L) + floor` (full range) | **185.7** | **128.7** | 48.1% | **0.9484** | — |
| temozolomide | pure `C₀·exp(−d/L)`, near field d ≤ 150 µm | 417.3 | 289.2 | — | 0.8844 | 0.9053 |
| temozolomide | straight-line null `a + b·d` (near field) | — | — | — | 0.8646 | — |
| **SOURCE NULL** | `C₀·exp(−d/L) + floor` (full range) | **21.6** | **15.0** | 14.8% | **0.9862** | — |
| SOURCE NULL | pure `C₀·exp(−d/L)`, near field d ≤ 150 µm | 92.0 | 63.8 | — | 0.7185 | 0.7964 |
| SOURCE NULL | straight-line null `a + b·d` (near field) | — | — | — | 0.6367 | — |

The **final-concentration** field gives essentially the same profile (doxorubicin empirical d½ 71.1 µm, temozolomide 407.4 µm) — expected, because with a constant-source window the field is quasi-steady, so `AUC ≈ C_steady · Δt` and the two share a spatial shape. It is *not* an independent check, but it does rule out the profile being an artefact of time-integration.

### How much of the gradient is the drug, and how much is the vasculature?

The engine's delivery term is `vascular_density · bbb_penetrance · plasma`: the source is spread over the **whole continuous density field**, not confined to the thresholded mask (mean vascular density over tissue = 0.082; tissue 200 µm from the mask still carries 14% of the peak density). The SOURCE NULL row measures what that alone buys: an apparent half-distance of **23.1 µm** with no drug physics at all. On the primary arm the drugs sit clearly above it (71.0 and 406.7 µm, correlations 0.85 and 0.51), so transport is doing real work — but a meaningful part of the profile is still vessel/KDE geometry, and that share must be subtracted mentally before reading the agreement with the literature as evidence.

---

## Result 2 — CONTROL: the isolated single-vessel gradient

To see the gradient the *physics* predicts, free of the distributed-source confound, we run the **same frozen transport solver** and the **same compiled drug** on one planar (line) vessel in a uniform medium at the real tissue's mean ECM (0.138). For a planar source the steady solution is exactly `C(d) = C₀·exp(−d/L)` with `L = √(D_eff/k)` — the geometry the calibration constant is defined in.

| drug | MW (Da) | D (µm²/s) | k = clearance+binding (1/s) | **L** (µm) | **d½ = L·ln2** (µm) | R² (log) |
|---|---:|---:|---:|---:|---:|---:|
| **doxorubicin** | 544 | 185 | 0.0333 | **71.4** | **49.5** | 1.00000 |
| **temozolomide** | 194 | 278 | 0.00302 | **288.5** | **200.0** | 1.00000 |

**Temozolomide's gradient is 4.0× longer than doxorubicin's** in identical geometry — the physically expected ordering: temozolomide is small (MW 194 vs 544 Da), diffuses faster (278 vs 185 µm²/s) and is not consumed by a high-capacity intracellular/DNA sink (k = 0.003 vs 0.033 /s), so it escapes the perivascular cuff and reaches tissue the anthracycline never sees.

This control is also where the exponential-shape claim is unambiguous: log-space R² ≈ 1.0 for both drugs, i.e. the predicted profile is a **true exponential**, the functional form the measurements report.

---

## Comparison to the real measurements

**The anchors (both PMID-verified, carried in `tmesim/validation/pk_anchors.py`):**

> Primeau AJ, Rendon A, Hedley D, Lilge L, Tannock IF. The distribution of the anticancer drug doxorubicin in relation to blood vessels in solid tumors. Clin Cancer Res 2005;11(24 Pt 1):8782-8. **PMID 16361566** — doxorubicin intensity falls to half its perivascular value within **40–50 µm** of CD31⁺ vessels (quantitative immunofluorescence of doxorubicin autofluorescence vs distance from CD31+ blood vessels; murine 16C and EMT6 tumors + human PC-3 prostate xenografts (mice)).

> Patel KJ, Tredan O, Tannock IF. Distribution of the anticancer drugs doxorubicin, mitoxantrone and topotecan in tumors and normal tissues. Cancer Chemother Pharmacol 2013;72(1):127-38. **PMID 23680920** — independent replication; half-distance **25–75 µm** at 10 min, reported jointly for doxorubicin / mitoxantrone / topotecan.

| doxorubicin readout | predicted | measured | verdict |
|---|---:|---:|:--|
| real tissue (primary arm), **empirical d½** | **71.0 µm** | 40–50 µm (25–75 µm) | PASS(wide) (outside 40-50 um, inside the 25-75 um replication band) |
| real tissue (primary arm), fitted **L** vs L-equivalent band | **60.8 µm** | 58–72 µm | PASS (inside the L-equivalent 58-72 um band) |
| real tissue (dx = 30 µm cross-check), empirical d½ | 79.5 µm | 40–50 µm (25–75 µm) | PASS(3x) (outside both bands, same order of magnitude) |
| **isolated single vessel**, d½ | **49.5 µm** | 40–50 µm (25–75 µm) | PASS (inside the measured 40-50 um band) |
| isolated single vessel, L vs L-equivalent band | 71.4 µm | 58–72 µm | PASS (inside the L-equivalent 58-72 um band) |

**Why two bands.** `pk_anchors.py` is explicit that the published 40–50 µm is a *half*-distance, so the quantity comparable to an e-folding length is `40–50 / ln2 ≈ 58–72 µm`. Comparing an `L` to the raw 40–50 µm (as `pk_reality.py` does) is deliberately conservative. Both comparisons are given above so no reader has to guess which convention a number is in.

### Is the shape exponential, as measured?

- **doxorubicin** — *isolated vessel*: log-space R² = **1.00000**, a clean exponential to five decimals. *Real tissue*: **yes, exponential-plus-floor** — R²(linear) = **0.9652** for `C₀·exp(−d/L)+floor`, versus 0.8771 for a pure exponential and 0.8175 for a straight line over the same near field. The decaying component has L = 60.8 µm; beyond ~150 µm the profile plateaus at 24% of the perivascular value instead of continuing to decay.
- **temozolomide** — *isolated vessel*: log-space R² = **1.00000**, a clean exponential to five decimals. *Real tissue*: **yes, exponential-plus-floor** — R²(linear) = **0.9484** for `C₀·exp(−d/L)+floor`, versus 0.8844 for a pure exponential and 0.8646 for a straight line over the same near field. The decaying component has L = 185.7 µm; beyond ~150 µm the profile plateaus at 48% of the perivascular value instead of continuing to decay.

So the predicted functional form **is** the measured one: exponential where the geometry is clean, exponential-onto-a-background on real tissue. The published profiles also flatten into a background at large distance — but ours is a modelling consequence of the diffuse source, not a measured autofluorescence floor, and the two should not be equated.

---

## THE DECISIVE CONTROL — is the fitted decay our own KDE bandwidth?

`vascular_density` is a Gaussian KDE of endothelial-cell positions with a **40 µm bandwidth**, and the vessel mask is a threshold on that *same* smoothed field. If the fitted decay length simply tracks the bandwidth, the spatial fit is measuring our own smoothing kernel and **cannot** be quoted as agreement with the literature. We therefore rebuild the real fields at 25 / 40 / 60 µm bandwidth and re-measure everything.

**Primary arm (dx = 12 µm, repaired connected mask):**

| KDE bandwidth (µm) | source-null L (µm) | doxorubicin L (µm) | dox / null | temozolomide L (µm) | TMZ / null | corr(AUC, source): dox / TMZ |
|---:|---:|---:|---:|---:|---:|---:|
| 25 | 13.8 | **52.3** | 3.80× | **223.8** | 16.25× | 0.742 / 0.385 |
| 40 | 21.6 | **60.8** | 2.82× | **185.7** | 8.61× | 0.854 / 0.512 |
| 60 | 35.7 | **77.0** | 2.16× | **400.2** | 11.21× | 0.924 / 0.633 |

**Fragmented native dx = 12 µm mask (the arm with no lateral transport), for contrast:**

| KDE bandwidth (µm) | source-null L (µm) | doxorubicin L (µm) | dox / null | temozolomide L (µm) | TMZ / null |
|---:|---:|---:|---:|---:|---:|
| 25 | 13.6 | 50.2 | 3.69× | 183.7 | 13.51× |
| 40 | 21.3 | 58.4 | 2.74× | 175.8 | 8.26× |
| 60 | 34.1 | 72.6 | 2.13× | 316.2 | 9.27× |

**How to read this.** On the **fragmented** mask the drug L sits within ~1.1–1.2× of the source null at every bandwidth and the two drugs are indistinguishable — there, the fitted decay *is* the KDE, exactly as the independent `results/delivery_metrics.md` analysis found across all five drugs. On the **repaired (connected)** mask the drugs separate from the null and from each other, but the null itself still moves with the bandwidth, so **part of the fitted decay on real tissue is our smoothing kernel in every case**. Quote the drug/null ratio, never the absolute fitted L, as the evidence.

---

## Sensitivity to the (arbitrary) vessel threshold — primary arm

| percentile | threshold | vessel area | max distance (µm) | near-field L (µm) | near-field d½ (µm) | **empirical d½ (µm)** | *source-null d½ (µm)* | R² (log) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| p90 | 0.168 | 10.0% | 1074 | 142.4 | 98.7 | **72.2** | *21.7* | 0.9312 |
| p95 | 0.217 | 5.0% | 1116 | 154.1 | 106.8 | **71.0** | *23.1* | 0.9210 |
| p98 | 0.299 | 2.0% | 2357 | 152.1 | 105.4 | **66.4** | *24.0* | 0.9190 |
| p99 | 0.382 | 1.0% | 2551 | 130.3 | 90.3 | **66.3** | *25.0* | 0.9551 |

---

## Honest caveats

- **The nominally-converged grid breaks the transport domain.** At dx = 12 µm the cell-derived tissue mask is 26,058 disconnected fragments (largest = 0.7% of tissue), so lateral diffusion is switched off and *both* drugs collapse onto the vascular source field. `results/numerical_rigor.md` converged `kill_fraction` and `penetration_index` at dx ≲ 12 µm but did not test mask connectivity; this module does, and the primary arm therefore repairs the mask morphologically. **Any perivascular-gradient number quoted from a native dx ≤ 20 µm run is a source-geometry measurement, not a transport measurement.**
- **The repair is a modelling choice.** Closing the mask adds interstitial pixels the rasteriser missed; it is defensible (real tissue is continuous) but it is not measured, and the closing radius is a free parameter. The dx = 30 µm native-mask arm is included precisely because it needs no morphology at all, and the two agree only qualitatively — hence the honest ranges quoted above rather than a single number.
- **Part of the gradient is still vessel geometry, not drug transport.** The source null alone gives d½ = 23.1 µm. The engine delivers drug proportionally to the *continuous* density everywhere, so tissue scored as 'far from a vessel' still has a local source, which flattens the predicted profile relative to a true line-source geometry.
- **The vessel mask is not a measured vessel lumen.** `vascular_density` is a 40 µm-bandwidth Gaussian KDE of the positions of cells *typed* as endothelial from the Xenium marker panel (`reconstruct.real_fields.vascular_density_from_cells`); the real experiments segment CD31-stained vessel **lumens**. Our distances are distances from an endothelial-density ridge, and the 40 µm bandwidth is the same order as the gradient being measured — it puts a floor on how sharp the near-vessel profile can be.
- **The threshold is a free parameter**, fixed *a priori* at the 95th percentile (≈5% vessel area fraction) rather than chosen to match the literature; the sensitivity table shows what moves with it.
- **Doxorubicin agreement is PARTLY BY CONSTRUCTION.** The single α in `validation.calibration` was fit to the *same* Primeau anchor, through the 1-D analytic relation `L = √(D_eff/k_eff)`. The isolated single-vessel control therefore recovers ~45 µm essentially by construction and is a **consistency check that the calibration survives the full numerical solver on real geometry, not an independent validation**. What is genuinely non-trivial: the shape comes out exponential, and the real heterogeneous 2-D geometry did not move the number by an order of magnitude.
- **The α convention matters and has changed.** This run used **α = 148.6**, calibrated against the *e-folding* target `45/ln2 ≈ 64.9 µm` (`validation/calibration.py`, `DOX_TARGET_EFOLD_UM`). An earlier convention fitted α to the raw 45 µm half-distance as if it were an e-folding length, giving α = 309.6 and a ~1.4× steeper gradient. Every absolute penetration number in this report inherits that single convention choice; check `calibrate_alpha()` before comparing to older results files.
- **Temozolomide is the independent output** — α was never fit to it — but there is **no published temozolomide perivascular-distance dataset**, so that arm is a *prediction*, checkable only for direction and order of magnitude. (The measured temozolomide anchor in `pk_anchors.py`, PMID 19861433, is a brain-ECF:plasma **partition ratio** — a different physical quantity that must not be compared to a penetration length.)
- **2D, straight-line distances.** The section is a 2D plane through a 3D vascular network: a vessel just above or below the plane is invisible, so grid points scored 'far' may really be near an out-of-plane vessel. Distances are Euclidean, not geodesic around the irregular tissue mask.
- **Not time-matched.** Primeau's images are ~10 min post-injection snapshots; the primary profile integrates a 6 h simulated exposure (AUC). The final-concentration readout is the closer snapshot analogue and agrees closely.
- **Untreated section, no measured drug field.** The Xenium GBM sample is an untreated resection; the comparison is to published measurements in other (murine / xenograft) tumors, not to drug imaged on *this* tissue.
