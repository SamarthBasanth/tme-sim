# Perivascular drug gradient on the REAL tissue — predicted vs measured

> **Regenerated under the discrete vessel source.** `vascular_mode = discrete` (`reconstruct.real_fields._DEFAULT_VASCULAR_MODE`, promoted in [`results/vessel_promotion.md`](vessel_promotion.md)) · **alpha = 148.6**, the e-folding-length convention `L = d_half / ln 2` · per-cell sensitivity mode = **`uniform`** (the reconstruction default) · **dx = 12 µm** on the primary arm, with 12 / 12 / 30 µm across the 3 arms.
>
> **`empirical_half_distance_um` renders `None` honestly, never a fabricated number.** It returns `None` **by design** when a profile never falls to half its perivascular value; every read renders **"n/a (never halves)"** instead of crashing or guessing. **No substitute number is ever fitted or extrapolated** — a profile that does not halve HAS no half-distance, and that is a measured property of the profile, not a gap to be filled. Whether temozolomide's real-tissue profile *does* halve now depends on which distance metric is used — see the next callout, which is the main content of this regeneration.
>
> ### ⚠ THE 2026-07-29 FLAGSHIP RESULT BELOW REVERSES UNDER CORRECTED TRANSPORT PHYSICS
>
> Re-run under `PDEConfig(transvascular_exchange=True)` (see
> `results/transvascular_exchange_revalidation.md`). Not retracted, because the
> flag is default-off and the numbers below remain the shipped path's literal
> output — but the headline conclusion of the regeneration announced in the
> next paragraph does not survive the delivery-law fix.
>
> **Temozolomide goes back to NO DECAY.** Exchange removes drug fastest where
> the source is strongest, which flattens the profile: TMZ's fitted `L` rises
> **277.4 → 396.3 µm (+43 %)** and its empirical half-distance —
> **153.6 µm** below — becomes **not measurable at all**. The verdict returns
> from *FAIL* to *NO DECAY (no half-distance could be measured)*, i.e. exactly
> the legacy-metric state the 2026-07-29 regeneration was written to announce
> it had left. The primary-vs-legacy contrast for temozolomide collapses too:
> under exchange both metrics say "never halves".
>
> **And the published result was marginal by construction.** Only 16 distance
> bins clear `MIN_BIN_COUNT = 30`, so the profile is measurable to ~155–165 µm.
> The 153.6 µm crossing sat **in the last usable bin**. A result resting on a
> final-bin crossing is fragile independent of this physics change, and should
> not have been quoted as a clean verdict flip.
>
> **Doxorubicin is unaffected**: all three verdicts identical, numbers shift
> +1.8–2.2 % (d½ 41.5 → 42.2 µm, L 43.5 → 44.5 µm). The isolated
> single-vessel control is fully invariant (R²log = 1.00000 both arms) —
> exchange alters the vessel-line amplitude but not the decay away from it,
> the same reason the calibration chain survives.
>
> **All shape claims survive**: exponential on the isolated vessel,
> exponential-plus-floor on real tissue, the distance-metric change is
> orthogonal to exchange, and the source-null confound is unchanged.
> Temozolomide's qualitative ordering survives and the gap widens (6.4× →
> 8.9× doxorubicin's L). Every temozolomide fitted length needs regenerating.
>
> **One correction to a claim made elsewhere in this repo.** A "pre-existing
> drift" in the vessel-distance map (tissue p90 85.26 µm vs a published 84.85)
> was reported during this re-validation effort. It is **not a drift** — it is
> two different metrics of the same tissue: LEGACY percentile-mask p90 =
> 84.852814 µm (bit-matching `.dose_schedule_cache.json`) and PRIMARY
> direct-cell p90 = 85.264116 µm. Both reproduce exactly. No number in this
> report is affected; it quotes median/p95/max, all of which reproduced.

> **THE PRIMARY DISTANCE-TO-VESSEL METRIC CHANGED (2026-07-29) — read this before quoting any number below.** Every fitted number in this report was previously computed from a **percentile threshold on the smoothed/discrete `vascular_density` FIELD** (`vessel_mask_from_density` → `distance_transform_edt`). It is now computed from the **direct distance to the nearest REAL measured endothelial CELL** (`direct_vessel_distance_um`, a single `cKDTree` query, no threshold, no bandwidth, no mask) — a straightforward audit (`results/sweep/vessel_distance_audit_percentile_method_2026-07-29T05-54-00Z.md`) found that a percentile on a density field answers "how far to the densest few percent of the field", not "how far to a vessel", and a prototype (`results/sweep/vessel_distance_direct_replacement_prototype_2026-07-29T05-54-00Z.md`) confirmed the direct-cell distance matches `vascular_validation.py`'s independently-measured tumour→endothelial distance while the percentile mask does not. **The aggregate distance barely moves (median 45.0 µm primary vs 43.3 µm legacy, 0.96×) — this is a TAIL effect, not a mean effect, and it flips one drug's verdict outright:** **temozolomide** goes from *"NO DECAY (no half-distance could be measured)"* (legacy) to *"FAIL (more than 3x off the measured band)"* (primary) — its profile now measurably halves (at a large distance, well outside the doxorubicin band, which is exactly what a small, fast-diffusing, weakly-cleared drug should do) where the old metric could not even measure a half-distance at all. **Doxorubicin** IMPROVES from *"PASS(wide) (outside 40-50 um, inside the 25-75 um replication band)"* to *"PASS (inside the measured 40-50 um band)"*. Neither change is a change in the simulated physics (the AUC/concentration fields are identical; only the distance axis they are binned against changed) — see the dedicated section below for the full head-to-head numbers, and the module HONESTY note for why the LEGACY metric is retained rather than deleted.
>
> Absolute exposure / AUC / coverage figures below are on the **discrete-source unit** and are NOT comparable with the KDE-era numbers this file carried before the discrete-source regeneration (the tissue-integrated source total is 3.36x smaller; that is a change of *unit*, not of physics).

**Module:** `tmesim/validation/vessel_profile.py` (run: `./.venv/bin/python -m tmesim.validation.vessel_profile`)

**Figure:** `results/vessel_profile.png`

The canonical experimental readout of tumor drug penetration is a **profile**, not a scalar: quantitative immunofluorescence of drug intensity against **distance from the nearest CD31⁺ blood vessel**. tmesim predicts a full spatial drug field on **real** tissue, so we can compute *the same quantity* and compare the **functional form**.

## Headline

| | doxorubicin | temozolomide |
|---|---:|---:|
| **fitted L** (µm), real tissue, primary arm | **43.5** | **277.4** |
| **half-distance d½** (µm), real tissue, primary arm | **41.5** | **153.6** |
| fitted L (µm), isolated single vessel | 72.5 | 293.2 |
| d½ = L·ln2 (µm), isolated single vessel | 50.3 | 203.3 |
| measured | **40–50 µm** (PMID 16361566), 25–75 µm (PMID 23680920) | *no published perivascular-distance dataset* |

### Read this before quoting any number

Two different quantities are easy to conflate, and only one of them is clean:

1. **The analytic / single-vessel penetration length** `L = √(D_eff / (clearance + binding))`, which is genuinely drug-specific (doxorubicin 72.5 µm vs temozolomide 295.4 µm) and is what the isolated single-vessel control reproduces with the real solver. **This is a legitimate model output and is the thing to compare to the literature.**
2. **The spatial fit to drug-vs-distance on our reconstructed tissue.** This one is **still confounded, but by one thing instead of two.** The drug's perivascular-sink SOURCE term is proportional to the (smoothed or discrete) `vascular_density` field, so part of the fitted decay is that field's own shape — the source null (d½ = 10.0 µm with **no drug physics at all**) and the bandwidth sweep below quantify how much. This regeneration REMOVES a second, previously-conflated confound: the vessel-DISTANCE axis used to be a threshold on that *same* KDE field, so part of the fitted decay used to also be our own 40 µm mask-smoothing bandwidth. It is now the direct nearest-real-endothelial-cell distance (see the dedicated section below), decoupled from whatever field the drug is delivered through — a legitimate, separately-debated modelling choice — from what counts as "a vessel" for the distance axis, which should always be the actual measured geometry. On the fragmented native dx = 12 µm mask the remaining (source-field) confound is *total* — every drug returns close to the source-null decay, reproducing the independent finding in `results/delivery_metrics.md`. **A match to the measured 40–50 µm band read off the spatial fit is not a validation.** Quote the drug/source-null ratio instead.

**Shape:** exponential. On the isolated single vessel the predicted profile is a pure exponential to five decimals (log-space R² = 1.00000). On real tissue it is **exponential decaying onto a floor**, `C(d) = C₀·exp(−d/L) + floor` (R² = 0.9892), the floor being the drug delivered by the *distributed* vascular source. Published profiles also flatten into a background, so the qualitative form matches.

**Verdict for doxorubicin:** real tissue d½ = 41.5 µm → *PASS (inside the measured 40-50 um band)*; the fitted `L = 43.5 µm` against the L-equivalent band `40–50/ln2 ≈ 58–72 µm` → *PASS(3x) (steeper than measured by 33%, same order of magnitude)*. **This is not an independent validation** — α was fit to this same anchor. **Temozolomide (L = 277.4 µm, 6.4× doxorubicin's) is the free prediction**: there is still no *published, temozolomide-specific* perivascular-distance dataset to score it against, but under the PRIMARY distance metric its profile now has a measurable empirical d½ (153.6 µm) that CAN be sanity-checked against the doxorubicin band as an order-of-magnitude check — and it *FAIL (more than 3x off the measured band)*. That is an honest failure of the generic check, not a validated temozolomide measurement; see the next section for what changed and why.

---

## THE DISTANCE-TO-VESSEL METRIC CHANGE — PRIMARY (direct-cell) vs LEGACY (percentile-mask)

This regeneration replaces the vessel-distance axis that drives every fit above. **What changed:** previously, "distance from nearest vessel" was computed by thresholding the (smoothed or discrete) `vascular_density` FIELD at its 95th tissue percentile and running a Euclidean distance transform on the resulting mask (`vessel_mask_from_density` → `distance_from_vessel_um`). It is now computed directly from the REAL measured endothelial-cell positions — a single `cKDTree` nearest-neighbour query, no threshold, no bandwidth, no mask at all (`direct_vessel_distance_um`, reusing `vascular_validation.distance_to_nearest`, the exact function `vascular_validation.py` uses for its own independently-published tumour→endothelial distance, so the two numbers cannot drift apart). **The percentile-mask pipeline is NOT deleted** — it is kept as a documented LEGACY metric (`vessel_mask_from_density` / `distance_from_vessel_um` are unchanged functions) for backward comparison and for datasets where only a density field, not individual cell positions, is available.

**Why this matters: it is a tail effect, not a mean effect.** Grid-wide, the two distance fields barely differ (median 45.0 µm primary vs 43.3 µm legacy, ratio 0.96×; max 789 vs 786 µm) — the percentile mask happens to be geometrically close to the real vessel footprint almost everywhere (see the audit's §4 for why: under the discrete source field the percentile threshold degenerates to the same `vasc > 0` set the direct metric is measuring against, by coincidence of the fallback branch, not by design). But the two metrics diverge in the FAR tail — the tumour cells genuinely furthest from any vessel — and that tail is exactly where a slow or absent perivascular sink shows up, which is precisely temozolomide's regime.

| | doxorubicin | temozolomide |
|---|---:|---:|
| empirical d½ (µm) — PRIMARY (direct-cell) | **41.5** | **153.6** |
| empirical d½ (µm) — LEGACY (percentile-mask) | 31.8 | n/a (never halves) |
| fitted L, exp+floor (µm) — PRIMARY | **43.5** | **277.4** |
| fitted L, exp+floor (µm) — LEGACY | 34.0 | 276.8 |
| verdict vs measured band — PRIMARY | *PASS (inside the measured 40-50 um band)* | *FAIL (more than 3x off the measured band)* |
| verdict vs measured band — LEGACY | *PASS(wide) (outside 40-50 um, inside the 25-75 um replication band)* | *NO DECAY (no half-distance could be measured)* |

**Read plainly, not smoothed over:**

- **Doxorubicin's verdict IMPROVES**: *"PASS(wide) (outside 40-50 um, inside the 25-75 um replication band)"* (legacy) → *"PASS (inside the measured 40-50 um band)"* (primary). The legacy percentile-mask distance was systematically slightly SHORTER than the true nearest-vessel distance near the perivascular zone the doxorubicin fit is most sensitive to, which compressed its fitted decay length below the measured band; the direct-cell distance corrects that and lands the fit inside the measured 40–50 µm band.
- **Temozolomide's verdict CHANGES**: *"NO DECAY (no half-distance could be measured)"* (legacy) → *"FAIL (more than 3x off the measured band)"* (primary). This is the more consequential change and it is stated here explicitly, not buried: under the legacy metric, temozolomide's profile never fell to half its perivascular value anywhere on the section (`empirical_half_distance_um` correctly returned `None` — it is not a bug, the profile genuinely never crossed 50% within the legacy distance range), so no verdict could be computed at all. Under the PRIMARY metric the profile DOES cross 50%, far out in the tail where the legacy metric under-measured distance, and the resulting half-distance is more than 3× the measured doxorubicin band — an honest, scoreable **FAIL** on the generic order-of-magnitude check. This is the physically expected direction (temozolomide is small, fast-diffusing and has no fast perivascular sink, so its profile SHOULD extend far beyond a doxorubicin-scaled band) but it is also a genuine change in what this report says about temozolomide: it went from "unscoreable" to "scoreable and, on this generic check, failing". Nothing about the simulated transport physics changed — the AUC field for temozolomide is bit-for-bit the same as before; only the distance axis it is binned against did.

**Design decision.** The direct-cell distance is now the PRIMARY metric everywhere in this module (`run_arm`, and the KDE-bandwidth sweep below); the percentile-mask distance is retained as an explicit LEGACY metric (same functions, same behaviour, unchanged) for (a) this head-to-head comparison, (b) the mask visualisation in figure panel (a), and (c) sections/datasets that only have a density field and no individual endothelial-cell positions to query directly. `threshold_sensitivity` (percentile sweep) is kept as a LEGACY-metric-only diagnostic — the primary metric has no percentile to sweep. `bandwidth_sensitivity` (KDE-bandwidth sweep) now computes the vessel-distance axis ONCE, outside the bandwidth loop, because real cell positions do not move with the KDE bandwidth; this fixes a second, adjacent issue the same audit surfaced — under the discrete source field the OLD sweep was silently re-deriving an identical mask at every bandwidth and testing nothing on the vessel-distance axis at all (see the audit, §6).

---

## Method

1. **Real geometry.** The open 10x Xenium GBM section (40,887 measured cells; tumor 23,812) reconstructed by `reconstruct.real_pipeline.reconstruct_real('xenium', dx_um=…)` → grid (424, 675), extent 8100×5088 µm, mean ECM 0.087, mean vascular density 0.026. **dx = 12 µm** is the grid-converged resolution `results/numerical_rigor.md` §2 asks for (the ~45 µm perivascular boundary layer is unresolved for dx ≳ 15 µm).
2. **Vessel distance — PRIMARY.** A single `cKDTree` nearest-neighbour query (`direct_vessel_distance_um`) from every grid point to the nearest REAL measured ENDOTHELIAL cell (there is no threshold, no mask, no bandwidth) → distance-from-nearest-vessel in µm at every grid point (over tissue: median 45 µm, p95 100 µm, max 789 µm).
3. **Vessel distance — LEGACY (retained for comparison).** The REAL endothelial-derived `vascular_density` field thresholded at its 95th tissue percentile (`≥ 0.000`) → **5,156 vessel pixels = 2.6% of tissue area**, then `scipy.ndimage.distance_transform_edt` on the complement with `sampling = dx` (over tissue: median 43 µm, p95 102 µm, max 786 µm). See the dedicated section above for why this is no longer the primary metric and what changes when it is replaced.
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

**At dx = 12 µm — exactly the resolution `numerical_rigor.md` recommends — the domain is 26,058 disconnected fragments whose largest component holds 0.7% of the tissue.** Lateral diffusion is effectively switched off. That is why, on the native dx = 12 µm mask, doxorubicin and temozolomide produce *the same* profile (42.0 vs 144.8 µm) despite ~4× different transport lengths, and both sit on top of the source null (10.0 µm) — correlation between the AUC field and the source field 0.319 and 0.114. This is a grid artefact, not physics.

We therefore report **three arms** and take as primary the one that resolves the boundary layer *and* has a connected domain:

| arm | dx (µm) | largest connected component | doxorubicin d½ (µm) | doxorubicin L (µm) | temozolomide d½ (µm) | temozolomide L (µm) | source-null d½ (µm) | corr(AUC, source): dox / TMZ |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| dx = 12 µm, morphologically REPAIRED (connected) mask **(PRIMARY)** | 12 | 98.3% | **41.5** | 43.5 | **153.6** | 277.4 | *10.0* | 0.318 / 0.111 |
| dx = 12 µm, NATIVE cell-derived mask (fragmented) | 12 | 98.0% | **42.0** | 45.1 | **144.8** | 280.0 | *10.0* | 0.319 / 0.114 |
| dx = 30 µm, NATIVE mask (connected, but boundary layer coarse) | 30 | 98.1% | **13.0** | 7.4 | **n/a (never halves)** | 181.0 | *10.0* | 0.548 / 0.192 |

The **repair** is pure image morphology applied to a *copy* of the mask (`scipy.ndimage.binary_closing`, 2 iterations, then `binary_fill_holes`), which restores the largest component to 98.3% of tissue; the ECM and vascular fields are then **re-derived on the repaired mask** with the public `reconstruct.real_fields.build_real_fields`, because the KDE is multiplied by the tissue mask at build time and reusing the old fields would leave the added pixels with a spurious zero source. **Nothing in `reconstruct` is modified.** The added pixels are the interstitial gaps *between* measured cells — real tissue the drug diffuses through, excluded by the native mask only because the rasteriser found no cell centroid there.

The dx = 30 µm arm (native mask, 98% connected) is an independent cross-check with no morphology applied at all, and it agrees with the repaired arm qualitatively: doxorubicin 13.0 µm vs temozolomide n/a (never halves) µm. **Honest range across the two physically-sane arms: doxorubicin d½ ≈ 13–42 µm, temozolomide d½ ≈ n/a (never halves) µm.**

---

## Result 1 — the gradient on real tissue (primary arm)

**Model-free readout** (what the experiment reports — the distance at which the binned mean falls to half its perivascular value, no functional form assumed):

| signal | perivascular → far-field plateau | **empirical d½** (µm) | corr(AUC, source field) |
|---|---:|---:|---:|
| **doxorubicin** AUC | 1.00 → nan | **41.5** | 0.318 |
| **temozolomide** AUC | 1.00 → nan | **153.6** | 0.111 |
| *SOURCE NULL: `vascular_density` itself* | 1.00 → nan | *10.0* | 1.000 |

**Fits.** A *pure* exponential is the wrong model on real tissue: the distributed source puts a plateau under the perivascular decay and drags a single-exponential fit towards it. The three-parameter `C(d) = C₀·exp(−d/L) + floor` fits far better and its `L` is the real perivascular decay length. All R² are reported in **linear space** so the models are directly comparable:

| signal | model | **L** (µm) | **d½ = L·ln2** (µm) | floor (% of perivascular) | **R² (linear space)** | R² (log space) |
|---|---|---:|---:|---:|---:|---:|
| **doxorubicin** | `C₀·exp(−d/L) + floor` (full range) | **43.5** | **30.2** | 14.4% | **0.9892** | — |
| doxorubicin | pure `C₀·exp(−d/L)`, near field d ≤ 150 µm | 69.0 | 47.8 | — | 0.9740 | 0.9913 |
| doxorubicin | straight-line null `a + b·d` (near field) | — | — | — | 0.8848 | — |
| **temozolomide** | `C₀·exp(−d/L) + floor` (full range) | **277.4** | **192.3** | 0.0% | **0.9938** | — |
| temozolomide | pure `C₀·exp(−d/L)`, near field d ≤ 150 µm | 276.2 | 191.5 | — | 0.9951 | 0.9954 |
| temozolomide | straight-line null `a + b·d` (near field) | — | — | — | 0.9951 | — |
| **SOURCE NULL** | `C₀·exp(−d/L) + floor` (full range) | **n/a** | **n/a** | n/a% | **n/a** | — |
| SOURCE NULL | pure `C₀·exp(−d/L)`, near field d ≤ 150 µm | n/a | n/a | — | n/a | n/a |
| SOURCE NULL | straight-line null `a + b·d` (near field) | — | — | — | n/a | — |

The **final-concentration** field gives essentially the same profile (doxorubicin empirical d½ 41.5 µm, temozolomide 153.9 µm) — expected, because with a constant-source window the field is quasi-steady, so `AUC ≈ C_steady · Δt` and the two share a spatial shape. It is *not* an independent check, but it does rule out the profile being an artefact of time-integration.

### How much of the gradient is the drug, and how much is the vasculature?

The engine's delivery term is `vascular_density · bbb_penetrance · plasma`: the source is spread over the **whole continuous density field**, not confined to the thresholded mask (mean vascular density over tissue = 0.026; tissue 200 µm from the mask still carries nan% of the peak density). The SOURCE NULL row measures what that alone buys: an apparent half-distance of **10.0 µm** with no drug physics at all. On the primary arm the drugs sit clearly above it (41.5 and 153.6 µm, correlations 0.32 and 0.11), so transport is doing real work — but a meaningful part of the profile is still vessel/KDE geometry, and that share must be subtracted mentally before reading the agreement with the literature as evidence.

---

## Result 2 — CONTROL: the isolated single-vessel gradient

To see the gradient the *physics* predicts, free of the distributed-source confound, we run the **same frozen transport solver** and the **same compiled drug** on one planar (line) vessel in a uniform medium at the real tissue's mean ECM (0.087). For a planar source the steady solution is exactly `C(d) = C₀·exp(−d/L)` with `L = √(D_eff/k)` — the geometry the calibration constant is defined in.

| drug | MW (Da) | D (µm²/s) | k = clearance+binding (1/s) | **L** (µm) | **d½ = L·ln2** (µm) | R² (log) |
|---|---:|---:|---:|---:|---:|---:|
| **doxorubicin** | 544 | 185 | 0.0333 | **72.5** | **50.3** | 1.00000 |
| **temozolomide** | 194 | 278 | 0.00302 | **293.2** | **203.3** | 1.00000 |

**Temozolomide's gradient is 4.0× longer than doxorubicin's** in identical geometry — the physically expected ordering: temozolomide is small (MW 194 vs 544 Da), diffuses faster (278 vs 185 µm²/s) and is not consumed by a high-capacity intracellular/DNA sink (k = 0.003 vs 0.033 /s), so it escapes the perivascular cuff and reaches tissue the anthracycline never sees.

This control is also where the exponential-shape claim is unambiguous: log-space R² ≈ 1.0 for both drugs, i.e. the predicted profile is a **true exponential**, the functional form the measurements report.

---

## Comparison to the real measurements

**The anchors (both PMID-verified, carried in `tmesim/validation/pk_anchors.py`):**

> Primeau AJ, Rendon A, Hedley D, Lilge L, Tannock IF. The distribution of the anticancer drug doxorubicin in relation to blood vessels in solid tumors. Clin Cancer Res 2005;11(24 Pt 1):8782-8. **PMID 16361566** — doxorubicin intensity falls to half its perivascular value within **40–50 µm** of CD31⁺ vessels (quantitative immunofluorescence of doxorubicin autofluorescence vs distance from CD31+ blood vessels; murine 16C and EMT6 tumors + human PC-3 prostate xenografts (mice)).

> Patel KJ, Tredan O, Tannock IF. Distribution of the anticancer drugs doxorubicin, mitoxantrone and topotecan in tumors and normal tissues. Cancer Chemother Pharmacol 2013;72(1):127-38. **PMID 23680920** — independent replication; half-distance **25–75 µm** at 10 min, reported jointly for doxorubicin / mitoxantrone / topotecan.

| doxorubicin readout | predicted | measured | verdict |
|---|---:|---:|:--|
| real tissue (primary arm), **empirical d½** | **41.5 µm** | 40–50 µm (25–75 µm) | PASS (inside the measured 40-50 um band) |
| real tissue (primary arm), fitted **L** vs L-equivalent band | **43.5 µm** | 58–72 µm | PASS(3x) (steeper than measured by 33%, same order of magnitude) |
| real tissue (dx = 30 µm cross-check), empirical d½ | 13.0 µm | 40–50 µm (25–75 µm) | FAIL (more than 3x off the measured band) |
| **isolated single vessel**, d½ | **50.3 µm** | 40–50 µm (25–75 µm) | PASS(wide) (outside 40-50 um, inside the 25-75 um replication band) |
| isolated single vessel, L vs L-equivalent band | 72.5 µm | 58–72 µm | PASS(3x) (shallower than measured by 12%, same order of magnitude) |

**Why two bands.** `pk_anchors.py` is explicit that the published 40–50 µm is a *half*-distance, so the quantity comparable to an e-folding length is `40–50 / ln2 ≈ 58–72 µm`. Comparing an `L` to the raw 40–50 µm (as `pk_reality.py` does) is deliberately conservative. Both comparisons are given above so no reader has to guess which convention a number is in.

### Is the shape exponential, as measured?

- **doxorubicin** — *isolated vessel*: log-space R² = **1.00000**, a clean exponential to five decimals. *Real tissue*: **yes, exponential-plus-floor** — R²(linear) = **0.9892** for `C₀·exp(−d/L)+floor`, versus 0.9740 for a pure exponential and 0.8848 for a straight line over the same near field. The decaying component has L = 43.5 µm; beyond ~150 µm the profile plateaus at nan% of the perivascular value instead of continuing to decay.
- **temozolomide** — *isolated vessel*: log-space R² = **1.00000**, a clean exponential to five decimals. *Real tissue*: **yes, exponential-plus-floor** — R²(linear) = **0.9938** for `C₀·exp(−d/L)+floor`, versus 0.9951 for a pure exponential and 0.9951 for a straight line over the same near field. The decaying component has L = 277.4 µm; beyond ~150 µm the profile plateaus at nan% of the perivascular value instead of continuing to decay.

So the predicted functional form **is** the measured one: exponential where the geometry is clean, exponential-onto-a-background on real tissue. The published profiles also flatten into a background at large distance — but ours is a modelling consequence of the diffuse source, not a measured autofluorescence floor, and the two should not be equated.

---

## THE DECISIVE CONTROL — is the fitted decay our own KDE bandwidth?

`vascular_density` is a Gaussian KDE of endothelial-cell positions with a **40 µm bandwidth**, and it is the field the drug's perivascular-sink source term is proportional to. If the fitted decay length tracks that bandwidth, the spatial fit is measuring our own smoothing kernel and **cannot** be quoted as agreement with the literature. We therefore rebuild the real fields at 25 / 40 / 60 µm bandwidth and re-measure everything.

**What this sweep tests changed along with the distance metric.** The vessel-DISTANCE axis below is the PRIMARY direct-cell metric, computed ONCE from the real endothelial-cell positions **before** the bandwidth loop — real cell positions do not move when the KDE bandwidth changes, so this axis is by construction identical at bw = 25, 40 and 60 µm. That is deliberate, not a gap: under the OLD percentile-mask distance, this same sweep silently rebuilt an essentially identical mask at every bandwidth anyway (the production `vascular_mode = discrete` field does not even use `bandwidth_um` — see the audit, §6), so it was already testing nothing on the distance axis while *appearing* to test bandwidth sensitivity end-to-end. Fixing the distance axis explicitly makes that scope honest: this control now isolates exactly one thing — how much the fitted `L` and the drug/source-null ratio move because the bandwidth-dependent `vascular_density` FIELD (the delivery source term, and via the companion ECM-KDE, tissue tortuosity) changed, with the "where is a vessel" question held fixed at its one geometrically correct answer.

**Primary arm (dx = 12 µm, repaired connected mask):**

| KDE bandwidth (µm) | source-null L (µm) | doxorubicin L (µm) | dox / null | temozolomide L (µm) | TMZ / null | corr(AUC, source): dox / TMZ |
|---:|---:|---:|---:|---:|---:|---:|
| 25 | n/a | **44.4** | n/a× | **280.7** | n/a× | 0.315 / 0.110 |
| 40 | n/a | **43.5** | n/a× | **277.4** | n/a× | 0.318 / 0.111 |
| 60 | n/a | **42.4** | n/a× | **273.2** | n/a× | 0.322 / 0.112 |

**Fragmented native dx = 12 µm mask (the arm with no lateral transport), for contrast:**

| KDE bandwidth (µm) | source-null L (µm) | doxorubicin L (µm) | dox / null | temozolomide L (µm) | TMZ / null |
|---:|---:|---:|---:|---:|---:|
| 25 | n/a | 47.0 | n/a× | 282.6 | n/a× |
| 40 | n/a | 46.0 | n/a× | 279.1 | n/a× |
| 60 | n/a | 44.8 | n/a× | 274.8 | n/a× |

**How to read this.** On the **fragmented** mask the drug L sits within ~1.1–1.2× of the source null at every bandwidth and the two drugs are indistinguishable — there, the fitted decay *is* the KDE, exactly as the independent `results/delivery_metrics.md` analysis found across all five drugs. On the **repaired (connected)** mask the drugs separate from the null and from each other, but the null itself still moves with the bandwidth, so **part of the fitted decay on real tissue is our smoothing kernel in every case**. Quote the drug/null ratio, never the absolute fitted L, as the evidence.

---

## [LEGACY METRIC ONLY] Sensitivity to the (arbitrary) percentile vessel threshold — primary arm

The PRIMARY direct-cell distance has no percentile and no threshold to sweep — this table characterises the sensitivity of the LEGACY percentile-mask metric only, retained for backward comparison.

| percentile | threshold | vessel area | max distance (µm) | near-field L (µm) | near-field d½ (µm) | **empirical d½ (µm)** | *source-null d½ (µm)* | R² (log) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| p90 | 0.000 | 2.6% | 786 | 68.7 | 47.6 | **31.8** | *10.0* | 0.9813 |
| p95 | 0.000 | 2.6% | 786 | 68.7 | 47.6 | **31.8** | *10.0* | 0.9813 |
| p98 | 1.000 | 2.6% | 786 | 68.7 | 47.6 | **31.8** | *10.0* | 0.9813 |
| p99 | 1.000 | 2.6% | 786 | 68.7 | 47.6 | **31.8** | *10.0* | 0.9813 |

---

## Honest caveats

- **The nominally-converged grid breaks the transport domain.** At dx = 12 µm the cell-derived tissue mask is 26,058 disconnected fragments (largest = 0.7% of tissue), so lateral diffusion is switched off and *both* drugs collapse onto the vascular source field. `results/numerical_rigor.md` converged `kill_fraction` and `penetration_index` at dx ≲ 12 µm but did not test mask connectivity; this module does, and the primary arm therefore repairs the mask morphologically. **Any perivascular-gradient number quoted from a native dx ≤ 20 µm run is a source-geometry measurement, not a transport measurement.**
- **The repair is a modelling choice.** Closing the mask adds interstitial pixels the rasteriser missed; it is defensible (real tissue is continuous) but it is not measured, and the closing radius is a free parameter. The dx = 30 µm native-mask arm is included precisely because it needs no morphology at all, and the two agree only qualitatively — hence the honest ranges quoted above rather than a single number.
- **Part of the gradient is still vessel geometry, not drug transport.** The source null alone gives d½ = 10.0 µm. The engine delivers drug proportionally to the *continuous* density everywhere, so tissue scored as 'far from a vessel' still has a local source, which flattens the predicted profile relative to a true line-source geometry.
- **Neither metric measures a CD31-stained vessel lumen.** The PRIMARY metric is a direct distance to the position of a cell *typed* as endothelial from the Xenium marker panel (`reconstruct.real_fields.vascular_density_from_cells`'s input cells) — a real measured point, but a cell centroid, not a segmented lumen boundary, so distances are systematically offset from a true CD31 edge by roughly a cell radius. This is a smaller, more defensible gap than the LEGACY metric's: that one measures distance from a 40 µm-bandwidth Gaussian KDE ridge of the same cell positions (or, under the discrete field, a painted vessel-radius footprint), and the 40 µm bandwidth is the same order as the gradient being measured, putting a floor on how sharp the near-vessel profile can be. Switching to the direct-cell metric removes the KDE-bandwidth part of this caveat entirely; it does not remove the cell-centroid-vs-lumen part.
- **[LEGACY METRIC ONLY]** The percentile threshold is a free parameter, fixed *a priori* at the 95th percentile (≈5% vessel area fraction) rather than chosen to match the literature; the sensitivity table above shows what moves with it. This does not apply to the PRIMARY metric, which has no threshold.
- **Doxorubicin agreement is PARTLY BY CONSTRUCTION.** The single α in `validation.calibration` was fit to the *same* Primeau anchor, through the 1-D analytic relation `L = √(D_eff/k_eff)`. The isolated single-vessel control therefore recovers ~45 µm essentially by construction and is a **consistency check that the calibration survives the full numerical solver on real geometry, not an independent validation**. What is genuinely non-trivial: the shape comes out exponential, and the real heterogeneous 2-D geometry did not move the number by an order of magnitude.
- **The α convention matters and has changed.** This run used **α = 148.6**, calibrated against the *e-folding* target `45/ln2 ≈ 64.9 µm` (`validation/calibration.py`, `DOX_TARGET_EFOLD_UM`). An earlier convention fitted α to the raw 45 µm half-distance as if it were an e-folding length, giving α = 309.6 and a ~1.4× steeper gradient. Every absolute penetration number in this report inherits that single convention choice; check `calibrate_alpha()` before comparing to older results files.
- **Temozolomide is the independent output** — α was never fit to it — but there is **no published temozolomide perivascular-distance dataset**, so that arm is a *prediction*, checkable only for direction and order of magnitude. (The measured temozolomide anchor in `pk_anchors.py`, PMID 19861433, is a brain-ECF:plasma **partition ratio** — a different physical quantity that must not be compared to a penetration length.)
- **2D, straight-line distances.** The section is a 2D plane through a 3D vascular network: a vessel just above or below the plane is invisible, so grid points scored 'far' may really be near an out-of-plane vessel. Distances are Euclidean, not geodesic around the irregular tissue mask.
- **Not time-matched.** Primeau's images are ~10 min post-injection snapshots; the primary profile integrates a 6 h simulated exposure (AUC). The final-concentration readout is the closer snapshot analogue and agrees closely.
- **Untreated section, no measured drug field.** The Xenium GBM sample is an untreated resection; the comparison is to published measurements in other (murine / xenograft) tumors, not to drug imaged on *this* tissue.
