# Perivascular drug gradient on the REAL tissue — predicted vs measured

> **Regenerated under the discrete vessel source.** `vascular_mode = discrete` (`reconstruct.real_fields._DEFAULT_VASCULAR_MODE`, promoted in [`results/vessel_promotion.md`](vessel_promotion.md)) · **alpha = 148.6**, the e-folding-length convention `L = d_half / ln 2` · per-cell sensitivity mode = **`uniform`** (the reconstruction default) · **dx = 12 µm** on the primary arm, with 12 / 12 / 30 µm across the 3 arms.
>
> **`empirical_half_distance_um` renders `None` honestly, never a fabricated number.** It returns `None` **by design** when a profile never falls to half its perivascular value; every read renders **"n/a (never halves)"** instead of crashing or guessing. **No substitute number is ever fitted or extrapolated** — a profile that does not halve HAS no half-distance, and that is a measured property of the profile, not a gap to be filled. Whether temozolomide's real-tissue profile *does* halve depends on BOTH which distance metric is used AND which delivery law is used — see the STATUS block immediately below (the delivery law, which is undecided) and the distance-metric callout after it. Those two questions are independent and are kept independent throughout.
>
> ### STATUS: this is a TWO-ARM document. Read the arm label on every number.
>
> **The delivery law is an OPEN decision that has not been made**, so this
> report carries both candidate laws side by side and labels every affected
> number with the arm it came from. Neither arm is "history".
>
> | arm | delivery law | status as of this writing |
> |---|---|---|
> | **exchange OFF** | one-way leak, `S = PS·Cp` | **the shipped default** — `tmesim/core/config.py`, `transvascular_exchange: bool = False`; verified at runtime, `PDEConfig().transvascular_exchange` → `False`. Every unlabelled number elsewhere in this repo is this arm. **Known to be physically wrong**: with no back-flux the steady state `PS·Cp/(clearance+binding)` is unbounded above the plasma level, and 5 of 11 GBM drugs sit above it (`results/transvascular_exchange_revalidation.md` §2). |
> | **exchange ON** | gradient-driven, `S = PS·(Cp − C)` | `PDEConfig(transvascular_exchange=True)`. Physically correct *form*, same `PS`, no new parameter; steady state `PS·Cp/(PS+clearance+binding) ≤ Cp` by construction. **Not the default**, and its absolute scale is conditional on an unjustified implicit `k₀ = 1 s⁻¹` (revalidation §6.2). |
>
> ### ⚠ A RESULT IS LOST UNDER THE CORRECTED LAW — stated plainly, not softened
>
> **Under exchange ON there is NO measurable temozolomide penetration depth on
> this section.** The profile never falls to half its perivascular value
> anywhere in the measurable range, so `empirical_half_distance_um` returns
> `None` and no half-distance exists to quote. The exchange-OFF value of
> **153.6 µm** has no exchange-ON counterpart — not a shifted number, an
> absent one. This is a **loss of a result**, and it is the whole point of
> this callout: if the corrected physics is adopted, the 2026-07-29 headline
> ("temozolomide's profile now measurably halves") does not survive, and the
> honest statement becomes "we cannot measure a temozolomide penetration depth
> here." The fitted decay length still exists and still rises,
> **L = 277.4 µm (OFF) → 396.3 µm (ON), +42.9 %**, but a fitted `L` is a model
> output, not the model-free readout the experiment reports.
>
> Mechanism: exchange removes drug fastest exactly where the source is
> strongest, which flattens the perivascular shoulder.
>
> **The loss is on the PRIMARY (direct-cell) distance metric, and it is total
> there.** Under exchange ON, temozolomide has **no** empirical half-distance on
> **any** of the three arms under the primary metric (dx12 repaired, dx12 native,
> dx30 native); under exchange OFF it had one on two of the three. Under the
> LEGACY percentile-mask metric the picture is mixed and must not be reported as
> agreement: on the primary arm both metrics say "never halves" in the ON arm,
> but on the dx12-native arm the legacy metric still halves (141.8 → 153.4 µm)
> and on the dx30 arm it still halves (203.8 → 251.4 µm). So the
> primary-vs-legacy temozolomide contrast that the distance-metric section below
> is built around collapses **on the primary arm** and does not behave uniformly
> across the others — which is itself a reason not to lean on that contrast.
>
> **And the OFF-arm result was marginal by construction.** Exactly **16**
> distance bins clear `MIN_BIN_COUNT = 30` on the primary arm (measured,
> identical in both delivery-law arms), with bin centres at 5, 15, … **155 µm**.
> The 153.6 µm crossing was interpolated **inside that final bin**. A result
> resting on a last-bin crossing is fragile independent of any physics change
> and should not have been quoted as a clean verdict flip.
>
> ### The three negative controls, which are what make the reversal credible
>
> 1. **Doxorubicin barely moves.** Primary arm: d½ 41.5 → 42.2 µm (**+1.76 %**),
>    L 43.5 → 44.5 µm (**+2.22 %**); LEGACY-metric d½ 31.8 → 32.4 µm (+1.91 %),
>    L 34.0 → 34.7 µm (+2.23 %). Across **all twelve** doxorubicin readouts
>    (3 arms × 2 distance metrics × {d½, L}) the OFF→ON shift spans
>    **+0.12 % to +2.64 %**, and **every doxorubicin verdict is unchanged in
>    every arm**. Compare temozolomide's fitted `L` over the matching set
>    (3 arms × 2 metrics): **+33.4 % to +122.7 %**, and on the primary metric at
>    dx = 12 µm, +39 % to +43 %.
> 2. **The isolated single-vessel control is invariant in shape** — and this
>    was tested non-vacuously (see the propagation note below): with a genuinely
>    non-zero `PS` on the vessel line, the perivascular amplitude falls by up to
>    **3.5×** (doxorubicin) and **7.7×** (temozolomide), while `L` holds at
>    **72.52 µm** (doxorubicin) and **293.23 → 293.25 µm**
>    (temozolomide, +0.01 %) with log-space R² = 1.00000 in every arm. Exchange
>    alters the vessel-line amplitude and not the decay away from it. **An
>    isolated-vessel control that stays invariant while the tissue result
>    reverses is direct evidence the reversal is physical, not a bug**: the
>    tissue reversal is driven by the *distributed* source, which is the only
>    thing the isolated geometry removes.
> 3. **Every shape claim survives.** Exponential on the isolated vessel
>    (R²log = 1.00000, both arms); exponential-plus-floor on real tissue in both
>    arms (temozolomide R² 0.9938 OFF / 0.9892 ON; doxorubicin 0.9892 / 0.9895).
>    The source null is **identical by construction** in both arms — it bins the
>    raw `vascular_density` field and solves no PDE — so the source-field
>    confound is unchanged. The distance-metric change is orthogonal to the
>    delivery law: both metrics are computed from the same simulated field in
>    each arm.
>
> Temozolomide's qualitative ordering survives and the gap widens (6.4× → 8.9×
> doxorubicin's fitted L).
>
> ### The silent-no-op trap, and how it was defeated
>
> `results/transvascular_exchange_revalidation.md` §6 lists `vessel_profile` as
> a module where `PDEConfig(transvascular_exchange=True)` is a **silent no-op**.
> That is true **of one function in the file, not the file**, and the difference
> decides which numbers above can be trusted
> (`tmesim/validation/propagation.py::DECLARATIONS`, enforced by
> `tests/test_propagation.py::test_vessel_profile_is_split_not_averaged`):
>
> | function | route to the solver | category |
> |---|---|---|
> | `analyze_drug` (all three tissue arms **and** the KDE-bandwidth sweep) | `sim.engine.simulate_response` | **PROPAGATED** — the flag really reaches the solver |
> | `isolated_vessel_profile` (the single-vessel control) | `simulate_drug_transport` **directly**, no `exchange_k` | **INERT** — the flag is a silent no-op |
>
> Every run behind this document was executed under
> `tmesim.validation.propagation.instrument`, which records whether `exchange_k`
> actually crossed the solver boundary, and every OFF/ON pair was compared with
> `propagation.digest`, a **bit-level** structural hash. Measured:
>
> * **Tissue arms — the comparison is real.** `exchange_k` reached the solver on
>   1/1 call in every ON run and 0/1 in every OFF run, and **the AUC fields are
>   not bit-identical** (primary arm, doxorubicin `ace1396266a6` → `02fb27369d66`;
>   temozolomide `e55ee7a84eb0` → `41301c90c4d0`, digests truncated).
> * **Isolated control — the trap fired, and was caught.** Setting the flag and
>   changing nothing else produced a **bit-identical** AUC field
>   (`dbe618f75210` for doxorubicin, `1d76c8f4d5e1` for temozolomide, both arms),
>   confirming the no-op directly rather than inferring it.
> * **A second trap under the first.** The engine-faithful injection
>   `exchange_k = exchange_coefficient(micro, drug)` is **identically zero** on
>   the isolated geometry, because `_uniform_micro` sets `vascular_density = 0`
>   and hand-paints a unit line source instead — so it too returns a bit-identical
>   field and would have produced a *vacuous* "invariant". The control was
>   therefore re-run with `PS` reconstructed from the geometry's own
>   normalisation (`source = PS·Cp` with peak 1.0) under **both** readings of the
>   undefined `Cp`: `PS = src` (Cp = 1; `PS = 1 s⁻¹` on the vessel line) and
>   `PS = src·bbb_penetrance` (Cp = 1/bbb). Both change the field (digests
>   `eabaafc8e963` / `f5c69791fb40` for doxorubicin) and neither changes `L`.
>   The invariance in control 2 above is that non-vacuous result.
>
> Runs: `experiments/vessel_profile_arms.py`. **The exchange-OFF arm was re-run,
> not assumed**, and it reproduced every published number it covers: all three
> tissue arms × both drugs × both distance metrics (primary arm d½ 41.5, L 43.5,
> R² 0.9892, LEGACY 31.8 / 34.0, corr 0.318, 16 usable bins; native arm 42.0 /
> 45.1 and 144.8 / 280.0; dx = 30 µm arm 13.0 / 7.4 and n/a / 181.0), plus both
> isolated single-vessel fits. So the two arms differ by the delivery law and by
> nothing else. **Scope limit:** the KDE-bandwidth sweep was re-run in the ON arm
> only; its OFF rows are carried from the previous generation, and the one
> overlap that could be checked (bw = 40 µm, repaired mask) reproduced the
> primary arm bit-for-bit in the ON arm.
>
> **One correction to a claim made elsewhere in this repo.** A "pre-existing
> drift" in the vessel-distance map (tissue p90 85.26 µm vs a published 84.85)
> was reported during this re-validation effort. It is **not a drift** — it is
> two different metrics of the same tissue: LEGACY percentile-mask p90 =
> 84.852814 µm (bit-matching `.dose_schedule_cache.json`) and PRIMARY
> direct-cell p90 = 85.264116 µm. Both reproduce exactly. No number in this
> report is affected; it quotes median/p95/max, all of which reproduced.
>
> **Regenerating this file from the module reproduces the exchange-OFF arm
> only.** `tmesim/validation/vessel_profile.py::write_report` has no exchange-ON
> path; the ON numbers and this callout are maintained by hand from
> `experiments/vessel_profile_arms.py` and would be dropped by a bare
> `python -m tmesim.validation.vessel_profile`.

> **THE PRIMARY DISTANCE-TO-VESSEL METRIC CHANGED (2026-07-29) — read this before quoting any number below.** Every fitted number in this report was previously computed from a **percentile threshold on the smoothed/discrete `vascular_density` FIELD** (`vessel_mask_from_density` → `distance_transform_edt`). It is now computed from the **direct distance to the nearest REAL measured endothelial CELL** (`direct_vessel_distance_um`, a single `cKDTree` query, no threshold, no bandwidth, no mask) — a straightforward audit (`results/sweep/vessel_distance_audit_percentile_method_2026-07-29T05-54-00Z.md`) found that a percentile on a density field answers "how far to the densest few percent of the field", not "how far to a vessel", and a prototype (`results/sweep/vessel_distance_direct_replacement_prototype_2026-07-29T05-54-00Z.md`) confirmed the direct-cell distance matches `vascular_validation.py`'s independently-measured tumour→endothelial distance while the percentile mask does not. **The aggregate distance barely moves (median 45.0 µm primary vs 43.3 µm legacy, 0.96×) — this is a TAIL effect, not a mean effect, and it flips one drug's verdict outright:** **temozolomide** goes from *"NO DECAY (no half-distance could be measured)"* (legacy) to *"FAIL (more than 3x off the measured band)"* (primary) — its profile now measurably halves (at a large distance, well outside the doxorubicin band, which is exactly what a small, fast-diffusing, weakly-cleared drug should do) where the old metric could not even measure a half-distance at all. **Doxorubicin** IMPROVES from *"PASS(wide) (outside 40-50 um, inside the 25-75 um replication band)"* to *"PASS (inside the measured 40-50 um band)"*. Neither change is a change in the simulated physics (the AUC/concentration fields are identical; only the distance axis they are binned against changed) — see the dedicated section below for the full head-to-head numbers, and the module HONESTY note for why the LEGACY metric is retained rather than deleted.
>
> **[ARM LABEL — the callout above is an exchange-OFF statement.]** "Its profile now measurably halves" is true under the shipped one-way delivery law and **false under the corrected gradient-driven law**, where temozolomide's profile does not halve on either distance metric on this arm. The doxorubicin half of that callout holds in both arms. See the STATUS block at the top.
>
> Absolute exposure / AUC / coverage figures below are on the **discrete-source unit** and are NOT comparable with the KDE-era numbers this file carried before the discrete-source regeneration (the tissue-integrated source total is 3.36x smaller; that is a change of *unit*, not of physics).

**Module:** `tmesim/validation/vessel_profile.py` (run: `./.venv/bin/python -m tmesim.validation.vessel_profile`)

**Figure:** `results/vessel_profile.png`

The canonical experimental readout of tumor drug penetration is a **profile**, not a scalar: quantitative immunofluorescence of drug intensity against **distance from the nearest CD31⁺ blood vessel**. tmesim predicts a full spatial drug field on **real** tissue, so we can compute *the same quantity* and compare the **functional form**.

## Headline

Both delivery-law arms, side by side. **Exchange OFF is the shipped default; exchange ON is the physically correct form.** Neither column is the "old" one.

| | doxorubicin, **exch OFF** | doxorubicin, **exch ON** | temozolomide, **exch OFF** | temozolomide, **exch ON** |
|---|---:|---:|---:|---:|
| **fitted L** (µm), real tissue, primary arm | **43.5** | **44.5** | **277.4** | **396.3** |
| **half-distance d½** (µm), real tissue, primary arm | **41.5** | **42.2** | **153.6** | **n/a (never halves)** |
| fitted L (µm), isolated single vessel | 72.52 | 72.52 | 293.23 | 293.24–293.25 |
| d½ = L·ln2 (µm), isolated single vessel | 50.27 | 50.27 | 203.25 | 203.26 |
| measured | **40–50 µm** (PMID 16361566), 25–75 µm (PMID 23680920) | ← same | *no published perivascular-distance dataset* | ← same |

The isolated single-vessel ON column spans the two readings of the control's undefined plasma normalisation (`PS = src` and `PS = src·bbb`); both are reported in full in the two-arms section below. **The single number that disappears rather than moves is temozolomide's empirical half-distance.**

### Read this before quoting any number

Two different quantities are easy to conflate, and only one of them is clean:

1. **The analytic / single-vessel penetration length** `L = √(D_eff / (clearance + binding))`, which is genuinely drug-specific (doxorubicin 72.5 µm vs temozolomide 295.4 µm) and is what the isolated single-vessel control reproduces with the real solver. **This is a legitimate model output and is the thing to compare to the literature.**
2. **The spatial fit to drug-vs-distance on our reconstructed tissue.** This one is **still confounded, but by one thing instead of two.** The drug's perivascular-sink SOURCE term is proportional to the (smoothed or discrete) `vascular_density` field, so part of the fitted decay is that field's own shape — the source null (d½ = 10.0 µm with **no drug physics at all**) and the bandwidth sweep below quantify how much. This regeneration REMOVES a second, previously-conflated confound: the vessel-DISTANCE axis used to be a threshold on that *same* KDE field, so part of the fitted decay used to also be our own 40 µm mask-smoothing bandwidth. It is now the direct nearest-real-endothelial-cell distance (see the dedicated section below), decoupled from whatever field the drug is delivered through — a legitimate, separately-debated modelling choice — from what counts as "a vessel" for the distance axis, which should always be the actual measured geometry. On the fragmented native dx = 12 µm mask the remaining (source-field) confound is *total* — every drug returns close to the source-null decay, reproducing the independent finding in `results/delivery_metrics.md`. **A match to the measured 40–50 µm band read off the spatial fit is not a validation.** Quote the drug/source-null ratio instead.

**Shape:** exponential, **in both delivery-law arms**. On the isolated single vessel the predicted profile is a pure exponential to five decimals (log-space R² = 1.00000, OFF and ON). On real tissue it is **exponential decaying onto a floor**, `C(d) = C₀·exp(−d/L) + floor` (doxorubicin R² = 0.9892 OFF / 0.9895 ON; temozolomide 0.9938 OFF / 0.9892 ON), the floor being the drug delivered by the *distributed* vascular source. Published profiles also flatten into a background, so the qualitative form matches. **No shape claim in this report depends on which delivery law is used.**

**Verdict for doxorubicin — identical in both arms.** Real tissue d½ = 41.5 µm (OFF) / 42.2 µm (ON) → *PASS (inside the measured 40-50 um band)* in both; the fitted `L = 43.5 µm` (OFF) / 44.5 µm (ON) against the L-equivalent band `40–50/ln2 ≈ 58–72 µm` → *PASS(3x) (steeper than measured by 33%, same order of magnitude)* OFF and *PASS(3x) (steeper than measured by 31%, same order of magnitude)* ON — same verdict class, both arms. **This is not an independent validation** — α was fit to this same anchor.

**Temozolomide is the free prediction, and this is where the two arms diverge.** There is no *published, temozolomide-specific* perivascular-distance dataset to score it against either way. Under **exchange OFF** (L = 277.4 µm, 6.4× doxorubicin's) the PRIMARY distance metric gives a measurable empirical d½ of 153.6 µm, which CAN be sanity-checked against the doxorubicin band as an order-of-magnitude check — and it *FAIL (more than 3x off the measured band)*: an honest failure of a generic check, not a validated temozolomide measurement. Under **exchange ON** (L = 396.3 µm, 8.9× doxorubicin's) **there is no empirical d½ at all** — the profile never falls to half its perivascular value, so the verdict is *NO DECAY (no half-distance could be measured)* and the generic check cannot be run. **Under the corrected delivery law, no temozolomide penetration depth is measurable on this section.** See the two-arms section for the full numbers and the distance-metric section for the separate, orthogonal question of which distance axis is used.

---

## THE DISTANCE-TO-VESSEL METRIC CHANGE — PRIMARY (direct-cell) vs LEGACY (percentile-mask)

This regeneration replaces the vessel-distance axis that drives every fit above. **What changed:** previously, "distance from nearest vessel" was computed by thresholding the (smoothed or discrete) `vascular_density` FIELD at its 95th tissue percentile and running a Euclidean distance transform on the resulting mask (`vessel_mask_from_density` → `distance_from_vessel_um`). It is now computed directly from the REAL measured endothelial-cell positions — a single `cKDTree` nearest-neighbour query, no threshold, no bandwidth, no mask at all (`direct_vessel_distance_um`, reusing `vascular_validation.distance_to_nearest`, the exact function `vascular_validation.py` uses for its own independently-published tumour→endothelial distance, so the two numbers cannot drift apart). **The percentile-mask pipeline is NOT deleted** — it is kept as a documented LEGACY metric (`vessel_mask_from_density` / `distance_from_vessel_um` are unchanged functions) for backward comparison and for datasets where only a density field, not individual cell positions, is available.

**Why this matters: it is a tail effect, not a mean effect.** Grid-wide, the two distance fields barely differ (median 45.0 µm primary vs 43.3 µm legacy, ratio 0.96×; max 789 vs 786 µm) — the percentile mask happens to be geometrically close to the real vessel footprint almost everywhere (see the audit's §4 for why: under the discrete source field the percentile threshold degenerates to the same `vasc > 0` set the direct metric is measuring against, by coincidence of the fallback branch, not by design). But the two metrics diverge in the FAR tail — the tumour cells genuinely furthest from any vessel — and that tail is exactly where a slow or absent perivascular sink shows up, which is precisely temozolomide's regime.

**The distance-metric question and the delivery-law question are orthogonal and are kept separate here.** Both metrics are post-hoc reductions of the *same* simulated field, so the metric comparison can be — and is — run inside each delivery-law arm. Both are shown, because the temozolomide contrast this section is built around exists in one arm and not the other.

| | doxorubicin, **exch OFF** | doxorubicin, **exch ON** | temozolomide, **exch OFF** | temozolomide, **exch ON** |
|---|---:|---:|---:|---:|
| empirical d½ (µm) — PRIMARY (direct-cell) | **41.5** | **42.2** | **153.6** | **n/a (never halves)** |
| empirical d½ (µm) — LEGACY (percentile-mask) | 31.8 | 32.4 | n/a (never halves) | n/a (never halves) |
| fitted L, exp+floor (µm) — PRIMARY | **43.5** | **44.5** | **277.4** | **396.3** |
| fitted L, exp+floor (µm) — LEGACY | 34.0 | 34.7 | 276.8 | 393.5 |
| verdict vs measured band — PRIMARY | *PASS (inside the measured 40-50 um band)* | *PASS (inside the measured 40-50 um band)* | *FAIL (more than 3x off the measured band)* | *NO DECAY (no half-distance could be measured)* |
| verdict vs measured band — LEGACY | *PASS(wide) (outside 40-50 um, inside the 25-75 um replication band)* | *PASS(wide) (outside 40-50 um, inside the 25-75 um replication band)* | *NO DECAY (no half-distance could be measured)* | *NO DECAY (no half-distance could be measured)* |

**On this (primary) arm the temozolomide PRIMARY-vs-LEGACY contrast exists only in the exchange-OFF column.** Under exchange ON both distance metrics return "never halves" here, so the metric change makes no difference to temozolomide on this arm and the paragraphs below describe an OFF-arm finding. (On the other two arms the ON-arm legacy metric *does* still halve — 153.4 µm at dx = 12 µm native, 251.4 µm at dx = 30 µm — while the primary metric does not; see the register near the end.) Doxorubicin's metric contrast is present and essentially identical in both arms.

**Read plainly, not smoothed over (exchange-OFF arm unless stated):**

- **Doxorubicin's verdict IMPROVES**: *"PASS(wide) (outside 40-50 um, inside the 25-75 um replication band)"* (legacy) → *"PASS (inside the measured 40-50 um band)"* (primary). The legacy percentile-mask distance was systematically slightly SHORTER than the true nearest-vessel distance near the perivascular zone the doxorubicin fit is most sensitive to, which compressed its fitted decay length below the measured band; the direct-cell distance corrects that and lands the fit inside the measured 40–50 µm band.
- **Temozolomide's verdict CHANGES**: *"NO DECAY (no half-distance could be measured)"* (legacy) → *"FAIL (more than 3x off the measured band)"* (primary). This is the more consequential change and it is stated here explicitly, not buried: under the legacy metric, temozolomide's profile never fell to half its perivascular value anywhere on the section (`empirical_half_distance_um` correctly returned `None` — it is not a bug, the profile genuinely never crossed 50% within the legacy distance range), so no verdict could be computed at all. Under the PRIMARY metric the profile DOES cross 50%, far out in the tail where the legacy metric under-measured distance, and the resulting half-distance is more than 3× the measured doxorubicin band — an honest, scoreable **FAIL** on the generic order-of-magnitude check. This is the physically expected direction (temozolomide is small, fast-diffusing and has no fast perivascular sink, so its profile SHOULD extend far beyond a doxorubicin-scaled band) but it is also a genuine change in what this report says about temozolomide: it went from "unscoreable" to "scoreable and, on this generic check, failing". Nothing about the simulated transport physics changed *within this arm* — the exchange-OFF AUC field for temozolomide is bit-for-bit the same as before; only the distance axis it is binned against did.
  **Under exchange ON this entire paragraph has no counterpart**: the PRIMARY metric also returns "never halves", so temozolomide is back to *unscoreable* on both metrics and the verdict flip described here does not occur. The corrected delivery law therefore *undoes* the one result this distance-metric regeneration was written to announce. Which of the two states is the right one to publish depends on the delivery-law decision, which has not been made.

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

**At dx = 12 µm — exactly the resolution `numerical_rigor.md` recommends — the domain is 26,058 disconnected fragments whose largest component holds 0.7% of the tissue.** Lateral diffusion is effectively switched off. That is why, on the native dx = 12 µm mask, doxorubicin and temozolomide produce *the same* profile (42.0 vs 144.8 µm under exchange OFF; 42.8 µm vs no half-crossing under exchange ON) despite ~4× different transport lengths, and both sit on top of the source null (10.0 µm) — correlation between the AUC field and the source field 0.319 and 0.114 (OFF), 0.324 and 0.135 (ON). This is a grid artefact, not physics, **and it is present in both delivery-law arms**.

We therefore report **three arms** and take as primary the one that resolves the boundary layer *and* has a connected domain:

**Exchange OFF (the shipped default):**

| arm | dx (µm) | largest connected component | doxorubicin d½ (µm) | doxorubicin L (µm) | temozolomide d½ (µm) | temozolomide L (µm) | source-null d½ (µm) | corr(AUC, source): dox / TMZ |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| dx = 12 µm, morphologically REPAIRED (connected) mask **(PRIMARY)** | 12 | 98.3% | **41.5** | 43.5 | **153.6** | 277.4 | *10.0* | 0.318 / 0.111 |
| dx = 12 µm, NATIVE cell-derived mask (fragmented) | 12 | 98.0% | **42.0** | 45.1 | **144.8** | 280.0 | *10.0* | 0.319 / 0.114 |
| dx = 30 µm, NATIVE mask (connected, but boundary layer coarse) | 30 | 98.1% | **13.0** | 7.4 | **n/a (never halves)** | 181.0 | *10.0* | 0.548 / 0.192 |

**Exchange ON (`PDEConfig(transvascular_exchange=True)`), same geometry, same α, same everything else:**

| arm | dx (µm) | largest connected component | doxorubicin d½ (µm) | doxorubicin L (µm) | temozolomide d½ (µm) | temozolomide L (µm) | source-null d½ (µm) | corr(AUC, source): dox / TMZ |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| dx = 12 µm, morphologically REPAIRED (connected) mask **(PRIMARY)** | 12 | 98.3% | **42.2** | 44.5 | **n/a (never halves)** | 396.3 | *10.0* | 0.323 / 0.135 |
| dx = 12 µm, NATIVE cell-derived mask (fragmented) | 12 | 98.0% | **42.8** | 46.3 | **n/a (never halves)** | 389.8 | *10.0* | 0.324 / 0.135 |
| dx = 30 µm, NATIVE mask (connected, but boundary layer coarse) | 30 | 98.1% | **13.0** | 7.4 | **n/a (never halves)** | 403.1 | *10.0* | 0.553 / 0.240 |

The source-null column is **identical by construction** in the two arms: it bins the raw `vascular_density` field and solves no PDE, so no delivery law can touch it. **Under exchange ON, temozolomide has no measurable half-distance on any of the three arms** — the readout that survives is the fitted `L`, which rises on every arm.

The **repair** is pure image morphology applied to a *copy* of the mask (`scipy.ndimage.binary_closing`, 2 iterations, then `binary_fill_holes`), which restores the largest component to 98.3% of tissue; the ECM and vascular fields are then **re-derived on the repaired mask** with the public `reconstruct.real_fields.build_real_fields`, because the KDE is multiplied by the tissue mask at build time and reusing the old fields would leave the added pixels with a spurious zero source. **Nothing in `reconstruct` is modified.** The added pixels are the interstitial gaps *between* measured cells — real tissue the drug diffuses through, excluded by the native mask only because the rasteriser found no cell centroid there.

The dx = 30 µm arm (native mask, 98% connected) is an independent cross-check with no morphology applied at all, and it agrees with the repaired arm qualitatively in **both** delivery-law arms: doxorubicin 13.0 µm vs temozolomide n/a (never halves) µm. **Honest range across the two physically-sane arms: doxorubicin d½ ≈ 13–42 µm (OFF) / 13–42 µm (ON), temozolomide d½ ≈ n/a (never halves) in both.**

One caveat on the dx = 30 µm arm's temozolomide `L`, which moves the most of any number in this report (181.0 → 403.1 µm, +123 %): that arm's exponential-plus-floor fit is the weakest in the report (R² 0.856 OFF / 0.827 ON, against ≥ 0.976 on the dx = 12 µm arms), so the large percentage is a fit on a poorly resolved boundary layer, not a sharper measurement. It is reported because it is what the run produced, not because it is the strongest evidence.

---

## Result 1 — the gradient on real tissue (primary arm)

**Model-free readout** (what the experiment reports — the distance at which the binned mean falls to half its perivascular value, no functional form assumed):

| signal | **empirical d½** (µm), exch OFF | **empirical d½** (µm), exch ON | corr(AUC, source), OFF | corr(AUC, source), ON |
|---|---:|---:|---:|---:|
| **doxorubicin** AUC | **41.5** | **42.2** | 0.318 | 0.323 |
| **temozolomide** AUC | **153.6** | **n/a (never halves)** | 0.111 | 0.135 |
| *SOURCE NULL: `vascular_density` itself* | *10.0* | *10.0 (identical by construction)* | 1.000 | 1.000 |

**This is the model-free row, and it is the one that breaks.** The empirical half-distance assumes no functional form — it is simply the distance at which the binned mean first falls to half its perivascular value — which is exactly why its disappearance under exchange ON cannot be argued away as a fitting artefact. There is no half-crossing to interpolate. The **final-concentration** field behaves identically (doxorubicin 41.5 → 42.3 µm; temozolomide 153.9 µm OFF → **no crossing** ON), so the loss is not an artefact of time-integration either.

Both arms use **16** usable distance bins reaching **155 µm**, so the measurable range is the same in each and the difference is not a change in what could be seen.

**Fits.** A *pure* exponential is the wrong model on real tissue: the distributed source puts a plateau under the perivascular decay and drags a single-exponential fit towards it. The three-parameter `C(d) = C₀·exp(−d/L) + floor` fits far better and its `L` is the real perivascular decay length. All R² are reported in **linear space** so the models are directly comparable:

**Exchange OFF (the shipped default):**

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

**Exchange ON (`PDEConfig(transvascular_exchange=True)`):**

| signal | model | **L** (µm) | **d½ = L·ln2** (µm) | floor (% of perivascular) | **R² (linear space)** | R² (log space) |
|---|---|---:|---:|---:|---:|---:|
| **doxorubicin** | `C₀·exp(−d/L) + floor` (full range) | **44.5** | **30.8** | 14.2% | **0.9895** | — |
| doxorubicin | pure `C₀·exp(−d/L)`, near field d ≤ 150 µm | 69.5 | 48.2 | — | 0.9757 | 0.9919 |
| doxorubicin | straight-line null `a + b·d` (near field) | — | — | — | 0.8883 | — |
| **temozolomide** | `C₀·exp(−d/L) + floor` (full range) | **396.3** | **274.7** | 0.0% | **0.9892** | — |
| temozolomide | pure `C₀·exp(−d/L)`, near field d ≤ 150 µm | 394.2 | 273.2 | — | 0.9923 | 0.9914 |
| temozolomide | straight-line null `a + b·d` (near field) | — | — | — | 0.9962 | — |
| **SOURCE NULL** | `C₀·exp(−d/L) + floor` (full range) | **n/a** | **n/a** | n/a% | **n/a** | — |
| SOURCE NULL | pure `C₀·exp(−d/L)`, near field d ≤ 150 µm | n/a | n/a | — | n/a | n/a |
| SOURCE NULL | straight-line null `a + b·d` (near field) | — | — | — | n/a | — |

**The `d½ = L·ln2` column is a MODEL-DERIVED number and must not be substituted for the missing empirical one.** Temozolomide's ON row reads 274.7 µm, but that is `L·ln2` for a fit whose profile never actually halves inside the measured range — it is an extrapolation of the fitted form far past the last usable bin (155 µm), which is why the model-free table above says "n/a" and this one does not. Do not quote 274.7 µm as a temozolomide penetration depth.

**Temozolomide's exchange-ON fit is also the one place where the straight-line null wins.** Its near-field R² is 0.9962 for a straight line against 0.9923 for a pure exponential — a profile that flat over 0–150 µm is barely distinguishable from linear, which is the same fact as "it never halves" seen from the fitting side. The doxorubicin near field beats its line null decisively in both arms (0.9757 vs 0.8883 ON; 0.9740 vs 0.8848 OFF).

The **final-concentration** field gives essentially the same profile in both arms (doxorubicin empirical d½ 41.5 µm OFF / 42.3 µm ON; temozolomide 153.9 µm OFF and **no crossing** ON) — expected, because with a constant-source window the field is quasi-steady, so `AUC ≈ C_steady · Δt` and the two share a spatial shape. It is *not* an independent check, but it does rule out the profile being an artefact of time-integration, in either arm.

### How much of the gradient is the drug, and how much is the vasculature?

The engine's delivery term is `vascular_density · bbb_penetrance · plasma`: the source is spread over the **whole continuous density field**, not confined to the thresholded mask (mean vascular density over tissue = 0.026; tissue 200 µm from the mask still carries nan% of the peak density). The SOURCE NULL row measures what that alone buys: an apparent half-distance of **10.0 µm** with no drug physics at all — **the same 10.0 µm in both delivery-law arms**, because the null solves no PDE. On the primary arm the drugs sit clearly above it in both arms (OFF: 41.5 and 153.6 µm, correlations 0.32 and 0.11; ON: 42.2 µm and *no half-crossing at all*, correlations 0.32 and 0.14), so transport is doing real work — but a meaningful part of the profile is still vessel/KDE geometry, and that share must be subtracted mentally before reading the agreement with the literature as evidence. Note that under exchange ON temozolomide cannot even be *compared* to the null on this readout, because it has no half-distance.

---

## Result 2 — CONTROL: the isolated single-vessel gradient

To see the gradient the *physics* predicts, free of the distributed-source confound, we run the **same frozen transport solver** and the **same compiled drug** on one planar (line) vessel in a uniform medium at the real tissue's mean ECM (0.087). For a planar source the steady solution is exactly `C(d) = C₀·exp(−d/L)` with `L = √(D_eff/k)` — the geometry the calibration constant is defined in.

**Exchange OFF (the shipped default):**

| drug | MW (Da) | D (µm²/s) | k = clearance+binding (1/s) | **L** (µm) | **d½ = L·ln2** (µm) | R² (log) |
|---|---:|---:|---:|---:|---:|---:|
| **doxorubicin** | 544 | 185 | 0.0333 | **72.5** | **50.3** | 1.00000 |
| **temozolomide** | 194 | 278 | 0.00302 | **293.2** | **203.3** | 1.00000 |

### This control under both delivery laws — and why the naive test of it is worthless

This is the one part of the module where `PDEConfig(transvascular_exchange=True)` is a **silent no-op** (`isolated_vessel_profile` calls `simulate_drug_transport` directly with no `exchange_k`; declared INERT in `tmesim/validation/propagation.py`). Worse, the obvious remedy is also empty here: `exchange_coefficient(micro, drug)` is **identically zero** on this geometry, because `_uniform_micro` sets `vascular_density = 0` and hand-paints a unit line source instead. Both of those produce a bit-identical AUC field and a vacuous "nothing changed".

So `PS` was rebuilt from the control's own normalisation — `source = PS·Cp` with a peak of 1.0 — under both readings of its undefined `Cp`, and the flag was injected as `exchange_k` explicitly:

| drug | arm | `PS` on the vessel line (1/s) | perivascular AUC | **L** (µm) | **d½ = L·ln2** (µm) | R² (log) | AUC field |
|---|---|---:|---:|---:|---:|---:|---|
| doxorubicin | exchange OFF (shipped) | — | 53 434.7 | **72.52** | **50.27** | 1.00000 | `dbe618f75210` |
| doxorubicin | flag set, `exchange_k` NOT passed | — | 53 434.7 | 72.52 | 50.27 | 1.00000 | `dbe618f75210` — **bit-identical, the no-op** |
| doxorubicin | `exchange_k = exchange_coefficient(micro, dp)` | **0.0000** | 53 434.7 | 72.52 | 50.27 | 1.00000 | `dbe618f75210` — **bit-identical, PS is zero here** |
| doxorubicin | ON, `PS = src` (Cp = 1) | 1.0000 | 15 379.5 | **72.52** | **50.27** | 1.00000 | `eabaafc8e963` — differs |
| doxorubicin | ON, `PS = src·bbb` (Cp = 1/bbb) | 0.0206 | 50 847.0 | **72.52** | **50.27** | 1.00000 | `f5c69791fb40` — differs |
| temozolomide | exchange OFF (shipped) | — | 144 156 | **293.23** | **203.25** | 1.00000 | `1d76c8f4d5e1` |
| temozolomide | flag set, `exchange_k` NOT passed | — | 144 156 | 293.23 | 203.25 | 1.00000 | `1d76c8f4d5e1` — **bit-identical, the no-op** |
| temozolomide | `exchange_k = exchange_coefficient(micro, dp)` | **0.0000** | 144 156 | 293.23 | 203.25 | 1.00000 | `1d76c8f4d5e1` — **bit-identical, PS is zero here** |
| temozolomide | ON, `PS = src` (Cp = 1) | 1.0000 | 18 781.6 | **293.25** | **203.26** | 1.00000 | `790be04daa49` — differs |
| temozolomide | ON, `PS = src·bbb` (Cp = 1/bbb) | 0.1254 | 78 478.6 | **293.24** | **203.26** | 0.99999 | `6f70d79e57ec` — differs |

**What this control establishes.** Exchange cuts the perivascular amplitude by up to **3.5×** (doxorubicin) and **7.7×** (temozolomide) and leaves the **decay length untouched** — doxorubicin `L` identical to five significant figures, temozolomide `L` moving by **+0.01 %**. That is exactly what the algebra predicts: `PS` is non-zero only on the vessel line, so away from it the reaction rate is still `clearance + binding` and `L = √(D_eff/k)` is unchanged. It is also why the α calibration chain survives the delivery-law change.

**This is the load-bearing negative control for the temozolomide reversal.** The isolated geometry differs from the real tissue in exactly one respect that matters here: it has no *distributed* source. The tissue result reverses and this one does not, which is direct evidence that the reversal is the distributed source interacting with the back-flux term — a physical mechanism — rather than a solver or harness artefact. Had this control also moved, a bug would have been the leading explanation.

**Temozolomide's gradient is 4.0× longer than doxorubicin's** in identical geometry — the physically expected ordering: temozolomide is small (MW 194 vs 544 Da), diffuses faster (278 vs 185 µm²/s) and is not consumed by a high-capacity intracellular/DNA sink (k = 0.003 vs 0.033 /s), so it escapes the perivascular cuff and reaches tissue the anthracycline never sees.

This control is also where the exponential-shape claim is unambiguous: log-space R² ≈ 1.0 for both drugs, i.e. the predicted profile is a **true exponential**, the functional form the measurements report.

---

## Comparison to the real measurements

**The anchors (both PMID-verified, carried in `tmesim/validation/pk_anchors.py`):**

> Primeau AJ, Rendon A, Hedley D, Lilge L, Tannock IF. The distribution of the anticancer drug doxorubicin in relation to blood vessels in solid tumors. Clin Cancer Res 2005;11(24 Pt 1):8782-8. **PMID 16361566** — doxorubicin intensity falls to half its perivascular value within **40–50 µm** of CD31⁺ vessels (quantitative immunofluorescence of doxorubicin autofluorescence vs distance from CD31+ blood vessels; murine 16C and EMT6 tumors + human PC-3 prostate xenografts (mice)).

> Patel KJ, Tredan O, Tannock IF. Distribution of the anticancer drugs doxorubicin, mitoxantrone and topotecan in tumors and normal tissues. Cancer Chemother Pharmacol 2013;72(1):127-38. **PMID 23680920** — independent replication; half-distance **25–75 µm** at 10 min, reported jointly for doxorubicin / mitoxantrone / topotecan.

Every verdict string below was produced by the repo's own `verdict_vs_measurement` / `verdict_L_vs_efold_band`, not written by hand.

| doxorubicin readout | predicted, **exch OFF** | predicted, **exch ON** | measured | verdict, **exch OFF** | verdict, **exch ON** |
|---|---:|---:|---:|:--|:--|
| real tissue (primary arm), **empirical d½** | **41.5 µm** | **42.2 µm** | 40–50 µm (25–75 µm) | PASS (inside the measured 40-50 um band) | PASS (inside the measured 40-50 um band) |
| real tissue (primary arm), fitted **L** vs L-equivalent band | **43.5 µm** | **44.5 µm** | 58–72 µm | PASS(3x) (steeper than measured by 33%, same order of magnitude) | PASS(3x) (steeper than measured by 31%, same order of magnitude) |
| real tissue (dx = 30 µm cross-check), empirical d½ | 13.0 µm | 13.0 µm | 40–50 µm (25–75 µm) | FAIL (more than 3x off the measured band) | FAIL (more than 3x off the measured band) |
| **isolated single vessel**, d½ | **50.3 µm** | **50.3 µm** | 40–50 µm (25–75 µm) | PASS(wide) (outside 40-50 um, inside the 25-75 um replication band) | PASS(wide) (outside 40-50 um, inside the 25-75 um replication band) |
| isolated single vessel, L vs L-equivalent band | 72.5 µm | 72.5 µm | 58–72 µm | PASS(3x) (shallower than measured by 12%, same order of magnitude) | PASS(3x) (shallower than measured by 12%, same order of magnitude) |

**Not one doxorubicin verdict changes.** This is the negative control that licenses reading the temozolomide change as real: the drug the delivery law was calibrated against is unmoved by it, so the temozolomide result is not a global rescaling of everything.

| temozolomide readout | predicted, **exch OFF** | predicted, **exch ON** | verdict, **exch OFF** | verdict, **exch ON** |
|---|---:|---:|:--|:--|
| real tissue (primary arm), **empirical d½**, PRIMARY metric | **153.6 µm** | **none — profile never halves** | FAIL (more than 3x off the measured band) | **NO DECAY (no half-distance could be measured)** |
| real tissue (dx = 12 µm native), empirical d½, PRIMARY metric | 144.8 µm | **none — profile never halves** | FAIL (more than 3x off the measured band) | **NO DECAY (no half-distance could be measured)** |
| real tissue (dx = 30 µm cross-check), empirical d½, PRIMARY metric | none — never halves | **none — profile never halves** | NO DECAY (no half-distance could be measured) | NO DECAY (no half-distance could be measured) |

There is no measured temozolomide perivascular-distance dataset under either delivery law, so these are scored only against the generic doxorubicin band as an order-of-magnitude check. **Under exchange ON that check cannot be run at all on the primary distance metric — on any of the three grid arms.**

**Why two bands.** `pk_anchors.py` is explicit that the published 40–50 µm is a *half*-distance, so the quantity comparable to an e-folding length is `40–50 / ln2 ≈ 58–72 µm`. Comparing an `L` to the raw 40–50 µm (as `pk_reality.py` does) is deliberately conservative. Both comparisons are given above so no reader has to guess which convention a number is in.

### Is the shape exponential, as measured?

**Answered in both delivery-law arms; the answer is the same in both.**

- **doxorubicin** — *isolated vessel*: log-space R² = **1.00000** in every arm (OFF, and ON under both `PS` readings), a clean exponential to five decimals. *Real tissue*: **yes, exponential-plus-floor** — R²(linear) = **0.9892 OFF / 0.9895 ON** for `C₀·exp(−d/L)+floor`, versus 0.9740 / 0.9757 for a pure exponential and 0.8848 / 0.8883 for a straight line over the same near field. The decaying component has L = 43.5 µm (OFF) / 44.5 µm (ON); beyond ~150 µm the profile plateaus at 14.4% (OFF) / 14.2% (ON) of the perivascular value instead of continuing to decay.
- **temozolomide** — *isolated vessel*: log-space R² = **1.00000** in every arm. *Real tissue*: **yes, exponential-plus-floor** — R²(linear) = **0.9938 OFF / 0.9892 ON**, versus 0.9951 / 0.9923 for a pure exponential and 0.9951 / **0.9962** for a straight line over the same near field. The decaying component has L = 277.4 µm (OFF) / 396.3 µm (ON); the fitted floor is ~0% of the perivascular value in both arms. **Note the ON straight-line number**: with the corrected law temozolomide's near field is fit *better* by a line than by an exponential, which is the fitting-side statement of "this profile no longer has a measurable decay".

So the predicted functional form **is** the measured one under either delivery law: exponential where the geometry is clean, exponential-onto-a-background on real tissue. No shape conclusion in this report depends on the delivery-law decision. The published profiles also flatten into a background at large distance — but ours is a modelling consequence of the diffuse source, not a measured autofluorescence floor, and the two should not be equated.

---

## THE DECISIVE CONTROL — is the fitted decay our own KDE bandwidth?

`vascular_density` is a Gaussian KDE of endothelial-cell positions with a **40 µm bandwidth**, and it is the field the drug's perivascular-sink source term is proportional to. If the fitted decay length tracks that bandwidth, the spatial fit is measuring our own smoothing kernel and **cannot** be quoted as agreement with the literature. We therefore rebuild the real fields at 25 / 40 / 60 µm bandwidth and re-measure everything.

**What this sweep tests changed along with the distance metric.** The vessel-DISTANCE axis below is the PRIMARY direct-cell metric, computed ONCE from the real endothelial-cell positions **before** the bandwidth loop — real cell positions do not move when the KDE bandwidth changes, so this axis is by construction identical at bw = 25, 40 and 60 µm. That is deliberate, not a gap: under the OLD percentile-mask distance, this same sweep silently rebuilt an essentially identical mask at every bandwidth anyway (the production `vascular_mode = discrete` field does not even use `bandwidth_um` — see the audit, §6), so it was already testing nothing on the distance axis while *appearing* to test bandwidth sensitivity end-to-end. Fixing the distance axis explicitly makes that scope honest: this control now isolates exactly one thing — how much the fitted `L` and the drug/source-null ratio move because the bandwidth-dependent `vascular_density` FIELD (the delivery source term, and via the companion ECM-KDE, tissue tortuosity) changed, with the "where is a vessel" question held fixed at its one geometrically correct answer.

This sweep runs through `analyze_drug`, so it is on the **PROPAGATED** side of the module and the delivery-law flag genuinely reaches it: `exchange_k` was confirmed to reach the solver on 1/1 call in every ON run. **Scope, stated so it is not over-read:** the exchange-ON columns below were re-run for this document; the exchange-OFF columns are carried unchanged from the previous generation and were *not* re-executed here. The one overlap that could be checked — bw = 40 µm on the repaired mask — reproduced the primary-arm ON run **bit-for-bit** (identical AUC-field digests, `02fb27369d66` and `41301c90c4d0`), which is an internal consistency check on the whole harness.

**Primary arm (dx = 12 µm, repaired connected mask):**

| KDE bandwidth (µm) | source-null L (µm) | dox L, **exch OFF** | dox L, **exch ON** | TMZ L, **exch OFF** | TMZ L, **exch ON** | TMZ empirical d½, ON | corr(AUC, source) ON: dox / TMZ |
|---:|---:|---:|---:|---:|---:|:--|---:|
| 25 | n/a | **44.4** | **45.4** | **280.7** | **400.4** | n/a (never halves) | 0.320 / 0.134 |
| 40 | n/a | **43.5** | **44.5** | **277.4** | **396.3** | n/a (never halves) | 0.323 / 0.135 |
| 60 | n/a | **42.4** | **43.3** | **273.2** | **391.4** | n/a (never halves) | 0.326 / 0.137 |

**Fragmented native dx = 12 µm mask (the arm with no lateral transport), for contrast:**

| KDE bandwidth (µm) | source-null L (µm) | dox L, **exch OFF** | dox L, **exch ON** | TMZ L, **exch OFF** | TMZ L, **exch ON** | TMZ empirical d½, ON |
|---:|---:|---:|---:|---:|---:|:--|
| 25 | n/a | 47.0 | 48.2 | 282.6 | 393.8 | n/a (never halves) |
| 40 | n/a | 46.0 | 47.2 | 279.1 | 389.7 | n/a (never halves) |
| 60 | n/a | 44.8 | 45.9 | 274.8 | 384.8 | n/a (never halves) |

**The bandwidth conclusion is delivery-law-independent, and the temozolomide loss is bandwidth-independent.** In both arms the fitted `L` moves only a few percent across a 25 → 60 µm bandwidth change (doxorubicin −4.5% OFF / −4.6% ON on the primary arm; temozolomide −2.7% OFF / −2.3% ON), so the "is the fitted decay our own smoothing kernel?" answer is unchanged by the delivery law. And temozolomide has **no** empirical half-distance at **any** of the six exchange-ON bandwidth/mask combinations — the loss is not an artefact of one KDE bandwidth.

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

## REGISTER — every temozolomide fitted length in this report, in both arms

The annotation this document replaces asked for exactly this: *"every temozolomide fitted length needs regenerating."* Every **exchange-ON** value below was re-run for this document (`experiments/vessel_profile_arms.py`), as was every **exchange-OFF** value except the six KDE-bandwidth rows, which are carried unchanged from the previous generation and are marked *(carried)*. Where an OFF value was re-run it reproduced the published number exactly.

| where | metric | **L, exch OFF** (µm) | **L, exch ON** (µm) | Δ | **empirical d½, OFF** (µm) | **empirical d½, ON** (µm) |
|---|---|---:|---:|---:|:--|:--|
| dx = 12 µm repaired **(PRIMARY ARM)** | PRIMARY (direct-cell) | **277.4** | **396.3** | **+42.9 %** | **153.6** | **n/a (never halves)** |
| dx = 12 µm repaired | LEGACY (percentile-mask) | 276.8 | 393.5 | +42.2 % | n/a (never halves) | n/a (never halves) |
| dx = 12 µm native | PRIMARY | 280.0 | 389.8 | +39.2 % | 144.8 | **n/a (never halves)** |
| dx = 12 µm native | LEGACY | 273.7 | 378.5 | +38.3 % | 141.8 | 153.4 |
| dx = 30 µm native | PRIMARY | 181.0 | 403.1 | +122.7 % | n/a (never halves) | n/a (never halves) |
| dx = 30 µm native | LEGACY | 235.9 | 314.6 | +33.4 % | 203.8 | 251.4 |
| KDE bw = 25 µm, repaired | PRIMARY | 280.7 *(carried)* | 400.4 | +42.6 % | — | n/a (never halves) |
| KDE bw = 40 µm, repaired | PRIMARY | 277.4 *(carried)* | 396.3 | +42.9 % | — | n/a (never halves) |
| KDE bw = 60 µm, repaired | PRIMARY | 273.2 *(carried)* | 391.4 | +43.3 % | — | n/a (never halves) |
| KDE bw = 25 µm, native | PRIMARY | 282.6 *(carried)* | 393.8 | +39.3 % | — | n/a (never halves) |
| KDE bw = 40 µm, native | PRIMARY | 279.1 *(carried)* | 389.7 | +39.6 % | — | n/a (never halves) |
| KDE bw = 60 µm, native | PRIMARY | 274.8 *(carried)* | 384.8 | +40.0 % | — | n/a (never halves) |
| **isolated single vessel** (no distributed source) | — | **293.23** | **293.24–293.25** | **+0.01 %** | 203.20 | 203.20 |
| near-field pure exponential, primary arm | PRIMARY | 276.2 | 394.2 | +42.7 % | — | — |

**Read the two extremes of that table together.** On real tissue, with a distributed source, temozolomide's fitted decay length rises by **38–43 %** on every dx = 12 µm configuration (and by 123 % on the poorly-resolved dx = 30 µm arm). On the isolated single vessel, with the distributed source removed and everything else held identical, it moves by **+0.01 %**. The delivery-law change acts on the *distributed* source and essentially nothing else — which is both the mechanism of the reversal and the reason to believe it is physics rather than a bug.

**And the column that matters most is the last one.** Under exchange ON, temozolomide's empirical half-distance is absent on the PRIMARY distance metric in **every one of the nine tissue runs** measured for it (3 arms + 3 bandwidths × 2 masks). Under exchange OFF it was present on two of the three arms (the bandwidth sweep reports only fitted `L`, in both arms, so the OFF d½ column is left blank there rather than guessed). **Under the corrected delivery law there is no temozolomide penetration depth to report from this section.**

---

## Honest caveats

- **The delivery law itself is undecided, and it changes the headline temozolomide result.** Everything above is reported in two arms because `PDEConfig.transvascular_exchange` is currently `False` (a one-way leak known to be physically wrong) and the correct gradient-driven form is available but not the default. The choice is not cosmetic here: it decides whether this report can quote a temozolomide penetration depth at all. **Do not quote any temozolomide number from this file without its arm label.**
- **The exchange-ON numbers inherit an uncalibrated scale.** `exchange_k` is `vascular_density · bbb_penetrance`, a dimensionless product used as a rate in s⁻¹ via an implicit, underived conversion factor of exactly 1 s⁻¹ (`results/transvascular_exchange_revalidation.md` §6.2). The *form* of the ON law is right and its influx and back-flux use the same `PS`, so the qualitative conclusion (exchange flattens the profile, temozolomide loses its half-crossing) is robust; but the ON magnitudes — including the 396.3 µm — would move under a calibrated `k₀`. The OFF arm has no such term and does not inherit this.
- **The isolated single-vessel control's ON arm required choosing a normalisation the control does not define.** Its source is hand-painted with peak 1.0 and its `vascular_density` is zero, so `Cp` — and therefore `PS = source/Cp` — is undefined. Both readings were run and both give the same invariant `L`, which is why the conclusion is stated as robust; but it is a reconstructed `PS`, not one the module itself produces.
- **The nominally-converged grid breaks the transport domain.** At dx = 12 µm the cell-derived tissue mask is 26,058 disconnected fragments (largest = 0.7% of tissue), so lateral diffusion is switched off and *both* drugs collapse onto the vascular source field. `results/numerical_rigor.md` converged `kill_fraction` and `penetration_index` at dx ≲ 12 µm but did not test mask connectivity; this module does, and the primary arm therefore repairs the mask morphologically. **Any perivascular-gradient number quoted from a native dx ≤ 20 µm run is a source-geometry measurement, not a transport measurement.**
- **The repair is a modelling choice.** Closing the mask adds interstitial pixels the rasteriser missed; it is defensible (real tissue is continuous) but it is not measured, and the closing radius is a free parameter. The dx = 30 µm native-mask arm is included precisely because it needs no morphology at all, and the two agree only qualitatively — hence the honest ranges quoted above rather than a single number.
- **Part of the gradient is still vessel geometry, not drug transport.** The source null alone gives d½ = 10.0 µm. The engine delivers drug proportionally to the *continuous* density everywhere, so tissue scored as 'far from a vessel' still has a local source, which flattens the predicted profile relative to a true line-source geometry.
- **Neither metric measures a CD31-stained vessel lumen.** The PRIMARY metric is a direct distance to the position of a cell *typed* as endothelial from the Xenium marker panel (`reconstruct.real_fields.vascular_density_from_cells`'s input cells) — a real measured point, but a cell centroid, not a segmented lumen boundary, so distances are systematically offset from a true CD31 edge by roughly a cell radius. This is a smaller, more defensible gap than the LEGACY metric's: that one measures distance from a 40 µm-bandwidth Gaussian KDE ridge of the same cell positions (or, under the discrete field, a painted vessel-radius footprint), and the 40 µm bandwidth is the same order as the gradient being measured, putting a floor on how sharp the near-vessel profile can be. Switching to the direct-cell metric removes the KDE-bandwidth part of this caveat entirely; it does not remove the cell-centroid-vs-lumen part.
- **[LEGACY METRIC ONLY]** The percentile threshold is a free parameter, fixed *a priori* at the 95th percentile (≈5% vessel area fraction) rather than chosen to match the literature; the sensitivity table above shows what moves with it. This does not apply to the PRIMARY metric, which has no threshold.
- **Doxorubicin agreement is PARTLY BY CONSTRUCTION.** The single α in `validation.calibration` was fit to the *same* Primeau anchor, through the 1-D analytic relation `L = √(D_eff/k_eff)`. The isolated single-vessel control therefore recovers ~45 µm essentially by construction and is a **consistency check that the calibration survives the full numerical solver on real geometry, not an independent validation**. What is genuinely non-trivial: the shape comes out exponential, and the real heterogeneous 2-D geometry did not move the number by an order of magnitude.
- **The α convention matters and has changed.** This run used **α = 148.6**, calibrated against the *e-folding* target `45/ln2 ≈ 64.9 µm` (`validation/calibration.py`, `DOX_TARGET_EFOLD_UM`). An earlier convention fitted α to the raw 45 µm half-distance as if it were an e-folding length, giving α = 309.6 and a ~1.4× steeper gradient. Every absolute penetration number in this report inherits that single convention choice; check `calibrate_alpha()` before comparing to older results files.
- **Temozolomide is the independent output** — α was never fit to it — but there is **no published temozolomide perivascular-distance dataset**, so the temozolomide result is a *prediction*, checkable only for direction and order of magnitude. Under the corrected delivery law it is not even checkable for magnitude on this section, because it has no measurable half-distance. (The measured temozolomide anchor in `pk_anchors.py`, PMID 19861433, is a brain-ECF:plasma **partition ratio** — a different physical quantity that must not be compared to a penetration length.)
- **2D, straight-line distances.** The section is a 2D plane through a 3D vascular network: a vessel just above or below the plane is invisible, so grid points scored 'far' may really be near an out-of-plane vessel. Distances are Euclidean, not geodesic around the irregular tissue mask.
- **Not time-matched.** Primeau's images are ~10 min post-injection snapshots; the primary profile integrates a 6 h simulated exposure (AUC). The final-concentration readout is the closer snapshot analogue and agrees closely.
- **Untreated section, no measured drug field.** The Xenium GBM sample is an untreated resection; the comparison is to published measurements in other (murine / xenograft) tumors, not to drug imaged on *this* tissue.
