# Does any per-object vessel-radius model generalise the vascular area fraction across sections? A genuine null.

**Scope, precisely.** `results/heldout_patient.md` §2.1 already showed that `vessel_geometry.py`'s data-derived radius (the fix documented in `results/vessel_geometry_fix.md`, which replaced the older per-endothelial-cell 7.5 µm painting) still misses the measured CD34 vascular area fraction (2.1 ± 0.3 %, Angelucci/Tamma C et al., *Cells* 2022;12(1):11, PMID 36611806) by 0.59× on the development section and 5.23× on the held-out second patient. This is **not** the already-closed spatial-clustering question (`results/vascular_validation.md` §6a, `results/sweep/vessel_object_model_2026-07-29T19-07-08Z.md`) — that asked whether vessels sit CSR-randomly. This asks a narrower question: **is there a non-circular, section-invariant function from a vessel object's own endothelial-cell data (`k`, its cells-per-cross-section) to its radius that reproduces the measured area fraction on more than one section?**

**Bottom line up front: no.** Six candidate radius functions (the two already built into `vessel_geometry.py`, plus four new ones motivated by the specific literature caveat the module's own docstring already flags) were tested on all three real GBM sections this repo has cached. Every one trades the two out-of-scope sections off against the development section in the **same direction**: whatever moves the held-out patient toward 2.1% moves the development section further **away** from it, because the two errors run in opposite signs (development section is 0.59× too low; both other sections are 5.2× too high) while every k-based radius model is monotonic in the same direction on all three. A direct search for the constant radius that *would* hit 2.1% on each section separately finds that the required radius is 7.28 µm (section 1), 3.12 µm (section 2) and 2.57 µm (section 3) — a perfect inverse rank correlation with each section's own (already measured, already known-defective) microvessel density (110.8 / 233.3 / 320.4 vessels/mm²). That inverse relationship is directionally consistent with a real, published GBM finding (MVD inversely correlated with vessel-size index — Halshtok Neiman et al., *Clin Neuroradiol* 2021, PMID 34940899), but operationalising it here would mean calibrating the vessel radius against this project's own vessel *count*, which `results/heldout_patient.md` §6.3 and `results/vessel_geometry_fix.md` §5 already identify, independently, as an unresolved, both-sections-fail segmentation artefact (MVD is 1.7–4.9× the published cohort range on every section tried). Conditioning one already-wrong number on another would not be a fix; it would be laundering one open defect through a second one. No code in `tmesim/` was changed.

## 1. What was already known, read first (not re-litigated)

* `tmesim/reconstruct/vessel_geometry.py`'s own docstring already documents that its `radius_mode='data'` step (`r = r₁·k`, circumferential tiling) is validated at `k=2` by an independent cross-check on the sections' own two-cell-object nuclear separations (implying `r(k=2) ≈ 5.7–6.4 µm` against the model's `6.23 µm`), **and** flags a genuine tension: Sargent SM et al. (*Vasc Biol* 2023;5:e230010, PMID 37582180) separately reports that *within the capillary zone*, two-EC capillaries are only ~0.4 µm wider (diameter) than one-EC capillaries — far less than strict tiling predicts — and states plainly that EC number "alone is insufficient to explain the full range of capillary diameters." The module already calls this "the single largest source of model uncertainty," not a solved problem.
* `results/vessel_geometry_fix.md` §5 already swept `radius_mode ∈ {data, data-k≤6, stereological}` on three sections (at slightly different endothelial-typing numbers than reported here — the marker-gating pipeline has moved on since that report, see §5 caveat below) and concluded **"no model dominates"**, attributing the residual disagreement explicitly to "segmentation, not mesh": endothelial-cell density (hence MVD, hence `k`) disagrees across sections by 4–27× against the human-EM ground truth, a defect in the typing/linking stage that this fix does not and should not paper over.
* This report's job was to check whether that "no model dominates" verdict still holds on the **current** pipeline (marker maps and typing have changed since that fix was written — the numbers below do not match `vessel_geometry_fix.md`'s table, e.g. section 2 is now 10.98% where that report had 7.14%), to test literature-motivated candidates that report did not try, and to look for an actual quantitative external anchor rather than assert one.

## 2. Hypothesis, stated before looking at the fix candidates

If the `r ∝ k` over-dispersion Sargent flags is the generalisation failure, then **flattening** the radius's dependence on `k` (sub-linear, capped, or literature-bounded) should shrink the area fraction more on sections with more multi-EC objects — and the held-out section (mean `k`=2.32) has far more multi-EC objects than the development section (mean `k`=1.36) — so flattening should close the gap on the held-out section *without* opening a new one on the development section, provided the development section's own `k` distribution is not itself already stretched thin by the same defect. This is falsifiable: if flattening `k` moves both sections in the *same wrong direction*, the defect is not `r(k)`'s functional form.

## 3. Method

All three real GBM Xenium sections cached in this repo were driven through the frozen production path — `reconstruct_real(source='xenium', sample=..., dx_um=12.0)` → `segment_vessels(..., link_distance_um=LINK_DISTANCE_UM)` (20 µm, unchanged) — exactly as `tmesim.validation.heldout_patient` does for sections 1–2. Section 3 (`data_cache/gse324268_xenium`, GEO GSE324268 `GBMtumor07`) is not part of the `heldout_patient` harness but is already used elsewhere in this repo (`validation/cell_state_impact.py`, `validation/endothelial_density.py`, `results/vessel_geometry_fix.md`) and was added here purely to get a third independent data point — n=2 cannot distinguish "this generalises" from "these two tumours differ," and n=3 barely can, but it is what is available.

For each section, the production `vessel_geometry.vessel_geometry()` / `vascular_area_fraction()` / `_subgrid_coverage()` were called **unmodified** for the two radius modes already built into the module (`'data'`, `'stereological'`, plus the already-supported `max_ec_per_cross_section` cap). Four further candidates were prototyped by constructing an alternative `geom` dict with a different `k → r` map and feeding it to the **same, unmodified** `_subgrid_coverage` painter, so only the radius assignment varies and the sub-grid rasterisation, tissue masking and area accounting are identical to production in every candidate:

| candidate | formula | rationale |
|---|---|---|
| C0 baseline | `r = r₁·k`, no cap | current production default (`radius_mode='data'`) |
| C1 stereological | `r = r̄ = 4.669 µm` for every object | current production alternative; no `k` dependence at all |
| C2 linear, cap k≤6 | `r = r₁·min(k,6)` | already-supported `max_ec_per_cross_section=6.0`, Sargent's observed ceiling |
| C3 linear, cap k≤2 | `r = r₁·min(k,2)` | restricts tiling to the true modal capillary range (1–2 EC) Sargent reports |
| C4 sqrt(k) | `r = r₁·√k` | area (not circumference) scales with EC count — a distinct, equally defensible physical assumption |
| C5 sqrt(k), cap k≤6 | `r = r₁·√min(k,6)` | C4 + the literature ceiling |
| C6 piecewise, literature-anchored | flat `+0.2 µm/EC` for `k≤2` (Sargent's measured ~0.4 µm diameter delta), resumes full tiling slope for `2<k≤6`, frozen beyond `k=6` | the most literally Sargent-faithful curve: flat where Sargent measured flat, tiling only where Sargent did not contradict it, no extrapolation past the observed range |

## 4. Results — area fraction (%), all three sections, against the 2.1 ± 0.3% target

| candidate | section 1 (dev) | section 2 (held-out) | section 3 (gse324268) |
|---|---:|---:|---:|
| **measured target** | **2.1** | **2.1** | **2.1** |
| C0 baseline (`data`, production default) | 1.246% (0.59×) | 10.982% (5.23×) | 10.938% (5.21×) |
| C1 stereological | 1.031% (0.49×) | 3.476% (1.66×) | 4.353% (2.07×) |
| C2 linear, cap k≤6 | 1.150% (0.55×) | 8.027% (3.82×) | 8.854% (4.22×) |
| C3 linear, cap k≤2 | 0.804% (0.38×) | 3.696% (1.76×) | 4.614% (2.20×) |
| C4 sqrt(k) | 0.775% (0.37×) | 4.451% (2.12×) | 5.007% (2.38×) |
| C5 sqrt(k), cap k≤6 | 0.759% (0.36×) | 3.946% (1.88×) | 4.638% (2.21×) |
| C6 piecewise, literature-anchored | 0.928% (0.44×) | 6.587% (3.14×) | 7.090% (3.38×) |

(C0's 1.246%/10.982% reproduce `results/heldout_patient.md`'s reported 1.25%/10.98% to three figures — this confirms the prototype's plumbing matches production exactly, not just approximately.)

**Every single candidate that closes the gap on section 2 and/or 3 (C1, C3, C4, C5) pulls section 1 further below target, not closer to it.** Section 1 starts at 0.59× and every flattening candidate pushes it lower still (down to 0.36–0.44×), because section 1's own `k` distribution (mean 1.36, only 7.3% of objects above `k=2`) already sits close to the capillary floor these candidates flatten toward — there is very little `k`-inflation on section 1 to remove. The candidates that help the two more-cellular sections do so specifically by shrinking the radius attached to their comparatively larger population of multi-EC objects (24.0% above `k=2` on section 2, 22.7% on section 3), which section 1 mostly does not have.

## 5. Is there a constant (or any monotone function of `k`) that could work? A direct search says no

For each section, the single per-object constant radius that makes `vascular_area_fraction` hit exactly 2.1% was found by bisection (using `radius_mode='constant'`, so this is not a new formula — it is asking "what would even the simplest possible model need"):

| section | tissue area (mm²) | MVD (vessels/mm²) | radius that hits 2.1% exactly (µm) |
|---|---:|---:|---:|
| 1 (dev) | 27.9 | 110.8 | **7.280** |
| 2 (held-out) | 89.5 | 233.3 | **3.119** |
| 3 (gse324268) | 98.5 | 320.4 | **2.574** |

The required radius falls monotonically as MVD rises — a **perfect inverse rank correlation (Spearman = −1.00, n=3)**. Applying any one section's "correct" radius to the others is not a near-miss, it is categorically wrong in both directions:

| radius calibrated on ↓ / applied to → | section 1 | section 2 | section 3 |
|---|---:|---:|---:|
| section 1's r*=7.28 µm | 2.13% (✓ by construction) | 6.28% (3.0×) | 7.92% (3.8×) |
| section 2's r*=3.12 µm | 0.55% (0.26×) | 2.10% (✓) | 2.63% (1.25×) |
| section 3's r*=2.57 µm | 0.42% (0.20×) | 1.68% (0.80×) | 2.10% (✓) |

No off-diagonal cell is within the published ±0.3 percentage-point tolerance. There is no radius, constant or `k`-dependent, that is simultaneously right for a sparsely-vascularised section and a densely-vascularised one, because what the data actually calls for is a radius that shrinks precisely where the (already known to be over-counted) vessel density rises.

## 6. Why this is not "found the fix, just condition on MVD" — the circularity, checked directly

The inverse MVD-vs-radius relationship found in §5 is not invented for this report: Halshtok Neiman et al. (*Clin Neuroradiol* 2021, PMID 34940899) find, in human cerebral neoplasms including GBM, that MRI vessel-size index rises with tumour volume (r=+0.502) while microvessel density falls with it (r=−0.531) — i.e. the same qualitative trade-off (denser vasculature, smaller-caliber vessels) this project's own three sections independently reproduce at r=−1.00. That direction is genuine tumour biology, not an artefact of this search. Two things nonetheless rule out building a fix on it here:

1. **This project's own MVD is not a trustworthy input.** `results/heldout_patient.md` §6.3 already calls MVD (110.8 / 233.3 vessels/mm² here, 65 [35–95] published) "a systematic property of the segmentation [over-splitting] rather than something the held-out patient revealed" — it fails on **both** already-tested sections, identically, for reasons unrelated to vessel geometry (the 20 µm linking distance and the endothelial marker gate, per `vessel_geometry_fix.md` §5's own diagnosis). Scaling the radius by a number that is itself already flagged as wrong, in an unknown and section-specific way, does not correct the area fraction — it recombines two uncontrolled errors and calls the result fixed. Whether the resulting area fraction would land near 2.1% by genuine cancellation or by coincidence cannot be told apart with the data in hand.
2. **n=3 cannot fit and validate a calibration in the same breath.** Any function `r(MVD)` capable of hitting three arbitrary targets exactly needs at least 2 free parameters against 3 data points — 1 residual degree of freedom, which is not a held-out test, it is curve-fitting with almost no slack. There is no fourth section available to check the fitted relationship against, which is exactly the overfitting this task was warned against forcing.

## 7. Is this fundamentally underdetermined? Yes — and here is the concrete missing signal

Xenium centroids give exactly two per-object numbers: how many endothelial nuclei fell in a segmented object (`N`), and that object's own 2D spatial extent (`E`). Both `k` (`= N / max(1, E/50µm)`) and every candidate radius above are functions of only those two numbers. What none of them carries is an actual **caliber signal** — a lumen boundary, a basement-membrane / perimeter marker, or mural-cell coverage, any of which would let radius be measured rather than inferred from nuclear count. Concretely, the current codebase already has a candidate for exactly that signal and has already declined to use it for an unrelated, structural reason: `tmesim/reconstruct/marker_panels.py` records that the 480-gene IO panel (section 2) measures the pericyte/mural programme directly (`RGS5` 12.5%, `PDGFRB` 9.5%, `NOTCH3` 24.9%, `ACTA2` 38.0% of cells) — pericyte/smooth-muscle coverage is a real, independent correlate of vessel caliber — but `CellType` has no `MURAL` member, so those markers are folded into `STROMAL`/`TUMOR` rather than kept as a vessel-wall signal, and (per that module's own note) is not on the 319-gene development panel at all, so it could not be a cross-section-consistent radius input even if `CellType` were extended. Absent a perimeter/caliber signal that (a) exists on every section and (b) is not itself a per-section-varying artefact, the vascular area fraction is underdetermined by 2D endothelial-cell-position data alone, exactly as the module's own docstring already half-admits ("the single largest source of model uncertainty").

## 8. Verdict

**Genuine, honest null.** No `k`-dependent or constant per-object radius model — six tried here, two already in production — generalises the vascular area fraction across more than one of the three real GBM sections this repo has cached, because the development section's error (0.59×, too low) and the other two sections' error (5.2×, too high) require the radius to move in opposite directions, and no section-invariant function of a vessel object's own endothelial-cell data does that. The one relationship that WOULD reconcile all three (radius inversely tied to each section's own MVD) is directionally supported by independent published GBM biology, but this project's own MVD is itself an unresolved, already-documented, both-sections-fail defect — conditioning a "fix" on it would be circular, not corrective, and there are not enough independent sections (n=3, need ≥1 held out from any fit) to validate such a calibration honestly. The path to an actual fix is a real caliber/perimeter signal (lumen boundary, basement-membrane marker, or a `MURAL` cell type built from the pericyte markers the IO panel already measures) that does not currently exist across all sections in a comparable form — not a better function of `k`.

**No code was changed.** `tmesim/reconstruct/vessel_geometry.py` and `tmesim/reconstruct/vessel_network.py` remain as they were; `results/heldout_patient.md`, `results/vessel_network_impact.md` and `results/vessel_geometry_fix.md` were read but not regenerated, because nothing they depend on moved.

## 9. Caveats

1. **n=3 sections, all GBM, all Xenium.** This cannot separate "no k-based radius model can ever work" from "these three tumours/panels happen to need opposite corrections" — but the internal consistency of the MVD-vs-required-radius trend across all three (not just the original two) is a stronger signal than n=2 alone would give.
2. **`vessel_geometry_fix.md`'s own per-section numbers are now stale.** Its §5 table (section 2 = 7.14% under `radius_mode='data'`) does not match this report's 10.98% (matching the current `heldout_patient.md`), because the marker-gating/typing pipeline has changed since that fix was written. This report's numbers are the current, live ones; `vessel_geometry_fix.md`'s qualitative conclusion ("no model dominates") is reconfirmed here on the current pipeline and with a third section added, but its specific numeric table should not be read as current.
3. **The MVD-vs-vessel-size literature citation (PMID 34940899) was read via search-result abstracts/snippets, not the full paper text**, so the correlation coefficients quoted (r=+0.502 tumour-volume-vs-vessel-size, r=−0.531 tumour-volume-vs-MVD) are secondhand from the abstract and should be verified against the primary text before being cited as a hard number anywhere load-bearing; here they are used only to support a qualitative direction, not a fitted slope.
4. **All candidates were scored against the vessel objects `segment_vessels` already produced at the frozen 20 µm linking distance** — i.e., this report holds the (separately, already-documented) MVD over-segmentation defect fixed and asks only whether the per-object radius can compensate for it. It cannot, and it should not be asked to: the two are genuinely separate defects at separate pipeline stages.

## 10. Reproduction

Standalone prototype, not committed under `tmesim/` (per this task's instructions — no `tmesim/` file was edited). Calls `reconstruct_real`, `segment_vessels`, `vessel_geometry.vessel_geometry`, `vessel_geometry._subgrid_coverage` and `vessel_geometry.vascular_area_fraction` verbatim; the only new code constructs alternative `geom` dicts (radius arrays) fed into the unmodified sub-grid painter.

```python
from tmesim.reconstruct.real_pipeline import reconstruct_real
from tmesim.reconstruct.vessel_network import LINK_DISTANCE_UM, segment_vessels
from tmesim.reconstruct.vessel_geometry import vessel_geometry, _subgrid_coverage

for sample in ("data_cache/xenium_gbm", "data_cache/xenium_gbm_io", "data_cache/gse324268_xenium"):
    micro = reconstruct_real(source="xenium", sample=sample, dx_um=12.0)
    vessels = segment_vessels(micro.cells, micro.domain, link_distance_um=LINK_DISTANCE_UM)
    geom = vessel_geometry(vessels, radius_mode="data")   # or 'stereological', or a custom k->r map
    cov, cell_area = _subgrid_coverage(micro.domain, geom)
    # area_fraction_pct = 100 * (cov * cell_area)[tissue].sum() / cell_area[tissue].sum()
```

Runtime: ~2 s (section 1, 40,887 cells) to ~30 s (sections 2–3, ~820–850k cells each) per section for reconstruction + segmentation; no PDE solve is needed for this check.
