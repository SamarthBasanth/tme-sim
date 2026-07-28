# Vessel geometry fix: a grid-independent, data-derived vascular area fraction

**Scope.** Vessel geometry only. Nothing in `tmesim/core`, `tmesim/physics`,
`tmesim/sim/engine.py` or any drug / PK / response-validation code was touched.
New module: `tmesim/reconstruct/vessel_geometry.py`. Modified:
`tmesim/reconstruct/vessel_network.py` (docstring + `segmentation_report` +
an opt-in `geometry=` keyword), `tests/test_vessel_network.py`.
Suite: **162 passed**.

**Sections.** `xenium_319` = `data_cache/xenium_gbm` (319-gene brain panel,
"section 1"); `xenium_io` = `data_cache/xenium_gbm_io` (480-gene IO panel,
"section 2"); `gse324268` = GEO GSE324268 `GBMtumor07` ("section 3", loaded
successfully). Segmentation is `segment_vessels` at the unchanged 20 µm linking
distance throughout; the vessel objects are identical before and after the fix.

**What the acceptance test is.** (a) **Grid-independence** — a quantity that
claims to describe tissue must not move when the mesh does. (b) **Fidelity to the
sample** — the field must reproduce *this* tumour's measured endothelial
distribution, including its unevenness. It is **not** "hit 2.1 %" or "land inside
35–95 vessels/mm²". Those are cohort means over tens of patients; glioblastoma
vasculature is abnormal and heterogeneous by definition, and that abnormality is
the disease, not an error. They are reported below as context only.

---

## 1. Reproducing the reported defect — and correcting its diagnosis

The brief was that the area fraction is "grid-imposed": at the production
`dx = 12 µm` one grid cell has an equivalent-disc radius of 6.8 µm, so the mesh
cannot represent the hardcoded 7.5 µm vessel, and the quantity therefore tracks
the mesh. **The first half of that is true; the conclusion is not.** Sweeping the
production sections over `dx = 30/20/12/8/6 µm` under the *unmodified* fixed-disc
implementation:

| dx (µm) | cell equiv. radius | s1 `xenium_319` | s2 `xenium_io` | s3 `gse324268` |
|---:|---:|---:|---:|---:|
| 30 | 16.93 | 2.563 % | 5.783 % | 17.047 % |
| 20 | 11.28 | 2.503 % | 5.811 % | 16.862 % |
| 12 |  6.77 | **2.499 %** | **5.791 %** | 16.865 % |
|  8 |  4.51 | 2.505 % | 5.775 % | 16.875 % |
|  6 |  3.39 | 2.497 % | 5.776 % | 16.893 % |
| **max/min** | | **1.027×** | **1.006×** | **1.011×** |

The area fraction barely moves. It was already inside the project's 1.05×
convention. The reason is a fact about lattices, not luck: point-sampling a fixed
physical region on a regular grid is an **unbiased** estimator of that region's
area at *any* pitch — the expected number of lattice centres inside a region of
area `A` is `A/dx²` whatever `dx` is. Coarsening the mesh adds variance, not bias.

**What *is* violently grid-dependent is the field's geometry**, and that is a real
defect:

| dx (µm) | painted pixels per vessel object (s1 / s2 / s3) | heterogeneity ratio, binary field (s1 / s2 / s3) |
|---:|---|---|
| 30 | 0.26 / 0.62 / 0.53 | 1.693 / 1.239 / 1.319 |
| 20 | 0.56 / 1.40 / 1.18 | 1.257 / 1.041 / 1.058 |
| 12 | 1.57 / 3.86 / 3.27 | 1.017 / 0.965 / 0.952 |
|  8 | 3.53 / 8.67 / 7.37 | 0.993 / 0.953 / 0.932 |
|  6 | 6.26 / 15.41 / 13.11 | 0.971 / 0.948 / 0.920 |
| **max/min** | **24.4× / 25.0× / 24.8×** | **1.744× / 1.307× / 1.433×** |

At `dx = 30 µm` **three quarters of section 1's vessels paint no pixel at all** —
they are absent from the delivery source. And the *unevenness* of the field
(§4) moves by up to 1.74× with the mesh, which for a patient-specific model is
the more damaging of the two.

### The defect that actually moves the number

The legacy painter puts a disc around **every endothelial cell**, so the painted
area is essentially `n_EC × π r²`. The area fraction therefore measures
**endothelial-cell density** — a property of the marker panel and the typing
threshold — rather than vessel density:

| | s1 | s2 | s3 |
|---|---:|---:|---:|
| tissue area | 27.9 mm² | 89.5 mm² | 98.5 mm² |
| endothelial cells | 3,984 (143 /mm²) | 34,179 (**382 /mm²**) | 103,817 (**1054 /mm²**) |
| vessel objects | 3,092 (**111 /mm²**) | 9,322 (**104 /mm²**) | 35,243 (358 /mm²) |
| endothelial cells per object | 1.29 | 3.67 | 2.95 |
| singleton objects | 84.8 % | 42.9 % | 52.8 % |
| legacy area fraction @ dx=12 | 2.50 % | 5.79 % | 16.87 % |

Sections 1 and 2 have **near-identical vessel densities (111 vs 104 /mm²)** but
**2.7× different endothelial densities** — and their reported area fractions
differed by exactly that ratio. A vessel profile lined by `k` endothelial nuclei
is one vessel, not `k` vessels. Section 1's previous "1.19× agreement" with 2.1 %
was the coincidence of `143 EC/mm² × π(7.5 µm)² ≈ 2.5 %`.

The second problem is that the 7.5 µm constant was justified partly *by the mesh*
— its own docstring called it "the smallest radius that reliably paints ≥1 grid
cell per vessel". A parameter chosen from the grid cannot be used to produce a
tissue measurement, whether or not the resulting number happens to be stable.

---

## 2. The fix — both (a) and (b), because they solve different halves

### (a) Sub-grid fractional-area painting

`vessel_geometry.subgrid_coverage_field` rasterises each footprint onto a lattice
whose pitch is a fixed **physical** length (`clip(r_min/4, 0.25, 1.0)` µm, set by
the vessels and never by `dx`), then bins the sub-pixels into the domain grid. A
cell returns the *fraction of its area* that is vessel. Consequences:

* a vessel far smaller than a grid cell deposits exactly its own area — the "no
  vessel below 6.8 µm is representable" limitation is gone;
* the painted total is computed **with no reference to `dx`**, so it is
  grid-independent by construction rather than by unbiasedness-on-average;
* the Monte-Carlo variance that produced the residual 1.03× disappears.

Two implementation details are load-bearing and are pinned by tests. Stencil
radii are grouped **logarithmically at 1 %** (a *relative* tolerance — absolute
rounding to whole sub-pixels costs half a sub-pixel and is worst on the smallest
vessels, which are the majority). And the **sub-pixel phase** of small discs is
kept rather than snapped: snapping is not merely imprecise, it is *biased*,
because every disc of a given radius then receives the identical lattice count and
the O(1) point-count error stops averaging out over the section. Measured against
a closed form on non-lattice-aligned discs, the painter is accurate to **≤1.2 %**
over `r = 3–20 µm` with a dx-spread of ≤1.0035×.

### (b) A per-vessel radius derived from the data

`vessel_geometry.vessel_geometry`. A segmented object's `N` endothelial cells are
distributed two ways: **around** the wall (`k` per cross-section) and **along** the
axis. The object's own in-plane extent `E` says how much of the split is axial:

```
k = N / max(1, E / 50 µm)          r = r₁ · k
footprint = everything within r of an axis segment of length max(0, E − 2r)
```

The two limits are the two things a thin section can show, and the formula
reproduces both without a switch: a **transverse cut** has `E ≈ 2r`, the axis
collapses to a point and the footprint is a single disc — *one* profile, not `N`;
an **in-plane run** has `E ≫ 2r`, `k` falls towards 1, and the footprint is a
stadium along the run.

**Every constant is measured; none is fitted, and none was chosen to hit 2.1 %.**

| constant | value | provenance |
|---|---:|---|
| endothelial cells per mm of vessel | 20 | **Human** cortex, H01 petavoxel EM: 4,604 endothelial cells lining ~230 mm of vasculature, "~20 per mm". Shapson-Coe 2024, PMID 38723085, [doi:10.1126/science.adk4858](https://doi.org/10.1126/science.adk4858) |
| ⇒ axial pitch, one cell | 50 µm | reciprocal of the above |
| vessel length density | 219 mm/mm³ | same dataset: 230 mm in 1.05 mm³ |
| vascular volume fraction | 1.5 % | same dataset: "collapsed blood vessels ~1.5 %" of neuropil volume |
| ⇒ **mean vessel radius r̄** | **4.669 µm** | *derived*: `VV = L_V·πr̄²` ⟹ `r̄ = sqrt(VV/(πL_V))`. No free parameter. |
| mean cells per cross-section k̄ | 1.5 | capillaries are "either one or two endothelial cells, with roughly equal proportions"; transition vessels two to six. Sargent 2023, PMID 37582180, [doi:10.1530/VB-23-0010](https://doi.org/10.1530/VB-23-0010) |
| ⇒ **r₁** (one cell per cross-section) | **3.113 µm** | *derived*: `k` cells tile the circumference so `r ∝ k`, hence `r₁ = r̄/k̄` |

Using one dataset (H01) for both the length density and the volume fraction is
deliberate: it makes the two mutually consistent instead of mixing cohorts.

**Two independent cross-checks on r₁,** neither used to set it:
* Sargent 2023 measures a one-cell capillary **lumen** diameter of ~2–5 µm, i.e. a
  lumen radius of 1.0–2.5 µm plus wall — the same scale as 3.11 µm.
* From the sections themselves: a two-cell profile is two cells interlocked around
  one lumen, so their nuclear centroids sit across a diameter. The measured median
  separation of two-cell objects is **11.3 µm** (s2) and **12.7 µm** (s3), implying
  `r(k=2) ≈ 5.7–6.4 µm` against this model's `2r₁ = 6.23 µm`. (Section 1's 16.8 µm
  is truncated against the 20 µm linking distance and is not usable.)

**Documented limitation of the `r ∝ k` step.** Sargent 2023 also reports that
*within* the capillary zone a two-cell capillary is only ~0.4 µm wider than a
one-cell capillary — far less than the factor 2 strict circumferential tiling
predicts — and states explicitly that endothelial number "alone is insufficient to
explain the full range of capillary diameters". `r ∝ k` therefore **over-disperses
radii at the capillary end**. This is the single largest source of model
uncertainty in the numbers below, and it is bracketed in §5, not hidden.

---

## 3. Grid-independence after the fix — the acceptance test

`vascular_area_fraction`, `radius_mode='data'` (the default), same sections, same
sweep:

| dx (µm) | s1 `xenium_319` | s2 `xenium_io` | s3 `gse324268` |
|---:|---:|---:|---:|
| 30 | 0.9916 % | 7.1409 % | 17.7868 % |
| 20 | 0.9906 % | 7.1399 % | 17.7908 % |
| 12 | 0.9908 % | 7.1425 % | 17.7906 % |
|  8 | 0.9908 % | 7.1429 % | 17.7908 % |
|  6 | 0.9907 % | 7.1425 % | 17.7904 % |
| **max/min** | **1.0010×** | **1.0004×** | **1.0002×** |
| *before (legacy)* | *1.027×* | *1.006×* | *1.011×* |

**PASS**, with three orders of magnitude of margin against the 1.05× convention
(≤0.1 % versus a 5 % allowance). The residual is the assignment of sub-pixels to
grid cells near the tissue-mask boundary; the painted *area* itself does not
depend on `dx` at all.

The **heterogeneity** of the field is grid-independent too, which the legacy field
was not:

| | s1 | s2 | s3 |
|---|---:|---:|---:|
| legacy binary, het. ratio spread over dx | 1.744× | 1.307× | 1.433× |
| **sub-grid, het. ratio spread over dx** | **1.014×** | **1.009×** | **1.024×** |

Pinned by `tests/test_vessel_network.py::test_area_fraction_is_grid_independent`
(and `..._preserves_measured_heterogeneity...` for the second table), which assert
≤1.05× across `dx = 30/20/12/8/6` and, for the area fraction, <1.005×.

---

## 4. Fidelity to the sample

Does the field put vasculature where the endothelial cells actually are, and leave
avascular tissue avascular? Endothelial cells and field mass are aggregated into
fixed 120/100 µm **physical** bins and compared. `heterogeneity ratio` =
CV(field per bin) / CV(measured endothelial count per bin); 1.0 means the measured
unevenness is carried exactly, well below 1.0 means it has been smoothed away.
At `dx = 12 µm`:

| section | field | het. ratio | Pearson r | Spearman ρ | tissue bins with no vessel (measured) |
|---|---|---:|---:|---:|---|
| s1 | KDE (historical) | 0.696 | 0.878 | — | **0.0 %** (34.6 %) |
| s1 | legacy binary | 1.017 | 0.933 | — | 34.6 % (34.6 %) |
| s1 | **sub-grid (fix)** | 2.315 | 0.786 | 0.922 | 31.8 % (34.6 %) |
| s2 | KDE | 0.724 | 0.906 | — | **0.5 %** (30.5 %) |
| s2 | legacy binary | 0.965 | 0.967 | — | 30.3 % (30.5 %) |
| s2 | **sub-grid (fix)** | 1.294 | 0.931 | 0.938 | 26.0 % (30.5 %) |
| s3 | KDE | 0.859 | 0.921 | — | **0.0 %** (0.9 %) |
| s3 | legacy binary | 0.952 | 0.953 | — | 0.9 % (0.9 %) |
| s3 | **sub-grid (fix)** | 1.656 | 0.908 | 0.920 | 0.7 % (0.9 %) |

* The **KDE fills in every avascular bin** (0.0–0.5 % empty against a measured
  34.6 / 30.5 / 0.9 %) and flattens the variance to 0.70–0.86 of the measured.
  This re-confirms the defect `vessel_network.py` was written to address.
* The **sub-grid field keeps avascular tissue avascular** (within 3, 5 and 0.2
  percentage points; the gap is footprints that legitimately straddle a bin edge)
  and correlates strongly with the measured cell counts (ρ = 0.92–0.94).
* Its heterogeneity ratio is **above** 1, i.e. it is *sharper* than the raw
  endothelial counts. That is the `r ∝ k` scaling: a bin holding one four-nucleus
  profile is given area ∝ 16, not ∝ 4. It is a modelling claim with a stated
  mechanism, and it is the step Sargent cautions against — see §5.

---

## 5. Per-section area fraction, reported without tuning

At `dx = 12 µm`, `radius_mode='data'`. **The published values are population
context, not targets.**

| section | legacy (fixed 7.5 µm, per cell) | **fixed (sub-grid, data-derived)** | published cohort 2.1 ± 0.3 % |
|---|---:|---:|---|
| s1 `xenium_319` | 2.50 % | **0.99 %** | 0.47× the cohort mean |
| s2 `xenium_io` | 5.79 % | **7.14 %** | 3.40× |
| s3 `gse324268` | 16.87 % | **17.79 %** | 8.47× |

*(Published: CD34⁺ vascular area fraction, human IDH-wt GBM, n = 30 — Tamma R,
Ribatti D et al., Cells 2022;12(1):11, PMID 36611806,
[doi:10.3390/cells12010011](https://doi.org/10.3390/cells12010011).)*

### Radius-model sensitivity (all data/literature-derived; none fitted)

| model | s1 | s2 | s3 |
|---|---:|---:|---:|
| `data` — `r = r₁·k`, **default** | 0.99 % | 7.14 % | 17.79 % |
| `data`, capped at Sargent's observed max `k ≤ 6` | 0.94 % | 5.35 % | 13.36 % |
| `stereological` — every object at the measured r̄ = 4.67 µm | 0.97 % | 1.98 % | 5.75 % |
| legacy 7.5 µm radius, but painted sub-grid | 2.15 % | 3.45 % | 10.42 % |

Fidelity of each radius model (het. ratio / Pearson / Spearman, dx = 12):

| model | s1 | s2 | s3 |
|---|---|---|---|
| `data` | 2.32 / 0.79 / 0.92 | 1.29 / 0.93 / 0.94 | 1.66 / 0.91 / 0.92 |
| `data`, `k ≤ 6` | 2.05 / 0.81 / 0.92 | 1.16 / 0.92 / 0.95 | 1.31 / 0.80 / 0.88 |
| `stereological` | 0.97 / 0.96 / 0.94 | 0.95 / 0.92 / 0.96 | 0.85 / 0.74 / 0.86 |

No model dominates: `stereological` tracks measured heterogeneity most closely on
s1 and s2 but has the weakest rank agreement on s3; `data` has the best Pearson on
the densest section. `data` is kept as the default because it is the per-vessel,
data-derived radius the fix is *for* — `stereological` is a single constant, even
if a measured one — but the spread of this table (a factor 3.6 on s2) is the
honest uncertainty on the absolute area fraction, and it is dominated by the
`r ∝ k` step, not by anything to do with the mesh.

### What the residual disagreement is: segmentation, not mesh

The three sections now disagree with each other by up to 18×, and only ~1.3× of
the change on any section came from the geometry. **This is a finding about the
endothelial call rate, not about the grid.** The evidence:

1. **Endothelial density is implausible against human EM ground truth.** H01
   (PMID 38723085) measures 4,604 endothelial cells in 1.05 mm³ of human cortex =
   4,385 /mm³, i.e. ~22 /mm² in a 5 µm section (~39 /mm² allowing an Abercrombie
   correction for a ~4 µm nuclear diameter). The reconstructions report **143 /
   382 / 1054 per mm²** — 4× to 27× that, before allowing anything for GBM
   hypervascularity.
2. **The sections disagree far more in cells than in vessels.** Endothelial
   density varies 7.4× across the three sections; vessel-object density varies
   3.4×; and sections 1 and 2 have essentially the same vessel density (111 vs
   104 /mm²) with 2.7× different cell densities. The panel, not the vasculature,
   is what changed.
3. **Objects with implausible cells-per-cross-section.** The `k` distribution runs
   to a maximum of 12 / 26 / 62 against Sargent's observed ceiling of 5–6. An
   object with `k = 62` (radius 193 µm on s3) is single-linkage percolation
   merging a whole vascular network into one object, not a wide vessel. Section 3
   also has 358 objects/mm² with 52.8 % singletons — simultaneously over-merging
   in dense regions and over-splitting in sparse ones.
4. **Vessel axis length per unit area** — 0.093 / 0.425 / 1.303 mm per mm² —
   spans 14× across sections, again tracking the cell count.

The 20 µm linking distance and the endothelial marker gate are the parameters that
now carry the between-section discrepancy. That is where the next investigation
belongs; it is out of scope here and this fix does not paper over it.

**And the caveat that limits everything absolute here.** A 2D section measures a
volume fraction only through a section-sampling depth. `L_V = n_A/(t·λ)` is
linear in the effective depth `t_eff`, and `t_eff` is not known: at the nominal
5 µm section thickness section 1 implies 1,430 mm/mm³ of vessel against H01's
219 mm/mm³, while an Abercrombie-corrected ~15 µm implies 477. Any absolute area
fraction from Xenium centroids inherits that factor-of-three ambiguity. What is
*not* ambiguous, and what this fix delivers, is that the number no longer depends
on the mesh and no longer depends on a constant chosen from it.

---

## 6. Effect on MVD and inter-vessel distance: **none, by construction**

| quantity | s1 | s2 | s3 | published cohort (context) |
|---|---:|---:|---:|---|
| MVD (objects/mm²) | 110.9 | 104.1 | 357.8 | 65 [35–95] |
| spacing from MVD (µm) | 95.0 | 98.0 | 52.9 | 124 [102.6–169] |
| inter-vessel NN median (µm) | 40.1 | 31.0 | 24.8 | — |
| median tumour→vessel (µm) | 43.3 | 39.6 | 17.5 | 58 [48.2–79.4] |

**Unchanged to machine precision before and after the fix.** These are computed
from cell *coordinates* and the tissue area; the vessel radius never enters any of
them, so correcting the radius cannot move them and must not be credited with
doing so. Pinned by
`tests/test_vessel_network.py::test_mvd_and_intervessel_distance_do_not_use_the_radius`,
which re-runs `segmentation_report` under three different radius settings and
asserts these three rows are bit-identical.

So the MVD and spacing figures stand where they were. Under the rescoped
criterion they are **not** defects: they are what this tissue measures. Sections 1
and 2 sit above the published cohort MVD range and section 3 far above it, with
correspondingly short inter-vessel spacing — a hypervascular, spatially clustered
architecture, which is what glioblastoma is. The part of them that *is* a
methodological question is the same one as in §5 — the 20 µm linking distance
simultaneously over-splits (84.8 % singletons on s1) and over-merges (a 1,142-cell
object on s3) — and it is unaffected by this fix either way.

---

## 7. What changed in the code

* **New** `tmesim/reconstruct/vessel_geometry.py` — measured constants with
  provenance, `vessel_geometry` (per-object radius/axis/footprint),
  `subgrid_coverage_field` (fractional-area painter), `vascular_area_fraction`,
  `vascular_heterogeneity` (fidelity-to-sample scoring). Runnable self-test:
  `./.venv/bin/python -m tmesim.reconstruct.vessel_geometry`.
* `vessel_network.segmentation_report` — headline `vessel_area_fraction_pct` is now
  the sub-grid data-derived value; the historical number is kept alongside as
  `legacy_binary_area_fraction_pct` so the size of the correction stays visible.
  Published rows gained `reference_kind` marking them as population context.
* `vessel_network.discrete_vascular_field` — new **opt-in** `geometry='subgrid'`.
  The default remains `'binary'` and is **bit-identical** to the historical field,
  so no already-published result moves silently. `vascular_density_from_cells`
  (and therefore every reconstruction) still routes through the binary default;
  promoting the sub-grid field to the production source is a deliberate re-run and
  is **not** done here.
* `VESSEL_RADIUS_UM = 7.5` and its docstring are retained but marked superseded,
  with the corrected diagnosis recorded in place of the old "the mesh cannot
  represent a smaller vessel" justification.
* **Tests** (`tests/test_vessel_network.py`, suite 162 passed): grid-independence
  of the area fraction across the dx sweep; sub-grid representation of vessels
  smaller than a cell; painter accuracy against a closed form on non-aligned
  discs; the radius being derived rather than assumed; one footprint per vessel
  object rather than per nucleus; the sub-grid field being opt-in with the default
  unchanged; MVD / inter-vessel distance being radius-free; heterogeneity
  preserved where the KDE loses it, and grid-independent.
