# The vascular area fraction and MVD failures: two defects, one of them a modality conflation

*Written after tracing the full chain from measured endothelial cell positions to the two
numbers `results/heldout_patient_validation.md` reports as FAIL, and re-measuring every
step on all three cached real GBM sections. Every number below was observed in this work.
No constant was tuned. No radius and no bandwidth was changed.*

> ### Headline
>
> **The 5.23× is not one defect, it is two, and the larger one is a modality conflation
> of exactly the kind that manufactured two fake regressions earlier in this project.**
>
> **The area-fraction target is the wrong quantity.** The model computes the geometric
> areal fraction of vessel *profiles* — lumen plus wall. The measured 2.1 ± 0.3 % it is
> scored against is, per the source's own Methods, *"the number of positively stained
> pixels on the total pixels in the image"* for a **CD34 immunostain** — endothelium only,
> lumen excluded — over **10 operator-selected 40× fields in the tumour centre** of 4 µm
> sections. That is not the same physical quantity, and the ± 0.3 is a **standard error**,
> not a standard deviation, so the 1.8–2.4 % pass band is a ±1 SE band on a cohort mean
> rather than a cohort range.
>
> **A commensurable anchor already exists, in a paper this module already cites.** Marien
> KM *et al.* (PMID 27583442) — the source of the 65 /mm² MVD number — also reports, for
> the same n = 20 GBM cohort, a **microvessel areal fraction of 4.70 % ± 2.38 %**, measured
> by stereological point counting in which a grid point counts as vascular when it falls
> on *"an endothelial cell **or a vessel lumen**"*. That is precisely what the model
> paints. Rescored against it, the held-out section's miss falls from **5.23× to 2.34×**
> and the development section's from 0.59× to 0.27×.
>
> **The model still fails, and the opposite-sign pattern survives the correction.** So the
> conflation was not the whole story: there is a second, real defect underneath it. It is
> the radius model. Because `r = r₁·k` and a transverse footprint has area `π r²`, the
> painted area scales as **k²** — quadratically in the per-object endothelial count — with
> a hard floor at `k = 1`. Section 1 sits on the floor (mean k = 1.36) and reads low;
> sections 2 and 3 sit well above it (2.32, 2.07) and read high. **One upstream cause, two
> transfer functions, opposite signs.**
>
> **And the mechanism that puts them on opposite sides of that floor is a single measured
> number.** The fixed 20 µm linking distance falls on **opposite sides of the two sections'
> median endothelial nearest-neighbour distance** — **22.9 µm** on the development section
> (so the rule mostly refuses to merge: 81.3 % singletons) and **11.4 µm** on the held-out
> section (so it merges avidly: 62.9 %). One rule, two qualitatively different behaviours.
>
> **MVD is a third thing, and it is over-counted by construction.** The measurement's
> connectivity rule is *continuity of stain*; the model's is *proximity of nuclei at
> 20 µm*. This package's own measured constant says consecutive endothelial nuclei on one
> vessel are **50 µm** apart (`EC_AXIAL_PITCH_UM`, PMID 38723085). Linking at 20 µm is
> provably below that, so one stained vessel is split into one object per nucleus, in the
> same direction on every section — the observed signature.
>
> **Tested directly, the merged-object convention is confirmed in kind but does not close
> the gap on its own.** The merge factor it would need is not a constant: **1.70× / 3.59× /
> 4.93×** across the three sections, tracking each section's endothelial areal density
> (153.9 / 645.2 / 756.7 EC per mm², a 4.9× spread) rather than any universal
> nuclei-per-cross-section figure (§7.3).
>
> **But the linking distance cannot simply be raised.** At 50 µm the held-out section
> **percolates**: a single object swallows **24.1 %** of all endothelial cells. Its
> apparently-in-band MVD of 37.3 /mm² at that distance is an artefact of that giant
> component, not a fix. **I therefore changed no constant.**
>
> **Blast radius: measured, and it is zero.** The production delivery source field is
> **bit-identical** across linking distances 10/20/50/80 µm, because
> `discrete_vascular_field(geometry='binary', amplitude='unit')` uses endothelial
> positions and never the cluster labels. Nothing in this audit touches any delivery, PK
> or kill number. That is reassuring for merge safety and damning for something else: the
> field that actually delivers drug still uses the per-endothelial-cell 7.5 µm painting
> that `vessel_geometry.py`'s own docstring calls **"Wrong #2"**.

---

## 1. The computation chain, end to end

Both reported quantities start from the same input and diverge at step 3.

```
CellPopulation (measured Xenium centroids, marker-gated CellType.ENDOTHELIAL)
  │
  ├─1─ vessel_network._endothelial_xy(cells)            -> xy (n_endo, 2), µm
  │
  ├─2─ vessel_network._single_linkage_labels(xy, eps=LINK_DISTANCE_UM=20 µm)
  │      cKDTree.sparse_distance_matrix -> connected_components
  │      == DBSCAN(eps=20, min_samples=1); no point is noise
  │      -> labels, n_vessels
  │
  ├─3a─ MVD  =  n_vessels / tissue_area_mm2
  │      tissue_area = domain.tissue_mask.sum() * dx² / 1e6
  │      -> segment_vessels()["mvd_per_mm2"]
  │      -> segmentation_report() row "mvd_per_mm2"     vs PUBLISHED 65 (35–95)
  │      (rows "intervessel_spacing_um" and "median_tumor_to_vessel_um" are
  │       DERIVED FROM THIS SAME NUMBER — 1/sqrt(MVD) — so they are not
  │       independent evidence about it)
  │
  └─3b─ per-object geometry, in vessel_geometry.vessel_geometry(radius_mode='data'):
         extent E  = max pairwise nuclear separation within the object   (from step 2)
         N         = endothelial nuclei in the object
         axial_units = max(1, E / EC_AXIAL_PITCH_UM)         EC_AXIAL_PITCH_UM = 50 µm
         k         = max(1, N / axial_units)                 (NO cap by default)
         r         = SINGLE_EC_VESSEL_RADIUS_UM · k          r₁ = 3.113 µm
         axis_len  = max(0, E − 2r)
         footprint = stadium: π r² + 2 r · axis_len
           │
           ├─ vessel_geometry._stadium_disc_centres -> union of discs
           ├─ vessel_geometry._subgrid_coverage     -> fractional coverage on a
           │    FIXED PHYSICAL sub-lattice (pitch = clip(r_min/4, 0.25, 1.0) µm),
           │    binned to the dx grid, masked to tissue
           └─ vessel_geometry.vascular_area_fraction -> area_fraction_pct
                = 100 · Σ(coverage · cell_area)|tissue / Σ(cell_area)|tissue
                -> segmentation_report() row "vascular_area_fraction_pct"
                                                        vs PUBLISHED 2.1 (1.8–2.4)
```

Two properties of this chain matter and are easy to miss:

* **The radius is a function of the nuclear COUNT, and the nuclear EXTENT is used only to
  decide how much of the count is axial.** The object's own extent `E` is an independent
  estimate of its radius — `vessel_geometry`'s docstring says so explicitly ("a transverse
  cut has `E ≈ 2r`") — and the code computes `E`, then discards it as a radius estimate.
* **`r₁ = 3.113 µm` is a hard floor.** Every object with `k = 1` gets exactly that radius,
  so a section whose objects are mostly singletons has its area fraction pinned at
  `n_endo · π r₁²  =  n_endo · 30.4 µm²` and cannot go lower whatever the tissue is doing.

### 1.1 What the production *field* uses — which is not this

`reconstruct.real_fields.vascular_density_from_cells` defaults to `mode='discrete'`, which
calls `segment_vessels` and then `discrete_vascular_field` **with its own defaults**:
`geometry='binary'`, `amplitude='unit'`, `radius_um=VESSEL_RADIUS_UM=7.5`. That path paints
a binary 7.5 µm disc **around every endothelial cell** and — measured, §6 — never reads the
cluster labels at all. So the sub-grid, data-derived geometry described above supplies the
**reported** area fraction only. It has never been promoted into the field that delivers
drug.

---

## 2. Commensurability, checked first — and it fails on both rows

### 2.1 The area fraction: the model and the measurement are not the same quantity

Retrieved from PubMed. The area-fraction anchor is Tamma R, Ingravallo G, Annese T,
d'Amati A, Lorusso L, Ribatti D, *"Tumor Microenvironment and Microvascular Density in
Human Glioblastoma"*, **Cells 2022;12(1):11**, [PMID 36611806](https://pubmed.ncbi.nlm.nih.gov/36611806/),
[DOI](https://doi.org/10.3390/cells12010011). Its Methods state, verbatim:

> "Serial histological sections of **4 um thickness** … Digital slides were inspected with
> Aperio ImageScope … at 20× magnification, and **10 fields with an equal area were
> selected for the analysis at 40× magnification in both the center of the tumor and in
> the healthy surrounding areas**. The protein expression was assessed with the **Positive
> Pixel Count** algorithm … and reported as **positivity percentage, defined as the number
> of positively stained pixels on the total pixels in the image**."

and its Results:

> "Morphometric analysis showed the significant increase of CD34+ into the tumor tissue
> (**CD34: 2.1 ± 0.3%**) as compared to the healthy tissue, CTRL (CD34: 0.58 ± 0.09%)."

Three separate mismatches, all confirmed from the source rather than inferred:

| | the model | the measurement |
|---|---|---|
| **what is counted** | geometric area of the vessel **profile** — lumen **plus** wall | **pixels positive for CD34 immunostain** — endothelium only; the lumen is unstained and contributes nothing |
| **sampling** | every vessel object in the whole tissue mask | **10 operator-selected 40× fields in the tumour centre** |
| **dispersion** | — | the paper's Statistical Analysis says data "are reported as **means ± SE**" and "scatter dot plots with mean and **standard error**". So ±0.3 is an **SE at n = 30**, implying a cohort **SD of 0.3·√30 = 1.64** |

The first row is the load-bearing one and it is the same error class as the
free-interstitial-vs-total-tissue pooling that produced two fake regressions in this
project's history: a stain-positivity fraction is **strictly less than** the vessel-profile
area fraction, by the lumen fraction, which grows with caliber and is not quantified here.

The third row means the declared pass band is wrong even on its own terms: **1.8–2.4 % is a
±1 SE interval on a cohort mean, not a range of tumours.** A ±1 SD cohort range would be
roughly 0.5–3.7 %.

> **I did not widen the band.** Doing so would convert the development section's 1.25 %
> from a miss into a pass, and a band that widens until the model fits inside it is exactly
> the rigging this project forbids — the statistical justification does not change how it
> would read. The defect is documented in `PUBLISHED` and the band is left bit-identical;
> re-scoping it is a reviewer's decision, not an agent's.

### 2.2 The commensurable anchor was already in the repository, uncited for this purpose

Marien KM *et al.*, *"Development and Validation of a Histological Method to Measure
Microvessel Density in Whole-Slide Images of Cancer Tissue"*, **PLoS One 2016;11(9):e0161496**,
[PMID 27583442](https://pubmed.ncbi.nlm.nih.gov/27583442/),
[DOI](https://doi.org/10.1371/journal.pone.0161496) — already cited by this module for the
65 /mm² MVD figure — reports **both** parameters for the same cohort. From its Results:

> "The mean microvessel density (Q_A) and its standard deviation for the four cancer types
> was: 112 ± 50 (RCC), 76 ± 21 (CRC), **65 ± 30 (GBM)**, and 43 ± 13 (OC) vessels per mm².
> **The mean areal fraction of the vessels (A_A) and its standard deviation for the four
> cancer types was: 11.92% ± 6.74% (RCC), 3.87% ± 1.18% (CRC), 4.70% ± 2.38% (GBM), and
> 2.99% ± 1.84% (OC).**"

and its counting rule for A_A, from Methods:

> "A grid point … was regarded as overlapping a microvessel when it fell on **an
> endothelial cell or a vessel lumen**."

**That is the model's quantity.** Lumen included, per unit tissue area, same specimens,
same unbiased whole-slide sampling design, same denominator convention as the MVD anchor.
It is a cohort mean ± 1 **SD** (n = 20), so 2.32–7.08 % is a genuine cohort range.

This is registered as a new `PUBLISHED` entry, `vessel_areal_fraction_pct`, and
`segmentation_report` now emits it as an additional row. The historical
`vascular_area_fraction_pct` row is left **numerically untouched** so nothing already
published moves.

### 2.3 Rescoring against the commensurable anchor

| section | model | vs CD34 stain-positivity 2.1 % (**not commensurable**) | vs vessel areal fraction 4.70 % (**commensurable**) |
|---|---:|---|---|
| 1 (development) | 1.246 % | 0.59× — LOW | **0.27× — LOW** (band 2.32–7.08) |
| 2 (held out) | 10.982 % | **5.23× — HIGH** | **2.34× — HIGH** |
| 3 (gse324268) | 10.938 % | 5.21× — HIGH | **2.33× — HIGH** |

**The headline 5.23× shrinks to 2.34× purely by comparing like with like.** Slightly over
half of the reported discrepancy on the held-out section was never a model error at all;
it was the lumen, the field selection and an SE read as an SD.

**And the model still fails.** Both sign errors survive, and they still run in opposite
directions. So §2 removes about half the magnitude and none of the structure. The
structure is §3.

### 2.4 MVD: the model cannot implement the measured convention

Marien's counting rule, verbatim from its Methods (the abstract states it more briefly):

> "every CD31-positive object, no matter how small, should be counted, **except suspected
> CD31-positive monocytes, macrophages and tumor cells**. Furthermore, **if CD31-positive
> objects were connected, they were considered a single object, while absence of staining
> defined two or more separate objects**."

The measurement's connectivity is **continuity of stain along the vessel wall**. CD31
stains the endothelial cell body, so two nuclei on one vessel are joined by continuous
staining and count once. The model's connectivity is **proximity of nuclear centroids at
`LINK_DISTANCE_UM` = 20 µm**.

These cannot agree, and the *direction* of the disagreement is fixed by a constant this
package already carries and already cites:

> `vessel_geometry.EC_AXIAL_PITCH_UM = 1000 / EC_PER_MM_VESSEL = 50 µm` — the spacing of
> consecutive endothelial nuclei along one vessel, from ~4,604 endothelial cells lining
> ~230 mm of vessel in the H01 human cortex reconstruction (PMID 38723085).

**20 µm < 50 µm.** The linking distance is provably below the nuclear pitch it would have
to span, so consecutive nuclei on a single stained vessel are split into separate objects
and MVD is over-counted **by construction, in the same direction on every section** —
precisely the signature `heldout_patient_validation.md` identified but attributed to
"over-splitting" without naming the constant that proves it.

The two constants are mutually contradictory within the same package: `vessel_network` links
at 20 µm while `vessel_geometry` simultaneously asserts that same-vessel nuclei are 50 µm
apart, and `vessel_geometry` uses the 50 µm figure to *divide* the nuclear count.

Marien names two further biases that both push the same way and that the model does not
control: CD31 "is not only expressed by endothelial cells, but also by other cells, such as
macrophages and platelets, which may result in an **overestimation of the number of
microvessels**" (the model's marker gate performs no such exclusion), and "depending on the
orientation of the vessel and the direction of the sectioning, **one vessel can appear as
separate shapes in the two-dimensional section**."

---

## 3. The second defect: the painted area is quadratic in the per-object nuclear count

This is what remains after §2, and it is the reason the two sections miss in opposite
directions.

For a **transverse** object (`E ≤ 2r`, which is the majority) `axis_len = 0` and the
footprint is a bare disc:

```
area  =  π r²  =  π (r₁ k)²  =  π r₁² · k²
```

The area is **quadratic in `k`**, the object's endothelial nuclei per cross-section, while
essentially every upstream error — how sensitive the marker gate is, how many nuclei the
linking rule pulls into one object — is **linear** in it. And `k` is uncapped by default
(`max_ec_per_cross_section=None`), so nothing bounds the amplification.

The consequence is a floor and a runaway on the same curve:

| | section 1 (dev) | section 2 (held out) | section 3 |
|---|---:|---:|---:|
| endothelial cells /mm² | **153.9** | **645.2** | **756.7** |
| median nearest-EC distance | **22.9 µm** (> 20 µm link) | **11.4 µm** (< link) | **12.9 µm** (< link) |
| singleton objects | **81.3 %** | 62.9 % | 58.5 % |
| MVD /mm² @ 20 µm | 110.8 | 233.3 | 320.4 |
| mean `k` @ 20 µm | **1.36** | **2.32** | **2.07** |
| area fraction | 1.246 % | 10.982 % | 10.938 % |
| vs commensurable 4.70 % | 0.27× LOW | 2.34× HIGH | 2.33× HIGH |

Section 1's objects are overwhelmingly singletons — because its median endothelial
nearest-neighbour distance (22.9 µm) is **larger** than the 20 µm linking distance — so its
area fraction is pinned near the `k = 1` floor of `n_endo · π r₁²` and **cannot** reach the
measured band no matter what the tissue looks like. Sections 2 and 3 have nearest-neighbour
distances (11.4, 12.9 µm) well **inside** the linking distance, so they sit above the floor
where the `k²` term takes over. §7.1 develops this.

**That is the answer to "why do the two failures have different shapes".** They are not two
independent problems. They are one upstream defect — the endothelial point pattern → vessel
object mapping of §2.4 — read out through two different transfer functions:

* **MVD is (approximately) linear** in that defect, so both sections miss the same way, and
  its miss is the honest size of the over-splitting: 1.70× / 3.59× / 4.93×.
* **The area fraction is quadratic with a floor**, so the same defect drives one section
  below the band and the others above it.

A single cause genuinely cannot produce both patterns — but a single cause *plus a
nonlinear read-out* can, and does.

### 3.1 The internal inconsistency, with no published number involved

`vessel_geometry`'s own model states that a transverse cut has `E ≈ 2r`. That makes the
object's measured extent an **independent** estimate of its radius, which the code computes
and then discards. Where `r₁·k` exceeds the extent-implied radius, the painted disc is wider
than the nuclei that generated it — geometrically impossible for one vessel profile.

`experiments/vascular_segmentation_audit.py` measures this per object, using the exact chord
relation for `k` nuclei tiling a circle (largest pairwise chord `2r·sin(π⌊k/2⌋/k)`, exact at
`k = 2`). Measured: **40.8 % of multi-cell objects on the held-out section** (35.0 % on
section 3, 13.8 % on section 1) are painted **wider than the nuclei that generated them**,
and the largest assigned radii — 88.5 µm and 129.2 µm — are 4.7× and 6.9× the largest brain
vessel *diameter* this same module cites. See §7.2.

This is a defect diagnosable **without reference to any measured target**, which is why it
is reported separately from §2.

---

## 4. Why the linking distance was NOT changed

§2.4 gives a principled, non-fitted candidate correction: link at the package's own measured
endothelial axial pitch, `EC_AXIAL_PITCH_UM = 50 µm`, so that consecutive nuclei on one
vessel merge as the measurement's convention requires. 50 µm was fixed by a cited human EM
measurement long before anyone looked at MVD; it is not a tuned parameter.

Measured, on all sections (`experiments/vascular_segmentation_audit.py`):

**Section 1 (development)**

| link (µm) | MVD /mm² | giant component (% of all EC) | mean `k` | area fraction |
|---:|---:|---:|---:|---:|
| 10 | 148.7 | 0.1 % | 1.04 | 0.510 % |
| 15 | 132.5 | 0.3 % | 1.16 | 0.746 % |
| **20 (production)** | **110.8** | **0.5 %** | **1.36** | **1.245 %** |
| 25 | 95.5 | 1.2 % | 1.51 | 1.692 % |
| 30 | 86.8 | 1.8 % | 1.61 | 2.022 % |
| 40 | 73.8 | 1.9 % | 1.75 | 2.554 % |
| **50 (= EC axial pitch)** | **63.3** | **1.9 %** | 1.85 | 2.997 % |

**Section 2 (held out)**

| link (µm) | MVD /mm² | giant component (% of all EC) | mean `k` | area fraction |
|---:|---:|---:|---:|---:|
| 10 | 482.2 | 0.0 % | 1.34 | 3.847 % |
| 15 | 316.6 | 0.2 % | 1.94 | 8.170 % |
| **20 (production)** | **233.3** | **0.3 %** | **2.32** | **10.982 %** |
| 25 | 179.5 | 0.4 % | 2.57 | 13.000 % |
| 30 | 136.4 | 0.5 % | 2.75 | 14.628 % |
| 40 | 73.9 | **6.4 %** | 2.97 | 16.714 % |
| **50 (= EC axial pitch)** | **37.3** | **24.1 %** | 3.02 | 17.207 % |

At 50 µm **both** sections' MVD land inside the published 35–95 band (63.3 and 37.3). That
looks like a fix, and it is not one:

1. **Section 2 percolates.** At 50 µm a *single* connected object contains **24.1 % of every
   endothelial cell in the section**. The MVD of 37.3 /mm² is produced by that giant
   component absorbing a quarter of the endothelium, not by correctly merging vessels. The
   number enters the band on its way through, not because it became right. Section 1, which
   is 4.5× less endothelial-dense, does not percolate (giant component 1.9 % at 50 µm) — so
   the correction is *section-dependent in kind*, not merely in degree.
2. **It makes the area fraction dramatically worse**, exactly as §3 predicts: merging raises
   `k`, and the area is quadratic in `k`. Section 2 goes 10.98 % → 17.21 %, i.e. from 2.34×
   to 3.66× the commensurable 4.70 %.
3. **Single linkage is the wrong operator for the job.** The measured convention is
   connectivity of a *stained tube*; single linkage on centroids has no notion of a tube and
   percolates as soon as the point density exceeds a threshold set by the linking distance.

So: the MVD defect is **real, isolated and explained**, and it is **not closeable by moving
the linking distance**. Implementing the measured convention needs an operator that can
represent an elongated stained object — an anisotropic or graph-based linking along a fitted
local vessel axis, or a caliber/perimeter signal — not a larger isotropic ball.

This also independently confirms the conclusion of
`results/sweep/vascular_area_fraction_generalization_2026-07-30T03-18-42Z.md`: the required
radius correlates inversely with MVD across sections because both are driven by the same
endothelial-density defect, and conditioning one on the other launders it.

### 4.1 An operational landmine found on the way

`segment_vessels` computes each object's extent with a dense pairwise distance matrix
(`np.linalg.norm(pts[:, None, :] - pts[None, :, :], axis=-1)`), which is **O(n²) in memory
for the largest single cluster**. At the production 20 µm this is harmless. Once a section
percolates it is not: the first run of this audit had to be killed at `link = 60 µm` on
section 2, where the giant component would have required a ~40 GB allocation. Anyone raising
`LINK_DISTANCE_UM`, or running these sections' denser neighbours, will hit it.

---

## 5. What was changed

**One file under `tmesim/`: `tmesim/reconstruct/vessel_network.py`.** The change is additive
and **numerically inert for every pre-existing row**.

1. **`PUBLISHED["vascular_area_fraction_pct"]`** — `what` corrected to describe what the
   source actually measured (CD34 **stain-positivity** percentage, Aperio Positive Pixel
   Count, **10 selected 40× tumour-centre fields**, 4 µm sections) instead of the previous
   "whole-field positive-pixel count"; first author corrected to Tamma R *et al.*; new
   `measures`, `commensurable_with_model=False`, `dispersion_is` and `caveat` fields
   recording the lumen-exclusion conflation and the SE-read-as-SD band defect. **`value`,
   `sd` and `range` are untouched**, so no verdict moves.
2. **`PUBLISHED["vessel_areal_fraction_pct"]`** — NEW. The commensurable anchor: 4.70 %
   ± 2.38 % (range 2.32–7.08), GBM n = 20, CD31 whole-slide SURS stereological point
   counting with the lumen included, PMID 27583442. Marked
   `commensurable_with_model=True`.
3. **`PUBLISHED["mvd_per_mm2"]`** — new `measures` field quoting the source's counting rule
   verbatim, and a `caveat` recording the convention mismatch of §2.4 including the
   20 µm-vs-50 µm contradiction. Value, sd and range untouched.
4. **`segmentation_report`** — every row now carries `measures`, `commensurable` and
   `caveat`, and one additional row `vessel_areal_fraction_pct` scores the **same**
   reconstructed area fraction against the commensurable anchor, so the size of the modality
   correction is always visible next to the historical comparison.

Nothing else under `tmesim/` was modified. In particular
`tmesim/validation/heldout_patient.py`, `tmesim/validation/endothelial_density.py`,
`tmesim/validation/kill_calibration.py`, `tmesim/validation/scorecard.py` and every existing
file under `results/` are untouched.

**What was deliberately NOT changed, and why**

* **`LINK_DISTANCE_UM` stays at 20 µm** — §4. The principled value percolates on the denser
  sections.
* **`SINGLE_EC_VESSEL_RADIUS_UM`, `VESSEL_RADIUS_UM`, the `r ∝ k` law, the sub-grid pitch,
  every bandwidth** — untouched. The generalisation sweep already showed no `k`-based radius
  model closes both sections, and §3 explains why: the defect is the *quadratic read-out of
  an upstream error*, not the value of a constant.
* **The 1.8–2.4 % band** — §2.1. Widening it to the statistically correct cohort range would
  flip section 1 to PASS. That must be a reviewed decision.

---

## 6. Blast radius — measured, not asserted

The brief flags that `PS` and `vascular_density` are perfectly confounded
(`results/transvascular_exchange_revalidation.md` §5), so any vessel-map change propagates
into every delivery number. **For the changes made here the propagation is exactly zero, and
for a linking-distance change it would also be zero.** Measured directly:

```
link= 10.0  n_vessels=  4149  field identical to production: True  max|diff|=0.000e+00
link= 20.0  n_vessels=  3091  field identical to production: True  max|diff|=0.000e+00
link= 50.0  n_vessels=  1765  field identical to production: True  max|diff|=0.000e+00
link= 80.0  n_vessels=   825  field identical to production: True  max|diff|=0.000e+00
```

The reason is structural: with the production defaults `geometry='binary'` and
`amplitude='unit'`, `discrete_vascular_field` builds a `cKDTree` over `vessels["xy"]` — *all*
endothelial cells — and asks only "is any endothelial cell within 7.5 µm of this grid-cell
centre". **The cluster labels are never read.** They are read only under
`amplitude='per_vessel_cells'`, which is not the default anywhere on the production path.
The production oxygen route is likewise unaffected: `OXYGEN_MODES` defaults to `'pde'`, and
the vessel-distance proxy defaults to `source='endothelial'`, not `'vessel_network'`.

So the vessel **segmentation** is, today, a pure reporting layer: it produces MVD and the
area fraction and feeds no physics.

**The uncomfortable corollary.** The delivery source is *not* the corrected sub-grid
geometry. It is still the binary 7.5 µm disc painted around **every endothelial cell** —
the construction `vessel_geometry.py`'s docstring calls **"Wrong #2 — the footprint is
painted per endothelial CELL, not per vessel … This is the defect that actually moves the
number."** That fix reached the *reported* area fraction and never reached the *field*. Any
future decision to promote `geometry='subgrid'` into production would be a genuine,
shape-changing perturbation of every delivery number — and not one the `PS` gauge can
absorb, because it changes *where* the source sits, not merely its scale. §7 quantifies how
different the two fields are on each section; that number is the real blast radius of the
fix this project has already written but never shipped.

---

## 7. Measured diagnostics

`experiments/vascular_segmentation_audit.py`, dx = 12 µm, production 20 µm linking, all
three cached real GBM sections.

| | **s1 (dev)** | **s2 (held out)** | **s3 (gse324268)** |
|---|---:|---:|---:|
| endothelial cells /mm² | **153.9** | **645.2** | **756.7** |
| MVD /mm² | 110.8 | 233.3 | 320.4 |
| EC per vessel object | 1.39 | 2.77 | 2.36 |
| singleton objects | **81.3 %** | 62.9 % | 58.5 % |
| **median nearest-EC distance (µm)** | **22.9** | **11.4** | **12.9** |
| EC with a neighbour ≤ 20 µm | **41.5 %** | 77.3 % | 75.2 % |
| mean `k` | 1.36 | 2.32 | 2.07 |
| p99 object radius (µm) | 18.7 | **36.6** | **30.2** |
| max object radius (µm) | 39.9 | **88.5** | **129.2** |
| area fraction | 1.245 % | 10.982 % | 10.938 % |
| painted area per EC (µm²) | 80.9 | 170.2 | 144.5 |
| multi-cell objects painted **wider than their own nuclei** | 13.8 % | **40.8 %** | 35.0 % |
| median `r_model / r_extent-implied` (multi-cell) | 0.76 | 0.92 | 0.88 |

### 7.1 The single most explanatory number in this table

> **The fixed 20 µm linking distance falls on OPPOSITE SIDES of the two sections' median
> endothelial nearest-neighbour distance: 22.9 µm on the development section and 11.4 µm on
> the held-out section.**

On section 1 the linking distance is *below* the typical nuclear spacing, so the rule mostly
**refuses to merge**: 81.3 % of objects are singletons and only 41.5 % of endothelial cells
have any neighbour within reach. `k` is pinned near 1, so the radius is pinned at the
`r₁ = 3.113 µm` floor and the area fraction cannot rise — **LOW**.

On sections 2 and 3 the linking distance is *above* the typical spacing, so the rule
**merges avidly**: mean `k` reaches 2.32 / 2.07 and the `k²` term of §3 takes over — **HIGH**.

One rule, two qualitatively different behaviours, decided by which side of 20 µm the
section's endothelial density puts it on. That is the mechanism behind the opposite signs.

### 7.2 Radii that no published brain vessel matches

`PUBLISHED["vessel_diameter_um"]` in this same module cites a maximum in-vivo human brain
caliber of **37.4 µm diameter** (cortical grey matter, PMID 19856441). The model assigns
single objects radii of up to **88.5 µm (s2)** and **129.2 µm (s3)** — diameters of 177 and
258 µm, i.e. **4.7× and 6.9× the largest caliber the module itself cites**. At the 99th
percentile the radii are 36.6 / 30.2 µm, i.e. diameters roughly **2× the largest cited
caliber**. Nothing bounds this because `max_ec_per_cross_section` defaults to `None`.

This is a defect visible **without any external areal target**: it is the module's own
vessel-diameter citation contradicting its own radius model.

**And the `k²` scaling has quietly undone the fix it was introduced to make.** The legacy
per-endothelial-cell painting deposits `π · 7.5² = 176.7 µm²` per endothelial cell. The
"fixed" sub-grid, per-*object* geometry deposits **170.2 µm² per endothelial cell on the
held-out section** — **96 % of the legacy per-cell figure**, and 144.5 µm² (82 %) on section
3. On section 1 it deposits 80.9 µm² (46 %). So on exactly the sections where
`vessel_geometry.py`'s "Wrong #2" argument said the per-cell painting inflates the area, the
replacement reproduces it almost exactly — because `area ∝ k²` while the number of objects
falls as `1/k`, so the two nearly cancel and the total goes back to scaling with the
**endothelial-cell count**. The fix removed the per-cell painting in form and reinstated its
scaling in substance.

### 7.3 The merged-object convention, tested directly

The hypothesis under test: Marien's rule merges all interconnected stained objects into one
vessel, so the model over-counts by the number of endothelial nuclei per counted vessel.

| | s1 | s2 | s3 |
|---|---:|---:|---:|
| model MVD / published 65 | **1.70×** | **3.59×** | **4.93×** |
| EC already merged per object at 20 µm | 1.39 | 2.77 | 2.36 |
| ⇒ EC per vessel required by the merged convention | **2.36** | **9.94** | **11.6** |

**Verdict on the lead: confirmed in kind and in direction, ruled out as a complete
quantitative explanation.** The convention mismatch is real (§2.4) and it correctly predicts
that every section misses the same way. But the merge factor it would need is **not a
constant**: it varies **2.9× across three sections** (1.70 / 3.59 / 4.93) and tracks each
section's endothelial areal density (153.9 / 645.2 / 756.7 EC/mm², a 4.9× spread) rather
than any universal nuclei-per-cross-section figure. A fixed "3–4 nuclei per cross-section"
correction would fit section 2 and miss sections 1 and 3 in opposite directions.

So MVD carries **two** superposed errors: the counting-convention mismatch (fixed direction,
affects all sections) and a section-to-section endothelial-density inconsistency that the
convention cannot absorb. The second is upstream of this module, in cell typing — Marien
excludes CD31+ monocytes, macrophages and tumour cells; the marker gate does not.

### 7.4 A framing correction worth recording

Two points on how these anchors are usually described in this project:

* **The MVD anchor is not Weidner hot-spot counting.** Marien's method is **systematic
  uniform random sampling**, and the paper explicitly contrasts the two: *"the hotspot
  method (Weidner's method) will quantify the strongest angiogenic areas of a tumor and the
  SURS method provides global information."* Hot-spot counting yields **higher** numbers, so
  scoring against a hot-spot convention would make the model's over-count look *smaller*
  than it is. The 65 /mm² is the global, unbiased figure, which is the harder and the
  correct bar for a whole-section reconstruction.
* **The two historical anchors are different assays on different cohorts**: CD34 (n=30,
  Tamma) for the area fraction, CD31 (n=20, Marien) for MVD. The commensurable anchor added
  in §2.2 is **CD31, n=20, same paper and same cohort as the MVD anchor**, so the two scored
  morphometry rows now come from one assay on one cohort instead of two — a consistency gain
  independent of the commensurability argument.

---

## 8. What could make this wrong

1. **The commensurable anchor is still a cohort mean over 20 patients**, and Marien's own
   paper reports A_A as the *less* reproducible of its two parameters — inter-ROI ICC 0.6 in
   GBM against 1.0 for the profile count, because point counting samples only a small area.
   A 2.34× miss survives that comfortably; a 1.3× miss would not.
2. **The lumen-exclusion argument is directional, not quantitative.** I have shown that CD34
   stain-positivity must under-report vessel-profile area, and therefore that 2.1 % is a
   lower bound on the model's quantity. I have **not** quantified by how much, and I did not
   try — doing so would require assuming a wall thickness and a DAB threshold, which is
   exactly the kind of free parameter this audit exists to avoid.
3. **n = 3 sections, all GBM, all Xenium, two panels.** The floor-vs-runaway account of §3
   rests on three points, and section 1 differs from sections 2 and 3 in gene panel as well
   as in endothelial density.
4. **The `k²` account assumes the transverse regime.** It is exact where `axis_len = 0`; for
   elongated objects the stadium term `2·r·axis_len` is linear in `r` and the scaling is
   intermediate. The transverse regime dominates at the production linking distance but not
   at 50 µm.
5. **`EC_AXIAL_PITCH_UM = 50 µm` is normal human cortex** (PMID 38723085), not GBM. GBM
   endothelium is hyperplastic and its nuclear pitch is plausibly shorter, which would
   narrow — but not close — the 20-vs-50 µm gap of §2.4. That the observed median nearest-
   endothelial-neighbour distance on the development section (26.4 µm, per
   `results/vascular_validation.md` §3.2) also exceeds the 20 µm linking distance is
   independent support that the gap is real.
6. **Nothing here validates the model's area fraction as correct.** It relocates about half
   the reported error to a measurement mismatch and explains the structure of the rest. Both
   rows still miss, and §4 shows the remaining defect is not closeable by any parameter this
   audit was willing to move.

---

## 9. Reproduction

```bash
./.venv/bin/python -m experiments.vascular_segmentation_audit     # no PDE solves
./.venv/bin/python -m pytest -q                                   # 754 passed
```

Deterministic. No PDE solve is required for anything in this report.

**Runtime warning, and it is the §4.1 landmine.** The per-section production diagnostics of
§7 cost ~40 s each. The *link sweep* is much more expensive on the two dense sections and
degrades sharply above 30 µm as the giant component grows: section 2 took 40 s at 40 µm and
245 s at 50 µm, and the sweep on section 3 above 30 µm did not complete within the budget of
this work and its 40/50 µm rows are therefore **not reported**. Every s3 figure quoted here
is from the production 20 µm segmentation or from the ≤ 30 µm sweep rows, which did
complete. The §4 conclusion rests on sections 1 and 2, where the full sweep is measured.

**Test result: 754 passed, 0 failed** — the declared baseline, unchanged. This is expected
rather than lucky: every pre-existing number emitted by `segmentation_report` is
bit-identical (§5), and §6 measures the production field to be bit-identical too.

**Files changed by this work:** `tmesim/reconstruct/vessel_network.py` (metadata + one
additional scored row; no pre-existing number moves), `experiments/vascular_segmentation_audit.py`
(new), and this file.

*Article metadata and full texts retrieved from PubMed.*
