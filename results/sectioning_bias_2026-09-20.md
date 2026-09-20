# Does taking a 2-D section of a 3-D vascular network bias this engine's vascular statistics?

**Date: 2026-09-20. Route R10 of `results/verification_routes_2026-09-20.md`.**

**§0 is a PRE-REGISTRATION**: it was written and saved before `experiments/sectioning_bias.py`
was run even once, and it fixes every threshold in advance, including the ones that were
then failed. **§1 onward were appended after the run** and report what it produced,
whichever way each threshold fell. Nothing in §0 was edited afterwards.

> ## VERDICT, in one paragraph
>
> **The sectioning bias does NOT explain the 3.59× MVD failure. REFUTED.** Over 283
> random 2-D sections of a real 3-D tumour vascular network, the largest profile count
> any single section produced was **1.99×** the ensemble mean; the fraction reaching
> 3.589× was **0 / 283**. The pre-registered statistic reached 1.769 against a bar of
> 2.871 — **T3: PARTIALLY ACCOUNTS, i.e. it does not account for it.**
>
> **For the 2.34× areal-fraction failure the pre-registered bar is formally MET
> (`R_total` = 3.302 vs 1.869, T3: ACCOUNTS FOR) — and the pass is weak, as §0.6 said in
> advance a pass under this deliberately generous test would be.** A single section
> reaches that miss 3.5 % of the time; the repository has *two* sections at 2.337× and
> 2.327×, a joint probability of 0.125 %, so luck does not explain the reproducibility.
> The same spread also brackets the third section's 0.27× miss in the opposite direction.
>
> **The finding that survives both halves: every vascular statistic tested FAILS the
> single-section spread threshold T2.** The areal fraction varies **8.3×** between the
> 5th and 95th percentile of random sections *of one unchanging tissue*. This engine has
> one section, and carries no such interval anywhere.
>
> **One systematic bias does exceed the T1 threshold:** a 5 µm section over-states vessel
> caliber by **1.63×**, and the engine's painted areal fraction is quadratic in radius.
>
> Delesse's principle holds to 3.9 % as it must, so there is **no** systematic
> areal-fraction bias to be recovered. Nothing was tuned; no shipped number moves.

`calibrate_alpha()` at START of this work = **22.223477711172563** (asserted equal to
22.223477711172563; assertion passed). Re-asserted at the end — see §8.

---

## 0. PRE-REGISTRATION

### 0.1 The hypothesis under test, stated so it can lose

The repository carries two vascular failures on the held-out GBM section
(`results/heldout_patient.md` §2, `results/vascular_segmentation_audit.md` §2.3):

| quantity | model, section 2 (held out) | anchor | ratio |
|---|---:|---|---:|
| `mvd_per_mm2` | **233.27** /mm² | Marien 2016 GBM **65 ± 30** /mm² (PMID 27583442) | **3.589× HIGH** |
| vascular areal fraction | **10.982 %** | Marien 2016 GBM **4.70 ± 2.38 %**, commensurable point count | **2.336× HIGH** |

The scouting hypothesis (R10): these are a **stereology/commensurability defect** — a
2-D section of a 3-D tube network does not report the 3-D truth, and nobody has
computed the size of the discrepancy on a network whose 3-D truth is known.

**I am testing this, not confirming it.** The hypothesis is in serious trouble before
any measurement, for a reason that must be stated up front and that constrains what a
"pass" could even mean:

> **Both sides of both failing comparisons are already 2-D section measurements.** The
> engine's MVD is objects per mm² on one 2-D Xenium section. Marien's Q_A is CD31
> profiles per mm² on one 2-D histology section. The engine's areal fraction is painted
> profile area per unit section area; Marien's A_A is stereological point counting on a
> section. **A bias that applies equally to both cancels identically out of the ratio and
> cannot explain one micron of either failure.** The only part of a sectioning bias that
> can survive into a model-vs-anchor ratio is the part that is *asymmetric* between the
> two protocols. There are exactly three candidate asymmetries, and they are what this
> route can actually test:
>
> * **(a) Section thickness.** Xenium is a nominal **5 µm** FFPE section
>   (`tmesim/validation/dimensionality.py:121`, `SECTION_THICKNESS_UM = 5.0`). Marien
>   2016 **does not state a section thickness anywhere in its text** — I fetched the full
>   open-access article and checked (§1.3). A thicker section over-projects: it shows
>   every vessel that passes anywhere through the slab, so profile counts and areal
>   fractions both rise with thickness.
> * **(b) Section orientation**, if the vasculature is anisotropic. A plane cutting an
>   anisotropic network across the grain reports something different from one cutting
>   along it.
> * **(c) Single-section luck.** The engine's number is **one** section. Marien's is a
>   cohort mean over 20 GBM patients × 15 systematically sampled ROIs. A single section
>   can deviate from a cohort mean by whatever the section-to-section spread allows, and
>   nothing averages that down on the model side.

So the pre-registered question is **not** "is there a sectioning bias" (there is, for
some statistics, by elementary stereology). It is: **can the asymmetric part of it reach
3.59× and 2.34×?**

### 0.2 The statistics I will compute — taken from the code, not assumed

Read from `tmesim/reconstruct/vessel_network.py::segment_vessels` (lines 388–535) and
`segmentation_report` (line 929), which are the functions that produce the failing rows:

| # | engine statistic | code definition | published anchor |
|---|---|---|---|
| S1 | `mvd_per_mm2` | `n_vessels / tissue_area_mm2`, where `n_vessels` is the number of single-linkage clusters of endothelial centroids at `LINK_DISTANCE_UM = 20` | 65 (35–95) /mm², Marien Q_A |
| S2 | `vascular_area_fraction_pct` | `100 · Σ(sub-grid coverage · cell area)|tissue / Σ(cell area)|tissue` from `vessel_geometry.vascular_area_fraction` — vessel **profile** area (lumen + wall) per unit section area | 4.70 ± 2.38 % (Marien A_A, commensurable) and 2.1 ± 0.3 % (Tamma CD34 stain positivity, **not** commensurable) |
| S3 | `intervessel_spacing_um` | `1000 / sqrt(mvd_per_mm2)` — **mechanically derived from S1**, not independent evidence | 124 (102.6–169) µm, itself derived from the same 65 /mm² |
| S4 | `median_tumor_to_vessel_um` | median over tumour centroids of the 2-D distance to the nearest **endothelial centroid** (`cKDTree(xy).query(tumor_xy, k=1)`) | 58 (48.2–79.4) µm, derived from 65 /mm² under a Poisson null |

S3 is excluded from the bias analysis as a statistic in its own right: it is a
deterministic function of S1 on both the model side and the anchor side, so its bias
factor is exactly `S1_bias^(−1/2)` and reporting it separately would double-count. It is
reported only as that derived consequence.

S4 is the one that matters physically — it is the length scale that sets everything in a
diffusion model — and it is the one where a 2-D section has an **unavoidable,
one-directional** bias: the nearest vessel to a point in 3-D may be out of the section
plane, so a 2-D measurement can only ever over-state the distance. On the model side S4
is genuinely measured (not derived from S1), so it is a real test.

I will additionally compute **S5, apparent vessel caliber** (equivalent-circle diameter
of the section profile) against the true 3-D diameter distribution, because the engine's
areal fraction is quadratic in an assumed radius and a radius bias is the second-order
route by which a sectioning effect could enter S2.

### 0.3 The 3-D reference quantities

On the true 3-D network, over a stated reference volume `V_ref`:

* `V_V` — vessel **volume fraction** = (volume of the union of the vessel cylinders) / `V_ref`.
  Computed as a union, not a sum of `π r² L`, so overlaps at branch nodes are not double counted.
* `L_V` — vessel **length density** = `Σ L_i / V_ref`.
* `d3D` — for uniformly sampled tissue points, the exact distance to the nearest vessel
  **surface**, `max(0, dist_to_axis − r)`, minimised over segments.
* the length-weighted 3-D diameter distribution.

### 0.4 How I will section

* The network is 104 straight cylindrical segments (see §1 — the scout said 106; that
  count includes the two header lines).
* A section plane is `{x : x·n = s}`. For the **primary** arm, `n` is drawn **isotropic
  uniform random (IUR)** — a normalised 3-D Gaussian, restricted to a hemisphere by
  symmetry — and `s` is drawn uniformly over the plane's full support across the
  reference box. IUR is the design under which the classical stereological identities are
  exact, so it is the design under which a deviation is attributable to the network
  rather than to me.
* **N = 400 IUR planes.** Separately, **120 planes per coordinate axis normal** (x, y, z)
  for the anisotropy arm (0.1b above), 360 in total.
* Section **thicknesses** `t ∈ {0, 4, 5, 15} µm`. `t = 0` is the mathematical plane; 4 µm
  is Tamma's stated histology thickness; 5 µm is the Xenium nominal; 15 µm is the
  Abercrombie-corrected effective depth `results/vessel_geometry_fix.md` already names as
  the upper end of the repo's own factor-of-three `t_eff` ambiguity. A section of
  thickness `t` is rendered as the **projection** of the vessels through the slab, at 1 µm
  sub-steps, which is what a transmitted-light histological section shows.
* Each section is **rasterised** at 1.0 µm pitch over the plane ∩ reference box, and every
  pixel is tested for membership in every (pruned) cylinder. This makes the areal fraction
  a genuine union — it is exactly Delesse point counting — and lets the profile count use
  Marien's own connectivity rule (*"if CD31-positive objects were connected, they were
  considered a single object"*) by taking **connected components of the in-plane vessel
  mask**, rather than counting axis crossings.
* For the classical identity check I additionally count **axis crossings** (segments whose
  centreline crosses the plane), which is the quantity the identity `Q_A = L_V/2` is about.
* Planes whose in-box tissue area is below **0.05 mm²** are excluded from the per-plane
  distributions (a field smaller than ~224 µm square is not a field a histologist would
  score), but are retained in the ratio-of-sums controls, which must stay unbiased.
* Seed **20260920**, fixed, single run, no re-draws.

### 0.5 IMPLEMENTATION CONTROLS — if these fail, the run is VOID, not the hypothesis

These are code verification, not results. Both are exact consequences of stereology for
an IUR design, independent of what the network looks like:

* **C1 (Delesse).** Ratio-of-sums areal fraction over all 400 IUR planes at `t = 0`,
  divided by the true `V_V`, must lie in **[0.95, 1.05]**.
* **C2 (Q_A = L_V/2).** Ratio-of-sums axis-crossing count per unit area over all 400 IUR
  planes, divided by `L_V / 2`, must lie in **[0.90, 1.10]**.

If either fails, my rasteriser or my sampler is wrong and no finding below is reportable.

### 0.6 PRE-REGISTERED THRESHOLDS

Let `B = median over planes of (2-D estimate / 3-D reference)`.

* **T1 — BIASED-UNSAFE.** A statistic is declared biased-unsafe if `B ∉ [0.667, 1.500]`,
  i.e. the systematic sectioning bias is at least 1.5× in either direction.
* **T2 — SPREAD-UNSAFE.** A statistic is declared spread-unsafe *for use on one section*
  if `P95 / P05` of the per-plane estimate exceeds **2.0**. This is the criterion that
  matters for this engine, which only ever has one section.

The connection test — the one that decides R10:

* **T3 — EXPLANATORY POWER.** For a known failure of size `F` (F = 3.589 for S1, F = 2.336
  for S2, both in the HIGH direction), define the **asymmetric residual** `R`, the largest
  ratio that sectioning alone can manufacture between *one* 2-D section (the model side)
  and the *ensemble mean* of 2-D sections (the anchor side, which is a 20-patient cohort
  mean and therefore not subject to single-section luck):

  ```
  R_single = P97.5( per-plane estimate ) / mean( per-plane estimate )     [same thickness]
  R_thick  = median( estimate at t = 5 µm ) / median( estimate at t = 4 µm )
  R_total  = R_single · R_thick        (the two asymmetries compounded, pro-hypothesis)
  ```

  Verdict, fixed now:
  * **ACCOUNTS FOR F** if `R_total ≥ F / 1.25` in the HIGH direction (≥ 2.871 for S1, ≥ 1.869 for S2).
  * **PARTIALLY ACCOUNTS** if `1.5 ≤ R_total < F / 1.25`.
  * **DOES NOT ACCOUNT** if `R_total < 1.5`, or if the direction is wrong.

  `R_total` is deliberately constructed to be **generous to the hypothesis**: it takes the
  97.5th percentile of single-section luck (not the median), compounds it with the full
  thickness effect, and ignores that the two asymmetries are not independent. A refutation
  under a generous test is a strong refutation; a pass under a generous test is a weak pass
  and will be labelled as one.

* **T4 — the additional, wider claim.** I will also report the **naive** bias `B` of each
  2-D statistic against its 3-D reference (T1/T2), because that is what R10 literally
  asked for and because it is the number that matters for any future claim of the form
  "this section tells you about the tissue". **T4 is reported whether or not T3 passes,
  and a T4 finding is never substituted for a failed T3.**

### 0.7 What a FAIL looks like, written before the run

If `R_total < 1.5` for S1 and S2 — which is what I expect, for the reason in §0.1 — then
**the sectioning hypothesis is REFUTED as an explanation for the 2.34× and the 3.59×**,
the two failures remain unexplained by stereology, and
`results/vascular_segmentation_audit.md`'s existing attribution (single-linkage
over-splitting at a 20 µm threshold below the 50 µm endothelial nuclear pitch, read out
through a radius model quadratic in the per-object nuclear count) stands unchallenged.
That is a real result: it stops a wrong explanation from taking hold, and it is reported
in the headline.

---

*(Sections 1 onward were written after the run.)*

---

## 1. What was fetched, and what it actually contains

### 1.1 The 3-D network — fetched by me, re-fetched by the script on every run

| | |
|---|---|
| page | `https://sites.arizona.edu/secomb/three-dimensional-network-data-rat-tumor-1998/` |
| file | `https://sites.arizona.edu/secomb/files/2024/01/rattum98-ffab1d1a873244ec.txt` |
| fetch | **HTTP 200, 8,049 bytes, 106 lines**, 2026-09-20T09:11:30Z |
| sha256 | `fc2bc2233268f27c2be55d9214cd49ba63d177e340f5d576453eaf5f9ab0c1e0` |

The page states, verbatim:

> "R3230AC mammary carcinoma in rat dorsal skin flap preparation."

> "Segments and nodes (end points of segments) are numbered. Each row in the network data
> file refers to one segment. Segment diameter, length, and Cartesian coordinates of nodes
> are given in microns. **The overall dimensions of the region are 550 x 520 x 230 microns.**"

and gives the citation: T.W. Secomb, R. Hsu, R.D. Braun, J.R. Ross, J.F. Gross and M.W.
Dewhirst, *"Theoretical simulation of oxygen transport to tumors by three-dimensional
networks of microvessels"*, in *Oxygen Transport to Tissue XX*, ed. A.G. Hudetz and D.F.
Bruley, Plenum, New York, 1998, pp. 629-634.

The file header, verbatim from the file:

```
 Segment   Node Node  diam. length        node from           node to
 No. name  from to                      x     y     z      x     y     z
```

### 1.2 One correction to the scouting survey, and one thing it got exactly right

**The scout recorded "106 segments". That is the LINE count.** The file is 2 header lines
plus **104 numbered segments**, numbered 1..104 contiguously. Everything below uses 104.
The error is small and does not affect the route; it is recorded because an uncorrected
count would have propagated into a density.

What the scout got right: the file is open, it is a real mapped tumour network, and it has
node coordinates *and* per-segment diameters, which is what this test needs. Verified
independently here: the straight-line distance between the two node coordinates of each
segment agrees with the segment's stated `length` column to a **maximum relative
difference of 2.25 × 10⁻⁴**, so each row really is a straight cylinder and the geometry is
internally consistent.

### 1.3 What I checked on the anchor side, and one thing that is not knowable

Marien 2016 (PMID 27583442) — the source of **both** failing anchors — was fetched in full
from PLOS ONE. Its two definitions, verbatim:

> "A grid point, which was designated by two perpendicular cross-lines, was regarded as
> overlapping a microvessel when it fell on an endothelial cell or a vessel lumen."

> "every CD31-positive object, no matter how small, should be counted, except suspected
> CD31-positive monocytes, macrophages and tumor cells. Furthermore, if CD31-positive
> objects were connected, they were considered a single object, while absence of staining
> defined two or more separate objects."

**Marien 2016 does not state a section thickness anywhere in its text.** I searched the
full fetched article. This matters, because section thickness is one of only three ways a
sectioning effect can survive into a model-vs-anchor ratio (§0.1). The 4 µm used below as
the histology comparator is the thickness quoted verbatim from the *other* (Tamma 2022,
PMID 36611806) anchor's Methods in `results/vascular_segmentation_audit.md` §2.1; I did
not re-fetch Tamma in this pass, MDPI returned 403 and PMC returned a bot wall, so that
one number is a repo-internal quotation rather than a first-hand one. The full thickness
sweep (0 / 4 / 5 / 15 µm) is reported so any other thickness can be substituted.

Marien also anticipates this route's question in its own Discussion:

> "Depending on the orientation of the vessel and the direction of the sectioning, one
> vessel can appear as separate shapes in the two-dimensional section."

That is a bias on the **measurement** side, in the same direction as the model's. It is
another reason the two sides do not simply cancel — and, as §4 shows, it is not nearly
large enough to matter.

---

## 2. The 3-D truth, and one property of this network that shapes everything else

| quantity | value |
|---|---|
| segments | **104** |
| total vessel centreline length | **7,465.7 µm** |
| reference volume (stated region) | 550 × 520 × 230 = **6.578 × 10⁷ µm³** |
| **vessel volume fraction `V_V`** (union, 1 µm voxels) | **1.4999 %** (naive Σπr²L: 1.5134 %) |
| **vessel length density `L_V`** | **113.49 mm/mm³** → stereological `Q_A = L_V/2` = **56.75 /mm²** |
| diameter | min 4.70, median 8.30, max 33.2 µm; **length-weighted mean 11.17 µm** |
| **median 3-D distance from tissue to nearest vessel SURFACE** | **42.40 µm** (p05 7.9, p95 99.3, max 164.4) |

**The network is strongly anisotropic, and this is the single most important caveat in
this document.** Length-weighted mean |unit-axis component|, against an isotropic
reference of 0.5 for each:

| axis | measured | isotropic |
|---|---:|---:|
| x | 0.588 | 0.5 |
| y | 0.659 | 0.5 |
| **z** | **0.108** | 0.5 |

The vessels lie almost entirely in the x–y plane. The node z-coordinates span 25–213 µm
against a stated region depth of 230 µm, and 39 of 104 segments have **exactly zero**
z-extent. This is what a **dorsal skin-flap window chamber** is: a thin, flat sheet of
tissue observed *en face* through a window, with far coarser resolution in depth than in
plane. Some of the flatness is the tumour; some of it is the preparation and the imaging
axis. **I cannot separate the two, and no conclusion below leans on the anisotropy being
biological.** Where a result depends on section orientation it is reported as a
sensitivity, never as the headline.

---

## 3. IMPLEMENTATION CONTROLS — PASS, and one change to how they were computed

| control | result | band | verdict |
|---|---:|---|---|
| **C1** Delesse, `A_A / V_V`, per normal | median **0.9981**, worst deviation **0.0045** | [0.95, 1.05] | **PASS** |
| **C2** crossing-integral recovery, per normal | median **0.9997**, worst deviation **0.0067** | [0.90, 1.10] | **PASS** |
| box-volume recovery, per normal | median **1.000000** | — | exact |
| C1 in the **as-pre-registered** random-plane form | **0.9706** | [0.95, 1.05] | PASS |
| C2 in the **as-pre-registered** random-plane form | **0.9623** | [0.90, 1.10] | PASS |

**Disclosure: the controls were computed in a stronger form than §0.5 specified, and the
pre-registered form is reported alongside so nothing is concealed.** §0.5 asked for a
ratio-of-sums over the 400 random planes. That construction confounds *is the rasteriser
correct* with *are 400 random offsets enough to average a high-variance integrand*, and on
this anisotropic network it moved by more than 10 % between a 60-plane and a 400-plane
draw. It was replaced by a deterministic per-normal Cavalieri sweep (`ds = 2 µm`, 16
normals), which is strictly stronger — the identity must hold for **every** normal, not
merely on average — and has no sampling noise at all. Both forms pass. **No finding
threshold (T1, T2, T3) was altered, before or after seeing any result.**

Two rasteriser defects were found and fixed *by* these controls, before any finding was
read, and both are recorded because they are the kind of thing that silently manufactures
a result:

1. **The reference ("tissue") area of an oblique section must not be estimated by counting
   pixel centres.** On this box that estimator runs up to **9 % low** for a diagonal
   normal, which propagated into both controls at ~11 % and would have voided the run for
   a reason having nothing to do with the network. It is now computed as the exact area of
   the convex polygon {plane} ∩ {box}, verified against the analytic identity
   ∫ A(s) ds = V_box to 10⁻⁶.
2. **Offsets sampled uniformly over a normal's support must be re-weighted by that support
   width** when planes of different normals are pooled, because the reference box is
   550 × 520 × 230 µm and its support width runs 230 → ~750 µm with the normal.

---

## 4. HEADLINE — the answer, both halves of it

> ### The sectioning bias is real, it is large in SPREAD, and it splits the two failures.
>
> **It cannot explain the 3.59× MVD miss. Not partially, not at the tail.** Over 283
> non-degenerate random sections of a real 3-D tumour network, the **largest** profile
> count any single section produced was **1.99×** the ensemble mean, and the fraction of
> sections at or above 3.589× was **0 / 283 = 0.0000**. The pre-registered
> asymmetric-residual statistic `R_total` reached **1.769** against a bar of **2.871**.
> **T3 verdict for MVD: PARTIALLY ACCOUNTS — which is to say, it does not account for it.**
> The 3.59× is not a stereology artefact and the repository's existing attribution —
> single-linkage over-splitting at a 20 µm threshold below the 50 µm endothelial nuclear
> pitch (`results/vascular_segmentation_audit.md` §2.4, §7.3) — stands unchallenged by
> this work.
>
> **For the 2.34× areal fraction it formally passes the pre-registered bar, and the pass
> is weak by construction and gets weaker on inspection.** `R_total` = **3.302** against a
> bar of **1.869**, so **T3 verdict: ACCOUNTS FOR**. But the honest reading of that number
> is: a single 2-D section lands at or above 2.337× the ensemble mean **3.53 %** of the
> time (10 / 283). Sectioning luck *can* produce the held-out section's miss. It is a
> 1-in-28 event. **And the repository has two sections at 2.337× and 2.327×** — two
> near-identical draws from a 3.5 % tail, joint probability **0.125 %**, which is not what
> luck looks like. So: **sectioning cannot be excluded as a contributor to the areal
> fraction miss, and it is not a sufficient explanation for its reproducibility.**
>
> **The finding that survives both halves, and is the one worth acting on: EVERY vascular
> statistic tested FAILS the single-section spread threshold T2, most of them by a wide
> margin.** The areal fraction varies **8.3×** between the 5th and 95th percentile of
> random sections *of one unchanging tissue*. This engine has one section. That spread is
> not an error bar it currently carries anywhere.

---

## 5. RESULTS

### 5.1 T1 — systematic bias (IUR sections, t = 5 µm, the Xenium thickness)

`B` = median over sections of (2-D estimate / 3-D reference). BIASED-UNSAFE outside
[0.667, 1.500].

| statistic | 3-D reference | 2-D median | **B** | T1 |
|---|---:|---:|---:|:--|
| areal fraction (%) | `V_V` = 1.500 | 1.740 | **1.160** | safe |
| profile count /mm², axis crossings | `L_V/2` = 56.75 | 53.51 | **0.943** | safe |
| profile count /mm², Marien connectivity rule | `L_V/2` = 56.75 | 50.36 | **0.887** | safe |
| **apparent vessel caliber (µm)** | 11.174 (length-wtd true) | 18.185 | **1.627** | **UNSAFE** |
| **distance to nearest vessel, 2-D vs 3-D, paired** | 1.000 | 1.198 | **1.198** | safe |

**Only one statistic carries a systematic bias past the threshold, and it is the vessel
caliber.** A 5 µm section over-states the true tube diameter by **1.63×** — 1.35× from the
oblique-cut ellipse alone at `t = 0`, rising to **2.10×** at a 15 µm effective depth. This
is the classical result and it is measured here rather than assumed. It matters for this
engine specifically, because the painted areal fraction is **quadratic** in an assumed
radius (`area ∝ π(r₁k)²`, `vessel_geometry.py`), so a caliber over-statement enters the
areal fraction squared.

**The Delesse result is worth stating precisely, because it is the part of the scout's
hypothesis that is simply false.** Over IUR planes at `t = 0` the **mean** areal fraction
is 1.4417 % against `V_V` = 1.4999 % — a ratio of **0.961**, unbiased to within sampling
error, exactly as Delesse's principle requires. **There is no systematic areal-fraction
bias to be recovered.** What there is instead is skew: the **median** section reads
0.8815 × `V_V`, because the distribution is strongly right-skewed — most sections miss the
vessels and a few catch them lengthwise.

### 5.2 T2 — single-section spread. **All five statistics FAIL.**

SPREAD-UNSAFE if P95/P05 across single sections exceeds 2.0.

| statistic (IUR, t = 5 µm) | P05 | median | P95 | **P95/P05** | CV | T2 |
|---|---:|---:|---:|---:|---:|:--|
| areal fraction (%) | 0.444 | 1.740 | 3.672 | **8.27** | 0.62 | **FAIL** |
| profile count /mm² (crossings) | 18.6 | 53.5 | 94.8 | **5.09** | 0.44 | **FAIL** |
| profile count /mm² (Marien rule) | 18.5 | 50.4 | 84.3 | **4.56** | 0.41 | **FAIL** |
| apparent caliber (µm) | 11.5 | 18.2 | 31.6 | **2.75** | 0.33 | **FAIL** |
| 2-D/3-D nearest-vessel ratio | 1.041 | 1.198 | 2.200 | **2.11** | — | **FAIL** |

**This is the result of the route.** Every vascular number this product reports is a draw
from a distribution this wide, and the product reports one draw with no interval on it.
For calibration, Marien's own cohort SDs imply CV 0.46 (Q_A) and 0.51 (A_A) *across
patients*; the CVs above — 0.41 and 0.62 — are **within a single unchanging tissue**,
from nothing but the choice of plane. Marien's paper independently reports A_A as the less
reproducible of its two parameters (inter-ROI ICC below 0.7), which is the same finding
from the measurement side.

`results/dimensionality.md` §1.4 already found the analogous thing for delivery — mean
2-D/3-D exposure 0.880–1.408, *"the scatter is which plane you cut"*. This route extends
that from a solved field to the **geometry statistics themselves**, on a real 3-D network
rather than on a synthetic one.

### 5.3 Section thickness — where the over-projection actually lands

Median over IUR sections, as a ratio to the mathematical plane (`t = 0`):

| statistic | t = 0 | t = 4 µm | t = 5 µm | t = 15 µm | t=5 / t=4 |
|---|---:|---:|---:|---:|---:|
| areal fraction (%) | 1.322 | 1.638 (1.24×) | 1.740 (1.32×) | 2.620 (**1.98×**) | **1.062** |
| profile count /mm² | 49.56 | 49.51 (1.00×) | 50.36 (1.02×) | 47.98 (0.97×) | **1.017** |
| apparent caliber (µm) | 15.13 | 17.65 (1.17×) | 18.19 (1.20×) | 23.51 (**1.55×**) | **1.030** |

Two things follow, and they pull in opposite directions:

* **Thickness moves the areal fraction a lot in absolute terms** — a 15 µm effective depth
  nearly doubles it. That is exactly the `t_eff` ambiguity `results/vessel_geometry_fix.md`
  already flagged (*"`L_V = n_A/(t·λ)` is linear in the effective depth `t_eff`, and
  `t_eff` is not known"*), now measured on a real network instead of argued.
* **But between the two thicknesses actually in play it is nearly nothing.** 5 µm vs 4 µm
  gives **1.062×** on the areal fraction and **1.017×** on the profile count. The
  thickness asymmetry — one of the only three ways a sectioning bias can survive into a
  model-vs-anchor ratio — **contributes at most 6 %** of either failure. It is not the
  explanation for anything.
* **The profile count is almost thickness-independent**, because over-projection adds new
  profiles at the same rate as it merges neighbouring ones under Marien's connectivity
  rule. That is a non-obvious result and it further weakens any thickness account of the
  MVD miss.

### 5.4 Orientation — large, and deliberately not used

Median over 120 planes per coordinate-axis normal, `t = 5 µm`:

| normal | areal fraction (%) | profile count /mm² |
|---|---:|---:|
| x | 1.797 | 66.9 |
| y | 1.947 | 71.1 |
| **z** (the en-face plane of the window chamber) | 1.781 | **10.5** |
| spread (max/min) | **1.09×** | **6.78×** |

Cutting this network across the grain versus along it changes the profile count by
**6.8×** — far more than the 3.59× under investigation. **I am not using this to explain
anything, for two reasons, and the reasons are the point:**

1. The anisotropy is at least partly an artefact of the dorsal skin-flap preparation (§2),
   so its magnitude is not transferable to an excised human GBM block.
2. Even if it were, the asymmetry only bites if the model's section and the anchor's
   sections are **systematically** differently oriented relative to the vasculature. Both
   are human GBM resection blocks cut without reference to vessel direction. There is no
   mechanism that would orient one and not the other, and inventing one would be exactly
   the sort of post-hoc rescue this route exists to avoid.

The honest use of this row is the opposite of an excuse: **it is the strongest single
demonstration that a profile count on one section is not a property of the tissue.**

### 5.5 The distance to the nearest vessel — the one that matters physically

This is the statistic that sets the length scale in a diffusion model, and it is measured
**paired**: for the same physical points, the in-plane distance to the nearest vessel
profile against the exact 3-D distance to the nearest vessel surface.

| | median | P05 | P95 |
|---|---:|---:|---:|
| true 3-D distance at those points (µm) | 41.9 | 32.6 | 58.8 |
| 2-D in-section distance (µm) | 60.7 | 44.0 | 104.8 |
| **ratio, per section** | **1.198** | 1.041 | 2.200 |

A 2-D section **over-states** how far tissue is from its blood supply, always — the
nearest vessel may be out of plane — by a median of **1.20×** and by up to **2.2×** on an
unlucky section. The bias is one-directional by construction and the measurement confirms
it: the minimum over 283 sections is 1.00, never below.

**This makes one of the repository's rows worse, not better, and that must be said.** The
held-out section reports `median_tumor_to_vessel_um` = **23.87 µm** against a
`PUBLISHED` 58 µm — already **0.41× LOW**. Correcting for the 1.20× sectioning inflation
pushes the implied true 3-D distance to roughly **20 µm**, i.e. **0.34×**. The sectioning
bias does not rescue that row; it deepens it. (Caveat, stated in the same breath: the
58 µm figure is itself *derived* from the 65 /mm² MVD under a Poisson null, not measured,
so the comparison is weak on both sides.)

---

## 6. THE CONNECTION — does this account for the 2.34× and the 3.59×?

### 6.1 The structural point, which was pre-registered before any measurement

**Both sides of both failing comparisons are 2-D section measurements** (§0.1). A bias
that applies equally to both cancels identically. The only parts that can survive are
thickness (§5.3: **≤ 1.06×**), orientation (§5.4: large on this network, but with no
mechanism to make it asymmetric between two human GBM blocks), and single-section luck
(§5.2: large). The pre-registered `R_total` compounds the first and the third, generously.

### 6.2 MVD, 3.589× HIGH — **the hypothesis is REFUTED**

| | |
|---|---:|
| known failure `F` | **3.589×** |
| `R_single` = P97.5 / mean | 1.739 |
| `R_thick` = median(t=5) / median(t=4) | 1.017 |
| **`R_total`** | **1.769** |
| bar to "ACCOUNTS FOR" (`F`/1.25) | 2.871 |
| **pre-registered verdict** | **PARTIALLY ACCOUNTS** |
| **P(one section ≥ 3.589× the ensemble mean)** | **0 / 283 = 0.0000** |
| largest ratio any single section reached | **1.986×** |

**Say it plainly: the sectioning bias does not explain the MVD failure.** The
pre-registered verdict word is "PARTIALLY", because `R_total` = 1.769 clears the 1.5 floor
— but the exceedance calculation makes the situation unambiguous: **no section, out of
283 random sections of a real tumour network, ever reached the miss**, and the ceiling of
what sectioning can do here (1.99×) is *below* even the mildest of the repository's three
MVD misses (1.70× / 3.59× / 4.93×) by more than the gap between them.

This closes a question the repository has carried since 2026-07-30 and removes a live
excuse. `results/vascular_segmentation_audit.md` §7.3 had already ruled out the
merged-object convention as a *complete* quantitative explanation while confirming it in
kind and direction; this route removes the remaining alternative. What is left is what
that audit named: the model's `n_vessels` is not a profile count at all — it is the number
of single-linkage clusters of endothelial **nuclei centroids** at 20 µm, a threshold
provably below the 50 µm axial nuclear pitch the same package asserts — plus a marker gate
that does not perform Marien's exclusion of CD31⁺ monocytes, macrophages and tumour cells.

### 6.3 Areal fraction, 2.336× HIGH — **a formal pass that should not be leaned on**

| | |
|---|---:|
| known failure `F` | **2.336×** |
| `R_single` = P97.5 / mean | 3.109 |
| `R_thick` | 1.062 |
| **`R_total`** | **3.302** |
| bar to "ACCOUNTS FOR" (`F`/1.25) | 1.869 |
| **pre-registered verdict** | **ACCOUNTS FOR** |
| P(one section ≥ 2.336× the ensemble mean) | **0.0353** (10 / 283) |
| P(**two** independent sections both ≥ 2.336×) | **0.00125** |
| largest ratio any single section reached | 3.927× |

The pre-registered threshold is met and is reported as met. **Three qualifications, all of
which were written into §0.6 before the run or follow directly from it:**

1. **`R_total` was built to be generous** — a 97.5th percentile compounded with a full
   thickness effect, with the pre-registration stating in advance that *"a pass under a
   generous test is a weak pass and will be labelled as one."* This is that label.
2. **It is a 1-in-28 event, and the repository has two of them.** Sections 2 and 3 sit at
   2.337× and 2.327× — two draws from a 3.5 % tail, landing within 0.4 % of each other.
   Joint probability under independent sectioning luck: **0.125 %**. Section-to-section
   luck does not reproduce itself to three significant figures. Something systematic is
   doing that, and it is not the plane.
3. **The low tail covers section 1 too, which is the tell.** Section 1 sits at **0.265×**
   the same anchor, and P(one section ≤ 1/2.336× the mean) = **0.124**. So the sectioning
   spread comfortably brackets *all three* of the repository's sections — 0.27×, 2.33×,
   2.34× — while explaining the pattern of none of them. A hypothesis that accommodates
   every observation equally well has not been tested by any of them. **That is the
   correct reading of this "pass", and it is the reason it is not being promoted into a
   recovered result.**

**What this does license, and it is not small:** the areal-fraction anchor comparison
**cannot be scored on one section to better than a factor of several**, and any future
statement of the form "this section's vascular areal fraction is X% and the literature
says Y%" must carry a single-section interval of roughly 0.14×–3.1× around the tissue's
own mean before any model error is discussed. No such interval exists anywhere in
`results/` today.

### 6.4 What changes, and what does not

**Nothing in this document changes any shipped number, any constant, or any verdict.**
No file under `tmesim/` was touched. `results/vascular_segmentation_audit.md` §6 measured
that the vessel segmentation is *"a pure reporting layer: it produces MVD and the area
fraction and feeds no physics"* (the delivery source field is bit-identical across linking
distances 10/20/50/80 µm), so the blast radius of everything above is confined to the two
reported rows and to how they are interpreted.

---

## 7. LIMITS — what this does NOT show

1. **One network. One tumour. One species.** A single R3230AC rat mammary carcinoma in a
   dorsal skin-flap window, n = 104 segments, 0.066 mm³. Nothing here is a cohort. The
   repository's flagship tissue is human GBM, and `results/landscape_pathology_3d.md`
   established — twice, adversarially — that **no 3-D vascular graph of human GBM exists
   publicly**. This route was run because the stereology of sectioning a tube network is a
   geometry question, not a biology question. That argument covers the *mechanism*; it
   does **not** license transferring the *magnitudes* (1.20×, 1.63×, 8.3× spread) to human
   GBM, whose vasculature is more heterogeneous, more tortuous and larger in caliber.
2. **The network is quasi-planar, and I cannot tell how much of that is the preparation.**
   Length-weighted mean |u_z| = 0.108 against an isotropic 0.5. A window chamber is a thin
   flat sheet imaged *en face*. Every orientation-dependent number in §5.4 inherits that
   and none of them is used to support a conclusion.
3. **A 3-D graph is not a Xenium section, and the difference is the whole of the MVD
   defect.** This network gives centreline geometry and diameters. It has **no endothelial
   nuclei**. So I can compute what a *faithful profile count* reports on a 2-D section —
   which is Marien's estimator — but I **cannot** emulate the engine's estimator, which
   clusters detected nuclear centroids at 20 µm. The bias measured here applies to the
   anchor's estimator and to any faithful profile count. It does not and cannot apply to
   the engine's, and that asymmetry is itself part of why the two numbers do not compare.
4. **Straight cylinders.** Each row is a straight segment of constant diameter. Real
   vessels taper and curve within a segment; curvature increases the chance that one
   vessel yields multiple profiles, which would push the profile count *up* and make the
   MVD refutation of §6.2 weaker. The magnitude is unquantified here.
5. **The reference volume is the stated bounding box**, 550 × 520 × 230 µm, which includes
   corners the mapped network does not reach (node z spans 25–213 µm). Using a tighter
   box would raise `V_V` and `L_V` proportionally and would raise the 2-D estimates by the
   same factor, so **ratios are unaffected**; absolute densities are not, and the absolute
   `V_V` = 1.50 % and `L_V` = 113 mm/mm³ should be read with that in mind.
6. **Marien's section thickness is unknown** (§1.3). The 4 µm comparator is borrowed from
   a different paper. The full 0/4/5/15 µm sweep is published in §5.3 so any other value
   can be substituted; at plausible histology thicknesses the conclusion does not move,
   because the profile count is nearly thickness-independent (§5.3).
7. **`P(exceedance)` is estimated from 283 planes**, so `0.0000` for MVD means "below
   roughly 1/283 ≈ 0.0035", not "impossible". The stronger statement is the one from the
   observed maximum: the largest single-section ratio in 283 draws was **1.99×**.
8. **This says nothing about whether the physics is right.** It is a geometry result about
   an estimator. It does not touch the transport law, the oxygen field, any drug, any
   compound-resolved quantity, or any comparison between compounds. No kill fraction, no
   clinical inference, no ranking of any kind is implied or licensed by anything above.

---

## 8. Reproduction, and the calibration guard

```
.venv/bin/python experiments/sectioning_bias.py                # main pass
.venv/bin/python experiments/sectioning_bias.py --exceedance   # supplementary, section 6
```

Seed **20260920**; 400 IUR planes (283 above the 0.05 mm² field floor), 120 planes per
coordinate-axis normal, thicknesses 0/4/5/15 µm, 1.0 µm raster pitch, 400,000 3-D sample
points. All numbers, including every distribution and both control forms, are in
`results/sectioning_bias_data.json`.

**`calibrate_alpha()` at START = `22.223477711172563`** — asserted equal to
22.223477711172563; assertion passed.
**`calibrate_alpha()` at END = `22.223477711172563`** — asserted equal to
22.223477711172563; assertion passed.
The assertion is in the script itself (`assert_alpha`), runs at both ends of both entry
points, and raises rather than warns.

**Files written by this work, and no others:** `experiments/sectioning_bias.py`,
`results/sectioning_bias_2026-09-20.md`, `results/sectioning_bias_data.json`. Nothing
under `tmesim/`, `app/`, `landing-3d/` or any pre-existing file in `results/` or
`archive/` was read-modified-written, and no constant anywhere was tuned.
