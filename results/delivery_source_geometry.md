# The delivery source runs on the geometry the code calls "Wrong #2" — and that is the right geometry for it

*Written after tracing both vascular representations from measured endothelial cell
positions to the PDE source, verifying the production field's link-distance invariance
by experiment on all four cached real sections, and measuring what a switch would move
on the nine published rows (seven of which carry a verdict). No radius, bandwidth or
linking distance was changed. The
only edit under `tmesim/` is a DOCSTRING in `reconstruct/vessel_geometry.py` recording
the finding; no constant, default or returned value moves, and `pytest -q` is unchanged
at the declared 754.*

> ### Headline
>
> **The audit's claim is confirmed and generalised: the production `vascular_density`
> is bit-identical across linking distances 10/20/50/80 µm on ALL FOUR cached real
> sections**, and the same harness *does* detect a change on the two paths that read
> cluster labels — so "identical" is a property of the default, not of the test.
>
> **But routing the delivery source through the corrected sub-grid geometry is the
> WRONG fix, and the reason is a physical law, not a tolerance.** The engine's delivery
> term is `vascular_density · bbb_penetrance · plasma` — a transvascular flux, whose
> geometric weight is the vessel **wall** area (`S/V`). The sub-grid field paints a
> **cross-section** area. With this module's own `r = r₁·k`, a cross-section area goes
> as `k²` while wall area goes as `k`, so the sub-grid field's source deposited per
> endothelial cell is **not constant — it spans 8.1× to 8.4× across `k` on every one of
> the four sections** (30.4 µm²/EC at `k = 1`, ~250 µm²/EC at `k ≥ 5`). The production
> binary field deposits a constant per endothelial cell, which — given this module's
> own two measured human constants — **is** the wall-area law.
>
> **So the two geometries are answers to two different questions, and the module's own
> docstring already scoped it that way**: the sub-grid geometry is "the default for the
> *reported* area fraction, which is the quantity being fixed". A stereological areal
> fraction is a cross-section quantity (Marien counts a grid point vascular when it
> falls on "an endothelial cell **or a vessel lumen**"). A leak source is not.
>
> **The areal-equivalence result the brief anticipated is confirmed, and it is the
> weaker half of the story.** Rasterised, per endothelial cell, in the order
> s1 / s2 / s3 / breast: the sub-grid field deposits **80.9 / 170.2 / 144.5 /
> 151.5 µm²** against the production field's effective **174.3 / 155.6 / 162.2 /
> 159.9 µm²** — ratios **0.464 / 1.094 / 0.891 / 0.947**. On the three dense sections
> the two fields carry the same total source to within 11 %.
> **But their supports are not the same**: the sub-grid field is nonzero on **2.0–2.2×
> as many grid cells**, Jaccard **0.44–0.53**, and it moves the median tissue-to-source
> distance by **−2.1 to −7.0 µm**. Same dose, different place. That is precisely the
> perturbation the `PS` gauge cannot absorb.
>
> **And on the module's own fidelity score the "fix" is worse.**
> `vessel_geometry.vascular_heterogeneity` exists to ask whether a field "puts its
> vasculature *where the cells actually are*". The production binary field scores a
> heterogeneity ratio of **0.94–1.01** (1.0 is perfect) and Pearson **r = 0.94–0.96**
> against the measured endothelial counts. The sub-grid field scores **1.45–2.26** and
> **r = 0.83–0.93** — it has *sharpened* the architecture by exactly the `k²`
> over-dispersion the audit diagnosed in §3.
>
> **`alpha` cannot absorb any of this, and the "refit" arm is degenerate — proved, not
> assumed.** `calibration.calibrate_alpha(target_efold_um, drug, ecm_density)` takes no
> microenvironment; its three inputs are `D_eff = 140.408947 µm²/s`,
> `k_passive = 6.308002e-05 /s` and `uptake_proxy = 2.237458e-04 /s`, none of which is
> a function of `vascular_density`. Re-solving it under a swapped source returns
> **148.60815765445545 — a difference of exactly 0.0**. Unlike the `D_tissue` sweep,
> this perturbation has **no gauge to hide in**.
>
> **Blast radius, measured over 27 transport solves: THREE of the seven SCORED
> verdicts flip.** GBM temozolomide **PASS → FAIL** (the flagship — the only drug with a
> matched patient-derived ex-vivo anchor), GBM cisplatin **PASS → FAIL**, and breast
> docetaxel **FAIL → PASS** — the project's *only* falsification, removed. Two more
> survive on margins cut **17.7×** and **13.7×**. The `binary` control arm reproduces all
> nine published ratios exactly, so nothing here is the harness. (Two of the nine
> published rows — breast paclitaxel and breast 5-fluorouracil — carry no verdict at
> all: one drug-matched study each, hence a zero-width band. Their ratios reproduce and
> move like the rest; see §5.0.)
>
> **And it is ~100 % SHAPE, not dose.** Holding the tissue-integrated source
> **bit-exactly** at the production total (`match_total_to`) restores mean tumour
> exposure to **1.02–1.03×** on GBM and 0.96–0.98× on breast, and the scored ratios stay
> at **0.36–0.48×** and **0.52–0.76×**. On GBM temozolomide the dose-matched arm is
> **worse** than the unmatched one. A gauge that exactly restores the dose removes about
> **1 %** of the movement on breast and **none** of it on GBM.
>
> **No `alpha` absorbs it either, even with the anchor abandoned.** Tuning `alpha` freely,
> the best the sub-grid source reaches is **0.2053** against a published 0.2301 — bought
> by pushing doxorubicin's penetration length to **129.5 µm, 1.80× above the top of its
> measured 57.7–72.1 µm band**.
>
> **The O(N²) landmine is real, there are TWO of them, and the audit named one.**
> Measured on all four sections: at 50 µm section 3's giant component reaches
> **67.2 %** of all endothelium, and `segment_vessels` would then allocate **~60 GB**
> in the extent step *and* **~60 GB** in the inter-vessel-distance step
> (`cKDTree.query(xy, k=counts.max()+1)`), which the audit did not name. At 80 µm
> section 3 needs **~132 GB** and **~89 GB**. **I declined to patch
> `vessel_network.py`** (it is outside my ownership) and instead shipped a verified
> drop-in in `experiments/delivery_source_geometry.py`, proved **bit-identical**
> (labels, counts and extents to 0.0 µm) against the shipping function on all four
> sections.
>
> **RECOMMENDATION: do not promote `geometry='subgrid'` into the delivery source.**
> Close the open question the audit left, and record the reason as a *law* mismatch
> rather than a tolerance. The 7.5 µm radius is a gauge constant that `PS` absorbs
> exactly; the `k²` scaling is not, because it changes vessels' weights *relative to
> each other* within a section.

---

## 1. The two chains, end to end

Both chains start from the same `CellPopulation` and the same single-linkage
segmentation, and they diverge at the painting step. Only one of them reaches the PDE.

### 1.1 Chain A — the REPORTED morphometry (sub-grid, per vessel OBJECT)

```
CellPopulation (measured Xenium centroids, marker-gated CellType.ENDOTHELIAL)
  │
  ├─ vessel_network._endothelial_xy(cells)                     -> xy (n_endo, 2) µm
  ├─ vessel_network._single_linkage_labels(xy, LINK_DISTANCE_UM = 20 µm)
  │     cKDTree.sparse_distance_matrix -> connected_components   -> labels, n_vessels
  ├─ vessel_network.segment_vessels(...)  -> {'xy','labels','cells_per_vessel',
  │                                           'extent_um', 'mvd_per_mm2', ...}
  │
  └─ vessel_geometry.vessel_geometry(vessels, radius_mode='data')
        axial_units = max(1, extent / EC_AXIAL_PITCH_UM)      EC_AXIAL_PITCH_UM = 50 µm
        k           = max(1, cells_per_vessel / axial_units)  (uncapped by default)
        r           = SINGLE_EC_VESSEL_RADIUS_UM · k          r₁ = 3.1125 µm
        axis_len    = max(0, extent − 2r)
        footprint   = stadium:  π r² + 2 r · axis_len
          │
          ├─ vessel_geometry._stadium_disc_centres      -> union of discs
          ├─ vessel_geometry._subgrid_coverage          -> fractional coverage on a
          │    FIXED PHYSICAL sub-lattice (pitch = clip(r_min/4, 0.25, 1.0) µm)
          ├─ vessel_geometry.subgrid_coverage_field     -> (ny, nx) in [0, 1]
          └─ vessel_geometry.vascular_area_fraction     -> area_fraction_pct
                └─ vessel_network.segmentation_report rows
                   'vascular_area_fraction_pct'  (vs 2.1 %, NOT commensurable)
                   'vessel_areal_fraction_pct'   (vs 4.70 %, commensurable)
```

**Consumers of chain A, exhaustively** (every `import` of `vessel_geometry` in the
repository):

| consumer | what it does with it | on a physics path? |
|---|---|---|
| `vessel_network.segmentation_report` | the two scored morphometry rows | no — writes a report |
| `vessel_network.discrete_vascular_field(geometry='subgrid')` | the opt-in field | **never the default anywhere** |
| `validation/endothelial_density.py` | `ec_per_cross_section` in the linking sweep | no — writes `results/endothelial_density.md` |
| `experiments/vascular_segmentation_audit.py` | the audit's diagnostics | no |
| `experiments/delivery_source_geometry.py` | this report | no |

**No solver, no engine path, and no published delivery, PK or kill number reads it.**

### 1.2 Chain B — the DELIVERY SOURCE (binary, per endothelial CELL)

```
reconstruct.real_pipeline._reconstruct_xenium
  └─ real_fields.vascular_density_from_cells(domain, cells)        mode='discrete'
       ├─ vessel_network.segment_vessels(cells, domain)            [labels computed]
       └─ vessel_network.discrete_vascular_field(domain, vessels)
            geometry='binary', amplitude='unit', radius_um = VESSEL_RADIUS_UM = 7.5
              tree = cKDTree(vessels['xy'])            <-- ALL endothelial cells
              hit  = tree.query_ball_point(grid_centres, r=7.5, return_length=True) > 0
              field = hit * tissue_mask                             values in {0.0, 1.0}
                                                       *** labels are NEVER read ***
  └─ micro.fields.vascular_density
       ├─ physics.oxygen.solve_oxygen        -> micro.fields.oxygen
       │     └─ cells.update_cell_states     -> HYPOXIC / QUIESCENT labels
       │           └─ kill_calibration.resistance_vector  -> the F_state term
       ├─ physics.pressure.solve_ifp         -> ifp, velocity  (a = vasc · Lp_eff)
       └─ sim.engine.delivery_source(work, drug, plasma)
             source = vascular_density · bbb_penetrance · plasma      (µM/s)
               └─ physics.diffusion.simulate_drug_transport -> AUC
                     └─ every exposure, coverage, penetration-depth and kill number
```

**Two properties of chain B decide everything below.**

* **`amplitude='unit'` makes the field an INDICATOR, not an area.** Its value is
  exactly 1.0 inside the disc, and `discrete_vascular_field`'s own docstring fixes the
  referent: *"1.0 means 'a perfused vessel wall leaking at the reference rate'"*. The
  tissue-integrated source is therefore proportional to the **painted area**, and the
  painted area is `≈ n_EC · π r²` minus overlap — **proportional to the endothelial
  cell count**.
* **`geometry='subgrid'` makes the field a FRACTIONAL AREA COVERAGE.** Its value is the
  fraction of the grid cell occupied by a vessel **cross-section**. The
  tissue-integrated source is then `Σ_objects (π r² + 2 r·axis_len)` — proportional to
  `Σ k²`, not to `n_EC`.

Those are two different physical laws wearing the same variable name, and §4.1 is about
which one a transvascular flux wants.

---

## 2. Verification by experiment — the production field does not read cluster labels

Not read off the source: measured, by varying `link_distance_um` and taking a SHA-256
of the field's bytes. **All four cached real sections**, `dx = 12 µm`.

| section | grid | n_endo | link 10 | link 20 | link 50 | link 80 | `default` invariant? |
|---|---|---:|---:|---:|---:|---:|---|
| s1 (dev, `xenium_gbm`) | 424×675 | 4,295 | 4,149 obj | 3,091 | 1,765 | 825 | **YES** |
| s2 (held out, `xenium_gbm_io`) | 905×956 | 57,733 | 43,154 | 20,875 | 3,337 | 285 | **YES** |
| s3 (`gse324268_xenium`) | 1196×880 | 74,515 | 62,637 | 31,551 | 2,347 | 73 | **YES** |
| breast (`xenium_breast_janesick/rep1`) | 458×628 | 13,677 | 11,125 | 5,756 | 1,626 | 438 | **YES** |

On every section the four `hash_default` digests are **the same 16-byte prefix**, the
field is `np.array_equal` to the one `reconstruct_real` built (`max|diff| = 0.0`
exactly), and the tissue-integrated source total is invariant (e.g. breast 15,191.0 at
every link).

**The harness can detect a difference — positive control.** The same loop hashes the
two paths that *do* read `vessels['labels']`, and both move at every step:

| section | `amplitude='per_vessel_cells'` invariant? | `geometry='subgrid'` invariant? |
|---|---|---|
| s1 | **NO** (4 distinct digests) | **NO** (4 distinct digests) |
| s2 | **NO** | **NO** |
| s3 | **NO** | **NO** |
| breast | **NO** | **NO** |

So the invariance is a property of the production **defaults**, exactly as the audit
said, and is now established on four sections instead of one.

**The mechanism, from the code, one line:** under `geometry='binary'` +
`amplitude='unit'` the function builds `cKDTree(vessels["xy"])` over *all* endothelial
cells and asks only *"is any endothelial cell within 7.5 µm of this grid-cell centre"*.
`vessels['labels']` is dereferenced only in the `amplitude='per_vessel_cells'` branch;
`vessels['cells_per_vessel']` and `vessels['extent_um']` only under
`geometry='subgrid'`. Neither is a production default anywhere.

**Consequence, restated precisely.** The entire vessel *segmentation* — the linking
distance, the cluster labels, MVD, the area fraction, everything
`results/vascular_segmentation_audit.md` re-measured — is today a **pure reporting
layer**. It feeds no physics. `segment_vessels` is still called on the production path,
but only its `xy` survives into the field.

---

## 3. What `vessel_geometry.py` says about itself

Quoted verbatim, because the module's own commentary is the primary evidence about
intent.

> **Wrong #1 — the radius is a constant, and the constant was chosen partly from the
> mesh.** ``VESSEL_RADIUS_UM``'s own docstring says it is "the smallest radius that
> reliably paints ≥1 grid cell per vessel while staying inside the published range". A
> parameter justified by the grid cannot then be used to produce a tissue measurement.

> **Wrong #2 — the footprint is painted per endothelial CELL, not per vessel.** This is
> the defect that actually moves the number. The painted area is essentially
> ``n_EC × π r²``, so the reported area fraction tracks the **endothelial-cell
> density**, which is a property of the marker panel and the typing threshold, not of
> the vasculature. On the two production sections the vessel densities are nearly
> identical (111 vs 104 objects/mm²) while the endothelial-cell densities differ by
> 2.7× (143 vs 382 EC/mm²) — and the reported area fractions differed by exactly that
> ratio (2.50 % vs 5.79 %). A vessel profile lined by ``k`` endothelial cells is
> **one** vessel, not ``k`` of them.

There is also a third numbered item, and it is a *retraction* rather than a defect:

> **NOT wrong (measured, see results/vessel_geometry_fix.md §1).** The binary
> centre-in-disc rule was suspected of making the area fraction an artefact of the mesh
> … That reasoning is wrong for an *aggregate* quantity. Point-sampling a fixed
> physical region on a regular lattice is an **unbiased** estimator of that region's
> area at any pitch … What *is* grid-dependent under the binary rule is **which vessels
> survive**: at ``dx = 30 µm`` the mean painted-pixel count is 0.26 per vessel object,
> so ~three quarters of the segmented vessels contribute no source at all and the
> field's *geometry* — not its integral — is a function of the mesh.

And the module states its own scope, which is the single most load-bearing sentence for
this report:

> Nothing in ``core`` / ``physics`` / ``sim.engine`` / ``molecular`` is touched, and the
> legacy binary field in ``vessel_network`` is left bit-identical so that no
> already-published result silently moves. The corrected geometry is opt-in at the field
> level (``discrete_vascular_field(..., geometry='subgrid')``) and **is the default for
> the *reported* area fraction, which is the quantity being fixed**.

Finally, the module documents the limitation that §4.1 turns into the decisive argument:

> **Documented limitation.** Sargent 2023 also reports that, *within* the capillary
> zone, two-endothelial-cell capillaries are only ~0.4 µm wider than one-cell
> capillaries — much less than the factor 2 that strict circumferential tiling predicts
> — and states explicitly that endothelial number "alone is insufficient to explain the
> full range of capillary diameters". The ``r ∝ k`` model therefore **over-disperses
> radii at the capillary end**.

**Reading these together.** "Wrong #2" is an argument about a **reported areal
fraction**, and against that quantity it is correct: an areal fraction must count each
vessel profile once. The module never claimed the *field* should switch — it says the
opposite. What the audit found is therefore not a fix that was written and forgotten;
it is a fix that was **correctly scoped** to the reporting layer, whose scope had never
been justified in the direction of the physics. §4 supplies that justification, and it
comes out on the side of not switching.

`VESSEL_RADIUS_UM`'s own justification has meanwhile been **rewritten** in
`vessel_network.py` and the mesh-derived sentence "Wrong #1" quotes is marked
`SUPERSEDED`. The radius now stands on cited in-vivo human morphometry — mean cortical
grey-matter caliber 37.4 ± 12.9 µm diameter (PMID 19856441), of which "7.5 µm is the
**low** end" — and is described as *"the outer wall radius of the perfused vessel
profile, i.e. the annulus of tissue in immediate contact with leaking endothelium"*.
That is a delivery-source referent, and it is not a mesh referent.

---

## 4. Correctness, feasibility and safety — three questions, three different answers

### 4.1 (a) CORRECTNESS — **NO.** The sub-grid geometry implements the wrong law for a flux

The engine's delivery term is

```
source = vascular_density · bbb_penetrance · plasma            (µM/s)
```

which is a **transvascular flux per unit tissue volume**. The physical form of such a
flux is `PS·(S/V)·(C_p − C_i)`: its geometric weight is the vessel **wall area per unit
tissue**, `S/V`. It is not a cross-sectional area and not a volume fraction.

**Endothelial-cell count is a direct wall-area proxy — by this module's own two measured
constants, in both regimes.**

* *Axially*: `EC_PER_MM_VESSEL = 20` (H01 human temporal cortex, ~4,604 endothelial
  cells lining ~230 mm of vessel, PMID 38723085) ⇒ one endothelial cell lines
  `EC_AXIAL_PITCH_UM = 50 µm` of wall. Wall area ∝ EC count.
* *Circumferentially*: the module's own tiling model, `2πr = k·w` ⇒ `r ∝ k`. A
  cross-section lined by `k` cells has perimeter `k` times that of a 1-cell vessel, so
  its wall area per unit length is also `∝ k` = ∝ EC count.

**Both arrangements give `S ∝ n_EC`.** So a source that deposits a *constant* per
endothelial cell implements the wall-area law, and a source that deposits `π r² ∝ k²`
per *object* (of which there are `∝ 1/k`) deposits `∝ k` **per endothelial cell** — it
implements a cross-section-area law.

Measured, per section, binned by `k`:

| `k` bin | s1 (dev) | s2 (held out) | s3 | breast |
|---|---:|---:|---:|---:|
| `k = 1` (singletons) | **30.4** | **30.4** | **30.4** | **30.4** |
| `2 ≤ k < 3` | 84.2 | 73.2 | 72.4 | 74.0 |
| `3 ≤ k < 5` | 152.7 | 133.3 | 138.2 | 137.7 |
| `k ≥ 5` | **245.4** | **255.7** | **251.7** | **256.4** |
| **spread (top bin / `k=1` bin)** | **8.06×** | **8.40×** | **8.27×** | **8.43×** |
| production binary, per EC, pre-overlap | 176.7 | 176.7 | 176.7 | 176.7 |

*(sub-grid footprint area deposited per endothelial cell, µm²; `k = 1` is exactly
`π r₁² = π · 3.1125² = 30.43 µm²` by construction.)*

> **An 8.1–8.4× spread in the per-endothelial-cell leak weight, reproduced to within
> 4 % across four independent sections, two tumour types and three gene panels.** The
> production field's spread is 1.00× by construction. A leak source that weights a
> `k = 7` vessel's endothelium 8.4× more heavily per cell than a `k = 1` vessel's is
> asserting that the same wall leaks 8.4× harder because the linking rule put more
> nuclei in one object — and `k` is an *uncapped* quantity that the audit already showed
> is driven by the section's endothelial density and by which side of the 20 µm linking
> distance its nearest-neighbour distance falls on.

**This is the `k²` defect of `vascular_segmentation_audit.md` §3, transferred verbatim
into the physics.** There it produced a reported area fraction that ran 0.27× low on one
section and 2.34× high on two others. In the delivery source it would produce the same
dispersion in *relative vessel leakiness within a single section* — which, unlike an
overall scale, nothing downstream can undo.

**Independent confirmation from the module's own fidelity score.**
`vessel_geometry.vascular_heterogeneity` was written to ask exactly this: does the field
"put its vasculature *where the cells actually are*", and does it "carry the measured
unevenness"? Scored against the raw endothelial positions at the documented 100 µm bin:

| section | het. ratio (binary) | het. ratio (sub-grid) | Pearson *r* (binary) | Pearson *r* (sub-grid) |
|---|---:|---:|---:|---:|
| s1 (dev) | **1.010** | 2.262 | **0.943** | 0.827 |
| s2 (held out) | **0.943** | 1.447 | **0.960** | 0.929 |
| s3 | **0.957** | 1.763 | **0.957** | 0.911 |
| breast | **0.949** | 1.481 | **0.962** | 0.915 |

A ratio of 1.0 means the field
carries exactly the measured unevenness. **The production field is at 0.94–1.01. The
"corrected" field is at 1.45–2.26 — it has sharpened the architecture, in the direction
`r ∝ k` over-disperses, on every section.** The module's own documented limitation
predicted this.

**What is genuinely wrong with the production field, and why it does not matter here.**
The binary field's 7.5 µm radius is arbitrary at the ~2× level, and it fixes the
constant of proportionality between EC count and source total. That constant is
**exactly the quantity `PS` absorbs**.
`results/transvascular_exchange_revalidation.md` §5, verbatim:

> **PS and `bbb_penetrance` are perfectly confounded.** Scaling `vascular_density` by
> f gives a BIT-IDENTICAL PS field to scaling `bbb_penetrance` by f. The model cannot
> distinguish a wrong vessel map from a wrong permeability.

So the 7.5 µm is a **gauge constant** — a uniform rescale of the source is, to the
model, the same object. The `k²` scaling is
not, because it re-weights vessels *relative to one another inside one section*, and no
global multiplier can undo a relative re-weighting.

> **The trade, stated plainly: switching would exchange a gauge defect that this project
> has already proved harmless for a physical-law defect that no gauge can remove.**

### 4.2 (b) FEASIBILITY — **YES, but with a caveat that cuts the same way**

**First, the areal comparison the brief flagged — confirmed, with one correction to how
it was stated.** The audit reported that the fixed geometry deposits **170.2 µm² per
endothelial cell** on the held-out section against the legacy **176.7 µm²**, i.e. 96 %,
and concluded that the fix "reinstated in substance the very scaling it was built to
remove". Both numbers reproduce here. But they are not the same kind of number: 176.7 is
the **analytic** `π · 7.5²`, while the field the PDE actually receives is
**de-duplicated on the grid**, so the legacy field's *effective* deposit is lower and
falls with endothelial density:

| | s1 (dev) | s2 (held out) | s3 | breast |
|---|---:|---:|---:|---:|
| endothelial cells /mm² | 153.9 | 645.2 | 756.7 | 354.6 |
| legacy, analytic `π r²` | 176.7 | 176.7 | 176.7 | 176.7 |
| **legacy, as rasterised into the field** | **174.3** | **155.6** | **162.2** | **159.9** |
| **sub-grid, as rasterised into the field** | **80.9** | **170.2** | **144.5** | **151.5** |
| **sub-grid / legacy (both rasterised)** | **0.464** | **1.094** | **0.891** | **0.947** |

Compared like with like, the sub-grid field deposits **9.4 % MORE** than the legacy
field on the held-out section, not 96 % of it. **The audit's conclusion survives the
correction and is slightly strengthened**: in areal terms the two fields are
interchangeable to within ~11 % on all three dense sections. The 0.464 on the sparse
development section is the `k = 1` floor of §4.1 — that section is 81 % singletons, so
almost every object gets the 30.4 µm² minimum.

> **So the brief's hypothesis "if the two representations are numerically near-identical
> in areal terms, switching may buy nothing at all" is CONFIRMED for the integral and
> REFUTED for the field.** Same total, different place, and a different law behind the
> difference.

Turning to feasibility proper: the sub-grid field is defined and non-degenerate
everywhere the binary field is, and then some:

| section | binary support (cells) | sub-grid support | ratio | binary support covered by sub-grid | Jaccard |
|---|---:|---:|---:|---:|---:|
| s1 (dev) | 5,199 | 9,520 | 1.83× | 98.4 % | 0.533 |
| s2 (held out) | 62,372 | 136,065 | 2.18× | 97.6 % | 0.443 |
| s3 | 83,940 | 170,491 | 2.03× | 96.9 % | 0.470 |
| breast | 15,191 | 32,329 | 2.13× | 97.7 % | 0.454 |

So the sub-grid support is very nearly a **superset** of the production support
(96.9–98.4 % containment), and the Jaccard of 0.44–0.53 is driven almost entirely by
the sub-grid field being twice as large, not by the two disagreeing about where vessels
are. Feasibility is not the obstacle.

**The caveat.** Singleton objects — **58.5 / 62.9 / 62.3 / 81.3 % of all objects**
(s3 / s2 / breast / s1) at the production linking
distance — receive `r = r₁ = 3.1125 µm`, a footprint of **30.4 µm²**, which at
`dx = 12 µm` is **0.21 of one grid cell**. Under the binary rule the same vessel paints
at least one full cell at amplitude 1.0. So the switch would attenuate the source at
isolated vessels by ~5× *relative to* clustered ones. That is the same 8.4× dispersion
of §4.1 seen from the other end, and it is why the two fields carry nearly the same
total while disagreeing about its distribution.

**On the "intentional smoothing" hypothesis the brief raised: it is supported, and it is
a resolution argument, not a physics argument.** `vessel_geometry.py` records that at
`dx = 30 µm` the binary rule paints **0.26 pixels per vessel object**, so ~¾ of vessels
contribute no source at all — the field's *geometry* is a function of the mesh even
though its *integral* is not. At the production `dx = 12 µm` one grid cell has an
equivalent-disc radius of **6.8 µm**, so `radius_um = 7.5 µm` is the smallest value that
reliably keeps every vessel in the source. The sub-grid painter removes that
representability limit by construction (a vessel far below the cell size still deposits
its own area). **That is a real advantage of the sub-grid painter — and it is
separable from the radius law.** See the recommendation in §7.2.

### 4.3 (c) SAFETY — **NO.** Measured in §5

Three of the seven scored verdicts flip, including the flagship (GBM temozolomide,
PASS → FAIL) and the project's only falsification (breast docetaxel, FAIL → PASS); two
more survive on margins cut by 13.7× and 17.7×. Holding the source total bit-exactly
fixed removes ~1 % of the movement on breast and none of it on GBM.

---

## 5. Blast radius — both arms, on the nine published rows (seven scored)

**27 transport solves**, `dx = 12 µm`, `dt = 10 s`, 6 h window, seed 0 — the published
grid, timestep, window, `PDConfig`, p99 bath-matching convention and frozen literature
bands, through `kill_calibration.decompose` / `ratio_test` unmodified. The swapped source
is installed *before* `prepare_tissue`, so oxygen, the cell-state labels and the IFP
field are all re-solved from it.

### 5.0 The control arm, first: this harness IS the published pipeline

| tissue | drug | published ratio | this run, `binary` arm |
|---|---|---:|---:|
| gbm | temozolomide | 0.2626 | **0.2626** |
| gbm | cisplatin | 0.2791 | **0.2791** |
| gbm | doxorubicin | 0.2077 | **0.2077** |
| breast | doxorubicin | 0.2301 | **0.2301** |
| breast | paclitaxel | 0.1503 | **0.1503** |
| breast | docetaxel | 0.1484 | **0.1484** |
| breast | epirubicin | 0.2301 | **0.2301** |
| breast | 5-fluorouracil | 0.3379 | **0.3379** |
| breast | cisplatin | 0.3085 | **0.3085** |

All nine reproduce to four decimals, and the derived margins reproduce
`results/p99_convention.md` §6's margin column exactly (+65.2 / +33.8 / +261.4 / +300.5
/ +265.1 / +877.0 / +1990.8 / +720.9 / +47.8 %). Anything that moves below is the source
geometry and nothing else.

> **Seven of these nine rows carry a verdict; two do not.** Breast paclitaxel and breast
> 5-fluorouracil each have exactly one drug-matched 3-D-vs-2-D study, so the min-max
> interval over that literature is a zero-width point that no continuous prediction can
> land inside. They are reported **NOT TESTABLE**, and they are *not* re-scored against
> the broad cross-tumour-type "solid" band — that band spans 62× and is the same
> non-drug-matched band four other panel drugs are explicitly refused a verdict on
> (`results/kill_calibration_breast.md` §3, `results/zero_width_bands.md`). Their two
> margin figures in the list above (+265.1 % and +720.9 %) are margins to that CONTEXT
> band and are not verdict margins. Their RATIOS are unaffected and are still reported
> and compared across arms below, because a ratio that moves is a real measurement of
> the blast radius whether or not a bar exists for it to cross. This document formerly
> scored both against the broad band and called the result a verdict; see
> `results/fallback_rule_propagation.md`.

### 5.1 The two alpha arms — the refit arm is DEGENERATE, and that is the finding

`results/dtissue_sensitivity.md` needed two arms because a change to `D` moves
`alpha`'s closed form. **A change to the source geometry does not**, and the reason is
structural rather than numerical:

```python
def calibrate_alpha(target_efold_um=DOX_TARGET_EFOLD_UM, drug="doxorubicin",
                    ecm_density=DEFAULT_ECM) -> float:
    ...
    k_target = d_eff / (float(target_efold_um) ** 2)
    alpha = (k_target - k_passive(dp)) / uptake_proxy(dp)
```

It takes **no `Microenvironment`, no `Domain` and no field**. Its three inputs are
`D_eff = 140.408947 µm²/s`, `k_passive = 6.308002e-05 /s` and
`uptake_proxy = 2.237458e-04 /s`, all compiled from the molecule; the fit target is the
analytic `L = sqrt(D_eff/k_eff)`, which contains no source term. Executed under the
swapped field it returns **148.60815765445545**, i.e. **`alpha_refit − alpha_frozen`
= 0.0 exactly** (`results/sweep/delivery_source_geometry_alpha.json`).

> **Both arms therefore run at the same alpha, and they are the same arm.** This is
> materially worse than the `D_tissue` case, not better: there, `alpha` re-solved and
> re-exported the change as a global multiplier
> (`results/alpha_identifiability.md` §1). Here `alpha` cannot even see the change.
> **The one calibrated constant in the model is blind to the geometry of the source it
> is calibrating a penetration length for.**

Because a degenerate arm proves nothing on its own, the second arm is instead the
**strictly stronger** gauge test: `subgrid_matched`, the sub-grid geometry with
`match_total_to` set to the production field, holding the tissue-integrated source
**bit-exactly** at the production total (breast 15191.0, GBM 5199.0 — verified in the
records). Since `PS` and `vascular_density` are perfectly confounded, a pure scale change
is absorbable *by construction*; **whatever survives the matched arm is what no global
gauge — `PS`, `bbb_penetrance`, `alpha` or plasma level — can remove.**

### 5.2 The scored verdicts

| tissue | drug | band | published / `binary` | `subgrid` | `subgrid_matched` | verdict |
|---|---|---|---:|---:|---:|---|
| gbm | **temozolomide** | 0.1590–0.4622 | 0.2626 (+65.2 %) | **0.1550** (−2.5 %) | **0.1255** (−21.1 %) | **PASS → FAIL** |
| gbm | **cisplatin** | 0.2086–1.5004 | 0.2791 (+33.8 %) | **0.1042** (−50.1 %) | **0.1179** (−43.5 %) | **PASS → FAIL** |
| gbm | doxorubicin | 0.0575–1.0000 | 0.2077 (+261.4 %) | 0.0660 (**+14.8 %**) | 0.0737 (+28.2 %) | PASS → PASS |
| breast | doxorubicin | 0.0575–1.0000 | 0.2301 (+300.5 %) | 0.1307 (+127.5 %) | 0.1326 (+130.7 %) | PASS → PASS |
| breast | paclitaxel | *none — n = 1, zero-width* | 0.1503 | 0.0763 | 0.0774 | **NOT TESTABLE** (no verdict in any arm) |
| breast | **docetaxel** | 0.0152–0.1463 | 0.1484 (**−1.48 % over the ceiling**) | **0.0752** | **0.0764** | **FAIL → PASS** |
| breast | epirubicin | 0.0110–1.0000 | 0.2301 (+1990.8 %) | 0.1307 (+1087.7 %) | 0.1326 (+1104.6 %) | PASS → PASS |
| breast | 5-fluorouracil | *none — n = 1, zero-width* | 0.3379 | 0.2530 | 0.2558 | **NOT TESTABLE** (no verdict in any arm) |
| breast | cisplatin | 0.2086–1.5004 | 0.3085 (+47.8 %) | 0.2160 (**+3.5 %**) | 0.2186 (+4.8 %) | PASS → PASS |

*(percentages are margin to the nearer band edge — the lower edge except for docetaxel,
which is scored against its ceiling. Paclitaxel and 5-fluorouracil have no band, so
they have no margin; their ratios move exactly as before and are shown unannotated.)*

> ### **Three of the seven scored verdicts flip, and two more of the survivors lose over 90 % of their margin.**
>
> * **GBM temozolomide: PASS → FAIL.** This is the flagship verdict — the *only* drug in
>   the project with a matched patient-derived ex-vivo anchor (GSE232173, PMID 37620410).
>   It falls from +65.2 % above its band floor to 2.5 % below it, and to 21.1 % below in
>   the dose-matched arm.
> * **GBM cisplatin: PASS → FAIL**, from +33.8 % to −50.1 %.
> * **Breast docetaxel: FAIL → PASS.** The project's **only falsification** is removed.
>   `results/p99_convention.md` §6.1 already ruled this verdict **UNDETERMINED at
>   `dx = 12`**, and warned that *"a proxy-derived correction that happens to remove the
>   project's only falsification is exactly the result to hold to the highest bar"*. A
>   source-geometry change that removes it — while breaking the flagship — is that
>   warning arriving on schedule.
> * **GBM doxorubicin** survives on a margin cut from +261.4 % to **+14.8 %**, and
>   **breast cisplatin** on one cut from +47.8 % to **+3.5 %** — a **13.7×** collapse.
> * **Only two of the seven scored verdicts** (breast doxorubicin, breast epirubicin)
>   keep a margin above 80 %. This bullet formerly read *"only four of the nine"* and
>   named paclitaxel and 5-fluorouracil among the survivors; both of those "margins"
>   were distances to the broad non-drug-matched context band, not to a bar either drug
>   is held to, so neither is a survivor of anything. Removing them makes this finding
>   worse, not better — the set of verdicts that comfortably survive the source-geometry
>   swap is half the size it was reported to be.

### 5.3 It is SHAPE, not dose — the matched arm settles it

| tissue | drug | mean exposure `subgrid`/`binary` | mean exposure `matched`/`binary` | scored ratio `subgrid`/`binary` | scored ratio `matched`/`binary` |
|---|---|---:|---:|---:|---:|
| gbm | temozolomide | 0.474 | **1.021** | 0.590 | **0.478** |
| gbm | cisplatin | 0.473 | **1.020** | 0.373 | **0.423** |
| gbm | doxorubicin | 0.477 | **1.027** | 0.318 | **0.355** |
| breast | doxorubicin | 0.913 | **0.964** | 0.568 | **0.576** |
| breast | paclitaxel | 0.918 | **0.969** | 0.507 | **0.515** |
| breast | docetaxel | 0.918 | **0.969** | 0.507 | **0.515** |
| breast | epirubicin | 0.913 | **0.964** | 0.568 | **0.576** |
| breast | 5-fluorouracil | 0.931 | **0.983** | 0.749 | **0.757** |
| breast | cisplatin | 0.924 | **0.975** | 0.700 | **0.709** |

**Read the GBM rows.** Holding the source total bit-exactly fixed restores the mean
tumour-cell exposure to **1.02–1.03×** of production — the dose is *fully* returned — and
the scored ratio is still at **0.36–0.48×**. For temozolomide the matched arm is **worse
than** the unmatched one (0.478 vs 0.590). On breast the matched arm returns 96–98 % of
the dose and retains **98.4–99.0 %** of the ratio movement.

> **The delivery blast radius is essentially 100 % shape. A gauge that exactly restores
> the dose removes ~1 % of it on breast and none of it on GBM.** That is the direct
> experimental refutation of "`PS` will absorb it", and it is why §4.1 frames the
> objection as a law rather than a constant.

### 5.4 The source change propagates through OXYGEN into the cell-state labels

Not a side effect worth burying: `prepare_tissue` re-solves oxygen from
`vascular_density` and re-labels every cell from it, so the state-resistance term `F_state`
moves too. Tumour-cell state counts (0 = proliferating, 1 = hypoxic, 2 = quiescent):

| panel | mean O₂ | proliferating | hypoxic | quiescent |
|---|---:|---:|---:|---:|
| gbm `binary` | 0.3521 | 2,391 | **0** | 21,421 |
| gbm `subgrid` | **0.1964** | 0 | **13,955** | 9,857 |
| gbm `subgrid_matched` | 0.2826 | 0 | 0 | **23,812** |
| breast `binary` | 0.3457 | 4,154 | 10,071 | 75,560 |
| breast `subgrid` | 0.3267 | 3,192 | **15,776** | 70,817 |
| breast `subgrid_matched` | 0.3328 | 3,729 | 14,544 | 71,512 |

On the sparse GBM section the raw sub-grid arm turns **58.6 % of tumour cells hypoxic
where the production field labels none at all**, because that section sits on the `k = 1`
floor and the sub-grid field carries only 0.464× the source (§4.2). On breast the hypoxic
fraction rises **57 %** (10,071 → 15,776). Correspondingly `F_state` moves 0.5502 →
0.4121 (GBM, and → 0.6758 for temozolomide, which carries the per-drug hypoxia override)
and 0.5063 → 0.4914 (breast).

**Absolute output moves too**: kill at the repo-nominal plasma level falls to
**0.40–0.55×** on GBM and **0.82–0.93×** on breast.

### 5.5 Could ANY alpha absorb it? Measured — no, and the attempt breaks the anchor

§5.1 proves the shipping calibration cannot move. This asks the strictly weaker question:
if one were willing to **abandon** the doxorubicin penetration anchor and tune `alpha`
freely, is there a value that restores breast doxorubicin's published 0.2301 under the
sub-grid source? Six more solves, `alpha` at 0.25× / 1× / 4× the published value:

| `alpha` | analytic `L` | vs the measured anchor | `binary` ratio | `subgrid` ratio | sub/bin |
|---:|---:|---|---:|---:|---:|
| 37.152 (0.25×) | **129.5 µm** | **1.80× above the 57.7–72.1 µm band** | 0.2993 | **0.2053** | 0.686 |
| **148.608 (published)** | **64.9 µm** | the fitted anchor | **0.2301** | **0.1307** | 0.568 |
| 594.433 (4×) | 32.5 µm | 0.56× the band floor | 0.1614 | 0.0827 | 0.512 |

**No swept `alpha` restores the published value.** The best the sub-grid source reaches is
**0.2053 at `alpha = 37.15`** — still **10.8 % below** the published 0.2301, and bought by
pushing the anchor drug's penetration length to **129.5 µm, exactly 2.0× the measured
64.92 µm target and 1.80× above the top of its 57.7–72.1 µm band** (Primeau 2005,
PMID 16361566, `d_half/ln2`). Restoring 0.2301 exactly would require `alpha` below 0.25×,
i.e. `L` well beyond 129.5 µm, which is further outside the band still.

**And the gap widens in the direction `alpha` would have to move to keep the anchor
honest**: the sub-grid/binary ratio runs 0.686 → 0.568 → 0.512 as `alpha` rises. So
`alpha` is not merely unable to absorb the geometry change; it trades against the one
measurement the model is calibrated to, and loses.

> **Both arms, concluded.** *Frozen*: `alpha = 148.60815765445545`, three verdicts flip.
> *Refit at the anchor*: the identical `alpha`, identical verdicts, difference exactly
> 0.0. *Refit against the delivery observable, anchor abandoned*: still cannot reach the
> published value, at the cost of a 2× anchor violation. There is no arm in which this
> change is absorbable.

### 5.6 Why the movement is drug-ordered

The scored ratio's sensitivity tracks penetration length: the short-`L` taxanes move most
(×0.507) and the long-`L` small molecules least (5-fluorouracil ×0.749, cisplatin
×0.700). That is the expected signature of a **source-geometry** perturbation — a solute
whose diffusion length greatly exceeds the inter-source spacing averages over the layout
and cannot see it; a short-`L` solute resolves it. It is also why the effect cannot be
tuned away drug-by-drug: the ordering is set by `L`, and `L` is fixed by the one
calibrated anchor.

---

## 6. The O(N²) memory hazard — two sites, and I declined to patch either

### 6.1 The two sites

`vessel_network.segment_vessels` contains **two** dense allocations that scale with the
largest cluster, not with the section:

```python
# site 1 — per-cluster extent (the one results/vascular_segmentation_audit.md §4.1 names)
d = np.linalg.norm(pts[:, None, :] - pts[None, :, :], axis=-1)     # (n, n, 2) then (n, n)
extent[k] = float(d.max())

# site 2 — inter-vessel nearest-neighbour distance (NOT previously named)
kmax = int(min(max(counts.max() + 1, 2), n_endo))
dists, idx = tree.query(xy, k=kmax)                                # (n_endo, kmax) x 2
```

Site 2 is the larger of the two on the denser sections, because `kmax` is set by the
**giant component's** size while the row count is the **whole section's** endothelium.

### 6.2 Measured, on all four sections

Giant-component fraction and the implied peak allocation (site 1 counts the `(n,n,2)`
difference plus the `(n,n)` norm; site 2 counts the float64 distances plus the int64
indices):

| section | link | giant comp. | % of all EC | site 1 | site 2 |
|---|---:|---:|---:|---:|---:|
| s1 (dev) | 80 µm | 154 | 3.6 % | <0.01 GB | 0.01 GB |
| breast | 60 µm | 4,019 | 29.4 % | 0.39 GB | 0.88 GB |
| breast | 80 µm | 7,108 | 52.0 % | 1.21 GB | 1.56 GB |
| s2 (held out) | 40 µm | 3,723 | 6.4 % | 0.33 GB | 3.44 GB |
| s2 (held out) | 50 µm | 13,891 | **24.1 %** | 4.63 GB | **12.83 GB** |
| s2 (held out) | 60 µm | 30,705 | 53.2 % | **22.63 GB** | **28.36 GB** |
| s2 (held out) | 80 µm | 54,297 | 94.0 % | **70.76 GB** | **50.16 GB** |
| s3 | 40 µm | 24,742 | 33.2 % | **14.69 GB** | **29.50 GB** |
| s3 | 50 µm | 50,040 | **67.2 %** | **60.10 GB** | **59.66 GB** |
| s3 | 80 µm | 74,272 | 99.7 % | **132.39 GB** | **88.55 GB** |

**At the production 20 µm every figure is under 0.7 GB and the hazard is invisible.**
The audit's report of a run killed at 60 µm on section 2 is reproduced: the two sites
together need ~51 GB there. **Section 3 is worse than section 2** and was not previously
characterised — it crosses 10 GB between 30 and 40 µm.

This also explains, mechanically, the audit's §9 runtime note that the section-3 sweep
"did not complete within the budget of this work" above 30 µm.

### 6.3 What I did instead of patching

`tmesim/reconstruct/vessel_network.py` is outside the files I own, and the brief names
it explicitly as freshly landed and not to be touched. **I therefore declined to patch
it**, and instead shipped a verified drop-in inside the file I do own:

* `experiments/delivery_source_geometry._max_pairwise_distance(pts)` — **exact** max
  pairwise distance in **O(n)** memory. A point set's diameter equals its convex hull's
  diameter, so only hull vertices need pairing; a chunked exact scan (256-row blocks) is
  the fallback for degenerate/collinear sets that Qhull rejects.
* `experiments/delivery_source_geometry._segment_safe(cells, domain, link)` — the subset
  of `segment_vessels`' output that the **fields** actually read (`xy`, `labels`,
  `cells_per_vessel`, `extent_um`), on the same single-linkage operator, skipping site 2
  entirely (no field reads `intervessel_nn`).

**Verified bit-identical** against the shipping `segment_vessels` at the production
20 µm on all four sections, as an in-run assertion rather than a claim:

```
s1_dev       labels_identical=True counts_identical=True extent_max_abs_diff_um=0.0 xy_identical=True
s2_heldout   labels_identical=True counts_identical=True extent_max_abs_diff_um=0.0 xy_identical=True
s3_gse324268 labels_identical=True counts_identical=True extent_max_abs_diff_um=0.0 xy_identical=True
breast       labels_identical=True counts_identical=True extent_max_abs_diff_um=0.0 xy_identical=True
```

**The minimal in-place patch, for whoever owns `vessel_network.py`.** Site 1 is a
three-line swap that changes no number (the hull's diameter *is* the set's diameter).
Site 2 needs a different operator, not a smaller one: the correct fix is a
`query_ball_point`-style expanding-radius search for the nearest *cross-label*
neighbour, and it is a behaviour-preserving but non-trivial rewrite. **Neither is
urgent**: nothing on the production path runs above 20 µm, and §2 has just established
that nothing on the production path reads the output at all.

---

## 7. Recommendation

### 7.1 Do NOT promote `geometry='subgrid'` into the delivery source. Close the question.

The open item `vascular_segmentation_audit.md` §6 left — *"the fix this project has
already written but never shipped"* — is answered, and the answer is that it should not
be shipped **into the source**. Not because the change is too disruptive (it is, but
that is not the argument), and not because the two geometries are equivalent (§4.2 shows
they are not), but because **they implement two different geometric laws and only one of
them is the law a transvascular flux obeys** (§4.1). The evidence is:

1. **A law argument from the module's own measured constants**, independent of any
   external target: wall area ∝ `k`, cross-section area ∝ `k²`, endothelial-cell count ∝
   wall area in both the axial and the circumferential regime.
2. **An 8.06–8.43× dispersion in per-endothelial-cell source weight**, reproduced to
   within 4 % on four sections, two tumour types and three gene panels, against 1.00× by
   construction for the production field.
3. **A worse score on the module's own fidelity metric** (`vascular_heterogeneity`
   1.45–2.26 vs 0.94–1.01, target 1.0) — the `r ∝ k` over-dispersion the module already
   documents as a limitation.
4. **A blast radius that is not absorbable by any gauge** (§5): three of the seven
   scored verdicts flip, ~100 % of the movement survives an exact dose match, and no `alpha` restores the
   published value even with the calibration anchor abandoned.

**The two flips in opposite directions are the decisive practical point.** A change that
removed the project's only FAIL would be suspicious on its own
(`p99_convention.md` §6.1 says so in terms). A change that removes it **and**
simultaneously breaks the one verdict backed by a matched patient-derived ex-vivo
measurement is not a fix under any reading.

**A "no change needed" result is the result here, and it is a real one:** it converts a
standing suspicion of a shipped-but-unrouted fix into a documented, measured scoping
decision, and it removes a change that would have flipped the project's only
falsification.

### 7.2 What IS worth taking from the sub-grid work — and what it would cost

The sub-grid **painter** and the sub-grid **radius law** are separable, and only the
second is objectionable. The painter (fractional-area rasterisation onto a fixed
*physical* sub-lattice) fixes a real defect that the binary rule has and that §4.2
confirms: at coarse `dx` the binary centre-in-disc rule silently deletes vessels
(0.26 painted pixels per object at `dx = 30 µm`), so the source's *geometry* is a
function of the mesh even though its integral is not.

A source that kept the wall-area law **and** gained the sub-grid painter would be one
footprint **per endothelial cell**, at a **fixed** radius, rasterised fractionally. That
is a strict improvement on both current options. **The current code cannot express it**:
`vessel_geometry` is per-*object* by construction, and `radius_mode='constant'` still
paints one stadium per object whose length is the object's own extent. Building it is a
new construction, not a flag.

**I did not build it**, deliberately. It would move every published exposure, coverage,
penetration-depth and kill number — the same blast radius §5 measures — and would need
its own revalidation and its own falsifiability pass. Recording it as the right next
step is the correct output of a diagnosis.

### 7.3 The one change I did make, and the one I declined

* **Made:** a scope section in `tmesim/reconstruct/vessel_geometry.py`'s module
  docstring recording that "Wrong #2" is an argument about an areal fraction and does
  not transfer to the delivery source, with the measured dispersion figures. It is
  **documentation only** — no constant, no default and no returned value changes, and
  `pytest -q` is unmoved at the declared 754 (§9). The purpose is that the next reader
  of the audit's §6 corollary finds the answer next to the code instead of re-deriving
  it.
* **Declined:** patching either O(N²) site in `vessel_network.py` (§6.3). Outside my
  ownership, not on any production path, and — per §2 — feeding nothing. A verified,
  bit-identical drop-in for site 1 and a written spec for site 2 are supplied instead.

### 7.4 What I did not touch, and why

`LINK_DISTANCE_UM`, `VESSEL_RADIUS_UM`, `SINGLE_EC_VESSEL_RADIUS_UM`, the sub-lattice
pitch, `DX_UM`, `alpha`, every `DRUG_LIBRARY` snapshot and every existing file under
`results/` are **bit-identical**. In particular the 8.4× dispersion of §4.1 would shrink
if `max_ec_per_cross_section` were capped at Sargent 2023's observed 5–6 cells per
cross-section — and capping it to make the delivery-source argument softer would be
exactly the rigging this project forbids. The cap is left at its default `None`, and the
dispersion is reported at the value the code actually ships.

---

## 8. What could make this wrong

1. **The wall-area argument assumes the leak is uniform along the wall.** If GBM
   vascular permeability is concentrated at large, structurally abnormal vessels — which
   is a real, documented feature of glioblastoma — then a super-linear weight in `k`
   would be defensible. **Nothing in this repo measures that**, `k` is not a measured
   caliber (§4.1), and adopting `k²` because it happens to be super-linear would be
   fitting a law to a hope. The honest statement is that the wall-area law is the
   *default* a flux boundary condition implies, and that departures from it need their
   own measurement.
2. **`amplitude='unit'` is itself unvalidated as an absolute.**
   `discrete_vascular_field`'s docstring says so explicitly, and
   `real_fields.VASCULAR_MODES` records that flipping between the KDE and the discrete
   field changes the unit of every absolute exposure by ~3.4×. This report compares two
   *geometries* under one amplitude convention; it does not validate that convention.
3. **The 8.1–8.4× dispersion is a property of `k`, which is a property of the linking
   distance.** At a different link distance the dispersion would differ in magnitude.
   Its *direction* would not: `area ∝ k²` against `wall ∝ k` holds for every `k > 1`.
4. **`vascular_heterogeneity` scores against endothelial-cell counts**, so a field built
   per endothelial cell has a structural advantage in it. The metric is not neutral
   between the two constructions. It is reported because it is the module's *own* stated
   fidelity criterion and because the sub-grid field's failure on it is large
   (1.45–2.26 against a target of 1.0), not because it settles the question by itself —
   §4.1's law argument does that, and it is independent of the metric.
5. **Four sections, two tumour types, one platform (Xenium), one `dx`.** The per-EC
   dispersion reproduces to within 4 % across all four, which is stronger than the audit
   had, but every section is 2-D and Xenium-typed.
6. **The GBM blast radius is inflated by the `k = 1` floor, and the matched arm is the
   number to quote.** The development section is 81 % singletons, so the raw sub-grid
   source there carries only 0.464× the production total (§4.2) and the oxygen field
   collapses (0.352 → 0.196), relabelling 58.6 % of tumour cells hypoxic. That is a real
   consequence of the switch, but it is a *dose* effect and a `PS` rescale would remove
   it. The verdict flips are quoted from the **dose-matched** arm as well as the raw one
   precisely so this cannot be the explanation — and both arms flip the same three
   verdicts.
7. **`match_total_to` matches the tissue-integrated total, not the per-drug delivered
   dose.** The residual 2–4 % in the matched arm's mean exposure is that difference plus
   the oxygen/state feedback, and it is reported rather than corrected.
8. **This report does not validate the production field as correct.** It establishes
   that the production field implements the right *law* with an arbitrary constant, and
   that the proposed replacement implements the wrong law. Both statements can be true
   while the field remains a poor absolute model of perfusion.

---

## 9. Reproduction

```bash
./.venv/bin/python -u -m experiments.delivery_source_geometry --stage verify   # ~80 s,  no PDE
./.venv/bin/python -u -m experiments.delivery_source_geometry --stage geom     # ~100 s, no PDE
./.venv/bin/python -u -m experiments.delivery_source_geometry --stage alpha    # ~2 s,   no PDE
./.venv/bin/python -u -m experiments.delivery_source_geometry --stage blast    # 27 solves, ~30 min
./.venv/bin/python -u -m experiments.delivery_source_geometry --stage absorb   # 6 solves,  ~6 min
./.venv/bin/python -m pytest -q
```

Deterministic; the `blast` stage checkpoints one JSON panel per line and resumes.
Raw records:

* `results/sweep/delivery_source_geometry_verify.json` — link sweep + field digests
* `results/sweep/delivery_source_geometry_geom.json` — the two geometries, four sections
* `results/sweep/delivery_source_geometry_alpha.json` — the degenerate refit arm
* `results/sweep/delivery_source_geometry_hazard.json` — giant components + allocations
* `results/sweep/delivery_source_geometry_blast.json` (+ `..._blast_dx12.jsonl`)
* `results/sweep/delivery_source_geometry_absorb.json` — the alpha-absorption sweep

**Runtime note, and it is §6's landmine.** The `verify` stage reaches an 80 µm linking
distance on sections where `segment_vessels` would need 50–132 GB. It completes in ~80 s
because it routes around both O(N²) sites (§6.3) — *not* because the hazard has been
fixed. Anyone re-running the audit's sweep through the shipping function will still hit
it.

**Files changed by this work:** `experiments/delivery_source_geometry.py` (new), this
file, and a **documentation-only** scope section added to
`tmesim/reconstruct/vessel_geometry.py`'s module docstring (§7.3). No constant, no
default, no returned value and no published number moves.

**Test result: `754 passed, 8 warnings` — the declared baseline, unchanged.** Expected
rather than lucky: the only edit under `tmesim/` is a module docstring, and §2 has
independently measured the production field to be bit-identical to the one every existing
result was generated against.
