# Cellular volume exclusion — the obstacle population the hindrance law never had

> # ⚠️ CORRECTED 2026-08-14 — THE MERGE RECOMMENDATION WAS SUPERSEDED: THE DEFAULT WAS FLIPPED, AND THE CHAIN IN §1 IS NO LONGER THE CHAIN THE ENGINE RUNS
>
> **This banner is the only edit made to this file.** No number, table or verdict
> below it has been changed or re-run. Everything below was true on
> **2026-08-10**.
>
> **What changed.** Commit `fd2a9e7` (2026-08-14, *"Flip two-route transport
> ON"*) promoted the permeable-cell law to the shipping default. Record:
> `results/two_route_promotion.md`.
>
> **Which statements are now false.**
>
> | where | statement | status |
> |---|---|---|
> | header | *"Nothing is wired into `compile_drug.py`… `DX_UM` and `alpha`'s published value are untouched."* | **FALSE.** `compile_drug.build_profile` now calls `two_route.two_route_diffusivity` on every compile, and `alpha` moved from `148.60815765445545` to `22.223477711172563`. |
> | §1 table | the five-step chain ending *"⇒ `DrugParams.diffusivity` = **184.75 µm²/s**"* for doxorubicin | **SUPERSEDED.** The chain gained a terminal two-route cell term. Read `results/two_route_promotion.md` for the shipping values. |
> | §1 | *"Where an obstacle volume fraction enters… **exactly three places, and all three count matrix**."* | **FALSE.** There is now a fourth, and it is the one this file argued was missing: a **cell** volume fraction, consumed by `two_route_diffusivity` at compile time. |
> | §2(a) | *"No function in the chain accepts a cellular fraction."* | **FALSE** — this was the defect, and it is closed. |
> | §7 | *"**MERGE the physics, do NOT flip the default.**"* | **SUPERSEDED, deliberately and with its gating conditions discharged.** The default was flipped four days later. This file's recommendation was correct on its date; `results/two_route_promotion.md` is the document that records why it was overtaken. |
>
> **What still stands.** The diagnosis this file exists for — that the hindrance
> law had no cellular obstacle population, and that cells are most of a tumour's
> volume — is **confirmed, not corrected**. So is §4.4's disclosure discipline
> and the finding that the correction **deepens rather than rescues** the panel's
> only FAIL; the promotion re-measured that and it deepened again. The cell
> volume fraction derived stereologically in §3.2 is still the derivation the
> shipping law uses.
>
> **One caveat carried forward.** This file's law enters through
> `ecm_mesh.tissue_hindrance_two_route`; the shipping path calls
> `two_route.two_route_diffusivity` from `compile_drug` directly. The wrapper is
> **not** the shipping call site — see `tmesim/molecular/ecm_mesh.py`'s module
> docstring.

**What this file is.** A structural defect in the transport physics, found
independently by two investigations approaching from opposite directions, plus
the physics that closes it, what the corrected law gives against measurement, and
the propagation of that correction through every published breast verdict.

Driver: `experiments/cellular_exclusion.py`. Code changed: `tmesim/molecular/ecm_mesh.py`
only (three new functions, all additive; the default behaviour of every existing
function is unchanged). Nothing is wired into `compile_drug.py`, no frozen
snapshot is touched, `DX_UM` and `alpha`'s published value are untouched.

> **The rule this file was written under.** No parameter was selected because of
> how close it lands to a measured target. The effective-medium law was chosen on
> two properties stated before any comparison (§3.1); the cell volume fraction was
> derived from the repo's own real section by a stereological identity (§3.2). §4.4
> discloses the value that *would* have hit the measured 12-fold exactly, so the
> reader can check that it is not the value used.

---

## 0. Bottom line

| | before | after (Maxwell cell exclusion) | measured |
|---|---:|---:|---|
| hindrance `D_free/D_eff`, tirapazamine, MCL-appropriate (no stroma) | **1.77×** | **6.61×** | **12×** (V79 MCL, PMID 9806526) |
| hindrance, full pipeline at the repo's ECM defaults | 1.77× | **11.73×** | — |
| `D_eff` (µm²/s), MCL-appropriate | **287.3** | **77.0** — *inside the band* | **74–130** (3 cell lines) |
| `D_eff` (µm²/s), full pipeline at repo ECM defaults | 287.3 | 43.4 | — |
| `D_tissue/D_free`, doxorubicin | 0.3997 | **0.0604** | 0.06–0.13 measured for doxorubicin in cellularised tumour tissue |
| `D_tissue` (cm²/s), doxorubicin | 1.404e-06 | **2.123e-07** | 2.9e-07 / 5.0e-07 / 6.7e-07 (three direct measurements) |

The defect is **confirmed**: no cellular volume fraction enters the hindrance law
anywhere (§2). One missing physical term moves both independent signatures in the
right direction by roughly the right amount, which is the strongest available
evidence that the diagnosis is the right one rather than two coincidences.

It does **not** land exactly on either target, and the two targets pull in
opposite directions. Against the MCL hindrance the correction reaches 6.6× of a
measured 12× (a 3.7× improvement on a 6.8× error, leaving 1.8×). Against the
direct tumour-tissue `D_tissue` measurements the correction **overshoots** — 2.12e-07
against a measured 2.9–6.7e-07 — because the pre-existing ECM legs were partly
standing in for cell-scale obstruction (§5.3). Both facts are reported; neither is
adjusted away.

**Blast radius (§5), both α arms, published `dx = 12 µm` grid, 12 solves:
exactly 1 of 12 verdicts moves — docetaxel FAIL → PASS, and only in `frozen`.**
In `recal`, the mode the codebase actually runs, **no verdict moves at all**, and
docetaxel moves *further* outside its band (1.5 % → 3.7 % outside). A 6.6×
diffusivity correction does not rescue the panel's only FAIL; it deepens it.

**Merge recommendation: MERGE the physics, do NOT flip the default.** §7.

---

## 1. The chain, end to end

Read out of live code by `experiments/cellular_exclusion.py --stage chain`, not
from memory. Doxorubicin, at the repo default `ecm_reference = 0.4`:

| step | code | value | **obstacle volume fraction consumed** |
|---|---|---|---|
| hydrodynamic radius | `stokes_einstein.hydrodynamic_radius_nm(mw)` | `r = 0.066·MW^(1/3)` = 0.5393 nm | **NONE** — a function of MW alone |
| free diffusivity | `stokes_einstein.free_diffusivity(r, 310 K)` | 351.31 µm²/s | **NONE** |
| micro-tortuosity | `brownian.simulate_effective_diffusivity` | × 0.6807 (300 nm disks) | `obstacle_fraction = ecm_reference = 0.4` — **ECM** |
| nm mesh size gate | `ecm_mesh.ogston_renkin_hindrance` | × 0.7726 (Ogston partition 0.9226 × Faxén drag 0.8375) | `fiber_volume_fraction = 0.05` — **ECM fibres** |
| **⇒ `DrugParams.diffusivity`** | `compile_drug.molecule_to_drugparams` | **184.75 µm²/s** | — |
| ECM tortuosity coupling | `physics.diffusion._effective_diffusivity` | × 0.7600 = `1 − 0.6·ecm_density` | `ecm_density` — **ECM** |
| **⇒ `D` in the PDE** | solver | **140.41 µm²/s = 1.404e-06 cm²/s** | — |

`D_tissue / D_free = 0.3997`. This reproduces `results/dtissue_sensitivity.md` §1
to four significant figures, from an independent code path, so the chain above is
the chain the engine actually runs.

**Where an obstacle volume fraction enters, and what it counts.** Exactly three
places, and all three count matrix:

1. `brownian.simulate_effective_diffusivity(free_D, obstacle_fraction=…)` — a 2-D
   Poisson field of 300 nm impenetrable disks. `compile_drug` passes
   `ecm_reference` into this slot.
2. `ecm_mesh.ogston_renkin_hindrance(…, fiber_volume_fraction=0.05)` — a ~2 nm
   collagen-scale fibre matrix.
3. `physics.diffusion.ALPHA_TORTUOSITY = 0.6` against the spatial `ecm_density`
   field.

---

## 2. The diagnosis, confirmed before it was fixed

**(a) No function in the chain accepts a cellular fraction.** Printed live
signatures:

```
ogston_renkin_hindrance('solute_radius_nm', 'fiber_radius_nm', 'fiber_volume_fraction')
simulate_effective_diffusivity('free_D_um2s', 'obstacle_fraction', 'particle_radius_nm', ...)
hydrodynamic_radius_nm('mw', 'is_antibody')
compile_drug.build_profile('molecule', 'ecm_reference', 'ic50', 'temperature_K', 'brownian_kw')
```

`--stage chain` asserts this rather than asserting it in prose: it fails the run
if any of those four signatures grows a cell-fraction parameter, so this claim
cannot silently rot.

**(b) The one fraction that does vary the answer is a matrix quantity, and it
excludes tumour cells by construction.** `ecm_reference` is fed from
`ecm_density`, which `reconstruct.real_fields.ecm_density_from_cells` builds as a
**KDE of `CellType.STROMAL` cells only**, used *"as a spatial proxy for
extracellular-matrix density"* (its own docstring). Tumour cells — 89,785 of the
167,780 cells on the breast section — contribute nothing to it. Doubling a
tumour's cell packing at fixed stromal content moves `D` by exactly zero.

**(c) The 300 nm obstacle population cannot be a mis-scaled cell term.**
`brownian.py`'s docstring calls its obstacles *"cells / ECM aggregates"*, but
`OBSTACLE_RADIUS_UM = 0.30`, and a tumour cell is 10–20 µm across — a factor of
~30 in radius, ~3×10⁴ in volume. More decisively, the *number* fed into that slot
is the ECM field, not a cell count. At most it is a partial, mis-scaled proxy;
it is not a cellular volume fraction.

**(d) The one place in the package that does carry a `cellular_volume_fraction`
is not in this chain.** `molecular/ion_trapping.py::cell_to_tissue_basis` takes
one — for a lysosomal partition, and `compile_drug` never imports
`ion_trapping`. `molecular/binding_site_barrier.py` *reads*
`drug_params.diffusivity` but never writes it.

**Verdict: CONFIRMED. Cellular volume exclusion is absent from the hindrance law,
and the two measured signatures need no other explanation.**

---

## 3. The physics added

### 3.1 The law, and why this law was chosen

`tmesim/molecular/ecm_mesh.py::cell_exclusion(cell_volume_fraction)` implements
**Maxwell's effective-medium result in the impermeable-inclusion limit**:

```
D_eff / D_0 = (1 − φ) / (1 + φ/2)
            = porosity (1 − φ)  ×  tortuosity 1/λ², λ² = 1 + φ/2
```

> Maxwell JC. *A Treatise on Electricity and Magnetism*, Vol. 1, Ch. IX, Art. 314
> (Clarendon Press, 1873) — the effective conductivity of a cluster of spheres of
> one conductivity embedded in an infinite medium of another, under a uniform
> far field; the oldest effective-property estimate in the literature. Diffusion
> in a two-phase medium solves the same Laplace problem as conduction, so the
> conductivity result transfers unchanged; setting the inclusion conductivity to
> zero gives the expression above.

**Two properties, both stated before any comparison to data, are why Maxwell and
not something else:**

1. **It is exact in the dilute limit** — derived from the single-sphere solution
   of Laplace's equation, with no adjustable parameter. Its dilute slope is the
   classic `1 − 3φ/2`, asserted in the module self-test.
2. **It is the rigorous upper bound.** Hashin & Shtrikman proved that for *any*
   isotropic two-phase composite whose inclusion phase has zero conductivity, the
   Maxwell expression is the **maximum** attainable effective conductivity.

> Hashin Z, Shtrikman S. "A variational approach to the theory of the effective
> magnetic permeability of multiphase materials." *J Appl Phys* 1962;33:3125–3131.
> [DOI](https://doi.org/10.1063/1.1728579)
>
> Two clauses of that paper's own abstract are what license its use here, and both
> were checked rather than assumed: *"The results are also valid for the dielectric
> constant, electric conductivity, heat conductivity, and **diffusivity** of such
> materials"* — so the transfer from permeability to a diffusion coefficient is the
> authors' own, not an analogy made in this file; and *"For the case of two-phase
> materials, the bounds derived are the most restrictive ones that can be given in
> terms of the phase permeabilities and volume fractions"* — so no sharper bound is
> obtainable from φ alone, which is all this model knows.

So Maxwell is the **least-hindering law the physics permits**. Real tissue can
only diffuse more slowly than this, never faster. That is the property you want
in a term you are adding to close a known gap: every error it makes is
conservative, and it cannot be accused of being picked because it hinders enough.

**Bruggeman is reported, never used.** `cell_exclusion` also returns
`(1 − φ)^{3/2}`, the differential/self-consistent effective medium — the standard
treatment at *high* inclusion fraction, where Maxwell's dilute derivation is out
of its range, and the "Bruggeman relation" (`λ² = ε^{-1/2}`) of the porous-media
literature.

> Bruggeman DAG. "Berechnung verschiedener physikalischer Konstanten von
> heterogenen Substanzen. I. Dielektrizitätskonstanten und Leitfähigkeiten der
> Mischkörper aus isotropen Substanzen." *Ann Phys (Leipzig)* 1935;416(7):636–664.
> [DOI](https://doi.org/10.1002/andp.19354160705)
>
> **Attribution caveat, checked rather than assumed.** The simple `(1 − φ)^{3/2}`
> form is a later *specialisation* universally called the "Bruggeman correlation";
> Bruggeman's own 1935 paper does not state it. See Tjaden B, Cooper SJ, Brett DJL,
> Kramer D, Shearing PR. "On the origin and application of the Bruggeman
> correlation for analysing transport phenomena in electrochemical systems."
> *Curr Opin Chem Eng* 2016;12:44–51.
> [DOI](https://doi.org/10.1016/j.coche.2016.02.006) — which states that Bruggeman
> "did not explicitly address the special case which leads to the simple analytic
> expression relating the porosity to the tortuosity factor." The form is used
> here as a bracket only, so the attribution caveat does not affect any result.

It always lies **below** the Maxwell value (asserted in the self-test, consistent
with the Hashin–Shtrikman bound), and at φ = 0.79 it is 1.6× more hindering. It is
carried through every table below **as the other end of a bracket**, so the reader
can see which side the truth lies on. Reporting a bracket and quoting the
conservative end is the opposite of tuning to a target.

**Composition with the existing ECM legs.**
`ecm_mesh.tissue_hindrance(solute_radius_nm, cell_volume_fraction, …)` returns
`ogston_renkin_hindrance × cell_exclusion`. The product form is sequential
homogenisation, licensed by scale separation: the fibre mesh has a ~6.9 nm pore
and cells are ~10 µm across — three orders of magnitude apart — so the mesh
homogenises into a local interstitial diffusivity that the cell-scale Maxwell
problem then takes as its matrix conductivity. This is the same composition rule
`compile_drug.build_profile` already uses to multiply its micro-tortuosity by the
mesh gate.

**Which basis, and why it matters.** `cell_exclusion` returns both the flux-basis
factor `(1−φ)/(1+φ/2)` (called `hindrance`) and the tortuosity-only factor
`1/(1+φ/2)` (called `tortuosity_factor`), because they differ by `(1−φ)` and
picking the wrong one is a 5× error at φ = 0.79. **The flux basis is the one used**,
on two grounds:

* it is what a slab-flux or tracer experiment measures and reports as "the tissue
  diffusion coefficient" — Hicks, Toley, Weinberg and Qian are all of that kind,
  so it is the commensurable quantity;
* the engine has **no interstitial porosity anywhere**. `physics/diffusion.py`
  solves `∂C/∂t = ∇·(D∇C) − kC + S` with no `α` in the accumulation term, the
  source deposits into the whole voxel, and the exposure metric integrates `C`
  over the whole tissue. The tortuosity-only factor would be the right one *if*
  the porosity appeared elsewhere; it does not, so dropping `(1−φ)` from `D` drops
  volume exclusion from the model entirely — which is the defect being fixed.

### 3.2 The cell volume fraction, derived from data

`ecm_mesh.cell_volume_fraction_from_section(total_cell_area, tissue_area)` applies
**Delesse's principle** (Delesse A, *Procédé mécanique pour déterminer la
composition des roches*, 1847 — the founding result of stereology, `A_A = V_V`):
for a randomly oriented plane through a
statistically isotropic two-phase medium, the expected **areal** fraction of a
phase on the section equals its **volume** fraction. A spatial-omics section that
carries a segmented area per cell therefore measures φ_cell directly — **with no
assumed cell radius and no assumed section thickness.**

Applied to the repo's own real Xenium breast section
(`data_cache/xenium_breast_janesick/rep1`, Janesick 2023, PMID 38114474) — the
section every published breast verdict is computed on — 167,780 segmented cells,
Σ`cell_area` = 3.249e+07 µm²:

| tissue-area definition | area (µm²) | φ_cell |
|---|---:|---:|
| bounding box | 4.118e+07 | **0.7892** |
| 50 µm occupancy mask | 3.823e+07 | 0.8500 |
| 100 µm occupancy mask | 4.152e+07 | 0.7826 |
| 200 µm occupancy mask | 4.256e+07 | 0.7635 |

**Adopted φ_cell = 0.7892.** The mask convention is worth about ±0.05 and the
adopted value sits mid-range.

**Three caveats travel with it and are not optional.**

1. **Segmentation convention.** Xenium's default boundary expands each nucleus
   outward until neighbours meet, so it *tiles* the section and over-states cell
   area where cells are sparse. φ_cell here is an **upper estimate**, and so is
   every hindrance derived from it.
2. **Delesse assumes isotropy.** Breast ducts are not isotropic; a section cut
   along a duct is not a random plane.
3. **Literature cross-check.** Tumour interstitial volume fractions measured by
   three independent methods (rat-tumour series 36–60 %; Kim 2004 MRI 27.5 ± 9.1 %,
   PMID 15334566; Krol 1999 inulin available volume 50 %, PMID 10463619 — all
   compiled in `results/tumour_physiology_pbpk_literature.md` §A) imply
   **φ_cell ≈ 0.40–0.72**, i.e. the data-derived 0.79 sits *above* the literature
   band. §4.3 therefore reports the whole curve, not one point.

---

## 4. What the corrected law gives against the measured targets

Tirapazamine, compiled through the **unmodified** pipeline exactly as
`experiments/anchor_prototype.py` does (same SMILES, MW 178.151 vs published
178.15, same `ecm_reference = 0.4`), so the "before" column is bit-identical to
`results/external_anchor_survey.md` §3.

### 4.1 Hindrance factor `D_free / D_eff` — dimensionless, absolute

| arm | Maxwell (used) | Bruggeman (bracket) |
|---|---:|---:|
| **BEFORE** (published survey §3.2) | **1.77×** | — |
| A. pipeline + cell term, repo ECM defaults | **11.73×** | 18.32× |
| B. **cell term alone** — MCL-appropriate: packed cells, little stroma | **6.61×** | 10.33× |
| **MEASURED**, V79 MCL vs culture medium (PMID 9806526) | **12.00×** | |

**Arm B is the row to quote against Hicks**, because an MCL is packed cells with
little matrix and the repo's ECM legs describe matrix the MCL does not have. It
reaches **6.61× of a measured 12×**. That is a **3.7× improvement** on the
published 1.77×, closing most of the 6.8× error and leaving **1.8×**.

Arm A — the whole existing chain with the cell term added — gives 11.73×, close
to 12. **Do not read that as the result.** Arm A includes an ECM contribution an
MCL does not have, so its proximity to the MCL measurement is partly the wrong
term doing the right work. The honest statement is that the corrected law brackets
the measurement between 6.6× and 11.7× depending on how much matrix you credit
the MCL with, against a before-value of 1.77× that is outside that bracket
entirely.

**A stated gap, not a silent substitution.** φ_cell here is **0.7892, measured on
a human breast tumour section — not on a multicellular layer.** A literature
search for a measured cell (or extracellular) volume fraction of a Wilson/Hicks
MCL was run and **returned nothing usable**; that search did not complete cleanly
and is recorded as an open item rather than as a null. So arm B applies a
breast-tumour packing density to a V79 monolayer stack. If an MCL is denser than
this section — which is the likely direction, since an MCL has no vasculature, no
necrosis and little stroma — then **arm B under-states the hindrance and the
residual 1.8× is an over-estimate**. Obtaining an MCL extracellular-space
measurement is the cheapest available way to sharpen this row, and is listed in
§7.

### 4.2 Effective diffusion coefficient — absolute, µm²/s

| arm | `D_eff` | vs measured 74–130 µm²/s |
|---|---:|---|
| **BEFORE** (published survey §3.1) | **287.3** | **OUTSIDE** (2.2–3.9× over, missing all three cell lines in the same direction) |
| A. pipeline + cell term, repo ECM defaults | 43.4 | OUTSIDE (1.7× under the lowest) |
| B. **cell term alone** (MCL-appropriate) | **77.0** | **INSIDE** |

Measured: V79-171b 74, MGH-U1 130 (PMID 9806526), SiHa 130 ± 5 (PMID 17869669).

Arm B lands **inside** the measured band, 4 % above its lower edge. This was not
aimed at: the same φ_cell and the same law produce 6.61× on the hindrance row,
which is *not* on target. The two rows are the same measurement expressed two
ways and they cannot both be hit by one φ, because the engine's `D_free`
(509.47 µm²/s for tirapazamine) is itself an unvalidated Stokes–Einstein estimate
at an interstitial viscosity of 1.2e-3 Pa·s. Landing inside on the absolute and
1.8× short on the ratio is the signature of a residual `D_free` error, not of a
tuned φ.

### 4.3 Sensitivity to φ_cell — the whole curve, since φ_cell is the input with the widest uncertainty

| φ_cell | Maxwell factor | arm B fold | arm B `D_eff` | arm A fold | arm A `D_eff` |
|---:|---:|---:|---:|---:|---:|
| 0.40 (bottom of literature band) | 0.5000 | 2.00× | 254.7 | 3.55× | 143.7 |
| 0.50 | 0.4000 | 2.50× | 203.8 | 4.43× | 114.9 |
| 0.60 | 0.3077 | 3.25× | 156.8 | 5.76× | 88.4 |
| 0.72 (top of literature band) | 0.2059 | 4.86× | 104.9 | 8.62× | 59.2 |
| 0.7635 (200 µm mask) | 0.1712 | 5.84× | 87.2 | 10.36× | 49.2 |
| **0.7892 (adopted)** | **0.1512** | **6.61×** | **77.0** | **11.73×** | **43.4** |
| 0.85 (50 µm mask) | 0.1053 | 9.50× | 53.6 | 16.84× | 30.2 |

Across the literature band alone (0.40–0.72) arm B spans 2.0–4.9×, i.e. the
conclusion "the correction is large and in the right direction" is robust to
φ_cell, while the specific value is not. **φ_cell is now the dominant uncertainty
in the hindrance law**, and it is a *measurable* one — which is a strictly better
position than the previous state, where the term was absent and no measurement
could constrain it.

### 4.4 The disclosure that shows this was not tuned

To hit the measured 12-fold *exactly* from the cell term alone would require
**φ_cell = 0.8800**. The value used is **0.7892**, derived from the section by
Delesse before the comparison was run, and it is 0.09 below that. The Bruggeman
bracket at the adopted φ gives 10.33× — also not 12. No φ, and no law, was
selected to reach the target.

---

## 5. Blast radius

### 5.1 How the correction reaches the engine

The cell term is **size-independent** — Maxwell's expression contains no solute
radius — so it arrives at the solver as a single multiplicative factor on every
drug's `DrugParams.diffusivity`:

```
scale = cell_exclusion(0.7892)["hindrance"] = 0.1512   (6.61-fold)
D_tissue(doxorubicin): 1.4041e-06  ->  2.1228e-07 cm^2/s
D_tissue / D_free:          0.3997  ->  0.0604
```

That is **exactly the perturbation shape `experiments/dtissue_sensitivity.py`
already sweeps**, so that harness is re-used unmodified, on the **published
`dx = 12 µm` grid**, in both α arms. The shipping-D records at `dx = 12` were
reused from `results/sweep/dtissue_confirm_dx12.jsonl`, so only the new D point
was solved. The harness's own assertion — that every exposure field is bit-distinct
across the sweep — passed, so the override is provably not a no-op.

### 5.2 The two arms, and why both are required

`alpha` is fitted in closed form so doxorubicin's `L = sqrt(D_eff/k_eff)` equals a
**measured** 64.92 µm penetration depth (Primeau 2005, PMID 16361566), which makes
`alpha` nearly proportional to `D`. `results/alpha_identifiability.md` establishes
that `alpha` is a gauge coordinate: it absorbs a change at the anchor and
re-exports it as a global multiplier on every other drug. So:

* **`frozen`** — D corrected, published `alpha = 148.6082` kept. `L` moves as `√D`;
  doxorubicin's `L` falls to **25.2 µm**, far below the measured 57.7–72.1 µm
  e-folding band. This arm **breaks a measurement the model currently reproduces.**
* **`recal`** — D corrected and `alpha` re-solved against the same measured anchor.
  `alpha` falls to **22.47**; `L` is pinned at 64.9 µm by construction. **This is
  what the codebase actually does**, because `calibrate_alpha()` is called at run
  time.

### 5.3 The correction overshoots the direct tumour-tissue measurements — and that is diagnostic

| | value |
|---|---:|
| corrected `D_tissue`, doxorubicin | **2.123e-07 cm²/s** |
| Toley 2013 LS174T 3-D cell mass (PMID 23860772) | 2.9e-07 |
| Weinberg 2007 rabbit VX2 liver carcinoma in vivo (PMID 17900740) | 5.0e-07 |
| Qian 2003 rat liver in vivo (PMID 12932647) | 6.7e-07 |
| geometric mean of the three | 4.60e-07 |

The corrected value sits **1.4× below the lowest** of the three direct
measurements and 2.2× below their geometric mean, having started 4.8× *above* the
lowest. It has crossed the band rather than landing in it.

The ratio tells the same story more sharply. `D_tissue/D_free` goes
**0.3997 → 0.0604**, and `results/dtissue_sensitivity.md` §8 records measured
`D_tissue/D_free` ratios in cellularised tumour tissue of **0.06–0.13 for
doxorubicin** and 0.045–0.083 for tirapazamine in multicellular layers. The
corrected ratio lands exactly on the bottom edge of the measured doxorubicin band
— a genuinely independent corroboration, since nothing about those ratios entered
the calculation.

**Why it overshoots, stated as a defect and not explained away.** Arm A stacks the
new cell term on top of two pre-existing terms that were *partly* standing in for
cell-scale obstruction: `brownian.py`'s 300 nm disk field driven by
`ecm_reference = 0.4`, whose own docstring calls its obstacles *"cells / ECM
aggregates"*, and `ALPHA_TORTUOSITY = 0.6` against `ecm_density`. Together those
contribute a factor 0.6807 × 0.7600 = 0.5173. If they were entirely matrix, arm A
is right and the residual is a `D_free` error. If they were partly cells, arm A
double-counts by up to ~1.9×. **This report cannot resolve which**, and it does
not guess: the size of the discrepancy is stated so the next piece of work can
resolve it against a measurement. It is the first item in the recommended order of
work (§7).

### 5.4 A resolution hazard the `frozen` arm creates and the `recal` arm does not

`results/grid_convergence_breast.md` established, and tested rather than asserted,
that the ratio's grid sensitivity is governed by the perivascular boundary layer
`L = sqrt(D_eff/k)` relative to `dx` — rank/log correlation between `L` and the
ladder spread **r = −0.930** (n = 6). In the `frozen` arm `L` moves as `√D`, i.e.
by `√0.1512 = 0.389`:

| drug | `L` published (µm) | `L/dx` @ dx=12 | `L` frozen-corrected | **`L/dx` corrected** | `L` recal-corrected |
|---|---:|---:|---:|---:|---:|
| 5-fluorouracil | 209.9 | 17.5 | 81.6 | 6.80 | ~209.9 (pinned) |
| cisplatin | 142.9 | 11.9 | 55.6 | 4.63 | ~142.9 (pinned) |
| doxorubicin | 64.9 | 5.4 | 25.2 | 2.10 | 64.9 (pinned by the anchor) |
| epirubicin | 64.9 | 5.4 | 25.2 | 2.10 | 64.9 (pinned) |
| paclitaxel | 28.8 | 2.4 | 11.2 | **0.93** | ~28.8 (pinned) |
| **docetaxel** | 28.2 | 2.4 | 11.0 | **0.91** | ~28.2 (pinned) |

**In the `frozen` arm the two taxanes fall below one grid cell per boundary layer
at the published `dx = 12 µm`.** They were already the two grid-fragile drugs — a
49.7 % and 48.5 % ratio spread across the ladder at `L/dx = 2.4` — and the frozen
correction pushes them to `L/dx < 1`, where the perivascular gradient the ratio is
built from is no longer resolved at all. **Frozen-arm taxane numbers at `dx = 12`
should not be treated as converged, and this file does not treat them as such.**

In the `recal` arm `L` is pinned by the measured penetration anchor and `L/dx` is
unchanged from the published run, so the grid error is the published grid error and
nothing new is introduced.

**What this file deliberately does NOT do with the slab error law.**
`results/boundary_layer_resolution.md` carries a controlled-slab error law and a
+8.3 % figure for docetaxel. That number is **not** used here to bound or
contextualise the movement on tissue: the p99-convention investigation has since
measured that **the tissue grid dependence for docetaxel has the opposite sign to
the controlled slab**, so the slab result does not transfer. It is cited as
slab-only, and the sign on tissue is under separate investigation by another
piece of work.

### 5.5 Which verdicts move

All 12 (mode × drug) points solved on the published `dx = 12 µm` grid, 458×628,
167,780 cells / 89,785 tumour, same `dt = 10 s`, same 6 h window, same seed 0, same
`PDConfig`, same p99 bath-matching convention, same frozen literature bands.
`alpha`: published **148.6082** in `frozen`, re-solved to **22.2281** in `recal`.

| mode | drug | ratio, shipping D | ratio, corrected D | move | verdict shipping | verdict corrected | **moved?** |
|---|---|---:|---:|---:|---|---|:--:|
| frozen | doxorubicin | 0.2301 | 0.1403 | −39.0 % | PASS | PASS | no |
| frozen | paclitaxel | 0.1503 | 0.0743 | −50.6 % | PASS | PASS | no |
| **frozen** | **docetaxel** | **0.1484** | **0.0727** | **−51.0 %** | **FAIL** | **PASS** | **YES** |
| frozen | epirubicin | 0.2301 | 0.1403 | −39.0 % | PASS | PASS | no |
| frozen | 5-fluorouracil | 0.3379 | 0.2556 | −24.4 % | PASS | PASS | no |
| frozen | cisplatin | 0.3085 | 0.2162 | −29.9 % | PASS | PASS | no |
| recal | doxorubicin | 0.2301 | 0.2330 | +1.25 % | PASS | PASS | no |
| recal | paclitaxel | 0.1503 | 0.1537 | +2.25 % | PASS | PASS | no |
| **recal** | **docetaxel** | **0.1484** | **0.1517** | **+2.19 %** | **FAIL** | **FAIL** | **no** |
| recal | epirubicin | 0.2301 | 0.2330 | +1.25 % | PASS | PASS | no |
| recal | 5-fluorouracil | 0.3379 | 0.3353 | −0.77 % | PASS | PASS | no |
| recal | cisplatin | 0.3085 | 0.3071 | −0.43 % | PASS | PASS | no |

**Exactly one verdict moves: docetaxel FAIL → PASS, in `frozen` only.**

**In `recal` — the mode the codebase actually runs — nothing moves, and docetaxel
moves the *wrong way for a rescue*.** Its ratio rises 0.1484 → 0.1517 against a band
top edge of 0.1463, i.e. from 1.5 % outside its band to **3.7 % outside**. A 6.6×
correction to the diffusivity does not rescue the panel's only FAIL; it deepens it.
That is the opposite of what a rigged correction produces, and it is the single
most useful line in this section.

**Why the two arms disagree so violently, and which to believe.** This is
`alpha` behaving exactly as `results/alpha_identifiability.md` says a gauge
coordinate behaves. In `frozen`, `D` falls 6.6× while the sink is held fixed, so
`L = sqrt(D/k)` collapses by 2.6× and every penetration-shaped quantity collapses
with it — hence −24 % to −51 %. In `recal`, `alpha` absorbs the change at the
anchor (148.6082 → 22.2281, a factor 0.1496 ≈ the 0.1512 applied to `D`) and
re-exports it as a near-exact compensation on every other drug, so `L` is pinned
and the scored ratio barely moves. **`recal` is the arm to believe**, because
`calibrate_alpha()` runs at compile time and because `frozen` silently breaks the
measured Primeau penetration anchor (§5.2).

**Two caveats on the `recal` movements, and they matter more than the movements.**

1. **They are smaller than the discretisation error already sitting under them.**
   `results/grid_convergence_breast.md` measures docetaxel's error at `dx = 12` as
   **−3.4 % against its Richardson limit**. The recal move is **+2.19 %**. A
   movement smaller than the known grid error is **not cleanly attributable to this
   physics change**, and this file does not claim it is. The *verdict* is safe —
   docetaxel is FAIL at `dx = 12` (0.1484), at `dx = 8` (0.1521), at the Richardson
   limit (0.1537) and at the fitted limit (0.1573), and the correction pushes it
   further out — but the *number* is not resolved to the precision quoted.
2. **The `frozen` movements are large and clean, but their grid is not.** −39 % to
   −51 % is an order of magnitude above the grid error, so the *direction and size*
   are attributable. But §5.4 shows the frozen arm drives both taxanes to
   `L/dx < 1`, so the frozen docetaxel flip is happening on a grid that no longer
   resolves the boundary layer producing it. **The one verdict that moves is the
   one whose movement is least trustworthy** — and it moves only in the arm the
   codebase does not run.


---

## 6. What this does and does not buy

**It buys correct attribution, not a better gradient.** `results/dtissue_sensitivity.md`
§8 already established the structural point and it applies verbatim here: because
`alpha` is fitted so that `D/k` reproduces a measured penetration depth, the
measurement constrains only the *ratio*. Lowering `D` 6.6× and letting `alpha`
fall 6.6× leaves the predicted gradient almost unchanged. What changes is the
claim the model makes about **how much drug is consumed perivascularly** — by
about 6.6× — and that claim is currently unmeasured.

**It buys a falsifiable input where there was an absent term.** Before: the
hindrance law had no cellular degree of freedom and no measurement could
constrain it. After: φ_cell is the dominant uncertainty, it is derivable from
every section the engine already ingests, and the law that consumes it is a
rigorous bound rather than a fit.

**It does not fix `D_free`.** §4.2 shows a residual ~1.8× on the hindrance ratio
that the cell term cannot absorb. `results/dtissue_sensitivity.md` §6 notes the
engine's `D_free` rests on an interstitial viscosity of 1.2e-3 Pa·s (~1.7× water
at 37 °C) that the repo's own `permeability.py` does not use. That is a separate,
unresolved question and this file does not touch it.

**It does not make the hindrance drug-specific.** Maxwell's factor is identical
for cisplatin and docetaxel. The real per-drug variation — Pruijn's 69-fold spread
across one congeneric series (`external_anchor_survey.md` §2.5) — is untouched by
this fix and remains the sharpest outstanding falsification of the diffusivity
model.

**It is spatially uniform.** One φ_cell for the whole section. A per-voxel
cellularity field is derivable from the same `cell_area` column and is the natural
next step, but it lives in `tmesim/reconstruct/`, which this work does not own.

**It will move threshold statistics much harder than it moves ratios.**
`results/dose_schedule.md`, regenerated independently while this work was running,
found temozolomide becomes **delivery-limited** once its transport is corrected:
IC50 coverage at 30× collapses **97.80 % → 3.43 %** at `dx = 12 µm`. That was a
delivery-law (back-flux) correction, not a diffusivity one, but it demonstrates the
mechanism — *coverage is a fraction-of-cells-above-a-threshold statistic and is
therefore far more sensitive than the p99-normalised 3-D:2-D ratio scored here,
which cancels a uniform exposure change exactly*. A further 6.6× reduction in
`D_eff` will move coverage-type quantities again, and in the same direction. That
file is not owned here and was not re-run; the point is flagged so nobody reads
§5's small `recal` ratio movements as evidence that the correction is harmless to
every published quantity. It is not.

---

## 7. Merge recommendation

**MERGE the physics. Do NOT flip the default in this change.**

**Merge now — `tmesim/molecular/ecm_mesh.py`.** `cell_exclusion`,
`cell_volume_fraction_from_section` and `tissue_hindrance` are **purely additive**.
No existing function's signature, default or return value changes;
`tissue_hindrance(r, 0.0)` reproduces `ogston_renkin_hindrance(r)` exactly (asserted).
Nothing calls them yet, so `DrugParams.diffusivity` is bit-identical, every frozen
`DRUG_LIBRARY` snapshot is valid, and every published number stands. Test suite:
754 passing, unchanged.

**Do not wire it into `compile_drug.build_profile` in this change**, for three
reasons that are about evidence, not about caution:

1. **The correction overshoots the direct measurements (§5.3), and the overshoot
   is diagnostic.** Two pre-existing ECM legs were partly standing in for
   cell-scale obstruction. Wiring the cell term in without deciding what those
   legs are for double-counts by a known factor. That decision touches
   `brownian.py`'s obstacle scale and `ALPHA_TORTUOSITY`, and it should be made
   deliberately.
2. **The `frozen` arm breaks a measured anchor.** Any merge that changes `D` must
   ship with `alpha` re-solved. `calibrate_alpha()` does that at run time, so the
   codebase would land in `recal` — but a reader who edits the diffusivity by hand
   would silently land in `frozen` and destroy the Primeau anchor. That trap
   should be closed first.
3. **φ_cell is currently a single global constant derived from one breast section
   with a known segmentation bias**, and it sits above the literature band. It
   should be a per-tissue field before it becomes a default.

**The order of work that follows from this file:**

1. Adopt the Hicks hindrance row (12×, PMID 9806526) as a registered scorecard
   FAIL, as `external_anchor_survey.md` §4 already recommends. It now measures a
   term the code contains rather than one it lacks.
2. Resolve the double-count: decide whether `brownian.py`'s 300 nm obstacle field
   and `ALPHA_TORTUOSITY = 0.6` are matrix terms or mis-scaled cell terms. §5.3
   gives the size of the discrepancy to resolve it against.
3. Build a per-voxel φ_cell field from the `cell_area` column already present in
   every Xenium cache (`tmesim/reconstruct/` — not owned here).
4. Obtain a measured cell or extracellular volume fraction for a Wilson/Hicks
   multicellular layer (§4.1). It is the one missing input that would turn the
   6.61×-vs-12× row from a bracket into a score, and it is a literature question,
   not a compute question.
5. Only then flip the default, in `recal` mode, and regenerate the affected
   reports — including the coverage-type quantities flagged in §6, which will move
   much harder than the ratios did.

**What the blast radius says about the eventual flip.** It is reassuring but not
permission. In `recal` — the mode a flip would land in — **no verdict moves and
every scored ratio moves under 2.3 %** (§5.5), so the published breast panel is not
resting on the missing cell term. But that robustness is itself the `alpha` gauge
absorbing the change (§5.2), not evidence that the transport is right; and §6
records that threshold statistics such as IC50 coverage are not protected by the
same cancellation. A flip should therefore be judged on the physics and on
§5.3's double-count, not on the verdict table staying still.
