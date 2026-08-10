# `OxygenConfig.supply` is a transvascular mass-transfer coefficient in s⁻¹, its literature value is ~2×10² s⁻¹, and the shipped 0.05 is ~4,000× too small

> **One-line answer.** `supply` is not a tuning knob of unknown type. It is the intravascular
> oxygen mass-transfer coefficient per unit vessel volume, in **s⁻¹**, and it is measurable.
> Matching it to Secomb's measured cerebral `K_TO₂,cap` gives **supply ≈ 194–251 s⁻¹** on the
> GBM section — the shipped `0.05` corresponds to a `K_TO₂,cap` of **0.00027** against a measured
> **1.213**, i.e. **4,570× too small**. The signature of the error is visible without any
> literature at all: under the shipped value **98.1 % of the entire blood-to-tissue oxygen drop
> happens across the vessel wall** and only 1.9 % across the tissue, which is the opposite of
> the measured ~45/55 split and the opposite of the Krogh picture the module's docstring cites.
>
> **The O2 variable is dimensional, not a convention** — but the *derivation of `supply` does not
> depend on that*, which is why this constant can be calibrated when `consumption` could not.
>
> **The uncomfortable part, stated plainly and reported separately in §6:** correcting `supply`
> ALONE — holding the uncalibrated `consumption` at 0.02 — drives mean tumour O2 to 0.9925 and
> `F_state` to exactly **1.0000**, which makes `results/kill_calibration.md`'s headline
> **2/3, with temozolomide failing**. That is the same 2/3 `results/oxygen_level_shift.md`
> measured under the historical KDE field, reached by a completely different route, and the
> predicted temozolomide ratio agrees with it **to sixteen significant figures**.
>
> **But `supply` cannot honestly be moved alone.** The same conversion machinery shows the
> model's oxygen SINK is ~140× too small on this section. Pinning the source with Secomb and
> the sink with Watabe's PET-measured glioma TMRO₂ produces, for the first time in this
> project, a GBM oxygen field whose **median tumour pO₂ (4.5–8.2 mmHg) brackets the measured
> human-GBM median (Lally 2006: 5.1 mmHg)** — and under *that* jointly calibrated
> configuration the ratio test is **3/3 PASS again**. The shipped configuration is two large
> compensating errors, and which verdict you get depends on whether you fix one or both.

**Status: derivation and recommendation. Nothing was changed.** `OxygenConfig.supply` is
still `0.05`, `consumption` is still `0.02`, and `tmesim/` is untouched by this work. The
only intervention anywhere is a harness-level wrapper around `physics.oxygen.solve_oxygen`
that hands it `dataclasses.replace(cfg, supply=S)` — the solve is unmodified library code and
`supply` is a public field of a public dataclass.

* Generator: `experiments/oxygen_supply_calibration.py`
* Raw: `results/sweep/oxsupply_algebra.json`, `oxsupply_conductance.json`, `oxsupply_ladder.json`,
  `oxsupply_joint.json`, `oxsupply_sinkmatch.json`, `oxsupply_sections.json`,
  `oxsupply_grid.json`, `oxsupply_killcal.json`
* Sections: `data_cache/xenium_gbm` (40,887 cells, 23,812 tumour, grid 424×675) and
  `data_cache/xenium_breast_janesick/rep1`, both at dx = 12 µm.

---

## 1. What `supply` IS, dimensionally — established before looking at any number

`tmesim/physics/oxygen.py` solves

```
D * lap(O) + supply * vasc * (1 - O) - consumption * rho_tum * O = 0
```

`D` is µm²/s and `lap` divides by µm², so `D*lap(O)` carries units of `O`/s. `vasc` is
dimensionless (clipped to [0,1]; under the shipped `mode='discrete'` it is **exactly 1.0
inside a vessel pixel and exactly 0 outside**, `reconstruct/vessel_network.py:539`), and
`(1-O)` is dimensionless. So by balance:

> **`supply` has units of s⁻¹. It is a first-order rate constant.**

Which one? The term `supply * vasc * (1-O)` is a linear transvascular exchange term — the
standard `(P·S/V)·(C_blood − C_tissue)` form. `P` is the wall permeability (µm/s) and `S/V`
the exchanging surface per unit volume (µm⁻¹), so their product is s⁻¹. Evaluating `S/V`
inside the vessel (surface per unit *vessel* volume, `2/R` for a cylinder):

> **`supply` = P_eff · (2/R) — the effective oxygen conductance from blood into tissue,
> per unit vessel volume, including whatever intraluminal resistance the blood itself
> contributes.**

Two structural consequences follow immediately, before any measurement:

**(a) `O = 1` is the BLOOD compartment, by construction.** The source vanishes only at
`O = 1`; a region of pure source with no sink relaxes to `O = 1` exactly. So whatever `O = 1`
means physically, it is the tension the vessel is equilibrating tissue *towards* — i.e.
intravascular pO₂. This matters in §2.

**(b) The classical Krogh formulation this module cites is the `supply → ∞` limit of what it
actually solves.** Krogh, and Grimes et al. 2014 (PMID 26064525) which the docstring cites for
the sink form, impose `O = 1` as a *Dirichlet boundary condition at the vessel wall*. The
distributed-source form here recovers that only when `supply` is large enough to pin the
vessel interior. The relevant comparison rates, all from the shipped constants
(`results/sweep/oxsupply_algebra.json`):

| rate | value (s⁻¹) | what it is |
|---|---:|---|
| `supply` (shipped) | **0.05** | the constant under investigation |
| `D / R²`, R = 7.5 µm | **35.56** | rate at which diffusion drains a vessel-sized region |
| `D / dx²`, dx = 12 µm | 13.89 | rate at which diffusion drains one grid cell |
| `consumption · ρ_tum`, ρ = 1 … 8 | 0.02 … 0.16 | the sink |

The source equilibration length `sqrt(D/supply)` is **200 µm** at the shipped value, against a
vessel radius of **7.5 µm**. The vessel is 27× smaller than the distance over which its own
source acts, so it cannot pin its own oxygen: **the model's vessels are not vessels, they are
weak volumetric sources sitting in a bath.**

### 1.1 The Krogh cross-check, before any literature value

`stage_algebra` builds a self-contained axisymmetric Krogh cylinder (source for r < 7.5 µm,
uniform sink beyond, no-flux at the tissue boundary) and compares the distributed-source
solution against the Dirichlet solution for the same geometry and sink:

| `supply` (s⁻¹) | ρ_tum = 1: O2 at wall | mean tissue O2 | **% of Dirichlet** | ρ_tum = 8: mean tissue O2 | **% of Dirichlet** |
|---:|---:|---:|---:|---:|---:|
| **0.05 (shipped)** | 0.0391 | 0.0381 | **3.9 %** | 0.0049 | **0.6 %** |
| 0.5 | 0.2887 | 0.2818 | 28.9 % | 0.0468 | 5.6 % |
| 5 | 0.7998 | 0.7807 | 80.0 % | 0.3081 | 36.8 % |
| 35.56 (= D/R²) | 0.9627 | 0.9396 | 96.3 % | 0.6612 | 79.0 % |
| 100 | 0.9841 | 0.9605 | 98.4 % | 0.7532 | 90.0 % |
| 1000 | 0.9964 | 0.9725 | 99.6 % | 0.8166 | 97.6 % |
| Dirichlet (O = 1 at wall) | 1.0000 | 0.9760 | 100 % | 0.8368 | 100 % |

**At the shipped value the model delivers 0.6–3.9 % of the oxygen the same geometry and the
same sink would deliver under the boundary condition its own documentation describes.** No
literature is needed for that statement.

---

## 2. Is the O2 variable dimensional, or a normalisation convention? — Both, and it does not matter for this derivation

This is the question the task flagged as changing what a correct answer looks like, so it is
answered explicitly rather than assumed.

**`O` is a normalised variable with a declared physical anchor.** `oxygen.py` defines
`REFERENCE_PO2_MMHG = 35.0` and `to_mmhg()` as a pure rescale, so `O` is `p / 35 mmHg` — a
dimensional pO₂ in disguise. The band edges at 0.2 / 0.4 are therefore 7 / 14 mmHg, which is a
*convention about where to draw a biological line*, not a physical constant.

**The anchor is in the wrong compartment, and that is a real (small) defect.** §1(a) showed
`O = 1` is the tension the vessel equilibrates towards — intravascular pO₂. `REFERENCE_PO2_MMHG`
is anchored to normal brain **tissue** pO₂ (Erecinska & Silver 2001, PMID 11718758), and the
module's own docstring argues for the tissue anchor on the grounds that "the PDE's `vasc` field
is a vascular-density source term feeding a tissue compartment, not a blood compartment". That
reasoning describes what the source *feeds*, not what the source *saturates at*. Secomb's
reference cortical values are `P̄_cap = 58.75 mmHg` and `P̄_tiss = 47.65 mmHg` — so a
compartment-correct anchor for `O = 1` would be ~55–60 mmHg, not 35, and every mmHg figure in
this report would scale up by ~1.7×. **This is flagged, not adopted:** re-anchoring moves the
band edges to 11.75 / 23.5 mmHg and is a separate decision from the one this report is about.

**Crucially, the derivation in §3 is independent of the anchor.** Writing the model's own
intravascular mass transfer coefficient out (§3.1) the reference pressure cancels exactly.
That is the structural reason `supply` is calibratable when `consumption` was not: the
project's two previous passes at `consumption`
(`results/oxygen_calibration_landscape.md` §4, `results/sweep/oxygen_ocr_literature_second_pass_*.md`)
both stalled on "the model has never defined what `O = 1` means in mmHg", and both were right
to — but that blocker does not apply to a *conductance*, which is a ratio of a flux to a
pressure difference and is blind to the scale on which both are expressed.

**So: this is situation (1), not situation (2). `supply` is a genuine physical constant with
a genuine measurable value, and it can be calibrated from public data.**

---

## 3. The literature value, with provenance

### 3.1 The closed form

Averaged over tissue, `solve_oxygen`'s source per unit tissue volume is
`supply · f_v · (1 − O_vessel)` in O/s, where `f_v` is the vessel-pixel **area** fraction of
tissue (which equals the volume fraction for an isotropic 3-D medium, by Delesse's principle).
Converting to a physical oxygen flux with `c = α·p` and `p = P_ref·O`:

```
M      = supply * f_v * alpha * P_ref * (1 - O_vessel)  =  supply * f_v * alpha * (P_ref - p_vessel)
dP_cap = P_ref - p_vessel
K_cap(model) = M / dP_cap  =  supply * f_v * alpha
```

The `(1 − O_vessel)` and `P_ref` both cancel. **`K_cap(model) = supply · f_v · α`** —
independent of `consumption`, of `diffusivity`, of the solved field, and of the mmHg anchor.
Inverting:

```
supply = K_TO2,cap / (alpha * f_v)
```

### 3.2 The measured coefficient — primary source, verified two ways

**Secomb TW. "The mass transfer coefficient for oxygen transport from blood to tissue in
cerebral cortex." *J Cereb Blood Flow Metab* 2020;40(8):1634–1646. PMID
[31423930](https://pubmed.ncbi.nlm.nih.gov/31423930/), PMC
[7370375](https://pmc.ncbi.nlm.nih.gov/articles/PMC7370375/), DOI
[10.1177/0271678X19870068](https://doi.org/10.1177/0271678X19870068).** (Retrieved via PubMed.)

Quoted verbatim from the full text: *"the predicted values are KTO2tiss = 0.997, KTO2cap =
1.213 and KTO2 = 0.547, all in units of cm3O2 · (100 cm3 · min · mmHg)−1"*, with
`KTO2 = JO2/(P̄cap − P̄tiss)` and `1/KTO2 = 1/KTO2tiss + 1/KTO2cap`, and steady-state
*"P̄cap = 58.75 mmHg, P̄wall = 53.72 mmHg and P̄tiss = 47.65 mmHg"*. Reference parameter set:
`rc = 2 µm`, `rt = 25 µm`, `Ncap = 509.3/mm²`, `HT = 0.281`, `D = 2410 µm²/s`,
`α = 3.89×10⁻⁵ cm³O₂·(cm³·mmHg)⁻¹`.

This coefficient set is **verified, not trusted**, two independent ways:

1. **Internal consistency.** `1/0.997 + 1/1.213 = 1.827 → 0.547` ✓, and
   `6.0 / (58.75 − 47.65) = 0.541 ≈ 0.547` ✓ against the stated reference CMRO₂.
2. **Reproduced from the paper's own analytic formula.** `stage_conductance` implements
   `K_tiss = 8πDαN_cap(1−δ²)²/(−4 lnδ − 3 + 4δ² − δ⁴)`, `δ = rc/rt`, evaluates it at the
   paper's own reference geometry, and **asserts** it lands on the reported 0.997. It returns
   **0.99714 — a relative error of 1.4×10⁻⁴.** The run fails loudly if it does not.

Why this source and not a wall-permeability figure: the paper states directly that the wall is
*not* the resistance — *"This is not the case for oxygen, which passes readily through cell
membranes, being a non-polar, lipid-soluble molecule… radial gradients of PO2 are present
within the microvessels and the surrounding tissue, and are not concentrated at the vessel
wall."* `K_TO₂,cap` is the right object because it is the *lumped* intravascular conductance
(red-cell spacing, plasma gap, wall) — the same lumped thing `supply` is. That intraluminal
resistance is real and non-negligible is the long-established Hellums result: Hellums JD,
*Microvasc Res* 1977;13(1):131–6, PMID [859450](https://pubmed.ncbi.nlm.nih.gov/859450/),
DOI [10.1016/0026-2862(77)90122-4](https://doi.org/10.1016/0026-2862(77)90122-4)
("the resistance in the capillary is a significant fraction of the overall resistance");
Baxley & Hellums, *Ann Biomed Eng* 1983;11(5):401–16, PMID
[6679691](https://pubmed.ncbi.nlm.nih.gov/6679691/), DOI
[10.1007/BF02584216](https://doi.org/10.1007/BF02584216); Hellums, Nair, Huang & Ohshima,
*Ann Biomed Eng* 1996;24(1):1–24, PMID [8669708](https://pubmed.ncbi.nlm.nih.gov/8669708/),
DOI [10.1007/BF02770991](https://doi.org/10.1007/BF02770991) (*"Intraluminal resistance to gas
transport … was neglected for a half-century following the early work of Krogh. In recent
years it has come to be understood that this neglect is seriously in error."*). Abstracts
retrieved and verified; the 1996 review's own tabulated mass-transfer coefficients were **not**
extracted (Springer full text is paywalled) and are flagged **UNVERIFIED** — they are not used
for any number here, only for the qualitative claim that intraluminal resistance is
non-negligible, which the abstracts state directly.

**O₂ solubility α.** Three independently-sourced values are carried as a band rather than one
picked: `3.0×10⁻⁵` (the classical Bunsen figure, 0.003 mL O₂/100 mL/mmHg for plasma at 37 °C),
`3.3×10⁻⁵` (already derived in-repo from Grimes et al. 2014's Ω constant,
`results/sweep/oxygen_ocr_literature_second_pass_*.md` Step 1), and `3.89×10⁻⁵` (Secomb's own).
All in cm³O₂·(cm³·mmHg)⁻¹. They span 30 %.

### 3.3 The number

| | GBM (`xenium_gbm`) | breast (`janesick/rep1`) |
|---|---:|---:|
| vessel-pixel area fraction of tissue `f_v` | 0.02683 | 0.05672 |
| **`supply` implied, α = 3.0 / 3.3 / 3.89 ×10⁻⁵** | **251.1 / 228.3 / 193.7** | **118.8 / 108.0 / 91.6** |
| cross-check: Secomb's own geometry (f_v = 0.0064) rescaled by radius (2 µm → 7.5 µm) | 280.8 / 255.3 / 216.5 | (same, section-independent) |
| **`K_cap` the SHIPPED `supply = 0.05` corresponds to** | **0.00024 – 0.00031** | 0.00051 – 0.00066 |
| ratio to Secomb's measured 1.213 | **1/4,570** | 1/2,160 |

> **Recommended value: `supply ≈ 2×10² s⁻¹`, with a defensible band of roughly
> 10² – 10³ s⁻¹.** The band's width comes from α (±30 %), from `K_TO₂,cap`'s own stated
> dependence on capillary density and tube hematocrit (tumour vessels have both lower and more
> heterogeneous — call it a factor of ~2 each way), and from the section-dependence below.

**Two honest caveats on the number, both of which matter:**

1. **`supply` is not section-invariant under this parameterisation.** `f_v` is a property of
   the section *and of the vessel-painting convention*, so the closed form gives 194 s⁻¹ on GBM
   and 92 s⁻¹ on breast — a 2.1× spread that is exactly the ratio of their painted vessel area
   fractions. The physically invariant quantity is the *conductance per unit tissue volume*
   `G = supply·f_v = K_TO₂,cap/α = 5.2–6.7 s⁻¹`. The radius-scaling route (which is
   section-independent by construction: 217–281 s⁻¹) agrees with the GBM closed form to 12 %
   and disagrees with the breast one by 2.4×, which locates the discrepancy in the breast
   section's painted vessel fraction rather than in the derivation.
2. **The 7.5 µm painting radius is already known to be wrong.** `results/vascular_segmentation_audit.md`
   establishes the production `vascular_density` is a union of 7.5 µm discs around endothelial
   *cells* (what `vessel_geometry.py` calls "Wrong #2"), giving 2.683 % coverage against the
   2.1 ± 0.3 % the module's own `PUBLISHED` block cites. `supply` derived through `f_v`
   inherits that error linearly. **Fixing the vessel geometry and fixing `supply` are the same
   task and should be done together.**

### 3.4 The signature, visible without any literature

`stage_ladder` records the two halves of the blood-to-tissue drop separately. Under the
shipped value the model puts essentially all of it across the wall:

| `supply` | O2 at vessel pixels | mean tumour O2 | **transvascular share of total drop** |
|---:|---:|---:|---:|
| **0.05 (shipped)** | 0.3613 | 0.3486 | **98.1 %** |
| 5 | 0.9822 | 0.9715 | 62.4 % |
| 35.56 | 0.9975 | 0.9897 | 24.8 % |
| **194 (derived)** | **0.9995** | **0.9925** | **6.3 %** |
| Dirichlet limit | 1.0000 | 0.9932 | 0.0 % |
| *Secomb, measured cortex* | *(P̄_wall/P̄_cap = 0.914)* | | **45.3 %** *(5.03 of 11.1 mmHg)* |

**The shipped configuration says the vessel wall is 50× more resistive than the tissue. The
measurement says they are comparable.** The derived value overshoots in the other direction
(6.3 % rather than 45 %) — because the model's *tissue* resistance is also too low, which is
§5's subject, not this one's.

---

## 4. What that value produces on the two real sections

`results/sweep/oxsupply_sections.json`, dx = 12 µm, everything else shipped:

| section | arm | O2 at vessels | **mean tumour O2** | mmHg (35-anchor) | median mmHg | H / Q / V | frac < 10 mmHg |
|---|---|---:|---:|---:|---:|---|---:|
| GBM | shipped 0.05 | 0.3613 | **0.3486** | 12.20 | 12.11 | 0.000 / 0.900 / 0.100 | 0.001 |
| GBM | **derived 194** | 0.9995 | **0.9925** | **34.74** | 34.76 | 0.000 / 0.000 / 1.000 | 0.000 |
| GBM | Dirichlet 10⁵ | 1.0000 | 0.9932 | 34.76 | 34.78 | 0.000 / 0.000 / 1.000 | 0.000 |
| breast | shipped 0.05 | 0.3345 | **0.2815** | 9.85 | 9.76 | 0.112 / 0.842 / 0.046 | 0.548 |
| breast | **derived 194** | 0.9994 | **0.9747** | **34.12** | 34.51 | 0.000 / 0.000 / 1.000 | 0.000 |
| breast | Dirichlet 10⁵ | 1.0000 | 0.9764 | 34.18 | 34.56 | 0.000 / 0.000 / 1.000 | 0.000 |

### 4.1 Against measured tumour pO₂ — the derived value ALONE makes the field worse

Human GBM pO₂, from the four studies already assembled in
`results/oxygen_calibration_landscape.md` §1 (not re-verified here; that report's own
UNVERIFIED flags carry over):

* Lally et al. 2006, PMID [17207315](https://pubmed.ncbi.nlm.nih.gov/17207315/) — **median
  tumour pO₂ = 5.1 mmHg**, n = 23 gliomas, Eppendorf histograph, intraoperative.
* Collingridge et al. 1999, PMID [10665789](https://pubmed.ncbi.nlm.nih.gov/10665789/) — both
  glioma and peritumoural brain contain regions **< 2.5 mmHg**.
* Beppu et al. 2002, PMID [12160140](https://pubmed.ncbi.nlm.nih.gov/12160140/) —
  intratumoral GBM ≈ **15.3 mmHg** (UNVERIFIED to the decimal in the source report).
* Nozue et al. 1997, PMID [9290690](https://pubmed.ncbi.nlm.nih.gov/9290690/) — interlab
  spread on an identical protocol was 1.5–5.6 mmHg, so treat any single value as a range.

**Verdict: `supply = 194` alone gives mean tumour pO₂ = 34.7 mmHg and a hypoxic fraction of
exactly zero, against a measured median of ~5 mmHg. That is worse than the shipped value's
12.2 mmHg, not better.** This is stated first and plainly because it is the finding that most
constrains what to do next: **correcting `supply` in isolation would be indefensible**, not
because of what it costs downstream, but because the resulting oxygen field is further from
the measurement than the one it replaces.

The reason is not that the derivation is wrong. It is that **the shipped `supply = 0.05` is a
compensating error**: it suppresses mean tumour O2 by 2.85× (0.9925 → 0.3486) through a
physically wrong mechanism — a fictitious wall resistance carrying 98 % of the drop — and
thereby lands the absolute pO₂ in a range that *looks* plausible. Remove the compensating
error without removing what it compensates for and the residual error (a sink ~140× too
small, §5) is exposed at full size.

**Breast**: no primary numeric median breast-carcinoma pO₂ was retrievable in this pass.
PubMed returned Vaupel's reviews (e.g. PMID
[9879340](https://pubmed.ncbi.nlm.nih.gov/9879340/),
[8798832](https://pubmed.ncbi.nlm.nih.gov/8798832/)) but no median value appears in the
retrievable abstracts and the full texts are not in PMC. **The breast pO₂ comparison is
therefore UNVERIFIED** and the breast numbers above are reported without a measured target.

---

## 5. What the same machinery says about the sink — because the two cannot be separated

The conversion in §3.1 also runs in the other direction, on the tissue side, and this is
where §4's residual goes.

**The model's oxygen sink is ~140× too small on this section.** At the Dirichlet limit (source
no longer the limiting factor), the physical consumption flux the shipped `consumption = 0.02`
produces is:

| | GBM | breast |
|---|---:|---:|
| tissue-mean `rho_tum` (cells per 144 µm² tile) | **0.1229** | 0.3352 |
| implied CMRO₂, cm³O₂·(100 cm³·min)⁻¹ | **0.0169** | 0.0454 |
| Watabe 2017 C6-glioma TMRO₂ (measured, ¹⁵O-PET) | 1.9 – 2.8 | — |
| Secomb cortical reference CMRO₂ | 6.0 | — |
| **shortfall vs Watabe** | **112 – 166×** | — |

Watabe T et al., *EJNMMI Res* 2017;7:16, PMID
[28210996](https://pubmed.ncbi.nlm.nih.gov/28210996/), DOI
[10.1186/s13550-017-0263-6](https://doi.org/10.1186/s13550-017-0263-6) — already assembled and
verified in `results/sweep/oxygen_ocr_literature_second_pass_*.md` §2b.

**This sharpens that report's own conclusion rather than contradicting it.** It converted the
same Watabe number and concluded `consumption` should be 12–25× higher (0.240–0.505 s⁻¹). That
conversion evaluated the sink at an *illustrative* `rho_tum` of 1–8, taken from the module's
synthetic self-test. **On the real section the tissue-mean `rho_tum` is 0.1229** — 8–65× lower
than the illustrative values — and the remaining factor is exactly that. The 12–25× and the
112–166× are the same calculation at two different densities; the second uses the section's
own measured density and is the one that applies here.

### 5.1 The jointly calibrated configuration — and it lands on the measured GBM median

`stage_sinkmatch` pins the source at the derived `supply = 194` and sweeps the sink
(`results/sweep/oxsupply_sinkmatch.json`):

| `consumption` | implied CMRO₂ | mean tumour pO₂ | **median** | p10 | p90 | HYPOXIC frac | frac < 10 mmHg | `F_state` |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.02 (shipped) | 0.017 | 34.74 | 34.76 | 34.56 | 34.89 | 0.0000 | 0.0000 | 1.0000 |
| 0.5 | 0.361 | 29.64 | 29.94 | 26.34 | 32.48 | 0.0000 | 0.0000 | 1.0000 |
| 2.5 | 1.171 | 19.25 | 19.21 | 12.43 | 25.96 | 0.0037 | 0.037 | 0.9177 |
| 5 | 1.684 | 13.84 | 13.32 | 6.89 | 21.31 | 0.1048 | 0.269 | 0.7108 |
| **10** | **2.227** ✔ | **9.15** | **8.19** | 3.22 | 16.23 | 0.4022 | 0.627 | **0.5249** |
| **20** | **2.748** ✔ | **5.65** | **4.53** | 1.30 | 11.34 | 0.7076 | 0.858 | **0.4167** |
| 50 | 3.361 | 2.76 | 1.82 | 0.34 | 6.20 | 0.9273 | 0.977 | 0.3655 |

✔ = implied CMRO₂ inside Watabe's measured 1.9–2.8 band.

> **The two rungs that reproduce the measured glioma oxygen consumption rate produce a
> median tumour pO₂ of 4.53 and 8.19 mmHg. Lally et al.'s measured human median is 5.1 mmHg.**

Nothing was fitted to that median — the source was pinned by Secomb's `K_TO₂,cap`, the sink by
Watabe's TMRO₂, and the pO₂ distribution fell out. It is the first time in this project that
the oxygen field has been anchored on both sides by independent measurements and landed on an
independent third measurement. It also, for the first time, predicts a **non-zero, non-absurd
hypoxic fraction (0.40–0.71)** for a tumour type defined partly by pseudopalisading necrosis —
resolving from the correct side the "0.0 is no more credible than 0.94" caveat the module
docstring has carried since the grid fixes.

**Do not over-read it.** Three real weaknesses: (i) `consumption ≈ 10–20 s⁻¹` is 500–1000×
the shipped default and 20–80× above even the higher of the two brackets the previous pass
derived, so it is a large extrapolation of that report's own band; (ii) the whole comparison
is on the 35 mmHg anchor, and §2 argues the compartment-correct anchor is ~1.7× higher, which
would shift the matched rung; (iii) the sink is still first-order in O, which
`results/oxygen_sensitivity_band.md` and Grimes et al. 2014 both establish is the wrong
functional form, and a Michaelis–Menten sink would reach the same CMRO₂ at a different
constant. **The claim here is that the source and sink are jointly identifiable and that the
jointly-calibrated point is consistent with three independent measurements — not that
`consumption = 10` is now calibrated.**

---

# 6. Consequences — reported separately, and deliberately last

Everything above is a statement about transport coefficients and measured oxygen tensions.
Nothing in it was chosen with reference to anything below. Sections 1–5 should be readable,
and disagreed with, without this section existing.

### 6.1 `F_state`

| configuration | mean tumour O2 | H / Q / V | `F_state` |
|---|---:|---|---:|
| shipped (`supply` 0.05, `consumption` 0.02) | 0.3486 | 0.000 / 0.900 / 0.100 | **0.5502** |
| historical KDE field (`results/oxygen_level_shift.md`) | 0.5945 | 0.000 / 0.000 / 1.000 | **1.0000** |
| **derived `supply` only (194, `consumption` 0.02)** | **0.9925** | 0.000 / 0.000 / 1.000 | **1.0000** |
| **jointly calibrated (194, 10)** | 0.2614 | 0.402 / … | **0.5249** |
| **jointly calibrated (194, 20)** | 0.1614 | 0.708 / … | **0.4167** |

Correcting the source alone reproduces the historical KDE's `F_state` **exactly** — 1.0000, every
tumour cell VIABLE. Correcting both lands at 0.42–0.52, remarkably close to the shipped 0.5502,
**but for a completely different reason**: the shipped 0.5502 is 90 % QUIESCENT and 0 % HYPOXIC,
whereas the jointly calibrated 0.42–0.52 is 40–71 % HYPOXIC. The scalar coincides; the biology
underneath it does not.

### 6.2 The GBM 3/3-vs-2/3 verdict — measured by re-running the unmodified generator

`validation.kill_calibration.run()` + `analyse()`, unmodified, `results/sweep/oxsupply_killcal.json`:

| arm | `supply` | `consumption` | **temozolomide** | lit band | verdict | cisplatin | doxorubicin | **count** |
|---|---:|---:|---:|---|---|---|---|---|
| published (shipped) | 0.05 | 0.02 | 0.2626 | 0.159 – 0.462 | PASS | PASS | PASS | **3/3** |
| historical KDE (`oxygen_level_shift.md`) | 0.05 | 0.02 | 0.46784 | 0.159 – 0.462 | **FAIL** | PASS | PASS | **2/3** |
| **derived `supply` only** | **194** | 0.02 | **0.46784** | 0.159 – 0.462 | **FAIL** | PASS (0.498) | PASS (0.369) | **2/3** |
| **jointly calibrated, CMRO₂ low** | **194** | **10** | **0.33218** | 0.159 – 0.462 | **PASS** | PASS (0.277) | PASS (0.221) | **3/3** |
| **jointly calibrated, CMRO₂ high** | **194** | **20** | **0.34240** | 0.159 – 0.462 | **PASS** | PASS (0.216) | PASS (0.170) | **3/3** |

**Three things worth saying about this table, in order of how much they matter.**

1. **Fixing `supply` alone costs the flagship result, and that is the honest finding.**
   2/3, temozolomide — the drug `kill_calibration.md` calls its matched anchor — failing.
   The derived-supply arm's predicted ratio is `0.46783514742473475`, which is the
   `oxygen_level_shift.md` KDE arm's value **to all sixteen digits**, because both routes drive
   `F_state` to exactly 1.0000 and the ratio is linear in it. Two unrelated interventions
   reaching a bit-identical number is a strong cross-check that this is one mechanism, not two.
2. **Fixing both restores 3/3 — but this must not be read as vindication.** It is a *different
   model* passing: the tumour is 40–71 % hypoxic instead of 90 % quiescent, absolute kills move,
   and `consumption = 10–20` is a much larger extrapolation than anything this project has
   adopted. It also gets *closer to its band edge on cisplatin*: at `consumption = 20` the
   predicted 0.21644 sits **3.7 % above** the band floor of 0.20864. A modestly different
   `consumption` would fail cisplatin instead of temozolomide.
3. **Under every configuration examined, the verdict is a function of an uncalibrated or
   newly-recalibrated oxygen constant.** The published 3/3 and the derived 2/3 and the joint
   3/3 are three readings of the same lever. `results/kill_calibration.md`'s framing that
   "nothing in that prediction was tuned" is accurate about the *drug* parameters and
   incomplete about the *oxygen* ones.

### 6.3 The grid pathology — the derived value eliminates it

`results/grid_convergence_gbm.md` §5 and `results/tail90_falsification.md` §6 both route their
non-monotone-`F_state` problem to the oxygen workstream. Re-running that ladder
(`results/sweep/oxsupply_grid.json`):

| dx (µm) | shipped: mean O2 | **shipped `F_state`** | derived: mean O2 | **derived `F_state`** |
|---:|---:|---:|---:|---:|
| 30 | 0.3364 | 0.5041 | 0.9728 | **1.0000** |
| 20 | 0.3446 | 0.5402 | 0.9862 | **1.0000** |
| 16 | 0.3478 | 0.5333 | 0.9903 | **1.0000** |
| 12 | 0.3486 | **0.5502** | 0.9925 | **1.0000** |
| 10 | 0.3461 | 0.5305 | 0.9930 | **1.0000** |
| 8 | 0.3482 | 0.5488 | 0.9933 | **1.0000** |
| **span** | | **0.0461, non-monotone** | | **0.0000** |

**Span exactly zero, at every rung.** This confirms the mechanism `oxygen_level_shift.md` §6
proposed and gives it its cause: `F_state`'s grid dependence is not a property of the section,
the mesh or the solver — it is what happens when an *under-supplied* field parks mean tumour O2
0.05 below an absolute band edge at 0.4, so that a few-percent grid perturbation reshuffles
thousands of cells across it. Restore the supply and the field moves clear of every band edge.

**Caveat, and it is not a small one:** the derived-supply arm achieves span zero by putting
100 % of tumour cells in one band, which is the same degenerate way the historical KDE achieved
it. Zero span here means "nothing can cross an edge because nothing is near one", not
"the discretisation converged". The jointly calibrated configuration, which sits at
`F_state` 0.42–0.52 with a genuine spread of cells across bands, was **not** run across the dx
ladder and its grid behaviour is **unknown**. That is the single most important untested
question this report leaves open.

---

## 7. Recommendation

**Do not change `supply` on its own.** §4 is decisive: it makes the oxygen field further from
measured GBM pO₂, and §6.2 shows it costs the flagship verdict — in that order of importance.

**Do treat `supply = 0.05` as established-wrong rather than merely uncalibrated.** Unlike
`consumption`, this constant has a closed-form physical identity, a measured target, and a
value that is wrong by ~3.7 orders of magnitude with an internal signature (§3.4: 98 % of the
oxygen drop across the vessel wall) that needs no external data to recognise.

**The defensible next step is a joint re-derivation, in this order:**

1. Fix the vessel geometry first (`results/vessel_geometry_fix.md`'s `geometry='subgrid'` path,
   already implemented and opt-in). `supply` derived through `f_v` inherits the 7.5 µm painting
   error linearly, so calibrating `supply` against a known-wrong `f_v` bakes that error in.
2. Re-derive `supply` from `K_TO₂,cap / (α · f_v)` on the corrected field.
3. Re-derive `consumption` against the *section's own measured* `rho_tum` (0.1229 on GBM),
   not an illustrative value — §5 shows this changes the answer by 8–65×.
4. Settle the sink's functional form (`results/oxygen_sensitivity_band.md`) before, not after,
   fixing the constant, since the two are coupled.
5. Re-anchor `REFERENCE_PO2_MMHG` to the capillary compartment (~55–60 mmHg) or document
   explicitly why the tissue anchor is kept on a variable whose `O = 1` is a blood tension (§2).
6. Only then re-run the dx ladder and the ratio test, and report the verdict that comes out.

**Every published `F_state`, absolute kill and kill-derived ratio on the GBM section should be
read as conditional on an oxygen source constant that is ~4,000× from its measured value.**
The claim classes `results/oxygen_level_shift.md` §5 lists are exactly the affected set; this
report changes the *size* of that exposure from "unidentified" to "identified and large", and
adds `results/oxygen_sensitivity_band.md` to it (that study swept the sink over 2,500× while
holding the source at a value now shown to be off by 4,000×).

---

## 8. Limitations

1. **One coefficient, one species, one tissue.** `K_TO₂,cap = 1.213` is mouse cerebral cortex.
   Tumour microvasculature differs in hematocrit, calibre and flow heterogeneity, all of which
   Secomb states the coefficient depends on. The ±2× band in §3.3 is an estimate of that, not a
   measurement of it.
2. **The Hellums 1996 tabulated coefficients are UNVERIFIED** (Springer paywall). They are used
   only for the qualitative claim that intraluminal resistance is non-negligible, which the
   retrieved abstracts state directly.
3. **Breast has no measured pO₂ target** (§4.1) — the comparison there is model-to-model only.
4. **`consumption = 10–20` is this report's extrapolation**, not an adopted value, and is a
   large step beyond the 0.240–0.505 the previous pass derived. §5 explains the difference
   (illustrative vs. measured `rho_tum`) but the previous pass's band was not re-derived here.
5. **The jointly calibrated configuration was not grid-converged** (§6.3 caveat) and was run on
   GBM only.
6. **The sink is still first-order in O throughout.** Every number here inherits that
   functional-form question, which is not this report's to settle.
7. **`F_state` here is computed exactly as `oxygen_level_shift.md` computes it**
   (`PDConfig` multipliers over `update_cell_states` bands), reproducing its published 0.5502
   at dx = 12 and its full shipped dx ladder to four decimals — that agreement is the
   cross-check that this is the same quantity.
8. **No code under `tmesim/` was modified, so `pytest` was not re-run.** The stated baseline of
   754 passing is unaffected by construction; the only files this work creates are
   `experiments/oxygen_supply_calibration.py`, this report, and `results/sweep/oxsupply_*.json`.

---

## 9. Reproduce

```
./.venv/bin/python experiments/oxygen_supply_calibration.py algebra conductance ladder \
                   joint sinkmatch sections grid killcal
```

Each stage checkpoints to `results/sweep/oxsupply_<stage>.json`. Total wall clock on this
machine is about 7 minutes, dominated by `killcal` (3 × ~48 s) and `grid`. `conductance`
**asserts** that the Krogh formula reproduces Secomb's own reported `K_TO₂,tiss = 0.997` and
aborts if it does not, so a silently-wrong coefficient set cannot propagate into the numbers.
