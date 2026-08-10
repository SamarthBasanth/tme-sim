# Interstitial fluid pressure: the engine's only absolute pressure quantity, scored

**Verdict in one line.** The engine's solved IFP is **NOT TESTABLE against the
human intracranial-tumour anchor** — not because the quantities fail to line up,
but because that anchor cannot discriminate a working model from an empty one —
and it **FAILS the equivalent breast anchor by ~40x**, low. The landing page
should stop displaying a numeric mmHg range.

Driver: `experiments/ifp_validation.py` · raw records:
`results/sweep/ifp_validation.json` · registered as scorecard group
`ifp-wick-in-needle` (3 rows: NOT TESTABLE, **FAIL**, NOT TESTABLE).

Nothing in `tmesim/physics/pressure.py` was read back, fitted or adjusted. The
module is imported and called exactly as it ships.

---

## 1. The anchor, verified against the primary source

According to PubMed, retrieved by PMID and read in full text from PMC:

> Boucher Y, Salehi H, Witwer B, Harsh GR, Jain RK. "Interstitial fluid pressure
> in intracranial tumours in patients and in rodents." *Br J Cancer*
> 1997;75(6):829-836. **PMID 9062403**, PMC2063404,
> [DOI](https://doi.org/10.1038/bjc.1997.148)

**Results, "Measurements in patients", verbatim:**

> "The IFP was measured in 11 patients with an astrocytoma, glioblastoma,
> meningioma or ganglioglioma (Table 3). The IFP in 10 of 11 tumours varied
> between 0.0 and 3.0 mmHg. In one astrocytoma the IFP was 9.5 mmHg. The mean
> IFP was 2.0 ± 2.5 mmHg."

### 1.1 What the ± is — the question that decides everything

**It is a STANDARD DEVIATION.** Not inferred, quoted. Statistical analysis
section, verbatim:

> "The data are presented as the mean ± standard deviation."

An SEM reading is independently impossible: SEM 2.5 at n = 11 implies SD 8.3,
which cannot be produced by 11 values of which 10 lie inside 0.0–3.0.

This is the single most important qualifier on any verdict here, and it is
exactly the situation the brief warned about. **A mean of 2.0 with an SD of 2.5,
on a quantity that is physically non-negative and that the solver clips at zero,
gives an effective admissible band of 0.0 to 4.5 mmHg.** The measured spread is
larger than the measured value. Almost nothing fails it.

### 1.2 Technique, tumour types, treatment state

* **Technique: wick-in-needle** (the paper cites Boucher et al 1991), taken
  intraoperatively "following the opening of the intracranial cavity and removal
  of the dura". Not micropuncture, not a catheter transducer.
* **Intracranial tumours specifically, not pooled** with the paper's peripheral
  or rodent work — Table 3 is a standalone human series.
* **All 11 patients were pre-treated and decompressed.** Dexamethasone 10 mg +
  mannitol 0.5 mg/kg + furosemide 20 mg i.v. before surgery, plus
  hyperventilation to PaCO₂ 25–32 mmHg. The authors' own abstract conclusion:
  *"antioedema therapy in patients with brain tumours is responsible in part for
  the low tumour IFP."* The same paper measures a 17–50 % fall in **rodent**
  tumour IFP from hypocapnia alone (27 ± 15 → 18 ± 10 mmHg). So 2.0 ± 2.5 mmHg
  is a **treated, decompressed, open-cranium** number, not an untreated
  intracranial-tumour baseline. The engine models no steroid, no osmotic agent,
  no diuretic and no decompression.

### 1.3 Table 3, transcribed in full

| Sex/age | Tumour type | Prior treatment | Tumour IFP | Surrounding brain |
|---|---|---|---:|---:|
| M/30 | R-Glioblastoma | X | 3.0 | 3.5 |
| M/39 | R-Glioblastoma | X | 1.5 | 2.5 |
| F/51 | Glioblastoma | – | 0.5 | 0.5 |
| F/48 | Astrocytoma | – | **9.5** | 1.5 |
| F/70 | Meningioma | – | 1.5 | – |
| M/48 | R-Glioblastoma | S,X | **0.0** | – |
| F/36 | R-Glioblastoma | S,X | 1.0 | 3.0 |
| M/49 | Oligo/Astrocytoma | – | 2.5 | – |
| F/30 | R-Astrocytoma | S,C,X | 2.5 | – |
| M/36 | Ganglioglioma | – | 0.5 | – |
| F/62 | R-Astrocytoma | X,C,S | 2.0 | – |

(R = recurrent; X = radiation; C = chemotherapy; S = surgery.)

Two things fall out of the table that the abstract does not give you:

* **A measured 0.0 mmHg.** The anchor's observed range literally reaches zero.
* **The glioblastoma-only subset** — the tissue the landing page actually shows —
  is n = 5, values 3.0 / 1.5 / 0.5 / 0.0 / 1.0, **mean 1.2 mmHg, range 0.0–3.0**.
  Lower still.

*Transcription check, reported rather than hidden:* the 11 values above give
mean 2.23, SD 2.59, whereas the paper states 2.0 ± 2.5. The discrepancy is small
and is most likely rounding or a single OCR digit in the PMC scan. **All scoring
below uses the paper's stated 2.0 ± 2.5**, which is also the more permissive
choice, so nothing here turns on it.

### 1.4 The cross-check that kills the anchor

According to PubMed:

> Less JR, Posner MC, Boucher Y, Borochovitz D, Wolmark N, Jain RK.
> "Interstitial hypertension in human breast and colorectal tumors."
> *Cancer Res* 1992;52(22):6371-4. **PMID 1423283** (no DOI in the PubMed record)

Same **wick-in-needle** technique, with **Boucher and Jain on both author
lists**. From its abstract, verbatim: *"IFP in normal tissues had a mean of
**2 ± 4 mm Hg (n = 11)**."*

**The intracranial-TUMOUR value (2.0 ± 2.5, n = 11) is statistically
indistinguishable from the same group's NORMAL-TISSUE value (2 ± 4, n = 11).**
An anchor that cannot separate tumour from normal tissue cannot separate this
engine from one that has no tumour in it. That is the reason the GBM row is
registered NOT TESTABLE rather than PASS.

---

## 2. What the engine's IFP field actually is

`tmesim/physics/pressure.py`, read and not modified.

**Equation** — steady-state screened Poisson (Helmholtz-type), 2-D:

    -div( K grad p ) + (a + d) * p = a * p_v

with, per grid node, `K` interstitial hydraulic conductivity, `a = vascular_density *
Lp_eff` the vessel-wall Starling source/coupling, `d = lambda_lymph * (1 - phi_tumour)`
lymphatic drainage suppressed inside the tumour, `phi` a Gaussian-smoothed
tumour-occupancy field normalised to max 1.

**Boundary conditions.** Dirichlet `p = 0` on all four domain borders **and on
every node outside the tissue mask**. The screening term provides the interior
closure. Numerics: 5-point finite volume, harmonic-mean face conductances,
`scipy.sparse` direct solve, then `clip(p, 0, None)`.

**Coefficients — all hardcoded, none fitted to anything:**

| symbol | value | units |
|---|---:|---|
| `_K0` | 1.2 | µm²/(mmHg·s) |
| `_K_ECM_ATTEN` | 0.6 | – |
| `_PV` | **25.0** | **mmHg** |
| `_LP_NORMAL` | 1.0e-3 | 1/(mmHg·s) |
| `_LEAK_GAIN` | 12.0 | – |
| `_LAMBDA_LYMPH` | 3.0e-2 | 1/(mmHg·s) |
| `_TUMOR_SIGMA` | 5.0 | **grid cells** (not µm — see §5.3) |

The module's own docstring: *"Parameters are order-of-magnitude biological
values, not patient-calibrated."* I grepped for any fit of these constants and
found none. **This makes any IFP score genuinely out-of-sample** — and equally
means nothing constrains the absolute scale.

**Units are genuinely mmHg, natively.** `_PV` is 25.0 mmHg and `K` carries mmHg
in its denominator, so `p` comes out of the solve dimensioned in mmHg with **no
conversion at the display layer**. The landing exporter passes `solve_ifp`'s
output straight through. Confirmed by reading both call sites. This is the rare
case where the engine emits a real dimensioned absolute — unlike the oxygen
field, whose mmHg↔dimensionless mapping was never derived.

---

## 3. The commensurable statistic

A **wick-in-needle** probe is a 23-gauge needle (0.64 mm OD) whose side hole,
packed with a nylon wick, equilibrates hydraulically with the *mobile
interstitial fluid* over a sub-millimetre to millimetre contact region, at **one
placement inside the tumour mass**.

So the commensurable model statistic is **the solved field's local mean over a
wick-sized window, at a tumour-bearing location, reported as a distribution over
possible placements**. Explicitly it is **NOT**:

* **not the domain mean** — a needle is not a slide average, and most of the
  section is not tumour;
* **not the field maximum** — a needle is not placed on the single hottest node;
* **not the landing page's 0.22–0.65 mmHg** — those are **iso-contour levels**,
  the 30 / 50 / 70 / 90 % fractions of a *coarse-grained maximum* (0.726 mmHg),
  chosen so four contours render legibly. They are a rendering artefact of
  `IFP_LEVEL_FRACTIONS = (0.30, 0.50, 0.70, 0.90)`, not the field's range. The
  field's actual range on the same section is **0 to 12.93 mmHg**. Scoring
  against the displayed numbers would be a category error, and I did not.

**Implementation.** Local mean over a Gaussian window of σ = 250 µm, computed by
normalised convolution over tissue nodes only (the same `coarse_grain` the
landing exporter uses, so p = 0 exterior nodes are not averaged in and no
spurious edge falloff is invented), restricted to nodes with tumour occupancy
φ ≥ 0.5, reported as the median over placements with IQR.

**Sensitivity to both free choices, reported rather than buried.** The exact
wick sampling volume is not derivable from the paper, and neither is the tumour
occupancy at which a surgeon would place the needle, so both were swept: σ = 125
/ 250 / 500 µm and φ ≥ 0.25 / 0.5 / 0.75. The window size moves the statistic by
under 1.4×. **The occupancy threshold matters more** — on GBM the median moves
0.126 → 0.462 → 0.646 mmHg across φ ≥ 0.25 / 0.5 / 0.75 (5.1×), on breast 0.307
→ 0.370 → 0.557 mmHg (1.8×) — so that choice is load-bearing and is stated
rather than assumed away.

**It does not change any verdict, because the FAIL survives the most generous
choice available.** Taking the single highest wick-scale value found *anywhere*
in the tumour, over every window size and every threshold — a statistic far more
favourable to the engine than any needle placement could be — gives **1.32 mmHg
(GBM)** and **0.94 mmHg (breast)**. The breast best case is still **6.4× below
the band floor** of 6 mmHg and **16× below** the measured mean. No defensible
commensurability choice rescues it.

---

## 4. The scores

All at the production spacing **dx = 12 µm**, the same grid the landing page
exports at.

### 4.1 GBM — `data_cache/xenium_gbm` (the landing-page section)

| | |
|---|---|
| predicted (wick-scale median) | **0.462 mmHg** (IQR 0.289–0.583, max 0.720) |
| measured | 2.0 ± 2.5 mmHg (mean ± **SD**, n = 11); GBM-only subset mean 1.2, range 0.0–3.0 |
| formally inside mean ± 1 SD? | yes (band effectively 0.0–4.5) |
| formally inside the observed 0.0–3.0? | yes |
| **registered verdict** | **NOT TESTABLE** |

**Why not PASS.** It is inside the band, and reporting that as a PASS would be
dishonest for three independent reasons:

1. **The bar cannot be failed by anything plausible.** Effective band 0.0–4.5
   mmHg on a non-negative quantity. A model that output identically zero would
   pass.
2. **The anchor does not separate tumour from normal tissue** (§1.4). Agreement
   with it is equally consistent with the engine having no tumour effect at all
   — which, on this section, is close to what it has.
3. **It is a treated, decompressed state** the engine does not model (§1.2).

And note the direction anyway: 0.462 mmHg is **4.3× below** the measured mean of
2.0 and **2.6× below** the GBM-only mean of 1.2. Even the easy test is being
passed from underneath.

### 4.2 Breast — `data_cache/xenium_breast_janesick/rep1`

This is the comparison that actually discriminates, because the breast
measurement is far from zero and ~7× its own study's normal-tissue value.

| | |
|---|---|
| predicted (wick-scale median) | **0.370 mmHg** (IQR 0.153–0.471, max 0.575) |
| measured (PMID 1423283) | **15 ± 9 mmHg** (mean ± SD), primary breast carcinoma |
| bar | inside 6–24 mmHg |
| **registered verdict** | **FAIL — 41× below the measured mean, 16× below the band floor** |

Corroborated independently:

> Nathanson SD, Nelson L. "Interstitial fluid pressure in breast cancer, benign
> breast conditions, and breast parenchyma." *Ann Surg Oncol* 1994;1(4):333-8.
> **PMID 7850532**, [DOI](https://doi.org/10.1007/BF03187139)

Invasive ductal carcinoma **29 ± 3 mmHg, n = 25** — and the paper writes
"29 ± 3 **(SE)**", so **that uncertainty is a standard error, not an SD**, and
must not be read as one. Against that value the engine is ~78× low. Its normal
parenchyma is −0.3 ± 0.1 mmHg.

This FAIL is consistent with every other absolute quantity this engine has been
scored on (D_eff 2.2–3.9× high, hindrance 7× weak, vascular area fraction 5.23×
off, MVD 3.59× off). It is a perfectly good result and it is the honest one.

### 4.3 Why it fails — the mechanism, which is worse than the number

`pressure.py`'s own docstring promises a **flat elevated plateau approaching
`_PV` = 25 mmHg** in the tumour core, with all the pressure drop squeezed into a
thin margin rim: the Jain effect. On real reconstructed tissue that is not what
happens.

| | GBM | breast |
|---|---:|---:|
| raw field max, in tissue | 12.93 mmHg | 12.43 mmHg |
| raw field **median** in tissue | **0.0036 mmHg** | **0.0055 mmHg** |
| raw field mean in tissue | 0.128 mmHg | 0.184 mmHg |
| fraction of tissue above 0.5 mmHg | 5.4 % | 9.3 % |
| fraction of tissue above 2 mmHg | 1.7 % | 2.7 % |
| fraction of tissue nodes with vascular density > 0 | 2.7 % | 5.7 % |

The solved field is **essentially zero over 94–97 % of the tissue with isolated
spikes to 10–20 mmHg at endothelial cells**. The vascular density that drives
the Starling source is zero at almost every node, so `a → 0` and `p → 0` there
regardless of tumour occupancy. **The engine reproduces roughly the right peak
magnitude and misses the plateau entirely** — which is the opposite of the
mechanism the module claims to demonstrate. The tumour-scale interstitial
hypertension that makes IFP matter for drug delivery is absent.

---

## 5. Grid sensitivity

Neither flagship panel is grid-converged (`results/grid_convergence_breast.md`,
`results/grid_convergence_gbm.md`), and IFP is solved on the same mesh. Four
spacings were run per section.

Wick-scale (σ = 250 µm) tumour median, mmHg:

| dx (µm) | GBM | breast | median L/dx (GBM) |
|---:|---:|---:|---:|
| 24 | 0.358 | 0.621 | 0.31 |
| 18 | 0.434 | 0.424 | 0.41 |
| **12 (production)** | **0.462** | **0.370** | 0.59 |
| 8 | 0.370 | 0.286 | 0.85 |

* **GBM: non-monotone, 1.29× spread**, rising then falling — the same
  non-monotone signature `grid_convergence_gbm.md` reports for the drug panel.
* **Breast: monotone decreasing, 2.17× spread**, and still falling **22.5 % over
  the last refinement** with no sign of a plateau. Extrapolated, the breast
  statistic is heading toward zero as dx → 0, which would make the FAIL *worse*,
  not better. **Refining the mesh does not rescue this number.**
* The raw field maximum falls monotonically everywhere (GBM 19.47 → 18.18 →
  12.93 → 10.81 mmHg from dx = 24 to 8), which is the signature of a
  point-source peak whose height is set by the cell size rather than by physics.

### 5.1 The mesh never resolves the operator

The screened-Poisson operator has screening length `L = sqrt(K / (a + d))`. From
the shipped coefficients, on real tissue:

| | GBM | breast |
|---|---:|---:|
| median L | **7.1 µm** | **6.2 µm** |
| L, p5–p99 | 6.6–9.8 µm | 5.5–9.1 µm |
| fraction of tissue with L < dx (at dx = 12) | **99.76 %** | **99.92 %** |
| median L / dx at production | 0.59 | 0.51 |

**The screening length is below every grid spacing the engine ever runs at
(8–30 µm), at more than 99.7 % of tissue nodes.** A screened-Poisson solve whose
screening length is under the mesh cannot represent the operator it discretises.

### 5.2 What the solve degenerates to

When the reaction term dominates the diagonal — median `(a+d)·dx²/K` = **2.87**
(GBM) and **3.80** (breast) at dx = 12 — each node relaxes toward its own *local
algebraic equilibrium* `a·p_v/(a+d)`, with the Laplacian barely coupling
neighbours. Measured directly: the solved field correlates with that pointwise
limit at **r = 0.90 (GBM) and r = 0.93 (breast)**.

**At production resolution the "solved pressure field" is closer to a pointwise
nonlinear map of the vascular-density KDE than to a transport solve.** The
margin pressure gradient that carries the Jain effect into the convection term
is not resolved. This is a real finding about the solver, independent of any
anchor.

### 5.3 A grid-dependence source inside the model itself

`_TUMOR_SIGMA = 5.0` is specified in **grid cells**, not micrometres. The
physical tumour-occupancy smoothing scale is therefore 40 / 60 / 90 / 120 µm at
dx = 8 / 12 / 18 / 24 µm. **The model changes when the mesh changes**, so the
ladder above is not a pure discretisation study — part of the movement is the
model being redefined. This is a defect in `pressure.py`, which I did not touch;
flagging it for whoever owns that file.

### 5.4 A metadata error in the exported asset

`landing-3d/public/models/ifp_contours.json` states the raw field's structure is
"per-vessel (screening length ~60 um)". The measured median screening length is
**7.1 µm** — the published figure is ~8.5× too large. The ~60 µm number more
plausibly describes the vascular-density KDE bandwidth, which is a different
quantity. I do not own `landing-3d/` or `experiments/export_landing_asset.py`;
flagging it for their owner.

---

## 6. Registration

Registered in `tmesim/validation/scorecard.py` as full-tier PDE-scope group
**`ifp-wick-in-needle`** (~15 s, needs both Xenium caches), three rows:

| row | status |
|---|---|
| GBM section IFP vs human intracranial-tumour IFP (PMID 9062403) | **NOT TESTABLE** — anchor cannot discriminate; full reason on the row |
| breast section IFP vs human primary-breast-carcinoma IFP (PMIDs 1423283, 7850532) | **FAIL** — 41× below the measured mean |
| IFP solve is under-resolved on every grid the engine runs at | **NOT TESTABLE** — no measurement; registered so the two rows above cannot be quoted without it |

Three new named conditionals: `ifp-uncalibrated` (every coefficient hardcoded
and unfitted — which is what makes the score out-of-sample *and* what leaves the
absolute scale unconstrained), `ifp-unresolved-mesh` (§5.1–5.2),
`ifp-grid-unconverged` (§5).

The NOT TESTABLE row exists specifically so nobody re-derives this question and
lands on an easy PASS. That is now closed.

---

## 7. Recommendation: the landing page should stop showing a numeric mmHg range

**Recommendation: remove the `0.22–0.65 mmHg` figure from section 02's label.
Keep the visual, keep "Solved", drop the number.** (I do not own `landing-3d/`
or `app/`; this is a recommendation to their owner, not a change.)

The reasoning, in order of severity:

1. **The engine's absolute IFP is wrong where it can be checked.** Against the
   only IFP measurement that discriminates, the engine is **41× low**. A
   specific mmHg figure presented as a credential is a claim the project cannot
   support.
2. **The displayed numbers are not what a reader will think they are.** A
   visitor reads "0.22–0.65 mmHg" as *the tumour's interstitial pressure range*.
   They are the 30 % and 90 % iso-contour levels of a 250 µm coarse-grained
   maximum. The field's real range on that section is 0–12.93 mmHg and its
   median is 0.0036 mmHg. Even if the physics were right, the label describes a
   rendering parameter, not a result.
3. **It is not grid-converged.** The underlying statistic moves 1.3× (GBM) to
   2.2× (breast) across the grid ladder. Re-exporting at a different dx would
   change the published number with no change in the science.
4. **The mesh does not resolve the operator** (§5.1). The specific digits are
   downstream of a discretisation that cannot represent the equation.

**What is still honest to show.** The *qualitative* asset is defensible and
worth keeping: the contours are genuinely the output of a real steady-state
solve on measured tissue geometry, the JSON already says
`is_calibrated_to_patient: false`, and the spatial *structure* — pressure
concentrated around perfused regions — is a real consequence of measured vessel
positions. The recommended label is something like **"Solved interstitial
pressure — relative, uncalibrated — 250 µm mean"**, with the mmHg values dropped
from both the label and the tooltip. If a number is wanted for credibility, the
honest one to show is not a pressure at all: it is that the field is a solve over
a measured geometry, not an illustration.

**What would make it displayable again.** A pressure model whose screening
length is resolved at production dx, whose tumour plateau survives on real
reconstructed tissue, and which lands inside 6–24 mmHg on the breast section.
None of that is a tuning exercise, and it must not become one: `_PV`, `_LEAK_GAIN`
and `_LAMBDA_LYMPH` must not be adjusted toward 15 mmHg, or the anchor stops
being an anchor. **Nothing was tuned in producing this report.**
