# The sink's functional form does NOT explain the missing hypoxia — it overshoots it, and the two forms' pO₂ supports do not overlap

**Date: 2026-09-20. Status: PRE-REGISTERED A/B, REFUTED on the registered threshold.**
Nothing under `tmesim/` was modified. No default was changed. No parameter was tuned
to any anchor. `OxygenConfig.sink_form` still ships as `"linear"` and `consumption`
still ships at `0.02`.

**Calibration guard, start of pass:**

```
$ .venv/bin/python experiments/oxygen_sink_form.py
START calibrate_alpha() = 22.223477711172563
ASSERT OK
```

*(The closing assertion is in §10. A transient failure of this guard, and its cause,
is reported in §9.1 rather than suppressed.)*

---

## 0. The headline, which is a refutation

`results/oxygen_calibration_landscape.md` §4 argued on paper, on 2026-07-28, that the
engine's first-order-everywhere oxygen sink *"will systematically under-predict
hypoxia"*, and that a saturating Michaelis–Menten sink is the literature-standard form.
The engine reports `HYPOXIC = 0.0000` on the flagship GBM section — the predicted
signature. `results/stale_blockers_2026-09-20.md` §4 ranked testing this **★1**, the
single highest-value verification available.

It was tested. **The pre-registered refutation branch fired.**

| | linear (ships) | michaelis_menten | measured (human glioma) |
|---|---:|---:|---|
| median pO₂, tumour cells | **12.0815 mmHg** | **1.0419 mmHg** | **5.1 mmHg** (Lally 2006) |
| **T1** — median in [1.5, 5.6] mmHg | **FAIL (HIGH)** | **FAIL (LOW)** | acceptance band = Nozue 1997 |
| P(pO₂ < 2.5 mmHg) | **0.0000** | 0.9780 | > 0 (Collingridge 1999) |
| **T2** — a sub-2.5 mmHg tail exists | **FAIL** | **PASS** | |
| minimum pO₂ anywhere a tumour cell sits | 10.1669 mmHg | 0.3370 mmHg | |
| maximum pO₂ anywhere a tumour cell sits | 14.7473 mmHg | **3.7775 mmHg** | |

**The two forms do not bracket the measurement in the ordinary sense — their supports
are disjoint.** The Michaelis–Menten field's **maximum** (3.7775 mmHg) lies *below* the
linear field's **minimum** (10.1669 mmHg). The measured cohort median of 5.1 mmHg falls
in the empty gap between them.

> **At the shipped `consumption`, neither sink form places a single tumour cell at the
> measured median pO₂.** Not one of 23,812, under either form.

So the functional form is **not** the explanation. Swapping it does not move the field
onto the measurement; it moves the field past the measurement and out the other side,
from a section with no hypoxia at all to a section that is uniformly more anoxic than
the measured median. **This closes, as a sufficient explanation, a hypothesis this
repository has carried on paper since 2026-07-28.**

**What the form *does* explain — and this is the one part that survives, §7.** The
linear sink's distribution is too *narrow*, by a factor no unit choice can fix: it
reaches only 84.15 % of the way down from its own median to its own minimum, where the
measured tissue must reach 49.02 %. **The linear sink cannot produce a sub-2.5 mmHg
region at any value of `REFERENCE_PO2_MMHG` whatsoever.** The Michaelis–Menten sink
can. So the July derivation was right about the *mechanism* (the form controls whether
a hypoxic tail is reachable at all) and wrong about the *consequence* (switching the
form does not reconcile the engine with the measurement).

---

## 1. A finding about the recommendation, before the finding about the physics

**The ★1 recommendation this task acted on was itself stale.** `results/stale_blockers_2026-09-20.md`
row 4 correctly identified that the capability shipped on 2026-07-29, one day after
`oxygen_calibration_landscape.md` declined the test. What it did not check is whether
the **test** had also been run. It had — the next day, in the same burst of work.

| what | where | when | status |
|---|---|---|---|
| Both sink forms, real flagship section, hypoxic fraction + mean O₂ | `results/sweep/oxygen_functional_form_prototype_2026-07-29T05-54-00Z.md` (Picard) | 2026-07-29 | **already run** |
| Same, production Newton solver, 11-point `consumption` grid | `results/oxygen_sensitivity_band.md` §3 | 2026-07-29 | **already run** |
| Linear arm's full per-cell mmHg percentile ladder + sub-2.5/sub-10 fractions | `results/sweep/oxsupply_sections.json` → `sections/gbm/shipped_0.05` | 2026-08-10 | **already computed, never published in any `.md`** |
| Linear arm's median mmHg across a dx ladder | `results/oxygen_joint_grid_test.md` §1.4 | 2026-08-15 | **already run** |
| **Any mmHg number under `sink_form="michaelis_menten"` on a real section** | — | — | **did not exist** |
| **A declared-in-advance PASS/FAIL band for oxygen** | — | — | **did not exist** (`grep "pre-regist" results/oxygen*.md` → 0 hits) |

This is the same defect `stale_blockers` was written to catch, one turn deeper: it
scanned for *declined tests whose capability now exists*, and did not scan for
*declined tests that were subsequently run*. The cost is the same in both directions —
work proposed that has already been done.

**The honest consequence for this document: the headline A/B was NOT pre-registerable,
and is not presented as though it were.** §3 declares, line by line, every number that
was read before measuring.

---

## 2. The two equations

Both are solved by `tmesim/physics/oxygen.py:417 solve_oxygen` on the same 5-point
finite-difference Laplacian with no-flux (Neumann) boundaries, `dx = 12 µm`.

**Common form.** Steady-state reaction–diffusion for dimensionless saturation
`O ∈ [0,1]`:

```
D·∇²O  +  supply·vasc·(1 − O)  −  sink(O)  =  0
```

with `D = 2000 µm²/s`, `supply = 0.05`, `vasc = micro.fields.vascular_density`, and
`ρ_tum` the tumour-cell **areal** density in cells per `TUMOR_REFERENCE_AREA_UM2 = 144 µm²`
(`oxygen.py:238 _tumor_density`), which is grid-spacing independent.

**Arm A — `sink_form="linear"` (the shipping default).**

```
sink(O) = consumption · ρ_tum · O                       consumption = 0.02
```

Linear in `O` everywhere, so the system collapses to one symmetric, diagonally
dominant sparse system solved directly, with no iteration:

```
(−D·L + diag(supply·vasc + consumption·ρ_tum)) @ O  =  supply·vasc
```

The sink vanishes as `O → 0`: **consumption switches itself off exactly where the
tissue is running out of oxygen.** That is the structural property under test.

**Arm B — `sink_form="michaelis_menten"` (opt-in; Grimes, Fletcher & Partridge 2014,
PMID 26064525).**

```
sink(O) = consumption · ρ_tum · O / (O + K)             K = mm_half_saturation = 1/35
```

`K = 1/35 ≈ 0.028571` is Grimes' base case `Km = 1 mmHg` expressed against this
project's own `O = 1 ≡ 35 mmHg` anchor. The sink saturates toward a **constant**
`consumption·ρ_tum` as `O` grows and, crucially, **does not switch off as `O → 0`** —
it falls only linearly through the `O/(O+K)` factor once `O ≲ K`. Nonlinear, so it is
solved by damped Newton–Raphson with an analytic Jacobian
`J = M + diag(reaction_bulk·K/(O+K)²)` and a backtracking line search
(`oxygen.py:313 _solve_michaelis_menten`), reusing the same sparse direct solver at
every step. It raises rather than returning a silently unconverged field.

**The mmHg map.** `oxygen.py:225 to_mmhg` is a pure linear rescale,
`pO₂ = O × REFERENCE_PO2_MMHG` with `REFERENCE_PO2_MMHG = 35.0` (`oxygen.py:215`),
anchored to normal non-tumour-bearing human brain **tissue** pO₂ (Erecinska & Silver
2001, PMID 11718758). There is no offset, which is what makes the shape statistic in §7
anchor-invariant. The engine's own hypoxia band `O < 0.2` is therefore `< 7 mmHg`.

*(This is, as far as the repository's own grep can tell, the first executable call site
of `to_mmhg` in committed code. Every prior mmHg conversion in `experiments/` multiplied
`REFERENCE_PO2_MMHG` by hand.)*

---

## 3. What was known before measuring — the disclosure that makes §4 honest

A pre-registration is worthless if the registrant already knows the answer. These
numbers were read **before** `experiments/oxygen_sink_form.py` was written, and are
recorded in that file's `PRIOR_INFORMATION` block:

| quantity | published value | source | read before measuring? |
|---|---:|---|---|
| linear hypoxic fraction @ `consumption=0.02` | 0.0000 | `oxygen_sensitivity_band.md` §3 | **yes** |
| MM hypoxic fraction @ `consumption=0.02` | 1.0000 | `oxygen_sensitivity_band.md` §3 | **yes** |
| linear mean O₂ (tumour) | 0.3486 | `oxygen_sensitivity_band.md` §3 | **yes** |
| MM mean O₂ (tumour) | 0.0347 | `oxygen_sensitivity_band.md` §3 | **yes** |
| linear median pO₂ | 12.107 mmHg | `oxygen_joint_grid_test.md` §1 | **yes** |
| linear full percentile ladder + P(<2.5) = 0.0 | p1 10.219 … p99 14.707 | `oxsupply_sections.json` | **yes** |

**Consequence, stated plainly.** The linear arm is **entirely spent**. Every quantity
this document reports for it was on disk beforehand, so the linear arm below is a
**reproduction**, not a measurement. Because its hypoxic fraction (`< 7 mmHg`) was
already published as exactly 0.0, its sub-2.5 and sub-5 mmHg fractions were *derivable*
as 0.0 and were likewise not open.

**What was genuinely open, and is therefore genuinely pre-registered:** the
Michaelis–Menten arm's **distribution in mmHg on a real reconstructed section** — its
median, its quantile ladder, its minimum and maximum, and how its 100 % sub-7 mmHg mass
splits across the 2.5 and 5 mmHg thresholds. No such number existed anywhere in the
repository. The only both-forms mmHg execution on record
(`stale_blockers_2026-09-20.md` §2.3) was on a synthetic 64×64 lattice and says so in
its own text.

---

## 4. The pre-registration

Frozen in `experiments/oxygen_sink_form.py` as the `PRE_REGISTRATION` dict **before the
file was ever executed**. Thresholds are stated in both directions so a null is as
interpretable as a hit.

**Section.** `reconstruct_real('xenium', dx_um=12.0)` — the flagship GBM Xenium section
at the grid-converged production spacing. Grid 424×675, 40,887 cells, **23,812 tumour
cells**, 193,752 tissue grid points. `dx` frozen at 12 µm in advance because
`oxygen_joint_grid_test.md` §1.4 shows this statistic crossing the Lally value as `dx`
moves.

**RunSpec.** Arm A `OxygenConfig()`. Arm B `OxygenConfig(sink_form="michaelis_menten")`.
Shared and unchanged in both arms: `diffusivity=2000.0`, `consumption=0.02`,
`supply=0.05`, `hypoxia_threshold=0.2`, `mm_half_saturation=0.02857142857142857`,
`mm_max_iter=50`, `mm_tol=1e-09`. **No parameter was tuned to any anchor.**

**Anchors — all four fetched from the primary PubMed record in this pass, not inherited.**

| id | source | fetched quote | role |
|---|---|---|---|
| **A1** | Lally BE *et al.*, *Cancer J* 2006;12(6):461-6, PMID [17207315](https://pubmed.ncbi.nlm.nih.gov/17207315/) | *"Tissue oxygen tension was measured intraoperatively using an Eppendorf pO2 Histograph. … we included 13 patients with high grade gliomas (HGG) and 10 with low grade gliomas (LGG). **The median tumor oxygen pressure for the entire cohort was 5.1 mmHg.**"* | the measured central value |
| **A2** | Collingridge DR *et al.*, *Radiother Oncol* 1999;53(2):127-31, PMID [10665789](https://pubmed.ncbi.nlm.nih.gov/10665789/) | *"This study quantifies the spatial distribution of pO2 in glioma and in the surrounding brain tissue. **Both glioma and peritumoural brain contain regions at oxygen tensions less than 2.5 mmHg.**"* | the tail must exist |
| **A3** | Nozue M *et al.*, *J Surg Oncol* 1997;66(1):30-8, PMID [9290690](https://pubmed.ncbi.nlm.nih.gov/9290690/) | *"**the median values of tumor pO2 varied from 1.5 mmHg to 5.6 mmHg among the laboratories**, and pO2 of normal subcutis also varied from 28 mmHg to 38 mmHg. … Caution needs to be exercised in using the absolute, median, or distribution of pO2 values"* | the acceptance band |
| **A4** | Pinel S *et al.* 2004, PMID [15093922](https://pubmed.ncbi.nlm.nih.gov/15093922/); Stüben G *et al.* 1998, PMID [9783885](https://pubmed.ncbi.nlm.nih.gov/9783885/) | *"Percentage of pO(2) values <2.5 mm Hg … 58.5 +/- 27.0% for GBM Nan1; … 91.5 +/- 8.3% … for U87"*; *"about 80% of all pO2 values were < or =2.5 mmHg"* | **supporting only** |

**A correction to the inherited reading of A1, found by re-fetching it.** The abstract
dichotomises *"the median pO2"* per patient and stratifies *"patients … with a median
tumor oxygen tensions either above or below"* the cut-points. **5.1 mmHg is a
cohort-level median OF PER-PATIENT MEDIANS, not the median of a pooled per-reading
histogram.** That is a different object from the engine's per-cell distribution over one
section, and §8 does not pretend otherwise. `verification_routes_2026-09-20.md` R9 and
`oxygen_supply_calibration.md` §4.1 both quote the 5.1 mmHg without this distinction.

**A4 is supporting only, and never decisive**, because the host vasculature is murine
and the tumour is a cell line. Re-fetching the primary records established something the
repository did not have on file: **no published per-reading fraction-below-threshold
exists for human *patient* glioma in the open abstracts.** Collingridge's abstract gives
no median, no *n* and no percentage; Rampling's (PMID [8005794](https://pubmed.ncbi.nlm.nih.gov/8005794/),
*"Six different tracks, each with a path length of 22 mm, were recorded per patient
representing 192 readings"*, *"Both the high grade intrinsic and metastatic brain tumors
showed significant regions of hypoxia"*) gives neither. **That is the binding limit on
how distributional this comparison can be**, and it is a property of the literature, not
of the engine.

**Thresholds.**

| id | rule | provenance |
|---|---|---|
| **T1** | An arm PASSES iff its per-tumour-cell median pO₂ ∈ **[1.5, 5.6] mmHg**. Direction of failure recorded. | A3's measured inter-laboratory span for an identical protocol on identical tissue — the measurement's own noise floor. It contains A1's 5.1 mmHg. |
| **T2** | An arm PASSES iff **P(pO₂ < 2.5 mmHg) > 0** among tumour cells. | A2. Deliberately the weakest possible form — one cell suffices — so failing it cannot be argued down. |
| **T3** | The form **EXPLAINS** the gap iff A fails T1 or T2 **and** B passes **both**. | — |
| **T4** | The forms are **INDISTINGUISHABLE** iff \|Δmedian\| < 0.5 mmHg **and** \|ΔP(<2.5)\| < 0.05. | 0.5 mmHg ≈ 12 % of the T1 band width. Declared so a null is interpretable. |
| **T5** | **REFUTED-BY-OVERSHOOT** iff B passes T2 but fails T1 **LOW** while A fails T1 **HIGH** — the forms bracket the measurement instead of one matching it. | Declared in advance *specifically* so this outcome could not later be rationalised as "the MM form is directionally right". If T5 fires, it leads the report. |

**Confounds, also declared in advance** (`CONFOUNDS` in the script): the source constant
is independently known to be ~4,000× off (§9.2); the mmHg anchor is itself contested
(§6); the median statistic is mesh-sensitive (`dx` frozen); and cell-weighted sampling
is not needle-track sampling (both reported, §5).

---

## 5. The measured distributions

pO₂ in mmHg via `to_mmhg`, `REFERENCE_PO2_MMHG = 35.0`. Per-tumour-cell (n = 23,812) is
the primary sample; the area-weighted tissue grid (n = 193,752) is reported alongside it
rather than instead of it, because an Eppendorf needle samples space, not cells.

### 5.1 Quantile ladders, tumour cells

| statistic | **linear** | **michaelis_menten** |
|---|---:|---:|
| min | **10.1669** | **0.3370** |
| p1 | 10.3358 | 0.4127 |
| p5 | 10.5605 | 0.5244 |
| p10 | 10.8156 | 0.6100 |
| p25 | 11.4856 | 0.8170 |
| **p50 (median)** | **12.0815** | **1.0419** |
| p75 | 12.5833 | 1.4238 |
| p90 | 13.8887 | 2.0358 |
| p95 | 14.1477 | 2.2701 |
| p99 | 14.5439 | 2.9244 |
| max | **14.7473** | **3.7775** |
| mean ± sd | 12.1548 ± 1.0264 | 1.1882 ± 0.5516 |

### 5.2 The same, area-weighted over the tissue grid

| statistic | **linear** | **michaelis_menten** |
|---|---:|---:|
| min | 10.1669 | 0.3370 |
| p10 | 11.0520 | 0.6742 |
| **p50** | **12.1535** | **1.0980** |
| p90 | 14.0316 | 2.1270 |
| max | 14.7680 | 3.7992 |
| mean ± sd | 12.2731 ± 1.0292 | 1.2601 ± 0.5724 |

The two weightings differ on the median by **0.60 %** under the linear form and
**5.39 %** under Michaelis–Menten — the nonlinear sink is the more sampling-sensitive of
the two, which is worth recording. Neither is close to the 11.04 mmHg difference the A/B
is about, and **neither changes any verdict**: the area-weighted Michaelis–Menten median
of 1.0980 mmHg is still below the T1 band, and the area-weighted linear median of
12.1535 mmHg is still above it. So **the cell-vs-needle sampling mismatch (confound C4)
does not drive any result here**, but it is not negligible under the MM arm and should
not be quoted as if it were.

### 5.3 Empirical CDF — the distribution-to-distribution comparison, as far as the literature permits

P(pO₂ < x), tumour cells:

| x (mmHg) | 0.5 | 1.0 | **2.5** | **5.0** | **7.0** | **10.0** | 14.0 | 20.0 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| **linear** | 0.0000 | 0.0000 | **0.0000** | **0.0000** | **0.0000** | **0.0000** | 0.9157 | 1.0000 |
| **michaelis_menten** | 0.0378 | 0.4498 | **0.9780** | **1.0000** | **1.0000** | **1.0000** | 1.0000 | 1.0000 |
| *measured constraint* | — | — | **> 0** (A2) | — | — | — | — | — |

**This is a four-point CDF comparison, not a two-sample distributional test, and the
reason is the literature, not the engine.** The primary human-glioma records publish
summary statistics, not histograms (§4). The one point on the measured CDF that is
firmly established for human patient glioma is `P(< 2.5) > 0`, and the one measured
central value is the cohort median of per-patient medians. Everything stronger would
require per-reading data that the open abstracts do not carry. Reporting a Kolmogorov–
Smirnov statistic against a distribution nobody has published would be an invented
number.

Against the **supporting** xenograft range (A4: 58.5 %, ~80 %, 91.5 % of readings below
2.5 mmHg), the Michaelis–Menten arm's **97.80 %** sits **above all three**, i.e. it is
more hypoxic than the most hypoxic published glioma xenograft control. The linear arm's
**0.00 %** sits below all three. Both overshoot the supporting anchor, in the same two
directions they overshoot the primary one.

### 5.4 Reproduction check against the already-published numbers

Not a new result — a check that this run is the same experiment.

| quantity | published | measured now | Δ |
|---|---:|---:|---:|
| linear hypoxic fraction | 0.0000 | 0.000000 | 0 |
| MM hypoxic fraction | 1.0000 | 1.000000 | 0 |
| linear mean O₂ (tumour) | 0.3486 | 0.347281 | −0.38 % |
| MM mean O₂ (tumour) | 0.0347 | 0.033950 | −2.16 % |
| linear median pO₂ | 12.107 mmHg | 12.0815 mmHg | −0.21 % |

The binary quantities reproduce exactly; the continuous ones sit 0.2–2.2 % low. **A
plausible cause, flagged as a hypothesis and not a verified attribution:** all three
published values predate the promotion of `discrete_vascular_field(geometry='area')` to
the shipping default in mid-August 2026 (`results/area_conservation_promotion.md` was
added 2026-08-15; `oxygen_joint_grid_test.md`'s own staleness banner is dated
2026-08-14), which changed the vessel source field these solves read. I did not test
that attribution by re-running under the old painting, and it is not needed
for any verdict here — the effect is two orders of magnitude smaller than the 11.04 mmHg
difference the A/B is about.

---

## 6. Controlled: the proof that only the sink form changed

Asserted at runtime, not claimed in prose.

**Config identity.** `dataclasses.asdict` diff between the two `OxygenConfig` objects:

```
CONFIG DIFF (must be exactly one field): {'sink_form': ('linear', 'michaelis_menten')}
  ASSERT OK -- the two arms differ in exactly one field.
```

**Input identity (CONTROL_3).** SHA-256 over everything `solve_oxygen` reads — domain
`nx`/`ny`/`dx`, the vascular density array, cell `x`, `y`, `cell_type`, and the derived
tumour areal density from `_tumor_density` — taken three times:

```
CONTROL_3 fingerprint BEFORE arm A: ea666e06b655e90a
CONTROL_3 fingerprint BETWEEN arms: ea666e06b655e90a
CONTROL_3 fingerprint AFTER arm B : ea666e06b655e90a
CONTROL_3 all three identical: True
```

One `Microenvironment` object is built once and handed to both solves. The inputs are
byte-identical, demonstrated rather than asserted.

**Determinism (CONTROL_1).** Arm A solved twice is bit-identical (`np.array_equal` →
`True`).

**Degeneracy (CONTROL_2) — the control that makes the comparison mean anything.** As
`K → ∞`, `c·K·O/(O+K) → c·O`: the Michaelis–Menten sink degenerates analytically to the
linear one. Run at `K = 10⁴` with `consumption = 0.02 × 10⁴ = 200`:

```
max|O_MM - O_linear| = 9.566e-06        CONTROL_2 PASS
```

The Newton solver reproduces the direct linear solve to 1×10⁻⁵ when the mathematics says
it must. **The 11.04 mmHg difference between the arms is therefore attributable to the
functional form and not to the solver.** Solve times: linear 5.67 s, Michaelis–Menten
39.10 s (Newton, converged, no `RuntimeError`).

---

## 7. The part of the July derivation that survives: an anchor-invariant shape test

**Flagged explicitly as POST-HOC. This was devised after seeing the two distributions
and is reported as a diagnostic, never as a registered test.** It is in the script's
`post_hoc_diagnostics` block under an `IS_POST_HOC_NOT_PREREGISTERED: True` key.

`to_mmhg` is a pure rescale with no offset, so the ratio `R = min / median` is
**invariant under every possible choice of `REFERENCE_PO2_MMHG`**. Combining the two
primary human anchors — A1's median of 5.1 mmHg with A2's requirement that regions below
2.5 mmHg exist — the measured tissue must reach `2.5 / 5.1 = 0.4902` of its own median.

| arm | R = min / median | can reach a sub-2.5 mmHg region at **some** anchor? |
|---|---:|---|
| **linear** | **0.8415** | **NO** |
| **michaelis_menten** | **0.3235** | **YES** |

**This is the one statement in this document that no choice of `REFERENCE_PO2_MMHG` can
change.** The linear sink's field is too *narrow*: across an entire GBM section it
spans only 10.1669–14.7473 mmHg, a factor of 1.45 from end to end. Rescaling moves that
window but cannot stretch it. **A first-order-everywhere sink structurally cannot produce
a hypoxic tail on this section at any unit convention** — which is precisely the
mechanism `oxygen_calibration_landscape.md` §4 derived on paper, now measured.

So the correct reading is a split verdict, and both halves matter:

* **On the mechanism — the July derivation is CONFIRMED.** The sink's functional form
  is what controls whether a hypoxic tail is reachable at all, and the shipped form
  makes it unreachable.
* **On the consequence — the July derivation is REFUTED.** *"Therefore switching the
  form will fix the hypoxia gap"* does not follow and is false. The form is necessary
  and nowhere near sufficient: it converts an error of one sign into a larger error of
  the other sign, and at the shipped `consumption` it never passes through the measured
  value at all.

A hypothesis that was carried for 54 days as a single proposition turns out to be two
propositions with opposite truth values. **Separating them is the result.**

---

## 8. Robustness: does the verdict depend on the contested mmHg anchor?

Declared in advance, because `results/oxygen_supply_calibration.md` §2 argues from
Secomb's reference cortical values (`P̄_cap = 58.75`, `P̄_tiss = 47.65` mmHg) that *"a
compartment-correct anchor for `O = 1` would be ~55–60 mmHg, not 35, and every mmHg
figure in this report would scale up by ~1.7×"*. The primary verdict is scored at the
**shipped, frozen** 35.0. The secondary pass re-scores at 57.5 mmHg (×1.643):

| arm | median @ 35.0 | T1 | median @ 57.5 | T1 @ 57.5 | `p_ref` it would NEED to hit 5.1 mmHg |
|---|---:|---|---:|---|---:|
| linear | 12.0815 | **FAIL (HIGH)** | 19.8482 | **FAIL (HIGH)** | **14.8 mmHg** |
| michaelis_menten | 1.0419 | **FAIL (LOW)** | 1.7116 | **PASS** | **171.3 mmHg** |

**Reported as a partial failure of robustness, not swapped.** The linear arm's T1 failure
is anchor-independent across the entire defensible range. **The Michaelis–Menten arm's
T1 failure is not** — at the contested capillary anchor its median enters the band. Any
future document quoting §0's table must carry this row with it.

**What is *not* rescued by the anchor, and is the stronger argument:** the last column.
For the linear field to place its median on the measured 5.1 mmHg, `O = 1` would have to
mean **14.8 mmHg** — below normal brain tissue pO₂, i.e. "fully oxygenated" would denote
a tension the literature calls hypoxic. For the Michaelis–Menten field, `O = 1` would
have to mean **171.3 mmHg** — roughly 1.7× arterial. The defensible range for this anchor
is about 35–60 mmHg. The two arms need 14.8 and 171.3, missing it low by 2.4× and high by
2.9×. **Neither field is rescuable by any physically defensible unit convention**, and
the T2 result is untouched by the anchor in any case (§7).

---

## 9. What this does not show

**9.1 — The calibration guard failed transiently at the start of this pass, and that is
reported rather than suppressed.** The first two invocations returned
`32.804763617621894` and `22.224082431732104` instead of the pinned
`22.223477711172563`. The cause was diagnosed without touching a file: another agent
held uncommitted edits in `tmesim/molecular/`, and the tree was read mid-write. Three
consecutive invocations then returned the pinned value, and both the START and END
asserts in the run of record passed. **This measurement is independent of that
territory**: tracing `sys.modules` shows the oxygen path imports only
`tmesim.molecular`, `tmesim.molecular.mgmt` and `tmesim.molecular.types` — not
`two_route`, the edited module — and `solve_oxygen` consumes no drug parameter and no
`alpha`. Nothing was reverted, stashed or checked out. But a reader should know the
guard was not green on first contact.

**9.2 — Both arms sit on a source constant independently known to be ~4,000× off.**
`results/oxygen_supply_calibration.md` derives `supply ≈ 194 /s` against the shipped
`0.05`, and states that *"Every published `F_state`, absolute kill and kill-derived
ratio on the GBM section should be read as conditional on an oxygen source constant that
is ~4,000× from its measured value."* **Both arms of this A/B inherit that.** So no
arm's absolute agreement or disagreement with a measured pO₂ is attributable to the sink
alone. What the A/B *does* isolate is the **difference** between the forms, because the
source is byte-identical across them (§6). This document does not license "the sink form
is the problem with the oxygen field" — it licenses "the sink form is not the
explanation for the missing hypoxia."

**9.3 — The two axes are coupled and only one was varied.** `consumption` was held at
0.02 in both arms, deliberately, so the comparison would be controlled. But
`oxygen_sensitivity_band.md` §3 already established that the forms reach saturation at
`consumption` values 5–10× apart, and §6 of that document states the axes are coupled.
**Nothing here tests whether some `(sink_form, consumption)` pair reproduces the
measurement.** A Michaelis–Menten sink at a lower `consumption` might well land in the
band; this run does not know, and finding such a pair by searching would be tuning a
parameter to an anchor, which is prohibited.

**9.4 — This is n = 1 in tissue.** One section, one patient, one platform. Seven
sections are already solved and on disk (`results/multisection_generalisation.json`) and
30 more are ingested; none was used here. Nothing in this document establishes that the
disjoint-support finding generalises.

**9.5 — `dx` was frozen, not swept.** 12 µm, chosen in advance. `oxygen_joint_grid_test.md`
§1.4 shows this statistic moving across the Lally value with mesh spacing for other
configurations. Whether the disjointness in §0 survives a `dx` ladder is untested.

**9.6 — The measured comparison is four points on a CDF, not a distribution test**, and
one of those points (A1) is a median of per-patient medians rather than a per-reading
median (§4). The engine's per-cell distribution over one section and Eppendorf's
per-reading distribution along needle tracks across 23 patients are different objects.
The disjointness in §0 is robust to this — a maximum below a measured median is a
mismatch under any reasonable pairing — but the *degree* of mismatch is not quantified
to better than that.

**9.7 — No claim here is clinical, diagnostic, predictive of response, or about any
individual.** No cross-compound ranking is made or implied; that claim class remains
withdrawn (2026-08-21) and absolutely. No absolute kill fraction is reported. Nothing
here is "validated". Nothing here recommends changing `sink_form`'s default, and this
document is not evidence that the Michaelis–Menten form is correct — §0 shows it
failing the registered median threshold in the opposite direction.

**9.8 — Two of the five anchors are not human patient glioma.** A3 (the acceptance band)
is a murine fibrosarcoma, used only as the measurement's own noise floor. A4 is human
glioma cell lines xenografted in mice, used only as supporting evidence and never
decisive. Both limits are stated at point of use.

---

## 10. Reproduction and closing calibration guard

```
$ ./.venv/bin/python experiments/oxygen_sink_form.py
```

Deterministic, no seeds; ~1 minute. Writes `results/oxygen_sink_form_data.json` with
every number above plus the full `PRIOR_INFORMATION`, `CONFOUNDS`, `ANCHORS` and
`PRE_REGISTRATION` blocks, the three fingerprints, and the post-hoc diagnostics.
Run at `git HEAD = f9d443d56f073b7b840daaf09c4549e74d89b2fc`.

```
END calibrate_alpha() = 22.223477711172563
ASSERT OK
```

Identical to the opening value in the banner. **No file under `tmesim/`, `app/`,
`docs/`, `scripts/` or `data_cache/` was modified in this pass, and no existing file in
`results/` or `archive/` was created, edited or touched.** This pass wrote exactly three
new files: `experiments/oxygen_sink_form.py`, this document, and
`results/oxygen_sink_form_data.json`.
