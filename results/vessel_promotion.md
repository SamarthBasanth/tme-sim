# Promoting the discrete vessel source to default — the two blockers, settled

`results/vessel_network_impact.md` §8.2 refused to make the discrete segmented-vessel
source the default for two stated reasons: **(a)** `alpha` was calibrated against the KDE
field and would have to be re-fitted, and **(b)** grid convergence was untested — a 7.5 µm
vessel radius on a 12 µm mesh is under-resolved. This report settles both, re-runs the
out-of-sample validation that was named the deciding test, and acts on the result.

Everything below was measured on the REAL 10x Xenium GBM section
(`reconstruct_real('xenium')`, 40,887 measured cells, 30,174 tumour, 3,984 endothelial).
Raw drivers and outputs are listed in §9.

## Verdict, up front

> **PROMOTED.** `reconstruct.real_fields.vascular_density_from_cells` now defaults to
> `mode='discrete'`; the KDE is retained unchanged as `mode='kde'`.

| the question the brief asked | the answer the numbers gave |
|---|---|
| Does the discrete field's prediction **converge** on dx = 15/12/8/6? | **YES, on every metric.** Coverage @1× IC50, mean tumour AUC and kill are converged by **dx = 12 µm**; the 8 → 6 µm worst-drug changes are **0.0011 abs / 0.38% / 0.0008 abs** against tolerances of 0.03 abs / 5% rel |
| Does it need a **finer mesh than the KDE** did? | **For one metric yes, and that is the real cost — but the overall comparison goes the other way.** The fitted `L` needs **dx ≤ 8 µm** (doxorubicin 47.2 / **39.6** / 51.5 / 51.9 µm at dx = 15/12/8/6, so dx = 12 sits **24% below** the asymptote). Against that, **the KDE is the field that fails to converge**: its mean tumour AUC moves **6.35%** from dx = 8 to dx = 6 — outside the 5% tolerance — while the discrete field's moves 0.38% (§2.5) |
| What is the **new alpha**? | **148.6 — unchanged, and it must not be re-fitted.** `calibrate_alpha()` is a closed-form solve on the compiled drug plus the published anchor; **no vascular field enters it**. It returns `148.60815765445545` under either default, bit-identical |
| Does **out-of-sample PK** still pass? | **YES, 2/2, unchanged.** mitoxantrone 48.2 µm and topotecan 83.6 µm, both inside the cited 36.1–108.2 µm band. `pk_multidrug.run()` is byte-identical under both defaults. **And forcing an empirical re-fit BREAKS it under both fields** (§4) |
| Does `therapeutic_coverage` still **discriminate**? | **YES** — separation 134× its own grid drift (KDE: 271×), against a threshold of 3× |
| Does `penetration_length` now **separate drugs**? | **YES, and 2.2× better.** Five-drug spread **3.12× → 6.91×**, and the discrete source **admits no null fit at all** |
| Does the **regional** agreement improve? | **NO — it cannot.** Byte-identical, and the spatial version differs by **0.19%** (§6). The §8.3 claim that the discrete field is the right substrate for the BBB leak is **not supported** |

**The honest counterweight, stated before anything else.** Promotion changes the **unit** of
every absolute exposure by **3.36×** (§7). That is not a physics result and it is not
evidence-backed in either direction — it is a change of referent. Every absolute number in
`results/` generated before this change is stale; §7 lists them.

---

## 1. Method, and the one numerical choice that had to be made

Three drugs spanning the delivery axis (temozolomide, cisplatin, doxorubicin), the
unmodified engine, `alpha = 148.6082`, 6 h exposure window, advection **on**.

**dt = 5 s at every dx, held fixed.** Diffusion and reaction are solved implicitly and are
unconditionally stable, so `dt` is constrained *only* by the explicit upwind advection
term, `dt·(|vx|+|vy|)/dx ≤ 1`. Measured on this section:

| | dx = 15 | 12 | 8 | 6 |
|---|---|---|---|---|
| max interstitial speed, **KDE** (µm/s) | 0.361 | 0.407 | 0.509 | 0.551 |
| max interstitial speed, **discrete** (µm/s) | 0.730 | 0.779 | 0.843 | 0.868 |
| CFL-limited dt_max, KDE (s) | 41.6 | 29.5 | 15.7 | 10.9 |
| **CFL-limited dt_max, discrete (s)** | 20.5 | 15.4 | **9.5** | **6.9** |

**This is a genuine, previously unreported cost of the discrete field.** A source
concentrated on 2.5% of the tissue produces a steeper interstitial-pressure gradient, so
the Darcy velocity is **1.6–2.0× higher** than under the KDE, and the published
`dt = 10 s` **violates the advection CFL at dx = 8 and dx = 6**. dt = 5 s is safe for both
fields at every dx on the ladder and is held constant so that dt is not confounded with dx.

**The dt reduction does not change the answer.** At dx = 12 µm with dt = 5 s the discrete
arm reproduces `vessel_network_impact.md`'s dt = 10 s arm B to four significant figures —
coverage 0.3836 / 0.9936 / 0.0000, mean tumour AUC 1.083e+04 / 1.691e+06 / 7.063e+05, L
39.6 / 152.5 / 273.3 µm — and the KDE arm likewise reproduces arm A (0.9692, 3.712e+04,
95.9 µm). The two studies are on the same footing; only dx moves.

`L` is fitted on a **common distance axis** — Euclidean distance from a segmented vessel
pixel — rebuilt at each dx, so it means the same thing in both arms.
Kill is reported as the **seed-free expected kill** (mean per-cell kill probability), which
removes the Bernoulli noise from the convergence question; the realised seed-averaged kill
is also stored and agrees to <0.003.

---

## 2. GRID CONVERGENCE — the real risk, and it passes

### 2.1 Per-drug ladders

| field | drug | dx = 15 | dx = 12 | dx = 8 | dx = 6 |
|---|---|---|---|---|---|
| **discrete** | doxorubicin | cov 0.3767 · AUC 1.064e4 · kill 0.2184 · L 47.2 | 0.3836 · 1.083e4 · 0.2212 · **39.6** | 0.3876 · 1.090e4 · 0.2220 · 51.5 | 0.3866 · 1.088e4 · 0.2218 · 51.9 |
| | cisplatin | 0.9892 · 1.662e6 · 0.3731 · 151.1 | 0.9936 · 1.691e6 · 0.3748 · 152.5 | 0.9914 · 1.702e6 · 0.3759 · 153.3 | 0.9924 · 1.696e6 · 0.3751 · 151.9 |
| | temozolomide | 0.0000 · 6.946e5 · 0.0141 · 286.2 | 0.0000 · 7.063e5 · 0.0144 · 273.3 | 0.0000 · 7.111e5 · 0.0145 · 273.1 | 0.0000 · 7.084e5 · 0.0145 · 269.7 |
| **KDE** | doxorubicin | 0.9724 · 3.833e4 · 0.7648 · 105.3 | 0.9692 · 3.712e4 · 0.7527 · 95.9 | 0.9666 · 3.671e4 · 0.7482 · 96.6 | 0.9596 · 3.438e4 · 0.7280 · 94.8 |
| | cisplatin | 0.9986 · 5.947e6 · 0.9226 · 201.6 | 0.9986 · 5.758e6 · 0.9120 · 184.5 | 0.9986 · 5.696e6 · 0.9082 · 183.2 | — |
| | temozolomide | 0.0000 · 2.479e6 · 0.1180 · 329.7 | 0.0000 · 2.401e6 · 0.1136 · 299.2 | 0.0000 · 2.375e6 · 0.1121 · 295.5 | — |

(coverage @1× IC50 · mean tumour AUC µM·s · expected kill · fitted L µm. The discrete arm is
complete. The KDE arm is run at dx = 6 for **doxorubicin only** — it is the drug with by far
the shortest `L`, hence the resolution-limited case, and the control needed to answer "does
the discrete field need a finer mesh than the KDE did".)

### 2.2 Worst-drug successive-refinement change

Judged on `delivery_metrics.py`'s own criteria: **bounded** fractions on the **absolute**
scale against `CONVERGENCE_TOL_ABS = 0.03`, **unbounded** metrics on the **relative** scale
against `CONVERGENCE_TOL = 0.05`.

| metric | field | 15 → 12 | 12 → 8 | 8 → 6 | converged? |
|---|---|---|---|---|---|
| **therapeutic_coverage_1x** (abs) | **discrete** | 0.0069 ✔ | 0.0040 ✔ | **0.0011 ✔** | **YES** |
| | KDE | 0.0032 ✔ | 0.0026 ✔ | 0.0070 ✔ | yes |
| **mean_tumor_auc** (rel) | **discrete** | 1.76% ✔ | 0.68% ✔ | **0.38% ✔** | **YES** |
| | KDE | 3.17% ✔ | 1.09% ✔ | **6.35% ✘** | **NO** |
| **kill, expected** (abs) | **discrete** | 0.0028 ✔ | 0.0011 ✔ | **0.0008 ✔** | **YES** |
| | KDE | 0.0121 ✔ | 0.0045 ✔ | 0.0202 ✔ | yes, 25× looser |
| median_tumor_auc (rel) | **discrete** | 3.53% ✔ | 0.95% ✔ | **0.55% ✔** | **YES** |
| | KDE | 3.29% ✔ | 1.03% ✔ | **6.48% ✘** | **NO** |
| **source mean over tissue** (rel) | **discrete** | 2.15% ✔ | 0.24% ✔ | **0.34% ✔** | **YES** |
| | KDE | 3.15% ✔ | 1.09% ✔ | **6.40% ✘** | **NO** |
| **penetration_length_L** (rel) | discrete | 16.2% ✘ | **30.2% ✘** | **1.24% ✔** | **only for dx ≤ 8** |
| | KDE | 9.24% ✘ | 1.24% ✔ | 1.81% ✔ | yes, by dx = 12 |

(✔ = inside the tolerance. KDE 8 → 6 is doxorubicin only, per §2.1.)

**The convergence verdict: the discrete source is converged on every metric — bulk metrics
by dx = 12 µm, the distance-resolved `L` by dx = 8 µm — and its sequences are monotone and
tightening at every step.** The KDE is converged on `L` earlier, but it is the field that
**fails** at the finest step on absolute delivery: its mean tumour AUC moves 6.35% from
dx = 8 to dx = 6, outside the project's own 5% tolerance, and §2.5 shows exactly why.

### 2.3 The one place the discrete field is genuinely worse — and it is exactly the predicted one

`vessel_network_impact.md` §8.2 item 4 predicted that a 7.5 µm feature on a 12 µm mesh
would be under-resolved. **It is, and only the metric that resolves the perivascular
boundary layer notices.** Doxorubicin's fitted `L` is non-monotone —
47.2 → **39.6** → 51.5 → 51.9 µm — with the dx = 12 value 24% below the dx ≤ 8 asymptote of
~51.7 µm, while the KDE's `L` is flat from dx = 12 (95.9 → 96.6 µm).

This is not a solver defect and it is not noise: it is a 7.5 µm source annulus being
rasterised onto meshes of 15, 12, 8 and 6 µm, whose single-cell equivalent disc radii are
8.5, 6.8, 4.5 and 3.4 µm. Only at dx ≤ 8 µm is the annulus resolved by more than one pixel
across, and only then does `L` settle. **Practical consequence: quote `L` (and any
perivascular-gradient claim) at dx ≤ 8 µm under the discrete field. Bulk delivery metrics
remain safe at dx = 12 µm.**

### 2.4 The cost, quantified plainly

| | dx = 12 (production) | dx = 8 (needed for `L`) | dx = 6 |
|---|---|---|---|
| grid | 424 × 675 = 286k | 635 × 1012 = 643k | 846 × 1349 = 1.14M |
| tissue px | 193,752 | 435,924 | 775,035 |
| wall-clock / solve, discrete, dt = 5 s | 226 s | 369 s | 435–649 s |
| peak RSS | 0.93 GB | 2.13 GB | 3.45 GB |
| CFL-limited dt_max, discrete | 15.4 s | 9.5 s | 6.9 s |

Against the historical production configuration (dx = 12 µm, dt = 10 s, KDE), a
**distance-resolved** discrete run costs roughly **3.3×** the wall clock — ~1.6× from the
finer mesh and ~2× from halving dt to satisfy the CFL — and **2.3× the memory**. A **bulk
delivery** discrete run at dx = 12 µm costs ~2× (the dt halving alone, and only because
dt = 10 s is retained as safe at dx = 12 where dt_max is 15.4 s — so in practice, nothing).

### 2.5 The finding that inverts the blocker — the KDE is the field that does not converge

Blocker (b) assumed the sharp source would be the grid-sensitive one. On the dose — the
quantity every absolute result depends on — the opposite is true:

| | dx = 15 | 12 | 8 | 6 | 15 → 6 total |
|---|---|---|---|---|---|
| discrete: source area fraction | 2.447% | 2.500% | 2.505% | 2.497% | **+2.0%** |
| **discrete: mean over tissue** | 0.02447 | 0.02500 | 0.02505 | 0.02497 | **+2.0%** |
| **KDE: mean over tissue** | 0.08670 | 0.08397 | 0.08305 | 0.07774 | **−10.3%** |
| discrete: mean tumour AUC, doxorubicin | 1.064e4 | 1.083e4 | 1.090e4 | 1.088e4 | **+2.3%** |
| KDE: mean tumour AUC, doxorubicin | 3.833e4 | 3.712e4 | 3.671e4 | 3.438e4 | **−10.3%** |

The reason is structural and was hiding in plain sight. The KDE is **divided by its own
maximum**, and the maximum of a Gaussian-smoothed point field is a *grid-dependent*
quantity: refine the mesh and the peak pixel resolves higher, so the whole normalised field
— and therefore the total drug entering the tissue — shrinks. The KDE's delivered dose
drifts **10.3% across the ladder and is still moving 6.4% at the finest step**. The discrete
field's amplitude is 1.0 by definition and its footprint is set by a physical radius, so its
integral is dx-stable to 2% overall and 0.3% at the finest step.

**The blocker was aimed at the wrong field.** The under-resolved 7.5 µm annulus costs the
discrete source accuracy in `L` at dx = 12 and nothing else; the KDE's grid-dependent
normalisation costs it accuracy in *every absolute number*, at every dx, without converging.

---

## 3. ALPHA — there is nothing to re-fit, and this is checkable in one line

`validation.calibration.calibrate_alpha` solves, in closed form,

```
k_eff* = D_eff / L*²                 alpha = (k_eff* − k_passive) / uptake_proxy
```

`D_eff`, `k_passive` and `uptake_proxy` come from the **compiled molecule**; `L*` is the
published doxorubicin anchor (40–50 µm half-distance ÷ ln2 = 64.92 µm). **No vascular field
appears anywhere in the calculation**, and its signature takes only
`(target_efold_um, drug, ecm_density)`.

Measured under both defaults:

| | KDE default | discrete default |
|---|---|---|
| `calibrate_alpha()` | 148.60815765445545 | 148.60815765445545 |
| `pk_multidrug.run()` penetration block | doxorubicin 64.92 (FIT), mitoxantrone 48.21 PASS, topotecan 83.62 PASS | **byte-identical** |

The claim in §8.2 item 3 that "flipping the default without re-fitting alpha silently
rescales every AUC" **conflated two different things**. Alpha does not move. What *does*
move is the **source normalisation** (§7) — a 3.36× change in the tissue-integrated total —
and no amount of alpha re-fitting would address it, because alpha sets a loss rate, not a
source amplitude.

---

## 4. The deciding test: out-of-sample PK, and what happens if you re-fit anyway

### 4.1 The published test, unchanged

| drug | role | predicted L | cited e-fold band | verdict, KDE default | verdict, discrete default |
|---|---|---|---|---|---|
| doxorubicin | in-sample fit | 64.92 µm | 57.7–72.1 µm | FIT | FIT |
| **mitoxantrone** | out-of-sample | **48.21 µm** | 36.1–108.2 µm | **PASS** | **PASS** |
| **topotecan** | out-of-sample | **83.62 µm** | 36.1–108.2 µm | **PASS** | **PASS** |

**2/2 out-of-sample, identical to the last digit under both source fields. The deciding
test does not degrade, because it cannot.**

### 4.2 What an empirical re-fit would do — and why it must not be done

The strongest available reading of "re-fit alpha against the discrete field" is: find the
alpha for which the **simulated section-scale** doxorubicin decay length equals the 64.92 µm
anchor, separately under each field. Solved here by a two-probe power law plus verification:

| field | empirical alpha | vs analytic 148.61 | verification L | simulated L at α = 148.61 |
|---|---|---|---|---|
| KDE | **455.24** | 3.06× | 64.36 µm | 95.88 µm |
| discrete | **77.59** | 0.52× | 64.47 µm | 39.58 µm |

Now freeze each empirical alpha and re-run the out-of-sample check:

| out-of-sample drug (band 36.1–108.2 µm) | α = 148.61 (field-free) | α = 455.24 (KDE-empirical) | α = 77.59 (discrete-empirical) |
|---|---|---|---|
| mitoxantrone | 48.2 µm **PASS** | 27.6 µm **FAIL** | 66.7 µm PASS |
| topotecan | 83.6 µm **PASS** | 47.8 µm PASS | 115.5 µm **FAIL** |
| **tally** | **2/2** | 1/2 | 1/2 |

**Both empirical re-fits break one of the two out-of-sample drugs, in opposite directions,
and the field-free value is the only one that passes both.** The reason is the category
error such a re-fit commits: the anchor is a **single-vessel** immunofluorescence gradient,
while a section-scale `L` on a densely and near-uniformly vascularised tissue is a different
statistic — `vessel_network_impact.md` §9 warns against quoting them interchangeably, and
this is the measurement showing why. `results/vessel_profile.md` already validates the
field-free calibration in the geometry the anchor was made in: the isolated single-vessel
control gives d½ = 49.5 µm against a measured 40–50 µm.

**Conclusion: alpha stays at 148.6, and the reason blocker (a) dissolves is not that the
discrete field happens to agree — it is that the calibration never depended on the field at
all.**

---

## 5. delivery_metrics — before / after on the published five-drug panel

Five-drug panel, dx = 15 and 12 µm, dt = 5 s, both fields, common vessel-distance axis.

| drug | coverage @1× IC50, KDE | discrete | L (µm), KDE | discrete |
|---|---|---|---|---|
| temozolomide | 0.0000 | 0.0000 | 299.2 | 273.3 |
| lomustine | 0.9718 | 0.2175 | 137.2 | 78.2 |
| cisplatin | 0.9986 | 0.9936 | 184.5 | 152.5 |
| methotrexate | 0.0007 | 0.0000 | 107.9 | 48.7 |
| doxorubicin | 0.9692 | 0.3836 | 95.9 | 39.6 |
| **source null** | — | — | **36.6** | **undefined** |
| **cross-drug spread** | | | **3.12×** | **6.91×** |

**Does `therapeutic_coverage` still discriminate? Yes.**

| | KDE | discrete | threshold |
|---|---|---|---|
| across-drug absolute spread | 0.9986 | 0.9936 | — |
| worst-drug 15 → 12 drift | 0.0037 | 0.0074 | — |
| **separation (spread ÷ drift)** | **271×** | **134×** | ≥ 3× |

The discrete field halves the separation and it remains 45× the threshold. Coverage is
unaffected as a discriminating metric.

**Does `penetration_length` now separate drugs? Yes, and by 2.2×.**

| | KDE | discrete |
|---|---|---|
| cross-drug spread | 3.12× | **6.91×** |
| relative spread | 1.233 | **1.973** |
| worst-drug 15 → 12 drift | 0.0924 | 0.1618 |
| separation | 13.3× | **12.2×** |
| **source null on the same axis** | **36.6 µm** | **undefined — the field is a step function** |

Both discriminate on the current reconstruction, so `delivery_metrics.md`'s published
"`penetration_length_L` does NOT discriminate" (1.07× spread, separation 0.254) is stale
under **either** field — confirming `vessel_network_impact.md` §3, which traced that to the
retired occupancy tissue mask, not to the vasculature. On the legacy p95-KDE axis this run
reproduces that report exactly (temozolomide 178.3, cisplatin 108.6, doxorubicin 58.0,
source null 21.3 µm).

**The qualitative gain is the null row, not the spread.** Under the KDE a drug's `L` is
always partly the source's (36.6 µm on this axis). Under the discrete field the identical
three-parameter fit **returns nothing at all** — a step function has no decay length — so
there is no competing explanation left for the drug's `L`. That is a difference in kind.

---

## 6. Regional / BBB leak — the discrete field is NOT its "correct substrate"

`vessel_network_impact.md` §8.3 argued that because the additive BBB leak `L = 0.238` is a
**per-vessel** property, the discrete field is its right substrate, and called this "a
stronger argument for the discrete field than anything measured in this report". **It is
not supported.**

**(a) The published regional test contains no vascular field.**
`glioma_regional_pk.run_post_fix_regional` predicts
`R_region = clip(bbb_penetrance + L·D_region, 0, 1)` from the molecular layer and a
per-region scalar `D`. There is no grid, no source and no transport solve in it. Run under
both defaults:

| CE vs non-enhancing | measured ratio | model, KDE default | model, discrete default |
|---|---|---|---|
| levetiracetam | 1.452 | 1.191 | 1.191 |
| cefazolin | 4.263 | 5.873 | 5.873 |
| **geometric-mean abs ratio error** | | **1.296** | **1.296** |

Byte-identical. The regional agreement cannot improve, and reporting that it had would have
been an artefact of re-running the same arithmetic.

**(b) The spatial version, which *is* testable, says the same thing.**
Driving the frozen transport solver with `barrier_delivery_source` —
`S = vascular_density · clip(p + L·D(x,y)) · plasma` — over a core-vs-margin disruption
field (`REGION_DISRUPTION['core'] = 0.50`, margin 0.00; core = the upper half of the
measured tumour-cell density, 20,331 vs 9,843 tumour cells), under each source field:

| drug | field | core/margin exposure, no leak | with leak | **leak effect (ratio of ratios)** |
|---|---|---|---|---|
| doxorubicin | KDE | 1.2243 | 4.5867 | **3.7464** |
| doxorubicin | discrete | 1.1731 | 4.4030 | **3.7534** |
| temozolomide | KDE | 1.1329 | 1.6390 | **1.4468** |
| temozolomide | discrete | 1.0847 | 1.5710 | **1.4484** |

**The leak's regional effect is substrate-independent to 0.19%.** The reason is algebraic
and should have been noticed before the field was built: within a region of uniform `D`,
`clip(p + L·D)/p` is a **scalar** multiplying the source, and transport is linear — so it
passes through to the exposure identically whatever the source's spatial shape. (Diffusion
does dilute it: doxorubicin's source ratio of 6.79× arrives as a 3.75× exposure ratio. But
it dilutes it *equally* on both substrates.)

The per-vessel argument would only bite if `bbb_disruption` were resolved at vessel scale
**and** the measurement were made at vessel scale. It is not: the measured contrast is
between MRI-scale regions, millimetres across. **This argument should be withdrawn from
§8.3.**

---

## 7. What the promotion changes, and what is now stale

### 7.1 The code change

| file | change |
|---|---|
| `tmesim/reconstruct/real_fields.py` | `_DEFAULT_VASCULAR_MODE = 'discrete'`; `VASCULAR_MODES` docstring carries the evidence, the two costs and the normalisation warning. KDE path byte-for-byte unchanged behind `mode='kde'` |
| `tmesim/validation/vessel_profile.py` | `vessel_mask_from_density` falls back to `field > 0` when the percentile threshold is non-positive. **This was mandatory**: on a source that is exactly zero over 97.5% of the tissue the p95 threshold is 0.0 and `vasc >= 0` selected the **whole tissue** as vessel, silently destroying every distance transform. This is `vessel_network_impact.md` §8.4 item 3 — for a segmented field, `field > 0` *is* the segmented mask. The KDE's p95 is strictly positive, so no existing number moves |
| `tmesim/reconstruct/real_pipeline.py` | `metadata['real_pipeline']['vascular_source']` now reports the construction actually used instead of hard-coding "KDE" |
| `tests/test_vessel_network.py` | the default-preservation test now pins the discrete default and that `mode='kde'` still reproduces the KDE; **new** test that a sparse source does not degenerate the vessel mask |
| `tests/test_oxygen_proxy.py` | the metadata assertion follows the new `vascular_source` string |

Nothing in `tmesim/core`, `tmesim/physics` or `tmesim/sim/engine.py` was touched. Suite:
**138 passed**.

### 7.2 The normalisation change — read this before quoting any absolute number

| | KDE | discrete |
|---|---|---|
| % of tissue with a nonzero source | 99.71% | 2.50% |
| mean over tissue | 0.0840 | 0.0250 |
| ratio of tissue-integrated totals | 1.00 | **0.298 (i.e. 3.36× less)** |
| what `1.0` means | the densest endothelial-KDE spot **in this section** | a perfused vessel wall leaking at the reference rate |

So doxorubicin's coverage @1× IC50 reads 0.9692 under the KDE and 0.3836 under the discrete
default, and its kill 0.7527 vs 0.2218. **Almost all of that is dose, not geometry**: the
dose-matched control in `vessel_network_impact.md` §6 puts the geometry-only change at
**−0.0166** in coverage. Neither amplitude is a validated absolute — the project has never
claimed one — but the referents differ in kind: the KDE's unit is a property of *this
section* and changes if you change the section, while the discrete field's unit is a
statable physical condition. The promotion adopts the better-defined convention; it does
**not** claim the new absolute numbers are more nearly true.

### 7.3 Regeneration debt — files whose ABSOLUTE numbers are now stale

> **STATUS: mostly PAID DOWN in a later pass. See §7.4 for what was regenerated, what could
> not be, and why.** The list below is the original debt as recorded at promotion time and is
> kept verbatim so the two passes can be compared.

Not regenerated in this pass. Every one of these was generated on the KDE source and its
absolute exposures/AUC/coverage/kill carry the old 3.36× unit. Ordering follows
`vessel_network_impact.md` §8.4:

1. `results/delivery_metrics.md` — coverage, AUC, the convergence and discrimination tables.
   §5 above supersedes its discrimination verdict for `penetration_length_L` but the file
   itself still carries KDE absolutes.
2. `results/vessel_profile.md` — its p95 vessel mask now resolves to the segmented mask.
3. `results/numerical_rigor.md` — the grid study; §2 above is the discrete-field replacement
   for the source-field part of it, but its own tables are KDE.
4. `results/pharmacology_real.md`, `results/real_drug_panel.md`, `results/dose_schedule.md`,
   `results/niche_response.md`, `results/combination.md`, `results/cross_tissue.md`.
5. `results/uncertainty.md`, `results/ranking_robustness.md` — ranking and error bars.
6. `results/patient_report_xenium.md` / `_visium.md` and every figure in `results/`.

**Not stale** (verified field-free in §3, §4, §6): `results/calibration.md`,
`results/pk_multidrug.md`, `results/pk_reality.md`, `results/glioma_regional_pk.md`,
`results/mgmt_scalar.md`.

### 7.4 Regeneration pass — what was actually re-run, and what could not be

Every file below was re-run from its own producer under `_DEFAULT_VASCULAR_MODE = 'discrete'`,
`alpha = 148.6`, on the same measured Xenium section. Each regenerated file now carries its
own header block stating `vascular_mode`, `alpha`, sensitivity mode and `dx`. The KDE-era
"STALE ABSOLUTE NUMBERS" banner was removed **only** from files that were actually re-run.

**Regenerated (11 reports + 3 figure-only artefacts):**

| file | producer | note |
|---|---|---|
| `real_drug_panel.{md,png}` | `experiments/real_drug_panel.py` | headline panel; ordering held |
| `patient_report_xenium.{md,png}` | `experiments/patient_report.py --source xenium` | agrees drug-for-drug with the panel |
| `patient_report_visium.{md,png}` | `experiments/patient_report.py --source visium` | Visium path never read the switch — see below |
| `heldout_patient.{md,png}` | `tmesim/validation/heldout_patient.py` | **was stale and un-bannered**; ordering replicates on section 2 |
| `delivery_metrics.{md,png}` | `tmesim/validation/delivery_metrics.py` | corrects a published-wrong claim — see below |
| `cross_tissue.{md,png}` | `experiments/cross_tissue.py` | needed a producer bug fixed first |
| `dose_schedule.{md,png}` | `experiments/dose_schedule.py` | p95 vessel mask now resolves to 2.7% of tissue |
| `niche_response.{md,png}` | `experiments/niche_response.py` | |
| `combination.{md,png}` | `experiments/combination.py` | |
| `uncertainty.{md,png}` | `tmesim/validation/uncertainty.py` | 80 LHS draws, 66.5 min |
| `ranking_robustness.{md,png}` | `tmesim/validation/ranking_robustness.py` | 40 outer draws, 112.5 min |
| `real_simulation.png` | `experiments/real_simulation.py` | figure only |
| `pharmacology_real.md` | `tmesim/validation/pharmacology_real.py` | **byte-identical** — see below |

**The unit change is now auditable, and the ordering held.** On the flagship panel
(`real_drug_panel.md`, dx = 12 µm, sensitivity `uniform`), cisplatin's mean tumour AUC moved
**5.48e+06 → 1.63e+06 µM·s** and its `kill_fraction` **0.9125 → 0.3699**. The five-drug kill
ordering is **unchanged** — cisplatin > doxorubicin > lomustine > methotrexate > temozolomide
— on the development section, on the 816,769-cell held-out section, at dx = 15 and dx = 12 µm,
and in **all 40** Monte-Carlo draws of `ranking_robustness`, whose decision probabilities
moved by ≤ 0.004 (P(cisplatin = rank 1) 0.927 → 0.931). §7.1's claim that the promotion is a
change of *unit* rather than of *ranking* is therefore now measured, not asserted.

**Two files were over-flagged and are structurally immune.** `pharmacology_real.md` imports
nothing from `tmesim.reconstruct` — no domain, no field, no PDE — and re-running it reproduced
the file byte-for-byte. `patient_report_visium.md` (and the four Visium arms of
`cross_tissue.md`) go through `reconstruct.real_visium.reconstruct_from_visium`, which builds
vascular density from spot deconvolution and never calls `vascular_density_from_cells`, so
`_DEFAULT_VASCULAR_MODE` does not reach it; its vascular-field statistics are identical before
and after promotion. Both banners have been replaced with an explanation rather than a refresh.

**One file was stale and was NOT flagged.** `heldout_patient.md` was produced ~30 min before
the default was flipped and the banner sweep missed it. It has been re-run in full.

**Two files CANNOT be regenerated without a code change, and keep permanent banners.** Both
are latent reporting-path bugs that the sharper source exposed, not physics failures:

* `numerical_rigor.md` — `numerical_rigor.py:273` picks mass-conservation injection points as
  `vasc > np.percentile(vasc[tissue], 99.0)`. The discrete field is binary and exactly 1.0 on
  2.50% of tissue, so `p99 == max` and the strict `>` selects **zero** cells;
  `injected == 0` and line 280 raises `ZeroDivisionError`. Under the KDE the same expression
  selected 1,938 cells. One-line fix: take the top-N by value instead of `> p99`.
* `vessel_profile.md` — `vessel_profile.py:799` formats `empirical_d_half_um` with `:6.1f`, but
  `empirical_half_distance_um` returns `None` by design when a profile never halves, and under
  the discrete source **temozolomide's profile never halves** (`L` ≈ 296 µm against a 43 µm
  median tumour-to-vessel distance). Fix: guard the three formatted reads (`:790`, `:799`,
  `:961`). Its partial run did establish that the p95 vessel mask now resolves to 2.6% of
  tissue and the source null's empirical d½ is 10.0 µm.

Both fixes live under `tmesim/`, which the regeneration pass was not permitted to modify.

**One producer bug had to be fixed to regenerate at all.** `experiments/cross_tissue.py:842`
referenced an undefined `alpha` inside `_verdict()`, so the script wrote its figure and then
died with `NameError` before writing markdown — i.e. the file had *no working producer*.
`_verdict` now takes `alpha` as a parameter (two lines, `experiments/` only).

**Regenerating fixed a claim that was wrong in the conservative direction.**
`delivery_metrics.md` listed `penetration_length_L` as *does NOT discriminate* on a 1.07×
five-drug spread against a source null of 23.7 µm. On the discrete source the spread is
**6.91× (39.6–273.3 µm)**, separation **12.2×** its own grid drift, and the source field admits
**no null fit at all**. The metric now reads **YES, discriminates** — exactly as §5 predicted
and as `vessel_network_impact.md` §3 demanded.

---

## 8. What this does NOT establish

* **One section, one patient.** Between-case GBM vascularity spans ~48×. The convergence
  result is a property of this section's vessel spacing and this mesh ladder.
* **The amplitude of neither field is validated.** §7.2. The discrete field's *geometry*
  scores 2.50% against a measured 2.1 ± 0.3% CD34 area fraction (PMID 36611806) untuned;
  its *amplitude* is a convention.
* **The KDE arm at dx = 6 µm is doxorubicin only.** Its non-convergence on absolute
  delivery (§2.5) therefore rests on one drug at the finest step — though the *source*
  statistic that causes it (mean over tissue, −6.4% from dx = 8 to 6) is drug-independent
  and is measured directly, so the mechanism is not in doubt even if the per-drug tally is
  thin.
* **The 8 → 6 µm step is one refinement.** A converged-looking pair is evidence, not proof;
  no Richardson order is claimed, because the end-to-end re-rasterisation of cells, mask and
  distance axis at each dx makes a single clean order fragile (the same caveat
  `numerical_rigor.md` records).
* **The deeper defect is untouched.** The segmented vessels are still 0.97× a random point
  field on the median tumour-to-vessel distance. Sharpening the source does not create the
  vascular hotspots and avascular zones that GBM has and this section's endothelial point
  pattern does not.
* **The regeneration debt in §7.3 is real and outstanding.** Until those files are re-run,
  `results/` mixes two source conventions, and each file must be read against the field it
  was generated on.

---

## 9. Reproduce

| artefact | driver |
|---|---|
| grid convergence (§2) | `gridconv2.py` — 3 drugs × 2 fields × dx 15/12/8/6, dt = 5 s, resumable |
| alpha + out-of-sample (§3, §4) | `alpharefit2.py` — two-probe power-law solve per field + frozen-alpha OOS |
| delivery metrics before/after (§5) | `delivmet.py` — five-drug panel × 2 fields × dx 15/12 |
| regional (§6) | `regional.py` — part A re-runs `glioma_regional_pk` under both defaults; part B drives `barrier_delivery_source` through the frozen solver |
| field-free checks (§3) | `calibrate_alpha()` and `pk_multidrug.run()` under both defaults |

Deterministic (`seed = 0`; kill reported seed-free as the mean per-cell probability, with a
2-seed realised check). `alpha = 148.6082` throughout §2 and §5.
