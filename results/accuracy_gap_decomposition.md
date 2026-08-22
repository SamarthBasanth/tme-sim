# Localising the "~14× too shallow" gradient gap — what it is a ratio of, and where the error actually lives

**Status: MEASURED, not argued.** Driver:
`experiments/accuracy_gap_decomposition.py`. Data:
`results/accuracy_gap_decomposition_data.json`. Target:
`results/kyle2004_penetration.md` §6.1 / C4.

**Nothing under `tmesim/`, `app/`, `docs/`, `landing-3d/`, `scripts/` or
`data_cache/` was modified. No default, constant or shipped parameter was
changed. `calibrate_alpha()` is asserted live at `22.223477711172563` in every
stage and holds. Nothing was fitted to any measurement.** Where a parameter is
scaled it is scaled on a **local `DrugParams` copy**, over a grid declared in the
script's constants before the run, the whole curve is reported, and no value is
adopted, written back or recommended.

---

## 0. Verdict up front

**The gap is real, it is in the optimistic direction, and it is about seven
times smaller than the way it has been quoted. It is also, on the product's own
readout, about five times smaller again.**

1. **"~14× too shallow" is a ratio of two DIMENSIONLESS ratios, not of two
   lengths.** It is Kyle 2004's measured near:far exposure ratio for
   daunorubicin (50) over the engine's modelled near:far exposure ratio
   (3.6457) = **13.715**. Both numbers are concentration ratios. Neither is a
   length.

2. **So the arithmetic that has been done on it is wrong.** "L goes like
   `sqrt(D_eff/k)`, so 14× in L is ~196× in `D_eff/k`" is a category error that
   overstates the problem by roughly 30×. Kyle's readout is **exponential** in
   the transport length. Run both sides through the same estimator and the
   length is **2.45× too long** (band **2.00–2.98×** once the paper's own SDs are
   propagated), i.e. `D_eff/k` is **6.0× too large**.

3. **Two thirds of even that is not a prediction failure.** Log-additively:
   **1.832× (67.5%) is inherited from doxorubicin**, the compound
   `calibrate_alpha()` is fit to, and **1.339× (32.5%) is daunorubicin-specific**.
   The doxorubicin term is a **conflict between two external anchors** and it is
   now formally hard: **Primeau gives `λ_dox` 57.71–72.13 µm, Kyle gives
   30.73–55.68 µm, and those bands are DISJOINT.** §5.

4. **It is a `k` error, not a `D_eff` error — and there is a landmine here.**
   Two anchors already in this repository, never combined, separate them:
   tirapazamine's `k` is **7.5–12.5× too small** and its `D_eff` is **1.5–2.6×
   too small as well, which pushes `λ` the OPPOSITE way and is currently masking
   part of the `k` error**. On the **MCL-appropriate arm** — the one this repo
   itself names for this preparation — `D_eff` lands **inside** its measured band
   and **`k` carries 100 % of a 2.10–3.61× length error**. **Correcting `D_eff`
   alone to its own scorecard row makes the perivascular gradient WORSE.** §6.

5. **The exponential amplification does NOT transfer to the product's own
   readout.** Measured on the real 10x Xenium GBM section: the map's
   half-distance responds to `λ` with an exponent of **0.92–1.22**, essentially
   linear, against **2.92** for Kyle's near:far ratio. **A 2.45× length error
   shows up on the map as ~2.2–3.0×, not 14×.** §7.

6. **The preparation-mismatch hypothesis was tested with the correction declared
   in advance, and it FAILED.** Setting arm A's ECM to what this repository
   declares an MCL to be makes the gap **worse**, 13.71× → **14.41×**. §4.2.

7. **And the largest single threat is that the anchor may not measure a length at
   all.** Measured: pure retarded diffusion with **zero loss** in Kyle's geometry
   reproduces near:far **13.6 at a retardation factor of 50 and 84.1 at 90** — a
   factor around 60–70 reproduces the measured 50 **with `λ = ∞`**. Kyle's own
   ~7× near-side accumulation for daunorubicin is direct evidence of exactly
   that store, and §9.2 of the Kyle file says the engine does not have one. §8.2.

**What to say instead of "the map is 14× too shallow":** *on the map's own
readout the penetration depth is about 2–3× too long; the length scale behind it
is ~2.5× too long; most of that is an unresolved conflict between two external
anchors rather than a prediction error; the term that is wrong is `k`, not
`D_eff`; and it is not yet established that the anchor measures a length scale
at all.*

---

## 1. Step 0 — is the 14× what it claims to be?

### 1.1 Name both quantities and their units

| | quantity | units |
|---|---|---|
| numerator | Kyle 2004 Table 1, drug-derived-fluorescence column: near-side signal over the first 30 µm divided by far-side signal over the last 30 µm of a **140 ± 15 µm** HCT-116 multicell layer, 1 h one-sided exposure at 10 µmol/L, far side sealed, frozen immediately after exposure | **dimensionless** |
| denominator | the engine's arm-A readout: mean 1 h **AUC of free interstitial** drug over the first 30 µm divided by mean AUC over the last 30 µm of a uniform 140 µm slab at **`dx = 2 µm`, `ecm = 0.0875`, planar source in column 0, no-flux far face, convection OFF**, through `simulate_response` | **dimensionless** |

**14× = 50 / 3.6457 = 13.715**, band from the measurement's own ±1 SD
**9.60–17.83**. A ratio of ratios.

**The related number, so the two are never confused.** §6's other headline,
"~3.7× compressed" (C3), is `(50/11)/(3.6457/2.9976) = 4.545/1.216 = ` **3.737**
— the *fold between two compounds*. **It is contained inside the 13.715, and the
factorisation is exact:**

```
13.7149  =  3.6696             x  3.7374
   14x      the DOXORUBICIN row     the compound under-separation
            residual, 11/2.9976     C3's "~3.7x compressed"
            -- alpha is FIT to
            this compound, so it
            is an ANCHOR CONFLICT
            (sec. 5), not a prediction
```

**C4 and C3 are not two independent failures. C4 is C3 multiplied by an anchor
conflict, and the anchor conflict is very slightly the larger of the two.**

### 1.2 Is the preparation comparable?

**Yes, and the concern raised about it is misdirected.** The 14× comes from
**arm A**, a uniform 140 µm slab with a one-sided source and a sealed far face,
which reproduces Kyle's apparatus deliberately. **The 10x Xenium GBM tissue
section plays no part in the 14×.** Arm B (the real section) produces the
*ordering* result and the C2 noise-floor result, not this one.

Two **inputs** to arm A are imported from an unrelated preparation. That is a
real defect and it is tested with a declared correction in §4.2 — where it makes
the gap worse.

### 1.3 Reproduce before arguing

Live, at the exact published configuration — **`dx = 2 µm`, `ecm = 0.0875`,
140 µm thick, 30 µm windows, 1 h, convection OFF, through `simulate_response`**:

| compound | live `R_auc` | §6.1 published | relative difference |
|---|---:|---:|---:|
| doxorubicin | 2.997645248639645 | 2.997645248639645 | **0.000e+00** |
| daunorubicin | 3.6456683362306537 | 3.6456683362306537 | **0.000e+00** |

Bit-exact. Everything below stands on it.

### 1.4 Is the model-side number computed the same way as the experiment-side number?

| mismatch | model | experiment | measured effect on `R` |
|---|---|---|---:|
| **estimator** | 1 h **AUC** (time-integrated) | a **snapshot** of tissue frozen immediately after exposure | model's AUC ratio is **1.042× larger** than its own snapshot ratio (daunorubicin) |
| **boundary** | volumetric source in column 0 | stirred reservoir at fixed concentration (Dirichlet) | model's snapshot ratio is **1.011× larger** than the Dirichlet steady state at the same `λ` |
| **measurand** | **free interstitial** drug | **total tissue** drug (fluorescence) | not measurable here — §8.3, signed but unquantified |

**`results/kyle2004_penetration.md` §6 states the boundary mismatch "biases the
modelled ratio downward relative to a clean fixed-concentration boundary". That
is refuted. Measured, it biases it UPWARD by 1.1 %.** §6 uses that sign to argue
that part of C4 is construction. It is not, and the small construction term
points the other way.

`results/spatial_noise_decomposition.md` found the distance **bin width** moves
the real-section half-distance by 18.37 µm. **On arm A there is no binning at
all** — two fixed 30 µm windows on a 1-D slab — so that hazard does not reach the
14×. The arm-A estimator sensitivities that do exist are in §8.1 and none is
capable of touching a factor of 14.

---

## 2. The transfer function, and why 14 is not 2.45

Kyle's geometry at steady state, with a first-order sink, a fixed-concentration
exposed face and a sealed far face, has `c(x) ∝ cosh((T−x)/λ)`, so

```
R(λ) = [ sinh(T/λ) − sinh((T−W)/λ) ] / sinh(W/λ)          T = 140 µm, W = 30 µm
```

**exponential in `T/λ`**. That single fact is the whole of §0.2.

**This is not an approximation of the engine — it IS the engine.** A diagnostic
`k`-scaling on a local `DrugParams` copy over the pre-declared grid 0.25×–32×
traces daunorubicin's `λ` from 128.81 µm to 11.39 µm, at the arm-A configuration
(`dx = 2 µm`, `ecm = 0.0875`, 140 µm, 30 µm windows, 1 h, convection OFF):

| `k` scale | `λ` (µm) | engine `R_auc` | engine `R_snapshot` | analytic `R(λ)` | snapshot / analytic |
|---:|---:|---:|---:|---:|---:|
| 0.25 | 128.81 | 1.623 | 1.508 | 1.499 | 1.0061 |
| 0.50 | 91.08 | 2.208 | 2.087 | 2.071 | 1.0078 |
| **1.00** | **64.40** | **3.646** | **3.499** | **3.460** | **1.0112** |
| 2.00 | 45.54 | 7.736 | 7.515 | 7.407 | 1.0145 |
| 4.00 | 32.20 | 22.667 | 22.215 | 21.850 | 1.0167 |
| 8.00 | 22.77 | 101.872 | 100.469 | 98.846 | 1.0164 |
| 16.00 | 16.10 | 819.270 | 811.480 | 802.380 | 1.0113 |
| 32.00 | 11.39 | 14724.447 | 14627.511 | 14653.061 | 0.9983 |

**The engine's snapshot ratio tracks the analytic transfer function to within
1.7 % across an 11-fold range in `λ`.** Inverting the measurement through `R(λ)`
is therefore legitimate, and it is the only way to put the two sides on one axis.

**Amplification.** `dlnR/dlnλ` is **−1.81** at the engine's `λ`, **−4.46** at
Kyle's, averaging **2.92** over the interval. `13.715 = 2.45^2.92`. **That is
where the 14 comes from, and it is a property of the readout, not of the model.**

### 2.1 The two sides on one axis

| | `R` (near:far) | implied `λ` (µm) |
|---|---:|---:|
| engine, daunorubicin, arm A | 3.6457 | **64.40** |
| **Kyle, daunorubicin** | **50 (35–65)** | **26.26 (24.80–28.54)** |
| engine, doxorubicin, arm A | 2.9976 | **72.49** |
| **Kyle, doxorubicin** | **11 (6–16)** | **39.56 (35.18–49.53)** |

**Daunorubicin: `λ` is 2.453× too long; `D_eff/k` is 6.02× too large.**

---

## 3. Ruled out, and ruled out by construction

**Source geometry and the vessel map contribute NOTHING to the 14×.** Arm A has
no vessel map, no `delivery_source` geometry and no reconstruction — it is a
uniform slab with a planar source occupying column 0. The candidate is excluded
before any measurement.

Measured anyway, because "excluded by construction" is how vacuous checks get
written: sweeping the source density 0.25 → 1.0 at the arm-A configuration moves
`R_auc` by **exactly 0.000e+00**. The readout is exactly invariant to source
level — an independent confirmation, on a different knob, of the sibling result
that the field is homogeneous of degree 1 in the entry term.

---

## 4. The decomposition

### 4.1 The ledger

Daunorubicin, arm A at `dx = 2 µm`, `ecm = 0.0875`, 140 µm, 30 µm windows, 1 h,
convection OFF. Multiplicative, exact to four decimal places:

| step | factor | running value |
|---|---:|---:|
| engine `R_auc`, as published | — | **3.6457** |
| AUC estimator → snapshot estimator | ×0.9597 | 3.4988 |
| source-column boundary → Dirichlet boundary | ×0.9890 | 3.4603 |
| **the `λ` error** | **×14.4495** | **50.0000** = measured |
| **product** | **13.7149** | = the "14×" |

**Read the first two rows carefully. They are worth 0.949× combined and they run
the WRONG WAY: they are not part of the gap, they conceal 5 % of it.** Correcting
the model's estimator to match Kyle's makes the discrepancy 14.29×, not smaller.

**Everything else — 14.4495 of the 13.7149 — is one term: `λ`.**

### 4.2 The preparation term, tested with the correction DECLARED IN ADVANCE — and it failed

Arm A gives an HCT-116 multicell layer the **GBM Xenium section's** mean ECM
(0.0875) and a **fixed** cell volume fraction `φ = 0.7892`. Both are wrong
inputs for this preparation. Promoting them into the ledger requires knowing what
the right inputs are, so both were **declared with their source before the run**,
and a prediction was registered.

**Declared: `ecm` for an HCT-116 MCC = 0.0.** Not a free choice.
`results/cellular_exclusion.md` §4.1: *"an MCL is packed cells with little matrix
and the repo's ECM legs describe matrix the MCL does not have"*, and
`scorecard.py`'s own arm-B note: *"the IMPERMEABLE-cell term alone, the
MCL-appropriate arm, since an MCL is packed cells with little matrix"*. This is
the repository's declared position about this preparation, applied.

**Registered prediction, written before the run:** *correcting `ecm` 0.0875 → 0.0
raises `D_eff` by 1.0554, raises `λ` by 1.0273, and therefore LOWERS `R`. The gap
gets WORSE by about 1.05×.*

**Result — the prediction held:**

| compound | shipped `ecm` 0.0875 | **corrected `ecm` 0.0** | fold on `R` | gap |
|---|---|---|---:|---|
| daunorubicin | `λ` 64.40 µm, `R` 3.646 | `λ` 66.16 µm, `R` 3.469 | **0.9517×** | **13.71× → 14.41×** |
| doxorubicin | `λ` 72.49 µm, `R` 2.998 | `λ` 74.47 µm, `R` 2.870 | **0.9574×** | 3.67× → 3.83× |

**Simulating the preparation Kyle actually used makes the gap WORSE, not better.
The ECM half of the preparation-mismatch hypothesis is refuted.**

The 1.24× that has been quoted for "`ecm` 0.4" is **not a preparation
correction** — `ecm = 0.4` is `calibrate_alpha()`'s `ecm_reference`, which is
*further* from an MCC than 0.0875 is, not closer. That number is a **separate
defect**: `results/kyle2004_penetration.md` §4's `L` column (dox 64.92, epi
64.92, dau 57.68, mit 48.33 µm) is computed at `ecm = 0.4` while arm A
immediately below it is solved at `ecm = 0.0875`, where the same compounds have
`L` = **72.49, 72.49, 64.40, 53.96 µm**. The file prints an `L` next to a run
that does not have that `L`, with no ECM stated in either place. Same defect
class as `results/codon_claim_surface_drift` — correct arithmetic, wrong object.

**`φ` is NOT corrected, and the reason matters more than the number.**
`results/cellular_exclusion.md` §4.1, verbatim: *"A literature search for a
measured cell (or extracellular) volume fraction of a Wilson/Hicks MCL was run
and **returned nothing usable**; that search did not complete cleanly and is
recorded as an open item rather than as a null."* **There is no value to correct
to, and none is invented here.** Only the *direction* is sourced, from the same
paragraph: *"If an MCL is denser than this section — which is the likely
direction, since an MCL has no vasculature, no necrosis and little stroma"*.

Reported as an open-ended signed sensitivity at the corrected `ecm = 0.0`, arm-A
configuration:

| `φ` | `λ` (µm) | analytic `R` | residual gap |
|---:|---:|---:|---:|
| 0.70 | 80.22 | 2.438 | 20.51× |
| **0.7892 (shipped)** | **66.16** | **3.298** | **15.16×** |
| 0.85 | 55.21 | 4.688 | 10.66× |
| 0.90 | 44.69 | 7.785 | 6.42× |
| **0.9647** | — | **50** | **1.00×** |

**That last row is the point, and it is a warning, not a result.** An **unmeasured**
parameter can absorb the **entire** residual. That is exactly why it must not be
allowed to absorb any of it until it is measured — moving `φ` to 0.9647 would be
tuning to the anchor, and it is recorded here only to show how cheap that
temptation is. It is not adopted, not written back, and must never be quoted as a
finding. **`φ` is the largest un-anchored lever in the model, and obtaining an
MCL extracellular-space measurement — which `results/cellular_exclusion.md` §7
already lists — is the only legitimate way to use it.**

Note the `φ` sweep holds `k` fixed, so it is a lower bound: denser packing should
raise `k` too. Note also that the two-route hindrance is at the impermeable
Maxwell limit (0.1512) for all four of these compounds (§5 of the Kyle file), so
`φ` is the *only* lever inside that term for this compound class.

### 4.3 Splitting `λ`

Log-additive in `λ`, so the terms are exact:

| term | fold in `λ` | share | what it is |
|---|---:|---:|---|
| **doxorubicin-inherited** | **1.832×** | **67.5 %** | the residual on the compound `calibrate_alpha()` is fit to. **Not a prediction — this term IS the Primeau-vs-Kyle anchor conflict** (§5). |
| **daunorubicin-specific** | **1.339×** | **32.5 %** | the genuine out-of-sample error. This is §6's "~3.7× compressed" (C3) expressed as a length. |
| total | 2.453× | 100 % | |

On Kyle's own axis, so the sizes are legible:

| | `λ` (µm) | `R` |
|---|---:|---:|
| engine as shipped | 64.40 | 3.460 |
| doxorubicin-inherited term removed | 35.15 | 16.041 |
| daunorubicin-specific term removed | 48.10 | 6.446 |
| both removed (= Kyle) | 26.26 | 50.001 |

**C4 is not an independent failure from C3. It is C3 sitting on top of an anchor
conflict, and the anchor conflict is the larger half.**

---

## 5. The two anchors cannot both be satisfied — a formal result

`calibrate_alpha()` fits α so that doxorubicin's `L` equals
`DOX_TARGET_EFOLD_UM` = 45 / ln2 = 64.92 µm at `ecm_reference = 0.4` — the
Primeau half-distance anchor converted to an e-folding length. Doxorubicin's `λ`
is therefore **pinned by construction** and cannot fail against Primeau;
`results/kyle2004_penetration.md` §2.7 correctly declines to score it.

It can fail against Kyle, and it does:

| anchor | doxorubicin `L` / `λ` |
|---|---|
| **Primeau 2005** (PMID 16361566), 40–50 µm half-distance, converted by the repo's own `half_distance_to_efold_length` | **57.71 – 72.13 µm** |
| **Kyle 2004** (PMID 15342419), near:far 11 ± 5, inverted through Kyle's own geometry with the paper's 140 ± 15 µm thickness SD propagated | **30.73 – 55.68 µm** |

**DISJOINT.** This is not a calibration being imperfect. **There is no value of α
— no value of any single parameter — that satisfies both anchors, because the
model gives each compound one `λ` and the two anchors demand two.** So either one
anchor is wrong, or **they measure different things**, which is live: Primeau is
an in-vivo vascularised xenograft with a decaying plasma input; Kyle is a 1 h
constant-concentration one-sided exposure of a packed multicell layer with a
sealed far face.

**Consequences, and they are permanent:**

* **"Fit α better" is not a fix.** It is the one move that provably cannot work.
  Any α trades Primeau agreement for Kyle agreement one-for-one.
* **67.5 % of the daunorubicin gap is this**, so at most 32.5 % of it is
  reachable by improving the transport model at all.
* `results/two_route_promotion.md` §8.3 records that *"the two external anchors
  now disagree about where the engine sits"* and calls it *"the sharpest open
  question the promotion leaves"*. **This is that disagreement, quantified, and
  it is now known to be formally unsatisfiable rather than merely
  uncomfortable.**
* The honest product statement: *absolute penetration depth is calibrated to one
  anchor; a second anchor in a different preparation disagrees by ~1.8× and the
  two cannot both be satisfied by this model.*

---

## 6. Is it `D_eff` or is it `k`? — and a landmine for the next reader

`λ = sqrt(D_eff/k)` **cannot** separate them. Neither can Kyle's data, nor any
half-distance, near:far ratio or penetration depth, because all of them are
functions of `λ` alone. The engine's own per-compound decomposition (§5 of the
Kyle file: `D` fold 1.0124 the wrong way, `k` fold 1.2826 the right way) is a
statement about the model's internals, not a measurement.

**One thing in this repository separates them, and it had never been assembled.**
Tirapazamine's `D_eff` and its `k` were **each measured**, in the same
multicellular-layer slab-flux experiment, fitted to a diffusion-reaction equation
of the same form as the engine's. Both are already separate rows in
`tmesim/validation/scorecard.py`. Compiled live at `ecm_reference = 0.4` with
α = 22.223477711172563:

| | engine | measured | fold |
|---|---:|---|---:|
| `D_eff` | 49.33 µm²/s | 74–130 µm²/s (V79-171b 74, MGH-U1 130, PMID 9806526; SiHa 130 ± 5, PMID 17869669) | **0.379–0.667×** (1.5–2.6× **too small**) |
| `k` | 6.9182e-4 /s | 5.1667e-3 – 8.6667e-3 /s (0.31 and 0.52 /min, PMID 9806526) | **0.0798–0.1339×** (7.5–12.5× **too small**) |
| **`λ = sqrt(D/k)`** | **267.04 µm** | **92.40 – 158.62 µm** | **1.68 – 2.89× too long** |

**That derived `λ` row brackets the 2.45× Kyle implies for daunorubicin — a
different compound class, a different laboratory, a different assay.** It is the
strongest corroboration available and it was sitting in two unrelated rows.

**The split, log-additive:**

| bound | `λ` fold | `k` share of magnitude | `D_eff` share | direction of the `D_eff` term |
|---|---:|---:|---:|---|
| vs low `D`, high `k` | 2.890× | **86.2 %** | 13.8 % | **opposite** |
| vs high `D`, low `k` | 1.683× | **67.5 %** | 32.5 % | **opposite** |

### 6.1 On the preparation Kyle actually used, `D_eff` is RIGHT and `k` carries all of it

`results/cellular_exclusion.md` §4.1 and `scorecard.py` both name the
**MCL-appropriate arm** as the impermeable-cell term **alone**, because an MCL is
packed cells with little matrix. Computed live on that arm:

| | engine, MCL-appropriate arm | measured | verdict |
|---|---:|---|---|
| `D_eff` | **77.01 µm²/s** | 74–130 µm²/s | **INSIDE the band** |
| `k` | 6.9182e-4 /s | 5.1667e-3 – 8.6667e-3 /s | 7.5–12.5× too small |
| **`λ`** | **333.64 µm** | 92.40 – 158.62 µm | **2.10 – 3.61× too long** |

**On the preparation this anchor is measured in, the `D_eff` term contributes
NOTHING to the length error and `k` carries 100 % of it — and the error is
LARGER than on the shipping arm, because the shipping arm's `D_eff` error was
cancelling part of it.**

### 6.2 ⚠ THE LANDMINE — read this before "fixing" a scorecard row

**`D_eff` and `k` are both too small, and they push `λ` in OPPOSITE directions.**
`D_eff` too small **shortens** `λ`; `k` too small **lengthens** it. `k`'s error
is much larger, so `λ` is too long — but `D_eff`'s error is **currently
cancelling 14–33 % of it**.

**Therefore: a well-meaning next reader who opens
`scorecard.py`'s "effective diffusivity (tirapazamine)" row, sees a clean
1.5–2.6× miss, corrects it in isolation, and re-runs, will make the perivascular
gradient error WORSE by 1.22–1.62× — `λ` multiplies by `sqrt(74/49.33) = 1.225`
at the low edge of the measured band and `sqrt(130/49.33) = 1.623` at the high
edge — and every gradient, penetration-depth and reach number in the product
with it.** That is not a reason to leave `D_eff`
wrong; a model should not be right by cancellation. It is a reason that
**`D_eff` and `k` must be moved together, and that anyone who sees the gradient
improve after touching `D_eff` alone should suspect `k` moved too.**

### 6.3 The caveat on the `k` row, carried in full

`scorecard.py` declares the `k` row `NOT_TESTABLE`, and its stated reason is
sound: the measured constant is **anoxic bioreductive metabolism of a
hypoxia-activated prodrug**, while the engine's `k` is a generic
binding-plus-clearance sink — *"the process is not the same process"*.

That objection is about **mechanism**, not about commensurability of the number:
both sit in the `−kC` slot of the same equation, both were fitted from the same
slab flux, and the **length** they jointly determine is precisely what a
slab-flux experiment measures well (`D` and `k` separately are only weakly
determined by it, which is why the bands are wide). **The derived `λ` row is
reported here as a DERIVED row with that caveat attached and must never be
quoted without it.** It also means the engine has no *mechanism* by which it
could have got tirapazamine right, so this is **consistency, not vindication**.

### 6.4 What would separate them prospectively, rather than by borrowing anchors

Two observables with different powers of `D` and `k`:

* **Time to steady state.** The slowest slab mode relaxes at `k + D(π/2T)²` — at
  the engine's daunorubicin values `5.03e-3 + 3.4e-3 /s`, i.e. the two terms are
  comparable, so a **time-resolved** near:far series separates them cleanly. Kyle
  exposed for a fixed 1 h only, so this dataset cannot do it — but it names the
  experiment: **repeat the MCC near:far at 15 min, 1 h and 6 h.** The engine's
  prediction is already on record and is sharp: its registered sweep moves `R`
  from 3.65 (1 h) to 3.52 (6 h), a 3.5 % change, because it is equilibrated. A
  store-dominated tissue would collapse. **This experiment also settles §8.2.**
* **Steady-state mass balance.** At steady state `∫S = k∫c` exactly, so **mean
  tissue concentration per unit delivered dose depends on `k` alone and is
  exactly independent of `D`**. Any tissue:plasma AUC (`Kp`) anchor is therefore
  a **pure `k` observable**. The repository already has the machinery
  (`experiments/breast_kp_anchor.py`, `results/breast_kp_validation.md`) and has
  never read it this way. **This is the cheapest unbought evidence for the
  question in this file.**

---

## 7. Does the `λ` argument survive the real vessel field?

**It survives, and the amplification does not — which is the most
decision-relevant result here.**

The slab argument is 1-D. On the real section, distance-to-nearest-vessel is not
a slab coordinate, sources overlap, and the readout is a model-free empirical
half-distance, not a near:far ratio. Two independent traces of `d½(λ)` on the
**10x Xenium human GBM section, hBrain 319-gene panel, `dx12_repaired` arm,
`dx = 12 µm`, convection ON, 6 h, `dt = 10 s`, 10 µm distance bins, AUC field
binned by `direct_vessel_distance_um`**:

**(a) The four shipped compounds** (`λ` varies through both `D` and `k`; tissue
built once and reused, so the geometry is byte-identical):

| compound | `λ` (µm) | `d½` (µm) | `d½/λ` |
|---|---:|---:|---:|
| doxorubicin | 72.49 | 45.19 | 0.623 |
| daunorubicin | 64.40 | 38.91 | 0.604 |
| mitoxantrone | 53.96 | 31.49 | 0.584 |

`dln d½ / dln λ` = **1.222**.

**(b) A diagnostic `k`-scaling on the same section** (`D` fixed, `λ` moved 4×;
three full solves). The `k`×1.0 point reproduces the published 38.9 µm at 38.99:

| `k` scale | `λ` (µm) | `λ/dx` | `d½` (µm) | `d½/λ` |
|---:|---:|---:|---:|---:|
| 1.0 | 64.40 | 5.37 | 38.99 | 0.605 |
| 4.0 | 32.20 | 2.68 | 20.62 | 0.640 |
| 16.0 | 16.10 | **1.34** | 14.18 | 0.881 |

`dln d½ / dln λ` = **0.919** over the resolved pair (`λ` 64.40 → 32.20), **0.729**
including the third point. **The third point is resolution-limited, not physics:
at `λ/dx = 1.34` the gradient is not resolved on a 12 µm grid and `d½` floors out
near the grid pitch.** It is reported rather than dropped, and it should not be
read as a transport result.

**Conclusion.** On the real vessel field the map's readout responds to `λ` with
an exponent of **0.92–1.22** — essentially **linear** — against **2.92** for
Kyle's near:far ratio.

* **The exponential amplification is a property of Kyle's near:far estimator and
  does NOT transfer to the readout the product ships.**
* **A 2.45× length error appears on the map as 2.2–3.0×, not 14×.**
* Concretely: doxorubicin's shipped `d½` is 45.2 µm; its `λ` error against Kyle
  is 1.83×, so if Kyle is right the honest map number is **≈ 22–26 µm**. That is
  the sentence a customer should hear — *"half the drug is gone by about 25 µm,
  not 45"* — and it is a different conversation from *"the map is 14× wrong"*.
* `d½/λ` ≈ **0.60** on the real vessel field against **ln2 = 0.693** for a 1-D
  exponential: the 2-D overlapping-source geometry **compresses** the
  half-distance by ~13 % relative to a slab. That is a real geometry effect and
  it is small.

**Caveat, stated because it is load-bearing.** Arm B runs with **convection ON**,
so `d½` is not a pure function of `λ` — a Péclet number depends on `D`
separately. Trace (a) varies `D` and `k` together; trace (b) varies `k` alone.
They bracket the exponent at 0.92–1.22 and both are near 1, which is the result;
neither on its own would be sufficient.

---

## 8. Threats to the comparison itself

These are not decomposition terms. **All of §4 is conditional on Kyle's MCC being
a steady-state, linear, first-order-sink diffusion-reaction system — the same
model class as the engine. §8.2 shows that is not established.**

### 8.1 The instrument is badly conditioned, and the paper's own SD proves it

`R ≈ exp((T−W)/λ)`, so `dlnR = dT/λ`. Kyle reports culture thickness as
**140 ± 15 µm (SD)**. Propagated through Kyle's own readout:

| at | `T` = 125 µm | `T` = 140 µm | `T` = 155 µm | spread over ±1 SD |
|---|---:|---:|---:|---:|
| the engine's `λ` (64.40 µm) | 2.774 | 3.460 | 4.335 | **1.56×** |
| Kyle's implied `λ` (26.26 µm) | 28.247 | 50.001 | 88.518 | **3.13×** |

**A ±15 µm thickness SD alone produces a 3.13× spread in the ratio at the
measurement's own operating point. The paper reports the ratio as 50 ± 15, a
30 % SD.** Either the thickness spread is between-batch and was controlled per
culture, or the reported ratio SD substantially understates the uncertainty on
the quantity being inverted. Nothing in the paper resolves it.

Propagating **both** reported SDs:

| | daunorubicin `λ` (µm) | error vs the engine's 64.40 µm |
|---|---|---|
| ratio SD only, `T` fixed at 140 µm | 24.80 – 28.54 | 2.26 – 2.60× |
| **both SDs propagated** | **21.61 – 32.16** | **2.00 – 2.98×** |

**Quote the second row.** The point estimate 2.45× is not the finding; the band
is.

Other arm-A estimator sensitivities, at the arm-A configuration, none capable of
touching a factor of 14:

| axis | swept over | `R_auc` |
|---|---|---|
| grid pitch `dx` | 1, 2, 4 µm | 3.625, 3.646, 3.758 (**1.04×**) |
| edge window | 15, 30, 45 µm | 4.222, 3.646, 3.169 (**1.33×**) |
| exposure window | 600 s, 1 h, 6 h | 4.823, 3.646, 3.522 (**1.37×**, almost all below 1 h) |

### 8.2 ⚠ THE 1 h SNAPSHOT MAY NOT MEASURE A LENGTH AT ALL

**This is the largest single threat in this file. It does not correct the
inversion — it can void it.**

The engine is **fully equilibrated** at 1 h: `1/k = 155 s`, `T²/D_eff = 732 s`,
exposure 3600 s — and its own registered 1 h → 6 h sweep moves `R` by only 3.5 %,
which is the confirmation rather than the assumption.

A tissue with an intracellular **store** is not. A linear equilibrium store
retards `D` and `k` by a factor `R_f`, leaving `λ` **unchanged** but multiplying
the time to equilibrium by `R_f`. Measured — **pure retarded diffusion, ZERO
loss, Kyle's geometry (Dirichlet at 0, sealed at 140 µm), snapshot at 3600 s**,
whose steady state is exactly uniform (`R = 1`):

| `R_f` | `D_apparent` (µm²/s) | `τ = Dt/(R_f T²)` | near:far at 1 h |
|---:|---:|---:|---:|
| 1 | 26.79 | 4.920 | 1.000 |
| 10 | 2.679 | 0.492 | 1.490 |
| 20 | 1.339 | 0.246 | 2.756 |
| **50** | **0.536** | **0.098** | **13.600** |
| **90** | **0.298** | **0.055** | **84.096** |
| 200 | 0.134 | 0.025 | 6289.360 |

**A retardation factor of about 60–70 reproduces the measured 50 with `λ = ∞`
and no sink whatsoever.** The measured ratio is jointly determined by a length
and by whether the preparation has equilibrated, and **a single 1 h timepoint
cannot tell them apart.**

Is a store that large plausible for this compound? **Kyle's own paper says so, in
the sentence `results/kyle2004_penetration.md` §2.5 already quotes:** *"DOX
appears to accumulate the least over the 1-hour exposure with EPI and DAU
reaching ~1.5 and ~7 times higher levels, respectively, on the near sides."*
That is a direct observation of exactly the store the retardation model needs, in
exactly the compound with the extreme ratio. §9.2 of that file states the engine
represents cellular uptake as a **loss**, not a reservoir, and that the
retardation term that would make it one was built, **refuted**
(`results/intracellular_retardation.md`), and ships OFF.

**Say this at full volume: it is not established that Kyle 2004's near:far ratio
measures a transport length scale.** If the MCC is store-retarded at 1 h then
§4's inversion to `λ = 26.26 µm` is **void**, the "14×" is a **time-constant**
statement misread as a **length** statement, and the correct response is not to
shorten `λ` but to add a store. **This is not resolvable from Kyle 2004 alone and
must not be asserted either way.** §6.4's time-resolved experiment settles it,
and it is the only thing that does.

### 8.3 Free interstitial drug versus total tissue drug — signed, unquantified

Kyle measures **drug-derived fluorescence** = total tissue drug. The engine
solves for **free interstitial** drug.

* If the store is **linear, spatially uniform and equilibrated**, total is a
  constant multiple of free and it **cancels exactly** from the ratio. This is
  what §4 assumes.
* If the store **saturates on the near side** — and `results/two_route_promotion.md`
  §8.6 records Choi 2011 (PMID 21225235) measuring exactly that for paclitaxel,
  38 % penetration at 1 µM against 100 % at 10 µM, while anthracyclines bind DNA
  at finite site density — the near-side signal is sublinear, the **measured ratio
  is biased DOWN**, and the true free-drug ratio is **larger than 50**. **This
  direction makes the gap worse.**
* If it is **not equilibrated**, §8.2 dominates everything.

**No number is offered and none is invented.** The signed direction is stated
because a term whose sign is known is worth reporting even when its size is not.

---

## 9. What was planted, and what bit

| plant | expected to be caught by | outcome |
|---|---|---|
| move `ALPHA_EXPECTED` to 22.2 | `_assert_alpha()` | **CAUGHT** — `AssertionError: calibrate_alpha() = 22.223477711172563, expected 22.2` |
| mutate the module constant `ECM_ARM_A` to 0.4, re-run the reproduction check | §1.3's relative-difference-vs-published check | **NOT CAUGHT.** Relative difference stayed at 0.000e+00. |
| evaluate the analytic transfer function at a `λ` 1.5× wrong | §2's engine-vs-analytic agreement | **CAUGHT** — agreement goes from 1.0537 to 1.8819 |

**The second plant is the interesting one and it is recorded rather than quietly
fixed.** The reproduction check was **not** vacuous — **the plant could not reach
it.** `slab_readouts()` bound `ecm`, `dx_um`, `thickness_um`, `exposure_s` and
`window_um` as **default arguments**, evaluated once at import, so mutating the
module constant afterwards changed nothing the function could see. Re-planted
through the argument path, the check bites on every axis:

| plant, through the argument path | `R_auc` | rel. diff vs published |
|---|---:|---:|
| clean (published config) | 3.6457 | 0.000e+00 |
| `ecm` 0.4 | 4.5374 | 2.446e-01 |
| thickness 100 µm | 2.0204 | 4.458e-01 |
| window 20 µm | 4.0029 | 9.799e-02 |
| exposure 600 s | 4.8231 | 3.230e-01 |
| `dx` 8 µm | 4.1476 | 1.377e-01 |

The trap has been removed — those parameters now resolve at call time — so a
future reader planting the obvious way gets a bite. **The defect this repository
keeps meeting is not "the check is wrong", it is "the plant never reached the
check". That happened here, on the first try, in a file written by someone who
was watching for it.**

The §4.2 preparation correction was additionally run as a **registered
prediction** rather than a plant: the direction and approximate size were written
into the script's constants before the stage was executed, and the run confirmed
them.

---

## 10. Ranked: what would actually move the number

### 1. Settle whether Kyle's 1 h MCC is equilibrated. **A GATE, not a fix.**
Until this is answered, §4's entire inversion may be measuring a time constant
(§8.2). Cost: a literature question first (does any MCC study report a *time
series* of the near:far ratio?), and a prediction the engine has already made
second (`R` moves 3.5 % from 1 h to 6 h; a store-retarded tissue collapses).
**Everything below is conditional on this and it should be answered first.**

### 2. Read the existing `Kp` anchors as a `k` anchor. **Fixable, cheap, highest value per unit work.**
Steady-state mass balance makes mean tissue concentration per unit dose depend on
`k` **alone**, exactly independent of `D` (§6.4). The machinery exists
(`experiments/breast_kp_anchor.py`, `results/breast_kp_validation.md`) and has
never been read this way. §6 says the error is a `k` error; **this is the only
route to confirming it prospectively rather than by borrowing tirapazamine's
rows.**

### 3. Measure an MCL extracellular-space fraction. **Fixable, and it is the only legitimate way to use the largest un-anchored lever.**
`φ = 0.7892` is a fixed constant applied to a GBM section, an MCC and everything
else, and it can absorb the **entire** residual (§4.2). The repository already
searched for a measured MCL value and found nothing usable. **Until it is
measured, `φ` must not move — moving it to reduce a residual is tuning to the
anchor.** Already listed in `results/cellular_exclusion.md` §7.

### 4. The Primeau-vs-Kyle anchor conflict. **NOT fixable — state it permanently.**
The bands are disjoint (§5). Both measurements may be right and the engine's
one-`λ`-per-compound model may be the thing that cannot hold both — different
preparations, different input functions. **This is a structural scope limit, not
a bug, and 67.5 % of the daunorubicin gap is it.** **"Fit α better" is provably
not a fix**, and attempting it would be tuning to an anchor.

### 5. A bound-drug compartment. **Item 1 in disguise; do not add it to close the 14×.**
`results/kyle2004_penetration.md` §10.3 already ranks this second. §8.2 sharpens
why: a store does not merely add near-side accumulation, it **changes the
dynamical regime of the measurement**. `results/intracellular_retardation.md`
refuted a retardation term against a *free-drug diffusion coefficient*, a
different measurand from Kyle's total tissue drug, which does not inherit that
refutation. **Add it, if at all, because item 1 shows the regime is wrong — never
because it closes a residual.**

### 6. `D_eff`. **Measurably wrong, and fixing it alone makes this number WORSE.**
See the landmine, §6.2. It should still be fixed — a model should not be right by
cancellation — but **`D_eff` and `k` must move together**, and nobody should
expect the perivascular gap to improve when `D_eff` is corrected.

### 7. The estimator and the construction. **Fixable, not worth doing for this purpose.**
Worth 0.949× combined (§4.1), in the direction that flatters the model. Fixing
them correctly makes the published discrepancy slightly **larger**. Do it for
correctness — §6 of the Kyle file states the boundary bias with the **wrong
sign**, and §4's `L` column is printed at an ECM the run below it does not use
(§4.2) — but not expecting movement.

### 8. Reconditioning the anchor. **NOT fixable — permanent.**
`R ≈ exp((T−W)/λ)` amplifies everything, including the source's own ±15 µm
thickness SD, by 3.13× (§8.1). **A near:far ratio across a slab is a poor
instrument for inferring a transport length, and that will be true of any future
paper reporting one.** Stop quoting the ratio error; quote the length error with
its propagated band, **2.00–2.98×**, and on the product's own readout quote
**2.2–3.0×** (§7).

---

## 11. Reproducing this

```
./.venv/bin/python experiments/accuracy_gap_decomposition.py --stage verify
./.venv/bin/python experiments/accuracy_gap_decomposition.py --stage anchors
./.venv/bin/python experiments/accuracy_gap_decomposition.py --stage prep
./.venv/bin/python experiments/accuracy_gap_decomposition.py --stage offslab
./.venv/bin/python experiments/accuracy_gap_decomposition.py --stage all
```

`verify`, `anchors` and `prep` are slab and analytic work, a few minutes each.
`offslab` runs three diagnostic solves on the 10x Xenium human GBM section at
`dx = 12 µm` and is the expensive stage. `--stage all` writes
`results/accuracy_gap_decomposition_data.json`.

**Suite state.** Nothing outside `experiments/accuracy_gap_decomposition.py`,
this file and `results/accuracy_gap_decomposition_data.json` was added or
modified. No assertion was deleted or weakened anywhere. **The full suite was not
run**: three other agents are editing `app/` and `docs/` concurrently, so a suite
result now would be exactly the kind of under-contention number this repository
has already been bitten by twice this week.

**Hand-offs outside this file's territory, by SYMBOL, not line number:**

* `results/kyle2004_penetration.md` §6 — the sentence stating that the
  source-term boundary *"biases the modelled ratio downward relative to a clean
  fixed-concentration boundary"*. Measured, the sign is **upward** (§1.4).
* `results/kyle2004_penetration.md` §4 — the `L (µm)` column. It is computed at
  `ecm = 0.4`; arm A below it is solved at `ecm = 0.0875`. Either state the ECM
  in the header or print the arm-A values (§4.2).
* `tmesim/validation/scorecard.py::_run_mcl_transport_coefficients`, the
  `effective diffusivity (tirapazamine)` row — it invites an isolated fix that
  would make every gradient number in the product worse. The landmine in §6.2
  belongs in that row's `note`.
