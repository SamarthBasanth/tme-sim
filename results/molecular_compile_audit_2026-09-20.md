# Molecular compile audit — the two-route law, re-derived

**Date:** 2026-09-20
**Scope:** `tmesim/molecular/two_route.py`, `tmesim/molecular/compile_drug.py`, and the
bring-your-own surface in `tmesim/molecular/byo.py` that delivers their output.
**Branch/commit at audit time:** `main` @ `f9d443d`.
**Nothing in this audit changed a shipping number.** Every plant listed in §7 was
reverted byte-identical and the suite was re-run green afterwards.

```
calibrate_alpha() START = 22.223477711172563   ASSERT OK
calibrate_alpha() END   = 22.223477711172563   ASSERT OK
```

---

## 0. Pre-registered thresholds

Written down before the scoring was run, and reported against as-is.

| # | check | threshold | verdict |
|---|---|---|---|
| P1 | independent recomputation of R² reproduces the recorded 0.718 | \|ΔR²\| ≤ 0.001 | **PASS** (0.718311) |
| P2 | the recorded `n` matches the scored `n` | exact | **FAILED** (recorded 71, scored 68) |
| P3 | residual vs molecular weight carries no systematic term | \|t\| < 2 on the partial coefficient | **FAILED** (t = −2.80) |
| P4 | the delivered surface flags a compound outside the scored envelope | a flag exists on each scored axis | **FAILED** (no logP or HBD flag exists) |
| P5 | an input the law cannot handle is refused, not defaulted | refusal on every non-physical input | **FAILED** (5 of 16 probes return NaN silently) |

---

## 1. The pipeline: what the law consumes and what it emits

`build_profile` is the only wiring point; `two_route=True` is the shipping default
(`compile_drug._TWO_ROUTE_DEFAULT`).

| # | stage | module / symbol | consumes | emits |
|---|---|---|---|---|
| 1 | hydrodynamic radius | `stokes_einstein.hydrodynamic_radius_nm` | `mw`, `is_antibody` | `radius_nm` |
| 2 | free diffusivity | `stokes_einstein.free_diffusivity` | `radius_nm`, `T=310 K`, η=1.2e-3 Pa·s | `free_D` µm²/s |
| 3 | micro-tortuosity | `brownian.simulate_effective_diffusivity` | `free_D`, `ecm_reference=0.4`, `radius_nm`, seed 0 | `micro_D` |
| 3b | ECM mesh gate | `ecm_mesh.ogston_renkin_hindrance` | `radius_nm` | `mesh_h` → **`interstitial_D = micro_D × mesh_h`** |
| 3c-i | membrane leg | `two_route.membrane_permeability_um_s` | **`logp`, `hbd`**, `radius_nm`, `T` | `log K = min(logP − 1.34·HBD, 3.0)`; `P = K·D_mem/d` µm/s |
| 3c-ii | cytoplasm leg | `two_route.cytoplasm_diffusivity_um2s` | `radius_nm`, `T`, Kao 0.27, η_water 6.9e-4 | `D_cyto` µm²/s |
| 3c-iii | effective medium | `two_route.two_route_factor` | `interstitial_D`, `D_cyto`, `P`, **φ = 0.7892**, `a = 7.8517 µm`, `K = 1` | `hindrance`, `transcellular_gain`, `effective_D_um2s` |
| 3c | **the law** | `two_route.two_route_diffusivity` | the three above | **`effective_D = interstitial_D × hindrance`** |
| 4–7 | unchanged legs | `permeability`, `binding`, `potency` | descriptors | BBB, free fraction, binding, clearance, IC50 |
| → | | `MolecularProfile` / `DrugParams` | | `.diffusivity` = the two-route value; `.interstitial_diffusivity` always carries the pre-promotion value |

**Solute inputs to the law, in full: `mw` (only through `radius_nm`), `logp`, `hbd`.**
Nothing else about the structure reaches the diffusivity. Two opt-in, default-off
terms exist and are bitwise no-ops when omitted: `efflux_ratio` (K → K/ER) and
`bound_to_free_ratio` (capacity only — it cannot move a conductance, and the module
asserts that rather than claiming it).

---

## 2. R² = 0.718, re-derived independently — **REPRODUCES**

Re-derived from `data_cache/pruijn2008/pruijn2008_tableS1.csv` (the pdfplumber
extraction with its own `PROVENANCE.md`), **not** from the hardcoded
`PRUIJN_2008_TABLE_S1` literal in `experiments/diffusivity_falsification.py`, and
scored through the shipping code path. The end-to-end
`build_profile(two_route=True)` compile was cross-checked against the direct
`two_route_diffusivity` call on a 10-compound sample: **10/10 bit-identical**.

```
pearson r (log10 pred vs log10 meas)          = 0.847532
r^2 (squared Pearson — what the repo reports) = 0.718311     <- recorded 0.718
R^2 = 1 - SSres/SStot, log10, 1:1 line        = 0.708926
R^2 = 1 - SSres/SStot, linear, 1:1 line       = 0.721813
```

**The number is real and it is not an inflated statistic.** The repo reports a
squared Pearson correlation, which is invariant to any affine rescale of the
prediction and can therefore flatter a biased model. Here it does not: the
coefficient of determination against the **1:1 line** is 0.709, within 0.01 of it,
because the geometric bias is 1.032× and the OLS slope is 1.117. The model is
close to unbiased, so the two statistics coincide. That is a genuine result.

Baselines on the same n:

| model | r² (squared Pearson) | R² vs 1:1 |
|---|---|---|
| MW-only engine (`two_route=False`) | 0.0221 | −2.239 |
| paracellular only (impermeable cells) | 0.0221 | −0.645 |
| **shipping two-route law** | **0.7183** | **+0.709** |
| the paper's own Eqn 1 (*fitted to this data*) | 0.8393 | +0.831 |

Supporting: ratio pred/meas spans 0.227–2.987, median 1.034, **98.5 % within 3×**,
geometric bias 1.032×. Measured spread 29.69×; model spread 8.34×.

### 2a. **FAILED (P2): the scored n is 68, not 71**

Three compounds — **SN28871, SN28874, SN28875** — carry no logP7.4 in either the
measured or the calculated column of Table S1. `experiments/two_route_transport.py`
excludes them and prints the honest count (`68/71 compounds have a logP7.4`). The
**claim surface does not.** Every one of these says 71:

- `tmesim/molecular/two_route.py:37` — "R^2 = 0.718 against 71 compounds"
- `tmesim/molecular/compile_drug.py:29, 87` — "R^2 = 0.718 on the same 71 compounds"
- `tmesim/molecular/compile_drug.py:261` — **the provenance string attached to every compiled compound**
- `tmesim/molecular/byo.py:48, 53, 229, 668, 868` — the BYO surface, five sites
- `tests/test_byo_compound.py:1151`, `tests/test_time_resolved_pk.py:108`,
  `tests/data/static_source_reference.json`

This is the same defect shape as the thirteen-day ranking sentence: the engine is
right, the documents are wrong, and the wrong number is stapled to every compile.
It is small — 68 vs 71 — but it is a count of the evidence, and it is overstated.

### 2b. Not checked by CI, by design

`tests/test_claim_guard.py:85` states it outright: *"Numbers are out of scope.
`R² = 0.718` cannot be checked this way; the only way to know it is still 0.718 is
to re-run the scoring, which is a CI cost."* The string-matching guards in
`tests/test_byo_compound.py` assert `"0.718" in ...`, which pins the **text** and
in doing so also pins the **wrong n**. This audit is the first independent
recomputation; it agrees.

### 2c. Minor, unreconciled, conservative

The recorded "~24 % of the measured spread" (`two_route.py:37`,
`compile_drug.py:30/87/261`, `results/two_route_promotion.md:17`) does not reproduce
exactly: 8.34/29.69 = **28.1 %**. `(8.3−1)/(29.7−1) = 25.4 %` is closer. The
recorded figure **understates** the model, so it errs safe; flagged, not corrected.

---

## 3. Residual structure — **FAILED (P3): there is a missing molecular-weight term**

The aggregate R² hides it. Regressing `log10(meas)` on `log10(pred)` and adding one
descriptor at a time:

| added term | coefficient | t(65) | R² | verdict |
|---|---|---|---|---|
| — (base) | slope 1.117 ± 0.086 | — | 0.7183 | |
| **+ MW (Da)** | **−0.000647 ± 0.000231** | **−2.80** | **0.7486** | **SIGNIFICANT** |
| + log₁₀ MW | −0.3566 ± 0.1494 | −2.39 | 0.7410 | SIGNIFICANT |
| + logP | +0.0255 ± 0.0243 | +1.05 | 0.7230 | ns |
| + HBD | −0.0464 ± 0.0279 | −1.66 | 0.7298 | ns |
| + MW + logP + HBD | — | −5.06 / +4.28 / −2.74 | **0.8136** | all three significant |

Raw residual correlation: `r(resid, MW) = +0.326`, `r(resid, HBD) = +0.257`,
`r(resid, logP) = −0.199`.

The drift is monotone and visible without any statistics — geometric mean of
pred/meas by MW quintile:

| MW band (Da) | n | geo(pred/meas) |
|---|---|---|
| 60.1 – 213.0 | 13 | 0.809× |
| 213.0 – 249.0 | 13 | 0.912× |
| 250.0 – 289.0 | 13 | 1.201× |
| 290.0 – 324.0 | 13 | 0.986× |
| 329.0 – 539.0 | 16 | 1.276× |

**A 1.58× systematic swing across the mass range, biased low at the light end and
high at the heavy end.** Mass enters this law at exactly one place — the
hydrodynamic radius, which feeds Stokes-Einstein in water and in the bilayer core —
so the residual says that single channel is the wrong shape in MW, not that the
membrane physics is wrong. By contrast the logP quintiles are flat to within
1.33×–0.83× with no monotone trend, i.e. the *chemistry* axis the promotion added
is behaving; the *mass* axis it inherited is not.

Adding all three descriptors linearly recovers R² = 0.814 against the law's 0.718.
That ~0.10 of R² sitting in a plain linear correction is the size of the missing
term. (It is not a proposal: fitting those coefficients to this dataset is exactly
what the paper's own Eqn 1 did, and doing it here would forfeit the blindness that
makes 0.718 worth anything.)

---

## 4. The extrapolation boundary — **FAILED (P4)**

### 4a. Where the envelope actually is

The law is **not fitted** to Pruijn — no coefficient was solved on it. But it was
**scored and selected** against it (`HBD_LOGP_DECREMENT` was chosen with the sweep
in hand), so the scored envelope is the envelope that licenses the R².

| axis | min | p5 | median | p95 | max |
|---|---|---|---|---|---|
| MW (Da) | 60.1 | 177 | 279 | 497 | **539** |
| logP7.4 | **−2.62** | −2.01 | +0.20 | +1.88 | **+2.96** |
| HBD | 0 | 0 | 1 | 3 | **6** |
| solute radius (nm) | 0.2585 | 0.3706 | 0.4313 | 0.5228 | 0.5371 |
| log K_mem | −10.66 | −5.30 | −1.55 | +0.54 | **+1.62** |

**0 of 68 benchmark compounds hit the `LOGK_MEMBRANE_MAX = 3.0` clamp.** The clamp
is therefore a free, exact out-of-envelope detector: any compound that clamps is
outside the set the law was scored on, and the code already computes the flag —
`two_route_diffusivity` returns `membrane_clamped` on every call. Nothing reads it.

### 4b. What the surface actually tells a user

`byo.py` has **one** envelope check and it is on MW alone:

```
MW_MIN = 50.0   MW_MAX = 2000.0            -> refusal
MW_CALIBRATED_MIN = 130.078                -> warning
MW_CALIBRATED_MAX = 853.918
```

- **There is no logP envelope constant in `byo.py`. There is no HBD envelope
  constant.** Confirmed by inspecting `byo.__all__`. Those are the two axes the
  two-route law actually discriminates on.
- The one warning that exists is keyed to the wrong thing: it is the span of the
  **two curated libraries** (5-fluorouracil → paclitaxel), and its own text says it
  is about "the MW → hydrodynamic-radius law". It is not the Pruijn envelope and
  does not claim to be. The two barely overlap: Pruijn runs 60.1–539 Da, the warning
  band 130.078–853.918 Da.

Compiled through the real `compile_user_compound`:

| compound | MW | logP | HBD | log K | clamped | D_eff | outside Pruijn on | byo warnings |
|---|---|---|---|---|---|---|---|---|
| 5-fluorouracil | 130.1 | −0.80 | 2 | −3.48 | no | 49.35 | — | 0 |
| paclitaxel | 853.9 | +3.74 | 4 | −1.62 | no | 36.58 | **MW, logP** | **0** |
| bexarotene | 334.5 | +5.79 | 1 | **+3.00** | **yes** | 200.25 | **logP, log K** | **0** |
| sucrose | 342.3 | −5.40 | 8 | −16.12 | no | 33.71 | **logP, HBD** | **0** |

**Paclitaxel is 1.58× above the top of the scored mass range and 0.78 log units
above the top of the scored logP range, and receives not one warning** — because
its MW happens to be the exact upper bound of the curated-library band.
Bexarotene is pinned at the partition clamp that no benchmark compound reaches,
and is delivered silently. Sucrose sits 2.8 log units below the scored logP floor
and 2 donors above the HBD ceiling, silently.

An extrapolated prediction is being handed over carrying the same
"R² = 0.718, blind and unfitted" provenance string as an interpolated one. That is
the quiet failure this check was written to find, and it is present.

---

## 5. φ = 0.7892 — sensitivity

`DEFAULT_CELL_VOLUME_FRACTION` enters every compiled compound. Sweeping it across
all 68 scored compounds (the anchor question is a separate agent's; this is leverage
only):

| φ | urea D | TPZ D | acridanone D | model spread | geo(D/D₀) over 68 | r² vs Pruijn |
|---|---|---|---|---|---|---|
| 0.0000 | 430.83 | 287.44 | 277.54 | 2.32× | 3.2968× | 0.0221 |
| 0.4000 | 215.41 | 144.37 | 258.30 | 2.77× | 2.0532× | 0.6832 |
| 0.6000 | 132.56 | 89.28 | 249.01 | 4.30× | 1.5064× | **0.7266** |
| 0.7200 | — | — | — | 6.25× | 1.187× | 0.7236 |
| **0.7892** | **65.12** | **44.40** | **240.41** | **8.34×** | **1.0000×** | **0.7183** |
| 0.9000 | 29.71 | 20.83 | 235.46 | 17.44× | 0.6804× | 0.7064 |
| 0.9500 | 14.60 | 10.77 | 233.25 | 33.57× | 0.5108× | 0.6986 |

**Local log-log sensitivity at φ = 0.7892 (±1 %), and it is not uniform:**

| compound | logP | HBD | d ln D / d ln φ |
|---|---|---|---|
| urea | −2.01 | 4 | **−4.028** |
| tirapazamine | −0.38 | 2 | **−3.933** |
| acridan-9(10H)-one | +2.96 | 1 | **−0.148** |

Three findings.

1. **φ is not a rescale; it is the spread.** A hydrophilic compound moves ~4 % per
   1 % of φ, a lipophilic one ~0.15 % — a 27× difference in leverage. The
   compound-to-compound spread the promotion rests on goes 2.77× at φ = 0.40 to
   8.34× at φ = 0.7892. At the bottom of the literature band the law would resolve
   barely more than the **falsified** MW-only model's 2.32×.
2. **φ was not tuned to the benchmark.** r² peaks at φ ≈ 0.60 (0.7266) and the
   shipping value scores *below* that (0.7183). A tuned φ would have landed on the
   optimum; this one did not. That is a real piece of evidence for the law's
   blindness and it survives the audit.
3. **The shipping φ sits above the literature band its own docstring cites.**
   `two_route.py:125` calls 0.7892 "An UPPER estimate" and records "Literature
   tumour band implies 0.40–0.72". 0.7892 > 0.72. Taking the band at face value, the
   compiled diffusivity carries a **1.19×–2.05× geometric uncertainty** from φ alone
   (per-compound 1.01×–3.31×), and ±10 % on φ alone gives 1.21×/0.78× geometric
   (per-compound up to 1.41×/0.61×). The docstring argues the same segmentation bias
   inflates `a` and so under-states the transcellular route, i.e. the bias is
   conservative for the *gain* — but it is **not** conservative for the spread, which
   is the quantity the promotion claims, and the docstring does not say so.

---

## 6. Sentinel / refusal paths — **FAILED (P5)**

### 6a. Non-physical inputs to `two_route_diffusivity`

11 of 16 probes refuse correctly with a specific `ValueError`. **Five return NaN
silently:**

| input | result |
|---|---|
| `cell_volume_fraction` = 1.0 / −0.1 / **nan** | REFUSED ✓ (all three) |
| `solute_radius_nm` = 0 / −1 | REFUSED ✓ |
| **`solute_radius_nm` = nan** | **NaN SILENTLY RETURNED** |
| `hbd` = −1 | REFUSED ✓ |
| **`hbd` = nan** | **NaN SILENTLY RETURNED** |
| **`logp` = nan** | **NaN SILENTLY RETURNED** |
| `logp` = ±inf | accepted (→ D = 210.7 / 40.06); defensible physics, undocumented |
| `interstitial_D_um2s` = 0 | REFUSED ✓ |
| **`interstitial_D_um2s` = nan** | **NaN SILENTLY RETURNED** |
| `temperature_K` = 0 | REFUSED ✓ |
| `cell_radius_um` = 0 | REFUSED ✓ |
| **`cell_radius_um` = nan** | **NaN SILENTLY RETURNED** |

**The cause is a single idiom.** Every guard is written `if x <= 0: raise`, which
NaN passes through because every comparison against NaN is False.
`cell_volume_fraction` is the one parameter written in the negated form
`if not (0.0 <= phi < 1.0)` — and it is the one parameter that catches NaN. The
module already contains the correct idiom in one place and the unsafe idiom in
five. Inverting those five guards is a two-character change per site and cannot
move any valid compile.

### 6b. The silent default in `compile_drug.build_profile` — the more serious one

```python
float(descriptors.get("logp", 0.0) or 0.0)
float(descriptors.get("hbd", 0)   or 0)
```

A molecule with **no logP and no HBD descriptor is not refused. It is compiled at
logP = 0, HBD = 0** — and those are not neutral values:

- `log K_mem = 0.00` sits at the **85th percentile** of the benchmark's membrane
  openness (median −1.55, p95 +0.54). A *missing* descriptor produces a membrane
  more open than 85 % of the compounds the law was scored on.
- At MW 263 Da: complete descriptors compile to **D = 45.68 µm²/s**; the same
  molecule with the descriptors absent compiles to **D = 189.34 µm²/s** — a
  **4.14× silent inflation**, delivered with the full `R² = 0.718` provenance.

The inconsistency is sharp and internal: `experiments/two_route_transport.py`
**excludes** the three no-logP Pruijn compounds from scoring on the grounds that the
law cannot score them, while `build_profile` **compiles them with a confident
number**. The experiment and the shipping path disagree about whether a missing
descriptor is handleable.

(For completeness: an explicit `logp=None` does **not** silently default — it raises
`TypeError` in `permeability.py:198`, outside this audit's territory. So the hole is
specifically the *absent-key* path, and only because of the `or` fallback.)

---

## 7. Plant-test log

Every plant was applied to the working tree, confirmed **REACHED** (the guard ran
and named its subject), watched go **RED**, then reverted byte-identical and
re-run **GREEN**. No `git` command was used at any point.

**Baseline:** `tests/test_two_route.py tests/test_molecular.py` → **88 passed**;
`tests/test_byo_compound.py -k Ordering` → **3 passed**.

### Plant A — the withdrawn ranking claim, re-inserted into `byo.CAN_CONCLUDE`

The exact defect that shipped for thirteen days. Guard: `TestOrderingIsNotOfferedAsEvidence`.

```
FAILED tests/test_byo_compound.py::TestOrderingIsNotOfferedAsEvidence::test_no_byo_string_offers_an_ordering_at_any_grade
1 failed, 2 passed, 81 deselected in 20.46s

E           byo.CAN_CONCLUDE[0]
E             'Use it to RANK'   -> instructs the reader to produce an ordering
E             'ranking is the claim' -> names an ordering as the product's claim
E             'Use it to RANK'   -> an imperative whose object is an ordering, with
E                                   no prohibition in the same clause
E           compile_user_compound[doxorubicin]['can_conclude'][0]      <- and the
E           compile_user_compound[extrapolated-mw]['can_conclude'][0]     payloads
```

It fires on the module constant **and on every compiled compound's payload**, which
is what makes it a real guard rather than a string check. Reverted → **3 passed**.
**The guard is green and unregressed.**

### Plant B — `DEFAULT_CELL_VOLUME_FRACTION` 0.7892 → 0.7000

```
7 failed, 81 passed in 28.41s
FAILED ...TestShippingDefault::test_alpha_is_pinned_at_its_repinned_value
FAILED ...TestShippingDefault::test_drugparams_diffusivity_is_bit_identical
FAILED ...TestShippingDefault::test_the_repin_is_exactly_the_two_route_factor
FAILED ...TestShippingDefault::test_alpha_follows_D_through_its_own_closed_form
FAILED ...TestEffluxComposition::test_an_open_membrane_is_strongly_moved_by_efflux
FAILED ...TestEffluxComposition::test_the_sourced_decrement_closes_paclitaxels_membrane
FAILED ...test_parameters_match_the_pre_registration
  assert two_route.DEFAULT_CELL_VOLUME_FRACTION == p["cell_volume_fraction"]["value"]
E assert 0.7 == 0.7892
```

φ is pinned against its pre-registration, and the α pin goes red with it —
confirming in code that α rides on φ, which is what §5 measures. Reverted → 88 passed.

### Plant C — `HBD_LOGP_DECREMENT` 1.34 → 1.0

1.0 is the value that **maximises** the Pruijn score (r² 0.725 vs 0.718) and was
rejected for being unsourced. Planting it is planting the tuning the rules forbid.

```
6 failed, 63 passed in 16.03s
FAILED ...test_the_promoted_parameter_matches_its_own_registration
FAILED ...TestShippingDefault::test_alpha_is_pinned_at_its_repinned_value
FAILED ...TestShippingDefault::test_drugparams_diffusivity_is_bit_identical
FAILED ...TestShippingDefault::test_the_repin_is_exactly_the_two_route_factor
FAILED ...TestShippingDefault::test_alpha_follows_D_through_its_own_closed_form
FAILED ...TestEffluxComposition::test_the_sourced_decrement_closes_paclitaxels_membrane
```

The anti-tuning guard bites specifically on the R²-optimal value. Reverted → 88 passed.

### Plant D — a NaN-refusal guard, run against **unmodified** shipping code

Inverted plant: the guard is the thing being planted, to prove the hole in §6a is
real and reached rather than hypothetical.

```
5 failed, 69 deselected in 0.87s
FAILED ...TestPlantD_NaNMustBeRefused::...[solute_radius_nm]
FAILED ...TestPlantD_NaNMustBeRefused::...[logp]
FAILED ...TestPlantD_NaNMustBeRefused::...[hbd]
FAILED ...TestPlantD_NaNMustBeRefused::...[interstitial_D_um2s]
FAILED ...TestPlantD_NaNMustBeRefused::...[cell_radius_um]

E AssertionError: cell_radius_um=NaN produced a NaN diffusivity SILENTLY -- no
E   exception, no sentinel. A compiled compound would carry NaN into DrugParams
E   with the full 'R^2 = 0.718' provenance string attached to it.
```

Red on all five, against shipping code. The guard was **not** left in the suite: a
failing test is not a deliverable, and fixing `two_route.py` mid-audit while other
agents are live is out of scope. Reverted → 88 passed. **§6a is the finding; the
fix is a recommendation, not an action taken here.**

### Post-plant state

```
tests/test_two_route.py tests/test_molecular.py   -> 88 passed
tests/test_byo_compound.py -k Ordering            ->  3 passed
git status --short tmesim/molecular/ tests/test_two_route.py tests/test_molecular.py
  (empty -- territory byte-identical to f9d443d)
calibrate_alpha() END = 22.223477711172563        -> ASSERT OK
```

---

## 8. What this audit does and does not establish

**Establishes.** R² = 0.718 is real, independently reproduced from the CSV
extraction through the shipping code path (0.718311), and is not an artefact of the
squared-Pearson statistic — the 1:1 coefficient of determination is 0.709. φ was not
tuned to the benchmark. The anti-tuning and anti-ranking guards are live and bite.

**Does not establish.** Nothing here makes any absolute diffusivity defensible; the
absolute anchors remain irreconcilable and untouched by this work. No
cross-compound ordering is offered, implied, or graded anywhere in this document —
§3 and §5 are statements about the **law's functional form** and about a
**parameter's leverage**, and neither is evidence for or against an ordering, which
remains withdrawn without qualification. The residual-structure result is a
diagnosis of a missing term, not a proposed term: fitting one to this dataset would
destroy the blindness that is the only reason 0.718 is worth anything.

**Open, in priority order.**

1. §6b — the absent-descriptor path compiles a 4.14× inflated diffusivity with full
   provenance. This is the one that can put a wrong number in front of a user today.
2. §4 — no logP/HBD envelope flag exists; `membrane_clamped` is already computed on
   every call and is a free exact detector that nothing reads.
3. §2a — the `n = 71` overstatement at eleven sites including the per-compile
   provenance string.
4. §6a — five NaN guards written in the unsafe direction, in a module that already
   contains the safe idiom.
5. §5 — the shipping φ sits above the literature band its own docstring cites, and
   the spread claim is a strong function of it.
