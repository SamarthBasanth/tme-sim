# The interval, as a first-class output

**Date:** 2026-09-20
**Code:** `tmesim/validation/uncertainty.py` §7, `app/server/uncertainty_api.py`, `tests/test_uncertainty_band.py`
**Frozen constant:** `calibrate_alpha() = 22.223477711172563` — asserted at the start and the end of every measurement in this file, and at the top and bottom of the test module.

---

## 0. Why

`results/kill_calibration.md` §5.1 says it plainly:

> **"Absolute kill is defensible as an INTERVAL, not as a number, and the interval is wide"** — temozolomide 0.0048–0.1442, a 30.1x span.

That sentence has, until today, lived only in a prose caveat inside a results document. The engine reports single-valued metrics on every surface a reader actually looks at: `metrics.kill_fraction` in `result.json`, in the API, in the UI panel. A reader who sees `kill_fraction = 0.1122` and nothing else has no way to know that the *characterised* uncertainty around it spans a factor of fifteen to thirty.

This is the interval built as an output rather than a footnote.

**It does not make absolute kill fraction defensible.** Every interval this feature emits carries `quotable: false` with the reason attached. `results/cell_state_impact.md` §5 recommends not quoting an absolute kill fraction *or* a kill ceiling as a result on this section at all — no cell-state route clears its own measurability bar there — and nothing here overturns that. **The interval is published as the SIZE OF THE REFUSAL.** The prohibition is the subject of the feature, not an obstacle to it: the feature's job is to make the prohibition visible at the point of use, with its magnitude attached and its two causes told apart.

---

## 1. What the band covers

Exactly two terms, because exactly two have MEASURED spans recorded in `results/`.

### 1.1 The BIOLOGICAL bar — intrinsic potency between patients

| | |
|---|---|
| **What it is** | the inter-culture spread of measured ex-vivo temozolomide sensitivity |
| **Measurement** | median viability **82.4 %** at 100 µM, inter-culture **IQR 62.1–97.1 %**, n = **40** patient-derived GBM cultures (PMID 37620410 / GEO GSE232173) |
| **Stored as** | IC50 endpoints **151.1–1882 µM** about the anchored **362.1 µM** |
| **Recorded in** | `results/kill_calibration.md` §4 |
| **Kind** | **BIOLOGICAL.** Patients really do differ by this much. |
| **Does a better model shrink it?** | **NO.** No amount of modelling narrows variation between people. |
| **What shrinks it** | measuring **this patient's own** sensitivity — an ex-vivo viability assay on this tumour, which replaces a between-patient spread with a within-assay one |

Its own caveat travels with it: the response variable is an operator-*reproducibility* table with %CV up to 42.9 %, so the bar is if anything optimistic.

### 1.2 The EPISTEMIC bar — which cell-state route

| | |
|---|---|
| **What it is** | three different routes to a cell-state label, on the **same cells** |
| **Measurement** | `F_state` = 0.5502 (`oxygen`, shipping default) · 0.6490 (`markers`) · 0.9550 (`vessel_distance`) |
| **Recorded in** | `results/cell_state_impact.md` §3.1/§3.2, propagated in `results/kill_calibration.md` §5.1 |
| **As a factor on kill** | **×1.0000 – ×1.7357** (= 0.9550 / 0.5502) |
| **Kind** | **EPISTEMIC.** *"Nothing physical changed between those columns; only the label on each cell did."* |
| **Does a better model shrink it?** | **YES.** It is our ignorance, and a measurement deletes it. |
| **What shrinks it** | Ki-67 IHC/IF on a serial section registered to the Xenium coordinates; **any** hypoxia-response gene (CA9 / ADM / NDRG1 / SLC2A1) on the panel; a measured pO₂ (`results/cell_state_impact.md` §5.1) |

It is **one-sided**: the shipping default sits at the *bottom* of the measured span, so this bar only ever widens the prediction upward.

### 1.3 Why conflating them would be worse than no feature

Merging the two into one bar tells a reader two false things at once: that a better model narrows a bar no model can narrow, and that a bar one serial-section stain would delete is a permanent property of the tumour. The second is the more expensive lie — it says "stop working on the model" about the one term that is entirely ours.

So the separation is **structural, not editorial**: `assert_kinds_separated()` runs over every payload before it is returned, every component carries its own `kind`, `kind_meaning`, `what_shrinks_it` and `source_document`, and `tests/test_uncertainty_band.py` plants a merged payload and watches the guard go red (§5, Plant 1).

### 1.4 A third quantity, deliberately NOT combined

`results/kill_calibration.md` §4 also records a *tighter* temozolomide bar: the bootstrap 95 % CI on the **cohort median** viability (10,000 resamples of n = 40) → IC50 **218–545 µM**. It answers a different question — *how well is the typical GBM's potency pinned down* — and it is **epistemic-statistical**: it shrinks as n grows. It is served in the registry (`COHORT_MEDIAN_SPAN`) and is never folded into either bar. Use the wide bar for a patient whose own sensitivity has not been measured; use this one when asking about the cohort.

---

## 2. What the band does NOT cover

None of the following is in any interval this feature emits, and each is returned with the payload so a consumer reads it without opening this file.

1. **Model-form error.** The interval is computed *inside* one fixed model structure. **Plausibly the largest term, and not in these bars.**
2. **The absence of outcome validation.** The Xenium section is an untreated resection. Nothing here has been scored against a measured response, so this is *not* a claim that the true value falls in the band 90 % of the time — or any % of the time.
3. **Delivery / transport uncertainty.** `alpha` is fit once on ONE doxorubicin penetration anchor and sets the perivascular sink for every compound. Its span is not measured, so it is not in the band. (§1–6 of `uncertainty.py` put an *expert prior* on it instead — a different object; see §6 below.)
4. **MGMT.** `results/mgmt_scalar.md` measures a further 2.9–5.0x scale-down of temozolomide kill on real sections. Not applied: a temozolomide band here is an **upper bound in the MGMT dimension** on top of its stated width.
5. **Section- and patient-level data uncertainty.** One section, one segmentation, no replicate.
6. **Cross-compound comparison.** Withdrawn 2026-08-21 and not reinstated here in any form. Every entry point takes exactly ONE run for ONE compound; `test_no_entry_point_takes_two_compounds` asserts it structurally.

---

## 3. Which metrics got a band, and which were refused

| metric | verdict | why |
|---|---|---|
| **`kill_fraction`** | **BANDED** (`quotable: false`) | both recorded spans apply and both propagate onto it |
| **`kill_ceiling`** | **BANDED**, epistemic only (`quotable: false`) | the Hill term is saturated, so intrinsic potency drops out entirely — the biological bar is **structurally zero**, not missing |
| `penetration_index` | **REFUSED — UNCHARACTERISED** | neither recorded span touches transport. The only interval this repo has for it is the posterior-predictive CI in §1–6 of `uncertainty.py`, propagated from **expert priors on analyst choices**, not measured spreads. It also carries a systematic **resolution bias** — a shift of every draw the same way, not a widening — which cannot be expressed as a band at all. |
| `penetration_length_um` | **REFUSED — UNCHARACTERISED** | same reason; a fit to the run's own AUC-vs-distance field, which neither span moves |
| `mean_tumor_auc` | **REFUSED — UNCHARACTERISED** | an exposure, set by transport and the PK schedule; no measured span exists |
| `therapeutic_coverage_1x` | **REFUSED — a TODO, not a data gap** | the potency span *does* move it (coverage is exposure in IC50 units), but coverage is a fraction above a threshold rather than a mean of a smooth per-cell function, and the band-edge re-evaluation has **never been checked against the engine for it**. An unchecked propagation is not a measured band. The check is cheap and this metric should get a band once it is done. |

Refused metrics still return their **point**, plus `status: UNCHARACTERISED` and the reason — which is the honest answer, and the one a UI needs in order to show a bare number without implying it is tight.

### 3.1 Per-compound: which bar each one is entitled to

| compound | biological bar | grade | verdict |
|---|---|---|---|
| **temozolomide** | IQR over **40 patient-derived cultures** | `measured_between_patient` | **BANDED** |
| **cisplatin** | IQR over **34 GDSC2 glioma monolayer CELL LINES** | `measured_between_cell_line` | **BANDED, population named** — it is a proxy for between-patient variation, not a measurement of it. `results/kill_calibration.md` §4: cisplatin *"gets no new anchor"*, its bar is *"the spread of the anchor source, not a re-anchoring"*. It must never be quoted as a between-patient bar. |
| **doxorubicin** | quoted literature range 0.1–1 µM | `quoted_literature_range` | **REFUSED.** Two round numbers a decade apart are not a characterised spread: no population, no n, no quantile, no dispersion statistic. Reporting it as an error bar would manufacture a distribution out of a citation. |
| anything else | none recorded | — | **REFUSED**, and the remaining epistemic-only band is labelled `BANDED_PARTIAL` with *"NARROWER than the characterised uncertainty, not equal to it"* |

The epistemic bar is refused outright off the flagship section (`xenium_319`): three sections on three panels is a sample size of three, and which routes are even admissible is decided by each panel's gene content.

---

## 4. Worked example

The stored production run **`r_00099`** — the shipping configuration, `xenium_319`, dx = 12 µm, temozolomide, `state_mode='oxygen'`, single 50 µM dose / 4 h half-life / 6 h window, 23,812 tumour cells. Recorded engine runtime **40.57 s**.

```
calibrate_alpha() = 22.223477711172563     [START ASSERT OK]

POINT  kill_fraction = 0.11221232991768859

reproduction gate (PRE-REGISTERED, 4.0 sigma):
    expected 0.111306947297789   Bernoulli sigma 0.0020137584436948327
    z = -0.4496                  PASSED

BIOLOGICAL  (potency, between-patient, n = 40):
    IC50 endpoints 151.1 - 1882.0 uM   (anchor 362.1)
    kill 0.034254 - 0.307841           span 8.987x
    a better model does NOT shrink this

EPISTEMIC   (cell-state route, 3 routes on the same cells):
    factor x1.0000 - x1.7357
    kill 0.112212 - 0.194771           span 1.736x
    a better model DOES shrink this
    F_state this run 0.5422 vs recorded oxygen-route 0.5502  (1.46% apart)

COMBINED band: 0.034254 - 0.534329   span 15.599x
               quotable = False
               covers = ['biological', 'epistemic']

kill_ceiling  point 0.542164   band 0.5502 - 0.9550 (1.736x)
              point_inside_band = False   [disclosed, see PR-3]

REFUSED (uncharacterised): mean_tumor_auc, penetration_index,
                           penetration_length_um, therapeutic_coverage_1x

calibrate_alpha() = 22.223477711172563     [END ASSERT OK]
```

### 4.1 Reading it

**A point of 0.1122 with a characterised band of 0.0343–0.5343.** A factor of **15.6**. Of that, a factor of **9.0 is biological** — it does not go away, and the only thing that narrows it is an ex-vivo assay on this patient's own tumour. A factor of **1.74 is ours** — one Ki-67 stain on a serial section would delete it.

### 4.2 The span is NOT a fixed multiplier, and this is a finding

`results/kill_calibration.md` §5.1 records 30.1x for temozolomide at the 2026-08-12 configuration. This run measures **15.6x**. Both are correct. The potency band is an IC50 band, and the kill it implies is a **Hill** function of exposure: at the higher exposure this run reaches, the potent edge of the band is closer to saturation and the span compresses (17.4x → 9.0x on the potency term alone). **A band quoted as a frozen ±% would be wrong for every run but the one it was measured on.**

That is why the recorded *IC50* endpoints are stored and the recorded *kill* endpoints are not. The IC50 band is a property of the measurement and does not move with the engine; the kill endpoints were measured before two physics changes (`d77d169` area-conserving vessel source, `fd2a9e7` two-route transport) that moved every kill in the repository. Re-emitting them would be precisely the claim-surface drift this project keeps having. **Every band is re-propagated through the run's own Hill on the run's own per-cell exposure and cell state, every time.**

### 4.3 The ceiling's point falls outside its own band — disclosed, not patched

This run's own `F_state` is **0.5422**; the recorded `oxygen`-route endpoint is **0.5502**, 1.46 % away (a different alpha era and a different oxygen field). The band endpoints are the *measurement* and are **not widened to contain the point**; the payload sets `point_inside_band: false` and says why. Moving a measured endpoint to swallow a run's value would be tuning a record to fit a run.

---

## 5. The guards, and the plants that fired them

Every guard was first exercised on the **healthy** object and asserted green, so that the red result is a reaction to the plant and not to a check that never reached its subject. (That is the defect this repository found seven times in one week; see the 2026-08 vacuous-green audit.)

**Baseline: 33 passed.** Each plant below was applied to the production code, the suite was run, and the file was then restored and verified byte-identical.

| # | planted defect | file | result |
|---|---|---|---|
| 1 | `assert_kinds_separated` made a no-op | `uncertainty.py` | **3 failed** |
| 2 | doxorubicin's quoted literature range promoted to `BANDABLE_GRADES` | `uncertainty.py` | **2 failed** |
| 3 | reproduction gate widened to 1e9 σ | `uncertainty.py` | **1 failed** |
| 4 | a transcribed IC50 endpoint drifted from its record (151.1 → 200.0) | `uncertainty.py` | **1 failed** |
| 5 | route stops refusing a missing exposure payload | `uncertainty_api.py` | **1 failed** |

```
========= PLANT 1: assert_kinds_separated made a no-op =========
E           Failed: DID NOT RAISE <class 'AssertionError'>
FAILED tests/test_uncertainty_band.py::test_plant_merging_the_two_kinds_is_rejected
FAILED tests/test_uncertainty_band.py::test_plant_dropping_the_kind_label_is_rejected
FAILED tests/test_uncertainty_band.py::test_plant_marking_an_absolute_kill_quotable_is_rejected
3 failed, 29 passed

========= PLANT 2: doxorubicin's quoted range promoted to a band =========
>           assert body["potency_spans"]["doxorubicin"]["bandable"] is False
E           assert True is False
FAILED ...::test_doxorubicin_quoted_range_is_refused_not_banded
FAILED ...::test_spans_route_serves_the_registry_and_the_distinction
2 failed, 30 passed

========= PLANT 3: reproduction gate widened =========
>       assert iv["reproduction"]["passed"] is False, iv["reproduction"]
E       AssertionError: {'abs_delta': 0.05788884631968594,
E                        'bernoulli_sigma': 0.0020598090002146275,
E                        'expected_kill': 0.11861453084849494,
E                        'gate_sigma': 1000000000.0, ...}
FAILED ...::test_plant_corrupted_exposure_refuses_the_band
1 failed, 31 passed

========= PLANT 4: a transcribed span drifts from its dated record =========
>       assert (tmz.lo, tmz.hi, tmz.point) == (151.1, 1882.0, 362.1)
E       assert (200.0, 1882.0, 362.1) == (151.1, 1882.0, 362.1)
E         At index 0 diff: 200.0 != 151.1
FAILED ...::test_potency_span_matches_the_dated_record
1 failed, 31 passed

========= PLANT 5: route stops refusing a missing exposure payload =========
FileNotFoundError: .../r_00099/field_drug_auc.bin
FAILED ...::test_run_route_refuses_a_run_whose_payloads_are_gone
1 failed, 31 passed

--- all reverted; BOTH FILES BYTE-IDENTICAL TO PRE-PLANT ---
33 passed
```

Beyond the guards above, the test module also plants the *defect the feature is for*: a compound with no recorded span, a compound whose span is a quoted range, and a section the span was not measured on — each must refuse rather than interpolate, and each is asserted to do so.

---

## 6. Pre-registered thresholds, and one that FAILED

Written before any measurement (`PR-1`…`PR-4`).

**PR-1 — reproduction gate: 4.0 σ. PASSED.**
The band re-evaluates the run's own kill at the band-edge IC50s, so the re-evaluation must first reproduce the run's own `kill_fraction`. Tolerance is the exact Bernoulli sampling scale of the engine's own draw, σ = √(Σ p(1−p)) / n. Measured on two independent paths: **z = −0.4496** (stored run `r_00099`, re-read from persisted payloads at dx = 12 µm) and **z = +0.3874** (in-process, engine's own float64 fields, dx = 80 µm). Threshold basis: Bernoulli sampling at n ≈ 2×10⁴ is normal to high accuracy; 4 σ is a value-independent convention.
*Disclosure:* an exploratory probe returning z = −0.4496 was run **before** this threshold was written down. 4.0 σ was not chosen to accommodate it — any threshold ≥ 1 would have.

**PR-2 — default-on iff the band costs < 1 % of the run it annotates. PASSED on every real configuration; FAILED on the smallest possible run.**

| path | measured | run it annotates | fraction | verdict |
|---|---:|---:|---:|---|
| in-process (arrays in hand), flagship dx = 12 µm | **1.75 ms** | 40.57 s | **0.0043 %** | PASS |
| in-process incl. array construction, dx = 80 µm real run | **4.6 ms** | 45.65 s | **0.010 %** | PASS |
| **in-process, smallest test run (dx = 80 µm, 30-min window)** | **2.7 ms** | **0.23 s** | **1.18 %** | **FAILED** |
| HTTP route, re-read from persisted payloads (warm) | **1.9 ms** | 40.57 s | 0.0047 % | PASS |
| HTTP route, first call in a process (cold) | **0.71 s** | 40.57 s | 1.76 % | **FAILED** |

The FAILED rows are reported, not tuned away. The band's cost is **flat** — it solves nothing, decodes three payloads in ~0.3 ms and does two numpy passes over 23,812 cells — so the *ratio* depends entirely on the denominator. It is default-on-cheap for every run anyone actually makes and is not for a run shorter than about half a second. The cold-route row is a one-time per-process cost: `calibration.calibrated_drugparams` recompiles the molecule at **0.55 s** and was 100 % of the route's cost, so `uncertainty_api._drug_constants` memoises the two floats it reads out of it; warm calls are 1.9 ms.

**Shipping decision, from the measurement:** the band is **DEFAULT-ON** on the in-process path (`finished_run_from_result` + `run_interval_report`, 4.6 ms on a 45 s run), and the HTTP route is an ordinary opt-in `GET` with its cold-start cost stated above.

**PR-3 — state-route consistency: DISCLOSED, deliberately NOT gated.**
No defensible tolerance exists. `results/cell_state_impact.md` §2.2 measures this route's VIABLE-fraction swing across mask × dx as **unbounded** (the band empties at one setting), so any threshold here would be invented. The run's own `F_state` is therefore *reported* beside the recorded 0.5502 — 0.5422 on `r_00099`, 1.46 % apart — with no pass/fail. §4.3 above is the consequence.

**PR-4 — `calibrate_alpha() == 22.223477711172563` at start and end. PASSED** (§4, and at the top and bottom of `tests/test_uncertainty_band.py`).

---

## 7. What this feature is NOT

* **Not a validated predictive interval.** There is no outcome label on this tissue. A 15x band around an unvalidated quantity is a statement about the model's internal consistency and the measured spread of two inputs — nothing more.
* **Not §1–6 of the same module.** Sections 1–6 of `tmesim/validation/uncertainty.py` are a Monte-Carlo propagation of **expert priors on analyst choices** (typing threshold, KDE bandwidth, alpha, tortuosity) and say so. §7 carries only **measured, recorded** spans. The two must not be read as the same object, and the module docstring says so at the top of the file.
* **Not a licence to quote a kill fraction.** `quotable: false`, always, with the reason in the payload.
* **Not a cross-compound anything.** One run, one compound, structurally.

---

## 8. Reproduction

```bash
.venv/bin/python -m pytest tests/test_uncertainty_band.py -q      # 33 passed
```

Route mount (one line, in `create_app`, beside the other routers):

```python
app.include_router(build_uncertainty_router(jobs))
```

Routes: `GET /api/uncertainty/spans` (the registry and the two kinds; no run needed) and `GET /api/runs/{run_id}/uncertainty` (that run's intervals and refusals, re-read from `geometry.bin` + `field_drug_auc.bin` + `field_oxygen.bin`; **nothing is re-simulated**).
