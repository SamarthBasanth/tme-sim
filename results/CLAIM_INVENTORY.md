# CLAIM INVENTORY — what Vivum claims, and what backs each claim

**Compiled 2026-08-12.** This file exists because `results/` holds ~80 documents,
many of which correct each other, and no single page answered *"what does this
product actually claim, and what backs it?"* This is that page.

**Precedence.** Where this file and any other document disagree, check the date
first: several `results/` files are internally correct but describe a
configuration that has since changed. The ordering rule used throughout is
**later measurement on the shipped configuration wins**, and every row below
names the file it is taken from. Three files are **mid-flight and owned by other
workstreams** — `results/cfl_resolution.md`, `results/oxygen_joint_grid_test.md`,
`results/kill_calibration.md` §3.1. Anything drawn from them is marked
**(PROVISIONAL)** and may move.

**Status vocabulary.**

* **DEFENSIBLE** — a competent reviewer would accept it as stated. It has been
  attacked and survived, and the narrow wording below is the wording that
  survived.
* **CONTINGENT** — the measurement is real but the verdict depends on a
  configuration choice that is itself unvalidated, or rests on one section, one
  paper, or a margin smaller than its own known error.
* **RETRACTED** — measured false, or superseded by a later measurement.

---

## 0. The answer, in one paragraph

Vivum is a research instrument whose **defensible output is the relative
comparison of compounds on one section of tissue**. That claim is narrow, it has
been attacked from several directions, and it has held. Essentially everything
absolute has failed against external measurement: interstitial pressure by 41x,
the oxygen supply constant by ~4,570x, tissue diffusivity and hindrance against
multicellular-layer anchors, vascular area fraction by 2.34x (5.23x as the
scorecard still scores it) and microvessel density by 3.59x on a held-out
patient. The transport physics itself is
**almost entirely unvalidated against external measurement**: in `vivum validate
--quick`, **13 of 13 scored rows never invoke the solver**; in the full run
**21 of 33** do not; and of the rows where a *solved field* is compared to an
*external measurement*, there are **three, and all three FAIL** — interstitial
pressure 41x low, 5-fluorouracil partition 12.8x high, methotrexate partition
10.5x low. **The narrow claim is real. The broad one is not.** Selling the broad
one would be the fourth time in this project that claims outran evidence.

---

## 1. The claim table

### 1a. Comparative outputs — the product's core

| # | Claim, narrowest true form | What backs it | What would falsify it | Status |
|---|---|---|---|---|
| C1 | **Scored ratio outputs are largely insulated against *common-mode, multiplicative* parameter error.** Correcting the transvascular delivery law changed absolute delivery by up to **9.58x** (5-FU) while the scored 3-D:2-D ratio moved **1.37x** — 84–86 % of the delivery error cancels for the three worst-affected drugs. | `breast_kill_revalidation.md` §5, one Xenium breast section (Janesick rep1), dx = 12 µm. Per-drug: 5-FU 9.58x→1.37x (13.8 % residual), cisplatin 8.86x→1.43x (16.4 %), cyclophosphamide 6.35x→1.45x (20.1 %). **Zero verdict flips.** | A *shape* perturbation rather than a scale one. Already found: `delivery_source_geometry.md` §5.3 shows a bit-exact dose match removes ~1 % of the movement on breast and **none** on GBM. | **DEFENSIBLE, narrow.** See note (i). |
| C2 | **Cross-drug ranking is invariant to common-mode perturbations** — cell-state route, oxygen route, α, and the exposure unit. | `cell_state_impact.md` §3.3 (all three state routes); `oxygen_level_shift.md` §4; `ranking_robustness.md` line 15 (a 3.36x change in the exposure unit moves decision probabilities by ≤ 0.004). | A drug-*specific* transport parameter. Already found — see C3. | **DEFENSIBLE, narrow.** |
| C3 | ~~Ranking is unconditionally invariant.~~ | — | — | **RETRACTED — measured false.** Correcting the delivery law changes the **top-1 GBM pick, cisplatin → vincristine**; top-3 changes; 4 of 45 pairs swap, two crossing gaps of 40.5 % and 35.6 % of the leader's kill against a Bernoulli sd ~0.003. P(cisplatin rank 1) **0.9277 → 0.3965**. `transvascular_exchange_revalidation.md` §4.1. Independently, a source-geometry change swaps cisplatin↔temozolomide at rank 1 in **both** raw and dose-matched arms (`delivery_source_geometry.md` §5.2). |
| C4 | **Therapeutic coverage is *exactly* invariant to the oxygen field.** | `oxygen_level_shift.md` §4 — verified with `==`, not `isclose`, across both vascular modes and 3 drugs. Structural: oxygen is not an input to the transport solve at all. | Nothing, for coverage — it is a structural identity. | **DEFENSIBLE.** |
| C5 | ~~*Ranking* is exactly invariant under oxygen perturbations.~~ Ranking is invariant *in measurement*, under a **near**-common multiplier. | Relative kill changes span **71.4–77.0 %** — a 5.6 pp spread across drugs against a 74 % common shift (`oxygen_level_shift.md` §4). `kill_fraction` is a covariance, not a scalar multiple. | A non-empty hypoxic band. `core/config.py:264` sets `TMZ_HYPOXIA_RESISTANCE = 0.80` against a global 0.35 — a genuinely **drug-dependent** oxygen→resistance map, inert today only because GBM HYPOXIC = 0.0000 at every dx. | **CONTINGENT.** The "drug-independent per-cell multiplier" premise is false in general; it is true only while the hypoxic band is empty. |

**Note (i) — how narrow C1 really is.** The three headline cancellation numbers
are not three independent probes: two are the breast 3-D:2-D ratio and one is the
GBM 3-D:2-D ratio, all of the same metric family, and all three cancel
*multiplicative scale* error, which the `plasma = c_bath / p99(e_unit)`
renormalisation removes **by construction** (`dtissue_sensitivity.md` §5 proves
this algebraically). Corrections to the numbers as commonly quoted:

* The "**3.07x diffusivity change, 99.6 % cancelling**" is mislabelled twice. The
  *input* is a **0.325x** change in `D_tissue` (1.404e-06 → 4.559e-07 cm²/s); 3.074x
  is the *exposure response*. **99.6 % is not printed anywhere** in
  `dtissue_sensitivity.md`; it is derived from the §7.3 row 3.074x → 1.004x, and
  applies to **doxorubicin and epirubicin only, `recal` mode only, dx = 12 only**.
  `grid_dependence_audit.md` §2.1 computes **97.8–98.6 %** over the full range and
  explicitly flags that 99.6 % is not in the source file. In **`frozen`** mode the
  same document measures **anti-cancellation** (surviving vs exposure 7.25–21.34).
* The cancellation exists in `recal` because **α is re-fit to the same doxorubicin
  anchor**, re-exporting a D change as a global multiplier. `dtissue_sensitivity.md`
  §8: *"D_tissue and the perivascular sink strength are not separately identified."*
  It is partly a gauge artefact, not purely robustness.
* The "**dt halving, 0.0003 %**" is from `grid_convergence_gbm.md` §7.3 (GBM, not
  breast) and its own caveat is explicit: it licenses **nothing** about absolute
  AUC or absolute kill, only that the *ratio* is dt-insensitive.
* The residual is not zero anywhere: 13.8 % at best, **~48 %** for the five
  low-permeability drugs. `breast_kill_revalidation.md` §5: *"'insulated' is not
  'immune'."*

### 1b. The 3-D:2-D kill panel — the main scored result

| # | Claim | Backing | Falsifier | Status |
|---|---|---|---|---|
| K1 | **GBM 3/3 PASS** at the shipped configuration. | `kill_calibration.md` §3: TMZ 0.2626 (band 0.1590–0.4622), cisplatin 0.2791 (0.2086–1.5004), doxorubicin 0.2077 (0.0575–1.0000). | Any change to the oxygen level or the source geometry. **Already found.** | **CONTINGENT.** |
| K2 | **Breast 5/6 PASS, docetaxel FAIL.** | `kill_calibration_breast.md` §3. | Grid, normaliser, source geometry. **All found** — and the shipped scorecard does not reproduce this count. | **CONTINGENT, and inconsistent with `vivum validate`. See note (iii).** |
| K3 | **Docetaxel's breast FAIL is the project's falsification instance.** | `kill_calibration_breast.md` §3: ratio **0.1484** vs band ceiling **0.1463** — **+1.5 % outside**. | Anything moving the ratio by >1.5 %. | **CONTINGENT → the shipped-convention verdict is UNDETERMINED.** See notes (ii) and (iii). |

**K1 contingency, from `kill_calibration.md` §3.1 (PROVISIONAL).** Panel counts
across five measured alternative configurations: **2/3, 2/3, 3/3, 3/3, 1/3.**

* **2/3** under the historical KDE oxygen field (shipped until 2026-07-28) — TMZ 0.46784 vs ceiling 0.46221, 1.2 % over.
* **2/3** under a literature-derived `supply` (194 s⁻¹) alone — TMZ **0.46784**, bit-identical to the row above. These two are independent in provenance but land on the **same lever** (oxygen level), so they are not two independent confirmations.
* **3/3** under joint `supply`+`consumption` calibration, but **cisplatin clears its floor by only 3.7 %**.
* **1/3** under `geometry='subgrid'` delivery source — TMZ *and* cisplatin drop.

`grid_convergence_gbm.md` adds that the 3/3 survives **on the width of the
literature bands, not on numerics**: movement/margin ratios 0.37x / 0.65x / 0.41x,
and dx = 12 sits at the **maximum of `F_state` over the whole ladder**, biasing all
three ratios high by +1.2 / +2.4 / +2.7 %.

**Note (ii) — docetaxel, stated precisely.** It is the only FAIL among the **9
published 3-D:2-D verdicts**, not the project's only falsification (see §1e). The
perturbations it survived are real:

* corrected delivery law: **1.5 % → 2.6 %** outside (`breast_kill_revalidation.md` §3 — *not* 3.7 %);
* 6.6x cellular-exclusion hindrance correction, **`recal` arm**: **1.5 % → 3.7 %** outside (`cellular_exclusion.md` §5.5) — but it flips **FAIL → PASS** in the `frozen` arm, and +2.19 % is smaller than the known −3.4 % grid error;
* `tail90` normaliser: **+1.5 % → +35.2 %** (`tail90_falsification.md` §7.1–7.2 — *not* `p99_convention.md`), a decisive FAIL — but `tail90` **has not been adopted**.

The perturbations it does **not** survive:

* **a 2x grid coarsening** — verdict changes between dx = 12 and dx = 24 (`dtissue_sensitivity.md` §7.2 + limitation 7); FAIL→PASS between dx = 12 and 16 (`grid_convergence_breast.md`). `grid_dependence_audit.md` §3.3 gives a **grid safety factor of 0.08**;
* the sub-grid delivery source removes it entirely, 0.1484 → 0.0752 (`delivery_source_geometry.md` §0).

`p99_convention.md` §6.1 therefore states the shipped-convention verdict outright:
*"The verdict is **UNDETERMINED at dx = 12 µm** — not FAIL, and not PASS."* That is
the current position under the convention the engine actually runs. **Do not
present docetaxel as a clean falsification instance.** The honest version is
stronger anyway: it is the only one of nine verdicts whose result is a function
of the grid at all, and the project found and published that itself.

**Note (iii) — NEW, found live on 2026-08-12: the shipped scorecard does not agree
with the breast report.** `vivum validate` (full, shipping physics) returns the
breast kill panel as **3 PASS / 3 FAIL**, not 5/1:

| drug | predicted | band the scorecard used | verdict |
|---|---|---|---|
| docetaxel | 0.1484 | 0.0152–0.1463 | FAIL |
| **5-fluorouracil** | 0.3379 | **0.1636–0.1636** | **FAIL** |
| **paclitaxel** | 0.1503 | **0.0412–0.0412** | **FAIL** |

The two new FAILs are scored against **zero-width bands** derived from a single
literature row (n = 1, "no spread to test against", as the row's own note says).
`kill_calibration_breast.md` §3 reaches 5/6 by **falling back to the broad
cross-tumour band** for exactly these two drugs; the scorecard does not fall back.

Neither is obviously wrong, but they are different bars and the discrepancy is
undisclosed. Note that a zero-width band on a continuous quantity is the exact
mirror of the defect for which the GBM IFP row was declared **NOT TESTABLE** — a
bar nothing can fail, versus a bar nothing can pass. By the scorecard's own stated
philosophy these two rows should probably be **NOT TESTABLE**, not FAIL. Until
that is settled: **do not state "breast 5/6 PASS, one FAIL" on any surface**, and
do not describe docetaxel as the only FAIL in the breast panel — the shipping tool
prints three.

### 1c. The spatial spread

| # | Claim | Backing | Falsifier | Status |
|---|---|---|---|---|
| S1 | Tumour cells in the macrophage-densest quintile retain **0.4291** of exposure vs **0.8592** in the sparsest — a raw structure ratio of **2.0025x**. | `immune_sink_impact.md` §5.1. One GBM section, one drug, one sink strength, dx = 20 µm. | — | **DEFENSIBLE as a raw statistic** ("reproduces to the digit"). |
| S2 | The excess attributable to **where these macrophages actually are** is **1.08–1.21x** over clustering-preserving nulls. | §5.3: clumpiness-matched permutation **1.6581** → excess **1.208x**; 400 µm block shuffle **1.8522** → excess **1.081x**. Observed sits above the **entire range** of both nulls, 20/20 draws. | A third spatially-aware null reaching 2.0025; or the effect failing to reproduce on a second section. | **CONTINGENT.** |
| S3 | ~~The comparator is a well-mixed model's 1.0.~~ | — | — | **RETRACTED** by §5.3 and by the measured well-mixed floor of **0.9343**. |

**S2's real caveats, all from `immune_sink_impact.md`.** (a) **p = 0.0476 is the
20-draw floor (1/21)** — the file says it *"must never be quoted as a resolved
p-value"*; there is **no confidence interval anywhere**. (b) **One section, one
drug, one strength.** (c) **82 % of the headline is a depth confound**: binning by
distance-to-tissue-edge with *zero* immune input already gives **1.6394** (§5.1a).
(d) Monotonicity across quintiles is **not** evidence — every control reproduces
it, including uniform scatter (§5.4). (e) The statistic is **refused in code** on
spot data (`require_measured_cell_positions`). (f) **The module ships default-off
and no shipped code path imports it** — this is a study, not a product output.

### 1d. Absolute quantities — the failure column

Every row here is scored against an independent published measurement. None is
defensible; they are listed so no one has to rediscover them.

| Quantity | Predicted | Measured | Error | Source |
|---|---|---|---|---|
| **Interstitial fluid pressure, breast** | 0.370 mmHg | 15 ± 9 mmHg (PMID 1423283) | **41x low** (16x below band floor of 6) | `ifp_validation.md` §4.2 |
| IFP, GBM | 0.462 mmHg | 2.0 ± 2.5 mmHg | 4.3x low — but **NOT TESTABLE**: a model outputting zero would pass | `ifp_validation.md` §4.1 |
| Oxygen `supply` constant | 0.05 | ≈194–251 s⁻¹ (PMID 31423930) | **~4,570x low** on GBM (**2,160x** on breast), partly compensated by a sink **112–166x** low | `oxygen_supply_calibration.md` §3.3, §5 |
| Effective diffusivity, MCL | 287.3 µm²/s | 74–130 µm²/s | FAIL | `anchor_registration.md` §1 |
| Hindrance `D_free/D_eff` | 1.77x → **6.61x** after cellular exclusion | **12x** (PMID 9806526) | **1.8x short** after correction | `cellular_exclusion.md` §0, §4.1 |
| Breast Kp, 5-fluorouracil | 7.7956 | 0.61 ± 0.11 (PMID 9205062) | **12.8x high** | `breast_kp_validation.md` |
| Breast Kp, methotrexate | 0.0571 | 0.60 ± 0.20 (PMID 9679959) | **10.5x low** — opposite direction, so no constant absorbs both | `breast_kp_validation.md` |
| Vascular area fraction, held-out patient | 10.982 % | **4.70 ± 2.38 %** (Marien 2016, PMID 27583442) | **2.34x** — but the **scorecard still scores 5.23x**, see below | `vascular_segmentation_audit.md` §2.3 |
| Microvessel density, held-out patient | 233.27 /mm² | 65 (35–95) /mm² | **3.59x** (s1 1.70x, s3 4.93x — the merge factor is not a constant) | `vascular_segmentation_audit.md` §7.3, `heldout_patient.md` §2 |

**The vascular-area row is a trap, and it cuts against us.** The previously quoted
**5.23x** error used an **incommensurable anchor** — Tamma 2022 (PMID 36611806)
measures CD34 stain-positivity 2.1 %, where the lumen contributes no stained
pixels while the model paints the whole vessel profile. Against the commensurable
Marien 2016 areal fraction the error is **2.34x**. But `vascular_segmentation_audit.md`
§5 states the historical anchor was left **untouched** and the commensurable one
**added as a new key**, and this is confirmed in code: `scorecard.py:1417` and
`:1740` still read `vascular_area_fraction_pct` (the 2.1 % CD34 value), while
`vessel_areal_fraction_pct` appears only in `reconstruct/vessel_network.py`. **So
`vivum validate` still prints 5.23x**, and `ifp_validation.md` §4.2 — written after
the audit — still cites 5.23x. Do not quote 2.34x as the product's number until
the scorecard scores it. Either way it is a FAIL.

### 1e. What the validation suite actually says

Run live on 2026-08-12, `vivum validate --quick`, shipping physics
(`transvascular_exchange=False`):

```
SUMMARY   PASS 7   FAIL 6   NOT TESTABLE 4   DEFERRED 8   (of 25 rows)
          13 of the 13 scored rows never invoke the transport PDE.
```

Full run, also executed live on 2026-08-12 (exit code 1), reproducing
`anchor_registration.md` §0 exactly: **PASS 13 / FAIL 20 / NOT TESTABLE 17 /
DEFERRED 0 over 50 rows; 33 scored, 21 of them never invoking the PDE.** So:

* the "**13 of 13**" figure is a **`--quick`** statement — in quick mode there are
  **zero** PDE-scoped scored rows. Row-by-row enumeration of the registry confirms
  it: 2 + 1 + 6 + 1 + 1 + 2 = 13 scored rows, every one `scope=analytic`;
* the full run does have **12 PDE-scoped scored rows**, so the solver is **not**
  wholly untested. They are 2 (`breast-transvascular`) + 1 (IFP breast) + 9
  (3 GBM kill drugs + 6 breast kill drugs);
* `anchor_registration.md` §7 says *"the only registered anchor that scores the
  solver against an external measurement on the strength of a solve remains
  `ifp-wick-in-needle`'s breast row"*. **That phrasing is too strong.**
  `breast-transvascular` does exactly that — free interstitial:plasma AUC vs human
  microdialysis, the PDE's own state variable, via `simulate_response`. The
  defensible narrow version is that **IFP is the only PDE anchor with zero fitted
  coefficients**; the kill and transvascular groups all cite the `alpha-fit`
  conditional.

**This correction makes the picture worse, not better.** Counting properly, there
are **three** scored rows in which a *solved field* meets an *external measurement*
— IFP breast, 5-FU breast Kp, methotrexate breast Kp — **and all three FAIL**
(41x low; 12.8x high; 10.5x low). The remaining 9 PDE rows compare a *ratio* to a
literature *band* and are conditional on the one fitted constant. The held-out
patient's 10 PDE solves are deliberately **not registered**, because they are
scored against no measured number at all — they compare two sections with each
other. Refinement makes the IFP failure worse, not better: breast IFP is monotone
decreasing across dx = 24/18/12/8 (0.621 → 0.424 → 0.370 → 0.286), a **2.17x**
spread (`ifp_validation.md` §5).

The full run's 12 PDE-scoped scored rows, verbatim: `breast-transvascular`
5-FU **FAIL** and methotrexate **FAIL**; `gbm-kill-3d2d` cisplatin/doxorubicin/
temozolomide **PASS ×3**; `breast-kill-3d2d` cisplatin/doxorubicin/epirubicin
**PASS ×3** and docetaxel/5-FU/paclitaxel **FAIL ×3**; `ifp-wick-in-needle` breast
**FAIL**. That is **6 PASS / 6 FAIL** on the solver-invoking rows, and every PASS
is a 3-D:2-D ratio conditional on the one fitted constant α.

**Roughly 18 distinct FAIL rows exist across `results/`** (breast Kp ×2, MCL ×2,
held-out geometry ×2, IFP, TMZ real-tissue d½, doxorubicin d½, regional PK ×2,
analytic penetration ×2, vascular heterogeneity, cross-platform composition,
Moncayo explant ×2, oxygen supply). Docetaxel is one falsification among many —
the project's honesty record is better than "one instance", and should be stated
that way.

### 1f. Numerics and resolution

| # | Claim | Status |
|---|---|---|
| N1 | Mass conserved to ~1e-13 on the real irregular mask; concentrations never negative. | **DEFENSIBLE** (`numerical_rigor.md`). |
| N2 | ~~dx = 12 µm is grid-converged.~~ | **RETRACTED.** Breast: docetaxel is **−3.4 %** from its Richardson limit, boundary layer `L` = 28.2 µm = **2.4 cells** at dx = 12, ladder spread **49.7 %** (`grid_convergence_breast.md`). GBM: the ladder is **non-monotone**, observed orders **0.96** and **−3.33**, GCI 5.44–6.42 %, power-law fits flagged **DEGENERATE** (`grid_convergence_gbm.md` §0). Neither panel is in an asymptotic regime. |
| N3 | ~~dx = 12 is the finest grid satisfying the advection CFL at dt = 10 s.~~ | **RETRACTED.** Measured CFL 0.663 at dx = 10; the GBM floor is **dx = 7 (CFL 0.955)** (`grid_convergence_gbm.md` §2.1; `cfl_resolution.md` **PROVISIONAL** re-measures dx = 8 at 0.939 and dx = 7 at 0.955 on the shipped configuration). |
| N4 | ~~`penetration_index` discriminates drugs.~~ | **RETRACTED** and retired — it spans 0.1946–0.1959 across a panel whose delivered exposure differs by 120x (`delivery_metrics.md`). |
| N5 | ~~The mesh is not the bottleneck once dx ≲ 12 µm.~~ | **RETRACTED by `numerical_rigor.md`'s own correction block** — the elasticities were mis-attributed (|E_pen| 0.747 belongs to `grid_dx_um`, not the endothelial threshold; |E_kill| 0.610 belongs to α) and the §2 grid trend had the **sign backwards**. |

### 1g. The negative claims — and these are worth as much as the positive ones

All **DEFENSIBLE**, all backed, and all should stay on every surface:

* **Vivum does not predict clinical response or patient outcome.** Backed by
  `DATA_HUNT.md`: a 39-agent, 164-dataset sweep found **no open dataset anywhere**
  pairing measured spatial drug distribution with measured treatment response in
  the same tumour, and **no glioma+temozolomide spatial drug map at all**.
* **Absolute kill fraction is not a defensible output.** Backed by §1d and by the
  0.3252 spread across state routes.
* **The spatial spread statistic is invalid on spot data** — refused in code.
* **The top-1 drug call is not decision-grade** — cisplatin 0.241 vs doxorubicin
  0.231, ~4 % apart; and P(rank 1) collapses 0.93 → 0.40 under the corrected
  delivery law. `ranking_robustness.md` §8: *"the reportable output of this
  pipeline for the top of the panel is a SET, not a winner."*

---

## 2. What a competent reviewer would accept today

**Would accept:**

1. The data provenance. All open, all cited, no synthetic data in the package.
2. **Relative comparison of compounds on one section**, presented as a ranking
   with a stated tie at the top and stated invariance limits.
3. Coverage's exact invariance to the oxygen field (a structural identity).
4. The solver's internal numerics — conservation, positivity, dt behaviour of ratios.
5. Every negative claim in §1g. These are the strongest part of the product.
6. The honesty apparatus itself: `vivum validate` prints failures three times,
   refuses to score a solver number that has not declared its delivery law, and
   exits 1 on any FAIL. This is unusually good and is a genuine selling point.

**Would not accept:**

1. **Any absolute number** — pressure, oxygen, diffusivity, partition coefficient,
   vascular morphometry, kill fraction. §1d is a wall of failures.
2. **"The transport physics is validated."** It is not. Three scored rows compare
   a solved field to an external measurement, and all three fail — in both
   directions, so no single constant can absorb them.
3. **The 3-D:2-D panel verdicts as robust.** GBM 3/3 → 2/3 or 1/3 under three
   perturbations; breast docetaxel is UNDETERMINED at the shipped resolution; and
   the shipped scorecard returns the breast panel as **3/6, not 5/6** (note iii).
4. **dx = 12 as a converged production grid.**
5. **The spatial spread as a headline finding** — one section, p at the 20-draw
   floor, 82 % of it reproduced by a depth confound with no immune content, from a
   module that ships default-off.
6. **Anything computed under the shipping delivery law as an absolute.** See §4.

---

## 3. Surface claims that are now stale

Reported only; `landing-3d/` and `app/` are not edited by this file.

### Landing page — `landing-3d/lib/copy.ts`

| Where | Says | Should say |
|---|---|---|
| Section 04 body | *"What these particular macrophages add is about **1.15x**"* | The measured excess is **1.081x and 1.208x** against two nulls — quote **1.08–1.21x**. |
| Guardrail comment, item 2 | *"~**1.13-1.20x**"* | Same. The low end is **1.081**, not 1.13. |
| Section 04 body | *"a null that preserves the macrophages' own clustering returns about **1.7**"* | 1.6581 rounds to 1.7, but the **second** primary null returns **1.8522** and is omitted. Quoting only the weaker null flatters the result. |
| Section 04 | attributes the residual spread to clustering | **82 % of it is reproduced by distance-to-tissue-edge alone, with no immune input** (`immune_sink_impact.md` §5.1a). This is a bigger caveat than the clustering one and is absent. |
| Section 04 | presents the statistic as a product finding | The module **ships default-off and no shipped code path imports it**. |
| Section 02 | *"Interstitial fluid pressure builds where vessels leak…"* | True as mechanism, but the solved IFP is **41x below** measured on breast. If the section shows a solved IFP field, it needs that caveat. |
| Section 05 | *"correcting the transvascular delivery law … changes the top pick on glioblastoma"* | **Correct and well done** — but it does not say the corrected law is **not the shipping default**. See §4. |

Section 04's quintile figures (0.43 / 0.86) and the 2.0x are **correct**. The
guardrails at the top of `copy.ts` are the best-calibrated text on any surface.

### README.md — the most stale surface

| Where | Says | Reality |
|---|---|---|
| Title | *"A **glioblastoma** … simulator"* | Breast is now a co-equal flagship with its own panel, anchors and FAIL. |
| Numerics row | *"Production dx = 12 µm is also the finest grid satisfying the advection CFL at dt = 10 s"* | **Retracted.** CFL 0.663 at dx = 10; floor is dx = 7. And `numerical_rigor.md` no longer contains the word CFL. |
| Numerics row | *"VALIDATED"* | Its own correction block withdraws the headline and the elasticity attributions (N5 above). |
| Correction 1 | *"`therapeutic_coverage` … grid-converged (±0.03)"* | Not established at dx = 12 for either panel (N2). |
| Correction 5 | *"The vascular field is a smooth KDE, not vessels"* | **Superseded.** `_DEFAULT_VASCULAR_MODE = "discrete"` since `ceb96c0`. The whole correction describes a configuration that no longer ships — and the promotion caused a **12x** growth in the oxygen route's downstream effect that the README never mentions (`oxygen_level_shift.md`). |
| Validation table | vascular area fraction *"10.98 % vs measured 2.1 %, 5.2×"* | **Not stale — this is what the scorecard still prints**, because the incommensurable CD34 anchor was never replaced. But it is a known-wrong comparison: against the commensurable Marien 2016 anchor the error is **2.34x**. The README should carry both numbers and say which one `vivum validate` uses. |
| Validation table | *"Transport physics vs measured tissue PK — **VALIDATED out-of-sample**"* | These are `analytic`-scope rows that **never invoke the solver**. Calling them transport-physics validation is precisely the conflation `scorecard.py` warns against: *"'Robust under exchange' and 'never evaluates the delivery law' are very different statements."* |
| Whole document | — | Does not mention IFP 41x, oxygen supply 4,570x, cellular exclusion, the breast docetaxel FAIL, the grid-convergence findings, or the delivery-law default. It points to `SESSION_AUDIT.md` as authoritative; that file is **2026-08-10** and at least nine of its statements are superseded (see `results/SESSION_AUDIT.md` cross-checks in §4). |

### Studio — `app/`

Mostly well-scoped: every "RANKING is invariant" string in `app/server/` is
qualified *"across all three cell-state routes"*, which is a common-mode
perturbation and therefore true (C2). Two gaps:

* **The Studio never mentions the transvascular delivery law.** The landing page
  carries this guardrail and the Studio does not, so a user who reaches the
  product without the marketing site never learns that the top-1 pick is
  law-dependent.
* `app/web/js/mock.js:471` — *"vascular_mode=kde did not converge under grid
  refinement in the held-out check; discrete vessels are the production default."*
  True on the default, but it implies discrete **is** converged. Neither is
  converged at dx = 12.

### Code comment

`tmesim/reconstruct/real_fields.py:146-149` asserts **in the present tense** that
`dt = 10 s` violates the CFL at `dx ≤ 8`. `cfl_resolution.md` §0 (PROVISIONAL)
measures 0.939 and 0.955 on the shipped field. The comment describes a superseded
reconstruction.

---

## 4. The finding that deserves its own section: the shipping delivery law is the known-wrong one

`tmesim/core/config.py:137` — `transvascular_exchange: bool = False`.

The default is the **one-way leak** `S = PS·Cp`, which has no path from
interstitium back to blood and whose steady state is unbounded above plasma. The
codebase's own comment says it is *"measurably wrong"*: it puts 5-FU at **7.80x**
plasma against a measured 0.61x and methotrexate at **0.057x** against a measured
0.60x — **wrong in opposite directions, so no single constant can absorb it.**
Five of eleven GBM drugs sit above plasma; max cellwise AUC/plasma-AUC is **72.6**.

The corrected law `S = PS·(Cp − C)` introduces **no new parameter**, is bounded by
plasma by construction, and moves 5-FU from 12.78x of measured to **~1.33x**. It
is off by default, for a stated and defensible reason — *"turning it on changes
every absolute delivery number this project has ever published, which is a
re-validation exercise, not a flag flip."*

That reasoning is sound as engineering. It is not sound as a claim surface. **Every
absolute number this product currently ships is produced under a delivery law the
codebase itself documents as measurably wrong**, and the corrected law reorders
the GBM panel. Any external communication that quotes a shipped absolute needs to
say this. Today, none does.

---

## 5. Corrections to the working belief set

Recorded because each was believed before this pass and is wrong or imprecise.

| Belief | Correction |
|---|---|
| Docetaxel's breast FAIL is defensible, the project's only falsification instance | **Both halves wrong.** `p99_convention.md` §6.1 downgrades it to **UNDETERMINED at dx = 12** under the shipped convention; it does not survive a 2x grid coarsening, and a sub-grid source removes it. And ~18 other FAIL rows exist. |
| It survived the delivery law at 1.5 % → 3.7 % outside | Two numbers conflated. Delivery law: **1.5 % → 2.6 %** (`breast_kill_revalidation.md` §3). The 3.7 % is the **hindrance** correction, `recal` arm only — it flips **FAIL → PASS** in `frozen`. |
| The +35.2 % normaliser result is in `p99_convention.md` | It is in **`tail90_falsification.md`** §7.1–7.2. `p99_convention.md` gives the tail90 ratio but not that margin. |
| A 3.07x diffusivity change with 99.6 % cancelling | The **input** is 0.325x, not 3.07x; **99.6 % is not printed** in the source; it is dox/epi, `recal`, dx = 12 only. Other scopes give 97.8–98.6 %, and `frozen` gives **anti**-cancellation. |
| Coverage *and ranking* are EXACTLY invariant under oxygen because it is a drug-independent multiplier | **Coverage yes, provably. Ranking no.** Relative kill changes span 71.4–77.0 % — near-common, not common. And oxygen→resistance **is** drug-dependent (`TMZ_HYPOXIA_RESISTANCE = 0.80` vs global 0.35), inert only while the hypoxic band is empty. |
| Spatial spread ~1.13–1.20x | **1.081x and 1.208x** — quote **1.08–1.21x**. Also: p is the 20-draw floor, no CI, one section, and 82 % of the raw ratio is a depth confound. |
| `vivum validate` reports 13 of 13 scored rows never invoke the solver, and the one anchor that does (IFP) fails | **First half true for `--quick` only** (verified live: 13/13, zero PDE-scoped scored rows; full run is 21 of 33). **Second half refuted**: `breast-transvascular` also scores a solved field against an external measurement (human microdialysis) and also fails, in both directions. There are **three** such rows and all three fail. IFP's genuine distinction is being the only PDE anchor with **no fitted coefficients**. |
| Vascular area fraction 2.34x off | **Number correct, framing wrong.** The commensurable anchor was **added, not substituted** — `scorecard.py:1417/:1740` still read the 2.1 % CD34 key, so `vivum validate` still prints **5.23x**, and `ifp_validation.md` §4.2 still cites it. |
| Oxygen supply 4,570x low | **Correct for GBM**; breast is **2,160x**, the source file's own title rounds to ~4,000x, and the ratio is of *conductances* (`K_cap`), not of the `supply` constant itself (~3,900x). The compensating sink is **112–166x**, not exactly 140x. |
| The GBM 3/3 is 2/3 under two *independent* perturbations | Correct in count, but the two 2/3 arms give **bit-identical** TMZ ratios (0.46784) because both are oxygen-**level** perturbations. Independent in provenance, not in mechanism. |
| Oxygen demotion reversed, solve out-predicts baseline by 1.26x, ρ ~ 0.15 | **Confirmed exactly** (0.1523279 / 0.1204202 = 1.26497), within-section, one patient, 816,769 cells. Add: the composite is carried by VEGFA alone — **HIF1A's ρ is +0.0140, the wrong sign** — and the PDE's effective diffusion length is **2–5x too long** vs Thomlinson & Gray. |
| Grid convergence at dx = 12 fails | **Confirmed**, with the numbers in N2. |
| Breast is 5/6 PASS with docetaxel the lone FAIL | **Not what the shipping tool prints.** Live full `vivum validate` returns **3 PASS / 3 FAIL** — 5-fluorouracil and paclitaxel FAIL against zero-width single-row bands the report replaces with a broad-band fallback. Undisclosed anywhere. See note (iii). |

---

## 6. What would most change this document

In rough order of how much each would move the honest assessment:

0. **Reconciling `vivum validate` with `kill_calibration_breast.md`** (note iii).
   The shipping tool and the flagship report disagree about how many drugs fail,
   right now, and that has to be settled before either count is quoted anywhere.
   It is also the cheapest item on this list.
1. **A second section, same tumour type.** Almost every defensible claim above is
   *n* = 1 in tissue. This is the cheapest, largest improvement available.
2. **Any external anchor that scores the solver on the strength of a solve and
   passes.** There is currently one, and it fails. `external_anchor_survey.md` §5
   reports that nothing public measures a spatially resolved free-interstitial
   concentration field in tumour tissue this engine could reconstruct — so this may
   require generating data, not finding it.
3. **Flipping `transvascular_exchange` to True** and re-validating, or publishing
   the caveat on every absolute. The current split is defensible engineering and
   indefensible communication.
4. **Resolving dx.** `L`/dx = 2.4 for the shortest-`L` drugs is not a resolved
   gradient, and it is the direct cause of the only verdict that is a function of
   the grid.
5. **Re-deriving `supply` for the discrete-vessel referent.** A 4,570x error and a
   ~140x compensating error is a coincidence the model should not be relying on.
