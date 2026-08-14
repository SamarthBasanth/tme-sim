# The two-route transport model is the default, and the claim is narrower than "fixed"

**Status: landed.** `tmesim/molecular/compile_drug.py` now applies the two-route
(paracellular + transcellular) cell term on every compile. Every
`DrugParams.diffusivity` in the project moved, the perivascular-sink constant
`alpha` moved with it, one external-anchor verdict moved **FAIL → PASS**, one
external-anchor reading moved the other way, and **0 of 12** breast-panel
verdicts moved — which was the registered prediction and was verified rather
than assumed.

Suite: **961 → 967 passed** (six net new tests; none removed). Eight assertions
fired and were re-pinned or inverted, none deleted; two bitwise reference cases were regenerated with their superseded
values kept. Nothing was fitted to the anchor dataset.

**The one-sentence claim.** Against 71 compounds measured by one method in one
cell line, blind and unfitted, the shipping transport model reaches
**R² = 0.718** — and still reproduces only about **24 %** of the measured spread,
from **one cell line**, and still cannot be trusted for absolute numbers.

---

## 0. What changed, in one table

| | superseded | shipping | factor |
|---|---:|---:|---:|
| `two_route` in `compile_drug` | absent | **ON** | — |
| `HBD_LOGP_DECREMENT` | 1.0 (convention) | **1.34** (sourced) | — |
| `calibrate_alpha()` | 148.60815765445545 | **22.223477711172563** | 0.14954 |
| doxorubicin `D` (µm²/s) | 184.74861430191334 | **27.92558262333494** | 0.15115 |
| Pruijn blind **R²** | 0.022 | **0.718** | 32.6× |
| Pruijn modelled spread | 2.3× | **8.3×** | (measured 29.7×) |
| Pruijn geometric bias | 3.40× high | **1.03× high** | — |
| Pruijn within 3× | 38 % | **99 %** | — |
| exact-mass pairs ordered | **0 of 31** | **24 of 31** | — |
| MCL hindrance factor | 1.7732-fold, **FAIL** | **10.3274-fold, PASS** | — |
| MCL absolute `D_eff` (TPZ) | 287.32 µm²/s, FAIL | **49.33 µm²/s, FAIL** | direction flipped |
| MCL k reading vs measured | 0.49–0.82× (would pass) | **0.080–0.134× (would fail)** | — |
| breast-panel verdicts moved | — | **0 of 12** | — |
| `kill_fraction`, real section `dx = 12` | 0.023980 | **0.066227** | 2.76× |
| `kill_fraction`, real section `dx = 80` | 0.026373 | **0.139384** | 5.29× |

Efflux, intracellular retardation, `physics/subgrid.py` and the vessel
`geometry` default are all **unchanged and still off**. §7 says why for each.

---

## 1. The decision, and the four gates it had to clear

`results/two_route_transport.md` §7 recommended *"MERGE the physics, do NOT flip
the default"* and gave four reasons. `results/two_route_efflux.md` §5.1 upgraded
that to *"YES, flip — but do the retardation term first"*.
`results/intracellular_retardation.md` then built the retardation term and
**refuted** it. That closes the sequence:

| gate | closed by | outcome |
|---|---|---|
| The one unsourced parameter sits at its own sweep optimum | `two_route_efflux.md` §3 | **CLOSED.** Sourced at 1.34 from open-access Abraham LSER coefficients; it is *not* the optimum and scores marginally **worse** on R². |
| Efflux is unresolved and paclitaxel is the largest number in the blast radius | `two_route_efflux.md` §2, §4.1 | **CLOSED.** HT-29 is functionally P-gp-negative (Wils 1994, PMID 7945455); 0 of 71 compounds carry an efflux assignment; and at the sourced decrement the paclitaxel outlier falls 4.94× → 1.60× *without any efflux term*. |
| The `K > 1` retardation term must be built and scored first | `intracellular_retardation.md` | **CLOSED, by refutation.** Built, verified against a direct transient solve to 0.01 %, and it makes the model decisively worse: R² 0.718 → 0.606, coverage 99 % → 0 %, matched-MW rank agreement 77.4 % → 51.6 % (a coin flip). Pruijn's own Discussion says why: their `Dmcl` is a **free-drug** coefficient. |
| Five of six worst residuals are lipophilic over-predictions, held only by the logP clamp | `two_route_efflux.md` §4.2 | **HALF FALSIFIED, HALF IMPROVED.** The clamp binds on **0 of 68** compounds, so it holds nothing; and at 1.34 the residuals become two-sided (3 over / 3 under) with coverage 96 % → 99 %. |

### 1.1 The evidence was re-verified from source before anything was changed

Not taken from the documents. Re-run first, on the unmodified tree:

| check | required | observed |
|---|---|---|
| `experiments/two_route_transport.py --stage prereg` | passes | **passes**, all 7 parameters `OK` |
| `calibrate_alpha()` with the flag off | `148.60815765445545` | **exactly equal** |
| `--stage hbd` at decrement 1.34 | R² = 0.718 | **0.718**, Spearman +0.821, geo bias 1.03×, within-3× 99 % |
| the sweep optimum is elsewhere | at 1.0 | **at 1.0, R² = 0.725** |

The stop condition was: if any of those failed to reproduce, stop and report.
None did.

---

## 2. What was actually changed

### 2.1 The wiring — one branch in `build_profile`

```python
interstitial_D = micro_D * mesh_h          # unchanged: Brownian × Ogston-Renkin
if two_route:                              # NEW, default True
    tr = two_route.two_route_diffusivity(
        interstitial_D, radius_nm, descriptors["logp"], descriptors["hbd"], ...)
    effective_D = interstitial_D * tr["hindrance"]
else:
    effective_D = interstitial_D           # bitwise the old value
```

`MolecularProfile` gained `interstitial_diffusivity` (always the pre-two-route
number, so no caller has to recompile to get it) and `two_route` (the law's own
output dict). Nothing else in the compile changed. `two_route=False` recovers
the superseded value **bitwise**, and that is asserted, not asserted-in-prose:
`tests/test_two_route.py::TestShippingDefault::test_the_superseded_value_is_still_exactly_recoverable`.

### 2.2 The pre-registration was EXTENDED, not bypassed

`results/sweep/tworoute_prereg.json` was frozen before the model was ever scored.
It is **untouched**. `two_route.HBD_LOGP_DECREMENT_PREREGISTERED = 1.0` keeps
live code checkable against it, and `--stage prereg` still asserts all seven
original parameters.

The promoted value has its own registration,
`results/sweep/tworoute_promotion_prereg.json`, which records the derivation
(3.621 × 0.37 = 1.340), both coefficient sources verbatim, the one remaining
modelling choice, the direction of its bias, and the terms deliberately left off.
`--stage prereg` asserts live code against **both** files, and
`test_the_promoted_parameter_matches_its_own_registration` checks the arithmetic
rather than the prose.

### 2.3 The test that existed to force this re-run did its job

`tests/test_two_route.py` asserted that `compile_drug.py` contained no reference
to `two_route` at all, *deliberately*, so that anyone wiring it in would be
forced to re-run the blast radius before the suite went green. It fired. It was
**inverted, not deleted** — `test_two_route_is_wired_into_compile_drug` now pins
the new contract (wired in, on by default, still switchable off).

---

## 3. The new alpha, and what actually moved it

**`calibrate_alpha()` = 22.223477711172563**, from 148.60815765445545.

`alpha` is a **gauge coordinate** (`results/alpha_identifiability.md`): it
absorbs a change at the calibration anchor and re-exports it as a global
multiplier. A change in its value is expected and is not by itself a problem.
What matters is that it moved for the stated reason, by the stated arithmetic.

**Doxorubicin has 6 hydrogen-bond donors, so its membrane is shut.** Its
transcellular gain is **1.0000002473**, i.e. the *new* physics contributes
2.5 × 10⁻⁷ of the change. What moved doxorubicin is the **impermeable-cell
Maxwell term** — the `P → 0` limit of the same law, `(1−φ)/(1+φ/2)` at
φ = 0.7892 — which is `0.15115445288971746`. The two-route factor is
`0.15115449027238376`. The difference between them is the entire transcellular
contribution at the anchor.

So: **the compound-specific physics that justifies the promotion is invisible at
the anchor.** The anchor moved by the compound-*independent* cell-exclusion term
that `results/cellular_exclusion.md` already analysed. What the two-route model
adds is the **differential** across the panel, and `alpha` cannot absorb that.

The alpha ratio is **0.149544**, not the D ratio **0.151154**, and the gap is
arithmetic rather than physical:

```
alpha = (D_eff·(1−0.6·φ_ECM)/L*²  −  k_passive) / uptake_proxy
```

`k_passive` (doxorubicin's clearance, `6.308e-5 /s`) is the genuinely-slow term
and is **not** amplified by `alpha`, so it is subtracted before the division and
does not scale. It was 0.19 % of `k_target` before and is 1.25 % now. That is the
whole discrepancy, and it is pinned in
`test_alpha_follows_D_through_its_own_closed_form`.

`L = sqrt(D_eff/k_eff)` still lands in the measured Primeau band by construction —
the fit target `64.92 µm` (PMID 16361566) did not move.

---

## 4. Regenerated: the blind score of the SHIPPING engine

`experiments/diffusivity_falsification.py`'s pre-flight now runs both ways, and
the second half is the wiring check:

```
PRE-FLIGHT 1: is the FALSIFIED model still a pure function of MW?
  tirapazamine from SMILES / from MW only : 287.3216 / 287.3483 um^2/s  (0.01%)
  -> descriptors do NOT reach it.  CONFIRMED
PRE-FLIGHT 2: does the SHIPPING model now use more than MW?
  tirapazamine, shipping / same MW, no logP, no donors : 49.3322 / 216.4827
  -> descriptors DO reach D_eff.  CONFIRMED
```

`experiments/two_route_transport.py --stage pruijn`, at the shipping decrement,
68 scoreable compounds:

| | shipping engine (superseded) | + impermeable cells | **shipping (two-route)** | measured |
|---|---:|---:|---:|---:|
| spread | 2.3× | 2.3× | **8.3×** | **29.7×** |
| Pearson r (log₁₀) | +0.149 | +0.149 | **+0.848** | — |
| **R²** | 0.022 | 0.022 | **0.718** | — |
| Spearman | +0.083 | +0.083 | **+0.821** | — |
| geometric bias | 3.40× high | 0.51× low | **1.03× high** | — |
| within 3× | 38 % | 68 % | **99 %** | — |
| ratio min / median / max | 0.76 / 3.65 / 17.84 | 0.11 / 0.55 / 2.70 | **0.23 / 1.03 / 2.99** | — |

**Matched molecular weight**, where a MW-only model is flat by construction:

| MW | n | measured fold | superseded | **shipping** |
|---:|---:|---:|---:|---:|
| 192 | 3 | 1.4× | 1.00× | **1.03×** |
| 213 | 3 | 2.5× | 1.00× | **1.10×** |
| 263 | 4 | 3.7× | 1.00× | **1.22×** |
| 289 | 3 | 1.4× | 1.00× | **1.68×** |
| 305 | 4 | 3.8× | 1.00× | **2.69×** |

**24 of 31** exact-mass pairs ordered correctly (77 %), against **0 of 31**.

**Tirapazamine**, the project's registered anchor, compiled from its own SMILES:

| | `D` (µm²/s) | vs measured 40.4 |
|---|---:|---|
| superseded | 287.4 | **7.1× over** |
| paracellular only | 43.4 | 1.08× |
| **shipping** | **44.4** | **1.10×** |

**φ_cell sensitivity, the widest-uncertainty input.** R² stays 0.683–0.727 across
0.40–0.85; the adopted 0.7892 is neither the R²-optimal (0.60) nor the
bias-optimal (0.85) value.

---

## 5. Every published number that changed

### 5.1 Compiled diffusivities, and the two-route factor per drug

`ecm_reference = 0.4`. The `vs dox` column is what survives the α gauge.

| drug | MW | logP | HBD | log K_mem | P (µm/s) | `D` superseded | **`D` shipping** | factor | gain | **vs dox** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| doxorubicin | 543.5 | 0.001 | 6 | −8.04 | 9.6e-7 | 184.7486 | **27.9256** | 0.15115 | 1.0000 | **1.0000** |
| epirubicin | 543.5 | 0.001 | 6 | −8.04 | 9.6e-7 | 184.7486 | **27.9256** | 0.15115 | 1.0000 | 1.0000 |
| docetaxel | 807.9 | 3.26 | 5 | −3.44 | 0.0335 | 156.2800 | **23.8623** | 0.15269 | 1.0102 | 1.0102 |
| 5-fluorouracil | 130.1 | −0.798 | 2 | −3.48 | 0.0565 | 323.8095 | **49.3493** | 0.15240 | 1.0083 | 1.0083 |
| cisplatin | 298.0 | 0.195 | 2 | −2.49 | 0.4214 | 236.1844 | **38.6696** | 0.16373 | 1.0832 | 1.0832 |
| **paclitaxel** | 853.9 | 3.736 | 4 | −1.62 | 2.155 | 151.3981 | **36.5813** | 0.24162 | **1.5985** | **1.5985** |
| temozolomide | 194.2 | −2.078 | 1 | −3.48 | 0.0567 | 277.9851 | **42.4245** | 0.15261 | 1.0097 | 1.0097 |

**Panel differential: 1.599×**, down from 4.939× at the pre-registered
convention. Sourcing the decrement is what shrank it: 4 donors × 1.34 = 5.36
subtracted from paclitaxel's Crippen logP 3.736 closes its membrane
(`P` 49.4 → 2.15 µm/s). The largest single number in the old blast radius was
substantially an artefact of an under-sourced donor penalty.

### 5.2 The breast panel — 0 of 12 verdicts move, verified

Two arms, per the convention of `results/mechanism_attribution.md` and
`results/vessel_profile.md`. Published `dx = 12 µm` grid, 458 × 628, 167,780
cells / 89,785 tumour, same `dt`, same 6 h window, same seed, same `PDConfig`,
same p99 bath-matching convention, same frozen literature bands. Each drug at
**its own** shipping `D`.

**Arm A — gauge held.** `alpha` frozen at the published 148.6082, only `D` moves.
Solved in `results/two_route_efflux.md` §4.3 and carried over: **1 of 6 moves,
docetaxel FAIL → PASS.** This arm silently breaks the measured Primeau anchor and
is not what the codebase runs; it is reported because it isolates the effect of
the `D` change from the gauge.

**Arm B — `recal`, which is what ships.** `alpha` re-solved on doxorubicin =
22.223477711172563, applied to every drug. This is the true wiring, and it is the
row `two_route_efflux.md` §5.2 Table B could only approximate.

| drug | ratio published | **ratio shipping** | band | verdict published | **verdict shipping** | moved? |
|---|---:|---:|---|---|---|:--:|
| doxorubicin | 0.2301 | **0.2330** | 0.0575 – 1.0 | PASS | **PASS** | no |
| paclitaxel | 0.1503 | **0.1747** | — | NOT TESTABLE | **NOT TESTABLE** | no |
| docetaxel | 0.1484 | **0.1520** | 0.0152 – 0.1463 | FAIL | **FAIL** | no |
| epirubicin | 0.2301 | **0.2330** | 0.0110 – 1.0 | PASS | **PASS** | no |
| 5-fluorouracil | 0.3379 | **0.3355** | — | NOT TESTABLE | **NOT TESTABLE** | no |
| cisplatin | 0.3085 | **0.3105** | 0.2086 – 1.5004 | PASS | **PASS** | no |

**0 of 6 scored rows move; 0 of 12 across both arms' scored points.** The
registered prediction in `two_route_efflux.md` §4.3 was *"0 of 12 in `recal`"* and
it landed. `frozen` and `recal` are numerically identical here (α 22.2235 vs
22.2234, a 5e-6 difference from the two code paths' rounding), which is the check
that `alpha_for_scale(1.0) == calibrate_alpha()`.

**Docetaxel's FAIL is deeper, not rescued.** Its ratio goes 0.1484 → 0.1520
against a band top edge of 0.1463: from 1.5 % outside to **3.9 % outside**. A
correction that deepens the panel's only failure is not a correction chosen to
flatter it.

**Paclitaxel moves most** (+16 %), which is expected — it is the only drug with a
materially different `vs dox` factor — and it is NOT TESTABLE in both arms, so it
carries no verdict either way.

### 5.3 The external transport anchors — one moved to PASS, one moved the other way

`vivum validate`, group `mcl-transport-coefficients`, tirapazamine.

| row | superseded | **shipping** | measured | verdict |
|---|---:|---:|---|---|
| **hindrance factor** `D_free/D_eff` | 1.7732-fold | **10.3274-fold** | 12-fold (V79-171b MCL, PMID 9806526) | **FAIL → PASS** |
| **absolute `D_eff`** | 287.32 µm²/s | **49.33 µm²/s** | 74–130 µm²/s (PMID 9806526, 17869669) | FAIL → FAIL, **direction flipped** |
| **consumption rate** `k`, calibrated | 4.2311e-3 /s | **6.9182e-4 /s** | 5.167e-3 – 8.667e-3 /s | NOT TESTABLE → NOT TESTABLE, **reading reversed** |

**The hindrance PASS was not bought by moving a bar.** The band is the repo's own
loosest declared factor (`pk_reality.ORDER_OF_MAGNITUDE_FACTOR` = 3×, i.e.
4–36-fold), it was chosen while the engine was failing it, and it has not been
touched. The pass is still 0.86× the central measured value. Pinned in
`test_mcl_verdicts_are_scored_on_the_engine_as_it_ships`, which asserts the band
is still exactly `(4.0, 36.0)`.

**The absolute row's flip matters more than it looks.** The engine was 2.21×
*above* the measured band; it is now 1.50× *below* it. Overshooting in the other
direction is progress in magnitude and is not a pass. Read against the hindrance
row, which now passes, the pattern is the signature
`results/cellular_exclusion.md` §4.2 named: **the hindrance ratio is now the
right size and the thing being hindered is not** — i.e. a residual free-diffusivity
error the two-route work does not touch.

**The consumption-rate row's disclosure reversed sign, and this is the most
dangerous thing in the change.** That row declines to score, and it said so with:
*"scored on k_calibrated this row would PASS the repo's 3× band. Declining to
score REMOVES a pass from the tally rather than removing a failure."* At the new
`alpha` that is **false**: k_calibrated is now 0.080–0.134× the measurement, i.e.
7.5–12.5× too small, outside the band. Declining to score now removes a **failure**.
A disclosure whose whole purpose is to say which way an omission cuts, and which
has quietly started cutting the other way, is worse than no disclosure. It is
corrected in `scorecard.py` and pinned by
`test_mcl_consumption_rate_declines_to_score_and_says_which_way`, which now
asserts the string `"removes a FAILURE"`.

**Note the coupling.** The Primeau fit pins the *ratio* `D/k`, so lowering `D`
lowers the fitted `k` by the same factor. The D-anchor improved and the k-anchor
degraded by roughly the same amount. Two independent external measurements now
**disagree** about where the engine sits. That is unresolved and is §8's sharpest
open item.

### 5.4 Absolute engine outputs move by 3–5×

Two bitwise reference cases in `tests/data/static_source_reference.json`:

| case | quantity | superseded | **shipping** | factor |
|---|---|---:|---:|---:|
| real Xenium, `dx = 80`, TMZ ×3 | `kill_fraction` | 0.026373 | **0.139384** | 5.29× |
| production, `dx = 12`, TMZ ×3 | `kill_fraction` | 0.023980 | **0.066227** | 2.76× |
| production, `dx = 12` | `penetration_index` | 0.495919 | **0.425139** | 0.857× |
| production, `dx = 12` | `plasma_mean` | — | — | **bit-identical** |

**Diffusivity fell 6.6× and exposure ROSE. That is the surprising part and it is
not a bug.** The Primeau anchor pins `D/k`, so the fitted loss rate
`k = clearance + alpha·binding` fell with `D`. A weaker sink leaves far more drug
in tissue even though it spreads more slowly.

`results/cellular_exclusion.md` §6 predicted exactly this asymmetry: the
p99-normalised ratio the breast panel scores cancels a common-mode exposure change
(0 of 12 verdicts moved), while **threshold statistics such as kill fraction and
IC50 coverage do not**. Anyone quoting an absolute `kill_fraction` from a run
predating this change is quoting a number that has moved by 3–5×. Absolute kill
fraction was already documented as a non-defensible output for independent
reasons, and this is a second one.

`plasma_mean` — pure PK, no spatial dependence, no diffusivity in it — is
bit-identical, which is the positive control on the mechanism.

---

## 6. The re-pins, and why the new values are right rather than merely different

Eight assertions fired. **None was deleted.** Every one was re-pinned with its
superseded value preserved, following the pattern of
`results/tumor_sigma_units.md` §5.

| test | what fired | disposition |
|---|---|---|
| `test_two_route.py::test_compile_drug_does_not_import_two_route` | the guard that existed to force this re-run | **inverted** into `test_two_route_is_wired_into_compile_drug`; new tests assert the old value is still bitwise recoverable and that the re-pin *equals the old pin times the law* |
| `test_two_route.py::test_alpha_is_unchanged` | 148.6082 → 22.2235 | re-pinned; a new test derives the new value from `calibration`'s own closed form and asserts the alpha ratio is **not** the D ratio |
| `test_two_route.py` diffusivity pin | 184.7486 → 27.9256 | re-pinned |
| `test_reality_layer.py::test_alpha_calibration_hits_anchor` | band `50 < α < 2000` | widened at the bottom to `5 < α`, with the reason recorded: the band guards the *order of magnitude of the diagnosis* (the structure-only binding prior is still >10× too small), and the assertion that actually pins correctness — that `L` lands in the measured band — is untouched |
| `test_app_server.py::test_result_carries_mandatory_provenance` | α in the Studio provenance | re-pinned to 22.2235 |
| `test_byo_compound.py::test_diffusivity_is_a_function_of_mw_alone` | the contract **inverted** | renamed `test_diffusivity_is_no_longer_a_function_of_mw_alone`; it now asserts isomers DIFFER, and a new test asserts a matched-mass pair with a real donor gap separates *in the physically required direction*. The honest caveat is pinned too: for that particular isomer pair the separation is ~6e-6, because Crippen barely distinguishes them |
| `test_byo_compound.py::test_the_diffusivity_falsification_is_disclosed_at_the_point_of_use` | the disclosure it pins became false | renamed `test_the_transport_limit_is_disclosed_at_the_point_of_use`; it now pins **all three halves** of the narrower claim — the new score, the old score it replaced, and the shortfall that remains — plus that the panel still carries a genuine hard stop, so softening the transport caveat did not soften the renderer's stop behaviour |
| `test_validate_cli.py` × 2 | the MCL verdict and the reversed disclosure | re-pinned; the band is asserted **unchanged** at (4.0, 36.0) so the PASS cannot have been bought by widening it |
| `test_time_resolved_pk.py` × 2 | bitwise references | regenerated by hand; superseded values kept under `_superseded_2026-08-13_mw_only_diffusivity`; cause, justification and failure pattern recorded in the JSON's `_repin_history` and in the test module docstring |

**Why the new bitwise values are right and not merely different.** The old values
encode a diffusivity law whose only solute input was molecular weight — a law
that was scored blind against 71 measured compounds and refuted at R² = 0.022,
and that was flat by construction exactly where the measurement moves most. The
new values are produced by a law scoring R² = 0.718 on the same data, which
brings the engine's dimensionless multicellular-layer hindrance inside its
measured band for the first time. These cases now run a transport coefficient
that has been tested externally and largely agrees with the measurement; they did
not before. That is a claim about the *coefficient*, not about the outputs: §5.4
is explicit that the outputs are not thereby correct.

**The failure pattern confirms the mechanism rather than correlating with it.**
Only a case whose drug comes from the molecular compile can move. Exactly those
did: `real_xenium_dx80_tmz_3dose` and `production_xenium_dx12_tmz_3dose`.
`default_1dose`, `weak_1dose`, `weak_3dose`, `weak_10dose_noconv` and
`solver_flat` all use hand-built `DrugParams` literals with no compile in the
path and are **byte-for-byte unchanged**;
`test_production_dx12_ifp_solve_matches_stored_reference_bitwise`, which has no
drug in it at all, passed untouched. This was the stated stop condition: had a
literal-drug case or the IFP pin moved, the mechanism would not have been what we
thought and the change would not have landed.

---

## 7. What was deliberately NOT changed, and why

| | state | why |
|---|---|---|
| **efflux** | OFF | `results/two_route_efflux.md`. HT-29 is functionally P-gp-negative (Wils 1994, PMID 7945455, verbatim: *"No functional expression of P-glycoprotein was found in the differentiated HT29-18-C1 cell line"*). 0 of 71 Pruijn compounds carry an efflux assignment — pinned by a test — so the term would change the blind score by **exactly zero**. And no primary numeric efflux ratio for paclitaxel, the one drug it was needed for, could be retrieved. |
| **intracellular retardation** | OFF | `results/intracellular_retardation.md` refuted it against this anchor: R² 0.718 → 0.606, geometric bias 1.03× → 0.098×, within-3× 99 % → 0 %, matched-MW rank agreement 77.4 % → **51.6 %**. Pruijn's `Dmcl` is a **free-drug** coefficient and their own Discussion says so. |
| **`tmesim/physics/subgrid.py`** | OFF, not bundled | `results/area_conservation_regeneration.md` found that once the vessel field is continuous the closure stops being a no-op and inflates the F_state span **325×**. Its default-off status is now load-bearing. |
| **vessel `geometry='area'`** | at its current default | Separately evidenced. Two consequential flips in one change makes attribution impossible if anything breaks. |
| **`OxygenConfig.supply`, band edges, `DX_UM`** | untouched | Out of scope by instruction and by prudence. |
| **`results/sweep/tworoute_prereg.json`** | untouched | It is a pre-registration. Editing it would destroy the only thing it is for. |

---

## 8. What is still not true about the transport model

This section is the point of the document. The claim that replaces
*"FALSIFIED, and the fix is default-off"* is **narrower than "fixed"**.

1. **It reproduces about a quarter of the measured spread.** 8.3× modelled
   against 29.7× measured. The model *under-separates* a congeneric series
   systematically — by roughly 4× — rather than separating it correctly. A user
   ranking similar-mass analogues should read the order as a hypothesis worth
   testing, not as a measurement. The shortfall was **predicted before the
   scoring run** (`tworoute_prereg.json`, registered prediction): the
   transcellular route saturates at a ceiling set by the measured
   cytoplasm/interstitium conductivity contrast, so a bounded model cannot
   reproduce an unbounded-looking spread.

2. **The entire test is ONE cell line.** Every one of the 71 measurements is
   HT-29. Pruijn's own paper reports its asymptote coefficients moving up to
   **100-fold** across SiHa, HCT8-Ea and HCT8-Ra. A model claiming to be physics
   rather than a fit should be scored on that movement, and it has not been.
   **This is the cheapest unbought evidence in the project.**

3. **Absolute numbers did not become trustworthy, and one anchor got worse.**
   Tirapazamine's compiled `D_eff` still misses its measured multicellular-layer
   band, now 1.50× low instead of 2.21× high. The calibrated consumption rate
   moved from within 2× of an independently measured rate constant to 7.5–12.5×
   below it. The two external anchors are coupled through the Primeau `D/k` fit
   and **now disagree about where the engine sits**. Nothing here resolves that,
   and it is the sharpest open question the promotion leaves.

4. **One modelling choice remains inside the sourced parameter.** 1.34 = 3.621 ×
   0.37 is derived from two measured coefficient sets, but 0.37 is the acidity of
   a *typical* donor group and the model consumes an integer donor **count**.
   Across the whole measured range of that choice R² stays 0.588–0.722, so
   nothing bets on it — but the choice has not disappeared. And 1.34 is an
   **upper** estimate for a bilayer: hexadecane is the harshest available
   membrane proxy, and no verified Abraham coefficient set for a
   DMPC/liposome-water system could be obtained.

5. **The descriptors feeding it are still wrong for at least one panel drug.**
   Crippen puts cisplatin at logP +0.20 against an experimental ≈ −2.2. The
   membrane leg inherits that error linearly, and at a decrement of 1.34 it is
   *more* sensitive to donor count than before, not less.

6. **The model is known-wrong for taxanes in a specific, measured way.** Choi
   2011 (PMID 21225235) measures paclitaxel at 38 % penetration through an HT-29
   multicellular layer at 1 µM against 100 % at 10 µM, in a P-gp-negative line.
   That is saturable intracellular binding. The model has no capacity term; the
   term built for it was refuted (§7); paclitaxel is on the shipped panel; and it
   is the drug whose diffusivity moved most (+60 % relative to the anchor).
   Paclitaxel is NOT TESTABLE on the panel, which limits the damage but does not
   remove it.

7. **The taxane grid-resolution hazard is real at this decrement.** Closing
   paclitaxel's membrane drops both taxanes' boundary layer per grid cell:
   docetaxel to `L/dx` ≈ 0.94 at `dx = 12 µm`, i.e. below one cell. No verdict
   turns on it (docetaxel is FAIL at every grid tested and this pushes it further
   out), but the taxane numbers are quoted on a grid that does not resolve the
   gradient producing them. A grid-convergence check at these scales was **not**
   run here and is owed.

8. **Every unresolved item in `results/cellular_exclusion.md` §7 is still
   unresolved**, and this model composes on top of all of them — in particular
   the `brownian.py` 300 nm obstacle field / `ALPHA_TORTUOSITY = 0.6` possible
   double-count, which sets the matrix conductivity the two-route law consumes.

9. **Nothing here says anything about a solved concentration field, and nothing
   here says anything about clinical response.** This scores one PDE coefficient.
   `results/external_anchor_survey.md` §5 stands unchanged: no public dataset
   measures a spatially resolved free-interstitial concentration field in tumour
   tissue that this engine could be scored against.

---

## 9. Claim surfaces updated

Every surface that said the transport model is FALSIFIED for similar-mass
compounds and that the fix is default-off and wired into nothing. Leaving any of
them would have been as bad as the overstatement they replaced.

| surface | change |
|---|---|
| `README.md` | the "read this before you rank a compound series" box, the external-anchor scorecard (two rows), the honest summary, and the "what would most improve this" list — item 0 is now *score a second cell line* |
| `landing-3d/lib/copy.ts` | guardrail 1b rewritten, keeping the superseded version verbatim so copy cannot drift back to **either** extreme; section 03 body copy |
| `app/server/report.py` | limit `diffusivity_mw_only_falsified` → `diffusivity_two_route_partial`, still **severity critical**, still positioned as limit 2c |
| `app/server/api.py` | the self-compiled-compound warning |
| `tmesim/cli.py` | CSV export: `diffusivity_is_falsified` → `diffusivity_model`, and `what_is_defensible` |
| `tmesim/molecular/byo.py` | the module honesty table, `WILL_NOT_TELL_YOU`, the Studio `diffusivity` field, and the hard-stop caveat |
| `tmesim/validation/scorecard.py` | the `mcl-transport-coefficients` docstring and all three rows, including the reversed disclosure |
| `app/web/js/ui/provenance.js`, `app/web/js/ui/spec.js` | the Studio's provenance panel and the self-compiled-compound tooltip |
| `docs/DATA_ACCESS_PLAN.md` | the "what was falsified" row, now recording that the falsified law was replaced rather than that it stands |
| `docs/QUICKSTART.md`, `docs/architecture.md`, `docs/MATHEMATICAL_MODEL.md` (§2.5 and the limiting-uncertainty list), `docs/molecular_layer.md` (§3 block, §4 grade table, §4 postscript, §5), `docs/REAL_BUILD_PLAN.md` | same substitution throughout |
| `tests/data/static_source_reference.json` | re-pin + `_repin_history` + superseded key |

`results/diffusivity_falsification.md`, `results/two_route_transport.md`,
`results/two_route_efflux.md`, `results/intracellular_retardation.md` and
`results/cellular_exclusion.md` are **left as written**. They are the historical
record of how this was established, every number in them is still exactly
reproducible with `two_route=False`, and this document is the landing record.

---

## 10. Reproducing this

```
./.venv/bin/python experiments/two_route_transport.py --stage prereg   # both registrations
./.venv/bin/python experiments/two_route_transport.py --stage limits   # the algebra
./.venv/bin/python experiments/two_route_transport.py --stage pruijn   # the blind score
./.venv/bin/python experiments/diffusivity_falsification.py            # the wiring pre-flight
./.venv/bin/python -m pytest tests/test_two_route.py -q                # the pins
./.venv/bin/python -m pytest -q                                        # 967 passed
```

The breast panel at the shipping configuration is
`dtissue_sensitivity.run(dx_um=12.0, modes=("frozen", "recal"),
explicit_points=[reference_transport_constants()["D_tissue_cm2_s"]])` — one
reference point, because `build_drug` scales each drug's **own** compiled `D`, so
scale = 1 is exactly "every drug at its own shipping diffusivity, at the new
alpha". 12 solves, ~10 minutes.
