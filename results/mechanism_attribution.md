# Mechanistic result visualisation — per-cell attribution

**What this adds.** A run used to return scalars: `kill_fraction`, therapeutic coverage,
mean tumour AUC. Those say *how much* happened. They cannot say **why**, and the kill model
already knows: `tmesim/cells/pkpd.py` computes a local exposure, a Hill term, a resistance
multiplier, a kill probability and a Bernoulli draw for every single cell, and discards all
five one line later. `tmesim/sim/mechanism.py` picks them back up and labels every tumour cell
with the reason it died, or the reason it did not.

It is a **read-out of the model's own arithmetic, not a second model**. `apply_drug_kill` now
takes an optional `diagnostics=` dict and fills it with the exact arrays it decided on — the
same objects, not recomputed copies — and every label is a comparison between two of them.
The run is bit-identical with the attribution on or off, which is asserted rather than claimed
(`tests/test_mechanism.py::test_attribution_is_bit_identical_to_a_run_without_it` compares
`kill_fraction`, `viable_tumor_fraction` and the whole AUC field).

---

## 1. Categories and thresholds

| code | label | limiter | condition |
|---|---|---|---|
| 0 | `KILLED` | — | the engine's own `killed` mask |
| 1 | `SURVIVED_UNREACHED` | delivery | `E/IC50 < unreached_ratio` |
| 2 | `SURVIVED_UNDERDOSED` | dose | `unreached_ratio <= E/IC50 < adequate_ratio` |
| 3 | `SURVIVED_RESISTANT` | resistance | exposure adequate, `m <= resistant_modulation` **and** `m < 1` |
| 4 | `SURVIVED_STOCHASTIC` | chance | both adequate, `p_kill >= stochastic_min_p_kill` |
| 5 | `UNATTRIBUTED` | — | none of the above. **Counted, never reassigned.** |
| 6 | `ALREADY_DEAD` | — | dead before the drug; counts in `kill_fraction`, not killed by it |
| 7 | `NOT_TUMOR` | — | non-tumour cell |

`m = sensitivity * resistance`, the exact multiplicative modulation the kill model applied.

**Two of the four thresholds are not choices, and two more are derived from the run's own
configuration rather than picked:**

| threshold | value | where it comes from |
|---|---|---|
| `adequate_ratio` | 1.0 | the *definition* of IC50 — `E == IC50` is exactly where the Hill term equals half its ceiling |
| `unreached_ratio` | 0.1 | the one conventional constant. At `hill = 1.2` it is 5.9% of the drug's achievable kill probability. Exposed as a parameter so it can be moved and the breakdown watched to move |
| `resistant_modulation` | **derived** 0.5 | `max(pd.hypoxia_resistance, pd.quiescent_resistance)` — the *mildest* resistance penalty the kill model itself applies |
| `stochastic_min_p_kill` | **derived** 0.25 | `f(adequate_ratio) * resistant_modulation` — the lowest `p_kill` the CHANCE branch can produce at `max_kill_rate == 1` |

The last derivation is what makes `UNATTRIBUTED` a live check rather than decoration: under the
production `PDConfig` it is **empty by construction**, and it fills up precisely when something
outside the bands suppressed the kill (a `max_kill_rate` ceiling below 1, a clip, a non-finite
exposure). It was worth having — see §4.

**Precedence.** Exposure bands bind before biology, because below IC50 no restoration of
sensitivity lifts the Hill term past half its ceiling. Cells that are *both* short of drug and
resistant are not hidden: `joint_limited.count` reports how many exposure-limited survivals
would still not have died on dose alone.

---

## 2. Bit-exact reconciliation

`mechanism.reconcile` rebuilds `kill_fraction` from the labels with the **same integer
arithmetic** `CellPopulation.viable_tumor_fraction` uses —
`1.0 - (n_tumor - (KILLED + ALREADY_DEAD)) / n_tumor` — and compares with `==`, not `allclose`.
`simulate_response` **raises** `MechanismReconciliationError` rather than returning a breakdown
that does not add up, and no payload is written for such a run.

Measured on the real Xenium 319 section, `dx = 12 um`, 6 h window, `dt = 10 s`, seed 0,
40,887 measured cells of which 23,812 are TUMOR:

| drug | `kill_fraction` (engine) | rebuilt from labels | `\|delta\|` |
|---|---|---|---|
| temozolomide | `0.01839408701495049` | `0.01839408701495049` | `0.0` |
| cisplatin | `0.43658659499412056` | `0.43658659499412056` | `0.0` |

Identical to the last bit of the mantissa, not to a tolerance.

> **These are the time-resolved-PK numbers, and that matters when comparing files.** The runs
> in this report use `time_resolved = true` with `plasma_peak = 50 µM`, `half_life = 4 h`
> (stated in §3 below). Other files in `results/` report the **static-source** config for the
> same drug and section — cisplatin `0.4361666386695784`, temozolomide `0.01835209138249627`.
> Both pairs were re-run directly from `reconstruct_real('xenium', dx_um=12)` and each
> reproduces to the last bit, so the ~0.1% gap (≈10 of 23,812 tumour cells) is a **config
> difference, not a reproducibility failure**.
>
> Worth naming, because "reconciles bit-exactly" appears throughout this report and could be
> misread: it is a claim that the *label array* rebuilds *this run's* `kill_fraction`. It is
> **not** a claim that two runs of different configs agree. A kill quoted without its PK mode
> is ambiguous.

---

## 3. The finding: temozolomide vs cisplatin on the flagship section

Both runs: `xenium_319`, `dx = 12 um`, 6 h, `plasma_peak = 50 uM`, `half_life = 4 h`,
`sensitivity_mode = uniform`, `time_resolved = true`, seed 0.

| | temozolomide | cisplatin |
|---|---|---|
| compiled IC50 | 606.05 uM | 24.55 uM |
| compiled `bbb_penetrance` | 0.1254 | 0.8577 |
| mean tumour AUC | 7.71e5 uM·s | 1.83e6 uM·s |
| therapeutic coverage 1× | **0.000** | **0.9955** |
| `kill_fraction` | 0.0184 | 0.4366 |
| tumour cells killed | 438 / 23,812 | 10,396 / 23,812 |
| **survival: DELIVERY** | **97.10 %** | 0.23 % |
| **survival: DOSE** | 2.90 % | 0.44 % |
| **survival: RESISTANCE** | 0.00 % | **96.42 %** |
| **survival: CHANCE** | 0.00 % | 2.91 % |
| UNATTRIBUTED | 0 | 0 |
| tumour E/IC50 — p05 / p50 / p95 / max | 0.037 / 0.056 / 0.090 / 0.157 | 1.79 / 3.19 / 5.84 / 14.37 |
| exposure shortfall at the median cell | **17.9×** | 0.31× (already past IC50) |

### 3.1 The expectation was wrong, and it is the model that says so

The brief expected temozolomide to be "overwhelmingly dose- or resistance-limited". **It is
not.** By the documented bands it is 97.1 % DELIVERY-limited, and the reason is not a marginal
threshold call: **no tumour cell anywhere on the section reaches even 16 % of IC50**
(`max E/IC50 = 0.157`), and the median cell is short by a factor of **17.9**. There is no
`unreached_ratio` in `[0.157, 1.0]` that would move a single cell out of the exposure-limited
regime and into RESISTANCE or CHANCE; only the split *between* UNREACHED and UNDERDOSED moves
with that constant, and at this population every cell is 6–27× below IC50.

`therapeutic_coverage_1x = 0.0` is the same statement in the pre-existing scalar. The
attribution is what makes it actionable: it says the 0.0 is not a rounding of "nearly enough",
it is 23,374 individual cells that never saw a pharmacologically meaningful concentration.

Cisplatin is the clean contrast and it inverts completely: coverage 0.9955, median cell at
3.19× IC50, and **96.4 % of its survivors are RESISTANCE-limited**. Delivery is solved; what
stops it is the tumour's own oxygen-driven state.

### 3.2 What "DELIVERY-limited" does and does not mean here

This must be read precisely, because the categories describe *the kill model*, not biology in
general:

* `UNREACHED` and `UNDERDOSED` are **two bands on one axis**, `E/IC50`. The kill model reaches
  `p_hill` through that ratio alone, so it **has no separate potency term** and cannot
  distinguish "too little drug arrived" from "the drug is too weak". Neither can this
  attribution, and it says so in every payload (`exposure_axis_note`).
* What *is* separable within a run: IC50 is a single per-drug constant, so the entire
  cell-to-cell **spread** in `E/IC50` is delivery, while the population's **position** relative
  to 1.0 is set jointly by delivery and potency. Temozolomide's 17.9× shortfall decomposes
  upstream into a 6.8× `bbb_penetrance` disadvantage (0.125 vs 0.858) and a 24.7× IC50
  disadvantage (606 vs 24.6 uM) against cisplatin — but that decomposition comes from the drug
  compilation, not from the per-cell kill, and this module does not assert it.

So the honest headline is: **temozolomide fails on this tissue because exposure never comes
within an order of magnitude of its potency**, and the delivery/dose split is a statement about
*how far short* it fell, not about which upstream term caused the shortfall.

### 3.3 The resistance state is nearly universal, and it is the section's known weak point

`joint_limited` for temozolomide is **21,073 of 23,374** exposure-limited survivors (90.2 %) —
that many are *also* at or below the resistance threshold. Raising the dose alone would not
have killed them. The same ratio for cisplatin is 98.9 % of its (tiny) exposure-limited set.

That is the engine's internal oxygen re-banding calling almost every tumour cell
HYPOXIC or QUIESCENT, and it lands directly on this section's documented limitation: the
319-gene panel **measures no hypoxia**, so no state route clears the bar
(`results/cell_state_impact.md` §5) and `kill_fraction` on this section carries an unresolved
uncertainty of up to 0.3252 absolute purely from the choice of state route. Cisplatin's
"96.4 % RESISTANCE-limited" therefore inherits that uncertainty in full, and the run's own
caveats say so. **The cross-drug contrast (delivery-limited vs resistance-limited) is the
defensible output; the absolute size of the resistance bucket is not.**

---

## 4. What the attribution refuses to do

`UNATTRIBUTED` was 0 on both real runs, as the derived thresholds guarantee under the
production `PDConfig`. Two guards exist so that it is a real check and not a formality:

* **A ceiling the bands do not model is not called "bad luck".** Drop `max_kill_rate` to 0.4
  and a cell at exactly the adequate boundary with no resistance has `p_kill = 0.2`, below the
  derived CHANCE floor of 0.25. It is labelled `UNATTRIBUTED` and counted, not filed under
  CHANCE. (`test_unexplained_survival_is_unattributed_not_reassigned`)
* **A defect found by this work.** `resistant_modulation` derives from the run's `PDConfig`, so
  the **counterfactual resistance-off arm** (`resistance_off: true`, used by the rank
  attribution) derives a threshold of 1.0 — under which `m <= threshold` is true for *every*
  cell, and 100 % of survivors in the one run where resistance was switched OFF would have been
  reported as RESISTANCE-limited with an empty CHANCE bucket: the exact inversion of the truth.
  `is_resistance_limited` therefore requires `m < 1.0` as well — "something was actually
  suppressed" — and the degenerate case is pinned by
  `test_resistance_off_arm_does_not_file_every_survivor_as_resistant`. The complementary test
  checks the guard does not evict QUIESCENT cells, which sit at `m == 0.5` exactly on the
  threshold.

---

## 5. Time-lapse

`PDEConfig.time_resolved` (`results/dose_schedule.md`) lets tissue concentration genuinely rise
**and fall**, so a time-lapse can show the drug clear rather than plateau.
`simulate_drug_transport(..., capture_frames=N)` keeps N snapshots of the **concentration**
field `C` — not the AUC field, which only accumulates and so can never show drug leaving —
evenly spaced over the window and always including `t = 0` and `t = total_time`. Capture is a
`copy` at selected steps: no arithmetic changes, and
`test_frame_capture_does_not_perturb_the_solve` asserts the AUC and final-concentration fields
are bit-identical with capture on and off.

Measured on the two runs above (24 frames, 939.1 s interval, 6 h window, 424×675 grid,
27.5 MB payload each), peak tissue concentration by frame:

```
temozolomide  0 → 144.3 → 142.4 → 136.4 → … → 57.8 → 55.2 uM
cisplatin     0 → 555.0 → 530.5 → 507.0 → … → 214.8 → 205.3 uM
```

Drug enters, spreads, and then clears — a 2.6× fall for temozolomide and a 2.7× fall for
cisplatin across the window. On the **default static** plasma source the same series rises
monotonically toward a quasi-steady state instead (temozolomide peaks at 96.1 uM at `t = T`),
because the delivery rate is held at its window average by construction. The run says so:
`RunResult.timelapse.plateau_note` carries that sentence whenever `time_resolved` is false, and
`/api/validate` warns before the run is submitted.

---

## 6. Reproduce

```bash
./.venv/bin/python -m app.server          # then POST the spec below, or:
./.venv/bin/python -m pytest tests/test_mechanism.py -q
```

```jsonc
{ "section": "xenium_319", "dx_um": 12.0, "drug": "temozolomide",
  "sensitivity_mode": "uniform", "vascular_mode": "discrete",
  "oxygen_mode": "pde", "state_mode": "oxygen", "mgmt": null,
  "pk": {"plasma_peak": 50.0, "half_life_s": 14400.0,
         "n_doses": 1, "dose_interval_s": 86400.0},
  "total_time_s": 21600.0,
  "attribute_mechanism": true, "timelapse_frames": 24, "time_resolved": true }
```

`GET /api/runs/{id}/mechanism` (TMSM, 45.8 kB for 40,887 cells) and
`GET /api/runs/{id}/timelapse` (TMSL, 27.5 MB for 24 frames) serve the two payloads; the wire
formats are frozen in `API_CONTRACT.md`.

**Not for clinical use.** This tissue is an untreated resection with no outcome label. Every
number here is a prediction of the physical model and has never been scored against a measured
treatment response.
