# Zero-width bands: a verdict rendered against a bar nothing can pass

**Measured and fixed 2026-08-12.** Found by `results/CLAIM_INVENTORY.md` note
(iii), which ran the shipping tool instead of reading about it.

**One line.** `vivum validate` scored the breast 3-D:2-D kill panel **3 PASS /
3 FAIL** while `results/kill_calibration_breast.md` published **5 PASS / 1
FAIL**; two of the scorecard's FAILs were rendered against **zero-width bands**
that no number can satisfy. Neither tally was right. The corrected panel is
**3 PASS / 1 FAIL / 6 NOT TESTABLE**, and the two paths can no longer disagree
because they no longer both decide.

---

## 1. Reproduction

`./.venv/bin/python -m tmesim.cli validate --only breast-kill-3d2d`, shipping
physics (`transvascular_exchange=False`), warm cache, before any change:

```
SUMMARY   PASS 3   FAIL 3   NOT TESTABLE 4   DEFERRED 0   (of 10 rows)
```

| drug | predicted ratio | band the scorecard scored against | verdict |
|---|---:|---|---|
| docetaxel | 0.1484 | 0.0152–0.1463 (folds 8.3–58x, n = 2) | FAIL |
| **5-fluorouracil** | 0.3379 | **0.1636–0.1636** (fold 7.5–7.5x, n = 1) | **FAIL** |
| **paclitaxel** | 0.1503 | **0.0412–0.0412** (fold 25–25x, n = 1) | **FAIL** |

It reproduced exactly as reported, to the printed digit. The two flagged bands
have **zero width**: their lower and upper edges are the same float. The set of
predictions that satisfies `lo <= x <= hi` when `lo == hi` has measure zero on a
continuous quantity, so the FAIL was guaranteed before the model ran.

## 2. Mechanism — confirmed, not assumed

The chain, top to bottom:

1. `experiments/kill_calibration_breast_report.DRUG_LIT_SCOPE` maps both
   paclitaxel and 5-fluorouracil to `{"kind": "drug_field"}`, i.e. *filter the
   shared `LIT_2D3D` table by its `drug` field*.
2. Counting that table by `drug` field gives: temozolomide 6, cisplatin 5,
   doxorubicin 5, **5-fluorouracil 1**, **paclitaxel 1**, mitomycin C 1,
   trabectedin 1. Both drugs have exactly **one** matched row —
   PMID 33395433 (Koch, colorectal, fold 7.5x) and PMID 41422245 (Mutahar,
   breast, fold 25x).
3. `kill_calibration.literature_bridge` reports the **range** of the matched
   folds, deliberately (`fold_min` / `fold_max`). Over one value those are the
   same value.
4. `kill_calibration.bridge_to_kill_ratio_band` maps `fold_min` → `ratio_hi`
   and `fold_max` → `ratio_lo`. Equal folds in, equal ratios out:
   `0.1636–0.1636` and `0.0412–0.0412`.
5. `kill_calibration.ratio_test` then asked `lo <= predicted <= hi` and, as it
   had to, said FAIL.

So the mechanism is exactly the n = 1 collapse `kill_calibration_breast.md` §3
already suspected: **a min-max over a single point is a point.** Nothing in the
transport model, the tissue, the potency or the grid is implicated.

### The two code paths that disagreed

`build_ratio_section` computed **both** candidate verdicts for every drug —
`test` against the drug-matched band and `test_solid` against the broad
cross-tumour-type `solid` band — and then left the *choice between them* to each
consumer:

* **the document generator** (`render_markdown`) picked `test_solid` whenever
  `band_is_drug_matched` was false, i.e. whenever n = 1 → paclitaxel and
  5-fluorouracil **PASS** → 5/6;
* **the scorecard** (`scorecard._run_breast_kill_3d2d`) picked `test`
  unconditionally whenever the drug had any drug-matched literature at all →
  the same two drugs **FAIL** → 3/6.

Same input object, same function, same run — two selection rules, written
twice, three years of documents quoting the one that was never the shipping
tool's. The defect is not that one of the two rules was wrong. It is that the
**selection was duplicated at all**.

## 3. Which behaviour is correct — and it is neither of them

**The scorecard's FAIL is indefensible.** A bar nothing can pass is the exact
mirror of the defect for which the GBM wick-in-needle IFP row is already
reported NOT TESTABLE: *"A verdict that cannot be failed is not a verdict"*
(`scorecard._run_ifp_wick_in_needle`). A verdict that cannot be passed is not
one either.

**The document's fallback is also wrong, and this is the harder call, because
the fallback is superficially reasonable and it is what four other experiment
drivers in this repo do.** Four reasons it does not survive:

1. **It is band-widening in disguise.** The drug's own evidence cannot test it,
   so a **62x-wide** band (0.4–25.0x, n = 20) is substituted, and both drugs
   pass *because of the substitution*. `results/vascular_segmentation_audit.md`
   §2.1 declined to widen a band it had just shown to be statistically wrong,
   for precisely this reason — *"a band that widens until the model fits inside
   it is exactly the rigging this project forbids"*, and re-scoping it is *"a
   reviewer's decision, not an agent's"*. The same line holds here, and this
   change does not re-scope anything: it declines to score.
2. **It contradicts this panel's own treatment of the same band.** The `solid`
   band is exactly what cyclophosphamide, tamoxifen, vinorelbine and palbociclib
   are shown against — and they are refused a verdict on it, in the scorecard's
   own words, "deliberately NOT scored against the broad, non-drug-matched
   band". One band cannot be unscorable for four drugs in a table and scorable
   for two others in the same table.
3. **It is partly circular.** Both single rows carry `scope=('solid',)` — they
   are *members* of the band being substituted. Worse for paclitaxel: its 25x
   point **is** the band's `fold_max`, i.e. its only datum defines the edge it
   is then tested against.
4. **It carries almost no information.** The solid band spans 0.4–25.0x, most
   of the published range of 3-D:2-D effects for cytotoxics. Passing it is close
   to unfalsifiable — the same weakness, in the opposite direction, that
   disqualified the IFP anchor.

**Decision: a single published study is a point estimate, not an interval, so a
drug with n = 1 drug-matched rows is NOT TESTABLE** — reported with the point
quoted, with the broad band quoted explicitly as context and explicitly not as a
verdict, and never scored. This is the same status the no-literature drugs get,
which is correct: in both cases the honest statement is *"no interval exists to
test against"*, and the difference between them (one measurement versus none) is
recorded in the reason string, not in the verdict.

Note that this decision **removes two PASSes** the project has been quoting and
adds no FAILs. It is not the reading that maximises passes; the reading that
maximises passes is the one being retired.

## 4. The fix — the two paths can no longer disagree

The defect was duplicated decision-making, so the fix is a single decision, not
a patched scorecard.

* **`kill_calibration.score_ratio_row`** (new) is now the **only** place that
  decides which band a 3-D:2-D row is scored against and whether it can be
  scored at all. It returns one object carrying `status`, `scored`, `band`,
  `test`, `n_rows`, `reason`, plus `context_band` / `context_test` that are
  computed on every row and are never the scored band.
* **`experiments/kill_calibration_breast_report.build_ratio_section`** attaches
  that object to each drug as `row["scored"]`. Both consumers — the markdown
  generator that writes `results/kill_calibration_breast.md` and
  `scorecard._run_breast_kill_3d2d` that `vivum validate` prints — now *render*
  `row["scored"]["status"]`. Neither derives a verdict. The legacy keys are kept
  but are themselves derived from the shared decision, so they cannot drift
  either.
* **`kill_calibration.ratio_test`** returns `NOT TESTABLE` with a stated reason,
  never `FAIL`, for any degenerate band — equal, crossed or non-finite edges —
  via the new `band_is_degenerate`. This is the belt: the guard lives in the
  only function that renders a 3-D:2-D verdict, so no caller anywhere in the
  repo can score against a zero-width band by forgetting to check.
* **`scorecard._kill_rows`** re-checks degeneracy and downgrades to NOT TESTABLE
  with a reason. This is the braces, and it covers the **GBM** path too, which
  has the same latent hazard (a drug with one `drug:` -scoped row would collapse
  identically) even though no GBM band is degenerate today.
* **`tests/test_validate_cli.py`** gains `TestZeroWidthBands` and
  `TestBreastKillPanelCannotDiverge`. The divergence guard is data-free: it
  feeds a synthetic result to the real `build_ratio_section`, renders the
  DOCUMENT's table cell via the newly extracted `ratio_table_row`, and asserts
  it equals the status `vivum validate` assigns for every drug in the panel. It
  also asserts no scored scorecard row carries a degenerate band, and — the
  point of the hard rule — that a single-study row stays NOT TESTABLE *even when
  the broad band would have passed it*.

### The generator hazard was live

`results/kill_calibration_breast.md` **is generated** by
`experiments/kill_calibration_breast_report.py`. Correcting the `.md` by hand
would have been reverted by the next re-run — the failure mode this project has
hit three times (immune-sink, oxygen-demotion, kill_calibration). The fix is
therefore in the generator, and a re-run now *emits* the corrected verdicts. The
checked-in `.md` was deliberately **not** edited by this work and is currently
stale; regenerating it is what it needs:

```
./.venv/bin/python experiments/kill_calibration_breast_report.py
```

Verified by running the generator with `OUT_MD` redirected to a scratch path. Its
§3 table, its headline block and its §1 literature bullets all now print
NOT TESTABLE for paclitaxel and 5-fluorouracil, and its verdict line reads
**"Of the 4 formally-testable drugs, 3/4 land inside their comparison band"**,
which is the scorecard's tally exactly.

## 5. Corrected tallies

Shipping physics (`transvascular_exchange=False`), dx = 12 µm, warm cache.

| tier | before | after |
|---|---|---|
| `--quick` | PASS 7, FAIL 6, NOT TESTABLE 4, DEFERRED 8 (25 rows) | **unchanged** — no kill anchor runs in quick |
| full | PASS 13, FAIL 20, NOT TESTABLE 17 (50 rows) | **PASS 13, FAIL 18, NOT TESTABLE 19** (50 rows) |

Two rows moved **FAIL → NOT TESTABLE**. No row moved to PASS. Scored rows in the
full tier drop from 33 to 31; the "21 of 33 scored rows never invoke the solver"
figure in `CLAIM_INVENTORY.md` §0 becomes **21 of 31**.

Per group, full tier:

| group | PASS | FAIL | NOT TESTABLE |
|---|---:|---:|---:|
| gbm-kill-3d2d | 3 | 0 | 0 |
| **breast-kill-3d2d** | **3** | **1** | **6** |

**The honest breast panel verdict: 3 PASS / 1 FAIL / 6 NOT TESTABLE.** Not 5/1,
not 3/3. Of ten panel drugs only four are testable at all: doxorubicin PASS,
cisplatin PASS, epirubicin PASS, **docetaxel FAIL**. Docetaxel is therefore
still the panel's only FAIL — but the sentence "breast is 5/6 PASS" is retired,
and so is "the shipping tool prints three FAILs".

## 6. What this does not fix

* **GBM is untouched: 3/3 PASS**, and every caveat on it stands
  (`kill_calibration.md` §3.1 measures 2/3, 2/3, 3/3, 3/3, 1/3 across five
  alternative configurations).
* **Docetaxel's FAIL is unchanged and remains UNDETERMINED at dx = 12 µm** under
  `p99_convention.md` §6.1. This work says nothing about it.
* **Four older experiment drivers still carry the retired fallback rule**
  in their own local copies of the selection logic —
  `experiments/dtissue_sensitivity.py`, `experiments/grid_convergence_breast.py`,
  `experiments/delivery_source_geometry.py` and `experiments/p99_convention.py`
  all compute `formal = t_specific if band_is_drug_matched else t_solid`
  (`experiments/oxygen_joint_grid_test.py` records the same flag but runs the
  GBM panel only, where every band has n >= 5, so it is unaffected). They were left
  alone deliberately: their published documents are frozen measurements, and
  editing the driver without regenerating the document would manufacture exactly
  the generator/document mismatch this file exists to close. **Any paclitaxel or
  5-fluorouracil "verdict" in those five documents is a broad-band substitution
  and is not a verdict under the shipping rule.** Their docetaxel numbers — the
  drug those documents are actually about — are unaffected, docetaxel having
  n = 2.
* **A second study for either drug would make it testable again.** The threshold
  is `kill_calibration.MIN_ROWS_FOR_BAND = 2`, and
  `test_single_study_drugs_are_not_testable_not_failed` pins the current n = 1
  so that new literature forces the tally to be revisited rather than silently
  absorbed.
