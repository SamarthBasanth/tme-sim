# Export honesty — what the two artefacts a user keeps now say

`results/first_user_readiness.md` walked a researcher's whole journey and found
one inversion worse than the rest:

> A researcher who compares the README to the exported report will find the
> report is the MORE CONFIDENT of the two.

The two things a user leaves with are the **generated report** (`report.html`,
the self-contained document) and the **exported CSV** (`vivum compare-drugs
--csv`). This file records what was wrong in those two artefacts, what changed,
and — as importantly — what did not.

Every defect below was **re-verified at source before it was fixed**. All six
were real. One number quoted in the readiness report was not (§7).

---

## 1. The report asserted a grid-convergence claim its own documentation retracts

**Verified.** `app/server/report.py`, limit `resolution`, shipped this on every
report:

> "This report ran at dx = 12 um. kill_fraction is grid-converged for dx <= 20
> um; penetration_index is NOT converged above ~15 um and is over-estimated on
> coarse grids ... (0.285 at dx = 60 um falling to 0.187 at dx = 10 um)."

at `"severity": "info"` — the lowest the renderer has — citing
`results/numerical_rigor.md` §2. Both halves are contradicted by the sources:

* That file's own **CORRECTION** block withdraws the headline ("*the mesh is not
  the bottleneck once dx ≲ 12 µm*"), records that its §4 elasticities were
  attributed to the wrong knobs, and states that **the §2 prose had the sign of
  the grid trend backwards**.
* The §2 table shows the penetration index **rising** with refinement — 0.1110
  at dx = 60 µm to 0.3062 at dx = 10 µm — so a coarse mesh **under**-states it.
  The quoted pair 0.285 / 0.187 appears **nowhere** in that column (0.2850 is
  the `kill_fraction` at dx = 15).
* `results/grid_convergence_breast.md`: docetaxel scores **0.1484** at the
  published dx = 12 µm against a Richardson dx → 0 estimate of **0.1537**, and
  its **verdict flips FAIL → PASS between dx = 12 and dx = 16 µm**.
* `results/grid_convergence_gbm.md`: the scored ratio is **non-monotone in dx**
  for all three drugs tested — not an asymptotic regime, so Richardson and
  power-law fits are outside their assumptions (observed orders **p = 0.96** and
  **p = −3.33**, GCI **5.44 %** and **6.42 %** on the finest pairs, both
  power-law fits flagged **DEGENERATE**).

**Changed.** The limit is not deleted — the dx a report ran at is exactly what a
limits block exists to state. It now says there is **no dx at which this
pipeline has been shown to be grid-converged**, quotes the retracted sentence as
history rather than as an assertion, carries the two contradicting ladders with
their numbers, states the corrected direction of the penetration trend, and is
**no longer filed as `info`** (it is `warn`, and its title reads *"Resolution —
the grid-convergence claim is WITHDRAWN"*).

Guarded against over-correcting the other way: the limit also says the solver
itself is **not** impeached — mass is conserved to ~4e-14 on the real irregular
mask and the prediction is dt-converged to <1e-2 for dt ≲ 20 s
(`numerical_rigor.md` §1 and §3, untouched by the correction).

The same retracted claim propagated to two more places in the Studio backend:

* `app/server/api.py` — `default_spec()`'s docstring called dx = 12 µm "the
  converged resolution"; the two `validate_spec` dx warnings asserted
  convergence at 20 µm and had the penetration direction backwards ("typically
  15-25 % high"). All three corrected. The constants `DX_KILL_CONVERGED_UM` /
  `DX_PENETRATION_CONVERGED_UM` are renamed to `DX_KILL_LADDER_COARSE_UM` /
  `DX_PENETRATION_LADDER_COARSE_UM` — they are the coarsest rungs a ladder was
  ever run at, which is all they were entitled to be — with the old names kept
  as aliases so nothing importing them breaks.
* **`app/server/jobs.py` still carries it**, in the per-run caveat at ~L707
  ("penetration_index is ... OVER-ESTIMATED on grids coarser than ...",
  "kill_fraction is grid-converged only for dx <= 20 um"). That file is not
  owned by this pass and four other agents are live in it, so it is **reported,
  not edited**. It is the last instance known to this pass.

## 2. The report hero named a single winner the same page calls law-dependent

**Verified.** The hero read *"Highest predicted kill on this specimen:
cisplatin"* and *"The ORDER is the claim"*, while limit 2b, added the day
before, correctly states the ranking is **LAW-DEPENDENT** — and cisplatin is
precisely the pick the corrected delivery law overturns (cisplatin →
vincristine, `results/transvascular_exchange_revalidation.md` §4.1).

**Changed — the headline is NOT deleted.** Ranking is the defensible output;
that is the project's claim. What was wrong was stating a winner as settled. The
hero now reads:

> **Top of this ranking, under `S = PS*Cp`: cisplatin** — `0.3699` · UPPER BOUND
> · LAW-DEPENDENT

and carries a second paragraph naming the alternative source term, the
permeability re-weighting (47× across this drug list), the measured top-1 change
**cisplatin → vincristine** (4 of 45 pairs swapped), and
`results/ranking_robustness.md` §8: *"the reportable output of this pipeline for
the top of the panel is a **SET**, not a winner"* — with P(cisplatin = rank 1) =
0.931 under that report's declared priors, falling to 0.587 in its stress
scenario. It closes: *this is the top of the ranking one law and one set of
priors produced, not this specimen's drug.*

The law named in the hero is the one **stored with that document** (read from
`PDEConfig` when its runs were made), falling back to the shipped default for a
document that predates the field. Nothing about the law is hard-coded in the
renderer.

The pre-existing honest lines ("a model output, not a recommendation", "the
ORDER is the claim; the LEVEL is not", the UPPER BOUND chip) are untouched.

## 3. A retired metric had a column; its replacement was computed and dropped

**Verified.** `penetration_index = mean(AUC)/p99(AUC)` is a **scale-invariant
shape statistic**: multiply the whole exposure field by any constant and it does
not move, so on real tissue it reports the section's vascular geometry, which
every drug shares. `results/delivery_metrics.md` §4 retired it for that reason —
it spans 0.3044–0.5107 across a panel whose delivered exposure spans 158×, and
its separation (spread / grid drift) is **15.2** against **67.7** for
`therapeutic_coverage` @1× IC50. Both exported artefacts demonstrated its
failure mode on their own front page: **temozolomide carries the highest
penetration index of the GBM panel while ranking last.** Meanwhile
`therapeutic_coverage_1x` was computed into every report row
(`report.py` ~L502), stored in `report.json`, and **never rendered**.

**Changed.** In the report's ranked panel and in `compare-drugs` (printed table
and CSV alike), the retired index has lost its column to `coverage @1x IC50`,
which was already there. Both surfaces explain the swap in place, with the
citation. Coverage is not exempt from anything else on the page: the report says
in the same breath that it is a delivery number, so it is produced under the
delivery law named above and carries the same mesh term as everything else.

The index is not carried under a renamed column in the CSV. It is retired; the
CSV is the artefact most likely to become a slide, and re-exporting a retired
metric under any name re-creates the defect. The provenance header says
explicitly that it is absent and why.

## 4. The exported CSV carried no provenance at all

**Verified.** Eight columns, no header comment, no sidecar. Absent: the section,
`dx`, the exposure window, `alpha`, the delivery law, the sensitivity mode, and
every word of "not for clinical use". Present at full float64 precision:
`kill_fraction`, the quantity the README says is not defensible.

**Changed.** Every CSV now opens with a `#` provenance block — inside the file,
so it cannot be separated from the numbers by copying it — carrying:

`not_for_clinical_use` · `generator` · `section_source` · `sample_path` ·
`patient_id` · `library` · `n_drugs` · `dx_um` · `grid` · `n_cells` ·
`n_tumor` · `exposure_window_s` · `plasma_peak_uM_override` · `seed` ·
`alpha` + `alpha_source` (the single doxorubicin anchor, PMID 16361566) ·
`oxygen_mode` · `state_mode` · `sensitivity_mode` +
`sensitivity_is_upper_bound` · `transvascular_exchange` · `delivery_law` +
`delivery_law_caveat` · `retired_metric_note` · `what_is_defensible`.

Every value is **read from the run** — the reconstruction report, the config
actually used, the calibrated alpha — never restated from memory.

The same block is written as a `<file>.provenance.json` sidecar for readers that
cannot skip comments, and the command prints how to read the file
(`pandas.read_csv(path, comment='#')`). The one cost: a naive `csv.DictReader`
that does not skip `#` lines now needs one line of filtering. That was judged
worth it — provenance that lives in a sidecar only is provenance that gets left
behind on the first copy-paste — and the tests read it the way a user would.

## 5. The delivery-law disclosure reached every Studio surface but not `compare-drugs`

**Verified.** `app/server/delivery_law.py` is used by `jobs.py`, `projects.py`,
`report.py`, `inspect.js`, `projects.js`, `provenance.js` and `reports.js`. It
was not used by `tmesim/cli.py` — the command the README calls the core use
case, and whose entire output is a ranking.

**Changed.** `cmd_compare_drugs` prints the disclosure before it simulates
anything, and repeats the one-line form under the table. It **imports the one
implementation** rather than restating the claim: that module is a runtime read
of `PDEConfig`, so it cannot go stale, and a second copy of the sentence is
exactly the thing that would drift. It imports nothing but `typing` at module
scope, so this does not pull in the optional `[app]` extra; if the `app` package
is genuinely absent the command says so loudly and still prints the raw
`PDEConfig.transvascular_exchange` value, rather than printing nothing.

## 6. Three different default exposure windows, and the ranking command printed none

**Verified.** CLI 1,800 s (`tmesim/cli.py::_FAST_TOTAL_TIME`), Studio 21,600 s
(`app/server/api.py::default_spec`), library 43,200 s
(`core.config.PDEConfig.total_time`). `cmd_simulate` and `cmd_molecular_run`
print their window; `cmd_compare_drugs` did not.

**Changed — printed, not harmonised.** The windows are deliberately not
silently unified: raising the CLI to the library default multiplies every
quickstart run's cost ~24×, and lowering the Studio's would change numbers it
has already published. Instead:

* `cmd_compare_drugs` prints `EXPOSURE WINDOW = <n> s (<h> h)`, says whether
  that came from its own default or from `--total-time`, names the library
  default read live from `PDEConfig`, and states that exposure, coverage and
  kill all scale with the window so two surfaces are not comparable unless the
  windows match. The window is also in the CSV provenance block.
* The Studio's window is now the named constant
  `app/server/api.py::STUDIO_DEFAULT_TOTAL_TIME_S`, documented at its definition
  as one of three, so the divergence is stated in code rather than discovered.
  The report's meta block already printed the window it used; that is unchanged.

## 7. One thing in the readiness report that did not check out

`results/first_user_readiness.md` states that correcting the delivery law makes
`P(cisplatin rank 1)` fall **0.9277 → 0.3965**, citing
`results/ranking_robustness.md`. **Those two numbers are not in that file.** It
reports P(cisplatin = rank 1) = **0.931** (90 % CI 0.929–0.934) under its §3
priors, and **0.587** in the stress scenario C. Its §8 conclusion — *"the
reportable output ... for the top of the panel is a SET, not a winner"* — does
stand, and so does the law-driven top-1 change (cisplatin → vincristine), which
is measured in `results/transvascular_exchange_revalidation.md` §4.1 rather than
in the robustness report. The hero quotes the numbers that exist, from the files
that contain them. The finding itself was correct; only that one citation was
not.

---

## What was deliberately NOT changed

* **No existing disclosure was weakened.** The delivery-law caveat, the CAVEATS
  panel hard stops, the provenance chips, the NOT FOR CLINICAL USE banners and
  the `top-ranked is LAW-DEPENDENT` line in the Studio's report list are all
  untouched. `app/web/js/ui/reports.js` needed no edit: it already carried the
  law, and there is no CSV export on the web surface.
* **No new absolute number.** Coverage replaces a retired metric; it does not
  add an absolute kill anywhere, and the UPPER BOUND chip and the interval limit
  still travel with the kill that was already rendered.
* **`penetration_index` is still stored** in `report.json` and in the run
  metrics. It is data; what was wrong was rendering it as a headline column.
* **`tmesim/reconstruct/`, `tmesim/molecular/`, `tmesim/validation/`,
  `tmesim/physics/`, `landing-3d/`, `README.md` and every existing results file**
  were not touched.
