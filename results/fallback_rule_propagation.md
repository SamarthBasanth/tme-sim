# The retired fallback rule had four more copies, and all four had published verdicts

*Baseline at completion: `pytest -q` = **775 passed** (the baseline moved from 772 to
775 under this work when the `_TUMOR_SIGMA` units fix landed — see §3.5). Written after
commit `2466ce1` ("Zero-width bands: neither the tool nor the document was
right") fixed the single-study band defect in `tmesim/validation/kill_calibration.py` and
the breast report generator, and disclosed — rather than silently patched — that four
older experiment drivers still carried the retired rule. This file confirms that
disclosure at source, routes all four through the shared decision, regenerates each of
their documents from the fixed driver, and records exactly what moved. **No PDE was
solved**, which is also why two of the four documents now carry a stated staleness
banner (§3.5).*

> ### Headline
>
> **The disclosure was accurate and complete. All four drivers genuinely carried the
> retired rule, all four had it reach a published document, and none of them was already
> correct.** In every case the affected drugs are the same two — breast **paclitaxel**
> (one drug-matched study, 25.0x, n = 1) and breast **5-fluorouracil** (one drug-matched
> row, 7.5x, n = 1) — whose min-max "band" is a zero-width point. Every GBM verdict is
> untouched (temozolomide n = 6, cisplatin n = 5, doxorubicin n = 5), as are docetaxel
> (n = 2), epirubicin (n = 2), doxorubicin (n = 5) and cisplatin (n = 5) on breast.
>
> **No ratio changed anywhere.** The band decision happens strictly after transport, so
> every regeneration re-scored a stored scalar against a frozen band. **100 verdict cells
> moved from `PASS` to `NOT TESTABLE`** across the four drivers' artefacts (dtissue sweep
> 40, dtissue confirm 8, dtissue refcheck 2, grid ladder 14, delivery blast 6, p99
> conventions 28, p99 bias 2); not one `PASS`, `FAIL` or numeric value moved otherwise.
>
> **One document's stated conclusion did rest on a substituted verdict**, and correcting
> it made that document's finding **worse**, not better:
> `delivery_source_geometry.md` claimed *"only four of the nine keep a margin above
> 80 %"* and named paclitaxel and 5-fluorouracil among the four survivors. Both
> "margins" were distances to the broad non-drug-matched context band. The true figure
> is **two of the seven scored verdicts**.
>
> **The document the brief flagged as the real risk — `p99_convention.md`, which ranks
> normalising conventions by how many verdicts each reproduces — survives intact, and
> for a checkable reason.** Its PASS-count row drops by **exactly 2 in every one of the
> twelve columns**, because paclitaxel and 5-fluorouracil scored PASS against the broad
> band under *every* convention. They were constant columns that discriminated between no
> two conventions. The ranking, the gaps, the `tail90` recommendation and the explicit
> rejection of `tail99` are all unchanged. The old table was wrong; it was not wrong in a
> way that was doing any work.

---

## 1. Confirming the disclosure at source

### 1.1 The rule, as each driver carried it

The retired rule is: *if the drug-matched literature band is not usable, score the row
against the broad cross-tumour-type `solid` band instead.* Since a drug-matched band was
deemed "usable" only when `n > 1`, the two single-study drugs fell through to a band
spanning **0.4x-25.0x fold** (62x wide, kill-ratio 0.0412-1.5004), which they both
passed comfortably.

| driver | line (pre-fix) | how it was spelled | reached a document? |
|---|---:|---|---|
| `experiments/dtissue_sensitivity.py` | 300 | `formal = t_specific if b["band_is_drug_matched"] else t_solid` | **YES** — `results/dtissue_sensitivity.md` |
| `experiments/grid_convergence_breast.py` | 204 | `formal = t_specific if b["band_is_drug_matched"] else t_solid` | **YES** — `results/grid_convergence_breast.md` |
| `experiments/delivery_source_geometry.py` | 687 | `formal = t_spec if b["band_is_drug_matched"] else t_solid` | **YES** — `results/delivery_source_geometry.md` |
| `experiments/p99_convention.py` | 771 | `used = b["band"] if b["band_is_drug_matched"] else b["band_solid"]` | **YES** — `results/p99_convention.md` |

The fourth is spelled differently — it substitutes the *band* rather than the *verdict* —
but it is the same decision and it produces the same substitution. **None of the four was
already correct**, so none was left alone.

`p99_convention.py` additionally reached the rule a second way: its `breast_bands()`
imported `build_bands()` from `dtissue_sensitivity.py` and then re-applied the
substitution itself. Fixing `build_bands` alone would not have fixed it.

### 1.2 Which drugs the rule actually touched

Read live from the frozen literature tables:

| drug | n | fold range | band status |
|---|---:|---|---|
| doxorubicin | 5 | 1.0-18.8x | interval — scored |
| cisplatin | 5 | 0.4-6.0x | interval — scored |
| docetaxel | 2 | 8.3-58.0x | interval — scored |
| epirubicin | 2 | 1.0-76.0x | interval — scored |
| **paclitaxel** | **1** | **25.0x (point)** | **zero-width — NOT TESTABLE** |
| **5-fluorouracil** | **1** | **7.5x (point)** | **zero-width — NOT TESTABLE** |
| temozolomide (GBM scope) | 6 | 2.0-5.0x | interval — scored |

The disclosure's claim that "their docetaxel numbers are unaffected (n = 2)" is
confirmed, and extends to epirubicin, doxorubicin, cisplatin and all three GBM drugs.

### 1.3 A related path that was checked and found already correct

`delivery_source_geometry._gbm_bands` and `p99_convention.gbm_bands` both carried a
*different* fallback — `if not b.get("n"): b = literature_bridge(scope="solid")`, i.e.
substitute the broad band when a drug has **no** drug-matched literature at all. On the
GBM panel this branch is dead: all three GBM drugs have n >= 5. It has been routed
through the shared decision anyway (where the same case is `NOT TESTABLE`), so it cannot
become live silently, but **no GBM verdict moved and no GBM number changed.**

---

## 2. What changed in the code

The decision is now imported, never re-derived. Each driver builds its bands, carries the
**bridges** alongside the mapped bands (the bridge is what knows `n`), and hands both to
`tmesim.validation.kill_calibration.score_ratio_row`, which returns the verdict. Nothing
under `tmesim/` was touched.

* `dtissue_sensitivity.py` — `build_bands` carries `specific_bridge` / `solid_bridge`;
  `score()` calls `score_ratio_row`; new `verdict_fields()` renders that one decision into
  the record schema; new `rescore_rows()` re-derives stored verdicts.
* `grid_convergence_breast.py` — same three changes. Its Richardson/power-law
  extrapolated verdicts now inherit the decision too: extrapolating an unscored drug's
  ratio to `dx -> 0` does not manufacture a band for it.
* `delivery_source_geometry.py` — `_run_panel` calls `score_ratio_row` for **both**
  tissues; new `rescore_panel()`; rows now persist `ic50_uM` / `hill` so a cached panel
  can be re-scored without a compile.
* `p99_convention.py` — `gbm_bands` / `breast_bands` return band *records* (band, context
  band, both bridges); `stage_d` and `stage_d_bias` both score through `_p99_score`, a
  one-line wrapper over `score_ratio_row`, so published, bias-corrected and
  convention-swept verdicts all come out of the same call.

**Stale caches cannot republish the retired rule.** Every driver re-derives stored
verdicts before any analysis reads them — on the resume path, on `--report-only`, and on
cached-panel replay. A checkpoint written last week is re-scored on load.

### 2.1 Why this cost no PDE solves, despite `tail90_falsification.md`

`results/tail90_falsification.md` §1 established that the tissue ladders' **per-cell
exposure vectors** were never written to disk — only summary statistics and a
sum-and-first fingerprint — so anything needing a conditional tail statistic at a
non-published rung requires a full re-solve. That was checked before starting, and
re-confirmed: `results/sweep/gridconv_breast_e_unit.npz` does not exist.

It does not bind here. **The band decision consumes one scalar per row, `ratio_exact`,
and that scalar is persisted at every rung, every sweep point and every arm.** The band
itself depends only on the frozen literature tables and on `(ic50, hill)`, neither of
which `dx`, `D_tissue` or the source geometry touches — every driver already says so in
its own `build_bands` docstring, and the regeneration below is the test of that claim.
Re-scoring is therefore exact:

| document | how it regenerated | solves |
|---|---|---:|
| `dtissue_sensitivity.md` | `--report-only` from `dtissue_sweep_dx24.json` | **0** |
| `grid_convergence_breast.md` | `--report-only` from `gridconv_records.jsonl` | **0** |
| `delivery_source_geometry.md` | `--stage blast` replayed all 6 panels from `delivery_source_geometry_blast_dx12.jsonl`; document hand-updated from the regenerated JSON | **0** |
| `p99_convention.md` | `--stage D` from the published run caches; document hand-updated from the regenerated `p99_verdicts.json` | **0** |

Two of the four documents (`delivery_source_geometry.md`, `p99_convention.md`) are
hand-written prose over machine-generated JSON rather than machine-emitted markdown. For
those, "regenerated from the fixed driver" means the JSON was regenerated by the fixed
driver and every affected figure in the prose was replaced from it. No number in either
document was typed from anywhere else.

---

## 3. What moved, document by document

**100 verdict cells moved, all `PASS -> NOT TESTABLE`, all paclitaxel or
5-fluorouracil. Zero ratios moved. Zero PASS/FAIL verdicts moved.**

### 3.1 `results/dtissue_sensitivity.md` — 50 cells; conclusion unaffected

Paclitaxel and 5-fluorouracil each carried a verdict at all 10 sweep points in both
`frozen` and `recal` modes (40 cells), plus the §7.1 zero-solve control against the
published cache (2) and the §7.3 published-grid confirmation set (8).

* §0 bottom line: *"6 verdict(s) PASS everywhere"* -> *"across the 4 of 6 swept drugs that
  carry a scorable band ... 4 verdict(s) PASS everywhere"*, with the two untestable drugs
  named and their exclusion explained.
* §3 per-drug headers for both drugs now state the zero-width band explicitly and quote
  the broad band as context, instead of presenting it as the bar.
* §4 flip map: the two rows previously read **ROBUST PASS**. A drug with no bar cannot be
  robust; they now read *"NOT TESTABLE at every D (no verdict, robust or otherwise)"*.
* §8 PASS counts: `6/6` and `5/6` -> `4/4` and `3/4`, with the untestable rows counted in
  no denominator.

**Does any conclusion rest on a moved verdict? No.** The document's three findings — that
`recal` absorbs a several-fold `D_tissue` change into the re-fitted sink, that the scored
ratio therefore barely moves, and that the measurement-backed diffusivity does **not**
rescue docetaxel's FAIL — are all about drugs whose verdicts did not move. The
recommendation is made on measurement provenance and was never a function of PASS count.

### 3.2 `results/grid_convergence_breast.md` — 14 cells (+4 extrapolated); conclusion unaffected

Both drugs carried a verdict at all 7 rungs, plus their Richardson and power-law
extrapolated verdicts.

* §0: *"5/6 verdicts are constant at every grid"* -> *"Of the 4 of 6 swept drugs that carry
  a scorable band, 3/4 verdicts are constant at every grid"*.
* §1 control table now matches the corrected published table (`NOT TESTABLE` on both
  sides), and still reports **verdicts ALL MATCH** — the control still passes.
* §3 headers for both drugs replaced with an explicit "NO TESTABLE BAND" statement.
* §4 summary: `published verdict` / `verdict @ dx=8` / `verdict at the converged value`
  all now `NOT TESTABLE` for the two drugs.

**Does any conclusion rest on a moved verdict? No.** The document exists to answer whether
docetaxel's FAIL is a grid artefact; it concludes FAIL at dx = 12, at dx = 8, under
Richardson extrapolation and under the whole-ladder fit. Docetaxel is n = 2 and untouched.
§5's boundary-layer correlation (`r = -0.930`, n = 6) is computed from decay lengths and
ratio spreads, not verdicts, and is unchanged.

**Incidental finding — a pre-existing generator/document drift, now closed.** Regenerating
this file also added §1b (the CFL question), §1c (no grid-error cancellation) and the
per-drug "successive refinement steps" lines, none of which the published document
carried: it had been generated by an earlier revision of `write_report`. This is the same
class of defect `2466ce1` was about, found in the same file for an unrelated reason.
One honest gap remains and is stated in the regenerated text rather than papered over:
§1b's `dt`-control points are **absent**, because `results/sweep/gridconv_dtcontrol.jsonl`
is not on disk and producing it needs solves.

### 3.3 `results/delivery_source_geometry.md` — 6 cells; ONE STATED CONCLUSION MOVED

Both drugs carried a verdict in all three source arms (`binary`, `subgrid`,
`subgrid_matched`).

* §5.2 verdict table: both rows' band column becomes *"none — n = 1, zero-width"*, their
  three per-arm margin percentages are removed (a margin to a band you are not held to is
  not a margin), and their verdict column becomes **NOT TESTABLE (no verdict in any
  arm)**. Their three ratios per drug are unchanged and still shown.
* The §5.2 callout heading: *"Three of the nine flip"* -> *"Three of the seven scored
  verdicts flip"*. The three flips themselves — GBM temozolomide PASS->FAIL, GBM cisplatin
  PASS->FAIL, breast docetaxel FAIL->PASS — are all n >= 2 drugs and are untouched.

**The conclusion that moved.** The final bullet read:

> *Only four of the nine (breast doxorubicin, paclitaxel, epirubicin, 5-fluorouracil) keep
> a margin above 80 %.*

Two of those four survivors were surviving against the broad context band. Corrected:

> *Only two of the seven scored verdicts (breast doxorubicin, breast epirubicin) keep a
> margin above 80 %.*

**The set of verdicts that comfortably survive the source-geometry swap is half the size
this document reported.** That strengthens its recommendation (do not promote
`geometry='subgrid'` into the delivery source) rather than weakening it — but it moved
because a substituted verdict was load-bearing for it, which is exactly what this exercise
was looking for.

The document's headline is not affected: it rests on the three flips, on the alpha arm
being provably degenerate, and on the shape-not-dose result in §5.3, none of which
involves paclitaxel or 5-fluorouracil as a verdict. §5.3's exposure and ratio folds for
both drugs are ratio quantities and are unchanged.

### 3.4 `results/p99_convention.md` — 30 cells; the comparison table survives, and here is why

This is the document the brief singled out, because it compares normalising conventions by
how many published verdicts each reproduces. Both drugs carried a verdict under **all
fourteen** conventions the driver sweeps (twelve of which the document tabulates), plus a
row each in the §6 bias table.

* §5 verdict tables: both breast rows become `—` (NOT TESTABLE) in every convention
  column, band column becomes *"none (n = 1)"*, ratios unchanged.
* §5 PASS-count row, **the load-bearing table**:

| | p50 | p90 | p95 | **p99** | p99.9 | max | **tail90** | tail95 | tail99 | mean | rms | e-wtd mean |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **was** (out of 9) | 7 | 8 | 8 | 8 | 8 | 8 | 8 | 8 | 9 | 7 | 7 | 7 |
| **now** (out of 7) | 5 | 6 | 6 | 6 | 6 | 6 | 6 | 6 | 7 | 5 | 5 | 5 |

**Every column drops by exactly 2.** Paclitaxel and 5-fluorouracil scored PASS against the
broad band under every single convention, so they contributed a constant `+2` to every
total and discriminated between no two conventions. The consequence is that the
document's comparison is unchanged in shape:

* `tail90`, the convention selected on numerical grounds, still reproduces the published
  verdict set exactly — now stated as **6/7 PASS, same 6, same 1 FAIL** (was 8/9).
* `tail99` is still the unique candidate that would clear the docetaxel FAIL — now
  **7/7** (was 9/9) — and is still rejected, on grid error and bootstrap variance.
* The moment conventions (mean / rms / L4 / e-weighted mean) still cost GBM temozolomide,
  the flagship verdict. `p99.9` and `max` still swap GBM cisplatin for breast docetaxel.
* No ranking, no gap, and no recommendation changes.

* §6 bias table: both rows lose their margin columns (they were **+265.1 %** and
  **+720.9 %** to the context band, presented as verdict margins) and their
  `PASS -> PASS / changed: no` cells become `NOT TESTABLE (no verdict to change) / n/a`.
  Their corrected ratios (0.1405 and 0.3158) are retained — the correction arithmetic is
  about the ratio and is unaffected.
* §6's answer changes denominator only: *"one of nine"* -> *"one of the seven scored
  verdicts"*. Breast docetaxel, at a **+1.5 %** margin, is still the only verdict the
  grid bias puts at risk, and §6.1's refusal to read that as "docetaxel passes" is
  untouched.
* §0 bullets 6, 7 and 8, and the "defensible statement" in §6.1, were updated to the
  corrected denominators. Where a sentence was about the nine tissue **ratios** rather
  than verdicts (the -2.43 % to +3.79 % spread, the sign-agreement count) the nine was
  left standing and the word "verdicts" corrected to "ratios", because all nine ratios
  are real.

**Does any conclusion rest on a moved verdict? No — but only because of a fact that had
to be checked rather than assumed.** Had paclitaxel or 5-fluorouracil flipped under any
convention, the comparison table would have been ranking conventions partly on
broad-band substitutions and would have needed to be rebuilt. They did not; they were
constant. That is now stated in the document, next to the table, with the old row quoted
so a reader can see the shift is uniform.

---

### 3.5 A second defect landed under this work, and two of these documents are stale because of it

While this was in progress, commit `5d6e8e4` ("Land the `_TUMOR_SIGMA` units fix")
changed the interstitial-pressure solve. The gaussian that smooths the tumour-occupancy
field had its width expressed in **grid cells** (`_TUMOR_SIGMA = 5.0`) and handed
straight to `ndimage.gaussian_filter`, whose `sigma` is in pixels — so the physical
smoothing length tracked the mesh (40 um at `dx = 8`, 60 at 12, 80 at 16, 120 at 24).
It is now a physical length, `_TUMOR_SIGMA_UM = 60.0`.

**`dx = 12 um` is the exact fixed point** — `60.0 / 12.0 == 5.0` in IEEE-754, so
`gaussian_filter` receives a bit-identical argument and every shipped `dx = 12` result
is unchanged, verified field-by-field in `results/tumor_sigma_units.md` §4.

That matters here because `phi` sets both the vessel-leak gain and the lymphatic
suppression in `solve_ifp`, whose Darcy velocity feeds the **convection term of every
transport solve** in these four drivers (`kill_calibration.prepare_tissue` ->
`unit_transport`, `convection=True`).

| document | grid | affected by `5d6e8e4`? |
|---|---|---|
| `delivery_source_geometry.md` | `dx = 12` only | **No** — exact fixed point |
| `p99_convention.md` (stage D) | `dx = 12` panels | **No** — exact fixed point |
| `dtissue_sensitivity.md` | sweep at **`dx = 24`**; §7.3 confirmation at `dx = 12` | **YES** for the sweep |
| `grid_convergence_breast.md` | ladder **8/10/16/18/24/30** plus `dx = 12` | **YES** for every rung but 12 |

`results/tumor_sigma_units.md` §3 measures the size of it: mean tissue IFP at `dx = 24`
moves **0.181899 -> 0.132363**, and the spread across the ladder collapses from
**1.89x to 1.09x**.

**Both affected documents were regenerated with zero solves, so their RATIOS predate the
fix and were not recomputed.** The band decision is exact regardless — it consumes only
the stored `ratio_exact` scalar — so every verdict reported in them is correct *for the
ratio shown*. The ratios themselves are stale at every grid but 12. Rather than let a
fresh regeneration timestamp imply otherwise, **a STALE TRANSPORT WARNING is now emitted
by both generators** and stands at the head of both documents, naming the affected rungs,
the fixed point, the measured IFP movement, and the solve cost of fixing it (120 and 42
transport solves respectively).

**The sharpest consequence is in `grid_convergence_breast.md`, and it is stated there.**
That document is a *grid-convergence study*, and the superseded constant made the model
change with the mesh — so its ladder was measuring discretisation error and a model
change at once, which is exactly the confound a convergence study cannot tolerate. Its
Richardson and power-law extrapolations to `dx -> 0` inherit that. **Its docetaxel
conclusion should be re-derived from a re-solved ladder before it is quoted again.** That
is a pre-existing defect that landed under this work rather than one this work
introduced, and re-solving it is out of scope for a no-solve band-rule propagation — but
it would have been dishonest to regenerate the file today and leave the point unmade.

---

## 4. What this does NOT fix, and what is now guaranteed

**Guaranteed.** There is exactly one place in the repository that decides whether a
3-D:2-D row is scored and against what: `kill_calibration.score_ratio_row`. It is called
by the scorecard, the breast report generator, and all four of these drivers. No consumer
re-derives it. `ratio_test` independently refuses to return `FAIL` on a degenerate band,
so a caller that forgot the guard still cannot manufacture a failure.

**Not fixed, and out of scope here.** Two of these four documents are hand-written prose
over generated JSON. Their numbers are now correct and were transcribed from the
regenerated artefacts, but nothing *enforces* that they stay in sync the way a
`write_report` does — the same structural weakness that let this defect persist in
`grid_convergence_breast.md`'s §1b/§1c for an unrelated reason (§3.2). That is a real
remaining exposure and is recorded here rather than claimed as closed.

**Not fixed: the underlying literature gap.** Paclitaxel and 5-fluorouracil are not
untestable because the model is deficient. They are untestable because exactly one study
each exists, and a single study is a point estimate. If a second independent 3-D-vs-2-D
measurement is published for either drug, both become testable immediately with no code
change — `score_ratio_row` reads `n` from the bridge. Until then, the honest output is a
ratio, a point, and no verdict.

---

## 5. Reproducing this

```
./.venv/bin/python -m experiments.dtissue_sensitivity     --report-only --dx 24   # emits the §3.5 staleness banner
./.venv/bin/python -m experiments.grid_convergence_breast --report-only              # emits the §3.5 staleness banner
./.venv/bin/python -m experiments.delivery_source_geometry --stage blast --dx 12
./.venv/bin/python -m experiments.p99_convention          --stage D
./.venv/bin/python -m pytest -q          # 775 passed
```

None of these solves a PDE. The third and fourth regenerate JSON; the corresponding
prose in `results/delivery_source_geometry.md` §5.0/§5.2 and `results/p99_convention.md`
§0/§5/§6 was updated from that JSON by hand.

See also: `results/zero_width_bands.md` (why a single study is a point estimate and not an
interval), `results/kill_calibration_breast.md` §3 (the corrected published panel),
`results/tail90_falsification.md` §1 (what is and is not persisted on disk).
