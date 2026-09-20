# POST-HOC REPORT OF AN EXISTING COMPUTATION — seven GBM sections, one compound, solved 2026-08-23 and read for the first time today

> ## CORRECTION (2026-09-20, same day) — two statements in this file are contradicted by its own tables
>
> Added by a later pass in the same session. **The tables are correct; two sentences that summarise them are not.** Appended rather than edited, which is this directory's rule: a dated record is corrected in view of the reader, never rewritten underneath them.
>
> **1 — "`F(0.1)` = 0 everywhere" is contradicted by this file's own §5.5 and §5.6.** §4.4 item 4 reads *"`F(0.1)` = 0 everywhere, single-cell and spot alike. Stable at the value that makes it useless."* Two tables below it print otherwise:
>
> | section | `F(0.1)` as printed here |
> |---|---:|
> | `visium_10x_parent_gbm` (§5.5, §5.6) | **0.001811** |
> | `xenium_319 @ dx 12` (§5.6) | **0.000498** |
>
> Six of eight entries are 0; two are not. The corrected statement is **"`F(0.1)` = 0 on six of the eight entries, and below 0.2% on the other two"** — which supports the same conclusion (the readout carries almost no information at this threshold) without asserting a uniformity the data does not show. The word doing the damage is *everywhere*: it converts a strong empirical pattern into a false universal, and a reader checking one row finds the file disagreeing with itself.
>
> **2 — "ranks first" overstates where `CLAIM_INVENTORY.md` §6 puts a second section.** §0 of this file says §6 ranks *"A second section, same tumour type"* **first** among everything that would change the project's honest self-assessment. §6's list is numbered **0–5**, and that item is number **1**. Item **0** is *"Reconciling `vivum validate` with `kill_calibration_breast.md`"*, which states of itself: *"It is also the cheapest item on this list."*
>
> So §6 ranks a second section **second**, behind an unresolved disagreement between the shipping tool and the flagship report about how many drugs fail. That does not weaken this file's reason to exist — the item is still ranked above everything else — but "first" is not what the source says, and the item actually ranked first is a live internal contradiction about a *count that is currently quoted*.
>
> Provenance, because it matters for how far this spreads: the "first" phrasing is inherited from `verification_routes_2026-09-20.md` R12 rather than invented here. That is the same claim-surface drift this project keeps finding — a phrasing propagating between documents faster than anyone re-checks it against the source.
>
> Nothing else in this file is disputed by this correction. Its central result — that the per-cell verdict takes three different values across sections under byte-identical thresholds, and that the flagship is the outlier — is unaffected by either item.


**Status: POST-HOC REPORT OF AN EXISTING COMPUTATION. NOT a pre-registered
test.** The computation this document describes was run on **2026-08-23** and
written to `results/multisection_generalisation.json`. I read that file on
**2026-09-20**, four weeks later. Nothing in this document was declared before
the numbers existed, and nothing in it may be read as though it had been.

That distinction is the whole reason for this header. **You cannot pre-register
a test against data that already exists.** This project withdrew a headline
claim on 2026-08-21 over exactly this species of error — a pre-declared anchor
failed and an unrelated floor was substituted after the fact
(`results/codon_transport_capacity_gap` thread; `results/anchor_registration.md`,
`results/capacity_ranking_robustness.md`). So no threshold is invented here and
then presented as fixed beforehand. §6 states what the natural pre-registration
*would* have been and whether the result *would* have passed it, and labels that
judgement as the retrospective thing it is.

**Calibration guard, start of pass:**

```
$ .venv/bin/python -c "from tmesim.validation.calibration import calibrate_alpha; print(repr(calibrate_alpha()))"
START calibrate_alpha() = 22.223477711172563
ASSERT OK
```

*(The closing assertion is in §9.)*

**Nothing under `tmesim/`, `app/`, `docs/`, `landing-3d/` or `experiments/` was
changed by this pass, and no existing file in `results/` or `archive/` was
edited. This file is the only thing created.** The re-runs in §3 were executed
against the working tree as other agents left it and wrote nothing into the
repository.

---

## 0. Why this file exists at all

`results/CLAIM_INVENTORY.md` §6 ranks, **first** among everything that would
change this project's honest self-assessment:

> **A second section, same tumour type.** Almost every defensible claim above is
> *n* = 1 in tissue. This is the cheapest, largest improvement available.

Seven sections were solved on 2026-08-23. The output has been on disk ever
since. As of this morning, **no file in `results/`, `docs/`, `README.md`,
`app/` or `landing-3d/` mentioned `multisection_generalisation`** — the only
hits in the entire tree were the two `experiments/` scripts that produce it and
`results/stale_blockers_2026-09-20.md`, the audit that found it:

```
$ grep -rl "multisection_generalisation" results/ docs/ README.md app/ landing-3d/ tmesim/ experiments/ tests/
results/stale_blockers_2026-09-20.md
experiments/multisection_generalisation.py
experiments/multisection_report.py
```

*(This file now makes a fourth hit. It is the first one in `results/` that is a
reading of the data rather than an audit noting that nobody had read it.)*

`results/verification_routes_2026-09-20.md` R12 — written **today**, by a
thread that did not know the runs existed — budgets **4–6 engineer-days** to
run what had already been run, and describes the project as "currently costing
four sections of replication for free."

### 0.1 One place it WAS written down, and the form it took there

The findings are not wholly unrecorded. The commit that added the JSON —
`bf15977`, 2026-08-24, *"The headline is capacity-sensitive, and the flagship
section is not representative"* — carries them in its **message**:

> One compound, one configuration (dx = 30, 6 h), seven sections: …
> The dominant limiter takes THREE different values across seven sections.
> Therapeutic coverage spans the entire range 0.0 to 1.0. p90/p10 spans 6x.

A commit message is a record, not a claim surface. Nothing reads it, no guard
checks it, and four weeks later `CLAIM_INVENTORY.md` §6 still ranks "a second
section" first without citing it.

**And note the form those three summary sentences take.** All three are computed
across **seven** sections — a set that includes the four whose per-cell spatial
statistics the shipped code refuses (§1.4). The same message states the refusal
correctly two paragraphs further down (*"The four Visium rows are SPOT data,
where per-cell spatial statistics are refused in code … They are recorded, not
tabulated as tissue readings, and the refusal is the result"*), so the prose was
right and the summary line was not. That seven-section **6×** range is exactly
the counterfactual quarantined in §5.4 of this document. The corresponding
single-cell figures are **1.2687×** (§4.2), and a dominant limiter taking **two**
values, not three (§4.1).

This is the claim-surface drift defect the repository already tracks
(`results/claim_surface_audit_2026-08-30.md`), arriving on the **understatement**
side, where no guard is watching.

---

## 1. What was actually run

### 1.1 The file

```
results/multisection_generalisation.json
  129,325 bytes   mtime 2026-08-23 03:43:53   sha256 d6b59ce56d2abae5…
  top-level keys: alpha, drug, runs, updated
  alpha   = 22.223477711172563          (the live value; asserted above)
  drug    = "temozolomide"              ONE compound
  updated = "2026-08-23 03:43:53"
  runs    = 8 entries, every one status = "ok"
  every entry carries alpha_at_run = 22.223477711172563
```

**Eight entries, seven sections.** The eighth entry is `xenium_319` re-solved at
a finer grid. It is a paired control, not an eighth specimen, and it is kept out
of every cross-section statement in this document.

### 1.2 The RunSpec — identical on all seven

Read out of the JSON, not retyped from the driver:

```json
{"section": "<varies>", "dx_um": 30.0, "drug": "temozolomide",
 "sensitivity_mode": "uniform", "vascular_mode": "discrete",
 "oxygen_mode": "pde", "state_mode": "oxygen", "mgmt": null,
 "pk": {"plasma_peak": 50.0, "half_life_s": 14400.0,
        "n_doses": 1, "dose_interval_s": 86400.0},
 "total_time_s": 21600.0, "attribute_mechanism": true, "timelapse_frames": 0}
```

A field-by-field diff of all eight specs against `xenium_319@dx30` returns
**exactly one differing key per entry**: `section` for six of them, `dx_um`
(30.0 → 12.0) for the eighth. Nothing else varies. `transvascular_exchange` is
`False`, the shipped default, on every run.

The driver is `experiments/multisection_generalisation.py`, and it calls
`app.server.jobs.execute_run` — the same function the Studio queue calls — so
these are the product's shipped readouts rather than a re-implementation of
them. All seven sections are **built-in registry sections that resolve by
filesystem path** under `data_cache/`.

### 1.3 What `status: "ok"` means, exactly

It means `execute_run` returned without raising. Nothing more. From the driver:
the `except Exception` branch writes `status: "FAILED"` with a traceback; the
success branch hard-codes `"status": "ok"` in `_extract`. **It is not a quality
verdict and not a statement that any number in the entry is meaningful.**

The pipeline's actual refusal channels are separate fields, and they are
recorded verbatim:

| refusal channel | fired on |
|---|---|
| `metric_reasons` (per-metric refusal) | **0 of 7** sections |
| `penetration_fit_error` | **0 of 7** — the exp+floor fit converged everywhere |
| `reach.null_reason` | **0 of 7** |
| `section_caveats` → `spatial_spread: NO` | **4 of 7** (every Visium section) |

Seven sections, verdicts that contradict each other (§4), and every per-metric
refusal channel silent. **Disagreement across tissue is not something this
pipeline can currently detect about itself.**

### 1.4 Which of the seven are refused, by name

`app/server/api.py::_SPOT_SOURCES = ("visium",)` and
`app/server/api.py::SPOT_SPATIAL_REFUSAL` attach a blocking `NO` caveat to every
spot-resolution section. Checked entry by entry:

| section | platform | `spatial_spread_valid` | blocking caveats |
|---|---|---|---|
| `xenium_319` | xenium | `True` | none |
| `xenium_io` | xenium | `True` | none |
| `gse324268` | xenium | `True` | none |
| `visium_10x_parent_gbm` | visium | **`False`** | `spatial_spread: NO` |
| `visium_gbm1` | visium | **`False`** | `spatial_spread: NO` |
| `visium_gbm2` | visium | **`False`** | `spatial_spread: NO` |
| `visium_gbm3` | visium | **`False`** | `spatial_spread: NO` |

The refusal is not a soft caveat. In the shipped constant's own words: on a
Visium section the pipeline "synthesises 12 pseudo-cells per spot at the spot
centroid with uniform jitter", and therefore **"THE DELIVERY MAP IS INVALID
HERE TOO, AND FOR THE SAME REASON … the falloff curve bins by distance to the
nearest endothelial cell …, the perivascular reference is the first bin of that
curve, and the delivery spread is a percentile ratio over the resulting field —
so all three measure the 12-pseudo-cells-per-spot synthesis, not this tissue's
vascular architecture. Read NOTHING on this section as a statement about where
its vessels are."**

**So the generalisation question can be asked of three sections, not seven.**
The four Visium sections produce numbers — the engine does not raise — and those
numbers appear in §5.4 **quarantined**, so the refusal is auditable. They are
not evidence about those tumours, and **no average anywhere in this document is
taken over them, or over them and the Xenium sections together.**

---

## 2. What is NOT comparable to what

The task this file answers asked for non-comparabilities to be reported rather
than averaged over. There are five, and the first three are first-order.

### 2.1 Four of seven sections are refused on the readouts that matter

§1.4. This is the largest single fact in the file.

### 2.2 `xenium_io` is typed by a DIFFERENT marker map — and the endothelial call is the source term of the entire delivery solve

`tmesim/reconstruct/marker_panels.py::select_marker_map` picks the map whose
markers the loaded panel actually measures. Run against the three cached panel
JSONs:

```
$ .venv/bin/python  (select_marker_map over each section's own gene_panel.json)
xenium_319   342 names  coverage {brain_319: 1.00, io_480: 0.50,  breast_313: 0.50 }  -> brain_319  ENDOTHELIAL = ['PECAM1']
xenium_io    503 names  coverage {brain_319: 0.50, io_480: 1.00,  breast_313: 0.714}  -> io_480     ENDOTHELIAL = ['VWF', 'FLT1', 'PLVAP']
gse324268    396 names  coverage {brain_319: 1.00, io_480: 0.786, breast_313: 0.643}  -> brain_319  ENDOTHELIAL = ['PECAM1']
```

`PECAM1` is **absent** from the IO panel (verified directly in
`data_cache/xenium_gbm_io/Xenium_GBM_gene_panel.json`), so `xenium_io`'s
endothelial class is called by a **three-gene** gate where the other two
sections use a **one-gene** gate. `vascular_mode = "discrete"` builds the
delivery source from exactly those endothelial cell positions, so **everything
downstream of the endothelial call — `mean_vascular_density`, `E_ref`, the
falloff, the delivery spread, the fitted penetration length — is produced by a
different measurement rule on `xenium_io` than on the other two.**

The module documents this as a known hazard in its own docstring: the
single-vs-multi-marker imbalance once "reversed its sign on this section", with
the endothelial fraction reading "4.18% here against 9.74/12.19% on the two
PECAM1 sections. That is the arity bug, and it is fixed in the SCORING
(`real_xenium.TYPING_AGGREGATION`), not in this map." The fix is claimed at the
scoring layer; the **gene sets remain different**, and no section in this
checkout carries both panels, so there is no measurement in this repository
showing the two maps agree on one tissue.

**Consequence, and it is a good one:** the two sections that disagree most
sharply — `xenium_319` (limiter `dose`) and `gse324268` (limiter `resistance`) —
are typed by the **same** map, `brain_319`, with the **same** one-gene
endothelial gate. **The verdict flip in §4 is therefore not a marker-panel
artefact.** The caveat lands on `xenium_io` specifically, and on any three-way
spread that includes it.

### 2.3 The eighth entry is a different grid, on one section only

`xenium_319@dx12` is the only entry at dx = 12 µm. It is a paired grid control
for the flagship and for nothing else. Neither section that carries the verdict
flip was ever solved at a second grid. See §7.1 — the existing renderer gets
this wrong in prose.

### 2.4 The tissue mask, and a 102× range in voxel count

| | Xenium sections | Visium sections |
|---|---|---|
| `tissue_mask_mode` | `geometric` | full bounding rectangle |
| `tissue_coverage_fraction` | 0.672 / 0.715 / 0.647 | 1.000 on all four |
| voxels at dx = 30 µm | 46,070 – 169,087 | 1,665 – 3,190 |

Every mask-area average — `mean_tissue_auc`, each falloff bin mean, the tissue
mean of anything — is an average over a different number of voxels on each row,
and the Visium rows have no tissue mask at all. Across the seven the voxel count
spans a factor of **102**.

### 2.5 The seven do not carry the same caveat set

`xenium_319` carries **eight** run caveats; `xenium_io` and `gse324268` carry
**seven**. The extra one on the flagship is section-specific and it matters:

> **ABSOLUTE KILL MUST NOT BE QUOTED FROM THIS SECTION UNDER ANY STATE ROUTE.**
> Its panel carries none of `['VEGFA', 'CA9', 'ADM', 'NDRG1', 'SLC2A1']`, so
> hypoxia is UNMEASURED here …

So the section whose `hypoxic_tumor_fraction` reads **0.00%** is also the one
section of the three whose panel cannot check a hypoxia call against anything.
That is not a reason to discard it. It is a reason not to read its 0% as a
measurement of that tumour.

---

## 3. Does it still reproduce? — four entries re-solved today

The file is four weeks old and two changes have landed on the run path since.
A stale number is worse than no number, so four entries were re-solved against
the working tree exactly as other agents left it this morning, and every field
diffed against the stored entry.

### 3.1 The result

| entry re-solved today | fields compared | identical | differ |
|---|---:|---:|---:|
| `xenium_319@dx30` (the flagship) | 48 | **48** | 0 |
| `visium_gbm2@dx30` (a spot section) | 48 | **48** | 0 |
| `gse324268@dx30` (851,509 cells) | 48 | **48** | 0 |
| `xenium_io@dx30` (816,769 cells) | 48 | 44 | **4** |

"Identical" is `==`, not `allclose`. It includes every float to the last digit,
the full seven-rung reach ladder, the full falloff arrays, the complete
`mech_counts` integer histogram, the survival breakdown and all thirteen
`drug_params` fields. `calibrate_alpha()` was asserted at
`22.223477711172563` before and after the re-run.

**Three of the four reproduce bit-for-bit. The transport solve, the
reconstruction, the oxygen solve, the kill model and the per-cell attribution
all produce the same numbers today that they produced on 2026-08-23.** The
stored file is live, not stale.

### 3.2 The one entry that does NOT reproduce, and why that is the good news

On `xenium_io@dx30`, four fields differ. All four are in the reach block, and
all four moved in the same direction — from **a number** to **a refusal**:

| field | stored 2026-08-23 | current engine |
|---|---|---|
| `reach_headline_F` | `0.0` | `None` |
| `reach_ladder_F` | `[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]` | `[]` |
| `E_ref` | `28.656596966167733` | `None` |
| `falloff_normalised_mean_auc` | 35 bins | 0 bins |

with a reason now attached where there was none:

> only 1 VIABLE tumour cell in this reconstruction, below the 30 this repository
> requires before a mean over a sample is reported
> (`vessel_profile.MIN_BIN_COUNT`, the same floor that gates the perivascular
> reference this ratio is divided by). F(x) = counts / 1, so every rung could
> only land on a multiple of 1 — that is a COUNT, not a fraction, and it is
> refused rather than reported with a caveat

Everything on that section **outside** the reach block still reproduces exactly:
`E_p90_over_p10 = 2.3211182478545296`, `penetration_length_um =
133.16164748093735`, `mech_dominant_limiter = 'resistance'`, `n_viable = 1`,
`n_tumor = 540533`, the entire mechanism histogram.

### 3.3 The two changes on the run path, and their dates

```
$ git log --since=2026-08-23 -- app/server/jobs.py app/server/reach.py tmesim/sim/ tmesim/physics/ tmesim/reconstruct/ tmesim/core/
```

**(a) `cab3c8e` — 2026-08-23 03:45:23, ninety seconds after this JSON was
written.** `app/server/reach.py` now refuses the ladder when `n_viable <
MIN_BIN_COUNT` (= 30). Its commit message:

> Found by running the shipped pipeline off the flagship section for the first
> time. On `xenium_io` — same platform, same disease, same configuration —
> `compute_reach` returned **n_viable = 1** and a clean-looking ladder with
> F(0.1) = 0.0 and no caveat anywhere on the panel. … THE FLOOR IS BORROWED,
> NOT INVENTED. `vessel_profile.MIN_BIN_COUNT = 30` is this repository's own
> existing answer to "how many samples before a mean means anything" … Picking
> a fresh number here, after seeing that `xenium_io` gave 1, would have been a
> threshold chosen knowing what the data does.

**This sweep already caused one shipped-code fix, ninety seconds after it
finished — and that fix is the only trace it left anywhere.** The sweep itself
was never written up. The defect it exposed could not have been found on the
flagship, whose same denominator is 2,019 cells.

**(b) `962b53a` — 2026-08-27.** `app/server/jobs.py` now attaches
`verdict_stability` to every run payload. For the compound this entire sweep is
built on:

```
$ .venv/bin/python -c "from app.server.jobs import verdict_stability; print(verdict_stability('temozolomide'))"
{'compound': 'temozolomide', 'stable': False, 'measured_on': 'xenium_319',
 'source': 'results/verdict_inversion_census.md',
 'note': 'temozolomide REVERSES its delivery/dose verdict when the missing
          transvascular exchange term is switched on, so the verdict above is
          one of the unstable ones.'}
```

That annotation did not exist when the sweep ran. It lands directly on §4 and
is carried there.

**Neither change touches the transport law, the reconstruction or the kill
model — which is exactly why 44 of 48 fields on `xenium_io` and 48 of 48 on the
other three are bit-identical.** Both changes are refusals and disclosures. The
engine got more honest, not different.

### 3.4 What this means for the stored file

**Seven of the eight stored entries are reproducible as they stand.** One entry
— `xenium_io@dx30` — has a reach block that **the shipped product would refuse
to emit today**. That block must not be quoted from the JSON. It is superseded
in code, and this document treats it as refused everywhere below.

---

## 4. What generalises, and what does not

The interesting output is the disagreement, not a mean. **No mean over seven
sections appears anywhere in this document, and none over three either.** The
four Visium sections contribute nothing to this section (§1.4).

### 4.1 The verdict the product exists to answer FLIPS

`tmesim/sim/mechanism.py` labels every tumour cell with why it survived and
reconciles the labels against `kill_fraction` bit-for-bit or raises.

| section | marker map | dominant limiter | its share |
|---|---|---|---:|
| `xenium_319` | `brain_319` | **`dose`** | 0.9985 |
| `xenium_io` | `io_480` | **`resistance`** | 0.6090 |
| `gse324268` | `brain_319` | **`resistance`** | 0.6629 |

**This is not a shift in magnitude. It is a different answer.** And it is not
explained away by §2.2: `xenium_319` and `gse324268` are typed by the *same*
marker map with the *same* one-gene endothelial gate, and they disagree.

**Where the flip comes from.** `mechanism.py`'s own documentation gives the kill
mapping as `m = sensitivity × resistance`, with `SURVIVED_RESISTANT` assigned
when exposure was adequate but `m <= resistant_modulation`, and
`resistant_modulation` "derived from the PDConfig the run actually used — it is
`max(hypoxia_resistance, quiescent_resistance)`". This sweep ran
`sensitivity_mode = "uniform"`, so `sensitivity ≡ 1` on every cell and the
**only** route into the resistance bucket is the cell-state term. And the
cell-state term tracks the field that differs most across these sections:

| | `xenium_319` | `xenium_io` | `gse324268` |
|---|---:|---:|---:|
| `hypoxic_tumor_fraction` | **0.00%** | 34.18% | 18.88% |
| `mean_oxygen` | 0.3507 | 0.2204 | 0.2874 |
| `mean_vascular_density` | 0.02664 | 0.09997 | 0.1228 |

**So the flip is real in the shipped output, and its proximate driver is the
least-validated field in the model.** The run's own caveats say so, in the
payload, unprompted:

* `oxygen (pde)`: "WEAK BUT BEST-IN-CLASS … against measured hypoxia markers on
  816,769 cells the solve scores rho **-0.1523**".
* `cell_state (oxygen)`: "**NOT A MEASUREMENT** … the tumour VIABLE fraction this
  produces swings across two tissue-mask conventions and three grid spacings,
  both of which are bookkeeping choices rather than measurements — 0.60%–98.19%
  (**165x**) …".
* and on `xenium_319` only, "hypoxia is UNMEASURED here" (§2.5).

Two further discounts on this flip, both from the shipped engine rather than
from me:

1. **The flagship's own split is over-determined.** `mech_joint_limited` on
   `xenium_319` is **19,596 cells — 92.69% of its exposure-limited survivals** —
   labelled `dose`-limited by the documented precedence rule while *also*
   sitting at or below the resistance threshold. Raising the dose alone would
   not have killed them. The clean-looking 99.85% `dose` share is not a clean
   counterfactual even on the one section it was measured on.
2. **Since 2026-08-27 the engine flags this very verdict as unstable for this
   very compound** (§3.3b): temozolomide `REVERSES` its delivery/dose verdict
   when the transvascular exchange term the transport law lacks is switched on.
   Measured on `xenium_319` only. So the flip in the table above is a
   disagreement between two verdicts that are each, separately, known to invert
   under a term that ships OFF.

**Thresholds were not tuned per section.** `mech_thresholds` is byte-identical
on every row: `{'unreached_ratio': 0.1, 'adequate_ratio': 1.0,
'resistant_modulation': 0.5, 'stochastic_min_p_kill': 0.25}`.

### 4.2 The one continuous readout that is close across tissue

| `E_p90/E_p10` (delivery spread) | `xenium_319` | `xenium_io` | `gse324268` | max/min |
|---|---:|---:|---:|---:|
| | 1.9936 | 2.3211 | 2.5293 | **1.2687×** |

`xenium_io` sits strictly between the other two, so dropping it for the
marker-map reason in §2.2 leaves the spread **unchanged at 1.2687×** over the
two same-map sections. This readout is also the one
`app/server/reach.py`'s `CONVENTION_NOTE` had already promoted to headline in
place of `F(0.1)`, on the flagship alone. **This is the first evidence that the
promotion picked a statistic that behaves the same way on more than one
specimen.** It is three specimens, and three specimens is three specimens.

### 4.3 Everything else moves, and several change interpretation

| readout | `xenium_319` | `xenium_io` | `gse324268` | |
|---|---:|---:|---:|---|
| fitted penetration length `L` (µm) | **277.69** | 133.16 | 103.18 | **2.691×**; flagship is the LONGEST |
| `therapeutic_coverage_1x` | **0.00%** | 72.71% | 75.89% | "no cell reaches IC50" is a fact about one section |
| `E_ref` perivascular reference (µM) | 7.802 | *refused today* | 35.621 | |
| falloff usable-bin support (µm) | 145 | *refused today* | 85 | no shared domain |
| `F(0.1)` reach headline | 0.0 | *refused today* | 0.0 | stable, and informationless |
| `F(1.0)` top ladder rung | 0.42051 | *refused today* | 0.0035057 | **119.9×** between the two that survive |
| `n_viable` (the reach denominator) | 2,019 | **1** | 43,643 | |
| `hypoxic_tumor_fraction` | 0.00% | 34.18% | 18.88% | a ratio is undefined against 0 |

`L` on the flagship is the **longest** of the three: the specimen this product
ships with is its most drug-penetrable one. `F(1.0)` — the ladder's only rung
that is not identically zero on the single-cell sections — differs by a factor
of **119.9** between the two sections on which it is still computed at all.

**`hypoxic_tumor_fraction` and `therapeutic_coverage_1x` have no max/min here on
purpose.** One of the three values is exactly 0 in each case, and a ratio
against 0 is not a number. §7.2 records that the existing renderer prints
"1.81×" and "1.04" for these two rows by silently dropping the zero — which
drops the flagship, the section that disagrees.

### 4.4 What DOES generalise, stated plainly

1. **The pipeline completes.** 7 of 7 runnable built-in sections, no crash, no
   metric refusal, wall clocks 6.4 s – 57.1 s. That is a real engineering
   result and it is the only sense in which "it generalises" holds without
   qualification.
2. **The per-cell attribution is internally consistent everywhere.**
   `mech_reconciliation.bit_exact` is `True` on **8 of 8** entries and
   `UNATTRIBUTED` is **0** on all eight. The *consistency* of the attribution
   generalises perfectly. Its *answer* does not — §4.1.
3. **The delivery spread stays within 1.2687× across three specimens** — §4.2.
4. **`F(0.1)` = 0 everywhere**, single-cell and spot alike. Stable at the value
   that makes it useless; the earlier decision to de-headline it was right and
   if anything understated.

### 4.5 Be explicit about n — seven sections is not seven patients

**The arm that carries the result is three sections, and the repository's own
provenance record says they are three different donors.**

| section | accession, as the registry records it | deposit |
|---|---|---|
| `xenium_319` | `Xenium_V1_FFPE_Human_Brain_Glioblastoma_With_Addon` (Xenium v1.3.0, FFPE; tissue from Avaden Biosciences) | 10x Genomics public |
| `xenium_io` | `Xenium_V1_Human_Brain_GBM_FFPE` (Xenium v2.0.0, FFPE) | 10x Genomics public |
| `gse324268` | GEO `GSM9572843`, `GBM_tumor07`, primary resected human GBM | GEO, Hölzel lab (Bonn) |

`docs/real_data/singlecell_spatial.md` states the separation of the two 10x
sections twice and in its own words — *"They are complementary sections of
**different patients**, so they cannot be merged cell-for-cell"* (:278) and
*"One section, one patient, and a different patient from the 319-gene
section"* (:415). `gse324268` is a separate deposit from a different
institution and is a third donor on any reading.

**The honest qualification:** neither 10x dataset ships a donor identifier in
this checkout's cached provenance, so two of the three separations rest on the
repository's documentation rather than on a field in the data. No section is a
serial section or a replicate of another — this project has been careful about
that distinction before, and declined to download 10x's
`Targeted_Visium_Human_Glioblastoma_Pan_Cancer` on the explicit grounds that it
is *"a targeted ~1000-gene panel of the SAME tissue block"* as an existing
sample (`docs/real_data/visium_gbm.md`).

**So: n = 3 sections from 3 donors, one compound, one configuration.** Every one
is glioblastoma, and every one is an untreated resection.

**The refused Visium arm is four sections and an unresolved contradiction.**
They contribute nothing to §4 (§1.4), but the record should say what it is:

* `docs/real_data/visium_gbm.md` says the three GSE194329 samples are *"three
  different tissues"*, each *"an independent primary human GBM tumor"*.
* `tmesim/data/real/visium.py:85` and the shipped section registry
  (`app/server/cache.py:571`, `tmesim/data/section_catalog.py:189`) label `GBM2`
  *"Glioblastoma, IDH-wildtype (**recurrent**)"* and *"Visium GBM · GBM2
  (recurrent)"*.

One of those two is wrong, **no patient identifier is recorded for any of the
four**, and the same GEO series is documented as containing a same-patient pair
(`GBM5_1`/`GBM5_2`, *"peri-tumoral pair"*). So "four sample ids" does not
license "four patients" for that arm. Recorded, not resolved — resolving it
would change nothing about §4.

**What this licenses on the claim surface.** The *n* = 1-in-tissue caveat that
`CLAIM_INVENTORY.md` attaches to nearly every surviving claim comes off only the
readouts in §4.4, and only in this form: *measured on three GBM sections from
three donors, one compound, one configuration, at one grid spacing*. Not
"replicated". Not with a confidence interval. **Three is not a sample, no
statistical claim is made anywhere in this document, and none can be.**

---

## 5. The tables

**Provenance of every table below.** They are the output of
`experiments/multisection_report.py`, which already existed — written
2026-08-23 alongside the sweep, and never run into `results/`. It interpolates
every figure from `results/multisection_generalisation.json`; nothing in it is
typed by hand. Regenerate with:

```
$ .venv/bin/python -c "import experiments.multisection_report as m; print(m.build(m.load()))"
```

That command was run today and its tables are quoted verbatim below. **Its
PROSE is not quoted, and §7 says why: three of its sentences are wrong, and its
default output path `results/multisection_generalisation.md` was deliberately
NOT written by this pass.**

### 5.1 What ran, and what it cost

| section | label | platform | cells | grid @dx | **voxels** | extent (µm) | tissue-mask mode | tissue coverage | vascular source | wall clock |
|---|---|---|---:|---|---:|---|---|---:|---|---:|
| `xenium_319` | Xenium GBM · 319-panel | xenium | 40,887 | 170×271 | 46,070 | 8130 × 5100 | geometric | 67.2% | DISCRETE segmented vessels from RE | 9.4 s |
| `xenium_io` | Xenium GBM · IO 480-panel | xenium | 816,769 | 363×383 | 139,029 | 11490 × 10890 | geometric | 71.5% | DISCRETE segmented vessels from RE | 48.3 s |
| `gse324268` | Xenium GBM · GSE324268 tumor 07 | xenium | 851,509 | 479×353 | 169,087 | 10590 × 14370 | geometric | 64.7% | DISCRETE segmented vessels from RE | 57.1 s |
| `visium_10x_parent_gbm` | Visium GBM · 10x parent section | visium | 41,560 | 57×50 | 2,850 | 5000 × 5700 | full rectangle | 100.0% | KDE of REAL ENDOTHELIAL spot fract | 7.6 s |
| `visium_gbm1` | Visium GBM · GBM1 (primary) | visium | 57,050 | 55×58 | 3,190 | 5800 × 5500 | full rectangle | 100.0% | KDE of REAL ENDOTHELIAL spot fract | 9.4 s |
| `visium_gbm2` | Visium GBM · GBM2 (recurrent) | visium | 34,707 | 54×55 | 2,970 | 5500 × 5400 | full rectangle | 100.0% | KDE of REAL ENDOTHELIAL spot fract | 8.3 s |
| `visium_gbm3` | Visium GBM · GBM3 (primary) | visium | 21,279 | 37×45 | 1,665 | 4500 × 3700 | full rectangle | 100.0% | KDE of REAL ENDOTHELIAL spot fract | 6.4 s |

Wall clocks were measured against three other agents running heavy simulations
concurrently on the same machine.

### 5.2 The per-section readouts — single-cell sections ONLY

The four Visium sections are absent from this table on purpose (§1.4).

| readout | **`xenium_319`** | **`xenium_io`** | **`gse324268`** | max/min |
|---|---|---|---|---|
| cells (total) | 40,887 | 816,769 | 851,509 | — |
| `n_tumor` | 23,812 | 540,533 | 515,887 | **22.7×** |
| `n_viable` (reach denominator) | 2,019 | 1 | 43,643 | **4.36e+04×** |
| `n_viable` / `n_tumor` | 8.4789% | 0.0002% | 8.4598% |  |
| cell fraction TUMOR | 58.2% | 66.2% | 60.6% | 1.14× |
| cell fraction ENDOTHELIAL | 10.5% | 7.1% | 8.8% | 1.49× |
| cell fraction IMMUNE | 19.5% | 11.4% | 12.6% | 1.71× |
| cell fraction STROMAL | 11.8% | 15.4% | 18.0% | 1.53× |
| mean vascular density | 0.02664 | 0.09997 | 0.1228 | **4.61×** |
| mean oxygen | 0.3507 | 0.2204 | 0.2874 | **1.59×** |
| hypoxic tumour fraction | 0.00% | 34.18% | 18.88% | **1.81×** |
| — *(order-of-magnitude only: the oxygen consumption constant is un-calibrated)* — |  |  |  |  |
| **`E_p90/E_p10`  (DELIVERY SPREAD — the headline)** | **1.994** | **2.321** | **2.529** | **1.27×** |
| `E_p10` (µM) | 135.1 | 453.2 | 522.1 | **3.86×** |
| `E_p50` (µM) | 185.2 | 714.6 | 802.2 | **4.33×** |
| `E_p90` (µM) | 269.4 | 1052 | 1321 | **4.9×** |
| **`F(0.1)` (reach headline)** | **0** | **0** | **0** | — |
| **`F(1.0)` (top ladder rung)** | **0.4205** | **0** | **0.003506** | — |
| `E_ref` (perivascular reference) | 7.802 uM | 28.66 uM | 35.62 uM | 4.57× |
| **fitted penetration length `L` (µm)** | **277.7** | **133.2** | **103.2** | **2.69×** |
| analytic penetration length (µm) | 284.2 | 247.3 | 249.5 | **1.15×** |
| falloff support (last usable bin, µm) | 145 | 345 | 85 | 4.06× |
| `therapeutic_coverage_1x` | 0.00% | 72.71% | 75.89% | — |
| median E / IC50 | 0.3055 | 1.179 | 1.324 | **4.33×** |
| mean tumour AUC (µM·s) | 4.225e+6 | 1.59e+7 | 1.859e+7 | **4.4×** |
| **dominant limiter** | **dose** | **resistance** | **resistance** | — |
| *kill_fraction (NOT a defensible output)* | *0.1121* | *0.3159* | *0.3481* | *—* |

**Two rows in that table carry a `max/min` that silently excludes a zero.**
`hypoxic tumour fraction` prints "1.81×" over values `0.00% / 34.18% / 18.88%`,
and the renderer's §7 prints "1.04" for `therapeutic_coverage_1x` over
`0.00% / 72.71% / 75.89%`. Both ratios are taken over the two non-zero entries
and drop the flagship — the section that disagrees. See §7.2. The correct
statement in both cases is that the ratio is undefined.

**The last row is the absolute kill fraction.** It is printed only because the
shipped payload carries it and silently dropping a number the product prints
would be its own kind of dishonesty. **Absolute kill fraction is never a
defensible output, no comparison anywhere in this document is built on it, and
the payload itself flags it** — `kill_fraction_not_defensible` is a field in
every entry.

### 5.3 The per-cell mechanism split — all seven, spot rows marked

`⚠︎spot` marks a section whose per-cell statistics are refused (§1.4). Those
four rows are shown so the refusal is auditable. **They are not evidence about
those tumours and nothing in §4 is computed from them.**

| section | KILLED | UNREACHED | UNDERDOSED | RESISTANT | STOCHASTIC | UNATTRIBUTED | **delivery** | **dose** | **resistance** | **chance** | **dominant** | reconciles bit-exact |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| `xenium_319` | 2,670 | 31 | 21,111 | 0 | 0 | 0 | 0.001466 | 0.9985 | 0 | 0 | **dose** | yes |
| `xenium_io` | 170,764 | 672 | 102,825 | 225,177 | 41,095 | 0 | 0.001817 | 0.2781 | 0.609 | 0.1111 | **resistance** | yes |
| `gse324268` | 179,598 | 33 | 86,613 | 222,918 | 26,725 | 0 | 9.813e-5 | 0.2576 | 0.6629 | 0.07947 | **resistance** | yes |
| `visium_10x_parent_gbm` ⚠︎spot | 19,452 | 0 | 6,835 | 0 | 8,147 | 0 | 0 | 0.4562 | 0 | 0.5438 | **chance** | yes |
| `visium_gbm1` ⚠︎spot | 19,939 | 0 | 16,369 | 0 | 7,771 | 0 | 0 | 0.6781 | 0 | 0.3219 | **dose** | yes |
| `visium_gbm2` ⚠︎spot | 15,678 | 0 | 0 | 0 | 2,364 | 0 | 0 | 0 | 0 | 1 | **chance** | yes |
| `visium_gbm3` ⚠︎spot | 10,855 | 0 | 0 | 0 | 1,780 | 0 | 0 | 0 | 0 | 1 | **chance** | yes |

`bit_exact` is `True` and `UNATTRIBUTED` is `0` on all eight entries.

### 5.4 The quarantine — what the four refused sections produced

Reproduced so the refusal in §1.4 is auditable rather than an assertion.
**None of these is a reading of those tumours.** On a spot section the
endothelial positions that source the delivery field are synthesised — 12
pseudo-cells per spot at the spot centroid with uniform jitter — so every
column below is a measurement of that synthesis.

| section | cells | `n_tumor` | `n_viable` | `E_p90/p10` | `F(0.1)` | `F(1.0)` | `L` (µm) | falloff support | dominant limiter | `coverage_1x` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|
| `visium_10x_parent_gbm` | 41,560 | 34,434 | 32,580 | 5.55 | 0.001811 | 0.8368 | 66.27 | 235 µm | chance | 69.52% |
| `visium_gbm1` | 57,050 | 44,079 | 31,348 | 9.845 | 0 | 0.8098 | 71.87 | 195 µm | dose | 51.53% |
| `visium_gbm2` | 34,707 | 18,042 | 18,042 | 1.62 | 0 | 0.5042 | 412.9 | 45 µm | chance | 100.00% |
| `visium_gbm3` | 21,279 | 12,635 | 12,635 | 1.818 | 0 | 0.5366 | 301.2 | 65 µm | chance | 100.00% |

For the record of how bad it would have been: had the refusal been ignored, the
delivery spread across all seven would have been quoted as **1.62–9.845
(6.08×)** and the penetration length as **66.27–412.9 µm (6.23×)**. Both ranges
are dominated by sections whose vascular geometry is synthetic. **That
seven-section range appears here once, as the counterfactual it is, and is
quoted nowhere else.**

### 5.5 The reach ladder `F(x)`

**The `xenium_io` row of this table is SUPERSEDED.** It is reproduced as the
file stores it; the shipped engine has refused that row since 2026-08-23
03:45:23 (§3.2), and it must not be quoted.

| section | `F(0.01)` | `F(0.02)` | `F(0.05)` | `F(0.1)` | `F(0.2)` | `F(0.5)` | `F(1)` | `n_viable` | `E_ref` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `xenium_319` | 0 | 0 | 0 | 0 | 0 | 0 | 0.4205 | 2,019 | 7.802 uM |
| `xenium_io` | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 **← ONE CELL** | 28.66 uM |
| `gse324268` | 0 | 0 | 0 | 0 | 0 | 0 | 0.003506 | 43,643 | 35.62 uM |
| `visium_10x_parent_gbm` ⚠︎spot | 0 | 0 | 0 | 0.001811 | 0.03511 | 0.4774 | 0.8368 | 32,580 | 56.8 uM |
| `visium_gbm1` ⚠︎spot | 0 | 0 | 0 | 0 | 0.0007656 | 0.3149 | 0.8098 | 31,348 | 46.57 uM |
| `visium_gbm2` ⚠︎spot | 0 | 0 | 0 | 0 | 0 | 0.006817 | 0.5042 | 18,042 | 109.5 uM |
| `visium_gbm3` ⚠︎spot | 0 | 0 | 0 | 0 | 0 | 0.0152 | 0.5366 | 12,635 | 103.1 uM |

### 5.6 The normalised falloff — the curves do not share a domain

A bin needs ≥ 30 grid points to be usable, so how far each curve extends is
itself a per-section measurement. **The `xenium_io` row is SUPERSEDED for the
same reason as §5.5.**

| section | 5 µm | 15 µm | 25 µm | 35 µm | 45 µm | 65 µm | 85 µm | 105 µm | 145 µm | last usable bin | n bins |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `xenium_319` | 1 | 0.951 | 0.906 | 0.893 | 0.869 | 0.813 | 0.742 | 0.71 | 0.522 | 145 µm | 15 |
| `xenium_io` | 1 | 0.964 | 0.913 | 0.844 | 0.776 | 0.637 | 0.529 | 0.452 | 0.228 | 345 µm | 35 |
| `gse324268` | 1 | 0.916 | 0.836 | 0.753 | 0.677 | 0.58 | 0.504 | *no usable bin* | *no usable bin* | 85 µm | 9 |
| `visium_10x_parent_gbm` ⚠︎spot | 1 | 0.985 | 0.894 | 0.761 | 0.659 | 0.563 | 0.524 | 0.509 | 0.397 | 235 µm | 22 |
| `visium_gbm1` ⚠︎spot | 1 | 0.99 | 0.981 | 0.856 | 0.723 | 0.537 | 0.452 | 0.36 | 0.319 | 195 µm | 20 |
| `visium_gbm2` ⚠︎spot | 1 | 1.01 | 0.981 | 0.958 | 0.849 | *no usable bin* | *no usable bin* | *no usable bin* | *no usable bin* | 45 µm | 5 |
| `visium_gbm3` ⚠︎spot | 1 | 0.993 | 0.999 | 0.959 | 0.86 | 0.84 | *no usable bin* | *no usable bin* | *no usable bin* | 65 µm | 7 |

### 5.7 The single paired grid control

Only `xenium_319` was solved at a second grid. Every other row's `@12` column is
empty because that run was never made — not because it agreed.

| section | `E_p90/p10` @30 | @12 | `L` (µm) @30 | @12 | limiter @30 | @12 | `F(1.0)` @30 | @12 | voxels @30 | @12 |
|---|---:|---:|---:|---:|---|---|---:|---:|---:|---:|
| `xenium_319` | 1.994 | 1.994 | 277.7 | 278.5 | dose | dose | 0.4205 | 0.373 | 46,070 | 286,200 |
| `xenium_io` | 2.321 | — | 133.2 | — | resistance | — | 0 | — | 139,029 | — |
| `gse324268` | 2.529 | — | 103.2 | — | resistance | — | 0.003506 | — | 169,087 | — |
| `visium_10x_parent_gbm` ⚠︎spot | 5.55 | — | 66.27 | — | chance | — | 0.8368 | — | 2,850 | — |
| `visium_gbm1` ⚠︎spot | 9.845 | — | 71.87 | — | dose | — | 0.8098 | — | 3,190 | — |
| `visium_gbm2` ⚠︎spot | 1.62 | — | 412.9 | — | chance | — | 0.5042 | — | 2,970 | — |
| `visium_gbm3` ⚠︎spot | 1.818 | — | 301.2 | — | chance | — | 0.5366 | — | 1,665 | — |

On the one section that has both, dx = 30 → 12 µm (a 6.2× increase in voxel
count, 46,070 → 286,200) moves the delivery spread by 1.0004×, `L` by 1.003×,
and `F(1.0)` from 0.42051 to 0.37301; the dominant limiter stays `dose`.
**Neither of the two sections that carry the verdict flip was solved at a second
grid, so whether that flip survives a change of mesh is UNTESTED.**

---

## 6. The post-hoc framing, stated once and properly

**The rule.** A pre-registration is a commitment made before the outcome can be
known. The outcome here has been on disk since 2026-08-23. **Nothing written
today can be a pre-registration of it**, and any threshold named now would be
named by someone who has already read the answers. This project withdrew a
headline claim on 2026-08-21 in precisely this shape — a pre-declared anchor
failed and a different floor was put in its place afterwards — which is why the
rule is stated here rather than assumed.

What is legitimate: to say what the natural pre-registration *would* have been,
and whether the result *would* have passed it, **labelled as the retrospective
judgement it is.** That is all §6.1 and §6.2 are.

### 6.1 The criteria that need no threshold — and therefore could not be tuned by having seen the data

These four are categorical. Their verdicts do not depend on any number chosen
today, which is exactly why they are the ones worth stating.

| | Criterion, as it would have been written on 2026-08-22 | Retrospective verdict |
|---|---|---|
| **P1** | The dominant limiter — the product's headline verdict — is the same on every single-cell section | **WOULD HAVE FAILED.** `dose` / `resistance` / `resistance` (§4.1) |
| **P2** | Every section either produces a readout or refuses it; none emits a number the pipeline cannot stand behind | **WOULD HAVE FAILED.** `xenium_io` emitted a full seven-rung ladder, an `E_ref` and a 35-bin falloff over a **one-cell** denominator, with `null_reason = None` (§3.2) |
| **P3** | The per-cell attribution reconciles bit-exactly with `kill_fraction` on every section | **WOULD HAVE PASSED.** 8 of 8, `UNATTRIBUTED = 0` everywhere |
| **P4** | The pipeline completes on every runnable built-in section without crashing | **WOULD HAVE PASSED.** 7 of 7 |

**This is a retrospective judgement, written on 2026-09-20 by someone who had
already read the file.** Its only claim to weight is that P1–P4 are
threshold-free, so no choice made today could have moved any of the four
verdicts.

### 6.2 The continuous readouts have NO pre-existing tolerance, and this document does not invent one

For `L` (**2.691×** across three sections), `E_ref`, the falloff supports
(145 / 85 µm on the two sections where they survive), `F(1.0)` (**119.9×**) and
the delivery spread (**1.2687×**), the question "did it replicate?" requires a
tolerance: *how close is close enough for two tumours to be said to agree?*

**No such tolerance exists anywhere in this repository, and none is created
here.** The two quantities that might look like candidates are not:

* `vessel_profile.MIN_BIN_COUNT = 30` is a **sample-size floor**, not an
  agreement tolerance.
* The flagship's own measured grid drift (§5.7) is a **measurement**, not a
  threshold — and it exists on one section out of three, so it cannot adjudicate
  a spread involving the other two.

So the spreads are reported and **the verdict is declined**. Naming a bar today,
with the numbers already in front of me, is the exact move that cost this
project a headline four weeks ago.

One observation that needs no threshold, because it is an ordering rather than a
comparison to a bar: **of the continuous readouts, the delivery spread's
1.2687× is the smallest by more than a factor of two** — the next-smallest is
`L` at 2.691×, and the largest is `F(1.0)` at 119.9×.

### 6.3 Two things about process that belong in the record

**(a) The closest thing to a blind statement that exists, and why it still is
not a pre-registration.** `results/verification_routes_2026-09-20.md` R12 was
written today by a thread that did not know these runs existed — the proof is
that it budgets 4–6 engineer-days to run them. Its "if it fails" clause was
therefore written **blind to this result**:

> a claim that holds on section 1 does not hold on sections 2–6: **this is the
> most valuable failure available** … a section-to-section failure would be the
> first from *tissue* variability, and it would be the honest bound on what one
> section licenses.

That is a statement of consequence written by an author who had not seen the
answer. **It is still not a pre-registration, because the data existed on disk
when it was written.** It is recorded as what it is.

**(b) The one threshold that WAS chosen after seeing this data, and the right
way it was handled.** `cab3c8e` set the reach refusal at `n_viable < 30` *after*
`xenium_io` had returned 1. Its own commit message names the hazard and the
remedy: the floor was **borrowed** from `vessel_profile.MIN_BIN_COUNT`, an
existing constant, because *"picking a fresh number here, after seeing that
`xenium_io` gave 1, would have been a threshold chosen knowing what the data
does."* That is the correct handling of this exact failure mode and it is worth
naming as precedent.

**(c) A process claim in the existing renderer that this document does not
restate.** `experiments/multisection_report.py` §7 asserts *"The rule below was
fixed before the dx = 12 µm arm had produced a single number."* That is a claim
about its author's process on 2026-08-23, not about the data, and **it is not
verifiable from anything in this checkout.** It is neither endorsed nor
contradicted here; it is simply not relied on.

---

## 7. Three defects in `experiments/multisection_report.py`, found today

The renderer's **tables** are pure interpolation and are quoted in §5 without
change. Its **prose** carries three errors. Naming them is the reason its
output was not written to `results/multisection_generalisation.md` by this pass.

### 7.1 It claims a grid control on seven sections that exists on one

Its §7.1 closes with: *"**The dominant limiter is unchanged by the grid on all 7
sections.** The cross-section flip in §4 is therefore not a mesh artefact — it
survives a 2.5× refinement in dx."*

Measured:

```
$ .venv/bin/python -c "import experiments.multisection_report as m; ..."
sections at dx30: [gse324268, visium_10x_parent_gbm, visium_gbm1, visium_gbm2,
                   visium_gbm3, xenium_319, xenium_io]
sections at dx12: ['xenium_319']
sections actually tested at both grids: ['xenium_319']
flips: []   ->  else-branch prints "all 7 sections"
```

The selection logic is right (`flips` only considers sections present at both
grids); the sentence it falls through to is wrong. **One section of seven was
solved at a second grid, and it is not one of the two that carry the flip.** The
conclusion drawn — that the flip is not a mesh artefact — **is not supported by
this data at all.** Whether the flip survives a change of mesh is untested.

This is the most consequential of the three: it is the sentence that would have
licensed reading the verdict flip as tissue.

### 7.2 Two `max/min` cells silently drop a zero — and the zero is always the flagship

`spread()` is documented as "max/min over the finite, strictly positive
entries", but the tables do not repeat that where a reader meets the number:

```
spread(hypoxic fraction over 3 Xenium) = 1.8106…   values [0.0, 0.34179…, 0.18877…]
spread(coverage_1x   over 3 Xenium) = 1.0436…      values [0.0, 0.72708…, 0.75885…]
```

So its §0 prints "**1.81×**" next to `0.00% · 34.18% · 18.88%`, and its §7
prints "**1.04**" for `therapeutic_coverage_1x` next to values that include
`0.00%`. In both cases the dropped entry is `xenium_319` — **the section that
disagrees.** The correct statement is that the ratio is undefined. §4.3 of this
document prints no ratio for either row.

### 7.3 One sentence about shipped code that stopped being true after ninety seconds

Its §2.2 says of the `n_viable == 1` hole: *"It is reported, not patched: `app/`
is outside this experiment's territory."* True at 2026-08-23 03:43:53. False at
2026-08-23 03:45:23, when `cab3c8e` patched it (§3.3a). **Publishing that
sentence today would put a false statement about shipped code into `results/`.**

### 7.4 What follows

`experiments/multisection_report.py` should not be run to its default `_OUT`
(`results/multisection_generalisation.md`) until §7.1–§7.3 are corrected. **This
pass did not run it to that path and did not create that file.** The renderer is
not this pass's territory to edit, and `results/` is a dated record — so this is
recorded here rather than fixed here.

---

## 8. What this does NOT show

1. **Nothing clinical.** Every section is an untreated resection with no outcome
   label. Every number is a prediction of a physical transport model and none
   has been scored against a measured treatment response. **Not for clinical
   use.** No response prediction, no diagnosis, no patient-level inference.
2. **Nothing about any compound but temozolomide**, and nothing comparing
   compounds. **Cross-compound ranking is withdrawn (2026-08-21), absolutely.**
   This is one compound across sections and nothing here revives that claim.
3. **Nothing from an absolute kill fraction.** It is in the payload, printed
   once in §5.2 with the payload's own `not_defensible` flag, and used for
   nothing.
4. **Nothing about spot-resolution tissue.** Four of seven sections are refused
   (§1.4). "Does it generalise across platform?" has **no answer** in this
   checkout — an absent answer, not a negative one.
5. **Nothing about any disease but glioblastoma.** All seven sections are GBM.
6. **No external anchor.** Nothing in this sweep is compared against a
   measurement of anything. It is the engine disagreeing with itself across
   tissue, which is informative about generality and says nothing about
   accuracy. It must never be presented as validation of the model.
7. **Nothing under the transport term the model lacks.** Every run has
   `transvascular_exchange = False`, the shipped default, and since 2026-08-27
   the engine itself flags temozolomide's delivery/dose verdict as one that
   **reverses** when that term is switched on (§3.3b, measured on `xenium_319`
   only). The §4.1 flip is a disagreement between two verdicts that are each
   separately known to invert.
8. **Nothing about grid convergence.** One section of seven has a second grid
   (§5.7), and the shipped caveat is unchanged: *there is no dx at which this
   pipeline has been shown to be grid-converged.*
9. **Nothing statistical.** *n* = 3 donors on the arm that carries the result
   (§4.5). No confidence interval, no p-value, no test. Three sections agreeing
   is three sections agreeing.
10. **Nothing about a section with no endothelial marker.** All seven carry one,
    so `vessel_profile.direct_vessel_distance_um`'s refusal path is **untested
    by these specimens, not passed by them.**
11. **Nothing that separates tissue from panel on `xenium_io`** (§2.2). That
    section's endothelial class — the source term of the whole delivery solve —
    is called by a different gene set from the other two.
12. **Nothing about the hypoxia field's correctness.** The verdict flip is
    driven by it, and the payload's own caveats rate it `rho = -0.1523` against
    measured markers and label the cell-state route **"NOT A MEASUREMENT"** with
    a 165× swing under bookkeeping choices (§4.1).

---

## 9. Closing calibration guard

```
$ .venv/bin/python -c "from tmesim.validation.calibration import calibrate_alpha; print(repr(calibrate_alpha()))"
END calibrate_alpha() = 22.223477711172563
ASSERT OK
```

Identical to the opening value. The same assertion was made and passed at the
start and the end of every re-run in §3, and every stored entry carries
`alpha_at_run = 22.223477711172563`.

**No code under `tmesim/`, `experiments/`, `app/`, `docs/`, `landing-3d/` or
`data_cache/` was modified in this pass, and no file in `results/` or `archive/`
other than this new one was created or touched.**
