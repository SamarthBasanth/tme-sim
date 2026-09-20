# Orphaned results index — work already done that nothing in the repo points at

**Date: 2026-09-20. Status: INDEX / REFERENCE AUDIT. No new computation, no new
measurement, no verdict of my own.** Every number below is read out of a file
that already existed when this pass started, and every one is attributed to the
file it came from. Nothing under `tmesim/`, `app/`, `docs/`, `experiments/`,
`landing-3d/`, `render.yaml` or `Dockerfile` was touched. **No existing file in
`results/` or `archive/` was edited, reworded or re-scored. This file is the
only thing created.**

**Calibration guard, start of pass:**

```
$ .venv/bin/python -c "from tmesim.validation.calibration import calibrate_alpha; a=calibrate_alpha(); assert a==22.223477711172563, a; print('alpha OK',a)"
alpha OK 22.223477711172563
```

*(Closing assertion in §6.)*

---

## 0. Why this file exists, and what it is not

This project's recurring defect is claim-surface drift, and it has so far been
tracked in one direction: documents claiming more than the engine supports.
[`results/claim_surface_audit_2026-08-30.md`](claim_surface_audit_2026-08-30.md)
and [`results/withdrawal_banner_log_2026-09-03.md`](withdrawal_banner_log_2026-09-03.md)
are both passes in that direction.

This is a pass in the other direction: **completed, reproducible work on disk
that no document points at.** The cost of that is different but real — it makes
the same work get done twice, and it makes the repository's own self-assessment
worse than its evidence.

**This index does not promote anything.** An orphaned result is not
automatically a good result. Three of the items below are *failures*, one is a
*contradiction the owner has explicitly reserved*, and the largest one is a
result whose headline is that the product's central verdict **does not hold off
its flagship section**. Surfacing them means surfacing that too.

---

## 1. The primary orphan — `results/multisection_generalisation.json`

### 1.1 What it is, and its provenance

| | |
|---|---|
| path | `results/multisection_generalisation.json` |
| size / mtime | 129,325 bytes, **2026-08-23 03:43:53** |
| internal stamp | `"updated": "2026-08-23 03:43:53"` |
| `"alpha"` | `22.223477711172563` — and `alpha_at_run` carries the same value on **all eight** entries |
| `"drug"` | `"temozolomide"` — **one compound, no comparison** |
| driver | `experiments/multisection_generalisation.py` (tracked) |
| renderer that was written for it | `experiments/multisection_report.py` (838 lines, tracked) |
| committed in | `bf15977`, 2026-08-24 23:34:53 −0500, 9 files / +7,867 lines |

Commit `bf15977` added the data and the drivers. **It added no `results/*.md`
for the multisection sweep at all.** The renderer's own `_OUT` target,
`results/multisection_generalisation.md`, was never produced.

**Eight entries: seven GBM sections solved at dx = 30 µm, plus one dx = 12 µm
replicate of the flagship.** All eight carry `status: "ok"`.

| entry | platform | `n_cells` |
|---|---|---:|
| `xenium_319@dx30` (flagship) | Xenium | 40,887 |
| `xenium_io@dx30` | Xenium | 816,769 |
| `gse324268@dx30` | Xenium | 851,509 |
| `visium_10x_parent_gbm@dx30` | Visium (spot) | 41,560 |
| `visium_gbm1@dx30` | Visium (spot) | 57,050 |
| `visium_gbm2@dx30` | Visium (spot) | 34,707 |
| `visium_gbm3@dx30` | Visium (spot) | 21,279 |
| `xenium_319@dx12` | Xenium | 40,887 |

### 1.2 The question it bears on, and where that ranking is written down

[`results/CLAIM_INVENTORY.md`](CLAIM_INVENTORY.md) **§6, "What would most change
this document"** (line 501), *"In rough order of how much each would move the
honest assessment"*:

> 1. **A second section, same tumour type.** Almost every defensible claim above
>    is *n* = 1 in tissue. This is the cheapest, largest improvement available.

**A precision note, because it matters for how loudly this may be stated.** That
list is numbered **0 through 5**, and item **0** — *"Reconciling `vivum validate`
with `kill_calibration_breast.md`"* — is printed above it. So "a second section"
is the **second entry in the list, numbered 1**.
[`results/verification_routes_2026-09-20.md`](verification_routes_2026-09-20.md)
R12 (line 811) nonetheless describes it as the item CLAIM_INVENTORY §6 *"ranks
**first** among everything that would change the honest assessment"*, and R12 is
itself ranked **4th** in that file's own priority table (line 954). **The repo
says "first" in one place and prints it second in another. I am not resolving
that; I am recording that the claim "the top-ranked open question" is the
repository's own phrasing and is not literally what §6's numbering shows.**

What is not ambiguous is that this is a named, standing, high-priority open item
in at least three documents, and that **four sections of replication were sitting
on disk unread for 28 days** —
[`results/verification_routes_2026-09-20.md`](verification_routes_2026-09-20.md)
line 840: *"it is currently costing four sections of replication for free."*

### 1.3 What it concludes — including what it fails

Read directly from the payload. **No mean over the seven appears here, and none
should: the sections disagree, and the disagreement is the finding.**

**(a) The per-cell verdict the product exists to answer takes three different
values across the seven sections.** `mech_dominant_limiter`, with
`mech_thresholds` byte-identical on every entry:

| section | limiter | share |
|---|---|---:|
| `xenium_319` (flagship) | **`dose`** | 0.9985 |
| `xenium_io` | **`resistance`** | 0.6090 |
| `gse324268` | **`resistance`** | 0.6629 |
| `visium_gbm1` ⚠︎spot | `dose` | 0.6781 |
| `visium_10x_parent_gbm` ⚠︎spot | `chance` | 0.5438 |
| `visium_gbm2` ⚠︎spot | `chance` | 1.0000 |
| `visium_gbm3` ⚠︎spot | `chance` | 1.0000 |

`xenium_319` and `gse324268` are typed by the **same** marker map and still
disagree. **The flagship is the only section reporting `dose`, and it is the
section the README, the landing copy and most of `results/` are written around.**

**(b) The delivery-spread ratio is the one continuous readout that stays close.**
`E_p90/E_p10` on the three single-cell sections: **1.9936, 2.3211, 2.5293** —
a max/min of **1.2687×** across three sections, which the registry records as
three different donors.

**(c) Almost everything else moves.** Fitted `penetration_length_um` on the same
three: **277.69, 133.16, 103.18 µm** (2.691×), with the **flagship the longest**
— i.e. the shipped specimen is the most drug-penetrable of the three.
`therapeutic_coverage_1x`: **0.00%, 72.71%, 75.89%**. `hypoxic_tumor_fraction`:
**0.00%, 34.18%, 18.88%**.

**(d) What it refuses.** Four of the seven sections are **Visium spot data** and
carry a `section_caveats` entry, written by the engine, that says per-cell
spatial statistics are invalid there *and* that **the delivery map is invalid
there too**, because the exposure field's source term is a density over
endothelial positions that are synthesised as 12 pseudo-cells per spot:

> *"Read NOTHING on this section as a statement about where its vessels are."*

**So "seven sections" is four refusals and a three-section arm.** Any use of the
number 7 that does not say so is an overstatement.

**(e) What it does not produce.** `kill_fraction` is present on every entry and
every entry also carries `kill_fraction_not_defensible`. It is not a defensible
output and is used for nothing here. `reach_headline_F` at `x = 0.1` is **0.0 on
six of the eight entries** — stable at the value that makes it useless.

### 1.4 Does anything supersede it? — **the withdrawn claim is NOT in it**

**This was the first thing checked, and it is the good news.**
`"drug": "temozolomide"`, one compound, on every entry. **There is no
cross-compound ordering anywhere in the file**, so the 2026-08-21 withdrawal
(`results/pk_ranking_robustness.md`, `results/capacity_ranking_robustness.md`;
banners logged in
[`withdrawal_banner_log_2026-09-03.md`](withdrawal_banner_log_2026-09-03.md))
does not touch it, and nothing in it revives that claim.

**Two things since 2026-08-23 do bear on it, and both are already recorded in
[`results/multisection_generalisation_2026-09-20.md`](multisection_generalisation_2026-09-20.md):**

1. **`cab3c8e`, 2026-08-23 03:45:23 — ninety seconds after the JSON was
   written.** `app/server/reach.py` now refuses the reach ladder when
   `n_viable < 30`. `xenium_io` has `n_viable = 1`. **Its stored reach block is
   superseded in code and must not be quoted from the JSON.** That sibling
   document re-solved four entries against today's tree and found 48/48 fields
   identical on three of them and 44/48 on `xenium_io`, the four differences all
   being *a number becoming a refusal*.
2. **`962b53a`, 2026-08-27.** `verdict_stability('temozolomide')` now returns
   `stable: False` — the compound this whole sweep is built on **reverses** its
   delivery/dose verdict when the missing transvascular exchange term is
   switched on (`results/verdict_inversion_census.md`, measured on `xenium_319`
   only). **The §1.3(a) flip is therefore a disagreement between two verdicts
   each separately known to invert under a term that ships OFF.**

### 1.5 Its current linkage state — the orphan was read 19 minutes before this pass

**This is the part of the brief that has already changed under it, and saying so
is the point of an index.** As of **2026-09-20 03:59:53**, the JSON is no longer
unread: [`results/multisection_generalisation_2026-09-20.md`](multisection_generalisation_2026-09-20.md)
(951 lines) is a full seven-section reading of it, correctly labelled
*"POST-HOC REPORT OF AN EXISTING COMPUTATION. NOT a pre-registered test."* It
is honest, it quarantines the four spot sections, it marks the `xenium_io` reach
block superseded, and its §8 lists twelve things the result does not show.

**But that document is itself an orphan.** A reference sweep of every text file
in the tree (excluding `.git`, `.venv`, `data_cache`) finds **no file anywhere
references it**. The same is true of three other files written today. So the
orphan has been *read*; it has not yet been *linked*.

**Two inaccuracies in that sibling document, recorded here because I may not
edit it:**

* Its **§4.4 item 4** states *"**`F(0.1)` = 0 everywhere**, single-cell and spot
  alike."* The payload disagrees on two of eight entries:
  `visium_10x_parent_gbm@dx30` has `reach_headline_F = 0.001810926949048496`
  and `xenium_319@dx12` has `0.00049800796812749`. Its own **§5.5 table prints
  the first of these correctly as `0.001811`**, so the prose summary overstates
  a stability its own table does not show. The direction of the error is small
  and toward tidiness, not toward the product.
* Its §0 quotes CLAIM_INVENTORY §6 as ranking the second-section item **first**.
  See §1.2 above: §6's list is numbered from 0 and that item is numbered 1. The
  phrasing is inherited from `verification_routes_2026-09-20.md` R12, not
  invented there.

### 1.6 Where it should be referenced from

None of these are edits I made; they are the places a reader would look and find
nothing. **`results/` and `archive/` are dated published records — the right
mechanism for the two `results/` entries below is a dated banner or a new dated
file, never a rewrite of the body.**

| surface | why it belongs there |
|---|---|
| `results/CLAIM_INVENTORY.md` §6 item 1, and the *n* = 1 rows in §1a/§1c | the item this data bears on; the *n* = 1 rows are the ones it tests |
| `results/verification_routes_2026-09-20.md` R12 | R12's "Data. Already on disk" paragraph is about exactly these sections |
| `README.md` | carries the flagship-section verdict as the product's headline; §1.3(a) is the direct qualifier |
| `results/heldout_patient.md`, `results/flagship_markers.md` | both argue about what transfers between sections; both already carry 2026-09-03 withdrawal banners for a different reason |
| the Studio panel / `app/` caveat surface | the run payload already ships `verdict_stability`; the section-dependence is the same class of disclosure |

---

## 2. The `p90_capacity_sweep` contradiction — exact state, **NOT RESOLVED HERE**

**The user has reserved this one. Nothing below changes either file. This is a
transcription of what each says.**

| | `results/p90_capacity_sweep.md` | `results/p90_capacity_sweep.json` |
|---|---|---|
| mtime | **2026-08-23 03:34:52** | **2026-08-24 23:33:50** |
| length | 186 lines | 20,680 bytes |
| says | **"Status: PRE-DECLARED, NOT YET RUN."** | a completed run, `elapsed_s` 1232.72 |
| ends | `*(Everything below this line is written after the sweep. The pre-declaration above is frozen.)*` — **and nothing follows it** | `"verdict": "THE HEADLINE IS CAPACITY-SENSITIVE AND SO IS THE CLAIM BUILT ON IT"` |

**Both files were added in the same commit, `bf15977`** — and that commit's own
message opens *"The headline is capacity-sensitive, and the flagship section is
not representative"* and prints the R-ladder in full. **So the verdict the `.md`
says does not exist yet is stated in the commit message that added the `.md`.**

What the JSON records, against the criteria the `.md` fixed in advance:

* **Reproduction check (§0.4, declared as a check, not a criterion): PASSES.**
  R = 1 measured `2.1106584343146793` / `3.942472241633886` against the published
  `2.111` / `3.942`.
* **C1 — value stability at R = 39: `fail: true`.** Temozolomide moves
  `2.1107 → 3.2969`, an absolute move of **1.1862 = 5.872× its own 0.202
  nuisance width**; doxorubicin moves `3.9425 → 4.9580`, **1.0155 = 4.494× its
  0.226**.
* **C2 — separation floor 1.74 over R = 3.2–4184: `fail: true`**, minimum
  **1.5039** at **R = 39**.
* **C3 — descriptive, no bar, declared as descriptive:** over the measured
  capacity range the ratio spans **2.2012 – 24,526.69** for one compound and
  **3.9854 – 226,391.07** for the other.
* All five planted defects bit (`plant.all_bit: true`); the two routes agree at
  convection OFF to `rel_diff` 5.77e−15 and **disagree at convection ON** by
  `rel_diff` 0.02673, which §0.5/§0.6 declared in advance as the planted trap.

**Two pre-declared predictions in §0.2 were wrong, and the JSON is the only
place that is visible:**

1. §0.2 predicted the R = 39 value would be *"at least 2x its R = 1 value"*.
   Measured **1.562×**. The commit message records this explicitly: *"the
   prediction was WRONG in the conservative direction and is recorded as wrong."*
2. §0.2 predicted the two compounds would **converge**, the separation *"to
   approach ~1 at R = 4184"*. The JSON records separation **rising** to
   **9.2304** at R = 4184. **Neither the `.md` nor the commit message records
   this second falsification anywhere.**

**One caveat that belongs on any future write-up of C2, stated and not acted
on.** C2 is a **two-compound separation**, and the mechanism that failed it — no
per-compound tissue-capacity term, uniform R breaking the result below the
best-available constant of 39 — is the *same* mechanism
`results/capacity_ranking_robustness.md` gave for the **2026-08-21 withdrawal of
cross-compound ranking**, which is absolute. C1, by contrast, is a
single-compound value-stability result and stands on its own. Separately:
`results/delivery_spread_noise.md`, the file C2's 1.74 floor is inherited from,
states a two-compound separation in its own headline and **does not carry a
2026-09-03 withdrawal banner** (`grep -c` for it in
`withdrawal_banner_log_2026-09-03.md` returns **0**; its only `withdraw` hit is
an internal note at line 92 withdrawing *one argument* for P90/P10, not the
metric). **Whether that file is inside or outside the withdrawal's scope is a
judgement I am flagging, not making.**

---

## 3. The rest of the orphan sweep

**Method.** Every file under `results/` and `archive/` (474 files) was searched
for by basename *and* by repo-relative path across every `.md`, `.py`, `.json`,
`.html`, `.ts/.tsx/.js`, `.yaml/.yml`, `.txt`, `.toml`, `.css`, `.mjs`, `.sh`
file in the tree, excluding `.git`, `.venv`, `node_modules`, `.pytest_cache`,
`data_cache` and `__pycache__`. **76 had zero references.** They are not all
equal, and the classification below is the useful part.

### 3.1 Genuinely unused finished documents

| file | date | what it is | where it belongs |
|---|---|---|---|
| [`multisection_generalisation_2026-09-20.md`](multisection_generalisation_2026-09-20.md) | 09-20 | the §1 reading of the seven-section sweep | CLAIM_INVENTORY §6, verification_routes R12, README |
| [`visium_section_level_recovery.md`](visium_section_level_recovery.md) | 09-20 | *"Is the spot-data refusal over-applied?"* — a jitter-reseed audit of the four Visium GBM sections, i.e. the direct follow-up to §1.3(d) | wherever the spot refusal is stated; it is the audit of that refusal |
| [`oxygen_sink_form_2026-09-20.md`](oxygen_sink_form_2026-09-20.md) | 09-20 | pre-registered A/B, **REFUTED on the registered threshold**; `sink_form` still ships `"linear"` | the oxygen-field caveat chain; it is a clean registered negative |
| [`molecular_compile_audit_2026-09-20.md`](molecular_compile_audit_2026-09-20.md) | 09-20 | re-derivation of the two-route law and the BYO surface | `MOLECULAR_CONTRACT.md`, `two_route_transport.md` |
| [`first_user_journey_2026-09-20.md`](first_user_journey_2026-09-20.md) | 09-20 | a real-browser newcomer walkthrough against a live `app.server` | `first_user_readiness.md` / `ship_readiness.md` lineage |
| [`security_audit.md`](security_audit.md) | 08-15 | read-only pass over `app/server`, `app/web`, `tmesim`, `landing-3d`, deps, git history | nothing references it; it is the only security pass on disk |
| [`studio_a11y_audit.md`](studio_a11y_audit.md) | 07-28 | keyboard/screen-reader pass driven headlessly against the real backend | same |
| [`landing_design.md`](landing_design.md), [`rebrand_landing.md`](rebrand_landing.md) | 08-15 | landing-page design and the Vivum→Codon rename pass | historically interesting; **note the name reverted to Vivum on 2026-08-22**, so `rebrand_landing.md` describes a rename that was undone |
| [`REALITY_CHECK_SUMMARY.md`](REALITY_CHECK_SUMMARY.md) | 07-26 | v2-era summary; calls the engine *"a first validation against measured reality"* | **do not promote.** Its framing predates every 2026-08 withdrawal and does not match the current claim surface |

### 3.2 Near-orphans worth naming

Referenced by exactly one file, where that one reference is fragile:

* [`REAL_VALIDATION.md`](REAL_VALIDATION.md) — reached only from
  `grid_dependence_audit.md`.
* [`withdrawal_banner_log_2026-09-03.md`](withdrawal_banner_log_2026-09-03.md) —
  **the log of the withdrawal banners is itself reachable only from
  `stale_blockers_2026-09-20.md`**, a file written 17 days later.
* [`export_honesty.md`](export_honesty.md), [`hypoxia_imputation_model.md`](hypoxia_imputation_model.md),
  [`conservation_audit.md`](conservation_audit.md),
  [`cell_state_multipliers_2026-09-20.md`](cell_state_multipliers_2026-09-20.md)
  — each reachable only from a single `app/` or `tests/` source file, not from
  any document.

### 3.3 Zero-reference files that are **not** orphaned work

Recorded so the count above is not read as 76 unused results:

* `results/sweep/*_partial.json` (9 files) — **restart checkpoints by design**;
  `area_conservation_regeneration.md` §§ and `subgrid_closure.md` describe the
  checkpoint protocol, and the completed non-`_partial` siblings *are*
  referenced.
* `results/sweep/descriptor_run_*.json` (5) — referenced **by pattern**:
  `descriptor_falsification.md:395` says *"Records land in
  `results/sweep/descriptor_run_<stage>.json`"*.
* `results/sweep/retard_*.json`, `tworoute_*.json`, `*_records_dx12.jsonl` — the
  `*_prereg.json` members of each family are cited by name in
  `intracellular_retardation.md` and `abraham_citation_erratum.md`; the raw run
  records are intermediates.
* `archive/synthetic_era/*.py` (15) — the pre-real-data era, correctly archived.
* Five `results/demo_*.png` / `molecular_simulation.png` — v1 demo art.

**One family is a real gap:** `results/sweep/gbmspace_batch_{A,B,C,D,followup}_2026-07-29T05-54-00Z.md`
(5 files) and nine `vasc_bandwidth_*` GBM sweep records have **no reference by
name or by pattern** anywhere, while their breast-tissue sibling
`vasc_bandwidth_breast_default_vs_15um_...md` *is* cited by
`vascular_validation_breast.md:174` and `experiments/vascular_validation_breast.py:678`.
The GBM half of that same sweep is unlinked.

---

## 4. What this index does NOT claim

1. **Nothing clinical.** Nothing here is a response prediction, a diagnosis, a
   treatment-planning output, a "digital twin" or personalised medicine. Every
   section indexed is an untreated resection with no outcome label. **Not for
   clinical use.**
2. **Nothing about compound ordering.** Cross-compound ranking is **withdrawn,
   2026-08-21, absolutely**, and nothing in this file revives it. §1's data is
   one compound. §2 reports the *failure* of a two-compound separation and flags
   its scope question rather than asserting any order.
3. **No absolute kill fraction.** The field exists in the §1 payload with the
   engine's own `kill_fraction_not_defensible` flag attached and is used for
   nothing.
4. **Nothing is externally anchored.** §1 is the engine disagreeing with itself
   across tissue. That is informative about generality and says **nothing about
   accuracy**, and it must never be presented as validation of the model.
5. **Nothing here re-scores or re-runs anything.** No solve was executed in this
   pass. Every figure is a transcription, and where two files disagree I have
   quoted both rather than picking one.
6. **§2 is not resolved.** The `.md`/`.json` contradiction is recorded exactly
   as it stands, by the owner's explicit instruction, and neither file was
   touched.
7. **Being orphaned is not a quality judgement.** §3.1 includes a refuted
   pre-registration and a document (`REALITY_CHECK_SUMMARY.md`) that should
   **not** be linked from anywhere current.

---

## 5. The one-line version

**The seven-section GBM sweep was solved on 2026-08-23, reproduces bit-for-bit
today on three of four re-checked entries, contains no withdrawn claim — and its
headline is that the product's per-cell verdict takes three different values
across seven sections, with the flagship the only one giving the answer the
repository is written around.** It was read for the first time on 2026-09-20,
into a document that nothing yet references. And `results/p90_capacity_sweep.md`
still says **"NOT YET RUN"** over a JSON that records a completed run failing
both of its own pre-declared criteria, added in the same commit as the `.md`, by
a commit message that states the verdict.

---

## 6. Closing calibration guard

```
$ .venv/bin/python -c "from tmesim.validation.calibration import calibrate_alpha; a=calibrate_alpha(); assert a==22.223477711172563, a; print('alpha OK',a)"
alpha OK 22.223477711172563
```

Identical to the opening value. **No file under `tmesim/`, `app/`, `docs/`,
`experiments/`, `landing-3d/` or `tests/` was modified in this pass, and no file
in `results/` or `archive/` other than this new one was created or touched.**
