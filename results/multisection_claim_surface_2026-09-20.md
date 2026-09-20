# CLAIM-SURFACE AUDIT against `results/multisection_generalisation.json` — what the product ships, and what three sections support

**Status: AN AUDIT OF SHIPPED TEXT, NOT A NEW COMPUTATION.** Nothing was
simulated for this file. It reads
[`results/multisection_generalisation.json`](multisection_generalisation.json)
(computed 2026-08-23, read for the first time 2026-09-20) and its report,
[`results/multisection_generalisation_2026-09-20.md`](multisection_generalisation_2026-09-20.md)
— **including that file's own correction block, which is where two of its
summary sentences are withdrawn** — and asks one question of every surface a
reader actually meets: *does anything the product claims, implies or relies on
require a generalisation across tissue sections that this result does not
support?*

**Calibration guard, start of pass:**

```
$ .venv/bin/python -c "from tmesim.validation.calibration import calibrate_alpha; a=calibrate_alpha(); assert a==22.223477711172563, a; print('alpha OK',a)"
alpha OK 22.223477711172563
```

*(The closing assertion is in §7.)*

---

## 0. The result this audit is against, and its bounds

Read the bounds before any verdict below, because every verdict is conditional
on them.

* **One compound.** Temozolomide. Nothing here touches any other compound, and
  **cross-compound ranking remains withdrawn (2026-08-21), absolutely.**
* **One configuration.** `dx = 30 µm`, 6 h window, `sensitivity_mode = uniform`,
  `vascular_mode = discrete`, `oxygen_mode = pde`, `state_mode = oxygen`,
  `transvascular_exchange = False`. `mech_thresholds` byte-identical on every
  row.
* **Three sections, not seven.** Four of the seven are Visium spot sections
  whose own payload attaches a blocking `spatial_spread: NO` caveat and states
  that the delivery map is invalid there. **This is a three-section arm plus
  four refusals.** No figure in this file is computed over seven sections, and
  the seven-section spreads that exist in the source are a quarantined
  counterfactual there.
* **n = 3 donors.** Three is not a sample. No statistical claim is made here and
  none is available.
* **No second grid where it matters.** Only `xenium_319` was solved at a second
  grid. **Neither section carrying the verdict flip was**, so whether the flip
  survives a change of mesh is UNTESTED.
* **No external anchor.** This is the engine disagreeing with itself across
  tissue. It is informative about generality and says nothing about accuracy,
  and it must never be presented as validation.

The three single-cell readouts this audit turns on, read out of the JSON:

| | `xenium_319` | `xenium_io` | `gse324268` |
|---|---:|---:|---:|
| dominant limiter | **`dose`** | **`resistance`** | **`resistance`** |
| its share | 0.998534 | 0.608967 | 0.662876 |
| `delivery` share | 0.001466 | 0.001817 | 0.000098 |
| `E_p90/E_p10` | 1.993572 | 2.321118 | 2.529296 |
| fitted `L` (µm) | 277.690 | 133.162 | 103.184 |
| `hypoxic_tumor_fraction` | 0.000000 | 0.341792 | 0.188772 |
| `therapeutic_coverage_1x` | 0.000000 | 0.727088 | 0.758856 |
| `median_E_over_ic50` | 0.305531 | 1.179037 | 1.323738 |
| `F(0.1)` | 0 | 0 | 0 |

Paired grid control, the only one that exists — `xenium_319` at dx = 30 against
dx = 12: `E_p90/E_p10` 1.993572 → 1.994409 (**1.0004×**), `L` 277.690 → 278.492
(**1.003×**), limiter `dose` → `dose`.

`xenium_319` and `gse324268` are typed by the **same** marker map (`brain_319`)
with the **same** one-gene `PECAM1` endothelial gate, and they disagree. The
flip is therefore **not** a marker-panel artefact. `xenium_io` is typed by
`io_480` and every three-way spread including it carries that caveat.

---

## 1. The audit table

`SUPPORTED` = the three-section arm backs it. `OVERSTATED` = it implies
cross-section generality this result contradicts. `UNDERSTATED` = this direction
is real and the copy hedges more than it needs to. `UNSETTLED` = named, not
resolved.

### 1a. `README.md` — OWNED, edited by this pass

| # | Claim, as it read | Verdict | Disposition |
|---|---|---|---|
| R1 | *"GBM's hypoxic fraction is 0.0000 at every grid spacing tested under the shipped configuration"* (standing limitations) | **OVERSTATED** | It is a property of `xenium_319`, not of GBM. The same shipped configuration at dx = 30 gives 34.18 % and 18.88 % on two other **registered** sections. Rewritten. |
| R2 | *"That drug-dependence is inert **today only because the GBM section's hypoxic fraction is 0.0000 at every grid spacing tested** — a non-empty hypoxic band would activate it"* (contingent, `TMZ_HYPOXIA_RESISTANCE`) | **OVERSTATED** | The trigger it names as hypothetical **has already fired** on two registered sections, and on both of them `resistance` is the dominant limiter. Rewritten to say the near-invariance is untested there. |
| R3 | *"(1) a **second section of the same tumour type** — almost every defensible claim above is n = 1 in tissue, and this is the cheapest large improvement available"* | **OVERSTATED (stale)** | It was bought on 2026-08-23 and the headline verdict did not survive it. Rewritten to say what was bought, what it returned, and what is still unbought. |
| R4 | *"did the drug never arrive, or did it arrive and stay too weak? … Vivum separates them, cell by cell"* (who it is for) | **OVERSTATED** | Presents the output as binary. `resistance` is a shipped verdict and it is the answer on 2 of the 3 sections. A third answer added, pointing at the new contingent entry. |
| R5 | *"the defensible output is a **ratio inside one map, between regions of one tissue section, for one compound**"* (honest summary) | **SUPPORTED** | Correctly scoped; the result strengthens it rather than contradicting it. Unchanged. |
| R6 | *"The spatial-attribution statistic is **refused in code** on spot data"* | **SUPPORTED** | Demonstrated on four Visium sections in this sweep, refusal firing as documented. Unchanged. |
| R7 | *"Absolute kill fraction is not a defensible output"* | **SUPPORTED** | The payload carries `kill_fraction_not_defensible` on every entry and no comparison anywhere is built on it. Unchanged. |
| R8 | *"One section per tumour, no replicates."* | **SUPPORTED** | Still true of every **scored** claim: no scored anchor has been re-run on the other two sections. Unchanged. |
| R9 | `E_p90/E_p10` behaviour across tissue | **UNDERSTATED — nothing said at all** | Added: 1.9936 / 2.3211 / 2.5293, spread 1.2687×, smallest continuous readout in the sweep by more than 2×, and 1.0004× across a 6.2× refinement in voxel count. |

### 1b. `landing-3d/` — OWNED, edited by this pass

| # | Claim | Verdict | Disposition |
|---|---|---|---|
| L1 | beat 06 body: *"get back a per-cell verdict — **delivery or potency** —"* | **OVERSTATED** | It named the one setting the model has never produced non-degenerately and omitted one it produces routinely. On the majority of this arm the verdict a user gets back is **neither half**. The phrase was DELETED (see §2 for why deleted rather than replaced). |
| L2 | beat 06 ledger, *"The 03 verdict · **INVERTS TWO WAYS**"* | **OVERSTATED (a count that is now wrong)** | A third route exists and it is the only one that is not a modelling choice: the section. Now `INVERTS THREE WAYS`, with the clause *"And so does another section of the same tumour."* |
| L3 | beat 03 body: *"**On the flagship section** 99.8% of surviving tumour cells are under-dosed…"* | **SUPPORTED** | Correctly scoped, and the figure reproduces: `mech_dominant_share` = 0.998534 on `xenium_319` at dx = 30 and 0.997824 at dx = 12. Unchanged. |
| L4 | beat 04: *"**In one section**, tumor cells in the macrophage-densest quintile…"* | **SUPPORTED** | Scoped, and untouched by this sweep — the immune sink term is in no shipped run. Unchanged. |
| L5 | beat 05 ledger, *"Spot platforms · NOT REPORTED"* | **SUPPORTED** | Now demonstrated on four sections. Unchanged. |
| L6 | beat 05 ledger, *"Absolute kill fraction · NOT DEFENSIBLE"* / *"Patient response · OUT OF SCOPE"* | **SUPPORTED** | Unchanged. |
| L7 | beat 06 ledger, *"“Arrives almost everywhere” · 2.2–3.0× TOO SHALLOW"* | **UNDERSTATED — and it stayed understated** | A second, independent reason the sentence flatters the drug is measured: of three sections the flagship is the most penetrable (277.69 against 133.16 and 103.18). The clause was written, measured, and **withdrawn for room**. Recorded in the file as a known omission and as the first thing that beat should buy. See §2. |
| L8 | guardrail *"A DELIVERY VERDICT HAS NEVER BEEN DEMONSTRATED NON-DEGENERATELY"* | **SUPPORTED — and now on more evidence** | `delivery` is dominant on **none** of the seven sections and never exceeds a 0.001817 share, while `dose` and `resistance` each dominate somewhere. Recorded in the guardrail; it is a reason to keep saying "two", not to promote the third. |
| L9 | `app/layout.tsx` `<meta description>` | **SUPPORTED** | No cross-section claim, no ordering, no "validated". The 2026-09-20 removal of the comparing-candidates clause is intact. **Unchanged — checked, not assumed.** |

### 1c. `app/server/` — NOT OWNED, reported in §3

| # | Claim | Verdict |
|---|---|---|
| S1 | `jobs.py` `verdict_stability()` — *"KEEPS its delivery/dose verdict"*, shipped into `extra_metrics` of **every** run document | **OVERSTATED** — the flag is measured on `xenium_319` only, carries `measured_on` but says nothing in its own prose, and the UI prints it verbatim on any section. |
| S2 | `report.py` — *"cells that survived because the drug never arrived … not because they resisted it. **This is the answer.**"* | **OVERSTATED** — "resisted it" is the dominant verdict on 2 of 3 sections. |
| S3 | `reach.py` `CONVENTION_NOTE` — the promotion of the delivery spread to headline, justified on the flagship alone | **UNDERSTATED** — three donors, 1.2687×, and 1.0004× across grids. This is the clearest understatement on any surface. |
| S4 | `reach.py` module docstring — *"**The headline** — `F(0.1)`, picked off that ladder"* vs `CONVENTION_NOTE`'s *"It is NO LONGER HEADLINED ANYWHERE"* | **OVERSTATED (internal contradiction)** — and the sweep makes it material: `F(0.1)` is 0 on all three single-cell sections. |
| S5 | `cross_section.py` `NO_CROSS_SECTION_NOISE_FLOOR` — *"it was measured on ONE section with the tissue held fixed"* | **UNDERSTATED** — still literally true of `delivery_spread_noise.md`, but a three-donor figure now exists and is not cited. |
| S6 | `cross_section.py` module design (one compound, section varies, spot refusal carried verbatim, one RunSpec signature per axis) | **SUPPORTED** — this module is the one surface built for exactly this question and it is correct. |

### 1d. `app/web/` — NOT OWNED, reported in §3

| # | Claim | Verdict |
|---|---|---|
| W1 | `js/ui/finding.js` — *"The claim above says drug is present essentially everywhere in the tumour and is essentially everywhere below the concentration that produces half its effect."* rendered **unconditionally** | **OVERSTATED — the worst instance found.** On a `resistance` section the claim above says no such thing. |
| W2 | `js/ui/finding.js` — *"THIS VERDICT HOLDS FOR …"* / *"THIS VERDICT REVERSES FOR …"* printed with no section qualifier | **OVERSTATED** — inherits S1. |
| W3 | Fixed binary headings — *"Delivery or dose"*, *"was it delivery, or was it dose?"*, *"the delivery-versus-potency split"* | **OVERSTATED** — fixed prose that renders identically on a `resistance` section. |
| W4 | `js/ui/doseinput.js` — *"for every tissue and every vessel geometry"* | **SUPPORTED, conditionally** — gated on `plasmaBoundHolds === true`, where it follows by arithmetic. This sweep ran the one-way law, so the branch does not apply to it. **Not a defect; checked and cleared.** |
| W5 | The per-cell data path (`LIMITER_LABEL`, `SURVIVAL_PHRASE`, `RESISTANCE-LIMITED`, dominant-share sentence) | **SUPPORTED** — emits all four outcomes correctly. Only the surrounding fixed copy is binary. |
| W6 | Tissue registration language | **SUPPORTED** — *"Nothing is uploaded from this browser for tissue."* The path convention holds across `app/web/`, `docs/` and the contract. |

### 1e. `docs/` — NOT OWNED, reported in §3

| # | Claim | Verdict |
|---|---|---|
| D1 | `architecture.md` — *"Earlier attribution of small predicted kills to hypoxic resistance was wrong — they are driven mainly by limited delivered exposure versus these drugs' IC50s."* | **OVERSTATED — the strongest single contradiction in `docs/`.** A general attribution claim, false on 2 of 3 sections. |
| D2 | `architecture.md` — *"The hypoxic fraction is grid-sensitive and only order-of-magnitude."* | **OVERSTATED by omission** — it attributes the spread to mesh. Holding dx fixed at 30 the spread is 0.00 / 34.18 / 18.88 %, i.e. **donor**, not mesh. |
| D3 | `MATHEMATICAL_MODEL.md` — *"the predicted hypoxic tumor fraction is 0.0 at every Δx"* | **OVERSTATED** — a property of one section, written as the model's behaviour. |
| D4 | `VERIFICATION_METHODOLOGY.md` — *"under the current shipping default the HYPOXIC state is empty at every Δx, so the hypoxia-resistance multiplier has no cells to act on"* | **OVERSTATED** — attributed to the shipping **default**, not to the section. It has cells to act on, on two registered sections. |
| D5 | `WORKSPACE_CONTRACT.md` — *"Tissue registers BY FILESYSTEM PATH. The word 'upload' never appears about tissue."* | **SUPPORTED** — the contract holds. |

---

## 2. The edits this pass made, and what they cost

### 2.1 `README.md`

1. **Who it is for.** Added the third answer: *"and reports a third answer,
   `resistance` (exposure was adequate and the cell survived anyway), when that
   is what the run says. On two of the three registered single-cell GBM
   sections, `resistance` is what it says."*
2. **What is contingent — a new first entry:** *"The per-cell verdict is
   SECTION-DEPENDENT, and the flagship is the outlier."* Carries the RunSpec,
   the three limiters and their shares, the same-marker-map point, the
   penetration lengths, and the bounds (three-section arm plus four refusals,
   n = 3, one compound, one configuration, no second grid on the flipping
   sections).
3. **What is contingent — a second new entry:** *"What SURVIVED the same test,
   and it is the delivery spread."* This is the understatement correction, and
   it is stated at exactly the strength the evidence carries: reproducibility of
   magnitude across three donors, **not** evidence that the statistic
   discriminates between tissues.
4. **`TMZ_HYPOXIA_RESISTANCE` inertness** narrowed to the flagship, with the
   two measured hypoxic fractions and an explicit *"untested there"*.
5. **Standing limitation on the hypoxic band** narrowed from GBM to
   `xenium_319`, with the two counter-values.
6. **"What would most improve this" item (1)** rewritten from an unbought
   improvement to a bought one that failed, naming what is still unbought.
7. **Data table:** added the third GBM donor (`GSM9572843` / `GSE324268` tumor
   07, 851,509 cells), which this page now quotes and previously did not list.
8. **"Where to look next":** a row pointing at the source report.

### 2.2 `landing-3d/lib/copy.ts`

**Rendered copy — three changes, and one of them is a deletion that paid for
the other two.**

* beat 06 body: `— delivery or potency —` **removed**.
* beat 06 ledger, `The 03 verdict`: state `INVERTS TWO WAYS` → `INVERTS THREE
  WAYS`; detail gains *"And so does another section of the same tumour."*
* beat 06 ledger, `“Arrives almost everywhere”`: **unchanged in the end** — see
  below.

**Why the deletion, and the measurement that forced it.** The obvious repair to
`delivery or potency` was to substitute the two settings the model has actually
produced — *"too weak a dose, or a cell its own state protects"*. That was
written and measured. It is 29 characters longer, and together with the two
ledger additions it **clipped beat 06 by 21.6px at 360x700 and 9.3px at
360x740** (`npm run audit:viewports`, exit 1, those two viewports named).
Baseline beat 06 headroom at 360x700 is **13.2px** against a ledger detail line
of ~13.8px — under one line. This file's own rule is that room may not be bought
out of a caveat, and that retracting an overclaim is the one edit that frees
space without costing disclosure. So the overclaim was retracted rather than
restated, and the room went to the verdict row.

**Measured after:** beat 06 worst case **692.8 of 700 at 360x700, 7.2px spare**,
down from 13.2px. The margin is now under half a detail line. The next edit that
needs a line on beat 06 has nothing left to sell.

**What was withdrawn for room, recorded as a known omission rather than
discovered later.** Two clauses, and **both of them cost us**, which is the
direction that matters:

1. On the verdict row: that of the three routes, the first two are modelling
   choices (a term the model lacks, a band threshold we chose) and only the
   third is tissue.
2. On the penetration row: *"So does the specimen — of three glioblastoma
   sections, 03's is the most penetrable."*

Both are recorded in `lib/copy.ts` at the point a future editor meets the rows,
with the instruction that they are the first two things beat 06 should buy. The
page therefore currently reads **softer than the evidence** on both, and that is
stated here rather than left to be found.

**Guardrail comments (no layout cost) added to `lib/copy.ts`:**

* A new head-of-file guardrail, *"THE PAGE'S LEAD IS READ OFF THE ONE SECTION OF
  THREE THAT GIVES THAT ANSWER, AND THAT SECTION IS THE MOST PENETRABLE OF THE
  THREE"*, carrying the numbers, what it licenses, what it forbids, the bounds,
  and — in its own section, so it cannot be dropped — what survived.
* An addition to the delivery-degeneracy guardrail recording that `delivery` is
  dominant on none of the seven sections and never exceeds 0.001817, which
  reinforces the guardrail from a new direction.

### 2.3 `landing-3d/app/layout.tsx`

**No edit.** The `<meta description>` was read in full and checked against this
result: it makes no cross-section claim, contains no ordering, and does not use
"validated". Recorded here as *checked and cleared* rather than silently
skipped, because this is the surface that once outlived a claim's removal by a
month.

### 2.4 Verification of the landing page

Run after the final state of the file:

| command | exit |
|---|---:|
| `npm run build` (includes `audit:origins`) | **0** |
| `npm run lint` | **0** |
| `npm run audit:viewports` | **0** — 11 viewports × 2 motion settings |
| `npm run audit:network` | **0** — zero cross-origin requests, scene fully loaded |

`npm run build` with `NEXT_PUBLIC_STUDIO_URL` unset exits **1** by design
(`lib/env.ts` refuses an unset or loopback value in a production build). That is
a pre-existing guard, not a regression; the builds above set it to a placeholder
public URL.

---

## 3. Recommended edits in files this pass does NOT own

Each carries the exact file, the exact current text and replacement text.
**None of these was applied.**

### 3.1 `app/server/jobs.py` — the sentence that reaches every run document

`verdict_stability()` is written into `extra_metrics` on every run. Its `note`
is compound-scoped prose about a section-scoped measurement, and `app/web`
prints the prose, not the `measured_on` key.

Current (`app/server/jobs.py`, in `verdict_stability`):

```python
        "note": (
            f"{key} KEEPS its delivery/dose verdict when the missing "
            "transvascular exchange term is switched on."
            if stable else
            f"{key} REVERSES its delivery/dose verdict when the missing "
            "transvascular exchange term is switched on, so the verdict above "
            "is one of the unstable ones."
        ),
```

Replacement:

```python
        "note": (
            f"{key} KEEPS its delivery/dose verdict when the missing "
            f"transvascular exchange term is switched on — MEASURED ON "
            f"{VERDICT_STABILITY_SECTION} ONLY, and this run may not be on "
            "that section. The verdict is also section-dependent independently "
            "of this term: at byte-identical thresholds the dominant limiter "
            "reads dose on xenium_319 and resistance on xenium_io and "
            "gse324268 (results/multisection_generalisation_2026-09-20.md "
            "sec. 4.1), so a HOLD here is a statement about the exchange term "
            "on one section, never about tissue."
            if stable else
            f"{key} REVERSES its delivery/dose verdict when the missing "
            f"transvascular exchange term is switched on, so the verdict above "
            f"is one of the unstable ones. MEASURED ON "
            f"{VERDICT_STABILITY_SECTION} ONLY. The verdict is also "
            "section-dependent independently of this term: at byte-identical "
            "thresholds the dominant limiter reads dose on xenium_319 and "
            "resistance on xenium_io and gse324268 "
            "(results/multisection_generalisation_2026-09-20.md sec. 4.1)."
        ),
```

Also recommended, same file, the comment above `VERDICT_STABILITY_SECTION`:

> Current: `#: Measured on the flagship only. Naming the section stops the flag being read` / `#: as a property of the compound anywhere.`
>
> Replacement: `#: Measured on the flagship only. Naming the section stops the flag being read` / `#: as a property of the compound anywhere — and, since 2026-09-20, it must also` / `#: stop it being read as a property of TISSUE: the dominant limiter itself is` / `#: section-dependent at byte-identical thresholds (dose on xenium_319,` / `#: resistance on xenium_io and gse324268,` / `#: results/multisection_generalisation_2026-09-20.md sec. 4.1), so this flag` / `#: scores ONE axis of instability on ONE section.`

### 3.2 `app/server/report.py` — *"This is the answer."*

Current (in the `SURVIVED_UNREACHED` lead):

```python
                '<span class="under">SURVIVED_UNREACHED — per-cell, on this '
                'measured tissue: cells that survived because the drug never '
                'arrived in quantity at their position, not because they '
                'resisted it. This is the answer. '
```

Replacement:

```python
                '<span class="under">SURVIVED_UNREACHED — per-cell, on this '
                'measured tissue: cells that survived because the drug never '
                'arrived in quantity at their position, not because they '
                'resisted it. This is the answer ON THIS SECTION. Resistance '
                'is a verdict this engine returns and it is the dominant one '
                'on two of the three glioblastoma sections measured at '
                'byte-identical thresholds '
                '(results/multisection_generalisation_2026-09-20.md sec. 4.1). '
```

### 3.3 `app/server/reach.py` — `CONVENTION_NOTE`, the understatement

This ships with every reach payload. The promotion of the delivery spread to
headline was justified **on the flagship alone**, and it now has three-donor
support. This is the one place where the correction runs in our favour, and a
document understating the code is read with the same trust as one overstating
it.

Current (tail of `CONVENTION_NOTE`):

```python
    "The headline is instead the delivery spread "
    "(P90/P10 of tumour-cell exposure), which needs no threshold at all and "
    "carries signal even where the map is flat. The full F(x) ladder is still "
    "returned, and is still the evidence — read the curve, not any one point "
    "on it."
```

Replacement:

```python
    "The headline is instead the delivery spread "
    "(P90/P10 of tumour-cell exposure), which needs no threshold at all and "
    "carries signal even where the map is flat. THAT PROMOTION WAS MADE ON THE "
    "FLAGSHIP ALONE AND HAS SINCE BEEN CHECKED ON MORE THAN ONE SPECIMEN: "
    "across three single-cell GBM sections from three donors, one compound, "
    "byte-identical thresholds, the delivery spread reads 1.9936 / 2.3211 / "
    "2.5293 — a spread of 1.2687x, the smallest continuous readout in that "
    "sweep by more than a factor of two — while F(0.1) is 0 on all three and "
    "F(1.0) differs by 119.9x between the two sections where it is computed. "
    "On the one section solved at two grids the delivery spread moves 1.0004x "
    "between dx = 30 and dx = 12 "
    "(results/multisection_generalisation_2026-09-20.md sec. 4.2, 5.7). READ "
    "THAT AT ITS ACTUAL STRENGTH: it is reproducibility of MAGNITUDE on three "
    "specimens, NOT evidence that this statistic discriminates between "
    "tissues, and three specimens is three specimens. The full F(x) ladder is "
    "still returned, and is still the evidence — read the curve, not any one "
    "point on it."
```

### 3.4 `app/server/reach.py` — the docstring that still calls `F(0.1)` the headline

The module docstring and `CONVENTION_NOTE` contradict each other, and the sweep
makes it material.

> Current, in *The quantities (exact definitions)*: `* **The headline** — ``F(0.1)``, picked off that ladder, not recomputed.`
>
> Replacement: `* **The ladder headline point** — ``F(0.1)``, picked off that ladder, not` / `  recomputed. It is NOT the product headline; :data:`CONVENTION_NOTE` says so` / `  and says why. It reads 0 on all three single-cell GBM sections measured at` / `  dx = 30 (results/multisection_generalisation_2026-09-20.md sec. 5.5), which` / `  is the de-headlining decision arriving as a measurement rather than a` / `  preference.`

The constant comment `#: The headline point of the ladder.` above `HEADLINE_X`
should gain `Of the LADDER, not of the product — see CONVENTION_NOTE.`

### 3.5 `app/server/cross_section.py` — `NO_CROSS_SECTION_NOISE_FLOOR`

The refusal to claim a noise floor is correct and should stay. What is now
available is a measured cross-donor span for one statistic, which the note does
not mention.

Append to `NO_CROSS_SECTION_NOISE_FLOOR`:

```python
    " ONE cross-donor span does now exist, for ONE statistic and ONE compound, "
    "and it is not a noise floor either: across three single-cell GBM sections "
    "the delivery spread E_p90/E_p10 reads 1.9936 / 2.3211 / 2.5293 (1.2687x) "
    "while the fitted penetration length spans 2.691x and F(1.0) spans 119.9x "
    "(results/multisection_generalisation_2026-09-20.md sec. 4.2, 4.3). Three "
    "donors is not a distribution and none of those numbers licenses a "
    "significance test here."
```

### 3.6 `app/web/js/ui/finding.js` — the unconditional refutation sentence

This is the most consequential `app/web` defect found: the block is appended
without being gated on the dominant limiter, so on a `resistance` section the
panel asserts the claim above says something it does not.

Current:

```js
      `Measure ${drug} directly in this same tissue — imaging mass `
      + 'spectrometry, or microdialysis on a matched block. The claim above '
      + 'says drug is present essentially everywhere in the tumour and is '
      + 'essentially everywhere below the concentration that produces half its '
      + 'effect. Finding tumour regions with no drug in them would refute the '
      + 'delivery half directly. ',
```

Replacement (gate the two middle sentences on `topKey === 'dose'`, the same key
the correctly-scoped sentence in this file already uses):

```js
      `Measure ${drug} directly in this same tissue — imaging mass `
      + 'spectrometry, or microdialysis on a matched block. ',
      topKey === 'dose'
        ? 'The claim above says drug is present essentially everywhere in this '
          + 'section and is essentially everywhere below the concentration that '
          + 'produces half its effect. Finding tumour regions with no drug in '
          + 'them would refute the delivery half directly. '
        : 'The claim above is about THIS section: a measured map disagreeing '
          + 'with it in this tissue would refute it. It is not a statement '
          + 'about tissue in general — at byte-identical thresholds this '
          + 'engine returns a different dominant limiter on other '
          + 'glioblastoma sections. ',
```

### 3.7 `app/web/js/ui/finding.js` — the verdict-stability sentence

Add the section scope the payload already carries. Current:

```js
             ? h('b', `THIS VERDICT HOLDS FOR ${String(vstab.compound).toUpperCase()} `
                 + '— but it reverses for 7 of the 12 shipped compounds. ')
```

Replacement:

```js
             ? h('b', `THIS VERDICT HOLDS FOR ${String(vstab.compound).toUpperCase()} `
                 + `ON ${String(vstab.measured_on || 'THE FLAGSHIP').toUpperCase()}, `
                 + 'WHERE IT WAS MEASURED — but it reverses for 7 of the 12 '
                 + 'shipped compounds, and the limiter is section-dependent '
                 + 'besides. ')
```

The same `measured_on` qualifier belongs on the `REVERSES` branch.

### 3.8 `app/web/` fixed binary headings

Three fixed strings render identically on a `resistance` section and should name
the third outcome. They are not false about the *question*, only about the
*answer set*; the minimal honest fix is on the one that asserts an output.

* `js/ui/reports.js` — *"One specimen and one compound in, one discrimination
  out: was it delivery, or was it dose?"* → *"One specimen and one compound in,
  one discrimination out: was it delivery, was it dose, or was it the cell's own
  state?"*
* `js/ui/inspect.js` rail heading *"Delivery or dose"* and `js/ui/doseinput.js`
  *"(the argument is under Delivery or dose)"* — these name a **panel**, not an
  output. **Leave them.** Renaming a navigation label to chase an output
  taxonomy would break the cross-references without improving a claim.
* `js/ui/provenance.js` — *"is the delivery-versus-potency split among
  survivors, above — a ratio inside one map"* → *"is the per-cell survival split
  among survivors, above — a ratio inside one map"*. This one renders on
  **every** run under the kill metric.

### 3.9 `docs/architecture.md`

Current:

```
**The hypoxic fraction is grid-sensitive and only order-of-magnitude.** The same
Visium tissue reports 100% hypoxic at dx = 100 µm and 0% at dx = 12 µm; the Xenium
section reported ≈94% coarse and ≈13% at dx = 12 µm. Earlier attribution of small
predicted kills to hypoxic resistance was wrong — they are driven mainly by limited
delivered exposure versus these drugs' IC50s.
```

Replacement:

```
**The hypoxic fraction is grid-sensitive AND donor-sensitive, and is only
order-of-magnitude.** The same Visium tissue reports 100% hypoxic at dx = 100 µm
and 0% at dx = 12 µm; the Xenium section reported ≈94% coarse and ≈13% at
dx = 12 µm. It is not only the mesh: holding dx fixed at 30 µm under the shipped
configuration, three registered single-cell GBM sections report 0.00%, 34.18%
and 18.88% (`results/multisection_generalisation_2026-09-20.md` §4.1). The
earlier attribution of small predicted kills to hypoxic resistance was wrong
**on the flagship section**, where delivered exposure against IC50 is the
binding constraint — but the attribution does not travel: on the other two
sections the dominant limiter is `resistance`, not `dose`, at byte-identical
thresholds.
```

### 3.10 `docs/MATHEMATICAL_MODEL.md`

> Current: `**the predicted hypoxic tumor fraction is $0.0$ at every $\Delta x$**`
>
> Replacement: `**the predicted hypoxic tumor fraction is $0.0$ at every $\Delta x$ ON THIS SECTION** (it is 34.18% and 18.88% on the two other registered single-cell GBM sections at $\Delta x = 30\ \mu m$, `results/multisection_generalisation_2026-09-20.md` §4.1)`

### 3.11 `docs/VERIFICATION_METHODOLOGY.md`

> Current: `under the current shipping default the HYPOXIC state is *empty at every Δx*, so the hypoxia-resistance multiplier has no cells to act on`
>
> Replacement: `under the current shipping default the HYPOXIC state is *empty at every Δx on `xenium_319`*, so the hypoxia-resistance multiplier has no cells to act on THERE. It is not a property of the default: at Δx = 30 µm the same configuration gives 34.18% on `xenium_io` and 18.88% on `gse324268`, and on both of those the dominant limiter is `resistance` (`results/multisection_generalisation_2026-09-20.md` §4.1)`

### 3.12 Two loose uses of "validated", unrelated to this result but found in the sweep

`docs/architecture.md` (*"lets the reconstruction and the pharmacology be
validated independently of each other"*) and `docs/real_data/pk_anchors.md`
(*"transport engine can be validated against reality more broadly"*). Both are
bare, both are the word `docs/outreach/README.md` already says to police, and
the policing was applied there but not here. Recommended: *"scored"* in both
places.

---

## 4. What this audit did NOT settle

1. **Whether the verdict flip survives a change of mesh.** Only `xenium_319` has
   a second grid, and it is not one of the sections that flips. This is the
   single cheapest thing that would change how much weight §1 carries, and no
   copy anywhere should imply it has been done.
2. **Whether the flip is tissue or panel on `xenium_io`.** Its endothelial class
   — the source term of the whole delivery solve — is called by a three-gene
   gate where the other two use one gene, and no section in this checkout
   carries both panels. The flip itself does **not** rest on this, because
   `xenium_319` and `gse324268` share a map and disagree; any three-way spread
   including `xenium_io` does.
3. **Whether `E_p90/E_p10` discriminates between tissues.** Three donors landing
   within 1.2687× is consistent with a well-behaved statistic and equally
   consistent with one that carries little tissue information. Nothing here
   separates those, and no copy may imply the first.
4. **Anything about any compound but temozolomide**, and nothing comparing
   compounds. Cross-compound ranking stays withdrawn.
5. **Anything about spot-resolution tissue.** Four of seven are refused. "Does
   it generalise across platform?" has no answer in this checkout — an absent
   answer, not a negative one.
6. **`results/cross_tissue.md`'s status.** It is cited by
   `app/server/cross_section.py` as "the measured precedent", it is one compound
   on five GBM microenvironments four of which are Visium, and it was computed
   at `alpha = 148.6` against the live 22.2235. It was not rewritten (it is a
   dated record) and it was not re-run. Whether it should still be cited as a
   precedent is left open and named.
7. **Whether any *scored* claim moves.** No anchor in `vivum validate` was
   re-run on the other two sections by this pass or by the sweep. Every
   "n = 1 in tissue" caveat on a scored row stands untouched.

---

## 5. The defect class this belongs to, stated once

This is the same defect the repository already tracks in
`results/claim_surface_audit_2026-08-30.md`: **the engine is sound and the
documents are wrong.** Every number in §0 came out of the shipped `execute_run`
path, every refusal fired as designed, and `mech_reconciliation.bit_exact` is
`True` on 8 of 8 entries. Nothing in the code misbehaved. What went wrong is
that a set of sentences across five surfaces describes a property of
`xenium_319` in the grammar of a property of GBM, of the shipping default, or of
the product.

The tally over the 35 claims in §1: **16 OVERSTATED, 4 UNDERSTATED, 15
SUPPORTED** (one of the supported rows, W4, was checked and cleared rather than
assumed). Two things are worth naming about the direction of the errors:

* **Four of the thirty-five run in our favour and were corrected upward** —
  `CONVENTION_NOTE`'s promotion of the delivery spread, `README.md` carrying no
  cross-donor figure at all, `cross_section.py`'s noise-floor note, and beat
  06's penetration row. That is the side the source report also found the
  2026-08-23 commit message erring on, and it is the side with no guard
  watching: nothing in this repository fails a build because a document is more
  pessimistic than the code.
* **The worst instance is not in a document at all.** It is
  `app/web/js/ui/finding.js`'s unconditional refutation sentence and
  `app/server/jobs.py`'s section-blind `note`, both of which reach a user
  through a run rather than through prose. A stale sentence in a `results/` file
  is read by whoever opens that file; a stale sentence in `extra_metrics` is
  read by everyone who runs anything.

---

## 6. Constraints honoured by this pass

* **Cross-compound ranking**: not restored, not implied, not referenced as a
  limit. No ordering between compounds appears anywhere above.
* **Absolute kill fraction**: quoted nowhere in this file, and nothing in §1 or
  §2 is built on it.
* **No clinical framing**: no response prediction, no diagnostic or
  treatment-planning language, no "digital twin", no personalised medicine, no
  AI/ML claim. "Validated" is used only when naming its misuse.
* **Tissue registers by filesystem path.** The word "upload" appears above only
  in the audit of whether it appears elsewhere. It does not, about tissue.
* **Every number here is traceable** to `results/multisection_generalisation.json`,
  to `results/multisection_generalisation_2026-09-20.md`, to a named source file
  and line, or to a command whose exit code is printed in §2.4. No number was
  invented and no citation was created.
* **No file in `results/` or `archive/` was rewritten.** This file is the only
  one created. Files edited by this pass: `README.md` and
  `landing-3d/lib/copy.ts`. `landing-3d/app/layout.tsx` was read and
  deliberately left unchanged. Nothing under `app/`, `docs/`, `tmesim/`,
  `experiments/`, `render.yaml` or `Dockerfile` was modified.

---

## 7. Closing calibration guard

```
$ .venv/bin/python -c "from tmesim.validation.calibration import calibrate_alpha; a=calibrate_alpha(); assert a==22.223477711172563, a; print('alpha OK',a)"
alpha OK 22.223477711172563
```

Identical to the opening value, and identical to the `alpha_at_run` carried by
every entry of the file this audit is against.
