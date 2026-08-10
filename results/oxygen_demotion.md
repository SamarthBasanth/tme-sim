# The oxygen demotion was wrong: a component deprecated on an inverted inequality

> **STATUS: WITHDRAWN AND REPLACED BY THIS CORRECTION RECORD.**
>
> This document used to argue for demoting the reaction-diffusion oxygen solve in
> favour of a distance-to-vessel surrogate. **Its central premise was an inverted
> inequality.** The solved field out-predicts the distance baseline; it does not
> lose to it. The demotion should never have happened.
>
> The document is **not deleted**. It is kept, retitled, and rewritten as the
> record of the mistake, because this is the most instructive failure in the
> project: *a component was deprecated on a misread inequality while every
> mechanical check in the same repository disagreed with the prose.*

---

## 1. What was claimed

The opening sentence of the superseded version, verbatim:

> `results/state_validation.md` §1.2 found the reaction-diffusion oxygen solve
> **falsified as value-add** — out-predicted against measured hypoxia markers by a
> one-line distance transform built from its own source cells, and by its own
> undiffused source. This report acts on that.

Its headline table, verbatim:

> | **Best predictor of measured VEGFA** | `vessel_distance_exact` at rho = -0.1647 (p = 0.00050) |
> | **The PDE** | rho = -0.1479 (p = 0.00050), rank 5 of 7 |

Its recommendation, verbatim:

> **The surrogate is better and cheaper on the one section where we can measure**
> — rho -0.1647 vs -0.1479, at 31x less compute at dx = 12 um. **There is no
> accuracy argument for the solve.**

And the source claim it rested on, from `results/state_validation.md` §1.2 as it
then read, verbatim:

> **`dist_endothelial` (-0.1204) beats the PDE (-0.1523) by -21%.**

On that basis a module was written — `tmesim/physics/oxygen_proxy.py` — whose
own docstring stated its purpose as supplying a proxy

> whose only defence is that it out-predicts the solve on the one section where
> we can measure.

## 2. Why it was wrong

**Every predictor in the §1.2 table is sign-aligned to oxygen.**
`state_validation.trivial_baselines` explicitly negates the distance
(`"values": -dist_endo`) and negates the cell count, and documents why:

> Each baseline is returned with the SIGN CONVENTION ALIGNED TO OXYGEN — i.e.
> larger means "more oxygenated" — so its Spearman rho against a hypoxia score is
> directly comparable to the simulated field's, no mental sign-flipping.

So a correct predictor gives a **negative** rho, and the **stronger** predictor is
the one with the **larger |rho|**. The §1.2 table's own header says exactly this.

**|-0.1523| > |-0.1204|.** The solved field is the stronger predictor.

The "-21%" figure in the superseded sentence was itself the tell:
0.1204 / 0.1523 = **0.79**, i.e. the baseline is **79% as strong as** the solved
field. That is a statement that the baseline is *weaker*. It was read as though a
negative percentage meant the PDE came out below the baseline.

### 2.1 Every arithmetic path in the repository always disagreed

This is the part worth internalising. The misreading was never anywhere except in
hand-written English. Four independent mechanical checks, all present in the repo
at the time, all disagreed with it:

| check | what it computes | what it returns |
|---|---|---|
| `state_validation.verdict_test1()` | `beats = abs(sim) > abs(best) * (1.0 + BASELINE_MARGIN)`, `BASELINE_MARGIN = 0.05` | **SUPPORTED**, reason string *"the solved field beats every trivial baseline"* |
| `validation/heldout_patient.py` | `abs(dist_rho) > abs(pde_rho)` → `"CONFIRMED" if ... else "REFUTED"` | **REFUTED** on the demotion claim, *"the solved field now wins"* |
| nested partials (§1.2) | partial Spearman both ways | simulated O2 retains **0.1227**; distance retains only **0.0792** |
| depth-partialled column (§1.2) | log-depth removed | **-0.1619** vs **-0.1185** — same ordering |

`verdict_test1` sits in **§1.6 of the very same generated document** whose §1.2
asserted the opposite. The document contradicted itself, in one file, from one
run, for as long as both sections existed. Nobody read down.

**The transferable lesson.** The defect was not in the measurement, the data, or
the statistics — all four numbers were correct throughout, and remain correct
today. It was that a *directional claim was written in prose* while the
corresponding *decision procedure was computed in code*, and the two were never
required to agree. A number rendered into a sentence by hand is an unchecked
assertion. See §7 for the fix.

## 3. What the arithmetic actually says

Re-derived from primitives at dx = 12 um on the IO section (816,769 cells,
480-gene panel), deliberately **not** by calling the reporting layer, so a bug in
the reporting layer could not hide in this check. Target: the measured hypoxia
score (mean z of depth-normalised VEGFA and HIF1A). All predictors sign-aligned to
oxygen.

| predictor | rho vs hypoxia score | rho vs VEGFA | rho with log-depth removed |
|---|---|---|---|
| **SIMULATED O2 (the PDE)** | **-0.1523** | **-0.2004** | **-0.1619** |
| dist_endothelial | -0.1204 | -0.1237 | -0.1185 |
| vascular_kde (the solver's own undiffused source) | -0.1111 | -0.1017 | -0.1119 |
| local_cell_density (sanity floor) | -0.0192 | -0.0031 | -0.0336 |

* **The PDE beats the distance baseline by 1.26x** on |rho| (0.1523279 / 0.1204202
  = 1.26497). *(`state_validation.md` §1.2 rounds this to "1.27x"; 1.26497 rounds
  to 1.26. The direction is not in doubt; only the second decimal place.)*
* **The diffusion step ADDS information rather than destroying it.** The solve
  takes a source field correlating at 0.1111 and returns one correlating at
  0.1523 — a 1.37x gain on |rho|. The superseded document's claim that "the
  diffusion step DESTROYS predictive information" was the same inversion applied
  to a second pair of numbers.
* **Nested partials agree.** Simulated O2 retains 0.1227 with distance held fixed;
  distance retains only 0.0792 with simulated O2 held fixed (the two are nested at
  rho = +0.3006). Distance-to-vessel is the redundant member of the pair, not the
  PDE.

### 3.1 A second, independent defect: this document's own table was stale

The numbers in §1 of the superseded version (**-0.1647** for the surrogate,
**-0.1479** for the PDE) do not merely have their ordering reversed. **They no
longer reproduce at all.** Two commits landed *after* this document was written
and changed the very inputs it compared:

* `4ad221f` — *"Classifier arity fix: endothelial spread 2.91x -> 1.49x across
  three real sections"*. This changed which cells are `CellType.ENDOTHELIAL`.
* `ceb96c0` — *"Discrete vessels PROMOTED to default"*. This changed the vascular
  source the PDE solves against.

`dist_endothelial` is built from the endothelial cell set; the PDE is solved from
the vascular source. Both moved. Re-running the same comparison under current code
at the same dx, on the same section, against the same target:

| route | rho vs VEGFA (as published here) | rho vs VEGFA (current code) | rho vs composite (current code) |
|---|---|---|---|
| **pde** | -0.1479 | **-0.2004** | **-0.1523** |
| vessel_network | -0.1521 | -0.1470 | -0.1343 |
| vessel_distance_exact | -0.1647 | -0.1237 | -0.1204 |
| dist_endothelial | -0.1647 | -0.1237 | -0.1204 |
| vessel_distance (grid-sampled) | -0.1636 | -0.1220 | -0.1183 |
| vascular_kde | -0.1456 | -0.1017 | -0.1111 |
| local_cell_density | -0.0031 | **-0.0031** | -0.0192 |

**Read the last row.** `local_cell_density` is the only predictor here that
depends on neither the endothelial typing nor the solve — it is a count of cells
within a radius. It is **unchanged to four decimal places**. Every row that
touches endothelial typing or the vascular source moved. That is a clean internal
control confirming the cause of the drift rather than assuming it.

The endothelial cell count is the direct evidence: this document reported
*"Segmentation discards 4,003 of **34,179** endothelial cells as singletons"*.
Current code finds **57,733** endothelial cells on the same section (and
segmentation keeps 44,612 points in 7,754 vessels of >= 2 cells, out of 20,875
objects, at a 20 um link distance).

So on this document's **own** primary target, VEGFA alone, under current code the
PDE leads the best surrogate by **1.62x** (0.2004 / 0.1237), and the nested test
on that target agrees: PDE retains **-0.1725** with distance held fixed, distance
retains **-0.0679** with the PDE held fixed.

Two further claims of the superseded §1 also invert under current code, and are
withdrawn:

* *"The 'more physical' vessel definition is WORSE"* (-0.1521 vs -0.1636). It is
  now **better**: `vessel_network` -0.1470 vs `vessel_distance` -0.1220.
* *"The vessel-distance surrogate wins ... a 1.11x larger effect."* It loses.

The identity the document asserted between `vessel_distance_exact` and
`dist_endothelial` **does** hold, and is confirmed exactly (both -0.12366370788383857):
`exp(-d/L)` and `-d` are both strictly monotone in `d`, and a monotone map cannot
change a Spearman rho.

## 4. The honest limits of the reversal — read this before quoting anything above

**The solved field winning does NOT mean the oxygen solve is validated.** A
correction that becomes an overclaim in the opposite direction is the same error
with its sign flipped. Specifically:

1. **rho ~ 0.15 is WEAK.** It explains a couple of percent of the rank variance in
   the marker. It cannot carry a claim about any individual cell, and it is not
   evidence that the absolute oxygen level is right — nothing here tests that.
   `OxygenConfig.consumption` has still never been fitted to a measured oxygen
   tension.
2. **This is a WITHIN-SECTION comparison, not out-of-sample.** Every number in §3
   comes from one section of one patient. It is a re-derivation on the data the
   claim was made from, not a held-out test.
3. **The composite is carried by VEGFA alone.** HIF1A contributes rho = **+0.0140**
   against the simulated field — the wrong sign and indistinguishable from zero.
   The "VEGFA + HIF1A hypoxia score" is, in practice, VEGFA.
4. **VEGFA is a hypoxia-RESPONSE gene, not an oxygen tension**, with its own
   kinetics and its own non-hypoxic inducers. The comparison is internal: it
   establishes that the solve is not redundant with the geometry it is built from,
   not that the solve is accurate about pO2.
5. **The margin is small and has no confidence interval.** 1.26x on a rho of 0.15,
   with the permutation test measuring each rho against zero rather than the
   predictors against each other.
6. **The circularity is inherited by both routes.** `FLT1` and `PLVAP` are IO-panel
   endothelial markers, so the cells defining "distance to a vessel" are themselves
   typed from the panel. This cannot explain the *gap* between the routes, but
   neither route's absolute rho is free of it.

**The defensible conclusion is: "the PDE is the best of four weak predictors, on
one section, within-sample."** It is not "the oxygen solve is validated". If in
doubt, understate it.

### 4.1 The one criticism of the PDE that SURVIVES

The reversal kills the *accuracy* case against the solve. It leaves the
*calibration* case fully intact, and that case is still open.

The PDE's own effective diffusion length is `sqrt(D / (consumption * rho_tum))`
(`oxygen_proxy.pde_diffusion_length_um`), which at the default `OxygenConfig` and
the measured tumour density is **~310 um on the IO section and ~800 um on the
brain section** — 2-5x the published oxygen-diffusion distance of ~150 um
(Thomlinson & Gray, Br J Cancer 1955;9(4):539-549, PMID 13304213,
doi:10.1038/bjc.1955.55).

*(Those two length scales are the module's own reported figures and were **not**
re-derived in this pass. They depend on the measured tumour areal density, which
the `4ad221f` reclassification may also have moved, so they should be re-measured
before being quoted as current. The argument does not depend on their exact
values — only on the length scale being materially larger than the published one,
which a factor of 2-5 comfortably allows.)*

A rank correlation cannot see this. Spearman rho is invariant under any monotone
rescaling, so the solve can win §3 comfortably and still have its diffusion length
wrong by a factor of 2-5. **Both are true at once.** What the reversal removes is
only the inference that this mis-scaling was *destroying predictive information* —
it demonstrably is not (0.1111 -> 0.1523).

Resolving it requires a measured oxygen tension — phosphorescence-quenching pO2 or
pimonidazole IHC on a serial section — not another rank statistic.

## 5. Verdict on `tmesim/physics/oxygen_proxy.py`

**KEEP, with the rationale corrected, and DEMOTE the recommendation that callers
switch to it.**

Being born of a false premise is not by itself disqualifying — the question is
whether the module has an independent justification. It has three, none of which
depend on the withdrawn claim:

1. **A baseline you keep is how a claim stays falsifiable.** The entire force of
   §1.2's SUPPORTED verdict is that the solve was made to beat a physics-free
   alternative built from its own inputs. That alternative must exist and be
   runnable, or the verdict is unfalsifiable decoration. Deleting the surrogate
   would *weaken* the corrected conclusion, not strengthen it. This is the
   strongest argument for keeping it, and it is the opposite of the argument that
   created it.
2. **Cost — the one claim the reversal leaves untouched.** The surrogate is a
   `cKDTree` build plus one nearest-neighbour query against a sparse direct solve
   on 865,180 unknowns at dx = 12 um. This was never an accuracy argument, so the
   inversion does not touch it. Measured this pass at dx = 12 um: vessel-distance
   grid **1.497 s**, at exact cell coordinates **1.051 s**, vessel-network variant
   **11.876 s**. *(The PDE solve was not timed in isolation in this pass; the
   superseded §2 runtime table — 11x to 31x — is pre-fix and should be re-measured
   before being quoted. The structural argument that a kd-tree query is
   asymptotically cheaper than a sparse solve stands on its own.)*
3. **`hypoxia_from_markers` is a MEASUREMENT and was never contingent on any of
   this.** It reads hypoxia off measured transcripts with no physics, and returns
   `None` with an explicit reason where the panel cannot measure it (the 319-gene
   brain panel carries no VEGFA, CA9, ADM, NDRG1 or SLC2A1, so on that section
   hypoxia is UNMEASURED and the interface says so). It is wired into the Studio's
   validation path. This is the most defensible thing in the module.

A fourth item outlived its motivation: `pde_diffusion_length_um` produced the
calibration finding in §4.1, which is real, independent of the rank comparison,
and still unresolved.

**What must change, and has:**

* The module docstring's title and premise are rewritten. It previously opened
  *"Oxygen SURROGATES that are not a solved oxygen tension — and are better at
  it"* and is now *"Cheap oxygen BASELINES, and the one route that is a
  measurement"*, carrying an explicit withdrawal warning at the top.
* Every in-module claim that the surrogate out-predicts the solve is replaced with
  the measured relation: it scores **0.79x** the solved field's |rho| — cheaper
  **and weaker**.
* The superseded recommendation *"use `oxygen_mode='vessel_distance'` for any NEW
  hypoxia claim"* is **withdrawn**. `'pde'` remains the default, and it is now the
  default on the evidence rather than in spite of it. `'vessel_distance'` is
  demoted to an explicitly-weaker cost-constrained option.

**What is NOT recommended:** deleting the module, or flipping the default. The
default never moved (`oxygen_mode='pde'`, bit-identical, asserted by
`tests/test_oxygen_proxy.py`), and that is the one piece of luck in this episode:
**no flagship result in `results/` was ever computed through the surrogate.** The
only surrogate-driven numbers anywhere were the sensitivity arm in the superseded
§3 of this document, which was explicitly a comparison and never fed anything
downstream. The damage was to the documentation and to the stated rationale of a
module — not to any published prediction.

## 6. What survives from the superseded document

Not everything in it was downstream of the inversion. These findings stand and are
worth keeping:

* **Therapeutic coverage is oxygen-invariant *by construction*.**
  `delivery_metrics.therapeutic_coverage` is a pure function of per-cell AUC and
  IC50, and the AUC field has no oxygen dependence anywhere in the pipeline. **No
  oxygen model of any kind — PDE, surrogate, or a perfect measurement — can change
  a delivery metric in this codebase.** Every coverage number in
  `results/delivery_metrics.md`, `results/dose_schedule.md` and
  `results/vessel_network_impact.md` is immune to this entire question. This is a
  structural argument and does not depend on any rho.
* **The oxygen route is a second-order knob** next to the tissue-mask convention,
  which swings the VIABLE fraction far harder.
* **`oxygen_mode='markers'` raising on a panel that cannot measure hypoxia** is the
  right interface behaviour and should be kept.
* **The `exp(-d/L)` mapping cannot be tuned to win**, because every strictly
  monotone map of distance has the identical Spearman rho. `L` sets band edges and
  nothing else, which is why it is anchored to published morphometry rather than
  fitted.

**What must be re-measured before it is quoted again:** the §2 runtime table and
the entire §3 downstream drug arm (state occupancies, kill deltas of -0.0039 /
-0.0262 / -0.0160) were computed pre-`4ad221f` and pre-`ceb96c0`. They are not
re-derived here. Their *qualitative* conclusion — that the oxygen route barely
moves predicted kill and cannot move coverage at all — is structural and survives;
their digits are stale.

## 7. The generator fix — so this cannot regenerate back to the defect

`results/state_validation.md` is a **generated** document. Correcting the markdown
without correcting `tmesim/validation/state_validation.py` would mean the next
regeneration silently reinstates the inversion — and a corrected document that
reverts on regeneration is *worse* than an uncorrected one, because it looks
reviewed.

The generator emitted the direction as hard-coded English:

```python
A(f"* **`dist_endothelial` ({de:+.4f}) beats the PDE ({so:+.4f}) by "
  f"{100*(abs(de)/abs(so)-1):.0f}%.** ...")
```

That sentence asserts an ordering **unconditionally**, whatever the numbers do.
Meanwhile `verdict_test1` in the same module computed `abs(sim) > abs(best)`.

**The fix is structural, not textual.** Every ordering claim is now derived from
the numbers through a new helper, `state_validation._ordering()`, which returns
the direction, the ratio, and the verb phrases. No sentence in the module states
an ordering the data did not produce. If the ordering ever legitimately flips, the
prose flips with it.

The generator also now **emits the correction record itself** (the `CORRECTION
(§1.2)` block), so the trail of the mistake is reproduced by regeneration rather
than being a hand-added note that the next run deletes.

**Verified, not assumed.** `write_report` was re-run against the independently
re-derived dx = 12 numbers and the output checked mechanically:

* With the real measured numbers: all nine superseded sentences are **absent**
  (*"is out-predicted by a physics-free distance-to-nearest-vessel baseline"*,
  *"beats the PDE"*, *"also beats the PDE"*,
  *"the diffusion step \*destroys\* predictive information"*,
  *"**The PDE is the redundant variable in that pair.**"*, *"The oxygen field is
  not falsified — it is out-competed"*, *"The model has converged onto an answer
  that a distance transform beats"*, *"and it does not matter"*, and the §1.3
  thresholding sentence), and the corrected text plus the `CORRECTION (§1.2)`
  block are **present**.
* With the two rho values artificially **swapped**, the prose flips back to the
  "out-competed" wording and the correction block disappears — confirming the text
  is following the arithmetic in both directions rather than being hard-coded the
  other way.

## 8. Sites that still carry the superseded claim

Every quote below was read from the file itself, not inherited from another
document's summary. Rows marked CORRECTED were fixed in this pass; the rest are
listed so they can be routed, not so they can be forgotten.

| site | what it says now | what it should say |
|---|---|---|
| `tmesim/reconstruct/real_pipeline.py` L106-108, **L234-235**, L821 | L234: `"FALSIFIED AS VALUE-ADD (results/state_validation.md §1.2, results/oxygen_demotion.md): out-predicted against measured "...` — this string travels into every `reconstruction_report()['oxygen']['caveat']` and thence into the Studio. L821: *"retained as the default despite being falsified as value-add"* | the solve out-predicts every trivial baseline (1.26x on \|rho\|) but at weak absolute rho, within-section. `'pde'` is the default **on** the evidence, not despite it |
| `tmesim/reconstruct/cell_state.py` L25 | *"`results/oxygen_demotion.md` §1 — the oxygen PDE is out-predicted as a ..."* | cite the corrected §3 ordering of this document |
| `tmesim/validation/oxygen_benchmark.py` **L707**, **L784** | L707 emits *"**falsified as value-add** — out-predicted against measured hypoxia markers "*. L784 hard-codes *"**The vessel-distance surrogate wins.**"* and then interpolates `winner = ranked[0]`. Its ranking (`sorted(..., key=lambda r: r["rho_vegfa"])`) is in fact **correct** while every rho is negative — most-negative first is largest-\|rho\| first — so under current numbers `ranked[0]` is now the **PDE**, and the sentence would emit *"The vessel-distance surrogate wins. `pde` at rho = -0.2004 against the PDE's -0.2004 — a 1.00x larger effect"* | give the prose the same `_ordering()` treatment applied to `state_validation.py`, and make the ranking robust to a positive rho by sorting on `abs()`. **Highest-priority remaining site**: it regenerates *this* document, so running it as-is overwrites this correction with the defect |
| ~~`docs/architecture.md`~~ | **no action needed.** Searched for `kd-tree`, `value-add`, `out-predict`, `redundant`, `demot` — **no instance of the claim survives in this file** | — |
| ~~`docs/VERIFICATION_METHODOLOGY.md` L36, L39-40, L1075~~ | **CORRECTED IN THIS PASS.** Defect #11 was listed as a successful falsification (*"Oxygen PDE falsified as value-add — beaten by a kd-tree query"*) | now marked **WITHDRAWN**, with the retraction and the reframing of the transferable pattern as *prose/arithmetic divergence* rather than trivial-baseline failure |
| `results/heldout_patient.md` (~L344) | **mostly already correct** — L238 reads *"a one-line `cKDTree` query ... scores -0.1204 — 0.79× the solved field ... so the solve does add predictive information here"*, and L336 says *"loses to the PDE"*. But L344 still concludes *"the hypoxia rows confirm that the reaction-diffusion solve is out-predicted by a distance transform"* | L344 alone contradicts L238 and L336 of its own document, and must be flipped. This file is the closest to correct of any site, because its generator adjudicates with `abs()` |
| `app/reports/rep_00001,3,4/report.html` | *"oxygen (pde): FALSIFIED AS VALUE-ADD ... (rho -0.1503 vs -0.1132)"* | **backwards by its own quoted numbers**: \|-0.1503\| > \|-0.1132\| |
| `tests/test_oxygen_proxy.py:257` | asserts `"FALSIFIED AS VALUE-ADD" in rep["oxygen"]["caveat"]` | pins the defective string in place; must be updated together with `real_pipeline.py` |

**Note the last two rows.** The Studio caveat is internally inconsistent by its
own quoted numbers, and a test actively pins the defective wording so that fixing
the caveat string will fail the suite until both are changed together. That is how
a wrong claim becomes load-bearing.

## 9. Reproduction

The correction in §3 was re-derived from primitives (loader -> `reconstruct_real`
-> `cKDTree` -> `spearman_full`) without calling the reporting layer, at
dx = 12 um on `data_cache/xenium_gbm_io`, 816,769 cells. The §3.1 refresh
additionally exercised `oxygen_proxy.oxygen_from_vessel_distance` and
`vessel_distance_um` under current code.

`tmesim/validation/oxygen_benchmark.py` has **not** been re-run to regenerate this
document, because its prose generator still hard-codes the superseded ordering
(§8). Running it as-is would overwrite this correction with the defect. It needs
the `_ordering()` treatment first.

**Data provenance.** 10x Genomics Xenium Human Brain Cancer (GBM) FFPE IO-panel
section, v2.0.0, 480 targets. CC BY 4.0, no login or DUA.
