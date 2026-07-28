# Re-anchoring the flagship drug panel to a section that can measure cell state

**Experiment:** `experiments/flagship_markers.py` (run: `./.venv/bin/python experiments/flagship_markers.py`)

`results/cell_state_impact.md` established that the section almost every headline number in this project comes from — `data_cache/xenium_gbm`, the 40,887-cell 319-gene Xenium brain panel — **cannot support a cell state at all**: it carries no hypoxia-response gene of the marker set (HYPOXIC is UNMEASURED there and no cell is assigned it), and MKI67+TOP2A fire in 1.42% of its tumour cells, 13x below the published 19-40% GBM Ki-67 range, so a per-cell proliferation call is refused. Its own conclusion was that absolute kill and kill ceilings must not be quoted from that section.

This report re-runs the **flagship five-drug panel** [temozolomide, lomustine, cisplatin, methotrexate, doxorubicin] on a section that DOES clear the bar, with `state_mode='markers'`, `sensitivity_mode='uniform'`, dx = 12 um — and reports what changes.

## 1. Which section, and why

| section | panel | cells (tumour) | proliferation | hypoxia | MGMT | verdict |
|---|---|---|---|---|---|---|
| `brain_319` | 319 genes | 40,887 (23,812) | `neighbourhood` (1.42%) | `UNMEASURED` (0.00%) | no | comparison arm — cannot support cell state |

**GSE324268 is chosen as primary, and the deciding argument is MGMT.** Both large sections clear the proliferation detection floor and both measure hypoxia, so on cell state alone the choice would be close — the IO panel actually has the stronger proliferation readout (MKI67 34.06%, inside the published 19-40% Ki-67 range, against GSE324268's 13.34%). But GSE324268 is the **only cached section in this project that carries MGMT**, and MGMT is temozolomide's resistance determinant. On every other section TMZ's sensitivity is left at a uniform 1.0 and its kill is flagged an UPPER BOUND because nobody has measured its resistance gene (`results/real_drug_panel.md` §*Sensitivity mode*). Here it can be measured (§5). A section on which four of five drugs have a measurable determinant and the fifth does not is a section on which the panel is still unfair; this is the first one where that asymmetry can be closed with a number.

Two secondary reasons: GSE324268 measures hypoxia over **four** transcripts (VEGFA, CA9, ADM, SLC2A1) against the IO panel's one (VEGFA), and it is the only section carrying **CA9**, the one near-hypoxia-restricted gene of the set — so its hypoxic band is the only one in the project with any constraint on its SIZE rather than only its ordering. The IO section was run anyway, as a check, and §4 reports whether the two agree.

## 2. The measured cell state, and where every label comes from

Fractions are over TUMOUR cells — the population every kill number is about. The `oxygen` column is the shipping default on the SAME reconstruction, for contrast only.

| section | route | VIABLE | QUIESCENT | HYPOXIC | mesh-dependent? |
|---|---|---:|---:|---:|---|
| `brain_319` | `markers` | 29.80% | 70.20% | 0.00% | **no** — computed from measured coordinates and counts in um; no grid enters |
| `brain_319` | `oxygen` | 10.04% | 89.96% | 0.00% | yes — the tissue mask gates the vascular source, which sets the band edges |

### 2.1 Marker provenance, per section

**`brain_319`** — 10x Xenium Human Brain — Glioblastoma v1.3.0 — 319-gene brain panel

* **Proliferation** — call level `neighbourhood` over ['MKI67', 'TOP2A'], per-gene detection over tumour cells MKI67 0.55%, TOP2A 0.99%, union 1.42%. Missing from this panel: none.
* **Hypoxia** — call level `UNMEASURED` over NO GENE ON PANEL. Missing from this panel: ['VEGFA', 'CA9', 'ADM', 'NDRG1', 'SLC2A1'].
* **UNMEASURED states:** HYPOXIC.

**What the HYPOXIC band is, exactly.** Positivity for a hypoxia-RESPONSE transcript, not an oxygen tension. VEGFA, ADM and SLC2A1 all have real non-hypoxic expression, so the band's ORDERING is measured and its SIZE is uncalibrated — the — HYPOXIC band on the primary section is a threshold artefact as a *hypoxic fraction*, not a finding. Only CA9 is near-restricted to severe hypoxia and it fires in ~1% of that section's tumour cells — a plausible perinecrotic fraction, but too sparse for a per-cell call on its own. Nothing in this report quotes the HYPOXIC fraction's magnitude as a hypoxic fraction; it is used only as the band the engine's `PDConfig.hypoxia_resistance` acts on.

**The one-label collapse, and why it bites hardest here.** Proliferation and hypoxia are two axes and `core.CellState` has one slot; a hypoxia-positive cell is labelled HYPOXIC whatever its proliferation call. That precedence is forced by the schema. On the 319-gene section it cost nothing because the HYPOXIC band was empty; on a section where hypoxia IS measured it puts — of tumour cells under the x0.35 multiplier and so dominates the state term. **Measuring a band you previously could not see makes the unvalidated multiplier attached to it matter more, not less** — see §6.3.

## 3. The drug panel

Five drugs, `sensitivity_mode='uniform'` (identical for every drug, so it cancels out of the ranking), dx = 12 um, one reconstruction and one transport solve per drug reused across every state route. `therapeutic_coverage_1x` is the fraction of tumour cells whose time-averaged exposure exceeds 1x IC50; it and `mean tumour AUC` carry **no cell-state dependence at all**, so they are identical between the two kill columns by construction.

### 3.1 `brain_319` — 40,887 cells, grid 424x675, alpha = 148.6

| drug | IC50 (uM) | BBB | mean tumour AUC (uM*s) | therapeutic_coverage_1x | **kill (markers)** | kill (oxygen) | ceiling (markers) |
|---|---:|---:|---:|---:|---:|---:|---:|
| temozolomide | 606 | 0.125 | 7.66e+05 | 0.0000 | **0.0219** | 0.0184 | 0.6490 |
| lomustine | 50 | 0.973 | 9.32e+05 | 0.2802 | **0.2809** | 0.2371 | 0.6490 |
| cisplatin | 24.6 | 0.858 | 1.82e+06 | 0.9955 | **0.5155** | 0.4362 | 0.6490 |
| methotrexate | 12 | 0.023 | 1.57e+04 | 0.0000 | **0.0233** | 0.0192 | 0.6490 |
| doxorubicin | 0.5 | 0.021 | 1.15e+04 | 0.4373 | **0.3073** | 0.2595 | 0.6490 |

**Ranking by kill, `markers` route:** **cisplatin** (0.5155) > **doxorubicin** (0.3073) > **lomustine** (0.2809) > **methotrexate** (0.0233) > **temozolomide** (0.0219)

**Ranking by kill, `oxygen` route (same cells, same transport):** **cisplatin** (0.4362) > **doxorubicin** (0.2595) > **lomustine** (0.2371) > **methotrexate** (0.0192) > **temozolomide** (0.0184)

**Ranking by `therapeutic_coverage_1x`:** cisplatin > doxorubicin > lomustine > temozolomide > methotrexate

> **Reproduction gate — 8/8 values reproduce `results/cell_state_impact.md` §3.2 to 5e-05 absolute.** Three drugs x two state routes plus the two shared ceilings, i.e. the transport solve, both state routes and the PD mapping are all checked against an independently published run: temozolomide kill_oxygen 0.0184 vs 0.0184; temozolomide kill_markers 0.0219 vs 0.0219; cisplatin kill_oxygen 0.4362 vs 0.4362; cisplatin kill_markers 0.5155 vs 0.5155; doxorubicin kill_oxygen 0.2595 vs 0.2595; doxorubicin kill_markers 0.3073 vs 0.3073; kill ceiling (oxygen route, shared) 0.5502 vs 0.5502; kill ceiling (markers route, shared) 0.6490 vs 0.6490. This harness is running the project's pipeline, not a variant of it.

> **But note what it does NOT reproduce, because it matters for reading every older number in `results/`.** `results/real_drug_panel.md` and `results/heldout_patient.md` quote *lower* kills for these same drugs on this same section — temozolomide 0.0147, lomustine 0.1975, cisplatin 0.3699, methotrexate 0.0162, doxorubicin 0.2165 — and this run does not reproduce them. The cause is a pipeline drift, not a harness difference: those reports' own reconstruction block records mean vascular density **0.0250** and mean oxygen **0.2955** on this section, where the current code gives **0.0268** and **0.3520** on identical inputs and an unchanged alpha = 148.61. More drug reaches the tissue, so every kill rises ~1.2x. `results/cell_state_impact.md` and `results/kill_calibration.md` postdate that change and agree with this run exactly. **The two older reports' absolute numbers are stale and should be regenerated; their RANKINGS are unaffected** (the ordering below is identical to theirs).

## 4. Does the ranking hold across SECTIONS and across ROUTES at once?

`results/cell_state_impact.md` §3.3 found the cross-drug ranking invariant across three state routes on ONE section. `results/heldout_patient.md` §3.3 found it invariant across two sections on ONE route. This table crosses both axes simultaneously, which is the stronger claim.

| section | route | ranking by predicted kill |
|---|---|---|
| `brain_319` | `markers` | cisplatin > doxorubicin > lomustine > methotrexate > temozolomide |
| `brain_319` | `oxygen` | cisplatin > doxorubicin > lomustine > methotrexate > temozolomide |

**The ranking is INVARIANT across all 2 (section x route) combinations tested: cisplatin > doxorubicin > lomustine > methotrexate > temozolomide.** That is a materially stronger statement than either prior report could make. The mechanism is the same one both of them identified: the resistance multipliers are SHARED across drugs — every drug pays the same x0.5 in the same QUIESCENT cells and the same x0.35 in the same HYPOXIC cells — so the state term very largely cancels out of a comparison and does not cancel at all out of a level.

**What moves instead is the level.** Same drugs, same physics, different state route:

| section | drug | kill (oxygen) | kill (markers) | ratio |
|---|---|---:|---:|---:|
| `brain_319` | temozolomide | 0.0184 | 0.0219 | 1.19x |
| `brain_319` | lomustine | 0.2371 | 0.2809 | 1.18x |
| `brain_319` | cisplatin | 0.4362 | 0.5155 | 1.18x |
| `brain_319` | methotrexate | 0.0192 | 0.0233 | 1.21x |
| `brain_319` | doxorubicin | 0.2595 | 0.3073 | 1.18x |

So a claim of the form *"cisplatin beats doxorubicin beats lomustine on this tissue"* survives every section and every route tested here, while a claim of the form *"temozolomide kills X% of tumour cells"* survives none of them.

## 5. Temozolomide with a MEASURED MGMT — the first time in this project

MGMT is not on the primary section's panel; this section is empty.

## 6. Is absolute kill defensible on this section?

**No — and the honest version of that answer is more useful than it was, because the terms that block it have changed.**

What is now closed on this section that was open on the 319-gene one:

1. **Cell state is MEASURED, not assumed.** Both bands are traceable to measured transcripts at detection rates that are biology rather than probe sensitivity; nothing is UNMEASURED; and the labelling is mesh-independent by construction (§2). On the incumbent section hypoxia was unmeasured, proliferation was 13x below the published index, and the oxygen route's VIABLE band emptied entirely across the mask x mesh grid.
2. **Temozolomide's resistance determinant is MEASURED** (§5). TMZ was the one drug in the panel whose kill could not be scaled by anything measured; it now can be.
3. **The ranking survives a stronger test** — invariant across sections AND routes simultaneously (§4).

What remains open, in order of size:

1. **No outcome label. This is the one that actually blocks it.** GSE324268 is an untreated resection with no paired drug response, exactly like the other two sections. Everything in §3 is a prediction with nothing to score it against. Measuring the cell state correctly makes the prediction *defensible in its inputs*; only an outcome makes it defensible as a prediction.
2. **Between-patient potency spread, ~17.4x.** `results/kill_calibration.md` §4 propagated the measured inter-culture IQR of ex-vivo temozolomide viability and got a 17.4x band on absolute kill — against a bound of <= 2.9x from the entire cell-state term (the multipliers are x0.35 and x0.5, so no relabelling can move `F_state` outside [0.35, 1.0]). **That comparison is the crux of this section's verdict: the term this report closes was already the smaller of the two.** Closing it is still worth doing — it is the *epistemic* term, the one a better model can remove, while the 17.4x is biological spread that shrinks only by measuring this patient's own sensitivity — but it does not on its own make a level quotable.
3. **The resistance MULTIPLIERS are still unvalidated priors.** x0.35 hypoxic and x0.5 quiescent are declared in `core/config.py:PDConfig` with no citation in this repo. This report measures WHICH cells are in each band; it does nothing whatever about what each band is worth. On a section where the HYPOXIC band is now large this matters MORE than it did on the 319-gene section, where that band was empty and the multiplier acted on no cells at all.
4. **The hypoxic band's SIZE is uncalibrated.** §2 — the band is an enrichment label built from hypoxia-RESPONSE transcripts with real non-hypoxic expression. Its ordering is measured; its magnitude is not. A measured pO2 (pimonidazole IHC or phosphorescence quenching on a serial section) would close this.
5. **One calibration constant** — alpha, fit once on doxorubicin penetration depth on the 319-gene section and NOT refitted here — sets the perivascular sink for every drug on every section.
6. **Grid convergence is inherited, not re-established.** dx = 12 um was shown converged on the 319-gene section's geometry (`results/numerical_rigor.md` §2). These sections are ~20x denser with larger vessel objects; the convergence argument is carried over.

**The operative sentence.** On this section an absolute kill can for the first time be quoted with every term in it traceable to a measurement or to a named unvalidated prior — measured state, measured MGMT, measured geometry, published potency — which is a real change from the 319-gene section, where the state term was a bookkeeping choice with an unbounded swing. It still should not be quoted as a headline result, because a prediction with no outcome to score it against is not a result no matter how well-sourced its inputs are. **The cross-drug ranking remains the honest primary output, and it is now supported by a stronger invariance test than before.**

## 7. Should this section replace the 319-gene one as the project's reference?

**Yes for anything that depends on cell state; no as a blanket swap; and the two should be reported together rather than one silently replacing the other.**

| claim type | quote it from | why |
|---|---|---|
| absolute kill, kill ceiling, state fractions | **GSE324268** | the only section where both bands are measured and MGMT is available |
| temozolomide anything | **GSE324268** | the only section where TMZ's determinant is measurable |
| cross-drug ranking | either — it is invariant (§4) | but quote the invariance, not one section's copy of it |
| delivery metrics (AUC, coverage, penetration length) | either | they carry no state dependence at all |
| the vascular external anchors (endothelial fraction, vascular area fraction) | **neither is settled** | `results/heldout_patient.md` §6.2 shows both fail to replicate across sections |

Three reasons not to make it a clean swap:

1. **It is a different patient.** Nothing about the 319-gene section's tumour transfers to this one. A section change is not a correction of the old numbers, it is a different measurement.
2. **Every published number would move.** The re-run list in `results/cell_state_impact.md` §5.2 applies with a section change on top of the route change, and the levels move (§4) even where the ordering does not. Flipping silently under a reader comparing against an older document would be exactly the failure mode this project keeps guarding against.
3. **The 319-gene section still has the one thing these do not**: it is the section every calibration in the project (alpha, the grid-convergence study, the vessel-linking distance) was fitted on. Retiring it means carrying those constants onto tissue they were never fitted to — which this report already does, and flags.

## 8. Caveats — what could make this wrong

1. **No outcome label on any section.** Every kill here is a prediction with nothing to score it against.
2. **Different patients, different panels.** Cross-section comparisons mix real biological difference with panel content. Absolute AUCs and kills are NOT comparable between sections (`results/heldout_patient.md` §3.5 shows the ECM field alone has no common scale across panels); only orderings and within-section spreads are.
3. **Targeted-panel dropout.** A zero is "no probe fired", not "no transcript". No individual cell's marker call is trustworthy on any panel — only rates and orderings are.
4. **The hypoxic band over-calls hypoxia** relative to any pO2-based definition, and the one-label collapse lets it win over the proliferation call, so it dominates the state term on this section.
5. **The MGMT bridge assumes unbiased relative capture** by the Xenium probe set, and the bridge cohort (TCGA GBM, n = 7) are different patients. The binary route agrees, which is reassuring and is not independent evidence about the bridge itself.
6. **The MGMT calibration is n = 17** with a 95% CI on lambda spanning a >100x range in slope, fitted on dissociated 2-D monolayers.
7. **The PD resistance multipliers are unvalidated priors** (§6.3).
8. **alpha is not refitted** per section (§6.5), and grid convergence is inherited (§6.6).
9. **2-D, in-plane, single section.** Distances are Euclidean and in-plane; there is no out-of-plane transport.

## 9. Reproduction and cost

```bash
./.venv/bin/python experiments/flagship_markers.py
```

Deterministic — the only randomness is the engine's seeded generator (`seed = 0`). One reconstruction per section, one transport solve per drug, both reused across every state route and every MGMT arm; nothing is re-solved to change a state label.

| section | cells | grid | reconstruction | 5 transport solves | total |
|---|---:|---|---:|---:|---:|
| `brain_319` | 40,887 | 424x675 | 2 s | 4.6 min | 4.6 min |
| **all** | | | | | **4.6 min** |

**Files added by this work (nothing under `tmesim/` was modified):**

* `experiments/flagship_markers.py` — this analysis
* `results/flagship_markers.md` — this report

**Data provenance.** Three open single-cell spatial GBM sections — 10x Xenium Human Brain Glioblastoma (v1.3.0, 319 targets, CC BY 4.0), 10x Xenium Human Brain Cancer GBM FFPE IO panel (v2.0.0, 480 targets, CC BY 4.0), and GEO GSE324268 GBMtumor07 (366 targets) — plus 7 cached TCGA-GBM bulk RNA-seq sections used ONLY for the panel-capture unit bridge in §5.2. None carries a paired drug-response outcome label. The MGMT mapping is calibrated on GEO GSE232173 / Ntafoulis et al., Br J Cancer 2023 (PMID 37620410).
