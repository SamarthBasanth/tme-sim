# Re-anchoring the flagship drug panel to a section that can measure cell state

**Experiment:** `experiments/flagship_markers.py` (run: `./.venv/bin/python experiments/flagship_markers.py`)

`results/cell_state_impact.md` established that the section almost every headline number in this project comes from — `data_cache/xenium_gbm`, the 40,887-cell 319-gene Xenium brain panel — **cannot support a cell state at all**: it carries no hypoxia-response gene of the marker set (HYPOXIC is UNMEASURED there and no cell is assigned it), and MKI67+TOP2A fire in 1.42% of its tumour cells, 13x below the published 19-40% GBM Ki-67 range, so a per-cell proliferation call is refused. Its own conclusion was that absolute kill and kill ceilings must not be quoted from that section.

This report re-runs the **flagship five-drug panel** [temozolomide, lomustine, cisplatin, methotrexate, doxorubicin] on a section that DOES clear the bar, with `state_mode='markers'`, `sensitivity_mode='uniform'`, dx = 12 um — and reports what changes.

## Headline

| | |
|---|---|
| **Primary section** | `data_cache/gse324268_xenium` — 851,509 cells (515,887 tumour), 366-gene panel |
| **Cell state** | MEASURED: proliferation `per_cell` over ['MKI67', 'TOP2A'], hypoxia `per_cell` over ['VEGFA', 'CA9', 'ADM', 'SLC2A1'] — **nothing UNMEASURED** |
| **Measured tumour state fractions** | VIABLE 13.43% / QUIESCENT 33.99% / HYPOXIC 52.58% |
| **Cross-drug ranking (markers)** | cisplatin > doxorubicin > lomustine > temozolomide > methotrexate |
| **Same ranking as the published 319-gene run?** | **NO** — see §4 |
| **Temozolomide, MGMT unmeasured (the old assumption)** | kill 0.1168 |
| **Temozolomide, MGMT MEASURED on this section** | kill 0.0485 (sensitivity 0.4175) |
| **Is absolute kill now defensible?** | see §6 — **no, but the reason has moved** |

> **This is a different patient on a different panel.** Nothing here transfers from or to the 319-gene section as biology, and neither section carries a drug-response outcome label. What travels is the ORDERING and the machinery; the levels do not.

## 1. Which section, and why

| section | panel | cells (tumour) | proliferation | hypoxia | MGMT | verdict |
|---|---|---|---|---|---|---|
| `gse324268_366` | 366 genes | 851,509 (515,887) | `per_cell` (31.15%) | `per_cell` (52.58%) | **YES** (1.72% of cells) | **PRIMARY** |
| `io_480` | 480 genes | 816,769 (540,533) | `per_cell` (34.06%) | `per_cell` (65.07%) | no | check |
| `brain_319` | 319 genes | 40,887 (23,812) | `neighbourhood` (1.42%) | `UNMEASURED` (0.00%) | no | comparison arm — cannot support cell state |

**GSE324268 is chosen as primary, and the deciding argument is MGMT.** Both large sections clear the proliferation detection floor and both measure hypoxia, so on cell state alone the choice would be close — the IO panel actually has the stronger proliferation readout (MKI67 34.06%, inside the published 19-40% Ki-67 range, against GSE324268's 13.34%). But GSE324268 is the **only cached section in this project that carries MGMT**, and MGMT is temozolomide's resistance determinant. On every other section TMZ's sensitivity is left at a uniform 1.0 and its kill is flagged an UPPER BOUND because nobody has measured its resistance gene (`results/real_drug_panel.md` §*Sensitivity mode*). Here it can be measured (§5). A section on which four of five drugs have a measurable determinant and the fifth does not is a section on which the panel is still unfair; this is the first one where that asymmetry can be closed with a number.

Two secondary reasons: GSE324268 measures hypoxia over **four** transcripts (VEGFA, CA9, ADM, SLC2A1) against the IO panel's one (VEGFA), and it is the only section carrying **CA9**, the one near-hypoxia-restricted gene of the set — so its hypoxic band is the only one in the project with any constraint on its SIZE rather than only its ordering. The IO section was run anyway, as a check, and §4 reports whether the two agree.

## 2. The measured cell state, and where every label comes from

Fractions are over TUMOUR cells — the population every kill number is about. The `oxygen` column is the shipping default on the SAME reconstruction, for contrast only.

| section | route | VIABLE | QUIESCENT | HYPOXIC | mesh-dependent? |
|---|---|---:|---:|---:|---|
| `gse324268_366` | `markers` | 13.43% | 33.99% | 52.58% | **no** — computed from measured coordinates and counts in um; no grid enters |
| `gse324268_366` | `oxygen` | 8.52% | 71.93% | 19.55% | yes — the tissue mask gates the vascular source, which sets the band edges |
| `io_480` | `markers` | 11.36% | 23.57% | 65.07% | **no** — computed from measured coordinates and counts in um; no grid enters |
| `io_480` | `oxygen` | 0.00% | 66.23% | 33.76% | yes — the tissue mask gates the vascular source, which sets the band edges |
| `brain_319` | `markers` | 29.80% | 70.20% | 0.00% | **no** — computed from measured coordinates and counts in um; no grid enters |
| `brain_319` | `oxygen` | 10.04% | 89.96% | 0.00% | yes — the tissue mask gates the vascular source, which sets the band edges |

### 2.1 Marker provenance, per section

**`gse324268_366`** — GEO GSE324268 GBMtumor07 — 366-gene brain + GBM add-on panel

* **Proliferation** — call level `per_cell` over ['MKI67', 'TOP2A'], per-gene detection over tumour cells MKI67 13.34%, TOP2A 25.72%, union 31.15%. Missing from this panel: none.
* **Hypoxia** — call level `per_cell` over ['VEGFA', 'CA9', 'ADM', 'SLC2A1'], per-gene detection over tumour cells VEGFA 45.52%, CA9 1.03%, ADM 10.33%, SLC2A1 12.63%, union 52.58%. Missing from this panel: ['NDRG1'].
* **UNMEASURED states:** **none** — every band on this section is traceable to a measured transcript.

**`io_480`** — 10x Xenium Human Brain Cancer (GBM) FFPE v2.0.0 — 480-gene IO panel

* **Proliferation** — call level `per_cell` over ['MKI67'], per-gene detection over tumour cells MKI67 34.06%, union 34.06%. Missing from this panel: ['TOP2A'].
* **Hypoxia** — call level `per_cell` over ['VEGFA'], per-gene detection over tumour cells VEGFA 65.07%, union 65.07%. Missing from this panel: ['CA9', 'ADM', 'NDRG1', 'SLC2A1'].
* **UNMEASURED states:** **none** — every band on this section is traceable to a measured transcript.

**`brain_319`** — 10x Xenium Human Brain — Glioblastoma v1.3.0 — 319-gene brain panel

* **Proliferation** — call level `neighbourhood` over ['MKI67', 'TOP2A'], per-gene detection over tumour cells MKI67 0.55%, TOP2A 0.99%, union 1.42%. Missing from this panel: none.
* **Hypoxia** — call level `UNMEASURED` over NO GENE ON PANEL. Missing from this panel: ['VEGFA', 'CA9', 'ADM', 'NDRG1', 'SLC2A1'].
* **UNMEASURED states:** HYPOXIC.

**What the HYPOXIC band is, exactly.** Positivity for a hypoxia-RESPONSE transcript, not an oxygen tension. VEGFA, ADM and SLC2A1 all have real non-hypoxic expression, so the band's ORDERING is measured and its SIZE is uncalibrated — the 52.6% HYPOXIC band on the primary section is a threshold artefact as a *hypoxic fraction*, not a finding. Only CA9 is near-restricted to severe hypoxia and it fires in 1.03% of that section's tumour cells — a plausible perinecrotic fraction, but too sparse for a per-cell call on its own. Nothing in this report quotes the HYPOXIC fraction's magnitude as a hypoxic fraction; it is used only as the band the engine's `PDConfig.hypoxia_resistance` acts on.

**The one-label collapse, and why it bites hardest here.** Proliferation and hypoxia are two axes and `core.CellState` has one slot; a hypoxia-positive cell is labelled HYPOXIC whatever its proliferation call. That precedence is forced by the schema. On the 319-gene section it cost nothing because the HYPOXIC band was empty; on a section where hypoxia IS measured it puts 52.6% of tumour cells under the x0.35 multiplier and so dominates the state term. **Measuring a band you previously could not see makes the unvalidated multiplier attached to it matter more, not less** — see §6.3.

## 3. The drug panel

Five drugs, `sensitivity_mode='uniform'` (identical for every drug, so it cancels out of the ranking), dx = 12 um, one reconstruction and one transport solve per drug reused across every state route. `therapeutic_coverage_1x` is the fraction of tumour cells whose time-averaged exposure exceeds 1x IC50; it and `mean tumour AUC` carry **no cell-state dependence at all**, so they are identical between the two kill columns by construction.

### 3.1 `gse324268_366` — 851,509 cells, grid 1196x880, alpha = 148.6

| drug | IC50 (uM) | BBB | mean tumour AUC (uM*s) | therapeutic_coverage_1x | **kill (markers)** | kill (oxygen) | ceiling (markers) |
|---|---:|---:|---:|---:|---:|---:|---:|
| temozolomide | 606 | 0.125 | 3.35e+06 | 0.0000 | **0.1168** | 0.0973 | 0.7249 |
| lomustine | 50 | 0.973 | 4.01e+06 | 0.9975 | **0.3885** | 0.4167 | 0.4883 |
| cisplatin | 24.6 | 0.858 | 7.88e+06 | 0.9999 | **0.4649** | 0.4920 | 0.4883 |
| methotrexate | 12 | 0.023 | 6.74e+04 | 0.0003 | **0.0774** | 0.0900 | 0.4883 |
| doxorubicin | 0.5 | 0.021 | 4.92e+04 | 0.9953 | **0.4027** | 0.4310 | 0.4883 |

**Ranking by kill, `markers` route:** **cisplatin** (0.4649) > **doxorubicin** (0.4027) > **lomustine** (0.3885) > **temozolomide** (0.1168) > **methotrexate** (0.0774)

**Ranking by kill, `oxygen` route (same cells, same transport):** **cisplatin** (0.4920) > **doxorubicin** (0.4310) > **lomustine** (0.4167) > **temozolomide** (0.0973) > **methotrexate** (0.0900)

**Ranking by `therapeutic_coverage_1x`:** cisplatin > lomustine > doxorubicin > methotrexate > temozolomide

### 3.2 `io_480` — 816,769 cells, grid 905x956, alpha = 148.6

| drug | IC50 (uM) | BBB | mean tumour AUC (uM*s) | therapeutic_coverage_1x | **kill (markers)** | kill (oxygen) | ceiling (markers) |
|---|---:|---:|---:|---:|---:|---:|---:|
| temozolomide | 606 | 0.125 | 2.87e+06 | 0.0000 | **0.1042** | 0.0797 | 0.7520 |
| lomustine | 50 | 0.973 | 3.46e+06 | 0.9702 | **0.3536** | 0.3486 | 0.4592 |
| cisplatin | 24.6 | 0.858 | 6.79e+06 | 0.9987 | **0.4330** | 0.4248 | 0.4592 |
| methotrexate | 12 | 0.023 | 5.8e+04 | 0.0000 | **0.0649** | 0.0657 | 0.4592 |
| doxorubicin | 0.5 | 0.021 | 4.22e+04 | 0.9691 | **0.3654** | 0.3603 | 0.4592 |

**Ranking by kill, `markers` route:** **cisplatin** (0.4330) > **doxorubicin** (0.3654) > **lomustine** (0.3536) > **temozolomide** (0.1042) > **methotrexate** (0.0649)

**Ranking by kill, `oxygen` route (same cells, same transport):** **cisplatin** (0.4248) > **doxorubicin** (0.3603) > **lomustine** (0.3486) > **temozolomide** (0.0797) > **methotrexate** (0.0657)

**Ranking by `therapeutic_coverage_1x`:** cisplatin > lomustine > doxorubicin > temozolomide > methotrexate

### 3.3 `brain_319` — 40,887 cells, grid 424x675, alpha = 148.6

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

> **But note what it does NOT reproduce, because it matters for reading every older number in `results/`.** `results/real_drug_panel.md` and `results/heldout_patient.md` quote *lower* kills for these same drugs on this same section — temozolomide 0.0147, lomustine 0.1975, cisplatin 0.3699, methotrexate 0.0162, doxorubicin 0.2165 — and this run does not reproduce them. The cause is a pipeline drift, not a harness difference: those reports' own reconstruction block records mean vascular density **0.0250** and mean oxygen **0.2955** on this section, where the current code gives **0.0268** and **0.3520** on identical inputs and an unchanged alpha = 148.61. More drug reaches the tissue, so every kill rises ~1.2x. `results/cell_state_impact.md` and `results/kill_calibration.md` postdate that change and agree with this run exactly. **The two older reports' absolute numbers are stale and should be regenerated; their KILL RANKINGS are unaffected BY THIS DRIFT** (the ordering below is identical to theirs). That is a claim about a ~1.2x common-mode rise in delivered drug, which moves every drug together and therefore cancels out of an ordering. It is not a claim that those reports' rankings are unconditional: a DRUG-SPECIFIC change to the delivery law does move them ([`results/transvascular_exchange_revalidation.md`](transvascular_exchange_revalidation.md) §4.1).

## 4. Does the ranking hold across SECTIONS and across ROUTES at once?

`results/cell_state_impact.md` §3.3 found the cross-drug ranking invariant across three state routes on ONE section. `results/heldout_patient.md` §3.3 found it invariant across two sections on ONE route. This table crosses both axes simultaneously, which is the stronger claim.

| section | route | ranking by predicted kill |
|---|---|---|
| `gse324268_366` | `markers` | cisplatin > doxorubicin > lomustine > temozolomide > methotrexate |
| `gse324268_366` | `oxygen` | cisplatin > doxorubicin > lomustine > temozolomide > methotrexate |
| `io_480` | `markers` | cisplatin > doxorubicin > lomustine > temozolomide > methotrexate |
| `io_480` | `oxygen` | cisplatin > doxorubicin > lomustine > temozolomide > methotrexate |
| `brain_319` | `markers` | cisplatin > doxorubicin > lomustine > methotrexate > temozolomide |
| `brain_319` | `oxygen` | cisplatin > doxorubicin > lomustine > methotrexate > temozolomide |

**The ranking is NOT invariant.** The distinct orderings observed are:

* cisplatin > doxorubicin > lomustine > methotrexate > temozolomide — `brain_319`/markers, `brain_319`/oxygen
* cisplatin > doxorubicin > lomustine > temozolomide > methotrexate — `gse324268_366`/markers, `gse324268_366`/oxygen, `io_480`/markers, `io_480`/oxygen

> ### Resolving this file's own contradiction, in favour of the table above
>
> **This file asserts BOTH "the ranking is NOT invariant" (here, §4) and "invariant
> across sections AND routes" (§6) and "cross-drug ranking | either — it is
> invariant (§4)" (§7). They cannot all stand. The measured table above is the
> authority and the other two are wrong as written** — they are unconditional prose
> in `experiments/flagship_markers.py` (lines 1293 and 1360) sitting above a §4
> paragraph that is emitted from a measured `if invariant:` branch, so the prose
> does not track its own evidence. Fixing this file without fixing that generator
> restores the contradiction on the next run.
>
> **What the six cells actually support**, stated so it cannot be over-read:
>
> * **The top 3 by kill — cisplatin > doxorubicin > lomustine — is invariant across
>   all six (section x route) combinations.** That is the real, and still useful,
>   result of §4.
> * **Ranks 4 and 5 are NOT invariant ACROSS SECTIONS**: `brain_319` gives
>   methotrexate > temozolomide, `gse324268_366` and `io_480` both give
>   temozolomide > methotrexate.
> * **Within any one section the ordering IS invariant across the two state
>   ROUTES**, all five places, which is the `results/cell_state_impact.md` §3.3
>   result reproduced on two further sections.
>
> So the axis matters: **route-invariant, top-3 section-invariant, tail
> section-dependent.** "Invariant across sections and routes" is true only of the
> top 3.
>
> **And none of the three axes covers the delivery law.** Section, route and the
> potency band are all COMMON-MODE — they move every drug together, which is
> exactly the cancellation mechanism the `if invariant:` branch above credits.
> Correcting the one-way vascular source `S = PS*Cp` to the gradient-driven
> `S = PS*(Cp - C)` is drug-specific: it re-weights each drug by its own
> `bbb_penetrance`, which spans **47x** across the GBM panel, and it changes the
> GBM **top-1** (cisplatin -> vincristine) — i.e. it breaks the very top-3 that
> survives everything in this file — changes the top-3 set on both panels, and
> swaps 4 of 45 GBM pairs, two of them across gaps of 40.5 % and 35.6 % of the
> leader's kill against a Bernoulli sd of ~0.003. Name the quantity too: every
> ranking in this file is by **kill**, and on breast the delivery (AUC) ranking is
> preserved exactly (rho = 1.000) while the kill ranking scrambles (rho = 0.830),
> because the Hill map saturates. Every exchange-arm figure is additionally
> conditional on `PS_SCALE_S_INV` (`tmesim/treatment/drugs.py`, default 1.0,
> calibrated 0.16122 1/s). See
> [`results/transvascular_exchange_revalidation.md`](transvascular_exchange_revalidation.md)
> §4.1 and §7.3.

**But the invariance is not equally safe at every adjacent pair, and the narrowest gap should be quoted with it.** The engine draws kills as Bernoulli trials over tumour cells, so a gap of order the sampling scale `sqrt(p(1-p)/n)` is not a resolved ordering:

| section | tightest adjacent pair (markers route) | gap | Bernoulli sigma | gap / sigma |
|---|---|---:|---:|---:|
| `gse324268_366` | doxorubicin vs lomustine | 0.0142 | 0.00068 | 20.8 |
| `io_480` | doxorubicin vs lomustine | 0.0117 | 0.00065 | 17.9 |
| `brain_319` | methotrexate vs temozolomide | 0.0014 | 0.00098 | 1.5 |

Every other adjacent pair is far wider. So the ordering is solid at the top — cisplatin's lead over the next drug is orders of magnitude above the sampling scale on every section — and the **methotrexate / temozolomide boundary at the bottom is the one place the ranking rests on a small margin**. §5 shows that measuring MGMT opens that particular gap decisively, which is the substantive reason to care about it here.

**`therapeutic_coverage_1x` still does not recover the kill ordering — on 3 of 3 sections.** It swaps doxorubicin and lomustine on both large sections (doxorubicin is delivery-starved but 1200x more potent, so it kills more while covering less) and it is pinned at 0.0000 for temozolomide and methotrexate wherever their exposure never reaches 1x IC50. `results/heldout_patient.md` §6.3 found exactly this on two sections and declined to call it a generalisation failure; a third section agreeing makes it a settled property of the metric rather than of any tissue.

**What moves instead is the level.** Same drugs, same physics, different state route:

| section | drug | kill (oxygen) | kill (markers) | ratio |
|---|---|---:|---:|---:|
| `gse324268_366` | temozolomide | 0.0973 | 0.1168 | 1.20x |
| `gse324268_366` | lomustine | 0.4167 | 0.3885 | 0.93x |
| `gse324268_366` | cisplatin | 0.4920 | 0.4649 | 0.94x |
| `gse324268_366` | methotrexate | 0.0900 | 0.0774 | 0.86x |
| `gse324268_366` | doxorubicin | 0.4310 | 0.4027 | 0.93x |
| `io_480` | temozolomide | 0.0797 | 0.1042 | 1.31x |
| `io_480` | lomustine | 0.3486 | 0.3536 | 1.01x |
| `io_480` | cisplatin | 0.4248 | 0.4330 | 1.02x |
| `io_480` | methotrexate | 0.0657 | 0.0649 | 0.99x |
| `io_480` | doxorubicin | 0.3603 | 0.3654 | 1.01x |
| `brain_319` | temozolomide | 0.0184 | 0.0219 | 1.19x |
| `brain_319` | lomustine | 0.2371 | 0.2809 | 1.18x |
| `brain_319` | cisplatin | 0.4362 | 0.5155 | 1.18x |
| `brain_319` | methotrexate | 0.0192 | 0.0233 | 1.21x |
| `brain_319` | doxorubicin | 0.2595 | 0.3073 | 1.18x |

Across the 15 (section x drug) cells above the state route moves the level by 0.86x-1.31x. Across SECTIONS it moves far more: the same drug's marker-route kill spans up to 5.3x (temozolomide, 0.1168 / 0.1042 / 0.0219).

**A finding worth naming: the state route matters LESS on the sections that can actually measure it.** On the two markers-capable sections the oxygen and marker routes give kills within 0.86x-1.02x of each other; on the 319-gene section, where the marker route has to fall back to a neighbourhood rule and the oxygen route's VIABLE band is a mask artefact, they differ by 1.18x-1.21x. The route disagreement on the incumbent section is not a measurement of how much cell state matters — it is a measurement of how badly that panel constrains it.

So a claim of the form *"cisplatin beats doxorubicin beats lomustine on this tissue"* survives every section and every route tested here, while a claim of the form *"temozolomide kills X% of tumour cells"* does not survive a change of section at all.

## 5. Temozolomide with a MEASURED MGMT — the first time in this project

Every prior temozolomide number in this project carries the flag *"UPPER BOUND — no intrinsic resistance is modelled for this drug, its determinant is not measurable on this panel"*. MGMT is absent from both 10x Xenium panels and has **no** panel surrogate (0/396 genes pass FDR, `results/gbm_tmz_exvivo.md` §B2). **GSE324268 carries MGMT.**

### 5.1 What is measured

| quantity | value |
|---|---:|
| MGMT transcripts in the section | 15,191 |
| cells with >= 1 MGMT transcript | 1.72% (1.38% of tumour cells) |
| panel pseudobulk | 96.36 CPM (log2(CPM+1) = 6.605) |

**MGMT is a per-patient scalar here, not a field, and that is a measurement statement rather than a modelling convenience:** it fires in 1.72% of cells at ~1 transcript each, so a per-cell MGMT map would be counting noise. This is exactly the resolution `tmesim/molecular/mgmt.py` was written for.

### 5.2 The unit bridge — and why the naive number is wrong

`mgmt_sensitivity` is calibrated in **whole-transcriptome** `log2(CPM+1)` (GSE232173 bulk RNA-seq, n = 17). A 366-gene targeted panel has a far smaller denominator, so its pseudobulk CPM is not that unit. Feeding the raw panel CPM in gives log2(CPM+1) = 6.605, which is **above the calibration ceiling 4.467** — an extrapolation, and a unit error rather than a biological finding.

The bridge is measured, not assumed: across the 7 cached **TCGA-GBM bulk** RNA-seq sections, the 366 genes of this panel carry 5.03% of all transcriptome counts (range 3.52-6.69%). Rescaling by that fraction converts a panel pseudobulk CPM into a bulk-equivalent one.

| route | MGMT level | TMZ sensitivity | in calibration range? |
|---|---|---:|---|
| **BRIDGED (primary)** | 4.85 CPM bulk-equivalent (log2 2.548) | **0.4175** | **yes** |
| bridged, over the capture range | 3.39-6.45 CPM | 0.3704-0.4809 | yes |
| BINARY `unmethylated` (independent check) | — (a status, not a level) | 0.3923 | n/a |
| NAIVE panel CPM (**rejected**) | 96.36 CPM (log2 6.605) | 0.1039 | **NO — unit error** |

**Two independent checks land on the bridged value.**

1. **The TCGA GBM cohort's own MGMT.** Those same 7 bulk sections have MGMT at a median 4.51 CPM (range 2.48-22.92). The bridged value for this Xenium section, 4.85 CPM, lands essentially on that median — i.e. a GBM's MGMT level comes out looking like a GBM's MGMT level, from a completely separate assay and cohort.
2. **The binary route agrees to 6%.** The continuous bridged sensitivity is 0.4175; the measured median-kill ratio between MGMT-expressing and MGMT-silent cultures — which needs no unit bridge at all — is 0.3923. Two routes through the calibration that share no arithmetic converge.

This section's tumour is therefore **MGMT-EXPRESSING** (bulk-equivalent 4.85 CPM against the calibration's 1.0 CPM expressing threshold), i.e. the resistant group. That is a measurement about this patient, and it is the opposite of the uniform-1.0 assumption every prior TMZ number in this project carries.

### 5.3 What it does to temozolomide

| MGMT input | sensitivity | TMZ kill (markers route) | x vs uniform | rank in the 5-drug panel |
|---|---:|---:|---:|---|
| none (uniform 1.0 — the old assumption) | 1.0000 | 0.1168 | 1.00x | #4 of 5 — **UPPER BOUND** |
| **MEASURED, bridged (primary)** | 0.4175 | 0.0485 | 0.42x | #5 of 5 |
| binary `unmethylated` | 0.3923 | 0.0455 | 0.39x | #5 of 5 |
| naive panel CPM (rejected) | 0.1039 | 0.0124 | 0.11x | #5 of 5 |
| bridged, 95% CI on lambda | 0.1769-0.9851 | 0.0208-0.1150 | — | #5 of 5 to #4 of 5 |

**Read the rank column carefully.** `results/mgmt_scalar.md` found that supplying a measured MGMT returns TMZ from an artificially-inflated #4 to #5 — but that promotion was specific to `determinant` mode, where TMZ's four competitors are scaled down by their own measured determinants and TMZ alone is not. This run is in `uniform` mode, where no drug is scaled, so TMZ starts last on delivery x potency alone and the measurement cannot demote it further. **What MGMT changes here is the MAGNITUDE and the LABEL, not the rank**: a number that was flagged *'no intrinsic resistance modelled'* now carries the resistance value used, its CI, and where it came from.

**It does, however, make the placement itself robust, and that IS the same correction `results/mgmt_scalar.md` describes.** Without MGMT, temozolomide (0.1168) sits below methotrexate (0.0774) by only -0.0394 — -88.2x the Bernoulli sampling scale on 515,887 tumour cells, i.e. a placement that is nearly a coin-flip. With the measured MGMT scalar the gap opens to 0.0289 (65x sigma). That mirrors the finding in `results/mgmt_scalar.md` exactly — there the measurement moved TMZ off an artificially-inflated #4 and opened a 0.5 pp gap to 7-9 pp; here it cannot move the rank because `uniform` mode never inflated it, but it converts a marginal last place into a decisive one.

**And it is not a small correction.** The measured scalar cuts temozolomide's predicted kill by 2.4x (0.1168 -> 0.0485). Every temozolomide kill previously published by this project — `results/real_drug_panel.md`, `results/dose_schedule.md`, `results/kill_calibration.md`, `results/heldout_patient.md` — is an upper bound in this dimension, exactly as `results/kill_calibration.md` limitation 8 says it is. This is the first run in the project where that bound is replaced by a measurement rather than restated as a caveat.

**What this does NOT establish.** The mapping is fitted on n = 17 paired cultures with a 95% CI on lambda spanning 0.006-0.680, i.e. a factor of >100 in slope; the direction is solid, the magnitude is approximate. The cultures are dissociated 2-D monolayers. The panel-to-bulk bridge assumes unbiased relative capture by the Xenium probe set, which is not established — the honest reading of the bridge is that it fixes an order-of-magnitude unit error, not that it calibrates a level. And the TCGA sections supplying the bridge are different patients.

## 6. Is absolute kill defensible on this section?

**No — and the honest version of that answer is more useful than it was, because the terms that block it have changed.**

What is now closed on this section that was open on the 319-gene one:

1. **Cell state is MEASURED, not assumed.** Both bands are traceable to measured transcripts at detection rates that are biology rather than probe sensitivity; nothing is UNMEASURED; and the labelling is mesh-independent by construction (§2). On the incumbent section hypoxia was unmeasured, proliferation was 13x below the published index, and the oxygen route's VIABLE band emptied entirely across the mask x mesh grid.
2. **Temozolomide's resistance determinant is MEASURED** (§5). TMZ was the one drug in the panel whose kill could not be scaled by anything measured; it now can be.
3. **The ranking survives a stronger test — but a narrower one than this line claims.** §4's measured table gives: the top 3 (cisplatin > doxorubicin > lomustine) invariant across all six section x route cells, all five places invariant across ROUTES within a section, and ranks 4-5 NOT invariant across sections. "Invariant across sections AND routes" is true only of the top 3. See the resolution box in §4 — and note that none of these axes is the delivery law, which does change the top-1 ([`results/transvascular_exchange_revalidation.md`](transvascular_exchange_revalidation.md) §4.1).

What remains open, in order of size:

1. **No outcome label. This is the one that actually blocks it.** GSE324268 is an untreated resection with no paired drug response, exactly like the other two sections. Everything in §3 is a prediction with nothing to score it against. Measuring the cell state correctly makes the prediction *defensible in its inputs*; only an outcome makes it defensible as a prediction.
2. **Between-patient potency spread, ~17.4x.** `results/kill_calibration.md` §4 propagated the measured inter-culture IQR of ex-vivo temozolomide viability and got a 17.4x band on absolute kill — against a bound of <= 2.9x from the entire cell-state term (the multipliers are x0.35 and x0.5, so no relabelling can move `F_state` outside [0.35, 1.0]). **That comparison is the crux of this section's verdict: the term this report closes was already the smaller of the two.** Closing it is still worth doing — it is the *epistemic* term, the one a better model can remove, while the 17.4x is biological spread that shrinks only by measuring this patient's own sensitivity — but it does not on its own make a level quotable.
3. **The resistance MULTIPLIERS are still unvalidated priors.** x0.35 hypoxic and x0.5 quiescent are declared in `core/config.py:PDConfig` with no citation in this repo. This report measures WHICH cells are in each band; it does nothing whatever about what each band is worth. On a section where the HYPOXIC band is now large this matters MORE than it did on the 319-gene section, where that band was empty and the multiplier acted on no cells at all.
4. **The hypoxic band's SIZE is uncalibrated.** §2 — the band is an enrichment label built from hypoxia-RESPONSE transcripts with real non-hypoxic expression. Its ordering is measured; its magnitude is not. A measured pO2 (pimonidazole IHC or phosphorescence quenching on a serial section) would close this.
5. **One calibration constant** — alpha, fit once on doxorubicin penetration depth on the 319-gene section and NOT refitted here — sets the perivascular sink for every drug on every section.
6. **Grid convergence is inherited, not re-established.** dx = 12 um was shown converged on the 319-gene section's geometry (`results/numerical_rigor.md` §2). These sections are ~20x denser with larger vessel objects; the convergence argument is carried over.

**The operative sentence.** On this section an absolute kill can for the first time be quoted with every term in it traceable to a measurement or to a named unvalidated prior — measured state, measured MGMT, measured geometry, published potency — which is a real change from the 319-gene section, where the state term was a bookkeeping choice with an unbounded swing. It still should not be quoted as a headline result, because a prediction with no outcome to score it against is not a result no matter how well-sourced its inputs are. **The cross-drug KILL ranking remains the honest primary output, and it is now supported by a stronger invariance test than before — stronger, not unconditional.** The test covers state routes and (for the top 3) sections. It does not cover the vascular delivery law, which changes the GBM top-1 ([`results/transvascular_exchange_revalidation.md`](transvascular_exchange_revalidation.md) §4.1).

## 7. Should this section replace the 319-gene one as the project's reference?

**Yes for anything that depends on cell state; no as a blanket swap; and the two should be reported together rather than one silently replacing the other.**

| claim type | quote it from | why |
|---|---|---|
| absolute kill, kill ceiling, state fractions | **GSE324268** | the only section where both bands are measured and MGMT is available |
| temozolomide anything | **GSE324268** | the only section where TMZ's determinant is measurable |
| cross-drug ranking | either for the TOP 3, which is section- and route-invariant (§4); **section-specific for ranks 4-5**, which differ between `brain_319` and the other two | quote the invariance with its axes attached — routes and (top-3 only) sections. It is NOT invariant across the delivery law, which changes the top-1 ([§4.1](transvascular_exchange_revalidation.md)) |
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
| `gse324268_366` | 851,509 | 1196x880 | 88 s | 33.0 min | 34.6 min |
| `io_480` | 816,769 | 905x956 | 24 s | 15.8 min | 16.3 min |
| `brain_319` | 40,887 | 424x675 | 2 s | 3.6 min | 3.7 min |
| **all** | | | | | **54.7 min** |

**Files added by this work (nothing under `tmesim/` was modified):**

* `experiments/flagship_markers.py` — this analysis
* `results/flagship_markers.md` — this report

**Data provenance.** Three open single-cell spatial GBM sections — 10x Xenium Human Brain Glioblastoma (v1.3.0, 319 targets, CC BY 4.0), 10x Xenium Human Brain Cancer GBM FFPE IO panel (v2.0.0, 480 targets, CC BY 4.0), and GEO GSE324268 GBMtumor07 (366 targets) — plus 7 cached TCGA-GBM bulk RNA-seq sections used ONLY for the panel-capture unit bridge in §5.2. None carries a paired drug-response outcome label. The MGMT mapping is calibrated on GEO GSE232173 / Ntafoulis et al., Br J Cancer 2023 (PMID 37620410).
