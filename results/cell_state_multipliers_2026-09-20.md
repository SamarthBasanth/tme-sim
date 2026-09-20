# The cell-state multipliers themselves

`results/kill_calibration.md` §5 says what this report is for, in the project's own words:
*"What the cell-state work must deliver is not primarily the fraction but a defence of the
multipliers themselves (x0.8 / x0.5), which are unvalidated priors in `PDConfig` and are NOT
bounded by anything above."*

This is that pass. It is **not** a fourth measurement of the QUIESCENT fraction —
`results/cell_state_impact.md` measured `F_state` three ways on the flagship section
(0.5502–0.9550) and this extends that rather than repeating it.

**Nothing was tuned, adopted or changed.** No multiplier value in `tmesim/` moved. Every number
below is a shipped value, a swept probe, or closed-form arithmetic on values already published
in `results/`. Every threshold was pre-registered before the sweeps ran, and the one that failed
is reported FAILED (§7.1).

## Headline

| | |
|---|---|
| Cell-state multipliers the kill path applies | **3** (plus VIABLE `x1` and DEAD ⇒ `0`, which are structural, not priors) |
| Defended by a citation | **0 of 3** |
| Defended by a documented modelling judgement anchored to real GBM data | **1 of 3** — temozolomide's `x0.80`; its honest form is the interval `[0.675, 1.000]` |
| **Undefended — unvalidated priors, no citation in this repo** | **2 of 3** — `quiescent_resistance = 0.50` and the global `hypoxia_resistance = 0.35` |
| Datasets in `data_cache/` that could supply a value for any of them | **0** — none pairs a drug response with a cell-state stratification at matched exposure |
| Leverage on the flagship section: `quiescent_resistance` over `[0,1]` | **10.1–11.8x** |
| Leverage on the flagship section: `hypoxia_resistance` over `[0,1]` | **exactly 1.000x** — the kill is bit-identical; the HYPOXIC band is empty |
| The same two on `io_480` (34.2 % hypoxic) and `gse324268_366` (18.9 %) | quiescent **4.1–6.6x**, hypoxia **1.3–2.0x** |
| Is the published **2.00x** relabelling bound right? | **Yes for temozolomide** — exactly, tightly, on all three sections. **No as a general statement:** cisplatin and doxorubicin are **2.857x**. |
| Does that bound constrain the multipliers? | **No. It restates them.** At `quiescent_resistance = 0.25` it is 4.00x; at 0.05, 20x. |
| Guards added | 4 in `tests/test_cell_state.py`; 3 plants, each confirmed reached, each red, each reverted (§7.2) |

> **The one-line finding.** The engine treats its two undefended state multipliers as
> interchangeable floats sitting side by side in one dataclass. They are not interchangeable. On
> the section every published kill number in this project was computed on, `quiescent_resistance`
> — which has **no measurement behind it anywhere, for any drug, in any tissue** — moves absolute
> kill by an order of magnitude, while `hypoxia_resistance`, the one the project has spent its
> entire literature-calibration effort on, moves it by **exactly nothing**, because the band it
> multiplies is empty. The effort went to the prior with zero leverage.

> **Read with the fraction work, not instead of it.** `results/cell_state_impact.md` asks which
> cells carry the label, and shows the three routes disagree about *which* while agreeing about
> *how many*. This asks what the label is worth once assigned. The two answers multiply.

## 1. The inventory — every cell-state multiplier the engine applies

Cell state reaches a kill decision through exactly **one** channel. That is asserted, not
assumed: `tests/test_cell_state.py::test_the_kill_path_applies_exactly_three_state_multipliers_and_no_fourth`
puts one cell of each state at identical exposure and identical sensitivity and requires the
per-cell kill probability to be the Hill term times exactly `1.0 / R_hypoxia / R_quiescence / 0`.
Independently, the transport field itself was re-solved with every multiplier set to 1.0 and
differenced against the shipped solve: `max |AUC_default - AUC_all-ones| = 0.0` exactly, on
every drug on every section (pre-registered H1). And the oxygen sink is proportional to local
tumour-cell *density* (`consumption * rho_tum * O`), not to cell state, so there is no
state → oxygen feedback and no further state multiplier upstream of the kill.

| # | multiplier | defined at | value | kind | what defends it today |
|---|---|---|---:|---|---|
| **M1** | `PDConfig.quiescent_resistance` | `core/config.py:239` | **0.50** | prior | **Nothing.** No citation in the repo; the project's own literature pass found no quantitative measurement of any kind, for any of its three drugs, in any tissue. |
| **M2** | `PDConfig.hypoxia_resistance` (global default; in force for cisplatin and doxorubicin) | `core/config.py:238` | **0.35** | prior | **Nothing.** No citation. Real GBM-specific cisplatin evidence exists but its magnitude was never extracted; the value was deliberately left unchanged rather than invented. |
| **M3** | `DrugParams.hypoxia_resistance` = `TMZ_HYPOXIA_RESISTANCE` (temozolomide only) | `core/config.py:317` | **0.80** | modelling judgement | A documented, reproducible conversion of a real GBM-cell-line IC50-fold range to this project's own operating point. Not a citation — the input fold `f` is a choice inside a literature range. |
| **M4** | VIABLE, implicit | `cells/pkpd.py` (`np.ones`) | 1.00 | structural | It *is* the reference state. Nothing to defend. |
| **M5** | DEAD ⇒ `p_kill = 0` | `cells/pkpd.py` (`alive` mask) | 0 | structural | "A killed cell cannot be killed again." Nothing to defend. |
| **M6** | QUIESCENT band top = `2.0 x hypoxia_threshold` | `cells/agents.py` | 2.0 | band edge | Its own docstring: *"a deliberately simple, documented choice"*. Undefended — but it selects **which** cells, i.e. the fraction, which is `cell_state_impact.md`'s subject, not this report's. |
| **M7** | `resistant_modulation = max(pd.hypoxia_resistance, pd.quiescent_resistance)` | `sim/mechanism.py:240` | 0.50 | derived | Derived, and derived **wrong for temozolomide** — §6.1. |

`report/attribution.py` adds no multiplier: its `o2resist` term is a counterfactual *simulation*
with `PDConfig(hypoxia_resistance=1, quiescent_resistance=1)`, i.e. both penalties switched off
together — which is exactly why it is named for the whole oxygen-driven resistance state rather
than for hypoxia alone. That is the honest construction and it is not double-counting anything.

### 1.1 M3's derivation, reproduced

`PDConfig`'s docstring converts an IC50 fold-shift `f` into an equivalent kill-probability
multiplier at fixed exposure via `R(r, f, h) = (1 + r^h) / (f^h + r^h)`, evaluated at
`r = E/IC50 = 0.056` (the median tumour-cell exposure ratio in `results/dose_schedule.md`) and
`h = 1.2`. Recomputed here, at that operating point:

| `f` (hypoxic IC50 fold-shift) | source | `R` at `r = 0.056`, `h = 1.2` |
|---|---|---:|
| **1.0** (no significant effect) | Musah-Eroje & Watson 2019, PMID 30694423 — the matched-oxygen **2-D** comparison, the geometry structurally closest to what `pkpd.py` computes | **1.000** |
| 1.1 | Wang *et al.* 2017, PMID 28801626 (U87) | 0.895 |
| **1.2** | Wang *et al.* 2017 (mid-range) — **the value shipped** | **0.808** |
| 1.3 | Wang *et al.* 2017 (U87 hypoxia-pretreated) | 0.736 |
| 1.4 | Wang *et al.* 2017 (GL261) | 0.675 |

The docstring's own claim — `R` = 0.806–0.811 across this project's documented dose band
`r = 0.034–0.078` — reproduces exactly (0.8062, 0.8083, 0.8106), and the override really is
applied on the real GDSC-calibrated compile path and not only to the synthetic library entry
(`validation/calibration.py:203`). **What the table shows is that the literature-consistent
interval for M3 is `R ∈ [0.675, 1.000]` — a 1.48x span that includes "no effect at all" — and
0.80 is a point chosen inside it.** M3's *direction and order of magnitude* are
literature-anchored for temozolomide; its *value* is a modelling judgement. That is a materially
better position than M1 or M2. It is still not a measurement.

## 2. What would defend each, and what `data_cache/` can bear on it

All three multipliers are the same quantity: **a ratio of kill probability between two cell
states at matched exposure.** Nothing else substitutes directly — that is the "metric mismatch"
the 2026-07-29 literature pass established, and it is why an IC50 fold cannot be plugged in
without first fixing an operating point.

Every dataset in `data_cache/` was checked against that requirement.

| dataset | what it is | can it bear on M1 / M2 / M3? |
|---|---|---|
| `gdsc_glioma_dose_response.csv` | GDSC2, 2-D glioma monolayer, `LN_IC50` / `AUC` per cell line | **No.** One normoxic, cycling condition per line; no oxygen arm, no cell-cycle arm. |
| `drugresponse/` (HGCC panel, `JOINED_GSE232173_TMZviability.csv`) | the ex-vivo 2-D TMZ viability anchor already used by `kill_calibration.md` §4 | **No.** Bulk viability per culture; no state stratification. |
| `hunt/geo_organoid_screens/GSE148842` | **single-cell** RNA-seq of acute human GBM slice cultures, 43 samples, **paired vehicle vs drug** arms | **The closest structural match in the cache, and still no.** It has the paired treated/untreated contrast on real human GBM tissue with per-cell proliferation transcripts. Two gaps: the drug arms are etoposide, panobinostat, RO4929097, tazemetostat, ipsenisib and Ana-12 — **none of this project's three drugs** — and scRNA-seq of survivors measures a composition shift, not a per-cell kill, because dying cells leave the assay. It could bear on the *direction* for a non-library cytotoxic. It cannot supply a value. |
| `hunt/radiopharm_pet/` (ACRIN-FMISO 6684) | FMISO PET hypoxia imaging paired with **clinical outcome** | **No**, twice over: it pairs hypoxia with outcome rather than with a matched-exposure kill, and using it that way would be clinical response prediction, which this project does not do. |
| `hunt/pharmacoscopy_exvivo/` (DARTT-1) | ex-vivo image-based drug screen, **AML**, markers CD3/CD34/CD33 | **No.** Wrong disease, wrong markers, no proliferation or hypoxia stratification. |
| `ivy_gap/` | anatomic-structure RNA-seq (pseudopalisading / perinecrotic zones) | **No.** Bears on *which* cells are hypoxic — the fraction — not on the multiplier. |
| the three Xenium sections | untreated resections | **No.** No drug was given; `kill_calibration.md` §4 says so in terms. |

**Stated plainly: no dataset in `data_cache/` can supply a value for any of the three
multipliers, because none of them pairs a drug response with a cell-state stratification at
matched exposure in glioma.** This is not a gap that better analysis of what is already
downloaded can close.

**What would close it, per multiplier.**

* **M1 (`quiescent_resistance`)** — the minimum viable experiment is small and specific: one GBM
  culture or slice, separated into cycling and G0 populations (EdU / Ki-67 sort, or a FUCCI
  reporter read *in situ*), exposed to **one** temozolomide concentration, viability read on both
  arms. That yields a survival ratio at a known baseline kill probability, which the
  metric-mismatch analysis identifies as the one convertible literature quantity. Until it
  exists, M1 has no magnitude.
* **M2 (`hypoxia_resistance`, cisplatin / doxorubicin)** — one supplementary-table lookup. Ahmed
  *et al.* 2018 (PMID 29492900, PMC5951876 Electronic Supplementary Material) already contains a
  statistically significant, GBM-specific cisplatin IC50 shift under hypoxia in U251 and SNB19;
  the 2026-07-29 pass identified it and did not retrieve it. Running that number through the same
  `R(r, f, h)` conversion M3 already uses would give cisplatin the treatment temozolomide has
  had. **This is the cheapest available improvement in this whole area and it is still
  outstanding** — though see §3 for why it is not the one that matters most.
* **M3 (temozolomide)** — already at the limit of what the literature supports. What remains is
  not more searching but an honest restatement: quote `[0.675, 1.000]`, not the point 0.80,
  wherever the value is load-bearing.

## 3. Leverage, per multiplier

The engine treats M1 and M2/M3 identically — two floats, side by side, carrying the same
"unvalidated prior" label. They do not deserve the same amount of worry, and this section
measures by how much they differ.

**Method.** One `sim.engine.simulate_response` per (section, drug) supplies the transport; the
labels are then held at the engine's own oxygen banding and **one multiplier at a time** is swept
over its full `[0, 1]` range, re-running the unmodified `cells.pkpd.apply_drug_kill` on a private
copy of the cells with the engine's own seeded generator. The swap path is gated: it must
reproduce `simulate_response`'s own `kill_fraction` **bit-for-bit** before any swept number is
believed. It did, for every drug on every section.

**Kill is exactly affine in each multiplier.** Pre-registered threshold: affine to `1e-9`
relative. Measured: maximum relative departure **3.8e-16** across all drugs and sections —
machine precision. The Hill clip is therefore inactive at every exposure this project reaches,
and the leverage of a state multiplier is *exactly* that band's share of the Hill mass:

```
    d(kill)   / d(m_state) = (1/N) * sum of p_hill_i over cells in that band
    d(ceiling)/ d(m_state) = that band's fraction of tumour cells
```

which also makes every ratio in this section **invariant to `alpha`, to the dose and to the
exposure field**. That matters here, because `calibrate_alpha()` now returns
22.223477711172563 where `cell_state_impact.md` ran at 148.608 (§6.3): the absolute kills below
are not comparable to that report's, and the leverage ratios are.

### 3.1 Leverage, three sections, dx = 12 um (flagship) and dx = 30 um (the two ~800k-cell sections)

| section | tumour cells | VIABLE / HYPOXIC / QUIESCENT | drug | **M1 span over `[0,1]`** | **M2/M3 span over `[0,1]`** |
|---|---:|---|---|---:|---:|
| `brain_319` dx 12 | 23,812 | 8.43 % / **0.00 %** / 91.57 % | temozolomide | **10.13x** | **1.000x** |
| `brain_319` dx 12 | | | cisplatin | **11.76x** | **1.000x** |
| `brain_319` dx 12 | | | doxorubicin | **11.45x** | **1.000x** |
| `io_480` dx 30 | 540,533 | 0.00 % / **34.18 %** / 65.82 % | temozolomide | 4.14x | 1.80x |
| `io_480` dx 30 | | | cisplatin | 6.57x | 2.03x |
| `gse324268_366` dx 30 | 515,887 | 8.46 % / **18.88 %** / 72.66 % | temozolomide | 4.22x | 1.32x |
| `gse324268_366` dx 30 | | | cisplatin | 5.82x | 1.42x |

Flagship absolute kills, for the record (alpha = 22.2235, repo-nominal dose):

| drug | multiplier swept | kill at `m = 0` | kill at `m = 1` | span | absolute move |
|---|---|---:|---:|---:|---:|
| temozolomide | `quiescent_resistance` | 0.019780 | 0.200319 | **10.13x** | 0.1805 |
| temozolomide | `hypoxia_resistance` (its `x0.80` override) | 0.112212 | 0.112212 | **1.000x** | **0.000000** |
| cisplatin | `quiescent_resistance` | 0.082395 | 0.969049 | **11.76x** | 0.8867 |
| cisplatin | `hypoxia_resistance` | 0.525743 | 0.525743 | **1.000x** | **0.000000** |
| doxorubicin | `quiescent_resistance` | 0.077020 | 0.881740 | **11.45x** | 0.8047 |
| doxorubicin | `hypoxia_resistance` | 0.476902 | 0.476902 | **1.000x** | **0.000000** |

**The hypoxia multiplier has exactly zero leverage on the flagship section — bit-identical kill
across its entire range, to the last digit, for all three drugs.** That is not a small effect
rounded away: the HYPOXIC band is empty under the shipping mask, so the multiplier multiplies
nothing. Pre-registered threshold `|kill(0) - kill(1)| < 1e-12`; measured **0.0 exactly**.

**And its leverage is a property of the patient, not of the engine.** It is 1.000x on the
flagship, 1.32–1.42x on `gse324268_366` and 1.80–2.03x on `io_480`, tracking that section's
HYPOXIC occupancy and nothing else. There is no configuration in which it overtakes M1.

Over a *literature-plausible* range rather than the full `[0, 1]` — for M1 the literature
supports only "`< 1`", so `[0.05, 1.00]` is as narrow as honesty allows:

| drug (flagship) | M1 over `[0.05, 1.00]` | M2/M3 over `[0.35, 1.00]` |
|---|---:|---:|
| temozolomide | **6.87x** | 1.000x |
| cisplatin | **7.63x** | 1.000x |
| doxorubicin | **7.54x** | 1.000x |

### 3.2 This inverts `kill_calibration.md` §5's ordering

That section compares the cell-state term against the potency anchor and concludes cell state is
"non-dominant": **2.00x** from relabelling, against **17.4x** from the measured inter-patient
potency spread and **2.8x** from the narrow cohort-median bar.

But 2.00x is the movement available from **relabelling**, with the multipliers held at their
shipped values. The movement available from the **multiplier** — the thing §5 itself names as
undefended — is **6.87x** for temozolomide on the same section. That is **2.5x larger than the
narrow potency bar it was being compared against**, and it sits between the two potency figures
rather than below both.

`kill_calibration.md` §5 is right that the *label* is non-dominant. It does not follow, and that
report does not claim, that *cell state* is non-dominant — and on the multiplier axis it is not.
§5's own note 1 (*"the bound above is conditional on them... this is the dependency to track"*)
is exactly the right instinct; this is the quantity it was asking for.

## 4. The boundedness claim, checked rather than inherited

`results/kill_calibration.md` §5: *"Because the engine applies cell state as a bounded multiplier
(x0.8 hypoxic for temozolomide specifically -- its own per-drug override, `core.config.DRUG_LIBRARY`
-- x0.5 quiescent, x1 viable), relabelling cells can move `F_state` only between 0.5 and 1.0. Over
the full 0-100 % quiescent sweep above, temozolomide's absolute 3-D kill moves **2.00x**."*

### 4.1 The general law, and why it did not actually need three sections

Because kill is **exactly** affine in each multiplier (§3, 3.8e-16 — the Hill clip is inactive at
every exposure this project reaches), the maximum movement available from *any* relabelling is

```
    max kill / min kill  =  1 / min(multipliers reachable by that drug)
```

exactly, and **tightly** — both endpoints are attained by labelling every cell. The `p_hill`
distribution cancels out of the ratio, so the bound does not depend on the section, the tissue,
the exposure field, the dose or `alpha`. It is arithmetic, not an empirical property, and the
right way to check it across sections is to confirm the affinity premise. It was measured on
three sections anyway, because the brief asked for more than the one it was derived on.

### 4.2 Measured, per drug, per section

Labels swept 0 → 100 % (exposure-ordered, the pessimistic assignment `kill_calibration.md` uses);
multipliers held at their shipped values.

| section | drug | label swept | **expectation ratio** | analytic `1/floor` | single-seed engine ratio |
|---|---|---|---:|---:|---:|
| `brain_319` dx 12 | temozolomide | QUIESCENT | **2.000000** | 2.000000 | 1.9541 |
| `brain_319` dx 12 | temozolomide | HYPOXIC | 1.250000 | (1/0.80) | 1.2406 |
| `brain_319` dx 12 | cisplatin | QUIESCENT | 2.000000 | 2.000000 | 2.0011 |
| `brain_319` dx 12 | cisplatin | **HYPOXIC** | **2.857143** | **2.857143** | 2.9098 |
| `brain_319` dx 12 | doxorubicin | QUIESCENT | 2.000000 | 2.000000 | 2.0103 |
| `brain_319` dx 12 | doxorubicin | **HYPOXIC** | **2.857143** | **2.857143** | 2.9377 |
| `io_480` dx 30 | temozolomide | QUIESCENT | **2.000000** | 2.000000 | 2.0031 |
| `io_480` dx 30 | temozolomide | HYPOXIC | 1.250000 | (1/0.80) | 1.2512 |
| `io_480` dx 30 | cisplatin | QUIESCENT | 2.000000 | 2.000000 | 2.0004 |
| `io_480` dx 30 | cisplatin | **HYPOXIC** | **2.857143** | **2.857143** | 2.8652 |
| `gse324268_366` dx 30 | temozolomide | QUIESCENT | **2.000000** | 2.000000 | 2.0033 |
| `gse324268_366` dx 30 | temozolomide | HYPOXIC | 1.250000 | (1/0.80) | 1.2510 |
| `gse324268_366` dx 30 | cisplatin | QUIESCENT | 2.000000 | 2.000000 | 1.9999 |
| `gse324268_366` dx 30 | cisplatin | **HYPOXIC** | **2.857143** | **2.857143** | 2.8692 |

**Verdict.**

1. **The 2.00x figure is SOUND, TIGHT and REPRODUCIBLE for temozolomide**, on all three sections,
   exactly. `kill_calibration.md` §5 is correct as scoped — it names temozolomide and it names the
   quiescent label, and both qualifications are load-bearing.
2. **It is not a property of the engine, and a reader who generalises "bounded by construction"
   will get the other two drugs wrong.** Cisplatin and doxorubicin have no per-drug override, so
   their reachable floor is the global `x0.35` and their relabelling bound is **2.857x — 43 %
   larger than 2.00x**. Both drugs sit in the same report's §3 PASS/FAIL panel. That figure
   appears nowhere in `results/`.
3. **The bound restates the unvalidated prior; it does not constrain it.** It is robust to M3:
   across the whole literature-consistent interval for temozolomide's hypoxia multiplier
   (`[0.675, 1.000]`, §1.1) both endpoints exceed 0.50, so the floor stays `quiescent_resistance`
   and the bound stays 2.00x. It is **not** robust to M1, the multiplier with nothing behind it:
   at `quiescent_resistance = 0.25` the bound is 4.00x, at 0.05 it is 20x. So "2.00x" is a
   restatement of a number nobody has measured, and `kill_calibration.md` §5's note 1 is exactly
   right to call it the dependency to track.

### 4.3 The bound is on the expectation. A single engine run breaks it about half the time.

The pre-registered soundness threshold (`measured ratio <= 1/floor x (1 + 1e-6)`) **FAILED on the
single-seed engine path** for cisplatin (2.9098 vs 2.8571, +1.8 %) and doxorubicin (2.9377,
+2.8 %). That is reported FAILED and the threshold was not moved. The cause was then measured
rather than asserted: both endpoints were re-evaluated over **64 seeds** on the same transport
result (flagship section, 23,812 tumour cells).

| drug | analytic bound | expectation ratio | seed 0 | 64-seed mean ± sd | range | seeds above the bound |
|---|---:|---:|---|---|---|---:|
| temozolomide | 2.000000 | **2.000000000** | 1.954117 (`z = +0.70` at the resistant endpoint) | **2.000012** ± 0.030602 | 1.9281 – 2.0898 | **35 / 64** |
| cisplatin | 2.857143 | **2.857142857** | 2.909836 (`z = −2.06`) | **2.854812** ± 0.025501 | 2.7965 – 2.9165 | **30 / 64** |
| doxorubicin | 2.857143 | **2.857142857** | 2.937736 (`z = −2.57`) | **2.856152** ± 0.027469 | 2.7967 – 2.9377 | **28 / 64** |

The excess is Bernoulli sampling, not a broken bound: the expectation lands on `1/floor` to nine
decimals, each 64-seed mean is within one standard error of it, and roughly half the seeds land
on each side — which is what a ratio whose expectation *equals* the bound must do.

**The practical consequence is worth stating, because it is the kind of thing that gets quoted
wrong.** `1/min(multiplier)` bounds the *expected* kill ratio. It does not bound a single seeded
`simulate_response` run, which exceeds it about half the time by up to ~3 % at this tumour-cell
count. A report quoting a relabelling bound should quote the expectation, or quote a single
realisation with its sampling scale attached.

## 5. What is undefended

Stated in the words the brief asked for, without softening.

**`PDConfig.quiescent_resistance = 0.50` is an unvalidated prior.** It has no citation in this
repo and no measurement behind it anywhere. The project's own literature pass
(`results/sweep/state_multiplier_literature_check_2026-07-29T05-54-00Z.md`) went looking with
four separate search strategies and reported the result plainly: *"no study was found that
reports a quantitative resistance ratio or fold-change between quiescent (G0-arrested) and
cycling cancer cells for temozolomide, cisplatin, or doxorubicin, in any tissue."* What exists is
a well-established mechanistic argument for the **direction** (`< 1`, because cell-cycle-dependent
cytotoxics preferentially kill cycling cells) and nothing at all about the **magnitude**. 0.50 is
not a measured value, not a fitted value, and not derived from one. And it is the multiplier that
carries **91.6 %** of tumour cells on the section every published kill number in this project was
computed on, moving absolute kill by **6.9–7.6x** across the range the literature actually
permits.

**`PDConfig.hypoxia_resistance = 0.35` is an unvalidated prior.** It has no citation in this repo.
It was left at 0.35 deliberately and honestly — the 2026-07-29 pass found real GBM-specific
cisplatin evidence whose magnitude it did not extract, and a real doxorubicin fold-change in the
wrong tissue under the wrong exposure protocol, and declined to invent a number from either. That
decision was correct. It does not make the value defended. It is in force for cisplatin and
doxorubicin — the two drugs whose evidence exists and has not been pulled.

**`TMZ_HYPOXIA_RESISTANCE = 0.80` is neither an unvalidated prior nor a measurement.** It is a
**modelling judgement with a documented, reproducible derivation** from real, GBM-cell-line,
hypoxia-versus-normoxia temozolomide IC50 data, converted Hill-consistently to this project's own
operating point. The derivation reproduces exactly (§1.1). What it is not is a lookup: the input
fold-shift `f = 1.2` is one point inside a literature range running from a clean null result
(`f = 1.0`, in the matched-oxygen 2-D comparison structurally closest to what `pkpd.py` computes)
to `f = 1.4`. **The honest form of M3 is the interval `R ∈ [0.675, 1.000]`, not the point 0.80.**

**Vocabulary that must not be used for any of the three.** None is "calibrated", "validated",
"fitted", or "literature-derived" — the last is true only of M3's *conversion procedure*, never of
the value that comes out of it, because that value depends on a choice of `f`. M1 and M2 are
unvalidated priors and should be called that every time they are quoted.

### 5.1 What this report does not do

It proposes no new value and adopts none. **No multiplier in `tmesim/` was changed.** Nothing here
says any of the three multipliers is *wrong* — only that two of them are undefended, and that
their leverage is wildly unequal. There is still no drug-response outcome label on any of the
three sections, so nothing in this repo can currently score a multiplier against reality.

## 6. Three claim-surface defects found on the way

None is fixed here: `tmesim/sim/`, `tmesim/validation/` and existing `results/` files are all
outside this agent's write territory, and `results/` files are dated records.

### 6.1 `sim/mechanism.py` resolves the hypoxia multiplier the way its own contract forbids

`cells/pkpd.py`'s docstring states the rule plainly: *"Any module that recomputes the
hypoxia/quiescence resistance vector directly instead of going through `apply_drug_kill` MUST call
[`effective_hypoxia_resistance`] instead of reading `pd.hypoxia_resistance` directly, or it will
silently ignore the per-drug override."*

`sim/mechanism.py:240-241` reads it directly:

```python
res = max(float(pd.hypoxia_resistance),
          float(pd.quiescent_resistance))
```

`resistant_modulation` is documented as *"the mildest resistance penalty the kill model itself
applies"*. For temozolomide the kill model applies `x0.80` to hypoxic cells, so the mildest
penalty is **0.80**; this derives **0.50**. The kill is unaffected — the kill path resolves
correctly — but the *explanation* is wrong in two places: a temozolomide hypoxic cell carrying
exactly the model's mildest penalty is **not** called `SURVIVED_RESISTANT`, and
`stochastic_min_p_kill`, derived from the same number, moves with it. The run's own provenance
string then prints `hypoxia_resistance=0.35` for a run that applied 0.80 (`mechanism.py:545`).

The defect is **latent on the flagship section** (empty HYPOXIC band → no cell is currently
mislabelled) and fires on any section with a populated HYPOXIC band — 34.18 % on `io_480`,
18.88 % on `gse324268_366` (§3.1). The effective value is already sitting in the diagnostics dict
`apply_drug_kill` hands it, under `diagnostics["hypoxia_resistance"]`.

### 6.2 `results/cell_state_impact.md` §3.1 contradicts its own §3.2, and the arithmetic proves it

§3.1's "resistance multiplier applied" column reads `x1.00 / x0.50 / x0.35` for **all three
drugs**. §3.2's ceiling column, on the same page, does not agree. Using §3.1's own
`vessel_distance` band fractions (VIABLE 0.9108 / QUIESCENT 0.0867 / HYPOXIC 0.0026):

```
  ceiling with r_hyp = 0.35 :  0.9108 + 0.0867*0.5 + 0.0026*0.35 = 0.95506
  ceiling with r_hyp = 0.80 :  0.9108 + 0.0867*0.5 + 0.0026*0.80 = 0.95623
  published                 :  cisplatin 0.9550, doxorubicin 0.9550, temozolomide 0.9561
```

Temozolomide's published ceiling could only have been produced with `x0.80`. The generator prints
the multiplier column from `PDConfig()` globals (`validation/cell_state_impact.py:1125-1126`)
while computing the ceiling drug-resolved (`_ceiling_decomposition(..., drug=dp)`), so the table
disagrees with itself by construction. It goes unnoticed on the flagship's default oxygen route
only because the HYPOXIC band is empty there, and the whole discrepancy hides in the 0.26 % of
cells the `vessel_distance` route labels hypoxic.

### 6.3 Two numbers those reports quote as current have moved

Recorded, not chased — both concern the *fraction* and the *dose*, not the multipliers.

* **Band occupancy.** `cell_state_impact.md` §2.1/§3.1 and `kill_calibration.md` §5 both quote the
  flagship section at dx = 12 um, `geometric` as **VIABLE 10.04 % / QUIESCENT 89.96 % / HYPOXIC
  0.00 %**. A fresh run today on the same path gives **8.43 % / 91.57 % / 0.00 %**. The bands come
  from the engine's own oxygen solve, which has no `alpha` dependence, so this is a change
  upstream of the kill, not a dose effect.
* **`alpha`.** `cell_state_impact.md` §3 quotes `alpha = 148.61`, as do `fstate_mechanism.md` §7
  and the `DRUG_LIBRARY` provenance note. `calibrate_alpha()` returns **22.223477711172563**
  today, asserted at the start and the end of every run behind this report (§7.4).
  **Absolute kills in this report are therefore not comparable to those reports' absolute kills.**
  Every leverage ratio and every bound in §3 and §4 is `alpha`-invariant — a consequence of the
  exact affinity result, not an assumption — so no conclusion here depends on which value is
  right.

## 7. Pre-registration, guards, reproduction

### 7.1 Every threshold, as registered and as scored

Registered before any sweep ran; none was moved afterwards.

| # | claim under test | pre-registered threshold | measured | verdict |
|---|---|---|---|---|
| H1 | the transport field carries no multiplier dependence | `max abs(AUC_default − AUC_all-ones) == 0.0` exactly | **0.0** (all drugs, all sections) | **PASS** |
| H2a | kill is affine in each multiplier (no Hill clipping) | affine to `1e-9` relative | **3.8e-16** worst case | **PASS** |
| H2b | an empty band gives its multiplier exactly zero leverage | `abs(kill(m=0) − kill(m=1)) < 1e-12` on `brain_319` | **0.0** exactly, all three drugs | **PASS** |
| H3a | relabelling cannot exceed `1 / min(reachable multiplier)` — **expectation path** | `<= 1/floor x (1 + 1e-6)` | equals `1/floor` to nine decimals (2.000000000, 1.250000, 2.857142857) | **PASS** |
| H3a' | the same, on the **single-seed engine path** | the same threshold | cisplatin **2.9098** and doxorubicin **2.9377** against a 2.8571 ceiling | **FAIL as registered** — diagnosed in §4.3 over 64 seeds; the threshold was NOT changed |
| H3b | the bound is tight (the endpoints reach it) | `>= 0.99 x` ceiling | reached exactly on the expectation path | **PASS** |
| H3c | the published **2.00x** generalises across drugs | holds for every drug, or the exceptions are named | **NOT GENERAL** — cisplatin and doxorubicin are 2.857x | **FAIL, named** |
| H4 | nothing is tuned | no multiplier value in `tmesim/` changes | `tmesim/` untouched (sha256 verified); the only tracked file changed is `tests/test_cell_state.py` | **PASS** |

### 7.2 The guards, and the plant test

Before this pass **nothing in the test suite pinned the multipliers, their count, their per-drug
resolution, or the boundedness statement derived from them** — three published numbers could have
been invalidated by a one-character edit with every check staying green. Four tests were added to
`tests/test_cell_state.py` (0.05 s, no real section needed):

1. `test_the_kill_path_applies_exactly_three_state_multipliers_and_no_fourth`
2. `test_the_shipped_multiplier_values_are_pinned_to_the_published_claim`
3. `test_hypoxia_multiplier_resolution_is_per_drug_not_global`
4. `test_relabelling_bound_is_one_over_the_smallest_reachable_multiplier`

Each was planted against. Every plant was confirmed to sit inside the executed code path *before*
the run, and reverted after, with `tmesim/cells/pkpd.py`'s sha256 checked back to its pristine
value each time. Plant 3 changes a value in `core/config.py`, which is outside this agent's write
territory, so it was applied as an out-of-tree pytest plugin — identical effect on every
`PDConfig()`, no repo file touched.

```
--- baseline: pristine tree, four guards green ---
c66ca83db6bd82c2edde1c2f079ebdeab4dc494a188ba7a78950b4b8cc4c6fa3  tmesim/cells/pkpd.py
....                                                                     [100%]
4 passed, 19 deselected in 0.05s

--- PLANT 1: a FOURTH state multiplier (VIABLE x0.95) in cells/pkpd.py ---
PLANT REACHED: the line is inside apply_drug_kill, executed by every kill call
E           AssertionError: VIABLE: kill probability is not p_hill x 1.0
E           assert np.float64(0.475) == 0.5 ± 1.0e-12
E             Obtained: 0.475
E             Expected: 0.5 ± 1.0e-12
E           assert np.float64(0.4) == 0.38 ± 1.0e-12
E             Obtained: 0.4
E             Expected: 0.38 ± 1.0e-12
E           assert 1.9 == 2.0 ± 2.0e-12
FAILED tests/test_cell_state.py::test_the_kill_path_applies_exactly_three_state_multipliers_and_no_fourth
FAILED tests/test_cell_state.py::test_hypoxia_multiplier_resolution_is_per_drug_not_global
FAILED tests/test_cell_state.py::test_relabelling_bound_is_one_over_the_smallest_reachable_multiplier
3 failed, 1 passed, 19 deselected in 0.22s
reverted -> c66ca83db6bd82c2edde1c2f079ebdeab4dc494a188ba7a78950b4b8cc4c6fa3

--- PLANT 2: read pd.hypoxia_resistance directly, ignoring the per-drug override ---
PLANT REACHED: the line is inside apply_drug_kill, executed by every kill call
E           assert 0.35 == 0.8
FAILED tests/test_cell_state.py::test_hypoxia_multiplier_resolution_is_per_drug_not_global
1 failed, 3 passed, 19 deselected in 0.28s
reverted -> c66ca83db6bd82c2edde1c2f079ebdeab4dc494a188ba7a78950b4b8cc4c6fa3

--- PLANT 3: the VALUE drifts, PDConfig.quiescent_resistance 0.50 -> 0.25 ---
    (out-of-tree pytest plugin: core/config.py is outside this agent's write
     territory; the effect on every PDConfig() is identical)
PLANT 3 applied: PDConfig.quiescent_resistance 0.50 -> 0.25
PLANT REACHED: PDConfig().quiescent_resistance = 0.25
E       AssertionError: quiescent_resistance moved; results/cell_state_multipliers_2026-09-20.md
        and the 2.00x bound in results/kill_calibration.md section 5 are now wrong
E       assert 0.25 == 0.5
E        +  where 0.25 = PDConfig(hypoxia_resistance=0.35, quiescent_resistance=0.25,
                                  max_kill_rate=1.0).quiescent_resistance
E           assert 4.0 == 2.0 ± 2.0e-12
E             Obtained: 4.0
E             Expected: 2.0 ± 2.0e-12
FAILED tests/test_cell_state.py::test_the_shipped_multiplier_values_are_pinned_to_the_published_claim
FAILED tests/test_cell_state.py::test_relabelling_bound_is_one_over_the_smallest_reachable_multiplier
plugin removed (no repo file was touched)

--- after all plants reverted ---
c66ca83db6bd82c2edde1c2f079ebdeab4dc494a188ba7a78950b4b8cc4c6fa3  tmesim/cells/pkpd.py
....                                                                     [100%]
4 passed, 19 deselected in 0.04s
```

Plant 3 is the one that matters most: it is not a hypothetical. It is precisely what a future
recalibration of the undefended M1 would look like, and it turns the published 2.00x into 4.00x
silently. The guard is what makes that visible.

Whole file, no plants: `23 passed in 22.08s`.

### 7.3 Reproduction

* Sections: `data_cache/xenium_gbm` (dx = 12 um), `data_cache/xenium_gbm_io` and
  `data_cache/gse324268_xenium` (dx = 30 um, for compute), `geometric` tissue mask, seed 0,
  `PDEConfig(**validation.delivery_metrics._PDE)`, drugs compiled through
  `validation.delivery_metrics.compile_real_drug(name, alpha)`.
* The swap path is gated against `sim.engine.simulate_response` bit-for-bit on every
  (section, drug) before any swept number is read.
* Nothing in `tmesim/` was modified. The only tracked files this work touches are
  `tests/test_cell_state.py` and this report.

### 7.4 The `calibrate_alpha()` assertions, verbatim

```
calibrate_alpha START = 22.223477711172563  ASSERT OK
...
calibrate_alpha END = 22.223477711172563  ASSERT OK
```

Printed and asserted at the start and the end of every run behind this report: the flagship
sweep, both ~800k-cell section sweeps, and the 64-seed endpoint study.

## 8. Where the next hour of cell-state work should go

1. **Measure M1, or stop quoting absolute kill.** It is the only multiplier with real leverage on
   the flagship section, it has no measurement behind it anywhere, and one small matched-exposure
   experiment (§2) would settle it. Everything else in this area is secondary to that.
2. **Pull the Ahmed 2018 supplementary cisplatin fold** and run it through M3's existing
   conversion. Cheap, bounded, and it retires M2 for one of the two drugs it governs. But note
   §3.1: on the flagship section this buys **zero** movement, and it only starts to matter on
   sections with a populated HYPOXIC band.
3. **Re-quote M3 as an interval.** `[0.675, 1.000]` wherever `x0.80` is load-bearing.
4. **Fix `sim/mechanism.py:240`** to call `effective_hypoxia_resistance` (§6.1). One line; the
   value is already in the diagnostics dict.
5. **Do not read "2.00x" as a property of the engine.** It is temozolomide's, it is the quiescent
   label's, it is an expectation not a realisation, and it moves with the very prior it is being
   used to reassure about.

---

### Sources

Literature in this report is second-hand from
`results/sweep/state_multiplier_literature_check_2026-07-29T05-54-00Z.md`, which performed the
primary search. Nothing new was retrieved here; the `R(r, f, h)` arithmetic in §1.1 is this
report's own recomputation of `PDConfig`'s published conversion.

* Musah-Eroje A & Watson S. *J Neurooncol* 2019;142(2):231-240. PMID **30694423** — the
  matched-oxygen 2-D null result.
* Wang H, Wan X, Xiong Y, *et al.* *Sci Rep* 2017;7(1):7965. PMID **28801626** — the 1.1–1.4x TMZ
  IC50 shifts that `TMZ_HYPOXIA_RESISTANCE` is derived from.
* Ahmed EM, Bandopadhyay G, Coyle B, Grabowska A. *Cell Oncol (Dordr)* 2018;41(3):319-328.
  PMID **29492900** — the un-extracted GBM cisplatin fold (§2).
* Rohwer N & Cramer T. *Drug Resist Updat* 2011;14(3):191-201. PMID **21466972** — mechanistic
  support for the *direction* of M1 only.
* Desoize B & Jardillier J. *Crit Rev Oncol Hematol* 2000;36(2-3):193-207. PMID **11033306** —
  names a high quiescent fraction as the leading multicellular-resistance mechanism.
