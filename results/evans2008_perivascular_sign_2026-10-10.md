# Evans 2008 vs the engine: does the oxygen field reproduce the measured per-vessel gradient SIGN SPLIT?

**Date: 2026-10-10. Status: PRE-REGISTERED. §1–§4 below were written and committed to this
file BEFORE a single engine number was produced.** Nothing under `tmesim/`, `app/`, `docs/`
or `landing-3d/` was modified. No default was changed. No parameter was tuned to any anchor.

**Calibration guard, start of pass:**

```
$ .venv/bin/python -c "from tmesim.validation.calibration import calibrate_alpha; a=calibrate_alpha(); assert a==22.223477711172563, a; print('alpha OK',a)"
alpha OK 22.223477711172563
```

*(Closing guard in §10.)*

---

## 1. What is being tested, and why this anchor

`results/data_scouting_threads_closed_2026-10-09.md` §6 promoted **Evans et al. 2008**
(PMID 19138031, PMC2629387, *Radiat Res* 170:677-690, DOI 10.1667/RR1207.1) to the top of
the scouting queue as *"the only item in either file that is adverse to a claim the engine
currently makes"* and *"the cheapest falsification available"*. This pass acts on it.

**The engine claim under test.** `tmesim/physics/oxygen.py` solves a steady-state
reaction-diffusion balance

```
D·lap(O) + supply·vasc·(1 − O) − consumption·rho_tum·O = 0
```

on the reconstruction grid, with `vasc` = the discrete segmented-vessel source
(`reconstruct.vessel_network.discrete_vascular_field`, nonzero only within
`VESSEL_RADIUS_UM = 7.5 µm` of a segmented vessel's endothelial cells and **exactly zero**
elsewhere), and `rho_tum` = the local tumour-cell areal density. `O = 1` is DEFINED to be
`REFERENCE_PO2_MMHG = 35.0 mmHg` (normal human brain tissue pO₂, Erecinska & Silver 2001,
PMID 11718758). `label_hypoxia` marks cells HYPOXIC where `O < 0.2` (≈ 7 mmHg), and
`hypoxic_tumor_fraction` (`tmesim/sim/engine.py:429`, via
`validation/delivery_metrics.py:574`) is the mean of that label over tumour cells. So the
whole hypoxia output of this engine is a functional of **one scalar field whose only spatial
driver is distance from a vessel, weighted by local tumour density.**

**The measurement.** Evans measured EF5 binding (%CRB) in 10 µm-wide radial bands around
**each individually identified CD31+ vessel** ("distance from each observed vessel", DFEOV)
in six human brain tumours, and fitted a line to the band means. Its conclusion:
*"neither distance to the nearest blood vessel nor distance from each observed blood vessel
provide reliable estimates of tissue pO(2)"*, because *"the presence of both positive and
negative slopes tended to average to near zero"*.

**Why this is adverse.** The engine's field cannot be anything other than a function of
distance-from-vessel. If the engine produces a reliable monotone fall of pO₂ away from
vessels where a human measurement finds the slope sign split, that is a statement about the
engine's **structure**, not about a parameter value.

---

## 2. The anchor numbers, verified from the primary source

Full text retrieved **2026-10-10** via `efetch?db=pmc&id=PMC2629387&retmode=xml`
(162,656 B), parsed with `experiments/evans2008/dump_tables.py` and `dump_text.py`. All five
tables are machine-readable. **Verification of the numbers reported in
`data_scouting_threads_closed_2026-10-09.md` §3.4 is recorded in §5 of this file, including
two corrections.** Nothing below is digitised from a figure; no figure was digitised
anywhere in this pass.

**Table 1 (the conversion ladder), verbatim:**

| Descriptor | %CRB | pO₂ (mmHg) |
|---|---|---|
| Physiological oxygenation | 1–3 | 100–20 |
| Mild hypoxia | 3–10 | 20–4 |
| Moderate hypoxia | 10–30 | 4–0.75 |
| Severe hypoxia | 30–100 | 0.75–0 |
| Anoxia | 100 | 0 |

**Table 4 (the countable statistic), verbatim — three tissue sections, vessel counts by
slope sign:**

| Specific oxygen level | left + | left − | middle + | middle − | right + | right − |
|---|---:|---:|---:|---:|---:|---:|
| Oxygenated | 0 | 0 | 0 | 0 | 5 | 4 |
| Mild hypoxia | 12 | 12 | 4 | 6 | 10 | 3 |
| Moderate hypoxia | 5 | 1 | 2 | 0 | 1 | 2 |
| Severe hypoxia | 1 | 0 | 0 | 0 | 0 | 1 |
| **Total vessels examined** | **18** | **13** | **6** | **6** | **16** | **10** |

**Pooled: 40 positive / 29 negative = 69 vessels, 0.5797 positive-slope fraction,
Clopper-Pearson 95% CI [0.4548, 0.6976].** Per section: 18/31 = 0.581 [0.391, 0.755];
6/12 = 0.500 [0.211, 0.789]; 16/26 = 0.615 [0.406, 0.798].

---

## 3. The mapping: EF5 %CRB → pO₂ → the engine's `O`. **The sign test needs no functional form.**

**Step 1 — the sign convention, from the paper's own words.** *"While a positive EF5 binding
gradient implies, as described above, a negative oxygen gradient (expected), the implication
of a negative EF5 binding gradient away from a blood vessel is a positive oxygen gradient
(unexpected)"*. So:

| Evans | implies | engine equivalent |
|---|---|---|
| **positive** EF5 slope (40 vessels) | pO₂ **falls** away from the vessel | `d⟨O⟩/dr < 0` — the EXPECTED diffusion gradient |
| **negative** EF5 slope (29 vessels) | pO₂ **rises** away from the vessel | `d⟨O⟩/dr > 0` — an ANTI-gradient |

**Step 2 — why the sign comparison is defensible without the conversion.** Table 1 is
strictly monotone *decreasing* in pO₂ across its entire tabulated range. `O = pO₂ / 35` is
a strictly increasing linear map. **The sign of a slope is invariant under any strictly
monotone transform of the dependent variable**, therefore

```
sign( d%CRB/dr )  =  − sign( dpO₂/dr )  =  − sign( dO/dr )
```

exactly, with **no assumption about the shape of the EF5–pO₂ curve, no calibration constant
and no use of `REFERENCE_PO2_MMHG`.** This is the reason the sign split — and not the slope
magnitude — is the primary statistic of this pass. It is the one quantity in Evans 2008 that
transfers to this engine without a modality conflation of the kind
`external_anchor_survey.md` §0.1 catalogues.

**Step 3 — the magnitude mapping, which IS assumption-bearing and is therefore REPORTED,
NOT SCORED.** Table 1's four finite class boundaries fit a power law
`pO₂ = 99.69 · %CRB^(−1.4272)` with relative residuals −0.3 %, +3.9 %, −6.8 %, +3.6 %
(fit computed from the four tabulated (%CRB, pO₂) pairs; the 100 %CRB → 0 mmHg point is
excluded because it is a limit, not a datum). Converting Evans' Table 3 slopes
(%CRB·µm⁻¹ × 1000) into mmHg·µm⁻¹ requires this fit **and** a choice of operating point,
so every magnitude number in §8 carries that caveat explicitly.

### What the mapping assumes — stated, not buried

1. **EF5 binding is inversely monotone in pO₂.** The paper's own premise, and it states
   *"EF5 binding has always been found to vary inversely with pO₂, with a 50–200-fold
   dynamic range"*. The sign test depends on this and nothing else.
2. **EF5 integrates; the engine does not.** Binding accumulates over the 16–24 h between
   IV dose (21 mg/kg) and resection. The engine's `O` is an *instantaneous steady state*.
   A time-average of a fluctuating pO₂ is not a steady state, and this is the mapping's
   largest assumption. **It cuts against the refutation**: temporal fluctuation in the
   tissue could produce measured anti-gradients that no instantaneous field would show.
   §9 states exactly how far this limits any verdict reached here.
3. **Table 1 is an in-vitro human-cell calibration applied to tissue**, after corrections
   the paper itself quantifies as spanning several orders of magnitude (Table 5: optics
   2.8×, AUC 1.9×, CRB 7.2×, competed-stain combined ~1000×).
4. **`O = 1 ≡ 35 mmHg` is a repo DECLARATION, not a measurement** (`oxygen.py`, "MMHG
   REFERENCE POINT"). It is used only in §8's magnitude arm and in the class assignment,
   never in the primary sign test.
5. **2D section vs 3D tissue.** Vessels outside the section plane exist in the tissue and
   do not exist in the engine. The paper says so itself: *"the tissue pO₂ will clearly be
   affected by vessels outside the plane and field of analysis"*. §9 bounds the verdict for
   this too.

---

## 4. PRE-REGISTRATION

### 4.1 The replication

For each section, the engine's own oxygen field is put through Evans' DFEOV algorithm:

1. `micro = reconstruct_real(source='xenium', sample=<cache>, dx_um=dx)` at shipped
   defaults (`oxygen_mode='pde'`, `state_mode='oxygen'`, `vascular_mode='discrete'`,
   `tissue_mask_mode='geometric'`).
2. The field measured is **the engine's own**: `solve_oxygen(micro, OxygenConfig())` is
   called explicitly and asserted **bit-identical** to `micro.fields.oxygen`. If that
   assertion fails, the pass reports the failure and produces no verdict.
3. `vessels = segment_vessels(micro.cells, micro.domain)` at the shipped
   `LINK_DISTANCE_UM = 20.0`.
4. Every tissue-mask grid-cell centre is assigned to its **nearest** vessel (nearest
   endothelial cell), which implements Evans' *"contours and bands are truncated when they
   intersect those from neighbouring vessels... pixels were assessed only once"*.
5. `r = (distance to that vessel's nearest endothelial cell) − VESSEL_RADIUS_UM`, i.e.
   distance from the painted vessel **periphery**, matching Evans' "from the vessel
   perimeter". Pixels with `0 ≤ r < 50 µm` are kept and binned into **5 bands of 10 µm**,
   band centres **5, 15, 25, 35, 45 µm** — the paper's actual analysis geometry (§5
   corrects the record on this).
6. Per vessel: band means of `O`, OLS slope against band centre, over **bands 1–5**
   (primary) and **bands 2–5** (the paper's second arm).
7. Vessel oxygenation class from the **band-1 mean**, through Table 1's pO₂ boundaries:
   physiological ≥ 20, mild [4, 20), moderate [0.75, 4), severe (0, 0.75), anoxia = 0 mmHg.
8. **Occupancy rule:** a vessel is analysable iff each of its 5 bands contains **≥ 3** grid-cell
   centres. The `≥ 1` rule is reported alongside as a sensitivity.

### 4.2 Sections, resolution, and arms — fixed in advance

**Sections: the three Xenium sections and only those** — `xenium_319`, `xenium_io`,
`gse324268`. The four Visium sections are **excluded by `app.server.api.SPOT_SPATIAL_REFUSAL`**:
their endothelial positions are synthesised (12 pseudo-cells per spot at the spot centroid),
so a per-vessel perivascular statistic computed on them is not a reading of that tissue.
Running them would manufacture agreement or disagreement out of pseudo-cells.

**Resolution: dx = 8 µm is the primary arm**, because `results/vessel_promotion.md` §2 is
this repo's own standing finding that *"distance-resolved perivascular claims still need
dx ≤ 8 µm"*. `dx = 12 µm` (production) and `dx = 5 µm` (finer) are sensitivity arms. A
section that cannot be reconstructed or solved at a given dx is recorded with its reason and
**never substituted**.

**Arms reported: `sink_form='linear'` (ships) is primary; `sink_form='michaelis_menten'`
(opt-in, `OxygenConfig(sink_form='michaelis_menten')`) is a second arm**, because
`results/oxygen_sink_form_2026-09-20.md` already established these two forms bracket the
measured pO₂ from opposite sides and do not overlap. Vessel subsets: **all segmented
vessels** is primary; **vessels with ≥ 4 endothelial cells** is a sensitivity arm (closer to
Evans' Table 2 mean vessel areas of 416–3385 µm², all larger than a single 7.5 µm-radius
disc's 177 µm²). Excluding small vessels is a choice that could move the answer, so it does
not get to be primary.

### 4.3 Primary statistic and the decision rule

**P_exp = (number of analysable vessels with `d⟨O⟩/dr < 0` over bands 1–5) / (number of
analysable vessels)**, pooled over the three Xenium sections at dx = 8 µm, linear sink, all
vessels. Evans' comparator is **0.5797, CP95 % CI [0.4548, 0.6976]**.

| branch | condition | what it means |
|---|---|---|
| **R1 — REFUTED, TOO CLEAN** | `P_exp ≥ 0.90` **and** `P_exp > 0.6976` | The engine manufactures a near-deterministic radial gradient where the measurement finds the sign split ~58/42. The engine's spatial oxygen claim is refuted on the one statistic that transfers without a conversion. **Leads the report.** |
| **R2 — CONSISTENT** | `0.4548 ≤ P_exp ≤ 0.6976` | The engine reproduces the measured sign heterogeneity. No refutation on this statistic, and the anchor is answered rather than dodged. |
| **R3 — DISCREPANT, NOT DETERMINISTIC** | `0.6976 < P_exp < 0.90` | A quantified discrepancy. **Explicitly NOT a refutation**, declared in advance so this outcome cannot later be upgraded into one. |
| **R4 — REFUTED, TOO NOISY** | `P_exp < 0.4548` | The engine is *more* disordered than the tissue — a different failure, and it would mean the anti-gradients come from the tumour-density sink rather than from anything Evans measured. |
| **R5 — COMPARISON VOID** | fewer than 20 analysable vessels on every section, or the bit-identity assertion in §4.1.2 fails | No verdict. The void is the finding. |
| **R6 — CLASS STRATIFICATION VOID** | every analysable engine vessel falls in a **single** Table 1 oxygenation class | Evans' Table 4 is stratified by class and spans four of them on one patient; if the engine spans one, the stratified comparison cannot be made. **Fires independently of R1–R4 and is reported alongside whichever of those fires.** |

**Predicted direction, committed in advance: R1 fires, and R6 fires with it.**

The reasoning, so that a hit is not retro-fitted. (a) `vasc ≥ 0` and the source term is
`supply·vasc·(1 − O)`; the solution obeys a maximum principle with `O ≤ 1`, so the term is
**non-negative everywhere — no vessel in this engine can ever be an oxygen sink.** (b) Under
`amplitude='unit'` every vessel carries the **same** source amplitude per unit of covered
area, so there is no per-vessel perfusion or oxygenation heterogeneity in the model at all.
(c) The sink is non-negative. For an isolated vessel the steady state is therefore monotone
decreasing in `r` by construction, and the only two mechanisms that can flip a sign are
tumour-density heterogeneity and a neighbouring vessel. (d) R6 is predicted because
`results/oxygen_sink_form_2026-09-20.md` §0 already reports the shipped linear sink confining
tumour-cell pO₂ to **10.1669–14.7473 mmHg** on the flagship section — a 4.6 mmHg total spread
that lies wholly inside Table 1's single "mild hypoxia" class (4–20 mmHg).

**If the predicted branch fires it is a refutation of the engine's spatial oxygen claim, not
a parameter miss, only to the extent that §9's two limits (temporal integration, 2D vs 3D)
permit.** Those limits are written here, before the run, precisely so that they cannot be
quietly dropped from the conclusion afterwards.

### 4.4 Hard constraints this pass runs under

* **No parameter is tuned to this anchor, in any arm.** Every `OxygenConfig` value used is
  either the shipped default or the documented opt-in `sink_form`.
* Nothing under `tmesim/`, `app/`, `docs/`, `landing-3d/` is edited.
* No existing file in `results/` or `archive/` is rewritten. This is one new dated file.
* No absolute kill fraction, no clinical or diagnostic framing, appears anywhere.
* No figure was digitised. Every Evans number is from a table in the PMC XML.

---

*Sections 5 onward were written after the run.*

---

## 5. Anchor verification. Three things verified, two corrected, one open item dissolved.

Source of record: `efetch?db=pmc&id=PMC2629387&retmode=xml`, retrieved 2026-10-10,
162,656 B. Tables printed by `experiments/evans2008/dump_tables.py`, body and figure
captions by `experiments/evans2008/dump_text.py`. **No figure was digitised.**

### 5.1 Verified exactly, as reported in `data_scouting_threads_closed_2026-10-09.md` §3.4

| claim in the record | verdict |
|---|---|
| Five machine-readable tables | **VERIFIED** — Tables 1–5 all parse |
| Table 1 ladder: 1–3 %CRB = 100–20 mmHg; 3–10 = 20–4; 10–30 = 4–0.75; 30–100 = 0.75–0; 100 = 0 | **VERIFIED verbatim** |
| Table 2 vessels/mm² = 244, 70, 16, 127, 24, 11 | **VERIFIED** (patients A–F) |
| Table 2 vessel area range 416 ± 878 … 3385 ± 5929 µm², **SD not SEM** | **VERIFIED** — footnote-free `± SD` in the header |
| Table 3 carries per-patient, per-class intercept **and** slope for bands 1–5 and 2–5 | **VERIFIED** — 13 populated rows, slope units `%CRB per µm × 1000` |
| Table 4 slope-sign counts 18 v 13, 6 v 6, 16 v 10 | **VERIFIED** |
| Table 5 correction factors: optics 2.8×, AUC 1.9×, CRB 7.2× | **VERIFIED**, and the record omitted the fourth column: competed-stain background and mean signal, *"Combined ∼1000"*. A ~1000× correction sits between the raw fluorescence and the %CRB scale. |
| Dose 21 mg/kg IV | **VERIFIED** |
| Six tumours: hemangiopericytoma, WHO Grade 2 glial, WHO Grade 3 glial, + glioblastomas | **VERIFIED**; precisely: A HPCTA, B oligoastrocytoma (G2), C anaplastic oligodendroglioma (G3), D and E *de novo* glioblastoma, F *recurrent* glioblastoma |

Additionally verified and not in the record: **20 tissue sections** were analysed across the
six patients (Table 2 `N` = 4, 3, 4, 3, 3, 3), after four slides were excluded (two >80 %
necrotic, two with antibody-microprecipitation artefact); the pixel is **1.75 µm** linear /
**3.1 µm²** areal.

### 5.2 CORRECTION 1 — Table 4 is **one patient**, not six, and the record does not say so

The record presents the sign counts as the paper's strongest asset without naming their
scope. They are **patient F only** — the recurrent glioblastoma — across its **three**
slides. Fig. 4's caption: *"The recurrent glioblastoma (patient F, lower panel) has a
relatively small number of vessels (n = 69, also see Fig. 5, Table 4)"*. Fig. 5's caption:
*"vessel gradients assessed on three slides from the recurrent glioblastoma from patient
F… Sixty-nine vessels were examined (also see Tables 2, 4)"*.

**Table 4 reconciles to patient F's Table 2 row exactly, class by class**, which settles the
scope without relying on the captions:

| class | Table 4, summed over the three sections and both signs | Table 2, patient F |
|---|---:|---:|
| oxygenated | 0+0 + 0+0 + 5+4 = **9** | **9** |
| mild hypoxia | 12+12 + 4+6 + 10+3 = **47** | **47** |
| moderate hypoxia | 5+1 + 2+0 + 1+2 = **11** | **11** |
| severe hypoxia | 1+0 + 0+0 + 0+1 = **2** | **2** |
| **total** | 31 + 12 + 26 = **69** | **69** |

No other patient's Table 2 row sums to 69, and patients A, B and C have zero vessels in any
hypoxic class, so the 47 / 11 / 2 counts can come from nowhere else.

**This matters for the comparison, in the direction that weakens the anchor**: the
countable statistic is *n* = 69 vessels in *one* recurrent glioblastoma, not 69 vessels
across six tumours. The broader *"positive and negative values averaged to near zero"*
claim does rest on all six patients (Table 3), but it is an average, not a count.

### 5.3 CORRECTION 2 — there is no 0–200 µm profile, so a listed open item dissolves

The record states the spatial axis is *"20 10-µm-wide bands"*, concludes that *"Table 3
tabulates only the first 0–50 µm; the full 200 µm profile is figure-only"*, and carries
*"Evans' 50–200 µm band profile (figure-only)"* forward as an open item.

The 20-band construction belongs to **Fig. 2, which is a feasibility test, not a binding
profile**: its Y axis is *"the percentage of vessels with at least some portion of each band
unencumbered by neighboring vessels"*. The analysis geometry is stated in Results:
*"**A maximum perivascular radius of 50 µm employing five equally spaced contours was
selected** since the use of five contours of 10 µm spacing compared to 10 contours of 5 µm
spacing did not change the EF5 binding gradients but did significantly reduce the
computation time"*, and *"Each vessel was thus characterized by average EF5 binding in the
five bands with average distances of 5, 15, 25, 35 and 45 µm"*.

So **no 0–200 µm EF5 profile exists in this paper, in a table or in a figure.** Table 3 is
not a truncation of a longer published profile; 0–50 µm is the whole published analysis.
The open item is not blocked — it has no referent. (Fig. 3's caption does mention one
qualitative ~200 µm observation about two specific vessels, which is an annotation on an
image, not a profile.)

### 5.4 One clarification that bears on the mapping

The record quotes EF5 as given *"approximately 24 h"* before surgery, which is the
abstract's wording. The Methods are wider: *"Sixteen to 24 h after completion of drug
administration, the tumor was resected"*. The binding therefore integrates pO₂ over a
**16–24 h** window, not a fixed one. This widens assumption 2 in §3.

---

## 6. The engine path, as found

`hypoxic_tumor_fraction` is produced by exactly one chain, confirmed by reading it:

`reconstruct.real_pipeline.reconstruct_real` → `reconstruct.real_fields.vascular_density_from_cells(mode='discrete')` → `reconstruct.vessel_network.segment_vessels` + `discrete_vascular_field` (discs of `VESSEL_RADIUS_UM = 7.5 µm` around measured endothelial cells, area-weighted, **exactly zero elsewhere**) → `physics.oxygen.solve_oxygen` (one sparse direct solve of the linear-sink steady state) → `physics.oxygen.label_hypoxia` / `cells.agents.update_cell_states` banding at `O < 0.2` → `validation.delivery_metrics.py:574` `hypoxic_tumor_fraction = mean(hypoxic)` over tumour cells → `sim/engine.py:429`.

**The field this pass measured is the engine's own.** On every run,
`fields.oxygen` was asserted **bit-identical** to `solve_oxygen(micro, OxygenConfig())`
(`pipeline_field_is_bit_identical_to_solve_oxygen = true`, 7/7 runs). Pre-registered branch
R5 did not fire.

**Two structural properties, verified at runtime rather than argued:**

1. **`O ≤ 1` everywhere on every run** (`field_o2_exceeds_one_anywhere = false`, 7/7). The
   source term is `supply·vasc·(1 − O)`, so it is **non-negative everywhere**:
   **no vessel in this engine can be an oxygen sink.** Evans' own proposed explanation for
   negative EF5 gradients — *"groups of vessels or vessel networks transporting blood into
   the tumor for nutrient supply (in-flowing) and out of the tumor for waste removal
   (out-flowing)"*, i.e. *"a significant fraction of them must be considered to be 'sinks'"*
   — is **not representable in this engine at any parameter value, in 2D or in 3D.**
2. **Every vessel has the same source amplitude per unit of covered area**
   (`amplitude='unit'`). There is no per-vessel perfusion, hematocrit or oxygenation
   heterogeneity in the model at all. Vessel-to-vessel differences enter only through
   geometry.

Property 1 is the same defect `results/transport_falsifier_scouting_2026-10-09.md` §(Chary &
Jain 1989, PMID 2748592) raises for the **advective** term — *"the engine treats segmented
vessels as a source with no corresponding sink — so the sign of the advective term near a
subset of vessels may be wrong"*. This pass confirms the identical structure in the
**oxygen** term and supplies the physiological reason it matters.

**Attribution correction to this pass's own brief.** The brief attributes the oxygen-sink
*"overshoots rather than explains"* finding to
`results/transport_falsifier_scouting_2026-10-09.md`. That file contains **no** mention of
oxygen (0 hits for `oxygen`, 0 for `overshoot`). The finding is
`results/oxygen_sink_form_2026-09-20.md`, whose title is *"The sink's functional form does
NOT explain the missing hypoxia — it overshoots it, and the two forms' pO₂ supports do not
overlap"*. Both files are connected here, each to what it actually says.

**Cross-check against that prior.** It reports the shipped linear sink confining tumour-cell
pO₂ to 10.1669–14.7473 mmHg and the Michaelis–Menten sink to 0.3370–3.7775 mmHg. This pass
measures the **whole-field** range on the same section as **10.162–14.768** and
**0.337–3.801** mmHg. The prior reproduces.

---

## 7. THE RESULT ON THE PRE-REGISTERED RULE: **R3 fired. The predicted branch did not.**

Primary arm as registered — three Xenium sections, dx = 8 µm, `sink_form='linear'`, all
segmented vessels, `min_count = 3`:

| section | vessels/mm² | analysable | expected sign (`dO/dr < 0`) | anti-gradient | **P_exp** | CP95 |
|---|---:|---:|---:|---:|---:|---|
| `xenium_319` | 110.8 | 2,805 | 2,520 | 285 | **0.8984** | [0.8866, 0.9093] |
| `xenium_io` | 233.5 | 7,308 | 5,887 | 1,421 | **0.8056** | [0.7963, 0.8146] |
| `gse324268` | 320.4 | 4,962 | 3,781 | 1,181 | **0.7620** | [0.7499, 0.7738] |
| **POOLED** | — | **15,075** | **12,188** | **2,887** | **0.8085** | **[0.8021, 0.8147]** |
| **Evans patient F (measured)** | **11** | 69 | 40 | 29 | **0.5797** | [0.4548, 0.6976] |

```
PRE-REGISTERED BRANCH: R3 DISCREPANT-NOT-DETERMINISTIC
  (0.6976 < 0.8085 < 0.90)
```

**The prediction committed in §4.3 was R1, and R1 did not fire.** `P_exp = 0.8085` is
clearly outside Evans' interval — the two CP95 intervals are **disjoint**, 0.8021 > 0.6976 —
so the discrepancy is not a sampling accident. But 0.8085 < 0.90, and §4.3 declared in
advance that this band is *"explicitly NOT a refutation… declared in advance so this outcome
cannot later be upgraded into one"*. **It is not upgraded. The scored verdict of this pass
is R3.**

**R6 also did not fire, and it missed by two vessels.** The rule required *every* analysable
engine vessel in a single Evans class. Pooled, **15,073 of 15,075 (99.987 %) are
`mild_hypoxia`**; `xenium_io` puts **2** vessels in `moderate_hypoxia`. On `xenium_319`
alone R6 *would* have fired (one class, 2,805/2,805). **The rule as written was
under-specified** — it did not say whether "every" meant per-section or pooled — and that is
a defect in this pass's own pre-registration, recorded rather than resolved in whichever
direction helps. Read pooled, R6 does not fire. Read per-section, it fires on one of three.
Either reading, the number to carry is **99.987 % in one class**.

### 7.1 Sensitivity arms, all pre-registered in §4.2

| arm | P_exp | branch it would fall in |
|---|---:|---|
| **primary** (dx 8, linear, all vessels, mc 3) | **0.8085** | **R3** |
| bands **2–5** (the paper's second arm), pooled | 0.7690 | R3 |
| `xenium_319` at dx = **12 µm** (production) | 0.8968 | R3 |
| `xenium_319` at dx = **8 µm** | 0.8984 | R3 |
| `xenium_319` at dx = **5 µm** | 0.8966 | R3 |
| `min_count = 1` instead of 3 (`xenium_319`, dx 8) | 0.8954 | R3 |
| vessels with **≥ 4 endothelial cells**, pooled | **0.9048** | R1 |
| `sink_form='michaelis_menten'` (opt-in; `xenium_319`, dx 8) | **0.9508** | R1 |

**The statistic is resolution-invariant** — 0.8968 / 0.8984 / 0.8966 at dx = 12 / 8 / 5 µm
on the same section, a spread of 0.002. So R3 is not a grid artefact, and the repo's own
`dx ≤ 8 µm` requirement for distance-resolved perivascular claims
(`results/vessel_promotion.md`) is satisfied without changing the answer.

**Two arms cross into R1, and neither gets to be the verdict.** The ≥ 4-endothelial-cell
subset (0.9048) is the arm closest to Evans' vessel sizes — §4.2 made it a sensitivity arm
precisely because *"excluding small vessels is a choice that could move the answer, so it
does not get to be primary"*, and it moved the answer. The Michaelis–Menten arm (0.9508) is
not the shipped default. Both are reported; neither is promoted.

### 7.2 The vessel-density mismatch, and why it does not rescue the engine

Evans' patient F — the only patient with a sign count — has **11 vessels/mm²**, the
*sparsest* of his six. The three engine sections carry **110.8, 233.5 and 320.4
vessels/mm²**: 10–29× denser, and `gse324268` is denser than Evans' densest patient
(A, 244/mm²). This is the largest like-for-like defect in the comparison and it is not
repairable from the data in hand.

It does **not** work in the engine's favour. `P_exp` falls monotonically as vessel density
rises across the three sections (110.8 → 0.8984, 233.5 → 0.8056, 320.4 → 0.7620); a linear
read of those three points extrapolates to **P_exp ≈ 0.96 at Evans' 11 vessels/mm²** — i.e.
*further* from the measurement, in R1 territory. The mechanism is visible: denser vessels
mean more neighbour interference, which is the engine's only route to an anti-gradient at
all. **Matching the density would be expected to make the discrepancy worse, not better.**
(Three points and a linear extrapolation over a 10× span is weak evidence and is offered
only as a direction, not a value.)

There is also a correspondence worth recording: Evans' Fig. 2 reports that at patient A's
244 vessels/mm² *"less than 10 % of vessels can be analyzed beyond 50 µm"*. Under this
pass's `min_count ≥ 3` rule the analysable share is 90.7 % / 35.0 % / 15.7 % of segmented
vessels at 110.8 / 233.5 / 320.4 vessels/mm² — the same steep fall with density, from an
independent implementation.

---

## 8. Reported, not scored: magnitude and class spread. **This is where the engine fails.**

Both comparisons below were designated **reported-not-scored** in §3 step 3 and §4.2, before
the run, because each depends on the Table 1 conversion or on the `O = 1 ≡ 35 mmHg`
declaration. Neither changes the §7 verdict. Both are large.

### 8.1 The engine's perivascular gradient is 19–54× shallower than Evans' *averaged* lines

Converting all 13 populated Table 3 rows through the Table 1 power law
(`experiments/evans2008/table3_magnitudes.py`, chain rule at each row's own intercept):

* **Patient A is excluded**: its intercept of 0.3 %CRB lies *outside* the 1–30 %CRB range
  the power law was fitted on, and extrapolating there returns 556 mmHg, which is
  physiologically impossible. Reporting that row would be fitting an artefact.
* The remaining **12 rows** give `|dpO₂/dr|` = **0.00028 to 0.825 mmHg/µm, median 0.0507**.
* **The conversion is violently nonlinear and that is a real limit on this arm**: the same
  %CRB slope maps to 0.0102 mmHg/µm at a 5.9 %CRB intercept and to 0.825 mmHg/µm at a
  1.0 %CRB intercept. These numbers are order-of-magnitude statements only.

The engine, same geometry, same 5 bands:

| section | median `|slope|` (mmHg/µm) | median `pO₂(band1) − pO₂(band5)` | max `|band1 − band5|` |
|---|---:|---:|---:|
| `xenium_319` | 0.000939 | **+0.036 mmHg** | 0.617 mmHg |
| `xenium_io` | 0.002616 | **+0.090 mmHg** | 1.774 mmHg |
| `gse324268` | 0.002348 | **+0.072 mmHg** | 0.915 mmHg |
| pooled | **0.00202** | **+0.080 mmHg** | — |

**The median total change in pO₂ across the entire 40 µm from band 1 to band 5 is
0.080 mmHg.** Against Evans' median converted slope the engine is **19–54× shallower**, and
Evans' Table 3 values are *averages over vessels*, about which the paper says explicitly:
*"the slopes for the individual vessels were much larger than those for the averages"*. So
19–54× is a **floor** on the discrepancy, not an estimate of it.

This is also resolution-invariant: median `|slope|` = 0.000904 / 0.000939 / 0.000967 at
dx = 12 / 8 / 5 µm on `xenium_319`.

### 8.2 The engine spans 2 of Evans' 4 perivascular classes, and 99.987 % of vessels sit in one

Using Evans' own Table 1 boundaries and the repo's own `REFERENCE_PO2_MMHG`:

| | oxygenated (≥20) | mild (4–20) | moderate (0.75–4) | severe (<0.75) | classes spanned |
|---|---:|---:|---:|---:|---:|
| **Evans patient F** (Table 2, *n* = 69) | **9** | **47** | **11** | **2** | **4** |
| engine pooled, linear sink (*n* = 15,075) | 0 | **15,073** | 2 | 0 | 2 |
| engine `xenium_319`, linear | 0 | 2,805 | 0 | 0 | 1 |
| engine `xenium_319`, Michaelis–Menten | 0 | 0 | 2,543 | 262 | 2 |

One recurrent glioblastoma produces perivascular band-1 values spanning from above 20 mmHg
to below 0.75 mmHg — a **≥27-fold** range — **within three slides of one tumour**. The
engine's band-1 values on 15,075 vessels across three tumours span **10.23–14.73 mmHg**
(`xenium_319`), a **1.44-fold** range. The shipped linear sink cannot place a vessel in
Evans' "oxygenated" class anywhere on any of the three sections, and the opt-in
Michaelis–Menten sink cannot place one there either — it moves the whole population down
two classes rather than spreading it.

**This comparison needs no slope and no power law.** It uses Table 1's class boundaries,
Table 2's class membership, and one declared constant. It is the cleanest adverse result in
this pass.

---

## 9. POST-HOC: why R3 and not R1. **The engine's anti-gradients are flatness, not mechanism.**

**Everything in this section is post-hoc. It was not pre-registered and it does not change
the §7 verdict** (`experiments/evans2008/posthoc_magnitude.py`). It exists because R3 can be
read as "the engine broadly agrees with Evans", and that reading is wrong.

A sign is only meaningful if the change it signs is resolvable. Pooled over the primary arm:

| resolvability floor on `|pO₂(band1) − pO₂(band5)|` | vessels kept | share | **P_exp** | CP95 |
|---|---:|---:|---:|---|
| ≥ 0 mmHg (**the primary arm**) | 15,075 | 100.0 % | **0.8085** | [0.8021, 0.8147] |
| ≥ 0.01 mmHg | 14,070 | 93.3 % | 0.8274 | [0.8210, 0.8336] |
| ≥ 0.05 mmHg | 9,815 | 65.1 % | 0.8689 | [0.8620, 0.8755] |
| ≥ 0.10 mmHg | 6,373 | 42.3 % | 0.8979 | [0.8902, 0.9052] |
| ≥ 0.25 mmHg | 1,725 | 11.4 % | 0.9565 | [0.9458, 0.9657] |
| ≥ 0.50 mmHg | 212 | 1.4 % | **0.9953** | [0.9740, 0.9999] |
| ≥ 1.00 mmHg | 7 | 0.05 % | 1.0000 | [0.5904, 1.0000] |

`P_exp` rises monotonically with resolvability to **0.9953 at ≥ 0.5 mmHg**. And the
anti-gradient vessels are systematically the flat ones: median `|slope|` in the
anti-gradient group is **2.2–2.8× smaller** than in the expected-sign group on every section
(0.000366 vs 0.001027; 0.001213 vs 0.003119; 0.001271 vs 0.002757 mmHg/µm).

**So the two 50/50-ish-looking results have opposite causes.**

* **Evans' near-zero average** comes from *large* slopes of both signs cancelling — the
  paper states it: *"The basis for these low average slopes was not that the slope for
  individual vessels was uniformly low. Rather, the presence of both positive and negative
  slopes tended to average to near zero"*, with Fig. 5 *"on a log scale"* to fit them all.
* **The engine's 19 % anti-gradients** come from a field that is nearly flat perivascularly:
  89.1 % of `xenium_319` vessels change by less than 0.1 mmHg across the whole 40 µm. Where
  the engine *does* produce a resolvable gradient it is 99.5 % the textbook sign.

R3 is therefore a **coincidence of summary statistics, not an agreement**. Recording it as
"the engine is broadly consistent with the measurement" would be precisely the defect
`codon-vacuous-green-checks` names: a check that passes without testing its subject.

---

## 10. What this refutes, what it does not, and what stays open

### 10.1 What it does NOT refute — stated first, because the pre-registered branch says so

* **The sign statistic does not refute the engine's spatial oxygen claim.** `P_exp = 0.8085`
  landed in the band §4.3 declared in advance to be a quantified discrepancy and *not* a
  refutation. The pass's own prediction (R1) was wrong, and the rule stands as written.
* **It does not establish that a 3D field would fail.** The engine is 2D; the measurement is
  a 2D section of a 3D tissue with vessels out of plane, and the paper says so:
  *"the tissue pO₂ will clearly be affected by vessels outside the plane and field of
  analysis"*. Some unknown part of Evans' 42 % anti-gradient fraction is out-of-plane
  vessels.
* **It does not establish that the measured anti-gradients would survive an instantaneous
  measurement.** EF5 integrates pO₂ over 16–24 h. A time-average of a fluctuating field is
  not a steady state, and this cuts *for* the engine: temporal fluctuation could produce
  measured anti-gradients that no instantaneous field would show. §3 assumption 2 was
  written before the run for exactly this reason.
* **It says nothing about drug transport.** `results/vessel_profile.md` measures a clean
  perivascular *drug* gradient in the same geometry, and O₂ is consumed orders of magnitude
  faster than a cytotoxic is cleared. Nothing here transfers to the delivery solver, and the
  engine must not borrow credibility in either direction.
* **No parameter was tuned, in any arm**, and no result here licenses tuning one. The 19–54×
  magnitude gap is not an invitation to move `consumption`:
  `results/oxygen_sink_form_2026-09-20.md` already showed that the one documented axis
  available overshoots the measured pO₂ from the other side.

### 10.2 What it does show

1. **The engine has no perivascular oxygen gradient worth the name.** Median total change
   across 40 µm from a vessel: **0.080 mmHg**. Resolution-invariant across dx = 12/8/5 µm.
   Its headline claim is not "a reliable monotone gradient" — it is **an almost flat field
   with a weak sign preference.** Both halves of that sentence are new.
2. **The engine cannot span the measured perivascular pO₂ range.** 2 of Evans' 4 classes,
   99.987 % of 15,075 vessels in one, against 4 classes in three slides of one tumour
   (§8.2). This needs no conversion beyond Table 1's own boundaries.
3. **A structural gap, verified at runtime, not argued.** `O ≤ 1` on 7/7 runs ⇒ the source
   term `supply·vasc·(1 − O)` is non-negative everywhere ⇒ **no vessel in this engine can
   be an oxygen sink**, and every vessel carries the same amplitude per unit area. Evans'
   proposed mechanism for negative gradients — in-flowing versus out-flowing vessels —
   **cannot be expressed in this model at any parameter value, in 2D or 3D.** This is the
   same structure `transport_falsifier_scouting_2026-10-09.md` flags for the advective term
   via Chary & Jain 1989. It is a model-form gap, not a parameter miss, and it is the one
   finding here that no amount of better data would change.
4. **The sign agreement is for the wrong reason (§9).** Post-hoc, and labelled as such.
5. **The anchor is cheaper and weaker than the record said.** Fully tabulated and free, as
   reported — but the countable statistic is one patient at 11 vessels/mm², against engine
   sections at 110.8–320.4, and the 0–200 µm profile the record lists as open does not
   exist.

### 10.3 Open, and named so it stays visible

* **Evans' per-vessel slope *values* are not published** — only Table 3's per-class averages
  and Table 4's signs. A magnitude test against *individual* vessels, which is the test the
  engine most clearly fails, would need Fig. 5 digitised (a log-axis scatter) or the
  authors' data. This pass digitised nothing.
* **No section in this repo matches patient F's 11 vessels/mm².** Until one does, the §7.2
  extrapolation is a direction, not a number.
* **Nothing here is registered in `tmesim/validation/scorecard.py`.** This pass added no
  scorecard entry, and the §8.2 class-spread comparison is the candidate that would be worth
  one — it is cheap, conversion-light and currently adverse.

---

## 11. Provenance and files

| file | role |
|---|---|
| `experiments/evans2008/dump_tables.py` | prints every `table-wrap` in the PMC XML |
| `experiments/evans2008/dump_text.py` | prints body sections and figure captions |
| `experiments/evans2008/perivascular_sign.py` | the DFEOV replication on the engine's field |
| `experiments/evans2008/table3_magnitudes.py` | Table 3 %CRB/µm → mmHg/µm via the Table 1 ladder |
| `experiments/evans2008/pool.py` | pools the sections and applies §4.3 as written |
| `experiments/evans2008/posthoc_magnitude.py` | §9 only, post-hoc, labelled |
| `experiments/evans2008/out/*.json` | every number in §7–§9 |

The PMC XML was downloaded to `/tmp/evans2008_src/` and parsed there; nothing was written to
`data_cache/` except read-only symlinks to the three already-present Xenium caches.

## 12. Closing calibration guard

```
$ .venv/bin/python -c "from tmesim.validation.calibration import calibrate_alpha; a=calibrate_alpha(); assert a==22.223477711172563, a; print('alpha OK',a)"
alpha OK 22.223477711172563
```

Asserted at the start and end of every driver in §11 as well, and unchanged throughout.
Nothing under `tmesim/`, `app/`, `docs/` or `landing-3d/` was modified; no existing file in
`results/` or `archive/` was edited; no parameter was tuned to any anchor; no figure was
digitised; no absolute kill fraction and no clinical framing appears anywhere above.

**Scored verdict of this pass: R3 — a quantified discrepancy on the sign statistic, not a
refutation. The pre-registered prediction was R1 and it was wrong. The adverse findings that
do stand are §8.2 (class spread), §10.2.1 (gradient magnitude) and §10.2.3 (vessels cannot
be sinks) — and only the last of those is beyond repair by data.**
