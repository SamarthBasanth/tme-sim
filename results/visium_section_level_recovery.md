# Is the spot-data refusal over-applied? — a jitter-reseed audit of the four Visium GBM sections

**Date: 2026-09-20. Status: MEASUREMENT.** Every number below is the output of a
command run in this pass. Nothing under `tmesim/`, `app/` or `landing-3d/` was
modified. No existing file in `results/` or `archive/` was altered. The driver is
`experiments/visium_section_level.py` (new, this pass); its raw record is a JSON
written outside the repository and the numbers quoted here are reproduced from
the driver's own stdout.

* commission: route **R12** of `results/verification_routes_2026-09-20.md` — *"Four
  of those Visium sections have been treated as unusable because the spatial-spread
  statistic is refused in code on spot data. That refusal is correct and applies to
  that statistic. It does not apply to section-level quantities."*
* the refusal under audit: `app/server/api.py::SPOT_SPATIAL_REFUSAL`
* `calibrate_alpha()` = **22.223477711172563** at the start and at the end of every
  stage (asserted by the driver; the run aborts if it moves).

**What this document does NOT do.** It produces no cross-compound ordering — that
claim is withdrawn (2026-08-21) and is untouched here; the driver runs ONE compound
(temozolomide) across sections and compares nothing between compounds. It produces
no absolute kill fraction as an output. It makes no clinical, diagnostic or
response claim. It does not weaken any refusal: every quantity that fails the
pre-registered criterion below stays refused, and two quantities that the
pre-registration wrongly cleared are moved BACK into the refused column in §3.2.

---

## 0. The answer

> ## The refusal's SUBJECT is right. Its stated REASON is wrong on this platform, and the recovery it is hiding is real but small.
>
> **1. Every quantity the refusal names is measurably contaminated, so none of them
> comes back.** Re-running the shipped pipeline with nothing changed but the
> pseudo-cell placement seed moves the fitted penetration length by up to **19.4 %**
> (`visium_gbm3`: 301.2 → 331.8 → 242.7 µm across three seeds), the worst-served
> exposure decile `E_p10` by up to **1.09 %**, and the headline delivery spread
> `E_p90/E_p10` by up to **0.60 %**. The refusal is not over-applied to its own list.
>
> **2. But the sentence that justifies refusing the exposure map is factually false
> for the Visium path.** The refusal says the exposure field's *"SOURCE is a density
> over ENDOTHELIAL CELL POSITIONS, and on a spot section those positions are
> synthesised"*. The Visium reconstructor never builds a source from cell positions.
> `reconstruct.real_fields.vascular_density_from_cells` was called **0 times** on all
> four Visium sections and **1 time** on the `xenium_319` positive control; the
> pipeline's own provenance string is `"KDE of REAL ENDOTHELIAL spot fraction"`
> against Xenium's `"DISCRETE segmented vessels from REAL ENDOTHELIAL cell
> positions"`; and the `vascular_density` array is **bitwise identical** across three
> jitter seeds on every section.
>
> **3. The verdict on the map survives anyway, by a route the text does not name.**
> With convection ON — the shipped default — the solved AUC field *is* jitter-
> dependent: **100 % of tissue pixels** change under a reseed (median 0.17 %, p99
> 1.76 %, max 2.82 %). The route is `physics/pressure.py`, which rasterises **TUMOUR**
> cell positions into a tumour-occupancy field to locate the leaky, poorly-drained
> region; the IFP moves by up to 2.625 and the Darcy velocity advects the drug. With
> `convection=False` the same field is **bitwise identical**. So the refusal reaches
> the right answer through a mechanism that is not the one it states, on the wrong
> cell type, with a ~2 % magnitude rather than the 5.7× distance-tail distortion
> `results/synthesised_positions.md` measured for the *distance* field.
>
> **4. Five section-level quantities do survive, and one of them is worth having.**
> `mean_vascular_density`, `mean_ecm_density` and `penetration_length_analytic_um`
> are **exactly** unchanged (0.0 deviation, bitwise-identical fields);
> `mean_tissue_auc` is unchanged to 3.0e-16 (1 ULP). But `tissue_coverage_fraction`
> is **1.0 on all four sections because the Visium reconstructor never sets a tissue
> mask at all** — it is the absence of a mask, not a measurement — and only
> **61.5–93.4 %** of the domain rectangle has a spot under it. And
> `mean_ecm_density` is **exactly 0.0** on `visium_gbm2`, because that section's
> deconvolved STROMAL fraction is identically zero at every one of its 2,891 spots.
>
> **5. The scout's claim that this is "free replication on the project's #1 priority"
> does not hold, and the free thing is somewhere else.** The #1 priority
> (`CLAIM_INVENTORY` §6.1) is that *almost every defensible claim is n = 1 in tissue*;
> those claims are C1, K1 and S1/S2, all per-cell or per-compound kill-panel claims
> on **Xenium** sections. None of the five recovered quantities backs any of them.
> What IS already paid for and unclaimed: `results/multisection_generalisation.json`
> has held all seven sections' payloads since **2026-08-23**, and
> `experiments/multisection_report.py` (838 lines, tracked) exists to render
> `results/multisection_generalisation.md` from it — **and that file has never
> existed**, in the working tree or anywhere in git history.

---

## 1. What the refusal actually refuses, and where it is enforced

### 1.1 The constant

`app/server/api.py:1594` — `_SPOT_SOURCES = ("visium",)`, matched case-insensitively
against a section record's `platform` (ingested) or `source` (built-in).

`app/server/api.py:1597` — `SPOT_SPATIAL_REFUSAL`, one string. It names **six**
quantities:

| # | quantity named in the text | mechanism the text gives |
|---|---|---|
| 1 | the per-cell attenuation spread | 12 pseudo-cells per spot at the centroid with uniform jitter |
| 2 | the immune-neighbourhood structure ratio | same |
| 3 | the delivery / exposure **map** | *"its SOURCE is a density over ENDOTHELIAL CELL POSITIONS"* |
| 4 | the falloff curve | bins by distance to the nearest endothelial cell (`vessel_profile.direct_vessel_distance_um`) |
| 5 | the perivascular reference | the first bin of that curve |
| 6 | the delivery spread | *"a percentile ratio over the resulting field"* |

and closes with *"Read NOTHING on this section as a statement about where its
vessels are."*

### 1.2 Every enforcement point, verified by reading the call sites

| where | what it does | does it BLOCK anything? |
|---|---|---|
| `api.py:1623 _with_caveats` | sets `spatial_spread_valid = False`, appends one `{capability:"spatial_spread", status:"NO", reason:SPOT_SPATIAL_REFUSAL}` row, sets `n_caveats_blocking` | **No.** Disclosure only |
| `api.py:1870 GET /api/sections` | the only call site of `_with_caveats` | No |
| `app/web/js/ui/spec.js:931` | renders `  [N refused]` on the `<option>` and the full reason in its `title` | **No** — the option is disabled only on `!available` or `runnable === false`; its own comment says *"A blocking caveat is not the same as unavailable: the section still runs"* |
| `app/web/js/ui/tissue.js:134,492,697` | restates it in the ingest panel before a path is typed | No |
| `app/web/js/ui/inspect.js:611-620` | a `refused` provenance band on a per-cell inspector | No |
| `tmesim/data/ingest.py:1592-1617` | a **second, near-identical copy** of the same text, for an ingested section's stored capability report | No |
| `tmesim/validation/immune_sink_impact.py:754 require_measured_cell_positions` | **raises `SpotDataError`** — the only hard refusal in the chain | Yes, but see below |

**The only hard refusal never fires in production.** `require_measured_cell_positions`
is reached from `immune_sink_impact._spatial_context`, and a grep over `app/` and
`tmesim/` finds **no** importer of that module outside its own file and comments —
consistent with `CLAIM_INVENTORY` §1c note (f), *"the module ships default-off and no
shipped code path imports it"*. So in the shipping product the refusal is
**entirely a disclosure**: a Visium section runs, the pipeline emits every number it
emits for a Xenium section, and the caveat text is what stands between those numbers
and a reader.

### 1.3 Why the reason is a good reason

It is. `reconstruct/real_visium.py:35` sets `_CELLS_PER_SPOT = 12`;
`_place_cells` (line 67) puts `round(frac * 12)` cells per type at the spot centroid
with `jitter_um = min(27.5, 0.45*dx)` inside a uniform disc. A ~55 µm spot is not a
cell, so a per-cell spatial statistic on those coordinates is a statistic about a
grid of spots wearing a cell's name. Nothing below contradicts that, and
`results/synthesised_positions.md` already measured the consequence for the
*distance* field (p99 vessel distance 5.67× distorted, tissue fraction beyond 150 µm
23× larger). **This document does not reopen that.**

---

## 2. The pre-registered criterion

Written into `experiments/visium_section_level.py`'s module docstring **before the
first run**, and reproduced here verbatim in substance.

A scalar emitted by the shipped pipeline on a Visium section is **RECOVERABLE** iff
it passes all three:

* **C1 — JITTER-INVARIANCE (measured).** Re-run the shipped pipeline with the
  pseudo-cell placement seed changed from the shipped `0` to `1` and `2`, changing
  nothing else. On **every one of the four sections**:
  * `RATIO` rule: `max_s |x_s − x_0| / |x_0| ≤ 0.01`;
  * `FRACTION` rule (quantities in [0,1], where a relative test near zero is unfair):
    `max_s |x_s − x_0| ≤ 0.01` absolute.
* **C2 — MULTIPLICITY-INDEPENDENCE (judged from code, declared in advance).** The
  quantity must be definable without `_CELLS_PER_SPOT = 12`. A cell COUNT is
  `12 × occupancy` and fails no matter how jitter-invariant it is. A statistic over
  the pseudo-cell population as a *fraction* of it is recorded as `WEIGHTED` — a real
  spot-level statistic weighted by the deconvolved fractions, but **not** a per-cell
  statistic — and `WEIGHTED` is treated as a C2 **fail**.
* **C3 — RESOLUTION HONESTY (judged from the definition, declared in advance).** No
  length scale below the ~100 µm Visium spot pitch may enter the definition. Passing
  C3 does not make the resolution limit go away.

**Pre-declared falsifier.** If `mean_vascular_density`, `mean_ecm_density` or
`mean_tissue_auc` moved by more than the C1 threshold under a reseed, then the
delivery source on a Visium section *is* downstream of the synthesised positions,
the refusal's stated mechanism is right as written, and the answer to the whole
audit is "no, correctly scoped". **That falsifier was not triggered** (§3.1).

**Why C1 alone is not enough, stated in advance.** Three seeds is a small sample, and
C1 tests only *"is this a function of the jitter DRAW"*. It cannot see a dependence
on the synthesis that is common to every draw — which is exactly what C2 exists for,
and exactly what caught `mean_oxygen` in §3.2.

---

## 3. What was measured

Configuration: temozolomide, `dx_um = 30` (the Visium loader raises it to the 100 µm
spot pitch), 6 h window, shipped defaults otherwise — the same configuration the
2026-08-23 sweep used, so its rows remain comparable. Seeds 0 / 1 / 2, where 0 is
the shipped behaviour.

**The probe is not vacuous.** On every section the reseed moved **100 %** of
pseudo-cells, median displacement **24.4–24.6 µm**.

**The pipeline is deterministic.** Two runs at the *same* seed on
`visium_10x_parent_gbm` returned `mean_tissue_auc`, `E_p50_uM`, `mean_oxygen`,
`kill_fraction` and `E_p90/E_p10` **bit-identical**, so every deviation below is the
jitter and nothing else.

### 3.1 Structural check, with a positive control

`reconstruct.real_fields.vascular_density_from_cells` — the function that builds a
vessel source *from endothelial cell positions* — was instrumented with a call
counter on the `real_pipeline` binding the reconstructor actually uses.

| section | `vascular_density_from_cells` calls (seeds 0/1/2) | pipeline's own `vascular_source` |
|---|---|---|
| `visium_10x_parent_gbm` | **0 / 0 / 0** | `KDE of REAL ENDOTHELIAL spot fraction` |
| `visium_gbm1` | **0 / 0 / 0** | `KDE of REAL ENDOTHELIAL spot fraction` |
| `visium_gbm2` | **0 / 0 / 0** | `KDE of REAL ENDOTHELIAL spot fraction` |
| `visium_gbm3` | **0 / 0 / 0** | `KDE of REAL ENDOTHELIAL spot fraction` |
| **`xenium_319` (positive control)** | **1 / 1 / 1** | `DISCRETE segmented vessels from REAL ENDOTHELIAL cell positions` |

The control proves the counter works. The Visium source is
`real_visium._field_from_fraction`: the per-spot ENDOTHELIAL fraction scattered at
the **measured** spot centroids and Gaussian-smoothed at a 100 µm bandwidth. No
pseudo-cell enters it. `app/server/cache.py:1296` says the same thing in its own
words when `vascular_mode='kde'` is requested on a Visium section — *"reconstruct_
from_visium builds vascular_density from the measured ENDOTHELIAL spot fraction and
exposes no mode switch"* — and the shipped `vascular_mode='discrete'` branch
re-derives nothing.

**The three arrays the transport solve consumes, compared BITWISE across seeds:**

| array | `10x_parent` | `gbm1` | `gbm2` | `gbm3` |
|---|---|---|---|---|
| `domain.tissue_mask` | identical | identical | identical | identical |
| `fields.vascular_density` | identical | identical | identical | identical |
| `fields.ecm_density` | identical | identical | identical | identical |
| `fields.oxygen` | **differs**, max abs 3.355e-03 | 3.508e-03 | 3.869e-03 | 3.462e-03 |

Oxygen differs because its sink is `consumption × rho_tum` and
`physics/oxygen.py::_tumor_density` bins **tumour-cell positions**. That is the
probe working, not a defect.

### 3.2 THE TWO PLACES THE PRE-REGISTRATION WAS WRONG

Recorded in full because a pre-registration that is quietly edited after the numbers
arrive is worthless. Both errors were in the **permissive** direction — both would
have handed back a quantity that should stay refused — and both were caught by the
data, not by re-reading the code out of virtue.

1. **`tissue_median_over_plasma_peak`.** Pre-registered as C2 = `grid`, on the
   strength of its name and of `jobs.py`'s comment calling it *"THE BOUND"*. Its
   measured C1 deviation came back as `5.8184e-03` — **exactly** `median_tumor_auc`'s,
   to every digit, on all four sections. That is not a coincidence:
   `app/server/jobs.py:1151` computes it as `E_p50_uM / plasma_peak`, and `E_p50_uM`
   is the 50th percentile of **per-tumour-cell** exposure, not of the field over the
   tissue mask. **Reclassified `weighted` → REFUSED.** The field name says "tissue"
   and the quantity is a cell statistic; that is worth fixing independently of this
   audit.
2. **`mean_oxygen`.** Pre-registered as C2 = `grid`, with a note flagging that its
   sink reads tumour-cell density — the note was right and the class was wrong.
   `_tumor_density` is a raw **count** per grid cell divided by `dx²`, so on a Visium
   section it is `12 × (spot tumour fraction) / dx² × ref_area` and scales linearly
   with `_CELLS_PER_SPOT`. **Reclassified `weighted` → REFUSED.** It passed C1 at
   1.5e-04, so it is not jitter-driven; it is multiplicity-driven, which is precisely
   the failure C2 was written to catch.

### 3.3 The enumerated list — what is refused, and why

`RATIO` rule unless marked; `max` is the worst of the four sections.

#### REFUSED — genuinely per-cell, or sub-pitch, or an artefact of the 12

| quantity | C1 max | fails | why |
|---|---:|---|---|
| **`penetration_length_um`** | **1.9416e-01** | C1, C3 | fitted against distance to the nearest endothelial PSEUDO-cell; per section 7.69 % / 5.69 % / 11.75 % / **19.42 %** |
| `E_p10_uM` | **1.0859e-02** | C1, C2 | 10th percentile of per-pseudo-cell exposure |
| `p10_E_over_ic50` | **1.0859e-02** | C1, C2 | same population |
| `E_p90_over_p10` | 6.0373e-03 | C2 | THE HEADLINE. A percentile ratio over the pseudo-cell population |
| `E_p90_uM` | 6.0196e-03 | C2 | same |
| `E_p50_uM` | 5.8184e-03 | C2 | same |
| `median_tumor_auc` | 5.8184e-03 | C2 | same |
| `tissue_median_over_plasma_peak` | 5.8184e-03 | C2 | **reclassified in §3.2** — it is `E_p50 / plasma_peak` |
| `low_o2_over_high_o2_exposure` | 5.1099e-03 | C2 | an exposure ratio between two pseudo-cell bands |
| `therapeutic_coverage_1x` (frac) | 2.4394e-03 | C2 | fraction of pseudo-tumour-cells above 1× IC50 |
| `n_viable_reach` | 1.5347e-03 | C2 | a pseudo-cell count |
| `mean_tumor_auc` | 1.2691e-03 | C2 | AUC sampled at pseudo-cell positions |
| `reach_headline_F` (frac) | 7.0426e-04 | C2, C3 | binned on the same sub-pitch vessel-distance transform |
| `viable_tumor_fraction` (frac) | 5.5426e-04 | C2 | state labels banded at pseudo-cell positions |
| `kill_fraction` (frac) | 5.5426e-04 | C2 | **never a defensible output at all**; recorded, never read |
| `mean_oxygen` | 1.5155e-04 | C2 | **reclassified in §3.2** — sink ∝ `_CELLS_PER_SPOT` |
| `n_cells`, `n_tumor` | 0.0 | C2 | `12 × occupancy` by construction; exactly jitter-invariant and exactly artefactual |

The last row is the clearest illustration of why C1 alone would have been a bad
criterion: `n_cells` is **perfectly** reproducible across seeds and is a pure
artefact of the number 12.

#### RECOVERABLE — passes C1, C2 and C3

| quantity | C1 max | note |
|---|---:|---|
| `mean_vascular_density` | **0.0** | field bitwise identical |
| `mean_ecm_density` | **0.0** | field bitwise identical |
| `penetration_length_analytic_um` | **0.0** | closed form; no field solve |
| `mean_tissue_auc` | 3.0048e-16 | 1 ULP; see the conservation caveat in §4.2 |
| `tissue_coverage_fraction` (frac) | 0.0 | **vacuous on this platform — see §4.1** |
| `plasma_peak_uM`, `plasma_peak_over_ic50` | 0.0 | PK inputs. No tissue content whatsoever; listed for completeness, not as a result |

### 3.4 The delivery map — the mechanism the refusal does not name

Measured on `visium_gbm1`, seeds 0 vs 1, `simulate_response` driven directly so the
convection switch could be flipped:

| | convection **ON** (shipped default) | convection **OFF** |
|---|---|---|
| AUC field bitwise identical across seeds | **No** | **Yes** |
| tissue pixels that change | **100.00 %** | 0.00 % |
| per-pixel relative change, p50 / p90 / p99 / max | 1.73e-03 / 8.53e-03 / 1.76e-02 / **2.82e-02** | 0 / 0 / 0 / 0 |
| `mean_tissue_auc` relative change | 4.04e-15 | 0.00e+00 |
| IFP field bitwise identical | No, max abs diff **2.625** | No, max abs diff 2.625 |

`PDEConfig.convection` is **True** by default. `tmesim/physics/pressure.py:103`
rasterises tumour cells into a tumour-occupancy field *"(used to locate the leaky,
poorly-drained tumor region)"*; that field moves under a reseed, the IFP moves with
it, and the Darcy velocity advects the drug. **So the exposure map on a Visium
section is jitter-dependent — through TUMOUR cells and the pressure solve, not
through the endothelial source term the refusal blames.**

The two facts sit together without contradiction: the map's *source* carries no
jitter (§3.1, bitwise-identical `vascular_density`), and the map's *advection* does.
The refusal's conclusion is right; its stated causal chain is not the one in the code.

---

## 4. The recovered quantities, on the four sections, with the caveats that apply

**Everything in this table is a ~100 µm-resolution quantity.** Nothing below the
Visium spot pitch was measured on these sections, and surviving a jitter test does
not change that. Every cell-type fraction feeding these numbers is a **deconvolution
product** — a model output scored against the Darmanis reference, not a
segmentation call.

| | `10x_parent` | `gbm1` | `gbm2` | `gbm3` |
|---|---:|---:|---:|---:|
| spots | 3,468 | 4,764 | 2,891 | 1,775 |
| grid at dx = 100 µm | 57 × 50 | 55 × 58 | 54 × 55 | 37 × 45 |
| **grid cells with a spot under them** | **78.1 %** | **93.4 %** | **61.5 %** | **68.7 %** |
| `tissue_coverage_fraction` as reported | 1.0 | 1.0 | 1.0 | 1.0 |
| `mean_vascular_density` | 0.13214017 | 0.11044153 | 0.30246834 | 0.32503647 |
| `mean_ecm_density` | 0.02985667 | 0.01516004 | **0.0** | 0.02541225 |
| `mean_tissue_auc` (µM·s) | 2.0160792e+07 | 1.6850204e+07 | 4.6147978e+07 | 4.9591226e+07 |
| `penetration_length_analytic_um` | 288.448 | 289.740 | 291.067 | 288.840 |
| deconvolved mean ENDOTHELIAL fraction | 0.0977 | 0.0579 | 0.2950 | 0.2372 |
| deconvolved mean STROMAL fraction | 0.0104 | 0.0033 | **0.00000** | 0.0008 |

### 4.1 `tissue_coverage_fraction` is vacuous here, and that matters for the row above it

`core/domain.py:41` defaults `tissue_mask` to **all-True**, and
`reconstruct_from_visium` never sets one. So `tissue_coverage_fraction = 1.0` on
every Visium section is *the absence of a mask*, not a property of the tissue — and
`_field_from_fraction`'s docstring claim that it *"confines to the tissue mask"* is
a no-op there. Measured directly: only **61.5–93.4 %** of each domain rectangle has a
spot under it, so **6.6–38.5 % of the area every "tissue" mean is taken over has no
measurement beneath it at all**, only KDE bleed from neighbouring spots. This is not
a reason to withdraw `mean_tissue_auc`; it is the reason the word "tissue" in its
name should be read as "bounding rectangle".

### 4.2 `mean_tissue_auc` survives for a reason that limits what it licenses

The AUC **field** is not jitter-invariant (§3.4) — 100 % of its pixels move. Its
tissue-mask **mean** is invariant to 4e-15 because advection redistributes drug
rather than adding or removing it, so the mean is conserved while the map under it
changes. **`mean_tissue_auc` therefore licenses nothing whatever about the map.** A
scalar passing C1 while the field it summarises fails is a real hazard in this
criterion and is recorded as such.

Separately: absolute exposure is already in `CLAIM_INVENTORY` §1d's failure column
(breast Kp 12.8× high for 5-FU, 10.5× **low** for methotrexate — opposite directions,
so no constant absorbs both). Replicating an absolute number that is already known to
be wrong on two independent anchors adds generality to a failure, not support to a
claim.

### 4.3 `mean_vascular_density` and `mean_ecm_density` are not densities

`real_visium._field_from_fraction` ends with `_normalize01`, a per-section min-max
rescale to [0,1]. So both numbers are **means of a field whose maximum is 1.0 by
construction on that same section** — shape statistics, not physical densities, and
their cross-section comparison is a comparison of mean-over-max ratios. Two further
facts bound them:

* the deconvolved ENDOTHELIAL fraction is **0.2950** on `gbm2` and **0.2372** on
  `gbm3`, against a measured human vascular areal fraction of **4.70 ± 2.38 %**
  (Marien 2016, PMID 27583442, as cited in `CLAIM_INVENTORY` §1d). A 24–30 %
  endothelial call is a deconvolution artefact by a factor of ~5–6;
* `gbm2`'s STROMAL fraction is **exactly 0.00000 at every one of its 2,891 spots**
  (sd 0, min 0, max 0). The ECM field is therefore identically zero, the
  `D_eff = D·(1 − 0.6·ecm)` tortuosity term is **inert** on that section, and
  `mean_ecm_density = 0.0` is a degenerate deconvolution, not a stroma-free tumour.
  This is visible downstream: `gbm2`'s `mean_tissue_auc` is the second-highest of
  the four and its (refused) fitted length is 412.9 µm.

### 4.4 `penetration_length_analytic_um` carries almost no section information

288.4 / 289.7 / 291.1 / 288.8 µm — a **0.9 %** spread across four sections, because
it is a closed form in `mean_ecm_density` and the compiled drug, and
`mean_ecm_density` is between 0.0 and 0.030 everywhere. It is exactly reproducible
and nearly constant, which are two different things and only the first is a virtue.

---

## 5. Verdict: is the refusal over-applied?

**On its list of quantities: no.** All six things `SPOT_SPATIAL_REFUSAL` names are
either measurably jitter-dependent (`penetration_length_um` at up to 19.4 %, the
exposure map at 100 % of pixels), sub-pitch by definition (the falloff, the
perivascular reference), or refused elsewhere and never computed in production (the
attenuation spread, the structure ratio). **Not one of them should come back.**

**On its stated reason: yes, in one specific and checkable way.** The clause *"its
SOURCE is a density over ENDOTHELIAL CELL POSITIONS"* is false for the Visium path,
and it is the **only** justification the text offers for refusing the exposure map.
A reviewer who checks it finds `vascular_density_from_cells` uncalled, a provenance
string that says "spot fraction", and a bitwise-identical vascular field — and is
then entitled to conclude the map survives, which it does not. **A refusal resting
on a false premise is a weaker refusal than one resting on a true one**, and this is
the project's own recurring defect class: the engine is sound, the document is wrong
(`results/CLAIM_INVENTORY.md`'s claim-surface-drift finding; `docs/VERIFICATION_
METHODOLOGY.md` §10b, *"a consistently-applied wrong convention validates itself"*).

**On section-level quantities: marginally, and the margin is small.** Five
quantities pass; two of the five (`tissue_coverage_fraction`, and the PK inputs) are
vacuous; two more (`mean_vascular_density`, `mean_ecm_density`) are normalised
deconvolution shape statistics with a ~5–6× artefact in one of their inputs and an
exactly-zero degenerate case in the other; one (`mean_tissue_auc`) is an absolute
exposure already failing two external anchors, and it licenses nothing about the map
it averages. R12's list of things the refusal *"does not apply to"* — *"the tissue
mask, the vascular density field, the ECM field, grid convergence, the invariance
measurements"* — is **correct as a statement about the refusal** and **optimistic as
a statement about what those quantities are worth on these four sections**.

**The strongest honest sentence available from this pass** is not about recovery at
all. It is: *"On four Visium GBM sections, re-running the shipped pipeline with
nothing changed but the seed that places pseudo-cells inside their spots moves the
fitted penetration length by up to 19.4 % and every pixel of the exposure map, while
leaving the vascular field, the ECM field and the tissue-mask mean exposure
bit-for-bit unchanged."* That is a measurement of what the spot refusal is actually
protecting against, and it is the first one this repository has that separates the
contaminated quantities from the clean ones **by experiment rather than by reading
the code**.

---

## 6. Proposed change to `app/server/api.py` — NOT APPLIED

`app/server/api.py` is the security boundary and another workstream is in it. The
following is a proposal for the file's owner. It is a change to **one string
constant** and to **nothing else**: no control flow, no new field, no removed field,
no change to `_SPOT_SOURCES`, `spatial_spread_valid`, `n_caveats_blocking`, or which
sections carry the caveat. Every currently-refused quantity stays refused.

**Why it is safe.** `tests/test_cross_section.py:249,261` compare the emitted reason
to the constant by identity and equality (`refusal is SPOT_SPATIAL_REFUSAL`,
`== SPOT_SPATIAL_REFUSAL`), never to a literal, so they pass unchanged.
`tests/test_ingest.py:842` asserts `"12" in reason` — that is the *ingest* copy, not
this one, and the replacement text keeps the "12 pseudo-cells" sentence regardless.
No JS reads the string's content; `spec.js`, `tissue.js` and `inspect.js` render it
verbatim or name the constant.

Replace lines **1604–1613** — from `"HERE TOO, AND FOR THE SAME REASON. The exposure
field is solved on the "` through `"where its vessels are. (An earlier version of
this sentence, written on "` — with the measured mechanism. The `-` block below is
verbatim from the **working-tree** `app/server/api.py` at the time of writing, and
the surrounding context lines are unchanged. (That file carries 132 lines of another
workstream's in-flight edits; `git diff HEAD -- app/server/api.py` was checked and
**none of them touch `SPOT_SPATIAL_REFUSAL`, `_SPOT_SOURCES` or `_with_caveats`** —
but the line numbers above are working-tree numbers, so match on the text, not on
the offsets.)

```diff
     "immune-neighbourhood structure ratio would measure that jitter, not the "
     "tissue (results/immune_sink_impact.md §5.6). THE DELIVERY MAP IS INVALID "
-    "HERE TOO, AND FOR THE SAME REASON. The exposure field is solved on the "
-    "grid, but its SOURCE is a density over ENDOTHELIAL CELL POSITIONS, and on "
-    "a spot section those positions are synthesised exactly like every other "
-    "cell's. The falloff curve bins by distance to the nearest endothelial "
-    "cell (validation.vessel_profile.direct_vessel_distance_um), the "
-    "perivascular reference is the first bin of that curve, and the delivery "
-    "spread is a percentile ratio over the resulting field - so all three "
-    "measure the 12-pseudo-cells-per-spot synthesis, not this tissue's "
-    "vascular architecture. Read NOTHING on this section as a statement about "
-    "where its vessels are. (An earlier version of this sentence, written on "
+    "HERE TOO, BUT NOT BY THE SAME ROUTE, AND THE ROUTE MATTERS. The delivery "
+    "SOURCE on a Visium section is NOT built from cell positions: "
+    "reconstruct.real_visium._field_from_fraction scatters the per-spot "
+    "ENDOTHELIAL fraction at the MEASURED spot centroids and smooths it at a "
+    "100 um bandwidth, and reconstruct.real_fields."
+    "vascular_density_from_cells is never called on this platform - measured, "
+    "0 calls on four Visium sections against 1 on xenium_319 "
+    "(results/visium_section_level_recovery.md §3.1). What DOES carry the "
+    "synthesis into the map is CONVECTION: physics.pressure rasterises TUMOUR "
+    "cell positions to locate the leaky, poorly-drained region, so the IFP and "
+    "its Darcy velocity move with the jitter and advect the drug. With "
+    "PDEConfig.convection at its shipped default of True, a reseed changes "
+    "100% of tissue pixels of the exposure field (median 0.17%, p99 1.76%, max "
+    "2.82%); with convection off the field is bitwise identical (ibid. §3.4). "
+    "The falloff curve bins by distance to the nearest endothelial cell "
+    "(validation.vessel_profile.direct_vessel_distance_um), the perivascular "
+    "reference is the first bin of that curve, and the delivery spread is a "
+    "percentile ratio over per-TUMOUR-CELL exposure sampled at those same "
+    "synthesised positions - so all three measure the "
+    "12-pseudo-cells-per-spot synthesis, not this tissue's vascular "
+    "architecture. Measured under a pure reseed: the fitted penetration length "
+    "moves by up to 19.4%, E_p10 by up to 1.09%, the delivery spread by up to "
+    "0.60% (ibid. §3.3). Read NOTHING on this section as a statement about "
+    "where its vessels are. WHAT THIS DOES NOT REFUSE, because it does not "
+    "read a cell position at all: the vascular and ECM FIELDS and the "
+    "tissue-mask MEAN of the exposure field, all three bitwise or 1-ULP "
+    "invariant under a reseed. They are ~100 um-resolution deconvolution "
+    "products, the reported tissue_coverage_fraction of 1.0 is the ABSENCE of "
+    "a tissue mask rather than full coverage, and the mean of a field says "
+    "nothing about the map under it (ibid. §4). (An earlier version of this "
+    "sentence, written on "
     "2026-08-21, claimed the exposure map and its falloff DID survive this "
     "refusal because they are solved on the grid. That was wrong: it confused "
     "where the PDE is discretised with where its source term comes from. The "
```

**Three notes for whoever applies it.**

1. The surviving parenthetical *"it confused where the PDE is discretised with where
   its source term comes from"* is now itself imprecise — the source term is exactly
   the part that turned out to be clean. It is left in place because it is a dated
   historical correction and this document's rule is that dated records are not
   rewritten; a second parenthetical would be the honest fix if one is wanted.
2. **`tmesim/data/ingest.py:1592-1617` carries a second, near-identical copy of the
   same false sentence** and needs the same correction. `api.py:1596`'s comment
   claims the refusal is *"Stated once here so the API, the picker and the report
   cannot drift apart"* — it is in fact stated twice in Python, and that is the drift
   risk the comment was written to prevent.
3. If the owner would rather not grow a long string, the minimum honest change is
   the two clauses in **bold** in §5: delete *"its SOURCE is a density over
   ENDOTHELIAL CELL POSITIONS"* and name convection instead. Everything else in the
   diff is elaboration.

---

## 7. What a recovery actually buys — the scout's claim, checked

R12 states that declining these sections *"is currently costing four sections of
replication for free"*, and the commission asks whether that is replication on the
project's own #1 priority.

**It is not, and the arithmetic is short.** `CLAIM_INVENTORY` §6 item 1 is *"A second
section, same tumour type. Almost every defensible claim above is n = 1 in tissue."*
The defensible claims are C1 (common-mode insulation, breast Xenium, a 3-D:2-D kill
ratio), C4 (a structural identity, which any input satisfies and which therefore
gains nothing from real tissue), K1 (the 3-D:2-D kill panel) and S1/S2 (the spatial
spread — refused on spot data, correctly). **None of the five recovered quantities
appears in the backing of any of them.** The recovered set cannot replicate a kill
panel, because kill is refused; it cannot replicate the spatial spread, because that
is the statistic the refusal is about; and it cannot replicate C1, because that is a
ratio of kill ratios between configurations.

**What the four sections DO buy, honestly stated:** a four-section demonstration that
the vascular field, the ECM field and the mean tissue exposure are exactly
reproducible under the synthesis — i.e. a **measurement of the refusal's own
boundary**, which is what §5 recommends reporting — plus the two degeneracies in §4.3
(`gbm2`'s identically-zero stroma, the 5–6× endothelial over-call) that are findings
about the **ingest and deconvolution** layer and were not previously written down.

**The thing that really is free, and is not this.** All seven sections were run on
**2026-08-23** by `experiments/multisection_generalisation.py`, and
`results/multisection_generalisation.json` (129 KB, 9 runs) has held the payloads
since. `experiments/multisection_report.py` — 838 lines, tracked, whose module
docstring says it *"writes the report… every number in the document is interpolated
from that file"* — was written to render `results/multisection_generalisation.md`.
**That file has never existed**: not in the working tree, and `git log --all --
results/multisection_generalisation.md` returns nothing. Its findings survive only
in commit `bf15977`'s message (*"the dominant limiter takes THREE different values
across seven sections… therapeutic coverage spans the entire range 0.0 to 1.0"*),
and `CLAIM_INVENTORY` §6 still ranks "a second section" first without citing it.
Rendering that report is a one-command job against data that is already on disk, and
it is a larger and cheaper win than anything in this document. **It is outside this
pass's write territory and is left for its owner.**

---

## 8. What would falsify this document

* **A Visium section whose `vascular_density` array is NOT bitwise identical across
  placement seeds.** That would mean §3.1's structural finding is section-specific
  rather than a property of the code path, and the refusal's original wording would
  be right somewhere.
* **A reseed sweep with more than three seeds in which `mean_vascular_density`,
  `mean_ecm_density` or `mean_tissue_auc` moves.** Three seeds is a small sample and
  §2 says so; the bitwise-identical *fields* are the stronger evidence, and they are
  what would have to fall.
* **A demonstration that `physics/pressure.py`'s tumour-occupancy field is not the
  route in §3.4** — e.g. the AUC field still moving under a reseed with
  `convection=False` on a different section. Measured here on `visium_gbm1` only.
* **A shipped code path that imports `validation.immune_sink_impact`.** §1.2's claim
  that the hard refusal never fires in production rests on a grep plus
  `CLAIM_INVENTORY` §1c(f); one importer would change the character of the whole
  refusal from a disclosure to a gate.
* **A section-level claim that the five recovered quantities do support and that
  §7 missed.** §7 is an argument from the claim table as it stands on 2026-09-20; a
  claim added later could change it.

---

`calibrate_alpha()` = **22.223477711172563**, asserted at the start of the run, after
each of the four sections, and at the end. Driver:
`experiments/visium_section_level.py`.
