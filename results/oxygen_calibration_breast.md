# Oxygen calibration on the breast reconstruction: does the same PDE look more or less credible on a second tumor type?

**Status: measurement + literature comparison. Zero core-engine changes.** This extends `results/oxygen_calibration_landscape.md` (glioma) and `results/cell_state_impact.md` (state-banding rigor) to the breast flagship (`results/breast_cancer_flagship.md`). It runs the existing, unmodified oxygen reaction-diffusion solve (`tmesim.physics.oxygen.solve_oxygen`, invoked internally by `reconstruct_real`) and the existing, unmodified cell-state banding (`tmesim.cells.agents.update_cell_states`) on the real breast Xenium reconstruction, converts the result to mmHg via the newly added `tmesim.physics.oxygen.to_mmhg()` / `REFERENCE_PO2_MMHG` anchor, and checks the answer against real, PMID-cited breast-tumor pO2 measurements (the Vaupel group's Eppendorf/computerized-pO2-histography electrode literature). No file under `tmesim/` was edited to produce this report — see §5.

## Headline

| | |
|---|---|
| **Reconstruction** | `reconstruct_real(source="xenium", sample="data_cache/xenium_breast_janesick/rep1", dx_um=12.0)` — 167,780 real cells (89,785 tumor), grid 458×628, tissue coverage 93.1% |
| **Oxygen field, whole tissue (masked mean, matching `breast_cancer_flagship.md`)** | mean O = 0.3457 |
| **Oxygen field, whole grid array** | mean 0.3502, median 0.3454, std 0.0891, range [0.1338, 0.5430] |
| **Oxygen at TUMOR cells (the population every downstream kill number is computed over)** | mean 0.2815, median 0.2788, std 0.0646, range [0.1338, 0.5422] |
| **Same, in mmHg via `to_mmhg()` (anchor `REFERENCE_PO2_MMHG = 35.0` mmHg)** | mean **9.85 mmHg**, median **9.76 mmHg**, std 2.26 mmHg, range **[4.68, 18.98] mmHg** |
| **Tumor-cell state fractions (`update_cell_states`, unmodified)** | **HYPOXIC 11.22%**, QUIESCENT 84.16%, VIABLE 4.63%, DEAD 0.00% |
| **Fraction of tumor cells below the literature's classic severe-hypoxia cutoff (2.5 mmHg)** | **0.0000%** (n = 0 of 89,785) |
| **Sink form / consumption observed at run time** | `sink_form = "linear"` (the original `consumption * rho_tum * O` form), `consumption = 0.02`, `supply = 0.05`, `diffusivity = 2000 µm²/s` — all `OxygenConfig` defaults, unchanged from what `results/oxygen_calibration_landscape.md` and `oxygen.py`'s own docstring describe for GBM. A concurrent task had already added a `sink_form` field and Michaelis–Menten option to `OxygenConfig` by the time this ran (see §5), but its default is `"linear"` and that is what this run took. |
| **Is breast's headline more plausible than GBM's 0.0% HYPOXIC?** | **Only superficially.** At the model's OWN threshold (7 mmHg, inherited unchanged from the dimensionless `t=0.2` band) breast reads a non-trivial 11.2% HYPOXIC. At the literature's actual severe-hypoxia cutoff (2.5 mmHg) breast is exactly as implausible as GBM: **0.0%**. See §4. |

---

## 1. The solved O2 field, and its mmHg conversion

Ran unmodified: `reconstruct_real` internally calls `physics.oxygen.solve_oxygen` with the default `OxygenConfig()` (`oxygen_mode='pde'`, the shipping default), storing the result in `micro.fields.oxygen`. Nothing about that call path was touched for this task.

| quantity | dimensionless O | mmHg (`to_mmhg`) |
|---|---:|---:|
| whole-grid mean | 0.3502 | 12.26 |
| whole-grid median | 0.3454 | 12.09 |
| whole-grid std | 0.0891 | 3.12 |
| whole-grid min / max | 0.1338 / 0.5430 | 4.68 / 19.00 |
| tumor-cell mean | 0.2815 | 9.85 |
| tumor-cell median | 0.2788 | 9.76 |
| tumor-cell std | 0.0646 | 2.26 |
| tumor-cell min / max | 0.1338 / 0.5422 | 4.68 / 18.98 |

Tumor-cell mmHg percentiles (the population every kill number in `results/breast_cancer_flagship.md` is computed over):

| percentile | mmHg | O |
|---|---:|---:|
| p1 | 5.37 | 0.1536 |
| p5 | 6.14 | 0.1754 |
| p10 | 6.88 | 0.1965 |
| p25 | 8.41 | 0.2404 |
| p50 (median) | 9.76 | 0.2788 |
| p75 | 11.22 | 0.3207 |
| p90 | 12.87 | 0.3678 |
| p95 | 13.91 | 0.3974 |
| p99 | 15.84 | 0.4525 |

**Note on the anchor.** `REFERENCE_PO2_MMHG = 35.0` mmHg is defined in `tmesim/physics/oxygen.py` from Erecinska & Silver 2001 (*Respir Physiol* 2001;128(3):263-276, PMID [11718758](https://pubmed.ncbi.nlm.nih.gov/11718758/)) — **normal, non-tumour-bearing human BRAIN tissue** pO2. It was chosen for GBM because `solve_oxygen`'s source term (`supply*vasc*(1-O)`) is a tissue-compartment quantity, and that reasoning is tissue-agnostic (it is about which physiological compartment `O=1` represents, not which organ). It was **not re-derived for breast tissue** — applying it here reuses a brain-tissue anchor on breast tissue. This matters quantitatively: Vaupel et al. 1991 (§3 below) measured a normal-breast median pO2 of 65 mmHg, roughly double the 35 mmHg brain anchor, so an anchor genuinely appropriate to breast physiology would very likely stretch this report's mmHg numbers upward across the board. **Every mmHg figure in this report should be read as "on the brain-tissue anchor, applied to breast tissue," not as a breast-specific calibration** — this is an honest limitation carried over unmodified from `oxygen.py`, not something this task introduced or resolved.

---

## 2. Cell-state fractions on breast tumor cells

`cells.agents.update_cell_states(micro, oxygen, cfg)` was called directly (unmodified function, `OxygenConfig()` defaults) — the same call `sim.engine.simulate_response` makes internally and the same one `tmesim/validation/oxygen_benchmark.py::_kill_from_oxygen` and `cell_state_impact.py` use for their own re-banding arms. `reconstruct_real`'s default `state_mode='oxygen'` does **not** call this itself (by design — see `reconstruction_report()['cell_state']`), so this is the same "install it yourself, on the frozen reconstruction" pattern those two GBM validation scripts already use.

| state | n (tumor cells) | fraction |
|---|---:|---:|
| HYPOXIC (`O < 0.2` = < 7.00 mmHg) | 10,071 | **11.2168%** |
| QUIESCENT (`0.2 ≤ O < 0.4` = [7.00, 14.00) mmHg) | 75,560 | **84.1566%** |
| VIABLE (`O ≥ 0.4` = ≥ 14.00 mmHg) | 4,154 | **4.6266%** |
| DEAD | 0 | 0.0000% |

n_tumor = 89,785. Band edges are set by `OxygenConfig.hypoxia_threshold = 0.2` (HYPOXIC below it, QUIESCENT up to `2×`, VIABLE above), unmodified from the GBM default, translated to mmHg by the same `to_mmhg()` anchor as §1.

**Contrast with GBM.** On the GBM flagship section, at the current default `consumption`, the PDE-route HYPOXIC tumor fraction is **0.0%** at every grid spacing (`tmesim/physics/oxygen.py`'s own docstring; corroborated by `results/cell_state_impact.md` §3.1, `results/oxygen_demotion.md` §3.1). Breast's 11.2% is not zero. Whether that makes it *more credible* is addressed directly in §4 — the short answer is that it depends entirely on which pO2 threshold is used to call a cell "hypoxic," and at the literature's own preferred threshold breast is equally at 0.0%.

---

## 3. Real breast-cancer tumor pO2 literature — PubMed search, PMID-verified

Searched PubMed directly (`breast cancer tumor oxygen tension Eppendorf polarographic`, `Vaupel breast cancer oxygenation pO2 tissue`, `breast cancer hypoxic fraction pO2 histography prognosis`). The Vaupel/Höckel/Mainz group's computerized-pO2-histography (Eppendorf-electrode-derived) program is exactly the literature the task named, and it is substantial and breast-specific — stronger, in fact, than what was available for glioma (`oxygen_calibration_landscape.md` found only 4 glioma-specific electrode studies; the Vaupel breast program alone returns 7+ directly relevant papers). Every number below is quoted from the article abstract/metadata returned by the PubMed tool; none was independently re-derived from a full-text table, so treat exact decimals with the same caution `oxygen_calibration_landscape.md` applied to its own UNVERIFIED-flagged numbers.

| Study | Tissue | Technique | Key reported pO2 |
|---|---|---|---|
| **Vaupel P, Schlenger K, Knoop C, Höckel M. "Oxygenation of human tumors: evaluation of tissue oxygen distribution in breast cancers by computerized O2 tension measurements." *Cancer Res* 1991;51(12):3316-3322.** PMID [2040005](https://pubmed.ncbi.nlm.nih.gov/2040005/) | Breast cancer (pathological stages T1–T4), fibrocystic disease, normal breast — awake patients, computerized electrode-motion pO2 histography | Direct polarographic pO2, systematic clinical series | **Normal breast median pO2 = 65 mmHg** (no readings below 12.5 mmHg). **Fibrocystic disease median = 67 mmHg** (no readings ≤ 5 mmHg). **Breast cancer, pooled T1–T4, median pO2 = 30 mmHg.** 6 of 15 breast cancers had tissue areas with pO2 between 0–2.5 mmHg (severe hypoxia/near-anoxia). No correlation found between oxygenation and stage/grade; tumor-to-tumor variability exceeded intra-tumor heterogeneity. |
| **Vaupel P, Mayer A, Briest S, Höckel M. "Oxygenation gain factor: a novel parameter characterizing the association between hemoglobin level and the oxygenation status of breast cancers." *Cancer Res* 2003;63(22):7634-7637.** PMID [14633681](https://pubmed.ncbi.nlm.nih.gov/14633681/) | Breast cancer patients, stratified by blood hemoglobin (cHb) | Same computerized pO2 histography | **Median breast-cancer pO2 = 3–15 mmHg** across the measured cHb range (8.5–14.7 g/dL) — i.e., breast-tumor median pO2 is anemia/perfusion-dependent, not a single number. **Normal tissue medians in the same patients: normal breast 52 mmHg, subcutis 51 mmHg, skeletal muscle 37 mmHg** (roughly constant, independent of cHb over the non-anemic range) — note this normal-breast figure differs from the 65 mmHg reported in the 1991 paper above, an inter-study spread of comparable size to the interlaboratory pO2 variability `oxygen_calibration_landscape.md` flagged for glioma (Nozue et al. 1997, PMID 9290690). |
| **Vaupel P, Schlenger KH, Hoeckel M, Okunieff P. "Oxygenation of mammary tumors: from isotransplanted rodent tumors to primary malignancies in patients." *Adv Exp Med Biol* 1992;316:361-371.** PMID [1288097](https://pubmed.ncbi.nlm.nih.gov/1288097/) | Rodent mammary tumor isografts, human breast-cancer xenografts, primary human breast cancer | Comparative pO2 histography across the model spectrum | Human primary lesions span the full range from "low-flow" isotransplants/xenografts (hypoxic/anoxic) to "high-flow" xenografts or normal tissue (well oxygenated) — perfusion efficiency, not stage/grade, is the dominant modulator of tumor oxygenation. |
| **Vaupel P, Höckel M, Mayer A. "Detection and characterization of tumor hypoxia using pO2 histography." *Antioxid Redox Signal* 2007;9(8):1221-1235.** PMID [17536958](https://pubmed.ncbi.nlm.nih.gov/17536958/), DOI [10.1089/ars.2007.1628](https://doi.org/10.1089/ars.2007.1628) | Meta-compilation of 125 clinical pO2-histography studies across **cervix, head-and-neck, and breast cancer** (not breast-exclusive — flagged below) | Review/meta-summary of the same electrode technique | **Overall median pO2 across the three tumor types = 10 mmHg; overall hypoxic fraction (pO2 ≤ 2.5 mmHg) ≈ 25%.** Tumor oxygenation is independent of clinical stage/histology/grade; metastases do not differ substantially from primaries. |
| **Vaupel P, Mayer A, Höckel M. "Impact of hemoglobin levels on tumor oxygenation: the higher, the better?" *Strahlenther Onkol* 2006;182(2):63-71.** PMID [16447012](https://pubmed.ncbi.nlm.nih.gov/16447012/), DOI [10.1007/s00066-006-1543-7](https://doi.org/10.1007/s00066-006-1543-7) | Review across tumor types incl. breast, cervix, head-and-neck | Review of pO2-vs-hemoglobin literature | Maximum tumor oxygenation (breast included) occurs at cHb ≈ 12–14 g/dL (women); both anemia and supraphysiological hemoglobin worsen oxygenation. Not a breast-exclusive pO2 number, cited for the anemia-dependence context relevant to §3's spread. |
| **Kallinowski F, Zander R, Hoeckel M, Vaupel P. "Tumor tissue oxygenation as evaluated by computerized-pO2-histography." *Int J Radiat Oncol Biol Phys* 1990;19(4):953-961.** PMID [2211264](https://pubmed.ncbi.nlm.nih.gov/2211264/), DOI [10.1016/0360-3016(90)90018-f](https://doi.org/10.1016/0360-3016(90)90018-f) | Isografted rodent tumors, xenografted human tumors (incl. breast), first patient measurements | Methods/validation paper for the pO2-histography system itself | Hypoxic and anoxic tissue areas frequently found in both isografted and xenografted tumors; extent depends on cell line, tumor size, implantation site, vascularity, perfusion. First patient measurements showed marked inter- and intra-tumor heterogeneity. |
| **Groebe K, Vaupel P. "Evaluation of oxygen diffusion distances in human breast cancer xenografts using tumor-specific in vivo data." *Int J Radiat Oncol Biol Phys* 1988;15(3):691-697.** PMID [3417489](https://pubmed.ncbi.nlm.nih.gov/3417489/), DOI [10.1016/0360-3016(88)90313-6](https://doi.org/10.1016/0360-3016(88)90313-6) | Human breast cancer xenografts | Diffusion-model computation from tumor-specific in vivo microvascular data | Radioresistance can develop at intercapillary distances **> 100 µm**; if intercapillary distance **> 140 µm**, hypoxia is present even at the arterial end of a microvessel. This is a breast-specific analogue of the Thomlinson & Gray 1955 distance (PMID [13304213](https://pubmed.ncbi.nlm.nih.gov/13304213/)) that `oxygen_demotion.md`'s vessel-distance surrogate already uses for GBM. |
| **Mueller-Klieser W, Schlenger KH, Walenta S, Gross M, Karbach U, Hoeckel M, Vaupel P. "Pathophysiological approaches to identifying tumor hypoxia in patients." *Radiother Oncol* 1991;20(Suppl 1):21-28.** PMID [2020765](https://pubmed.ncbi.nlm.nih.gov/2020765/), DOI [10.1016/0167-8140(91)90182-g](https://doi.org/10.1016/0167-8140(91)90182-g) | Cervix and breast cancer patients | pO2 histography, cryospectrophotometry, bioluminescence ATP imaging (multi-modal) | pO2 values consistently lower in tumor than in adjacent normal tissue; great inter- and intra-tumoral heterogeneity even among same-stage/grade tumors. |

**What was searched for but NOT found (an honest gap, exactly parallel to GBM's §2 in `oxygen_calibration_landscape.md`).** A dedicated PubMed search for a breast-tissue-level oxygen consumption rate (OCR) in diffusion-model-comparable units — the breast analogue of Watabe et al. 2017's rat-glioma ¹⁵O-PET TMRO2 (PMID 28210996) that `oxygen.py`'s docstring now cites for GBM — returned **zero results** (`breast cancer tumor metabolic rate of oxygen PET 15O`). This gap is real and should not be papered over: neither tumor type in this codebase has a directly usable, tissue-level, absolute OCR number to calibrate `OxygenConfig.consumption` against. `OxygenConfig.consumption = 0.02` was run **unmodified and un-recalibrated** for this breast result, carrying forward every uncertainty `oxygen.py`'s own docstring already documents for GBM (the Watabe-vs-Diepart/Milotti bracket spans roughly two orders of magnitude around this same default).

---

## 4. Is breast's non-zero HYPOXIC fraction more credible than GBM's 0.0%? — an honest, threshold-sensitive answer

**Not simply "yes."** The headline 11.2% HYPOXIC fraction (§2) is non-zero and therefore looks, at first glance, less obviously wrong than GBM's flat 0.0% — a tumor type is not defined by necrosis the way GBM is, so a non-trivial hypoxic fraction is a priori more plausible for breast than for GBM even before comparing to literature. But that comparison is being made against the model's **own** threshold (7 mmHg), which was inherited unchanged from the GBM dimensionless band (`t=0.2`) and never independently justified for breast tissue. Re-running the same tumor-cell mmHg distribution against a ladder of thresholds actually used in the literature above gives a different picture:

| mmHg cutoff | source of the cutoff | fraction of breast tumor cells below it (this model) |
|---:|---|---:|
| **2.5 mmHg** | Vaupel 2007's own "hypoxic fraction" definition (PMID 17536958) — the SAME severe-hypoxia cutoff `oxygen_calibration_landscape.md` used for GBM (via Collingridge 1999) | **0.0000%** (0 of 89,785 tumor cells) |
| 5.0 mmHg | roughly the lower end of the 1991 6/15-tumor "0–2.5 mmHg" hypoxic-area band, generously widened | 0.2952% (265 cells) |
| **7.0 mmHg** | **the model's own HYPOXIC threshold** (`t=0.2` translated via `to_mmhg`) | 11.2168% (10,071 cells) — the headline number in §2 |
| 10.0 mmHg | Vaupel 2007's pooled cervix/H&N/breast median pO2 (PMID 17536958); also the broader "hypoxic" cutoff `oxygen_calibration_landscape.md` used for GBM (Lally 2006) | 54.8254% (49,225 cells) |
| 14.0 mmHg | the model's own VIABLE threshold (`t=0.4`) | 95.3734% (85,631 cells) |

**Read plainly: at the literature's own preferred severe-hypoxia cutoff (2.5 mmHg), this model predicts EXACTLY ZERO hypoxic breast tumor cells — the identical implausible zero `oxygen.py`'s docstring already documents for GBM.** The model's solved oxygen field never gets colder than 4.68 mmHg anywhere a tumor cell sits (§1), so no threshold at or below ~4.7 mmHg can ever be crossed by construction, regardless of what fraction is reported at a looser cutoff. The 11.2% headline is real (it is what `update_cell_states` actually computes with the shipping default `hypoxia_threshold=0.2`), but it is a statement about a 7 mmHg cutoff the model invented for itself, not a reproduction of the ~2.5 mmHg cutoff the Vaupel literature (and GBM's own landscape review) treats as the operative definition of "hypoxic." Vaupel's own numbers make substantial sub-2.5-mmHg regions routine in real breast tumors (6 of 15 tumors in the 1991 series, ~25% pooled hypoxic fraction in the 2007 compilation) — so a floor of 4.68 mmHg that is never crossed is, by this project's own established standard for judging such things (`oxygen_calibration_landscape.md` §4's "0.0% is not credible" argument), **equally not credible for breast as it is for GBM.**

**So: same failure mode, different visibility.** GBM's problem is visible immediately because its own threshold reports 0.0%. Breast's identical structural problem — the solved field never reaches genuinely severe hypoxia — is masked by the fact that the model's threshold (7 mmHg) happens to sit above the floor (4.68 mmHg) it never crosses, producing a non-zero, superficially plausible-looking 11.2%. Whether 11.2% (at 7 mmHg) or 0.0% (at 2.5 mmHg) is "the" model's HYPOXIC fraction depends entirely on which of these two numbers is quoted, and neither has been independently validated against a breast-specific measurement — this project has no matched pO2 measurement on the actual Xenium Rep1 section, only population-level literature from different patients. **The honest conclusion is not "breast looks better than GBM" — it is "breast's headline number happens to fall on the more forgiving side of an unresolved threshold-and-calibration problem that afflicts both tumor types identically."** This is the same `OxygenConfig.consumption`-uncalibrated-level problem `oxygen.py`'s own docstring names for GBM, observed here to reproduce on a second, independently reconstructed real tumor.

---

## 5. Reproduction and zero core-engine changes

**Script (not committed under `tmesim/`):** a standalone script that calls `tmesim.reconstruct.real_pipeline.reconstruct_real`, `tmesim.cells.agents.update_cell_states`, and `tmesim.physics.oxygen.to_mmhg` — all pre-existing, unmodified functions — and only reads/reports their outputs. No new keyword arguments, no monkeypatching, no edits to any file under `tmesim/`.

```python
micro = reconstruct_real(source="xenium", sample="data_cache/xenium_breast_janesick/rep1", dx_um=12.0)
oxygen = np.asarray(micro.fields.oxygen)                      # solved by reconstruct_real internally
mmhg = tmesim.physics.oxygen.to_mmhg(oxygen)                   # unmodified helper
tmesim.cells.agents.update_cell_states(micro, oxygen, OxygenConfig())  # unmodified, called explicitly
```

**Confirmed zero edits under `tmesim/` from this task.** `git status`/`git diff` at the time of writing shows `tmesim/core/config.py` and `tmesim/physics/oxygen.py` modified — both are the concurrent Michaelis–Menten-sink / consumption-sensitivity-band task flagged in this task's own briefing, not this report's work. Specifically observed: that concurrent task had already landed `OxygenConfig.sink_form` (default `"linear"`, unchanged behavior) and `mm_half_saturation` / `mm_max_iter` / `mm_tol` fields by the time this script ran, plus the `REFERENCE_PO2_MMHG` / `to_mmhg()` addition this report depends on. This run used every default (`sink_form="linear"`, `consumption=0.02`, `supply=0.05`, `diffusivity=2000`), i.e., the **original linear sink**, not the new Michaelis–Menten option — stated here per the task's instruction to cite what was actually observed rather than assume. No line under `tmesim/` was added, removed, or edited by this task.

**Data provenance.** `data_cache/xenium_breast_janesick/rep1` — Janesick et al. 2023 (PMID [38114474](https://pubmed.ncbi.nlm.nih.gov/38114474/)), 10x Xenium FFPE Human Breast Cancer Rep1, 313-gene targeted panel, the same reconstruction `results/breast_cancer_flagship.md` uses. This section is untreated resected tissue with no paired pO2 or outcome measurement — every number in this report is a model prediction checked against POPULATION-level literature from different patients, not a validated per-section measurement.

**Caveats carried forward, stated plainly:**
1. `REFERENCE_PO2_MMHG = 35 mmHg` is a normal-BRAIN-tissue anchor (Erecinska & Silver 2001) applied here to breast tissue without re-derivation; Vaupel's own normal-breast pO2 (52–65 mmHg depending on cohort) is roughly double it (§1).
2. `OxygenConfig.consumption` remains uncalibrated against any measured pO2 for either tumor type, and no breast-specific OCR literature exists to calibrate it against (§3) — the same gap `oxygen_calibration_landscape.md` found for glioma.
3. The HYPOXIC/QUIESCENT/VIABLE band edges (7 / 14 mmHg) are untouched GBM defaults, not re-derived for breast tissue.
4. `reconstruction_report()['patient_id']` still returns the hardcoded string `"xenium-gbm-real"` on this breast section — a known, previously-documented cosmetic label bug (`results/breast_cancer_flagship.md` §1b), irrelevant to any number reported here but worth flagging so a reader is not confused by it.
5. n = 1 breast section, no outcome label, no paired pO2 measurement — everything here is model-vs-population-literature, not model-vs-ground-truth.
