# Second-pass OCR literature search: did a differently-angled search find what `oxygen_calibration_landscape.md` missed?

**Status: research only. No code changed.** This is a targeted, second-pass follow-up to `results/oxygen_calibration_landscape.md` §2, which searched hard for a glioma-specific tissue-level oxygen consumption rate (OCR) in units usable by `OxygenConfig.consumption` (`tmesim/core/config.py:76-82`, wired into `tmesim/physics/oxygen.py`) and came up short of anything glioma-specific and directly usable — only cell-line Seahorse (2D monolayer) data and Grimes et al. 2014 (PMID 26064525), a non-glioma theoretical paper. This pass approaches the same question from three different angles the first pass did not use: (1) measured-pO2-gradient-plus-fitted-diffusion-model papers for any brain tumor or normal brain tissue, (2) CMRO2 (cerebral metabolic rate of oxygen) as a standard clinical quantity, both normal and glioma-specific (PET-based), and (3) an explicit unit-by-unit conversion attempt into `OxygenConfig.consumption`'s exact units, now read directly from `tmesim/physics/oxygen.py`.

## Units of `OxygenConfig.consumption`, confirmed from source

`tmesim/physics/oxygen.py:181-239` solves the steady state of

```
D * lap(O) + supply * vasc * (1 - O) - consumption * rho_tum * O = 0
```

`O` is dimensionless relative saturation in `[0,1]`. `vasc` (`[0,1]`, clipped) and `rho_tum` (`_tumor_density`, tumor cells per `TUMOR_REFERENCE_AREA_UM2 = 144 um^2` reference tile, a dimensionless occupancy typically O(1)-O(8) per the module's own self-test) are both dimensionless. `D * lap(O)` has units `O / s` (since `D` is `um^2/s` and `lap` divides by `um^2`), so by balance in the PDE, **`consumption` (and `supply`) must have units of `s^-1`** — a first-order rate constant, exactly as the class docstring implies ("per-tumor-cell consumption rate (relative)") but never states in SI-comparable terms. Default value: `consumption = 0.02` (`tmesim/core/config.py:80`). This confirmed unit (`s^-1`, first-order-in-`O` sink) is the target for everything below, and it matters: it rules in one class of literature number (first-order kinetic constants) and rules out, or forces a conversion for, another (zeroth-order volumetric fluxes like Grimes' `a_max` or raw CMRO2).

A parallel, independent piece of work in this same results tree, `results/sweep/oxygen_mmhg_reference_proposal_2026-07-29T05-54-00Z.md`, proposes `REFERENCE_PO2_MMHG = 35.0` mmHg (Erecinska & Silver 2001, PMID 11718758) as the value `O = 1` should be defined to represent — i.e. it resolves the "what does O=1 mean in mmHg" blocker that `oxygen_calibration_landscape.md` §4 flagged as unresolved. That proposal is used below (not re-derived here) wherever a concrete mmHg anchor is needed for a worked conversion, and it is cited, not endorsed independently — this document did not re-verify Erecinska & Silver 2001.

---

## 1. Measured pO2 gradient + fitted diffusion-consumption model — searched again, three real hits, none glioma

This is the exact type of paper the first pass said it could not find for glioma: a study that measures a pO2 gradient at known depth/distance increments and back-calculates a tissue-level consumption rate by fitting a diffusion-reaction model to the gradient's *shape*. A second, differently-angled search (PubMed via `Krogh cylinder`, `reaction-diffusion... consumption rate`, and general web search, rather than the first pass's glioma-first framing) found three real instances of this method. **None is glioma or brain-tumor tissue**, but two are directly comparable tissue types (retina, and a well-perfused solid tumor), and they let the "not found for glioma" conclusion be checked from outside rather than just re-asserted.

### 1a. Dewhirst, Secomb, Ong, Hsu, Gross 1994 — the closest methodological match found

Dewhirst MW, Secomb TW, Ong ET, Hsu R, Gross JF. "Determination of local oxygen consumption rates in tumors." *Cancer Res* 1994;54(13):3333-3336. PMID [8012945](https://pubmed.ncbi.nlm.nih.gov/8012945/).

Oxygen microelectrodes measured PO2 profiles across microregions (87-286 um wide) of R3230AC rat mammary adenocarcinoma in a dorsal skin-flap preparation, with adjacent microvessels mapped. **Measured PO2 values were compared against theoretical (Krogh-type) diffusion simulations and local consumption rates were deduced by the fit** — precisely the method the task asked to search for. Six profiles gave consumption rates **0.83-2.22 cm^3 O2/100 g/min, mean (+/- SD) 1.52 +/- 0.51 cm^3 O2/100 g/min**.

Not glioma, not even brain tissue (mammary adenocarcinoma, rat, extracranial skin-flap window). But it is the real thing methodologically — a directly-fitted, tissue-level, diffusion-model-derived OCR — and its magnitude (see §3 conversion) turns out to land in the same order of magnitude as the PET-based glioma number found in §2, which is a useful independent cross-check even though neither is glioma-*and*-diffusion-fit simultaneously.

### 1b. Zuckerman, Cheasty, Wang 1993 — same method, neural tissue (retina), not glioma

Zuckerman R, Cheasty JE, Wang Y. "Optical mapping of inner retinal tissue PO2." *Curr Eye Res* 1993;12(9):809-825. PMID [8261793](https://pubmed.ncbi.nlm.nih.gov/8261793/), DOI [10.3109/02713689309020386](https://doi.org/10.3109/02713689309020386).

An optical (fluorescence-quenching) oxygen-mapping technique, not an electrode, but methodologically the same idea: **tissue pO2 gradients around arterioles were fit to a Krogh cylinder diffusion model**, yielding an O2 consumption rate of **5.48 +/- 0.50 (SEM) x 10^-3 mL O2/mL tissue/min** for light-adapted rat inner retina, with the Krogh cylinder radius fit at 144.7 +/- 5.5 um. The paper explicitly notes its retinal pO2 histogram (mean 25.0 mmHg) "showed remarkable correspondence to those determined invasively in brain by others, using O2 microcathodes, possibly suggesting a similarity in the underlying capillary architectures of the two neural tissues" — a secondary claim (not independently verified against the cited "others" in this pass), but it is the only hit in this search that explicitly draws a neural-tissue-to-brain correspondence for a fitted OCR value. Still not glioma, not tumor tissue at all.

### 1c. Milotti, Stella, Chignola 2017 — not a new measurement, but the best functional-form match found, with a citation trail to real in vivo tumor numbers already in `s^-1`

Milotti E, Stella S, Chignola R. "Pulsation-limited oxygen diffusion in the tumour microenvironment." *Sci Rep* 2017;7:39762. DOI [10.1038/srep39762](https://doi.org/10.1038/srep39762).

This is a theoretical/simulation paper (general solid tumors, not glioma), but it is unusually valuable for this specific codebase for a structural reason the first-pass report flagged as a problem, not just a units gap: **its reaction-diffusion model uses a first-order-in-concentration sink term, `f(Phi) ~ gamma * Phi`** — i.e. exactly tme-sim's own `-consumption * rho_tum * O` functional form (linear in the oxygen variable, not the zeroth-order/Michaelis-Menten form Grimes 2014 uses and that the first-pass report flagged as the literature default). Confirmed by direct fetch of the paper text (`pmc.ncbi.nlm.nih.gov/articles/PMC5206636/`): "the reaction term f(Phi(r,t), r, t) ~ gamma*Phi," with gamma treated as a first-order rate constant.

Critically, the paper cites a **real, measured, in vivo range for this exact rate constant, already in `s^-1`**: "the rates of oxygen consumption in different areas of in vivo tumours have been elegantly and precisely measured, and they have been shown to vary in the range **1.66 x 10^-4 - 5 x 10^-3 s^-1 (mean value 2.16 x 10^-3 s^-1)**" — citing three papers by Diepart et al. (refs 28-30 in Milotti 2017):

- Diepart C, Jordan BF, Gallez B. "A New EPR Oximetry Protocol to Estimate the Tissue Oxygen Consumption In Vivo." *Radiat Res* 2009;172(2):220-225. PMID [19630526](https://pubmed.ncbi.nlm.nih.gov/19630526/), DOI [10.1667/RR1448.1](https://doi.org/10.1667/RR1448.1).
- Diepart C et al. "Comparison of methods for measuring oxygen consumption in tumor cells in vitro." *Anal Biochem* 2010;396(2):250-256. (PMID not retrieved in this pass; not the same record as the two below.)
- Diepart C, Magat J, Jordan BF, Gallez B. "In vivo mapping of tumor oxygen consumption using 19F MRI relaxometry." *NMR Biomed* 2011;24(5):458-463. PMID [20891023](https://pubmed.ncbi.nlm.nih.gov/20891023/), DOI [10.1002/nbm.1604](https://doi.org/10.1002/nbm.1604).

These Diepart papers derive the `s^-1` rate constant from the *kinetics of pO2 return to baseline after a carbogen-breathing oxygen-saturation challenge* (EPR oximetry and, in the 2011 paper, 19F MRI relaxometry mapping) in mouse liver and muscle/tumor tissue — not from a spatial pO2-gradient fit, so it is a different measurement paradigm than §1a/1b (temporal kinetics vs. spatial gradient fitting), but it estimates the same physical quantity a first-order sink term needs, and does so in exactly `s^-1`. **Not glioma-specific** (mouse tumor models in the Diepart papers are not stated as brain tumors in the abstracts retrieved) and **not independently re-derived here** (taken from Milotti 2017's citation, not fetched from the Diepart primary texts directly) — flagged UNVERIFIED-to-primary-source, same convention `oxygen_calibration_landscape.md` uses elsewhere.

**No fourth hit**: repeated PubMed searches combining `glioma`/`glioblastoma`/`brain tumor` with `Krogh`, `reaction-diffusion`, `microelectrode... diffusion... consumption... fit`, and `xenograft... oxygen... diffusion... model... consumption` returned **zero results** for the glioma-specific combination in every query tried (see query log below). This reproduces, from a different angle, the first pass's core finding: **a paper that fits a diffusion-consumption model to a measured pO2 gradient specifically in glioma or any brain tumor tissue was not found in this pass either.**

---

## 2. CMRO2 — the genuinely new thing this pass found that the first pass missed entirely

`oxygen_calibration_landscape.md` does not mention CMRO2 anywhere. This is the one place this second pass materially outperforms the first: CMRO2 is a standard, decades-old, PET-quantifiable clinical metric, and **glioma-specific CMRO2 does exist via 15O-PET**, including one paper close enough to a glioma tissue model that its number is directly usable as an order-of-magnitude anchor.

### 2a. Normal human brain CMRO2 — well established

Powers WJ, Grubb RL, Darriet D, Raichle ME. "Cerebral blood flow and cerebral metabolic rate of oxygen requirements for cerebral function and viability in humans." *J Cereb Blood Flow Metab* 1985;5(4):600-608. PMID [3877067](https://pubmed.ncbi.nlm.nih.gov/3877067/), DOI [10.1038/jcbfm.1985.89](https://doi.org/10.1038/jcbfm.1985.89). PET study, 50 subjects. Reports the **minimum** rCMRO2 compatible with normal neurological function as **1.3 mL O2/100 g-min**; this is a floor for viability, not the resting-normal point value.

The commonly cited **resting-normal** point value, per general web search of the PET-methods literature (not independently PubMed-verified paper-by-paper in this pass; flagged accordingly): **~130-150 umol O2/100 g/min whole-brain** (one MRI-CMRO2 methods paper reports 132.1 +/- 20.0 umol/100 g/min in young healthy subjects, "in agreement with literature reports using PET"), or equivalently, in the older PET-literature units directly comparable to the glioma number in §2b, **roughly 3.0-3.8 mL O2/100 g/min for gray matter**. Both unit conventions (umol/100g/min and mL/100g/min) appear in the literature; 1 mL O2 (STP) = 44.6 umol, so 3.0-3.8 mL/100g/min = 134-170 umol/100g/min — internally consistent with the ~130-150 umol figure above.

### 2b. Glioma-specific CMRO2 (called TMRO2, "tumoral" metabolic rate of oxygen) via 15O-PET — a real, glioma-specific, tissue-level hit

Watabe T, Kanai Y, Ikeda H, Horitsugi G, Matsunaga K, Kato H, Isohashi K, Abe K, Shimosegawa E, Hatazawa J. "Quantitative evaluation of oxygen metabolism in the intratumoral hypoxia: 18F-fluoromisonidazole and 15O-labelled gases inhalation PET." *EJNMMI Res* 2017;7:16. PMID [28210996](https://pubmed.ncbi.nlm.nih.gov/28210996/), PMC [5313496](https://pmc.ncbi.nlm.nih.gov/articles/PMC5313496/), DOI [10.1186/s13550-017-0263-6](https://doi.org/10.1186/s13550-017-0263-6).

Ten Fischer rats bearing **C6 glioma** (a standard rat glioma cell line, the same tissue type used throughout the glioma-modeling literature) were imaged with both 18F-FMISO PET (hypoxia marker) and the **steady-state 15O-labelled-gas inhalation method, methodologically identical to human CMRO2 PET** (the paper explicitly follows the standard Frackowiak steady-state-inhalation quantification used for human CMRO2). Tumor regions were stratified into three hypoxia bands by FMISO uptake (Low/Medium/High SUV) and TMRO2 (tumoral metabolic rate of oxygen, exact same physical quantity and units as CMRO2) was reported for each:

| FMISO group | TBF (mL/100mL/min) | OEF (%) | **TMRO2 (mL O2/100 mL/min)** | TBV (mL/100mL) |
|---|---:|---:|---:|---:|
| Low (non-hypoxic) | 55 +/- 30 | 33 +/- 14 | **2.8 +/- 1.3** | 5.7 +/- 2.1 |
| Medium (moderately hypoxic) | 32 +/- 17 | 36 +/- 17 | **1.9 +/- 1.0** | 4.3 +/- 1.9 |
| High (severely hypoxic) | 30 +/- 15 | 41 +/- 16 | **2.1 +/- 1.1** | 3.9 +/- 1.2 |

Hypoxic tumor regions showed significantly lower TBF, TMRO2, and TBV than non-hypoxic regions, with OEF rising in the severe-hypoxia group (increased extraction from reduced supply) — internally consistent physiology. **This is a real, PET-measured, glioma-specific (C6 rat glioma), tissue-level oxygen consumption number, in the exact same units (mL O2/100 mL tissue/min) as normal-brain CMRO2**, which makes a like-for-like comparison possible: glioma TMRO2 (1.9-2.8) sits noticeably below the normal-brain range (~3.0-3.8, §2a) — consistent with the "luxury perfusion" / reduced-oxygen-demand phenomenon described elsewhere in the glioma-PET literature (e.g., the comparison group in Mineura K, Sasajima T, Kowada M, Shishido F, Uemura K, "[Positron emission tomography (PET) study in patients with meningiomas]" *No To Shinkei* 1990;42(2):145-151, PMID [2357416](https://pubmed.ncbi.nlm.nih.gov/2357416/) — Japanese-language, abstract only, reports comparison malignant-glioma and contralateral-gray-matter rCMRO2 values but the glioma-specific number itself was not extractable from the retrieved abstract in this pass, flagged UNVERIFIED).

**Caveat, stated plainly:** rat C6 xenograft, not human GBM tissue; PET-derived (a whole-voxel functional average over a ~1mm-diameter sphere of interest per the paper's own methodology), not a spatial-gradient diffusion-model fit like §1a/1b. It answers a different (but arguably more directly useful) question than the §1 papers: not "what's the local point consumption rate implied by a measured gradient's shape," but "what's the actual measured volumetric O2 flux this glioma tissue is consuming, in the same units clinical CMRO2 uses." Both are legitimate routes to the same underlying physical quantity.

---

## 3. Converting CMRO2/TMRO2 into `OxygenConfig.consumption`'s units — worked, with every assumption shown

This is the exercise the task asked for explicitly: "work out (show your steps) whether CMRO2 can be converted into the same units `tmesim/physics/oxygen.py`'s `consumption` constant uses."

**Step 1 — establish a physical O2-solubility constant.** The PDE's sink term `consumption * rho_tum * O` is a *fractional saturation loss rate*; to compare it to a physical O2 flux (mL O2/mL tissue/time) you need to know how much O2 (mL) a unit of tissue holds at `O = 1`. This requires Henry's-law tissue O2 solubility, a real physical constant, not a free modeling choice. Grimes et al. 2014 (`oxygen_calibration_landscape.md` §2) already supplies this implicitly: its conversion constant `Omega = 3.0318 x 10^7 mmHg*kg/m^3` relates consumption `a(p)` (m^3 O2/kg/s) to a pressure-drop rate via `D*lap(p) = Omega*a(p)`, i.e. `1/Omega` is the tissue O2 solubility in m^3 O2 per kg tissue per mmHg. Converting units (tissue density ~1000 kg/m^3 ~ 1 g/mL): `1/Omega = 3.299 x 10^-8 m^3/(kg*mmHg) = 3.3 x 10^-5 mL O2/(mL tissue * mmHg)`. This independently-derived value from Grimes' own constant is used below as `alpha = 3.3 x 10^-5 mL O2/(mL tissue * mmHg)` — a genuine cross-check, not an assumption invented for this document, since it falls directly out of a number already vetted in the first-pass report.

**Step 2 — pick a reference pO2 for `O = 1`.** Using the parallel proposal's `REFERENCE_PO2_MMHG = 35.0` mmHg (`results/sweep/oxygen_mmhg_reference_proposal_2026-07-29T05-54-00Z.md`, cited to Erecinska & Silver 2001, PMID 11718758 — not independently re-verified here), the tissue O2 "capacity" at full saturation is:

```
C_ref = alpha * REFERENCE_PO2_MMHG = 3.3e-5 * 35 = 1.155e-3 mL O2 / mL tissue
```

**Step 3 — relate the sink term to a physical flux.** The PDE's sink `consumption * rho_tum * O` (units `O/s`) becomes a physical O2 consumption flux (mL O2/mL tissue/s) when multiplied by `C_ref`:

```
physical_flux = consumption * rho_tum * O * C_ref
=>  consumption = physical_flux / (rho_tum * O * C_ref)
```

**Step 4 — plug in Watabe 2017's glioma TMRO2.** Convert `2.8 mL O2/100 mL/min` (non-hypoxic C6 glioma) to SI-comparable units: `2.8 / 100 / 60 = 4.667e-4 mL O2/mL tissue/s`.

Two illustrative parameter choices (both plausible per `tmesim/physics/oxygen.py`'s own self-test docstring, which uses `rho_tum ~ 8` for a hypoxic core and describes `rho_tum ~ 1` as "roughly one tumor cell per cell-sized tile"):

| Scenario | rho_tum | O (assumed) | TMRO2 used | Implied `consumption` (s^-1) |
|---|---:|---:|---:|---:|
| Sparse tumor, well-oxygenated | 1 | 0.8 | 2.8 (non-hypoxic) | **0.505** |
| Dense hypoxic core | 8 | 0.15 | 2.0 (hypoxic, ~avg of Medium/High) | **0.240** |

**Step 5 — compare to the current default and to Diepart/Milotti's already-`s^-1` number.** `OxygenConfig.consumption = 0.02` (`tmesim/core/config.py:80`). Both scenarios above land **12-25x above** the current default — i.e., this conversion route, taken at face value, suggests the current default may be *too low* relative to measured glioma O2 demand.

But the Diepart-derived, already-unit-matched first-order rate constant from Milotti 2017 (§1c) — which needs **no** `C_ref`/`REFERENCE_PO2_MMHG` conversion at all, since it is already a fitted `s^-1` first-order rate constant analogous to tme-sim's lumped `consumption * rho_tum` product — gives **1.66e-4 to 5e-3 s^-1 (mean 2.16e-3 s^-1)**. Comparing this to the current default's lumped product (`consumption * rho_tum = 0.02 * rho_tum`, i.e. `0.02` to `0.16` s^-1 for `rho_tum` in `[1, 8]`) puts the **current default 4x-100x above** the Diepart/Milotti range instead.

**This is the honest, load-bearing finding of this section: two independent, literature-grounded conversion routes disagree with each other about which direction the current default is wrong, by roughly an order of magnitude each way.** The Watabe/CMRO2 route (§2b, spatially-averaged PET flux in glioma tissue specifically) says the default is too low; the Diepart/Milotti route (§1c, temporal-kinetics-derived first-order rate constant, general tumor, already unit-matched) says it may be too high. Both are real numbers from real, cited sources; neither should be treated as *the* answer. The disagreement itself is informative — it shows the uncertainty here is not just "we don't have a number," it's "the numbers we now have, even restricted to ones in defensible units, don't converge," which is a stronger and more specific form of the first-pass report's "un-calibrated" conclusion, not a resolution of it.

**What this section does NOT do, per scope:** it does not recommend a new value for `consumption`, and it does not resolve the zeroth-order-vs-first-order sink functional-form question `oxygen_calibration_landscape.md` §2 already flagged as a prerequisite modeling decision (Watabe's PET-measured TMRO2 is itself a zeroth-order-style *bulk flux*, not evidence about whether the *local* kinetics are first-order — using it inside a first-order sink term the way §3 does is itself the same category of modeling assumption the first-pass report warned against making silently). This section shows the conversion is *computable*, not that any single resulting number is *right*.

---

## 4. Verdict

Restating the task's own framing directly:

**Did this pass find a better number than `oxygen_calibration_landscape.md` already has, or does its "not found" conclusion hold up under a second, differently-angled search?**

**Both, in different places — and that is the honest answer, not a hedge:**

1. **The specific "not found" conclusion — a diffusion-model fit to a measured pO2 gradient, done specifically in glioma or any brain tumor tissue — holds up.** A second search, using entirely different search terms (Krogh cylinder, reaction-diffusion, CMRO2, cerebral metabolic rate of oxygen, rather than the first pass's hypoxia/pO2-threshold framing) still did not find one. Two real instances of that exact method exist in the literature (Dewhirst/Secomb/Ong/Hsu/Gross 1994, mammary adenocarcinoma; Zuckerman/Cheasty/Wang 1993, retina) — confirming the method is real and used elsewhere, which makes its absence for glioma specifically more notable, not less. This is now evidence from two independent search passes that the gap is genuine rather than an artifact of under-searching.

2. **A materially better, glioma-specific number was found by approaching from an angle the first pass never used (CMRO2/PET) — Watabe et al. 2017's PET-measured TMRO2 in C6 glioma, 1.9-2.8 mL O2/100 mL/min, stratified by hypoxia grade.** This is glioma-specific, tissue-level, in a standard and independently cross-checkable unit (directly comparable to the equally well-established normal-brain CMRO2 figure, ~3.0-3.8 mL/100g/min), and — unlike the Seahorse cell-line numbers the first pass found — does not require an unknown 2D-to-3D cell-packing conversion, because it is already a bulk tissue-volume measurement. It is a genuine addition to the evidence base this project has.

3. **The unit conversion into `OxygenConfig.consumption` is now demonstrably *computable*, which the first pass could not show because it never confirmed the constant's units or found a p_ref proposal to convert against — but doing the computation surfaces a *new*, sharper problem, not a solution.** Two defensible, differently-sourced conversions of "how big should `consumption` be" disagree with each other by roughly an order of magnitude in opposite directions relative to the current default (§3). That is a stronger, more specific form of "un-calibrated" than the first pass could state, precisely because this pass had real candidate numbers to disagree with each other. The first-pass report's structural finding (zeroth-order vs. first-order sink mismatch, `p_ref` never defined) is not resolved by this document — it is reinforced, now with numbers on both sides of it that don't agree.

**Net honest statement for whoever calibrates this next:** do not treat either §2b or the Diepart/Milotti number as a target to hit. Treat their disagreement as a real uncertainty band, and note that resolving it will need the zeroth-vs-first-order sink question settled first (per `oxygen_calibration_landscape.md` §2, unaddressed here) — converting a bulk PET flux into a *local, first-order-in-O* rate constant the way §3 does is itself already smuggling in the "consumption is proportional to O everywhere" assumption tme-sim currently makes and the literature's zeroth-order default disputes.

---

## Appendix: search queries run in this pass (for reproducibility)

PubMed (via MCP `mcp__claude_ai_PubMed__search_articles`): `oxygen diffusion distance brain tumor reaction-diffusion model fit consumption rate` (0 hits); `cerebral metabolic rate of oxygen CMRO2 normal brain value` (104 hits, top 15 reviewed); `Krogh cylinder brain oxygen consumption model microelectrode` (1 hit -> Zuckerman 1993, retina); `glioma oxygen consumption rate reaction-diffusion model` (0 hits); `15O-PET CMRO2 glioma measurement` (0 hits); `glioma CMRO2 PET oxygen metabolism 15O2` (0 hits); `glioblastoma cerebral oxygen metabolism positron emission tomography` (39 hits, top 15 reviewed); `tumor oxygen tension microelectrode diffusion model consumption rate fit` (0 hits); `brain tissue oxygen consumption rate mathematical model vasculature` (25 hits, top 15 reviewed); `cerebral cortex oxygen tension Krogh model diffusion microelectrode consumption` (0 hits); `Secomb tumor oxygen diffusion model consumption rate vasculature simulation` (1 hit -> not relevant on inspection); `tumor spheroid oxygen consumption rate diffusion model fit pO2 gradient` (0 hits); `glioblastoma xenograft oxygen consumption rate reaction diffusion Krogh` (0 hits); `Determination of local oxygen consumption rates in tumors microelectrode` (1 hit, wrong paper by keyword match); `Groebe Mueller-Klieser oxygen consumption rate tumor spheroid diffusion coefficient` (0 hits); `Mineura glioma positron emission tomography oxygen metabolism` (15 hits reviewed -> PMID 2357416); `malignant glioma regional cerebral oxygen metabolism PET rCMRO2` (11 hits reviewed); `Grote Thews cerebral cortex oxygen pressure fields Krogh model` (0 hits); `glioma xenograft oxygen microelectrode pO2 profile diffusion model consumption fit brain` (0 hits); `Diepart Gallez EPR oximetry tissue oxygen consumption in vivo` (2 hits -> PMID 19630526, 20891023); `Diepart in vivo mapping tumor oxygen consumption 19F MRI relaxometry` (1 hit -> PMID 20891023).

WebSearch: `glioma tissue oxygen consumption rate diffusion model fit pO2 gradient mL O2 per 100g per min`; `"Determination of local oxygen consumption rates in tumors" tumor microelectrode 1.52 cm3 O2 100g min` (-> PMID 8012945); `normal human brain CMRO2 standard value umol/100g/min PET`; `glioma CMRO2 15O-PET reduced oxygen metabolism tumor value ml/100g/min pubmed` (-> Watabe 2017); `"C6 glioma" OR glioblastoma xenograft brain oxygen tension microelectrode Krogh cylinder model fitted consumption rate` (no new hits).

WebFetch: `pmc.ncbi.nlm.nih.gov/articles/PMC5206636/` (Milotti 2017, twice, for citation and reference-list detail); `pmc.ncbi.nlm.nih.gov/articles/PMC5313496/` (Watabe 2017, twice, for citation and value extraction); `researchgate.net/publication/15186190...` (403 Forbidden, could not fetch; PMID/details obtained via PubMed metadata lookup instead).
