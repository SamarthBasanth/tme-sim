# Imaging-derived TME parameters: a landscape survey for tme-sim

**Status: research only. No code changed.** This document surveys what clinical and preclinical imaging modalities quantify about the tumour microenvironment (TME), maps those outputs onto the fields tme-sim's transport PDE needs (vascular source, ECM/tortuosity, oxygen, interstitial fluid pressure (IFP), geometry), and answers the four questions posed: what exists as standard-of-care for a GBM patient, which single modality would help most, whether the resolution gap is bridgeable, and what public datasets exist. Every factual claim below is cited; anything I could not independently confirm is marked **UNVERIFIED**.

tme-sim's own numbers, for reference throughout (from this repo, not an external source): production grid spacing `dx = 12 µm` (docs/architecture.md:20); the calibrated perivascular delivery boundary layer is ~45 µm thick; the vascular source field is currently a KDE of `ENDOTHELIAL` (PECAM1+/FLT1+) cell positions from a single 2D, 5 µm-thick Xenium section (docs/MATHEMATICAL_MODEL.md:1363).

**Real section geometry** — the three actual reconstructed sections, all at `dx = 12 µm`:

| Section | Dataset | Grid | Extent | Bounding area | Tissue area |
|---|---|---|---|---|---|
| **S1** flagship, 319-panel | `data_cache/xenium_gbm` | (424, 675) | **8100 × 5088 µm = 8.1 × 5.09 mm** | 41.2 mm² | **27.9 mm²** (67.70% coverage) |
| **S2** held-out, 480-panel IO | `data_cache/xenium_gbm_io` | (905, 956) | **11472 × 10860 µm = 11.5 × 10.9 mm** | 124.6 mm² | **89.5 mm²** (71.83% coverage) |
| **S3** GEO GSE324268 GBMtumor07, 366-panel | `gse324268_366` | 1196 × 880 | ≈14.4 × 10.6 mm | **151.1 mm²** (recorded) | **UNVERIFIED** — not recorded in the files I read |

Sources: `docs/REPRODUCE.md:362` and `:372`; `results/delivery_metrics.md:23`; `results/heldout_patient.md:66-72` (extents, tissue areas, coverage); `results/cell_state_impact.md:81` (S3 grid); `results/classifier_fix.md:112` (S3 bounding area 151.1 mm²).

> **Correction notice.** An earlier draft of this report sized the domain at ~1.44 mm from the `120×120 at Δx = 15 µm` figure in docs/MATHEMATICAL_MODEL.md:93. **That is a synthetic scale illustration, not the production configuration**, and every consequence I drew from it was wrong. See §6 for the documentation defect that caused it. Two caveats on the corrected numbers themselves: (a) the extents I was handed in review (`[8120.0, 5100.0]`, `[10400.0, 7200.0]`) do not exactly match the repo — S1 is recorded as `(8100.0, 5088.0)` (a trivial 0.25% difference, immaterial), but the S3 figure `[10400, 7200]` = 74.9 mm² is **half** the 151.1 mm² recorded in `results/classifier_fix.md:112`; I have used the repo-recorded values throughout and flag the discrepancy for resolution. (b) I nearly compounded the original error by reading the `0.2870` column of `results/cell_state_impact.md:81` as S3's tissue coverage — it is **mean O₂**, per the table header at line 67. S3's tissue area is genuinely unrecorded in the files I read, so S3 appears below on bounding area only.

---

## 1. Modality-by-modality survey

### 1.1 DCE-MRI (Tofts / extended Tofts model) — the central candidate

**What it measures.** Dynamic contrast-enhanced MRI tracks the T1-shortening effect of a gadolinium-based contrast agent (e.g., Gd-DTPA) as it moves from blood plasma into the extravascular extracellular space (EES) and back. Fitting the tissue concentration-time curve to a two-compartment kinetic model yields:

| Parameter | Physical meaning | Units | Range |
|---|---|---|---|
| **Ktrans** | Volume transfer constant: rate of contrast delivery from plasma to EES per unit tissue volume. **Flow-limited regime:** Ktrans = plasma flow per unit tissue volume. **Permeability-limited regime:** Ktrans = permeability–surface-area product (PS) per unit tissue volume. | min⁻¹ | tumour-dependent (see below) |
| **ve** | Volume fraction of EES (interstitial space) per unit tissue volume — literally the "room" available for the tracer, i.e. interstitial volume fraction | dimensionless, 0–1 | typically 0.1–0.5 in tumour |
| **vp** | Fractional plasma volume per unit tissue volume | dimensionless, 0–1 | usually small (<0.1) except highly vascular lesions |
| **kep** | Efflux rate constant from EES back to plasma, = Ktrans/ve | min⁻¹ | derived |

This is the original Tofts standardized nomenclature, proposed as an international standard: Tofts PS et al., "Estimating kinetic parameters from dynamic contrast-enhanced T1-weighted MRI of a diffusable tracer: standardized quantities and symbols," *J Magn Reson Imaging* 1999;10(3):223–32, PMID [10508281](https://pubmed.ncbi.nlm.nih.gov/10508281/), DOI [10.1002/(SICI)1522-2586(199909)10:3<223::AID-JMRI2>3.0.CO;2-S](https://doi.org/10.1002/(SICI)1522-2586(199909)10:3%3C223::AID-JMRI2%3E3.0.CO;2-S).

**Why this matters for tme-sim specifically.** Ktrans is dimensionally and physically extremely close to what our PDE calls the vascular delivery source. In the current model (docs/MATHEMATICAL_MODEL.md), the source term is `S = s · ρ_vasc · (C_blood − C_tissue)`-type flux with rate constant `s` (1/s, calibrated to 0.05) multiplying a [0,1] vascular density field. Ktrans (min⁻¹, convert to s⁻¹ by /60) is literally a measured per-voxel version of that same rate constant — under permeability-limited conditions Ktrans **is** PS per unit tissue volume, i.e. exactly the coefficient that belongs in front of the concentration-gradient source term in a Kety–Schmidt-style two-compartment PDE. ve is a direct, per-voxel measurement of the interstitial volume fraction the model currently treats as a nuisance/calibrated constant. This is the single most direct correspondence between an imaging output and a term in our governing equations — see §3 for whether it is usable at our resolution.

The distributed-parameter model further separates blood flow F from permeability-surface-area product PS, both with units of cm/min·(cm²/100 g tissue) = mL/min/100 g, identical dimensions to flow. Source: [mriquestions.com Ktrans and permeability](https://mriquestions.com/k-trans--permeability.html); Sourbron & Buckley review, PMC [3752406](https://pmc.ncbi.nlm.nih.gov/articles/PMC3752406/).

**Quantitative values in GBM.** Reported baseline Ktrans in glioblastoma contrast-enhancing tumour: median ~0.109 min⁻¹ (range 0.050–0.205 min⁻¹) in one longitudinal cohort (PMC [11119591](https://pmc.ncbi.nlm.nih.gov/articles/PMC11119591/)); another cohort reports median 0.14 min⁻¹ (range 0.06–0.26 min⁻¹); a threshold of Ktrans ≈ 0.08 min⁻¹ has been used to separate recurrent tumour from radiation necrosis. **UNVERIFIED as a single consensus number** — values vary substantially by acquisition protocol, AIF choice, and model (standard vs. extended Tofts vs. distributed-parameter). Treat 0.05–0.25 min⁻¹ as the working GBM enhancing-tumour range.

**Resolution / acquisition.** Clinical DCE-MRI voxels are **1.25 × 1.25 × 3.0 mm** in a conventional protocol (4.0 s scan interval) and **1.1 mm isotropic** in a high-resolution research protocol (1.6 s scan interval), per a glioma reproducibility study ([*Sci Rep* 2021, PMC 8636480](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8636480/)). Note the conventional protocol is strongly anisotropic — slices are 2.4× thicker than the in-plane edge — because temporal and spatial resolution trade off directly in a dynamic acquisition. It is **3D** (multi-slice or volumetric) but each voxel's kinetic parameters are a spatial average over that mm-scale volume. Exact per-voxel dimensions and their consequences are worked through in §3.1–3.4.

**Clinical status.** DCE-MRI perfusion is **common in neuro-oncology research and increasingly in trials**, but is not yet uniform standard-of-care at every centre the way T1/T2/FLAIR/post-contrast T1 is; DSC-MRI is more entrenched as the routine perfusion sequence in many neuro-oncology protocols (see §1.2). Practical Neurology overview: [Brain Tumor Imaging: Advanced MRI Techniques](https://practicalneurology.com/diseases-diagnoses/imaging-testing/brain-tumor-imaging/31602/).

**File format / software.** DICOM (raw dynamic series); parametric maps (Ktrans, ve, vp, kep) are typically derived post-hoc and may be saved back as DICOM secondary captures or NIfTI depending on the pipeline. Common software: OHIF/Osirix plugins, IB DCE, Tissue4D (Siemens), FireVoxel ([firevoxel.org](https://firevoxel.org/docs/html/userguide/dce_models.html)), MIStar, in-house MATLAB/Python (DCEMRI.jl, PyDCE), and commercial platforms (Olea Sphere).

### 1.2 DSC-MRI: rCBV, rCBF, MTT

**What it measures.** Dynamic susceptibility contrast MRI tracks the T2*/T2 signal drop during first-pass of a gadolinium bolus (susceptibility effect rather than T1-relaxivity effect). Standard derived maps: relative cerebral blood volume (**rCBV**), relative cerebral blood flow (**rCBF**), mean transit time (**MTT**), plus leakage-corrected K2 maps, percentage signal recovery (PSR), and peak height (PH).

**Units.** rCBV/rCBF are usually reported as *relative* (ratio to normal-appearing white matter), dimensionless; absolute CBF can be reported in mL/100 g/min with arterial-input-function calibration. MTT in seconds.

**Glioma grading performance.** rCBV is the most sensitive of the DSC parameters for grading: pooled sensitivity 0.92, specificity 0.81, AUC 0.93 in a meta-analysis (PMC [7444047](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7444047/), Bai et al. — **UNVERIFIED exact PMID/citation for meta-analysis, values as reported by search aggregator, cross-check against primary text advised**). A widely cited rCBV threshold of **1.75** separates high- from low-grade glioma with sensitivity 95.0%, specificity 57.5% (AJNR, Law et al. 2003, PMID [14625221](https://pubmed.ncbi.nlm.nih.gov/14625221/), PMC [8148904](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8148904/)); other cohorts report cutoffs from 1.98 to 2.91, reflecting protocol dependence. MTT does **not** reliably separate high- from low-grade glioma (same source).

**3D / resolution.** DSC-MRI is a 2D-multislice or 3D EPI acquisition; in-plane resolution typically ~1.5–2 mm, slice thickness 3–5 mm; effectively 3D coverage of the tumour but coarse through-plane resolution.

**Clinical status: standard-of-care.** DSC-MRI perfusion (rCBV in particular) is **the most established, most routinely acquired perfusion sequence in neuro-oncology** — used for grading, distinguishing tumour recurrence from pseudoprogression/radiation necrosis, and treatment response.

**File format / software.** DICOM source; rCBV/rCBF/MTT maps computed and stored as DICOM or NIfTI via NordicICE, Olea Sphere, IB Neuro, FuncTool, or open-source `dsc_mri_toolbox`/`ANTsX`.

### 1.3 ASL perfusion

**What it measures.** Magnetically labels arterial blood water as an endogenous tracer (no exogenous contrast) and measures its washout into tissue — quantitative cerebral blood flow (CBF) in absolute units.

**Units.** mL/100 g/min (absolute), or relative tumour blood flow (rTBF) normalized to normal-appearing tissue.

**Glioma application.** ASL-derived rCBFmax/CBFmax discriminate glioma progression from treatment-related change with sensitivity/specificity 0.88/0.83 and 0.93/0.84 respectively (systematic review/meta-analysis, PMC [10576519](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10576519/)). ASL CBF correlates with histopathological microvascular density in MRI-guided biopsies (Neuroradiology, PMID [28236048](https://pubmed.ncbi.nlm.nih.gov/28236048/) — direct histology validation of an imaging perfusion metric, notable for our purposes).

**Advantage over DCE/DSC:** no contrast agent needed (usable in renal-impaired patients, repeatable without cumulative gadolinium), quantitative and less sensitive to BBB disruption artefact than Gd-based methods.

**3D / resolution.** 3D pseudo-continuous ASL (pCASL) is now the recommended clinical implementation; typical voxel size ~3–4 mm isotropic (coarser than DCE/DSC due to inherently low SNR of the labeling approach).

**Clinical status.** **Common in research and increasingly in trials**; not yet universal standard-of-care but growing, particularly as a gadolinium-free adjunct to DSC.

**File format / software.** DICOM; CBF maps computed via vendor sequences (Siemens, GE, Philips all ship product ASL) or open tools (ASL toolbox in FSL/BASIL).

### 1.4 Diffusion MRI: ADC and IVIM

**ADC (apparent diffusion coefficient).** Measures the rate of water self-diffusion in tissue, in units of **mm²/s** (typically reported ×10⁻³ or ×10⁻⁶ mm²/s). ADC is reduced where cell density is high because intra/extracellular membranes restrict diffusion — a genuine, if noisy, cellularity/tortuosity proxy. Reported inverse correlation between ADC and histological cellularity in GBM: Spearman's ρ = −0.48 (p = 0.002) (PMC [4965093](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4965093/) / PMID [27467557](https://pubmed.ncbi.nlm.nih.gov/27467557/)). Mean ADC ≈ 1122 ± 337 ×10⁻⁶ mm²/s corresponding to mean cellularity ≈ 3309 ± 1915 cells/mm² in one 37-patient GBM cohort (same source). The correlation is described as **weak-to-moderate**: ADC is confounded by necrosis and edema, both of which also elevate extracellular volume fraction, so it is a lumped proxy for "extracellular volume fraction" rather than a clean tortuosity measurement.

**Does ADC constrain our ECM/tortuosity field?** Partially and indirectly. ADC is inversely related to extracellular volume fraction, which is mechanistically close to what governs effective diffusivity in a porous-medium/tortuosity formulation (D_eff = D_free / τ², with τ related to extracellular volume fraction via Ogston/Mackie-Meares-type relations). But ADC folds together cellularity, necrosis, and edema into one scalar, so recovering a clean ECM tortuosity field from ADC alone requires additional assumptions or complementary structural data (e.g., necrosis segmentation from T2/FLAIR to disambiguate).

**IVIM (intravoxel incoherent motion).** A bi-exponential extension of DWI, using multiple b-values to separate pure diffusion (D, mm²/s) from a pseudo-diffusion/perfusion term (D*, mm²/s, capillary microcirculation) and perfusion fraction (f, dimensionless 0–1, fraction of signal attributable to microvascular flow). Original description: Le Bihan et al. 1988, *Radiology* PMID (foundational, pre-PubMed-indexed-by-this-tool; widely cited as Radiology 1988;168:497-505). In glioma, perfusion fraction f is significantly larger in high-grade (17.5 ± 6.3%) vs. low-grade glioma (5.8 ± 3.8%) (PMC [4677415](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4677415/)). IVIM is valuable specifically **because** it decomposes the ADC signal into a diffusion (tortuosity/cellularity-like) component and a perfusion-like component, partially resolving the confound above — but it does not use a gadolinium bolus, so its "perfusion" term is a much cruder, unitless proxy than DCE/DSC-derived flow.

**3D / resolution.** DWI/IVIM is acquired as 3D multi-slice EPI; typical clinical voxel 2×2×4–5 mm (highly anisotropic, thick slices); IVIM specifically needs multiple b-value acquisitions, extending scan time and often trading spatial resolution for SNR.

**Clinical status.** Standard DWI/ADC is **standard-of-care** in essentially every brain tumour MRI protocol worldwide. IVIM specifically is **research-only / occasional trial use** — it requires a non-standard multi-b-value protocol most clinical centres do not run by default.

**File format / software.** DICOM; ADC maps are computed on-scanner (standard product feature) and archived as DICOM; IVIM parameter maps (D, D*, f) require offline fitting (in-house Python/MATLAB, or research packages) and are less commonly archived as standard series.

### 1.5 Hypoxia imaging: PET tracers and OE-MRI

**18F-FMISO PET.** The most extensively validated hypoxia PET tracer; 2-nitroimidazole trapped in hypoxic (but viable) cells via bioreductive metabolism. Quantified as **tumour-to-blood ratio (TBR)** or tumour-to-normal-tissue ratio (TNR), unitless. A TBR/T:B threshold of **≥1.2** is the standard cutoff defining "hypoxic volume" (HV) in multiple GBM studies (PMID [19091885](https://pubmed.ncbi.nlm.nih.gov/19091885/), J Nucl Med 2008;50:36-44); at 90 min post-injection, mean TBR = 1.72 ± 0.22 across 9 GBM subjects in a dynamic-imaging study (PMID [31520369](https://pubmed.ncbi.nlm.nih.gov/31520369/), DOI [10.1007/s13246-019-00797-8](https://doi.org/10.1007/s13246-019-00797-8)). Hypoxic volume and mean T/B correlate with tumour aggressiveness metrics derived from biomathematical growth models fit to MRI (PMID [19366800](https://pubmed.ncbi.nlm.nih.gov/19366800/), DOI [10.1158/0008-5472.CAN-08-3884](https://doi.org/10.1158/0008-5472.CAN-08-3884)). FMISO-PET has also been used to generate absolute oxygen tension (pO₂) maps in GBM (*Sci Rep* 2017, [10.1038/s41598-017-08646-y](https://www.nature.com/articles/s41598-017-08646-y)) — this is the closest imaging analogue to our oxygen field's target quantity.

**18F-FAZA PET.** Faster clearance than FMISO, giving lower background/tumour uptake ratio but potentially better contrast at earlier imaging times; less validated than FMISO in GBM specifically (review: [Clin Transl Imaging 2020](https://link.springer.com/article/10.1007/s40336-020-00358-0)).

**64Cu-ATSM PET.** Copper-based bioreductive hypoxia tracer; validated against pimonidazole immunohistochemistry in canine tumours (PMC [3403947](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3403947/)); a 62Cu generator variant does not require an on-site cyclotron, improving practical access relative to F-18 tracers requiring nearby cyclotron production. **GBM-specific validation UNVERIFIED in this search pass** — most ATSM literature I found is in other tumour types; treat ATSM-in-GBM as plausible but not confirmed here.

**Units/spatial resolution.** PET is inherently unitless-ratio or, with full kinetic modeling, an absolute rate constant (mL/min/g-type for the K_i from Patlak analysis). PET spatial resolution is coarse: clinical PET/CT or PET/MR scanners resolve **~4–6 mm** at best (worse than any MRI modality here), and it is fully 3D volumetric.

**OE-MRI (oxygen-enhanced MRI).** Measures ΔR1 (change in T1 relaxation rate) induced by breathing 100% O2 vs. air; well-oxygenated tissue shows a positive ΔR1 (dissolved paramagnetic O2 in plasma/interstitial fluid shortens T1), hypoxic tissue does not respond. Validated preclinically against direct pO2 probes and immunohistochemistry (*Cancer Res* 2016, [10.1158/0008-5472.CAN-15-1789](https://aacrjournals.org/cancerres/article/76/4/787/616114/Oxygen-Enhanced-MRI-Accurately-Identifies), oncotarget preclinical glioma model PMC [5458254](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5458254/)). A prospective GBM trial (the MANGO study) is investigating OE- and BOLD-MRI as hypoxia biomarkers specifically in glioblastoma (PMC [10773871](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10773871/)). OE-MRI requires no radiotracer and no ionizing radiation, and is fully 3D at standard MRI resolution (~1–3 mm), a real practical advantage over PET — but it measures a relaxivity change proxy for oxygenation, not a directly quantitative pO2 map the way FMISO-derived pO2 maps attempt to be.

**Clinical status.** All hypoxia-imaging methods here (FMISO/FAZA/ATSM PET, OE-MRI) are **research-only in GBM** — none is part of routine standard-of-care neuro-oncology imaging protocols; they appear in trials and single/multi-centre research studies, not in the default diagnostic package.

**File format / software.** PET: DICOM, reconstructed on scanner console (GE/Siemens/Philips), kinetic modeling and parametric TBR maps typically computed offline (PMOD, in-house MATLAB/Python, QIBA-compliant pipelines). OE-MRI: DICOM source, ΔR1 maps computed offline (custom pipelines; no single dominant commercial package).

### 1.6 Vessel size imaging (VSI) / vessel architectural imaging (VAI)

**What it measures.** Exploits the differential sensitivity of gradient-echo (ΔR2*) vs. spin-echo (ΔR2) signal decay during first-pass DSC bolus passage to large vs. small vessels respectively. The ratio ΔR2*/ΔR2, combined with rCBV, yields an estimate of mean **microvessel radius/diameter** and microvessel density, without needing histology. VAI extends this by tracking the relative timing shift of the ΔR2*/ΔR2 peaks to classify vessel type (arteriole/capillary/venule-dominated) and estimate oxygen saturation changes.

**Units.** Vessel radius in µm (typically reported diameter 2–50 µm range accessible by the technique); vessel density as a relative index.

**Validation.** MRI-derived VSI correlates with histology-derived vessel size index in high-grade human glioma (r = 0.49, p = 0.01) using image-guided biopsy co-registration (*Sci Rep* 2019, PMID [30742782](https://pubmed.ncbi.nlm.nih.gov/30742782/) — **note: this PMID from my search actually resolves to a Jain lab tumour-vessel-normalization review, not the VSI validation paper; the VSI validation paper itself is PMID [30808879](https://pubmed.ncbi.nlm.nih.gov/30808879/), *Sci Rep* 2019, "Validation of vessel size imaging (VSI) in high-grade human gliomas using magnetic resonance imaging, image-guided biopsies, and quantitative immunohistochemistry" — cross-check recommended**). This is one of the only imaging techniques directly histology-validated against vessel morphometry in glioma, making it notable for calibrating/validating our vascular field even where it can't supply it directly.

**3D / resolution.** Requires the same DSC bolus acquisition (dual-echo GE/SE EPI), so shares DSC's ~1.5–2 mm in-plane / 3–5 mm slice resolution; fully 3D coverage.

**Clinical status.** **Research-only.** Requires a non-standard dual-echo DSC sequence not part of routine clinical protocols; used in a handful of specialized neuro-oncology imaging research centres.

**File format / software.** DICOM source; VSI/VAI maps computed with in-house or specialized (e.g., Aarhus/Aalborg group's msmpr toolbox — **UNVERIFIED name**) research pipelines; no widely deployed commercial product identified in this search.

### 1.7 Non-invasive interstitial fluid pressure (IFP) imaging

This is the imaging capability most directly relevant to tme-sim's `physics/pressure.py::solve_ifp` module, and also the least mature.

**Ground truth invasive method (for context):** wick-in-needle (WIN) manometry, direct puncture. Normal interstitial pressure is roughly −2 to 0 mmHg; tumour IFP is commonly 10× higher, with reported values in human tumours up to **60 mmHg** (review, PMC [4965107](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4965107/)). This IFP elevation and its consequences for drug delivery is the "Jain effect" our model already implements — see Jain's own reviews: Martin, Seano & Jain, *Annu Rev Physiol* 2019, PMID [30742782](https://pubmed.ncbi.nlm.nih.gov/30742782/), DOI [10.1146/annurev-physiol-020518-114700](https://doi.org/10.1146/annurev-physiol-020518-114700); Jain, Martin & Stylianopoulos, *Annu Rev Biomed Eng* 2014, PMID [25014786](https://pubmed.ncbi.nlm.nih.gov/25014786/), DOI [10.1146/annurev-bioeng-071813-105259](https://doi.org/10.1146/annurev-bioeng-071813-105259); Goel et al., *Physiol Rev* 2011, PMID [21742796](https://pubmed.ncbi.nlm.nih.gov/21742796/), DOI [10.1152/physrev.00038.2010](https://doi.org/10.1152/physrev.00038.2010).

**Non-invasive imaging estimates of IFP — glioma-specific.** Elmghirbi et al. built a DCE-MRI-based estimator of tumour IFP (TIFP) and hydraulic conductivity in a rat U251 glioma model, validated directly against wick-in-needle: R = 0.76 (p < 0.0001) for MRI-TIFP vs. WIN-TIFP in controls, R = 0.93 (p = 0.014) in bevacizumab-treated animals; estimated hydraulic conductivity K ≈ (2.3 ± 3.1) × 10⁻⁶ mm/mmHg·s (PMID [29524243](https://pubmed.ncbi.nlm.nih.gov/29524243/), DOI [10.1002/mrm.27163](https://doi.org/10.1002/mrm.27163)). A follow-up study by the same group used DCE-MRI in the same U251 model to track acute IFP, interstitial flux, and vascular parameter changes after bevacizumab (PMID [33817893](https://pubmed.ncbi.nlm.nih.gov/33817893/), DOI [10.1002/nbm.4516](https://doi.org/10.1002/nbm.4516)) — this is essentially the exact quantity our `solve_ifp` module computes theoretically, measured empirically via imaging, in glioma, with histology-comparable validation. **This is preclinical rat-model work, not yet demonstrated in human GBM patients (UNVERIFIED for human translation).**

**MR elastography (MRE).** Estimates tissue stiffness (shear modulus) via shear-wave propagation imaging; used as an input to biomechanical models that predict IFP and solid stress, rather than measuring IFP directly (*Sci Rep* 2020, [10.1038/s41598-020-62367-3](https://www.nature.com/articles/s41598-020-62367-3)). Complementary to DCE-based IFP estimation, not a competing direct measurement.

**Clinical status.** **Research-only, and not yet validated in human GBM** as far as this search could establish — this is the least mature of all the modalities surveyed for direct clinical translation, despite being conceptually the best match to one of our five required fields.

**File format / software.** DICOM source (standard DCE-MRI acquisition); IFP/hydraulic-conductivity maps are the output of bespoke research pipelines (Patlak/Logan-plot-based fluid-mechanical or empirical multivariate models), not vendor or open-source standard software.

### 1.8 Preclinical 3D vasculature: micro-CT, synchrotron, photoacoustic

**Micro-CT.** Non-destructive 3D X-ray tomography; minimum voxel size **1–5 µm** with contrast agent perfusion (Microfil, barium-based casting agents); without contrast agent, effective vessel detection limit is closer to **~50 µm** (PMC [7674347](https://pmc.ncbi.nlm.nih.gov/articles/PMC7674347/)). Fully 3D, ex vivo/preclinical (mouse/rat brain, xenograft) — this resolution (1–5 µm) is actually **finer than our own 12 µm production grid**, i.e. micro-CT-derived vasculature could in principle directly source our vascular field at native or finer resolution, in 3D, for a mouse xenograft model (not a human patient).

**Synchrotron radiation micro-CT / phase-contrast tomography (SRµCT/SRXPCT).** Sub-micron to ~1.4 µm resolution capillary imaging in select ROIs using synchrotron X-ray sources; contrast-agent-free phase-contrast variants have been applied to rodent brain tumour models (arXiv [1409.5453](https://arxiv.org/pdf/1409.5453); PMC [4983401](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4983401/)). Highest-resolution 3D vasculature imaging surveyed here, but requires synchrotron beamtime access — a handful of facilities worldwide, not a routine lab technique.

**Photoacoustic imaging / PACT.** Label-free (endogenous haemoglobin contrast) or targeted-probe 3D imaging of vasculature in live small animals; sub-100 µm resolution to ~10 mm depth with all-optical scanners (PMID [22612139](https://pubmed.ncbi.nlm.nih.gov/22612139/)); 3D PACT systems achieve isotropic ~380 µm resolution over 12.8×12.8×12.8 mm volumes (PMC [11256897](https://pmc.ncbi.nlm.nih.gov/articles/PMC11256897/)). Notable for being **in vivo, longitudinal, and label-free**, unlike micro-CT/synchrotron/light-sheet, which are ex vivo/terminal — enabling tracking of the same tumour's vasculature over time.

**Clinical status.** All three are **preclinical/research-only** (small-animal imaging); photoacoustic imaging has some early clinical translation for superficial human tissue (breast, skin) but not brain (skull is a major acoustic barrier).

**File format / software.** Micro-CT: DICOM or vendor-native (Bruker SkyScan, Zeiss); reconstruction/segmentation via Amira, Dragonfly, or open-source (ImageJ/Fiji + BoneJ/plugins). Synchrotron: raw sinogram data plus reconstructed TIFF/HDF5 stacks, facility-specific pipelines. Photoacoustic: mostly proprietary/in-house reconstruction (k-Wave open-source toolbox is common for simulation/reconstruction research).

### 1.9 Tissue clearing + light-sheet microscopy (iDISCO, CLARITY, CUBIC)

**What it measures.** Chemically renders whole organs/tumours/explants optically transparent (lipid removal + refractive-index matching), then images intact 3D vascular networks (typically CD31/PECAM1-immunolabeled, matching our own Xenium endothelial marker) via light-sheet fluorescence microscopy — no physical sectioning, so no registration artefacts.

**Methods.** iDISCO/iDISCO+ and 3DISCO (organic-solvent-based, better for adult/dense tissue), CLARITY (hydrogel-based aqueous clearing), CUBIC (aqueous, good for whole mouse brain), often combined (iDISCO+CUBIC) for optimal clearing plus labeling (academic.oup.com Biology of Reproduction review [10.1093/biolre/ioaa009], PMC [10551892](https://pmc.ncbi.nlm.nih.gov/articles/PMC10551892/)). Applied specifically to brain tumour microenvironment 3D quantification: *Angiogenesis* 2018, [10.1007/s10456-017-9565-6](https://link.springer.com/article/10.1007/s10456-017-9565-6) ("Optical clearing and fluorescence deep-tissue imaging for 3D quantitative analysis of the brain tumor microenvironment").

**Resolution / 3D.** True cellular-to-subcellular resolution (~1 µm lateral, a few µm axial depending on light-sheet objective), genuinely volumetric, no interpolation between planes required (unlike serial-section reconstruction). This is the only preclinical/explant technique surveyed here that gives **native 3D vascular topology at sub-cellular resolution with no reconstruction artefact** — finer than tme-sim's own 12 µm grid.

**Limitation.** Ex vivo, terminal (tissue must be fixed and cleared), and — critically — clearing protocols typically work on **explants, organoids, or intact small organs (mouse brain), not large human surgical specimens** at this resolution; applicability to human GBM resection specimens at full-tumour scale is limited by clearing depth/penetration and antibody penetration into large dense tissue, though whole human organ clearing has been demonstrated at coarser resolution for e.g. dental pulp (PMC [6659648](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6659648/)) and lymphoma nodes (PMC [10750216](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10750216/)).

**Clinical status.** **Research-only**, PDX/mouse-model and small human specimen (biopsy/explant) use.

**File format / software.** Raw light-sheet output: proprietary (e.g., Zeiss Lightsheet Z.1/7, LaVision Ultramicroscope) multi-terabyte stacks, typically converted to OME-TIFF/HDF5/N5/Zarr for downstream analysis; segmentation/vascular network extraction via Imaris, Vesselucida, or open-source (VesSAP, TubeMap from the iDISCO consortium).

### 1.10 CODA (Johns Hopkins) and serial-histology 3D reconstruction

**What it is.** CODA (not an acronym) is a pipeline for non-linear registration of serial histological (H&E) whole-slide sections plus deep-learning semantic segmentation, producing quantitative 3D reconstructions of multi-centimetre-scale tissue volumes at cellular resolution. Kiemen AL et al., "CODA: quantitative 3D reconstruction of large tissues at cellular resolution," *Nat Methods* 2022;19(11):1490-1499, PMID [36280719](https://pubmed.ncbi.nlm.nih.gov/36280719/), DOI [10.1038/s41592-022-01650-9](https://doi.org/10.1038/s41592-022-01650-9). Demonstrated on pancreas, skin, lung, and liver, segmenting up to nine microanatomical structures **including vasculature** to 96.6% accuracy. Digitization at 0.5 µm/pixel lateral resolution; sections analysed at a reduced 2 µm working resolution; sections taken every third slide (source-dependent axial spacing, effectively tens-of-µm axial sampling **UNVERIFIED precise number, reported inconsistently across secondary sources as "12 µm" in one search snippet, which is not internally consistent with "every third section" of typical ~5 µm sections — treat the axial (z) sampling interval as on the order of 10–20 µm and verify against the primary paper before using**).

**Applications beyond original pancreas paper.** Skin/scalp 3D reconstruction of dermis vasculature and immune cell density (*Commun Biol* 2023, [10.1038/s42003-023-04991-z](https://www.nature.com/articles/s42003-023-04991-z)); AI-driven 3D mapping of microvascular invasion in hepatocellular carcinoma (JHU thesis repository, [jscholarship.library.jhu.edu](https://jscholarship.library.jhu.edu/items/27878869-d383-485a-8f30-7dd0687fc3ba) — **UNVERIFIED, thesis-level source not peer-reviewed**); a CODA-inspired pipeline for whole-heart 3D reconstruction of myocarditis with integrated spatial transcriptomics (PMID [41401444](https://pubmed.ncbi.nlm.nih.gov/41401444/)); a more recent CODA extension mapping 3D immune remodeling in pancreatic precancer (PMID [42325993](https://pubmed.ncbi.nlm.nih.gov/42325993/)). **No GBM-specific application of CODA was found in this search** — it is a general-purpose method, plausible for GBM resection specimens, but not yet demonstrated there as far as I can confirm.

**Why this is directly relevant to tme-sim.** CODA is doing, on serial H&E, essentially the same kind of thing tme-sim's `reconstruct/real_xenium.py` does on a single Xenium section — except CODA explicitly stacks and registers many sections to recover a true 3D tissue volume, whereas our current pipeline uses one 5 µm 2D section and treats it as representative. CODA-style 3D serial reconstruction of a GBM specimen, if performed with an endothelial IHC stain (CD31/CD34) on interleaved sections, would give a genuinely 3D vascular density field at cellular resolution from the *same tissue block* — a direct, feasible upgrade path that requires no new imaging technology, only more sections and a registration pipeline.

**Resolution / 3D.** Sub-cellular lateral (0.5–2 µm), tens-of-µm axial section spacing — genuinely 3D but coarser through-plane than in-plane, and (unlike light-sheet) requires physical sectioning with attendant registration/distortion correction.

**Clinical status.** **Research-only**, requires serial-sectioning workflow well beyond standard pathology practice (routine GBM histology is typically a handful of representative sections, not hundreds of serial slides).

**File format / software.** Whole-slide images (typically Aperio SVS, or similar pyramidal TIFF formats) as input; CODA is an open-source MATLAB/Python pipeline (code released with the *Nat Methods* paper); output 3D label volumes typically saved as multi-page TIFF/NIfTI/custom formats for downstream morphometric analysis.

---

## 2. Q1 — What already exists in the imaging archive for a GBM patient at a major centre?

| Tier | Modality / parameter | Basis |
|---|---|---|
| **Standard-of-care** (essentially every GBM MRI protocol) | T1, T1-post-Gd, T2, FLAIR, DWI/**ADC** | Universal brain tumour protocol component |
| **Standard-of-care** (most major neuro-oncology centres) | **DSC-MRI: rCBV, rCBF, MTT** | Entrenched routine perfusion sequence for grading/pseudoprogression workup; strong grading evidence base (§1.2) |
| **Common in trials / many centres, not universal** | **DCE-MRI: Ktrans, ve, vp, kep** | Increasingly used, especially for BBB-permeability and treatment-response questions, but not yet as uniformly protocolized as DSC |
| **Common in trials / growing** | **ASL perfusion** | Gadolinium-free CBF, growing adoption particularly where repeated contrast dosing is a concern |
| **Research-only** | **IVIM** (needs non-standard multi-b-value DWI) | Requires protocol beyond standard DWI |
| **Research-only** | **Hypoxia PET** (FMISO/FAZA/ATSM), **OE-MRI** | No routine hypoxia-imaging standard of care exists in GBM anywhere as far as this search found |
| **Research-only** | **VSI/VAI** (needs dual-echo GE/SE DSC) | Specialized sequence, handful of academic centres |
| **Research-only, and not yet human-translated** | **Non-invasive IFP estimation** | Demonstrated in rat glioma models only (§1.7) |
| **Not applicable to living patients** | Micro-CT, synchrotron, photoacoustic, tissue clearing/light-sheet, CODA | All ex vivo/preclinical or resection-specimen techniques |

**Bottom line for Q1:** for a real GBM patient scanned at a major centre, the imaging archive will essentially always contain **T1/T1-Gd/T2/FLAIR/DWI-ADC and, at most such centres, DSC-MRI perfusion maps (rCBV/rCBF/MTT)**. DCE-MRI and ASL are common but not guaranteed. Everything with a direct quantitative link to vascular permeability-surface-area (beyond the coarse rCBV proxy), oxygenation, or IFP is either trial-only or has never left the animal-model stage.

## 3. Q3 — The resolution mismatch (addressed first, because it gates everything else)

### 3.1 The imaging voxel dimensions I assume

Anisotropy matters enormously here, so in-plane and through-plane are stated separately, with citations:

| Modality | In-plane | Slice thickness | In-plane footprint | Source |
|---|---|---|---|---|
| **DSC-MRI**, conventional 3T | 1.875 × 1.875 mm (128² matrix) | **5 mm**, 15–20 slices | 3.52 mm² | [PMC 5766430](https://pmc.ncbi.nlm.nih.gov/articles/PMC5766430/); ASFNR recommendations [AJNR 36(6):e41](https://www.ajnr.org/content/36/6/e41), [PMC 5074767](https://pmc.ncbi.nlm.nih.gov/articles/PMC5074767/); HGG consensus protocol [PMC 7523451](https://pmc.ncbi.nlm.nih.gov/articles/PMC7523451/) |
| **DCE-MRI**, conventional | 1.25 × 1.25 mm | **3.0 mm** | 1.56 mm² | [*Sci Rep* 2021, PMC 8636480](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8636480/) |
| **DCE-MRI**, high-resolution research | 1.1 × 1.1 mm | **1.1 mm** (isotropic) | 1.21 mm² | same |
| **Hypoxia PET** (FMISO/FAZA/ATSM) | ~4–6 mm | ~4–6 mm | ~25 mm² (at 5 mm) | §1.5 |

Note the anisotropy is severe for the two routine sequences: DSC slices are **2.7× thicker** than the in-plane voxel edge, conventional DCE **2.4×**. Only the high-resolution research DCE protocol is isotropic.

### 3.2 The per-edge ratio — the number that governs the source-field question

- DSC in-plane 1.875 mm / 12 µm = **156×**
- DCE conventional 1.25 mm / 12 µm = **104×**
- DCE high-res 1.1 mm / 12 µm = **92×**

So one clinical imaging voxel edge is **92–156× our grid spacing**. Because our simulation is a **2D** grid (not a 3D volume — an important correction to my earlier draft, which wrongly cubed this ratio), the right areal statement is how many of our 2D cells sit under one imaging voxel *footprint*:

- DSC: 3.52 mm² / (12 µm)² = **~24,400 sim cells per voxel footprint**
- DCE conventional: **~10,900**
- DCE high-res: **~8,400**

**This is the number that kills the source-field idea, and it is unchanged by the domain-size correction.** ~10⁴ of our finite-volume cells lie under a single imaging voxel. A Ktrans/rCBV/ADC value carries exactly zero information about where, among those ~10⁴ cells, the individual capillaries sit — which is precisely what our 12 µm grid and 45 µm calibrated perivascular boundary layer exist to resolve (docs/architecture.md:20-21). A Ktrans map cannot distinguish "a few closely-spaced high-permeability microvessels" from "many evenly-spread lower-permeability microvessels" producing the same averaged transfer constant; both give one voxel value, but they produce very different micro-scale penetration profiles — the exact quantity tme-sim computes.

### 3.3 The z-mismatch — a separate point, and it is worse than the in-plane one

A Xenium section is **5 µm thick**. The imaging slice it would be compared against is 1.1–5 mm thick:

- vs. DSC (5 mm slice): **1000×**
- vs. conventional DCE (3 mm slice): **600×**
- vs. high-res DCE (1.1 mm slice): **220×**

Two consequences, and the second is the one that matters:

1. The measured Ktrans/rCBV at any in-plane location is an average over a tissue column **220–1000× thicker** than the section actually samples. Tumour composition varies substantially over 3–5 mm in z (rim → core → necrosis), so this is a large, irreducible noise term on any section-to-voxel comparison — not a registration error that better methods could remove, but a genuine physical mismatch of sampled volumes.
2. **The entire section lies inside a single imaging slice in z, regardless of slice thickness.** So the pairing is intrinsically 2D — a section plane against the in-plane voxel grid of *one* slice. There is no third dimension to the comparison. This is why the voxel counts below are computed as in-plane areas, not volumes.

### 3.4 How many imaging voxels actually cover our real sections

This is where the corrected domain size changes the picture. Voxel counts = section area / in-plane footprint:

| Section | Area basis | DSC (3.52 mm²) | DCE conv (1.56 mm²) | DCE hi-res (1.21 mm²) | PET (~25 mm²) |
|---|---|---|---|---|---|
| **S1** flagship | bounding 41.2 mm² | 12 | 26 | 34 | 1.6 |
| **S1** flagship | **tissue 27.9 mm²** | **8** | **18** | **23** | 1.1 |
| **S2** IO | bounding 124.6 mm² | 35 | 80 | 103 | 5.0 |
| **S2** IO | **tissue 89.5 mm²** | **25** | **57** | **74** | 3.6 |
| **S3** GSE324268 | bounding 151.1 mm² | 43 | 97 | 125 | 6.0 |

**The corrected finding:** our real sections span roughly **8–34 perfusion voxels** (S1, depending on protocol and whether tissue or bounding area is used) up to **25–125** (S2, S3) — not "less than one," which is what I wrongly claimed from the synthetic 120×120 figure. That is a materially different regime, and §3.5 works out what it does and does not buy. For hypoxia PET the correction is smaller but still real: S1 covers ~1–2 PET voxels rather than "a fraction of one," and even S3 only ~6 — PET remains far too coarse for any within-section regional analysis.

Preclinical micro-CT (1–5 µm) and light-sheet/CODA (0.5–2 µm) remain the exceptions: **finer than or comparable to** our 12 µm grid, and genuinely capable of supplying source-field-resolution vasculature — but only for mouse/xenograft tissue or resected human specimens, never a living patient.

### 3.5 Verdict on Q3

The source-field conclusion **stands unchanged and should not be softened**: at 92–156× per edge and ~10⁴ sim cells per voxel footprint, no clinical imaging modality can supply a spatially resolved vascular, ECM, oxygen, or IFP source field at 12 µm. Painting one averaged value across ~10⁴ grid cells discards the heterogeneity the model exists to resolve; it would be strictly worse than our current measured single-cell-resolution endothelial KDE, imperfect and 2D as that is. This is not fixable by interpolation — the sub-voxel information is not in the data.

What the corrected domain size changes is the **consequence**, not the verdict. With 10–125 voxels covering a section rather than a fraction of one, imaging is no longer limited to supplying a single scalar. It could in principle supply a *regional map* across the specimen — which turns "calibration against one number" into a genuinely falsifiable multi-point correlation test. Whether that test is actually runnable is §3.6, and the answer is more constrained than the voxel count alone suggests.

### 3.6 Q2b — Is a WITHIN-SECTION regional validation actually possible?

**The proposed experiment.** Partition the section into the ~10–125 tiles corresponding to imaging voxels; compute simulated regional mean exposure (or reconstructed vascular density) per tile; correlate against measured Ktrans (or rCBV) in the matching voxel. A positive correlation would be real, external, falsifiable evidence that the reconstruction's vascular heterogeneity tracks independently measured perfusion. That is a much stronger test than a single scalar comparison.

**Power.** Using Fisher-z (α = 0.05 two-sided, 80% power; N = 7.85/(atanh r)² + 3):

| N (paired voxels) | Min. detectable r @80% power | Power to detect **r = 0.5** |
|---|---|---|
| 8 (S1 tissue, DSC) | 0.85 | **23%** |
| 18 (S1 tissue, DCE conv) | 0.62 | **57%** |
| 23 (S1 tissue, DCE hi-res) | 0.56 | **69%** |
| 25 (S2 tissue, DSC) | 0.54 | **73%** |
| 57 (S2 tissue, DCE conv) | 0.36 | **98%** |
| 74 (S2 tissue, DCE hi-res) | 0.32 | **>99%** |
| 97 (S3 bbox, DCE conv) | 0.29 | **>99%** |

**What effect size should we expect?** The published imaging-vs-histology correlations in glioma are the honest benchmark, and they cluster around r ≈ 0.5: VSI-MRI vs. histological vessel size **r = 0.49** (PMID [30808879](https://pubmed.ncbi.nlm.nih.gov/30808879/)), ADC vs. histological cellularity **ρ = −0.48** ([PMC 4965093](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4965093/)). These are studies by groups who do this professionally, with purpose-built registration. r ≈ 0.5 is the realistic ceiling, not a pessimistic floor.

Against that benchmark: **the flagship section S1 alone is underpowered** (23–69% power depending on protocol — it would miss a true r = 0.5 between a third and three-quarters of the time). S2 and S3 at DCE resolution are nominally well powered (98%+). So the answer is protocol- and section-dependent, and the flagship is the *worst* case, not the best.

**The spatial-autocorrelation caveat, which eats much of that power.** Nominal N overstates effective N here, badly. Both fields are spatially smooth: tumour regional structure (necrotic core / enhancing rim / infiltrative margin) has a correlation length of several mm, and the imaging point-spread function plus partial-volume effects smooth the measured side further. Adjacent 1.25 mm voxels are therefore strongly dependent, and effective degrees of freedom could plausibly be 2–4× below nominal — which would knock S2's N = 57 down to ~14–28 effective, back into the marginal band. **UNVERIFIED as a quantitative factor** (it depends on the actual variogram, which nobody has measured for this pairing), but directionally certain. The honest mitigation is to treat *specimens/patients* as the independent replicate unit and pool across them, rather than leaning on within-section voxel count.

**The registration problem — this is the real blocker, not the power.** The experiment requires knowing where each section tile sits in the pre-operative imaging volume. That demands: (a) the tissue block's 3D position and orientation in the imaging frame, and (b) non-rigid correction for brain shift at craniotomy, formalin/paraffin shrinkage, and sectioning/mounting distortion (folding, tearing, arbitrary rotation and flipping).

Two established approaches, with very different applicability:

- **Point-biopsy via stereotactic neuronavigation — mature, and this is what the field actually does.** The biopsy coordinate is recorded at the moment of sampling and mapped back to pre-op MRI. This is exactly the method behind the VSI validation (PMID [30808879](https://pubmed.ncbi.nlm.nih.gov/30808879/)) and the ASL-vs-microvascular-density study (PMID [28236048](https://pubmed.ncbi.nlm.nih.gov/28236048/)). It works, it is published, and its realized correlations are the r ≈ 0.49 quoted above. **But it yields one paired point per biopsy — typically a handful per patient — not a 10–125-tile map from a single section.**
- **Whole-mount specimen registration — mature in prostate, hard-to-infeasible in GBM.** Methods like RAPSODI and statistical-deformation-model registration ([arXiv 1907.00324](https://arxiv.org/pdf/1907.00324); *Med Phys* 2024 [10.1002/mp.16753](https://aapm.onlinelibrary.wiley.com/doi/10.1002/mp.16753); *Sci Rep* 2021 [10.1038/s41598-021-85028-5](https://www.nature.com/articles/s41598-021-85028-5)) reconstruct serial whole-mount histology in 3D and register it to MRI, often aided by patient-specific 3D-printed sectioning molds. These are prostate-dominated **because the prostate is removed intact as a whole organ with preserved geometry.** GBM resection is typically piecemeal debulking (often ultrasonic aspiration), so there is usually **no intact specimen geometry to register at all.** This is the structural blocker, and it is not a methods-maturity problem that will be solved by a better algorithm.

**And the first-order blocker, which precedes both:** none of our three sections has any paired MRI. S1 and S2 are open 10x demo sections and S3 is GEO GSE324268 — public specimen data with no accompanying patient imaging (`results/cell_state_impact.md:278` notes none carries a paired drug-response label either). **The within-section experiment cannot be run on the data we currently hold, at any voxel count.** It requires prospective collection — a patient with pre-op DCE-MRI, an en-bloc or neuronavigation-localized resection sample, and spatial transcriptomics on that sample — or finding an existing cohort where all three coexist. I did not find such a public cohort (§5).

**Honest bottom line on Q2b:** the power arithmetic is *not* the binding constraint — pooling across specimens solves that. The binding constraints are (i) no paired imaging exists for our current sections, and (ii) piecemeal GBM resection makes whole-section-to-voxel-grid registration infeasible in the general case, forcing a fallback to the point-biopsy design (a few localized samples per patient, pooled across many patients). That fallback design is real, established, and worth pursuing — but it is a prospective multi-patient study, not an analysis we can run on what we have.

### 3.7 Where clinical imaging genuinely helps, restated

1. **Regional multiplier / boundary condition** at mm scale — 10–125 tiles per section is enough to set a spatially varying envelope, not just a global constant. This is a real upgrade over what I claimed in the first draft.
2. **A multi-point falsifiable calibration test** — conditional on solving the provenance and registration problems in §3.6.
3. **Biopsy targeting** — DCE/DSC/hypoxia-PET maps identify which region (high-Ktrans rim, hypoxic core) to sample for the Xenium/CODA/light-sheet data that actually supplies source-field resolution.
4. **Preclinical micro-CT/light-sheet/CODA** remain the only route to genuine source-field-resolution 3D vasculature — ex vivo or animal models only.

## 4. Q2 — Which single modality would most improve the model, and by how much?

Given the resolution finding above, "improve the model" has to mean something achievable — not replacing the source field, but improving what imaging *can* contribute (§3, points 1–3).

**Argument for DCE-MRI (Ktrans/ve/vp) as the single most valuable addition, used as a regional constraint rather than a source field:**

1. It is the only modality on this list whose derived parameter has an essentially 1:1 dimensional and physical correspondence to a term already in our PDE (Ktrans ↔ our vascular delivery rate constant × density field; ve ↔ our interstitial volume fraction). No other modality maps this cleanly onto our equations — DSC-MRI's rCBV is a blood *volume* fraction proxy (closer to our ρ_vasc density field conceptually, but without the permeability information Ktrans carries), and hypoxia PET/OE-MRI map onto the oxygen field only indirectly via a resistance/consumption relationship we would have to additionally model.
2. Unlike VSI/VAI, IVIM, or IFP-imaging, DCE-MRI is **common in trials at exactly the population tme-sim would need for validation** (GBM patients with paired outcome data) — see the TCIA collections in §5 (QIN-GBM-TR, RIDER-NEURO-MRI both have DCE-MRI). This means DCE-derived Ktrans/ve can be obtained today, from public data, for real GBM patients, whereas IFP-imaging is not yet human-validated at all (§1.7) and micro-CT/light-sheet require new specimen acquisition.
3. Concretely, what it would buy: a **regional rescaling of the vascular source rate constant `s`** (docs/MATHEMATICAL_MODEL.md's `s = 0.05 /s` calibrated constant) by measured Ktrans — turning a spatially uniform calibrated constant into a spatially varying field anchored to real per-patient measurement at the mm scale the imaging actually resolves. **Corrected from my first draft:** because a real section spans 18–125 DCE voxels (§3.4), this is not a 2–4 zone piecewise approximation — it is a genuine mm-scale map across the specimen, which is a substantially better upgrade than I originally credited.
4. It would supply an external validation *test*, not merely a number. My first draft framed this as "does our patch-averaged flux land in the observed Ktrans range (~0.05–0.25 min⁻¹)" — a single scalar sanity check that is nearly unfalsifiable given how wide that range is. With 18–125 paired tiles per section, the test becomes: **does simulated regional exposure correlate with measured Ktrans across voxels within the same specimen?** That is a real, falsifiable, directional hypothesis with a defined null. This is the single biggest change the corrected domain size makes to the report's conclusions.

**Does the correction change the Q2 answer? — checking the review's proposal rather than agreeing with it.** The proposal was that 10–40 voxels strengthens the regional-calibration case from a single scalar to a real falsifiable test. **I agree with that specific claim, and §3.6's power arithmetic supports it** — at N = 57–125 (S2, S3 at DCE resolution) the test has >98% power against the realistic r ≈ 0.5 effect size. So the DCE-MRI recommendation is *strengthened*, and I have upgraded points 3 and 4 accordingly.

**But I do not agree that voxel count was ever the binding constraint, and I want to be explicit about that rather than nod along.** Three things temper it, in increasing order of severity:

- **The flagship section specifically is underpowered** (23–69% power, §3.6). If the intended demonstration is on S1 — the section every other flagship result in this repo uses — the answer is that 8–23 tissue voxels is *not* enough against r ≈ 0.5. The well-powered configurations are S2 and S3, which are not the flagship. Anyone planning this experiment should plan it on S2/S3 or on pooled specimens, not on S1.
- **Spatial autocorrelation plausibly cuts effective N by 2–4×**, which returns even S2 to the marginal band (§3.6).
- **Decisively: none of our three sections has paired MRI, and piecemeal GBM resection blocks whole-section registration in the general case** (§3.6). So the experiment is not runnable today at *any* voxel count, and the realistic fallback — stereotactic point biopsies pooled across patients — collapses back to a few points per patient anyway, making within-section voxel count largely moot.

**By how much would this help — honestly, and revised upward but bounded:** it is still a calibration/validation improvement rather than a mechanistic one; it does not change what the model resolves at the micro-scale. But it is a *better* calibration improvement than I first said: a multi-point falsifiable regional test rather than a single-scalar range check, and a mm-scale source-rate map rather than a 2–4 zone approximation. The cost of that upgrade is not analytical — it is prospective specimen collection with imaging provenance, which is the expensive part. And the modality that would help the *micro-scale source field itself* remains non-clinical: ex vivo light-sheet or CODA-style 3D vascular reconstruction of resected tissue (§1.9, §1.10).

## 5. Q4 — Public datasets: TCIA collections pairing imaging with tissue/outcome

| Collection | Content | Patients/studies | Format | URL |
|---|---|---|---|---|
| **IvyGAP** (Ivy Glioblastoma Atlas Project) | MRI/CT (pre/post-surgery/follow-up) linked via de-identified specimen numbers to Ivy GAP's cellular-resolution ISH histology and RNA-seq by anatomic structure (leading edge, infiltrating tumor, cellular tumor, microvascular proliferation, pseudopalisading-around-necrosis) — the same anatomic zones tme-sim already uses for Ivy GAP correlation | 39 patients, 390 studies, 139.98 GB | DICOM | [cancerimagingarchive.net/collection/ivygap](https://www.cancerimagingarchive.net/collection/ivygap/); manifest: `GC_manifest_IvyGAP_20260326.csv`. **Restricted license required** (facial-reconstruction risk in some images) — email help@cancerimagingarchive.net |
| **TCGA-GBM** | Radiology imaging (structural MRI) matched to TCGA genomic/clinical data (clinical/genomic data lives in GDC, not TCIA) | Large TCGA GBM cohort | DICOM | [cancerimagingarchive.net/collection/tcga-gbm](https://www.cancerimagingarchive.net/collection/tcga-gbm/) |
| **BraTS-TCGA-GBM** | Segmentation labels + radiomic features for pre-op TCGA-GBM scans | Subset of TCGA-GBM | DICOM + derived segmentations | [wiki.cancerimagingarchive.net BraTS-TCGA-GBM](https://wiki.cancerimagingarchive.net/pages/viewpage.action?pageId=24282666) |
| **QIN-BRAIN-DSC-MRI** ("Glioma DSC-MRI Perfusion Data with Standard Imaging and ROIs") | Post-contrast T1w + co-registered DSC-MRI, with binary ROIs (lesion, AIF, normal white matter, normal cortex, whole brain) | 49 subjects, 52 studies, 30.76 GB | DICOM | [cancerimagingarchive.net/collection/qin-brain-dsc-mri](https://www.cancerimagingarchive.net/collection/qin-brain-dsc-mri/); DOI [10.7937/K9/TCIA.2016.5DI84Js8](https://doi.org/10.7937/K9/TCIA.2016.5DI84Js8) |
| **RIDER-NEURO-MRI** | Test-retest (repeat scans ~2 days apart, same 1.5T scanner) DCE-MRI with T1-mapping, DTI, contrast-enhanced 3D FLASH/FLAIR — 19 recurrent GBM patients | 19 patients, 108 studies, 368 series, 70,220 images, 7.79 GB | DICOM | [cancerimagingarchive.net/collection/rider-neuro-mri](https://www.cancerimagingarchive.net/collection/rider-neuro-mri/); DOI [10.7937/K9/TCIA.2015.VOSN3HN1](https://doi.org/10.7937/K9/TCIA.2015.VOSN3HN1) |
| **QIN GBM Treatment Response** (QIN-GBM-TR) | Repeated DCE-MRI (2 baseline scans 2–5 days apart, post-surgery pre-therapy) for perfusion repeatability studies | Reported in *Neuro-Oncol Adv* study using this cohort | DICOM | Referenced via [academic.oup.com/noa/article/3/1/vdab174](https://academic.oup.com/noa/article/3/1/vdab174/6433056) (PMID/DOI **UNVERIFIED** in this pass — locate directly on TCIA portal under "QIN GBM Treatment Response") |

**Also flagged but not deeply verified this pass:** RIDER Neuro MRI overlaps with "test-retest reproducibility" use cases rather than outcome-paired use cases; none of the above directly pairs **hypoxia-PET or IFP-imaging** with GBM outcome at scale as a public downloadable dataset — I found no equivalent public FMISO-PET-GBM or IFP-GBM TCIA collection in this search pass. If one exists it is **UNVERIFIED / not found**.

A recent review specifically surveying public glioma MRI datasets is a good starting map for further discovery: "Exploring Adult Glioma through MRI: A Review of Publicly Available Datasets to Guide Efficient Image Analysis," PMC [11773385](https://pmc.ncbi.nlm.nih.gov/articles/PMC11773385/), arXiv preprint [2409.00109](https://arxiv.org/html/2409.00109v1).

---

## 6. Documentation defect that caused the error in the first draft

This should be fixed in the repo, because it will mislead the next reader exactly as it misled me.

**The defect:** `docs/MATHEMATICAL_MODEL.md:93` states "The tissue domain is ~1.8 mm across (120×120 at Δx = 15 µm)". In context this sits inside a scale argument for why all-atom MD is infeasible, so it is a synthetic illustration — but it is presented in the declarative voice ("The tissue domain **is**"), and it is **the only explicit domain-extent-in-millimetres figure anywhere in the mathematical model document.** A reader looking to answer "how physically large is the simulated domain?" finds 1.8 mm and finds nothing contradicting it. `docs/MATHEMATICAL_MODEL.md:579` reinforces the impression by quoting `_test_runtime_120` asserting "a 120×120 grid," and `:575` refers to "hour-long windows on a 120×120 grid" as if that were the production case.

The true production extents (8.1 × 5.09 mm, 11.5 × 10.9 mm, ≈14.4 × 10.6 mm) live only in `docs/REPRODUCE.md:362,372`, `results/delivery_metrics.md:23`, and `results/heldout_patient.md:72` — i.e. in the reproduce guide and results files, never in the architecture or mathematical-model documents that a reader consults for model structure. The 120×120 figure is off by **4.5×–8× per edge and 20×–75× in area** from every real section.

**Suggested fixes:**
1. At `docs/MATHEMATICAL_MODEL.md:93`, mark the 120×120 figure explicitly as a synthetic scale illustration and add the real extents alongside it.
2. Add the real section geometry table (grid, extent, tissue area for all three sections) to `docs/architecture.md` next to the existing "Production resolution is `dx = 12 µm`" statement at line 20 — that is where a reader looks for production configuration, and it currently states resolution without ever stating extent.
3. Minor, same family: `docs/architecture.md:339` documents the signature as `Domain(nx, ny, dx=10.0, ...)`, whose default `dx = 10.0` contradicts the "Production resolution is `dx = 12 µm`" statement 319 lines above. Worth a parenthetical noting the default is not the production value.
4. Consider resolving the S3 extent discrepancy noted in the correction notice at the top of this document (`[10400, 7200]` = 74.9 mm² as given in review vs. 151.1 mm² recorded in `results/classifier_fix.md:112`).

---

## Closing verdict

**Q3 is the crux, and the honest answer is that it kills the "imaging as a source field" version of this idea — for a reason that survives the correction to my domain-size arithmetic.** Clinical DCE-MRI/DSC-MRI in-plane voxels are **92–156× coarser per edge** than our 12 µm grid, putting **~8,400–24,400 of our finite-volume cells under a single imaging voxel footprint**; the through-plane mismatch against a 5 µm section is a further **220–1000×**. No clinical imaging modality surveyed here can supply a spatially resolved vascular, ECM, oxygen, or IFP **source field** at the scale our PDE operates on. Painting one averaged Ktrans/rCBV/ADC/TBR value across ~10⁴ of our grid cells is not adding spatial information — it is discarding the heterogeneity our model exists to resolve, dressed up as a data upgrade.

**What I got wrong in the first draft, and what it changes.** I sized our domain at ~1.4 mm from a synthetic 120×120 illustration and concluded our whole simulated patch was smaller than one perfusion voxel. That was false. The real sections are 8.1 × 5.09 mm, 11.5 × 10.9 mm, and ≈14.4 × 10.6 mm, spanning **18–125 DCE-MRI voxels** (8–43 at DSC resolution). This does not rescue the source-field idea — the per-edge ratio is untouched — but it materially upgrades the *validation* case: imaging can supply a mm-scale regional map across a specimen, not a single scalar, which converts calibration from an unfalsifiable range check into a real correlation test with a defined null and (on the larger sections) >98% power against the r ≈ 0.5 effect size the field actually achieves.

**What now binds is provenance, not resolution.** The regional-validation experiment cannot be run on our current data at any voxel count: none of our three sections has paired MRI, and piecemeal GBM resection means there is usually no intact specimen geometry to register against an imaging volume — which forces the fallback to stereotactic point biopsies pooled across many patients, an established design (that is how the VSI and ASL validation studies got their r ≈ 0.49) but a prospective multi-patient study rather than an analysis of what we hold. Note also that the flagship section S1, which every other result in this repo is built on, is the *worst* candidate for this test at 8–23 tissue voxels and 23–69% power; S2 and S3 are the viable ones.

The one part of the imaging landscape that *does* resolve at or finer than our grid — micro-CT (1–5 µm), synchrotron phase-contrast (~1.4 µm), and tissue-clearing/light-sheet or CODA-style serial reconstruction (0.5–2 µm) — achieves that only ex vivo or in animal models. That reframes the real opportunity as "a better way to build the next-generation explant/preclinical microenvironment reconstruction," not "unlock the existing clinical imaging archive." If the goal is a source field at our resolution, the honest next step is CODA-style 3D serial-section reconstruction or light-sheet imaging of a resected GBM specimen with an endothelial stain — not pulling Ktrans off a PACS.
