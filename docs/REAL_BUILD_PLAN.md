# The Real Build — Patient Microenvironment Simulation (no synthetic)

## Principle (hard rules)
1. **No synthetic ground truth, no invented patients.** Every input traces to a real measurement. The synthetic cohort generator and all synthetic-ground-truth experiments are RETIRED from the pipeline.
2. **Whole microenvironment, not just tumor.** Model the niche: malignant cells + their states, astrocytes/oligodendrocytes, CAF/reactive stroma, tumor-associated macrophages & T cells, endothelium/vasculature, ECM, and the peritumoral / infiltrative-margin / necrotic regions.
3. **Rigorous & honest.** Validate reconstruction against held-out real data; validate physics against real measured tissue pharmacokinetics; be explicit about what is genuinely closed vs. access-gated.

## The one honest constraint
Building the microenvironment from **real** data is fully achievable now from open sources (spatial transcriptomics, single-cell atlases, Ivy GAP, TCGA). The **gold-standard validation** — real per-patient tissue paired with that patient's *measured drug response* — is largely behind controlled access (HTAN/dbGaP + patient-derived-organoid screens). We will (a) get every open real dataset, (b) use real organoid/PDX drug-response where openly available, and (c) document the exact access path for the rest. No synthetic stand-ins.

## Layers

### L1 — Real data (foundation)
- **Spatial architecture**: 10x Visium GBM (open), 10x Xenium single-cell-resolution GBM (if open), Ivy GAP anatomic structures (leading edge, cellular tumor, microvascular proliferation, perinecrotic zone).
- **Single-cell reference atlas**: an open GBM scRNA-seq atlas (e.g. Neftel 2019 / GSE131928-class) to deconvolve spatial spots into real cell-type fractions and malignant states.
- **Vasculature / hypoxia / ECM signals**: from real markers (endothelial PECAM1/CD31; hypoxia CA9/VEGFA/HIF targets; ECM collagen/fibronectin programs) mapped spatially.
- **Drug molecules**: real structures + measured GDSC IC50 (already built — keep).
- **Outcomes for validation**: real tissue-PK penetration anchors (multiple drugs); open patient-derived-organoid / PDX GBM drug screens where available.

### L2 — Reconstruction (real sample -> real microenvironment)
- Spot/cell deconvolution -> real cell-type composition per location.
- Place real cell populations at real coordinates; build real vascular network, ECM density, and an oxygen field inferred from real perfusion/hypoxia signal.
- Whole-niche assembler including margin/normal/necrotic regions.

### L3 — Simulation (keep the physics engine, extend to the niche)
- Real physics engine (convection-diffusion-reaction transport, oxygen, interstitial pressure, Hill kill) — keep. **Not "already validated": that phrasing stood here and is retracted.** The solver's input coefficients are externally testable and the engine currently misses them — the MW-only diffusivity law is falsified against 71 measured compounds, and the multicellular-layer `D_eff` and hindrance anchors both FAIL. What is calibrated is one penetration anchor (α, one doxorubicin measurement).
- Extend response to niche cell types (immune/stromal modulation) mechanistically, parameterized from real biology only.

### L4 — Validation (rigorous, real)
- Reconstruction fidelity vs held-out real spatial data.
- Physics vs multiple real measured tissue-PK anchors (not just doxorubicin).
- Drug response vs real organoid/PDX data where obtainable; else documented access gap.

## Concurrent agent plan
**Wave 1 — REAL DATA ACQUISITION (launching now, parallel):** hunt + download + characterize each real dataset; report exact format + what is open vs gated.
1. GBM single-cell reference atlas (deconvolution reference)
2. 10x Visium GBM sections (multiple, open)
3. 10x Xenium / single-cell spatial GBM (open, if any)
4. Ivy GAP full anatomic + expression
5. Real drug-response outcomes (organoid/PDX/HCMI/GEO) + GDSC
6. Real tissue-PK penetration compendium (multi-drug, literature)

**Wave 2 — REAL RECONSTRUCTION (after data shapes known):** deconvolution; vasculature/ECM/oxygen from real markers; whole-niche assembler.

**Wave 3 — SIMULATION + VALIDATION:** niche-extended engine; reconstruction validation; multi-drug PK validation; response validation.

**Wave 4 — INTEGRATION:** remove synthetic modules from the pipeline; real-only CLI/tests/docs.

Wave 1 gates everything — we build the real pipeline against the *real* data formats, not assumptions.
