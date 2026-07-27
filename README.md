# tmesim

**A glioblastoma tumor-microenvironment drug-transport simulator built entirely on
real measured data.**

`tmesim` reconstructs a patient's glioblastoma **microenvironment** — not just the
tumor, but the whole niche: malignant cells and their states, stroma, immune
infiltrate, endothelium, ECM, oxygen — from real open spatial-transcriptomics
tissue, then runs mechanistic transport and pharmacodynamics through that
reconstruction to predict where a drug actually goes and what it can kill.

Everything upstream of the prediction is measured. Nothing is synthetic.

> **The one thing this does not do:** predict a real patient's drug response and
> check it. No open dataset pairs spatial GBM tissue with that patient's measured
> response, so the response axis is **unvalidated**. That is a data-access gap, not
> an engineering gap, and it is stated everywhere in this repo rather than buried.
> See [What is validated, and what is not](#what-is-validated-and-what-is-not).

---

## The real data foundation

All open, all downloaded, all traceable to a citation. Provenance for each is in
[`docs/real_data/`](docs/real_data/).

| Layer | Real dataset | Scale |
|---|---|---|
| Single-cell spatial | **10x Xenium** GBM section (CC BY 4.0) | **40,887 cells**, true µm coordinates, 319-gene panel |
| Spot spatial | 4 Visium GBM sections (10x public + **GSE194329**) | 1,775–4,764 spots each, whole transcriptome |
| Cell reference | Darmanis 2017, GEO **GSE84465** | 23,460 genes × 3,589 annotated cells, 7 cell types |
| Niche architecture | **Ivy GAP** (Allen Institute, Puchalski 2018) | 270 samples across 7 anatomic structures |
| Drug potency | **GDSC** glioma + **HGCC** (Johansson 2020, PMID 32668248) | 262 drugs × 85 patient-derived cultures |
| Physics ground truth | 7 **PMID-verified** tissue-PK anchors, 6 drugs | penetration lengths + brain:plasma ratios |

**Synthetic data is gone.** `tmesim/data/synthetic.py` and the fabricated-geometry
loader are deleted from the package. The synthetic-validation era is archived under
[`archive/synthetic_era/`](archive/synthetic_era/README.md), and its honest findings
— including what synthetic validation *could not* establish — are preserved in
[`results/STAGE0_FINDINGS.md`](results/STAGE0_FINDINGS.md) and
[`results/DRUG_SELECTION_RIGOR.md`](results/DRUG_SELECTION_RIGOR.md).

---

## Install

Python ≥ 3.9.

```bash
git clone <this-repo> tme-sim
cd tme-sim
python3 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -e .
```

Pulls numpy, scipy, matplotlib, pandas, scikit-image, networkx, scikit-learn,
requests, tifffile (see `pyproject.toml`).

---

## Quickstart

```bash
# What the CLI offers.
python -m tmesim.cli --help

# The 10 library GBM drugs with structure-derived descriptors (no data needed, ~0.4 s).
python -m tmesim.cli molecules

# Download the auto-acquirable open datasets into ./data_cache/  (~2 min).
python -m tmesim.cli acquire

# Reconstruct the REAL whole-niche microenvironment from the Xenium GBM section
# at the production resolution (~17 s).
python -m tmesim.cli reconstruct --source xenium --dx 12

# Reconstruct a real Visium section (native ~100 µm spot pitch, ~44 s).
python -m tmesim.cli reconstruct --source visium

# Reconstruct + compile a real drug + simulate transport and kill (~1 min 53 s).
python -m tmesim.cli simulate --source xenium --dx 12 --drug temozolomide

# Same, driven from the molecular library (structure -> DrugParams -> engine).
python -m tmesim.cli molecular-run --drug doxorubicin --source xenium --dx 12
```

`reconstruct` prints the niche report (cell-type fractions, mean vascular/ECM/
oxygen, the ECM marker transcripts used, and the per-cell sensitivity model in
force) and writes a figure. `simulate` adds the compiled drug parameters and the
transport→kill result.

The Xenium/Visium/GSE84465/HGCC caches are open but too large for `acquire`; the
exact `curl` commands are in [`docs/REPRODUCE.md`](docs/REPRODUCE.md) §2b.

### Reproduce everything

```bash
./scripts/reproduce_real.sh              # core pipeline, ~1 h 45 m
./scripts/reproduce_real.sh --with-deep  # + niche / combination / dose / patient reports, ~3 h
./scripts/reproduce_real.sh --dry-run    # print every command, run nothing
```

One script: acquire (if missing) → reconstruct → flagship simulation → drug panel
→ cross-tissue → validations. Deterministic (seed 0). Setup, data acquisition,
expected outputs and runtimes: [`docs/REPRODUCE.md`](docs/REPRODUCE.md).

---

## Corrections — read these before quoting any number

These are real errors that were found and fixed. They are listed here rather than
in a footnote because several of them retract previously published headlines.

### 1. Production resolution is dx = 12 µm, and `penetration_index` does not rank drugs

Two separate problems with the same metric:

* **Not converged at coarse grids.** The calibrated perivascular boundary layer is
  ~45 µm thick, so any grid with dx ≳ 45 µm cannot resolve it and reports an
  artificially flat (high) index. On the Xenium section the index falls
  0.396 (dx 30) → 0.311 (20) → 0.204 (15) → **0.196 (12)** → 0.197 (10). The dx = 30 µm
  values used in earlier runs sit ~2× high and are superseded.
* **Scale-invariant even when converged.** `penetration_index = mean(AUC) / p99(AUC)`
  is a *ratio*: multiply the whole exposure field by any constant and it does not
  move. At dx = 12 µm it therefore reports the **vascular geometry**, which is
  identical for all drugs on one section — it collapses to **0.194–0.196 across the
  entire five-drug panel** even though their delivered exposure differs by **~122×**.

  **Compare drugs by mean tissue AUC or therapeutic (IC50) coverage instead.** The
  earlier ranking "temozolomide penetrates best, doxorubicin barely reaches the
  tumor bulk" is **retracted**. See [`results/real_drug_panel.md`](results/real_drug_panel.md)
  and [`results/delivery_metrics.md`](results/delivery_metrics.md).

`kill_fraction`, by contrast, *is* grid-converged (it plateaus by dx ≲ 20 µm and
moves 0.0035 relative between the two finest grids). See
[`results/numerical_rigor.md`](results/numerical_rigor.md) §2.

### 2. The hypoxic fraction is grid-sensitive and only order-of-magnitude

The same Visium tissue reports **100% hypoxic at dx = 100 µm and 0% at dx = 12 µm**;
the Xenium section reported ≈94% on a coarse grid and ≈13% at dx = 12 µm. The
earlier attribution of small predicted kills to hypoxic resistance was wrong — at
~13% hypoxic they are driven mainly by limited delivered exposure relative to these
drugs' IC50s. Treat the hypoxic fraction as an order-of-magnitude figure, not a
measurement.

### 3. α uses the E-FOLD convention — the earlier 309.6 was a unit error

One mechanistic constant, α, is fit **once** to doxorubicin's real perivascular
anchor (Primeau 2005, PMID 16361566) and then held fixed for every other drug. The
anchor is published as a **half-distance** `d_half`; the simulator's penetration
length is an **e-folding length**, and the two are related by

```
d_half = L · ln 2        =>        L = d_half / ln 2 ≈ 1.443 · d_half
```

Fitting α to the raw half-distance skips that conversion and inflates α by exactly
`1/ln²2 = 2.081×`. That is what produced the old **α = 309.6**; the correct value is
**α ≈ 148.6** (`calibrate_alpha()` → 148.608, target e-fold length 64.92 µm). Any
report still quoting 309.6 is stale. See
[`docs/MATHEMATICAL_MODEL.md`](docs/MATHEMATICAL_MODEL.md) §7–8 and
[`results/calibration.md`](results/calibration.md).

### 4. Per-cell sensitivity is drug-specific only where a real marker was measured

`reconstruct.real_pipeline.apply_drug_sensitivity` sets each cell's intrinsic
sensitivity from a **measured determinant transcript on the 319-gene Xenium panel**,
and *only* where the panel actually carries one:

| drug class | determinant | basis |
|---|---|---|
| doxorubicin, etoposide, mitoxantrone (TOP2 poisons) | **TOP2A** | cytotoxicity scales with TOP2A protein level (PMID 18574145) |
| cisplatin, lomustine, carmustine, methotrexate, vincristine | **MKI67** | proxy for cell-cycle dependence |
| **temozolomide**, procarbazine, bevacizumab | **none — uniform 1.0** | TMZ's canonical determinant is **MGMT**, which is *not on this panel* and whose methylation is not measurable from these transcripts |

Uniform 1.0 is a documented no-resistance upper bound, not a measurement. The
default reconstruction is drug-agnostic (uniform 1.0 for every drug), so drug
*rankings* are unaffected by it; **absolute** kill fractions are conditional on it,
and the reconstruction report says so on every run.

---

## What is validated, and what is not

### Validated against real data

* **Reconstruction fidelity.** Xenium marker phenotyping vs Visium NNLS
  deconvolution — two independent platforms and two independent estimators —
  agree at **Pearson 0.997**; endothelial fraction agrees to **0.03%** (9.74% vs
  9.77%). Niche labels are marker-faithful under 10,000-permutation tests
  (PECAM1 in perivascular d = 1.29, p = 1e-4; VEGFA in necrotic d = 1.35, p = 1e-4).
  Deconvolution is stable (±1.0% per-spot std; held-out pseudobulk recovery
  r = 0.935). A real bug was found and fixed here: a marker-argmax classifier
  inflated ENDOTHELIAL to 20.6%; the threshold-gated replacement gives 9.7%,
  corroborated independently by the Visium 9.77%.
  → [`results/reconstruction_fidelity.md`](results/reconstruction_fidelity.md)
* **Transport physics vs real measured tissue PK.** With α frozen from
  doxorubicin, the model **predicts out-of-sample** penetration for mitoxantrone
  (48.2 µm) and topotecan (83.6 µm, a chemically distinct class) — both inside
  their real cited bands. → [`results/pk_multidrug.md`](results/pk_multidrug.md)
* **Numerics.** Mass is conserved to ~1e-13 on the real irregular tissue mask,
  concentrations stay non-negative, the implicit scheme is stable across a 32× dt
  range, and `kill_fraction` is grid-converged.
  → [`results/numerical_rigor.md`](results/numerical_rigor.md)
* **Potency axis corroborated real-vs-real.** GDSC and HGCC — two independent real
  screens, no model in the loop — agree at **Spearman +0.83 (p = 2.4e-5)**. This is
  also *why* the molecular layer uses measured IC50 rather than a structure→potency
  QSAR: that prediction is weak and confounded.
  → [`results/pharmacology_real.md`](results/pharmacology_real.md)

### NOT validated

* **Predicted drug response vs real patient outcomes.** Every open spatial GBM
  dataset is an untreated resection with no paired measured response; every open
  drug-response dataset (GDSC, HGCC) is dissociated cells with no
  microenvironment. **No open dataset bridges the two.** Closing the loop requires
  HTAN GBM spatial multi-omics (dbGaP + HTAN DUA) paired with NCI PDMR / organoid
  biobank screens (MTA). The pipeline is built to ingest exactly that.
* **BBB-penetrance ordering is only partially concordant** with real brain:plasma
  ratios (Spearman 0.50): methotrexate is correctly least penetrant, but
  cisplatin/temozolomide invert — the known limit of structure-only rules with no
  efflux or carrier transport.
* **Absolute kill fractions** are conditional on the sensitivity model above, on
  the single-anchor α, and on 2D geometry. Rankings are far more trustworthy than
  absolutes; ratios more than either.

Full summary: [`results/REAL_VALIDATION.md`](results/REAL_VALIDATION.md).
Every equation, with its provenance: [`docs/MATHEMATICAL_MODEL.md`](docs/MATHEMATICAL_MODEL.md).

---

## Module map

```
data/real/  ─►  reconstruct/real_*  ─►  Microenvironment  ─►  physics/ + cells/  ─►  sim/  ─►  validation/  ─►  viz/
(real bytes)     (real niche)            (core object)         (transport + kill)     (engine)   (vs real PK)      (figures)
                                              ▲
                                    molecular/ (structure → DrugParams)
```

| Package | Responsibility |
|---|---|
| `tmesim.core` | **Frozen** data model: `Domain`, `Microenvironment`, `CellPopulation`, `TissueFields`, `DrugParams`, the config dataclasses and the `*Result` containers. Everything builds against these. |
| `tmesim.data.real` | Loaders for the real caches: `singlecell_spatial` (Xenium), `visium` (4 sections), `scref_gbm` (GSE84465), `ivygap`, `drug_response` (HGCC + GDSC). |
| `tmesim.data` | `acquire.py` — best-effort live download (TCGA-GBM, Ivy GAP, GDSC); the only module permitted network access. `spatial_ingest.py` — the single-sample Visium reader the `real/visium` loader generalises. |
| `tmesim.molecular` | Structure → physics. Verified structures/SMILES (`library`), descriptors, Stokes–Einstein + Brownian + Ogston–Renkin mesh diffusivity, permeability/BBB, measured GDSC potency, the saturable binding-site barrier, and `compile_drug` which emits a `DrugParams`. |
| `tmesim.reconstruct` | Real tissue → real niche. `real_xenium` (threshold-gated marker typing, per-cell drug sensitivity), `real_visium` (NNLS deconvolution vs the Darmanis atlas), `real_fields` (endothelial-KDE vasculature, matrix-transcript ECM), `region_label` (Ivy GAP niche zones), `real_pipeline` (`reconstruct_real`, `apply_drug_sensitivity`, `reconstruction_report`). |
| `tmesim.physics` | `diffusion` (convection–diffusion–reaction transport → `TransportResult`), `oxygen` (steady-state O₂ + hypoxia banding), `pressure` (interstitial fluid pressure and the Jain-effect velocity). |
| `tmesim.cells` | `agents` (oxygen → cell state), `pkpd` (local drug AUC → Hill kill probability, modulated by sensitivity and state). |
| `tmesim.treatment` | `drugs` — drug lookup, plasma PK profile, and the spatial delivery source (`vascular_density × bbb_penetrance × plasma`). |
| `tmesim.sim` | `engine.simulate_response` — the deterministic per-tissue integrator. |
| `tmesim.validation` | The honesty layer. `calibration` (the single α), `pk_anchors` / `pk_multidrug` / `pk_reality` / `glioma_regional_pk` (real measured PK), `reconstruction_fidelity`, `vascular_validation`, `cellstate_validation`, `delivery_metrics` (the defensible replacements for `penetration_index`), `vessel_profile`, `numerical_rigor`, `uncertainty`, `pharmacology_real`. |
| `tmesim.viz` | `plots`, `molecular_plots`, `explainer`, `patient_report`. |
| `tmesim.cli` | argparse entry point: `reconstruct`, `simulate`, `molecules`, `molecular-run`, `acquire`. |
| `experiments/` | Runnable real-data studies: flagship simulation, drug panel, cross-tissue, niche response, combination, dose/schedule, per-patient report. |

Data-flow and physics detail: [`docs/architecture.md`](docs/architecture.md). The
frozen public interface: [`ARCHITECTURE_CONTRACT.md`](ARCHITECTURE_CONTRACT.md);
where the two disagree, the contract wins. The molecular layer's own contract:
[`MOLECULAR_CONTRACT.md`](MOLECULAR_CONTRACT.md).

---

## Why glioblastoma

GBM is the disease where the microenvironment is not a second-order correction but
the crux:

* **Delivery is the binding constraint.** The blood–brain barrier admits only a
  small, spatially heterogeneous fraction of most systemic agents.
* **There is no clean margin.** Malignant cells infiltrate normal brain well past
  any resectable boundary, and the drug has to physically reach them.
* **Transport-limited regions are the resistant ones.** Sparse, tortuous
  vasculature drives hypoxia *and* poor delivery in the same place, while elevated
  interstitial pressure opposes inward flux.

The real panel makes the point concretely: on one real tumor, **doxorubicin is 50×
more potent than cisplatin** (IC50 0.5 vs 24.6 µM) and still loses, because its
delivered exposure is ~120× lower (BBB 0.051, largest MW, consumed within ~45 µm of
a vessel). **Delivery, not potency, decides the outcome** — and that is exactly the
term a non-spatial model has no way to represent.

---

## Honest limitations

* **2D.** A whole-slide section is naturally 2D, but out-of-plane vascular
  connectivity is invisible: a vessel just above the plane makes a point look
  "far" when it is not.
* **One anchor sets the delivery scale.** α is fit to a single doxorubicin
  measurement whose own uncertainty is plausibly ±50–100%. It has the *smallest*
  local elasticity of the tuned knobs (|E_pen| ≈ 0.11) but the widest prior range,
  so it remains a leading global uncertainty.
* **Reconstruction gates dominate the local sensitivity.** The endothelial typing
  threshold (|E_pen| ≈ 0.52) and the KDE bandwidth (|E_kill| ≈ 0.44) move the
  prediction more than the solver does. The limiting factor is the biology and the
  reconstruction assumptions, not the numerics.
* **Fields are inferred, not measured.** Vascular density, ECM and oxygen are
  reconstructed from real markers and real perfusion signal — faithful to the data,
  but model fields. `vascular_density` is max-normalised, so it carries geometry,
  not an absolute perfused fraction.
* **Reduced biology.** No adaptive resistance, no immune dynamics, no efflux or
  carrier transport, and systemic PK is a simplified sustained plasma level.
* **Research prototype.** Not a clinical decision tool. No diagnostic, regulatory
  or treatment-selection use.

---

## Where to look next

| Question | File |
|---|---|
| How do I reproduce all of this? | [`docs/REPRODUCE.md`](docs/REPRODUCE.md) |
| What is the model, mathematically? | [`docs/MATHEMATICAL_MODEL.md`](docs/MATHEMATICAL_MODEL.md) |
| How do the pieces fit together? | [`docs/architecture.md`](docs/architecture.md) |
| What has actually been validated? | [`results/REAL_VALIDATION.md`](results/REAL_VALIDATION.md) |
| Where did each dataset come from? | [`docs/real_data/`](docs/real_data/) |
| Is the solver trustworthy? | [`results/numerical_rigor.md`](results/numerical_rigor.md) |
| How should I compare drugs? | [`results/delivery_metrics.md`](results/delivery_metrics.md) |
| What did the synthetic era show? | [`archive/synthetic_era/README.md`](archive/synthetic_era/README.md) |
