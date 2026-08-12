# tmesim

**A glioblastoma tumour-microenvironment drug-transport simulator built entirely on
real measured data — and adversarially audited.**

`tmesim` reconstructs a patient's glioblastoma **microenvironment** — not just the
tumour, but the whole niche: malignant cells and their states, stroma, immune
infiltrate, endothelium, ECM, oxygen — from real open spatial-transcriptomics
tissue, then runs mechanistic transport and pharmacodynamics through that
reconstruction to predict where a drug actually goes and what it can kill.

Everything upstream of the prediction is measured. Nothing is synthetic.

> ### Read this before you read anything else
>
> A rigour audit ([`results/SESSION_AUDIT.md`](results/SESSION_AUDIT.md)) attacked
> this model and found **nine distinct defects**, several of them mutually
> cancelling. Every headline number published before that audit was wrong in at
> least one way. The engine and the reconstruction were sound; the *claims* were
> not.
>
> `SESSION_AUDIT.md` is the authoritative statement of current state — which
> defects are fixed, which are diagnosed but outstanding, and which results are
> retracted. **This README summarises it; where the two differ, the audit wins.**

---

## The real data foundation

All open, all downloaded, all traceable to a citation. Per-dataset provenance is
in [`docs/real_data/`](docs/real_data/); the complete re-download record for the
gitignored cache is [`docs/DATA_INVENTORY.md`](docs/DATA_INVENTORY.md).

| Layer | Real dataset | Scale |
|---|---|---|
| Single-cell spatial (primary) | **10x Xenium** GBM FFPE, Human Brain 319-gene panel (CC BY 4.0) | **40,887 cells**, true µm coordinates, 8.1 × 5.1 mm |
| Single-cell spatial (deep) | **10x Xenium** GBM FFPE, Immuno-Oncology **480-gene** panel (CC BY 4.0) | **816,769 cells**, 11.5 × 10.8 mm — 20× the cells, 3× the area |
| Spot spatial | 4 Visium GBM sections (10x public parent + **GSE194329** GBM1/2/3) | 1,775–4,764 spots each, whole transcriptome |
| Cell reference | Darmanis 2017, GEO **GSE84465** | 23,460 genes × 3,589 annotated cells, 7 cell types |
| Niche architecture | **Ivy GAP** (Allen Institute, Puchalski 2018) | 270 RNA-seq samples across 7 anatomic structures |
| Drug potency | **GDSC2** glioma + **HGCC** (Johansson 2020, PMID 32668248) | 262 drugs × 85 patient-derived cultures |
| Physics ground truth | **PMID-verified** tissue-PK anchors, 6 drugs | penetration lengths + brain:plasma ratios |
| External validation corpus | the **data hunt** caches under `data_cache/hunt/` (~5.3 GB, 40 directories) | 164 datasets surveyed, 143 open — see [`results/DATA_HUNT.md`](results/DATA_HUNT.md) |

The two Xenium panels are **complementary, not redundant**: they share only 60
genes. The brain panel carries PECAM1/SOX2/OLIG2; the IO panel loses those but
gains CD3E, FOXP3, PDCD1, CD274, **VEGFA** and **HIF1A**. Marker maps are selected
by panel content (`reconstruct.marker_panels`), and the 319-gene behaviour is
unchanged bit-for-bit.

**Synthetic data is gone.** `tmesim/data/synthetic.py` and the fabricated-geometry
loader are deleted from the package; the reproducibility harness fails loudly if
either reappears. The synthetic-validation era is archived under
[`archive/synthetic_era/`](archive/synthetic_era/README.md), with its honest
findings — including what synthetic validation *could not* establish — preserved in
[`results/STAGE0_FINDINGS.md`](results/STAGE0_FINDINGS.md) and
[`results/DRUG_SELECTION_RIGOR.md`](results/DRUG_SELECTION_RIGOR.md).

---

## Install

Python ≥ 3.9 (tested on CPython 3.9.6, macOS Apple Silicon).

```bash
git clone <this-repo> tme-sim
cd tme-sim
python3 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -e .
pip install rdkit openpyxl statsmodels xlrd pyxlsb rdata   # see docs/REPRODUCE.md §1
```

`pip install -e .` pulls numpy, scipy, matplotlib, pandas, scikit-image, networkx,
scikit-learn, requests and tifffile. The second line is **not** in
`pyproject.toml` and is needed for the molecular descriptors (rdkit) and the
external-validation readers; see [`docs/REPRODUCE.md`](docs/REPRODUCE.md) §1.

---

## Quickstart

Every command below was executed in this repo before being documented; the
runtimes are measured, not estimated. New to this repo?
[`docs/QUICKSTART.md`](docs/QUICKSTART.md) is a shorter, install-to-first-result
walkthrough (CLI + Python API), aimed at someone who has never seen this
codebase's module layout.

`pip install -e .` also registers a `vivum` console script — `vivum --help` is
the exact same tool as `python -m tmesim.cli --help` below, just installed
rather than run as a module; either form works everywhere on this page.

```bash
# What the CLI offers.
python -m tmesim.cli --help

# The 10 library GBM drugs with structure-derived descriptors (no data, ~0.5 s).
python -m tmesim.cli molecules

# Compare several drugs on ONE reconstructed section (new: compare-drugs).
python -m tmesim.cli compare-drugs --drugs temozolomide,cisplatin,doxorubicin

# Launch the Studio web UI (needs `pip install -e '.[app]'`).
python -m tmesim.cli studio

# Download the auto-acquirable open datasets into ./data_cache/  (~2 min, ~90 MB).
python -m tmesim.cli acquire

# Reconstruct the REAL whole-niche microenvironment from the 319-panel Xenium GBM
# section at the PRODUCTION resolution (~6 s, 40,887 cells).
python -m tmesim.cli reconstruct --source xenium --dx 12

# Same pipeline on the 480-gene Immuno-Oncology section (~45 s, 816,769 cells).
python -m tmesim.cli reconstruct --source xenium --sample data_cache/xenium_gbm_io --dx 12

# Reconstruct a real Visium section (native ~100 µm spot pitch, ~44 s).
python -m tmesim.cli reconstruct --source visium

# Reconstruct + compile a real drug + simulate transport and kill (~2 min).
python -m tmesim.cli simulate --source xenium --dx 12 --drug temozolomide

# Same, driven from the molecular library (structure -> DrugParams -> engine).
python -m tmesim.cli molecular-run --drug doxorubicin --source xenium --dx 12

# Metric known-answer checks, no simulation (~1 s). Proves penetration_index is
# scale-invariant and therapeutic_coverage is not.
python -m tmesim.validation.delivery_metrics --selftest
```

`reconstruct` prints the niche report (cell-type fractions, mean vascular/ECM/
oxygen, the ECM marker transcripts used, the tissue-mask mode, and the per-cell
sensitivity model in force) and writes a figure. `simulate` adds the compiled drug
parameters and the transport→kill result.

The Xenium, Visium, GSE84465 and HGCC caches are open but too large for
`acquire`; the exact `curl` commands are in
[`docs/REPRODUCE.md`](docs/REPRODUCE.md) §2b, and the full backup record —
including every `data_cache/hunt/` directory — is
[`docs/DATA_INVENTORY.md`](docs/DATA_INVENTORY.md).

### Reproduce everything

```bash
./scripts/reproduce_real.sh              # core pipeline
./scripts/reproduce_real.sh --with-deep  # + niche / combination / dose / patient reports
./scripts/reproduce_real.sh --list       # print the stages and exit
./scripts/reproduce_real.sh --dry-run    # print every command, run nothing
./scripts/reproduce_real.sh --only panel # run one stage
```

One script: preflight → data (acquire if missing) → reconstruct → flagship
simulation → drug panel → cross-tissue → validation → external validation.
Deterministic (seed 0). Setup, data acquisition, expected outputs and runtimes:
[`docs/REPRODUCE.md`](docs/REPRODUCE.md).

---

## Five corrections — read these before quoting any number

These are real errors that were found and fixed, or real confounds that were
found and are still open. They are here rather than in a footnote because several
of them **retract previously published headlines**.

### 1. Production resolution is dx = 12 µm, and `penetration_index` is RETIRED

Two separate problems with the same metric.

* **Not converged at coarse grids.** The calibrated perivascular boundary layer is
  ~45 µm thick, so any grid with dx ≳ 45 µm cannot resolve it and reports an
  artificially flat (high) index. Values quoted at dx = 30 µm or dx = 100 µm are
  superseded.
* **Scale-invariant even when converged — it cannot discriminate drugs at all.**
  `penetration_index = mean(AUC) / p99(AUC)` is a *ratio*: multiply the whole
  exposure field by any constant and it does not move. At dx = 12 µm it therefore
  reports the **vascular geometry**, which is identical for every drug on one
  section. Across the five-drug panel it spans **0.1946–0.1959 — a spread of
  0.0014 — while those same drugs' delivered exposure differs by 120×.** Its grid
  drift (0.0183) is ~13× larger than its entire across-drug signal.

  **Use `therapeutic_coverage` at 1× IC50 instead** — the fraction of tumour cells
  whose time-averaged exposure `E = AUC/T` exceeds the drug's measured IC50. It is
  dimensionally sound, grid-converged (±0.03), and spans **0.311–0.925** on the
  same panel (separation 21.9 vs 0.08). Quote it with the 0.5×/2× IC50 band
  beside it and `median_tumor_auc` (µM·s) as the supporting absolute.

  The earlier ranking "temozolomide penetrates best, doxorubicin barely reaches
  the tumour bulk" is **retracted**, as is the cross-tissue "kill spans 2% → 24%"
  headline (that spread mixed dx = 12 µm with dx = 100 µm runs).
  → [`results/delivery_metrics.md`](results/delivery_metrics.md),
  [`results/real_drug_panel.md`](results/real_drug_panel.md)

`kill_fraction`, by contrast, **is** grid-converged (plateaus by dx ≲ 20 µm).

### 2. The hypoxic fraction is grid-sensitive and only order-of-magnitude

The Xenium section reported ≈94% hypoxic at dx = 30 µm, ≈13% at dx = 12 µm and
≈7% at dx = 8 µm — not converged; the same Visium tissue reports 100% at
dx = 100 µm and 0% at dx = 12 µm. Two grid-dependence bugs (a raw-count oxygen
sink scaling as dx², and a grid-dependent tissue mask) were **partially cancelling
each other**, which made the output look more stable than it was. Both are now
fixed — the tissue mask is `geometric` (the union of 30 µm physical discs around
measured centroids, so the domain does not move with `dx`) — and the earlier
attribution of small predicted kills to hypoxic resistance was **wrong**: they are
driven mainly by limited delivered exposure relative to these drugs' IC50s. Treat
any hypoxic fraction as an order-of-magnitude figure, not a measurement.

### 3. α uses the E-FOLD convention — the earlier 309.6 was a unit error

One mechanistic constant, α, is fit **once** to doxorubicin's real perivascular
anchor (Primeau 2005, PMID 16361566) and then held fixed for every other drug and
every tissue. The anchor is published as a **half-distance** `d_half`; the
simulator's penetration length is an **e-folding length**:

```
d_half = L · ln 2        =>        L = d_half / ln 2 ≈ 1.443 · d_half
```

Fitting α to the raw half-distance skips that conversion and inflates α by exactly
`1/ln²2 = 2.081×`. That produced the old **α = 309.6**; the correct value is
**α ≈ 148.6** (`calibrate_alpha()` → 148.608, target e-fold length 64.92 µm).
`calibrate_alpha` now takes an e-folding target and a regression test guards it, so
the mistake cannot recur silently. **Any report still quoting 309.6 is stale, and
absolute µM·s figures from before the fix are not comparable to current ones.**
→ [`docs/MATHEMATICAL_MODEL.md`](docs/MATHEMATICAL_MODEL.md) §7–8,
[`results/calibration.md`](results/calibration.md)

### 4. Per-cell sensitivity is drug-specific from real markers — except for temozolomide

`reconstruct.real_pipeline.apply_drug_sensitivity` sets each cell's intrinsic
sensitivity from a **measured determinant transcript**, and only where the panel
actually carries one:

| drug class | determinant | basis |
|---|---|---|
| doxorubicin, etoposide, mitoxantrone (TOP2 poisons) | **TOP2A** | cytotoxicity scales with TOP2A level (PMID 18574145) |
| cisplatin, lomustine, carmustine, methotrexate, vincristine | **MKI67** | proxy for cell-cycle dependence |
| **temozolomide**, procarbazine, bevacizumab | **none — uniform 1.0** | TMZ's canonical determinant is **MGMT**, which is absent from **both** Xenium panels (319-gene *and* 480-gene), and methylation is not measurable from these transcripts at all |

Uniform 1.0 is a documented no-resistance upper bound, not a measurement.
MGMT's absence is not merely unmodelled but **unmodellable from this modality** —
so a determinant-mode cross-drug comparison quietly advantages temozolomide,
because its resistance determinant is the one that is *unmeasured*. That unfairness
is a **diagnosed, unfixed defect** (audit #8). The reconstruction default is
therefore **uniform 1.0 for every drug**, which is the only fair cross-drug mode:
sensitivity cancels out of any ranking, and absolute kills are an equal upper bound.

### 5. The vascular field is a smooth KDE, not vessels — this confounds every vessel-relative claim

`vascular_density` is a **Gaussian KDE (bandwidth ~40 µm) of measured endothelial
cell positions**, max-normalised. It is the source term for every drug
(`delivery = vascular_density × bbb_penetrance × plasma`), and being a blur rather
than a vessel network it produces three converging artifacts:

* fitted penetration length is **26–28 µm for all five drugs** against a
  **source-null of 23.7 µm** — the decay is measuring KDE bandwidth, not drug
  transport ([`results/delivery_metrics.md`](results/delivery_metrics.md));
* the exposure profile has a **flat floor at ~12% of peak**, so no tissue is truly
  drug-starved and delivery is never the binding constraint
  ([`results/dose_schedule.md`](results/dose_schedule.md));
* "perivascular regions get more drug" is **circular** — the region label and the
  drug source share the same endothelial signal
  ([`results/niche_response.md`](results/niche_response.md)).

Independently, `validation.vascular_validation` finds the arrangement is *nearly
structureless*: median tumour-to-vessel distance is within 1% of a complete
spatial randomness null at the same density, whereas real GBM is strongly
clustered. The model therefore **under-represents the poorly-delivered tail**.

**Treat every vessel-relative spatial claim in this repo as confounded.**
Replacing the KDE with a discrete segmented vessel network (literature basis:
~20 endothelial cells per mm of vessel, 1–2 EC per capillary cross-section —
Shapson-Coe 2024 *Science*, Sargent 2023) is the single highest-value structural
fix outstanding. `tmesim/reconstruct/vessel_network.py` exists as the beginning of
that work.

---

## What is validated, and what is not

| Axis | Status | Evidence |
|---|---|---|
| **Reconstruction fidelity** | **MIXED — one claim retracted, one recovered, two still fail** | Originally reported as VALIDATED (Xenium vs Visium, Pearson **0.997**). Re-run after a classifier-arity fix (`4ad221f`) found the number that mattered — Pearson r on the harder non-tumour subvector — collapsed **0.945 → 0.016** (rank agreement ρ = 1.00 → 0.40): **RETRACTED**, → [`results/revalidation.md`](results/revalidation.md). Separately, a held-out second patient (`results/heldout_patient.md`) originally found the endothelial-fraction-vs-human-volume-EM agreement did NOT replicate (4.18% vs the measured 8.0%) — but net of later fixes (discrete-vessel promotion, the same classifier-arity fix) it now **DOES replicate** when re-run (7.07% vs 8.0%), a genuine improvement not yet reflected in `results/SESSION_AUDIT.md`. The SAME held-out test found the **vascular area fraction** claim still fails badly (10.98% vs measured 2.1%, 5.2×) and **MVD** is not reproduced on either section — both real, current, unresolved gaps. → [`results/heldout_patient.md`](results/heldout_patient.md) is the most current single source for all of this. |
| **Numerics** | **VALIDATED** | Mass conserved to ~1e-13 on the real irregular mask; concentrations never negative; dt-converged across a 32× range; `kill_fraction` grid-converged. Production dx = 12 µm is also the finest grid satisfying the advection CFL at dt = 10 s. → [`results/numerical_rigor.md`](results/numerical_rigor.md) |
| **Transport physics vs measured tissue PK** | **VALIDATED out-of-sample** | With α frozen from doxorubicin, the model predicts mitoxantrone (48.2 µm) and chemically-distinct topotecan (83.6 µm), both inside their real cited bands. → [`results/pk_multidrug.md`](results/pk_multidrug.md) |
| **Potency** | **CORROBORATED real-vs-real** | GDSC and HGCC — two independent real screens, no model in the loop — agree at Spearman **+0.83** (p = 2.4e-5). This is *why* the molecular layer uses measured IC50 rather than a structure→potency QSAR. → [`results/pharmacology_real.md`](results/pharmacology_real.md) |
| **Delivery / partition magnitude, in human glioma** | **VALIDATED, with an error budget** | Scored against **125 human glioma tissue samples across 42 patients** plus 27 PMID-anchored regional PK measurements: 5/6 drugs within 10×, 4/6 within 3×; aggregate bias 2.02×, but **1.01× — essentially unbiased — excluding one outlier**. Cross-drug ordering: measured penetrant-vs-excluded separation **4.66×**, model predicts **3.65×** from structure alone. → [`results/glioma_regional_pk.md`](results/glioma_regional_pk.md) |
| **Regional (CE vs non-enhancing) ordering** | **FAILS** | Measured contrast-enhancing vs non-enhancing separation is **1.45×–4.26×** (p < 0.0001, within-patient permutation); the model predicts exactly **1.00×**, because there is **no barrier-integrity term anywhere in it**. The validation derived the fix from the data: a multiplicative barrier cannot work (3.2× spread across drugs) but an **additive leak is drug-independent to within 1.3×** (L = 0.211 / 0.230 / 0.272 across three independent studies), giving `source = vascular_density · (bbb_penetrance + L · bbb_disruption) · plasma`, **L ≈ 0.24** — the signature of size-nonselective paracellular opening. That term now exists in `tmesim/reconstruct/barrier.py` and is reachable via `treatment.drugs.delivery_source(..., bbb_disruption=...)`, but it is **opt-in and OFF by default** — passing `None` reproduces the pre-fix source exactly, so every shipped result still predicts 1.00×. |
| **Spatial / vessel-relative claims** | **CONFOUNDED** | By our own vascular smoothing — correction 5 above. |
| **The vascular delivery law itself** | **DEFECT FOUND; CORRECTED FORM IMPLEMENTED, OPT-IN, OFF BY DEFAULT** | The shipped source is a **one-way leak**, `S = PS·Cp`, with no path back from interstitium to blood, so steady-state tissue concentration is unbounded above the plasma level — and on the GBM flagship **5 of 11 drugs are predicted above plasma** (cisplatin **2.71×**, carmustine 2.54×, lomustine 1.39×, procarbazine 1.23×, temozolomide 1.14×), which passive delivery cannot do. The gradient-driven form `S = PS·(Cp − C)` is the same source minus a first-order loss at the same rate, needs **no new parameter**, and reproduces its closed-form steady state to **6.3e-15**; it is reachable as `PDEConfig(transvascular_exchange=True)` and is **off by default**, so every number shipped in this repo is still the one-way arm and is bit-for-bit unchanged. Turning it on moves absolute kills (temozolomide **−56%**, cisplatin **−77%**, doxorubicin **−3.5%**) **and changes the kill ranking** — see the two ranking bullets below. Separately, the PS field carries an implicit, uncalibrated rate scale of exactly 1 s⁻¹ (`tmesim/treatment/drugs.py::PS_SCALE_S_INV`, default **1.0**, calibrated value **0.16122 1/s**), so every exchange-on number anywhere in this repo is conditional on that scale as well. → [`results/transvascular_exchange_revalidation.md`](results/transvascular_exchange_revalidation.md) |
| **Passive transport as a sufficient explanation of tissue partition** | **FALSIFIED, diagnostically** | On 6 matched human GIST tumour/normal imatinib MSI pairs the transport-only model gives 7.9× geometric-mean fold error vs 13.0× for the no-partitioning null — it beats the null by only ~1.6× and 4 of 5 liver pairs lie *below* the model's structural floor, so **no penetration length can reproduce them**. The missing physics is active hepatic uptake and lysosomal ion-trapping, neither of which is a diffusion process. → [`results/gist_imatinib.md`](results/gist_imatinib.md) |
| **Predicted drug RESPONSE vs real patient outcomes** | **CANNOT BE VALIDATED with open data** | A 39-agent, 164-dataset sweep established that **no open dataset anywhere pairs measured spatial drug distribution with measured treatment response in the same tumour** — not one — and **no glioma + temozolomide spatial drug map exists at all** (METASPACE: zero). The experiment that would close the loop (the Agar/White adavosertib GBM platform) has **no accession**; it is a months-long collaboration request. → [`results/DATA_HUNT.md`](results/DATA_HUNT.md) |

### Also honest about

* **`bbb_penetrance` is a structure-only score with no efflux or carrier
  transport**, so it is only partially concordant with measured brain:plasma
  ratios. The audit's "a mean cannot express a veto" defect **has since been
  fixed**: crossing the barrier is a *conjunction* (small enough AND lipophilic
  enough AND not too polar AND not too H-bonding), so the four desirabilities are
  now combined by a **geometric** mean, where any single near-zero desirability
  drives the product to zero. It adds no tunable parameter — for equal
  desirabilities the two means coincide. On the 6-drug human-glioma scorecard the
  worst-case error falls from **65.6× to 6.9×** (the eribulin outlier: 0.334 →
  0.035 against a measured 0.0051) and all 6 drugs move inside the
  order-of-magnitude band. It makes two drugs *worse* — methotrexate and cefazolin
  carry disqualifying descriptors too — and the module says so.
* **The top-1 drug call is not decision-grade, and it does not survive a change
  of delivery law.** On the shipped one-way source, cisplatin (0.241) and
  doxorubicin (0.231) differ by ~4%, well inside plausible parameter
  uncertainty. Under gradient-driven exchange the pair **reverses** — cisplatin
  0.2311 vs doxorubicin 0.2547 at dx = 12 µm, and 0.2335 vs 0.2547 at dx = 8 µm,
  so it is not a mesh artefact — and on the full 11-drug GBM panel top-1 moves
  **cisplatin → vincristine** with the top-3 set changing too. The Monte-Carlo
  ranking panel that certified this call goes **P(cisplatin rank 1) 0.9277 →
  0.3965** and P(doxorubicin rank 1) 0.0430 → 0.5586, which destroys that
  module's own pre-declared decision-grade bar. Report the top of the panel as a
  **set**, not a winner. →
  [`results/transvascular_exchange_revalidation.md`](results/transvascular_exchange_revalidation.md)
  §4.1, [`results/ranking_robustness.md`](results/ranking_robustness.md)
* **"Ranking is invariant" is false as written, and a ranking claim has to name
  the quantity it is about.** DELIVERY-rank invariance does not imply KILL-rank
  invariance: on the breast panel the AUC ranking is *perfectly* preserved
  across the two delivery laws (Spearman ρ = **1.000**) while the kill ranking
  scrambles (ρ = **0.830**), because the Hill map saturates — drugs on the
  plateau lose almost no kill when AUC drops, drugs on the steep limb lose most
  of theirs. On GBM the kill-rank correlation is high (ρ = +0.9273, Kendall
  +0.8222) and the top-1 still changes; 4 of 45 pairs swap, two of them across
  gaps of 40.5% and 35.6% of the leader's kill against a Bernoulli sd of ~0.003.
  No constant absorbs this: the attenuation is
  `(clearance+binding)/(PS+clearance+binding)` and PS is driven by
  `bbb_penetrance`, which spans **47×** across the GBM panel, so it is a
  **differential re-weighting by permeability**, not a rescale. The ranking's
  historical robustness (to α, to the 3.36× KDE→discrete source change, to the
  cell-state route) does not generalise, because those perturbations were
  common-mode and this one is drug-specific. →
  [`results/transvascular_exchange_revalidation.md`](results/transvascular_exchange_revalidation.md)
  §4.1
* **Absolute kill fractions** are conditional on the sensitivity model, on the
  single-anchor α, on 2D geometry, and on the delivery law (temozolomide −56%,
  cisplatin −77%, doxorubicin −3.5% under exchange). Ratios and same-tissue
  comparisons are more trustworthy than absolutes — but see the two bullets
  above for exactly where a ranking still breaks.

Authoritative current state: [`results/SESSION_AUDIT.md`](results/SESSION_AUDIT.md).
Every equation with its provenance:
[`docs/MATHEMATICAL_MODEL.md`](docs/MATHEMATICAL_MODEL.md).

---

## Module map

```
data/real/  ─►  reconstruct/real_*  ─►  Microenvironment  ─►  physics/ + cells/  ─►  sim/  ─►  validation/  ─►  viz/
(real bytes)     (real niche)            (core object)         (transport + kill)     (engine)   (vs real data)   (figures)
                                              ▲
                                    molecular/ (structure → DrugParams)
```

| Package | Responsibility |
|---|---|
| `tmesim.core` | **Frozen** data model: `Domain`, `Microenvironment`, `CellPopulation`, `TissueFields`, `DrugParams`, the config dataclasses and the `*Result` containers. Everything builds against these. |
| `tmesim.data.real` | Loaders for the real caches: `singlecell_spatial` (both Xenium panels — select with `load(path=...)`), `visium` (4 sections), `scref_gbm` (GSE84465), `ivygap`, `drug_response` (HGCC + GDSC2). |
| `tmesim.data` | `acquire.py` — best-effort live download (TCGA-GBM clinical + expression, Ivy GAP, GDSC2 glioma); the only module permitted network access. `spatial_ingest.py` — the single-sample Visium reader that `real/visium` generalises. |
| `tmesim.molecular` | Structure → physics. Verified structures/SMILES (`library`), RDKit-or-fallback `descriptors`, Stokes–Einstein radius, Brownian + Ogston–Renkin mesh hindrance, `permeability`/BBB, measured GDSC `potency`, the saturable `binding_site_barrier`, and `compile_drug` which emits a `DrugParams`. |
| `tmesim.reconstruct` | Real tissue → real niche. `real_xenium` (threshold-gated marker typing, geometric tissue mask, per-cell sensitivity), `marker_panels` (panel-aware marker selection: `brain_319` vs `io_480`), `real_visium` (NNLS deconvolution vs the Darmanis atlas), `deconvolve`, `real_fields` (endothelial-KDE vasculature, matrix-transcript ECM, oxygen hand-off), `region_label` (Ivy GAP niche zones), `vessel_network` (the discrete-vessel replacement for the KDE, correction 5), `real_pipeline` (`reconstruct_real`, `apply_drug_sensitivity`, `reconstruction_report`). |
| `tmesim.physics` | `diffusion` (convection–diffusion–reaction transport → `TransportResult`), `oxygen` (steady-state O₂ + hypoxia banding), `pressure` (interstitial fluid pressure and the Jain-effect velocity). |
| `tmesim.cells` | `agents` (oxygen → cell state), `pkpd` (local drug AUC → Hill kill probability, modulated by sensitivity and state). |
| `tmesim.treatment` | `drugs` — drug lookup, plasma PK profile, and the spatial delivery source (`vascular_density × bbb_penetrance × plasma`, a one-way leak by default; `exchange_coefficient` supplies the gradient-driven `PS·(Cp − C)` form used when `PDEConfig(transvascular_exchange=True)`), plus `PS_SCALE_S_INV`, the transvascular rate scale that is implicitly 1 s⁻¹ today. |
| `tmesim.sim` | `engine.simulate_response` — the deterministic per-tissue integrator. |
| `tmesim.validation` | The honesty layer, and where this project spends most of its code. **Internal:** `calibration` (the single α), `numerical_rigor`, `uncertainty`, `reconstruction_fidelity`, `vascular_validation`, `cellstate_validation`, `delivery_metrics` (the defensible replacements for `penetration_index`), `vessel_profile`, `ranking_robustness`. **Against real external measurements:** `pk_anchors` / `pk_multidrug` / `pk_reality` (tissue PK), `glioma_regional_pk` (human glioma partition), `gist_imatinib` (human GIST MSI), `gbm_tmz_exvivo` (GSE232173 — glioma + TMZ + measured ex-vivo response), `moncayo_explant` (HGSOC MALDI-MSI + GeoMx), `pharmacology_real` (GDSC vs HGCC). **Under construction:** `cellstate_validation`, `ranking_robustness`, `vessel_network_impact` — these have no entry point yet and produce nothing; see [`docs/REPRODUCE.md`](docs/REPRODUCE.md) §3. **Retired:** `experiment`, `baselines`, `metrics` — the synthetic-cohort harness. |
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

The real panel makes the point concretely. Ranking by **AUC/IC50** reproduces the
simulated kill ordering exactly: delivery spans **120×**, potency spans **1212×**,
and the *product* decides. Doxorubicin is **49× more potent** than cisplatin, is
**120× less delivered**, and — on the shipped one-way delivery law — loses. That
trade-off is exactly the term a non-spatial model has no way to represent.

It is also exactly why the *outcome* of that contest is not a fixed fact. The
delivery axis is the fragile one: correcting the delivery law compresses the GBM
panel's delivery spread from **386× to 90×** by re-weighting every drug by its
own permeability, and this particular contest reverses (cisplatin 0.2311 vs
doxorubicin 0.2547). The mechanism — delivery × potency, resolved in space — is
the durable claim; the winner is not. →
[`results/transvascular_exchange_revalidation.md`](results/transvascular_exchange_revalidation.md)
§4.1

---

## Honest limitations

* **The response axis is unvalidated and cannot currently be validated.** See the
  table above. This is a data-access gap, not an engineering gap.
* **The vasculature is a blur, not a vessel network.** Correction 5. This is the
  largest structural error in the model.
* **No barrier-integrity term.** Correction table, regional axis. The fix is
  specified (L ≈ 0.24) and not implemented.
* **2D.** A whole-slide section is naturally 2D, but out-of-plane vascular
  connectivity is invisible: a vessel just above the plane makes a point look
  "far" when it is not. Distances are Euclidean, not geodesic.
* **One anchor sets the delivery scale.** α is fit to a single doxorubicin
  measurement whose own uncertainty is plausibly ±50–100%. It has the *smallest*
  local elasticity of the tuned knobs (|E_pen| ≈ 0.11) but the widest prior range,
  so it remains a leading global uncertainty.
* **Reconstruction gates dominate the local sensitivity.** The endothelial typing
  threshold (|E_pen| ≈ 0.52) and the KDE bandwidth (|E_kill| ≈ 0.44) move the
  prediction more than the solver does. The biology and the reconstruction
  assumptions are the limiting factor, not the numerics.
* **Fields are inferred, not measured.** Vascular density, ECM and oxygen are
  reconstructed from real markers and real perfusion signal — faithful to the data,
  but model fields. `vascular_density` is max-normalised, so it carries geometry,
  not an absolute perfused fraction, and is not comparable across tissues.
* **Reduced biology.** No adaptive resistance, no immune dynamics, no efflux or
  carrier transport (ABCB1/ABCG2 are on neither Xenium panel either), and systemic
  PK is a simplified sustained plasma level.
* **One section per tumour, no replicates.** Where a claim rests on a single
  section, the write-ups say so.
* **Research prototype.** Not a clinical decision tool. No diagnostic, regulatory
  or treatment-selection use.

The design principle behind all of this: **the modules that could flatter the
result and the modules that check it are separate, and the checking modules are
allowed to fail loudly.** Every retraction in this repo — the ENDOTHELIAL argmax
bug, the penetration-index ranking, the cross-tissue spread, the α unit error, the
hypoxia attribution — was found that way, not by the tests. Fifty-five tests were
green throughout.

---

## Where to look next

| Question | File |
|---|---|
| **What is the current state, honestly?** | [`results/SESSION_AUDIT.md`](results/SESSION_AUDIT.md) |
| I'm new here — how do I install this and get a first result? | [`docs/QUICKSTART.md`](docs/QUICKSTART.md) |
| I want to script drug/dose/sweep experiments, MATLAB-style | [`PROTOCOL_LANGUAGE.md`](PROTOCOL_LANGUAGE.md) (spec) + Studio's "protocol" tab (editor) |
| I want my team to reach Studio from somewhere other than my own machine | [`docs/DEPLOYMENT.md`](docs/DEPLOYMENT.md) |
| How do I reproduce all of this? | [`docs/REPRODUCE.md`](docs/REPRODUCE.md) |
| Where did every byte of data come from, and how do I re-download it? | [`docs/DATA_INVENTORY.md`](docs/DATA_INVENTORY.md) |
| What is the model, mathematically? | [`docs/MATHEMATICAL_MODEL.md`](docs/MATHEMATICAL_MODEL.md) |
| How do the pieces fit together? | [`docs/architecture.md`](docs/architecture.md) |
| What open data exists to validate this, and what does not? | [`results/DATA_HUNT.md`](results/DATA_HUNT.md) |
| How should I compare drugs? | [`results/delivery_metrics.md`](results/delivery_metrics.md) |
| Does the drug ranking survive a corrected delivery law? | [`results/transvascular_exchange_revalidation.md`](results/transvascular_exchange_revalidation.md) — no, not the KILL ranking |
| Is the solver trustworthy? | [`results/numerical_rigor.md`](results/numerical_rigor.md) |
| Is the vasculature real? | [`results/vascular_validation.md`](results/vascular_validation.md) |
| What did the synthetic era show? | [`archive/synthetic_era/README.md`](archive/synthetic_era/README.md) |
