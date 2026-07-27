# Reproducing the real results

Everything in `results/` is regenerated from **real, open, measured data** by one
script:

```bash
./scripts/reproduce_real.sh
```

There is no synthetic path. `tmesim/data/synthetic.py` and the fabricated-geometry
loader are deleted; the synthetic-validation era is retired to
`archive/synthetic_era/` (its honest findings are preserved in
`results/STAGE0_FINDINGS.md` and `results/DRUG_SELECTION_RIGOR.md`).

Every run is deterministic (global seed 0). Re-running reproduces the same numbers
and the same figures.

---

## 1. Environment

Python ≥ 3.9. Tested on macOS (Apple Silicon) with CPython 3.9.6.

```bash
git clone <this-repo> tme-sim
cd tme-sim
python3 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -e .
```

That installs numpy, scipy, matplotlib, pandas, scikit-image, networkx,
scikit-learn, requests and tifffile (see `pyproject.toml`).

Verify:

```bash
python -m tmesim.cli --help
python -m tmesim.cli molecules        # ~0.4 s, no data needed
```

`molecules` lists the 10 library GBM drugs with structure-derived descriptors. If
RDKit is not installed the descriptor `source` column reads `fallback_table`
instead of `rdkit` for some rows; nothing else changes.

The harness defaults to `./.venv/bin/python`. Override with `PY`:

```bash
PY=$(which python) ./scripts/reproduce_real.sh
```

---

## 2. Data acquisition

Two classes of input, both fully open (no login, no DUA).

### 2a. Automatic — `python -m tmesim.cli acquire`

```bash
python -m tmesim.cli acquire                   # writes into ./data_cache/
python -m tmesim.cli acquire --out data_cache --timeout 60
```

Best-effort live download; it never crashes, it reports per-source status.
Approximately 2 minutes, ~90 MB. It fetches:

| source | what | lands in |
|---|---|---|
| `tcga_gbm_clinical` | 80 TCGA-GBM clinical cases (GDC REST API) | `data_cache/tcga_gbm_clinical.json` |
| `tcga_gbm_expression` | 8 open STAR gene-count files | `data_cache/tcga_expression/` |
| `ivy_gap` | Ivy GAP gene-expression matrix zip (270 RNA-seq samples) | `data_cache/ivy_gap/gene_expression_matrix.zip` |
| `gdsc_glioma` | 12,001 glioma (GBM/LGG) GDSC2 dose-response rows | `data_cache/gdsc_glioma_dose_response.csv` |

The Ivy GAP zip must then be **extracted** where the loader looks — the harness
does this for you, or do it by hand:

```bash
mkdir -p data_cache/ivygap
unzip -o data_cache/ivy_gap/gene_expression_matrix.zip -d data_cache/ivygap
# fpkm_table.csv, rows-genes.csv, columns-samples.csv, README.txt
```

Four extra Ivy GAP summary tables (per-sample structure ontology, sub-block
histology fractions, ISH energies, per-donor molecular subtype/MGMT/survival) are
not in the zip:

```bash
for n in rna_seq_samples_details sub_block_details gene_expression_details tumor_details; do
  curl -L -o "data_cache/ivygap/$n.csv" \
    "https://glioblastoma.alleninstitute.org/api/v2/gbm/$n.csv"
done
```

### 2b. Manual — the big open files

These are open but too large for the best-effort acquire step, so they are fetched
once with `curl`. All URLs below were re-verified live (HTTP 200, `content-length`
matching the byte counts recorded in `docs/real_data/`).

**10x Xenium GBM section** — 40,887 single cells, true µm coordinates, 319-gene
panel (CC BY 4.0). ~16.6 MB. This is the primary tissue for every headline result.

```bash
mkdir -p data_cache/xenium_gbm
BASE=https://cf.10xgenomics.com/samples/xenium/1.3.0/Xenium_V1_FFPE_Human_Brain_Glioblastoma_With_Addon/Xenium_V1_FFPE_Human_Brain_Glioblastoma_With_Addon
curl -L -o data_cache/xenium_gbm/Xenium_GBM_cells.csv.gz               "${BASE}_cells.csv.gz"                # 2,167,253 B
curl -L -o data_cache/xenium_gbm/Xenium_GBM_cell_feature_matrix.tar.gz "${BASE}_cell_feature_matrix.tar.gz"  # 14,279,458 B
curl -L -o data_cache/xenium_gbm/Xenium_GBM_gene_panel.json            "${BASE}_gene_panel.json"             # 152,134 B
```

**10x Visium GBM — 10x public "parent" section** (whole transcriptome). ~63 MB.

```bash
mkdir -p data_cache/visium_gbm && cd data_cache/visium_gbm
BASE=https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_Glioblastoma/Parent_Visium_Human_Glioblastoma
curl -L -O "${BASE}_filtered_feature_bc_matrix.tar.gz"    # 54,134,945 B
curl -L -O "${BASE}_spatial.tar.gz"                       #  9,134,896 B
tar xzf Parent_Visium_Human_Glioblastoma_filtered_feature_bc_matrix.tar.gz
tar xzf Parent_Visium_Human_Glioblastoma_spatial.tar.gz
cd ../..
# -> data_cache/visium_gbm/{filtered_feature_bc_matrix,spatial}
```

**Visium GBM1 / GBM2 / GBM3 — GEO GSE194329** (three further primary human GBM
sections; needed for the cross-tissue stage). ~700 MB of tarballs; keep only
`spatial/` and `filtered_feature_bc_matrix/` from each.

```bash
cd data_cache/visium_gbm
curl -L -O https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5833nnn/GSM5833533/suppl/GSM5833533_GBM1_spaceranger_out.tar.gz  # 258,814,702 B
curl -L -O https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5833nnn/GSM5833534/suppl/GSM5833534_GBM2_spaceranger_out.tar.gz  # ~189 MB
curl -L -O https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5833nnn/GSM5833535/suppl/GSM5833535_GBM3_spaceranger_out.tar.gz  # ~252 MB
for s in GBM1 GBM2 GBM3; do mkdir -p "$s"; tar xzf GSM5833*_${s}_spaceranger_out.tar.gz -C "$s" --strip-components=1; done
cd ../..
# -> data_cache/visium_gbm/GBM{1,2,3}/{filtered_feature_bc_matrix,spatial}
```

**Darmanis GBM scRNA-seq atlas — GEO GSE84465** (the deconvolution reference:
23,460 genes × 3,589 annotated cells, 7 cell types). ~21 MB.

```bash
mkdir -p data_cache/scref_gbm
curl -L -o data_cache/scref_gbm/GSE84465_GBM_All_data.csv.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84465/suppl/GSE84465_GBM_All_data.csv.gz        # 20,519,629 B
curl -L -o data_cache/scref_gbm/GSE84465_series_matrix.txt.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84465/matrix/GSE84465_series_matrix.txt.gz      #    145,533 B
```

**HGCC patient-derived drug screen** (Johansson 2020, PMID 32668248 — 262 drugs ×
85 patient-derived GBM cultures). ~0.8 MB.

```bash
mkdir -p data_cache/drugresponse
curl -L -o data_cache/drugresponse/hgcc_drugresponse.json http://hgcc.se/portal/json/drugsummary.json
curl -L -o data_cache/drugresponse/hgcc_bootstrap.json    http://hgcc.se/portal/json/bootstrap.json
```

Full provenance for each dataset — sizes, licences, what was deliberately *not*
downloaded, and the access path for the gated sources — is in `docs/real_data/`:
`singlecell_spatial.md`, `visium_gbm.md`, `scref_gbm.md`, `ivygap.md`,
`drug_response.md`, `pk_anchors.md`.

### 2c. What cannot be downloaded

Nothing in the open tier pairs **spatial GBM tissue** with that patient's
**measured drug response**. Closing that loop needs HTAN GBM spatial multi-omics
(dbGaP + HTAN DUA) plus NCI PDMR / organoid-biobank screens (MTA). The pipeline is
built to ingest exactly that; see `docs/real_data/drug_response.md` for the ranked
access path.

---

## 3. Running the harness

```bash
./scripts/reproduce_real.sh                 # core pipeline
./scripts/reproduce_real.sh --with-deep     # + niche / combination / dose / patient reports
./scripts/reproduce_real.sh --list          # print the stages
./scripts/reproduce_real.sh --dry-run       # print every command, run nothing
./scripts/reproduce_real.sh --only panel    # one stage
./scripts/reproduce_real.sh --from cross    # resume from a stage onwards
```

It runs `set -euo pipefail` — the first failure stops the run. Preflight always
executes (it is a few seconds) so a broken venv or a missing cache fails fast
rather than 40 minutes in.

### Stages

| # | stage | what it runs | writes |
|---|---|---|---|
| 0 | `preflight` | interpreter + imports; asserts no synthetic module reappeared | — |
| 1 | `data` | `tmesim.cli acquire` if the auto sources are missing; extracts the Ivy GAP zip; hard-fails with an exact file list if a manual dataset is absent | `data_cache/` |
| 2 | `reconstruct` | `cli reconstruct` on Xenium @ dx=12 µm and Visium @ native pitch; `validation.reconstruction_fidelity` | `reconstruct_*.png`, `reconstruction_fidelity.md` |
| 3 | `flagship` | `experiments/real_simulation.py`; `cli simulate` smoke check | `real_simulation.png`, `simulate_xenium_*.png` |
| 4 | `panel` | `experiments/real_drug_panel.py`; `validation.delivery_metrics` | `real_drug_panel.{md,png}`, `delivery_metrics.{md,png}` |
| 5 | `cross` | `experiments/cross_tissue.py` (5 tumors at matched dx=12 µm + the dx=100 µm control) | `cross_tissue.{md,png}` |
| 6 | `validate` | `calibration`, `pk_multidrug`, `pk_reality`, `pharmacology_real`, `vessel_profile`, `numerical_rigor`, `uncertainty` | `calibration.md`, `pk_multidrug.md`, `pk_reality.md`, `pharmacology_real.md`, `vessel_profile.{md,png}`, `numerical_rigor.md`, `uncertainty.{md,png}` |
| — | `deep` (opt-in) | `niche_response`, `combination`, `dose_schedule`, `patient_report` ×2 | `niche_response.{md,png}`, `combination.{md,png}`, `dose_schedule.{md,png}`, `patient_report_{xenium,visium}.{md,png}` |

---

## 4. Expected runtimes

Wall-clock on an Apple Silicon laptop, single process, warm data cache.

| step | runtime | source of the number |
|---|---:|---|
| `cli molecules` | 0.4 s | measured |
| `cli acquire` | ~2 min | measured (4/4 sources OK) |
| `cli reconstruct --source xenium --dx 30` | 4.5 s | measured |
| `cli reconstruct --source xenium --dx 12` | 17 s | measured |
| `cli reconstruct --source visium` | 44 s | measured |
| `cli simulate --source xenium --dx 12 --drug temozolomide` | 1 min 53 s | measured |
| `validation.reconstruction_fidelity` | ~3 min | estimate |
| `experiments/real_simulation.py` | ~4 min | estimate |
| `experiments/real_drug_panel.py` | ~10 min | estimate (5 drugs × the measured single-drug cost) |
| `validation.delivery_metrics` | 11 min | reported in `results/delivery_metrics.md` |
| `experiments/cross_tissue.py` | 19.5 min | reported in `results/cross_tissue.md` |
| `validation.calibration` / `pk_multidrug` / `pk_reality` / `pharmacology_real` | ~1–3 min each | estimate |
| `validation.vessel_profile` | ~6 min | estimate |
| `validation.numerical_rigor` | ~45 min | estimate (7-point dx ladder + 6-point dt ladder + 8 sensitivity re-runs) |
| `validation.uncertainty` | 24.4 min | reported in `results/uncertainty.md` |
| `experiments/patient_report.py --source xenium` | 11 min 18 s | reported in `results/patient_report_xenium.md` |
| **core pipeline (stages 0–6)** | **~1 h 45 m** | sum of the above |
| **`--with-deep`** | **~3 h** | sum of the above |

Steps marked *estimate* have not been individually stopwatch-timed; the harness
prints the real elapsed seconds for each step as it goes, so the first full run
gives you exact numbers for your machine.

---

## 5. Expected outputs

After a full core run, `results/` holds:

**Reconstruction** — `reconstruct_xenium.png`, `reconstruct_visium.png`,
`reconstruction_fidelity.md`. Expect TUMOR 73.8% / STROMAL 7.9% / IMMUNE 8.5% /
ENDOTHELIAL 9.7% on a 8130 × 5100 µm grid, 40,887 cells, and cross-platform
concordance vs the Visium NNLS deconvolution at Pearson ≈ 0.997.

**Flagship + panel** — `real_simulation.png`, `real_drug_panel.{md,png}`,
`delivery_metrics.{md,png}`.

**Cross-tissue** — `cross_tissue.{md,png}`, including the retraction of the earlier
dx-confounded "kill spans 2% → 24%" headline.

**Validation** — `calibration.md` (α ≈ 148.6, e-fold convention), `pk_multidrug.md`
(2 out-of-sample drugs inside their real cited bands), `pk_reality.md`,
`pharmacology_real.md` (GDSC vs HGCC Spearman +0.83), `vessel_profile.{md,png}`,
`numerical_rigor.md`, `uncertainty.{md,png}`.

### Three corrections you should see reflected

1. **Production resolution is dx = 12 µm.** `penetration_index` is not converged at
   coarser grids (it is over-stated ~2× at dx = 30 µm), *and* even at dx = 12 µm it
   is a scale-invariant shape statistic of the vascular field that does not
   discriminate drugs — it collapses to 0.194–0.196 across the whole five-drug
   panel. Compare drugs by **mean tissue AUC** or **therapeutic (IC50) coverage**
   instead. See `results/delivery_metrics.md`.
2. **Hypoxic fraction is grid-sensitive** and only order-of-magnitude. The same
   Visium tissue reports 100% hypoxic at dx = 100 µm and 0% at dx = 12 µm. Do not
   quote it as a measurement.
3. **α uses the e-fold convention** `L = d_half / ln 2`. Fitting to the raw
   half-distance inflates α by exactly `1/ln²2 = 2.081×`; the earlier **309.6 was a
   convention error** and the correct value is **≈ 148.6**. A regenerated report
   still quoting 309.6 is stale.

### Reproducibility notes

* Global seed 0; `--seed` on the CLI changes it. Transport and reconstruction are
  fully deterministic — only the per-cell kill draw uses the RNG, and the
  audits average it over seeds 0–3.
* `data_cache/` and `results/regime_pp*/` are git-ignored, so a fresh clone always
  starts from the acquisition step.
* Do not compare numbers across grid resolutions. Anything quoted at dx = 30 µm or
  dx = 100 µm is superseded by the dx = 12 µm value.

---

## 6. Troubleshooting

| symptom | cause / fix |
|---|---|
| `python interpreter not found at .../.venv/bin/python` | create the venv (§1) or set `PY=$(which python)` |
| `tmesim is not importable` | `pip install -e .` inside the venv |
| `MISSING: data_cache/xenium_gbm/...` | run the manual `curl` commands in §2b |
| `only N/3 GSE194329 Visium sections cached` (warning) | the cross-tissue stage needs GBM1/2/3; fetch them in §2b |
| Ivy GAP loader raises "extract gene_expression_matrix.zip into data_cache/ivygap/" | run the `unzip` in §2a |
| `NotOpenSSLWarning` from urllib3 | harmless; system LibreSSL vs urllib3 v2 |
| RDKit SMILES parse error for `ANTIBODY:bevacizumab` | expected — bevacizumab is an antibody with no SMILES; it falls back to a tabulated descriptor row |
| a step is much slower than the table | the runtimes assume a warm `data_cache/`; the first Xenium/Visium read also extracts tarballs |
