# Reproducing the real results

Everything in `results/` is regenerated from **real, open, measured data** by one
script:

```bash
./scripts/reproduce_real.sh
```

There is no synthetic path. `tmesim/data/synthetic.py` and the fabricated-geometry
loader are deleted, and the harness **hard-fails in preflight** if either
reappears. The synthetic-validation era is retired to `archive/synthetic_era/`
(its honest findings are preserved in `results/STAGE0_FINDINGS.md` and
`results/DRUG_SELECTION_RIGOR.md`).

Every run is deterministic (global seed 0). Re-running reproduces the same numbers
and the same figures.

**Companion documents.** [`DATA_INVENTORY.md`](DATA_INVENTORY.md) is the backup
record for the gitignored `data_cache/` — every dataset, its accession, its
licence, and the exact command to re-download it.
[`../results/SESSION_AUDIT.md`](../results/SESSION_AUDIT.md) is the authoritative
statement of what the regenerated numbers are entitled to claim.

---

## 1. Environment

Python ≥ 3.9. Verified on macOS (Apple Silicon), CPython **3.9.6**.

```bash
git clone <this-repo> tme-sim
cd tme-sim
python3 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install --upgrade pip          # PEP 660: pip >= 21.3, and macOS 14 ships 21.2.4
pip install -c constraints.txt -e .
```

Without the pip upgrade, `-e .` fails on the system Python with
`ERROR: File "setup.py" or "setup.cfg" not found` — there is no `setup.py` here,
so an editable install needs PEP 660 support. Measured on a clean clone,
2026-08-15.

That installs numpy, scipy, matplotlib, pandas, scikit-image, networkx,
scikit-learn, requests and tifffile (see `pyproject.toml`).

### 1a. `constraints.txt` — why the `-c` is not optional

**`pyproject.toml` declares bounds. `constraints.txt` declares versions.** If you
intend to *reproduce a number* rather than merely run the code, you need the
second one.

This matters more here than in most projects, because this repo's references are
**bit-exact**. `tests/test_time_resolved_pk.py` compares against a stored
`tests/data/static_source_reference.json` with no tolerance at all, and `results/`
quotes values to 10+ significant figures. A change of ~1e-14 in a solver AUC — the
size of a different vectorised reduction order in a newer numpy — is amplified
through the nonlinear kill model into a visible difference.

Measured, 2026-08: resolving the dependency set fresh gave

| | fresh resolve | reference environment (`.venv`, and what every number in `results/` was made with) |
|---|---|---|
| numpy | 2.4.6 | **2.0.2** |
| scipy | 1.17.1 | **1.13.1** |
| pandas | 3.0.5 | **2.3.3** |
| scikit-image | 0.26.0 | **0.24.0** |
| scikit-learn | 1.9.0 | **1.6.1** |
| matplotlib | 3.11.1 | **3.9.4** |

**5 of 6 bitwise-reference failures disappeared once the versions were pinned back
to the reference set** — e.g. a solver AUC of `20642672.153042123` against a stored
`20642672.15304199`. Without the `-c`, two installs a month apart are not
numerically identical.

What is in the file, and what is deliberately not:

* **Pinned exactly** — anything that can move a float: `numpy`, `scipy`,
  `scikit-image`, `scikit-learn`, `pandas`, `networkx`, `h5py`, `tifffile`,
  `imagecodecs`, `rdkit`, `pyarrow`, `openpyxl`, `statsmodels`, plus `joblib` and
  `threadpoolctl` (they set scikit-learn's chunk size and the BLAS thread count, and
  a different thread count is a different summation tree — same numpy, different
  last ULP). `matplotlib` is pinned too, but only for byte-identical figures; **no
  number in `results/` is computed by matplotlib**, so a mismatch there is never a
  reason to distrust a value.
* **Not pinned** — `pytest` (a test runner cannot change a float),
  `fastapi`/`uvicorn`/`httpx` (they serve JSON over numbers computed elsewhere, and
  release often enough that pinning them would make the file rot faster than the
  numerics it protects), and every transitive/platform-specific wheel
  (`pillow`, `fonttools`, `certifi`, …). That last exclusion is the point: a raw
  `pip freeze` is unusable on a second platform six weeks later.

A constraints file only *bounds* what pip installs — it never installs anything by
itself. So `-c constraints.txt` is safe to add to any `pip install` in this repo,
including the optional extras below, and packages listed there that you do not want
simply stay absent.

**To regenerate it** (after a deliberate, verified upgrade — not casually):

```bash
./.venv/bin/python -m pip freeze          # the reference environment
```

Hand-copy the new versions of the packages already listed in `constraints.txt` into
it. **Do not replace the file with the freeze output.** Then re-run the bit-exact
tests, and if a stored reference now differs, *that is the finding* — investigate
and record it. Do not update the reference value to match the new library.

```bash
./.venv/bin/python -m pytest tests/test_time_resolved_pk.py tests/test_physics.py
```

`Dockerfile` and `.github/workflows/docker.yml` both install with `-c
constraints.txt`, so the image and the informational cross-platform reference job
are pinned the same way. That is what makes a red in that job mean *genuine
linux/amd64 float divergence* rather than version drift.

### Extras that `pyproject.toml` does not declare

Three packages are imported by the codebase but are **missing from
`pyproject.toml`**. Install them or the corresponding modules degrade or fail:

```bash
pip install -c constraints.txt -e '.[validation]'
```

| package | imported by | what happens without it |
|---|---|---|
| `rdkit` | `molecular.descriptors`, `validation.moncayo_explant` | **now a CORE dependency, no longer optional** — declared in `pyproject.toml` since ribociclib and palbociclib were added, because neither has a `_FALLBACK_TABLE` entry and both compile to `mw=0.0` without it |
| `openpyxl` | `validation.gbm_tmz_exvivo` | cannot read the `.xlsx` supplements |
| `statsmodels` | `validation.moncayo_explant` | the mixed-effects fit cannot run |

Reading some `data_cache/hunt/` formats additionally needs `xlrd` (legacy `.xls`),
`pyxlsb` (the Istituto Mario Negri `.xlsb` MSI workbooks) and `rdata` (the CTRPv2
`.rds` PharmacoSet). Install them only if you are working with Tier B data:

```bash
pip install -c constraints.txt xlrd pyxlsb rdata
```

The harness preflight reports which extras are present, so you find out in
10 seconds rather than 40 minutes in.

### Verify the install

```bash
python -m tmesim.cli --help
python -m tmesim.cli molecules                                  # ~0.5 s, no data needed
python -m tmesim.cli sections                                   # 0 of 7 runnable on a clone
python -m tmesim.validation.delivery_metrics --selftest         # ~3 s, no data needed
```

`molecules` lists the 10 library GBM drugs. `--selftest` runs known-answer checks
on the metric definitions and prints the arithmetic proof that `penetration_index`
is scale-invariant over a 10⁶× dose change while `therapeutic_coverage` moves from
0.0004 to 1.0000.

The harness defaults to `./.venv/bin/python`. Override with `PY`:

```bash
PY=$(which python) ./scripts/reproduce_real.sh
```

---

## 2. Data acquisition

`data_cache/` is gitignored and holds **6.1 GB**. A fresh clone always starts from
acquisition. Three classes of input:

| class | how | size | needed for |
|---|---|---:|---|
| **AUTO** | `python -m tmesim.cli acquire` | ~90 MB | GDSC2 potency, Ivy GAP niches, TCGA background |
| **AUTO (sections)** | `acquire --sections xenium_319` (default) / `--sections all` | 16.6 MB / 304 MB | the tissue itself — reconstruct, simulate, Studio |
| **MANUAL** | `curl` recipes below | ~830 MB | Everything else in the core pipeline |
| **HUNT** | `data_cache/hunt/` | 5.3 GB | The external-validation stage only |

**The distinction that decides whether anything runs at all:** the four tabular
AUTO sources parameterise *drugs*. A tissue *section* is what a drug is
simulated on, and it lives in its own cache directory. Fetching all four tabular
sources and no section leaves every section unavailable — which is exactly what
`acquire` used to do, silently. `python -m tmesim.cli sections` is the
authoritative list of which sections are runnable right now and, for each one
that is not, the literal command that fixes it.

**The complete record — every accession, licence, byte count and re-download
command, including all 38 `data_cache/hunt/` directories — is
[`DATA_INVENTORY.md`](DATA_INVENTORY.md).** This section is the short version.

### 2a. Automatic — `python -m tmesim.cli acquire`

```bash
python -m tmesim.cli acquire                   # writes into ./data_cache/
python -m tmesim.cli acquire --list            # print the plan, download nothing
python -m tmesim.cli acquire --sections all    # + every auto-fetchable section
python -m tmesim.cli acquire --sections none   # tabular sources only
```

Best-effort live download; it never crashes, it reports per-source status, and
it ends by printing how many sections are runnable and what to run next. It
fetches:

| source | what | lands in |
|---|---|---|
| `tcga_gbm_clinical` | 80 TCGA-GBM clinical cases (GDC REST API) | `data_cache/tcga_gbm_clinical.json` |
| `tcga_gbm_expression` | 8 open STAR gene-count files | `data_cache/tcga_expression/` |
| `ivy_gap` | Ivy GAP gene-expression matrix zip (270 RNA-seq samples) | `data_cache/ivy_gap/gene_expression_matrix.zip` |
| `gdsc_glioma` | 12,001 glioma (GBM/LGG) GDSC2 dose-response rows | `data_cache/gdsc_glioma_dose_response.csv` |
| `section:xenium_319` | the 319-panel Xenium GBM section, 3 files, 16,598,845 B | `data_cache/xenium_gbm/` |

**Measured 2026-08-15 on a clean clone (macOS, CPython 3.9.6): 1 m 50 s (5/5
sources OK) and 3 m 52 s (4/5 — the Ivy GAP host returned HTTP 504, which
`acquire` reported and stepped over). Both runs left 1 of 7 sections
runnable and a first simulation ran afterwards.** The
previous claim in this file — "17 s" — did not survive re-measurement, and a
third cold-start run measured 3 m 40 s, so read this as *minutes, dominated by
the GDC endpoints*, not as a number to plan around. Re-run `acquire` freely: every file
already on disk at its expected byte count is skipped, so a repeat run is cheap
and a partial failure resumes at file granularity.

The Ivy GAP zip must then be **extracted** where the loader looks — the harness
does this for you, or by hand:

```bash
mkdir -p data_cache/ivygap
unzip -o data_cache/ivy_gap/gene_expression_matrix.zip -d data_cache/ivygap
# -> fpkm_table.csv, rows-genes.csv, columns-samples.csv, README.txt
```

Four extra Ivy GAP summary tables are not in the zip:

```bash
for n in rna_seq_samples_details sub_block_details gene_expression_details tumor_details; do
  curl -L -o "data_cache/ivygap/$n.csv" \
    "https://glioblastoma.alleninstitute.org/api/v2/gbm/$n.csv"
done
```

> **⚠** `cancerrxgene.org` returned HTTP 410 on 2026-07-27. Do not cite the GDSC
> web portal as a live URL; `acquire` uses the Sanger COG bucket
> (`https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/`), which still serves.

### 2b. Manual — the big open files

All open (no login, no DUA). Everything in this subsection can also be fetched
by hand, and the two Xenium sections no longer need to be: `acquire --sections
all` downloads them from these exact URLs with these exact byte counts. The
commands are kept here because they are the ground truth the fetcher is built
from, and because a `curl` needs neither the package installed nor the venv
active.

**Xenium GBM, 319-gene brain panel** — 40,887 cells, the primary tissue for every
headline result. ~17 MB. **Auto: `acquire --sections xenium_319` (the default).**

```bash
mkdir -p data_cache/xenium_gbm
BASE=https://cf.10xgenomics.com/samples/xenium/1.3.0/Xenium_V1_FFPE_Human_Brain_Glioblastoma_With_Addon/Xenium_V1_FFPE_Human_Brain_Glioblastoma_With_Addon
curl -L -o data_cache/xenium_gbm/Xenium_GBM_cells.csv.gz               "${BASE}_cells.csv.gz"                # 2,167,253 B
curl -L -o data_cache/xenium_gbm/Xenium_GBM_cell_feature_matrix.tar.gz "${BASE}_cell_feature_matrix.tar.gz"  # 14,279,458 B
curl -L -o data_cache/xenium_gbm/Xenium_GBM_gene_panel.json            "${BASE}_gene_panel.json"             # 152,134 B
```

**Xenium GBM, 480-gene Immuno-Oncology panel** — 816,769 cells (20× more), 11.5 ×
10.8 mm. ~288 MB. Optional: the harness warns and skips its step if absent.
**Auto: `acquire --sections xenium_io`.**

```bash
mkdir -p data_cache/xenium_gbm_io
BASE=https://cf.10xgenomics.com/samples/xenium/2.0.0/Xenium_V1_Human_Brain_GBM_FFPE/Xenium_V1_Human_Brain_GBM_FFPE
curl -L -o data_cache/xenium_gbm_io/Xenium_GBM_cells.csv.gz               "${BASE}_cells.csv.gz"                #  44,907,408 B
curl -L -o data_cache/xenium_gbm_io/Xenium_GBM_cell_feature_matrix.tar.gz "${BASE}_cell_feature_matrix.tar.gz"  # 242,459,483 B
curl -L -o data_cache/xenium_gbm_io/Xenium_GBM_gene_panel.json            "${BASE}_gene_panel.json"             #     256,404 B
```

Both sections use **identical local basenames** — the *directory* is the only
selector. Point the pipeline at either with `--sample`:

```bash
python -m tmesim.cli reconstruct --source xenium --dx 12                                    # 319-panel
python -m tmesim.cli reconstruct --source xenium --sample data_cache/xenium_gbm_io --dx 12  # 480-panel
```

**Visium GBM — 10x public "parent" section.** ~63 MB.

```bash
mkdir -p data_cache/visium_gbm && cd data_cache/visium_gbm
BASE=https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_Glioblastoma/Parent_Visium_Human_Glioblastoma
curl -L -O "${BASE}_filtered_feature_bc_matrix.tar.gz"    # 54,134,945 B
curl -L -O "${BASE}_spatial.tar.gz"                       #  9,134,896 B
tar xzf Parent_Visium_Human_Glioblastoma_filtered_feature_bc_matrix.tar.gz
tar xzf Parent_Visium_Human_Glioblastoma_spatial.tar.gz
cd ../..
```

**Visium GBM1 / GBM2 / GBM3 — GEO GSE194329** (needed for the cross-tissue stage).
~700 MB of tarballs.

```bash
cd data_cache/visium_gbm
curl -L -O https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5833nnn/GSM5833533/suppl/GSM5833533_GBM1_spaceranger_out.tar.gz  # 258,814,702 B
curl -L -O https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5833nnn/GSM5833534/suppl/GSM5833534_GBM2_spaceranger_out.tar.gz  # ~189 MB
curl -L -O https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5833nnn/GSM5833535/suppl/GSM5833535_GBM3_spaceranger_out.tar.gz  # ~252 MB
for s in GBM1 GBM2 GBM3; do mkdir -p "$s"; tar xzf GSM5833*_${s}_spaceranger_out.tar.gz -C "$s" --strip-components=1; done
cd ../..
```

**Darmanis GBM scRNA-seq atlas — GEO GSE84465** (the deconvolution reference;
without it the Visium path cannot run at all). ~21 MB.

```bash
mkdir -p data_cache/scref_gbm
curl -L -o data_cache/scref_gbm/GSE84465_GBM_All_data.csv.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84465/suppl/GSE84465_GBM_All_data.csv.gz   # 20,519,629 B
curl -L -o data_cache/scref_gbm/GSE84465_series_matrix.txt.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84465/matrix/GSE84465_series_matrix.txt.gz #    145,533 B
```

**HGCC patient-derived drug screen** (Johansson 2020, PMID 32668248). ~0.8 MB.

```bash
mkdir -p data_cache/drugresponse
curl -L -o data_cache/drugresponse/hgcc_drugresponse.json http://hgcc.se/portal/json/drugsummary.json
curl -L -o data_cache/drugresponse/hgcc_bootstrap.json    http://hgcc.se/portal/json/bootstrap.json
```

### 2c. `data_cache/hunt/` — the external-validation corpus

**5.3 GB, 38 directories, ~1,950 files**, produced by the 39-agent data hunt
(`experiments/data_hunt_workflow.js` → `results/DATA_HUNT.md`, 164 datasets
surveyed, 143 open). It is **not** reproducible with a single command and there is
no manifest for 30 of the 38 directories.

Only **four** directories are read by code today:

| directory | read by | needed to run |
|---|---|---|
| `hunt/autoradiography_microdialysis/` | `validation.glioma_regional_pk` | the `external` stage |
| `hunt/maldi_msi_drug/` + `hunt/verified_0/` | `validation.gist_imatinib` | the `external` stage |
| `hunt/geomx_dsp/` | `validation.moncayo_explant` | the `external` stage |
| `hunt/geo_organoid_screens/` | `validation.gbm_tmz_exvivo` | the `external` stage |

`hunt/maldi_msi_drug/` ships its own harvest scripts
(`harvest_gist_imatinib.py`, `harvest_gist_quantitative.py`) against the METASPACE
GraphQL endpoint — that directory is self-reproducing. The others are not.

**Two files are hand-curated and cannot be re-fetched with a `curl`:**
`hunt/autoradiography_microdialysis/parsed/*.csv` (extracted from seven papers)
and `hunt/geomx_dsp/NatCommun2026_SourceData_SupFig2_MSI_perpixel.csv.gz`
(extracted from sheet `SupF. 2 (MSI)` of a 130 MB journal Source Data workbook
that is in *neither* the Zenodo nor the GEO deposit). **Back these two up
separately** — see `DATA_INVENTORY.md` §B1, §B2 and "Known gaps".

Without the hunt caches, run everything except the external stage:

```bash
./scripts/reproduce_real.sh --from preflight   # then Ctrl-C, or:
for s in reconstruct flagship panel cross validate; do ./scripts/reproduce_real.sh --only $s; done
```

### 2d. What cannot be downloaded

Nothing in the open tier pairs **spatial GBM tissue** with that patient's
**measured drug response**. A 164-dataset sweep found not one. There is also **no
glioma + temozolomide spatial drug map anywhere** (METASPACE: zero). The
experiment that would close the loop — the Agar/White adavosertib GBM MALDI-MSI
platform, drug concentration per pixel in phase-1 patient tissue with molecular
response — has **no accession**; it is a months-long collaboration request.
Ranked access paths for the gated tier, and the confirmed dead ends, are in
`results/DATA_HUNT.md` §4.

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
executes (about 10 s) so a broken venv, a missing extra or a resurrected synthetic
module fails fast. Each step prints its **real elapsed seconds**, so the first
full run gives exact numbers for your machine.

### Stages

| # | stage | what it runs | writes |
|---|---|---|---|
| 0 | `preflight` | interpreter + imports + optional-extras report; asserts no synthetic module reappeared; delivery-metric known-answer self-test | — |
| 1 | `data` | `cli acquire` if the auto sources are missing; extracts the Ivy GAP zip; hard-fails with an exact file list if a manual dataset is absent; warns on the IO panel, the 3 GEO Visium sections and the hunt caches | `data_cache/` |
| 2 | `reconstruct` | `cli reconstruct` on Xenium 319-panel @ dx=12 µm, Xenium 480-panel IO @ dx=12 µm (skipped if absent), Visium @ native pitch; `validation.reconstruction_fidelity` | `reconstruct_*.png`, `reconstruction_fidelity.md` |
| 3 | `flagship` | `experiments/real_simulation.py`; `cli simulate` smoke check | `real_simulation.png`, `simulate_xenium_*.png` |
| 4 | `panel` | `experiments/real_drug_panel.py`; `validation.delivery_metrics` | `real_drug_panel.{md,png}`, `delivery_metrics.{md,png}` |
| 5 | `cross` | `experiments/cross_tissue.py` (5 tumours at matched dx=12 µm + the dx=100 µm control) | `cross_tissue.{md,png}` |
| 6 | `validate` | `calibration`, `pk_multidrug`, `pharmacology_real`, **`vascular_validation`**, `vessel_profile`, `numerical_rigor`, `uncertainty`, `pk_anchors` (prints only) | `calibration.md`, `pk_multidrug.md`, `pharmacology_real.md`, `vascular_validation.{md,png}`, `vessel_profile.{md,png}`, `numerical_rigor.md`, `uncertainty.{md,png}` |
| 7 | `external` | **`glioma_regional_pk`**, **`gist_imatinib`**, **`gbm_tmz_exvivo`**, **`moncayo_explant`**, `reconstruct.vessel_network` self-test | `glioma_regional_pk.{md,png}`, `gist_imatinib.{md,png}`, `gbm_tmz_exvivo.{md,png}`, `moncayo_explant.{md,png}` |
| — | `deep` (opt-in) | `niche_response`, `combination`, `dose_schedule`, `patient_report` ×2 | `niche_response.{md,png}`, `combination.{md,png}`, `dose_schedule.{md,png}`, `patient_report_{xenium,visium}.{md,png}` |

### Steps the harness deliberately does NOT run

The harness prints a `[skip]` line with the reason for each. **Do not add them
back** — their docstrings advertise `python -m …` commands that are currently
**no-ops**: the modules have no entry point and exit without writing anything.

| module | state |
|---|---|
| `validation.cellstate_validation` | `run()`/`make_figure()` exist; no `write_report`, no `__main__` |
| `validation.ranking_robustness` | **RUNNABLE — this row previously said "no `__main__`" and instructed readers to report cisplatin vs doxorubicin as a TIE. Both are superseded.** It landed, it has a `__main__`, and the measurement *contradicted* that prior: P(cisplatin = rank 1) = 0.9277 under the report's declared priors. **That in turn has been overturned** — those priors are all **common-mode**, and a **drug-specific** transport parameter breaks it: correcting the delivery law drops P(cisplatin rank 1) to **0.3965** and changes the top-1 GBM pick to vincristine (`results/transvascular_exchange_revalidation.md` §4.1). Read `results/ranking_robustness.md` — and its §8, which concludes the reportable output for the top of the panel is a **SET, not a winner** |
| `validation.vessel_network_impact` | **RUNNABLE** — has `__main__`; `results/vessel_network_impact.{md,png}` exist. (Previously listed here as having no `__main__`; that was wrong.) |
| *(hypoxia vs VEGFA/HIF1A)* | `marker_panels` points at `validation/hypoxia_validation.py`, which does not exist |
| `validation.pk_reality` | `__main__` only **prints**. `results/pk_reality.md` has **no producer** — re-running will not refresh it, and the on-disk copy is stale |

**Before trusting a docstring here, check for `if __name__ == "__main__":`.**
Several of these modules advertise a `python -m …` invocation and an output file
in their docstring while having neither. `moncayo_explant` and `gbm_tmz_exvivo`
were in this list until they gained entry points; the harness now runs both.

---

## 4. Expected runtimes

Wall-clock on an Apple Silicon laptop, single process, warm data cache.
**"measured" means stopwatch-timed in this repository.**

| step | runtime | source |
|---|---:|---|
| `cli --help` | < 1 s | measured |
| `cli molecules` | ~0.5 s | measured |
| `validation.delivery_metrics --selftest` | 3 s | measured |
| `cli acquire` | **17 s** (4/4 sources OK) | measured |
| `cli reconstruct --source xenium --dx 12` | **6 s** | measured |
| `cli reconstruct --source xenium --sample data_cache/xenium_gbm_io --dx 12` | **45 s** (816,769 cells) | measured |
| `cli reconstruct --source visium` | ~44 s | measured (earlier run) |
| `cli simulate --source xenium --dx 12 --drug temozolomide` | **3 min 58 s** | measured |
| `validation.pk_anchors` | ~5 s | measured |
| `validation.vascular_validation` | 65 s | reported in `results/vascular_validation.md` |
| `reconstruct.vessel_network` (self-test) | ~2 min | estimate — runs and passes |
| `validation.reconstruction_fidelity` | ~3 min | estimate |
| `experiments/real_simulation.py` | ~4 min | estimate |
| `experiments/real_drug_panel.py` | **19.6 min** | measured (10 drugs × 2 sensitivity modes; panel expanded from 5 to 10 drugs — see `results/real_drug_panel.md`) |
| `validation.delivery_metrics` (full) | 11 min | reported in `results/delivery_metrics.md` |
| `experiments/cross_tissue.py` | 19.5 min | reported in `results/cross_tissue.md` |
| `validation.calibration` / `pk_multidrug` / `pharmacology_real` | ~1–3 min each | estimate |
| `validation.vessel_profile` | ~6 min | estimate |
| `validation.numerical_rigor` | ~45 min | estimate (7-point dx ladder + 6-point dt ladder + 8 sensitivity re-runs) |
| `validation.uncertainty` | 24.4 min | reported in `results/uncertainty.md` |
| `validation.glioma_regional_pk` | ~5 min | estimate (20,000-permutation within-patient test) |
| `validation.gist_imatinib` | ~5 min | estimate |
| `validation.gbm_tmz_exvivo` | ~3 min | estimate |
| `validation.moncayo_explant` | ~15 min | estimate (376,736 MSI pixels + a per-slice diffusion solve) |
| `experiments/patient_report.py --source xenium` | 11 min 18 s | reported in `results/patient_report_xenium.md` |
| **core pipeline (stages 0–7)** | **~2 h 15 m** | sum of the above |
| **`--with-deep`** | **~3 h 45 m** | sum of the above |

---

## 5. Expected outputs

**Reconstruction** — `reconstruct_xenium.png`, `reconstruct_visium.png`,
`reconstruction_fidelity.md`. On the 319-panel Xenium section at dx = 12 µm expect
exactly:

> ### ⚠️ STALE BLOCK — the composition line below is PRE-CLASSIFIER-FIX, and the concordance claim under it is RETRACTED
>
> **The `TUMOR 0.7380 …` line is the *old* classifier's output and the shipped code no
> longer produces it.** The classifier-arity fix (`4ad221f`, `results/classifier_fix.md`)
> replaced the `mean_z` aggregation, which systematically over-scored classes defined by
> few markers. Two independent files agree on what the current defaults produce on this
> section — `results/revalidation.md` §4a and `results/heldout_patient.md` §1:
>
> | | pre-fix (printed below) | current defaults |
> |---|---:|---:|
> | TUMOR | 0.7380 | **0.5824** |
> | STROMAL | 0.0792 | **0.1178** |
> | IMMUNE | 0.0854 | **0.1948** |
> | ENDOTHELIAL | 0.0974 | **0.1050** |
> | mean vascular density | 0.0840 | **0.0268** (discrete source, `results/heldout_patient.md` §1) |
> | mean oxygen (tissue) | 0.5755 | **0.3520** (`results/heldout_patient.md` §1) |
>
> `n_cells`, `grid`, `extent_um` and the mask settings are unchanged. The block is left
> verbatim rather than rewritten because no run was performed for this edit; treat the
> right-hand column as the expectation and the block as a historical print.

```
n_cells 40887   grid (424, 675)   extent_um (8100.0, 5088.0)
TUMOR 0.7380  STROMAL 0.0792  IMMUNE 0.0854  ENDOTHELIAL 0.0974  NECROTIC 0.0000
mean_vascular_density 0.0840   mean_ecm_density 0.0759   mean_oxygen 0.5755
tissue_mask_mode geometric   tissue_disc_radius_um 30.0   tissue_coverage_fraction 0.6770
```

**Do not expect — and do not quote — cross-platform concordance at Pearson ≈ 0.997.**
That claim is **RETRACTED** (`results/revalidation.md` §4; `results/reconstruction_fidelity.md`
carries the retraction banner). Under current defaults the Xenium-vs-Visium 4-vector Pearson
is **0.980**, the **non-tumour subvector collapses 0.945 → 0.016** — statistically
indistinguishable from no agreement — composition rank agreement falls from an exact match
(ρ = 1.00) to **ρ = 0.40**, and mean absolute composition difference worsens **4.44% → 12.21%**.
Reverting the classifier alone reproduces the old numbers to four decimal places, which is how
the project knows the 0.997 was agreement between a *biased* estimator and another estimator
rather than a real cross-platform result. **Cross-platform composition concordance is an OPEN
question in this project, not evidence.** The endothelial column is the one that barely moved
(9.74% → 10.50% against Visium's 9.77%); the other three moved by 4–11 percentage points.

On the 480-panel IO section at dx = 12 µm expect `n_cells 816769` and
`extent_um (11472.0, 10860.0)`. **The composition figures this paragraph used to
give — TUMOR 0.7909 / STROMAL 0.0998 / IMMUNE 0.0674 / ENDOTHELIAL 0.0418 — are
also pre-classifier-fix.** `results/heldout_patient.md` §1 reports the current
defaults on this section as TUMOR **66.18%** / STROMAL **15.37%** / IMMUNE **11.39%** /
ENDOTHELIAL **7.07%**. The two panels' endothelial gates are still *not* directly
comparable — the IO panel has no PECAM1 and gates on VWF/FLT1/PLVAP instead — so
every composition difference between the sections is jointly a difference in biology
and a difference in what is measurable.

Against the one external anchor the reconstruction has (human cortex volume EM,
**8.0%** endothelial, PMID 38723085), **both** sections now land close: 10.50% on the
319-panel development section and 7.07% on the 480-panel held-out section, which
`results/heldout_patient.md` scores as **REPLICATES**. Read that as a property of
those sections and those marker gates, not as a validated property of the classifier —
the owning file says so explicitly, and the endothelial fraction only replicates after
two rounds of correcting the *comparison methodology itself*
(`results/reconstruction_validation_roadmap.md` §5).

**Flagship + panel** — `real_simulation.png`, `real_drug_panel.{md,png}`,
`delivery_metrics.{md,png}`.

**Cross-tissue** — `cross_tissue.{md,png}`, including the retraction of the earlier
dx-confounded "kill spans 2% → 24%" headline.

**Internal validation** — `calibration.md` (α ≈ 148.6, e-fold convention),
`pk_multidrug.md` (2 out-of-sample drugs inside their real cited bands),
`pharmacology_real.md` (GDSC vs HGCC Spearman +0.83),
`vascular_validation.{md,png}` (the vasculature **FAILS** the clustering test),
`vessel_profile.{md,png}`, `numerical_rigor.md`, `uncertainty.{md,png}`.

**External validation** — `glioma_regional_pk.{md,png}` (partition magnitude
passes at 1.01× bias excluding one outlier; cross-drug ordering passes at 3.65×
predicted vs 4.66× measured; **regional ordering FAILS** at exactly 1.00×),
`gist_imatinib.{md,png}` (transport-only **falsified** as a sufficient
explanation: 7.9× fold error vs a 13.0× null),
`gbm_tmz_exvivo.{md,png}` (GSE232173 — right disease, right drug, measured ATP
viability in patient-derived GSC cultures; potency axis only, no drug
concentration and no microenvironment) and `moncayo_explant.{md,png}`.

### Five corrections you should see reflected

1. **Production resolution is dx = 12 µm, and `penetration_index` is RETIRED.**
   Not converged at coarser grids, *and* even converged it is a scale-invariant
   shape statistic of the vascular field: it spans **0.1946–0.1959 (spread 0.0014)
   across the five-drug panel while their delivered exposure spans 120×**, and its
   own grid drift (0.0183) is ~13× larger than that entire signal. Use
   **`therapeutic_coverage` @ 1× IC50** (spans 0.311–0.925, separation 21.9,
   ±0.03 error bar), quoted with the 0.5×/2× band and `median_tumor_auc` (µM·s).
   **The CLI's `simulate` still prints a `penetration_index` line — ignore it.**
2. **Hypoxic fraction is grid-sensitive** and only order-of-magnitude. The same
   Visium tissue reports 100% at dx = 100 µm and 0% at dx = 12 µm. Two
   grid-dependence bugs were partially cancelling; both are fixed and the tissue
   mask is now `geometric` (dx-independent). Do not quote a hypoxic fraction as a
   measurement.
3. **α uses the e-fold convention** `L = d_half / ln 2`. Fitting the raw
   half-distance inflates α by exactly `1/ln²2 = 2.081×`; the earlier **309.6 was a
   convention error** and the correct value is **≈ 148.6**. A regenerated report
   still quoting 309.6 is stale, and absolute µM·s figures from before the fix are
   not comparable to current ones.
4. **Per-cell sensitivity is drug-specific from real markers — except for
   temozolomide.** TOP2A for TOP2 poisons (doxorubicin, etoposide, mitoxantrone);
   MKI67 as a cell-cycle proxy for cisplatin/lomustine/carmustine/methotrexate/
   vincristine; **uniform 1.0 for temozolomide, because MGMT is absent from BOTH
   Xenium panels** and its methylation is not measurable from these transcripts at
   all. The reconstruction default is uniform 1.0 for *every* drug — the only fair
   cross-drug mode, since a determinant-mode comparison quietly advantages TMZ.
5. **This correction is itself SUPERSEDED — read it as history.** It used to say
   the vascular field is a smooth 40 µm KDE, that fitted penetration length is
   26–28 µm against a source-null of 23.7 µm, that the exposure profile has a
   flat floor at ~12 % of peak, and that the discrete replacement is "built and
   self-tested but dormant". **All four are now wrong.** `_DEFAULT_VASCULAR_MODE
   = "discrete"` (`reconstruct/real_fields.py:179`) — segmented vessels are the
   shipped default, not a blur. The flat floor was re-diagnosed as diffusion
   physics rather than a KDE artefact, and the length degeneracy was traced to
   the tissue mask (`docs/VERIFICATION_METHODOLOGY.md`). What survives: **neither
   vascular mode is grid-converged at dx = 12 µm**, and `therapeutic_coverage` is
   robust because it is computed over the fixed set of real measured cells, not
   over vessel-relative geometry.

### Reproducibility notes

* Global seed 0; `--seed` on the CLI changes it. Transport and reconstruction are
  fully deterministic — only the per-cell kill draw uses the RNG, and the audits
  average it over seeds 0–3.
* `data_cache/` and `results/regime_pp*/` are git-ignored, so a fresh clone always
  starts from the acquisition step.
* Do not compare numbers across grid resolutions. Anything quoted at dx = 30 µm or
  dx = 100 µm is superseded by the dx = 12 µm value.
* The CLI `--dx` **default is 30.0**, not 12. Always pass `--dx 12` for anything
  quotable.
* The additive-leak barrier term (`reconstruct.barrier`, L ≈ 0.24) is **opt-in and
  off by default**, so a regenerated run still predicts a 1.00× regional ratio.

---

## 6. Troubleshooting

| symptom | cause / fix |
|---|---|
| `python interpreter not found at .../.venv/bin/python` | create the venv (§1) or set `PY=$(which python)` |
| `tmesim is not importable` | `pip install -c constraints.txt -e .` inside the venv |
| preflight prints `optional extras MISSING` | `pip install -c constraints.txt -e '.[validation]'` (§1) |
| `MISSING: data_cache/xenium_gbm/...` | run the manual `curl` commands in §2b |
| `WARN: data_cache/xenium_gbm_io/ absent` | optional; the IO-panel reconstruct step is skipped. Fetch it in §2b if you want it |
| `WARN: only N/3 GSE194329 Visium sections cached` | the cross-tissue stage needs GBM1/2/3; fetch them in §2b |
| `WARN: N/3 required data_cache/hunt/ caches present` | the `external` stage will fail; see `DATA_INVENTORY.md` Tier B, or skip that stage |
| Ivy GAP loader raises "extract gene_expression_matrix.zip into data_cache/ivygap/" | run the `unzip` in §2a |
| `python -m tmesim.validation.<name>` exits silently, writes nothing | the module has no `__main__` — see §3, "Steps the harness deliberately does NOT run". Currently only `cellstate_validation` (`ranking_robustness` and `vessel_network_impact` have both landed and are runnable) |
| `results/pk_reality.md` never updates | expected — no producer exists in the current `pk_reality.py`; the file is stale |
| `NotOpenSSLWarning` from urllib3 | harmless; system LibreSSL vs urllib3 v2 |
| RDKit `SMILES Parse Error: ANTIBODY:bevacizumab` | expected — bevacizumab is an antibody with no SMILES; it falls back to a tabulated descriptor row |
| a step is much slower than the table | the runtimes assume a warm `data_cache/`; the first Xenium/Visium read also extracts tarballs |
| numbers disagree with a `results/*.md` you remember | check `results/SESSION_AUDIT.md` first — nine defects were fixed and several headlines retracted |
