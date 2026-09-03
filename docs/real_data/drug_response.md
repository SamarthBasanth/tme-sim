# Real drug-response OUTCOME data — the ground-truth landscape (open vs gated)

*Last hunted live: 2026-07-26. Loader: `tmesim/data/real/drug_response.py`
(`from tmesim.data.real import load_hgcc, load_gdsc`). Cache:
`data_cache/drugresponse/`.*

This is the **outcome / ground-truth** layer: real *measured* drug response that
the simulation must eventually be scored against. The job was to find REAL,
OPEN, PAIRED (molecular/spatial + measured response) glioblastoma drug-response
data, and to be brutally honest about what is open vs access-gated.

## TL;DR (the honest headline)

- **A genuinely open, patient-derived, quantitative, molecularly-paired GBM
  drug-response dataset exists and is now cached here: the HGCC atlas
  (Johansson et al. 2020, *Cell Reports*).** 262 drugs × 85 patient-derived GBM
  cultures, each with molecular subtype + overall survival, and matched
  expression/CNA/mutation openly available. **This is the best open first
  benchmark.**
- **BUT it is a *dissociated-cell* model — there is NO tissue microenvironment
  and NO spatial context.** Same limitation as GDSC. HGCC is closer to patient
  biology than GDSC (low-passage patient-derived vs immortalised lines), but it
  is still cells in a dish.
- **The gold standard tmesim actually needs — a patient's spatial tissue paired
  with that patient's measured drug response — does not exist in any source this
  project has found, open or gated.** No open dataset closes that loop. This
  matches `docs/real_data_status.md` section C. Until 2026-09-03 this bullet named
  HTAN/dbGaP and a PDMR MTA as the gated route holding it; both were checked and
  neither does — see `docs/DATA_ACCESS_PLAN.md` §6–§7.
- GBM **organoid** papers (Jacob 2020, Golebiewska 2020) have OPEN molecular
  data but their drug readouts are **small and largely qualitative/figure-level
  (no clean per-model IC50/AUC matrix)** — useful as sanity checks, not as a
  quantitative training/validation matrix.

## Ranked source map

Ranked by (openness × how paired-with-molecular × usefulness as a quantitative
outcome benchmark). Open = direct download, no login/DUA.

| # | Source | Model type | Open? | Response readout | Paired molecular | TME? | Scale | Accession / URL |
|---|---|---|---|---|---|---|---|---|
| 1 | **HGCC — Johansson 2020, *Cell Reports*** | **patient-derived** GBM cultures | **OPEN** ✅ | normalised survival vs control (AUC-like; **lower=more sensitive**) | **OPEN**: subtype+survival (bundled); expr/CNA/mutation (portal + GEO) | ✗ dissociated | **262 drugs × 85 cultures** | PMID 32668248, doi:10.1016/j.celrep.2020.107897, http://hgcc.se/portal |
| 2 | **GDSC2 glioma (Sanger 8.5)** | immortalised lines (GBM+LGG) | **OPEN** ✅ | **LN_IC50 + AUC** (fitted) | OPEN via Cell Model Passports (generic lines) | ✗ dissociated | 286 drugs × 51 lines (34 GBM) | https://www.cancerrxgene.org |
| 3 | Jacob 2020 GBO biobank (*Cell*) | patient-derived organoids | OPEN molecular; **drug QUALITATIVE** | KI67% / size in figures; genotype→drug correlations (no numeric matrix) | OPEN: mutation panel (Table S1), exome AF (Table S3), RNA-seq GSE141947 | partial (organoid) | 4 drugs, ~8–70 GBOs | PMID 31883794, doi:10.1016/j.cell.2019.11.036, GSE141947 |
| 4 | Golebiewska 2020 PDOX (*Acta Neuropathol*) | patient-derived orthotopic xenografts + organoids | molecular OPEN; **drug figure-level**; models on request | IC50+AUC (CellTiter-Glo) targeted panels + one 41-drug array | OPEN: aCGH GSE137959, methyl GSE137845, expr GSE134470, scRNA GSE128195, seq PRJNA627814 | ✅ in-vivo stroma/vessels | ~9 drugs targeted; 41-drug on 1 model | PMID 33009951, doi:10.1007/s00401-020-02226-7 |
| 5 | HCMI GBM models (GDC) | next-gen patient-derived models | OPEN molecular; **NO drug screen** | none systematic | OPEN: WGS/WXS/RNA-seq/methylation + clinical (GDC open tier) | ✗ | dozens of GBM models | https://hcmi-searchable-catalog.nci.nih.gov/ |
| 6 | NCI PDMR | PDX / organoids | molecular open; drug partial; **models via MTA** | in-vivo/ex-vivo drug response (subset public) | OPEN: WES+RNA-seq public DB | ✅ PDX in-vivo | GBM subset | https://pdmr.cancer.gov |
| 7 | **HTAN GBM spatial multi-omics** | **patient tissue** (CODEX/MIBI/Visium+scRNA) | **GATED** 🔒 | no per-sample measured kill; SoC + outcome | GATED: dbGaP + HTAN DUA | ✅ **true spatial TME** | cohort | https://data.humantumoratlas.org |

Notes on the "partial" cells:
- **HGCC molecular is fully open** but the *large* matrices (expr/CNA/mutation)
  are left un-cached by default (URLs below) to keep the repo cache small; the
  subtype + survival needed to stratify response are bundled and cached.
- **Jacob Table S6 is a red herring for outcomes**: it holds GSEA feature lists
  + a "GBO age" sheet, **not** a numeric drug-response table. The actual drug
  effect (KI67%) is reported only in Figure 6. Verified by parsing the file.
- **Golebiewska drug data** is presented in figures/curves; there is no single
  downloadable per-model × per-drug IC50/AUC supplementary matrix, and PDOX
  models are shared "on reasonable request" (EuroPDX / PDXFinder), so it is not
  a turnkey open matrix despite fully-open omics.

## What was actually downloaded (`data_cache/drugresponse/`)

| File | Bytes | Source | Contents (verified by parsing) |
|---|---|---|---|
| `hgcc_drugresponse.json` | 257,719 | HGCC portal `json/drugsummary.json` | drug-response matrix `{drug: {culture: value}}`, **262 drugs × 85 cultures**, value = normalised survival vs control |
| `hgcc_bootstrap.json` | 507,685 | HGCC portal `json/bootstrap.json` | bundled `clinical` (101 cultures: subtype/sex/age/survival/vitalstat), `druginfo`, `genes` |
| `jacob2020_TableS1.xlsx` | 44,039 | Cell mmc1 | patient info + **CPD mutation panel** (58 GBOs × genes) |
| `jacob2020_TableS2.xlsx` | 18,373 | Cell mmc2 | biobanked GBO characterisation |
| `jacob2020_TableS3.xlsx` | 320,708 | Cell mmc3 | RNA-seq/exome QC + **per-GBO exome variant allele frequencies** |
| `jacob2020_TableS4/5.xlsx` | 36k/15k | Cell mmc4/5 | scRNA metrics; xenograft engraftment |
| `jacob2020_TableS6.xlsx` | 110,219 | Cell mmc6 | GSEA feature lists (NOT a numeric drug matrix — see note) |
| `johansson2020_mmc2.xlsx` | 15,769 | Cell Reports Table S1 | 100 HGCC cultures: data generated + subtype/sex/age/survival |
| `johansson2020_mmc3.xlsx` | 171,259 | Cell Reports Table S2 | 263-compound library metadata + per-compound **average AUC** + CV-R |
| `johansson2020_mmc4–7.xlsx` | ~10k each | Cell Reports Tables | MOA enrichment; TP53 variants; qPCR assays; validation dose ranges |

Total cached ≈ 1.5 MB. Everything above was pulled with `curl` from public,
no-login endpoints (Elsevier `ars.els-cdn.com` supplementary + the HGCC portal
JSON) and verified by opening/parsing (openpyxl / json) — see the run log.

### HGCC matched omics (open, not cached by default — for later molecular pairing)
- `http://hgcc.se/portal/json/expr.json` (~28 MB, expression)
- `http://hgcc.se/portal/json/cna.json` (~29 MB, copy number)
- `http://hgcc.se/portal/json/mutation.json` (~11 MB) / `muttable.json` (~44 MB)
- GEO mirrors: **GSE72218** (expression), **GSE72209** (CNA) — from the original
  HGCC resource paper (Xie et al. 2015, *EBioMedicine*, PMID 26629530).

## The single best open benchmark

**HGCC (Johansson et al. 2020).** Rationale:
1. **Quantitative** per-(culture, drug) sensitivity readout (not qualitative).
2. **Patient-derived** low-passage cultures spanning all transcriptional
   subtypes — biologically closer to real tumours than GDSC's immortalised lines.
3. **Molecularly paired and open**: subtype + survival bundled; full
   expression/CNA/mutation openly downloadable for the SAME cultures.
4. **Turnkey**: cached and loadable now.

Loader usage:
```python
from tmesim.data.real import load_hgcc, load_gdsc
hgcc = load_hgcc()          # 14,958 rows: culture, drug, response(+subtype,survival,...)
gdsc = load_gdsc(gbm_only=True)  # 8,384 rows: cell_line, drug, ln_ic50, auc, ...
```
`load_hgcc()` returns response as normalised survival vs control (**lower = more
sensitive**; e.g. Bortezomib mean ≈ 0.42; classic cytotoxics Dactinomycin /
anthracyclines rank most potent — a sanity check that the direction is right).

Use **GDSC (#2)** as the complementary quantitative reference (true LN_IC50/AUC,
cross-validated against HGCC on shared drug classes).

## The honest gap to gold-standard paired validation

The best open data (HGCC, GDSC) is **dissociated cells with no microenvironment
and no spatial coordinates.** tmesim's whole premise is that the *spatial
microenvironment* changes drug response. Therefore:

- HGCC/GDSC can validate the **pharmacology / Hill-kill and molecular→sensitivity**
  parts of the engine (does predicted sensitivity track measured sensitivity
  across real patient-derived models and subtypes?). This is a real, non-trivial
  check and it is achievable now.
- They **cannot** validate the microenvironment contribution — there is no
  measured "response of this spatially-resolved tissue to this drug" in any open
  source.

> **CORRECTED 2026-09-03 — the two gated routes below do not hold the data this
> section assumes they do, and neither is worth an application.** Both were checked
> in [`../DATA_ACCESS_PLAN.md`](../DATA_ACCESS_PLAN.md) §6–§7 and
> [`../outreach/03_controlled_access_skeleton.md`](../outreach/03_controlled_access_skeleton.md)
> §2: **HTAN has no GBM atlas at all**, and **NCI PDMR has no glioblastoma models and
> no public drug-response data**. Those files carry the counts and the dates; they are
> not restated here. Steps 1–3 below are kept as the record of what this document
> assumed, not as instructions to follow — do not open either application on them.

**What this section assumed would close the gold-standard loop:**
1. **HTAN GBM** spatial multi-omics (tissue microenvironment) — access step:
   register at https://data.humantumoratlas.org, obtain **dbGaP** approval +
   sign the **HTAN Data Use Agreement** for the imaging/sequence levels.
2. A **measured per-model drug response** for the same/similar patients — access
   step: **NCI PDMR** model + drug-response request + **MTA** at
   https://pdmr.cancer.gov, and/or a GBM **organoid biobank** (Jacob GBOs via the
   Penn/JHU group; Golebiewska PDOX via EuroPDX/PDXFinder "on request").
3. Validation = reconstruct the HTAN tissue, predict response with tmesim, and
   score against the paired PDMR/organoid measured kill from the matched model.

The defensible claim is: **tmesim can be benchmarked against real, open,
patient-derived pharmacology (HGCC + GDSC), and the spatial-microenvironment
contribution to response is UNTESTED — no open dataset pairs a measured spatial drug
distribution with a measured response in the same tumour.** Until 2026-09-03 that
last clause said the contribution was qualitatively validated pending gated data,
which claimed a grade of validation for a link nothing has scored, and pointed at two
routes that do not hold the data. No synthetic stand-ins were used or proposed.
