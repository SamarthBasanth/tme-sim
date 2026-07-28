# Zenodo 10.5281/zenodo.8226867 (Jonas Lab / Kibur Medical glioma IMD) — deposit check

**Verdict: NOT MACHINE-READABLE.** Neither file in the Zenodo deposit, nor the paper, nor its
preprint contains a machine-readable table of per-reservoir/per-site drug concentration paired with
per-site response (caspase-3, pH2AX). The two deposited files are real data, but neither one is the
pairing this project needs: one is bulk pathway-enrichment output (not spatial, not per-site), the
other is spatially-resolved but is endogenous-metabolite MALDI-MSI from a single tissue section of a
single patient (no drug channel, no IHC response channel, no reservoir-position field). The
concentration-vs-distance and response-vs-distance data that *does* exist (Figs. 2–5 of the paper) is
presented only as plotted points on figures, in every version of the manuscript checked, including a
CC-BY-4.0 preprint of the same study. The paper's own data-availability statement is unambiguous: what
isn't in the paper/Supplementary Materials is "available on request... not publicly available due to
patient privacy."

This closes out the one open item flagged in `results/landscape_response_platforms.md` §4 point 3 — the
Jonas Lab/Kibur GBM IMD is downgraded from "glioma-specific near-miss, pending verification" to
"glioma-specific gated asset," on par with the Agar/White MALDI-MSI dataset already in
`results/DATA_HUNT.md`.

---

## 1. Deposit identity and version check

- Record fetched: **Zenodo 8226867** (https://zenodo.org/records/8226867, API:
  https://zenodo.org/api/records/8226867). Also resolves from the DOI https://doi.org/10.5281/zenodo.8226867.
- Concept DOI: 10.5281/zenodo.8226866. A search over `conceptrecid:8226866` returns exactly **one**
  version — record 8226867 is both the first and the latest (`"is_last": true`). No version history to
  worry about.
- Title: "Supplement data for: Intratumoral drug-releasing microdevices allow in situ high throughput
  pharmaco phenotyping in patients with gliomas"
- Creators: Oliver Jonas, Pierpaolo Peruzzi (Brigham and Women's Hospital)
- License: CC-BY-4.0. Access: open, no registration. Published 2023-09-04. 87 downloads at time of check.
- A second search (`q=creators.name:"Jonas" AND (microdevice OR glioma OR IMD)`) turned up no other Jonas
  Lab / Kibur Medical Zenodo deposit — this is the only one.

## 2. Full file listing

| File | Size | Format | MD5 |
|---|---|---|---|
| `Transcriptomic data - Nanostring.xlsx` | 110,314 bytes (110 KB) | Excel 2007+ (xlsx), 4 sheets | d516c972e752a3f442f8298fed8b5d5c |
| `MALDI data - glutathione and PEG.h5ad` | 25,056,046 bytes (25 MB) | AnnData / HDF5 (h5ad) | 16cfea7877fb9af188b2ea564672285d |

Direct download URLs (fetched and inspected both, in full, below):
- https://zenodo.org/records/8226867/files/Transcriptomic%20data%20-%20Nanostring.xlsx?download=1
- https://zenodo.org/records/8226867/files/MALDI%20data%20-%20glutathione%20and%20PEG.h5ad?download=1

## 3. What is actually inside each file

### 3.1 `Transcriptomic data - Nanostring.xlsx`

Opened directly (parsed the xlsx's internal XML — `xl/workbook.xml` + `xl/worksheets/*.xml` +
`xl/sharedStrings.xml` — since no `openpyxl`/`pandas` were available in the sandbox; this is a faithful
read of the raw spreadsheet, not a guess). **Four sheets, one per drug**: `Abemaciclib`, `Osimertinib`,
`Everolimus`, `Lapatinib` (575, 359, 1150, and 408 rows respectively, each 3 columns, dimension
`A1:C575` etc.).

Each sheet is a **Gene Set Enrichment Analysis (GSEA) pathway table** — Reactome/MSigDB-style pathway
name, normalized enrichment score, and p-value, e.g. (Abemaciclib sheet, rows 1–5):

```
[None, 'Normalized enrichment score', 'pValue']
['Cell Cycle', '-2.6988261869999999', '1.4549699999999999E-4']
['Cell Cycle Checkpoints', '-2.6345744209999999', '1.5671500000000001E-4']
['G2/M Checkpoints', '-2.4485839110000001', '3.24044E-4']
['Homologous DNA Pairing and Strand Exchange', '-2.3451819949999999', '5.2119500000000001E-4']
```

This is **bulk, per-drug pathway enrichment**, almost certainly pooled across whatever reservoirs/patients
received that drug (the sheet names have no patient or reservoir identifier, no distance, no dose, no
coordinate). There is no per-site row anywhere in this file. It is machine-readable, but it answers a
different question ("which pathways move with drug X") than the one this project needs ("what happened
at position Z given concentration Y at position Z").

### 3.2 `MALDI data - glutathione and PEG.h5ad`

Opened directly with `h5py` (h5ad is HDF5; `anndata`/`pandas` were not available, so this is a raw
group/dataset walk of the AnnData structure — X/obs/var — which is sufficient to read shapes, dtypes,
coordinate ranges and the full var annotation list).

Structure: `X` (sparse CSR, 8,466 pixels × 1,171 m/z features), `obs` (8,466 rows — one per pixel),
`var` (1,171 rows — one per detected m/z feature).

**`obs` fields are genuinely spatial**: `imgx`/`imgy` (pixel-grid coordinates, range 0–107 × 0–97),
`realx`/`realy` (µm-scale coordinates, range 359–466 × −230 to −133), `x`/`y` (integer versions of the
same), `run` and `section` (categorical — but **only one category each**: run =
`20200708_OllieBrain02_L5-L1_NEDC_50um`, section = `L4`). **This h5ad is a single MALDI-MSI raster from a
single tissue section of a single patient/run — not a cohort-wide, multi-patient, multi-reservoir table.**
50 µm pixel size is legible from the run name.

**`var` fields**: an m/z value plus a `metabolite` annotation string. Of 1,171 features, only **73 are
annotated to a named compound** (the rest are `unknown-N`). The 73 named compounds are entirely
**endogenous small molecules**: amino acids (alanine, serine, cysteine, glutamine, glutamate…), TCA-cycle
intermediates (succinate, malate, citrate/isocitrate, α-ketoglutarate…), nucleotides (AMP/ADP/ATP,
GMP/GDP, UMP/UDP/UTP, CMP…), lipids (C16:0, C16:1, C18:0–2, C20:4, several unlabeled `LipidN`), and
**glutathione** (var index 675) — matching the file's stated scope. **No PEG entry was found among the 73
named metabolites** despite the filename (PEG is the IMD's drug-carrier polymer; it may be present only
among the 1,098 `unknown-N` unannotated m/z peaks, or the annotation the paper used for "PEG" simply
didn't survive into this deposit's `metabolite` field — either way it isn't retrievable by name from this
file). **None of the four therapeutic drugs (temozolomide, lomustine, irinotecan, carboplatin, lapatinib,
osimertinib, abemaciclib, everolimus, doxorubicin) appear anywhere in the 73 named entries.** This file
does not contain a drug-concentration channel at all — it is an endogenous-metabolism/oxidative-stress
readout (glutathione is a redox marker, plausible as a response-adjacent proxy), not the dose axis, and
it carries no caspase-3 or pH2AX values (those are IHC/immunofluorescence readouts, a completely
different modality never represented in this h5ad).

**Bottom line on this file**: it is real, spatially-registered, per-pixel, machine-readable data — but
(a) it covers only one section from one of six patients, (b) it has no drug-concentration variable, and
(c) it has no IHC response variable. Even if reservoir coordinates were known and overlaid on this
pixel grid, the file alone cannot answer "concentration at Z vs. response at Z" because neither
concentration nor caspase-3/pH2AX is a column in it.

## 4. The actual concentration-vs-distance-vs-response data: figures only, confirmed two ways

The paper's central pairing — drug autofluorescence/MALDI-TOF concentration by distance from each
reservoir (Fig. 2), and pH2AX/cleaved-caspase-3 by distance and by patient (Figs. 3–5) — was checked two
ways:

**(a) Published version** (PMID 37672566 / PMC10754230 / DOI 10.1126/scitranslmed.adi0069). The article
is **not Open Access** (Europe PMC: `"isOpenAccess":"N"`; PMC's own OA service returns
`idIsNotOpenAccess` for PMC10754230). The single supplementary file is a combined 27.4 MB docx
(`NIHMS1945207-supplement-Supplementary_information.docx`) gated behind PMC's proof-of-work JS challenge
and not retrievable by direct fetch in this pass; ScienceDirect/science.org's own supplement endpoint
returned 403. The paper's Data Availability statement, quoted verbatim from PMC:

> "All data associated with this study are present in the paper or the Supplementary Materials.
> Additional data files related to transcriptomics and mass spectrometry that are not presented in the
> manuscript or in the supplementary materials section are available at 10.5281/zenodo.8226867."

i.e. the Zenodo deposit is explicitly scoped to be an *overflow* bucket for transcriptomics/MS files that
didn't fit in the manuscript — not a replacement for the concentration/response figures, which by this
statement's own logic must live "in the paper or the Supplementary Materials" (figures, not necessarily
tables).

**(b) A CC-BY-4.0 preprint of the same trial** (Research Square, DOI 10.21203/rs.3.rs-1647435/v2, posted
2022-05-26, same NCT04135807, same 6 patients, same authors, hosted at
https://kiburmed.com/wp-content/uploads/2022/07/Peruzzi-et-al-submitted-to-Nat-Med.pdf — fetched and read
in full, all 34 pages). This preprint is openly licensed and was readable end-to-end, unlike the paywalled
STM version, and it removes any doubt about what the underlying data actually is:

- Fig. 2a/c ("Drug release profiles from each patient for Doxorubicin and Lapatinib") is a **line plot**
  of local drug dose (µM) vs. distance (µm), one colored line per patient — a GraphPad Prism chart, not a
  table.
- Fig. 3a/b and Fig. 4b (pH2AX/CC3 vs. patient, and vs. distance) are likewise **scatter/bar plots with
  error bars**, explicitly stated to be produced in "GraphPad Prism version 9.0."
- Fig. 5 (per-drug pH2AX index per patient, 6 drugs) is again bar/scatter plots only.
- The preprint's own **Data Availability** section, verbatim:

  > "Any data which are not immediately available in the manuscript are available on request from the
  > corresponding author (PP.) Those data are not publicly available due to patient privacy."

No table anywhere in the 34-page preprint (methods, results, extended-data tables, or reference list)
contains a numeric per-reservoir concentration or per-region pH2AX/CC3 value with a distance or
coordinate identifier. Extended Data Table 1 (patient clinical/molecular characteristics) and Table 1
(trial retrieval/feasibility statistics) are the *only* two tables in the document, and neither is the
dose-response table — both are metadata about the trial, not per-site measurements. These are almost
certainly what became "table S1" and "table S2" in the final published version's in-text citations
(`[table S1]`, `[table S2]` — referenced in the PMC page but only reachable inside the gated docx), which
is consistent with them being trial/patient metadata rather than the concentration×response dataset.

## 5. Does anything let you match a "site" to a spatial location?

- The Nanostring xlsx: no. No patient ID, no reservoir ID, no distance, no coordinate — just pathway name
  × drug (aggregated).
- The MALDI h5ad: yes for pixel position within its one section (`imgx`/`imgy`/`realx`/`realy`), but no
  reservoir-ID or distance-from-reservoir field, no drug channel, and no response (IHC) channel — so even
  though it has coordinates, there is nothing to pair them with inside this file.
- The paper (any version): reservoir/distance identifiers exist implicitly as the x-axis of Fig. 2/3/4/5,
  but the underlying numbers were never tabulated in any version I could access — they exist only as
  pixel positions on a rendered chart, which would require manual digitization (e.g. WebPlotDigitizer)
  to recover approximate values, and even then only mean/SD per patient-drug, not the full per-pixel
  scatter.

## 6. Second/related deposit or dataset check

- Zenodo search across the concept record and across `creators.name:"Jonas"` found no other deposit from
  this group.
- Web search for a follow-up Kibur Medical / Jonas Lab public dataset (post-2023) found only: (a) a 2024
  AACR conference press release ("Kibur Medical Presents Preliminary Clinical Data on Intratumoral
  Microdosing of High Grade Gliomas") — no data deposit, abstract/press-release only
  (https://www.prweb.com/releases/Kibur_Medical_Presents_Preliminary_Clinical_Data_on_Intratumoral_Microdosing_of_High_Grade_Gliomas_at_the_American_Academy_of_Cancer_Research_in_Orlando_FL/prweb19277027.htm);
  (b) the Jonas Lab publications page (https://jonaslab.bwh.harvard.edu/publications/), which lists no
  glioma follow-up beyond the 2023 STM paper as of this check, though the preprint text (line ~473)
  states a **follow-up clinical trial with a two-staged procedure (72 h indwelling time) "is currently
  underway"** — worth monitoring, but nothing public yet.
- No second Zenodo/GEO/other accession for this platform in glioma was found anywhere in this pass.

## 7. What would actually be needed

To get the real per-reservoir concentration × per-region caspase-3/pH2AX pairing (not what's in this
Zenodo deposit), in order of cheapest first:

1. **Read the gated supplementary docx directly**, if institutional/journal access is available:
   `https://pmc.ncbi.nlm.nih.gov/articles/PMC10754230/bin/NIHMS1945207-supplement-Supplementary_information.docx`
   (blocked here by PMC's proof-of-work challenge + non-OA status) or via Science.org's supplement page
   for DOI 10.1126/scitranslmed.adi0069 (blocked here by 403/paywall). Given the preprint's identical
   figure set and identical data-availability language, there is no strong reason to expect this docx
   contains a numeric table that the preprint's figures don't — but it hasn't been directly opened, so
   this is the one remaining unverified path, not a closed one.
2. **Digitize the figures** (WebPlotDigitizer or similar) from Fig. 2 (concentration vs. distance, per
   patient, per drug) and Fig. 3b (pH2AX/CC3 vs. distance, per patient) in either the preprint PDF (fully
   accessible, CC-BY-4.0) or the STM version — this recovers approximate mean/SD curves, not per-pixel
   raw values, but would be sufficient for a first-pass distance-binned comparison against a transport
   model.
3. **Contact the corresponding authors directly** — the data-availability statement in both manuscript
   versions explicitly names this as the intended route: Pierpaolo Peruzzi (pperuzzi@bwh.harvard.edu) and
   Oliver Jonas (ojonas@bwh.harvard.edu), Brigham and Women's Hospital. Ask specifically for the
   per-reservoir raw values underlying Fig. 2 (concentration vs. distance, all 6 patients × 2 drugs) and
   Fig. 3b/5 (pH2AX/CC3 vs. distance/patient, all 6 patients × 8 drugs), in tabular form, noting patient
   privacy is unlikely to be a blocker since these are de-identified tissue-level measurements, not
   clinical records — the stated privacy caveat appears aimed at clinical/demographic fields (Extended
   Data Table 1), not the pharmaco-phenotyping arrays.

## 8. Sources

- Zenodo record: https://zenodo.org/records/8226867 · API: https://zenodo.org/api/records/8226867 ·
  version search: https://zenodo.org/api/records/?q=conceptrecid:8226866&all_versions=true
- Files: https://zenodo.org/records/8226867/files/Transcriptomic%20data%20-%20Nanostring.xlsx?download=1 ·
  https://zenodo.org/records/8226867/files/MALDI%20data%20-%20glutathione%20and%20PEG.h5ad?download=1
- Published paper: PMID 37672566, PMC10754230, DOI 10.1126/scitranslmed.adi0069 —
  https://pmc.ncbi.nlm.nih.gov/articles/PMC10754230/
- PMC OA-status check: https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id=PMC10754230 (returns
  `idIsNotOpenAccess`)
- Europe PMC record/supplementary-file check: https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10754230/supplementaryFiles
  (returns `"Article with id PMC10754230 is not open access"`), search:
  https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=PMC10754230&format=json
- CC-BY-4.0 preprint (same trial, fully read): Research Square DOI 10.21203/rs.3.rs-1647435/v2, PDF hosted
  at https://kiburmed.com/wp-content/uploads/2022/07/Peruzzi-et-al-submitted-to-Nat-Med.pdf
- Kibur Medical AACR 2024 press release:
  https://www.prweb.com/releases/Kibur_Medical_Presents_Preliminary_Clinical_Data_on_Intratumoral_Microdosing_of_High_Grade_Gliomas_at_the_American_Academy_of_Cancer_Research_in_Orlando_FL/prweb19277027.htm
- Jonas Lab publications page: https://jonaslab.bwh.harvard.edu/publications/
