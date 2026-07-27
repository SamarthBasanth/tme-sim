# Curated & derived data — IRREPLACEABLE, kept in git

Everything else under `data_cache/` is gitignored because it can be re-downloaded
(see `docs/DATA_INVENTORY.md` for exact URLs and `curl` commands).

**These files cannot be.** They were hand-extracted from paper PDFs/supplements or
built as derived joins during the data hunt. Losing them means redoing manual work.

| File | What it is | Why irreplaceable |
|---|---|---|
| `autoradiography_microdialysis/parsed/regional_intratumoral_PK_anchors_literature.csv` | 27 PMID/DOI-anchored human intratumoral PK measurements | hand-curated from primary literature |
| `autoradiography_microdialysis/parsed/oh2026_neurooncol_TableS3_per_sample_tumor_drug_ratios.csv` | 125 samples / 42 human glioma cases, per-sample tumour:plasma ratios with region flags | extracted from a paper supplement |
| `autoradiography_microdialysis/parsed/jackson2018_*` | TMZ microdialysis raw table | extracted from supplement |
| `geo_organoid_screens/JOINED_GSE232173_TMZviability.csv` | GBM TMZ ex-vivo viability joined to expression | derived join (verified arithmetically correct; note it sources the operator-reproducibility table, not GEO) |
| `maldi_msi_drug/gist_imatinib_*_index.csv` | index over 110 MSI ion-image arrays | derived during harvest |
| `verified_0/gist_*clinical*.csv` | GIST spatial↔clinical join (55/55, zero orphans) | derived join |
| `pdx_encyclopedia/mayo_vaubel2020_TableS*` | 37 GBM PDX lines, per-arm survival + MGMT | paper supplements |
| `slice_organotypic/Mann2023_*` | organotypic slice dose-response | extracted from PDF text layer — note `TableS3_patient_tumors_UNALIGNED.csv` is flagged UNRESOLVED (blanks lost in extraction); do not use the patient rows |

Provenance for each is in `docs/DATA_INVENTORY.md`.
