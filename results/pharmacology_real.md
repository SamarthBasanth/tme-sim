# Real-outcome validation: the drug-POTENCY axis

> **Vascular-source promotion: NOT AFFECTED, and this was verified, not assumed.** This
> file was carrying the blanket "STALE ABSOLUTE NUMBERS" banner applied after
> [`results/vessel_promotion.md`](vessel_promotion.md) promoted the discrete
> segmented-vessel source (`vascular_mode = discrete`). The banner was **over-broad here**:
> `tmesim/validation/pharmacology_real.py` imports nothing from `tmesim.reconstruct` —
> no domain, no vascular field, no PDE — it reads HGCC and GDSC tables and molecular
> descriptors only. Re-running it under the discrete default reproduced this file
> **byte-for-byte** (diff empty apart from the removed banner), so the banner has been
> removed rather than refreshed. No `alpha`, no `dx` and no sensitivity mode enter this
> report.

First test of tmesim's molecular / pharmacology layer against a **real, measured** outcome: does structure-derived chemistry (and independently, measured GDSC potency) predict **real drug sensitivity in patient-derived glioblastoma cultures (HGCC, Johansson et al. 2020, *Cell Reports*)**?

> **SCOPE / HONESTY.** This validates the **potency axis only**. Both HGCC and GDSC are *dissociated cells* -- no tissue architecture, no vasculature, no immune/stromal niche. Nothing here supports (or tests) tmesim's spatial-microenvironment premise. Potency is also notoriously hard to predict from structure alone -- which is exactly why `tmesim.molecular.potency` uses *measured* GDSC IC50s rather than a QSAR.

## 1. Coverage (HGCC drug name -> verified structure)

- HGCC drugs total: **262** (262-compound repurposing-style screen x 85 patient-derived cultures)
- Mapped to an RDKit-verified SMILES (parses + MW within 3% of literature): **69** (26%)
- Unmapped: **193** -- obscure repurposing compounds, non-oncology agents, and code-name entries (e.g. `F4`, `G3`, `ZD`, `musc`, `PRM-123`); structures were **not fabricated** for these.
- Mapping favours the built `tmesim.molecular.library` first, then a curated MW-checked add-on table for recognised agents.

### Most-sensitive real HGCC drugs (lower mean response = more potent)

| rank | drug | mean response | n cultures |
|---|---|---|---|
| 1 | Dactinomycin | 0.210 | 45 |
| 2 | Idarubicin | 0.230 | 45 |
| 3 | Daunorubicin hydrochloride | 0.252 | 45 |
| 4 | Mitoxantrone dihydrochloride | 0.310 | 84 |
| 5 | Plicamycin | 0.314 | 85 |
| 6 | Ikarugamycin | 0.316 | 85 |
| 7 | Doxorubicin | 0.329 | 84 |
| 8 | ToM | 0.363 | 85 |
| 9 | Antibiotic A-23187 | 0.397 | 45 |
| 10 | Bortezomib | 0.421 | 85 |

## 2. Descriptor model vs. shuffle baseline (cross-validated)

8 physicochemical descriptors ['mw', 'logp', 'tpsa', 'hbd', 'hba', 'rotatable_bonds', 'aromatic_rings', 'formal_charge'] -> HGCC mean sensitivity, on **n = 69** mapped drugs. Out-of-fold (5-fold KFold) scores, with an empirical label-permutation null (300 shuffles):

| model | CV Spearman | CV R^2 | perm p-value | null Spearman mean / p95 |
|---|---|---|---|---|
| RandomForest | 0.421 | 0.414 | 0.003 | -0.024 / 0.218 |
| Ridge (std) | 0.391 | 0.277 | 0.007 | -0.032 / 0.214 |

A mean predictor has CV R^2 ~ 0 by construction; a negative CV R^2 means the descriptor model does **worse than predicting the average**. The permutation p-value is the honest test: does the model beat randomly shuffled labels?

## 3. Real-vs-real concordance (GDSC LN_IC50 vs HGCC sensitivity)

Two INDEPENDENT wet-lab screens, matched on drug name for the overlap (**n = 18** drugs). Both encode 'lower = more sensitive', so agreement is a POSITIVE Spearman -- **no model in the loop**.

- **Spearman rho = 0.827** (p = 2.35e-05), HGCC mean response vs GDSC median LN_IC50.

- `potency.resolve_ic50` marks these library drugs 'measured' from the GDSC glioma screen: ['temozolomide', 'methotrexate', 'cisplatin'] (the module's official path; we extend the same median-LN_IC50 aggregation to all 18 overlapping drugs for power).

| drug | HGCC mean response | GDSC median LN_IC50 |
|---|---|---|
| Dactinomycin | 0.210 | -4.167 |
| Mitoxantrone dihydrochloride | 0.310 | 0.900 |
| Bortezomib | 0.421 | -5.192 |
| MG-132 | 0.618 | -1.458 |
| Vinblastine Sulfate | 0.626 | -3.364 |
| Paclitaxel | 0.647 | -3.566 |
| Teniposide | 0.647 | 0.628 |
| Dasatinib | 0.656 | 0.484 |
| Vinorelbinetartrate | 0.657 | -2.592 |
| Rapamycin | 0.722 | -2.490 |
| Lapatinib | 0.743 | 3.868 |
| Sorafenib tosylate | 0.749 | 2.988 |
| Crizotinib | 0.758 | 2.948 |
| Gefitinib | 0.783 | 3.758 |
| Axitinib | 0.871 | 3.219 |
| Tamoxifen citrate | 0.916 | 4.526 |
| Nilotinib | 0.920 | 3.303 |
| Temozolomide | 1.053 | 6.407 |

## Verdict (honest)

- **Structure -> real potency: weak-but-real.** Cross-validated descriptor model on real HGCC sensitivity: RandomForest Spearman 0.421 / R^2 0.414 (permutation p = 0.003); Ridge 0.391 / 0.277 (p = 0.007). There IS a statistically real signal -- but read it carefully. It is almost certainly a **drug-CLASS confound, not a genuine within-class QSAR**: the sensitive end of this repurposing library is dominated by large, structurally distinctive natural-product cytotoxics (dactinomycin, the anthracyclines, mitoxantrone, taxanes, vinca alkaloids) that bulk descriptors trivially separate from small, inert non-oncology agents (acetaminophen, propranolol, the antihistamines). The dominant feature is simply molecular weight (univariate Spearman ~ -0.53: bigger = more sensitive), with polarity/H-bonding (TPSA, HBA, HBD ~ -0.4 to -0.47) following; aromatic-ring count carries essentially no signal (~+0.08). The model is recognising 'is this a big complex cytotoxic scaffold?', which correlates with potency here -- it is NOT predicting potency from structure within a class. True potency still depends on target engagement, which is why `tmesim.molecular.potency` relies on *measured* GDSC IC50s rather than a descriptor QSAR.

- **Real-vs-real (the more compelling check): the two independent screens agree positively** (Spearman 0.827 over 18 shared drugs). Two dissociated-cell datasets agreeing on potency ordering is the strongest real signal here -- and it is a check on the DATA, not on any tmesim model.

- **Bottom line.** The pharmacology layer's *measured* potency (GDSC) is the trustworthy component -- and it is strongly corroborated by a fully independent patient-derived screen. Structure-only descriptors show a real but confounded (class-level, not within-class) signal and should not be trusted as a potency predictor. None of this validates the spatial microenvironment model -- that requires patient tissue paired with that patient's measured response, which remains access-gated (see `tmesim.data.real.drug_response.SOURCES`).
