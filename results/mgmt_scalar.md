# MGMT as a per-patient scalar for temozolomide sensitivity

**The uniform-1.0 temozolomide assumption is now optional, and the way out is a measured number.** `results/gbm_tmz_exvivo.md` established that MGMT expression predicts measured ex-vivo TMZ response (Spearman +0.516, p = 0.034, n = 17), that this is ~27 % of the between-model rank variance and an 81-point viability spread, and that **no** Xenium panel gene surrogates for it (0/396 pass FDR; panel base rate 14.1 % vs 12.3 % genome-wide, i.e. no enrichment) and GDSC cannot resolve it (rho = -0.068, p = 0.71). The fix is not a panel proxy. **MGMT does not need to be spatial — it is a per-patient scalar**, recoverable as Visium pseudobulk and routinely measured in the clinic.

Code: `tmesim/molecular/mgmt.py` (the mapping), `tmesim/validation/mgmt_scalar.py` (this analysis), and one optional `mgmt=` argument on `tmesim.reconstruct.real_pipeline.apply_drug_sensitivity`. **Omitting `mgmt` reproduces the previous behaviour exactly** — asserted numerically in `tests/test_mgmt_scalar.py`.

## 1. The calibration

Fitted on the **n = 17** patient-derived GBM stem-cell cultures that have BOTH a published viability number (BJC2023_MOESM2 Supplemental Table 2) and paired RNA-seq (GEO GSE232173). Response = **mean % cell viability at 100 uM TMZ** vs DMSO, 6-day ATP assay — *lower = more sensitive*. Predictor = `log2(MGMT CPM + 1)` from the raw counts.

Observed range: MGMT 0.00–21.12 CPM (`log2(CPM+1)` 0.000–4.467); viability 23.2–104.2 %.

### Primary mapping — exponential in kill space (the default)

`sensitivity` enters the engine **multiplicatively** (`cells/pkpd.py`: `p_kill = p_hill * sensitivity * R_hypoxia * R_quiescence`), so the mapping is fitted multiplicatively. With `kill = 1 - viability/100`, OLS of `ln(kill)` on `x = log2(CPM+1)` gives

```
s(x) = exp(-lambda * x),   lambda = 0.3429 per log2 unit
  95 % CI (t, 15 df) : 0.0059 – 0.6798
  bootstrap 95 % CI (10k)  : 0.0125 – 0.6448
  p = 0.0466   r2 = 0.239   n = 17
```

Anchored so **s = 1.0 at MGMT = 0 CPM**. That anchor is the whole point: the model's existing uniform-1.0 TMZ default becomes the *MGMT-null special case* of the mapping, so supplying no MGMT changes nothing and supplying a measured one can only scale kill **down**. Bounded in (0, 1], floored at 0.05 (a resistant tumour is not immortal: TMZ's N7-meG / N3-meA adducts are not MGMT substrates, and no zero-kill culture was observed in n = 40).

**Direction is the Esteller/Hegi biology**: Spearman(x, viability) = +0.516 (p = 0.0339) and Spearman(x, kill) = -0.530 (p = 0.0286) — **high MGMT -> repaired O6-meG -> resistant -> LOW sensitivity**.

Two of the 17 cultures measured >= 99 % viability (one at 104.2 %, i.e. more cells than the DMSO control). That is assay noise, not negative kill, so they are pinned at 0.01 kill for the log fit rather than dropped.

### Secondary mapping — linear in viability (sensitivity analysis only)

The report's own OLS, `viability% = 65.72 + 6.754 * x` (slope 95 % CI 1.47–12.04, p = 0.0156, r2 = 0.331), inverted the same way. Available as `model='linear'`. It is **not** the default because it predicts **zero kill at x = 5.08** and goes negative above that — and three of our four Visium sections sit at x >= 4.0, uncomfortably close. An unbounded multiplier is not usable in the engine.

### Does the functional form earn its keep? Leave-one-out

| model | LOO MAE (percentage points of viability) |
|---|---:|
| constant (what the model does today) | 16.98 |
| linear in viability | 13.24 |
| **exponential in kill (the default)** | **14.52** |

Both one-covariate models beat the constant out-of-sample. The linear one fits viability slightly better; the exponential one is the one that is **usable as a bounded multiplier**, and it costs 1.28 pp of LOO MAE to get that. That trade is stated, not hidden.

### Binary input (methylated / unmethylated)

For a clinically reported promoter status the **measured group contrast** is used rather than the regression: median kill 0.2944 (MGMT-silent, n = 9) vs 0.1155 (MGMT-expressing, n = 8), Mann-Whitney **p = 0.0152**, Cohen d = 1.30.

| status | meaning | sensitivity | 95 % CI |
|---|---|---:|---|
| `methylated` | promoter methylated -> MGMT silenced -> sensitive | 1.0000 | anchor |
| `unmethylated` | MGMT expressed -> resistant | 0.3923 | 0.097 – 0.629 (bootstrap) |

**Independent cross-check.** Averaging the *continuous* `s(x)` within the same two groups gives 0.942 (silent) and 0.309 (expressing), a ratio of 0.328 against the measured 0.392. Two different routes through the same data agree.

### Prefer continuous expression over binary methylation — and why

In this same cohort (MOESM3 `CoxPH_tissue_005`, n = 56 parental tumours) MGMT **expression** was prognostic for overall survival (**HR 1.85, p = 6.1e-05**, q = 0.0079), externally replicated in GLIOTRAIN (n = 70, HR 1.21, p = 0.025) — while MGMT **promoter methylation**, the routine clinical biomarker, was **not significant for OS in the same cohort (p = 0.16)**. Methylation is a thresholded, noisy proxy for the quantity that actually matters, which is how much MGMT protein the tumour can make. So `mgmt_sensitivity(level)` is the primary interface and `mgmt_sensitivity(status=...)` is the fallback.

### The mapping, tabulated

| MGMT (CPM) | log2(CPM+1) | sensitivity (exponential) | 95 % CI | sensitivity (linear) |
|---:|---:|---:|---|---:|
| 0.0 | 0.000 | 1.0000 | 1.000 – 1.000 | 1.0000 |
| 1.0 | 1.000 | 0.7097 | 0.507 – 0.994 | 0.8030 |
| 3.0 | 2.000 | 0.5037 | 0.257 – 0.988 | 0.6060 |
| 7.0 | 3.000 | 0.3575 | 0.130 – 0.982 | 0.4090 |
| 15.0 | 4.000 | 0.2538 | 0.066 – 0.977 | 0.2120 |
| 21.0 | 4.459 | 0.2168 | 0.050 – 0.974 | 0.1215 |
| 25.0 | 4.700 | 0.1996 | 0.050 – 0.973 | 0.0740 |

## 2. MGMT recovered from our own REAL Visium sections (pseudobulk)

Visium is whole-transcriptome, so MGMT *is* in the reference. Summed over all in-tissue spots it yields one honest number per section.

| section | accession | spots | % spots MGMT+ | MGMT UMI | pseudobulk CPM | log2(CPM+1) | -> TMZ sensitivity | 95 % CI |
|---|---|---:|---:|---:|---:|---:|---:|---|
| `GBM2` | GSM5833534 (GBM2) | 2,891 | 5.3 % | 165 | 7.00 | 3.000 | **0.3575** | 0.130 – 0.982 |
| `GBM1` | GSM5833533 (GBM1) | 4,764 | 10.6 % | 562 | 15.33 | 4.029 | **0.2512** | 0.065 – 0.977 |
| `GBM3` | GSM5833535 (GBM3) | 1,775 | 29.4 % | 665 | 16.61 | 4.138 | **0.2420** | 0.060 – 0.976 |
| `10x_parent_gbm` | Parent_Visium_Human_Glioblastoma | 3,468 | 25.8 % | 1119 | 25.52 | 4.729 | **0.1976** | 0.050 – 0.973 |

**Plausible per-patient MGMT range from our own data: sensitivity 0.198 – 0.357** — i.e. these four real GBM sections imply TMZ sensitivity between 20 % and 36 % of the MGMT-silent case. None of them is MGMT-silent.

**This is pseudobulk, NOT a resolvable spatial field — stated plainly.** MGMT is detected in only 5.3–29.4 % of spots and each positive spot carries ~1–2 UMI, i.e. it sits at the detection floor. A per-spot 'MGMT map' built from that would be counting noise dressed up as biology. One scalar per section is the honest product — and a scalar is exactly what the determinant needs to be.

Note also that these sections are a different assay (3' spot-based Visium) from the bulk RNA-seq the mapping was calibrated on. CPM is nominally the same unit and the values (7.0–25.5 CPM) land inside/just above the calibration range (0–21.1 CPM), but the two are not identically calibrated. The `10x_parent_gbm` value (25.5 CPM, log2 4.73) is an **extrapolation** past the calibration ceiling (4.47) and is flagged as such by `mgmt_sensitivity(..., detail=True)['extrapolated']`.

## 3. Impact on the REAL Xenium tumour at the converged dx = 12 um

*(Not run in this pass — re-run `./.venv/bin/python -m tmesim.validation.mgmt_scalar` without `--quick`.)*

## 4. Honest caveats

1. **n = 17.** The 95 % CI on lambda is 0.0059–0.6798 — the lower bound is nearly zero (no effect) and the upper bound is ~2× the point estimate. The *direction* is solid and independently replicated; the *magnitude* is approximate. Every sensitivity this module returns carries that CI in `detail=True`.
2. **The response variable is noisy by construction.** It is the operator-*reproducibility* table (2–3 operators per culture over multiple years, %CV up to 43 %, median 8.3 %), not the primary screen — the primary n = 66 screen was published as figures only. Measurement noise in y attenuates the fitted slope, so this mapping is if anything conservative.
3. **Ex-vivo 2-D cultures are not spatial tissue.** These are dissociated GSC monolayers: no architecture, no vasculature, no diffusion gradient, no immune or neuronal compartment (the authors' own GSEA shows those programmes are lost in culture). **Nothing spatial is calibrated here**, and this dataset cannot test the transport half of the model at all. The scalar multiplies the *pharmacodynamic* term only; the delivery physics is untouched.
4. **Visium MGMT is pseudobulk, not a spatial field.** 5.3–29.4 % of spots positive at ~1–2 UMI each. It is used as one number per section and labelled as such everywhere it appears.
5. **Cross-assay.** Calibrated on bulk RNA-seq CPM of GSC cultures, applied to 3' spot-based Visium pseudobulk CPM of tissue. Same nominal unit, different chemistry and different cellular composition (Visium includes stroma, immune and normal brain; the cultures are pure tumour). One of the four sections extrapolates past the calibration ceiling.
6. **A clinically measured MGMT is the proper input.** The Visium numbers here demonstrate that a real per-patient value is recoverable from our own data; they are not a claim that Visium is the right assay. In a real deployment the input is the patient's own MGMT result — ideally quantitative expression, not the binary methylation call, for the reason in §1.
7. **Temozolomide only.** The nitrosoureas (lomustine, carmustine) and procarbazine are also O6-alkylators whose lesions MGMT reverses, but no measured dose-response was fitted for them, so `mgmt=` is a recorded no-op for those drugs rather than a silent extrapolation.
8. **Still no outcome label on the Xenium section.** It is an untreated resection. Adding a measured determinant makes the prediction *less* unconditional; it does not make it a validated response.

## References

* Ntafoulis I, Kleijn A, Ju J, *et al.* Br J Cancer 2023;129(8):1327-1338. PMID **37620410**. GEO **GSE232173**.
* Hegi ME, *et al.* NEJM 2005;352:997-1003. PMID **15758010** (MGMT methylation and TMZ benefit).
* Esteller M, *et al.* NEJM 2000;343:1350-1354. PMID **11070098** (MGMT silencing and alkylator response).
* Stupp R, *et al.* NEJM 2005;352:987-996. PMID **15758009**.
* Prior work in this repo: `results/gbm_tmz_exvivo.md` (the measurement), `results/real_drug_panel.md` (the panel this impact section re-runs), `results/numerical_rigor.md` (why dx = 12 um).
