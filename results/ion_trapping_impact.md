# Tissue-differential lysosomal ion trapping — implemented, parameterised from measurement, and it does NOT rescue GIST

**Verdict up front.** This is the one candidate the chain had not actually excluded. It was implemented with the retardation law already validated against Crank's analytical solution, and its single free parameter was taken from the literature rather than from the target. The measurement demands a normal:tumour retardation ratio of **11.73x**; the measured rat liver lysosomal volume fraction of **1.25% of tissue volume** (PMID 33148688) can produce at most **4.47x**, and only in the physically unattainable limit of a tumour cell with no lysosomes at all. That is a **2.62x shortfall**, the same kind of bound argument that excluded hepatic uptake. Over **n = 6** matched GIST pairs the geometric-mean fold error goes from **7.92x** to **3.00x** — but only at that impossible setting, and it costs the intestinal control (1.11x -> 4.87x); at the lowest lysosomal volume fraction ever MEASURED in a cancer cell it is 4.87x, and at the highest it is 16.38x — worse than doing nothing. **This is the third excluded mechanism.**

| scorecard | before | after ion trapping | verdict |
|---|---|---|---|
| GIST tumour:normal, n = 6 pairs (geomean fold error) | **7.92x** | **3.00x** | improves, but ONLY at `phi_lys(tumour) = 0`, which does not exist |
| .. at the lowest MEASURED cancer-cell value (MCF7) | 7.92x | 4.87x | still 4.9x wrong |
| .. at the highest MEASURED cancer-cell value (MDA-MB-231) | 7.92x | 16.38x | **WORSE — the sign flips** |
| .. tumour given its host organ's measured value | 7.92x | 7.89x | the term cancels exactly, as Moncayo requires |
| GIST intestinal control (the one pair the model gets right) | 1.11x | 4.87x | **BROKEN** at the same setting that improves the liver pairs |
| human glioma partition, n = 6 drugs (GM abs fold error) | **1.995x** | **1.624x** | improves — the first mechanism in this chain that helps anything |
| .. applying R regardless of assay type (physically wrong) | 1.995x | 2.144x | worse |
| .. within 3x | 4/6 | 5/6 | |
| .. adavosertib specifically | 0.235x | 1.241x | the drug the uptake run named, and the only one that moves |

Six findings, in order of how much they are worth:

1. **The required lysosomal volume fraction is outside the measured range, and that is the result.** With the tumour granted zero lysosomes — an impossibility, chosen because it is the most generous case — the hepatocyte compartment would have to be **3.86% of liver volume**. The measured value is **1.25%**, a **3.09x** overshoot. If the tumour is instead given its host organ's content, the demand rises to **18.5%** — 3.5x more than **spleen**, the most lysosome-rich tissue in the body at 5.3%.
2. **The mechanism has a hard ceiling and the target is above it.** Henderson-Hasselbalch trapping saturates at `10^(pH_cyt - pH_lys)` = 316.2, so the whole question is the volume fraction — and at the measured one the largest tumour:normal asymmetry obtainable is 4.47x against 11.73x demanded. Ruling a mechanism out by a literature ceiling without ever fitting it is the same move that excluded hepatic uptake, now applied to its replacement.
3. **The sign of the tumour:normal contrast is not established, and at the only measured tumour values it can run backwards.** The mechanism needs normal liver to be MORE lysosome-rich than GIST tissue. No matched cancer-versus-normal lysosomal volume-fraction pair exists on any basis; the two studies that do report volumes point the other way (Glunde 2003 PMID 14965446, Kundra & Dean 1990 PMID 1693341); and at MDA-MB-231's measured 3.67% of cell volume the predicted ratio moves AWAY from the measurement, 7.92x -> 16.38x.
4. **The one configuration that helps GIST breaks the intestinal control.** Gut `phi_lys` (1.22%) and liver `phi_lys` (1.25%) are barely different, so a mechanism driven by acidic-compartment content cannot move the five liver pairs without moving the intestinal pair by nearly the same factor. It does: 1.11x -> 4.87x. The uptake run set that as the falsifiable side-condition for a real mechanism, and this one fails it.
5. **On glioma it is the first thing in this chain that helps — and the help is one drug wide.** GM absolute fold error 1.995x -> 1.624x, 4/6 -> 5/6 within 3x, driven entirely by adavosertib (0.235x -> 1.241x), the exact drug the uptake run identified as the one uptake could not fix. It is not evidence: n = 6, one mover, sign-test p = 1.00, and adavosertib's pKa is a structure class rather than a measurement.
6. **A structural observation that decides the placement question.** Ion trapping raises TOTAL tissue drug and leaves unbound EXTRACELLULAR drug alone, because a lysosome is a reservoir and at equilibrium a reservoir does not shift the free concentration outside the cell. That is the same `C_total = R * C_free` identity Moncayo proved. Applied honestly, it means the term cannot touch the three microdialysis drugs on the glioma scorecard — including eribulin, the largest single error there. Applied dishonestly (to everything), the two weak bases move by exactly opposite amounts in log space and the geometric mean barely shifts. §6.

## 1. Nothing outside the two new files was touched

The efflux run added one default-off keyword to `treatment.drugs.delivery_source`. The uptake run added nothing outside its two files. **This run does the same**: no shipped code path imports `molecular.ion_trapping`, so the term is off by construction rather than by default. All three claims are verified numerically.

- **120/120 delivery-source comparisons bit-for-bit identical** (`np.array_equal`, not `allclose`) across 4 drug names x 5 penetrances x 3 plasma concentrations x 2 branches. All identical: **True**.
- **60/60 inertness checks return exactly 1.0** — every drug in every tumour tissue (none has a measured `phi_lys`), and every non-base in every tissue that does. All exactly 1.0: **True**.
- Shipped files modified: **0**.

And the checks are not vacuous:

- the factor bites where `phi_lys` IS measured (imatinib in normal liver, R > 4): **True**
- the factor is genuinely tissue-differential: **True**
- and it cancels EXACTLY when the two tissues share a `phi_lys`, which is Moncayo's proof restated as a unit test: **True**

### The formulation is the validated one, asserted rather than claimed

`molecular/ion_trapping.py` restates `pH_partition_ratio` and `retardation_factor` rather than importing them, because `molecular` must not depend on `validation`. Restating can become drifting, so the two are compared against `validation/moncayo_explant.py`'s functions over 77 (pKa, phi) combinations with `==`, not a tolerance: **77/77 identical**, all identical **True**, pH constants identical **True**. That stepper was checked against Crank's analytical solution to 0.53%.

## 2. The parameter — where `phi_lys` comes from, and where it does not

This is the crux. A tissue-differential retardation has exactly one input per tissue, and if that input is chosen to fit the target then the whole exercise is circular. So every value carries a PMID and the provenance of each is stated, including the ones that undermine it.

### Measured, fraction of TISSUE volume

| tissue | `phi_lys` | measured? | provenance |
|---|---:|---|---|
| `generic` | -- | NO | **no measurement exists** |
| `normal_liver` | 1.250% | yes | **DIRECTLY MEASURED** |
| `liver_metastasis` | -- | NO | **no measurement exists** |
| `normal_intestine` | 1.220% | yes | scaled from liver by BMP |
| `intestinal_tumour` | -- | NO | **no measurement exists** |
| `normal_brain` | 1.360% | yes | scaled from liver by BMP |
| `glioma` | -- | NO | **no measurement exists** |
| `kidney` | 1.730% | yes | scaled from liver by BMP |
| `spleen` | 5.270% | yes | scaled from liver by BMP |
| `muscle` | 0.110% | yes | scaled from liver by BMP |
| `adipose` | 0.027% | yes | scaled from liver by BMP |
| `lung` | 1.510% | yes | scaled from liver by BMP |
| `heart` | 0.320% | yes | scaled from liver by BMP |
| `skin` | 0.160% | yes | scaled from liver by BMP |

Source: **PMID 33148688** (Schmitt 2020 *Drug Metab Dispos* 49:53-61, doi:10.1124/dmd.120.000161), Tables 1-2, rat, from the lysosomal marker lipid bis(monoacylglycero)phosphate. Range 0.03% (adipose) to 5.3% (spleen).

Four caveats, all of which cut against the mechanism this table was assembled to test:

1. **Only LIVER is a measurement.** The paper's own Table 2 footnote: *"Liver was taken as reference, whose lysosomal content was determined experimentally. Other tissues were calculated relative to the liver by bis(monoacylglycero)phosphate tissue concentrations."* These are one measurement plus a lipid-ratio extrapolation, not fourteen independent numbers. Fortunately liver is the tissue this run needs.
2. **The liver value is rat, and it is derived rather than read off.** `(1 - f_EW) x fV_lys(hepatocyte) = 0.839 x 1.49% = 1.25%`. The hepatocyte cell-basis input 1.49% +/- 0.08% is **PMID 30409837** (Schmitt 2019 *DMD* 47:49-57, doi:10.1124/dmd.118.084541), LysoTracker Red spot volumetry over >350,000 endo-/lysosomes in plated rat hepatocytes from 6 rats. The independent EM cross-check — **PMID 833203** (Blouin, Bolender & Weibel 1977 *J Cell Biol* 72:441, doi:10.1083/jcb.72.2.441, n = 5, point counting) — gives hepatocyte dense bodies 0.82% +/- 0.16% of CELL volume plus autophagic vacuoles 0.23%, i.e. 1.1% of cell = 0.92% of tissue. **The headline uses the LARGER of the two**, because LysoTracker labels all acidic compartment while EM dense bodies do not, and because it is the choice generous to the mechanism. The pessimistic arm is in §7.
3. **No normal-human liver value is retrievable at all.** NO NORMAL-HUMAN LIVER LYSOSOMAL VOLUME FRACTION IS RETRIEVABLE. Every 'human' value in the PBPK literature is rat, a textbook, or a fit: Liu & Jusko 2020 (PMID 33277347, doi:10.1124/jpet.120.000385) states 'rat lysosomal fractions were employed' and its liver 9.55% is FITTED, not anatomy; Assmus 2017 (PMID 28823852, doi:10.1016/j.ejps.2017.08.014) tabulates liver 1.12% of tissue but its caption says 'in rat' and its values are Blouin 1977's; Kazmi 2013 (PMID 23378628) 0.7% traces to a review whose own work is rat; Hallifax & Houston 2007's 1% cites a textbook; the 0.5% that circulates via Trapp 2008 (doi:10.1007/s00249-008-0338-4) is a footnoted ASSUMPTION chosen 'to get a volume ratio lysosome to cytosol of 1:200 (De Duve et al. 1974)'. The only retrievable human hepatocyte number, 3.1% of CYTOPLASM (Yu 1988, PMID 3169726, doi:10.1016/s0046-8177(88)80149-7), is from postmortem livers after shock or head trauma, includes lipofuscin granules, and that paper's own conclusion is that trauma EXPANDS hepatocyte lysosomes -- so it is not a healthy-human parameter. Rohr 1976 (PMID 823698) and Roessner 1978 (PMID 654842) certainly contain normal-human values and are not obtainable online.
4. **Beware the basis.** `% of TISSUE`, `% of CELL` and `% of CYTOPLASM` differ by ~20% and ~40% respectively and are freely interchanged in the PBPK literature. Only tissue-basis values are used directly here; `ion_trapping.cell_to_tissue_basis` does every conversion and reproduces the published liver number to 5e-5.

### NOT measured: tumour tissue, on any basis

This is the honest gap and it is not a small one.

**`liver_metastasis`** — No lysosomal volume fraction has ever been measured in GIST tissue or in GIST-T1 cells, on any basis. Both GIST/lysosome papers -- Burger 2015 PMID 26108972 doi:10.1124/mol.114.097451 and Sasaki 2021 PMID 33826667 doi:10.1371/journal.pone.0249650 -- are pharmacological and colocalisation studies with no volumetry. A Europe PMC full-text search for the exact phrase 'lysosomal volume fraction' returns 6 papers corpus-wide and none concerns GIST. Beware the review PMID 36899844, which states '70-80% of IM was found to accumulate in lysosomes ... in IM-sensitive GIST-1 cells' citing Burger 2015 -- that figure is not in Burger 2015 and must not be propagated. Callers must supply an explicit assumption; see CANCER_CELL_PHI_LYS.

**`glioma`** — The one on-target study exists and its numbers are not retrievable: Forsby, Griph & Ivarsson 1983, PMID 6137128, doi:10.1111/j.1699-0463.1983.tb02766.x, stereological volume density of secondary lysosomes in a human normal glial line vs a human malignant glioma line, on a CYTOPLASM basis. Only the direction survives in the abstract ('the volume density of lysosomes was significantly higher in the post-confluent glioma cells'); the paper has no PMCID, is not open access, is not in the Wayback Machine, and has 3 citations, so the four numbers were never quoted onward. This is an interlibrary-loan problem, not a search problem.

The closest thing to a tumour parameter that exists is the only modern set of measured cancer-cell lysosomal volume fractions, on a CELL basis:

| cell line | `phi_lys`, % of cell volume | on a tissue basis (x 0.839) |
|---|---:|---:|
| MCF7 | 0.500% | 0.419% |
| T47D | 0.780% | 0.654% |
| MDA-MB-468 | 0.835% | 0.701% |
| MDA-MB-231 | 3.670% | 3.079% |

PMID 33172973 (Collins 2021 J Pharmacol Exp Ther, doi:10.1124/jpet.120.000309): lysosomal percent volume of total cellular volume, human breast cancer lines, live-cell confocal z-stacks; 0.50% +/- 0.19% (MCF7) to 3.67% +/- 1.09% (MDA-MB-231).

Three caveats travel with it and all three cut against using it: the panel is human **breast** lines only, there is **no normal comparator in the paper**, and the internal spread is 7-fold. The paper's own point is that lysosomotropic exposure inflates the parameter 1.5-3x within 24 h — which matters directly, because imatinib is itself lysosomotropic and the GIST specimens were taken from patients on imatinib.

**And the direction of the cancer-versus-normal contrast — which is the sign of the whole mechanism — is not established.** NO MATCHED CANCER-VS-NORMAL LYSOSOMAL VOLUME-FRACTION PAIR EXISTS, on any basis. The widely repeated 'cancer cells have more lysosomes' claim rests on LAMP1 westerns, EM appearance and cathepsin activity, not volumetry, and the two studies that do report volumes point the other way: Glunde 2003 (PMID 14965446, doi:10.1016/S1476-5586(03)80037-4) gives non-malignant MCF-12A 16.8 um^3/cell against metastatic MDA-MB-231 11.2 and MDA-MB-435 11.6, i.e. the NORMAL line is higher; and Kundra & Dean 1990 (PMID 1693341, doi:10.1016/0014-4827(90)90261-8) found 'number of lysosomes and their total volume were also much larger in 3T3 cells' than in SV40-transformed SV3T3, i.e. transformation REDUCED it. Perera 2015 (PMID 26168401, doi:10.1038/nature14587) gives the strongest cancer-vs-normal contrast (15.0 +/- 1.9 vs 1.0 +/- 0.9 lysosomes per cell in PDA vs normal pancreas) but that is a COUNT, not a volume fraction.

## 3. The chemistry — pKa from structure, not from a wish list

Basic centres are found from structure with the SMARTS `moncayo_explant.py` already uses (an sp3 nitrogen that is not an amide/sulfonamide/phosphoramide N, not attached to an aromatic ring, not already charged), and the pKa is read off the same curated substructure-class table. That table cannot be tuned drug by drug, which is the point of using it. Every SMILES is asserted against its PubChem InChIKey before any chemistry runs.

| drug | InChIKey | basic centres | class | pKa (structure) | pKa (literature) | pKa used | trapping ratio |
|---|---|---:|---|---:|---:|---:|---:|
| imatinib | `KTUFNOKKBVMGRW-UHFFFAOYSA-N` | 2 | tertiary_amine | 9.8 | 8.07 | 8.07 | 278.8 |
| adavosertib | `BKWJAKQVGHWELA-UHFFFAOYSA-N` | 1 | tertiary_amine | 9.8 | -- | 9.80 | 315.4 |
| eribulin | `UFNVPOGXISZXJD-JBQZKEIOSA-N` | 1 | primary_amine | 10.0 | -- | 10.00 | 315.7 |
| levetiracetam | `HPHUVLMMVZITSG-LURJTMIESA-N` | 0 | none | 1.0 | -- | 1.00 | 1.0 |
| cefazolin | `MLYYVTUWGNIJIB-BXKDBHETSA-N` | 0 | none | 1.0 | -- | 1.00 | 1.0 |
| temozolomide | `BPEGJWRSRHCHSN-UHFFFAOYSA-N` | 0 | none | 1.0 | -- | 1.00 | 1.0 |
| methotrexate | `FBOZXECLQNJBKD-ZDUSSCGKSA-N` | 0 | none | 1.0 | -- | 1.00 | 1.0 |

Three things worth reading off that table:

- **Only three of seven drugs are weak bases at all**: imatinib (piperazine), adavosertib (4-methylpiperazine) and eribulin (primary amine). Temozolomide, cefazolin, levetiracetam and methotrexate have no basic centre that passes the filter, so the term is exactly 1.0 for them in every tissue. That is a real constraint on how much this mechanism could ever explain.
- **Imatinib uses 8.07, the repo's own value.** That is `IMATINIB_PKA_BASE` from `validation/gist_imatinib.py`, whose comment records the literature range as 7.7-8.07. The structure class would give 9.8, which produces a LARGER trapping ratio — using the repo's smaller literature value is the conservative choice against the mechanism, and both are run in §7.
- **Adavosertib's pKa could not be found anywhere** — not in the primary literature, not in vendor datasheets, not in the phase I PK papers. The structure-class 9.8 is used and swept over 6.0-10.4 in §6 rather than defended, because the glioma result rests on this drug alone.

Note the ceiling this imposes. The trapping ratio saturates at `10^(pH_cyt - pH_lys)` = 316.2, so imatinib's 278.8 is already 88% of the maximum any base can reach. There is no stronger base to appeal to; the only remaining lever is the volume fraction, and §4 is about that lever failing.

## 4. Scorecard 1 — GIST tumour:normal, n = 6 matched pairs

The prediction is the same reaction-diffusion delivery ratio the published run scored — same 1-D Krogh slab, same intervessel half-spacings, same penetration length `L = 35.87` um from the doxorubicin-anchored `alpha` fitted once to an unrelated measurement — with one tissue-differential retardation on each side:

```
  ratio = (R_tumour * E_tumour) / (R_normal * E_normal)
  R     = 1 + phi_lys * (trapping_ratio - 1)
```

`E_tumour = 0.5572`, `E_liver = 0.9614`, `E_intestine = 0.8644` — all unchanged. The ONLY new quantity is `phi_lys`.

### 4.1 Five configurations, generous ones first

| configuration | `phi_lys` tumour | `R_normal` | `R_tumour` | scale `R_t/R_n` | predicted liver ratio | geomean fold error (n = 6) | within 3x |
|---|---:|---:|---:|---:|---:|---:|---:|
| `off` | -- | 1.000 | 1.000 | 1.0000 | 0.5796 | **7.92x** | 1/6 |
| `limit` | 0.000% | 4.472 | 1.000 | 0.2236 | 0.1296 | **3.00x** | 4/6 |
| `collins_low` | 0.419% | 4.472 | 2.165 | 0.4842 | 0.2806 | **4.87x** | 2/6 |
| `collins_high` | 3.079% | 4.472 | 9.553 | 2.1361 | 1.2380 | **16.38x** | 1/6 |
| `host` | 1.250% | 4.472 | 4.472 | 1.0000 | 0.5796 | **7.89x** | 1/6 |

Measured liver geomean: **0.0494**. The no-partitioning null is 12.96x, so 7.92x was already beating it and 3.00x beats it further — but only at a value of `phi_lys(tumour)` that no cell has.

- **`off`** — mechanism off -- the shipped model
- **`limit`** — measured normal phi_lys, tumour phi_lys = 0 (physically unattainable best case; an upper bound, not a prediction)
- **`collins_low`** — tumour phi_lys from the LOWEST measured cancer line, MCF7 0.50% of cell volume -> 0.42% of tissue (PMID 33172973)
- **`collins_high`** — tumour phi_lys from the HIGHEST measured cancer line, MDA-MB-231 3.67% of cell volume -> 3.08% of tissue (PMID 33172973)
- **`host`** — tumour tissue given its host organ's measured value -- the term cancels, which is Moncayo's proof reappearing

### 4.2 The whole response surface, so nothing is hidden by a choice

No tumour value is defensible, so none is asserted. The sweep runs from the impossible zero, through the measured cancer-cell range, to the host-organ value where the term cancels.

| `phi_lys` tumour | `R_tumour` | scale | predicted liver | geomean fold error | within 3x | intestinal control | nominal floor | liver pairs still below it |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.000% | 1.000 | 0.2236 | 0.1296 | 3.00x | 4/6 | 4.87x | 0.0466 | 3/5 |
| 0.100% | 1.278 | 0.2857 | 0.1656 | 3.42x | 1/6 | 3.81x | 0.0595 | 4/5 |
| 0.200% | 1.556 | 0.3478 | 0.2016 | 3.90x | 1/6 | 3.13x | 0.0725 | 4/5 |
| 0.419% **(MCF7)** | 2.165 | 0.4842 | 0.2806 | 4.87x | 2/6 | 2.25x | 0.1009 | 4/5 |
| 0.600% | 2.667 | 0.5963 | 0.3456 | 5.59x | 2/6 | 1.83x | 0.1242 | 4/5 |
| 0.800% | 3.222 | 0.7205 | 0.4176 | 6.34x | 2/6 | 1.51x | 0.1501 | 5/5 |
| 1.000% | 3.778 | 0.8447 | 0.4896 | 7.05x | 1/6 | 1.29x | 0.1760 | 5/5 |
| 1.250% **(= liver)** | 4.472 | 1.0000 | 0.5796 | 7.89x | 1/6 | 1.09x | 0.2083 | 5/5 |
| 1.500% | 5.166 | 1.1553 | 0.6695 | 8.86x | 1/6 | 1.06x | 0.2407 | 5/5 |
| 2.000% | 6.555 | 1.4658 | 0.8495 | 11.24x | 1/6 | 1.34x | 0.3054 | 5/5 |
| 2.500% | 7.944 | 1.7764 | 1.0295 | 13.62x | 1/6 | 1.63x | 0.3701 | 5/5 |
| 3.079% **(MDA-MB-231)** | 9.553 | 2.1361 | 1.2380 | 16.38x | 1/6 | 1.96x | 0.4450 | 5/5 |
| 4.000% | 12.110 | 2.7080 | 1.5695 | 20.77x | 1/6 | 2.48x | 0.5642 | 5/5 |
| 5.270% | 15.638 | 3.4969 | 2.0266 | 26.82x | 0/6 | 3.21x | 0.7285 | 5/5 |

**The surface is monotone.** The best the mechanism ever does is 3.00x at `phi_lys(tumour) = 0`, and it degrades from there. There is no interior optimum to find and nothing to tune toward.

### 4.3 The structural floor — broken, but not where it matters

`results/gist_imatinib.md` proved the transport-only model has a hard floor at any penetration length: `E_t/E_n >= a_liver/a_tumour`. With a retardation on, that becomes `>= (a_liver/a_tumour) * (R_t/R_n)`.

| configuration | nominal floor | liver pairs below it | extreme-geometry floor | liver pairs below it |
|---|---:|---:|---:|---:|
| `off` | 0.2083 | 5/5 | 0.0833 | 4/5 |
| `limit` | 0.0466 | 3/5 | 0.0186 | 0/5 |
| `collins_low` | 0.1009 | 4/5 | 0.0403 | 1/5 |
| `collins_high` | 0.4450 | 5/5 | 0.1780 | 5/5 |
| `host` | 0.2083 | 5/5 | 0.0833 | 4/5 |

So the honest answer to "is the floor broken?" is **partly, and only in the impossible case**. At `phi_lys(tumour) = 0` the extreme-geometry floor drops from 0.0833 to 0.0186 and all five liver ratios become reachable — the previously-unreachable 4/5 becomes 0/5. But the **nominal** floor only falls to 0.0466 and 3/5 liver ratios are still below it. At the lowest MEASURED cancer-cell value the nominal floor is 0.1009 with 4/5 still below, and at the highest it is above the starting floor. **Breaking a floor by assuming a tumour has no lysosomes is not breaking it.**

### 4.4 Per patient, at the most generous setting

| patient | organ | measured | predicted (before) | predicted (limit) | fold before | fold at limit |
|---|---|---:|---:|---:|---:|---:|
| A | Liver | 0.0231 | 0.5796 | 0.1296 | 25.06x | 5.60x |
| AA | Intestine | 0.7159 | 0.6446 | 0.1469 | 1.11x | 4.87x |
| J | Liver | 0.1418 | 0.5796 | 0.1296 | 4.09x | 1.09x |
| K | Liver | 0.0439 | 0.5796 | 0.1296 | 13.20x | 2.95x |
| L | Liver | 0.0433 | 0.5796 | 0.1296 | 13.39x | 2.99x |
| U | Liver | 0.0473 | 0.5796 | 0.1296 | 12.26x | 2.74x |

**Patient AA is the control and it is where the mechanism gives itself away.** It is an intestinal primary against normal intestine, the one pair the shipped model already predicts to 1.11x. Gut `phi_lys` (1.22%) is within 3% of liver's (1.25%), so a mechanism made of acidic-compartment content applies there with almost the same strength and the control moves to 4.87x. Control survives: **False**. The uptake run set this as the falsifiable side-condition — *"a term that fixed the liver pairs by breaking the intestinal pair would be a fitted contrast, not a mechanism"* — and ion trapping fails it, where a hepatocyte-restricted transporter at least would not have.

## 5. Required versus measured — the bound argument

This is the section that decides the question, and it is set up so that **every discretionary choice favours the mechanism**:

- lysosomal lumen pH **4.7** rather than the 5.3 that Schmitt 2020's own PBPK model uses — a lower lumen pH means stronger trapping;
- the LysoTracker-based liver `phi_lys` of 1.25% rather than the EM-based 0.92%;
- tumour `phi_lys` set to exactly **zero**, a tumour cell containing no acidic compartment at all.

| quantity | value |
|---|---:|
| imatinib basic pKa (the repo's own value) | 8.07 |
| lysosome:cytosol trapping ratio at pH 4.7 / 7.2 | 278.8 |
| .. and its saturating ceiling for ANY base | 316.2 |
| `R_normal / R_tumour` the measurement demands | **11.73** |
| implied hepatocyte `phi_lys` (tumour = 0) | **3.86%** |
| MEASURED rat liver `phi_lys` (PMID 33148688) | **1.25%** |
| required / measured | **3.09x** |
| **is the required value inside the measured range?** | **False** |
| largest `R_n/R_t` the measured value can produce | **4.47** |
| **shortfall** | **2.62x** |
| implied hepatocyte `phi_lys` if tumour = host organ | 18.5% |
| highest `phi_lys` of ANY tissue in the body (spleen) | 5.3% |
| tumour `phi_lys` needed if liver had SPLEEN's content | 0.12% |
| lowest MEASURED cancer-cell `phi_lys`, tissue basis | 0.42% |

**Read the last four rows together, because that is where the argument closes.** Taken alone, the required 3.86% is not off the biological scale — spleen reaches 5.3%. But it is required of **liver**, whose measured value is 1.25%, and no tissue other than spleen exceeds 1.73%. And even granting liver the spleen's content, the tumour would then have to sit at 0.12% — below skeletal muscle (0.11%) and 3.5x below the lowest cancer cell line anyone has measured. Both ends have to be wrong at once, in opposite directions, for the mechanism to work.

**Excluded by its own literature: True.**

### What this corrects in `results/uptake_impact.md` §10

That section sized the same quantity and got 3.86% — the same number this run gets — and concluded the mechanism *"cannot be dismissed on magnitude"* because 3.86% *"sits inside the range `moncayo_explant.py` already sweeps"*. **That was the wrong comparator.** The Moncayo sweep (0.5-5%) was a plausibility range chosen to bracket an unknown, not a measurement, and its own reference value of 1.4% was assembled as `0.7 cellularity x 0.02 lysosomal fraction` rather than measured. Compared against what was actually measured in liver, 3.86% is a 3.09x overshoot. The redirect was right that Moncayo's *dispersion* proof does not apply to a between-tissue ratio; it was wrong that the magnitude was reachable.

## 6. Scorecard 2 — the six-drug measured human glioma partition

The BEFORE column is `glioma_regional_pk.run_scorecard(aggregation="geometric", apply_leak=True)` called live. AFTER multiplies the whole predicted partition by `R`, and that placement differs on purpose from how efflux and uptake were placed: those are membrane transporters and act on the route a molecule takes, whereas ion trapping is an intracellular equilibrium acting on whatever free drug is present however it arrived, including drug that came through a broken tight junction. There is no transcellular-only version of this mechanism.

The old ceiling of 1.0 does not apply either, for the reason `reconstruct.barrier.effective_penetrance` gives in its own docstring: *"a partition above 1 requires tissue binding or active uptake"*. Lysosomal sequestration IS tissue binding, so the total-drug ceiling rises to `R`. `phi_lys(glioma)` = 1.36%, assumed equal to the MEASURED normal-brain value (1.36% of tissue, PMID 33148688); no glioma measurement is retrievable.

### 6.1 The placement question is physics, and it is decidable

A lysosome is an intracellular reservoir. At equilibrium it raises TOTAL tissue drug by `R` and leaves the unbound EXTRACELLULAR concentration alone — the same `C_total = R * C_free` identity Moncayo proved. So a tissue-homogenate LC-MS/MS measurement sees `R` and an intracerebral microdialysis probe does not. Three of the six drugs on this scorecard are microdialysis, including eribulin, which is both a strong base and the largest single error on the board.

| drug | measurement | sees total drug? | `R` from chemistry | `R` applied | predicted before | predicted after | fold before | fold after |
|---|---|---|---:|---:|---:|---:|---:|---:|
| eribulin | `intracerebral_microdialysis` | **no** | 5.280 | 1.000 | 0.0353 | 0.0353 | 6.944x | 6.944x |
| cefazolin | `tissue_homogenate_LCMS` | yes | 1.000 | 1.000 | 0.1526 | 0.1526 | 1.020x | 1.020x |
| methotrexate | `intracerebral_microdialysis` | **no** | 1.000 | 1.000 | 0.1415 | 0.1415 | 0.795x | 0.795x |
| temozolomide | `intracerebral_microdialysis` | **no** | 1.000 | 1.000 | 0.1594 | 0.1594 | 0.817x | 0.817x |
| levetiracetam | `tissue_homogenate_LCMS` | yes | 1.000 | 1.000 | 0.9466 | 0.9466 | 1.358x | 1.358x |
| adavosertib | `tissue_homogenate_LCMS` | yes | 5.276 | 5.276 | 0.5785 | 3.0522 | 0.235x | 1.241x |

| aggregate | before | after (method-aware) | after (applied to everything) |
|---|---:|---:|---:|
| geometric-mean ABS fold error | 1.9946x | **1.6244x** | 2.1436x |
| geometric-mean signed fold error | 1.0662x | 1.4067x | 1.8564x |
| within 3x | 4/6 | 5/6 | 5/6 |
| within 10x | 6/6 | 6/6 | 5/6 |
| worst drug | 6.94x | 6.94x | 36.67x |
| improved? | -- | **True** | False |

**A structural note on the wrong placement, because it is not obvious.** Under `all_tissues` both weak bases get the same `R` (their pKa values both saturate the trapping ratio), eribulin is over-predicted and adavosertib under-predicted, so the two `|log10 fold|` values move by `+log10(R)` and `-log10(R)` and the geometric mean is *exactly* unchanged — until adavosertib crosses 1.0 and starts moving back the other way. That is why several rows in the §7 sensitivity table read 1.995x on both sides. The apparent robustness is an artefact of two errors cancelling, not a result.

### 6.2 It is one drug wide, and that drug's pKa is not measured

2/6 drugs are weak bases and 1/6 move under the method-aware placement. Sign test over |log10 fold error|: 1 better, 0 worse, p = 1.00. with 1 of 6 drugs moving at all, no test can reach significance at any effect size; reported to make that explicit rather than to claim a result.

| drug dropped | GM abs before | GM abs after |
|---|---:|---:|
| eribulin | 1.5542x | 1.2148x |
| cefazolin | 2.2808x | 1.7828x |
| methotrexate | 2.1872x | 1.7096x |
| temozolomide | 2.1992x | 1.7190x |
| levetiracetam | 2.1541x | 1.6837x |
| adavosertib | 1.7143x | 1.7143x |

Adavosertib's pKa could not be found in any retrievable source, so it is swept rather than defended:

| assumed pKa | `R(glioma)` | adavosertib predicted | adavosertib fold | GM abs fold error (n = 6) |
|---:|---:|---:|---:|---:|
| 6.0 | 1.254 | 0.726 | 0.295x | 1.9206x |
| 6.5 | 1.713 | 0.991 | 0.403x | 1.8234x |
| 7.0 | 2.659 | 1.538 | 0.625x | 1.6946x |
| 7.5 | 3.856 | 2.230 | 0.907x | 1.5928x |
| 8.0 | 4.701 | 2.719 | 1.105x | 1.5934x |
| 8.5 | 5.082 | 2.940 | 1.195x | 1.6143x |
| 9.0 | 5.220 | 3.020 | 1.228x | 1.6215x |
| 9.8 | 5.276 | 3.052 | 1.241x | 1.6244x |
| 10.4 | 5.284 | 3.057 | 1.243x | 1.6248x |

The improvement is robust to the pKa above about 7.0 and disappears below 6.5. An aniline-attached 4-methylpiperazine is basic in that range beyond reasonable doubt, so the direction survives — but note that a pKa of about 7.7 would land adavosertib exactly on its measured value, which is a fact to be reported and not a parameter to be chosen.

## 7. Does the verdict depend on anything we chose?

Four axes, each of which could have been chosen the other way. The GIST rows are all evaluated at `phi_lys(tumour) = 0`, i.e. at the mechanism's best possible case.

| lumen pH | imatinib pKa | liver `phi_lys` | `R_liver` | GM fold error at the limit | required `phi_lys` | required / measured | shortfall | still excluded? |
|---:|---:|---|---:|---:|---:|---:|---:|---|
| 4.7 | 8.07 (literature) | 1.25% (schmitt_lysotracker) | 4.47 | 3.01x | 3.9% | 3.09x | 2.62x | **True** |
| 4.7 | 8.07 (literature) | 0.92% (blouin_em) | 3.56 | 3.39x | 3.9% | 4.18x | 3.29x | **True** |
| 4.7 | 9.80 (structure_class) | 1.25% (schmitt_lysotracker) | 4.93 | 2.91x | 3.4% | 2.73x | 2.38x | **True** |
| 4.7 | 9.80 (structure_class) | 0.92% (blouin_em) | 3.90 | 3.19x | 3.4% | 3.70x | 3.01x | **True** |
| 5.3 | 8.07 (literature) | 1.25% (schmitt_lysotracker) | 1.86 | 5.23x | 15.5% | 12.42x | 6.29x | **True** |
| 5.3 | 8.07 (literature) | 0.92% (blouin_em) | 1.64 | 5.70x | 15.5% | 16.82x | 7.16x | **True** |
| 5.3 | 9.80 (structure_class) | 1.25% (schmitt_lysotracker) | 1.98 | 5.02x | 13.7% | 10.97x | 5.93x | **True** |
| 5.3 | 9.80 (structure_class) | 0.92% (blouin_em) | 1.72 | 5.51x | 13.7% | 14.86x | 6.81x | **True** |

**Excluded on 8/8 settings.** The headline configuration is the one most generous to the mechanism; every other choice makes the shortfall larger, and using the lumen pH that Schmitt's own model uses makes it roughly 2.4x larger again.

| lumen pH | placement | `phi_lys(glioma)` | GM abs before | GM abs after | improved? | within 3x | adavosertib fold |
|---:|---|---|---:|---:|---|---:|---:|
| 4.7 | `method_aware` | 1.36% (normal_brain) | 1.995x | 1.624x | True | 5/6 | 1.241x |
| 4.7 | `method_aware` | 0.42% (collins_low) | 1.995x | 1.734x | True | 5/6 | 0.545x |
| 4.7 | `method_aware` | 3.08% (collins_high) | 1.995x | 1.827x | True | 5/6 | 2.512x |
| 4.7 | `all_tissues` | 1.36% (normal_brain) | 1.995x | 2.144x | False | 5/6 | 1.241x |
| 4.7 | `all_tissues` | 0.42% (collins_low) | 1.995x | 1.995x | False | 5/6 | 0.545x |
| 4.7 | `all_tissues` | 3.08% (collins_high) | 1.995x | 2.712x | False | 5/6 | 2.512x |
| 5.3 | `method_aware` | 1.36% (normal_brain) | 1.995x | 1.768x | True | 5/6 | 0.485x |
| 5.3 | `method_aware` | 0.42% (collins_low) | 1.995x | 1.902x | True | 4/6 | 0.312x |
| 5.3 | `method_aware` | 3.08% (collins_high) | 1.995x | 1.626x | True | 5/6 | 0.802x |
| 5.3 | `all_tissues` | 1.36% (normal_brain) | 1.995x | 1.995x | False | 5/6 | 0.485x |
| 5.3 | `all_tissues` | 0.42% (collins_low) | 1.995x | 1.995x | False | 4/6 | 0.312x |
| 5.3 | `all_tissues` | 3.08% (collins_high) | 1.995x | 1.995x | False | 5/6 | 0.802x |

The glioma improvement survives every setting of the method-aware placement and fails every setting of the physically wrong one. That is the right way round — the placement that is justified a priori is the one that helps — but it should still be read as one drug moving at n = 6, not as a validated mechanism.

## 8. Verdict and recommendation

### Recommendation: **default OFF**, and the module ships that way.

Nothing in the shipped pipeline imports `molecular.ion_trapping`, and the term returns exactly 1.0 for every drug in every tumour tissue because none has a measured lysosomal volume fraction. It is available to any caller who supplies one explicitly, as a documented mechanism with a stated parameter gap rather than an endorsed default.

### What this run establishes

1. **The chain's own correction was half right.** Moncayo's uniform-`R` proof genuinely does not apply to a between-tissue ratio — this run confirms that algebraically and as a unit test, and the term does move the GIST ratio where a scalar efflux term provably could not. But the magnitude claim in `uptake_impact.md` §10 compared the required 3.86% against a *swept* range rather than a *measured* one. Against the measurement it is a 3.09x overshoot.
2. **Ion trapping is excluded by a bound, exactly as hepatic uptake was.** Demanded 11.73x, available 4.47x at the measured liver value in a physically unattainable best case, 2.62x short. Uptake was 7.8x short. Two mechanisms have now been ruled out by literature ceilings without either being fitted.
3. **The sign is not supported and can invert.** No matched cancer-versus-normal lysosomal volume-fraction measurement exists; the two studies reporting volumes point the other way; and at the highest measured cancer-cell value the GIST prediction moves away from the data.
4. **The falsifiable side-condition fails.** The intestinal control moves from 1.11x to 4.87x at the same setting that improves the liver pairs, because gut and liver `phi_lys` differ by 3%.
5. **On glioma it helps, and it is the first thing in this chain that does.** That is worth recording precisely because the two previous runs were negatives on both boards. It rests on one drug, an unmeasured pKa and an assumed glioma `phi_lys`, and it is not claimed as evidence.
6. **The placement question has a physical answer, and it changes the glioma verdict.** Total-drug assays see `R`; extracellular microdialysis does not. Honouring that is what makes the glioma result positive, and ignoring it makes it negative.

### What this run does NOT establish

- **That lysosomal sequestration of imatinib is not real.** It plainly is: PMID 26108972 demonstrates it in GIST-T1 with bafilomycin A1, NH4Cl and LysoTracker colocalisation. What is shown is narrower — the TISSUE-LEVEL volume fractions the mechanism would need to explain the GIST tumour:normal contrast are not the ones that were measured.
- **Nothing at adequate power.** Six matched pairs with one independent tissue contrast; six drugs of which one moves.
- **Nothing about the tumour side, because nothing is measured there.** Every tumour number in this report is either an assumption or a breast cancer cell line standing in for a GIST.
- **Nothing about spatial structure.** This is a per-tissue constant. The Moncayo target needs a spatial FIELD of acidic-compartment content and that remains untested and unmeasured.

### Where the chain stands now

| candidate | status |
|---|---|
| passive reaction-diffusion transport | **falsified** for tissue partition — below the structural floor at any L (`gist_imatinib.md`) |
| scalar lysosomal ion-trapping, WITHIN a section | **falsified** — uniform R gives no dispersion at equilibrium (`moncayo_explant.md`) |
| active efflux (P-gp / BCRP) | **tested, does not help**; provably inert on a tissue ratio; wrong sign in liver (`efflux_impact.md`) |
| carrier-mediated uptake (OATP / OCT / OAT) | **tested, does not help**; named carrier falsified; class excluded by a literature bound (`uptake_impact.md`) |
| **tissue-differential lysosomal ion trapping** | **tested, does not rescue GIST**; excluded by a measured-parameter bound (2.62x short); breaks the intestinal control; sign unsupported. **Helps the glioma scorecard by one drug** (this run) |
| spatial acidic-compartment FIELD (within-tumour structure) | **OPEN**, unchanged; no measurement in usable units exists |

**Three mechanisms have now been named, implemented, cited and scored, and none of them closes the GIST residual.** Two were excluded by literature ceilings and one by a scalar-cancellation proof. That is a stronger conclusion than any of the three negatives on its own: the residual in the GIST tumour:normal contrast is not accounted for by any transporter or sequestration term that today's literature can parameterise, and the next step is **measurement rather than modelling**.

### The three measurements that would settle it

Named specifically, because "more data" is not a plan:

1. **A lysosomal volume fraction for GIST tissue, or for GIST-T1 cells, on a stated basis.** It does not exist. Collins 2021's method (live-cell confocal z-stacks normalised to measured cell volume) applied to GIST-T1 and to matched normal hepatocytes would decide the sign and the magnitude in one experiment.
2. **Forsby, Griph & Ivarsson 1983 (PMID 6137128).** Stereological lysosomal volume density in a human normal glial line versus a human malignant glioma line — exactly the glioma parameter this run had to assume. The paper is not open access, has no PMCID, was never archived and has three citations, so the four numbers were never quoted onward. This is an interlibrary-loan request, not a research programme.
3. **A normal-human liver lysosomal volume fraction.** Every value in the PBPK literature is rat, a textbook or a fit. Rohr 1976 (PMID 823698) and Roessner 1978 (PMID 654842) contain them and are not obtainable online.

---

*Generated by `tmesim/validation/ion_trapping_impact.py`. Additive only: no file under `core/`, `physics/`, `molecular/` (other than the new `molecular/ion_trapping.py`), `reconstruct/`, `sim/` or `treatment/` was modified, and the bit-exactness of the untouched delivery path is verified in §1.*
