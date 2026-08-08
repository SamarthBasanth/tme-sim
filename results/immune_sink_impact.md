# Immune (TAM) drug sequestration — implemented as a SPATIAL sink, and the spatial structure is the whole result

> ## CORRECTION NOTICE — the headline spatial number is smaller than this report originally claimed
>
> An audit of §5 found that the 2.0025x structure ratio was being quoted against a null that is too weak to be informative, and that a large part of the spread is a transport gradient with depth into tissue rather than immune biology. **The honest excess attributable to where the macrophages actually sit is ~1.08-1.21x, not 2.0x.** The effect is still there and it still clears an honest null — but the excess over that null is 0.08-0.21, against the 1.0 excess that "2.0x versus a well-mixed 1.0" advertised.
>
> The nulls, the floors, the depth confound and the CV diagnostic below are all computed by `tmesim/validation/immune_sink_impact.py`, so re-running the command in §9 reproduces this correction rather than reverting it. What changed, in one place:
>
> | claim as originally published | corrected |
> |---|---|
> | structure ratio 2.0025x "against a well-mixed model's exactly 1.0" | observed 2.0025x is real and reproduces to the digit. But the defensible comparator is a **spatially-aware** null, not a well-mixed one: clumpiness-matched permutation **1.6581**, 400 um block shuffle **1.8522**. **Excess = 1.208x and 1.081x.** |
> | the label-permutation null (mean 1.4145) is the evidence | that null under-clumps the sink field it is being compared against — CV **0.83** vs the observed **1.341** — so it loses almost automatically. It is **demoted** to a floor check (§5.3). |
> | the quintile gradient is an immune effect | binning the same tumour cells by **distance to tissue edge alone**, with no reference to immune density at all, already gives SR **1.6394** (§5.1a). |
> | monotonicity across all five quintiles is corroboration | **every** control is monotone, including uniform-in-tissue scatter. Monotonicity is not evidence (§5.4). |
> | p90/p10 spread 2.4078x against a well-mixed 1.0 | the well-mixed floor for that statistic is **1.4281**, not 1.0 (§5.1). |
> | (not previously stated) | **the statistic is INVALID on Visium/spot data**; the generator now REFUSES to compute it there (§5.6). |
>
> Nothing in §1-§4 or §6 changes. Those sections were not implicated by the audit.

**Verdict up front.** The immune compartment was pharmacologically inert: 11.4-19.5% of the segmented cells across the three real sections occupied space and consumed nothing. It now has a mechanism — a first-order drug loss proportional to the local MEASURED immune-cell area fraction — and the honest findings are three. **(i) The literature supports a DIRECTION and not a MAGNITUDE**, and for free doxorubicin the one direct measurement is NEGATIVE, so the magnitude is swept and never fitted. **(ii) The reduction is spatially structured rather than uniform, but the structure attributable to immune geometry is ~1.08-1.21x, not the 2.0x this report previously led with.** The raw quintile contrast is 2.0025x (densest quintile keeps 0.4291 of its exposure against 0.8592 in the sparsest), but most of that is reproduced by nulls that keep the macrophages' clumpiness (1.6581 and 1.8522) and by a depth-into-tissue gradient with no immune content at all (1.6394); against honest spatially-aware nulls the observed value still sits above their entire ranges (p = 0.0476, the 20-draw floor), which is what remains defensible. **(iii) It does not improve the measured human-glioma anchor**, so it ships DEFAULT-OFF.

| question | answer |
|---|---|
| drugs with a MEASURED free-drug TAM uptake ratio | **0 of 3** |
| mean tumour AUC, doxorubicin, S=0 -> S=20 | 11504.2 -> 5938.0 uM*s (0.516x) |
| therapeutic coverage at 1x IC50, same range | 0.4373 -> 0.1022 |
| predicted kill, same range | 0.2595 -> 0.1676 |
| **cross-drug RANKING across the whole sweep** | **INVARIANT** |
| spatial structure ratio, raw (sparsest / densest immune quintile) | 2.0025x |
| **same ratio, excess over an honest spatially-aware null** | **1.08-1.21x**, p = 0.0476 (= the 20-permutation floor) |
| same ratio from tissue-EDGE DISTANCE alone, no immune input | 1.6394x — see §5.1a |
| validity on Visium / spot transcriptomics | **INVALID — refused by the generator (§5.6)** |
| human glioma scorecard, GM abs fold error | 1.995x -> 1.995x at best swept S (0.0) |
| **recommendation** | **DEFAULT OFF** — see §7 |

---

## 1. Nothing outside the two new files was touched, and OFF is bit-identical

The escalation set by `results/uptake_impact.md` and `results/ion_trapping_impact.md` is followed exactly: **no shipped code path imports `tmesim.molecular.immune_sink`**, so the term is off BY CONSTRUCTION rather than by default. Files created: `tmesim/molecular/immune_sink.py`, `tmesim/validation/immune_sink_impact.py`, `tests/test_immune_sink.py`, this report. Files modified elsewhere: **0**.

A spatially-varying sink cannot be expressed by rescaling the AUC field afterwards — drug diffuses INTO a depleted region from its neighbours, and that non-locality is the entire microenvironment effect being measured. So the transport is RE-SOLVED. The shipped solver is not modified and not copied: its three numerical primitives (`_effective_diffusivity`, `_build_diffusion_operator`, `_advective_divergence`) are imported unchanged, exactly as `validation/moncayo_explant.py::bath_diffusion` already does, and the single functional difference is one line:

```python
k = np.full((ny, nx), drug.clearance + drug.binding) + sink_rate
```

Because `x + 0.0 == x` exactly in IEEE-754, the OFF path must be BIT-identical, and that is asserted rather than assumed:

- **36/36 solver comparisons bit-for-bit identical** (`.tobytes()` equality on `auc` and `final_concentration`, `np.array_equal` on `times` and the concentration history, `==` on the penetration index; never `allclose`) across 3 (clearance, binding) pairs x convection on/off x 2 timesteps x {`sink_rate=None`, `sink_rate=zeros`} plus the time-resolved-source branch. All identical: **True**.
- Every drug's LIBRARY strength is exactly 0.0, including a name not in the table: **True**.

And on the REAL flagship section the gate is end-to-end rather than synthetic: the `S = 0` arm of §4 is required to reproduce `sim.engine.simulate_response` bit-for-bit or the run aborts.

| drug | S=0 AUC bit-identical to the engine | S=0 kill identical | engine kill |
|---|---|---|---|
| doxorubicin | **True** | **True** | 0.259491 |
| temozolomide | **True** | **True** | 0.018352 |
| cisplatin | **True** | **True** | 0.436167 |

The checks are not vacuous — a positive sink genuinely bites, only ever lowers exposure, and is spatially non-uniform:

- a positive sink changes the AUC field: **True**
- it never RAISES exposure anywhere: **True**
- and the attenuation genuinely VARIES across tissue (range 0.0024, which would be 0 for a well-mixed sink): **True**

## 2. The literature — a DIRECTION, not a MAGNITUDE, and one direct NEGATIVE

This was searched before it was coded. The result is uncomfortable and it is the most important section of the report, because everything downstream is conditional on it.

| finding | what was measured | direction | PMID / DOI |
|---|---|---|---|
| Miller et al., Nat Commun 2015;6:8692 — **the best-quantified case** | single-cell imaging of a PLGA-b-PEG **nanoparticle** carrying a fluorescent Pt(IV) pro-drug. TAMs were **4% of cells but accumulated 30% of total intratumoral nanoparticle**; payload was **>3x elevated within ~15 um of a TAM** with gradients resolved over **<25 um**; macrophage depletion **LOWERED** intratumoral Pt 24% (5.0 -> 3.8 uM) and made tumours grow **faster** | **DEPOT — the opposite of a sink** | 26503691 / 10.1038/ncomms9692 |
| Islam et al., Nanotheranostics 2022;6(4):451-464 — **the direct free-doxorubicin measurement** | TAM-rich WEHI-164 fibrosarcoma. Liposomal formulations "shifted cellular drug uptake to TAM and to monocytic MDSC"; **"free doxorubicin cellular drug uptake was below detectable levels"** | **NEGATIVE for free drug** | 36105861 / 10.7150/ntno.75045 |
| Mayer et al., J Pharmacol Exp Ther 1997;280(3):1406-1414 | two murine fibrosarcomas differing **4-fold in TAM content**; TAM-rich tumours took 1.5-2x more liposomal doxorubicin, but CD11b+ cells correlated **poorly** with doxorubicin fluorescence and the authors concluded TAM uptake does **not** explain it (attributing it to TAM-driven vascular permeability) | **NEGATIVE for TAM uptake** | 9067330 |
| Zhitomirsky & Assaraf, Oncotarget 2015;6(2):1143-1156; Drug Resist Updat 2016;24:23-33 | doxorubicin is a hydrophobic weak base (amino-sugar pKa ~8.3) passively ion-trapped into acidic organelles; sequestered drug is held away from its nuclear target | **sink, mechanistic, no TAM-specific ratio** | 25544758 / 10.18632/oncotarget.2732; 26830313 / 10.1016/j.drup.2015.11.004 |
| Aghai & Tokes, Cancer Chemother Pharmacol 1990;25(4):295-298 | doxorubicin compartmentalisation across normal human marrow lineages: the nuclei of **macrophages and activated phagocytes were the LEAST doxorubicin-accessible of any lineage**, i.e. in macrophages the drug is held in the cytoplasmic/vesicular compartment | **sink, single-cell, no ratio** | 2295115 / 10.1007/BF00684889 |

Read the first two rows together, because they are the ones that decide the framing. **The well-quantified TAM drug-handling literature is about CARRIERS, not free drug, and where it is quantified the sign is a DEPOT.** The one experiment that measured free doxorubicin in TAMs found it below detection. So a TAM SINK for free cytotoxics rests on mechanism (lysosomal ion trapping in the most lysosome-rich cell in the niche) plus a single-cell compartmentalisation observation from 1990 — a direction, with no number attached.

**Consequently `immune_sink_strength_for(drug)` returns exactly 0.0 for every drug**, the same identity `uptake_factor` and `ion_trapping` adopted, and the magnitude below is a SWEPT parameter. Per drug:

| drug | measured ratio | direction the literature supports | why |
|---|---|---|---|
| doxorubicin | **NONE** | sink (mechanistic, unquantified for free drug) | hydrophobic weak base (amino-sugar pKa ~8.3); passive ion-trapping into acidic organelles, of which macrophages carry the largest complement. In normal human marrow the n... |
| temozolomide | **NONE** | no evidence in either direction | small (194 Da), neutral, non-basic, freely diffusing, and it is a PRO-DRUG that decomposes spontaneously at physiological pH rather than being taken up and held. There is... |
| cisplatin | **NONE** | depot, not sink, and only when nano-encapsulated | a small neutral/aquated platinum species with no basic centre, so no ion-trapping route. The one quantified TAM result for a platinum agent is Miller 2015's Pt(IV) PRO-DR... |

Note the contrast structure this gives, which is the point of including the other two drugs: **doxorubicin** has a mechanism and no number, **cisplatin**'s only quantified TAM result is a nanoparticle DEPOT whose sign is opposite, and **temozolomide** — small, neutral, non-basic, a pro-drug that decomposes spontaneously — has no proposed route at all. A sink applied equally to all three is therefore already known to be wrong; it is applied equally anyway, because the alternative is to invent three different unmeasured numbers.

### 2.1 The functional form and its units

```
  dC/dt = div(D_eff grad C) - div(vC) + S - (clearance + binding) C
                                          - k_sink(x, y) C        <-- NEW

  k_sink(x, y) = k_max * phi_imm(x, y)                 [1/s]
```

| symbol | meaning | units | where it comes from |
|---|---|---|---|
| `phi_imm` | local IMMUNE-cell area fraction | dimensionless, [0,1] | **MEASURED** — Gaussian KDE of the real immune-cell centroids (cells/um^2) times the per-cell cross-section `pi r^2` (um^2/cell) |
| `r` | immune-cell cross-sectional radius | um | DECLARED 7.5 um (macrophage/microglial soma ~15 um across) |
| bandwidth | footprint over which one macrophage depletes drug | um | DECLARED 25 um, **literature-anchored**: Miller 2015 resolved TAM payload gradients over <25 um (PMID 26503691) |
| `k_max` | loss constant of a voxel packed ENTIRELY with macrophages | 1/s | **SWEPT**, never fitted |
| `S` | `k_max / (clearance + binding)` | dimensionless | the sweep axis: "how many times the drug's own tissue sink does a fully-macrophage voxel add" |

Only the PRODUCT `k_max * phi_imm` is physically meaningful, so `r` and `k_max` trade off exactly against each other — which is precisely why `k_max` is swept rather than asserted, and why no number here should be read as a calibrated rate. **Zero parameters were fitted to any scorecard.**

## 3. The measured immune field, and where the literature anchor lands on it

Immune-cell area fraction on all three real sections at dx = 12 um, from the frozen `reconstruct_real` default path (matched-null classifier).

| section | cells | IMMUNE | IMMUNE % | grid | mean phi | median | p90 | p99 | max |
|---|---|---|---|---|---|---|---|---|---|
| `brain_319` | 40,887 | 7,965 | 19.48% | 424x675 | **0.0495** | 0.0282 | 0.1258 | 0.2962 | 1.0000 |
| `io_480` | 816,769 | 93,004 | 11.39% | 905x956 | **0.1833** | 0.1530 | 0.3417 | 0.6923 | 1.0000 |
| `gse324268_366` | 851,509 | 107,417 | 12.61% | 1196x880 | **0.1925** | 0.1802 | 0.3110 | 0.4919 | 1.0000 |

The distribution is strongly right-skewed on every section — the p99 is 3-6x the mean — which is the reconstruction saying the immune compartment is CLUMPED, not spread. That clumping is what §5 tests, and it is independently corroborated: `results/classifier_fix.md` §4 measures nearest-neighbour IMMUNE label enrichment of 1.41-1.83x above random on these same sections under the shipping classifier.

**Converting the one literature anchor into this model's units.** Miller 2015 gives a fraction `f` of intratumoral drug that ends up inside macrophages (0.30 from the 4%-of-cells/30%-of-drug partition; 0.24 from the depletion experiment). Two first-order losses in parallel split the drug in proportion to their rate constants, so `f = <k_sink>/(<k_sink> + k_base)` and therefore `S = f / ((1 - f) * <phi_imm>)`. Run through each section's OWN measured field:

| section | mean phi | S at f = 0.24 (depletion) | S at f = 0.30 (partition) |
|---|---|---|---|
| `brain_319` | 0.0495 | **6.38** | **8.65** |
| `io_480` | 0.1833 | **1.72** | **2.34** |
| `gse324268_366` | 0.1925 | **1.64** | **2.23** |

So the anchor lands at **S ~ 6.4-8.7** on the flagship section, which is why the sweep runs 0.0, 1.0, 2.0, 5.0, 10.0, 20.0 — it BRACKETS the anchor rather than sitting on it. **This is a bracketing calculation, not a calibration.** The anchor is a polymeric nanoparticle in a murine fibrosarcoma; the drugs here are free small molecules in human glioma; and in the source experiment the effect ran in the DEPOT direction. It is used only to answer "how large could a diverted fraction plausibly be at all", and the answer is that even the most TAM-avid formulation ever measured diverted 24-30%.

## 4. The sweep on the real flagship section

**Section.** `data_cache/xenium_gbm` at dx = 12 um, grid 424x675, 23,812 tumour cells. Perivascular-sink alpha = 148.61, calibrated once on doxorubicin exactly as every flagship run does. ONE reconstruction and ONE oxygen + IFP prologue are shared by every arm, so the comparison isolates the sink. `_PDE` is imported from `validation/delivery_metrics.py` so this and every flagship delivery number solve the identical problem.

**Why this section and not all three.** Every absolute AUC / coverage / kill number in `results/` was produced on this section at this dx, so it is the only one against which a delta is interpretable. The immune FIELD is measured on all three (§3), and it is worth reading §3 alongside this table: `io_480` and `gse324268_366` carry mean immune area fractions 3.9x the flagship's, because they are far denser sections, so the same S would produce a substantially larger effect there. The levels below are flagship-specific; the SHAPE of the S-dependence is not.

| drug | IC50 (uM) | S | diverted fraction f | mean tumour AUC (uM*s) | mean exposure (uM) | coverage @1x IC50 | predicted kill | penetration index |
|---|---|---|---|---|---|---|---|---|
| doxorubicin | 0.5 | **0** | 0.000 | 11504.2 | 0.5326 | 0.4373 | 0.2595 | 0.3044 |
| doxorubicin | 0.5 | 1 | 0.047 | 10752.4 | 0.4978 | 0.3942 | 0.2507 | 0.3163 |
| doxorubicin | 0.5 | 2 | 0.090 | 10153.6 | 0.4701 | 0.3568 | 0.2417 | 0.3224 |
| doxorubicin | 0.5 | 5 | 0.198 | 8854.2 | 0.4099 | 0.2689 | 0.2231 | 0.3273 |
| doxorubicin | 0.5 | 10 | 0.331 | 7490.4 | 0.3468 | 0.1847 | 0.1991 | 0.3218 |
| doxorubicin | 0.5 | 20 | 0.498 | 5938.0 | 0.2749 | 0.1022 | 0.1676 | 0.3025 |
| temozolomide | 606.047 | **0** | 0.000 | 765899.6 | 35.4583 | 0.0000 | 0.0184 | 0.5107 |
| temozolomide | 606.047 | 1 | 0.047 | 722992.1 | 33.4719 | 0.0000 | 0.0172 | 0.5291 |
| temozolomide | 606.047 | 2 | 0.090 | 686692.7 | 31.7913 | 0.0000 | 0.0159 | 0.5468 |
| temozolomide | 606.047 | 5 | 0.198 | 603401.8 | 27.9353 | 0.0000 | 0.0139 | 0.5971 |
| temozolomide | 606.047 | 10 | 0.331 | 511609.4 | 23.6856 | 0.0000 | 0.0117 | 0.6067 |
| temozolomide | 606.047 | 20 | 0.498 | 403457.0 | 18.6786 | 0.0000 | 0.0086 | 0.5793 |
| cisplatin | 24.5503 | **0** | 0.000 | 1821653.3 | 84.3358 | 0.9955 | 0.4362 | 0.4283 |
| cisplatin | 24.5503 | 1 | 0.047 | 1712977.8 | 79.3045 | 0.9949 | 0.4305 | 0.4459 |
| cisplatin | 24.5503 | 2 | 0.090 | 1622919.2 | 75.1352 | 0.9942 | 0.4246 | 0.4608 |
| cisplatin | 24.5503 | 5 | 0.198 | 1420895.3 | 65.7822 | 0.9890 | 0.4081 | 0.4827 |
| cisplatin | 24.5503 | 10 | 0.331 | 1203019.2 | 55.6953 | 0.9628 | 0.3837 | 0.4891 |
| cisplatin | 24.5503 | 20 | 0.498 | 950039.8 | 43.9833 | 0.8383 | 0.3461 | 0.4726 |

Relative to OFF, per drug, across the swept range:

| drug | AUC at S_max / AUC at S=0 | coverage change (pp) | kill change (pp) |
|---|---|---|---|
| doxorubicin | **0.5162x** | -33.52 | -9.19 |
| temozolomide | **0.5268x** | +0.00 | -0.98 |
| cisplatin | **0.5215x** | -15.71 | -9.00 |

**Read the three drugs against each other — the asymmetry is the useful part, and it is not about the sink.** The AUC column is nearly identical across drugs, because the sink is applied with the same `S` to all three and the transport is close to linear in it. What differs is how that same fractional exposure loss TRANSLATES, and the translator is each drug's position on its own Hill curve:

| drug | mean exposure / IC50 at S=0 | AUC retained at S_max | **kill retained** at S_max | coverage @1x IC50, S=0 -> S_max |
|---|---|---|---|---|
| doxorubicin | 1.0652 | 0.5162x | **0.6457x** | 0.4373 -> 0.1022 |
| temozolomide | 0.0585 | 0.5268x | **0.4668x** | 0.0000 -> 0.0000 |
| cisplatin | 3.4352 | 0.5215x | **0.7936x** | 0.9955 -> 0.8383 |

A drug sitting FAR BELOW its IC50 is on the steep, near-proportional part of the Hill curve, so a given fractional exposure loss costs it nearly the same fraction of its kill. A drug sitting AT or ABOVE its IC50 for much of the tissue is partly saturated, so the same fractional exposure loss costs it proportionally less kill — but costs it a great deal of COVERAGE, because coverage is a threshold count and the threshold is exactly where its cells are. Temozolomide is the extreme of the first case (its coverage at 1x IC50 is 0.0000 at every strength — 606 uM is simply unreachable in this tissue, so the metric cannot move) and doxorubicin of the second (coverage collapses while kill is comparatively defended). **The practical consequence: which metric is sensitive to a delivery insult depends on the drug, and reporting only one of AUC / coverage / kill would have hidden two-thirds of this effect.**

### 4.1 The cross-drug ranking

| S | ranking by predicted kill | by coverage @1x IC50 | by mean tumour AUC |
|---|---|---|---|
| 0 | cisplatin > doxorubicin > temozolomide | cisplatin > doxorubicin > temozolomide | cisplatin > temozolomide > doxorubicin |
| 1 | cisplatin > doxorubicin > temozolomide | cisplatin > doxorubicin > temozolomide | cisplatin > temozolomide > doxorubicin |
| 2 | cisplatin > doxorubicin > temozolomide | cisplatin > doxorubicin > temozolomide | cisplatin > temozolomide > doxorubicin |
| 5 | cisplatin > doxorubicin > temozolomide | cisplatin > doxorubicin > temozolomide | cisplatin > temozolomide > doxorubicin |
| 10 | cisplatin > doxorubicin > temozolomide | cisplatin > doxorubicin > temozolomide | cisplatin > temozolomide > doxorubicin |
| 20 | cisplatin > doxorubicin > temozolomide | cisplatin > doxorubicin > temozolomide | cisplatin > temozolomide > doxorubicin |

**The ranking by predicted kill is INVARIANT across the entire swept range.** That is the expected structure and it is the same one `results/cell_state_impact.md` §3.3 found for the cell-state route: the sink is applied with the SAME strength to all three drugs, so it very largely cancels out of a comparison and does not cancel at all out of a level. A claim of the form *"cisplatin beats doxorubicin beats temozolomide on this tissue"* survives every sequestration strength tested; a claim about an absolute kill fraction does not.

## 5. THE SPATIAL QUESTION — is the reduction structured, or just smaller?

This is the section the whole module exists for, and it is the section the audit corrected. **A uniform reduction is pharmacologically boring: a compartmental model reproduces it exactly by lowering one rate constant.** The question is whether tumour cells sitting in macrophage-dense neighbourhoods are measurably worse off than tumour cells that are not — because that is a statement only a spatial model can make.

**Read the answer first, then the working.** The raw quintile contrast is 2.0025x and reproduces to the digit. But three things have to be subtracted from it before it means what this report originally claimed it meant: a depth-into-tissue transport gradient that carries 1.6394x on its own with no immune input (§5.1a), a null that under-clumps the sink field and therefore loses without contest (§5.3), and a monotonicity argument that every control also satisfies (§5.4). What survives all three is an excess of **1.08-1.21x** over honest spatially-aware nulls, above their entire ranges at p = 0.0476 (the 20-draw floor).

**Configuration.** `doxorubicin` at S = 10 (the sweep point nearest the Miller anchor), flagship section `data_cache/xenium_gbm`, dx = 20 um, transvascular exchange off, 23,812 tumour cells and 7,965 immune cells. The statistic is the per-tumour-cell **attenuation** `a_i = AUC_sink,i / AUC_base,i`, from the SAME prologue, source field and velocity field — the only difference between the two solves is `k_sink`.

The observed statistic is computed twice by two independent code paths in this file (`spatial_structure` for §5.5's mesh check and `_spatial_context` for the nulls below); they agree to 0.000000 absolute, which is the internal consistency check on everything that follows.

**A declared cost decision, and it is priced.** Every null draw is a full transport re-solve, so the nulls run at dx = 20 um rather than the flagship 12 um (~2.8x fewer unknowns). Because this project treats mesh-dependence as a first-class failure mode, the SAME observed statistic is also computed at the flagship spacing and both are reported in §5.5 — so the reader can see directly whether the structure is a property of the tissue or of the grid.

### 5.1 Tumour-cell exposure binned by LOCAL IMMUNE DENSITY

| quintile of local phi_imm | n tumour cells | mean phi_imm | mean AUC, sink OFF | mean AUC, sink ON | **mean attenuation** |
|---|---|---|---|---|---|
| 1 (sparsest) | 4,763 | 0.0007 | 8398.1 | 7379.7 | **0.8592** |
| 2 (.) | 4,762 | 0.0107 | 10031.7 | 7848.5 | **0.7708** |
| 3 (.) | 4,762 | 0.0338 | 10642.0 | 7545.7 | **0.6937** |
| 4 (.) | 4,762 | 0.0675 | 11806.0 | 6966.1 | **0.5557** |
| 5 (densest) | 4,763 | 0.1586 | 15941.5 | 7006.8 | **0.4291** |

**Structure ratio** (mean attenuation in the sparsest quintile / mean attenuation in the densest) = **2.0025x**. This number is correct and reproduces to the digit; what it is compared against is where this report previously went wrong. Read §5.1a, §5.3 and §5.4 before quoting it.

**The two floors, MEASURED rather than asserted.** The claim "a well-mixed model gives exactly 1.0" is not assumed here: a genuinely uniform sink field of the same tissue-mean magnitude is built, re-solved, and binned identically.

| statistic | observed | well-mixed floor, MEASURED | what the excess over the floor actually is |
|---|---|---|---|
| structure ratio (sparsest/densest quintile) | **2.0025** | **0.9343** | the analytic 1.0000 is a FAIR comparator for this binned statistic — the measured floor sits slightly below it, so quoting 1.0 is if anything mildly conservative |
| per-cell attenuation spread (p90/p10) | **2.4078** | **1.4281** | **the analytic 1.0 is NOT the right floor here.** 2.4078 must be read against 1.4281, i.e. an excess of 1.686x |

The p90/p10 figure of 2.4078x (over a full range of 0.0000 to 0.9915, against a section-mean attenuation of 0.6617) was previously quoted as though its floor were 1.0. **It is 1.4281.** Any external use of the p90/p10 number must carry that floor with it.

The well-mixed floor is the easy comparator, and it is not the informative one. A uniform sink is not a serious rival hypothesis — the serious rivals are "the macrophages are clumped, but not clumped *where the data says*" (§5.3) and "this is depth into tissue, not immune density" (§5.1a).

#### 5.1a A large part of the spread is DEPTH INTO TISSUE, not immune biology

This is the largest single correction in this report and it was not tested in the original version.

Take the same 23,812 tumour cells, the same two solves, and bin them by **distance to the tissue edge alone** — no reference to immune density anywhere in the binning:

| binning variable | sparsest/shallowest bin | densest/deepest bin | **structure ratio** |
|---|---|---|---|
| local immune density `phi_imm` (table above) | 0.8592 | 0.4291 | **2.0025** |
| **distance to tissue edge, no immune input at all** | **0.7850** | **0.4789** | **1.6394** |

So a variable containing **zero immune information** already reaches 1.6394 of the headline 2.0025 — 82% of the ratio. The two binning variables are correlated: **corr(phi_imm, edge-distance) = 0.372** (Spearman 0.581).

**The depth gradient comes from the measured tissue, not from the kernel.** The segmented cell-density envelope rises **2.5x** from the outer 25 um band into the interior. Cells deep in the section are simply in denser tissue, which carries both more macrophages and different transport. The uniform-in-tissue scatter control (§5.4) passes the identical KDE and reaches only 1.2126, so the KDE is not manufacturing the gradient by itself.

**And the published label-permutation null cannot remove this confound, because it PRESERVES it.** Permuting the immune label among measured cell positions keeps every immune cell on a real cell, so it inherits the same density envelope: the permuted fields carry corr(phi_imm, edge-distance) = **0.354** against the observed 0.372. The null reproduces the depth gradient it should be controlling for.

**How much survives after stratifying it out.** Splitting the tumour cells into five edge-distance strata and computing the immune-density structure ratio *within* each:

| edge-distance stratum | 1 (shallowest) | 2 | 3 | 4 | 5 (deepest) |
|---|---|---|---|---|---|
| within-stratum structure ratio | **1.33** | **1.56** | **1.76** | **2.03** | **1.56** |

An immune-density contrast persists at every depth, so the effect is not purely depth. But it is 1.33-2.03x within strata rather than 2.0025x overall, and it is not monotone in depth — which is what one expects when two partly-collinear variables are being separated on a single section.

### 5.2 Exposure deficit versus DISTANCE to the nearest macrophage

The same effect read the other way round, and on the length scale Miller 2015 measured (>3x payload elevation within ~15 um of a TAM; gradients resolved over <25 um):

| distance to nearest immune cell | n tumour cells | mean AUC, sink OFF | mean attenuation |
|---|---|---|---|
| 0-10 um | 360 | 25332.1 | **0.4493** |
| 10-20 um | 4,859 | 13780.3 | **0.5144** |
| 20-30 um | 5,737 | 11920.6 | **0.5795** |
| 30-50 um | 6,055 | 10820.7 | **0.6797** |
| 50-100 um | 5,588 | 9233.7 | **0.8206** |
| 100-inf um | 1,213 | 7431.1 | **0.8814** |

**This table is NOT independent corroboration of Miller's length scale, and must not be read as such.** The kernel bandwidth was CHOSEN as 25 um from that same measurement (§2.1), so the model was told what length scale to have and the fact that it has it is arithmetic, not evidence.

**Nor is the monotonicity of this column evidence of anything** — §5.4 shows every control reproduces monotonicity, including controls with no measured immune structure in them at all. The original version of this report leaned on that monotonicity; it should not have.

What the table does show that is not circular is that the deficit does not saturate until well beyond the bandwidth (drug also diffuses in from further away, which is exactly the non-local effect an AUC rescale could not have produced), and that the worst-off cells are a small minority — 360 of 23,812 tumour cells lie within 10 um of a macrophage. Note also that the sink-OFF AUC column falls with distance too: cells near macrophages are in denser, better-perfused tissue and started with MORE drug, so the attenuation column is the right comparison and the raw exposure column would point the other way. That same tissue-density envelope is the confound quantified in §5.1a.

### 5.3 The nulls — the published one was too weak, and here are the honest ones

**What was wrong.** The original version of this report tested the observed structure against ONE null: permuting the IMMUNE LABEL among the measured cell positions. That null does hold the tissue outline, the cell layout, the vessel geometry, the ECM field and the total sink magnitude fixed — but it destroys the immune compartment's clumping entirely, and clumping is the very quantity the structure ratio measures.

**Measured, not argued.** The sink field the label-permutation null produces has coefficient of variation **0.832** (range 0.813-0.856). The observed sink field has CV **1.341**. A null whose field is 38% less variable than the observed field will lose a variance-driven contest almost regardless of whether any real spatial biology exists — it is not being asked a fair question.

And the dependence is direct rather than incidental. Sweeping the clump size of the permutation, holding the immune COUNT and therefore the total sink magnitude fixed:

| clump size (cells) | 1 | 5 | 20 | 50 |
|---|---|---|---|---|
| sink-field CV | 0.83 | 1.29 | 1.84 | 2.03 |
| resulting structure ratio | **1.41** | **1.67** | **2.24** | **2.67** |

**The structure ratio very largely tracks sink-field heterogeneity.** Any null that does not match the observed CV of 1.341 is measuring CV, not geometry.

#### The primary comparators

Two nulls fix this. Both keep **every immune cell on a real measured cell position** — nothing is scattered into empty space, and nothing leaves the tissue mask:

- **clumpiness-matched permutation** — the immune label moves in clumps of 5 mutually-nearest MEASURED cells instead of one cell at a time, so the CV asymmetry above is largely removed (1.291 against the observed 1.341) and only the *location* of the clumps varies.
- **400 um block shuffle** — the section is tiled into 400 um blocks and each block's observed immune COUNT is reassigned to another block, which draws that many of its own measured cells. The block-scale count heterogeneity of the real section is inherited exactly and only relocated.

20 draws each, each one a full transport re-solve.

| null | rank | null mean SR | null range | sink-field CV | observed | one-sided p | **observed / null mean** |
|---|---|---|---|---|---|---|---|
| **clumpiness-matched permutation** | **PRIMARY** | **1.6581** | 1.6053-1.6854 | 1.291 | 2.0025 | 0.0476 | **1.208x** |
| **400 um block shuffle** | **PRIMARY** | **1.8522** | 1.7483-1.9106 | 1.130 | 2.0025 | 0.0476 | **1.081x** |
| label permutation (as originally published) | DEMOTED — floor check only | 1.4145 | 1.3879-1.4448 | 0.832 | 2.0025 | 0.0476 | 1.416x |
| jitter reshuffle (re-draw the Visium jitter only) | control, see below | 2.0177 | 2.0060-2.0317 | 1.329 | 2.0025 | 1.0000 | 0.992x |
| spot emulation (the whole Visium synthesis) | control, see below | 1.9820 | 1.9698-1.9899 | 1.413 | 2.0025 | 0.0476 | 1.010x |
| torus shift | **UNUSABLE on this mask** — see below | -- | -- | 1.780 | 2.0025 | -- | -- |

**The honest number is the last column: 2.0025 / 1.6581 = 1.208x and 2.0025 / 1.8522 = 1.081x.** That — an excess of roughly 1.08 to 1.21x — is what is attributable to where the macrophages actually are, once a null is allowed to be as clumped as the data. It is not 2.0x. This report previously said 2.0x, and that was quoting a raw ratio against a floor of 1.0 as though the floor were the alternative hypothesis.

**What still holds, and it is not nothing.** The observed 2.0025 sits **above the entire range of both primary nulls** (1.6053-1.6854 and 1.7483-1.9106), at p = 0.0476. A clumpiness-matched null and a block-shuffle null both fail to reach the observed value in 20 out of 20 draws. Putting the same number of equally-clumped macrophages *somewhere else* on real measured cell positions produces a measurably weaker exposure gradient than putting them where the data says they are. A 1.08-1.21x effect that survives an honest null is worth more than a 2.0x effect resting on a weak one.

#### The platform controls — what happens if the positional detail goes away

**Jitter reshuffle — re-draw the Visium jitter and nothing else.** Every measured immune cell is displaced by a uniform-in-disc offset of at most 27.5 um, the exact amplitude `reconstruct/real_visium.py::_place_cells` applies (`jitter_um = min(27.5, 0.45 * dx_um)`). Nothing about which neighbourhood a macrophage is in changes. On a SPOT reconstruction this control is a no-op by construction — the positions were generated that way. Here, 20 draws give **2.0177** (range 2.0060-2.0317, CV 1.329 against the observed 1.341): it **reproduces the observed value** (2.0177 against 2.0025, 1.008x). Positional detail below the jitter amplitude therefore carries essentially none of the statistic: the structure ratio is built from immune density at scales COARSER than ~27.5 um, not from it. This control does not discriminate, and nothing is rested on it.

**Spot emulation — the whole Visium synthesis.** Snapping every measured immune cell to its 100 um lattice centroid and re-jittering it (literally what `_place_cells` does, keeping the per-spot immune COUNT and destroying everything finer) gives **1.9820** at CV 1.413: it **reproduces the observed value** (1.9820 against 2.0025, 0.990x). Positional detail below the 100 um spot pitch therefore carries essentially none of the statistic: the structure ratio is built from immune density at scales COARSER than the 100 um spot pitch, not from it. This control does not discriminate, and nothing is rested on it.

**This does NOT reproduce a claim the §5 audit made, and the audit's version should be dropped rather than repeated.** The audit reported a jitter reshuffle at 1.4188 on this section and read its collapse as the decisive evidence separating Xenium from Visium. Re-implemented here — displace every measured immune cell by at most 27.5 um, change nothing else — the control returns 2.0177, i.e. the observed value, and the full spot synthesis returns 1.9820 if run as well. Either the audit's control was a different operation or its number does not reproduce; either way **this report does not rest the Xenium result on a jitter control**, and no downstream document should. What the Xenium result rests on is the two primary nulls in the table above, which do not involve the platform question at all. The scope limit in §5.6 rests on PROVENANCE — spot platforms never localise an immune cell — and not on this control.

#### Two methodological notes that must travel with these numbers

1. **Torus / Moran-style nulls are unusable on this tissue mask and nothing rests on them.** Circularly shifting the immune point pattern lands only **66.5%** of shifted cells back inside tissue (over 20 random shifts; range 52.6-87.5%). The surviving field has mean `phi_imm` at **0.65x** the observed and CV inflated to **1.33x** the observed — i.e. the shift simultaneously weakens and sharpens the sink, which biases the structure ratio upward for entirely mechanical reasons. **This is measured here as a FEASIBILITY diagnostic and the torus null is never scored.** The clumpiness-matched permutation and the block shuffle are the defensible nulls, precisely because both keep every immune cell on a real measured position inside the mask.
2. **p = 0.0476 is the 20-draw floor** (1/21). It appears in every row above for that reason and **must never be quoted as a resolved p-value**. It means only "not exceeded in 20 draws". The informative columns are the null range and the observed/null ratio.

**One control that is not about structure at all.** Mean attenuation is a magnitude, not a structure statistic: observed 0.6617 against a label-permutation mean of 0.6091 and a clumpiness-matched mean of 0.6338. The same number of macrophages removes roughly the same total amount of drug wherever they sit. A compartmental model gets that quantity right by construction; what it cannot express is **which** tumour cells lose the drug. That distinction is still the point of the module — it is the *size* of the spatial term that has been corrected, not its existence.

### 5.4 Monotonicity across the quintiles is NOT evidence, and was previously treated as though it were

The original report pointed to the attenuation falling smoothly across all five immune quintiles (0.8592 -> 0.4291) as corroboration. It is not. Every control tested is also monotone across all five quintiles:

| field | quintile 1 | quintile 5 | monotone across all 5? | ratio |
|---|---|---|---|---|
| **observed** | 0.8592 | 0.4291 | yes | 2.0025 |
| label permutation (as originally published) | 0.7045 | 0.4980 | **yes** | 1.41 |
| clumpiness-matched permutation | 0.8042 | 0.4850 | **yes** | 1.66 |
| 400 um block shuffle | 0.8242 | 0.4451 | **yes** | 1.85 |
| jitter reshuffle (re-draw the Visium jitter only) | 0.8612 | 0.4268 | **yes** | 2.02 |
| spot emulation (the whole Visium synthesis) | 0.8651 | 0.4365 | **yes** | 1.98 |
| uniform-in-tissue scatter | 0.6990 | 0.5764 | **yes** | 1.21 |

Uniform-in-tissue scatter contains **no measured immune structure whatsoever** and is still monotone. Monotonicity is a generic consequence of binning a smooth field by its own local value — it would be surprising only in its absence. **Only the RATIO carries information, and only against a null that matches the observed field's clumpiness (§5.3).** Any downstream document quoting "monotone across all five quintiles" as support should drop that phrasing.

### 5.5 Is the structure a property of the tissue, or of the mesh?

`results/state_validation.md` and `results/cell_state_impact.md` both turn on the same question, so it is asked here rather than assumed. The identical observed statistic, on the identical section and drug at the identical strength, at two grid spacings:

| dx (um) | structure ratio (sparsest/densest quintile) | attenuation p90/p10 | mean attenuation |
|---|---|---|---|
| 12 (flagship) | **1.9683** | 2.2505 | 0.6741 |
| 20 (nulls run here) | **2.0025** | 2.4078 | 0.6617 |

Swing across the mesh: **1.017x**. That is small enough that the structure is a property of the measured cell geometry rather than of the discretisation, which is the bar this repo applies to every spatial claim. Note that mesh-stability is a necessary check, not a sufficient one: a discretisation-stable number can still be a depth artefact (§5.1a) or be quoted against the wrong null (§5.3), and here it was both.

### 5.6 SCOPE LIMIT — this statistic is INVALID on Visium / spot transcriptomics, and the generator now REFUSES to compute it there

This is a hard limit, not a caution. **The structure ratio computed on Visium or any other spot-based platform is meaningless, because the platform does not measure immune-cell positions at all.**

**The mechanism, read off the reconstruction code.** `tmesim/reconstruct/real_visium.py::_place_cells` builds the cell population from a deconvolved spot table by synthesising a fixed **12 pseudo-cells per ~100 um spot** (`_CELLS_PER_SPOT = 12`), apportioned across types by the deconvolved fractions, and placing all of them **at the spot centroid with uniform-in-disc jitter of at most 27.5 um** (`jitter_um = min(27.5, 0.45 * dx_um)`). There is no measured position for any individual immune cell anywhere in that pipeline. The only spatial information present is the per-spot immune FRACTION on a 100 um lattice; everything finer is drawn from a random number generator.

**Why every plausible null is therefore a no-op.** The sink field is a 25 um-bandwidth KDE. Sub-spot jitter of <=27.5 um is at or below that bandwidth, so re-drawing it — or reshuffling labels within a spot — changes the KDE almost not at all. Any control that preserves the per-spot immune count is, by construction, comparing a field to itself.

**The refusal is in the code, not in this prose.** `validation/immune_sink_impact.py::require_measured_cell_positions` raises `SpotDataError` for any reconstruction whose own metadata says the positions were synthesised from spots (`metadata['real_pipeline']['source']`, `metadata['n_spots']`); the shared context builder every null family runs through calls it first, so no null can bypass it, and an unstamped reconstruction of unknown provenance is refused too. It reads the reconstruction's metadata rather than the sample PATH, so a renamed section cannot get past it. `tests/test_immune_sink_nulls.py` pins the refusal.

**And the cost of ignoring the limit is measured here rather than asserted.** Passing the measured single-cell immune positions through the identical spot synthesis (§5.3, spot emulation) returns **1.9820** against the observed 2.0025 — 0.990x. **A spot-resolution field reproduces essentially the entire claimed effect.** Read that as a warning and not as a reassurance: it means the structure ratio is carried by immune density at the spot scale and coarser, so on a platform where the per-spot immune content is a DECONVOLUTION ESTIMATE rather than a count of localised cells, the statistic will report whatever the deconvolution wrote — at full apparent strength, and with no null available on that platform able to tell the two apart, because every control that preserves per-spot counts is a no-op there. Any use of a Visium-derived structure ratio in a deck, a landing page or a customer-facing claim is unsupported.

**This does not retract the single-cell result, and it does not rest on a platform control either.** What that result rests on is the two primary nulls of §5.3, which never leave the measured section: an equally-clumped rearrangement of the same macrophages onto other measured positions does not reach the observed value in any draw. The scope limit here is about PROVENANCE — a spot platform never localises an immune cell, so there is nothing on it for such a null to rearrange — and the generator enforces it by refusing to run.

## 6. Does it improve the measured human-glioma anchor? NO.

`results/glioma_regional_pk.md` is the only measured scorecard this term can reach: 6 drugs with measured human glioma tissue:plasma partitions, currently at a geometric-mean absolute fold error of 1.995x. At quasi-steady state a parallel loss multiplies the partition by `(1 - f)`:

```
  before:  R = clip( bbb_penetrance + L * D_sample, 0, 1 )
  after :  R = clip( bbb_penetrance + L * D_sample, 0, 1 ) * (1 - f)
```

**Modality is not optional here.** Microdialysis measures FREE extracellular drug, which is exactly what a sink into macrophages removes, so the factor applies. Tissue homogenate measures TOTAL drug INCLUDING what is inside the macrophages, so sequestration does not lower it at all and the factor must not be applied. This is the same free-versus-total bookkeeping `results/ion_trapping_impact.md` §6 established, used here with the opposite sign. Both variants are shown because the wrong one is what makes the choice load-bearing.

Applied at each section's own mean phi_imm = 0.0495:

| S | diverted fraction f | modality-aware GM abs fold | (physically wrong) applied to all GM abs fold | within 3x | within 10x | worst |
|---|---|---|---|---|---|---|
| **0 (OFF)** | 0.000 | **1.995x** | 1.995x | 4/6 | 6/6 | 6.94x |
| 1 | 0.047 | **2.011x** | 2.014x | 4/6 | 6/6 | 6.62x |
| 2 | 0.090 | **2.026x** | 2.045x | 4/6 | 6/6 | 6.32x |
| 5 | 0.198 | **2.069x** | 2.133x | 4/6 | 6/6 | 5.57x |
| 10 | 0.331 | **2.133x** | 2.340x | 4/6 | 6/6 | 4.64x |
| 20 | 0.498 | **2.237x** | 2.831x | 4/6 | 6/6 | 4.25x |

**No swept strength improves the anchor.** The best modality-aware point (S = 0) gives 1.995x against a baseline of 1.995x.

**Why it cannot help, structurally.** The sink can only push predictions DOWN, and it can only reach the three microdialysis drugs. Of those, the model currently OVER-predicts only eribulin; it already UNDER-predicts temozolomide and methotrexate. So any strength large enough to help eribulin drags the other two further below their measurements, and the geometric mean over the set gets worse before eribulin's 6.9x is dented. This is a one-directional term applied to a scorecard whose errors point both ways — the same structural mismatch `results/uptake_impact.md` §3 documented for carrier-mediated uptake, mirrored.

Per-drug at the strength nearest the Miller anchor (S = 10, modality-aware):

| drug | modality | measured | pred BEFORE | fold BEFORE | pred AFTER | fold AFTER | sink applied? |
|---|---|---|---|---|---|---|---|
| eribulin | intracerebral_microdialysis | 0.0051 | 0.0353 | 6.94x | 0.0236 | **4.64x** | **yes** |
| cefazolin | tissue_homogenate_LCMS | 0.1496 | 0.1526 | 1.02x | 0.1526 | **1.02x** | no (total-drug assay) |
| methotrexate | intracerebral_microdialysis | 0.1780 | 0.1415 | 0.79x | 0.0946 | **0.53x** | **yes** |
| temozolomide | intracerebral_microdialysis | 0.1951 | 0.1594 | 0.82x | 0.1066 | **0.55x** | **yes** |
| levetiracetam | tissue_homogenate_LCMS | 0.6972 | 0.9466 | 1.36x | 0.9466 | **1.36x** | no (total-drug assay) |
| adavosertib | tissue_homogenate_LCMS | 2.4600 | 0.5785 | 0.24x | 0.5785 | **0.24x** | no (total-drug assay) |

## 7. RECOMMENDATION — DEFAULT OFF

**Leave the immune sink default-off (it is off by construction — no shipped path imports it), and quote its spatial result as a capability rather than a calibrated prediction.** In order of weight:

1. **The magnitude is unmeasured and the literature's own sign is contested.** The best-quantified TAM drug-handling experiment in the literature (Miller 2015) found a DEPOT that made tumours respond BETTER, and the one direct free-doxorubicin measurement (Islam 2022) found TAM uptake below detection. Turning on a sink by default would assert a direction the strongest available data does not support for free drug.
2. **It does not improve the measured anchor**, and §6 explains why it structurally cannot: it is a one-directional term on a scorecard whose errors point both ways. This is the same outcome the efflux, uptake and ion-trapping runs reached, and the same conclusion follows.
3. **The spatial result is real but is SMALLER than this report originally claimed, and must be quoted at its corrected size.** The raw quintile contrast is 2.0025x, but the defensible statement is: **an excess of 1.08-1.21x over spatially-aware nulls that match the observed clumpiness (§5.3), above the entire range of both in 20 draws.** Most of the raw 2.0025 is accounted for by things that are not immune geometry: depth into tissue alone reaches 1.6394x with no immune input whatsoever (§5.1a), and simply matching the observed clumpiness reaches 1.6581x (§5.3). The correct external phrasing is *"tumour cells in macrophage-dense neighbourhoods lose measurably more drug than an equally-clumped random arrangement of the same macrophages would predict, by ~1.08-1.21x"* — not *"a 2.0x spread that a compartmental model gives as 1.0"*. It remains a genuine demonstration that geometry carries information a compartmental model cannot express; it is not a 2.0x demonstration. **And it is single-cell-only: on spot data the statistic is invalid and the generator refuses to compute it (§5.6).**
4. **The ranking is safe either way.** The cross-drug ranking is invariant across the entire swept range, so no comparative conclusion in `results/` is at risk from this term being off — or from it being turned on.

### 7.1 What measurement would settle it

One experiment closes the magnitude question: **free doxorubicin (or its fluorescent analogue) dosed into a TAM-containing tumour, then flow-sorted CD163+/TREM2+ macrophages versus tumour cells assayed for drug per cell**, giving the per-cell enrichment ratio that does not currently exist. Miller 2015 did exactly this for the CARRIER; nobody has published it for the free drug. Failing that, single-cell intravital imaging of the free drug rather than the vehicle would do. Until then this term has a literature-supported SIGN and an unmeasured SIZE, and every number in §4 must be quoted with the strength it was produced at.

## 8. Caveats — what could make this wrong

1. **The magnitude is swept, not measured.** Nothing in §4 or §5 is a prediction of how much drug real TAMs remove. Only the SHAPE of the spatial effect is independent of the magnitude.
2. **The sink is irreversible, and the best evidence says it should not be.** This module models uptake into macrophages as a pure loss. Miller 2015's central finding is the opposite: TAMs released payload back to neighbouring tumour cells over <25 um. A re-release term would partly or wholly cancel the sink and could flip its sign. **That term is not implemented**, and it is the single largest structural omission here.
3. **The immune label is a marker call, not a TAM call.** The classifier gates PTPRC/CD68/AIF1/CD163 (and CD163/TREM2 on the IO panel); it does not distinguish macrophages from lymphocytes, and only macrophages have the phagocytic/lysosomal machinery the mechanism rests on. `results/classifier_fix.md` §6 further notes the matched-null rule mildly OVER-estimates multi-marker class fractions, so `phi_imm` is likely a mild over-estimate.
4. **`r` = 7.5 um and the 25 um bandwidth are declared constants.** They set the absolute scale of `phi_imm` and trade off exactly against `k_max`. The bandwidth is anchored to a measured gradient length, the radius is not anchored to anything measured on these sections.
5. **2D.** A section is a slice; an area fraction is not a volume fraction, and a macrophage's true drug catchment is a 3D neighbourhood. This applies to every field in the project and is not specific to this term, but it bites harder here because the mechanism is explicitly about local geometry.
6. **Three sections, three patients, three panels, and no outcome label anywhere.** Nothing here says which strength is right, because nothing here can be scored against a measured drug response.
7. **The spatial structure ratio is partly confounded with depth into tissue, and one section cannot fully separate them.** §5.1a measures corr(phi_imm, edge-distance) = 0.372 on a cell-density envelope that rises 2.5x from the tissue rim inward. Stratifying gives 1.33-2.03x within strata, but with 23,812 tumour cells on ONE section there is no way to decide how much of the residual is immune geometry and how much is residual collinearity. A second single-cell section with a different immune-versus-depth arrangement is the cheapest test that would resolve it, and it has not been run.
8. **Everything in §5 is one section, one drug, one strength.** The corrected 1.08-1.21x is a single measurement with a p-value pinned at its 20-draw floor (0.0476). It is not an interval and no confidence should be attached to the second decimal place.
9. **The statistic does not transfer to spot platforms** (§5.6), so the module's flagship spatial claim is available only where single-cell centroids are measured. That restricts it to Xenium/CosMx-class data and excludes the Visium sections in `data_cache/`; the generator raises rather than reporting a number there.
10. **The scorecard test rests on assumption A1** (tissue:plasma partition = effective penetrance) and on a quasi-steady-state reading of it. Both are the scorecard's own standing assumptions, inherited unchanged.

## 9. Reproduction

```bash
./.venv/bin/python -m tmesim.validation.immune_sink_impact
```

Deterministic — the engine's seeded generator (`seed = 0`) and the nulls' own seeds (20260728, 20260729, 20260730, 20260731, 20260732) are the only randomness. Runtime is dominated by transport solves: 21 on the 424x675 flagship grid (§4, including the three engine reference runs), 22 for the demoted label-permutation null and 108 more at dx = 20 um for the honest nulls, the floors and the CV sweep (§5.1-§5.4), and 2 at dx = 12 um (§5.5). Each stage is cached to `results/immune_sink_impact_cache.pkl` as it completes, so an interrupted run resumes rather than restarting; delete that file to force a full recompute. **Every §5 number above — the nulls, the CV diagnostic, the measured floors, the edge-distance confound, the monotonicity controls and the platform controls — is produced by this command.** `write_report` refuses to write this file at all if the honest nulls are missing, so the corrected claim cannot be silently reverted by a regeneration.

**Files added by this work:**

* `tmesim/molecular/immune_sink.py` — new; the field, the swept magnitude, the literature table and the re-solve
* `tmesim/validation/immune_sink_impact.py` — new; this analysis, including the honest nulls of §5.3 and the spot-data refusal of §5.6
* `tests/test_immune_sink.py` — new; pins the bit-identity and the non-vacuity
* `tests/test_immune_sink_nulls.py` — new; pins that the primary nulls exceed the label permutation, that the CV diagnostic is emitted, and that the spot-data guard fires
* `results/immune_sink_impact.md` — this report

*Generated by `tmesim/validation/immune_sink_impact.py`. Additive only: no file under `core/`, `physics/`, `sim/`, `treatment/`, `reconstruct/` or `molecular/` (other than the new `molecular/immune_sink.py`) was modified, and the bit-exactness of the untouched transport path is verified in §1 both synthetically and end-to-end on the real section.*
