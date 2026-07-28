# Carrier-mediated uptake (OATP / OCT / OAT) — implemented, cited, and it does NOT help either

**Verdict up front: this is the mechanism both prior runs named, it was implemented with zero free parameters and one unambiguous functional form, and it fails on both measured scorecards. It is left DEFAULT-OFF. But this is not the same negative twice — unlike the efflux run, this one produces a positive DISPROOF of the specific mechanism `results/gist_imatinib.md` named, a quantitative bound showing the whole hepatic-uptake class is insufficient by the literature's own numbers, and a redirect to a candidate the chain eliminated using an argument that does not apply to it.**

| scorecard | before | after uptake | verdict |
|---|---|---|---|
| GIST tumour:normal, n = 6 pairs (geomean fold error) | **7.92x** | **7.92x** | IDENTICAL — the factor resolves to 1.0 |
| .. at the literature's UPPER BOUND on hepatic uptake (1.5x) | 7.92x | 5.65x | still 7.8x short of what the data demands |
| human glioma partition, n = 6 drugs (GM abs fold error) | **1.995x** | **2.006x** | WORSE |
| .. within 3x | 4/6 | 4/6 | |
| .. within 10x | 6/6 | 6/6 | |
| .. adavosertib specifically | 0.235x | 0.235x | unchanged — the number does not exist |

Five findings, in order of how much they are worth:

1. **The mechanism GIST named is not merely unquantified — it is falsified.** `results/gist_imatinib.md` named *"active hepatic uptake (OCT1/OATP)"*. Verified rather than assumed, the OCT1 half is contradicted by five OCT1-expressing cell lines and oocytes with a working positive control (MPP+ up 8-fold in the same oocytes), by shRNA knockdown, and — decisively — by **Oct1(-/-) mice whose HEPATIC imatinib is indistinguishable from wild-type** (PMID 24352644). Imatinib is a potent OCT1 *inhibitor* (IC50 0.095 uM), which is exactly the property that makes the inhibitor-difference assays underlying the clinical "OCT-1 activity" literature read as though it were a substrate.
2. **And the surviving hepatic-uptake candidates are quantitatively excluded, by their own literature.** The measurement demands `U_normal` = **11.7x**. The largest value any published source supports is **1.5x** — a **7.8x shortfall**. Even granting that ceiling, the geomean fold error only reaches 5.65x and 4/5 liver ratios still lie below the reachable floor. This is the same *kind* of argument that falsified the transport model, now applied to its proposed replacement.
3. **On glioma the term is blocked on the one drug that mattered, twice over.** Adavosertib is a properly controlled OATP1A2 substrate (PMID 28928160, HEK293 vs vector P = 0.001, abolished by 100 uM sulfobromophthalein P = 0.005, with OATP2B1 and OAT3 negative in the same experiment) — but the uptake fold appears **only as bars in Figure 4A**, with no value in the text, tables or abstract. And the one number that does exist, Li 2017's human GBM `Kp,uu` of 3.2, is circular against the scorecard target *and* is itself back-fitted: that paper's PBPK uptake clearance of 44 L/h was "assigned based on sensitivity analysis", i.e. tuned so the model reproduces 3.2. Passive plus efflux alone give 0.78-0.92.
4. **A widely repeated mechanistic story does not survive reading the source.** "The acidic tumour microenvironment activates OATP1A2 uptake" is asserted in `results/efflux_impact.md` on Li 2017's authority. In fact **no pH experiment was ever run on OATP1A2** in that paper — all of its pH-dependence work is on MDCKII efflux, and its own mechanism is *reduced efflux* at low pH (ER 9.2 -> 3.7 for ABCB1), not activated uptake. The PBPK model explicitly holds uptake constant. This run corrects the record.
5. **The redirect is the most useful output, and it is a correction to the chain's own reasoning.** `results/moncayo_explant.md` eliminated ion-trapping by proving a uniform retardation `R` gives `C = R*C_bath` everywhere, so it cannot produce *spatial dispersion within a section*. That proof is correct — **and it does not apply to a ratio BETWEEN two tissues**, which is what GIST measures. Meanwhile the group that dismantled the OCT1 story concluded that **lysosomal sequestration** determines intracellular imatinib, tested in GIST-T1 among other lines (PMID 26108972). Sizing it with imatinib's own pKa: the required hepatocyte lysosomal volume fraction is **3.9%**, inside the range `moncayo_explant.py` already sweeps. §11.

## 1. Nothing outside the two new files was touched

The efflux run added one default-off keyword to `treatment.drugs.delivery_source`. **This run adds nothing at all** outside `tmesim/molecular/uptake.py` and `tmesim/validation/uptake_impact.py`: no shipped code path imports `molecular.uptake`, so the term is off by construction rather than by default. Both halves are verified numerically.

- **120/120 delivery-source comparisons bit-for-bit identical** (`np.array_equal`, not `allclose`) across 4 drug names x 5 penetrances x 3 plasma concentrations x 2 branches. All identical: **True**.
- **100/100 composition comparisons bit-for-bit identical** to `reconstruct.barrier.effective_penetrance` wherever the factor is exactly 1.0 — ceiling included. All identical: **True**.
- Shipped files modified: **0**.

And the checks are not vacuous. Each uses a drug/tissue pair that actually carries a number — doxorubicin in glioma, OATP1A2, 2.0-fold, PMID 24554572 — because a non-vacuity check run on a pair whose factor is 1.0 would pass trivially:

- the factor bites where the carrier IS expressed: **True**
- the composition is genuinely tissue-differential (glioma > normal liver for the same drug): **True**
- the free-drug ceiling really does rise above 1.0: **True**

That last one is a real capability even though nothing on the scorecards gets to use it. `reconstruct.barrier.effective_penetrance` clips penetrance at 1.0 and its docstring has always said why: *"a partition above 1 requires tissue binding or active uptake, neither of which this term models."* This module models active uptake, so it raises the ceiling to the uptake factor — not a free choice but the steady-state bound `Kp,uu -> (P + CL_up)/P = U`. Adavosertib's measured 2.46 sits above the old ceiling and is unreachable without it. The machinery is in place; the literature will not supply the number to drive it.

## 2. The functional form — and why there is only one of it

Passive permeability `P` in both directions plus an active influx clearance `CL_up` in parallel gives a steady-state unbound partition

```
  f_uptake = Kp,uu = (P + CL_up)/P = 1 + b_up = UR
```

The efflux module had to choose between `1/ER` and `2/(1+ER)` — forms that differ by 2x at high ER — and justify the choice on n = 2. **That ambiguity does not arise here.** All three ways an uptake ratio is measured collapse to the same number:

| observable | value |
|---|---|
| cell accumulation fold, transfected / vector control | `1 + CL_up/P` |
| directional `Papp(A->B) / Papp(B->A)` | `1 + CL_up/P` |
| `Papp(A->B)` vs the passive-only monolayer | `1 + CL_up/P` |

Derived from the three-compartment mass balance and **verified numerically** in `molecular/uptake.py`'s self-test, not asserted. An in vivo tissue `Kp,uu` is the same quantity read in vivo, so it needs no conversion either. One form, no form-shopping — one fewer discretionary choice than the efflux term had, which is why the sensitivity table in §10 has no form axis.

Carriers acting in parallel on one membrane have additive clearances, so `UR_combined = sum_i UR_i - (n - 1)`, the exact mirror of the efflux module's parallel-pump rule.

**Zero free parameters.** Every ratio in the table is a published measurement with a PMID. Nothing is fitted to either scorecard. The one scaling factor that WAS available — converting transporter protein abundance in fmol/ug into a transport clearance — is refused for exactly the reason the efflux module refused its `SF`: the only data we could calibrate it on are the two scorecards under test.

## 3. The tissue expression map — what makes this term non-scalar

A transporter contributes in a tissue if and only if that tissue is listed here. This is the property efflux lacked: a per-drug scalar cancels out of a within-drug tissue ratio, `(f E_t)/(f E_n) = E_t/E_n`, whereas a per-(drug, tissue) factor scales it by `U_t/U_n`. Both identities are verified in the module's self-test.

| transporter | gene | expressed in | source |
|---|---|---|---|
| OATP1B1 | SLCO1B1 | normal_liver | PMID 10644574 (Abe 2000): LIVER-EXCLUSIVE, hepatocyte basolateral (sinusoidal) membrane; cloning PMID 10358072. |
| OATP1B3 | SLCO1B3 | normal_liver | PMID 10779507 (Konig 2000): liver-restricted in NORMAL tissue, hepatocyte basolateral, perivenous-zonated. |
| OCT1 | SLC22A1 | normal_liver | PMID 9187257 (Zhang 1997): HUMAN OCT1 is liver-predominant. Species caveat stated in the source itself: rodent Oct1 is also kidney, intestine and colon, so mouse data do not transfer tissue-for-tissue. |
| OATP1A2 | SLCO1A2 | normal_brain, glioma, kidney | PMID 15632119 (Lee 2005): human BRAIN CAPILLARY ENDOTHELIUM, luminal membrane, plus renal distal nephron. In glioma specifically PMID 16357150 (Bronger 2005 Cancer Res), 61 human glioma specimens: protein in the luminal membrane of the endothelial cells forming the BBB and the blood-TUMOUR barrier, but NOT in the glioma cells themselves. |
| OAT1 | SLC22A6 | kidney | PMID 11912245 (Motohashi 2002): kidney proximal tubule, basolateral membrane, cortex. |
| OAT3 | SLC22A8 | kidney | PMID 11912245 (Motohashi 2002): kidney proximal tubule, basolateral membrane, cortex. |
| OCT2 | SLC22A2 | kidney | PMID 11912245 (Motohashi 2002): kidney proximal tubule basolateral membrane, cortex and medullary ray; the same paper reports hOCT1 protein UNDETECTABLE in kidney, which is the human OCT1/OCT2 division of labour. |
| RFC | SLC19A1 | normal_liver, liver_metastasis, normal_intestine, intestinal_tumour, normal_brain, glioma, kidney | PMID 12144527 (Whetstine 2002 Biochem J): hRFC transcripts 'ubiquitously expressed' across a 76-tissue array, highest in placenta and lowest in skeletal muscle. |
| CTR1 | SLC31A1 | **(none)** | *no verifiable localisation citation — see notes* |
| SLC22A16 | SLC22A16 | **(none)** | *no verifiable localisation citation — see notes* |
| GLUT3 | SLC2A3 | normal_brain, glioma | PMID 1505464 (Nagamatsu 1992): neuron-predominant in human cerebral cortex by double labelling. Glioblastoma / brain-tumour-initiating cells PMID 23995067 (Flavahan 2013 Nat Neurosci). |

Four honest caveats, all of which cut against the result this map was built to produce:

1. **Binary is a modelling choice.** The alternative — scaling transport by measured protein abundance — needs a transport-per-fmol conversion fittable only against the datasets under test. The binary map is a literature fact; a scaled map would be a fitted parameter.
2. **`OATP1B3` is NOT tumour-absent in general.** It is liver-restricted in *normal* tissue (PMID 10779507) but abundantly expressed in gastric, colon and pancreatic cancer (PMID 11375950). It is coded as absent from the GIST liver metastasis because GIST is a mesenchymal KIT-driven tumour rather than an adenocarcinoma — an inference from lineage, since **no measurement of any uptake carrier in GIST tissue exists** (PubMed `SLCO1A2 AND GIST` returns zero; the one GIST transporter review, PMID 21163869, covers only ABCB1/ABCC1).
3. **`CTR1` and `SLC22A16` are encoded with NO tissues**, and therefore contribute nothing. CTR1 is universally described as the cisplatin uptake route but no verifiable primary human tissue-distribution citation could be retrieved; SLC22A16 carries the cleanest doxorubicin uptake Km anywhere (5.2 uM, PMID 15963465) but was cloned as testis-specific (PMID 12089149). Both are recorded and neutralised rather than silently omitted.
4. **`OATP1A2` is deliberately NOT given the intestine**, and this exclusion is not neutral. OATP1A2 is widely asserted to sit in the enterocyte apical membrane and to mediate duodenal absorption of imatinib (PMID 21508937), but no verified human intestinal localisation citation could be retrieved, and the module forbids an uncited tissue claim. Including it would move the GIST intestinal control pair — the one pair the shipped model already gets right. Sources also disagree on brain: PMID 16357150 places OATP1A2 in glioma BBB/blood-tumour-barrier endothelium and explicitly **not** in glioma cells (which is exactly where a penetrance multiplier belongs), while PMID 25132355 localises it to neurons and puts OATP2B1 in the capillaries instead.

## 4. The literature-sourced substrate table

18 drugs and 11 carriers: all ten in `molecular.library.GBM_DRUGS`, the four glioma-scorecard drugs outside it, imatinib (GIST), and the three PARP inhibitors (Moncayo) — the same 18 the efflux table covers, so the two mechanisms are scored over an identical library. Every non-`unknown` cell carries a PMID.

- confirmed substrate of at least one carrier: **8**
- .. WITH a usable published number: **4**
- .. WITHOUT one (factor forced to exactly 1.0): **4**
- .. whose only number is withheld as circular: **1**
- measured non-substrate: **3**
- unknown status throughout: **7**

**This is the same binding constraint the efflux run hit, and it is worse here.** There, 7 of 18 drugs carried a usable published ratio. Here it is 4 of 18 — and the field's reason is structural rather than accidental: uptake-transporter papers report **Km/Vmax and inhibitor IC50s**, which characterise the carrier, while a penetrance multiplier needs a **ratio against passive permeability in the same system**. Methotrexate is the extreme case: a confirmed substrate of six carriers spanning liver, kidney, BBB and every tissue, with abundant kinetics, and not one usable fold.

| drug | normal_liver | liver_metastasis | normal_intestine | intestinal_tumour | normal_brain | glioma | kidney | reason (glioma) |
|---|---|---|---|---|---|---|---|---|
| imatinib | 1.00 | 1.00 | 1.00 | 1.00 | **1.50** | **1.50** | **1.50** | applied |
| adavosertib | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | only_circular_observable |
| temozolomide | 1.00 | 1.00 | 1.00 | 1.00 | **1.70** | **1.70** | 1.00 | applied |
| lomustine | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | non_substrate |
| carmustine | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | non_substrate |
| procarbazine | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | unknown_substrate_status |
| methotrexate | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | substrate_but_no_published_ratio |
| doxorubicin | 1.00 | 1.00 | 1.00 | 1.00 | **2.00** | **2.00** | **2.00** | applied |
| cisplatin | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | no_transporter_in_this_tissue |
| etoposide | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | unknown_substrate_status |
| vincristine | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | no_transporter_in_this_tissue |
| bevacizumab | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | antibody |
| eribulin | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | unknown_substrate_status |
| cefazolin | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | no_transporter_in_this_tissue |
| levetiracetam | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | non_substrate |
| niraparib | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | unknown_substrate_status |
| rucaparib | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | unknown_substrate_status |
| olaparib | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | unknown_substrate_status |

### Citations, per assignment

**imatinib**
  - `OCT1` **no**: PMID 24352644 (Nies 2014 Clin Cancer Res): NO increase over vector control in FIVE OCT1-expressing lines (HEK, MDCK, V79, CHO) and in oocytes, with the MPP+ positive control up 8-fold in the SAME oocytes; decynium-22 and prazosin failed to inhibit; shRNA knockdown no effect; and Oct1(-/-) mice show NO change in plasma or HEPATIC imatinib. Corroborated by PMID 18483382 (Hu/Sparreboom 2008): oocytes not a substrate (P = 0.21), HEK293 only 1.20-fold, and reduced-function SLC22A1 variants had no effect on oral clearance in 74 patients (P = 0.99). Also PMID 30734690 (Blanc Mettral 2019): not a substrate of OCT1, OCT2 or OCT3.
    - *THIS ROW FALSIFIES THE MECHANISM results/gist_imatinib.md NAMED, and it is the single most consequential entry in the module. The large clinical 'OCT-1 activity' literature (PMID 16597591 White 2006, PMID 17761829 White 2007, PMID 17568400 Wang 2008, PMID 18559609 Hiwase 2008 with a prazosin-inhibitable IUR of 6.38 ng/200,000 cells) is real and clinically predictive — but every one of those studies measures an INHIBITOR-DIFFERENCE in non-transfected cells, never transport by OCT1. PMID 26108972 (Burger/Sparreboom 2015 Mol Pharmacol) shows in HEK293, K562, SD-1 and GIST-T1 that the prazosin-sensitive component is LYSOSOMAL SEQUESTRATION, not SLC22A1. And imatinib is a potent OCT1 INHIBITOR (IC50 0.095 uM, PMID 24352644), which is exactly the property that makes inhibitor-difference assays read as if it were a substrate. Inhibitor is not substrate; that distinction is the whole difference between this row being 'no' and being 'yes'.*
  - `OATP1B1` **unknown**
    - *GENUINELY DISPUTED between two direct transfected-cell studies, so coded unknown rather than resolved in either direction. FOR: PMID 23340295 (Zimmerman 2013 Clin Cancer Res), HEK293 stably transfected, 0.1 uM, 30 min — imatinib among the TKIs with statistically significant uptake over vector in BOTH OATP1B1 and OATP1B3. AGAINST: PMID 18483382 (Hu 2008), oocytes at 0.2 uM — OATP1B1 negative while OATP1B3 was positive.*
  - `OATP1B3` **yes**: PMID 23340295 (Zimmerman 2013 Clin Cancer Res): HEK293 stably transfected with OATP1B3, [3H]imatinib 0.1 uM, 30 min, statistically significant uptake over vector control; the same paper states imatinib was 'previously identified as a substrate of OATP1B3'. Concordant with PMID 18483382 (Hu 2008), oocytes 0.2 uM.
    - *CONFIRMED SUBSTRATE, NO USABLE NUMBER — the exact situation doxorubicin was in for the efflux table. Zimmerman publishes imatinib's uptake only as bars in Figure 1A with asterisks; no fold appears in the text, the tables or the figure caption, and no Km/Vmax exists for imatinib on any OATP. Reading a value off the figure would be a guess, so the factor is 1.0. What the figure DOES license is an upper BOUND, and the report uses it as such rather than as a value: sorafenib was the highest TKI in that panel at ~1.5-fold, so imatinib is necessarily below 1.5-fold. A bound that falsifies is legitimate where a bound used as a value would not be.*
  - `OATP1A2` **yes**: PMID 21508937 (Eechoute/Sparreboom 2011 Clin Pharmacol Ther): Xenopus oocytes, [3H]imatinib, uptake over water-injected control; rosuvastatin-inhibitable at both pH values tested. Independently PMID 21633340 (Yamakawa 2011 Clin Pharmacol Ther): significantly enhanced uptake in OATP1A2-transfected HEK293 (P = 0.002), naringin-inhibitable in HEK293, Caco-2 and K562.
    - in vitro uptake ratio **1.50** — Xenopus oocyte uptake of [3H]imatinib vs water-injected control at pH 7.4 (1.5-fold, P = 0.0008); rises to 2.4-fold at pH 5.0 (P = 0.0001), i.e. the carrier is pH-sensitive and the physiological-pH value is recorded here as the conservative one — PMID 21508937 (Eechoute 2011 Clin Pharmacol Ther)
    - *THE ONLY NUMERIC UPTAKE RATIO ANYWHERE FOR IMATINIB — and it is for a transporter that is NOT hepatic, so it cannot rescue the GIST prediction. OATP1A2 sits in brain capillary endothelium and renal distal nephron (PMID 15632119), not in hepatocytes. Two further warnings travel with it. (1) THE CLINICAL TRANSLATION FAILED, in the source paper itself: no association with OATP1A2 variants in 94 patients, and rosuvastatin had NO effect on imatinib PK in 12 GIST patients. The authors' own conclusion is that 'this transporter by itself is unlikely to contribute substantially to the absorption profiles of imatinib in humans'. (2) The 1.5-fold is an oocyte number at 0.2-1 uM.*
  *Notes:* THE ANCHOR CASE, AND THE MECHANISM IT WAS SUPPOSED TO SUPPLY IS FALSIFIED. results/gist_imatinib.md named 'active hepatic uptake (OCT1/OATP)' as the missing physics. Verified rather than assumed, the OCT1 half is contradicted by five transfected cell lines, oocytes with a working positive control, shRNA knockdown and Oct1(-/-) mice whose HEPATIC imatinib is indistinguishable from wild-type (PMID 24352644). The OATP half survives as substrate status but carries no publishable number and is bounded below 1.5-fold, against a required ~11.7-fold. Imatinib's real intracellular-accumulation mechanism, per the group that dismantled the OCT1 story, is LYSOSOMAL SEQUESTRATION (PMID 26108972) — tested in GIST-T1 among other lines. That is the redirection this run produces and it is developed in results/uptake_impact.md. Note also that no measurement of ANY uptake carrier in GIST tumour tissue exists: PubMed 'SLCO1A2 AND GIST' returns zero, and the one GIST transporter review (PMID 21163869) covers only ABCB1/ABCC1. Human liver:plasma and human tumour Kp,uu for imatinib are both UNKNOWN; the available figures are mouse (liver:plasma AUC 2.29, PMID 21635257), and Oct1(-/-) mice show that number is not OCT1-dependent.

**adavosertib**
  - `OATP1A2` **yes**: PMID 28928160 (Li 2017 Clin Cancer Res): cellular uptake of AZD1775 in HEK293 overexpressing OATP1A2 significantly higher than vector control (P = 0.001), reduced to control level by 100 uM sulfobromophthalein (P = 0.005), with an estrone-3-sulfate positive control (P < 0.001). OATP2B1 and OAT3 negative in the same experiment.
    - in vivo Kp,uu **3.20** (human glioblastoma tumour) — unbound tumour-to-plasma concentration ratio, median 3.2 (range 1.3-24.4) in n = 20 patients, phase 0 single oral dose with resection at 2-26 h; fu MEASURED in the patients' own specimens (fu,brain 0.067 +/- 0.013, fu,plasma 0.182 +/- 0.024), Kp,total median 9.1 — PMID 28928160 (Li 2017 Clin Cancer Res) — **WITHHELD as circular against glioma_regional_pk**
    - *A PROPERLY CONTROLLED POSITIVE WITH NO PUBLISHABLE NUMBER. The uptake fold appears ONLY as bars in Figure 4A; no value is given in the text, the tables or the abstract, and no Km/Vmax is reported. PubMed 'AZD1775 AND OATP1A2' returns exactly one record, this one, and it publishes only P-values. Under this module's rule a figure is not a number, so uptake_ratio stays None — the same treatment doxorubicin got in the efflux table. THE Kp,uu IS WITHHELD AS CIRCULAR: it is a human glioblastoma tumour-to-plasma ratio, i.e. the same quantity in the same tissue as the scorecard target. Worse, it is not even an independent measurement of the mechanism: Li 2017's PBPK uptake clearance of 44 L/h was 'assigned based on sensitivity analysis', i.e. TUNED so that the model reproduces 3.2 (passive plus efflux alone give Kp,uu 0.78-0.92). And no pH experiment was ever run on OATP1A2 — all the pH work in that paper is on MDCKII efflux — so the widely repeated 'acidic tumour microenvironment activates OATP1A2 uptake' story is a PBPK inference, not a measurement.*
  *Notes:* THE PRE-REGISTERED TARGET OF THIS RUN, AND IT CANNOT BE SCORED. The efflux run identified adavosertib as the drug where efflux was most clearly the wrong mechanism and named uptake as the right one. Uptake IS the right mechanism — the substrate status is demonstrated with inhibitor and positive controls — but the model needs a number and no number exists that is both published and non-circular. Two independent human cohorts agree the drug accumulates (Li 2017 Kp,uu median 3.2, n = 20; Lee 2025 PMID 39820282 intratumoral:plasma 4.18 +/- 3.36 contrast-enhancing and 0.74 +/- 0.63 non-enhancing), and note that the non-enhancing value BELOW 1.0 means the accumulation is a blood-TUMOUR-barrier phenomenon in leaky core, not a BBB phenomenon. Contradicting mouse data are on record and unresolved: brain:plasma 1.0 in wild-type mice (PMID 29147815) and normal-brain:blood 0.04-0.05 (PMID 25609063), a 9-25x interspecies gap neither human paper explains. Placement note: PMID 16357150 puts OATP1A2 in the glioma BBB/blood-tumour BARRIER endothelium and explicitly NOT in glioma cells, which is exactly where a penetrance multiplier acts — so the composition used here is supported rather than assumed.

**temozolomide**
  - `GLUT3` **yes**: PMID 40813374 (Diao 2025 Cell Death Discov): GLUT3 knockdown in U343-MG reduced TMZ uptake by 54.3%; stable GLUT3 overexpression in U251-MG increased uptake ~70% vs empty-vector control, LC-MS/MS with a 2H3-TMZ internal standard; site-directed mutagenesis of 8 residues; mouse xenograft confirmation.
    - in vitro uptake ratio **1.70** — cellular TMZ uptake, stable GLUT3-overexpressing U251-MG vs stable empty-vector control (+70%), LC-MS/MS with 2H3-TMZ internal standard — the canonical transfected-versus-vector form of the observable — PMID 40813374 (Diao 2025 Cell Death Discov)
    - *THE WEAKEST NUMERIC ROW IN THE TABLE, INCLUDED WITH ITS PROBLEMS ON DISPLAY, and leave-one-out tested in the report. (1) Single 2025 paper, no independent replication found. (2) The companion affinity data are internally suspect: the paper reports GLUT3 glucose Kd 1-10 uM against a well-established glucose Km of ~1.4 mM, a ~1000x discrepancy, so the MST numbers should not be trusted — the uptake numbers are a different measurement and are the only thing used here. (3) Kd is binding, not transport; no Km/Vmax for TMZ transport exists anywhere. (4) The knockdown arm implies a larger endogenous ratio of 1/(1-0.543) = 2.19; the smaller overexpression figure is used because it matches the canonical observable and is the conservative one. (5) The paper's own framing claim that TMZ uptake 'is regulated by Wnt signalling' could not be traced to a primary source and is not relied on.*
  *Notes:* Corrects a common assumption in the other direction from the efflux table. TMZ is routinely described as entering cells purely by passive diffusion; a 2025 report contradicts that with a carrier. Everything else is UNKNOWN rather than negative: no OATP, OCT, OAT, MCT, ENT, PEPT or OCTN assay for TMZ exists. LAT1 is coded nowhere because the only LAT1 result (PMID 38052257) concerns an ENGINEERED amino-acid prodrug, not temozolomide. The active species MTIC is generated intracellularly by spontaneous, pH-driven ring opening (PMID 8049205) rather than by a localised enzyme, and no MTIC transport study exists. The human CSF anchor, AUC_CSF/AUC_plasma = 20% in n = 35 (PMID 15173079), is total drug and is not attributable to any transporter, so it is NOT used here.

**lomustine**
  - `OATP1B1` **no**: PMID 557367 (Begleiter 1977 Cancer Res): 'uptake of CCNU also occurs by passive diffusion' — temperature-independent, Na-independent, unaffected by NaF/NaCN/2,4-DNP, and NOT self-inhibited by excess unlabelled CCNU or BCNU, i.e. non-saturable. A non-saturable, energy-independent uptake excludes carrier mediation generically.
    - *Recorded against OATP1B1 as the representative carrier because the source excludes carrier mediation as a class rather than testing named transporters. TRAP IN THIS PAPER: total RADIOACTIVITY uptake for ring-labelled CCNU WAS temperature-dependent, which the authors attribute to drug decomposition preceding uptake. Anyone citing 'temperature-dependent CCNU uptake' has misread it.*
  *Notes:* An improvement on the efflux table, where lomustine was 'unknown' in both cells: here there is a genuine, well-designed MEASURED NEGATIVE from 1977 that excludes carrier-mediated uptake as a class. Corroborated by PMID 6227421 (intracellular CCNU constant at 85 pmol/uL in L1210 over 10-60 min). Discard the DrugBank-derived claim that lomustine is a major CYP2D6 substrate: it is unverifiable and contradicts the GLEOSTINE label, which states the elimination pathways have not been characterised.

**carmustine**
  - `OATP1B1` **no**: PMID 557367 (Begleiter 1977 Cancer Res): 'uptake of intact BCNU was by passive diffusion'; cell/medium distribution ratio 0.2-0.6, i.e. BELOW 1 and therefore not concentrative, and no saturation. Replicated by PMID 1692254: steady-state cell:medium 0.52 +/- 0.12 in SK-MG-1 human glioma.
    - *As for lomustine, recorded against a representative carrier because the source excludes carrier mediation as a class. A cell:medium ratio below 1 is a particularly clean negative for an UPTAKE term: a concentrative carrier cannot produce it.*
  *Notes:* CRITICAL MIS-ATTRIBUTION RISK, AVOIDED AND RECORDED. PMID 8341272 (Noe 1993 Mol Pharmacol) reports the only Km/Vmax in the entire nitrosourea literature — Km 2.39 +/- 0.37 mM, Vmax 236 +/- 53 pmol/uL/s, epinephrine Ki 163 +/- 15 uM, via the extraneuronal monoamine transporter (EMT/OCT3, SLC22A3) — but it is for SarCNU, an EXPERIMENTAL ANALOGUE, not carmustine. The same paper explicitly contrasts SarCNU with BCNU, 'which enters cells primarily by passive diffusion'. Importing that Km/Vmax would attribute a different molecule's carrier to this one. (Even within the SarCNU series it is contested: PMID 11181897 found no correlation between EMT expression and SarCNU efficacy across nine xenografts.) Also do not miscite PMID 6631721: MeCCNU/semustine impairing renal PAH transport in rats is transporter INJURY, not substrate status.

**procarbazine**
  - no carrier assignment: see notes.
  *Notes:* GENUINELY UNKNOWN for every carrier, and honestly so. The single on-topic primary source is PMID 708471 (Lam 1978 Biochem Pharmacol), 'On the mechanism of uptake of procarbazine by L5178Y lymphoblasts in vitro' — PubMed carries no abstract and the full text could not be retrieved, so its conclusion cannot be stated. The paper exists and is exactly on topic; that is all that can be said. Do NOT infer 'passive by analogy' with the nitrosoureas: the same group published carrier-mediated transport for melphalan in the same cell line (PMID 837375). Procarbazine is small, lipophilic and equilibrates rapidly with CSF, which is consistent with passive diffusion — but that is physicochemical REASONING, not a transport measurement, and it does not exclude a carrier. Modelling caveat: procarbazine is a prodrug requiring CYP/MAO activation (PMID 2908840, PMID 1316811) and transport of the active azoxy metabolites has never been studied at all. Its MAO INHIBITION is enzyme pharmacology and must not be read as transporter substrate status.

**methotrexate**
  - `OATP1A2` **yes**: PMID 16702441 (Badagnani 2006 JPET): saturable OATP1A2-mediated MTX transport in Xenopus oocytes, Km 457 +/- 118 uM, Vmax 17.5 +/- 4.9 pmol/oocyte/60 min; acidic pH stimulated MTX uptake up to 7-fold.
    - *A MEASURED Km IS NOT AN UPTAKE RATIO, and the distinction is load-bearing. This module's factor is `UR = 1 + CL_up/P`, which needs the transporter clearance AND the passive permeability in the same system. Km/Vmax alone gives neither, and no transfected-versus-vector fold is published for MTX on OATP1A2. The 'up to 7-fold' figure is a pH RATIO within the transfected cells, not a transfected-over-control ratio, and using it would be a category error. (The 30-fold value in that abstract is for the estrone-3-sulfate probe, not MTX — a trap.) Factor stays 1.0.*
  - `OATP1B1` **yes**: PMID 26474710 (Kim 2015 Mol Pharm): MTX uptake in HEK293-OATP1B1, characterised via inhibitor IC50s (rifampicin 0.9 uM at 1B1, 0.3 uM at 1B3; telmisartan 6.7 / 7.9 uM).
    - *Substrate status is implied by the assay design rather than by a reported Km or fold; no usable ratio. Population-PK support exists but is NOT transport evidence: SLCO1B1 variants explain 10.7% of MTX clearance variance in childhood ALL (PMID 22147369).*
  - `OATP1B3` **yes**: PMID 11375950 (Abe 2001 Gastroenterology): saturable, dose-dependent MTX transport in Xenopus oocytes; OATP1B3/LST-2 transfection increased MTX sensitivity. No Km reported.
  - `OAT1` **yes**: PMID 12130730 (Takeda 2002 JPET): hOAT1-mediated MTX transport, Km 553.8 uM.
  - `OAT3` **yes**: PMID 12130730 (Takeda 2002 JPET): hOAT3-mediated MTX transport, Km 21.1 uM (hOAT4 17.8 uM).
  - `RFC` **yes**: PMID 36071163 (Dang 2022 Nature): cryo-EM structure of human RFC WITH METHOTREXATE BOUND; 'major importer of folates ... as well as chemotherapeutic agents such as methotrexate'.
    - *No Km could be verified from any primary source read for this table — the commonly quoted low-micromolar value did not appear in any retrievable paper, so it is marked unverified rather than hard-coded. In any case RFC is UBIQUITOUS and therefore cancels out of any tumour:normal ratio (see the RFC transporter note).*
  *Notes:* THE BEST-COVERED DRUG IN THE TABLE AND STILL NO USABLE FACTOR. Methotrexate is a confirmed substrate of six distinct carriers spanning liver, kidney, BBB and every tissue — and not one of them yields the quantity this module needs. Km/Vmax and inhibitor IC50s are abundant; transfected-over-control uptake FOLDS are absent. That is a precise statement of the coverage gap: the field measures affinity and inhibition, while a penetrance multiplier needs a ratio against passive permeability in the same system. PCFT/SLC46A1 is left out entirely: 27 PubMed records link MTX and PCFT but no primary transport paper could be verified, and an unverified carrier is not encoded. Also do not cite PMID 20448812 for OATP1B1/1B3 MTX transport — its probe is fluorescein-methotrexate, a conjugate, not MTX.

**doxorubicin**
  - `OATP1A2` **yes**: PMID 24554572 (Durmus 2014 Int J Cancer): doxorubicin uptake in OATP1A2-overexpressing HEK293 vs control. Independently confirmed in a HeLa vaccinia single-transporter system, PMID 27777271 (Lee 2016 Mol Pharmacol), where four SLCO1A2 variants showed impaired transport.
    - in vitro uptake ratio **2.00** — ~2-fold cellular uptake in OATP1A2-transfected HEK293 over vector control — the canonical transfected-versus-control observable, in a single-transporter background — PMID 24554572 (Durmus 2014 Int J Cancer)
    - *A QUIETLY IMPORTANT ROW. Doxorubicin was the anchor case the efflux table could NOT put a number on — the textbook P-gp substrate and the molecule this repo's perivascular constant alpha was fitted to, left at exactly 1.0 because no defensible efflux ratio existed. Here it is one of the few drugs that DOES carry a usable published uptake ratio. Doxorubicin is not on either scorecard, so this changes no score; it is recorded because it shows the coverage gap is drug-specific and mechanism-specific rather than uniform.*
  - `OATP1B1` **unknown**
    - *GENUINELY ASSAY-FORMAT DEPENDENT, so coded unknown. NEGATIVE in vitro: both PMID 24554572 and PMID 27777271 found no doxorubicin uptake by OATP1B1 in HEK293 / HeLa single-transporter systems. POSITIVE in the SAME paper (PMID 27777271) in MDCKII TRANSWELL assays: 'OATP1B1, OATP1B3, and OATP1A2, either alone or in combination with MDR1, significantly transported doxorubicin'. POSITIVE in vivo: liver-specific transgenic OATP1B1 partially rescued hepatic uptake in Oatp1a/1b-/- mice (PMID 24554572), where Oatp1a/1b deficiency lowered the liver-to-plasma ratio up to 4.1-fold.*
  - `OATP1B3` **unknown**
  - `SLC22A16` **yes**: PMID 15963465 (Okabe 2005 BBRC): saturable, dose-dependent doxorubicin transport by SLC22A16 (OCT6/CT2) in Xenopus oocytes, apparent Km 5.2 +/- 0.4 uM; Jurkat transfectants significantly more doxorubicin-sensitive at 2 uM.
    - *Contributes nothing because SLC22A16 is encoded with no tissues (testis-specific, PMID 12089149) — see that transporter's note. The Km is real and is the cleanest doxorubicin uptake kinetic anywhere, which is exactly why it is recorded and neutralised rather than silently omitted.*
  - `OCT1` **unknown**
    - *PubMed returns ZERO records for doxorubicin uptake by SLC22A1 or SLC22A3. This is a true null, not a negative.*
  *Notes:* The mirror image of its efflux row. Note the in vivo mouse numbers are attributable to Oatp1a/1b as a family, not to any single human orthologue: Oatp1a/1b-/- gave plasma AUC 1.3-fold higher and liver AUC 1.4-fold lower with a liver-to-plasma ratio up to 4.1-fold lower (PMID 24554572), and Slco1a/1b-/- gave plasma AUC 1.7-fold higher with a 2.3-fold lower liver-to-plasma ratio and 40% reduced clearance (PMID 27777271). Those are mouse-family numbers and are NOT used as a human single-transporter uptake ratio.

**cisplatin**
  - `CTR1` **yes**: PMID 12370430 (Ishida 2002 PNAS): yeast CTR1 deletion confers resistance with reduced intracellular cisplatin; mouse cells lacking one or two mCtr1 alleles show resistance and decreased accumulation in parallel with gene dosage.
    - *ACCUMULATION AND RESISTANCE, NOT A DIRECT FLUX MEASUREMENT, and no Km/Vmax exists. Combined with the unverified tissue distribution (see the CTR1 transporter note) this row contributes exactly nothing, which is the honest outcome for a mechanism everybody cites and nobody has localised in a retrievable human paper.*
  - `OCT2` **yes**: PMID 16914559 (Yonezawa 2006 JPET): HEK293 transient transfection — both cytotoxicity and CELLULAR ACCUMULATION of cisplatin enhanced by hOCT2; OCTN1 and OCTN2 explicitly negative, and carboplatin and nedaplatin were not substrates of any of the seven carriers tested (a useful internal control).
    - *Renal, so it contributes nothing in brain or liver. In vivo corroboration: OCT1/2 double-knockout mice show no ototoxicity and only mild nephrotoxicity, and cimetidine protects wild-type mice (PMID 20110413). Pharmacogenetics (association only, not transport): SLC22A2 rs316019 protective against ototoxicity.*
  - `OCT1` **unknown**
    - *CONTRADICTORY and left unresolved. PMID 16914559 reports accumulation 'weakly' enhanced by hOCT1, i.e. a weak substrate; PMID 16314463 (Ciarimboli 2005) found cisplatin at 100 uM inhibited transport via hOCT2-HEK293 but NOT hOCT1-HEK293, and no effect in human hepatocyte couplets. Note also that Ciarimboli measured inhibition of a fluorescent cation probe plus apoptosis rather than cisplatin flux, so its claims are inferential in both directions.*
  *Notes:* Expectation verified rather than assumed, as in the efflux table, and with the same conclusion from the other side: platinum handling runs through copper transporters and glutathione, not through the OATP/OCT carriers this module models. The net effect on both scorecards is exactly 1.0, and cisplatin is on neither.

**etoposide**
  - `OATP1B1` **unknown**
    - *DISPUTED WITHIN A SINGLE PAPER, AND THE ABSTRACT IS THE PART THAT IS WRONG. PMID 26663398 (Drenberg 2016 Clin Pharmacol Ther) abstract asserts 'cytarabine, daunorubicin, etoposide, and mitoxantrone are substrates for OATP1B1 and the mouse ortholog Oatp1b2', but its own Results state verbatim that intracellular uptake was significantly enhanced over vector control 'with the exception of VP16' — VP16 being etoposide. Read the number, not the claim. The paper's Discussion then cites OTHERS for etoposide's OATP1B1 substrate status, not its own null arm.*
  - `OATP1B3` **unknown**
    - *Tested only as an INHIBITOR: Ki 13.5 uM against OATP1B3-mediated [3H]CCK-8 (PMID 25618019). No transport of etoposide was demonstrated, so this must not be upgraded to substrate.*
  *Notes:* A SPECIES SPLIT, recorded because it is exactly the kind of thing that gets flattened into a wrong 'yes'. In PMID 26663398 the MOUSE ortholog Oatp1b2 did significantly transport etoposide while human OATP1B1 did not reach significance in the same experiment; in vivo, Oatp1b2-/- mice showed 1.2-2.4-fold reduced liver accumulation and 1.3-1.5-fold higher plasma AUC. Monkey OATP1B1/1B3/2B1 also transport it (PMID 30652318). The human single-transporter evidence is the null arm. Etoposide is on neither scorecard.

**vincristine**
  - `OATP1B1` **yes**: PMID 28828541 (Lancaster 2017 Pharm Res): OATP1B1-transfected CHO cells showed significantly increased vincristine uptake versus wild-type, reduced by OATP inhibitors.
    - *REAL BUT QUANTITATIVELY MINOR, and the source itself says so: in suspended human hepatocytes 'only 10% of the total vincristine uptake process could be attributed to OATP-mediated transport'. So even if a fold existed, the correct factor would be near 1.1, not near the value a transfected CHO line would suggest. No fold, no Km, no Vmax is published; single 2017 paper; heterologous non-hepatic background.*
  - `OATP1B3` **yes**: PMID 28828541 (Lancaster 2017 Pharm Res): same experiment, OATP1B3-transfected CHO cells.
    - *Also an INHIBITOR of OATP1B3, Ki 30.6 uM (PMID 25618019) — not substrate evidence. And beware the vinca substitution trap that the efflux table also flags: much of the canonical literature is VINBLASTINE (Ki 10.2 uM at OATP1B1, 15.9 uM at OATP1B3), which is not interchangeable with vincristine for a quantitative term.*
  *Notes:* Hepatic disposition of vincristine is dominated by a benzbromarone-sensitive EFFLUX mechanism (possibly MRP3) and by MDR1/P-gp, not by uptake. 'Vincristine is an OATP1B1 substrate' is true and should not carry a large uptake term in any PBPK model. Factor 1.0 for want of a number, which here happens to also be roughly the right answer.

**bevacizumab**
  - no carrier assignment: see notes.
  *Notes:* NOT APPLICABLE rather than measured-negative, and the module short-circuits antibodies to 1.0 before this row is read. Solute carriers translocate small molecules through a substrate-binding pore; a ~149 kDa IgG1 cannot be an SLC substrate. That is mechanistic REASONING, not an experimental negative — nobody has run a bevacizumab OATP1B1 assay and nobody should. The relevant transport route is FcRn (FCGRT), a RECEPTOR and not an SLC: IgG catabolism is 10x faster in FcRn-deficient mice (PMID 8643606), it mediates bidirectional transcytosis across human epithelium (PMID 10510331), and it is present in endothelium, macrophages, enterocytes and choroid plexus epithelium (PMID 28402755, PMID 12067234). IMPORTANT NUANCE THAT CUTS AGAINST THE USUAL NARRATIVE: PMID 19636712 measured a blood-corrected brain:plasma AUC ratio of 0.0022 for IgG1 in wild-type mice and 0.0021 in FcRn-KO (P = 0.33), concluding FcRn is NOT responsible for poor brain exposure — so FcRn should not be cited as the cause of an antibody's BBB exclusion. This repo already handles bevacizumab in a size/permeability term (bbb_penetrance = 0 via the is_antibody and MW gates), which is the right place for it.

**eribulin**
  - `OATP1B1` **unknown**
  - `OATP1B3` **unknown**
  - `OCT1` **unknown**
  - `OAT1` **unknown**
  - `OAT3` **unknown**
  *Notes:* CODED UNKNOWN, NOT 'NO', AND THE PRECEDENT IS THIS REPO'S OWN. The HALAVEN US Prescribing Information section 12.3 states that eribulin at clinically relevant concentrations 'is not a substrate of ... OATP1B1, OATP1B3 ... OCT1, OCT2 ... OAT1, OAT3 ... MATE1 ... or BCRP', and that it IS a substrate of P-gp. That is a label claim with NO PMID and an unpublished underlying study — exactly the situation in which results/efflux_impact.md coded eribulin's BCRP status 'unknown' rather than 'no'. The same standard is applied here rather than a more convenient one. Note also the label's qualifier 'at clinically relevant concentrations', which leaves higher concentrations untested. OATP1A2, LAT1, MCT, ENT, PEPT and OCTN have never been tested. Either way the factor is 1.0, so eribulin's 6.94x over-prediction — the largest single error on the glioma scorecard — is untouched by this mechanism.

**cefazolin**
  - `OAT3` **yes**: PMID 14984259 (Sakurai 2004 Pharm Res): 'the uptake study using HEK293 transfectants revealed that cefazolin and phenolsulfonphthalein were transported by hOAT3'. Corroborated by PMID 16098483 (Ueo 2005 Biochem Pharmacol): accumulation significantly greater in HEK-hOAT3 than control.
    - *No cefazolin-specific fold could be verified (paywalled); companion cephalosporins in the same panel ran 35-50x over control, which is suggestive but is a DIFFERENT MOLECULE'S number and is not imported. No Km for cefazolin at hOAT3 exists. Renal anyway, so it contributes nothing in brain.*
  - `OAT1` **no**: PMID 16098483 (Ueo 2005 Biochem Pharmacol): uptake of these antibiotics by HEK-hOAT1 was 'within two-fold of that by control cells'; the open-access review PMC3153278 (PMID 22069672) reads the primary data and states hOAT1-dependent uptake was non-significant for cefazolin and that 'cefazolin is also not an hOAT1 substrate'.
    - *THIS CORRECTS results/efflux_impact.md, WHICH RECORDED CEFAZOLIN AS 'a confirmed OAT1 substrate (PMID 15618660, Ki 72 uM vs PAH)'. Two errors are compounded in that note and both are fixed here. (1) PMID 15618660 is RAT Oat1, and the human orthologue does not transport cefazolin — a species split that a 2025 review (PMID 40563504) already mis-propagates into a human OAT1 substrate table. (2) A Ki is INHIBITION of PAH transport, not substrate status. The rat result stands on its own terms (rOAT1 significantly stimulated cefazolin uptake, IC50 72 +/- 8 uM vs PAH) and is recorded as rat.*
  - `OATP1B1` **yes**: PMID 18974612 (Nakakariya 2008 Drug Metab Pharmacokinet): OATP1B1-mediated cefazolin transport, Km 20.8 +/- 3.6 mM, Vmax 31.7 pmol/hr/oocyte.
    - *MILLIMOLAR, i.e. two orders of magnitude above clinical cefazolin plasma concentrations (~100-200 uM total), and the value is EXTRAPOLATED BEYOND THE TESTED RANGE (0.1-10 mM). Recorded for completeness; almost certainly not clinically meaningful. Factor 1.0 regardless, since a Km is not an uptake ratio.*
  - `OATP1B3` **yes**: PMID 18974612 (Nakakariya 2008 Drug Metab Pharmacokinet): OATP1B3-mediated cefazolin transport, Km 3.89 +/- 2.21 mM, Vmax 10.5 pmol/hr/oocyte.
    - *Millimolar, as OATP1B1. Same caveat.*
  *Notes:* PEPT1 is deliberately coded nowhere and the reason is a good one: cefazolin lacks the alpha-amino group that PEPT1 substrates carry, and this is verified rather than assumed — GlySar uptake was 'unaffected by cefazolin' (PMID 22490229), with no intracellular pH drop in Caco-2 (PMID 15618692) and a negative rabbit in vivo result (PMID 11765153), while cefadroxil and cephradine (which DO have the alpha-amino group) are substrates. PEPT2 has never been directly tested. HUMAN IN VIVO CAVEAT: probenecid prolongs cefazolin serum half-life by 20-40% (PMID 4758839, n = 12, replicated), but probenecid is non-selective across OAT1, OAT3, MRP2 and MRP4, so the interaction does not localise to OAT3 — and PMID 4758839 itself concludes cefazolin is 'cleared primarily by the glomerulus', a flagged contradiction with the modern OAT3-secretion picture. An expression-phenotype correlation exists (renal hOAT3 mRNA vs cefazolin elimination rate constant, r = 0.536, PMID 14984259) but no genotype-PK association does. Everything relevant here is renal; in glioma cefazolin's factor is 1.0.

**levetiracetam**
  - `OCT1` **no**: PMID 35955563 (Sachkova 2022 Int J Mol Sci): OCT1-HEK293 vs mock uptake ratio = 1.00 against the authors' own substrate threshold of > 1.5.
    - *A clean measured negative in a properly controlled screen.*
  *Notes:* THE ORPHAN-CARRIER HYPOTHESIS WAS TESTED HEAD-ON AND FAILED, which is a stronger result than 'unknown'. The proton-coupled organic-cation antiporter is the long-postulated unidentified influx carrier at the BBB; PMID 35955563 screened 72 drugs in hCMEC/D3 and levetiracetam gave antiport ratios of 1.17 +/- 0.28 and 1.19 +/- 0.24 against the authors' criterion that ratios near 1.0 indicate NON-substrates. No saturable brain uptake exists in any species: serum and CSF rise linearly with dose in rat (PMID 10210031) and brain-ECF entry is linear with no regional heterogeneity (PMID 11454660). CORRECTION TO A COMMON PREMISE, including one this repo could easily have inherited: levetiracetam's passive permeability is MODEST, not high — PS = 0.015 mL/min/g, 21x lower than brivaracetam (PMID 26663401) — so its slow brain entry is low passive permeability, not a missing carrier. Two traps avoided: the only Michaelis-Menten constant in the levetiracetam literature (Km 435 uM, PMID 14530892) belongs to a blood ESTERASE, not a transporter; and the metabolite ucb L057 is actively secreted and probenecid-sensitive while the parent is not, so 'levetiracetam is actively secreted' conflates the two. Every other carrier is UNKNOWN (never tested), which is a genuine null rather than a negative.

**niraparib**
  - `OATP1B1` **no**: PMID 34730366 (Martins/Loos/Schinkel 2021 Mol Pharm): Oatp1a/1b deletion did not significantly affect niraparib oral bioavailability or LIVER DISTRIBUTION in vivo (it did reduce liver uptake of the M1 metabolite). Concordant with the ZEJULA US PI section 12.3 and the EMA EPAR (EMA/648982/2017), both of which state niraparib is not a substrate of OATP1B1, OATP1B3 or OCT1 — label claims with no PMID, cited here only as corroboration of the peer-reviewed in vivo negative.
    - *An in vivo hepatic negative in the transporter-knockout design is a stronger negative than a monolayer screen, and it is peer-reviewed, which is why this row is 'no' where eribulin's is 'unknown'.*
  - `OATP1B3` **unknown**
  - `OCT1` **unknown**
  *Notes:* Distinguish substrate from inhibitor, which matters for this drug: niraparib is a potent MATE1/MATE2-K INHIBITOR (IC50 0.18 uM / <= 0.14 uM, ZEJULA label) and that is not transport. UNRESOLVED CONTRADICTION recorded rather than smoothed: the FDA label says niraparib 'did not inhibit ... OCT1' while the EMA EPAR says it 'weakly inhibits OCT1 with an IC50 = 34.4 uM' and advises caution with metformin. The measured 6-7x brain:plasma rise in Abcb1a/1b-KO mice (PMID 34730366) is LOSS OF EFFLUX, not gain of uptake, and belongs to the efflux module. No SLCO1B1/SLC22A1 genotype-exposure study exists for any PARP inhibitor.

**rucaparib**
  - `OATP1B1` **no**: PMID 32129697 (Liao 2020 Xenobiotica): 'rucaparib was a substrate for P-gp and BCRP, but not for OATP1B1, OATP1B3, OAT1, OAT3, or OCT2.' Peer-reviewed, unlike the label-only negatives elsewhere in this table.
  - `OATP1B3` **no**: PMID 32129697 (Liao 2020 Xenobiotica): same experiment.
  - `OAT1` **no**: PMID 32129697 (Liao 2020 Xenobiotica): same experiment.
  - `OAT3` **no**: PMID 32129697 (Liao 2020 Xenobiotica): same experiment.
  - `OCT1` **unknown**
    - *TESTED ONLY AS AN INHIBITOR, NEVER AS A SUBSTRATE. Rucaparib inhibits OCT1 with IC50 4.3 uM (PMID 32129697), which is not evidence either way about transport.*
  *Notes:* THE DRUG MOST LIKELY TO BE MISCODED IN THIS WHOLE TABLE, because it is a potent inhibitor of the very carriers it is not a substrate of: MATE2-K IC50 0.19 uM, MATE1 0.63 uM, OCT1 4.3 uM, OCT2 31 uM (PMID 32129697). UNRESOLVED REGULATORY CONTRADICTION recorded: the RUBRACA US PI states flatly that 'rucaparib inhibited OATP1B1, OATP1B3, OAT1, OAT3 ...' while the EMA EPAR states that at clinical exposures it did NOT inhibit them; PMID 32129697 explains the divergence (only 42-66% inhibition at 300 uM, far above clinical exposure) but the US PI omits the concentration qualifier. The clinical creatinine rise is attributed to MATE1/MATE2-K INHIBITION blocking renal creatinine secretion — again inhibition, not substrate status.

**olaparib**
  - `OATP1B1` **unknown**
    - *THE SINGLE MOST CONSEQUENTIAL UNRESOLVED ITEM IN THE PARP SET, and the contradiction is INTERNAL TO ONE DOCUMENT — the EMA Lynparza EPAR (EMA/CHMP/789139/2014) section 2.3.3. Study 13ASTRUKP7S2, in HEK293 overexpressing OATP1B1/1B3/OCT1 at 1-100 uM with working positive controls, found 'no evidence that olaparib was a substrate for any of OATP1B1, OATP1B3 or OCT1'. Study KMX042 in the same document suggested possible OATP1B1 involvement in hepatic uptake, but the assessor records that it 'was assessed qualitatively and not quantitatively' and that 'no conclusion can be drawn due to the variability of the test system'. Coding this 'no' would over-state the evidence.*
  - `OATP1B3` **unknown**
  - `OCT1` **unknown**
  *Notes:* You cannot read uptake-substrate status off the FDA label for this drug: the LYNPARZA US PI section 12.3 lists olaparib only as an INHIBITOR of BCRP, OATP1B1, OCT1, OCT2, OAT3, MATE1 and MATE2-K, and as a substrate+inhibitor of P-gp, without stating uptake-substrate status at all. PMID 27684210 (McCormick & Swaisland 2016 Xenobiotica) says only that olaparib 'was found to be a substrate of MDR1 but not of several other transporters' without enumerating them. Inhibitor IC50s, for completeness and NOT as substrate evidence: OATP1B1 20.3-27.1 uM, OCT1 37.9 uM, OCT2 19.9 uM, OAT3 18.4 uM, MATE1 5.50 uM.

### Corrections this pass makes to earlier reports in this repo

Three, all recorded because a repo that only ever adds findings is not checking itself:

1. **`results/efflux_impact.md` states adavosertib "ACCUMULATES in human glioma because OATP1A2-mediated uptake dominates in the acidic tumour microenvironment".** The accumulation is real and well measured. The *explanation* is not: no pH experiment was ever run on OATP1A2 in the source paper, its pH work is entirely on MDCKII efflux, its own mechanism is reduced efflux at low pH, and its PBPK model holds uptake clearance constant at a value fitted to reproduce the answer.
2. **`results/efflux_impact.md` records cefazolin as "a confirmed OAT1 substrate (PMID 15618660, Ki 72 uM vs PAH)".** Two errors compounded: that PMID is **rat** Oat1, and human OAT1 does not transport cefazolin (uptake within two-fold of control, PMID 16098483); and a Ki is *inhibition* of PAH transport, not substrate status. Cefazolin's real carrier is OAT3 (PMID 14984259).
3. **`results/gist_imatinib.md` names "active hepatic uptake (OCT1/OATP)" as the missing physics.** The OCT1 half is falsified (§5.1). The OATP half survives as substrate status but is quantitatively insufficient (§5.2).

## 5. Scorecard 1 — GIST tumour:normal, the target uptake was named for

The prediction is the published one with one factor added per side:

```
  before:  R = E_tumour / E_normal
  after :  R = (U_tumour * E_tumour) / (U_normal * E_normal)
```

`E` is unchanged (1-D Krogh slab, L = 35.9 um from the doxorubicin-anchored `alpha` fitted once to an unrelated measurement), the intervessel half-spacings are unchanged. The ONLY new quantity is `U`.

| organ | side | tissue | U | carriers applied | reason |
|---|---|---|---|---|---|
| Liver | tumour | liver_metastasis | **1.000** | -- | no_transporter_in_this_tissue |
| Liver | normal | normal_liver | **1.000** | -- | substrate_but_no_published_ratio |
| Intestine | tumour | intestinal_tumour | **1.000** | -- | no_transporter_in_this_tissue |
| Intestine | normal | normal_intestine | **1.000** | -- | no_transporter_in_this_tissue |

**`U_normal` for liver resolves to 1.000**, so the prediction is unchanged at 7.92x. The reason is not that imatinib fails to be a hepatic uptake substrate — OATP1B3 is a confirmed one — but that the only place its uptake fold was ever published is Figure 1A of PMID 23340295, as bars with asterisks and no number in the text, tables or caption. This is the doxorubicin situation from `results/efflux_impact.md` repeated exactly, on the anchor drug of the other scorecard.

| patient | organ | measured | pred before | fold before | pred after | fold after |
|---|---|---|---|---|---|---|
| A | Liver | **0.0231** | 0.5796 | 25.06x | 0.5796 | 25.06x |
| AA | Intestine | **0.7159** | 0.6446 | 1.11x | 0.6446 | 1.11x |
| J | Liver | **0.1418** | 0.5796 | 4.09x | 0.5796 | 4.09x |
| K | Liver | **0.0439** | 0.5796 | 13.20x | 0.5796 | 13.20x |
| L | Liver | **0.0433** | 0.5796 | 13.39x | 0.5796 | 13.39x |
| U | Liver | **0.0473** | 0.5796 | 12.26x | 0.5796 | 12.26x |

### 5.1 The named mechanism is falsified, not merely unquantified

`results/gist_imatinib.md` concluded: *"The dominant term is not the tumour being hard to reach — it is normal liver actively concentrating a weak base via hepatic uptake transporters and lysosomal trapping."* The uptake half of that sentence was checked rather than adopted, and it does not survive.

| design | result | source |
|---|---|---|
| 5 OCT1-expressing lines (HEK, MDCK, V79, CHO) + oocytes | **no increase over vector control**; MPP+ positive control up 8-fold in the SAME oocytes | PMID 24352644 |
| decynium-22 and prazosin on OCT1 cells | **failed to inhibit** imatinib uptake | PMID 24352644 |
| shRNA OCT1 knockdown | **no effect** | PMID 24352644 |
| **Oct1(-/-) mice, plasma and HEPATIC imatinib** | **no change vs wild-type** | PMID 24352644 |
| oocytes, 0.2 uM | **not a substrate, P = 0.21**; HEK293 only 1.20-fold | PMID 18483382 |
| reduced-function SLC22A1 variants, oral clearance, n = 74 | **no effect, P = 0.99** | PMID 18483382 |
| oocytes + HEK | **not a substrate of OCT1, OCT2 or OCT3** | PMID 30734690 |

The large clinical "OCT-1 activity" literature (PMID 16597591, 17761829, 17568400; PMID 18559609 gives a prazosin-inhibitable IUR of 6.38 ng/200,000 cells) is real and genuinely predicts response — but every one of those studies measures an **inhibitor-difference in non-transfected cells**, never transport by OCT1. PMID 26108972 shows in HEK293, K562, SD-1 and **GIST-T1** that the prazosin-sensitive component is **lysosomal sequestration**, not SLC22A1. And imatinib inhibits OCT1 with IC50 0.095 uM — the precise property that makes an inhibitor-difference assay read like substrate evidence. **Inhibitor is not substrate**, and that single distinction is the difference between this row being "no" and being "yes".

### 5.2 And the rest of the class is excluded by its own numbers

A factor of 1.0 for want of a published number leaves open the possibility that the mechanism is right and merely unmeasured. So the sharper question is asked: **what happens at the largest value any published source can support?**

- Upper bound used: **1.5x**. Source: PMID 23340295 (Zimmerman 2013 Clin Cancer Res), Figure 1A: imatinib below the panel maximum of ~1.5-fold (sorafenib); concordant with the only numeric imatinib uptake ratio published for any carrier, 1.5-fold on OATP1A2 at pH 7.4, PMID 21508937 (Eechoute 2011 Clin Pharmacol Ther).
- Value the measurement DEMANDS: **11.73x**.
- **Shortfall: 7.8x.**

| U_normal | predicted liver ratio | geomean fold error (n = 6) | within 3x |
|---|---|---|---|
| 1.00 | 0.5796 | 7.92x | 1/6 |
| 1.20 | 0.4830 | 6.80x | 1/6 |
| 1.50 **<- literature ceiling** | 0.3864 | 5.65x | 2/6 |
| 2.00 | 0.2898 | 4.44x | 2/6 |
| 3.00 | 0.1932 | 3.17x | 2/6 |
| 5.00 | 0.1159 | 2.21x | 5/6 |
| 8.00 | 0.0724 | 1.75x | 5/6 |
| 11.73 **<- what the data demands** | 0.0494 | 1.45x | 6/6 |
| 15.00 | 0.0386 | 1.48x | 5/6 |
| 20.00 | 0.0290 | 1.71x | 5/6 |

Measured liver geomean: **0.0494**. At the literature ceiling the prediction is 0.3864 — still 7.8x too high — and the structural floor only falls from 0.208 to 0.139, leaving 4/5 measured liver ratios still below it.

**Excluded by its own literature: True.** This is the same *class* of argument that made the original GIST run valuable — there, no penetration length could reach the measurement; here, no published hepatic uptake ratio can. A mechanism can be ruled out by a bound without ever being implemented at a value, and doing so is cheaper and more honest than fitting one.

### 5.3 Does the GIST verdict survive the choices we made?

The published run showed the intervessel half-spacings are the largest single lever on the transport prediction. Every row re-solves the whole prediction.

| axis | setting | U_normal | GM before | GM after | improved? |
|---|---|---|---|---|---|
| transporter_source | prefer=in_vitro, allow_circular=False | 1.00 | 7.92x | 7.92x | no |
| transporter_source | prefer=in_vitro, allow_circular=True | 1.00 | 7.92x | 7.92x | no |
| transporter_source | prefer=in_vivo, allow_circular=False | 1.00 | 7.92x | 7.92x | no |
| transporter_source | prefer=in_vivo, allow_circular=True | 1.00 | 7.92x | 7.92x | no |
| geometry | baseline geometry | 1.00 | 7.92x | 7.92x | no |
| geometry | tumour half-spacing 30 um (band end) | 1.00 | 11.22x | 11.22x | no |
| geometry | tumour half-spacing 120 um (band end) | 1.00 | 5.22x | 5.22x | no |
| geometry | normal-liver half-spacing 10 um (band end) | 1.00 | 7.83x | 7.83x | no |
| geometry | normal-liver half-spacing 15 um (band end) | 1.00 | 8.03x | 8.03x | no |
| geometry | penetration length L = 5 um | 1.00 | 3.67x | 3.67x | no |
| geometry | penetration length L = 400 um (uncalibrated) | 1.00 | 12.87x | 12.87x | no |

**0 of 11 configurations improve the GIST prediction.** The verdict is not an artefact of the geometry, the source preference or the circularity guard.

### 5.4 What this does NOT establish, at n = 6

The published run's caveat applies unchanged: **n = 6 matched pairs, five of them the same tissue pairing**, i.e. roughly ONE independent discriminating comparison. That was enough to falsify a transport model and it is enough to bound a replacement, but:

- The upper-bound argument inherits every assumption of the delivery model it sits on top of — the Krogh geometry, the literature half-spacings, and the doxorubicin-anchored `alpha`. It shows hepatic uptake cannot rescue *this* model, not that hepatic uptake is unimportant.
- Ion suppression differs between liver parenchyma and GIST stroma and does not cancel in the ratio. It is an uncontrolled confounder inflating the measured contrast, and therefore inflating the required `U_normal`.
- `U_tumour = 1.0` is an inference from lineage and from HCC downregulation (§8), not a measurement on these sections.

## 6. Scorecard 2 — the six-drug measured human glioma partition

BEFORE is not pasted: it is `glioma_regional_pk.run_scorecard(aggregation="geometric", apply_leak=True)` called live. AFTER changes one thing:

```
  before:  P_eff = clip( bbb_penetrance     + L * D_sample , 0, 1 )
  after :  P_eff = clip( bbb_penetrance * U + L * D_sample , 0, U )
```

| drug | measured | struct. penetrance | disrupted frac | **U** | why | pred before | fold before | pred after | fold after |
|---|---|---|---|---|---|---|---|---|---|
| eribulin | 0.0051 | 0.0353 | 0.00 | **1.0000** | unknown_substrate_status | 0.0353 | 6.944x | 0.0353 | 6.944x |
| cefazolin | 0.1496 | 0.0488 | 0.44 | **1.0000** | no_transporter_in_this_tissue | 0.1526 | 1.020x | 0.1526 | 1.020x |
| methotrexate | 0.1780 | 0.0225 | 0.50 | **1.0000** | substrate_but_no_published_ratio | 0.1415 | 0.795x | 0.1415 | 0.795x |
| temozolomide | 0.1951 | 0.1254 | 0.14 | **1.7000** | applied | 0.1594 | 0.817x | 0.2471 | 1.267x |
| levetiracetam | 0.6972 | 0.8399 | 0.45 | **1.0000** | non_substrate | 0.9466 | 1.358x | 0.9466 | 1.358x |
| adavosertib | 2.4600 | 0.4595 | 0.50 | **1.0000** | only_circular_observable | 0.5785 | 0.235x | 0.5785 | 0.235x |

| aggregate | before | after | direction |
|---|---|---|---|
| GM \|fold error\| | **1.995x** | **2.006x** | WORSE |
| GM signed bias | 1.066x | 1.147x | |
| within 3x | 4/6 | 4/6 | |
| within 10x | 6/6 | 6/6 | |
| worst \|fold error\| | 6.94x | 6.94x | |

**1 of 6 drugs moves at all.** The `why` column gives which of the eight distinct reasons produced each 1.0, and the distinctions matter: `non_substrate` (measured not to be one) is a finding, `substrate_but_no_published_ratio` is a coverage gap, `no_transporter_in_this_tissue` is the mechanism working correctly, and `only_circular_observable` is a number we refused.

### 6.1 Adavosertib — the pre-registered target, blocked twice

The efflux run identified adavosertib as the case where efflux was most clearly the wrong mechanism and named uptake as the right one. This run was written to test that. It cannot.

**The substrate status verifies.** Li et al. did not merely hypothesise it: uptake of AZD1775 in HEK293 overexpressing OATP1A2 was significantly higher than vector control (P = 0.001) and was reduced to control level by 100 uM sulfobromophthalein (P = 0.005), with an estrone-3-sulfate positive control and with OATP2B1 and OAT3 negative in the same experiment. That is a properly controlled positive, and PMID 16357150 puts OATP1A2 in the glioma blood-tumour-barrier endothelium — exactly where a penetrance multiplier acts.

**Block 1: the number does not exist.** The uptake fold appears only as bars in Figure 4A. No value in the text, tables or abstract; no Km, no Vmax; PubMed `AZD1775 AND OATP1A2` returns exactly one record and it publishes only P-values. Under this repo's standing rule — the one that forced doxorubicin to 1.0 in `results/efflux_impact.md` — a figure is not a number.

**Block 2: the only number that does exist is circular, and back-fitted.** See §7.

So adavosertib's factor is 1.00 and its fold error is unchanged at 0.235x. **The mechanism is right, the measurement is right, and the model can use neither.**

### 6.2 Temozolomide is the only mover, and it moves the wrong way

TMZ carries the table's other usable ratio — GLUT3, 1.70x (PMID 40813374) — and it makes the prediction worse, 0.817x -> 1.267x, taking a drug that was slightly under-predicted to slightly over-predicted. This is the temozolomide-and-efflux pattern in reverse: there, a genuine but tiny 1.5-fold efflux term hurt because TMZ was already under-predicted; here a 1.70-fold uptake term overshoots.

It is also the weakest row in the table and is reported as such: a single unreplicated 2025 paper, whose companion affinity data are internally inconsistent by ~1000x (GLUT3 glucose Kd reported at 1-10 uM against an established Km of ~1.4 mM), and whose knockdown arm implies a larger ratio (2.19) than the overexpression arm used here (1.70). Leave-one-out below shows the scorecard verdict does not depend on it either way.

### 6.3 Which drug is doing what

| drug dropped | GM \|fold\| before | GM \|fold\| after | uptake helps without it? |
|---|---|---|---|
| eribulin | 1.554x | 1.565x | no |
| cefazolin | 2.281x | 2.297x | no |
| methotrexate | 2.187x | 2.202x | no |
| temozolomide | 2.199x | 2.199x | no |
| levetiracetam | 2.154x | 2.169x | no |
| adavosertib | 1.714x | 1.726x | no |

### 6.4 Is the change even distinguishable at n = 6?

- drugs that move at all: **1 of 6**; of those, 1 worsen and 0 improve
- mean change in `|log10 fold error|`: **+0.0025** (positive = worse)
- exact sign test: p = 1.000; Wilcoxon signed-rank: p = 0.317

**with 1 of 6 drugs moving at all, neither test can reach significance at any effect size; both are reported to make that explicit.** No p-value here should be read as evidence for anything. What the scorecard establishes is not that uptake hurts; it is that the term cannot be *exercised* on five of the six drugs, and why each one is blocked.

## 7. The circularity guard, and exactly what it costs

A number one has decided not to use should be printed, not hidden.

| drug | carrier | tissue | withheld value | quantity | circular against | factor default | factor if allowed |
|---|---|---|---|---|---|---|---|
| adavosertib | OATP1A2 | normal_brain | **3.20** | unbound tumour-to-plasma concentration ratio, median 3.2 (range 1.3-24.4) in n = 20 patien... | glioma_regional_pk | 1.00 | 3.20 |
| adavosertib | OATP1A2 | glioma | **3.20** | unbound tumour-to-plasma concentration ratio, median 3.2 (range 1.3-24.4) in n = 20 patien... | glioma_regional_pk | 1.00 | 3.20 |
| adavosertib | OATP1A2 | kidney | **3.20** | unbound tumour-to-plasma concentration ratio, median 3.2 (range 1.3-24.4) in n = 20 patien... | glioma_regional_pk | 1.00 | 3.20 |

**Why this is circular and not merely "in vivo".** The efflux module used in vivo mouse knockout fold changes without complaint, and rightly: a knockout contrast in a mouse is a mechanistic input, measured against a genetic control, in a different species from the target. Li 2017's `Kp,uu` is none of those things. It is a human glioblastoma tumour-to-plasma drug concentration ratio — the same quantity, in the same tissue, in the same disease as the scorecard target.

The two studies ARE independent (Li 2017: Karmanos/Barrow, n = 20, unbound ratio, median 3.2 range 1.3-24.4; Lee 2025 PMID 39820282: Adult Brain Tumor Consortium, total homogenate ratio 4.18 +/- 3.36 enhancing and 0.74 +/- 0.63 non-enhancing), so this is not double-counting one dataset. Enabling it would be a real cross-cohort replication. **It would not be a prediction from mechanism**, which is the only thing this repo scores.

And there is a second, independent reason to refuse it. Li 2017's PBPK uptake clearance of **44 L/h was "assigned based on sensitivity analysis"** — i.e. tuned so the model reproduces 3.2. Passive plus efflux alone give Kp,uu 0.78-0.92 (max 1.03 at sensitivity extremes). So the value is not an independent measurement of uptake at all; it is the residual required to close the mass balance, and importing it would import a fitted parameter under the name of a literature value.

For completeness, the number it would have bought: GM abs fold error 1.995x -> **1.695x**, adavosertib to 0.646x. **That is an upper bound obtained by using the answer, and it is not claimed as a result.**

Note what this implies about the headline. If the guard were dropped, the glioma scorecard would improve — and the improvement would be worthless. That is the whole reason the guard is a code path with a reason code rather than a paragraph of good intentions.

## 8. The mechanism's premise, tested out of sample

The GIST argument rests on a claim the GIST dataset cannot test, because no transporter stain accompanies the deposit: **hepatic uptake carriers are present in normal liver and effectively absent from a tumour sitting in it.** The hardest possible case for that claim is hepatocellular carcinoma — a tumour that actually derives from hepatocytes. If the premise failed even there, it would certainly fail for a GIST metastasis.

| gene | carrier | tumour vs adjacent non-tumour liver | source |
|---|---|---|---|
| SLC22A1 | OCT1 | protein significantly DOWN in HCC vs paired adjacent normal liver, n = 71 paired, tissue-microarray IHC, P < 0.0001 | PMID 22196450 (Schaeffeler 2011 Genome Med) |
| SLC22A1 | OCT1 | 39 advanced HCC biopsies: absent 10%, weak 38%; and of the tumours that DO express it only ~1/3 have it at the plasma membrane (15% PM vs 36% cytosolic) — mislocalisation, not just loss of message | PMID 28178663 (Geier/Marin 2017 Oncotarget, TRANSFER) |
| SLC22A1 | OCT1 | aberrant splice variants that fail to reach the plasma membrane; >= 1 inactivating variant in 48% of HCC and 40% of cholangiocarcinoma | PMID 23532667 (Herraez/Marin 2013 Hepatology) |
| SLCO1B3 | OATP1B3 | decreased in 60% of 13 paired HCC / adjacent non-tumour samples, mRNA and protein, via HNF3-beta repression | PMID 14739090 (Vavricka 2004 J Hepatol) |
| SLCO1B1 | OATP1B1 | CONTRADICTS THE CLEAN STORY: OATP1B1 (OATP-C) was NOT significantly decreased in the same 13 paired samples | PMID 14739090 (Vavricka 2004 J Hepatol) |

All carriers downregulated: **False** — and the dissent is recorded in the table rather than dropped: OATP1B1 was **not** significantly decreased in the same paired samples that found OATP1B3 down in 60% (PMID 14739090). The premise is directionally supported and not uniform. No number from these papers enters the model, and HCC is not GIST.

## 9. Does the glioma verdict depend on anything we chose?

Every discretionary choice, crossed: source preference, the circularity guard, where the factor multiplies, and which ceiling applies. Sixteen configurations. There is no functional-form axis because the form is unique (§2) — the efflux run needed eight configurations partly to cover a form ambiguity that does not exist here.

| prefer | circular | placement | ceiling | GM before | GM after | 3x | 10x | worst after | ada after | improved? |
|---|---|---|---|---|---|---|---|---|---|---|
| in_vitro | False | transcellular | uptake | 1.995x | 2.006x | 4->4 | 6->6 | 6.94x | 0.235x | no |
| in_vitro | False | transcellular | unity | 1.995x | 2.006x | 4->4 | 6->6 | 6.94x | 0.235x | no |
| in_vitro | False | total | uptake | 1.995x | 2.037x | 4->4 | 6->6 | 6.94x | 0.235x | no |
| in_vitro | False | total | unity | 1.995x | 2.037x | 4->4 | 6->6 | 6.94x | 0.235x | no |
| in_vitro | True | transcellular | uptake | 1.995x | 1.695x | 4->5 | 6->6 | 6.94x | 0.646x | yes |
| in_vitro | True | transcellular | unity | 1.995x | 1.831x | 4->5 | 6->6 | 6.94x | 0.407x | yes |
| in_vitro | True | total | uptake | 1.995x | 1.678x | 4->5 | 6->6 | 6.94x | 0.752x | yes |
| in_vitro | True | total | unity | 1.995x | 1.859x | 4->5 | 6->6 | 6.94x | 0.407x | yes |
| in_vivo | False | transcellular | uptake | 1.995x | 2.006x | 4->4 | 6->6 | 6.94x | 0.235x | no |
| in_vivo | False | transcellular | unity | 1.995x | 2.006x | 4->4 | 6->6 | 6.94x | 0.235x | no |
| in_vivo | False | total | uptake | 1.995x | 2.037x | 4->4 | 6->6 | 6.94x | 0.235x | no |
| in_vivo | False | total | unity | 1.995x | 2.037x | 4->4 | 6->6 | 6.94x | 0.235x | no |
| in_vivo | True | transcellular | uptake | 1.995x | 1.695x | 4->5 | 6->6 | 6.94x | 0.646x | yes |
| in_vivo | True | transcellular | unity | 1.995x | 1.831x | 4->5 | 6->6 | 6.94x | 0.407x | yes |
| in_vivo | True | total | uptake | 1.995x | 1.678x | 4->5 | 6->6 | 6.94x | 0.752x | yes |
| in_vivo | True | total | unity | 1.995x | 1.859x | 4->5 | 6->6 | 6.94x | 0.407x | yes |

**0 of 8 non-circular configurations improve the glioma scorecard** (8 of 16 including the circular ones). The verdict is not an artefact of the placement, the ceiling or the source preference.

## 10. The redirect — and a correction to this chain's own reasoning

Two mechanisms have now been implemented, cited and scored against these datasets, and neither worked. That would be a dead end if the elimination argument that got us here were sound. It is not, in one specific place, and that is the most useful thing this run found.

**What Moncayo proved.** `results/moncayo_explant.md` §4.3: with `phi_lys` uniform in space, the retardation `R` is uniform, so at equilibrium `C_total = R * C_bath` **everywhere** and the dispersion goes to zero however large `R` becomes. Correct, and this module does not dispute it.

**What it does not prove.** That is a statement about dispersion *within one tissue*. The GIST target is a ratio *between two tissues*, and two tissues with different acidic-compartment content have different `R`:

```
  C_tumour / C_normal = (R_t * E_t) / (R_n * E_n)      does NOT cancel
```

This is the same algebra that lets a tissue-differential uptake term move the ratio where a scalar efflux term could not. **The chain eliminated ion-trapping for GIST using a proof about a different quantity.** Both the efflux report (§5.1) and this task's framing inherited that error.

**And the evidence points there.** The same group whose experiments falsified imatinib/OCT1 concluded that **lysosomal sequestration, not SLC22A1, determines intracellular imatinib** — in HEK293, K562, SD-1 and, to the point, **GIST-T1** (PMID 26108972 (Burger/Sparreboom 2015 Mol Pharmacol)). The prazosin-sensitive signal underpinning the entire clinical "OCT-1 activity" literature is, in that work, lysosomal. So the mechanism that explains why the uptake story looked true is the mechanism the chain discarded.

### Sizing it — is the required asymmetry reachable, or absurd?

Not a model, not a prediction, and nothing fitted. This uses the Henderson-Hasselbalch trapping ratio already implemented in `validation/moncayo_explant.py` and imatinib's own basic pKa already used by `validation/gist_imatinib.py`. An absurd answer would rule the mechanism out cheaply; a plausible one identifies the next piece of work.

| quantity | value |
|---|---|
| imatinib basic pKa (already in the repo) | 8.07 |
| lysosome:cytosol trapping ratio at pH 4.7 / 7.2 | **278.8** |
| `R_normal/R_tumour` required by the measurement | **11.73** |
| implied hepatocyte lysosomal volume fraction (tumour ~0) | **3.86%** |
| `phi_lys` the Moncayo module uses for generic tissue | 1.4% |
| largest `phi_lys` the Moncayo module swept | 5% |
| **reachable inside that swept range?** | **True** |

The required 3.9% sits inside the range `moncayo_explant.py` already sweeps and is high but not absurd for a hepatocyte, which is among the most lysosome-rich cells in the body. Note carefully what this is and is not: **it is a plausibility check on a single number, not evidence.** It says the candidate cannot be dismissed on magnitude, which is exactly the test the hepatic-uptake candidate just failed in §5.2. That asymmetry — one candidate excluded by a bound, the other not — is the entire content of this section.

**Two things would be needed to turn it into a test**, and neither is available today: a measured acidic-compartment or lysosomal-content contrast between hepatocytes and GIST tissue (a literature pass for this project established that no measured mg/g or mol% acidic-phospholipid content exists for tumour tissue, and no published PBPK model carries tumour-specific lysosomal parameters), and an independent tumour-in-liver dataset so the test is not run on the same six pairs that generated the hypothesis.

## 11. Verdict and recommendation

### Recommendation: **default OFF, and more strongly than efflux.**

Efflux was left off after being tested and found wrong. Uptake is left off after being tested and found **inapplicable**: on 16 of 18 drugs the literature supplies no usable ratio, on the GIST anchor the named carrier is falsified and the class is bounded out, and on the glioma anchor the one usable number would be circular. Nothing in the shipped pipeline calls `molecular.uptake` and nothing was changed to accommodate it; it is available to any caller who imports it, as a documented mechanism rather than an endorsed default.

**Do not read this as "uptake transport is unimportant."** It plainly is important — a 4.1-fold liver-to-plasma change in Oatp1a/1b-/- mice for doxorubicin is real. What is shown is narrower and firmer: a literature-sourced, zero-free-parameter uptake term cannot be built for these drugs today, and where it can be built it does not improve these two measured human datasets.

### What this run establishes

1. **Uptake is implementable with zero free parameters and, unlike efflux, one unambiguous functional form** — all three published observables collapse to `UR = 1 + CL_up/P`, verified from the mass balance.
2. **A tissue-differential factor is structurally capable of moving a tumour:normal ratio**, where the efflux scalar provably was not. Both identities are verified to machine precision. The architecture is correct; the inputs are missing.
3. **The mechanism `results/gist_imatinib.md` named is falsified** by five cell lines, oocytes with a positive control, shRNA knockdown, and Oct1(-/-) mice with unchanged hepatic imatinib.
4. **The remaining hepatic-uptake class is excluded by a bound**: 11.7x demanded, 1.5x available, 7.8x short. Ruling a mechanism out by a literature ceiling, without ever fitting it, is the cheapest honest form of falsification available and it is the main methodological contribution here.
5. **Literature coverage is the binding constraint for the second run running, and now the reason is structural.** Uptake papers report Km/Vmax and inhibitor IC50s; a penetrance multiplier needs a ratio against passive permeability in the same system. Those are different measurements and the field mostly publishes the first.
6. **Three factual corrections to earlier reports in this repo** (§4), including one — the acidic-TME/OATP1A2 story — that this repo had been repeating on a misreading of its own cited source.

### What this run does NOT establish

- **Nothing at adequate power.** Six matched pairs with one independent tissue contrast; six drugs of which one moves.
- **Not that adavosertib's uptake mechanism is wrong.** It is documented, controlled, in the right direction, and unquantified. If the value behind Figure 4A of PMID 28928160 were published, the glioma half of this test could actually be run — one number, from one experiment already performed.
- **Nothing about the tumour side of the GIST contrast.** `U_tumour = 1.0` is an inference from lineage and HCC data, not a measurement.
- **Nothing about lysosomal sequestration.** §10 sizes it; it does not test it. A plausible magnitude is not evidence.

### Where the chain stands now

| candidate | status |
|---|---|
| passive reaction-diffusion transport | **falsified** for tissue partition — below the structural floor at any L (`gist_imatinib.md`) |
| scalar lysosomal ion-trapping, WITHIN a section | **falsified** — uniform R gives no dispersion at equilibrium (`moncayo_explant.md`) |
| active efflux (P-gp / BCRP) | **tested, does not help**; provably inert on a tissue ratio; wrong sign in liver (`efflux_impact.md`) |
| carrier-mediated uptake (OATP / OCT / OAT) | **tested, does not help**; named carrier falsified; class excluded by a literature bound (this run) |
| **tissue-differential acidic-compartment / lysosomal content** | **OPEN, and wrongly eliminated earlier.** Direct evidence in GIST-T1 (PMID 26108972); required magnitude plausible (§10); no measurement exists |
| spatial transporter or acidic-compartment FIELD (within-tumour structure) | **OPEN**, unchanged from `moncayo_explant.md`; no measurement in usable units exists |

The last two are different objects and should not be conflated again. A per-tissue *constant* is enough to move a tumour:normal ratio and cannot produce structure inside a tumour; a spatial *field* is needed for the Moncayo target and is not needed for the GIST one. This run's contribution to that distinction is that the GIST-shaped version of ion-trapping was never actually ruled out.

---

*Generated by `tmesim/validation/uptake_impact.py`. Additive only: no file under `core/`, `physics/`, `molecular/` (other than the new `molecular/uptake.py`), `reconstruct/`, `sim/` or `treatment/` was modified, and the bit-exactness of the untouched delivery path is verified in §1.*
