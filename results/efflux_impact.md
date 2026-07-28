# Active efflux (P-gp / BCRP) — implemented, cited, and it does NOT help

**Verdict up front: implemented as a zero-free-parameter, fully cited mechanism, scored against both measured human datasets, and left DEFAULT-OFF because it makes the one scorecard it can act on worse and cannot act on the other at all.**

| scorecard | before | after efflux | verdict |
|---|---|---|---|
| human glioma partition, n = 6 drugs (GM abs fold error) | **1.995x** | **2.424x** | WORSE |
| .. within 3x | 4/6 | 4/6 | unchanged |
| .. within 10x | 6/6 | 5/6 | WORSE |
| .. worst case | 6.94x | 16.23x | WORSE |
| GIST tumour:normal, n = 6 pairs (geomean fold error) | **7.92x** | **7.92x** | IDENTICAL — the term cancels exactly |

Three findings, in order of how much they are worth:

1. **A scalar per-drug efflux factor is structurally incapable of moving a tumour:normal ratio.** It cancels: `R = (f·E_t)/(f·E_n) = E_t/E_n`. Verified to machine precision on all six GIST pairs (§5). This is the SAME class of structural limit the Moncayo explant run established for scalar ion-trapping, arrived at independently, and it means the GIST diagnosis cannot be answered by any per-drug property — only by a per-TISSUE transporter field.
2. **On the glioma scorecard the term is contradicted by a direct human measurement of the same drug in the same tissue.** Adavosertib carries the best-evidenced efflux ratio in the table (14.1, two assay types, two groups), and the same paper reports human GBM tumour Kp,uu of median 3.2 — adavosertib ACCUMULATES in human glioma because OATP1A2 UPTAKE dominates in the acidic tumour microenvironment. Applying the efflux term drives it from 0.24x to 0.062x.
3. **Literature coverage is the binding constraint, not the physics.** Of 18 drugs, 7 carry a usable published number and 5 are CONFIRMED substrates with no usable ratio anywhere — including doxorubicin, the textbook P-gp substrate and the molecule this repo's perivascular constant was fitted to. A mechanism that silently does nothing for half the library is not ready to be on by default.

## 1. The one touched function is bit-identical when off

`treatment.drugs.delivery_source` gained a single keyword, `efflux=False`. `core`, `physics` and `sim.engine` are untouched.

- **180/180 comparisons bit-for-bit identical** (`np.array_equal`, not `allclose`) across 4 drug names x 5 penetrances x 3 plasma concentrations x 3 branches (implicit default, explicit `efflux=False`, and the barrier branch).
- All identical: **True**.
- And the check is not vacuous: with `efflux=True` on a strong substrate the source does change (True).

## 2. The literature-sourced substrate table

18 drugs: all ten in `molecular.library.GBM_DRUGS`, the four glioma-scorecard drugs outside it, imatinib (GIST), and the three PARP inhibitors (Moncayo). Every non-`unknown` cell carries a PMID. `unknown` means searched and not found, and is treated downstream as a non-substrate — but is reported separately, because *measured not to be a substrate* and *nobody has looked* are different states of knowledge.

- P-gp (ABCB1): **10 yes / 3 no / 5 unknown**
- BCRP (ABCG2): **8 yes / 5 no / 5 unknown**
- Confirmed substrate WITH a usable published number: **7**
- Confirmed substrate WITHOUT one: **5** (factor forced to exactly 1.0 — we do not invent ratios)
- Measured non-substrate: **6**

| drug | P-gp | BCRP | in vitro ER | in vivo KO fold | ER used | **factor** | source | status |
|---|---|---|---|---|---|---|---|---|
| temozolomide | yes | yes | -- | 1.5 | 1.50 | **0.6667** | in_vivo | applied |
| lomustine | unknown | unknown | -- | -- | -- | **1.0000** | -- | unknown_substrate_status |
| carmustine | unknown | unknown | -- | -- | -- | **1.0000** | -- | unknown_substrate_status |
| procarbazine | unknown | unknown | -- | -- | -- | **1.0000** | -- | unknown_substrate_status |
| methotrexate | no | yes | -- | -- | -- | **1.0000** | -- | substrate_but_no_published_er |
| doxorubicin | yes | no | -- | -- | -- | **1.0000** | -- | substrate_but_no_published_er |
| cisplatin | no | no | -- | -- | -- | **1.0000** | -- | non_substrate |
| etoposide | yes | unknown | -- | -- | -- | **1.0000** | -- | substrate_but_no_published_er |
| vincristine | yes | no | -- | -- | -- | **1.0000** | -- | substrate_but_no_published_er |
| bevacizumab | no | no | -- | -- | -- | **1.0000** | -- | non_substrate |
| levetiracetam | unknown | no | -- | -- | -- | **1.0000** | -- | unknown_substrate_status |
| cefazolin | unknown | yes | -- | -- | -- | **1.0000** | -- | substrate_but_no_published_er |
| adavosertib | yes | yes | 14.1 | 25.0 | 14.10 | **0.0709** | in_vitro | applied |
| eribulin | yes | unknown | -- | 30.0 | 30.00 | **0.0333** | in_vivo | applied |
| imatinib | yes | yes | 39.0 | 28.0 | 39.00 | **0.0256** | in_vitro | applied |
| niraparib | yes | yes | -- | 6.5 | 6.50 | **0.1538** | in_vivo | applied |
| rucaparib | yes | yes | -- | 14.6 | 14.60 | **0.0685** | in_vivo | applied |
| olaparib | yes | yes | -- | 8.1 | 8.10 | **0.1235** | in_vivo | applied |

### Citations, per assignment

**temozolomide** — P-gp `yes`: PMID 29852323 (de Gooijer 2018 Neoplasia): WEAK. No clear transport in LLC-Mdr1a monolayers, but Abcb1a/b-/- mouse brain accumulated ~20% more TMZ than WT.
  BCRP `yes`: PMID 29852323: TMZ is a GOOD Bcrp1 substrate — clear vectorial transport in MDCKII-Bcrp1, abolished by elacridar 5 uM.
  In vivo knockout fold: 1.5. in vivo brain-to-plasma fold change, Abcb1a/b;Abcg2-/- vs WT FVB mice, 50 mg/kg i.v.; reproduced pharmacologically by elacridar. PMID 29852323 (de Gooijer 2018 Neoplasia), p < 0.001, plasma clearance unchanged
  *Notes:* NO in vitro ER exists and none can: the authors state TMZ's ~30 min half-life at pH 7.4 makes Papp in cm/s uncalculable. The combined knockout effect is only 1.5x (1.2x for each pump alone) — i.e. TMZ is a real but very weak efflux substrate, which corrects the common assumption that it is not one at all. Goldwirt 2014 (PMID 24867782) asserts ABCB1 transport in its abstract but reports a mdr1a-/- vs +/+ brain:plasma ratio of 1.1, i.e. no effect; read its number, not its claim.

**lomustine** — P-gp `unknown` (no source found)
  BCRP `unknown` (no source found)
  *Notes:* No transport assay of any kind exists — no bidirectional monolayer, no transfected-cell transport study, no Mdr1a/b-/- or Bcrp1-/- mouse PK. The one relevant primary result is a RESISTANCE readout, not a transport measurement: PMID 36252775 (Radtke 2022 Pharmacol Res) found that CRISPR ABCB1 knockout in glioblastoma cells improved the temozolomide response and altered the carmustine response but did NOT alter the lomustine response. That is suggestive of non-substrate status and is recorded here, but a viability assay is not a substrate assay, so the status stays UNKNOWN rather than being upgraded to 'no'.

**carmustine** — P-gp `unknown` (no source found)
  BCRP `unknown` (no source found)
  *Notes:* No primary substrate study exists. Searches for carmustine/BCNU against P-glycoprotein, ABCB1, ABCG2/BCRP, efflux transporters and knockout-mouse PK returned zero substrate assays; every carmustine+MDR1 hit is a gene-therapy chemoprotection study in which MDR1 is transduced to protect marrow against a DIFFERENT drug while MGMT protects against BCNU (e.g. PMID 10933194). Mechanistically BCNU is small, lipophilic and hydrolyses within minutes, which makes it both a poor efflux candidate and hard to assay — but that is reasoning, not a citation.

**procarbazine** — P-gp `unknown` (no source found)
  BCRP `unknown` (no source found)
  *Notes:* No primary literature found. The broadest productive query ('procarbazine transport') returns ~13 records, none addressing ABC-transporter substrate status. Circumstantially procarbazine is small (MW 221), lipophilic, and used in the CNS-directed PCV regimen, which is consistent with low efflux liability — an inference, not data, and not a basis for a numeric term.

**methotrexate** — P-gp `no`: PMID 38981917 (Volpe 2024, FDA): MTX efflux ratio measured in C2BBe1 wild-type vs single/double CRISPR knockouts falls in BCRP-KO and MRP2-KO cells but NOT in P-gp-KO cells.
  BCRP `yes`: PMID 14500392 (Volk & Schneider 2003 Cancer Res): ATP-dependent vesicular transport by WILD-TYPE (R482) BCRP, Km 680 uM, Vmax 2400 pmol/mg/min; PMID 12874005 (Chen 2003): Km 1.34 mM, Vmax 687 pmol/mg/min. MTX-Glu2/Glu3 also transported. In vivo confirmation PMID 15313923, PMID 19996279.
  *Notes:* A confirmed BCRP substrate with NO usable number for this purpose. The available in vivo figure (1.8-1.9x, PMID 15313923) is a PLASMA CLEARANCE change in Bcrp1-/- mice, not a brain- or tissue-to-plasma penetration fold change, so it does not measure the quantity this module needs and is deliberately NOT used as fold_ko. Note also the mirror-image of doxorubicin: the R482T/R482G drug-selected BCRP mutants LOSE MTX transport while gaining doxorubicin transport. MTX is additionally an OAT/OATP/RFC substrate, which is out of this module's scope and is part of why an efflux-only term is incomplete for it. One contradicting P-gp claim exists (PMID 30977678) from Caco-2 accumulation rather than bidirectional flux; not weighted.

**doxorubicin** — P-gp `yes`: PMID 3472246 (Ueda 1987 PNAS): MDR1 cDNA transfection confers doxorubicin resistance — the canonical demonstration. In vivo: PMID 10681369 (van Asperen 2000), biliary secretion of unchanged doxorubicin falls 13.3% -> 2.4% of dose in mdr1a-/- mice (~5.5x).
  BCRP `no`: PMID 11559526 (Honjo 2001) and PMID 14612912 (Robey 2003): WILD-TYPE R482 human BCRP does NOT transport doxorubicin; only the R482T/R482G drug-selected mutants do (mutant transfectants 13-71x resistant vs 3-4x for wild-type).
  *Notes:* THE ANCHOR CASE, AND IT CARRIES NO NUMBER. Doxorubicin is the textbook P-gp substrate and the drug this repo's perivascular-uptake constant alpha was originally fitted to, yet no defensible single in vitro ER could be attributed: PMID 39669004 (Kotze 2024) ran it bidirectionally in MDCK-II-MDR1 and states ER >= 2, but publishes ER only as pH curves in figures with values in unretrieved supplementary material. Recording 'ER ~ 2' from a figure would be a guess, so the factor stays 1.0. The widespread belief that doxorubicin is a BCRP substrate traces entirely to mutant cell lines and is coded here as 'no' for wild-type BCRP.

**cisplatin** — P-gp `no`: PMID 22792399 (Stordal 2012 PLoS One): elacridar does not reverse cisplatin resistance in a P-gp overexpressor; PMID 24403508 (Doherty 2014): P-gp-overexpressing KB-8-5-11 cells are COLLATERALLY SENSITIVE to cisplatin; PMID 16534613 (Takara 2006): cisplatin-selected HeLa show DECREASED MDR1 expression/activity.
  BCRP `no`: PMID 18820913, PMID 19638344 — negative in BCRP substrate screens.
  *Notes:* Expectation verified rather than assumed: platinum resistance runs through copper transporters (CTR1/ATP7A/ATP7B) and glutathione/GST, not ABC efflux. One trap avoided: PMID 15139851 (Alqawi 2004) shows cisplatin competes for a photoaffinity label on ABCG2 — that is BINDING, not transport, and must not be coded as substrate evidence.

**etoposide** — P-gp `yes`: PMID 20028753 (Lagas 2010 Clin Cancer Res): P-gp restricts oral (re)uptake and mediates gut-wall excretion; PMID 12649196 (Allen 2003): GF120918 in WT mice raised oral etoposide plasma exposure 4-5x; PMID 11245453: mdr1a/1b;mrp1 triple-KO mice are 3-5x more etoposide-sensitive.
  BCRP `unknown` (no source found)
  *Notes:* Genuinely split by species/system, hence BCRP unknown: mouse Bcrp1 transports etoposide across MDCKII monolayers and confers resistance reversible by Ko143 and GF120918 (PMID 12649196), while human ABCG2 in HEK293 conferred no appreciable etoposide resistance (PMID 14612912). NO NUMBER IS RECORDED. A candidate ER of 5.38 (C2BBe1 wild-type, falling to 3.65 in MDR1-KO and 1.24 in MDR1/MRP2 double-KO, PMID 38981917) could not be verified against the source table, and in any case that wild-type ratio is dominated by MRP2/ABCC2 rather than P-gp — etoposide is primarily an MRP2 substrate, so attributing its whole efflux to P-gp would over-state the term this module models.

**vincristine** — P-gp `yes`: PMID 20501632 (Wang 2010 Neuro-Oncology); PMID 11245453 (Johnson 2001): mdr1a/1b;mrp1 triple-KO mice show 128x increased vincristine toxicity; PMID 11489261.
  BCRP `no`: PMID 18820913 — negative in a BCRP substrate screen.
  *Notes:* Substrate status is not in doubt but the QUANTITATIVE brain evidence is weaker than its reputation. The only brain-specific figure found (PMID 20501632) is a 1.3-3.6x normal-brain partition increase in a TRIPLE Pgp+Mrp1 knockout, so it is not attributable to P-gp and BCRP alone and is not recorded as fold_ko. No numeric monolayer ER was located. Beware also that much of the canonical 'vinca + P-gp' literature is VINBLASTINE (e.g. PMID 7910522), which is not interchangeable with vincristine for a quantitative term.

**bevacizumab** — P-gp `no`: PMID 23421687 (Desai 2013 ACS Chem Neurosci): ER = (P_pass + Vmax/Km)/P_pass is undefined as P_pass -> 0; PMID 27572805 (Pardridge 2016): therapeutic antibodies do not cross the BBB and must be re-engineered for receptor-mediated transcytosis.
  BCRP `no`: PMID 23421687; PMID 27572805 — same argument.
  *Notes:* NOT APPLICABLE rather than measured-negative. ABC efflux transporters act on solutes that partition into the lipid bilayer; a ~149 kDa IgG has essentially no passive transcellular permeability, so the efflux ratio is ill-posed, not merely unmeasured. Bevacizumab's barrier restriction belongs in a size/permeability term, and this repo already handles it there (bbb_penetrance = 0 via the is_antibody and MW >= 1000 gates). efflux_factor short-circuits antibodies to 1.0 before this row is read.

**levetiracetam** — P-gp `unknown` (no source found)
  BCRP `no`: PMID 25645391 (Roemermann 2015 Neuropharmacology): no transport in MDCKII-Bcrp1 or MDCKII-BCRP concentration-equilibrium assays at therapeutic concentrations, in a paper whose lamotrigine positive control DID transport — a well-controlled negative.
  *Notes:* P-gp status is genuinely DISPUTED, and the dispute looks like a species and assay-format artefact, so 'unknown' is the honest code. FOR: PMID 18824002 (Luna-Tortos 2008) found directional transport by human P-gp in a concentration-equilibrium transport assay — used precisely because levetiracetam's high passive permeability masks low-affinity efflux in a conventional bidirectional assay — though it was also MK571-inhibitable, so MRPs contribute. AGAINST: PMID 17045309 (Baltes 2006) found transport by mouse mdr1a but NOT human MDR1; PMID 15120740 (Potschka 2004) found verapamil and probenecid did not raise cortical levetiracetam in rat microdialysis; PMID 17916018 rated it only a weak canine P-gp substrate. No numeric ER exists either way — a CETA reports AUC differences, not a ratio.

**cefazolin** — P-gp `unknown` (no source found)
  BCRP `yes`: PMID 18339814 (Kato 2008 Drug Metab Dispos): ATP-dependent uptake of cefazolin into RAT Bcrp-expressing Sf9 membrane vesicles.
  *Notes:* Two caveats make this row weak in both cells. P-gp: there is NO substrate study at all — four PubMed hits for cefazolin + P-glycoprotein, none testing substrate status — hence unknown, not 'no'. BCRP: the positive is RAT Bcrp, measured as vesicular uptake in a biliary-excretion context, not a human-BCRP monolayer ratio, and human BCRP was not tested. No ER exists, so the factor is 1.0 regardless. Cefazolin's real transporter biology is organic-anion, not ABC-efflux: a confirmed OAT1 substrate (PMID 15618660, Ki 72 uM vs PAH) and an MRP4 substrate (PMID 17344354, Km 80 uM, kidney:plasma up 2.7x in Mrp4-/- mice).

**adavosertib** — P-gp `yes`: PMID 28928160 (Li 2017 Clin Cancer Res); independently confirmed PMID 29147815 (de Gooijer 2018 Invest New Drugs).
  BCRP `yes`: PMID 28928160; PMID 29147815.
  In vitro ER: P-gp 10.6, BCRP 4.5, combined 14.10. bidirectional MDCKII-ABCB1 (ER 10.6) and MDCKII-ABCG2 (ER 4.5) monolayers at basolateral pH 7.4; abolished by cyclosporine A 10 uM and Ko143 1 uM respectively; loperamide/gefitinib positive controls, 85-115% mass recovery. PMID 28928160 (Li 2017 Clin Cancer Res)
  In vivo knockout fold: 25.0. in vivo brain-to-plasma ratio, wild-type vs Abcb1a/b;Abcg2-/- mice. PMID 29147815 (de Gooijer 2018 Invest New Drugs)
  *Notes:* The best-evidenced drug in the table — two assay types, two independent groups — and, for our purposes, the most dangerous. Two warnings. (1) The ERs are pH-DEPENDENT: ABCB1 falls 9.2 -> 3.7 and ABCG2 4.7 -> 2.6 as basolateral pH drops 7.4 -> 6.5, so the acidic tumour microenvironment roughly halves both. The pH 7.4 values are recorded, which is the efflux-unfavourable end for a tumour. (2) DIRECTLY CONTRADICTED IN OUR OWN TISSUE: despite dual efflux, Li 2017 measured human GBM tumour Kp,uu at a MEDIAN OF 3.2 (range 1.3-24.4) — i.e. adavosertib ACCUMULATES in human glioma, because OATP1A2-mediated uptake dominates in the acidic tumour microenvironment. The mouse and human results genuinely disagree and the authors say so explicitly. Applying this row to human glioma is contradicted by a direct human measurement of the same drug in the same tissue.

**eribulin** — P-gp `yes`: PMID 21162698 (Taur 2010 Xenobiotica).
  BCRP `unknown` (no source found)
  In vivo knockout fold: 30.0. in vivo brain penetration, abcb1a-deficient CF-1 mice vs CF-1 wild-type (oral bioavailability also 62.3% vs 7.6%). PMID 21162698 (Taur 2010 Xenobiotica)
  *Notes:* The strongest single efflux number in the table and the drug this whole exercise was motivated by. Caveats: it rests on ONE primary paper (an Eisai in-house study) and the CF-1 model knocks out abcb1a ONLY, not abcb1a/1b. No in vitro ER exists — the Caco-2 work in the same paper used eribulin as an INHIBITOR of digoxin efflux (IC50 > 10 uM), which is inhibitor data and must not be conflated with substrate data. BCRP is coded unknown, not 'no', because the only source for the negative is the FDA Halaven label section 12.3, whose underlying in vitro study is unpublished and carries no PMID. Counterpoint worth carrying: the current label reports no clinically meaningful DDI with ketoconazole or rifampin, so the 30x mouse brain effect does not translate to a clinical IV interaction.

**imatinib** — P-gp `yes`: PMID 12604685 (Dai 2003 JPET): B->A flux 39x greater than A->B in MDCKII-MDR1, reversed by zosuquidar.
  BCRP `yes`: PMID 15251980 (Burger 2004 Blood): established in BCRP-overexpressing lines with mitoxantrone competition and Ko143 reversal. IN VITRO ONLY — see notes.
  In vitro ER: P-gp 39.0, BCRP --, combined 39.00. bidirectional MDCKII-MDR1, Papp(B->A)/Papp(A->B) = 39 (parental MDCKII background 8, so the transfectant-specific component is real but the background is not negligible). PMID 12604685 (Dai 2003 JPET)
  In vivo knockout fold: 28.0. in vivo brain-to-plasma, Mdr1a/1b;Bcrp TRIPLE knockout vs WT. PMID 19225039 (Bihorel/Oberoi 2009 Drug Metab Dispos)
  *Notes:* P-gp is unambiguous and replicated: 11x brain:plasma at 1 h after oral dosing and ~7x after IV in mdr1a/b-/- mice (PMID 12604685), 4x in Mdr1a/1b-/- (PMID 19225039), 5.5x by in situ brain perfusion (PMID 17696988), 6.4x IV (PMID 17380257). BCRP is WEAKER than its reputation and the two windows disagree: in vitro yes, but in vivo Bcrp(-/-) mice show brain:plasma comparable to wild-type (PMID 19225039) and Bcrp1(-/-) uptake matches WT at non-saturating concentrations (PMID 17696988). The 28x figure is the P-gp/Bcrp SYNERGY in the triple knockout and must not be attributed to BCRP alone. NOTE FOR THE GIST VALIDATION: none of this helps there. A scalar per-drug factor cancels exactly out of a tumour:normal ratio, and hepatic P-gp/BCRP sit on the CANALICULAR membrane exporting imatinib into bile, which lowers NORMAL liver and pushes the predicted ratio the wrong way.

**niraparib** — P-gp `yes`: PMID 34730366 (Martins/Loos/Schinkel 2021 Mol Pharm): human ABCB1 transported niraparib moderately in MDCK-II cells; oral elacridar in WT mice reproduced the knockout phenotype.
  BCRP `yes`: PMID 34730366 (in vitro, mouse Abcg2); PMID 41998207 (Macdonald 2026 Br J Cancer, polarised-cell efflux assays). IN VITRO ONLY: single Abcg2-/- mice were indistinguishable from wild-type, so BCRP is not BBB-limiting for niraparib in vivo.
  In vivo knockout fold: 6.5. in vivo brain-to-plasma fold change, Abcb1a/1b-/- and Abcb1a/1b;Abcg2-/- vs WT mice (reported as a 6- to 7-fold range; 6.5 is the midpoint of the authors' own range, not a fitted value). PMID 34730366 (Martins 2021 Mol Pharm)
  *Notes:* No numeric in vitro transwell ER is published for ANY of the three PARP inhibitors in accessible sources — the source papers report transport qualitatively or place ratios in paywalled supplements. Niraparib is the brain-penetrant member of the class: measured brain:plasma AUC 0.29-0.33 vs ~0.03 for olaparib in two xenograft models (PMID 30647846), corroborated by NHP and mouse imaging (PMID 41023116). That ~10x measured separation is a useful external check on this module and it does not pass it — see results/efflux_impact.md.

**rucaparib** — P-gp `yes`: PMID 24962512 (Durmus/Schinkel 2015 Pharm Res); PMID 26438157 (Parrish/Elmquist/Sarkaria 2015 Mol Cancer Ther): intracellular rucaparib reduced to 32.0+/-7.4% of control in MDR1 cells, reversed by zosuquidar.
  BCRP `yes`: PMID 24962512 (very efficiently transported by mouse Abcg2); PMID 26438157 (intracellular rucaparib 14.5+/-5.2% of control in BCRP1 cells, reversed by Ko143).
  In vivo knockout fold: 14.6. steady-state osmotic-pump brain-to-plasma, Mdr1a/1b-/-Bcrp1-/- triple knockout (1.61+/-0.25) vs WT (0.11+/-0.08). PMID 26438157 (Parrish 2015 Mol Cancer Ther), abstract values
  *Notes:* INTERNAL INCONSISTENCY IN THE SOURCE, recorded rather than smoothed: the abstract gives 1.61 vs 0.11 (14.6x), the Results section of the same paper gives 2.0+/-0.37 vs 0.12+/-0.11 (~17x), and the Discussion says 'almost 20-fold'. The abstract value is used because it is the most conservative and the most citable; the true uncertainty on this row is at least 14.6-20x. Two independent groups agree on dual P-gp/BCRP substrate status.

**olaparib** — P-gp `yes`: PMID 35341895 (Song 2022 Eur J Pharm Sci): brain-to-plasma UNBOUND ratio rose 5.6x in Abcb1a/1b-/- mice vs WT; elacridar reproduced the effect.
  BCRP `yes`: PMID 30396915 (Song 2018 JPET): efficiently transported by human ABCG2 in vitro — the component is masked by endogenous canine ABCB1 in MDCKII and only appears once that is suppressed (corrected ER rose 3x under verapamil, 5x under shRNA knockdown). IN VIVO MINOR: single Abcg2-/- was indistinguishable from WT.
  In vivo knockout fold: 8.1. in vivo brain-to-plasma UNBOUND concentration ratio (Kp,uu), Abcb1a/1b;Abcg2-/- vs WT mice. PMID 35341895 (Song 2022 Eur J Pharm Sci)
  *Notes:* Abcb1 dominates (5.6x alone) and Abcg2 adds little in vivo (8.1x combined) despite efficient human ABCG2 transport in vitro — a clean example of why in vitro substrate status does not imply an in vivo barrier effect. Note this quantity is a Kp,uu ratio, i.e. already unbound-corrected, which is exactly the quantity FORM_KPUU models.

## 3. The functional form, and why it is not a knob

Both candidate forms parameterise the same quantity, the efflux burden `b = CL_efflux / P_passive`, and both give `ER = 1 + b`:

- **`kpuu` (default): `f = 1/ER`.** Single-barrier unbound-partition relation — efflux clearance in parallel with passive permeability, which is the derivation behind the published definition `ER = (P_pass + Vmax/Km)/P_pass` itself (PMID 23421687). The `+1` is already inside ER, so dividing by ER does not double-count.
- **`monolayer`: `f = 2/(1+ER)`.** Resolves the cell interior (apical / cell / basolateral, pump at the apical face); exactly 2x more permissive at high ER.

Both return exactly 1.0 at ER = 1, which is why the plain `1/(1+ER)` starting form was rejected: it would silently halve every molecule with no efflux at all.

`kpuu` is the default on three grounds, all pre-committed: it is the form with published empirical support against measured brain `Kp,uu` (r^2 = 0.813 for a high-expressing MDCK-MDR1 line, PMID 24009309; transwell-to-knockout agreement over 32 CNS compounds, PMID 17962372; `Kp,uu` inversely tracking ER above ER ~ 3, PMID 35156331); it does not double-count the `+1`; and it wins the one internal consistency test available. **An in vivo knockout fold change is the same model read differently** — deleting the pump sets `b = 0`, so `f = 1/fold_KO` in EITHER form — which is why using whichever observable a drug happens to have is not form-shopping.

### The consistency test (n = 2 — a tie-breaker, not a calibration)

Two drugs carry both observables. Under `kpuu` an in vitro ER predicts `fold_KO = ER`; under `monolayer`, `(1+ER)/2`.

| drug | in vitro ER | measured KO fold | `kpuu` predicts | `monolayer` predicts |
|---|---|---|---|---|
| adavosertib | 14.1 | 25.0 | 14.1 (1.77x off) | 7.5 (3.31x off) |
| imatinib | 39.0 | 28.0 | 39.0 (0.72x off) | 20.0 (1.40x off) |

Geometric-mean absolute discrepancy: **1.57x for `kpuu`, 2.15x for `monolayer`**. n = 2. Reported because it points the same way as the published correlations, not because two drugs settle anything.

**A scaling factor was available and was refused.** The best-supported practical formulation in the literature is `Kp,uu = 1/(1 + SF·(ER-1))` with an assay- and species-specific `SF` calibrated to neuro-PK data (PMID 30683809; human BBB P-gp is 2.3x LOWER than mouse mdr1a and human BCRP 1.9x higher, PMID 21291474). There is exactly one dataset here we could calibrate `SF` on — the six-drug scorecard this term is about to be tested against — so fitting it would have made the test circular. The shipped term has **zero free parameters** and is correspondingly allowed to fail. It did.

## 4. Scorecard 1 — the six-drug measured human glioma partition

BEFORE is not pasted: it is `glioma_regional_pk.run_scorecard(aggregation="geometric", apply_leak=True)` called live, i.e. the model exactly as it ships. AFTER changes one thing:

```
  before:  P_eff = clip( bbb_penetrance              + L * D_sample , 0, 1 )
  after :  P_eff = clip( bbb_penetrance * f_efflux   + L * D_sample , 0, 1 )
```

Efflux multiplies the TRANSCELLULAR term only — a disrupted tight junction is a paracellular shunt around the pump. (The physically-wrong alternative, multiplying the whole bracket, is run in §7 and changes nothing.)

| drug | measured | modality | struct. penetrance | disrupted-sample frac | **f_efflux** | why | pred before | fold before | pred after | **fold after** |
|---|---|---|---|---|---|---|---|---|---|---|
| eribulin | 0.0051 | intracerebral_microdialysis | 0.0353 | 0.00 | **0.0333** | applied | 0.0353 | 6.944x | 0.0012 | **0.231x** |
| cefazolin | 0.1496 | tissue_homogenate_LCMS | 0.0488 | 0.44 | **1.0000** | substrate_but_no_published_er | 0.1526 | 1.020x | 0.1526 | **1.020x** |
| methotrexate | 0.1780 | intracerebral_microdialysis | 0.0225 | 0.50 | **1.0000** | substrate_but_no_published_er | 0.1415 | 0.795x | 0.1415 | **0.795x** |
| temozolomide | 0.1951 | intracerebral_microdialysis | 0.1254 | 0.14 | **0.6667** | applied | 0.1594 | 0.817x | 0.1176 | **0.603x** |
| levetiracetam | 0.6972 | tissue_homogenate_LCMS | 0.8399 | 0.45 | **1.0000** | unknown_substrate_status | 0.9466 | 1.358x | 0.9466 | **1.358x** |
| adavosertib | 2.4600 | tissue_homogenate_LCMS | 0.4595 | 0.50 | **0.0709** | applied | 0.5785 | 0.235x | 0.1516 | **0.062x** |

| aggregate | before | after | direction |
|---|---|---|---|
| GM \|fold error\| | **1.995x** | **2.424x** | WORSE |
| GM signed bias | 1.066x | 0.460x | |
| within 3x | 4/6 | 4/6 | |
| within 10x | 6/6 | 5/6 | WORSE |
| worst \|fold error\| | 6.94x | 16.23x | WORSE |
| GM \|fold error\|, ex-adavosertib | 1.714x | 1.657x | better |

**3 of 6 drugs move at all.** The other 3 get a factor of exactly 1.0 — not because they are non-substrates, but because the literature has no usable number for them (methotrexate and cefazolin are confirmed BCRP substrates with no applicable ratio; levetiracetam's P-gp status is genuinely disputed).

### Which drug is doing this

| drug dropped | GM \|fold\| before | GM \|fold\| after | efflux helps without it? |
|---|---|---|---|
| eribulin | 1.554x | 2.159x | no |
| cefazolin | 2.281x | 2.881x | no |
| methotrexate | 2.187x | 2.763x | no |
| temozolomide | 2.199x | 2.614x | no |
| levetiracetam | 2.154x | 2.721x | no |
| adavosertib | 1.714x | 1.657x | YES |

The term only helps once adavosertib is removed — and even then it buys 1.714x -> 1.657x, a 3% reduction in error on five drugs. **Dropping the drug that a mechanism gets wrong, and then reporting that the mechanism works, is exactly the move this repo exists not to make.** The ex-adavosertib figure is reported because adavosertib has an independently documented defect (its measured 2.46 is a tissue-homogenate ratio above the model's structural ceiling of 1.0), not because it licenses the headline.

### Per drug, honestly

- **eribulin** (improves): factor 0.0333 from ER 30.0 (in_vivo). 6.94x -> 0.231x.
- **temozolomide** (WORSENS): factor 0.6667 from ER 1.5 (in_vivo). 0.82x -> 0.603x.
- **adavosertib** (WORSENS): factor 0.0709 from ER 14.1 (in_vitro). 0.24x -> 0.062x.

  - *eribulin* is the drug this exercise was for, and efflux does move it the right way — but past the target. It was over-predicted 6.94x; with a 30x mouse knockout factor it is now under-predicted, so the error shrinks without the prediction becoming good, and it still fails the 3x band.
  - *temozolomide* is a genuine surprise from the literature pass: TMZ IS a weak efflux substrate (a good Bcrp1 substrate, PMID 29852323), not a non-substrate as commonly assumed. But the measured combined-knockout effect is only 1.5x, and TMZ was already slightly under-predicted, so even that tiny factor hurts.
  - *adavosertib* is the whole story. See §4.1.

### 4.1 Why adavosertib breaks it — and why that is a finding, not bad luck

Adavosertib has the best-evidenced efflux data in the table: in vitro ER 10.6 (MDCKII-ABCB1) and 4.5 (MDCKII-ABCG2) with proper inhibitor controls and 85-115% mass recovery (PMID 28928160), independently confirmed with a 25x knockout brain fold change by a second group (PMID 29147815). If the efflux term is going to work anywhere, it should work here.

It fails, and **the same paper says why**: Li et al. measured human GBM tumour `Kp,uu` at a **median of 3.2 (range 1.3-24.4)**. Adavosertib does not merely reach human glioma despite dual efflux — it *accumulates* there, because OATP1A2-mediated **uptake** dominates in the acidic tumour microenvironment. The authors state the mouse and human results disagree. Two further details point the same way: the ERs are pH-dependent and roughly halve as basolateral pH falls 7.4 -> 6.5 (i.e. tumour pH suppresses the efflux the model is applying), and the measured value the scorecard grades against (2.46) is itself above 1.0.

So the model over-corrects a drug whose real behaviour in the target tissue is the opposite of the mechanism being added. That is not a tuning failure; it is the efflux term being **the wrong mechanism for this drug in this tissue**.

### 4.2 Is the change even distinguishable at n = 6?

- drugs that move at all: **3 of 6**; of those, 2 worsen and 1 improve
- mean change in `|log10 fold error|`: **+0.0846** (positive = worse)
- exact sign test over the movers: p = 1.000; Wilcoxon signed-rank: p = 0.593

**with 3 of 6 drugs moving at all, neither test can reach significance at any effect size; both are reported to make that explicit.** The honest reading is that six drugs cannot establish that efflux helps OR that it hurts. What they can do is show that it does not deliver the large, obvious improvement that would justify enabling an untested term — and that the one drug with the best efflux data moves the furthest in the wrong direction for a reason the literature already documented.

## 5. Scorecard 2 — GIST tumour:normal, where the term cannot even apply

`results/gist_imatinib.md` diagnosed the missing mechanism as being on the **normal-liver** side: active hepatic uptake (OCT1/OATP) plus lysosomal trapping concentrating imatinib in normal liver. The question here is whether efflux, the tractable half of that diagnosis, moves the prediction toward the measured 0.023-0.142.

Imatinib's factor is **0.0256** (ER 39, in_vitro, PMID 12604685) — one of the strongest in the table. It changes nothing.

### 5.1 The structural result

The prediction is a RATIO of the same drug in two tissues:

```
  R = (f · E_tumour) / (f · E_normal) = E_tumour / E_normal
```

A scalar per-drug factor cancels. Verified to machine precision on all 6 matched pairs: **True**.

| model | geomean fold error | worst | within 3x | predicted liver ratio |
|---|---|---|---|---|
| M1 transport-only (published) | **7.92x** | 25.06x | 1/6 | 0.5796 |
| M1 + scalar efflux, both sides | **7.92x** | 25.06x | 1/6 | 0.5796 |
| efflux on NORMAL only (the real hepatic biology) | **298.19x** | 977.31x | 0/6 | 22.6026 |
| efflux on TUMOUR only (unsupported) | **5.10x** | 43.31x | 3/6 | 0.0149 |

**This is the same class of result the Moncayo explant run produced, reached independently.** There, a spatially uniform ion-trapping retardation `R` was shown to give `C = R·C_bath` everywhere at equilibrium, so it could buy a time constant but never a spatial pattern. Here, a drug-uniform efflux factor `f` gives `f·E_t / f·E_n`, so it can shift a drug's absolute penetrance but never a tissue contrast. Two different mechanisms, the same algebra: **a scalar cannot produce a differential.**

### 5.2 And the tissue-differential version has the wrong sign

Suppose we abandoned the scalar and gave tumour and normal their own factors. Two things follow, and they are both bad.

- **The real hepatic biology moves the prediction AWAY from the measurement.** P-gp and BCRP in liver sit on the CANALICULAR membrane and export imatinib into BILE — i.e. they REMOVE drug from normal hepatocytes. Applying efflux to the normal side only takes the liver prediction from 0.580 to 22.603 against a measured geometric mean of 0.049, and the fold error from 7.92x to 298.19x. Worse, and worse for the right reason: efflux on the normal side is the opposite of the *uptake* the GIST run named.
- **The version that WOULD help is the one the biology contradicts.** Efflux on the tumour side only gives 0.0149 and a fold error of 5.10x — numerically better than the published 7.92x. But it asserts that a GIST liver metastasis effluxes imatinib while the adjacent normal liver does not, and normal liver is among the highest ABCB1/ABCG2-expressing tissues in the body. The required tumour/normal factor ratio to hit the measurement exactly is **0.0853**, i.e. an efflux-ratio asymmetry of **11.7x in favour of the tumour** — in the direction opposite to known transporter expression. Buying a better number by asserting inverted biology is not progress, and we are not doing it.

**Conclusion for GIST: unchanged.** The transport model remains structurally falsified for tissue partition, and efflux is confirmed NOT to be the missing mechanism. The GIST report's own diagnosis stands: what is missing is active hepatic UPTAKE (OCT1/OATP), which is a fundamentally different object — it raises normal-tissue concentration rather than lowering it, and it is tissue-specific rather than drug-scalar.

## 6. An external check the term also fails

Neither scorecard contains a PARP inhibitor, and the following contrast played no part in building the table — so this is the one genuinely out-of-sample test available.

Measured (Sun et al. 2018 Oncotarget, PMID 30647846 (brain:plasma AUC ratio, two xenograft models)): niraparib brain:plasma AUC 0.31 vs olaparib 0.03 — a **10.3x** separation, and one of the better-known facts about this drug class.

| drug | in vivo KO fold | efflux factor | measured brain:plasma |
|---|---|---|---|
| niraparib | 6.5 | 0.1538 | 0.31 |
| olaparib | 8.1 | 0.1235 | 0.03 |
| rucaparib | 14.6 | 0.0685 | -- |

- Ordering correct: **True** (niraparib > olaparib).
- Predicted separation: **1.25x** against a measured 10.3x — the term UNDER-separates by **8.3x**.

So on the one external contrast available, efflux gets the direction right and the magnitude wrong by an order of magnitude — the same pattern as the Moncayo run's compound-specificity test (`ordering yes, magnitude no`). A mechanism that reproduces orderings but not magnitudes cannot be used to produce quantitative predictions, which is what `delivery_source` does.

## 7. Does the verdict depend on anything we chose?

Every discretionary choice in the term, crossed: functional form, which published observable wins when a drug has both, and where the factor multiplies.

| form | prefer | placement | GM \|fold\| before | GM \|fold\| after | 3x | 10x | worst after | improved? |
|---|---|---|---|---|---|---|---|---|
| kpuu | in_vitro | transcellular | 1.995x | 2.424x | 4->4 | 6->5 | 16.23x | no |
| kpuu | in_vitro | total | 1.995x | 3.065x | 4->4 | 6->5 | 59.96x | no |
| kpuu | in_vivo | transcellular | 1.995x | 2.464x | 4->4 | 6->5 | 17.91x | no |
| kpuu | in_vivo | total | 1.995x | 3.372x | 4->4 | 6->5 | 106.31x | no |
| monolayer | in_vitro | transcellular | 1.995x | 2.355x | 4->4 | 6->5 | 13.68x | no |
| monolayer | in_vitro | total | 1.995x | 2.762x | 4->4 | 6->5 | 32.11x | no |
| monolayer | in_vivo | transcellular | 1.995x | 2.464x | 4->4 | 6->5 | 17.91x | no |
| monolayer | in_vivo | total | 1.995x | 3.372x | 4->4 | 6->5 | 106.31x | no |

**0 of 8 configurations improve the scorecard.** The verdict is not an artefact of the form, the source preference, or the placement. The most permissive configuration available (`monolayer` form, which halves the penalty at high ER) still does not recover the shipped model's accuracy.

Two rows of this table are a free correctness check rather than a sensitivity: the two `in_vivo` rows are IDENTICAL across the two forms. That is exactly what the derivation in §3 requires — with `prefer=in_vivo` every moving drug resolves through a knockout fold change, and `f = 1/fold_KO` in both forms — so if those rows had differed, the two forms would not be two windows onto one model.

## 8. Verdict and recommendation

### Recommendation: **default OFF**, and the module ships that way.

`treatment.drugs.delivery_source(..., efflux=False)` is the default and is bit-identical to the previous code (§1). `molecular.efflux` is fully implemented, fully cited, self-tested, and available to any caller who passes `efflux=True` — as a documented, available mechanism, not as an endorsed one.

### What this run establishes

1. **Efflux is implementable from the literature with zero free parameters.** 18 drugs assigned, every claim cited, 7 carrying a usable published ratio. No value was fitted to any dataset in this repo.
2. **It does not improve the one measured scorecard it can act on**, in any of the 8 configurations tested: GM absolute fold error 1.99x -> 2.42x, order-of-magnitude band 6/6 -> 5/6, worst case 6.9x -> 16.2x.
3. **It cannot act on the other scorecard at all**, and the proof is algebraic rather than empirical: a scalar per-drug factor cancels out of a within-drug tissue ratio. This independently reproduces the structural lesson of the Moncayo explant run.
4. **The failure is diagnostic, and it names the same mechanism GIST named.** Adavosertib fails because human glioma UPTAKE (OATP1A2) beats efflux in the measured tissue. GIST fails because normal-liver UPTAKE (OCT1/OATP) is the missing term and efflux has the wrong sign there. Two independent datasets, two different tissues, the same answer: **the missing mechanism is uptake transport, not efflux.** That is a sharper statement than either report could make alone, and it is the useful output of this run.

### What this run does NOT establish

- **Not that efflux is unimportant in vivo.** It plainly is — a 30x mouse brain knockout effect for eribulin is real. What is shown is that a scalar per-drug efflux term does not improve THESE six measured human partition values, and cannot touch the GIST target at all.
- **Nothing at adequate power.** Six drugs, three of which move; two GIST placements that differ by an algebraic identity rather than by data. No p-value here should be read as evidence for anything.
- **Nothing about the drugs with no number.** Half the confirmed substrates in the table get a factor of 1.0 because no usable ratio is published, so this test is silent about them — including doxorubicin, the classic P-gp substrate.
- **Nothing about efflux inhibition, resistance, or intracellular accumulation.** This is the delivery axis only.

### What would change the verdict

A **per-tissue transporter field** — the object that both this run and the GIST run converge on needing, because it is what would make the term non-scalar and therefore able to produce a tissue contrast at all.

This is better placed than the spatial acidic-compartment field the Moncayo run needed, for which no measurement in usable units exists at all. Absolute BBB transporter protein abundance HAS been quantified by targeted proteomics (PMID 21291474: human BBB P-gp 6.06 fmol/ug, BCRP 8.14 fmol/ug), which is a real parameterisation in real units — the thing ion-trapping lacked. But two specific obstacles are worth stating so nobody starts this expecting it to be easy:

1. **Not from the spatial data currently in this repo.** `reconstruct.real_xenium` documents its own panel audit, and ABCB1 and ABCG2 are explicitly on the ABSENT list for the 319-gene Xenium Human Brain panel. The Visium path carries deconvolved cell-type fractions, not a transporter transcript field. So a spatially-resolved transporter map would need a dataset this repo does not currently load.
2. **Transcript is not transport.** Even with expression in hand, converting mRNA to a transport clearance requires exactly the kind of scaling factor refused in §3 — which is fittable only against the same six drugs, and therefore circular for the same reason.

**And efflux is the lower-value of the two candidates.** Both measured datasets point at UPTAKE, not efflux: OATP1A2 in human glioma (which is why adavosertib accumulates where this term predicts it should be excluded) and OCT1/OATP in normal liver (which is what GIST named, and which efflux has the wrong sign for). An uptake term is the better next piece of work, and this run's main contribution is having ruled out the alternative on measured data rather than on intuition.

---

*Generated by `tmesim/validation/efflux_impact.py`. Additive only: `core`, `physics`, `sim.engine`, `molecular.permeability` and `reconstruct.barrier` are untouched; `treatment.drugs.delivery_source` gained one default-off keyword whose bit-exactness when off is verified in §1.*
