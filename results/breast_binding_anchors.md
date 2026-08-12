# Independent anchors for the binding sink — a literature survey, four new candidates, and the null that survives

**Purpose.** `results/breast_tissue_permeability.md` localised the remaining
breast-anchor error to the alpha-scaled loss term `k = clearance + binding`,
and `results/binding_sink_refuted.md` refuted that term in magnitude and in
direction. Both stop at the same place: **`alpha` must not be fit to the breast
anchor**, so what is needed is an *independent measurement* of tissue binding
or clearance. This document is the search for one.

**What this document adds.** It does not repeat
`results/binding_sink_refuted.md`. That file already carries Blakeley's `k30`
values, Mader's 5-FU `N.A.` null, the Grantab/Tunggal head-to-head penetration
comparison, Ekström's rat `Kp,uu ~ 1`, the BBB `Kp,uu` pair, and the Sane
efflux knockout. Everything below is **new to this repo**: four candidate
anchors not previously cited, one re-reading of a source already cited, and a
sharpened statement of the null.

Literature retrieved via PubMed. Entries marked **[verified here]** were read
directly from the PubMed abstract or PMC full text during this survey. Entries
marked **[abstract-only]** had their body text checked for and found
unreachable — every retrieval route from this container is recorded in §8, and
the tier is not upgraded on hope.

**Nothing here was simulated.** This container has no `data_cache/`, so no
engine run was performed. §9 specifies exactly what to run when one is
available, with its prediction written down first.

---

## 1. The quantity being anchored, in the model's units

At steady state under transvascular exchange, `Kp = PS/(PS + k)`, so any
measurement that constrains `k` must be a **first-order loss rate from the free
interstitial pool**, in 1/s. Three numbers frame everything below.

| quantity | 5-fluorouracil | methotrexate | ratio |
|---|---:|---:|---:|
| model `k = clearance + binding` (alpha = 148.608) | 5.586e-3 /s | 2.670e-2 /s | **4.78x** |
| implied interstitial half-life | 2.1 min | 26 s | — |
| `k` required to hit the measured anchor at measured PS | 6.39e-4 – 2.35e-3 /s | 6.67e-4 – 2.45e-3 /s | **1.04x** |

The required band is `k = PS*(1 - Kp)/Kp` evaluated at the measured breast PS
range 1.0e-3 – 3.67e-3 /s (DCE-MRI, `results/breast_tissue_permeability.md` §1)
and the measured anchors `Kp` = 0.61 (5-FU, PMID 9205062) and 0.60 (MTX,
PMID 9679959). **Reproducing both anchors requires a flat `k`. The model has a
4.78x spread.** That is the question every candidate below is being asked to
answer, and only that question.

## 2. Candidate 1 — methotrexate polyglutamate retention in human BREAST cancer cells

**This is the strongest new candidate, and it is breast-specific.** It is also
the exact mechanism the task brief flagged as MTX's real, large, drug-specific
binding term.

Jolivet J, Schilsky RL, Bailey BD, Drake JC, Chabner BA. "Synthesis, retention,
and biological activity of methotrexate polyglutamates in cultured human breast
cancer cells." *J Clin Invest* 1982;70(2):351-60. **PMID 6178756**, PMC371243,
DOI 10.1172/JCI110624. **[verified here]** Companion, same lab, same finding:
Schilsky RL *et al.*, *Adv Exp Med Biol* 1983;163:247-57, **PMID 6193688**,
DOI 10.1007/978-1-4757-5241-0_19 **[verified here]**.

Modality: HPLC/gel-filtration speciation of intracellular MTX and its
γ-polyglutamates, cultured human breast carcinoma lines, 24 h at 2 µM MTX, then
24 h washout into drug-free medium.

| cell line | polyglutamate as % of intracellular drug | polyglutamate amount | intracellular drug retained after 24 h washout |
|---|---:|---:|---:|
| MCF-7 | 55.4 % | 51.9 nmol/g | **22 %** |
| ZR-75-B | 87.6 % | 62.4 nmol/g | **51 %** |
| MDA-231 | 32 % | 4.06 nmol/g | recovers fully; no persistent inhibition at 2 µM |

**Independent?** Yes, completely. Different modality (cell-culture HPLC),
different decade, different laboratory, and no partition or penetration
quantity from this project enters it.

**Converted to the model's units, it does three things.**

1. **It gives a measured release rate for the MTX intracellular pool.**
   Reading the washout as first order, 22 % remaining at 24 h is
   **1.75e-5 /s** (t½ 11.0 h) and 51 % remaining is **7.79e-6 /s** (t½ 24.7 h).
   The model's `k` for MTX is 2.670e-2 /s — **1500x to 3400x faster**.
2. **It settles the compartment, quantitatively.** This rate is the *release*
   rate of an intracellular sink, not the *removal* rate from the interstitium,
   so it cannot be substituted for `k` directly. But it makes
   `results/binding_sink_refuted.md` §2's compartment argument numerical rather
   than rhetorical: MTX's famous retention is a slow intracellular pool with a
   half-life of **11-25 hours**, which is 2-4x longer than the engine's entire
   6-hour simulation window. Over that window it is a nearly static reservoir,
   not a fast sink.
3. **It shows the term is capacity-limited, and the model's is not.** The paper
   states the retained polyglutamate is present "in excess of dihydrofolate
   reductase binding capacity" — i.e. the authors measured the saturable pool
   and found it exceeded. `binding.py`'s term is an unbounded first-order sink
   that never saturates. Whatever the right magnitude, the functional form is
   wrong for this mechanism.

**The finding that most damages a per-drug binding constant.** MDA-231
accumulates **12-15x less** polyglutamate than MCF-7/ZR-75-B (4.06 vs 51.9 and
62.4 nmol/g), and the paper says outright this "could only partially be
explained by decreased ability of these cells to take up free drug." MTX tissue
binding is therefore **not a property of the molecule** — it varies by more
than an order of magnitude across three human breast lines. `binding.py`
computes it from logP and hydrogen-bond-donor count, which are identical in all
three.

**Caveat carried, not dropped.** The mass basis of "nmol/g" could not be
verified: PMC371243 exposes the abstract only, its body being a scanned PDF
with no text layer. Reading it as nmol per ml of cell water gives an
intracellular:medium ratio of ~47x (MCF-7), ~36x (ZR-75-B) and ~6x (MDA-231)
after 24 h, but **that conversion is unverified and must not be used as a
partition coefficient** until the printed table is read.

## 3. Candidate 2 — 5-fluorouracil measured in human BREAST tumour tissue at near steady state

Furuya Y, Yamamoto K, Kohno N, Yamamoto M, Saitoh Y. "Nocturnal 5-fluorouracil
infusion to patients with breast cancer prior to surgery." *Cancer Lett*
1997;116(1):79-84. **PMID 9177461**, DOI 10.1016/S0304-3835(97)00163-8
**[verified here]**.

Modality: HPLC on surgical specimens, n = 13 primary breast cancers, after a
10-hour constant i.v. infusion of 300 mg/m² 5-FU (2100 h to 0700 h, volumetric
pump) — a genuine near-steady-state exposure, not a bolus tail.

| tissue | 5-FU concentration |
|---|---:|
| breast tumour | **25.0 ± 4.1 ng/g** (mean ± SEM; range 6-49), n = 13 |
| normal breast parenchyma | **below detection, < 3 ng/g** |
| adipose | **below detection, < 3 ng/g** |
| ER-positive tumours | 14.4 ± 4.5 ng/g |
| ER-negative tumours | 31.8 ± 5.0 ng/g |

25.0 ng/g is **0.192 µM** at MW 130.08 and unit density.

**Independent?** Yes. Different cohort, different country, different modality
(tissue HPLC, not microdialysis) from the Müller anchors.

**What it does and does not give.** It is the only measured *absolute* 5-FU
concentration in human breast tumour this survey found, and it is at steady
state, which is the regime `Kp = PS/(PS + k)` describes. **It cannot be turned
into a `Kp` here**, because the paper reports no matched plasma concentration
and this survey will not supply one from an assumed clearance — that would
manufacture the number the repo's convention forbids. Recorded as a measured
tissue level and as the highest-value *future* anchor: one plasma 5-FU
measurement during an identical infusion converts it into a second, fully
independent breast `Kp`, by a modality that does not share microdialysis's
recovery-calibration assumptions.

**A second, weaker use, available now.** Tumour ≥ 8x normal breast parenchyma
and adipose (25 vs < 3 ng/g) is a measured within-patient contrast that any
corrected model must still produce. The engine currently has no mechanism that
would distinguish them other than vascular density.

**The ER split is a warning, recorded rather than used.** A 2.2x difference
between ER-positive and ER-negative tumours in the same study, at the same
dose, means breast `k` (or breast PS) carries real biological structure that a
single per-drug constant cannot represent. It is the tissue-level analogue of
the MDA-231 result in §2.

## 4. Candidate 3 — methotrexate interstitial exposure equals FREE plasma, re-read

`results/binding_sink_refuted.md` §1 uses Blakeley JO *et al.*, *J Neurooncol*
2009;91(1):51-8, **PMID 18787762**, PMC4029102, DOI 10.1007/s11060-008-9678-2
for its `k30` values. The full text carries a second quantity that file does not
use **[verified here, PMC full text]**:

> "Considering that MTX is 54 %-59 % protein bound in human plasma it appears
> that interstitial concentrations achieved in regions of BBB disruption are
> effectively limited by the free fraction of drug in systemic blood."

With AUC(ECF)/AUC(plasma) = 0.28-0.31 in contrast-enhancing tumour and
`fu,plasma` = 0.41-0.46, the **unbound** partition is

    Kp,uu = 0.28/0.46 to 0.31/0.41  =  0.61 to 0.76

and therefore `k/PS = 0.32 to 0.64`, i.e. `PS/(PS+k) = 0.61-0.76`.

**Compare the breast anchor: 0.60.** Two tumours, two organs, two laboratories,
two decades apart, same drug, same answer to within the measurement spread.
The model puts MTX at 0.0523-0.0571.

Also from the same full text, and not previously recorded here: peak ECF
concentration was **27-30 %** of peak plasma in enhancing tumour, ECF
concentrations exceeded 2 µM for **20-26 hours** in both enhancing and
non-enhancing regions, and in enhancing tumour the ECF terminal slope was
**plasma-limited** ("similar to the Cp-time profiles"), which is why that half
of Blakeley's `k30` spread is not a pure tissue constant.

**Independent?** Yes of the breast anchors; **not** of `binding_sink_refuted.md`
§1, which is the same paper. Counted once, as a re-reading, not as a fifth
source.

**Caveat that must travel with it.** This is brain. Importing a brain number
into breast is the precise error `results/breast_kp_validation.md` §3.1 exists
to record. It is used here only as a *ratio* (`Kp,uu`), which is the form least
contaminated by barrier physiology, and only as corroboration of direction.

## 5. Candidate 4 — measured Kp for BOTH drugs in the SAME mice (lead, abstract-tier)

Poulin P, Chen YH, Ding X, Gould SE, Hop CECA, Messick K, Oeh J, Liederer BM.
"Prediction of drug distribution in subcutaneous xenografts of human tumor cell
lines and healthy tissues in mouse." *J Pharm Sci* 2015;104(4):1508-21.
**PMID 25615572**, DOI 10.1002/jps.24336. **[abstract-only]**

`results/binding_sink_refuted.md` §6 already names this as the single
highest-value retrieval action. This survey confirms it is closed from here
(§8) and adds **one statement that is usable from the abstract alone**:

> the tissue-composition model — which predicts `Kp` from **lipid binding and
> plasma-protein binding only**, with no drug-specific sink — landed within
> **twofold of measured in 77 %** and **threefold in 93 %** of cases across six
> antineoplastics (including MTX and 5-FU) in three human tumour xenografts and
> four healthy tissues; the abstract names its failures explicitly — **brain
> for every drug, liver for MTX, and every xenograft for gemcitabine.**

**Neither MTX nor 5-FU in tumour appears in that exception list.** So a
mechanism containing no drug-specific binding sink reproduced both drugs'
measured tumour partitioning to within 2-3x. That is evidence against a large
drug-to-drug sink term, obtained without reading the tables — but it bounds the
*model class*, not the ratio, and it is a mouse subcutaneous xenograft. It does
not close the null in §7.

## 6. Candidate 5 — mouse BREAST tumour methotrexate microdialysis with fitted rate constants (lead, closed)

Alsulimani HH, Kim J, Sani SN. "Microdialysis-directed intra-tumor
pharmacokinetic modeling of methotrexate in mice and humans." *J Pharm Pharm
Sci* 2016;19(2):239-51. **PMID 27518173**, DOI 10.18433/jpps.v19i2.26074.
**[abstract-only]**

From the abstract **[verified here]**: microdialysis of free MTX in the tumours
of **breast tumour-bearing mice**, fitted to a hybrid PBPK model with tumour
blood flow and volume and plasma as a forcing function, which "produced various
transfer rate constants between plasma and tumor", then scaled to predict human
tumour concentrations.

**This is the closest published match to what the model needs that this survey
found**, and `results/binding_sink_refuted.md` does not cite it: a *breast*
tumour, *free* (microdialysate) drug, and fitted first-order plasma↔tumour
transfer constants — the same compartmental structure as `k13`/`k30`, in the
right tissue. Its numeric parameters are behind the University of Alberta
journal host, which is egress-blocked from this container (§8).

**Ranked action: retrieve this alongside Poulin.** Poulin gives the ratio across
both drugs in one experiment; Alsulimani gives an absolute rate for MTX in
breast. Together they would answer §1's question outright.

## 7. Which way the evidence points

**Away from 4.78x, and toward flat or inverted — but the direct measurement
still does not exist.**

Pointing away from 4.78x (i.e. supporting a flat `k`):

* §2. MTX's extra retention is measured to be **intracellular**, and slow —
  release t½ 11-25 h, 1500-3400x slower than the model's MTX `k`. It is in the
  wrong compartment and on the wrong timescale to be a fast interstitial sink.
* §2. MTX binding varies **12-15x across three human breast lines** with
  identical logP and HBD count, so it is not a molecular constant at all.
* §4. MTX's measured tumour `Kp,uu` is **0.61-0.76**, matching the 0.60 breast
  anchor and 5-FU's 0.61. Near-equilibration with free plasma, for both drugs.
* §5. A no-sink partition model reproduced **both** drugs' measured tumour `Kp`
  to within 2-3x, and neither is listed among its failures.

Pointing toward a spread:

* **Nothing found in this survey.** The one real mechanism that would justify
  MTX > 5-FU — polyglutamation — was found, measured, in breast, and it turns
  out to sit in the wrong compartment (§2). This was the task brief's leading
  candidate for a drug-specific term and it did not survive contact with the
  measurement.

**The null that survives, and it is the same one
`results/binding_sink_refuted.md` §5 reported.** No study has measured a
first-order *interstitial* loss rate for both drugs in the same tissue by the
same modality. Methotrexate has one (Blakeley, brain, n = 4, spanning 350x);
5-fluorouracil has none (Mader 2003 reports `t½` = N.A. in tumour, subcutaneous
tissue and plasma alike). **A measured rate ratio does not exist.** This survey
did not find one and does not claim to have closed it. What it changes is the
weight of the surrounding evidence: four new lines, three of them breast, all
pointing the same way, and the strongest remaining argument *for* a
drug-specific MTX sink now measured and disqualified on compartment grounds.

**Stated as a claim with its own limits:** the evidence supports a roughly flat
`k` across these two drugs *strongly enough that the 4.78x spread should be
treated as an artefact of `binding.py`'s structural heuristic*, and *not
strongly enough to assign a value to `k`.* The value must still come from §5 or
§6, or from the experiment in §9.

## 8. A consistency check that is NOT a new constraint

Worth recording because it looks like a contradiction and is not.

The 19F-MRS literature reports a tumour pool of free 5-FU with a half-life far
longer than blood: **0.4-2.1 h** in human tumours (Wolf W *et al.*, *PNAS*
1990;87(1):492-6, **PMID 2296605**, PMC53290, DOI 10.1073/pnas.87.1.492
**[verified here]**); **0.33-1.3 h** in 6 of 11 tumours including breast
(Presant CA *et al.*, *J Clin Oncol* 1990;8(11):1868-73, **PMID 2230874**, DOI
10.1200/JCO.1990.8.11.1868 **[verified here]**); **51 of 99** patients with
breast, colorectal and other tumours above a 20-minute threshold (Wolf W *et
al.*, *NMR Biomed* 1998;11(7):380-7, **PMID 9859944** **[verified here]**);
**19 of 57** above it (Presant CA *et al.*, *Lancet* 1994;343(8907):1184-7,
**PMID 7909867**, DOI 10.1016/S0140-6736(94)92399-X **[verified here]**); blood
t½ 5-15 min.

Naively this is an upper bound on `PS + k` of **5.78e-4 /s** at the 20-minute
threshold, down to **9.17e-5 /s** at 2.1 h — which would put the model's 5-FU
`k` alone 10x to 61x too fast. **Do not use it that way.**
`results/binding_sink_refuted.md` §2.1 already disqualifies the substitution:
19F-MRS measures *total* tissue drug, and for the binding component it carries
the wrong sign.

The correct reading is a consistency check, and it passes. At the measured DCE
PS of 1.0e-3 – 3.67e-3 /s, a *freely exchanging* interstitial pool
re-equilibrates with plasma on a timescale of `ln2/(PS/Kp)` = **1.9 to 7.0
minutes**. Any 19F signal persisting 20-126 minutes therefore *cannot* be the
freely exchanging pool — it must be intracellular. Which is exactly what §2
found independently for the other drug. The two measurements agree about where
the retained drug is; neither tells us how fast it leaves the interstitium.

## 9. The exact experiment to run when `data_cache` is available

No engine run was performed for this document. Three runs, in this order. Each
prediction below is **written down before the run**, per the convention
`results/breast_tissue_permeability.md` established.

**E1 — the flat-`k` counterfactual, zero free parameters.**
Driver: a new `experiments/breast_binding_ratio.py`, modelled on
`experiments/breast_tissue_permeability.py` (override compiled `DrugParams`;
touch nothing under `tmesim/`). Substitute the **experimental** logP for both
drugs — MTX −1.85, 5-FU −0.89 (Hansch/Leo/Hoekman via PubChem CID 126941 and
3385, already recorded in `results/binding_sink_refuted.md` §3) — recompute
`binding.tissue_binding_rate`, keep `alpha = 148.608` **untouched**, keep PS at
the measured breast value, exchange ON, same section, `dx = 12 um`, 6 h,
`dt = 10 s`, seed 0. This is a measured-input substitution, **not** a fit to the
anchor.
*Prediction:* the model's own `k` ratio collapses from 4.83x to 1.51x, the two
drugs' predicted `Kp` converge to within ~1.5x — and **neither lands in its
band**, because the *level* of `alpha`, not its drug-to-drug spread, is still
10-40x too high. If both drugs land in band, this prediction is wrong and the
spread was the whole error.

**E2 — the washout test, an out-of-sample falsification of `k` that never
touches the breast `Kp` anchor.**
Run 5-FU and MTX to quasi-steady state, then set `plasma_conc = 0` and fit the
decay of `mean_tumor_auc`'s underlying mean interstitial concentration.
*Prediction:* t½ = 2.1 min (5-FU) and 26 s (MTX), reproducing `ln2/k`
analytically inside the engine, against Blakeley's measured MTX ECF terminal
t½ of 3.15-14.51 h. This scores the magnitude claim *in the solver* rather than
on paper, and it is the run that would catch a silent no-op — instrument it
with `tmesim/validation/propagation.py` (commit 639358d), because a driver
that calls `simulate_drug_transport` directly and forgets `exchange_k` returns
a bit-identical, entirely plausible wrong answer.

**E3 — the requirement ladder, reported as a requirement and NOT written back.**
Sweep `k` for both drugs over [1e-4, 3e-2] /s at fixed measured PS and report
which `k` puts each drug in its band. *Prediction:* a **single flat** `k` in
6.4e-4 – 2.4e-3 /s puts **both** drugs in band simultaneously — the analytic
solve in §1 says it must, and E3 tests whether the full spatial solve agrees
with the analytic steady state on real tissue.
**The resulting value must not be written into `binding.py` or into `alpha`.**
That would be fitting a third parameter to the two data points that diagnosed
the first two, which is the thing this whole line of work has refused to do
four times. E3 is a diagnostic. The value of `k` still has to come from §5 or
§6.

## 10. What would falsify the conclusion in §7

* **Poulin (§5) or Alsulimani (§6), read in full, showing a measured
  tumour `Kp` or transfer-rate ratio near 4.78x for MTX:5-FU.** That single
  table would overturn §7 outright. It is why both are ranked as the top
  retrieval actions rather than treated as settled.
* **A measured interstitial elimination rate for 5-FU in any solid tumour**,
  showing 5-FU's `k` to be 4-5x *below* MTX's. Mader 2003's `N.A.` says this may
  not be achievable for 5-FU by microdialysis at clinical doses; a different
  modality would be needed.
* **Evidence that MTX polyglutamation depletes the interstitium fast enough to
  act as a 6-hour sink** — i.e. a measured *uptake* clearance, not the *release*
  rate §2 supplies. This is the one measurement that would rescue the model's
  direction, and this survey did not find it for either drug.
* **E1 landing both drugs in band.** That would mean the spread, not the level,
  was the entire error, and §7's "treat the level separately" framing is wrong.

## 11. Honest nulls — what was searched and returned nothing

Searched via PubMed and returned zero usable results:

* `fu,tissue` / fraction unbound in tissue homogenate, by equilibrium dialysis
  or any other method, for **either** drug, in **any** tissue, in **any**
  species. This confirms and does not merely repeat
  `results/binding_sink_refuted.md` §5 — it was re-searched independently here
  with the same outcome. **There is no measured tissue-binding constant for
  5-fluorouracil or methotrexate.** This is the single biggest gap.
* Measured intratumoural methotrexate elimination half-life from serial biopsy
  in any human solid tumour. The MTX intratumoural literature that exists is
  brain (PMID 18787762 microdialysis; PMID 15170895 a single cystic ependymoma
  with AUC(tumour)/AUC(serum) = 1.95; PMID 1652112 intratumoural infusion in
  glioblastoma) — all CNS, none breast, none a clean tissue washout.
* 5-fluorouracil tumour microdialysis with a resolvable elimination phase in
  any animal model.
* Methotrexate polyglutamate accumulation measured in **human breast tumour
  tissue** (as opposed to cell lines, §2, or leukaemic blasts — PMID 30849662
  measures MTX-PG accumulation against BCRP/MRP4 expression in childhood ALL,
  which is the wrong tissue and the wrong compartment for this purpose).
* A primary measurement of 5-fluorouracil plasma protein binding. The ~8-12 %
  figure is widely quoted; this survey did not locate a primary source for it
  and has therefore **not** used it to compute a 5-FU `Kp,uu`, which is why §4
  reports that quantity for methotrexate only.

Retrieval routes attempted and blocked by the container's egress proxy, so that
the trail is not re-walked: `journals.library.ualberta.ca` (Alsulimani full
text), `onlinelibrary.wiley.com` and `sciencedirect.com` (Poulin full text),
`europepmc.org` and `www.ebi.ac.uk` (both papers, and Blakeley's Table 2),
`pmc.ncbi.nlm.nih.gov` over HTTP. PMC content reachable **only** through the
PubMed MCP `get_full_text_article` tool, which returns body text but **strips
tables** — this is why Blakeley's `k30` values in
`results/binding_sink_refuted.md` §1 could be confirmed in prose but not
re-derived here, and why PMC371243 (§2) yielded abstract only.

## 12. Status

**The binding sink still has no independent measured value, and the breast
transport claim remains UNVALIDATED.** What has changed is that the search
space is now mapped: the mechanism most likely to justify the model's 4.78x has
been found, measured in breast, and disqualified on compartment and timescale
grounds; three further measured lines point the other way; and the two
documents that would settle it numerically are identified, ranked, and blocked
only by access.

`alpha` was not touched. No parameter was changed. No number in this document
was fitted to anything.
