# The two blockers that kept the two-route model default-off

> # ⚠️ CORRECTED 2026-08-14 — THE BLOCKERS WERE DISCHARGED AND THE DEFAULT WAS FLIPPED. THE TITLE NOW DESCRIBES A CLOSED STATE.
>
> **This banner is the only edit made to this file.** No number, table or verdict
> below it has been changed or re-run. Everything below was true on
> **2026-08-13**.
>
> **What changed.** This file produced the evidence for a flip decision and
> explicitly *"does not flip anything."* Commit `fd2a9e7` (2026-08-14, *"Flip
> two-route transport ON"*) took that decision. Record:
> `results/two_route_promotion.md`.
>
> **Specifically false now:** the header's *"Nothing is wired into
> `compile_drug.py`, `DrugParams.diffusivity` is bit-identical, and
> `calibrate_alpha()` still returns **148.60815765445545**."* All three clauses
> are false of the shipping engine; `alpha` is **22.223477711172563**. The title
> — *"the two blockers that kept the two-route model default-off"* — should be
> read in the past tense.
>
> **What still stands, and is the reason this file matters.** The efflux
> composition is an algebraic identity with no free parameter; the H-bond-donor
> decrement is taken from published Abraham coefficients and **scores worse on R²
> than the value it replaces**, and the sourced value was adopted anyway. That
> discipline is unchanged and was re-verified at the promotion.
>
> **Two qualifications added later, which this file could not have carried.**
> 1. `results/descriptor_falsification.md` finds the published blind R² was
>    computed with the **paper's own measured descriptors**, not the ones the
>    shipping pipeline computes. Every R² below — 0.718, 0.725, and the
>    0.588–0.722 sweep range — inherits that qualification.
> 2. The H-bond-donor count is the variable this file's §3 sweeps. A units
>    question has since been raised about which donor count the engine feeds that
>    law at runtime. It is **under investigation and unresolved**; nothing here is
>    corrected for it. Do not treat §3.3's robustness range as settled until it
>    closes.

**What this file is.** `results/two_route_transport.md` added a transcellular
route that moved the blind Pruijn score from R² = 0.022 to 0.725, and then
recommended against flipping the default for two stated reasons: an unresolved
efflux question centred on paclitaxel (§5.1 there), and one parameter with no
primary source sitting at its own sweep optimum (§4.4 there). This file resolves
both, and produces the evidence for a flip decision. **It does not flip anything.**

Driver: `experiments/two_route_efflux.py`. Code changed: one additive function
and one default-`None` keyword in `tmesim/molecular/two_route.py`, plus nine
tests in `tests/test_two_route.py`. `tmesim/molecular/efflux.py` is **unchanged**
— see §2.4 for why, which is itself a finding. Nothing is wired into
`compile_drug.py`, `DrugParams.diffusivity` is bit-identical, and
`calibrate_alpha()` still returns **148.60815765445545**.

> **The rule this file was written under, unchanged from its parent.** Nothing is
> fitted to Pruijn. The efflux composition is an algebraic identity with **no free
> parameter** (§2.1, verified to 3.7e-16). The H-bond-donor decrement is now taken
> from published Abraham coefficients and **scores worse on R² than the value it
> replaces** (§3.3) — which is the point. The pre-registered
> `results/sweep/tworoute_prereg.json` is untouched and its checking stage still
> passes.

---

## 0. Bottom line

**Blocker 1 — efflux. The hypothesis was right about the mechanism and wrong
about the fix.**

* Efflux composes with the two-route law through an **exact identity**, not a
  modelling choice: a cell with a pump has its *interior* resistance multiplied
  by the efflux ratio and its *membrane* resistance untouched, so efflux enters
  as `K → K/ER` at a parameter the model already has. **Zero new free
  parameters.** Verified against a direct slab solve to 3.7e-16 over 420 points.
* The composition confirms the hypothesis exactly: at ER = 39 paclitaxel's
  transcellular gain falls **4.94× → 1.14×**, while **doxorubicin and epirubicin
  do not move at any efflux ratio, including 10⁴**, because their membranes are
  shut. Under the previous transport model the term was identically zero for
  *every* drug. That is the precise sense in which efflux had nothing to act on.
* **But it is not scoreable against Pruijn, for three independent sourced
  reasons (§2.2), and the strongest is new: HT-29 is functionally
  P-glycoprotein-negative.**
* **And it cannot be applied to paclitaxel anyway.** No primary numeric efflux
  ratio for paclitaxel could be retrieved (§2.4). Under the efflux module's own
  house rule — never invent a ratio — paclitaxel's row would carry no number, so
  the term is 1.0 for the one drug it was needed for.

**Blocker 2 — the H-bond-donor decrement. Resolved by sourcing it.**

| `hbd_decrement` | provenance | R² | Spearman | geo bias | within 3× |
|---:|---|---:|---:|---:|---:|
| 0.000 | parameter-free, nothing unsourced | 0.345 | +0.703 | 2.20× | 66 % |
| 1.000 | the pre-registered **convention** | **0.725** | +0.842 | 1.25× | 96 % |
| **1.340** | **SOURCED** — Abraham octanol/hexadecane | **0.718** | +0.821 | **1.03×** | **99 %** |

The sourced value is **not** the sweep optimum, scores marginally worse on R²
(0.718 vs 0.725), and is **better** on geometric bias (1.03× vs 1.25×) and
within-3× coverage (99 % vs 96 %). Across the entire measured range of the one
remaining choice it never falls below R² = 0.588.

**The two blockers turn out to be the same blocker.** At the sourced decrement
the paclitaxel outlier that motivated the efflux question **largely disappears
without any efflux term**: the panel differential falls **4.94× → 1.60×** and
paclitaxel's transcellular gain **4.94× → 1.60×**. The 4.94× was substantially an
artefact of an under-sourced donor penalty, not evidence of missing efflux.

**Blast radius: exactly 1 of 12 verdicts moves — docetaxel FAIL → PASS, in
`frozen` only. In `recal`, the mode the codebase actually runs, none moves.** The
same pattern as the parent. The sourced decrement changes the blast radius'
magnitude, not its verdicts (§4.3).

**Recommendation: §5. Short version — yes, flip, at 1.34, in `recal`, efflux off,
but add the retardation term first.**

---

## 1. What was asked, and what each answer rests on

| question | answer | rests on |
|---|---|---|
| Is efflux commensurable with Pruijn? | **No** | 3 independent reasons, §2.2 |
| Does two-route + efflux fix paclitaxel? | mechanism yes, **application no** | §2.1, §2.4 |
| What does it do to the Pruijn score? | **exactly nothing** — 0 of 71 compounds carry a ratio | §2.2 |
| How was the HBD parameter resolved? | **sourced** at 1.34, scores slightly worse | §3 |
| Flip? | see §5 | §4 |

---

## 2. Blocker 1 — efflux

### 2.1 How efflux composes, and why it is an identity rather than a choice

The question is not whether to multiply something by `1/ER`. It is **where the
pump acts**. Take one cell as a slab crossed by the transport direction: an entry
membrane with passive permeability `P` and outward active clearance `E`, an
interior of conductivity `σ_c` over path length `L`, and an exit membrane with
the same `P` and `E`. The three steady-state flux statements are

```
J = P*(C1 - Ca) - E*Ca          [entry face: passive in, pump out]
J = (σ_c/L)*(Ca - Cb)           [the interior]
J = P*(Cb - C2) + E*Cb          [exit face]
```

Eliminating `Ca` and `Cb` gives, in closed form,

```
(C1 - C2)/J  =  2/P  +  (L/σ_c) * ER,        ER = 1 + E/P
```

**Efflux multiplies the interior resistance by the efflux ratio and leaves the
membrane resistance untouched.** In this module's series form
`1/σ_c* = 1/σ_c + 1/(a·P)` that is exactly `σ_c → σ_c/ER`, and since
`σ_c = K·D_cyto` it is exactly **`K → K/ER` with `P` unchanged**. Efflux
therefore enters at a parameter the model already has and introduces **no new
degree of freedom**.

> `--stage algebra` solves the slab directly at 420 `(P, σ_c, L, E)` points and
> checks all three flux statements and the closed form. Max relative error
> **3.7e-16**. Asserted again in
> `tests/test_two_route.py::TestEffluxComposition`.

**This is the same object `efflux.FORM_KPUU` already models.** That module's
default form is `f = Kp,uu = P_pass/(P_pass + CL_eff) = 1/ER` — an *unbound
partition ratio across one membrane at steady state*. `K` here is also an unbound
partition ratio across one membrane (cytoplasm over interstitium). Composing them
is multiplying a partition by a partition. The derivation above is
compartment-agnostic.

**The falsifiable consequence, and it is the interesting one.** The interior
resistance is what efflux multiplies, so where the membrane is *shut* the
interior contributes nothing to the series and efflux changes the tissue
diffusivity by **exactly nothing**:

| probe | membrane resistance | gain at ER = 1 | gain at ER = 10⁴ | moved? |
|---|---:|---:|---:|:--:|
| doxorubicin-like (logP 0.0, HBD 6) | 1.99e5 | 1.0000 | 1.0000 | **no** |
| 5-FU-like (logP −0.8, HBD 2) | 125.8 | 1.0478 | 1.0006 | yes |
| cisplatin-like (logP 0.2, HBD 2) | 12.58 | 1.3368 | 1.0005 | yes |
| docetaxel-like (logP 3.26, HBD 5) | 10.95 | 1.2745 | 1.0003 | yes |
| **paclitaxel-like (logP 3.74, HBD 4)** | **0.363** | **3.3055** | **1.0003** | **yes** |

Under the previous transport model there were no cell inclusions at all, so `K`
never entered the calculation and this term was **identically zero for every
drug**. That is the precise sense in which efflux had nothing to act on until a
transcellular route existed — and it is a consequence of the algebra, not a
rationalisation of one. **The hypothesis in the brief was correct.**

### 2.2 Is efflux commensurable with the Pruijn dataset? No — three reasons

Asked and answered *before* scoring anything, because scoring first and asking
afterwards is how a dataset gets quietly fitted.

**(i) Coverage is zero.** `--stage commensurability`: **0 of the 71** Pruijn
compounds carries any assignment in `efflux.EFFLUX_TABLE`. They are 67
benzotriazine di-N-oxide analogues of tirapazamine from the Auckland Cancer
Society Research Centre plus NLCQ-1, urea, mannitol and acridan-9(10H)-one. A
literature pass found **no published efflux ratio for tirapazamine, SR-4233 or
any SN-numbered compound**; `tirapazamine + P-glycoprotein` returns exactly one
PubMed record and it is a medicinal-chemistry paper using tirapazamine only as a
cytotoxicity reference. So `efflux_factor()` returns exactly 1.0 for every
compound in the dataset and applying it would change the score **by exactly
zero** — not by a little, by zero. There is no test here. This fact is pinned by
`test_none_of_the_pruijn_compounds_carries_an_efflux_ratio` so it cannot rot
silently.

**(ii) HT-29 is functionally P-glycoprotein-negative.** This is the load-bearing
one and it is new. Pruijn's measurement is diffusion through an **HT-29**
multicellular layer. Three independent primary sources:

> **Wils P, Phung-Ba V, Warnery A, Lechardeur D, Raeissi S, Hidalgo IJ, Scherman
> D.** *Biochem Pharmacol* 1994;48(7):1528–1530, PMID **7945455**,
> [DOI](https://doi.org/10.1016/0006-2952(94)90580-0) — *"No functional expression
> of P-glycoprotein was found in the differentiated HT29-18-C1 cell line"*,
> measured in the same experiment in which Caco-2 P-gp was *"very high"*.

> **Lorke DE et al.** *J Nucl Med* 2001;42(4):646–654, PMID **11337555** —
> parental HT-29 cells *"completely devoid of JSB-1 staining"* (anti-P-gp);
> rhodamine 123 *"accumulated strongly"* in them.

> **Rimet O, Mirrione A, Barra Y.** *Biochem Biophys Res Commun* 1999;259(1):43–49,
> PMID **10334913**, [DOI](https://doi.org/10.1006/bbrc.1999.0711) — HT29-D4
> *"constitutively expresses a very low level of the MDR1 gene product"*; MDR1 had
> to be transfected to obtain resistance.

Every numeric ratio in `EFFLUX_TABLE` comes from a **transfected P-gp
overexpressor** (MDCKII-ABCB1, MDCKII-MDR1) or a **transporter-knockout mouse**.
Applying such a number to a P-gp-negative calibration tissue would charge Pruijn
for a transporter it does not have.

*The mirror-image caveat, recorded because it cuts the other way:* HT-29 **does**
endogenously express functional **BCRP/ABCG2** and MRP1 — Yanase 2004 (PMID
15367706, *"HT-29 cells, which endogenously express BCRP"*), Jendzelovský 2009
(PMID 20024169), Perego 2001 (PMID 11507048, BCRP is the axis HT-29 amplifies
under selection, *"independent of expression of P-glycoprotein or MRP"*). So a
**BCRP** term would be commensurable with Pruijn in principle. None of the 71
compounds has a BCRP ratio either, so this remains unactionable — but it is the
one efflux route that is not ruled out by the tissue, and it is currently
unmodelled. *No primary source quantifying HT-29 ABCB1 or ABCG2 against Caco-2 or
MDCKII-MDR1 on a common scale was found; marked UNVERIFIED.*

**(iii) The sign is not even settled, and Pruijn's own authors say so.** From the
paper's Conclusions, verbatim: *"Active transporters such as ABC pumps are known
to modify extravascular drug transport in tissue, and can enhance penetration of
MCLs by reducing sequestration in cells and increasing paracellular transport."*
They argue efflux would **increase** MCL penetration. That is the opposite sign
to §2.1's steady-state conductance result — and both can be right, because
Pruijn's measurement is a **transient** and reducing intracellular *capacity*
speeds a transient front while reducing *conductance* slows it. The two-route
model is a steady-state conductivity model with **no capacity term**, which
`results/two_route_transport.md` §6 already lists as known-missing physics (item
1, the `K > 1` retardation basis). So the model cannot represent the effect
Pruijn's authors describe, and scoring a one-sided efflux term against their
transient would be scoring the wrong sign.

> The paper's own composition, confirmed verbatim from the publisher page:
> *"67 structurally diverse benzotriazine di-N-oxides (analogues of the
> hypoxia-activated prodrug tirapazamine) plus four miscellaneous compounds"*.
> Its predictive model uses **MW, logP₇.₄, H-bond donors and H-bond acceptors
> only**, at R² = 0.907, with **no transporter term of any kind**.

**Conclusion: efflux is not commensurable with Pruijn, and it was therefore
evaluated on the tissue panel instead**, as the brief directed.

### 2.3 On the tissue panel: the mechanism works

`--stage panel`, breast panel, ER swept because the published ratios come from
transfected overexpressors and the claim that a breast tumour cell pumps that
hard is not in evidence. **Transcellular gain:**

| ER | doxorubicin | paclitaxel | docetaxel | epirubicin | 5-FU | cisplatin | panel differential |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1.0000 | **4.9395** | 1.4646 | 1.0000 | 1.0393 | 1.3738 | **4.939×** |
| 3 | 1.0000 | 2.6470 | 1.3983 | 1.0000 | 1.0386 | 1.3264 | 2.647× |
| 10 | 1.0000 | 1.5424 | 1.2657 | 1.0000 | 1.0366 | 1.2262 | 1.542× |
| 20 | 1.0000 | 1.2770 | 1.1801 | 1.0000 | 1.0341 | 1.1572 | 1.277× |
| **39** | 1.0000 | **1.1435** | 1.1117 | 1.0000 | 1.0302 | 1.0995 | **1.144×** |
| 100 | 1.0000 | 1.0564 | 1.0503 | 1.0000 | 1.0220 | 1.0457 | 1.056× |

**Read the doxorubicin and epirubicin columns.** They are 1.0000 at *every*
efflux ratio including 100. Both are confirmed P-gp substrates. The term cannot
move them because their membranes are shut, exactly as §2.1 predicted. This is
the sharpest available demonstration that the composition is a physical statement
and not a knob: it declines to act on the two drugs a knob would have been most
tempted to act on, one of which is the α calibration anchor.

So the answer to *"does two-route + efflux reduce paclitaxel's gain toward
something defensible?"* is **yes, mechanistically** — 4.94× → 1.14× at ER = 39.

### 2.4 And it still cannot be applied — the literature will not supply the number

`efflux.py`'s house rule is that a ratio is recorded only when a specific number
can be attributed to a specific paper. Against that rule:

| drug | P-gp | BCRP | usable numeric ratio? | source |
|---|---|---|---|---|
| **paclitaxel** | **yes** | unknown | **NO** | see below |
| docetaxel | yes | unknown | **yes, ER ≈ 20** | Wils 1994, PMID 7945455 |
| epirubicin | yes | unknown | no | Lo & Huang 2000, PMID 10677583 — direction only |
| 5-fluorouracil | **unknown** — no positive evidence found | unknown | no | — |
| doxorubicin | yes | no | no (already recorded) | — |
| cisplatin | no | no | n/a, measured non-substrate | — |

**Paclitaxel, the drug the whole blocker is about, has no usable number.**

* *P-gp substrate status is not in doubt* — Sparreboom 1997 (PMID 9050899),
  Takeuchi 2006 (PMID 16779700), Taub 2005 (PMID 16093365).
* *In vitro monolayer ER:* every candidate is behind a paywalled table or figure
  (Takeuchi 2006, Taub 2005, Hugger 2002 PMID 12210046, Mitra 2011 PMID
  21538355). The one open-access source, Cai 2022 (PMID 35509628, PMC9060980),
  reports ERs only *under inhibition* and never states the paclitaxel-alone
  control ER in text; its reported `Papp` values are also implausibly large for
  paclitaxel. **Unusable.** Mitra 2011 does establish verbatim that paclitaxel
  efflux is strongly **pH-dependent**, so any single ER is valid only at its
  measured pH.
* *In vivo knockout:* Kemper 2003 (PMID 12855665) reports an **11-fold** increase
  in paclitaxel **brain uptake** in P-gp knockout mice. That is a brain *level*,
  not `Kp(KO)/Kp(WT)`; plasma AUC also rises 1.3–2× in the same knockouts
  (Sparreboom 1997; Lagas 2006, PMID 17062689), so the ratio fold is smaller and
  **deriving it would be inventing a number**. This is the same rejection the
  module already applies to methotrexate's 1.8–1.9× plasma-clearance figure.
* Lagas 2006 additionally shows **MRP2/ABCC2 is a co-equal determinant** of
  paclitaxel disposition i.v., so a P-gp-only term is under-specified for it.

So paclitaxel's row would read `substrate_but_no_published_er`, factor **1.0**,
and the panel differential would be **unchanged at 4.94×**. Adding docetaxel's
sourced ER ≈ 20 changes the differential not at all, because the differential is
set by paclitaxel over doxorubicin.

**`tmesim/molecular/efflux.py` was therefore left unchanged.** The four candidate
rows above are fully cited and ready, but adding them would change
`substrate_summary()` from 18 drugs to 22 and silently stale the counts quoted in
three places in `results/efflux_impact.md`, which is generated by
`tmesim/validation/efflux_impact.py` — a file this work does not own. The
evaluation above was done through the module's own sanctioned
`efflux_ratio=` override, which its docstring provides for exactly this case.
**Recommended follow-up: add the four rows in a change that can re-run
`efflux_impact.py`.**

### 2.5 The better-sourced diagnosis for paclitaxel is not efflux

The literature pass turned up a **tissue-scale measurement in Pruijn's own cell
line**, for two of the panel drugs:

> **Choi MS, Kim SH, Kuh HJ.** *Oncol Rep* 2011;25(3):863–870, PMID **21225235**,
> [DOI](https://doi.org/10.3892/or.2011.1138). Paclitaxel and 5-FU through
> multicellular layers of **HT-29** and DLD-1. Paclitaxel penetration was **38 %
> at 24 h at 1 µM** in HT-29, **80 % at 48 h at 1 µM**, but **100 % at 10 µM by
> 48 h**. *"PTX showed limited penetration in both MCLs in contrast to the full
> penetration seen by 5-FU."*

Three things follow, and they point away from efflux:

1. **HT-29 is P-gp-negative (§2.2 ii), so the paclitaxel restriction measured
   here cannot be P-gp efflux.**
2. The restriction is **saturable** — 38 % at 1 µM against 100 % at 10 µM. That
   is the signature of **binding/sequestration**, not of a pump, and paclitaxel
   binds tubulin avidly.
3. Saturable intracellular binding is precisely the **`K > 1` retardation term
   that `results/two_route_transport.md` §6 already lists as missing** (known
   missing physics, item 1). The model represents a partitioning solute's extra
   *conductance* and not its extra *capacity*.

Corroborating, in the same direction: Grantab, Sivananthan & Tannock, *Cancer Res*
2006;66(2):1033–1039, PMID 16424039, attribute differential paclitaxel MCL
penetration to **cell packing density and adhesion, not efflux**; Al-Abd 2008
(PMID 18271941) reproduces the HT-29 > DLD-1 ordering. The one paper that *does*
implicate P-gp at tissue scale — Martin 2003, PMID 14562035, BODIPY-taxol
accumulating in MCF7(WT) but not P-gp-expressing MCF7(Adr) spheroids, restored by
nicardipine — uses a fluorescent conjugate in a drug-selected line and yields
**no extractable numeric fold change**.

**So the paclitaxel row remains the model's most questionable output, the
contradiction flagged in the parent §5.1 is NOT resolved by efflux, and the
best-sourced explanation is a capacity term this model does not have.** That is a
negative result for efflux and it is the honest one.

---

## 3. Blocker 2 — the H-bond-donor decrement

### 3.1 What was and was not retrievable

**Not retrievable, and named as such:**

* **Abraham MH, Chadha HS, Whiting GS, Mitchell RC.** *J Pharm Sci*
  1994;83(8):1085–1100, PMID **7983591**,
  [DOI](https://doi.org/10.1002/jps.2600830806) — the canonical Δlog P citation,
  and the one `results/two_route_transport.md` §7 named as the right source.
  **HARD-CLOSED.** Unpaywall (`oa_status: closed`, empty `oa_locations`),
  OpenAlex (`is_oa: false`) and Semantic Scholar (`CLOSED`) independently report
  no open-access copy anywhere; no secondary source quoting its regression
  equation verbatim was located. Its abstract is verified and states *"An analysis
  of 288 delta log P values shows that solute hydrogen-bond acidity is the major
  factor"* — but carries **no coefficients**. Its fitted coefficients are recorded
  **UNVERIFIED** and are **not used**.
* **Seiler P.** *Eur J Med Chem* 1974;9:473–479 — the `I_H` increments. Predates
  DOI assignment; no full text, no scan, and no secondary source quoting the
  numeric table. **Zero verified values.** Any `I_H` number in circulation should
  be treated as unsourced until someone photographs the page. Recorded
  **UNVERIFIED**.

**Retrievable, verbatim, open-access:** the two constituent Abraham LSER
coefficient sets, from which the Δlog P donor coefficient follows by subtraction.

### 3.2 The sourced value, derived in one line

Abraham's LSER is `log P = c + e·E + s·S + a·A + b·B + v·V`, with `A` the solute's
hydrogen-bond **acidity** (donor strength). The octanol → alkane correction is the
**difference** of two such equations, so its donor term is `(a_oct − a_alk)·A`:

| system | `a` (donor / acidity coefficient) |
|---|---:|
| water → wet octan-1-ol | **+0.034** |
| water → hexadecane | **−3.587** |
| **Δlog P donor coefficient** | **3.621 per unit A** |

> Both read verbatim from the open-access wet-solvent water-to-organic log P
> coefficient table of **Acree WE Jr, Grubbs LM, Abraham MH**, "Prediction of
> Partition Coefficients and Permeability of Drug Molecules in Biological
> Systems with Abraham Model Solute Descriptors…", in *Toxicity and Drug
> Testing*, IntechOpen 2012 (octan-1-ol wet `a` = 0.034; hexadecane
> `a` = −3.587; cyclohexane `a` = −3.740). Independently confirmed in Table 1 of
> **Abraham MH, Acree WE Jr**, "Descriptors for Pentane-2,4-dione and Its
> Derivatives", *J Solution Chem* 2017;46(8):1625–1638,
> [DOI](https://doi.org/10.1007/s10953-017-0667-y), PMC5579141, which carries
> both values.
>
> **Not the source**, despite earlier revisions of this file: *ACS Omega*
> 2019;4(2):2883–2892, [DOI](https://doi.org/10.1021/acsomega.8b03242),
> PMC6648601 — its Table 1 has no hexadecane row. Earlier revisions also
> cross-checked against an unnamed "Acree/Abraham compilation" quoting hexane
> −3.578 and heptane −3.571; the IntechOpen table prints **−3.599** and
> **−3.317**, so those two values did not reproduce and the cross-check has been
> withdrawn. Neither enters any calculation. See
> `results/abraham_citation_erratum.md`.

The model's input is an **integer donor count** (`NumHDonors`), not `A`, so a
typical `A` per donor group is required. Measured `A` values, verified in the
open-access tables of PMC5830486 and PMC3477586:

| donor group | `A` | ⇒ decrement = 3.621·A |
|---|---:|---:|
| aliphatic amine | 0.16 | 0.579 |
| aniline NH₂ | 0.26 | 0.941 |
| **aliphatic OH** | **0.37** | **1.340** |
| benzyl alcohol | 0.39 | 1.412 |
| benzamide NH₂ | 0.49 | 1.774 |
| phenol OH | 0.60 | 2.173 |
| carboxylic acid | 0.62 | 2.245 |

**What this is and is not.** It is a value **derived in one line from two
independently measured coefficient sets**, not a directly measured per-donor
decrement — no such quantity exists in any primary source that could be reached.
**One modelling choice remains and is stated rather than buried:** which donor
group is typical. This converts one unsourced number into a sourced formula with
one documented choice; it does not make the choice disappear. The whole measured
range is swept in §3.3.

**Direction of bias, and it runs against this number.** Hexadecane is the
*harshest available* membrane proxy. A real bilayer's interfacial region (ester
carbonyls, phosphate) supplies hydrogen-bond acceptors, so the donor penalty for a
**liposome/water** partition should be **smaller** than 3.621·A. No Abraham
coefficient set for a DMPC/liposome-water system could be verified — flagged
**UNVERIFIED** rather than corrected for. It means **1.34 is an upper estimate for
a bilayer**, and that the pre-registered 1.0 sits *inside* the plausible band
rather than below it.

**Two housekeeping notes that matter.** (i) `Lipinski.NumHDonors` counts an NH₂ as
**one** donor while `NHOHCount` counts it as **two**; the `A` values above are per
**group**, so they pair with `NumHDonors`. The engine uses `hbd` = `NumHDonors`,
which is the correct pairing. (ii) A flat per-donor count is physically wrong for
chelating motifs — the ClogPalk parameter file (Kenny, Montanari & Prokopczyk,
*J Comput Aided Mol Des* 2013;27:389–402; parameter file read in full) carries
explicit intramolecular-H-bond terms such as −0.7 for catechol and −4.0 for
8-hydroxyquinoline. Group increments are the better long-term answer; they are not
drop-in replacements here because they are decrements on `log P_alk` relative to a
saturated hydrocarbon of equal surface area, not on Δlog P(oct−alk).

### 3.3 Scored — and the sourced value is not the optimum

`--stage hbd`. 68 compounds, nothing fitted, same harness as the parent file.

| decrement | provenance | spread | **R²** | Spearman | geo bias | within 3× |
|---:|---|---:|---:|---:|---:|---:|
| 0.000 | **parameter-free**, nothing unsourced | 6.4× | 0.345 | +0.703 | 2.20× | 66 % |
| 0.579 | sourced, A = 0.16 aliphatic amine (low end) | 6.5× | 0.598 | +0.843 | 1.66× | 88 % |
| 0.941 | sourced, A = 0.26 aniline NH₂ | 6.6× | 0.722 | +0.840 | 1.30× | 96 % |
| 1.000 | the pre-registered **convention** | 7.1× | **0.725** | +0.842 | 1.25× | 96 % |
| **1.340** | **SOURCED, A = 0.37 aliphatic OH (central)** | 8.3× | **0.718** | +0.821 | **1.03×** | **99 %** |
| 1.412 | sourced, A = 0.39 benzyl alcohol | 8.4× | 0.716 | +0.817 | 1.00× | 99 % |
| 1.774 | sourced, A = 0.49 benzamide NH₂ | 8.5× | 0.679 | +0.814 | 0.86× | 97 % |
| 2.173 | sourced, A = 0.60 phenol OH | 8.7× | 0.603 | +0.791 | 0.77× | 94 % |
| 2.245 | sourced, A = 0.62 carboxylic acid (high end) | 8.8× | 0.588 | +0.788 | 0.76× | 94 % |

**The sourced value scores worse on R² than the convention it replaces, and that
is the point.** 0.718 against 0.725 — a 1 % difference — and it is **better** on
geometric bias (1.03× against 1.25×, i.e. essentially unbiased) and on within-3×
coverage (99 % against 96 %). The R² optimum remains at 1.000, so the disclosure
in `results/two_route_transport.md` §4.4 stands unchanged and unweakened; what has
changed is that the *adopted* value no longer sits on it.

**The robustness statement that matters commercially.** Across the **entire
measured range** of the one remaining choice — every donor group from aliphatic
amine to carboxylic acid, 0.579 to 2.245 — R² stays between **0.588 and 0.722**,
against **0.022** for the shipping engine. The result does not depend on getting
the choice right.

**The parameter-free floor is unchanged and still available:** decrement 0 gives
R² = 0.345 and Spearman +0.703 with nothing unsourced anywhere. §5 says which to
prefer.

---

## 4. Blast radius

### 4.1 The correction, three arms

`--stage blast`, published `dx = 12 µm` breast grid, both α arms, each drug at its
own two-route `D`.

| drug | scale @ 1.0 | vs dox | **scale @ 1.34 (sourced)** | **vs dox** | scale @ 0.0 | vs dox |
|---|---:|---:|---:|---:|---:|---:|
| doxorubicin | 0.1512 | 1.0000 | 0.1512 | 1.0000 | 0.7949 | 1.0000 |
| **paclitaxel** | 0.7466 | **4.9393** | **0.2416** | **1.5985** | 0.9501 | 1.1952 |
| docetaxel | 0.2214 | 1.4645 | 0.1527 | 1.0102 | 0.9403 | 1.1829 |
| epirubicin | 0.1512 | 1.0000 | 0.1512 | 1.0000 | 0.7949 | 1.0000 |
| 5-fluorouracil | 0.1571 | 1.0392 | 0.1524 | 1.0083 | 0.4749 | 0.5974 |
| cisplatin | 0.2077 | 1.3737 | 0.1637 | 1.0832 | 0.8048 | 1.0125 |
| **panel differential** | | **4.939×** | | **1.599×** | | **2.001×** |

**This is the finding that couples the two blockers.** At the sourced decrement
the paclitaxel outlier that motivated the entire efflux question **largely
disappears without any efflux term**: 4 donors × 1.34 = 5.36 subtracted from
Crippen logP 3.74 puts `log K_mem` at −1.62 instead of −0.26, membrane
permeability at 2.16 µm/s instead of 49.4, and the membrane resistance ratio at
8.39 instead of 0.366 — i.e. **paclitaxel's membrane is no longer nearly open**.
The 4.94× was substantially an artefact of an under-sourced donor penalty.

Note also what the **parameter-free** arm does, because it is the argument against
it: at decrement 0, doxorubicin's own membrane opens (`log K_mem` = 0, `P` = 106
µm/s) and it acquires a **5.26× transcellular gain**. Doxorubicin is a
notoriously poorly-permeating, strongly tissue-binding anthracycline and is the α
calibration anchor. **The parameter-free variant buys "nothing unsourced" at the
price of physics that is wrong about the anchor drug.**

### 4.2 The parent's *other* stated blocker also moves, and one of its premises is false

`results/two_route_transport.md` §7 reason 2 was: *"Five of the six worst
residuals are lipophilic over-predictions (§4.3). The top of the logP range is
held only by the clamp at 3, and that clamp was adopted from a different function
for a different purpose."* Both halves change.

**The residual structure stops being one-sided.**

| | worst-six composition | ratio range | within 3× |
|---|---|---:|---:|
| decrement 1.0 | **5 over-predictions, 1 under** | 0.41 – 3.98 | 96 % |
| **decrement 1.34 (sourced)** | **3 over, 3 under** | 0.23 – 2.99 | **99 %** |

At the sourced value the systematic lipophilic over-prediction the parent flagged
is gone; the worst residual is now an *under*-prediction (SN26434, 0.23×) and the
error is two-sided. Coverage improves 96 % → 99 %. The worst single ratio does get
worse in the under direction (0.41× → 0.23×), and that is the honest cost.

**The clamp premise is false, and it is false in a way that matters for the
choice between variants.** `LOGK_MEMBRANE_MAX = 3.0` binds on **0 of 68** Pruijn
compounds at decrement 1.0 *and* at 1.34 — the most lipophilic compound in the set
reaches `log K` ≈ 1.96. So the clamp is doing no work on the anchor dataset and
cannot be what holds the top of the range. But on the **breast panel**:

| decrement | compounds pinned at the clamp |
|---|---|
| 1.00 | none |
| **1.34 (sourced)** | **none** |
| 0.00 (parameter-free) | **paclitaxel (3.74) and docetaxel (3.26) — both clamped** |

**Only the parameter-free variant puts drugs on that clamp, and it puts both
taxanes there.** Their entire membrane permeability is then set by a constant
adopted from a different function for a different purpose — which is precisely
the objection the parent raised. This is a second, independent argument against
the parameter-free variant, alongside its opening of doxorubicin's membrane
(§4.1).

### 4.3 Verdicts

**The doxorubicin-anchored *uniform* correction is unchanged from the parent, and
that part needs no new solve.** Doxorubicin's membrane is shut at both decrements,
so its two-route scale is **0.151158551795757 at 1.0 and 0.15115449027238376 at
1.34** — a relative difference of **2.7e-5**, four orders of magnitude inside the
1e-3 tolerance the harness matches `D` points on. The uniform-correction column is
therefore numerically the parent's, and its published verdicts carry over
verbatim: **1 of 12 moves in `frozen` (docetaxel FAIL → PASS), 0 of 12 in
`recal`** (`results/two_route_transport.md` §5.2 Table A). Epirubicin is identical
to doxorubicin throughout.

What genuinely changes is **Table B, each drug at its own `D`** — and it changes
in the direction of *less* movement, because every non-anchor drug's scale moves
**closer** to doxorubicin's:

| drug | vs dox @ 1.0 | **vs dox @ 1.34** | change |
|---|---:|---:|---|
| paclitaxel | 4.9393 | **1.5985** | 3.1× closer to the anchor |
| docetaxel | 1.4645 | **1.0102** | 1.4× closer |
| cisplatin | 1.3737 | **1.0832** | 1.3× closer |
| 5-fluorouracil | 1.0392 | **1.0083** | 1.0× closer |
| epirubicin | 1.0000 | 1.0000 | unchanged |

**The parent's Table B found exactly one verdict movement — docetaxel FAIL → PASS,
in `frozen` only — and it was driven by docetaxel's scale (0.2214) sitting
*above* the uniform one (0.1512).** At 1.34 docetaxel's scale is **0.1527**, i.e.
essentially *on* the uniform value (1.0102× the anchor). Its per-drug row
therefore collapses onto the uniform row, and the parent's own monotonicity
argument applies unchanged: the ratio is monotonically increasing in `D` at fixed
α, so a scale that is 1 % above the uniform one gives a ratio 1 % above the
uniform one. **Docetaxel's `frozen` FAIL → PASS movement, the only verdict that
moved in the parent, is expected to persist; its `recal` FAIL is expected to
persist as a FAIL, marginally further outside its band.**

**Solved.** `--stage blast --only 1.34`, published `dx = 12 µm` grid, 458×628,
167,780 cells / 89,785 tumour, same `dt`, same 6 h window, same seed, same
`PDConfig`, same p99 bath-matching convention, same frozen literature bands. α is
the published **148.6082** in `frozen` and re-solved in `recal`. Four rows are
published **NOT TESTABLE** (paclitaxel and 5-fluorouracil, both arms) — single-study
bands — leaving 8 scored verdicts.

| mode | drug | scale | ratio ship | ratio corr | verdict → verdict | α | moved? |
|---|---|---:|---:|---:|---|---:|:--:|
| frozen | doxorubicin | 0.1512 | 0.2301 | 0.1403 | PASS → PASS | 148.6082 | no |
| frozen | paclitaxel | 0.2416 | 0.1503 | 0.0914 | NOT TESTABLE | 148.6082 | no |
| **frozen** | **docetaxel** | 0.1527 | **0.1484** | **0.0731** | **FAIL → PASS** | 148.6082 | **YES** |
| frozen | epirubicin | 0.1512 | 0.2301 | 0.1403 | PASS → PASS | 148.6082 | no |
| frozen | 5-fluorouracil | 0.1524 | 0.3379 | 0.2559 | NOT TESTABLE | 148.6082 | no |
| frozen | cisplatin | 0.1637 | 0.3085 | 0.2200 | PASS → PASS | 148.6082 | no |
| recal | doxorubicin | 0.1512 | 0.2301 | 0.2330 | PASS → PASS | 22.2234 | no |
| recal | paclitaxel | 0.2416 | 0.1503 | 0.1522 | NOT TESTABLE | 35.6934 | no |
| **recal** | **docetaxel** | 0.1527 | **0.1484** | **0.1517** | **FAIL → FAIL** | 22.4521 | **no** |
| recal | epirubicin | 0.1512 | 0.2301 | 0.2330 | PASS → PASS | 22.2234 | no |
| recal | 5-fluorouracil | 0.1524 | 0.3379 | 0.3353 | NOT TESTABLE | 22.4093 | no |
| recal | cisplatin | 0.1637 | 0.3085 | 0.3073 | PASS → PASS | 24.0953 | no |

**Exactly 1 of 12 verdicts moves — docetaxel FAIL → PASS, and only in `frozen`.
In `recal`, the mode the codebase actually runs, no verdict moves.** That is the
same verdict pattern as the parent at decrement 1.0, and the analytic prediction
above landed: docetaxel's `frozen` movement persisted and its `recal` FAIL
persisted. **The sourced decrement changes the blast radius' magnitude, not its
verdicts.**

Note also that the `recal` α values are now much tighter across the panel —
22.22 to 35.69, against the parent's 22.22 to 110.88 — because every drug's scale
has moved toward the anchor's. The panel is closer to being describable by a
single α, which is what a smaller differential means.

**The one cost, and it is a real one.** `results/cellular_exclusion.md` §5.4
flagged that a uniform correction drove both taxanes below one grid cell per
boundary layer at `dx = 12 µm`, and `results/two_route_transport.md` §5.1
recorded that the two-route correction was **gentler** on exactly those two drugs
because their membranes were open. Closing paclitaxel's membrane gives that
hazard back:

| correction (frozen arm) | paclitaxel `L/dx` | docetaxel `L/dx` |
|---|---:|---:|
| cell-exclusion, uniform 0.1512 | 0.93 | 0.93 |
| two-route @ decrement 1.0 | **2.07** | **1.13** |
| **two-route @ decrement 1.34 (sourced)** | **1.18** | **0.94** |
| two-route @ decrement 0.0 | 2.34 | 2.33 |

**At the sourced value docetaxel falls back below one grid cell per boundary
layer (0.94), and paclitaxel drops from 2.07 to 1.18.** The parent's "one
resolution hazard is *better* here" claim does **not** survive the sourcing. This
does not change a verdict — docetaxel is FAIL at every grid tested and the
movement is in `frozen` only — but it means the `frozen` numbers for both taxanes
are quoted on a grid that no longer resolves the gradient producing them, and a
flip should be accompanied by a grid-convergence check at these scales rather
than assuming `results/grid_convergence_breast.md` still covers them.

---

## 5. Recommendation

**The decision is the reader's. This section produces the evidence and a
recommendation, and flips nothing.**

### 5.1 What I recommend

**Recommendation: YES, flip — but do the retardation term first.**

Concretely, the configuration I recommend flipping *to*, when it is flipped:

| | value | why |
|---|---|---|
| `two_route` | **ON** | §5.2 |
| `hbd_decrement` | **1.34**, re-registered | §3, §5.4 last bullet |
| efflux | **OFF** | §2, and now for much stronger reasons |
| α mode | **`recal`** | `frozen` silently breaks the measured Primeau anchor |

**But not in this change, and not before the `K > 1` retardation term
(§5.3 item 1).** The parent's five-item queue is now a two-step sequence — add
the capacity term, then flip — rather than the five it was. If the retardation
work stalls, flipping anyway is defensible provided paclitaxel's row carries the
§2.5 caveat explicitly; that is a judgement call and it belongs to the reader,
not to this file.

The parent §7 gave **four reasons not to flip**. Two are closed or materially
weakened, two are untouched:

| parent §7 reason not to flip | status after this work |
|---|---|
| 1. The one unsourced parameter sits at its own sweep optimum | **CLOSED.** §3. Sourced at 1.34 from open-access Abraham coefficients; it is *not* the optimum and scores marginally worse on R². |
| 2. Five of six worst residuals are lipophilic over-predictions, and the top of the logP range is held only by the clamp at 3 | **HALF FALSIFIED, HALF IMPROVED.** §4.2 — the clamp binds on **0 of 68** Pruijn compounds, so it holds nothing; and at 1.34 the residuals become two-sided (3 over / 3 under) with coverage 96 % → 99 %. |
| 3. Every unresolved item in `cellular_exclusion.md` §7 is still unresolved | **OPEN, unchanged.** |
| 4. The descriptor layer is wrong by ~2.4 log units for cisplatin | **OPEN, unchanged, and owned elsewhere.** |

The parent §7 also gave a **five-step order of work**, of which step 5 was the
flip itself. Step 1 (obtain the Abraham conversion) is **done**. Step 2 (the
retardation basis) is **not done but is now much better specified** — §2.5 gives
it a measured, same-cell-line, saturable target it did not have. Steps 3 (second
cell line) and 4 (the `brownian.py` / `ALPHA_TORTUOSITY` double-count) are
**untouched**.

Separately, and not on either list: the efflux contradiction that the parent §5.1
called *"the largest single number in the blast radius, and … the reason §7
recommends against flipping the default"* is **resolved as far as it can be** —
see §5.2 item 3.

### 5.2 Why flip, in order of weight

1. **The unsourced parameter is gone, and its replacement is not flattering.**
   1.34 is derived in one line from two independently measured, open-access
   Abraham coefficient sets. It scores **worse** on R² than the convention it
   replaces (0.718 vs 0.725) and better on bias (1.03× vs 1.25×) and coverage
   (99 % vs 96 %). The single most attackable fact about the parent file — that
   the adopted value sat on its own sweep optimum — no longer holds for the
   adopted value.
2. **The blast radius shrinks by a factor of three.** The panel differential
   falls **4.94× → 1.60×** and the paclitaxel outlier that was *"the largest
   single number in the blast radius, and the reason §7 recommends against
   flipping the default"* falls with it (§4.1). It was substantially an artefact
   of an under-sourced donor penalty.
3. **The efflux question is closed rather than pending.** Not by adding efflux,
   but by establishing with three sourced reasons that it is not testable against
   the anchor and not applicable to the drug that needed it (§2.2, §2.4). A
   pending unknown has become a documented boundary.
4. **The result is robust to the one remaining judgement.** Across the entire
   measured range of donor acidity — 0.579 to 2.245 — R² stays 0.588–0.722
   against a baseline of 0.022 (§3.3). Flipping does not bet on getting the
   choice right.
5. **`recal`, not `frozen`.** `recal` is the mode the codebase actually runs, and
   `frozen` silently breaks the measured Primeau penetration anchor. This is
   carried over unchanged from the parent §5.2 and is not re-litigated here.

### 5.3 Why not in the same change, and what would change my mind

1. **The model is now known to be wrong for taxanes in a specific, measured
   way.** §2.5 is not a hedge: Choi 2011 measures paclitaxel at 38 % penetration
   through an HT-29 multicellular layer at 1 µM against 100 % at 10 µM, in a
   P-gp-negative line. That is saturable intracellular binding, the model has no
   capacity term, and paclitaxel is on the panel. Shipping a default that is
   known-wrong for a named drug against a named measurement is a different act
   from shipping it as opt-in. **Add the `K > 1` retardation basis first, or
   flip with paclitaxel's row explicitly carrying that caveat.**
2. **One cell line, still.** Pruijn's own paper reports its asymptote
   coefficients moving up to 100-fold across SiHa, HCT8-Ea and HCT8-Ra. Nothing
   in this work tested that movement. A model claiming to be physics rather than
   a fit should be scored on a second cell line before it becomes the default,
   and that remains the cheapest unbought piece of evidence.
3. **The descriptor layer still feeds it a logP that is ~2.4 log units wrong for
   cisplatin**, and at a decrement of 1.34 the membrane leg is *more* sensitive
   to donor count than before, not less.
4. **The taxane grid-resolution hazard comes back (§4.3, last table).** Sourcing
   the decrement closes paclitaxel's membrane, which drops both taxanes' `L/dx`
   — docetaxel to **0.94**, back below one cell per boundary layer. No verdict
   turns on it, but a flip should carry a grid-convergence check at these scales
   rather than inheriting `results/grid_convergence_breast.md`.

**What would change my mind toward flipping immediately:** a second-cell-line
score that holds up. **What would change my mind against flipping at all:** a
verified DMPC/liposome-water Abraham coefficient set showing the bilayer donor
penalty is far below 3.621·A, which would move the sourced value materially below
1.0 and reopen §3.

### 5.4 What I recommend *against*

* **Against shipping the parameter-free variant as the default.** It is the
  tempting answer — nothing unsourced anywhere — and it is worse physics on two
  independent counts. It opens doxorubicin's membrane to a 5.26× transcellular
  gain (§4.1), and doxorubicin is both a notoriously poorly-permeating
  anthracycline and the α calibration anchor; and it puts **both taxanes exactly
  on the `logK` clamp of 3** (§4.2), which is the very constant the parent
  objected to. Keep it as the disclosed floor it already is (R² = 0.345), not as
  the shipped default.
* **Against turning efflux on**, for the reasons in §2. Its default-off status is
  unchanged, but the *reason* is now much stronger than `results/efflux_impact.md`
  had: not merely that it failed to improve a scorecard, but that it is
  structurally unable to move the α anchor (§2.1), untestable against the
  transport anchor (§2.2), and unsourceable for the drug it was needed for
  (§2.4).
* **Against changing `HBD_LOGP_DECREMENT` in `two_route.py` without
  re-registering it.** It is deliberately left at the pre-registered 1.0 so that
  `results/sweep/tworoute_prereg.json` and its checking stage stay valid. 1.34
  lives in `experiments/two_route_efflux.py::HBD_ABRAHAM` and is passed
  explicitly. A flip should move it in a change that re-registers it.

### 5.5 The follow-ups this produces, in value order

1. **Add the `K > 1` retardation/capacity basis** and score it against Choi 2011
   (PMID 21225235), which supplies a saturable concentration dependence in the
   anchor's own cell line. This is now the top item and it has a target.
2. **Score against a second cell line.** Unchanged from the parent, still the
   cheapest unbought evidence.
3. **Add the four cited efflux rows** (docetaxel ER ≈ 20, paclitaxel/epirubicin
   substrate-only, 5-FU unknown) in a change that can re-run
   `tmesim/validation/efflux_impact.py` — see §2.4 for why they were not added
   here. All four are recorded verbatim in
   `results/sweep/tworoute_efflux_prereg.json`.
4. **Retrieve a DMPC/liposome-water Abraham coefficient set.** It is the one
   literature item that could move the sourced decrement materially, and it would
   move it *down*.
5. **Consider a BCRP/ABCG2 term.** It is the one efflux route HT-29 does not rule
   out (§2.2), and it is currently unmodelled.

---

## 6. What this licenses saying, and what it does not

**Established:**

* Efflux composes with the two-route law through an algebraic identity with zero
  free parameters, verified to 3.7e-16, and the composition makes a falsifiable
  prediction (a shut membrane is immune to efflux) that it satisfies.
* Efflux is **not commensurable** with the Pruijn dataset, for three independent
  reasons, the strongest of which is a sourced fact about the calibration tissue:
  HT-29 is functionally P-gp-negative.
* The H-bond-donor decrement is now **sourced**, from open-access Abraham
  coefficients, at 1.34; it is not the sweep optimum and it scores marginally
  worse on R² and better on bias and coverage.
* At the sourced decrement the panel differential falls from 4.94× to 1.60×.

**Not established, and must not be claimed:**

* **That paclitaxel is now explained.** It is not. §2.5 relocates the diagnosis
  from efflux to saturable intracellular binding, on the strength of a
  tissue-scale measurement in the right cell line — but the model has no capacity
  term, so it cannot represent that, and paclitaxel remains its most questionable
  row. Its band is single-study and it is published NOT TESTABLE.
* **That efflux is refuted.** It is *unmeasurable here*. The mechanism is real
  and the composition is correct; the literature will not supply the input for
  the drug that needs it, and the anchor dataset cannot test it.
* **That 1.34 is a measured constant.** It is a one-line derivation from two
  measured coefficient sets, with one stated remaining choice, and hexadecane
  makes it an upper estimate for a bilayer.
* **That the engine's diffusivity is now correct.** Everything in
  `results/two_route_transport.md` §6 still binds, unchanged.

**Known missing physics, updated:**

1. **The `K > 1` retardation/capacity basis is now the top item, not the second.**
   §2.5 gives it a measured, same-cell-line, saturable signature to be tested
   against (Choi 2011, PMID 21225235: 38 % at 1 µM vs 100 % at 10 µM through
   HT-29 MCLs). That is a sharper target than the parent file had.
2. **A BCRP/ABCG2 term is the one efflux route the anchor tissue does not rule
   out** (§2.2). Unactionable today for want of ratios on the compound set.
3. Everything else in `results/two_route_transport.md` §6 and
   `results/cellular_exclusion.md` §7 carries over unchanged.

---

## 7. Reproducing this

```
./.venv/bin/python experiments/two_route_efflux.py --stage algebra          # <1 s
./.venv/bin/python experiments/two_route_efflux.py --stage commensurability # <1 s
./.venv/bin/python experiments/two_route_efflux.py --stage panel            # ~5 s
./.venv/bin/python experiments/two_route_efflux.py --stage hbd              # ~60 s
./.venv/bin/python experiments/two_route_efflux.py --stage blast            # long
./.venv/bin/python -m tmesim.molecular.two_route                            # self-test
./.venv/bin/python -m pytest tests/test_two_route.py -q
```
