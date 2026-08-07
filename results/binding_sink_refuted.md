# The binding sink is wrong in magnitude and in direction

**The third defect, now characterised.** `results/breast_tissue_permeability.md`
localised the remaining breast-anchor error to the alpha-scaled binding sink
`k = clearance + binding`. This tests that term against measured data and
against the code's own stated bounds. It fails both.

**No parameter was changed.** This is a diagnosis, and the fix it points to is
not one that can be made from inside the model.

## 1. The magnitude is refuted, by a wide margin

The model's loss rates imply interstitial half-lives of:

| drug | model `k` | implied interstitial t½ |
|---|---:|---|
| 5-fluorouracil | 5.586e-3 /s | **2.1 minutes** |
| methotrexate | 2.670e-2 /s | **26 seconds** |

Against measured human tumour extracellular-fluid kinetics:

* **Methotrexate ECF elimination, measured directly.** Blakeley JO *et al.*,
  *J Neurooncol* 2009;91(1):51-8, PMID 18787762, PMC4029102 — human, n=4,
  microdialysis probes in glioma. The paper's own model is
  `dC_ecf/dt = (k13*V1*Cp - k30*V3*C_ecf)/V3` and defines `k30` as "the rate
  constant for elimination of MTX from brain ECF" — first-order, from the
  interstitium, free drug. **This is exactly the model's quantity.** Values:
  2.42, 17.50, 0.52, 0.05 /h, i.e. 6.72e-4 to 1.39e-5 /s, geometric mean
  **2.85e-4 /s**. Terminal ECF half-lives **3.15-14.51 hours**.
* **5-fluorouracil interstitial half-life is not resolvable.** Mader RM *et
  al.*, *Br J Cancer* 2003;88(5):782-7, PMC2376334 — the one group that tried
  reports `t½` as **N.A. in plasma, tumour AND subcutaneous tissue**.

A 26-second interstitial half-life is not a plausible physical quantity for a
drug whose measured tumour-ECF terminal half-lives are 3 to 14 hours.

**Stated precisely, because the measured spread is enormous.** Blakeley's
`k30` spans 350-fold across only 4 patients, so no single measured value is
well constrained. The defensible claim is that the model's methotrexate `k`
(2.670e-2 /s) sits **above the entire measured range** (1.39e-5 to 4.86e-3
/s) — 5.4x above its top and **94x above its geometric mean**. Two further
caveats belong with that: for the two patients with disrupted barrier the ECF
terminal slope is *plasma-limited* ("the t½,z of the C_ecf-time curves was
similar to the C_p-time profiles"), so it is not a pure tissue constant; and
this is brain, not breast.

## 2. The direction is refuted too

The model has methotrexate losing interstitial drug **4.78x faster** than
5-fluorouracil. Four independent lines say otherwise:

1. **The only head-to-head measurement.** Grantab R *et al.*, *Cancer Res*
   2006;66(2):1033-9, PMID 16424039 — multicellular layers, 6 h penetration.
   Normalising each drug to its OWN cell-free control (which removes the
   free-solution diffusivity difference): 15.1 % vs 12.8 % in one line, 40.2 %
   vs 73.6 % in the other. Methotrexate is at parity or **better**, never 4.78x
   worse. Corroborated by Tunggal JK *et al.*, *Clin Cancer Res* 1999;5:1583-6,
   PMID 10389947: "The rate of penetration of 5-FU through MCLs is comparable
   to that of methotrexate."
2. **Methotrexate shows essentially NO net interstitial loss in normal
   peripheral tissue.** Ekström PO *et al.*, *Cancer Chemother Pharmacol*
   1996;37(5):394-400, PMID 8599860 — rat, steady-state infusion, with BOTH
   zero-flow and no-net-flux microdialysis (the two rigorous recovery methods):
   "the level of unbound MTX in the circulation is fairly similar to the
   extracellular concentrations in the muscle and liver", i.e. **Kp,uu ~ 1**.
3. **At the blood-brain barrier the two drugs are indistinguishable.**
   5-FU Kp,uu,brain 0.13 (PMID 25547867) vs methotrexate 0.10-0.17
   (PMID 30951819, 36669926).
4. **Classic methotrexate PBPK models muscle as pure ECF distribution — no
   cellular uptake, no binding at all** (NCBI Bookshelf NBK218077, verified in
   full text).

**The mechanism the model is reaching for is the wrong compartment.**
Methotrexate's famous tissue retention is *intracellular* — polyglutamation
and tight DHFR binding, retained for days. That slows efflux from cells; it
does not accelerate removal from the interstitium. And where methotrexate
genuinely IS lower, the measured cause is **active efflux transport**, not
binding: Sane R *et al.*, *DMD* 2014;42(4):537-40, PMID 24464805, showed an
Abcg2/Abcc4 double knockout roughly doubled brain Kp, 0.11 → 0.23.

A first-order *interstitial* loss term is being asked to carry a mechanism it
does not physically represent.

### 2.1 The compartment convention inverts the sign

This is the sharpest form of the argument, and it cuts the other way from the
naive reading of the literature.

`k` is loss from the **free interstitial** pool. Under that convention,
intracellular trapping is a *loss from the interstitium* and should make `k`
**larger**, not smaller.

* **5-fluorouracil traps.** ~50 % of human tumours show a 19F-MRS tumour
  half-life above the 20-minute trapping threshold (51 of 99 patients, Wolf
  1998 PMID 9859944; 19 of 57, Presant 1994 PMID 7909867), by intracellular
  anabolism to fluoronucleotides.
* **Methotrexate does not trap — it tracks plasma.** Reported independently in
  **five** studies: rat C6 glioma (PMID 10468023), rabbit VX2 muscle r = 0.97
  (PMID 9487560), a human case (PMID 8996531), human glioma ECF "closely
  paralleled the C_p-time profile" (PMID 18787762), and rat muscle/liver
  Kp,uu ~ 1 (PMID 8599860).

So the drug with the fast interstitial sink should be **5-fluorouracil**. The
model has methotrexate 4.78x faster. **Read with the correct compartment
convention, the measured evidence does not merely fail to support the model's
ratio — it inverts it.**

This also disqualifies a tempting shortcut: 19F-MRS half-lives cannot be
substituted for `k`, because they measure *total* tissue drug and, for the
binding component, carry the wrong sign.

## 3. Where the 4.78x actually comes from, in the code

`binding = alpha * tissue_binding_rate(descriptors)`, and

    score = (logP - 2.5)/2.5 + (HBD - 2)/4   ->  sigmoid  ->  [1e-5, 5e-3] /s

| drug | logP (RDKit) | HBD | logP term | HBD term | score |
|---|---:|---:|---:|---:|---:|
| 5-fluorouracil | -0.798 | 2 | -1.319 | **0.000** | -1.319 |
| methotrexate | +0.268 | 5 | -0.893 | **+0.750** | -0.143 |

Structure-derived ratio **4.83x**, matching the residual almost exactly.
**64 % of the gap is the hydrogen-bond-donor term alone.**

Two things are wrong with that:

* **The heuristic's sign is inverted for ionised acids.** It asserts
  "lipophilic + H-bond donors stick to tissue". Methotrexate's high donor count
  comes from its two **carboxylates** — the very groups that make it a dianion
  at pH 7.4 and therefore *prevent* tissue partitioning.
* **The lipophilicity is wrong.** RDKit gives methotrexate logP = +0.268; the
  experimental value is **-1.85** (Hansch/Leo/Hoekman, via PubChem CID 126941;
  5-FU is -0.89 from the same source). Substituting *measured* logP for both
  drugs — not fitting to any anchor — collapses the model's own formula from
  **4.83x to 1.51x**, against a measured 1.02x.

## 4. The code already contradicts itself

`molecular/binding.py` documents its output band as a "physically plausible"
[1e-5, 5e-3] /s. Every raw value obeys it. After the `alpha = 148.6`
amplification, **18 of 19 calibrated library drugs exceed that documented
ceiling**, by up to **29.8x**:

    docetaxel     1.491e-1 /s   29.8x        paclitaxel   1.387e-1 /s   27.7x
    vinorelbine   1.235e-1 /s   24.7x        tamoxifen    1.230e-1 /s   24.6x
    vincristine   8.906e-2 /s   17.8x        irinotecan   6.012e-2 /s   12.0x
    methotrexate  2.663e-2 /s    5.3x

This is an internal inconsistency, checkable without any literature, and it
agrees in sign and rough magnitude with the external refutation in §1.

## 5. What is NOT established

**The 5-FU:methotrexate loss-rate RATIO has no independent measured
constraint.** Methotrexate has exactly one measured interstitial elimination
constant (brain, n=4, spanning 350x); 5-fluorouracil has **none** — Mader 2003
is the only study that attempted it and reported `N.A.` in tumour,
subcutaneous tissue *and* plasma alike, because interstitial 5-FU
(median c_max 0.26 µg/ml) was too low to fit a terminal slope. No `fu,tissue`
value exists for either drug in any tissue in any species.

Because Blakeley's `k30` spans 350-fold, **both** 5-FU rodent tumour values
fall *inside* methotrexate's measured range. A 4.78x ratio is therefore
neither supported nor excluded by measurement — it is **unconstrained**. That
is the reportable null, and it is the honest ceiling on how strong §2's
inversion argument can be made: the direction is contradicted mechanistically
and by the one head-to-head penetration measurement, but no measured *rate
ratio* exists to contradict it numerically.

**There has never been a head-to-head rate comparison.** The one group that
measured both drugs with the same modality in the same tissue (Müller/Vienna,
human breast) published AUC ratios only, and both full texts are
paywalled — AACR returned HTTP 403 on every attempt, across multiple sessions
and access routes.

**A correction to this project's own earlier reasoning.**
`results/breast_tissue_permeability.md` argued that because the two measured
anchors are flat (0.61 vs 0.60) and measured permeability is flat (1.02x),
`k` must be flat too. That inference is **weaker than it was stated**. An
interstitial:plasma AUC ratio is a NET exposure quantity and cannot be
decomposed into influx and efflux — two drugs with different loss rates can
share one. And methotrexate's 0.60 carries **CV = 100 %**, so it is a much
softer anchor than 5-FU's 0.61. The defensible claims here are about the
**magnitude** and **direction** of `k` (§1, §2), not about the anchors proving
equality.

## 6. What would settle it

1. **Poulin P *et al.*, *J Pharm Sci* 2015;104(4):1508-21, PMID 25615572** —
   measured steady-state Kp for BOTH drugs in the SAME mice, including three
   tumour xenografts. The only source with a controlled comparison. Closed
   access; ~10 retrieval routes exhausted. **One reprint request away**, and it
   is the highest-value single action remaining on this line.
2. Müller M *et al.*, *Breast Cancer Res Treat* 2000;60(3):211-7, PMID
   10930108 — the Vienna group's own synthesis, which by design contains both
   anchor cohorts' interstitial concentration-time profiles. Also closed.
3. A purpose-built tumour microdialysis study with a resolvable washout phase
   — which, given Mader's N.A. result, may not be achievable for 5-FU at all.

**Do not fit `alpha` to the breast anchor in the meantime.** It would be a
third parameter absorbing the error of the first two on two data points, and
it would break the out-of-sample mitoxantrone test `alpha` currently passes.
