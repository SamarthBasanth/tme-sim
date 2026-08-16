# DRAFT — NOT SENT. FOR HUMAN REVIEW, EDITING AND SENDING.

> **Status: unsent draft.** No contact has been made with Prof. Agar, Prof. Santagata, Prof.
> White, or anyone at BWH, DFCI, HMS or MIT. They do not know this project exists and have
> agreed to nothing. Fill every `[BRACKET]`. Do not invent an affiliation or a referral.
>
> **Priority 3** of the access plan. Highest scientific ceiling, lowest probability, longest
> timeline. Send it anyway — the cost is one email.

---

## Target

The paper: Lopez BGC, Kohale IN, Du Z, …, **Santagata S**, Sarkaria J, **White FM**,
**Agar NYR**. *Multimodal platform for assessing drug distribution and response in clinical
trials.* **Neuro-Oncology** 2022;24(1):64–77. **PMID 34383057, PMC8730776.**

| | |
|---|---|
| **Nathalie Y. R. Agar** | Director, Surgical Molecular Imaging Laboratory (SMIL), BWH Neurosurgery. Lab site <https://agarlab.bwh.harvard.edu/> advertises **a web contact form only** — no email or phone is printed on it. The paper's corresponding line gives `nathalie_agar@dfci.harvard.edu`. |
| **Sandro Santagata** | HMS Systems Biology page <https://sysbio.med.harvard.edu/sandros> prints `ssantagata@partners.org` and assistant `cbuzzell@bwh.harvard.edu`. The paper prints `ssantagata@bics.bwh.harvard.edu`. **`santagata.bwh.harvard.edu` is currently down.** |
| **Forest M. White** (MIT, Koch Institute) | Owner of the phosphoproteomics arm and therefore of the PXD018782 TMT design map. Relevant for the softest ask below. |

**Verify all three addresses on the day you send.** Given the ambiguity, the SMIL web form plus
one email to the paper's corresponding address is the safest belt-and-braces route.

**Structural facts to know before writing** (detail in `../DATA_ACCESS_PLAN.md` §4):

- **The PSOC is defunct.** U54 CA210180 ran 2016-08-29 to 2022-07-31; `psoc.mit.edu` is
  NXDOMAIN; NCI lists it under Past Projects. **Do not address them via the PSOC** — it dates
  the sender by four years.
- **SMIL is not a fee-for-service core.** No rate card, no service menu. This is a
  collaboration request or nothing.
- The patient specimens are from **NCT01849146**, the ABTC phase-1 adavosertib + RT + TMZ study
  with a recurrent-GBM surgical window-of-opportunity arm. **Phase-1 patient tissue means IRB
  and a DUA, and a DUA requires a counterparty institution.** If the sender has no affiliation,
  say so early; it changes what is realistic and they will work it out anyway.

## Ask the small thing first

The instinct is to ask for the human MALDI-MSI. Don't, not in email one. Two smaller asks are
far more answerable, involve no patient tissue, and establish that the sender is a serious
person:

1. **The TMT channel → experimental arm design map for PXD018782.** PRIDE PXD018782 is public
   and CC0, but our verification pass found it contains **no quantified phosphosite response** —
   all 80 files are raw spectra or peptide identifications, with zero reporter-ion, abundance
   or ratio fields. Without the channel-to-arm mapping the deposit is close to uninterpretable.
   This is a one-file, no-sensitivity request that would help everyone who downloads it, and it
   is a legitimate thing to raise about a public deposit.
2. **The PDX arm of the MALDI-MSI** (GBM22 / GBM84). Not human-subject data, therefore not
   IRB-encumbered, and it already carries the headline result we most want to test — 1–19 µM
   measured adavosertib with intratumoral CV 20.9–85.6%.

The phase-1 patient MSI and the t-CyCIF montages ("uploaded to Omero" with no accession, no
URL, and nothing findable on Minerva, cycif.org/data or HTAN) can follow a reply.

---

## Subject line options

- `PXD018782 — TMT channel/arm mapping, and a question about the MALDI-MSI`
- `Testing a GBM drug-transport model against measured intratumoral adavosertib`

---

## Draft message

Dear Prof. Agar and Prof. Santagata,

I am writing about the 2022 *Neuro-Oncology* multimodal platform paper (Lopez et al.,
PMID 34383057) with one small request and one larger question.

I maintain an open-source mechanistic model of drug transport through patient-specific
glioblastoma tissue: it reconstructs the microenvironment from spatial transcriptomics — cell
types, vasculature, ECM, oxygen, on real micron coordinates — and solves reaction–diffusion
transport and pharmacodynamics through it. I have spent several months trying to find any
dataset in the world that can test its central claim, and yours is the only one that can.

**Where the model stands, plainly.** It is a research prototype. It is not clinically validated
and nothing in it should inform a treatment decision. An adversarial audit of my own work found
eleven distinct defects, several mutually cancelling and invisible to a green test suite; every
headline number I published before that audit was wrong in at least one way. One component was
falsified outright and retracted: a penetration metric that could not discriminate between drugs
by construction. A second retraction has since been **retracted in its turn** — I had demoted my
own oxygen PDE for losing to plain distance-to-nearest-vessel, and that was an inverted
inequality: the predictors are sign-aligned, so |−0.1523| > |−0.1204| means the solve *wins*, by
1.26×. I record it because it is the most instructive failure in the project — a component was
deprecated on a misread inequality while every mechanical check in the same repository disagreed
with the prose — and because the corrected reading is still modest: at ρ ≈ 0.15, within-section,
the solve is the best of four weak predictors, not a validated oxygen field, and its diffusion
length is still 2–5× the published value because its consumption term has never been fitted.
What survives is the delivery axis, and it is scored: against 125 human glioma
tissue samples across 42 patients plus 27 PMID-anchored regional PK measurements, predicted
tumour:plasma partition lands within about a **2× geometric-mean fold error**, 6/6 drugs inside
10×, with the BBB-penetrant vs BBB-excluded separation reproduced from molecular structure alone
(3.65× predicted, 4.66× measured). To be exact about what that scores: it is a **structure-only
permeability coefficient**, and those rows never invoke the transport solver. Of the checks that
do compare a solved field to an external measurement there are three, and all three fail. I would
not want "the delivery axis is scored" read as "the concentration field is validated." One
failure was informative enough to fix the model: the
measured contrast-enhancing vs non-enhancing ratio is 1.45–4.26× and we predicted exactly 1.00×,
because there was no barrier-integrity term anywhere. The data itself showed the fix — a
multiplicative barrier factor cannot work (3.2× spread across drugs) but an **additive leak is
drug-independent to within 1.3×** across three independent studies (L = 0.211 / 0.230 / 0.272),
which is the signature of size-nonselective paracellular opening.

**What I cannot test at all is response.** I ran a systematic sweep of 164 datasets. There is
no open dataset anywhere that pairs a measured spatial drug distribution with a measured
treatment response in the same tumour, and there is no glioma + temozolomide spatial drug map
in existence — METASPACE returns zero. Your platform is the experiment I would design if I
could design one: measured adavosertib per pixel, in the right disease, with co-registered
phosphoproteomics and t-CyCIF, in phase-1 patient specimens.

**The small request first.** PXD018782 is public and I have downloaded it, but as deposited it
contains raw spectra and peptide identifications with no reporter-ion abundances or ratios, and
no mapping from TMT channel to experimental arm. Would it be possible to share that design map?
It is a single small file, involves no patient data, and would make an already-public deposit
usable by anyone who downloads it.

**The larger question.** The paper states the MALDI-MSI data are available upon request. I
would be glad to discuss the phase-1 specimens if that is ever appropriate, but I recognise
that is human-subject data from NCT01849146 and would require an IRB pathway and a data-use
agreement — [**PICK ONE, HONESTLY:** *and I should say up front that I currently have no
institutional affiliation, which I understand may make a DUA impractical; if so, I would rather
know now.* / *and I would be working through [INSTITUTION], whose office can execute a DUA.*]

So the more tractable version of the question: **would you consider sharing the PDX arm
(GBM22 / GBM84)?** It is not human-subject data, and it already carries the measurement I most
want to test myself against — the 1–19 µM range and the 20.9–85.6% intratumoral CV. A model
that cannot reproduce that CV distribution from tissue architecture alone is wrong, and I would
like to find that out.

**What I would offer.** Before receiving anything I would commit, in a timestamped
version-controlled file, to the predicted intratumoral CV, the predicted spatial
autocorrelation of exposure, and the predicted relationship between local exposure and local
response — so that whatever follows is a prediction and not a fit. I would share the
reconstruction pipeline and the validation harness (permutation tests, leave-one-out
cross-validation and an explicit uncertainty Monte Carlo, built to report where a model fails),
and I would **publish the result if the model fails**, which I would rather state now than
later. Authorship entirely as you prefer, including none.

If the answer is no, or no reply, I completely understand — I am aware this is an unsolicited
email about a four-year-old paper. Thank you for the paper regardless; the reported CV range is
the number I have been designing against.

Best regards,

[YOUR FULL NAME]
[AFFILIATION, or: "Independent researcher — no institutional affiliation"]
[EMAIL]
[REPOSITORY / PROJECT LINK, if public]

---

## Reviewer notes

- **Do not** mention the PSOC, or `psoc.mit.edu`, or suggest routing through it. It ended in
  July 2022.
- **Do not** call SMIL a core facility or ask for a quote. There is no service route.
- **Do not** claim to have read the t-CyCIF data or the MSI. We have neither. The only thing
  downloaded is PXD018782.
- **The affiliation sentence must be honest and must be chosen deliberately.** An unaffiliated
  sender who elides the question wastes everyone's time — a DUA needs a counterparty legal
  entity. Saying it plainly is better received than discovering it in month three.
- **Consider co-sending with Forest White (MIT)** on the TMT design-map ask specifically — it
  is his lab's arm of the paper, and he is a corresponding author.
- **Sequencing.** `../DATA_ACCESS_PLAN.md` §4.3 recommends sending this *after* an institutional
  collaboration exists (Target 1), because it can then be made jointly and has an institution
  behind it. If sent cold anyway, treat gate G2 (8 weeks, no reply → stop) as binding rather
  than following up repeatedly.
