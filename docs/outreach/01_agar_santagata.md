# DRAFT — NOT SENT. FOR HUMAN REVIEW, EDITING AND SENDING.

> **Status: unsent draft.** No contact has been made with Prof. Agar, Prof. Santagata, Dr Lee,
> Prof. Wen, Prof. White, or anyone at BWH, DFCI, HMS, MIT or the ABTC. They do not know this
> project exists and have agreed to nothing. Fill every `[BRACKET]`. Do not invent an
> affiliation or a referral.
>
> **Priority 3** of the access plan. Highest scientific ceiling, lowest probability, longest
> timeline. Send it anyway — the cost is one email.
>
> **Fact-check pass 2026-08-12** (against PubMed records, not search snippets). Three changes:
> the primary citation is now the **2025 trial report**, not the 2022 platform paper; the
> stale `3.65× vs 4.66×` separation figure was removed; the barrier-leak paragraph was
> rewritten because the fix it describes as needed is now implemented. See "What changed".

---

## Target

**Primary paper (new — cite this one first).** Lee EQ, Alexander BM, Romo CG, Supko JG,
**Agar NYR**, Talebi Z, Ahluwalia MS, Desai AS, Dietrich J, Kaley TJ, Peereboom DM, Gantchev J,
Baquer G, **Santagata S**, Takebe N, Desideri S, Fisher JD, Sims M, Ye X, Ligon KL, Nabors LB,
Grossman SA, **Wen PY**. *Phase I Study of Adavosertib with Radiotherapy and Temozolomide in
Newly Diagnosed Glioblastoma and Intratumoral Drug Levels in Recurrent Glioblastoma.*
**Clin Cancer Res** 2025;31(6):983–992. **PMID 39820282, PMC11913579**,
doi:10.1158/1078-0432.CCR-24-2311.

This is the ABTC trial itself — the source of the phase-1 specimens. It publishes the MALDI-MSI
drug quantitation, and its data-availability statement is the current one. Verified verbatim
from PMC on 2026-08-12:

> "The data that support the findings of this study are available on request from the
> corresponding author. The data are not publicly available due to privacy restrictions."

**Method paper (still cite, second).** Lopez BGC, Kohale IN, Du Z, …, Lee EQ, **Santagata S**,
Sarkaria J, **White FM**, **Agar NYR**. *Multimodal platform for assessing drug distribution and
response in clinical trials.* **Neuro-Oncology** 2022;24(1):64–77. **PMID 34383057,
PMC8730776**, doi:10.1093/neuonc/noab197. Author list and pagination verified against the
PubMed record 2026-08-12. Its statement — "The MALDI-MSI data are available upon request" —
still stands unamended three and a half years on.

| | |
|---|---|
| **Nathalie Y. R. Agar** | BWH Neurosurgery / DFCI. Director, Surgical Molecular Imaging Laboratory (SMIL). Lab site <https://agarlab.bwh.harvard.edu/> advertises **a web contact form only**. The 2022 paper's corresponding line gives `nathalie_agar@dfci.harvard.edu`; **independently corroborated** by the PubMed affiliation record for PMID 40295505 (2025), which prints `Nathalie_Agar@dfci.harvard.edu`. |
| **Sandro Santagata** | BWH Pathology / DFCI. HMS Systems Biology page <https://sysbio.med.harvard.edu/sandros> prints `ssantagata@partners.org` and assistant `cbuzzell@bwh.harvard.edu`. The 2022 paper prints `ssantagata@bics.bwh.harvard.edu`. **`santagata.bwh.harvard.edu` is currently down.** |
| **Corresponding author of the 2025 trial paper** | `[READ OFF THE JOURNAL PAGE — NOT RECORDED HERE]`. The data-availability statement points at "the corresponding author" and PubMed does not carry an email for this record. Eudocia Q. Lee (first author) and Patrick Y. Wen (last author) are both DFCI/BWH; **do not guess which**. Open the AACR article page and read it. |
| **Forest M. White** (MIT, Koch Institute) | Senior author of the phosphoproteomics arm of the 2022 paper and therefore the owner of the PXD018782 TMT design map. Not an author on the 2025 trial paper. *Whether he is a designated corresponding author on the 2022 paper is unverified here* — check before co-addressing. |

**Verify every address on the day you send.**

**Structural facts to know before writing** (detail in `../DATA_ACCESS_PLAN.md` §4):

- **The PSOC is defunct.** U54 CA210180 ran 2016-08-29 to 2022-07-31; `psoc.mit.edu` is
  NXDOMAIN; NCI lists it under Past Projects. **Do not address them via the PSOC.**
- **SMIL is not a fee-for-service core.** No rate card, no service menu. Collaboration or
  nothing.
- **The trial is closed.** NCT01849146 enrolled August 2013 – May 2020 and was terminated
  2020-05-05; 71 patients total, of whom 12 went through the intratumoral-drug-distribution
  arm. Nothing further will be collected. Whatever exists, exists.
- **Phase-1 patient tissue means IRB and a DUA, and a DUA requires a counterparty
  institution.** If the sender has no affiliation, say so early.

## Ask the small thing first

Not the human MALDI-MSI. Two smaller asks are far more answerable, involve no patient tissue,
and establish that the sender is a serious person:

1. **The TMT channel → experimental arm design map for PXD018782.** PRIDE PXD018782 is public
   and CC0, but our verification pass found it contains **no quantified phosphosite response** —
   all 80 files are raw spectra or peptide identifications, with zero reporter-ion, abundance
   or ratio fields. Without the channel-to-arm mapping the deposit is close to
   uninterpretable. One file, no patient data, and it helps everyone who downloads it.
2. **The PDX arm of the MALDI-MSI** (GBM22 / GBM84). Not human-subject data, therefore not
   IRB-encumbered, and it already carries the headline number we most want to be tested
   against — 1–19 µM measured adavosertib with intratumoral CV 20.9–85.6% (2022 paper, GBM84).

The phase-1 patient MSI and the t-CyCIF montages ("uploaded to Omero" with no accession, no
URL, nothing findable on Minerva, cycif.org/data or HTAN) can follow a reply.

---

## Subject line options

- `PXD018782 — TMT channel/arm mapping, and a question about the adavosertib MSI`
- `A GBM transport model against your measured intratumoral adavosertib (CCR 2025)`

---

## Draft message

Dear Prof. Agar and Prof. Santagata,

I am writing about the intratumoral drug-distribution work in your 2025 *Clinical Cancer
Research* paper (Lee et al., PMID 39820282) and the 2022 *Neuro-Oncology* platform paper it
builds on, with one small request and one larger question.

I maintain an open-source mechanistic model of drug transport through patient-specific
glioblastoma tissue: it reconstructs the microenvironment from spatial transcriptomics — cell
types, vasculature, ECM, oxygen, on real micron coordinates — and solves reaction–diffusion
transport and pharmacodynamics through it. It is a research prototype, not clinically
validated, and nothing in it should inform a treatment decision. An adversarial audit of my own
work found eleven distinct defects and retracted two components outright — an oxygen PDE beaten
by plain distance-to-nearest-vessel, and a penetration metric that could not discriminate
between drugs by construction. What survives is the delivery axis, and it is scored: against
125 human glioma tissue samples across 42 patients plus 27 PMID-anchored regional PK
measurements, predicted tumour:plasma partition lands within about a **2× geometric-mean fold
error**, 6/6 drugs inside 10×. Adavosertib is one of those six, from your trial's own reported
tissue PK.

**One number in your 2025 paper is why I am writing.** Your two modalities disagree on region:
LC-MS/MS on homogenates gives a mean intratumoral:plasma ratio of 4.18 in contrast-enhancing
tissue against 0.74 in non-enhancing, while the MALDI-MSI shows "no apparent correlation
between drug concentration and contrast-enhancing versus non-contrast enhancing tumour" across
the ten specimens imaged. That difference is the exact quantity my model produces. Its
regional term is an **additive** barrier leak — derived from the observation that a
multiplicative barrier factor cannot be made drug-independent (3.2× spread across drugs) while
an additive one is drug-independent to within 1.3× across three independent studies — and on
the drugs where regional measurements exist it predicts contrast-enhancing:non-enhancing ratios
of 1.19× and 5.87× against measured 1.45× and 4.26×. I would like to know whether it survives
contact with a within-section measurement, which is something no open dataset in the world can
tell me.

**The small request first.** PXD018782 is public and I have downloaded it, but as deposited it
contains raw spectra and peptide identifications with no reporter-ion abundances or ratios, and
no mapping from TMT channel to experimental arm. Would it be possible to share that design
map? It is a single small file, involves no patient data, and would make an already-public
deposit usable by anyone who downloads it.

**The larger question.** Both papers state the MALDI-MSI data are available on request. I
recognise the phase-1 specimens are human-subject data from a closed trial and would require an
IRB pathway and a data-use agreement — [**PICK ONE, HONESTLY:** *and I should say up front that
I currently have no institutional affiliation, which I understand may make a DUA impractical;
if so, I would rather know now.* / *and I would be working through [INSTITUTION], whose office
can execute a DUA.*]

So the tractable version: **would you consider sharing the PDX arm (GBM22 / GBM84)?** It is not
human-subject data, and it carries the measurement I most want to be tested against — the
1–19 µM range and the 20.9–85.6% intratumoral CV. A model that cannot reproduce that CV
distribution from tissue architecture alone is wrong, and I would like to find that out.

**What I would offer.** Before receiving anything I would commit, in a timestamped
version-controlled file, to the predicted intratumoral CV, the predicted spatial
autocorrelation of exposure, and the predicted relationship between local exposure and local
response — so that whatever follows is a prediction and not a fit. I would share the
reconstruction pipeline and the validation harness (permutation tests, leave-one-out
cross-validation, an explicit uncertainty Monte Carlo, built to report where a model fails),
and I would **publish the result if the model fails**, which I would rather state now than
later. Authorship entirely as you prefer, including none.

If the answer is no, or no reply, I understand completely. Thank you for the papers regardless;
the reported CV range is the number I have been designing against.

Best regards,

[YOUR FULL NAME]
[AFFILIATION, or: "Independent researcher — no institutional affiliation"]
[EMAIL]
[REPOSITORY / PROJECT LINK, if public]

---

## What changed in the 2026-08-12 fact-check

| was | now | why |
|---|---|---|
| Led with the 2022 *Neuro-Oncology* paper; called it "a four-year-old paper" | Leads with **Clin Cancer Res 2025;31(6):983–992, PMID 39820282** | The trial report published 2025-03-17, after this draft was written. It is the source of the phase-1 specimens, it publishes the MSI quantitation, and its data-availability statement is the current one. Writing about a 2022 paper when a 2025 one exists reads as not having checked. |
| "reproduces the measured BBB-penetrant vs BBB-excluded separation (3.65× predicted, 4.66× measured) from molecular structure alone" | **removed** | `results/glioma_regional_pk.md` §P2: 3.65× is the **pre-fix** value. Post-fix the model predicts **17.20×** against measured 4.66× pooled / 9.3× intact-barrier — an over-separation of 1.8×, not a match. `SESSION_AUDIT.md` §1 still carries the pre-fix figure in its table; it is superseded. **Do not quote 3.65× as a current result anywhere.** |
| "there was no barrier-integrity term anywhere … the data itself showed the fix" (written as a gap) | The leak is **implemented** (`tmesim/reconstruct/barrier.py`, L = 0.238) and its current output is quoted: 1.19× / 5.87× against measured 1.45× / 4.26× | Describing an implemented fix as an open gap understates the work and would be corrected by anyone who read the repo. |
| — | Added the MSI-vs-LC-MS/MS regional discrepancy from the 2025 paper as the hook | It is a specific, checkable, expert observation about *their* data, and it is the one quantity this model exists to produce. It replaces a generic "yours is the only dataset that can test it" appeal. |

**Verified unchanged:** author lists, PMIDs, PMCIDs, journal, volume, issue and pagination for
both papers (PubMed, 2026-08-12); Agar's DFCI address, corroborated by a second, more recent
record; the 1–19 µM / 20.9–85.6% CV figures as belonging to the **PDX** arm of the 2022 paper.

## Reviewer notes

- **Do not** mention the PSOC, or `psoc.mit.edu`, or suggest routing through it. It ended in
  July 2022.
- **Do not** call SMIL a core facility or ask for a quote. There is no service route.
- **Do not** claim to have read the t-CyCIF data or the MSI. We have neither. The only thing
  downloaded is PXD018782.
- **Quote their numbers exactly or not at all.** The 4.18 / 0.74 ratios are the 200 mg surgical
  cohort; the 425 mg cohort gives 7.30 / 2.62. The MSI concentration range is 970 nM–7,270 nM
  across 10 specimens from 6 patients. Getting these wrong to the people who measured them
  ends the conversation.
- **The affiliation sentence must be honest and chosen deliberately.** A DUA needs a
  counterparty legal entity. Saying so plainly is better received than discovering it in month
  three.
- **The pre-registration this email promises does not exist yet.** There is no pre-registration
  file in `results/`. Write it before sending, or delete the offer. Promising a timestamped
  prediction and then producing it after seeing the data is the one thing that would make this
  worse than not writing at all.
- **Consider co-sending with Forest White (MIT)** on the TMT design-map ask specifically — it
  is his lab's arm of the 2022 paper. Verify his corresponding-author status first.
- **Sequencing.** `../DATA_ACCESS_PLAN.md` §4.3 recommends sending this *after* an institutional
  collaboration exists (Target 1). If sent cold anyway, treat gate G2 (8 weeks, no reply → stop)
  as binding rather than following up repeatedly.
- **Same institution, second asset.** The Jonas Lab / Kibur Medical implantable-microdevice
  glioma study (Peruzzi & Jonas, BWH; PMID 37672566) is the *other* glioma dataset that pairs
  per-site drug concentration with per-site response and is not machine-readable
  (`results/zenodo_kibur_deposit_check.md`). Both live at BWH. **Do not bundle them into one
  email** — but if a BWH collaboration ever opens, it is the obvious second ask.
