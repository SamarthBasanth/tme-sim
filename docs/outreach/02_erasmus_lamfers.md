# DRAFT — NOT SENT. FOR HUMAN REVIEW, EDITING AND SENDING.

> **Status: unsent draft.** No contact has been made with Dr Lamfers, the GLIOscreen group,
> or anyone at Erasmus MC. They do not know this project exists and have agreed to nothing.
> Fill every `[BRACKET]`. Do not invent an affiliation. Verify the address before sending.
>
> **Priority 1** of the access plan — this is the decisive experiment
> (`results/ranking_robustness.md` §8, item 1).
>
> **Fact-check pass 2026-08-12.** Two substantive corrections: a **retracted** validation
> number was being quoted as current, and the cisplatin/doxorubicin motivation was two
> revisions out of date. See "What changed".

---

## Target

**Dr Martine L. M. Lamfers** — PI, GLIOscreen, Department of Neurosurgery / Brain Tumor Center,
Erasmus MC Cancer Institute, Rotterdam.

- Address: **m.lamfers@erasmusmc.nl** — read directly off the PubMed affiliation record for
  PMID 37620410 (re-verified 2026-08-12), and consistent with her Erasmus MC group page and
  with PMID 40425813. **Re-verify on the day you send.**
- Group page: <https://www.erasmusmc.nl/en/research/groups/glioscreen>

**Do not** address correspondence about this to Ioannis Ntafoulis at Erasmus MC — public
profiles indicate he has since moved to CHDR Leiden. He is the first author and worth citing;
he is not the contact.

**Consider copying, but only if the first reply is positive** — copying four people on a cold
email reads as a mass mailing:

- **Prof. Johan M. Kros** (Pathology & Clinical Bioinformatics) — a co-author on PMID 37620410
  (initials **JM**, verified against the PubMed author list 2026-08-12) *and* the senior figure
  on their GeoMx DSP brain work with Dr Dana Mustafa. The single most natural bridge between
  the drug screen and a spatial assay.
- **Dr Miao-Ping Chien** (Molecular Genetics / Oncode) — a co-author on their 2025 GSC paper
  (PMID 40425813, verified) and a spatial single-cell methods developer.
- **Prof. Sieger Leenstra** — senior/surgical PI on the group (co-author on both papers).

---

## Subject line options

- `Reanalysis of GSE232173, and one experiment it points to`
- `Ex-vivo drug sensitivity + spatial profiling on the same GBM specimen — a question`

The first is probably best: it leads with something we did with *their* data, not with a
request.

---

## Draft message

Dear Dr Lamfers,

I have been working on an open-source mechanistic model of drug transport through
patient-specific glioblastoma tissue, and I have spent the last few weeks scoring it against
your 2023 *British Journal of Cancer* ex-vivo screen (Ntafoulis et al., PMID 37620410,
GSE232173). That reanalysis pointed at one experiment your group is, as far as I can find,
uniquely placed to run, and I would like to ask whether it is of any interest.

**What I did with your data, first, so this is concrete.** I took our potency layer's
temozolomide IC50 — extrapolated from GDSC2 across 51 glioma lines, with the honest caveat that
all 51 fitted IC50s sit above GDSC's 30 µM assay ceiling — recovered each line's Hill slope
from its own reported (IC50, AUC) pair, and evaluated at 100 µM. Zero parameters fitted to your
data. Predicted median viability 78.1%; your 40 cultures measured 82.4% (Mann-Whitney p = 0.27).
Separately, I re-derived the MGMT association from the raw GSE232173 counts and got
Spearman ρ = +0.516 (p = 0.034, n = 17 with a published viability number and paired RNA-seq),
rising to +0.735 on the 15 cultures in your tissue-paired cohort — an independent replication
of your ρ = +0.71 through a different pipeline. I also found and flagged the labelling
discrepancy in MOESM2 Supplemental Table 5 (the PFS rows read `UM 8.55 / M 5.92`, which
contradicts the Results text) rather than silently resolving it. The write-up is at
[`results/gbm_tmz_exvivo.md`, LINK].

**Where the model stands, plainly.** It is a research prototype. It is not clinically
validated, it is not a device, and nothing in it should inform a treatment decision. An
adversarial audit of my own work found **eleven distinct defects**, several mutually cancelling
and invisible to a passing test suite. Two components were **falsified and retracted** — the
oxygen PDE, beaten by plain distance-to-nearest-vessel, and a penetration metric that turned
out to be definitionally unable to discriminate between drugs. A third claim, cross-platform
agreement of the tissue reconstruction, was **also retracted** after a classifier fix: the
number that mattered collapsed and that comparison is now an open question rather than a pass.
What survives is the delivery axis: scored against 125 human glioma tissue samples across 42
patients plus 27 PMID-anchored regional PK measurements, the model predicts tumour:plasma
partition to about a **2× geometric-mean fold error**, with 6/6 drugs inside 10×. The pattern
is consistent and not flattering — every transport-level validation has held, and every
reconstruction-level one has been retracted.

**What I cannot do, and why I am writing.** The response half has never been tested. I ran a
systematic sweep of 164 datasets, and there is no open dataset anywhere — in any tumour type —
that pairs a measured spatial drug distribution with a measured treatment response in the same
tumour. More specifically to your field: as far as I can establish, **no published study in
glioma has run an ex-vivo drug sensitivity screen and a spatial molecular assay on the same
resection.** The closest is Shekarian et al. (*Sci Adv* 2022, PMID 35776791), which is CODEX on
treated explants rather than a dose-response screen. Yu et al. (*Nat Commun* 2025,
PMID 40295505) took stereotactic GBM cores through spatial transcriptomics, t-CyCIF and
MALDI-MSI and derived models from the same biopsy — and then wrote, in their own discussion,
that having models matched to multi-omics from the same biopsy "would allow for functional
integration with multi-omics results." They did not do it.

**The ask.** For one GBM resection already entering your GSC pipeline, would it be feasible to:

1. retain an adjacent FFPE or OCT block from the same tumour region;
2. add **cisplatin** and **doxorubicin** as full dose-response curves to whatever panel you run
   on the derived culture (I do not know whether they are already among the 20 agents screened
   in Posthoorn-Verheul et al., PMID 40425813 — if they are, the incremental cost is nil); and
3. run any spatial assay on the retained block — Xenium, Visium, GeoMx or CosMx all work for
   me. I understand this capability exists within Erasmus MC, and mention it only because it
   means the question is about pairing two things you already do rather than acquiring
   something new.

Those two drugs are not arbitrary. Our panel ranking is decided by delivery × potency, and
cisplatin and doxorubicin sit at opposite ends of both axes — cisplatin has roughly 157× the
delivered exposure and 49× the IC50. Propagating our reconstruction, calibration and potency
priors through 10,240 complete re-rankings puts cisplatin first with probability 0.931 on the
transport law we ship. But when we re-run it with the physically correct gradient-driven
vascular source (`S = PS(Cp − C)` rather than a one-way leak), that becomes 0.397 and
**doxorubicin takes rank 1 at 0.559** — the order reverses. So the reportable output today is a
set, not a winner, and the ambiguity has two separable halves. The delivery half needs an
absolute permeability measurement I do not have. The potency half is collapsible: our own
sensitivity analysis says the dominant remaining term is "the median glioma cell line" standing
in for "this tumour", and a per-patient IC50 pair for those two drugs replaces a population
spread with a measurement error. Doxorubicin's prior is the widest of any drug on our panel for
a mundane reason — it is simply absent from the GDSC glioma screen.

A retrospective version would also be valuable and cheaper: the two extra curves on cultures
derived from specimens for which a spatial section already exists in your archive.

**What I would provide.** Before seeing any of your data I would commit, in a timestamped and
version-controlled file, to the predicted drug ordering, the predicted per-region exposure, and
the predicted MGMT effect size — so the comparison is a prediction rather than a fit. I would
run the reconstruction and transport model on your section at no cost, and hand over the code
and the validation harness (permutation tests, leave-one-out cross-validation, an explicit
uncertainty Monte Carlo) whether or not the result is favourable. I will **publish the result
if the model fails**; I would rather agree that now than have it become awkward later.

On authorship: entirely as you prefer, including a supplementary analysis, an acknowledgement,
or nothing at all. My interest is in whether the model is right.

**Administratively**, my understanding is that secondary use of coded residual tissue is not
WMO-liable and would go through your NWBTC and an MDTA rather than a full METC submission, and
that in practice I would be receiving derived tables — IC50s and a spatial expression matrix —
rather than any material. I would follow whatever route you specify. [**IF APPLICABLE — DELETE
IF NOT:** I should flag now that this work [is / may become] commercial, which I understand
routes through the Erasmus MC TTO rather than the standard non-commercial MTA.]

If this is not something you have capacity for, a one-line no is genuinely useful and I will
not follow up. If it is interesting but the framing is wrong, I would welcome being told so.

With thanks for making GSE232173 and the supplementary tables public — the reanalysis above was
only possible because you did.

Best regards,

[YOUR FULL NAME]
[AFFILIATION, or: "Independent researcher — no institutional affiliation"]
[EMAIL]
[REPOSITORY / PROJECT LINK, if public]

---

## What changed in the 2026-08-12 fact-check

| was | now | why |
|---|---|---|
| "The reconstruction agrees across two independent spatial platforms at **Pearson 0.997**." | Sentence deleted; replaced with an explicit statement that this claim was **retracted** | `results/revalidation.md` §223: "Anything downstream quoting 'Pearson 0.997 cross-platform concordance' is now [wrong]." After the classifier-arity fix the 4-vector fell to 0.980 and the non-tumour subvector — the one that mattered — collapsed to **0.016**. The 0.997 was agreement between two biased estimators. This was the single most serious error in the outreach set: a retracted number, presented as a headline validation, in an email going out over the sender's name. |
| "reproduces the measured BBB-penetrant vs BBB-excluded separation (3.65× predicted, 4.66× measured) from molecular structure alone" | **removed** | Pre-fix value. The current model predicts 17.20× against measured 4.66× pooled — an over-separation, not a match (`results/glioma_regional_pk.md` §P2). |
| "They currently differ by about 4% in predicted kill … the honest output today is a tie" | Replaced with the current position: P(cisplatin rank 1) = 0.931 on the shipped law, **0.397 with doxorubicin at 0.559** under the corrected gradient-driven source | `results/ranking_robustness.md` §0 and its conditional-retraction notice. The 4% figure predates the discrete-vessel promotion *and* the delivery-law correction. The new framing is both current and a better argument: the ask is motivated by a stated order reversal, not by a small margin. |
| "the 20-agent panel of Posthoorn-Verheul et al." | "among the 20 agents screened in Posthoorn-Verheul et al." | Verified against the abstract: "enabled the screening of 20 anti-cancer agents … for 16 out of 18 HGG samples." It is a screen of 20 agents, not necessarily a fixed named panel. |
| "I am aware this capability exists within Erasmus MC, including in Prof. Kros's group" | "I understand this capability exists within Erasmus MC" | The Kros/GeoMx mention stays in the cc-list rationale where it belongs; naming a colleague's lab in the body reads as telling a PI how to run her institution. |

**Verified unchanged and correct:** PMID 37620410 (Ntafoulis I, Kleijn A, Ju J, …, Lamfers MLM;
*Br J Cancer* 2023;129(8):1327–1338); the corresponding address `m.lamfers@erasmusmc.nl` read
off the record itself; Kros **JM**; Chien MP on PMID 40425813; PMID 35776791 (Shekarian T et al.,
*Sci Adv* 2022;8(26):eabn9440, Hutter G senior, Basel); PMID 40295505 (Yu KKH et al.,
*Nat Commun* 2025;16:3957) and the quoted sentence, confirmed verbatim: "Having models that are
matched to multi-omics data from the same needle biopsy would allow for functional integration
with multi-omics results and further mechanistic insights into resistance."; all GSE232173
reanalysis figures (78.1 / 82.4 / p = 0.27 / ρ = +0.516 / +0.735 / MOESM2 `UM 8.55 / M 5.92`);
157× and 49×.

## Reviewer notes

- **Length.** Still long for a cold email. If cutting, keep: the GSE232173 reanalysis numbers
  (paragraph 2), the three-point ask, the order-reversal paragraph, and the pre-registration
  offer. Cut the audit paragraph to one sentence and link out.
- **The pre-registration this email promises does not exist yet.** There is no pre-registration
  file in `results/`. Write it before sending, or delete the offer.
- **Attachments.** Attach nothing on the first email. Link to one page if the repo is public.
- **Do not claim** the specific instrument at PARTS. Reporting of a 2024 10x spatial
  acquisition is secondary-source and unverified; the PARTS page itself lists no instruments.
- **Do not overstate the gap claim.** "As far as I can establish" and "as far as I can find" are
  load-bearing. Full-text search indexes open access only; a paywalled Methods section could
  contain a counterexample.
- **The commercial disclosure is not optional** if it applies. See the biobank regulation note
  in `../DATA_ACCESS_PLAN.md` §2.4.
- **Fallbacks if there is no reply in 6 weeks** (gate G1): Weller/Snijder at Zurich, whose
  Tumor Profiler programme splits every biopsy between imaging mass cytometry and pharmacoscopy
  by design and has the exact capability, just never pointed at glioma — **Weller M and
  Snijder B are both verified consortium authors** on PMID 40425842 (*Nat Med* 2025;31(7):2430)
  and Snijder B on PMID 33482122 (*Cancer Cell* 2021;39(3):288); Golebiewska at LIH; Hutter at
  Basel (PMID 35776791).
- **A second candidate has appeared since this draft was written.** The Bonn group behind
  GSE324268 published *Core2Edge* on 2026-07-28 (PMID 42521807, *Nat Protoc*): human GBM
  organoids implanted into **human** organotypic brain slices, explicitly offered as a
  drug-screening platform, from a group that has already deposited a Xenium section of a
  primary GBM. It is ex vivo and therefore has no vasculature — it cannot test perfusion-driven
  delivery — but it is a real second address for the same question. See
  `04_holzel_gse324268.md`.
