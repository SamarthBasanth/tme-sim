# DRAFT — NOT SENT. FOR HUMAN REVIEW, EDITING AND SENDING.

> **Status: unsent draft.** No contact has been made with Dr Lamfers, the GLIOscreen group,
> or anyone at Erasmus MC. They do not know this project exists and have agreed to nothing.
> Fill every `[BRACKET]`. Do not invent an affiliation. Verify the address before sending.
>
> **Priority 1** of the access plan — this is the decisive experiment
> (`results/ranking_robustness.md` §8).

---

## Target

**Dr Martine L. M. Lamfers** — PI, GLIOscreen, Department of Neurosurgery / Brain Tumor Center,
Erasmus MC Cancer Institute, Rotterdam.

- Address: **m.lamfers@erasmusmc.nl** — the corresponding-author address printed on
  PMID 37620410, and consistent with her Erasmus MC group page and PMID 40425813.
  **Re-verify on the day you send.**
- Group page: <https://www.erasmusmc.nl/en/research/groups/glioscreen>

**Do not** address correspondence about this to Ioannis Ntafoulis at Erasmus MC — public
profiles indicate he has since moved to CHDR Leiden. He is the first author and worth citing;
he is not the contact.

**Consider copying, but only if the first reply is positive** — copying four people on a cold
email reads as a mass mailing:

- **Prof. Johan M. Kros** (Pathology & Clinical Bioinformatics) — a co-author on PMID 37620410
  *and* the senior figure on their GeoMx DSP brain work with Dr Dana Mustafa. The single most
  natural bridge between the drug screen and a spatial assay.
- **Dr Miao-Ping Chien** (Molecular Genetics / Oncode) — a co-author on their own 2025 GSC
  paper and a spatial single-cell methods developer.
- **Prof. Sieger Leenstra** — senior/surgical PI on the group.

---

## Subject line options

- `Ex-vivo drug sensitivity + spatial profiling on the same GBM specimen — a question`
- `Reanalysis of GSE232173, and one experiment it points to`
- `A spatial transport model scored against your TMZ ex-vivo data`

The second is probably best: it leads with something we did with *their* data, not with a
request.

---

## Draft message

Dear Dr Lamfers,

I have been working on an open-source mechanistic model of drug transport through
patient-specific glioblastoma tissue, and I have spent the last few weeks scoring it against
your 2023 *British Journal of Cancer* ex-vivo screen (Ntafoulis et al., PMID 37620410,
GSE232173). I am writing because that reanalysis pointed clearly at one experiment your group
is, as far as I can find, uniquely placed to run — and I would like to ask whether it is of any
interest.

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

**Where the model actually stands, stated plainly.** It is a research prototype. It is not
clinically validated, it is not a device, and nothing in it should inform a treatment decision.
An adversarial audit of my own work found **eleven distinct defects**, several of them mutually
cancelling and invisible to a passing test suite; every headline number I had published before
that audit was wrong in at least one way. One component was **falsified and retracted**: a
penetration metric that turned out to be definitionally unable to discriminate between drugs. A
second retraction has since been **withdrawn in its turn** — I had demoted my own oxygen PDE for
losing to plain distance-to-nearest-vessel, and that comparison was an inverted inequality. The
predictors are sign-aligned to oxygen, so |−0.1523| > |−0.1204| means the solve out-predicts the
baseline, by 1.26×. I mention it because it is the most instructive failure in the project: a
component was deprecated on a misread inequality while every mechanical check in the same
repository disagreed with the prose. The corrected reading stays modest — at ρ ≈ 0.15,
within-section, it is the best of four weak predictors rather than a validated oxygen field, and
its diffusion length remains 2–5× the published value because the consumption term has never been
fitted to a measured tension. What survives is the
delivery axis: scored against 125 human glioma tissue samples across 42 patients plus 27
PMID-anchored regional PK measurements, the model predicts tumour:plasma partition to about a
**2× geometric-mean fold error**, with 6/6 drugs inside 10×, and reproduces the measured
BBB-penetrant vs BBB-excluded separation (3.65× predicted, 4.66× measured) from molecular
structure alone. The reconstruction agrees across two independent spatial platforms at
Pearson 0.997.

**What I cannot do, and why I am writing.** The response half has never been tested. I ran a
systematic sweep of 164 datasets, and there is no open dataset anywhere — in any tumour type —
that pairs a measured spatial drug distribution with a measured treatment response in the same
tumour. More specifically to your field: as far as I can establish, **no published study in
glioma has run an ex-vivo drug sensitivity screen and a spatial molecular assay on the same
resection.** The closest is Shekarian et al. (*Sci Adv* 2022, PMID 35776791), which is CODEX on
treated explants rather than a dose-response screen, and two 2025 reviews cite it as the only
example of the design. Yu et al. (*Nat Commun* 2025, PMID 40295505) took stereotactic GBM cores
through spatial transcriptomics, t-CyCIF and MALDI-MSI and derived models from the same biopsy
— and then wrote, in their own discussion, that having models matched to multi-omics from the
same biopsy "would allow for functional integration with multi-omics results." They did not do
it.

**The ask.** For one GBM resection already entering your GSC pipeline, would it be feasible to:

1. retain an adjacent FFPE or OCT block from the same tumour region;
2. add **cisplatin** and **doxorubicin** as full dose-response curves to whatever panel you run
   on the derived culture (I do not know whether they are already in the 20-agent panel of
   Posthoorn-Verheul et al., PMID 40425813 — if they are, the incremental cost is nil); and
3. run any spatial assay on the retained block — Xenium, Visium, GeoMx or CosMx all work for
   me. I am aware this capability exists within Erasmus MC, including in Prof. Kros's group,
   and I mention it only because it means the question is about pairing two things you already
   do rather than acquiring something new.

Those two drugs are not arbitrary. Our panel ranking is decided by delivery × potency, and
cisplatin and doxorubicin sit at opposite ends of both axes — cisplatin has roughly 157× the
delivered exposure and 49× the IC50. They currently differ by about 4% in predicted kill, which
is well inside our own parameter uncertainty, so the honest output today is a tie rather than a
recommendation. A Monte Carlo over that uncertainty shows the dominant remaining term is
"the median glioma cell line" standing in for "this tumour." A per-patient IC50 pair for those
two collapses it. Doxorubicin's uncertainty is particularly inflated because it is simply
absent from the GDSC glioma screen.

A retrospective version would also be valuable and cheaper: the two extra curves on cultures
derived from specimens for which a spatial section already exists in your archive.

**What I would provide.** Before seeing any of your data I would commit, in a timestamped and
version-controlled file, to the predicted drug ordering, the predicted per-region exposure, and
the predicted MGMT effect size — so the comparison is a prediction rather than a fit. I would
run the reconstruction and transport model on your section at no cost, and hand over the code
and the validation harness (a set of scoring modules with permutation tests, leave-one-out
cross-validation and an explicit uncertainty Monte Carlo) whether or not the result is
favourable. I would like to say explicitly that **I will publish the result if the model
fails**; I would rather agree that now than have it become awkward later.

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

## Reviewer notes

- **Length.** This is long for a cold email. If cutting, keep: the GSE232173 reanalysis numbers
  (paragraph 2), the three-point ask, and the pre-registration offer. Cut the audit paragraph
  down to one sentence and link out.
- **Attachments.** Attach nothing on the first email. Link to one page if the repo is public.
- **The Kros mention is deliberate but risky.** It shows homework and makes the ask concrete,
  but it can read as telling a PI how to run her institution. If it feels presumptuous, replace
  with: "I understand spatial profiling capability exists at Erasmus MC, though I do not know
  how it is organised."
- **Do not claim** the specific instrument at PARTS. Reporting of a 2024 10x spatial
  acquisition is secondary-source and unverified; the PARTS page itself lists no instruments.
- **Do not overstate the gap claim.** "As far as I can establish" and "as far as I can find" are
  load-bearing. Full-text search indexes open access only; a paywalled Methods section could
  contain a counterexample.
- **The commercial disclosure is not optional** if it applies. See the biobank regulation note
  in `../DATA_ACCESS_PLAN.md` §2.4.
- **Fallbacks if there is no reply in 6 weeks** (gate G1): Weller/Snijder at Zurich, whose
  Tumor Profiler programme (PMID 40425842, PMID 33482122) splits every biopsy between imaging
  mass cytometry and pharmacoscopy by design and has the exact capability, just never pointed
  at glioma; Golebiewska at LIH; Hutter at Basel (PMID 35776791).
