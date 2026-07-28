# DRAFT — NOT SENT. FOR HUMAN REVIEW, EDITING AND SENDING.

> **Status: unsent draft.** No contact has been made with the depositors of GSE324268. Fill
> every `[BRACKET]`. **The submitter's name and email must be read off the GEO record itself
> before sending** — see "Target" below; they have not been verified here.
>
> **Priority 2** of the access plan, and the smallest, highest-probability ask in the set.

---

## Why this one is worth sending

`results/gbm_tmz_exvivo.md` established that MGMT — the determinant of temozolomide
sensitivity — is on **neither** Xenium panel we use (319-gene brain, 480-gene IO; 739-gene
union, MGMT absent), that **no panel gene surrogates for it** (0/396 pass both FDR bars; the
panel's MGMT-correlation base rate is 14.1% against 12.3% genome-wide, i.e. no enrichment at
all), and that GDSC cannot resolve it either (ρ = −0.068, p = 0.71). The conclusion was that
the uniform-TMZ-sensitivity assumption is currently unavoidable, at a measured cost of ~27% of
the between-model rank variance in ex-vivo response.

**GSE324268 changes that, and it is a free download.** Verified from the deposit's own panel
JSON and expression matrix:

- Panel = Xenium Human Brain **v1.1** (266 genes) **+ a custom 100-gene add-on named
  "Glioblastoma Panel"** (design_id `ECQNGH`, `hBrain_100g`) = **366 targets**.
- **MGMT (`ENSG00000170430`) is one of the 100 custom genes, with 8 probe sets.**
- `GSM9572843` `GBM_tumor07` — primary resected human GBM, **851,509 cells**, with **MGMT =
  15,191 counts across 14,675 cells (1.7%)**. For scale in the same sample: CXCL12 = 18,856
  counts, EGFR = 13.4M, GFAP = 15.7M.
- `GSM9572842` — a patient-derived GBM organoid on a brain slice.
- Full Xenium outputs per sample (transcripts / cells / boundaries parquet,
  `cell_feature_matrix.h5`, morphology OME-TIFF). **Free direct download from NCBI FTP, no
  login, no DUA.**
- The add-on also carries ATRX, CDKN2A, CDK4, CA9 and CXCL12 — several of which are directly
  useful to us.

So **no email is needed to get the data.** The email is for the two things the deposit does not
contain, both cheap for them to answer.

---

## Target

<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE324268> — Hölzel lab, Bonn (public
March 2026).

**Read the submitter name and contact email off the GEO record before sending.** They are
listed there. **Do not guess an address, and do not assume the lab head is the submitter.**
Cross-check against the associated publication if one is linked by the time you send — if the
paper is out, the corresponding author is the better target and the paper should be cited in
the first line.

---

## Subject line

`GSE324268 — MGMT in the Xenium GBM add-on panel: two short questions`

---

## Draft message

Dear [NAME],

Thank you for depositing GSE324268 with the full Xenium outputs and the panel JSON — the custom
100-gene add-on is, as far as I can find, the only public human glioblastoma spatial dataset
that measures **MGMT** with dedicated probes, and it is directly useful to me.

Briefly, why: I maintain an open-source mechanistic model of drug transport through
patient-specific glioblastoma tissue, reconstructed from spatial transcriptomics. It is a
research prototype — not clinically validated, and nothing in it should inform a treatment
decision. Its delivery predictions are scored against 125 human glioma tissue samples across 42
patients (roughly 2× geometric-mean fold error on tumour:plasma partition), but its
temozolomide predictions carry a known, quantified hole: MGMT is absent from every standard
Xenium panel, no panel gene surrogates for it (I screened all 396 panel genes measurable in
GSE232173 — zero pass an FDR bar against both MGMT and measured ex-vivo response), and GDSC's
TMZ arm has no dynamic range to resolve it. So per-cell TMZ sensitivity in my model is
currently a constant, and I have measured what that costs: about 27% of the between-model rank
variance in measured ex-vivo TMZ response, across an 81-percentage-point viability spread.

Your panel is the first thing I have found that could close that. Two questions, both small:

1. **Is the MGMT promoter-methylation status of `GBM_tumor07` (GSM9572843) known, and could you
   share it?** My analysis of the Erasmus MC ex-vivo cohort (Ntafoulis et al., PMID 37620410)
   suggests MGMT *expression* is the better predictor — in that cohort expression was prognostic
   for OS at HR 1.85, p = 6.1e-05, while promoter methylation was not significant at p = 0.16 —
   but having both on one specimen would let me check the expression↔methylation relationship
   directly on tissue where I also have 851,509 cells of spatial context. Any other routine
   molecular annotation for that specimen (IDH, 1p/19q, EGFR amplification) would be equally
   welcome, as would the WHO grade and whether it was primary or recurrent.

2. **How did you choose the 100 add-on genes, and would you be willing to share the rationale
   or the design file?** I am weighing whether to specify a custom panel myself, and the
   reasoning behind an existing GBM-specific add-on is worth more to me than the gene list
   alone. I would also be interested in whether the MGMT probes behaved as expected — at 1.7%
   per-cell detection in your data I suspect MGMT is only honestly usable aggregated over
   regions or populations rather than per cell, and I would value your view before I build
   anything on it.

If it is useful, I am happy to share back whatever I derive from the section — the
reconstruction (typed cells, vascular field, ECM and oxygen on your coordinates) and the
transport solve are open source, and I would send you the outputs and the code regardless of
whether the result is favourable. No obligation either way; the deposit as it stands is already
more than most.

Best regards,

[YOUR FULL NAME]
[AFFILIATION, or: "Independent researcher — no institutional affiliation"]
[EMAIL]
[REPOSITORY / PROJECT LINK, if public]

---

## Reviewer notes

- **Download and reproduce the numbers before sending.** The 851,509 cells and 15,191 MGMT
  counts were verified by an automated check, not by hand. Quoting a number back to the people
  who generated it is a good way to demonstrate seriousness and a bad way to be wrong. Confirm
  them locally first.
- **Check whether the associated paper has published** by the time you send. If it has, cite it
  in the first sentence and write to the corresponding author instead.
- **Do not ask for tissue or for a new run.** The whole value of this contact is that it is a
  two-question, five-minute reply. Keep it that way.
- **Do not overstate the "only public dataset" claim.** GSE305393 (CosMx 6K, human glioma TMA,
  18 patients) also carries MGMT, and the Visium HD human brain cancer dataset carries it in an
  anaplastic oligodendroglioma. The accurate claim is: the only public **Xenium** human **GBM**
  dataset with dedicated MGMT probes. "As far as I can find" is doing real work in that
  sentence.
- **Carry the caveat.** Every one of these platforms measures MGMT transcript abundance, not
  promoter methylation, which is the clinically actionable biomarker. Never let that slip out
  of a claim.
