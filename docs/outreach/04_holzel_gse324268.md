# DRAFT — NOT SENT. FOR HUMAN REVIEW, EDITING AND SENDING.

> **Status: unsent draft.** No contact has been made with the depositors of GSE324268. Fill
> every `[BRACKET]`. **The submitter's name and contact email must be read off the GEO record
> itself before sending** — see "Target". They have not been verified here (NCBI is not
> reachable from the session that wrote this).
>
> **Priority 2** of the access plan, and the smallest, highest-probability ask in the set.
>
> **Fact-check and rewrite pass 2026-08-12.** The dataset has since been **downloaded,
> checksum-verified and analysed** (`results/mgmt_spatial.md`), which changes what this email
> should say: it now reports a finding about their panel rather than a suspicion. A candidate
> associated publication has also appeared. See "What changed".

---

## Why this one is worth sending

`results/gbm_tmz_exvivo.md` established that MGMT — the determinant of temozolomide
sensitivity — is on **neither** Xenium panel we use (319-gene brain, 480-gene IO; 739-gene
union, MGMT absent), that **no panel gene surrogates for it** (0/396 pass both FDR bars; the
panel's MGMT-correlation base rate is 14.1% against 12.3% genome-wide, i.e. no enrichment at
all), and that GDSC cannot resolve it either (ρ = −0.068, p = 0.71). The cost of the resulting
uniform-TMZ-sensitivity assumption was measured at ~27% of the between-model rank variance in
ex-vivo response.

**GSE324268 was the free download that let that be tested, and it has now been tested.**
Verified from the deposit's own panel JSON and expression matrix, with sha256s recorded in
`results/mgmt_spatial.md` §1:

- Panel = Xenium Human Brain **hBrain_v1.1** (266 targets) **+ a custom 100-target add-on named
  "Glioblastoma Panel"** (design_id `ECQNGH`, `hBrain_100g`, designed 2024-01-08). 400 target
  entries in the spec; **366 appear as `Gene Expression` features** in the shipped matrix.
- **MGMT (`ENSG00000170430`) is one of the add-on targets**, `gene_coverage = 8` — eight probe
  pairs on **one** gene target, i.e. one column in the count matrix, not eight features.
- `GSM9572843` `GBMtumor07` — primary resected human GBM, **851,509 cells**, with **MGMT =
  15,191 counts across 14,675 cells (1.72%)**. Both figures reproduce exactly from the HDF5.
- `GSM9572842` `GBM1organoid` — a patient-derived GBM organoid on a brain slice.
- Full Xenium outputs per sample (transcripts / cells / boundaries parquet,
  `cell_feature_matrix.h5`, morphology OME-TIFF). **Free direct download from NCBI FTP, no
  login, no DUA.**

**And the analysis produced a result worth telling them about.** MGMT is genuinely detected —
6.6× the mean negative-control probe on a per-probe basis — but it is **not spatially
resolvable at this depth**, and the deposit's own controls are what say so:

- max 4 MGMT counts in any single cell; 1.19% of tumour cells positive. No per-cell MGMT call
  exists at this depth and none was made.
- Moran's I = +0.00152 (z = +2.31, p = 0.011) — statistically non-zero, but **46 of 47
  abundance-matched genes are more clustered**.
- At 400 µm bins, MGMT dispersion φ = 1.382 against the pooled negative-control probes' 1.374;
  **thinned to the controls' exact count depth MGMT falls to φ = 1.211 ± 0.064, below the
  controls**.
- This is not an underpowered null: the same test on the same cells returns φ up to 154 and
  Moran's I up to +0.175 for real genes at the same abundance. The assay can see structure
  here. MGMT has none to see.
- Consequence for us: a per-specimen MGMT **scalar** is not a limitation to be apologised for,
  it is the empirically correct object. Swapping the scalar for spatially-resolved MGMT moves
  predicted TMZ kill by at most 0.69 percentage points, and a counts-shuffled spatial null
  reproduces 46–89% of that shift.

So **no email is needed to get the data, and none is needed to use it.** The email is for the
one thing the deposit does not contain and cannot contain — a promoter-methylation call — plus
one question about panel design, and an offer to hand them a QC observation about their own
probes.

---

## Target

<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE324268> — series title, as recorded by
the loader that downloaded it: *"Spatial transcriptomics of human glioblastoma and of a
GBM-organoid-on-brain-slice model"*. Hölzel lab, Bonn; public 2026-03.

**Read the submitter name and contact email off the GEO record before sending. Do not guess an
address, and do not assume the lab head is the submitter.**

**Candidate associated publication, found 2026-08-12 — check before you write.** A Bonn group
published *Core2Edge: a human glioblastoma organoid–brain slice model capturing infiltration and
transcriptional heterogeneity from core to single-cell dispersion*, **Nat Protoc**, online
**2026-07-28**, **PMID 42521807**, doi:10.1038/s41596-026-01412-3. Melhem A, Pregler BEF, …,
**Hölzel M**, …, Potthoff AL, **Schneider M**. The corresponding address printed in the PubMed
affiliation record is `matthias.schneider@ukbonn.de` (Department of Neurosurgery / Brain Tumor
Translational Research Group, University Hospital Bonn).

**This is a lead, not a confirmed link.** It is the same institution, the same lab head among
the authors, and the same organoid-on-brain-slice model that `GSM9572842` describes — but
whether GEO lists it as the series' citation was **not verifiable from the session that wrote
this note**, because NCBI was unreachable. Open the GEO record. If the paper is linked there,
cite it in the first line and write to whoever GEO names, cross-checked against the paper's
corresponding author.

---

## Subject line

`GSE324268 — MGMT in the Xenium GBM add-on panel: two questions and one observation`

---

## Draft message

Dear [NAME],

Thank you for depositing GSE324268 with the full Xenium outputs and the panel JSON. The custom
100-target add-on is, as far as I can find, the only public human glioblastoma **Xenium**
dataset with dedicated MGMT probes, and I have now downloaded and analysed the `GBMtumor07`
section. Two short questions, and one observation about your probes that may be worth having.

Briefly, why: I maintain an open-source mechanistic model of drug transport through
patient-specific glioblastoma tissue, reconstructed from spatial transcriptomics. It is a
research prototype — not clinically validated, and nothing in it should inform a treatment
decision. Its delivery predictions are scored against 125 human glioma tissue samples across 42
patients (roughly 2× geometric-mean fold error on tumour:plasma partition), but its
temozolomide predictions carried a known hole: MGMT is absent from every standard Xenium panel,
no panel gene surrogates for it (I screened all 396 panel genes measurable in GSE232173 — zero
pass an FDR bar against both MGMT and measured ex-vivo response), and GDSC's TMZ arm has no
dynamic range to resolve it. Your section is the first thing I found that could test whether
that hole is closable spatially.

**The observation.** MGMT is genuinely detected in `GBMtumor07` — 15,191 counts across 14,675
of 851,509 cells, which is 6.6× the mean negative-control probe on a per-probe basis. But I
could not make it spatially usable, and your deposit's own controls are what settled it. Moran's
I is +0.00152 (p = 0.011) with 46 of 47 abundance-matched genes more clustered; at 400 µm bins
MGMT's dispersion is 1.382 against 1.374 for the pooled negative-control probes, and thinned to
the controls' exact count depth it falls to 1.211 ± 0.064 — below them. The test is not
underpowered: real genes at the same abundance in the same cells return dispersions up to 154.
So the honest reading is that MGMT here is a per-specimen scalar, not a field, and I have kept
it as one. I mention this because it is a statement about the probes rather than about the
biology, and if you are designing further add-ons it may be worth knowing before someone builds
a per-cell MGMT call on eight probe pairs.

**Question 1 — is the MGMT promoter-methylation status of `GBMtumor07` (GSM9572843) known, and
could you share it?** Given the above, a per-specimen scalar is exactly the object I need. My
analysis of the Erasmus MC ex-vivo cohort (Ntafoulis et al., PMID 37620410) suggests MGMT
*expression* is the better predictor — in that cohort expression was prognostic for OS at
HR 1.85, p = 6.1e-05, while promoter methylation was not significant at p = 0.16 — but having
both on one specimen would let me check the expression↔methylation relationship directly on
tissue where I also have 851,509 cells of spatial context. Any other routine molecular
annotation for that specimen (IDH, 1p/19q, EGFR amplification) would be equally welcome, as
would the WHO grade and whether it was primary or recurrent.

**Question 2 — how did you choose the 100 add-on targets, and would you be willing to share the
rationale or the design file?** I am weighing whether to specify a custom panel myself, and the
reasoning behind an existing GBM-specific add-on is worth more to me than the gene list alone.
I would also value your view on whether the MGMT probes behaved as you expected.

If it is useful, I am happy to share back whatever I derive from the section — the
reconstruction (typed cells, vascular field, ECM and oxygen on your coordinates), the transport
solve, and the MGMT analysis above are all open source, and I would send you the outputs and
the code regardless of whether the result is favourable. No obligation either way; the deposit
as it stands is already more than most.

Best regards,

[YOUR FULL NAME]
[AFFILIATION, or: "Independent researcher — no institutional affiliation"]
[EMAIL]
[REPOSITORY / PROJECT LINK, if public]

---

## What changed in the 2026-08-12 pass

| was | now | why |
|---|---|---|
| Sample called `GBM_tumor07` | **`GBMtumor07`** (and `GBM1organoid` for GSM9572842) | The actual sample titles, read off the downloaded files. Quoting a sample name back wrong to the people who named it is a small error with a large signal. |
| "MGMT … with 8 probe sets" | "`gene_coverage = 8` — eight probe **pairs** on **one** gene target, one column in the matrix" | `results/mgmt_spatial.md` §1 records this as an explicit correction to the earlier brief. "8 probe sets" implies 8 features. |
| "Panel = … v1.1 (266 genes) + 100 = **366 targets**" | 400 target entries in the spec, **366** `Gene Expression` features in the matrix | Both numbers are real and they are not the same number. |
| "at 1.7% per-cell detection I **suspect** MGMT is only honestly usable aggregated over regions" | The measured result: not usable per cell **and not usable regionally either** — below the deposit's own negative controls at matched depth | The suspicion has been tested. Sending a suspicion when the answer is on disk wastes the recipient's attention; sending the answer is a contribution to their QC. |
| Reviewer note: "Download and reproduce the numbers before sending — verified by an automated check, not by hand" | Note removed; the download, sha256 verification and reproduction are done and cited | `results/mgmt_spatial.md` §1–§2. |
| "Check whether the associated paper has published" (open instruction) | A named candidate: **PMID 42521807**, *Nat Protoc*, 2026-07-28, with its corresponding address and an explicit warning that the GEO link is unverified | Published after this draft was written and after the 164-dataset sweep. |

**Verified unchanged:** every number in the "why this is worth sending" block against
`results/mgmt_spatial.md` and `results/gbm_tmz_exvivo.md` — 0/396, 14.1% vs 12.3%,
ρ = −0.068 / p = 0.71, ~27% rank variance, 81-point viability spread, HR 1.85 / p = 6.1e-05,
methylation p = 0.16 (n = 56), design_id `ECQNGH`, `ENSG00000170430`, 851,509 cells,
15,191 counts, 14,675 cells.

**One discrepancy left open.** The loader docstring
(`tmesim/data/real/mgmt_spatial.py`) describes `GSM9572842` as an organoid on a **mouse** brain
slice, while the series title says only "brain-slice model" and the Bonn group's *Core2Edge*
protocol is explicitly a **human** brain-slice system. Not resolvable without the GEO record.
It does not affect anything used here — that sample is deliberately not loaded — but **do not
describe the organoid sample's host species in the email.**

## Reviewer notes

- **Check whether the paper is linked on GEO** by the time you send. If it is, cite it in the
  first sentence and write to the corresponding author instead of the submitter.
- **Do not ask for tissue or for a new run.** The whole value of this contact is that it is a
  two-question, five-minute reply. Keep it that way.
- **Do not overstate the "only public dataset" claim.** GSE305393 (CosMx 6K, human glioma TMA,
  18 patients) also carries MGMT, and the Visium HD human brain cancer dataset carries it in an
  anaplastic oligodendroglioma. The accurate claim is: the only public **Xenium** human **GBM**
  dataset with dedicated MGMT probes. "As far as I can find" is doing real work there.
- **Carry the caveat.** Every one of these platforms measures MGMT transcript abundance, not
  promoter methylation, which is the clinically actionable biomarker. Never let that slip out
  of a claim.
- **The observation paragraph is the reason this email is worth their time.** If it is cut for
  length, cut question 2 instead. A finding about their own probes, offered for nothing, is
  what makes an unsolicited email from a stranger answerable.
- **A larger, later ask exists and is deliberately not made here.** *Core2Edge* (PMID 42521807)
  is a human GBM-organoid-in-human-brain-slice platform its authors present as a drug-screening
  system, from the group that deposited a Xenium section of a primary GBM. That is one
  institution holding both halves of the paired experiment `02_erasmus_lamfers.md` asks
  Rotterdam for — minus the vasculature, which is not a small minus for a transport model.
  **Do not raise it in the first email.**
