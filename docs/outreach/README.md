# Outreach drafts — NOT SENT, NOT APPROVED

> ## ⚠️ EVERY FILE IN THIS DIRECTORY IS AN UNSENT DRAFT FOR HUMAN REVIEW.
>
> Nothing here has been sent. No contact has been made with any person or group named
> in these files. Nobody named here has agreed to anything, been approached, or knows
> this project exists.
>
> **A human must read, edit, and personally send each message.** Do not automate this.

**Last fact-checked: 2026-08-12** — every citation, author list, institution, accession and
quoted result in these drafts was re-verified against the primary record (PubMed article
records, not search snippets; the repo's own results files; the downloaded bytes where they
exist). Corrections are logged in a "What changed" table at the foot of each draft. The drafts
are **send-ready and unsent**; what stands between them and being sent is a human decision and
the `[BRACKET]` fields, not further research.

## Before sending anything

1. **Fill every `[SQUARE BRACKET]` placeholder.** They mark facts only the sender knows —
   name, affiliation (or the honest absence of one), links, availability. **Do not invent an
   affiliation, a funder, a title, or a prior introduction.** "Independent researcher, no
   institutional affiliation" is a legitimate and sayable thing; a fabricated one ends the
   conversation permanently and correctly.
2. **Verify the contact details yourself the day you send.** Every address in these drafts
   came from a public source (a paper's corresponding-author line read off the record, an
   institutional faculty page) and is recorded with its source here and in
   `docs/DATA_ACCESS_PLAN.md`. Addresses rot, people move, and one of the lab websites
   referenced is currently down. Where no address is known, the draft says so and leaves a
   marked placeholder — **do not fill one in by guessing a format.**
3. **Write the pre-registration before you send.** Drafts 01 and 02 both promise a timestamped,
   version-controlled prediction made *before* any data arrives. **No such file exists in
   `results/` as of 2026-08-12.** Either write it first or delete the offer. Producing it after
   seeing the data is the one failure that would be worse than never writing.
4. **Re-read the honest-position paragraph in each draft and confirm it is still true.** It
   states: 11 defects found and fixed; delivery validated in human glioma to ~2× geometric-mean
   error; components falsified and retracted; response prediction unvalidatable with open data;
   research prototype, not clinically validated.
5. **Decide the commercial question before you send.** Erasmus MC's biobank regulation routes
   any commercial involvement through their TTO rather than the standard non-commercial MTA
   (§8.3.3). If this project has or may have a commercial dimension, say so in the first email.
6. **Send the short version.** Drafts 01 and 02 each carry a **Short version** (three or four
   paragraphs) followed by a **Long version**. The short one is the one to send; the long one
   exists because its numbers are what a reply will ask for. Drafts 03 and 04 are already
   short. A busy PI reads the first three lines: keep (a) what we built, (b) what we already
   did with *their* data, (c) the one specific thing we are asking for.

## Numbers that must never appear in an outgoing email

These were retracted or superseded inside the repo and were found being quoted as current
during the 2026-08-12 pass. Both had already reached a draft.

| do not write | why | what is true instead |
|---|---|---|
| "cross-platform reconstruction concordance, **Pearson 0.997**" | **Retracted.** `results/revalidation.md` — after the classifier-arity fix the 4-vector fell to 0.980 and the non-tumour subvector to **0.016**. The 0.997 was two biased estimators agreeing. | Cross-platform concordance is an **open question**. Every reconstruction-level external validation this project made has been retracted or fails on a held-out section; every transport-level one has held. Say that — it is more credible than the number was. |
| "penetrant-vs-excluded separation **3.65×** predicted vs 4.66× measured" | **Pre-fix.** `results/glioma_regional_pk.md` §P2 — post-fix the model predicts **17.20×**, an over-separation of 1.8×. `results/SESSION_AUDIT.md` §1 still shows 3.65× in a table; that table is superseded. | Quote the aggregate instead: ~2× geometric-mean fold error, 6/6 drugs within 10×. Or quote the regional ratios: 1.19× / 5.87× predicted against 1.45× / 4.26× measured. |
| "cisplatin and doxorubicin differ by about **4%** in predicted kill" | Predates both the discrete-vessel promotion and the delivery-law correction. | P(cisplatin rank 1) = 0.931 on the shipped one-way law; **0.397, with doxorubicin at 0.559, under the corrected gradient-driven source.** The order reverses. `results/ranking_robustness.md` §0. |

## Standing rules these drafts follow

- Lead with what has been **built and validated**, not with what we want.
- State the stage plainly: **research prototype, not clinically validated.** No patient has
  been or will be treated on its output.
- State what we offer in return: the validation harness, the reconstruction pipeline, a
  **pre-registered prediction**, and co-authorship on their terms — including none.
- Commit in writing to publishing a negative result. Say it in the first email so it is not a
  surprise in month six.
- Ask for the **smallest useful thing** first. The large ask can follow a reply.
- No fabricated credentials, affiliation, funding, prior contact, or endorsement.

## Contents

| File | To | Ask | Priority |
|---|---|---|---|
| `01_agar_santagata.md` | Nathalie Agar / Sandro Santagata, BWH-DFCI | The adavosertib GBM MALDI-MSI ("available on request" in both the 2022 method paper and the **2025 trial report**), starting with the non-patient PDX arm and the PXD018782 TMT design map | 3 |
| `02_erasmus_lamfers.md` | Martine Lamfers, GLIOscreen, Erasmus MC | Cisplatin + doxorubicin dose-response on a GSC culture whose parent specimen is also spatially profiled — the decisive experiment | **1** |
| `03_controlled_access_skeleton.md` | dbGaP / HTAN / PDMR / TCIA | Generic application skeleton with the specific data elements needed. **Not filable — blocked on institutional affiliation (gate G4).** | 4 |
| `04_holzel_gse324268.md` | The GSE324268 depositors (Hölzel lab, Bonn) | MGMT promoter-methylation status for `GBMtumor07`, the custom-panel design rationale, and a QC observation about their MGMT probes offered in return. Smallest, highest-probability ask in the set. | 2 |

**Not yet drafted, and it should be:** the Jonas Lab / Kibur Medical implantable-microdevice
glioma study (Peruzzi & Jonas, BWH; PMID 37672566; `results/zenodo_kibur_deposit_check.md`).
It pairs per-reservoir drug concentration with per-site caspase-3 / pH2AX response in human
glioma — and publishes it only as plotted points, with the underlying table "available on
request". It is the second glioma-specific gated asset, at the same institution as target 1.

Rationale, timelines, costs and the full access mechanics are in
[`../DATA_ACCESS_PLAN.md`](../DATA_ACCESS_PLAN.md).
