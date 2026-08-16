# Outreach drafts — NOT SENT, NOT APPROVED

> ## ⚠️ EVERY FILE IN THIS DIRECTORY IS AN UNSENT DRAFT FOR HUMAN REVIEW.
>
> Nothing here has been sent. No contact has been made with any person or group named
> in these files. Nobody named here has agreed to anything, been approached, or knows
> this project exists.
>
> **A human must read, edit, and personally send each message.** Do not automate this.

## Before sending anything

1. **Fill every `[SQUARE BRACKET]` placeholder.** They mark facts only the sender knows —
   name, affiliation (or the honest absence of one), links, availability. **Do not invent an
   affiliation, a funder, a title, or a prior introduction.** "Independent researcher, no
   institutional affiliation" is a legitimate and sayable thing; a fabricated one ends the
   conversation permanently and correctly.
2. **Verify the contact details yourself the day you send.** Every address in these drafts
   came from a public source (a paper's corresponding-author line, an institutional faculty
   page) and is recorded with its source in `docs/DATA_ACCESS_PLAN.md`. Addresses rot, people
   move, and one of the lab websites referenced is currently down. Check.
3. **Re-read the honest-position paragraph in each draft and confirm it is still true.** It
   states: 11 defects found and fixed; a structure-only permeability coefficient **scored** in
   human glioma to ~2× geometric-mean error (**not** "delivery validated" — those rows never
   invoke the transport solver, and the three checks that do compare a solved field to an
   external measurement all fail); two components falsified; response prediction unvalidatable
   with open data; research prototype, not clinically validated. If any of that has changed,
   update the draft — do not let an outdated claim go out over your name.

   > **This step has already failed once, on 2026-08-16.** `02_erasmus_lamfers.md` was still
   > asserting *"the reconstruction agrees across two independent spatial platforms at
   > Pearson 0.997"* after that number had been retracted in `results/revalidation.md` — the
   > non-tumour subvector had collapsed **0.945 → 0.016**. The draft is corrected and carries
   > a correction banner. **The word to police in these drafts is "validated."** This project's
   > documented failure mode is that the engine stays sound while the claim surfaces drift, and
   > an outreach letter is the surface where drift reaches a named external person.
4. **Decide the commercial question before you send.** Erasmus MC's biobank regulation routes
   any commercial involvement through their TTO rather than the standard non-commercial MTA
   (§8.3.3). If this project has or may have a commercial dimension, say so in the first email.
   Discovering it later is the fastest way to lose a collaboration.
5. **Shorten.** These drafts are longer than a first cold email should be. A busy PI reads the
   first three lines. Cut anything that is not (a) what we built, (b) what we already did with
   *their* data, (c) the one specific thing we are asking for.

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
| `01_agar_santagata.md` | Nathalie Agar / Sandro Santagata, BWH | The adavosertib GBM MALDI-MSI ("available upon request"), starting with the non-patient PDX arm and the PXD018782 TMT design map | 3 |
| `02_erasmus_lamfers.md` | Martine Lamfers, GLIOscreen, Erasmus MC | Cisplatin + doxorubicin dose-response on a GSC culture whose parent specimen is also spatially profiled — the decisive experiment | **1** |
| `03_controlled_access_skeleton.md` | dbGaP / HTAN / PDMR / TCIA | Generic application skeleton with the specific data elements needed. **Not sendable — blocked on institutional affiliation.** | 4 |
| `04_holzel_gse324268.md` | The GSE324268 depositors (Hölzel lab, Bonn) | MGMT promoter-methylation status for `GBM_tumor07`, plus the custom-panel design rationale. Smallest, highest-probability ask in the set. | 2 |

Rationale, timelines, costs and the full access mechanics are in
[`../DATA_ACCESS_PLAN.md`](../DATA_ACCESS_PLAN.md).
