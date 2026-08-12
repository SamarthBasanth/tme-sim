# Pre-registration — per-patient cisplatin/doxorubicin IC50 on a spatially profiled GBM specimen

**Status: FILED, UNSENT.** This is the pre-registration that
[`docs/outreach/02_erasmus_lamfers.md`](../docs/outreach/02_erasmus_lamfers.md)
promises ("a timestamped, pre-registered prediction filed before I see anything
of yours"). That email must not be sent until this file is committed, because
the offer is the substance of the ask. Filed 2026-08-12; git history is the
timestamp.

**Not yet reviewed by a human.** Drafted from the repo's own published numbers.
The predictions below are falsifiable as written, which is the point; a human
should confirm they are willing to be held to them before the email goes out.

---

## 1. What is being registered, and why it is worth registering

The proposed experiment pairs one GBM resection's spatial section with a
per-patient cisplatin and doxorubicin IC50 pair from the derived culture.

The reason to pre-register rather than compare afterwards is specific: **this
project's two candidate delivery laws predict OPPOSITE drug orderings on the
same panel.** That makes the measurement a discriminating test between them
rather than a goodness-of-fit exercise — but only if the prediction is fixed
first, because with two arms available after the fact, one of them always
"agrees".

| | shipped one-way law `S = PS·Cp` | corrected exchange law `S = PS·(Cp − C)` |
|---|---|---|
| P(cisplatin rank 1) | **0.9277** | **0.3965** |
| P(doxorubicin rank 1) | **0.0430** | **0.5586** |
| kill, cisplatin (dx = 12 µm) | — | **0.2311** |
| kill, doxorubicin (dx = 12 µm) | — | **0.2547** |
| modal top-1 | cisplatin | doxorubicin |

Source: [`transvascular_exchange_revalidation.md`](transvascular_exchange_revalidation.md)
§4.1, from 10,240 complete re-rankings. `ranking_robustness.md` §0 publishes the
one-way pair at 0.931 / 0.041, a separate run of the same panel; the ordering
claim registered below does not depend on which of the two is used, and no point
value is registered.

These two drugs were not chosen for convenience. They sit at opposite ends of
both axes the ranking is decided on — roughly **157× apart in delivered exposure
and 49× in IC50** — which is why the ordering is sensitive to the delivery law
at all.

## 2. Registered prediction 1 — drug ordering

**Registered now, before any of the collaborator's data exists.**

Measured IC50 is an **input**, not a comparator: it replaces the population
prior ("the median glioma cell line" standing in for "this tumour"), which the
project's own sensitivity analysis names as the dominant remaining term.
Doxorubicin's prior is the widest on the panel for the mundane reason that it is
absent from the GDSC glioma screen. So the registered object is a **decision
rule**, not a number:

> **P1.** Given the measured pair (IC50_cisplatin, IC50_doxorubicin) from this
> patient's culture, the model is run on this patient's section under BOTH
> delivery laws, with no other parameter changed, and each arm emits a kill
> ordering. **The one-way arm is predicted to rank cisplatin first; the exchange
> arm is predicted to rank doxorubicin first.** If the measured response
> ordering is available, it agrees with at most one of them.

**Before unblinding, and committed to this file as an amendment:** the
**discriminating window** — the region of (IC50_cisplatin, IC50_doxorubicin)
space in which the two arms disagree. It is computable from this patient's
section and the panel's fixed delivery ratio without any response data, and it
must be committed *before* the IC50s are read. If the measured pair falls
outside that window, both arms agree and **the experiment is uninformative about
the delivery law** — that outcome is registered here in advance so it cannot
later be reported as a pass.

**Falsification.**
- Measured ordering matches the exchange arm and not the one-way arm → evidence
  for the corrected law, and a direct argument for flipping the default.
- Measured ordering matches the one-way arm → the corrected law makes a wrong
  prediction on real response despite being the physically correct form, which
  is the more interesting outcome and must be published as such.
- Measured pair outside the discriminating window → uninformative, reported as
  uninformative.
- Neither arm matches → both wrong; the ranking claim does not survive at all.

## 3. Registered predictions 2 and 3 — two-stage, and why

Per-region exposure and the MGMT effect size cannot be given numeric values now,
because both depend on the collaborator's tissue, which does not exist in this
project's hands. Registering a number for them today would be fabrication.

They are registered as a **two-stage procedure**, which is the standard remedy
and is binding:

> **P2 (per-region exposure).** On receipt of the spatial section and before any
> response or IC50 data is unblinded, the reconstruction and transport model are
> run and the predicted per-region exposure committed to this file as a dated
> amendment, with the run id.
>
> **P3 (MGMT effect size).** Same protocol. Note the standing limitation this
> project has already measured: MGMT is absent from every Xenium panel
> ([`DATA_HUNT.md`](DATA_HUNT.md)), and in the one deposit where regional MGMT
> was checked the probes failed QC against the deposit's own negative controls
> ([`mgmt_spatial.md`](mgmt_spatial.md)). If the collaborator's assay carries
> MGMT with adequate power, P3 is registerable; if it does not, **P3 is
> withdrawn rather than weakened**, and the email should not promise it.

**Blinding.** Stage-2 amendments must be committed before the response data is
looked at. The git history is the evidence; a commit after unblinding is not a
prediction and must not be presented as one.

## 4. Conditions the predictions carry

- **Every exchange-arm number is conditional on `PS_SCALE_S_INV = 1.0`**
  (`tmesim/treatment/drugs.py`), the shipped default and an unjustified one; the
  calibrated value is 0.16122 1/s. A calibrated scale changes the exchange arm's
  magnitudes. The registered claim is the **ordering**, which is what the
  calibration is not expected to move — and that expectation is itself
  falsifiable and should be checked in stage 2.
- The exchange arm's kill values are quoted at dx = 12 µm. The ordering persists
  at dx = 8 µm (cisplatin 0.2335 vs doxorubicin 0.2547), so it is not a mesh
  artefact.
- The comparison is **kill ordering**, not delivery/AUC ordering. These are not
  interchangeable: on breast, AUC ranking is perfectly preserved (ρ = 1.000)
  while kill ranking scrambles (ρ = 0.830), because the Hill map saturates.

## 5. Moves ruled out in advance

Registered here so that agreement obtained by any of them is not reportable:

1. **Refitting `alpha` to this data.** It is fit once, on GBM doxorubicin
   penetration depth, and passes an out-of-sample mitoxantrone test. Refitting
   would absorb the error of two other parameters and break that test. This
   project has declined to do it four times; declining a fifth is registered.
2. **Per-drug adjustment of `bbb_penetrance`** to recover an ordering.
3. **Changing `sensitivity_mode`** (uniform vs determinant) after seeing the
   data.
4. **Switching the compared quantity** from kill ordering to delivery/AUC
   ordering after the fact — see §4.
5. **Selecting the arm after the fact.** Both arms are reported, always, with
   this file cited.

## 6. Publication commitment

The email states: "I will publish the result if the model fails." That is
registered here too. A failure of P1 in either direction is a result about the
delivery law and gets written up in `results/` with the same treatment every
other falsified prediction in this repo has received — struck through and kept,
not deleted.

## 7. Cross-references

- The ask: [`docs/outreach/02_erasmus_lamfers.md`](../docs/outreach/02_erasmus_lamfers.md)
- The order reversal: [`transvascular_exchange_revalidation.md`](transvascular_exchange_revalidation.md) §4.1
- The retracted decision-grade claim: [`ranking_robustness.md`](ranking_robustness.md) §0
- Why the response axis needs an external collaborator at all:
  [`SESSION_AUDIT.md`](SESSION_AUDIT.md) §5, [`DATA_HUNT.md`](DATA_HUNT.md)
