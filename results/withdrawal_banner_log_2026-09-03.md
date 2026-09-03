# Withdrawal-banner log — 2026-09-03

**What this pass did, and the only thing it did.** Ten files in `results/` carried
a cross-compound ordering as a live-looking claim with no notice anywhere in the
file that the claim had been withdrawn on **2026-08-21**. Each has had a dated
banner **prepended immediately after its H1**. No existing body text was altered,
reworded, re-scored or deleted anywhere. Verified per file: `git diff --numstat`
shows `+34 / -0`, a single hunk `@@ -2,0 +3,34 @@`, and the H1 plus every byte
from the old line 3 onward hashes identically to the committed blob at `7742351`.

The banner pattern is this project's own, established by
[`results/numerical_rigor.md`](numerical_rigor.md)'s correction block and by
[`results/CLAIM_INVENTORY.md`](CLAIM_INVENTORY.md)'s 2026-08-14 block. The
wording model is `README.md`'s withdrawal box.

**Source of the finding:**
[`results/claim_surface_audit_2026-08-30.md`](claim_surface_audit_2026-08-30.md),
which named three of these ten. The other seven came from an independent sweep of
`results/*.md` described in §3 below.

---

## 1. The reasons cited in every banner

Both are measurements in this repo, not arguments, and both are quoted from the
files that made them:

* [`results/pk_ranking_robustness.md`](pk_ranking_robustness.md) — per-compound
  plasma exposure across the retained compounds spans **113.1x**; the published
  order first breaks at a spread of **1.53x**.
* [`results/capacity_ranking_robustness.md`](capacity_ranking_robustness.md) —
  there is **no per-compound tissue-capacity term** in the engine at all, and
  under a *uniform* capacity the published order still breaks, at **R = 36.3**
  (GBM 10-drug) and **R = 5.28** (breast 10-drug), both below the best available
  constant of **39**.

No banner restates the withdrawn sentence, and no banner reprints a compound
ordering.

## 2. The three files named by the audit

| file | what was found | verified how |
|---|---|---|
| [`ranking_robustness.md`](ranking_robustness.md) | whole file; its subject is the withdrawn claim | `grep -c` for `withdraw` / `CORRECTION` / `2026-08-21` each returned **0** |
| [`alpha_identifiability.md`](alpha_identifiability.md) | line 19 names a cross-compound ordering as the one surviving output | `grep -c -i withdraw` returned **0**; line 19 read directly |
| [`CLAIM_INVENTORY.md`](CLAIM_INVENTORY.md) | §2 item 2, at lines **337–338**, lists the withdrawn claim under *"what a competent reviewer would accept today"* | lines read directly; rows **C2** (line 91) and **C3** (line 128) *were* annotated on 2026-08-21, so the miss is confirmed as a miss, not an omission of the whole file |

## 3. The seven found by sweep

`results/*.md` was swept two ways — for an explicit compound-ordering chain
(`X > Y` over the shipping panels' compound names) and for files with ten or more
mentions of *rank* / *ordering* — then every hit was read to decide whether the
withdrawn ordering is **load-bearing** in it, and each was checked for an existing
notice with `grep -c -i withdraw`.

| file | why it is in the same state | notice before this pass |
|---|---|---|
| [`real_drug_panel.md`](real_drug_panel.md) | its headline deliverable is a ten-compound kill-ordered table, plus §*Sensitivity mode*, which exists to argue that ordering is fair | none (its 3 hits are `retracted`, about the dx = 30 µm penetration figures) |
| [`patient_report_xenium.md`](patient_report_xenium.md) | §2 is a ranked virtual compound panel; the report opens by promising one ranked panel per section | none |
| [`patient_report_visium.md`](patient_report_visium.md) | same, spot-level, at dx = 100 µm | none |
| [`flagship_markers.md`](flagship_markers.md) | headline reports the cross-compound ordering and names the ordering as the quantity that transfers between sections | none |
| [`heldout_patient.md`](heldout_patient.md) | §3 asks whether the ordering holds on a second patient; §6 and §6.5 report ordering results as generalising | none |
| [`DRUG_SELECTION_RIGOR.md`](DRUG_SELECTION_RIGOR.md) | its entire question is per-patient compound selection, scored as top-1 accuracy and regret | none |
| [`bring_your_own_compound.md`](bring_your_own_compound.md) | §6.1's **CAN** clause grants permission to compare a supplied compound against others and read the ordering | none for the ranking; its 2026-08-14 banner is about the diffusivity law |

`CLAIM_INVENTORY.md` and `bring_your_own_compound.md` each now carry **two**
banners. Each new banner says so, and says the older one is about something else.

## 4. Considered and deliberately left alone

None of these was touched. Each is recorded with the reason, so the decision can
be reversed by someone who disagrees with where the line was drawn.

**(a) The withdrawal's own evidence, or documents that already record it.**
[`pk_ranking_robustness.md`](pk_ranking_robustness.md),
[`capacity_ranking_robustness.md`](capacity_ranking_robustness.md),
[`uptake_proxy_double_count.md`](uptake_proxy_double_count.md) (written
2026-08-21 *after* the withdrawal; its subject is the ranking failure and it
reinforces it), [`axis_inventory.md`](axis_inventory.md) (line 18 already records
that an unswept axis broke the ranking claim on 2026-08-21), and
[`transvascular_exchange_revalidation.md`](transvascular_exchange_revalidation.md)
(§4.1 is titled for the ranking-invariance claim being false). Bannering these
would misdescribe them.

**(b) Files where an ordering appears as an invariance probe, not as the
document's output.** [`cell_state_impact.md`](cell_state_impact.md),
[`oxygen_level_shift.md`](oxygen_level_shift.md),
[`immune_sink_impact.md`](immune_sink_impact.md). Their subjects are the
cell-state route, the oxygen route and the TAM sink. `README.md`'s C2 row already
states that the subject of those invariances no longer exists and that a true
invariance is not a licence to order compounds.

**(c) External-anchor comparisons against a measured ordering.**
[`kyle2004_penetration.md`](kyle2004_penetration.md) (already carries a
2026-08-22 correction), [`moncayo_explant.md`](moncayo_explant.md),
[`pk_multidrug.md`](pk_multidrug.md) Test 2. These score a model quantity against
a published measurement; they do not offer an ordering as a product output.

**(d) Same-compound orderings, which the withdrawal does not reach.**
[`regimen.md`](regimen.md) (temozolomide schedules ranked by cumulative dose),
[`uncertainty.md`](uncertainty.md) (a *parameter* sensitivity ranking).

**(e) The strongest borderline case, flagged rather than bannered.**
[`kill_calibration.md`](kill_calibration.md) has **no** withdrawal notice, and two
headline blockquotes (around lines 34 and 46) assert the cross-compound ranking as
invariant and as the simulator's primary output. That is a live assertion of the
withdrawn claim. It was left alone because the *subject* of the file is the
absolute-kill 3-D:2-D ratio test, not the ranking, and this pass was scoped to
documents whose subject was abolished. **It is the first place to look if this
pass is extended.**

**(f) Named as out of scope by this pass's own brief.**
[`delivery_metrics.md`](delivery_metrics.md) — the audit flags it for a stale
`alpha` (it prints 148.6; live `calibrate_alpha()` is `22.223477711172563`) and
for a `README.md` line that routed compound-comparison questions to it. That
`README.md` line has moved since the audit was written, and the file's subject is
the delivery metrics rather than the ranking, so it is out of scope here. The
stale `alpha` is a separate defect and is still open.
[`heldout_patient_validation.md`](heldout_patient_validation.md) — discusses one
ordering change, but its subject is what the held-out anchor does and does not
test.

## 5. `p90_capacity_sweep.md` — reported, not touched

Deliberately left exactly as found, as instructed. What is there:

* [`results/p90_capacity_sweep.md`](p90_capacity_sweep.md) line 3 reads
  **"Status: PRE-DECLARED, NOT YET RUN."** (file mtime 2026-08-23 03:34).
* `results/p90_capacity_sweep.json` (mtime **2026-08-24 23:33**, i.e. *later*)
  contains a completed Arm A run: `arm_A.analysis.C1.fail = True`,
  `arm_A.analysis.C2.fail = True`, and
  `arm_A.analysis.verdict = "THE HEADLINE IS CAPACITY-SENSITIVE AND SO IS THE
  CLAIM BUILT ON IT"` — the fail/fail row of the verdict table the document itself
  pre-declared at §0.3.

So the sweep was run, it landed on the worst of the four pre-declared outcomes,
and the document still tells a reader it has not been run. **That is a live
decision for the repository owner** — the failing verdict concerns
`E_p90/E_p10`, the shipped headline delivery number, and publishing it is a
substantive call, not a cleanup. No banner was added and no byte of either file
was changed.
