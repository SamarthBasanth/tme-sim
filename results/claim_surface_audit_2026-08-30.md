# Claim-surface audit — 2026-08-30

**Status: RUN AND MEASURED.** Read-only. Nothing under `tmesim/`, `app/`, `docs/`,
`landing-3d/`, `experiments/`, `scripts/` or any pre-existing file in `results/`
or `archive/` was created, edited or deleted by this audit. This file is the only
file written.

**Scope.** The *claim surface* — text a reader **outside** this repository can
reach: `README.md`, `results/*.md`, `docs/**`, `landing-3d/**` copy, and
user-visible strings in `app/server/report.py` and `app/web/js/**`. Auditing
against, not around, the changes that landed **today**: `app/server/account.py`
(`GET /api/account`, `GET /api/limits`), the `pruned_payloads` record,
`app/web/js/ui/account.js`, and `scripts/prune_run_payloads.py`.

**The defect being hunted** is this project's named recurring one:
*claim-surface drift* — the engine is sound, the DOCUMENTS are wrong, and
`results/` is the claim surface and has no CI.

---

## 0. Verdict up front

**The engine and the Studio UI are, on nearly every axis checked, the most
current and best-caveated surfaces in the repository. The drift is in the
documents that sit outside them — and in one `results/` file whose own data
contradicts its own status line.**

| category | count | of which `docs/` | what it means |
|---|---:|---:|---|
| **(d) PROHIBITED outright** | **39** | 29 | text asserting or licensing something the project has banned |
| **(c) CONTRADICTED by current code / filesystem** | **86** | 75 | measurably false against the tree as it stands today |
| **(b) SUPERSEDED by a later dated result** | **62** | 44 | true when written, false now because a later measurement landed |
| **(a) STILL TRUE (notable passes)** | **30+** | — | recorded because a pass that is not recorded gets re-litigated |
| **borderline — recorded, not counted** | **3** | — | judgment calls, argued in place |
| **(e) UNVERIFIED** | **11 classes** | — | I could not check these; they are neither pass nor fail |

**Total violations (b + c + d): 187 distinct claim instances, 148 of them in
`docs/`.**

Counts are **claim instances**, not table rows — a row covering nine loose uses of
*"validated"* counts nine. Two honesty notes on this table:

* **An earlier draft of this header said 8 / 13 / 17.** It was written before the
  `app/web/**` and `docs/**` sweeps landed and was under-counted on every line.
  The correction is left visible rather than silently applied, since a tidied
  count is the exact species of defect this audit exists to find.
* **`docs/` was audited by delegation and most of its rows are not personally
  re-derived.** Rows in §2.3 marked **[D]** come from that sweep. I spot-checked
  the five most serious and all five held — **and one delegated reading was wrong
  in my favour**: `docs/DEPLOYMENT.md:743` was relayed to me as saying *"Six."*
  and actually says *"Ten."* (my own count of `TMESIM_*` variables: **12**). The
  verdict survives the correction; the quote did not. Treat unspot-checked [D]
  rows as high-confidence but second-hand, and the `docs/` sub-totals as
  approximate to within a few instances.

The single most consequential finding is not in any of the files the task named.
It is `results/p90_capacity_sweep.json`, whose recorded verdict is
`"THE HEADLINE IS CAPACITY-SENSITIVE AND SO IS THE CLAIM BUILT ON IT"` — while
`results/p90_capacity_sweep.md`, the published record beside it, still reads
**"Status: PRE-DECLARED, NOT YET RUN."**

---

## 1. Methodology

### 1.1 How the surface was enumerated

Mechanically, by path, not by memory:

```
$ ls results/*.md | wc -l
     144
$ find docs -type f | wc -l
      25
```

plus `README.md`, the `landing-3d/` copy tree (`lib/copy.ts`, `lib/timeline.ts`,
`app/page.tsx`, `app/layout.tsx`, `components/ui/*.tsx`, `components/canvas/*.tsx`,
`components/scroll/*.tsx`, `landing-3d/README.md`, excluding `node_modules/` and
`.next/`), `app/server/report.py`, and all 59 files under `app/web/`
(`js/**`, `index.html`).

Two files that are formally *inside* `results/` but are 10.9 MB frozen HTML
snapshots of the whole Studio — `results/studio_export.html` and
`results/studio_export_artifact.html` — were treated as claim surface and audited
by content search rather than by reading.

### 1.2 An environment trap that fired, and nearly produced a false negative

`grep` on this machine is **`ugrep`**, and it does **not** word-split an unquoted
shell variable holding a file list. The first prohibited-term sweep of this audit
was written as `grep -rin "$pat" $SURF ... 2>/dev/null` and returned **zero hits
for every pattern including `validated`, `digital twin` and `upload`** — which
would have been reported as a clean surface. The `2>/dev/null` hid the cause:

```
$ SURF="README.md docs app/server/report.py app/web/js results/CLAIM_INVENTORY.md"
$ grep -rin "validated" $SURF | head -3
ugrep: warning: README.md docs app/server/report.py app/web/js results/CLAIM_INVENTORY.md: No such file or directory
```

Every sweep in this audit was re-run with **explicit literal paths**. This is
recorded because it is the same failure mode the audit is about: a confident
result pointed at the wrong object.

`python` does not exist in this environment; every Python invocation below is
`.venv/bin/python`. No command in this audit was judged by grepping its output
for a success word.

### 1.3 Verdict vocabulary

* **(a) STILL TRUE** — checked against code, the filesystem, or a later dated
  `results/` file, and it holds.
* **(b) SUPERSEDED** — true when written; a later dated result measured
  otherwise. The claim's *date* is the defect, not its arithmetic.
* **(c) CONTRADICTED** — measurably false against the tree as it stands today.
* **(d) PROHIBITED** — asserts or licenses something on the banned list
  regardless of whether it is true.
* **(e) UNVERIFIED** — I could not check it. **Recorded as unverified, never as
  either a pass or a fail.** §5 enumerates these.

Text that *records a withdrawal* or *narrates a historical failure* is **not** a
violation and is not counted as one. That distinction is applied throughout;
where it is close, the row says so.

### 1.4 Ground truth established first, from code

Every verdict below rests on these, which were measured before any document was
read:

| fact | evidence |
|---|---|
| The product ships the **STEADY-STATE** reach map | `app/server/jobs.py:1064` calls `compute_reach(..., drug=dp, convection=PDE_CONVECTION)`; `app/server/reach.py:747` — *"``drug`` GIVEN -> ``window_mode = "steady_state"``"* |
| In the steady problem a capacity term `R` is **ABSENT**, not small | `app/server/reach.py:391-396` — *"it multiplies the TIME DERIVATIVE AND NOTHING ELSE, and the steady problem has no time derivative, so R is absent from it by inspection — not small, not bounded, ABSENT"* |
| The steady solve is **48.6× cheaper to solve** but makes a run **4.2 % SLOWER** | `app/server/reach.py:36-43` — *"0.90 s against 43.70 s… **This is not a performance win for the product and must not be quoted as one.**… adds 1.854 s, **+4.2 %**"* |
| Cross-compound ranking is **withdrawn absolutely** | `app/server/reach.py:5-7`; `pyproject.toml:90`; `app/server/api.py:2115-2118` (*"the response carries no ``ranking`` key"*); `tmesim/protocol/executor.py:1019-1021` |
| Live `calibrate_alpha()` = **22.223477711172563** | `$ .venv/bin/python -c "from tmesim.validation.calibration import calibrate_alpha; print(repr(calibrate_alpha()))"` → `22.223477711172563` |
| **54** runs pruned today, **19** still hold `timelapse.bin`, **96** runs total | `$ grep -l pruned_payloads app/runs/*/result.json \| wc -l` → 54; `$ ls app/runs/*/timelapse.bin \| wc -l` → 19; `$ ls app/runs/r_*/result.json \| wc -l` → 96 |

---

## 2. The full table

### 2.1 `README.md`

| file:line | claim | verdict | evidence |
|---|---|---|---|
| `README.md:509-511` | *"Vivum is a research instrument whose defensible output is the **relative comparison of compounds on one tissue section**"* | **(d) PROHIBITED** | This is, word for word, the sentence `README.md:28-29` says was **withdrawn on 2026-08-21**: *"This sentence used to read 'the relative comparison of COMPOUNDS on one tissue section', and that claim was withdrawn on 2026-08-21."* The retraction is at the top of the page; the retracted sentence is still live 480 lines below it, in the section headed **"The honest summary"**. `README.md:34-36`: *"The prohibition is absolute rather than conditional."* |
| `README.md:513-515` | *"the diffusivity law reproduces about a quarter of the measured compound-to-compound spread, **so a ranking over similar-mass analogues is weak evidence rather than a measurement**"* | **(d) PROHIBITED** | Grades a cross-compound ranking as *weak evidence*, i.e. as evidence. The prohibition admits no grade: `landing-3d/lib/copy.ts:146-149` — *"No ordering, ranked list, shortlist, 'top pick', 'which compound wins', pairwise call or cross-compound comparison may appear… in any form, hedged or not, positive or negative — including as a limit."* |
| `README.md:542-543` | *"the only item on this list that changes what a ***ranking*** means"* | **(d) PROHIBITED** | Frames the project's top improvement priority around a product output that does not exist. `pyproject.toml:90`: *"`compare-drugs` was withdrawn 2026-08-21 with the cross-drug ranking claim."* |
| `README.md:645` | *"\| How should I compare drugs? \| [`results/delivery_metrics.md`](results/delivery_metrics.md) \|"* | **(d) PROHIBITED** | A navigation table row offering cross-compound comparison as a supported question, pointing at a document (§2.2) that is itself superseded on two axes. |
| `README.md:56-58` | *"it reaches **R² = 0.718, Spearman +0.821**, with **99 %** of compounds within 3× of measurement, a geometric bias of **1.03×**, and **24 of 31** exact-mass pairs ordered correctly"* | **(b) SUPERSEDED** | `results/hbd_convention_audit.md:3-7` (2026-08-21): *"the published `R² = 0.718` describes neither the shipping engine nor the correct physics."* `:267-270`: *"the residual 26 rows leave R² anywhere in **0.498–0.839**, an interval 0.34 wide… `0.718` is **neither refuted nor confirmed**."* README states 0.718 flat and cites neither `hbd_convention_audit.md` nor `descriptor_falsification.md`. |
| `README.md:438` | *"8.3× spread, **R² = 0.718** as shipped"* | **(b) SUPERSEDED** | Same. |
| `README.md:539` | *"the two-route promotion (R² 0.022 → 0.718, blind)"* | **(b) SUPERSEDED** | Same. |
| `README.md:639` | *"`results/CLAIM_INVENTORY.md` — **the single settled answer**, and the source for most of this page"* | **(b) SUPERSEDED** | `CLAIM_INVENTORY.md:41` — *"**Compiled 2026-08-12.**"* Its last commit is 2026-08-21 (`git log -- results/CLAIM_INVENTORY.md`). `grep -n 'steady\|reach map\|verdict_inversion\|08-22\|08-23\|08-25' results/CLAIM_INVENTORY.md` returns **one** hit, at `:425`, unrelated. It knows nothing of the 2026-08-22 steady-state cutover, the Aug-22 accuracy/inversion work, the Aug-23 p90 pre-registration or the Aug-25 verdict-inversion census. |
| `README.md` **as a whole** | the page presents itself as the project's current claim register | **(b) SUPERSEDED** | Date-stamping every `results/*.md` README cites, by that file's last commit: the **newest is 2026-08-21**. README cites **none** of `steady_state_reach.md`, `reach_invariances_measure.md`, `permeability_split.md`, `verdict_inversion_census.md`, `p90_capacity_sweep.md`, `accuracy_gap_decomposition.md`, `kyle2004_penetration.md`, `spatial_noise_decomposition.md`, `mechanism_attribution.md`, `delivery_spread_noise.md`, `hbd_convention_audit.md`, `dx_sensitivity_steady.md`. Nine days of published measurement, including the cutover that changed what object the product computes, are absent from the outermost surface. |
| `README.md` — steady-state cutover | absent entirely | **(b) SUPERSEDED** | No occurrence of `window_mode`, `steady_state` or the reach map's mode anywhere in README. The product's headline map changed object on 2026-08-22 (`app/server/reach.py:21-23`: *"WHICH OBJECT: the map is a STEADY-STATE object (2026-08-22)"*) and README does not mention it. |
| `README.md:191` | *"**Studio** — the browser UI (reconstruct, **compare**, protocol editor, caveats panel)"* | **(a) STILL TRUE — with a note** | The surviving `compare` is same-compound cross-specimen, not cross-compound: `app/server/api.py:2107-2112` — *"Compare the matched runs of a project — the SAME compound across specimens… that ranking was withdrawn… and the response carries no `ranking` key."* Defensible, but the bare word "compare" beside a withdrawn ranking is a foreseeable misread. **Not counted as a violation.** |
| `README.md:216` | *"The sample directory is read where it lies and **never uploaded anywhere**"* | **(a) STILL TRUE** | Uses the word "upload" about tissue, but to deny it, and it states the path-based model correctly. `docs/DEPLOYMENT.md:74-76`: *"Sections are registered by path, not by upload… There is no HTTP upload route for tissue."* **Not a violation.** |
| `README.md:87-95` | the *"What Vivum does not do"* block | **(a) STILL TRUE** | Every prohibited category is negated explicitly and correctly. |
| `README.md:1-11, 20-48` | the top-of-page withdrawal block | **(a) STILL TRUE** | Accurate, dated, and cites its two sources. It is the part of README that has kept up. |

### 2.2 `results/*.md`

| file:line | claim | verdict | evidence |
|---|---|---|---|
| `results/p90_capacity_sweep.md:3` | *"**Status: PRE-DECLARED, NOT YET RUN.**"* — and §2 is empty below *"Everything below this line is written after the sweep."* | **(b) SUPERSEDED — the most serious finding in this audit** | The sweep **was run**. `results/p90_capacity_sweep.json` is dated **2026-08-24 23:33**, the `.md` **2026-08-23 03:34**. Read with `.venv/bin/python`: `arm_A.analysis.verdict` = **`"THE HEADLINE IS CAPACITY-SENSITIVE AND SO IS THE CLAIM BUILT ON IT"`**; `C1.fail` = `true` (temozolomide moves **1.1862** from R=1 to R=39 against its own nuisance width 0.202 — **5.87×** that width; doxorubicin **1.0155** against 0.226 — **4.49×**); `C2.fail` = `true` (`min_separation` **1.5039** at R=39, against the pre-declared floor **1.74**); `elapsed_s` 1232.7; `alpha` asserted at `22.223477711172563`; the R=1 reproduction check passes (`2.1107` vs published `2.111`; `3.9425` vs `3.942`). The pre-registration's own verdict rule (`.md` §0.3) maps fail/fail to *"the headline is capacity-sensitive and so is the claim built on it"*. **The subject is `E_p90/E_p10`, which is a live shipped headline** — `app/web/js/ui/reach.js:638-646` renders `num(p90 / p10, 2) + 'x'` labelled **`DELIVERY SPREAD`**. |
| `results/CLAIM_INVENTORY.md:337-338` | *"2. **Relative comparison of compounds on one section**, presented as a ranking with a stated tie at the top and stated invariance limits."* — listed under **"## 2. What a competent reviewer would accept today → Would accept:"** | **(d) PROHIBITED** | An unhedged endorsement of cross-compound ranking as an acceptable product output, in the file `README.md:639` calls *"the single settled answer."* The same file's rows C2 (`:91`) and C3 (`:92`) were updated on 2026-08-21 to record the falsification; §2 line 337 was missed. `git log --format=%ad --date=short -- results/CLAIM_INVENTORY.md` → `2026-08-21, 2026-08-21, 2026-08-12`. |
| `results/kyle2004_penetration.md:740-742` | *"**The engine's absolute perivascular gradient is ~14× too shallow**… (daunorubicin 3.65 against 50 ± 15)."* — in the file's own **conclusion** | **(b) SUPERSEDED — by a correction 214 lines above it in the same file** | `:504-526`, headed **"CORRECTED 2026-08-22 — none of it is construction, and the '14×' is not a length error"**: *"the transport length `λ` is **2.45× too long**… The 14× is that error amplified by `dlnR/dlnλ ≈ 2.9`."* And *"**on the product's own readout the amplification does not transfer.** On the Xenium GBM section the map's half-distance responds to `λ` with exponent **0.92–1.22**… A 2.45× `λ` error is **~2.2–3.0× on the map**, not 14×."* The uncorrected 14× also survives at `:50`, `:89`, `:103`, `:496`, `:667`. A correction that the file's own summary does not adopt is not a correction the next reader will see. |
| `results/alpha_identifiability.md:19` | *"The cross-drug **ranking is the one thing that holds**"* | **(d) PROHIBITED / (b) SUPERSEDED** | No withdrawal banner anywhere in the file (`grep -c -iE 'WITHDRAWN\|no longer ranks' results/alpha_identifiability.md` → 0). Last commit 2026-08-08, i.e. before the 2026-08-21 withdrawal. |
| `results/ranking_robustness.md` — whole file, `:1` title *"Ranking robustness — is the per-patient drug ranking DECISION-GRADE?"* | presents cross-drug ranking as the product's output | **(d) PROHIBITED / (b) SUPERSEDED** | No withdrawal banner (`grep -c` → 0). Last commit 2026-07-28. It is cited approvingly from `README.md:299`, `README.md:401`, `docs/architecture.md:299` and `results/CLAIM_INVENTORY.md:91`,`:327`. A reader arriving from any of those four sees a document whose subject was abolished, with nothing on the page saying so. |
| `results/delivery_metrics.md:9` | *"`calibrate_alpha()` → **α = 148.6**, fit ONCE on doxorubicin"* | **(c) CONTRADICTED by current code** | Live `calibrate_alpha()` = `22.223477711172563` (§1.4). Aggravating: `README.md:645` sends every reader asking *"How should I compare drugs?"* to this file. No withdrawal banner (`grep -c` → 0). |
| `results/studio_export.html`, `results/studio_export_artifact.html` (10.9 MB each, 2026-08-22 08:59) | frozen full-Studio snapshots published as records | **(b) SUPERSEDED** | Both contain `DELIVERY SPREAD` (1 hit each) and the singular verdict-inversion wording — *"on a shipped compound at shipped parameters IT INVERTS"* — and **zero** occurrences of `verdict_inversion_census` or `p90_capacity_sweep`. `results/verdict_inversion_census.md:0` (2026-08-25) supersedes the singular framing: *"It reverses for **7 of the 12** scored compounds — 58.3 % of the shipped panel"*, temozolomide included, with **57.2 %** of its cells relabelled. The live UI has been updated (`app/web/js/ui/finding.js:345-347`); these two exports have not and cannot be, being dated records. |
| `results/dx_sensitivity_steady.md:163` | check **C** cites *"the published container `kill_fraction` `0.11221232991768859` **and grid `286,200`** (`docs/DEPLOYMENT.md`)"* | **(c) CONTRADICTED — half the citation** | `grep -n '0.1122123' docs/DEPLOYMENT.md` → `272:` ✓. `grep -n '286,200\|286200' docs/DEPLOYMENT.md` → **no match**. The grid figure lives at `app/server/api.py:163`, `results/oxygen_sensitivity_band.md:27` and `results/steady_state_reach.md:133`. Low severity; recorded because a pre-registration's guard is only as good as its anchor. |
| `results/dx_sensitivity_steady.md:3` | *"**Status: PRE-REGISTRATION ONLY — written before the sweep was run.**"* | **(a) STILL TRUE** | Created today 00:57. `results/dx_sensitivity_steady_data.json` does **not exist** (`ls` → *no matches found*); the driver `experiments/dx_sensitivity_steady.py` does (32,881 bytes, 00:56). This is an **honest open pre-registration** — the correct state — and is listed here only so it is not confused with the p90 row above, which is the same shape gone wrong. |
| `results/verdict_inversion_census.md` (2026-08-25) | *"the flagship compound, on the flagship section, **in the configuration the landing page and every results file are written around**, flips from dose to delivery with 57.2 % of its cells changing label"* | **(a) STILL TRUE — and it is an indictment of the rest of this table** | The file names its own blast radius. `grep -rn 'verdict_inversion_census'` across the claim surface returns hits **only** in `app/web/js/ui/finding.js:375` and `app/server/jobs.py:199,239`. Neither `README.md`, nor `docs/**`, nor `landing-3d/` cites it. |
| `results/steady_state_reach.md:26-27` | *"**0.90 s** against **43.70 s**… **48.6×**"* / *"cost added to a RUN: **+1.854 s**, i.e. **+4.2 %**. The run still marches; **see §2 before quoting the 48.6×**"* | **(a) STILL TRUE** | Matches `app/server/reach.py:36-43` exactly. The "60× cheaper, 4.4 s vs ~250 s" figures do not appear anywhere on the claim surface (swept: `grep -rn '60x\|60×\|4\.4 s\|~250\|250 s'` over README, docs, app, results/*.md, landing-3d/lib, API_CONTRACT.md — every hit is an unrelated number). **That retraction held.** |
| `results/reach_invariances_measure.md:616-621` | the *"~10× grid sensitivity"* retraction (measured 0.6 %) | **(a) STILL TRUE** | Propagated into `app/web/js/ui/reach.js:614-628` as a dated correction block, and into `results/p90_capacity_sweep.md:22-26` and `results/spatial_noise_decomposition.md:97`. The fabricated 10× claim is not asserted anywhere live. **That retraction held.** |
| `results/hbd_convention_audit.md` (2026-08-21) | *"the published `R² = 0.718` describes neither the shipping engine nor the correct physics"* | **(a) STILL TRUE — and orphaned** | `grep -rn 'hbd_convention_audit' README.md docs app/web/js app/server landing-3d/lib` → hits **only** at `landing-3d/lib/copy.ts:87` and `:416`. The landing page is the **only** surface carrying this correction. See §2.4 and §2.5. |

### 2.3 `docs/**` — all 25 files

**Coverage note.** The seven priority files plus `VERIFICATION_METHODOLOGY.md`
and the 17 remaining files were swept exhaustively. Rows marked **[D]** are
**delegated findings I did not personally re-derive**; rows without the marker I
verified myself. Five of the most serious [D] rows were spot-checked and are
reported here as verified — one of them with a correction to the delegated
reading, recorded in place.

**The structural finding: `docs/` does not know the product's current output.**
`grep -rn -iE 'reach map|compute_reach|window_mode|under.reach|E_ref'` over the
whole tree returns only oxygen-solver hits. The steady-state reach map — which
`app/web/js/ui/reach.js:1` calls *"the discrimination this product now sells"* —
is documented **nowhere in `docs/`**, nor is the capacity term `R`, nor the
+4.2 % run cost. Corollary, and it is a genuine pass: **the "60× cheaper" framing
is absent from `docs/`, not wrong.** Verified: `grep -rn -E '60x|48\.6|0\.90 s|43\.70|4\.4 s|250 s' docs/` returns only false positives from "148.6".

#### Prohibited

| file:line | claim | verdict | evidence |
|---|---|---|---|
| **`docs/DEPLOYMENT.md:15-17`** | *"**Vivum reconstructs a tissue section and ranks compounds against it**; it performs no diagnosis and produces nothing that is intended to select therapy for a patient."* | **(d) PROHIBITED** | **Verified verbatim.** The withdrawn claim stated as *what the product does*, in the opening scope paragraph of the deployment document — inside the very sentence whose job is to bound scope honestly. An operator reading only §0 to decide whether they may deploy this is told the product ranks compounds. |
| **`docs/VERIFICATION_METHODOLOGY.md:413-417`** | *"Ranked by `penetration_index` the panel comes out `cisplatin > lomustine > temozolomide > methotrexate > doxorubicin`. Ranked by the replacement metric it comes out `cisplatin > doxorubicin > lomustine > methotrexate > temozolomide` — which is the *identical* order the simulated kill fraction produces from an independent computation."* | **(d) PROHIBITED** | **Verified verbatim.** Two complete five-compound orderings printed in full, and one cross-validated against another. This is not a hedge or a fossil — it is the output itself. [D] adds that the ordering is also stale: `results/delivery_metrics.md:98-104` now gives TMZ 0.5107 > cisplatin 0.4283 > lomustine 0.3743 > MTX 0.3273 > dox 0.3044. |
| **`docs/architecture.md:700`** | *"**Read a ranking over similar-mass analogues as weak evidence.**"* | **(d) PROHIBITED** | **Verified verbatim.** The identical construction to `README.md:513-515` and `tmesim/molecular/byo.py:222-223` — the same hedge, on a third surface. `docs/QUICKSTART.md:266` says the opposite: *"no longer offered as a ranking by this tool."* |
| `docs/architecture.md:497-500` | *"Compare drugs by `mean_tumor_auc` … or `therapeutic_coverage` … **genuinely drug-discriminating**"* | **(d) PROHIBITED** [D] | Cross-compound comparison offered as the recommended workflow. |
| `docs/REPRODUCE.md:555` | recommends `therapeutic_coverage` *because it discriminates between compounds* | **(d) PROHIBITED** [D] | Same shape. |
| `docs/REPRODUCE.md:540` | *"**cross-drug ordering passes** at 3.65× predicted vs 4.66× measured"* | **(d) PROHIBITED + (b) SUPERSEDED** [D] | A cross-compound ordering reported as a PASS. Also stale: the shipping geometric-mean model predicts **17.20×** (`results/glioma_regional_pk.md:241`); 3.65× is the retired arithmetic-mean value. |
| `docs/DATA_ACCESS_PLAN.md:25, :27, :80-84, :92-99, :112, :626-630` | *"cross-drug penetrant-vs-excluded separation 3.65× predicted"*; *"the delivery×potency **kill ordering** … survives"*; pre-registers *"the predicted **drug ordering**"* and *"the predicted **cisplatin/doxorubicin ordering**"*; *"the **top-1 drug call**"*; *"cisplatin has ~157× the mean AUC and ~49× the IC50"* | **(d) PROHIBITED ×6** [D] | Six distinct prohibited constructions in one outward-facing document, several of them pre-registered as future acceptance criteria — i.e. commitments to re-assert the withdrawn claim. |
| `docs/DATA_ACCESS_PLAN.md:585` | *"a **validated** screen"* | **(d) PROHIBITED (wording)** [D] | Loose affirmative "validated". |
| `docs/outreach/02_erasmus_lamfers.md:154-157` | *"**Our panel ranking is decided by delivery × potency** … cisplatin has roughly 157× the delivered exposure and 49× the IC50 … they currently differ by about **4 % in predicted kill**"* | **(d) PROHIBITED ×2** [D] | A ranking asserted as current, **plus** an absolute kill difference as an output — two prohibitions in one sentence, in a draft addressed to an external collaborator. [D] adds the 4 % is unsourced anywhere in the repo and contradicted by both files carrying the quantity (`ranking_robustness.md:202` +0.1588; `real_drug_panel.md` 0.4362 vs 0.2595). |
| `docs/outreach/02:167`, `docs/real_data/pk_anchors.md:129-134, :144-149`, `docs/real_data/drug_response.md:116-119` | pre-register *"the predicted **drug ordering**"* / *"the sim should reproduce at least this **ordering**"* / a cross-model sensitivity ordering, as acceptance criteria | **(d) PROHIBITED ×4** [D] | `pk_anchors.md`'s criterion additionally **fails** on its own terms — `results/pk_multidrug.md` Test 2 reports ρ = 0.50. |
| `docs/QUICKSTART.md:11`, `docs/UI_ARCHITECTURE.md:323`, `docs/VERIFICATION_METHODOLOGY.md:798`, `docs/molecular_layer.md:398, :439`, `docs/real_data/ivygap.md:111-112`, `docs/real_data/drug_response.md:135-138`, `docs/outreach/README.md:47`, `docs/outreach/03:134` | loose affirmative *"validated"* | **(d) PROHIBITED (wording) ×9** [D] | `outreach/README.md:47` is the sharpest: it instructs drafters to *"lead with what has been built and **validated**"* — **thirteen lines after `:34` names "validated" as the word to police**. |
| `docs/outreach/03_controlled_access_skeleton.md:112-113` | tells a data-access committee the model predicts *"the resulting **local response**"* | **(d) PROHIBITED** [D] | Clinical-response framing, six lines before the same document concedes response predictions have never been tested. |

#### Contradicted by current code or the filesystem

| file:line | claim | verdict | evidence |
|---|---|---|---|
| **`docs/DEPLOYMENT.md:743`** | *"**Ten.** There is no config file; **these are the whole surface**, and only five are security-relevant."* | **(c) CONTRADICTED** | **Verified, and the delegated reading corrected**: the sweep reported this line as saying *"Six."* — it says **"Ten."** My own count: `grep -rhoE 'TMESIM_[A-Z_]+' tmesim/ app/ scripts/ \| sort -u` → **12**. Undocumented: `TMESIM_MAX_GRID_CELLS` (caps grid size and produces a client-facing refusal, `api.py:175,723,2276`), `TMESIM_BILLING_PROVIDER` (`account.py`, shipped **today**), `TMESIM_STUDIO_API_KEY_SOURCE`, `TMESIM_BAKED`, `TMESIM_HOST_MEMORY_MB`. A doc that says "these are the whole surface" is making a completeness claim, which is the kind that rots silently. |
| `docs/QUICKSTART.md:705` | `# kill_fraction: 0.0198` for the doc's own stated recipe (dx_um=80, seed 0, total_time 1800) | **(c) CONTRADICTED** [D] | Delegate ran the exact recipe: **0.0608**, a **3.07×** miss. A reproduction recipe that does not reproduce is the highest-value class of doc defect, because the reader executes it. |
| `docs/QUICKSTART.md:577` | `run history /path/to/tme-sim/app/server/.runs` | **(c) CONTRADICTED** [D] | `runs_dir()` → `/Users/samarthbasanth/tme-sim/app/runs`. |
| `docs/QUICKSTART.md:742` | *"the breast section and HGCC … `vivum sections` prints the exact commands for each"* | **(c) CONTRADICTED** [D] | `SECTION_RECIPES` holds 7 entries, **all GBM**; `vivum sections` prints nothing for either. |
| `docs/QUICKSTART.md:479, :534` | *"**13 of the 13 rows it scores never invoke the transport PDE**"* / *"Only the full run produces `PDE` rows at all"* | **(c) CONTRADICTED** [D] | `AnchorGroup(key='immune-sink', tier='quick', scope='PDE')` is in the live registry. Note this same figure is quoted approvingly at `README.md:466` and `docs/DATA_ACCESS_PLAN.md:25`. |
| `docs/REPRODUCE.md:435` | *"cli acquire — **17 s** (4/4 sources OK) — measured"* | **(c) CONTRADICTED — internally** [D] | `:216` of the **same file**: *"The previous claim in this file — '17 s' — **did not survive re-measurement**"*, giving 1 m 50 s / 3 m 40 s / 3 m 52 s. The §4 table was never updated. |
| `docs/REPRODUCE.md:154` | *"`molecules` lists the **10** library GBM drugs"* | **(c) CONTRADICTED** [D] | `len(GBM_DRUGS)` = **12**. |
| `docs/REPRODUCE.md:169`, `:329` | `data_cache/` *"holds **6.1 GB**"*; Tier B *"**5.3 GB**, 38 directories, ~1,950 files"* | **(c) CONTRADICTED** [D] | Measured **3.9 GB**; `data_cache/hunt` **360 MB / 24 dirs / 1,217 files**. |
| `docs/REPRODUCE.md:552` | *"across the **five-drug** panel"* | **(c) CONTRADICTED — internally** [D] | `:441` of the same file: *"panel expanded from **5 to 10** drugs."* |
| `docs/architecture.md:558` | *"The fix is `reconstruct.vessel_network` (**built, dormant**)"* | **(c) CONTRADICTED — internally** [D] | `:224-228` of the same file: *"**THIS IS THE SHIPPED DEFAULT** … that has been false since commit `ceb96c0`."* |
| `docs/architecture.md:664` | *"**13 of 13 scored rows never invoke the solver**"* | **(c) CONTRADICTED** [D] | Same `immune-sink` evidence. |
| `docs/MATHEMATICAL_MODEL.md:364-372, :366, :375` | §2.6 `P = Kp·Dm/d` as a single membrane leg; `Kp` clamped to `[1e-3, 1e5]`; `η_m = 0.10 Pa·s` | **(c) CONTRADICTED ×3** [D] | `permeability.py:216-220` is a **series** resistance with a 30 µm unstirred water layer, documented nowhere; `:66-67` `_LOGP_KP_MAX = 3.0` → `[1e-3, 1e3]`; `:48` `_ETA_MEMBRANE_PA_S = **1.0**`. |
| `docs/MATHEMATICAL_MODEL.md:934` | IFP smoothed with *"σ = **5 grid cells**"* | **(c) CONTRADICTED** [D] | `pressure.py` `_TUMOR_SIGMA_UM = 60.0` — a physical length, dx-independent. |
| `docs/MATHEMATICAL_MODEL.md:1394-1396` | `ρ_vasc` = *"**KDE** of measured ENDOTHELIAL cell positions"* | **(c) CONTRADICTED** [D] | `_DEFAULT_VASCULAR_MODE = "discrete"`. |
| `docs/MATHEMATICAL_MODEL.md:148, :1064` | MW 150 kDa → *"r_s ≈ **5.3** nm"*; `mitoxantrone` listed in the determinant table | **(c) CONTRADICTED — minor ×2** [D] | Computed **5.506 nm**; `get_molecule('mitoxantrone')` → unknown drug, so `vivum simulate --drug mitoxantrone` fails. |
| `docs/UI_ARCHITECTURE.md:200-202` | the READ rail section inventory | **(c) CONTRADICTED — incomplete** [D] | Omits `finding.js` (*"THE FIRST SCREEN"*), `reach.js`, `verdict.js`, `compound.js`, `account.js` and five more. The document describes a UI two rewrites old. |
| `docs/VERIFICATION_METHODOLOGY.md` — 14 rows incl. `:9`, `:27`, `:59`, `:409-410`, `:411-412`, `:425-428`, `:1112-1113`, `:1150-1153`, `:1158-1161` | *"Every number below is traceable"*; `penetration_index` *"definitionally unable to discriminate"* and *"withdrawn"*; *"the grid drift is 13× the entire across-drug spread"*; *"No tumour cell in any arm receives zero drug"* | **(c) CONTRADICTED ×14** [D] | Three load-bearing numbers trace to nothing (`33.4 µm`, the synthetic `0.998 → 0.000` ladder, the codebase counts). `penetration_index` is still computed and printed (`cli.py:239,252,1100,1106`) and `results/delivery_metrics.md:130` scores it **discriminating, separation 15.2**. The drift inequality has **inverted** — spread 0.206 now exceeds drift 0.0185 by 11×. `vessel_network_impact.md:239-247` records 0.07–0.15 % of cells at exactly zero. `:1158-1161` attributes the 4e-14 conservation to the oxygen solve; `numerical_rigor.md:23` attributes it to the **drug-transport** operator. |
| `docs/DATA_INVENTORY.md` — ~14 rows | *"6.1 GB"*; *"gitignored"*; *"every dataset on disk"*; Tier A *"~0.9 GB"*; A1/A2 marked **MANUAL**; `data_cache/ivygap/` *"36 MB + 89 MB"*; *"8 open STAR files"*; `verified_*` *"371 MB"* | **(c) CONTRADICTED ×14** [D] | **3.87 GB**; `git ls-files data_cache` returns **9 tracked files**; **13 datasets omitted** — including `xenium_breast_janesick/`, which `scorecard.py:2686,2707,2728` declares a hard `needs=` and which has **no URL, no bytes and no recipe anywhere in the file**; Tier A **3.51 GB**; `section_catalog.py:88,112` both **AUTO**; `data_cache/ivygap/` **does not exist** and `ivygap.load_structures()` raises `FileNotFoundError`; **7** STAR files; **69 MB**. |
| `docs/DATA_ACCESS_PLAN.md:75`, `:23` | *"`tmesim/validation/` (**24 modules**)"*; *"the synthetic-data era is deleted"* | **(c) CONTRADICTED ×2** [D] | **40 files**; `core/config.py:267-274` still ships hand-set "-like" library entries. |
| `docs/data_acquisition_plan.md` — whole file | documents `tmesim/data/loaders.py`, `load_real_cohort()` producing *"131 real records"*, and *"a **synthetic 2-D layout**"* as the current pipeline | **(c) CONTRADICTED — the most out-of-date file in the tree, and the only stale one with no banner** [D] | Neither module exists, and `scripts/reproduce_real.sh:177-180` treats their absence as a **pass condition**. `:58-66, :73-75, :121-126` describe the exact pipeline `REAL_BUILD_PLAN.md:4` retires. `:84-88`, `:104-112` on HTAN/PDMR are contradicted by the later-dated `DATA_ACCESS_PLAN.md:650-651`. |
| `docs/molecular_layer.md:161-166, :195-196, :250-260, :468-472` | `GBM_DRUGS` list; *"the self-test verifies … antibody ~1–10 µm²/s"*; the `build_profile` signature; *"computed from the molecule's **mass** alone"* | **(c) CONTRADICTED ×4** [D] | Omits `ribociclib` (12, not 11); `stokes_einstein.py:169-179` carries an explicit HONESTY NOTE **refusing** that figure and asserts `10 ≤ D ≤ 60` (computed 34.37); the documented signature **omits `two_route: bool = True`**, the shipping default-ON law; and the mass-only sentence is **self-contradicted three paragraphs earlier by its own §3 box**. |
| `docs/real_data/ivygap.md:23, :97-99, :103-106, :116-118, :131-139` | documents eight files in `data_cache/ivygap/` and derives mean fractional areas, subtype/MGMT counts and a marker signature table from them | **(c) CONTRADICTED — documents a cache that does not exist** [D] | Directory absent; only `ivy_gap/gene_expression_matrix.zip` (4 members). Four of the eight files are **neither on disk nor in the zip**, so all three derived tables are **unobtainable from anything the repo holds**, and **every line of the `:131-139` loader code block raises `FileNotFoundError` today**. |
| `docs/real_data/singlecell_spatial.md:304-316, :379-381` | *"**Known downstream gap (NOT fixed here)** … only two classes … **left open**"*; *"**No Xenium panel in this project measures MGMT**"* | **(c) CONTRADICTED — all four clauses** [D] | `marker_panels.py:137 IO_PANEL_MARKER_MAP` implements exactly the substitutes named unavailable, `select_marker_map()` picks it at coverage **1.00**, and `results/classifier_fix.md:33,193` scores the IO section with **four** classes. And `data_cache/gse324268_xenium/` is a 366-gene panel with **8 MGMT probe sets** — the dataset `docs/outreach/04` is entirely built on. |
| `docs/real_data/real_data_status.md:80` | `eff_D = 278 µm²/s` | **(c) CONTRADICTED — and it defeats the file's own banner** [D] | Live compile gives **42.42**. The `:5-8` banner flags four known-wrong numbers and **misses this fifth**, which makes `:84` (*"Re-running … reproduces this end to end"*) false. |
| `docs/real_data/pk_anchors.md:92`, `docs/real_data/drug_response.md:47, :124-131` | HTAN GBM spatial / NCI PDMR as the gated gold standard | **(c) CONTRADICTED** [D] | `DATA_ACCESS_PLAN.md:650-651` and the file's own sibling `outreach/03:34-35` (*"Nothing to request today" / "**Do not apply**"*). |
| `docs/REAL_BUILD_PLAN.md:6, :16` | *"Validate reconstruction against held-out real data"* as achieved; VEGFA/HIF1A as available markers | **(c) CONTRADICTED ×2** [D] | `results/heldout_patient.md:104` calls the held-out section *"the clearest single failure of generalisation in this report"*; both genes are **absent from the 319 panel**. |
| `docs/outreach/01:100-106`, and *"~2× geometric-mean fold error"* in **all four** outreach drafts | *"the model was fixed"*; the 1.995× headline | **(c) CONTRADICTED** [D] | 1.995× requires the FIX-2 barrier leak, which `treatment/drugs.py:141,159` ships **OFF** and which **no `app/` call site passes**. The structure-only default-path number is **3.451×**. Four external-facing drafts quote a number the shipped default does not produce. |

#### Superseded

| file:line | claim | verdict | evidence |
|---|---|---|---|
| `docs/architecture.md:695`, `docs/MATHEMATICAL_MODEL.md:351`, `docs/REAL_BUILD_PLAN.md:26`, `docs/DATA_ACCESS_PLAN.md:28`, `docs/molecular_layer.md:370, :397, :480`, `docs/QUICKSTART.md:255` | **R² = 0.718** asserted flat | **(b) SUPERSEDED ×8** | `results/hbd_convention_audit.md:3-7`, `:267-270` — interval **0.498–0.839**. QUICKSTART alone carries the `descriptor_falsification` axis (`:274-292`, R² → ≈0.27) but not the hbd axis. |
| `docs/architecture.md:639, :641`, `docs/MATHEMATICAL_MODEL.md:1327, :1334, :1342`, `docs/REPRODUCE.md:533, :565`, `docs/VERIFICATION_METHODOLOGY.md:862, :878` | **α = 148.6** | **(c) CONTRADICTED ×9** | Live `calibrate_alpha()` = **22.223477711172563**. See §2.7. |
| `docs/architecture.md:299, :467, :705`, `docs/REPRODUCE.md:413`, `docs/MATHEMATICAL_MODEL.md:707` | *"Report a **SET**, not a winner"*; *"drug **rankings** are unaffected"*; *"preserves … **relative ranking**"* | **(b) SUPERSEDED ×5** | *"Report a SET"* was the 2026-08-09 position; post-2026-08-21 there is no reportable top-of-panel at all, set or winner. |
| `docs/architecture.md:490-492, :538, :539, :540-541`, `docs/REPRODUCE.md:552-555`, `docs/VERIFICATION_METHODOLOGY.md:331-343, :400-407, :422-428, :442-445, :455-457, :482-488, :1071-1073` | the entire `delivery_metrics` number set — `penetration_index` 0.194–0.196, fitted L 26–28 µm, source null 23.7 µm, coverage 0.311–0.925, separation 21.9, 120× exposure spread, 30,174 cells, 9,216 re-rankings | **(b) SUPERSEDED ~20 rows** [D] | `results/delivery_metrics.md` now gives 0.3044–0.5107; L **280.0 / 78.0 / 150.2 / 53.4 / 45.1 µm**; source null **nan**; coverage **0–0.9955**; separation **67.7**; **158×**; **23,812** cells; **10,240** re-rankings. Two documents are quoting an entire superseded results table between them. |
| `docs/VERIFICATION_METHODOLOGY.md:20-21` | *"55 tests green at the time; the suite now collects **114**"* | **(b) SUPERSEDED** [D] | `pytest --collect-only -q` → **1430 collected**; a filtered run gave **1296 passed, 0 failed**. The doc understates by ~12.5×. |
| `docs/VERIFICATION_METHODOLOGY.md:560-578, :835-854, :940` | the MKI67 story (2.22 % vs 0.444 %, 5.0×, p = 0.0078); *"a **165×** swing"*; VIABLE 0.60–98.19 %; HYPOXIC 13.43 % | **(b) SUPERSEDED** [D] | `results/state_validation.md:303-305` now reads **nan% / 0 of 0 / inf× / n = 0 vs 23,812**; `:270-277` VIABLE **0.00–10.04 %**, HYPOXIC **0.00 % in all six rows**. |
| `docs/DEPLOYMENT.md:328-329` | *"`data_cache/` (**5.3 GB**), `app/runs/` (**647 MB**)"* | **(b) SUPERSEDED — by today's prune** | Measured today: `data_cache` **3.9 GB**, `app/runs` **709 MB** (down from ~1.49 GB; the prune freed **782,579,352 B**). |
| `docs/QUICKSTART.md:238`, `docs/MATHEMATICAL_MODEL.md:828, :1131` | *"`penetration_index` … **the two summary numbers**"*; *"$A_{ref}$ is the **converged** production grid"*; *"the penetration index is used as the **primary signal**"* | **(b) SUPERSEDED — internally ×3** [D] | Each is contradicted later **in its own file**: `REPRODUCE.md:549` (*"RETIRED … **ignore it**"*), `MATHEMATICAL_MODEL.md:1639-1641` (*"there is **no Δx** at which this pipeline has been shown to be grid-converged"*), `:1642` (*"retired regardless"*). |
| `docs/QUICKSTART.md:395, :457` | *"`--quick` covers **8 of the 16** registered anchor groups"*; the 18-failure inventory | **(b) SUPERSEDED** [D] | `build_registry()` → **17** groups, **9** quick. The inventory omits `brain-slice-tissue-capacity` (*"execution record, **FAIL**"*) and `ifp-wick-in-needle` (*"breast **FAILS ~40x**"*). |
| `docs/real_data/pk_anchors.md:46-48, :136-140` | *"the sim predicts **L ≈ 796 µm**, still ~11-14× too large"*; *"the current engine **over-predicts**"* | **(b) SUPERSEDED** [D] | `results/pk_multidrug.md` gives doxorubicin **64.9 µm**, inside the band; out-of-sample is **2/2 PASS**. |
| `docs/DATA_ACCESS_PLAN.md:27, :112-118, :310-315, :599-612` | *"the closed-form α calibration **survives**"*; the 0.241 / 0.231 / 0.582 figures; four action items | **(b) SUPERSEDED** [D] | `results/alpha_identifiability.md`: α is *"a **gauge coordinate**, not a measurement"*. The three figures **appear nowhere in `results/ranking_robustness.md`**, which reports 0.931. The action items are already complete. |
| `docs/outreach/02:125-128`, `docs/outreach/01:95-97`, `02:109-110` | *"the transport-level checks **have held** … that asymmetry is the honest summary"*; *"3.65× predicted"* | **(b) SUPERSEDED** [D] | `pk_ranking_robustness.md` §0 and `capacity_ranking_robustness.md` §0; the shipping model gives **17.20×**. |

#### Notable passes in `docs/`

| file:line | claim | verdict | evidence |
|---|---|---|---|
| `docs/QUICKSTART.md:219-225, :693, :711` | the CLI transcript — `effective_D 42.425 · ic50 606.047 / 1.200 · binding 4.428e-04 · clearance 5.794e-05 · bbb 0.125`; `cells: 40887 tumor: 23812 grid: (65, 103)`; IC50 to 16 digits | **(a) STILL TRUE** [D] | Recompiled live: **every digit** matches, including `606.0466722391783` and `0.12537158358322562`. |
| `docs/REPRODUCE.md:56-71, :361-368` | 15 pinned package versions; 11 manual curl byte counts | **(a) STILL TRUE** [D] | `constraints.txt` and `pip list` match on all 15; every byte count matches `section_catalog.py` and `ls -la`. |
| `docs/REPRODUCE.md:150-157` | the `delivery_metrics --selftest` expected output | **(a) STILL TRUE** [D] | Ran it — output matches `:156-157` exactly. |
| `docs/MATHEMATICAL_MODEL.md:1659-1684` | the 40-row appendix mapping every equation to a `module::function` | **(a) STILL TRUE** [D] | Programmatic `hasattr` over all 40 → **all refs resolve**. |
| `docs/MATHEMATICAL_MODEL.md:66-68, :773-776, :927-931, :1040-1052` | `DrugParams`, `OxygenConfig`, `pressure.py` and `PDConfig` default constants | **(a) STILL TRUE ×19** [D] | Every one matches the live dataclass. |
| `docs/UI_ARCHITECTURE.md:56-58, :288-295, :339, :373-375` | 9 height tokens, 6 font sizes, 8 spacing steps, `--radius: 2px`, the `overflow: hidden` contract | **(a) STILL TRUE ×25** [D] | `tokens.css` and `base.css` — exact, comments included. |
| `docs/DEPLOYMENT.md:74-83` | *"Sections are registered **by path, not by upload** … no HTTP upload route for tissue"*; *"the server makes **no outbound network requests**"* | **(a) STILL TRUE** | `grep -rn 'requests\.\|urllib\|httpx\.\|aiohttp' app/server/*.py` → **exit 1**; fonts vendored. |
| `docs/RUNBOOK_DEPLOY.md:142-144, :176-184, :283-285` | 7/7 sections available; 38 sections / 38 available; `/geometry` **368,003 B**, `/field/*` **1,144,828 B** | **(a) STILL TRUE** [D] | Byte counts exact (1,144,828 = 424×675×4+28); 7 built-in + 31 ingested = 38. |
| `docs/real_data/visium_gbm.md`, `docs/real_data/scref_gbm.md` | all 20 per-sample numbers; 23,460 × 3,589 and seven cell-type counts | **(a) STILL TRUE — zero defects in both files** [D] | Fully reproduced. |
| `docs/outreach/04_holzel_gse324268.md` | 851,509 cells; MGMT 15,191 counts / 14,675 cells / 1.72 %; design_id `ECQNGH`; 8 probe sets; 366 targets | **(a) STILL TRUE — the best-evidenced document in the tree** [D] | Every number reproduces exactly from the local `.h5`, discharging the file's own `:117-120` instruction to verify before sending. |


### 2.4 `landing-3d/**` copy

| file:line | claim | verdict | evidence |
|---|---|---|---|
| `landing-3d/app/layout.tsx:171` | the page `<meta name="description">` ends *"Built for pharma, biotech, and research teams **comparing candidate therapies**."* | **(d) PROHIBITED** | Cross-compound comparison framing, on the highest-visibility surface the page has — search results and link unfurls. Directly against the guardrail written into the sibling file: `landing-3d/lib/copy.ts:146-149` — *"No ordering, ranked list, shortlist, 'top pick', 'which compound wins', pairwise call or **cross-compound comparison** may appear on this page **in any form, hedged or not, positive or negative — including as a limit**."* `layout.tsx` does not import `copy.ts`, so nothing mechanically couples the guardrail to the string it governs. |
| `landing-3d/lib/copy.ts:1030, 1032` | ledger row `state: "14× TOO SHALLOW"`, detail *"about 14× shallower than measured, and the error flatters the drug"* | **(b) SUPERSEDED** | `results/kyle2004_penetration.md:504-526` (2026-08-22), cited by `copy.ts:1011-1017` as this row's own source: the 14× *"is not a length error"*; λ is **2.45×** too long; *"A 2.45× λ error is **~2.2–3.0× on the map**, not 14×."* The landing page prints the amplified estimator figure flat. Aggravating: the source file's own conclusion (`:740`) re-asserts 14×, so the landing page is faithfully copying a number its source corrected but did not adopt. |
| `landing-3d/lib/copy.ts:902-905` | ledger row `state: "0.498–0.839, NOT 0.718"` with *"0.718 was scored on the paper's donor counts, in a convention the model does not use"* | **(a) STILL TRUE — and it is the only surface that carries this** | Matches `results/hbd_convention_audit.md:3-7`, `:267-270` exactly. Every other surface (§2.1, §2.3, §2.5) still asserts 0.718 flat. |
| `landing-3d/lib/copy.ts:564`, `:899` | *"how much drug the tissue can store… cannot move a converged field at all — checked to one part in a billion… It does move what a finite exposure window accumulates"* | **(a) STILL TRUE** | Matches `app/server/reach.py:391-396` and `results/capacity_ranking_robustness.md` §2b. The steady/finite distinction is carried correctly, and the two halves travel together in both places they appear. |
| `landing-3d/lib/copy.ts:920-929` | `Absolute kill fraction · NOT DEFENSIBLE`; `Patient response · OUT OF SCOPE` / *"Not predicted. Not for diagnostic or treatment-planning use."* | **(a) STILL TRUE** | Every prohibited category negated explicitly on the rendered surface. |
| `landing-3d/lib/copy.ts:967` | *"the solve is a single process on your own machine, and it sends no tissue or compound data anywhere"* — no occurrence of "upload" anywhere in rendered copy | **(a) STILL TRUE** | Narrowed to a software property; the withdrawn firewall claim did not return. |
| `landing-3d/lib/copy.ts:688` | *"A bench assay on that section could refute it; **a ranking never could**."* | **borderline — recorded, not counted** | The only occurrence of the word "ranking" in rendered copy. It is a contrast to a thing the product does not do, not an offer — but `copy.ts:148-149` bans the word *"positive or negative — including as a limit"*, and `:149-151` says a fossil row *"invites a reader to believe the feature is still there and merely caveated."* Judgment call. |
| `landing-3d/components/ui/MetricOverlay.tsx:77-78` | the panel's summary row renders `Spread` → **`2.0x`**; `1.08–1.21×` appears nowhere on the panel | **(b) SUPERSEDED — soft** | `copy.ts:908-911` registers the verdict as *"1.08–1.21×, NOT 2.0×"*. The instrument panel gives typographic primacy to the retracted headline and leaves the correction to a ledger in a different beat. All three numbers on the panel are real and the footnote is honest; the emphasis is not. |
| `landing-3d/components/ui/MetricOverlay.tsx:33` (comment) | *"~1.13-1.20x, not 2.0x"* | **(c) CONTRADICTED** | `copy.ts:443-444` and the rendered copy both say **1.08–1.21×**. Comment-only, no visitor sees it, but two files that must agree do not. |
| `landing-3d/README.md:21` | *"This project uses `next/font` to automatically optimize and load **Geist**"* | **(c) CONTRADICTED** | Unmodified `create-next-app` boilerplate. `app/layout.tsx:78-80` records Geist as *deliberately not chosen*; `:125-152` loads Archivo + Atkinson Hyperlegible Next + Atkinson Hyperlegible Mono. Carries no product claim, so low severity — but it is a false statement in a file named README. |
| `landing-3d/lib/timeline.ts:18`, `:299`; `components/canvas/Markers.tsx:27`; `app/globals.css:11`; `components/canvas/ConcentrationField.tsx:77` | *"AI Analysis"*, *"Treatment Planning"*, *"AI flagged region"* | **(a) STILL TRUE** | All are **comments recording their own removal**. `lib/timeline.ts:22-29`'s `SECTIONS` labels reach no DOM (`ProgressRail` takes a `count`, not labels). **Not violations.** |
| `landing-3d/` rendered numbers | 41×, 99.8 %, 12.8×, 10.5×, 0.43/0.86/2.0×, 82 %, 1.7, 40,887, 8.09 × 5.07 mm, 58.2/19.5/11.8/10.5 %, 250 µm, 3.2–4184, n = 33, 68 compounds | **(e) UNVERIFIED individually** | Extracted and catalogued (96 visitor-visible rows) but not each re-derived. Two were checked and hold: `40,887` (`tmesim/reconstruct/real_xenium.py:4`), `3.2–4184 mL/g, n = 33` (`app/server/reach.py`, `results/brain_slice_anchor.md`). No speed or cost claim appears anywhere on the page — swept for `48.6`, `60x`, `4.4s`, `250s`, `1.74`, `1.56`, `4.77`, `GMFE`: **zero hits in rendered copy**. |

### 2.5 `app/server/report.py` and `app/web/js/**`

| file:line | claim | verdict | evidence |
|---|---|---|---|
| `app/web/js/ui/provenance.js:92` | *"Scored BLIND and UNFITTED against 71 compounds… **R^2 = 0.718**, 99% within 3x of measurement, and 24 of 31 exact-mass pairs ordered correctly"* | **(b) SUPERSEDED** | `results/hbd_convention_audit.md:3-7`, `:267-270`. This is the Studio's live provenance panel — the place a user goes to ask what backs the number on screen. It does cite `descriptor_falsification.md` and does say *"ordering COMPOUNDS by it, or by anything else here, is withdrawn"* (both good); it does not carry the hbd interval. |
| `app/web/js/ui/spec.js:892` | *"R^2 = 0.718 blind against 71 measured compounds"* | **(b) SUPERSEDED** | Same. |
| `app/server/report.py:335` | *"reaches R^2 = 0.718 with Spearman +0.821, 99 % of compounds within…"* | **(b) SUPERSEDED** | Same — and this is the **generated therapy report**, the artefact a user exports and sends onward. |
| `app/server/report.py:353` | *"produces it is scored, blind and unfitted, at R^2 = 0.718 across…"* | **(b) SUPERSEDED** | Same. |
| `app/server/api.py:798` | *"reaches R^2 = 0.718 blind"* | **(b) SUPERSEDED** | Same (API-served copy). |
| `app/web/js/mock.js:546` | `alpha: 148.6 + rng() * 9` | **(c) CONTRADICTED** | Live α = 22.223477711172563. Mock/demo data, shown to a user in mock mode. Low severity, but it is a number on a screen. |
| `app/web/js/ui/reach.js:638-646` | the **`DELIVERY SPREAD`** headline — `num(p90 / p10, 2) + 'x'`, *"the best-perfused tenth of tumour cells receives this many times the exposure of the worst-perfused tenth… **no dose, no potency, no threshold to choose**"* | **(b) SUPERSEDED — and this is the #1 live claim** | `grep -n -i capacity app/web/js/ui/reach.js` returns capacity caveats at `:850`, `:863`, `:874` (the **falloff** curve) and `:1341`, `:1353`, `:1360` (the **ladder**). **None is attached to `E_p90/E_p10`.** Its own pre-registered capacity sweep failed both criteria on 2026-08-24 (§2.2) and the published `.md` still says the sweep has not been run. Nothing in `README.md`, `docs/`, `landing-3d/` or the UI cites `p90_capacity_sweep`. |
| `app/web/js/ui/reach.js:614-628` | the dated `CORRECTION 2026-08-22` retracting the *"~10× grid sensitivity"* claim to a measured 0.6 % | **(a) STILL TRUE** | Matches `results/reach_invariances_measure.md:616-621`. |
| `app/web/js/ui/reach.js:850-886`, `:1339-1370` | three-arm, mode-aware capacity caveats (`UNRECORDED` / `STEADY` / finite-window) | **(a) STILL TRUE** | Matches `app/server/reach.py:380-470`. The steady-state cutover, the capacity absence and the finite-window residual are all carried correctly, per run, with the `window_mode`-absent case refused rather than defaulted. **This is the best-maintained claim surface in the repository.** |
| `app/web/js/ui/finding.js:345-376` | *"58% of the panel (7 of 12; bevacizumab excluded, its source is identically zero). Up to 57% of tumour cells change their per-cell label."* + the stable-compound rule + *"See results/verdict_inversion_census.md"* | **(a) STILL TRUE** | Matches `results/verdict_inversion_census.md` §0 and §2 exactly, including the per-run branch for a compound not in the census (`:351-353`: *"absence from that census is not evidence of stability"*). |
| `app/web/js/ui/welcome.js:52` | *"Nothing here is a clinical or treatment-planning tool."* | **(a) STILL TRUE** | Correct negation. |
| `app/web/js/main.js:1123-1127` | comment: *"A run without a time-lapse is ordinary; **the reason already rides in `RunResult.timelapse.reason`** and is shown in the run's provenance panel."* — the `catch` block swallows the failure and calls `scrubber.clear()` | **(c) CONTRADICTED as of today** | For the 51 runs measured in §2.6, `RunResult.timelapse.reason` is `null`. `main.js:1118` gates on `info.available === false`, which is `false` here because `available` is `true`, so the code proceeds, receives **410 Gone**, and discards the server's correct explanation without surfacing it. The 410's `reason` field — the one thing `scripts/prune_run_payloads.py` was written to produce — never reaches the user. |
| `app/server/report.py:2640-2643` | *"…on the panel this document used to rank it **moved the top-1 pick outright, cisplatin → vincristine**, swapping 4 of 45 pairs."* | **(d) PROHIBITED** | The **generated therapy report** — the artefact a user exports and forwards — prints a named cross-compound top-1 ordering. The **same document's own limit 2**, `report.py:224-228`, refuses to do exactly this: *"which two compounds traded places **is NOT restated here, because naming them would print the very ordering this page has withdrawn**."* The hero hedges it as *"quoted here as EVIDENCE OF THE LAW'S LEVERAGE and not as a comparison"* (`:2643-2645`), but the guardrail admits no hedge — `landing-3d/lib/copy.ts:148`: *"in any form, hedged or not."* One document, two incompatible policies, ~2,400 lines apart. |
| `app/server/report.py:2661-2665` | *"with **P(cisplatin = rank 1) at 0.931** under that report's declared priors falling to **0.587** in its stress scenario"* | **(d) PROHIBITED — lesser** | Framed as history (*"This sentence used to appear here as a caveat on a RANKING"*, `:2661`) and closed with *"The ranking is gone"* (`:2666`), which is the correct narration pattern. It still puts a named compound at rank 1 on the rendered page. Weaker than the row above; recorded, and counted once with it in §0 rather than twice. |
| `tmesim/molecular/byo.py:222-223` — rendered verbatim in `app/web/js/ui/compound.js` (`bulletList`, ~`:307`) | *"**Read the ORDER over a mass-matched series as weak evidence**, not as a measurement, and read no absolute diffusivity from it at all."* | **(d) PROHIBITED** | The same construction as `README.md:513-515`: grades a cross-compound ordering as *weak evidence*, i.e. as evidence. Directly against `app/web/js/ui/provenance.js:92`, rendered in the same Studio — *"ordering COMPOUNDS by it, or by anything else here, **is withdrawn**"* — and against `app/web/js/ui/compound.js:352` / `spec.js:886-888` (*"do not set it beside another compound"*). `byo.py:648` compounds it: *"read this before you rank on it."* |
| `app/web/js/ui/compound.js:139-141` | *"In particular: **diffusivity is a function of molecular weight alone**, and that has been measured against 71 compounds and refuted, so a series of similar-mass analogues is the one case this model is known not to resolve."* | **(c) CONTRADICTED by current code** | The engine ships the **two-route** law, not MW-only. `tmesim/molecular/two_route.py` is on the compile path; `compile_drug.build_profile` routes `logp` and `hbd` into `effective_D` (`results/CLAIM_INVENTORY.md:11-12`). Contradicted by `provenance.js:92`, `report.py:335-337`, `spec.js:891-892` — and decisively by `tmesim/molecular/byo.py:216-219`, whose text renders **inside this same panel, directly beneath this paragraph**: *"a two-route cell term… R^2 = 0.718… The MW-only model this replaced scored R^2 = 0.022."* The panel states two contradictory facts about its own compiler, a few lines apart. The paragraph describes the **falsified predecessor** as the shipping model. |
| `app/web/js/ui/inspect.js:386-397` renders `penetration_length_um`; `:382` (comment) records *"measures as ~14x too shallow (results/kyle2004_penetration.md sec. 6)"* | the rendered length carries no such qualifier | **(b) SUPERSEDED** | Two defects at once: the qualifier lives only in a comment no user sees, and the qualifier itself is the figure `results/kyle2004_penetration.md:504-526` corrected on 2026-08-22 (λ is 2.45× long; ~2.2–3.0× on the map). Same class as `landing-3d/lib/copy.ts:1030`. |
| `app/web/js/ui/history.js:227` | *"compounds this product ranked spans **113x**, while the published ordering survives a spread of only 1.53x"* | **(c) CONTRADICTED — minor** | Every other surface says **113.1×**: `app/server/report.py:219`, `app/web/js/ui/provenance.js:82`, `app/web/js/protocol/stub.js:55-56`, `app/web/js/editor/templates.js:126`, `README.md:31`. One surface rounding a shared constant differently is how two numbers become two claims. |
| `app/web/js/ui/account.js:21-23` (module header) | the panel *"repeats the server's own word for what those numbers are"* | **(c) CONTRADICTED — introduced today** | It does not. `account.js:404-407` **re-authors** the inventory sentence in the client; `usage.kind` (`app/server/account.py`, `"current inventory, not a meter"`) is never read. `usage.reconstruction_cache` is likewise served and never rendered. A caveat authored in two places is corrected in one — the defect `app/web/js/ui/verdict.js` exists to prevent, reintroduced on the day the panel shipped. |
| `app/web/js/mock.js:616` | `run.runtime_s = run._dur / 1000 * 4.1` — a synthesised value rendered in the `runtime` field | **(c) CONTRADICTED — minor** | Mock mode is fenced loudly (`app/web/css/components.css:1431`: `content: "SYNTHETIC · NO BACKEND · js/mock.js · NOT DATA"`), which is why this is minor rather than serious. It is listed because a fabricated runtime and a fabricated α (`mock.js:546`) sit in the same file and only the file, not the fields, is labelled. |
| `app/web/js/ui/finding.js:360-362` | *"On this section the stable ones are **VINCRISTINE, DOXORUBICIN, METHOTREXATE and ETOPOSIDE**."* | **borderline — recorded, not counted** | The subject is verdict **stability** under a missing term, not efficacy, and it is a faithful rendering of `results/verdict_inversion_census.md` §4. But it is a named multi-compound list on the claim surface, where `report.py:224-228` establishes that naming compounds in a set is itself the harm. Judgment call. |
| `app/web/js/ui/account.js:42`, `:137`, `:150` | tab label *"Plan & billing"*; subtitle *"access, API and **billing** for this deployment"*; button tooltip *"…what this server refuses, and **billing**"* | **borderline — recorded, not counted** | All three name billing as something this deployment *has*, and all three are read **before** the disclaimer at `:437` is reachable. A user who never opens the tab has seen *"Plan & billing"* and nothing else. The tab body is unmissable once entered, so this is framing, not a false claim. |
| `app/web/js/ui/tissue.js:272` | *"Nothing is uploaded from this browser for tissue."* — with `:269-271` scoping that the protocol panel *does* upload drug files | **(a) STILL TRUE** | The "upload" prohibition is honoured and the one legitimate exception is disclosed rather than elided. `welcome.js:107` and `tissue.js:193` agree. |
| `app/web/js/ui/spec.js:845-860` | *"this repository records exactly one end-to-end Studio run time — **47.13 s**… That is one machine and one cold start, **not a benchmark**"* and, for any other spec, *"the number does not transfer… **this repo records no measured relation between dx and run time**, so none is estimated here"* | **(a) STILL TRUE** | Matches `docs/DEPLOYMENT.md:267-280`, including that file's own record of retracting *"in 4 s on 40,887 cells"*. **No rendered string anywhere claims the steady-state cutover made runs faster or cheaper** — swept across all 59 `app/web/` files and `report.py`. The 4.2 %-slower result is honoured by omission. |
| `app/server/report.py:144-149`, `:2129`, `:2595-2598`, `:2896-2898` | `NOT_FOR_CLINICAL_USE` (*"Do not use it to select, rank, exclude, or dose any therapy for any patient"*); banner top and bottom; the four chips `RELATIVE EXPOSURE ONLY / UPPER BOUND / LAW-DEPENDENT / NOT A RECOMMENDATION`; `<title>` carrying NOT FOR CLINICAL USE | **(a) STILL TRUE** | The exported report is fenced on four independent surfaces including its own filename-visible title. |
| user-visible strings in `app/web/js/**` and `app/server/report.py` | | **(a) STILL TRUE — as a class** | Exhaustively swept across all 59 `app/web/` files and `report.py`: **no** `digital twin`, **no** personalised/personalized medicine, **no** AI/ML/deep-learning/neural copy, **no** positive standalone use of *"validated"* (every rendered instance is a negation), **no** "upload" applied to tissue. `report.py` carries ~60 hardcoded measurements; three mechanical weak points found — a hardcoded `45.0` divisor (`report.py:2409`), a hardcoded `10` µm bin fallback (`reach.js:1252`), and the `113x`/`113.1x` split above. |

### 2.6 Today's changes — `account.py`, `account.js`, `pruned_payloads`, the prune script

**The billing question, answered directly: no, `account.js`'s billing copy cannot
be read as a live billing capability.** It is the strongest disclosure surface
added to this repository in the period audited.

| file:line | claim | verdict | evidence |
|---|---|---|---|
| `app/web/js/ui/account.js:436-439` | `div.caveat.crit` rendered **first**, before any table: *"**BILLING IS NOT ENABLED ON THIS DEPLOYMENT.** No charge is made, no payment method is stored, and no usage is metered."* | **(a) STILL TRUE** | Backed by `app/server/account.py::_billing()`: `"enabled": False` unconditionally; `plan`, `seats`, `period` all `None`; `metered_quantities` and `invoices` empty lists. The `enabled` flag is not derived from the provider variable, so naming a provider cannot flip it. |
| `app/web/js/ui/account.js:441-451` | the billing table: `status`, `provider`, `plan`, `seats`, `billing period`, `metered quantities`, `invoices` | **(a) STILL TRUE** | Every row renders `NULL_GLYPH` for an absent value rather than a plausible default. An empty slot under an unmissable "not enabled" banner is not a capability claim. |
| `app/web/js/ui/account.js:453-458` | *"**A PROVIDER IS NAMED BUT NOTHING IS WIRED TO IT.** TMESIM_BILLING_PROVIDER is set on this host, and no implementation is mounted. Nothing is charged."* | **(a) STILL TRUE** | Matches `account.py::_billing()`'s `"configured-but-not-implemented"` branch. The one state that could be mistaken for working billing is the one called out loudest. |
| `app/web/js/ui/account.js:461-469` | *"What has to exist before usage can be billed"* + *"These are **preconditions, not a schedule**."* | **(a) STILL TRUE** | Matches `account.py`'s `blocked_on`: per-caller identity, and a durable usage record. Neither is presented as forthcoming. |
| `app/web/js/ui/account.js:403-407` | *"**THIS IS AN INVENTORY, NOT A METER.** These are the things that exist on this server right now. Deleting a run makes the number go down. Nothing here accumulates, nothing is billed from it."* — rendered **before** the numbers | **(a) STILL TRUE** | Matches `account.py`'s `"kind": "current inventory, not a meter"`. The framing precedes the figures, which is the only ordering under which it works. |
| `app/web/js/ui/account.js:419-422` | *"**NO RUN CAN BE ATTRIBUTED TO A PERSON.**"* | **(a) STILL TRUE** | Matches `account.py`'s `attribution` string. |
| `app/web/js/ui/account.js:137` | panel subtitle: *"access, API and **billing** for this deployment"* | **borderline — recorded, not counted** | A subtitle naming billing as something the panel is *for*. Refuted immediately and loudly by the tab itself. Recorded only because a subtitle is read before a tab body. |
| `app/server/api.py:2377-2389` | the `410 Gone` path: *"this run **DID** produce {label}; the file was deleted on {at} to reclaim disk ({bytes} bytes)… **re-run the RunSpec in this record**"* | **(a) STILL TRUE** | Correct, and it reads `pruned_payloads` **before** the fall-through branches, exactly as `scripts/prune_run_payloads.py`'s docstring says it must. `API_CONTRACT.md:237-245` documents it accurately. The server-side half of today's work is right. |
| `scripts/prune_run_payloads.py` — citation protection | *"Any run id named anywhere in `results/`, `archive/`, `docs/` or `README.md`… a reader following a citation must find what the citation points at."* | **(a) STILL TRUE** | Measured: 54 runs pruned, 19 still holding `timelapse.bin`, 96 runs total. The `--keep-recent` and cited-run protections both held. |
| **51 of 54 pruned runs** — `app/runs/*/result.json` | the run document simultaneously carries `pruned_payloads.timelapse = {at: "2026-08-30", bytes: …}` **and** `timelapse = {"available": true, "n_frames": 24, "reason": null, "payload": "GET /api/runs/{id}/timelapse (TMSL)"}` | **(c) CONTRADICTED by the filesystem, introduced today** | Measured with `.venv/bin/python` over all 96 run documents: **54** carry a pruned `timelapse`; of those, **51** still report `timelapse.available === true`. Examples: `r_00011`, `r_00012`, `r_00013` — each `n_frames: 24`, `reason: null`, pruned `2026-08-30`. `app/server/jobs.py:1541-1544` (`result_view`) copies `self.body` wholesale into the response, and both `timelapse` and `pruned_payloads` are in the persist allowlist (`jobs.py:1650-1666`). `result_view`'s `annotate` path backfills `delivery_law` and `claim_set` but **nothing reconciles `timelapse.available` against `pruned_payloads`**. So `GET /api/runs/{id}` asserts, in one payload, both that the time-lapse is available and that it was deleted. **This is the answer to "check that nothing now claims those payloads are present": something does, on 51 runs, in the API's own response body.** |

### 2.7 Cross-cutting: one constant, nine stale citations

`calibrate_alpha()` returns **22.223477711172563**. Measured live:

```
$ .venv/bin/python -c "from tmesim.validation.calibration import calibrate_alpha; print(repr(calibrate_alpha()))"
22.223477711172563
```

What moved is the constant, not the physics: the e-fold *target* is unchanged and
the calibration still lands it (`L_dox` = 64.921 µm), and the out-of-sample
predictions survive the change (mitoxantrone **48.3 µm** against a published
48.2). The two-route promotion shrank doxorubicin's `D_eff`, and
α = (D_eff/L² − k_cl)/k_b scales with it.

Stale at: `docs/MATHEMATICAL_MODEL.md:1327, :1334, :1342`;
`docs/architecture.md:639, :641`; `docs/REPRODUCE.md:533, :565`;
`docs/VERIFICATION_METHODOLOGY.md:862, :878`. Also
`results/delivery_metrics.md:9` (§2.2), `API_CONTRACT.md:101, :745`, and
`app/web/js/mock.js:546`.

**The sharp edge:** several of these files carry a correction notice reading
*"any report still quoting **309.6** is stale"* (`docs/architecture.md:641`,
`docs/REPRODUCE.md:566`). The notice is now itself one revision behind — a
staleness warning that has gone stale is worse than none, because it signals the
number has been checked.

---

## 3. The most dangerous LIVE claims, ranked

Ranked by *what a reader would do wrongly because of it*, not by how wrong it is.

**1. `results/p90_capacity_sweep.md:3` — a pre-registered FAIL/FAIL on the
product's headline delivery number, unpublished, while the record says the sweep
has not been run.**
The `.md` says *"Status: PRE-DECLARED, NOT YET RUN"* and its §2 is empty.
`results/p90_capacity_sweep.json` (2026-08-24 23:33) records
`verdict: "THE HEADLINE IS CAPACITY-SENSITIVE AND SO IS THE CLAIM BUILT ON IT"`,
`C1.fail: true` (temozolomide moves **5.87×** its own full nuisance-grid width;
doxorubicin **4.49×**), `C2.fail: true` (separation **1.5039** at R=39 against
the pre-declared floor 1.74). The subject, `E_p90/E_p10`, is rendered live as
`DELIVERY SPREAD` at `app/web/js/ui/reach.js:638`, and that block carries no
capacity caveat while the falloff and the ladder beside it carry three each.
This is the exact shape of the 2026-08-21 capacity-anchor failure the project
already withdrew a claim over — an anchor that failed, then went quiet.

**2. `docs/DEPLOYMENT.md:15-17` — the withdrawn claim, stated as what the product
does, in the scope paragraph of the deployment document.**
*"**Vivum reconstructs a tissue section and ranks compounds against it**; it
performs no diagnosis and produces nothing that is intended to select therapy for
a patient."* Verified verbatim. This ranks second because of **who reads it**: §0
of `DEPLOYMENT.md` is what an operator, a security reviewer or an institutional
data-access committee reads to decide whether they may run this at all — and the
sentence bounding the product's scope names ranking as the thing it does. The
prohibition is absolute (`README.md:34-36`).

**3. `docs/VERIFICATION_METHODOLOGY.md:413-417` — two complete cross-compound
orderings, printed in full, one cross-validated against the other.**
*"Ranked by `penetration_index` the panel comes out `cisplatin > lomustine >
temozolomide > methotrexate > doxorubicin`. Ranked by the replacement metric it
comes out `cisplatin > doxorubicin > lomustine > methotrexate > temozolomide` —
which is the *identical* order the simulated kill fraction produces."* Verified
verbatim. Every other ranking violation on this surface is a hedge, a fossil or a
framing; this one **is the output**. It is also stale — `results/delivery_metrics.md:98-104`
now orders TMZ first — so it is simultaneously prohibited and wrong.

**4. `README.md:509-511` — the withdrawn sentence, live, under "The honest
summary".**
*"Vivum is a research instrument whose defensible output is the relative
comparison of compounds on one tissue section."* `README.md:28-29` says this
sentence was withdrawn on 2026-08-21. A reader who skips to the summary — which
is what a summary is for — gets the withdrawn claim with no marker, reinforced at
`:513-515` and by the navigation row at `:645`.

**5. `results/CLAIM_INVENTORY.md:337-338` — cross-compound ranking listed under
"What a competent reviewer would accept today".**
Unhedged, in the file `README.md:639` calls *"the single settled answer."* The
same file's rows C2 and C3 were updated on 2026-08-21 to record the
falsification; §2 was not.

**6. `app/server/report.py:2640-2643` — the generated therapy report prints the
ordering its own limit refuses to print.**
*"it moved the top-1 pick outright, **cisplatin → vincristine**."* Limit 2 of the
same document (`:224-228`) declines to name those compounds *"because naming them
would print the very ordering this page has withdrawn."* The therapy report is
the artefact a user **exports and sends onward** — it travels further from this
repository than any other surface audited.

**7. `app/web/js/ui/compound.js:139-141` — the Studio tells a user the shipping
diffusivity model is the one that was falsified.**
*"diffusivity is a function of molecular weight alone, and that has been measured
against 71 compounds and refuted."* The engine ships the two-route law, and the
correct text renders **in the same panel, directly beneath this paragraph**. This
is the bring-your-own-compound path — read by exactly the user about to hand the
model a structure.

*Runners-up, in order:* (8) the **51 pruned runs** whose API payload claims
`timelapse.available: true` while carrying its own deletion record (§2.6) — one
of two defects introduced **today**; (9) **`docs/data_acquisition_plan.md`**,
which documents a module that does not exist and a synthetic pipeline that was
retired, and is **the only stale file in `docs/` with no banner**; (10)
**`docs/real_data/ivygap.md`**, which derives three tables from a cache that is
not on disk and whose loader block raises `FileNotFoundError` on every line;
(11) **`R² = 0.718` asserted flat on twenty-two surfaces** including the therapy
report and the Studio provenance panel, when
`results/hbd_convention_audit.md:3-7` says it *"describes neither the shipping
engine nor the correct physics"* and only the landing page carries the
correction; (12) **`landing-3d/app/layout.tsx:171`** — *"teams comparing
candidate therapies"* in the `<meta description>`; (13) **α = 148.6 on nine
documentation surfaces** against a live `22.223477711172563` (§2.7); (14) the
**four outreach drafts** quoting a *"~2× geometric-mean fold error"* that the
shipped default path does not produce (it gives 3.451×).

---

## 3b. Two operational findings that are not claims

**The prune script's citation protection has a false negative on ranges, and it
already fired.** `scripts/prune_run_payloads.py`'s `RUN_ID_RE = rb"r_\d{5}"`
cannot see a hyphenated range. `docs/RUNBOOK_DEPLOY.md:369` cites
`` `r_00027`-`r_00029` ``. Measured:

| run | `pruned_payloads` | `timelapse.bin` |
|---|---|---|
| `r_00027` | 0 | present |
| **`r_00028`** | **1** | **gone** |
| `r_00029` | 0 | present |

The two endpoints matched the regex and were protected; **the middle of the cited
range did not and was pruned.** The script's own docstring states the rule it
broke: *"a reader following a citation must find what the citation points at."*
Any published record citing a run **range** rather than an explicit list is
unprotected.

**Another process was editing `tmesim/` during this audit.** Transient `# PLANT`
markers were observed at `tmesim/molecular/permeability.py:325`,
`tmesim/reconstruct/barrier.py:227`, `tmesim/treatment/drugs.py:261` and
`app/server/cache.py:358` around 00:59 and had reverted by 01:03. Re-checked at
the end: `grep -rn 'PLANT'` over all four returns **exit 1**, and
`app/server/cache.py:351-362` holds the correct `commonpath`-based `_within`, not
`return True`. A filtered suite run finished **1296 passed, 0 failed**. No
planted defect is live. Recorded because every code-read verdict in this file was
taken during that window.

---

## 4. The shape of the drift

Three patterns, stated because the individual rows do not show them:

**The Studio is ahead of the documents.** `app/server/reach.py`,
`app/web/js/ui/reach.js` and `app/web/js/ui/finding.js` carry the steady-state
cutover, the capacity absence, the 4.2 %-slower correction, the retracted 10×
claim and the 7-of-12 verdict inversion — all correctly, all dated, several with
per-run branches for the case where a record cannot support the sentence.
`README.md` stops at 2026-08-21 and `results/CLAIM_INVENTORY.md` at 2026-08-12.
The engine is not the problem and has not been.

**Corrections propagate to exactly one surface and stop.**
`hbd_convention_audit.md` reached only `landing-3d/lib/copy.ts`.
`verdict_inversion_census.md` reached only `finding.js` and `jobs.py`.
`p90_capacity_sweep.json` reached nothing at all. There is no mechanism that
notices a new `results/` file and asks which existing claims it moves.

**Two files corrected themselves and then re-asserted the uncorrected version
lower down.** `results/kyle2004_penetration.md` corrects 14× at `:504-526` and
re-asserts it in its own conclusion at `:740`; `results/CLAIM_INVENTORY.md`
records the ranking falsification at `:91-92` and endorses ranking at `:337`.
In both cases the surviving downstream copy — the landing page's `14× TOO
SHALLOW`, README's *"the single settled answer"* — took the wrong half.

---

## 5. What this audit did NOT cover

Stated at full extent, because an audit's silence is the thing most likely to be
misread as a pass.

1. **No external measurement was checked against its source.** Every PMID-anchored
   number on the surface — 15 ± 9 mmHg (PMID 1423283), 74–130 µm²/s (PMID 9806526,
   17869669), 12× hindrance, 0.61 ± 0.11 and 0.60 ± 0.20 breast Kp, 4.70 ± 2.38 %
   vascular area, 65 (35–95) /mm² MVD, 194–251 s⁻¹ oxygen supply, the Pruijn 2008
   29.7× spread, the Primeau 40–50 µm anchor — is **UNVERIFIED**. I did not obtain
   a single paper. Where a row above is marked SUPERSEDED it is superseded by a
   *later file in this repository*, never by my reading of a citation.

2. **Nothing was re-run.** `vivum validate` was not executed in either mode, so
   `README.md:466-471`'s tallies (`PASS 7 · FAIL 6 · NOT TESTABLE 4 · DEFERRED 8`
   of 25; `PASS 13 · FAIL 18 · NOT TESTABLE 19` of 50) and
   `README.md:388`'s `PASS 3 FAIL 1 NOT TESTABLE 6` are **UNVERIFIED**. The full
   scorecard is documented at 15–25 minutes and the kill panels considerably more.
   The only live measurement taken was `calibrate_alpha()`.

3. **No simulation output was reproduced.** The 3-D:2-D kill ratios
   (0.2518 / 0.2739 / 0.2081), the grid ladders, the Monte-Carlo rank
   probabilities (0.9277 → 0.3965), the 99.8 % dose-limited share, the 2.0025
   immune spread, the 41× IFP gap — all **UNVERIFIED**. I read them; I did not
   compute them.

4. **`archive/` was not audited.** `archive/synthetic_era/` (28 entries) is a
   dated published record and was outside the enumerated scope. It is reachable
   from `README.md:270` and is therefore claim surface that nobody checked today.

5. **Four root-level contract documents were not audited**:
   `API_CONTRACT.md` (48 KB), `ARCHITECTURE_CONTRACT.md`, `MOLECULAR_CONTRACT.md`,
   `PROTOCOL_LANGUAGE.md`. They were outside the enumerated scope. Two incidental
   hits were recorded in passing (`API_CONTRACT.md:101` and `:745` carry
   `"alpha": 148.6` in example payloads; `API_CONTRACT.md:237-245` documents the
   410 path correctly) but no systematic pass was made. `README.md:621` says
   *"where it and `docs/architecture.md` disagree, the contract wins"* — so an
   unaudited contract outranks an audited document.

6. **`results/*.md` was swept, not read.** All 144 files were swept mechanically
   for prohibited terms, for the specific retracted numbers named in the brief
   (`60x`, `4.4 s`, `~250 s`, `1.74`, `1.56`, `10x grid`), for assertive ranking
   language, and for currency against later dated files. **Roughly 20 were read in
   substance.** A claim that is wrong in a way none of those sweeps would catch —
   a mis-stated method, an inverted inequality, a number transposed inside a table
   — would not appear in this audit. `results/numerical_rigor.md`'s own history
   shows this project has shipped exactly that kind of defect before.

7. **No figure was opened.** `results/*.png` (roughly 40 files, several megabytes
   each) carry axis labels, captions and plotted values that are claim surface.
   None was viewed.

8. **The landing page was read, not rendered.** No `npm run build`, no browser.
   Conditional rendering, the reduced-motion path's divergence from the cinematic
   path, and anything injected at runtime are **UNVERIFIED**. The claim that both
   experiences render `SECTION_COPY` identically was taken from source structure.

9. **The pruned-run defect was reasoned from source, not exercised.** The
   contradiction in §2.6 was established from `app/runs/*/result.json` on disk
   plus `jobs.py:1541-1544` and `jobs.py:1650-1666`. No server was started and no
   `GET /api/runs/{id}` was issued. The *documents* are proven contradictory; the
   served payload is inferred from the code path that serves them.

10. **Most of the `docs/` findings are delegated and second-hand.** All 25
    files were swept, so §2.3 is a *coverage*-complete pass — but I personally
    re-derived only the five most serious rows plus the α and R² classes. Rows
    marked **[D]** rest on another agent's reading. One relayed quote was wrong
    (`DEPLOYMENT.md:743`, corrected in place), which is the base rate to assume
    for the rest: the verdicts are high-confidence, the exact wording of an
    unspot-checked [D] row is not. Anything in §2.3 that would trigger an edit
    should be re-read at its line before acting.

11. **Test coverage of the claim surface was not assessed.** `results/` has no CI
    — that is the premise of this audit — but I did not check whether `tests/`
    pins any of the strings above, which would at least catch a future edit. The
    one relevant data point: `tests/test_app_server.py:1894-1935` exercises
    `pruned_payloads`, including asserting it survives a restart. Note what it
    does **not** assert — that `timelapse.available` is reconciled against it,
    which is the §2.6 defect. A filtered suite run during this audit finished
    **1296 passed, 0 failed**. The suite is green and the claim surface has 187
    violations; those two facts are not in tension, they are the finding.

---

## 6. Reproducing this audit

Every command in this file was run from `/Users/samarthbasanth/tme-sim` with
`.venv/bin/python`. The load-bearing ones:

```bash
# The headline finding.
.venv/bin/python -c "
import json; d=json.load(open('results/p90_capacity_sweep.json'))
a=d['arm_A']['analysis']
print(a['verdict']); print('C1 fail', a['C1']['fail']); print('C2', a['C2'])"

# The live calibration constant, against nine documents that assert 148.6.
.venv/bin/python -c "
from tmesim.validation.calibration import calibrate_alpha; print(repr(calibrate_alpha()))"

# The pruned-run contradiction.
.venv/bin/python -c "
import json,glob
n=b=0
for p in glob.glob('app/runs/*/result.json'):
    d=json.load(open(p))
    if 'timelapse' in (d.get('pruned_payloads') or {}):
        n+=1
        if (d.get('timelapse') or {}).get('available') is True: b+=1
print(n,'pruned;',b,'still report available:true')"

# Which results/ files README actually cites, by date.
for f in $(grep -o 'results/[a-zA-Z0-9_]*\.md' README.md | sort -u); do
  echo "$(git log -1 --format=%ad --date=short -- $f)  $f"; done | sort -r | head

# NOTE: grep here is ugrep and will not word-split $VAR file lists.
# Always pass explicit literal paths.
```
