# The descriptor layer that ships has never been scored

**What this file is.** `results/two_route_promotion.md` publishes **R² = 0.718**
for the two-route transport law against Pruijn 2008, and that number is in the
README, `docs/QUICKSTART.md` and the marketing copy. It was produced by feeding
the law **the paper's own measured `logP₇.₄` and `HD` columns**
(`experiments/two_route_transport.py:246-249`). The shipping engine feeds the
same law RDKit's `Crippen.MolLogP` and `Lipinski.NumHDonors`
(`tmesim/molecular/compile_drug.py:190-197`). Those are not the same inputs, and
the second set has never been scored against anything.

This file measures the gap and **repairs nothing**. No constant was moved, no
descriptor method was changed, `lipophilicity.resolve_logp` was **not** applied.

Driver: `experiments/descriptor_falsification.py`. Code changed:
none in `tmesim/`. Docs changed: `docs/QUICKSTART.md` (a stale line, §6).

---

## 0. Bottom line

| question | answer |
|---|---|
| Can the shipping descriptor layer be run on the Pruijn set at all? | **On 5 of 71 compounds.** The other 66 are ACSRC `SN` code numbers with no public structure. |
| How wrong is the shipping `logP` on the 5 it *can* run on? | **RMS 1.62 log units**, bias +0.51, worst +3.28. |
| Is the `HBD` column even the same quantity? | **No.** Table S1 counts donor *hydrogens*; the engine counts donor *heavy atoms*. Verified 5/5 vs 3/5. |
| Is ionisation handled anywhere on the transport path? | **No.** No pKa, no pH, no logD. Crippen is a neutral-species logP. |
| How much of R² = 0.718 survives an error that size? | **≈ 0.27** (0.266 ± 0.082), or **≈ 0.22** once the donor-count convention is included. |

**The honest headline: very little of the published 0.718 is a claim about the
shipping engine.** It is a claim about the *physics*, evaluated on descriptors
supplied by the measurement. Propagating the descriptor error this file actually
measured through the same scoring run collapses R² from **0.718 to ≈ 0.27** — a
loss of about **62 %** of it.

**But not to zero, and that matters.** The MW-only model this replaced scored
**0.022**. A projected 0.22–0.27 is still an order of magnitude above that, so
the transcellular physics is not erased by bad descriptors — it is *masked* by
them. The two failures are separable, which is the whole point of running this
before the second-cell-line test.

**What this file could NOT do, stated up front.** The clean version of this
experiment — run the engine's descriptors on all 71 compounds — **is not
runnable with any data reachable from this repo.** §1. Everything below is what
is left after that, and it is weaker than what was asked for.

---

## 1. The hard limit: 66 of 71 compounds have no structure

`descriptors.compute_descriptors` takes a SMILES. Pruijn 2008 does not publish
SMILES — Table S1 publishes *drawings*. 66 of the 71 compounds are Auckland
Cancer Society Research Centre internal code numbers (`SN24280`, `SN28287`, …),
unpublished benzotriazine di-N-oxide analogues. PubChem PUG-REST name lookup
returns `PUGREST.NotFound` for them (checked: `SN28287`, `SN29047`, `SN24280`).

They were **deliberately not hand-drawn from the PDF images**. A transcribed
structure is an invented input, and inventing the input is the exact defect this
file exists to detect.

Five compounds are public and were resolved from PubChem, each verified against
Table S1's own MW column (`--stage structures` **asserts** agreement to 1 Da and
fails the run otherwise):

| compound | Table S1 MW | RDKit MW | \|Δ\| | PubChem CID |
|---|---:|---:|---:|---:|
| Urea | 60.1 | 60.06 | 0.04 | 1176 |
| Mannitol | 182.0 | 182.17 | 0.17 | 6251 |
| SN24280 (TPZ) | 178.0 | 178.15 | 0.15 | 36400 |
| acridan-9(10H)-one | 195.0 | 195.22 | 0.22 | 2015 |
| NLCQ-1 | 332.0 | 331.76 | 0.24 | 398939 |

The TPZ SMILES is byte-identical to the repo's own `TPZ_SMILES`, asserted, so
the anchor column is the one the project already publishes. NLCQ-1's PubChem
record is the **hydrochloride** (MW 368.2); the `.Cl` was dropped because Table
S1's MW is the free base, and the MW assertion is what establishes that this was
the right call rather than a convenience.

> **Consequence.** No amount of re-scoring can turn 0.718 into a score of the
> shipping engine, because on 93 % of this dataset the shipping engine produces
> no answer at all. That is a property of the dataset, not of effort.

---

## 2. Per-compound descriptor error — the deliverable

`--stage descriptors`. Engine `logP` = RDKit Crippen `MolLogP`, the **neutral**
octanol/water partition coefficient. Paper `logP₇.₄` = a **measured distribution
coefficient at pH 7.4**.

| compound | Table S1 `logP₇.₄` | Crippen `logP` | **Δ logP** | S1 `HD` | engine `HBD` | Δ HBD |
|---|---:|---:|---:|---:|---:|---:|
| Mannitol | −2.62 | −3.585 | **−0.965** | 6 | 6 | 0 |
| NLCQ-1 | +0.22 | **+3.495** | **+3.275** | 1 | 1 | 0 |
| SN24280 (TPZ) | −0.38 | −0.916 | −0.536 | 2 | **1** | **−1** |
| Urea | −2.01 | −0.976 | **+1.034** | 4 | **2** | **−2** |
| acridan-9(10H)-one | +2.96 | +2.681 | −0.279 | 1 | 1 | 0 |

**mean +0.506, RMS 1.618, MAE 1.218, range [−0.965, +3.275].**

### 2.1 The failures are a class, not a scatter

Grouped by structure the five separate cleanly, and the boundaries are chemical:

* **Neutral aromatics** (acridone, TPZ): |Δ| ≤ 0.54. Crippen's home ground. Both
  slightly *too polar* — the harmless direction.
* **Small neutral polar** (urea, mannitol): |Δ| ≈ 1.0, **both signs**.
  Atom-contribution failure on hydrogen-bonding polyols and amides. This is the
  same class `lipophilicity.py` already documents for methotrexate (+2.1) and
  cisplatin (≈ +2.4), so it is a *third and fourth* instance of a known defect,
  not a new one.
* **Ionisable base** (NLCQ-1): **Δ = +3.28**, the largest error in the set by
  more than 3×. The 4-aminoquinoline is protonated at pH 7.4, so its `logD₇.₄`
  sits far below its neutral `logP`. **One-signed and always positive.**

**The actionable form: the engine is safe on neutral aromatics and dangerous on
bases.** Any compound with a basic nitrogen — which is most oncology
small-molecule chemistry — inherits a one-directional, multi-log-unit
over-estimate of lipophilicity, which opens the transcellular route too far and
**over-predicts D**.

### 2.2 Ionisation is not handled. At all.

Asked of the code, not of the docs (`_ionisation_audit`, which reads the shipping
call path with `inspect.getsource`):

```
descriptors.compute_descriptors -> 'logp' is Crippen.MolLogP:        True
compile_drug passes descriptors['logp'] straight to two_route:       True
pH / pKa / logD / ionisation tokens anywhere on that path:           NONE
```

`logP` and `logP₇.₄` are **different quantities**, and the engine computes the
first while every published transport number was scored against the second. This
is a **systematic discrepancy with a fixed sign**, not noise: for any compound
ionised at pH 7.4, `logD₇.₄ < logP`, so the engine reports it more lipophilic
than it is.

The repo does carry a pKa table — `tmesim/molecular/ion_trapping.py` — but it
serves lysosomal trapping and **is not on this path**. The capability exists in
the codebase and the transport model does not use it.

### 2.3 The `HD` column is not the `hbd` the engine computes

This is the finding I did not expect. RDKit exposes two donor counts and Table S1
tracks the one the engine **does not** use:

| compound | Table S1 `HD` | `Lipinski.NumHDonors` (engine) | `Descriptors.NHOHCount` |
|---|---:|---:|---:|
| Mannitol | 6 | 6 | 6 |
| NLCQ-1 | 1 | 1 | 1 |
| SN24280 (TPZ) | 2 | **1** | 2 |
| Urea | 4 | **2** | 4 |
| acridan-9(10H)-one | 1 | 1 | 1 |
| **agreement with Table S1** | — | **3 / 5** | **5 / 5** |

`NumHDonors` counts donor **heavy atoms** (an `–NH₂` is 1). `NHOHCount` counts
donor **hydrogens** (an `–NH₂` is 2). On this evidence Pruijn's `HD` is the
hydrogen count. So the published R² was scored on one quantity and the engine
consumes a different one — and for every `–NH₂` or `–NH₃⁺` the engine's count is
**one lower**.

That is not cosmetic. `two_route` multiplies this count by
`HBD_LOGP_DECREMENT = 1.34` log units per donor. **One unit of miscount is 1.34
log units of membrane partition coefficient — a factor of 22 in `K_mem`.**

> This bears on the `1.34` question and I am not going to pretend otherwise, so
> here is the careful version. `1.34` is sourced from Abraham LSER coefficients
> and is *not* the value that maximises the Pruijn score (1.0 is, at 0.725 vs
> 0.718), and `tests/test_two_route.py:558` exists to stop someone "improving"
> it. **Nothing here argues for touching 1.34.** What it argues is that the
> *count* 1.34 multiplies may be the wrong quantity, and that is a separate
> question from the coefficient's value — one that can be settled by reading
> Pruijn's Table S1 convention rather than by fitting anything. Deliberately left
> for a human. §7.

### 2.4 End-to-end, on the five it can run on

Shipping compile (`build_profile`, two-route ON, engine descriptors, engine MW)
against the measured HT29 `D_mcl`, alongside the same compound's row in the
published Pruijn run:

| compound | measured | **shipping** | ratio | published run (paper descriptors) | ratio |
|---|---:|---:|---:|---:|---:|
| Mannitol | 21.1 | 43.1 | 2.04× | 43.1 | 2.04× |
| **NLCQ-1** | 60.6 | **200.6** | **3.31×** | 82.2 | **1.36×** |
| SN24280 (TPZ) | 40.4 | 49.3 | 1.22× | 44.4 | 1.10× |
| Urea | 45.0 | 65.5 | 1.46× | 65.1 | 1.45× |
| acridan-9(10H)-one | 317.0 | 239.5 | 0.76× | 240.4 | 0.76× |

The one ionisable compound is the one that breaks: **1.36× → 3.31×** on the
descriptor substitution alone. Four of five barely move, which is consistent with
§2.1 — the engine is fine until it meets a base.

n = 5 is far too few for an R², which is why ratios are reported instead.

---

## 3. The large-n arm came back NULL, and it is reported as null

`--stage substitute`. The plan was: 49 compounds carry **both** a measured and a
"calculated" `logP₇.₄` in Table S1, so score the identical physics twice and
change only the descriptor source.

It does not work, for a reason worth recording. From the paper's Experimental
section: the calculated column is **ACD LogP/LogD prediction software v8.0**,
and the authors used *"ACD/LogP System Training and Accuracy Extender"* to train
it **on 98 compounds of this very series** before predicting the rest. The
consequence is measurable in the data:

* **22 of 49** calculated values are **exactly equal** to the measured value.
* **26 more** differ by ≤ 0.15 log units.
* **Exactly one** (TPZ, +0.71) differs by more than that.
* RMS "error": **0.107 log units** — an order of magnitude better than any
  general-purpose logD predictor, because it is not a prediction.

Scored: R² **0.6756 → 0.6744**, a change of **−0.0012** on n = 49.

**That number must not be quoted as evidence that computed descriptors are
harmless.** It is evidence that this column is not a computed descriptor. The
large-n version of the question **cannot be answered with data in hand.**

---

## 4. Propagation: how much of 0.718 survives

`--stage propagate`. Baseline reproduced live: **R² = 0.718** on n = 68,
asserted to match the published figure to 3.1e-4 before anything else runs.
Physics and every constant held fixed; only the descriptors move.

### 4.1 Systematic `logP` offset — the shape ionisation error takes

| Δ logP | R² | Spearman | spread | geo bias | within 3× |
|---:|---:|---:|---:|---:|---:|
| −2.00 | 0.322 | +0.494 | 6.39× | 0.58× | 72 % |
| −1.00 | 0.588 | +0.751 | 8.10× | 0.72× | 94 % |
| **0.00** | **0.718** | **+0.821** | **8.34×** | **1.03×** | **99 %** |
| +1.00 | 0.673 | +0.852 | 7.24× | 1.61× | 91 % |
| +2.00 | 0.455 | +0.770 | 6.57× | 2.24× | 66 % |

A uniform offset is the *gentlest* form of descriptor error — it is largely
absorbed into the geometric bias, and Spearman actually **improves** at +1.0.
R² is remarkably tolerant of it. Note the implication: **a single global logP
shift cannot be diagnosed from R² alone**; it shows up in the bias column.

### 4.2 Random `logP` scatter — the shape estimator error takes

200 replicates per row, seed 20260821.

| σ (log units) | mean R² | sd | p05 | p95 | mean ρ |
|---:|---:|---:|---:|---:|---:|
| 0.00 | 0.718 | — | — | — | +0.821 |
| 0.25 | 0.692 | 0.022 | 0.656 | 0.731 | +0.806 |
| 0.50 | 0.627 | 0.045 | 0.563 | 0.696 | +0.772 |
| 0.75 | 0.541 | 0.064 | 0.447 | 0.649 | +0.720 |
| **1.00** | **0.447** | 0.071 | 0.320 | 0.551 | +0.653 |
| 1.50 | 0.303 | 0.083 | 0.166 | 0.430 | +0.538 |
| 2.00 | 0.229 | 0.085 | 0.099 | 0.382 | +0.461 |

**Unbiased scatter is far more destructive than systematic bias.** σ = 1.0 log
unit — a *modest* expectation for a generic cLogP estimator — already halves R².
This is the single most useful row in the file for the sibling experiment: a
second cell line that scores ≈ 0.45 is fully consistent with correct physics and
σ ≈ 1 descriptors, and says nothing about the physics on its own.

### 4.3 The donor-count convention

Which compounds carry an `–NH₂` is unknowable without the 66 missing structures,
so this is a **bracket, not a claim**:

| HD shift | R² | Spearman | spread | geo bias | within 3× |
|---:|---:|---:|---:|---:|---:|
| **0 (published)** | **0.718** | +0.821 | 8.34× | 1.03× | 99 % |
| −1 for every compound | **0.569** | +0.836 | 6.55× | 1.77× | 78 % |
| −2 for every compound | 0.362 | +0.718 | 6.39× | 2.16× | 66 % |

The true value lies between 0.569 and 0.718 depending on how many of the 66
carry a primary amine. TPZ — the parent of 67 of the 71 — **does**, so it is not
a rare substructure in this series.

### 4.4 Combined, at the magnitudes this file measured

Bias and σ are taken from §2, not chosen. 200 replicates.

| scenario | bias | σ | ΔHD | mean R² | sd |
|---|---:|---:|---:|---:|---:|
| Crippen error as measured (n = 5) | +0.51 | 1.72 | 0 | **0.266** | 0.082 |
| … + engine donor convention | +0.51 | 1.72 | −1 | **0.215** | 0.085 |
| Crippen error excluding the base (n = 4) | −0.19 | 0.86 | 0 | **0.495** | 0.067 |
| … + engine donor convention | −0.19 | 0.86 | −1 | **0.402** | 0.072 |

**Read the range, not the point estimate: 0.22 – 0.50, against a published
0.718.** The spread between the top and bottom rows is entirely "does the
compound set contain ionisable bases or not", which is the §2.1 result restated.

These rows rest on a **five-compound** estimate. That is the weakest number in
this file. It is used because it is the only direct measurement of the shipping
descriptor against this dataset that exists — and its weakness *is* the finding.

---

## 5. What this means for the second-cell-line test

The sibling experiment predicts a second cell line with unchanged constants. The
de-confounding rule that falls out of §4:

| observed on the second line | consistent with |
|---|---|
| R² ≳ 0.6 | physics **and** descriptors both fine |
| R² ≈ 0.45 – 0.60 | descriptor scatter σ ≈ 0.5–1.0 alone. **Physics unindicted.** |
| R² ≈ 0.20 – 0.45 | descriptor error at the magnitude measured here, alone. **Physics still unindicted.** |
| R² < 0.20 | below anything descriptor error alone produces here — the first regime that actually implicates the physics |
| Spearman holds up while R² collapses | systematic bias (ionisation), not scatter — check the geometric-bias column |
| geometric bias moves, R² does not | uniform logP offset (§4.1) |

**Provided the second-line test uses the paper's descriptor columns, as this one
did, none of the above applies to it** — it is then a clean test of the physics
and of the cell-line-transferability of the constants. **If it uses engine
descriptors, it is confounded**, and the table above is the decoder.

The recommendation is therefore concrete: **run the second-line test on the
measured descriptor columns**, exactly as the first was, and treat the
descriptor question as this separate, unresolved thing.

---

## 6. `docs/QUICKSTART.md` corrected

Lines 214–215 said the two-route model was *"default-OFF and wired into
nothing"*. That has been false since `fd2a9e7`; `compile_drug.py:90` reads
`_TWO_ROUTE_DEFAULT = True` and `build_profile` applies the term on every
compile (verified in code, not from the changelog). The line also **contradicted
line 199 of the same file**, eleven lines above, which describes the same model
as the source of the shipping `D_eff` column.

Corrected, and the descriptor caveat from this file added alongside it, because a
reader who has just been told `D_eff` scores R² = 0.718 should be told in the
same breath that the number was earned on inputs the engine does not compute.

---

## 7. What I did not do, and what is left for a human

**Not done, deliberately:**

1. **`lipophilicity.resolve_logp` was not applied.** It exists, has no callers,
   and applying it once moved `alpha` from 148.6 to 68.6 and broke two bitwise
   pins. Nothing here was made to depend on it.
2. **`HBD_LOGP_DECREMENT` was not touched.** §2.3 is about the *count*, not the
   coefficient.
3. **No structure was transcribed from the Table S1 drawings.**
4. **No existing results file was edited.** Two inconsistencies found in passing
   are reported rather than fixed (below).

**Decisions that are now a human's:**

* **Does §2.3 justify applying `lipophilicity.resolve_logp`?** My read: *it is
  orthogonal.* `resolve_logp` fixes measured-vs-Crippen for **11 library drugs**
  and does nothing for an arbitrary user compound, which is where §2.1's failure
  class actually bites. It would improve the library drugs and leave the
  structural defect — no pKa anywhere — untouched. **The thing §2 argues for is
  a logD₇.₄ estimate, not a better logP table.** But it is a real improvement on
  a real defect and the case for it is stronger after this file than before it.
* **Should the engine's `hbd` be `NHOHCount` instead of `NumHDonors`?** §2.3 says
  the published score was earned on the hydrogen count. This is a
  `tmesim/molecular/**` change, outside my territory, and it moves `alpha`. It is
  also the cheapest single thing in this file.
* **Should the published 0.718 be re-labelled?** It is not wrong, but everywhere
  it appears it should say *"scored on measured descriptors"*. §6 does this for
  QUICKSTART only.

**Two inconsistencies found in passing, not fixed (not my territory):**

* `results/two_route_promotion.md` §4 gives tirapazamine *"compiled from its own
  SMILES"* as **44.4 µm²/s (1.10×)**. 44.4 is the value from the Pruijn row,
  which uses the **paper's** descriptors. The actual SMILES compile is
  **49.3 µm²/s (1.22×)** — the same document quotes `49.3322` eleven lines
  above, in the pre-flight block. Both numbers reproduced here (§2.4).
* `data_cache/pruijn2008/pruijn2008_tableS1.csv` spells one compound
  `Acridan-9(10H)-one` where `experiments/diffusivity_falsification.py` spells it
  `acridan-9(10H)-one`. Any case-sensitive join silently drops it (this file
  matches case-insensitively and says so). Owned by a sibling.

---

## 8. Reproducing this

```bash
./.venv/bin/python experiments/descriptor_falsification.py --stage structures
./.venv/bin/python experiments/descriptor_falsification.py --stage descriptors
./.venv/bin/python experiments/descriptor_falsification.py --stage substitute
./.venv/bin/python experiments/descriptor_falsification.py --stage propagate
```

Records land in `results/sweep/descriptor_run_<stage>.json`. `--stage substitute`
**asserts** that the live baseline still reproduces R² = 0.718 before it reports
anything, so this file cannot quietly rot against a moved engine. Monte-Carlo
rows are seeded (`20260821`) and reproduce exactly.
