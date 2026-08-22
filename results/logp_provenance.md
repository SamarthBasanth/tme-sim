# logP provenance audit and blast-radius test

> # ⚠️ CORRECTED 2026-08-14 — §1 FINDING 2 AND §2.3 ARE FALSE OF THE SHIPPING ENGINE. logP NOW ENTERS TRANSPORT DIRECTLY.
>
> **This banner is the only edit made to this file.** No number, table or
> recommendation below it has been changed, re-run or corrected. Everything below
> was true of the engine as it stood on **2026-08-08**, the date of the audit.
>
> **What changed.** Commit `fd2a9e7` (2026-08-14, *"Flip two-route transport ON:
> the first physics change that adds capability"*) made the two-route cell term
> the shipping default. `tmesim/molecular/compile_drug.py::build_profile` now
> passes `descriptors["logp"]` and `descriptors["hbd"]` into
> `two_route.two_route_diffusivity` and multiplies its hindrance factor into
> `effective_D` on **every compile**. The dated record is
> `results/two_route_promotion.md`.
>
> **Which specific statements are now false.**
>
> | where | statement | status |
> |---|---|---|
> | §1, finding 2 | *"logP cannot move transport directly… logP is not an input to it at any point."* | **FALSE as of `fd2a9e7`.** `logp` is an input to `effective_diffusivity`. |
> | §1, finding 2 | *"Its only routes to an output are the first-order loss rate (`binding`, `clearance`) and the delivery scalar (`bbb_penetrance`)."* | **INCOMPLETE.** There is now a fourth route, and it is the transport coefficient itself. |
> | §2.3 | *"`No logP anywhere in that chain.`"* and the five-line `effective_D` chain above it | **SUPERSEDED.** The chain gained a terminal `× two-route cell term` step that consumes `logp` and `hbd`. |
> | §2.3 | *"doxorubicin and epirubicin… both compile to `D = 184.75 µm²/s`"* | **The equality survives, the value does not.** Crippen is stereochemistry-blind, so dox and epi still receive identical `logp`/`hbd` and therefore still compile to the same `D`. That `D` is no longer 184.75 µm²/s. |
> | §2.4 | *"Because logP touches no expensive quantity, the sweep does not need to re-run the Brownian simulation per scenario."* | **NO LONGER VALID as a method.** `experiments/logp_sensitivity.py` caches `d_eff` per drug and re-evaluates only the closed-form terms; under the shipping engine `d_eff` is itself a function of the swept variable. Every ratio in §5 is therefore a **lower bound on the blast radius, not a bound on it**. |
>
> **What this does to §5.3, the file's dominant finding.** §5.3 measures that
> correcting doxorubicin's logP from Crippen `0.001` to the measured `1.27` moves
> `alpha` from `148.6082` to `68.6368` (0.4619×) and rescales every other drug's
> penetration length by ~1.47×. **That measurement was made on the MW-only
> engine and does not transfer.** Under the shipping engine the same edit moves
> the anchor's own `D_eff` as well as its loss rate, so both sides of the
> Primeau `D/k` fit move together and the ~1.47× global rescale cannot be assumed
> to hold in magnitude or in sign. `alpha` has separately moved to
> **22.223477711172563** at the promotion (`results/two_route_promotion.md`),
> so `148.6082` is no longer the starting point either. **§5.3 must be re-measured
> before it is acted on.** Nothing in §5.3 has been deleted; it is the record of
> what the pre-`fd2a9e7` engine did.
>
> **What is NOT affected, and still stands.** §3 (the literature audit: which
> logP values are measured, refuted, laundered or propagated) is a claim about
> the *literature*, not about the engine, and is untouched by this. §6's
> per-drug recommendations on *provenance* — adopt cisplatin's `-2.53`, store
> palbociclib as pH-resolved logD, do not adopt tamoxifen's 7.1 or epirubicin's
> 1.27 — stand unchanged. What does not stand is the sequencing argument in §6.2,
> which rests on §5's magnitudes.
>
> **Why this matters operationally.** `tmesim/molecular/lipophilicity.py::resolve_logp`
> exists with measured values and, as of this banner, **still has no caller** in
> `tmesim/` or `app/` — the only references outside the module itself are in
> `experiments/descriptor_falsification.py`. Anyone reading this file's finding 2
> would conclude that wiring it up is transport-inert. **It is not.** Wiring
> `resolve_logp` in is now a change to the transport coefficient of every
> compound whose logP it supersedes, and needs a revalidation run, not a
> provenance review.

**Scope.** Where octanol/water logP is used in tmesim, what the repo's current values
actually are, how they compare to a literature audit of nine drugs, and how far a
simulated output moves if a logP is wrong.

**Files touched by this audit:** `results/logp_provenance.md` (this file) and
`experiments/logp_sensitivity.py` (new). Nothing under `tmesim/` was modified.
All proposed code changes are in section 6 as a diff to be applied by someone else.

---

## 1. Bottom line

Four findings, in order of how much they should change what anyone does next.

1. **`DRUG_LIBRARY` in `tmesim/core/config.py` has no logP field.** There is no
   literature logP anywhere in the repo on `main`. **Every logP tmesim uses is an
   RDKit Crippen `MolLogP` prediction computed from SMILES at runtime.** The audit's
   nine-drug table therefore does not describe values this repo holds — it describes
   values this repo *never adopted*. The audit's headline hazards (epirubicin's
   propagated 1.27, palbociclib's refuted 0.99, tamoxifen's laundered 7.1) are **not
   present in `main`**. That is the good news.

2. **logP cannot move transport directly.** `DrugParams.diffusivity` is a pure
   function of hydrodynamic radius, which is a pure function of MW. logP is not an
   input to it at any point. Its *only* routes to an output are the first-order loss
   rate (`binding`, `clearance`) and the delivery scalar (`bbb_penetrance`). This
   bounds the blast radius a priori and is why the sweep in section 5 was cheap —
   each drug is compiled once and every logP scenario is then closed-form.

3. **But the two surviving routes are not small, and they multiply.** Loss rate sets
   the penetration length `L`; `bbb_penetrance` is a *linear* multiplier on the
   delivery source. Drug delivered per unit vessel length scales as `bbb x L`. Over
   the plausible error ranges this composite moves by up to **6.6x** (tamoxifen) and
   **2.2x** (palbociclib's L alone). This is not a null result.

4. **The single most consequential logP in the whole repo is doxorubicin's, and it is
   not one of the nine audited drugs.** Doxorubicin is the sole anchor for the one
   free calibration constant `alpha`. Correcting doxorubicin's logP from the Crippen
   0.001 to the measured 1.27 changes `alpha` by 0.462x and **rescales every other
   drug's penetration length by ~1.47x**, globally, without any of those drugs' own
   logP changing. A per-drug sweep that did not hold this separate would misattribute
   a global recalibration as nine local effects.

---

## 2. The call graph

Established by grep + direct read of the modules, and confirmed numerically in
section 4 by reproducing the frozen `DRUG_LIBRARY` rows from Crippen logP alone.

### 2.1 Where logP is born

There is exactly one producer:

```
tmesim/molecular/descriptors.py:85
    "logp": round(Crippen.MolLogP(mol), 3),
```

reached via `descriptors.compute_descriptors(smiles)` <- `descriptors.from_smiles(name, smiles)`.
Both drug libraries call it:

- `tmesim/molecular/library.py`         -> `GBM_DRUGS` (SMILES) -> `from_smiles`
- `tmesim/molecular/breast_library.py`  -> `BREAST_DRUGS` (SMILES) -> `from_smiles`

The only other source is `descriptors._FALLBACK_TABLE` (lines 51-65), used **only**
when RDKit cannot parse the SMILES. Its comment states outright that its values were
"pre-computed from RDKit" — so the fallback is also a Crippen prediction, not
literature. On `main` the fallback fires for exactly one entry, `bevacizumab`
(`logp: -1.0`, an antibody with no small-molecule graph). Every real small molecule
resolves with `source == "rdkit"` (verified: all 10 GBM + all 8 breast drugs).

### 2.2 Where logP goes

```
descriptors["logp"]
│
├─ permeability.membrane_permeability(descriptors)      [permeability.py:198-200]
│     Kp = 10 ** clamp(logp, -3, +3)                    (_LOGP_KP_MAX = 3.0)
│     -> MolecularProfile.membrane_permeability
│     -> DEAD END w.r.t. the engine. Not a DrugParams field, never reaches the PDE.
│        Consumed only by reporting/validation:
│          validation/glioma_regional_pk.py, validation/gist_imatinib.py,
│          validation/moncayo_explant.py
│
├─ permeability.bbb_penetrance(descriptors)             [permeability.py:295-324]
│     d_logp = exp(-(logp - 2.5)^2 / (2 * 1.7^2))       Gaussian, centre 2.5, width 1.7
│     geometric-mean desirability with d_mw, d_tpsa, d_hbd
│     -> DrugParams.bbb_penetrance                      [compile_drug.py:289]
│     -> LIVE. treatment/drugs.py:152
│          S(x,y) = vascular_density(x,y) * bbb_penetrance * plasma_conc
│        i.e. a LINEAR multiplier on the delivery source term of the transport PDE.
│        Also enters reconstruct/barrier.py:405 (effective_penetrance).
│
├─ binding.plasma_free_fraction(descriptors)            [binding.py:88-91]
│     ff = 1 / (1 + 10**((logp - 3.0)/1.5))
│     -> MolecularProfile.free_fraction
│     -> DEAD END w.r.t. the engine. Not a DrugParams field. Reporting only
│        (validation/gist_imatinib.py:400).
│
├─ binding.tissue_binding_rate(descriptors)             [binding.py:116-120]
│     score = (logp - 2.5)/2.5 + (hbd - 2)/4 ; rate = geom_interp(1e-5, 5e-3, sigmoid(score))
│     -> DrugParams.binding                             [compile_drug.py:287]
│     -> LIVE, and AMPLIFIED. validation/calibration.py:210 uses this exact field as
│        `uptake_proxy`, then calibrated_drugparams sets
│          binding := alpha * uptake_proxy(drug)
│        with alpha ~= 148.6 fit on doxorubicin.
│
└─ binding.clearance_rate(descriptors)                  [binding.py:142-146]
      score = 0.4*(logp - 2)/3 - 0.4*(mw - 300)/700 ; rate = geom_interp(1e-5, 1e-3, sigmoid(score))
      -> DrugParams.clearance                           [compile_drug.py:288]
      -> LIVE, un-amplified (calibration.py:226 k_passive).
```

### 2.3 What logP does NOT feed

`DrugParams.diffusivity` — the transport coefficient — is set at
`compile_drug.py:284` from `profile.effective_diffusivity`, which is built at
`compile_drug.py:119-139` from:

```
radius_hydrodynamic_nm  <- descriptors.hydrodynamic_radius_nm(mw, is_antibody)   [MW only]
free_D                  <- stokes_einstein.free_diffusivity(radius, T)           [radius only]
micro_D                 <- brownian.simulate_effective_diffusivity(free_D, phi)  [free_D, radius]
mesh_h                  <- ecm_mesh.ogston_renkin_hindrance(radius)              [radius only]
effective_D             = micro_D * mesh_h
```

**No logP anywhere in that chain.** Confirmed numerically: doxorubicin and
epirubicin, which have bit-identical descriptors, both compile to `D = 184.75 um^2/s`,
and every logP scenario in section 5 leaves `D` untouched by construction.

Also unaffected by logP: `ic50` (from `molecular/potency.py` / GDSC2 measurements),
`hill` (constant 1.2), `molecular_weight`.

### 2.4 The consequence for the sensitivity test

Because logP touches no expensive quantity, the sweep does not need to re-run the
Brownian simulation per scenario. `experiments/logp_sensitivity.py` compiles each
drug **once**, caches `d_eff`, and re-evaluates `tissue_binding_rate`,
`clearance_rate` and `bbb_penetrance` in closed form for every logP value. The full
12-drug, 40-scenario panel runs in the time of 12 compiles.

---

## 3. Per-drug audit

### 3.1 What the repo actually holds

There is no recorded logP provenance in the repo for any drug, because there is no
stored logP. The provenance for all nine is identical and is structural:
**RDKit Crippen `MolLogP`, computed at runtime from a PubChem SMILES.** The SMILES
themselves are verified (`breast_library.py:112-124` re-checks every one against the
PubChem molecular weight to <0.1 g/mol); it is the *logP derived from them* that is a
prediction.

| Drug | Repo logP (Crippen) | Where the SMILES lives | Repo provenance | Audit verdict on the *literature* value | Repo exposure |
|---|---|---|---|---|---|
| Cisplatin | **+0.195** | `library.py` GBM_DRUGS | RDKit Crippen prediction | **FOUND: -2.53 +/- 0.28** (Screnci 2000, PMID 10732773, Table 3; shake-flask n=5) — but octanol/**0.9% saline**, total Pt by ICP-MS | **Worst absolute error in the panel: +2.73 log units.** Crippen is a topological additive model and has no parameters for a square-planar Pt(II) complex; it should not be trusted here at all |
| Paclitaxel | **+3.736** | `breast_library.py` | RDKit Crippen prediction | NOT FOUND | No measured value to compare against; unverifiable in either direction |
| Docetaxel | **+3.260** | `breast_library.py` | RDKit Crippen prediction | NOT FOUND. 4.26 downgraded mid-audit (two EMA EPARs = same ASMF boilerplate quoted twice, no method); Sanofi Taxotere SDS conflicts at 3.2 | Crippen 3.260 happens to sit between the two unsourced numbers. Coincidence, not corroboration |
| Epirubicin | **+0.001** | `breast_library.py` | RDKit Crippen prediction | NOT FOUND — 1.27 is **probable propagation** from doxorubicin's HSDB/Hansch field (epirubicin is the 4'-epimer) | **The propagation never happened in `main`.** But see 3.2 — a *new* artifact is about to be created |
| Vinorelbine | **+4.754** | `breast_library.py` | RDKit Crippen prediction | NOT FOUND — logD only | Unverifiable. The real OECD 107 logD exists but its pH-resolved values were not retrieved, so nothing can be substituted without new retrieval |
| Tamoxifen | **+5.996** | `breast_library.py` | RDKit Crippen prediction | NOT FOUND — 7.1 is a **laundered prediction** (verified against the OPERA/PHYSPROP files); 7.87 +/- 0.75 (El-Kattan 2001, PMID 11250108) is a **calculator output**, not a measurement | **Prediction vs prediction.** 5.996 is not "wrong relative to 7.1"; both are model outputs. The honest statement is that tamoxifen's logP is *unknown*, with a ~2-log-unit spread across predictors. Open lead: Takacs-Novak & Avdeef, PMID 8877846 — genuine interlab shake-flask study of 23 drugs, paywalled, **could not determine whether tamoxifen is among them.** Open lead, not a negative |
| Irinotecan | **+4.091** | `library.py` GBM_DRUGS | RDKit Crippen prediction | NOT FOUND — logD only | Unverifiable; same situation as vinorelbine |
| Palbociclib | **+2.966** | `breast_library.py` (structure reused by `library.py`) | RDKit Crippen prediction | **NOT FOUND — 0.99 actively REFUTED** (absent from Health Canada IBRANCE monograph rev. 2024-09-04 s.13 and FDA label 207103s008; DrugBank's citation does not check out). Measured logD **does** exist: EMA Ibrance EPAR EMA/652627/2016 p.28, OECD 107, pH5 0.228 / pH7 1.11 / pH9 2.26 | **The refuted 0.99 is not in the repo — but 2.966 is *above the entire measured logD range*, including the pH 9 endpoint.** Highest-confidence correctable error in the panel |
| Ribociclib | **+2.796** | `library.py` GBM_DRUGS | RDKit Crippen prediction | NOT FOUND — logD only | Unverifiable; same situation as vinorelbine |

Additional drugs in the same code path, included because one of them is the
calibration anchor:

| Drug | Repo logP (Crippen) | Note |
|---|---|---|
| **Doxorubicin** | **+0.001** | **THE ANCHOR.** Holds the canonical HSDB field `log Kow = 1.27` citing Hansch p.187 (CID 31703). Crippen is off by 1.27 log units on the one drug that sets `alpha` for everything |
| Temozolomide | -2.078 | |
| Cyclophosphamide | +1.884 | |
| 5-fluorouracil | -0.798 | |
| Methotrexate | +0.268 | |

### 3.2 Two artifacts to flag explicitly

**(a) The epirubicin/doxorubicin stereoisomer trap.** Crippen `MolLogP` is a
topological atom-contribution model; it is blind to stereochemistry. Doxorubicin and
epirubicin differ only in the 4'-OH configuration of the sugar, so RDKit returns
**bit-identical descriptor dicts** for them — verified:

```
mw 543.525 | logp 0.001 | tpsa 206.07 | hbd 6 | hba 12 | rot 5 | arom 2 | charge 0 | r_h 0.5386
```

for both. `breast_library.py:28-32` markets epirubicin as "a real same-class,
same-tissue contrast to the panel's doxorubicin anchor." **On `main` that contrast
does not exist at the descriptor level** — every logP-derived quantity is identical
for the two drugs, and the entire epirubicin-vs-doxorubicin difference in any panel
output comes from IC50 alone.

This matters *more*, not less, once a measured-logP table is introduced. A table that
contains doxorubicin's real 1.27 but not epirubicin's (correct, since epirubicin has
no Hansch-backed field of its own) makes two stereoisomers with identical computed
descriptors receive logP **1.27 vs 0.001** — a 1.27-log-unit split created purely by
table membership, not by chemistry. Section 5 shows that split alone moves
epirubicin's `L` by 0.68x. Whatever is done, dox and epi must be handled as a pair
with a written justification, not independently.

**(b) A measured-logP table already exists on another branch, and it contains the
anchor.** The in-progress worktree at
`.claude/worktrees/agent-ac4e848d0655b5856/tmesim/molecular/lipophilicity.py`
defines `EXPERIMENTAL_LOGP` with 10 entries: `temozolomide, lomustine, carmustine,
procarbazine, methotrexate, doxorubicin, etoposide, vincristine, 5-fluorouracil,
cyclophosphamide`. Two observations, both load-bearing:

- **None of the nine audited drugs is in it.** Palbociclib is absent (confirmed by
  key enumeration), which independently corroborates that the branch's palbociclib
  logP of +2.966 is the RDKit Crippen value carried through unchanged — it is
  numerically identical to `main`'s. The refuted 0.99 was correctly never adopted.
- **Doxorubicin IS in it, at 1.27.** That is the anchor. Merging that table therefore
  triggers the global `alpha` recalibration quantified in section 5.3 — a ~1.47x
  rescale of every drug's penetration length — as a *side effect* of a change that
  reads like a per-drug provenance fix. This should be surfaced in that branch's own
  writeup before it merges.

---

## 4. Verification that the audited quantities are the live ones

`experiments/logp_sensitivity.py` reproduces the frozen `DRUG_LIBRARY` snapshots from
Crippen logP alone, which proves the call graph in section 2 is the real one and not
a plausible-looking read:

```
palbociclib-like   alpha*binding 0.044354404561 vs 0.044354404561 | clearance 1.0525e-04 vs 1.0525e-04 | bbb 0.578756 vs 0.578756  -> MATCH
ribociclib-like    alpha*binding 0.039932267072 vs 0.039932267072 | clearance 1.0342e-04 vs 1.0342e-04 | bbb 0.799518 vs 0.799518  -> MATCH
irinotecan-like    alpha*binding 0.060124566987 vs 0.060120000000 | clearance 1.1414e-04 vs 1.1410e-04 | bbb 0.150914 vs 0.151000  -> MATCH
```

(palbociclib/ribociclib match to 12 decimal places; irinotecan's row was written to 4
significant figures and matches to that precision.) The recomputed `alpha` is
`148.608158`, identical to the `148.60815765445545` frozen in `config.py`'s comment.

So: `DRUG_LIBRARY["palbociclib-like"].binding = 0.0443544` **is**
`148.608 x tissue_binding_rate(logp=2.966, hbd=2)`. The Crippen prediction is baked
into the committed library values.

One exception worth noting: the three Stage-0 rows (`temozolomide-like`,
`doxorubicin-like`, `cisplatin-like`) carry hand-set numbers
(e.g. `cisplatin-like` `bbb_penetrance=0.15`, `binding` left at the 5e-4 default),
**not** structure-derived ones. Those three rows are immune to any logP change. Every
other row, and every call to `calibration.calibrated_drugparams(...)` — which is what
`experiments/real_drug_panel.py` and the breast panel actually use — is not.

---

## 5. Sensitivity: how far does the output move?

Harness: `experiments/logp_sensitivity.py`. Metric: the analytic reaction-diffusion
penetration length `L = sqrt(D_eff / k_eff)`, `k_eff = clearance + alpha*binding` —
the same quantity `validation.calibration` is calibrated against (doxorubicin,
Primeau et al. 2005, PMID 16361566; 40-50 um half-distance -> ~64.9 um e-folding
length). Plus `bbb_penetrance`, the linear delivery multiplier. ECM density 0.4.

### 5.1 The alpha caveat, resolved

`calibrate_alpha` fits on **doxorubicin only**. Therefore, for all nine audited drugs
(none of which is doxorubicin), **the alpha-fixed and alpha-refit columns are
identical by construction** — refitting alpha under a perturbation to a non-anchor
drug returns the same alpha. Both columns are printed anyway so this is visible
rather than assumed; they agree to every printed digit for all nine. They diverge
only for doxorubicin itself, where the refit column returns exactly 1.000x (the drug
snaps back to its fit target — that is the fit, not a prediction).

### 5.2 Per-drug results (alpha fixed = alpha refit, for all nine)

Baselines are the repo's current Crippen values.

| Drug | Baseline (Crippen) | Scenario | Tier | L (um) | L ratio | bbb ratio | **bbb x L** |
|---|---|---|---|---|---|---|---|
| Cisplatin | +0.195, L=142.93, bbb=0.858 | **-2.53** | measured | 238.94 | 1.672x | 0.203x | **0.339x** |
| | | -2.25 (+1 SD) | measured | 230.15 | 1.610x | 0.259x | 0.418x |
| | | -2.81 (-1 SD) | measured | 247.31 | 1.730x | 0.159x | 0.275x |
| Palbociclib | +2.966, L=58.49, bbb=0.579 | **0.228** (pH 5) | measured logD | 130.52 | **2.232x** | 0.647x | 1.443x |
| | | **1.11** (pH 7) | measured logD | 102.79 | 1.758x | 0.868x | 1.525x |
| | | **2.26** (pH 9) | measured logD | 72.77 | 1.244x | 1.013x | 1.260x |
| | | 0.99 | **refuted** | 106.37 | 1.819x | 0.842x | 1.531x |
| Tamoxifen | +5.996, L=36.50, bbb=0.602 | 7.1 | prediction | 28.31 | 0.776x | 0.446x | 0.346x |
| | | 7.87 | prediction | 24.53 | 0.672x | 0.226x | **0.152x** |
| Epirubicin | +0.001, L=64.92, bbb=0.021 | 1.27 | **propagation** | 44.15 | 0.680x | 1.031x | 0.701x |
| Docetaxel | +3.260, L=28.22, bbb=0.019 | 4.26 | prediction | 22.77 | 0.807x | 0.994x | 0.802x |
| | | 3.2 | prediction | 28.63 | 1.014x | 1.000x | 1.015x |
| Paclitaxel | +3.736, L=28.79 | +/-1 log unit | unsourced | 23.07-37.55 | 0.80-1.30x | ~1.00x | ~0.80-1.30x |
| Vinorelbine | +4.754, L=31.24 | +/-1 log unit | unsourced | 24.76-41.11 | 0.79-1.32x | 0.88-1.10x | ~0.70-1.45x |
| Irinotecan | +4.091, L=47.50 | +/-1 log unit | unsourced | 35.74-64.54 | 0.75-1.36x | 0.72-1.21x | ~0.54-1.64x |
| Ribociclib | +2.796, L=62.00 | +/-1 log unit | unsourced | 45.86-84.40 | 0.74-1.36x | 0.92-0.98x | ~0.68-1.33x |

**Magnitudes, stated plainly:**

- A **1 log unit** logP error moves `L` by roughly **0.75x-1.36x** for every drug —
  a consistent ~1.8x spread across the 2-log-unit window. The mapping is smooth and
  monotone, with no cliff.
- The **largest single move from a non-anchor logP change is palbociclib at 2.23x**
  in `L` (Crippen 2.966 -> measured pH 5 logD 0.228). Even the mildest measured
  endpoint (pH 9, 2.26) moves it 1.24x. **Every measured logD endpoint moves
  palbociclib in the same direction — the repo currently under-predicts its
  penetration.**
- The **largest composite (`bbb x L`) move is tamoxifen at 0.152x**, i.e. a **6.6x
  reduction** in delivered drug — but this is a *prediction-to-prediction* comparison
  and must not be reported as a correction. It is best read as: tamoxifen's simulated
  output is uncertain by ~6.6x because nobody knows its logP, not because the repo's
  value is wrong. Tamoxifen sits at logP 5.996 on the far tail of the
  `bbb_penetrance` Gaussian (centre 2.5, width 1.7), which is why it is the most
  logP-sensitive drug in the panel by a wide margin.
- **Cisplatin's measured value moves the composite to 0.339x** (a ~3x reduction),
  driven mostly by `bbb_penetrance` collapsing 0.858 -> 0.174. Note this cuts against
  the repo's own defence at `config.py:414` that structure-derived `bbb_penetrance`
  cross-checks well against measurement — that defence was made for irinotecan, and
  cisplatin is the counter-example.
- **Docetaxel is nearly insensitive** (0.80-1.01x composite) because it sits deep in
  the region where `bbb_penetrance` is already ~0.019 and flat. The 4.26-vs-3.2
  literature conflict is, for this simulator, close to immaterial. That is a genuine
  null and should stop anyone spending more time chasing that citation.

### 5.3 The anchor re-export (the finding that dominates the rest)

Correcting **doxorubicin** from Crippen 0.001 to the measured HSDB/Hansch **1.27**:

```
alpha 148.6082 -> 68.6368   (0.4619x)
doxorubicin's OWN L: 64.921 -> 64.921 um   (pinned to the fit target by construction -- NOT a prediction)
```

Every other drug moves, with **no change to its own logP**:

| Drug | L before | L after | ratio |
|---|---|---|---|
| cisplatin | 142.93 | 209.26 | 1.464x |
| paclitaxel | 28.79 | 42.35 | 1.471x |
| docetaxel | 28.22 | 41.51 | 1.471x |
| epirubicin | 64.92 | 95.42 | 1.470x |
| vinorelbine | 31.24 | 45.94 | 1.471x |
| tamoxifen | 36.50 | 53.66 | 1.470x |
| irinotecan | 47.50 | 69.82 | 1.470x |
| palbociclib | 58.49 | 85.94 | 1.469x |
| ribociclib | 62.00 | 91.09 | 1.469x |

A near-uniform **~1.47x** global rescale. This is larger than most of the per-drug
effects in section 5.2 and it is invisible to a per-drug audit. It is also the change
most likely to actually happen, since doxorubicin's 1.27 is one of the few genuinely
Hansch-backed values and is already sitting in the other branch's table.

**Recommendation on this point specifically:** doxorubicin's logP correction should
be a deliberate, separately-reviewed recalibration with its own revalidation run
(`calibration.revalidate` against the mitoxantrone out-of-sample anchor, PMID
23680920), not a line in a provenance cleanup. Any `results/` file quoting an
absolute penetration length is invalidated by it.

---

## 6. Recommendation, per drug

The right mechanism is **not** to add a logP field to `DRUG_LIBRARY` — logP is not a
`DrugParams` quantity, and `config.py` is correctly free of it. The right mechanism is
a resolver at the descriptor layer that supersedes the Crippen value where a defensible
measurement exists and records the tier when it does not. That resolver **already
exists on the other agent's branch** (`tmesim/molecular/lipophilicity.py`,
`resolve_logp`). The recommendation is therefore to extend that table rather than
build a second mechanism.

| Drug | Recommendation | Rationale |
|---|---|---|
| **Cisplatin** | **Store as measured, with a hard caveat** — `-2.53 +/- 0.28` | The one real measurement in the panel. Crippen's +0.195 is off by 2.73 log units and Crippen has no business modelling a Pt(II) complex. Caveat must be stored, not just cited: octanol/**0.9% saline**, total Pt by ICP-MS — octanol/saline, not textbook octanol/water |
| **Palbociclib** | **Store as logD-with-pH** — 0.228@pH5, 1.11@pH7, 2.26@pH9 | Only drug where a schema change is clearly justified by the data. Crippen 2.966 is above the entire measured range. Use pH 7 (1.11) as the physiological default if a scalar is required, and record it as logD, never logP |
| **Ribociclib** | **Set to null / unknown; retrieve the logD** | Real OECD 107 logD exists but its pH-resolved values were not retrieved in this audit. Do not substitute anything until they are. Sensitivity is ~1.3x, so this is low-priority |
| **Vinorelbine** | **Set to null / unknown; retrieve the logD** | Same as ribociclib |
| **Irinotecan** | **Set to null / unknown; retrieve the logD** | Same as ribociclib. Note it is also a prodrug (SN-38), already flagged at `config.py:410` |
| **Tamoxifen** | **Keep Crippen 5.996, with an explicit "unknown, +/-2 log units" provenance note. Do NOT adopt 7.1 or 7.87** | Both alternatives are predictions. Replacing one prediction with another while relabelling it "experimental" is the exact failure this audit exists to stop. **This is the highest-value open lead in the audit** (6.6x composite sensitivity): PMID 8877846, Takacs-Novak & Avdeef — a genuine interlab shake-flask/potentiometric study of 23 diverse drugs, paywalled, and **whether tamoxifen is among the 23 was not determined.** Resolving that one question is worth more than the rest of this list combined |
| **Epirubicin** | **Keep Crippen 0.001, with a note; do NOT adopt 1.27** | 1.27 is doxorubicin's value. Adopting it would be propagation. But the note must record the stereoisomer trap (section 3.2a): dox and epi have bit-identical RDKit descriptors, so if doxorubicin is corrected to 1.27 and epirubicin is not, the resulting 1.27-log-unit split is an artifact of table membership. Handle the pair together |
| **Paclitaxel** | **Keep Crippen 3.736, with a note** | Nothing found. Sensitivity is modest (~1.3x) |
| **Docetaxel** | **Keep Crippen 3.260, with a note; do NOT adopt 4.26 or 3.2** | Both unsourced/methodless and mutually contradictory. Sensitivity is the lowest in the panel (0.80-1.01x composite) — **this one genuinely does not matter.** Close the thread |
| **Doxorubicin** *(not in the nine, but decisive)* | **Replace with the measured 1.27 — as a separate, revalidated change** | Genuinely Hansch-backed. But it is the sole `alpha` anchor; see 5.3. Requires its own revalidation run and invalidates absolute penetration numbers elsewhere in `results/` |

### 6.1 Proposed diff

**Not applied.** `tmesim/` was not touched by this audit. Apply by hand.

This targets `tmesim/molecular/lipophilicity.py` **on the other agent's branch**
(`.claude/worktrees/agent-ac4e848d0655b5856/`), extending its existing
`EXPERIMENTAL_LOGP` table, which currently contains none of the nine audited drugs.
It deliberately adds **only** cisplatin (the one measurement) and leaves palbociclib
to a separate logD-aware field, because storing a logD in a table named `logp` is how
this class of error starts.

```diff
--- a/tmesim/molecular/lipophilicity.py
+++ b/tmesim/molecular/lipophilicity.py
@@ EXPERIMENTAL_LOGP
+    # ---- Added by the logP provenance audit (results/logp_provenance.md) ----
+    # The ONLY drug in the nine-drug breast/CDK audit with a real measured
+    # octanol-phase partition coefficient. Every other one of the nine has an
+    # EMPTY Hansch/HSDB route -- see results/logp_provenance.md section 3.1.
+    "cisplatin": {
+        "logp": -2.53,
+        "sd": 0.28,
+        "n": 5,
+        "method": "shake-flask (Albert 1979 method)",
+        "source": (
+            "Screnci D et al., Br J Cancer 2000;82(4):966-972, PMID 10732773, "
+            "PMC2374394, Table 3 (NOT Table 1, which is a neurotoxicity ranking)."
+        ),
+        "confidence": "medium",
+        "note": (
+            "CAVEAT, DO NOT DROP: the aqueous phase was 0.9% NaCl SALINE, not "
+            "pure water, and platinum was quantified as TOTAL Pt by ICP-MS. "
+            "This is an octanol/SALINE partition coefficient, not a textbook "
+            "octanol/water log Kow. It is still far closer to reality than the "
+            "RDKit Crippen value it supersedes (+0.195): Crippen is a "
+            "topological atom-contribution model with no parameterisation for a "
+            "square-planar Pt(II) complex. Superseding moves cisplatin's "
+            "simulated penetration length 142.93 -> 238.94 um (1.67x) and its "
+            "bbb_penetrance 0.858 -> 0.174 (0.20x); net delivered drug 0.34x. "
+            "See results/logp_provenance.md section 5.2."
+        ),
+    },

@@ (new module-level table, adjacent to EXPERIMENTAL_LOGP)
+# ---------------------------------------------------------------------------
+# Measured logD (OECD 107) -- a DISTINCT quantity from logP. Stored separately
+# and keyed by pH ON PURPOSE: silently collapsing a pH-dependent logD into a
+# field named "logp" is precisely the error this table exists to prevent.
+# Consumers must opt in explicitly and state the pH they are assuming.
+# ---------------------------------------------------------------------------
+EXPERIMENTAL_LOGD: Dict[str, Dict] = {
+    "palbociclib": {
+        "by_ph": {5.0: 0.228, 7.0: 1.11, 9.0: 2.26},
+        "method": "OECD 107 (shake-flask)",
+        "source": "EMA Ibrance EPAR, EMA/652627/2016, p.28.",
+        "confidence": "medium",
+        "note": (
+            "REFUTATION ON RECORD: the widely-propagated palbociclib logP of "
+            "0.99 could NOT be found in the document it is attributed to. The "
+            "Health Canada IBRANCE monograph (rev. 2024-09-04) section 13 "
+            "states only pKa 7.4 and 3.9 plus pH-dependent solubility -- no "
+            "partition coefficient and no '0.99'. FDA label 207103s008 is "
+            "identical. DrugBank's citation does not check out against the "
+            "document it names. None of the three measured logD values equals "
+            "0.99. Do not reintroduce it from any secondary aggregator.\n"
+            "The repo's RDKit Crippen value (+2.966) is ABOVE the entire "
+            "measured range including the pH 9 endpoint. Superseding with the "
+            "pH 7 value (1.11) moves palbociclib's penetration length "
+            "58.49 -> 102.79 um (1.76x); the pH 5 value gives 2.23x, the "
+            "largest non-anchor logP effect in the audited panel."
+        ),
+    },
+}
```

Plus a provenance-note-only change for the four drugs that stay on Crippen, so the
"unknown" is recorded rather than implied:

```diff
--- a/tmesim/molecular/breast_library.py
+++ b/tmesim/molecular/breast_library.py
@@ module docstring
+LIPOPHILICITY PROVENANCE (results/logp_provenance.md, audited)
+-------------------------------------------------------------
+Every logP this module produces is an RDKit Crippen MolLogP PREDICTION derived
+from the SMILES below. The SMILES are verified against PubChem; the logP is
+not, and for these drugs it CANNOT be: a literature audit found the
+Hansch/HSDB measured-log-Kow route EMPTY for all of paclitaxel, docetaxel,
+epirubicin, vinorelbine, tamoxifen and palbociclib.
+
+  paclitaxel   Crippen +3.736  no measured value found anywhere.
+  docetaxel    Crippen +3.260  4.26 (EMA) is methodless ASMF boilerplate quoted
+                               twice; Sanofi's own Taxotere SDS conflicts at
+                               3.2. Neither adopted. Sensitivity is the LOWEST
+                               in the panel (~1.0x) -- this does not matter.
+  epirubicin   Crippen +0.001  the commonly-quoted 1.27 is DOXORUBICIN's HSDB
+                               value (Hansch p.187); epirubicin is its
+                               4'-epimer and has no field of its own. NOT
+                               adopted -- that would be propagation. NOTE:
+                               Crippen is stereochemistry-blind, so epirubicin
+                               and doxorubicin have BIT-IDENTICAL descriptors
+                               here. Correcting doxorubicin without epirubicin
+                               creates a 1.27-log-unit split between two
+                               stereoisomers that is an artifact of table
+                               membership, not chemistry. Handle as a pair.
+  vinorelbine  Crippen +4.754  logD (OECD 107) exists; values not retrieved.
+  tamoxifen    Crippen +5.996  7.1 is a PREDICTION laundered into an
+                               "experimental" PHYSPROP/OPERA field (verified
+                               from the source files); 7.87 +/- 0.75
+                               (El-Kattan 2001, PMID 11250108) is a calculator
+                               output -- the same table reports nicardipine
+                               HYDROCHLORIDE logP -0.99, impossible for a salt.
+                               NEITHER is a measurement; neither adopted.
+                               OPEN LEAD: PMID 8877846 (Takacs-Novak & Avdeef)
+                               is a genuine interlab shake-flask study of 23
+                               drugs; full text paywalled and whether tamoxifen
+                               is among them was NOT determined. This is an
+                               open lead, not a negative result -- and it is
+                               the highest-value one in the audit, because
+                               tamoxifen is the panel's most logP-sensitive
+                               drug (up to 6.6x in delivered drug).
+  palbociclib  Crippen +2.966  the propagated 0.99 is REFUTED (see
+                               lipophilicity.EXPERIMENTAL_LOGD). Real OECD 107
+                               logD exists: 0.228/1.11/2.26 at pH 5/7/9.
+                               Crippen sits above that whole range.
```

### 6.2 Sequencing

1. **Do first, cheap and safe:** land the provenance notes (6.1 second diff). They
   change no number and stop the next reader from assuming these values are measured.
2. **Do next:** add cisplatin's measured value and palbociclib's `EXPERIMENTAL_LOGD`.
   Both are well-evidenced and neither touches `alpha`.
3. **Do separately, with revalidation:** doxorubicin -> 1.27. This is a global
   recalibration (~1.47x on every penetration length), not a provenance fix. Requires
   `calibration.revalidate` against the mitoxantrone anchor (PMID 23680920) and a
   sweep of `results/` for invalidated absolute lengths. Decide epirubicin in the same
   change, not before or after.
4. **Retrieve, do not guess:** the OECD 107 logD values for ribociclib, vinorelbine,
   irinotecan; and PMID 8877846's compound list for tamoxifen.
5. **Close:** docetaxel. The literature conflict is real and, for this simulator,
   immaterial.

---

## 7. Reproducing this

```
./.venv/bin/python experiments/logp_sensitivity.py
```

Compiles 10 drugs once each (the only expensive step; ~1 Brownian sim panel per drug),
then evaluates all logP scenarios in closed form. Prints the `DRUG_LIBRARY`
reproduction check, per-drug scenarios with both alpha columns, the composite
`bbb x L` table, and the anchor re-export table. Every alternative logP in the script
is tagged with its tier (`measured` / `measured_logD` / `refuted` / `propagation` /
`prediction` / `unsourced`) and its citation, in the `SCENARIOS` dict.

**Discipline note.** No number in this document was invented. Every literature value
traces to the audit table supplied with its own PMID/document reference; every
simulator number is output from the run above. Where a value exists but was not
retrieved (the three logD sets, PMID 8877846's compound list), it is recorded as an
open lead and left blank rather than filled in — and no prediction is restated
anywhere here as a measurement.
