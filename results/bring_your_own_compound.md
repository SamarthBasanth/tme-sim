# Bring your own compound — what a self-compiled molecule is, and what it is not

**What this file is.** The shipped drug picker has thirteen curated molecules (the twelve
in `molecular.library.GBM_DRUGS` plus the cross-library palbociclib). A researcher's
own series is the reason they are here, so the Studio now accepts a SMILES and compiles
it into the same `core.DrugParams` a library drug is compiled into. This file records
what that path actually produces, how trustworthy each field is, what the compiler
refuses and why, and the exact claim a user is permitted to make from the result.

**Code:** `tmesim/molecular/byo.py` (new), `app/server/compounds.py` (new),
`app/web/js/ui/compound.js` (new), `tests/test_byo_compound.py` (new, 80 tests).
Wiring only in `app/server/api.py`, `app/server/jobs.py`, `app/web/js/api.js`,
`app/web/js/main.js`, `app/web/js/ui/spec.js`, and one function exported from
`app/web/js/ui/inspect.js`. Nothing in `tmesim/core`, `tmesim/physics`,
`tmesim/validation` or `tmesim/reconstruct` was modified.

---

## 0. Bottom line

1. **A user compound is compiled by the SAME chain as a library drug, with the same
   `alpha`** — pinned by an exact-equality test against `calibrated_drugparams ->
   apply_barrier_to_drugparams` (§4). Two compounds calibrated differently cannot be
   compared at all, so this is the precondition for everything else here.

2. **8 of 11 compiled fields are predictions.** Only `mw` and `tpsa` are exact functions
   of the graph the user typed. `ic50` is the user's own number or a flagged placeholder.
   Everything else — logP, hydrodynamic radius, free and tissue diffusivity, BBB
   penetrance, binding, clearance, Hill — is model output. The Studio leads with that
   ratio rather than burying it.

3. **The transport ranking is FALSIFIED for analogues of similar molecular weight** —
   which is the usual shape of a real compound series, and therefore the sharpest thing
   in this file. Diffusivity is a function of MW alone; `results/diffusivity_falsification.md`
   tested that against 71 measured compounds and refuted it. §3.

4. **Absolute kill, AUC and penetration depth remain non-defensible**, exactly as for a
   library drug. Supplying the molecule yourself changes nothing about that; the UI says
   so at the point of use, in the caveats panel, in the picker option, in the RunSpec
   validator's warnings, and on the finished run.

5. **logP is an RDKit Crippen prediction — for a user compound and for every curated drug
   alike.** There is exactly one logP route in this tree, so a self-compiled compound is
   on the same footing as a library one here. A name-keyed *measured* logP tier has been
   proposed for this repo more than once; `byo._check_logp_tier` is an inert backstop that
   refuses a user compound the day one lands, because such a lookup is safe for a vetted
   library and unsafe for a form field. §3.2.

---

## 1. The path, end to end

```
SMILES ──► descriptors.compute_descriptors        (RDKit; RAISES on bad input)
       ──► descriptors.from_smiles                (adds the MW-derived radius)
       ──► compile_drug.build_profile             (Stokes-Einstein + Brownian + ECM mesh)
       ──► compile_drug.molecule_to_drugparams    (-> core.DrugParams)
       ──► calibration.calibrated_drugparams(dp, alpha)   (loss terms <- the one fit)
       ──► binding_site_barrier.apply_barrier_to_drugparams
       ──► core.DrugParams  ── the object the engine integrates
```

`byo.compile_user_compound` composes those steps; it reimplements none of them. The last
two are the ones a naive "just call `molecule_to_drugparams`" implementation would omit,
and omitting them would put a user compound on a different loss scale from every library
drug — i.e. it would silently destroy the ordering, which is the only output either has.

### 1.1 The pre-existing path could not be exposed as-is

`descriptors.from_smiles` is deliberately crash-proof: given an unparseable SMILES and an
unknown name it returns a `Molecule` with `mw = 0.0`, every descriptor zero, and
`source='fallback_table'`. Downstream that becomes a `DrugParams` with a degenerate
diffusivity rather than an error. That is correct for a batch pipeline that must not die
mid-panel and **exactly wrong for a form field**. `compile_user_compound` therefore calls
`compute_descriptors` (which does raise) and refuses on the fallback path outright. Both
directions are pinned by
`test_the_from_smiles_fallback_can_never_stand_in_for_a_user_compound`.

---

## 2. Every compiled field and what it really is

Worked example: aspirin (`CC(=O)Oc1ccccc1C(=O)O`) named `ACME-101`, IC50 12.5 µM supplied.

| field | value | unit | Studio chip | status |
|---|---:|---|---|---|
| `mw` | 180.159 | g/mol | `computed` | exact from the graph |
| `tpsa` | 63.6 | Å² | `computed` | exact from the graph |
| `ic50` | 12.5 | µM | `literature` | the user's own measurement, on their cells |
| `logp` | 1.310 | — | `assumed` | **RDKit Crippen prediction** |
| `radius_hydrodynamic_nm` | 0.3728 | nm | `assumed` | **from MW alone** |
| `free_diffusivity` | 507.6 | µm²/s | `assumed` | Stokes-Einstein on the predicted radius |
| `diffusivity` | 286.1 | µm²/s | `assumed` | **from MW alone**, end to end |
| `bbb_penetrance` | 0.9596 | frac | `assumed` | CNS-MPO-style rules |
| `binding` | 1.127e-02 | 1/s | `assumed` | overwritten by `alpha × uptake_proxy` |
| `clearance` | 9.73e-05 | 1/s | `assumed` | weak logP/MW heuristic |
| `hill` | 1.2 | — | `assumed` | a constant, for every compound |

`binding` is quoted at this tree's `alpha` (148.608), refit from doxorubicin on every
process start. It is not a constant of the repo — during this work a change to
doxorubicin's logP moved it to 68.637 and back. Every compound's `binding` moves with it
*together*, which is exactly why an ordering survives such a change and a level does not.

**Nothing is `measured`.** No instrument touched this compound in this build. The chip
vocabulary is the Studio's existing one (`app/web/js/ui/provenance.js`), reused rather
than reinvented, so a chip means the same thing here as it does on a metric.

### 2.1 Diffusivity is a function of molecular weight and nothing else — see §3

`r = 0.066 · MW^(1/3)` → Stokes-Einstein → Brownian micro-tortuosity → Ogston-Renkin mesh
gate. Shape, logP, charge and H-bonding never enter it. Asserted, not asserted-about:
n-butylbenzene and *tert*-butylbenzene (C₁₀H₁₄ both ways) compile to **bit-identical**
diffusivity while their Crippen logP differs
(`test_diffusivity_is_a_function_of_mw_alone`). Two isomers in a user's series are
indistinguishable to the transport model by construction.

### 2.2 The hindrance law is still short of measurement

`results/cellular_exclusion.md` reaches 6.6× of a measured 12× hindrance — a residual
1.8× gap that no compound in this repo closes. A user compound inherits that gap exactly.
It is printed on the `diffusivity` row, with the file name, so the operator can go read
it.

### 2.3 IC50 is the user's or it is absent

Structure cannot give potency without docking or an assay. Supplied → used exactly as
given and labelled the user's own (`literature` chip: real, but measured on other
material). Absent → a flagged 25 µM placeholder, and the compile emits a caveat saying
that any kill fraction or coverage number derived from it **must not** be read as
potency, and that a ranking including that compound is a **delivery** ranking only.

---

## 3. The transport ranking is falsified for similar-mass analogues — the finding that matters most

A user's series is almost always **congeneric**: a scaffold with substituent variations,
so the compounds cluster in molecular weight. That is exactly the case the transport model
has been measured against and **refuted**.

`results/diffusivity_falsification.md` tested `MolecularProfile.effective_diffusivity`
against 71 compounds whose tissue diffusion coefficients were measured in HT29
multicellular layers by one method in one cell line (Pruijn FB, Patel K, Hay MP, Wilson
WR, Hicks KO, *Aust J Chem* 2008;61(9):687–693, [DOI](https://doi.org/10.1071/CH08240)) —
the same quantity in the same equation, so commensurability is not in question:

| | measured | this engine |
|---|---:|---:|
| spread across the 71 compounds | **29.7×** | **2.4×** |
| correlation with the measurement | — | **R² = 0.028** |
| spread at matched MW, n=4 at 263 Da | **3.7×** | **1.00×** |
| spread at matched MW, n=4 at 305 Da | **3.8×** | **1.00×** |

Molecular weight — diffusivity's **only** molecular input here — explains **2.8 %** of the
variance in the measured values. At matched MW the engine is flat *by construction*, and
the measurement is not.

**This is not a precision caveat, and it is not hedged as one.** It appears in three
places a user cannot miss: on the `diffusivity` row itself (with the numbers and the file
name), in the CANNOT-conclude list, and as a hard-stop caveat worded to trip the Studio's
`MUST NOT` renderer. `test_the_diffusivity_falsification_is_disclosed_at_the_point_of_use`
pins all three.

**What it costs the product claim.** The CAN list no longer offers an unqualified
"rank your compounds". It offers an ordering *of this model's inputs*, explicitly flagged
as "not a validated prediction", and it promotes the one thing the compiler does most
reliably: telling a user that two of their compounds are **indistinguishable to it**. For
a congeneric series that is frequently the honest answer, and it is now the headline of
the CAN list rather than a footnote.

A difference driven by an **IC50 the user measured** stands on much firmer ground than one
driven by transport, and the field table separates the two so that distinction is visible.

### 3.1 Crippen is stereochemistry-blind

Doxorubicin and epirubicin differ at one stereocentre and produce **bit-identical**
descriptor dicts — hence identical diffusivity, BBB penetrance and binding
(`test_crippen_is_stereochemistry_blind`). **If a series varies by stereochemistry, this
model cannot rank it at all.** Combined with §3, the two failure modes cover most of what
a medicinal chemist actually varies within a series.

### 3.2 Every logP is a Crippen prediction, and a measured tier must not be name-keyed

`results/logp_provenance.md` §1 established that every logP in this repo is an RDKit
Crippen prediction and that the Hansch/experimental route is empty for all nine drugs it
audited. **That is still true in this tree**, so a self-compiled compound is on the *same*
logP footing as a curated drug — no asymmetry is claimed, because none exists.

A name-keyed *measured* logP tier has been proposed for this repo more than once. It is
safe for a curated library (somebody checked that the name and the structure are the same
drug) and **unsafe for a form field**: a user could type any structure, call it `aspirin`,
and receive aspirin's published logP — which then propagates into `plasma_free_fraction`,
`tissue_binding_rate`, `clearance_rate`, `membrane_permeability` **and** `bbb_penetrance`
before anything else in the pipeline runs. `byo._check_logp_tier` refuses a user compound
carrying `logp_tier == 'experimental'` (`name_matches_measured_logp`), with a rename as the
remedy. It costs one dict lookup and is inert today;
`test_a_name_keyed_measured_logp_would_be_refused` is what makes it non-inert the day such
a tier lands.

---

## 4. Parity with the library, asserted by exact equality

`test_same_structure_same_alpha_gives_identical_transport` compiles lomustine's own
SMILES under the name `MY-LOMUSTINE` and demands `diffusivity`, `molecular_weight`,
`binding`, `clearance`, `bbb_penetrance`, `ic50` and `hill` all equal — **exactly**, not
within a tolerance — to
`apply_barrier_to_drugparams(calibrated_drugparams("lomustine", alpha))`.

Lomustine is chosen because it has no by-name special case (no TMZ hypoxia override, no
cross-library reuse), so the only remaining difference between the curated compile and the
user compile is the name — which must not matter.

`test_logp_is_the_same_tier_for_a_user_compound_and_a_library_drug` guards the assumption
underneath that: it asserts no library drug carries a `logp_tier`, so if a name-keyed
measured tier ever lands, this parity test's premise fails loudly rather than silently.

`test_calibrated_drugparams_passthrough_contract_holds` separately pins the contract this
composition depends on: `calibration._compile` documents that a `DrugParams` passes
through unchanged, so `calibrated_drugparams(dp, alpha)` applies `binding = alpha ×
uptake_proxy(dp)` and `clearance = k_passive(dp)`. If that ever changes, the test fails
here with a clear reason instead of quietly putting user compounds on a different loss
scale.

---

## 5. What is refused, and why each refusal exists

Each raises `CompoundRefused(code, message, field)` → HTTP 422 with the sentence on the
field at fault. Never a stack trace.

| input | code | why it is refused rather than compiled |
|---|---|---|
| `(blank)` | `smiles_empty` | nothing to compile |
| `not-a-smiles`, `InChI=...`, `c1ccccc` | `smiles_unparseable` | `from_smiles` would return a zeroed molecule |
| >600 chars | `smiles_too_long` | a biologic or a polymer |
| `EVQLVESGGGLVQPGG...` | `smiles_is_a_sequence` | a one-letter amino-acid paste; names the real mistake instead of "unparseable" |
| `CC(=O)[O-].[Na+]` | `smiles_multi_component` | MW would include the counter-ion, and diffusivity is MW alone — the salt would travel bonded to the drug |
| `[Gd]`, `[U]` | `crippen_unparameterised_atom` | Crippen does not fail on an unknown element, it assigns contribution **exactly 0.0** — a confident logP computed as if those atoms were absent |
| `C`, `CCO` (MW < 50) | `mw_below_range` | the `0.066` radius prefactor was set at ~200 g/mol |
| MW > 2000 | `mw_above_range` | outside the band `compile_drug`'s own bridge self-test asserts |
| MW ≥ 10 000 or `is_antibody` | `biologic_out_of_range` | the radius law switches to a globular-protein fit and every descriptor loses meaning; the library's one antibody (bevacizumab) is **hand-entered protein values, not a compiled structure**, which is precisely what a SMILES cannot reproduce |
| name `temozolomide`, `paclitaxel`, `mitoxantrone`, … | `name_collision` | a run record must say unambiguously which structure it simulated |
| a molecule carrying `logp_tier='experimental'` | `name_matches_measured_logp` | §3.2 — inert in this tree; the guard that fires the day a name-keyed measured tier lands |
| `../etc/passwd`, `has/slash` | `name_charset` | it is also a path component in the store |
| IC50 `0`, `-3`, `NaN`, `banana`, `5e6` | `ic50_not_positive` / `ic50_not_finite` / `ic50_not_a_number` / `ic50_absurd` | the Hill term divides by it; `5e6 µM` is 1 M, i.e. a unit error |
| any non-finite compiled parameter | `compile_nonfinite` | parses, but produces no usable physics |

**In range but extrapolating is a warning, not a refusal.** The two curated libraries span
130.078 g/mol (5-fluorouracil) to 853.918 (paclitaxel). Inside `[50, 2000]` but outside
that span, the MW→radius law is extrapolating; the compile succeeds and says so.

---

## 6. What the Studio lets a user do, and what it tells them

* **`Setup → Define → Your compound`** (collapsed by default): Name, SMILES, optional
  IC50, optional note. **Compile** is a dry run that stores nothing; **Add to picker**
  stores it and it becomes a legal `RunSpec.drug`.
* **The headline** above the provenance table: *"8 of 11 compiled fields are PREDICTIONS
  from the structure you typed, not measurements."*
* **One row per field** with its value, its provenance chip and its method sentence. The
  predicted rows carry a warning rule down the left edge.
* **Caveats** rendered by the same component the result panel uses, so the hard stops
  (`NOT DEFENSIBLE`, `MUST NOT`) get the same stop glyph and tinted field they get on a
  finished run. Four are emitted on every compile: the transport falsification (§3), the
  Crippen logP prediction, the self-compiled-compound blanket limit, and the IC50's
  origin.
* **CAN / CANNOT conclude**, verbatim from `byo.CAN_CONCLUDE` / `CANNOT_CONCLUDE`, so the
  wording in the module and the wording on screen cannot drift.
* **The drug picker option** is marked `[YOURS · 8/11 fields predicted]` with the full
  caveat on hover — because by the time a compound is chosen the compile screen may be
  closed.
* **`POST /api/validate`** returns a warning naming the compound as self-compiled, telling
  the operator to read the ordering rather than the levels, and carrying the falsification
  numbers (2.4× vs 29.7×, R² = 0.028) so the qualification travels with the permission.
* **The finished run** carries `drug_params.structure_source = "user_smiles"`,
  `provenance.user_compound` (SMILES, the whole field table, the can/cannot lists) and
  the compile-time caveats in `provenance.caveats`, where the existing caveats panel
  renders them alongside the standing limits.
* **Stub mode refuses to compile at all.** If the backend is not answering, a fabricated
  `DrugParams` would look exactly like a real one; the panel says so instead.

### 6.1 The permitted claim, in full

**CAN:** compare this compound against others on the **same section** with the **same
RunSpec** and read the **ordering**, not the levels (`sensitivity_mode='uniform'` makes the
per-cell sensitivity identical for every compound, so it cancels out) — understanding that
this is an ordering of *this model's inputs* and **not a validated prediction**; see
*which* property separates two compounds, and weigh a difference driven by an IC50 you
measured far above one driven by transport; **detect that two of your compounds are
indistinguishable to this model**, which is the most reliable thing the compiler does and
is often the honest answer for a congeneric series.

**CANNOT:** read its absolute kill fraction, absolute AUC/exposure, or absolute
penetration depth in µm; **rank a set of analogues of similar molecular weight on
transport** (§3 — measured and refuted, not merely uncertain); conclude anything about a
stereoisomer, salt form, prodrug or formulation; conclude anything clinical.

---

## 7. Storage

`app/compounds/<id>/compound.json`, written atomically, mirroring `app/server/projects.py`.
**Only the inputs are persisted** — name, SMILES, supplied IC50, note. Every compiled
number is recomputed on read, so a physics change propagates instead of leaving a stale
parameter set on disk that still looks current. Pinned by
`test_only_the_inputs_are_persisted`.

---

## 8. Tests

`tests/test_byo_compound.py` — **80 tests**, in four groups:

1. **Refusals** — `TestRefusals` (35). Every code in §5, that each message is a complete
   sentence over 40 characters, that the field named is one the form has, that no stack
   trace reaches the operator, and that the crash-proof `from_smiles` fallback can never
   stand in for a user compound.
2. **The honesty claims asserted numerically** — `TestProvenanceIsTrue` (15). logP *is*
   the Crippen number and the recorded tier agrees; Crippen *is* stereochemistry-blind;
   diffusivity *is* a function of MW alone; the caveats contain the exact tokens the
   Studio's hard-stop renderer keys on; can/cannot are stated and disjoint.
3. **Parity with the library** — `TestParityWithTheLibrary` (5). Exact equality on a
   `rdkit_crippen`-tier drug, the `calibrated_drugparams` pass-through contract, the
   measured-logP refusal, and the uncalibrated-preview disclosure.
4. **The Studio wiring** — `TestStore` (10), `TestEndpoints` (10),
   `TestRunSpecValidation` (3), `TestEndToEndRun` (2). Store lifecycle, restart survival,
   the in-process-only compile memo, the endpoints and their 422s, the drug catalogue, the
   RunSpec validator, and two **real end-to-end runs** on `xenium_319`: one user compound
   producing a completed run whose result says what it is, and one user compound ranked
   against temozolomide on the same section.
