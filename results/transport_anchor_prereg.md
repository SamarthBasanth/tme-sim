# Pre-registration: the next out-of-sample test of the transport law

**Status: REGISTERED, NOT RUN.** Nothing in this document is a result. No number
below is an outcome. Every number is either a threshold fixed in advance, or a
value already published elsewhere in this repository and cited to its file and
line.

**Three tests are registered here — one live, one blocked, one awaiting data.**

| | test | status |
|---|---|---|
| **B** | **Rat brain-slice tissue partitioning** — measured `Vu,brain` across six open-access papers sharing the Uppsala slice method, against the engine's implied unbound tissue distribution volume | **PRIMARY, live** |
| **A** | **A second cell line** — Pruijn's SiHa / HCT8-Ea / HCT8-Ra with shipping constants unchanged | **BLOCKED**, registered anyway (§5.1) |
| **C** | **Intracellular partitioning** — `Kp,uu,cell` against `cell_partition = 1.0` | **REGISTERED, NOT RUNNABLE** — 4 scoreable brain values (§3.8). The intracellular accumulation term **remains unconstrained**, and that is the honest status |

**This file was renamed twice.** Commissioned as `results/pruijn_prequel.md`, a
pre-registration for scoring the engine against the 71 HT-29 compounds of Pruijn
2008. That test cannot mean what it was commissioned to mean. §0 says why, and §0
is the most important section here for any later reader.

---

## 0. The 71 HT-29 compounds are spent

### 0.1 The sentence a future reader must not be able to blur

> **The 71 HT-29 compounds of Pruijn 2008 are the data the shipping transport law
> was selected on. Agreement with them is a consistency check, never external
> validation, and no claim that the model generalises may cite them.**

### 0.2 The evidence, from this repository's own record

| commit | file | what it did |
|---|---|---|
| `d3b9955` | `results/diffusivity_falsification.md` | Scored the then-shipping molecular-weight-only `D_eff` blind against the 71 and **falsified it** — R² = 0.022, modelled spread 2.3× against a measured 29.7×, flat by construction at matched MW. |
| `fd2a9e7` | `results/two_route_promotion.md` | Built the two-route cell term **because of that failure**, scored it on the **same 71** at R² = 0.718, and flipped it on as the shipping default. `calibrate_alpha()` moved 148.60815765445545 → 22.223477711172563. |

Three architectural decisions were made against this dataset:

1. **The membrane leg exists because of it** — §7 of the falsification measured
   that the spread tracks hydrogen-bond donors (R² = 0.305) and logP₇.₄
   (R² = 0.290), not MW (R² = 0.028). The transcellular route was built to consume
   exactly those two descriptors.
2. **Intracellular retardation is off because of it** —
   `results/intracellular_retardation.md` built the term and **refuted** it on
   these compounds (R² 0.718 → 0.606; matched-MW rank agreement 77.4% → 51.6%).
3. **Efflux is off partly because of it** — 0 of 71 compounds carry an efflux
   assignment, so the term would move the score by exactly zero.

A dataset that decided which mechanisms exist in a model is not available to
validate that model. The mitigations were real — `tworoute_prereg.json` was frozen
before scoring, carries a binding *"DO NOT FIT TO PRUIJN"*, and the one parameter
that moved went **away** from the Pruijn optimum — but a pre-registration protects
against parameter tuning; it does not un-spend a dataset that selected structure.

### 0.3 What Pruijn-HT-29 may still be used for

**Spent is not useless.** It remains a **regression guard**, a **consistency
check**, and the home for the **matched-MW degeneracy diagnostic**. If ever
registered in `codon validate`, its status must be `NOT TESTABLE`, reason
*in-sample — the transport law's structure was selected on this dataset*, exactly
as the doxorubicin penetration anchor is already declined in
`scorecard.py::_run_perivascular_penetration`. **One exception:** a component never
scored on these compounds is not spent with respect to them — §6.

### 0.4 Two corrections to the record, and a counting ambiguity

Recorded here because a later reader should not have to rediscover them.

1. **"Asymptote coefficients move up to 100-fold" is wrong as a summary.** The
   phrase is the authors' but refers to HCT8-Ea's increase at the high-logP₇.₄
   asymptote, not to a spread in fitted coefficients. The actual spread in `a` is
   **1.35 log₁₀ ≈ 22×**. This document uses 22×.
2. **SiHa is weaker than it appears.** Its `c` is printed identical to HT-29
   *including the standard error*, which strongly implies `c` was never re-fitted
   for SiHa and only `a` moved.
3. **The paper's own compound counts do not reconcile.** Table S1 lists **71**;
   Table 1 states **n = 69**; three rows are marked *"in table but not in graph"*;
   three carry no logP at all; 71 − 3 = **68**, which is the count the repo's
   scoring run uses. **The paper never says which compounds were dropped.** Any
   scoring against these 71 inherits that ambiguity, and it must be stated in any
   report that quotes an n.

---

## 1. Which test to run first — the judgement, argued

The coordinator asked which failure I would rather discover first. The practical
question is settled by A being blocked, but the answer decides whether an email to
Pruijn's corresponding author is worth sending, so here it is.

### 1.1 The two failures are not the same kind of thing

* **A failing** — the law is a description of HT-29 rather than a transport model.
  This is **evaluative**: it tells you the model does not transfer, without telling
  you which term is wrong.
* **B failing** — the engine's tissue drug-capacity assumption, never checked
  against any measurement, is in the wrong place. This is **diagnostic**: it names
  the term. (Test **C**, if its data arrives, does the same for the intracellular
  partition `K` specifically.)

### 1.2 My answer: B's failure is the more valuable discovery, and I would want it first

Four reasons, and one concession.

**(i) B can reveal that A would have passed for the wrong reason; A cannot return
the favour.** The membrane leg was selected on HT-29 **with every tissue-capacity
term frozen** — `cell_partition` at 1.0 and retardation off. If tissue capacity is
in fact strongly compound-dependent and correlated with lipophilicity — the
direction basic lipophilic compounds are expected to go, and **6 of B's 9 compounds
are bases** — then the logP/HBD membrane term has been **absorbing that
dependence**. R² = 0.718 would then be partly a proxy fit: the right answer reached
through the wrong term. That
is `docs/VERIFICATION_METHODOLOGY.md` **pattern 1, compensating errors**, described
there as *"the most dangerous defect class in scientific code."* A model with two
errors that cancel can transfer across cell lines and still be mechanistically
wrong. **A would pass it. B is the only one of the two that can see it.**

**(ii) B has more discriminating power, and now it has n as well.** On the set
actually scored (§3.4 E6, cyclic peptides excluded) the measured `Vu,brain` spread
is **1308-fold, 3.1 decades** against the **29.7×** of the dataset the law was built
on — **44× wider**, at **n = 35 rat values**. Under
`docs/VERIFICATION_METHODOLOGY.md` §5, discriminating power is signal against noise
floor, and B now has both a wider signal and — for the first time in this project —
**a measured noise floor to compare it against** (§3.4a).

**(iii) The failures have asymmetric follow-on cost.** A failing sends you back to
"find another dataset" — the state we are already in. B failing **names the fix**:
constrain the tissue-capacity coefficient. B passing **removes a live objection**
that currently caveats every transport number the engine produces.

**(iv) B is available and A is blocked.** Secondary, but real.

**The concession, stated because it cuts against my answer.** A tests
**generalisation**, which is the property actually claimed. B tests **one input
coefficient**, in a different organ and a different species. A's PASS would be
worth considerably more as evidence *for* the claim than B's PASS. The honest
formulation:

> **B is the better failure to discover. A is the better pass to obtain.**

Since the purpose is to find out whether the model is wrong rather than to
accumulate support, the more informative failure comes first — **B**.

### 1.3 The practical consequence

**The email to Pruijn is worth sending regardless of the odds.** A's PASS is the
only one of the two that could ever support a generalisation claim, and no quantity
of B substitutes for it. B tells us whether the mechanism is right; only A tells us
whether it transfers. Buying A is buying a different thing, not a better version of
the same thing.

---

## 2. Common ground: what these tests do and do not touch

### 2.1 The scope fence, applied to both

Both tests score **input coefficients of the transport PDE, not the solver and not
any solved field**. By `results/anchor_registration.md` §2's own convention that
makes them `analytic`, the same category as the quick anchors that never invoke the
transport PDE. `scorecard.py` states the rule directly:

> *"These modules must be reported as UN-TESTABLE, not as unaffected. 'Robust under
> exchange' and 'never evaluates the delivery law' are very different statements and
> must not be pooled in a coverage count."*

Registered consequences for both tests:

* scorecard `scope` = `SCOPE_ANALYTIC`, **never** `SCOPE_PDE`;
* conditional `coefficient-not-field` attaches, which states the row *"can FALSIFY
  the hindrance law, and it cannot confirm or refute the solver"*;
* **no result from either test may support** any claim about a solved concentration
  field, vessel spacing, transvascular delivery, interstitial convection, **any
  absolute kill fraction** (never a defensible output of this engine), or anything
  clinical or diagnostic;
* a run must record, via a call counter on `physics.diffusion.simulate_drug_transport`
  and `sim.engine.simulate_response`, that the solver was **not entered** — the
  instrumentation pattern already used in
  `results/transvascular_exchange_revalidation.md` §6.1.

### 2.2 The entry point

```
tmesim.molecular.descriptors.from_smiles(name, smiles)     -> Molecule
tmesim.molecular.compile_drug.build_profile(mol, ecm_reference=0.4)
                                                            -> MolecularProfile
```

Test B additionally reads `tmesim/molecular/two_route.py`'s `cell_partition`
parameter, shipped at **1.0**.

### 2.3 FROZEN in every scored arm of both tests

| quantity | shipped value | where |
|---|---|---|
| `cell_partition` (K) | **1.0** | `two_route.py:357`, `tworoute_prereg.json` |
| `HBD_LOGP_DECREMENT` | **1.34** | `two_route.py:245` |
| `cell_volume_fraction` | **0.7892** | `tworoute_prereg.json` |
| `cell_radius_um` | **7.8517** | `tworoute_prereg.json` |
| `cytoplasm_D_over_water` | **0.27** | `tworoute_prereg.json` |
| `membrane_viscosity_pa_s` / `bilayer_thickness_nm` / `logk_membrane_max` | **1.0 / 4.0 / 3.0** | `tworoute_prereg.json`, `two_route.LOGK_MEMBRANE_MAX` |
| `ecm_reference` | **0.4** | `anchor_prototype.py:84` |
| `_TWO_ROUTE_DEFAULT` | **True** | `compile_drug.py:90` |
| `alpha` | **22.223477711172563** | not in either test's loop; see §9.4 |

> **DO NOT FIT TO THE ANCHOR.** No parameter, exclusion, threshold or descriptor
> choice may be changed by looking at a held-out residual.

---

## 3. TEST B — rat brain-slice tissue partitioning *(PRIMARY)*

### 3.1 The dataset, and the cascade that defines the scored set

**Source.** Six open-access primary papers sharing the **identical Uppsala brain-slice
method**. Method homogeneity is what makes them poolable at all; it is a property of
the sources and must be re-verified per paper at GATE 1, not assumed from this
sentence.

**Position as retrieved:**

| quantity | scoreable values | range (mL·g brain⁻¹) | spread |
|---|---:|---|---:|
| `Vu,brain`, all species | **52** | — | — |
| `Vu,brain`, **rat** | **35** | **0.436 → 4184** | **9596×** |
| `Vu,brain`, rat, **excluding cyclic peptides** | — | — | **1308×** |
| `Kp,uu,cell`, brain | **4** | — | 54× |

**The scored set is defined by a cascade, applied in this order**, each step
property-based and fixed before any residual exists:

```
all Vu,brain values                       (52, all species)
  -> E2  rat only                         (35)
  -> E1  drop values imported from another study      (oxycodone)
  -> E5  drop rows whose Vu,brain was DERIVED from a model, not measured
  -> E6  drop cyclic peptides
  => THE SCORED SET
```

**The exact retained n after every step, and the per-compound table, are written
into the frozen JSON (§7.3) before the measured column is read.** This document
fixes the *rules*; it does not pre-state a final count it has not verified.

> **The headline spread quoted in any report is the spread of the set actually
> scored — approximately 1308× on the rat, non-peptide set — and NOT the 9596×
> figure.** 9596× is the flattering number and it describes a set that includes
> compounds the model is registered as unable to handle (§3.4 E6). Quoting it would
> be advertising discriminating power the test does not have.

Even at 1308× this is **44× wider** than the 29.7× the transport law was built
against, at **n well above** the §8.1 full-verdict threshold. Every number in this
subsection is a published measurement or a count of published measurements — not a
prediction, and not an outcome of anything registered here.

### 3.2 What is under test, and what is NOT

**This dataset is `Vu,brain`.** That narrows the claim, and the narrowing is not
cosmetic:

> **B constrains the TISSUE-PARTITIONING side — how much drug a gram of tissue
> holds per unit unbound interstitial concentration. It does NOT constrain the
> intracellular partition `cell_partition` (K) in `two_route.py`.**

`Vu,brain` is a *lumped* quantity: it aggregates intracellular partitioning,
non-specific tissue binding and lysosomal trapping into one number. `K` is one
component of that lump. Scoring `Vu,brain` and reporting it as a constraint on `K`
would be the modality pooling `results/external_anchor_survey.md` §0 says this
project has *"manufactured a false regression twice"* by committing.

Test **C** (§3.8) is where `K` would be constrained — and §3.8 records that it
**cannot currently be run**.

### 3.3 GATE 1 — does the engine express this quantity at all?

**This gate must clear before any number is compared, and NOT TESTABLE is a likely
and entirely legitimate outcome.** It is registered as likely *now* so that
discovering it later cannot be presented as a surprise or worked around.

`Vu,brain` is a **volume of distribution in mL·g⁻¹**. The engine's nearest
quantities are:

| candidate | what it is | commensurable? |
|---|---|---|
| `binding.tissue_binding_rate` | a first-order **rate**, 1/s | **No — wrong dimension.** A rate is not a volume |
| `binding.plasma_free_fraction` | dimensionless, **plasma** compartment | No — wrong compartment |
| `two_route.cell_partition` (K) | dimensionless intracellular partition | a **component**, not the lump (§3.2) |
| `ion_trapping` retardation factor `R` | dimensionless tissue accumulation | closest in **role**; must be checked for definition and compartment |

**The gate.** All four must hold, each verified against the primary text and
recorded per item:

1. the engine quantity is a **ratio of tissue content to unbound interstitial
   concentration**, in the same direction;
2. its **dimensions reduce to mL·g⁻¹** under a stated, pre-registered tissue-density
   assumption, or it is dimensionless and the measured value is converted by a
   stated assumption — **the assumption is written down before the comparison, not
   after**;
3. both are **unbound-referenced** on the denominator;
4. both describe **equilibrium partitioning**, not a steady state maintained by
   active transport.

**If the only commensurable engine quantity is a rate**, then B is
**NOT TESTABLE**, with the reason *the model cannot express the measured quantity*
— the status `scorecard.py::_run_ps_scale` already carries for `PS_SCALE_S_INV`,
where the finding was that *"read as rates the drugs are mutually incoherent … read
as partitions they agree to ~1.04×. That is the signature of a units error, not a
calibration error."* **That would itself be a real and reportable finding**: that
the engine has no representation of tissue drug capacity in commensurable units.
It is not a failed exercise and must not be written up as one.

### 3.4 Exclusions, all decided on properties before any residual exists

**E1 — values imported from another study are EXCLUDED.** Oxycodone's SD-rat value
carries the footnote *"Data from Fridén et al., 2011"* — imported, not measured in
the paper supplying its neighbours. Pooling across studies imports an inter-study
methodological difference nothing in the table bounds, and
`results/external_anchor_survey.md` §0 is explicit that *"an anchor is a value PLUS
its tissue, its marker or modality, its counting/sampling convention, and its
uncertainty"*. **Fixed here, before anyone knows whether it is a hit or a miss.** At
n = 35 this is now a small fraction of the set rather than a ninth of it, so the
cost of the conservative call has fallen while the reasoning is unchanged. The
retained-with-oxycodone arm is computed and reported as a labelled sensitivity and
**is not THE result**.

**E2 — non-rat species are EXCLUDED, and species must never be pooled.** The scored
set is **rat only** (35 of the 52). Pooling species would inflate apparent n while
adding a species confound. Recorded explicitly so no later reader merges to 52. A
mouse arm may be scored **separately**, reported **separately**, and never combined.

**E3 — `Kp,uu,brain` is EXCLUDED.** A brain slice has **no intact blood–brain
barrier**, which is precisely what makes `Vu,brain` a clean tissue-partitioning
measurement. `Kp,uu,brain` is an *in vivo* measure **across** the BBB and is
confounded by P-glycoprotein efflux, which Codon does not model
(`permeability.membrane_permeability` provenance: *"passive only; ignores active
transport / efflux"*; efflux ships default-off). The acquisition agent has captured
these values in a file named `kp_uu_brain_DO_NOT_SCORE.csv` with this reasoning
written into it. **That naming is the correct engineering response** — the guard
lives next to the data rather than in a document someone may not read.

**E4 — ionisation is NOT used as an exclusion here, and the reason is registered.**
X4 in §7.2 would exclude compounds substantially ionised at pH 7.4, which would
remove much of the signal — bases are where tissue partitioning is largest. Removing
them would leave a test of nothing. The set is instead **stratified** by ion class
and the stratified reading is pre-declared in §3.7 (B-S2). This is a deliberate,
reasoned departure from a generic criterion, made before any residual is visible and
recorded rather than applied silently.

**E5 — rows whose `Vu,brain` was DERIVED FROM A MODEL are EXCLUDED. This is the most
important exclusion in the document.**

Fridén 2014's 47-row table is **not 47 measurements**. Its footnote h records that
rows carrying an `fu,brain` value had `Vu,brain` **calculated from a pH-partition
model**. The acquisition agent split them **23 measured / 24 derived** and marked the
derived rows unscoreable.

> **Scoring against a derived value would be scoring Codon against somebody else's
> pH-partition model.**

And the hazard is worse than generic circularity, because of *which* model it is: a
pH-partition model is **structurally the same kind of model** as Codon's own
ion-trapping module — `docs/VERIFICATION_METHODOLOGY.md` §8 Case B is literally about
Henderson–Hasselbalch pH partition and its saturation ceiling. Agreement between the
two would be **model-to-model concordance wearing the costume of a measurement**,
and it would be *more* likely on exactly the basic, lipophilic compounds where the
question matters most (§3.10). This is `docs/VERIFICATION_METHODOLOGY.md` §3
circularity — *"a property of the provenance graph, not of the numbers"* — occurring
**inside a primary source**, one level deeper than a derived-data hazard is normally
looked for.

**The count discrepancy is registered, not resolved.** The paper's body text says
**31/16** where its own table supports **23/24**. Both sum to 47, so the disagreement
is about **which rows are derived**, not how many rows exist — which makes it a
classification ambiguity rather than a transcription error, and means the
conservative split cannot be checked by arithmetic alone. **The conservative reading
(23 measured) is adopted.** If anyone later reconciles the discrepancy, the
reconciliation must be **documented and committed before any rescoring**, never
after seeing how the extra rows score.

**E6 — cyclic peptides are EXCLUDED, on physical grounds, decided now.**

They are the compounds that widen the set from 1308× to 9596×, so excluding them
**costs this test its most impressive headline number**. It is still the right call,
for three structural reasons — none of which is about how they score:

1. **They sit in the model's unvalidated size gap.**
   `stokes_einstein.hydrodynamic_radius_nm(mw, is_antibody)` has an empirical
   small-molecule scaling and a separate antibody branch, and §7.2 X2 already records
   that *between the two it is unvalidated*. Cyclic peptides are squarely in that gap.
2. **Their permeability is governed by a mechanism the model does not contain.**
   Cyclic-peptide membrane transport is dominated by **conformational shielding and
   intramolecular hydrogen bonding** — the molecule hides its own donors. The
   engine's membrane leg computes permeability from a Crippen logP and an **integer
   donor count**, which is precisely the descriptor that conformational shielding
   invalidates. A donor that is intramolecularly satisfied is counted as exposed.
3. **Crippen is out-of-domain for them**, in the same way §7.2 X3 records for
   square-planar Pt(II): a topological additive model has no representation of
   conformation-dependent polarity.

**The with-peptides arm is computed and reported as a labelled sensitivity**, and if
it disagrees with the headline that disagreement is reported prominently. But the
verdict is carried by the set the model is physically entitled to be scored on, and
**the headline spread quoted is 1308×, not 9596×**.

### 3.4a The measured reproducibility floor — the lower bound on any threshold

For the first time in this project, the assay's own noise is **measured rather than
assumed**, from compounds appearing in **two independent rat papers**:

| compound | value A | value B | ratio |
|---|---:|---:|---:|
| paroxetine | 714 | 876 | **1.23×** |
| diazepam | 17.8 | 20 | **1.12×** |

> **DECLARED FLOOR: 1.25×.** Rounded **up** from the worst observed ratio (1.23×),
> because rounding up is conservative for a noise floor — it declares more noise, not
> less. **Scoring Codon tighter than 1.25× on rat `Vu,brain` is scoring against
> assay noise, not against physics.**

Consequences, binding:

1. **No threshold anywhere in this document may be set below 1.25×**, now or by
   future amendment. THE TEST's 3.0× bar clears it comfortably; this is recorded so
   that any future proposal to *tighten* the bar has a floor it cannot pass.
2. **No agreement tighter than 1.25× may be described as agreement.** If the engine
   lands within 1.25× of a measured value, the defensible sentence is *"within the
   assay's own reproducibility"*, not *"accurate to 1.2×"*.
3. **The floor is reported beside every fold-error figure**, so a reader can see how
   much of any residual is assay noise. Two replicate pairs is a thin basis and the
   report must say so — this is a floor estimated from **n = 2 pairs**, and it is
   quoted as an order-of-magnitude bound, not a calibrated uncertainty.

**The mouse cross-strain overlaps are NOT used, and the acquisition agent's judgement
is upheld.** Those overlaps agree to under 2%, with a **bit-identical oxycodone
value** — and the agent judged that implausibly good for independent measurement,
concluding it more likely reflects **value reuse between papers than independent
replication**, and flagged it **unusable rather than reporting the flattering
number**. That is the correct call and it is exactly
`docs/VERIFICATION_METHODOLOGY.md` §7, *refusing your own positives*: a sub-2%
reproducibility floor would have made every subsequent comparison look far more
demanding and far more impressive, and it would have been an artefact of
transcription. **A bit-identical value across two papers is evidence of copying, not
of precision.** The 1.25× rat floor stands as the declared figure.

### 3.5 THE TEST — and why this statistic

**First, the fact that decides which statistic is even available.** With everything
in §2.3 frozen, the engine's implied tissue partition is the **same number for every
compound** — it has no per-compound tissue-capacity term (retardation was built and
**refuted**, `results/intracellular_retardation.md`; `cell_partition` is fixed at
1.0). A prediction with zero variance has **no ranks**.

> **Spearman ρ is therefore UNDEFINED here, not merely underpowered.** It must be
> reported as **NOT COMPUTABLE**, never as 0.00, and never as a failure. This is a
> property of the model, not of the sample size, and n = 35 does not change it.

**Second, why a bare fold-error bar would be the wrong test.** Across a 1308-fold
measured spread, *any* constant predictor carries a large unsigned geometric-mean
fold error **by construction** — the spread alone forces it, however well-centred the
constant is. A bar stated directly on `GMFE(engine)` would be near-unreachable
regardless of centring, and would conflate two different failures: *the assumed value
is in the wrong place* and *the model has no per-compound variation*. The second is
already known and is registered in the paragraph above. A bar that cannot
distinguish them measures nothing.

**The test is therefore baseline-relative — the house's own pattern.**

> **THE TEST (B).** Let
> ```
> GMFE(p)      = exp( mean_i | ln( p / Vu_i ) | )      over the retained compounds
> GMFE_floor   = min over all constants c of GMFE(c)    (the best any constant could do)
> THE TEST     : GMFE(engine) / GMFE_floor  <=  3.0
> ```
> The engine's assumed value must be **within 3×, in geometric-mean-fold-error
> terms, of the best constant anyone could have chosen** on this measurement.

**What this does.** It quotients out the irreducible cost of having no per-compound
term — which is `GMFE_floor`, and which is arithmetic about the measurement rather
than a finding about Codon — and isolates the only question the data can answer about
a constant: **is it in the right place?** It is the same structure as
`scorecard.py::_run_gist_imatinib`'s *"beat the no-partitioning null by >= 3x GM fold
error"*, with the oracle constant serving as the floor rather than the rival. **3.0**
is the repository's bar (`BAND_TIGHT`, `ORDER_OF_MAGNITUDE_FACTOR`,
`delivery_metrics.md` §3b), adopted unchanged, and it has already produced a
published FAIL. It also clears the §3.4a reproducibility floor of 1.25× comfortably,
so the bar is not scoring against assay noise.

**Registered in advance: THE TEST is FORGIVING, and a PASS is narrow.** Because
`GMFE_floor` grows with the measured spread, dividing by it on a 1308-fold set leaves
a permissive bar. **A PASS means the engine's constant is not grossly miscentred. It
does not mean the engine is accurate, and it must never be reported as though it
did.** The centring ratio (B-S4) is reported precisely so a reader can see how much
room the bar leaves. This is stated here, before the run, rather than discovered as a
caveat afterwards.

### 3.5a `GMFE_floor` is a headline finding whatever the verdict

> **`GMFE_floor` is reported in the results document's opening summary, in the
> scorecard row's `note`, and in every abstract or summary of this work — under a
> PASS, a FAIL, and a NOT TESTABLE alike. It may not be omitted from any summary,
> and it may not be relegated to a methods section.**

The reason it must be protected this way: **since the engine's `Vu,brain` prediction
is constant, `GMFE_floor` is the entire story.** It is the best error any constant
could achieve on real measured tissue partitioning — which is exactly **the price the
engine pays for having no per-compound tissue-capacity term**, quantified against
data for the first time. It answers the question a user actually has: *what would a
per-compound capacity term be worth?*

A PASS on THE TEST would say only that the constant is well-placed **within** that
floor. Reporting the PASS without the floor would let a permissive bar (§3.5) stand
in for an accuracy claim the test never made — the precise failure mode this project
has retracted results for. **A PASS is the outcome most likely to bury this number,
which is why the rule is written before the result exists.**

Two further properties, both registered now:

* **`GMFE_floor` is a property of the measurement alone.** It does not depend on any
  engine constant, any descriptor, or any code version. It therefore **does not move
  if a bug fix lands between registration and execution** (§9.5) — the headline
  finding of test B is immune to the `HD` units defect and to any recalibration of
  `alpha`.
* **It is reported with a bootstrap interval** over compounds, on the same resamples
  as §3.6, so its own precision is visible.

### 3.6 Power — and whether the §8.1 override is still needed

**Asked and answered, in both directions, as required.**

**On power: the override is no longer needed, and I say so plainly.** At **n = 35**
the scored set clears §8.1's `n ≥ 20` full-verdict tier outright. The bootstrap-CI
rule was introduced when the set was n = 8 and the generic ladder would have returned
NOT TESTABLE; that justification has now lapsed. Recording the lapse rather than
quietly keeping a rule whose stated reason expired.

**On sufficiency: the bootstrap-CI rule is RETAINED, for a different and still-live
reason.** §8.1's ladder is calibrated by *minimum detectable ρ* — a quantity that is
**undefined for a constant predictor** (§3.5). The ladder can tell us n is adequate;
it cannot tell us whether a *ratio of two geometric-mean fold errors* is resolved,
because it was never designed to. So:

> **B's verdict stands only if the bootstrap 90 % confidence interval on
> `GMFE(engine) / GMFE_floor` lies ENTIRELY on one side of 3.0.** If the interval
> straddles 3.0, the verdict is **UNDERPOWERED → NOT TESTABLE**, whichever side the
> point estimate falls on.

This is `results/ranking_robustness.md`'s own convention — *"probability whose 90 %
credible interval lies entirely above 0.90, and the runner-up's lies entirely below
0.10"* — reused unchanged. **The distinction being drawn is between sample size and
statistical resolution**: n = 35 makes a straddle much less likely than n = 8 did,
but "unlikely" is not "impossible", and the rule costs nothing and can still return
NOT TESTABLE.

Bootstrap over compounds, BCa, ≥ 10 000 resamples, seed recorded. Point estimate and
interval both printed, for THE TEST ratio and for `GMFE_floor` separately.

### 3.7 Reported alongside — ranked, and none may replace THE TEST

| # | measure | what it shows |
|---|---|---|
| **B-S1** | **`GMFE_floor`** and the constant achieving it, with bootstrap interval | **the headline finding (§3.5a)** — the price of having no per-compound tissue term, measured |
| **B-S2** | **centring ratio**: engine's constant ÷ the `GMFE`-minimising constant | how much room THE TEST's forgiving bar (§3.5) actually leaves. **Diagnostic, not an acceptance criterion** |
| **B-S3** | **stratified by ion class**: bases vs non-bases | §3.4 E4. The engine has no ionisation term in this path, so a miss concentrated in bases is the **expected direction**. A stratum carries a verdict only if it independently clears §8.1; otherwise diagnostic. **Strata may never be used to select a subset for the headline** |
| **B-S4** | **direction** of the bias — systematically low or high? | a one-directional offset is a correctable error; a symmetric spread is an irreducible one |
| **B-S5** | within-3× and within-10× coverage | the repo's two declared bands |
| **B-S6** | worst-case fold error, compound named | stops a good mean hiding a catastrophic miss |
| **B-S7** | **the 1.25× reproducibility floor**, printed beside every fold-error figure | §3.4a. Lets a reader see how much of any residual is assay noise |
| **B-S8** | Spearman ρ **only if GATE 1 yields a varying prediction** | otherwise **NOT COMPUTABLE** (§3.5) |
| **B-S9** | Spearman ρ between **measured `Vu,brain`** and **Crippen logP** | **the compensating-error probe** of §1.2(i). Two measured/computed quantities with no Codon prediction in the loop — **hypothesis-generating only, never a result about the engine**, and now further weakened by §6's descriptor finding |
| **B-S10** | the with-oxycodone, with-cyclic-peptide, and with-derived-rows arms | the sensitivity arms of §3.4. Reported; **never** the verdict |

**No secondary may be promoted to THE TEST after the run.** Reporting several
statistics without a pre-named primary is equivalent to no pre-registration at all.

### 3.8 TEST C — `Kp,uu,cell`: registered, and NOT RUNNABLE

**Position:** only **4 scoreable brain `Kp,uu,cell` values** were obtained, spanning
54×.

> **Four points cannot discriminate a partition coefficient, and test C is therefore
> NOT RUNNABLE. `cell_partition = 1.0` — a hardcoded, load-bearing input to the
> transcellular route — REMAINS UNCONSTRAINED BY ANY MEASUREMENT. That is the honest
> status and it is to be reported as such wherever the two-route law is described.**

This is deliberately recorded as an **open limitation rather than a gap in a to-do
list**. `docs/VERIFICATION_METHODOLOGY.md` §7 distinguishes FALSIFIED / PARTIALLY
SUPPORTED / **UNDERPOWERED** / NOT TESTABLE, and warns that *"collapsing
'underpowered' into either 'supported' or 'refuted' is the most common way an honest
analysis becomes a dishonest one."* Four points is underpowered, and scoring them
anyway would manufacture a verdict.

**Registered so it cannot drift:**

1. **C is not merged into B, ever** — B is a lumped tissue quantity, C an
   intracellular partition (§3.2). Pooling them would restate exactly the error §3.2
   exists to prevent.
2. **A PASS on B does not constrain `K`**, and the §3.9 sentence says so explicitly.
3. **C becomes runnable at n ≥ 20** scoreable values meeting §3.3's commensurability
   gate — the §8.1 threshold, which does govern C because C's prediction is also
   constant and the same reasoning as §3.6 applies. Below that it stays registered
   and unrun.
4. **If C is later run, THE TEST for C is the same baseline-relative form** as §3.5 —
   the shipping constant within **3.0×** of `GMFE_floor` — with the same bootstrap
   rule, the same headline treatment of `GMFE_floor`, and its own row and verdict.
5. **The 4 values are not quoted as evidence in either direction**, including as
   reassurance that `K = 1.0` is "about right".

### 3.9 What a PASS on B licenses — the exact sentence

This project's recurring failure is a real result restated one notch broader than it
earned. The following is the **only** claim a PASS supports, and it is to be quoted
verbatim rather than paraphrased:

> **On a PASS:** *"Across N structurally diverse compounds whose unbound brain
> tissue distribution volume was measured in rat brain slices by a single shared
> method, the tissue partition the engine assumes lies within 3× — in
> geometric-mean-fold-error terms — of the best single constant obtainable on that
> measurement. The best any constant could achieve is GMFE_floor = F, which is the
> error the engine incurs for having no per-compound tissue-capacity term at all.
> This constrains one input coefficient of the transport PDE. The transport solver
> was not invoked, no concentration field was computed, no tumour tissue was
> involved, no per-compound variation was predicted or tested, and nothing here is
> evidence that the engine's tissue predictions are validated."*

**`GMFE_floor` is inside the PASS sentence deliberately** (§3.5a), so that the
number a PASS is most likely to bury cannot be dropped when the sentence is quoted.

**Not claimable on a PASS:** that the transport model is validated; that `D_eff` is
accurate; that the engine is accurate on tissue partitioning (§3.5 — the bar is
forgiving); that the two-route law generalises; that `cell_partition` is constrained
(§3.2, §3.8 — it remains unconstrained); that anything about human tumour follows
(§3.10).

### 3.10 What a PASS licenses about human tumour — less than it will feel like

My honest view, written before the result:

> **A PASS licenses only that the assumed tissue partition is in a plausible range
> for mammalian tissue generally. It licenses nothing about tumour.**

1. **Different species** — rat.
2. **Different organ** — brain is atypically lipid-rich, and a partition is exactly
   the quantity lipid content moves.
3. **Different pH environment, and this is the sharp one.** Tumour interstitium is
   acidic. For a weak base, tissue partitioning is a pH-partition quantity, so these
   values are taken under a pH gradient that is **not** the tumour one — and bases
   dominate both the set and its upper tail. **The anchor is weakest exactly where
   the quantity matters most**, and that sentence should be quoted whenever this
   anchor is cited. The B-S3 ion-class stratification exists to make the size of
   this exposure visible rather than arguable.
4. **The engine's ion-trapping module has a documented structural ceiling** —
   `docs/VERIFICATION_METHODOLOGY.md` §8 Case B records that the
   Henderson–Hasselbalch trapping ratio saturates for any sufficiently strong base,
   which is a resolution limit no parameterisation removes. A brain-slice partition
   does not constrain it.

---

## 4. Outcome buckets and withdrawal — test B

### 4.1 Buckets

**B-A — CONSTRAINED.** GATE 1 clears; THE TEST passes with the bootstrap interval
entirely below 3.0. *Claimable:* §3.9's sentence, verbatim. Nothing more.

**B-B — MISCALIBRATED.** GATE 1 clears; the interval lies entirely **above** 3.0.
*Reading:* the engine's assumed tissue partition is in the wrong place by more than
the bar, and every transport number inherits it. **FAIL** — it must appear in
`codon validate` and make it exit 1. Required next action is to constrain the
coefficient, **not** to widen the bar. B-S4 decides the shape of the fix: a
one-directional offset is a single corrected value; a symmetric spread needs a
per-compound term and a full blast-radius re-run.

**B-C — UNDERPOWERED.** The interval straddles 3.0. **NOT TESTABLE** (§3.6). Not a
soft pass; never pooled with one.

**B-D — NOT COMMENSURABLE.** GATE 1 fails. **NOT TESTABLE**, reason *the model
cannot express the measured quantity* — and §3.3 registers that this is a real and
reportable finding rather than a failed exercise.

### 4.2 What would make us abandon the claim entirely — test B

**The claim B can kill:**

> The engine's treatment of tissue drug capacity is a physically defensible
> approximation, rather than a placeholder whose value was never checked against
> anything.

**Any ONE of the following withdraws it — not narrows it.**

**BW1 — the assumed value is outside the measured range entirely.** The engine's
implied `Vu,brain` lies outside `[min, max]` of the retained measurements. The
engine's tissue capacity would then be set to a value **no compound in a
structurally diverse set of 35 exhibits**. This is a reachability argument
(`docs/VERIFICATION_METHODOLOGY.md` §8), which produced a definitive negative on
n = 6 in `results/gist_imatinib.md`, and it is decisive because it is not a
statistical claim at all — one comparison, no resampling.

*Registered as a property of the data, not an outcome:* the rat range spans
**0.436 → 4184 mL·g⁻¹**, so it straddles unity. Any engine constant inside that
interval clears BW1, which makes BW1 a genuine but **not easily tripped** condition
on this set — recorded so that clearing it is not later mistaken for a result.

**BW2 — the tissue capacity term is a proxy that the membrane leg has absorbed.**
B-S9 shows measured `Vu,brain` correlates with Crippen logP at |ρ| ≥ 0.7
(permutation p < 0.05, with the minimum detectable value printed),
**and** supplying a logP-dependent tissue partition materially degrades the HT-29
re-score relative to the shipping constant. That combination is the signature of
**compensating errors** — the membrane term absorbing tissue capacity's
compound-dependence, so that supplying the true dependence breaks the cancellation.
`docs/VERIFICATION_METHODOLOGY.md` §1: *"When a correct fix makes the output look
worse, that is information, not a regression."* Under BW2 the mechanism claim is
withdrawn **even though the HT-29 fit is unchanged — the fit was never the
evidence.**

*(The HT-29 re-score inside BW2 is legitimate despite §0: it is a
mechanism-attribution ablation on spent data, reported under that heading. It can
only ever **withdraw** a claim, never support one.)*

**BW3 — the model cannot express the measured quantity, and the reason is a
units/role error rather than an absence.** GATE 1 reveals that a quantity the engine
*consumes as* a partition is *defined as* a rate, or vice versa — the
`_run_ps_scale` signature. That is a **units/convention defect**: the parameter must
be withdrawn and re-derived, **not re-fitted**, because as
`docs/VERIFICATION_METHODOLOGY.md` §10b records, a consistently-applied wrong
convention validates itself.

### 4.3 What does NOT trigger withdrawal

* **Failure to predict per-compound variation** — arithmetic, registered in §3.5,
  and not evidence about anything.
* **A miss concentrated in the base stratum** — B-S3; the engine has no ionisation
  term in this path and that is a **known**, already-documented limitation, not a
  new finding. It bounds what the anchor constrains; it does not withdraw a claim.
* **A `Kp,uu,cell` mismatch** — not in this dataset (§3.2); test C's business.
* **`Kp,uu,brain` disagreement** — excluded by §3.4 E3; must not be computed.
* **A mouse-column disagreement** — excluded by §3.4 E2; scored separately if at all.
* **A single outlier compound** — reported as B-S6, investigated, **not removed**
  (§7.1). At n = 35 a single removal is easier to rationalise as harmless, which
  makes the prohibition more important rather than less.
* **A large `GMFE_floor`** — that is the measured price of the missing per-compound
  term (§3.5a), and it is a finding, not a failure of the engine to clear a bar.
* **Too few compounds** — bucket B-C. Absence of a test is not a falsification.

### 4.4 Is B falsifiable? — the audit

| condition | reachable? | why |
|---|---|---|
| BW1 | **yes** | the engine's implied value either lies inside `[0.436, 4184]` or it does not; one comparison, no statistics |
| BW2 | **yes** | both halves computable in the same run; the HT-29 re-score re-runs existing code |
| BW3 | **yes** | the precedent exists in this repo (`ps_scale`) and produced a NOT TESTABLE |
| THE TEST failing | **yes** | `GMFE_floor` is a property of the measurement alone; the engine's ratio to it is unconstrained above |

**B is falsifiable.** Every condition computes from the same run that would produce a
positive result, at no extra cost, and none requires a judgement call.

**The honest weak point is no longer n.** At 35 rat values the sample clears §8.1
outright. The two remaining soft spots are named instead: **GATE 1** may find no
commensurable engine quantity at all (§3.3), and **THE TEST is forgiving by
construction** on a wide-spread set (§3.5) — which is why a PASS is bounded by the
§3.9 sentence and why `GMFE_floor` is protected as a headline (§3.5a).

---

## 5. TEST A — a second cell line *(BLOCKED, registered anyway)*

### 5.1 Why a blocked test is still worth registering

The acquisition agent has the primary source — both PDFs, 71/71 rows extracted and
cross-checked field-by-field against the repo's existing literal with **zero
discrepancies** — and the finding is definitive:

> **Pruijn 2008 does not publish per-compound `Dmcl` for SiHa, HCT8-Ea or
> HCT8-Ra.** The nine-compound cassette was run in all four lines — Table 1's
> per-line fits could not exist otherwise — but only the HT-29 arm is tabulated.
> Not in Table S1 (caption: *"for HT29 MCLs"*), not as raw points in any figure.
> Fig. 4 is HT-29 only; Fig. 5 is *simulations*, not data. **There is nothing to
> digitize.**

What exists per line is urea `Dmcl` in four lines (one compound, four points, 20.5×
spread) and two *fitted* coefficients with the rest held at HT-29 values. **The
per-line MCL thicknesses are back-computed from the urea `Dmcl`, so using them to
predict urea is circular** — that circularity is registered here as a forbidden
input.

**Registering it now is the point.** Fixing the criterion while the data is
unobtainable is the strongest possible pre-registration: it cannot be shaped by the
data because the data does not exist to shape it. If the cassette arms ever arrive,
the analysis is already committed.

**Unblocking condition:** the corresponding author (`f.pruijn@auckland.ac.nz`) is
the only route to the per-line cassette arms; the Accessory Publication states the
compound set *"is available from the author"*. §1.3 argues that email is worth
sending regardless of the odds.

### 5.2 The criterion, fixed now

If per-compound `Dmcl` for a second line is obtained:

> **THE TEST (A).** On the held-out line, with §2.3 frozen, arm **M**
> (`build_profile(mol)`) must beat arm **B** (`build_profile(mol, two_route=False)`
> — the falsified MW-only predecessor, bitwise recoverable and pinned by
> `tests/test_two_route.py::TestShippingDefault::test_the_superseded_value_is_still_exactly_recoverable`)
> on unsigned geometric-mean fold error by **≥ 3.0×**.
>
> ```
> GMFE(arm) = exp( mean_i | ln( D_pred,i / D_meas,i ) | )
> THE TEST  :  GMFE(B) / GMFE(M)  >=  3.0
> ```

**Unsigned**, because `docs/VERIFICATION_METHODOLOGY.md` §7c: *"Never accept a
signed/symmetric error metric as an acceptance criterion"* — signed measures let
errors cancel and therefore defend compensating errors. **3.0** because it is the
repository's own bar and **has already produced a published FAIL**
(`scorecard.py::_run_gist_imatinib`). A third arm **N** — every compound at the
measured geometric mean — is an oracle: **beating N is weak evidence, failing to
beat N is decisive.**

Secondary measures: within-3× coverage; Spearman ρ with permutation p and the
minimum detectable ρ; modelled ÷ measured spread; signed geometric bias
(**diagnostic only, never an acceptance criterion**); matched-mass fold-separation
and pairs ordered correctly; worst-case fold error with the compound named.

### 5.2a Which descriptors test A is scored on — settled, not left open

`results/descriptor_falsification.md` has landed and its recommendation is adopted
**unconditionally**:

> **THE TEST (A) is scored on the paper's own measured descriptor columns
> (`logP7.4`, `HD`), exactly as the original scoring did.** The shipping-descriptor
> arm is computed and reported **separately**, and never carries the verdict.

**Why.** A is a test of whether the *physics* transfers to a tissue it was not
selected on. Feeding it descriptors that are themselves a known, large error source
would confound the two, and a poor result would be uninterpretable — it could equally
mean the transport law does not transfer or that Crippen mis-estimated a logP. Test A
answers one question; the descriptor pipeline is a separate question with its own
document. This is the same separation `results/two_route_promotion.md` §5 draws when
it insists the physics provenance be named on every row.

**The consequence, stated plainly:** A's PASS would be a statement about the two-route
**law**, not about the **shipping pipeline**. Both must be named in any report, and
the shipping-pipeline number is the one a user actually gets (§6).

### 5.2b Pre-declared R² regimes — so a descriptor loss is not misread as a physics failure

`results/descriptor_falsification.md` supplies an interpretation key, and it is
encoded here **before** any A result exists, because without it any descriptor-fed R²
in the 0.2–0.5 band would be argued about after the fact:

| observed R² on a descriptor-fed arm | pre-declared reading |
|---|---|
| **≥ 0.45** | **Fully consistent with correct physics plus σ ≈ 1 descriptor noise.** Does **NOT** implicate the transport law. May not be reported as a physics failure |
| **0.20 – 0.45** | **AMBIGUOUS.** Descriptor noise and physics error are not separable at this level. The verdict must come from the measured-descriptor arm (§5.2a); the descriptor-fed arm is reported without a physics interpretation |
| **< 0.20** | **The first regime that implicates the physics itself.** Only here may a descriptor-fed arm be cited as evidence against the transport law — and even then the measured-descriptor arm decides |

**This band is binding in both directions.** It forbids reading a mid-range R² as
vindication of the physics just as firmly as it forbids reading one as refutation.

**The caveat that travels with it, and it is load-bearing.** The projection behind
these bands rests on a **five-compound** descriptor-error estimate, because **66 of
the 71 Pruijn compounds are ACSRC `SN` codes with no resolvable structure anywhere**.
`descriptor_falsification.md` reports that weakness as one of its own findings, and it
propagates here: **the regime boundaries above are order-of-magnitude guides, not
calibrated thresholds**, and any report using them must say so. They are precise
enough to stop a category error and not precise enough to settle a close call.

### 5.3 "Constants unchanged", operationalised

> **Geometry may be measured. The answer may not be.**

**Admissible per-line inputs** — `cell_volume_fraction` and `cell_radius_um` only,
and only if all four hold: sourced from a measurement **independent of the `Dmcl`
values**; not derived from or fitted to any diffusion measurement; written into the
frozen JSON before the measured column is read; reported per line with its source.
**If any fails, both stay frozen at their shipped values in the headline arm**; a
per-line arm may be reported as **exploratory** and may not carry the verdict.

**Forbidden — supplying any of these voids test A:** the line's asymptote
coefficients (that is the answer in the paper's parameterisation); its QSAR fit or
any coefficient of it; any per-line scale factor chosen to reduce residuals; any
`Dmcl` value from the held-out line in any aggregate form; **the per-line MCL
thickness, when urea is among the scored compounds** (§5.1 — it is back-computed
from urea's own `Dmcl`); and selecting which line to score after seeing performance.

**Pre-flight, required before any scoring**, following the pattern already in
`experiments/diffusivity_falsification.py`:

```
PRE-FLIGHT A: is the shipping compile INVARIANT to cell-line identity?
  -> the two D_eff must be BITWISE IDENTICAL
PRE-FLIGHT B: if per-line geometry is admitted under 5.3, does it reach D_eff?
  -> the two D_eff must DIFFER
```

**Pre-flight A is expected to pass, and that is the substance of §5.4.**

### 5.4 The between-line level: arithmetic, not evidence

With shipping defaults `build_profile` takes a `Molecule` and nothing else — no
cell line, no tissue, no geometry. The prediction for a compound is the **same
number for every line**. The model is **flat by construction between lines**.

> **A between-line level miss, on its own, is NOT new evidence against the law, and
> closing it by supplying per-line φ_cell is NOT new evidence for it.** The
> informative question is not *did it miss* but *could it have hit*.

**The reachability procedure, which IS informative:** sweep φ_cell across its
published admissible range **0.40–0.85** (`results/cellular_exclusion.md` §4.3)
**under the shipping two-route law** — *the published sweep in that file is the
impermeable-cell arm B, not the shipping permeable-cell law, so the bound must be
recomputed and may not be quoted from that table* — and obtain `R_max`, the largest
between-line `D_eff` ratio cell packing alone can produce.

| outcome | reading |
|---|---|
| measured movement **inside** the reachable set | packing is a *sufficient* account in principle; whether the *right* φ_cell does it is a separate, weaker question and must be labelled so |
| measured movement **outside** at every admissible φ_cell | **no packing reproduces it.** A mechanism is missing — withdrawal condition **AW5** |

**And the quantity must match.** The 22× figure of §0.4 is a spread in *fitted
asymptote coefficients*, which is not a per-compound `Dmcl` spread. If the only
available between-line quantity is a coefficient the engine does not compute, that
comparison is **NOT TESTABLE**, not approximated.

### 5.5 Buckets and withdrawal — test A

**A-A — TRANSFERS, NARROWLY.** THE TEST passes, within-3× ≥ 2/3, ρ > 0 at
permutation p < 0.05. *Claimable:* the law's advantage over molecular weight
transfers to a line it was not selected on. Nothing more.

**A-B — ORDERING RIGHT, MAGNITUDE WRONG.** *The most likely outcome, pre-classified
so it cannot be spun.* ρ ≥ +0.5 at p < 0.05 but THE TEST fails or coverage < 2/3.
*Reportable as* **PARTIAL, and the claim narrows**: it becomes a **ranking** claim;
compiled `D_eff` **may not be quoted as a number** for an unseen tissue (every
surface that emits one — `results/bring_your_own_compound.md`, `tmesim/molecular/byo.py`
— must carry that); it is recorded **FAIL on THE TEST** with ρ in the `note` field,
not a pass on a substitute metric; and it is reported with the symmetry the
methodology demands — the oxygen PDE tops its baseline table on ρ *and* has a
diffusion length 2–5× the published value, and *"a reader who takes only the first
has mislearned the chapter."* **Not success, not withdrawal.**

**A-C — LEVEL RIGHT, ORDERING ABSENT.** THE TEST passes but ρ ≈ 0, or M fails to
beat N. The cell term is acting as a **constant rescale**; mechanism claim
**FAIL**.

**A-D / A-E — fails outright / NOT TESTABLE.**

**Withdrawal conditions for A. Any ONE withdraws the claim** that the law tracks
measured diffusivity better than MW in a tissue it was not selected on:

* **AW1** — `GMFE(B)/GMFE(M) < 1.0`: the falsified predecessor is not beaten.
  *(Registered asymmetry: 1.0–3.0 is a failure of THE TEST and a narrowing; below
  1.0 is withdrawal.)*
* **AW2** — `GMFE(N) ≤ GMFE(M)`: an oracle constant beats the model.
* **AW3** — ρ ≤ −0.3 at permutation p < 0.05: the ordering **inverts**. A user is
  systematically misled rather than merely uninformed. Additionally triggers a
  formal proposal to flip `_TWO_ROUTE_DEFAULT` back to `False`, since ordering is
  what the promotion was justified on.
* **AW4** — **the sharpest.** At matched molecular weight the model's predicted
  fold-separation is ≤ **1.05×** (flat — the MW-only signature) while the
  measurement separates by ≥ **2×** in the same mass groups. This is the exact
  diagnostic that killed the predecessor, which scored **0 of 31** on it. It is a
  *within-group* comparison, so it survives small n and does not depend on absolute
  level. If the membrane leg discriminates on HT-29 and is flat elsewhere, the
  discrimination was a property of the fitting set.
* **AW5** — structural unreachability at every admissible φ_cell (§5.4).

**Does NOT trigger withdrawal:** a between-line level miss (§5.4); under-separation
at the registered magnitude (§6.3); a further absolute `D_eff` miss (already
conceded — `two_route_promotion.md` §8.3 records tirapazamine at 1.50× low, having
been 2.21× high); a single outlier; too few compounds.

---

## 6. The descriptor confound — it sits upstream of both tests

### 6.1 The finding

`experiments/two_route_transport.py:246-247`, verbatim:

> *"logP7.4 and HD are the paper's OWN descriptor columns -- inputs, never the
> scored quantity."*

Every published Pruijn score — R² = 0.718, Spearman +0.821, geometric bias 1.03×,
within-3× 99% — was obtained by feeding the law **Pruijn's measured logP₇.₄ and
donor columns**. The engine at §2.2's entry point computes its own: **RDKit Crippen
`MolLogP`** and an RDKit donor count.

> **The R² = 0.718 headline is a score for the two-route *law* supplied with
> measured descriptors. It is not a score for the *pipeline* a user invokes via
> `build_profile`, and that pipeline has never been scored against any measured
> diffusion coefficient.**

Because the descriptor layer was never *selected* against these data, measuring it
on HT-29 is legitimate and is **not** spent (§0.3).

**That measurement has now landed** — `results/descriptor_falsification.md`. Its
projection for the shipping descriptor pipeline is **R² ≈ 0.22–0.50, central ≈ 0.27**,
against the published 0.718 — **roughly a 62 % loss**. Two things follow, and both
are already binding above:

1. **Test A is scored on the paper's measured descriptor columns** (§5.2a), so the
   physics question and the descriptor question stay separate.
2. **A mid-range R² on a descriptor-fed arm may not be read as a physics failure**
   (§5.2b), because ≈ 0.45 is fully consistent with correct physics plus σ ≈ 1
   descriptors.

**And the estimate's own weakness is part of its content.** It rests on **five
compounds**, because **66 of the 71 Pruijn compounds are ACSRC `SN` codes with no
resolvable structure anywhere**. That is why §5.2b's boundaries are declared as
order-of-magnitude guides rather than calibrated thresholds. A five-point error model
extrapolated across 71 compounds is a real result and a thin one simultaneously, and
this document treats it as both.

### 6.2 Why it is not a footnote

`results/logp_provenance.md` §1: *"Every logP tmesim uses is an RDKit Crippen
`MolLogP` prediction computed from SMILES at runtime."* Correcting doxorubicin's
logP from Crippen 0.001 to the measured 1.27 moves `alpha` by **0.462×**. Crippen
computes a **neutral octanol/water logP** while Pruijn's column is a **measured
logP₇.₄**, a distribution coefficient at pH 7.4 — `two_route_transport.py:492-495`
prints exactly this warning. Crippen is **stereochemistry-blind** and has no
parameters for coordination complexes (cisplatin +0.195 against a measured −2.53).

**And the existing bound on this hazard is stale.** `logp_provenance.md` §1 finding
2 states *"logP cannot move transport directly … logP is not an input to it at any
point."* True on 2026-08-08; the two-route promotion landed 2026-08-14 and routes
`logp` and `hbd` **directly into `effective_diffusivity`** (`compile_drug.py` step
3c). The file carries **no correction notice and no mention of `two_route`** —
verified by grep, not assumed. This is `docs/VERIFICATION_METHODOLOGY.md` §10c
recurring: a results file gone stale against a changing codebase. **The blast
radius of a wrong logP on `D_eff` is not bounded by any existing document.**

### 6.3 How the sibling's result changes each test

| sibling finds | effect on **B** | effect on **A** |
|---|---|---|
| **descriptors reproduce R² ≈ 0.718** | THE TEST unaffected. **B-S7 and withdrawal condition BW2 become more interpretable**, since a `Vu,brain`-vs-logP correlation then reflects real chemistry rather than descriptor noise | A's arms are clean; THE TEST is scored on the shipping descriptors and the result stands on its own |
| **descriptors degrade the score materially** | **THE TEST is still unaffected** — it compares a constant against a measured distribution and consumes no descriptor. But **BW2 is weakened to hypothesis-generating only**, because a correlation with a *noisy* logP cannot establish that the membrane leg absorbed tissue capacity. BW2 may then be *raised* but not *acted on* | **A is confounded.** A must then be scored in a 2×2 — {measured, Crippen} × {HT-29, held-out line} — with THE TEST on the *measured-descriptor* cell to isolate cell-line transfer, and the Crippen/held-out cell reported separately as the honest user-facing number |

**The attribution rule, pre-registered:**

```
descriptor-attributable degradation  =  (HT-29, measured)  ->  (HT-29, Crippen)
line-attributable degradation        =  (HT-29, measured)  ->  (held-out, measured)
```

**If the descriptor-attributable degradation exceeds the line-attributable
degradation, the shipping-pipeline result is confounded by its own descriptor
layer**, is reported **NOT TESTABLE** as a statement about cell-line transfer, and
the required next action is descriptor sourcing — not a change to the transport
law. Conversely, if the HT-29/Crippen cell degrades while the held-out/measured
cell does not, the finding is about `descriptors.py` and `two_route.py` must be
explicitly exonerated.

**A sixth withdrawal condition, applying to the pipeline rather than the law:**

* **W6** — the discrimination does not survive the shipping descriptors and
  attribution assigns the loss to descriptors. This does **not** withdraw the claim
  about the *law*, but it **withdraws the claim about the shipping pipeline**, and
  the two must be separated in every surface quoting R² = 0.718 to a user — because
  that number was never produced by the code path the user runs.

### 6.4 The registered prior expectation on under-separation

On HT-29 the shipping model reproduces 8.3× against a measured 29.7× — about
**24%**, an under-separation of roughly **4×** — and `two_route_promotion.md` §8.1
records that this was *predicted before the scoring run*. **Registered now:** we
**expect** under-separation of broadly this magnitude again. Therefore
under-separation near ~4× is **not** evidence against the model **and not** evidence
for it; only a materially larger shortfall is informative, and AW4 states where it
becomes so. Writing this down in advance is what stops the same number reading as
"consistent with expectation" in a good result and "a known limitation" in a bad
one.

---

### 6.5 The `HD` convention mismatch — a second, DIRECTIONAL descriptor defect

Distinct from §6.1–6.4, and not the same kind of problem. That one is **noise**; this
one is **systematic and signed**.

**The defect, as verified by the coordinator:** Pruijn's `HD` column is donor
**hydrogens** (`NHOHCount`). The engine computes `Lipinski.NumHDonors`, donor **heavy
atoms**. They agree for hydroxyls and diverge for NH₂ and primary amines — urea is
**paper 4, engine 2**.

> **The two-route law was SCORED against one quantity and is FED another at runtime.**

**Why this is not simply "the engine is wrong".** `HBD_LOGP_DECREMENT = 1.34` was
derived as `3.621 × 0.37`, where 0.37 is described as the acidity of a *typical donor
group* — language that points at the **heavy-atom** convention — while the R² = 0.718
that validated it was computed under the **hydrogen** convention. So the constant may
be defined for one convention and was validated under the other.

That is precisely the pattern `docs/VERIFICATION_METHODOLOGY.md` §10b records for the
α discrepancy of exactly `1/ln²2`, where *"both sides of the comparison used the wrong
convention, which is why it passed"*, and which it generalises as:

> **"A consistently-applied wrong convention validates itself."**

**This document takes no position on which convention is correct.** That belongs to
the running investigation, which owns it; nothing has been changed. What is registered
here is what the answer does to these tests.

**Consequences, registered in advance:**

| test | exposure |
|---|---|
| **B** | **IMMUNE.** B's prediction is a constant that consumes no donor count, and `GMFE_floor` is a property of the measurement alone (§3.5a). Neither THE TEST nor the headline finding moves under any resolution of this defect |
| **A** | **EXPOSED, and it changes how A must be run.** §5.2a specifies A's primary arm on the paper's measured columns — **including `HD`**. If the constant is defined for heavy-atom donors, that arm feeds the wrong convention into it |

**The rule for A, fixed now:**

1. If the investigation establishes which convention `HBD_LOGP_DECREMENT` is defined
   for, A's primary arm uses **that** convention, converting the paper's column if
   needed, with the conversion stated.
2. If it cannot be established before A runs, **A is run under BOTH conventions and
   both are reported.** If they agree on the bucket (§5.5), the verdict stands and the
   convention question is noted as immaterial to it. **If they disagree, A's verdict
   is NOT TESTABLE pending the investigation** — not the more favourable of the two.
3. The propagation route is named so its size is not guessed at: the decrement is
   applied **per donor**, so a miscount of `k` donors shifts the membrane partition by
   `1.34 × k` log units. **The magnitude of the resulting blast radius is for the
   investigation to measure, not for this document to assert.**

**It may move `alpha` and the Pruijn R².** §9.5 governs what happens if it does.

---

## 7. Exclusions and anti-tamper *(both tests)*

### 7.1 The rule that makes the rest meaningful

> **Exclusion is decided from compound properties and source metadata ALONE, and
> frozen before any predicted value is joined to any measured value. No compound
> may be excluded because of how it scores. No criterion may be added, widened or
> reinterpreted after any residual is visible.**

`docs/VERIFICATION_METHODOLOGY.md` §6: goalpost-moving is *"invisible after the fact
— including to the person doing it."* Dropping outliers after seeing residuals is
the most efficient way to make either test worthless and would be indistinguishable
from honest curation in the written record.

### 7.2 Criteria

| # | criterion | model limitation forcing it |
|---|---|---|
| **X1** | no authoritative structure resolves to a SMILES | every descriptor is computed from SMILES |
| **X2** | MW above the small-molecule ceiling declared in §8.3 | `stokes_einstein.hydrodynamic_radius_nm` uses an empirical small-molecule scaling with a separate antibody branch; between them it is unvalidated |
| **X3** | metal coordination complex (e.g. square-planar Pt(II)) | **verified**: Crippen *"is a topological additive model and has no parameters for a square-planar Pt(II) complex"* (`logp_provenance.md`) — out-of-domain, not merely uncertain |
| **X4** | permanently charged, or substantially ionised at the assay pH | the compile has **no ionisation term**; the membrane leg consumes a neutral-species logP |
| **X5** | documented substrate of a transporter expressed in the assay system | active transport is not modelled; efflux ships default-off |
| **X6** | measured value is a bound (`>`, `<`), censored, or non-numeric | nothing to score against |
| **X7** | no descriptor available for the arm being scored | the same rule the HT-29 run already applies (`two_route_transport.py:289-291`), applied **per arm**, with the count reported per arm |
| **X8** *(B only)* | the measured quantity is `Kp,uu,brain`; or the value is imported from another study; or the row is a non-rat species | §3.4 E1–E3. `Kp,uu,brain` is efflux-confounded; imported values carry an unbounded inter-study difference; species may not be pooled |
| **X9** | the reported value was **DERIVED FROM A MODEL** rather than measured | §3.4 E5. Scoring against it scores Codon against somebody else's model — and where that model is a pH-partition model, against a structural cousin of Codon's own ion-trapping module. **Applies to every dataset, not only B**: any future anchor must be audited for imputed rows at the footnote level, not the table level |
| **X10** | cyclic peptide, or any solute whose membrane transport is governed by conformational shielding / intramolecular hydrogen bonding | §3.4 E6. Sits in the unvalidated gap between the small-molecule and antibody radius branches, and its permeability is set by a mechanism the engine's integer-donor-count membrane leg cannot represent |

**X5 requires work before the run.** The HT-29 score was protected by a fact about
the tissue: HT-29 is functionally P-gp-negative (Wils 1994, PMID 7945455) and 0 of
71 compounds carry an efflux assignment. **No other system carries that
protection.** Transporter status must be established from primary sources
**before** scoring; if it cannot be, that is a declared confound in the row, not a
silent assumption. Carried forward: HT-29 **does** express BCRP/ABCG2, that was
never quantified, and BCRP is unmodelled.

### 7.3 Mechanical protection

1. Retained and excluded lists with criterion codes, the named system, the
   admissibility decisions, and every numeric bar in this document are written to
   `results/sweep/transport_anchor_prereg.json` and **hashed**.
2. Committed **before** the scoring script may read any measured column.
3. Asserted against live code and the live compound list at a `--stage prereg` gate,
   following `experiments/two_route_transport.py --stage prereg` against
   `results/sweep/tworoute_prereg.json`.
4. **Both counts always reported, per arm**: how many the source contains, how many
   survived. For anything scored against Pruijn's 71, the §0.4 counting ambiguity
   must be restated.
5. An exclusion discovered *after* scoring that should have applied is reported as a
   **second, clearly labelled sensitivity arm**; the pre-registered arm remains the
   headline; both printed.

*(This document does not create that JSON. Writing it is the first task of whoever
runs a test, before measured values are read.)*

---

## 8. Power, partial data, and slots to fill

### 8.1 Power

**This ladder governs tests A and C. It does NOT govern test B** — §3.6 replaces it
with a bootstrap-interval rule, in advance and with the reason, because the ladder
was written for correlation-based tests and B's prediction has no ranks (§3.5).

| n retained | status |
|---|---|
| **n ≥ 20** | full verdict |
| **10 ≤ n < 20** | verdict reported, **flagged UNDERPOWERED**, with the **minimum detectable effect** at α = 0.05, power = 0.8 printed beside it |
| **n < 10** | **NOT TESTABLE**. No verdict, no number quoted as validation |

*"A test that cannot resolve the real effect produces 'no evidence', not 'evidence
of no effect'."* At 10 ≤ n < 20 the minimum detectable effect is printed
unconditionally, whichever way the result goes.

**Structural requirement:** the retained set must span **at least two distinct
chemical scaffolds or series**. A congeneric series alone tests MW-degeneracy but
cannot exercise the membrane leg, which is where the two-route law does its work.

### 8.2 The dataset must separate its own compounds

The measured spread across the retained set must exceed **3×** the median reported
measurement uncertainty, else **NOT TESTABLE**. `docs/VERIFICATION_METHODOLOGY.md`
§5 applied to the measurement:

> A dataset that does not itself separate its compounds cannot test a model's
> separation, and a passing score on one would be measuring the assay.

### 8.3 Partial data

1. **The retained subset must be defined by a rule independent of the measured
   value** — never by which rows were easiest to extract *after* their values were
   visible. The rule is stated verbatim in the report.
2. **If only summary statistics are obtainable**, only the spread comparison is
   computable; everything else is **NOT TESTABLE** and must be reported as such
   rather than approximated from a digitised figure without saying so.
3. **If n < 10 after §7**, report the shortfall and this document. **Do not
   substitute the HT-29 score and describe it as validation** (§0.1). Do not
   substitute a system chosen later on the basis of performance.
4. **Partial retrieval is a headline number**: *"n of N compounds"*, with N as
   published.

### 8.4 Slots to fill before any prediction is computed

Everything in §2.3, §3.2–§3.6, §4.2, §5.2–§5.4, §6.3 and §7.1 is **fixed**. These
may be filled:

| slot | needed for |
|---|---|
| GATE 1 verification against primary text, item by item; **the engine quantity chosen as commensurable, and the density/conversion assumption used** | §3.3 |
| re-verification of all nine Table 6 rows, the oxycodone footnote, and the ion-class assignments against the primary text | §3.1, §3.4 E1 |
| assay pH, buffer and slice-incubation conditions for the brain-slice measurement | X4, §3.10 reason 3 |
| the numeric **MW ceiling** | X2 |
| transporter status of each assay system, from primary sources | X5 |
| **n** retained after each cascade step (§3.1), per arm | §8.1, §3.1 |
| which donor convention `HBD_LOGP_DECREMENT` is defined for | §6.5 — governs whether test A runs one arm or two |
| the engine state (commit SHA + every §2.3 constant) at execution | §9.5 |
| the cyclic-peptide row list, so E6 is auditable rather than asserted | §3.4 E6 |
| the Fridén 2014 measured/derived row assignment, per row | §3.4 E5 |
| whether the source publishes **measurement uncertainty** | §8.2; if absent the loosest declared band applies and no tighter one may be invented |
| scorecard group key and PMIDs | §8.5 |

### 8.5 Registration in `codon validate`

Registered as an `AnchorGroup` with `scope = SCOPE_ANALYTIC` (**never**
`SCOPE_PDE`), carrying `coefficient-not-field`, plus `mcl-no-vasculature` for A. A
NOT TESTABLE outcome is registered too, with its reason. **Failures are never
suppressed and never softened**; an anchor that fails must appear and must make
`codon validate` exit 1.

---

## 9. Provenance, the ratchet, and anti-tamper

### 9.1 Commit discipline

1. **This file must be committed before any scoring script runs.** Its commit
   timestamp preceding the first run is the only evidence these bars were set in
   advance.
2. **The numeric bars must be mirrored into
   `results/sweep/transport_anchor_prereg.json`** and asserted against live code at a
   `--stage prereg` gate. `docs/VERIFICATION_METHODOLOGY.md` §6: *"Write the
   acceptance criterion into the analysis module, not the write-up."* This document
   is the write-up; the JSON is the binding copy. **This document does not create
   it.**
3. **Amendments** are recorded in the log below with the date, the reason, and an
   explicit statement of **whether any result was visible at the time**. An
   amendment made with results visible does not void the test, but it converts the
   bar from pre-registered to post-hoc and the report must say so.

### 9.2 The ratchet — the honest part

This claim has already been narrowed twice under pressure. *"The transport model
predicts tissue diffusivity"* became, after `diffusivity_falsification.md`, *"the
MW-only model is falsified"*. That became, after `two_route_promotion.md` §8, a
claim explicitly *"narrower than 'fixed'"* — conceding absolute accuracy, conceding
~76% of the measured spread, conceding that the test is one cell line. **Each
narrowing that follows a failure leaves a claim easier to satisfy than the one it
replaced.** Relative and ordering claims are cheap; absolute ones are not.

> If the claim is narrowed **again** after either test, the narrowing must be
> recorded in this file with the date, the result that caused it, and the sentence
> it replaced — so the sequence of retreats is visible in one place rather than
> distributed across seven documents.

A claim that has retreated four times and never been withdrawn is not a robust
claim; it is an unfalsifiable one not yet recognised as such. The only defence is to
make the retreats countable.

### 9.3 What would make either test unfalsifiable, and the guards

Both become unreachable *in practice* — while remaining reachable in principle,
which is worse, because the document would still look rigorous — if either is
allowed:

1. **The dataset or cell line is chosen after seeing how the model does on it.**
   Guard: the system is named in the frozen JSON before any prediction is computed;
   if several are available, **all are scored and all are reported**, with the
   primary named in advance. Reporting a subset is equivalent to choosing after the
   fact.
2. **Exclusion criteria move after residuals are visible.** Guard: §7.

### 9.4 `alpha` is not in either loop, and must not be laundered

`alpha` was re-solved at promotion (148.60815765445545 → 22.223477711172563). Both
tests score coefficients upstream of it, so `alpha` is not in the loop — a feature,
and it keeps them clean. But any statement derived from a penetration length
`L = sqrt(D_eff/k)` **is** `alpha`-conditional and inherits the `alpha-fit`
conditional. **A clean coefficient result must not be used to launder an
`alpha`-dependent one.**

### 9.5 What happens if a bug fix lands between registration and execution

This is live, not hypothetical: the `HD` convention defect (§6.5) is under
investigation and may move `alpha` and the Pruijn R² before these tests run. The
policy is fixed **now**, while nobody knows which way a fix would move a result.

**My position, which agrees with the coordinator's and adds three riders:**

> **The registration STANDS. The run is reported against BOTH engine states. The
> verdict is carried by the shipping engine AT EXECUTION TIME.**

**Why the registration stands.** A threshold is a statement about *what would count as
success*, not about *which commit produced the number*. A fix changes the model's
output; it does not change what a defensible result looks like. Re-deriving a
threshold after a fix — having seen which way the fix moved things — is post-hoc
regardless of how principled the re-derivation sounds. §6 of the methodology:
goalpost-moving is *"invisible after the fact — including to the person doing it."*

**Why the verdict follows the shipping engine.** That is what a user gets. The pre-fix
arm is reported as *the state this registration was written against*, so a reader can
see whether the fix moved the verdict. **If the fix changes the bucket, that fact is
the headline of the report, not a footnote** — a defect whose repair flips a verdict
is the most informative thing such a run can produce.

**Rider 1 — both arms must be recomputed in the SAME engine state.** THE TEST is a
ratio: arm M against arm B for test A, and `GMFE(engine)` against `GMFE_floor` for
test B. **Never compare a post-fix M against a pre-fix B.** The pre-fix numbers are
already published and will be sitting there, correctly formatted and tempting to
reuse; reusing them would silently convert a controlled comparison into an
uncontrolled one. This is `docs/VERIFICATION_METHODOLOGY.md` §9's matched-control
discipline extended across *time* rather than across arms, and it is the single most
likely way this run gets quietly corrupted.

**Rider 2 — the fix must be justified independently of these tests' residuals.** Its
justification must be committed **before** the run. A fix proposed *because* of a
residual observed here is a post-hoc model change, and it spends this anchor exactly
as the two-route term spent Pruijn (§0). **In that case the test must be
re-registered and re-run on fresh data**, and the original result reported as what it
was: the observation that motivated a model change. §0's lesson, applied
prospectively rather than discovered afterwards.

**Rider 3 — no mixing mid-run.** If a fix lands while a run is in progress, the run is
**void** and repeated end-to-end in a single engine state. Every run records the
commit SHA and the values of every constant in §2.3.

**One stability property, and it is worth stating.** `GMFE_floor` — test B's headline
finding (§3.5a) — is computed from the measurement alone and contains no engine
quantity. **It does not move under any bug fix, any recalibration of `alpha`, or any
resolution of the `HD` convention question.** Test B's most important number is
therefore robust to everything currently in flight, which is a good reason to run B
first and an argument for protecting that number as a headline.

### Amendment log

*(none — registered, not yet run)*

---

## Appendix — sources for every factual statement above

Nothing here is a result of either registered test, because neither has been run.

| statement | source |
|---|---|
| MW-only model falsified on 71 HT-29 compounds; R² = 0.022; spread 2.3× vs 29.7×; 0 of 31 exact-mass pairs | `results/diffusivity_falsification.md`; commit `d3b9955` |
| Two-route law built in response, scored on the same 71 at R² = 0.718, promoted to default | `results/two_route_promotion.md`; commit `fd2a9e7` |
| Spread 8.3× vs 29.7×; within-3× 38% → 99%; Spearman +0.821; bias 1.03×; 24 of 31 pairs | `results/two_route_promotion.md` §0, §4 |
| `alpha` 148.60815765445545 → 22.223477711172563 | `results/two_route_promotion.md` §0 |
| Spread tracks HBD (R² = 0.305) and logP₇.₄ (R² = 0.290), not MW (R² = 0.028) | `results/diffusivity_falsification.md` §7 |
| Intracellular retardation refuted on these compounds | `results/intracellular_retardation.md` |
| Efflux off; HT-29 functionally P-gp-negative; HT-29 does express BCRP | `results/two_route_efflux.md`; Wils 1994, PMID 7945455 |
| Prereg frozen pre-scoring; "DO NOT FIT TO PRUIJN"; decrement moved away from the optimum | `results/sweep/tworoute_prereg.json`, `results/sweep/tworoute_promotion_prereg.json` |
| Under-separation predicted in advance; ~24% of spread; one cell line | `results/two_route_promotion.md` §8.1, §8.2 |
| Absolute `D_eff` still misses its band (1.50× low, previously 2.21× high) | `results/two_route_promotion.md` §8.3 |
| `cell_partition` = K, cytoplasm/interstitium equilibrium partition, "an unbound partition ratio across one membrane", shipped at 1.0 | `tmesim/molecular/two_route.py:357, 389-391, 485-487` |
| Pruijn 2008 tabulates per-compound `Dmcl` for HT-29 only; per-line thicknesses back-computed from urea; asymptote `a` spread 1.35 log₁₀ ≈ 22×; SiHa `c` printed identical to HT-29 including SE | acquisition agent, primary PDFs, 71/71 rows cross-checked against the repo literal with zero discrepancies |
| Count ambiguity: Table S1 = 71, Table 1 n = 69, three "in table but not in graph", three with no logP, 71−3 = 68 | acquisition agent, primary source |
| **All prior Pruijn scoring used the paper's own logP₇.₄ and HD columns, not Crippen** | `experiments/two_route_transport.py:246-247` (verbatim) |
| Crippen neutral logP vs measured logP₇.₄ warning; the 68-compound exclusion rule | `experiments/two_route_transport.py:492-495, 289-291` |
| Every logP is Crippen; doxorubicin correction moves `alpha` by 0.462×; no Pt(II) parameters; stereochemistry-blind | `results/logp_provenance.md` §1, §3 |
| `logp_provenance.md` §1 finding 2 is stale against the shipping engine | file dated 2026-08-08 vs promotion 2026-08-14; no `two_route` reference, no correction notice (verified by grep) |
| φ_cell admissible range 0.40–0.85 published as a sweep (arm B, impermeable-cell) | `results/cellular_exclusion.md` §4.3 |
| Frozen parameter values | `tmesim/molecular/two_route.py:245,254,357`; `results/sweep/tworoute_prereg.json`; `experiments/anchor_prototype.py:84` |
| Bars of 3.0 already in use and already producing a FAIL | `scorecard.py::_run_gist_imatinib`; `pk_reality.ORDER_OF_MAGNITUDE_FACTOR`; `glioma_regional_pk.BAND_TIGHT`; `delivery_metrics.md` §3b |
| PASS / FAIL / NOT TESTABLE / DEFERRED semantics; failures never suppressed | `tmesim/validation/scorecard.py` module docstring |
| `coefficient-not-field`, `mcl-no-vasculature` conditionals; the `analytic`-vs-`PDE` rule | `tmesim/validation/scorecard.py::CONDITIONALS`; `results/anchor_registration.md` §2 |
| The units-error signature precedent | `tmesim/validation/scorecard.py::_run_ps_scale` |
| Anchor definition; two prior modality-pooling regressions | `results/external_anchor_survey.md` §0 |
| Compensating errors as the most dangerous defect class; reachability; rank correlation cannot see calibration; signed metrics defend cancellation; pre-declared bars; stale results files | `docs/VERIFICATION_METHODOLOGY.md` §1, §2, §5, §6, §7, §8, §10c |
| 52 scoreable `Vu,brain` values, 35 in rat, range 0.436 → 4184 mL·g⁻¹ = 9596×; 1308× excluding cyclic peptides | six open-access primary papers sharing the Uppsala slice method, via the acquisition agent. **Every row to be re-verified against primary text at GATE 1** |
| oxycodone's *"Data from Fridén et al., 2011"* footnote | Loryan 2013 Table 6 — PMID 23336814, PMC3602653 — **excluded, §3.4 E1** |
| Replicate pairs: paroxetine 714 / 876 (1.23×), diazepam 17.8 / 20 (1.12×) → declared floor **1.25×** | two independent rat papers — §3.4a |
| Mouse cross-strain overlaps agreeing to < 2 % with a bit-identical oxycodone | acquisition agent — judged **implausibly good for independent measurement, likely value reuse; flagged UNUSABLE and not quoted** (§3.4a) |
| Fridén 2014 footnote h: rows carrying an `fu,brain` value had `Vu,brain` **calculated from a pH-partition model**; agent split 23 measured / 24 derived; body text says 31/16 | Fridén 2014 — **derived rows excluded, §3.4 E5**. Both splits sum to 47, so the disagreement is over *which* rows, not how many |
| `Kp,uu,cell`: 4 scoreable brain values, 54× spread | acquisition agent — **test C NOT RUNNABLE, §3.8** |
| `Kp,uu,brain` captured to `kp_uu_brain_DO_NOT_SCORE.csv` | acquisition agent — §3.4 E3 |
| Shipping descriptor pipeline projected at R² ≈ 0.22–0.50, central ≈ 0.27 vs published 0.718 (~62 % loss), on a **five-compound** error estimate because 66 of 71 are ACSRC `SN` codes with no resolvable structure | `results/descriptor_falsification.md` — §6.1, §5.2b |
| R² ≈ 0.45 consistent with correct physics + σ ≈ 1 descriptors; < 0.20 the first regime implicating physics | `results/descriptor_falsification.md` — encoded as §5.2b, declared order-of-magnitude guides |
| Pruijn `HD` = donor hydrogens (`NHOHCount`); engine computes `Lipinski.NumHDonors` = donor heavy atoms; urea paper 4 / engine 2 | coordinator, verified directly — §6.5. Investigation in progress; **nothing changed** |
| `HBD_LOGP_DECREMENT = 1.34 = 3.621 × 0.37`, 0.37 the acidity of a *typical donor group* | `results/sweep/tworoute_promotion_prereg.json` — the convention ambiguity of §6.5 |
| "A consistently-applied wrong convention validates itself"; the `1/ln²2` α error | `docs/VERIFICATION_METHODOLOGY.md` §10b |
