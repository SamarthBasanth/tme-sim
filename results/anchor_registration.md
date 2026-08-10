# Anchor registration — the Hicks multicellular-layer coefficients

**What this file records.** Three external anchors identified by
`results/external_anchor_survey.md` and re-scored by `results/cellular_exclusion.md`
were sitting unregistered: nothing in `vivum validate` scored them, so they could
not fail in CI and did not appear in any tally. They are now registered as the
anchor group **`mcl-transport-coefficients`** in `tmesim/validation/scorecard.py`.
This file states what was registered, the scope justification for each row, the
conditionals each row inherits, an audit of the previously-registered
`ifp-wick-in-needle` group, the new `vivum validate` tallies, and what is still
unregistered and why.

**The rule this registration was written under.** Nothing was tuned. No verdict
was softened. The one decision that mattered — the scope tag — was made against
the project's own interest, and §2 explains why.

---

## 0. Bottom line

| | before | after |
|---|---|---|
| `vivum validate --quick` | PASS 7, FAIL 4, NOT TESTABLE 3, DEFERRED 8 (22 rows) | PASS 7, **FAIL 6**, **NOT TESTABLE 4**, DEFERRED 8 (**25 rows**) |
| the quick-tier disclosure line | *"11 of the 11 scored rows never invoke the transport PDE"* | ***"13 of the 13 scored rows never invoke the transport PDE"*** |
| `vivum validate` (full) | PASS 13, FAIL 18, NOT TESTABLE 16 (47 rows); 19 of 31 scored rows not PDE | PASS 13, **FAIL 20**, **NOT TESTABLE 17** (**50 rows**); **21 of 33** scored rows not PDE |

**Two failures were added and zero passes.** The registration made the scorecard
strictly worse-looking, which is the correct direction for an honest one.

**The headline gap did NOT close, and this registration could not close it.** The
"N of N never invoke the PDE" line moved from 11/11 to 13/13 — the *count* went up
by two and the *fraction* stayed at 100%. These three anchors measure the transport
PDE's **input coefficients**, not any field it solves, so tagging them `PDE` would
have moved that line by two rows without a single extra solve. That would have
improved the number by defeating the disclosure the number exists to make. §2 is
the full argument.

---

## 1. What was registered

Anchor group `mcl-transport-coefficients`, quick tier (no tissue, ~1.6 s), scope
`analytic`. Three rows, registered separately because the engine can be wrong
about any one of them without being wrong about the others — and is.

Measured constants are **imported from `experiments/anchor_prototype.py`**, which
already owns them and traces each to its PMID; no measured number is re-typed into
the scorecard. Predictions come from the unmodified molecular pipeline
(`descriptors.from_smiles` → `compile_drug.build_profile`) at the repo's own
`ecm_reference = 0.4`.

### Sources, verified against PubMed

According to PubMed, both records were retrieved and every value below is quoted
from the retrieved abstract, not inferred:

> Hicks KO, Fleming Y, Siim BG, Koch CJ, Wilson WR. "Extravascular diffusion of
> tirapazamine: effect of metabolic consumption assessed using the multicellular
> layer model." *Int J Radiat Oncol Biol Phys* 1998;42(3):641-9. **PMID 9806526**,
> [DOI](https://doi.org/10.1016/s0360-3016(98)00268-5)
>
> The record states verbatim: *"The flux of TPZ through well-oxygenated MCL
> (equilibrated with 95% O2) was well fitted as Fickian diffusion without
> reaction, with a D(M) of 7.4 x 10(-7) cm2s(-1) (12-fold lower than in culture
> medium) for V79 and 1.3 x 10(-6) cm2s(-1) for MGH-U1 MCL. Flux of TPZ was
> suppressed under anoxia, and fitting the data required inclusion of a reaction
> term with a rate constant for metabolic consumption of TPZ of 0.52 min(-1) for
> V79 and 0.31 min(-1) for MGH-U1 MCL."*

> Hicks KO, Myint H, Patterson AV, Pruijn FB, Siim BG, Patel K, Wilson WR.
> "Oxygen dependence and extravascular transport of hypoxia-activated prodrugs:
> comparison of the dinitrobenzamide mustard PR-104A and tirapazamine."
> *Int J Radiat Oncol Biol Phys* 2007;69(2):560-71. **PMID 17869669**,
> [DOI](https://doi.org/10.1016/j.ijrobp.2007.05.049)
>
> The record states verbatim: *"The diffusion coefficient of PR-104A in
> multicellular layers (4.42 +/- 0.15 x 10(-7) cm2 s(-1)) was lower than that of
> tirapazamine (1.30 +/- 0.05 x 10(-6) cm2 s(-1))."*

Every number used below appears in one of those two quotations. Nothing was taken
from a title, and no uncertainty was invented where the record states none.

**One detail worth carrying, because it is not in either summary table.** The
diffusion coefficients were fitted to flux through **well-oxygenated** MCLs, where
the data were *"well fitted as Fickian diffusion without reaction"* — so `D` is
measured in the *absence* of the reaction term. The consumption constant was
fitted separately, **under anoxia**. `D` and `k` therefore come from two different
oxygenation states of the same system, which is a reason to score them as two rows
and never as a single implied penetration length.

### The three rows

| row | measured | engine, as it ships | verdict | bar |
|---|---|---|---|---|
| effective diffusivity | **74–130 µm²/s**, three MCL cell lines (V79-171b 74, MGH-U1 130, PMID 9806526; SiHa 130 ± 5, PMID 17869669) | **287.3 µm²/s** (3.88× the lowest, 2.21× the highest) | **FAIL** | inside the measured 74–130 µm²/s range |
| hindrance `D_free/D_eff` | **12-fold** vs culture medium (V79-171b, PMID 9806526; no published uncertainty) | **1.77-fold** (0.15× measured) | **FAIL** | 3× either side of 12, i.e. **4.0–36.0-fold** |
| first-order consumption `k` | **8.667e-03 /s** (0.52 /min, V79-171b) and **5.167e-03 /s** (0.31 /min, MGH-U1) | two readings, 43× apart, straddling the target | **NOT TESTABLE** | n/a |

#### Row 1 — effective diffusivity: **FAIL**

287.3 µm²/s against a measured 74–130. The three measurements agree with **each
other** to 1.8× while the engine misses all three **in the same direction** — a
systematic bias, not scatter. It does pass the repo's loosest declared band
(`pk_reality.ORDER_OF_MAGNITUDE_FACTOR = 3.0`, i.e. ≤ 390 µm²/s), and the row says
so; that is not the bar adopted, because `external_anchor_survey.md` §4
recommended this anchor **as a FAIL row** and that recommendation is the bar.

#### Row 2 — hindrance: **FAIL, and it lands SHORT rather than passing on improvement**

This is the sharpest row available anywhere in the registry: a pure dimensionless
ratio of one drug's diffusivity in tissue versus in free medium. No units, no
compartment translation, no cohort matching, no counting convention.

The measured 12-fold carries **no published uncertainty**, so a band had to come
from somewhere. Rather than invent one, the row uses the repo's own loosest
declared band — the 3× order-of-magnitude factor `pk_reality` already ships —
giving 4.0–36.0-fold. **That is the most generous bar available and the shipping
engine still misses it**, at 1.77-fold.

The row records the cellular-exclusion arm explicitly and records that it is
**short**:

* shipped: **1.77×**
* with the merged-but-default-off cell term, arm B (cell term alone, the
  MCL-appropriate arm): **6.6157×** — a 3.7× improvement that closes most of the
  6.8× error and **still lands 1.81× short of the measured 12×**
* arm A (full existing chain plus the cell term): 11.73×, close to 12 — **and not
  the result**, because arm A credits an MCL with ECM the MCL does not have. Its
  proximity to the measurement is partly the wrong term doing the right work
  (`cellular_exclusion.md` §4.1).

Improving is not passing, and the row does not record it as one. (`cellular_exclusion.md`
§0 prints arm B as `6.61×`, rounded down from 6.6157; the scorecard prints the
exact value in `values` and says so in the note, so the two cannot read as a
discrepancy.)

#### Row 3 — consumption rate: **NOT TESTABLE**, and the direction of that call is stated

Two independent reasons, either sufficient:

1. **The process is not the same process.** The measured constant is anoxic
   bioreductive metabolism of a hypoxia-activated prodrug — a real, drug-specific
   consumption reaction the engine's panel drugs do not have — while the engine's
   `k` is a generic binding-plus-clearance sink. Same units, same slot in the same
   equation, **different physical process**. `external_anchor_survey.md` §2.1, the
   owning document, declares it a measured **order-of-magnitude bound** and
   explicitly *not* a like-for-like score.
2. **The engine has no single number to score.** `k_raw = clearance + binding =
   9.747e-05 /s` is **0.0112–0.0189×** measured (53–89× too small). `k_calibrated =
   clearance + α·binding = 4.231e-03 /s` at α = 148.6082 is **0.488–0.819×**
   measured. The two differ by 43× and **straddle** the measurement. Nothing about
   this anchor determines which reading is the engine's `k`.

**Which way this call costs, stated in the row itself.** Scored on `k_calibrated` —
the value the solver actually uses — this row would **PASS** the repo's 3× band.
Declining to score it **removes a pass from the tally, not a failure**. That is the
opposite of softening, and the row says so in its `reason` so a reader cannot
mistake it for a hidden FAIL.

The FAIL-shaped half is not suppressed either: the reason records that the
uncalibrated sink is 53–89× too small, which independently reproduces
`results/pk_reality.md`'s ~300× diagnosis on a *different drug*, in a *different
system*, against a *directly measured rate constant* rather than a penetration
depth.

---

## 2. The scope decision — the crux

**All three rows are tagged `analytic`. None is PDE evidence.**

The transport PDE is `∂C/∂t = ∇·(D∇C) − kC + S`. `D_eff` and `k` are its **input
coefficients**. The engine produces them in closed form from molecular descriptors:

```
SMILES -> descriptors.from_smiles -> compile_drug.build_profile -> MolecularProfile
          (.effective_diffusivity, .hindrance, .clearance, .tissue_binding_rate)
```

No solve occurs anywhere in that chain. The prediction is read off a dataclass.

This is exactly the situation `external_anchor_survey.md` warns about when it notes
that logP-type quantities never touch the solver, and it is the same category
`glioma_regional_pk` and `gist_imatinib` are already registered under. So:

> **The honest answer is that these coefficients are analytic inputs, not solver
> outputs.** They can falsify the hindrance law. They cannot confirm or refute the
> solver.

**Why this mattered more than usual here.** Tagging these three rows `PDE` would
have moved the summary's *"N of N scored rows never invoke the transport PDE"* line
from 11/11 to 11/13 in the quick tier — an apparent improvement in the project's
single most-quoted honesty statistic — **purchased entirely by mislabelling, with
no new solve and no new evidence about the solver**. That would be worse than
leaving the line where it was, so the line was left where it was. It now reads
13/13.

**What these rows ARE worth, stated positively.** `results/alpha_identifiability.md`
establishes that the project has essentially **one** measurement constraining
either PDE coefficient, and that one (the doxorubicin perivascular half-distance,
PMID 16361566) is in-sample by construction because α is solved to it in closed
form. These three rows are the first **out-of-sample, absolute, same-equation**
constraints on the coefficients themselves. Tirapazamine is not on the engine's
panel and nothing in the repo was ever fitted to a TPZ measurement. That is a real
gain; it is just not a gain in solver evidence.

---

## 3. Verification that the solver is never entered

The `analytic` claim was not asserted, it was measured. Two independent instruments
were run over `_run_mcl_transport_coefficients()`:

1. **Attribute wrappers** on `tmesim.physics.diffusion.simulate_drug_transport` and
   `tmesim.sim.engine.simulate_response`.
2. **`sys.setprofile`** over the whole runner, recording the source file of every
   Python function actually entered — which catches indirect and renamed calls that
   an attribute wrapper would miss.

Result:

```
rows produced: 3
  FAIL         scope=analytic  physics=n/a — transport solver never invoked for this row
  FAIL         scope=analytic  physics=n/a — transport solver never invoked for this row
  NOT TESTABLE scope=analytic  physics=n/a — transport solver never invoked for this row

attribute-wrapper call counts: {'simulate_drug_transport': 0, 'simulate_response': 0}
functions entered inside physics/diffusion.py or sim/engine.py: NONE

VERDICT: SOLVER NEVER ENTERED — scope=analytic is correct
```

The profiler half of that check is now a **standing regression test**
(`tests/test_validate_cli.py::TestRegistryHygiene::test_mcl_coefficient_anchors_are_analytic_and_never_enter_the_solver`),
so the scope tag cannot be quietly upgraded later: any future change that makes
these rows enter the solver, or that relabels them `PDE` without doing so, fails
the suite.

---

## 4. Commensurability — what a multicellular layer does and does not license

This project has manufactured a false regression twice by pooling incompatible
modalities, so the boundary is stated in the row itself (conditional
`mcl-no-vasculature`), not in a footnote here.

**A multicellular layer is a stack of cells grown on a porous support: packed
cells, little stroma, NO VASCULATURE, no necrosis, no interstitial convection and
no plasma source term.** Hicks measured flux across that slab and fitted it to a
diffusion-reaction equation **of the same form the engine solves**. That is
precisely what makes `D` and `k` commensurable **as coefficients**: same units,
same functional role, same equation, no unit conversion and no compartment
translation.

**What it licenses.** Scoring the engine's compiled per-drug `D_eff`, and the
hindrance law that produces it, against a directly measured tissue diffusion
coefficient. The hindrance row in particular is informative *because* an MCL is
almost pure packed cells: it asks whether the hindrance law accounts for cellular
volume exclusion at all, and the answer as shipped is that it does not.

**What it does NOT license.** Any claim about a vascularised tissue section.
Nothing in an MCL constrains vessel spacing, the transvascular source term,
IFP-driven convection, or the spatial concentration field the solver produces on
real reconstructed geometry. An MCL has no vasculature, so it cannot corroborate or
refute a single thing the engine does with vasculature. Reading these rows as
evidence that the engine's tissue predictions are validated would be the exact
modality pooling `external_anchor_survey.md` §0 catalogues three live instances of.

**A second, narrower mismatch, disclosed rather than smoothed over.** V79-171b is a
Chinese hamster line and MGH-U1 and SiHa are human; the measured band pools three
cell lines across two publications. The band is used as a band — 74–130 µm²/s —
and no per-line verdict is claimed.

---

## 5. Conditionals recorded

Three new keys were added to `CONDITIONALS` and are printed in full by
`vivum validate` whenever any row cites them.

### `cell-exclusion-default-off` — cited by rows 1 and 2

States that `molecular/ecm_mesh.py::cell_exclusion` is **merged but not wired in**:
`cellular_exclusion.md` §7 recommends merging the physics and explicitly *not*
flipping the default, nothing in `compile_drug` calls it, and
`DrugParams.diffusivity` is bit-identical to the pre-merge value. **Every verdict
in this group is therefore scored on the engine AS IT SHIPS, with the term OFF.**
The cell-term-on arm is reported alongside each verdict and never supplies one.

The three caveats the conditional carries, all in the row rather than elsewhere:

* **(a) The wrong system's packing density.** φ_cell = 0.7892 is measured by
  Delesse's principle on the repo's real Xenium **breast tumour** section and then
  applied to a **V79 hamster multicellular layer**. A literature search for an MCL
  extracellular-volume measurement returned nothing usable
  (`cellular_exclusion.md` §4.1). If an MCL is denser than the breast section —
  the likely direction, since an MCL has no vasculature, no necrosis and little
  stroma — then arm B **under-states** the hindrance and the residual 1.8× is an
  over-estimate.
* **(b) φ_cell is an upper estimate.** Xenium's default boundary segmentation
  expands each nucleus until neighbours meet, so it tiles the section and
  over-states cell area where cells are sparse. Every hindrance derived from it is
  an upper estimate too.
* **(c) The ~1.9× double-count.** Stacking the cell term on the existing chain
  double-counts by up to ~1.9×, because `brownian.py`'s 300 nm obstacle field
  (whose own docstring calls its obstacles *"cells / ECM aggregates"*) and
  `ALPHA_TORTUOSITY = 0.6` were partly standing in for cell-scale obstruction.
  `cellular_exclusion.md` §5.3 states in terms that it **cannot resolve which**.

**Why this is a conditional and not a blocker.** None of (a)–(c) affects a single
verdict in this group, because no verdict is computed from the cell-term-on arm.
They are recorded because the arm's *numbers* are quoted in the rows, and a quoted
number without its caveats is how the survey's §0.1 failures happened.

### `mcl-no-vasculature` — cited by all three rows

The commensurability bound of §4, in full, in the row.

### `coefficient-not-field` — cited by all three rows

States that the row scores an **input coefficient**, not a solved field; that the
solver is never entered; that this was verified by call counter (§3); and the
consequence — *"it can FALSIFY the hindrance law, and it cannot confirm or refute
the solver."*

Row 3 additionally cites the pre-existing **`alpha-fit`**, because its
`k_calibrated` reading is downstream of the one fitted constant.

`analytic-only` was deliberately **not** cited by these rows, even though the
physics claim in it is accurate. Its text attributes the verification to
call-counter instrumentation in
`results/transvascular_exchange_revalidation.md` §6.1, and that instrumentation
covered `glioma_regional_pk` and `gist_imatinib` — not these rows. Citing it would
have attached a provenance claim these rows have not earned, so
`coefficient-not-field` carries their own verification instead.

---

## 6. Audit of the existing `ifp-wick-in-needle` registration

`results/ifp_validation.md` registered `ifp-wick-in-needle` before this work. It is
the first genuinely PDE-invoking registered anchor, and it holds up.

**Registration is complete and correct.** Three rows, all `scope=PDE`, all carrying
a declared delivery law (so none is downgraded by `_stamp_physics`), all with
PMIDs, and the group is gated on both real caches.

* GBM vs Boucher 1997 (PMID 9062403) — **NOT TESTABLE**, on the strongest available
  grounds: the anchor cannot discriminate. Correct, and the reason is the right
  one — an anchor whose band a model producing *nothing* would pass, and which
  cannot separate tumour from normal tissue (Less 1992, 2 ± 4 mmHg, same technique,
  overlapping authors).
* Breast vs Less 1992 (PMID 1423283) — **FAIL**, the discriminating row. Confirmed
  live: it is one of the 12 PDE-scoped scored rows in the full run, and the only
  PDE FAIL outside the two kill/transvascular groups.
* The numerical-resolution row — **NOT TESTABLE**, registered so the two above
  cannot be quoted without it. Correct.

**One defect found and fixed.** The GBM row cited `ifp-uncalibrated`,
`ifp-unresolved-mesh` and `untreated-section` but **not** `ifp-grid-unconverged` —
while that conditional's own text describes the **GBM** ladder explicitly
(0.358 → 0.434 → 0.462 → 0.370 mmHg across dx = 24/18/12/8 µm, non-monotone, 1.29×
spread). The median the GBM row publishes is one point on that ladder. A row that
declines to score still **publishes a number**, and that number is no more
grid-converged than the breast one. `ifp-grid-unconverged` was added to the GBM
row.

This is a strictly-more-disclosure change: it alters no verdict, no bar and no
predicted value, and it makes one more caveat print.

**One thing deliberately left alone.** The breast row's note reports ~40× low
against Less 1992 and ~78× low against Nathanson 1994 (PMID 7850532), and flags
that Nathanson's uncertainty is a **standard error, not an SD** and must not be
read as one. That is right and was not touched.

---

## 7. New tallies

Both tiers re-run after registration. `transvascular_exchange=False` (the shipping
law) in both.

### `vivum validate --quick`

| | before | after |
|---|---:|---:|
| PASS | 7 | **7** |
| FAIL | 4 | **6** |
| NOT TESTABLE | 3 | **4** |
| DEFERRED | 8 | 8 |
| total rows | 22 | **25** |
| scored rows | 11 | **13** |
| scored rows that never invoke the PDE | 11 | **13** |

> *"13 of the 13 scored rows never invoke the transport PDE. They cannot be cited
> as evidence about the solver."*

### `vivum validate` (full)

| | before | after |
|---|---:|---:|
| PASS | 13 | **13** |
| FAIL | 18 | **20** |
| NOT TESTABLE | 16 | **17** |
| DEFERRED | 0 | 0 |
| total rows | 47 | **50** |
| scored rows | 31 | **33** |
| scored rows that never invoke the PDE | 19 | **21** |

Exit code 1 in both tiers, before and after.

### Did the "never invoke the PDE" line move?

**Yes, and in the wrong direction — which is the honest outcome.**

* **Quick tier: 11 of 11 → 13 of 11+2 = 13 of 13.** The fraction is unchanged at
  **100%**. Two more scored rows exist and neither of them touches the solver.
* **Full tier: 19 of 31 → 21 of 33.** The fraction moved from 61.3% to 63.6%, i.e.
  the proportion of scored rows that are *not* solver evidence went **up**.

**This registration did not close the validation hole. It measured two more things
that are not the hole.** The count of PDE-scoped scored rows is unchanged at 12 in
the full run and 0 in the quick run. The only registered anchor that scores the
solver against an external measurement on the strength of a solve remains
`ifp-wick-in-needle`'s breast row, and it FAILS.

**What would actually move the line**, from `external_anchor_survey.md` §5: nothing
public was found that measures a **spatially resolved free-interstitial
concentration field** in tumour tissue this engine could reconstruct. Microdialysis
gives one point per probe with no coordinates; MCLs give slab flux with no
vasculature; perivascular immunofluorescence exists only for autofluorescent or
antibody-detectable drugs. Registering more coefficient anchors — however good —
cannot change that ratio, and this file should not be read as progress toward it.

---

## 8. What remains unregistered, and why

Each of these is now recorded in `scorecard.OUT_OF_SCOPE`, which `vivum validate
--list` prints, so the omissions are auditable rather than invisible.

| candidate | why not registered |
|---|---|
| **Pruijn 2008 `D_mcl` for 71 compounds** (doi:10.1071/CH08240), and the bounded version in PMID 16912264 (2.7–187 µm²/s across 16 TPZ analogues in HT29 MCLs) | **The highest-value unregistered anchor identified, blocked on retrieval and not on judgement.** The engine derives `D_eff` from a Stokes–Einstein radius that is a function of **MW alone** times a hindrance depending on the ECM fraction and **not on the molecule**, so across a congeneric series it must predict a near-constant `D_eff` while the measurement spans **69-fold**. That is a decisive falsification of the diffusivity model needing only the compound list. The publisher returns HTTP 403 to automated retrieval and the record is not in PubMed or PMC, so the compound-by-compound table is **UNVERIFIED** and no row can be built without inventing numbers. |
| **Capecitabine tumour:plasma 1.76, tumour:s.c. 1.11 in human breast lesions** (PMID 12618890) | Commensurable (free-interstitial over **ultrafiltered** plasma) and recommended for adoption by the survey §4, but **no module computes it**. It would reuse `experiments/breast_kp_anchor.py`'s solve path against a new anchor table — new code in a file this registration does not own. Registering a row whose runner does not exist is not possible. **When built, score the parent drug only**: DFCR, DFUR and 5-FU are metabolites generated intracellularly en route, so their interstitial concentration carries a local generation term the engine does not model. This one *would* be a PDE-scoped row. |
| **Tunggal 1999 (PMID 10389947) and Cowan 2001 (PMID 11149410) MCL flux ratios**, on four of the engine's own panel drugs | **Not actionable without new geometry, and the omission is expensive.** "% of the rate through the blank Teflon membrane" is a ratio of steady-state fluxes through a **1-D slab with reaction**; the engine solves 2-D reconstructed tissue with vessel sources and has no 1-D slab configuration. The measurement says doxorubicin and mitoxantrone are 3–5× slower than methotrexate while the calibrated panel puts doxorubicin at `L` = 64.9 µm and methotrexate at 75.2 µm — essentially equal, i.e. a **directional disagreement on the engine's own drugs** that nothing currently scores. Cowan's concentration invariance from 10 µM to 3 mM is separately a direct test of the engine's first-order-sink linearity assumption, and it **confirms** it over 300-fold. |
| PR-104A `D_mcl` = 44.2 ± 1.5 µm²/s (PMID 17869669) | Verified and available, but PR-104A was not compiled through the pipeline in the prototype this registration drives, so there is no engine number to score against it. Adding it would be a second drug on the same anchor class rather than a new class of evidence. Recorded here so it is not re-found as novel. |
| Jørgensen 2024 (PMID 38332155), Grossman 2013 (PMID 23649683), Ekstrøm 1997, Brana 2014 (PMID 25213162) | Rejected by the survey §4 on substance and not re-litigated: mixed denominators with no tumour; no plasma leg; n = 1 with no recovery correction; unretrievable. |

**One open item that would sharpen a registered row rather than add one.** A
measured cell (or extracellular) volume fraction for a Wilson/Hicks multicellular
layer would turn the hindrance row's cell-term-on arm from a bracket into a score,
by removing caveat (a) of §5. `cellular_exclusion.md` §7 item 4 lists it as a
literature question, not a compute question, and it is still open.

---

## 9. What changed

* `tmesim/validation/scorecard.py` — three new `CONDITIONALS` keys
  (`mcl-no-vasculature`, `cell-exclusion-default-off`, `coefficient-not-field`);
  the constant `PHI_CELL_BREAST_SECTION` (documented, used only to report the
  cell-term-on arm, and determining no verdict); the runner
  `_run_mcl_transport_coefficients`; the `mcl-transport-coefficients` registry
  entry; `ifp-grid-unconverged` added to the IFP GBM row; three new
  `OUT_OF_SCOPE` entries.
* `tests/test_validate_cli.py` — three regression tests, the load-bearing one being
  the profiler check that these rows never enter the solver.
* `results/anchor_registration.md` — this file.

Nothing under `tmesim/molecular/`, `tmesim/physics/`, `tmesim/reconstruct/`,
`app/`, `landing-3d/` or `docs/` was touched, and no existing results file was
modified. Test suite: **757 passing** (754 baseline + 3 new).
