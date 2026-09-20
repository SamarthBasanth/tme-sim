# Verification routes — what this engine can be scored against, today, without a clinical outcome

**Date: 2026-09-20. Status: SURVEY. Nothing here is a result.** No number below is an
outcome of running anything. Every number is either (a) a published measurement I
fetched in this pass, with the URL and the sentence I am relying on recorded in §6,
or (b) a value already published elsewhere in this repository and cited to its file.
No code was changed, no anchor was registered, no existing file was edited.

**Scope of the commission.** *"Scour far and wide for verification processes, without
being overly cynical."* The second clause is the harder half. This project has twice
talked itself out of a legitimate result — once by moving a goalpost after a
pre-registered anchor failed, once by writing off sections whose per-cell statistic is
refused in code for a reason that did not apply to the question being asked. Both
directions are errors. So this document is deliberately asymmetric in one respect: it
spends most of its length on **what can be done**, and confines *"what cannot be
verified without clinical data"* to a short §5, because that section is already the
most thoroughly written part of this repository (`results/DATA_HUNT.md` §6,
`docs/DATA_ACCESS_PLAN.md`, `results/landscape_response_platforms.md`) and does not
need a fifth telling.

**What this document does not license.** Nothing here produces or supports a
cross-compound ordering — withdrawn 2026-08-21, absolutely, and every route below is
either a *within-compound* comparison, a *within-tissue* comparison, or a coefficient
scored against its own measurement. Nothing here produces an absolute kill fraction, a
clinical-response prediction, or a diagnostic or treatment-planning output. Several
routes are compound-resolved; **compound-resolved is not cross-compound ranking**, and
a PASS on any of them would not make it one.

**Counts.** 14 routes. **10 verified obtainable today** by me, in this pass, with a
fetched URL or a live API query (R1, R2, R3, R4, R5, R6, R9, R10, R11, R12). **2
partially obtainable** — the primary record exists and is open, but the specific
numbers I would score against sit behind a paywall or in a figure (R7, R8). **2 open
but weak on the delivery axis** and ranked accordingly (R13, R14). **11 candidate
routes ruled out as already dead**, listed in §1.2.

---

## 0. How to read this

Each route carries the same six fields, in the same order:

| field | what it means |
|---|---|
| **Tests** | the one quantity under test, stated narrowly enough to fail |
| **Engine predicts** | what the shipped configuration currently emits for it, cited to file or to a `results/` document |
| **Data** | the measurement, with accession/DOI/PMID, and whether *I* obtained it or only located it |
| **Effort** | engineer-days, where one day = a day of one competent engineer with the repo already checked out and `.venv` working. Excludes writing the report |
| **If it fails** | what a FAIL would tell you. **A route where failure teaches nothing is not worth running** and is marked as such |
| **If it passes** | the exact sentence a PASS licenses — written now, before the run, so it cannot be inflated afterwards |

**One rule inherited from `results/transport_anchor_prereg.md` §0.1 applies to every
route below:** a dataset that selected a mechanism is spent as a test of that
mechanism. The 71 HT-29 compounds of Pruijn 2008 selected the two-route law's
structure. **No route below scores against them**, and any route that touches
multicellular layers is explicitly built on a *different* laboratory, cell line and
compound set for exactly this reason.

---

## 1. The baseline — so that nothing here is re-proposed

### 1.1 What is already scored

17 anchor groups are registered in `tmesim/validation/scorecard.py::build_registry()`.
Five touch the PDE. The picture, from `results/CLAIM_INVENTORY.md` §1e and
`results/anchor_registration.md` §0:

* **Three scored rows compare a *solved field* to an *external measurement*, and all
  three FAIL** — breast IFP 41× low, 5-FU breast Kp 12.8× high, methotrexate breast Kp
  10.5× low. The last two fail in *opposite directions*, so no single constant absorbs
  them.
* The remaining nine PDE rows compare a *ratio* to a literature band and are all
  conditional on the one fitted constant, `alpha = 22.223477711172563`.
* `results/kyle2004_penetration.md` is the only document that scores a **solved spatial
  field** through the shipped solver path. 3 of 6 criteria pass; the magnitude criterion
  fails, and the honest form of that failure is that the transport length λ is
  **2.45× too long** (`D_eff/k` ~6× too large).
* The pre-registered tissue-capacity test **FAILED** (`results/brain_slice_anchor.md`):
  the engine's implied `Vu,brain` is **1.000 for every compound**, below all 33
  measurements, ratio 12.05 against a pre-set bar of 3.0. `GMFE_floor = 4.77×` is the
  finding that survives any recalibration.

### 1.2 Confirmed dead — do not re-propose

Carried forward from `results/DATA_HUNT.md` §4, `results/external_anchor_survey.md` §4,
`results/kyle2004_penetration.md` §2, `results/zenodo_kibur_deposit_check.md` and
`docs/DATA_ACCESS_PLAN.md`, plus three I re-verified myself in this pass:

| ruled out | why |
|---|---|
| **Glioma + temozolomide spatial drug map** | Does not exist, open or gated. **Re-verified by me, 2026-09-20**, by direct query of the METASPACE GraphQL API: `countDatasets` = **20,417** total; `simpleQuery` returns **temozolomide 0, doxorubicin 0, paclitaxel 0, olaparib 0** and **imatinib 2** (the GIST set the repo already uses). The positive control proves the query works. This upgrades the repo's `[REPORTED]` zero-TMZ finding to first-hand |
| **Agar/White adavosertib GBM MALDI-MSI** | No accession; *"available upon request"*; the PSOC is defunct |
| **Jonas/Kibur glioma microdevice, as a machine-readable pairing** | `results/zenodo_kibur_deposit_check.md` opened both Zenodo files: one is bulk GSEA, the other a single-section endogenous-metabolite MSI raster with no drug channel and no IHC channel. *(The figures are a separate, live route — see R8. The deposit is dead; the figures are not.)* |
| **Pruijn 2008 HT-29, as validation** | Spent — it selected the transport law's structure |
| **Grantab 2006 (PMID 16424039), as a cross-drug penetration comparison** | Rejected by `results/kyle2004_penetration.md` §2 on detectability-driven dose asymmetry. **The rejection is correct and I am not overturning it.** *(Its within-drug packing contrast is a different experiment — see R4, which is built on the open 2012 sibling, not on this paper.)* |
| **Tunggal 1999 (PMID 10389947); Trédan 2009 (PMID 19176397)** | Rejected on substance: self-contradictory concentrations, no SD/SEM, dose asymmetry, inseparable parent/metabolite signal |
| **Brana 2014 (PMID 25213162)** | Unretrievable |
| **3-D vascular network graph of human GBM** | `results/landscape_pathology_3d.md` searched twice, adversarially, and confirmed the negative: HiP-CT (0 brain-tumour specimens, 0 vessel segmentations), IDR (0/147), Brain Image Library (0/14,033). *(A rat **tumour** network does exist and was not checked — see R10.)* |
| **HTAN, PDMR, dbGaP-gated TCIA** | Gates require eRA Commons + a tenure-track PI + an institutional Signing Official. Unavailable to an unaffiliated individual. HTAN has no GBM atlas; PDMR has zero high-grade glioma |
| **QIN-GBM-Treatment-Response** | Both scans are pre-therapy. Contains no response data despite the name |
| **Clinical DCE/DSC-MRI as a source field** | `results/landscape_imaging_derived.md`: 92–156× coarser per edge than the 12 µm grid, and no section we hold has paired MRI |

### 1.3 Three blockers that are now stale — this is where the new routes come from

**These are the most useful paragraphs in this document.** Three separate `results/`
files declined a class of test on capability grounds. All three capabilities now exist.
Nobody appears to have gone back.

1. **"The engine has no 1-D slab geometry."** — `results/external_anchor_survey.md`
   §4, declining Tunggal/Cowan-class multicellular-layer tests as *"COMMENSURABLE —
   PARTIAL. NOT ACTIONABLE without new code."* `results/anchor_registration.md` repeats
   it: *"Not actionable without new geometry, and the omission is expensive."*
   **Stale.** `results/kyle2004_penetration.md` (2026-08-22) solves a **140 µm slab
   through `simulate_response`** in its arm A, 28 slab solves including a grid sweep.
   The geometry exists.

2. **"The frozen public solver has no Dirichlet facility, so it structurally cannot
   represent a bath."** — `tmesim/validation/moncayo_explant.py:536-545`, and true of
   `physics/diffusion.py`. **Stale as a blocker on the *analysis* layer.**
   `moncayo_explant.py::bath_diffusion` (line 528) already implements a Dirichlet bath
   over an arbitrary boolean mask by importing the engine's own
   `_effective_diffusivity` and `_build_diffusion_operator` unchanged and adding only
   the boundary condition, with a local-equilibrium retardation factor, backward Euler
   and a cached LU. It is **verified against Crank's plane-sheet solution eq. 4.18 to a
   worst relative error of 0.53 %** (`verify_against_analytic_slab()`, line 608).
   Its signature is `bath_diffusion(mask, dx, drug, ecm, R, c_bath, total_time, ...)` —
   a boolean mask and a `DrugParams`. No transcriptomics, no cells, no vessels.
   **Every in-vitro transport measurement in the literature that uses a donor/receiver
   diffusion chamber is runnable against this today**, with the geometry registered by
   filesystem path like any other tissue.

3. **"The model's oxygen field is dimensionless and has never defined what `O = 1`
   corresponds to in mmHg — that missing mapping is the binding constraint."** —
   `results/oxygen_calibration_landscape.md`, final section. **Stale.**
   `tmesim/physics/oxygen.py` now declares `REFERENCE_PO2_MMHG` with `O = 1 ≡ 35 mmHg`
   in its docstring (lines 30-47) and ships a `to_mmhg` conversion. The mapping is a
   *declaration* rather than a calibration, which is a real caveat — but a declared
   mapping is scoreable and a missing one is not. **The measured human glioma pO₂
   literature that document catalogued is now comparable to a solved field.**

Separately: `results/CLAIM_INVENTORY.md` §6 item 0 (reconcile `vivum validate` with
`kill_calibration_breast.md` over zero-width bands) **is closed** —
`results/zero_width_bands.md` records the fix in both the tool and the generator. Do
not re-propose it.

---

## 2. The routes

### Class A — code verification. No external data required.

These are verification in the narrow V&V sense: *is the code solving the equations it
claims to solve?* — as distinct from validation, *are those the right equations?* For a
product positioned as a verification trust-layer this class is not a consolation prize;
it is the only class where the answer is unambiguous, cheap, and cannot be argued with.

---

#### R1 — Method of Manufactured Solutions on the shipped finite-volume operator

**Tests.** The **observed order of accuracy** of `physics/diffusion.py`'s discretisation
— the diffusion operator, the explicit upwind advection term, the reaction term, and
the no-flux treatment at the irregular tissue/non-tissue boundary — against its
*designed* order, on a problem whose exact solution is known by construction.

> **⚠ ALREADY IN PROGRESS — DO NOT DUPLICATE. Found in the working tree while this
> document was being written.** `experiments/mms_transport_verification.py` exists,
> untracked, 702 lines, mtime **2026-09-20 03:12**, and is owned by a concurrent
> workstream, not by this survey. Its module docstring states the same diagnosis this
> route reaches independently — *"A Richardson slope on a QoI tells you the answer is
> still moving with the mesh; it cannot tell you whether the PDE solver solves the PDE
> it claims to. Nothing in the tree verified that."* — and records the formal orders it
> will test against (2 in space for diffusion, 1 for donor-cell advection, 1 in time).
> **No `results/` document records an MMS outcome**, so the route below still describes
> an open question; what has changed is that the *implementation* effort is already
> spent by someone else. Read R1's "Effort" as **0 days of new work, plus whatever it
> costs to register the outcome**, and read its rank accordingly. Nothing in this
> document was written by, or should be attributed to, that workstream, and this file
> did not touch it.

**Engine predicts.** No prediction exists: no `results/` document records an observed
order of accuracy for the transport discretisation. What exists
instead is a set of grid-refinement symptoms that MMS is the standard instrument for
diagnosing: `results/grid_convergence_gbm.md` §0 reports the GBM ladder as
**non-monotone**, observed orders **0.96** and **−3.33**, GCI 5.44–6.42 %, and
power-law fits flagged **DEGENERATE**; `results/grid_convergence_breast.md` reports a
ladder spread of 49.7 % with docetaxel **−3.4 %** from its Richardson limit. Two claims
were retracted on the back of this (`CLAIM_INVENTORY` N2, N3), and the only verdict in
the whole product that is a function of the grid — breast docetaxel, currently
**UNDETERMINED at dx = 12 µm** — is downstream of it.

**Data.** None. The manufactured solution is chosen by the analyst; the source term is
whatever makes it exact. Standard reference: Salari K, Knupp P, *Code Verification by
the Method of Manufactured Solutions*, Sandia National Laboratories SAND2000-1444, June
2000 — **fetched**, https://www.osti.gov/servlets/purl/759450 (10,006,067 bytes). Its
abstract, verbatim: *"The strength of the MMS procedure is that it can identify any
coding mistake that affects the order-of-accuracy of the numerical method."*

**Effort.** 2–3 engineer-days. One day to build the manufactured source term and the
refinement driver; one to run the ladder on a rectangle and on a real irregular mask;
one to write up. The DCT Helmholtz solver in `validation/dimensionality.py:530` gives a
second, independent implementation of the same stencil for free.

**If it fails** (observed order below design on some term): **this is the single most
diagnostic outcome available anywhere in this document.** It would mean the
non-monotone ladders are a *discretisation defect*, not under-resolved geometry — and
it would name the term, because MMS is run term by term. The prime suspects are the
first-order upwind advection (which caps the whole scheme at first order wherever
convection is active, and convection is **ON by default**, `PDEConfig.convection =
True`) and the no-flux treatment at the irregular mask boundary. A failure here
retroactively qualifies **every** grid-refinement statement in `results/`.

**If it passes** (observed order matches design on every term, on both a rectangle and a
real mask): the licensed sentence is exactly — *"The shipped transport discretisation
attains its designed order of accuracy under the method of manufactured solutions on
the domains the product actually solves."* That sentence says **nothing** about whether
the physics is right, and must never be quoted next to one that does. It does mean the
non-monotone ladders are a geometry/resolution problem, which redirects the docetaxel
verdict work to the grid rather than the solver.

---

#### R2 — Extend the analytic-solution battery beyond Crank 4.18

**Tests.** Agreement of the solver against closed-form solutions on problems it can
already express: the 1-D slab with first-order reaction (**this is the exact solution
Hicks and Wilson fit their MCL flux data to**, so it is the analytic backbone of R4–R6);
the steady Krogh cylinder; the 2-D line source `K₀(r/L)` and 3-D point source
`exp(−r/L)/r`.

**Engine predicts.** Partly precedented and partly untested.
`moncayo_explant.py::verify_against_analytic_slab()` already checks the bath geometry
against Crank eq. 4.18 at **0.53 %** worst relative error.
`dimensionality.py::analytic_comparison` (line 280) already establishes that
`L = sqrt(D/k)` is identical in 2-D and 3-D and that only the prefactor differs, and
`control_transverse_is_exact` (line 643) agrees to ~1e-16. **What is missing is the
slab-with-reaction case** — precisely the one every MCL route needs.

**Data.** None.

**Effort.** 1–2 engineer-days, and about half of it is already written.

**If it fails.** Same diagnostic value as R1, narrower scope, lower cost. If the
slab-with-reaction case does not reproduce, R4–R6 cannot be interpreted at all, because
any disagreement with the measurement would be unattributable between the physics and
the solve.

**If it passes.** *"The solver reproduces the closed-form solutions of the problems it
is about to be scored on, to a stated tolerance."* This is a **precondition** for R4–R6,
not an independent result, and should be reported as one.

---

#### R3 — Cross-code verification against BioFVM

**Tests.** Whether two independently written codes, given the same domain, the same
`D`, the same decay rate and the same source distribution, produce the same field.

**Engine predicts.** No prediction. This has never been done. The engine has exactly
**one** independent implementation of its own stencil — the DCT Helmholtz solver — and
that one shares the repository, the author and the conventions, so it is not
independent in the sense that matters.

**Data.** BioFVM, Ghaffarizadeh A, Friedman SH, Macklin P, *Bioinformatics*
2016;32(8):1256-8, PMID **26656933**, DOI 10.1093/bioinformatics/btv730 — abstract
**fetched**. The sentences relied on: *"Solutions are first-order accurate in time and
second-order accurate in space"* and *"available for download at
http://BioFVM.MathCancer.org and http://BioFVM.sf.net under the Apache License
(v2.0)."* Open licence, no gate.

**Effort.** 3–4 engineer-days, most of it build-and-marshal rather than science.
BioFVM is C++/OpenMP; the comparison has to be driven through a file interchange.

**If it fails** (the two fields disagree beyond the combined discretisation error):
a genuine finding, and one that R1 alone could miss — MMS catches order-of-accuracy
mistakes, not sign conventions, unit slips or boundary-condition mismatches that are
*consistently* wrong. `results/hbd_convention_audit.md` and
`docs/VERIFICATION_METHODOLOGY.md` §10b (*"A consistently-applied wrong convention
validates itself"*) are exactly this failure class, and this repository has shipped it
before.

**If it passes.** *"The shipped transport solve agrees with an independent open-source
diffusive transport solver on a common benchmark, to a stated tolerance."* Modest, but
it is the kind of sentence a buyer of a verification layer actually wants, and no
competitor in this space has it either.

---

### Class B — in-vitro tissue transport. No vasculature, no clinical outcome, runnable on `bath_diffusion` today.

The unifying argument: the engine's transport law has been scored **as a coefficient**
(R² = 0.718 on 71 spent compounds) and, once, **as a solved field** (Kyle 2004, 3 of 6
criteria). These three routes score it as a solved field, in three different
laboratories, on three cell lines and three compound sets that the law was never fit to.
None of them requires a vessel, an outcome, or a patient.

---

#### R4 — Grantab & Tannock 2012: a measured packing-density perturbation ★

**Tests.** The **cell-volume-fraction term** `φ` in the two-route hindrance law. The
engine computes `hindrance = (1 + 2βφ)/(1 − βφ)` (Maxwell/Hasselman-Johnson effective
medium, `tmesim/molecular/two_route.py::two_route_factor`, lines 420-447) with
**φ = 0.7892**, a number derived from *this repository's own Xenium breast section* by
Delesse's principle (`two_route.py:125`) and **never scored against any external
measurement**. This route perturbs φ by a measured amount and asks whether the predicted
change in penetration matches the measured change.

**Engine predicts.** `hindrance` is an explicit, differentiable function of φ, so the
prediction is computable in closed form before any solve, and then again as a solved
field through `bath_diffusion` on a slab of the measured thickness. **The engine
currently predicts nothing for this experiment because it has never been asked.** Its
neighbouring verdicts: MCL absolute `D_eff` **FAIL** (49.33 µm²/s against a 74–130
µm²/s band), MCL hindrance ratio **PASS** (10.3274× against a measured 12×) — *two
anchors on the same physics that disagree about where the engine sits*, which is
`results/two_route_promotion.md`'s own unresolved finding and which this route is
positioned to adjudicate.

**Data. Open, fetched in full, numbers in hand.** Grantab RH, Tannock IF, *BMC Cancer*
2012;12:214, PMID **22672469**, PMC **3407510**, DOI 10.1186/1471-2407-12-214.
Fully open access. What I extracted directly from the article text:

* **Packing density, measured as % nuclear area, before and after a non-toxic
  perturbation**: EGTA 10 mM / 90 min reduced it *"from 50 ± 7.3% to 40 ± 1.6% and from
  40 ± 6.8% to 31 ± 4.5%, in MCC derived from HCT-8Ea and HCT-8E11 cells
  respectively (p = 2.2 × 10⁻⁹ and p = 8.1 × 10⁻⁶)"*, and the paper explicitly checks
  that the change is *"uniform throughout the MCC and not limited to the MCC surface
  (p = 0.78 and p = 0.82)"* — which matters, because a surface-only change would not be
  a bulk φ change.
* **Penetration change, same drug, same dose, same MCC**: *"EGTA pre-treatment improved
  penetration of 5-fluorouracil and gemcitabine (measured at 6 h) by ~3.3-4.4 fold in
  MCC derived from HCT-8Ea cells and by ~1.5-1.9 fold in MCC derived from HCT-8E11
  cells"*, and *"EGTA also led to ~2-fold increase in penetration of doxorubicin
  (p = 0.02 for HCT-8E11, p = 0.002 for HCT-8Ea)"*.
* **A second, mechanistically independent perturbation**: *"Bortezomib pre-treatment
  improved penetration of 5-fluorouracil by ~1.7-fold and gemcitabine by ~3-fold in MCC
  derived from HCT-8Ea cells and increased penetration of both drugs by ~2-fold in MCC
  derived from HCT-8E11 cells at 6 h."*
* **The observable is dimensionless and internally standardised**: penetration is
  reported as `C/C∞`, *"where C is the measured drug concentration and C∞ is the
  calculated drug concentration at equilibrium between the two drug compartments"*, with
  *"[14C]-sucrose … at a concentration of 3 μM in all experiments as an internal
  standard"*.

**Why this one survives the objection that killed Grantab 2006.**
`results/kyle2004_penetration.md` §2 rejected Grantab 2006 because its 15.4-fold
concentration spread was *"explicitly driven by detectability"* — a **cross-drug**
comparison at unequal doses. It then said, correctly, *"It is the right dataset for a
different experiment."* **This is that different experiment.** Every comparison here is
*the same drug, at the same dose, in the same cell line, with the same internal
standard, before and after a measured change in packing density*. Dose asymmetry
cancels identically out of a within-drug ratio. **Note also that Grantab 2006 itself is
not obtainable**: Unpaywall reports `"oa_status":"closed"`, `"has_repository_copy":
false`, `"oa_locations":[]` for DOI 10.1158/0008-5472.CAN-05-3077, and
aacrjournals.org returns a Cloudflare challenge. This route uses the 2012 paper only.

**Effort.** 3–4 engineer-days. A slab mask of the measured thickness, `bath_diffusion`
with a two-compartment donor/receiver boundary, φ perturbed by the measured nuclear-area
change under a **pre-registered** conversion from % nuclear area to cell volume fraction
(this conversion is the route's main methodological risk and must be fixed before any
residual is visible), predicted `C/C∞(6 h)` ratio scored against 3.3–4.4× / 1.5–1.9× /
~2× / ~1.7× / ~3×. Eight independent measured ratios across two drugs, two cell lines
and two perturbations.

**If it fails.** It names the term, and the term is one the product leans on hard.
φ = 0.7892 enters *every* compiled compound's `effective_D`. A failure in **magnitude**
says the effective-medium law has the wrong sensitivity to cell packing; a failure in
**direction** would be much worse and would falsify the law outright. Either way this is
the first external constraint of any kind on a number derived from the repo's own
tissue — and `docs/VERIFICATION_METHODOLOGY.md` §1 calls the compensating-error class
(*"the most dangerous defect class in scientific code"*) exactly what an unconstrained
φ invites, since φ and the membrane leg were fit and frozen together.

**If it passes.** *"The two-route hindrance law reproduces the measured change in
multicellular-layer drug penetration produced by a measured change in cellular packing
density, within a pre-declared band, on two drugs, two sublines and two independent
perturbations, in a laboratory and cell line the law was not built on."* That is the
strongest single sentence available to this product from any route in this document,
and it is a **transport-physics** sentence, not a response sentence.

---

#### R5 — Hong, Wilson & Hicks 2019: HCT116 MCL coefficients, with an extracellular-only marker ★

**Tests.** Two things at once, and the second is the reason to run it.
(i) Absolute `D_eff` for two compounds the law has never seen, in a cell line it has
never seen. (ii) **The extracellular leg in isolation.** The paper measures mannitol —
a pure paracellular marker with no transcellular route — in the *same* MCLs as the drug.
The engine's extracellular leg (Stokes–Einstein → Brownian micro-tortuosity →
Ogston–Renkin fibre gate → `D_eff = D·(1 − 0.6·ecm)` with **`ALPHA_TORTUOSITY = 0.6`
carrying no citation at all**, `physics/diffusion.py:75`) is currently constrained by
**nothing**. A marker with no membrane leg constrains it and only it.

**Engine predicts.** For a compound with no membrane permeability the two-route factor
reduces **exactly** to `(1 − φ)/(1 + φ/2)`, so the mannitol prediction is a closed-form
consequence of φ and the ECM tortuosity constant, with zero free parameters. The engine
has never emitted it.

**Data. Open, fetched in full, numbers in hand.** Hong CR, Wilson WR, Hicks KO,
*Neoplasia* 2019;21(2):159-…, PMID **30591421**, PMC **6314220**, DOI
10.1016/j.neo.2018.11.009. Quoting the article text directly:

* *"tissue (MCL) diffusion coefficient in the absence of metabolism, giving
  **D = (1.82 ± 0.07) x 10⁻⁷ cm² s⁻¹**"* for TH-302 (evofosfamide) in HCT116 MCLs
  — i.e. **18.2 µm²/s**.
* *"D for IPM = **(1.33 ± 0.04) x 10⁻⁷ cm²/s**, similar to the extracellular marker
  mannitol [47] in HCT116 MCLs (**D = (1.46 ± 0.02) x 10⁻⁷ cm²/s**) and thus consistent
  with paracellular (extracellular) diffusion of IPM."* — **13.3 and 14.6 µm²/s**.
* *"The first order rate constant for TH-302 metabolism was determined from transport
  through anoxic HCT116 MCLs (**kmet,0 = 0.0115 ± 0.05 s⁻¹**)"*.
* MCL thickness *"**139 ± 2 μm** based on 14C-urea flux"*, three MCLs per condition,
  with full time-courses in donor and receiver compartments (Figures 2 and 3) and
  *"(All model parameters are listed in Table S1)"*.

**This also closes an open item the repo has carried since 2026-08-10.**
`results/anchor_registration.md` lists as still open: *"A measured cell (or
extracellular) volume fraction for a Wilson/Hicks multicellular layer… it is still
open."* Mannitol and ¹⁴C-urea in the same MCLs are that measurement, read through an
Ogston/Mackie–Meares relation — which the engine already implements
(`molecular/ecm_mesh.ogston_renkin_hindrance`, lines 153-190).

**Effort.** 2–3 engineer-days for the coefficient arm (it slots straight into the
existing `mcl-transport-coefficients` group); +2 for the solved-field arm scoring the
receiver-compartment time-course through `bath_diffusion`.

**If it fails.** The two-anchor contradiction the repo already owns —
absolute `D_eff` FAIL, hindrance ratio PASS — gets a third, independent reading in a
third cell line, and the mannitol row **attributes** it: if mannitol is right and TH-302
is wrong, the membrane leg is at fault; if mannitol is wrong, the extracellular leg and
`ALPHA_TORTUOSITY = 0.6` are. **The current anchor set cannot make that attribution at
all.** That is the whole value of this route.

**If it passes.** *"The engine's extracellular transport leg reproduces the measured
multicellular-layer diffusion coefficient of a pure paracellular marker, and its full
law reproduces the measured coefficient of an out-of-sample compound in an
out-of-sample cell line, within a pre-declared band."* Note the two halves are
separately reportable and must be — a pass on one and a fail on the other is the most
informative outcome and must not be summarised as "mixed".

---

#### R6 — Phillips 1998: a fourth laboratory, and a *time-course*, not a coefficient

**Tests.** The solved field, as a function of time, in an independent laboratory
(Bradford, not Auckland, not Toronto) and a fourth cell line (DLD-1). The observable is
the **half-time of appearance in the receiver compartment**, which depends on `D`, on
metabolic consumption, and on the slab thickness jointly — so it is a stiffer test than
a fitted coefficient.

**Engine predicts.** Nothing yet. Its nearest emitted quantity is the analytic
penetration length `L = sqrt(D_eff/k_eff)` and the fitted `alpha`.

**Data. Open, fetched.** Phillips RM, Loadman PM, Cronin BP, *Br J Cancer*
1998;77(12):2112-9, PMID **9649122**, PMC **2150429**, DOI 10.1038/bjc.1998.355. From
the abstract, verbatim: *"Both … tirapazamine and … NITP rapidly penetrated DLD-1
multicell layers (**50.9 +/- 12.1 microm thick**) with **t(1/2) values of 1.36 and 2.38
h** respectively, whereas the rate of penetration of … EO9 and doxorubicin through
multicell layers was significantly slower (**t(1/2) = 4.62 and 13.1 h**)."*

**A caveat that must be stated in the same breath.** Doxorubicin is the compound
`alpha` is fit to (Primeau 2005), so **the doxorubicin row here is in-sample and must
be registered NOT TESTABLE**, exactly as `scorecard.py::_run_perivascular_penetration`
already declines doxorubicin's own penetration anchor. Three of four rows are scoreable;
tirapazamine additionally overlaps the Hicks anchor set, so the two genuinely
independent rows are **NITP and EO9**. Small n, stated up front.

**Effort.** 2 engineer-days, once R2's slab-with-reaction case is verified.

**If it fails.** Lower diagnostic yield than R4 or R5 — the failure would be
attributable to `D`, to `k`, or to the thickness estimate, and the paper reports a
±24 % thickness spread. Its value is **replication across laboratories**, not
attribution. A model that matched Auckland and missed Bradford would be describing a
laboratory, not tissue.

**If it passes.** *"Predicted multicellular-layer breakthrough half-times agree with
measurement on two out-of-sample compounds in a fourth laboratory."* n = 2. Say so.

---

### Class C — transport in real tissue

---

#### R7 — Extracellular-space diffusion in human glioma, by real-time iontophoresis

**Tests.** `ALPHA_TORTUOSITY = 0.6` and the extracellular leg, **in the actual flagship
tissue** — human high-grade astrocytoma — using a probe that has no transcellular route
and therefore isolates the geometric hindrance of the extracellular space.

**Engine predicts.** `D_eff = D·(1 − 0.6·ecm)`, i.e. a hindrance factor of `1/(1 − 0.6·ecm)`
with a floor at `max(1e-6, 1e-3·D)` (`physics/diffusion.py:78-88`). The measured
quantity is `λ² = D_free/D_apparent`, dimensionally the same object.
`results/uncertainty.md` already notes the correspondence and calls 0.6 *"a MODELLING
CONSTANT, never fit to this section"* — which is the point.

**Data. Located and abstract fetched; the numbers I actually need are not.** Vargová L,
Homola A, Zámecník J, Tichý M, Beneš V, Syková E, *Glia* 2003;42(1):77-88, PMID
**12594739**, DOI 10.1002/glia.10204. From the fetched abstract, verbatim: *"ECS volume
fraction alpha … tortuosity lambda (lambda(2) = free/apparent diffusion coefficient) and
nonspecific uptake k' were studied by the real-time tetramethylammonium method in acute
slices of human tissue… The average values of alpha, lambda, and k' in control cortex
were **0.24, 1.55, and 3.66 x 10⁻³ s⁻¹**, respectively… **In cellular regions of
high-grade astrocytomas (WHO grade III and IV), alpha and lambda were further
increased, and k' was significantly larger than in controls.**"*

**The honest obtainability statement.** The **control** values are in the abstract and
are usable now. The **grade III/IV** values — which are the ones that matter, because
the flagship is GBM — are in the paper's tables, and *Glia* 2003 is not open access. I
did not obtain them in this pass. **This route is therefore PARTIALLY OBTAINABLE**, and
it must not be ranked as though the GBM numbers were in hand. The cheap next step is a
single retrieval attempt (Internet Archive / interlibrary), which the repo has
successfully used before — `results/diffusivity_falsification.md` recovered Pruijn 2008
that way after diagnosing its 403 as *"a Cloudflare bot challenge… not a paywall"*.

**Effort.** 0.5 day to attempt retrieval; 1–2 days to score if it succeeds. If it does
not succeed, **stop** — do not substitute a review's summary for the primary table.

**If it fails.** Strongly diagnostic. `ALPHA_TORTUOSITY` and the MCL-derived hindrance
are the only two constraints on the extracellular leg, and one is uncited. A third
reading, in the right tissue and in situ, would either reconcile the existing
disagreement (MCL 12× vs the engine's 10.33×, against a brain-ECS λ² that is likely
nearer 3–5×) or prove that **the MCL hindrance and the in-tissue hindrance are simply
different numbers** — which would be a finding about the whole MCL-anchoring strategy,
not just about this engine.

**If it passes.** *"The engine's extracellular hindrance is consistent with
tortuosity measured in situ in human high-grade astrocytoma, within a pre-declared
band."* One paper, one method, one laboratory.

---

#### R8 — In-human glioma drug concentration versus distance, digitised from a CC-BY figure

**Tests.** The **transport length** of a drug in living human glioma — the one quantity
the repo has three mutually inconsistent readings of and no arbiter.

**Engine predicts.** `results/kyle2004_penetration.md` (corrections of 2026-08-22)
records the state of the contradiction precisely: the engine's λ is **2.45× too long**
(range 2.00–2.98×), `D_eff/k` ~6× too large, and the two existing external readings are
**disjoint** — *"Primeau implies λ_dox 57.71–72.13 µm, Kyle implies 30.73–55.68 µm."*
`alpha` is fitted to the Primeau end of that disjoint pair. A third, independent,
**in-human, in-vivo** reading is the only thing that adjudicates.

**Data. The record is open; the numbers are in a figure.** Peruzzi P, Jonas O, et al.,
*Sci Transl Med* 2023, PMID **37672566**, PMC **10754230**, DOI 10.1126/scitranslmed.adi0069
— **not** open access. But `results/zenodo_kibur_deposit_check.md` established, by
reading all 34 pages, that a **CC-BY-4.0 preprint of the same trial, same 6 patients,
same NCT04135807** is fully readable at Research Square DOI 10.21203/rs.3.rs-1647435/v2,
mirrored at https://kiburmed.com/wp-content/uploads/2022/07/Peruzzi-et-al-submitted-to-Nat-Med.pdf,
and that its Fig. 2a/c is *"a **line plot** of local drug dose (µM) vs. distance (µm),
one colored line per patient"* for **doxorubicin and lapatinib**. That same document
lists digitisation as its own route 2 and judges it *"sufficient for a first-pass
distance-binned comparison against a transport model."*

**The norm this route has to satisfy.** `results/DATA_HUNT.md` §4 says, of a paper whose
concentrations existed only as five binned colour ranges, *"A figure is not a dataset."*
That rule is right and this route does not violate it, but it only survives under three
conditions, all of which must be written into the pre-registration:
1. **Score a shape parameter only** — the exponential decay length fitted to the
   digitised profile. Never an absolute µM. The repo's own rule that ion intensity is
   not concentration applies with equal force to a digitised axis.
2. **Publish the digitisation error**, estimated by re-digitising blind and by
   round-tripping a known axis tick.
3. **Label every downstream number as digitised**, in the row's `note`, permanently.

**The in-sample trap, stated before the run.** Doxorubicin is the compound `alpha` is
fit to. The **doxorubicin arm of this route is in-sample and must be reported NOT
TESTABLE**, exactly as the Primeau row already is. **Lapatinib is the out-of-sample
arm** and is the whole test. If the lapatinib profile cannot be digitised cleanly, this
route returns nothing and should be abandoned rather than rescued with doxorubicin.

**Effort.** 2–3 engineer-days. Digitisation is hours; the pre-registration and the error
analysis are the work.

**If it fails.** It would say the engine's transport length is wrong in living human
brain tumour — which is the claim that matters most and is currently supported by
xenograft and in-vitro measurements only. It would also, independently of the engine,
tell you which of Primeau and Kyle the in-human value sits nearer, which is a
publishable observation about the *anchors* regardless of what the model does.

**If it passes.** *"The engine's predicted transport length for an out-of-sample
compound is consistent with a digitised, in-human, in-vivo concentration-versus-distance
profile, to a stated digitisation error, n = 6 patients, one compound."* Every clause of
that sentence is load-bearing.

---

### Class D — the solved oxygen field

---

#### R9 — Eppendorf pO₂ histography in human glioma versus the solved oxygen field ★

**Tests.** The solved oxygen field, **as a field**, converted to mmHg through the
declared reference — against measured intraoperative pO₂ distributions in human glioma.
This is a PDE-scoped row, and the repository currently has almost none that are both
scored and discriminating.

**Engine predicts.** `solve_oxygen` (`physics/oxygen.py:417`) solves
`D∇²O + supply·vasc·(1 − O) − consumption·ρ_tum·O = 0` with **`supply = 0.05` and
`consumption = 0.02`, both explicitly uncalibrated** (`core/config.py:145-151`), and
`O = 1 ≡ 35 mmHg` declared. The hypoxia band is `O < 0.2`, i.e. **< 7 mmHg**. The engine
currently reports **GBM HYPOXIC = 0.0000 at every grid spacing**
(`results/CLAIM_INVENTORY.md` C5) — *zero* hypoxic cells.

**Data. Open, three independent studies, all fetched.**
* Lally BE, Rockwell S, Fischer DB, Collingridge DR, Piepmeier JM, Knisely JP,
  *Cancer J* 2006;12(6):461-6, PMID **17207315**: *"Tissue oxygen tension was measured
  intraoperatively using an Eppendorf pO2 Histograph… we included 13 patients with high
  grade gliomas (HGG) and 10 with low grade gliomas (LGG). **The median tumor oxygen
  pressure for the entire cohort was 5.1 mmHg.**"*
* Collingridge DR, Piepmeier JM, Rockwell S, Knisely JP, *Radiother Oncol*
  1999;53(2):127-31, PMID **10665789**: *"**Both glioma and peritumoural brain contain
  regions at oxygen tensions less than 2.5 mmHg.**"*
* Rampling R, Cruickshank G, Lewis AD, Fitzsimmons SA, Workman P, *Int J Radiat Oncol
  Biol Phys* 1994;29(3):427-31, PMID **8005794**: *"**Six different tracks, each with a
  path length of 22 mm, were recorded per patient representing 192 readings**… Both the
  high grade intrinsic and metastatic brain tumors showed **significant regions of
  hypoxia**."*

**The measurement is a distribution, not a point, and that is the point.** Eppendorf
histography returns a histogram of readings along a needle track. The right comparison
is **distribution to distribution** — the engine's per-cell pO₂ distribution over a
section against the measured per-reading distribution — scored on the median and on the
fraction below 2.5 and 10 mmHg, with the inter-laboratory variability floor that
`results/oxygen_calibration_landscape.md` already documented (Nozue 1997, PMID 9290690:
median tumour pO₂ varied 1.5–5.6 mmHg for the *same* technique on the *same* model
across six laboratories) as the noise floor the bar must clear. A model predicting a
hypoxic fraction of exactly zero fails against a measurement whose *median* sits below
the model's own hypoxia threshold, and it fails by a margin no laboratory variability
can absorb.

**Effort.** 2–3 engineer-days. The solver, the conversion and the section are all in
place; the work is the pre-registered distributional statistic and the band.

**If it fails.** Very likely, and **that is not a reason to skip it — it is the reason
to run it.** Three things follow immediately and none of them requires the test to pass.
(i) It converts the oxygen story from a *coefficient* failure (supply 4,570× low,
compensated by a sink 112–166× low — a coincidence `CLAIM_INVENTORY` §6 item 5 already
flags as one the model should not rely on) into a *field* failure, which is a different
and more serious statement. (ii) It tests the compensation directly: two errors that
cancel in the mean will not cancel in the tails, so the *shape* of the predicted
histogram is the discriminator between "mis-scaled" and "structurally wrong". (iii)
`results/oxygen_calibration_landscape.md` already identified the structural suspect —
the engine's sink is **first-order in O everywhere** (`−c·ρ_tum·O`) where the literature
treats oxygen consumption as effectively **zeroth-order** until near-anoxia, and a
first-order-everywhere sink *"will systematically under-predict hypoxia"*. **A zero
hypoxic fraction is precisely the predicted signature of that defect.** A FAIL here
would be the first measurement confirming a mechanism the repo derived on paper.

**If it passes.** *"The solved oxygen field reproduces the measured distribution of
intraoperative tissue pO₂ in human glioma, on the median and on the sub-2.5 mmHg
fraction, within a band set by the published inter-laboratory variability."* This would
be a surprise, and if it happens the first thing to check is whether the declared
`O = 1 ≡ 35 mmHg` is doing the work rather than the physics.

---

### Class E — vascular geometry and stereology

---

#### R10 — A real 3-D tumour microvascular network, and the 2-D sectioning bias ★

**Tests.** Whether the engine's two worst geometry failures are a *segmentation* error
or a *stereology* error. The reconstruction reports **MVD 233.27 /mm²** against a
measured 65 (35–95) /mm² — **3.59×** — and **vascular area fraction 10.982 %** against
Marien 2016's commensurable 4.70 ± 2.38 % — **2.34×**
(`results/vascular_segmentation_audit.md`, `results/heldout_patient.md`). Both measured
anchors are **2-D histological counts on thin sections**. The reconstruction is also
2-D. Nobody has computed what a 2-D section of a *known* 3-D tumour vasculature reports,
relative to the truth.

**Engine predicts.** `results/vessel_geometry_fix.md` §§145-311 already worries about
exactly this — it invokes Delesse's principle and the stereological relation
`L_V = n_A/(t·λ)`, and notes that using one dataset for both length density and volume
fraction is a weakness. **It has never had a 3-D tumour network to test the relation
against.** It used H01 (normal human cortex from an epilepsy resection), which
`results/external_anchor_survey.md` §0.1 separately flags as the wrong tissue.

**Data. Open, downloadable, and I verified it.** Timothy Secomb's group publishes the
mapped 3-D geometry of a rat tumour microvascular network:
`https://sites.arizona.edu/secomb/files/2024/01/rattum98-ffab1d1a873244ec.txt` —
**fetched 2026-09-20, HTTP 200, 8,049 bytes, 106 lines.** Header, verbatim from the
file: `Segment / Node Node / diam. length / node from x y z / node to x y z`. Page text:
*"R3230AC mammary carcinoma in rat dorsal skin flap preparation"*, region dimensions
*"550 x 520 x 230 microns"*, all values in microns. Citation on the page: T.W. Secomb,
R. Hsu, R.D. Braun, J.R. Ross, J.F. Gross, M.W. Dewhirst, *"Theoretical simulation of
oxygen transport to tumors by three-dimensional networks of microvessels"*, in *Oxygen
Transport to Tissue XX*, Plenum, New York, 1998, pp. 629-634.

**Why this was not found before.** `results/landscape_pathology_3d.md` searched hard,
twice, adversarially, for a 3-D vascular graph — and searched HiP-CT, IDR, the Brain
Image Library, VesselGraph, Bernabeu and Stolz. **It searched for a human GBM network,
and correctly concluded none exists.** It did not check a physiologist's lab page for a
*rat tumour* network with node/edge/diameter geometry. The species and model are the
honest limitation; they do not affect the **stereological** question, which is a
geometry question about sectioning, not a biology question about rats.

**What to compute.** Three things, all cheap and all exact.
(i) Take repeated virtual 2-D sections of the known 3-D network at the section
thicknesses the anchors used, and compute the *measured* MVD and areal fraction a
histologist would report, against the *true* 3-D length density and volume fraction.
That ratio **is the sectioning bias**, computed rather than assumed.
(ii) Compute the true 3-D distance-to-nearest-vessel distribution — the length scale
that sets everything downstream in a diffusion model — and compare its shape to the
reconstruction's.
(iii) Independently, a second use: the network is a legitimate *source geometry* for a
delivery solve, which is what Hicks 2006 (PMID **16912264**, abstract fetched) used a
mapped network for, reporting *"a strong correlation between model-predicted and
measured hypoxic cell kill (R² = 0.89) but a poor correlation when the model did not
include extravascular transport (R² = 0.32)"*. **Do not pursue (iii) as a kill test** —
that crosses into response. Pursue it only as a geometry sanity check.

**Effort.** 1–2 engineer-days. An 8 KB text file, a parser, and a sectioning loop.
**This is the cheapest substantive route in this document.**

**If it fails** — i.e. the computed sectioning bias does *not* account for the 2.34×/3.59×:
the errors are real segmentation errors and the reconstruction over-counts vessels. That
closes a question the repo has been unable to settle and removes a live excuse.

**If it "passes"** — i.e. the bias does account for a substantial part of them: **this is
a genuine recovery of a result the project currently scores as a failure**, and it is
exactly the kind of thing the brief warns against missing. It would not make the anchors
pass; it would mean the anchor comparison is **not commensurable**, which is the same
defect the repository already identified and fixed once for the Tamma CD34 row
(`vascular_segmentation_audit.md` §5 — where the incommensurable anchor was left in
place and `vivum validate` still prints 5.23×). A second instance of the same defect, in
the same anchor group, found by computation rather than by argument, would be a strong
result. **Report the bias with an uncertainty interval over sectioning planes**, and do
not quote a single number.

---

### Class F — human tumour interstitial pharmacokinetics

---

#### R11 — Extend the microdialysis anchor beyond n = 2, and test the time-resolved path

**Tests.** `breast-transvascular` — one of only **three** scored rows in the whole
product where a solved field meets an external measurement. It currently has **two
compounds**, and both FAIL, in opposite directions.

**Engine predicts.** 5-FU **7.7956** against a measured 0.61 ± 0.11 (**12.8× high**);
methotrexate **0.0571** against 0.60 ± 0.20 (**10.5× low**)
(`results/breast_kp_validation.md`). `results/breast_transport_anchors.md` §6 names the
two structural causes — *"a blood-**brain**-barrier permeability applied to breast
tissue, and a one-way delivery source with no back-flux"* — and the second of those is
the `transvascular_exchange = False` default (`core/config.py:137`) that the codebase's
own comment calls *"measurably wrong"*.

**Data. Open, two studies, both fetched.**
* **Carboplatin.** Konings IR, Engels FK, Sleijfer S, Verweij J, Wiemer EA, Loos WJ,
  *Cancer Chemother Pharmacol* 2009;64(3):509-16, PMID **19096847**, PMC **2691802**:
  *"Plasma and microdialysate samples from tumor and adipose normal tissues were
  collected up to 47 h after dosing in **eight** carboplatin-treated patients with an
  accessible (sub)cutaneous tumor… Pharmacokinetics were evaluable in tumor tissue in
  **6/8** patients… with **exposure ratios of tissue versus plasma ranging from 0.64 to
  1.46**."*
* **5-FU, time-resolved.** Konings IR, Sleijfer S, Mathijssen RH, et al., *Cancer
  Chemother Pharmacol* 2011;67(5):1055-62, PMID **20652702**, PMC **3082021**: *"In
  **seven** patients … treated with a continuous 5-fluorouracil infusion, plasma and
  microdialysate samples from tumor and normal adipose tissue were collected over a
  period of **5 days**… Concentration-time curves of unbound 5-fluorouracil were lower
  in tumor tissue compared to the curves in plasma, but **exposure ratios of tumor
  tissue versus plasma increased during the 5-day infusion period.**"*

**Two distinct tests, and the second is the more interesting.**
(a) **Carboplatin** adds a platinum compound and a *normal-tissue paired control* to an
anchor group that has n = 2. It also adds a measurement whose ratio brackets 1.0, which
matters because the engine's one-way source has *"a steady state unbounded above
plasma"* and puts five of eleven GBM drugs above plasma — a measurement centred near 1.0
is a sharper discriminator than one centred near 0.6.
(b) **5-FU over 5 days** is the only external measurement located anywhere in this pass
that can test `PDEConfig.time_resolved`, which is **OFF by default** (`config.py:87`)
and therefore untested by every anchor in the registry. The measured signature — a
tumour:plasma ratio that **rises** during a constant infusion — is a dynamic,
dimensionless, shape observable that no static partition coefficient can express.

**The caveat that must be stated.** Both cohorts are *"(sub)cutaneous"* tumours of mixed
type, not breast primaries. They are **not** substitutes for the Müller breast rows and
must be registered as their own rows with their own tissue field, not pooled —
`results/external_anchor_survey.md` §0 records that this project *"manufactured a false
regression twice"* by pooling across modality and tissue.

**Effort.** 3 days for (a). 5–6 days for (b), because the time-resolved path has to be
exercised and the plasma driver has to be built from the reported infusion schedule.

**If it fails.** (a) tells you whether the 12.8×/10.5× straddle is a two-compound
accident or a general property — which is currently unknown and is the difference
between "two bad rows" and "the delivery law is wrong". (b) tells you whether the engine
can produce a *rising* tissue:plasma ratio at all. It probably cannot: a linear,
time-invariant transport model with a constant source and constant clearance approaches
a constant ratio. **If the reachable set of the model excludes a rising ratio, that is a
reachability refutation in the sense of `docs/VERIFICATION_METHODOLOGY.md` §8 — derivable
analytically, costing almost nothing, and unarguable.** Derive it before running the
solver.

**If it passes.** *"The engine's predicted tumour interstitial:plasma exposure ratio
agrees with human microdialysis measurement for an additional compound, in
subcutaneous tumour, n = 6."* And separately, for (b), *"the engine reproduces the
measured direction of change in that ratio over a 5-day infusion."*

---

### Class G — internal generalisation

---

#### R12 — Replicate the defensible claims on the sections already on disk ★

**Tests.** Whether the claims that currently survive attack survive it **on a second and
third section**. This is not an external anchor and must never be presented as one. It
is the item `results/CLAIM_INVENTORY.md` §6 ranks **first** among everything that would
change the honest assessment: *"**A second section, same tumour type.** Almost every
defensible claim above is n = 1 in tissue. This is the cheapest, largest improvement
available."*

**Engine predicts.** Every C-row and S-row in `CLAIM_INVENTORY` §1a/§1c is n = 1 in
tissue. The spatial-spread finding is *"one section, one drug, one strength"*, with
p at the 20-draw floor and no confidence interval.

**Data. Already on disk, registered by filesystem path, zero acquisition cost.** Per
`docs/DATA_INVENTORY.md` and this session's re-inspection: `data_cache/visium_gbm/`
(**4 spot-level sections**: the 10x parent plus GSE194329 GBM1/2/3),
`data_cache/greenwald_gbm_visium/` (**2 sections**, `GBM_ZH1007inf` and
`GBM_ZH1007nec`), `data_cache/gse324268_xenium/` (**1 GBM section**, 851,509 cells, and
it **carries MGMT**), `data_cache/gbmspace_visium/` (**30 sample directories**, 1.1 GB,
ingestion currently blocked on an `.h5ad` reader), plus the two Xenium GBM sections and
the Janesick breast pair already in use.

**The specific error this route exists to correct.** Four of those Visium sections have
been treated as unusable because the spatial-spread statistic is **refused in code** on
spot data (`require_measured_cell_positions`). **That refusal is correct and applies to
that statistic.** It does not apply to section-level quantities — the tissue mask, the
vascular density field, the ECM field, grid convergence, the invariance measurements, or
any ratio computed over a section rather than over per-cell neighbourhoods. Declining
the whole section because one statistic is refused is the second of the two errors the
brief names, and it is currently costing four sections of replication for free.

**Effort.** 4–6 engineer-days for the four Visium sections and the two Greenwald
sections (no new loader needed). Add ~3 for an `.h5ad` reader, which unlocks 30 more.

**If it fails** — i.e. a claim that holds on section 1 does not hold on sections 2–6:
**this is the most valuable failure available**, because it applies to the claims the
product is actually sold on rather than to a coefficient. Two of this project's
retractions came from configuration changes; a section-to-section failure would be the
first from *tissue* variability, and it would be the honest bound on what one section
licenses.

**If it passes.** Every surviving claim's narrowest true form gains *"replicated on N
sections from M patients"*, and the n = 1 caveat — currently attached to nearly
everything — comes off the ones that earn it. **No new anchor, no new external
measurement, no new physics.** Purely an upgrade in the generality of what is already
established, which is exactly why it is cheap and exactly why it is undervalued.

---

### Class H — open, but weak. Listed so the ranking is complete.

---

#### R13 — Interstitial fluid pressure: the *shape*, not the magnitude

**Tests.** The **spatial profile** of the solved IFP field — specifically the distance
from the tissue margin over which pressure rises to its interior plateau. This is a
length, and it is independent of the magnitude that already failed.

**Engine predicts.** A screened-Poisson field whose decay scale is
`sqrt(K/(a + d))`. **The scorecard already computes this number** — `screening_length_um`
with a median for GBM and for breast, plus the fraction of tissue where it falls below
the mesh (`tmesim/validation/scorecard.py:2493-2540`, `tmesim/physics/pressure.py:173`).
Its resolution row already flags that the screening length is below every grid spacing.
I have not run it and am not quoting its value.

**Data. Open, fetched.** Boucher Y, Baxter LT, Jain RK, *Cancer Res* 1990;50(15):4478-84,
PMID **2369726**: *"The majority (**86%**) of the pressure profiles demonstrated a large
gradient in the periphery leading to a plateau of almost uniform pressure in the deeper
layers of the tumors. Within isolated tumors, pressures reached plateau values at a
distance of **0.2 to 1.1 mm from the surface**… These results demonstrate for the first
time that the IFP is elevated throughout the tumor and drops precipitously to normal
values in the tumor's periphery."* n = 9 + 6 + 7 tumours, micropipette servo-null.

**Effort.** 1 engineer-day. The field is already solved on every run and the length is
already computed.

**If it fails.** Almost certain, and informative anyway. An IFP boundary layer that is
orders of magnitude shorter than the measured 200–1100 µm would mean the pressure field
is not merely **mis-scaled by 41×** but **structurally wrong** — and since IFP drives
the advection term, which is **ON by default**, it would establish that convection in
this engine is inert for a reason that should worry someone. That is a materially
different statement from the existing magnitude FAIL, and it is one sentence away from
being computable.

**Why it is ranked in Class H despite all that.** Because the outcome is close to
foreseeable from numbers the repository already holds, so the *information gained* is
smaller than the diagnostic language above suggests — and because the anchor is rat
tissue-isolated and subcutaneous tumours, not human glioma. It should be run as a
half-day rider on something else, not as a headline.

---

#### R14 — CIVO / GSE235475: graded intratumoral exposure in intact human tumour

**Tests.** The radius of a drug-effect zone produced by a point microdose in intact
human tumour, against a predicted diffusion radius over the same interval.

**Engine predicts.** A point-source diffusion problem the engine can express directly
(`bath_diffusion` with a point Dirichlet, or the existing analytic point-source form in
`dimensionality.py::analytic_comparison`). No current prediction.

**Data. Open, confirmed by me via the NCBI eSummary API, 2026-09-20.** GEO
**GSE235475**, *"Trackable Intratumor Microdosing and Spatial Profiling Provide Early
Insights into Activity of Investigational Agents in the Intact Tumor Microenvironment"*,
*Homo sapiens*, released 2023-06-28. From the record's summary, verbatim: *"12 patients
with head and neck carcinoma (HNC). Patients scheduled for tumor resection received
percutaneous intratumor injections of subasumstat and vehicle control 1 to 4 days prior
to surgery, **resulting in spatially localized and graded regions of drug exposure
(~1000-2000 microns in diameter)**. Drug-exposed (**n=214**) and unexposed regions
(**n=140**) were compared by GeoMx Digital Spatial Profiler, with evaluation at single
cell resolution in a subset of these by CosMx."*

**Effort.** 4–5 engineer-days.

**Why it is weak, stated plainly.** **The exposure axis is a pharmacodynamic readout,
not a measured drug concentration.** The ~1000–2000 µm diameter is the size of the
region in which an *effect* was detectable, which depends on an unknown and
compound-specific detection threshold as well as on transport. Converting a detectable-
effect radius into a transport parameter requires assuming that threshold. That makes a
disagreement unattributable between transport, potency and assay sensitivity — and
`docs/VERIFICATION_METHODOLOGY.md`'s verdict vocabulary would make this **NOT TESTABLE**
rather than a soft result.

**If it fails.** Teaches little, for the reason above. **This route is included for
completeness and to record why it is not recommended**, not as a candidate.

**If it passes.** Teaches little for the same reason. A route whose pass and fail are
both uninformative is not worth running, and this document says so rather than padding
the list.

---

## 3. The ranking — information gained per unit effort

Effort in engineer-days. "Attribution" means: would a failure **name the wrong term**,
or only tell you that something is wrong?

| rank | route | effort | obtainable? | scores a *solved field*? | attributes a failure? | information / effort |
|---:|---|---:|---|---|---|---|
| **1** | **R4 — Grantab 2012 packing-density perturbation** | 3–4 | **verified open** (PMC3407510) | yes | **yes — names φ** | **highest** |
| **2** | **R10 — Secomb 3-D tumour network, sectioning bias** | **1–2** | **verified, 8 KB downloaded** | no (geometry) | **yes — separates segmentation from stereology** | **highest** |
| **3** | **R1 — Method of manufactured solutions** | **0 new** (implementation already in the tree, owned elsewhere — see the banner on R1) | **no data needed** | n/a (code verification) | **yes — term by term** | very high |
| **4** | **R12 — Replicate on sections already on disk** | 4–6 | **on disk** | yes | no (generality, not attribution) | very high |
| **5** | **R5 — Hong 2019 HCT116 MCL + mannitol** | 2–5 | **verified open** (PMC6314220) | both arms | **yes — extracellular vs membrane leg** | high |
| **6** | **R9 — Eppendorf pO₂ histography vs solved field** | 2–3 | **verified open** (3 PMIDs) | **yes** | yes — sink functional form | high |
| **7** | **R11a — carboplatin microdialysis** | 3 | **verified open** (PMC2691802) | **yes** | partial | medium-high |
| **8** | **R11b — 5-FU 5-day time-course** | 5–6 | **verified open** (PMC3082021) | **yes** | **yes — reachability, derivable first** | medium-high |
| **9** | **R2 — analytic battery extension** | 1–2 | **no data needed** | n/a | yes | medium (and a *precondition* for 5, 6) |
| **10** | **R6 — Phillips 1998 DLD-1 half-times** | 2 | **verified open** (PMC2150429) | yes | no — replication only | medium |
| **11** | **R8 — digitised in-human profile** | 2–3 | **partial** — CC-BY figure, not a table | yes | yes, but n = 6 and digitised | medium |
| **12** | **R3 — BioFVM cross-code** | 3–4 | **verified open** (Apache 2.0) | n/a | yes — convention errors | medium |
| **13** | **R13 — IFP shape** | **1** | **verified open** (PMID 2369726) | yes | yes, but outcome near-foreseeable | low-medium |
| **14** | **R7 — human glioma ECS tortuosity** | 0.5 + 1–2 | **partial** — GBM values paywalled | no (coefficient) | yes | **conditional** — high if retrieved, zero if not |
| — | **R14 — CIVO GSE235475** | 4–5 | verified open | no | **no** | **not recommended** |

**Three notes on the ranking, because a table hides them.**

* **R2 is a dependency, not a competitor.** R5 and R6 cannot be interpreted until the
  slab-with-reaction analytic case is verified. Budget it inside them.
* **R7's rank is conditional by construction.** Half a day decides whether it is a
  strong route or no route. Spend the half day before ranking it against anything.
* **R10 is ranked second on effort, not on depth.** It is 1–2 days and it can recover a
  result currently scored as a failure. Nothing else on this list has that ratio.

---

## 4. The one I would run first, and why

> ### **R4 — the Grantab & Tannock 2012 packing-density perturbation.**

Five reasons, and one concession.

**(i) It is the only route that constrains a parameter the product cannot function
without and that nothing currently constrains.** `φ = 0.7892` enters the compiled
`effective_D` of *every* compound the engine has ever run, and it was derived from this
repository's own Xenium breast section by Delesse's principle — **a number the product
measured on its own tissue and then used as physics everywhere else.** The brain-slice
test already established that the engine has no per-compound tissue-capacity term. φ is
the *geometric* analogue of the same gap, and it is the one that is still testable.

**(ii) It is a within-drug ratio, which is the only class of claim this project has
successfully defended.** Every absolute has failed — pressure, oxygen, diffusivity,
partition, morphometry. Every survivor has been a ratio or an invariance. R4's
observable is a ratio of `C/C∞` at fixed drug, fixed dose, fixed line, fixed internal
standard. Dose asymmetry, detection efficiency, and absolute calibration cancel
identically. **This route is shaped like the things this engine has been able to get
right**, which is a reason to expect it to be informative rather than merely damning.

**(iii) It has eight independent measured ratios across two mechanisms.** 5-FU and
gemcitabine and doxorubicin, two sublines, EGTA and bortezomib. The bortezomib arm is
mechanistically independent of the EGTA arm — one is a calcium-dependent adhesion
disruption, the other proteasome inhibition — so a model that matched one and missed the
other would be caught. Most anchors in this repository are n = 1 to n = 6 with a single
mechanism.

**(iv) The data is open, fetched, and the numbers are already in this document.** No
gate, no MTA, no request, no institutional affiliation, no digitisation. PMC3407510,
fully open, and §2/R4 above quotes the packing densities, the penetration fold-changes,
the statistical significance and the definition of the observable verbatim. The 2006
sibling — which *is* paywalled, which the repo *did* correctly reject, and which is the
paper everyone reaches for — is not needed.

**(v) The capability to run it already exists and was declared missing.**
`bath_diffusion` gives the donor/receiver Dirichlet geometry, verified against Crank to
0.53 %. The repository declined this entire class of test on the grounds that it had no
slab geometry; it has had one since 2026-08-22. **Running R4 is partly an act of
noticing**, which is the cheapest kind of progress available.

**The concession.** R4 tests a *coefficient's sensitivity*, exercised through a solved
slab — not a solved field in patient tissue. It will not move the product's central
weakness, which is that only three scored rows compare a solved field to an external
measurement and all three fail. **R9 and R11 are the routes that attack that weakness
directly, and both should follow.** But R4 is the one that produces the most defensible
new sentence per day spent, and it produces it in the class of claim — transport
physics, scored as a ratio — that this project has actually earned the right to make.

**Run order I would actually commit to:** **R10** (1–2 days, may recover a result) →
**R2's slab-with-reaction case** (half a day, a dependency) → **R4** (3–4 days, the
headline) → **R9** and **R12** in parallel. **R1 is not in this order because it is
already being built by someone else** (see the banner on R1); the only thing to do
there is register its outcome when it lands, and it unblocks every grid statement in
`results/` when it does. The sequence above is roughly two-and-a-half weeks for one
engineer and it would change five rows of `CLAIM_INVENTORY`.

---

## 5. What cannot be verified without clinical data

Short, because it is already the best-documented part of this repository, and
proportionate, because it is not the point of this document.

1. **Response.** Predicted kill, coverage, or any ordering of therapeutic benefit cannot
   be scored against anything. No open dataset anywhere pairs a measured spatial drug
   distribution with a measured treatment response in the same tumour
   (`results/DATA_HUNT.md`, 39 agents, 164 datasets). This is not a gap this engine
   could close with better physics; it is the same categorical gap that keeps phase 1–3
   trials necessary after mature computational chemistry.
2. **Absolute kill fraction.** Withdrawn, permanently, and correctly.
3. **Cross-compound ordering.** Withdrawn 2026-08-21, absolutely. **None of the 14
   routes above creates a path back to it**, and a PASS on all 14 would not.
4. **The delivery → kill link, in space.** You can validate spatial delivery (ex vivo,
   non-glioma, small n) or response (glioma, no drug concentration, assumed uniform
   protocol). No open dataset does both in the same tissue.
5. **Glioma + temozolomide, spatially.** Does not exist. Re-verified first-hand against
   METASPACE in this pass (§1.2).
6. **Anything requiring dbGaP, an MTA, or an institutional Signing Official.**
   Structurally unavailable to an unaffiliated individual. Treat as absent, not as
   pending.

**The one thing worth saying about this list that is not a limitation.** Items 1, 4 and
5 constrain the *response* axis. **They do not constrain the transport axis, and 12 of
the 14 routes above live entirely on the transport axis.** That asymmetry is the honest
shape of this product: the physics can be verified further than it has been, and the
response cannot be verified at all. Both halves of that sentence should appear on every
surface, in that order.

---

## 6. Provenance — every external source, fetched in this pass

All fetched 2026-09-20 unless stated. Where I quote a sentence, I fetched the document
containing it. Where I did not fetch the numbers, I say so in the route.

| source | how | what I relied on |
|---|---|---|
| Vargová 2003, *Glia* 42(1):77-88, PMID 12594739, DOI 10.1002/glia.10204 | NCBI eutils efetch, abstract | control cortex α 0.24, λ 1.55, k′ 3.66e-3 s⁻¹; α and λ further increased in WHO III/IV. **Grade-specific values NOT obtained** |
| Grantab & Tannock 2012, *BMC Cancer* 12:214, PMID 22672469, PMC3407510 | full text, `curl` to pmc.ncbi.nlm.nih.gov, HTTP 200, 195,531 bytes | packing density 50±7.3→40±1.6 % and 40±6.8→31±4.5 %; EGTA penetration ×3.3–4.4 / ×1.5–1.9 (5-FU, gemcitabine), doxorubicin ~×2 (p=0.02, p=0.002); bortezomib 1 µM ×1.7 / ×3 / ×2; observable `C/C∞` at 6 h with [14C]-sucrose internal standard; in vivo doxorubicin distribution ~2× greater in loosely packed xenografts |
| Grantab 2006, *Cancer Res* 66(2):1033-9, PMID 16424039 | eutils abstract; Unpaywall API on DOI 10.1158/0008-5472.CAN-05-3077 | `"oa_status":"closed"`, `"has_repository_copy":false`, `"oa_locations":[]`. aacrjournals.org returns a Cloudflare challenge (HTTP 403). **Not obtainable** |
| Hong, Wilson & Hicks 2019, *Neoplasia* 21(2), PMID 30591421, PMC6314220, DOI 10.1016/j.neo.2018.11.009 | full text, `curl` to pmc.ncbi.nlm.nih.gov, HTTP 200, 224,138 bytes | D(TH-302) = (1.82±0.07)e-7 cm²/s; D(IPM) = (1.33±0.04)e-7; D(mannitol, HCT116 MCL) = (1.46±0.02)e-7; kmet,0 = 0.0115±0.05 s⁻¹; thickness 139±2 µm by ¹⁴C-urea; Table S1 holds all parameters |
| Phillips, Loadman & Cronin 1998, *Br J Cancer* 77(12):2112-9, PMID 9649122, PMC2150429 | eutils efetch, abstract | DLD-1 MCL 50.9±12.1 µm; t½ = 1.36 h (tirapazamine), 2.38 h (NITP), 4.62 h (EO9), 13.1 h (doxorubicin) |
| Hicks et al. 1998, *IJROBP* 42(3):641-9, PMID 9806526 | eutils efetch, abstract | V79 D_M = 7.4e-7 cm²/s, *"12-fold lower than in culture medium"*; MGH-U1 1.3e-6; k_met anoxic 0.52 and 0.31 min⁻¹. **Already the repo's anchor — fetched to confirm what the 12× refers to** |
| Hicks et al. 2006, *JNCI* 98(16):1118-28, PMID 16912264 | eutils efetch, abstract | D range 0.027e-6 to 1.87e-6 cm²/s across 16 TPZ analogs; R² = 0.89 with extravascular transport vs 0.32 without |
| Primeau et al. 2005, *Clin Cancer Res* 11(24):8782-8, PMID 16361566 | eutils efetch, abstract | doxorubicin halves at 40–50 µm from vessels; vessel-to-hypoxia distance 90–140 µm. **Already the `alpha` fit — fetched to confirm the in-sample status** |
| Lally et al. 2006, *Cancer J* 12(6):461-6, PMID 17207315 | eutils efetch, abstract | Eppendorf pO₂ Histograph, 23 analysed patients (13 HGG, 10 LGG), **median tumour pO₂ 5.1 mmHg** |
| Collingridge et al. 1999, *Radiother Oncol* 53(2):127-31, PMID 10665789 | eutils efetch, abstract | *"Both glioma and peritumoural brain contain regions at oxygen tensions less than 2.5 mmHg"* |
| Rampling et al. 1994, *IJROBP* 29(3):427-31, PMID 8005794 | eutils efetch, abstract | 20 patients, 6 tracks × 22 mm, 192 readings, *"significant regions of hypoxia"* |
| Boucher, Baxter & Jain 1990, *Cancer Res* 50(15):4478-84, PMID 2369726 | eutils efetch, abstract | 86 % of profiles show a peripheral gradient → interior plateau; plateau reached **0.2–1.1 mm from the surface**; n = 9 + 6 + 7 |
| Konings et al. 2009, *Cancer Chemother Pharmacol* 64(3):509-16, PMID 19096847, PMC2691802 | eutils efetch, abstract | carboplatin, 8 patients, evaluable in tumour 6/8, **tissue:plasma exposure ratios 0.64–1.46** |
| Konings et al. 2011, *Cancer Chemother Pharmacol* 67(5):1055-62, PMID 20652702, PMC3082021 | eutils efetch, abstract | 5-FU, 7 patients, 5-day continuous infusion, tumour:plasma ratio **increases during the infusion**; circadian rhythmicity in tumour ECF |
| Secomb rat tumour network | `curl`, https://sites.arizona.edu/secomb/files/2024/01/rattum98-ffab1d1a873244ec.txt — **HTTP 200, 8,049 bytes, 106 lines**; page https://sites.arizona.edu/secomb/three-dimensional-network-data-rat-tumor-1998/ | segment/node/diameter/length/xyz columns in microns; region 550 × 520 × 230 µm; R3230AC mammary carcinoma, rat dorsal skin flap; cite Secomb, Hsu, Braun, Ross, Gross & Dewhirst, *Oxygen Transport to Tissue XX*, Plenum 1998, 629-634 |
| METASPACE | live GraphQL query to https://metaspace2020.eu/graphql | `countDatasets` = **20,417**; `simpleQuery`: imatinib **2**, olaparib **0**, temozolomide **0**, doxorubicin **0**, paclitaxel **0** |
| GEO GSE235475 (CIVO) | NCBI eutils esummary, db=gds | 12 HNC patients, subasumstat, exposure regions **~1000–2000 µm diameter**, **214** exposed / **140** unexposed GeoMx ROIs, CosMx subset, released 2023-06-28 |
| BioFVM | eutils efetch, abstract; PMID 26656933, DOI 10.1093/bioinformatics/btv730 | *"first-order accurate in time and second-order accurate in space"*; Apache License v2.0 |
| Salari & Knupp 2000, SAND2000-1444 | `curl` https://www.osti.gov/servlets/purl/759450 — HTTP 200, 10,006,067 bytes; title page and abstract extracted | *"The strength of the MMS procedure is that it can identify any coding mistake that affects the order-of-accuracy of the numerical method."* |
| Peruzzi/Jonas 2023 preprint | **not re-fetched in this pass**; relied on `results/zenodo_kibur_deposit_check.md`, which fetched and read all 34 pages of the CC-BY-4.0 Research Square preprint (DOI 10.21203/rs.3.rs-1647435/v2) | Fig. 2a/c is a line plot of local drug dose (µM) vs distance (µm), one line per patient, doxorubicin and lapatinib |

---

## 7. What this document is not

* **Not a result.** Nothing has been run. Every "engine predicts" is cited to an existing
  `results/` document or to a file and line, not computed here.
* **Not a registration.** No row was added to `tmesim/validation/scorecard.py`. Each
  route that gets run needs its own pre-registration with the bar fixed first, in the
  form `results/transport_anchor_prereg.md` established — and the lesson that document's
  execution taught applies to all of them: **a bar is only credible if you can point at
  something it failed, and the required response to a failure is to constrain the
  coefficient, not to widen the bar.**
* **Not exhaustive on the negative side.** §5 is short by design. §1.2 carries the dead
  ends forward rather than re-deriving them.
* **Not a claim that any of this makes the product work.** Twelve of fourteen routes
  test transport physics. The response axis remains unverifiable, and every route above
  could pass without changing that by one word.
