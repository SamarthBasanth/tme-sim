# Scoring the engine's doxorubicin cell-membrane permeability against Degerstedt 2024 — a PRE-REGISTERED test

**Status: SCORED. Pre-registered before computation. Nothing in `tmesim/` was changed; `alpha` was not touched.**

Commissioned to *"score the transvascular permeability term against a measurement."*
**That framing is corrected in §1 before anything is scored: the engine has no
transvascular wall permeability, and the quantity the scouting pass actually found
a measurement for is a CELL-membrane permeability.** The test below is the
cell-membrane one, which is the comparison that exists.

Calibration guard, start and end of pass:

```
$ .venv/bin/python -c "from tmesim.validation.calibration import calibrate_alpha; a=calibrate_alpha(); assert a==22.223477711172563, a; print('alpha OK',a)"
alpha OK 22.223477711172563
```

Artefacts (all new; no existing file in `results/` or `archive/` was edited):

| file | what it is |
|---|---|
| `experiments/degerstedt_permeability_score.py` | the two-stage script; `--stage prereg` refuses to overwrite a frozen pre-registration, `--stage score` refuses to run without one |
| `experiments/degerstedt_prereg_2026-10-10.json` | the **frozen pre-registration**, `sha256 = bd1c6447f31c55ed3edeb938b6e03e8426cd0652dc67cc96a08a9b00332d2c7b` |
| `experiments/degerstedt_score_2026-10-10.json` | every number below, machine-readable, carrying that same sha256 |

---

## 0. Headline

1. **The measurement is commensurable.** The paper's `P_cell` is defined by
   `P = k_cell·V_cell/(A_cell·60)` with `k_cell` in min⁻¹ and `V/A = r/3` — i.e.
   the exact inversion of a membrane-flux law — and the engine's `P` enters as an
   interfacial conductance `σ_c/(a·P)`. Both are µm/s, both are flux per unit
   concentration difference across a cell membrane, both use a cell `V/A` as the
   length scale. **No unit or definition mismatch manufactures this discrepancy.**
   Independent check: Eq. 2 reproduces the published `9.00 × 10⁻⁴ µm/s` from the
   published `k_cell = 0.0163 min⁻¹` at an implied cell radius of **9.94 µm**, a
   physically sensible HepG2 radius.

2. **The scouting pass's 935× reproduces: 934.17×.** But it is **one of two**
   shipped doxorubicin "membrane permeability" numbers, and **they straddle the
   measurement.** `two_route.membrane_permeability_um_s` is 934× BELOW it;
   `permeability.membrane_permeability` is 18,971× ABOVE it. The two engine terms
   disagree with **each other by 1.77 × 10⁷ (7.25 decades)**. No one-sided
   statement of the form "the engine's permeability is X-fold off" is available.
   **The internal inconsistency is the finding, not the gap.**

3. **The 934× gap has no consequence for any shipped doxorubicin transport
   number.** Substituting the measured permeability into the shipped two-route law
   moves `D_eff` from **27.92558 → 27.93203 µm²/s, +0.023 %**. The scouting pass's
   description of it as *"the most consequential number in the file"* is not
   supportable on `D_eff` grounds.

4. **The scouting pass's α reframing is NOT SUPPORTED.** The claim was that α
   patches *"a measurably wrong permeability"* rather than a binding capacity.
   Converted into the fastest first-order sink a membrane-limited process can
   sustain at the engine's own geometry, the measured permeability supplies
   **1.287 × 10⁻³ s⁻¹** against α's calibrated **5.035 × 10⁻³ s⁻¹** — between
   **5.4 % and 25.6 %** of it, depending on concentration basis. Closing the
   permeability gap *completely* still leaves α short by **3.9×**. The 934×
   permeability gap and the 17.6× α amplification are neither the same magnitude
   nor the same code path.

5. **What the measurement establishes instead is sharper than what was claimed,
   and it also runs against the engine.** Because the measured permeability is a
   *ceiling* on uptake rate and not merely a smaller number, and because **a
   binding capacity sets how long a flux can be absorbed, not how fast it
   arrives**, α = 22.22 asserts a perivascular cellular-uptake sink
   **3.9× to 18.6× above the maximum rate a measured doxorubicin cell
   permeability permits.** `calibration.py`'s "high-capacity intracellular / DNA
   binding" attribution cannot rescue that, because capacity is the wrong
   operator for a rate bound. This is a **new, independent, measured upper bound
   that α violates** — not the reframing the scouting pass asserted, and not a
   vindication of the existing docstring either.

6. **The engine is internally inconsistent about doxorubicin by 162× before α and
   3608× after it.** The engine's own `P_mem = 9.63 × 10⁻⁷ µm/s` can sustain a
   cellular-uptake sink of at most **1.378 × 10⁻⁶ s⁻¹**, while the same compile's
   own `binding` is **2.237 × 10⁻⁴ s⁻¹** (162×) and **4.972 × 10⁻³ s⁻¹** after α
   (3608×). Two shipped terms describing the same physical act — doxorubicin
   entering a cell — disagree by three orders of magnitude *inside one compile*,
   and neither is tested against the other.

---

## 1. The commission's framing is wrong, and the correction has to come first

The task was to score **"the transvascular permeability term."** There is no such
term to score.

`tmesim/treatment/transvascular.py` states it in its own docstring, verbatim:

> A genuine transvascular permeability is `PS [1/s] = P [um/s] * (S/V) [1/um]`:
> a wall permeability times a vessel surface-to-volume density. **The repository
> has neither factor.** It has `bbb_penetrance`, a dimensionless brain:plasma
> PARTITION coefficient that three separate modules state is not a permeability,
> multiplied by an undocumented 1 s⁻¹. And `S/V` needs a vessel radius that this
> repo does not carry as a measured quantity.

Confirmed live: the shipped entry term is
`PS = vascular_density × bbb_penetrance × PS_SCALE_S_INV` with
`PS_SCALE_S_INV = 1.0`. `implied_permeability_um_s` can back-compute a wall
permeability, but only by the caller supplying a vessel radius the repo does not
have; at `R = 3/5/10 µm` it returns `0.0309 / 0.0514 / 0.1028 µm/s` for
doxorubicin, and its own docstring says **"Do not quote the output of this
function as the engine's permeability."** It is honoured here: those three numbers
are not scored.

And **Degerstedt's number is the wrong anchor for a wall permeability anyway.** It
is a permeability across a *tumour cell's* membrane, measured from intracellular
accumulation inside a gel. It is the right anchor for the engine's **cell**
membrane permeability, which is what the scouting pass actually compared it to,
and what is scored below.

*(The paper does quote, as its ref [65], a doxorubicin permeability of
`32 ± 23 × 10⁻⁴ µm/s` across an MDCK monolayer lining a cylindrical microvessel —
the closest thing in it to a wall permeability. **That value is NOT verified at
its own source** and is used below only to widen one sensitivity band. It is
marked UNVERIFIED everywhere it appears.)*

---

## 2. The engine's actual permeability paths and their real values

There are **two** functions in the shipped package that both return "a passive
membrane permeability for doxorubicin in µm/s", and they are both live. Values
re-derived here from a live compile, not inherited:

| | path | doxorubicin value | what it feeds |
|---|---|---|---|
| **A** | `tmesim/molecular/two_route.py::membrane_permeability_um_s` | **9.634241 × 10⁻⁷ µm/s** | `two_route_factor`'s interfacial resistance `σ_c/(a·P)` → the tissue-scale diffusivity hindrance → `DrugParams.diffusivity` |
| **B** | `tmesim/molecular/permeability.py::membrane_permeability` | **17.074278 µm/s** | `MolecularProfile.membrane_permeability`; **reported, not solved with** — it does not enter the PDE. Read by `validation/moncayo_explant.py` and `validation/glioma_regional_pk.py` |

**Ratio B/A = 1.772 × 10⁷.**

### 2.1 Which one the comparison is about, and why

**A.** It is the one that reaches the transport solver. `compile_drug.build_profile`
runs with `two_route=True` by default and sets
`effective_D = interstitial_D × two_route_out["hindrance"]`; that becomes
`DrugParams.diffusivity`. **B** is carried on the profile and reported, but
`molecule_to_drugparams` maps only `profile.effective_diffusivity`,
`tissue_binding_rate`, `clearance` and `bbb_penetrance` into `DrugParams` — **B
never enters the PDE.** So A is the shipped permeability in the sense that
matters, and A is what the scouting pass quoted.

**B is still scored (TEST 2), because its existence is what destroys the
one-sided claim.**

### 2.2 Why A is so small — the mechanism, not a mystery

```
log10 K_mem = min(logP − 1.34·n_HBD,  3.0)
```

For doxorubicin as the engine compiles it: Crippen `logP = 0.001`, `n_HBD = 6`.

```
log K_mem = 0.001 − 1.34×6 = −8.039      (NOT clamped)
```

Six hydrogen-bond donors at 1.34 log units each is an **8-decade** partition
penalty. `A = K_mem × D_mem/d` then lands at 9.63 × 10⁻⁷ µm/s. The number is a
direct, auditable consequence of `HBD_LOGP_DECREMENT = 1.34` applied to a
6-donor molecule, and `two_route.py` already documents that constant as an
**UPPER estimate** for a bilayer (hexadecane is the harshest membrane proxy; a
real bilayer's interface supplies H-bond acceptors). The direction of that
admitted bias is exactly the direction of the gap measured below.

The engine's two values are also bracketed by a band **the repository itself
already quotes**. `permeability.py`'s own docstring: *"real Caco-2 Papp for these
drugs is 0.001–1 µm/s, i.e. 1–4 orders below what this returns."* Against that
self-cited band:

| | value | vs the repo's own stated real band |
|---|---|---|
| **A** `two_route` | 9.634 × 10⁻⁷ µm/s = 9.6 × 10⁻¹¹ cm/s | **1038× BELOW** the bottom |
| **measured** | 9.00 × 10⁻⁴ µm/s = 9.0 × 10⁻⁸ cm/s | **0.90× the bottom** — essentially at it |
| **B** `permeability.py` | 17.074 µm/s = 1.7 × 10⁻³ cm/s | **17× ABOVE** the top |

The measurement lands where the repository's own prose says real drug
permeabilities land. **Both shipped numbers do not.** That is a third,
in-repo-only corroboration that the measurement is the right order and the engine
is not.

---

## 3. The measurement, as I verified it

**Source.** Degerstedt O, O'Callaghan P, Clavero AL, Gråsjö J, Eriksson O,
Sjögren E, Hansson P, Heindryckx F, Kreuger J, Lennernäs H. "Quantitative
imaging of doxorubicin diffusion and cellular uptake in biomimetic gels with
human liver tumor cells." *Drug Deliv Transl Res* 2024;**14(4)**:970–983.
**PMID 37824040, PMCID PMC10927899**, doi:10.1007/s13346-023-01445-1.

**Licence, verbatim from the full text:** *"This article is licensed under a
Creative Commons Attribution 4.0 International License, which permits use,
sharing, adaptation, distribution and reproduction in any medium or format."*
(`https://creativecommons.org/licenses/by/4.0/`). Retrieved as JATS XML from
Europe PMC and read directly, not from an abstract.

### 3.1 What was actually measured — this is the part that decides everything

The scouting pass's number is confirmed, and the chain behind it is now explicit.

**The raw observable** is intracellular doxorubicin fluorescence, converted to
concentration against a calibration curve, in 250 µm zones along the gel. Methods,
verbatim:

> The average cellular DOX fluorescence intensities for each zone were converted
> to DOX concentrations using calibration curves … By linear regression, the
> intracellular uptake of DOX (µM / min) was determined as the slope of the
> concentration vs time relation. The obtained intracellular uptake for each zone
> was in turn then plotted over the gel concentration for the corresponding zone
> (µM) at three timepoints (60, 120 and 180 min) resulting in a linear
> relationship between intracellular uptake and concentration. **By linear
> regression, the overall transport rate constant of DOX into the cell
> (k_cell, min⁻¹) was determined as the slope of the intracellular uptake vs
> concentration relation.**

**The definition**, Eq. 2 verbatim:

> `P_cell = (k_cell × V_cell) / (A_cell × 60)`
> where V_cell and A_cell are the volume and area of the cells respectively.
> **Approximating the cells are spherical, the ratio V_cell / A_cell can be
> simplified to r/3 where r is the cell radius.**

**The values**, verbatim:

> … resulted in an overall transport rate constant for DOX (k_cell) of
> **0.0163 ± 0.001 min⁻¹** and consequently an apparent cell permeability
> (P_cell) of **9.00 ± 0.74 × 10⁻⁴ µm/s**. The lower density of HepG2 cells
> resulted in slightly higher maximum cellular uptake rate … even if the k_cell
> and P_cell were lower, with **0.0120 ± 0.002 min⁻¹** and
> **6.49 ± 1.25 × 10⁻⁴ µm/s**, respectively.

**System:** HepG2 human hepatoma cells at 2 × 10⁶ (high density) or 1 × 10⁶ (low
density) cells/mL, suspended in a "biomimetic cirrhotic" gel of rat-tail type I
collagen 5 mg/mL + type I-S fibrinogen at pH 6.86, in a 3D-printed-mould chip
with a 20 µM doxorubicin donor. **Not GBM, not brain, not in vivo.**

**It is an APPARENT permeability** — the paper says so in the name, and it is
derived from one-way intracellular accumulation with no back-flux term. It
therefore lumps active transport, binding and partition. Doxorubicin is a P-gp
substrate and a DNA intercalator. The engine's `P_mem` is a *passive* permeability.
**These are not strictly the same object.**

**Direction of its own error, stated by the authors:** fluorescence self-quenches
above ~25 µM, a level reached within 90 min in the zone nearest the donor, and
DNA binding quenches further. Verbatim: *"This all suggests that we are more
likely to underestimate rather than overestimate the determined DOX apparent cell
permeability."* **So 9.00 × 10⁻⁴ µm/s is a LOWER bound**, and every ceiling
derived from it below is a lower bound on a ceiling.

### 3.2 Is it commensurable at all? — TEST 0, decided by reading, before computing

This was the gate: a unit or definition mismatch here would manufacture a 10³
discrepancy out of nothing.

| requirement | measurement | engine (`two_route`) | verdict |
|---|---|---|---|
| units µm/s | stated µm/s | stated µm/s | ✔ |
| defined as flux per unit Δconcentration across a **cell** membrane | Eq. 2 is the exact inversion of `V·dC_in/dt = P·A·ΔC` | `biot = σ_c/(a·P)`; `1/P` is an interfacial resistance in series with the inclusion interior (Hasselman & Johnson 1987 Eqn 8) | ✔ |
| length scale is a cell `V/A` | `V/A = r/3`, `r` the cell radius | `a = DEFAULT_CELL_RADIUS_UM = 7.8517 µm`, entering only through `σ_c/(a·P)` | ✔ |
| the `/60` is a unit conversion, not a model factor | `k_cell` is min⁻¹; `/60` → s⁻¹ | n/a | ✔ |

**Independent arithmetic check that no factor was lost.** Inverting Eq. 2 on the
published pair:

```
r = P_cell × 60 × 3 / k_cell = 9.00e-4 × 60 × 3 / 0.0163 = 9.9387 µm
```

A HepG2 radius of **9.94 µm** (≈ 20 µm diameter) is physically correct. Eq. 2
therefore closes on its own published inputs at a sensible radius. **There is no
hidden factor of 60, 10³ or anything else.**

> **TEST 0: PASS. The comparison IS commensurable.** This answer was fixed before
> any ratio was computed, and it is the one that most of the finding depends on.

---

## 4. The pre-registration

Frozen to `experiments/degerstedt_prereg_2026-10-10.json`,
`sha256 = bd1c6447f31c55ed3edeb938b6e03e8426cd0652dc67cc96a08a9b00332d2c7b`,
before any of §5 was computed. `--stage score` will not run without it and
records that hash in its output. The script refuses to overwrite it.

**Disclosure, stated in the pre-registration itself: THIS IS NOT A BLIND
PRE-REGISTRATION.** The engine values and the measured value were both in hand
when the criteria were frozen, and a hand conversion of `P_cell` into a sink rate
had already been sketched. What the freeze buys is that **the decision rules and
predicted bands cannot be moved after the numbers land** — not blindness. A reader
who wants a blind test should treat TEST 1 and TEST 2 as bookkeeping and TEST 3
and TEST 4 as the registered predictions, since neither of those had been computed
with the real shipped inputs at freeze time.

| test | registered prediction | registered pass criterion |
|---|---|---|
| **0** commensurability | — (decided by reading) | all four rows of §3.2 hold, else **everything below is VOID and the answer is "not commensurable"** |
| **1** the ratio | engine **LOWER**, `R1 > 1` | `400 ≤ R1 ≤ 2000` |
| **2** which engine term | engine **HIGHER**, `R2 > 10⁴` | registered for the record (value in hand), **plus** a registered consequence: if `R1 > 1` **and** `R2 > 1`, the two shipped terms straddle the measurement and the headline **must** be the internal inconsistency, not a one-sided gap |
| **3** consequence for the shipped number | `D_eff` moves **UP, negligibly**; reason given in advance: `σ_c/(a·P) = 2.18 × 10⁷` at shipped `P`, so a 934× rise leaves it near 2.3 × 10⁴ ≫ 1 | `|ΔD_eff|/D_eff < 0.02`. A pass **forbids** calling the gap "the most consequential number" on `D_eff` grounds and forces the α argument to stand alone |
| **4** the α test | outcome **B**, `k_ceiling` in `[8 × 10⁻⁴, 1.5 × 10⁻³] s⁻¹`, **3–7× BELOW** `k_eff(α)` | **A** `k_ceiling ≥ k_eff(α)` → reframing **SUPPORTED**; **B** `k_struct ≤ k_ceiling < k_eff(α)` → **partial overlap**, assert neither framing; **C** `k_ceiling < k_struct` → reframing **REFUTED** |
| **5** transvascular term | no such term exists | definitional correction, established by reading code — **not** a test |

TEST 4's construction, fixed in advance: membrane-limited uptake cannot remove
extracellular solute faster than `P × (S/V)`. At the engine's own geometry
`S/V_tissue = 3φ/a`. Both concentration bases are reported;
`k_ceiling_interstitial = k_ceiling_tissue/(1−φ)` is the **permissive** end and
**the verdict uses it**.

Also registered in advance, as the reading of outcome B: *"α asserts a
perivascular sink ABOVE the ceiling that a measured doxorubicin cell permeability
can sustain. A binding CAPACITY cannot rescue that: capacity sets how long a flux
can be absorbed, permeability sets the flux. So B is NOT 'the same hole counted
twice'."*

---

## 5. Results

Engine geometry used throughout, all shipped values:
`φ_cell = 0.7892`, `a = 7.8517 µm`, `S/V = 3φ/a = 0.301540 µm⁻¹`.

### TEST 1 — the ratio. **PASS, as predicted.**

| | doxorubicin cell permeability |
|---|---|
| engine, shipped (`two_route`) | **9.634241 × 10⁻⁷ µm/s** |
| Degerstedt, measured (HD HepG2) | **9.00 × 10⁻⁴ µm/s** |
| **R1 = measured / engine** | **934.17×** — engine LOWER |

Inside the registered `[400, 2000]` band. The scouting pass's 935× reproduces.

### TEST 2 — which engine term. **PASS, and it changes the headline.**

| | value | vs measured |
|---|---|---|
| `permeability.membrane_permeability` | **17.074278 µm/s** | **18,971× ABOVE** |
| `two_route.membrane_permeability_um_s` | **9.634241 × 10⁻⁷ µm/s** | **934× BELOW** |
| **internal spread** | **1.772 × 10⁷ (7.25 decades)** | — |

`straddle = true`. **The registered consequence fires.** The engine does not have
"a permeability that is 935× wrong." It has two permeabilities for the same
molecule that differ by 7.25 decades, with the measurement between them. Which
number you quote decides the **sign** of the disagreement. Any sentence of the
form "the engine's doxorubicin permeability is X-fold off" is under-specified.

*(`permeability.py`'s own docstring already confesses the direction of its own
error — "real Caco-2 Papp for these drugs is 0.001–1 µm/s, i.e. 1–4 orders below
what this returns." It was never scored against a number. It is now: 18,971×.)*

### TEST 3 — consequence for the shipped transport number. **PASS: none.**

Measured `P` substituted into the shipped law as a probe; every other input at its
shipped value; nothing in `tmesim/` changed.

| | `D_eff` (µm²/s) | `σ_c/(a·P)` | transcellular gain |
|---|---|---|---|
| shipped, `P = 9.634 × 10⁻⁷` | **27.925583** | 2.1808 × 10⁷ | 1.0000002 |
| probe, `P = 9.00 × 10⁻⁴` (measured) | **27.932027** | 2.3345 × 10⁴ | 1.0002310 |
| ceiling, `P → ∞` | 169.017126 | 0 | 6.0524 |

**Relative change: +0.023 %.** Well inside the registered `< 2 %`.

The membrane resistance ratio falls by 934× exactly as predicted and is *still*
2.3 × 10⁴ ≫ 1, so the cells stay effectively impermeable. The transcellular route
*does* have headroom — the `P → ∞` ceiling is 6.05× — but the measured
permeability reaches **0.004 %** of the way into it. Supplementary: `D_eff` would
need `P ≈ 0.397 µm/s` to rise even 10 %, which is **442× the measured value** and
**4.1 × 10⁵×** the engine's.

> **The registered consequence fires: the 934× gap has NO effect on any shipped
> doxorubicin transport number**, and the scouting pass's characterisation of it as
> *"the most consequential number in the file"* is not supportable on `D_eff`
> grounds. The α argument now has to stand entirely on its own.

### TEST 4 — the α test. **Outcome B, exactly as predicted and inside the predicted band.**

Maximum first-order removal rate a membrane-limited process can sustain at the
engine's own `φ` and `a`:

| `P` source | `P` (µm/s) | `k_ceiling` tissue basis (s⁻¹) | `k_ceiling` interstitial basis (s⁻¹) |
|---|---|---|---|
| **engine's own** `P_mem` | 9.634 × 10⁻⁷ | 2.905 × 10⁻⁷ | 1.378 × 10⁻⁶ |
| **measured, LD HepG2** | 6.49 × 10⁻⁴ | 1.957 × 10⁻⁴ | 9.284 × 10⁻⁴ |
| **measured, HD HepG2** | **9.00 × 10⁻⁴** | **2.714 × 10⁻⁴** | **1.287 × 10⁻³** |
| MDCK microvessel, **UNVERIFIED at source** | 32 × 10⁻⁴ | 9.649 × 10⁻⁴ | 4.577 × 10⁻³ |

Against the engine's loss rates for the same compile:

| | value (s⁻¹) |
|---|---|
| `k_struct` = `clearance + binding`, uncalibrated | 2.868 × 10⁻⁴ |
| `k_eff(α)` = `clearance + α·uptake_proxy`, α = 22.2235 | **5.035 × 10⁻³** |
| α's amplification on `k` | **17.56×** |

`k_struct = 2.868 × 10⁻⁴ ≤ k_ceiling = 1.287 × 10⁻³ < k_eff(α) = 5.035 × 10⁻³`

> **Outcome B — PARTIAL OVERLAP.** Predicted B. `k_ceiling = 1.287 × 10⁻³ s⁻¹`
> fell inside the registered `[8 × 10⁻⁴, 1.5 × 10⁻³]` band, and the shortfall
> **3.91×** fell inside the registered 3–7×.

**The decisive ratios:**

| statement | factor |
|---|---|
| `k_eff(α)` / ceiling, **interstitial** (permissive) basis | **3.91×** |
| `k_eff(α)` / ceiling, tissue-average basis | **18.55×** |
| `k_eff(α)` / ceiling from the **UNVERIFIED** MDCK value | **1.10×** |
| engine's own `binding` / ceiling from engine's own `P_mem` | **162.4×** |
| engine's **calibrated** `binding` / ceiling from engine's own `P_mem` | **3608×** |

**Using the engine's own `a = 7.8517 µm` is the conservative choice.** At the
paper's own implied `r = 9.9387 µm` the ceiling drops to `1.017 × 10⁻³ s⁻¹` and
α's excess *rises* to 4.95×. The verdict's direction is not an artefact of whose
cell radius is used.

---

## 6. What this does and does not imply for α

### 6.1 The reframing the scouting pass asserted is NOT SUPPORTED

The scouting pass wrote:

> An **external, measured** doxorubicin permeability 935× above the engine's is
> independent evidence that the missing ~300× lives in the **membrane/uptake
> leg** — a named, shipped parameter — rather than in an unmeasurable "binding
> capacity". That reframes α from a mechanistic constant into a patch for a
> measurably wrong permeability.

**Three separate reasons that does not follow.**

**(a) The magnitudes do not match, in either direction.** The permeability gap is
934×. α's amplification of `k` is 17.56×. If the permeability gap *were* the α
hole, closing it would close α. It does not: substituting the measured
permeability raises the sustainable sink ceiling to `1.287 × 10⁻³ s⁻¹`, which is
still **3.91× short** of α's `5.035 × 10⁻³ s⁻¹`. **A 934× repair cannot fill a
17.6× hole and still leave it 3.9× open unless the two were never the same
quantity** — and they were not: the 934× is a *permeability* error, the 3.91× is
what remains of a *rate* deficit after that error is fully repaired.

**(b) They are structurally independent in the code, not one hole seen twice.**
α multiplies `uptake_proxy`, which is `drug.binding`, which is
`molecular.binding.tissue_binding_rate(logP, HBD)`. That function **never reads
`P_mem`**, and `P_mem` never reaches `binding`. `P_mem` reaches the solver only
through the diffusivity hindrance (TEST 3), and `binding` reaches it only through
`k = clearance + binding` at `physics/diffusion.py:315`. **Repairing `P_mem` would
not change `binding` by one digit** — TEST 3 measured exactly that: 0.023 %.
There is no code path on which the two could be double-counted.

**(c) The direction of the remaining error is wrong for the claim.** The claim
needs the measurement to *supply* α's missing sink. It supplies between 5.4 % and
25.6 % of it and **caps** the rest.

### 6.2 What the measurement establishes instead — and it is worse for α, not better

The measured permeability is not merely a smaller number than α needs. **It is a
ceiling**, because a first-order sink `k·C` acting on extracellular solute
requires that solute to physically cross a membrane at rate `k·C`, and the
fastest a membrane of permeability `P` can pass it is `P·(S/V)·C`. The ceiling
holds at steady state and holds for any intracellular capacity whatsoever.

**And this is the point `calibration.py` cannot answer.** Its docstring attributes
the missing sink to *"rapid cellular uptake plus high-capacity intracellular /
DNA binding"* and names α as encoding *"the intracellular binding CAPACITY
(bound:free ratio) the structure-only heuristic cannot know."*

**A capacity cannot raise a flux ceiling.** Capacity determines how long a given
influx can be absorbed before the cell saturates and back-flux begins; it does
nothing to the influx rate itself. So the half of `calibration.py`'s attribution
that could in principle reach `5.035 × 10⁻³ s⁻¹` is *"rapid cellular uptake"* —
and that is precisely the quantity Degerstedt measured, at
`9.00 × 10⁻⁴ µm/s`, capping the sink at `1.287 × 10⁻³ s⁻¹`.

> **α = 22.2235 asserts a perivascular cellular-uptake sink 3.91× (permissive
> basis) to 18.55× (tissue basis) above the maximum rate that a measured
> doxorubicin cell permeability can sustain at the engine's own cell volume
> fraction and radius. The "high-capacity binding" attribution cannot rescue
> this, because capacity is the wrong operator for a rate bound.**

This is a genuinely new constraint on α, and it is independent of every anchor α
has been scored against before. It is **not** the reframing that was claimed, and
it is **not** a vindication of the existing docstring either.

**Two further internal facts make it harder to dismiss.**

First, `results/uptake_proxy_double_count.md` already established — measured on
the shipped solver, not argued — that `uptake_proxy` is **wired into an
irreversible sink, not a capacity**: multiplying `binding` by 10 moves the
steady-state concentration to 0.1089×, while a capacity of 10 leaves it at
0.999999999998×. So the term α scales really is a removal rate in s⁻¹, and the
flux ceiling really does bind it. That document's own conclusion — *"the engine
holds a per-compound quantity that its own docstring calls 'a capacity ordering'
and applies it as an irreversible sink … a mis-assignment"* — is corroborated
here from an entirely external direction, and now has a **number** attached to the
mis-assignment: 3.91×.

Second, the engine's own two terms were **already** inconsistent before any
external measurement arrived. Its own `P_mem` supports at most
`1.378 × 10⁻⁶ s⁻¹`; its own uncalibrated `binding` is `2.237 × 10⁻⁴ s⁻¹`
(**162×** more) and its calibrated `binding` is `4.972 × 10⁻³ s⁻¹` (**3608×**
more). **Degerstedt resolves which of the two is closer to right on the
permeability side** — `P_mem` is 934× too small, `binding` is within 1.2× of
what the measured permeability supports before α — **and then leaves α's 17.6×
standing on nothing.**

### 6.3 The honest reading, and the escape route that has to be named

**Where it lands:** the measured permeability almost exactly reproduces the
engine's *uncalibrated* structural sink (`2.714 × 10⁻⁴` tissue-basis against
`k_struct = 2.868 × 10⁻⁴ s⁻¹`, a ratio of 0.946). The structure-only prior that
`pk_reality.md` called ~300× too slow turns out to be **the right order for
membrane-limited doxorubicin uptake.** What is unaccounted for is the amplification
α adds on top of it.

**The escape route, stated so it can be attacked.** The ceiling bounds the
**cellular-uptake** component of the sink. If part of α·`uptake_proxy` represents
**extracellular** sequestration — doxorubicin binding to collagen and ECM protein,
which `tissue_binding_rate`'s logP/HBD form does not distinguish from
intracellular binding — then that part is not bounded by a cell permeability and
§6.2's factor shrinks. Two things narrow the route: `calibration.py` names the
mechanism as *intracellular* (*"cell lipids, protein and DNA"*, *"rapid cellular
uptake"*), and reversible extracellular binding is a **capacity**, which
`uptake_proxy_double_count.md` §1.1 proves cannot move a steady state at all,
while an *irreversible* extracellular sink would be degradation — which
`calibration.py` explicitly excludes (*"The missing physics is not slow metabolic
clearance"*). The route is narrow but it is real, and nothing here closes it.

**What it would take to close it:** a measured doxorubicin ECM/collagen binding
rate in s⁻¹, or an independent perivascular sink rate measured directly rather
than inferred from a penetration length. Neither is in this repository.

---

## 7. Limitations — what this cannot establish whatever the numbers said

* **HepG2 hepatoma in a collagen/fibrinogen gel is not GBM in brain.** Nothing
  here transfers a doxorubicin number to another compound or another tissue.
  Cross-compound ranking remains **WITHDRAWN (2026-08-21)** and nothing above
  reinstates it.
* **`P_cell` is APPARENT.** It lumps active transport, binding and partition;
  doxorubicin is a P-gp substrate and a DNA intercalator; the engine's `P_mem` is
  passive. A 934× gap is far wider than that caveat covers, but the caveat is real
  and it was registered before scoring. **For TEST 4 the lumping runs the
  *favourable* way** — an apparent permeability that already contains active
  uptake is the *right* quantity to bound a total uptake sink with, and it is
  still 3.91× short.
* **The measurement is a LOWER bound** by the authors' own statement about
  self-quenching, so every ceiling above is a lower bound on a ceiling. The
  UNVERIFIED MDCK value (32 × 10⁻⁴ µm/s) puts α only **1.10×** above the ceiling
  rather than 3.91×. **The direction of §6.2 survives the entire available
  measured range; the margin does not.** At the top of that range α sits
  essentially *at* the ceiling, which is a weaker claim than 3.91× and must be
  quoted alongside it.
* **No absolute simulator output is validated or invalidated here.** The subject
  is one shipped parameter and one derived rate. No absolute kill fraction was
  computed; none is defensible.
* **TEST 3 is doxorubicin-specific.** `P_mem` is a strong function of HBD count;
  for a low-HBD lipophilic compound the transcellular route does open and the
  membrane leg does matter. Nothing above says the two-route membrane leg is
  inconsequential in general — only that it is inconsequential for doxorubicin,
  which is the one compound α is fitted to.
* **The pre-registration is not blind** (§4). TEST 1 and TEST 2 were bookkeeping.
* **`calibration.py`'s docstring is also stale on a second point** noticed in
  passing and not pursued: it describes α as *"a ~150× amplification"* and
  `_verdict_text` still formats that phrase, while the live α is 22.22 and its
  amplification on `k` is 17.56×. That is the `results/`-is-the-claim-surface
  defect appearing inside `tmesim/` prose. **Not fixed here** — this pass does not
  edit `tmesim/`.

---

## 8. Closing calibration guard

```
$ .venv/bin/python -c "from tmesim.validation.calibration import calibrate_alpha; a=calibrate_alpha(); assert a==22.223477711172563, a; print('alpha OK',a)"
alpha OK 22.223477711172563
```

No parameter was tuned. α was not changed. No file under `tmesim/`, `app/`,
`docs/` or `landing-3d/` was touched. No existing file in `results/` or
`archive/` was edited.
