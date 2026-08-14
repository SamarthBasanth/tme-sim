# tmesim.molecular — The Molecular-Scale Layer

This document describes the molecular layer that sits *upstream* of the tissue
engine documented in [`architecture.md`](architecture.md). Its job is narrow and
specific: turn a real drug **molecule** (a structure, given as SMILES) into the
handful of physical parameters the tissue-scale transport engine needs — a
`core.DrugParams` — using physics-based estimators and a genuine particle
simulation rather than hand-set numbers.

It is a companion to [`../MOLECULAR_CONTRACT.md`](../MOLECULAR_CONTRACT.md), the
frozen API contract for this layer; where the two ever disagree, the contract
wins. The tissue engine and `tmesim.core` are **frozen** — the molecular layer
only *produces* a `DrugParams`; it does not touch the engine.

Conventions (matching the contract): length in **micrometers (µm)**, or **nm**
where noted; time in **seconds**; diffusivity in **µm²/s**; concentration in
**µM**; temperature in **K**. Body temperature = 310 K. Interstitial-fluid
viscosity ≈ 1.2×10⁻³ Pa·s (water at 37 °C).

> **Implementation status.** This block previously said (dated 2026-07) that only
> `types.py` and `stokes_einstein.py` were implemented and that `descriptors`,
> `library`, `brownian`, `permeability`, `binding` and `compile_drug` were merely
> "specified by the contract". **That is long out of date — all of them are
> implemented and on the shipped path.** The layer is now ~22 modules, and §2
> below documents only a subset of them: also present are `efflux.py`,
> `ion_trapping.py`, `two_route.py`, `immune_sink.py`, `uptake.py`,
> `lipophilicity.py`, `mgmt.py`, `byo.py`, `ecm_mesh.py`,
> `binding_site_barrier.py`, `breast_library.py`, `potency_breast.py` and
> `potency.py`. Several are opt-in mechanistic corrections that ship default-OFF.
> Treat §2 as a partial reference and the source as authoritative.

---

## 1. The multiscale idea

### 1.1 Why you cannot just run molecular dynamics on a tissue

The honest temptation is to say "we simulate the molecules, so simulate *all*
the molecules." You cannot. Atomistic molecular dynamics (MD) integrates Newton's
equations for every atom, and both the number of atoms and the number of
timesteps in a microenvironment are astronomically out of reach.

**Atom count.** Liquid water holds ≈ 3.3×10²² molecules/cm³, i.e. ≈ 1×10¹¹ atoms
per µm³. A *single* ~(20 µm)³ voxel — the size of one cell — therefore already
contains roughly

```
(20 µm)³ × 1×10¹¹ atoms/µm³  ≈  8×10³  µm³ × 1×10¹¹  ≈  10¹⁵ atoms.
```

The tissue engine's domain is not one cell — it is ~1.8 mm across (a 120×120 grid
at dx = 15 µm). A full microenvironment slab is thus **~10²⁰–10²¹ atoms**. Even
the single-voxel figure, ~10¹⁵ atoms, is already six orders of magnitude beyond
what routine all-atom MD handles (millions of atoms).

**Timestep count.** MD is bound to a ~1–2 fs (10⁻¹⁵ s) timestep by the fastest
bond vibrations. Drug transport and pharmacodynamics play out over **seconds to
hours** (10³–10⁴ s). Even reaching a mere 1 ms of simulated time takes

```
10⁻³ s / 10⁻¹⁵ s  =  10¹²  timesteps,
```

and the biologically relevant window is 10¹⁸–10¹⁹ timesteps.

**The product is the point.** ~10¹⁵ atoms × ~10¹² timesteps ≈ **10²⁷ atom-steps**
just to reach 1 ms in one voxel. At an optimistic 10⁹ atom-steps/s, that is ~10¹⁸
seconds — on the order of **10¹¹ years**, longer than the age of the universe, for
a single voxel and a millisecond. Whole-tissue, whole-PK atomistic MD is not
"expensive," it is categorically not computable. Anyone who claims to "simulate
the molecular dynamics of a tumor" is either redefining the words or overselling.

### 1.2 The two honest methods used instead

The molecular layer does **not** attempt whole-tissue MD. It uses exactly two
methods, and keeps them clearly separated:

1. **Physicochemical / cheminformatics estimation from structure (RDKit).**
   Cheap, structure-derived descriptors (MW, logP, TPSA, H-bond donors/acceptors,
   etc.) feed well-established empirical relations for permeability, blood–brain
   barrier (BBB) penetrance, protein binding, and clearance. These are
   *estimates* — literature-grade rules of thumb, not first-principles physics.

2. **Genuine overdamped Langevin / Brownian dynamics at the micro-scale.**
   To get an *effective* diffusivity through crowded ECM/cell tissue, the layer
   runs an actual particle simulation: many tracer particles random-walk through a
   periodic box packed with excluded-volume obstacles, and the effective
   diffusivity is read off from their mean-squared displacement. This is real,
   dimensionally-correct transport physics at a scale where it *is* computable
   (hundreds of point tracers, thousands of steps, a few seconds of wall time).

Target potency (IC50) is the deliberate exception — it is **not** reliably
derivable from structure without docking or an assay, and is treated as an input
or a clearly-flagged low-confidence estimate (see §4–§5).

### 1.3 Data flow

```
   drug structure            micro-scale physics              tissue engine
   (SMILES / name)           + cheminformatics                (FROZEN core)

  ┌───────────────┐   descriptors   ┌───────────────────┐   ┌────────────────┐
  │  RDKit         │───────────────►│ stokes_einstein    │   │                 │
  │  descriptors   │   MW, logP,     │   free_diffusivity │   │  physics/       │
  │  + library     │   TPSA, HBD…    ├───────────────────┤   │  diffusion,     │
  └───────────────┘                 │ brownian (particle │   │  pkpd …         │
          │                          │   sim) hindrance    │   │                 │
          │                          ├───────────────────┤   │  consumes       │
          │                          │ permeability / bbb │──►│  DrugParams     │
          └─────────────────────────►│ binding / clearance│   │  unchanged      │
                                     └─────────┬─────────┘   └────────────────┘
                                               │ compile_drug.build_profile
                                               ▼
                                     ┌───────────────────┐
                                     │  MolecularProfile  │  (rich, w/ provenance)
                                     │        └──►  core.DrugParams  (bridge)
                                     └───────────────────┘

   structure ──►  molecular sim ──►  DrugParams  ──►  tissue engine
```

`MolecularProfile` carries the full physics output plus a `provenance` dict
(method + confidence per quantity); `molecule_to_drugparams` narrows it to the
six numbers the frozen engine consumes.

---

## 2. Module reference

All modules live under `tmesim/molecular/`. Signatures are summarized from the
contract; shared types are in [`types.py`](../tmesim/molecular/types.py) and must
not be modified.

### `types.py` — shared data objects

- **`Molecule(name, smiles, descriptors: dict, source)`** with `.get(key,
  default)`. Descriptor keys: `mw` (g/mol), `logp`, `tpsa` (Å²), `hbd`, `hba`,
  `rotatable_bonds`, `aromatic_rings`, `formal_charge`, `radius_hydrodynamic_nm`
  (nm), `is_antibody`. `source` is `"rdkit"` or `"fallback_table"`.
- **`MolecularProfile(...)`** with `.summary()` — the full per-drug, per-context
  physics output: `free_diffusivity`, `effective_diffusivity`, `hindrance`,
  `hindrance_curve`, `membrane_permeability`, `bbb_penetrance`, `free_fraction`,
  `tissue_binding_rate`, `clearance`, `ic50`, `ecm_reference`, and `provenance`.

### `descriptors.py` — structure → numbers

- **`compute_descriptors(smiles) -> dict`** — RDKit computes `mw, logp, tpsa,
  hbd, hba, rotatable_bonds, aromatic_rings, formal_charge, is_antibody`
  (`False` for small molecules). Raises `ValueError` on an unparseable SMILES.
- **`from_smiles(name, smiles) -> Molecule`** — builds a `Molecule`, also filling
  `radius_hydrodynamic_nm` from MW (via the `stokes_einstein` relation). **Degrades
  gracefully**: if RDKit fails to parse, it fills from a small internal descriptor
  table keyed by `name` and sets `source="fallback_table"`.

**Physics/relation:** none of its own — it is the cheminformatics front door.
**In:** SMILES (+ name). **Out:** a populated `Molecule`. **Units:** as listed
above.

### `library.py` — the GBM drug set

- **`GBM_DRUGS: dict[str, str]`** — name → SMILES for real GBM-relevant agents:
  temozolomide, lomustine (CCNU), carmustine (BCNU), procarbazine, methotrexate,
  doxorubicin, cisplatin, etoposide, vincristine, irinotecan (a topoisomerase-I
  inhibitor prodrug; note its compiled DrugParams describe the parent molecule,
  not its active metabolite SN-38 — see `molecular/potency.py`), plus one
  macromolecule marker (`bevacizumab`, `is_antibody=True`, MW ≈ 149 000) to
  anchor the *size extreme*.
- **`get_molecule(name) -> Molecule`** (KeyError-safe listing of valid names),
  **`all_molecules() -> list[Molecule]`** — both go through
  `descriptors.from_smiles`.

**In:** a drug name. **Out:** `Molecule`(s). The antibody entry is deliberate: it
forces the whole pipeline (radius, diffusivity, permeability, BBB) to span from a
~0.4 nm small molecule to a ~5 nm antibody so the physics can be seen to behave
across three orders of magnitude in size.

### `stokes_einstein.py` — first-principles free diffusivity *(implemented)*

- **`hydrodynamic_radius_nm(mw, is_antibody=False) -> float`** — empirical MW →
  Stokes radius, two documented regimes:
  - Small molecules: `r_nm ≈ 0.066 · mw^(1/3)` (mass≈volume, radius≈volume^⅓;
    e.g. MW 300 → r ≈ 0.44 nm).
  - Globular proteins/antibodies: `r_nm ≈ 0.0515 · mw^0.392` (e.g. IgG MW 150 000
    → r ≈ 5.3 nm). The protein branch is forced for `is_antibody` **and** as a
    safety net for any `mw ≥ 5000`.
- **`free_diffusivity(radius_nm, temperature_K=310.0, viscosity_pa_s=1.2e-3) ->
  float`** — the **Stokes–Einstein** relation

  ```
  D = kB · T / (6 · π · η · r)
  ```

  with `kB = 1.380649×10⁻²³ J/K`. Computed in SI (radius nm → m) then reported in
  **µm²/s** (× 10¹², since 1 m² = 10¹² µm²). Dimensional check is carried in the
  docstring. The self-test verifies `D ∝ 1/r` and `D ∝ 1/η` exactly, and the
  small-molecule (~400–500 µm²/s) vs antibody (~1–10 µm²/s) separation.

**Physics:** genuine first-principles hydrodynamics (a sphere in a viscous
continuum). **In:** MW → radius → **Out:** dilute-fluid diffusivity, µm²/s. This
is the *free* diffusivity — no crowding yet.

### `brownian.py` — effective diffusivity by particle simulation

- **`simulate_effective_diffusivity(free_D_um2s, obstacle_fraction,
  particle_radius_nm=0.5, box_um=20.0, n_particles=400, steps=4000, dt_s=None,
  seed=0) -> dict`** — overdamped Langevin/Brownian dynamics of many
  non-interacting tracer particles through a periodic box packed with randomly
  placed spherical obstacles at volume fraction `obstacle_fraction`
  (excluded-volume: moves into obstacles are rejected/reflected). Tracks
  mean-squared displacement (MSD) on **unwrapped** positions and returns
  `{effective_D, hindrance = effective_D/free_D, msd, times, obstacle_fraction}`.
  See §3.
- **`hindrance_vs_ecm(free_D_um2s, ecm_fractions=(0.0…0.6), **kw) -> dict`** —
  sweeps obstacle fraction, returning `{ecm_fraction: effective_D}` for the
  `hindrance_curve`.

**Physics:** genuine Brownian dynamics. **In:** free `D` (µm²/s) + obstacle
volume fraction. **Out:** effective `D` (µm²/s) and dimensionless hindrance ∈ [0,1].

### `permeability.py` — passive membrane / BBB crossing

- **`membrane_permeability(descriptors, temperature_K=310.0) -> float`** — passive
  permeability (µm/s) via a **solubility–diffusion** model: `P ≈ Kp · Dm / d`,
  with the membrane partition coefficient `Kp` from logP, the in-membrane
  diffusivity `Dm` from size, and bilayer thickness `d ≈ 4 nm`.
- **`bbb_penetrance(descriptors) -> float`** — a `[0,1]` estimate from MW, logP,
  TPSA, HBD (physicochemical BBB rules; antibodies ≈ 0). Contractual ordering
  sanity: **temozolomide (crosses) > cisplatin / doxorubicin (poor) > antibody
  (≈ 0)**.

**Physics:** a mechanistic *model* (solubility–diffusion) parameterized by
cheminformatic descriptors — real physics form, estimated inputs. **In:**
descriptors. **Out:** `P` (µm/s), `bbb_penetrance` (dimensionless [0,1]).

### `binding.py` — protein binding, retention, clearance

- **`plasma_free_fraction(descriptors) -> float`** — unbound fraction ∈ [0,1]
  (lower for very lipophilic drugs).
- **`tissue_binding_rate(descriptors) -> float`** — first-order tissue
  retention/binding, 1/s (feeds `DrugParams.binding`).
- **`clearance_rate(descriptors) -> float`** — first-order clearance/metabolism,
  1/s (feeds `DrugParams.clearance`).

**Physics:** cheminformatic QSAR-style estimates, not derived from structure by
first principles. **In:** descriptors. **Out:** dimensionless free fraction; rates
in 1/s.

### `compile_drug.py` — the bridge

- **`build_profile(molecule, ecm_reference=0.4, ic50=None, temperature_K=310.0,
  brownian_kw=None) -> MolecularProfile`** — runs the full pipeline: Stokes–Einstein
  `free_diffusivity` → Brownian `hindrance` at `ecm_reference` (plus a small
  `hindrance_curve`) → `permeability` + `bbb_penetrance` →
  `binding`/`free_fraction`/`clearance`. Uses a provided `ic50`, else a
  **flagged rough estimate** recorded as low-confidence in `provenance`.
- **`molecule_to_drugparams(molecule_or_profile, ecm_reference=0.4, ic50=None,
  **kw) -> core.DrugParams`** — narrows the profile to the frozen engine's
  fields: `diffusivity = effective_diffusivity`, `molecular_weight = mw`, `ic50`,
  `bbb_penetrance`, `binding = tissue_binding_rate`, `clearance`; `hill` defaults
  to 1.2.

**In:** a `Molecule` (+ optional IC50 / ECM context). **Out:** a
`MolecularProfile` and a drop-in `core.DrugParams`. This is the seam where the
molecular layer hands off to the frozen tissue engine.

---

## 3. The Brownian-dynamics method

This is the layer's one piece of genuine, from-scratch transport physics, so it
is worth stating precisely.

### 3.1 Overdamped Langevin dynamics

At the scale of a solute in interstitial fluid, inertia is irrelevant: the
Reynolds number is tiny and momentum relaxes far faster than position changes.
The correct description is the **overdamped Langevin** (Brownian) equation, in
which the friction term dominates and the acceleration term is dropped:

```
dx/dt  =  (1/γ) F(x)  +  √(2 D) · ξ(t)
```

with γ the drag, `F` any external/excluded-volume force, `D` the free diffusivity
from Stokes–Einstein, and `ξ(t)` unit-variance Gaussian white noise. For free
tracers (`F = 0`) this discretizes to a Gaussian random walk whose per-step
displacement in each dimension has standard deviation

```
σ  =  √(2 · D · dt).
```

This √(2 D dt) step is the heart of the simulation. It is the discrete statement
that variance grows linearly with time at rate 2D per dimension — i.e. genuine
diffusion, not a ballistic (∝ t) or anomalous walk.

### 3.2 Excluded-volume obstacles

Crowding enters as geometry. The periodic box is packed with randomly placed
spherical obstacles at volume fraction `obstacle_fraction` (the ECM/cell fraction
the drug cannot enter). Each proposed particle move is checked against the
obstacles; a move that would place the tracer inside an obstacle is
rejected/reflected (excluded volume). Tracers therefore detour around obstacles,
travelling a longer, more convoluted path for the same net displacement — the
microscopic origin of hindered diffusion.

### 3.3 MSD → D estimator

The observable is the **mean-squared displacement** computed on *unwrapped*
positions (so periodic wrapping does not artificially cap displacement):

```
MSD(t)  =  ⟨ |x(t) − x(0)|² ⟩   averaged over all tracers.
```

For 2D diffusion, `MSD(t) = 4 · D · t`, so the effective diffusivity is the slope
of MSD vs time divided by 4:

```
effective_D  =  slope( MSD vs t ) / 4        (2D).
```

The contract's self-test checks that this recovers real diffusion: (1) MSD linear
in time at long times (diffusive, not ballistic or sub-diffusive); (2) with **zero
obstacles**, `effective_D ≈ free_D` to within a few percent — the estimator is
unbiased; (3) hindrance decreases monotonically with obstacle fraction; and
(4) a dilute-limit sanity check, `D_eff/D₀ ≈ 1 − a·φ` for small volume fraction
φ, with the fitted slope reported.

### 3.4 Hindrance / tortuosity, and how it feeds the tissue engine

The dimensionless ratio

```
hindrance  =  effective_D / free_D  ∈  [0, 1]
```

is the fraction of the dilute-fluid diffusivity that survives crowding. It is the
inverse-square of the classical **tortuosity** λ (`hindrance = 1/λ²`): tortuosity
measures how much longer the actual diffusive path is than the straight-line
distance because the solute must wind around impenetrable matrix and cells. A
hindrance of 0.5 means the drug diffuses at half its free rate through that
tissue.

This is exactly the quantity the tissue engine needs. The frozen engine's
transport PDE ([`architecture.md` §4.1](architecture.md)) uses an effective
diffusivity `D_eff = D · (1 − α · ecm_density)` — a *linear phenomenological*
tortuosity knob. The molecular layer replaces the hand-picked `D` with a
**mechanistically derived** `effective_diffusivity`: Stokes–Einstein gives the
free `D` from the molecule's hydrodynamic radius, and the Brownian simulation at
`ecm_reference` gives the crowding hindrance for that radius in *that* ECM
density. The result flows into `DrugParams.diffusivity`.

> **This chain has been measured against 71 compounds. Read this before quoting
> anything in §4.** The radius is `0.066 · mw**(1/3)`, so the Stokes–Einstein and
> Brownian legs take **molecular weight** as their only solute input. Until
> 2026-08-13 that was the whole chain, and it was **FALSIFIED**: scored blind and
> unfitted against Pruijn 2008 (71 compounds, HT29 multicellular layers, one
> method, one cell line) the measurement spans **29.7×** where the engine spanned
> **2.4×**, at **R² = 0.022**; at *identical* molecular weight the engine was flat
> **by construction — exactly 1.00×** — where the measurement spans **3.7×** (four
> compounds at 263 Da) and **3.8×** (four at 305 Da), so it ordered **0 of 31**
> exact-mass pairs.
>
> The two dominant predictors of the measurement are hydrogen-bond donor count
> (R² = 0.305) and logP₇.₄ (R² = 0.290) — both computed here and, at the time,
> discarded. A **two-route** cell term (paracellular between cells + transcellular
> through membranes) brings them in through a membrane permeability, and it is now
> applied on every compile. On the same 71 compounds, blind and unfitted, the
> composed chain reaches **R² = 0.718, Spearman +0.821**, with **99 %** of
> compounds within 3× of measurement, a geometric bias of **1.03×**, and **24 of
> 31** exact-mass pairs ordered correctly.
>
> **"Mechanistically derived" is now true of the form *and* of the content — but
> the content is only about a quarter right.** The model reproduces **8.3×** of
> the measured **29.7×** spread, so a similar-mass series is *under-separated*
> systematically rather than not separated at all; the entire test is **one cell
> line**; and absolute transport did not improve — tirapazamine's compiled `D_eff`
> now misses its measured multicellular-layer band by 1.5× (49.3 µm²/s against a
> measured 74–130), having previously missed it by 2.2× the other way. The
> practical consequence for a user: **read a ranking over similar-mass analogues
> as a hypothesis worth testing, not as a transport measurement.**
> → [`../results/two_route_promotion.md`](../results/two_route_promotion.md),
> [`../results/diffusivity_falsification.md`](../results/diffusivity_falsification.md)

---

## 4. What is real vs. estimated

The layer mixes three very different epistemic grades of number. Keeping them
separate is the whole honesty argument, so they are labeled explicitly in
`MolecularProfile.provenance` and summarized here.

| Quantity | Module / function | Grade | Basis |
|---|---|---|---|
| Free diffusivity `free_diffusivity` | `stokes_einstein.free_diffusivity` | **First-principles physics** | Stokes–Einstein `D = kB T / (6πηr)` — exact hydrodynamics for a sphere in a continuum. |
| Effective diffusivity / hindrance `effective_diffusivity`, `hindrance` | `brownian.simulate_effective_diffusivity` × `ecm_mesh.ogston_renkin_hindrance` × `two_route.two_route_diffusivity` | **PARTIALLY SUPPORTED against measurement** | Overdamped Langevin dynamics with excluded-volume obstacles, an Ogston–Renkin fibre-mesh gate, and a two-route permeable-cell effective-medium term (Hasselman & Johnson 1987). Scored blind on 71 measured compounds: **R² = 0.718**, 99 % within 3×, 24 of 31 exact-mass pairs ordered. The MW-only predecessor scored R² = 0.022 and ordered 0 of 31. Still recovers only 8.3× of a measured 29.7× spread, from one cell line, and the ABSOLUTE value still misses its measured band by 1.5×. |
| Hydrodynamic radius `radius_hydrodynamic_nm` | `stokes_einstein.hydrodynamic_radius_nm` | **Empirical scaling** | Documented MW→radius power laws (small-molecule cube-root; globular-protein fit). Well-validated but empirical. |
| Membrane permeability `membrane_permeability` | `permeability.membrane_permeability` | **Cheminformatic model** | Solubility–diffusion model (`P ≈ Kp·Dm/d`) — real physical *form*, descriptor-estimated inputs. |
| BBB penetrance `bbb_penetrance` | `permeability.bbb_penetrance` | **Cheminformatic estimate** | Physicochemical BBB rules on MW/logP/TPSA/HBD. Structure-only; ordering-calibrated, not magnitude-calibrated. |
| Free fraction / binding / clearance `free_fraction`, `tissue_binding_rate`, `clearance` | `binding.*` | **Cheminformatic estimate** | QSAR-style descriptor heuristics. Directionally reasonable, not fit to data. |
| Potency `ic50` | provided, or `compile_drug` flagged estimate | **Input / low-confidence** | Not derivable from structure without docking or assay. Treated as an input; any structure-only value is explicitly flagged low-confidence in `provenance`. |

The takeaway is deliberately unglamorous, and **measurement has since made it
less glamorous still**. This paragraph used to say that two quantities — free
diffusivity and the Brownian hindrance behind effective diffusivity — "are real
physics". The mechanisms are; the composed prediction is not. `effective_diffusivity`
is the **one quantity in this table that has been scored against a large external
dataset, and it failed** — while the permeability/BBB/binding rows it was being
favourably compared against were never tested that hard. So: the grades above
rank *derivational purity*, and derivational purity turned out not to predict
which row survives contact with data. Read the PARTIALLY SUPPORTED row as the
load-bearing one — it is the only row that has been scored against a large
external dataset at all, it was FALSIFIED on first contact, and what replaced it
is measured to be right about a quarter of the spread rather than right. Potency
is essentially an input. We do not claim otherwise.

---

## 5. Honest limitations

These are load-bearing caveats. A due-diligence reader should weigh them before
treating any molecular-layer number as calibrated.

- **Structure-only BBB under-predicts real CNS drugs.** Physicochemical BBB rules
  score passive permeability from MW/logP/TPSA/HBD. They will *under*-predict
  drugs that get into the brain by mechanisms the descriptors cannot see —
  **temozolomide** is the canonical example: it is a small, moderately polar
  molecule that crosses well clinically, partly aided by its chemistry and dosing,
  and a naive rules score can rank it lower than its real CNS exposure warrants.
  The contract pins *ordering* (TMZ > cisplatin/doxorubicin > antibody), not
  absolute magnitude.
- **Active transport and efflux: implemented, opt-in, default-OFF.** This
  previously read "there is no P-glycoprotein (P-gp) or other efflux/influx
  transporter model", and that was **wrong**. `tmesim/molecular/efflux.py`
  implements P-gp / BCRP from literature efflux ratios, and
  `tmesim/molecular/ion_trapping.py` implements ion trapping; both are scored
  (`results/efflux_impact.md`, `results/ion_trapping_impact.md`). Default-off is
  the *validated* recommendation, not an admission of absence: turning efflux on
  does not improve the measured anchor, and ABCB1/ABCG2 are on the explicit
  ABSENT list for both Xenium panels, so nothing spatial constrains it here.
  **As shipped (default-off) the consequence below still holds**: a passive-only
  solve will *over*-predict tissue exposure for P-gp substrates and does not
  represent transporter-mediated resistance.
- **The micro-sim is 2D.** The Brownian dynamics run in a 2D box (MSD = 4Dt).
  Real crowding, percolation thresholds, and obstacle connectivity differ between
  2D and 3D; the 2D hindrance is a defensible proxy, not the true 3D value. (This
  mirrors, and compounds with, the tissue engine's own 2D dodge in
  [`architecture.md` §6](architecture.md).)
- **No metabolism networks.** Clearance is a single first-order rate, not a
  metabolic pathway. Active metabolites, saturable enzymes (CYP kinetics), and
  drug–drug interactions are out of scope.
- **Potency is not derivable from structure.** IC50 is an input or a flagged rough
  estimate. Deriving it properly needs molecular docking against the actual target
  or an experimental assay — neither is in this layer. Treat any structure-only
  IC50 as a placeholder, not a prediction.

### The key point

Even with all of the above, the layer buys something real. Before it, the tissue
engine's per-drug parameters (diffusivity, BBB penetrance, binding, clearance)
were **hand-set** order-of-magnitude guesses. After it, those same parameters are
derived from the drug's structure — free diffusivity from Stokes–Einstein,
effective diffusivity from a real Brownian simulation, the rest from documented
cheminformatic models with explicit provenance. That is a real upgrade from
"assumed" to "computed".

**But not, for diffusivity, to "computed from the molecule."** That phrase stood
here and it is measurably wrong: `effective_diffusivity` is computed from the
molecule's **mass** alone, and at matched mass it is flat by construction while
the measurement moves 3.7× (see the box in §3). The upgrade is real for the
descriptor-consuming quantities; for transport it bought a better-founded
*form* around the same single input.

It is **not**, however, a substitute for validation. Mechanistic ≠ correct — and
this is no longer a hypothetical. This paragraph used to call the comparison
against measured tissue PK "the necessary next step". **That step has been
taken, and the layer lost — then partially recovered.** `effective_diffusivity`
was scored against 71 measured compounds and falsified (§3); the two-route model
that replaced it scores R² = 0.718 on the same data and brought the
multicellular-layer **hindrance** anchor inside its measured band for the first
time (FAIL → PASS, 1.77-fold → 10.33-fold against a measured 12-fold). The
multicellular-layer **absolute `D_eff`** anchor still FAILS in `vivum validate`,
with its miss now on the low side rather than the high; and `bbb_penetrance`
still fails its rank ordering against measured brain:plasma ratios. What remains genuinely untested is
the *solved field*: nothing public measures a spatially resolved free-interstitial
concentration in tumour tissue. So the honest position is not "unvalidated"
but "the input coefficients are testable and currently fail in identifiable,
structural ways", which is a sharper and less comfortable statement than the one
this paragraph used to make. Treat these numbers as less than a defensible,
first-principles-
anchored estimate.
