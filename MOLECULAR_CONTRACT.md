# tmesim.molecular — Molecular-Layer Contract (v1)

**Read fully before coding.** This adds a molecular-scale layer that turns a real
drug MOLECULE (structure) into the physical parameters the tissue engine needs.
The tissue engine and `tmesim.core` are FROZEN — do not modify them. The molecular
layer's final output is a `core.DrugParams` (drop-in for the existing engine) plus
a richer `MolecularProfile`.

## The multiscale idea (and the honesty line)
molecular structure  ──►  physics/particle sim  ──►  DrugParams  ──►  (existing tissue engine)

- This is NOT atomistic MD of a tissue (not computable). Two honest methods only:
  1. **Physicochemical / cheminformatics estimation** from structure (RDKit).
  2. **Genuine overdamped Langevin / Brownian dynamics** at the micro-scale to get
     an EFFECTIVE diffusivity in a crowded (ECM/cell) medium.
- Target potency (IC50) is NOT reliably derivable from structure without docking/
  assay — treat as an input or a clearly-flagged rough estimate. Be honest in
  `provenance`.

Units: length um (or nm where noted), time s, diffusivity um^2/s, concentration uM,
temperature K. Body temp = 310 K. Interstitial fluid viscosity ~ 1.2e-3 Pa·s.

## Shared types (already implemented in `tmesim.molecular.types`, do not modify)
- `Molecule(name, smiles, descriptors: dict, source)` — `.get(key, default)`.
  descriptor keys: mw, logp, tpsa, hbd, hba, rotatable_bonds, aromatic_rings,
  formal_charge, radius_hydrodynamic_nm, is_antibody.
- `MolecularProfile(molecule, free_diffusivity, effective_diffusivity, hindrance,
  membrane_permeability, bbb_penetrance, free_fraction, tissue_binding_rate,
  clearance, ic50, ecm_reference, hindrance_curve, provenance)` — `.summary()`.
- Bridge target: `tmesim.core.DrugParams(name, diffusivity, molecular_weight, ic50,
  hill, binding, clearance, bbb_penetrance)`.

RDKit IS installed (verified). Prefer it; but every module must DEGRADE GRACEFULLY
to a curated fallback if a SMILES fails to parse.

## Module APIs (implement exactly)

### `tmesim/molecular/descriptors.py`
```python
def compute_descriptors(smiles: str) -> dict:
    """RDKit -> {mw, logp, tpsa, hbd, hba, rotatable_bonds, aromatic_rings,
    formal_charge, is_antibody(False for small mols)}. Raises ValueError on bad SMILES."""
def from_smiles(name: str, smiles: str) -> Molecule:
    """Build a Molecule; also compute radius_hydrodynamic_nm from mw (see stokes_einstein
    helper or a local MW->radius relation) and store it in descriptors. Fallback: if
    RDKit parse fails, fill from a small internal descriptor table keyed by name and set
    source='fallback_table'."""
```

### `tmesim/molecular/library.py`
```python
GBM_DRUGS: dict[str, str]   # name -> SMILES, real GBM-relevant agents:
  # temozolomide, lomustine (CCNU), carmustine (BCNU), procarbazine, methotrexate,
  # doxorubicin, cisplatin, etoposide, vincristine. Include one macromolecule marker
  # e.g. 'bevacizumab' as an antibody (is_antibody=True, mw~149000) to show the size extreme.
def get_molecule(name: str) -> Molecule   # via descriptors.from_smiles; KeyError-safe listing names
def all_molecules() -> list[Molecule]
```

### `tmesim/molecular/stokes_einstein.py`
```python
def hydrodynamic_radius_nm(mw: float, is_antibody: bool=False) -> float:
    """Empirical MW->hydrodynamic radius (protein/small-molecule scaling). Document the relation."""
def free_diffusivity(radius_nm: float, temperature_K: float=310.0,
                     viscosity_pa_s: float=1.2e-3) -> float:
    """Stokes-Einstein D = kB*T/(6*pi*eta*r), returned in um^2/s. Document unit conversion."""
```

### `tmesim/molecular/brownian.py`  (GENUINE PARTICLE PHYSICS — be rigorous)
```python
def simulate_effective_diffusivity(free_D_um2s: float, obstacle_fraction: float,
                                   particle_radius_nm: float=0.5, box_um: float=20.0,
                                   n_particles: int=400, steps: int=4000, dt_s: float=None,
                                   seed: int=0) -> dict:
    """Overdamped Langevin / Brownian dynamics: many non-interacting tracer particles
    random-walk (step std = sqrt(2*free_D*dt)) through a periodic box packed with
    randomly-placed spherical obstacles at volume fraction `obstacle_fraction`
    (excluded volume: reject/reflect moves into obstacles). Track mean-squared
    displacement (MSD) using UNWRAPPED positions; effective_D = slope of MSD vs time /4
    (2D) . Return {effective_D, hindrance=effective_D/free_D, msd, times, obstacle_fraction}."""
def hindrance_vs_ecm(free_D_um2s: float, ecm_fractions=(0.0,0.1,0.2,0.3,0.4,0.5,0.6),
                     **kw) -> dict:
    """Run the sim across ECM/obstacle fractions; return {ecm_fraction: effective_D}."""
```
Self-test MUST verify: (1) MSD linear in time (diffusive, not ballistic/sub-diffusive
at long time); (2) with zero obstacles effective_D ~= free_D (within a few %);
(3) hindrance decreases monotonically with obstacle_fraction; (4) dilute-limit sanity:
for small fraction, hindrance approx matches a known relation (e.g. D_eff/D0 ~ 1 - a*phi;
report the fitted slope). Keep runtime to a few seconds (vectorize over particles).

### `tmesim/molecular/permeability.py`
```python
def membrane_permeability(descriptors: dict, temperature_K: float=310.0) -> float:
    """Passive permeability (um/s) via a solubility-diffusion model: P ~ Kp * Dm / d,
    partition Kp from logP, membrane diffusion Dm from size, d~4nm bilayer. Document."""
def bbb_penetrance(descriptors: dict) -> float:
    """[0,1] estimate from MW, logP, TPSA, HBD (physicochemical BBB rules; antibodies ~0).
    ORDERING must be sane: temozolomide (crosses) > cisplatin/doxorubicin (poor) > antibody(~0)."""
```

### `tmesim/molecular/binding.py`
```python
def plasma_free_fraction(descriptors: dict) -> float   # [0,1], lower for very lipophilic
def tissue_binding_rate(descriptors: dict) -> float     # 1/s retention (DrugParams.binding)
def clearance_rate(descriptors: dict) -> float          # 1/s metabolism/clearance
```

### `tmesim/molecular/compile_drug.py`  (THE BRIDGE — wave 2)
```python
def build_profile(molecule: Molecule, ecm_reference: float=0.4, ic50: float|None=None,
                  temperature_K: float=310.0, brownian_kw: dict|None=None) -> MolecularProfile:
    """Run the full molecular pipeline: SE free_diffusivity -> brownian hindrance at
    ecm_reference (and a small hindrance_curve) -> permeability + bbb_penetrance ->
    binding/free_fraction/clearance. ic50: use provided; else a flagged rough estimate
    (document as low-confidence in provenance). Assemble MolecularProfile with provenance."""
def molecule_to_drugparams(molecule_or_profile, ecm_reference: float=0.4,
                           ic50: float|None=None, **kw) -> DrugParams:
    """Return a core.DrugParams: diffusivity=effective_diffusivity, molecular_weight=mw,
    ic50, bbb_penetrance, binding=tissue_binding_rate, clearance. hill default 1.2."""
```

## Rules for every agent
- Create only your files; never modify `tmesim/core/*`, the tissue engine, or `types.py`.
- venv: `/Users/samarthbasanth/tme-sim/.venv/bin/python` (editable install; RDKit installed).
- Physics must be dimensionally correct with units in docstrings. Deterministic (seed).
- Vectorize (esp. brownian: step all particles at once). One call in a few seconds.
- Add `if __name__=='__main__':` self-test and RUN it with venv python; report results.
- If a tool call hits a network/API error, wait and retry (wifi has been intermittent).
