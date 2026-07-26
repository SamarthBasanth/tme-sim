"""The multiscale BRIDGE: compile a real drug MOLECULE into the tissue engine's
``core.DrugParams``.

This module runs the whole molecular pipeline end-to-end and assembles the two
outputs the contract asks for:

  molecular structure
      -> Stokes-Einstein free diffusivity          (stokes_einstein.py)
      -> Brownian-dynamics hindrance in crowded ECM (brownian.py)
      -> passive permeability + BBB penetrance      (permeability.py)
      -> plasma/tissue binding + clearance          (binding.py)
      => MolecularProfile  (rich, with provenance)
      => core.DrugParams   (drop-in for the FROZEN tissue engine)

Nothing here modifies ``tmesim.core``, the tissue engine, or ``types.py``.

Honesty line (carried from MOLECULAR_CONTRACT.md): diffusivity/permeability/BBB
are physics- and cheminformatics-based estimates from structure; potency (IC50)
is the weakest link and is treated as an INPUT or a clearly-flagged low-confidence
default. Every quantity in the profile carries a provenance string (method +
confidence) so each number can be defended.

Units: length um (nm where noted), time s, diffusivity um^2/s, concentration uM,
temperature K.
"""
from __future__ import annotations

from typing import Optional, Union

from ..core import DrugParams
from . import binding, permeability, stokes_einstein
from . import brownian
from . import ecm_mesh, potency
from .types import Molecule, MolecularProfile

# ---------------------------------------------------------------------------
# Simulation-cost defaults. A full compile runs one Brownian sim at the ECM
# reference PLUS a 4-point coarse hindrance curve = 5 sims. We keep each sim
# small and deterministic so a full compile finishes in a few seconds. Callers
# can override any of these via `brownian_kw`.
# ---------------------------------------------------------------------------
# Main effective-D sim at the ECM reference (the number that feeds DrugParams).
_MAIN_SIM_KW = {"n_particles": 300, "steps": 2500, "seed": 0}
# Coarse curve: cheaper still (fewer particles/steps) since it is only for plots.
_CURVE_SIM_KW = {"n_particles": 200, "steps": 1500, "seed": 0}
# Fractions sampled for the coarse hindrance curve (brackets the ECM reference).
_CURVE_FRACTIONS = (0.0, 0.2, 0.4, 0.6)

# Rough default potency when none is supplied. Flagged LOW confidence: structure
# alone cannot give potency without docking/assay. 25 uM matches the order of
# magnitude of the engine's temozolomide-like default (core.DRUG_LIBRARY).
_DEFAULT_IC50_UM = 25.0

# Default Hill coefficient handed to DrugParams (matches the engine's TMZ-like row).
_DEFAULT_HILL = 1.2


def build_profile(
    molecule: Molecule,
    ecm_reference: float = 0.4,
    ic50: Optional[float] = None,
    temperature_K: float = 310.0,
    brownian_kw: Optional[dict] = None,
) -> MolecularProfile:
    """Run the full molecular pipeline for one drug in one ECM context.

    Pipeline
    --------
    1. hydrodynamic radius: prefer the precomputed ``radius_hydrodynamic_nm``
       descriptor; else derive it from MW via
       ``stokes_einstein.hydrodynamic_radius_nm(mw, is_antibody)``.
    2. free diffusivity: ``stokes_einstein.free_diffusivity(radius, T)`` (um^2/s).
    3. effective diffusivity: a genuine 2D overdamped Brownian-dynamics sim through
       an obstacle field at ``ecm_reference`` volume fraction
       (``brownian.simulate_effective_diffusivity``). The tracer radius is set to
       the molecule's hydrodynamic radius so larger drugs are hindered more.
       ``hindrance = effective_D / free_D`` in [0, 1].
    4. a COARSE ``hindrance_curve`` over obstacle fractions (0, 0.2, 0.4, 0.6),
       run at reduced particle/step counts so the whole compile stays fast.
    5. passive ``membrane_permeability`` (um/s) and ``bbb_penetrance`` [0,1]
       (``permeability``).
    6. ``free_fraction`` [0,1], ``tissue_binding_rate`` (1/s), ``clearance`` (1/s)
       (``binding``).
    7. IC50: use ``ic50`` if provided; else the flagged low-confidence default.

    Parameters
    ----------
    molecule : Molecule
        Drug with physicochemical descriptors (from molecular.descriptors).
    ecm_reference : float
        ECM/obstacle volume fraction the effective diffusivity is computed at.
    ic50 : float or None
        Potency in uM. If None, a rough default is used and flagged in provenance.
    temperature_K : float
        Absolute temperature (default 310 K = body temperature).
    brownian_kw : dict or None
        Extra kwargs forwarded to the Brownian sims (e.g. n_particles, steps,
        seed, box_um). Overrides the internal cost defaults.

    Returns
    -------
    MolecularProfile
        Full physics output with a provenance entry for every quantity.
    """
    descriptors = molecule.descriptors
    mw = float(molecule.get("mw", 0.0) or 0.0)
    is_antibody = bool(molecule.get("is_antibody", False))

    # --- 1. hydrodynamic radius (nm) ---
    radius_nm = molecule.get("radius_hydrodynamic_nm")
    if radius_nm:
        radius_nm = float(radius_nm)
        radius_src = "molecule.radius_hydrodynamic_nm (precomputed MW-scaling)"
    else:
        radius_nm = stokes_einstein.hydrodynamic_radius_nm(mw, is_antibody)
        radius_src = "stokes_einstein.hydrodynamic_radius_nm(mw, is_antibody)"

    # --- 2. free (dilute-fluid) diffusivity, um^2/s ---
    free_D = stokes_einstein.free_diffusivity(radius_nm, temperature_K=temperature_K)

    # --- 3. effective diffusivity via Brownian dynamics at ecm_reference ---
    # Tracer radius = molecule hydrodynamic radius (bigger drug -> bigger excluded
    # region -> more hindrance), unless the caller overrides it.
    main_kw = dict(_MAIN_SIM_KW)
    main_kw.setdefault("particle_radius_nm", radius_nm)
    if brownian_kw:
        main_kw.update(brownian_kw)
    sim = brownian.simulate_effective_diffusivity(free_D, ecm_reference, **main_kw)
    micro_D = float(sim["effective_D"])  # tortuosity from ~300 nm obstacles

    # --- 3b. nm-scale ECM mesh size gate (Ogston-Renkin fiber matrix) ---
    # The 300 nm micro-obstacles capture tortuosity but barely resolve drug SIZE
    # (a 0.4 nm drug and a 5.5 nm antibody see nearly the same field). The ECM is
    # really a ~few-nm fiber mesh, where solute radius strongly gates diffusion.
    # We multiply the micro-tortuosity by the fiber-matrix effective hindrance so
    # size differences (esp. small molecule vs antibody) show up mechanistically.
    mesh = ecm_mesh.ogston_renkin_hindrance(radius_nm)
    mesh_h = float(mesh["effective_hindrance"])
    effective_D = micro_D * mesh_h
    hindrance = effective_D / free_D if free_D > 0 else 0.0

    # --- 4. coarse hindrance curve (reduced cost) ---
    curve_kw = dict(_CURVE_SIM_KW)
    curve_kw.setdefault("particle_radius_nm", radius_nm)
    if brownian_kw:
        curve_kw.update(brownian_kw)
    hindrance_curve = brownian.hindrance_vs_ecm(
        free_D, ecm_fractions=_CURVE_FRACTIONS, **curve_kw
    )

    # --- 5. passive permeability + BBB penetrance ---
    perm = permeability.membrane_permeability(descriptors, temperature_K=temperature_K)
    bbb = permeability.bbb_penetrance(descriptors)

    # --- 6. plasma/tissue binding + clearance ---
    free_fraction = binding.plasma_free_fraction(descriptors)
    tissue_binding = binding.tissue_binding_rate(descriptors)
    clearance = binding.clearance_rate(descriptors)

    # --- 7. IC50: input > measured (GDSC) > literature > flagged default ---
    if ic50 is not None:
        ic50_val = float(ic50)
        ic50_prov = "user-provided (input); confidence: as supplied"
    else:
        try:
            res = potency.resolve_ic50(molecule.name)
            ic50_val = float(res["ic50_uM"])
            ic50_prov = (
                f"{res.get('source')} ({res.get('confidence')})"
                + (f", n_cell_lines={res['n_cell_lines']}" if res.get("n_cell_lines") else "")
                + f"; {res.get('note', '')}".rstrip("; ")
            )
        except Exception as exc:  # potency lookup must never break a compile
            ic50_val = _DEFAULT_IC50_UM
            ic50_prov = (
                f"rough default {_DEFAULT_IC50_UM} uM; confidence: LOW -- potency "
                f"lookup failed ({exc!r}); structure alone cannot give potency."
            )

    provenance = {
        "hydrodynamic_radius_nm": f"{radius_src}; confidence: medium (empirical MW scaling)",
        "free_diffusivity": (
            "Stokes-Einstein D=kB*T/(6*pi*eta*r), interstitial fluid eta=1.2e-3 Pa.s; "
            "confidence: high (first-principles physics given the radius)"
        ),
        "effective_diffusivity": (
            f"micro-tortuosity from 2D Brownian dynamics at phi={ecm_reference} "
            f"(brownian, fit_r2={sim.get('fit_r2', float('nan')):.4f}, "
            f"D={micro_D:.1f}) x nm-scale ECM fiber-matrix size gate "
            f"(ecm_mesh.ogston_renkin_hindrance, mesh_h={mesh_h:.3f}); confidence: "
            "medium (genuine particle sim + Ogston-Renkin fiber theory; ECM "
            "geometry is a stylized proxy)"
        ),
        "hindrance": (
            "effective_D/free_D from the Brownian sim; confidence: medium"
        ),
        "hindrance_curve": (
            f"brownian.hindrance_vs_ecm at phi={_CURVE_FRACTIONS} "
            f"(reduced n_particles/steps for speed); confidence: medium (coarse)"
        ),
        "membrane_permeability": (
            "solubility-diffusion P=Kp*Dm/d from logP + size (permeability."
            "membrane_permeability); confidence: low-medium (passive only; ignores "
            "active transport / efflux)"
        ),
        "bbb_penetrance": (
            "CNS-MPO-style physicochemical rules on MW/logP/TPSA/HBD (permeability."
            "bbb_penetrance); confidence: low (structure-only; known to under-predict "
            "temozolomide, and cannot capture active transport)"
        ),
        "free_fraction": (
            "logistic in logP (albumin-binding proxy; binding.plasma_free_fraction); "
            "confidence: low (direction + rough magnitude only)"
        ),
        "tissue_binding_rate": (
            "logP + HBD stickiness heuristic (binding.tissue_binding_rate); "
            "confidence: low (order-of-magnitude prior)"
        ),
        "clearance": (
            "mild logP/MW heuristic (binding.clearance_rate); confidence: low "
            "(weak prior, not a metabolism prediction)"
        ),
        "ic50": ic50_prov,
    }

    return MolecularProfile(
        molecule=molecule,
        free_diffusivity=free_D,
        effective_diffusivity=effective_D,
        hindrance=hindrance,
        membrane_permeability=perm,
        bbb_penetrance=bbb,
        free_fraction=free_fraction,
        tissue_binding_rate=tissue_binding,
        clearance=clearance,
        ic50=ic50_val,
        ecm_reference=ecm_reference,
        hindrance_curve=hindrance_curve,
        provenance=provenance,
    )


def molecule_to_drugparams(
    molecule_or_profile: Union[Molecule, MolecularProfile],
    ecm_reference: float = 0.4,
    ic50: Optional[float] = None,
    **kw,
) -> DrugParams:
    """Compile a Molecule (or an already-built MolecularProfile) into DrugParams.

    If given a ``Molecule``, ``build_profile`` is run first (forwarding
    ``ecm_reference``, ``ic50``, and any extra kwargs such as ``temperature_K`` or
    ``brownian_kw``). If given a ``MolecularProfile``, it is used directly (and a
    supplied ``ic50`` overrides the profile's).

    Mapping to the FROZEN engine's DrugParams:
        diffusivity      = profile.effective_diffusivity   (hindered tissue D)
        molecular_weight = molecule MW
        ic50             = supplied ic50, else profile.ic50
        hill             = 1.2 (default)
        binding          = profile.tissue_binding_rate
        clearance        = profile.clearance
        bbb_penetrance   = profile.bbb_penetrance
    """
    if isinstance(molecule_or_profile, MolecularProfile):
        profile = molecule_or_profile
    elif isinstance(molecule_or_profile, Molecule):
        profile = build_profile(
            molecule_or_profile, ecm_reference=ecm_reference, ic50=ic50, **kw
        )
    else:
        raise TypeError(
            "molecule_or_profile must be a Molecule or MolecularProfile, got "
            f"{type(molecule_or_profile).__name__}"
        )

    molecule = profile.molecule
    mw = float(molecule.get("mw", 0.0) or 0.0)
    ic50_val = float(ic50) if ic50 is not None else float(profile.ic50)

    return DrugParams(
        name=molecule.name,
        diffusivity=profile.effective_diffusivity,
        molecular_weight=mw,
        ic50=ic50_val,
        hill=_DEFAULT_HILL,
        binding=profile.tissue_binding_rate,
        clearance=profile.clearance,
        bbb_penetrance=profile.bbb_penetrance,
    )


if __name__ == "__main__":
    import math
    import time as _time

    from . import library
    # Confirm the bridge target is the real, frozen engine type.
    from ..core import DrugParams as CoreDrugParams

    t0 = _time.time()
    print("=" * 70)
    print("compile_drug BRIDGE self-test: molecule -> MolecularProfile -> DrugParams")
    print("=" * 70)

    tmz_mol = library.get_molecule("temozolomide")
    dox_mol = library.get_molecule("doxorubicin")

    tmz = build_profile(tmz_mol)
    dox = build_profile(dox_mol)

    for name, prof in (("temozolomide", tmz), ("doxorubicin", dox)):
        print(f"\n--- {name} ---")
        print("summary:", prof.summary())
        print(f"  free_D     = {prof.free_diffusivity:8.2f} um^2/s")
        print(f"  eff_D      = {prof.effective_diffusivity:8.2f} um^2/s")
        print(f"  hindrance  = {prof.hindrance:8.3f}")
        print(f"  bbb        = {prof.bbb_penetrance:8.3f}")
        print(f"  ic50       = {prof.ic50} uM")
        print(f"  ic50 prov  : {prof.provenance['ic50']}")
        print(f"  curve      : "
              + ", ".join(f"phi={p}:{d:.1f}" for p, d in prof.hindrance_curve.items()))

    # --- MECHANISTIC assertions -------------------------------------------
    # Doxorubicin (MW ~544) is larger than temozolomide (MW ~194): larger drug
    # => smaller radius-inverse free_D AND heavier ECM hindrance => LOWER eff_D.
    assert dox_mol.get("mw") > tmz_mol.get("mw"), "doxorubicin should be heavier than TMZ"
    assert dox.free_diffusivity < tmz.free_diffusivity, (
        f"larger dox free_D {dox.free_diffusivity:.2f} should be < TMZ "
        f"{tmz.free_diffusivity:.2f}")
    assert dox.effective_diffusivity < tmz.effective_diffusivity, (
        f"larger dox eff_D {dox.effective_diffusivity:.2f} should be < TMZ "
        f"{tmz.effective_diffusivity:.2f}")
    # BBB: doxorubicin (large, very polar, many donors) crosses no better than TMZ.
    assert dox.bbb_penetrance <= tmz.bbb_penetrance, (
        f"dox bbb {dox.bbb_penetrance:.3f} should be <= TMZ {tmz.bbb_penetrance:.3f}")
    print("\n[mechanistic] dox: LOWER free_D, LOWER eff_D, LOWER-or-equal bbb than TMZ  OK")

    # --- Valid, finite, sane DrugParams from the engine's real type --------
    def _finite(x):
        return isinstance(x, (int, float)) and math.isfinite(x)

    for name, mol, prof in (("temozolomide", tmz_mol, tmz), ("doxorubicin", dox_mol, dox)):
        dp = molecule_to_drugparams(mol)
        # It must be the engine's frozen DrugParams type.
        assert isinstance(dp, CoreDrugParams), f"{name}: not a core.DrugParams"
        print(f"\nDrugParams({name}) = {dp}")
        for fld in ("diffusivity", "molecular_weight", "ic50", "hill",
                    "binding", "clearance", "bbb_penetrance"):
            v = getattr(dp, fld)
            assert _finite(v), f"{name}.{fld} not finite: {v}"
        # sane physical ranges
        assert 0.1 <= dp.diffusivity <= 1000.0, f"{name} diffusivity out of range: {dp.diffusivity}"
        assert 50.0 <= dp.molecular_weight <= 2000.0, f"{name} MW out of range: {dp.molecular_weight}"
        assert dp.ic50 > 0.0, f"{name} ic50 must be positive: {dp.ic50}"
        assert 0.0 <= dp.bbb_penetrance <= 1.0, f"{name} bbb out of [0,1]: {dp.bbb_penetrance}"
        assert 1e-5 <= dp.binding <= 5e-3, f"{name} binding out of band: {dp.binding}"
        assert 1e-5 <= dp.clearance <= 1e-3, f"{name} clearance out of band: {dp.clearance}"
        assert dp.name == name, f"name mismatch: {dp.name} != {name}"

    # Consistency: passing a prebuilt profile yields the same eff_D-derived diffusivity.
    dp_from_profile = molecule_to_drugparams(tmz)
    assert abs(dp_from_profile.diffusivity - tmz.effective_diffusivity) < 1e-9
    # Passing a Molecule builds a profile internally and also works.
    dp_from_mol = molecule_to_drugparams(tmz_mol)
    assert isinstance(dp_from_mol, CoreDrugParams)

    print(f"\nAll bridge self-tests passed. Runtime {_time.time() - t0:.2f} s")
