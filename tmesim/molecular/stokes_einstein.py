"""Stokes-Einstein free diffusivity and MW -> hydrodynamic-radius scaling.

This module gives the molecular layer its *dilute-fluid* diffusivity: the
diffusivity a molecule would have drifting through pure interstitial fluid with
no crowding. The Brownian-dynamics module (`brownian.py`) later hinders this
value to get the effective diffusivity through crowded ECM.

Two independent empirical relations live here:

1. MW -> hydrodynamic radius  (`hydrodynamic_radius_nm`)
2. radius -> Stokes-Einstein diffusivity  (`free_diffusivity`)

Everything is dimensionally explicit; every unit conversion is documented in the
docstrings and in inline comments so the numbers can be defended.

Units convention (matches the molecular contract):
  length      : nm in this module's radius API, converted to m internally for SI;
                diffusivity reported in um^2/s (tissue-engine convention).
  temperature : K   (body temperature default 310 K)
  viscosity   : Pa.s (interstitial fluid ~ 1.2e-3 Pa.s, close to water at 37 C)
"""
from __future__ import annotations

import math

# Boltzmann constant, SI: J/K = kg.m^2.s^-2.K^-1
_KB = 1.380649e-23  # J/K (CODATA exact)


def hydrodynamic_radius_nm(mw: float, is_antibody: bool = False) -> float:
    """Empirical molecular-weight -> hydrodynamic (Stokes) radius, in nm.

    Two documented regimes, both widely used empirical fits:

    Small molecules / general solutes
    ---------------------------------
        r_nm ~= 0.066 * mw**(1/3)          (mw in g/mol, r in nm)

    This is the classic "mass ~ volume, radius ~ volume^(1/3)" scaling. Assuming a
    roughly constant molecular density, volume scales linearly with MW, so the
    radius scales as MW^(1/3). The prefactor 0.066 nm.(g/mol)^-1/3 reproduces the
    ~0.4-0.5 nm radii measured for typical small drug molecules
    (e.g. mw=300 -> r ~= 0.066 * 300**(1/3) ~= 0.44 nm), consistent with
    small-molecule Stokes radii from the literature.

    Proteins / antibodies (globular)
    --------------------------------
        r_nm ~= 0.0515 * mw**0.392         (mw in g/mol, r in nm)

    A standard globular-protein hydrodynamic-radius fit (a power law with exponent
    ~0.39, close to but below 1/3 because proteins pack denser and are not perfect
    spheres). It gives, e.g., a ~150 kDa antibody (mw=150000):
        r ~= 0.0515 * 150000**0.392 ~= 5.3 nm,
    matching the well-known ~5 nm hydrodynamic radius of IgG antibodies.

    The `is_antibody` flag forces the protein branch. As a safety net, any very
    large molecule (mw >= 5000 g/mol) also uses the protein branch, because the
    small-molecule cube-root fit under-predicts radii badly for macromolecules.

    Parameters
    ----------
    mw : float
        Molecular weight in g/mol (== Da).
    is_antibody : bool
        If True, use the globular-protein scaling regardless of mw.

    Returns
    -------
    float
        Hydrodynamic radius in nm.
    """
    if mw <= 0:
        raise ValueError(f"mw must be positive, got {mw!r}")

    use_protein_scaling = is_antibody or mw >= 5000.0
    if use_protein_scaling:
        # Globular-protein power law: r[nm] = 0.0515 * mw^0.392
        return 0.0515 * mw ** 0.392

    # Small-molecule cube-root law: r[nm] = 0.066 * mw^(1/3)
    return 0.066 * mw ** (1.0 / 3.0)


def free_diffusivity(
    radius_nm: float,
    temperature_K: float = 310.0,
    viscosity_pa_s: float = 1.2e-3,
) -> float:
    """Stokes-Einstein free diffusivity of a sphere, returned in um^2/s.

    Physics
    -------
        D = kB * T / (6 * pi * eta * r)

    with
        kB  = 1.380649e-23 J/K      (Boltzmann constant)
        T   = temperature in K
        eta = dynamic viscosity in Pa.s
        r   = hydrodynamic radius in m

    Dimensional check (SI):
        [kB*T]        = J = kg.m^2.s^-2
        [6*pi*eta*r]  = (Pa.s).m = (kg.m^-1.s^-2 . s).m = kg.s^-1
        [D]           = (kg.m^2.s^-2) / (kg.s^-1) = m^2/s   OK.

    Unit conversions, step by step
    ------------------------------
      1. radius:  r_m  = radius_nm * 1e-9         (nm -> m)
      2. compute D_SI = kB*T / (6*pi*eta*r_m)     (result in m^2/s)
      3. report:  D_um2s = D_SI * 1e12            (m^2/s -> um^2/s)
                  because 1 m^2 = (1e6 um)^2 = 1e12 um^2, and seconds are unchanged.

    Parameters
    ----------
    radius_nm : float
        Hydrodynamic radius in nm.
    temperature_K : float
        Absolute temperature in K (default 310 K = body temperature).
    viscosity_pa_s : float
        Dynamic viscosity in Pa.s (default 1.2e-3, interstitial fluid at 37 C).

    Returns
    -------
    float
        Free (dilute-fluid) diffusivity in um^2/s.
    """
    if radius_nm <= 0:
        raise ValueError(f"radius_nm must be positive, got {radius_nm!r}")
    if temperature_K <= 0:
        raise ValueError(f"temperature_K must be positive, got {temperature_K!r}")
    if viscosity_pa_s <= 0:
        raise ValueError(f"viscosity_pa_s must be positive, got {viscosity_pa_s!r}")

    # Step 1: nm -> m
    r_m = radius_nm * 1e-9

    # Step 2: Stokes-Einstein in SI -> m^2/s
    d_si = _KB * temperature_K / (6.0 * math.pi * viscosity_pa_s * r_m)

    # Step 3: m^2/s -> um^2/s  (x 1e12)
    return d_si * 1e12


if __name__ == "__main__":
    print("=== stokes_einstein self-test ===")

    # --- Small molecule (~300 g/mol) ---
    mw_small = 300.0
    r_small = hydrodynamic_radius_nm(mw_small)
    D_small = free_diffusivity(r_small)
    print(f"small molecule  mw={mw_small:>8.0f} g/mol  "
          f"r={r_small:.3f} nm  D={D_small:.1f} um^2/s")

    # radius sanity: ~0.44 nm as documented
    assert 0.40 <= r_small <= 0.50, f"small-molecule radius off: {r_small}"
    # D sanity vs literature small-molecule tissue D ~ hundreds of um^2/s
    assert 300.0 <= D_small <= 600.0, f"small-molecule D out of range: {D_small}"

    # --- Antibody (~150000 g/mol) ---
    mw_ab = 150000.0
    r_ab = hydrodynamic_radius_nm(mw_ab, is_antibody=True)
    D_ab = free_diffusivity(r_ab)
    print(f"antibody        mw={mw_ab:>8.0f} g/mol  "
          f"r={r_ab:.3f} nm  D={D_ab:.2f} um^2/s")

    # Antibody hydrodynamic radius ~5-6 nm (well-known IgG value).
    assert 4.0 <= r_ab <= 6.5, f"antibody radius off: {r_ab}"

    # HONESTY NOTE ON THE ANTIBODY FREE-D TARGET
    # ------------------------------------------
    # The contract's "~1-10 um^2/s" for an antibody is the HINDERED tissue
    # diffusivity produced downstream by brownian.py, NOT this free (dilute-fluid)
    # Stokes-Einstein value. The physically correct FREE diffusivity of a ~150 kDa
    # IgG at 310 K in 1.2e-3 Pa.s is ~34 um^2/s: literature reports IgG D0 ~= 4.0e-7
    # cm^2/s (~40 um^2/s) in water at 20 C, which rescales (D ~ T/eta) to ~35 um^2/s
    # here. We refuse to shrink the radius to fake a smaller free-D. What the free
    # value must show is a clear, order-of-magnitude-scale drop vs small molecules;
    # the additional drop to single-digit um^2/s comes from ECM crowding (hindrance).
    assert 10.0 <= D_ab <= 60.0, f"antibody free D out of physical range: {D_ab}"

    # --- Size separation: small molecule diffuses far faster than the antibody. ---
    # For FREE Stokes-Einstein, D_small/D_ab == r_ab/r_small (~12x here). A 30x+
    # separation only appears after ECM hindrance hits the large antibody harder.
    ratio = D_small / D_ab
    print(f"D_small / D_antibody = {ratio:.1f}x  (free-D size separation)")
    assert ratio >= 8.0, f"expected small-mol D >> antibody D, got {ratio:.1f}x"

    # --- D scales as 1/r (monotonic decreasing in radius) ---
    radii = [0.4, 0.8, 1.6, 3.2, 6.4]
    Ds_r = [free_diffusivity(r) for r in radii]
    assert all(Ds_r[i] > Ds_r[i + 1] for i in range(len(Ds_r) - 1)), \
        "D must decrease monotonically with radius"
    # verify inverse proportionality: D*r is constant
    products = [D * r for D, r in zip(Ds_r, radii)]
    assert max(products) - min(products) < 1e-6 * max(products), \
        "D*r should be constant (D ~ 1/r)"
    print(f"D~1/r check: D*r products all equal ({products[0]:.4f}) OK")

    # --- D scales as 1/eta (monotonic decreasing in viscosity) ---
    etas = [0.8e-3, 1.0e-3, 1.2e-3, 2.0e-3, 4.0e-3]
    Ds_eta = [free_diffusivity(1.0, viscosity_pa_s=e) for e in etas]
    assert all(Ds_eta[i] > Ds_eta[i + 1] for i in range(len(Ds_eta) - 1)), \
        "D must decrease monotonically with viscosity"
    products_eta = [D * e for D, e in zip(Ds_eta, etas)]
    assert max(products_eta) - min(products_eta) < 1e-6 * max(products_eta), \
        "D*eta should be constant (D ~ 1/eta)"
    print(f"D~1/eta check: D*eta products all equal ({products_eta[0]:.6f}) OK")

    print("ALL CHECKS PASSED")
