"""Treatment / dosing module: drug lookup, systemic pharmacokinetics, and the
spatial delivery source that couples blood concentration into the tissue PDE.

Responsibilities (see ARCHITECTURE_CONTRACT.md):

- ``get_drug``       — fetch a :class:`DrugParams` from ``core.DRUG_LIBRARY``.
- ``plasma_profile`` — the vascular (blood) concentration C_plasma(t) in uM, a
  multi-dose exponential-decay 1-compartment model driven by :class:`PKConfig`.
- ``delivery_source``— the (ny, nx) instantaneous drug delivery-rate source that
  ``physics.diffusion`` adds to the transport PDE. It is a simple linear leak
  proxy: vessels leak drug into tissue in proportion to the local vascular
  density, the drug's blood-brain-barrier penetrance, and the current blood
  concentration.

Units: concentration uM, time seconds. All functions are vectorized and
deterministic (no randomness).
"""
from __future__ import annotations

import numpy as np

from tmesim.core import DRUG_LIBRARY, DrugParams, Microenvironment, PKConfig


def get_drug(name: str) -> DrugParams:
    """Return the :class:`DrugParams` for ``name`` from ``core.DRUG_LIBRARY``.

    Parameters
    ----------
    name : str
        Key into the drug library (e.g. ``"temozolomide-like"``).

    Returns
    -------
    DrugParams
        The parameter set for the requested drug.

    Raises
    ------
    KeyError
        If ``name`` is not in the library. The message lists the available
        names so the caller can correct a typo quickly.
    """
    try:
        return DRUG_LIBRARY[name]
    except KeyError:
        available = ", ".join(sorted(DRUG_LIBRARY)) or "(library is empty)"
        raise KeyError(
            f"Unknown drug {name!r}. Available drugs: {available}."
        ) from None


def plasma_profile(pk: PKConfig, times: np.ndarray) -> np.ndarray:
    """Vascular (blood) drug concentration in uM at each requested time.

    A 1-compartment multi-dose model. Each of ``pk.n_doses`` doses is
    administered at ``k * pk.dose_interval`` seconds (k = 0, 1, ...), instantly
    reaching ``pk.plasma_peak`` uM and then decaying with first-order
    elimination at the drug's plasma ``pk.half_life``. The total plasma
    concentration is the linear superposition of the doses given so far::

        C(t) = plasma_peak * sum_k  H(t - k*dt_dose) * exp(-lambda*(t - k*dt_dose))
        lambda = ln(2) / half_life

    where ``H`` is the Heaviside step (a dose contributes nothing before it is
    given). Superposition means overlapping doses accumulate, so a second dose
    given before the first has washed out raises the peak above ``plasma_peak``.

    Parameters
    ----------
    pk : PKConfig
        Pharmacokinetic settings (peak, half-life, dose count, interval).
    times : np.ndarray
        Times in seconds at which to evaluate the profile. Any shape; the
        result matches it.

    Returns
    -------
    np.ndarray
        Plasma concentration in uM, same shape as ``times``, all >= 0.
    """
    times = np.asarray(times, dtype=np.float64)
    lam = np.log(2.0) / pk.half_life  # elimination rate constant, 1/s

    # Dose administration times: (n_doses,) column so it broadcasts over `times`.
    n_doses = max(int(pk.n_doses), 0)
    if n_doses == 0:
        return np.zeros_like(times)
    dose_times = np.arange(n_doses, dtype=np.float64) * pk.dose_interval

    # elapsed[k, ...] = time since dose k at each requested time.
    elapsed = times[None, ...] - dose_times.reshape((n_doses,) + (1,) * times.ndim)

    # Each dose contributes only after it is given (elapsed >= 0).
    active = elapsed >= 0.0
    contrib = np.where(active, pk.plasma_peak * np.exp(-lam * np.where(active, elapsed, 0.0)), 0.0)

    return contrib.sum(axis=0)


def delivery_source(
    micro: Microenvironment,
    drug: DrugParams,
    plasma_conc: float,
    bbb_disruption: "np.ndarray | float | None" = None,
    efflux: bool = False,
) -> np.ndarray:
    """Instantaneous (ny, nx) drug delivery-rate source into the tissue.

    A simple linear vascular-leak proxy: at each grid point the rate at which
    drug enters the interstitium is proportional to the local perfused vessel
    fraction, the drug's blood-brain-barrier penetrance, and the current blood
    concentration::

        S(x, y) = vascular_density(x, y) * bbb_penetrance * plasma_conc

    (A fuller model would use ``(plasma_conc - tissue_conc)`` as the driving
    gradient; the transport PDE already removes drug via clearance/binding, so
    this linear-leak form keeps the treatment module decoupled from the tissue
    state and is what the contract specifies.)

    Regional barrier integrity (optional, OFF by default)
    -----------------------------------------------------
    ``bbb_penetrance`` is a pure function of the drug's 2D descriptors, so the
    form above predicts a contrast-enhancing-vs-non-enhancing fold change of
    EXACTLY 1.000 for every drug — against measured within-patient fold changes
    of 1.45x (levetiracetam) and 4.26x (cefazolin), p < 0.0001. Passing
    ``bbb_disruption`` switches on the measured-derived additive paracellular
    leak from :mod:`tmesim.reconstruct.barrier`::

        S(x, y) = vascular_density(x, y)
                  * clip(bbb_penetrance + L * bbb_disruption(x, y), 0, 1)
                  * plasma_conc

    with ``L = 0.238`` fitted across three independent drug/assay combinations
    (see :data:`tmesim.reconstruct.barrier.L_PARACELLULAR_LEAK`). An ADDITIVE
    form is used because a multiplicative one is drug-dependent (3.2x spread
    across drugs vs 1.3x for additive) and therefore cannot be a property of the
    barrier.

    **Backward compatibility is exact.** ``bbb_disruption=None`` (the default)
    short-circuits before any barrier code runs and returns the original
    expression bit-for-bit; ``bbb_disruption=0`` likewise returns the original
    values exactly, since ``x + L*0.0 == x`` in IEEE-754 and the clip is a no-op
    for ``bbb_penetrance`` already in [0, 1]. Both are covered by tests.

    Active efflux (optional, OFF by default, and OFF is the RECOMMENDATION)
    ----------------------------------------------------------------------
    ``efflux=True`` multiplies the TRANSCELLULAR penetrance by a literature-
    sourced P-gp/BCRP factor from :mod:`tmesim.molecular.efflux`, looked up by
    ``drug.name``::

        S = vascular_density
            * clip(bbb_penetrance * efflux_factor(drug.name)
                   + L * bbb_disruption, 0, 1)
            * plasma_conc

    Efflux multiplies only the transcellular term because a disrupted tight
    junction is a paracellular shunt AROUND the pump. Drugs that are not in the
    efflux table, or that have no published efflux ratio, get a factor of exactly
    1.0, so ``efflux=True`` is a no-op for most of the library.

    **This defaults to False on evidence, not on caution.** Re-scoring the
    six-drug measured human glioma partition scorecard with the term active makes
    it WORSE (geometric-mean absolute fold error 1.99x -> 2.43x, order-of-
    magnitude band 6/6 -> 5/6), and the term is structurally incapable of moving
    the GIST tumour:normal target at all because a scalar per-drug factor cancels
    out of a within-drug tissue ratio. See ``results/efflux_impact.md`` for the
    full before/after and the reasons. The mechanism is implemented, cited and
    available; it is not enabled because it does not validate.

    Engine usage
    ------------
    The transport engine calls ``plasma_profile`` once to get C_plasma(t) over
    the dosing course, then applies this source each timestep. Because the field
    is linear in ``plasma_conc``, the engine may equivalently build the unit
    field once with ``plasma_conc=1.0`` and multiply it by C_plasma(t) at every
    step instead of rebuilding it — the two are identical. This stays true with
    the barrier term on: it modulates the penetrance, not the plasma coupling.

    Parameters
    ----------
    micro : Microenvironment
        Supplies ``micro.fields.vascular_density`` (ny, nx), the vessel map.
    drug : DrugParams
        Supplies ``bbb_penetrance``.
    plasma_conc : float
        Blood concentration in uM at this instant (from ``plasma_profile``).
    bbb_disruption : ndarray (ny, nx), float, or None
        Per-location barrier disruption in [0, 1] (0 = intact, 1 = fully
        disrupted). ``None`` (default) disables the barrier term entirely and
        reproduces the pre-fix source exactly. Build a field with
        :func:`tmesim.reconstruct.barrier.bbb_disruption_field`.
    efflux : bool
        ``False`` (default) leaves the function bit-for-bit as it was. ``True``
        applies the P-gp/BCRP penetrance factor for ``drug.name`` from
        :mod:`tmesim.molecular.efflux`. Default-off is the validated
        recommendation, not a placeholder — see above.

    Returns
    -------
    np.ndarray
        (ny, nx) non-negative delivery-rate source in uM/s (per unit leak),
        same shape as the vascular density field.
    """
    vascular = np.asarray(micro.fields.vascular_density, dtype=np.float64)
    penetrance = drug.bbb_penetrance
    if efflux:
        # Imported lazily so the default path touches no extra module.
        from tmesim.molecular.efflux import efflux_factor
        penetrance = penetrance * efflux_factor(drug.name)
    if bbb_disruption is None:
        return vascular * penetrance * float(plasma_conc)

    # Imported lazily: reconstruct depends on treatment, not the other way round.
    from tmesim.reconstruct.barrier import effective_penetrance

    p_eff = effective_penetrance(penetrance, bbb_disruption)
    return vascular * p_eff * float(plasma_conc)


if __name__ == "__main__":
    # --- Self-test: exercise the PK model and the delivery source ------------
    import numpy as np

    # np.trapz was renamed to np.trapezoid in newer numpy; support both.
    _trapz = getattr(np, "trapezoid", np.trapz)

    from tmesim.core import (
        CellPopulation,
        CellState,
        CellType,
        Domain,
        Microenvironment,
        TissueFields,
    )

    # 1. get_drug: valid lookup + helpful KeyError -------------------------
    tmz = get_drug("temozolomide-like")
    assert isinstance(tmz, DrugParams) and tmz.name == "temozolomide-like"
    try:
        get_drug("aspirin")
        raise AssertionError("expected KeyError for unknown drug")
    except KeyError as e:
        assert "Available drugs" in str(e) and "temozolomide-like" in str(e)

    # 2. plasma_profile: single-dose decay --------------------------------
    pk1 = PKConfig(plasma_peak=50.0, half_life=3600.0, n_doses=1,
                   dose_interval=3600.0 * 24)
    t = np.linspace(0.0, 3600.0 * 8, 2000)
    c1 = plasma_profile(pk1, t)
    assert c1.shape == t.shape
    assert np.all(c1 >= 0.0)
    assert np.isclose(c1[0], 50.0), c1[0]          # peak at t=0
    assert np.all(np.diff(c1) <= 1e-9)             # monotonically decays
    # half-life check: concentration halves after one half_life
    c_half = plasma_profile(pk1, np.array([3600.0]))[0]
    assert np.isclose(c_half, 25.0, rtol=1e-3), c_half

    # 3. plasma_profile: superposition of multiple doses ------------------
    pk3 = PKConfig(plasma_peak=50.0, half_life=3600.0, n_doses=3,
                   dose_interval=3600.0 * 6)
    c3 = plasma_profile(pk3, t)
    # Just after a dose, concentration must jump up vs. just before it.
    just_before = plasma_profile(pk3, np.array([3600.0 * 6 - 1.0]))[0]
    just_after = plasma_profile(pk3, np.array([3600.0 * 6 + 1.0]))[0]
    assert just_after > just_before, (just_before, just_after)
    # The second dose's contribution stacks on the tail of the first, so the
    # peak of the 3-dose profile exceeds the single-dose peak of 50.
    assert c3.max() > 50.0, c3.max()
    # Superposition sanity: 3-dose profile == sum of three shifted single doses.
    manual = np.zeros_like(t)
    for k in range(3):
        manual += plasma_profile(
            PKConfig(plasma_peak=50.0, half_life=3600.0, n_doses=1),
            t - k * 3600.0 * 6,
        )
    assert np.allclose(c3, manual, atol=1e-9)
    # More doses deliver more total drug (area under curve).
    assert _trapz(c3, t) > _trapz(c1, t)

    # 4. delivery_source: scales with vasculature and bbb_penetrance ------
    domain = Domain(nx=20, ny=16, dx=15.0)
    vasc = np.zeros(domain.shape)
    vasc[:, :10] = 0.2   # left half: sparse vessels
    vasc[:, 10:] = 0.8   # right half: dense vessels
    fields = TissueFields(ecm_density=domain.ones() * 0.5, vascular_density=vasc)
    cells = CellPopulation(
        x=np.array([10.0]), y=np.array([10.0]),
        cell_type=np.array([CellType.TUMOR]), state=np.array([CellState.VIABLE]),
        sensitivity=np.array([0.5]),
    )
    micro = Microenvironment(domain=domain, fields=fields, cells=cells)

    drug_lo = DrugParams(name="lowbbb", bbb_penetrance=0.1)
    drug_hi = DrugParams(name="highbbb", bbb_penetrance=0.5)

    src = delivery_source(micro, drug_hi, plasma_conc=10.0)
    assert src.shape == domain.shape
    assert np.all(src >= 0.0)
    # Higher where vascular_density is higher.
    assert src[:, 10:].mean() > src[:, :10].mean()
    # Scales linearly with bbb_penetrance (0.5 / 0.1 == 5x).
    src_lo = delivery_source(micro, drug_lo, plasma_conc=10.0)
    assert np.allclose(src, src_lo * 5.0), (src.mean(), src_lo.mean())
    # Linear in plasma_conc: unit field * plasma == full field (engine shortcut).
    unit = delivery_source(micro, drug_hi, plasma_conc=1.0)
    assert np.allclose(unit * 10.0, src)
    # Zero blood concentration -> zero source.
    assert np.all(delivery_source(micro, drug_hi, plasma_conc=0.0) == 0.0)

    print("drugs.py self-test PASSED")
    print(f"  single-dose peak={c1[0]:.2f} uM, halves to {c_half:.2f} uM at 1 half-life")
    print(f"  3-dose superposed peak={c3.max():.2f} uM (> single-dose 50)")
    print(f"  AUC 3-dose={_trapz(c3, t):.1f} vs 1-dose={_trapz(c1, t):.1f} uM*s")
    print(f"  delivery_source dense/sparse vessel mean ratio="
          f"{src[:, 10:].mean() / src[:, :10].mean():.2f} (vasc 0.8/0.2)")
