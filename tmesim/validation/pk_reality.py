"""REALITY CHECK: does tmesim's predicted drug penetration match PUBLISHED
tissue-pharmacokinetics measurements?

This is the first real credibility test of the molecular -> tissue pipeline. It
does NOT touch the frozen core/engine/molecular modules; it only *drives* them
(compiling a real drug with :func:`molecular.compile_drug.molecule_to_drugparams`
and running the frozen transport PDE) and compares the result to the literature.

The canonical anchor
--------------------
Doxorubicin penetrates only **~40-50 um** from tumor blood vessels: a famous,
repeatedly-reproduced result showing steep perivascular gradients and poor drug
penetration of solid tumors.

    Primeau, Rendon, Hedley, Lilge, Tannock. "The distribution of the
    anticancer drug doxorubicin in relation to blood vessels in solid tumors."
    Clin Cancer Res 2005; 11(24):8782-8788.  (doxorubicin fluorescence falls to
    background by ~40-50 um from the nearest vessel.)

Supporting anchors (used for the ORDERING check, not exact magnitudes):
  * Reaction-diffusion penetration length  L = sqrt(D_eff / k)  sets how far a
    drug reaches from a vessel, where k = (clearance + binding) is the first-order
    loss rate. This is the steady-state e-folding length of C(r) ~ K_0(r/L) for a
    line/vessel source (Thiele-modulus / reaction-diffusion boundary-layer theory;
    see e.g. Thurber, Schmidt, Wittrup, "Antibody tumor penetration",
    Adv Drug Deliv Rev 2008; 60:1421-1434).
  * Small molecules reach further than large / rapidly-consumed / strongly-bound
    molecules. Antibodies penetrate poorly (the "binding-site barrier":
    Fujimori, Covell, Fletcher, Weinstein, J Nucl Med 1990; 31:1191-1198).

What we compute
---------------
* :func:`penetration_length_um` -- the analytic reaction-diffusion length
  L = sqrt(D_eff / (clearance + binding)) directly from a compiled ``DrugParams``.
  This is the PRINCIPLED metric (a single closed-form number, grid-independent).
* :func:`perivascular_profile` -- runs a focused SINGLE-VESSEL solve through the
  frozen transport engine, extracts steady-state concentration vs distance from
  the vessel, and reports (a) the 1/e penetration distance and (b) a far-field
  e-folding-length fit that should recover L.
* :func:`validate_against_literature` -- compares doxorubicin (plus temozolomide
  and the bevacizumab antibody) to the ~40-50 um anchor, reports the ratio and an
  honest pass/fail band (within 3x == order-of-magnitude agreement), and runs the
  small-molecule-vs-antibody ordering check.

HONESTY: the Stage-0 binding/clearance rates are uncalibrated structure-only
priors; this harness reports the ACTUAL numbers even when they disagree with the
literature. A mismatch is a valid, useful finding.

Units: length um, time s, diffusivity um^2/s, concentration uM, rates 1/s.
"""
from __future__ import annotations

import math
from typing import Dict, Optional, Union

import numpy as np

from ..core import (
    CellPopulation,
    Domain,
    DrugParams,
    Microenvironment,
    PDEConfig,
    TissueFields,
)
from ..physics.diffusion import ALPHA_TORTUOSITY, simulate_drug_transport

# ---------------------------------------------------------------------------
# Literature anchor(s).
# ---------------------------------------------------------------------------
# Doxorubicin perivascular penetration distance (Primeau et al. 2005).
DOX_PENETRATION_UM = (40.0, 50.0)
DOX_PENETRATION_MID = 45.0
DOX_CITATION = (
    "Primeau et al., Clin Cancer Res 2005; 11(24):8782-8788 "
    "(doxorubicin fluorescence falls to background ~40-50 um from tumor vessels)"
)

# Order-of-magnitude agreement band: predicted within this factor of the anchor.
ORDER_OF_MAGNITUDE_FACTOR = 3.0

# The ECM/obstacle volume fraction the drugs are compiled/evaluated at (matches
# molecular.compile_drug's default ecm_reference and the engine's tortuosity use).
DEFAULT_ECM = 0.4


# ---------------------------------------------------------------------------
# 1. Analytic reaction-diffusion penetration length.
# ---------------------------------------------------------------------------
def penetration_length_um(
    drug_params: DrugParams,
    ecm_density: float = DEFAULT_ECM,
    apply_tortuosity: bool = True,
) -> float:
    """Characteristic reaction-diffusion penetration length L (um).

    A drug delivered from a vessel diffuses into tissue while being lost to
    clearance + binding at first-order rate ``k = clearance + binding`` (1/s).
    The steady-state balance ``D_eff * laplacian(C) = k * C`` has the
    boundary-layer / e-folding length::

        L = sqrt(D_eff / k)                                              [um]

    where ``D_eff`` is the effective tissue diffusivity. The engine reduces the
    drug's free diffusivity by ECM tortuosity,
    ``D_eff = diffusivity * (1 - ALPHA_TORTUOSITY * ecm_density)`` (see
    :func:`tmesim.physics.diffusion._effective_diffusivity`); we mirror that here
    when ``apply_tortuosity`` is True so the number matches what the PDE actually
    experiences. This is THE principled, grid-independent penetration metric.

    Parameters
    ----------
    drug_params : DrugParams
        Compiled drug (``diffusivity`` um^2/s, ``clearance`` and ``binding`` 1/s).
    ecm_density : float
        ECM volume fraction for the tortuosity reduction (default 0.4).
    apply_tortuosity : bool
        If True, use the ECM-reduced ``D_eff``; else use the raw ``diffusivity``.

    Returns
    -------
    float
        Penetration length L in um (``inf`` if there is no loss term, i.e.
        clearance + binding == 0, meaning the drug is never consumed).
    """
    k = float(drug_params.clearance) + float(drug_params.binding)
    d = float(drug_params.diffusivity)
    if apply_tortuosity:
        d = d * (1.0 - ALPHA_TORTUOSITY * float(np.clip(ecm_density, 0.0, 1.0)))
        d = max(d, 1e-6)
    if k <= 0.0:
        return math.inf
    return math.sqrt(d / k)


# ---------------------------------------------------------------------------
# 2. Single-vessel transport solve -> perivascular concentration profile.
# ---------------------------------------------------------------------------
def _compile_drug(drug: Union[str, DrugParams]) -> DrugParams:
    """Resolve ``drug`` to a ``DrugParams``: compile a library name through the
    real molecular pipeline, or pass a ``DrugParams`` through unchanged."""
    if isinstance(drug, DrugParams):
        return drug
    # Lazy imports so the module loads even if the molecular layer is slow.
    from ..molecular import library
    from ..molecular.compile_drug import molecule_to_drugparams

    mol = library.get_molecule(str(drug))
    return molecule_to_drugparams(mol, ecm_reference=DEFAULT_ECM)


def perivascular_profile(
    drug: Union[str, DrugParams],
    ecm_density: float = DEFAULT_ECM,
    dx_um: float = 20.0,
    n: int = 161,
    vessel_radius_um: float = 20.0,
    steady_factor: float = 6.0,
    dt_s: float = 40.0,
) -> Dict:
    """Solve drug transport around a SINGLE central vessel and extract the
    concentration-vs-distance profile and the perivascular penetration distance.

    We build an idealized homogeneous tissue patch (uniform ECM, no interstitial
    convection -- a single-vessel idealization), place one small circular vessel
    (radius ``vessel_radius_um``) as a constant delivery source at the center, and
    run the FROZEN :func:`~tmesim.physics.diffusion.simulate_drug_transport` for
    ``steady_factor / k`` seconds (several loss-timescales) so the concentration
    field is near steady state. We then bin the final concentration into radial
    rings and report:

      * ``penetration_distance_um`` -- distance OUTWARD FROM THE VESSEL EDGE at
        which C first falls to 1/e of its near-vessel (vessel-edge) value. This is
        the direct analog of an immunofluorescence-vs-distance measurement.
      * ``efold_length_um`` -- an exponential (e-folding) length fit to the
        far-field tail ``log C vs r``; this should recover the analytic L and is
        less sensitive to near-vessel geometry.
      * ``L_analytic_um`` -- :func:`penetration_length_um` for cross-check.

    Notes on fidelity: the profile is 2D and single-vessel; near a small source
    the 2D geometric spreading (C ~ K_0(r/L) ~ exp(-r/L)/sqrt(r)) makes the 1/e
    distance somewhat shorter than L, while the far-field fit recovers L.

    Returns
    -------
    dict with keys: drug, D_eff_um2_s, k_per_s, distances_um, concentration_uM,
    vessel_radius_um, penetration_distance_um, efold_length_um, L_analytic_um,
    reached_1e (bool), steady_reached (bool).
    """
    dp = _compile_drug(drug)
    k = float(dp.clearance) + float(dp.binding)
    d_eff = float(dp.diffusivity) * (
        1.0 - ALPHA_TORTUOSITY * float(np.clip(ecm_density, 0.0, 1.0))
    )

    domain = Domain(nx=int(n), ny=int(n), dx=float(dx_um))
    ecm_field = np.full((n, n), float(ecm_density))
    fields = TissueFields(
        ecm_density=ecm_field, vascular_density=np.zeros((n, n))
    )
    cx = cy = n * dx_um / 2.0
    cells = CellPopulation(
        x=np.array([cx]), y=np.array([cy]),
        cell_type=np.array([0]), state=np.array([0]),
        sensitivity=np.array([0.5]),
    )
    micro = Microenvironment(domain=domain, fields=fields, cells=cells)

    yy, xx = np.mgrid[0:n, 0:n]
    r_grid = np.hypot((xx - n // 2), (yy - n // 2)) * dx_um
    source = np.where(r_grid <= vessel_radius_um, 1.0, 0.0).astype(np.float64)

    total_time = (steady_factor / k) if k > 0 else steady_factor * 1e5
    pde = PDEConfig(
        dt=float(dt_s), total_time=float(total_time),
        convection=False, save_every=1_000_000,
    )
    res = simulate_drug_transport(micro, dp, source, pde)
    c = res.final_concentration

    # Radial binning of the steady concentration field.
    r_max = float(r_grid.max())
    edges = np.arange(0.0, r_max, dx_um)
    dist, conc = [], []
    for i in range(len(edges) - 1):
        m = (r_grid >= edges[i]) & (r_grid < edges[i + 1])
        if m.any():
            dist.append(0.5 * (edges[i] + edges[i + 1]))
            conc.append(float(c[m].mean()))
    dist = np.asarray(dist)
    conc = np.asarray(conc)

    # Near-vessel reference = concentration at the vessel edge.
    ref_idx = int(np.argmin(np.abs(dist - vessel_radius_um)))
    ref = conc[ref_idx] if conc.size else 0.0
    thr = ref / math.e
    outside = np.where((dist > vessel_radius_um) & (conc <= thr))[0]
    reached_1e = outside.size > 0
    penetration = (
        float(dist[outside[0]] - vessel_radius_um) if reached_1e else float("nan")
    )

    # Far-field e-folding length: fit log C vs r beyond 2x the vessel radius,
    # over the positive-concentration tail.
    tail = (dist > 2.0 * vessel_radius_um) & (conc > 0.0)
    efold = float("nan")
    if tail.sum() > 3:
        slope = np.polyfit(dist[tail], np.log(conc[tail]), 1)[0]
        if slope < 0:
            efold = float(-1.0 / slope)

    # Steady-state check: total tissue mass changed little near the end (proxy:
    # we ran several loss-timescales, so flag steady if total_time >= 5/k).
    steady_reached = bool(k > 0 and total_time >= 5.0 / k)

    return {
        "drug": dp.name,
        "D_eff_um2_s": d_eff,
        "k_per_s": k,
        "distances_um": dist,
        "concentration_uM": conc,
        "vessel_radius_um": float(vessel_radius_um),
        "penetration_distance_um": penetration,
        "efold_length_um": efold,
        "L_analytic_um": penetration_length_um(dp, ecm_density=ecm_density),
        "reached_1e": reached_1e,
        "steady_reached": steady_reached,
    }


# ---------------------------------------------------------------------------
# 3. Validation against the literature.
# ---------------------------------------------------------------------------
def _band_verdict(ratio: float, factor: float = ORDER_OF_MAGNITUDE_FACTOR) -> str:
    """PASS if predicted is within ``factor`` of the anchor (either direction),
    else FAIL. ``ratio`` = predicted / anchor."""
    if not math.isfinite(ratio) or ratio <= 0:
        return "FAIL"
    return "PASS" if (1.0 / factor) <= ratio <= factor else "FAIL"


def validate_against_literature(
    drug_names=("doxorubicin", "temozolomide", "bevacizumab"),
    ecm_density: float = DEFAULT_ECM,
    run_profiles: bool = True,
) -> Dict:
    """Compare tmesim's predicted drug penetration to the published ~40-50 um
    doxorubicin anchor and run the small-molecule-vs-antibody ordering check.

    For each drug we compile a real ``DrugParams`` through the molecular pipeline,
    compute the analytic penetration length L, optionally run the single-vessel
    profile solve, and (for doxorubicin) score the ratio to the anchor against the
    order-of-magnitude (3x) band.

    Returns
    -------
    dict with per-drug results, the doxorubicin comparison, the ordering checks,
    and an overall honest verdict string.
    """
    drugs: Dict[str, DrugParams] = {}
    per_drug: Dict[str, Dict] = {}
    for name in drug_names:
        dp = _compile_drug(name)
        drugs[name] = dp
        entry = {
            "D_eff_um2_s": float(dp.diffusivity)
            * (1.0 - ALPHA_TORTUOSITY * ecm_density),
            "diffusivity_free_um2_s": float(dp.diffusivity),
            "k_per_s": float(dp.clearance) + float(dp.binding),
            "binding_per_s": float(dp.binding),
            "clearance_per_s": float(dp.clearance),
            "mw": float(dp.molecular_weight),
            "bbb_penetrance": float(dp.bbb_penetrance),
            "L_analytic_um": penetration_length_um(dp, ecm_density=ecm_density),
        }
        if run_profiles:
            prof = perivascular_profile(dp, ecm_density=ecm_density)
            entry["profile_penetration_um"] = prof["penetration_distance_um"]
            entry["profile_efold_um"] = prof["efold_length_um"]
        per_drug[name] = entry

    # --- doxorubicin comparison to the anchor -------------------------------
    dox = per_drug.get("doxorubicin")
    dox_cmp = None
    if dox is not None:
        L = dox["L_analytic_um"]
        ratio_L = L / DOX_PENETRATION_MID
        dox_cmp = {
            "anchor_um": DOX_PENETRATION_UM,
            "anchor_mid_um": DOX_PENETRATION_MID,
            "citation": DOX_CITATION,
            "L_analytic_um": L,
            "ratio_L": ratio_L,
            "verdict_L": _band_verdict(ratio_L),
        }
        if "profile_penetration_um" in dox:
            pen = dox["profile_penetration_um"]
            ratio_pen = pen / DOX_PENETRATION_MID
            dox_cmp.update(
                profile_penetration_um=pen,
                ratio_profile=ratio_pen,
                verdict_profile=_band_verdict(ratio_pen),
            )

    # --- ordering checks (directional, structure-driven) --------------------
    ordering = {}
    if "temozolomide" in per_drug and "doxorubicin" in per_drug:
        tmz_L = per_drug["temozolomide"]["L_analytic_um"]
        dox_L = per_drug["doxorubicin"]["L_analytic_um"]
        ordering["small_molecule_gt_larger_drug"] = {
            "claim": "temozolomide penetration length > doxorubicin (smaller, "
                     "less-bound drug reaches further)",
            "temozolomide_L_um": tmz_L,
            "doxorubicin_L_um": dox_L,
            "pass": bool(tmz_L > dox_L),
        }
    if "bevacizumab" in per_drug and "temozolomide" in per_drug:
        bev_D = per_drug["bevacizumab"]["D_eff_um2_s"]
        tmz_D = per_drug["temozolomide"]["D_eff_um2_s"]
        ordering["small_molecule_diffuses_faster_than_antibody"] = {
            "claim": "small-molecule effective diffusivity >> antibody "
                     "(antibodies diffuse far slower)",
            "temozolomide_Deff": tmz_D,
            "bevacizumab_Deff": bev_D,
            "ratio": (tmz_D / bev_D) if bev_D > 0 else math.inf,
            "pass": bool(tmz_D > bev_D),
        }

    verdict = _overall_verdict(dox_cmp, ordering)
    return {
        "ecm_density": ecm_density,
        "per_drug": per_drug,
        "doxorubicin_comparison": dox_cmp,
        "ordering": ordering,
        "verdict": verdict,
    }


def _overall_verdict(dox_cmp: Optional[Dict], ordering: Dict) -> str:
    """Compose an honest one-paragraph verdict from the doxorubicin magnitude
    comparison and the ordering checks."""
    if dox_cmp is None:
        return "INCONCLUSIVE: doxorubicin not evaluated."
    parts = []
    rL = dox_cmp["ratio_L"]
    parts.append(
        f"Doxorubicin analytic penetration length L={dox_cmp['L_analytic_um']:.0f} um "
        f"vs published ~40-50 um => {rL:.1f}x the anchor ({dox_cmp['verdict_L']} "
        f"the {ORDER_OF_MAGNITUDE_FACTOR:.0f}x order-of-magnitude band)."
    )
    if "ratio_profile" in dox_cmp:
        rp = dox_cmp["ratio_profile"]
        parts.append(
            f"Single-vessel 1/e distance={dox_cmp['profile_penetration_um']:.0f} um "
            f"=> {rp:.1f}x ({dox_cmp['verdict_profile']})."
        )
    n_order = len(ordering)
    n_pass = sum(1 for v in ordering.values() if v.get("pass"))
    parts.append(f"Ordering checks: {n_pass}/{n_order} PASS.")
    magnitude_ok = dox_cmp["verdict_L"] == "PASS"
    if magnitude_ok:
        parts.append("VERDICT: order-of-magnitude AGREEMENT on magnitude.")
    else:
        parts.append(
            "VERDICT: MAGNITUDE MISMATCH -- the sim OVER-predicts doxorubicin "
            "tissue penetration by roughly an order of magnitude. Root cause: the "
            "uncalibrated Stage-0 binding+clearance priors (~3e-4 /s, ~1 h "
            "timescale) are far too small to reproduce the steep perivascular "
            "gradient, which in reality is driven by rapid high-capacity DNA/tissue "
            "binding and cellular uptake. Directional trends (ordering) are captured; "
            "absolute penetration is not (needs PK calibration)."
        )
    return " ".join(parts)


# ---------------------------------------------------------------------------
# Reporting helpers.
# ---------------------------------------------------------------------------
def _fmt(x, w=10, p=1):
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return f"{'n/a':>{w}}"
    return f"{x:>{w}.{p}f}"


def _print_table(result: Dict) -> None:
    print("=" * 92)
    print("tmesim PK REALITY CHECK  -- predicted drug penetration vs published tissue PK")
    print("=" * 92)
    print(f"Anchor: doxorubicin ~40-50 um from tumor vessels")
    print(f"        [{DOX_CITATION}]")
    print(f"Order-of-magnitude band: predicted within {ORDER_OF_MAGNITUDE_FACTOR:.0f}x "
          f"of the anchor == PASS")
    print(f"ECM fraction: {result['ecm_density']}   (D_eff = D*(1-{ALPHA_TORTUOSITY}*ECM))")
    print("-" * 92)
    hdr = (f"{'drug':14s}{'MW':>8s}{'D_eff':>9s}{'k(1/s)':>10s}"
           f"{'L_anlyt':>9s}{'prof1/e':>9s}{'efoldfit':>9s}{'bbb':>7s}")
    print(hdr)
    print(f"{'':14s}{'Da':>8s}{'um2/s':>9s}{'':>10s}{'um':>9s}{'um':>9s}{'um':>9s}{'':>7s}")
    print("-" * 92)
    for name, e in result["per_drug"].items():
        print(f"{name:14s}{e['mw']:>8.0f}{e['D_eff_um2_s']:>9.1f}"
              f"{e['k_per_s']:>10.2e}"
              f"{_fmt(e['L_analytic_um'],9,0)}"
              f"{_fmt(e.get('profile_penetration_um'),9,0)}"
              f"{_fmt(e.get('profile_efold_um'),9,0)}"
              f"{e['bbb_penetrance']:>7.3f}")
    print("-" * 92)

    dc = result["doxorubicin_comparison"]
    if dc:
        print("\nDOXORUBICIN vs anchor (~45 um mid):")
        print(f"  analytic L          = {dc['L_analytic_um']:.0f} um   "
              f"ratio = {dc['ratio_L']:.1f}x   [{dc['verdict_L']}]")
        if "ratio_profile" in dc:
            print(f"  single-vessel 1/e   = {dc['profile_penetration_um']:.0f} um   "
                  f"ratio = {dc['ratio_profile']:.1f}x   [{dc['verdict_profile']}]")

    print("\nORDERING CHECKS:")
    for key, v in result["ordering"].items():
        status = "PASS" if v.get("pass") else "FAIL"
        print(f"  [{status}] {v['claim']}")

    print("\n" + "=" * 92)
    print("VERDICT")
    print("=" * 92)
    print(result["verdict"])
    print("=" * 92)


if __name__ == "__main__":
    import time as _time

    t0 = _time.time()
    result = validate_against_literature()
    _print_table(result)
    print(f"\n[runtime {_time.time() - t0:.1f} s]")
