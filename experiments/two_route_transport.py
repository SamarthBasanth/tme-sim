"""TWO-ROUTE TISSUE DIFFUSIVITY -- paracellular + transcellular -- scored ONCE,
BLIND, against the 71-compound Pruijn 2008 dataset.

What this script does, in order
-------------------------------
1. ``--stage prereg``  Prints every parameter of the model and its provenance,
   and ASSERTS that each matches ``results/sweep/tworoute_prereg.json``, which
   was written before the model was scored. If anyone tunes a constant to the
   measurement afterwards, this stage fails.
2. ``--stage limits``  Checks the algebra rather than asserting it in prose: the
   ``P -> 0`` limit reproduces ``ecm_mesh.cell_exclusion`` to machine precision,
   the ``phi -> 0`` and zero-contrast limits give 1, and the Maxwell-with-series-
   resistance form is identical to Hasselman & Johnson (1987) Eqn 8.
3. ``--stage pruijn``  Scores the model against all 71 measured tissue diffusion
   coefficients: overall spread, correlation, the matched-molecular-weight test,
   and tirapazamine. Also runs the DISCLOSURE sweep over the one parameter with
   no primary source, so the reader can see the adopted value is not the one that
   maximises agreement.
4. ``--stage drugs``   What the model does to the three drugs the project cares
   about -- temozolomide, cisplatin, doxorubicin -- plus tirapazamine.
5. ``--stage blast``   Blast radius on the published dx = 12 um breast grid, both
   alpha arms, re-using ``experiments/dtissue_sensitivity.py`` unmodified.

THE RULE THIS FILE WAS WRITTEN UNDER
------------------------------------
**Nothing is fitted to Pruijn.** The Pruijn dataset is the cleanest out-of-sample
transport anchor this project has, and fitting any constant to it would destroy
that permanently. The model's structure comes from the effective-medium
literature and its parameters from the repo's own section, from Kao 1993, and
from constants already published in ``permeability.py``. The paper's own Eqn-1
QSAR is deliberately NOT adopted: it is fitted to HT29 and the paper itself shows
its asymptotes moving up to 100-fold across cell lines, so adopting it would make
the engine in-sample against the only clean transport test it has.

The 71 measured values and their provenance live in
``experiments/diffusivity_falsification.py`` and are imported from there, not
re-typed.

Run:
    ./.venv/bin/python experiments/two_route_transport.py --stage prereg
    ./.venv/bin/python experiments/two_route_transport.py --stage limits
    ./.venv/bin/python experiments/two_route_transport.py --stage pruijn
    ./.venv/bin/python experiments/two_route_transport.py --stage drugs
    ./.venv/bin/python experiments/two_route_transport.py --stage blast
"""
from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import sys
import time
from typing import Dict, List, Optional

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from tmesim.molecular import compile_drug, descriptors, ecm_mesh  # noqa: E402
from tmesim.molecular import stokes_einstein, two_route  # noqa: E402
from tmesim.molecular.breast_library import get_molecule_breast  # noqa: E402
from tmesim.molecular.library import get_molecule  # noqa: E402
from tmesim.molecular.types import Molecule  # noqa: E402
from tmesim.validation.pk_reality import DEFAULT_ECM  # noqa: E402

from experiments.diffusivity_falsification import (  # noqa: E402
    PRUIJN_2008_TABLE_S1, TABLE_UNITS_TO_UM2_S, _pearson, _spearman, _ols, fold,
)

OUT_DIR = os.path.join(_ROOT, "results", "sweep")
PREREG = os.path.join(OUT_DIR, "tworoute_prereg.json")

#: Tirapazamine, exactly as experiments/anchor_prototype.py and
#: experiments/cellular_exclusion.py compile it, so "before" columns match.
TPZ_SMILES = "C1=CC=C2C(=C1)[N+](=NC(=[N+]2[O-])N)[O-]"
TPZ_MEASURED_HT29_UM2_S = 40.4   # Pruijn 2008 Table S1, "SN24280 (TPZ)"


def _hr(c: str = "=") -> None:
    print(c * 78)


# =========================================================================== #
# STAGE 1 -- the pre-registration
# =========================================================================== #
def stage_prereg() -> dict:
    _hr()
    print("STAGE 1 -- PRE-REGISTERED PARAMETERS (frozen before any scoring)")
    _hr()
    with open(PREREG) as fh:
        reg = json.load(fh)
    p = reg["parameters"]

    live = {
        "cell_volume_fraction": two_route.DEFAULT_CELL_VOLUME_FRACTION,
        "cell_radius_um": two_route.DEFAULT_CELL_RADIUS_UM,
        "cytoplasm_D_over_water": two_route.CYTOPLASM_D_OVER_WATER,
        "membrane_viscosity_pa_s": two_route.ETA_MEMBRANE_PA_S,
        "bilayer_thickness_nm": two_route.BILAYER_THICKNESS_NM,
        "logk_membrane_max": two_route.LOGK_MEMBRANE_MAX,
        "hbd_logp_decrement": two_route.HBD_LOGP_DECREMENT,
    }
    print(f"  {'parameter':28s} {'live value':>14s} {'pre-registered':>16s}  match")
    for k, v in live.items():
        reg_v = p[k]["value"]
        ok = abs(float(v) - float(reg_v)) <= 1e-12 * max(1.0, abs(float(reg_v)))
        print(f"  {k:28s} {v:>14.6g} {reg_v:>16.6g}  {'OK' if ok else 'DRIFTED'}")
        assert ok, (
            f"{k} has DRIFTED from the pre-registration: live {v} vs {reg_v}. "
            "A parameter changed after the model was scored is a fit to the "
            "measurement, whatever it is called.")
    # permeability.py's own constants must still be the ones adopted.
    from tmesim.molecular import permeability as _perm
    assert two_route.ETA_MEMBRANE_PA_S == _perm._ETA_MEMBRANE_PA_S
    assert two_route.BILAYER_THICKNESS_NM == _perm._BILAYER_THICKNESS_NM
    assert two_route.LOGK_MEMBRANE_MAX == _perm._LOGP_KP_MAX
    assert two_route.ETA_WATER_PA_S == _perm._ETA_WATER_PA_S
    print("\n  membrane constants still identical to permeability.py's published"
          " values: OK")
    print(f"\n  The ONE parameter with no primary source for its magnitude:")
    print(f"    hbd_logp_decrement = {two_route.HBD_LOGP_DECREMENT} log units/donor"
          f" -- {p['hbd_logp_decrement']['status']}")
    print("    Its structure is forced by Walter & Gutknecht 1986 (PMID 3735402)")
    print("    and Conradi 1991/1992 (PMID 1808606 / 1614980). Stage 3 sweeps it")
    print("    and discloses the value that would maximise the Pruijn score.")
    return {"prereg": reg, "live": live}


# =========================================================================== #
# STAGE 2 -- the limits, checked in code
# =========================================================================== #
def stage_limits() -> dict:
    _hr()
    print("STAGE 2 -- LIMIT CHECKS (algebra, not prose)")
    _hr()
    phi = two_route.DEFAULT_CELL_VOLUME_FRACTION
    sigma_m, d_cyto = 265.0, 200.0

    shipped = ecm_mesh.cell_exclusion(phi)["hindrance"]
    zero_p = two_route.two_route_factor(phi, sigma_m, d_cyto, 0.0)["hindrance"]
    print(f"  [A] P -> 0  vs  ecm_mesh.cell_exclusion({phi})['hindrance']")
    print(f"      two_route : {zero_p!r}")
    print(f"      shipped   : {shipped!r}")
    print(f"      |diff|    : {abs(zero_p - shipped):.3g}")
    assert abs(zero_p - shipped) < 1e-15
    print("      -> the two-route law is a STRICT GENERALISATION of the shipped one.")

    z = two_route.two_route_factor(0.0, sigma_m, d_cyto, 1.0)["hindrance"]
    n = two_route.two_route_factor(phi, sigma_m, sigma_m, float("inf"))["hindrance"]
    print(f"\n  [B] phi -> 0 : {z!r}   (must be 1)")
    print(f"  [C] zero contrast (D_cell == D_interstitium, P -> inf) : {n!r}")
    assert abs(z - 1.0) < 1e-15 and abs(n - 1.0) < 1e-15

    print("\n  [D] identity with Hasselman & Johnson (1987) J Composite Materials")
    print("      21:508-515, doi:10.1177/002199838702100602, Eqn 8:")
    worst = 0.0
    for p_ in (1e-4, 1e-2, 1.0, 100.0, 1e4):
        a = two_route.two_route_factor(phi, sigma_m, d_cyto, p_)["hindrance"]
        b = two_route.hasselman_johnson_hindrance(
            phi, sigma_m, d_cyto, two_route.DEFAULT_CELL_RADIUS_UM, p_)
        worst = max(worst, abs(a - b))
        print(f"      P={p_:>9.4g} um/s   ours {a:.14f}   HJ87 {b:.14f}")
    assert worst < 1e-14
    print(f"      max |difference| over the range = {worst:.3g}")

    ceil_ = two_route.two_route_ceiling(phi, sigma_m, d_cyto)
    print(f"\n  [E] REGISTERED CEILING at P -> inf: hindrance {ceil_['hindrance']:.4f}"
          f", i.e. at most {ceil_['transcellular_gain']:.2f}x the impermeable value.")
    print("      Set by the cytoplasm/interstitium conductivity contrast, which is")
    print("      measured (Kao 1993), not tunable. The model CANNOT exceed this.")
    return {"impermeable": shipped, "hj87_max_diff": worst,
            "ceiling": ceil_["hindrance"], "ceiling_gain": ceil_["transcellular_gain"]}


# =========================================================================== #
# STAGE 3 -- Pruijn, scored once, blind
# =========================================================================== #
def _engine_interstitial_D(mw: float) -> float:
    """Today's engine D_eff for a molecular weight -- the UNMODIFIED pipeline.

    This is the matrix (interstitial) conductivity the two-route law consumes:
    free Stokes-Einstein x Brownian micro-tortuosity x Ogston-Renkin mesh gate,
    with no cells anywhere. Compiled through `compile_drug.build_profile` exactly
    as `experiments/diffusivity_falsification.py` does, so the "before" column
    here is bit-comparable with the published falsification.
    """
    return compile_drug.build_profile(
        Molecule(name=f"mw{mw}", smiles="", descriptors={"mw": float(mw)})
    ).effective_diffusivity


def _score(rows: List[dict], key: str) -> dict:
    meas = [r["meas"] for r in rows]
    pred = [r[key] for r in rows]
    lm = [math.log10(v) for v in meas]
    lp = [math.log10(v) for v in pred]
    r = _pearson(lp, lm)
    ratios = sorted(p / m for p, m in zip(pred, meas))
    n = len(ratios)
    med = ratios[n // 2] if n % 2 else (ratios[n // 2 - 1] + ratios[n // 2]) / 2
    return {
        "n": n, "min": min(pred), "max": max(pred), "fold": fold(pred),
        "pearson_log10": r, "r2": r * r, "spearman": _spearman(pred, meas),
        "ratio_min": min(ratios), "ratio_median": med, "ratio_max": max(ratios),
        "within_3x": sum(1 for x in ratios if 1 / 3 <= x <= 3) / n,
        "geo_bias": 10 ** (sum(math.log10(x) for x in ratios) / n),
    }


def _build_rows(hbd_decrement: float, cache: Optional[Dict[float, float]] = None,
                phi: Optional[float] = None) -> List[dict]:
    """One row per Pruijn compound. logP7.4 and HD are the paper's OWN descriptor
    columns -- inputs, never the scored quantity. Compounds with no logP7.4 in
    either the measured or the calculated column are carried with the engine
    value only and excluded from the two-route scoring (n is reported)."""
    cache = cache if cache is not None else {}
    phi = two_route.DEFAULT_CELL_VOLUME_FRACTION if phi is None else phi
    rows = []
    for (name, mw, hd, ha, logp, src, d_meas, d_sem, d_paper) in PRUIJN_2008_TABLE_S1:
        if mw not in cache:
            cache[mw] = _engine_interstitial_D(mw)
        eng = cache[mw]
        radius = stokes_einstein.hydrodynamic_radius_nm(mw)
        row = {
            "name": name, "mw": mw, "hd": hd, "logp": logp, "src": src,
            "meas": d_meas * TABLE_UNITS_TO_UM2_S,
            "paper": d_paper * TABLE_UNITS_TO_UM2_S,
            "eng": eng, "radius_nm": radius,
        }
        if logp is not None:
            tr = two_route.two_route_diffusivity(
                eng, radius, logp, hd, cell_volume_fraction=phi,
                hbd_decrement=hbd_decrement)
            row.update({
                "two_route": tr["effective_D_um2s"],
                "para_only": tr["paracellular_D_um2s"],
                "P_um_s": tr["membrane_permeability_um_s"],
                "gain": tr["transcellular_gain"],
                "log_k": tr["log_k_membrane"],
            })
        rows.append(row)
    return rows


def stage_pruijn(verbose: bool = True) -> dict:
    _hr()
    print("STAGE 3 -- SCORED AGAINST PRUIJN 2008 (n=71, HT29 MCLs). BLIND, ONCE.")
    _hr()
    cache: Dict[float, float] = {}
    t0 = time.time()
    all_rows = _build_rows(two_route.HBD_LOGP_DECREMENT, cache)
    print(f"  compiled {len(cache)} distinct molecular weights through the "
          f"unmodified pipeline ({time.time()-t0:.1f}s)")
    rows = [r for r in all_rows if "two_route" in r]
    print(f"  {len(rows)}/{len(all_rows)} compounds have a logP7.4 in Table S1; "
          f"the other {len(all_rows)-len(rows)} are excluded (no logP -> no "
          "membrane leg). The engine-only columns still use all 71.")

    if verbose:
        print("\n  PER-COMPOUND (um^2/s), sorted by measured:")
        print(f"  {'compound':<20}{'MW':>7}{'HD':>4}{'logP':>7}{'logK':>7}"
              f"{'P(um/s)':>10}{'meas':>9}{'engine':>9}{'para':>8}{'2-route':>9}"
              f"{'2r/meas':>9}")
        for r in sorted(rows, key=lambda r: r["meas"]):
            print(f"  {r['name']:<20}{r['mw']:>7.0f}{r['hd']:>4}{r['logp']:>7.2f}"
                  f"{r['log_k']:>7.2f}{r['P_um_s']:>10.3g}{r['meas']:>9.1f}"
                  f"{r['eng']:>9.1f}{r['para_only']:>8.1f}{r['two_route']:>9.1f}"
                  f"{r['two_route']/r['meas']:>9.2f}")

    out = {}
    for key, label in (("eng", "engine (shipping, MW-only)"),
                       ("para_only", "paracellular only (impermeable cells)"),
                       ("two_route", "TWO-ROUTE (para + transcellular)")):
        out[key] = _score(rows, key)
    meas = [r["meas"] for r in rows]
    out["measured"] = {"min": min(meas), "max": max(meas), "fold": fold(meas)}
    paper = [r["paper"] for r in rows]

    print("\n  RESULT 1 -- SPREAD (the statistic the falsification turned on)")
    print(f"  {'model':40s} {'min':>8s} {'max':>8s} {'fold':>8s}")
    print(f"  {'MEASURED Dmcl':40s} {min(meas):8.1f} {max(meas):8.1f} "
          f"{fold(meas):7.1f}x")
    for key, label in (("eng", "engine (shipping, MW-only)"),
                       ("para_only", "paracellular only (impermeable cells)"),
                       ("two_route", "TWO-ROUTE (para + transcellular)")):
        s = out[key]
        print(f"  {label:40s} {s['min']:8.1f} {s['max']:8.1f} {s['fold']:7.1f}x")
    print(f"  {'(the paper OWN Eqn 1, not adopted)':40s} {min(paper):8.1f} "
          f"{max(paper):8.1f} {fold(paper):7.1f}x")

    print("\n  RESULT 2 -- CORRELATION WITH MEASUREMENT (log10)")
    print(f"  {'model':40s} {'Pearson r':>10s} {'R^2':>8s} {'Spearman':>10s}")
    for key, label in (("eng", "engine (shipping, MW-only)"),
                       ("para_only", "paracellular only"),
                       ("two_route", "TWO-ROUTE")):
        s = out[key]
        print(f"  {label:40s} {s['pearson_log10']:>+10.3f} {s['r2']:>8.3f} "
              f"{s['spearman']:>+10.3f}")

    print("\n  RESULT 3 -- ABSOLUTE ACCURACY")
    print(f"  {'model':40s} {'min':>8s} {'median':>8s} {'max':>8s} "
          f"{'<3x':>6s} {'geo bias':>9s}")
    for key, label in (("eng", "engine (shipping, MW-only)"),
                       ("para_only", "paracellular only"),
                       ("two_route", "TWO-ROUTE")):
        s = out[key]
        print(f"  {label:40s} {s['ratio_min']:8.2f} {s['ratio_median']:8.2f} "
              f"{s['ratio_max']:8.2f} {s['within_3x']:>5.0%} {s['geo_bias']:8.2f}x")

    # ---------------- the kill shot: matched molecular weight ---------------- #
    print("\n  RESULT 4 -- THE MATCHED-MW TEST. Can the model produce ANY spread")
    print("  where a MW-only model is flat by construction?")
    by_mw: Dict[float, List[dict]] = {}
    for r in rows:
        by_mw.setdefault(r["mw"], []).append(r)
    print(f"\n  {'MW':>6}{'n':>4}{'measured range':>20}{'meas fold':>11}"
          f"{'engine':>9}{'para':>8}{'2-route':>10}")
    ties = []
    for mw, sub in sorted(by_mw.items()):
        if len(sub) < 3:
            continue
        m_ = [r["meas"] for r in sub]
        rec = {"mw": mw, "n": len(sub), "meas_fold": fold(m_),
               "meas_lo": min(m_), "meas_hi": max(m_),
               "eng_fold": fold([r["eng"] for r in sub]),
               "para_fold": fold([r["para_only"] for r in sub]),
               "two_route_fold": fold([r["two_route"] for r in sub]),
               "names": [r["name"] for r in sub]}
        ties.append(rec)
        print(f"  {mw:>6.0f}{len(sub):>4}"
              f"{f'{min(m_):.1f} - {max(m_):.1f}':>20}{rec['meas_fold']:>10.1f}x"
              f"{rec['eng_fold']:>8.2f}x{rec['para_fold']:>7.2f}x"
              f"{rec['two_route_fold']:>9.2f}x")
    out["matched_mw"] = ties

    # within-tie rank agreement, the sharpest form of the test
    conc = disc = 0
    for mw, sub in by_mw.items():
        if len(sub) < 2:
            continue
        for i in range(len(sub)):
            for j in range(i + 1, len(sub)):
                a, b = sub[i], sub[j]
                if a["meas"] == b["meas"] or a["two_route"] == b["two_route"]:
                    continue
                if (a["meas"] - b["meas"]) * (a["two_route"] - b["two_route"]) > 0:
                    conc += 1
                else:
                    disc += 1
    out["matched_mw_pairs"] = {"concordant": conc, "discordant": disc,
                               "total": conc + disc}
    print(f"\n  Within-tie pair ordering (all exact-MW pairs, n={conc+disc}): "
          f"{conc} concordant, {disc} discordant "
          f"({conc/(conc+disc):.0%} correct). A MW-only model scores 0/0 -- it "
          "cannot order them at all.")

    # ---------------- tirapazamine ------------------------------------------ #
    tpz = next(r for r in rows if r["name"].startswith("SN24280"))
    out["tirapazamine"] = tpz
    print(f"\n  RESULT 5 -- TIRAPAZAMINE (the project's registered anchor)")
    print(f"    measured, HT29 MCL (Pruijn Table S1) : {tpz['meas']:.1f} um^2/s")
    print(f"    engine, shipping                     : {tpz['eng']:.1f} um^2/s "
          f"({tpz['eng']/tpz['meas']:.1f}x over)")
    print(f"    paracellular only (impermeable cells): {tpz['para_only']:.1f} um^2/s "
          f"({tpz['para_only']/tpz['meas']:.2f}x)")
    print(f"    TWO-ROUTE                            : {tpz['two_route']:.1f} um^2/s "
          f"({tpz['two_route']/tpz['meas']:.2f}x)")
    print(f"    membrane P = {tpz['P_um_s']:.4g} um/s, transcellular gain "
          f"{tpz['gain']:.3f}x")

    # ---------------- disclosure sweep -------------------------------------- #
    print("\n  RESULT 6 -- DISCLOSURE SWEEP over the ONE unsourced parameter.")
    print("  If the adopted value were tuned to this dataset it would be the one")
    print("  that maximises R^2. It is printed alongside so you can check.")
    print(f"\n  {'hbd_decrement':>14}{'spread':>10}{'R^2':>9}{'Spearman':>10}"
          f"{'geo bias':>10}{'<3x':>7}")
    sweep = []
    for dec in (0.0, 0.5, 1.0, 1.5, 2.0, 3.0):
        rr = [r for r in _build_rows(dec, cache) if "two_route" in r]
        s = _score(rr, "two_route")
        s["hbd_decrement"] = dec
        sweep.append(s)
        mark = "  <-- ADOPTED" if abs(dec - two_route.HBD_LOGP_DECREMENT) < 1e-12 else ""
        print(f"  {dec:>14.1f}{s['fold']:>9.1f}x{s['r2']:>9.3f}"
              f"{s['spearman']:>+10.3f}{s['geo_bias']:>9.2f}x{s['within_3x']:>6.0%}"
              f"{mark}")
    best = max(sweep, key=lambda s: s["r2"])
    out["hbd_sweep"] = sweep
    out["hbd_sweep_best_r2"] = best["hbd_decrement"]
    print(f"\n  Best R^2 in the sweep is at hbd_decrement = "
          f"{best['hbd_decrement']:.1f} (R^2 = {best['r2']:.3f}). "
          f"ADOPTED = {two_route.HBD_LOGP_DECREMENT:.1f}.")

    # ---------------- phi sensitivity --------------------------------------- #
    print("\n  RESULT 7 -- sensitivity to phi_cell, the widest-uncertainty input")
    print(f"  {'phi_cell':>9}{'spread':>10}{'R^2':>9}{'geo bias':>10}")
    phis = []
    for ph in (0.40, 0.50, 0.60, 0.72, 0.7635, 0.7892, 0.85):
        rr = [r for r in _build_rows(two_route.HBD_LOGP_DECREMENT, cache, phi=ph)
              if "two_route" in r]
        s = _score(rr, "two_route")
        s["phi"] = ph
        phis.append(s)
        print(f"  {ph:>9.4f}{s['fold']:>9.1f}x{s['r2']:>9.3f}{s['geo_bias']:>9.2f}x")
    out["phi_sweep"] = phis
    out["rows"] = rows
    return out


# =========================================================================== #
# STAGE 4 -- the drugs the project actually cares about
# =========================================================================== #
def stage_drugs() -> dict:
    _hr()
    print("STAGE 4 -- THE DRUGS THE PROJECT CARES ABOUT")
    _hr()
    names = ["temozolomide", "cisplatin", "doxorubicin"]
    mols = {n: get_molecule(n) for n in names}
    mols["tirapazamine"] = descriptors.from_smiles("tirapazamine", TPZ_SMILES)

    print(f"  {'drug':<14}{'MW':>8}{'logP':>7}{'HBD':>5}{'logK':>7}{'P(um/s)':>11}"
          f"{'engine D':>10}{'para':>9}{'2-route':>10}{'gain':>7}{'vs eng':>8}")
    out = []
    for name, mol in mols.items():
        prof = compile_drug.build_profile(mol, ecm_reference=DEFAULT_ECM)
        eng = float(prof.effective_diffusivity)
        radius = float(mol.get("radius_hydrodynamic_nm")
                       or stokes_einstein.hydrodynamic_radius_nm(
                           float(mol.get("mw")), bool(mol.get("is_antibody"))))
        tr = two_route.two_route_diffusivity(
            eng, radius, float(mol.get("logp", 0.0)), float(mol.get("hbd", 0)))
        rec = {"drug": name, "mw": float(mol.get("mw")),
               "logp": float(mol.get("logp", 0.0)), "hbd": float(mol.get("hbd", 0)),
               "log_k": tr["log_k_membrane"], "P_um_s": tr["membrane_permeability_um_s"],
               "engine_D_um2s": eng, "para_D_um2s": tr["paracellular_D_um2s"],
               "two_route_D_um2s": tr["effective_D_um2s"],
               "gain": tr["transcellular_gain"],
               "scale_vs_engine": tr["hindrance"]}
        out.append(rec)
        print(f"  {name:<14}{rec['mw']:>8.2f}{rec['logp']:>7.2f}{rec['hbd']:>5.0f}"
              f"{rec['log_k']:>7.2f}{rec['P_um_s']:>11.4g}{eng:>10.1f}"
              f"{rec['para_D_um2s']:>9.1f}{rec['two_route_D_um2s']:>10.1f}"
              f"{rec['gain']:>7.2f}{rec['scale_vs_engine']:>8.4f}")

    tpz = next(r for r in out if r["drug"] == "tirapazamine")
    print(f"\n  tirapazamine vs its MEASURED HT29 value "
          f"({TPZ_MEASURED_HT29_UM2_S:.1f} um^2/s):")
    print(f"    engine {tpz['engine_D_um2s']:.1f} "
          f"({tpz['engine_D_um2s']/TPZ_MEASURED_HT29_UM2_S:.1f}x over) -> "
          f"two-route {tpz['two_route_D_um2s']:.1f} "
          f"({tpz['two_route_D_um2s']/TPZ_MEASURED_HT29_UM2_S:.2f}x)")

    print("\n  NOTE, and it is not a small one: the engine's `logp` is RDKit's")
    print("  Crippen estimate of the NEUTRAL octanol/water logP, while Pruijn's")
    print("  column is a MEASURED logP7.4 (a distribution coefficient at pH 7.4).")
    print("  For cisplatin in particular Crippen gives +0.20 against an")
    print("  experimental ~-2.2 (permeability.py records this), so the membrane leg")
    print("  inherits that error wholesale. This model does not fix descriptors.")
    return {"drugs": out}


# =========================================================================== #
# STAGE 5 -- blast radius
# =========================================================================== #
def stage_blast(dx_um: float = 12.0, per_drug: bool = True) -> dict:
    _hr()
    print("STAGE 5 -- BLAST RADIUS, published dx = 12 um breast grid, BOTH arms")
    _hr()
    import experiments.dtissue_sensitivity as dts

    ref = dts.reference_transport_constants()
    d_ship = float(ref["D_tissue_cm2_s"])

    # per-drug scales: unlike the cell-exclusion term, the two-route factor is
    # COMPOUND-SPECIFIC, so it is not one global multiplier.
    scales = {}
    for name in dts.TESTABLE:
        # The scored panel straddles two libraries (the breast panel plus the
        # general one), exactly as `compile_breast_drug` does.
        try:
            mol = get_molecule_breast(name)
        except KeyError:
            mol = get_molecule(name)
        prof = compile_drug.build_profile(mol, ecm_reference=DEFAULT_ECM)
        radius = float(mol.get("radius_hydrodynamic_nm")
                       or stokes_einstein.hydrodynamic_radius_nm(
                           float(mol.get("mw")), bool(mol.get("is_antibody"))))
        tr = two_route.two_route_diffusivity(
            float(prof.effective_diffusivity), radius,
            float(mol.get("logp", 0.0)), float(mol.get("hbd", 0)))
        scales[name] = float(tr["hindrance"])

    s_dox = scales["doxorubicin"]
    print(f"  shipping D_tissue (doxorubicin) = {d_ship:.4e} cm^2/s")
    print(f"\n  {'drug':16s} {'two-route scale':>16s} {'fold':>8s} "
          f"{'vs doxorubicin':>15s}")
    for k, v in scales.items():
        print(f"  {k:16s} {v:16.4f} {1/v:7.2f}x {v/s_dox:14.4f}x")
    spread = max(scales.values()) / min(scales.values())
    print(f"\n  spread of the correction ACROSS the scored panel: {spread:.3f}x")
    print("  alpha is fitted on doxorubicin, so in `recal` the doxorubicin scale")
    print("  is absorbed exactly and what survives is the RATIO column above.")

    d_new = d_ship * s_dox
    print(f"\n  doxorubicin-anchored corrected D_tissue = {d_new:.4e} cm^2/s")

    os.makedirs(OUT_DIR, exist_ok=True)
    ckpt = os.path.join(OUT_DIR, f"tworoute_records_dx{int(dx_um)}.jsonl")
    seed = os.path.join(OUT_DIR, "dtissue_confirm_dx12.jsonl")
    if not os.path.exists(ckpt) and os.path.exists(seed) and int(dx_um) == 12:
        shutil.copyfile(seed, ckpt)
        print(f"  seeded {ckpt} from the existing dx=12 records")

    points = [d_new]
    if per_drug:
        for k, v in scales.items():
            p = d_ship * v
            if all(abs(p - q) > 1e-6 * abs(q) for q in points):
                points.append(p)
    out_json = os.path.join(OUT_DIR, f"tworoute_sweep_dx{int(dx_um)}.json")
    res = dts.run(dx_um=dx_um, modes=("frozen", "recal"), out_json=out_json,
                  checkpoint=ckpt, explicit_points=points)

    def pick(mode, drug, d):
        cands = [r for r in res["records"] if r["mode"] == mode and r["drug"] == drug]
        if not cands:
            return None
        best = min(cands, key=lambda r: abs(r["D_tissue_cm2_s"] - d))
        rel = abs(best["D_tissue_cm2_s"] - d) / abs(d)
        assert rel < 1e-3, (mode, drug, d, best["D_tissue_cm2_s"], rel)
        return best

    print("\n  TABLE A -- doxorubicin-anchored UNIFORM correction (directly")
    print("  comparable with results/cellular_exclusion.md 5.5)")
    print(f"  {'mode':7s} {'drug':16s} {'ratio ship':>11s} {'ratio corr':>11s} "
          f"{'verdict ship':>13s} {'verdict corr':>13s} {'moved':>6s}")
    moves_uniform = []
    for mode in ("frozen", "recal"):
        for drug in res["drugs"]:
            a, b = pick(mode, drug, d_ship), pick(mode, drug, d_new)
            if a is None or b is None:
                continue
            moved = a["verdict"] != b["verdict"]
            moves_uniform.append({"mode": mode, "drug": drug,
                                  "ratio_ship": a["ratio_exact"],
                                  "ratio_corr": b["ratio_exact"],
                                  "verdict_ship": a["verdict"],
                                  "verdict_corr": b["verdict"], "moved": moved,
                                  "alpha_ship": a["alpha"], "alpha_corr": b["alpha"]})
            print(f"  {mode:7s} {drug:16s} {a['ratio_exact']:11.4f} "
                  f"{b['ratio_exact']:11.4f} {a['verdict']:>13s} "
                  f"{b['verdict']:>13s} {'YES' if moved else 'no':>6s}")
    n_uniform = sum(1 for m in moves_uniform if m["moved"])
    print(f"  verdicts that MOVE under the uniform correction: {n_uniform} of "
          f"{len(moves_uniform)}")

    moves_perdrug = []
    if per_drug:
        print("\n  TABLE B -- each drug at ITS OWN two-route D (what a real wiring")
        print("  would do). `frozen` is exact. `recal` reads alpha re-solved at that")
        print("  drug's own scale rather than at doxorubicin's; the two differ by")
        print("  the alpha column, which is printed so the approximation is visible.")
        print(f"  {'mode':7s} {'drug':16s} {'scale':>8s} {'ratio ship':>11s} "
              f"{'ratio corr':>11s} {'verdict':>9s}->{'verdict':>9s} {'alpha':>9s} "
              f"{'moved':>6s}")
        for mode in ("frozen", "recal"):
            for drug in res["drugs"]:
                a = pick(mode, drug, d_ship)
                b = pick(mode, drug, d_ship * scales[drug])
                if a is None or b is None:
                    continue
                moved = a["verdict"] != b["verdict"]
                moves_perdrug.append({"mode": mode, "drug": drug,
                                      "scale": scales[drug],
                                      "ratio_ship": a["ratio_exact"],
                                      "ratio_corr": b["ratio_exact"],
                                      "verdict_ship": a["verdict"],
                                      "verdict_corr": b["verdict"], "moved": moved,
                                      "alpha_corr": b["alpha"]})
                print(f"  {mode:7s} {drug:16s} {scales[drug]:8.4f} "
                      f"{a['ratio_exact']:11.4f} {b['ratio_exact']:11.4f} "
                      f"{a['verdict']:>9s}->{b['verdict']:>9s} {b['alpha']:>9.4f} "
                      f"{'YES' if moved else 'no':>6s}")
        n_pd = sum(1 for m in moves_perdrug if m["moved"])
        print(f"  verdicts that MOVE under the per-drug correction: {n_pd} of "
              f"{len(moves_perdrug)}")

    return {"scales": scales, "scale_spread": spread, "D_ship_cm2_s": d_ship,
            "D_dox_corrected_cm2_s": d_new,
            "alpha_published": float(ref["alpha_published"]),
            "moves_uniform": moves_uniform, "moves_perdrug": moves_perdrug,
            "out_json": out_json}


# =========================================================================== #
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all",
                    choices=("all", "prereg", "limits", "pruijn", "drugs", "blast"))
    ap.add_argument("--dx", type=float, default=12.0)
    ap.add_argument("--no-per-drug", action="store_true")
    a = ap.parse_args()

    t0 = time.time()
    out = {}
    if a.stage in ("all", "prereg"):
        out["prereg"] = stage_prereg()
    if a.stage in ("all", "limits"):
        out["limits"] = stage_limits()
    if a.stage in ("all", "pruijn"):
        out["pruijn"] = stage_pruijn()
    if a.stage in ("all", "drugs"):
        out["drugs"] = stage_drugs()
    if a.stage in ("all", "blast"):
        out["blast"] = stage_blast(dx_um=a.dx, per_drug=not a.no_per_drug)

    os.makedirs(OUT_DIR, exist_ok=True)
    # NOTE: deliberately NOT f"tworoute_{stage}.json" -- that would collide with,
    # and overwrite, the pre-registration file this run is supposed to be checked
    # against. (It did exactly that once.)
    path = os.path.join(OUT_DIR, f"tworoute_run_{a.stage}.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(f"\nwrote {path}   ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
