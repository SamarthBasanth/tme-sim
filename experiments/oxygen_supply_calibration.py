"""DERIVING `OxygenConfig.supply` FROM THE PDE IT APPEARS IN.

`results/oxygen_level_shift.md` established that `supply = 0.05` was never
re-derived after commit `ceb96c0` changed what it multiplies, and that a factor
of 2 on it walks `F_state` across essentially its whole attainable range. This
script does the derivation the way a transport coefficient should be derived:
identify what the constant IS dimensionally, find what the physics literature
says that quantity is, and then MEASURE what the implied value does.

NOTHING IS CHANGED. `tmesim/` is untouched; `OxygenConfig.supply` keeps its
shipped default of 0.05 everywhere. The only intervention is a harness-level
wrapper around `physics.oxygen.solve_oxygen` that replaces the `cfg` it is
handed with `dataclasses.replace(cfg, supply=S)`. The solve itself is
unmodified library code and `supply` is a public, documented field of a public
dataclass, so this is the same class of intervention `oxygen_level_shift.py`
used for `mode='kde'`.

WHY THE FIELD-SCALING TRICK DOES NOT WORK HERE (read before "simplifying" this)
------------------------------------------------------------------------------
`oxygen_level_shift.py`'s `cliff` stage swept the source by scaling the vascular
FIELD, which is equivalent to scaling `supply` -- but only DOWNWARD.
`solve_oxygen` does `vasc = np.clip(vasc, 0.0, 1.0)` before forming the source,
and the shipped discrete field is already exactly 1.0 inside vessels, so any
upward scale factor is silently clipped away (that report's §3.4 measured the
resulting oxygen field to be bit-identical). This investigation is entirely
about values of `supply` ABOVE the shipped one, so it must move `supply` itself.

STAGES (each checkpoints to results/sweep/oxsupply_<stage>.json):

  algebra   -- pure arithmetic + a self-contained Krogh-cylinder reference
               solve. Establishes the dimensional identity of `supply`, the
               three competing rates it must be compared against, and how far
               the shipped value is from the Dirichlet (pinned-vessel) limit
               that the classical Krogh formulation assumes. No section needed.
  conductance -- THE DERIVATION. Closed form for the model's own intravascular
               mass transfer coefficient, inverted against Secomb 2020's
               measured cerebral K_TO2,cap. Self-tests the coefficient set by
               reproducing that paper's own reported K_TO2,tiss from its own
               analytic formula, and ABORTS if it does not.
  joint     -- derived vs shipped `supply` crossed with the literature band for
               `consumption`, because the two are not separately identifiable.
  sinkmatch -- source pinned at the derived `supply`, sink swept until the
               implied CMRO2 matches Watabe 2017's PET-measured glioma TMRO2;
               reports the resulting pO2 distribution against the measured
               human-GBM pO2 literature.
  ladder    -- `supply` swept over 8 decades on the real GBM section at
               dx = 12 um. Records the two halves of the blood-to-tissue
               oxygen drop separately (transvascular vs. tissue-diffusive),
               mean tumour O2, the three state bands and `F_state`.
  sections  -- shipped vs derived `supply` on BOTH real sections (GBM brain and
               breast), same quantities.
  grid      -- the dx ladder from `results/grid_convergence_gbm.md` re-run at
               the derived `supply`, to test whether the non-monotone `F_state`
               survives it.
  killcal   -- the unmodified `validation.kill_calibration.run()`/`analyse()`
               under the derived `supply`. Reported LAST and separately, per
               the task's ordering rule.

Run: ./.venv/bin/python experiments/oxygen_supply_calibration.py <stage> ...
"""
from __future__ import annotations

import dataclasses
import json
import os
import sys
import time

import numpy as np

_T0 = time.time()
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SWEEP = os.path.join(ROOT, "results", "sweep")

BRAIN = "data_cache/xenium_gbm"
BREAST = "data_cache/xenium_breast_janesick/rep1"
DX = 12.0

#: The value this report derives (see results/oxygen_supply_calibration.md).
#: NOT written to any config -- used only to label the arms below.
#: supply = K_TO2_cap / (alpha * f_vessel); see :func:`stage_conductance`.
DERIVED_SUPPLY = 194.0

# --------------------------------------------------------------------------- #
# Literature constants used by the closed-form derivation. Every one of these
# is a primary-source value; see results/oxygen_supply_calibration.md for the
# citation and the UNVERIFIED flags.
# --------------------------------------------------------------------------- #

#: Intravascular (capillary) oxygen mass transfer coefficient per unit TISSUE
#: volume, mouse cerebral cortex reference parameter set. Units
#: cm^3 O2 / (100 cm^3 tissue . min . mmHg).
#: Secomb TW, "The mass transfer coefficient for oxygen transport from blood to
#: tissue in cerebral cortex", J Cereb Blood Flow Metab 2020;40(8):1634-1646.
#: PMID 31423930, DOI 10.1177/0271678X19870068.
K_TO2_CAP = 1.213
#: Extravascular (tissue) partner coefficient from the same reference set.
K_TO2_TISS = 0.997
#: Secomb's reference geometry: capillary radius, tissue cylinder radius,
#: capillary length density, tissue O2 diffusivity and O2 solubility.
SECOMB = {
    "rc_um": 2.0, "rt_um": 25.0, "Ncap_per_mm2": 509.3, "HT": 0.281,
    "D_um2_s": 2410.0, "alpha_cm3O2_per_cm3_mmHg": 3.89e-5,
    "CMRO2_cm3O2_per_100cm3_min": 6.0,
}
#: O2 solubility band, cm^3 O2 / (cm^3 tissue . mmHg). Low end is the classical
#: Bunsen-coefficient figure for plasma at 37 C (0.003 mL O2/100 mL/mmHg); the
#: mid value is the one already derived in-repo from Grimes et al. 2014's Omega
#: (results/sweep/oxygen_ocr_literature_second_pass_*.md Step 1); the high end
#: is Secomb's own.
ALPHA_BAND = (3.0e-5, 3.3e-5, 3.89e-5)


def log(msg: str) -> None:
    print(f"[{time.time() - _T0:7.1f}s] {msg}", flush=True)


def save(stage: str, payload: dict) -> str:
    os.makedirs(SWEEP, exist_ok=True)
    path = os.path.join(SWEEP, f"oxsupply_{stage}.json")
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2, default=float)
    log(f"checkpoint -> {os.path.relpath(path, ROOT)}")
    return path


# --------------------------------------------------------------------------- #
# The one intervention.
# --------------------------------------------------------------------------- #
_PATCH_TARGETS = (
    "tmesim.physics.oxygen",
    "tmesim.reconstruct.real_fields",
    "tmesim.sim.engine",
    "tmesim.validation.kill_calibration",
    "tmesim.validation.oxygen_benchmark",
    "tmesim.validation.immune_sink_impact",
)


class force_supply:
    """Rebind `solve_oxygen` to a wrapper that overrides `cfg.supply` only.

    Every module that calls `solve_oxygen` imported the name at import time, so
    rebinding it inside `physics.oxygen` alone would not be seen. The wrapper
    forwards to the UNMODIFIED function with `dataclasses.replace(cfg,
    supply=S)`; no numerical code is replaced, and `supply=None` is a no-op
    passthrough so the shipped arm runs the identical code path.
    """

    def __init__(self, supply):
        self.supply = None if supply is None else float(supply)
        self._saved = {}

    def __enter__(self):
        import importlib

        from tmesim.physics.oxygen import solve_oxygen as _orig

        S = self.supply

        def _wrapped(micro, cfg):
            if S is not None:
                cfg = dataclasses.replace(cfg, supply=S)
            return _orig(micro, cfg)

        fn = _orig if S is None else _wrapped
        for mod_name in _PATCH_TARGETS:
            mod = importlib.import_module(mod_name)
            if hasattr(mod, "solve_oxygen"):
                self._saved[mod_name] = mod.solve_oxygen
                mod.solve_oxygen = fn
        return self

    def __exit__(self, *exc):
        import importlib

        for mod_name, saved in self._saved.items():
            setattr(importlib.import_module(mod_name), "solve_oxygen", saved)
        return False


# --------------------------------------------------------------------------- #
# Measurement shared by every section stage.
# --------------------------------------------------------------------------- #
def measure(micro) -> dict:
    """Everything this investigation reads off one reconstructed section.

    The two numbers that matter for the DERIVATION (as opposed to the
    consequences) are `o2_at_vessel_pixels` and the drop decomposition: the
    model's blood compartment is O = 1 BY CONSTRUCTION (the source term is
    `supply*vasc*(1-O)`, which vanishes only at O = 1), so
    `1 - o2_at_vessel_pixels` is the share of the total blood-to-tissue oxygen
    drop the model attributes to the transvascular step, and
    `o2_at_vessel_pixels - mean_o2_tumour` is the share it attributes to
    diffusion through tissue. The classical Krogh picture puts essentially all
    of it in the second term.
    """
    from tmesim.cells.agents import update_cell_states
    from tmesim.core import (CellState, CellType, Microenvironment,
                             OxygenConfig, PDConfig)
    from tmesim.physics.oxygen import REFERENCE_PO2_MMHG

    pdc = PDConfig()
    mult = {"HYPOXIC": pdc.hypoxia_resistance,
            "QUIESCENT": pdc.quiescent_resistance, "VIABLE": 1.0}

    oxy = np.asarray(micro.fields.oxygen, dtype=np.float64)
    vasc = np.asarray(micro.fields.vascular_density, dtype=np.float64)
    tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
    vessel = (vasc > 0) & tissue

    work = Microenvironment(
        domain=micro.domain, fields=micro.fields, cells=micro.cells.copy(),
        patient_id=micro.patient_id, metadata=micro.metadata)
    update_cell_states(work, oxy, OxygenConfig())
    st = np.asarray(work.cells.state)
    tum = np.asarray(micro.cells.cell_type) == int(CellType.TUMOR)
    n_tum = int(tum.sum())
    fr = {s.name: float((st[tum] == int(s)).sum()) / n_tum
          for s in (CellState.HYPOXIC, CellState.QUIESCENT, CellState.VIABLE)}

    col, row = micro.domain.to_index(np.asarray(micro.cells.x)[tum],
                                     np.asarray(micro.cells.y)[tum])
    o_tum = oxy[row, col]

    o_vess = float(oxy[vessel].mean()) if vessel.any() else float("nan")
    o_tis = float(oxy[tissue].mean())
    mean_tum = float(o_tum.mean())

    return {
        "n_tumour_cells": n_tum,
        "grid": [int(v) for v in micro.domain.shape],
        "dx_um": float(micro.domain.dx),
        "vessel_pixel_fraction_of_tissue": float(vessel.sum() / max(tissue.sum(), 1)),
        "mean_vascular_tissue": float(vasc[tissue].mean()),
        "o2_at_vessel_pixels": o_vess,
        "o2_max_anywhere": float(oxy.max()),
        "mean_o2_tissue": o_tis,
        "mean_o2_tumour": mean_tum,
        "mean_o2_tumour_mmhg": mean_tum * REFERENCE_PO2_MMHG,
        "o2_at_vessel_pixels_mmhg": o_vess * REFERENCE_PO2_MMHG,
        # Drop decomposition (see docstring).
        "drop_transvascular": 1.0 - o_vess,
        "drop_tissue_diffusive": o_vess - mean_tum,
        "transvascular_share_of_drop": float(
            (1.0 - o_vess) / max(1.0 - mean_tum, 1e-12)),
        "band_fractions": fr,
        "f_state": float(sum(fr[k] * mult[k] for k in fr)),
        "hypoxic_frac_tumour": fr["HYPOXIC"],
        "tumour_frac_below_2p5mmhg": float(
            (o_tum * REFERENCE_PO2_MMHG < 2.5).mean()),
        "tumour_frac_below_10mmhg": float(
            (o_tum * REFERENCE_PO2_MMHG < 10.0).mean()),
        "tumour_o2_mmhg_percentiles": {
            str(p): float(np.percentile(o_tum * REFERENCE_PO2_MMHG, p))
            for p in (1, 5, 10, 25, 50, 75, 90, 99)},
    }


# --------------------------------------------------------------------------- #
# Stage: algebra -- the derivation itself, no section required
# --------------------------------------------------------------------------- #
def stage_algebra() -> dict:
    """Dimensional identity of `supply`, the competing rates, and the Krogh limit.

    Part A is arithmetic on the shipped constants. Part B is a self-contained
    axisymmetric Krogh-cylinder solve (built here, not imported) that answers
    the one question the arithmetic cannot: at what `supply` does the
    distributed-source formulation this module uses become numerically
    indistinguishable from the Dirichlet (O = 1 at the vessel wall) boundary
    condition the classical Krogh/Grimes formulation imposes?
    """
    import scipy.sparse as sp
    from scipy.sparse.linalg import spsolve

    from tmesim.core import OxygenConfig
    from tmesim.reconstruct.vessel_network import VESSEL_RADIUS_UM

    cfg = OxygenConfig()
    D = float(cfg.diffusivity)
    R = float(VESSEL_RADIUS_UM)

    out = {
        "shipped": {"supply_s^-1": float(cfg.supply),
                    "consumption_s^-1": float(cfg.consumption),
                    "diffusivity_um2_s": D,
                    "vessel_radius_um": R},
        "competing_rates_s^-1": {
            # Rate at which diffusion empties a vessel-sized region.
            "D_over_R2": D / R ** 2,
            # Rate at which diffusion empties one grid cell, at each dx used.
            "D_over_dx2_dx8": D / 8.0 ** 2,
            "D_over_dx2_dx12": D / 12.0 ** 2,
            "D_over_dx2_dx30": D / 30.0 ** 2,
            # The sink, at rho_tum = 1 and rho_tum = 8 (module self-test range).
            "sink_rho1": float(cfg.consumption) * 1.0,
            "sink_rho8": float(cfg.consumption) * 8.0,
        },
        "lengths_um": {
            # sqrt(D/supply): how far oxygen diffuses before the source
            # equilibrates it. Must be << R for the vessel to pin its own O2.
            "source_equilibration_length_shipped": float(np.sqrt(D / cfg.supply)),
            "sink_diffusion_length_rho1": float(np.sqrt(D / (cfg.consumption * 1.0))),
            "sink_diffusion_length_rho8": float(np.sqrt(D / (cfg.consumption * 8.0))),
            "vessel_radius": R,
        },
    }
    out["damkohler_shipped"] = {
        "supply_over_D_over_R2": float(cfg.supply / (D / R ** 2)),
        "supply_over_D_over_dx2_dx12": float(cfg.supply / (D / 12.0 ** 2)),
    }
    # Wall-only permeability estimate: P = D_wall/t for an endothelium of
    # thickness t, S/V = 2/R for a cylinder, supply = P * 2/R.
    out["wall_only_estimate"] = {
        f"t_{t}um": {"P_um_s": D / t, "supply_s^-1": (D / t) * 2.0 / R}
        for t in (0.2, 0.5, 1.0, 2.0)
    }

    # ---------------- Part B: Krogh-cylinder convergence to Dirichlet -------- #
    # Axisymmetric steady state on r in [0, R_t], uniform sink Q for r > R:
    #   (1/r) d/dr ( r D dO/dr ) + supply*1{r<R}*(1-O) - Q*1{r>R}*O = 0
    # Neumann at r = R_t (Krogh's no-flux tissue boundary) and regularity at 0.
    # Compared against the Dirichlet problem O(R) = 1 on r in [R, R_t].
    def krogh(supply, Q, R_t=60.0, R_v=R, n=4000):
        dr = R_t / n
        r = (np.arange(n) + 0.5) * dr           # cell centres
        rf = np.arange(1, n) * dr               # interior faces
        a_p = np.zeros(n)
        rows, cols, vals = [], [], []
        for k in range(n - 1):
            c = D * rf[k] / (dr * dr * r[k])
            rows.append(k); cols.append(k + 1); vals.append(c)
            a_p[k] += c
            c2 = D * rf[k] / (dr * dr * r[k + 1])
            rows.append(k + 1); cols.append(k); vals.append(c2)
            a_p[k + 1] += c2
        src = np.where(r < R_v, float(supply), 0.0)
        snk = np.where(r >= R_v, float(Q), 0.0)
        rows.extend(range(n)); cols.extend(range(n))
        vals.extend(-(a_p + src + snk))
        A = sp.coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsc()
        o = spsolve(A, -src)
        return r, np.clip(o, 0.0, 1.0)

    def krogh_dirichlet(Q, R_t=60.0, R_v=R, n=4000):
        # Same discretisation, but O pinned to 1 for r < R_v.
        BIG = 1e12
        return krogh(BIG, Q, R_t=R_t, R_v=R_v, n=n)

    krogh_out = {}
    for rho in (1.0, 8.0):
        Q = float(cfg.consumption) * rho
        r_ref, o_ref = krogh_dirichlet(Q)
        tis = r_ref >= R
        ref_mean = float(np.average(o_ref[tis], weights=r_ref[tis]))
        rows = []
        for s in (0.05, 0.5, 5.0, 35.56, 100.0, 356.0, 1e3, 1e4, 1e5):
            r_, o_ = krogh(s, Q)
            m = float(np.average(o_[tis], weights=r_[tis]))
            rows.append({
                "supply": float(s),
                "o2_at_vessel_axis": float(o_[0]),
                "o2_at_vessel_wall": float(o_[np.searchsorted(r_, R)]),
                "mean_tissue_o2": m,
                "fraction_of_dirichlet_mean": m / ref_mean,
                "relative_error_vs_dirichlet": abs(m - ref_mean) / ref_mean,
            })
        krogh_out[f"rho_tum={rho:g}"] = {
            "sink_Q_s^-1": Q,
            "dirichlet_reference_mean_tissue_o2": ref_mean,
            "dirichlet_reference_wall_o2": float(o_ref[np.searchsorted(r_ref, R)]),
            "sweep": rows,
        }
        log(f"algebra/krogh rho={rho:g}: Dirichlet mean tissue O2 = {ref_mean:.4f}; "
            f"shipped supply=0.05 gives {rows[0]['mean_tissue_o2']:.4f} "
            f"({100*rows[0]['fraction_of_dirichlet_mean']:.1f}% of it)")
    out["krogh_cylinder"] = krogh_out
    return out


# --------------------------------------------------------------------------- #
# Stage: conductance -- the closed-form derivation
# --------------------------------------------------------------------------- #
def stage_conductance() -> dict:
    """Derive `supply` in closed form, then check the model's TISSUE side too.

    THE DERIVATION. `solve_oxygen`'s source, per unit tissue volume averaged
    over tissue, is ``supply * f_v * (1 - O_vessel)`` in units of O/s, where
    ``f_v`` is the vessel-pixel area fraction of tissue. Converting O to a
    physical oxygen volume via ``c = alpha * p`` and ``p = P_ref * O`` makes
    the transvascular flux per unit tissue volume

        M = supply * f_v * alpha * P_ref * (1 - O_vessel)
          = supply * f_v * alpha * (P_ref - p_vessel)

    and the transvascular partial-pressure drop is exactly
    ``dP_cap = P_ref - p_vessel``. So the model's own intravascular mass
    transfer coefficient is

        K_cap(model) = M / dP_cap = supply * f_v * alpha

    -- the ``(1 - O_vessel)`` cancels, so this is INDEPENDENT of `consumption`,
    of `diffusivity`, and of the solved field. It is a pure statement about the
    coefficient. Setting it equal to Secomb's measured ``K_TO2_cap`` inverts to

        supply = K_TO2_cap / (alpha * f_v)

    which is the value this report derives. Note the one uncomfortable
    consequence, which is reported rather than hidden: ``f_v`` is a property of
    the SECTION and of the vessel-painting convention, so `supply` as a global
    config scalar is only section-invariant to the extent the painted vessel
    area fraction is.

    THE TISSUE-SIDE CHECK. The same conversion applied to the Dirichlet-limit
    solve gives ``K_tiss(model) = M / dP_tiss``, which can be compared against
    Secomb's ``K_TO2_tiss``. That one DOES depend on `consumption`, on the
    geometry and on the grid, so it is a check on the rest of the model rather
    than part of the derivation -- and it is reported as such.
    """
    import numpy as np

    from tmesim.core import OxygenConfig
    from tmesim.physics.oxygen import REFERENCE_PO2_MMHG, _tumor_density
    from tmesim.reconstruct.real_pipeline import reconstruct_real
    from tmesim.reconstruct.vessel_network import VESSEL_RADIUS_UM

    def k_tiss_krogh(rc_um, rt_um, D_um2_s, alpha):
        """Secomb 2020 eq. for the Krogh-cylinder TISSUE mass transfer coefficient.

            K_tiss = 8 pi D alpha Ncap (1-d^2)^2 / (-4 ln d - 3 + 4 d^2 - d^4)

        with d = rc/rt and Ncap = 1/(pi rt^2). Returned in the paper's units,
        cm^3 O2 / (100 cm^3 . min . mmHg). This function is SELF-TESTED below
        against the paper's own reported 0.997 for its own reference geometry,
        which is how the fetched coefficient is verified rather than trusted.
        """
        d = rc_um / rt_um
        g = -4.0 * np.log(d) - 3.0 + 4.0 * d ** 2 - d ** 4
        n_cap_per_cm2 = 1.0 / (np.pi * (rt_um * 1e-4) ** 2)
        D_cm2_s = D_um2_s * 1e-8
        return (8.0 * np.pi * D_cm2_s * alpha * n_cap_per_cm2
                * (1.0 - d ** 2) ** 2 / g) * 100.0 * 60.0

    cfg = OxygenConfig()
    _self_test = float(k_tiss_krogh(SECOMB["rc_um"], SECOMB["rt_um"],
                                    SECOMB["D_um2_s"],
                                    SECOMB["alpha_cm3O2_per_cm3_mmHg"]))
    if abs(_self_test - K_TO2_TISS) / K_TO2_TISS > 0.01:
        raise AssertionError(
            f"Krogh K_tiss formula does not reproduce Secomb's reported "
            f"{K_TO2_TISS}: got {_self_test}")
    log(f"conductance: Krogh K_tiss formula reproduces Secomb's reported "
        f"{K_TO2_TISS} as {_self_test:.4f} -- coefficient set verified")

    out = {
        "krogh_formula_self_test": {
            "secomb_reported_K_tiss": K_TO2_TISS,
            "formula_value": _self_test,
            "relative_error": float(abs(_self_test - K_TO2_TISS) / K_TO2_TISS),
        },
        "K_TO2_cap_cm3O2_per_100cm3_min_mmHg": K_TO2_CAP,
        "K_TO2_tiss_cm3O2_per_100cm3_min_mmHg": K_TO2_TISS,
        "secomb_reference": dict(SECOMB),
        "reference_po2_mmhg": float(REFERENCE_PO2_MMHG),
        "vessel_radius_um": float(VESSEL_RADIUS_UM),
    }

    # Secomb's own capillary volume fraction, two independent ways.
    f_v_secomb_krogh = (SECOMB["rc_um"] / SECOMB["rt_um"]) ** 2
    f_v_secomb_density = (SECOMB["Ncap_per_mm2"] / 1e6) * np.pi * SECOMB["rc_um"] ** 2
    out["secomb_capillary_volume_fraction"] = {
        "from_krogh_radii": float(f_v_secomb_krogh),
        "from_length_density": float(f_v_secomb_density),
    }

    # K in per-second units.
    k_cap_s = K_TO2_CAP / 100.0 / 60.0     # cm^3O2/(cm^3 . s . mmHg)
    out["K_TO2_cap_cm3O2_per_cm3_s_mmHg"] = float(k_cap_s)
    # Transvascular conductance per unit tissue volume, in the model's O units.
    out["G_transvascular_s^-1_by_alpha"] = {
        f"alpha={a:.2e}": float(k_cap_s / a) for a in ALPHA_BAND}
    # supply implied by Secomb's OWN geometry (a per-vessel-volume rate).
    out["supply_from_secomb_own_geometry_s^-1"] = {
        f"alpha={a:.2e}": float(k_cap_s / (a * f_v_secomb_krogh))
        for a in ALPHA_BAND}

    out["sections"] = {}
    for label, sample in (("gbm", BRAIN), ("breast", BREAST)):
        log(f"conductance: reconstructing {label} (shipped config)")
        micro = reconstruct_real("xenium", sample=sample, dx_um=DX)
        vasc = np.asarray(micro.fields.vascular_density, dtype=np.float64)
        tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
        f_v = float((np.clip(vasc, 0, 1) * tissue).sum() / max(tissue.sum(), 1))
        rho = _tumor_density(micro)
        rho_bar = float(rho[tissue].mean())

        sec = {
            "vessel_area_fraction_of_tissue": f_v,
            "mean_rho_tum_per_144um2": rho_bar,
            "supply_implied_s^-1": {
                f"alpha={a:.2e}": float(k_cap_s / (a * f_v)) for a in ALPHA_BAND},
            # Same thing expressed the other way: what K_cap does the SHIPPED
            # supply correspond to?
            "K_cap_of_shipped_supply": {
                f"alpha={a:.2e}": float(cfg.supply * f_v * a * 100.0 * 60.0)
                for a in ALPHA_BAND},
        }
        # Theory value of the model's OWN tissue conductance: an equivalent
        # Krogh cylinder with the model's vessel radius and the section's own
        # vessel area fraction (rt = rc / sqrt(f_v)).
        rt_equiv = float(VESSEL_RADIUS_UM / np.sqrt(f_v))
        sec["krogh_equivalent"] = {
            "rc_um": float(VESSEL_RADIUS_UM),
            "rt_um": rt_equiv,
            "delta_rc_over_rt": float(VESSEL_RADIUS_UM / rt_equiv),
            "K_tiss_theory": {
                f"alpha={a:.2e}": float(
                    k_tiss_krogh(VESSEL_RADIUS_UM, rt_equiv, cfg.diffusivity, a))
                for a in ALPHA_BAND},
        }
        # Scaling check: Secomb's per-vessel-volume supply rescaled by radius
        # (supply = 2*P/R scales as 1/R).
        sec["supply_from_radius_scaling_s^-1"] = {
            f"alpha={a:.2e}": float(
                (k_cap_s / (a * f_v_secomb_krogh)) * SECOMB["rc_um"] / VESSEL_RADIUS_UM)
            for a in ALPHA_BAND}

        # --- tissue side, from the Dirichlet-limit solve on this section ---
        with force_supply(1e5):
            micro_d = reconstruct_real("xenium", sample=sample, dx_um=DX)
        m = measure(micro_d)
        alpha_mid = ALPHA_BAND[1]
        # Physical consumption flux per unit tissue volume implied by the model.
        oxy_d = np.asarray(micro_d.fields.oxygen, dtype=np.float64)
        M_model_O_per_s = float(
            (cfg.consumption * rho * oxy_d)[tissue].mean())      # O/s
        M_phys = M_model_O_per_s * alpha_mid * REFERENCE_PO2_MMHG  # cm3O2/cm3/s
        # Both halves averaged over the SAME support (tissue), so K_tiss is a
        # self-consistent conductance rather than a mixed-support ratio.
        dP_tiss = (m["o2_at_vessel_pixels"] - m["mean_o2_tissue"]) * REFERENCE_PO2_MMHG
        sec["dirichlet_limit"] = {
            "mean_o2_tissue": m["mean_o2_tissue"],
            "mean_o2_tumour": m["mean_o2_tumour"],
            "o2_at_vessel_pixels": m["o2_at_vessel_pixels"],
            "drop_tissue_diffusive_mmhg": float(dP_tiss),
            "implied_CMRO2_cm3O2_per_100cm3_min": float(M_phys * 100.0 * 60.0),
            "K_tiss_model_cm3O2_per_100cm3_min_mmHg": (
                float(M_phys / dP_tiss * 100.0 * 60.0) if dP_tiss > 0 else None),
            "K_tiss_model_over_secomb": (
                float(M_phys / dP_tiss * 100.0 * 60.0 / K_TO2_TISS)
                if dP_tiss > 0 else None),
        }
        out["sections"][label] = sec
        log(f"  {label}: f_v={f_v:.5f}  supply implied "
            f"{sec['supply_implied_s^-1'][f'alpha={ALPHA_BAND[1]:.2e}']:.1f} s^-1  "
            f"| shipped supply corresponds to K_cap="
            f"{sec['K_cap_of_shipped_supply'][f'alpha={ALPHA_BAND[1]:.2e}']:.5f} "
            f"(Secomb {K_TO2_CAP})")
        log(f"  {label}: Dirichlet-limit implied CMRO2 = "
            f"{sec['dirichlet_limit']['implied_CMRO2_cm3O2_per_100cm3_min']:.3f} "
            f"cm3O2/100cm3/min (Secomb reference 6.0; Watabe C6 glioma 1.9-2.8); "
            f"K_tiss(model)/K_tiss(Secomb) = "
            f"{sec['dirichlet_limit']['K_tiss_model_over_secomb']}")
        save("conductance_partial", out)
    return out


# --------------------------------------------------------------------------- #
# Stage: joint -- supply and consumption together
# --------------------------------------------------------------------------- #
def stage_joint() -> dict:
    """Cross the derived `supply` with the literature band for `consumption`.

    `supply` and `consumption` are not independently identifiable from the
    solved field alone -- only their ratio sets the level. This project has
    already bracketed `consumption` (results/oxygen_sensitivity_band.md;
    results/sweep/oxygen_ocr_literature_second_pass_*.md, which converts
    Watabe et al. 2017's PET-measured C6-glioma TMRO2 into 0.240-0.505 s^-1,
    12-25x the shipped 0.02). This stage asks the question neither study could:
    with the SOURCE side pinned by its own physics, does the sink band produce
    a GBM oxygen field that matches the measured human GBM pO2 literature?
    """
    import dataclasses as _dc

    from tmesim.core import OxygenConfig
    from tmesim.reconstruct.real_pipeline import reconstruct_real

    out = {"sample": BRAIN, "dx_um": DX, "cells": {}}
    for sup_label, sup in (("shipped_0.05", 0.05), ("derived_194", DERIVED_SUPPLY)):
        for con_label, con in (("shipped_0.02", 0.02),
                               ("watabe_lo_0.240", 0.240),
                               ("watabe_hi_0.505", 0.505),
                               ("diepart_hi_0.005", 0.005)):
            log(f"joint: supply={sup_label} consumption={con_label}")
            with _force_cfg(supply=sup, consumption=con):
                micro = reconstruct_real("xenium", sample=BRAIN, dx_um=DX)
            rec = measure(micro)
            rec["supply"] = float(sup)
            rec["consumption"] = float(con)
            out["cells"][f"{sup_label}|{con_label}"] = rec
            log(f"  meanO2(tum)={rec['mean_o2_tumour']:.4f} "
                f"({rec['mean_o2_tumour_mmhg']:.2f} mmHg)  "
                f"hypox={rec['hypoxic_frac_tumour']:.4f}  "
                f"<10mmHg={rec['tumour_frac_below_10mmhg']:.4f}  "
                f"F_state={rec['f_state']:.4f}")
            save("joint_partial", out)
    _ = OxygenConfig, _dc
    return out


def stage_sinkmatch() -> dict:
    """With the SOURCE pinned by its own physics, what sink matches measured CMRO2?

    `stage_conductance` measures the oxygen consumption flux the shipped
    `consumption` actually produces on the real section, in the same units as
    clinical CMRO2, and finds it is ~140x below Watabe et al. 2017's
    PET-measured C6-glioma TMRO2. That comparison is only meaningful once the
    source is not starving the tissue, so it is re-done here at the derived
    supply, sweeping `consumption` until the implied CMRO2 lands in the
    measured band -- and then asking the question that actually matters: does
    the resulting field reproduce the measured human-GBM pO2 distribution?
    """
    from tmesim.core import OxygenConfig
    from tmesim.physics.oxygen import REFERENCE_PO2_MMHG, _tumor_density
    from tmesim.reconstruct.real_pipeline import reconstruct_real

    alpha = ALPHA_BAND[1]
    out = {"sample": BRAIN, "dx_um": DX, "supply": DERIVED_SUPPLY,
           "alpha_used": alpha,
           "watabe_TMRO2_band_cm3O2_per_100cm3_min": [1.9, 2.8],
           "secomb_cortex_reference_cm3O2_per_100cm3_min": 6.0,
           "rungs": {}}
    for con in (0.02, 0.1, 0.5, 1.0, 2.5, 5.0, 10.0, 20.0, 50.0):
        log(f"sinkmatch: consumption={con:g} at supply={DERIVED_SUPPLY:g}")
        with _force_cfg(supply=DERIVED_SUPPLY, consumption=con):
            micro = reconstruct_real("xenium", sample=BRAIN, dx_um=DX)
        rec = measure(micro)
        tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
        rho = _tumor_density(micro)
        oxy = np.asarray(micro.fields.oxygen, dtype=np.float64)
        M = float((con * rho * oxy)[tissue].mean()) * alpha * REFERENCE_PO2_MMHG
        rec["consumption"] = float(con)
        rec["implied_CMRO2_cm3O2_per_100cm3_min"] = float(M * 100.0 * 60.0)
        out["rungs"][f"c{con:g}"] = rec
        log(f"  CMRO2={rec['implied_CMRO2_cm3O2_per_100cm3_min']:.3f}  "
            f"meanO2(tum)={rec['mean_o2_tumour_mmhg']:.2f} mmHg  "
            f"median={rec['tumour_o2_mmhg_percentiles']['50']:.2f} mmHg  "
            f"hypox={rec['hypoxic_frac_tumour']:.4f}  "
            f"<10mmHg={rec['tumour_frac_below_10mmhg']:.4f}  "
            f"F_state={rec['f_state']:.4f}")
        save("sinkmatch_partial", out)
    _ = OxygenConfig
    return out


class _force_cfg(force_supply):
    """Like :class:`force_supply` but overrides `consumption` too."""

    def __init__(self, supply=None, consumption=None):
        super().__init__(supply)
        self.consumption = None if consumption is None else float(consumption)

    def __enter__(self):
        import importlib

        from tmesim.physics.oxygen import solve_oxygen as _orig

        S, C = self.supply, self.consumption

        def _wrapped(micro, cfg):
            kw = {}
            if S is not None:
                kw["supply"] = S
            if C is not None:
                kw["consumption"] = C
            if kw:
                cfg = dataclasses.replace(cfg, **kw)
            return _orig(micro, cfg)

        for mod_name in _PATCH_TARGETS:
            mod = importlib.import_module(mod_name)
            if hasattr(mod, "solve_oxygen"):
                self._saved[mod_name] = mod.solve_oxygen
                mod.solve_oxygen = _wrapped
        return self


# --------------------------------------------------------------------------- #
# Stage: ladder -- supply swept on the real GBM section
# --------------------------------------------------------------------------- #
SUPPLY_LADDER = (0.005, 0.05, 0.5, 1.0, 5.0, 10.0, 35.56, 100.0, 356.0,
                 1000.0, 1e4, 1e5)


def stage_ladder() -> dict:
    from tmesim.reconstruct.real_pipeline import reconstruct_real

    out = {"sample": BRAIN, "dx_um": DX, "rungs": {}}
    for s in SUPPLY_LADDER:
        log(f"ladder: supply = {s:g}")
        with force_supply(s):
            micro = reconstruct_real("xenium", sample=BRAIN, dx_um=DX)
        rec = measure(micro)
        rec["supply"] = float(s)
        out["rungs"][f"s{s:g}"] = rec
        log(f"  supply={s:<10g} O2@vessel={rec['o2_at_vessel_pixels']:.4f}  "
            f"meanO2(tum)={rec['mean_o2_tumour']:.4f}  "
            f"H/Q/V {rec['band_fractions']['HYPOXIC']:.3f}/"
            f"{rec['band_fractions']['QUIESCENT']:.3f}/"
            f"{rec['band_fractions']['VIABLE']:.3f}  "
            f"F_state={rec['f_state']:.4f}")
        save("ladder_partial", out)
    return out


# --------------------------------------------------------------------------- #
# Stage: sections -- both real sections, shipped vs derived
# --------------------------------------------------------------------------- #
def stage_sections() -> dict:
    from tmesim.reconstruct.real_pipeline import reconstruct_real

    out = {"dx_um": DX, "derived_supply": DERIVED_SUPPLY, "sections": {}}
    for label, sample in (("gbm", BRAIN), ("breast", BREAST)):
        out["sections"][label] = {}
        for arm, s in (("shipped_0.05", None), ("derived", DERIVED_SUPPLY),
                       ("dirichlet_limit_1e5", 1e5)):
            log(f"sections: {label} / {arm}")
            with force_supply(s):
                micro = reconstruct_real("xenium", sample=sample, dx_um=DX)
            rec = measure(micro)
            rec["supply"] = float(s) if s is not None else 0.05
            out["sections"][label][arm] = rec
            log(f"  {label}/{arm}: O2@vessel={rec['o2_at_vessel_pixels']:.4f} "
                f"meanO2(tum)={rec['mean_o2_tumour']:.4f} "
                f"({rec['mean_o2_tumour_mmhg']:.2f} mmHg) "
                f"F_state={rec['f_state']:.4f} "
                f"hypox={rec['hypoxic_frac_tumour']:.4f}")
            save("sections_partial", out)
    return out


# --------------------------------------------------------------------------- #
# Stage: grid -- does the non-monotone F_state survive the derived supply?
# --------------------------------------------------------------------------- #
def stage_grid() -> dict:
    from tmesim.reconstruct.real_pipeline import reconstruct_real

    out = {"sample": BRAIN, "derived_supply": DERIVED_SUPPLY, "ladder": {}}
    for dx in (30.0, 20.0, 16.0, 12.0, 10.0, 8.0):
        for arm, s in (("shipped", None), ("derived", DERIVED_SUPPLY)):
            log(f"grid: dx={dx:g} arm={arm}")
            with force_supply(s):
                micro = reconstruct_real("xenium", sample=BRAIN, dx_um=dx)
            rec = measure(micro)
            rec["arm"] = arm
            out["ladder"][f"dx{dx:g}_{arm}"] = rec
            log(f"  dx={dx:g} {arm}: meanO2(tum)={rec['mean_o2_tumour']:.4f} "
                f"F_state={rec['f_state']:.4f}")
            save("grid_partial", out)

    for arm in ("shipped", "derived"):
        vals = [out["ladder"][k]["f_state"] for k in out["ladder"]
                if k.endswith("_" + arm)]
        out[f"f_state_span_{arm}"] = float(max(vals) - min(vals))
        out[f"f_state_values_{arm}"] = vals
    return out


# --------------------------------------------------------------------------- #
# Stage: killcal -- the consequence, reported last
# --------------------------------------------------------------------------- #
def stage_killcal() -> dict:
    from tmesim.validation import kill_calibration as kc

    out = {"arms": {}}
    # (label, supply, consumption). consumption=None leaves the shipped 0.02.
    # The last two arms are the JOINTLY calibrated configurations from
    # `stage_sinkmatch` -- the source pinned by Secomb's K_TO2_cap and the sink
    # pinned by Watabe's PET-measured C6-glioma TMRO2 band.
    for arm, s, c in (("derived_supply_only", DERIVED_SUPPLY, None),
                      ("joint_cmro2_lo", DERIVED_SUPPLY, 10.0),
                      ("joint_cmro2_hi", DERIVED_SUPPLY, 20.0)):
        log(f"killcal: kill_calibration.run() under supply={s} consumption={c}")
        t0 = time.time()
        with _force_cfg(supply=s, consumption=c):
            res = kc.run()
            analysis = kc.analyse(res)
        rows = {}
        for r in res.get("drugs", res.get("rows", [])):
            if not isinstance(r, dict):
                continue
            name = r.get("drug") or r.get("name")
            rows[name] = {k: r.get(k) for k in
                          ("drug", "ratio_exact", "f_penetration", "f_state",
                           "implied_ic50_fold", "verdict") if k in r}
        out["arms"][arm] = {
            "supply": float(s),
            "consumption": (float(c) if c is not None else 0.02),
            "rows": rows,
            "state_breakdown": res.get("state"),
            "verdicts": _verdicts(analysis),
            "wall_s": time.time() - t0,
        }
        log(f"killcal {arm}: {out['arms'][arm]['verdicts']}")
        save("killcal_partial", out)
    return out


def _verdicts(analysis) -> dict:
    found = {}

    def walk(o, path=""):
        if isinstance(o, dict):
            for k, v in o.items():
                if k in ("verdict", "ratio_exact", "f_state", "f_penetration",
                         "predicted_ratio", "implied_ic50_fold", "band_lo",
                         "band_hi", "lit_ratio_lo", "lit_ratio_hi"):
                    found[f"{path}.{k}"] = v
                else:
                    walk(v, f"{path}.{k}")
        elif isinstance(o, (list, tuple)):
            for i, v in enumerate(o):
                walk(v, f"{path}[{i}]")

    walk(analysis)
    return found


STAGES = {
    "algebra": stage_algebra,
    "conductance": stage_conductance,
    "joint": stage_joint,
    "sinkmatch": stage_sinkmatch,
    "ladder": stage_ladder,
    "sections": stage_sections,
    "grid": stage_grid,
    "killcal": stage_killcal,
}


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv == ["all"]:
        argv = ["algebra", "conductance", "ladder", "joint", "sections",
                "grid", "killcal"]
    for stage in argv:
        if stage not in STAGES:
            print(f"unknown stage {stage!r}; expected {sorted(STAGES)}")
            return 2
        log(f"=== STAGE {stage} ===")
        save(stage, STAGES[stage]())
        log(f"=== STAGE {stage} done ===")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
