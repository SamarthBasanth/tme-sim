"""Is p99 bath-matching the RIGHT normalising convention, or just A convention?

The finding this script exists to adjudicate
--------------------------------------------
``tmesim/validation/kill_calibration.decompose`` sets the plasma level so that
the **99th percentile** of the per-tumour-cell exposure field sits exactly at the
bath concentration, and only then applies the Hill map. Three independent
investigations converged on the consequence:

* ``results/boundary_layer_resolution.md`` sec. 4.4/5.3 -- with exchange OFF the
  domain MEAN of the exposure field has EXACTLY zero grid error (a discrete
  conservation identity, max |rel err| 2.1e-11 over 180 points), but the
  ``mean/p99`` shape statistic is off by **+8.3 % at dx = 12 um** and refinement
  walks it down only to +5.4 % (dx = 8), +2.2 % (4), +1.1 % (2).
* ``results/dtissue_sensitivity.md`` sec. 5/7.3 -- a 3.07x change in mean tumour
  exposure moves the scored ratio by 0.997-1.009x.
* ``results/grid_convergence_gbm.md`` sec. 7 -- halving dt moves the exposure
  field by 1.2e-04 and the scored ratio by 0.0003 %.

All three files give the same explanation: p99 bath-matching removes
multiplicative SCALE error exactly, so only SHAPE error survives.

What this script adds
---------------------
STAGE A -- **the algebra, verified numerically rather than asserted.** The scored
    ratio is shown to be a functional of the NORMALISED SHAPE VECTOR
    ``u = e / N(e)`` alone, for any normaliser ``N`` that is positively
    homogeneous of degree 1. Exact scale invariance is therefore a property of
    the whole homogeneous family, **not of p99**. Measured to machine precision
    over 12 orders of magnitude of scaling, on the real published exposure
    fields, for every candidate convention.

STAGE B -- **grid error of every candidate convention** on the SAME controlled
    slab that found the +8.3 % (``experiments/boundary_layer_resolution.py``,
    which has a closed-form steady state). Nothing is re-derived: this module
    imports ``slab_exact``, ``slab_numeric`` and ``slab_exact_reference`` from
    that script so the numbers are commensurable with the study that raised the
    question.

STAGE C -- **sensitivity to the diffusivity perturbation**, on the real breast
    section at the PUBLISHED dx = 12 um grid, reusing
    ``experiments/dtissue_sensitivity.py``'s own ``build_drug`` / ``e_unit_for``
    so the perturbation is the one that file measured. The reference arm is the
    published run cache and is asserted against the published record.

STAGE D -- **verdict consequences**, computed LAST and reported separately, so
    that no convention can be selected on its PASS count. Nine published rows
    (3 GBM + 6 breast) recomputed under every convention; SEVEN of them carry a
    verdict. Breast paclitaxel and breast 5-fluorouracil have one drug-matched
    study each, so their min-max band is a zero-width point and they are NOT
    TESTABLE under every convention. Their ratios are still swept -- the
    scale-invariance result is a property of the ratio and needs no band -- but
    they are counted in no PASS denominator and are never re-scored against the
    broad cross-tumour-type band. The verdict for every row comes from
    ``kill_calibration.score_ratio_row``; this module does not decide it. See
    ``results/zero_width_bands.md`` and
    ``results/fallback_rule_propagation.md``.

STAGE E -- **estimator stability**: bootstrap and subsample behaviour of each
    normaliser on the real per-cell exposure distributions. An order statistic in
    the extreme tail is not a consistent estimator of anything; the mean is.

STAGE F -- **what is the bath analogue, physically?** Locates the cells that
    actually sit in a perfused voxel on each real section and reports where their
    mean exposure falls on the percentile scale, against where each candidate
    convention falls. Prices the spheroid side of the same analogy from geometry.

STAGE G -- **does the slab number transfer to tissue?** The slab and both
    published tissue dx ladders all solve dx = 12 AND dx = 8, so the same dx step
    can be compared on both. Tests the claim in
    ``results/boundary_layer_resolution.md`` sec. 8 that the +8.3 % applies to
    "every 3-D:2-D verdict in this repo". **It does not** -- the slab error is
    drug-independent by construction, the tissue error is drug-dependent, and the
    two shortest-L drugs move in the opposite direction.

HARD RULE OBSERVED: no convention is selected by its verdict count. Stages B, C,
E, F and G are the selection evidence; stage D is a consequence report.

THE KNOWN TRAP (``tmesim/validation/propagation.py``)
-----------------------------------------------------
``kill_calibration.unit_transport`` calls ``simulate_drug_transport`` DIRECTLY,
so ``PDEConfig`` flags are silent no-ops there and ``kill_calibration.py`` is
declared ``INERT``. **This script never passes ``exchange_k`` itself** and so is
NOT a wiring point: stage B's exchange-ON arm is produced by calling
``boundary_layer_resolution.slab_numeric(exchange=True)``, and that module is
already a declared ``WiringPoint`` with the assertion it requires. Stage C runs
the shipped exchange-OFF path, which is what every published verdict ran on.

Run:
    ./.venv/bin/python experiments/p99_convention.py --stage A
    ./.venv/bin/python experiments/p99_convention.py --stage B
    ./.venv/bin/python experiments/p99_convention.py --stage C   # ~40 min
    ./.venv/bin/python experiments/p99_convention.py --stage D
    ./.venv/bin/python experiments/p99_convention.py --stage E
    ./.venv/bin/python experiments/p99_convention.py --stage F
    ./.venv/bin/python experiments/p99_convention.py --stage G   # needs stage B
Writes:
    results/sweep/p99_algebra.json
    results/sweep/p99_gridslab.json
    results/sweep/p99_dtissue_dx12.json   (+ .jsonl checkpoint)
    results/sweep/p99_verdicts.json
    results/sweep/p99_stability.json
    results/sweep/p99_perivascular.json
    results/sweep/p99_transfer.json
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, ".."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

OUT_DIR = os.path.join(_ROOT, "results", "sweep")
_T0 = time.time()


def log(msg: str) -> None:
    print(f"[p99 {time.strftime('%H:%M:%S')} +{time.time()-_T0:7.1f}s] {msg}",
          flush=True)


# --------------------------------------------------------------------------- #
# THE CANDIDATE CONVENTIONS
# --------------------------------------------------------------------------- #
#: Every candidate normaliser ``N``. A convention is admissible as a
#: "bath-matching" rule only if it is POSITIVELY HOMOGENEOUS OF DEGREE 1, i.e.
#: ``N(lambda*e) = lambda*N(e)`` for every ``lambda > 0`` -- that is precisely
#: the algebraic condition that makes the scored ratio invariant to a
#: multiplicative error in the exposure field (stage A proves it). ``fixed`` is
#: deliberately included and is NOT homogeneous: it is the negative control that
#: shows what the protection is worth.
#:
#: ``kind`` groups them by estimator class, which is what stage E is about:
#:   order   -- an order statistic (percentile / max)
#:   moment  -- a smooth integral functional of the field
#:   none    -- no normalisation at all
NORMALISERS: Tuple[Dict, ...] = (
    dict(tag="p50", kind="order", homogeneous=True,
         label="median", fn=lambda e: float(np.percentile(e, 50.0))),
    dict(tag="p90", kind="order", homogeneous=True,
         label="p90", fn=lambda e: float(np.percentile(e, 90.0))),
    dict(tag="p95", kind="order", homogeneous=True,
         label="p95", fn=lambda e: float(np.percentile(e, 95.0))),
    dict(tag="p99", kind="order", homogeneous=True,
         label="p99 (SHIPPED)", fn=lambda e: float(np.percentile(e, 99.0))),
    dict(tag="p99.9", kind="order", homogeneous=True,
         label="p99.9", fn=lambda e: float(np.percentile(e, 99.9))),
    dict(tag="max", kind="order", homogeneous=True,
         label="max", fn=lambda e: float(np.max(e))),
    dict(tag="mean", kind="moment", homogeneous=True,
         label="mean", fn=lambda e: float(np.mean(e))),
    dict(tag="rms", kind="moment", homogeneous=True,
         label="rms (L2)", fn=lambda e: float(np.sqrt(np.mean(np.square(e))))),
    dict(tag="L4", kind="moment", homogeneous=True,
         label="L4 norm", fn=lambda e: float(np.mean(np.power(e, 4.0)) ** 0.25)),
    dict(tag="sauter", kind="moment", homogeneous=True,
         label="exposure-weighted mean  mean(e^2)/mean(e)",
         fn=lambda e: float(np.mean(np.square(e)) / max(np.mean(e), 1e-300))),
    # --- conditional tail means (CVaR / expected shortfall) ----------------
    # The physical object the bath-matching convention is reaching for is a
    # SHELL of well-perfused cells, not one quantile: a spheroid's outer layer
    # is a set of cells, all of them near the bath. A conditional tail mean is
    # that object -- it keeps the upper-tail meaning of a high percentile but
    # is an INTEGRAL over the tail rather than a single order statistic, so it
    # is neither set by one grid cell nor a step function of the grid.
    dict(tag="tail99", kind="tail", homogeneous=True,
         label="mean of the top 1 % (CVaR 99)",
         fn=lambda e: float(np.mean(e[e >= np.percentile(e, 99.0)]))),
    dict(tag="tail95", kind="tail", homogeneous=True,
         label="mean of the top 5 % (CVaR 95)",
         fn=lambda e: float(np.mean(e[e >= np.percentile(e, 95.0)]))),
    dict(tag="tail90", kind="tail", homogeneous=True,
         label="mean of the top 10 % (CVaR 90)",
         fn=lambda e: float(np.mean(e[e >= np.percentile(e, 90.0)]))),
    dict(tag="fixed", kind="none", homogeneous=False,
         label="fixed physical reference (NO normalisation)", fn=None),
)
NORM_BY_TAG: Dict[str, Dict] = {n["tag"]: n for n in NORMALISERS}
HOMOGENEOUS_TAGS: Tuple[str, ...] = tuple(
    n["tag"] for n in NORMALISERS if n["homogeneous"])
SHIPPED_TAG = "p99"


def normalise(e: np.ndarray, tag: str, fixed_ref: Optional[float] = None
              ) -> float:
    """The reference exposure ``N(e)`` under convention ``tag``.

    ``fixed`` needs an externally supplied constant: the whole point of that
    convention is that the reference does NOT come from the field being scored,
    so it must be frozen once (here: the published run's own p99) and then
    reused unchanged at every perturbed point.
    """
    if tag == "fixed":
        if fixed_ref is None:
            raise ValueError("the 'fixed' convention needs a frozen reference")
        return float(fixed_ref)
    return float(NORM_BY_TAG[tag]["fn"](np.asarray(e, dtype=np.float64)))


# --------------------------------------------------------------------------- #
# The scored functional, generalised over the normaliser
# --------------------------------------------------------------------------- #
def decompose_general(e_unit, sens, resist, c_bath: float, ic50: float,
                      hill: float, max_kill: float = 1.0, *, tag: str = "p99",
                      fixed_ref: Optional[float] = None) -> Dict:
    """``kill_calibration.decompose`` with the normaliser left free.

    This is a GENERALISATION, not a re-implementation: with ``tag='p99'`` it is
    asserted bit-identical to the shipped ``decompose`` in
    :func:`assert_matches_shipped`, which every stage calls before using it.
    The engine's own Hill map (``kill_calibration.p_hill``) is used unmodified.
    """
    from tmesim.validation.kill_calibration import p_hill

    e_unit = np.asarray(e_unit, dtype=np.float64)
    sens = np.asarray(sens, dtype=np.float64)
    resist = np.asarray(resist, dtype=np.float64)

    ref_unit = normalise(e_unit, tag, fixed_ref)
    if ref_unit <= 0.0:
        raise ValueError(f"reference exposure is zero under {tag}")
    plasma = float(c_bath) / ref_unit
    E = plasma * e_unit

    p3 = p_hill(E, ic50, hill, max_kill)
    p2 = float(p_hill(c_bath, ic50, hill, max_kill))

    kill_3d = float(np.mean(p3 * sens * resist))
    kill_2d = float(p2 * np.mean(sens))
    return {
        "normaliser": tag,
        "ref_unit_uM_per_uM_plasma": ref_unit,
        "plasma_required_uM": plasma,
        "f_penetration": float(np.mean(p3) / p2) if p2 > 0 else float("nan"),
        "f_state": float(np.mean(resist)),
        "ratio_exact": float(kill_3d / kill_2d) if kill_2d > 0 else float("nan"),
        "kill_3d_model": kill_3d,
        "kill_2d_model": kill_2d,
        "E_median_over_bath": float(np.median(E) / c_bath),
    }


def assert_matches_shipped(e_unit, sens, resist, c_bath, ic50, hill,
                           max_kill) -> Dict:
    """Refuse to proceed unless ``decompose_general(tag='p99')`` reproduces the
    SHIPPED ``kill_calibration.decompose`` bit-for-bit.

    Without this guard every number below would be a comparison against this
    file's own arithmetic rather than against the published convention.
    """
    from tmesim.validation.kill_calibration import decompose

    ship = decompose(e_unit, sens, resist, c_bath, ic50, hill, max_kill)
    mine = decompose_general(e_unit, sens, resist, c_bath, ic50, hill, max_kill,
                             tag="p99")
    out = {}
    for k in ("ref_unit_uM_per_uM_plasma", "plasma_required_uM",
              "f_penetration", "f_state", "ratio_exact", "kill_3d_model",
              "kill_2d_model"):
        a, b = float(ship[k]), float(mine[k])
        assert a == b or (a == 0.0 and b == 0.0), (
            f"decompose_general(p99) does not reproduce the shipped decompose "
            f"on {k}: {a!r} vs {b!r}. Refusing to continue -- every alternative "
            f"convention below would be scored against the wrong baseline.")
        out[k] = a
    return out


# --------------------------------------------------------------------------- #
# Loading the published, cached exposure fields
# --------------------------------------------------------------------------- #
GBM_CACHE = os.path.join(_ROOT, "data_cache", "kill_calibration_runs.json")
BREAST_CACHE = os.path.join(_ROOT, "data_cache",
                            "kill_calibration_breast_runs.json")
#: The six drugs ``results/kill_calibration_breast.md`` calls formally testable.
BREAST_TESTABLE = ("doxorubicin", "paclitaxel", "docetaxel", "epirubicin",
                   "5-fluorouracil", "cisplatin")
GBM_DRUGS = ("temozolomide", "cisplatin", "doxorubicin")


def load_gbm_panel() -> Dict:
    """The published GBM panel: cached e_unit + a fresh engine prologue.

    The transport solves are read from the run cache written by the published
    ``kill_calibration.run`` (the same file that module itself resumes from);
    the cell-state labels are re-derived by running the engine's OWN prologue,
    which is deterministic and involves no PDE transport.
    """
    from tmesim.core import PDConfig, PDEConfig, PKConfig, SimConfig
    from tmesim.reconstruct.real_pipeline import (apply_drug_sensitivity,
                                                  reconstruct_real)
    from tmesim.validation.calibration import calibrate_alpha
    from tmesim.validation.kill_calibration import (BATH_MULTIPLES, DT_S,
                                                    DX_UM, EXVIVO_TMZ_UM,
                                                    NOMINAL_HALF_LIFE_S,
                                                    NOMINAL_PEAK_UM,
                                                    SAVE_EVERY, SEED, WINDOW_S,
                                                    compile_drug,
                                                    prepare_tissue,
                                                    resistance_vector)
    from tmesim.core import Microenvironment

    cache = json.load(open(GBM_CACHE))
    log("GBM: reconstructing the real Xenium GBM section at dx = 12 um ...")
    micro = reconstruct_real(source="xenium", dx_um=float(DX_UM))
    pde = PDEConfig(dt=DT_S, total_time=WINDOW_S, convection=True,
                    save_every=SAVE_EVERY)
    pk = PKConfig(plasma_peak=NOMINAL_PEAK_UM, half_life=NOMINAL_HALF_LIFE_S,
                  n_doses=1, dose_interval=3600.0 * 24)
    pd_cfg = PDConfig()
    cfg = SimConfig(pde=pde, pk=pk, pd=pd_cfg, seed=SEED)
    work, labeled, _vel = prepare_tissue(micro, cfg)
    tumor = labeled.tumor_mask()
    state_t = np.asarray(labeled.state)[tumor]
    alpha = float(calibrate_alpha())

    out: Dict = {"tissue": "gbm", "dx_um": float(DX_UM), "alpha": alpha,
                 "n_tumor": int(tumor.sum()), "drugs": {},
                 "_work": work, "_labeled": labeled, "_tumor": tumor}
    for name in GBM_DRUGS:
        drug = compile_drug(name, alpha)
        key = f"{name}|{DX_UM}|{DT_S}|{WINDOW_S}|{SAVE_EVERY}"
        e_unit = np.asarray(cache[key]["e_unit"], dtype=np.float64)
        m2 = Microenvironment(domain=work.domain, fields=work.fields,
                              cells=labeled.copy(), patient_id=work.patient_id,
                              metadata=work.metadata)
        apply_drug_sensitivity(m2, name, mode="uniform")
        sens = np.asarray(m2.cells.sensitivity, dtype=np.float64)[tumor]
        c_bath = (float(EXVIVO_TMZ_UM) if name == "temozolomide"
                  else 1.0 * float(drug.ic50))
        out["drugs"][name] = {
            "e_unit": e_unit, "sens": sens,
            "resist": resistance_vector(state_t, pd_cfg, drug=drug),
            "c_bath_uM": c_bath, "ic50_uM": float(drug.ic50),
            "hill": float(drug.hill), "max_kill": float(pd_cfg.max_kill_rate),
            "c_bath_label": ("100 uM (ex-vivo assay)" if name == "temozolomide"
                             else "1xIC50"),
        }
    _ = BATH_MULTIPLES
    return out


def load_breast_panel() -> Dict:
    """The published breast panel: cached e_unit + a fresh engine prologue."""
    from tmesim.core import PDConfig, PDEConfig, PKConfig, SimConfig
    from tmesim.reconstruct.real_pipeline import reconstruct_real
    from tmesim.validation.calibration import calibrate_alpha
    from tmesim.validation.kill_calibration import (prepare_tissue,
                                                    resistance_vector)

    from experiments.kill_calibration_breast import (BREAST_SAMPLE_PATH, DT_S,
                                                     DX_UM,
                                                     NOMINAL_HALF_LIFE_S,
                                                     NOMINAL_PEAK_UM,
                                                     SAVE_EVERY, SEED,
                                                     WINDOW_S,
                                                     apply_drug_sensitivity_on_copy,
                                                     compile_breast_drug)

    cache = json.load(open(BREAST_CACHE))
    log("breast: reconstructing the real Xenium breast section at dx = 12 um ...")
    micro = reconstruct_real(source="xenium", sample=BREAST_SAMPLE_PATH,
                             dx_um=float(DX_UM))
    pde = PDEConfig(dt=DT_S, total_time=WINDOW_S, convection=True,
                    save_every=SAVE_EVERY)
    pk = PKConfig(plasma_peak=NOMINAL_PEAK_UM, half_life=NOMINAL_HALF_LIFE_S,
                  n_doses=1, dose_interval=3600.0 * 24)
    pd_cfg = PDConfig()
    cfg = SimConfig(pde=pde, pk=pk, pd=pd_cfg, seed=SEED)
    work, labeled, velocity = prepare_tissue(micro, cfg)
    tumor = labeled.tumor_mask()
    state_t = np.asarray(labeled.state)[tumor]
    alpha = float(calibrate_alpha())

    out: Dict = {"tissue": "breast", "dx_um": float(DX_UM), "alpha": alpha,
                 "n_tumor": int(tumor.sum()), "drugs": {},
                 "_work": work, "_labeled": labeled, "_tumor": tumor,
                 "_velocity": velocity, "_pde": pde, "_pd": pd_cfg,
                 "_state_t": state_t}
    for name in BREAST_TESTABLE:
        drug, _mol, _info = compile_breast_drug(name, alpha)
        key = (f"{name}|{DX_UM}|{DT_S}|{WINDOW_S}|{SAVE_EVERY}"
               f"|ic50={float(drug.ic50):.6g}")
        e_unit = np.asarray(cache[key]["e_unit"], dtype=np.float64)
        m2 = apply_drug_sensitivity_on_copy(work, labeled, name)
        sens = np.asarray(m2.cells.sensitivity, dtype=np.float64)[tumor]
        out["drugs"][name] = {
            "e_unit": e_unit, "sens": sens,
            "resist": resistance_vector(state_t, pd_cfg, drug=drug),
            "c_bath_uM": 1.0 * float(drug.ic50), "ic50_uM": float(drug.ic50),
            "hill": float(drug.hill), "max_kill": float(pd_cfg.max_kill_rate),
            "c_bath_label": "1xIC50",
        }
    _ = SEED
    return out


# --------------------------------------------------------------------------- #
# STAGE A -- the algebra
# --------------------------------------------------------------------------- #
#: Multiplicative field errors to test scale-invariance against. Deliberately
#: spanning 12 orders of magnitude: if the invariance is exact it costs nothing
#: to demand it over an absurd range, and if it is only approximate this finds
#: out immediately.
SCALE_FACTORS: Tuple[float, ...] = (1e-6, 1e-3, 0.1, 0.5, 2.0, 10.0, 1e3, 1e6)


def stage_a(panels: Sequence[Dict]) -> Dict:
    """What p99 bath-matching preserves and destroys, measured not asserted.

    A1  exact invariance under a multiplicative (SCALE) field error, for every
        degree-1-homogeneous normaliser, and the failure of the same test for
        the non-homogeneous ``fixed`` control.
    A2  the response to a pure SHAPE error of matched size, so the two can be
        compared on the same axis.
    A3  the response to a SPATIALLY VARYING multiplicative error -- the class
        that a grid refinement actually produces, and which is neither pure
        scale nor pure shape.
    """
    out: Dict = {"scale_factors": list(SCALE_FACTORS), "panels": []}
    rng = np.random.default_rng(20260809)

    for panel in panels:
        prec: Dict = {"tissue": panel["tissue"], "drugs": []}
        for name, d in panel["drugs"].items():
            e, sens, resist = d["e_unit"], d["sens"], d["resist"]
            cb, ic, hh, mk = (d["c_bath_uM"], d["ic50_uM"], d["hill"],
                              d["max_kill"])
            base_check = assert_matches_shipped(e, sens, resist, cb, ic, hh, mk)
            fixed_ref = float(np.percentile(e, 99.0))

            row: Dict = {"drug": name, "c_bath_uM": cb,
                         "shipped_ratio": base_check["ratio_exact"],
                         "conventions": {}}
            for tag in [n["tag"] for n in NORMALISERS]:
                base = decompose_general(e, sens, resist, cb, ic, hh, mk,
                                         tag=tag, fixed_ref=fixed_ref)
                # --- A1: pure multiplicative SCALE error -------------------
                worst_scale = 0.0
                for s in SCALE_FACTORS:
                    pert = decompose_general(e * s, sens, resist, cb, ic, hh,
                                             mk, tag=tag, fixed_ref=fixed_ref)
                    rel = abs(pert["ratio_exact"] / base["ratio_exact"] - 1.0)
                    worst_scale = max(worst_scale, rel)
                # --- A2: pure SHAPE error, size-matched --------------------
                # A monotone power-law warp e -> e^(1+eps), renormalised to hold
                # the MEAN fixed, so it changes only the distribution's shape.
                eps = 0.02
                warped = np.power(e / e.mean(), 1.0 + eps)
                warped *= e.mean() / warped.mean()
                shape_resp = abs(
                    decompose_general(warped, sens, resist, cb, ic, hh, mk,
                                      tag=tag, fixed_ref=fixed_ref
                                      )["ratio_exact"] / base["ratio_exact"]
                    - 1.0)
                # --- A3: spatially varying multiplicative error ------------
                # 2 % iid multiplicative noise: the class a grid change makes.
                noise = np.exp(rng.normal(0.0, 0.02, size=e.size))
                noisy_resp = abs(
                    decompose_general(e * noise, sens, resist, cb, ic, hh, mk,
                                      tag=tag, fixed_ref=fixed_ref
                                      )["ratio_exact"] / base["ratio_exact"]
                    - 1.0)
                row["conventions"][tag] = {
                    "ref_unit": base["ref_unit_uM_per_uM_plasma"],
                    "ratio_exact": base["ratio_exact"],
                    "f_penetration": base["f_penetration"],
                    "worst_rel_change_under_pure_scale": worst_scale,
                    "rel_change_under_2pct_shape_warp": shape_resp,
                    "rel_change_under_2pct_iid_multiplicative": noisy_resp,
                }
            prec["drugs"].append(row)
        out["panels"].append(prec)
    return out


# --------------------------------------------------------------------------- #
# STAGE B -- grid error on the controlled slab with the exact steady solution
# --------------------------------------------------------------------------- #
#: The dx ladder. dx = 12 um is the shipped grid; the rest bracket it so the
#: convergence ORDER of each convention can be read off, not just its value.
SLAB_DX: Tuple[float, ...] = (1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0, 24.0,
                              30.0)
#: Dimensionless bath multiples ``c_bath / IC50`` at which the Hill-mapped
#: functional is scored. 1.0 is the breast panel's reporting concentration;
#: 0.165 = 100/606 is temozolomide's (the ex-vivo assay concentration over its
#: library IC50). 10.0 probes the saturated end.
BATH_MULTIPLES_SLAB: Tuple[float, ...] = (0.165, 1.0, 10.0)
SLAB_HILL = 1.2


#: Number of synthetic "tumour cells" laid down on the slab for the
#: population-sampled variant of stage B (see :func:`stage_b`). Fixed physical
#: positions, INDEPENDENT of dx, because that is what the real pipeline has:
#: ``kill_calibration`` reads the AUC field at ``domain.to_index(cells.x,
#: cells.y)``, i.e. it nearest-neighbour samples the grid field at a cell set
#: that does not move when the grid changes.
SLAB_N_CELLS = 20_000


def _functionals(c: np.ndarray, hill: float,
                 multiples: Sequence[float]) -> Dict[str, Dict[str, float]]:
    """Every convention's scored functionals on one concentration profile.

    ``G_lin`` = ``mean(C)/N(C)`` is the low-dose / linear-Hill limit of
    ``f_penetration`` and, for ``N = p99``, is EXACTLY the shipped
    ``TransportResult.penetration_index``. ``fpen_<m>`` is the full Hill-mapped
    ``f_penetration`` at ``c_bath = m*IC50``, which is the quantity
    ``decompose`` actually scores.
    """
    from tmesim.validation.kill_calibration import p_hill

    c = np.asarray(c, dtype=np.float64)
    out: Dict[str, Dict[str, float]] = {}
    for n in NORMALISERS:
        if not n["homogeneous"]:
            continue          # 'fixed' has no meaning on a bare profile
        ref = float(n["fn"](c))
        vals = {"ref": ref, "G_lin": float(c.mean() / ref)}
        u = c / ref
        for m in multiples:
            vals[f"fpen_{m:g}"] = float(
                np.mean(p_hill(m * u, 1.0, hill)) / p_hill(m, 1.0, hill))
        out[n["tag"]] = vals
    return out


def stage_b(geometries: Sequence[str] = ("vessel", "wide_aligned"),
            delta_um: float = 15.0,
            dx_list: Sequence[float] = SLAB_DX) -> Dict:
    """Grid error of every convention on the closed-form two-region slab.

    Imported wholesale from ``experiments/boundary_layer_resolution.py`` -- the
    slab, its exact steady state, the cell-averaged reference and the shipped
    solver call are all that script's, so these numbers sit on the same axis as
    the +8.3 %/+17.8 % it reported.

    The exact reference for each functional is its CONTINUUM value, evaluated on
    a 4e5-point sampling of the closed form (the same construction
    ``blayer.slab_exact_penetration`` uses for the shipped penetration index),
    because that is the dx -> 0 limit the discrete statistic converges to.
    """
    import experiments.boundary_layer_resolution as blayer

    d = 140.41              # um^2/s, the shipping D_tissue
    l_out = 64.92           # um, the measured doxorubicin e-fold anchor
    kappa = d / l_out ** 2
    lx, cp = 480.0, 1.0
    ps = d / (delta_um ** 2) - kappa
    geo = {t: a for t, a in blayer.SLAB_GEOMETRIES}

    out: Dict = {"D_um2_s": d, "L_out_um": l_out, "kappa_per_s": kappa,
                 "Lx_um": lx, "delta_um": delta_um, "PS_per_s": ps,
                 "hill": SLAB_HILL, "bath_multiples": list(BATH_MULTIPLES_SLAB),
                 "geometries": {}}

    for tag in geometries:
        a_um = geo[tag]
        series: Dict = {"a_um": a_um, "exact": {}, "points": []}
        # Fixed physical "cell" positions -- dx-independent, exactly as the
        # tissue pipeline's cell set is.
        x_cells = (np.arange(SLAB_N_CELLS) + 0.5) * (lx / SLAB_N_CELLS)
        series["n_cells_sampled"] = int(SLAB_N_CELLS)

        for law, exchange in (("off", False), ("on", True)):
            k_in = (ps + kappa) if exchange else kappa
            xs = np.linspace(0.0, lx, 400_001)
            c_exact = blayer.slab_exact(xs, d=d, k_in=k_in, k_out=kappa,
                                        s_in=ps * cp, a=a_um, lx=lx)
            series["exact"][law] = _functionals(c_exact, SLAB_HILL,
                                                BATH_MULTIPLES_SLAB)
            # the population-sampled reference uses the SAME cell positions as
            # the numeric arm, so the finite-N sampling cancels exactly and the
            # only thing measured is the grid.
            c_exact_cells = blayer.slab_exact(x_cells, d=d, k_in=k_in,
                                              k_out=kappa, s_in=ps * cp,
                                              a=a_um, lx=lx)
            series["exact"][law + "_cells"] = _functionals(
                c_exact_cells, SLAB_HILL, BATH_MULTIPLES_SLAB)
        for dx in dx_list:
            if abs(round(lx / dx) * dx - lx) > 1e-9:
                continue
            rec: Dict = {"dx_um": float(dx), "a_over_dx": a_um / dx,
                         "interface_aligned":
                             abs(round(a_um / dx) * dx - a_um) < 1e-9,
                         "laws": {}}
            idx = np.clip((x_cells / dx).astype(int), 0,
                          int(round(lx / dx)) - 1)
            for law, exchange in (("off", False), ("on", True)):
                num = blayer.slab_numeric(d=d, ps=ps, kappa=kappa, cp=cp,
                                          a=a_um, lx=lx, dx=dx,
                                          exchange=exchange)
                for suffix, field in (("", num["c"]), ("_cells", num["c"][idx])):
                    got = _functionals(field, SLAB_HILL, BATH_MULTIPLES_SLAB)
                    ref = series["exact"][law + suffix]
                    errs: Dict = {}
                    for ntag, vals in got.items():
                        ex = ref[ntag]
                        errs[ntag] = {
                            k: (vals[k] - ex[k]) / ex[k] if ex[k] != 0 else
                            float("nan") for k in vals}
                    rec["laws"][law + suffix] = {
                        "values": got, "rel_err": errs,
                        "steady_drift": num["steady_drift"],
                        "mean": num["mean"]}
            series["points"].append(rec)
            log(f"  slab {tag} dx={dx:g}: "
                + " ".join(
                    f"{t}:{rec['laws']['off']['rel_err'][t]['G_lin']*100:+.2f}%"
                    for t in ("mean", "p95", "p99", "max")))
        out["geometries"][tag] = series
    return out


# --------------------------------------------------------------------------- #
# STAGE C -- the diffusivity perturbation on the published breast grid
# --------------------------------------------------------------------------- #
#: The two decisive D values ``results/dtissue_sensitivity.md`` sec. 7.3 re-ran on
#: the PUBLISHED dx = 12 um grid, and the alphas that file recorded for them.
#: Nothing here is a new choice: this is that file's own reduced point set.
D_SHIPPING_CM2S = 1.404089e-06
D_PERTURBED_CM2S = 4.559e-07


def stage_c(checkpoint: Optional[str] = None) -> Dict:
    """How does each convention respond to a real parameter perturbation?

    The reference arm is the PUBLISHED exposure cache (asserted, per drug,
    against the ``_e_unit_sum`` / ``_e_unit_first`` fingerprints recorded in
    ``results/sweep/dtissue_confirm_dx12.json``, so a silently stale cache
    cannot pass). Only the perturbed arms are solved: 6 drugs x 2 modes.

    ``frozen`` moves D with the published alpha and is the near-pure SHAPE
    perturbation; ``recal`` re-solves alpha against the same measured
    doxorubicin penetration anchor, which is what the codebase actually does,
    and is the near-pure SCALE perturbation (3.07x in mean exposure). Between
    them they bracket the two error classes stage A separates.
    """
    from tmesim.validation.kill_calibration import DT_S as _dt  # noqa: F401

    from experiments.dtissue_sensitivity import (build_drug, e_unit_for,
                                                 sweep_points)

    panel = load_breast_panel()
    work, labeled, tumor = panel["_work"], panel["_labeled"], panel["_tumor"]
    velocity, pde, pd_cfg = panel["_velocity"], panel["_pde"], panel["_pd"]
    alpha_pub = panel["alpha"]
    _ = sweep_points

    confirm = json.load(open(os.path.join(OUT_DIR,
                                          "dtissue_confirm_dx12.json")))
    fingerprint = {(r["mode"], round(r["D_tissue_cm2_s"], 12), r["drug"]):
                   (r["_e_unit_sum"], r["_e_unit_first"])
                   for r in confirm["records"]}

    out: Dict = {"dx_um": panel["dx_um"], "alpha_published": alpha_pub,
                 "D_shipping_cm2_s": D_SHIPPING_CM2S,
                 "D_perturbed_cm2_s": D_PERTURBED_CM2S,
                 "n_tumor": panel["n_tumor"], "arms": []}
    scale = D_PERTURBED_CM2S / D_SHIPPING_CM2S

    done: Dict = {}
    if checkpoint and os.path.exists(checkpoint):
        for line in open(checkpoint):
            try:
                r = json.loads(line)
            except Exception:
                continue
            done[(r["mode"], r["drug"])] = np.asarray(r["e_unit"],
                                                      dtype=np.float64)
        log(f"stage C: resumed {len(done)} perturbed solves from checkpoint")

    for mode in ("frozen", "recal"):
        for name in BREAST_TESTABLE:
            d = panel["drugs"][name]
            e_ref = d["e_unit"]
            # --- guard: is the cached reference arm the PUBLISHED one? ------
            fp = fingerprint.get((mode, round(D_SHIPPING_CM2S, 12), name))
            assert fp is not None, (
                f"no published dx = 12 record for {mode}/{name} at the shipping "
                f"D -- cannot verify the reference arm")
            s_got, f_got = float(e_ref.sum()), float(e_ref[0])
            assert (abs(s_got - fp[0]) <= 1e-6 * abs(fp[0])
                    and abs(f_got - fp[1]) <= 1e-12 * max(abs(fp[1]), 1e-30)), (
                f"{mode}/{name}: the cached reference exposure field does not "
                f"match the published dx = 12 record "
                f"(sum {s_got!r} vs {fp[0]!r}). Refusing to continue.")

            if (mode, name) in done:
                e_pert = done[(mode, name)]
                drug_p, _m, _i, alpha_p = build_drug(name, scale, mode,
                                                     alpha_pub)
            else:
                t0 = time.time()
                drug_p, _m, _i, alpha_p = build_drug(name, scale, mode,
                                                     alpha_pub)
                e_pert = e_unit_for(work, labeled, tumor, drug_p, pde, velocity)
                log(f"stage C: {mode}/{name} perturbed solve "
                    f"{time.time()-t0:.0f}s (alpha {alpha_p:.4f})")
                if checkpoint:
                    with open(checkpoint, "a") as fh:
                        fh.write(json.dumps({"mode": mode, "drug": name,
                                             "alpha": alpha_p,
                                             "e_unit": e_pert.tolist()}) + "\n")
            assert not np.array_equal(e_ref, e_pert), (
                f"{mode}/{name}: the perturbed exposure field is BIT-IDENTICAL "
                f"to the reference -- the D override did not reach the solver.")

            cb, ic, hh, mk = (d["c_bath_uM"], d["ic50_uM"], d["hill"],
                              d["max_kill"])
            assert_matches_shipped(e_ref, d["sens"], d["resist"], cb, ic, hh,
                                   mk)
            fixed_ref = float(np.percentile(e_ref, 99.0))
            row: Dict = {
                "mode": mode, "drug": name, "alpha_perturbed": alpha_p,
                "exposure_mean_fold": float(e_pert.mean() / e_ref.mean()),
                "exposure_p99_fold": float(np.percentile(e_pert, 99.0)
                                           / np.percentile(e_ref, 99.0)),
                "conventions": {},
            }
            for tag in [n["tag"] for n in NORMALISERS]:
                a = decompose_general(e_ref, d["sens"], d["resist"], cb, ic, hh,
                                      mk, tag=tag, fixed_ref=fixed_ref)
                b = decompose_general(e_pert, d["sens"], d["resist"], cb, ic,
                                      hh, mk, tag=tag, fixed_ref=fixed_ref)
                fold = b["ratio_exact"] / a["ratio_exact"]
                row["conventions"][tag] = {
                    "ratio_reference": a["ratio_exact"],
                    "ratio_perturbed": b["ratio_exact"],
                    "ratio_fold": float(fold),
                    "abs_ln_fold": float(abs(math.log(fold))),
                }
            out["arms"].append(row)
            log(f"stage C: {mode}/{name} exposure x"
                f"{row['exposure_mean_fold']:.3f}  p99 ratio fold "
                f"{row['conventions']['p99']['ratio_fold']:.4f}")
    return out


# --------------------------------------------------------------------------- #
# STAGE D -- verdict consequences (computed LAST, reported separately)
# --------------------------------------------------------------------------- #
def _band_entry(spec: Optional[Dict], solid: Dict, c_bath: float, ic50: float,
                hill: float) -> Dict:
    """One drug's band record: the drug-matched band, the broad context band,
    and BOTH bridges, so the consumer can hand them to
    ``kill_calibration.score_ratio_row`` rather than choosing between them here.

    Nothing in this module decides which band scores a row. It used to
    (``used = band if band_is_drug_matched else band_solid``), and that copy of
    the retired fallback rule is what put two broad-band substitutions into this
    document's §5 and §6 tables. See ``results/fallback_rule_propagation.md``.
    """
    from tmesim.validation.kill_calibration import bridge_to_kill_ratio_band
    band = (bridge_to_kill_ratio_band(spec, c_bath, ic50, hill)
            if spec else None)
    band_solid = bridge_to_kill_ratio_band(solid, c_bath, ic50, hill)
    return {"band": band, "band_solid": band_solid,
            "specific_bridge": spec, "solid_bridge": solid,
            "specific_n": int(spec["n"]) if spec else 0}


def gbm_bands(panel: Dict) -> Dict:
    """The published GBM literature bands, from the frozen literature tables.

    All three GBM drugs have n >= 5 drug-matched rows, so none is touched by the
    single-study rule -- but the decision is still made in one place.
    """
    from tmesim.validation.kill_calibration import literature_bridge
    solid = literature_bridge(scope="solid")
    out = {}
    for name, d in panel["drugs"].items():
        b = literature_bridge(scope=f"drug:{name}")
        spec = b if b.get("n") else None
        out[name] = _band_entry(spec, solid, d["c_bath_uM"], d["ic50_uM"],
                                d["hill"])
    return out


def breast_bands(panel: Dict) -> Dict:
    """The published breast literature bands, via the breast report's own
    band builder. Both the drug-matched band and the broad context band are
    carried; which (if either) scores a row is ``score_ratio_row``'s call."""
    from experiments.dtissue_sensitivity import build_bands
    ic50_hill = {n: (d["ic50_uM"], d["hill"])
                 for n, d in panel["drugs"].items()}
    raw = build_bands(list(panel["drugs"]), ic50_hill)
    return {name: {"band": b["band"], "band_solid": b["band_solid"],
                   "specific_bridge": b["specific_bridge"],
                   "solid_bridge": b["solid_bridge"],
                   "specific_n": b["specific_n"]}
            for name, b in raw.items()}


def _p99_score(predicted_ratio: float, b: Dict) -> str:
    """The verdict for one predicted ratio against one band record.

    A one-line wrapper over ``kill_calibration.score_ratio_row`` so that every
    verdict string in this file -- published, bias-corrected, or convention-swept
    -- comes out of the same decision. Returns the status only; callers here
    never need the rest.
    """
    from tmesim.validation.kill_calibration import score_ratio_row
    return score_ratio_row(predicted_ratio,
                           specific_bridge=b["specific_bridge"],
                           specific_band=b["band"],
                           context_bridge=b["solid_bridge"],
                           context_band=b["band_solid"])["status"]


def _scope_label(b: Dict) -> str:
    if b["specific_bridge"] is None:
        return "NONE drug-matched — not scored"
    if b["specific_n"] < 2:
        return f"single study (n = {b['specific_n']}) — zero-width, not scored"
    return f"drug-matched (n = {b['specific_n']})"


def stage_d(panels: Sequence[Tuple[Dict, Dict]]) -> Dict:
    """Recompute all nine published verdicts under every convention.

    Deliberately the LAST stage, and its output is reported in a separate
    section of the write-up, because selecting a normalising convention by how
    many PASSes it produces would be exactly the kind of rigging this repo
    exists not to do.
    """
    from tmesim.validation.kill_calibration import (kill_ratio_to_ic50_fold,
                                                    score_ratio_row)
    out: Dict = {"panels": []}
    for panel, bands in panels:
        prec: Dict = {"tissue": panel["tissue"], "drugs": []}
        for name, d in panel["drugs"].items():
            e, sens, resist = d["e_unit"], d["sens"], d["resist"]
            cb, ic, hh, mk = (d["c_bath_uM"], d["ic50_uM"], d["hill"],
                              d["max_kill"])
            assert_matches_shipped(e, sens, resist, cb, ic, hh, mk)
            fixed_ref = float(np.percentile(e, 99.0))
            b = bands[name]
            scorable = b["specific_bridge"] is not None and b["specific_n"] >= 2
            row: Dict = {"drug": name, "c_bath_label": d["c_bath_label"],
                         "band_scope": _scope_label(b),
                         "scored": bool(scorable),
                         "specific_n": b["specific_n"],
                         "band_ratio_lo": (b["band"]["ratio_lo"] if scorable
                                           else None),
                         "band_ratio_hi": (b["band"]["ratio_hi"] if scorable
                                           else None),
                         "band_fold_lo": (b["band"]["fold_lo"] if scorable
                                          else None),
                         "band_fold_hi": (b["band"]["fold_hi"] if scorable
                                          else None),
                         "context_ratio_lo": b["band_solid"]["ratio_lo"],
                         "context_ratio_hi": b["band_solid"]["ratio_hi"],
                         "conventions": {}}
            for tag in [n["tag"] for n in NORMALISERS]:
                dec = decompose_general(e, sens, resist, cb, ic, hh, mk,
                                        tag=tag, fixed_ref=fixed_ref)
                s = score_ratio_row(dec["ratio_exact"],
                                    specific_bridge=b["specific_bridge"],
                                    specific_band=b["band"],
                                    context_bridge=b["solid_bridge"],
                                    context_band=b["band_solid"])
                row["conventions"][tag] = {
                    "ratio_exact": dec["ratio_exact"],
                    "f_penetration": dec["f_penetration"],
                    "implied_ic50_fold": kill_ratio_to_ic50_fold(
                        dec["ratio_exact"], cb, ic, hh),
                    "verdict": s["status"],
                    "scored": bool(s["scored"]),
                    "direction": (s["test"]["direction"] if s["scored"]
                                  else "not scored — no testable band"),
                    # CONTEXT ONLY: what the broad non-drug-matched band would
                    # have said. Never counted as a verdict.
                    "context_verdict": s["context_test"]["verdict"],
                }
            prec["drugs"].append(row)
        out["panels"].append(prec)
    return out


def measured_bias_at_dx(dx: float = 12.0, geometry: str = "vessel",
                        law: str = "off", tag: str = SHIPPED_TAG
                        ) -> Dict[float, float]:
    """The stage-B measured grid bias of ``f_penetration``, per bath multiple.

    Read from ``results/sweep/p99_gridslab.json`` rather than quoted, so the
    verdict-risk arithmetic below cannot drift from the measurement. Returns
    ``{c_bath/IC50 -> relative error of f_penetration at that dx}``.
    """
    slab = json.load(open(os.path.join(OUT_DIR, "p99_gridslab.json")))
    pts = slab["geometries"][geometry]["points"]
    p = next(q for q in pts if abs(q["dx_um"] - dx) < 1e-9)
    err = p["laws"][law]["rel_err"][tag]
    return {float(k[5:]): float(v) for k, v in err.items()
            if k.startswith("fpen_")}


def stage_d_bias(panels: Sequence[Tuple[Dict, Dict]],
                 biases: Sequence[float] = (-0.083, -0.054, -0.022, -0.011)
                 ) -> Dict:
    """What the measured grid bias does to the nine verdicts, if corrected.

    ``results/boundary_layer_resolution.md`` sec. 4.4/6 measures the shape
    statistic ``mean/p99`` as reading **+8.35 % HIGH at dx = 12 um** against an
    exact reference (+5.4 % at dx = 8, +2.2 % at 4, +1.1 % at 2). Correcting a
    +8.3 % over-read means multiplying the scored ratio by ``1/1.083``, i.e. the
    entries below with ``bias = -0.083``. The correction is applied to the
    scored RATIO because in the low-dose limit ``f_penetration`` IS that shape
    statistic; the mapping is exact only there, which is stated as a limit
    rather than hidden.
    """
    from tmesim.validation.kill_calibration import kill_ratio_to_ic50_fold
    out: Dict = {"biases": list(biases), "panels": []}
    #: The stage-B measurement, per dx, keyed by bath multiple. Applied per drug
    #: at ITS OWN reporting concentration, because the Hill map damps the bias:
    #: the +8.3 % quoted repo-wide is the low-dose / linear limit.
    try:
        measured = {dx: measured_bias_at_dx(dx) for dx in (12.0, 8.0, 4.0, 2.0)}
    except Exception as exc:                       # stage B not run yet
        log(f"stage D: measured grid bias unavailable ({exc}); "
            f"only the flat corrections are reported")
        measured = {}
    out["measured_fpen_bias_by_dx_and_bath_multiple"] = {
        str(k): v for k, v in measured.items()}
    for panel, bands in panels:
        prec: Dict = {"tissue": panel["tissue"], "drugs": []}
        for name, d in panel["drugs"].items():
            e, sens, resist = d["e_unit"], d["sens"], d["resist"]
            cb, ic, hh, mk = (d["c_bath_uM"], d["ic50_uM"], d["hill"],
                              d["max_kill"])
            b = bands[name]
            base = decompose_general(e, sens, resist, cb, ic, hh, mk,
                                     tag=SHIPPED_TAG)["ratio_exact"]
            scorable = b["specific_bridge"] is not None and b["specific_n"] >= 2
            band_lo = b["band"]["ratio_lo"] if scorable else None
            band_hi = b["band"]["ratio_hi"] if scorable else None

            def _score(x):
                return _p99_score(x, b)

            row: Dict = {"drug": name, "ratio_published": base,
                         "verdict_published": _score(base),
                         "scored": bool(scorable),
                         "specific_n": b["specific_n"],
                         "band_scope": _scope_label(b),
                         "band_ratio_lo": band_lo,
                         "band_ratio_hi": band_hi,
                         # A margin is a distance to the bar this drug is held
                         # to. An unscored drug has no bar, so it has no margin;
                         # quoting its distance to the broad context band as one
                         # is exactly the substitution this file used to make.
                         "margin_to_lo_pct": ((base / band_lo - 1.0) * 100
                                              if scorable and band_lo > 0
                                              else None),
                         "margin_to_hi_pct": ((base / band_hi - 1.0) * 100
                                              if scorable and band_hi > 0
                                              else None),
                         "corrected": {}}
            for eps in biases:
                r = base * (1.0 + eps)
                v = _score(r)
                row["corrected"][f"{eps:+.3f}"] = {
                    "ratio": r, "verdict": v,
                    "implied_ic50_fold": kill_ratio_to_ic50_fold(r, cb, ic, hh),
                    "changed": v != row["verdict_published"],
                }
            # --- the MEASURED, concentration-matched correction -------------
            m = float(cb) / float(ic)
            row["bath_multiple"] = m
            row["measured_correction"] = {}
            for dx, table in measured.items():
                mm = min(table, key=lambda k: abs(math.log(k / m)))
                bias = table[mm]
                r = base / (1.0 + bias)      # undo an over-read of size `bias`
                v = _score(r)
                row["measured_correction"][f"dx{dx:g}"] = {
                    "nearest_bath_multiple_measured": mm,
                    "measured_fpen_bias": bias,
                    "ratio_corrected": r, "verdict": v,
                    "implied_ic50_fold": kill_ratio_to_ic50_fold(r, cb, ic, hh),
                    "changed": v != row["verdict_published"],
                }
            prec["drugs"].append(row)
        out["panels"].append(prec)
    return out


# --------------------------------------------------------------------------- #
# STAGE E -- estimator stability
# --------------------------------------------------------------------------- #
BOOT_N = 400
BOOT_SEED = 20260809


def stage_e(panels: Sequence[Dict], boot_n: int = BOOT_N) -> Dict:
    """Is the normaliser a stable estimator on the cell population it reads?

    Two independent probes, neither of which involves the grid:

    ``bootstrap``  resample the tumour cells with replacement. This is the
        sampling variability the normaliser inherits from the fact that the
        cell set is one finite realisation of a tissue.
    ``halfsample`` score on a random HALF of the cells. An order statistic deep
        in the tail is estimated from ``(1-q)*n`` cells, so halving n halves
        that support; the mean does not care.
    """
    rng = np.random.default_rng(BOOT_SEED)
    out: Dict = {"boot_n": int(boot_n), "panels": []}
    for panel in panels:
        prec: Dict = {"tissue": panel["tissue"], "n_tumor": panel["n_tumor"],
                      "drugs": []}
        for name, d in panel["drugs"].items():
            e, sens, resist = d["e_unit"], d["sens"], d["resist"]
            cb, ic, hh, mk = (d["c_bath_uM"], d["ic50_uM"], d["hill"],
                              d["max_kill"])
            n = e.size
            fixed_ref = float(np.percentile(e, 99.0))
            idx_boot = [rng.integers(0, n, size=n) for _ in range(boot_n)]
            idx_half = [rng.permutation(n)[: n // 2] for _ in range(boot_n)]
            row: Dict = {"drug": name, "n_cells": int(n),
                         "n_cells_above_p99": int((e > np.percentile(e, 99.0)).sum()),
                         "n_cells_above_p999": int((e > np.percentile(e, 99.9)).sum()),
                         "conventions": {}}
            for tag in [nm["tag"] for nm in NORMALISERS]:
                base = decompose_general(e, sens, resist, cb, ic, hh, mk,
                                         tag=tag, fixed_ref=fixed_ref
                                         )["ratio_exact"]
                for label, idxs in (("bootstrap", idx_boot),
                                    ("halfsample", idx_half)):
                    vals = np.array([
                        decompose_general(e[i], sens[i], resist[i], cb, ic, hh,
                                          mk, tag=tag, fixed_ref=fixed_ref
                                          )["ratio_exact"] for i in idxs])
                    row["conventions"].setdefault(tag, {})[label] = {
                        "mean": float(vals.mean()),
                        "cv_pct": float(100.0 * vals.std(ddof=1)
                                        / max(abs(vals.mean()), 1e-300)),
                        "bias_vs_full_pct": float(
                            100.0 * (vals.mean() / base - 1.0)),
                    }
                row["conventions"][tag]["ratio_full"] = float(base)
            prec["drugs"].append(row)
            log(f"stage E: {panel['tissue']}/{name} done")
        out["panels"].append(prec)
    return out


# --------------------------------------------------------------------------- #
# STAGE F -- what is p99 actually standing in for?
# --------------------------------------------------------------------------- #
def stage_f(panels: Sequence[Dict]) -> Dict:
    """Locate the PHYSICAL bath analogue on the model's own exposure distribution.

    The convention's stated justification is that "the best-perfused cells in
    tissue are the analogue" of a spheroid's outer surface. That is a structural
    claim, and it is checkable: the structural analogue of the outer shell is
    the set of tumour cells sitting in a PERFUSED voxel (the shipped vascular
    field is binary, so 'perfused' is exactly defined). This computes where that
    set's mean exposure falls on the whole distribution's percentile scale, and
    compares it with where each candidate convention falls.

    It also prices the SPHEROID side of the same analogy from geometry alone:
    the fraction of a spheroid's cells lying within one cell diameter of its
    surface is ``1 - (1 - d/R)^3``, which is what "the outer surface sits at the
    bath" refers to in the papers the bands come from.
    """
    out: Dict = {"spheroid_shell_fractions": [], "panels": []}
    for r_um in (50.0, 100.0, 150.0, 250.0):
        for d_um in (10.0, 15.0, 20.0):
            out["spheroid_shell_fractions"].append({
                "spheroid_radius_um": r_um, "cell_diameter_um": d_um,
                "outer_shell_volume_fraction":
                    float(1.0 - max(0.0, 1.0 - d_um / r_um) ** 3),
            })
    for panel in panels:
        work, labeled = panel.get("_work"), panel.get("_labeled")
        tumor = panel.get("_tumor")
        prec: Dict = {"tissue": panel["tissue"], "drugs": []}
        if work is None:
            out["panels"].append(prec)
            continue
        vasc = np.asarray(work.fields.vascular_density, dtype=np.float64)
        col, row = work.domain.to_index(labeled.x, labeled.y)
        perfused_cell = (vasc[row, col] > 0.0)[tumor]
        prec["n_tumor"] = int(tumor.sum())
        prec["n_tumor_in_perfused_voxel"] = int(perfused_cell.sum())
        prec["perfused_tumour_cell_fraction"] = float(
            perfused_cell.sum() / max(1, tumor.sum()))
        for name, d in panel["drugs"].items():
            e = d["e_unit"]
            row_out: Dict = {"drug": name}
            if perfused_cell.any():
                ref = float(e[perfused_cell].mean())
                row_out["perivascular_mean_exposure"] = ref
                row_out["its_percentile_on_the_full_distribution"] = float(
                    100.0 * (e < ref).mean())
            for n in NORMALISERS:
                if not n["homogeneous"]:
                    continue
                v = float(n["fn"](e))
                row_out[f"pct_of_{n['tag']}"] = float(100.0 * (e < v).mean())
            prec["drugs"].append(row_out)
        out["panels"].append(prec)
    return out


# --------------------------------------------------------------------------- #
# STAGE G -- does the slab's grid error transfer to tissue?
# --------------------------------------------------------------------------- #
def stage_g() -> Dict:
    """Compare the slab's predicted dx dependence with the MEASURED tissue one.

    ``results/boundary_layer_resolution.md`` sec. 8 states that the +8.3 % applies
    to "any p99-normalised or distribution-shaped quantity -- which includes every
    3-D:2-D verdict in this repo". That is a transferable claim and it is
    testable, because the slab and both tissue ladders all solve dx = 12 AND
    dx = 8. The commensurable quantity is therefore the SAME dx step on both:
    ``ratio(12)/ratio(8) - 1``.

    The tissue numbers are PARSED from the two published ladder reports (owned by
    sibling workstreams; nothing is recomputed here). The slab number comes from
    this file's own stage B.
    """
    import re

    slab = json.load(open(os.path.join(OUT_DIR, "p99_gridslab.json")))
    pts = {p["dx_um"]: p for p in slab["geometries"]["vessel"]["points"]}

    def slab_step(dx_hi: float, dx_lo: float, m: float, tag: str) -> float:
        """Predicted ``ratio(dx_hi)/ratio(dx_lo) - 1`` from the two grid errors."""
        key = f"fpen_{m:g}"
        a = pts[dx_hi]["laws"]["off"]["rel_err"][tag][key]
        b = pts[dx_lo]["laws"]["off"]["rel_err"][tag][key]
        return (1.0 + a) / (1.0 + b) - 1.0

    out: Dict = {"comparison": "ratio(dx=12) / ratio(dx=8) - 1",
                 "slab_geometry": "vessel (a = 7.5 um), exchange OFF",
                 "rows": []}

    # --- GBM: parse the per-drug ladder tables --------------------------------
    gbm_md = open(os.path.join(RESULTS_DIR_MD, "grid_convergence_gbm.md")).read()
    for drug in GBM_DRUGS:
        m = re.search(rf"### {re.escape(drug)} —.*?\n(.*?)\nStep-by-step",
                      gbm_md, re.S)
        if m is None:
            log(f"stage G: no ladder table for gbm/{drug} -- skipped")
            continue
        vals = {}
        for line in m.group(1).split("\n"):
            if not line.startswith("| ") or "---" in line or "dx" in line:
                continue
            c = [x.strip() for x in line.strip("|").split("|")]
            vals[float(c[0].split()[0])] = float(c[3].replace("*", ""))
        mult = 0.165 if drug == "temozolomide" else 1.0
        out["rows"].append(_transfer_row("gbm", drug, vals[12.0], vals[8.0],
                                         mult, slab_step(12.0, 8.0, mult, "p99")))

    # --- breast: parse the summary table -------------------------------------
    br_md = open(os.path.join(RESULTS_DIR_MD,
                              "grid_convergence_breast.md")).read()
    m = re.search(r"\| drug \| ratio @ dx=12 \| ratio @ finest \|.*?\n\|[-:| ]+\|\n(.*?)\n\n",
                  br_md, re.S)
    if m is None:
        log("stage G: breast summary table not found -- breast rows skipped")
    else:
        for line in m.group(1).split("\n"):
            if not line.startswith("| "):
                continue
            c = [x.strip() for x in line.strip("|").split("|")]
            out["rows"].append(_transfer_row(
                "breast", c[0], float(c[1]), float(c[2]), 1.0,
                slab_step(12.0, 8.0, 1.0, "p99")))

    tis = [r["tissue_step"] for r in out["rows"]]
    out["summary"] = {
        "n": len(tis),
        "tissue_step_min": min(tis), "tissue_step_max": max(tis),
        "tissue_spread_pct_points": 100.0 * (max(tis) - min(tis)),
        "slab_step_at_1xIC50": slab_step(12.0, 8.0, 1.0, "p99"),
        "n_sign_agree": sum(1 for r in out["rows"] if r["sign_agrees"]),
        "n_sign_opposite": sum(1 for r in out["rows"] if not r["sign_agrees"]),
        "opposite_drugs": [r["drug"] for r in out["rows"]
                           if not r["sign_agrees"]],
        "verdict": (
            "The slab error is drug-INDEPENDENT by construction; the tissue error "
            "is drug-dependent and not one-signed. The slab figure is a correct "
            "measurement of the scoring statistic's error in isolation and is "
            "valid for RANKING conventions, but it does not transfer as a bias "
            "correction to a tissue verdict."),
    }
    log(f"stage G: tissue step spans {100*min(tis):+.2f}% to {100*max(tis):+.2f}%, "
        f"slab predicts {100*out['summary']['slab_step_at_1xIC50']:+.2f}%; "
        f"{out['summary']['n_sign_opposite']} of {len(tis)} disagree in sign "
        f"({', '.join(out['summary']['opposite_drugs']) or 'none'})")
    return out


def _transfer_row(tissue: str, drug: str, r12: float, r8: float, mult: float,
                  slab: float) -> Dict:
    step = r12 / r8 - 1.0
    return {"tissue": tissue, "drug": drug.strip(), "bath_multiple": mult,
            "ratio_dx12": r12, "ratio_dx8": r8, "tissue_step": step,
            "slab_step_predicted": slab,
            "sign_agrees": bool(step * slab > 0.0),
            "slab_over_tissue": (slab / step) if step != 0 else float("inf")}


RESULTS_DIR_MD = os.path.join(_ROOT, "results")


# --------------------------------------------------------------------------- #
def _json_safe(obj):
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items()
                if not (isinstance(k, str) and k.startswith("_"))}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return obj


def _dump(obj, name: str) -> None:
    os.makedirs(OUT_DIR, exist_ok=True)
    path = os.path.join(OUT_DIR, name)
    with open(path, "w") as fh:
        json.dump(_json_safe(obj), fh, indent=2)
    log(f"written -> results/sweep/{name}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="A",
                    choices=["A", "B", "C", "D", "E", "F", "G", "ADEF"])
    ap.add_argument("--boot", type=int, default=BOOT_N)
    args = ap.parse_args()

    if args.stage == "B":
        _dump(stage_b(), "p99_gridslab.json")
        return

    if args.stage == "G":
        _dump(stage_g(), "p99_transfer.json")
        return

    if args.stage == "C":
        _dump(stage_c(checkpoint=os.path.join(OUT_DIR,
                                              "p99_dtissue_dx12.jsonl")),
              "p99_dtissue_dx12.json")
        return

    gbm = load_gbm_panel()
    log(f"GBM panel loaded: {gbm['n_tumor']} tumour cells")
    breast = load_breast_panel()
    log(f"breast panel loaded: {breast['n_tumor']} tumour cells")
    panels = [gbm, breast]

    if args.stage in ("A", "ADEF"):
        _dump(stage_a(panels), "p99_algebra.json")
    if args.stage in ("D", "ADEF"):
        pairs = [(gbm, gbm_bands(gbm)), (breast, breast_bands(breast))]
        _dump({"verdicts": stage_d(pairs), "bias": stage_d_bias(pairs)},
              "p99_verdicts.json")
    if args.stage in ("F", "ADEF"):
        _dump(stage_f(panels), "p99_perivascular.json")
    if args.stage in ("E", "ADEF"):
        _dump(stage_e(panels, args.boot), "p99_stability.json")


if __name__ == "__main__":
    main()
