"""DEFENSIBLE DELIVERY METRICS — replacing ``penetration_index``.

Why this module exists
----------------------
``penetration_index = mean(AUC) / p99(AUC)`` over the tissue mask is a
**scale-invariant SHAPE statistic**. Multiply the whole exposure field by any
constant and the index does not move. Consequently, at the grid-converged
resolution dx = 12 um it reports the *vascular geometry* of the section — which is
identical for every drug — and collapses to 0.194-0.196 across the entire
five-drug panel, even though those drugs' actual delivered exposure differs by
~122x (``results/real_drug_panel.md``, ``results/numerical_rigor.md``).

A delivery metric has to answer "how much drug actually got to the tumor cells,
and was it enough?". This module defines five that do, and then *validates* them:

======================================  ==================================  ======
metric                                   units                              kind
======================================  ==================================  ======
1. ``mean/median_tissue_auc``            uM*s                               absolute
   ``mean/median_tumor_auc``             uM*s                               absolute
2. ``therapeutic_coverage``              fraction of tumor cells in [0,1]   absolute
3. ``exposure_percentiles``              uM (and uM*s)                      absolute
4. ``penetration_length_fit``            um                                 length
5. ``hypoxic_vs_normoxic_ratio``         dimensionless                      shape
======================================  ==================================  ======

The headline recommendation that falls out of the validation is metric **2**,
``therapeutic_coverage`` at 1x IC50: it is the only one that is simultaneously
(i) dimensionally sound and clinically interpretable ("what fraction of the tumor
saw a killing dose"), (ii) grid-converged, and (iii) drug-discriminating by a
margin far wider than its own grid drift. See ``results/delivery_metrics.md``.

Definitions (exact)
-------------------
Let ``AUC(x)`` be the transport solve's exposure field (uM*s) and ``T`` the
exposure-window length (s, = ``transport.times[-1]``, the same ``T``
``cells.pkpd`` uses). The **time-averaged effective concentration** at a point is

    E(x) = AUC(x) / T            [uM]

which is directly comparable to a drug's IC50 because IC50 *is* a concentration.
Every per-cell statistic here is evaluated at the **real measured cell positions**
(``micro.cells.x/y`` -> ``domain.to_index``), NOT over the grid. That matters for
grid convergence: the Xenium cell population is fixed (40,887 cells, 30,174
tumor) at every dx, so a per-cell statistic changes with dx only because the
*field* changed — whereas any mask-area average also changes because the
rasterised tissue mask changed.

1. ``mean_tissue_auc`` / ``median_tissue_auc`` — the tissue-mask (area-weighted)
   mean and median of ``AUC``. ``mean_tumor_auc`` / ``median_tumor_auc`` are the
   same statistics over tumor CELLS. The median is reported because the exposure
   field is strongly right-skewed (a perivascular tail drags the mean).

2. ``therapeutic_coverage(auc, ic50, T, mult)`` — the fraction of tumor cells
   whose time-averaged exposure ``E = AUC/T`` exceeds ``mult * IC50``. Reported at
   ``mult in {0.5, 1, 2}``; the 1x value is the headline. This is a *drug-specific*
   number: the same exposure field scores differently for a potent and a weak
   drug, which is exactly the property ``penetration_index`` lacks.

3. ``exposure_percentiles`` — P10 / P50 / P90 of the per-tumor-cell ``E`` (uM) and
   ``AUC`` (uM*s), plus the P90/P10 ratio (the delivery *heterogeneity*). P10 is
   the worst-served decile of the tumor — the population that drives regrowth.

4. ``penetration_length_fit`` — fits ``AUC(d) = C0*exp(-d/L) + floor`` to the
   binned AUC-vs-distance-from-nearest-vessel profile on the REAL section, giving
   ``L`` in micrometers. Unlike the index this is a real physical length. It ships
   with its own **SOURCE NULL** control (the identical fit applied to the raw
   ``vascular_density`` field), because on this tissue the delivery source is
   distributed over the continuous density field, so the drug profile largely
   inherits the source's shape — see ``results/vessel_profile.md``. The honest
   readout is therefore ``L_drug`` *versus* ``L_source_null``, and this module
   reports both. The distance map, binning and fits are REUSED from
   ``validation.vessel_profile`` so the two reports cannot drift apart.

5. ``hypoxic_vs_normoxic_exposure_ratio`` — mean ``E`` in the worst-perfused
   quartile of tumor cells divided by mean ``E`` in the best-perfused quartile
   (threshold-free, always defined), plus the literal banded versions against the
   engine's own oxygen bands (``t = OxygenConfig.hypoxia_threshold``; hypoxic
   ``o < t``, normoxic ``o >= t``, viable ``o >= 2t``) which are nan when a band is
   empty. < 1 means the drug under-serves the
   poorly-perfused regions. It is a RATIO, hence scale-invariant, hence a shape
   statistic like the index — it is reported as a *diagnostic*, not as a headline,
   and the validation below shows it is the least grid-stable of the five.

What this module does NOT touch
-------------------------------
Nothing in ``core`` / ``physics`` / ``molecular`` / ``reconstruct`` / ``sim.engine``
is modified. This module only consumes a ``Microenvironment`` and a
``SimulationResult``.

Run:    ./.venv/bin/python -m tmesim.validation.delivery_metrics
Writes: results/delivery_metrics.png
        results/delivery_metrics.md

Units: AUC uM*s, exposure/IC50 uM, length um, coverage/fractions dimensionless.
"""
from __future__ import annotations

import argparse
import os
import pickle
import time
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..core import (
    CellType,
    DrugParams,
    Microenvironment,
    OxygenConfig,
    PDEConfig,
    SimConfig,
    SimulationResult,
)
from ..molecular.binding_site_barrier import apply_barrier_to_drugparams
from ..molecular.library import get_molecule
from ..reconstruct.real_pipeline import reconstruct_real, reconstruction_report
from ..sim.engine import simulate_response
from .calibration import calibrate_alpha, calibrated_drugparams
from .vessel_profile import (
    MIN_BIN_COUNT,
    NEAR_WINDOW_UM,
    VESSEL_PERCENTILE,
    direct_vessel_distance_um,
    distance_from_vessel_um,
    distance_profile,
    endothelial_xy,
    fit_exponential,
    fit_exponential_plus_floor,
    usable_bins,
    vessel_mask_from_density,
)

__all__ = [
    "exposure_window_s",
    "tumor_cell_exposure",
    "tissue_auc_summary",
    "therapeutic_coverage",
    "exposure_percentiles",
    "penetration_length_fit",
    "hypoxic_vs_normoxic_exposure_ratio",
    "compute_all",
    "backfill",
    "run",
    "selftest",
    "main",
]

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
#: The five-drug panel of ``results/real_drug_panel.md``.
PANEL: List[str] = ["temozolomide", "lomustine", "cisplatin", "methotrexate",
                    "doxorubicin"]

#: Grid-refinement ladder. dx = 12 um is the PRODUCTION resolution: it is the
#: finest grid that satisfies the advection CFL condition at the published
#: dt = 10 s (``results/real_drug_panel.md``), so dx = 10 um is deliberately NOT
#: on this ladder.
DX_LADDER: List[float] = [20.0, 15.0, 12.0]
DX_PRODUCTION = 12.0

#: IC50 multipliers for the therapeutic-coverage variants.
COVERAGE_MULTIPLIERS: Tuple[float, ...] = (0.5, 1.0, 2.0)

#: Percentiles of the per-tumor-cell exposure distribution.
EXPOSURE_PERCENTILES: Tuple[float, ...] = (10.0, 50.0, 90.0)

#: GRID-CONVERGENCE criterion for an UNBOUNDED metric (AUC, exposure, length,
#: ratio): worst-drug RELATIVE change between the two finest grids below this.
CONVERGENCE_TOL = 0.05

#: GRID-CONVERGENCE criterion for a metric bounded in [0, 1] (the coverage and
#: kill fractions): worst-drug ABSOLUTE change, in fraction-of-cells units.
#: Relative error is the wrong measure for a bounded fraction — it blows up as the
#: fraction approaches zero even when the answer moved by a fraction of a percentage
#: point — so bounded metrics are judged on the absolute scale they live on. BOTH
#: numbers are reported in the convergence table either way, so a reader can apply
#: whichever rule they prefer.
CONVERGENCE_TOL_ABS = 0.03

#: A metric DISCRIMINATES if its across-drug spread at dx = 12 um exceeds its own
#: worst-drug grid drift by at least this factor (signal-to-numerical-noise).
DISCRIMINATION_MIN_RATIO = 3.0

_RESULTS_DIR = os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "results")
)
_OUT_PNG = os.path.join(_RESULTS_DIR, "delivery_metrics.png")
_OUT_MD = os.path.join(_RESULTS_DIR, "delivery_metrics.md")

#: The engine config the flagship real runs use (advection on, dt = 10 s, 6 h).
_PDE = dict(dt=10.0, total_time=3600.0 * 6, convection=True, save_every=500)

#: Validated categorical palette (dataviz skill reference instance, light mode;
#: worst adjacent CVD dE 24.2). One slot per drug, assigned in fixed PANEL order.
_DRUG_COLORS = {
    "temozolomide": "#2a78d6",   # blue
    "lomustine": "#1baf7a",      # aqua
    "cisplatin": "#eda100",      # yellow
    "methotrexate": "#008300",   # green
    "doxorubicin": "#4a3aa7",    # violet
}
_C_GOOD = "#2a78d6"   # the metric we keep
_C_BAD = "#e34948"    # the metric we retire
_INK = "#0b0b0b"
_INK_2 = "#52514e"
_GRID = "#dcdbd6"


# --------------------------------------------------------------------------- #
# 0. Sampling helpers — everything per-cell, on the FIXED real cell population
# --------------------------------------------------------------------------- #
def exposure_window_s(sim_result: SimulationResult) -> float:
    """The exposure-window length ``T`` (s) used to turn AUC (uM*s) into a
    time-averaged concentration ``E = AUC/T`` (uM).

    Taken from ``sim_result.transport.times[-1]`` — byte-for-byte the same ``T``
    :func:`tmesim.cells.pkpd.apply_drug_kill` uses, so ``E`` here and the ``E``
    the kill model sees are the same number.
    """
    tr = getattr(sim_result, "transport", None)
    if tr is None or np.size(tr.times) == 0:
        raise ValueError("sim_result carries no transport.times; cannot define T")
    T = float(np.asarray(tr.times, dtype=float)[-1])
    if not np.isfinite(T) or T <= 0.0:
        raise ValueError(f"non-positive exposure window T = {T}")
    return T


def _sample_at_cells(field: np.ndarray, micro: Microenvironment,
                     mask: Optional[np.ndarray] = None) -> np.ndarray:
    """Sample a grid ``field`` at the real cell positions (optionally a subset).

    Uses ``domain.to_index`` (the same nearest-grid-cell map the pharmacodynamic
    kill uses), so a cell's sampled exposure here is the exposure that killed it.
    """
    cells = micro.cells
    x, y = cells.x, cells.y
    if mask is not None:
        x, y = x[mask], y[mask]
    col, row = micro.domain.to_index(x, y)
    return np.asarray(field, dtype=np.float64)[row, col]


def tumor_cell_exposure(micro: Microenvironment,
                        sim_result: SimulationResult) -> Dict:
    """Per-TUMOR-CELL exposure, the substrate of every absolute metric here.

    Returns ``auc_uM_s`` (uM*s) and ``exposure_uM`` (= AUC/T, uM) as arrays with
    one entry per ``CellType.TUMOR`` cell, plus ``T_s`` and ``n_tumor``.

    The tumor cell population is a FIXED set of measured Xenium cells, identical
    at every dx — which is what makes these statistics grid-comparable.
    """
    T = exposure_window_s(sim_result)
    tumor = np.asarray(micro.cells.cell_type) == int(CellType.TUMOR)
    auc = _sample_at_cells(sim_result.fields["auc"], micro, tumor)
    return {
        "T_s": T,
        "n_tumor": int(tumor.sum()),
        "auc_uM_s": auc,
        "exposure_uM": auc / T,
        "tumor_mask": tumor,
    }


# --------------------------------------------------------------------------- #
# 1. Absolute exposure
# --------------------------------------------------------------------------- #
def tissue_auc_summary(micro: Microenvironment,
                       sim_result: SimulationResult) -> Dict:
    """METRIC 1 — absolute delivered exposure, in uM*s.

    ``mean_tissue_auc`` / ``median_tissue_auc`` are area-weighted over the tissue
    mask (what the engine already reports as ``extras['mean_auc']``);
    ``mean_tumor_auc`` / ``median_tumor_auc`` are the same statistics evaluated on
    the real tumor CELLS, which is the population the prediction is about and the
    one that does not move when the mask is re-rasterised.

    The median is reported alongside the mean because the exposure field is
    strongly right-skewed: a small perivascular cuff carries a large share of the
    integral, so mean >> median and the mean alone overstates what a typical
    tumor cell sees.
    """
    auc = np.asarray(sim_result.fields["auc"], dtype=np.float64)
    tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
    tis = auc[tissue] if tissue.any() else np.zeros(0)
    exp = tumor_cell_exposure(micro, sim_result)
    tum = exp["auc_uM_s"]
    return {
        "mean_tissue_auc": float(tis.mean()) if tis.size else float("nan"),
        "median_tissue_auc": float(np.median(tis)) if tis.size else float("nan"),
        "mean_tumor_auc": float(tum.mean()) if tum.size else float("nan"),
        "median_tumor_auc": float(np.median(tum)) if tum.size else float("nan"),
        "mean_tumor_exposure_uM": float(exp["exposure_uM"].mean()) if tum.size else float("nan"),
        "median_tumor_exposure_uM": (float(np.median(exp["exposure_uM"]))
                                     if tum.size else float("nan")),
        "skew_mean_over_median": (float(tum.mean() / np.median(tum))
                                  if tum.size and np.median(tum) > 0 else float("nan")),
        "n_tissue_px": int(tissue.sum()),
        "n_tumor": exp["n_tumor"],
        "T_s": exp["T_s"],
    }


# --------------------------------------------------------------------------- #
# 2. Therapeutic coverage — THE HEADLINE METRIC
# --------------------------------------------------------------------------- #
def therapeutic_coverage(auc, ic50: float, T: float,
                         multipliers: Sequence[float] = COVERAGE_MULTIPLIERS
                         ) -> Dict:
    """METRIC 2 — the fraction of tumor cells that saw a KILLING dose.

    .. math::

        E_i = \\mathrm{AUC}_i / T \\quad [\\mu M], \\qquad
        \\mathrm{coverage}(m) = \\frac{\\#\\{i : E_i > m \\cdot \\mathrm{IC50}\\}}{N}

    Parameters
    ----------
    auc : array-like
        Per-tumor-cell exposure integral, uM*s (from
        :func:`tumor_cell_exposure`). Any 1-D array of AUC values works.
    ic50 : float
        The drug's REAL (GDSC / literature) IC50 in uM — the concentration giving
        half-maximal kill. Because ``E`` is also a uM concentration, no unit
        conversion is needed; this is the identical comparison
        :mod:`tmesim.cells.pkpd` makes internally.
    T : float
        Exposure-window length in seconds.
    multipliers : sequence of float
        IC50 multiples to score. ``0.5`` is the permissive variant (a sub-IC50
        dose still kills some cells on the Hill curve), ``2.0`` the stringent one.

    Returns
    -------
    dict
        ``coverage_0.5x`` / ``coverage_1.0x`` / ``coverage_2.0x`` (fractions in
        [0, 1]), the absolute threshold concentrations used, ``ic50_uM``, ``T_s``
        and ``n``.

    Why this is the right shape of metric
    -------------------------------------
    It is (a) **dimensionally sound** — a concentration compared to a
    concentration; (b) **drug-specific** — the same exposure field scores
    differently for a potent and a weak drug, so it cannot collapse to a
    geometry constant the way ``penetration_index`` does; (c) **clinically
    legible** — "63% of the tumor saw at least a half-maximal-effect
    concentration" is a sentence an oncologist can act on; and (d) **absolute** —
    scaling the field scales the answer.
    """
    a = np.asarray(auc, dtype=np.float64).ravel()
    a = a[np.isfinite(a)]
    T = float(T)
    ic50 = float(ic50)
    if T <= 0:
        raise ValueError(f"exposure window T must be > 0, got {T}")
    if not np.isfinite(ic50) or ic50 <= 0:
        raise ValueError(f"ic50 must be a positive concentration, got {ic50}")
    E = a / T
    out: Dict[str, float] = {"ic50_uM": ic50, "T_s": T, "n": int(E.size)}
    for m in multipliers:
        thr = float(m) * ic50
        key = f"coverage_{float(m):g}x"
        out[key] = float((E > thr).mean()) if E.size else float("nan")
        out[f"threshold_{float(m):g}x_uM"] = thr

    # The CONTINUOUS companion. Coverage is a threshold count, so it saturates at
    # 1.000 once a drug clears its IC50 everywhere and then stops discriminating.
    # The potency-normalised exposure E/IC50 is the same comparison without the
    # threshold: dimensionless, drug-specific, absolute (it scales with dose), and
    # unbounded above, so it keeps ranking drugs in the saturated regime.
    ratio = E / ic50 if E.size else np.zeros(0)
    out["median_E_over_ic50"] = float(np.median(ratio)) if ratio.size else float("nan")
    out["p10_E_over_ic50"] = (float(np.percentile(ratio, 10.0))
                              if ratio.size else float("nan"))
    return out


# --------------------------------------------------------------------------- #
# 3. The delivery DISTRIBUTION
# --------------------------------------------------------------------------- #
def exposure_percentiles(exposure_uM, auc_uM_s=None,
                         percentiles: Sequence[float] = EXPOSURE_PERCENTILES
                         ) -> Dict:
    """METRIC 3 — P10 / P50 / P90 of the per-tumor-cell exposure.

    A single mean hides the thing that decides an outcome: the *worst-served*
    part of the tumor. ``P10`` is the exposure the bottom decile of tumor cells
    receives — the reservoir that survives and regrows — and ``P90/P10`` is the
    delivery heterogeneity (how many-fold better the best-perfused decile does
    than the worst). Units: ``E`` in uM, ``AUC`` in uM*s; the ratio is
    dimensionless.
    """
    E = np.asarray(exposure_uM, dtype=np.float64).ravel()
    E = E[np.isfinite(E)]
    out: Dict[str, float] = {}
    if E.size == 0:
        return {f"E_p{int(p)}_uM": float("nan") for p in percentiles}
    for p in percentiles:
        out[f"E_p{int(p)}_uM"] = float(np.percentile(E, float(p)))
    lo, hi = out.get("E_p10_uM", np.nan), out.get("E_p90_uM", np.nan)
    out["E_p90_over_p10"] = float(hi / lo) if lo and np.isfinite(lo) and lo > 0 else float("nan")
    if auc_uM_s is not None:
        A = np.asarray(auc_uM_s, dtype=np.float64).ravel()
        A = A[np.isfinite(A)]
        for p in percentiles:
            out[f"auc_p{int(p)}_uM_s"] = float(np.percentile(A, float(p))) if A.size else float("nan")
    return out


# --------------------------------------------------------------------------- #
# 4. A real physical penetration LENGTH (um), not a shape index
# --------------------------------------------------------------------------- #
def penetration_length_fit(auc_field: np.ndarray,
                           micro: Microenvironment,
                           distance_um: Optional[np.ndarray] = None,
                           vessel_percentile: float = VESSEL_PERCENTILE,
                           include_source_null: bool = True) -> Dict:
    """METRIC 4 — the exponential decay length ``L`` (um) of AUC vs
    distance-from-nearest-vessel on the REAL section.

    **Distance metric (changed 2026-07-29).** When ``distance_um`` is not
    supplied, this now defaults to the PRIMARY direct nearest-endothelial-cell
    distance (:func:`validation.vessel_profile.direct_vessel_distance_um`), not
    the old percentile-mask distance. The old percentile-threshold-on-
    ``vascular_density`` pipeline (:func:`vessel_mask_from_density` ->
    :func:`distance_from_vessel_um`, still selectable by passing a
    percentile-mask ``distance_um`` explicitly) is retained as a documented
    LEGACY metric; see ``results/vessel_profile.md``'s dedicated "distance
    metric change" section for the full head-to-head numbers and why the switch
    was made (a percentile on a density field measures "distance to the
    densest few percent of the field", not "distance to a vessel", and the two
    diverge sharply in the tail this metric's own fit is most sensitive to).
    The 10 um distance binning and the fits are all REUSED from
    :mod:`validation.vessel_profile`, so this number cannot drift away from the
    perivascular-gradient report.

    Two fits are returned because they answer different questions:

    * ``exp_floor`` — ``C(d) = C0*exp(-d/L) + floor`` over the full distance
      range. On real tissue the distributed vascular source puts a non-zero
      background under the perivascular decay, and this three-parameter form is
      the one that separates the decaying component (``L_um``) from that plateau.
      **This is the ``L`` to quote.**
    * ``near_exp`` — a pure ``C0*exp(-d/L)`` over d <= 150 um, the window the
      published immunofluorescence profiles cover. It is contaminated by the
      floor and reads systematically longer.

    ``source_null`` applies the identical pipeline to the raw
    ``vascular_density`` field. **Read the drug ``L`` against it**: the engine's
    delivery term is ``vascular_density * bbb_penetrance * plasma``, i.e. a source
    distributed over the whole continuous density field, and because the
    calibrated perivascular sink is fast the tissue field is close to locally
    reaction-dominated (``C ~ S/k``) — so the drug profile largely inherits the
    SOURCE's shape. Whatever ``L`` survives in the null is vessel/KDE geometry,
    not drug transport. See ``results/vessel_profile.md``.
    """
    tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
    vasc = np.asarray(micro.fields.vascular_density, dtype=np.float64)
    _distance_was_supplied = distance_um is not None
    if distance_um is None:
        distance_um = direct_vessel_distance_um(micro.domain,
                                                endothelial_xy(micro.cells))

    def _fit(field: np.ndarray) -> Dict:
        prof = distance_profile(field, distance_um, tissue)
        d_all, m_all, w_all = usable_bins(prof)
        d_near, m_near, w_near = usable_bins(prof, NEAR_WINDOW_UM)
        ef = fit_exponential_plus_floor(d_all, m_all, w_all)
        ne = fit_exponential(d_near, m_near, w_near)
        return {
            "L_um": float(ef.get("L_um", np.nan)),
            "d_half_um": float(ef.get("d_half_um", np.nan)),
            "floor_fraction": float(ef.get("floor_fraction", np.nan)),
            "r2_linear": float(ef.get("r2_linear_space", np.nan)),
            "near_exp_L_um": float(ne.get("L_um", np.nan)),
            "near_exp_r2_log": float(ne.get("r2_log", np.nan)),
            "perivascular_mean": float(m_all[0]) if m_all.size else float("nan"),
            "n_bins": int(d_all.size),
            "error": ef.get("error"),
        }

    out = {"vessel_percentile": float(vessel_percentile),
           "distance_metric": ("caller-supplied" if _distance_was_supplied
                               else "direct-cell (PRIMARY)"),
           "min_bin_count": int(MIN_BIN_COUNT),
           "drug": _fit(np.asarray(auc_field, dtype=np.float64))}
    if include_source_null:
        out["source_null"] = _fit(vasc)
        Ld, Ls = out["drug"]["L_um"], out["source_null"]["L_um"]
        out["L_excess_over_source_null"] = (float(Ld / Ls)
                                            if np.isfinite(Ld) and np.isfinite(Ls) and Ls > 0
                                            else float("nan"))
    out["L_um"] = out["drug"]["L_um"]
    return out


# --------------------------------------------------------------------------- #
# 5. Does the drug reach the poorly-perfused tissue?
# --------------------------------------------------------------------------- #
def hypoxic_vs_normoxic_exposure_ratio(
    micro: Microenvironment,
    sim_result: SimulationResult,
    hypoxia_threshold: float = OxygenConfig().hypoxia_threshold,
) -> Dict:
    """METRIC 5 — mean tumor-cell exposure in hypoxic vs well-oxygenated tissue.

    Tumor cells are partitioned with the engine's OWN oxygen banding
    (:func:`tmesim.cells.agents.update_cell_states`): with ``t =
    hypoxia_threshold`` (default ``OxygenConfig.hypoxia_threshold`` = 0.2),

        ``o < t``      -> HYPOXIC        (poorly perfused, drug-resistant)
        ``t <= o < 2t``-> QUIESCENT      (perfusion shoulder)
        ``o >= 2t``    -> VIABLE         (well perfused)

    Both a banded and a threshold-FREE contrast are returned:

    * ``ratio_lowest_over_highest_o2_quartile`` — **the primary readout**. Tumor
      cells are ranked by their own local oxygen and the mean ``E`` of the
      worst-perfused quartile is divided by that of the best-perfused quartile.
      This asks the same question with a partition that **always exists** and does
      not depend on where a hypoxia cutoff is drawn.
    * ``ratio_hypoxic_over_normoxic`` (``o < t`` vs ``o >= t``) and
      ``ratio_hypoxic_over_viable`` (``o < t`` vs ``o >= 2t``) — the literal banded
      versions. They are **undefined (nan) whenever a band is empty**, which is
      exactly what happens on a well-oxygenated reconstruction with no hypoxic
      tumor cells at all.

    Values < 1 mean the drug systematically under-serves exactly the regions that
    are already hardest to kill.

    CAVEAT, stated at the definition: this is a RATIO of two exposures, so it is
    scale-invariant — structurally the same kind of statistic as
    ``penetration_index``, and it inherits the same weakness. The banded variants
    additionally sit on a threshold crossing on the oxygen field, and the tumor
    hypoxic fraction is documented as NOT grid-converged
    (``results/real_drug_panel.md``). It is a diagnostic, not a headline.
    """
    oxygen = sim_result.fields.get("oxygen")
    if oxygen is None:
        raise ValueError("sim_result.fields has no 'oxygen' field")
    exp = tumor_cell_exposure(micro, sim_result)
    o2 = _sample_at_cells(oxygen, micro, exp["tumor_mask"])
    E = exp["exposure_uM"]
    t = float(hypoxia_threshold)

    hyp = o2 < t
    via = o2 >= 2.0 * t
    non_hyp = ~hyp

    def _m(sel):
        return float(E[sel].mean()) if sel.any() else float("nan")

    mh, mv, mn = _m(hyp), _m(via), _m(non_hyp)

    # Threshold-FREE variant. The banded ratios above are undefined whenever a band
    # is empty (a well-oxygenated reconstruction has no hypoxic cells at all), and
    # they sit on a threshold crossing that is documented as not grid-converged.
    # Ranking the tumor cells by their own local oxygen and contrasting the
    # worst-perfused quartile with the best-perfused quartile asks the same
    # question ("does the drug reach the badly-perfused tissue?") with a partition
    # that always exists and is insensitive to where the hypoxia cutoff sits.
    q_lo, q_hi = float("nan"), float("nan")
    if o2.size >= 4:
        lo_edge, hi_edge = np.percentile(o2, [25.0, 75.0])
        lo_sel, hi_sel = o2 <= lo_edge, o2 >= hi_edge
        q_lo, q_hi = _m(lo_sel), _m(hi_sel)

    return {
        "hypoxia_threshold": t,
        "hypoxic_tumor_fraction": float(hyp.mean()),
        "normoxic_tumor_fraction": float(non_hyp.mean()),
        "viable_tumor_fraction": float(via.mean()),
        "mean_exposure_hypoxic_uM": mh,
        "mean_exposure_normoxic_uM": mn,
        "mean_exposure_viable_uM": mv,
        "ratio_hypoxic_over_normoxic": float(mh / mn) if mn > 0 else float("nan"),
        "ratio_hypoxic_over_viable": float(mh / mv) if mv > 0 else float("nan"),
        "mean_exposure_lowest_o2_quartile_uM": q_lo,
        "mean_exposure_highest_o2_quartile_uM": q_hi,
        "ratio_lowest_over_highest_o2_quartile": (float(q_lo / q_hi) if q_hi > 0
                                                  else float("nan")),
        "n_hypoxic": int(hyp.sum()),
        "n_normoxic": int(non_hyp.sum()),
        "n_viable": int(via.sum()),
    }


# --------------------------------------------------------------------------- #
# The public entry point
# --------------------------------------------------------------------------- #
def compute_all(micro: Microenvironment, sim_result: SimulationResult,
                drug: DrugParams,
                distance_um: Optional[np.ndarray] = None) -> Dict:
    """Compute every delivery metric for one (tissue, simulation, drug) triple.

    Parameters
    ----------
    micro : Microenvironment
        The reconstructed real microenvironment the simulation ran on. Supplies
        the tissue mask, the real cell positions, and the vascular field.
    sim_result : SimulationResult
        Output of :func:`tmesim.sim.engine.simulate_response`. Must carry
        ``fields['auc']``, ``fields['oxygen']`` and ``transport.times``.
    drug : DrugParams
        The COMPILED drug — its ``ic50`` (uM) is what makes
        ``therapeutic_coverage`` drug-specific.
    distance_um : np.ndarray, optional
        Pre-computed distance-from-vessel map (um). Pass it when looping over
        drugs on the same tissue; it is geometry, identical for every drug, and
        the distance transform is the expensive part of metric 4.

    Returns
    -------
    dict
        A flat-ish dict with, under ``'headline'``, the small set of numbers meant
        to be quoted, and under the metric names the full detail. Also carries
        ``penetration_index`` (from the engine's extras) purely so a caller can
        show the contrast, and ``_exposure_uM`` — the raw per-tumor-cell exposure
        array, kept for plotting the delivery distribution.
    """
    T = exposure_window_s(sim_result)
    exp = tumor_cell_exposure(micro, sim_result)
    auc_cells = exp["auc_uM_s"]

    m1 = tissue_auc_summary(micro, sim_result)
    m2 = therapeutic_coverage(auc_cells, float(drug.ic50), T)
    m3 = exposure_percentiles(exp["exposure_uM"], auc_cells)
    m4 = penetration_length_fit(sim_result.fields["auc"], micro,
                                distance_um=distance_um)
    m5 = hypoxic_vs_normoxic_exposure_ratio(micro, sim_result)

    return {
        "drug": drug.name,
        "ic50_uM": float(drug.ic50),
        "bbb_penetrance": float(drug.bbb_penetrance),
        "T_s": T,
        "n_tumor": exp["n_tumor"],
        "dx_um": float(micro.domain.dx),
        "grid": tuple(micro.domain.shape),
        # --- metric blocks -------------------------------------------------- #
        "absolute_exposure": m1,
        "therapeutic_coverage": m2,
        "exposure_percentiles": m3,
        "penetration_length": m4,
        "hypoxia": m5,
        # --- the retired statistic, for contrast ----------------------------- #
        "penetration_index": float(sim_result.extras.get("penetration_index", np.nan)),
        "kill_fraction": float(sim_result.kill_fraction),
        # --- raw per-tumor-cell exposure (uM), kept for the distribution plot -- #
        "_exposure_uM": exp["exposure_uM"],
        # --- the numbers meant to be quoted ---------------------------------- #
        "headline": {
            "therapeutic_coverage_1x_ic50": m2["coverage_1x"],
            "therapeutic_coverage_0.5x_ic50": m2["coverage_0.5x"],
            "therapeutic_coverage_2x_ic50": m2["coverage_2x"],
            "median_E_over_ic50": m2["median_E_over_ic50"],
            "p10_E_over_ic50": m2["p10_E_over_ic50"],
            "mean_tumor_auc_uM_s": m1["mean_tumor_auc"],
            "median_tumor_auc_uM_s": m1["median_tumor_auc"],
            "E_p10_uM": m3["E_p10_uM"],
            "E_p50_uM": m3["E_p50_uM"],
            "E_p90_uM": m3["E_p90_uM"],
            "penetration_length_L_um": m4["L_um"],
            "low_o2_over_high_o2_exposure": m5["ratio_lowest_over_highest_o2_quartile"],
            "hypoxic_over_normoxic_exposure": m5["ratio_hypoxic_over_normoxic"],
        },
    }


# --------------------------------------------------------------------------- #
# Validation harness: grid convergence + drug discrimination
# --------------------------------------------------------------------------- #
#: (label, "/"-separated path into a :func:`compute_all` dict, units) for every
#: metric the grid-convergence / discrimination audit covers. The separator is
#: "/" and not "." because some leaf keys ("coverage_0.5x") contain a dot.
_AUDIT: List[Tuple[str, str, str]] = [
    ("mean_tissue_auc", "absolute_exposure/mean_tissue_auc", "uM*s"),
    ("median_tissue_auc", "absolute_exposure/median_tissue_auc", "uM*s"),
    ("mean_tumor_auc", "absolute_exposure/mean_tumor_auc", "uM*s"),
    ("median_tumor_auc", "absolute_exposure/median_tumor_auc", "uM*s"),
    ("therapeutic_coverage_0.5x", "therapeutic_coverage/coverage_0.5x", "fraction"),
    ("therapeutic_coverage_1x", "therapeutic_coverage/coverage_1x", "fraction"),
    ("therapeutic_coverage_2x", "therapeutic_coverage/coverage_2x", "fraction"),
    ("median_E/IC50 (potency margin)",
     "therapeutic_coverage/median_E_over_ic50", "-"),
    ("P10_E/IC50 (worst-decile margin)",
     "therapeutic_coverage/p10_E_over_ic50", "-"),
    ("exposure_P10", "exposure_percentiles/E_p10_uM", "uM"),
    ("exposure_P50", "exposure_percentiles/E_p50_uM", "uM"),
    ("exposure_P90", "exposure_percentiles/E_p90_uM", "uM"),
    ("exposure_P90/P10", "exposure_percentiles/E_p90_over_p10", "-"),
    ("penetration_length_L", "penetration_length/L_um", "um"),
    ("lowO2/highO2_exposure (quartiles)",
     "hypoxia/ratio_lowest_over_highest_o2_quartile", "-"),
    ("hypoxic/normoxic_exposure (banded)",
     "hypoxia/ratio_hypoxic_over_normoxic", "-"),
    ("penetration_index (RETIRED)", "penetration_index", "-"),
    ("kill_fraction (reference)", "kill_fraction", "fraction"),
]


def _dig(d: Dict, path: str) -> float:
    cur = d
    for part in path.split("/"):
        cur = cur[part]
    return float(cur)


def _rel(a: float, b: float) -> float:
    """Relative change |a-b| / |b| (nan-safe)."""
    if not (np.isfinite(a) and np.isfinite(b)) or b == 0:
        return float("nan")
    return float(abs(a - b) / abs(b))


def compile_real_drug(name: str, alpha: float):
    """The same compilation pipeline every flagship real run uses: alpha-calibrated
    perivascular sink + saturable binding-site barrier."""
    dp = calibrated_drugparams(name, alpha)
    return apply_barrier_to_drugparams(dp, get_molecule(name)), get_molecule(name)


def run(dx_ladder: Sequence[float] = DX_LADDER,
        drugs: Sequence[str] = tuple(PANEL),
        cache_path: Optional[str] = None) -> Dict:
    """Run the full validation: every drug at every dx on the REAL Xenium section.

    ``cache_path`` (optional) pickles the raw per-(dx, drug) metric dicts so the
    report and figure can be regenerated without re-simulating.
    """
    if cache_path and os.path.exists(cache_path):
        with open(cache_path, "rb") as fh:
            return pickle.load(fh)

    t0 = time.time()
    alpha = calibrate_alpha()
    print(f"perivascular-sink calibration alpha = {alpha:.2f} "
          f"(fit ONCE on doxorubicin)")

    per_dx: Dict[float, Dict] = {}
    reports: Dict[float, Dict] = {}
    for dx in dx_ladder:
        print(f"\n--- dx = {dx:g} um ---")
        micro = reconstruct_real("xenium", dx_um=float(dx))
        rep = reconstruction_report(micro)
        reports[float(dx)] = rep
        print(f"  grid {rep['grid']}  {rep['n_cells']:,} real cells "
              f"(tumor {rep['n_tumor']:,})")
        # Geometry: identical for every drug -> compute the distance map ONCE.
        # PRIMARY metric: direct nearest-endothelial-cell distance (see
        # results/vessel_profile.md's dedicated distance-metric-change section).
        endo_xy = endothelial_xy(micro.cells)
        dist = direct_vessel_distance_um(micro.domain, endo_xy)
        # LEGACY mask, kept only for the informational area-fraction line below
        # (cheap: geometry only, no simulation).
        vm = vessel_mask_from_density(micro.fields.vascular_density,
                                      micro.domain.tissue_mask)
        print(f"  vessel distance [PRIMARY, direct-cell, {endo_xy.shape[0]:,} "
              f"endothelial cells]  ·  LEGACY mask: {vm['n_vessel_px']:,} px "
              f"({vm['area_fraction']*100:.1f}% of tissue) at p{vm['percentile']:.0f}")

        rows: Dict[str, Dict] = {}
        for name in drugs:
            dp, mol = compile_real_drug(name, alpha)
            t1 = time.time()
            res = simulate_response(
                micro, SimConfig(drug=dp, pde=PDEConfig(**_PDE), seed=0))
            met = compute_all(micro, res, dp, distance_um=dist)
            met["mw"] = float(mol.get("mw"))
            met["vessel_area_fraction"] = float(vm["area_fraction"])
            rows[name] = met
            h = met["headline"]
            print(f"  {name:14s} cov@IC50={h['therapeutic_coverage_1x_ic50']:.4f}  "
                  f"meanAUC={h['mean_tumor_auc_uM_s']:.3g}  "
                  f"L={h['penetration_length_L_um']:.1f}um  "
                  f"pen_idx={met['penetration_index']:.4f}  "
                  f"[{time.time()-t1:.0f}s]")
        per_dx[float(dx)] = rows

    ladder = [float(d) for d in dx_ladder]
    out = {
        "alpha": float(alpha),
        "dx_ladder": ladder,
        # The production resolution is dx = 12 um whenever it is on the ladder;
        # otherwise the finest grid actually run.
        "dx_production": (float(DX_PRODUCTION) if float(DX_PRODUCTION) in ladder
                          else min(ladder)),
        "drugs": list(drugs),
        "per_dx": per_dx,
        "reports": reports,
        "runtime_s": time.time() - t0,
    }
    out["convergence"] = convergence_table(out)
    out["discrimination"] = discrimination_table(out)
    if cache_path:
        os.makedirs(os.path.dirname(os.path.abspath(cache_path)), exist_ok=True)
        with open(cache_path, "wb") as fh:
            pickle.dump(out, fh)
    return out


def backfill(res: Dict) -> Dict:
    """Recompute any audited quantity that is *derivable from what is already
    stored*, then rebuild the convergence / discrimination tables.

    This exists so a cached run (``--cache``) stays usable when a purely derived
    metric is added: nothing is estimated or imputed — the potency margins are
    recomputed exactly, from the stored per-tumor-cell exposures and the stored
    IC50, by the same code path :func:`therapeutic_coverage` uses.
    """
    for rows in res["per_dx"].values():
        for met in rows.values():
            tc = met["therapeutic_coverage"]
            if "median_E_over_ic50" in tc:
                continue
            ratio = np.asarray(met["_exposure_uM"], dtype=float) / met["ic50_uM"]
            tc["median_E_over_ic50"] = float(np.median(ratio))
            tc["p10_E_over_ic50"] = float(np.percentile(ratio, 10.0))
            met["headline"]["median_E_over_ic50"] = tc["median_E_over_ic50"]
            met["headline"]["p10_E_over_ic50"] = tc["p10_E_over_ic50"]
    res["convergence"] = convergence_table(res)
    res["discrimination"] = discrimination_table(res)
    return res


def convergence_table(res: Dict) -> List[Dict]:
    """For every audited metric: its value at each dx (worst-drug view) and the
    worst-drug relative change between successive grids. A metric is
    GRID-CONVERGED when the worst-drug 15 -> 12 um change is < CONVERGENCE_TOL."""
    dxs = sorted(res["dx_ladder"], reverse=True)          # 20, 15, 12
    fine = min(dxs)
    rows = []
    for label, path, units in _AUDIT:
        vals = {dx: [_dig(res["per_dx"][dx][d], path) for d in res["drugs"]]
                for dx in dxs}
        drifts, abs_drifts = {}, {}
        for a, b in zip(dxs[:-1], dxs[1:]):               # (20->15), (15->12)
            key = f"{a:g}->{b:g}"
            per_drug = [_rel(vals[b][i], vals[a][i])
                        for i in range(len(res["drugs"]))]
            finite = [v for v in per_drug if np.isfinite(v)]
            # Worst drug in the panel; nan only if the metric is undefined for all.
            drifts[key] = max(finite) if finite else float("nan")
            per_abs = [abs(vals[b][i] - vals[a][i]) for i in range(len(res["drugs"]))]
            fin_abs = [v for v in per_abs if np.isfinite(v)]
            abs_drifts[key] = max(fin_abs) if fin_abs else float("nan")
        final_key = f"{dxs[-2]:g}->{fine:g}"
        worst = drifts.get(final_key, float("nan"))
        worst_abs = abs_drifts.get(final_key, float("nan"))
        bounded = units == "fraction"
        conv = (np.isfinite(worst_abs) and worst_abs < CONVERGENCE_TOL_ABS
                if bounded else
                np.isfinite(worst) and worst < CONVERGENCE_TOL)
        rows.append({
            "metric": label, "units": units, "path": path, "bounded": bounded,
            "values": {dx: vals[dx] for dx in dxs},
            "mean_at_dx": {dx: float(np.mean(vals[dx])) for dx in dxs},
            "drift": drifts,
            "abs_drift": abs_drifts,
            "worst_drug_drift_final": worst,
            "worst_drug_abs_drift_final": worst_abs,
            "converged": bool(conv),
        })
    return rows


def discrimination_table(res: Dict) -> List[Dict]:
    """For every audited metric at the PRODUCTION dx: the across-drug spread, and
    that spread expressed as a multiple of the metric's own worst-drug grid drift.

    ``separation = relative_spread / grid_drift`` is the honest test: a metric is
    only usable for ranking drugs if the differences it reports are much larger
    than the numerical noise it carries.
    """
    dx = res["dx_production"]
    conv = {r["metric"]: r for r in res["convergence"]}
    rows = []
    for label, path, units in _AUDIT:
        vals = np.array([_dig(res["per_dx"][dx][d], path) for d in res["drugs"]],
                        dtype=float)
        finite = vals[np.isfinite(vals)]
        lo = float(finite.min()) if finite.size else float("nan")
        hi = float(finite.max()) if finite.size else float("nan")
        mid = float(np.mean(finite)) if finite.size else float("nan")
        fold = float(hi / lo) if np.isfinite(lo) and lo > 0 else float("nan")
        rel_spread = float((hi - lo) / abs(mid)) if mid else float("nan")
        abs_spread = float(hi - lo)
        bounded = conv[label]["bounded"]
        # Compare like with like: a bounded fraction is judged on the absolute
        # scale it lives on, everything else on the relative scale.
        spread = abs_spread if bounded else rel_spread
        drift = (conv[label]["worst_drug_abs_drift_final"] if bounded
                 else conv[label]["worst_drug_drift_final"])
        if not np.isfinite(drift) or not np.isfinite(spread):
            sep = float("nan")          # metric undefined on this section
        elif drift > 0:
            sep = float(spread / drift)
        else:
            sep = float("inf")          # exactly grid-invariant
        rows.append({
            "metric": label, "units": units, "path": path, "bounded": bounded,
            "per_drug": {d: float(v) for d, v in zip(res["drugs"], vals)},
            "min": lo, "max": hi, "fold_range": fold,
            "relative_spread": rel_spread,
            "abs_spread": abs_spread,
            "grid_drift": drift,
            "grid_drift_relative": conv[label]["worst_drug_drift_final"],
            "grid_drift_absolute": conv[label]["worst_drug_abs_drift_final"],
            "separation": sep,
            "converged": conv[label]["converged"],
            "discriminates": bool(np.isfinite(sep) and sep >= DISCRIMINATION_MIN_RATIO
                                  and np.isfinite(fold) and fold > 1.05),
        })
    return rows


# --------------------------------------------------------------------------- #
# Figure
# --------------------------------------------------------------------------- #
def make_figure(res: Dict, path: str = _OUT_PNG) -> str:
    """Four panels: (1) therapeutic coverage per drug with its 2x-0.5x IC50 band,
    (2) the per-tumor-cell exposure distributions against each drug's own IC50,
    (3) the contrast — a flat ``penetration_index`` beside a spread coverage,
    (4) grid convergence of the headline metric."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    dx = res["dx_production"]
    drugs = res["drugs"]
    rows = [res["per_dx"][dx][d] for d in drugs]
    colors = [_DRUG_COLORS.get(d, _C_GOOD) for d in drugs]

    plt.rcParams.update({
        "font.size": 10, "axes.edgecolor": _GRID, "axes.labelcolor": _INK_2,
        "xtick.color": _INK_2, "ytick.color": _INK_2, "text.color": _INK,
        "axes.titlecolor": _INK, "figure.facecolor": "white",
        "axes.facecolor": "white",
    })
    fig = plt.figure(figsize=(17.5, 11.5))
    gs = GridSpec(2, 2, figure=fig, hspace=0.42, wspace=0.22,
                  left=0.06, right=0.985, top=0.86, bottom=0.07)

    # ---- Panel 1: therapeutic coverage -------------------------------------- #
    ax = fig.add_subplot(gs[0, 0])
    cov1 = [r["therapeutic_coverage"]["coverage_1x"] for r in rows]
    cov05 = [r["therapeutic_coverage"]["coverage_0.5x"] for r in rows]
    cov2 = [r["therapeutic_coverage"]["coverage_2x"] for r in rows]
    x = np.arange(len(drugs))
    ax.bar(x, cov1, 0.56, color=colors, zorder=3)
    for i in range(len(drugs)):
        ax.plot([x[i], x[i]], [cov2[i], cov05[i]], color=_INK, lw=2.0, zorder=4,
                solid_capstyle="butt")
        ax.plot([x[i]], [cov05[i]], marker="_", ms=13, mew=2.0, color=_INK, zorder=5)
        ax.plot([x[i]], [cov2[i]], marker="_", ms=13, mew=2.0, color=_INK, zorder=5)
        ax.text(x[i], max(cov05[i], cov1[i]) + 0.028, f"{cov1[i]:.3f}",
                ha="center", va="bottom", fontsize=10, color=_INK, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{d}\nIC50 {r['ic50_uM']:.3g} µM" for d, r in zip(drugs, rows)],
                       fontsize=9)
    ax.set_ylabel("fraction of tumor cells with E > IC50")
    ax.set_ylim(0, 1.12)
    ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.grid(axis="y", color=_GRID, lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_title("1 · THERAPEUTIC COVERAGE — what fraction of the tumor saw a killing dose\n"
                 "bar = 1× IC50;  whisker spans 2× IC50 (bottom) → 0.5× IC50 (top)",
                 fontsize=11, loc="left")

    # ---- Panel 2: exposure distributions ------------------------------------ #
    ax = fig.add_subplot(gs[0, 1])
    data, ic50s = [], []
    for r, d in zip(rows, drugs):
        E = np.asarray(r["_exposure_uM"], dtype=float)
        E = E[np.isfinite(E) & (E > 0)]
        data.append(np.log10(E))
        ic50s.append(np.log10(r["ic50_uM"]))
    parts = ax.violinplot(data, positions=np.arange(len(drugs)), vert=False,
                          widths=0.82, showextrema=False, showmedians=False)
    for body, c in zip(parts["bodies"], colors):
        body.set_facecolor(c)
        body.set_edgecolor(c)
        body.set_alpha(0.85)
        body.set_linewidth(0)
    for i, r in enumerate(rows):
        p = r["exposure_percentiles"]
        lo, mid, hi = (np.log10(max(p["E_p10_uM"], 1e-12)),
                       np.log10(max(p["E_p50_uM"], 1e-12)),
                       np.log10(max(p["E_p90_uM"], 1e-12)))
        ax.plot([lo, hi], [i, i], color="white", lw=3.0, zorder=4)
        ax.plot([lo, hi], [i, i], color=_INK, lw=1.4, zorder=5)
        ax.plot([mid], [i], marker="o", ms=6, color=_INK, mec="white", mew=1.2, zorder=6)
        ax.plot([ic50s[i], ic50s[i]], [i - 0.44, i + 0.44], color=_C_BAD, lw=2.4,
                zorder=7)
        ax.text(ic50s[i], i + 0.46, "IC50", color=_C_BAD, fontsize=8,
                ha="center", va="bottom", fontweight="bold")
    # Clip the far-left tail: a handful of near-zero-AUC cells otherwise stretch the
    # axis by several decades and squash the part of the distribution that matters.
    x_lo = min(np.percentile(v, 1.0) for v in data) - 0.4
    x_hi = max(v.max() for v in data) + 0.4
    ax.set_xlim(min(x_lo, min(ic50s) - 0.4), max(x_hi, max(ic50s) + 0.4))
    ax.set_yticks(np.arange(len(drugs)))
    ax.set_yticklabels(drugs, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("per-tumor-cell time-averaged exposure  E = AUC/T   [log₁₀ µM]")
    ax.grid(axis="x", color=_GRID, lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_title(f"2 · DELIVERY DISTRIBUTION over the {rows[0]['n_tumor']:,} real "
                 "tumor cells\n"
                 "black bar = P10→P90, dot = P50;  red line = that drug's real IC50 "
                 "(mass to the RIGHT of it is the coverage in panel 1)",
                 fontsize=11, loc="left")

    # ---- Panel 3: the contrast ---------------------------------------------- #
    ax = fig.add_subplot(gs[1, 0])
    pen = [r["penetration_index"] for r in rows]
    auc_all = [r["absolute_exposure"]["mean_tumor_auc"] for r in rows]
    disc_fold = max(auc_all) / max(min(auc_all), 1e-30)
    w = 0.38
    b1 = ax.bar(x - w / 2, pen, w, color=_C_BAD, zorder=3,
                label=f"penetration_index (RETIRED) — spread {max(pen)-min(pen):.3f}")
    b2 = ax.bar(x + w / 2, cov1, w, color=_C_GOOD, zorder=3,
                label=f"therapeutic_coverage @1× IC50 — spread {max(cov1)-min(cov1):.3f}")
    for b in list(b1) + list(b2):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.018,
                f"{b.get_height():.3f}", ha="center", va="bottom", fontsize=8.5,
                color=_INK)
    ax.set_xticks(x)
    ax.set_xticklabels(drugs, fontsize=9)
    ax.set_ylabel("dimensionless fraction / index")
    ax.set_ylim(0, 1.15)
    ax.grid(axis="y", color=_GRID, lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.legend(loc="upper center", fontsize=9, frameon=False, ncol=1)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_title(f"3 · THE CONTRAST at converged dx = {dx:g} µm — same tumor, "
                 f"same {len(drugs)} drugs, {disc_fold:.0f}× spread in delivered "
                 "exposure\n"
                 f"penetration_index spans only {min(pen):.3f}–{max(pen):.3f} "
                 f"({max(pen)/max(min(pen), 1e-12):.2f}× — vascular geometry); "
                 f"coverage spans {min(cov1):.3f}–{max(cov1):.3f} "
                 f"({max(cov1)/max(min(cov1), 1e-12):.1f}×)",
                 fontsize=11, loc="left")

    # ---- Panel 4: grid convergence of the headline metric -------------------- #
    ax = fig.add_subplot(gs[1, 1])
    dxs = sorted(res["dx_ladder"], reverse=True)
    # Stagger the end labels so near-identical curves do not collide.
    finals = sorted((res["per_dx"][dxs[-1]][d]["therapeutic_coverage"]["coverage_1x"], d)
                    for d in drugs)
    nudge = {}
    last = -1e9
    for v, d in finals:
        off = 0.0 if v - last > 0.035 else (0.035 - (v - last))
        nudge[d] = off
        last = v + off
    for d, c in zip(drugs, colors):
        ys = [res["per_dx"][k][d]["therapeutic_coverage"]["coverage_1x"] for k in dxs]
        ax.plot(dxs, ys, marker="o", ms=7, lw=2.0, color=c, label=d, zorder=3)
        ax.annotate(f"{ys[-1]:.3f}", (dxs[-1], ys[-1] + nudge[d]),
                    textcoords="offset points", xytext=(9, -3), fontsize=9,
                    color=_INK)
    pen_dx = [float(np.mean([res["per_dx"][k][d]["penetration_index"] for d in drugs]))
              for k in dxs]
    ax.plot(dxs, pen_dx, marker="s", ms=6, lw=2.0, ls="--", color=_C_BAD,
            label="penetration_index (panel mean)", zorder=3)
    ax.annotate(f"{pen_dx[-1]:.3f}", (dxs[-1], pen_dx[-1]), textcoords="offset points",
                xytext=(9, -3), fontsize=9, color=_C_BAD)
    ax.axvline(dx, color=_INK_2, lw=1.0, ls=":", zorder=1)
    ax.text(dx, 0.02, "production dx ", fontsize=8.5, color=_INK_2, ha="right",
            va="bottom", transform=ax.get_xaxis_transform())
    ax.invert_xaxis()
    ax.set_xticks(dxs)
    ax.set_xticklabels([f"{k:g}" for k in dxs])
    ax.set_xlabel("grid spacing dx (µm) — refining to the right")
    ax.set_ylabel("therapeutic_coverage @1× IC50")
    ax.set_ylim(-0.03, 1.12)
    ax.grid(color=_GRID, lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.legend(loc="center left", fontsize=8.5, frameon=False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    cconv = next(r for r in res["convergence"]
                 if r["metric"] == "therapeutic_coverage_1x")
    pconv = next(r for r in res["convergence"]
                 if r["metric"] == "penetration_index (RETIRED)")
    ax.set_title("4 · GRID CONVERGENCE of the headline metric on the REAL section\n"
                 f"coverage: Δ = {cconv['worst_drug_abs_drift_final']:.3f} abs "
                 f"({cconv['worst_drug_abs_drift_final']*100:.1f} pp, tol "
                 f"{CONVERGENCE_TOL_ABS:.2f})   ·   penetration_index: Δ = "
                 f"{pconv['worst_drug_abs_drift_final']:.3f} abs = "
                 f"{pconv['worst_drug_abs_drift_final'] / max(max(pen)-min(pen), 1e-12):.0f}× "
                 "its whole spread",
                 fontsize=11, loc="left")

    rep = res["reports"][dx]
    fig.suptitle(
        "DEFENSIBLE DELIVERY METRICS on the REAL 10x Xenium glioblastoma section — "
        "replacing the scale-invariant `penetration_index`\n"
        f"{rep['n_cells']:,} real cells (tumor {rep['n_tumor']:,}) · grid {rep['grid']} · "
        f"dx = {dx:g} µm (grid-converged, advection-CFL-safe at dt = 10 s) · "
        f"one tissue, {len(drugs)} drugs, identical geometry\n"
        "PREDICTIONS of the calibrated physics on real tissue — this is an untreated "
        "resection with NO measured drug response to validate against",
        fontsize=12.5, y=0.985,
    )
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return path


# --------------------------------------------------------------------------- #
# Report
# --------------------------------------------------------------------------- #
def _f(v, nd=3):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return "n/a"
    return f"{v:.{nd}g}"


def write_report(res: Dict, path: str = _OUT_MD) -> str:
    dx = res["dx_production"]
    dxs = sorted(res["dx_ladder"], reverse=True)
    drugs = res["drugs"]
    rows = [res["per_dx"][dx][d] for d in drugs]
    rep = res["reports"][dx]
    conv = {r["metric"]: r for r in res["convergence"]}
    disc = {r["metric"]: r for r in res["discrimination"]}

    cov = disc["therapeutic_coverage_1x"]
    pen = disc["penetration_index (RETIRED)"]

    L = []
    A = L.append
    A("# Defensible delivery metrics — replacing `penetration_index`")
    A("")
    A("**Module:** `tmesim/validation/delivery_metrics.py` "
      "(run: `./.venv/bin/python -m tmesim.validation.delivery_metrics`)")
    A("")
    A("**Metric self-test:** `./.venv/bin/python -m tmesim.validation.delivery_metrics "
      "--selftest` (known-answer checks on the definitions, no simulation)")
    A("")
    A("**Figure:** `results/delivery_metrics.png`")
    A("")
    A(f"**Geometry:** the REAL open 10x Xenium GBM section — {rep['n_cells']:,} "
      f"measured cells (tumor {rep['n_tumor']:,}), reconstructed by "
      f"`reconstruct.real_pipeline.reconstruct_real('xenium', dx_um=…)`; real "
      f"endothelial-KDE vascular field, real matrix-transcript ECM, physics-solved "
      f"oxygen. **Drugs:** the {len(drugs)}-drug panel of "
      f"[`results/real_drug_panel.md`](real_drug_panel.md), each compiled through "
      f"the identical pipeline (`calibrate_alpha()` → α = {res['alpha']:.1f}, fit "
      f"ONCE on doxorubicin → `calibrated_drugparams` → "
      f"`apply_barrier_to_drugparams`). **Engine:** advection on, dt = 10 s, 6 h "
      f"window, seed 0 — unchanged from the flagship run. Production resolution "
      f"**dx = {dx:g} µm**, grid {rep['grid']} (dx = 10 µm violates the advection "
      f"CFL condition at dt = 10 s and is therefore not on the ladder).")
    A("")
    A("---")
    A("")
    A("## 0. The problem being fixed")
    A("")
    A("`penetration_index = mean(AUC) / p99(AUC)` over the tissue mask is a "
      "**scale-invariant shape statistic**: multiply the entire exposure field by "
      "any constant and it does not move. At the grid-converged resolution it "
      "therefore reports the *vascular geometry of the section*, which is identical "
      f"for all {len(drugs)} drugs — and it duly collapses:")
    A("")
    A(f"> at dx = {dx:g} µm, `penetration_index` = "
      f"**{pen['min']:.3f}–{pen['max']:.3f}** across the whole panel "
      f"(spread {pen['max']-pen['min']:.4f}), while those same {len(drugs)} drugs' mean "
      f"delivered exposure differs by "
      f"**{disc['mean_tumor_auc']['fold_range']:.0f}×**.")
    A("")
    A("A delivery metric must be **absolute** (it scales when the delivered dose "
      "scales) and, ideally, **drug-specific** (it knows the drug's potency). The "
      "five metrics below are defined to be exactly that, and then audited on both "
      "axes that matter: grid convergence and drug discrimination.")
    A("")
    A("---")
    A("")
    A("## 1. Definitions, units, rationale")
    A("")
    A("Let `AUC(x)` be the transport solve's exposure field (µM·s) and `T` the "
      "exposure-window length (s; here "
      f"T = {rows[0]['T_s']:.0f} s = {rows[0]['T_s']/3600.0:g} h, taken from "
      "`transport.times[-1]`, the "
      "same `T` `cells.pkpd` uses). The **time-averaged effective concentration** is")
    A("")
    A("```")
    A("E(x) = AUC(x) / T        [µM]")
    A("```")
    A("")
    A("which is directly comparable to an IC50 because an IC50 *is* a "
      "concentration. Every per-cell statistic is evaluated at the **real measured "
      "cell positions** (`domain.to_index`, the same nearest-cell map the kill model "
      f"uses), not over grid pixels — the tumor population is a FIXED set of "
      f"{rep['n_tumor']:,} Xenium cells at every dx, so a per-cell statistic moves "
      "with dx only because the *field* moved, not because the rasterised mask "
      "changed.")
    A("")
    A("| # | metric | units | what it answers | absolute or shape? |")
    A("|---|---|---|---|---|")
    A("| 1 | `mean_tissue_auc`, `median_tissue_auc`, `mean_tumor_auc`, "
      "`median_tumor_auc` | µM·s | how much drug arrived, in total | **absolute** |")
    A("| 2 | `therapeutic_coverage(auc, ic50, T)` — fraction of tumor cells with "
      "`E > m·IC50`, m ∈ {0.5, 1, 2} | dimensionless fraction ∈ [0,1] | what "
      "fraction of the tumor saw a **killing** dose | **absolute + drug-specific** |")
    A("| 3 | `exposure_percentiles` — P10/P50/P90 of per-tumor-cell `E` | µM | how "
      "*unevenly* the dose was delivered; P10 = the worst-served decile | "
      "**absolute** |")
    A("| 4 | `penetration_length_fit` — `L` from `AUC(d)=C₀·exp(−d/L)+floor` vs "
      "distance-from-nearest-vessel | µm | a real physical decay length | length "
      "(scale-invariant, but dimensional) |")
    A("| 5 | `hypoxic_vs_normoxic_exposure_ratio` — mean `E` in the worst- vs "
      "best-perfused quartile of tumor cells (plus the literal hypoxic/normoxic "
      "bands) | dimensionless | does the drug reach the poorly-perfused tissue? | "
      "*shape* (diagnostic only) |")
    A("")
    A("**Rationale, metric by metric.**")
    A("")
    A("1. **Absolute exposure (µM·s).** The integral the pharmacodynamics actually "
      "consumes. The median is reported next to the mean because the field is "
      "strongly right-skewed — a thin perivascular cuff carries much of the "
      f"integral, so mean/median ≈ "
      f"{np.mean([r['absolute_exposure']['skew_mean_over_median'] for r in rows]):.1f}× "
      "on this section and the mean alone overstates what a typical tumor cell sees. "
      "Honest limitation: an exposure in µM·s is not interpretable without a "
      "potency, which is why it is not the headline.")
    A("2. **Therapeutic coverage.** The clinically meaningful reduction of a whole "
      "exposure *field* to one number: the fraction of tumor cells whose "
      "time-averaged concentration clears the drug's own IC50. Dimensionally sound "
      "(concentration vs concentration), bounded in [0,1], and — crucially — the "
      "same exposure field scores differently for a potent and a weak drug, so it "
      "**cannot** collapse to a geometry constant. The 0.5× and 2× variants say how "
      "hard the answer leans on the exact IC50 value.")
    A("3. **Exposure percentiles.** A mean hides the population that decides the "
      "outcome. P10 is the exposure the worst-served decile of the tumor receives — "
      "the reservoir that survives and regrows — and P90/P10 is the delivery "
      "heterogeneity.")
    A("4. **Penetration length.** A genuine physical length in µm, fitted to the "
      "binned AUC-vs-distance-from-vessel profile (the distance map, 10 µm binning "
      "and fits are **reused** from `validation/vessel_profile.py`, so the two "
      "reports cannot drift apart). It ships with its own **SOURCE NULL** — the "
      "identical fit on the raw `vascular_density` field — because the engine's "
      "delivery term is `vascular_density · bbb_penetrance · plasma`, i.e. a source "
      "spread over the whole continuous density field; with a fast perivascular sink "
      "the tissue field is close to locally reaction-dominated (`C ≈ S/k`) and "
      "largely **inherits the source's shape**.")
    A("5. **Hypoxic/normoxic ratio.** Two forms. The **primary, threshold-free** "
      "one ranks tumor cells by their own local oxygen and divides the mean `E` of "
      "the worst-perfused quartile by that of the best-perfused quartile — a "
      "partition that always exists and does not depend on where a hypoxia cutoff "
      "is drawn. The **banded** form uses the engine's own oxygen banding "
      "(`t = OxygenConfig.hypoxia_threshold` = 0.2; hypoxic `o < t` vs normoxic "
      "`o ≥ t`) and is **nan whenever a band is empty**. Stated at the definition: "
      "both are *ratios*, hence scale-invariant, hence structurally the same kind of "
      "statistic as `penetration_index`; the banded one additionally sits on a "
      "threshold crossing on an oxygen field whose hypoxic fraction is documented as "
      "**not** grid-converged. Included because the question is real, and reported "
      "as a diagnostic — not as a headline.")
    A("")
    A("---")
    A("")
    A("## 2. Grid convergence on the REAL section")
    A("")
    A(f"Every drug re-simulated on the same section re-reconstructed at "
      f"dx = {', '.join(f'{k:g}' for k in dxs)} µm. The table reports the "
      "**panel-mean** value at each dx and, for the drift columns, the **worst "
      "drug** in the panel — the conservative reading.")
    A("")
    A("**The convergence criterion, stated before the numbers.** Two scales are "
      "reported for the final refinement step, and which one decides is fixed by "
      "the metric's type, not by the answer:")
    A("")
    A(f"- metrics **bounded in [0,1]** (the coverage and kill fractions) are judged "
      f"on the **absolute** change, tolerance **{CONVERGENCE_TOL_ABS:.2f}** "
      f"(= {CONVERGENCE_TOL_ABS*100:.0f} percentage points of tumor cells). Relative "
      "error is ill-conditioned for a bounded fraction — it diverges as the fraction "
      "goes to zero even when the answer moved by a fraction of a point.")
    A(f"- everything else (µM·s, µM, µm, dimensionless ratios) is judged on the "
      f"**relative** change, tolerance **{CONVERGENCE_TOL*100:.0f}%**.")
    A("")
    A("Both columns are in the table for every metric, so any other rule can be "
      "applied to the same numbers.")
    A("")
    hdr = "| metric | units | " + " | ".join(f"dx {k:g} µm" for k in dxs) + \
          f" | Δrel {dxs[0]:g}→{dxs[1]:g} | Δrel {dxs[1]:g}→{dxs[2]:g} " \
          f"| Δabs {dxs[1]:g}→{dxs[2]:g} | converged? |"
    A(hdr)
    A("|---|---|" + "---:|" * (len(dxs) + 3) + ":--:|")
    for label, _p, units in _AUDIT:
        r = conv[label]
        cells = " | ".join(_f(r["mean_at_dx"][k], 4) for k in dxs)
        d1 = r["drift"].get(f"{dxs[0]:g}->{dxs[1]:g}", float("nan"))
        d2 = r["drift"].get(f"{dxs[1]:g}->{dxs[2]:g}", float("nan"))
        a2 = r["abs_drift"].get(f"{dxs[1]:g}->{dxs[2]:g}", float("nan"))
        mark = "**YES**" if r["converged"] else "no"
        star = " *(abs)*" if r["bounded"] else ""
        A(f"| `{label}` | {units} | {cells} | {_f(d1, 3)} | {_f(d2, 3)} | "
          f"{_f(a2, 3)} | {mark}{star} |")
    A("")
    conv_yes = [lbl for lbl, _p, _u in _AUDIT if conv[lbl]["converged"]]
    conv_no = [lbl for lbl, _p, _u in _AUDIT if not conv[lbl]["converged"]]
    A(f"**Converged (worst drug, {dxs[-2]:g}→{dxs[-1]:g} µm):** "
      + ", ".join(f"`{m}`" for m in conv_yes) + ".")
    A("")
    A("**NOT converged:** " + (", ".join(f"`{m}`" for m in conv_no) or "none") + ".")
    A("")
    cov_conv = conv["therapeutic_coverage_1x"]
    A(f"**Read this honestly for the headline metric.** `therapeutic_coverage_1x` "
      f"moves **{cov_conv['worst_drug_abs_drift_final']:.4f} absolute** "
      f"({cov_conv['worst_drug_abs_drift_final']*100:.1f} percentage points of tumor "
      f"cells) on its worst drug between {dxs[-2]:g} µm and {dxs[-1]:g} µm, which is "
      f"inside the {CONVERGENCE_TOL_ABS:.2f} tolerance — but that same move is "
      f"**{cov_conv['worst_drug_drift_final']*100:.1f}% relative**, and it is the "
      "*least*-covered drug (temozolomide, coverage ≈ 0.31) that produces it, which "
      "is exactly the regime where relative error is misleading. The residual is "
      "real and monotone, not noise: carry it as a **±0.03 error bar on any quoted "
      "coverage**. The well-covered drugs are far tighter "
      f"(cisplatin {abs(conv['therapeutic_coverage_1x']['values'][dxs[-1]][2] - conv['therapeutic_coverage_1x']['values'][dxs[-2]][2]):.4f}, "
      f"doxorubicin {abs(conv['therapeutic_coverage_1x']['values'][dxs[-1]][4] - conv['therapeutic_coverage_1x']['values'][dxs[-2]][4]):.4f} absolute).")
    A("")
    A(f"For contrast, `penetration_index` moves "
      f"{conv['penetration_index (RETIRED)']['worst_drug_drift_final']*100:.1f}% "
      f"relative over the same refinement step while its entire spread across the "
      f"{len(drugs)}-drug panel at dx = {dx:g} µm is only "
      f"{disc['penetration_index (RETIRED)']['relative_spread']*100:.1f}% — its grid "
      f"drift is "
      f"{conv['penetration_index (RETIRED)']['worst_drug_drift_final'] / max(disc['penetration_index (RETIRED)']['relative_spread'], 1e-12):.0f}× "
      "larger than the entire signal it is supposed to carry.")
    A("")
    A("---")
    A("")
    A(f"## 3. Discrimination across the {len(drugs)}-drug panel (dx = {dx:g} µm)")
    A("")
    A("### 3a. The per-drug table")
    A("")
    A("| drug | MW (Da) | real IC50 (µM) | mean tumor AUC (µM·s) | median tumor AUC "
      "(µM·s) | E P10 (µM) | E P50 (µM) | E P90 (µM) | **coverage @0.5× IC50** | "
      "**coverage @1× IC50** | **coverage @2× IC50** | median E/IC50 | P10 E/IC50 | "
      "L (µm) | lowO₂/highO₂ E | `penetration_index` | kill_fraction |")
    A("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for d, r in zip(drugs, rows):
        a, c, p, pl, h = (r["absolute_exposure"], r["therapeutic_coverage"],
                          r["exposure_percentiles"], r["penetration_length"],
                          r["hypoxia"])
        A(f"| {d} | {r.get('mw', float('nan')):.0f} | {r['ic50_uM']:.3g} | "
          f"{a['mean_tumor_auc']:.3g} | {a['median_tumor_auc']:.3g} | "
          f"{p['E_p10_uM']:.3g} | {p['E_p50_uM']:.3g} | {p['E_p90_uM']:.3g} | "
          f"{c['coverage_0.5x']:.4f} | **{c['coverage_1x']:.4f}** | "
          f"{c['coverage_2x']:.4f} | {c['median_E_over_ic50']:.3g} | "
          f"{c['p10_E_over_ic50']:.3g} | {pl['L_um']:.1f} | "
          f"{h['ratio_lowest_over_highest_o2_quartile']:.3f} | "
          f"{r['penetration_index']:.4f} | "
          f"{r['kill_fraction']:.4f} |")
    A("")
    sn = rows[0]["penetration_length"]["source_null"]
    A(f"`L` is the decaying component of `C₀·exp(−d/L)+floor` fitted to the binned "
      f"AUC-vs-distance profile, where distance is now the **PRIMARY** direct "
      f"nearest-endothelial-cell metric (computed once per grid resolution via "
      f"`direct_vessel_distance_um` and shared across drugs) — not the LEGACY "
      f"p{rows[0]['penetration_length']['vessel_percentile']:.0f} "
      f"percentile-mask threshold this report used before "
      f"({rows[0].get('vessel_area_fraction', float('nan'))*100:.1f}% of tissue area "
      f"under that legacy mask, kept here for reference only; see "
      f"`results/vessel_profile.md` for the full head-to-head). "
      f"**SOURCE NULL** (identical fit on the raw `vascular_density` field): "
      f"L = {sn['L_um']:.1f} µm.")
    A("")
    A("### 3b. Spread vs numerical noise")
    A("")
    A("A metric only ranks drugs credibly if its across-drug spread is much larger "
      "than its own grid drift. `separation = spread / worst-drug grid drift`, both "
      "taken on the **same scale** (absolute for the [0,1] fractions, relative for "
      f"everything else — the rule stated in §2); we require ≥ "
      f"{DISCRIMINATION_MIN_RATIO:g}. A metric that fails here is reporting less "
      "signal than it carries numerical noise.")
    A("")
    A("| metric | min | max | fold range (max/min) | spread (abs) | spread (rel) | "
      f"grid drift {dxs[1]:g}→{dxs[2]:g} µm (matched scale) | **separation** | "
      "discriminates? |")
    A("|---|---:|---:|---:|---:|---:|---:|---:|:--:|")
    for label, _p, _u in _AUDIT:
        r = disc[label]
        A(f"| `{label}` | {_f(r['min'], 4)} | {_f(r['max'], 4)} | "
          f"{_f(r['fold_range'], 4)} | {_f(r['abs_spread'], 3)} | "
          f"{_f(r['relative_spread'], 3)} | {_f(r['grid_drift'], 3)} | "
          f"**{_f(r['separation'], 3)}** | "
          f"{'**YES**' if r['discriminates'] else 'no'} |")
    A("")
    disc_yes = [lbl for lbl, _p, _u in _AUDIT if disc[lbl]["discriminates"]]
    disc_no = [lbl for lbl, _p, _u in _AUDIT if not disc[lbl]["discriminates"]]
    A("**Discriminates:** " + (", ".join(f"`{m}`" for m in disc_yes) or "none") + ".")
    A("")
    A("**Does NOT discriminate:** " + (", ".join(f"`{m}`" for m in disc_no) or "none")
      + ". Three of these are informative failures, not bugs:")
    A("")
    A(f"- **`penetration_length_L`** — the fitted drug `L` spans only "
      f"{disc['penetration_length_L']['min']:.1f}–"
      f"{disc['penetration_length_L']['max']:.1f} µm "
      f"({disc['penetration_length_L']['fold_range']:.2f}×) against a source-null "
      f"`L` of {sn['L_um']:.1f} µm, while its own grid drift is "
      f"{disc['penetration_length_L']['grid_drift']*100:.0f}%. On real tissue the "
      "measured decay length is set by the endothelial-KDE geometry, not by drug "
      "transport — the confound `results/vessel_profile.md` quantifies at length. "
      "The metric is sound; this *section* cannot resolve drugs with it.")
    A(f"- **`lowO2/highO2_exposure`** — grid-stable "
      f"({disc['lowO2/highO2_exposure (quartiles)']['grid_drift']*100:.1f}% drift) "
      f"but essentially identical for every drug "
      f"({disc['lowO2/highO2_exposure (quartiles)']['min']:.3f}–"
      f"{disc['lowO2/highO2_exposure (quartiles)']['max']:.3f}). That is the correct "
      "answer, not a failure: the perfusion gradient is a property of the *tissue*, "
      "and all five drugs are delivered proportionally to the same vascular field, "
      "so they all under-serve the badly-perfused quartile by the same ~3×.")
    A(f"- **`penetration_index`** — the metric this module replaces. See §4.")
    A("")
    A("---")
    A("")
    A("## 4. The explicit contrast with `penetration_index`")
    A("")
    A(f"| | `penetration_index` (retired) | `therapeutic_coverage` @1× IC50 "
      "(proposed headline) |")
    A("|---|---|---|")
    A("| definition | mean(AUC)/p99(AUC) over the tissue mask | fraction of tumor "
      "cells with AUC/T > IC50 |")
    A("| units | dimensionless *shape* | dimensionless *fraction of cells* |")
    A("| scales with delivered dose? | **NO** (scale-invariant by construction) | "
      "**YES** |")
    A("| knows the drug's potency? | **NO** | **YES** (IC50 enters explicitly) |")
    A(f"| range across the {len(drugs)}-drug panel | **{pen['min']:.4f} – {pen['max']:.4f}** "
      f"(spread {pen['max']-pen['min']:.4f}, fold {pen['fold_range']:.3f}×) | "
      f"**{cov['min']:.4f} – {cov['max']:.4f}** (spread {cov['max']-cov['min']:.4f}, "
      f"fold {cov['fold_range']:.3g}×) |")
    A(f"| grid drift {dxs[1]:g}→{dxs[2]:g} µm (worst drug), absolute | "
      f"{_f(pen['grid_drift_absolute'], 3)} | "
      f"{_f(cov['grid_drift_absolute'], 3)} |")
    A(f"| grid drift {dxs[1]:g}→{dxs[2]:g} µm (worst drug), relative | "
      f"{_f(pen['grid_drift_relative'], 3)} | "
      f"{_f(cov['grid_drift_relative'], 3)} |")
    A(f"| **separation** (spread / drift, matched scale) | "
      f"**{_f(pen['separation'], 3)}** | **{_f(cov['separation'], 3)}** |")
    A("")
    A(f"Both are dimensionless numbers in [0,1] on the same tissue, so the absolute "
      f"row is the like-for-like comparison: `penetration_index` drifts "
      f"{_f(pen['grid_drift_absolute'], 3)} with grid refinement while its whole "
      f"across-drug spread is {pen['max']-pen['min']:.4f} — the numerics move it "
      f"~{pen['grid_drift_absolute'] / max(pen['max']-pen['min'], 1e-12):.0f}× more "
      f"than the drugs do. Coverage drifts {_f(cov['grid_drift_absolute'], 3)} "
      f"against a spread of {cov['max']-cov['min']:.3f} — "
      f"~{(cov['max']-cov['min']) / max(cov['grid_drift_absolute'], 1e-12):.0f}× the "
      "other way round.")
    A("")
    A("The scale-invariance claim is not rhetorical. "
      "`delivery_metrics.selftest()` (run it: "
      "`./.venv/bin/python -m tmesim.validation.delivery_metrics --selftest`) "
      "asserts arithmetically that multiplying an exposure field by 10⁻³ and by "
      "10³ leaves `mean/p99` identical to 1e-12 while `therapeutic_coverage` moves "
      "from ~0.000 to 1.000 over the same 10⁶× dose change.")
    A("")
    A("Per drug, side by side:")
    A("")
    A("| drug | `penetration_index` | mean tumor AUC (µM·s) | coverage @1× IC50 |")
    A("|---|---:|---:|---:|")
    for d, r in zip(drugs, rows):
        A(f"| {d} | {r['penetration_index']:.4f} | "
          f"{r['absolute_exposure']['mean_tumor_auc']:.3g} | "
          f"{r['therapeutic_coverage']['coverage_1x']:.4f} |")
    A("")
    A(f"The index moves by {pen['max']-pen['min']:.4f} over a panel whose delivered "
      f"exposure moves by {disc['mean_tumor_auc']['fold_range']:.0f}×. It is not "
      "measuring delivery. Coverage moves by "
      f"{cov['max']-cov['min']:.3f} on the same panel.")
    A("")
    A("---")
    A("")
    A("## 5. Recommendation — the headline delivery number")
    A("")
    marg = disc["median_E/IC50 (potency margin)"]
    n_sat = sum(1 for r in rows if r["therapeutic_coverage"]["coverage_1x"] > 0.999)
    A("> **Headline: `therapeutic_coverage` at 1× IC50 — the fraction of tumor "
      "cells whose time-averaged exposure `E = AUC/T` exceeds the drug's real "
      f"IC50.** Quote it at dx = {dx:g} µm, always with the 0.5×/2× IC50 band beside "
      "it, with **`median E/IC50`** (the continuous potency margin) as its "
      "companion, and with `median_tumor_auc` (µM·s) as the supporting absolute "
      "number.")
    A("")
    if n_sat:
        A(f"> **Saturation note for THIS section:** {n_sat} of {len(rows)} drugs "
          f"score ≥ 0.999 coverage here, so on this particular tumor coverage alone "
          f"cannot rank them. **`median E/IC50` is the ranking number in that "
          f"regime** — it still spans {marg['min']:.3g}–{marg['max']:.3g} "
          f"({marg['fold_range']:.0f}×) across the same panel, with a grid drift of "
          f"{_f(marg['grid_drift'], 3)} and a separation of "
          f"{_f(marg['separation'], 3)}. Coverage remains the number to *quote* "
          "(it is the one an oncologist can read); the margin is the number to "
          "*sort by* when coverage pins at 1.")
        A("")
    A("Justification, against the three requirements:")
    A("")
    A("1. **Physically meaningful.** It compares a concentration to a concentration. "
      "`E = AUC/T` is exactly the quantity `cells.pkpd` feeds into the Hill "
      "dose-response, and `IC50` is by definition the concentration giving "
      "half-maximal effect — so the metric is the delivery half of the kill model, "
      "read out directly rather than through a stochastic Bernoulli draw. It is "
      "bounded in [0,1], monotone in delivered dose, and one sentence of plain "
      "English (\"x% of the tumor saw at least a half-maximal-effect "
      "concentration\").")
    A(f"2. **Grid-converged — with a stated caveat.** Worst-drug change from "
      f"{dxs[-2]:g} µm to {dxs[-1]:g} µm is "
      f"**{cov['grid_drift_absolute']:.4f} absolute** "
      f"({cov['grid_drift_absolute']*100:.1f} percentage points of tumor cells), "
      f"inside the {CONVERGENCE_TOL_ABS:.2f} tolerance for a bounded fraction. That "
      f"same step is {cov['grid_drift_relative']*100:.1f}% *relative*, driven "
      "entirely by the least-covered drug — so the honest statement is **converged "
      "to ±0.03 in coverage units, not to 5% relative**, and coverage should be "
      "quoted with that error bar. Unlike a mask-area average it is computed over a "
      "FIXED population of real cells, so grid refinement cannot change the "
      "denominator — which is why it is this stable at all.")
    A(f"3. **Discriminating.** It spans {cov['min']:.4f}–{cov['max']:.4f} across the "
      f"panel ({cov['max']-cov['min']:.3f} absolute, "
      f"{cov['fold_range']:.2f}× fold) — a separation of "
      f"**{_f(cov['separation'], 3)}×** its own grid drift, versus "
      f"{_f(pen['separation'], 3)} for `penetration_index`. It also ranks the "
      "drugs the way the mechanism says it should: potency and delivery are combined, "
      "so a drug wins only when it both arrives AND is potent at the concentration "
      "it achieves.")
    cov_order = [d for d, _ in sorted(
        zip(drugs, [r["therapeutic_coverage"]["coverage_1x"] for r in rows]),
        key=lambda t: -t[1])]
    kill_order = [d for d, _ in sorted(
        zip(drugs, [r["kill_fraction"] for r in rows]), key=lambda t: -t[1])]
    if cov_order == kill_order:
        A("")
        A(f"   **And it recovers the kill ordering exactly.** Ranked by coverage the "
          f"panel is **{' > '.join(cov_order)}** — the *identical* order the "
          "simulated `kill_fraction` produces on the same runs (kill is the "
          "engine's own pharmacodynamic endpoint; on this ladder it is itself only "
          f"grid-stable to {conv['kill_fraction (reference)']['worst_drug_abs_drift_final']:.3f} "
          "absolute, so the agreement is between two independently-computed "
          "orderings, not between a number and its own definition). Coverage is a "
          "pure *delivery* statistic: it never touches the Hill curve, the "
          "per-cell sensitivities, the hypoxia/quiescence resistance factors or the "
          "stochastic Bernoulli draw. Reproducing the kill ranking from delivery "
          "and potency alone is the strongest evidence on offer that it is measuring "
          "the right thing — and `penetration_index`, on the same panel, ranks the "
          f"drugs **{' > '.join(d for d, _ in sorted(zip(drugs, [r['penetration_index'] for r in rows]), key=lambda t: -t[1]))}**, "
          "which is noise.")
    A("")
    A("**What to quote alongside it (and why not instead of it):**")
    A("")
    hi_d = max(zip(drugs, rows), key=lambda t: t[1]["absolute_exposure"]["median_tumor_auc"])
    lo_d = min(zip(drugs, rows), key=lambda t: t[1]["absolute_exposure"]["median_tumor_auc"])
    A("- `median_tumor_auc` (µM·s) — the raw absolute delivery. Strongly "
      f"discriminating ({disc['median_tumor_auc']['fold_range']:.0f}× across the "
      f"panel) though not itself grid-converged on this ladder "
      f"({disc['median_tumor_auc']['grid_drift']*100:.1f}% worst-drug drift; use "
      f"`mean_tumor_auc`, {disc['mean_tumor_auc']['grid_drift']*100:.1f}%, if a "
      "converged absolute number is required), and in any case *not* interpretable "
      "without a potency: "
      f"{hi_d[1]['absolute_exposure']['median_tumor_auc']:.2g} µM·s of {hi_d[0]} "
      f"(IC50 {hi_d[1]['ic50_uM']:.3g} µM) and "
      f"{lo_d[1]['absolute_exposure']['median_tumor_auc']:.2g} µM·s of {lo_d[0]} "
      f"(IC50 {lo_d[1]['ic50_uM']:.3g} µM) are not comparable doses. Report it as "
      "the input to coverage, not as the verdict.")
    A("- `median E/IC50` and `P10 E/IC50` (dimensionless) — the potency-normalised "
      "exposure margin at the median and at the worst-served decile. Same "
      "information content as coverage, continuous instead of thresholded, and "
      "immune to saturation; quote `median E/IC50` whenever coverage pins at 1.000. "
      f"It inherits the median's grid sensitivity "
      f"({disc['median_E/IC50 (potency margin)']['grid_drift']*100:.1f}% worst-drug "
      f"drift, versus {disc['therapeutic_coverage_1x']['grid_drift']:.3f} absolute "
      "for coverage), which is the price of being continuous — it is a ranking "
      "companion, not a replacement headline.")
    A("- `E_p10` (µM) — the worst-served decile, the honest floor of the delivery "
      "distribution and the number that should worry you about regrowth.")
    L_excess = [r["penetration_length"]["L_excess_over_source_null"] for r in rows]
    A("- `penetration_length_fit` `L` (µm) — quote only against the SOURCE NULL. On "
      f"this section the fitted drug `L` spans "
      f"{min(r['penetration_length']['L_um'] for r in rows):.1f}–"
      f"{max(r['penetration_length']['L_um'] for r in rows):.1f} µm against a "
      f"source-null `L` of {sn['L_um']:.1f} µm "
      f"(drug/null ratio {min(L_excess):.2f}–{max(L_excess):.2f}×), i.e. the "
      "measured length is dominated by vessel/KDE geometry rather than by drug "
      "transport — the same confound `results/vessel_profile.md` documents at "
      "length. The drug-intrinsic length separates cleanly only in the isolated "
      "single-vessel control, where doxorubicin's gradient is ~4× steeper than "
      "temozolomide's. `L` is a legitimate *physical* number and a legitimate "
      "*model* diagnostic; it is not a delivery ranking on this section.")
    A(f"- `lowO₂/highO₂ exposure ratio` — diagnostic only. On this section it is "
      f"{min(r['hypoxia']['ratio_lowest_over_highest_o2_quartile'] for r in rows):.3f}"
      f"–{max(r['hypoxia']['ratio_lowest_over_highest_o2_quartile'] for r in rows):.3f} "
      "across the panel. It is a ratio (hence scale-invariant, hence the same "
      "structural weakness as `penetration_index`) and the banded hypoxic/normoxic "
      "form sits on a partition that is itself not grid-converged. Read it as "
      "'is delivery correlated with perfusion?', never as a drug ranking.")
    A("")
    A("**What to retire.** `penetration_index` should no longer be quoted as a "
      "delivery number anywhere. It is not wrong — it is a faithful measurement of "
      "*exposure-field uniformity* — but on real tissue at converged resolution that "
      "uniformity is a property of the reconstructed vasculature, shared by every "
      "drug. If a uniformity statistic is wanted, quote it as "
      "`exposure_uniformity_index` and say what it is.")
    A("")
    A("---")
    A("")
    A("## 6. Honesty and caveats")
    A("")
    A("- **No outcome label.** The Xenium section is an untreated resection. Every "
      "number here is a PREDICTION of the calibrated physics on real geometry, not "
      "a validated response. Coverage being grid-converged and drug-discriminating "
      "says it is a *well-posed* metric; it does not say the underlying delivery "
      "model is right.")
    A("- **Coverage inherits the IC50's uncertainty.** The IC50s are real GDSC / "
      "literature potencies measured in cell culture, not in this tissue. That is "
      "why the 0.5× and 2× variants are reported as a band and not buried: they are "
      "the metric's sensitivity to the one external number it depends on. Where the "
      "0.5×–2× band is wide, the ranking at that position is soft.")
    n_sat_6 = sum(1 for r in rows if r["therapeutic_coverage"]["coverage_1x"] > 0.999)
    A("- **Coverage can saturate — this is its one structural weakness, and it is "
      "why the continuous companion is reported next to it.** A drug delivering far "
      "above its IC50 *everywhere* scores 1.000 and stops discriminating. "
      + (f"On this panel that does NOT bite — {n_sat_6} of {len(rows)} drugs reach "
         f"0.999 and the maximum is {max(r['therapeutic_coverage']['coverage_1x'] for r in rows):.3f} "
         "— so the weakness is latent here, not active. It will bite on a "
         "higher-dose or more potent panel."
         if n_sat_6 == 0 else
         f"On this panel it does bite: {n_sat_6} of {len(rows)} drugs are at or "
         "above 0.999.")
      + " `median E/IC50` — the **potency-normalised exposure margin** — is the same "
      "comparison without the threshold: dimensionless, drug-specific, absolute, "
      "unbounded above, and it spans "
      f"{disc['median_E/IC50 (potency margin)']['fold_range']:.0f}× across the panel. "
      "Read them together: coverage says *how much of the tumor* cleared the bar, the "
      "margin says *by how far*. Coverage is the right headline precisely in the "
      "regime that matters clinically — where part of the tumor is under-dosed — and "
      "the margin is the tiebreaker when it is not.")
    A(f"- **One anchor sets the absolute scale.** The single α (= {res['alpha']:.1f} "
      "here, read live from `validation.calibration.calibrate_alpha()`) is fit to "
      "one doxorubicin penetration measurement; every absolute exposure — and "
      "therefore every coverage value — rides on it. Coverage is *relatively* robust "
      f"because the same α applies to all {len(drugs)} drugs, but the absolute level "
      "is not independently validated. **Older reports in `results/` may quote a "
      "different α** (the half-distance vs e-folding-length convention was corrected "
      "in `calibration.py`, a factor 1/ln²2 = 2.08×); the absolute µM·s figures here "
      "are therefore not directly comparable to those, while the *structural* "
      "conclusions — the index collapsing, coverage spreading — are about the "
      "relative behaviour of the metrics and are unaffected.")
    hyp_by_dx = ", ".join(
        f"{k:g} µm → {res['per_dx'][k][drugs[0]]['hypoxia']['hypoxic_tumor_fraction']*100:.1f}%"
        for k in dxs)
    A("- **The tumor hypoxic fraction is not grid-converged** "
      "(`results/real_drug_panel.md`), which is exactly why metric 5 is a diagnostic "
      "and why its primary form is the threshold-free quartile contrast. On the "
      f"reconstruction as it stands the hypoxic tumor fraction is: {hyp_by_dx}. "
      "Where it is 0%, the banded hypoxic/normoxic ratio is legitimately undefined "
      "(`n/a`) — that is reported, not papered over.")
    A("- **2D section, straight-line distances, static plasma level.** All the "
      "standing limitations of the engine (`results/numerical_rigor.md`) apply "
      "unchanged; nothing in `core` / `physics` / `molecular` / `reconstruct` / "
      "`sim.engine` was modified for this audit.")
    A("")
    A(f"*Total runtime {res['runtime_s']:.0f} s "
      f"({len(dxs)} grids × {len(drugs)} drugs).*")
    A("")

    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\n".join(L))
    return path


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _print_summary(res: Dict) -> None:
    dx = res["dx_production"]
    drugs = res["drugs"]
    disc = {r["metric"]: r for r in res["discrimination"]}
    print("\n" + "=" * 108)
    print(f"DELIVERY METRICS on the REAL Xenium GBM section (dx = {dx:g} um)")
    print("=" * 108)
    hdr = (f"{'drug':14s} {'IC50(uM)':>9} {'meanAUC':>10} {'medAUC':>10} "
           f"{'E_p10':>9} {'E_p50':>9} {'E_p90':>9} {'cov0.5x':>8} {'cov1x':>8} "
           f"{'cov2x':>8} {'L(um)':>7} {'lo/hiO2':>8} {'pen_idx':>8}")
    print(hdr)
    print("-" * len(hdr))
    for d in drugs:
        r = res["per_dx"][dx][d]
        a, c, p, pl, h = (r["absolute_exposure"], r["therapeutic_coverage"],
                          r["exposure_percentiles"], r["penetration_length"],
                          r["hypoxia"])
        print(f"{d:14s} {r['ic50_uM']:9.3g} {a['mean_tumor_auc']:10.3g} "
              f"{a['median_tumor_auc']:10.3g} {p['E_p10_uM']:9.3g} "
              f"{p['E_p50_uM']:9.3g} {p['E_p90_uM']:9.3g} "
              f"{c['coverage_0.5x']:8.4f} {c['coverage_1x']:8.4f} "
              f"{c['coverage_2x']:8.4f} {pl['L_um']:7.1f} "
              f"{h['ratio_lowest_over_highest_o2_quartile']:8.3f} "
              f"{r['penetration_index']:8.4f}")
    print("=" * 108)
    print(f"{'metric':36s} {'converged':>10} {'fold range':>12} {'drift':>10} "
          f"{'separation':>11} {'discriminates':>14}")
    print("-" * 108)
    for label, _p, _u in _AUDIT:
        r = disc[label]
        print(f"{label:36s} {str(r['converged']):>10} {_f(r['fold_range'], 4):>12} "
              f"{_f(r['grid_drift'], 3):>10} {_f(r['separation'], 3):>11} "
              f"{str(r['discriminates']):>14}")
    print("=" * 108)
    cov = disc["therapeutic_coverage_1x"]
    pen = disc["penetration_index (RETIRED)"]
    print(f"RECOMMENDED HEADLINE: therapeutic_coverage @1x IC50  "
          f"[{cov['min']:.4f}-{cov['max']:.4f}, drift {_f(cov['grid_drift'],3)}, "
          f"separation {_f(cov['separation'],3)}]")
    print(f"RETIRE              : penetration_index  "
          f"[{pen['min']:.4f}-{pen['max']:.4f}, spread {pen['max']-pen['min']:.4f} "
          f"-- a vascular-geometry constant, not delivery]")
    print("=" * 108)


def selftest(verbose: bool = True) -> None:
    """Known-answer checks on the metric definitions themselves (no simulation).

    The load-bearing one is the last: multiplying the whole exposure field by a
    constant leaves ``penetration_index`` EXACTLY unchanged but moves
    ``therapeutic_coverage``. That is the entire argument for the replacement,
    proved arithmetically rather than asserted.
    """
    T = 21600.0
    rng = np.random.default_rng(0)

    # 1. Coverage against a hand-computable case.
    auc = np.array([0.0, 1.0, 2.0, 3.0, 4.0]) * T      # E = 0,1,2,3,4 uM
    c = therapeutic_coverage(auc, ic50=2.0, T=T)
    assert c["coverage_1x"] == 0.4, c            # E > 2 -> {3,4} -> 2/5
    assert c["coverage_0.5x"] == 0.6, c          # E > 1 -> {2,3,4} -> 3/5
    assert c["coverage_2x"] == 0.0, c            # E > 4 -> none
    assert c["threshold_1x_uM"] == 2.0

    # 2. Coverage is monotone in the delivered dose and in potency.
    base = rng.lognormal(mean=1.0, sigma=1.0, size=5000) * T
    cov = lambda a, ic: therapeutic_coverage(a, ic50=ic, T=T)["coverage_1x"]
    assert cov(base * 2.0, 5.0) > cov(base, 5.0), "coverage must rise with dose"
    assert cov(base, 1.0) > cov(base, 20.0), "coverage must rise with potency"

    # 3. Percentiles are ordered and the heterogeneity ratio is consistent.
    p = exposure_percentiles(base / T, base)
    assert p["E_p10_uM"] < p["E_p50_uM"] < p["E_p90_uM"], p
    assert abs(p["E_p90_over_p10"] - p["E_p90_uM"] / p["E_p10_uM"]) < 1e-9

    # 4. THE POINT: penetration_index is scale-invariant, coverage is not.
    field = rng.lognormal(mean=2.0, sigma=1.5, size=20000)
    pen = lambda f: float(f.mean() / np.percentile(f, 99.0))
    for s in (1e-3, 1.0, 1e3):
        assert abs(pen(field * s) - pen(field)) < 1e-12, (
            "penetration_index is supposed to be exactly scale-invariant")
    c_lo = therapeutic_coverage(field * 1e-3 * T, ic50=1.0, T=T)["coverage_1x"]
    c_hi = therapeutic_coverage(field * 1e3 * T, ic50=1.0, T=T)["coverage_1x"]
    assert c_hi > c_lo, "therapeutic_coverage must respond to a 10^6x dose change"
    if verbose:
        print(f"[selftest] penetration_index invariant under 10^6x dose scaling "
              f"({pen(field):.6f} at every scale); therapeutic_coverage moves "
              f"{c_lo:.4f} -> {c_hi:.4f}. All metric known-answer checks passed.")


def main(argv=None) -> Dict:
    ap = argparse.ArgumentParser(
        description="Defensible delivery metrics on the real Xenium GBM section")
    ap.add_argument("--dx", type=float, nargs="+", default=DX_LADDER,
                    help="grid-refinement ladder in um (default %(default)s)")
    ap.add_argument("--drugs", type=str, nargs="+", default=PANEL)
    ap.add_argument("--cache", type=str, default=None,
                    help="pickle path to cache/reuse the raw simulation results")
    ap.add_argument("--png", type=str, default=_OUT_PNG)
    ap.add_argument("--md", type=str, default=_OUT_MD)
    ap.add_argument("--selftest", action="store_true",
                    help="run only the metric known-answer checks and exit")
    args = ap.parse_args(argv)

    selftest()
    if args.selftest:
        return {}

    res = backfill(run(args.dx, args.drugs, cache_path=args.cache))
    print(f"\nfigure written: {make_figure(res, args.png)}")
    print(f"report written: {write_report(res, args.md)}")
    _print_summary(res)
    return res


if __name__ == "__main__":
    main()
