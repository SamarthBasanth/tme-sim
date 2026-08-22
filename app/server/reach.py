"""THE SCALE-FREE REACH MAP — how much of the tumour the drug never reaches.

Why this module exists
----------------------
Codon has stopped ranking compounds. The cross-compound ordering was measured
and does not hold (``results/pk_ranking_robustness.md``,
``results/capacity_ranking_robustness.md``), so the product is the SINGLE-compound
map: one real tissue, one compound, where the drug reaches and how much of the
tumour it never reaches.

The replacement headline therefore has to be **SCALE-FREE**, and that constraint
is the entire design.
``validation.delivery_metrics.therapeutic_coverage`` cannot be it. Its own
docstring lists "(d) **absolute** — scaling the field scales the answer" as a
property, and the two quantities it reads — the IC50 and the plasma
concentration — are exactly the two the robustness studies measured as uncertain
(per-compound plasma exposure across the retained GBM compounds spans
**113.1x**, ``results/pk_ranking_robustness.md`` sec. 0). Shipping an absolute
number as the headline would rebuild the same failure one compound at a time.

What IS protected is a **RATIO BETWEEN TWO POINTS OF ONE MAP**, and only to the
extent stated here. Two perturbations were swept and measured:

* **Plasma / dose scale — EXACT, for the finite-window AUC field.** The delivery
  source is one multiplication by the plasma level and the transport PDE is
  linear in ``C``, so ``AUC`` is exactly proportional to it. Measured at the
  shipping configuration (real Xenium GBM section, dx = 12 um, dt = 10 s,
  convection ON, 6 h window) over plasma factors 0.005 to 26 — a 5,200x span —
  the worst relative deviation of the AUC field from an exact rescale anywhere
  was **1.41e-14** (``results/pk_ranking_robustness.md`` sec. 1.2). A factor that
  multiplies every point of the map divides out of a ratio between points.
* **Tissue capacity — EXACT for a CONVERGED gradient shape, NOT for a finite
  window.** ``results/capacity_ranking_robustness.md`` sec. 2b states it at that
  strength and no more: *"steady-state and equilibrium-referenced outputs —
  penetration length, interstitial:plasma ratio, converged gradient shape — are
  EXACTLY invariant to tissue capacity (verified to 1e-9). The finite-window kill
  ranking is not."* A capacity term ``R`` enters as ``dC/dt = (1/R)[...]``, i.e.
  a pure rescale of TIME, so ``C_R(x, t) = C_1(x, t/R)`` and over a finite window
  ``AUC_R(x, T) = R * AUC_1(x, T/R)``. The spatial shape is then the shape at a
  DIFFERENT effective window: the between-point ratio is exactly
  capacity-invariant only where the field has converged, and carries a residual
  where it has not. This product's window is finite (6 h by default), so the
  reach ratio sits in the "converged gradient shape" family only to the extent
  the field has converged, and that is not asserted here.

**Neither invariance protects against error in the transport law itself**, whose
descriptor path is unscored (``results/descriptor_falsification.md``). And the
two bullets above are a LIST OF THE AXES THAT WERE SWEPT, not a coverage claim:
silence about an unswept axis here is silence, not a pass.

The quantities (exact definitions)
----------------------------------
Let ``AUC(x)`` be the transport solve's exposure field (uM*s).

* **Perivascular reference exposure** ``E_ref`` (uM*s) — the mean AUC over tissue
  in the FIRST usable bin of the distance-from-nearest-endothelial-cell profile.
  Bins are ``vessel_profile.BIN_UM`` = 10 um wide and a bin is usable only with at
  least ``vessel_profile.MIN_BIN_COUNT`` = 30 grid points in it. This is ALREADY
  computed as ``delivery_metrics.penetration_length_fit(...)["drug"]
  ["perivascular_mean"]`` and that value is REUSED here verbatim — never
  recomputed — so the reach number and the published penetration length can
  never drift apart.
* **Reach ratio** ``r_i = AUC_i / E_ref`` for tumour cell ``i``. Dimensionless.
  Numerator and denominator are both linear in any plasma scale factor, so it is
  invariant to dose; both are read off the same field at the same window, so a
  capacity term cancels from it exactly wherever that field has converged (see
  above for where it does not).
* **Under-reached fraction** ``F(x)`` — the fraction of VIABLE tumour cells with
  ``r_i < x``. "Viable" is defined EXACTLY as
  ``delivery_metrics.hypoxic_vs_normoxic_exposure_ratio`` defines it — local
  oxygen ``o >= 2 * OxygenConfig.hypoxia_threshold`` at the cell's own grid
  point, the engine's own banding (``tmesim.cells.agents.update_cell_states``) —
  so the two cannot disagree about which cells they are talking about.
* **The curve** — ``F(x)`` over the fixed log ladder :data:`REACH_LADDER`.
* **The headline** — ``F(0.1)``, picked off that ladder, not recomputed.
* **The scale-free falloff** — the same distance bins that produced ``E_ref``,
  each bin's mean AUC divided by ``E_ref``: exposure vs distance from the nearest
  vessel, normalised, with the bin distances and per-bin counts.

**The 0.1 threshold is a REPORTING CONVENTION, not a physical constant.** It is
returned in the payload as :data:`CONVENTION_NOTE` so a reader meets the sentence
where they meet the number, and the whole ladder is returned beside the headline
so the choice stays visible. That is deliberate: ``1x IC50`` acquired an
authority it had not earned by being quoted alone, and ``1/10`` must not repeat
it.

REFUSAL, not defaulting
-----------------------
Whenever a quantity is undefined — no tumour cells, no viable tumour cells, no
usable distance bin, ``E_ref`` non-finite or non-positive — this returns
``null_reason`` and **no number**. Nothing is substituted, extrapolated or
imputed. This is the repo's standing convention (see
``vessel_profile.D_HALF_NEVER`` for the same rule applied to a half-distance
that does not exist).

Nothing in ``tmesim/`` is modified by this module; it only consumes a
``Microenvironment``, a ``SimulationResult`` and the already-computed
``penetration_length_fit`` block.

Units: AUC / E_ref uM*s, distance um, reach ratio and F(x) dimensionless.
"""
from __future__ import annotations

from typing import Dict, Optional, Sequence

import numpy as np

from tmesim.core import Microenvironment, OxygenConfig, SimulationResult
from tmesim.validation.delivery_metrics import (
    penetration_length_fit,
    tumor_cell_exposure,
)
# Imported deliberately rather than re-implemented: this is the same
# nearest-grid-cell sampling map (``domain.to_index``) that
# ``hypoxic_vs_normoxic_exposure_ratio`` uses to read a tumour cell's local
# oxygen, and the "viable" band here has to be the identical set of cells.
from tmesim.validation.delivery_metrics import _sample_at_cells
from tmesim.validation.vessel_profile import (
    BIN_UM,
    MIN_BIN_COUNT,
    distance_profile,
    usable_bins,
)

__all__ = [
    "REACH_LADDER",
    "HEADLINE_X",
    "CONVENTION_NOTE",
    "SCALE_INVARIANCE_NOTE",
    "compute_reach",
    "null_reach",
]

#: The fixed log ladder of reach ratios at which ``F(x)`` is reported. Returned
#: in full, always, so the headline threshold is visibly one point on a curve
#: rather than a constant of nature.
REACH_LADDER: Sequence[float] = (0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0)

#: The headline point of the ladder. It is a member of :data:`REACH_LADDER` and
#: the headline is read OFF the ladder, never recomputed.
HEADLINE_X: float = 0.1

#: Shipped with the number, every time. The threshold is a choice.
CONVENTION_NOTE: str = (
    "x = 0.1 is a REPORTING CONVENTION, not a physical constant. Nothing in the "
    "transport physics, the tissue or the compound marks it out, and no "
    "measurement in this repo selects it. It is NO LONGER HEADLINED ANYWHERE, "
    "and the reason is worth stating because it is a result rather than a "
    "preference: measured on the flagship Xenium GBM section at dx = 12 um, "
    "F(0.1) is 0.0005 — the compound reaches essentially the whole tumour, and "
    "a headline that reads 0.05% on the specimen this product ships with is not "
    "a headline. Retuning the threshold to a value that reads better was "
    "declined: any threshold chosen AFTER seeing what the field does is "
    "post-hoc selection wearing a convention's clothes, which is exactly how "
    "'1x IC50' came to be defended. The headline is instead the delivery spread "
    "(P90/P10 of tumour-cell exposure), which needs no threshold at all and "
    "carries signal even where the map is flat. The full F(x) ladder is still "
    "returned, and is still the evidence — read the curve, not any one point "
    "on it."
)

#: What the ratio IS protected against, at the strength each was measured, and
#: what it is not protected against. Returned in the payload.
SCALE_INVARIANCE_NOTE: str = (
    "The reach ratio is a ratio between two points of ONE map, so a factor that "
    "multiplies every point of the map divides out of it. Two axes were swept. "
    "(1) PLASMA/DOSE SCALE: exact for the finite-window AUC field — worst "
    "relative deviation from an exact rescale 1.41e-14 over plasma factors "
    "0.005 to 26 (a 5,200x span), measured at the shipping configuration "
    "(results/pk_ranking_robustness.md sec. 1.2). "
    "(2) TISSUE CAPACITY: the physics of the cited invariance is CORRECT and the "
    "REGIME CLAIM WAS NOT — read that distinction before reading the source. "
    "capacity_ranking_robustness.md sec. 2b says capacity is exact for a "
    "CONVERGED gradient; that is true. What was wrong, in an earlier version of "
    "this very note, was the assumption that this product's configuration IS a "
    "converged gradient. It is not, and the next reader of that document will "
    "make the same inference unless the gap between its claim and our window is "
    "written down. So what follows is NOT an invariance — it is a BOUNDED "
    "RESIDUAL, "
    "and the bound differs by an order of magnitude between the two readings "
    "below. A capacity term acts as a rescale of time "
    "(AUC_R(x,T) = R * AUC_1(x,T/R)), so it is exact only for a CONVERGED "
    "gradient shape (1e-9, results/capacity_ranking_robustness.md sec. 2b) and "
    "THIS WINDOW IS NOT CONVERGED: measured on the flagship Xenium GBM section, "
    "the gap to the infinite-window limit is 0.01302 for temozolomide at the "
    "shipping 6 h and decays as 281/T with no plateau, so reaching 0.001 would "
    "need about 78 h. Consequences, measured: the UNDER-REACHED LADDER moves by "
    "at most 0.0011 over capacity R = 1 to 39 and is flat to R = 60, so it is "
    "usable with the bound stated as R <~ 60 — but over the full measured "
    "capacity range (3.2 to 4184 mL/g, n = 33) it moves up to 0.80, which is "
    "WORSE than the compound-ordering failure this product withdrew. The "
    "NORMALISED FALLOFF moves 0.3772 absolute at R = 39 on a curve whose whole "
    "range is 1.0 to 0.51, with the fitted penetration length collapsing 279.7 "
    "to 76.7 um: the SHAPE claim inherits the same disease as the ranking and "
    "MUST NOT be presented as capacity-invariant. "
    "(3) TRANSVASCULAR PENETRANCE: exactly invariant, verified bit-identical at "
    "bbb_penetrance 0.1 and 1.0, because it is a pure global multiplier on the "
    "source. That is a LIMIT AS MUCH AS A PROPERTY - this map cannot see a "
    "permeability barrier at all, and depends on the compound only through "
    "L = sqrt(D_eff/k). Do not read it as a barrier or brain-penetration "
    "finding. "
    "Neither invariance protects against error in the transport law itself, "
    "whose descriptor path is unscored (results/descriptor_falsification.md). "
    "This is the list of axes that were SWEPT, not a coverage claim: an axis "
    "not named here was not tested."
)

#: How ``E_ref`` is obtained. Recorded in the payload so a consumer can see that
#: the reach number and the penetration length read the same bin.
E_REF_SOURCE: str = (
    "delivery_metrics.penetration_length_fit(...)['drug']['perivascular_mean'] "
    f"— mean AUC over tissue in the first distance bin ({BIN_UM:g} um wide, "
    f">= {MIN_BIN_COUNT} grid points) of the distance-from-nearest-"
    "endothelial-cell profile, REUSED verbatim, never recomputed"
)

#: How "viable" is defined, in the payload, in the words of the metric it copies.
VIABLE_DEFINITION: str = (
    "a TUMOR cell whose local oxygen at its own grid point satisfies "
    "o >= 2 * OxygenConfig.hypoxia_threshold — the engine's own banding "
    "(tmesim.cells.agents.update_cell_states), identical to the 'viable' band of "
    "delivery_metrics.hypoxic_vs_normoxic_exposure_ratio, so the two cannot "
    "disagree about which cells they describe"
)


def null_reach(reason: str, **known) -> Dict:
    """A refusal. Every quantity that could not be computed is absent, not
    defaulted; anything that genuinely IS known (cell counts, the convention
    text) is still reported."""
    out: Dict = {
        "null_reason": reason,
        "headline_under_reached_fraction": None,
        "headline_x": float(HEADLINE_X),
        "curve": None,
        "falloff": None,
        "E_ref_uM_s": None,
        "E_ref_source": E_REF_SOURCE,
        "n_viable": None,
        "n_tumor": None,
        "viable_definition": VIABLE_DEFINITION,
        "hypoxia_threshold": None,
        "convention": CONVENTION_NOTE,
        "scale_invariance": SCALE_INVARIANCE_NOTE,
        "penetration_fit_error": None,
        "units": _UNITS,
    }
    out.update(known)
    return out


_UNITS = {
    "E_ref": "uM*s",
    "reach_ratio": "dimensionless (AUC_i / E_ref)",
    "under_reached_fraction": "fraction of VIABLE tumour cells, in [0, 1]",
    "distance_um": "um",
    "normalised_mean_auc": "dimensionless (bin mean AUC / E_ref)",
}


def compute_reach(
    micro: Microenvironment,
    sim_result: SimulationResult,
    distance_um: Optional[np.ndarray] = None,
    penetration: Optional[Dict] = None,
    ladder: Sequence[float] = REACH_LADDER,
    hypoxia_threshold: float = OxygenConfig().hypoxia_threshold,
) -> Dict:
    """The scale-free reach map for one (tissue, simulation) pair.

    Parameters
    ----------
    micro : Microenvironment
        The reconstructed real microenvironment the simulation ran on.
    sim_result : SimulationResult
        Output of :func:`tmesim.sim.engine.simulate_response`; must carry
        ``fields['auc']`` and ``fields['oxygen']``.
    distance_um : np.ndarray, optional
        Pre-computed distance-from-nearest-vessel map (um). The distance
        transform is the expensive part of this computation and it is pure
        geometry, so a caller that has already built one (for
        :func:`delivery_metrics.compute_all`, say) should pass the SAME array
        here rather than paying for it twice.
    penetration : dict, optional
        An already-computed :func:`delivery_metrics.penetration_length_fit`
        result. Pass it when the caller has one — ``E_ref`` is read out of
        ``penetration['drug']['perivascular_mean']`` either way, so passing it
        only avoids re-running the fit. When omitted the fit is run here with
        ``include_source_null=False`` (the source null plays no part in this
        metric and would double the binning work); ``['drug']
        ['perivascular_mean']`` is the identical number either way.
    ladder : sequence of float
        Reach ratios at which ``F`` is evaluated. Defaults to
        :data:`REACH_LADDER`; the headline is picked off it.
    hypoxia_threshold : float
        Passed through to the "viable" band, defaulting to
        ``OxygenConfig.hypoxia_threshold`` — the same default
        :func:`delivery_metrics.hypoxic_vs_normoxic_exposure_ratio` uses.

    Returns
    -------
    dict
        ``headline_under_reached_fraction`` (``F(0.1)``), ``headline_x``,
        ``curve`` (``{'x': [...], 'under_reached_fraction': [...]}``),
        ``falloff`` (``{'distance_um', 'normalised_mean_auc', 'count',
        'bin_um', 'min_bin_count'}``), ``E_ref_uM_s`` + ``E_ref_source``,
        ``n_viable``, ``n_tumor``, ``viable_definition``, ``convention``,
        ``scale_invariance``, ``penetration_fit_error``, ``units`` and
        ``null_reason``.

        ``null_reason`` is ``None`` on success and a sentence otherwise, and
        when it is set every derived number is ``None``. See the module
        docstring: this refuses rather than defaults.

    Note on ``penetration_fit_error``
    ---------------------------------
    A failure of the exponential-plus-floor fit is REPORTED (as
    ``penetration_fit_error``) but does not null this result, because no
    quantity here reads a fitted parameter: ``E_ref`` is a bin mean, and the
    ladder and the falloff are ratios of measured bin/cell values. Nulling on it
    would be discarding a number that is defined, which is the opposite of the
    refusal rule.
    """
    x = np.asarray(list(ladder), dtype=np.float64)

    # --- 1. the tumour cells and their exposure ---------------------------- #
    exp = tumor_cell_exposure(micro, sim_result)
    n_tumor = int(exp["n_tumor"])
    auc_cells = np.asarray(exp["auc_uM_s"], dtype=np.float64)
    if n_tumor == 0:
        return null_reach(
            "no TUMOR cells in this reconstruction, so there is no tumour for "
            "the drug to reach or miss", n_tumor=0)

    # --- 2. the viable band, defined exactly as metric 5 defines it --------- #
    oxygen = sim_result.fields.get("oxygen")
    if oxygen is None:
        return null_reach(
            "sim_result.fields carries no 'oxygen' field, so the VIABLE band "
            "(o >= 2 * hypoxia_threshold) that this fraction is taken over "
            "cannot be identified; no substitute population is used",
            n_tumor=n_tumor)
    o2 = _sample_at_cells(oxygen, micro, exp["tumor_mask"])
    viable = np.asarray(o2, dtype=np.float64) >= 2.0 * float(hypoxia_threshold)
    n_viable = int(viable.sum())
    if n_viable == 0:
        return null_reach(
            "no VIABLE tumour cells in this reconstruction (no tumour cell has "
            f"local oxygen >= 2 * hypoxia_threshold = "
            f"{2.0 * float(hypoxia_threshold):g}), so a fraction OF the viable "
            "population is undefined; it is not fallen back to all tumour cells",
            n_tumor=n_tumor, n_viable=0,
            hypoxia_threshold=float(hypoxia_threshold))

    # --- 3. the perivascular reference, REUSED not recomputed --------------- #
    # The distance transform is the expensive step and it is pure geometry, so
    # it is resolved ONCE here and handed to both the fit and the falloff.
    if distance_um is None:
        distance_um = _distance_for(micro)
    if penetration is None:
        penetration = penetration_length_fit(
            sim_result.fields["auc"], micro, distance_um=distance_um,
            include_source_null=False)
    fit_error = penetration["drug"].get("error")
    e_ref = float(penetration["drug"]["perivascular_mean"])
    if not np.isfinite(e_ref) or e_ref <= 0.0:
        return null_reach(
            f"the perivascular reference exposure E_ref is {e_ref!r}: the "
            "distance-from-vessel profile has no usable first bin (a bin needs "
            f">= {MIN_BIN_COUNT} grid points) or its mean AUC is not positive, "
            "so every reach ratio would be a division by it",
            n_tumor=n_tumor, n_viable=n_viable,
            hypoxia_threshold=float(hypoxia_threshold),
            penetration_fit_error=fit_error)

    # --- 4. the curve ------------------------------------------------------- #
    # F(x) = fraction of VIABLE tumour cells with r_i < x. Evaluated by binary
    # search on the sorted ratios, which makes monotonicity in x exact rather
    # than incidental, and keeps the strict "<" of the definition.
    r = auc_cells[viable] / e_ref
    r_sorted = np.sort(r)
    counts = np.searchsorted(r_sorted, x, side="left")
    F = counts.astype(np.float64) / float(n_viable)

    on_ladder = np.flatnonzero(x == float(HEADLINE_X))
    if on_ladder.size == 0:
        return null_reach(
            f"the headline threshold x = {HEADLINE_X:g} is not on the ladder "
            f"{[float(v) for v in x]} this call was given, and it is read OFF "
            "the ladder rather than recomputed, so no headline is reported",
            n_tumor=n_tumor, n_viable=n_viable, E_ref_uM_s=e_ref,
            hypoxia_threshold=float(hypoxia_threshold),
            penetration_fit_error=fit_error)
    headline_idx = int(on_ladder[0])

    # --- 5. the scale-free falloff ----------------------------------------- #
    # The SAME bins that produced E_ref: same distance map, same binning, same
    # usable-bin rule, so falloff[0] is 1.0 by construction.
    tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
    prof = distance_profile(np.asarray(sim_result.fields["auc"],
                                       dtype=np.float64),
                            distance_um, tissue)
    d_bins, m_bins, n_bins = usable_bins(prof)
    if m_bins.size == 0 or float(m_bins[0]) != e_ref:
        # The falloff must be the same profile E_ref came from. If it is not,
        # something upstream has changed the binning between the two calls and
        # the honest answer is to say so, not to publish two profiles that
        # disagree.
        return null_reach(
            "the distance profile rebuilt here does not reproduce the "
            "perivascular mean that delivery_metrics.penetration_length_fit "
            f"reported (E_ref = {e_ref!r}); the falloff and the reference "
            "would come from different binnings, so neither is reported",
            n_tumor=n_tumor, n_viable=n_viable, E_ref_uM_s=e_ref,
            hypoxia_threshold=float(hypoxia_threshold),
            penetration_fit_error=fit_error)

    falloff = m_bins / e_ref

    return {
        "null_reason": None,
        "headline_under_reached_fraction": float(F[headline_idx]),
        "headline_x": float(HEADLINE_X),
        "curve": {
            "x": [float(v) for v in x],
            "under_reached_fraction": [float(v) for v in F],
        },
        "falloff": {
            "distance_um": [float(v) for v in d_bins],
            "normalised_mean_auc": [float(v) for v in falloff],
            "count": [int(v) for v in n_bins],
            "bin_um": float(prof["bin_um"]),
            "min_bin_count": int(MIN_BIN_COUNT),
        },
        "E_ref_uM_s": e_ref,
        "E_ref_source": E_REF_SOURCE,
        "n_viable": n_viable,
        "n_tumor": n_tumor,
        "viable_definition": VIABLE_DEFINITION,
        "hypoxia_threshold": float(hypoxia_threshold),
        "convention": CONVENTION_NOTE,
        "scale_invariance": SCALE_INVARIANCE_NOTE,
        "penetration_fit_error": fit_error,
        "units": _UNITS,
    }


def _distance_for(micro: Microenvironment) -> np.ndarray:
    """The PRIMARY distance-from-nearest-vessel map for this tissue.

    Only reached when a caller supplied no ``distance_um``; it is the same
    direct nearest-endothelial-cell distance
    (:func:`vessel_profile.direct_vessel_distance_um`) that
    :func:`delivery_metrics.penetration_length_fit` defaults to, so the two
    profiles are the same profile.
    """
    from tmesim.validation.vessel_profile import (
        direct_vessel_distance_um,
        endothelial_xy,
    )

    return direct_vessel_distance_um(micro.domain, endothelial_xy(micro.cells))
