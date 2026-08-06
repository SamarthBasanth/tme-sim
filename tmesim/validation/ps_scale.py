"""The transvascular rate scale ``k_0`` — the missing unit in the delivery law.

THE DEFECT THIS EXISTS TO MAKE VISIBLE
---------------------------------------
``treatment.drugs.delivery_source`` computes::

    S = vascular_density * bbb_penetrance * plasma_conc

and both ``sim/engine.py`` and ``treatment/drugs.py`` declare the result to be
in ``uM/s``. It is not. Dimensionally::

    vascular_density [dimensionless volume fraction]
        x bbb_penetrance [dimensionless PARTITION coefficient]
        x plasma [uM]                                        =  uM,  NOT uM/s

``bbb_penetrance`` is not a permeability and the repo says so in three places
(``core/config.py``'s field docstring, ``reconstruct/barrier.py`` — "It is NOT
a permeability with physical units" — and ``molecular/permeability.py``, whose
own calibration is against measured tissue:plasma RATIOS). The engine has been
obtaining ``s^-1`` from an implicit, undocumented conversion factor of exactly
``1 s^-1`` that nothing anywhere derives.

This module makes that factor explicit, gives it a name, and calibrates it
against independently measured brain permeability-surface-area products.

WHY IT WAS INVISIBLE
--------------------
The GBM transport anchor is a penetration DEPTH: ``L = sqrt(D_eff / k)`` with
``k = clearance + binding``. ``PS`` does not appear in that expression at all,
so the doxorubicin-fit/mitoxantrone-predicted chain could not have caught a
scale error of ANY magnitude. See
``results/transvascular_exchange_revalidation.md`` section 6.2.

WHAT A SINGLE ``k_0`` CAN AND CANNOT DO
----------------------------------------
It can remove the systematic BIAS. It provably cannot remove the drug-to-drug
SPREAD: the two measured anchors below imply ``k_0`` values that disagree by
8.70x, so the best any single constant achieves is to split that disagreement,
leaving a residual of exactly ``sqrt(8.70) = 2.95x``. That residual is not a
fitting failure — it is the engine's cross-drug permeability ordering being
wrong, and it needs a better ``bbb_penetrance``, not a better ``k_0``.

Note the ordering error runs in BOTH directions and so cannot be rescaled
away: in brain the engine UNDER-separates temozolomide:methotrexate (5.57x
against a measured 48.5x), while in breast it OVER-separates
methotrexate:5-fluorouracil (28.1x against a measured 1.02x).

HONESTY ON THE FIT
------------------
``k_0`` is fit ONLY to independently measured brain PS values. It is NOT fit
to the breast microdialysis anchor (``results/breast_kp_validation.md``), to
any kill fraction, or to any outcome. Fitting it to the thing it is then
scored against would be circular, and this project's own standard forbids it.
The breast anchor is therefore still available as an out-of-sample check.
"""
from __future__ import annotations

import math
from typing import Dict

__all__ = [
    "MEASURED_BRAIN_PS",
    "ENGINE_TISSUE_MEAN_PS",
    "calibrate_ps_scale",
    "scorecard",
]

#: Independently MEASURED blood-to-brain permeability-surface-area products,
#: in 1/s (equivalently mL/g/s; brain density ~1.04 g/mL makes the distinction
#: a 4% effect, far below the discrepancies at stake here).
#:
#: These are in vivo rodent measurements with an INTACT barrier. They are the
#: only quantity commensurable with what the engine's PS field claims to be.
MEASURED_BRAIN_PS: Dict[str, Dict] = {
    "methotrexate": {
        "ps_s_inv": 3.3e-5,
        "note": (
            "PA = 3.3e-5 /s intact (rising to 28e-5 after osmotic barrier "
            "opening). Ohno K, Fredericks WR, Rapoport SI, Surg Neurol "
            "1979;12(4):323-8, PMID 524251. Corroborated independently by "
            "Shapiro WR et al., Cancer Res 1988;48:694-701, PMID 3335031 "
            "(K1 = 1.3-1.9 uL/g/min in rat cortex = 2.2-3.2e-5 /s), a "
            "different lab and a different method (quantitative "
            "autoradiography). Note Shapiro also measured methotrexate "
            "in an actual C6 glioma at 1.3-8.2 uL/g/min, so the 'tumour "
            "vasculature is leakier' defence does not rescue the engine's "
            "value either."
        ),
    },
    "temozolomide": {
        "ps_s_inv": 96e-3 / 60.0,   # 96 uL/g/min -> 1.6e-3 mL/g/s
        "note": (
            "PS = 96 +/- 22 uL/g/min, rat in situ brain perfusion, n=24, "
            "3H-labelled. Smith QR et al., Fluids Barriers CNS 2024, PMID "
            "39690422. WEAKER than the methotrexate anchor: a single primary "
            "source, and in situ saline perfusion rather than in vivo, so "
            "plasma protein binding and efflux are not intact. Correcting "
            "for free fraction (methotrexate ~50% bound, temozolomide ~15%) "
            "shrinks the cross-drug disagreement from 8.70x to ~5.1x but "
            "does not remove it."
        ),
    },
}

#: The engine's own tissue-mean PS on the real GBM flagship section
#: (``data_cache/xenium_gbm``, dx = 12 um, oxygen_mode='pde',
#: state_mode='oxygen'), i.e. ``mean(vascular_density * bbb_penetrance)`` over
#: the tissue mask. Measured directly; tissue-mean vascular density is 0.0268.
ENGINE_TISSUE_MEAN_PS: Dict[str, float] = {
    "methotrexate": 6.039e-4,
    "temozolomide": 3.364e-3,
    # Recorded for context, not used in the fit -- no measured brain PS for
    # doxorubicin could be obtained. Smith 2024 explicitly EXCLUDED it for
    # "marked active efflux"; free doxorubicin in normal rat brain is below
    # detection (Wohlfart 2011, PMID 21616104). Qualitatively this value is
    # the worst of the three, but it cannot be scored.
    "doxorubicin": 5.519e-4,
}


def calibrate_ps_scale() -> Dict:
    """Fit the single scale ``k_0`` that best reconciles engine PS with measured.

    Geometric-mean fit over the anchored drugs, which is the right estimator
    for a multiplicative discrepancy and makes the residual symmetric in
    log-space by construction.

    Returns a dict with ``k0``, the per-drug implied scales, the disagreement
    between them, and the residual fold error before and after.
    """
    implied = {
        d: MEASURED_BRAIN_PS[d]["ps_s_inv"] / ENGINE_TISSUE_MEAN_PS[d]
        for d in MEASURED_BRAIN_PS
    }
    logs = [math.log(v) for v in implied.values()]
    k0 = math.exp(sum(logs) / len(logs))

    def gm_fold(scale: float) -> float:
        folds = [
            abs(math.log(ENGINE_TISSUE_MEAN_PS[d] * scale
                         / MEASURED_BRAIN_PS[d]["ps_s_inv"]))
            for d in MEASURED_BRAIN_PS
        ]
        return math.exp(sum(folds) / len(folds))

    lo, hi = min(implied.values()), max(implied.values())
    return {
        "k0": k0,
        "implied_per_drug": implied,
        "disagreement_fold": hi / lo,
        "gm_fold_before": gm_fold(1.0),
        "gm_fold_after": gm_fold(k0),
        # The residual is EXACTLY the half-spread: no single constant can do
        # better, and this identity is the proof that the remaining error is
        # an ordering problem in bbb_penetrance, not a scale problem.
        "irreducible_residual": math.sqrt(hi / lo),
        "n_anchors": len(implied),
    }


def scorecard() -> str:
    """Human-readable before/after, for reports and for the CLI."""
    c = calibrate_ps_scale()
    out = [
        "Transvascular rate scale k_0 -- calibration against measured brain PS",
        "=" * 72,
        f"{'drug':16s} {'engine PS':>12s} {'measured':>12s} "
        f"{'k0=1':>10s} {'k0 fit':>10s}",
    ]
    for d in MEASURED_BRAIN_PS:
        e = ENGINE_TISSUE_MEAN_PS[d]
        m = MEASURED_BRAIN_PS[d]["ps_s_inv"]
        out.append(f"{d:16s} {e:12.4e} {m:12.4e} "
                   f"{e / m:9.2f}x {e * c['k0'] / m:9.2f}x")
    out += [
        "",
        f"  fitted k_0                 = {c['k0']:.5f} 1/s",
        f"  per-drug scales disagree by  {c['disagreement_fold']:.2f}x",
        f"  GM |fold| error   {c['gm_fold_before']:.3f}x -> {c['gm_fold_after']:.3f}x",
        f"  irreducible residual       = {c['irreducible_residual']:.3f}x "
        f"(= sqrt of the disagreement, by construction)",
        "",
        "  k_0 fixes the SCALE. It provably cannot fix the drug-to-drug",
        "  SPREAD -- that needs a better bbb_penetrance, not a better k_0.",
    ]
    return "\n".join(out)


if __name__ == "__main__":
    print(scorecard())
