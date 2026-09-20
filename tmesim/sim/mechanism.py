"""Per-cell mechanism attribution: WHY each tumour cell died, or did not.

A run returns scalars — ``kill_fraction``, coverage, AUC. Those say *how much*
happened; they cannot say *why*. The kill model already knows the answer for
every single cell and throws it away one line after computing it. This module
picks it back up.

What it is (and is not)
-----------------------
This is a **read-out of the kill model's own arithmetic**, not a second model.
:func:`tmesim.cells.pkpd.apply_drug_kill` fills a ``diagnostics`` dict with the
exact arrays it decided on — the same objects, not recomputed copies — and every
label here is a comparison between two of those numbers. Nothing is inferred, no
mechanism the engine does not model is invented, and a cell whose fate does not
follow from those numbers is labelled ``UNATTRIBUTED`` and counted, never
rounded into the nearest plausible bucket.

The kill mapping being attributed (:mod:`tmesim.cells.pkpd`) is::

    E       = AUC / T                                   uM   (local exposure)
    ratio   = E / IC50                                  --   (exposure vs potency)
    f       = ratio^h / (1 + ratio^h)                   --   (Hill fraction of ceiling)
    p_hill  = max_kill_rate * f
    m       = sensitivity * resistance                  --   (the cell's own biology)
    p_kill  = clip(p_hill * m, 0, 1)
    dead    = alive AND (u < p_kill),  u ~ U(0,1)

so exactly three things can keep a cell alive: ``f`` was small (not enough drug
arrived), ``m`` was small (the cell's biology suppressed the kill), or the
Bernoulli draw ``u`` came up above a ``p_kill`` that was not small. The
categories below are those three, with the first split by how far short the
exposure fell.

Categories
----------
``KILLED``
    The engine's own ``killed`` mask. Not re-derived — the mask itself.

``SURVIVED_UNREACHED`` — **DELIVERY-limited**
    ``ratio < unreached_ratio``. The drug never arrived. At the default
    ``unreached_ratio = 0.1`` and a typical Hill coefficient ``h = 1.2``, the
    Hill term is ``f = 0.059``: even a perfectly sensitive, fully viable cell
    could not have reached 6% of the drug's achievable kill probability. What is
    wrong here is transport, not biology and not dose at the tumour cell.

``SURVIVED_UNDERDOSED`` — **DOSE-limited**
    ``unreached_ratio <= ratio < adequate_ratio``. Drug arrived, but below the
    concentration at which it is half-effective. ``adequate_ratio = 1.0`` is not
    a tuned constant: ``ratio == 1`` *is* ``E == IC50``, the pharmacological
    definition of the half-effect point, and the exact point where ``f = 0.5``.

``SURVIVED_RESISTANT`` — **RESISTANCE-limited**
    Exposure was adequate (``ratio >= adequate_ratio``) but
    ``m <= resistant_modulation``. The default ``resistant_modulation`` is
    **derived from the PDConfig the run actually used** — it is
    ``max(hypoxia_resistance, quiescent_resistance)``, i.e. the *mildest*
    penalty the kill model itself applies. A cell is called resistant exactly
    when its own sensitivity and state suppressed its kill probability by at
    least as much as the weakest resistance term the model has. On the
    production ``PDConfig`` that is 0.5.

``SURVIVED_STOCHASTIC`` — **CHANCE**
    Exposure adequate, biology adequate, ``p_kill >= stochastic_min_p_kill`` —
    and the cell still drew above it. Its neighbour with identical inputs died.
    The default floor is likewise derived, not chosen: it is
    ``f(adequate_ratio) * resistant_modulation``, the smallest ``p_kill`` this
    branch can produce when ``max_kill_rate == 1``.

``UNATTRIBUTED``
    A survival the four categories above do not explain. With the default
    (derived) thresholds and ``max_kill_rate == 1`` this set is **empty by
    construction**, which is what makes it a live check rather than decoration:
    it fills up precisely when something outside the bands suppressed the kill —
    a ``max_kill_rate`` ceiling below 1, a clip, or a non-finite exposure. Its
    count is always reported.

``ALREADY_DEAD``
    Tumour cells that were ``DEAD`` before the drug was applied. They are
    counted in ``kill_fraction`` (which is ``1 - alive_tumour/original_tumour``)
    but the drug did not kill them, so they are held separately rather than
    inflating ``KILLED``.

``NOT_TUMOR``
    Non-tumour cells. The payload is aligned to the geometry payload's cell
    order, which is the whole population, so they need a code. ``kill_fraction``
    is a tumour-only quantity and they take no part in the breakdown.

Precedence, and why survivals can be over-determined
----------------------------------------------------
A cell can be short of drug *and* resistant. Only one label fits in a byte, so
the order above is applied as written: **exposure bands first**. The reason is
counterfactual leverage — below ``IC50`` no amount of restored sensitivity can
lift the Hill term above 0.5, so the binding constraint is the exposure. The
cells for which both hold are not hidden: ``joint_limited`` in the aggregate
reports how many ``UNREACHED`` / ``UNDERDOSED`` cells are *also* below the
resistance threshold, i.e. how many would still not have died on dose alone.

Reconciliation
--------------
:func:`reconcile` recomputes ``kill_fraction`` from the labels using the *same*
integer arithmetic ``CellPopulation.viable_tumor_fraction`` uses, and compares
with ``==``. Not ``allclose`` — bit-for-bit. If an attribution cannot reproduce
the number it claims to explain, it is lying about the run, and
:func:`attribute` raises :class:`MechanismReconciliationError` rather than
returning a plausible-looking breakdown.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Sequence

import numpy as np

from ..core import CellType, PDConfig

__all__ = [
    "LABELS",
    "LABEL_CODE",
    "LABEL_LIMITER",
    "KILLED",
    "SURVIVED_UNREACHED",
    "SURVIVED_UNDERDOSED",
    "SURVIVED_RESISTANT",
    "SURVIVED_STOCHASTIC",
    "UNATTRIBUTED",
    "ALREADY_DEAD",
    "NOT_TUMOR",
    "MechanismThresholds",
    "MechanismReconciliationError",
    "classify",
    "reconcile",
    "attribute",
]

# --------------------------------------------------------------------------- #
# Label codes. These are the wire values of the u8 mechanism payload; append
# only, never renumber.
# --------------------------------------------------------------------------- #
KILLED = 0
SURVIVED_UNREACHED = 1
SURVIVED_UNDERDOSED = 2
SURVIVED_RESISTANT = 3
SURVIVED_STOCHASTIC = 4
UNATTRIBUTED = 5
ALREADY_DEAD = 6
NOT_TUMOR = 7

LABELS = (
    "KILLED",
    "SURVIVED_UNREACHED",
    "SURVIVED_UNDERDOSED",
    "SURVIVED_RESISTANT",
    "SURVIVED_STOCHASTIC",
    "UNATTRIBUTED",
    "ALREADY_DEAD",
    "NOT_TUMOR",
)

LABEL_CODE = {name: i for i, name in enumerate(LABELS)}

#: The limiting factor each survival category names, for the aggregate roll-up.
LABEL_LIMITER = {
    "SURVIVED_UNREACHED": "delivery",
    "SURVIVED_UNDERDOSED": "dose",
    "SURVIVED_RESISTANT": "resistance",
    "SURVIVED_STOCHASTIC": "chance",
    "UNATTRIBUTED": "unattributed",
}

#: Short, category-by-category text the UI shows next to the legend.
LABEL_DESCRIPTION = {
    "KILLED": "killed by the drug (the engine's own Bernoulli outcome)",
    "SURVIVED_UNREACHED":
        "DELIVERY-limited: local exposure E/IC50 below the unreached band — the "
        "drug never arrived in a pharmacologically meaningful amount",
    "SURVIVED_UNDERDOSED":
        "DOSE-limited: drug arrived but E stayed below IC50, so the Hill term "
        "never reached half of its ceiling",
    "SURVIVED_RESISTANT":
        "RESISTANCE-limited: exposure was adequate, but sensitivity x state "
        "suppressed the kill probability by at least the model's own mildest "
        "resistance penalty",
    "SURVIVED_STOCHASTIC":
        "CHANCE: exposure and biology both adequate; the cell drew above its "
        "kill probability",
    "UNATTRIBUTED":
        "not explained by the bands above — reported, never reassigned",
    "ALREADY_DEAD":
        "tumour cell already DEAD before the drug was applied; counted in "
        "kill_fraction but not killed by the drug",
    "NOT_TUMOR":
        "non-tumour cell; kill_fraction is a tumour-only quantity",
}


class MechanismReconciliationError(AssertionError):
    """The labels do not reproduce the engine's own ``kill_fraction``."""


@dataclass(frozen=True)
class MechanismThresholds:
    """The bands that turn the kill model's numbers into a label.

    Every value is a parameter, and the two that would otherwise be arbitrary
    default to ``None`` and are **derived from the run's own configuration** by
    :meth:`resolve` rather than picked.

    unreached_ratio : float
        ``E / IC50`` below which the survival is called DELIVERY-limited.
        Default 0.1 — at a typical ``hill = 1.2`` that is 5.9% of the drug's
        achievable kill probability, i.e. the drug is absent for practical
        purposes. This is the one genuinely conventional constant here; it is
        exposed so a caller can move it and see the breakdown move.
    adequate_ratio : float
        ``E / IC50`` at or above which exposure is called adequate. Default 1.0,
        which is not a choice: it is the definition of IC50 and the exact point
        at which the Hill term equals half its ceiling.
    resistant_modulation : float, optional
        ``sensitivity * resistance`` at or below which an adequately-exposed
        survival is called RESISTANCE-limited. ``None`` (default) derives
        ``max(pd.hypoxia_resistance, pd.quiescent_resistance)`` — the mildest
        resistance penalty the kill model itself applies.
    stochastic_min_p_kill : float, optional
        The smallest ``p_kill`` that may be called CHANCE. ``None`` (default)
        derives ``f(adequate_ratio) * resistant_modulation``, the lowest value
        the STOCHASTIC branch can produce at ``max_kill_rate == 1``, so that
        UNATTRIBUTED is empty by construction unless the model applied a
        suppression the bands do not describe.
    """

    unreached_ratio: float = 0.1
    adequate_ratio: float = 1.0
    resistant_modulation: Optional[float] = None
    stochastic_min_p_kill: Optional[float] = None
    #: Where the hypoxia half of `resistant_modulation` came from. A real
    #: field rather than an attribute set after construction, because this is
    #: a frozen dataclass -- and it SHOULD be frozen: thresholds that can be
    #: mutated after resolution are thresholds a later line can move.
    hypoxia_source: Optional[str] = None

    def resolve(self, pd: PDConfig, hill: float,
                drug: Optional["DrugParams"] = None) -> "MechanismThresholds":
        """Return a copy with every ``None`` replaced by its derived value.

        ``drug`` IS LOAD-BEARING AND WAS MISSING UNTIL 2026-09-20.

        The derived resistance threshold has to be the mildest penalty the KILL
        MODEL actually applies, or the classifier draws its boundary somewhere
        the model never goes. This read ``pd.hypoxia_resistance`` directly --
        the global default -- while ``cells/pkpd.py`` applies a PER-DRUG
        override when the drug carries one, and that module's docstring says in
        so many words that any module recomputing this path MUST call
        ``effective_hypoxia_resistance`` instead of reading the default. This
        one did not.

        MEASURED on run r_00099 (temozolomide, flagship): the threshold was
        recorded as 0.50 from ``max(0.35, 0.50)``, while the model applies 0.80
        to that compound, so the correct derivation is ``max(0.80, 0.50)``.

        Latent on the flagship, which has no hypoxic cells at all -- but it
        fires wherever the hypoxic band is occupied, and two independent audits
        found it on the same day from opposite directions.

        ``drug=None`` keeps the old behaviour and SAYS SO in the resolved
        thresholds, rather than being silently indistinguishable from a drug
        that happens to set no override.
        """
        res = self.resistant_modulation
        hyp_source = None
        if res is None:
            if drug is not None:
                from ..cells.pkpd import effective_hypoxia_resistance

                hyp = float(effective_hypoxia_resistance(drug, pd))
                hyp_source = "drug override via cells.pkpd"
            else:
                hyp = float(pd.hypoxia_resistance)
                hyp_source = ("PDConfig default -- no drug was supplied, so a "
                              "per-drug override could not be seen")
            res = max(hyp, float(pd.quiescent_resistance))
        floor = self.stochastic_min_p_kill
        if floor is None:
            r = float(self.adequate_ratio)
            f = (r ** float(hill)) / (1.0 + r ** float(hill)) if r > 0 else 0.0
            floor = f * res
        out = MechanismThresholds(
            unreached_ratio=float(self.unreached_ratio),
            adequate_ratio=float(self.adequate_ratio),
            resistant_modulation=float(res),
            stochastic_min_p_kill=float(floor),
            # Recorded so it reaches the run document: a reader who sees 0.50
            # for a compound the model penalises at 0.80 should be able to
            # tell whether the drug was consulted, without re-deriving it.
            hypoxia_source=hyp_source,
        )
        return out

    def as_dict(self) -> Dict:
        return {
            "unreached_ratio": self.unreached_ratio,
            "adequate_ratio": self.adequate_ratio,
            "resistant_modulation": self.resistant_modulation,
            "stochastic_min_p_kill": self.stochastic_min_p_kill,
            "hypoxia_source": getattr(self, "hypoxia_source", None),
        }


# --------------------------------------------------------------------------- #
# Classification
# --------------------------------------------------------------------------- #
def is_resistance_limited(modulation: np.ndarray,
                          thresholds: MechanismThresholds) -> np.ndarray:
    """Did this cell's own biology actually hold the kill down?

    Two conditions, and the second is not redundant:

    * ``m <= resistant_modulation`` — the suppression is at least as strong as
      the mildest resistance penalty the kill model applies, and
    * ``m < 1.0`` — there *is* a suppression.

    Without the second clause the label degenerates exactly where it matters
    most. The counterfactual arm of the rank attribution runs with
    ``PDConfig(hypoxia_resistance=1, quiescent_resistance=1)``; the derived
    threshold is then ``max(1, 1) = 1.0``, every unresisted cell satisfies
    ``m <= 1.0``, and every adequately-exposed survivor would be filed as
    RESISTANCE-limited in the one run where resistance was switched OFF — while
    CHANCE, which is what those survivals actually are, would be empty. The
    guard says the obvious thing out loud: a cell whose modulation is 1.0 had
    nothing suppressed, so nothing about its biology can be the reason it lived.
    """
    m = np.asarray(modulation, dtype=np.float64)
    return (m <= float(thresholds.resistant_modulation)) & (m < 1.0)


def classify(
    cell_type: Sequence[int],
    diagnostics: Dict,
    thresholds: MechanismThresholds,
) -> np.ndarray:
    """Label every cell with the code that explains its fate.

    Parameters
    ----------
    cell_type : array-like of int
        ``tmesim.core.CellType`` codes, in the population's own order — the same
        order the geometry payload is written in.
    diagnostics : dict
        The dict :func:`tmesim.cells.pkpd.apply_drug_kill` filled in. Must carry
        ``alive``, ``killed``, ``ratio``, ``modulation`` and ``p_kill``.
    thresholds : MechanismThresholds
        Already :meth:`~MechanismThresholds.resolve` d — this function does not
        guess a PDConfig.

    Returns
    -------
    np.ndarray of uint8, length ``n_cells``.
    """
    for key in ("alive", "killed", "ratio", "modulation", "p_kill"):
        if key not in diagnostics:
            raise KeyError(
                f"diagnostics is missing {key!r}: pass diagnostics= to "
                f"apply_drug_kill so the attribution reads the kill model's own "
                f"arrays instead of re-deriving them"
            )
    if thresholds.resistant_modulation is None or \
            thresholds.stochastic_min_p_kill is None:
        raise ValueError(
            "thresholds still carry None; call MechanismThresholds.resolve(pd, "
            "hill) first so the derived values are explicit"
        )

    ctype = np.asarray(cell_type)
    n = ctype.size
    alive = np.asarray(diagnostics["alive"], dtype=bool)
    killed = np.asarray(diagnostics["killed"], dtype=bool)
    ratio = np.asarray(diagnostics["ratio"], dtype=np.float64)
    modulation = np.asarray(diagnostics["modulation"], dtype=np.float64)
    p_kill = np.asarray(diagnostics["p_kill"], dtype=np.float64)
    for name, arr in (("alive", alive), ("killed", killed), ("ratio", ratio),
                      ("modulation", modulation), ("p_kill", p_kill)):
        if arr.size != n:
            raise ValueError(
                f"diagnostics[{name!r}] has {arr.size} entries but there are "
                f"{n} cells — the attribution is not aligned to the population"
            )

    labels = np.full(n, np.uint8(NOT_TUMOR), dtype=np.uint8)
    tumor = ctype == int(CellType.TUMOR)
    if not tumor.any():
        return labels

    # Order matters and is documented in the module docstring: exposure bands
    # bind before biology, because below IC50 no restoration of sensitivity can
    # lift the Hill term past half its ceiling.
    usable = np.isfinite(ratio) & np.isfinite(modulation) & np.isfinite(p_kill)
    survived = tumor & alive & ~killed

    labels[tumor & ~alive] = np.uint8(ALREADY_DEAD)
    labels[tumor & killed] = np.uint8(KILLED)

    labels[survived & ~usable] = np.uint8(UNATTRIBUTED)
    rest = survived & usable
    unreached = rest & (ratio < thresholds.unreached_ratio)
    labels[unreached] = np.uint8(SURVIVED_UNREACHED)
    rest = rest & ~unreached
    underdosed = rest & (ratio < thresholds.adequate_ratio)
    labels[underdosed] = np.uint8(SURVIVED_UNDERDOSED)
    rest = rest & ~underdosed
    resistant = rest & is_resistance_limited(modulation, thresholds)
    labels[resistant] = np.uint8(SURVIVED_RESISTANT)
    rest = rest & ~resistant
    stochastic = rest & (p_kill >= thresholds.stochastic_min_p_kill)
    labels[stochastic] = np.uint8(SURVIVED_STOCHASTIC)
    labels[rest & ~stochastic] = np.uint8(UNATTRIBUTED)

    return labels


# --------------------------------------------------------------------------- #
# Reconciliation
# --------------------------------------------------------------------------- #
def reconcile(labels: np.ndarray, kill_fraction: float) -> Dict:
    """Rebuild ``kill_fraction`` from the labels and compare bit-for-bit.

    The arithmetic mirrors
    :meth:`tmesim.core.CellPopulation.viable_tumor_fraction` exactly — integer
    counts, one true division, one subtraction from ``1.0`` — so an equal
    result is an equality of IEEE-754 doubles, not of rounded summaries.
    """
    labels = np.asarray(labels, dtype=np.uint8)
    n_tumor = int((labels != NOT_TUMOR).sum())
    n_killed = int((labels == KILLED).sum())
    n_already_dead = int((labels == ALREADY_DEAD).sum())
    n_dead = n_killed + n_already_dead

    if n_tumor == 0:
        viable = 0.0
    else:
        viable = (n_tumor - n_dead) / n_tumor
    rebuilt = 1.0 - viable

    engine = float(kill_fraction)
    exact = bool(rebuilt == engine)
    return {
        "n_tumor": n_tumor,
        "n_killed": n_killed,
        "n_already_dead": n_already_dead,
        "kill_fraction_from_labels": float(rebuilt),
        "kill_fraction_from_engine": engine,
        "bit_exact": exact,
        "abs_difference": abs(rebuilt - engine),
        "method": (
            "1.0 - (n_tumor - (KILLED + ALREADY_DEAD)) / n_tumor, the same "
            "integer arithmetic CellPopulation.viable_tumor_fraction uses, "
            "compared with == (not allclose)"
        ),
    }


# --------------------------------------------------------------------------- #
# The public entry point
# --------------------------------------------------------------------------- #
def attribute(
    cell_type: Sequence[int],
    diagnostics: Dict,
    kill_fraction: float,
    pd: Optional[PDConfig] = None,
    thresholds: Optional[MechanismThresholds] = None,
    drug: Optional["DrugParams"] = None,
) -> Dict:
    """Classify every cell, aggregate the survivals, and verify the total.

    Returns a dict with ``labels`` (uint8, aligned to the population order),
    ``counts``, ``survival_breakdown``, ``joint_limited``, ``thresholds`` and
    ``reconciliation``.

    Raises
    ------
    MechanismReconciliationError
        If the labels do not reproduce ``kill_fraction`` bit-for-bit, or if any
        tumour cell was left without a label.
    """
    pd = pd if pd is not None else PDConfig()
    hill = float(diagnostics.get("hill", 1.0))
    resolved = (thresholds or MechanismThresholds()).resolve(
        pd, hill, drug=drug)

    labels = classify(cell_type, diagnostics, resolved)
    rec = reconcile(labels, kill_fraction)
    if not rec["bit_exact"]:
        raise MechanismReconciliationError(
            f"mechanism attribution does not reconcile with the engine: labels "
            f"give kill_fraction={rec['kill_fraction_from_labels']!r} but the "
            f"engine reported {rec['kill_fraction_from_engine']!r} "
            f"(|delta| = {rec['abs_difference']:.3e}). The attribution would be "
            f"describing a run that did not happen."
        )

    counts = {name: int((labels == code).sum())
              for code, name in enumerate(LABELS)}
    n_tumor = rec["n_tumor"]
    n_survived = sum(counts[k] for k in LABEL_LIMITER)

    survival_breakdown = {}
    for name, limiter in LABEL_LIMITER.items():
        survival_breakdown[limiter] = (
            counts[name] / n_survived if n_survived else 0.0)

    # Over-determined survivals: short of drug AND below the resistance
    # threshold. They are labelled by the exposure band (documented precedence);
    # this is the honest count of how many of those would ALSO have survived a
    # dose increase.
    ratio = np.asarray(diagnostics["ratio"], dtype=np.float64)
    modulation = np.asarray(diagnostics["modulation"], dtype=np.float64)
    exposure_limited = (labels == SURVIVED_UNREACHED) | (labels == SURVIVED_UNDERDOSED)
    also_resistant = exposure_limited & is_resistance_limited(modulation, resolved)

    tumor_ratio = ratio[np.asarray(cell_type) == int(CellType.TUMOR)]
    finite_ratio = tumor_ratio[np.isfinite(tumor_ratio)]
    median_ratio = float(np.median(finite_ratio)) if finite_ratio.size else None
    shortfall = (
        float(resolved.adequate_ratio) / median_ratio
        if median_ratio and median_ratio > 0 else None
    )

    return {
        "labels": labels,
        "label_names": list(LABELS),
        "label_descriptions": dict(LABEL_DESCRIPTION),
        "counts": counts,
        "n_cells": int(labels.size),
        "n_tumor": n_tumor,
        "n_survived": int(n_survived),
        "killed_fraction_of_tumor": (rec["n_killed"] / n_tumor) if n_tumor else 0.0,
        "survival_breakdown": survival_breakdown,
        "survival_breakdown_note": (
            "fractions of SURVIVING tumour cells, so they sum to 1. "
            "'delivery' = SURVIVED_UNREACHED, 'dose' = SURVIVED_UNDERDOSED, "
            "'resistance' = SURVIVED_RESISTANT, 'chance' = SURVIVED_STOCHASTIC."
        ),
        "joint_limited": {
            "count": int(also_resistant.sum()),
            "fraction_of_exposure_limited": (
                float(also_resistant.sum() / exposure_limited.sum())
                if exposure_limited.any() else 0.0),
            "note": (
                "exposure-limited survivals that are ALSO at or below the "
                "resistance threshold. They are labelled by their exposure band "
                "(documented precedence: below IC50 the exposure is the binding "
                "constraint), but raising the dose alone would not have killed "
                "them. Reported so the delivery/dose share is not read as a "
                "clean counterfactual."
            ),
        },
        "exposure_ratio_percentiles_tumor": {
            "p05": float(np.percentile(finite_ratio, 5)) if finite_ratio.size else None,
            "p50": float(np.percentile(finite_ratio, 50)) if finite_ratio.size else None,
            "p95": float(np.percentile(finite_ratio, 95)) if finite_ratio.size else None,
            "max": float(finite_ratio.max()) if finite_ratio.size else None,
        },
        "exposure_shortfall_fold": shortfall,
        "exposure_shortfall_note": (
            "adequate_ratio / median(E/IC50) over tumour cells: the fold change "
            "in exposure the median tumour cell would need to reach IC50. Above "
            "1 the median cell is short of drug; below 1 it is already past the "
            "half-effect point and something else is holding the kill down."
        ),
        # This is the honesty note the split above requires, and it is a
        # property of the KILL MODEL, not of this module: p_hill depends on E
        # and IC50 only through their ratio, so within one run the model has no
        # separate potency term to attribute to.
        "exposure_axis_note": (
            "UNREACHED and UNDERDOSED are two bands on ONE axis, E/IC50, not two "
            "independent mechanisms. The kill model reaches p_hill through that "
            "ratio alone, so it cannot separate 'too little drug arrived' from "
            "'the drug is too weak' — and neither can this attribution. What IS "
            "separable within a run: IC50 is a single per-drug constant, so the "
            "entire CELL-TO-CELL spread in E/IC50 is delivery, while the "
            "population's overall position relative to 1.0 is set jointly by "
            "delivery and potency. Read the split as a statement about how far "
            "short the exposure fell, not about which upstream term caused it."
        ),
        "thresholds": resolved.as_dict(),
        "threshold_derivation": {
            "unreached_ratio": (
                "convention, exposed as a parameter: E/IC50 below which the "
                "Hill term is a few percent of its ceiling"),
            "adequate_ratio": (
                "definition of IC50 — E == IC50 is exactly where the Hill term "
                "equals half its ceiling"),
            "resistant_modulation": (
                "DERIVED from this run's PDConfig as max(hypoxia_resistance="
                f"{float(pd.hypoxia_resistance)}, quiescent_resistance="
                f"{float(pd.quiescent_resistance)}) — the mildest resistance "
                "penalty the kill model itself applies. A cell is called "
                "RESISTANCE-limited only if it is at or below this AND its "
                "modulation is below 1.0, i.e. something was actually "
                "suppressed"
                + (". This run switched resistance OFF, so the threshold is 1.0 "
                   "and the second condition is what carries the label: no cell "
                   "is resistance-limited here, by construction and correctly."
                   if float(resolved.resistant_modulation) >= 1.0 else "")),
            "stochastic_min_p_kill": (
                "DERIVED as f(adequate_ratio) * resistant_modulation at "
                f"hill={hill}, the lowest p_kill the CHANCE branch can produce "
                "at max_kill_rate = 1, so UNATTRIBUTED is empty by construction "
                "unless the model suppressed the kill in a way the bands do not "
                "describe"),
        },
        "pd": {
            "hypoxia_resistance": float(pd.hypoxia_resistance),
            "quiescent_resistance": float(pd.quiescent_resistance),
            "max_kill_rate": float(pd.max_kill_rate),
        },
        "drug": {
            "ic50_uM": float(diagnostics.get("E50", float("nan"))),
            "hill": hill,
            "exposure_window_s": float(diagnostics.get("duration_s", float("nan"))),
        },
        "reconciliation": rec,
    }
