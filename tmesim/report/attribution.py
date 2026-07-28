"""The computed, additive, falsifiable rank attribution.

Lifted VERBATIM out of ``experiments/patient_report.py`` so that the CLI report
and the Studio's ``POST /api/reports`` share one implementation. The arithmetic
is unchanged; only its address is.

The decomposition
-----------------
Predicted kill is driven by exposure relative to potency, then discounted by the
tumor's own oxygen-driven resistance (hypoxic AND quiescent cells). Each drug's
standing is decomposed into three **exactly additive log10 terms** ("decades"),
each centred on the panel mean so they read as advantage (+) / penalty (-)
versus the panel:

    C_tumor_i      = mean AUC sampled at this tumor's TUMOR-cell positions / T
    drive_i        = C_tumor_i / IC50_i                    (dimensionless)

    delivery_i     = log10 C_tumor_i - mean_j log10 C_tumor_j
    potency_i      = -( log10 IC50_i - mean_j log10 IC50_j )
    o2resist_i     = log10 rho_i     - mean_j log10 rho_j

      where rho_i = kill_i / kill_i(resistance OFF) is measured by a SECOND,
      otherwise-identical simulation of the same drug on the same tissue with
      ``PDConfig(hypoxia_resistance=1, quiescent_resistance=1)`` -- i.e. the
      resistance term is a counterfactual simulation, not a hand-wave. Because
      BOTH penalties are switched off together, the term is named for the whole
      oxygen-driven resistance state rather than for hypoxia alone.

    delivery_i + potency_i  ==  log10 drive_i - mean_j log10 drive_j   (exact)

The "dominant driver" for a drug is the term of largest magnitude, and the
reason for each ADJACENT RANK GAP is the term that differs most between the two
drugs. :func:`attribute` also returns the Spearman rank correlation between the
summed attribution and the actual simulated kill ranking, so the attribution is
itself falsifiable rather than decorative.
"""
from __future__ import annotations

from typing import Dict, List, Sequence

import numpy as np

__all__ = [
    "TERMS",
    "TERM_LABEL",
    "DEFAULT_WINDOW_S",
    "attribute",
    "gap_reason",
    "pct_rank",
    "spearman",
]

#: The three additive terms, in the order they are reported.
TERMS = ("delivery_decades", "potency_decades", "hypoxia_decades")

TERM_LABEL = {
    "delivery_decades": "DELIVERY",
    "potency_decades": "POTENCY",
    # The counterfactual switches off BOTH the hypoxic and the quiescent kill
    # penalty, so the honest name for this term is the whole oxygen-driven
    # resistance state, not hypoxia alone.
    "hypoxia_decades": "O2-RESISTANCE",
}

#: The exposure window every drug in a panel must share for the decomposition
#: to be comparable (6 h, ``experiments/patient_report.py``).
DEFAULT_WINDOW_S = 3600.0 * 6


def pct_rank(values: Sequence[float], i: int,
             higher_is_better: bool = True) -> float:
    """Percentile of item ``i`` within the panel (0 = worst, 100 = best)."""
    v = np.asarray(values, dtype=np.float64)
    n = v.size
    if n < 2:
        return 100.0
    order = np.argsort(v if higher_is_better else -v)
    rank = int(np.where(order == i)[0][0])
    return 100.0 * rank / (n - 1)


def spearman(a: Sequence[float], b: Sequence[float]) -> float:
    from scipy.stats import rankdata

    ra, rb = rankdata(a), rankdata(b)
    if ra.std() == 0 or rb.std() == 0:
        return float("nan")
    return float(np.corrcoef(ra, rb)[0, 1])


def attribute(rows: List[Dict], window_s: float = DEFAULT_WINDOW_S) -> Dict:
    """Attach the exact additive log10 attribution to every row, in place.

    ``delivery + potency`` reproduces the drug's log10 exposure-to-potency drive
    relative to the panel mean *exactly*; ``hypoxia`` is measured by the
    counterfactual resistance-off simulation. Returns the method's own
    faithfulness check.

    Each row must carry ``name``, ``tumor_auc`` (mean AUC in uM*s at the TUMOR
    cell positions), ``ic50`` (uM), ``kill`` and ``kill_no_resistance``.
    """
    eps = 1e-30
    c_tumor = np.array([max(r["tumor_auc"], eps) / window_s for r in rows])
    ic50 = np.array([max(r["ic50"], eps) for r in rows])
    rho = np.array([
        (r["kill"] / r["kill_no_resistance"]) if r["kill_no_resistance"] > 0 else 1.0
        for r in rows])
    rho = np.clip(rho, eps, None)

    log_c, log_ic50, log_rho = np.log10(c_tumor), np.log10(ic50), np.log10(rho)
    delivery = log_c - log_c.mean()
    potency = -(log_ic50 - log_ic50.mean())
    hypoxia = log_rho - log_rho.mean()

    kills = [r["kill"] for r in rows]
    for i, r in enumerate(rows):
        r["c_tumor_uM"] = float(c_tumor[i])
        r["drive"] = float(c_tumor[i] / ic50[i])
        r["rho"] = float(rho[i])
        r["hypoxia_cost"] = float(1.0 - rho[i])
        r["delivery_decades"] = float(delivery[i])
        r["potency_decades"] = float(potency[i])
        r["hypoxia_decades"] = float(hypoxia[i])
        r["total_decades"] = float(delivery[i] + potency[i] + hypoxia[i])
        r["exposure_pct"] = pct_rank(c_tumor, i, higher_is_better=True)
        r["potency_pct"] = pct_rank(ic50, i, higher_is_better=False)
        r["hypoxia_pct"] = pct_rank(rho, i, higher_is_better=True)
        r["penetration_pct"] = pct_rank([x["L_um"] for x in rows], i, True)
        dom = max(TERMS, key=lambda k: abs(r[k]))
        # Name the term AND its direction — "DELIVERY" alone is ambiguous
        # between "wins on delivery" and "is crippled by delivery".
        r["dominant"] = (f"{TERM_LABEL[dom]} "
                         f"({'advantage' if r[dom] >= 0 else 'penalty'})")
        r["dominant_term"] = TERM_LABEL[dom]
        r["dominant_key"] = dom

    totals = [r["total_decades"] for r in rows]
    return {
        "spearman_total_vs_kill": spearman(totals, kills),
        "spearman_delivery_vs_kill": spearman(
            [r["delivery_decades"] for r in rows], kills),
        "spearman_potency_vs_kill": spearman(
            [r["potency_decades"] for r in rows], kills),
        "panel_mean_hypoxia_cost": float(np.mean([r["hypoxia_cost"] for r in rows])),
        "window_s": float(window_s),
    }


def gap_reason(better: Dict, worse: Dict, *, from_below: bool = False) -> str:
    """Which attribution term actually explains an adjacent rank gap.

    Deltas are always computed as ``better - worse`` (the higher-ranked drug
    minus the lower-ranked one); ``from_below`` only changes how the sentence is
    phrased, so the last-placed drug is described relative to the drug ABOVE it
    rather than nonsensically against itself.
    """
    deltas = {k: better[k] - worse[k] for k in TERMS}
    key = max(deltas, key=lambda k: abs(deltas[k]))
    fold = 10.0 ** abs(deltas[key])
    if from_below:
        # Restate every delta from the LOWER-ranked drug's point of view, so the
        # signs in one sentence all share a single reference drug.
        others = ", ".join(f"{TERM_LABEL[k].lower()} {-deltas[k]:+.2f}"
                           for k in TERMS if k != key)
        return (
            f"Versus **{better['name']}** immediately above it, the term that "
            f"differs most is **{TERM_LABEL[key]}** ({-deltas[key]:+.2f} "
            f"decades for this drug ≈ {fold:.1f}× "
            f"{'against it' if deltas[key] > 0 else 'in its favour'}); its "
            f"other terms stand at {others} relative to that drug."
        )
    others = ", ".join(
        f"{TERM_LABEL[k].lower()} {deltas[k]:+.2f}" for k in TERMS if k != key)
    return (
        f"Versus **{worse['name']}** immediately below it, the term that differs "
        f"most is **{TERM_LABEL[key]}** ({deltas[key]:+.2f} decades ≈ "
        f"{fold:.1f}× {'in its favour' if deltas[key] > 0 else 'against it'}); "
        f"the other terms differ by {others}."
    )
