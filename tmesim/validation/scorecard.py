"""THE EXTERNAL-ANCHOR SCORECARD — one runnable inventory of every place this
model is scored against an INDEPENDENT MEASURED number.

This module is the data layer behind ``vivum validate`` (``tmesim.cli``). It
does not implement a single new comparison. Every number it prints is produced
by calling the shipping validation module that already owns that comparison,
unmodified, and every verdict is that module's own declared bar. What this
module adds is a *registry*: the anchors were previously spread across
``tmesim/validation/*.py`` and ``experiments/*.py`` with no single entry point,
so "how does the model currently score against reality" could not be answered
in one command.

What counts as an anchor here
-----------------------------
A row appears below **only** if it compares a model-produced number to a
quantitative measurement published by someone else, carrying a PMID (or a PMC
id where no PMID exists). Explicitly EXCLUDED, and listed in
:data:`OUT_OF_SCOPE` with reasons:

* internal-consistency, grid-convergence and mass-conservation checks
  (``numerical_rigor``, ``dimensionality``'s analytic identities),
* ablations that re-score an anchor already listed here under a modified model
  (``efflux_impact``, ``ion_trapping_impact``, ``uptake_impact``,
  ``immune_sink_impact.glioma_scorecard`` — all of them re-run
  :func:`glioma_regional_pk.run_scorecard`, so counting them would count the
  same six measurements up to four times),
* real-vs-real corroboration with no model in the loop
  (``pharmacology_real``'s GDSC/HGCC concordance),
* dataset landscape surveys (``results/landscape_*.md``,
  ``results/breast_transport_anchors.md``, ``results/oxygen_calibration_landscape.md``).

The four statuses, and why NOT TESTABLE is load-bearing
-------------------------------------------------------
``PASS``
    The model ran and its prediction landed inside the band the owning module
    declares. The bar is named in every row (``bar`` field) because different
    modules declare different bars, and a scorecard that hides which bar it
    used is worthless.
``FAIL``
    The model ran and missed. **Failures are never suppressed and never
    softened.** ``vivum validate`` exits non-zero when any row is FAIL, so this
    is usable as a CI gate.
``NOT TESTABLE``
    No verdict is possible, for a stated ``reason``. This is not a soft pass and
    must never be pooled with one. Four distinct reasons occur, and conflating
    any of them with PASS is exactly the error this project retracts things for:

    * *in-sample* — the number is the FIT, not a test of it. The doxorubicin
      penetration anchor is solved for in closed form, so the model reproduces
      it by construction and could not have failed.
    * *no measured target exists* — the literature supports a DIRECTION but
      never measured the MAGNITUDE, so there is nothing to score against
      (the immune/TAM sink).
    * *required data absent* — the anchor needs bytes under ``data_cache/``
      that this checkout does not have.
    * *the model cannot express the measured quantity* — see below.
``DEFERRED``
    Skipped by ``--quick``, which runs only the anchors that need no real
    tissue. A DEFERRED row is **not** a pass; anchors that FAIL in the full run
    appear as DEFERRED under ``--quick``, which is why the command prints a
    banner saying so.

ANALYTIC anchors do not test the PDE
------------------------------------
Each row carries a ``scope``. ``PDE`` means the transport solver actually ran.
``analytic`` means the prediction is a closed-form expression over molecular
descriptors and the solver was never invoked — verified for
``glioma_regional_pk`` and ``gist_imatinib`` by call-counter instrumentation in
``results/transvascular_exchange_revalidation.md`` §6.1, which states the rule
directly: *"These modules must be reported as UN-TESTABLE, not as unaffected.
'Robust under exchange' and 'never evaluates the delivery law' are very
different statements and must not be pooled in a coverage count."* An analytic
row still carries a real PASS/FAIL against its measurement — it just cannot be
cited as evidence about the solver, and the printed scorecard says so.

Conditionals
------------
Several anchors are conditional on something uncalibrated or on something the
model was itself fitted to. Rather than bury that, each row lists the keys of
the caveats it inherits (:data:`CONDITIONALS`), and the renderer prints every
cited caveat in full. The two that reach furthest:

* ``alpha-fit`` — every penetration and exposure number downstream of the one
  fitted constant ``alpha`` (``validation.calibration.calibrate_alpha``).
* ``ps-scale`` — every transvascular-exchange number is conditional on
  ``tmesim/treatment/drugs.py::PS_SCALE_S_INV = 1.0``, an implicit, underived
  dimensionless-to-rate conversion (``results/transvascular_exchange_revalidation.md``
  §6.2).

Nothing under ``tmesim/`` outside this file, and nothing under ``results/``, is
written or modified by running this module.
"""
from __future__ import annotations

import os
import sys
import time
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple

# --------------------------------------------------------------------------- #
# Statuses and scopes
# --------------------------------------------------------------------------- #
PASS = "PASS"
FAIL = "FAIL"
NOT_TESTABLE = "NOT TESTABLE"
DEFERRED = "DEFERRED"

#: Ordered worst-first, for sorting and for the summary line.
STATUS_ORDER = (FAIL, NOT_TESTABLE, DEFERRED, PASS)

#: ``scope``: did the transport PDE actually run to produce this prediction?
SCOPE_PDE = "PDE"
SCOPE_ANALYTIC = "analytic"
SCOPE_GEOMETRY = "geometry"

#: Tiers. ``quick`` anchors need no real tissue and finish in seconds.
TIER_QUICK = "quick"
TIER_FULL = "full"

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))


def _repo_path(*parts: str) -> str:
    return os.path.join(REPO_ROOT, *parts)


class NotTestable(Exception):
    """Raised by a runner when no verdict is possible.

    The message becomes the row's ``reason``. Raising this is the ONLY
    sanctioned way for a runner to decline to produce a verdict — returning an
    empty row list, or a PASS, would silently convert "unreachable" into
    "passing".
    """


# --------------------------------------------------------------------------- #
# Named caveats, printed in full whenever any row cites them
# --------------------------------------------------------------------------- #
CONDITIONALS: Dict[str, str] = {
    "alpha-fit": (
        "Inherits alpha, the single free constant fit ONCE in closed form to the "
        "doxorubicin perivascular anchor (PMID 16361566) by "
        "validation.calibration.calibrate_alpha. Anything downstream of alpha is "
        "conditional on that one fit and is not independent of it."
    ),
    "ps-scale": (
        "Conditional on an UNCALIBRATED unit scale. tmesim/treatment/drugs.py::"
        "PS_SCALE_S_INV = 1.0 is an implicit, undocumented, underived conversion "
        "from a dimensionless product to a rate in 1/s; nothing in the repo "
        "derives it. See results/transvascular_exchange_revalidation.md sec. 6.2. "
        "Every transvascular-exchange-ON number carries this."
    ),
    "leak-fit": (
        "The additive paracellular leak constant "
        "reconstruct.barrier.L_PARACELLULAR_LEAK = 0.238 was FITTED to the same "
        "contrast-enhancing-vs-non-enhancing contrasts these rows are scored "
        "against, so the leak-on score is not fully out-of-sample "
        "(results/glioma_regional_pk.md)."
    ),
    "analytic-only": (
        "Computed in closed form from molecular descriptors; the transport PDE is "
        "never invoked. Verified by call-counter instrumentation in "
        "results/transvascular_exchange_revalidation.md sec. 6.1. This row cannot "
        "confirm or refute the solver."
    ),
    "gdsc-potency": (
        "Intrinsic potency comes from the GDSC cell-line screen, not from this "
        "patient's tissue; the 3-D:2-D ratio therefore constrains the PRODUCT of "
        "delivery and cell-state resistance, not either factor alone "
        "(results/kill_calibration.md)."
    ),
    "band-wider-than-assay": (
        "The literature band is wider than the assay spread it is meant to "
        "correct (the 2-D reference's own IQR is ~19x, PMID 34794398), so a pass "
        "here is a weak pass by construction (results/kill_calibration.md)."
    ),
    "untreated-section": (
        "The open spatial section is an UNTREATED resection. There is no paired "
        "measured drug response on this tissue; the anchor is a cohort-level "
        "literature value, not this patient's outcome."
    ),
    "brain-permeability-in-breast": (
        "PS is bbb_penetrance, a blood-BRAIN-barrier term, applied to breast "
        "tissue. results/transvascular_exchange_revalidation.md sec. 5 states the "
        "fix corrected the FORM of the delivery law, not the VALUE of the "
        "permeability."
    ),
}


# --------------------------------------------------------------------------- #
# Row / group records
# --------------------------------------------------------------------------- #
@dataclass
class Row:
    """One model-vs-measurement comparison."""

    group: str
    label: str
    status: str
    scope: str
    quantity: str          # what model quantity vs what measured quantity
    predicted: str
    measured: str
    pmids: Tuple[str, ...] = ()
    bar: str = ""          # the declared pass bar, named
    note: str = ""
    reason: str = ""       # REQUIRED when status is NOT TESTABLE / DEFERRED
    conditional: Tuple[str, ...] = ()
    values: Dict = field(default_factory=dict)   # machine-readable numbers

    def to_dict(self) -> Dict:
        return {
            "group": self.group,
            "label": self.label,
            "status": self.status,
            "scope": self.scope,
            "quantity": self.quantity,
            "predicted": self.predicted,
            "measured": self.measured,
            "pmids": list(self.pmids),
            "bar": self.bar,
            "note": self.note,
            "reason": self.reason,
            "conditional": list(self.conditional),
            "values": dict(self.values),
        }


@dataclass(frozen=True)
class AnchorGroup:
    """One registered anchor: metadata plus the callable that scores it."""

    key: str
    title: str
    tissue: str
    tier: str
    scope: str
    source: str                     # the module/script that owns the comparison
    report: str                     # the results/*.md that publishes it
    runner: Callable[[], List[Row]]
    needs: Tuple[str, ...] = ()     # repo-relative paths that must exist
    est: str = ""                   # rough runtime, for --list


def classify(predicted: float, lo: float, hi: float) -> str:
    """PASS if ``lo <= predicted <= hi``, else FAIL. Bounds are inclusive.

    Deliberately the only place a band comparison is spelled out, so the
    PASS/FAIL rule is testable in isolation without any tissue or PDE.
    """
    p, a, b = float(predicted), float(lo), float(hi)
    if a > b:
        a, b = b, a
    return PASS if (a <= p <= b) else FAIL


# --------------------------------------------------------------------------- #
# QUICK anchors — analytic / literature-table only, no real tissue, no PDE
# --------------------------------------------------------------------------- #
def _run_perivascular_penetration() -> List[Row]:
    """Perivascular penetration length vs measured half-distances.

    Drives ``validation.pk_multidrug.run_penetration`` (which itself drives
    ``calibration.calibrate_alpha``). Doxorubicin is the FIT and is reported as
    NOT TESTABLE; mitoxantrone and topotecan are out-of-sample predictions with
    alpha frozen.
    """
    from .pk_multidrug import run_penetration

    res = run_penetration()
    rows: List[Row] = []
    for r in res["rows"]:
        lo, hi = r["lequiv_band_um"]
        hlo, hhi = r["half_band_um"]
        measured = (f"{hlo:.0f}-{hhi:.0f} um half-distance "
                    f"(= {lo:.0f}-{hi:.0f} um e-fold, d_half/ln2)")
        predicted = f"L = {r['predicted_L_um']:.1f} um"
        values = {
            "predicted_L_um": float(r["predicted_L_um"]),
            "efold_band_um": [float(lo), float(hi)],
            "half_distance_band_um": [float(hlo), float(hhi)],
            "ratio_to_band_mid": float(r["ratio_to_lequiv_mid"]),
            "module_verdict": r["verdict"],
            "alpha": float(res["alpha"]),
        }
        if r["in_sample"]:
            rows.append(Row(
                group="gbm-penetration", label=r["drug"], status=NOT_TESTABLE,
                scope=SCOPE_ANALYTIC,
                quantity="analytic penetration length L = sqrt(D_eff/k_eff) vs "
                         "measured perivascular half-distance (immunofluorescence "
                         "vs distance-from-vessel)",
                predicted=predicted, measured=measured, pmids=("16361566",),
                bar="n/a — this IS the fit",
                reason="IN-SAMPLE. alpha is solved in closed form so that L equals "
                       "this target exactly; the model reproduces it by "
                       "construction and could not have failed. Reporting it as a "
                       "PASS would be circular.",
                conditional=("analytic-only",), values=values,
            ))
            continue
        # Out-of-sample. Use the owning module's own declared verdict rather than
        # inventing a stricter bar here: run_penetration passes on the published
        # band, or on the 3x order-of-magnitude band pk_reality declares.
        verdict = r["verdict"]
        status = PASS if verdict in ("PASS", "PASS(3x)") else FAIL
        note = r.get("class", "")
        if verdict == "PASS(3x)":
            note = ("OUTSIDE the published band; passes only on the 3x "
                    "order-of-magnitude bar (pk_reality.ORDER_OF_MAGNITUDE_FACTOR). "
                    + note)
        rows.append(Row(
            group="gbm-penetration", label=r["drug"], status=status,
            scope=SCOPE_ANALYTIC,
            quantity="OUT-OF-SAMPLE analytic penetration length L = "
                     "sqrt(D_eff/k_eff), alpha frozen, vs measured perivascular "
                     "half-distance",
            predicted=predicted + f"  ({r['ratio_to_lequiv_mid']:.2f}x band mid)",
            measured=measured, pmids=(str(r["pmid"]),),
            bar="inside the published e-fold band, else within 3x of its midpoint",
            note=note, conditional=("alpha-fit", "analytic-only"), values=values,
        ))
    return rows


def _run_bbb_ordering() -> List[Row]:
    """Structure-only ``bbb_penetrance`` ORDERING vs measured CSF/tissue:plasma
    ordering. Ranks only — a 0-1 score is not a partition ratio, and the owning
    module refuses to compare their magnitudes."""
    from .pk_multidrug import run_partition

    res = run_partition()
    real_order = " < ".join(res["real_order"])
    sim_order = " < ".join(res["sim_order"])
    disc = res["discordant_pairs"]
    status = PASS if disc == 0 else FAIL
    note = ""
    if disc:
        pairs = "; ".join(f"{a} vs {b}" for a, b in res["discordant"])
        note = f"discordant pair(s): {pairs}"
    pmids = tuple(str(r["pmid"]) for r in res["rows"])
    return [Row(
        group="bbb-ordering", label="MTX/cisplatin/TMZ", status=status,
        scope=SCOPE_ANALYTIC,
        quantity="rank order of structure-only bbb_penetrance vs rank order of "
                 "measured brain-or-CSF:plasma partition ratio (RANKS ONLY — the "
                 "dimensionless score is not a partition ratio)",
        predicted=f"model order {sim_order} (Kendall tau {res['kendall_tau']:.2f}, "
                  f"Spearman rho {res['spearman_rho']:.2f})",
        measured=f"measured order {real_order}",
        pmids=pmids, bar="every pairwise order concordant",
        note=note, conditional=("analytic-only",),
        values={"kendall_tau": float(res["kendall_tau"]),
                "spearman_rho": float(res["spearman_rho"]),
                "concordant_pairs": int(res["concordant_pairs"]),
                "discordant_pairs": int(disc),
                "real_order": list(res["real_order"]),
                "model_order": list(res["sim_order"])},
    )]


def _run_glioma_partition() -> List[Row]:
    """Six-drug human glioma tissue:plasma partition scorecard.

    ``apply_leak=True`` is the model AS IT SHIPS (the additive paracellular
    leak is on). Scored on the module's own TIGHT band, ``BAND_TIGHT = 3.0``,
    which is the column ``results/glioma_regional_pk.md`` prints PASS/FAIL in;
    the wider 10x order-of-magnitude outcome is carried in the note.
    """
    from .glioma_regional_pk import BAND_OOM, BAND_TIGHT, run_scorecard

    sc = run_scorecard(aggregation="geometric", apply_leak=True)
    rows: List[Row] = []
    for r in sc["rows"]:
        status = PASS if r["pass_3x"] else FAIL
        oom = "inside 10x" if r["pass_10x"] else "OUTSIDE 10x too"
        free = "microdialysis (free drug)" if r["free_drug_comparator"] \
            else "tissue homogenate (TOTAL drug, incl. intracellular)"
        note = f"{oom}; n={r['n_measured']}; comparator = {free}"
        if not r["free_drug_comparator"]:
            note += (". Homogenate measures total drug and can exceed 1 by "
                     "definition, while the engine computes free interstitial — "
                     "the comparison is not modality-matched")
        rows.append(Row(
            group="glioma-partition", label=r["drug"], status=status,
            scope=SCOPE_ANALYTIC,
            quantity="predicted tissue:plasma partition (assumption A1: "
                     "R = clip(bbb_penetrance + L*disruption, 0, 1)) vs MEASURED "
                     "human glioma tissue:plasma AUC ratio",
            predicted=f"{r['predicted_partition']:.4f}  "
                      f"({r['fold_error']:.2f}x measured)",
            measured=f"{r['measured_mean']:.4f} "
                     f"(range {r['measured_min']:.4f}-{r['measured_max']:.4f})",
            pmids=tuple(str(r["pmids"]).split("|")),
            bar=f"fold error within {BAND_TIGHT:g}x (module BAND_TIGHT); "
                f"{BAND_OOM:g}x reported alongside",
            note=note, conditional=("analytic-only", "leak-fit"),
            values={"predicted_partition": float(r["predicted_partition"]),
                    "measured_mean": float(r["measured_mean"]),
                    "fold_error": float(r["fold_error"]),
                    "pass_3x": bool(r["pass_3x"]), "pass_10x": bool(r["pass_10x"]),
                    "n_measured": int(r["n_measured"]), "method": str(r["method"])},
        ))
    return rows


def _run_tmz_exvivo() -> List[Row]:
    """GDSC-extrapolated temozolomide viability at 100 uM vs MEASURED ex-vivo
    ATP viability in patient-derived glioma stem-cell cultures (n=40)."""
    from .gbm_tmz_exvivo import PMID_NTAFOULIS, test_a_potency_axis, verify_provenance

    prov = verify_provenance()
    a = test_a_potency_axis(prov)
    # The module reports Mann-Whitney as a ``(U, p)`` tuple. Prefer the n=40
    # cohort (parsed from the primary .docx supplement); fall back to the
    # audited n=17 derived CSV, which is what a checkout without python-docx
    # gets.
    mw = a.get("mw_pred_vs_meas40") or a.get("mw_pred_vs_meas17")
    if mw is None:
        raise NotTestable(
            "gbm_tmz_exvivo.test_a_potency_axis returned no Mann-Whitney "
            "comparison against the measured cohort")
    p = float(mw[1]) if isinstance(mw, (tuple, list)) else float(mw)
    pred, meas = a["pred_median"], a["meas40_median"]
    fold = pred / meas if meas else float("inf")
    # The declared bar is a distributional one: the paper's measured viability
    # distribution and the model's extrapolated one must not be separable.
    status = PASS if float(p) >= 0.05 else FAIL
    return [Row(
        group="gbm-tmz-exvivo", label="temozolomide", status=status,
        scope=SCOPE_ANALYTIC,
        quantity="GDSC Hill curve extrapolated to the assay's 100 uM dose "
                 "(3.3x above the GDSC concentration ceiling) vs MEASURED ATP "
                 "viability in patient-derived GSC cultures",
        predicted=f"{pred:.1f}% viability (median, IQR "
                  f"{a['pred_iqr'][0]:.1f}-{a['pred_iqr'][1]:.1f})",
        measured=f"{meas:.1f}% viability (median, n={len(a['meas40'])})",
        pmids=(str(PMID_NTAFOULIS),),
        bar="Mann-Whitney p >= 0.05 (distributions not separable)",
        note=f"Mann-Whitney p = {float(p):.3f}; predicted/measured = {fold:.2f}x. "
             f"Genuinely out-of-range: GDSC's own dose window tops out below the "
             f"assay concentration.",
        conditional=("gdsc-potency",),
        values={"pred_median_pct": float(pred), "meas_median_pct": float(meas),
                "fold": float(fold), "mannwhitney_p": float(p),
                "n_measured": int(len(a["meas40"]))},
    )]


def _run_gist_imatinib() -> List[Row]:
    """Analytic Krogh-slab tumour:normal imatinib partition vs MEASURED
    MALDI-MSI tumour:normal ion ratios in human GIST (n=6 matched patients).

    Published as FALSIFIED: the transport-only model barely beats a
    no-partitioning null. This row exists so the falsification stays visible.
    """
    import tmesim.validation.gist_imatinib as G

    df = G.load_join()
    char = G.characterize(df)
    targets = G.measured_targets(char)
    null = G.plasma_null(char)
    compiled = G.compile_imatinib()
    pred = G.predict(compiled)
    sc = G.score(targets, compiled, pred, null)
    model_gm = float(sc["M1"]["geomean_fold_error"])
    null_gm = float(sc["null_unity"]["geomean_fold_error"])
    rho = float(sc["M1"]["spearman_rho"])
    # Declared bar (results/gist_imatinib.md): a transport-only explanation has
    # to beat the no-partitioning null by more than a factor of 3 to count as
    # sufficient. It beats it by ~1.6x, which is why the report calls the model
    # FALSIFIED as a sufficient explanation rather than merely imprecise.
    beats_by = null_gm / model_gm if model_gm > 0 else float("nan")
    status = PASS if beats_by >= 3.0 else FAIL
    return [Row(
        group="gist-imatinib", label="imatinib", status=status,
        scope=SCOPE_ANALYTIC,
        quantity="analytic 1-D Krogh-slab delivery efficiency x tissue "
                 "partitioning -> tumour:normal imatinib ratio, vs MEASURED "
                 "MALDI-MSI tumour:normal ion ratio in human GIST",
        predicted=f"GM |fold error| {model_gm:.2f}x over n={int(sc['n'])} patients "
                  f"(Spearman rho {rho:.3f})",
        measured=f"no-partitioning null gives {null_gm:.2f}x "
                 f"(model beats null by only {beats_by:.2f}x)",
        pmids=("PMC6650609",),
        bar="beat the no-partitioning null by >= 3x GM fold error",
        note="FALSIFIED as a sufficient explanation (results/gist_imatinib.md): "
             "4 of 5 liver ratios lie BELOW the model's structural floor, so no "
             "penetration length can reproduce them. n=6 with 5 of 6 sharing the "
             "same tissue pairing, i.e. ~one independent comparison.",
        conditional=("analytic-only",),
        values={"model_gm_abs_fold": model_gm, "null_gm_abs_fold": null_gm,
                "beats_null_by": beats_by, "n": int(sc["n"]),
                "spearman_rho": rho},
    )]


def _run_ps_scale() -> List[Row]:
    """The implicit unit scale ``PS_SCALE_S_INV`` against measured intact-BBB
    permeability-surface-area products.

    Reported NOT TESTABLE, and the reason is the whole point: the engine's
    ``bbb_penetrance`` is a dimensionless PARTITION, while the delivery law uses
    it as a RATE in 1/s. The two readings are not the same quantity, so no
    single scale can be scored — and the drug-to-drug incoherence proves it.
    """
    from .ps_scale import calibrate_ps_scale

    d = calibrate_ps_scale()
    per = ", ".join(f"{k} {v:.3g}" for k, v in d["implied_per_drug"].items())
    return [Row(
        group="ps-scale", label="PS_SCALE_S_INV", status=NOT_TESTABLE,
        scope=SCOPE_ANALYTIC,
        quantity="engine volume-averaged permeability-surface-area product "
                 "(1/s) vs MEASURED intact-blood-brain-barrier PA product",
        predicted=f"a single k0 = {d['k0']:.4g} 1/s would be needed; the "
                  f"per-drug values disagree by {d['disagreement_fold']:.2f}x "
                  f"({per})",
        measured=f"n={int(d['n_anchors'])} measured brain PA anchors",
        pmids=("524251", "39690422"),
        bar="n/a",
        reason="THE MODEL CANNOT EXPRESS THE MEASURED QUANTITY. bbb_penetrance is "
               "a dimensionless brain:plasma PARTITION, but the delivery law "
               "consumes it as a first-order RATE via an implicit, underived "
               "conversion of exactly 1 1/s (PS_SCALE_S_INV). Read as rates the "
               "drugs are mutually incoherent by "
               f"{d['disagreement_fold']:.1f}x; read as partitions they agree to "
               "~1.04x. That is the signature of a units error, not a "
               "calibration error, and no single k0 can remove it (irreducible "
               f"residual {d['irreducible_residual']:.2f}x). See "
               "results/transvascular_exchange_revalidation.md sec. 6.2.",
        conditional=("ps-scale", "analytic-only"),
        values={k: (float(v) if isinstance(v, (int, float)) else v)
                for k, v in d.items() if k != "implied_per_drug"},
    )]


def _run_immune_sink() -> List[Row]:
    """The immune (TAM) drug sink.

    NOT TESTABLE, and this is a literature fact rather than a missing file: no
    drug in the panel has a MEASURED free-drug TAM uptake ratio. The cited work
    supports a DIRECTION only, and the single best-quantified case is a depot,
    i.e. the OPPOSITE sign. The sink strength is therefore swept and never
    fitted, and it ships default-off.
    """
    return [Row(
        group="immune-sink", label="TAM sequestration", status=NOT_TESTABLE,
        scope=SCOPE_PDE,
        quantity="first-order drug loss proportional to the MEASURED local "
                 "immune-cell area fraction vs a measured free-drug "
                 "macrophage:interstitium uptake ratio",
        predicted="swept over the strength range, never fitted",
        measured="NONE — 0 of 3 panel drugs has a measured free-drug TAM uptake "
                 "ratio",
        pmids=("26503691",),
        bar="n/a",
        reason="NO MEASURED TARGET EXISTS. The literature supports a DIRECTION "
               "and not a MAGNITUDE, and for free doxorubicin the one direct "
               "measurement is NEGATIVE (a depot, the opposite sign). The "
               "module's permutation test (structure ratio 2.0025x, p = 0.0476) "
               "is an INTERNAL label-permutation null on the model's own field, "
               "not a comparison against any measured number, so it is not an "
               "external anchor and is not scored here. The sink does not improve "
               "the measured human-glioma anchor and ships DEFAULT-OFF "
               "(results/immune_sink_impact.md).",
        conditional=(),
        values={"n_drugs_with_measured_uptake_ratio": 0, "n_drugs": 3},
    )]


# --------------------------------------------------------------------------- #
# FULL anchors — need real tissue under data_cache/, and mostly the PDE
# --------------------------------------------------------------------------- #
def _import_experiment(name: str):
    """Import a driver under ``experiments/`` as a module.

    ``experiments`` is an implicit namespace package once the repository root is
    importable, which is exactly how ``experiments/kill_calibration_breast_report.py``
    already imports its sibling. Importing (rather than executing) a driver runs
    only its module body — the ``if __name__ == "__main__"`` block, and therefore
    every write into ``results/``, is not reached.
    """
    if REPO_ROOT not in sys.path:
        sys.path.insert(0, REPO_ROOT)
    import importlib

    return importlib.import_module(f"experiments.{name}")


def _run_breast_transvascular(exchange: bool = False) -> List[Row]:
    """Predicted interstitial:plasma AUC ratio on the real breast Xenium section
    vs MEASURED in-vivo human microdialysis.

    This is the project's flagship external transport falsification. The engine
    quantity and the measured quantity are the same physical thing: the PDE's
    state variable is free interstitial concentration, which is what a
    microdialysis probe samples.

    Only the two ANCHORED drugs are simulated (the driver's own panel run also
    covers ~11 unanchored drugs for spread, which is a different question and is
    not needed for a scorecard).
    """
    import numpy as np

    from ..core import PDConfig, PDEConfig, PKConfig, SimConfig
    from ..reconstruct.real_pipeline import apply_drug_sensitivity, reconstruct_real
    from ..sim.engine import simulate_response
    from . import delivery_metrics

    bkp = _import_experiment("breast_kp_anchor")
    section = _repo_path(bkp.SECTION_PATH)
    if not os.path.isdir(section):
        raise NotTestable(f"requires the real breast Xenium cache at "
                          f"{bkp.SECTION_PATH} (absent)")

    micro = reconstruct_real(source="xenium", sample=section, dx_um=bkp.DX_UM,
                             oxygen_mode="pde", state_mode="oxygen")
    sens0 = np.array(micro.cells.sensitivity, copy=True)
    state0 = np.array(micro.cells.state, copy=True)
    meta0 = dict(micro.metadata)

    rows: List[Row] = []
    for name, (mean, se, n, pmid) in sorted(bkp.ANCHOR.items()):
        # Reset the shared reconstruction so runs cannot contaminate each other
        # (the driver does exactly this).
        micro.cells.sensitivity = np.array(sens0, copy=True)
        micro.cells.state = np.array(state0, copy=True)
        micro.metadata = dict(meta0)
        apply_drug_sensitivity(micro, name, mode="uniform")

        dp, mol = bkp.compile_drug(name)
        cfg = SimConfig(
            drug=dp, pk=PKConfig(),
            pde=PDEConfig(dt=10.0, total_time=bkp.TOTAL_TIME_S, convection=True,
                          transvascular_exchange=bool(exchange)),
            pd=PDConfig(), seed=0,
        )
        res = simulate_response(micro, cfg)
        summ = delivery_metrics.tissue_auc_summary(micro, res)
        plasma_auc = float(res.extras["plasma_mean"]) * bkp.TOTAL_TIME_S
        ratio = summ["mean_tumor_auc"] / plasma_auc

        lo, hi = mean - se, mean + se
        status = classify(ratio, lo, hi)
        cond = ["alpha-fit", "untreated-section"]
        law = ("one-way leak S = PS*Cp (repo default)" if not exchange
               else "gradient-driven exchange S = PS*(Cp - C)")
        if exchange:
            cond += ["ps-scale", "brain-permeability-in-breast"]
        rows.append(Row(
            group="breast-transvascular", label=name, status=status,
            scope=SCOPE_PDE,
            quantity="mean_tumor_auc / (plasma_mean * total_time) — free "
                     "interstitial drug AUC at tumour cells over plasma AUC — vs "
                     "MEASURED in-vivo microdialysis interstitial:plasma AUC in "
                     "primary human breast tumour",
            predicted=f"{ratio:.4f}  ({ratio / mean:.2f}x measured)",
            measured=f"{mean:.2f} +/- {se:.2f} (mean +/- SE, n={n})",
            pmids=(str(pmid),),
            bar="inside measured mean +/- SE",
            note=f"delivery law: {law}",
            conditional=tuple(cond),
            values={"ratio": float(ratio), "measured_mean": float(mean),
                    "measured_se": float(se), "n": int(n),
                    "fold_of_measured": float(ratio / mean),
                    "transvascular_exchange": bool(exchange)},
        ))
    return rows


def _kill_rows(group: str, tissue_note: str, tests: Dict, per_drug: Dict,
               conditional: Sequence[str]) -> List[Row]:
    """Shared renderer for both 3-D:2-D ratio tests (GBM and breast).

    ``tests`` is ``kill_calibration.analyse``'s ``A["tests"]`` shape, i.e. the
    output of the frozen ``ratio_test``; nothing is re-scored here.
    """
    rows: List[Row] = []
    for name in sorted(tests):
        t = tests[name]
        d = per_drug.get(name, {})
        status = PASS if t["verdict"] == PASS else FAIL
        matched = t.get("scope", "")
        note = f"literature scope used: {matched}"
        if matched == "solid":
            note += (" — NO drug-matched 3-D:2-D band exists, so this falls back "
                     "to the broad cross-tumour-type cytotoxic band")
        if t["verdict"] != PASS:
            note += f". {t['direction']}"
        note += f". {tissue_note}"
        rows.append(Row(
            group=group, label=name, status=status, scope=SCOPE_PDE,
            quantity="predicted 3-D:2-D kill ratio on the real section "
                     "(delivery x cell-state resistance) vs MEASURED "
                     "spheroid-vs-monolayer IC50 fold shift",
            predicted=f"ratio {t['predicted_ratio']:.4f} "
                      f"(implied IC50 fold {t.get('implied_fold', float('nan')):.2f}x)",
            measured=f"ratio band {t['lit_ratio_lo']:.4f}-{t['lit_ratio_hi']:.4f} "
                     f"(measured fold {t['fold_lo']:.2g}-{t['fold_hi']:.2g}x)",
            pmids=tuple(_bridge_pmids(matched)),
            bar="predicted ratio inside the literature-derived ratio band "
                "(kill_calibration.ratio_test)",
            note=note, conditional=tuple(conditional),
            values={"predicted_ratio": float(t["predicted_ratio"]),
                    "lit_ratio_lo": float(t["lit_ratio_lo"]),
                    "lit_ratio_hi": float(t["lit_ratio_hi"]),
                    "implied_fold": float(t.get("implied_fold", float("nan"))),
                    "scope": matched,
                    "ic50_uM": float(d.get("ic50_uM", float("nan")))},
        ))
    return rows


def _bridge_pmids(scope: str) -> Tuple[str, ...]:
    """PMIDs of the literature rows the named bridge scope aggregates."""
    from .kill_calibration import LIT_2D3D

    out: List[str] = []
    for r in LIT_2D3D:
        if scope in r.get("scope", ()) and r.get("fold") is not None:
            p = str(r.get("pmid", ""))
            if p and p not in out:
                out.append(p)
    return tuple(out)


def _run_gbm_kill_3d2d() -> List[Row]:
    """GBM 3-D:2-D kill ratio on the real Xenium GBM section."""
    from . import kill_calibration as KC

    if not os.path.isdir(_repo_path("data_cache", "xenium_gbm")):
        raise NotTestable("requires the real Xenium GBM cache at "
                          "data_cache/xenium_gbm (absent)")
    # verify_engine=False: the full-engine cross-check is a self-consistency
    # check of the fast path, not an external anchor, and it is the single most
    # expensive call in the module.
    res = KC.run(use_cache=True, verbose=False, verify_engine=False)
    A = KC.analyse(res)
    return _kill_rows(
        "gbm-kill-3d2d",
        "Real untreated Xenium GBM section, dx = 12 um, 6 h window.",
        A["tests"], res["per_drug"],
        ("alpha-fit", "gdsc-potency", "band-wider-than-assay", "untreated-section"),
    )


def _run_breast_kill_3d2d() -> List[Row]:
    """Breast 3-D:2-D kill ratio on the real Janesick Xenium breast section.

    Drugs with NO published 3-D-vs-2-D comparison in any tumour type are
    reported NOT TESTABLE rather than scored against a non-drug-matched band —
    an honest null, exactly as ``results/kill_calibration_breast.md`` reports it.
    """
    kcb = _import_experiment("kill_calibration_breast")
    rep = _import_experiment("kill_calibration_breast_report")

    if not os.path.isdir(kcb.BREAST_SAMPLE_PATH):
        raise NotTestable("requires the real breast Xenium cache at "
                          "data_cache/xenium_breast_janesick/rep1 (absent)")
    res = kcb.run(verbose=False, verify_engine=False, use_cache=True)
    section = rep.build_ratio_section(res)

    rows: List[Row] = []
    for name in sorted(section):
        row = section[name]
        h = row["headline"]
        if not row.get("testable"):
            rows.append(Row(
                group="breast-kill-3d2d", label=name, status=NOT_TESTABLE,
                scope=SCOPE_PDE,
                quantity="predicted 3-D:2-D kill ratio vs MEASURED "
                         "spheroid-vs-monolayer IC50 fold shift",
                predicted=f"ratio {h['ratio_exact']:.4f}",
                measured="NONE",
                pmids=(), bar="n/a",
                reason="NO MEASURED TARGET EXISTS. A structured PubMed search "
                       "found ZERO published 3-D-vs-2-D spheroid/monolayer "
                       "potency comparisons for this drug in ANY tumour type. "
                       "Its predicted ratio is reported for completeness and is "
                       "deliberately NOT scored against the broad, "
                       "non-drug-matched band "
                       "(results/kill_calibration_breast.md).",
                conditional=("alpha-fit", "gdsc-potency"),
                values={"predicted_ratio": float(h["ratio_exact"])},
            ))
            continue
        t = row["test"]
        t = dict(t)
        t.setdefault("scope", f"drug:{name}")
        t["fold_lo"] = row["band"]["fold_lo"]
        t["fold_hi"] = row["band"]["fold_hi"]
        t["implied_fold"] = h.get("implied_ic50_fold", float("nan"))
        rows.extend(_kill_rows(
            "breast-kill-3d2d",
            "Real untreated Janesick Xenium breast section (PMID 38114474), "
            "dx = 12 um, 6 h window."
            + ("" if row.get("band_is_drug_matched")
               else " Single-row band (n=1): no spread to test against."),
            {name: t}, res["per_drug"],
            ("alpha-fit", "gdsc-potency", "untreated-section"),
        ))
    return rows


def _run_vascular_morphometry() -> List[Row]:
    """Reconstructed vascular geometry vs measured human GBM morphometry.

    Reconstruction only — no drug, no PDE — so ``scope`` is ``geometry``. These
    anchors constrain the tissue the transport runs ON, which is a real and
    separate claim from the transport itself.
    """
    from . import vascular_validation as VV

    if not os.path.isdir(_repo_path("data_cache", "xenium_gbm")):
        raise NotTestable("requires the real Xenium GBM cache at "
                          "data_cache/xenium_gbm (absent)")
    from ..reconstruct.real_pipeline import reconstruct_real

    micro = reconstruct_real(source="xenium", dx_um=VV.DX_UM)
    res = VV.measure(micro)
    prim = res["primary"]
    kde_abs = {a["threshold"]: a["frac_geom_mask"] for a in res["kde"]["absolute"]}
    R = VV.LIT_RANGES

    rows: List[Row] = []

    def _add(label, value, band_key, quantity, pred_s, meas_s, pmids, note="",
             values=None):
        lo, hi = R[band_key]
        rows.append(Row(
            group="vascular-morphometry", label=label,
            status=classify(value, lo, hi), scope=SCOPE_GEOMETRY,
            quantity=quantity, predicted=pred_s, measured=meas_s,
            pmids=pmids,
            bar=f"inside the module's frozen LIT_RANGES['{band_key}'] "
                f"= {lo:.4g}..{hi:.4g}",
            note=note,
            values=dict(values or {}, value=float(value),
                        band=[float(lo), float(hi)]),
        ))

    # 1. microvessel density — the primary anchor.
    mvd = float(prim["mvd_per_mm2"])
    _add("microvessel density", mvd, "mvd_gbm_per_mm2",
         "reconstructed vessel-cluster density (per mm^2) vs MEASURED "
         "whole-slide CD31 microvessel density in human GBM",
         f"{mvd:.1f} /mm^2",
         f"{VV.GBM_MVD_MEAN:.0f} +/- {VV.GBM_MVD_SD:.0f} /mm^2 (n=20 patients)",
         ("27583442",),
         "not reproduced on the held-out patient either "
         "(results/heldout_patient.md)")

    # 2. Intercapillary distance, the one DIRECTLY MEASURED spacing anchor
    #    (Voronoi, not derived from an MVD via a Poisson assumption).
    nn = float(prim["vessel_nn_median"])
    _add("intercapillary distance (directly measured)", nn,
         "intercapillary_yoshii_glioma_um",
         "median reconstructed vessel-to-nearest-vessel distance (um) vs "
         "MEASURED Voronoi intercapillary distance in the proliferating area "
         "of human glioma",
         f"{nn:.1f} um",
         "77.8 +/- 15.9 um (n=22 patients; normal white matter 123 +/- 4.2)",
         ("3359450",),
         "the only spacing anchor here that is a direct measurement rather than "
         "a Poisson consequence of a published density")

    # 3. Tumour cells beyond the diffusion-limited cord radius.
    frac_pct = 100.0 * float(prim["tumor_to_vessel_frac_gt_100"])
    lo, hi = R["frac_gt_100_pct"]
    _add("tumour cells >100 um from a vessel", frac_pct, "frac_gt_100_pct",
         "fraction of reconstructed tumour cells further than 100 um from the "
         "nearest vessel vs the Poisson fraction implied by the MEASURED GBM "
         "microvessel density",
         f"{frac_pct:.1f}%", f"{lo:.1f}-{hi:.1f}% implied by 65 +/- 30 /mm^2",
         ("27583442", "13304213"),
         "100 um is the load-bearing length because it is the lower edge of "
         "Thomlinson & Gray's 100-180 um diffusion-limited viable-cord band")

    # 4. Vascular area fraction — NOT TESTABLE, and the owning module says so.
    #    The reconstructed value depends on which KDE threshold is chosen, and
    #    the threshold is not independently determined, so agreement with the
    #    measured 2.1% is a choice rather than a result.
    area_pct = 100.0 * float(kde_abs[0.3])
    lo, hi = R["vascular_area_pct"]
    rows.append(Row(
        group="vascular-morphometry", label="vascular area fraction",
        status=NOT_TESTABLE, scope=SCOPE_GEOMETRY,
        quantity="reconstructed vascular area fraction (%) vs MEASURED "
                 "whole-section CD34+ area fraction in IDH-wildtype GBM",
        predicted=f"{area_pct:.2f}% (KDE threshold 0.3; the 0.2/0.5 thresholds "
                  f"give {100 * kde_abs[0.2]:.2f}% / {100 * kde_abs[0.5]:.2f}%)",
        measured=f"{lo:.1f}-{hi:.1f}% (2.1 +/- 0.3%, n=30)",
        pmids=("36611806",), bar="n/a",
        reason="THRESHOLD-CHOSEN, NOT A VALIDATION. The reconstructed value is "
               "set by which KDE threshold is picked, and nothing determines "
               "that threshold independently, so landing on the measured value "
               "is a choice rather than a result. results/vascular_validation.md "
               "annotates this row with exactly that wording. It also DOES NOT "
               "REPLICATE on the held-out patient (5.23x, "
               "results/heldout_patient.md).",
        values={"area_pct_at_0.3": area_pct,
                "area_pct_at_0.5": 100 * float(kde_abs[0.5]),
                "area_pct_at_0.2": 100 * float(kde_abs[0.2]),
                "band": [float(lo), float(hi)]},
    ))
    return rows


# --------------------------------------------------------------------------- #
# The registry
# --------------------------------------------------------------------------- #
def build_registry(exchange: bool = False) -> List[AnchorGroup]:
    """The full anchor registry, quick tier first.

    ``exchange`` selects the delivery law for the breast transvascular anchor:
    ``False`` (the repo default) is the shipping one-way leak; ``True`` is the
    opt-in gradient-driven exchange, whose every number is conditional on the
    uncalibrated ``PS_SCALE_S_INV`` (see :data:`CONDITIONALS`).
    """
    return [
        AnchorGroup(
            key="gbm-penetration",
            title="GBM perivascular penetration depth (fit + 2 out-of-sample)",
            tissue="GBM", tier=TIER_QUICK, scope=SCOPE_ANALYTIC,
            # run_penetration drives calibration.calibrate_alpha and covers a
            # superset of calibration.revalidate's comparison: the same frozen
            # alpha, the same mitoxantrone anchor, plus topotecan (a chemically
            # DISTINCT out-of-sample drug, so a more independent test than
            # revalidate's single same-class one).
            source="tmesim/validation/pk_multidrug.py::run_penetration "
                   "(drives calibration.calibrate_alpha; same anchors as "
                   "calibration.revalidate, plus topotecan)",
            report="results/calibration.md, results/pk_multidrug.md",
            runner=_run_perivascular_penetration, est="~10 s",
        ),
        AnchorGroup(
            key="bbb-ordering",
            title="bbb_penetrance ordering vs measured brain/CSF:plasma ordering",
            tissue="GBM", tier=TIER_QUICK, scope=SCOPE_ANALYTIC,
            source="tmesim/validation/pk_multidrug.py::run_partition",
            report="results/pk_multidrug.md",
            runner=_run_bbb_ordering, est="~5 s",
        ),
        AnchorGroup(
            key="glioma-partition",
            title="Human glioma tissue:plasma partition, 6 measured drugs",
            tissue="GBM", tier=TIER_QUICK, scope=SCOPE_ANALYTIC,
            source="tmesim/validation/glioma_regional_pk.py::run_scorecard",
            report="results/glioma_regional_pk.md",
            runner=_run_glioma_partition, est="~2 s",
            needs=("data_cache/hunt/autoradiography_microdialysis/parsed/"
                   "regional_intratumoral_PK_anchors_literature.csv",
                   "data_cache/hunt/autoradiography_microdialysis/parsed/"
                   "oh2026_neurooncol_TableS3_per_sample_tumor_drug_ratios.csv"),
        ),
        AnchorGroup(
            key="gbm-tmz-exvivo",
            title="Temozolomide potency extrapolated out of range vs measured "
                  "ex-vivo viability",
            tissue="GBM", tier=TIER_QUICK, scope=SCOPE_ANALYTIC,
            source="tmesim/validation/gbm_tmz_exvivo.py::test_a_potency_axis",
            report="results/gbm_tmz_exvivo.md",
            runner=_run_tmz_exvivo, est="~10 s",
            needs=("data_cache/gdsc_glioma_dose_response.csv",
                   "data_cache/hunt/geo_organoid_screens"),
        ),
        AnchorGroup(
            key="gist-imatinib",
            title="GIST imatinib tumour:normal partition vs measured MALDI-MSI",
            tissue="GIST", tier=TIER_QUICK, scope=SCOPE_ANALYTIC,
            source="tmesim/validation/gist_imatinib.py::score",
            report="results/gist_imatinib.md",
            runner=_run_gist_imatinib, est="~15 s",
            needs=("data_cache/hunt/maldi_msi_drug",
                   "data_cache/hunt/verified_0/gist_spatial_plus_clinical.csv"),
        ),
        AnchorGroup(
            key="ps-scale",
            title="Implicit unit scale PS_SCALE_S_INV vs measured brain PA product",
            tissue="GBM", tier=TIER_QUICK, scope=SCOPE_ANALYTIC,
            source="tmesim/validation/ps_scale.py::calibrate_ps_scale",
            report="results/transvascular_exchange_revalidation.md sec. 6.2",
            runner=_run_ps_scale, est="<1 s",
        ),
        AnchorGroup(
            key="immune-sink",
            title="Immune (TAM) drug sequestration magnitude",
            tissue="GBM", tier=TIER_QUICK, scope=SCOPE_PDE,
            source="tmesim/validation/immune_sink_impact.py",
            report="results/immune_sink_impact.md",
            runner=_run_immune_sink, est="<1 s",
        ),
        AnchorGroup(
            key="breast-transvascular",
            title="Breast interstitial:plasma AUC vs measured human microdialysis",
            tissue="breast", tier=TIER_FULL, scope=SCOPE_PDE,
            source="experiments/breast_kp_anchor.py",
            report="results/breast_kp_validation.md",
            runner=lambda: _run_breast_transvascular(exchange=exchange),
            est="~4-8 min (reconstruction + 2 PDE solves, no cache)",
            needs=("data_cache/xenium_breast_janesick/rep1",),
        ),
        AnchorGroup(
            key="gbm-kill-3d2d",
            title="GBM 3-D:2-D kill ratio vs measured spheroid-vs-monolayer shift",
            tissue="GBM", tier=TIER_FULL, scope=SCOPE_PDE,
            source="tmesim/validation/kill_calibration.py::run + analyse",
            report="results/kill_calibration.md",
            runner=_run_gbm_kill_3d2d,
            est="~2-4 min warm (data_cache/kill_calibration_runs.json), longer cold",
            needs=("data_cache/xenium_gbm", "data_cache/gdsc_glioma_dose_response.csv"),
        ),
        AnchorGroup(
            key="breast-kill-3d2d",
            title="Breast 3-D:2-D kill ratio vs measured spheroid-vs-monolayer shift",
            tissue="breast", tier=TIER_FULL, scope=SCOPE_PDE,
            source="experiments/kill_calibration_breast.py "
                   "(+ kill_calibration_breast_report.build_ratio_section)",
            report="results/kill_calibration_breast.md",
            runner=_run_breast_kill_3d2d,
            est="~3-6 min warm (data_cache/kill_calibration_breast_runs.json)",
            needs=("data_cache/xenium_breast_janesick/rep1",
                   "data_cache/gdsc_breast_dose_response.csv"),
        ),
        AnchorGroup(
            key="vascular-morphometry",
            title="Reconstructed vascular geometry vs measured human GBM morphometry",
            tissue="GBM", tier=TIER_FULL, scope=SCOPE_GEOMETRY,
            source="tmesim/validation/vascular_validation.py::measure",
            report="results/vascular_validation.md",
            runner=_run_vascular_morphometry, est="~30 s",
            needs=("data_cache/xenium_gbm",),
        ),
    ]


#: Validations that were inventoried and DELIBERATELY not registered, with the
#: reason. Printed by ``vivum validate --list`` so the omissions are auditable
#: rather than invisible.
OUT_OF_SCOPE: Tuple[Tuple[str, str], ...] = (
    ("tmesim/validation/efflux_impact.py, ion_trapping_impact.py, "
     "uptake_impact.py, immune_sink_impact.py::glioma_scorecard",
     "ablations: each re-runs glioma_regional_pk.run_scorecard under a modified "
     "model. Registering them would count the SAME six measurements up to four "
     "times and inflate the anchor count."),
    ("tmesim/validation/pharmacology_real.py (GDSC vs HGCC concordance, "
     "rho = 0.827, PMID 32668248)",
     "real-vs-real corroboration with no model in the loop: it validates the "
     "potency INPUT, not any model output."),
    ("tmesim/validation/vessel_profile.py (doxorubicin d_half on real tissue "
     "vs PMID 16361566)",
     "NOT independent: alpha was fit to this same anchor, so the agreement is "
     "in-sample. The module's own report says so and directs you to the "
     "drug/source-null ratio instead."),
    ("tmesim/validation/numerical_rigor.py, dimensionality.py analytic identities, "
     "reconstruction_fidelity.py",
     "internal consistency (mass conservation, grid/dt convergence, "
     "cross-platform agreement). No external measured number."),
    ("tmesim/validation/mgmt_scalar.py",
     "a FIT (lambda = 0.3429 on n=17), not a test of one."),
    ("tmesim/validation/state_validation.py, cellstate_validation.py, "
     "oxygen_demotion",
     "model-vs-marker correlation WITHIN the same section; no external published "
     "value is being reproduced."),
    ("tmesim/validation/moncayo_explant.py (PMID 41844641)",
     "structurally capped: 89.9% of measured variance is within-slice "
     "ROI-to-ROI and the ROI table carries no X/Y coordinates, so a "
     "slice-resolved prediction cannot reach it. Not scoreable pass/fail."),
    ("tmesim/validation/heldout_patient.py (PMIDs 38723085, 36611806, 27583442)",
     "a real held-out external-anchor test, but ~18 min and 10 PDE solves — too "
     "slow for a suite meant to be run routinely. Run it directly."),
    ("results/oxygen_calibration_breast.md / oxygen_calibration_landscape.md",
     "the mmHg <-> dimensionless-O mapping was never derived for these tissues "
     "(REFERENCE_PO2_MMHG = 35.0 is a normal-human-BRAIN anchor), so nothing is "
     "scoreable. Reported in those files as an honest null."),
    ("results/breast_tissue_permeability.md",
     "a falsified FIX to an anchor already registered here "
     "(breast-transvascular), not a separate measurement."),
)


# --------------------------------------------------------------------------- #
# Suite driver
# --------------------------------------------------------------------------- #
def _missing(needs: Sequence[str]) -> List[str]:
    return [p for p in needs if not os.path.exists(_repo_path(p))]


def run_suite(quick: bool = False, only: Optional[Sequence[str]] = None,
              exchange: bool = False,
              on_progress: Optional[Callable[[AnchorGroup], None]] = None,
              ) -> Dict:
    """Run the registered anchors and return a machine-readable scorecard.

    Parameters
    ----------
    quick
        Run only the ``quick``-tier anchors (analytic / literature-table, no
        real tissue, seconds). Full-tier anchors become ``DEFERRED`` rows with a
        stated reason — never PASS, and never silently dropped.
    only
        Restrict to these group keys.
    exchange
        Delivery law for the breast transvascular anchor (see
        :func:`build_registry`).
    on_progress
        Called with each group just before it runs.
    """
    groups = build_registry(exchange=exchange)
    if only:
        wanted = {k.strip() for k in only if k.strip()}
        unknown = wanted - {g.key for g in groups}
        if unknown:
            raise KeyError(f"unknown anchor key(s): {sorted(unknown)}. "
                           f"Known: {sorted(g.key for g in groups)}")
        groups = [g for g in groups if g.key in wanted]

    rows: List[Row] = []
    group_meta: List[Dict] = []
    t_start = time.time()
    for g in groups:
        if on_progress is not None:
            on_progress(g)
        t0 = time.time()
        deferred = quick and g.tier != TIER_QUICK
        missing = _missing(g.needs)
        if deferred:
            grows = [Row(
                group=g.key, label="(whole anchor)", status=DEFERRED,
                scope=g.scope, quantity=g.title, predicted="-", measured="-",
                reason=f"DEFERRED by --quick: this anchor needs real tissue from "
                       f"data_cache/ and {g.est}. It is NOT a pass. Re-run "
                       f"without --quick.",
            )]
        elif missing:
            grows = [Row(
                group=g.key, label="(whole anchor)", status=NOT_TESTABLE,
                scope=g.scope, quantity=g.title, predicted="-", measured="-",
                reason="REQUIRED DATA ABSENT: " + ", ".join(missing),
            )]
        else:
            try:
                grows = list(g.runner())
            except NotTestable as exc:
                grows = [Row(
                    group=g.key, label="(whole anchor)", status=NOT_TESTABLE,
                    scope=g.scope, quantity=g.title, predicted="-", measured="-",
                    reason=str(exc),
                )]
            except Exception as exc:  # noqa: BLE001 - a crash is not a pass
                grows = [Row(
                    group=g.key, label="(whole anchor)", status=NOT_TESTABLE,
                    scope=g.scope, quantity=g.title, predicted="-", measured="-",
                    reason=f"ANCHOR CRASHED and therefore produced no verdict: "
                           f"{type(exc).__name__}: {exc}",
                )]
        rows.extend(grows)
        group_meta.append({
            "key": g.key, "title": g.title, "tissue": g.tissue, "tier": g.tier,
            "scope": g.scope, "source": g.source, "report": g.report,
            "seconds": round(time.time() - t0, 2), "n_rows": len(grows),
        })

    counts = {s: sum(1 for r in rows if r.status == s) for s in STATUS_ORDER}
    cited = sorted({c for r in rows for c in r.conditional})
    n_scored = counts[PASS] + counts[FAIL]
    n_analytic_scored = sum(
        1 for r in rows
        if r.status in (PASS, FAIL) and r.scope != SCOPE_PDE)
    return {
        "mode": "quick" if quick else "full",
        "exchange": bool(exchange),
        "counts": counts,
        "n_rows": len(rows),
        "n_scored": n_scored,
        "n_scored_not_pde": n_analytic_scored,
        "exit_code": 1 if counts[FAIL] else 0,
        "seconds": round(time.time() - t_start, 2),
        "groups": group_meta,
        "rows": [r.to_dict() for r in rows],
        "conditionals": {k: CONDITIONALS[k] for k in cited},
    }
