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
``geometry`` means the row scores the reconstructed TISSUE — vessel density,
spacing, area fraction, cellularity — with no drug and no solver anywhere in
it; those anchors constrain what the transport runs ON, which is a real claim
but a different one. ``analytic`` means the prediction is a closed-form
expression over molecular descriptors and the solver was never invoked —
verified for
``glioma_regional_pk`` and ``gist_imatinib`` by call-counter instrumentation in
``results/transvascular_exchange_revalidation.md`` §6.1, which states the rule
directly: *"These modules must be reported as UN-TESTABLE, not as unaffected.
'Robust under exchange' and 'never evaluates the delivery law' are very
different statements and must not be pooled in a coverage count."* An analytic
row still carries a real PASS/FAIL against its measurement — it just cannot be
cited as evidence about the solver, and the printed scorecard says so.

Physics provenance
------------------
Every row carries a ``physics`` string naming the transvascular delivery law
the number was produced under, and the suite result carries a ``physics``
block for the run as a whole. This is not decoration: the breast
interstitial:plasma anchor differs by roughly **10x** between
``transvascular_exchange=False`` (one-way leak, ``S = PS*Cp``, what the engine
ships) and ``True`` (gradient-driven, ``S = PS*(Cp - C)``) — 5-FU lands at
12.78x of measured under the former and near 1.33x under the latter. A
scorecard that printed one number without naming the law would be ambiguous by
an order of magnitude, which is precisely the failure mode this project keeps
finding in its own results files.

Three rules enforce it:

* the shipping default is READ from ``PDEConfig`` (:func:`repo_default_exchange`)
  rather than hardcoded, so it cannot desynchronise;
* the breast runner reads the law back off the config that was actually solved
  and reports a FAIL provenance row if it differs from the one requested;
* :func:`_stamp_physics` DOWNGRADES any ``PDE``-scoped row that carries a
  verdict without declaring a law to ``NOT TESTABLE``. An unlabelled solver
  number is uninterpretable, not merely imprecise, and must never pass quietly.

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
    #: The physics the number was produced under. NEVER left blank on a row
    #: whose scope is PDE — see :func:`physics_label`.
    physics: str = ""

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
            "physics": self.physics,
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


# --------------------------------------------------------------------------- #
# Physics provenance
# --------------------------------------------------------------------------- #
#: The two transvascular delivery laws, and what each one IS. The breast
#: interstitial:plasma anchor differs by ~10x between them (5-FU lands at
#: 12.78x measured under the one-way law and ~1.33x under exchange), so a
#: scorecard that prints one number without naming the law is ambiguous by an
#: order of magnitude. That is the same failure this project keeps finding in
#: its own results files, so every PDE row carries its law explicitly.
PHYSICS_LAWS: Dict[bool, str] = {
    False: "transvascular_exchange=False — one-way leak, S = PS*Cp "
           "(source never depends on local tissue concentration)",
    True: "transvascular_exchange=True — gradient-driven exchange, "
          "S = PS*(Cp - C) (opt-in; every number CONDITIONAL on PS_SCALE_S_INV)",
}

#: Rows whose prediction never went through the solver. Not "unaffected by the
#: delivery law" — simply not a statement about it at all.
PHYSICS_NOT_APPLICABLE = "n/a — transport solver never invoked for this row"


def repo_default_exchange() -> bool:
    """The delivery law the engine SHIPS with, read from PDEConfig itself.

    Read rather than hardcoded so that a change to the shipping default cannot
    silently desynchronise from what this scorecard claims the default is.
    """
    from ..core import PDEConfig

    return bool(PDEConfig().transvascular_exchange)


def physics_label(exchange: Optional[bool]) -> str:
    """Human-readable provenance for the physics a row ran under.

    ``None`` means the solver was not invoked at all.
    """
    if exchange is None:
        return PHYSICS_NOT_APPLICABLE
    law = PHYSICS_LAWS[bool(exchange)]
    if bool(exchange) != repo_default_exchange():
        law += "  [NOT THE SHIPPING DEFAULT]"
    return law


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
               "module's permutation test is an INTERNAL null on the model's "
               "own field, not a comparison against any measured number, so it "
               "is not an external anchor and is not scored here whatever its "
               "value. Do not quote its magnitude from memory: "
               "results/immune_sink_impact.md carries a CORRECTION NOTICE "
               "demoting the original label-permutation null (that null "
               "under-clumps the field it is compared against) and reporting a "
               "much smaller excess against spatially-aware nulls, part of it "
               "attributable to a transport gradient with depth into tissue "
               "rather than to immune biology. The sink does not improve the "
               "measured human-glioma anchor and ships DEFAULT-OFF.",
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
        # Read the delivery law back OFF THE CONFIG THAT WAS ACTUALLY SOLVED,
        # rather than trusting the flag we passed in. This anchor moves by
        # ~10x between the two laws, so a silent desync between "what the
        # scorecard says it ran" and "what ran" would be worth an order of
        # magnitude — the exact failure this command exists to catch.
        ran_exchange = bool(cfg.pde.transvascular_exchange)
        if ran_exchange != bool(exchange):
            rows.append(Row(
                group="breast-transvascular", label=name, status=FAIL,
                scope=SCOPE_PDE,
                quantity="physics-provenance check: requested delivery law vs "
                         "the law the solved config actually carried",
                predicted=f"solved with transvascular_exchange={ran_exchange}",
                measured=f"requested transvascular_exchange={bool(exchange)}",
                pmids=(str(pmid),),
                bar="the solved config must carry the requested delivery law",
                note="NOT a comparison against measurement — a provenance "
                     "failure. The anchor differs by ~10x between the two laws, "
                     "so a number produced under the wrong one is not merely "
                     "imprecise, it is mislabelled. Refusing to score.",
                physics=physics_label(ran_exchange),
            ))
            continue

        res = simulate_response(micro, cfg)
        summ = delivery_metrics.tissue_auc_summary(micro, res)
        plasma_auc = float(res.extras["plasma_mean"]) * bkp.TOTAL_TIME_S
        ratio = summ["mean_tumor_auc"] / plasma_auc

        lo, hi = mean - se, mean + se
        status = classify(ratio, lo, hi)
        cond = ["alpha-fit", "untreated-section"]
        law = ("one-way leak S = PS*Cp (repo default)" if not ran_exchange
               else "gradient-driven exchange S = PS*(Cp - C)")
        if ran_exchange:
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
            note=f"delivery law: {law}. This anchor moves by roughly 10x "
                 f"between the two laws, so the number above is only "
                 f"interpretable together with the physics line.",
            conditional=tuple(cond),
            physics=physics_label(ran_exchange),
            values={"ratio": float(ratio), "measured_mean": float(mean),
                    "measured_se": float(se), "n": int(n),
                    "fold_of_measured": float(ratio / mean),
                    "transvascular_exchange": ran_exchange,
                    "repo_default_transvascular_exchange":
                        repo_default_exchange()},
        ))
    return rows


def _kill_rows(group: str, tissue_note: str, tests: Dict, per_drug: Dict,
               conditional: Sequence[str],
               pmids: Optional[Dict[str, Sequence[str]]] = None) -> List[Row]:
    """Shared renderer for both 3-D:2-D ratio tests (GBM and breast).

    ``tests`` is ``kill_calibration.analyse``'s ``A["tests"]`` shape, i.e. the
    output of the frozen ``ratio_test``; nothing is re-scored here.

    ``pmids`` overrides the scope-derived PMID lookup per drug. The breast
    report builds several of its bands from report-local row tuples that carry
    no ``scope``, so the GBM path's scope lookup cannot find them; without the
    override those rows would print "PMID n/a", which for an external-anchor
    scorecard is not a cosmetic defect.
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
            pmids=tuple((pmids or {}).get(name) or _bridge_pmids(matched)),
            # The kill drivers never override PDEConfig, so they run under
            # whatever the engine ships. Read it rather than assume it.
            physics=physics_label(repo_default_exchange()),
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
            pmids={name: _rows_pmids(row.get("specific_bridge"))},
        ))
    return rows


def _rows_pmids(bridge: Optional[Dict]) -> Tuple[str, ...]:
    """PMIDs of the literature rows a ``kill_calibration.literature_bridge``
    result aggregated, in order and de-duplicated."""
    out: List[str] = []
    for r in (bridge or {}).get("rows", ()):
        p = str(r.get("pmid", "") or "")
        if p and p not in out:
            out.append(p)
    return tuple(out)


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

    # 5. The SAME measured CD34+ area fraction, but against the THRESHOLD-FREE
    #    reconstruction. This is a genuinely scoreable row and row 4 is not, so
    #    they must not be pooled: vessel_geometry.vascular_area_fraction paints
    #    a sub-grid vessel footprint from data-derived radii (20 EC/mm, PMID
    #    38723085; mean degree 1.5, PMID 37582180 — literature inputs that are
    #    independent of the 2.1% target), with no threshold to pick and no
    #    hardcoded radius. segmentation_report scores it against its own frozen
    #    PUBLISHED range, so the bar below is that module's, not this one's.
    from ..reconstruct.vessel_network import PUBLISHED, segment_vessels
    from ..reconstruct.vessel_network import segmentation_report

    vessels = segment_vessels(micro.cells, micro.domain)
    seg = segmentation_report(vessels, micro.domain)
    af_row = next(r for r in seg["rows"]
                  if r["key"] == "vascular_area_fraction_pct")
    slo, shi = af_row["published_range"]
    rows.append(Row(
        group="vascular-morphometry",
        label="vascular area fraction (threshold-free)",
        status=PASS if af_row["verdict"] == "PASS" else FAIL,
        scope=SCOPE_GEOMETRY,
        quantity="grid-independent, data-derived sub-grid vascular area "
                 "fraction (%) vs MEASURED whole-section CD34+ area fraction "
                 "in IDH-wildtype GBM",
        predicted=f"{af_row['reconstructed']:.2f}%  "
                  f"({af_row['ratio']:.2f}x measured; module verdict "
                  f"{af_row['verdict']}). Legacy fixed-radius binary field gave "
                  f"{seg['legacy_binary_area_fraction_pct']:.2f}%",
        measured=f"{af_row['published']:.1f}% (range {slo:.1f}-{shi:.1f}%, n=30)",
        pmids=(str(af_row["pmid"]),),
        bar=f"inside the module's frozen PUBLISHED range "
            f"{slo:.1f}-{shi:.1f}% (vessel_network.segmentation_report)",
        note="This row and the KDE row above score the SAME measured number "
             "against two different reconstructions and must not be pooled: "
             "this one has no threshold to choose, so it can fail. The vessel "
             "radii come from literature constants (20 endothelial nuclei per "
             "mm of vessel, PMID 38723085; mean vessel degree 1.5, PMID "
             "37582180) that were not tuned to the 2.1% target. "
             "results/vessel_geometry_fix.md sec. 5 reports 0.99 / 7.14 / "
             "17.79% across the three available GBM sections, i.e. it does NOT "
             "replicate across sections. " + str(PUBLISHED[
                 "vascular_area_fraction_pct"]["source"]),
        conditional=(),
        values={"area_pct": float(af_row["reconstructed"]),
                "ratio_to_measured": float(af_row["ratio"]),
                "band": [float(slo), float(shi)],
                "module_verdict": str(af_row["verdict"]),
                "legacy_binary_area_pct": float(
                    seg["legacy_binary_area_fraction_pct"])},
    ))
    return rows


def _run_breast_vascular_morphometry() -> List[Row]:
    """Reconstructed BREAST vasculature vs measured human invasive-ductal-
    carcinoma morphometry — the breast twin of :func:`_run_vascular_morphometry`.

    Different tissue, different publications, so this duplicates no GBM row.
    It is also the group with the largest number of honest NOT TESTABLE rows:
    ``experiments/vascular_validation_breast.py`` searched for and did not find
    a breast counterpart to GBM's whole-slide MVD anchor (PMID 27583442) or its
    directly-measured intercapillary-distance anchor (PMID 3359450), and says so.
    """
    from . import vascular_validation as vv

    vvb = _import_experiment("vascular_validation_breast")
    if not os.path.isdir(_repo_path(vvb.BREAST_SAMPLE)):
        raise NotTestable(f"requires the real breast Xenium cache at "
                          f"{vvb.BREAST_SAMPLE} (absent)")
    res = vv.measure(vvb.reconstruct_breast())
    prim, het, kde = res["primary"], res["heterogeneity"], res["kde"]
    rows: List[Row] = []

    # 1. Microvessel density vs Weidner's hotspot counts.
    mvd = float(prim["mvd_per_mm2"])
    lo = vvb.w2mm2(vvb.WEIDNER_NOMET_MEAN)
    hi = vvb.w2mm2(vvb.WEIDNER_MET_MEAN)
    rows.append(Row(
        group="breast-vascular-morphometry", label="microvessel density",
        status=classify(mvd, lo, hi), scope=SCOPE_GEOMETRY,
        quantity="reconstructed vessel-cluster density (per mm^2, 20 um "
                 "single-linkage) vs MEASURED factor-VIII hotspot microvessel "
                 "counts in human invasive ductal breast carcinoma",
        predicted=f"{mvd:.0f} /mm^2  ({mvd / ((lo + hi) / 2):.2f}x band mid)",
        measured=f"{lo:.0f}-{hi:.0f} /mm^2 (45-101 vessels per 200x field, "
                 f"node-negative to node-positive means)",
        pmids=("1701519",),
        bar=f"inside {lo:.0f}-{hi:.0f} /mm^2",
        note="WEAKER EVIDENCE THAN THE GBM EQUIVALENT, twice over. Weidner "
             "counts the single most vascular HOTSPOT, so the published band is "
             "already biased high and exceeding it is worse than it looks. And "
             "the per-field-to-per-mm^2 conversion uses a 0.74 mm^2 field area "
             "INFERRED from a different paper's microscope (PMID 7539869, "
             "bladder), not stated in the breast studies themselves.",
        values={"mvd_per_mm2": mvd, "band": [float(lo), float(hi)],
                "weidner_field_mm2": float(vvb.WEIDNER_FIELD_MM2)},
    ))

    # 2. Vascular area fraction vs Haisan's digital CD34 relative endothelial
    #    area. Same KDE-threshold defect as the GBM row: NOT TESTABLE.
    #    Keyed by threshold, not by position, so a reordering upstream cannot
    #    silently retarget this row at a different threshold.
    kde_abs = {a["threshold"]: a["frac_geom_mask"] for a in kde["absolute"]}
    area_pct = 100.0 * float(kde_abs[0.3])
    rows.append(Row(
        group="breast-vascular-morphometry", label="vascular area fraction",
        status=NOT_TESTABLE, scope=SCOPE_GEOMETRY,
        quantity="reconstructed vascular area fraction (%) vs MEASURED "
                 "whole-slide digital CD34 relative endothelial area in human "
                 "invasive ductal carcinoma (n=50)",
        predicted=f"{area_pct:.2f}% (KDE field threshold 0.3)",
        measured="0.91-4.20% across compartments (stroma 0.91, parenchyma 1.95, "
                 "invasive front 4.20)",
        pmids=("24073397",), bar="n/a",
        reason="THRESHOLD-CHOSEN, NOT A VALIDATION — the same defect as the GBM "
               "row. The value is set by which KDE field threshold is picked and "
               "nothing determines that threshold independently, so landing near "
               "the measured range would be a choice rather than a result. "
               "results/vascular_validation_breast.md sec. 3.5 gives this row no "
               "validation weight. Unlike GBM there is also no threshold-free "
               "breast counterpart registered here.",
        values={"area_pct_at_0.3": area_pct,
                "measured_range_pct": [0.91, 4.20]},
    ))

    # 3. Tile-to-tile heterogeneity. A measured breast number EXISTS but it is
    #    a within-tumour COMPARTMENT contrast, not the between-tile spread the
    #    model produces, so the two are not the same quantity.
    p90p10 = float(het["p90_over_p10"])
    rows.append(Row(
        group="breast-vascular-morphometry", label="vascular heterogeneity",
        status=NOT_TESTABLE, scope=SCOPE_GEOMETRY,
        quantity="p90/p10 endothelial-density spread across 500 um tiles vs "
                 "MEASURED CD34 area contrast between pathologist-defined "
                 "tumour compartments",
        predicted=f"{p90p10:.1f}x p90/p10 across 500 um tiles",
        measured="4.61x invasive front vs tumour-associated stroma (a "
                 "three-compartment contrast, n=50)",
        pmids=("24073397",), bar="n/a",
        reason="NOT THE SAME QUANTITY. The measured 4.61x is a contrast between "
               "three pathologist-defined COMPARTMENTS within a tumour; the "
               "model's number is an unlabelled tile-to-tile spread. No "
               "between-case breast CD34 area-fraction range analogous to GBM's "
               "Sipos 2024 (PMID 38928515) was found despite a dedicated search, "
               "so there is nothing to score the spread against.",
        values={"p90_over_p10": p90p10, "measured_compartment_contrast": 4.61},
    ))

    # 4 & 5. The two anchors GBM has and breast provably does not.
    rows.append(Row(
        group="breast-vascular-morphometry",
        label="intercapillary distance", status=NOT_TESTABLE,
        scope=SCOPE_GEOMETRY,
        quantity="median reconstructed vessel-to-nearest-vessel distance (um) "
                 "vs a measured intercapillary distance in breast IDC",
        predicted=f"{float(prim['vessel_nn_median']):.1f} um",
        measured="NONE",
        pmids=(), bar="n/a",
        reason="NO MEASURED TARGET EXISTS. A dedicated search found no Voronoi "
               "or otherwise directly-measured intercapillary-distance study for "
               "breast IDC. GBM has one (Yoshii & Sugiyama, PMID 3359450); "
               "breast does not, and reporting the GBM value here would be "
               "importing a brain measurement into breast tissue.",
        values={"vessel_nn_median_um": float(prim["vessel_nn_median"])},
    ))
    rows.append(Row(
        group="breast-vascular-morphometry",
        label="tumour cells >100 um from a vessel", status=NOT_TESTABLE,
        scope=SCOPE_GEOMETRY,
        quantity="fraction of reconstructed tumour cells further than 100 um "
                 "from the nearest vessel vs a measured breast equivalent",
        predicted=f"{100.0 * float(prim['tumor_to_vessel_frac_gt_100']):.1f}%",
        measured="NONE",
        pmids=(), bar="n/a",
        reason="NO MEASURED TARGET EXISTS. The GBM row derives its band from "
               "Marien's whole-slide unbiased MVD (PMID 27583442), whose full "
               "text was retrieved and confirmed to cover colorectal, "
               "glioblastoma, ovarian and renal carcinoma and explicitly NOT "
               "breast — the authors themselves flag breast as undone work. No "
               "follow-up exists.",
        values={"frac_gt_100_pct":
                100.0 * float(prim["tumor_to_vessel_frac_gt_100"])},
    ))
    return rows


def _run_endothelial_cellularity() -> List[Row]:
    """Reconstructed section CELLULARITY and endothelial fraction vs measured
    human histology — a different axis from the vessel-geometry anchors above.

    Only section s1 is run (``quick=True`` inside the owning module), because
    the three sections cost minutes and s1 is the section every other anchor
    here uses.

    Of the module's comparisons, exactly one is modality-matched enough to
    score: nuclei per mm^2 measured on 5 um H&E against a 5 um Xenium section.
    The endothelial-fraction comparison is reported NOT TESTABLE, which is the
    module's own conclusion — its entire point is that the apparent over-call
    factorises into definitional terms that carry no endothelial information.
    """
    from . import endothelial_density as ED

    if not os.path.isdir(_repo_path("data_cache", "xenium_gbm")):
        raise NotTestable("requires the real Xenium GBM cache at "
                          "data_cache/xenium_gbm (absent)")
    res = ED.run(sections=("s1",), do_regions=False)
    sec = res["sections"]["s1"]
    if sec.get("status") != "ok":
        raise NotTestable(f"section s1 not evaluated: {sec.get('status')}")
    agg = sec["by_aggregation"][res["production_aggregation"]]
    rec = agg["reconciliation"]

    obs_cells = float(rec["obs_cells_per_mm2"])
    _mean, clo, chi = ED.GBM_CELLULARITY["hgg_nuclei_per_mm2"]
    ctrl = ED.GBM_CELLULARITY["control_nuclei_per_mm2"][0]
    rows = [Row(
        group="endothelial-cellularity", label="section cellularity",
        status=classify(obs_cells, clo, chi), scope=SCOPE_GEOMETRY,
        quantity="reconstructed cells per mm^2 of tissue on a 5 um Xenium "
                 "section vs MEASURED nuclei per mm^2 on 5 um H&E sections of "
                 "human high-grade glioma",
        predicted=f"{obs_cells:.0f} cells/mm^2  "
                  f"({obs_cells / float(_mean):.2f}x the measured mean)",
        measured=f"{float(_mean):.0f} /mm^2 (range {clo:.0f}-{chi:.0f}; "
                 f"matched non-tumour control {ctrl:.0f} /mm^2)",
        pmids=(str(ED.GBM_CELLULARITY["pmid"]),),
        bar=f"inside the measured {clo:.0f}-{chi:.0f} /mm^2 range",
        note="MODALITY-MATCHED, which is why it is scored: both sides are "
             "2-D nuclear profiles on a 5 um section, so no volume-to-section "
             "depth conversion enters. Corroborated to within 5% by an "
             "independent 37-patient, 561-biopsy series "
             f"(4557 +/- 3169 /mm^2, PMID "
             f"{ED.GBM_CELLULARITY_CORROBORATION['pmid']}). Two caveats, both "
             "of which push the SAME way as the observed miss and neither of "
             "which excuses it: the measured fields were selected for MAXIMUM "
             "cellularity in both arms, and Xenium retains only DAPI nuclei "
             "that segment and receive assigned transcripts, so in-situ "
             "under-detection relative to H&E nuclear counting is expected. "
             "The row is scored anyway because the section geometry genuinely "
             "matches; what it bounds is the reconstruction's cell recall, "
             "not the transport model.",
        values={"obs_cells_per_mm2": obs_cells,
                "measured_mean": float(_mean),
                "band": [float(clo), float(chi)],
                "cellularity_ratio_vs_h01": float(rec["cellularity_ratio"])},
    ), Row(
        group="endothelial-cellularity", label="endothelial cell fraction",
        status=NOT_TESTABLE, scope=SCOPE_GEOMETRY,
        quantity="reconstructed endothelial percentage of all cells vs the "
                 "MEASURED endothelial fraction of a 1.05 mm^3 volume-EM "
                 "reconstruction of normal-appearing human temporal cortex",
        predicted=f"{float(rec['obs_ec_per_mm2']):.0f} endothelial/mm^2; "
                  f"fraction ratio {float(rec['fraction_ratio']):.2f}x the "
                  f"reference, of which {1.0 / float(agg['purity']['purity']):.2f}x "
                  f"is estimated classifier over-call",
        measured=f"{float(rec['h01']['ec_fraction_pct']):.2f}% of cells "
                 f"(4,604 endothelial of 57,180, n=1 volume)",
        pmids=(str(ED.H01["pmid"]), str(ED.CORTEX_STEREOLOGY["pmid"])),
        bar="n/a",
        reason="THE TWO SIDES DO NOT COUNT THE SAME THING, and the reference is "
               "the wrong tissue. One unit of the measurement is a whole cell "
               "traced in 3-D through 1.05 mm^3; one unit of the model is a 2-D "
               "polygon on a ~5 um section, so the comparison needs the "
               "profile-counting relation N_A = N_V*(t+h) and an effective "
               "sampling depth that is estimated, not measured. The reference "
               "is also NORMAL temporal cortex, n=1, with no published spread, "
               "while the section is tumour. endothelial_density.py deliberately "
               "reports a factorised BUDGET here rather than a verdict, because "
               f"the raw ratio ({float(rec['apparent_ratio_naive']):.1f}x) is "
               f"dominated by definitional depth and cellularity terms that "
               "carry no endothelial information at all.",
        values={"obs_ec_per_mm2": float(rec["obs_ec_per_mm2"]),
                "h01_ec_fraction_pct": float(rec["h01"]["ec_fraction_pct"]),
                "fraction_ratio": float(rec["fraction_ratio"]),
                "apparent_ratio_naive": float(rec["apparent_ratio_naive"]),
                "apparent_ratio_depth_corrected": float(
                    rec["apparent_ratio_depth_corrected"]),
                "classifier_purity": float(agg["purity"]["purity"])},
    )]
    return rows


def _run_heldout_patient_geometry() -> List[Row]:
    """The SAME published human-GBM anchors, on the section nothing was fitted to.

    Every other GBM row in this registry is scored on ``data_cache/xenium_gbm``
    -- the section the marker gate, the 30 um tissue-mask radius, the 20 um
    vessel linking distance and the 7.5 um vessel radius were all chosen while
    looking at. This group scores the identical published values on
    ``data_cache/xenium_gbm_io``: a DIFFERENT patient on a DIFFERENT 480-gene
    panel (the two panels share 60 genes), reconstructed by the same frozen
    production path with no constant refitted. It is the only out-of-sample
    TISSUE in the registry, and it is structurally the same move as
    ``breast-vascular-morphometry`` -- same anchors, different specimen.

    Only the GEOMETRY half of ``heldout_patient.py`` is registered, and that is
    a deliberate, load-bearing restriction. The module's other half runs ten
    transport-PDE solves, but those solves are scored against **no measured
    number at all** -- they compare the two sections with each other. Nothing in
    this group invokes the transport PDE, and the group must not be cited as
    evidence about the solver. See :data:`OUT_OF_SCOPE` for the parts left out.

    Runtime is ~3-5 min and is dominated by loading and reconstructing 816,769
    cells, NOT by any solve. The whole module was previously left unregistered
    on a runtime estimate that belonged to the drug panel.
    """
    from ..reconstruct.vessel_network import (
        LINK_DISTANCE_UM, VESSEL_RADIUS_UM, segment_vessels, segmentation_report)
    from . import heldout_patient as HP

    spec = next(s for s in HP.SECTIONS if not s["dev"])
    sample = str(spec["sample"])
    if not os.path.isdir(_repo_path(sample)):
        raise NotTestable(f"requires the held-out Xenium GBM cache at {sample} "
                          "(absent)")
    micro = HP.reconstruct_section(sample, HP.DX_UM)
    vessels = segment_vessels(micro.cells, micro.domain,
                              link_distance_um=LINK_DISTANCE_UM)
    seg = segmentation_report(vessels, micro.domain, radius_um=VESSEL_RADIUS_UM)
    by_key = {r["key"]: r for r in seg["rows"]}
    comp = HP.composition_block(micro, HP._retry(
        HP.load_singlecell_spatial, path=sample, dense=False))

    #: Both scored rows use segmentation_report's OWN frozen PUBLISHED bands and
    #: PMIDs, exactly as the section-1 threshold-free area-fraction row above
    #: does, so the bar is the owning module's and not this one's.
    held_out = ("HELD OUT: this section was never inspected while any "
                "reconstruction constant was chosen. The 20 um linking distance "
                "and the vessel radius model are section-1 / literature values, "
                "used here unchanged and NOT refitted.")
    cohort = ("The published value is a population cohort mean, which is "
              "context for a single patient rather than a target; "
              "vessel_network.segmentation_report says so on every row it "
              "emits, and that caveat travels with these two.")

    rows: List[Row] = []

    af = by_key["vascular_area_fraction_pct"]
    alo, ahi = af["published_range"]
    rows.append(Row(
        group="heldout-patient-geometry",
        label="vascular area fraction (held-out section)",
        status=classify(float(af["reconstructed"]), alo, ahi),
        scope=SCOPE_GEOMETRY,
        quantity="reconstructed threshold-free sub-grid vascular area fraction "
                 "(%) on the HELD-OUT section vs MEASURED whole-section CD34+ "
                 "area fraction in IDH-wildtype human GBM",
        predicted=f"{float(af['reconstructed']):.2f}%  "
                  f"({float(af['ratio']):.2f}x the measured mean)",
        measured=f"{float(af['published']):.2f}% (range {alo:.1f}-{ahi:.1f}%, "
                 f"n=30)",
        pmids=(str(af["pmid"]),),
        bar=f"inside vessel_network.PUBLISHED range {alo:.1f}-{ahi:.1f}%",
        note=(held_out + " THE LOAD-BEARING ROW: it is the only anchor here "
              "that tests the vascular FIELD rather than the point pattern, and "
              "it has no threshold to pick. Section 1 lands at 1.25% (0.59x, "
              "LOW); this section lands high, so the two sections miss the same "
              "published value in OPPOSITE directions -- which is why "
              "results/sweep/vascular_area_fraction_generalization_"
              "2026-07-30T03-18-42Z.md could find no section-invariant radius "
              "model that fixes both. " + cohort),
        values={"area_fraction_pct": float(af["reconstructed"]),
                "band": [float(alo), float(ahi)],
                "ratio_vs_published": float(af["ratio"]),
                "section": sample},
    ))

    mv = by_key["mvd_per_mm2"]
    mlo, mhi = mv["published_range"]
    rows.append(Row(
        group="heldout-patient-geometry",
        label="microvessel density (held-out section)",
        status=classify(float(mv["reconstructed"]), mlo, mhi),
        scope=SCOPE_GEOMETRY,
        quantity="reconstructed vessel-cluster density (per mm^2, 20 um "
                 "single-linkage) on the HELD-OUT section vs MEASURED "
                 "whole-slide CD31 microvessel density in human GBM",
        predicted=f"{float(mv['reconstructed']):.0f} /mm^2  "
                  f"({float(mv['ratio']):.2f}x the measured mean)",
        measured=f"{float(mv['published']):.0f} /mm^2 (range {mlo:.0f}-"
                 f"{mhi:.0f}, n=20 patients)",
        pmids=(str(mv["pmid"]),),
        bar=f"inside vessel_network.PUBLISHED range {mlo:.0f}-{mhi:.0f} /mm^2",
        note=(held_out + " This is NOT a held-out failure on its own: section 1 "
              "misses the same band in the same direction (111 /mm^2), so the "
              "two sections AGREE with each other while both disagree with the "
              "literature. That points at the segmentation over-splitting "
              "endothelial cells into vessel objects, not at the held-out "
              "patient. " + cohort),
        values={"mvd_per_mm2": float(mv["reconstructed"]),
                "band": [float(mlo), float(mhi)],
                "ratio_vs_published": float(mv["ratio"]),
                "n_vessels": int(vessels["n_vessels"]),
                "n_endothelial_cells": int(vessels["n_endothelial_cells"]),
                "section": sample},
    ))

    rows.append(Row(
        group="heldout-patient-geometry",
        label="endothelial cell fraction (held-out section)",
        status=NOT_TESTABLE, scope=SCOPE_GEOMETRY,
        quantity="reconstructed endothelial percentage of all cells on the "
                 "HELD-OUT section vs the MEASURED endothelial fraction of a "
                 "1.05 mm^3 volume-EM reconstruction of normal-appearing human "
                 "temporal cortex",
        predicted=f"{float(comp['endothelial_pct']):.2f}% of cells "
                  f"({float(comp['endothelial_fold_vs_volume_em']):.2f}x the "
                  f"reference)",
        measured="8.05% of cells (4,604 endothelial of 57,180, n=1 volume)",
        pmids=("38723085",), bar="n/a",
        reason="NOT TESTABLE FOR TWO INDEPENDENT REASONS, and heldout_patient.py "
               "scores it PASS against a 1.5x tolerance anyway -- which is why "
               "it is carried here explicitly rather than left to that report. "
               "(1) SAME DEFECT AS THE SECTION-1 ROW: endothelial-cellularity "
               "above already reports this exact comparison NOT TESTABLE, "
               "because one unit of the measurement is a whole cell traced in "
               "3-D through 1.05 mm^3 and one unit of the model is a 2-D polygon "
               "on a ~5 um section, and the reference is normal cortex, n=1. "
               "Nothing about a second section repairs that. (2) IT IS NOT "
               "BLIND: results/classifier_fix.md sec. 0 selected the production "
               "typing rule (matched_null_z) on a cross-SECTION consistency "
               "criterion computed WITH this section in it, and that selection "
               "is what moved this number from 4.18% (0.52x -- outside the 1.5x "
               "tolerance) to its present value. The rule was not optimised "
               "against the 8.0% anchor, but it was optimised toward agreement "
               "with a section that is anchored near it, so the resulting "
               "agreement is not out-of-sample evidence.",
        values={"endothelial_pct": float(comp["endothelial_pct"]),
                "fold_vs_volume_em": float(
                    comp["endothelial_fold_vs_volume_em"]),
                "endothelial_per_mm2": float(comp["endothelial_per_mm2"]),
                "cell_density_per_mm2": float(comp["cell_density_per_mm2"]),
                "pre_fix_endothelial_pct": 4.18,
                "section": sample},
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
            runner=_run_vascular_morphometry, est="~1 min",
            needs=("data_cache/xenium_gbm",),
        ),
        AnchorGroup(
            key="breast-vascular-morphometry",
            title="Reconstructed breast vasculature vs measured human IDC "
                  "morphometry",
            tissue="breast", tier=TIER_FULL, scope=SCOPE_GEOMETRY,
            source="experiments/vascular_validation_breast.py "
                   "(drives validation.vascular_validation.measure verbatim)",
            report="results/vascular_validation_breast.md",
            runner=_run_breast_vascular_morphometry, est="~20 s",
            needs=("data_cache/xenium_breast_janesick/rep1",),
        ),
        AnchorGroup(
            key="endothelial-cellularity",
            title="Section cellularity and endothelial fraction vs measured "
                  "human histology",
            tissue="GBM", tier=TIER_FULL, scope=SCOPE_GEOMETRY,
            source="tmesim/validation/endothelial_density.py::run "
                   "(section s1 only) -> definitional_reconciliation",
            report="results/endothelial_density.md",
            runner=_run_endothelial_cellularity, est="~1-2 min (section s1 only)",
            needs=("data_cache/xenium_gbm",),
        ),
        AnchorGroup(
            key="heldout-patient-geometry",
            title="HELD-OUT second GBM patient: reconstructed geometry vs the "
                  "same measured human GBM histology",
            tissue="GBM", tier=TIER_FULL, scope=SCOPE_GEOMETRY,
            source="tmesim/validation/heldout_patient.py::reconstruct_section + "
                   "vascular_block / composition_block, section 2 "
                   "(data_cache/xenium_gbm_io) -- GEOMETRY ONLY, the module's "
                   "10 PDE solves are NOT registered",
            report="results/heldout_patient.md, "
                   "results/heldout_patient_validation.md",
            runner=_run_heldout_patient_geometry,
            est="~2-3 min (measured 126 s on a loaded machine: load + "
                "reconstruct 816,769 cells + vessel segmentation; NO PDE solve)",
            needs=("data_cache/xenium_gbm_io",),
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
    ("tmesim/validation/heldout_patient.py -- the DRUG PANEL and HYPOXIA halves "
     "(its geometry half IS now registered, as heldout-patient-geometry)",
     "SUPERSEDES an earlier note that said this module was omitted on runtime "
     "grounds only. That was wrong twice over. The geometry half costs ~4 min "
     "and no PDE solve, so runtime never justified hiding it, and it is now "
     "registered. The two halves that stay out are out on substance: (a) the "
     "five-drug panel runs 10 transport-PDE solves but is scored against NO "
     "measured number -- it compares the two sections with EACH OTHER, so it "
     "yields no external anchor, and it is downstream of alpha besides; (b) the "
     "VEGFA/HIF1A hypoxia check is model-vs-marker WITHIN one section, the same "
     "category state_validation.py is excluded under, and it is a re-derivation "
     "on the very section that finding was made on, so it is not out-of-sample. "
     "Worth running directly anyway: computing that comparison arithmetically "
     "is what shows results/state_validation.md sec. 1.2 has its inequality "
     "backwards (|-0.1523| > |-0.1204|, so the solved field OUT-predicts the "
     "distance baseline it is said to lose to)."),
    ("results/oxygen_calibration_breast.md / oxygen_calibration_landscape.md",
     "the mmHg <-> dimensionless-O mapping was never derived for these tissues "
     "(REFERENCE_PO2_MMHG = 35.0 is a normal-human-BRAIN anchor), so nothing is "
     "scoreable. Reported in those files as an honest null."),
    ("results/breast_tissue_permeability.md",
     "a falsified FIX to an anchor already registered here "
     "(breast-transvascular), not a separate measurement."),
    ("results/binding_sink_refuted.md (PMID 18787762 human glioma ECF "
     "elimination k30, 3.15-14.51 h half-life)",
     "a genuine model-parameter-vs-measurement falsification, but NO module "
     "computes it: the document is hand-written with no generator in the repo, "
     "so there is nothing to run. It is the diagnosis attached to the "
     "breast-transvascular failure rather than a separate anchor."),
    ("tmesim/validation/pk_reality.py",
     "the pre-calibration predecessor of pk_multidrug.run_penetration on the "
     "SAME drug and the SAME anchor (PMID 16361566); superseded, and "
     "registering it would double-count."),
    ("tmesim/validation/pk_anchors.py",
     "a data compendium and unit converters. It imports nothing from the model "
     "and runs no simulation, so it produces no prediction to score."),
    ("tmesim/validation/revalidation.py",
     "re-runs the glioma / pk / gist / reconstruction entry points under "
     "current code. Every number in it is already registered or already "
     "excluded above."),
    ("tmesim/validation/vessel_network_impact.py",
     "an A/B ablation of the KDE vs discrete vessel source; it reuses "
     "values and PMIDs verbatim from results/vascular_validation.md. Its one "
     "novel number is the threshold-free area fraction, which is registered "
     "here under vascular-morphometry and attributed to its real owner, "
     "reconstruct/vessel_network.py::segmentation_report."),
    ("tmesim/validation/oxygen_benchmark.py",
     "this IS results/oxygen_demotion.md — within-section marker correlation, "
     "already excluded above."),
    ("tmesim/validation/mgmt_spatial.py, classifier_fix.py",
     "spatial autocorrelation against the assay's own negative-control probes, "
     "and a bug fix scored on cross-section CONSISTENCY. classifier_fix's own "
     "report states the volume-EM anchor 'was never optimised against' and is "
     "not its target."),
    ("tmesim/validation/cell_state_impact.py",
     "compares three state-labelling routes to each other. Its Ki-67 and "
     "Thomlinson-Gray numbers are probe-sensitivity checks on the DATA and a "
     "model INPUT respectively, not scored model outputs."),
    ("tmesim/validation/uncertainty.py, ranking_robustness.py, propagation.py, "
     "delivery_metrics.py",
     "Monte-Carlo error bars, order-stability under parameter uncertainty, "
     "flag-propagation instrumentation, and metric self-validation against a "
     "source null. No external measured number is scored in any of them."),
    ("tmesim/validation/experiment.py, baselines.py, metrics.py "
     "(results/stage0_result.json, results/STAGE0_FINDINGS.md)",
     "the Stage-0 kill-test harness. The only run on disk is a synthetic "
     "fixture cohort from the archived synthetic era, not real tissue."),
    ("experiments/cross_tissue.py, real_drug_panel.py, breast_drug_panel.py, "
     "dose_schedule.py, combination.py, niche_response.py, regimen.py, "
     "patient_report.py, flagship_markers.py",
     "forward-prediction drivers on real geometry. Their PMIDs are data "
     "provenance (Janesick PMID 38114474) or plasma-PK INPUTS (Stupp PMID "
     "15758009 and ~20 others in regimen.py); no model output is compared to "
     "any of them."),
    ("results/gbm_3d_registration.md, results/registration_3d_breast_poc.md",
     "registration scored by Chamfer/Hausdorff/density correlation against its "
     "own unregistered baseline. The PMIDs are dataset provenance."),
    ("results/hypoxia_imputation_model.md, results/oxygen_sensitivity_band.md, "
     "results/mechanism_attribution.md, results/vessel_promotion.md",
     "held-out-cell ML imputation within one platform; a parameter sweep plus "
     "two-solver numerical agreement; a read-out of the kill model's own "
     "arithmetic; and a default-promotion decision that re-runs an already "
     "registered check. None introduces a new measurement."),
    ("results/tumour_physiology_pbpk_literature.md",
     "explicitly 'reference material, not wired into any model' — a survey."),
)


# --------------------------------------------------------------------------- #
# Suite driver
# --------------------------------------------------------------------------- #
def _missing(needs: Sequence[str]) -> List[str]:
    return [p for p in needs if not os.path.exists(_repo_path(p))]


def _stamp_physics(rows: Sequence[Row]) -> List[Row]:
    """Give every row a physics-provenance string, and refuse to score a
    solver number that did not declare one.

    A row whose ``scope`` is ``PDE`` and that carries a PASS/FAIL verdict is a
    claim about the transport model, and the breast anchor shows those claims
    move by ~10x depending on which delivery law was active. A row like that
    with no declared law is therefore not a weaker result, it is an
    UNINTERPRETABLE one, so it is downgraded to NOT TESTABLE rather than left
    to pass quietly. Rows that never invoked the solver — analytic, geometry,
    and any PDE-scoped row that produced no number — are stamped
    :data:`PHYSICS_NOT_APPLICABLE`, which is a statement that the question does
    not arise, not that the physics agreed.
    """
    out: List[Row] = []
    for r in rows:
        if r.physics:
            out.append(r)
            continue
        produced_a_number = r.status in (PASS, FAIL)
        if r.scope == SCOPE_PDE and produced_a_number:
            r.status = NOT_TESTABLE
            r.reason = (
                "PHYSICS PROVENANCE MISSING. This row ran the transport solver "
                "and produced a verdict, but did not declare which "
                "transvascular delivery law it ran under. The breast "
                "interstitial:plasma anchor moves by roughly 10x between "
                "transvascular_exchange=False and True, so an undeclared "
                "number cannot be interpreted at all — it is not a weaker "
                "result, it is an unlabelled one. Refusing to score it. "
                + (r.reason or "")
            ).strip()
            r.bar = ""
            r.physics = "UNDECLARED — this is the defect described above"
        else:
            r.physics = PHYSICS_NOT_APPLICABLE
        out.append(r)
    return out


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
        rows.extend(_stamp_physics(grows))
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
    default_exchange = repo_default_exchange()
    return {
        "mode": "quick" if quick else "full",
        "exchange": bool(exchange),
        # Physics provenance for the whole run, so a stored JSON scorecard can
        # never be compared against a later one without the delivery law being
        # visible on both sides.
        "physics": {
            "transvascular_exchange_requested": bool(exchange),
            "transvascular_exchange_repo_default": default_exchange,
            "is_shipping_default": bool(exchange) == default_exchange,
            "law": PHYSICS_LAWS[bool(exchange)],
            "note": "The breast interstitial:plasma anchor differs by roughly "
                    "10x between the two laws. Any comparison of two "
                    "scorecards that does not first check this field is "
                    "meaningless.",
        },
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
