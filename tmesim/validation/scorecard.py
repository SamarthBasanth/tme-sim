"""THE EXTERNAL-ANCHOR SCORECARD — one runnable inventory of every place this
model is scored against an INDEPENDENT MEASURED number.

This module is the data layer behind ``codon validate`` (``tmesim.cli``). It
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
    softened.** ``codon validate`` exits non-zero when any row is FAIL, so this
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
    "ifp-uncalibrated": (
        "Every coefficient in tmesim/physics/pressure.py is a hardcoded "
        "order-of-magnitude constant, by its own docstring 'not "
        "patient-calibrated': _K0 = 1.2, _K_ECM_ATTEN = 0.6, _PV = 25.0 mmHg, "
        "_LP_NORMAL = 1e-3, _LEAK_GAIN = 12.0, _LAMBDA_LYMPH = 3e-2. None was "
        "fitted to any measurement, which makes any IFP score genuinely "
        "out-of-sample — and also means nothing constrains the absolute scale. "
        "_PV = 25.0 in particular is an upper bound on the predicted core "
        "plateau and was NOT chosen against any of the anchors scored here."
    ),
    "ifp-unresolved-mesh": (
        "The screened-Poisson screening length sqrt(K/(a+d)) implied by those "
        "coefficients is ~6-7 um on both real sections, BELOW every grid "
        "spacing the engine runs at (8-30 um), at >99.7% of tissue nodes. The "
        "reaction term dominates the diagonal, so the discrete solution "
        "collapses toward the pointwise algebraic limit a*p_v/(a+d) (measured "
        "correlation r = 0.90 GBM / 0.93 breast) and the interior Laplacian "
        "barely couples neighbours. Any IFP number is conditional on a mesh "
        "that does not resolve the operator being solved "
        "(results/ifp_validation.md sec. 5)."
    ),
    "ifp-grid-unconverged": (
        "The wick-scale IFP statistic is not grid-converged. Across dx = 24, "
        "18, 12, 8 um it moves 0.358 -> 0.434 -> 0.462 -> 0.370 mmHg on GBM "
        "(non-monotone, 1.29x spread) and 0.621 -> 0.424 -> 0.370 -> 0.286 "
        "mmHg on breast (monotone, 2.17x spread, still falling 22.5% over the "
        "last refinement with no plateau). Two contributions are known: the "
        "vascular-density KDE is rebuilt per grid, and _TUMOR_SIGMA = 5.0 in "
        "pressure.py is specified in GRID CELLS, so the physical tumour "
        "smoothing scale is 40/60/90/120 um at dx = 8/12/18/24 um and the "
        "model itself changes with the mesh."
    ),
    "mcl-no-vasculature": (
        "COMMENSURABILITY BOUND, and it is the load-bearing caveat on this "
        "group. A multicellular layer is a stack of cells grown on a porous "
        "support: packed cells, little stroma, NO VASCULATURE, no necrosis, no "
        "interstitial convection and no plasma source term. Hicks measured flux "
        "across that slab and fitted it to a diffusion-reaction equation of the "
        "SAME FORM the engine solves, which is exactly what makes D and k "
        "commensurable AS COEFFICIENTS. What that licenses: scoring the "
        "engine's compiled per-drug D_eff, and the hindrance law that produces "
        "it, against a directly measured tissue diffusion coefficient in the "
        "same units and the same functional role. What it does NOT license: any "
        "claim about a vascularised tissue section. Nothing in an MCL "
        "constrains vessel spacing, the transvascular source term, IFP-driven "
        "convection, or the spatial concentration field the solver produces on "
        "real reconstructed geometry. Reading an MCL row as evidence that the "
        "engine's tissue predictions are validated would be the modality "
        "pooling results/external_anchor_survey.md sec. 0 catalogues three live "
        "instances of, and which this project has manufactured a false "
        "regression from twice."
    ),
    "cell-exclusion-default-off": (
        "The cellular volume-exclusion term "
        "(molecular/ecm_mesh.py::cell_exclusion, Maxwell's "
        "impermeable-inclusion limit) is MERGED BUT NOT WIRED IN. "
        "results/cellular_exclusion.md sec. 7 recommends merging the physics and "
        "explicitly NOT flipping the default; nothing in compile_drug calls it, "
        "so DrugParams.diffusivity is bit-identical to the pre-merge value. "
        "Every VERDICT in this group is therefore scored on the engine AS IT "
        "SHIPS, with the term OFF. The cell-term-ON arm is reported alongside "
        "each verdict and carries three caveats of its own, none of which is a "
        "footnote: (a) phi_cell = 0.7892 is measured by Delesse's principle on "
        "the repo's real Xenium BREAST TUMOUR section and then applied to a V79 "
        "hamster multicellular layer — a literature search for an MCL "
        "extracellular-volume measurement returned nothing usable, so the "
        "packing density is imported from the wrong system "
        "(cellular_exclusion.md sec. 4.1); (b) Xenium's boundary segmentation "
        "tiles the section, so phi_cell is an UPPER estimate and every "
        "hindrance derived from it is too; (c) stacking the cell term on the "
        "existing chain DOUBLE-COUNTS by up to ~1.9x, because brownian.py's "
        "300 nm obstacle field and ALPHA_TORTUOSITY = 0.6 were partly standing "
        "in for cell-scale obstruction, and cellular_exclusion.md sec. 5.3 "
        "states in terms that it cannot resolve which."
    ),
    "coefficient-not-field": (
        "This row scores an INPUT COEFFICIENT of the transport PDE, not any "
        "solved field. Two provenance routes carry this tag and both were "
        "verified rather than asserted, by a call counter on "
        "physics.diffusion.simulate_drug_transport and "
        "sim.engine.simulate_response: (i) the MCL coefficient rows, whose "
        "prediction is read straight off MolecularProfile after "
        "descriptors.from_smiles -> compile_drug.build_profile (counters 0/0, "
        "plus a sys.setprofile sweep, recorded in "
        "results/anchor_registration.md sec. 3); (ii) the brain-slice "
        "tissue-capacity row, whose prediction is the two-route capacity term "
        "s_eff evaluated at the shipping constants (counters 0/0, recorded in "
        "results/brain_slice_anchor.md sec. 0). It is a stronger row than an "
        "unanchored coefficient and a weaker one than a field test: it can "
        "FALSIFY the coefficient law it scores, and it cannot confirm or "
        "refute the solver."
    ),
    "capacity-normalised-to-one": (
        "THE ENGINE'S TISSUE CAPACITY IS 1.0 BY CONSTRUCTION, NOT BY CHOICE, "
        "and that is the load-bearing caveat on the brain-slice row. The "
        "quantity scored is the two-route capacity term written down in "
        "molecular/two_route.py's own derivation, s_eff = (1 - phi) + "
        "phi*K*(1 + K_b): matrix (interstitial) capacity normalised to 1, "
        "inclusion (cell) capacity K*(1 + K_b), combined by an exact volume "
        "average. At the shipping constants phi = 0.7892, K = cell_partition "
        "= 1.0 and K_b absent (compile_drug.build_profile never passes "
        "bound_to_free_ratio), so s_eff = 1.000000 EXACTLY, for every "
        "compound. Nobody set that value: it is the normalisation implied by "
        "publishing a CONDUCTANCE and declining to publish a capacity "
        "(two_route.py: 'this module's published output is a CONDUCTANCE, "
        "computed against a convention in which the no-binding capacity is "
        "1'). EVERY dimensionless partition/capacity quantity in the shipping "
        "engine equals exactly 1.0 -- cell_partition, "
        "two_route.retardation_factor, the unshipped ion_trapping "
        "retardation R, and s_eff -- so the verdict does not depend on which "
        "of them a reader prefers. TWO CONSEQUENCES FOR THE FIX, and neither "
        "knob is clean. (1) K (cell_partition) is NOT a capacity-only knob: "
        "it also sets the inclusion conductivity in two_route_factor "
        "(sigma_c = cell_partition * cytoplasm_D), so raising it to lift the "
        "capacity would move D_eff and re-open the Pruijn HT-29 fit -- and "
        "prereg sec. 3.2 independently forbids using Vu,brain, a LUMPED "
        "quantity, to set K, which is one component of the lump. THIS ROW "
        "DOES NOT LICENSE MOVING cell_partition AND MUST NOT BE CITED AS IF "
        "IT DID; K remains unconstrained by any measurement. (2) K_b "
        "(intracellular retardation) IS the per-compound capacity term the "
        "engine lacks -- and it was built and REFUTED on Pruijn HT-29 "
        "(results/intracellular_retardation.md: R^2 0.718 -> 0.606; "
        "matched-MW rank agreement 77.4% -> 51.6%). The term this failure "
        "points at and the term the earlier work rejected are the SAME TERM, "
        "and the two datasets point in opposite directions. That tension is "
        "unresolved and is a model decision, not a scorecard one."
    ),
    "brain-slice-not-tumour": (
        "COMMENSURABILITY BOUND (prereg results/transport_anchor_prereg.md "
        "sec. 3.10). The measurement is an unbound distribution volume in RAT "
        "BRAIN SLICES: 300 um fresh slices in oxygenated artificial ECF, "
        "Vi-corrected, at the assay's own pH -- no vasculature, no plasma "
        "source term, no tumour, and not human. What it licenses: scoring the "
        "engine's implied tissue drug CAPACITY, an input coefficient, in the "
        "same units and the same functional role. What it does NOT license: "
        "any statement about the transport solver (never entered here), any "
        "solved concentration field, vessel spacing, transvascular delivery, "
        "interstitial convection, any absolute kill fraction, D_eff accuracy, "
        "or human tumour tissue. The anchor is weakest exactly where the "
        "quantity matters most: tumour interstitium is acidic and rat brain "
        "is not, and the ion-trapping direction that follows is precisely the "
        "term the engine does not ship."
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


#: Cell volume fraction of the repo's own real Xenium breast section, derived by
#: Delesse's principle (areal fraction = volume fraction on a random plane) from
#: 167,780 segmented cells: sum(cell_area) 3.249e+07 um^2 over a 4.118e+07 um^2
#: bounding box. Published in results/cellular_exclusion.md sec. 3.2 and produced
#: by experiments/cellular_exclusion.py --stage phi. It is NOT a code default —
#: nothing in tmesim/ reads it — and it is used here ONLY to report the
#: cell-term-ON arm alongside each verdict. NO VERDICT IN THIS GROUP DEPENDS ON
#: IT: every status below is computed from the shipping compiled coefficient.
PHI_CELL_BREAST_SECTION = 0.7892


def _run_mcl_transport_coefficients() -> List[Row]:
    """The engine's two compiled transport coefficients vs Hicks' MEASURED
    multicellular-layer diffusion coefficient, hindrance factor and consumption
    rate.

    Why this group is worth a registry slot
    ---------------------------------------
    The transport PDE has exactly two per-drug coefficients,
    ``dC/dt = D_eff*lap(C) - k*C + S``, and Hicks fitted flux through a
    multicellular layer **to a diffusion-reaction model of the same form** to
    obtain both. Same units, same functional role, same equation: no unit
    conversion, no compartment translation, no modality bridge to get wrong.
    Tirapazamine is not on the engine's panel and nothing in the repo was ever
    fitted to a TPZ measurement, so all three rows are fully out-of-sample.

    Why the scope is ``analytic`` and not ``PDE`` — the crux
    -------------------------------------------------------
    These rows score the solver's INPUTS, not its OUTPUTS. The prediction comes
    from ``descriptors.from_smiles`` -> ``compile_drug.build_profile`` and is
    read off ``MolecularProfile``; no solve happens. Registering them as PDE
    evidence would move the summary's "N of N scored rows never invoke the
    transport PDE" line for free and by mislabelling, which is worse than
    leaving that line where it is. They are ``analytic`` for the same reason
    ``glioma_regional_pk`` is, and they carry ``coefficient-not-field`` to say
    so in the row rather than in a footnote.

    Three separate rows, because they can pass and fail independently
    ----------------------------------------------------------------
    The diffusivity is absolute in um^2/s, the hindrance is a dimensionless
    ratio against the same drug's free diffusivity, and the consumption rate is
    an absolute 1/s. The engine can be wrong about any one without being wrong
    about the others — and in fact is: it misses the absolute band and the
    dimensionless ratio by different factors and in a way that
    results/cellular_exclusion.md sec. 4.2 reads as a residual ``D_free`` error.

    Scored ON THE ENGINE AS IT SHIPS
    --------------------------------
    The verdicts below are the shipping compiled numbers, and since
    ``results/two_route_promotion.md`` the shipping compile INCLUDES the cell
    term: ``compile_drug.build_profile`` applies the two-route (paracellular +
    transcellular) law by default. That promotion moved a verdict in this group,
    which is exactly what a scorecard exists to notice:

    * the dimensionless **hindrance** row goes **FAIL -> PASS** — 1.77-fold
      shipping before, **10.33-fold** now, against a measured 12-fold;
    * the absolute **effective diffusivity** row stays **FAIL** but FLIPS
      DIRECTION — 287.3 um^2/s (2.21x above the measured band) before,
      **49.3 um^2/s** (1.50x below it) now;
    * the **consumption-rate** row stays NOT TESTABLE, but the reading it
      declines to score changed sign: at the old alpha, scoring it would have
      manufactured a PASS; at the new alpha it would manufacture a FAIL. That
      sentence is corrected in the row rather than left to rot, because a
      "declining to score removes a pass" claim that has become "removes a
      failure" is the most damaging kind of stale disclosure.

    The pre-promotion arm is computed live with ``two_route=False`` and reported
    in every row's ``values`` under ``superseded_``, never as the verdict.

    NOTHING IS TUNED HERE. The measured constants are imported from
    ``experiments/anchor_prototype.py``, which owns them and traces each to its
    PMID; the predictions come from the unmodified molecular pipeline at the
    repo's own ``ecm_reference = 0.4``.
    """
    from ..molecular import compile_drug, descriptors, ecm_mesh
    from .pk_reality import ORDER_OF_MAGNITUDE_FACTOR

    AP = _import_experiment("anchor_prototype")

    mol = descriptors.from_smiles("tirapazamine", AP.TPZ_SMILES)
    prof = compile_drug.build_profile(mol, ecm_reference=AP.ECM_REFERENCE)
    # The PRE-PROMOTION arm, computed live rather than quoted from a document,
    # so the before/after in every note below cannot drift from the code.
    prof_old = compile_drug.build_profile(mol, ecm_reference=AP.ECM_REFERENCE,
                                          two_route=False)

    d_ship = float(prof.effective_diffusivity)
    d_free = float(prof.free_diffusivity)
    fold_ship = 1.0 / float(prof.hindrance)
    d_super = float(prof_old.effective_diffusivity)
    fold_super = 1.0 / float(prof_old.hindrance)

    # The IMPERMEABLE-cell term on its own, computed live rather than quoted.
    cell = ecm_mesh.cell_exclusion(PHI_CELL_BREAST_SECTION)
    cell_scale = float(cell["hindrance"])
    # Arm B ("cell term ALONE") is the MCL-appropriate arm: an MCL is packed
    # cells with little matrix, and the repo's ECM legs describe matrix the MCL
    # does not have. Arm A is the whole existing chain with the cell term added,
    # which since the promotion is what the engine SHIPS -- the only difference
    # from d_ship is that the shipping term is the PERMEABLE-cell generalisation,
    # and for tirapazamine the transcellular route adds only 2.2 %.
    d_armB, fold_armB = d_free * cell_scale, 1.0 / cell_scale
    d_armA = d_super * cell_scale
    fold_armA = d_free / d_armA

    lo = min(v[0] for v in AP.D_MCL_MEASURED_UM2_S.values())
    hi = max(v[0] for v in AP.D_MCL_MEASURED_UM2_S.values())
    d_pmids = tuple(dict.fromkeys(v[2] for v in AP.D_MCL_MEASURED_UM2_S.values()))

    mcl_arm = (
        f"PRE-PROMOTION arm, reported for comparison and NOT scored: with "
        f"two_route=False the same compile gives D_eff {d_super:.1f} um^2/s and "
        f"{fold_super:.2f}-fold hindrance, which is what every number in "
        f"results/diffusivity_falsification.md and results/cellular_exclusion.md "
        f"was computed against. Two further arms, also unscored: arm B (the "
        f"IMPERMEABLE-cell term alone, the MCL-appropriate arm, since an MCL is "
        f"packed cells with little matrix) gives {d_armB:.1f} um^2/s and "
        f"{fold_armB:.2f}-fold; arm A (the old chain plus the impermeable term) "
        f"gives {d_armA:.1f} um^2/s and {fold_armA:.2f}-fold. Arm A's proximity "
        f"to the measurement is partly the WRONG TERM DOING THE RIGHT WORK — it "
        f"credits an MCL with ECM the MCL does not have "
        f"(results/cellular_exclusion.md sec. 4.1), and that objection now "
        f"applies to the SHIPPING number too, because the shipping chain is arm "
        f"A with a permeable membrane. Arm B's fold is {fold_armB:.4f} exactly; "
        f"cellular_exclusion.md sec. 0 prints it as 6.61x, rounded down.")

    rows: List[Row] = []

    # ---------------------------------------------------------------- row 1 #
    rows.append(Row(
        group="mcl-transport-coefficients",
        label="effective diffusivity (tirapazamine)",
        status=classify(d_ship, lo, hi), scope=SCOPE_ANALYTIC,
        quantity="compiled MolecularProfile.effective_diffusivity (um^2/s) — "
                 "the D_eff the transport PDE would be handed — vs the MEASURED "
                 "tissue diffusion coefficient of the same drug in "
                 "multicellular layers, obtained by fitting slab flux to a "
                 "diffusion-reaction model of the same form",
        predicted=f"{d_ship:.1f} um^2/s as shipped "
                  f"({d_ship / lo:.2f}x the lowest measured cell line, "
                  f"{d_ship / hi:.2f}x the highest); {d_super:.1f} um^2/s before "
                  f"the two-route promotion, i.e. the miss FLIPPED DIRECTION "
                  f"from {d_super / hi:.2f}x ABOVE the band to {d_ship / lo:.2f}x "
                  f"BELOW it",
        measured=f"{lo:.0f}-{hi:.0f} um^2/s across three MCL cell lines "
                 f"(V79-171b 74, MGH-U1 130, PMID 9806526; SiHa 130 +/- 5, "
                 f"PMID 17869669)",
        pmids=d_pmids,
        bar=f"inside the measured {lo:.0f}-{hi:.0f} um^2/s range. "
            f"results/external_anchor_survey.md sec. 4 recommends this anchor "
            f"as a FAIL row and that is the bar adopted",
        note=f"THE MISS FLIPPED DIRECTION AND DID NOT CLOSE. The three "
             f"measurements agree with EACH OTHER to 1.8x. Before the two-route "
             f"promotion the engine was above all three ({d_super:.1f} um^2/s); "
             f"it is now below all three ({d_ship:.1f} um^2/s), by a smaller "
             f"factor. Overshooting a band in the other direction is progress "
             f"in magnitude and is NOT a pass, and this row is not recorded as "
             f"one. It does pass the repo's loosest declared band "
             f"({ORDER_OF_MAGNITUDE_FACTOR:.0f}x, "
             f"pk_reality.ORDER_OF_MAGNITUDE_FACTOR), which is not the bar "
             f"here. READ THIS ROW AGAINST THE HINDRANCE ROW BELOW, WHICH NOW "
             f"PASSES: the dimensionless ratio lands inside its band while the "
             f"absolute value does not, which is the signature of a residual "
             f"FREE-diffusivity error — the hindrance is right and the thing it "
             f"is hindering is not (results/cellular_exclusion.md sec. 4.2). "
             f"{mcl_arm}",
        conditional=("mcl-no-vasculature", "cell-exclusion-default-off",
                     "coefficient-not-field"),
        values={"shipped_D_eff_um2s": d_ship,
                "superseded_mw_only_D_eff_um2s": d_super,
                "free_diffusivity_um2s": d_free,
                "measured_band_um2s": [float(lo), float(hi)],
                "cell_term_on_armA_D_eff_um2s": d_armA,
                "cell_term_on_armB_D_eff_um2s": d_armB,
                "phi_cell": PHI_CELL_BREAST_SECTION,
                "ecm_reference": float(AP.ECM_REFERENCE)},
    ))

    # ---------------------------------------------------------------- row 2 #
    # The measured 12-fold carries NO published uncertainty (the abstract says
    # "12-fold lower than in culture medium" and nothing more), so a band has to
    # come from somewhere. Rather than invent one, use the repo's own loosest
    # declared band — the 3x order-of-magnitude factor pk_reality already ships.
    # That is the most generous bar available, and the shipping engine still
    # misses it.
    h_lo, h_hi = (AP.MEASURED_HINDRANCE_FOLD_V79 / ORDER_OF_MAGNITUDE_FACTOR,
                  AP.MEASURED_HINDRANCE_FOLD_V79 * ORDER_OF_MAGNITUDE_FACTOR)
    rows.append(Row(
        group="mcl-transport-coefficients",
        label="hindrance factor D_free/D_eff (tirapazamine)",
        status=classify(fold_ship, h_lo, h_hi), scope=SCOPE_ANALYTIC,
        quantity="compiled 1/MolecularProfile.hindrance — the dimensionless "
                 "fold reduction of the same drug's diffusivity in tissue "
                 "versus in free medium — vs the MEASURED fold reduction of "
                 "tirapazamine's diffusion coefficient in a V79-171b "
                 "multicellular layer relative to culture medium",
        predicted=f"{fold_ship:.2f}-fold as shipped "
                  f"({fold_ship / AP.MEASURED_HINDRANCE_FOLD_V79:.2f}x "
                  f"measured); {fold_super:.2f}-fold "
                  f"({fold_super / AP.MEASURED_HINDRANCE_FOLD_V79:.2f}x measured) "
                  f"before the two-route promotion",
        measured=f"{AP.MEASURED_HINDRANCE_FOLD_V79:.0f}-fold lower than in "
                 f"culture medium (V79-171b MCL; no uncertainty published)",
        pmids=("9806526",),
        bar=f"within the repo's own loosest declared band, "
            f"pk_reality.ORDER_OF_MAGNITUDE_FACTOR = "
            f"{ORDER_OF_MAGNITUDE_FACTOR:.0f}x either side of the measured "
            f"{AP.MEASURED_HINDRANCE_FOLD_V79:.0f}-fold, i.e. "
            f"{h_lo:.1f}-{h_hi:.1f}-fold. The measured value carries no "
            f"published uncertainty, so inventing a tighter band would be "
            f"inventing a number; this is the most GENEROUS bar available, it "
            f"was chosen while the engine was FAILING it, and it has not been "
            f"widened or moved now that the engine passes",
        note=f"THE SHARPEST ROW IN THIS GROUP, AND THE ONE THE TWO-ROUTE "
             f"PROMOTION MOVED: no units, no compartment translation, no cohort "
             f"matching — a pure dimensionless ratio of one drug's diffusivity "
             f"in tissue versus in free medium. It went "
             f"{fold_super:.2f}x -> {fold_ship:.2f}x against a measured "
             f"{AP.MEASURED_HINDRANCE_FOLD_V79:.0f}x, i.e. "
             f"FAIL -> PASS ({fold_ship / AP.MEASURED_HINDRANCE_FOLD_V79:.2f}x "
             f"measured, inside a band of {h_lo:.0f}-{h_hi:.0f}). "
             f"WHAT THAT DOES AND DOES NOT MEAN. The diagnosis this row carried "
             f"was structural and correct: the old engine's only obstacle "
             f"population was ECM, so cellular volume exclusion — most of a "
             f"tumour's volume — was absent from the hindrance law, and an MCL "
             f"is almost pure packed cells, which is exactly why it exposed "
             f"that. The shipping law now has cells in it and the row passes. "
             f"But the band is 3x wide either side of a measurement that "
             f"carries NO published uncertainty, the pass is "
             f"{AP.MEASURED_HINDRANCE_FOLD_V79 / fold_ship:.2f}x short of the "
             f"central value, and the ABSOLUTE diffusivity row above still "
             f"fails — so what has been established is that the hindrance RATIO "
             f"is now the right size, not that the transport model is right. "
             f"{mcl_arm}",
        conditional=("mcl-no-vasculature", "cell-exclusion-default-off",
                     "coefficient-not-field"),
        values={"shipped_fold": fold_ship,
                "superseded_mw_only_fold": fold_super,
                "measured_fold": float(AP.MEASURED_HINDRANCE_FOLD_V79),
                "band_fold": [float(h_lo), float(h_hi)],
                "cell_term_on_armA_fold": fold_armA,
                "cell_term_on_armB_fold": fold_armB,
                "armB_shortfall_x": float(
                    AP.MEASURED_HINDRANCE_FOLD_V79 / fold_armB),
                "phi_cell": PHI_CELL_BREAST_SECTION},
    ))

    # ---------------------------------------------------------------- row 3 #
    k_raw = float(prof.clearance + prof.tissue_binding_rate)
    from .calibration import calibrate_alpha

    alpha = float(calibrate_alpha())
    k_cal = float(prof.clearance + alpha * prof.tissue_binding_rate)
    k_meas = {line: v[0] for line, v in AP.K_CONSUMPTION_MEASURED_PER_S.items()}
    k_lo, k_hi = min(k_meas.values()), max(k_meas.values())
    rows.append(Row(
        group="mcl-transport-coefficients",
        label="first-order tissue consumption rate k (tirapazamine)",
        status=NOT_TESTABLE, scope=SCOPE_ANALYTIC,
        quantity="the engine's first-order interstitial loss rate k (1/s) — "
                 "the -kC sink of the transport PDE — vs the MEASURED anoxic "
                 "metabolic consumption rate constant of tirapazamine fitted "
                 "in the same position of the same diffusion-reaction equation",
        predicted=f"TWO readings that differ by {k_cal / k_raw:.0f}x and land "
                  f"on OPPOSITE sides of the measurement: k_raw = clearance + "
                  f"binding = {k_raw:.3e} /s ({k_raw / k_hi:.4f}-"
                  f"{k_raw / k_lo:.4f}x measured), and k_calibrated = clearance "
                  f"+ alpha*binding = {k_cal:.3e} /s at alpha = {alpha:.4f} "
                  f"({k_cal / k_hi:.3f}-{k_cal / k_lo:.3f}x measured)",
        measured=f"{k_hi:.3e} /s (0.52 /min, V79-171b) and {k_lo:.3e} /s "
                 f"(0.31 /min, MGH-U1)",
        pmids=("9806526",),
        bar="n/a",
        reason=(
            "PARTIALLY COMMENSURABLE ONLY, and the owning survey says so before "
            "any number is quoted. TWO INDEPENDENT REASONS, either sufficient. "
            "(1) THE PROCESS IS NOT THE SAME PROCESS. The measured constant is "
            "ANOXIC BIOREDUCTIVE METABOLISM of a hypoxia-activated prodrug — a "
            "real, drug-specific consumption reaction that the engine's panel "
            "drugs do not have — while the engine's k is a generic "
            "binding-plus-clearance sink. Same units, same slot in the same "
            "equation, DIFFERENT PHYSICAL PROCESS. "
            "results/external_anchor_survey.md sec. 2.1 declares it a measured "
            "ORDER-OF-MAGNITUDE BOUND on how large a real interstitial "
            "consumption rate can be, and explicitly NOT a like-for-like "
            "score. (2) THE ENGINE HAS NO SINGLE NUMBER TO SCORE. k_raw and "
            f"k_calibrated differ by {k_cal / k_raw:.0f}x and straddle the "
            "measurement: the uncalibrated sink is 53-89x TOO SMALL (which "
            "independently reproduces results/pk_reality.md's ~300x diagnosis "
            "on a different drug, in a different system, against a directly "
            "measured rate constant rather than a penetration depth — that "
            "FAIL-shaped fact is recorded here, not suppressed), while the "
            "alpha-amplified sink lands within 1.2-2.0x. Picking the flattering "
            "reading would manufacture a PASS and picking the other would "
            "manufacture a FAIL, and nothing about this anchor determines the "
            "choice. NOTE WHICH WAY THIS COSTS — AND THAT IT REVERSED. Until "
            "results/two_route_promotion.md, scoring this row on k_calibrated "
            "would have PASSED the repo's 3x band (0.49-0.82x measured at "
            "alpha = 148.6082), so declining to score removed a PASS from the "
            "tally. The promotion cut alpha to 22.2235 — because doxorubicin's "
            "compiled D fell 6.6x and alpha is fitted from L = sqrt(D/k) "
            "against the same measured anchor — and k_calibrated fell with it. "
            f"Scored on k_calibrated TODAY this row would FAIL: "
            f"{k_cal / k_hi:.3f}-{k_cal / k_lo:.3f}x measured, i.e. "
            f"{k_lo / k_cal:.1f}-{k_hi / k_cal:.1f}x TOO SMALL, outside the 3x "
            "band. So declining to score now removes a FAILURE from the tally, "
            "not a pass. The decision not to score is unchanged and rests on "
            "reason (1), which is about commensurability and is untouched by "
            "any of this; what changed is which way the omission flatters us, "
            "and that is stated here rather than left as a stale claim."),
        note=f"WHAT MAY HONESTLY BE SAID, AND IT IS NOW WEAKER THAN IT WAS. "
             f"Before the two-route promotion the magnitude alpha implied for a "
             f"small-molecule interstitial loss rate landed within 2x of one "
             f"independently measured rate constant of the same class, in a fit "
             f"that saw no tirapazamine, no MCL and no rate constant — weaker "
             f"than validation, stronger than nothing "
             f"(results/external_anchor_survey.md sec. 3.3). At the new alpha it "
             f"lands {k_lo / k_cal:.1f}-{k_hi / k_cal:.1f}x below that "
             f"measurement instead. THIS IS THE COST SIDE OF THE PROMOTION AND "
             f"IT IS NOT COSMETIC: the same change that brought the "
             f"dimensionless hindrance row inside its band pushed this reading "
             f"out of its own. The two anchors are coupled — the Primeau "
             f"penetration fit pins the RATIO D/k, so lowering D lowers k by the "
             f"same factor — and they now disagree about where the engine sits. "
             f"That disagreement is unresolved and is the sharpest open question "
             f"the promotion leaves (results/two_route_promotion.md sec. 8).",
        conditional=("mcl-no-vasculature", "alpha-fit", "coefficient-not-field"),
        values={"k_raw_per_s": k_raw, "k_calibrated_per_s": k_cal,
                "alpha": alpha,
                "measured_k_per_s": {k: float(v) for k, v in k_meas.items()},
                "raw_ratio_to_measured": [k_raw / k_hi, k_raw / k_lo],
                "calibrated_ratio_to_measured": [k_cal / k_hi, k_cal / k_lo]},
    ))
    return rows


def _run_brain_slice_capacity() -> List[Row]:
    """The engine's implied tissue drug CAPACITY vs measured rat brain-slice
    ``Vu,brain`` — TEST B of ``results/transport_anchor_prereg.md``.

    THE CRITERION WAS FIXED BEFORE THE NUMBER WAS COMPUTED
    ------------------------------------------------------
    This is the repo's first genuinely pre-registered external test. The plan
    (``results/transport_anchor_prereg.md``) was written and committed by one
    agent; the dataset cascade, the statistic, the 3.0 bar, the bootstrap and
    the withdrawal conditions were all fixed in it; a different agent then
    executed it and published ``results/brain_slice_anchor.md``. Every bar
    below is transcribed from the plan by
    ``experiments/brain_slice_anchor.py``, which also mirrors them into the
    frozen ``results/sweep/transport_anchor_prereg.json`` (prereg sec. 7.3).
    Nothing here may be re-chosen now that the result is known.

    THE HEADLINE IS ENGINE-INDEPENDENT, AND IT IS NOT DROPPABLE
    -----------------------------------------------------------
    ``GMFE_floor`` — the best unsigned geometric-mean fold error ANY constant
    predictor could achieve on this measurement — is the number prereg
    sec. 3.5a insists be published whatever the verdict. It contains no engine
    quantity, no descriptor and no code version, so it does not move under any
    recalibration of ``alpha``, any resolution of the ``HD``-convention
    question, or any bug fix. It is what quantifies the missing per-compound
    term: a re-centred constant would pass THE TEST by construction and would
    still be wrong by ``GMFE_floor`` on average, so re-centring must never be
    reported as a fix. The floor is therefore carried in ``measured``, in the
    ``note`` and in ``values``, and a test pins that it cannot be dropped.

    SCOPE IS ``analytic``, AND THAT IS NOT A FORMALITY
    --------------------------------------------------
    This row scores an INPUT COEFFICIENT — a capacity term compiled in closed
    form — not a solved field. Prereg sec. 8.5 registers it as
    ``SCOPE_ANALYTIC`` and says **never** ``SCOPE_PDE``. The executor verified
    the solver was not entered (both call counters 0), and the runner
    re-verifies it live on every run and records the counts in ``values``.
    Tagging it ``PDE`` would move the summary's "N of N scored rows never
    invoke the transport PDE" line without a single extra solve, which is the
    disclosure that line exists to make.

    NOTHING IS TUNED HERE
    ---------------------
    The measured values come from ``data_cache/brain_slice/vu_brain.csv``,
    whose ``INCLUSION_RULE.md`` was frozen before extraction; the cascade,
    statistics and engine quantity all come from
    ``experiments/brain_slice_anchor.py`` unmodified. The verdict is recomputed
    LIVE rather than quoted, so if the engine's capacity term ever moves, this
    row moves with it instead of going stale.
    """
    import numpy as np

    B = _import_experiment("brain_slice_anchor")

    all_rows = B._load_rows()
    scored = B.cascade(all_rows)["scored"]
    v = B._values(scored)

    # prereg sec. 2.3 / 9.5 rider 3. A drift here does NOT soften the verdict:
    # the test is recomputed live either way, and the drift is shouted in the
    # note. Converting a FAIL to NOT TESTABLE on a constant change is exactly
    # the "check that stopped checking" pattern this row exists to avoid.
    try:
        B.assert_frozen_constants()
        frozen_drift = ""
    except AssertionError as exc:          # pragma: no cover - drift is a defect
        frozen_drift = str(exc)

    # GATE 1 (prereg sec. 3.3), with the solver guard the executor used.
    with B._SolverGuard() as guard:
        gate = B.engine_tissue_capacity()
    solver_calls = dict(guard.counts)

    pred = float(gate["implied_vu_brain_ml_per_g"])
    floor, c_opt = B.gmfe_floor(v)
    _f, f_lo, f_hi, _fz0, _fa, _ = B.bca(
        v, lambda x: B.gmfe_floor(x)[0], B.BOOTSTRAP_RESAMPLES,
        B.BOOTSTRAP_SEED, B.BOOTSTRAP_CI)
    g_eng = B.gmfe(pred, v)
    ratio = g_eng / floor
    _r, r_lo, r_hi, z0, acc, _ = B.bca(
        v, lambda x: B.gmfe(pred, x) / B.gmfe_floor(x)[0],
        B.BOOTSTRAP_RESAMPLES, B.BOOTSTRAP_SEED, B.BOOTSTRAP_CI)

    bar_x = float(B.THE_TEST_BAR)
    # prereg sec. 3.6: the interval must lie ENTIRELY on one side of the bar. A
    # straddle is UNDERPOWERED (bucket B-C) and is NOT a soft pass.
    if r_lo > bar_x:
        status, bucket = FAIL, "B-B MISCALIBRATED"
    elif r_hi < bar_x:
        status, bucket = PASS, "B-A CONSTRAINED"
    else:
        status, bucket = NOT_TESTABLE, "B-C UNDERPOWERED"

    lo_m, hi_m = float(v.min()), float(v.max())
    bw1 = not (lo_m <= pred <= hi_m)
    below = int(np.sum(pred < v))
    within3 = B._coverage(pred, v, 3.0) * 100.0
    within10 = B._coverage(pred, v, 10.0) * 100.0
    pmids = tuple(dict.fromkeys(r["source_pmid"] for r in scored))

    note = (
        f"GMFE_floor = {floor:.2f}x IS THE HEADLINE, AND IT IS "
        f"ENGINE-INDEPENDENT (prereg sec. 3.5a). {floor:.2f}x "
        f"(90 % BCa CI {f_lo:.2f}-{f_hi:.2f}) is the best unsigned "
        f"geometric-mean fold error ANY constant predictor could achieve on "
        f"these {len(v)} measurements, attained at the constant "
        f"{c_opt:.0f} mL.g brain^-1. It contains no engine quantity and does "
        f"not move under any recalibration of alpha, any resolution of the "
        f"HD-convention question, or any bug fix. Read plainly: the engine "
        f"has NO per-compound tissue-capacity term at all, and this is what "
        f"that costs. A single re-centred constant would cut GMFE(engine) "
        f"from {g_eng:.2f}x to {floor:.2f}x and would then PASS THE TEST BY "
        f"CONSTRUCTION while still being wrong by {floor:.2f}x on average -- "
        f"nearly {floor / B.REPRODUCIBILITY_FLOOR:.0f}x the assay's own "
        f"{B.REPRODUCIBILITY_FLOOR:.2f}x reproducibility floor. RE-CENTRING "
        f"IS THEREFORE NOT A FIX AND MUST NOT BE REPORTED AS ONE; the "
        f"headroom a better CONSTANT is competing for is zero, and the "
        f"headroom a per-compound term is competing for is "
        f"{floor:.2f}x down to {B.REPRODUCIBILITY_FLOOR:.2f}x. "
        f"DIRECTION: the miss is one-directional -- the engine is below "
        f"{below} of {len(v)} measured values (B-S4), within-3x coverage "
        f"{within3:.1f} %, within-10x {within10:.1f} %, worst case "
        f"{float(np.max(np.maximum(pred / v, v / pred))):.0f}x -- and the "
        f"non-base stratum's ratio is 6.58, so prereg sec. 4.3's "
        f"'concentrated in the base stratum' exemption was checked and does "
        f"NOT apply. WITHDRAWAL CONDITION BW1 IS "
        f"{'TRIGGERED' if bw1 else 'NOT TRIGGERED'}: the engine's implied "
        f"value {pred:.3f} lies "
        f"{'OUTSIDE' if bw1 else 'inside'} the retained measured range "
        f"[{lo_m:g}, {hi_m:g}]. The claim BW1 withdraws, verbatim from prereg "
        f"sec. 4.2, is 'the engine's treatment of tissue drug capacity is a "
        f"physically defensible approximation, rather than a placeholder "
        f"whose value was never checked against anything' -- WITHDRAWN, not "
        f"narrowed. "
        f"TWO JUDGEMENT CALLS ARE RECORDED RATHER THAN RESOLVED, both flagged "
        f"by the executor and both carried forward here. (1) GATE 1's ENGINE "
        f"QUANTITY. prereg sec. 3.3 lists four candidates and rejects three; "
        f"the quantity actually scored, the two-route capacity term s_eff, is "
        f"not on that list, and sec. 8.4 hands the slot to the executor. THE "
        f"ALTERNATIVE READING, STATED SO THE DECISION IS VISIBLE AND "
        f"REVERSIBLE: a reader could hold that a capacity fixed at 1 by a "
        f"NORMALISATION CONVENTION is not a quantity the engine 'expresses' "
        f"at all, in which case GATE 1 fails and this bucket is B-D, NOT "
        f"TESTABLE, reason 'the model cannot express the measured quantity'. "
        f"THIS ROW REGISTERS THE LESS FLATTERING READING -- FAIL -- because "
        f"s_eff = (1-phi) + phi*K*(1+K_b) is a real expression whose "
        f"parameters happen to be unity and it passes all four of sec. 3.3's "
        f"gate items on its own terms. Both readings are non-passes and "
        f"neither supports the claim, and the NUMBER is unaffected either "
        f"way: cell_partition, two_route.retardation_factor, the unshipped "
        f"ion_trapping R and s_eff are ALL exactly 1.0 as shipped. (2) PREREG "
        f"sec. 4.2 IS INTERNALLY INCONSISTENT ABOUT BW1 AND NEEDS AMENDMENT "
        f"UNDER sec. 9.5 -- which an implementer may not do, so it is flagged "
        f"here and not touched. BW1's operative sentence says 'outside "
        f"[min, max] of the RETAINED measurements' (retained min "
        f"{lo_m:g} -> trips), while its illustrative paragraph quotes the "
        f"pre-E6 range 0.436-4184, which straddles unity (-> would not trip). "
        f"The executor applied the operative sentence because 'retained' has "
        f"a single fixed meaning in that document. NOTHING TURNS ON IT: THE "
        f"TEST fails by {ratio / bar_x:.1f}x under either reading, and BW1 is "
        f"a separate condition from THE TEST. "
        f"The engine prediction has ZERO VARIANCE -- one distinct capacity "
        f"state over 210 descriptor combinations spanning logP -4..8, HBD "
        f"0..16 and radius 0.15-3.0 nm -- so Spearman rho is NOT COMPUTABLE "
        f"(prereg sec. 3.5 registered this in advance); it is not 0.00 and "
        f"not a separate failure. Sensitivity: at rho_tissue = 1.04 g.mL^-1 "
        f"the engine value is "
        f"{gate['implied_vu_brain_ml_per_g_at_1p04']:.3f} and every verdict "
        f"is unchanged. Every sensitivity arm the executor ran also fails "
        f"(E1/E5/E6 relaxed: ratios 10.94, 12.17, 10.45), so no exclusion is "
        f"doing load-bearing work; the separate, never-pooled mouse arm sits "
        f"at 3.11 and carries no verdict (prereg sec. 3.4 E2). Reproducibility "
        f"floor {B.REPRODUCIBILITY_FLOOR:.2f}x (n = 2 replicate pairs; an "
        f"order-of-magnitude bound, not a calibrated uncertainty) prints "
        f"beside every fold-error figure above, per prereg sec. 3.4a."
    )
    if frozen_drift:                        # pragma: no cover - drift is a defect
        note = ("!! prereg sec. 2.3 FROZEN CONSTANT DRIFT: " + frozen_drift +
                " The verdict above was recomputed LIVE against the moved "
                "engine and is NOT softened, but this row is no longer the "
                "pre-registered execution and prereg sec. 9.5 requires the "
                "coordinator to re-register it. " + note)
    if any(solver_calls.values()):          # pragma: no cover - a scope defect
        note = (f"!! SCOPE TAG DEFECT: the transport solver WAS entered while "
                f"producing this row ({solver_calls}), so scope=analytic is "
                f"no longer true of it. The verdict stands; the tag does not. "
                + note)

    return [Row(
        group="brain-slice-tissue-capacity",
        label="implied tissue capacity Vu,brain (33 compounds, rat brain slice)",
        status=status, scope=SCOPE_ANALYTIC,
        quantity="the engine's implied unbound brain tissue distribution "
                 "volume Vu,brain (mL.g brain^-1) -- the two-route capacity "
                 "term s_eff = (1-phi) + phi*K*(1+K_b) from "
                 "molecular/two_route.py at the shipping constants, converted "
                 "at a density of 1.00 g.mL^-1 declared IN ADVANCE -- vs "
                 "MEASURED rat brain-slice Vu,brain, the same ratio of tissue "
                 "content to unbound interstitial concentration, measured by "
                 "the equilibrium brain-slice method",
        predicted=f"{pred:.3f} mL.g brain^-1, THE SAME NUMBER FOR EVERY "
                  f"COMPOUND (s_eff = "
                  f"{gate['s_eff_ml_per_ml_tissue']:.6f} exactly). "
                  f"GMFE(engine) = {g_eng:.2f}x, i.e. ratio to the floor = "
                  f"{ratio:.3f} against a bar of {bar_x:.1f} "
                  f"[90 % BCa CI {r_lo:.2f}-{r_hi:.2f}, entirely "
                  f"{'ABOVE' if r_lo > bar_x else 'below'} the bar] -- bucket "
                  f"{bucket}. The engine sits {c_opt / pred:.1f}x below the "
                  f"best available constant and below "
                  f"{below}/{len(v)} measurements",
        measured=f"{len(v)} retained rows, {lo_m:g}-{hi_m:g} mL.g brain^-1 "
                 f"(spread {hi_m / lo_m:.0f}x), rat, equilibrium brain slice, "
                 f"from open-access primary sources. GMFE_floor = "
                 f"{floor:.2f}x (90 % BCa CI {f_lo:.2f}-{f_hi:.2f}) at the "
                 f"optimal constant {c_opt:.0f} mL.g brain^-1 -- the best ANY "
                 f"constant predictor could do, and the number that "
                 f"quantifies the missing per-compound term",
        pmids=pmids,
        bar=f"PRE-REGISTERED IN ADVANCE, results/transport_anchor_prereg.md "
            f"sec. 3.5: GMFE(engine) / GMFE_floor <= {bar_x:.1f}, with the "
            f"90 % BCa interval ({B.BOOTSTRAP_RESAMPLES} resamples over "
            f"compounds, seed {B.BOOTSTRAP_SEED}) required to lie ENTIRELY on "
            f"one side of the bar (sec. 3.6); a straddle is UNDERPOWERED and "
            f"NOT a soft pass. The bar is baseline-relative and forgiving BY "
            f"CONSTRUCTION on a set this wide, and it was fixed and committed "
            f"before any prediction was computed -- it has not been widened, "
            f"and prereg sec. 4.1 says the required next action on a miss is "
            f"to constrain the coefficient, NOT to widen the bar",
        reason=("" if status != NOT_TESTABLE else
                f"UNDERPOWERED (bucket B-C, prereg sec. 3.6): the 90 % BCa "
                f"interval [{r_lo:.2f}, {r_hi:.2f}] straddles the bar "
                f"{bar_x:.1f}, so the sign of the comparison is not "
                f"established. This is NOT a soft pass and must never be "
                f"pooled with one."),
        note=note,
        conditional=("coefficient-not-field", "capacity-normalised-to-one",
                     "brain-slice-not-tumour"),
        values={
            "verdict_bucket": bucket,
            "engine_vu_brain_ml_per_g": pred,
            "s_eff": float(gate["s_eff_ml_per_ml_tissue"]),
            "gmfe_engine_x": float(g_eng),
            "gmfe_floor_x": float(floor),
            "gmfe_floor_ci90": [float(f_lo), float(f_hi)],
            "gmfe_floor_optimal_constant_ml_per_g": float(c_opt),
            "ratio": float(ratio),
            "ratio_ci90": [float(r_lo), float(r_hi)],
            "the_test_bar": bar_x,
            "n": int(v.size),
            "measured_range_ml_per_g": [lo_m, hi_m],
            "reproducibility_floor_x": float(B.REPRODUCIBILITY_FLOOR),
            "bw1_reachability_triggered": bool(bw1),
            "n_below_engine": below,
            "within_3x_pct": float(within3),
            "within_10x_pct": float(within10),
            "spearman_rho": "NOT COMPUTABLE -- zero-variance prediction",
            "density_g_per_ml": float(B.TISSUE_DENSITY_G_PER_ML),
            "bca": {"resamples": int(B.BOOTSTRAP_RESAMPLES),
                    "seed": int(B.BOOTSTRAP_SEED),
                    "level": float(B.BOOTSTRAP_CI),
                    "z0": float(z0), "a": float(acc)},
            "solver_calls": solver_calls,
            "frozen_constants_drift": frozen_drift,
            "prereg": "results/transport_anchor_prereg.md",
            "prereg_frozen_json": "results/sweep/transport_anchor_prereg.json",
            "execution_record": "results/brain_slice_anchor.md",
        },
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
    from .kill_calibration import band_is_degenerate

    rows: List[Row] = []
    for name in sorted(tests):
        t = tests[name]
        d = per_drug.get(name, {})
        # A degenerate (zero-width) band is refused a verdict here as well as in
        # kill_calibration.ratio_test, so that NO path into this renderer — GBM,
        # breast, or one added later — can print PASS/FAIL against a bar nothing
        # can pass. Same treatment the GBM wick-in-needle IFP row gets for the
        # mirror defect, a bar nothing can fail.
        degenerate = (bool(t.get("degenerate"))
                      or t.get("verdict") == NOT_TESTABLE
                      or band_is_degenerate(t["lit_ratio_lo"], t["lit_ratio_hi"]))
        status = (NOT_TESTABLE if degenerate
                  else (PASS if t["verdict"] == PASS else FAIL))
        matched = t.get("scope", "")
        note = f"literature scope used: {matched}"
        if matched == "solid":
            note += (" — NO drug-matched 3-D:2-D band exists, so this falls back "
                     "to the broad cross-tumour-type cytotoxic band")
        if t["verdict"] not in (PASS, NOT_TESTABLE):
            note += f". {t['direction']}"
        note += f". {tissue_note}"
        rows.append(Row(
            group=group, label=name, status=status, scope=SCOPE_PDE,
            reason=(t.get("reason") or
                    "ZERO-WIDTH literature band: no prediction can land inside "
                    "it, so no verdict is rendered.") if degenerate else "",
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

    **This function decides nothing.** Which band a drug is scored against, and
    whether it can be scored at all, is decided once by
    ``kill_calibration.score_ratio_row`` and carried on
    ``build_ratio_section(...)[drug]["scored"]``; the report generator that
    writes ``results/kill_calibration_breast.md`` renders the same object. That
    is deliberate: this function and that generator used to derive the verdict
    separately and disagreed about this exact panel — the report published
    5 PASS / 1 FAIL while ``codon validate`` printed 3 PASS / 3 FAIL, two of the
    FAILs against zero-width bands. See ``results/zero_width_bands.md``.

    Two kinds of NOT TESTABLE reach this scorecard, and they are different
    facts about the literature:

    * **No published 3-D-vs-2-D comparison in any tumour type** (cyclophosphamide,
      tamoxifen, vinorelbine, palbociclib) — no measured target exists.
    * **Exactly one published comparison** (paclitaxel, 5-fluorouracil) — a
      measured target exists but it is a POINT, and a min-max over one point is a
      zero-width interval. Not re-scored against the broad non-drug-matched band,
      which is the same band the four drugs above are refused a verdict on.
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
        s = row["scored"]
        if not s["scored"]:
            band = s["band"]
            ctx = s["context_test"]
            measured = "NONE"
            if band is not None:
                measured = (f"a single measured point, fold "
                            f"{band['fold_central']:.4g}x -> zero-width ratio "
                            f"band {band['ratio_lo']:.4f}-{band['ratio_hi']:.4f}")
            reason = s["reason"]
            if ctx is not None:
                reason += (f" For context only, and NOT a verdict: against the "
                           f"broad non-drug-matched 'solid' band "
                           f"({s['context_band']['fold_lo']:.4g}-"
                           f"{s['context_band']['fold_hi']:.4g}x, n="
                           f"{s['context_band']['n_rows']}) this ratio would "
                           f"score {ctx['verdict']}.")
            rows.append(Row(
                group="breast-kill-3d2d", label=name, status=NOT_TESTABLE,
                scope=SCOPE_PDE,
                quantity="predicted 3-D:2-D kill ratio vs MEASURED "
                         "spheroid-vs-monolayer IC50 fold shift",
                predicted=f"ratio {h['ratio_exact']:.4f}",
                measured=measured,
                # The predicted ratio IS a solver output — the row is
                # verdictless because no bar exists, not because nothing ran —
                # so the delivery law is declared rather than stamped "not
                # applicable".
                physics=physics_label(repo_default_exchange()),
                pmids=tuple(_rows_pmids(row.get("specific_bridge"))), bar="n/a",
                reason=reason,
                conditional=("alpha-fit", "gdsc-potency"),
                values={"predicted_ratio": float(h["ratio_exact"]),
                        "n_literature_rows": int(s["n_rows"])},
            ))
            continue
        t = dict(s["test"])
        t.setdefault("scope", f"drug:{name}")
        t["fold_lo"] = s["band"]["fold_lo"]
        t["fold_hi"] = s["band"]["fold_hi"]
        t["implied_fold"] = h.get("implied_ic50_fold", float("nan"))
        rows.extend(_kill_rows(
            "breast-kill-3d2d",
            "Real untreated Janesick Xenium breast section (PMID 38114474), "
            f"dx = 12 um, 6 h window. Drug-matched band, n={s['n_rows']}.",
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


def _run_ifp_wick_in_needle() -> List[Row]:
    """Solved interstitial fluid pressure vs MEASURED human wick-in-needle IFP.

    The engine solves IFP natively in mmHg (``physics/pressure.py``: ``_PV =
    25.0`` mmHg, ``K`` in um^2/(mmHg*s)), with no conversion at the display
    layer, so it is one of the very few ABSOLUTE dimensioned quantities the
    engine emits. Until now nothing external scored it, while the landing page
    displayed a specific mmHg range taken from it.

    Three rows, and the split between them is the whole point.

    1. **GBM vs Boucher 1997 (PMID 9062403) -- NOT TESTABLE.** Not because the
       units or the technique fail to line up (they do line up), but because
       the ANCHOR CANNOT DISCRIMINATE. The measured intracranial-tumour value is
       2.0 +/- 2.5 mmHg where the 2.5 is a **standard deviation** (the paper's
       Statistical analysis section says so verbatim), n = 11, on a
       non-negative quantity whose solver clips at zero -- so the admissible
       band is effectively 0 to 4.5 mmHg and a model that produced NOTHING
       would pass it. The killer is the cross-check: Less 1992 (PMID 1423283),
       same wick-in-needle technique with Boucher and Jain on both author
       lists, measured NORMAL tissue at 2 +/- 4 mmHg, n = 11 -- statistically
       indistinguishable from the intracranial TUMOUR figure. An anchor that
       cannot separate tumour from normal tissue cannot separate this engine
       from one with no tumour in it. Registering that as a PASS would be the
       weakest row on the scorecard while reading like one of the strongest.
    2. **Breast vs Less 1992 (PMID 1423283) -- the DISCRIMINATING row, and it
       FAILS.** Primary breast carcinoma IFP is 15 +/- 9 mmHg (SD), far from
       zero and far from normal tissue, so a wrong model cannot pass by
       producing nothing. The engine lands ~40x low.
    3. **A numerical-resolution row**, which is why row 2 fails the way it
       does: the screened-Poisson operator's screening length
       ``sqrt(K/(a+d))`` is ~6-7 um on both real sections, BELOW every grid
       spacing the engine ever runs at, so > 99.7% of tissue nodes are
       unresolved and the field collapses toward the pointwise algebraic limit
       ``a*p_v/(a+d)``.

    NOTHING IS TUNED HERE. ``experiments/ifp_validation.py`` imports
    ``physics/pressure.py`` and calls it as it ships; no coefficient is read
    back, fitted or adjusted toward any measured value.
    """
    for need in ("data_cache/xenium_gbm", "data_cache/xenium_breast_janesick/rep1"):
        if not os.path.isdir(_repo_path(*need.split("/"))):
            raise NotTestable(f"requires the real Xenium cache at {need} (absent)")
    IV = _import_experiment("ifp_validation")

    gbm = IV.measure_section("gbm", 12.0)
    breast = IV.measure_section("breast", 12.0)
    sc = IV.score([gbm, breast])
    A, B = IV.ANCHOR, IV.BREAST_ANCHOR

    def _wick(rec):
        return rec["wick"]["sigma_250um"]["tumour_phi_ge_0.5"]

    g, b = _wick(gbm), _wick(breast)
    commensurability = (
        "COMMENSURABILITY. A wick-in-needle probe is a 23-gauge needle whose "
        "side hole equilibrates with the mobile interstitial fluid over a "
        "sub-millimetre to millimetre contact region, at ONE placement inside "
        "the tumour mass. The commensurable model statistic is therefore the "
        "solved field's local mean over a wick-sized window at tumour-bearing "
        "locations -- NOT the domain mean, NOT the field maximum, and NOT the "
        "landing page's iso-contour levels, which are fractions of a "
        "coarse-grained maximum chosen for rendering. Reported here as the "
        "median over placements at a 250 um Gaussian window over nodes with "
        "tumour occupancy >= 0.5; the 125 um and 500 um windows are in "
        "results/sweep/ifp_validation.json and move it by under 1.4x.")
    unresolved = (
        f"The screened-Poisson screening length sqrt(K/(a+d)) has median "
        f"{gbm['resolution']['screening_length_um']['p50']:.1f} um (GBM) and "
        f"{breast['resolution']['screening_length_um']['p50']:.1f} um (breast), "
        f"BELOW the 12 um production spacing at "
        f"{100 * gbm['resolution']['frac_tissue_screening_below_dx']:.2f}% and "
        f"{100 * breast['resolution']['frac_tissue_screening_below_dx']:.2f}% "
        f"of tissue nodes respectively.")

    rows: List[Row] = []

    rows.append(Row(
        group="ifp-wick-in-needle",
        label="GBM section IFP vs human intracranial-tumour IFP",
        status=NOT_TESTABLE, scope=SCOPE_PDE,
        quantity="solved interstitial fluid pressure (mmHg) at wick-scale "
                 "placements in the tumour mass of the real Xenium GBM section "
                 "vs MEASURED intraoperative wick-in-needle IFP in human "
                 "intracranial tumours",
        predicted=f"{g['p50_mmhg']:.3f} mmHg (median over placements; IQR "
                  f"{g['p25_mmhg']:.3f}-{g['p75_mmhg']:.3f}, max "
                  f"{g['max_mmhg']:.3f})",
        measured=f"{A['mean_mmhg']} +/- {A['spread_mmhg']} mmHg "
                 f"(mean +/- SD, n={A['n']}); 10 of 11 tumours lay in 0.0-3.0, "
                 f"one astrocytoma at 9.5; glioblastoma-only subset n=5 mean "
                 f"{sc['anchor_recomputed']['glioblastoma_only_mean_mmhg']:.1f}, "
                 f"range 0.0-3.0",
        pmids=(A["pmid"], B["pmid"]),
        bar="n/a — no bar is declared, because no bar on this anchor "
            "discriminates",
        reason=(
            "THE ANCHOR CANNOT DISCRIMINATE, for three independent reasons, and "
            "the engine's number would sit inside it either way. "
            "(1) UNCERTAINTY. The paper's Statistical analysis section states "
            "verbatim that 'the data are presented as the mean +/- standard "
            "deviation', so 2.5 is an SD over 11 patients, not an SEM. On a "
            "non-negative quantity that the solver clips at zero, mean +/- 1 SD "
            "is effectively 0.0-4.5 mmHg, and the paper's own per-patient table "
            "has 10 of 11 tumours in 0.0-3.0 INCLUDING A MEASURED 0.0. Any "
            "engine output from nothing up to 4.5 mmHg passes. A verdict that "
            "cannot be failed is not a verdict. "
            "(2) IT DOES NOT SEPARATE TUMOUR FROM NORMAL TISSUE. Less 1992 "
            "(PMID 1423283) measured normal human tissue at 2 +/- 4 mmHg, n=11, "
            "by the SAME wick-in-needle technique with Boucher and Jain on both "
            "author lists. The intracranial-tumour figure is statistically "
            "indistinguishable from it, so agreement with it is not evidence "
            "that the engine reproduces tumour interstitial hypertension -- it "
            "is equally consistent with the engine having no tumour effect at "
            "all, which on this section is close to what it has. "
            "(3) TREATED, DECOMPRESSED STATE. All 11 patients received "
            "dexamethasone 10 mg + mannitol 0.5 mg/kg + furosemide 20 mg i.v. "
            "and were hyperventilated to PaCO2 25-32 mmHg, and IFP was taken "
            "after the cranium was opened and the dura removed. The authors' "
            "own conclusion is that 'antioedema therapy in patients with brain "
            "tumours is responsible in part for the low tumour IFP', and the "
            "same paper measures a 17-50% fall in rodent tumour IFP from "
            "hypocapnia alone (27 +/- 15 -> 18 +/- 10 mmHg). The engine models "
            "no steroid, no osmotic agent, no diuretic and no decompression, so "
            "even a matching number would be matching a different physiological "
            "state. " + commensurability),
        note=("Not a null result: the DISCRIMINATING version of this same "
              "comparison is the breast row below, which the engine FAILS by "
              "~40x. Read the two together. " + unresolved),
        # ``ifp-grid-unconverged`` belongs on THIS row too, not only on the
        # breast one below. The caveat's own text describes the GBM ladder
        # (0.358 -> 0.434 -> 0.462 -> 0.370 mmHg, non-monotone, 1.29x spread),
        # and the median quoted above is one point on it. A row that declines
        # to score still PUBLISHES a number, and that number is no more
        # grid-converged than the breast one.
        conditional=("ifp-uncalibrated", "ifp-unresolved-mesh",
                     "ifp-grid-unconverged", "untreated-section"),
        values={"wick250_tumour_median_mmhg": g["p50_mmhg"],
                "wick250_tumour_iqr_mmhg": [g["p25_mmhg"], g["p75_mmhg"]],
                "raw_field_max_mmhg": gbm["raw"]["max_mmhg"],
                "raw_field_median_in_tissue_mmhg": gbm["raw"]["p50_mmhg"],
                "raw_field_mean_in_tissue_mmhg": gbm["raw"]["mean_in_tissue_mmhg"],
                "anchor_mean_mmhg": A["mean_mmhg"],
                "anchor_sd_mmhg": A["spread_mmhg"],
                "normal_tissue_reference_mmhg": [2.0, 4.0],
                "dx_um": 12.0},
        physics=physics_label(repo_default_exchange()),
    ))

    blo, bhi = B["mean_mmhg"] - B["spread_mmhg"], B["mean_mmhg"] + B["spread_mmhg"]
    rows.append(Row(
        group="ifp-wick-in-needle",
        label="breast section IFP vs human primary-breast-carcinoma IFP",
        status=classify(b["p50_mmhg"], blo, bhi), scope=SCOPE_PDE,
        quantity="solved interstitial fluid pressure (mmHg) at wick-scale "
                 "placements in the tumour mass of the real Xenium breast "
                 "section vs MEASURED wick-in-needle IFP in human primary "
                 "breast carcinoma",
        predicted=f"{b['p50_mmhg']:.3f} mmHg (median over placements; IQR "
                  f"{b['p25_mmhg']:.3f}-{b['p75_mmhg']:.3f}, max "
                  f"{b['max_mmhg']:.3f}) — "
                  f"{sc['breast_vs_breast_anchor']['fold_low_vs_mean']:.0f}x "
                  f"BELOW the measured mean, "
                  f"{sc['breast_vs_breast_anchor']['fold_low_vs_band_floor']:.0f}x "
                  f"below the band floor",
        measured=f"{B['mean_mmhg']:.0f} +/- {B['spread_mmhg']:.0f} mmHg "
                 f"(mean +/- SD) in primary breast carcinoma; all neoplasms "
                 f"{B['all_neoplasms_mmhg']}; normal tissue "
                 f"{B['normal_tissue_mmhg']}",
        pmids=(B["pmid"], "7850532"),
        bar=f"inside the measured mean +/- 1 SD = {blo:.0f}-{bhi:.0f} mmHg",
        note=("THE LOAD-BEARING ROW OF THIS GROUP, and it is the one the "
              "intracranial anchor could not be. 15 +/- 9 mmHg is far from zero "
              "and 7x the same study's normal-tissue value, so unlike the brain "
              "figure it cannot be passed by producing nothing — and the engine "
              "does not pass it. Corroborated independently by Nathanson 1994 "
              "(PMID 7850532), 29 +/- 3 mmHg in invasive ductal carcinoma, n=25 "
              "— whose uncertainty is a STANDARD ERROR, not an SD, and must not "
              "be read as one; against that value the engine is ~78x low. "
              "The direction is consistent with every other absolute quantity "
              "this engine has been scored on. Mechanistically the miss is not "
              "mysterious: pressure.py's own docstring predicts a flat tumour "
              "plateau approaching _PV = 25 mmHg, but on real reconstructed "
              "tissue the vascular density that drives the Starling source is "
              "zero at 94.3% of breast tissue nodes (97.3% on GBM), so the "
              "solved field is near-zero almost everywhere (median "
              f"{breast['raw']['p50_mmhg']:.4f} mmHg in tissue) with isolated "
              f"spikes to {breast['raw']['max_mmhg']:.1f} mmHg at endothelial "
              "cells. The engine reproduces the peak magnitude and misses the "
              "plateau, which is the opposite of the Jain effect it claims. " +
              commensurability),
        conditional=("ifp-uncalibrated", "ifp-unresolved-mesh",
                     "ifp-grid-unconverged", "untreated-section"),
        values={"wick250_tumour_median_mmhg": b["p50_mmhg"],
                "wick250_tumour_iqr_mmhg": [b["p25_mmhg"], b["p75_mmhg"]],
                "band_mmhg": [blo, bhi],
                "fold_below_measured_mean":
                    sc["breast_vs_breast_anchor"]["fold_low_vs_mean"],
                "raw_field_max_mmhg": breast["raw"]["max_mmhg"],
                "raw_field_median_in_tissue_mmhg": breast["raw"]["p50_mmhg"],
                "dx_um": 12.0},
        physics=physics_label(repo_default_exchange()),
    ))

    rows.append(Row(
        group="ifp-wick-in-needle",
        label="IFP solve is under-resolved on every grid the engine runs at",
        status=NOT_TESTABLE, scope=SCOPE_PDE,
        quantity="screened-Poisson screening length sqrt(K/(a+d)) implied by "
                 "physics/pressure.py's own shipped coefficients vs the grid "
                 "spacings the engine actually solves on (8-30 um)",
        predicted=f"median screening length "
                  f"{gbm['resolution']['screening_length_um']['p50']:.1f} um "
                  f"(GBM) / "
                  f"{breast['resolution']['screening_length_um']['p50']:.1f} um "
                  f"(breast); L/dx = "
                  f"{gbm['resolution']['median_L_over_dx']:.2f} / "
                  f"{breast['resolution']['median_L_over_dx']:.2f} at the 12 um "
                  f"production spacing",
        measured="n/a — this row is scored against the discretisation, not "
                 "against a measurement",
        pmids=(), bar="n/a",
        reason=(
            "NO MEASUREMENT IS INVOLVED; this row is registered so that the two "
            "rows above cannot be quoted without it. " + unresolved + " A "
            "screened-Poisson solve whose screening length is below the mesh "
            "cannot represent the operator it discretises: the reaction term "
            "dominates the diagonal (median (a+d)*dx^2/K = "
            f"{gbm['resolution']['reaction_dominance_median']:.1f} on GBM, "
            f"{breast['resolution']['reaction_dominance_median']:.1f} on "
            "breast), so each node relaxes toward its own LOCAL algebraic "
            "equilibrium a*p_v/(a+d) with the Laplacian barely coupling "
            "neighbours. Measured directly: the solved field correlates with "
            "that pointwise algebraic limit at r = "
            f"{gbm['resolution']['corr_solved_vs_local_equilibrium']:.2f} (GBM) "
            f"and {breast['resolution']['corr_solved_vs_local_equilibrium']:.2f} "
            "(breast). At production dx the 'solved pressure field' is "
            "therefore closer to a pointwise nonlinear map of the vascular-"
            "density KDE than to a transport solve, and the margin pressure "
            "gradient that carries the Jain effect into the convection term is "
            "not resolved. Refining the mesh does not fix it and is not free: "
            "the field is not grid-converged either (see the "
            "ifp-grid-unconverged conditional), and pressure.py's tumour-"
            "occupancy smoothing _TUMOR_SIGMA = 5.0 is specified in GRID CELLS, "
            "so the physical smoothing scale is itself a function of dx (60 um "
            "at dx=12, 120 um at dx=24) and the model changes as the mesh "
            "changes."),
        note="Grid ladder, wick-scale tumour median: GBM 0.358 -> 0.434 -> "
             "0.462 -> 0.370 mmHg and breast 0.621 -> 0.424 -> 0.370 -> 0.286 "
             "mmHg at dx = 24, 18, 12, 8 um. Full ladder in "
             "results/sweep/ifp_validation.json.",
        conditional=("ifp-unresolved-mesh", "ifp-grid-unconverged"),
        values={"screening_length_um_gbm":
                    gbm["resolution"]["screening_length_um"],
                "screening_length_um_breast":
                    breast["resolution"]["screening_length_um"],
                "frac_tissue_unresolved_gbm":
                    gbm["resolution"]["frac_tissue_screening_below_dx"],
                "frac_tissue_unresolved_breast":
                    breast["resolution"]["frac_tissue_screening_below_dx"],
                "corr_solved_vs_local_equilibrium_gbm":
                    gbm["resolution"]["corr_solved_vs_local_equilibrium"],
                "corr_solved_vs_local_equilibrium_breast":
                    breast["resolution"]["corr_solved_vs_local_equilibrium"]},
        physics=physics_label(repo_default_exchange()),
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
            key="mcl-transport-coefficients",
            title="Compiled transport coefficients vs measured multicellular-"
                  "layer D, hindrance and consumption rate (Hicks 1998/2007)",
            tissue="MCL", tier=TIER_QUICK, scope=SCOPE_ANALYTIC,
            source="experiments/anchor_prototype.py (owns the measured "
                   "constants, each traced to its PMID) driving the unmodified "
                   "molecular pipeline descriptors.from_smiles -> "
                   "compile_drug.build_profile, which since "
                   "results/two_route_promotion.md APPLIES the two-route cell "
                   "term -- so the verdicts are the cell-term-ON numbers. The "
                   "pre-promotion arm (two_route=False) and the "
                   "impermeable-cell arms are computed live and score nothing",
            report="results/two_route_promotion.md, "
                   "results/external_anchor_survey.md sec. 2.1 + 3, "
                   "results/cellular_exclusion.md, "
                   "results/anchor_registration.md",
            runner=_run_mcl_transport_coefficients, est="~10 s",
        ),
        AnchorGroup(
            key="brain-slice-tissue-capacity",
            title="PRE-REGISTERED (TEST B): implied tissue drug capacity "
                  "Vu,brain vs measured rat brain-slice Vu,brain, 33 compounds",
            tissue="rat brain", tier=TIER_QUICK, scope=SCOPE_ANALYTIC,
            source="experiments/brain_slice_anchor.py (owns the sec. 3.1 "
                   "cascade, the statistics of sec. 3.5 and the GATE 1 engine "
                   "quantity; every bar transcribed from the pre-registration "
                   "and mirrored into results/sweep/transport_anchor_prereg.json "
                   "before any prediction was computed) driving the unmodified "
                   "molecular.two_route capacity term. The transport solver is "
                   "NOT entered -- call counters re-checked live on every run",
            report="results/brain_slice_anchor.md (execution record, FAIL), "
                   "results/transport_anchor_prereg.md (the PRE-REGISTRATION: "
                   "sec. 3.5 THE TEST, sec. 3.5a GMFE_floor, sec. 4.1 buckets, "
                   "sec. 4.2 withdrawal, sec. 8.5 this registration), "
                   "data_cache/brain_slice/PROVENANCE.md",
            runner=_run_brain_slice_capacity, est="~3 s",
            needs=("data_cache/brain_slice/vu_brain.csv",
                   "data_cache/brain_slice/INCLUSION_RULE.md"),
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
        AnchorGroup(
            key="ifp-wick-in-needle",
            title="Solved interstitial fluid pressure vs measured human "
                  "wick-in-needle IFP (GBM not testable, breast FAILS ~40x)",
            tissue="GBM + breast", tier=TIER_FULL, scope=SCOPE_PDE,
            source="experiments/ifp_validation.py::measure_section + score "
                   "(calls tmesim.physics.pressure.solve_ifp exactly as it "
                   "ships; no coefficient is read back or fitted)",
            report="results/ifp_validation.md",
            runner=_run_ifp_wick_in_needle, est="~15 s (2 sections, no drug)",
            needs=("data_cache/xenium_gbm",
                   "data_cache/xenium_breast_janesick/rep1"),
        ),
    ]


#: Validations that were inventoried and DELIBERATELY not registered, with the
#: reason. Printed by ``codon validate --list`` so the omissions are auditable
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
    ("Pruijn 2008 (doi:10.1071/CH08240) — D_mcl for 71 compounds spanning "
     "~69-fold, and the bounded version in PMID 16912264 (2.7-187 um^2/s "
     "across 16 TPZ analogues in HT29 MCLs)",
     "THE HIGHEST-VALUE UNREGISTERED ANCHOR IDENTIFIED, and it is blocked on "
     "retrieval, not on judgement. The engine derives D_eff from a "
     "Stokes-Einstein radius that is a function of MW ALONE times a hindrance "
     "that depends on the ECM fraction and not on the molecule, so across a "
     "congeneric series it must predict a near-constant D_eff while the "
     "measurement spans 69-fold — a decisive falsification of the diffusivity "
     "model that needs only the compound list. The publisher returns HTTP 403 "
     "to automated retrieval and the record is not in PubMed or PMC, so the "
     "compound-by-compound table is UNVERIFIED and no row can be built from it "
     "without inventing numbers (results/external_anchor_survey.md sec. 2.5)."),
    ("Capecitabine tumour:plasma 1.76 and tumour:s.c. 1.11 in human breast "
     "lesions (PMID 12618890)",
     "COMMENSURABLE and recommended for adoption by "
     "results/external_anchor_survey.md sec. 4, but NO MODULE COMPUTES IT: it "
     "would reuse experiments/breast_kp_anchor.py's solve path against a new "
     "anchor table, which is new code in a file this registration does not "
     "own. Registering a row whose runner does not exist is not possible; this "
     "entry records the omission so it stays visible. When it is built, score "
     "the PARENT DRUG ONLY — DFCR, DFUR and 5-FU are metabolites generated "
     "intracellularly en route, so their interstitial concentration carries a "
     "local generation term the engine does not model."),
    ("TEST C of results/transport_anchor_prereg.md (sec. 3.8) -- Kp,uu,cell, "
     "data_cache/brain_slice/kp_uu_cell.csv -- and the two other columns that "
     "sit beside the registered brain-slice anchor: "
     "data_cache/brain_slice/kp_uu_brain_DO_NOT_SCORE.csv, and the MOUSE arm "
     "of vu_brain.csv",
     "REGISTERED IN THE PRE-REGISTRATION AND DELIBERATELY NOT SCORED, each "
     "for a reason fixed before the data were read. TEST C is NOT RUNNABLE: "
     "only 4 scoreable Kp,uu,cell values survive open-access retrieval, below "
     "any power the plan would accept, and prereg sec. 3.8 registers it as "
     "such rather than running it small (it is also the one column that WOULD "
     "constrain cell_partition, which therefore remains unconstrained by any "
     "measurement). Kp,uu,brain is excluded by prereg sec. 3.4 E3 and its file "
     "is named DO_NOT_SCORE: it is a blood-brain-BARRIER ratio, not a tissue "
     "capacity, so scoring it against a capacity term would be the "
     "units/role error prereg sec. 4.2 BW3 exists to catch. The mouse arm is "
     "excluded by E2 -- species are never pooled -- and is reported "
     "separately in results/brain_slice_anchor.md B-S10 (ratio 3.11, no "
     "verdict); it is NOT evidence that the rat result is marginal."),
    ("Tunggal 1999 (PMID 10389947) and Cowan 2001 (PMID 11149410) "
     "multicellular-layer flux ratios on FOUR OF THE ENGINE'S OWN PANEL DRUGS",
     "NOT ACTIONABLE WITHOUT NEW GEOMETRY, and the omission is expensive. "
     "'% of the rate through the blank Teflon membrane' is a ratio of "
     "steady-state fluxes through a 1-D slab with reaction; the engine solves "
     "2-D reconstructed tissue with vessel sources and has no 1-D slab "
     "configuration. The measurement says doxorubicin and mitoxantrone are "
     "3-5x slower than methotrexate while the calibrated panel puts "
     "doxorubicin at L = 64.9 um and methotrexate at 75.2 um — essentially "
     "equal, i.e. a DIRECTIONAL DISAGREEMENT on the engine's own drugs that "
     "nothing currently scores. Cowan's concentration invariance from 10 uM to "
     "3 mM is separately a direct test of the engine's first-order-sink "
     "linearity assumption, and it CONFIRMS it over 300-fold."),
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
