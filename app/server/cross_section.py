"""``GET /api/cross-section`` — ONE COMPOUND, MANY SECTIONS.

The comparison this product can defend, and the one it cannot
--------------------------------------------------------------
Ordering two DIFFERENT COMPOUNDS by predicted effect was withdrawn on
2026-08-21 and is not restored here in any form — not behind a flag, not as a
limit, not as evidence. ``results/pk_ranking_robustness.md``: per-compound
plasma exposure spans **113.1x** across the compounds the GBM panel ranks while
the published ordering first breaks at a spread of only **1.53x**.
``results/capacity_ranking_robustness.md``: the breast top-1 pick changes under
this repo's own structural capacity rule.

The SAME compound on two pieces of tissue is a different question and a
legitimate one. That is the only question this module asks:

    the compound is FIXED, the section VARIES, and nothing else may vary.

``results/cross_tissue.md`` is the measured precedent — one compound
(temozolomide) on five real GBM microenvironments — and it is also the warning.
Its own headline was RETRACTED once because the five tissues had been run at
different grid resolutions, so the published spread "mixed real biology with
numerical error". That retraction is why this module groups runs by
:func:`app.server.projects.spec_signature` and will place only ONE signature on
an axis: same mesh, same modes, same PK, same window, same arm — everything
except the tissue.

What is on the axis
-------------------
The per-section quantity is the **reach ladder** ``F(x)`` from
:mod:`app.server.reach`: the fraction of VIABLE tumour cells whose exposure is
below ``x`` times that section's own perivascular reference. It is a ratio
between two points of ONE map, so it is the *relative* within-tissue statement
that survived the withdrawal — *"Within a single compound on a single tissue,
the RELATIVE exposure between regions is the defensible output"*
(``app.server.projects.RELATIVE_EXPOSURE_CLAIM``) — rather than an absolute
level that does not travel between sections.

``F(0.1)`` is NOT headlined here, and that is deliberate rather than an
oversight: :data:`app.server.reach.CONVENTION_NOTE` records that the threshold
is a reporting convention nothing in the physics marks out, and that retuning it
after seeing the field would be "post-hoc selection wearing a convention's
clothes". The whole ladder is shown, and the verdict below is ORDINAL — it reads
the curve, never one point on it.

THE DISAGREEMENT IS THE OUTPUT
------------------------------
A compound whose reach verdict flips between two patients' tissue is the most
informative thing this engine currently says, so the flip is computed rather
than left for a reader to spot. Two independent detectors, and neither invents a
threshold:

1. **Crossing.** For a pair of sections, the sign of ``F_a(x) - F_b(x)`` at every
   rung both sections can supply. Constant sign at every rung -> the order holds.
   Mixed -> the curves CROSS and there is no ordering of those two tissues at
   all. Ordinal, threshold-free, and a statement about two SECTIONS — never
   about two compounds.

2. **Fraction against count.** ``F`` is a fraction, and its denominator is that
   section's own viable-tumour-cell count. ``results/heldout_patient.md`` §1.1 is
   the archetype and it is why this detector exists at all: between two real GBM
   sections the endothelial *fraction* fell to 0.67x while the endothelial
   *areal density* rose 4.19x, "these two rows point in opposite directions, and
   that is the finding", because section 2 simply had 6.2x the cell density. A
   fraction that falls while the count it is a fraction OF rises is the same
   defect waiting to happen here, so both are shown side by side and a
   direction disagreement is reported in its own block.

What this module REFUSES rather than renders
--------------------------------------------
* **Spot-resolution sections contribute no per-cell statistic.**
  ``app/server/api.py::_SPOT_SOURCES`` and
  ``app/server/api.py::SPOT_SPATIAL_REFUSAL`` are imported and carried verbatim,
  never paraphrased and never re-derived. That refusal states in terms that the
  DELIVERY MAP is invalid on a Visium section too, and for the same reason: the
  exposure field's SOURCE is a density over endothelial cell POSITIONS, and
  those positions are synthesised — 12 pseudo-cells per spot at the spot
  centroid with uniform jitter. Every quantity this view puts on its axis is
  built from that map (see :data:`SPOT_REFUSED_QUANTITIES`), so a spot section
  appears in the view with the refusal IN PLACE OF the numbers.
* **Comparability is a backend verdict, never a display choice.** The absolute
  metrics reuse :func:`app.server.projects._kill_comparability` — the same
  function ``GET /api/projects/{id}/compare`` refuses with — rather than a second
  copy of the rule. Where it says ``comparable: false`` this module emits no
  numbers for that row at all, only the refusal and its blockers.
* **Two signatures never share an axis** (``results/cross_tissue.md``'s
  retraction, above), and neither do two window modes: a steady-state map and a
  finite-window map answer different questions
  (:data:`app.server.reach.STEADY_STATE_NOTE`).
* **Two compounds never share an axis.** :func:`assert_one_compound_axis` walks
  the FINISHED payload and raises on a second compound anywhere in it. It audits
  the output rather than the input precisely because the input filter is the
  thing a future edit would loosen.

Units and denominators
----------------------
``F(x)`` is dimensionless and is a fraction **of that section's own viable
tumour cells**, a count that differs between sections and is itself banded by an
oxygen field whose absolute level has never been calibrated
(``results/cross_tissue.md`` †: a ±40% change in the lumped consumption constant
swings the hypoxic tumour fraction between ~0.0 and ~0.6). Every fraction on the
axis therefore travels with its own denominator, and that caveat travels with
every fraction.
"""
from __future__ import annotations

import time
from typing import Callable, Dict, List, Optional, Sequence, Tuple

from fastapi import APIRouter, HTTPException, Query

from .cache import lookup_section, redact_exception
from .jobs import jsonable
# Imported from the module that already owns these rules. `_drug_of`,
# `_kill_comparability` and `_HEADLINE_METRICS` are private to
# `app/server/projects.py`, and importing a private name is a coupling this
# file accepts DELIBERATELY: the alternative is a second copy of the
# comparability rule, and a second copy of a refusal is how a refusal comes to
# disagree with itself. `tests/test_cross_section.py` pins every one of these
# imports, so a rename in projects.py goes red in the test suite rather than at
# request time.
from .projects import (  # noqa: PLC2701 - see the comment above
    RELATIVE_EXPOSURE_FULL,
    ProjectError,
    _HEADLINE_METRICS,
    _drug_of,
    _kill_comparability,
    capability_report,
    spec_signature,
)
from .reach import HEADLINE_X, REACH_LADDER

__all__ = [
    "CrossSectionError",
    "ONE_COMPOUND_AXIS",
    "WITHDRAWN_CROSS_COMPOUND_AXIS",
    "SPOT_REFUSED_QUANTITIES",
    "DENOMINATOR_CAVEAT",
    "NO_CROSS_SECTION_NOISE_FLOOR",
    "assert_one_compound_axis",
    "available_compounds",
    "build_cross_section_router",
    "cross_section_payload",
]


class CrossSectionError(Exception):
    """A refusal from this module, mapped to an HTTP status by the router."""

    def __init__(self, message: str, status: int = 400,
                 detail: Optional[Dict] = None):
        super().__init__(message)
        self.message = message
        self.status = int(status)
        self.detail = detail or {}


# --------------------------------------------------------------------------- #
# The one sentence this view is, and the one it refuses to become
# --------------------------------------------------------------------------- #
#: What the axis IS. Returned in the payload so the browser renders this
#: sentence rather than inventing its own description of the comparison.
ONE_COMPOUND_AXIS = (
    "ONE COMPOUND, MANY SECTIONS. The compound is held FIXED and the SECTION "
    "varies; that is the only comparison drawn here. Two sections on one axis "
    "is the legitimate question — the same molecule meeting two different "
    "pieces of real tissue. Two COMPOUNDS on one axis is the withdrawn one, "
    "whatever the axis is labelled."
)

#: The refusal, raised by :func:`assert_one_compound_axis`. It names the defect
#: as the BUILDER's, because by the time this text is produced the request was
#: legitimate and the payload was not.
WITHDRAWN_CROSS_COMPOUND_AXIS = (
    "REFUSED: more than one compound reached a single axis of this view. "
    "Ordering two different compounds by predicted effect was WITHDRAWN on "
    "2026-08-21, absolutely — no flag, no mode, no caveat and no axis label "
    "restores it, because one grid on one colour scale IS a ranking however "
    "the compounds got there. The two measurements that withdrew it: "
    "results/pk_ranking_robustness.md — per-compound plasma exposure spans "
    "113.1x across the compounds the GBM panel ranks while the published "
    "ordering first breaks at a spread of 1.53x; "
    "results/capacity_ranking_robustness.md — the breast top-1 pick changes "
    "under this repo's own structural capacity rule. The same rule is enforced "
    "at two other layers and this is the third: "
    "app/server/projects.py::_kill_comparability raises on a second compound, "
    "and tmesim/protocol/executor.py::_check_one_compound plans zero runs for a "
    "two-compound program. This payload was BUILT and then REFUSED by its own "
    "audit, so the defect is in the builder, not in the request."
)

#: Every quantity on this axis is read off the delivery map, and on a spot
#: section that map is invalid — stated by ``app/server/api.py::
#: SPOT_SPATIAL_REFUSAL`` in terms, which is imported and carried verbatim
#: rather than summarised. Each entry is (quantity, the sentence of the refusal
#: that covers it).
SPOT_REFUSED_QUANTITIES: Tuple[Tuple[str, str], ...] = (
    ("reach ladder F(x)",
     "F(x) counts VIABLE TUMOUR CELLS whose own exposure falls below x times "
     "the perivascular reference — a per-cell statistic read off the delivery "
     "map, whose source is a density over endothelial cell POSITIONS that this "
     "section synthesises."),
    ("perivascular reference exposure E_ref",
     "E_ref is the first usable bin of the distance-from-nearest-endothelial-"
     "cell profile, and the refusal names that profile "
     "(validation.vessel_profile.direct_vessel_distance_um) explicitly."),
    ("normalised falloff (exposure vs distance from the nearest vessel)",
     "the falloff bins by distance to the nearest endothelial cell, so on a "
     "spot section it measures the 12-pseudo-cells-per-spot synthesis rather "
     "than this tissue's vascular architecture."),
    ("penetration_length_um",
     "the penetration length is FIT on that same distance profile, so it "
     "inherits the same defect; E_ref is reused from that fit verbatim "
     "(app/server/reach.py) so the two cannot be refused separately."),
)

#: Travels with every fraction on the axis. Two independent reasons, both
#: measured, and neither is a reason to withhold the number — they are reasons
#: it must never be read without its denominator beside it.
DENOMINATOR_CAVEAT = (
    "F(x) is a FRACTION OF THIS SECTION'S OWN VIABLE TUMOUR CELLS, and that "
    "denominator is not common across sections. Two things follow and both are "
    "measured. (1) A fraction and the count beneath it can move in OPPOSITE "
    "directions when the denominator changes: results/heldout_patient.md §1.1 "
    "measured exactly that on two real GBM sections — the endothelial FRACTION "
    "fell to 0.67x while the endothelial AREAL DENSITY rose 4.19x, because the "
    "second section held 6.2x the cell density — and its own conclusion is "
    "'these two rows point in opposite directions, and that is the finding'. "
    "So the count is shown beside every fraction here, and a direction "
    "disagreement is reported rather than smoothed over. (2) WHICH cells are "
    "'viable' is banded by the oxygen field, whose absolute level has never "
    "been calibrated: results/cross_tissue.md records that a ±40% change in "
    "the lumped consumption constant swings the real-section hypoxic tumour "
    "fraction between ~0.0 and ~0.6, and instructs that such fractions be read "
    "as 'a lot / a little' and never differenced. The denominators of two "
    "sections are therefore not on a common footing, which is the second "
    "reason the verdict below is ORDINAL rather than a difference of numbers."
)

#: Said once, at the point where a reader would otherwise supply it themselves.
NO_CROSS_SECTION_NOISE_FLOOR = (
    "NO CROSS-SECTION NOISE FLOOR HAS BEEN MEASURED FOR THIS LADDER. A crossing "
    "reported below is a statement about THESE runs, not a measured effect, and "
    "nothing here is a significance test. The nearest measurement this repo "
    "owns is results/delivery_spread_noise.md, which moved one compound's "
    "delivery spread by <= 10% across a 2.5x change in grid resolution, mask "
    "repair on/off and a 4x change in exposure window — and it was measured on "
    "ONE section with the tissue held fixed, so it bounds ANALYSIS choices on "
    "one tissue, not the difference between two tissues. That document's own "
    "'not established' list says it does not generalise to other tissues; it is "
    "cited here as the size of a nuisance span and for nothing else."
    " ONE cross-donor span does now exist, for ONE statistic and ONE compound, "
    "and it is not a noise floor either: across three single-cell GBM sections "
    "the delivery spread E_p90/E_p10 reads 1.9936 / 2.3211 / 2.5293 (1.2687x) "
    "while the fitted penetration length spans 2.691x and F(1.0) spans 119.9x "
    "(results/multisection_generalisation_2026-09-20.md sec. 4.2, 4.3). Three "
    "donors is not a distribution and none of those numbers licenses a "
    "significance test here."
)

#: Why one signature. ``results/cross_tissue.md``'s retraction, in one sentence.
ONE_SIGNATURE_ONLY = (
    "Only ONE RunSpec signature is ever placed on this axis. "
    "results/cross_tissue.md retracted its own published cross-tissue spread "
    "because the sections had been run at different grid resolutions, so the "
    "reported range 'mixed real biology with numerical error'; the grid is "
    "matched here by construction, and runs at any other signature are listed "
    "as available sets rather than silently merged in."
)

#: What the ladder verdict is, and is not.
ORDINAL_VERDICT_NOTE = (
    "The verdict is ORDINAL: it reads the SIGN of the difference between two "
    "sections at each rung of the ladder and reports whether that sign is "
    "constant. It deliberately does not difference the fractions, does not "
    "pick a rung and does not pick a threshold — x = 0.1 is a reporting "
    "convention nothing in the physics marks out (app/server/reach.py::"
    "CONVENTION_NOTE) and a headline chosen after seeing the curve is post-hoc "
    "selection. Read the curve, not any one point on it."
)


# --------------------------------------------------------------------------- #
# Spot sections — imported, never assumed
# --------------------------------------------------------------------------- #
def _spot_rules() -> Tuple[Sequence[str], str]:
    """``(_SPOT_SOURCES, SPOT_SPATIAL_REFUSAL)`` from ``app/server/api.py``.

    Imported INSIDE the function on purpose. ``api.py`` mounts this module's
    router from inside ``create_app``; a module-level import back into ``api``
    would make the two files' import order load-bearing for no gain, and this
    call happens once per request against an already-imported module.

    Never re-implemented: which platforms are spot platforms, and the sentence
    that refuses them, have exactly one home. If ``api.py`` adds a platform,
    this view stops publishing per-cell statistics for it on the same commit.
    """
    from .api import SPOT_SPATIAL_REFUSAL, _SPOT_SOURCES

    return _SPOT_SOURCES, SPOT_SPATIAL_REFUSAL


def _platform_of(record: Optional[Dict], capability: Optional[Dict]) -> str:
    """The platform string, read the way ``api.py::_with_caveats`` reads it."""
    for src in (record or {}, capability or {}):
        platform = str(src.get("platform") or src.get("source") or "").lower()
        if platform:
            return platform
    return ""


def _spot_verdict(record: Optional[Dict],
                  capability: Optional[Dict]) -> Dict:
    """Does this section's platform refuse per-cell statistics?

    ``spatial_spread_valid`` follows ``api.py``'s tri-state exactly: ``None``
    for a platform we could not name — "we did not establish this" is a
    different claim from "this is valid" — and a section we cannot name is NOT
    treated as refused, because inventing a refusal is as dishonest as
    inventing a number. It is reported as unknown and its numbers carry that.
    """
    spot_sources, refusal = _spot_rules()
    platform = _platform_of(record, capability)
    is_spot = platform in spot_sources
    return {
        "platform": platform or None,
        "per_cell_statistics_refused": bool(is_spot),
        "spatial_spread_valid": (None if not platform else not is_spot),
        "refusal": refusal if is_spot else None,
        "refused_quantities": (
            [{"quantity": q, "why": why} for q, why in SPOT_REFUSED_QUANTITIES]
            if is_spot else []),
    }


# --------------------------------------------------------------------------- #
# The guard
# --------------------------------------------------------------------------- #
#: Keys anywhere in a payload that name a compound. Walked generically so a
#: field added later cannot smuggle a second compound past this audit by being
#: new.
_COMPOUND_KEYS = ("compound", "drug", "compounds", "drugs")


def _compounds_in(node, found: set) -> set:
    """Every compound named anywhere beneath ``node``, however deep."""
    if isinstance(node, dict):
        for key, value in node.items():
            if key in _COMPOUND_KEYS:
                if isinstance(value, str) and value.strip():
                    found.add(value.strip().lower())
                elif isinstance(value, (list, tuple)):
                    for item in value:
                        if isinstance(item, str) and item.strip():
                            found.add(item.strip().lower())
            _compounds_in(value, found)
    elif isinstance(node, (list, tuple)):
        for item in node:
            _compounds_in(item, found)
    return found


def assert_one_compound_axis(payload: Dict) -> str:
    """Audit a FINISHED payload: exactly one compound is named anywhere in it.

    This is the mechanical guard, and it is mechanical in two senses. It walks
    the payload RECURSIVELY rather than checking the fields this file happens
    to write today, so a field added next month is audited the day it is added.
    And it runs on the OUTPUT rather than the input, because the input filter —
    "select runs whose drug is X" — is precisely the line a future edit would
    relax, and an audit that trusts the filter it is guarding guards nothing.

    Returns the single compound. Raises :class:`CrossSectionError` with
    :data:`WITHDRAWN_CROSS_COMPOUND_AXIS` on anything else, as a 500: a caller
    cannot cause this, so it is this server's defect, not their request.
    """
    found = sorted(_compounds_in(payload, set()))
    if len(found) == 1:
        return found[0]
    raise CrossSectionError(
        f"{WITHDRAWN_CROSS_COMPOUND_AXIS} Compounds found on the axis: "
        f"{found or '(none)'}.",
        status=500,
        detail={"compounds": found, "axis": ONE_COMPOUND_AXIS},
    )


# --------------------------------------------------------------------------- #
# Selecting the runs
# --------------------------------------------------------------------------- #
def _normalise_compound(raw) -> str:
    """One compound name, or a refusal.

    A comma- or space-separated list is refused with the withdrawal text rather
    than with "unknown compound": a caller who typed two names is asking for the
    output that was withdrawn, and telling them their compound does not exist
    would be a lie about the reason.
    """
    if not isinstance(raw, str) or not raw.strip():
        raise CrossSectionError(
            "compound: required. This view holds ONE compound fixed and varies "
            "the section, so it cannot be built without being told which "
            "compound. GET /api/cross-section/compounds lists the compounds "
            "this server has run on two or more sections, as a SET.",
            status=400)
    name = raw.strip().lower()
    parts = [p for p in name.replace(",", " ").split() if p]
    if len(parts) > 1:
        raise CrossSectionError(
            f"{WITHDRAWN_CROSS_COMPOUND_AXIS} Requested: {parts}. Ask for one "
            "of them; there is no mode of this route that answers for two.",
            status=400,
            detail={"compounds": sorted(set(parts))})
    return parts[0]


def _run_views(jobs) -> List[Dict]:
    """Every finished run this server holds, as ``RunResult`` views."""
    out: List[Dict] = []
    for job in jobs.list():
        if getattr(job, "status", None) != "done":
            continue
        try:
            out.append(jsonable(job.result_view()))
        except Exception:  # noqa: BLE001 - one unreadable run must not hide the rest
            continue
    return out


def available_compounds(jobs) -> Dict:
    """Which compounds this server can build a cross-section view for.

    A SET, alphabetical, never ordered by any metric and never carrying a
    "best" mark (``docs/WORKSPACE_CONTRACT.md`` §6). The only number beside a
    compound is how many SECTIONS it has been run on, which is a count of
    tissue and not a statement about the compound.
    """
    per_compound: Dict[str, Dict[str, set]] = {}
    for run in _run_views(jobs):
        drug = _drug_of(run)
        section = (run.get("spec") or {}).get("section")
        if not drug or not section:
            continue
        slot = per_compound.setdefault(
            drug, {"sections": set(), "signatures": set(), "runs": set()})
        slot["sections"].add(str(section))
        slot["signatures"].add(spec_signature(run.get("spec") or {}))
        slot["runs"].add(run.get("run_id"))
    rows = [
        {
            "compound": name,
            "n_sections": len(slot["sections"]),
            "n_runs": len(slot["runs"]),
            "n_signatures": len(slot["signatures"]),
            "sections": sorted(slot["sections"]),
            # A view needs two sections to be a comparison at all. One section
            # is still listed, with the reason, rather than hidden: "not
            # offered" and "not present" are different facts.
            "comparable_sections": len(slot["sections"]) >= 2,
            "reason": (None if len(slot["sections"]) >= 2 else
                       "run on one section only, and one section is not a "
                       "cross-section comparison"),
        }
        for name, slot in sorted(per_compound.items())
    ]
    return {
        "compounds": rows,
        "ordering": (
            "ALPHABETICAL, and deliberately so. This list is a SET of what has "
            "been run. It is never ordered by any metric and carries no 'best' "
            "mark: " + WITHDRAWN_CROSS_COMPOUND_AXIS),
        "axis": ONE_COMPOUND_AXIS,
        "generated_at": time.time(),
    }


# --------------------------------------------------------------------------- #
# The reach ladder across sections
# --------------------------------------------------------------------------- #
def _reach_block(run: Dict) -> Optional[Dict]:
    block = run.get("reach")
    return block if isinstance(block, dict) else None


def _ladder_of(reach: Optional[Dict]) -> Optional[List[float]]:
    curve = (reach or {}).get("curve") or {}
    xs = curve.get("x")
    if not isinstance(xs, list) or not xs:
        return None
    try:
        return [float(x) for x in xs]
    except (TypeError, ValueError):
        return None


def _fractions_of(reach: Optional[Dict]) -> Optional[List[Optional[float]]]:
    curve = (reach or {}).get("curve") or {}
    values = curve.get("under_reached_fraction")
    if not isinstance(values, list) or not values:
        return None
    out: List[Optional[float]] = []
    for v in values:
        out.append(float(v) if isinstance(v, (int, float))
                   and not isinstance(v, bool) else None)
    return out


def _section_reach_state(section: Dict) -> Dict:
    """Can this section put a number on the ladder, and if not, why not.

    Four refusals, in the order a reader needs them: the platform refusal
    first (it is a statement about the tissue, not about this run), then the
    run's own null_reason, then a run that predates the reach map, then a
    malformed curve. Nothing is defaulted and nothing is interpolated.
    """
    if section["spot"]["per_cell_statistics_refused"]:
        return {
            "usable": False,
            "reason": section["spot"]["refusal"],
            "kind": "platform",
        }
    reach = section.get("reach")
    if reach is None:
        return {
            "usable": False,
            "kind": "absent",
            "reason": (
                "this run carries no reach map. It predates the map, or it was "
                "a reconstruction-only run with no compound to take a ratio "
                "for. Re-run this RunSpec to produce one; nothing is inferred "
                "for it here."),
        }
    if reach.get("null_reason"):
        return {"usable": False, "kind": "null", "reason": str(reach["null_reason"])}
    if _ladder_of(reach) is None or _fractions_of(reach) is None:
        return {
            "usable": False,
            "kind": "malformed",
            "reason": ("this run's reach map carries no readable F(x) curve, so "
                       "it contributes no rung. It is not treated as zero."),
        }
    return {"usable": True, "kind": None, "reason": None}


def _ladder_rows(usable: List[Dict]) -> Tuple[List[Dict], Optional[List[float]]]:
    """One row per rung: each section's fraction and the count beneath it.

    The ladder is each RUN'S OWN ``curve.x``. Two sections whose ladders differ
    are refused by the caller rather than interpolated onto a common grid —
    interpolating would manufacture a number at an ``x`` no run reported.
    """
    if not usable:
        return [], None
    ladder = _ladder_of(usable[0]["reach"])
    rows: List[Dict] = []
    for i, x in enumerate(ladder or []):
        values: Dict[str, Optional[float]] = {}
        counts: Dict[str, Optional[int]] = {}
        denominators: Dict[str, Optional[int]] = {}
        for section in usable:
            fractions = _fractions_of(section["reach"]) or []
            frac = fractions[i] if i < len(fractions) else None
            n_viable = section.get("n_viable")
            values[section["section"]] = frac
            denominators[section["section"]] = n_viable
            counts[section["section"]] = (
                int(round(frac * n_viable))
                if frac is not None and isinstance(n_viable, int) else None)
        rows.append({
            "x": float(x),
            "is_convention_point": bool(abs(float(x) - float(HEADLINE_X)) < 1e-12),
            "under_reached_fraction": values,
            # DERIVED, and labelled: fraction * denominator, rounded. It is
            # published because a fraction without the count beneath it is the
            # exact shape of results/heldout_patient.md §1.1's defect — not
            # because the engine reported a count.
            "under_reached_cells": counts,
            "of_viable_tumor_cells": denominators,
            "cells_are_derived": (
                "count = round(F(x) * viable tumour cells in that section). "
                "Derived here from the two published numbers, not reported by "
                "the run."),
        })
    return rows, ladder


def _pairwise_stability(rows: List[Dict], usable: List[Dict]) -> Dict:
    """Where the answer is stable across tissue, and where it is not."""
    ids = [s["section"] for s in usable]
    pairs: List[Dict] = []
    for i, a in enumerate(ids):
        for b in ids[i + 1:]:
            signs: List[Dict] = []
            count_signs: List[Dict] = []
            for row in rows:
                fa = row["under_reached_fraction"].get(a)
                fb = row["under_reached_fraction"].get(b)
                if fa is None or fb is None:
                    continue
                sign = 0 if fa == fb else (1 if fa > fb else -1)
                signs.append({"x": row["x"], "sign": sign,
                              "values": {a: fa, b: fb}})
                ca = row["under_reached_cells"].get(a)
                cb = row["under_reached_cells"].get(b)
                if ca is not None and cb is not None:
                    count_signs.append({
                        "x": row["x"],
                        "sign": 0 if ca == cb else (1 if ca > cb else -1),
                        "counts": {a: ca, b: cb},
                    })
            nonzero = {s["sign"] for s in signs if s["sign"] != 0}
            flips = len(nonzero) > 1
            # Where the FRACTION and the COUNT disagree about which section is
            # worse — results/heldout_patient.md §1.1's pattern, at a rung.
            by_x = {s["x"]: s for s in count_signs}
            divergences = [
                {"x": s["x"], "fraction_sign": s["sign"],
                 "count_sign": by_x[s["x"]]["sign"],
                 "fractions": s["values"], "counts": by_x[s["x"]]["counts"]}
                for s in signs
                if s["x"] in by_x and s["sign"] != 0
                and by_x[s["x"]]["sign"] != 0
                and s["sign"] != by_x[s["x"]]["sign"]
            ]
            pairs.append({
                "sections": [a, b],
                "n_rungs_compared": len(signs),
                "order_flips": flips,
                "flip_rungs": [s["x"] for s in signs if s["sign"] != 0] if flips else [],
                "per_rung": signs,
                "higher_throughout": (
                    None if flips or not nonzero
                    else (a if nonzero == {1} else b)),
                "fraction_count_divergences": divergences,
                "verdict": (
                    "ORDER FLIPS ALONG THE LADDER — these two tissues cannot be "
                    "ordered for this compound at all"
                    if flips else
                    "ORDER HOLDS AT EVERY RUNG COMPARED"
                    if nonzero else
                    "IDENTICAL AT EVERY RUNG COMPARED"),
            })
    any_flip = any(p["order_flips"] for p in pairs)
    any_divergence = any(p["fraction_count_divergences"] for p in pairs)
    if len(usable) < 2:
        verdict, level = "NOT A COMPARISON", "info"
        text = (
            f"{len(usable)} section can put a number on this ladder, and a "
            "comparison needs two. Section 1 above says why each of the others "
            "cannot.")
    elif any_flip:
        verdict, level = "THE ANSWER IS NOT STABLE ACROSS TISSUE", "crit"
        text = (
            "At least one pair of sections CROSSES: which tissue is more "
            "under-reached depends on which rung of the ladder you read. There "
            "is no ordering of those two tissues for this compound. This is the "
            "most informative thing this view says — it is the disagreement, "
            "not a failure of the run.")
    else:
        verdict, level = "THE ORDER HOLDS AT EVERY RUNG COMPARED", "ok"
        text = (
            "Every pair of sections keeps the same sign at every rung both "
            "could supply: the curves do not cross. That is an ORDINAL "
            "statement about these tissues for this compound, at this "
            "signature, and it is not a measured effect size.")
    return {
        "verdict": verdict,
        "level": level,
        "text": text,
        "pairs": pairs,
        "any_order_flips": any_flip,
        "any_fraction_count_divergence": any_divergence,
        "how_read": ORDINAL_VERDICT_NOTE,
        "noise_floor": NO_CROSS_SECTION_NOISE_FLOOR,
    }


def _reach_comparability(usable: List[Dict], refused: List[Dict]) -> Dict:
    """Blockers and caveats for the ladder itself.

    Separate from :func:`app.server.projects._kill_comparability`, which is
    about an ABSOLUTE level. This one is about a within-tissue RATIO, so it
    blocks on different things — but it is written in the same shape, and the
    absolute block below still defers to that function rather than to this one.
    """
    blockers: List[str] = []
    caveats: List[str] = []

    for section in refused:
        blockers.append(
            f"{section['label'] or section['section']} ({section['section']}): "
            f"{section['reach_state']['reason']}")

    modes = {s["window_mode"] for s in usable if s.get("window_mode")}
    if len(modes) > 1:
        blockers.append(
            f"the sections' reach maps are not the same OBJECT: window_mode "
            f"{sorted(modes)} appear together. A steady-state map is the reach "
            "this tissue's geometry permits and a finite-window map is the "
            "reach one dosing window achieved; they answer different questions "
            "and are not placed on one axis (app/server/reach.py::"
            "STEADY_STATE_NOTE, ::FINITE_WINDOW_NOTE).")

    ladders = {tuple(_ladder_of(s["reach"]) or ()) for s in usable}
    if len(ladders) > 1:
        blockers.append(
            "the sections report different F(x) ladders, so no rung is shared. "
            "Nothing is interpolated onto a common grid here: that would "
            "manufacture a value at an x no run reported.")

    if len(usable) < 2 and not blockers:
        blockers.append(
            "fewer than two sections can supply this ladder, so there is "
            "nothing to compare.")

    panels = {s["section"]: (s.get("capability") or {}).get("panel")
              for s in usable}
    if len(set(panels.values())) > 1:
        caveats.append(
            f"the sections are on DIFFERENT panels ({panels}). The vessel map "
            "every reach ratio is referenced to is built from whichever "
            "endothelial transcripts the panel carries, and "
            "results/heldout_patient.md §1.2 measured exactly how much that "
            "matters: a 1-marker gate (PECAM1, 24.57% detection) against a "
            "3-marker mean-z gate (VWF 4.15%, FLT1 9.42%, PLVAP 4.36%) pushes "
            "the endothelial call down for a purely technical reason. Part of "
            "any difference below is a difference in what was measurable.")

    caveats.append(DENOMINATOR_CAVEAT)
    caveats.append(RELATIVE_EXPOSURE_FULL)
    for key in ("scale_invariance", "window_note", "convention"):
        for section in usable:
            text = (section.get("reach") or {}).get(key)
            if isinstance(text, str) and text and text not in caveats:
                caveats.append(text)
    caveats.append(
        "NOT FOR CLINICAL USE. No section here carries a drug-response outcome "
        "label, so no reach curve on this axis has ever been scored against a "
        "measured response, and this is not a statement about a patient.")

    return {"comparable": not blockers, "blockers": blockers,
            "caveats": caveats}


# --------------------------------------------------------------------------- #
# The absolute metrics — gated by the EXISTING backend verdict
# --------------------------------------------------------------------------- #
def _absolute_block(usable_and_refused: List[Dict], metric: str) -> Dict:
    """One row for ``metric``: this compound's own value on each section.

    The verdict is :func:`app.server.projects._kill_comparability`, unchanged
    and unwrapped — the same function ``GET /api/projects/{id}/compare`` refuses
    with. Where it says ``comparable: false`` this block carries NO values at
    all, only the refusal and its blockers: showing two absolute numbers side by
    side under a "not comparable" chip is still showing two absolute numbers
    side by side, and the instruction this view was built to is to render the
    refusal IN PLACE OF the number.
    """
    entries = [
        {
            # `_kill_comparability` reads `specimen.label`, `specimen.section`
            # and `specimen.id`. A section is not a specimen, but the rule it
            # applies is about the PANEL and the RUN, both of which are present
            # here — so the rule is reused rather than re-derived.
            "specimen": {"id": s["section"], "label": s["label"],
                         "section": s["section"]},
            "capability": s["capability"],
            "run": s["run"],
        }
        for s in usable_and_refused
    ]
    if len(entries) < 2:
        return {
            "metric": metric,
            "comparable": False,
            "blockers": ["fewer than two sections ran this compound at this "
                         "signature, so there is no cross-section value to "
                         "read."],
            "caveats": [],
            "values": {},
            "state_measured": {},
            "withheld": True,
        }
    verdict = _kill_comparability(entries)
    values = {}
    for entry in entries:
        raw = (entry["run"].get("metrics") or {}).get(metric)
        values[entry["specimen"]["id"]] = (
            float(raw) if isinstance(raw, (int, float))
            and not isinstance(raw, bool) else None)

    spot_withheld = {}
    if metric == "penetration_length_um":
        # The one headline metric that is FIT on the vessel-distance profile,
        # which SPOT_SPATIAL_REFUSAL names explicitly. A spot section does not
        # contribute it even when the comparability verdict is otherwise clear.
        for section in usable_and_refused:
            if section["spot"]["per_cell_statistics_refused"]:
                spot_withheld[section["section"]] = section["spot"]["refusal"]
                values[section["section"]] = None

    comparable = bool(verdict["comparable"]) and not spot_withheld
    return {
        "metric": metric,
        "comparable": comparable,
        # NEVER a headline: absolute kill is not one anywhere in this product,
        # and this block sits below the ladder for that reason.
        "is_headline": False,
        "blockers": list(verdict["blockers"]) + (
            [f"{sid}: {why}" for sid, why in sorted(spot_withheld.items())]),
        "caveats": list(verdict["caveats"]),
        "state_measured": verdict["state_measured"],
        # Withheld wholesale when not comparable. The refusal takes the place
        # of the numbers rather than sitting beside them.
        "values": values if comparable else {},
        "withheld": not comparable,
        "withheld_reason": (
            None if comparable else
            "the backend verdict on this pair is comparable=false, so the "
            "values are not rendered at all. A number shown under a refusal is "
            "still a number shown."),
    }


# --------------------------------------------------------------------------- #
# The payload
# --------------------------------------------------------------------------- #
def cross_section_payload(
    jobs,
    *,
    compound: str,
    signature: Optional[str] = None,
    metric: str = "kill_fraction",
    capability: Callable[[str], Dict] = capability_report,
    section_record: Callable[[str], Optional[Dict]] = lookup_section,
) -> Dict:
    """ONE compound across every section that ran it at one RunSpec signature.

    ``capability`` and ``section_record`` are injected for the same reason
    ``compare_project`` takes ``resolve_run``: the tests build an axis out of
    sections that do not have to exist on this checkout's disk.
    """
    if metric not in _HEADLINE_METRICS:
        raise CrossSectionError(
            f"metric: {metric!r} is not one of {list(_HEADLINE_METRICS)}. It "
            "selects which absolute column is shown BELOW the ladder and can "
            "never change the ladder, the verdict or which sections appear.",
            status=400)
    name = _normalise_compound(compound)

    runs = [r for r in _run_views(jobs) if _drug_of(r) == name]
    problems: List[str] = []

    # --- group by signature: one experiment per axis ----------------------- #
    groups: Dict[str, Dict[str, Dict]] = {}
    for run in runs:
        section = (run.get("spec") or {}).get("section")
        if not section:
            problems.append(
                f"run {run.get('run_id')} names no section and cannot be placed "
                "on an axis of sections")
            continue
        sig = spec_signature(run.get("spec") or {})
        slot = groups.setdefault(sig, {})
        previous = slot.get(str(section))
        # Newest wins, and the earlier one is reported rather than dropped in
        # silence — the same rule compare_project applies.
        if previous is None:
            slot[str(section)] = run
            continue
        # More than one run of this compound on this section at this
        # signature. Newest wins, and the loser is REPORTED either way —
        # `jobs.list()` is newest-first, so the run being displaced is usually
        # the arriving one, and recording the problem only in the replacing
        # branch would have kept that (commoner) case silent.
        newer, older = ((run, previous)
                        if (run.get("created_at") or 0) > (previous.get("created_at") or 0)
                        else (previous, run))
        slot[str(section)] = newer
        problems.append(
            f"section {section} has more than one run of {name} at signature "
            f"{sig}; the newest ({newer.get('run_id')}) is on the axis and "
            f"{older.get('run_id')} is not")

    sets = sorted(
        ({"signature": sig, "n_sections": len(slot),
          "sections": sorted(slot), "runs": {s: r.get("run_id")
                                             for s, r in sorted(slot.items())}}
         for sig, slot in groups.items()),
        key=lambda s: (s["n_sections"], s["signature"]), reverse=True)

    chosen = signature or (sets[0]["signature"] if sets else None)
    if signature is not None and signature not in groups:
        raise CrossSectionError(
            f"signature: {signature!r} has no runs of {name!r} on this server. "
            f"Signatures available for it: {[s['signature'] for s in sets]}",
            status=404)
    for entry in sets:
        entry["is_primary"] = entry["signature"] == chosen

    selected = groups.get(chosen or "", {})

    # --- one row per section ---------------------------------------------- #
    sections: List[Dict] = []
    for section_id, run in sorted(selected.items()):
        cap: Dict = {}
        cap_error: Optional[str] = None
        try:
            cap = capability(section_id)
        except ProjectError as exc:
            # A DELIBERATE refusal — most often a run whose section is no longer
            # registered on this server. `cache.redact_exception`'s own
            # docstring exempts the capability report from redaction because it
            # is written for a reader, so it is carried verbatim. The section
            # KEEPS ITS ROW: a run that cannot be described is still a run that
            # exists, and dropping it would quietly shrink the axis.
            cap_error = str(exc)
            problems.append(f"section {section_id}: {cap_error}")
            # A panel we cannot read is a panel that has NOT cleared the
            # absolute-kill bar, and the reason is the refusal itself rather
            # than a None that reads as "no reason given".
            cap = {"label": section_id, "available": False,
                   "absolute_kill_supported": False,
                   "absolute_kill_reason": cap_error}
        except Exception as exc:  # noqa: BLE001 - a section must still be listed
            cap_error = redact_exception(
                exc, f"capability for section {section_id!r} could not be read",
                brief=True)
            problems.append(cap_error)
            cap = {"label": section_id, "available": False,
                   "absolute_kill_supported": False,
                   "absolute_kill_reason": cap_error}
        try:
            record = section_record(section_id)
        except Exception:  # noqa: BLE001
            record = None
        reach = _reach_block(run)
        spot = _spot_verdict(record, cap)
        labels: List[str] = []
        if not cap.get("absolute_kill_supported"):
            labels.append("ABSOLUTE KILL NOT SUPPORTED BY THIS PANEL")
        unsupported = sorted(
            n for n, blk in (cap.get("state_routes") or {}).items()
            if blk.get("supported") is False)
        if unsupported:
            labels.append("UNSUPPORTED STATE ROUTE(S): " + ", ".join(unsupported))
        if spot["per_cell_statistics_refused"]:
            labels.append("PER-CELL STATISTICS REFUSED ON THIS PLATFORM")
        if not cap.get("available", True):
            labels.append("SECTION NOT LOADABLE IN THIS CHECKOUT")
        if cap_error:
            labels.append("CAPABILITY UNKNOWN: " + cap_error)
        row = {
            "section": section_id,
            "label": cap.get("label") or (record or {}).get("label") or section_id,
            # THE RUN'S OWN compound, read back off its RunSpec — never the
            # compound that was REQUESTED. This is what makes the audit bite:
            # if the selection above ever admitted a run of another compound,
            # this row would carry that compound's name and
            # `assert_one_compound_axis` would refuse the whole payload. A row
            # stamped with the requested name would agree with itself no matter
            # which run produced it, which is a guard that cannot fail.
            "compound": _drug_of(run),
            # Echoed for the same reason, and shown in the view: the axis is
            # only one experiment if these agree.
            "spec": {k: v for k, v in (run.get("spec") or {}).items()
                     if k in ("section", "drug", "dx_um", "state_mode",
                              "sensitivity_mode", "vascular_mode",
                              "oxygen_mode", "exposure_window_s")},
            "run_id": run.get("run_id"),
            "created_at": run.get("created_at"),
            "run": run,
            "capability": cap,
            "panel": cap.get("panel"),
            "labels": labels,
            "unsupported_state_routes": unsupported,
            "absolute_kill_supported": bool(cap.get("absolute_kill_supported")),
            "absolute_kill_reason": cap.get("absolute_kill_reason"),
            "spot": spot,
            "reach": reach,
            "window_mode": (reach or {}).get("window_mode"),
            "n_viable": (reach or {}).get("n_viable"),
            "n_tumor": (reach or {}).get("n_tumor"),
            "E_ref": (None if spot["per_cell_statistics_refused"]
                      else (reach or {}).get("E_ref_uM_s")),
            "E_ref_units": (reach or {}).get("E_ref_units"),
        }
        row["reach_state"] = _section_reach_state(row)
        sections.append(row)

    # THE GUARD, first pass. The final audit below is the one that cannot be
    # walked past, but it runs after `_absolute_block`, and by then
    # `projects._kill_comparability` would already have raised its own
    # cross-compound refusal one layer down. Both refusals are correct; this
    # one is the one that names THIS view's axis, so it goes first. A
    # projection is audited rather than the raw rows: the run bodies hanging
    # off them carry compiled drug blocks, and an audit should read the axis,
    # not everything the axis was built from.
    assert_one_compound_axis({"sections": [
        {"compound": s["compound"], "spec": s["spec"]} for s in sections]}
        if sections else {"compound": name})

    usable = [s for s in sections if s["reach_state"]["usable"]]
    refused = [s for s in sections if not s["reach_state"]["usable"]]
    rows, ladder = _ladder_rows(usable)
    verdict = _reach_comparability(usable, refused)
    stability = _pairwise_stability(rows, usable) if verdict["comparable"] else {
        "verdict": "REFUSED",
        "level": "crit",
        "text": ("The ladder is not comparable across these sections, so no "
                 "stability verdict is computed. The blockers are the output."),
        "pairs": [],
        "any_order_flips": False,
        "any_fraction_count_divergence": False,
        "how_read": ORDINAL_VERDICT_NOTE,
        "noise_floor": NO_CROSS_SECTION_NOISE_FLOOR,
    }

    payload = {
        "compound": name,
        "axis": {
            "fixed": "compound",
            "varies": "section",
            "compound": name,
            "n_sections": len(sections),
            "statement": ONE_COMPOUND_AXIS,
        },
        "withdrawn": WITHDRAWN_CROSS_COMPOUND_AXIS,
        "signature": chosen,
        "signature_note": ONE_SIGNATURE_ONLY,
        "matched_sets": sets,
        "metric": metric,
        "metric_means": (
            "which absolute column is shown below the ladder. It selects a "
            f"column, not a conclusion, and every row carries {list(_HEADLINE_METRICS)} "
            "on the run behind it regardless."),
        "sections": [
            {k: v for k, v in s.items() if k not in ("run", "reach")}
            for s in sections
        ],
        "reach": {
            "ladder": ladder or list(REACH_LADDER),
            "ladder_is": ("each RUN'S OWN curve.x, not a ladder this view "
                          "chose. Sections whose ladders differ are refused, "
                          "never interpolated onto a common grid."),
            "rows": rows,
            "sections_on_axis": [s["section"] for s in usable],
            "sections_refused": [
                {"section": s["section"], "label": s["label"],
                 "kind": s["reach_state"]["kind"],
                 "reason": s["reach_state"]["reason"]}
                for s in refused
            ],
            "comparable": verdict["comparable"],
            "blockers": verdict["blockers"],
            "caveats": verdict["caveats"],
            "headline_x_is_a_convention": HEADLINE_X,
        },
        "stability": stability,
        "absolute": _absolute_block(sections, metric),
        "problems": problems,
        "generated_at": time.time(),
        "reading_order": [
            "sections — what each panel and platform can support, and the "
            "denominator each fraction below is a fraction OF",
            "reach — one compound's F(x) ladder on each section, refusals in "
            "place of the numbers they refuse",
            "stability — where the answer holds across tissue and where it "
            "flips. The flip is the output",
            "absolute — this compound's own level per section, withheld "
            "wholesale unless the backend verdict says comparable",
        ],
    }
    # THE GUARD. Last thing before the payload leaves: see
    # `assert_one_compound_axis` for why it audits the output.
    assert_one_compound_axis(payload)
    return payload


# --------------------------------------------------------------------------- #
# Router
# --------------------------------------------------------------------------- #
def build_cross_section_router(*, jobs) -> APIRouter:
    """Mount ``GET /api/cross-section`` and ``/api/cross-section/compounds``.

    Gated by the same middleware every other ``/api/`` route is: ``create_app``
    tests ``get_route_path(request.scope).startswith("/api/")``, so a router
    mounted there is authenticated by construction.
    """
    router = APIRouter()

    def _refuse(exc: CrossSectionError):
        raise HTTPException(status_code=exc.status,
                            detail={"error": exc.message, **exc.detail})

    @router.get("/api/cross-section/compounds")
    def get_cross_section_compounds() -> Dict:
        """The compounds this server has run, as a SET. Never ordered."""
        return jsonable(available_compounds(jobs))

    @router.get("/api/cross-section")
    def get_cross_section(
        # Optional at the HTTP layer so a caller who names no compound gets
        # this module's own 400 — which says what the view is and how to find
        # a compound — rather than FastAPI's 422 field error.
        compound: Optional[str] = Query(None,
                                        description="the ONE compound held fixed"),
        signature: Optional[str] = Query(None),
        metric: str = Query("kill_fraction"),
    ) -> Dict:
        """ONE compound, many sections — where the answer holds and where it
        flips. Two different compounds are refused at 400, and a payload that
        somehow carried two would be refused at 500 by its own audit."""
        try:
            return jsonable(cross_section_payload(
                jobs, compound=compound, signature=signature, metric=metric))
        except CrossSectionError as exc:
            _refuse(exc)
        except ProjectError as exc:
            # `_kill_comparability` raises this when it is handed more than one
            # compound. It is the SAME refusal as the audit's, one layer down,
            # and it is surfaced rather than becoming a bare 500.
            raise HTTPException(status_code=exc.status or 500,
                                detail={"error": str(exc.message)})

    return router
