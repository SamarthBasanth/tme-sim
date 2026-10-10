"""``POST /api/preview`` -- a COARSE SOLVE shown fast, not a model.

WHAT THIS IS
------------
The same solver the Studio queue runs, on a 2.5x coarser grid
(``dx = 30 um`` instead of the shipped ``dx = 12 um``), returning ONE readout:
the delivery spread ``E_p90/E_p10``. It is ~6.8x cheaper than the real run and
its disagreement with it has been MEASURED, pair by pair, on real sections.

WHY A COARSE SOLVE AND NOT A LEARNED SURROGATE
----------------------------------------------
``results/surrogate_feasibility_2026-10-09.md`` fitted ridge and gradient
boosting to 921 real solver runs and scored them leave-one-DONOR-out. On the
delivery spread the best model landed at MAE **0.4375** against an across-donor
range of **0.5357238964376265** -- 82 % of the readout's entire dynamic range.
The same solver on the coarser grid lands inside **0.0909 absolute**, which is
~4.8x better than the best model's MEAN error. So the thing a surrogate had to
beat was never "no preview"; it was this, and it lost.

The coarse solve also has three properties no fitted model has, and they matter
more than the error ratio:

* **Its output is the engine's.** It carries the engine's own caveats and the
  engine's own refusals, so it cannot disagree with the real run about what is
  defensible. A model can return a confident number on a section where the
  engine refuses to produce one at all.
* **It has no training set, so it cannot go stale.** The surrogate would be
  bound to one engine state (``alpha = 22.223477711172563``, the one-way
  delivery law, ``PDE_DT_S = 10.0``, convection ON) and would keep returning
  numbers after any of those moved. This repository's withdrawal rate is not
  hypothetical: the grid-convergence headline and the cross-compound ranking
  went in one week.
* **It makes no AI/ML claim**, so nothing here goes near the prohibition
  ``tests/test_home_page.py`` enforces on user-visible surfaces.

WHAT THIS REFUSES, AND WHY EACH REFUSAL IS A MEASUREMENT
--------------------------------------------------------
Every refusal below is a number someone took, not a policy someone preferred.
See :data:`REFUSED_READOUTS` and :data:`SPOT_PREVIEW_REFUSAL`.

1. **Spot (Visium) sections.** Not "because the statistics are weak there" --
   because ``dx_um`` IS NOT A LEVER THERE AT ALL.
   ``tmesim/reconstruct/real_pipeline.py::_reconstruct_visium`` reads
   ``visium_dx = float(dx_um) if dx_um and dx_um >= 50.0 else 100.0``, so a
   preview at 30 um and a real run at 12 um are BOTH clamped to 100 um and are
   the IDENTICAL SOLVE at the IDENTICAL COST. A "preview" there is a no-op that
   would report 0.000 % disagreement because it compared a run with itself --
   which is exactly how the feasibility study's first pass came to report a
   perfect result on 18 of 30 pairs (its sec. 9.1). Offering it would be the
   ``results/codon_vacuous_green_checks`` pattern with a progress bar on it.

2. **The fitted penetration length.** Measured disagreement up to **16.0405 %**
   / **17.3158 um** across the 12 genuine pairs -- the same order as the
   **19.4 %** this repository already records for that readout under a pure
   change of the pseudo-cell placement seed (``api.py::SPOT_SPATIAL_REFUSAL``).
   A 16 % error on a fitted length is not a preview of the answer, it is a
   different answer.

3. **The penetration index**, and this one is a DEVIATION from the recommendation
   in sec. 8 of the feasibility record, taken on evidence that record did not
   assemble. Its ``dx = 30`` vs ``dx = 12`` behaviour is UNRECONCILED across
   three of this repository's own measurements of nominally the same Xenium
   section, in two directions and over two orders of magnitude:

   ===========================================  ==========  ==========  =========
   source                                       dx = 30     dx = 12     apart
   ===========================================  ==========  ==========  =========
   ``results/numerical_rigor.md`` sec. 2        0.2247      0.3035      26.0 %
   ``results/cross_tissue.md`` ("Reading")      0.396       0.196       102 %
   ``experiments/.../coarse_baseline.json``     0.310702    0.313970    1.04 %
   ===========================================  ==========  ==========  =========

   (The third row is doxorubicin, the compound whose ~45 um penetration length
   ``numerical_rigor.md`` sec. 2 names as the reason the coarse grid cannot
   resolve the boundary layer. The first two rows rise and fall respectively.)

   **And the disagreement sits in the rung a preview is made of.** Row 3
   against row 1, same section, same compound, same two grids: 3.4 % apart at
   ``dx = 12`` (0.313970 against 0.3035) and 38.3 % apart at ``dx = 30``
   (0.310702 against 0.2247). The records agree where the product solves and
   diverge by nearly 1.4x where the preview would stand, so this is not a
   disagreement off to one side of the feature -- it is underneath it.

   ``numerical_rigor.md``'s own CORRECTION item 4 states the consequence:
   **"Neither ladder should be quoted as the solver's order of accuracy until
   that is settled."** There is therefore NO measured agreement figure to put
   beside a preview of this readout, and a preview without one is just a wrong
   number shown quickly. The shipped per-run RESOLUTION caveat
   (``jobs.py``) independently names ``dx`` as "the LARGEST single driver of the
   penetration index measured on that section (|elasticity| 0.747, 3.1x the
   next knob)" -- previewing a readout across a change in its own largest
   driver needs a settled figure, and this one has three.

4. **Any absolute level.** ``median_auc_over_peak`` agrees to 0.3215 % at worst
   -- the best agreement of any readout measured -- and is still refused,
   because agreement is not defensibility. It is an ABSOLUTE exposure level and
   ``jobs.py`` holds that tissue-over-plasma above 1 "IS THE LAW'S KNOWN
   UNBOUNDEDNESS MADE VISIBLE". Showing the product's least defensible number
   faster does not make it quotable. Anything whose name contains ``kill`` is
   refused on the same principle and is additionally never a defensible output.

5. **A pair nobody measured.** The agreement is for ``dx = 30`` against
   ``dx = 12`` on two Xenium sections and six compounds -- 12 pairs. Outside
   that envelope this module REFUSES rather than widening the tolerance. A
   preview's error bar is a measurement, and the way to grow the envelope is to
   run more pairs into ``coarse_baseline.json``, which this module reads at
   request time. It does not carry a second copy of those figures: if the data
   file is absent, THE FEATURE IS ABSENT, because the number and the statement
   of how far it may be from the truth ship together or not at all.

A PREVIEW IS NOT A RUN
----------------------
This module calls :func:`app.server.jobs.execute_run` directly and **discards
``_payloads``**. Nothing is enqueued, nothing is persisted, no ``run_id`` is
minted. So a preview number cannot be opened in the Read rail, attached to a
Specimen, exported, put in a report, or counted on a project card -- the four
routes by which a fast number becomes a quoted one. The payload says
``quotable: false`` and :func:`_audit_payload` refuses to return a body that
has lost either that flag or its agreement block.
"""
from __future__ import annotations

import json
import math
import os
import statistics as st
import time
from typing import Dict, List, Optional, Sequence, Tuple

from fastapi import APIRouter, HTTPException, Request

__all__ = [
    "AGREEMENT_SOURCE",
    "MEASURED_COMPOUNDS",
    "MEASURED_SECTIONS",
    "PREVIEW_DX_UM",
    "PREVIEW_LABEL",
    "PREVIEWABLE_READOUT",
    "REFUSED_READOUTS",
    "SPOT_PREVIEW_REFUSAL",
    "TARGET_DX_UM",
    "PreviewError",
    "build_preview_router",
    "measured_agreement",
    "preview_payload",
]

# --------------------------------------------------------------------------- #
# The one pair that was measured
# --------------------------------------------------------------------------- #
#: The grid the preview solves on.
PREVIEW_DX_UM: float = 30.0

#: The grid the REAL run solves on, and the only one the agreement covers.
#: NOT a claim that 12 um is converged -- ``results/CLAIM_INVENTORY.md`` sec. 1f
#: RETRACTS that claim, the GBM ladder is non-monotone with observed orders 0.96
#: and -3.33, and "neither panel is in an asymptotic regime". There is no
#: converged answer for either grid to approximate, and the preview does not
#: need one: the question is NOT "is this the true value" but "does the number
#: shown now match the number this product will print when the real run
#: finishes".
TARGET_DX_UM: float = 12.0

#: The sections the 12 genuine pairs were measured on. Both Xenium, and both
#: are sections where ``dx_um`` actually changes the grid
#: (170x271 -> 424x675 and 184x252 -> 458x628).
MEASURED_SECTIONS: Tuple[str, ...] = ("xenium_319", "xenium_breast")

#: The compounds those pairs covered. A SET, alphabetical, never ordered --
#: cross-compound ranking was withdrawn 2026-08-21 and nothing here restores
#: it. The compound is an input to the solve and never a thing being compared.
MEASURED_COMPOUNDS: Tuple[str, ...] = (
    "cisplatin", "doxorubicin", "lomustine", "methotrexate",
    "temozolomide", "vincristine",
)

#: Where the agreement figures come from, read at request time rather than
#: copied into this file. ``results/`` is this repository's claim surface and
#: has no CI; a second copy of a measured number in source is the defect this
#: project keeps finding (``results/codon_claim_surface_drift``).
AGREEMENT_SOURCE = os.path.join(
    "experiments", "surrogate_feasibility", "coarse_baseline.json")

#: The ONE readout a coarse solve may preview. ``app/server/reach.py``'s
#: ``CONVENTION_NOTE`` names the delivery spread (P90/P10 of tumour-cell
#: exposure) as the headline in place of ``F(0.1)``, and it is a RATIO of two
#: percentiles of one field -- no dose, no potency, no threshold to choose.
PREVIEWABLE_READOUT: str = "E_p90_over_p10"

#: Shown beside the number, every time, at no smaller than the number's own
#: size. The UI does not compose this sentence; it renders this one.
PREVIEW_LABEL: str = (
    "PREVIEW - A COARSER SOLVE, NOT THIS RUN'S RESULT. This is the same solver "
    "on a 2.5x coarser grid (dx = 30 um against this run's dx = 12 um), shown "
    "because it costs about a sixth as much. The real run is still going and "
    "its number is the one this product stands behind. DO NOT QUOTE THIS "
    "NUMBER: it is not in any result.json, it has no run_id, and it cannot be "
    "exported or reported."
)

#: Readout -> why a coarse solve may not preview it. Keyed by the name a caller
#: would ask for; every one of these is a number, not a preference.
REFUSED_READOUTS: Dict[str, str] = {
    "penetration_length_um": (
        "THE FITTED PENETRATION LENGTH IS NOT PREVIEWABLE. Measured across the "
        "12 genuine dx = 30 / dx = 12 pairs it disagrees by up to 16.0405 % "
        "(17.3158 um), against a median of 1.8053 %. That is the same order as "
        "the 19.4 % this repository already records for this readout under a "
        "pure change of the pseudo-cell placement seed, with nothing about the "
        "tissue changed (app/server/api.py::SPOT_SPATIAL_REFUSAL). A 16 % "
        "error on a fitted length is not a preview of the answer, it is a "
        "different answer. Wait for the run; it is the only number for this."
    ),
    "penetration_index": (
        "THE PENETRATION INDEX HAS NO SETTLED dx = 30 / dx = 12 AGREEMENT, so "
        "there is no figure to show beside a preview of it. Three measurements "
        "in this repository on nominally the same Xenium section disagree in "
        "two directions: results/numerical_rigor.md sec. 2 has it RISING "
        "0.2247 -> 0.3035 (26.0 % apart), results/cross_tissue.md has it "
        "FALLING 0.396 -> 0.196 (102 % apart), and "
        "experiments/surrogate_feasibility/coarse_baseline.json has it FLAT at "
        "0.310702 -> 0.313970 (1.04 % apart). AND THE DISAGREEMENT SITS IN THE "
        "RUNG A PREVIEW IS MADE OF: on the same section, the same compound and "
        "the same two grids, the stored pairs agree with that ladder to 3.4 % "
        "at dx = 12 um (0.313970 against 0.3035) and differ by 38.3 % at "
        "dx = 30 um (0.310702 against 0.2247). The coarse rung is the preview, "
        "so this is not a disagreement off to one side of the feature - it is "
        "underneath it. numerical_rigor.md's own CORRECTION item 4 states the "
        "consequence: neither ladder should be quoted until that is settled. "
        "The shipped RESOLUTION caveat on every run independently names dx as "
        "the LARGEST single driver of this index (|elasticity| 0.747, 3.1x the "
        "next knob). Settle the coarse rung first; a preview cannot be built "
        "on a disagreement this wide in exactly the place it would stand."
    ),
    "median_auc_over_peak": (
        "AN ABSOLUTE EXPOSURE LEVEL IS NOT PREVIEWED HERE, AND THE REASON IS "
        "NOT ITS AGREEMENT. It agrees best of anything measured - 0.3215 % at "
        "worst across the 12 pairs - and is still refused, because agreement "
        "is not defensibility. app/server/jobs.py holds that a tissue median "
        "over plasma peak above 1 IS THE LAW'S KNOWN UNBOUNDEDNESS MADE "
        "VISIBLE, never the tissue concentrating the drug, and "
        "results/CLAIM_INVENTORY.md sec. 2 lists any absolute number among "
        "what a reviewer would not accept. Showing the product's least "
        "defensible number faster does not make it quotable. WHAT TO DO: read "
        "the delivery spread instead - it is a RATIO between two points of one "
        "map, it is what this product headlines, and it is the one readout "
        "this route previews. If you need the absolute level itself, run the "
        "real solve and read it there with the caveats the run carries; it is "
        "not more defensible at 12 um than at 30 um, and this route declines "
        "to imply that it is."
    ),
    "reach_headline_F": (
        "THE REACH LADDER IS NOT PREVIEWED, AND IT DOES NOT NEED TO BE. "
        "app/server/reach.py computes it at steady state, so a self-normalised "
        "ratio there cannot depend on the schedule: measured fold change "
        "exactly 1.0 across a 4x dose, a 16x half-life and a 6x window, and "
        "bit-identical rungs over 4 checked pairs. It is also 27.5 % exactly "
        "0.0 on the 825 finite values collected, and refused outright on one of "
        "the four sections where per-cell statistics are defensible. "
        "reach.py::CONVENTION_NOTE is explicit that x = 0.1 is NO LONGER "
        "HEADLINED ANYWHERE. Read the curve on the finished run, not a point on "
        "a preview of it."
    ),
}

#: The spot refusal, stated in terms of the LINE OF CODE that makes it true
#: rather than in terms of the statistics being weak. Both are the case; only
#: this one makes the preview a NO-OP rather than an approximation.
SPOT_PREVIEW_REFUSAL: str = (
    "NO PREVIEW ON A SPOT SECTION, BECAUSE dx_um IS NOT A LEVER HERE. "
    "tmesim/reconstruct/real_pipeline.py::_reconstruct_visium reads "
    "`visium_dx = float(dx_um) if dx_um and dx_um >= 50.0 else 100.0`, so a "
    "30 um preview and a 12 um run are BOTH clamped to 100 um: the identical "
    "grid, the identical solve, the identical cost. There is no faster answer "
    "to show. A preview offered here would report 0.000 % disagreement with "
    "the real run because it had compared a run with ITSELF - which is how the "
    "study behind this feature first came to report a perfect result on 18 of "
    "its 30 pairs, and is the results/codon_vacuous_green_checks pattern "
    "exactly. Separately and independently, the delivery spread is not a "
    "reading of this tissue at all here "
    "(app/server/api.py::SPOT_SPATIAL_REFUSAL). Run the real solve; on this "
    "section it is the only solve there is."
)


class PreviewError(Exception):
    """A refusal, with the HTTP status the route should answer with."""

    def __init__(self, message: str, *, status: int = 422, **detail):
        super().__init__(message)
        self.message = message
        self.status = status
        self.detail = detail


# --------------------------------------------------------------------------- #
# The agreement, derived from the measured pairs
# --------------------------------------------------------------------------- #
def _pct(sorted_vals: Sequence[float], q: float) -> float:
    """Nearest-rank percentile.

    The SAME convention as ``experiments/surrogate_feasibility/derived.py``,
    deliberately: linear interpolation over 12 values puts the p90 of the
    delivery spread at 1.1757 % where the published record says 1.2001 %, and a
    feature whose error bar differs from the record it cites in the second
    decimal is a claim-surface drift waiting to be found.
    """
    if not sorted_vals:
        return float("nan")
    i = min(int(round(q * (len(sorted_vals) - 1))), len(sorted_vals) - 1)
    return sorted_vals[i]


def _f(v) -> Optional[float]:
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) else None


def _agreement_path() -> str:
    from .cache import repo_root

    return os.path.join(repo_root(), AGREEMENT_SOURCE)


def measured_agreement() -> Optional[Dict]:
    """The preview's error bar, recomputed from the measured pairs.

    ``None`` when the data file is absent, and the caller must then REFUSE the
    preview rather than fall back to a literal. That is the whole discipline of
    this function: the number and the statement of how far it may be from the
    real solve are one artefact, and a server that has lost the second does not
    get to serve the first.

    **Only pairs where the grid ACTUALLY DIFFERS are counted.** A pair whose
    two runs produced the same grid is a run compared with itself, its error is
    0.000 % by construction, and mixing those in is what made the study's first
    pass report a perfect result (sec. 9.1 of the record). The filter is on the
    measured grids, not on the platform name, so it cannot be fooled by a
    section whose platform string is wrong.
    """
    path = _agreement_path()
    try:
        with open(path) as fh:
            cb = json.load(fh)
    except (OSError, ValueError):
        return None

    pairs = cb.get("pairs") or []
    real, tautologies = [], 0
    for p in pairs:
        c, f = p.get("coarse") or {}, p.get("fine") or {}
        if c.get("status") != "ok" or f.get("status") != "ok":
            continue
        if c.get("grid") == f.get("grid"):
            tautologies += 1
            continue
        real.append(p)
    if not real:
        return None

    rel: List[float] = []
    abs_: List[float] = []
    for p in real:
        cv = _f((p.get("coarse") or {}).get(PREVIEWABLE_READOUT))
        fv = _f((p.get("fine") or {}).get(PREVIEWABLE_READOUT))
        if cv is None or fv is None or fv == 0:
            continue
        rel.append(abs(cv - fv) / abs(fv))
        abs_.append(abs(cv - fv))
    if not rel:
        return None
    rel.sort()
    abs_.sort()

    cw = sorted(_f(p["coarse"].get("wall_clock_s")) or 0.0 for p in real)
    fw = sorted(_f(p["fine"].get("wall_clock_s")) or 0.0 for p in real)
    loads = [x for x in
             ([_f(p["coarse"].get("loadavg_1m")) for p in real]
              + [_f(p["fine"].get("loadavg_1m")) for p in real])
             if x is not None]

    return {
        "readout": PREVIEWABLE_READOUT,
        "preview_dx_um": _f(cb.get("preview_dx_um")) or PREVIEW_DX_UM,
        "target_dx_um": _f(cb.get("shipped_dx_um")) or TARGET_DX_UM,
        "n_pairs": len(rel),
        "n_sections": len({p["section"] for p in real}),
        "n_compounds": len({p["drug"] for p in real}),
        "sections": sorted({p["section"] for p in real}),
        # A SET, alphabetical. Never an ordering, and no metric beside it.
        "compounds": sorted({p["drug"] for p in real}),
        "median_rel_err_pct": st.median(rel) * 100.0,
        "p90_rel_err_pct": _pct(rel, 0.90) * 100.0,
        "max_rel_err_pct": rel[-1] * 100.0,
        "median_abs_err": st.median(abs_),
        "max_abs_err": abs_[-1],
        # The two comparators that say whether the error bar is small. Both are
        # WITHIN-compound or across-DONOR figures; neither is a ranking.
        "across_donor_range": 0.5357238964376265,
        "across_donor_range_source":
            "results/surrogate_feasibility_2026-10-09.md sec. 5.2 - the "
            "delivery spread's full range over the 4 donors where per-cell "
            "spatial statistics are defensible (1.9936 / 2.3211 / 2.3787 / "
            "2.5293)",
        "nuisance_band_width": 0.202,
        "nuisance_band_source":
            "results/delivery_spread_noise.md sec. 0 - temozolomide's own span "
            "over the full nuisance grid (1.969 - 2.171). A WITHIN-compound "
            "figure; the compound separation that record also reports is not "
            "quoted here and no ordering between compounds is implied",
        # Cost. `results/performance_profile_2026-09-20.md` sec. 0.1-0.2 measures
        # FIVE different wall-clocks for one RunSpec on unchanged code and
        # concludes that what is stable is the SHAPE, not the absolutes -- so
        # the RATIO is the quotable part and the load it was taken under ships
        # with it.
        "coarse_median_s": st.median(cw),
        "target_median_s": st.median(fw),
        "cost_ratio": (st.median(fw) / st.median(cw)) if st.median(cw) else None,
        "loadavg_1m_min": min(loads) if loads else None,
        "loadavg_1m_max": max(loads) if loads else None,
        "n_cpu": sorted({p["coarse"].get("n_cpu") for p in real
                         if p["coarse"].get("n_cpu")}),
        # Reported, not hidden: how many stored pairs were discarded as
        # tautologies. A reader can check that the envelope is 12 and not 30.
        "pairs_discarded_same_grid": tautologies,
        "source": AGREEMENT_SOURCE,
        "record": "results/surrogate_feasibility_2026-10-09.md sec. 6",
    }


# --------------------------------------------------------------------------- #
# Refusals
# --------------------------------------------------------------------------- #
def _spot_sources() -> Tuple[str, ...]:
    """``_SPOT_SOURCES`` from ``app/server/api.py``, imported not copied.

    Deferred because ``api.create_app`` mounts this router, so a module-level
    import of ``api`` here is a cycle.
    """
    from .api import _SPOT_SOURCES

    return tuple(_SPOT_SOURCES)


def _section_row(section: str) -> Optional[Dict]:
    from .cache import list_sections

    return next((r for r in list_sections() if r.get("id") == section), None)


def _check_readout(requested) -> str:
    """Which readout the caller asked for, or a refusal naming the measurement.

    A caller who names nothing gets the previewable one. A caller who names the
    previewable one gets it. Everything else is refused with its own figure --
    including names this module has never heard of, because an unknown readout
    is one nobody measured a dx agreement for.
    """
    if requested is None:
        return PREVIEWABLE_READOUT
    name = str(requested).strip()
    if name == PREVIEWABLE_READOUT:
        return name
    if name in REFUSED_READOUTS:
        raise PreviewError(REFUSED_READOUTS[name], readout=name,
                           previewable=[PREVIEWABLE_READOUT])
    if "kill" in name.lower():
        raise PreviewError(
            "ABSOLUTE KILL FRACTION IS NEVER A DEFENSIBLE OUTPUT of this "
            "pipeline and is not previewed, approximated or estimated here at "
            "any resolution. This is not a statement about the coarse grid; it "
            "is the product's position on the quantity, and a faster route to "
            "it is not an improvement. WHAT TO DO: read the delivery spread "
            "instead - a within-map ratio, which is the comparison this "
            "product defends.",
            readout=name, previewable=[PREVIEWABLE_READOUT])
    raise PreviewError(
        f"NO MEASURED dx = {PREVIEW_DX_UM:g} / dx = {TARGET_DX_UM:g} "
        f"AGREEMENT EXISTS FOR {name!r}, so there is nothing to show beside a "
        "preview of it and this route will not invent one. The only readout a "
        f"coarse solve previews here is {PREVIEWABLE_READOUT!r} (the delivery "
        "spread), because it is the only one whose disagreement with the real "
        "run has been measured pair by pair. Measure the pairs and the "
        "envelope grows; until then, wait for the run.",
        readout=name, previewable=[PREVIEWABLE_READOUT])


def _check_section(section: str) -> Dict:
    if not section or not isinstance(section, str):
        raise PreviewError(
            "A preview needs a section. Name one from GET /api/sections.",
            field="section")
    row = _section_row(section)
    if row is None:
        raise PreviewError(
            f"Unknown section {section!r}. GET /api/sections lists what this "
            "server can reconstruct.", status=404, section=section)
    platform = str(row.get("platform") or row.get("source") or "").lower()
    if platform in _spot_sources():
        raise PreviewError(SPOT_PREVIEW_REFUSAL, status=422, section=section,
                           platform=platform,
                           clamped_dx_um=100.0, remedy="run the real solve")
    if section not in MEASURED_SECTIONS:
        raise PreviewError(
            f"NO PREVIEW ON {section!r}: its dx = {PREVIEW_DX_UM:g} / "
            f"dx = {TARGET_DX_UM:g} agreement has not been measured. The 12 "
            "pairs behind this feature cover "
            + ", ".join(MEASURED_SECTIONS)
            + " and nothing else, and this route will not stretch a figure "
            "taken on two sections over a third. The engine still runs here - "
            "POST /api/runs is unaffected. To extend the preview to this "
            "section, solve it at both resolutions for the compounds you care "
            "about and add the pairs to " + AGREEMENT_SOURCE + ".",
            section=section, measured_sections=list(MEASURED_SECTIONS))
    return row


def _check_available(row: Dict) -> None:
    """Can this server actually build the reconstruction?

    DELIBERATELY LAST, after every envelope check. The envelope checks are pure
    policy — they need no tissue and no solve — and they are the ones a caller
    most needs to hear, so a server with no ``data_cache/`` mounted still
    answers "this readout is not previewable" rather than masking it behind
    "this section is not here". It also keeps those refusals testable on a
    machine without the tissue data, which is where a refusal's wording is
    most likely to be edited.
    """
    if row.get("available", True):
        return
    raise PreviewError(
        f"Section {row.get('id')!r} is in the catalogue but is not "
        "reconstructable on this server: "
        + str(row.get("reason") or "no reason recorded")
        + ". A preview cannot be cheaper than a reconstruction that cannot be "
        "built. WHAT TO DO: populate the section's cache directory as "
        "`vivum sections` describes; POST /api/runs cannot run here either, so "
        "this is not a limitation of the preview.",
        status=409, section=row.get("id"))


def _check_compound(drug) -> str:
    if not drug or not isinstance(drug, str):
        raise PreviewError(
            "A preview needs a compound. Name one from GET /api/drugs.",
            field="drug")
    if drug not in MEASURED_COMPOUNDS:
        raise PreviewError(
            f"NO PREVIEW FOR {drug!r}: its dx = {PREVIEW_DX_UM:g} / "
            f"dx = {TARGET_DX_UM:g} agreement has not been measured. The "
            "measured set is "
            + ", ".join(MEASURED_COMPOUNDS)
            + " - listed alphabetically as a SET, with no ordering between "
            "them implied or available. A user-compiled compound "
            "(results/bring_your_own_compound.md) is a shipped surface and is "
            "deliberately NOT given a borrowed error bar here. Run it at both "
            "resolutions and add the pair to " + AGREEMENT_SOURCE + ".",
            drug=drug, measured_compounds=list(MEASURED_COMPOUNDS))
    return drug


def _check_target_dx(dx) -> float:
    target = _f(dx)
    if target is None or target <= 0:
        raise PreviewError(
            "A preview needs the dx_um the REAL run will use, so it can state "
            "which measured pair it is quoting.", field="dx_um")
    if target == PREVIEW_DX_UM:
        raise PreviewError(
            f"THIS RUN IS ALREADY AT dx = {PREVIEW_DX_UM:g} um, which is the "
            "preview's own resolution, so a preview would be the same solve at "
            "the same cost - a progress bar in front of nothing. Submit the "
            "run.",
            status=409, dx_um=target)
    if target != TARGET_DX_UM:
        raise PreviewError(
            f"NO MEASURED PAIR FOR dx = {PREVIEW_DX_UM:g} um against "
            f"dx = {target:g} um. The 12 pairs behind this feature compare "
            f"dx = {PREVIEW_DX_UM:g} with dx = {TARGET_DX_UM:g} and nothing "
            "else, so for this run there is no figure saying how far a preview "
            "could be from it - and a preview that cannot state that distance "
            "is just a wrong number shown quickly. This is not a tolerance "
            "that can be widened by judgement; it is a pair that has to be "
            "measured. WHAT TO DO: either set dx_um to "
            f"{TARGET_DX_UM:g} um - the resolution this product ships and the "
            "one the pairs cover - and preview that, or submit the run as it "
            "stands and wait for it; POST /api/runs is unaffected and will "
            "answer for this dx. Note that GET /api/defaults coarsens dx_um on "
            "a memory-constrained host (api.py::_default_dx_um returns 21 um at "
            "a 100,000-cell ceiling), and on such a host there is no preview "
            "because there is no measurement of that pair - which is a "
            "measurement to take, not a tolerance to widen.",
            dx_um=target, preview_dx_um=PREVIEW_DX_UM,
            measured_target_dx_um=TARGET_DX_UM)
    return target


# --------------------------------------------------------------------------- #
# The solve
# --------------------------------------------------------------------------- #
class _PreviewCtx:
    """The progress channel ``execute_run`` expects, with nowhere to store.

    A :class:`app.server.jobs._Ctx` writes stage, progress and a log into a
    ``Job``, which ``JobManager`` persists. This one writes into a list that is
    dropped when the request ends, because a preview must not become a row in
    the run store -- that is the route by which a fast number becomes a quoted
    one.
    """

    def __init__(self) -> None:
        self.stages: List[str] = []
        self.logs: List[str] = []

    def log(self, msg: str) -> None:
        self.logs.append(str(msg))

    def stage(self, name: str) -> None:
        self.stages.append(str(name))

    def check_cancelled(self) -> None:
        return None


def _coarse_spec(spec: Dict) -> Dict:
    """The caller's RunSpec with ONE field changed, and the extras switched off.

    The preview differs from the real run in ``dx_um`` ALONE among the fields
    that reach the physics. The two fields that are forced are both pure
    read-outs with no effect on the solve: ``timelapse_frames = 0`` (frames the
    preview will not show) and ``attribute_mechanism = False`` (a per-cell
    attribution the preview does not render). ``jobs.py`` records that a run is
    bit-identical with mechanism attribution on or off.

    ``resistance_off`` is NOT forced: it is a counterfactual ARM of a report,
    and silently flipping an arm would make the preview a preview of a
    different question.
    """
    out = dict(spec)
    out["dx_um"] = PREVIEW_DX_UM
    out["timelapse_frames"] = 0
    out["attribute_mechanism"] = False
    out.pop("kind", None)  # a TPL protocol program is not a RunSpec
    return out


#: Keys that must never appear anywhere in a preview payload, at any depth.
#: Checked on the FINISHED body rather than trusted to the input filter,
#: because the filter is what a future edit relaxes. The same discipline
#: ``app/server/cross_section.py`` applies to its own payload.
_FORBIDDEN_SUBSTRINGS: Tuple[str, ...] = (
    "kill", "penetration_length", "penetration_index", "auc", "coverage",
    "E_p90_uM", "E_p10_uM", "E_p50_uM",
)


def _audit_payload(body: Dict) -> Dict:
    """Refuse to return a body that has lost its label or its error bar.

    Three things make this payload honest and all three are checked here, after
    the fact, rather than assumed from the code above:

    1. the number is labelled a preview and marked not quotable;
    2. the measured agreement travels WITH it;
    3. no refused readout and no absolute level rode along.

    A failure is a 500 on purpose. A preview that has quietly shed its caveat
    is worse than no preview, so this route would rather break loudly than
    serve one.
    """
    if body.get("is_preview") is not True or body.get("quotable") is not False:
        raise PreviewError(
            "preview payload audit FAILED: a preview body must carry "
            "is_preview=true and quotable=false.", status=500)
    if not body.get("preview_label"):
        raise PreviewError(
            "preview payload audit FAILED: the preview label is missing, and "
            "an unlabelled fast number is the one thing this route exists to "
            "prevent.", status=500)
    ag = body.get("agreement") or {}
    for k in ("max_rel_err_pct", "max_abs_err", "n_pairs", "median_rel_err_pct"):
        if _f(ag.get(k)) is None:
            raise PreviewError(
                "preview payload audit FAILED: the measured agreement is "
                f"missing {k!r}. The number and the statement of how far it "
                "may be from the real solve ship together or not at all.",
                status=500)
    if body.get("value") is not None and not ag:
        raise PreviewError(
            "preview payload audit FAILED: a preview value with no agreement "
            "block.", status=500)

    def walk(node, path="") -> None:
        if isinstance(node, dict):
            for k, v in node.items():
                low = str(k).lower()
                for bad in _FORBIDDEN_SUBSTRINGS:
                    if bad.lower() in low:
                        raise PreviewError(
                            "preview payload audit FAILED: key "
                            f"{path + '.' + str(k)!r} names {bad!r}. A preview "
                            "carries ONE readout, the delivery spread, and no "
                            "absolute level and no refused readout may ride "
                            "along with it.", status=500)
                walk(v, path + "." + str(k))
        elif isinstance(node, list):
            for i, v in enumerate(node):
                walk(v, f"{path}[{i}]")

    # The refusal PROSE legitimately names the readouts it refuses, so the
    # key-name walk runs over the payload with those blocks removed.
    audited = {k: v for k, v in body.items()
               if k not in ("refused", "preview_label", "agreement")}
    walk(audited)
    return body


def preview_payload(jobs, spec: Dict) -> Dict:
    """One coarse solve, with its measured agreement, labelled as a preview."""
    from .jobs import execute_run, jsonable

    if not isinstance(spec, dict):
        raise PreviewError("the body must be a RunSpec object.", status=400)

    agreement = measured_agreement()
    if agreement is None:
        raise PreviewError(
            "THE PREVIEW IS UNAVAILABLE ON THIS SERVER because the measured "
            "agreement between the coarse and the shipped grid is not present: "
            + AGREEMENT_SOURCE + " is missing, unreadable, or holds no pair "
            "whose two runs used different grids. This route deliberately "
            "carries no second copy of those figures, so when the measurement "
            "is gone the feature is gone rather than the error bar being gone. "
            "POST /api/runs is unaffected - run the real solve.",
            status=503, source=AGREEMENT_SOURCE)

    readout = _check_readout(spec.get("readout"))
    section = _check_section(str(spec.get("section") or ""))
    drug = _check_compound(spec.get("drug"))
    target_dx = _check_target_dx(spec.get("dx_um"))
    # Last, and see its docstring: a server with no tissue mounted should still
    # be able to tell a caller that a readout is not previewable.
    _check_available(section)

    coarse = _coarse_spec(spec)
    t0 = time.time()
    ctx = _PreviewCtx()
    try:
        result = execute_run(coarse, ctx, cache=_recon_cache(jobs))
    except Exception as exc:  # noqa: BLE001
        from .cache import redact_exception

        # Hoisted into its own statement rather than concatenated inline. The
        # text is identical either way; what changes is that the concatenation
        # below no longer MENTIONS `exc`, and
        # tests/test_response_redaction.py's sweep flags any `+` expression
        # that names a caught exception -- correctly, since it cannot tell
        # from the shape whether the redactor was applied or bypassed. The
        # same hoist was applied to cache.py::state_dir_problem on
        # 2026-09-20 for the same reason, so this is the house idiom and not
        # a one-off dodge of a guard.
        detail = redact_exception(exc, "preview coarse solve", brief=True)

        raise PreviewError(
            "The coarse solve did not finish, so there is no preview: "
            + detail
            + ". This says nothing about whether the real run will succeed - "
            "POST /api/runs is a different solve on a different grid.",
            status=500, stage=(ctx.stages[-1] if ctx.stages else None))
    wall = time.time() - t0

    # NOTHING IS PERSISTED. `_payloads` holds the geometry and field blobs a
    # real run writes to disk; dropping it here is what keeps a preview out of
    # the run store, out of /api/runs, and therefore out of every surface that
    # could quote it.
    result.pop("_payloads", None)

    xm = result.get("extra_metrics") or {}
    rep = result.get("report") or {}
    p90, p10 = _f(xm.get("E_p90_uM")), _f(xm.get("E_p10_uM"))
    value = (p90 / p10) if (p90 is not None and p10 is not None and p10 > 0) else None
    null_reason = None
    if value is None:
        null_reason = (
            "The coarse solve produced no delivery spread on this section: one "
            "of the two exposure percentiles is absent or the lower one is "
            "zero. The engine declining to produce a number is carried through "
            "rather than filled in - which is the point of previewing with the "
            "engine instead of with a model, since a model would have returned "
            "something here."
        )

    body: Dict = {
        # --- what this is, stated before the number ----------------------- #
        "is_preview": True,
        "quotable": False,
        "preview_label": PREVIEW_LABEL,
        "readout": readout,
        "readout_name": "delivery spread  P90 / P10",
        "value": value,
        "unit": "x",
        "null_reason": null_reason,
        # The denominator travels with the fold-change (WORKSPACE_CONTRACT
        # sec. 6: every number carries its unit and its denominator).
        "n_tumor": rep.get("n_tumor"),
        "derivation":
            "extra_metrics.E_p90_uM / extra_metrics.E_p10_uM from a dx = "
            f"{PREVIEW_DX_UM:g} um solve, over the tumour cell set. A ratio "
            "between two points of ONE map: within-map heterogeneity only, "
            "never a permeability result, never an absolute level, and never "
            "comparable across compounds.",
        # --- how far it may be from the real run -------------------------- #
        "agreement": agreement,
        # --- what it is a preview OF -------------------------------------- #
        "preview_dx_um": PREVIEW_DX_UM,
        "target_dx_um": target_dx,
        "preview_grid": rep.get("grid"),
        "preview_wall_clock_s": wall,
        "section": section.get("id"),
        "drug": drug,
        # --- what it refused to preview, so the absence is visible -------- #
        "refused": [{"readout": k, "why": v}
                    for k, v in sorted(REFUSED_READOUTS.items())],
        # The engine's own caveats, verbatim, from the same place a real run's
        # provenance block holds them. A preview that dropped them would be
        # claiming MORE than the run it previews -- including the RESOLUTION
        # caveat, which on this body reads `dx = 30 um` and so states the
        # preview's own coarseness in the engine's words rather than this
        # module's. ``claim_set`` rides along for the reason jobs.py puts it
        # beside them: whatever lifts the caveats out of a record carries the
        # one field that says how old they are.
        "caveats": list((result.get("provenance") or {}).get("caveats") or []),
        "claim_set": (result.get("provenance") or {}).get("claim_set"),
        "superseded_by":
            "the real run. Submit POST /api/runs with your own spec; this "
            "preview has no run_id, is stored nowhere, and will be replaced "
            "rather than compared against.",
    }
    return _audit_payload(jsonable(body))


def _recon_cache(jobs):
    """The reconstruction cache the real queue uses, when it has one.

    A preview at dx = 30 um builds a DIFFERENT cache entry from a run at
    dx = 12 um (``make_key`` includes dx), so sharing the cache cannot make the
    two solves contaminate each other -- it only means the preview's
    reconstruction is reused by the next preview.
    """
    from .cache import RECON_CACHE

    return getattr(jobs, "cache", None) or RECON_CACHE


# --------------------------------------------------------------------------- #
# Router
# --------------------------------------------------------------------------- #
def build_preview_router(*, jobs) -> APIRouter:
    """Mount ``POST /api/preview`` and ``GET /api/preview/agreement``.

    Gated by the same middleware every other ``/api/`` route is: ``create_app``
    tests ``get_route_path(request.scope).startswith("/api/")``, so a router
    mounted there is authenticated by construction.
    """
    router = APIRouter()

    @router.get("/api/preview/agreement")
    def get_preview_agreement() -> Dict:
        """The error bar, and the refusals, WITHOUT running anything.

        Served separately so a surface can show what a preview would and would
        not tell you before spending the 8.7 s -- and so the absence of a
        measurement is readable rather than only discoverable by being refused.
        """
        from .jobs import jsonable

        ag = measured_agreement()
        return jsonable({
            "available": ag is not None,
            "agreement": ag,
            "unavailable_reason": None if ag is not None else (
                "the measured agreement is not present on this server: "
                + AGREEMENT_SOURCE + " is missing, unreadable, or holds no "
                "pair whose two runs used different grids"),
            "preview_dx_um": PREVIEW_DX_UM,
            "target_dx_um": TARGET_DX_UM,
            "previewable": [PREVIEWABLE_READOUT],
            "preview_label": PREVIEW_LABEL,
            "measured_sections": list(MEASURED_SECTIONS),
            "measured_compounds": list(MEASURED_COMPOUNDS),
            "spot_refusal": SPOT_PREVIEW_REFUSAL,
            "refused": [{"readout": k, "why": v}
                        for k, v in sorted(REFUSED_READOUTS.items())],
        })

    @router.post("/api/preview")
    async def post_preview(request: Request) -> Dict:
        """One coarse solve of the delivery spread, with its measured error bar.

        Synchronous on purpose. The median coarse solve on the measured
        sections is 8.6678 s against 58.7711 s for the real run, and a preview
        that needed its own queue, its own poll and its own cancel would be a
        second run-shaped object a user could mistake for a run.
        """
        from .api import _json_body

        try:
            spec = await _json_body(request)
        except HTTPException:
            raise
        except Exception:  # noqa: BLE001
            raise HTTPException(status_code=400,
                                detail={"error": "the body is not valid JSON."})
        try:
            return preview_payload(jobs, spec)
        except PreviewError as exc:
            raise HTTPException(status_code=exc.status,
                                detail={"error": exc.message, **exc.detail})

    return router
