"""Cell states from MEASURED markers — and an honest refusal where they cannot be.

WHY THIS MODULE EXISTS
----------------------
The shipping cell-state labelling is a banding of the SIMULATED oxygen field
(``cells.agents.update_cell_states``: ``o < t`` HYPOXIC, ``t <= o < 2t``
QUIESCENT, else VIABLE), and those labels gate drug kill through
``PDConfig.hypoxia_resistance`` (0.35) and ``PDConfig.quiescent_resistance``
(0.50). Three independent analyses in ``results/`` have now shown that chain
cannot carry an absolute number:

* ``results/dose_schedule.md`` §1.3 — for temozolomide, whose per-cell
  sensitivity is uniform 1.0 on the 319-gene panel, **~89% of the kill ceiling
  comes from the QUIESCENT band alone**. The ceiling is `sensitivity x
  resistance`, so it IS the state labelling.
* ``results/state_validation.md`` §2.2 — the model's tumour VIABLE fraction runs
  **0.60% to 98.19% (165x)** across two tissue-mask conventions and three grid
  spacings. Neither the mask nor the mesh is a measurement. Under the shipping
  ``geometric`` mask the HYPOXIC band is EMPTY at every dx, so
  ``hypoxia_resistance`` acts on no cells at all. (That table is KDE-era;
  re-measured on the same six-configuration grid under the current DISCRETE
  vessel source, ``results/cell_state_impact.md`` §2.2 gets 0.00%-10.04% — the
  VIABLE band itself goes empty, so the swing is unbounded. The conclusion does
  not depend on which source is in force.)
* ``results/oxygen_demotion.md`` §1 — the oxygen PDE's own effective diffusion
  length (310-801 um) is 2-5x the measured oxygen diffusion distance (~150 um)
  because ``OxygenConfig.consumption`` has never been fitted to a measured oxygen
  tension. (This bullet also used to say the PDE was out-predicted as a hypoxia
  predictor by a one-line distance-to-nearest-vessel transform. **That was an
  inverted inequality and is withdrawn** — the solve scores rho -0.1523 against
  the baseline's -0.1204, and wins. The length-scale criticism above is
  INDEPENDENT of that and survives it: a Spearman rho is invariant to monotone
  rescaling, so the solve can win the correlation and still have the wrong length
  scale. The banding this module does depends on the length scale, not on the
  rank correlation, so nothing here changes.)

So the state label — the single largest term in every absolute kill number this
project reports — is currently a bookkeeping choice. This module supplies the
alternative that is a MEASUREMENT: proliferative / quiescent from measured
proliferation transcripts, hypoxic from measured hypoxia-response transcripts,
with an explicit **UNMEASURED** verdict wherever the panel cannot support the
call.

WHAT IT WILL AND WILL NOT DO
----------------------------
1. **It will refuse.** Where a panel carries no gene of a marker set, that state
   is declared UNMEASURED and **no cell is assigned it**. The 319-gene Xenium
   Human Brain panel carries no VEGFA / CA9 / ADM / NDRG1 / SLC2A1, so on that
   section HYPOXIC is unmeasured — it has not been simulated away and must not be
   reported as if it had.
2. **It will not fabricate a per-cell call out of dropout.** On the 319-gene
   panel MKI67 fires in 0.35% of all cells (0.55% of tumour cells) and TOP2A in
   1.08% (0.99% of tumour cells), for a union of 1.42% of tumour cells, against
   published GBM Ki-67 labelling indices of 19-40% (:data:`LITERATURE_KI67`) —
   13x below the lower bound on the union, and up to 114x below on MKI67 alone. A per-cell "this cell is cycling" call at that detection rate is a
   probe-sensitivity readout, not biology. When the union detection of the
   specific markers falls below :data:`PER_CELL_DETECTION_FLOOR` the module
   switches to a **NEIGHBOURHOOD** call (see below) and says so in the coverage
   report; it never silently produces per-cell labels from a 0.35% detection
   rate.
3. **It will report what it used.** Every return carries a coverage block naming
   each requested marker, whether the panel measures it, its measured detection
   rate, the CALL LEVEL that was therefore chosen, and the list of states that
   are UNMEASURED on that panel.

THE TWO CALL LEVELS
-------------------
``'per_cell'``
    A cell is positive iff at least one of the present markers has >= 1 count in
    that cell. This is the only convention-free binarisation available: the
    threshold is "a probe fired", not a number someone picked. Admitted only when
    the union detection rate reaches :data:`PER_CELL_DETECTION_FLOOR`.

``'neighbourhood'``
    The per-cell call is dropout-limited, so positivity is aggregated over a
    fixed PHYSICAL window: for each cell, the fraction of cells within
    :data:`NEIGHBOURHOOD_RADIUS_UM` that are marker-positive. A cell is called
    proliferative iff its local positive fraction strictly EXCEEDS the
    section-wide positive rate, i.e. iff it sits in a neighbourhood ENRICHED for
    measured proliferation. What this label means, exactly: *the tissue around
    this cell carries more measured proliferation signal than this section's
    average*. It is **not** a per-cell Ki-67 call and the resulting fraction is
    **not** a proliferative index — dropout of this severity destroys the level,
    and only the spatial ordering survives. The window is in micrometres and is
    computed from the measured cell coordinates, so it is **completely
    independent of the simulation mesh** — which is the property the oxygen route
    does not have. It is NOT independent of the window radius, though: the
    binarisation is flagged ``'CONVENTION'`` in the coverage block and the radius
    sweep is reported (see :data:`NEIGHBOURHOOD_RADIUS_UM`). A reader who needs a
    level rather than an ordering should treat this axis as UNMEASURED and use
    the continuous ``per_cell['proliferation_neighbourhood']['fraction']``
    covariate instead.

MARKER CHOICES, AND THEIR SPECIFICITY
-------------------------------------
Proliferation (:data:`PROLIFERATION_MARKERS`): ``MKI67``, ``TOP2A``, ``PCNA``.
Only ``MKI67`` and ``TOP2A`` (:data:`CELL_CYCLE_SPECIFIC_MARKERS`) are used for
the positive call by default. ``PCNA`` is measured, reported, and EXCLUDED from
the default rule on purpose: PCNA protein has a long half-life and the transcript
is expressed in non-cycling and DNA-repairing cells, so a PCNA+ call over-counts
proliferation. It fires in 49.5% of tumour cells on the 319-gene panel against
MKI67's 0.55% — a 90x gap, on the same cells, that is specificity and not
biology. Pass
``include_low_specificity=True`` to fold it in and see what it does.

Hypoxia (:data:`HYPOXIA_MARKERS`): ``VEGFA``, ``CA9``, ``ADM``, ``NDRG1``,
``SLC2A1``. These are hypoxia-**RESPONSE** transcripts — the output of a
transcriptional programme with its own kinetics and its own non-hypoxic inducers
— not an oxygen tension. ``HIF1A`` is deliberately absent: oxygen sensing acts on
the HIF-1alpha PROTEIN via O2-dependent prolyl hydroxylation and VHL-mediated
proteolysis (Ivan et al., Science 2001;292(5516):464-468, PMID 11292862,
doi:10.1126/science.1059817; Jaakkola et al., Science 2001;292(5516):468-472,
PMID 11292861, doi:10.1126/science.1059796), so its mRNA has no obligation to
track oxygen, and ``results/state_validation.md`` §1.1 measured it at rank
249/480 (rho = +0.0085) against the simulated field.

**THE HONEST LIMIT ON THE HYPOXIC BAND.** A detection-based hypoxia call is
traceable to a measurement but is NOT calibrated to a hypoxic fraction. VEGFA
fires in 63.4% of IO-panel cells and ADM/SLC2A1 in 12-15% of GSE324268 cells;
these genes have real non-hypoxic expression, so "positive for a hypoxia-response
transcript" over-calls hypoxia relative to any pO2-based definition. Of the set
only ``CA9`` (:data:`HYPOXIA_SPECIFIC_MARKERS`) is near-restricted to severe
hypoxia. The module therefore reports the HYPOXIC fraction as an ENRICHMENT
label with an UNCALIBRATED SIZE, and
``hypoxia_markers=HYPOXIA_SPECIFIC_MARKERS`` is available as the conservative
variant. Nothing here calibrates an absolute
hypoxic fraction, and nothing here should be read as doing so.

THE ONE-LABEL COLLAPSE
----------------------
Proliferation and hypoxia are two axes; ``core.CellState`` has one label. The
collapse used here is **HYPOXIC takes precedence**: a hypoxia-marker-positive
cell is labelled HYPOXIC whatever its proliferation call, matching the oxygen
route's band ordering (the lowest oxygen band wins). That precedence is a
convention forced by the schema, it is stated in the returned provenance, and it
is the only place in this module where a choice is made that is not a marker
threshold.

WHAT IS NOT TOUCHED
-------------------
``core``, ``sim.engine``, ``physics`` and ``cells.pkpd`` are unmodified. The
shipping state route is still the oxygen banding
(``reconstruct.real_pipeline.reconstruct_real(state_mode='oxygen')``), and
``tests/test_cell_state.py`` asserts the default is bit-identical.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy import sparse as _sp
from scipy.spatial import cKDTree

from ..core.microenvironment import CellState, CellType

__all__ = [
    "CELL_CYCLE_SPECIFIC_MARKERS",
    "HYPOXIA_MARKERS",
    "HYPOXIA_SPECIFIC_MARKERS",
    "LITERATURE_KI67",
    "LOW_SPECIFICITY_MARKERS",
    "NEIGHBOURHOOD_MIN_CELLS",
    "NEIGHBOURHOOD_RADIUS_UM",
    "PER_CELL_DETECTION_FLOOR",
    "PROLIFERATION_MARKERS",
    "STATE_MODES",
    "marker_detection",
    "neighbourhood_positive_fraction",
    "panel_coverage",
    "state_fractions_from_states",
    "state_from_markers",
    "state_from_markers_micro",
]

# ---------------------------------------------------------------------------
# Marker sets
# ---------------------------------------------------------------------------

#: Measured proliferation transcripts looked for, in preference order.
PROLIFERATION_MARKERS: Tuple[str, ...] = ("MKI67", "TOP2A", "PCNA")

#: The subset used for the POSITIVE CALL by default. Ki-67 (MKI67) is expressed
#: in G1/S/G2/M and absent in G0 — the definition of the labelling index — and
#: TOP2A peaks in G2/M. Both are cell-cycle restricted.
CELL_CYCLE_SPECIFIC_MARKERS: Tuple[str, ...] = ("MKI67", "TOP2A")

#: Measured, reported, and excluded from the default positive call. See the
#: module docstring: PCNA is expressed in non-cycling and DNA-repairing cells and
#: its detection rate (45.6% on the 319-gene panel) is 130x MKI67's on the same
#: cells, which is a specificity gap rather than a biological one.
LOW_SPECIFICITY_MARKERS: Tuple[str, ...] = ("PCNA",)

#: Measured hypoxia-RESPONSE transcripts looked for, in preference order. NOT an
#: oxygen tension; HIF1A deliberately excluded (see the module docstring).
HYPOXIA_MARKERS: Tuple[str, ...] = ("VEGFA", "CA9", "ADM", "NDRG1", "SLC2A1")

#: The near-hypoxia-restricted subset, offered as the conservative variant. CA9
#: is a direct HIF-1 target whose expression in tumours is largely confined to
#: severely hypoxic / perinecrotic regions; the rest of
#: :data:`HYPOXIA_MARKERS` have substantial normoxic expression.
HYPOXIA_SPECIFIC_MARKERS: Tuple[str, ...] = ("CA9",)

#: Published GBM Ki-67 labelling indices, as the external scale against which a
#: measured per-cell detection rate is judged admissible. All PubMed-verified in
#: ``results/state_validation.md`` §2.1:
#:   19% — Fudaba et al., World Neurosurg 2021;149:e469-e480 (PMID 33567368)
#:   27% — Alkhaibary et al., Hematol Oncol Stem Cell Ther 2018;12(2):82-88 (PMID 30552865)
#:   30% — Zheng et al., Zhonghua Bing Li Xue Za Zhi 2021;50(3):229-235 (PMID 33677887)
#:   40% — Caramanti et al., Cureus 2022;14(4):e24246 (PMID 35602791)
#: The reference standard is itself imprecise (Nielsen et al., Diagn Pathol
#: 2018;13:38, PMID 29885671: inter-observer kappa 0.32 whole-slide, 0.04
#: hot-spot, n=235 gliomas) — but not by two orders of magnitude.
LITERATURE_KI67: Tuple[float, float] = (0.19, 0.40)

# ---------------------------------------------------------------------------
# Decision constants
# ---------------------------------------------------------------------------

#: Minimum union detection rate for a PER-CELL call to be admitted.
#:
#: Anchored, not picked: the published GBM Ki-67 range starts at 19%
#: (:data:`LITERATURE_KI67`), and 0.05 is that lower bound divided by ~4. A panel
#: detecting the marker in >= 5% of cells is within a factor of ~4-8 of the
#: published biology and a per-cell call is a defensible (if conservative)
#: readout of it; a panel detecting it in 0.35% is 54x below the lower bound and
#: a per-cell call there is reporting probe sensitivity. The factor of 4 is a
#: declared tolerance — it is the one number in this module that is a judgement —
#: and both the rate and the resulting call level are always returned so a reader
#: can apply a different one.
PER_CELL_DETECTION_FLOOR = 0.05

#: Radius (um) of the neighbourhood window used when the per-cell call is
#: inadmissible.
#:
#: **Chosen by a stated rule, and its sensitivity is reported rather than hidden.**
#: The rule: the smallest swept radius at which >95% of windows meet the
#: :data:`NEIGHBOURHOOD_MIN_CELLS` floor on the SPARSEST of the three cached
#: sections (the 319-gene brain section, 23,812 tumour cells over ~41 mm^2).
#: Measured window occupancy there, tumour population:
#:
#:     R (um)    median cells/window    windows below the floor
#:       50              8                       67.8%
#:       75             16                       10.8%
#:      100             27                        0.9%
#:      150             58                        0.2%
#:
#: 50 um — the radius ``validation.state_validation.DENSITY_RADIUS_UM`` uses for
#: its local-crowding baseline — is therefore DEGENERATE on that section and is
#: not the default. 100 um is also the right order for a "local niche": it is
#: within a factor of 1.5 of the Thomlinson & Gray oxygen-diffusion distance
#: (150 um, PMID 13304213), the length scale over which tumour tissue is
#: physiologically homogeneous.
#:
#: **The radius is a convention and it moves the answer.** On the brain section
#: the resulting enriched (proliferative) tumour fraction runs 4.83% -> 42.50%
#: over R = 50 -> 150 um (8.8x), because as R -> 0 only the positive cells themselves are
#: enriched (so the fraction tends to the detection rate) and as R -> infinity
#: every window tends to the section rate (so the fraction tends to 0). That
#: sweep is reported in ``results/cell_state_impact.md`` and is the reason the
#: neighbourhood LABEL is flagged ``binarisation='CONVENTION'``: its ORDERING is
#: measured, its FRACTION is not.
#:
#: The window is in MICROMETRES on the measured coordinates, so — unlike the
#: oxygen route — it does not depend on the simulation mesh at all.
NEIGHBOURHOOD_RADIUS_UM = 100.0

#: Minimum cells in a window for its positive fraction to be used. Below this the
#: window's fraction is replaced by the section-wide rate — the least informative
#: value available, and therefore the honest filler (the same convention
#: ``physics.oxygen_proxy._nadaraya_watson`` uses for unsupported grid cells).
NEIGHBOURHOOD_MIN_CELLS = 10

#: The state routes this project can take. ``'oxygen'`` is and remains the
#: shipping default (see ``reconstruct.real_pipeline.reconstruct_real``).
#:
#: ``'vessel_distance'`` is implemented in
#: :func:`reconstruct.real_pipeline._install_states` rather than here: it bands
#: ``physics.oxygen_proxy.oxygen_from_vessel_distance`` through the UNMODIFIED
#: ``cells.agents.update_cell_states``, so it reads no transcript and belongs to
#: neither this module's marker machinery nor the default's PDE. It is listed here
#: because this tuple is the one enumeration of the state routes, and
#: ``results/cell_state_impact.md`` scores all three against the same bar.
#:
#: A FOURTH route, ML-imputed hypoxia, exists but is DELIBERATELY NOT a member of
#: this tuple. ``STATE_MODES`` is mirrored, and asserted equal, in
#: ``app/server/cache.py`` and ``tmesim/protocol/ast.py`` (three independent
#: copies kept in lockstep by ``tests/test_app_server.py`` and
#: ``tests/test_protocol.py``) — extending the DSL/server surface for a route
#: this session validated on THREE real Visium sections was judged out of
#: proportion to that evidence base. ``tmesim.ml.panel_imputation`` and
#: ``reconstruct.real_pipeline.install_ml_imputed_hypoxia`` implement the same
#: opt-in / refuse-unless-requested / distinct-provenance-tag discipline this
#: tuple's routes use, as an explicit function a caller applies to an ALREADY
#: reconstructed ``Microenvironment`` rather than as a ``state_mode=`` string.
#: See ``results/panel_gap_imputation_sweep.md``.
STATE_MODES: Tuple[str, ...] = ("oxygen", "vessel_distance", "markers")


# ---------------------------------------------------------------------------
# Panel introspection
# ---------------------------------------------------------------------------
def _counts_column(expression, genes: np.ndarray, gene: str) -> Optional[np.ndarray]:
    """Raw per-cell counts for one panel gene, or ``None`` if the panel lacks it."""
    hit = np.flatnonzero(genes == gene)
    if hit.size == 0:
        return None
    j = int(hit[0])
    if _sp.issparse(expression):
        return np.asarray(expression[:, [j]].todense()).ravel().astype(np.float64)
    return np.asarray(expression)[:, j].astype(np.float64)


def marker_detection(
    expression,
    genes: Sequence[str],
    marker_set: Sequence[str],
    subset: Optional[np.ndarray] = None,
) -> Dict[str, object]:
    """Which of ``marker_set`` the panel measures, and how often each fires.

    Parameters
    ----------
    expression : sparse matrix or ndarray, (n_cells, n_genes)
        Raw counts.
    genes : sequence of str
        Panel gene symbols aligned to ``expression``'s columns.
    marker_set : sequence of str
        Markers to look for.
    subset : boolean ndarray (n_cells,), optional
        Restrict the detection rates to this population (e.g. tumour cells).

    Returns
    -------
    dict
        ``present`` / ``missing`` (lists), ``detection`` (gene -> nonzero
        fraction), ``union_detection`` (fraction positive for >= 1 present
        marker), ``positive`` (boolean array over the FULL population, ``None``
        when nothing is present) and ``n_scored``.
    """
    genes = np.asarray(genes, dtype=str)
    n = int(expression.shape[0])
    sel = np.ones(n, dtype=bool) if subset is None else np.asarray(subset, dtype=bool)
    if sel.shape != (n,):
        raise ValueError(f"subset shape {sel.shape} != ({n},)")

    present: List[str] = []
    missing: List[str] = []
    detection: Dict[str, float] = {}
    positive = np.zeros(n, dtype=bool)
    for g in marker_set:
        raw = _counts_column(expression, genes, str(g))
        if raw is None:
            missing.append(str(g))
            continue
        present.append(str(g))
        hit = raw > 0
        detection[str(g)] = float(hit[sel].mean()) if sel.any() else 0.0
        positive |= hit

    return {
        "requested": [str(g) for g in marker_set],
        "present": present,
        "missing": missing,
        "detection": detection,
        "union_detection": (float(positive[sel].mean()) if (present and sel.any())
                            else 0.0),
        "positive": (positive if present else None),
        "n_scored": int(sel.sum()),
    }


def panel_coverage(
    expression,
    genes: Sequence[str],
    subset: Optional[np.ndarray] = None,
    proliferation_markers: Sequence[str] = CELL_CYCLE_SPECIFIC_MARKERS,
    hypoxia_markers: Sequence[str] = HYPOXIA_MARKERS,
    per_cell_detection_floor: float = PER_CELL_DETECTION_FLOOR,
) -> Dict[str, object]:
    """Marker coverage + the CALL LEVEL each state axis therefore gets.

    Pure introspection: no states are produced. Use it to ask "can this panel
    support a marker-derived state at all?" before paying for a reconstruction.
    """
    prol = marker_detection(expression, genes, proliferation_markers, subset)
    hyp = marker_detection(expression, genes, hypoxia_markers, subset)
    prol_level, prol_reason = _call_level(
        prol, "proliferation", float(per_cell_detection_floor)
    )
    hyp_level, hyp_reason = _call_level(
        hyp, "hypoxia", float(per_cell_detection_floor)
    )
    unmeasured = []
    if prol_level == "UNMEASURED":
        unmeasured += ["VIABLE(proliferative)", "QUIESCENT"]
    if hyp_level == "UNMEASURED":
        unmeasured.append("HYPOXIC")
    return {
        "n_panel_genes": int(np.asarray(genes, dtype=str).size),
        "n_cells": int(expression.shape[0]),
        "n_scored": prol["n_scored"],
        "per_cell_detection_floor": float(per_cell_detection_floor),
        "proliferation": {**{k: v for k, v in prol.items() if k != "positive"},
                          "call_level": prol_level, "reason": prol_reason,
                          "binarisation": _BINARISATION[prol_level]},
        "hypoxia": {**{k: v for k, v in hyp.items() if k != "positive"},
                    "call_level": hyp_level, "reason": hyp_reason,
                    "binarisation": _BINARISATION[hyp_level]},
        "unmeasured_states": unmeasured,
    }


#: What the binary label costs at each call level. ``'per_cell'`` costs nothing —
#: the threshold is "a probe fired". ``'neighbourhood'`` costs a window radius,
#: and the resulting FRACTION moves with it (see :data:`NEIGHBOURHOOD_RADIUS_UM`),
#: so only the ordering is measured.
_BINARISATION: Dict[str, Optional[str]] = {
    "per_cell": (
        "MEASURED — a cell is positive iff a probe fired in it. No expression "
        "threshold, no quantile, no fitted constant."
    ),
    "neighbourhood": (
        "CONVENTION — enrichment above the section-wide positive rate inside a "
        "fixed-radius window. The ORDERING is measured; the resulting FRACTION "
        "depends on the window radius (4.83% -> 42.50%, 8.8x, over R = 50 -> 150 "
        "um on the 319-gene section) and is NOT a measured index. Treat the level as "
        "unmeasured and use the continuous neighbourhood fraction instead."
    ),
    "UNMEASURED": None,
}


def _window_summary(block: Optional[Dict[str, object]]) -> Optional[Dict[str, object]]:
    """Window diagnostics for a neighbourhood axis, minus the per-cell array."""
    if block is None:
        return None
    return {k: v for k, v in block.items() if k != "fraction"}


def _call_level(block: Dict[str, object], axis: str, floor: float) -> Tuple[str, str]:
    """Decide ``'per_cell'`` / ``'neighbourhood'`` / ``'UNMEASURED'`` and say why."""
    present = list(block["present"])  # type: ignore[arg-type]
    missing = list(block["missing"])  # type: ignore[arg-type]
    u = float(block["union_detection"])  # type: ignore[arg-type]
    if not present:
        return "UNMEASURED", (
            f"no {axis} marker of {block['requested']} is on this panel (all "
            f"absent: {missing}). {axis.capitalize()} is UNMEASURED on this "
            "section — it has not been simulated away and must not be reported as "
            "if it had."
        )
    if u >= floor:
        return "per_cell", (
            f"union detection {u * 100:.2f}% over {present} is at or above the "
            f"{floor * 100:.0f}% floor, so a per-cell 'a probe fired' call is "
            "admitted. No expression threshold is chosen anywhere."
        )
    return "neighbourhood", (
        f"union detection {u * 100:.3f}% over {present} is BELOW the "
        f"{floor * 100:.0f}% floor "
        + (f"(published GBM Ki-67 index {LITERATURE_KI67[0] * 100:.0f}-"
           f"{LITERATURE_KI67[1] * 100:.0f}%, i.e. this panel is "
           f"{LITERATURE_KI67[0] / max(u, 1e-12):.0f}x below the lower bound), "
           if axis == "proliferation" else "")
        + "so a PER-CELL call would be reporting probe sensitivity rather than "
          "biology. Aggregated to a neighbourhood positive fraction instead; the "
          "resulting fraction is an ENRICHMENT ordering, NOT an index."
    )


# ---------------------------------------------------------------------------
# The neighbourhood aggregation
# ---------------------------------------------------------------------------
def neighbourhood_positive_fraction(
    x: np.ndarray,
    y: np.ndarray,
    positive: np.ndarray,
    population: Optional[np.ndarray] = None,
    radius_um: float = NEIGHBOURHOOD_RADIUS_UM,
    min_cells: int = NEIGHBOURHOOD_MIN_CELLS,
) -> Dict[str, object]:
    """Local marker-positive fraction in a fixed physical window around each cell.

    For every cell, counts the cells of ``population`` within ``radius_um`` and
    how many of them are ``positive``, and returns the ratio. Windows holding
    fewer than ``min_cells`` cells get the section-wide rate instead (the least
    informative value available).

    The window is in MICROMETRES on the measured coordinates, so this quantity is
    completely independent of the simulation grid spacing — which is exactly the
    property the oxygen-banded route lacks.

    Returns ``{'fraction', 'n_window', 'n_positive_window', 'section_rate',
    'n_unsupported', 'radius_um', 'min_cells'}``.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    positive = np.asarray(positive, dtype=bool)
    n = x.size
    if y.size != n or positive.size != n:
        raise ValueError("x, y and positive must be the same length")
    pop = np.ones(n, dtype=bool) if population is None else np.asarray(population, bool)
    if pop.shape != (n,):
        raise ValueError(f"population shape {pop.shape} != ({n},)")
    r = float(radius_um)
    if not r > 0.0:
        raise ValueError(f"radius_um must be > 0; got {r}")

    pts = np.column_stack([x, y])
    tree_pop = cKDTree(pts[pop])
    n_win = np.asarray(
        tree_pop.query_ball_point(pts, r=r, return_length=True, workers=-1),
        dtype=np.float64,
    )
    pos_pop = pop & positive
    if pos_pop.any():
        tree_pos = cKDTree(pts[pos_pop])
        n_pos = np.asarray(
            tree_pos.query_ball_point(pts, r=r, return_length=True, workers=-1),
            dtype=np.float64,
        )
    else:
        n_pos = np.zeros(n, dtype=np.float64)

    section_rate = float(positive[pop].mean()) if pop.any() else 0.0
    frac = np.full(n, section_rate, dtype=np.float64)
    ok = n_win >= int(min_cells)
    frac[ok] = n_pos[ok] / n_win[ok]
    return {
        "fraction": frac,
        "n_window": n_win,
        "n_positive_window": n_pos,
        "section_rate": section_rate,
        "n_unsupported": int((~ok).sum()),
        "radius_um": r,
        "min_cells": int(min_cells),
    }


# ---------------------------------------------------------------------------
# The state assignment
# ---------------------------------------------------------------------------
def state_from_markers(
    expression,
    genes: Sequence[str],
    x: Optional[np.ndarray] = None,
    y: Optional[np.ndarray] = None,
    cell_type: Optional[np.ndarray] = None,
    state: Optional[np.ndarray] = None,
    proliferation_markers: Sequence[str] = CELL_CYCLE_SPECIFIC_MARKERS,
    hypoxia_markers: Sequence[str] = HYPOXIA_MARKERS,
    include_low_specificity: bool = False,
    scoring_population: str = "tumor",
    radius_um: float = NEIGHBOURHOOD_RADIUS_UM,
    min_cells: int = NEIGHBOURHOOD_MIN_CELLS,
    per_cell_detection_floor: float = PER_CELL_DETECTION_FLOOR,
) -> Dict[str, object]:
    """Assign VIABLE / QUIESCENT / HYPOXIC from MEASURED markers, or say you can't.

    Parameters
    ----------
    expression : sparse matrix or ndarray, (n_cells, n_genes)
        Raw per-cell counts, row-aligned to ``x`` / ``y`` / ``cell_type``.
    genes : sequence of str
        Panel gene symbols aligned to ``expression``'s columns.
    x, y : ndarray (n_cells,), optional
        Measured centroids in um. REQUIRED whenever an axis falls to the
        neighbourhood call; if they are absent that axis is declared UNMEASURED
        rather than guessed.
    cell_type : ndarray (n_cells,), optional
        ``core.CellType`` codes. Used only to choose the scoring population.
    state : ndarray (n_cells,), optional
        Existing per-cell states. DEAD cells are never relabelled (no
        resurrection), matching ``cells.agents.update_cell_states``.
    proliferation_markers, hypoxia_markers : sequence of str
        Marker sets. Defaults are :data:`CELL_CYCLE_SPECIFIC_MARKERS` and
        :data:`HYPOXIA_MARKERS`.
    include_low_specificity : bool
        Fold :data:`LOW_SPECIFICITY_MARKERS` (PCNA) into the proliferation call.
        Default ``False`` — see the module docstring for why.
    scoring_population : {'tumor', 'all'}
        Which cells the detection rates and the neighbourhood denominators are
        computed over. ``'tumor'`` (default) needs ``cell_type``; it falls back
        to ``'all'`` with a recorded note if none is supplied.
    radius_um, min_cells
        Neighbourhood window; see :func:`neighbourhood_positive_fraction`.
    per_cell_detection_floor : float
        See :data:`PER_CELL_DETECTION_FLOOR`.

    Returns
    -------
    dict
        ``states``
            ``int8`` array of :class:`core.CellState` codes, length n_cells.
            Cells the panel cannot speak about stay VIABLE and the coverage block
            names the state that is unmeasured — read ``states`` only together
            with ``coverage``.
        ``coverage``
            Per-axis: requested / present / missing markers, per-gene detection,
            union detection, the CALL LEVEL and the reason for it, plus
            ``unmeasured_states``.
        ``fractions``
            State fractions over the scoring population and over tumour cells.
        ``per_cell``
            The intermediate arrays (positivity, neighbourhood fractions).
        ``provenance``
            One-paragraph statement of exactly what the labels mean, including
            the HYPOXIC-takes-precedence collapse and the uncalibrated-size
            caveat.

    Raises
    ------
    ValueError
        On a row-count mismatch between ``expression`` and the coordinates — a
        misalignment would attach one cell's transcripts to another's position,
        so it is a hard error rather than a warning.
    """
    genes_arr = np.asarray(genes, dtype=str)
    n = int(expression.shape[0])
    for name, arr in (("x", x), ("y", y), ("cell_type", cell_type), ("state", state)):
        if arr is not None and int(np.asarray(arr).size) != n:
            raise ValueError(
                f"{name} has {np.asarray(arr).size} entries but expression has {n} "
                "rows; the reconstruction and the expression matrix are not aligned"
            )

    # --- scoring population --------------------------------------------------
    pop_mode = str(scoring_population).strip().lower()
    pop_note = None
    if pop_mode == "tumor" and cell_type is not None:
        pop = np.asarray(cell_type) == int(CellType.TUMOR)
        pop_label = "TUMOR cells"
    else:
        if pop_mode == "tumor":
            pop_note = ("scoring_population='tumor' was requested but no cell_type "
                        "was supplied; scored over ALL cells instead")
        pop = np.ones(n, dtype=bool)
        pop_label = "all cells"
    if not pop.any():  # degenerate section — score everything rather than nothing
        pop = np.ones(n, dtype=bool)
        pop_label = "all cells (no cell of the requested type)"

    prol_set = list(proliferation_markers)
    if include_low_specificity:
        prol_set += [g for g in LOW_SPECIFICITY_MARKERS if g not in prol_set]

    prol = marker_detection(expression, genes_arr, prol_set, pop)
    hyp = marker_detection(expression, genes_arr, hypoxia_markers, pop)
    prol_level, prol_reason = _call_level(
        prol, "proliferation", float(per_cell_detection_floor)
    )
    hyp_level, hyp_reason = _call_level(
        hyp, "hypoxia", float(per_cell_detection_floor)
    )

    have_xy = x is not None and y is not None
    nbh: Dict[str, Optional[Dict[str, object]]] = {"proliferation": None,
                                                   "hypoxia": None}

    def _resolve(block, level, reason, axis):
        """Turn a call level into a boolean per-cell label (or ``None``)."""
        if level == "UNMEASURED":
            return None, level, reason
        pos = np.asarray(block["positive"], dtype=bool)
        if level == "per_cell":
            return pos, level, reason
        if not have_xy:
            return None, "UNMEASURED", (
                reason + " The neighbourhood aggregation needs measured cell "
                "coordinates and none were supplied, so this axis is UNMEASURED."
            )
        agg = neighbourhood_positive_fraction(
            x, y, pos, population=pop, radius_um=radius_um, min_cells=min_cells
        )
        nw = np.asarray(agg["n_window"])
        nbh[axis] = {
            k: v for k, v in agg.items()
            if k not in ("fraction", "n_window", "n_positive_window")
        }
        nbh[axis].update({  # type: ignore[union-attr]
            "fraction": agg["fraction"],
            "median_cells_per_window": float(np.median(nw)) if nw.size else 0.0,
            "frac_windows_below_floor": (
                float((nw < int(min_cells)).mean()) if nw.size else 0.0),
        })
        enriched = np.asarray(agg["fraction"]) > float(agg["section_rate"])
        return enriched, level, reason

    prol_call, prol_level, prol_reason = _resolve(
        prol, prol_level, prol_reason, "proliferation"
    )
    hyp_call, hyp_level, hyp_reason = _resolve(hyp, hyp_level, hyp_reason, "hypoxia")

    # --- assemble the single CellState label ---------------------------------
    states = np.full(n, int(CellState.VIABLE), dtype=np.int8)
    if prol_call is not None:
        states[~prol_call] = int(CellState.QUIESCENT)
    if hyp_call is not None:
        states[hyp_call] = int(CellState.HYPOXIC)  # HYPOXIC takes precedence
    if state is not None:  # never resurrect a DEAD cell
        prior = np.asarray(state)
        dead = prior == int(CellState.DEAD)
        states[dead] = int(CellState.DEAD)

    unmeasured: List[str] = []
    if prol_level == "UNMEASURED":
        unmeasured += ["VIABLE(proliferative)", "QUIESCENT"]
    if hyp_level == "UNMEASURED":
        unmeasured.append("HYPOXIC")

    coverage = {
        "n_panel_genes": int(genes_arr.size),
        "n_cells": n,
        "scoring_population": pop_label,
        "n_scored": int(pop.sum()),
        "scoring_population_note": pop_note,
        "per_cell_detection_floor": float(per_cell_detection_floor),
        "include_low_specificity": bool(include_low_specificity),
        "proliferation": {**{k: v for k, v in prol.items() if k != "positive"},
                          "call_level": prol_level, "reason": prol_reason,
                          "binarisation": _BINARISATION[prol_level],
                          "window": _window_summary(nbh["proliferation"])},
        "hypoxia": {**{k: v for k, v in hyp.items() if k != "positive"},
                    "call_level": hyp_level, "reason": hyp_reason,
                    "binarisation": _BINARISATION[hyp_level],
                    "window": _window_summary(nbh["hypoxia"])},
        "unmeasured_states": unmeasured,
        "neighbourhood_radius_um": float(radius_um),
        "neighbourhood_min_cells": int(min_cells),
    }

    tumor = (np.asarray(cell_type) == int(CellType.TUMOR)
             if cell_type is not None else pop)
    return {
        "states": states,
        "coverage": coverage,
        "fractions": {
            "scored": state_fractions_from_states(states, pop),
            "tumor": state_fractions_from_states(states, tumor),
            "all": state_fractions_from_states(states, None),
        },
        "per_cell": {
            "proliferation_positive": prol["positive"],
            "hypoxia_positive": hyp["positive"],
            "proliferation_call": prol_call,
            "hypoxia_call": hyp_call,
            "proliferation_neighbourhood": nbh["proliferation"],
            "hypoxia_neighbourhood": nbh["hypoxia"],
        },
        "provenance": _provenance(coverage),
        "kind": "MEASURED marker-derived cell state",
    }


def _provenance(coverage: Dict[str, object]) -> str:
    """The paragraph that has to travel with any number derived from these states."""
    p = coverage["proliferation"]  # type: ignore[index]
    h = coverage["hypoxia"]  # type: ignore[index]
    bits = [
        "MEASURED marker-derived cell states. "
        f"Proliferation: {p['call_level']} over {p['present'] or 'NOTHING'} "  # type: ignore[index]
        f"(union detection {float(p['union_detection']) * 100:.3f}%). "  # type: ignore[index]
        f"Hypoxia: {h['call_level']} over {h['present'] or 'NOTHING'} "  # type: ignore[index]
        f"(union detection {float(h['union_detection']) * 100:.3f}%)."  # type: ignore[index]
    ]
    if coverage["unmeasured_states"]:
        bits.append(
            "UNMEASURED on this panel: "
            f"{', '.join(coverage['unmeasured_states'])}. "  # type: ignore[index]
            "No cell was assigned an unmeasured state; those cells default to "
            "VIABLE and that default is a placeholder, NOT a finding."
        )
    if p["call_level"] == "neighbourhood":  # type: ignore[index]
        bits.append(
            "The proliferation label is a NEIGHBOURHOOD ENRICHMENT call over a "
            f"{coverage['neighbourhood_radius_um']:.0f} um window, not a per-cell "
            "Ki-67 call; its FRACTION is not a proliferative index."
        )
    if h["call_level"] != "UNMEASURED":  # type: ignore[index]
        bits.append(
            "The HYPOXIC label is positivity for a hypoxia-RESPONSE transcript, "
            "not an oxygen tension, and its SIZE is uncalibrated — these genes "
            "have real non-hypoxic expression, so the band over-calls hypoxia "
            "relative to any pO2-based definition."
        )
        if p["call_level"] != "UNMEASURED":  # type: ignore[index]
            bits.append(
                "Two axes were collapsed into one CellState label with HYPOXIC "
                "taking precedence over the proliferation call."
            )
    return " ".join(bits)


def state_fractions_from_states(
    states: np.ndarray, mask: Optional[np.ndarray] = None
) -> Dict[str, float]:
    """Fraction of ``states`` (optionally restricted to ``mask``) in each state.

    Same keys and convention as ``cells.agents.state_fractions`` so the routes are
    directly comparable, but it takes a plain array rather than a
    ``Microenvironment`` — the marker route produces states before any grid
    exists.
    """
    st = np.asarray(states)
    sel = np.ones(st.size, dtype=bool) if mask is None else np.asarray(mask, bool)
    names = {CellState.VIABLE: "viable", CellState.HYPOXIC: "hypoxic",
             CellState.QUIESCENT: "quiescent", CellState.DEAD: "dead"}
    k = int(sel.sum())
    if k == 0:
        return {v: 0.0 for v in names.values()}
    sub = st[sel]
    return {name: float((sub == int(s)).sum()) / k for s, name in names.items()}


def state_from_markers_micro(micro, expression, genes: Sequence[str], **kwargs):
    """:func:`state_from_markers` driven by a reconstructed ``Microenvironment``.

    Convenience wrapper that pulls ``x`` / ``y`` / ``cell_type`` / ``state`` off
    ``micro.cells``. ``expression`` rows must be index-aligned to ``micro.cells``
    — checked, never assumed. Does NOT modify ``micro``; the caller decides
    whether to install the returned states.
    """
    cells = micro.cells
    if int(expression.shape[0]) != int(cells.n):
        raise ValueError(
            f"expression rows ({expression.shape[0]}) != cells ({cells.n}); the "
            "reconstruction and the expression matrix are not aligned"
        )
    return state_from_markers(
        expression, genes, x=cells.x, y=cells.y, cell_type=cells.cell_type,
        state=cells.state, **kwargs,
    )


if __name__ == "__main__":  # pragma: no cover - exercised on the real cached sections
    from ..data.real.singlecell_spatial import load as _load

    print("=== cell_state self-test on the REAL cached Xenium sections ===")
    for path in ("data_cache/xenium_gbm", "data_cache/xenium_gbm_io",
                 "data_cache/gse324268_xenium"):
        d = _load(path=path, dense=False)
        xy = np.asarray(d["cell_xy_um"], dtype=np.float64)
        cov = panel_coverage(d["expression"], d["genes"])
        print(f"\n{path}  ({cov['n_panel_genes']} genes, {cov['n_cells']:,} cells)")
        for axis in ("proliferation", "hypoxia"):
            b = cov[axis]
            print(f"  {axis:14s} present={b['present']} missing={b['missing']}")
            print(f"  {'':14s} detection={ {k: round(v, 5) for k, v in b['detection'].items()} }")
            print(f"  {'':14s} -> {b['call_level']}")
        print(f"  UNMEASURED states: {cov['unmeasured_states'] or 'none'}")

        out = state_from_markers(
            d["expression"], d["genes"], x=xy[:, 0], y=xy[:, 1],
            scoring_population="all",
        )
        print(f"  fractions (all cells): "
              + ", ".join(f"{k}={v:.4f}" for k, v in out["fractions"]["all"].items()))
        assert out["states"].shape == (d["expression"].shape[0],)
        # Determinism.
        again = state_from_markers(d["expression"], d["genes"],
                                   x=xy[:, 0], y=xy[:, 1],
                                   scoring_population="all")
        assert np.array_equal(out["states"], again["states"])

    print("\ncell_state self-test PASSED")
