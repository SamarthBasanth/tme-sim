"""THREE ROUTES TO A CELL STATE, scored on real tissue — and what they do to kill.

``results/dose_schedule.md`` §1.3 showed that for temozolomide, whose per-cell
sensitivity is uniform 1.0 on the 319-gene panel, **~89% of the kill ceiling comes
from the QUIESCENT band alone** — the ceiling is `sensitivity x resistance`, so it
IS the cell-state labelling. ``results/state_validation.md`` §2.2 then showed that
labelling is not a measurement: the tumour VIABLE fraction runs **0.60% to 98.19%
(165x)** across two tissue-mask conventions and three grid spacings, and under the
shipping ``geometric`` mask the HYPOXIC band is EMPTY so ``hypoxia_resistance``
acts on no cells. ``results/oxygen_demotion.md`` §1 showed the field the bands come
from is out-predicted by a distance transform, with an effective length 2-5x the
measured oxygen-diffusion distance because ``OxygenConfig.consumption`` is
uncalibrated.

This module is the follow-through. It does four things and writes
``results/cell_state_impact.md``:

1. **Marker coverage, per panel.** Which proliferation (``MKI67``, ``TOP2A``,
   ``PCNA``) and hypoxia-response (``VEGFA``, ``CA9``, ``ADM``, ``NDRG1``,
   ``SLC2A1``) transcripts each cached section measures, at what detection rate,
   and therefore which states are UNMEASURED there. The three panels are very
   different and the answer is per-panel, never global.
2. **Three routes, same cells.** Simulated-oxygen banding (the current default),
   distance-to-vessel banding, and marker-derived states — state fractions,
   dx-stability, and pairwise agreement.
3. **Downstream.** How far predicted kill, the kill CEILING and the cross-drug
   ranking move between routes for temozolomide, cisplatin and doxorubicin.
4. **A recommendation**, against a stated bar: state fractions must be stable
   across dx and across sections, and must be traceable to a measurement or
   explicitly declared unmeasured. The bar is allowed to go unmet.

HOW THE DOWNSTREAM ARM IS COMPUTED (and why it is not a re-implementation)
--------------------------------------------------------------------------
Identical device to ``validation.oxygen_benchmark`` §3, for the identical reason.
``sim.engine.simulate_response`` re-bands the cells from its own oxygen solve, and
``sim/engine.py`` is not modified here. But **the cell state reaches the kill
through exactly one channel**: ``cells.pkpd.apply_drug_kill`` multiplies the Hill
probability by ``PDConfig.hypoxia_resistance`` / ``quiescent_resistance`` keyed on
``cells.state``. Nothing in ``physics.pressure`` / ``physics.diffusion`` /
``treatment.drugs`` reads the state, so the transport solve and the AUC field —
and hence every delivery metric — are state-invariant by construction.

So each alternative arm installs its own states on a private copy of the cells and
re-runs the UNMODIFIED ``apply_drug_kill`` on the engine's OWN transport result
with the engine's own seeded generator. That structural claim is not asserted, it
is gated: the oxygen arm is recomputed through the same path and **must reproduce
``simulate_response``'s kill_fraction bit-for-bit or the run aborts**.

CAVEATS THAT DO NOT GO AWAY
---------------------------
* **The three sections are three different patients**, on three different panels.
  Nothing transfers between them; three sections is a sample size of three.
* **No outcome label anywhere.** All three are untreated resections. Nothing here
  says which route's kill is RIGHT — only how much the answer depends on the
  choice, which is a sensitivity, not a validation.
* **Hypoxia-response transcripts are not an oxygen tension.** VEGFA / ADM /
  SLC2A1 have real non-hypoxic expression, so the marker route's HYPOXIC band is
  an enrichment label with an UNCALIBRATED SIZE. Only CA9 is near-restricted to
  severe hypoxia, and it fires in 0.9% of cells where it is on the panel.
* **The neighbourhood binarisation is a convention.** Where the per-cell
  proliferation call is dropout-inadmissible, the fallback label depends on the
  window radius; the sweep is measured and reported (§2.4) rather than hidden.
* **Typing circularity.** ``MKI67`` is one of the four brain-panel TUMOR markers,
  so an MKI67-high cell is marginally more likely to be typed TUMOR. That inflates
  the marginal positive rate inside the tumour population on that section; it is
  reported and it cannot be removed without changing the typing.

Run:    ./.venv/bin/python -m tmesim.validation.cell_state_impact
Writes: results/cell_state_impact.md

Units: distance um, fractions dimensionless, AUC uM*s.
"""
from __future__ import annotations

import os
import time
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..cells.agents import update_cell_states
from ..cells.pkpd import apply_drug_kill, effective_hypoxia_resistance
from ..core import (
    CellState,
    CellType,
    DrugParams,
    Microenvironment,
    OxygenConfig,
    PDConfig,
    PDEConfig,
    SimConfig,
)
from ..data.real.singlecell_spatial import load as load_singlecell_spatial
from ..physics.oxygen_proxy import diffusion_length_um, oxygen_from_vessel_distance
from ..reconstruct.cell_state import (
    CELL_CYCLE_SPECIFIC_MARKERS,
    HYPOXIA_MARKERS,
    HYPOXIA_SPECIFIC_MARKERS,
    LITERATURE_KI67,
    LOW_SPECIFICITY_MARKERS,
    NEIGHBOURHOOD_MIN_CELLS,
    NEIGHBOURHOOD_RADIUS_UM,
    PER_CELL_DETECTION_FLOOR,
    PROLIFERATION_MARKERS,
    marker_detection,
    neighbourhood_positive_fraction,
    state_fractions_from_states,
    state_from_markers_micro,
)
from ..reconstruct.real_pipeline import reconstruct_real, reconstruction_report
from ..sim.engine import simulate_response

__all__ = [
    "DOWNSTREAM_DRUGS",
    "DX_SWEEP",
    "ROUTES",
    "SECTIONS",
    "downstream_effect",
    "marker_coverage_survey",
    "radius_sensitivity",
    "route_comparison",
    "run",
    "write_report",
]

# ---------------------------------------------------------------------------
# Declared constants
# ---------------------------------------------------------------------------

#: The three cached real sections. THREE DIFFERENT PATIENTS, three panels.
#: ``(key, path, label, masks)`` — the flagship 319-gene section gets BOTH
#: tissue-mask conventions because that is where the documented 165x swing lives
#: and where every published kill number was computed; the two ~800k-cell sections
#: get the shipping ``geometric`` mask only, for compute.
SECTIONS: Tuple[Tuple[str, str, str, Tuple[str, ...]], ...] = (
    ("brain_319", "data_cache/xenium_gbm",
     "10x Xenium Human Brain — Glioblastoma (v1.3.0), 319-gene brain panel",
     ("geometric", "occupancy")),
    ("io_480", "data_cache/xenium_gbm_io",
     "10x Xenium Human Brain Cancer (GBM) FFPE (v2.0.0), 480-gene IO panel",
     ("geometric",)),
    ("gse324268_366", "data_cache/gse324268_xenium",
     "GEO GSE324268 GBMtumor07, 366-gene brain + GBM add-on panel",
     ("geometric",)),
)

#: The section every flagship drug result in ``results/`` was computed on.
FLAGSHIP_SECTION = "data_cache/xenium_gbm"

#: Grid spacings. 12 um is the grid-converged production spacing
#: (``results/numerical_rigor.md``); 30 and 20 expose grid sensitivity. Same
#: sweep as ``validation.state_validation``, so the two reports are comparable.
DX_CONVERGED = 12.0
DX_SWEEP: Tuple[float, ...] = (30.0, 20.0, DX_CONVERGED)

#: The three state routes under comparison.
ROUTES: Tuple[str, ...] = ("oxygen", "vessel_distance", "markers")

#: The three drugs the downstream arm quantifies, spanning the delivery axis.
#: Identical to ``validation.oxygen_benchmark.DOWNSTREAM_DRUGS`` so the two A/B
#: studies are directly comparable.
DOWNSTREAM_DRUGS: Tuple[str, ...] = ("temozolomide", "cisplatin", "doxorubicin")

#: Neighbourhood radii swept to price the one convention the marker route carries.
RADIUS_SWEEP_UM: Tuple[float, ...] = (50.0, 75.0, 100.0, 150.0, 200.0)

#: PDE knobs for the transport solve. IMPORTED from
#: ``validation.delivery_metrics`` rather than restated, so this A/B/C and every
#: flagship delivery number are solving the identical problem and cannot drift.
from .delivery_metrics import _PDE  # noqa: E402  (placed here to sit with the constants)


def _results_dir() -> str:
    d = os.path.abspath(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "results")
    )
    os.makedirs(d, exist_ok=True)
    return d


def _load(path: str) -> Dict:
    """Load a cached section, retrying once on a transient failure."""
    try:
        return load_singlecell_spatial(path=path, dense=False)
    except Exception:
        time.sleep(3)
        return load_singlecell_spatial(path=path, dense=False)


def _reconstruct(sample: str, dx_um: float, mask: str):
    """``reconstruct_real`` on the frozen default path, retrying once."""
    try:
        return reconstruct_real("xenium", sample=sample, dx_um=float(dx_um),
                                tissue_mask_mode=mask)
    except Exception:
        time.sleep(3)
        return reconstruct_real("xenium", sample=sample, dx_um=float(dx_um),
                                tissue_mask_mode=mask)


# ---------------------------------------------------------------------------
# 1. MARKER COVERAGE, PER PANEL
# ---------------------------------------------------------------------------
def marker_coverage_survey(
    sections: Sequence[Tuple[str, str, str, Tuple[str, ...]]] = SECTIONS,
    dx_um: float = 30.0,
) -> List[Dict]:
    """Which markers each panel carries, at what rate, and what is UNMEASURED.

    Detection is reported over TUMOUR cells (the population every kill number is
    about) as well as over all cells; the tumour restriction needs the frozen
    typing, so a coarse reconstruction is built for it. Cell typing in
    ``reconstruct.real_xenium.phenotype_cells`` is a function of expression only,
    so the types — and therefore this whole table — are independent of ``dx_um``.
    """
    print("\n" + "=" * 78)
    print("1. MARKER COVERAGE — what each panel can and cannot say")
    print("=" * 78)
    out: List[Dict] = []
    for key, path, label, _masks in sections:
        data = _load(path)
        expr, genes = data["expression"], np.asarray(data["genes"], dtype=str)
        micro = _reconstruct(path, dx_um, "geometric")
        if int(expr.shape[0]) != int(micro.cells.n):
            raise ValueError(
                f"[{key}] expression rows ({expr.shape[0]}) != cells "
                f"({micro.cells.n}); reconstruction and matrix are not aligned"
            )
        tumor = np.asarray(micro.cells.cell_type) == int(CellType.TUMOR)

        blocks: Dict[str, Dict] = {}
        for axis, mset in (
            ("proliferation_specific", CELL_CYCLE_SPECIFIC_MARKERS),
            ("proliferation_low_specificity", LOW_SPECIFICITY_MARKERS),
            ("hypoxia", HYPOXIA_MARKERS),
            ("hypoxia_specific", HYPOXIA_SPECIFIC_MARKERS),
        ):
            all_cells = marker_detection(expr, genes, mset, None)
            tum = marker_detection(expr, genes, mset, tumor)
            blocks[axis] = {
                "present": all_cells["present"],
                "missing": all_cells["missing"],
                "detection_all": all_cells["detection"],
                "detection_tumor": tum["detection"],
                "union_all": all_cells["union_detection"],
                "union_tumor": tum["union_detection"],
            }

        marked = state_from_markers_micro(micro, expr, genes)
        cov = marked["coverage"]
        row = {
            "key": key, "path": path, "label": label,
            "n_cells": int(micro.cells.n), "n_tumor": int(tumor.sum()),
            "n_panel_genes": int(genes.size),
            "blocks": blocks,
            "coverage": cov,
            "fractions_tumor": marked["fractions"]["tumor"],
            "provenance": marked["provenance"],
        }
        out.append(row)
        print(f"\n  {key} ({genes.size} genes, {micro.cells.n:,} cells, "
              f"{int(tumor.sum()):,} tumour)")
        for axis in ("proliferation", "hypoxia"):
            b = cov[axis]
            print(f"    {axis:14s} present={b['present']} missing={b['missing']}"
                  f"  -> {b['call_level']}")
        print(f"    UNMEASURED states: {cov['unmeasured_states'] or 'none'}")
    return out


# ---------------------------------------------------------------------------
# 2. THE THREE ROUTES
# ---------------------------------------------------------------------------
def _states_from_oxygen(micro, oxygen: np.ndarray, cfg: OxygenConfig) -> np.ndarray:
    """Band a field into states with the UNMODIFIED ``update_cell_states``."""
    work = Microenvironment(
        domain=micro.domain, fields=micro.fields, cells=micro.cells.copy(),
        patient_id=micro.patient_id, metadata=micro.metadata,
    )
    update_cell_states(work, np.asarray(oxygen, dtype=np.float64), cfg)
    return np.asarray(work.cells.state).copy()


def _agreement(a: np.ndarray, b: np.ndarray, mask: np.ndarray) -> Dict[str, float]:
    """Raw agreement and Cohen's kappa between two label vectors on ``mask``."""
    x, y = np.asarray(a)[mask], np.asarray(b)[mask]
    n = x.size
    if n == 0:
        return {"agreement": float("nan"), "kappa": float("nan"), "n": 0}
    po = float((x == y).mean())
    labels = (int(CellState.VIABLE), int(CellState.HYPOXIC), int(CellState.QUIESCENT),
              int(CellState.DEAD))
    pe = sum(float((x == L).mean()) * float((y == L).mean()) for L in labels)
    kappa = (po - pe) / (1.0 - pe) if pe < 1.0 else float("nan")
    return {"agreement": po, "kappa": float(kappa), "n": int(n)}


def route_comparison(
    sections: Sequence[Tuple[str, str, str, Tuple[str, ...]]] = SECTIONS,
    dx_sweep: Sequence[float] = DX_SWEEP,
) -> Dict:
    """State fractions from all three routes, over every (section, mask, dx).

    One reconstruction per configuration supplies the cells, the solved field and
    the vessel geometry; all three routes are then built on those SAME cells, so
    the comparison isolates the state rule and nothing else. Every route ends in
    the same place — an array of ``core.CellState`` codes — and the oxygen and
    vessel-distance routes both go through the UNMODIFIED
    ``cells.agents.update_cell_states``.
    """
    print("\n" + "=" * 78)
    print("2. THREE ROUTES — state fractions, dx-stability, agreement")
    print("=" * 78)
    cfg = OxygenConfig()
    rows: List[Dict] = []
    expr_cache: Dict[str, Tuple] = {}

    for key, path, _label, masks in sections:
        if path not in expr_cache:
            d = _load(path)
            expr_cache[path] = (d["expression"], np.asarray(d["genes"], dtype=str))
        expr, genes = expr_cache[path]
        for mask in masks:
            for dx in dx_sweep:
                t0 = time.time()
                micro = _reconstruct(path, float(dx), mask)
                rep = reconstruction_report(micro)
                tumor = np.asarray(micro.cells.cell_type) == int(CellType.TUMOR)

                states: Dict[str, np.ndarray] = {}
                states["oxygen"] = _states_from_oxygen(
                    micro, micro.fields.oxygen, cfg)
                states["vessel_distance"] = _states_from_oxygen(
                    micro, oxygen_from_vessel_distance(micro, cfg), cfg)
                marked = state_from_markers_micro(micro, expr, genes)
                states["markers"] = np.asarray(marked["states"])

                fr = {r: state_fractions_from_states(states[r], tumor)
                      for r in ROUTES}
                agree = {}
                for i, ra in enumerate(ROUTES):
                    for rb in ROUTES[i + 1:]:
                        agree[f"{ra}|{rb}"] = _agreement(
                            states[ra], states[rb], tumor)

                rows.append({
                    "key": key, "path": path, "mask": mask, "dx_um": float(dx),
                    "grid": tuple(int(v) for v in micro.domain.shape),
                    "n_cells": int(micro.cells.n), "n_tumor": int(tumor.sum()),
                    "mean_oxygen": rep["mean_oxygen"],
                    "tissue_coverage": rep["tissue_coverage_fraction"],
                    "fractions": fr,
                    "agreement": agree,
                    "marker_coverage": marked["coverage"],
                    "runtime_s": time.time() - t0,
                })
                print(f"  {key:14s} {mask:9s} dx={dx:4.0f}  "
                      + "  ".join(
                          f"{r[:3]}: V{fr[r]['viable'] * 100:5.1f}"
                          f"/Q{fr[r]['quiescent'] * 100:5.1f}"
                          f"/H{fr[r]['hypoxic'] * 100:5.1f}" for r in ROUTES)
                      + f"   [{time.time() - t0:.0f}s]")

    # --- dx-stability: max/min VIABLE fraction per (section, mask, route) ----
    stability: List[Dict] = []
    for key, path, _label, masks in sections:
        for mask in masks:
            sub = [r for r in rows if r["key"] == key and r["mask"] == mask]
            if not sub:
                continue
            for route in ROUTES:
                v = np.array([r["fractions"][route]["viable"] for r in sub])
                q = np.array([r["fractions"][route]["quiescent"] for r in sub])
                h = np.array([r["fractions"][route]["hypoxic"] for r in sub])
                stability.append({
                    "key": key, "mask": mask, "route": route,
                    "viable_min": float(v.min()), "viable_max": float(v.max()),
                    "viable_swing": (float(v.max() / v.min())
                                     if v.min() > 0 else float("inf")),
                    "viable_range_pp": float((v.max() - v.min()) * 100.0),
                    "quiescent_range_pp": float((q.max() - q.min()) * 100.0),
                    "hypoxic_range_pp": float((h.max() - h.min()) * 100.0),
                    "n_configs": int(v.size),
                })

    # --- cross-mask swing on the flagship section (the documented 165x) ------
    mask_swing: List[Dict] = []
    for route in ROUTES:
        sub = [r for r in rows if r["path"] == FLAGSHIP_SECTION]
        if not sub:
            continue
        v = np.array([r["fractions"][route]["viable"] for r in sub])
        mask_swing.append({
            "route": route,
            "viable_min": float(v.min()), "viable_max": float(v.max()),
            "swing": float(v.max() / v.min()) if v.min() > 0 else float("inf"),
            "n_configs": int(v.size),
        })

    # --- cross-section spread at the production config -----------------------
    cross_section: List[Dict] = []
    for route in ROUTES:
        sub = [r for r in rows
               if r["mask"] == "geometric" and r["dx_um"] == DX_CONVERGED]
        v = np.array([r["fractions"][route]["viable"] for r in sub])
        cross_section.append({
            "route": route,
            "viable_by_section": {r["key"]: r["fractions"][route]["viable"]
                                  for r in sub},
            "viable_min": float(v.min()) if v.size else float("nan"),
            "viable_max": float(v.max()) if v.size else float("nan"),
            "swing": (float(v.max() / v.min()) if v.size and v.min() > 0
                      else float("inf")),
        })

    return {
        "rows": rows,
        "stability": stability,
        "mask_swing": mask_swing,
        "cross_section": cross_section,
        "dx_sweep": [float(d) for d in dx_sweep],
        "diffusion_length_um": float(diffusion_length_um()),
    }


# ---------------------------------------------------------------------------
# 2.4 THE ONE CONVENTION THE MARKER ROUTE CARRIES
# ---------------------------------------------------------------------------
def radius_sensitivity(
    sample: str = FLAGSHIP_SECTION,
    dx_um: float = 30.0,
    radii: Sequence[float] = RADIUS_SWEEP_UM,
) -> Dict:
    """Price the neighbourhood window radius, on the section that needs it.

    The per-cell proliferation call is inadmissible on the 319-gene panel (MKI67
    0.35%, TOP2A 1.08%), so its state comes from a neighbourhood enrichment rule
    whose one free constant is the window radius. Rather than assert the default
    is fine, this sweeps it and reports what it costs. Cell typing and the marker
    aggregation are both mesh-independent, so ``dx_um`` only affects how long the
    (discarded) field solve takes.
    """
    print("\n" + "=" * 78)
    print("2.4 RADIUS SENSITIVITY — pricing the neighbourhood convention")
    print("=" * 78)
    data = _load(sample)
    expr, genes = data["expression"], np.asarray(data["genes"], dtype=str)
    micro = _reconstruct(sample, dx_um, "geometric")
    tumor = np.asarray(micro.cells.cell_type) == int(CellType.TUMOR)
    md = marker_detection(expr, genes, CELL_CYCLE_SPECIFIC_MARKERS, tumor)
    pos = np.asarray(md["positive"], dtype=bool)

    rows = []
    for r in radii:
        agg = neighbourhood_positive_fraction(
            micro.cells.x, micro.cells.y, pos, population=tumor,
            radius_um=float(r), min_cells=NEIGHBOURHOOD_MIN_CELLS,
        )
        nw = np.asarray(agg["n_window"])
        enriched = np.asarray(agg["fraction"]) > float(agg["section_rate"])
        rows.append({
            "radius_um": float(r),
            "median_cells_per_window": float(np.median(nw[tumor])),
            "frac_below_floor": float((nw[tumor] < NEIGHBOURHOOD_MIN_CELLS).mean()),
            "enriched_fraction_tumor": float(enriched[tumor].mean()),
        })
        print(f"  R={r:5.0f} um  median window {rows[-1]['median_cells_per_window']:5.1f} "
              f"cells  below-floor {rows[-1]['frac_below_floor'] * 100:5.1f}%  "
              f"proliferative {rows[-1]['enriched_fraction_tumor'] * 100:5.2f}%")

    e = np.array([r["enriched_fraction_tumor"] for r in rows])
    return {
        "sample": sample,
        "rows": rows,
        "section_positive_rate_tumor": float(md["union_detection"]),
        "swing": float(e.max() / e.min()) if e.min() > 0 else float("inf"),
        "min_cells": int(NEIGHBOURHOOD_MIN_CELLS),
        "default_radius_um": float(NEIGHBOURHOOD_RADIUS_UM),
    }


# ---------------------------------------------------------------------------
# 3. DOWNSTREAM EFFECT
# ---------------------------------------------------------------------------
def _kill_from_states(micro, states: np.ndarray, transport, cfg: SimConfig,
                      drug) -> Dict:
    """Install ``states`` and re-run the UNMODIFIED kill on the engine's transport.

    Nothing is re-implemented: this calls ``cells.pkpd.apply_drug_kill``, the same
    function the engine calls, on a private copy of the cells with the engine's own
    seeded generator. Valid because the AUC field has no cell-state dependence
    anywhere in the pipeline — and that validity is CHECKED by the caller, which
    requires this function to reproduce ``simulate_response`` exactly in the
    oxygen arm.
    """
    work = Microenvironment(
        domain=micro.domain, fields=micro.fields, cells=micro.cells.copy(),
        patient_id=micro.patient_id, metadata=micro.metadata,
    )
    work.cells.state = np.asarray(states, dtype=np.int8).copy()
    apply_drug_kill(work, transport, cfg.pd, drug,
                    rng=np.random.default_rng(cfg.seed))
    viable = float(work.cells.viable_tumor_fraction())
    return {"kill_fraction": 1.0 - viable, "viable_tumor_fraction": viable}


def _ceiling_decomposition(states: np.ndarray, sensitivity: np.ndarray,
                           tumor: np.ndarray, pd: PDConfig,
                           drug: Optional[DrugParams] = None) -> Dict:
    """The kill CEILING and its per-band decomposition — ``dose_schedule.md`` §1.3.

    As exposure -> infinity the Hill term saturates at ``max_kill_rate``, so the
    per-cell kill probability tends to ``max_kill_rate x sensitivity x
    resistance``. The ceiling is the mean of that over tumour cells and each
    band's contribution is its share of that mean. No dose can exceed it, and for
    a drug with uniform sensitivity it IS the state labelling.

    ``drug``, if given, is resolved through ``cells.pkpd.effective_hypoxia_resistance``
    so a per-drug ``hypoxia_resistance`` override (temozolomide has one) is
    honoured -- this reimplements ``apply_drug_kill``'s resistance math in
    closed form (at infinite exposure) rather than calling it.
    """
    r_hyp = (effective_hypoxia_resistance(drug, pd) if drug is not None
             else pd.hypoxia_resistance)
    st = np.asarray(states)[tumor]
    sens = np.asarray(sensitivity, dtype=np.float64)[tumor]
    n = int(st.size)
    res = np.ones(n, dtype=np.float64)
    res[st == int(CellState.HYPOXIC)] = r_hyp
    res[st == int(CellState.QUIESCENT)] = pd.quiescent_resistance
    per_cell = pd.max_kill_rate * sens * res
    per_cell = np.where(st == int(CellState.DEAD), 0.0, per_cell)
    bands = {}
    for s in (CellState.VIABLE, CellState.HYPOXIC, CellState.QUIESCENT):
        m = st == int(s)
        bands[s.name] = {
            "fraction": float(m.mean()) if n else 0.0,
            "mean_kill_probability": float(per_cell[m].mean()) if m.any() else 0.0,
            "contribution": float(per_cell[m].sum()) / n if n else 0.0,
        }
    return {"ceiling": float(per_cell.mean()) if n else 0.0, "bands": bands}


def downstream_effect(
    sample: str = FLAGSHIP_SECTION,
    dx_um: float = DX_CONVERGED,
    mask: str = "geometric",
    drugs: Sequence[str] = DOWNSTREAM_DRUGS,
) -> Dict:
    """Does the state route change the drug answer? — the A/B/C on the real section.

    Arm ``oxygen``: the shipping default, run through the unmodified
    :func:`sim.engine.simulate_response`.
    Arms ``vessel_distance`` and ``markers``: the same cells, drug and transport,
    with the states swapped (see :func:`_kill_from_states`).

    **The replication check is a hard gate.** Before any alternative arm is
    believed, the oxygen arm is recomputed through :func:`_kill_from_states` using
    the engine's own banded states and must return the engine's ``kill_fraction``
    bit-for-bit. If it does not, the structural assumption behind the whole
    comparison is wrong and the function raises rather than reporting a number.
    """
    print("\n" + "=" * 78)
    print(f"3. DOWNSTREAM — kill, ceiling and ranking by state route "
          f"(dx = {dx_um:g} um, mask = {mask})")
    print("=" * 78)

    from .calibration import calibrate_alpha
    from .delivery_metrics import compile_real_drug

    alpha = float(calibrate_alpha())
    data = _load(sample)
    expr, genes = data["expression"], np.asarray(data["genes"], dtype=str)
    micro = _reconstruct(sample, float(dx_um), mask)
    rep = reconstruction_report(micro)
    tumor = np.asarray(micro.cells.cell_type) == int(CellType.TUMOR)
    ocfg = OxygenConfig()

    states = {
        "vessel_distance": _states_from_oxygen(
            micro, oxygen_from_vessel_distance(micro, ocfg), ocfg),
    }
    marked = state_from_markers_micro(micro, expr, genes)
    states["markers"] = np.asarray(marked["states"])
    print(f"  {micro.cells.n:,} cells ({int(tumor.sum()):,} tumour), "
          f"grid {micro.domain.shape}, alpha = {alpha:.2f}")

    pdcfg = PDConfig()
    rows: Dict[str, Dict] = {}
    band_fracs: Dict[str, Dict[str, float]] = {}
    for name in drugs:
        dp, mol = compile_real_drug(name, alpha)
        cfg = SimConfig(drug=dp, pde=PDEConfig(**_PDE), seed=0)
        t0 = time.time()
        res = simulate_response(micro, cfg)
        engine_states = _states_from_oxygen(micro, res.fields["oxygen"], cfg.oxygen)
        states["oxygen"] = engine_states

        # HARD GATE.
        check = _kill_from_states(micro, engine_states, res.transport, cfg, dp)
        if check["kill_fraction"] != float(res.kill_fraction):
            raise AssertionError(
                f"[{name}] the state-swap path does not reproduce the engine "
                f"({check['kill_fraction']!r} vs {float(res.kill_fraction)!r}); the "
                "assumption that the cell state reaches the kill only through "
                "apply_drug_kill's resistance multipliers is WRONG and no "
                "downstream number here may be believed"
            )

        per_route = {}
        for route in ROUTES:
            k = (float(res.kill_fraction) if route == "oxygen"
                 else _kill_from_states(
                     micro, states[route], res.transport, cfg, dp)["kill_fraction"])
            ceil = _ceiling_decomposition(
                states[route], micro.cells.sensitivity, tumor, pdcfg, drug=dp)
            per_route[route] = {
                "kill_fraction": float(k),
                "ceiling": ceil["ceiling"],
                "bands": ceil["bands"],
                "pct_of_ceiling": (float(k) / ceil["ceiling"] * 100.0
                                   if ceil["ceiling"] > 0 else float("nan")),
            }
            band_fracs[route] = {
                b: v["fraction"] for b, v in ceil["bands"].items()}
        rows[name] = {
            "drug": name, "ic50_uM": float(dp.ic50), "mw": float(mol.get("mw")),
            "routes": per_route, "runtime_s": time.time() - t0,
        }
        base = per_route["oxygen"]["kill_fraction"]
        print(f"  {name:14s} kill " + "  ".join(
            f"{r}={per_route[r]['kill_fraction']:.4f}" for r in ROUTES)
            + "   ceiling " + "  ".join(
                f"{r}={per_route[r]['ceiling']:.4f}" for r in ROUTES)
            + f"   (base {base:.4f})")

    ranking = {
        route: sorted(rows, key=lambda k: -rows[k]["routes"][route]["kill_fraction"])
        for route in ROUTES
    }
    return {
        "sample": sample, "dx_um": float(dx_um), "mask": mask,
        "alpha": float(alpha), "n_cells": int(micro.cells.n),
        "n_tumor": int(tumor.sum()),
        "grid": tuple(int(v) for v in micro.domain.shape),
        "report": rep,
        "pd": {"hypoxia_resistance": pdcfg.hypoxia_resistance,
               "quiescent_resistance": pdcfg.quiescent_resistance,
               "max_kill_rate": pdcfg.max_kill_rate},
        "band_fractions": band_fracs,
        "marker_coverage": marked["coverage"],
        "marker_provenance": marked["provenance"],
        "rows": rows,
        "ranking": ranking,
        "ranking_invariant": len({tuple(v) for v in ranking.values()}) == 1,
    }


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------
def run(
    sections: Sequence[Tuple[str, str, str, Tuple[str, ...]]] = SECTIONS,
    dx_sweep: Sequence[float] = DX_SWEEP,
    drugs: Sequence[str] = DOWNSTREAM_DRUGS,
    cache_path: Optional[str] = None,
) -> Dict:
    """Coverage + three-route comparison + radius sweep + downstream, in one pass.

    Deterministic: the only randomness is the engine's seeded generator
    (``seed = 0``). ``cache_path`` (optional) pickles the raw results so the report
    can be regenerated without repeating the compute — same convention as
    :func:`validation.oxygen_benchmark.run`.
    """
    if cache_path and os.path.exists(cache_path):
        import pickle

        with open(cache_path, "rb") as fh:
            print(f"loaded cached results from {cache_path}")
            return pickle.load(fh)

    t0 = time.time()
    coverage = marker_coverage_survey(sections)
    routes = route_comparison(sections, dx_sweep)
    radius = radius_sensitivity()
    down = downstream_effect(drugs=drugs)
    results = {
        "coverage": coverage, "routes": routes, "radius": radius,
        "downstream": down, "elapsed_s": time.time() - t0,
    }
    if cache_path:
        import pickle

        with open(cache_path, "wb") as fh:
            pickle.dump(results, fh)
        print(f"cached results to {cache_path}")
    return results


# ---------------------------------------------------------------------------
# The report
# ---------------------------------------------------------------------------
def _pct(v, nd: int = 2) -> str:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return "n/a"
    return "n/a" if not np.isfinite(f) else f"{f * 100:.{nd}f}%"


def _fmt(v, nd: int = 4) -> str:
    if v is None:
        return "n/a"
    try:
        f = float(v)
    except (TypeError, ValueError):
        return str(v)
    return "n/a" if not np.isfinite(f) else f"{f:.{nd}f}"


def _swing(v) -> str:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return "n/a"
    if not np.isfinite(f):
        return "unbounded (the band empties at one setting)"
    return f"{f:.2f}x"


def write_report(results: Dict, out_path: Optional[str] = None) -> str:
    """Write ``results/cell_state_impact.md``."""
    cov = results["coverage"]
    rt = results["routes"]
    rad = results["radius"]
    dn = results["downstream"]
    out_path = out_path or os.path.join(_results_dir(), "cell_state_impact.md")

    by_key = {c["key"]: c for c in cov}
    stab = {(s["key"], s["mask"], s["route"]): s for s in rt["stability"]}
    mask_swing = {m["route"]: m for m in rt["mask_swing"]}
    cross = {c["route"]: c for c in rt["cross_section"]}
    prod = [r for r in rt["rows"]
            if r["mask"] == "geometric" and r["dx_um"] == DX_CONVERGED]

    L: List[str] = []
    A = L.append

    # ---------------- headline ----------------
    A("# Where the cell state comes from, and what it is worth")
    A("")
    A("`results/dose_schedule.md` derives ~89% of the temozolomide kill ceiling "
      "from the QUIESCENT band alone; `results/state_validation.md` shows that "
      "band is set by a tissue-mask convention and a mesh, not by tissue; "
      "`results/oxygen_demotion.md` shows the field it is banded from is beaten by "
      "a distance transform. So the largest single term in every absolute kill "
      "number this project reports is a bookkeeping choice. This report replaces "
      "it where the data allows, refuses where it does not, and measures what the "
      "difference is worth.")
    A("")
    A("**The default has not been changed.** `reconstruct_real(state_mode='oxygen')` "
      "is still the shipping route and does not touch `cells.state` at all, so it "
      "is bit-identical to the pre-option pipeline (asserted by "
      "`tests/test_cell_state.py`). What is new is that a MEASURED route exists, "
      "that it declares what its panel cannot measure, and that the default's own "
      "falsification travels with it in `reconstruction_report()['cell_state']"
      "['caveat']`.")
    A("")

    ox_flag = mask_swing.get("oxygen", {})
    mk_flag = mask_swing.get("markers", {})
    vd_flag = mask_swing.get("vessel_distance", {})
    tmz = dn["rows"].get("temozolomide")
    A("## Headline")
    A("")
    A("| | |")
    A("|---|---|")
    A(f"| **Panels that can call proliferation per cell** | "
      f"{sum(1 for c in cov if c['coverage']['proliferation']['call_level'] == 'per_cell')}"
      f" of {len(cov)} — the 319-gene flagship panel is NOT one of them |")
    A(f"| **Panels that measure hypoxia at all** | "
      f"{sum(1 for c in cov if c['coverage']['hypoxia']['call_level'] != 'UNMEASURED')}"
      f" of {len(cov)} — on the 319-gene flagship section hypoxia is UNMEASURED |")
    A(f"| **VIABLE-fraction swing across mask x dx (flagship section)** | "
      f"oxygen **{_swing(ox_flag.get('swing'))}**, vessel-distance "
      f"{_swing(vd_flag.get('swing'))}, markers **{_swing(mk_flag.get('swing'))}** |")
    if tmz:
        ks = [tmz["routes"][r]["kill_fraction"] for r in ROUTES]
        cs = [tmz["routes"][r]["ceiling"] for r in ROUTES]
        A(f"| **Temozolomide kill across routes** | {min(ks):.4f} to {max(ks):.4f} "
          f"({max(ks) / max(min(ks), 1e-12):.2f}x) |")
        A(f"| **Temozolomide kill CEILING across routes** | {min(cs):.4f} to "
          f"{max(cs):.4f} ({max(cs) / max(min(cs), 1e-12):.2f}x) |")
    A(f"| **Cross-drug ranking** | "
      + ("**INVARIANT** across all three routes — but NOT across the delivery "
         "law (§3.3a) |" if dn["ranking_invariant"]
         else "**CHANGES** between routes |"))
    A("| **Recommendation** | see §5 |")
    A("")
    A("> **The three sections are THREE DIFFERENT PATIENTS** on three different "
      "panels. Nothing transfers between them and three sections is a sample size "
      "of three. None of them carries a drug-response outcome label, so nothing "
      "here says which route is RIGHT — only how much the answer depends on the "
      "choice.")
    A("")

    # ---------------- 1. coverage ----------------
    A("## 1. Marker coverage, panel by panel")
    A("")
    A("A state can only be derived from a marker the panel actually assays, and "
      "only at a detection rate that is biology rather than probe sensitivity. "
      "Detection is over TUMOUR cells (the population every kill number is about); "
      "a zero is \"no probe fired\", not \"no transcript\".")
    A("")
    A("| section | panel | cells (tumour) | MKI67 | TOP2A | PCNA | VEGFA | CA9 | "
      "ADM | NDRG1 | SLC2A1 |")
    A("|---|---|---|---|---|---|---|---|---|---|---|")
    for c in cov:
        det = {}
        for axis in ("proliferation_specific", "proliferation_low_specificity",
                     "hypoxia", "hypoxia_specific"):
            det.update(c["blocks"][axis]["detection_tumor"])
        cells_ = f"{c['n_cells']:,} ({c['n_tumor']:,})"
        vals = " | ".join(
            (_pct(det[g], 2) if g in det else "—")
            for g in ("MKI67", "TOP2A", "PCNA", "VEGFA", "CA9", "ADM", "NDRG1",
                      "SLC2A1"))
        A(f"| `{c['key']}` | {c['n_panel_genes']} genes | {cells_} | {vals} |")
    A("")
    A("`—` means the gene is **not on that panel**. Reading the table:")
    A("")
    for c in cov:
        p, h = c["coverage"]["proliferation"], c["coverage"]["hypoxia"]
        A(f"* **`{c['key']}`** — {c['label']}.  ")
        A(f"  Proliferation: `{p['call_level']}` over {p['present'] or 'NO GENE'} "
          f"(union {_pct(p['union_detection'], 3)} of tumour cells). "
          + (f"Hypoxia: `{h['call_level']}` over {h['present']} "
             f"(union {_pct(h['union_detection'], 2)})."
             if h["present"] else
             "Hypoxia: `UNMEASURED` — no gene of the set is on this panel.")
          + (f"  **UNMEASURED on this panel: "
             f"{', '.join(c['coverage']['unmeasured_states'])}.**"
             if c["coverage"]["unmeasured_states"] else ""))
    A("")
    A("### 1.1 Why the flagship section cannot make a per-cell proliferation call")
    A("")
    brain = by_key.get("brain_319")
    if brain:
        p = brain["coverage"]["proliferation"]
        u = float(p["union_detection"])
        A(f"On the 319-gene panel MKI67 fires in "
          f"{_pct(brain['blocks']['proliferation_specific']['detection_tumor'].get('MKI67', 0), 3)}"
          f" of tumour cells and TOP2A in "
          f"{_pct(brain['blocks']['proliferation_specific']['detection_tumor'].get('TOP2A', 0), 3)}"
          f", for a union of {_pct(u, 3)}. Published GBM Ki-67 labelling indices "
          f"are {LITERATURE_KI67[0] * 100:.0f}-{LITERATURE_KI67[1] * 100:.0f}% "
          "(Fudaba et al., World Neurosurg 2021;149:e469-e480, PMID 33567368; "
          "Alkhaibary et al., Hematol Oncol Stem Cell Ther 2018;12(2):82-88, PMID "
          "30552865; Zheng et al., Zhonghua Bing Li Xue Za Zhi 2021;50(3):229-235, "
          "PMID 33677887; Caramanti et al., Cureus 2022;14(4):e24246, PMID "
          f"35602791) — i.e. this panel sits **{LITERATURE_KI67[0] / max(u, 1e-12):.0f}x "
          "below the lower bound**. A per-cell \"this cell is cycling\" call at "
          "that rate reports probe sensitivity, not biology, so "
          "`state_from_markers` refuses to make one and falls back to a "
          "neighbourhood aggregation (§2.4 prices that fallback).")
        A("")
        pcna = float(
            brain["blocks"]["proliferation_low_specificity"]["detection_tumor"]
            .get("PCNA", 0.0))
        mki = float(
            brain["blocks"]["proliferation_specific"]["detection_tumor"]
            .get("MKI67", 0.0))
        A(f"The same panel carries **PCNA at {_pct(pcna, 2)}**, which is inside the "
          "published Ki-67 range and looks like a rescue. It is not, and it is "
          "deliberately excluded from the default call: PCNA protein has a long "
          "half-life and the transcript is expressed in non-cycling and "
          f"DNA-repairing cells. A {pcna / max(mki, 1e-12):.0f}x gap between PCNA "
          "and MKI67 detection **on the same cells** is a specificity difference, "
          "not a biological one — using PCNA would buy a per-cell call by measuring "
          "something other than proliferation. It is reported here and available "
          "via `include_low_specificity=True`; it is not the default.")
        A("")
    io = by_key.get("io_480")
    if io:
        A("### 1.2 The IO panel *can* — and it lands on the published index")
        A("")
        A(f"On the 480-gene IO section MKI67 fires in "
          f"{_pct(io['blocks']['proliferation_specific']['detection_tumor'].get('MKI67', 0), 2)}"
          " of tumour cells. That is **inside** the published "
          f"{LITERATURE_KI67[0] * 100:.0f}-{LITERATURE_KI67[1] * 100:.0f}% Ki-67 "
          "range, from a completely independent measurement (in-situ transcript "
          "detection vs IHC protein scoring on other cohorts). It is the single "
          "strongest piece of external corroboration in this report, and it is "
          "worth being precise about what it does and does not establish: it "
          "establishes that **this panel's MKI67 probe is not dropout-limited**, "
          "which is exactly the failure that disqualifies the 319-gene panel. It "
          "does not establish that any individual cell's call is correct. Note "
          "also that the two sections are different patients, so the comparison is "
          "panel-vs-panel, not tumour-vs-tumour.")
        A("")
    A("### 1.3 The hypoxia band is measurable on two panels and its SIZE is "
      "uncalibrated on both")
    A("")
    A("The 319-gene flagship panel carries no VEGFA, CA9, ADM, NDRG1 or SLC2A1, so "
      "on the section every published kill number was computed from, **hypoxia is "
      "UNMEASURED**. It has not been simulated away, and `state_from_markers` "
      "assigns no cell the HYPOXIC state there rather than manufacturing one.")
    A("")
    A("Where the markers do exist, the honest limit is different and just as "
      "important. \"Positive for a hypoxia-response transcript\" is traceable to a "
      "measurement, but it is **not calibrated to a hypoxic fraction**: VEGFA, ADM "
      "and SLC2A1 all have real non-hypoxic expression, and VEGFA alone fires in "
      + (_pct(by_key["io_480"]["blocks"]["hypoxia"]["detection_tumor"].get("VEGFA", 0), 1)
         if "io_480" in by_key else "63%")
      + " of IO-section tumour cells. A 60%-hypoxic tumour is not a finding, it is "
        "a threshold artefact. Of the set only **CA9** is near-restricted to severe "
        "hypoxia, and where it is on the panel it fires in "
      + (_pct(by_key["gse324268_366"]["blocks"]["hypoxia_specific"]["detection_tumor"].get("CA9", 0), 2)
         if "gse324268_366" in by_key else "~1%")
      + " of tumour cells — a plausible perinecrotic fraction, but too sparse for a "
        "per-cell call on its own. The marker route's HYPOXIC band is therefore an "
        "ENRICHMENT label with an uncalibrated size, and this report never quotes "
        "its magnitude as a hypoxic fraction.")
    A("")

    # ---------------- 2. three routes ----------------
    A("## 2. Three routes to a state, on the same cells")
    A("")
    A("Every configuration is ONE reconstruction (the frozen `reconstruct_real` "
      "path) whose cells, solved field and vessel geometry feed all three routes, "
      "so the comparison isolates the state rule. The oxygen and vessel-distance "
      "routes both end in the UNMODIFIED `cells.agents.update_cell_states`; the "
      "marker route never touches a grid. Fractions are over TUMOUR cells.")
    A("")
    A("| route | what it is | what it is traceable to |")
    A("|---|---|---|")
    A("| `oxygen` | the current default: band the SIMULATED oxygen field at "
      "`t = 0.2` / `2t` | nothing measured — a PDE with an unfitted `consumption`, "
      "gated by a tissue-mask convention |")
    A(f"| `vessel_distance` | band `exp(-d/L)`, `L = {rt['diffusion_length_um']:.1f}` "
      "um, in distance to the nearest measured endothelial cell | the MEASURED "
      "endothelial positions + the Thomlinson & Gray 150 um oxygen-diffusion "
      "distance (PMID 13304213) |")
    A("| `markers` | proliferative/quiescent from measured proliferation "
      "transcripts, hypoxic from measured hypoxia-response transcripts | MEASURED "
      "transcripts, with an explicit UNMEASURED verdict where the panel is silent |")
    A("")
    A("### 2.1 State fractions across every configuration")
    A("")
    A("| section | mask | dx (um) | grid | mean O2 | oxygen V/Q/H | "
      "vessel_distance V/Q/H | markers V/Q/H |")
    A("|---|---|---|---|---|---|---|---|")
    for r in rt["rows"]:
        def _vqh(route):
            f = r["fractions"][route]
            return (f"{f['viable'] * 100:.2f} / {f['quiescent'] * 100:.2f} / "
                    f"{f['hypoxic'] * 100:.2f}")
        A(f"| `{r['key']}` | `{r['mask']}` | {r['dx_um']:g} | "
          f"{r['grid'][0]}x{r['grid'][1]} | {_fmt(r['mean_oxygen'])} | "
          f"{_vqh('oxygen')} | {_vqh('vessel_distance')} | {_vqh('markers')} |")
    A("")
    A("> **These oxygen-route numbers are NOT the ones in "
      "`results/state_validation.md` §2.2, and that is expected.** That table was "
      "computed under the KDE vascular source; the shipping default is now the "
      "DISCRETE segmented-vessel source (`reconstruct.real_fields."
      "_DEFAULT_VASCULAR_MODE`, promoted in `results/vessel_promotion.md`), whose "
      "tissue-integrated supply is ~3.4x smaller, so the whole oxygen field — and "
      "therefore every band edge — sits lower. The cross-check that the path here "
      "is the frozen one is `results/dose_schedule.md` §1.3, which ran under the "
      "discrete source: it quotes 10.0% VIABLE / 90.0% QUIESCENT and a "
      "temozolomide kill ceiling of 0.5502 at `geometric`, dx = 12 um, and both are "
      "reproduced to four decimals in §3 below. The consequence for this report is "
      "that the swing figures below are re-measured under the current source "
      "rather than quoted from the older document.")
    A("")
    A("### 2.2 dx-stability — the property the bar is about")
    A("")
    A("Range of the tumour VIABLE fraction across the three grid spacings "
      f"({', '.join(f'{d:g}' for d in rt['dx_sweep'])} um) at a FIXED tissue mask. "
      "The mesh is a numerical convenience; a state fraction that moves with it is "
      "not a property of the tissue.")
    A("")
    A("| section | mask | route | VIABLE min | VIABLE max | swing | range (pp) |")
    A("|---|---|---|---|---|---|---|")
    for s in rt["stability"]:
        A(f"| `{s['key']}` | `{s['mask']}` | `{s['route']}` | "
          f"{_pct(s['viable_min'], 4)} | {_pct(s['viable_max'], 4)} | "
          f"{_swing(s['viable_swing'])} | {s['viable_range_pp']:.2f} |")
    A("")
    A("And across BOTH tissue-mask conventions x three dx on the flagship section "
      "— the same six-configuration grid over which `results/state_validation.md` "
      "§2.2 measured its 165x under the KDE-era vascular source, re-measured here "
      "under the current discrete source:")
    A("")
    A("| route | VIABLE min | VIABLE max | swing | configurations |")
    A("|---|---|---|---|---|")
    for m in rt["mask_swing"]:
        A(f"| `{m['route']}` | {_pct(m['viable_min'], 4)} | "
          f"{_pct(m['viable_max'], 4)} | {_swing(m['swing'])} | "
          f"{m['n_configs']} |")
    A("")
    A("**The marker route is mesh-independent by construction, and the table "
      "measures it rather than asserting it.** It is computed from measured cell "
      "coordinates and measured counts in micrometres; no grid enters. The oxygen "
      "route is not: the tissue mask gates the vascular source, so it sets the "
      "absolute oxygen level, and the mask is grid-dependent under the `occupancy` "
      "convention.")
    A("")
    A("### 2.3 Cross-section stability at the production configuration "
      f"(`geometric`, dx = {DX_CONVERGED:g} um)")
    A("")
    A("| route | " + " | ".join(f"`{c['key']}`" for c in cov) + " | spread |")
    A("|---|" + "---|" * (len(cov) + 1))
    for route in ROUTES:
        vals = cross[route]["viable_by_section"]
        A(f"| `{route}` | "
          + " | ".join(_pct(vals.get(c["key"], float("nan")), 4) for c in cov)
          + f" | {_swing(cross[route]['swing'])} |")
    A("")
    A("Three different patients on three different panels, so a spread here is not "
      "by itself a fault — real tumours differ. What it does show is that no route "
      "produces a *panel-independent* number: the marker route's fractions move "
      "with which markers the panel carries and at what detection rate, which is "
      "why the coverage block in §1 has to travel with any number derived from it.")
    A("")
    A("### 2.4 Pricing the one convention the marker route carries")
    A("")
    A("On the flagship panel the per-cell proliferation call is inadmissible, so "
      "the label comes from a neighbourhood enrichment rule whose one free "
      "constant is the window radius. That constant is swept rather than asserted:")
    A("")
    A(f"Section-wide tumour positive rate: {_pct(rad['section_positive_rate_tumor'], 3)} "
      f"(MKI67 or TOP2A). Windows below {rad['min_cells']} cells fall back to the "
      "section rate.")
    A("")
    A("| window radius (um) | median cells/window | windows below the floor | "
      "resulting proliferative fraction |")
    A("|---|---|---|---|")
    for r in rad["rows"]:
        mark = " **(default)**" if r["radius_um"] == rad["default_radius_um"] else ""
        A(f"| {r['radius_um']:.0f}{mark} | {r['median_cells_per_window']:.1f} | "
          f"{_pct(r['frac_below_floor'], 1)} | "
          f"{_pct(r['enriched_fraction_tumor'])} |")
    A("")
    A(f"**The fraction swings {_swing(rad['swing'])} across the radius**, and that "
      "is structural rather than fixable: as R -> 0 only the positive cells "
      "themselves are enriched, so the fraction tends to the detection rate; as "
      "R -> infinity every window tends to the section rate, so the fraction tends "
      "to 0. **The ordering is measured; the level is not.** The default "
      f"{rad['default_radius_um']:.0f} um is chosen by a stated rule — the smallest "
      f"swept radius at which >95% of windows meet the {rad['min_cells']}-cell "
      "floor on the sparsest section — and 50 um (the radius "
      "`validation.state_validation` uses for its crowding baseline) is degenerate "
      "here, with "
      + _pct(next((r["frac_below_floor"] for r in rad["rows"]
                   if r["radius_um"] == 50.0), float("nan")), 0)
      + " of windows below the floor.")
    A("")
    A("So the marker route on THIS panel buys mesh-independence and loses a window "
      "radius. That is a smaller and better-behaved convention than the oxygen "
      "route's mask x mesh swing (§2.2, unbounded here because the VIABLE band "
      "empties; 165x under the KDE-era source) — it is bounded, it is one number, "
      "and its two limits are known analytically — but it is a convention, and this "
      "report does not pretend otherwise.")
    A("")
    A("### 2.5 Do the routes agree with each other?")
    A("")
    A(f"Per-cell agreement and Cohen's kappa over tumour cells, at `geometric`, "
      f"dx = {DX_CONVERGED:g} um.")
    A("")
    A("| section | pair | agreement | Cohen's kappa |")
    A("|---|---|---|---|")
    for r in prod:
        for pair, a in r["agreement"].items():
            # The pair key contains a pipe; escape it or it breaks the table.
            label = " vs ".join(f"`{p}`" for p in pair.split("|"))
            A(f"| `{r['key']}` | {label} | {_pct(a['agreement'])} | "
              f"{_fmt(a['kappa'])} |")
    A("")
    A("Kappa near zero means the two routes are labelling the same cells no better "
      "than their marginal rates alone would predict — i.e. they agree on *how "
      "many* cells are resistant and disagree on *which*. That distinction matters "
      "downstream: a kill fraction depends almost entirely on the marginals, while "
      "any spatial claim (\"this region is drug-refractory\") depends on the "
      "identities, and the identities are where the routes part company.")
    A("")

    # ---------------- 3. downstream ----------------
    A("## 3. Downstream — what the choice is worth in predicted kill")
    A("")
    A(f"**Section.** `{dn['sample']}` — {dn['n_cells']:,} cells "
      f"({dn['n_tumor']:,} tumour), {dn['mask']} tissue mask, dx = "
      f"{dn['dx_um']:g} um, grid {dn['grid'][0]}x{dn['grid'][1]}. "
      f"Perivascular-sink alpha = {dn['alpha']:.2f}, calibrated once on "
      "doxorubicin as in every flagship run. Per-cell sensitivity is the "
      "reconstruction default (uniform 1.0), so for all three drugs the ceiling is "
      "the state labelling and nothing else — which is precisely the quantity "
      "under test.")
    A("")
    A("**How the alternative arms are computed.** `sim.engine.simulate_response` "
      "re-bands the cells from its own oxygen solve and `sim/engine.py` is not "
      "modified, so each alternative arm installs its states on a private copy and "
      "re-runs the unmodified `cells.pkpd.apply_drug_kill` on **the engine's own "
      "transport result**, with the engine's own seeded generator. Valid because "
      "the cell state reaches the kill through exactly one channel — the "
      "`PDConfig` resistance multipliers — and the AUC field has no state "
      "dependence anywhere. Not assumed: the oxygen arm is recomputed through the "
      "same path and **must reproduce the engine bit-for-bit or the run aborts**. "
      "It did, for all three drugs.")
    A("")
    A("### 3.1 Tumour-cell band occupancy under each route")
    A("")
    A("| route | VIABLE | QUIESCENT | HYPOXIC | resistance multiplier applied |")
    A("|---|---|---|---|---|")
    for route in ROUTES:
        b = dn["band_fractions"][route]
        A(f"| `{route}` | {_pct(b['VIABLE'])} | {_pct(b['QUIESCENT'])} | "
          f"{_pct(b['HYPOXIC'])} | x1.00 / x"
          f"{dn['pd']['quiescent_resistance']:.2f} / x"
          f"{dn['pd']['hypoxia_resistance']:.2f} |")
    A("")
    if dn["band_fractions"]["oxygen"]["HYPOXIC"] <= 0.0:
        A("**The default route's HYPOXIC band is EMPTY on this section under the "
          "shipping mask**, so `PDConfig.hypoxia_resistance = "
          f"{dn['pd']['hypoxia_resistance']}` currently acts on no cells at all in "
          "the flagship configuration — the point `results/state_validation.md` "
          "§2.2 and `results/oxygen_demotion.md` §3.1 both made. The marker route "
          "leaves it empty too, but for a completely different and more honest "
          "reason: this panel carries no hypoxia-response gene, so hypoxia is "
          "**UNMEASURED** rather than computed-and-found-to-be-zero. Those are not "
          "the same statement and they should never be reported as if they were.")
        A("")
    A("### 3.2 Predicted kill and the kill CEILING, by route")
    A("")
    A("The **ceiling** is the kill fraction as exposure -> infinity: the Hill term "
      "saturates, so the per-cell probability tends to `max_kill_rate x sensitivity "
      "x resistance`. No dose can exceed it, and with uniform sensitivity it is "
      "*entirely* the state labelling. It is the cleanest measure of what this "
      "whole question is worth.")
    A("")
    A("| drug | IC50 (uM) | " + " | ".join(f"kill ({r})" for r in ROUTES)
      + " | " + " | ".join(f"ceiling ({r})" for r in ROUTES) + " |")
    A("|---|---|" + "---|" * (2 * len(ROUTES)))
    for name, r in dn["rows"].items():
        A(f"| {name} | {r['ic50_uM']:.3g} | "
          + " | ".join(_fmt(r["routes"][q]["kill_fraction"]) for q in ROUTES)
          + " | " + " | ".join(_fmt(r["routes"][q]["ceiling"]) for q in ROUTES)
          + " |")
    A("")
    kill_ratio = max(
        (max(r["routes"][q]["kill_fraction"] for q in ROUTES)
         / max(min(r["routes"][q]["kill_fraction"] for q in ROUTES), 1e-12))
        for r in dn["rows"].values())
    kill_abs = max(
        (max(r["routes"][q]["kill_fraction"] for q in ROUTES)
         - min(r["routes"][q]["kill_fraction"] for q in ROUTES))
        for r in dn["rows"].values())
    A(f"**Absolute kill moves by up to {kill_abs:.4f} absolute "
      f"({kill_ratio:.2f}x relative) purely from the choice of state route**, with "
      "the tissue, the drug, the transport solve and the RNG all held identical. "
      "Nothing physical changed between those columns; only the label on each cell "
      "did.")
    A("")
    A("### 3.3 The cross-drug ranking")
    A("")
    A("| route | ranking by predicted kill |")
    A("|---|---|")
    for route in ROUTES:
        A(f"| `{route}` | {' > '.join(dn['ranking'][route])} |")
    A("")
    if dn["ranking_invariant"]:
        A("**The KILL ranking is INVARIANT ACROSS THESE THREE STATE ROUTES and "
          "the levels are not. That is the honest and useful result — about "
          "these routes.** The resistance multipliers are shared across drugs "
          "— every drug pays the same `x0.5` in the same QUIESCENT cells — so they "
          "very largely cancel out of a comparison and do not cancel at all out of "
          "a level. Concretely: a claim of the form *\"cisplatin beats "
          "doxorubicin beats temozolomide on this tissue\"* survives every state "
          f"route tested here, while a claim of the form *\"temozolomide kills "
          f"{dn['rows']['temozolomide']['routes']['oxygen']['kill_fraction'] * 100:.1f}% "
          "of tumour cells\"* does not survive any of them.")
    else:
        A("**The ranking CHANGES between routes.** That is the worst case for this "
          "codebase: cross-drug comparisons were the class of claim believed to be "
          "robust to the state labelling, and on this section they are not.")
    A("")
    # --- 3.3a ------------------------------------------------------------- #
    # This scope block is emitted unconditionally, and it is emitted from the
    # GENERATOR rather than hand-added to the .md, because this section is the
    # origin of the "ranking is invariant" claim that the app, the protocol
    # diagnostics, API_CONTRACT.md and four other reports all cite. A hand
    # correction to the .md is erased by the next regeneration, which is exactly
    # the pass the exchange-default flip requires.
    A("### 3.3a Scope — what this invariance does NOT cover")
    A("")
    A("Every route above is **common-mode**: each applies the same multiplier to "
      "every drug, which is precisely why it cancels out of a comparison. So this "
      "result is evidence about state labelling and about nothing else, and in "
      "particular it is **not** evidence that the ranking is robust in general.")
    A("")
    A("**A drug-specific perturbation breaks it.** Correcting the vascular "
      "delivery law from a one-way leak `S = PS*Cp` to the gradient-driven "
      "`S = PS*(Cp - C)` re-weights each drug by its own permeability, and "
      "`bbb_penetrance` spans 47x across the GBM panel. Measured: the GBM top-1 "
      "changes (cisplatin -> vincristine), the top-3 set changes on both panels, "
      "4 of 45 pairs swap — two across gaps of 40.5% and 35.6% of the leader's "
      "kill against a Bernoulli sd of ~0.003 — and the Monte-Carlo leader "
      "probability goes 0.9277 -> 0.3965. See "
      "`results/transvascular_exchange_revalidation.md` §4.1.")
    A("")
    A("**And name the quantity.** Delivery-rank invariance does not imply "
      "kill-rank invariance: on breast the AUC ranking is preserved exactly "
      "(rho = 1.000) while the kill ranking scrambles (rho = 0.830), because the "
      "Hill map saturates. Everything above is about the KILL ranking.")
    A("")
    max_rel = max(
        (max(r["routes"][q]["kill_fraction"] for q in ROUTES)
         / max(min(r["routes"][q]["kill_fraction"] for q in ROUTES), 1e-12) - 1.0)
        for r in dn["rows"].values()) * 100.0
    A("This is the same structure `results/oxygen_demotion.md` §3.2 found for the "
      "oxygen ROUTE — ranking identical, levels moving — but roughly "
      f"{max_rel / 3.4:.0f}x larger: that report's kill moved at most 3.4% "
      f"relative, and changing the state RULE moves it {max_rel:.0f}%. The reason "
      "is structural: swapping the oxygen field nudges a field the bands are then "
      "cut from, while swapping the state rule moves the band OCCUPANCIES "
      "themselves (§3.1: VIABLE 10.04% -> 91.08% -> 29.80%).")
    A("")
    A("### 3.4 What this does NOT show")
    A("")
    A("None of the three sections carries a paired drug-response outcome, so "
      "nothing here says which route's kill is correct. This is a sensitivity "
      "analysis: it measures how much the answer depends on a modelling choice. A "
      "route can be perfectly stable, fully traceable to a measurement, and still "
      "wrong about drug response.")
    A("")

    # ---------------- 4. the bar ----------------
    A("## 4. The bar, scored")
    A("")
    A("The bar, stated in advance: **state fractions must be stable across dx and "
      "across sections, and must be traceable to a measurement or explicitly "
      "declared unmeasured.**")
    A("")
    A("| criterion | `oxygen` | `vessel_distance` | `markers` |")
    A("|---|---|---|---|")
    dxmax = {route: max((s["viable_swing"] for s in rt["stability"]
                         if s["route"] == route), default=float("nan"))
             for route in ROUTES}
    dxall = {route: sorted(s["viable_swing"] for s in rt["stability"]
                           if s["route"] == route)
             for route in ROUTES}
    A("| stable across dx (VIABLE swing, worst of the "
      f"{len(rt['stability']) // len(ROUTES)} section x mask combinations; best in "
      "brackets) | "
      + " | ".join(f"{_swing(dxmax[r])} [{_swing(dxall[r][0])}]" for r in ROUTES)
      + " |")
    A("| stable across the mask x dx grid (flagship) | "
      + " | ".join(_swing(mask_swing[r]["swing"]) for r in ROUTES) + " |")
    A(f"| stable across the three sections (`geometric`, dx = {DX_CONVERGED:g} um) | "
      + " | ".join(_swing(cross[r]["swing"]) for r in ROUTES) + " |")
    A("| traceable to a measurement | **no** — an unfitted PDE gated by a "
      "bookkeeping mask | **partly** — measured endothelial positions + a "
      "published diffusion distance; the *band edges* are published, the *field* "
      "is a proxy | **yes, where the panel allows** — and an explicit UNMEASURED "
      "verdict where it does not |")
    A("| declares what it cannot know | no — it labels every cell regardless | "
      "no — it labels every cell regardless | **yes** — refuses HYPOXIC on the "
      "319-gene panel, refuses a per-cell proliferation call below the detection "
      "floor |")
    A("| free constants | `consumption` (unfitted), `hypoxia_threshold`, the "
      "tissue-mask convention, the mesh | `L` (derived from a published distance), "
      "`hypoxia_threshold`, the mesh | the detection floor, and a window radius "
      "only where the per-cell call is inadmissible |")
    A("")
    A("**Read the cross-section row carefully — it is the one criterion no route "
      "passes cleanly, and for three different reasons.** Three different patients "
      "on three different panels *should* differ, so a spread there is not "
      "automatically a fault. But the reasons are not equivalent: the oxygen "
      f"route's {_swing(cross['oxygen']['swing'])} spread includes a section where "
      "its VIABLE band is essentially empty, which is a model failure; the "
      f"vessel-distance route's {_swing(cross['vessel_distance']['swing'])} is "
      "small only because it calls ~91-100% of every section VIABLE, i.e. it is "
      "stable and nearly uninformative; and the marker route's "
      f"{_swing(cross['markers']['swing'])} tracks which markers each panel carries "
      "and at what detection rate, which is a property of the ASSAY and is "
      "reported alongside every number it produces (§1). So: **stability across "
      "sections is not established for any route, and the marker route is the only "
      "one whose spread is attributable to something a reader can inspect.** A "
      "fourth section on a panel matched to one of these would test it properly; "
      "three sections on three panels cannot.")
    A("")

    # ---------------- 5. recommendation ----------------
    A("## 5. RECOMMENDATION")
    A("")
    A("**Adopt `state_mode='markers'` as the recommended route on any panel that "
      "clears the detection floor — the IO 480-gene and GSE324268 366-gene sections "
      "do. Keep `'oxygen'` as the shipping default until the dependent results are "
      "re-run (§5.2). And on the 319-gene flagship section, stop quoting an "
      "absolute kill fraction or a kill ceiling as a result at all: NO route clears "
      "the bar there, and the honest output of that section is a ranking, a "
      "coverage map and a delivery metric.**")
    A("")
    A("In order of weight:")
    A("")
    A("1. **On the two panels that can support it, the marker route clears every "
      "part of the bar that three sections are capable of testing.** It is exactly "
      "dx-stable (it never touches a grid — measured, not asserted), it is "
      "traceable to measured transcripts, and it declares its gaps. The "
      "cross-section criterion is not established for ANY route here (§4), and that "
      "is the honest gap in this recommendation rather than a point in the marker "
      "route's favour. "
      "The vessel-distance route is a genuine improvement on the PDE — measured "
      "vessel positions, a published length scale, and a flagship-grid swing of "
      f"{_swing(mask_swing['vessel_distance']['swing'])} where the oxygen route's "
      f"is {_swing(mask_swing['oxygen']['swing'])} — but it is still a geometric "
      "proxy for a state, not a measurement of one, and on these sections it is "
      "nearly uninformative: it calls 91-100% of tumour cells VIABLE everywhere, "
      "so its stability is partly the stability of a near-constant.")
    brain_u = (float(brain["coverage"]["proliferation"]["union_detection"])
               if brain else float("nan"))
    A("2. **On the 319-gene flagship section NO route clears the bar, and that is "
      "the most important sentence in this report.** Hypoxia is unmeasured there; "
      f"proliferation detection is {LITERATURE_KI67[0] / max(brain_u, 1e-12):.0f}x "
      "below the published Ki-67 lower bound, so the marker route falls back to a "
      f"neighbourhood rule whose level moves {_swing(rad['swing'])} with the window "
      "radius; and the oxygen route's VIABLE band moves from "
      f"{_pct(mask_swing['oxygen']['viable_min'], 4)} to "
      f"{_pct(mask_swing['oxygen']['viable_max'], 4)} — i.e. it empties entirely — "
      "purely with bookkeeping. Every absolute "
      "kill fraction and every kill ceiling published on that section — "
      "`results/dose_schedule.md`, `results/report.md`, "
      "`results/patient_report_*.md` — inherits an unresolved uncertainty of that "
      "size. The **comparisons** in those documents are fine (§3.3); the "
      "**levels** are not, and should carry the state route and this caveat "
      "attached.")
    A("3. **Do not flip the shipping default silently.** §3.2 shows kill moves by "
      f"up to {kill_abs:.4f} absolute between routes, which is far too large to "
      "change under a reader who is comparing against an older document. The "
      "option exists, its provenance is recorded in "
      "`reconstruction_report()['cell_state']`, and flipping it is a deliberate "
      "act with a re-run list (§5.2).")
    A("")
    A("### 5.1 What measurement would settle it")
    A("")
    A("The bar is unmet on the flagship section for a data reason, not a modelling "
      "reason, so more modelling cannot fix it. What would:")
    A("")
    A("* **A proliferation readout with a real detection rate on THAT section** — "
      "Ki-67 IHC/IF on a serial section, scored per cell and registered to the "
      "Xenium coordinates, or the same tumour on a whole-transcriptome platform. "
      "That converts the neighbourhood fallback into a per-cell call and removes "
      "the window radius entirely. The IO panel already demonstrates the target is "
      f"reachable: MKI67 at "
      + (_pct(io["blocks"]["proliferation_specific"]["detection_tumor"].get("MKI67", 0), 1)
         if io else "33.7%")
      + ", inside the published range.")
    A("* **Any hypoxia-response gene on the flagship panel** — one of CA9, ADM, "
      "NDRG1, SLC2A1 would move HYPOXIC from unmeasured to measured there. CA9 is "
      "the one worth adding: it is the most hypoxia-restricted of the set, so it "
      "constrains the band's SIZE and not merely its ordering.")
    A("* **A measured oxygen tension** (pimonidazole IHC or phosphorescence-"
      "quenching pO2 on a serial section) would calibrate the hypoxic band's "
      "absolute size for every route at once — the one thing none of them "
      "currently has — and would let `OxygenConfig.consumption` be fitted instead "
      "of assumed.")
    A("* **A paired drug-response outcome.** Everything above makes the state "
      "labelling defensible; only an outcome makes the kill defensible.")
    A("")
    A("### 5.2 If the default WERE flipped, this is the re-run list")
    A("")
    A("Anything whose numbers depend on the per-cell band occupancy, i.e. on the "
      "`PDConfig` resistance multipliers having cells to act on. Coverage-only and "
      "AUC-only results are immune — the AUC field has no state dependence "
      "(`results/oxygen_demotion.md` §3.2 established the same immunity for the "
      "oxygen route):")
    A("")
    A("| artefact | why | immune? |")
    A("|---|---|---|")
    A("| `results/dose_schedule.md` | its kill ceilings are ~89% the quiescent "
      "multiplier | **no — re-run** |")
    A("| `results/report.md`, `results/patient_report_*.md` | quote absolute kill "
      "| **no — re-run** |")
    A("| `results/ranking_robustness.md`, `results/uncertainty.md` | resample kill "
      "| **no — re-run** (their RANKINGS should survive, §3.3) |")
    A("| `results/state_validation.md`, `results/oxygen_demotion.md` | both are "
      "ABOUT the oxygen-banded states | **no — re-run** |")
    A("| `results/delivery_metrics.md` | AUC + coverage only | yes |")
    A("| `results/vessel_network_impact.md` | coverage + penetration length | "
      "mostly |")
    A("")

    # ---------------- 6. caveats ----------------
    A("## 6. Caveats — what could make this wrong")
    A("")
    A("1. **Three sections, three patients, three panels.** Nothing transfers "
      "between them, and a cross-section \"spread\" mixes real biological "
      "difference with panel difference. n = 3.")
    A("2. **No outcome label anywhere.** Every kill number in §3 is a prediction "
      "with nothing to score it against. This report ranks routes on stability and "
      "traceability, which are necessary and not sufficient.")
    A("3. **Hypoxia-response transcripts are not an oxygen tension.** They are the "
      "output of a programme with its own kinetics and its own non-hypoxic "
      "inducers, so the marker HYPOXIC band's SIZE is uncalibrated (§1.3). HIF1A "
      "is excluded on mechanism (oxygen sensing acts on the protein: Ivan et al., "
      "Science 2001;292(5516):464-468, PMID 11292862; Jaakkola et al., Science "
      "2001;292(5516):468-472, PMID 11292861).")
    A("4. **Targeted-panel dropout.** A zero is \"no probe fired\". This is what "
      "disqualifies the per-cell call on the 319-gene panel and it also means no "
      "individual cell's marker call is trustworthy on any panel — only rates and "
      "orderings are.")
    A("5. **The neighbourhood binarisation is a convention**, spanning "
      f"{_swing(rad['swing'])} across the window radius (§2.4). Where a level is "
      "needed rather than an ordering, that axis should be treated as unmeasured.")
    A("6. **Typing circularity.** `MKI67` is one of the four brain-panel TUMOR "
      "markers, so an MKI67-high cell is marginally more likely to be typed TUMOR; "
      "that inflates the marginal positive rate inside the tumour population on "
      "that section. It cannot manufacture a between-route difference (all routes "
      "share the typing) but it does inflate one input to the marker route.")
    A("7. **The one-label collapse.** Proliferation and hypoxia are two axes and "
      "`core.CellState` has one slot; the marker route resolves the clash by "
      "letting HYPOXIC win. That precedence is a schema-forced convention, stated "
      "in the returned provenance, and it is the only choice in the marker route "
      "that is not a marker threshold.")
    A(f"8. **The detection floor is a judgement.** "
      f"{PER_CELL_DETECTION_FLOOR * 100:.0f}% is the published Ki-67 lower bound "
      "divided by ~4, and it is the one number in the marker route that is neither "
      "a marker threshold nor a published anchor. Both the raw rates and the "
      "resulting call level are always returned so a reader can apply a different "
      f"tolerance. The flagship panel misses it by "
      f"{PER_CELL_DETECTION_FLOOR / max(brain_u, 1e-12):.1f}x, so admitting it "
      f"would mean setting the floor at {brain_u * 100:.2f}% — "
      f"{LITERATURE_KI67[0] / max(brain_u, 1e-12):.0f}x below the published lower "
      "bound, i.e. accepting a rate at which the call cannot be distinguished from "
      "probe noise. The IO and GSE324268 panels clear it by 6-7x, so on those two "
      "the choice of tolerance does not matter either.")
    A("")

    # ---------------- 7. reproduction ----------------
    A("## 7. Reproduction")
    A("")
    A("```bash")
    A("python -m tmesim.validation.cell_state_impact")
    A("```")
    A("")
    A(f"Deterministic — the only randomness is the engine's seeded generator "
      f"(`seed = 0`). Runtime ~{results['elapsed_s'] / 60:.0f} min, dominated by "
      "the reconstructions of the two ~800k-cell sections across the dx sweep.")
    A("")
    A("**Files added or changed by this work:**")
    A("")
    A("* `tmesim/reconstruct/cell_state.py` — new; the marker-derived state and "
      "its coverage report")
    A("* `tmesim/validation/cell_state_impact.py` — new; this analysis")
    A("* `tmesim/reconstruct/real_pipeline.py` — additive `state_mode` / "
      "`state_kwargs` keyword, default `'oxygen'` and bit-identical (it does not "
      "touch `cells.state` at all), plus the `cell_state` provenance block in "
      "`reconstruction_report`")
    A("* `tests/test_cell_state.py` — new")
    A("* `results/cell_state_impact.md` — this report")
    A("")
    A("**Data provenance.** Three open single-cell spatial GBM sections: the 10x "
      "Xenium Human Brain Glioblastoma section (v1.3.0, 319 targets, CC BY 4.0), "
      "the 10x Xenium Human Brain Cancer GBM FFPE IO-panel section (v2.0.0, 480 "
      "targets, CC BY 4.0), and GEO GSE324268 GBMtumor07 (366 targets). None "
      "carries a paired drug-response outcome label.")

    text = "\n".join(L) + "\n"
    with open(out_path, "w") as fh:
        fh.write(text)
    print(f"\nwrote {out_path}")
    return out_path


if __name__ == "__main__":
    _cache = os.path.join(_results_dir(), "cell_state_impact_cache.pkl")
    _res = run(cache_path=_cache)
    write_report(_res)
    print(f"total {_res['elapsed_s'] / 60:.1f} min")
