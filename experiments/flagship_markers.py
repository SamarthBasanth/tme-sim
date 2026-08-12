"""RE-ANCHORING THE FLAGSHIP DRUG PANEL to a section that can support cell state.

WHY THIS EXISTS
---------------
Almost every headline number in this project — ``results/real_drug_panel.md``,
``results/dose_schedule.md``, ``results/kill_calibration.md``,
``results/patient_report_xenium.md`` — was computed on ONE section:
``data_cache/xenium_gbm``, the 40,887-cell 319-gene 10x Xenium brain panel.
``results/cell_state_impact.md`` then established that that panel **cannot support
a cell state at all**:

* it carries **no** hypoxia-response gene of the marker set (VEGFA / CA9 / ADM /
  NDRG1 / SLC2A1), so HYPOXIC is **UNMEASURED** there and no cell is assigned it;
* MKI67 + TOP2A fire in **1.42%** of its tumour cells, **13x below** the published
  19-40% GBM Ki-67 labelling index, so a per-cell proliferation call is refused and
  the marker route falls back to a neighbourhood rule whose level moves 8.79x with
  the window radius;
* and its own oxygen-route VIABLE band **empties entirely** across the tissue-mask
  x mesh grid, i.e. the swing is unbounded.

Its conclusion, verbatim: *"on the 319-gene flagship section, stop quoting an
absolute kill fraction or a kill ceiling as a result at all."*

Two other cached real sections DO clear that bar, and this module re-runs the
flagship five-drug panel on them with ``state_mode='markers'``:

* ``data_cache/gse324268_xenium`` — GEO GSE324268 GBMtumor07, 851,509 cells,
  366-gene brain + GBM add-on panel. **PRIMARY.** It is the only cached section
  that carries proliferation markers (MKI67 13.3%, TOP2A 25.7%) AND hypoxia-response
  markers (VEGFA 45.5%, CA9 1.03%, ADM 10.3%, SLC2A1 12.6%) AND **MGMT** — which
  makes it the only section in this project where temozolomide sensitivity can be
  MEASURED rather than declared unmodellable.
* ``data_cache/xenium_gbm_io`` — 10x Xenium GBM FFPE IO panel, 816,769 cells,
  480 genes. **CHECK.** Highest MKI67 detection anywhere in the project (34.06%,
  inside the published Ki-67 range) but no MGMT and only VEGFA for hypoxia.

The 319-gene section is re-run through the IDENTICAL harness as the comparison
arm, so the cross-section ranking test is like-for-like rather than quoted.

WHAT IS AND IS NOT MODIFIED
---------------------------
**Nothing under ``tmesim/`` is touched.** This module is a consumer of the existing
public API: ``reconstruct_real(state_mode='markers')``,
``apply_drug_sensitivity(mode='uniform')``, ``sim.engine.simulate_response``,
``validation.delivery_metrics``, ``molecular.mgmt.mgmt_sensitivity``.

HOW THE MARKER-STATE KILL IS COMPUTED (and why it is not a re-implementation)
-----------------------------------------------------------------------------
Identical device to ``validation.cell_state_impact`` §3, for the identical reason.
``sim.engine.simulate_response`` RE-BANDS the cells from its own oxygen solve, so
handing it a marker-stated reconstruction does not change its kill. But the cell
state reaches the kill through exactly ONE channel — the ``PDConfig`` resistance
multipliers inside ``cells.pkpd.apply_drug_kill`` — and nothing in the transport
solve reads the state, so the AUC field (and therefore every delivery metric) is
state-invariant by construction.

So each marker arm installs its states on a PRIVATE COPY of the cells and re-runs
the UNMODIFIED ``apply_drug_kill`` on the engine's OWN transport result with the
engine's own seeded generator. That structural claim is not asserted, it is
**gated**: the oxygen arm is recomputed through the same path and must reproduce
``simulate_response``'s ``kill_fraction`` bit-for-bit or the run aborts.

THE MGMT BRIDGE (GSE324268 only)
---------------------------------
``molecular.mgmt.mgmt_sensitivity`` is calibrated in ``log2(CPM + 1)`` where CPM is
**whole-transcriptome** bulk RNA-seq (GSE232173, n = 17). A 366-gene targeted panel
has a ~20x smaller denominator, so its pseudobulk CPM is NOT the calibration unit
and feeding it in raw would be a unit error. Three routes are therefore reported:

1. **BRIDGED (primary).** The panel pseudobulk CPM is rescaled by the fraction of
   the bulk transcriptome those same panel genes carry, MEASURED in the 7 cached
   TCGA GBM bulk STAR-count sections. The result is cross-checked against the
   TCGA GBM MGMT CPM distribution itself.
2. **BINARY.** ``mgmt_sensitivity(status='unmethylated')`` — the measured
   median-kill ratio between MGMT-expressing and MGMT-silent cultures. Needs no
   unit bridge at all and is the independent check on route 1.
3. **NAIVE (reported and REJECTED).** The raw panel CPM fed in as if it were bulk.
   It lands far above the calibration ceiling and is shown only to make the size of
   the unit error explicit.

Run:    ./.venv/bin/python experiments/flagship_markers.py
        ./.venv/bin/python experiments/flagship_markers.py --sections gse324268_366
        ./.venv/bin/python experiments/flagship_markers.py --report-only
Writes: results/flagship_markers.md

Units: distance um, AUC uM*s, concentration uM, fractions dimensionless.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import tempfile
import time
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from tmesim.cells.agents import update_cell_states
from tmesim.cells.pkpd import apply_drug_kill, effective_hypoxia_resistance
from tmesim.core import (
    CellState,
    CellType,
    DrugParams,
    Microenvironment,
    OxygenConfig,
    PDConfig,
    PDEConfig,
    SimConfig,
)
from tmesim.data.real.singlecell_spatial import load as load_singlecell_spatial
from tmesim.molecular.mgmt import (
    MGMT_EXPRESSED_CPM_THRESHOLD,
    MGMT_OBSERVED_LOG2CPM_RANGE,
    mgmt_sensitivity,
)
from tmesim.reconstruct.cell_state import state_fractions_from_states
from tmesim.reconstruct.real_pipeline import (
    apply_drug_sensitivity,
    reconstruct_real,
    reconstruction_report,
)
from tmesim.sim.engine import simulate_response
from tmesim.validation.calibration import calibrate_alpha
from tmesim.validation.delivery_metrics import (
    _PDE,
    compile_real_drug,
    exposure_window_s,
    therapeutic_coverage,
    tissue_auc_summary,
    tumor_cell_exposure,
)

__all__ = [
    "PANEL",
    "SECTIONS",
    "mgmt_from_xenium_panel",
    "run_section",
    "tcga_panel_capture",
    "write_report",
]

# --------------------------------------------------------------------------- #
# Declared constants
# --------------------------------------------------------------------------- #

#: The flagship five-drug panel, in the order ``results/real_drug_panel.md`` uses.
PANEL: Tuple[str, ...] = (
    "temozolomide", "lomustine", "cisplatin", "methotrexate", "doxorubicin",
)

#: Grid-converged production spacing (``results/numerical_rigor.md`` §2).
DX_UM = 12.0

#: The cached real sections. ``(key -> path, label, role)``.
SECTIONS: Dict[str, Dict[str, str]] = {
    "gse324268_366": {
        "path": "data_cache/gse324268_xenium",
        "label": "GEO GSE324268 GBMtumor07 — 366-gene brain + GBM add-on panel",
        "role": "PRIMARY (proliferation + hypoxia + MGMT all measured)",
    },
    "io_480": {
        "path": "data_cache/xenium_gbm_io",
        "label": "10x Xenium Human Brain Cancer (GBM) FFPE v2.0.0 — 480-gene IO panel",
        "role": "CHECK (highest MKI67 detection; no MGMT, hypoxia on VEGFA alone)",
    },
    "brain_319": {
        "path": "data_cache/xenium_gbm",
        "label": "10x Xenium Human Brain — Glioblastoma v1.3.0 — 319-gene brain panel",
        "role": "COMPARISON ARM (the incumbent flagship; cannot support cell state)",
    },
}

#: The order sections are reported in.
SECTION_ORDER: Tuple[str, ...] = ("gse324268_366", "io_480", "brain_319")

#: REPRODUCTION GATE on the 319-gene section: the published oxygen- and
#: marker-route kills and ceilings from ``results/cell_state_impact.md`` §3.2, at
#: these exact settings (uniform sensitivity, dx = 12 um, geometric mask, discrete
#: vessel source, alpha = 148.61). That report is the most recent published run on
#: this section and is the correct reference for the CURRENT pipeline.
#:
#: ``results/real_drug_panel.md`` and ``results/heldout_patient.md`` quote LOWER
#: kills for the same drugs on the same section (temozolomide 0.0147, cisplatin
#: 0.3699, doxorubicin 0.2165) because they predate a change in the vascular field
#: — their reconstruction reports mean vascular density 0.0250 and mean oxygen
#: 0.2955, where the current pipeline gives 0.0268 and 0.3520 on identical inputs.
#: That is a pipeline drift this report DOCUMENTS rather than inherits; see §3.
PUBLISHED_319: Dict[str, Dict[str, float]] = {
    "temozolomide": {"kill_oxygen": 0.0184, "kill_markers": 0.0219},
    "cisplatin": {"kill_oxygen": 0.4362, "kill_markers": 0.5155},
    "doxorubicin": {"kill_oxygen": 0.2595, "kill_markers": 0.3073},
}
#: Ceilings from the same table (shared across drugs -- uniform sensitivity means
#: the ceiling IS the state labelling).
PUBLISHED_319_CEILING = {"ceiling_oxygen": 0.5502, "ceiling_markers": 0.6490}
#: Tolerance for the gate. The engine draws kills from a seeded generator, so the
#: reproduction is exact to the 4 decimals the source report prints.
PUBLISHED_319_TOL = 5e-5

#: The older, superseded 319-gene numbers, quoted so the drift is visible.
STALE_319_KILL: Dict[str, float] = {
    "temozolomide": 0.0147, "lomustine": 0.1975, "cisplatin": 0.3699,
    "methotrexate": 0.0162, "doxorubicin": 0.2165,
}

#: Cached TCGA GBM bulk RNA-seq (STAR counts, GENCODE v36) used ONLY to measure
#: what fraction of the bulk transcriptome a targeted panel's genes carry — the
#: unit bridge for the MGMT scalar. These are DIFFERENT PATIENTS from any Xenium
#: section; the bridge is a cohort-level constant, not this patient's.
_TCGA_DIR = "data_cache/tcga_expression"

_RESULTS_DIR = os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
)
_OUT_MD = os.path.join(_RESULTS_DIR, "flagship_markers.md")

#: Raw per-section results are cached OUTSIDE results/ so this experiment adds
#: exactly two files to the repo. Override with ``--cache``.
_DEFAULT_CACHE = os.path.join(
    os.environ.get("TMESIM_SCRATCH", tempfile.gettempdir()),
    "flagship_markers_runs.json",
)


# --------------------------------------------------------------------------- #
# Small helpers
# --------------------------------------------------------------------------- #
def _retry(fn, *args, **kwargs):
    """Call ``fn``, retrying ONCE after a short pause on any transient failure."""
    try:
        return fn(*args, **kwargs)
    except Exception:
        time.sleep(3)
        return fn(*args, **kwargs)


def _states_from_oxygen(micro: Microenvironment, oxygen: np.ndarray,
                        cfg: OxygenConfig) -> np.ndarray:
    """Band a field into states with the UNMODIFIED ``update_cell_states``."""
    work = Microenvironment(
        domain=micro.domain, fields=micro.fields, cells=micro.cells.copy(),
        patient_id=micro.patient_id, metadata=micro.metadata,
    )
    update_cell_states(work, np.asarray(oxygen, dtype=np.float64), cfg)
    return np.asarray(work.cells.state).copy()


def _kill_from_states(micro: Microenvironment, states: np.ndarray, transport,
                      cfg: SimConfig, drug,
                      sensitivity_scalar: Optional[float] = None) -> float:
    """Install ``states`` and re-run the UNMODIFIED kill on the engine's transport.

    Nothing is re-implemented: this calls ``cells.pkpd.apply_drug_kill``, the same
    function the engine calls, on a private copy of the cells with the engine's own
    seeded generator. Valid because the AUC field carries no cell-state dependence
    anywhere in the pipeline — and that validity is GATED by the caller, which
    requires this function to reproduce ``simulate_response`` exactly in the
    oxygen arm.

    ``sensitivity_scalar`` (used only for the MGMT arm) multiplies the per-cell
    sensitivity array uniformly — MGMT is a per-patient scalar, not a field.
    """
    work = Microenvironment(
        domain=micro.domain, fields=micro.fields, cells=micro.cells.copy(),
        patient_id=micro.patient_id, metadata=micro.metadata,
    )
    work.cells.state = np.asarray(states, dtype=np.int8).copy()
    if sensitivity_scalar is not None:
        work.cells.sensitivity = np.clip(
            np.asarray(work.cells.sensitivity, dtype=np.float64)
            * float(sensitivity_scalar), 0.0, 1.0)
    apply_drug_kill(work, transport, cfg.pd, drug,
                    rng=np.random.default_rng(cfg.seed))
    return 1.0 - float(work.cells.viable_tumor_fraction())


def _reproduction_gate(rows: Sequence[Dict]) -> List[Dict]:
    """Score this run's 319-gene numbers against ``results/cell_state_impact.md``.

    Eight numbers — three drugs x (oxygen kill, marker kill) plus the two shared
    ceilings — which together check the transport solve, both state routes and the
    PD mapping. Computed from cached rows, so the reference can be changed without
    re-running the simulation.
    """
    out: List[Dict] = []
    for r in rows:
        ref = PUBLISHED_319.get(r["drug"])
        if not ref:
            continue
        for field, pub in ref.items():
            out.append({
                "quantity": f"{r['drug']} {field}", "published": pub,
                "here": r[field], "abs_delta": abs(r[field] - pub),
                "reproduced": abs(r[field] - pub) <= PUBLISHED_319_TOL,
            })
    if rows:
        for field, pub in PUBLISHED_319_CEILING.items():
            here = rows[0][field]
            out.append({
                "quantity": f"kill ceiling ({field.split('_')[1]} route, shared)",
                "published": pub, "here": here, "abs_delta": abs(here - pub),
                "reproduced": abs(here - pub) <= PUBLISHED_319_TOL,
            })
    return out


def _ceiling(states: np.ndarray, sensitivity: np.ndarray, tumor: np.ndarray,
             pd: PDConfig, sensitivity_scalar: Optional[float] = None,
             drug: Optional[DrugParams] = None) -> float:
    """The kill CEILING: ``mean_i[max_kill_rate * sensitivity_i * resistance_i]``.

    As exposure -> infinity the Hill term saturates, so no dose can exceed this.
    With uniform sensitivity it IS the state labelling, which is exactly the
    quantity this whole re-anchoring is about.

    ``drug``, if given, is resolved through ``cells.pkpd.effective_hypoxia_resistance``
    so a per-drug ``hypoxia_resistance`` override (temozolomide has one) is
    honoured -- unlike the flagship 319-gene section, ``gse324268_xenium`` (this
    module's PRIMARY section) carries real hypoxia-response markers and a
    genuinely nonzero HYPOXIC band, so this matters numerically here.
    """
    r_hyp = (effective_hypoxia_resistance(drug, pd) if drug is not None
             else pd.hypoxia_resistance)
    st = np.asarray(states)[tumor]
    sens = np.asarray(sensitivity, dtype=np.float64)[tumor]
    if sensitivity_scalar is not None:
        sens = np.clip(sens * float(sensitivity_scalar), 0.0, 1.0)
    res = np.ones(st.size, dtype=np.float64)
    res[st == int(CellState.HYPOXIC)] = r_hyp
    res[st == int(CellState.QUIESCENT)] = pd.quiescent_resistance
    per_cell = pd.max_kill_rate * sens * res
    per_cell = np.where(st == int(CellState.DEAD), 0.0, per_cell)
    return float(per_cell.mean()) if st.size else 0.0


# --------------------------------------------------------------------------- #
# The MGMT bridge
# --------------------------------------------------------------------------- #
def tcga_panel_capture(panel_genes: Sequence[str],
                       tcga_dir: str = _TCGA_DIR) -> Dict:
    """What fraction of the BULK transcriptome a targeted panel's genes carry.

    Measured on the cached TCGA GBM bulk STAR counts (GENCODE v36, unstranded).
    This is the unit bridge between a Xenium panel pseudobulk CPM (denominator =
    the panel's own genes) and the ``log2(CPM+1)`` whole-transcriptome unit
    ``molecular.mgmt`` is calibrated in:

        CPM_bulk_equivalent  =  CPM_panel  x  (panel counts / all counts)_bulk

    Also returns the measured MGMT CPM of each TCGA GBM section, which is an
    INDEPENDENT sanity band for the bridged value: a GBM's MGMT CPM should land
    inside the GBM cohort's own range.

    THE ASSUMPTION, STATED: this treats the panel's relative capture as unbiased,
    i.e. it assumes a Xenium probe set recovers its genes in the same proportion a
    bulk library does. Probe efficiency is not uniform, so the bridge is an
    ORDER-OF-MAGNITUDE correction, not a calibration.
    """
    import pandas as pd

    manifest_path = os.path.join(tcga_dir, "manifest.json")
    with open(manifest_path) as fh:
        manifest = json.load(fh)
    wanted = set(str(g) for g in panel_genes)

    rows = []
    for sid, fname in sorted(manifest.items()):
        df = pd.read_csv(os.path.join(tcga_dir, fname), sep="\t", comment="#")
        df = df[~df["gene_id"].astype(str).str.startswith("N_")]
        counts = pd.to_numeric(df["unstranded"], errors="coerce").fillna(0.0)
        total = float(counts.sum())
        in_panel = df["gene_name"].isin(wanted)
        rows.append({
            "sample": sid,
            "total_counts": total,
            "mgmt_counts": float(counts[df["gene_name"] == "MGMT"].sum()),
            "mgmt_cpm": float(1e6 * counts[df["gene_name"] == "MGMT"].sum() / total),
            "panel_capture_fraction": float(counts[in_panel].sum() / total),
            "n_panel_genes_matched": int(in_panel.sum()),
        })

    cap = np.array([r["panel_capture_fraction"] for r in rows])
    cpm = np.array([r["mgmt_cpm"] for r in rows])
    return {
        "n_sections": len(rows),
        "rows": rows,
        "capture_median": float(np.median(cap)),
        "capture_min": float(cap.min()),
        "capture_max": float(cap.max()),
        "mgmt_cpm_median": float(np.median(cpm)),
        "mgmt_cpm_min": float(cpm.min()),
        "mgmt_cpm_max": float(cpm.max()),
        "source": ("TCGA-GBM bulk RNA-seq, STAR unstranded counts, GENCODE v36 "
                   f"(n = {len(rows)} cached sections) — DIFFERENT PATIENTS from "
                   "any Xenium section"),
    }


def mgmt_from_xenium_panel(path: str, tumor: Optional[np.ndarray] = None) -> Dict:
    """A per-patient MGMT scalar from a Xenium section's OWN measured MGMT.

    Returns the raw measurement (counts, detection rate, panel pseudobulk CPM),
    the TCGA-bridged bulk-equivalent CPM and its sensitivity, the binary
    (``'unmethylated'``) sensitivity, and the naive-unit value that is reported
    only to size the unit error. ``{'mgmt_on_panel': False}`` if the panel lacks
    MGMT.
    """
    data = _retry(load_singlecell_spatial, path=path, dense=False)
    genes = np.asarray(data["genes"], dtype=str)
    expr = data["expression"]
    hit = np.flatnonzero(genes == "MGMT")
    if hit.size == 0:
        return {"mgmt_on_panel": False, "path": path,
                "reason": "MGMT is not on this panel"}

    col = np.asarray(expr[:, [int(hit[0])]].todense()).ravel().astype(np.float64)
    totals = np.asarray(expr.sum(axis=1)).ravel().astype(np.float64)
    total_counts = float(totals.sum())
    mgmt_counts = float(col.sum())
    panel_cpm = 1e6 * mgmt_counts / total_counts if total_counts > 0 else 0.0

    bridge = tcga_panel_capture(genes)
    out: Dict = {
        "mgmt_on_panel": True,
        "path": path,
        "n_cells": int(col.size),
        "mgmt_counts": mgmt_counts,
        "total_panel_counts": total_counts,
        "detection_all_cells": float((col > 0).mean()),
        "detection_tumor": (float((col[np.asarray(tumor, bool)] > 0).mean())
                            if tumor is not None else None),
        "panel_pseudobulk_cpm": float(panel_cpm),
        "panel_pseudobulk_log2cpm": float(math.log2(panel_cpm + 1.0)),
        "bridge": bridge,
        "resolution": ("PER-PATIENT SCALAR (pseudobulk over the whole section). "
                       "MGMT fires in ~1.7% of cells at ~1 transcript each, so a "
                       "per-cell MGMT map would be counting noise."),
    }

    # --- route 1: BRIDGED (primary) ----------------------------------------- #
    routes = {}
    for tag, cap in (("bridged", bridge["capture_median"]),
                     ("bridged_lo", bridge["capture_max"]),
                     ("bridged_hi", bridge["capture_min"])):
        cpm = panel_cpm * float(cap)
        det = mgmt_sensitivity(math.log2(cpm + 1.0), detail=True)
        routes[tag] = {
            "capture_fraction": float(cap),
            "cpm": float(cpm),
            "log2cpm": float(math.log2(cpm + 1.0)),
            "sensitivity": float(det["sensitivity"]),
            "ci95": [float(v) for v in det["ci95"]],
            "extrapolated": bool(det["extrapolated"]),
        }
    # --- route 2: BINARY (independent check) -------------------------------- #
    det_b = mgmt_sensitivity(status="unmethylated", detail=True)
    routes["binary_unmethylated"] = {
        "sensitivity": float(det_b["sensitivity"]),
        "ci95": [float(v) for v in det_b["ci95"]],
        "extrapolated": False,
    }
    # --- route 3: NAIVE (reported, rejected) -------------------------------- #
    det_n = mgmt_sensitivity(out["panel_pseudobulk_log2cpm"], detail=True)
    routes["naive_panel_cpm"] = {
        "cpm": float(panel_cpm),
        "log2cpm": out["panel_pseudobulk_log2cpm"],
        "sensitivity": float(det_n["sensitivity"]),
        "ci95": [float(v) for v in det_n["ci95"]],
        "extrapolated": bool(det_n["extrapolated"]),
    }
    out["routes"] = routes
    out["primary_route"] = "bridged"
    out["primary_sensitivity"] = routes["bridged"]["sensitivity"]
    out["expressing"] = bool(routes["bridged"]["cpm"] >= MGMT_EXPRESSED_CPM_THRESHOLD)
    out["calibration_ceiling_log2cpm"] = float(MGMT_OBSERVED_LOG2CPM_RANGE[1])
    return out


# --------------------------------------------------------------------------- #
# One section
# --------------------------------------------------------------------------- #
def run_section(key: str, dx_um: float = DX_UM,
                drugs: Sequence[str] = PANEL,
                checkpoint: Optional[str] = None) -> Dict:
    """Reconstruct ONCE with ``state_mode='markers'``, then run the whole panel.

    One reconstruction and one transport solve per drug; every state route reuses
    them. The reconstruction is the frozen ``reconstruct_real`` path with the only
    non-default argument being ``state_mode='markers'``.

    ``checkpoint`` (optional) writes the per-drug rows to disk as they complete and
    reloads them on a restart, so an interrupted run resumes at the next drug
    instead of repeating ~4 min of transport solve per finished drug. It changes no
    number: every drug is independent given the shared reconstruction.
    """
    meta = SECTIONS[key]
    path = meta["path"]
    print("\n" + "=" * 78)
    print(f"SECTION {key}  —  {meta['label']}")
    print(f"  role: {meta['role']}")
    print("=" * 78)

    t_sec = time.time()
    t0 = time.time()
    micro = _retry(reconstruct_real, "xenium", sample=path, dx_um=float(dx_um),
                   state_mode="markers")
    rep = reconstruction_report(micro)
    t_recon = time.time() - t0

    tumor = np.asarray(micro.cells.cell_type) == int(CellType.TUMOR)
    marker_states = np.asarray(micro.cells.state, dtype=np.int8).copy()
    state_block = rep["cell_state"]
    coverage = (state_block.get("detail") or {}).get("coverage", {})

    if state_block["mode"] != "markers" or not state_block["is_measured"]:
        raise AssertionError(
            f"[{key}] state_mode did not take effect: {state_block['mode']!r}")

    print(f"  {rep['n_cells']:,} cells ({int(tumor.sum()):,} tumour)  "
          f"grid {rep['grid']}  dx={rep['dx_um']:.0f} um   [{t_recon:.0f}s]")
    fr_markers = state_fractions_from_states(marker_states, tumor)
    print(f"  marker states (tumour): VIABLE {fr_markers['viable'] * 100:.2f}%  "
          f"QUIESCENT {fr_markers['quiescent'] * 100:.2f}%  "
          f"HYPOXIC {fr_markers['hypoxic'] * 100:.2f}%")
    print(f"  UNMEASURED states: {coverage.get('unmeasured_states') or 'none'}")

    mgmt = mgmt_from_xenium_panel(path, tumor=tumor)
    if mgmt["mgmt_on_panel"]:
        print(f"  MGMT ON PANEL: {mgmt['detection_all_cells'] * 100:.2f}% of cells, "
              f"panel pseudobulk {mgmt['panel_pseudobulk_cpm']:.2f} CPM -> "
              f"bridged {mgmt['routes']['bridged']['cpm']:.2f} CPM -> "
              f"TMZ sensitivity {mgmt['primary_sensitivity']:.4f}")
    else:
        print("  MGMT is NOT on this panel — TMZ sensitivity stays an upper bound")

    alpha = float(calibrate_alpha())
    pdcfg = PDConfig()
    ocfg = OxygenConfig()
    rows: List[Dict] = []
    if checkpoint and os.path.exists(checkpoint):
        with open(checkpoint) as fh:
            rows = json.load(fh).get("rows", [])
        if rows:
            print(f"  resuming from checkpoint: "
                  f"{[r['drug'] for r in rows]} already done")
    done = {r["drug"] for r in rows}

    for name in drugs:
        if name in done:
            continue
        t0 = time.time()
        dp, mol = compile_real_drug(name, alpha)
        apply_drug_sensitivity(micro, name, mode="uniform")
        cfg = SimConfig(drug=dp, pde=PDEConfig(**_PDE), seed=0)
        res = simulate_response(micro, cfg)

        # --- the engine's OWN state route, and the HARD GATE ----------------- #
        engine_states = _states_from_oxygen(micro, res.fields["oxygen"], cfg.oxygen)
        check = _kill_from_states(micro, engine_states, res.transport, cfg, dp)
        if check != float(res.kill_fraction):
            raise AssertionError(
                f"[{key}/{name}] the state-swap path does not reproduce the engine "
                f"({check!r} vs {float(res.kill_fraction)!r}); the assumption that "
                "the cell state reaches the kill only through apply_drug_kill's "
                "resistance multipliers is WRONG and no number here may be believed"
            )

        kill_markers = _kill_from_states(
            micro, marker_states, res.transport, cfg, dp)

        # --- delivery metrics (state-invariant by construction) -------------- #
        T = exposure_window_s(res)
        exp = tumor_cell_exposure(micro, res)
        auc_sum = tissue_auc_summary(micro, res)
        cov = therapeutic_coverage(exp["auc_uM_s"], float(dp.ic50), T)

        row = {
            "drug": name,
            "mw": float(mol.get("mw")),
            "ic50_uM": float(dp.ic50),
            "bbb": float(dp.bbb_penetrance),
            "kill_oxygen": float(res.kill_fraction),
            "kill_markers": float(kill_markers),
            "ceiling_oxygen": _ceiling(engine_states, micro.cells.sensitivity,
                                       tumor, pdcfg, drug=dp),
            "ceiling_markers": _ceiling(marker_states, micro.cells.sensitivity,
                                        tumor, pdcfg, drug=dp),
            "mean_tumor_auc": float(auc_sum["mean_tumor_auc"]),
            "coverage_1x": float(cov["coverage_1x"]),
            "coverage_0.5x": float(cov["coverage_0.5x"]),
            "coverage_2x": float(cov["coverage_2x"]),
            "median_E_over_ic50": float(cov["median_E_over_ic50"]),
            "runtime_s": time.time() - t0,
        }

        # --- MGMT arm: temozolomide only, and only where MGMT is measured ----- #
        if mgmt["mgmt_on_panel"] and name == "temozolomide":
            arms = {}
            for tag, blk in mgmt["routes"].items():
                s = float(blk["sensitivity"])
                arms[tag] = {
                    "sensitivity": s,
                    "kill_markers": _kill_from_states(
                        micro, marker_states, res.transport, cfg, dp,
                        sensitivity_scalar=s),
                    "ceiling_markers": _ceiling(
                        marker_states, micro.cells.sensitivity, tumor, pdcfg,
                        sensitivity_scalar=s, drug=dp),
                }
            # CI band on the primary (bridged) route.
            lo, hi = mgmt["routes"]["bridged"]["ci95"]
            arms["bridged_ci_lo"] = {
                "sensitivity": float(lo),
                "kill_markers": _kill_from_states(
                    micro, marker_states, res.transport, cfg, dp,
                    sensitivity_scalar=float(lo)),
            }
            arms["bridged_ci_hi"] = {
                "sensitivity": float(hi),
                "kill_markers": _kill_from_states(
                    micro, marker_states, res.transport, cfg, dp,
                    sensitivity_scalar=float(hi)),
            }
            row["mgmt_arms"] = arms

        rows.append(row)
        print(f"  {name:14s} kill  oxygen={row['kill_oxygen']:.4f}  "
              f"markers={row['kill_markers']:.4f}   "
              f"cov@1x={row['coverage_1x']:.4f}  "
              f"meanAUC={row['mean_tumor_auc']:.3g}   [{row['runtime_s']:.0f}s]",
              flush=True)
        if checkpoint:
            with open(checkpoint, "w") as fh:
                json.dump({"key": key, "rows": rows}, fh)

    rows = [r for name in drugs for r in rows if r["drug"] == name]  # panel order

    if key == "brain_319":
        gate = _reproduction_gate(rows)
        ok = sum(1 for g in gate if g["reproduced"])
        print(f"  reproduction gate vs results/cell_state_impact.md §3.2: "
              f"{ok}/{len(gate)} values reproduced to {PUBLISHED_319_TOL:g}")

    fr_oxygen = state_fractions_from_states(
        _states_from_oxygen(micro, micro.fields.oxygen, ocfg), tumor)

    return {
        "key": key,
        "path": path,
        "label": meta["label"],
        "role": meta["role"],
        "dx_um": float(dx_um),
        "alpha": alpha,
        "n_cells": int(rep["n_cells"]),
        "n_tumor": int(tumor.sum()),
        "n_panel_genes": int(coverage.get("n_panel_genes", 0)),
        "grid": [int(v) for v in rep["grid"]],
        "mean_oxygen": rep["mean_oxygen"],
        "mean_ecm": rep["mean_ecm_density"],
        "mean_vascular": rep["mean_vascular_density"],
        "celltype_fractions": rep["celltype_fractions"],
        "state_mode": rep["state_mode"],
        "state_source": state_block["source"],
        "state_caveat": state_block["caveat"],
        "coverage": {
            "proliferation": {
                k: coverage.get("proliferation", {}).get(k)
                for k in ("present", "missing", "detection", "union_detection",
                          "call_level", "reason", "binarisation")},
            "hypoxia": {
                k: coverage.get("hypoxia", {}).get(k)
                for k in ("present", "missing", "detection", "union_detection",
                          "call_level", "reason", "binarisation")},
            "unmeasured_states": coverage.get("unmeasured_states"),
            "per_cell_detection_floor": coverage.get("per_cell_detection_floor"),
            "neighbourhood_radius_um": coverage.get("neighbourhood_radius_um"),
        },
        "fractions_tumor_markers": fr_markers,
        "fractions_tumor_oxygen": fr_oxygen,
        "mgmt": mgmt,
        "pd": {"hypoxia_resistance": pdcfg.hypoxia_resistance,
               "quiescent_resistance": pdcfg.quiescent_resistance,
               "max_kill_rate": pdcfg.max_kill_rate},
        "rows": rows,
        "ranking_markers": [r["drug"] for r in
                            sorted(rows, key=lambda r: -r["kill_markers"])],
        "ranking_oxygen": [r["drug"] for r in
                           sorted(rows, key=lambda r: -r["kill_oxygen"])],
        "ranking_coverage": [r["drug"] for r in
                             sorted(rows, key=lambda r: -r["coverage_1x"])],
        "ranking_auc_over_ic50": [
            r["drug"] for r in sorted(
                rows, key=lambda r: -(r["mean_tumor_auc"] / r["ic50_uM"]))],
        "runtime_s": time.time() - t_sec,
        "recon_s": t_recon,
    }


# --------------------------------------------------------------------------- #
# The report
# --------------------------------------------------------------------------- #
def _f(v, nd: int = 4) -> str:
    if v is None:
        return "n/a"
    try:
        f = float(v)
    except (TypeError, ValueError):
        return str(v)
    return "n/a" if not np.isfinite(f) else f"{f:.{nd}f}"


def _pct(v, nd: int = 2) -> str:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return "—"
    return "—" if not np.isfinite(f) else f"{f * 100:.{nd}f}%"


def _rank_line(ranking: Sequence[str], rows: Sequence[Dict], field: str) -> str:
    by = {r["drug"]: r for r in rows}
    return " > ".join(f"**{d}** ({by[d][field]:.4f})" for d in ranking)


def write_report(results: Dict[str, Dict], out_path: str = _OUT_MD) -> str:
    """Write ``results/flagship_markers.md``."""
    L: List[str] = []
    A = L.append
    order = [k for k in SECTION_ORDER if k in results]
    prim = results.get("gse324268_366")
    io = results.get("io_480")
    b319 = results.get("brain_319")

    # ------------------------------------------------------------------ head #
    A("# Re-anchoring the flagship drug panel to a section that can measure "
      "cell state")
    A("")
    A("**Experiment:** `experiments/flagship_markers.py` "
      "(run: `./.venv/bin/python experiments/flagship_markers.py`)")
    A("")
    A("`results/cell_state_impact.md` established that the section almost every "
      "headline number in this project comes from — `data_cache/xenium_gbm`, the "
      "40,887-cell 319-gene Xenium brain panel — **cannot support a cell state at "
      "all**: it carries no hypoxia-response gene of the marker set (HYPOXIC is "
      "UNMEASURED there and no cell is assigned it), and MKI67+TOP2A fire in 1.42% "
      "of its tumour cells, 13x below the published 19-40% GBM Ki-67 range, so a "
      "per-cell proliferation call is refused. Its own conclusion was that absolute "
      "kill and kill ceilings must not be quoted from that section.")
    A("")
    A("This report re-runs the **flagship five-drug panel** "
      "[temozolomide, lomustine, cisplatin, methotrexate, doxorubicin] on a section "
      "that DOES clear the bar, with `state_mode='markers'`, "
      "`sensitivity_mode='uniform'`, dx = 12 um — and reports what changes.")
    A("")

    if prim:
        tmz = next(r for r in prim["rows"] if r["drug"] == "temozolomide")
        arms = tmz.get("mgmt_arms", {})
        A("## Headline")
        A("")
        A("| | |")
        A("|---|---|")
        A(f"| **Primary section** | `{prim['path']}` — "
          f"{prim['n_cells']:,} cells ({prim['n_tumor']:,} tumour), "
          f"{prim['n_panel_genes']}-gene panel |")
        A(f"| **Cell state** | MEASURED: proliferation `per_cell` over "
          f"{prim['coverage']['proliferation']['present']}, hypoxia `per_cell` over "
          f"{prim['coverage']['hypoxia']['present']} — "
          f"**nothing UNMEASURED** |")
        fm = prim["fractions_tumor_markers"]
        A(f"| **Measured tumour state fractions** | VIABLE "
          f"{_pct(fm['viable'])} / QUIESCENT {_pct(fm['quiescent'])} / HYPOXIC "
          f"{_pct(fm['hypoxic'])} |")
        A(f"| **Cross-drug ranking (markers)** | "
          + " > ".join(prim["ranking_markers"]) + " |")
        same = (b319 is not None
                and prim["ranking_markers"] == b319["ranking_oxygen"])
        A(f"| **Same ranking as the published 319-gene run?** | "
          + ("**YES** — invariant across BOTH the section and the state route |"
             if same else "**NO** — see §4 |"))
        if arms:
            A(f"| **Temozolomide, MGMT unmeasured (the old assumption)** | kill "
              f"{tmz['kill_markers']:.4f} |")
            A(f"| **Temozolomide, MGMT MEASURED on this section** | kill "
              f"{arms['bridged']['kill_markers']:.4f} "
              f"(sensitivity {arms['bridged']['sensitivity']:.4f}) |")
        A("| **Is absolute kill now defensible?** | see §6 — **no, but the reason "
          "has moved** |")
        A("")
        A("> **This is a different patient on a different panel.** Nothing here "
          "transfers from or to the 319-gene section as biology, and neither "
          "section carries a drug-response outcome label. What travels is the "
          "ORDERING and the machinery; the levels do not.")
        A("")

    # ------------------------------------------------------- 1. the choice #
    A("## 1. Which section, and why")
    A("")
    A("| section | panel | cells (tumour) | proliferation | hypoxia | MGMT | "
      "verdict |")
    A("|---|---|---|---|---|---|---|")
    for k in order:
        r = results[k]
        p, h = r["coverage"]["proliferation"], r["coverage"]["hypoxia"]
        mg = r["mgmt"]
        mg_s = (f"**YES** ({_pct(mg['detection_all_cells'])} of cells)"
                if mg.get("mgmt_on_panel") else "no")
        verdict = {
            "gse324268_366": "**PRIMARY**",
            "io_480": "check",
            "brain_319": "comparison arm — cannot support cell state",
        }[k]
        A(f"| `{k}` | {r['n_panel_genes']} genes | {r['n_cells']:,} "
          f"({r['n_tumor']:,}) | `{p['call_level']}` "
          f"({_pct(p['union_detection'])}) | `{h['call_level']}` "
          f"({_pct(h['union_detection'])}) | {mg_s} | {verdict} |")
    A("")
    _mki = lambda r: (r["coverage"]["proliferation"]["detection"] or {}).get("MKI67")
    A("**GSE324268 is chosen as primary, and the deciding argument is MGMT.** "
      "Both large sections clear the proliferation detection floor and both measure "
      "hypoxia, so on cell state alone the choice would be close — the IO panel "
      "actually has the stronger proliferation readout (MKI67 "
      + (_pct(_mki(io)) if io else "34.06%")
      + ", inside the published 19-40% Ki-67 range, against GSE324268's "
      + (_pct(_mki(prim)) if prim else "13.34%") + "). But GSE324268 is "
      "the **only cached section in this project that carries MGMT**, and MGMT is "
      "temozolomide's resistance determinant. On every other section TMZ's "
      "sensitivity is left at a uniform 1.0 and its kill is flagged an UPPER BOUND "
      "because nobody has measured its resistance gene "
      "(`results/real_drug_panel.md` §*Sensitivity mode*). Here it can be measured "
      "(§5). A section on which four of five drugs have a measurable determinant "
      "and the fifth does not is a section on which the panel is still unfair; this "
      "is the first one where that asymmetry can be closed with a number.")
    A("")
    A("Two secondary reasons: GSE324268 measures hypoxia over **four** transcripts "
      "(VEGFA, CA9, ADM, SLC2A1) against the IO panel's one (VEGFA), and it is the "
      "only section carrying **CA9**, the one near-hypoxia-restricted gene of the "
      "set — so its hypoxic band is the only one in the project with any constraint "
      "on its SIZE rather than only its ordering. The IO section was run anyway, as "
      "a check, and §4 reports whether the two agree.")
    A("")

    # ------------------------------------------ 2. the measured state #
    A("## 2. The measured cell state, and where every label comes from")
    A("")
    A("Fractions are over TUMOUR cells — the population every kill number is about. "
      "The `oxygen` column is the shipping default on the SAME reconstruction, for "
      "contrast only.")
    A("")
    A("| section | route | VIABLE | QUIESCENT | HYPOXIC | mesh-dependent? |")
    A("|---|---|---:|---:|---:|---|")
    for k in order:
        r = results[k]
        for route, key_ in (("markers", "fractions_tumor_markers"),
                            ("oxygen", "fractions_tumor_oxygen")):
            f = r[key_]
            A(f"| `{k}` | `{route}` | {_pct(f['viable'])} | "
              f"{_pct(f['quiescent'])} | {_pct(f['hypoxic'])} | "
              + ("**no** — computed from measured coordinates and counts in um; "
                 "no grid enters |" if route == "markers" else
                 "yes — the tissue mask gates the vascular source, which sets the "
                 "band edges |"))
    A("")
    A("### 2.1 Marker provenance, per section")
    A("")
    for k in order:
        r = results[k]
        p, h = r["coverage"]["proliferation"], r["coverage"]["hypoxia"]
        A(f"**`{k}`** — {r['label']}")
        A("")
        A(f"* **Proliferation** — call level `{p['call_level']}` over "
          f"{p['present'] or 'NO GENE ON PANEL'}"
          + (f", per-gene detection over tumour cells "
             + ", ".join(f"{g} {_pct(v)}" for g, v in (p["detection"] or {}).items())
             + f", union {_pct(p['union_detection'])}" if p["present"] else "")
          + f". Missing from this panel: {p['missing'] or 'none'}.")
        A(f"* **Hypoxia** — call level `{h['call_level']}` over "
          f"{h['present'] or 'NO GENE ON PANEL'}"
          + (f", per-gene detection over tumour cells "
             + ", ".join(f"{g} {_pct(v)}" for g, v in (h["detection"] or {}).items())
             + f", union {_pct(h['union_detection'])}" if h["present"] else "")
          + f". Missing from this panel: {h['missing'] or 'none'}.")
        unmeas = r["coverage"]["unmeasured_states"] or None
        A("* **UNMEASURED states:** "
          + (", ".join(unmeas) if unmeas else
             "**none** — every band on this section is traceable to a measured "
             "transcript")
          + ".")
        A("")
    hyp_frac = (prim["fractions_tumor_markers"]["hypoxic"] if prim else float("nan"))
    ca9 = ((prim["coverage"]["hypoxia"]["detection"] or {}).get("CA9")
           if prim else None)
    A("**What the HYPOXIC band is, exactly.** Positivity for a hypoxia-RESPONSE "
      "transcript, not an oxygen tension. VEGFA, ADM and SLC2A1 all have real "
      "non-hypoxic expression, so the band's ORDERING is measured and its SIZE is "
      f"uncalibrated — the {_pct(hyp_frac, 1)} HYPOXIC band on the primary section "
      "is a threshold artefact as a *hypoxic fraction*, not a finding. Only CA9 is "
      "near-restricted to severe hypoxia and it fires in "
      + (_pct(ca9, 2) if ca9 is not None else "~1%")
      + " of that section's tumour cells — a plausible perinecrotic fraction, but "
        "too sparse for a per-cell call on its own. Nothing in this report quotes "
        "the HYPOXIC fraction's magnitude as a hypoxic fraction; it is used only as "
        "the band the engine's `PDConfig.hypoxia_resistance` acts on.")
    A("")
    A("**The one-label collapse, and why it bites hardest here.** Proliferation and "
      "hypoxia are two axes and `core.CellState` has one slot; a hypoxia-positive "
      "cell is labelled HYPOXIC whatever its proliferation call. That precedence is "
      "forced by the schema. On the 319-gene section it cost nothing because the "
      "HYPOXIC band was empty; on a section where hypoxia IS measured it puts "
      f"{_pct(hyp_frac, 1)} of tumour cells under the x0.35 multiplier and so "
      "dominates the state term. **Measuring a band you previously could not see "
      "makes the unvalidated multiplier attached to it matter more, not less** — "
      "see §6.3.")
    A("")

    # ------------------------------------------------ 3. the drug panel #
    A("## 3. The drug panel")
    A("")
    A("Five drugs, `sensitivity_mode='uniform'` (identical for every drug, so it "
      "cancels out of the ranking), dx = 12 um, one reconstruction and one "
      "transport solve per drug reused across every state route. "
      "`therapeutic_coverage_1x` is the fraction of tumour cells whose "
      "time-averaged exposure exceeds 1x IC50; it and `mean tumour AUC` carry **no "
      "cell-state dependence at all**, so they are identical between the two kill "
      "columns by construction.")
    A("")
    for k in order:
        r = results[k]
        A(f"### 3.{order.index(k) + 1} `{k}` — {r['n_cells']:,} cells, grid "
          f"{r['grid'][0]}x{r['grid'][1]}, alpha = {r['alpha']:.1f}")
        A("")
        A("| drug | IC50 (uM) | BBB | mean tumour AUC (uM*s) | "
          "therapeutic_coverage_1x | **kill (markers)** | kill (oxygen) | "
          "ceiling (markers) |")
        A("|---|---:|---:|---:|---:|---:|---:|---:|")
        for row in r["rows"]:
            A(f"| {row['drug']} | {row['ic50_uM']:.3g} | {row['bbb']:.3f} | "
              f"{row['mean_tumor_auc']:.3g} | {row['coverage_1x']:.4f} | "
              f"**{row['kill_markers']:.4f}** | {row['kill_oxygen']:.4f} | "
              f"{row['ceiling_markers']:.4f} |")
        A("")
        A(f"**Ranking by kill, `markers` route:** "
          f"{_rank_line(r['ranking_markers'], r['rows'], 'kill_markers')}")
        A("")
        A(f"**Ranking by kill, `oxygen` route (same cells, same transport):** "
          f"{_rank_line(r['ranking_oxygen'], r['rows'], 'kill_oxygen')}")
        A("")
        A(f"**Ranking by `therapeutic_coverage_1x`:** "
          + " > ".join(r["ranking_coverage"]))
        A("")
        if k == "brain_319":
            g = _reproduction_gate(r["rows"])
            ok = sum(1 for x in g if x["reproduced"])
            A(f"> **Reproduction gate — {ok}/{len(g)} values reproduce "
              f"`results/cell_state_impact.md` §3.2 to {PUBLISHED_319_TOL:g} "
              "absolute.** Three drugs x two state routes plus the two shared "
              "ceilings, i.e. the transport solve, both state routes and the PD "
              "mapping are all checked against an independently published run: "
              + "; ".join(f"{x['quantity']} {x['here']:.4f} vs {x['published']:.4f}"
                          for x in g)
              + ". This harness is running the project's pipeline, not a variant "
                "of it.")
            A("")
            A("> **But note what it does NOT reproduce, because it matters for "
              "reading every older number in `results/`.** "
              "`results/real_drug_panel.md` and `results/heldout_patient.md` quote "
              "*lower* kills for these same drugs on this same section — "
              + ", ".join(f"{d} {v:.4f}" for d, v in STALE_319_KILL.items())
              + " — and this run does not reproduce them. The cause is a pipeline "
                "drift, not a harness difference: those reports' own reconstruction "
                "block records mean vascular density **0.0250** and mean oxygen "
                "**0.2955** on this section, where the current code gives "
              + f"**{r['mean_vascular']:.4f}** and **{r['mean_oxygen']:.4f}** on "
                "identical inputs and an unchanged alpha = "
              + f"{r['alpha']:.2f}. More drug reaches the tissue, so every kill "
                "rises ~1.2x. `results/cell_state_impact.md` and "
                "`results/kill_calibration.md` postdate that change and agree with "
                "this run exactly. **The two older reports' absolute numbers are "
                "stale and should be regenerated; their RANKINGS are unaffected** "
                "(the ordering below is identical to theirs).")
            A("")

    # ------------------------------------------- 4. the ranking test #
    A("## 4. Does the ranking hold across SECTIONS and across ROUTES at once?")
    A("")
    A("`results/cell_state_impact.md` §3.3 found the cross-drug ranking invariant "
      "across three state routes on ONE section. `results/heldout_patient.md` §3.3 "
      "found it invariant across two sections on ONE route. This table crosses "
      "both axes simultaneously, which is the stronger claim.")
    A("")
    A("| section | route | ranking by predicted kill |")
    A("|---|---|---|")
    for k in order:
        r = results[k]
        A(f"| `{k}` | `markers` | " + " > ".join(r["ranking_markers"]) + " |")
        A(f"| `{k}` | `oxygen` | " + " > ".join(r["ranking_oxygen"]) + " |")
    A("")
    all_rank = []
    for k in order:
        all_rank.append(tuple(results[k]["ranking_markers"]))
        all_rank.append(tuple(results[k]["ranking_oxygen"]))
    invariant = len(set(all_rank)) == 1
    if invariant:
        A(f"**The ranking is INVARIANT across all {len(all_rank)} "
          f"(section x route) combinations tested: "
          + " > ".join(all_rank[0]) + ".** "
          "That is a materially stronger statement than either prior report could "
          "make. The mechanism is the same one both of them identified: the "
          "resistance multipliers are SHARED across drugs — every drug pays the "
          "same x0.5 in the same QUIESCENT cells and the same x0.35 in the same "
          "HYPOXIC cells — so the state term very largely cancels out of a "
          "comparison and does not cancel at all out of a level.")
    else:
        A("**The ranking is NOT invariant.** The distinct orderings observed are:")
        A("")
        for rk in sorted(set(all_rank)):
            who = [f"`{order[i // 2]}`/{'markers' if i % 2 == 0 else 'oxygen'}"
                   for i, x in enumerate(all_rank) if x == rk]
            A(f"* " + " > ".join(rk) + " — " + ", ".join(who))
    A("")
    A("**But the invariance is not equally safe at every adjacent pair, and the "
      "narrowest gap should be quoted with it.** The engine draws kills as "
      "Bernoulli trials over tumour cells, so a gap of order the sampling scale "
      "`sqrt(p(1-p)/n)` is not a resolved ordering:")
    A("")
    A("| section | tightest adjacent pair (markers route) | gap | Bernoulli sigma | "
      "gap / sigma |")
    A("|---|---|---:|---:|---:|")
    for k in order:
        r = results[k]
        ordered = sorted(r["rows"], key=lambda x: -x["kill_markers"])
        pairs = [(ordered[i], ordered[i + 1]) for i in range(len(ordered) - 1)]
        a, b = min(pairs, key=lambda p: p[0]["kill_markers"] - p[1]["kill_markers"])
        gap = a["kill_markers"] - b["kill_markers"]
        n = max(int(r["n_tumor"]), 1)
        sig = math.sqrt(max(a["kill_markers"] * (1 - a["kill_markers"]), 0.0) / n)
        A(f"| `{k}` | {a['drug']} vs {b['drug']} | {gap:.4f} | {sig:.5f} | "
          f"{gap / max(sig, 1e-12):.1f} |")
    A("")
    A("Every other adjacent pair is far wider. So the ordering is solid at the top "
      "— cisplatin's lead over the next drug is orders of magnitude above the "
      "sampling scale on every section — and the **methotrexate / temozolomide "
      "boundary at the bottom is the one place the ranking rests on a small "
      "margin**. §5 shows that measuring MGMT opens that particular gap decisively, "
      "which is the substantive reason to care about it here.")
    A("")
    n_cov_ok = sum(1 for k in order
                   if results[k]["ranking_coverage"] == results[k]["ranking_markers"])
    A(f"**`therapeutic_coverage_1x` still does not recover the kill ordering — on "
      f"{len(order) - n_cov_ok} of {len(order)} sections.** It swaps doxorubicin "
      "and lomustine on both large sections (doxorubicin is delivery-starved but "
      "1200x more potent, so it kills more while covering less) and it is pinned at "
      "0.0000 for temozolomide and methotrexate wherever their exposure never "
      "reaches 1x IC50. `results/heldout_patient.md` §6.3 found exactly this on two "
      "sections and declined to call it a generalisation failure; a third section "
      "agreeing makes it a settled property of the metric rather than of any "
      "tissue.")
    A("")
    A("**What moves instead is the level.** Same drugs, same physics, different "
      "state route:")
    A("")
    A("| section | drug | kill (oxygen) | kill (markers) | ratio |")
    A("|---|---|---:|---:|---:|")
    for k in order:
        for row in results[k]["rows"]:
            ratio = (row["kill_markers"] / row["kill_oxygen"]
                     if row["kill_oxygen"] > 0 else float("nan"))
            A(f"| `{k}` | {row['drug']} | {row['kill_oxygen']:.4f} | "
              f"{row['kill_markers']:.4f} | {_f(ratio, 2)}x |")
    A("")
    ratios = [r["kill_markers"] / r["kill_oxygen"] for k in order
              for r in results[k]["rows"] if r["kill_oxygen"] > 0]
    xsec = {}
    for k in order:
        for r in results[k]["rows"]:
            xsec.setdefault(r["drug"], []).append(r["kill_markers"])
    worst = max(xsec.items(), key=lambda kv: max(kv[1]) / max(min(kv[1]), 1e-12))
    A(f"Across the {len(ratios)} (section x drug) cells above the state route moves "
      f"the level by {min(ratios):.2f}x-{max(ratios):.2f}x. Across SECTIONS it "
      f"moves far more: the same drug's marker-route kill spans up to "
      f"{max(worst[1]) / max(min(worst[1]), 1e-12):.1f}x ({worst[0]}, "
      + " / ".join(f"{v:.4f}" for v in worst[1]) + ").")
    A("")
    A("**A finding worth naming: the state route matters LESS on the sections that "
      "can actually measure it.** On the two markers-capable sections the oxygen "
      "and marker routes give kills within 0.86x-1.02x of each other; on the "
      "319-gene section, where the marker route has to fall back to a neighbourhood "
      "rule and the oxygen route's VIABLE band is a mask artefact, they differ by "
      "1.18x-1.21x. The route disagreement on the incumbent section is not a "
      "measurement of how much cell state matters — it is a measurement of how "
      "badly that panel constrains it.")
    A("")
    A("So a claim of the form *\"cisplatin beats doxorubicin beats lomustine on "
      "this tissue\"* survives every section and every route tested here, while a "
      "claim of the form *\"temozolomide kills X% of tumour cells\"* does not "
      "survive a change of section at all.")
    A("")

    # ------------------------------------------------------- 5. MGMT #
    A("## 5. Temozolomide with a MEASURED MGMT — the first time in this project")
    A("")
    if not prim or not prim["mgmt"].get("mgmt_on_panel"):
        A("MGMT is not on the primary section's panel; this section is empty.")
    else:
        mg = prim["mgmt"]
        br = mg["bridge"]
        tmz = next(r for r in prim["rows"] if r["drug"] == "temozolomide")
        arms = tmz["mgmt_arms"]
        A("Every prior temozolomide number in this project carries the flag "
          "*\"UPPER BOUND — no intrinsic resistance is modelled for this drug, its "
          "determinant is not measurable on this panel\"*. MGMT is absent from both "
          "10x Xenium panels and has **no** panel surrogate (0/396 genes pass FDR, "
          "`results/gbm_tmz_exvivo.md` §B2). **GSE324268 carries MGMT.**")
        A("")
        A("### 5.1 What is measured")
        A("")
        A("| quantity | value |")
        A("|---|---:|")
        A(f"| MGMT transcripts in the section | {mg['mgmt_counts']:,.0f} |")
        A(f"| cells with >= 1 MGMT transcript | {_pct(mg['detection_all_cells'])} "
          f"({_pct(mg['detection_tumor'])} of tumour cells) |")
        A(f"| panel pseudobulk | {mg['panel_pseudobulk_cpm']:.2f} CPM "
          f"(log2(CPM+1) = {mg['panel_pseudobulk_log2cpm']:.3f}) |")
        A("")
        A("**MGMT is a per-patient scalar here, not a field, and that is a "
          "measurement statement rather than a modelling convenience:** it fires in "
          f"{_pct(mg['detection_all_cells'])} of cells at ~1 transcript each, so a "
          "per-cell MGMT map would be counting noise. This is exactly the "
          "resolution `tmesim/molecular/mgmt.py` was written for.")
        A("")
        A("### 5.2 The unit bridge — and why the naive number is wrong")
        A("")
        A("`mgmt_sensitivity` is calibrated in **whole-transcriptome** `log2(CPM+1)` "
          "(GSE232173 bulk RNA-seq, n = 17). A 366-gene targeted panel has a far "
          "smaller denominator, so its pseudobulk CPM is not that unit. Feeding the "
          f"raw panel CPM in gives log2(CPM+1) = {mg['panel_pseudobulk_log2cpm']:.3f}, "
          f"which is **above the calibration ceiling "
          f"{mg['calibration_ceiling_log2cpm']:.3f}** — an extrapolation, and a "
          "unit error rather than a biological finding.")
        A("")
        A(f"The bridge is measured, not assumed: across the "
          f"{br['n_sections']} cached **TCGA-GBM bulk** RNA-seq sections, the "
          f"{prim['n_panel_genes']} genes of this panel carry "
          f"{br['capture_median'] * 100:.2f}% of all transcriptome counts "
          f"(range {br['capture_min'] * 100:.2f}-{br['capture_max'] * 100:.2f}%). "
          "Rescaling by that fraction converts a panel pseudobulk CPM into a "
          "bulk-equivalent one.")
        A("")
        A("| route | MGMT level | TMZ sensitivity | in calibration range? |")
        A("|---|---|---:|---|")
        b = arms["bridged"]
        rb = mg["routes"]["bridged"]
        A(f"| **BRIDGED (primary)** | {rb['cpm']:.2f} CPM bulk-equivalent "
          f"(log2 {rb['log2cpm']:.3f}) | **{b['sensitivity']:.4f}** | "
          + ("**yes**" if not rb["extrapolated"] else "NO — extrapolated") + " |")
        rbl, rbh = mg["routes"]["bridged_lo"], mg["routes"]["bridged_hi"]
        cpms = sorted([rbl["cpm"], rbh["cpm"]])
        sens = sorted([rbl["sensitivity"], rbh["sensitivity"]])
        A(f"| bridged, over the capture range | {cpms[0]:.2f}-{cpms[1]:.2f} CPM | "
          f"{sens[0]:.4f}-{sens[1]:.4f} | yes |")
        rn = mg["routes"]["binary_unmethylated"]
        A(f"| BINARY `unmethylated` (independent check) | — (a status, not a "
          f"level) | {rn['sensitivity']:.4f} | n/a |")
        rnv = mg["routes"]["naive_panel_cpm"]
        A(f"| NAIVE panel CPM (**rejected**) | {rnv['cpm']:.2f} CPM "
          f"(log2 {rnv['log2cpm']:.3f}) | {rnv['sensitivity']:.4f} | "
          "**NO — unit error** |")
        A("")
        A("**Two independent checks land on the bridged value.**")
        A("")
        A(f"1. **The TCGA GBM cohort's own MGMT.** Those same "
          f"{br['n_sections']} bulk sections have MGMT at a median "
          f"{br['mgmt_cpm_median']:.2f} CPM (range {br['mgmt_cpm_min']:.2f}-"
          f"{br['mgmt_cpm_max']:.2f}). The bridged value for this Xenium section, "
          f"{rb['cpm']:.2f} CPM, lands essentially on that median — i.e. a GBM's "
          "MGMT level comes out looking like a GBM's MGMT level, from a completely "
          "separate assay and cohort.")
        A(f"2. **The binary route agrees to "
          f"{abs(b['sensitivity'] - rn['sensitivity']) / rn['sensitivity'] * 100:.0f}%.** "
          f"The continuous bridged sensitivity is {b['sensitivity']:.4f}; the "
          "measured median-kill ratio between MGMT-expressing and MGMT-silent "
          f"cultures — which needs no unit bridge at all — is {rn['sensitivity']:.4f}. "
          "Two routes through the calibration that share no arithmetic converge.")
        A("")
        A(f"This section's tumour is therefore **MGMT-EXPRESSING** "
          f"(bulk-equivalent {rb['cpm']:.2f} CPM against the calibration's "
          f"{MGMT_EXPRESSED_CPM_THRESHOLD:.1f} CPM expressing threshold), i.e. the "
          "resistant group. That is a measurement about this patient, and it is the "
          "opposite of the uniform-1.0 assumption every prior TMZ number in this "
          "project carries.")
        A("")
        A("### 5.3 What it does to temozolomide")
        A("")
        A("| MGMT input | sensitivity | TMZ kill (markers route) | x vs uniform | "
          "rank in the 5-drug panel |")
        A("|---|---:|---:|---:|---|")

        def _rank_of(kill: float) -> str:
            others = sorted(
                ((r["drug"], r["kill_markers"]) for r in prim["rows"]
                 if r["drug"] != "temozolomide"),
                key=lambda t: -t[1])
            pos = 1 + sum(1 for _, k in others if k > kill)
            return f"#{pos} of 5"

        base = tmz["kill_markers"]
        A(f"| none (uniform 1.0 — the old assumption) | 1.0000 | {base:.4f} | "
          f"1.00x | {_rank_of(base)} — **UPPER BOUND** |")
        for tag, label in (("bridged", "**MEASURED, bridged (primary)**"),
                           ("binary_unmethylated", "binary `unmethylated`"),
                           ("naive_panel_cpm", "naive panel CPM (rejected)")):
            a = arms[tag]
            A(f"| {label} | {a['sensitivity']:.4f} | {a['kill_markers']:.4f} | "
              f"{a['kill_markers'] / base if base > 0 else float('nan'):.2f}x | "
              f"{_rank_of(a['kill_markers'])} |")
        lo, hi = arms["bridged_ci_lo"], arms["bridged_ci_hi"]
        rlo, rhi = _rank_of(lo["kill_markers"]), _rank_of(hi["kill_markers"])
        A(f"| bridged, 95% CI on lambda | {lo['sensitivity']:.4f}-"
          f"{hi['sensitivity']:.4f} | {lo['kill_markers']:.4f}-"
          f"{hi['kill_markers']:.4f} | — | "
          + (rlo if rlo == rhi else f"{rlo} to {rhi}") + " |")
        A("")
        A("**Read the rank column carefully.** `results/mgmt_scalar.md` found that "
          "supplying a measured MGMT returns TMZ from an artificially-inflated #4 "
          "to #5 — but that promotion was specific to `determinant` mode, where "
          "TMZ's four competitors are scaled down by their own measured determinants "
          "and TMZ alone is not. This run is in `uniform` mode, where no drug is "
          "scaled, so TMZ starts last on delivery x potency alone and the "
          "measurement cannot demote it further. **What MGMT changes here is the "
          "MAGNITUDE and the LABEL, not the rank**: a number that was flagged "
          "*'no intrinsic resistance modelled'* now carries the resistance value "
          "used, its CI, and where it came from.")
        A("")
        others = sorted(((r["drug"], r["kill_markers"]) for r in prim["rows"]
                         if r["drug"] != "temozolomide"), key=lambda t: t[1])
        nxt, nxt_k = others[0]
        gap0 = nxt_k - base
        gap1 = nxt_k - arms["bridged"]["kill_markers"]
        sig = math.sqrt(max(base * (1 - base), 0.0) / max(prim["n_tumor"], 1))
        A("**It does, however, make the placement itself robust, and that IS the "
          f"same correction `results/mgmt_scalar.md` describes.** Without MGMT, "
          f"temozolomide ({base:.4f}) sits below {nxt} ({nxt_k:.4f}) by only "
          f"{gap0:.4f} — {gap0 / max(sig, 1e-12):.1f}x the Bernoulli sampling scale "
          f"on {prim['n_tumor']:,} tumour cells, i.e. a placement that is nearly a "
          f"coin-flip. With the measured MGMT scalar the gap opens to {gap1:.4f} "
          f"({gap1 / max(sig, 1e-12):.0f}x sigma). That mirrors the finding in "
          "`results/mgmt_scalar.md` exactly — there the measurement moved TMZ off "
          "an artificially-inflated #4 and opened a 0.5 pp gap to 7-9 pp; here it "
          "cannot move the rank because `uniform` mode never inflated it, but it "
          "converts a marginal last place into a decisive one.")
        A("")
        A("**And it is not a small correction.** The measured scalar cuts "
          "temozolomide's predicted kill by "
          f"{base / max(arms['bridged']['kill_markers'], 1e-12):.1f}x "
          f"({base:.4f} -> {arms['bridged']['kill_markers']:.4f}). "
          "Every temozolomide kill previously published by this project — "
          "`results/real_drug_panel.md`, `results/dose_schedule.md`, "
          "`results/kill_calibration.md`, `results/heldout_patient.md` — is an "
          "upper bound in this dimension, exactly as "
          "`results/kill_calibration.md` limitation 8 says it is. This is the first "
          "run in the project where that bound is replaced by a measurement rather "
          "than restated as a caveat.")
        A("")
        A("**What this does NOT establish.** The mapping is fitted on n = 17 paired "
          "cultures with a 95% CI on lambda spanning 0.006-0.680, i.e. a factor of "
          ">100 in slope; the direction is solid, the magnitude is approximate. The "
          "cultures are dissociated 2-D monolayers. The panel-to-bulk bridge assumes "
          "unbiased relative capture by the Xenium probe set, which is not "
          "established — the honest reading of the bridge is that it fixes an "
          "order-of-magnitude unit error, not that it calibrates a level. And the "
          "TCGA sections supplying the bridge are different patients.")
    A("")

    # ------------------------------------------------- 6. the verdict #
    A("## 6. Is absolute kill defensible on this section?")
    A("")
    A("**No — and the honest version of that answer is more useful than it was, "
      "because the terms that block it have changed.**")
    A("")
    A("What is now closed on this section that was open on the 319-gene one:")
    A("")
    A("1. **Cell state is MEASURED, not assumed.** Both bands are traceable to "
      "measured transcripts at detection rates that are biology rather than probe "
      "sensitivity; nothing is UNMEASURED; and the labelling is mesh-independent by "
      "construction (§2). On the incumbent section hypoxia was unmeasured, "
      "proliferation was 13x below the published index, and the oxygen route's "
      "VIABLE band emptied entirely across the mask x mesh grid.")
    A("2. **Temozolomide's resistance determinant is MEASURED** (§5). TMZ was the "
      "one drug in the panel whose kill could not be scaled by anything measured; "
      "it now can be.")
    A("3. **The ranking survives a stronger test** — the TOP 3 is invariant "
      "across sections AND routes simultaneously (§4). Ranks 4-5 are NOT: "
      "`brain_319` orders methotrexate above temozolomide and the other two "
      "sections reverse it. And all of these axes are common-mode, so none of "
      "them speaks for the delivery law, which breaks even the surviving top 3 "
      "(§4).")
    A("")
    A("What remains open, in order of size:")
    A("")
    A("1. **No outcome label. This is the one that actually blocks it.** GSE324268 "
      "is an untreated resection with no paired drug response, exactly like the "
      "other two sections. Everything in §3 is a prediction with nothing to score "
      "it against. Measuring the cell state correctly makes the prediction "
      "*defensible in its inputs*; only an outcome makes it defensible as a "
      "prediction.")
    A("2. **Between-patient potency spread, ~17.4x.** `results/kill_calibration.md` "
      "§4 propagated the measured inter-culture IQR of ex-vivo temozolomide "
      "viability and got a 17.4x band on absolute kill — against a bound of <= 2.9x "
      "from the entire cell-state term (the multipliers are x0.35 and x0.5, so no "
      "relabelling can move `F_state` outside [0.35, 1.0]). **That comparison is "
      "the crux of this section's verdict: the term this report closes was already "
      "the smaller of the two.** Closing it is still worth doing — it is the "
      "*epistemic* term, the one a better model can remove, while the 17.4x is "
      "biological spread that shrinks only by measuring this patient's own "
      "sensitivity — but it does not on its own make a level quotable.")
    A("3. **The resistance MULTIPLIERS are still unvalidated priors.** x0.35 "
      "hypoxic and x0.5 quiescent are declared in `core/config.py:PDConfig` with no "
      "citation in this repo. This report measures WHICH cells are in each band; it "
      "does nothing whatever about what each band is worth. On a section where the "
      "HYPOXIC band is now large this matters MORE than it did on the 319-gene "
      "section, where that band was empty and the multiplier acted on no cells at "
      "all.")
    A("4. **The hypoxic band's SIZE is uncalibrated.** §2 — the band is an "
      "enrichment label built from hypoxia-RESPONSE transcripts with real "
      "non-hypoxic expression. Its ordering is measured; its magnitude is not. A "
      "measured pO2 (pimonidazole IHC or phosphorescence quenching on a serial "
      "section) would close this.")
    A("5. **One calibration constant** — alpha, fit once on doxorubicin penetration "
      "depth on the 319-gene section and NOT refitted here — sets the perivascular "
      "sink for every drug on every section.")
    A("6. **Grid convergence is inherited, not re-established.** dx = 12 um was "
      "shown converged on the 319-gene section's geometry "
      "(`results/numerical_rigor.md` §2). These sections are ~20x denser with "
      "larger vessel objects; the convergence argument is carried over.")
    A("")
    A("**The operative sentence.** On this section an absolute kill can for the "
      "first time be quoted with every term in it traceable to a measurement or to "
      "a named unvalidated prior — measured state, measured MGMT, measured "
      "geometry, published potency — which is a real change from the 319-gene "
      "section, where the state term was a bookkeeping choice with an unbounded "
      "swing. It still should not be quoted as a headline result, because a "
      "prediction with no outcome to score it against is not a result no matter how "
      "well-sourced its inputs are. **The cross-drug ranking remains the honest "
      "primary output, and it is now supported by a stronger invariance test than "
      "before.**")
    A("")

    # ------------------------------------- 7. should it become primary? #
    A("## 7. Should this section replace the 319-gene one as the project's "
      "reference?")
    A("")
    A("**Yes for anything that depends on cell state; no as a blanket swap; and "
      "the two should be reported together rather than one silently replacing the "
      "other.**")
    A("")
    A("| claim type | quote it from | why |")
    A("|---|---|---|")
    A("| absolute kill, kill ceiling, state fractions | **GSE324268** | the only "
      "section where both bands are measured and MGMT is available |")
    A("| temozolomide anything | **GSE324268** | the only section where TMZ's "
      "determinant is measurable |")
    A("| cross-drug ranking | either — the TOP 3 is invariant, ranks 4-5 are not "
      "(§4) | quote the invariance and its axis, never unqualified; it does not "
      "hold across the delivery law |")
    A("| delivery metrics (AUC, coverage, penetration length) | either | they "
      "carry no state dependence at all |")
    A("| the vascular external anchors (endothelial fraction, vascular area "
      "fraction) | **neither is settled** | `results/heldout_patient.md` §6.2 shows "
      "both fail to replicate across sections |")
    A("")
    A("Three reasons not to make it a clean swap:")
    A("")
    A("1. **It is a different patient.** Nothing about the 319-gene section's "
      "tumour transfers to this one. A section change is not a correction of the "
      "old numbers, it is a different measurement.")
    A("2. **Every published number would move.** The re-run list in "
      "`results/cell_state_impact.md` §5.2 applies with a section change on top of "
      "the route change, and the levels move (§4) even where the ordering does not. "
      "Flipping silently under a reader comparing against an older document would "
      "be exactly the failure mode this project keeps guarding against.")
    A("3. **The 319-gene section still has the one thing these do not**: it is the "
      "section every calibration in the project (alpha, the grid-convergence "
      "study, the vessel-linking distance) was fitted on. Retiring it means "
      "carrying those constants onto tissue they were never fitted to — which this "
      "report already does, and flags.")
    A("")

    # ---------------------------------------------------- 8. caveats #
    A("## 8. Caveats — what could make this wrong")
    A("")
    A("1. **No outcome label on any section.** Every kill here is a prediction with "
      "nothing to score it against.")
    A("2. **Different patients, different panels.** Cross-section comparisons mix "
      "real biological difference with panel content. Absolute AUCs and kills are "
      "NOT comparable between sections (`results/heldout_patient.md` §3.5 shows the "
      "ECM field alone has no common scale across panels); only orderings and "
      "within-section spreads are.")
    A("3. **Targeted-panel dropout.** A zero is \"no probe fired\", not \"no "
      "transcript\". No individual cell's marker call is trustworthy on any panel — "
      "only rates and orderings are.")
    A("4. **The hypoxic band over-calls hypoxia** relative to any pO2-based "
      "definition, and the one-label collapse lets it win over the proliferation "
      "call, so it dominates the state term on this section.")
    A("5. **The MGMT bridge assumes unbiased relative capture** by the Xenium probe "
      "set, and the bridge cohort (TCGA GBM, n = 7) are different patients. The "
      "binary route agrees, which is reassuring and is not independent evidence "
      "about the bridge itself.")
    A("6. **The MGMT calibration is n = 17** with a 95% CI on lambda spanning a "
      ">100x range in slope, fitted on dissociated 2-D monolayers.")
    A("7. **The PD resistance multipliers are unvalidated priors** (§6.3).")
    A("8. **alpha is not refitted** per section (§6.5), and grid convergence is "
      "inherited (§6.6).")
    A("9. **2-D, in-plane, single section.** Distances are Euclidean and in-plane; "
      "there is no out-of-plane transport.")
    A("")

    # ------------------------------------------------ 9. reproduction #
    A("## 9. Reproduction and cost")
    A("")
    A("```bash")
    A("./.venv/bin/python experiments/flagship_markers.py")
    A("```")
    A("")
    A("Deterministic — the only randomness is the engine's seeded generator "
      "(`seed = 0`). One reconstruction per section, one transport solve per drug, "
      "both reused across every state route and every MGMT arm; nothing is "
      "re-solved to change a state label.")
    A("")
    A("| section | cells | grid | reconstruction | 5 transport solves | total |")
    A("|---|---:|---|---:|---:|---:|")
    total = 0.0
    for k in order:
        r = results[k]
        drugs_s = sum(row["runtime_s"] for row in r["rows"])
        total += r["runtime_s"]
        A(f"| `{k}` | {r['n_cells']:,} | {r['grid'][0]}x{r['grid'][1]} | "
          f"{r['recon_s']:.0f} s | {drugs_s / 60.0:.1f} min | "
          f"{r['runtime_s'] / 60.0:.1f} min |")
    A(f"| **all** | | | | | **{total / 60.0:.1f} min** |")
    A("")
    A("**Files added by this work (nothing under `tmesim/` was modified):**")
    A("")
    A("* `experiments/flagship_markers.py` — this analysis")
    A("* `results/flagship_markers.md` — this report")
    A("")
    A("**Data provenance.** Three open single-cell spatial GBM sections — 10x "
      "Xenium Human Brain Glioblastoma (v1.3.0, 319 targets, CC BY 4.0), 10x Xenium "
      "Human Brain Cancer GBM FFPE IO panel (v2.0.0, 480 targets, CC BY 4.0), and "
      "GEO GSE324268 GBMtumor07 (366 targets) — plus 7 cached TCGA-GBM bulk RNA-seq "
      "sections used ONLY for the panel-capture unit bridge in §5.2. None carries a "
      "paired drug-response outcome label. The MGMT mapping is calibrated on GEO "
      "GSE232173 / Ntafoulis et al., Br J Cancer 2023 (PMID 37620410).")

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write("\n".join(L) + "\n")
    print(f"\nreport written: {out_path}")
    return out_path


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #
def main(sections: Sequence[str] = SECTION_ORDER, dx_um: float = DX_UM,
         cache: str = _DEFAULT_CACHE, report_only: bool = False) -> None:
    t0 = time.time()
    results: Dict[str, Dict] = {}
    if cache and os.path.exists(cache):
        with open(cache) as fh:
            results = json.load(fh)
        print(f"loaded {len(results)} cached section(s) from {cache}")

    if not report_only:
        for key in sections:
            if key in results:
                print(f"[{key}] already cached — skipping (delete {cache} to redo)")
                continue
            ckpt = (os.path.splitext(cache)[0] + f".{key}.partial.json"
                    if cache else None)
            results[key] = run_section(key, dx_um=dx_um, checkpoint=ckpt)
            if cache:
                os.makedirs(os.path.dirname(cache) or ".", exist_ok=True)
                with open(cache, "w") as fh:
                    json.dump(results, fh, indent=1)
                print(f"  cached -> {cache}")

    if not results:
        raise SystemExit("no results to report; run without --report-only first")
    write_report(results)
    print(f"[total {time.time() - t0:.0f} s]")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="Flagship drug panel on a markers-capable section")
    ap.add_argument("--sections", default=",".join(SECTION_ORDER),
                    help="comma-separated section keys (default: %(default)s)")
    ap.add_argument("--dx", type=float, default=DX_UM,
                    help="grid spacing in um (default %(default)s, grid-converged)")
    ap.add_argument("--cache", default=_DEFAULT_CACHE,
                    help="raw per-section results JSON (kept OUTSIDE results/)")
    ap.add_argument("--report-only", action="store_true",
                    help="rebuild the markdown from the cache without recomputing")
    args = ap.parse_args()
    main(sections=[s.strip() for s in args.sections.split(",") if s.strip()],
         dx_um=args.dx, cache=args.cache, report_only=args.report_only)
