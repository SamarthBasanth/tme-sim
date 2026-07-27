"""Whole-niche REAL microenvironment assembler (the flagship reconstruction).

This is the single entry point that turns a REAL spatial-transcriptomics sample
into an assembled, physics-ready :class:`~tmesim.core.Microenvironment` with NO
synthetic data anywhere in the path (see ``docs/REAL_BUILD_PLAN.md``). It wires
together the Wave-2 real modules:

  * cells + phenotypes  — :func:`reconstruct.real_xenium.reconstruct_from_xenium`
    (single-cell Xenium, threshold-gated marker typing) or
    :func:`reconstruct.real_visium.reconstruct_from_visium` (spot deconvolution);
  * REAL vascular field — KDE of measured ENDOTHELIAL cells
    (:func:`reconstruct.real_fields.vascular_density_from_cells`);
  * REAL ECM field      — for Xenium, an EXPRESSION-WEIGHTED KDE of the panel's
    real collagen / matrix transcripts (COL1A2, POSTN, VCAN, DCN, LAMA2, THBS1,
    ...) via :func:`reconstruct.real_fields.ecm_density_from_positions`; for
    Visium (no per-cell transcripts here) the stromal-fraction proxy that
    ``reconstruct_from_visium`` already builds;
  * REAL-perfusion oxygen — the unmodified physics solver driven by the real
    vascular field (:func:`reconstruct.real_fields.add_oxygen`).

Every field derives from a measured byte on disk. The per-cell ``sensitivity``
carried through from the Xenium reconstruction is PROVISIONAL (a marker-scaled
placeholder, not a measured drug response) — this open section has no outcome
label. Nothing here modifies core / engine / physics / molecular.
"""
from __future__ import annotations

import time
from typing import Dict, List, Optional

import numpy as np

from ..core import (
    CellPopulation,
    CellType,
    Microenvironment,
    OxygenConfig,
    TissueFields,
)
from ..data.real.singlecell_spatial import load as load_singlecell_spatial
from .real_fields import (
    add_oxygen,
    ecm_density_from_cells,
    ecm_density_from_positions,
    vascular_density_from_cells,
)
from .real_visium import reconstruct_from_visium
from .real_xenium import reconstruct_from_xenium

# REAL extracellular-matrix / matricellular transcripts carried on the Xenium
# Human Brain panel. The ECM field is an expression-weighted KDE of these (only
# the ones actually present in a given panel are used at runtime).
ECM_MARKERS: List[str] = [
    "COL1A2", "COL12A1", "COL25A1", "POSTN", "VCAN",
    "DCN", "FBLN1", "LAMA2", "THBS1", "SPON1", "TGFBI",
]

_DEFAULT_VISIUM_SAMPLE = "10x_parent_gbm"


def _load_xenium_retry(sample: Optional[str]):
    """Load the real Xenium data, retrying once on any transient/IO error."""
    kwargs = {"dense": False}
    if sample is not None:
        kwargs["path"] = sample
    try:
        return load_singlecell_spatial(**kwargs)
    except Exception:
        time.sleep(3)
        return load_singlecell_spatial(**kwargs)


def _ecm_weights_from_expression(expression, genes: np.ndarray):
    """Per-cell REAL ECM abundance = summed log1p counts of the panel's matrix
    genes. Returns ``(weights (n,), used_markers)``; weights are all-zero (and the
    caller should fall back to the stromal proxy) if no ECM marker is on the panel.
    """
    genes = np.asarray(genes, dtype=str)
    gene_index = {g: j for j, g in enumerate(genes)}
    used = [g for g in ECM_MARKERS if g in gene_index]
    if not used:
        return np.zeros(expression.shape[0], dtype=np.float64), used
    cols = [gene_index[g] for g in used]
    from scipy import sparse

    if sparse.issparse(expression):
        sub = expression[:, cols].toarray()
    else:
        sub = np.asarray(expression)[:, cols]
    weights = np.log1p(sub.astype(np.float64)).sum(axis=1)
    return weights, used


def _reconstruct_xenium(sample: Optional[str], dx_um: float) -> Microenvironment:
    """Assemble the real single-cell Xenium microenvironment with REAL fields.

    Cells + threshold-gated phenotypes + domain come from
    :func:`reconstruct_from_xenium` (the canonical, fixed-typing cell
    reconstructor). We then rebuild the physical fields as REAL fields: vascular
    from the measured endothelial cells, ECM from the real matrix-transcript KDE,
    and oxygen from the physics solver driven by that real vascular field.
    """
    micro = reconstruct_from_xenium(sample_path=sample, dx_um=dx_um)
    domain, cells = micro.domain, micro.cells

    # REAL vascular density: KDE of the measured ENDOTHELIAL cell positions.
    vascular = vascular_density_from_cells(domain, cells)

    # REAL ECM: expression-weighted KDE of the panel's collagen / matrix genes.
    # Load once more (aligned row-for-row to the deterministic reconstruction) to
    # read the per-cell matrix-transcript abundance the CellPopulation drops.
    data = _load_xenium_retry(sample)
    weights, ecm_markers = _ecm_weights_from_expression(
        data["expression"], data["genes"]
    )
    if weights.size == cells.n and float(weights.sum()) > 0.0:
        ecm = ecm_density_from_positions(domain, cells.x, cells.y, weights=weights)
        ecm_source = (
            "expression-weighted KDE of REAL matrix transcripts "
            f"({', '.join(ecm_markers)})"
        )
    else:
        # Panel carries no ECM transcript (or row mismatch): stromal-cell proxy.
        ecm = ecm_density_from_cells(domain, cells)
        ecm_source = "STROMAL-cell density proxy (no ECM transcript available)"

    micro.fields = TissueFields(
        ecm_density=ecm, vascular_density=vascular, oxygen=None
    )
    add_oxygen(micro, OxygenConfig())

    micro.metadata = dict(micro.metadata)
    micro.metadata["real_pipeline"] = {
        "source": "xenium",
        "vascular_source": "KDE of REAL ENDOTHELIAL cell positions",
        "ecm_source": ecm_source,
        "ecm_markers_used": ecm_markers,
        "oxygen_source": "physics.solve_oxygen driven by the real vascular field",
    }
    return micro


def _reconstruct_visium(sample: Optional[str], dx_um: float) -> Microenvironment:
    """Assemble the real spot-level Visium microenvironment and add oxygen.

    ``reconstruct_from_visium`` already builds vascular density (from the real
    ENDOTHELIAL spot fraction) and an ECM proxy (from the real STROMAL fraction);
    per-cell collagen transcripts are not available at spot resolution, so the
    stromal proxy is the honest ECM here. We only add the real-perfusion oxygen.
    """
    sample_id = sample or _DEFAULT_VISIUM_SAMPLE
    # Visium spot pitch is ~100 um; a fine (xenium-style) dx would over-resolve
    # the spots, so fall back to the loader default when dx is left small.
    visium_dx = float(dx_um) if dx_um and dx_um >= 50.0 else 100.0
    micro = reconstruct_from_visium(sample_id, dx_um=visium_dx)
    if micro.fields.oxygen is None:
        add_oxygen(micro, OxygenConfig())
    micro.metadata = dict(micro.metadata)
    micro.metadata["real_pipeline"] = {
        "source": "visium",
        "vascular_source": "KDE of REAL ENDOTHELIAL spot fraction",
        "ecm_source": "STROMAL-fraction proxy (spot resolution; no per-cell ECM transcript)",
        "ecm_markers_used": [],
        "oxygen_source": "physics.solve_oxygen driven by the real vascular field",
    }
    return micro


def reconstruct_real(
    source: str = "xenium",
    sample: Optional[str] = None,
    dx_um: float = 30.0,
) -> Microenvironment:
    """Reconstruct a whole-niche REAL microenvironment from open spatial data.

    Parameters
    ----------
    source : {'xenium', 'visium'}
        'xenium' -> single-cell-resolution 10x Xenium GBM section (threshold-gated
        marker typing, expression-weighted ECM). 'visium' -> spot-level 10x Visium
        GBM (reference-based deconvolution, stromal-proxy ECM).
    sample : str or None
        Xenium: cache directory (None -> loader default). Visium: sample id (None
        -> ``'10x_parent_gbm'``).
    dx_um : float
        Continuum grid spacing (um). Default 30 um suits the single-cell Xenium
        section; for Visium a value < 50 falls back to the 100 um spot pitch.

    Returns
    -------
    Microenvironment
        Real cells at real positions, REAL vascular / ECM / oxygen fields, ready
        for :func:`tmesim.sim.engine.simulate_response`.
    """
    src = str(source).strip().lower()
    if src == "xenium":
        return _reconstruct_xenium(sample, dx_um)
    if src == "visium":
        return _reconstruct_visium(sample, dx_um)
    raise ValueError(f"unknown source {source!r}; expected 'xenium' or 'visium'")


def reconstruction_report(micro: Microenvironment) -> Dict:
    """Summarise a reconstructed real microenvironment: cell-type fractions, mean
    vascular / ECM / oxygen (over the tissue mask), and the grid."""
    ct = np.asarray(micro.cells.cell_type)
    n = int(micro.cells.n)
    counts = {c.name: int((ct == int(c)).sum()) for c in CellType}
    fractions = {
        c.name: (counts[c.name] / n if n else 0.0) for c in CellType
    }

    tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)

    def _mean(field: Optional[np.ndarray]):
        if field is None:
            return None
        field = np.asarray(field, dtype=np.float64)
        return float(field[tissue].mean()) if tissue.any() else float(field.mean())

    f = micro.fields
    return {
        "patient_id": micro.patient_id,
        "source": micro.metadata.get("real_pipeline", {}).get("source"),
        "n_cells": n,
        "n_tumor": int((ct == int(CellType.TUMOR)).sum()),
        "grid": tuple(micro.domain.shape),
        "dx_um": float(micro.domain.dx),
        "extent_um": tuple(float(v) for v in micro.domain.extent_um),
        "celltype_counts": counts,
        "celltype_fractions": fractions,
        "mean_vascular_density": _mean(f.vascular_density),
        "mean_ecm_density": _mean(f.ecm_density),
        "mean_oxygen": _mean(f.oxygen),
        "ecm_markers_used": micro.metadata.get("real_pipeline", {}).get(
            "ecm_markers_used", []
        ),
    }


if __name__ == "__main__":
    micro = reconstruct_real("xenium")
    report = reconstruction_report(micro)
    print("=== REAL whole-niche microenvironment (xenium) ===")
    for k, v in report.items():
        print(f"  {k}: {v}")
