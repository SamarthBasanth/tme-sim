"""Reconstruct a REAL cell-resolved microenvironment from the Xenium GBM section.

This is the flagship *real* reconstruction (see ``docs/REAL_BUILD_PLAN.md``): it
turns the real 40,887-cell 10x Xenium single-cell-resolution glioblastoma
section into a real, cell-resolved :class:`tmesim.core.Microenvironment`. There
is NO synthetic data anywhere in this path — every cell position and every
expression count is a measured byte on disk
(``tmesim.data.real.singlecell_spatial.load``).

What is real vs. provisional
----------------------------
REAL (measured):
  * per-cell (x, y) centroids in micrometers;
  * per-cell targeted-panel counts for the 319-gene Xenium Human Brain panel;
  * the marker-gene phenotype of each cell (derived only from that real signal).
PROVISIONAL (a documented modelling choice, NOT a real outcome):
  * the per-cell ``sensitivity`` in [0,1] — a placeholder scaled from the tumour
    marker signal. It is NOT a measured drug response; this open section carries
    no paired outcome label. Do not treat it as ground truth.
  * ``ecm_density`` / ``vascular_density`` — Gaussian-smoothed spatial densities
    of the stromal / endothelial cell phenotypes (a proxy field, not a stain).

Phenotyping (marker -> CellType) — THRESHOLD-GATED
--------------------------------------------------
Markers are log1p-transformed then z-scored across all cells (so each marker
contributes on a comparable scale), and each type's score is the *mean* z-score
over its present markers (mean, not sum, so types with different marker counts
compare on a comparable scale). A cell is assigned a type only when that type's
score both (a) clears a POSITIVE threshold :data:`TYPING_THRESHOLD` (z >= 1.0, a
marker signal at least ~1 SD above the cohort mean — a genuine positive call, not
merely noise) AND (b) is the maximum among the gated types. Cells that clear no
gate default to ``TUMOR`` (the malignant background of a GBM section).

Why the gate (and not a bare arg-max). A bare arg-max INFLATES the single-marker
``ENDOTHELIAL`` class to ~20%: its score is one raw PECAM1 z-score (high variance)
while the multi-marker types are *means* pulled toward zero by averaging, so
endothelium wins any near-zero tie. Requiring a real positive signal (z >= 1.0)
breaks that single-vs-multi-marker imbalance and pulls the endothelial fraction
down to a biologically realistic ~10% (consistent with the Visium deconvolution's
~9.8%). The map (only markers actually present in the 319-gene panel are used):

  TUMOR       <- EGFR, SOX2, OLIG2, MKI67        (malignant glioma)
  IMMUNE      <- PTPRC, CD68, AIF1, CD163        (TAM / lymphoid)
  ENDOTHELIAL <- PECAM1                          (vessel-lining)
  STROMAL     <- MOG (oligodendrocyte),          (normal / stroma niche)
                 AQP4, GJA1 (astrocyte)

NECROTIC is a valid :class:`CellType` but is NOT assigned here: the panel
carries no necrosis marker, so no cell is phenotyped necrotic from expression
(reported honestly rather than invented).
"""
from __future__ import annotations

import os
from typing import Dict, List

import numpy as np
from scipy import sparse
from scipy.ndimage import gaussian_filter

from ..core.domain import Domain
from ..core.microenvironment import (
    CellPopulation,
    CellState,
    CellType,
    Microenvironment,
    TissueFields,
)
from ..data.real.singlecell_spatial import load as load_singlecell_spatial

# ---------------------------------------------------------------------------
# Marker -> CellType map (documented; only present markers are used at runtime).
# ---------------------------------------------------------------------------
MARKER_MAP: Dict[CellType, List[str]] = {
    CellType.TUMOR: ["EGFR", "SOX2", "OLIG2", "MKI67"],
    CellType.IMMUNE: ["PTPRC", "CD68", "AIF1", "CD163"],
    CellType.ENDOTHELIAL: ["PECAM1"],
    CellType.STROMAL: ["MOG", "AQP4", "GJA1"],
}

# Bandwidth (um) for the Gaussian density fields (ECM / vascular). ~one cell
# diameter scale of spatial smoothing so a field varies over tissue niches.
_DENSITY_BANDWIDTH_UM = 50.0

# Positive marker-signal gate (per-cell mean z-score) a type must clear to be
# assigned. z >= 1.0 == the marker signal is at least ~1 SD above the cohort
# mean (a genuine positive call). This is what breaks the single-vs-multi-marker
# imbalance that a bare arg-max suffers (see the module docstring).
TYPING_THRESHOLD = 1.0


def _zscore_columns(dense: np.ndarray) -> np.ndarray:
    """log1p then per-column (per-gene) z-score across cells. Returns float64."""
    x = np.log1p(dense.astype(np.float64))
    mean = x.mean(axis=0, keepdims=True)
    std = x.std(axis=0, keepdims=True)
    std[std == 0.0] = 1.0  # constant genes contribute 0, not NaN
    return (x - mean) / std


def phenotype_cells(
    expression, genes: np.ndarray, threshold: float = TYPING_THRESHOLD
) -> tuple[np.ndarray, Dict[CellType, List[str]], List[str]]:
    """Assign a CellType to each cell by a THRESHOLD-GATED marker z-score rule.

    Each type's score is the mean log1p-z-score of its present markers. A cell is
    assigned the type with the highest score AMONG the types whose score clears a
    positive ``threshold`` (a genuine marker-positive call); cells that clear no
    gate default to ``TUMOR`` (the malignant GBM background). The positive gate is
    what breaks the single-vs-multi-marker imbalance a bare arg-max produces (which
    inflates the single-marker ENDOTHELIAL class); see the module docstring.

    Parameters
    ----------
    expression : scipy.sparse or ndarray (n_cells, n_genes)
        Real per-cell targeted-panel counts.
    genes : ndarray of str
        Column labels for ``expression``.
    threshold : float
        Positive mean-z gate a type must clear to be assignable (default
        :data:`TYPING_THRESHOLD` = 1.0, i.e. ~1 SD above the cohort mean).

    Returns
    -------
    cell_type : ndarray (n_cells,) int8 of CellType values.
    used_map : dict CellType -> list of markers actually found in the panel.
    missing : list of markers requested but absent from the panel.
    """
    genes = np.asarray(genes, dtype=str)
    gene_index = {g: j for j, g in enumerate(genes)}

    used_map: Dict[CellType, List[str]] = {}
    missing: List[str] = []
    needed_cols: List[int] = []
    for ct, markers in MARKER_MAP.items():
        present = [m for m in markers if m in gene_index]
        absent = [m for m in markers if m not in gene_index]
        missing.extend(absent)
        used_map[ct] = present
        needed_cols.extend(gene_index[m] for m in present)

    # Densify only the marker columns (small; keeps memory light for 40k cells).
    needed_cols = sorted(set(needed_cols))
    if sparse.issparse(expression):
        sub = expression[:, needed_cols].toarray()
    else:
        sub = np.asarray(expression)[:, needed_cols]
    z = _zscore_columns(sub)  # (n_cells, len(needed_cols))
    col_pos = {c: k for k, c in enumerate(needed_cols)}

    n_cells = z.shape[0]
    types = list(used_map.keys())
    scores = np.full((n_cells, len(types)), -np.inf, dtype=np.float64)
    for t_idx, ct in enumerate(types):
        cols = [col_pos[gene_index[m]] for m in used_map[ct]]
        if cols:
            scores[:, t_idx] = z[:, cols].mean(axis=1)

    # Threshold gate: a type is only assignable to a cell when its mean z-score
    # clears the positive threshold. Non-clearing scores are set to -inf so they
    # can never win the arg-max.
    gated = np.where(scores >= float(threshold), scores, -np.inf)
    winner = gated.argmax(axis=1)
    cell_type = np.array([int(types[i]) for i in winner], dtype=np.int8)

    # Cells that clear NO gate (all types below threshold) default to the
    # malignant background: TUMOR.
    ungated = ~np.isfinite(gated.max(axis=1))
    cell_type[ungated] = int(CellType.TUMOR)
    return cell_type, used_map, sorted(set(missing))


def _provisional_sensitivity(
    expression, genes: np.ndarray, cell_type: np.ndarray
) -> np.ndarray:
    """PROVISIONAL per-cell sensitivity in [0,1] scaled from the tumour signal.

    NOT a measured drug response. Defined as the min-max normalised mean z-score
    of the TUMOR markers (cells with stronger malignant signal -> nearer 1).
    """
    genes = np.asarray(genes, dtype=str)
    gene_index = {g: j for j, g in enumerate(genes)}
    cols = [gene_index[m] for m in MARKER_MAP[CellType.TUMOR] if m in gene_index]
    if sparse.issparse(expression):
        sub = expression[:, cols].toarray()
    else:
        sub = np.asarray(expression)[:, cols]
    z = _zscore_columns(sub).mean(axis=1)
    lo, hi = np.percentile(z, 1), np.percentile(z, 99)
    if hi <= lo:
        return np.full(len(z), 0.5, dtype=np.float64)
    return np.clip((z - lo) / (hi - lo), 0.0, 1.0)


def _density_field(
    x_um: np.ndarray,
    y_um: np.ndarray,
    domain: Domain,
    bandwidth_um: float = _DENSITY_BANDWIDTH_UM,
) -> np.ndarray:
    """Gaussian kernel-density estimate of point (x,y) on the grid, in [0,1].

    Points are binned onto the (ny, nx) grid then convolved with a Gaussian of
    the given bandwidth (a discretised Gaussian KDE), and rescaled by its max.
    """
    grid = domain.zeros()
    if len(x_um) > 0:
        col, row = domain.to_index(x_um, y_um)
        np.add.at(grid, (row, col), 1.0)
    sigma = max(bandwidth_um / domain.dx, 1e-6)
    dens = gaussian_filter(grid, sigma=sigma, mode="constant")
    m = dens.max()
    if m > 0:
        dens = dens / m
    return dens


def reconstruct_from_xenium(
    sample_path: str | None = None, dx_um: float = 30.0
) -> Microenvironment:
    """Reconstruct the real cell-resolved GBM microenvironment from Xenium.

    Parameters
    ----------
    sample_path : str, optional
        Directory holding the cached Xenium files. ``None`` uses the loader
        default (``data_cache/xenium_gbm``).
    dx_um : float
        Grid spacing in micrometers for the continuum Domain / fields.

    Returns
    -------
    Microenvironment
        Real cells at real (x, y) positions with marker-derived CellTypes, a
        Domain covering the real tissue extent, and minimal continuum fields
        (vascular from endothelial density, ECM from stromal density; oxygen
        left None for the physics module to fill).
    """
    # --- load real data (retry once on transient failure) -------------------
    kwargs = {"dense": False}
    if sample_path is not None:
        kwargs["path"] = sample_path
    try:
        data = load_singlecell_spatial(**kwargs)
    except Exception:
        data = load_singlecell_spatial(**kwargs)

    xy = np.asarray(data["cell_xy_um"], dtype=np.float64)
    expression = data["expression"]
    genes = np.asarray(data["genes"], dtype=str)

    # --- phenotype every cell from real marker signal -----------------------
    cell_type, used_map, missing = phenotype_cells(expression, genes)
    sensitivity = _provisional_sensitivity(expression, genes, cell_type)

    # --- Domain covering the real extent (shift so min coord -> 0) ----------
    x = xy[:, 0] - xy[:, 0].min()
    y = xy[:, 1] - xy[:, 1].min()
    nx = int(np.ceil(x.max() / dx_um)) + 1
    ny = int(np.ceil(y.max() / dx_um)) + 1

    # tissue_mask = grid cells that contain at least one real cell.
    tissue_mask = np.zeros((ny, nx), dtype=bool)
    col = np.clip((x / dx_um).astype(int), 0, nx - 1)
    row = np.clip((y / dx_um).astype(int), 0, ny - 1)
    tissue_mask[row, col] = True

    domain = Domain(nx=nx, ny=ny, dx=dx_um, tissue_mask=tissue_mask)

    # --- CellPopulation: real positions + types, all VIABLE -----------------
    state = np.full(len(x), int(CellState.VIABLE), dtype=np.int8)
    cells = CellPopulation(
        x=x, y=y, cell_type=cell_type, state=state, sensitivity=sensitivity
    )

    # --- minimal continuum fields from real cell densities ------------------
    endo = cell_type == int(CellType.ENDOTHELIAL)
    strom = cell_type == int(CellType.STROMAL)
    vascular_density = _density_field(x[endo], y[endo], domain)
    ecm_density = _density_field(x[strom], y[strom], domain)
    fields = TissueFields(
        ecm_density=ecm_density,
        vascular_density=vascular_density,
        oxygen=None,  # filled later by the physics module
    )

    # --- type distribution for metadata / reporting -------------------------
    dist = {ct.name: int((cell_type == int(ct)).sum()) for ct in CellType}
    metadata = {
        "assay": data["source"]["dataset"],
        "panel": data["source"]["panel"],
        "n_panel_genes": int(len(genes)),
        "marker_map": {ct.name: used_map.get(ct, []) for ct in MARKER_MAP},
        "markers_missing_from_panel": missing,
        "celltype_distribution": dist,
        "sensitivity_note": (
            "PROVISIONAL: min-max scaled TUMOR marker z-score, NOT a measured "
            "drug response (this open section has no outcome label)."
        ),
        "fields_note": (
            "vascular_density = Gaussian KDE of ENDOTHELIAL cells; "
            "ecm_density = Gaussian KDE of STROMAL cells; oxygen left None."
        ),
        "has_outcome_label": bool(data["source"].get("has_outcome_label", False)),
    }

    return Microenvironment(
        domain=domain,
        fields=fields,
        cells=cells,
        patient_id="xenium-gbm-real",
        metadata=metadata,
    )


if __name__ == "__main__":
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    micro = reconstruct_from_xenium()

    print("=== REAL Xenium GBM microenvironment ===")
    print("summary:")
    for k, v in micro.summary().items():
        print(f"  {k}: {v}")

    print(f"\ngrid dims (ny, nx): {micro.domain.shape}  dx={micro.domain.dx} um")
    print(f"extent (um):        {micro.domain.extent_um}")

    print("\nmarker -> type map (present markers):")
    for ct_name, markers in micro.metadata["marker_map"].items():
        print(f"  {ct_name:<12}: {markers}")
    if micro.metadata["markers_missing_from_panel"]:
        print("markers requested but ABSENT from panel:",
              micro.metadata["markers_missing_from_panel"])

    ct = micro.cells.cell_type
    n = micro.cells.n
    print(f"\ncell-type distribution among {n:,} real cells:")
    for c in CellType:
        k = int((ct == int(c)).sum())
        print(f"  {c.name:<12}: {k:>7,}  ({100.0 * k / n:5.2f}%)")

    # --- scatter colored by CellType --------------------------------------
    colors = {
        CellType.TUMOR: "#d62728",
        CellType.STROMAL: "#2ca02c",
        CellType.IMMUNE: "#1f77b4",
        CellType.ENDOTHELIAL: "#ff7f0e",
        CellType.NECROTIC: "#7f7f7f",
    }
    fig, ax = plt.subplots(figsize=(10, 10))
    for c in CellType:
        m = ct == int(c)
        if m.any():
            ax.scatter(
                micro.cells.x[m], micro.cells.y[m], s=2, c=colors[c],
                label=f"{c.name} ({int(m.sum()):,})", linewidths=0,
            )
    ax.set_aspect("equal")
    ax.invert_yaxis()  # image convention
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.set_title("Real Xenium GBM microenvironment — cells by phenotype")
    ax.legend(markerscale=4, loc="upper right", fontsize=8)

    out_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "..", "results"
    )
    out_dir = os.path.abspath(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "real_xenium_microenvironment.png")
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    print(f"\nfigure written: {out_path}")
