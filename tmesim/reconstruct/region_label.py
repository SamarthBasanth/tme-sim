"""Niche-region labeling from REAL Ivy GAP per-region expression signatures.

Assign each region of a reconstructed microenvironment (a Visium spot or an
aggregated Xenium neighborhood) to a genuine GBM **anatomic zone** —

    margin        (leading edge / infiltrating tumour: the tumour/normal-brain
                   interface)
    core          (cellular tumour: dense proliferative bulk)
    perivascular  (microvascular proliferation / hyperplastic blood vessels)
    necrotic      (pseudopalisading + perinecrotic / hypoxic zone)

by scoring the location's local expression against the **real** Ivy GAP
per-anatomic-structure RNA-seq signatures
(:func:`tmesim.data.real.ivygap.structure_signatures`, grouped by ``region``).
No synthetic data enters this module: the signatures are measured FPKM from
laser-microdissected GBM structures (Puchalski et al. 2018, *Science*), and the
local expression to classify is measured Visium / Xenium counts.

Pipeline
--------
1. :func:`region_signatures` builds a real per-region mean-FPKM signature over a
   region-discriminative gene panel (data-driven top-variance genes across the
   four zones, unioned with canonical GBM markers).
2. :func:`label_regions` scores each location against every region signature by
   rank (Spearman) correlation over the genes SHARED between the panel and the
   local assay, and assigns the best-correlating region.
3. :func:`region_field` rasterizes the per-location labels into smooth
   per-region fraction fields on a :class:`~tmesim.core.Domain` grid, for the
   whole-niche assembler and for plotting.

The scoring is cross-platform (FPKM signature vs. UMI/targeted-count local
expression): both are log1p-transformed and rank-correlated, which is invariant
to the monotone scale difference between the two assays.

Only numpy / scipy / matplotlib / stdlib are used. Core is not modified.
"""
from __future__ import annotations

import time
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
from scipy.stats import rankdata

from ..core.domain import Domain
from ..data.real import ivygap

# Canonical, literature-grounded GBM zone markers. Always folded into the
# signature panel (when present in Ivy GAP) so the region biology is anchored,
# and used by the honest marker-separation validation in ``__main__``.
REGION_MARKERS: Dict[str, List[str]] = {
    # Leading-edge / infiltrating margin looks the most like normal brain:
    # oligodendrocyte + neuronal / synaptic programs.
    "margin": ["MBP", "PLP1", "MOBP", "MAG", "SNAP25", "SYT1", "GABRA1", "RBFOX3"],
    # Cellular-tumour core: classic GBM proliferative / progenitor program.
    "core": ["EGFR", "SOX2", "OLIG2", "MKI67", "CDK4", "PDGFRA", "DLL3"],
    # Perivascular / microvascular proliferation: endothelial + angiogenic.
    "perivascular": ["PECAM1", "VWF", "CD34", "CLDN5", "ANGPT2", "ESM1", "FLT1"],
    # Necrotic / perinecrotic hypoxic zone: HIF targets.
    "necrotic": ["CA9", "VEGFA", "ADM", "NDRG1", "SLC2A1", "HILPDA", "IGFBP3"],
}

# Region column order kept everywhere for stable plots / fields.
REGION_ORDER: List[str] = ["margin", "core", "perivascular", "necrotic"]

_UNLABELED = "unlabeled"


# ---------------------------------------------------------------------------
# small utilities
# ---------------------------------------------------------------------------
def _retry(fn, *args, **kwargs):
    """Call ``fn`` and retry ONCE on a transient/OS/network-style error.

    The Ivy GAP / Visium loaders read a local cache, so this only ever fires on
    a genuine I/O hiccup; a ``FileNotFoundError`` (missing cache) is NOT
    retried — it is a hard configuration error.
    """
    try:
        return fn(*args, **kwargs)
    except FileNotFoundError:
        raise
    except (OSError, ConnectionError, TimeoutError):
        time.sleep(1.0)
        return fn(*args, **kwargs)


def _densify(expr) -> np.ndarray:
    """Return a dense float64 (n_loc, n_genes) view of a dense or sparse matrix."""
    if hasattr(expr, "toarray"):
        return expr.toarray().astype(np.float64)
    return np.asarray(expr, dtype=np.float64)


# ---------------------------------------------------------------------------
# 1. real per-region signatures
# ---------------------------------------------------------------------------
def region_signatures(
    gene_set: Optional[Sequence[str]] = None,
    n_genes: int = 300,
    min_fpkm: float = 5.0,
    include_markers: bool = True,
    reference_histology_only: bool = False,
) -> Tuple[List[str], np.ndarray, List[str]]:
    """Build the REAL Ivy GAP per-region mean-expression signature matrix.

    Groups the Ivy GAP anatomic-structure FPKM by microenvironment ``region``
    (margin / core / perivascular / necrotic) and reduces it to a
    region-discriminative gene panel.

    Parameters
    ----------
    gene_set : sequence of str, optional
        Explicit gene panel to use. If ``None`` (default) the panel is chosen
        data-drivenly: the ``n_genes`` genes with the largest spread across the
        four region means (in log space), restricted to genes expressed at
        ``>= min_fpkm`` in at least one region, unioned with the canonical
        :data:`REGION_MARKERS` (when ``include_markers``).
    n_genes : int
        Size of the data-driven variable-gene panel (ignored if ``gene_set``).
    min_fpkm : float
        Minimum peak-across-regions FPKM for a gene to be eligible (noise floor).
    include_markers : bool
        Fold the canonical :data:`REGION_MARKERS` into the panel.
    reference_histology_only : bool
        Passed through to Ivy GAP. ``False`` (default) so the CSC-cluster niches
        populate the perivascular (CThbv) and necrotic (CTpnz) region means.

    Returns
    -------
    region_names : list[str]
        The region labels, in :data:`REGION_ORDER` (rows of ``signature_matrix``).
    signature_matrix : np.ndarray, shape (n_regions, n_panel_genes)
        Mean FPKM per region over the panel (raw FPKM; the scorer log-transforms).
    gene_names : list[str]
        Panel gene symbols, aligned to the columns of ``signature_matrix``.
    """
    sig = _retry(
        ivygap.structure_signatures,
        by="region",
        reference_histology_only=reference_histology_only,
    )  # genes x regions DataFrame
    # Keep the canonical region order, dropping any region absent from the data.
    regions = [r for r in REGION_ORDER if r in sig.columns]
    sig = sig[regions]

    if gene_set is not None:
        panel = [g for g in gene_set if g in sig.index]
    else:
        expressed = sig[sig.max(axis=1) >= float(min_fpkm)]
        # Spread across regions in log space = how region-discriminative a gene is.
        log_e = np.log1p(expressed.values)
        spread = log_e.max(axis=1) - log_e.min(axis=1)
        order = np.argsort(spread)[::-1]
        top = list(expressed.index[order[: int(n_genes)]])
        panel = list(top)
        if include_markers:
            have = set(panel)
            for genes in REGION_MARKERS.values():
                for g in genes:
                    if g in sig.index and g not in have:
                        panel.append(g)
                        have.add(g)

    if not panel:
        raise ValueError("region_signatures: empty gene panel (check gene_set / cache)")

    sub = sig.loc[panel]
    signature_matrix = sub.values.T.astype(np.float64)  # regions x genes
    return regions, signature_matrix, list(sub.index)


# ---------------------------------------------------------------------------
# 2. label each location by correlation to the region signatures
# ---------------------------------------------------------------------------
def _spearman_matrix(X: np.ndarray, S: np.ndarray) -> np.ndarray:
    """Row-wise Spearman correlation between every row of X and every row of S.

    X : (n_loc, g), S : (n_regions, g). Returns (n_loc, n_regions).
    Rank each row across genes, then Pearson-correlate the ranks. Rows of X with
    no variance (e.g. all-zero over the shared panel) yield NaN and are handled
    by the caller.
    """
    Xr = rankdata(X, axis=1).astype(np.float64)
    Sr = rankdata(S, axis=1).astype(np.float64)
    Xc = Xr - Xr.mean(axis=1, keepdims=True)
    Sc = Sr - Sr.mean(axis=1, keepdims=True)
    Xn = np.linalg.norm(Xc, axis=1, keepdims=True)
    Sn = np.linalg.norm(Sc, axis=1, keepdims=True)
    # errstate spans the matmul too: some BLAS backends raise a spurious
    # divide/overflow FP flag on gemm even though the (finite) result is fine.
    with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
        corr = (Xc @ Sc.T) / (Xn * Sn.T)
    return corr


def label_regions(
    local_expr,
    local_coords: Optional[np.ndarray],
    gene_names: Sequence[str],
    signatures: Optional[Tuple[List[str], np.ndarray, List[str]]] = None,
    return_scores: bool = False,
) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray, List[str], List[str]]]:
    """Assign each location to its best-correlating Ivy GAP anatomic region.

    Parameters
    ----------
    local_expr : ndarray or scipy sparse, shape (n_loc, n_genes)
        Local measured expression — Visium spot counts or aggregated Xenium
        neighborhood counts — with locations as rows and genes as columns.
    local_coords : ndarray (n_loc, 2), optional
        Location coordinates in micrometers. Not used for scoring (each location
        is scored independently); accepted for API symmetry and downstream use.
    gene_names : sequence of str
        Gene symbols aligned to the columns of ``local_expr``.
    signatures : (region_names, signature_matrix, gene_names), optional
        Precomputed :func:`region_signatures` output; computed if ``None``.
    return_scores : bool
        If True also return the per-location correlation matrix and names.

    Returns
    -------
    region_per_location : ndarray (n_loc,) of str
        The assigned region for each location (one of :data:`REGION_ORDER`, or
        ``"unlabeled"`` for a location with no signal over the shared panel).
    (if ``return_scores``) scores : ndarray (n_loc, n_regions), region_names :
        list[str], shared_genes : list[str].
    """
    if signatures is None:
        signatures = region_signatures()
    region_names, sig_mat, sig_genes = signatures

    # Shared genes between the local assay and the signature panel.
    col_of = {g: j for j, g in enumerate(gene_names)}
    shared = [g for g in sig_genes if g in col_of]
    if len(shared) < 3:
        raise ValueError(
            f"only {len(shared)} shared genes between local assay and Ivy GAP "
            "signature panel — too few to correlate. Check gene symbols."
        )
    local_cols = np.array([col_of[g] for g in shared], dtype=int)
    sig_cols = np.array([sig_genes.index(g) for g in shared], dtype=int)

    X = _densify(local_expr)[:, local_cols]  # (n_loc, n_shared)
    S = sig_mat[:, sig_cols]                 # (n_regions, n_shared)

    # Log1p both sides (Spearman is rank-invariant to this, but it keeps the
    # transform explicit and matches the Pearson fallback if ever swapped in).
    X = np.log1p(np.clip(X, 0.0, None))
    S = np.log1p(np.clip(S, 0.0, None))

    scores = _spearman_matrix(X, S)  # (n_loc, n_regions)

    labels = np.full(scores.shape[0], _UNLABELED, dtype=object)
    finite = np.isfinite(scores).all(axis=1)
    best = np.argmax(np.where(np.isfinite(scores), scores, -np.inf), axis=1)
    reg_arr = np.array(region_names, dtype=object)
    labels[finite] = reg_arr[best[finite]]
    labels = labels.astype(str)

    if return_scores:
        return labels, scores, list(region_names), shared
    return labels


# ---------------------------------------------------------------------------
# 3. rasterize labels into per-region fraction fields on the grid
# ---------------------------------------------------------------------------
def region_field(
    domain: Domain,
    coords_um: np.ndarray,
    region_labels: Sequence[str],
    bandwidth_um: float = 100.0,
) -> Dict[str, np.ndarray]:
    """Rasterize per-location region labels into smooth per-region fraction fields.

    For each region, bins the locations assigned to it onto the grid, smooths
    with a Gaussian (bandwidth ~ the spot pitch), and divides by the smoothed
    total of all LABELED locations so every returned field is a spatial fraction
    in ``[0, 1]``; across regions the fields sum to ~1 where any labeled location
    has support and to 0 elsewhere.

    Parameters
    ----------
    domain : Domain
        Target grid; fields are returned with shape ``domain.shape`` (ny, nx).
    coords_um : ndarray (n_loc, 2)
        Location coordinates in micrometers, aligned to ``region_labels``.
    region_labels : sequence of str
        Region per location from :func:`label_regions`.
    bandwidth_um : float
        Gaussian smoothing bandwidth in micrometers (default 100 um ~ the Visium
        centre-to-centre spot pitch).

    Returns
    -------
    dict[str, np.ndarray]
        Region name -> (ny, nx) fraction field in [0, 1]. Keys cover every region
        in :data:`REGION_ORDER` that appears in ``region_labels``.
    """
    from scipy.ndimage import gaussian_filter

    coords_um = np.asarray(coords_um, dtype=np.float64)
    labels = np.asarray(region_labels, dtype=object).astype(str)
    if coords_um.shape[0] != labels.shape[0]:
        raise ValueError(
            f"coords_um ({coords_um.shape[0]}) and region_labels "
            f"({labels.shape[0]}) length mismatch"
        )

    sigma = float(bandwidth_um) / float(domain.dx)
    col, row = domain.to_index(coords_um[:, 0], coords_um[:, 1])

    present = [r for r in REGION_ORDER if r in set(labels)]
    # Smoothed count of all labeled locations = the normalizing denominator.
    labeled = labels != _UNLABELED
    total = domain.zeros()
    if labeled.any():
        np.add.at(total, (row[labeled], col[labeled]), 1.0)
    total_s = gaussian_filter(total, sigma=sigma)

    eps = 1e-12
    fields: Dict[str, np.ndarray] = {}
    for reg in present:
        m = labels == reg
        g = domain.zeros()
        if m.any():
            np.add.at(g, (row[m], col[m]), 1.0)
        g_s = gaussian_filter(g, sigma=sigma)
        frac = np.where(total_s > eps, g_s / (total_s + eps), 0.0)
        frac = np.where(domain.tissue_mask, frac, 0.0)
        fields[reg] = np.clip(frac, 0.0, 1.0)
    return fields


# ---------------------------------------------------------------------------
# demo / self-test on REAL data
# ---------------------------------------------------------------------------
def _domain_from_coords(coords_um: np.ndarray, dx: float = 50.0) -> Domain:
    """A grid that covers the coordinate cloud with a small margin."""
    w = float(coords_um[:, 0].max())
    h = float(coords_um[:, 1].max())
    nx = int(np.ceil((w + dx) / dx))
    ny = int(np.ceil((h + dx) / dx))
    return Domain(nx=nx, ny=ny, dx=dx)


def _cp10k_log(expr, genes: Sequence[str], want: Sequence[str]) -> Dict[str, np.ndarray]:
    """Per-spot counts-per-10k, log1p, for a few marker genes (for validation)."""
    X = _densify(expr)
    tot = X.sum(axis=1, keepdims=True)
    tot[tot == 0] = 1.0
    norm = np.log1p(X / tot * 1e4)
    col = {g: j for j, g in enumerate(genes)}
    return {g: norm[:, col[g]] for g in want if g in col}


if __name__ == "__main__":
    import os

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from ..data.real import visium

    SAMPLE = "10x_parent_gbm"
    print(f"== Niche-region labeling on REAL Visium sample {SAMPLE!r} ==\n")

    # --- real signatures ----------------------------------------------------
    region_names, sig_mat, panel = region_signatures()
    print(f"Ivy GAP region signature: {len(region_names)} regions "
          f"({', '.join(region_names)}) x {len(panel)}-gene panel")

    # --- real local expression ----------------------------------------------
    d = _retry(visium.load_sample, SAMPLE)
    genes = d["genes"]
    coords = d["coords_um"]
    print(f"Visium: {d['n_spots']} spots x {d['n_genes']} genes "
          f"({d['accession']})")

    labels, scores, rnames, shared = label_regions(
        d["expression"], coords, genes,
        signatures=(region_names, sig_mat, panel),
        return_scores=True,
    )
    print(f"shared genes (panel ∩ Visium): {len(shared)} / {len(panel)}")

    # --- region distribution -------------------------------------------------
    print("\nRegion label distribution:")
    uniq, counts = np.unique(labels, return_counts=True)
    dist = dict(zip(uniq.tolist(), counts.tolist()))
    for r in REGION_ORDER + [_UNLABELED]:
        if r in dist:
            print(f"  {r:12s} {dist[r]:5d}  ({100*dist[r]/len(labels):5.1f} %)")

    # --- marker validation ---------------------------------------------------
    print("\nMarker validation (mean per-spot log CP10k by ASSIGNED region):")
    val_genes = ["PECAM1", "CA9", "VEGFA", "MKI67", "MBP"]
    mk = _cp10k_log(d["expression"], genes, val_genes)
    header = "  region       " + "".join(f"{g:>9s}" for g in mk)
    print(header)
    for r in REGION_ORDER:
        m = labels == r
        if not m.any():
            continue
        row = "".join(f"{mk[g][m].mean():9.3f}" for g in mk)
        print(f"  {r:12s}" + row)
    # explicit checks the task asked for
    peri = labels == "perivascular"
    nec = labels == "necrotic"
    if "PECAM1" in mk and peri.any():
        rest = ~peri
        print(f"\n  PECAM1 in perivascular={mk['PECAM1'][peri].mean():.3f} "
              f"vs rest={mk['PECAM1'][rest].mean():.3f}  "
              f"-> {'HIGHER (validated)' if mk['PECAM1'][peri].mean() > mk['PECAM1'][rest].mean() else 'NOT higher'}")
    for g in ("CA9", "VEGFA"):
        if g in mk and nec.any():
            rest = ~nec
            print(f"  {g:6s} in necrotic    ={mk[g][nec].mean():.3f} "
                  f"vs rest={mk[g][rest].mean():.3f}  "
                  f"-> {'HIGHER (validated)' if mk[g][nec].mean() > mk[g][rest].mean() else 'NOT higher'}")

    # --- rasterize + save a plot --------------------------------------------
    dom = _domain_from_coords(coords, dx=50.0)
    fields = region_field(dom, coords, labels)
    print(f"\nregion_field on {dom.nx}x{dom.ny} grid (dx={dom.dx} um): "
          f"{list(fields.keys())}")

    os.makedirs("results", exist_ok=True)
    colors = {
        "margin": "#4c9f70", "core": "#d1495b",
        "perivascular": "#3b6ea5", "necrotic": "#8a5a2b",
        _UNLABELED: "#bbbbbb",
    }
    fig, axes = plt.subplots(1, 2, figsize=(15, 6.5))
    ax = axes[0]
    for r in REGION_ORDER + [_UNLABELED]:
        m = labels == r
        if m.any():
            ax.scatter(coords[m, 0], coords[m, 1], s=9, c=colors[r],
                       label=f"{r} ({m.sum()})", edgecolors="none")
    ax.set_aspect("equal"); ax.invert_yaxis()
    ax.set_title(f"Ivy GAP region label per spot — {SAMPLE}")
    ax.set_xlabel("x (um)"); ax.set_ylabel("y (um)")
    ax.legend(markerscale=2, fontsize=8, loc="upper right")

    ax = axes[1]
    peri_f = fields.get("perivascular")
    nec_f = fields.get("necrotic")
    ext = [0, dom.extent_um[0], dom.extent_um[1], 0]
    if nec_f is not None:
        ax.imshow(nec_f, cmap="Oranges", alpha=0.85, extent=ext, origin="upper")
    if peri_f is not None:
        ax.imshow(peri_f, cmap="Blues", alpha=0.5, extent=ext, origin="upper")
    ax.set_aspect("equal")
    ax.set_title("Perivascular (blue) & necrotic (orange) fraction fields")
    ax.set_xlabel("x (um)"); ax.set_ylabel("y (um)")

    fig.tight_layout()
    out = os.path.join("results", "real_regions.png")
    fig.savefig(out, dpi=130)
    print(f"\nsaved {out}")
    print("\nregion_label.py ran on REAL data OK")
