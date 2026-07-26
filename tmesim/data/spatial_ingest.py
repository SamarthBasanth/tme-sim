"""Ingest a REAL spatial-transcriptomics glioblastoma sample into a tmesim
:class:`~tmesim.core.microenvironment.Microenvironment`.

This closes the loop the rest of the package only *approximated*: instead of a
synthetic 2-D layout parameterised by bulk features (see
``tmesim/data/synthetic.py`` and ``docs/data_acquisition_plan.md`` section B),
we place agents at the **measured spatial coordinates** of a 10x Visium spot
array and phenotype each spot from its **measured marker-gene expression**, then
hand the result to the *unmodified* ``reconstruct`` package to build the ECM /
vascular continuum fields the physics + sim stages consume.

What is real here vs. approximated
----------------------------------
REAL (measured bytes on disk):
  * spot (x, y) coordinates from ``spatial/tissue_positions_list.csv`` converted
    to micrometers via the Visium ``scalefactors_json.json``;
  * per-spot expression of a GBM marker panel from the MEX sparse matrix
    (``matrix.mtx.gz`` + ``features.tsv.gz`` + ``barcodes.tsv.gz``).

APPROXIMATED (documented, honest):
  * one Visium spot is a ~55 um well pooling ~1-10 cells, NOT a single cell — we
    model each spot as ONE representative agent at the spot centroid;
  * phenotype = argmax over a small hand-curated marker panel (no reference-based
    deconvolution);
  * per-agent ``sensitivity`` is a tumor/proliferation-signal heuristic (mirrors
    ``reconstruct.segmentation._sensitivity_from_tumor``) because open spatial
    data carries **no drug-response label** (the crux — see
    ``docs/real_data_status.md``);
  * ECM + vasculature come from ``reconstruct`` estimators, not a measured vessel
    stain.

The dataset (open, directly downloadable, no login):
  10x Genomics "Human Glioblastoma, Whole Transcriptome Analysis" (Visium v1,
  Space Ranger 1.2.0). Files under
  ``cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_Glioblastoma/``.
"""
from __future__ import annotations

import gzip
import json
import os
import tarfile
from typing import Dict, List, Optional, Tuple, Union

import numpy as np

# ---------------------------------------------------------------------------
# Dataset provenance (verified reachable 2026-07-26; all HTTP 200)
# ---------------------------------------------------------------------------
VISIUM_GBM_BASE = (
    "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/"
    "Parent_Visium_Human_Glioblastoma"
)
VISIUM_GBM_FILES = {
    # name -> (url suffix, approx content-length bytes verified via HEAD)
    "spatial": (
        f"{VISIUM_GBM_BASE}/Parent_Visium_Human_Glioblastoma_spatial.tar.gz",
        9_134_896,
    ),
    "filtered_matrix": (
        f"{VISIUM_GBM_BASE}/"
        "Parent_Visium_Human_Glioblastoma_filtered_feature_bc_matrix.tar.gz",
        54_134_945,
    ),
}
VISIUM_GBM_WEB_SUMMARY = (
    f"{VISIUM_GBM_BASE}/Parent_Visium_Human_Glioblastoma_web_summary.html"
)
VISIUM_GBM_LANDING = (
    "https://www.10xgenomics.com/datasets/"
    "human-glioblastoma-whole-transcriptome-analysis-1-standard-1-2-0"
)

# Visium v1 geometry: 55 um spot diameter. um-per-pixel = 55 / spot_diameter_px.
_SPOT_DIAMETER_UM = 55.0

# GBM marker panel: class -> gene symbols. Deliberately small and defensible.
#   TUMOR       : glioma stemness / classical driver / proliferation
#   STROMAL     : reactive astrocyte (GFAP)
#   IMMUNE      : TAM / microglia (CD68, AIF1=IBA1) + T cell (CD3E)
#   ENDOTHELIAL : vessel lining (PECAM1=CD31) + angiogenesis (VEGFA)
MARKER_PANEL: Dict[str, List[str]] = {
    "tumor": ["SOX2", "OLIG2", "EGFR", "MKI67"],
    "stroma": ["GFAP"],
    "immune": ["CD68", "AIF1", "CD3E"],
    "endothelial": ["PECAM1", "VEGFA"],
}
# Panel class -> CellType, resolved through reconstruct.segmentation aliases.
_CLASS_ORDER = ("tumor", "stroma", "immune", "endothelial")


def _log(msg: str) -> None:
    print(f"[spatial_ingest] {msg}", flush=True)


# ===========================================================================
# 1. Locate / fetch the open dataset
# ===========================================================================
def find_open_spatial_gbm(
    cache_dir: str = "data_cache/visium_gbm",
    timeout: int = 180,
    download: bool = True,
) -> Dict:
    """Document and (best-effort) fetch the best OPEN spatial GBM dataset.

    Returns a status dict describing the 10x Visium Human Glioblastoma sample,
    its verified direct-download URLs, whether the bytes are present locally, and
    — critically — that it carries **no paired drug-response outcome label**.

    Never raises: a network failure degrades to ``status='documented_only'`` with
    the URLs still returned, so the caller/reader can fetch manually.

    Parameters
    ----------
    cache_dir : where the tarballs are cached / extracted.
    timeout : per-request timeout (s) for the download attempt.
    download : if False, only report what is already cached (no network).
    """
    os.makedirs(cache_dir, exist_ok=True)
    status: Dict = {
        "dataset": "10x Genomics Human Glioblastoma (Visium v1, Space Ranger 1.2.0)",
        "assay": "10x Visium spatial gene expression (whole transcriptome)",
        "organism": "Homo sapiens, glioblastoma (GBM) surgical section",
        "landing_page": VISIUM_GBM_LANDING,
        "web_summary": VISIUM_GBM_WEB_SUMMARY,
        "urls": {k: v[0] for k, v in VISIUM_GBM_FILES.items()},
        "expected_bytes": {k: v[1] for k, v in VISIUM_GBM_FILES.items()},
        "has_outcome_label": False,
        "outcome_note": (
            "OPEN spatial GBM data has NO paired drug-response / RECIST outcome. "
            "This sample is a single untreated resection; it supports geometry "
            "reconstruction + a mechanistic prediction, NOT supervised validation."
        ),
        "cached": {},
        "status": "documented_only",
    }

    def _cached_ok(path: str, min_bytes: int) -> bool:
        return os.path.exists(path) and os.path.getsize(path) >= min_bytes

    for key, (url, exp) in VISIUM_GBM_FILES.items():
        dest = os.path.join(cache_dir, os.path.basename(url))
        # 0.9*exp tolerance so a slightly different mirror size still counts.
        if _cached_ok(dest, int(0.9 * exp)):
            status["cached"][key] = {"path": dest, "bytes": os.path.getsize(dest)}
            continue
        if not download:
            continue
        ok = _try_download(url, dest, timeout)
        if ok and _cached_ok(dest, int(0.9 * exp)):
            status["cached"][key] = {"path": dest, "bytes": os.path.getsize(dest)}

    n_cached = len(status["cached"])
    if n_cached == len(VISIUM_GBM_FILES):
        status["status"] = "downloaded"
    elif n_cached > 0:
        status["status"] = "partial"
    status["cache_dir"] = cache_dir
    return status


def _try_download(url: str, dest: str, timeout: int) -> bool:
    """Stream ``url`` to ``dest`` with one retry. Returns True on success. Never
    raises (network is allowed to be flaky)."""
    try:
        import requests  # local import so the module loads without requests
    except Exception as exc:  # noqa: BLE001
        _log(f"requests unavailable ({exc}); cannot download {url}")
        return False
    for attempt in range(2):
        try:
            with requests.get(url, timeout=timeout, stream=True) as r:
                r.raise_for_status()
                tmp = dest + ".part"
                with open(tmp, "wb") as fh:
                    for chunk in r.iter_content(chunk_size=1 << 16):
                        if chunk:
                            fh.write(chunk)
                os.replace(tmp, dest)
            _log(f"downloaded {os.path.basename(dest)} ({os.path.getsize(dest)} B)")
            return True
        except Exception as exc:  # noqa: BLE001
            _log(f"download attempt {attempt + 1} failed for {url}: {exc}")
    return False


# ===========================================================================
# 2. Read the real Visium bytes (coordinates + expression)
# ===========================================================================
def _ensure_extracted(cache_dir: str) -> Tuple[str, str]:
    """Ensure ``spatial/`` and ``filtered_feature_bc_matrix/`` are extracted in
    ``cache_dir``. Extracts the cached tarballs if needed. Returns the two dir
    paths. Raises FileNotFoundError if neither cache nor extracted dir exists."""
    spatial_dir = os.path.join(cache_dir, "spatial")
    matrix_dir = os.path.join(cache_dir, "filtered_feature_bc_matrix")

    for sub, tar_key in ((spatial_dir, "spatial"), (matrix_dir, "filtered_matrix")):
        if os.path.isdir(sub):
            continue
        tar_name = os.path.basename(VISIUM_GBM_FILES[tar_key][0])
        tar_path = os.path.join(cache_dir, tar_name)
        if not os.path.exists(tar_path):
            raise FileNotFoundError(
                f"neither {sub} nor {tar_path} present — run "
                f"find_open_spatial_gbm() first to fetch the dataset."
            )
        _log(f"extracting {tar_name} ...")
        with tarfile.open(tar_path, "r:gz") as tf:
            tf.extractall(cache_dir)
    return spatial_dir, matrix_dir


def _read_positions(spatial_dir: str) -> "Dict[str, Tuple[float, float]]":
    """Read ``tissue_positions_list.csv`` -> {barcode: (pxl_col_x, pxl_row_y)}
    for in-tissue spots only. Columns (no header) are:
    barcode, in_tissue, array_row, array_col, pxl_row_in_fullres, pxl_col_in_fullres.
    """
    path = os.path.join(spatial_dir, "tissue_positions_list.csv")
    coords: Dict[str, Tuple[float, float]] = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split(",")
            if len(parts) < 6:
                continue
            barcode, in_tissue = parts[0], parts[1]
            if in_tissue != "1":
                continue
            pxl_row = float(parts[4])  # y
            pxl_col = float(parts[5])  # x
            coords[barcode] = (pxl_col, pxl_row)
    return coords


def _read_scalefactor_umpp(spatial_dir: str) -> float:
    """um-per-pixel from ``scalefactors_json.json`` (55 um / spot_diameter_px)."""
    path = os.path.join(spatial_dir, "scalefactors_json.json")
    with open(path) as fh:
        sf = json.load(fh)
    diam = float(sf["spot_diameter_fullres"])
    return _SPOT_DIAMETER_UM / diam


def _read_marker_expression(
    matrix_dir: str, panel_genes: List[str]
) -> "Tuple[List[str], np.ndarray, np.ndarray, List[str]]":
    """Read the MEX sparse matrix and return
    ``(barcodes, marker_expr, total_counts, present_genes)`` where

      barcodes      : list[str], length n_spots (matrix column order)
      marker_expr   : (n_markers, n_spots) dense counts for the panel genes found
      total_counts  : (n_spots,) total UMI per spot (for CPM normalization)
      present_genes : the subset of ``panel_genes`` actually in the matrix, in
                      the same row order as ``marker_expr``

    Only the panel-gene rows are densified, so memory stays tiny.
    """
    from scipy.io import mmread

    # features.tsv.gz: ensembl_id \t symbol \t type  -> row index of each symbol.
    feat_path = os.path.join(matrix_dir, "features.tsv.gz")
    symbol_to_row: Dict[str, int] = {}
    with gzip.open(feat_path, "rt") as fh:
        for i, line in enumerate(fh):
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 2:
                symbol_to_row.setdefault(cols[1], i)  # first occurrence wins

    bc_path = os.path.join(matrix_dir, "barcodes.tsv.gz")
    with gzip.open(bc_path, "rt") as fh:
        barcodes = [ln.strip() for ln in fh if ln.strip()]

    mtx_path = os.path.join(matrix_dir, "matrix.mtx.gz")
    with gzip.open(mtx_path, "rb") as fh:
        mat = mmread(fh).tocsr()  # genes x spots

    total_counts = np.asarray(mat.sum(axis=0)).ravel().astype(np.float64)

    present = [g for g in panel_genes if g in symbol_to_row]
    if present:
        rows = [symbol_to_row[g] for g in present]
        marker_expr = mat[rows, :].toarray().astype(np.float64)  # (n_markers, n_spots)
    else:
        marker_expr = np.zeros((0, len(barcodes)), dtype=np.float64)
    return barcodes, marker_expr, total_counts, present  # type: ignore[return-value]


# ===========================================================================
# 3. Phenotype spots from marker genes -> CellPopulation
# ===========================================================================
def _minmax(a: np.ndarray) -> np.ndarray:
    lo, hi = float(a.min()), float(a.max())
    if hi - lo <= 0.0:
        return np.zeros_like(a)
    return (a - lo) / (hi - lo)


def _phenotype_spots(
    present_genes: List[str],
    marker_expr: np.ndarray,
    total_counts: np.ndarray,
) -> "Tuple[np.ndarray, np.ndarray, Dict[str, float]]":
    """CPM+log1p normalize marker expression, score each panel class per spot, and
    assign each spot a CellType (argmax class) + a tumor-signal-derived sensitivity.

    Returns ``(cell_type, sensitivity, class_means)``.
    """
    from tmesim.core import CellType

    n_spots = total_counts.shape[0]
    safe_tot = np.where(total_counts > 0, total_counts, 1.0)
    # CPM-like (per-10k) + log1p, then per-gene min-max so genes with different
    # dynamic range compete fairly (mirrors segmentation._phenotype's per-channel
    # normalization, but on genes instead of image channels).
    norm = np.empty_like(marker_expr)
    for j in range(marker_expr.shape[0]):
        cpm = marker_expr[j] / safe_tot * 1e4
        norm[j] = _minmax(np.log1p(cpm))

    gene_index = {g: i for i, g in enumerate(present_genes)}
    class_score = np.zeros((len(_CLASS_ORDER), n_spots), dtype=np.float64)
    class_means: Dict[str, float] = {}
    for ci, cls in enumerate(_CLASS_ORDER):
        genes = [g for g in MARKER_PANEL[cls] if g in gene_index]
        if genes:
            rows = [gene_index[g] for g in genes]
            class_score[ci] = norm[rows, :].mean(axis=0)
        class_means[cls] = float(class_score[ci].mean())

    # Phenotype = the class each spot is most RELATIVELY enriched for. We z-score
    # each class signal across spots before argmax, otherwise a single ubiquitously
    # high marker dominates: in GBM, GFAP (astrocyte) is broadly expressed — both
    # reactive astrocytes AND astrocytic tumor cells make it — so a raw argmax
    # collapses ~all spots to "stroma". Relative enrichment is a standard coarse
    # dominant-module assignment. It is NOT reference-based deconvolution (which
    # would need scanpy/cell2location, unavailable here) — see docs/real_data_status.md.
    z = (class_score - class_score.mean(axis=1, keepdims=True)) / (
        class_score.std(axis=1, keepdims=True) + 1e-9
    )
    class_to_celltype = {
        "tumor": int(CellType.TUMOR),
        "stroma": int(CellType.STROMAL),
        "immune": int(CellType.IMMUNE),
        "endothelial": int(CellType.ENDOTHELIAL),
    }
    winner = np.argmax(z, axis=0)
    cell_type = np.array(
        [class_to_celltype[_CLASS_ORDER[w]] for w in winner], dtype=np.int8
    )

    # Sensitivity from the ABSOLUTE tumor-class signal (proliferative,
    # drug-targetable), mirroring reconstruct.segmentation._sensitivity_from_tumor:
    #   tumor cells   -> 0.3 + 0.6*sig ;  non-tumor -> 0.1 + 0.2*sig.
    tumor_sig = np.clip(class_score[_CLASS_ORDER.index("tumor")], 0.0, 1.0)
    is_tumor = cell_type == int(CellType.TUMOR)
    sensitivity = np.where(is_tumor, 0.3 + 0.6 * tumor_sig, 0.1 + 0.2 * tumor_sig)
    sensitivity = np.clip(sensitivity, 0.0, 1.0)
    return cell_type, sensitivity, class_means


# ===========================================================================
# 4. Assemble the Microenvironment via the (unmodified) reconstruct package
# ===========================================================================
def microenvironment_from_spatial(
    path_or_data: Union[str, Dict],
    target_grid: int = 150,
    seed: int = 0,
    patient_id: str = "Visium-GBM",
):
    """Convert a 10x Visium GBM sample into a tmesim ``Microenvironment``.

    Pipeline
    --------
    1. read real spot coordinates (px) + Visium scalefactor -> micrometers;
    2. read real marker-gene expression from the MEX matrix;
    3. phenotype each spot (one representative agent per spot) by marker argmax,
       assign tumor-signal-derived sensitivity;
    4. build a ``Domain`` sized to the tissue extent;
    5. build the vascular-density field with the **real** reconstruct estimator
       (``reconstruct.vasculature.build_vascular_network`` with a tumor mask
       derived from the tumor-spot density) and the ECM field with
       ``reconstruct.spatial.build_tissue_fields``;
    6. bundle via ``reconstruct.assemble.assemble_microenvironment``.

    Parameters
    ----------
    path_or_data : str | dict
        A cache directory containing the 10x tarballs (or extracted ``spatial/``
        and ``filtered_feature_bc_matrix/``), OR a pre-loaded dict with keys
        ``coords_px`` {barcode:(x,y)}, ``umpp`` float, ``barcodes`` list,
        ``marker_expr`` (n_markers,n_spots), ``total_counts`` (n_spots,),
        ``present_genes`` list (mainly for testing / offline use).
    target_grid : approximate max grid dimension (sets ``dx``).
    seed : RNG seed for the vasculature reconstruction.
    patient_id : identifier stored on the microenvironment.

    Returns
    -------
    Microenvironment
    """
    from tmesim.core import CellPopulation, CellState, CellType, Domain
    from tmesim.reconstruct.assemble import assemble_microenvironment
    from tmesim.reconstruct.spatial import build_tissue_fields
    from tmesim.reconstruct.vasculature import build_vascular_network

    # --- gather raw inputs (real bytes or a supplied dict) -------------------
    if isinstance(path_or_data, dict):
        d = path_or_data
        coords_px = d["coords_px"]
        umpp = float(d["umpp"])
        barcodes = list(d["barcodes"])
        marker_expr = np.asarray(d["marker_expr"], dtype=np.float64)
        total_counts = np.asarray(d["total_counts"], dtype=np.float64)
        present_genes = list(d["present_genes"])
        source = "in-memory dict"
    else:
        cache_dir = str(path_or_data)
        spatial_dir, matrix_dir = _ensure_extracted(cache_dir)
        coords_px = _read_positions(spatial_dir)
        umpp = _read_scalefactor_umpp(spatial_dir)
        panel_flat = [g for genes in MARKER_PANEL.values() for g in genes]
        barcodes, marker_expr, total_counts, present_genes = _read_marker_expression(
            matrix_dir, panel_flat
        )
        source = cache_dir

    # --- align spots present in BOTH the matrix and the position list --------
    keep_idx = [i for i, bc in enumerate(barcodes) if bc in coords_px]
    if not keep_idx:
        raise ValueError("no barcodes shared between matrix and tissue positions")
    keep_idx = np.asarray(keep_idx, dtype=int)
    xs_px = np.array([coords_px[barcodes[i]][0] for i in keep_idx], dtype=np.float64)
    ys_px = np.array([coords_px[barcodes[i]][1] for i in keep_idx], dtype=np.float64)
    marker_expr = marker_expr[:, keep_idx] if marker_expr.size else marker_expr
    total_counts = total_counts[keep_idx]

    # --- px -> um, shift to a [0, extent] origin -----------------------------
    x_um = (xs_px - xs_px.min()) * umpp
    y_um = (ys_px - ys_px.min()) * umpp
    width_um = float(x_um.max()) if x_um.size else 1.0
    height_um = float(y_um.max()) if y_um.size else 1.0

    # --- domain sized to the tissue; dx from the target grid -----------------
    dx = max(max(width_um, height_um) / max(target_grid, 1), 5.0)
    nx = int(np.ceil(width_um / dx)) + 2
    ny = int(np.ceil(height_um / dx)) + 2
    domain = Domain(nx=nx, ny=ny, dx=dx)
    ext_w, ext_h = domain.extent_um
    # keep every agent strictly inside the physical extent (assemble asserts this)
    x_um = np.clip(x_um, 0.0, ext_w - 1e-6)
    y_um = np.clip(y_um, 0.0, ext_h - 1e-6)

    # --- phenotype spots -----------------------------------------------------
    cell_type, sensitivity, class_means = _phenotype_spots(
        present_genes, marker_expr, total_counts
    )
    state = np.full(cell_type.shape[0], int(CellState.VIABLE), dtype=np.int8)
    cells = CellPopulation(
        x=x_um, y=y_um, cell_type=cell_type, state=state, sensitivity=sensitivity
    )

    # --- reconstruct continuum fields (REAL reconstruct package) -------------
    # tumor mask for vascular suppression: smoothed tumor-agent density on the grid.
    from scipy.ndimage import gaussian_filter

    tum = cells.tumor_mask()
    tumor_mask = np.zeros(domain.shape, dtype=bool)
    if tum.any():
        grid = domain.zeros()
        col, row = domain.to_index(cells.x[tum], cells.y[tum])
        np.add.at(grid, (row, col), 1.0)
        dens = gaussian_filter(grid, sigma=max(1.0, 60.0 / dx))
        if dens.max() > 0:
            dens = dens / dens.max()
        tumor_mask = dens > 0.30

    vascular = build_vascular_network(domain, seed=seed, tumor_mask=tumor_mask)
    fields = build_tissue_fields(domain, cells, vascular_density=vascular)

    metadata: Dict = {
        "source": "10x Visium Human Glioblastoma (Space Ranger 1.2.0)",
        "provenance_url": VISIUM_GBM_LANDING,
        "ingest_path": source,
        "n_spots": int(cells.n),
        "umpp": umpp,
        "dx_um": dx,
        "marker_panel": MARKER_PANEL,
        "present_genes": present_genes,
        "class_mean_signal": class_means,
        "has_outcome_label": False,
        "spatial_is_real": True,
        "phenotype_is_deconvolved": False,
        "note": (
            "Coordinates + marker expression are MEASURED; one agent per Visium "
            "spot (~1-10 cells pooled); phenotype = marker argmax (not "
            "deconvolution); sensitivity is a heuristic (no drug-response label "
            "exists in open data)."
        ),
    }

    micro = assemble_microenvironment(
        domain, cells, fields, patient_id=patient_id, metadata=metadata
    )
    return micro


# ===========================================================================
# 5. Exercise everything
# ===========================================================================
if __name__ == "__main__":
    import time
    from collections import Counter

    from tmesim.core import CellType

    print("=" * 72)
    print("spatial_ingest: REAL 10x Visium GBM -> reconstruct -> sim")
    print("=" * 72)

    # --- 1. locate / fetch the open dataset ---
    status = find_open_spatial_gbm()
    print("\n[find_open_spatial_gbm]")
    print(f"  dataset : {status['dataset']}")
    print(f"  status  : {status['status']}")
    print(f"  outcome label present : {status['has_outcome_label']}")
    for k, v in status.get("cached", {}).items():
        print(f"  cached  : {k:16s} {v['bytes']:>12,d} B  {v['path']}")
    for k, u in status["urls"].items():
        print(f"  url     : {k:16s} {u}")

    if status["status"] not in ("downloaded", "partial"):
        print("\nNo local bytes; cannot run the real pipeline. URLs above are open.")
        raise SystemExit(0)

    # --- 2. reconstruct a REAL microenvironment ---
    t0 = time.time()
    micro = microenvironment_from_spatial(status["cache_dir"], target_grid=150, seed=0)
    print(f"\n[microenvironment_from_spatial]  ({time.time() - t0:.1f}s)")
    print("  summary:")
    for k, v in micro.summary().items():
        print(f"    {k}: {v}")
    dist = Counter(int(t) for t in micro.cells.cell_type)
    print("  phenotype distribution:",
          {CellType(k).name: v for k, v in sorted(dist.items())})
    print("  class mean signal    :",
          {k: round(v, 3) for k, v in micro.metadata["class_mean_signal"].items()})
    tsens = micro.cells.sensitivity[micro.cells.tumor_mask()]
    if tsens.size:
        print(f"  tumor sensitivity    : mean={tsens.mean():.3f} "
              f"[{tsens.min():.3f}, {tsens.max():.3f}]")

    # --- 3. compile a REAL drug and run the sim on the real-derived geometry ---
    from tmesim.core import PDEConfig, SimConfig
    from tmesim.molecular import compile_drug, library
    from tmesim.sim.engine import simulate_response

    tmz_mol = library.get_molecule("temozolomide")
    drug = compile_drug.molecule_to_drugparams(tmz_mol)  # structure -> DrugParams
    print(f"\n[compiled drug] {drug.name}: eff_D={drug.diffusivity:.1f} um^2/s, "
          f"ic50={drug.ic50} uM, bbb={drug.bbb_penetrance:.3f}")

    pde = PDEConfig(dt=10.0, total_time=3600.0 * 6, convection=True, save_every=200)
    cfg = SimConfig(drug=drug, pde=pde, seed=0)
    t1 = time.time()
    result = simulate_response(micro, cfg)
    print(f"\n[simulate_response]  ({time.time() - t1:.1f}s)  "
          "-- real-derived geometry ACCEPTED by the engine")
    print(f"  predicted_response (kill frac) = {result.predicted_response:.4f}")
    print(f"  viable_tumor_fraction          = {result.viable_tumor_fraction:.4f}")
    print(f"  penetration_index              = {result.extras['penetration_index']:.4f}")
    print(f"  mean_oxygen                    = {result.extras['mean_oxygen']:.4f}")
    print(f"  hypoxic_tumor_fraction         = {result.extras['hypoxic_tumor_fraction']:.4f}")
    print(f"  n_tumor                        = {int(result.extras['n_tumor'])}")

    assert 0.0 <= result.predicted_response <= 1.0
    for key in ("oxygen", "ifp", "auc", "vascular_density"):
        assert result.fields[key].shape == micro.domain.shape
    print("\nEND-TO-END OK: real Visium coords+expression -> reconstruct -> sim.")
