"""Loader for GEO **GSE324268** — a human primary GBM Xenium section that
**does** measure MGMT.

Why this module exists
----------------------
``results/gbm_tmz_exvivo.md`` §B2 and ``results/mgmt_scalar.md`` §2 concluded
that MGMT is not spatially resolvable, on the strength of two facts:

  1. MGMT is absent from both Xenium panels this repo held (319-gene brain,
     480-gene immuno-oncology; 739-gene union, ``MGMT present = False``); and
  2. in the four Visium sections MGMT sits at the detection floor
     (5.3-29.4 % of spots, ~1-2 UMI each), so only a pseudobulk scalar is
     honest.

Fact 1 is a statement about **those two panels**, and it was over-generalised
into "MGMT cannot be resolved spatially". This dataset falsifies the
generalisation: it is a human primary GBM Xenium run whose custom add-on panel
carries a dedicated 8-probe MGMT target.

Dataset (OPEN, free direct download, no login, no DUA)
------------------------------------------------------
GEO **GSE324268** — "Spatial transcriptomics of human glioblastoma and of a
GBM-organoid-on-brain-slice model" (Hölzel lab, Bonn; public 2026-03).

This module loads the **human primary tumour** sample only:

  * **GSM9572843** — ``GBMtumor07``, 851,509 segmented cells.

Panel (read from the shipped ``GSE324268_gene_panel.json.gz``, not assumed):

  * base ``Xenium Human Brain Gene Expression`` **hBrain_v1.1**, 266 targets;
  * add-on ``hBrain_100g`` design_id **ECQNGH**, described in the file as
    ``"Glioblastoma Panel "``, 100 targets, designed 2024-01-08;
  * 400 target entries, of which **366** appear as ``Gene Expression``
    features in the shipped count matrix (the remainder resolve to unassigned
    / deprecated codewords).

**MGMT** (``ENSG00000170430``) is one of the add-on targets, with
``info.gene_coverage = 8`` — eight probe pairs on one gene target, i.e. ONE
column in the feature matrix, not eight.

The second sample in the series (GSM9572842, ``GBM1organoid``) is a
GBM-organoid-on-mouse-brain-slice model and is deliberately NOT loaded here:
it is not human primary tissue.

Files this module expects in ``path`` (defaults to
``data_cache/gse324268_xenium``)
--------------------------------------------------------------------------
Downloaded by :func:`download` from the NCBI GEO FTP mirror; sizes and
sha256s in :data:`GSE324268_FILES` were verified against the bytes on disk.
The multi-GB morphology OME-TIFFs and the 3.3 GB transcript table are
deliberately NOT fetched — the cell x gene matrix plus the cell table is the
smallest complete unit for a cell-resolved reconstruction.

  * ``GSM9572843_GBMtumor07_cell_feature_matrix.h5``  (10x CellRanger-style
    HDF5, CSC, 541 features x 851,509 cells)
  * ``GSM9572843_GBMtumor07_cells.parquet.gz``  (per-cell centroids + QC)
  * ``GSE324268_gene_panel.json.gz``  (the panel spec quoted above)

Reading the parquet needs ``pyarrow``; :func:`download` therefore also writes a
``.csv.gz`` transcode of the cell table so the rest of the module (and the
existing Xenium loader) needs only pandas. :func:`load` prefers the CSV and
falls back to the parquet.

Two entry points
----------------
:func:`load`
    Returns the SAME dict contract as
    :func:`tmesim.data.real.singlecell_spatial.load` — ``cell_xy_um``,
    ``expression`` (cells x genes), ``genes``, ``cell_ids``, ``cell_meta``,
    ``coordinate_units``, ``source`` — read straight from the HDF5, which is
    fast (~10 s) and needs no intermediate files.

:func:`materialize_xenium_cache`
    Writes the 10x MEX layout (``cell_feature_matrix/{matrix.mtx.gz,
    features.tsv.gz,barcodes.tsv.gz}``, ``Xenium_GBM_cells.csv.gz``,
    ``Xenium_GBM_gene_panel.json``) that the EXISTING loader and
    ``reconstruct.real_xenium.reconstruct_from_xenium`` already consume, so
    this section can be reconstructed by the unmodified pipeline. Idempotent.
    This is the reason no file under ``tmesim/reconstruct/`` is touched.

Panel selection is also unchanged and needs no new registry entry: this panel
carries **all twelve** markers of ``marker_panels.BRAIN_PANEL_MARKER_MAP``
(EGFR, SOX2, OLIG2, MKI67, PTPRC, CD68, AIF1, CD163, PECAM1, MOG, AQP4, GJA1),
so ``select_marker_map`` scores it brain_319 = 1.00 and picks the existing
brain map. Verified in :func:`panel_marker_coverage`.

Honest limits
-------------
* **Different patient** from every other section in this repo. Nothing here is
  a paired measurement on tissue used elsewhere.
* **No outcome label.** GSE324268 is an untreated resection; there is no drug
  response, no PFS/OS and no MGMT promoter-methylation call in the deposit.
* Xenium measures MGMT **transcript**, not promoter methylation. See
  ``results/mgmt_spatial.md`` for why expression is nonetheless the quantity
  the model wants.
"""
from __future__ import annotations

import gzip
import hashlib
import io
import json
import os
import shutil
import sys
import urllib.request
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import sparse

# ---------------------------------------------------------------------------
# Provenance. URLs verified HTTP 200 on 2026-07-28; byte sizes and sha256s are
# of the files actually downloaded and used by ``results/mgmt_spatial.md``.
# ---------------------------------------------------------------------------
GSE324268_ACCESSION = "GSE324268"
GSE324268_SAMPLE = "GSM9572843"
GSE324268_LANDING = (
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE324268"
)
_GEO_SAMPLE_BASE = (
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9572nnn/GSM9572843/suppl"
)
_GEO_SERIES_BASE = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE324nnn/GSE324268/suppl"
)

#: local basename -> (url, verified bytes, verified sha256)
GSE324268_FILES: Dict[str, Tuple[str, int, str]] = {
    "GSM9572843_GBMtumor07_cell_feature_matrix.h5": (
        f"{_GEO_SAMPLE_BASE}/GSM9572843_GBMtumor07_cell_feature_matrix.h5",
        61_566_910,
        "816934dbae8e6b41b4cd3072092817c0387be88ef5dd77b186fea47edce43b5c",
    ),
    "GSM9572843_GBMtumor07_cells.parquet.gz": (
        f"{_GEO_SAMPLE_BASE}/GSM9572843_GBMtumor07_cells.parquet.gz",
        15_533_090,
        "571bfd1d6914d939c4d00c7fca7e19e9399604138883eb92326c7f4148ca283c",
    ),
    "GSE324268_gene_panel.json.gz": (
        f"{_GEO_SERIES_BASE}/GSE324268_gene_panel.json.gz",
        7_263,
        "b2784a12b606dfbd0672c24d8398af9df699680503809e4b84651fa991aff6e9",
    ),
}

#: Files in the deposit that are deliberately NOT downloaded, and why.
GSE324268_SKIPPED = {
    "GSM9572843_GBMtumor07_morphology.ome.tif.gz": "19.96 GB image bundle — not needed for a cell-resolved reconstruction",
    "GSM9572843_GBMtumor07_transcripts.parquet.gz": "3.34 GB per-transcript table — the cell x gene matrix is the aggregated form of it",
    "GSM9572843_GBMtumor07_cell_boundaries.parquet.gz": "74 MB polygon table — centroids suffice; the pipeline consumes centroids",
    "GSM9572843_GBMtumor07_nucleus_boundaries.parquet.gz": "66 MB polygon table — as above",
    "GSM9572842_GBM1organoid_*": "organoid-on-brain-slice model, not human primary tissue",
}

DEFAULT_CACHE = os.path.join("data_cache", "gse324268_xenium")

_H5 = "GSM9572843_GBMtumor07_cell_feature_matrix.h5"
_CELLS_PARQUET = "GSM9572843_GBMtumor07_cells.parquet.gz"
_CELLS_CSV = "GSM9572843_GBMtumor07_cells.csv.gz"
_PANEL_GZ = "GSE324268_gene_panel.json.gz"

# 10x MEX layout written by ``materialize_xenium_cache``; basenames must match
# ``tmesim.data.real.singlecell_spatial`` exactly.
_MEX_DIR = "cell_feature_matrix"
_MEX_CELLS = "Xenium_GBM_cells.csv.gz"
_MEX_PANEL = "Xenium_GBM_gene_panel.json"

_GENE_TYPE = "Gene Expression"

#: The gene this whole module exists for, and its Ensembl id as shipped.
MGMT_SYMBOL = "MGMT"
MGMT_ENSEMBL = "ENSG00000170430"


# ---------------------------------------------------------------------------
# download
# ---------------------------------------------------------------------------
def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(path: str = DEFAULT_CACHE, verify: bool = True) -> Dict[str, Dict]:
    """Fetch the three GSE324268 files into ``path`` (idempotent).

    Retries once on a transient network error, which is the documented policy
    for every network path in this repo.

    Returns a per-file record of bytes and sha256, and whether they match the
    values frozen in :data:`GSE324268_FILES`.
    """
    os.makedirs(path, exist_ok=True)
    out: Dict[str, Dict] = {}
    for name, (url, want_bytes, want_sha) in GSE324268_FILES.items():
        dest = os.path.join(path, name)
        if not os.path.exists(dest):
            try:
                _fetch(url, dest)
            except Exception:
                _fetch(url, dest)  # retry once
        got_bytes = os.path.getsize(dest)
        got_sha = _sha256(dest) if verify else None
        out[name] = {
            "url": url,
            "bytes": got_bytes,
            "bytes_match": got_bytes == want_bytes,
            "sha256": got_sha,
            "sha256_match": (got_sha == want_sha) if verify else None,
        }
    # Transcode the parquet cell table to csv.gz so nothing downstream needs
    # pyarrow. Only done once.
    csv_path = os.path.join(path, _CELLS_CSV)
    if not os.path.exists(csv_path):
        df = _read_cells_parquet(os.path.join(path, _CELLS_PARQUET))
        df.to_csv(csv_path, index=False, compression="gzip")
    out[_CELLS_CSV] = {
        "url": "derived (pyarrow transcode of the .parquet.gz)",
        "bytes": os.path.getsize(csv_path),
        "bytes_match": None,
        "sha256": _sha256(csv_path) if verify else None,
        "sha256_match": None,
    }
    return out


def _fetch(url: str, dest: str) -> None:
    tmp = dest + ".part"
    with urllib.request.urlopen(url, timeout=120) as r, open(tmp, "wb") as fh:
        shutil.copyfileobj(r, fh, length=1 << 20)
    os.replace(tmp, dest)


def _read_cells_parquet(parquet_gz: str) -> pd.DataFrame:
    with gzip.open(parquet_gz, "rb") as fh:
        return pd.read_parquet(io.BytesIO(fh.read()))


# ---------------------------------------------------------------------------
# panel
# ---------------------------------------------------------------------------
def panel_spec(path: str = DEFAULT_CACHE) -> Dict:
    """Return the ``payload`` block of the shipped gene-panel JSON."""
    with gzip.open(os.path.join(path, _PANEL_GZ), "rt") as fh:
        return json.load(fh)["payload"]


def mgmt_panel_entry(path: str = DEFAULT_CACHE) -> Dict:
    """The MGMT target record from the panel file — probe count and provenance.

    Verifying the claim "MGMT is on this panel" from the panel spec itself,
    rather than from a README, is the point of this function.
    """
    for t in panel_spec(path)["targets"]:
        data = t.get("type", {}).get("data", {})
        if data.get("name") == MGMT_SYMBOL:
            return {
                "symbol": data.get("name"),
                "ensembl": data.get("id"),
                "gene_coverage": t.get("info", {}).get("gene_coverage"),
                "codewords": t.get("codewords"),
                "source_panel": t.get("source", {}).get("identity", {}),
            }
    raise KeyError(f"{MGMT_SYMBOL} not present in the panel spec at {path}")


def panel_marker_coverage(path: str = DEFAULT_CACHE) -> Dict:
    """How the EXISTING panel registry scores this panel.

    Returned so ``results/mgmt_spatial.md`` can state, from a measurement
    rather than an assertion, that no new entry in
    ``reconstruct.marker_panels.PANEL_REGISTRY`` is required and that the two
    previously-cached sections are therefore untouched.
    """
    from ...reconstruct.marker_panels import (
        BRAIN_PANEL_MARKER_MAP,
        select_marker_map,
    )

    genes = load_genes(path)
    mm, name, cov = select_marker_map(genes)
    present = set(genes.tolist())
    return {
        "selected_panel": name,
        "coverage": cov,
        "is_brain_map_object": mm is BRAIN_PANEL_MARKER_MAP,
        "markers_present": {
            ct.name: [m for m in ms if m in present] for ct, ms in mm.items()
        },
        "markers_absent": sorted(
            {m for ms in mm.values() for m in ms if m not in present}
        ),
    }


# ---------------------------------------------------------------------------
# load
# ---------------------------------------------------------------------------
def _h5_features(path: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    import h5py

    with h5py.File(os.path.join(path, _H5), "r") as f:
        m = f["matrix"]
        name = np.array([s.decode() for s in m["features/name"][:]])
        ident = np.array([s.decode() for s in m["features/id"][:]])
        ftype = np.array([s.decode() for s in m["features/feature_type"][:]])
    return name, ident, ftype


def load_genes(path: str = DEFAULT_CACHE) -> np.ndarray:
    """Just the ``Gene Expression`` symbols — cheap, no matrix read."""
    name, _ident, ftype = _h5_features(path)
    return name[ftype == _GENE_TYPE]


def _source(path: str, n_genes: int) -> Dict[str, object]:
    return {
        "dataset": (
            "GEO GSE324268 — human primary glioblastoma, 10x Xenium In Situ "
            "(Hölzel lab, Bonn; public 2026-03)"
        ),
        "accession": f"{GSE324268_ACCESSION} / {GSE324268_SAMPLE}",
        "sample": "GBMtumor07",
        "assay": "Xenium In Situ (single-cell-resolution spatial)",
        "panel": (
            "Xenium Human Brain Gene Expression base (hBrain_v1.1, 266 targets) "
            "+ custom add-on 'Glioblastoma Panel' (hBrain_100g, design_id "
            f"ECQNGH, 100 targets) = {n_genes} Gene Expression features. "
            "MGMT (ENSG00000170430) is on the add-on with 8 probe pairs."
        ),
        "license": "GEO public deposit — free direct download, no login, no DUA",
        "landing_page": GSE324268_LANDING,
        "has_outcome_label": False,
        "measures_mgmt": True,
        "cache_dir": path,
        "n_panel_genes": int(n_genes),
    }


def load(
    path: str = DEFAULT_CACHE,
    genes_only: bool = True,
    dense: bool = False,
) -> Dict[str, object]:
    """Load GSE324268 ``GBMtumor07`` in the standard Xenium dict contract.

    Identical keys, orientation and units to
    :func:`tmesim.data.real.singlecell_spatial.load`, so any consumer of that
    loader can consume this one.

    Parameters
    ----------
    path : str
        Cache directory (see :func:`download`).
    genes_only : bool
        Restrict to the 366 ``Gene Expression`` features, dropping the 168
        negative-control / unassigned / deprecated codeword features. Pass
        ``False`` to keep them — the negative controls are the background
        reference the MGMT detection claim is scored against.
    dense : bool
        Densify ``expression``. **851,509 x 366 float32 is ~1.2 GB**; the
        default CSR is strongly preferred.

    Raises
    ------
    FileNotFoundError
        If the cache is absent. Call :func:`download` first.
    """
    import h5py

    h5_path = os.path.join(path, _H5)
    if not os.path.exists(h5_path):
        raise FileNotFoundError(
            f"Missing {h5_path}. Run "
            "`python -m tmesim.data.real.mgmt_spatial download` first."
        )

    with h5py.File(h5_path, "r") as f:
        m = f["matrix"]
        n_feat, n_cell = (int(v) for v in m["shape"][:])
        X = sparse.csc_matrix(
            (m["data"][:], m["indices"][:], m["indptr"][:]),
            shape=(n_feat, n_cell),
        )
        barcodes = np.array([s.decode() for s in m["barcodes"][:]])
        name = np.array([s.decode() for s in m["features/name"][:]])
        ftype = np.array([s.decode() for s in m["features/feature_type"][:]])

    cells = _load_cell_table(path)
    if not (len(cells) == n_cell == len(barcodes)):
        raise ValueError(
            f"cell-count mismatch: cells table={len(cells)}, "
            f"matrix cols={n_cell}, barcodes={len(barcodes)}"
        )
    if not (cells["cell_id"].to_numpy().astype(str) == barcodes).all():
        order = {cid: j for j, cid in enumerate(barcodes)}
        X = X[:, cells["cell_id"].map(order).to_numpy()]

    keep = (ftype == _GENE_TYPE) if genes_only else np.ones(n_feat, bool)
    expr = X[keep, :].T.tocsr()  # cells x features
    genes = name[keep]
    if dense:
        expr = expr.toarray().astype(np.float32)

    return {
        "cell_xy_um": cells[["x_centroid", "y_centroid"]].to_numpy(np.float64),
        "expression": expr,
        "genes": genes,
        "feature_types": ftype[keep],
        "cell_ids": cells["cell_id"].to_numpy().astype(str),
        "cell_meta": cells,
        "coordinate_units": "micrometer",
        "source": _source(path, int((ftype == _GENE_TYPE).sum())),
    }


def _load_cell_table(path: str) -> pd.DataFrame:
    csv_path = os.path.join(path, _CELLS_CSV)
    if os.path.exists(csv_path):
        return pd.read_csv(csv_path)
    return _read_cells_parquet(os.path.join(path, _CELLS_PARQUET))


# ---------------------------------------------------------------------------
# MGMT extraction
# ---------------------------------------------------------------------------
def mgmt_per_cell(path: str = DEFAULT_CACHE) -> Dict[str, object]:
    """Per-cell MGMT counts plus everything needed to judge them as a signal.

    Returns
    -------
    dict with
        mgmt : ndarray (n_cells,) int  — MGMT counts per cell
        gene_total : ndarray (n_cells,) — total ``Gene Expression`` counts
        neg_probe_total / neg_codeword_total : ndarray (n_cells,)
        xy : ndarray (n_cells, 2) float — centroids in um
        cell_area / nucleus_area : ndarray (n_cells,)
        gene_totals : DataFrame — per-gene section totals (background ranking)
        neg_probe_totals / neg_codeword_totals : ndarray — per-FEATURE totals
        n_mgmt_probes : int — probe pairs on the MGMT target (from the panel)
    """
    d = load(path, genes_only=False)
    expr = d["expression"]           # cells x all features, CSR
    genes = np.asarray(d["genes"])   # all feature names
    ftype = np.asarray(d["feature_types"])

    j = np.where((genes == MGMT_SYMBOL) & (ftype == _GENE_TYPE))[0]
    if len(j) != 1:
        raise KeyError(
            f"expected exactly one {MGMT_SYMBOL} Gene Expression feature, "
            f"found {len(j)}"
        )
    mgmt = np.asarray(expr[:, j[0]].todense()).ravel().astype(np.int32)

    is_gene = ftype == _GENE_TYPE
    is_np = ftype == "Negative Control Probe"
    is_nc = ftype == "Negative Control Codeword"

    per_cell = lambda mask: np.asarray(  # noqa: E731
        expr[:, np.where(mask)[0]].sum(axis=1)
    ).ravel()
    per_feature = lambda mask: np.asarray(  # noqa: E731
        expr[:, np.where(mask)[0]].sum(axis=0)
    ).ravel()

    gene_tot = per_feature(is_gene)
    meta = d["cell_meta"]
    return {
        "mgmt": mgmt,
        "gene_total": per_cell(is_gene),
        "neg_probe_total": per_cell(is_np),
        "neg_codeword_total": per_cell(is_nc),
        "xy": d["cell_xy_um"],
        "cell_area": meta["cell_area"].to_numpy(np.float64),
        "nucleus_area": meta["nucleus_area"].to_numpy(np.float64),
        "gene_totals": pd.DataFrame(
            {"gene": genes[is_gene], "total": gene_tot}
        ).sort_values("total", ascending=False, ignore_index=True),
        "neg_probe_totals": per_feature(is_np),
        "neg_codeword_totals": per_feature(is_nc),
        "n_neg_probes": int(is_np.sum()),
        "n_neg_codewords": int(is_nc.sum()),
        "n_mgmt_probes": int(mgmt_panel_entry(path)["gene_coverage"]),
        "expression": expr,
        "all_genes": genes,
        "feature_types": ftype,
        "source": d["source"],
    }


# ---------------------------------------------------------------------------
# MEX materialisation, so the UNMODIFIED reconstruct pipeline can read it
# ---------------------------------------------------------------------------
def materialize_xenium_cache(path: str = DEFAULT_CACHE, force: bool = False) -> Dict:
    """Write the 10x MEX layout the existing Xenium loader expects.

    ``tmesim.data.real.singlecell_spatial.load`` (and therefore
    ``reconstruct.real_xenium.reconstruct_from_xenium``) reads
    ``Xenium_GBM_cells.csv.gz``, ``Xenium_GBM_gene_panel.json`` and
    ``cell_feature_matrix/{matrix.mtx.gz,features.tsv.gz,barcodes.tsv.gz}``.
    This writes exactly those, in place, from the GEO HDF5 — so the section can
    be reconstructed **by the unmodified pipeline** and no file under
    ``tmesim/reconstruct/`` needs to change.

    All 541 features are written (not just the 366 genes) so the on-disk MEX is
    a faithful copy of the deposit; the loader's own ``genes_only`` filter then
    behaves exactly as it does for the other two sections.

    Idempotent: returns immediately if the layout is already present.
    """
    import h5py
    from scipy.io import mmwrite

    mex_dir = os.path.join(path, _MEX_DIR)
    mtx = os.path.join(mex_dir, "matrix.mtx.gz")
    done = (
        os.path.exists(mtx)
        and os.path.exists(os.path.join(path, _MEX_CELLS))
        and os.path.exists(os.path.join(path, _MEX_PANEL))
    )
    if done and not force:
        return {"materialized": False, "reason": "already present", "path": path}

    os.makedirs(mex_dir, exist_ok=True)

    with h5py.File(os.path.join(path, _H5), "r") as f:
        m = f["matrix"]
        n_feat, n_cell = (int(v) for v in m["shape"][:])
        X = sparse.csc_matrix(
            (m["data"][:], m["indices"][:], m["indptr"][:]),
            shape=(n_feat, n_cell),
        )
        barcodes = [s.decode() for s in m["barcodes"][:]]
        name = [s.decode() for s in m["features/name"][:]]
        ident = [s.decode() for s in m["features/id"][:]]
        ftype = [s.decode() for s in m["features/feature_type"][:]]

    buf = io.BytesIO()
    mmwrite(buf, X.tocoo(), field="integer")
    with gzip.open(mtx, "wb", compresslevel=1) as fh:
        fh.write(buf.getvalue())
    del buf

    with gzip.open(os.path.join(mex_dir, "features.tsv.gz"), "wt") as fh:
        for i, s, t in zip(ident, name, ftype):
            fh.write(f"{i}\t{s}\t{t}\n")
    with gzip.open(os.path.join(mex_dir, "barcodes.tsv.gz"), "wt") as fh:
        fh.write("\n".join(barcodes) + "\n")

    cells = _load_cell_table(path)
    cells.to_csv(os.path.join(path, _MEX_CELLS), index=False, compression="gzip")

    with gzip.open(os.path.join(path, _PANEL_GZ), "rt") as fh:
        panel = json.load(fh)
    with open(os.path.join(path, _MEX_PANEL), "w") as fh:
        json.dump(panel, fh)

    return {
        "materialized": True,
        "path": path,
        "n_features": n_feat,
        "n_cells": n_cell,
        "nnz": int(X.nnz),
        "matrix_bytes": os.path.getsize(mtx),
    }


def reconstruct(path: str = DEFAULT_CACHE, dx_um: float = 30.0, **kwargs):
    """Reconstruct this section with the **unmodified** existing pipeline.

    Materialises the MEX layout if needed, then calls
    ``reconstruct.real_xenium.reconstruct_from_xenium(sample_path=path)``.
    Nothing in ``reconstruct/`` is patched, subclassed or shadowed.
    """
    from ...reconstruct.real_xenium import reconstruct_from_xenium

    materialize_xenium_cache(path)
    return reconstruct_from_xenium(sample_path=path, dx_um=dx_um, **kwargs)


if __name__ == "__main__":  # pragma: no cover - smoke test
    if len(sys.argv) > 1 and sys.argv[1] == "download":
        for k, v in download().items():
            print(f"{k}\n   {v['bytes']:,} bytes  sha256 {str(v['sha256'])[:16]}  "
                  f"match={v['sha256_match']}")
        sys.exit(0)

    e = mgmt_panel_entry()
    print("MGMT panel entry:", e)
    d = load()
    print(f"cells        : {d['cell_xy_um'].shape[0]:,}")
    print(f"gene targets : {len(d['genes'])}")
    print(f"panel        : {d['source']['panel']}")
    print("marker coverage:", panel_marker_coverage())
