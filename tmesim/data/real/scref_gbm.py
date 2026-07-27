"""Load the REAL GBM single-cell RNA-seq reference atlas.

Source: Darmanis et al. 2017, *Cell Reports* 21(5):1399-1410,
"Single-Cell RNAseq analysis of diffuse neoplastic infiltrating cells at the
migrating front of human glioblastoma."  GEO accession **GSE84465**
(public / open access; NCBI GEO data are freely redistributable).

This is a full-length (SmartSeq2, plate-based) scRNA-seq atlas of 4 primary
human GBM tumours, sampling both the tumour core and the peritumoral
("periphery") brain, with per-cell cell-type annotations published by the
authors in the GEO series-matrix metadata.  It serves as the **deconvolution
reference**: a genes x cells expression matrix plus a real cell-type label per
cell, so spatial (Visium) spots can be resolved into real cell-type fractions
and malignant vs. non-malignant content.

Cached files (see ``docs/real_data/scref_gbm.md``):
  * ``GSE84465_GBM_All_data.csv.gz``     - genes x cells raw HTSeq counts.
  * ``GSE84465_series_matrix.txt.gz``    - per-cell annotations (cell type,
    tissue [Tumor/Periphery], patient id, neoplastic flag, tSNE cluster).

The two are joined on the cell barcode ``<plate_id>.<well>`` which is the column
name in the counts matrix (e.g. ``1001000173.G8``).

Only scipy/pandas/stdlib are used; no scanpy/anndata required.
"""
from __future__ import annotations

import gzip
import os
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# Default cache location (repo-root/data_cache/scref_gbm).
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", ".."))
_DEFAULT_CACHE = os.path.join(_REPO_ROOT, "data_cache", "scref_gbm")

_COUNTS_FILE = "GSE84465_GBM_All_data.csv.gz"
_META_FILE = "GSE84465_series_matrix.txt.gz"

# HTSeq summary rows that trail the real gene rows in the counts CSV; not genes.
_NON_GENE_ROWS = {
    "no_feature",
    "ambiguous",
    "too_low_aQual",
    "not_aligned",
    "alignment_not_unique",
}


def _parse_series_matrix(path: str) -> pd.DataFrame:
    """Parse per-cell annotations from a GEO series-matrix .txt.gz.

    Returns a DataFrame indexed by the ``<plate_id>.<well>`` cell barcode with
    columns: cell_type, tissue, patient_id, neoplastic, tsne_cluster, plate_id,
    well, gsm.
    """
    titles: List[str] = []
    gsms: List[str] = []
    # characteristic key -> list of per-sample values
    chars: Dict[str, List[str]] = {}

    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                titles = [c.strip('"') for c in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                gsms = [c.strip('"') for c in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                vals = [c.strip('"') for c in line.rstrip("\n").split("\t")[1:]]
                # Each characteristics line is "key: value"; key is stable per line.
                key = None
                for v in vals:
                    if ": " in v:
                        key = v.split(": ", 1)[0].strip()
                        break
                if key is None:
                    continue
                parsed = [v.split(": ", 1)[1] if ": " in v else "" for v in vals]
                # A GSE can repeat the same key across lines; keep them distinct.
                base = key
                n = 1
                while key in chars:
                    n += 1
                    key = f"{base}_{n}"
                chars[key] = parsed

    n_samples = len(titles)
    plate = chars.get("plate id", [""] * n_samples)
    well = chars.get("well", [""] * n_samples)
    barcodes = [f"{p}.{w}" for p, w in zip(plate, well)]

    df = pd.DataFrame(
        {
            "cell_type": chars.get("cell type", [""] * n_samples),
            "tissue": chars.get("tissue", [""] * n_samples),
            "patient_id": chars.get("patient id", [""] * n_samples),
            "neoplastic": chars.get("neoplastic", [""] * n_samples),
            "tsne_cluster": chars.get("tsne cluster", [""] * n_samples),
            "plate_id": plate,
            "well": well,
            "gsm": gsms if len(gsms) == n_samples else [""] * n_samples,
        },
        index=barcodes,
    )
    df.index.name = "cell_id"
    return df


def load_reference(
    cache_dir: str = _DEFAULT_CACHE,
    max_cells: Optional[int] = None,
    random_state: int = 0,
) -> Tuple[pd.DataFrame, pd.DataFrame, np.ndarray]:
    """Load the Darmanis GBM scRNA-seq reference from the cache.

    Parameters
    ----------
    cache_dir : str
        Directory holding the two cached GSE84465 files.
    max_cells : int, optional
        If given, return a random subsample of this many cells (stratified is
        not applied; a documented uniform subsample). ``None`` returns all.
    random_state : int
        Seed for the optional subsample.

    Returns
    -------
    expression : pandas.DataFrame
        genes x cells raw counts (index = gene symbol, columns = cell barcode).
    cell_labels : pandas.DataFrame
        per-cell annotations aligned to ``expression.columns`` (cell_type,
        tissue, patient_id, neoplastic, tsne_cluster, ...).
    gene_names : numpy.ndarray
        1-D array of gene symbols (== ``expression.index``).

    Raises
    ------
    FileNotFoundError
        If the cache files are missing (run the acquisition step first).
    """
    counts_path = os.path.join(cache_dir, _COUNTS_FILE)
    meta_path = os.path.join(cache_dir, _META_FILE)
    for p in (counts_path, meta_path):
        if not os.path.exists(p):
            raise FileNotFoundError(
                f"Missing cached reference file: {p}. "
                "Download GSE84465 into data_cache/scref_gbm/ first "
                "(see docs/real_data/scref_gbm.md)."
            )

    # Counts: first column is the gene symbol (unnamed in the header). Despite
    # the ``.csv`` name the file is whitespace-delimited.
    expr = pd.read_csv(counts_path, sep=r"\s+", index_col=0)
    expr.index = expr.index.astype(str)
    # Drop trailing HTSeq summary rows so only real genes remain.
    expr = expr.loc[~expr.index.isin(_NON_GENE_ROWS)]

    meta = _parse_series_matrix(meta_path)

    # Align annotations to the counts columns (intersection; both are the same
    # <plate>.<well> barcodes in practice).
    common = [c for c in expr.columns if c in meta.index]
    expr = expr[common]
    labels = meta.loc[common].copy()

    if max_cells is not None and max_cells < expr.shape[1]:
        rng = np.random.default_rng(random_state)
        keep = rng.choice(expr.shape[1], size=max_cells, replace=False)
        keep.sort()
        cols = expr.columns[keep]
        expr = expr[cols]
        labels = labels.loc[cols]

    gene_names = np.asarray(expr.index, dtype=object)
    return expr, labels, gene_names


if __name__ == "__main__":
    X, y, genes = load_reference()
    print(f"expression: {X.shape[0]} genes x {X.shape[1]} cells")
    print(f"cell-type vocabulary: {sorted(y['cell_type'].unique())}")
    print(y["cell_type"].value_counts().to_string())
    print(f"tissue split: {dict(y['tissue'].value_counts())}")
    print(f"patients: {sorted(y['patient_id'].unique())}")
