"""Spot deconvolution: real Visium spots -> real cell-type fractions.

This turns each REAL 10x Visium spot's expression profile into a vector of real
cell-type fractions, using the REAL Darmanis GBM single-cell atlas
(GSE84465, ``tmesim.data.real.scref_gbm``) as the deconvolution reference.
There is NO synthetic data anywhere in this module.

Method (deliberately simple + honest):

1. Load the Darmanis genes x cells reference and its per-cell ``cell_type``
   labels. Collapse the 7 published labels onto the coarse ``CellType`` classes
   of the core model:

       Neoplastic                               -> TUMOR
       Immune cell                              -> IMMUNE
       Vascular                                 -> ENDOTHELIAL
       Astrocyte / Oligodendrocyte / OPC / Neuron -> STROMAL

2. Restrict both reference and spot to the SHARED gene set (symbols present in
   both the Visium panel and the reference), library-size-normalize each single
   cell over those shared genes, then average per collapsed class to get one
   mean-expression *signature* per cell type.
3. For each spot, solve ``spot_expr ~= signatures @ w`` with non-negative least
   squares (:func:`scipy.optimize.nnls`) and normalize ``w`` to sum to 1 -> the
   spot's cell-type fractions.

Public API
----------
``build_signatures(shared_genes)``  -> (signatures[G x K], type_names)
``deconvolve_spot(spot_expr, signatures)`` -> fractions[K] (sum 1)
``deconvolve_sample(sample_id)``    -> (coords_um[S,2], fractions[S,K], type_names)
"""
from __future__ import annotations

from typing import List, Sequence, Tuple

import numpy as np
from scipy.optimize import nnls

from ..core.microenvironment import CellType
from ..data.real.scref_gbm import load_reference
from ..data.real.visium import load_sample

# ---------------------------------------------------------------------------
# Darmanis label -> coarse CellType mapping
# ---------------------------------------------------------------------------
# Keys are lower-cased label substrings so the (mis-spelled) "Astocyte" label in
# the GEO metadata is still caught.  Order of the canonical output classes:
_TYPE_ORDER: List[CellType] = [
    CellType.TUMOR,
    CellType.STROMAL,
    CellType.IMMUNE,
    CellType.ENDOTHELIAL,
]


def map_label(label: str) -> "CellType | None":
    """Map a Darmanis ``cell_type`` string onto a coarse core ``CellType``.

    Returns ``None`` for labels that carry no assignment (none currently, but
    keeps the collapse robust to unexpected vocabulary).
    """
    s = str(label).strip().lower()
    if "neoplastic" in s:
        return CellType.TUMOR
    if "immune" in s:
        return CellType.IMMUNE
    if "vascular" in s or "endothel" in s:
        return CellType.ENDOTHELIAL
    if any(k in s for k in ("astocyte", "astrocyte", "oligodendrocyte", "opc", "neuron")):
        return CellType.STROMAL
    return None


def _normalize_columns_to_unit_sum(mat: np.ndarray) -> np.ndarray:
    """Library-size normalize each COLUMN (cell) to sum 1 over the given genes.

    Cells with zero counts over the shared genes are left as all-zero.
    """
    colsums = mat.sum(axis=0)
    out = np.zeros_like(mat, dtype=np.float64)
    nz = colsums > 0
    out[:, nz] = mat[:, nz] / colsums[nz]
    return out


def build_signatures(
    shared_genes: Sequence[str],
    reference: "Tuple | None" = None,
) -> Tuple[np.ndarray, List[str]]:
    """Build per-cell-type mean-expression signatures over ``shared_genes``.

    Parameters
    ----------
    shared_genes : ordered gene symbols. The returned signature rows are aligned
        to this exact order. Symbols not found in the reference contribute an
        all-zero row (they carry no discriminative information).
    reference : optional pre-loaded ``(expr_df, labels_df, gene_names)`` tuple as
        returned by :func:`tmesim.data.real.scref_gbm.load_reference`; loaded
        on demand if omitted.

    Returns
    -------
    signatures : (n_shared_genes, K) float array, one column per cell type.
        Each column is the mean of the library-size-normalized single-cell
        profiles of that class over ``shared_genes`` (columns are on a
        comparable scale, so NNLS weights read as fractions).
    type_names : list[str] of length K — the ``CellType`` names, in the order of
        the signature columns. Only classes actually present in the reference
        are included.
    """
    if reference is None:
        expr_df, labels_df, _ = load_reference()
    else:
        expr_df, labels_df, _ = reference

    # Collapse labels -> coarse CellType (drop unmapped cells).
    mapped = labels_df["cell_type"].map(lambda v: map_label(v))
    keep = mapped.notna()
    mapped = mapped[keep]
    cells_kept = labels_df.index[keep]

    # Restrict reference to shared genes, preserving requested row order. Genes
    # absent from the reference get an all-zero row.
    shared_genes = list(shared_genes)
    ref_index = expr_df.index
    # First occurrence wins for any duplicate symbol in the reference index.
    present = ref_index[~ref_index.duplicated()]
    present_set = set(present)
    ref_sub = expr_df.loc[~ref_index.duplicated()]  # dedupe rows
    ref_sub = ref_sub[cells_kept]                    # keep mapped cells

    G = len(shared_genes)
    row_of = {g: i for i, g in enumerate(shared_genes)}
    # Dense reference submatrix over the shared/found genes.
    found = [g for g in shared_genes if g in present_set]
    found_rows = np.array([row_of[g] for g in found], dtype=int)
    ref_found = ref_sub.loc[found].to_numpy(dtype=np.float64)  # (n_found, n_cells)

    # Library-size normalize each cell over the shared genes (comparable columns).
    ref_norm = _normalize_columns_to_unit_sum(ref_found)

    type_names: List[str] = []
    cols: List[np.ndarray] = []
    mapped_vals = mapped.to_numpy()
    for ct in _TYPE_ORDER:
        sel = mapped_vals == ct
        if not sel.any():
            continue
        mean_found = ref_norm[:, sel].mean(axis=1)  # (n_found,)
        full = np.zeros(G, dtype=np.float64)
        full[found_rows] = mean_found
        cols.append(full)
        type_names.append(ct.name)

    signatures = np.column_stack(cols) if cols else np.zeros((G, 0))
    return signatures, type_names


def deconvolve_spot(spot_expr: np.ndarray, signatures: np.ndarray) -> np.ndarray:
    """Deconvolve one spot's expression into cell-type fractions via NNLS.

    Parameters
    ----------
    spot_expr : (G,) expression vector over the shared gene set (raw counts are
        fine — absolute scale does not affect the normalized fractions).
    signatures : (G, K) signature matrix from :func:`build_signatures`.

    Returns
    -------
    fractions : (K,) non-negative array summing to 1. If NNLS returns an
        all-zero solution (e.g. an empty spot), a uniform vector is returned.
    """
    b = np.asarray(spot_expr, dtype=np.float64).ravel()
    # Library-size normalize the spot to the signature scale (unit sum). This
    # does not change the normalized fractions but keeps NNLS well-conditioned.
    bsum = b.sum()
    if bsum > 0:
        b = b / bsum
    # scipy 1.13's nnls emits spurious FP warnings from its internal residual-norm
    # matmul even when the returned solution is correct; silence just those.
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        w, _ = nnls(signatures, b)
    total = w.sum()
    if total <= 0:
        K = signatures.shape[1]
        return np.full(K, 1.0 / K) if K else w
    return w / total


def deconvolve_sample(
    sample_id: str,
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """Deconvolve every spot of a REAL Visium sample into cell-type fractions.

    Returns
    -------
    coords_um : (S, 2) spot centroids in micrometers (from the loader).
    fractions : (S, K) per-spot cell-type fractions (each row sums to 1).
    type_names : list[str] of length K, aligned to the fraction columns.
    """
    sample = load_sample(sample_id)
    coords_um = np.asarray(sample["coords_um"], dtype=np.float64)
    genes = list(sample["genes"])
    expr = sample["expression"].tocsr()  # (S, G_visium)

    # Shared gene set = Visium symbols also present in the reference.
    expr_df, labels_df, _ = load_reference()
    ref_syms = set(expr_df.index[~expr_df.index.duplicated()])
    # Dedupe Visium symbols (first occurrence) to a well-defined column mapping.
    seen = set()
    shared_genes: List[str] = []
    shared_cols: List[int] = []
    for j, g in enumerate(genes):
        if g in ref_syms and g not in seen:
            seen.add(g)
            shared_genes.append(g)
            shared_cols.append(j)
    if not shared_genes:
        raise ValueError(f"sample {sample_id!r}: no genes shared with reference")

    signatures, type_names = build_signatures(
        shared_genes, reference=(expr_df, labels_df, None)
    )

    # Dense spot x shared-gene matrix, then per-spot NNLS.
    expr_sub = expr[:, np.asarray(shared_cols, dtype=int)].toarray().astype(np.float64)
    S = expr_sub.shape[0]
    K = signatures.shape[1]
    fractions = np.zeros((S, K), dtype=np.float64)
    for i in range(S):
        fractions[i] = deconvolve_spot(expr_sub[i], signatures)

    return coords_um, fractions, type_names


if __name__ == "__main__":
    sid = "10x_parent_gbm"
    coords, fr, names = deconvolve_sample(sid)
    print(f"[{sid}] spots={fr.shape[0]}  types={names}")
    print("mean cell-type fractions across spots:")
    for name, m in zip(names, fr.mean(axis=0)):
        print(f"  {name:12s} {m:.3f}")
