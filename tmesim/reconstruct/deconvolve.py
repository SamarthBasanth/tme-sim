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

from collections import OrderedDict
from typing import Dict, List, Sequence, Tuple

import numpy as np
from scipy.optimize import nnls

from ..core.microenvironment import CellType
from ..data.real.scref_gbm import (_OVERRIDE_ENV, load_reference,
                                   memory_headroom, reference_fingerprint,
                                   reference_gene_symbols)
from ..data.real.visium import load_sample

# --------------------------------------------------------------------------- #
# What a Visium deconvolution actually costs, and the refusal built on it
# --------------------------------------------------------------------------- #
# MEASURED 2026-09-20, darwin 24.6.0 / CPython 3.9 / pandas 2.3.3 / numpy
# 2.0.2, peak RSS rise of one cold `deconvolve_sample(sample_id)` per process
# (``resource.getrusage(RUSAGE_SELF).ru_maxrss``, which is BYTES on darwin and
# KiB on Linux). Three or more cold runs per section; the MAXIMUM observed is
# recorded. Full table and method in results/scref_memory_2026-09-20.md.
#
# THE POINT OF THIS TABLE. Fixing the reference loader was not enough and the
# numbers say so plainly: `load_reference()` went from a measured 1,961 MiB
# peak to 610 MiB, and the figures below barely moved, because the dominant
# term is `build_signatures` densifying a ~19,000 x 3,589 submatrix to float64
# TWICE. A Visium section needs ~3 GB on this machine whatever the reference
# loader does, so a 512 MB container cannot serve one at any `dx_um`. That is
# what the refusal below has to say, and it must not imply otherwise.
MEASURED_VISIUM_PEAK_MIB: Dict[str, float] = {
    "GBM1": 3454.8,
    "GBM2": 3151.8,
    "GBM3": 3156.1,
    "10x_parent_gbm": 3146.1,
}

#: The CHEAPEST Visium deconvolution ever measured here, across every section
#: and every run (GBM1, post-memo run: 2,924.0 MiB). Used as the threshold for
#: a sample with no measurement of its own -- an ingested section, say. It is
#: deliberately the minimum and not the mean: this is the figure below which
#: NOTHING measured here has ever fitted, so refusing under it never refuses a
#: run that was observed to succeed.
MEASURED_VISIUM_MIN_MIB = 2924.0


class DeconvolutionTooLargeForHost(MemoryError):
    """A Visium deconvolution will not fit in this host's declared memory.

    ``MemoryError`` subclass for the same reason as
    ``scref_gbm.ReferenceTooLargeForHost``: it IS the memory error, raised
    before the allocation rather than by the OOM killer afterwards.

    NOT YET WIRED TO THE CLIENT. ``app/server/jobs.py`` renders a failed run's
    error verbatim only for classes named in ``RUN_REFUSALS``; this one is not
    in that set, so a visitor currently gets the redacted sentence plus a
    reference id while the text below goes to stderr and
    ``GET /api/runs/{id}/traceback``. The container is saved either way -- the
    refusal happens before the allocation -- but the reader is not served.
    Adding this name and ``ReferenceTooLargeForHost`` to ``RUN_REFUSALS``
    closes it; see results/scref_memory_2026-09-20.md §4.5.
    """


def check_visium_memory_budget(sample_id: str) -> None:
    """Refuse a Visium deconvolution this host cannot hold. See the table above.

    This is the refusal that ``TMESIM_MAX_GRID_CELLS`` cannot make. That
    ceiling counts VOXELS, and every megabyte measured above is spent before a
    grid exists -- on the Visium matrix, the scRNA reference and the signature
    densification. Shrinking ``dx_um`` does not move any of it, so a visitor
    who picks a Visium section walks straight past the grid ceiling.

    Refuses only against a DECLARED limit (cgroup v2 -> v1 ->
    ``TMESIM_HOST_MEMORY_MB``, read by ``app/server/account.py::memory_limit``).
    Undeclared memory raises nothing.
    """
    h = memory_headroom()
    if h["override"] or h["headroom_bytes"] is None:
        return
    mib = 1024.0 * 1024.0
    need_mib = MEASURED_VISIUM_PEAK_MIB.get(sample_id, MEASURED_VISIUM_MIN_MIB)
    have_mib = h["headroom_bytes"] / mib
    if need_mib <= have_mib:
        return

    measured = (f"MEASURED at {need_mib:,.0f} MiB for this section"
                if sample_id in MEASURED_VISIUM_PEAK_MIB else
                f"not measured for this section; the cheapest Visium section "
                f"measured here needed {need_mib:,.0f} MiB")
    raise DeconvolutionTooLargeForHost(
        f"Refusing to deconvolve Visium sample {sample_id!r}: peak memory is "
        f"{measured}, and this host has {have_mib:,.0f} MiB of headroom "
        f"({h['limit_bytes'] / mib:,.0f} MiB limit, source {h['source']}, minus "
        f"{h['rss_bytes'] / mib:,.0f} MiB already resident)."
        f"\n\nThis is NOT a grid-size problem. TMESIM_MAX_GRID_CELLS counts "
        f"voxels and every byte above is spent before a grid exists, so no "
        f"choice of dx_um avoids it. The measured cost is dominated by "
        f"densifying the scRNA reference to float64 in build_signatures, not "
        f"by the requested resolution."
        f"\n\nOptions: (1) choose a Xenium section -- the Xenium path does not "
        f"load this reference at all, so none of the cost above applies to it "
        f"(its own cost was not measured in this pass); "
        f"(2) run Visium on an instance with more memory; or (3) set "
        f"{_OVERRIDE_ENV}=1 to attempt it anyway and accept that the container "
        f"may be OOM-killed. Measurements: "
        f"results/scref_memory_2026-09-20.md."
    )

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


# --------------------------------------------------------------------------- #
# Signature memoisation
# --------------------------------------------------------------------------- #
# WHAT IS CACHED AND WHY IT IS THIS AND NOT THE REFERENCE.
#
# ``load_reference()`` peaks at a MEASURED 452 MiB and stays 322.7 MiB resident
# (results/scref_memory_2026-09-20.md). Memoising IT would pin that 322.7 MiB
# for the life of a process that deploys to a 512 MB container -- the fix that
# looks obvious and is wrong. What the Visium path actually wants out of that
# 452 MiB is this: a ``(G_shared, K)`` float64 matrix with K = 4, i.e. a few
# hundred kB. That is what is cached.
#
# KEY. ``(reference_fingerprint(), tuple(shared_genes))``. The fingerprint is
# (path, size, mtime_ns) of both cached files, so replacing the reference
# invalidates rather than serves a stale signature. ``shared_genes`` is part of
# the key because it HAS to be: ``_normalize_columns_to_unit_sum`` divides each
# cell by its total over the shared genes, so the signature for gene set A is
# not a row-slice of the signature for gene set B. A per-section cache is the
# largest correct unit, not a conservative choice.
#
# COPY-ON-READ. Callers get a fresh array every time. Handing out the cached
# object would make one caller's in-place edit everyone else's result, which is
# a correctness bug rather than an optimisation; ~640 kB a copy is not worth
# arguing about.
#
# A caller that supplies its own ``reference=`` (the held-out split in
# ``tmesim/validation/reconstruction_fidelity.py``) is NEVER cached: that frame
# is a subset chosen by the caller and the fingerprint cannot describe it.
_SIGNATURE_CACHE: "OrderedDict[Tuple, Tuple[np.ndarray, List[str]]]" = OrderedDict()

#: Enough for every built-in Visium section plus room; each entry is ~640 kB.
_SIGNATURE_CACHE_MAX = 8


def signature_cache_clear() -> None:
    """Drop every memoised signature. For tests and for an operator REPL."""
    _SIGNATURE_CACHE.clear()


def signature_cache_info() -> Dict:
    """``{'entries': n, 'max': m}`` -- what is currently memoised."""
    return {"entries": len(_SIGNATURE_CACHE), "max": _SIGNATURE_CACHE_MAX}


def build_signatures(
    shared_genes: Sequence[str],
    reference: "Tuple | None" = None,
) -> Tuple[np.ndarray, List[str]]:
    """Build per-cell-type mean-expression signatures over ``shared_genes``.

    Memoised when ``reference`` is omitted; see the block above for the key,
    why the reference frame itself is deliberately not the thing cached, and
    why the result is copied on the way out.

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
    if reference is not None:
        return _build_signatures_uncached(shared_genes, reference)

    key = (reference_fingerprint(), tuple(str(g) for g in shared_genes))
    hit = _SIGNATURE_CACHE.get(key)
    if hit is None:
        # The reference is loaded, consumed and dropped inside this call, so
        # its 452 MiB peak is not held alongside anything else.
        sig, names = _build_signatures_uncached(shared_genes, load_reference())
        _SIGNATURE_CACHE[key] = (sig, names)
        while len(_SIGNATURE_CACHE) > _SIGNATURE_CACHE_MAX:
            _SIGNATURE_CACHE.popitem(last=False)
        hit = (sig, names)
    else:
        _SIGNATURE_CACHE.move_to_end(key)
    sig, names = hit
    return sig.copy(), list(names)


def _build_signatures_uncached(
    shared_genes: Sequence[str],
    reference: Tuple,
) -> Tuple[np.ndarray, List[str]]:
    """:func:`build_signatures` with the memo removed. Unchanged arithmetic."""
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

    Raises
    ------
    DeconvolutionTooLargeForHost
        Before loading anything, if this host declares a memory limit smaller
        than the measured peak of this section. See
        :func:`check_visium_memory_budget`.
    """
    # REFUSE FIRST. `load_sample` alone was measured at +400 to +620 MiB, so
    # the check has to precede it or it is checking after the fact.
    check_visium_memory_budget(sample_id)

    sample = load_sample(sample_id)
    coords_um = np.asarray(sample["coords_um"], dtype=np.float64)
    genes = list(sample["genes"])
    expr = sample["expression"].tocsr()  # (S, G_visium)

    # Shared gene set = Visium symbols also present in the reference.
    #
    # This used to be `load_reference()` -- a MEASURED +1,444 MiB parse of a
    # 23,460 x 3,589 counts matrix -- to build a set of gene symbols and
    # nothing else, and it then held that frame resident while
    # build_signatures made three more copies of it. `reference_gene_symbols()`
    # reads only the symbol column: +30 MiB, and verified array_equal to the
    # index the full parse returns.
    ref_syms = set(reference_gene_symbols())
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

    # No `reference=` argument: build_signatures loads, consumes and frees the
    # reference inside the call (and memoises the small result), instead of
    # this frame staying alive across it.
    signatures, type_names = build_signatures(shared_genes)

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
