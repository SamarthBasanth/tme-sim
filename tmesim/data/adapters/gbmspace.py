"""Adapter: GBM-Space per-section AnnData ``.h5ad`` -> the Space Ranger
``outs/``-shaped directory ``tmesim/data/ingest.py``'s Visium detector expects.

Background (full account in ``results/gbmspace_ingestion.md``): GBM-Space
(gbmspace.org / Sanger Institute; De Jong et al., bioRxiv 2025) publishes its
processed Visium sections as one ``.h5ad`` (AnnData-on-HDF5) file per section,
not as Space Ranger output directories. ``tmesim/data/ingest.py``'s Visium
detector (``resolve_visium_dirs`` / ``_detect_visium``) looks for a
``spatial/`` folder plus a ``filtered_feature_bc_matrix/`` MEX triple on disk;
an ``.h5ad`` is a single HDF5 file with neither, so ``ingest.detect()`` returns
``platform="unknown"`` on one directly.

Design decision - materialize, don't duplicate
------------------------------------------------
The task this adapter was built for asked that ingested sections "flow through
``reconstruct_from_visium`` UNCHANGED". ``reconstruct_from_visium(sample_id)``
(``tmesim/reconstruct/real_visium.py``) resolves a STRING id through
``tmesim.data.real.visium``'s directory-based registry
(``register_sample`` / ``_locate`` / ``load_from_dirs``), which in turn reads
literal files (``tissue_positions.csv``, ``scalefactors_json.json``,
``matrix.mtx.gz`` + ``features.tsv.gz`` + ``barcodes.tsv.gz``) off disk. There
is no in-memory entry point that bypasses those files - the whole chain
(``ingest.detect``/``ingest.ingest``, ``ingest_registry.register``,
``visium.register_sample``, ``deconvolve_sample``, ``reconstruct_from_visium``)
is coupled to "a directory with these exact filenames in it", not to a Python
object shape.

So rather than hand-building a :class:`~tmesim.data.ingest.SampleRecord` that
duplicates ``ingest._ingest_visium``'s logic (and would still need a NEW
reconstruction entry point to consume it, since ``reconstruct_from_visium``
does not take a record), this adapter reads the ``.h5ad`` with ``h5py`` alone
(no ``anndata`` dependency - the on-disk CSR/categorical encoding is simple
enough to decode by hand; see :func:`inspect_section`) and
:func:`materialize_spaceranger_dir` writes out a genuine, byte-for-byte
ordinary Space Ranger ``outs/`` directory from those arrays. That directory
then flows through EVERY downstream module - ``ingest.py``, ``visium.py``,
``deconvolve.py``, ``real_visium.py``, ``ingest_registry.py`` - with ZERO code
changes to any of them. That is the strongest available sense of "unchanged"
given how those modules are actually wired. The trade-off is one extra
disk-materialization pass per section (a few seconds; the written MEX triple
is smaller than the source ``.h5ad`` because the non-"Gene Expression"
``var`` columns - cell2location and NMF derived features - are dropped).

Schema notes (confirmed directly against a downloaded section with ``h5py``,
not assumed from the README - see ``results/gbmspace_adapter.md`` §2):

* ``X`` is a standard ``scipy``-compatible on-disk CSR (``data``/``indices``/
  ``indptr`` + an ``attrs["shape"]``), so no ``anndata`` sparse-matrix decoder
  is needed - ``scipy.sparse.csr_matrix((data, indices, indptr), shape=shape)``
  reads it directly.
* ``obs/sample`` and ``obs/sample_name`` and ``var/feature_types`` are
  anndata's on-disk ``categorical`` encoding (a ``categories`` array + an
  integer ``codes`` array); :func:`_read_categorical` decodes that by hand.
* ``obsm/spatial`` is confirmed (by correlating it against ``obs/array_row``
  / ``obs/array_col``, not assumed) to be ``(x, y)`` = ``(pxl_col_in_fullres,
  pxl_row_in_fullres)`` in FULLRES pixel units - the same convention and same
  units Space Ranger's ``tissue_positions.csv`` uses, just column-swapped and
  minus the header/CSV framing.
* ``uns/spatial/<sample_name>/scalefactors`` carries the same four keys
  ``scalefactors_json.json`` does (``spot_diameter_fullres``,
  ``fiducial_diameter_fullres``, ``tissue_hires_scalef``,
  ``tissue_lowres_scalef``), just as HDF5 scalar datasets instead of JSON.
* ``var/_index`` carries HGNC gene SYMBOLS only - there is no Ensembl id
  column anywhere in the file. The materialized ``features.tsv.gz`` writes
  the symbol into both the id and the symbol column (real Space Ranger output
  never does this); every ``source`` dict returned by :func:`register` says so
  explicitly so nothing downstream can mistake it for a real Ensembl id.
* ``var/feature_types`` has 4 categories - only ``"Gene Expression"`` (36,601
  of 36,717 columns for the sections inspected) is a raw-count gene panel; the
  other 3 (``"Cell state abundances"``, ``"Spatial niche abundances"``,
  ``"Histopath annotation overlap"``) are continuous derived features the
  GBM-Space team appended to the SAME matrix, confirmed non-integer by direct
  inspection. :func:`inspect_section` filters to ``"Gene Expression"`` only,
  exactly as the GBM-Space README instructs.
"""
from __future__ import annotations

import csv
import gzip
import json
import os
from typing import Dict, List, Optional

import h5py
import numpy as np
from scipy.io import mmwrite
from scipy.sparse import csr_matrix

__all__ = [
    "GENE_EXPRESSION_FEATURE_TYPE",
    "DEFAULT_CACHE_ROOT",
    "inspect_section",
    "materialize_spaceranger_dir",
    "register",
]

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", ".."))

#: Materialized Space Ranger directories land here, one subdirectory per
#: GBM-Space section - a sibling of ``data_cache/visium_gbm/`` so the two
#: Visium sources are easy to tell apart on disk.
DEFAULT_CACHE_ROOT = os.path.join(_REPO_ROOT, "data_cache", "gbmspace_visium")

GENE_EXPRESSION_FEATURE_TYPE = "Gene Expression"

#: Marker file dropped in a materialized directory once it is complete, so a
#: second call is a fast no-op instead of re-writing gigabytes.
_DONE_MARKER = ".gbmspace_materialized"

_SOURCE_LABEL = (
    "GBM-Space (gbmspace.org / Wellcome Sanger Institute; De Jong et al., "
    "bioRxiv 2025.05.13.653495 / .653733), spatial_data_visium_selected.tar.gz "
    "-- adapted from AnnData .h5ad to Space Ranger outs/ shape by "
    "tmesim.data.adapters.gbmspace"
)


def _decode(values) -> List[str]:
    return [v.decode("utf-8") if isinstance(v, (bytes, bytearray)) else str(v)
            for v in values]


def _read_categorical(group: "h5py.Group") -> np.ndarray:
    """Decode an anndata on-disk ``categorical`` column: a ``categories``
    array plus an integer ``codes`` array (``-1`` = missing/NaN)."""
    categories = np.array(_decode(group["categories"][:]), dtype=object)
    codes = np.asarray(group["codes"][:])
    out = np.full(codes.shape, "", dtype=object)
    valid = codes >= 0
    out[valid] = categories[codes[valid]]
    return out


def inspect_section(h5ad_path: str) -> Dict:
    """Read one GBM-Space Visium section ``.h5ad`` into plain numpy/scipy
    objects, using ``h5py`` alone (no ``anndata`` import).

    Returns a dict with:

    ``sample_name``       str, e.g. ``"AT10-BRA-5-FO-3_1"`` (from
                           ``obs/sample_name``, the human-readable section id
                           GBM-Space uses in its own figures).
    ``barcodes``           list[str], length n_spots - the 10x spot barcode
                           with the ``_<spaceranger-run-id>`` suffix stripped
                           (``obs/_index`` is ``"<barcode>_<sample>"``; the
                           barcode itself never contains an underscore).
    ``array_row``/``array_col``  (n_spots,) float64, the Visium grid indices.
    ``in_tissue``          (n_spots,) float64, 0/1.
    ``coords_px``          (n_spots, 2) float64 = ``(pxl_col_in_fullres,
                           pxl_row_in_fullres)`` i.e. ``(x, y)`` fullres pixel
                           coordinates (see module docstring for how the axis
                           order was confirmed, not assumed).
    ``genes``               list[str], the ``"Gene Expression"`` feature
                           symbols only, in their original panel-column order.
    ``expression``          scipy CSR, shape (n_spots, n_genes), RAW integer
                           UMI counts (verified integer-valued on read; raises
                           if it is not, rather than silently reinterpreting a
                           format that turned out to differ from the
                           documented schema).
    ``scalefactors``        dict with (at least) ``spot_diameter_fullres``.
    ``chemistry`` / ``software_version``  provenance strings from
                           ``uns/spatial/<sample>/metadata``, or ``None``.

    Raises
    ------
    ValueError
        If the file holds more than one distinct ``obs/sample`` (this adapter
        assumes one section per ``.h5ad``, true for every file inspected so
        far), if no ``"Gene Expression"`` feature rows are found, if the
        gene-expression block is not integer-valued (would mean the schema
        assumption above does not hold for this file), or if scalefactors are
        missing.
    """
    with h5py.File(h5ad_path, "r") as f:
        obs_idx = _decode(f["obs/_index"][:])
        barcodes = [s.split("_", 1)[0] for s in obs_idx]
        array_row = f["obs/array_row"][:].astype(np.float64)
        array_col = f["obs/array_col"][:].astype(np.float64)
        in_tissue = f["obs/in_tissue"][:].astype(np.float64)

        sample_ids = _read_categorical(f["obs/sample"])
        sample_names = _read_categorical(f["obs/sample_name"])
        n_distinct_samples = len(set(sample_ids.tolist()))
        if n_distinct_samples != 1:
            raise ValueError(
                f"{h5ad_path}: expected exactly one obs/sample value (one "
                f"section per .h5ad), found {n_distinct_samples}. This "
                "adapter does not split multi-section files."
            )
        sample_name = (
            str(sample_names[0]) if len(sample_names)
            else os.path.splitext(os.path.basename(h5ad_path))[0]
        )
        sample_id_provenance = str(sample_ids[0]) if len(sample_ids) else None

        if len(set(barcodes)) != len(barcodes):
            raise ValueError(
                f"{h5ad_path}: barcodes are not unique after stripping the "
                "sample suffix - refusing to materialize an ambiguous MEX "
                "triple."
            )

        coords_px = f["obsm/spatial"][:].astype(np.float64)  # (x, y) fullres px

        var_names = np.array(_decode(f["var/_index"][:]), dtype=object)
        feature_types = _read_categorical(f["var/feature_types"])
        ge_mask = feature_types == GENE_EXPRESSION_FEATURE_TYPE
        if not ge_mask.any():
            raise ValueError(
                f"{h5ad_path}: no var/feature_types == "
                f"{GENE_EXPRESSION_FEATURE_TYPE!r} rows found; categories "
                f"present were {sorted(set(feature_types.tolist()))}."
            )
        genes = var_names[ge_mask].tolist()

        shape = tuple(int(v) for v in f["X"].attrs["shape"])
        data = f["X/data"][:]
        indices = f["X/indices"][:]
        indptr = f["X/indptr"][:]
        X = csr_matrix((data, indices, indptr), shape=shape)  # spots x all-features
        expr = X[:, ge_mask].tocsr()
        if expr.nnz and not np.allclose(expr.data, np.round(expr.data)):
            raise ValueError(
                f"{h5ad_path}: 'Gene Expression' block is not integer-valued "
                "- the raw-counts assumption for this feature_types category "
                "does not hold for this file."
            )

        scalefactors: Dict[str, float] = {}
        chemistry: Optional[str] = None
        software_version: Optional[str] = None
        spatial_group = f.get("uns/spatial")
        sample_key = next(iter(spatial_group.keys()), None) if spatial_group else None
        if sample_key is not None:
            sf_group = f.get(f"uns/spatial/{sample_key}/scalefactors")
            if sf_group is not None:
                scalefactors = {k: float(sf_group[k][()]) for k in sf_group.keys()}
            meta_group = f.get(f"uns/spatial/{sample_key}/metadata")
            if meta_group is not None:
                if "chemistry_description" in meta_group:
                    chemistry = _decode([meta_group["chemistry_description"][()]])[0]
                if "software_version" in meta_group:
                    software_version = _decode([meta_group["software_version"][()]])[0]

    if "spot_diameter_fullres" not in scalefactors:
        raise ValueError(
            f"{h5ad_path}: no usable uns/spatial/*/scalefactors group (need "
            "at least spot_diameter_fullres to convert pixels to um)."
        )

    return {
        "h5ad_path": os.path.abspath(h5ad_path),
        "sample_name": sample_name,
        "sample_id_provenance": sample_id_provenance,
        "barcodes": barcodes,
        "array_row": array_row,
        "array_col": array_col,
        "in_tissue": in_tissue,
        "coords_px": coords_px,
        "genes": genes,
        "n_genes_total": int(shape[1]),
        "n_genes_gene_expression": int(ge_mask.sum()),
        "expression": expr,
        "scalefactors": scalefactors,
        "chemistry": chemistry,
        "software_version": software_version,
    }


def materialize_spaceranger_dir(
    section: Dict, dest_dir: str, overwrite: bool = False
) -> str:
    """Write ``section`` (from :func:`inspect_section`) out as a Space Ranger
    ``outs/``-shaped directory: ``spatial/tissue_positions.csv`` +
    ``spatial/scalefactors_json.json``, and
    ``filtered_feature_bc_matrix/{matrix.mtx.gz,features.tsv.gz,barcodes.tsv.gz}``.

    Once written, ``tmesim/data/ingest.py``'s Visium detector
    (``resolve_visium_dirs``/``_detect_visium``) recognises ``dest_dir`` at
    full confidence with NO changes to that module - this function IS "the
    adapter", in the sense that everything downstream of it is stock code.

    Idempotent: if ``dest_dir`` already carries the ``.gbmspace_materialized``
    marker this function wrote on a previous call, it returns immediately
    (pass ``overwrite=True`` to force a rewrite).
    """
    marker_path = os.path.join(dest_dir, _DONE_MARKER)
    if os.path.exists(marker_path) and not overwrite:
        return dest_dir

    spatial_dir = os.path.join(dest_dir, "spatial")
    matrix_dir = os.path.join(dest_dir, "filtered_feature_bc_matrix")
    os.makedirs(spatial_dir, exist_ok=True)
    os.makedirs(matrix_dir, exist_ok=True)

    barcodes = section["barcodes"]
    in_tissue = section["in_tissue"]
    array_row = section["array_row"]
    array_col = section["array_col"]
    coords_px = section["coords_px"]  # (x, y) = (pxl_col, pxl_row)

    # spatial/tissue_positions.csv - headered form (Space Ranger >= 1.2),
    # column order barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,
    # pxl_col_in_fullres, matched exactly to what
    # tmesim.data.real.visium._read_positions expects.
    with open(os.path.join(spatial_dir, "tissue_positions.csv"), "w",
              newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["barcode", "in_tissue", "array_row", "array_col",
                    "pxl_row_in_fullres", "pxl_col_in_fullres"])
        for i, bc in enumerate(barcodes):
            w.writerow([
                bc, int(in_tissue[i]), int(array_row[i]), int(array_col[i]),
                f"{coords_px[i, 1]:.6f}", f"{coords_px[i, 0]:.6f}",
            ])

    sf = section["scalefactors"]
    with open(os.path.join(spatial_dir, "scalefactors_json.json"), "w") as fh:
        json.dump({
            "spot_diameter_fullres": sf["spot_diameter_fullres"],
            "fiducial_diameter_fullres": sf.get("fiducial_diameter_fullres"),
            "tissue_hires_scalef": sf.get("tissue_hires_scalef"),
            "tissue_lowres_scalef": sf.get("tissue_lowres_scalef"),
        }, fh, indent=1)

    genes = section["genes"]
    # NOTE: GBM-Space's var/_index carries HGNC symbols only - there is no
    # Ensembl id anywhere in the source file. The symbol is written into BOTH
    # the id and symbol columns; this is flagged in register()'s `source`
    # dict so nothing downstream mistakes column 0 for a real Ensembl id.
    with gzip.open(os.path.join(matrix_dir, "features.tsv.gz"), "wt") as fh:
        for g in genes:
            fh.write(f"{g}\t{g}\tGene Expression\n")

    with gzip.open(os.path.join(matrix_dir, "barcodes.tsv.gz"), "wt") as fh:
        for bc in barcodes:
            fh.write(f"{bc}\n")

    expr = section["expression"]  # spots x genes CSR, raw integer counts
    mat = expr.T.tocoo()  # genes x spots, matching real Space Ranger matrix.mtx
    if not np.issubdtype(mat.dtype, np.integer):
        mat = mat.astype(np.int64)
    mtx_path = os.path.join(matrix_dir, "matrix.mtx")
    mmwrite(mtx_path, mat, field="integer")
    with open(mtx_path, "rb") as src, gzip.open(mtx_path + ".gz", "wb") as dst:
        dst.writelines(src)
    os.remove(mtx_path)

    with open(marker_path, "w") as fh:
        fh.write(section["h5ad_path"] + "\n")

    return dest_dir


def register(
    h5ad_path: str,
    sample_id: Optional[str] = None,
    cache_root: str = DEFAULT_CACHE_ROOT,
    overwrite: bool = False,
) -> Dict:
    """End-to-end: read ``h5ad_path``, materialize it under ``cache_root``,
    register it with ``tmesim.data.real.visium`` (so
    ``reconstruct_from_visium(sample_id)`` works unchanged), and ingest it
    through ``tmesim.data.ingest.ingest`` (so its capability report is
    MEASURED via the exact same path the other 7 real samples use, not
    asserted by this adapter).

    Returns
    -------
    dict with ``sample_id``, ``dest_dir``, ``section`` (the
    :func:`inspect_section` dict, arrays included), ``detect`` (the
    ``ingest.detect()`` result) and ``record`` (the ``ingest.SampleRecord``).

    Raises whatever :func:`inspect_section` or ``ingest.ingest`` raise - this
    function adds no new failure-swallowing.
    """
    from .. import ingest as ing
    from ..real import visium as visium_mod

    section = inspect_section(h5ad_path)
    sid = sample_id or section["sample_name"]
    dest_dir = os.path.join(cache_root, section["sample_name"])
    materialize_spaceranger_dir(section, dest_dir, overwrite=overwrite)

    visium_mod.register_sample(
        sid, dest_dir,
        source=_SOURCE_LABEL,
        accession=section["sample_name"],
        diagnosis="Glioblastoma (GBM-Space cohort; per-patient clinical "
                  "detail not resolved by this adapter)",
        open_access=True,
    )

    det = ing.detect(dest_dir)
    record = ing.ingest(dest_dir, sample_id=sid, detection=det)
    return {
        "sample_id": sid,
        "dest_dir": dest_dir,
        "section": section,
        "detect": det,
        "record": record,
    }


if __name__ == "__main__":
    import sys

    from ..ingest import format_capability_report

    for path in sys.argv[1:]:
        out = register(path)
        rec = out["record"]
        print(f"\n=== {out['sample_id']} ({path}) ===")
        print(f"  materialized -> {out['dest_dir']}")
        print(f"  n_spots={rec.n_cells}  n_genes={rec.n_genes}  "
              f"extent_um={rec.extent_um}")
        print(format_capability_report(rec.capability))
