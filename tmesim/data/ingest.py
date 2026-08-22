"""Generic sample ingestion: point it at a directory, it works out what is there.

Why this module exists
----------------------
Before it, a section was a hardcoded row. ``app/server/cache.py`` enumerated
three Xenium cache directories and four Visium sample ids **by path**;
``data/real/singlecell_spatial.py`` required the three files to be named
``Xenium_GBM_cells.csv.gz`` / ``Xenium_GBM_cell_feature_matrix.tar.gz`` /
``Xenium_GBM_gene_panel.json``; ``data/real/visium.py`` would only load a sample
whose id was in its ``_SAMPLES`` dict. A new patient specimen could not be loaded
without editing code. That is the difference between a demo and a tool.

This module is the generic path. It takes ANY directory, decides what it is,
loads it through the EXISTING loaders wherever the layout matches (so a built-in
section ingests to bit-identical arrays — asserted in
``tests/test_ingest.py``), validates it, and — the part that matters more than
the parsing — reports what the specimen can and cannot support.

The capability report is the point
----------------------------------
``results/cell_state_impact.md`` established that the 319-gene flagship panel
fails BOTH cell-state axes: it carries no hypoxia-response transcript at all, and
its proliferation markers fire in 1.42% of tumour cells, 13x below the published
19-40% Ki-67 range and below the 5% floor at which a per-cell call is allowed.
A user who does not know that gets a number that means nothing. So the ingester
answers the question BEFORE anything is run:

    HYPOXIA           UNMEASURED - panel carries none of VEGFA, CA9, ADM, NDRG1,
                      SLC2A1. HYPOXIC is assigned to NO cell. That is not the
                      same statement as "computed and found to be zero".
    state_mode=markers UNAVAILABLE - proliferation detection 1.42% is below the
                      5.0% per-cell floor; the module would fall back to a
                      neighbourhood rule whose level moves 8.79x with the window.

Two levels of evidence, always distinguished:

* **panel presence** - is the gene on the panel at all? Cheap (read
  ``features.tsv.gz``), so ``detect``/``scan`` can answer it with no matrix read.
* **detection rate** - what fraction of cells actually express it? Needs the
  count matrix. :func:`ingest` always measures it; :func:`scan_report` measures
  it only when asked, and says ``UNMEASURED`` rather than guessing when it has
  not.

Layouts recognised
------------------
``xenium``
    A Xenium output bundle: a per-cell centroid table (``cells.csv.gz`` /
    ``*_cells.csv.gz`` / the repo's ``Xenium_GBM_cells.csv.gz``), a
    ``cell_feature_matrix`` MEX triple either extracted as a directory or still
    inside ``*cell_feature_matrix.tar.gz``, and optionally ``gene_panel.json``
    (from which provenance is read).
``visium``
    Space Ranger output: a ``spatial/`` directory carrying
    ``tissue_positions[_list].csv`` + ``scalefactors_json.json``, and a
    ``filtered_feature_bc_matrix`` (or ``raw_...``) MEX triple. Both may be at
    the root, one level down, or inside the sibling tarballs Space Ranger ships.
``mex``
    A bare MEX triple (``matrix.mtx`` + ``features.tsv``/``genes.tsv`` +
    ``barcodes.tsv``, gzipped or not) PLUS a coordinates table - any CSV/TSV with
    a barcode/cell-id column and an x/y pair. This is the escape hatch for a
    platform we have never seen.
``multi_sample``
    The directory is a container: it holds no sample of its own, but two or more
    of its children are samples. Reported with the child list rather than
    guessed at.

Nothing is guessed. :func:`detect` returns a confidence and an explicit list of
what is missing, and :func:`ingest` refuses with a reason rather than building a
broken :class:`~tmesim.core.microenvironment.Microenvironment`.
"""
from __future__ import annotations

import csv
import gzip
import io
import json
import os
import tarfile
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

__all__ = [
    "IngestError",
    "SampleRecord",
    "PLATFORMS",
    "detect",
    "ingest",
    "scan_report",
    "capability_report",
    "format_capability_report",
    "resolve_xenium_files",
    "resolve_visium_dirs",
    "panel_genes",
    "validate_record",
    "MIN_CELLS",
    "MIN_EXTENT_UM",
    "MAX_EXTENT_UM",
]


class IngestError(ValueError):
    """Raised when a directory cannot be ingested, with the reason why."""


#: Platforms :func:`detect` can name. ``"unknown"`` is also a possible verdict.
PLATFORMS: Tuple[str, ...] = ("xenium", "visium", "mex")

# --------------------------------------------------------------------------- #
# Validation thresholds. Every one of these is a REFUSAL boundary, so each is
# justified against the specimens we hold rather than picked round.
# --------------------------------------------------------------------------- #

#: Fewer cells than this is not a tissue section. The smallest section in this
#: repo is the 10x parent Visium sample at 3,468 in-tissue spots; the smallest
#: plausible ROI crop is a few hundred. 50 refuses only degenerate input.
MIN_CELLS = 50

#: A section spanning less than this in either axis cannot be gridded at the
#: production dx (12-30 um) - fewer than ~4 voxels across.
MIN_EXTENT_UM = 100.0

#: 10 cm. Larger than any histology section, so an extent above it means the
#: coordinates are not micrometers (pixels and nanometres both land here).
MAX_EXTENT_UM = 100_000.0

#: Below this a coordinate column is more likely to be millimetres than um: the
#: Xenium sections here span 5,100-24,000 um, Visium 6,500-8,900 um.
_SUSPECT_SMALL_EXTENT_UM = 500.0


# --------------------------------------------------------------------------- #
# Filename resolution (stdlib only - imported by the loaders themselves)
# --------------------------------------------------------------------------- #
_MEX_MATRIX = ("matrix.mtx.gz", "matrix.mtx")
_MEX_FEATURES = ("features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv")
_MEX_BARCODES = ("barcodes.tsv.gz", "barcodes.tsv")

_XENIUM_CELL_COLUMNS = ("x_centroid", "y_centroid")

#: Column-name candidates for a generic coordinates table, lowercased.
_X_COLUMNS = ("x_centroid", "x", "x_um", "x_micron", "x_location", "pxl_col_in_fullres",
              "centroid_x", "x_position", "col")
_Y_COLUMNS = ("y_centroid", "y", "y_um", "y_micron", "y_location", "pxl_row_in_fullres",
              "centroid_y", "y_position", "row")
_ID_COLUMNS = ("cell_id", "barcode", "cell", "id", "cellid", "spot_id", "index")


def _listdir(path: str) -> List[str]:
    try:
        return sorted(os.listdir(path))
    except OSError:
        return []


def _is_mex_dir(path: str) -> bool:
    if not os.path.isdir(path):
        return False
    names = set(_listdir(path))
    return (any(n in names for n in _MEX_MATRIX)
            and any(n in names for n in _MEX_FEATURES)
            and any(n in names for n in _MEX_BARCODES))


def _mex_missing(path: str) -> List[str]:
    """Which of the three MEX components are absent from ``path``."""
    names = set(_listdir(path))
    missing = []
    for label, options in (("matrix.mtx[.gz]", _MEX_MATRIX),
                           ("features.tsv[.gz]", _MEX_FEATURES),
                           ("barcodes.tsv[.gz]", _MEX_BARCODES)):
        if not any(n in names for n in options):
            missing.append(label)
    return missing


def resolve_xenium_files(path: str) -> Dict[str, Optional[str]]:
    """Locate the components of a Xenium bundle in ``path``.

    Returns ``{cells, matrix_dir, matrix_tar, panel, transcripts}``, each an
    absolute-or-as-given path or ``None``.

    RESOLUTION ORDER IS LOAD-BEARING. The repo-canonical basenames
    (``Xenium_GBM_cells.csv.gz`` etc.) are matched FIRST, so every cache already
    in ``data_cache/`` resolves to exactly the files
    ``data.real.singlecell_spatial.load`` used before this module existed and its
    output is unchanged. ``data_cache/gse324268_xenium`` holds BOTH spellings
    (a canonical copy beside the original ``GSM9572843_GBMtumor07_*``), and the
    canonical copy is the one the published numbers were measured from.
    """
    names = _listdir(path)
    lower = {n.lower(): n for n in names}

    def pick(*exact: str, suffixes: Sequence[str] = (),
             contains: Sequence[str] = (), reject: Sequence[str] = ()) -> Optional[str]:
        for want in exact:
            if want.lower() in lower:
                return os.path.join(path, lower[want.lower()])
        for n in names:
            nl = n.lower()
            if any(r in nl for r in reject):
                continue
            if any(nl.endswith(s) for s in suffixes) or any(c in nl for c in contains):
                if os.path.isfile(os.path.join(path, n)):
                    return os.path.join(path, n)
        return None

    cells = pick(
        "Xenium_GBM_cells.csv.gz", "cells.csv.gz", "cells.csv",
        suffixes=("_cells.csv.gz", "_cells.csv"),
    )

    matrix_dir = None
    for cand in ("cell_feature_matrix", "filtered_feature_bc_matrix",
                 "raw_feature_bc_matrix"):
        p = os.path.join(path, cand)
        if _is_mex_dir(p):
            matrix_dir = p
            break
    if matrix_dir is None and _is_mex_dir(path):
        matrix_dir = path

    matrix_tar = pick(
        "Xenium_GBM_cell_feature_matrix.tar.gz", "cell_feature_matrix.tar.gz",
        suffixes=("cell_feature_matrix.tar.gz", "cell_feature_matrix.tar"),
    )

    panel = pick(
        "Xenium_GBM_gene_panel.json", "gene_panel.json", "gene_panel.json.gz",
        suffixes=("_gene_panel.json", "_gene_panel.json.gz"),
    )

    transcripts = pick(
        "transcripts.csv.gz", "transcripts.parquet",
        suffixes=("_transcripts.csv.gz",),
    )

    return {
        "cells": cells,
        "matrix_dir": matrix_dir,
        "matrix_tar": matrix_tar,
        "panel": panel,
        "transcripts": transcripts,
    }


def resolve_visium_dirs(path: str) -> Dict[str, Optional[str]]:
    """Locate ``(spatial_dir, matrix_dir)`` of a Space Ranger output in ``path``.

    THE TWO DIRECTORIES MUST BE SIBLINGS. Resolving them independently is what
    turns a container into a phantom sample: ``data_cache/`` contains a
    ``spatial/`` three levels down under ``visium_gbm/`` and a MEX triple three
    levels down under ``gse324268_xenium/``, and an independent search happily
    pairs those two unrelated specimens into one "sample". So the search is for
    a single PARENT holding both - ``path`` itself first, then a child named
    ``outs`` (the Space Ranger convention). Nothing deeper: a directory that
    merely CONTAINS samples is a container, and :func:`detect` reports it as one
    rather than picking a child arbitrarily.

    Also returns ``spatial_tar`` / ``matrix_tar`` when only the tarballs are
    present, so :func:`detect` can say "present but not extracted".
    """
    spatial = None
    matrix = None
    roots = [path]
    for name in _listdir(path):
        if name.lower() in ("outs", "output"):
            p = os.path.join(path, name)
            if os.path.isdir(p):
                roots.append(p)
    for root in roots:
        sp = next((os.path.join(root, n) for n in _listdir(root)
                   if _is_spatial_dir(os.path.join(root, n))), None)
        mx = next((os.path.join(root, n) for n in _listdir(root)
                   if n != "spatial" and _is_mex_dir(os.path.join(root, n))), None)
        if sp is not None or mx is not None:
            spatial, matrix = sp, mx
        if sp is not None and mx is not None:
            break

    names = _listdir(path)
    spatial_tar = next(
        (os.path.join(path, n) for n in names if n.lower().endswith("spatial.tar.gz")),
        None)
    matrix_tar = next(
        (os.path.join(path, n) for n in names
         if n.lower().endswith("feature_bc_matrix.tar.gz")), None)

    return {"spatial_dir": spatial, "matrix_dir": matrix,
            "spatial_tar": spatial_tar, "matrix_tar": matrix_tar}


def _dirs_at(root: str, depth: int) -> List[str]:
    """Directories exactly ``depth`` levels under ``root`` (``depth=0`` -> root)."""
    if depth == 0:
        return [root] if os.path.isdir(root) else []
    out: List[str] = []
    for parent in _dirs_at(root, depth - 1):
        for name in _listdir(parent):
            p = os.path.join(parent, name)
            if os.path.isdir(p):
                out.append(p)
    return out


def _is_spatial_dir(path: str) -> bool:
    if not os.path.isdir(path):
        return False
    names = {n.lower() for n in _listdir(path)}
    has_pos = any(n.startswith("tissue_positions") and ".csv" in n for n in names)
    return has_pos and "scalefactors_json.json" in names


# --------------------------------------------------------------------------- #
# detect()
# --------------------------------------------------------------------------- #
def detect(path: str) -> Dict:
    """Work out what kind of spatial sample ``path`` holds.

    Returns
    -------
    dict with keys
        ``path``        the directory inspected (absolute).
        ``platform``    one of :data:`PLATFORMS`, ``"multi_sample"`` or
                        ``"unknown"``.
        ``layout``      the specific on-disk shape, e.g.
                        ``"xenium_bundle_extracted"``, ``"xenium_bundle_tar"``,
                        ``"spaceranger_outs"``, ``"spaceranger_tar"``,
                        ``"mex_plus_coordinates"``, ``"container"``, ``None``.
        ``confidence``  0.0-1.0. 1.0 means every component of that layout is
                        present and named as expected; it drops for each
                        component that had to be inferred or is missing.
        ``files``       resolved component paths (``None`` where absent).
        ``issues``      what is missing or ambiguous, in plain words. NOT
                        empty just because confidence is high.
        ``notes``       non-blocking observations.
        ``contained_samples``  child directories that are themselves samples.
        ``ingestable``  whether :func:`ingest` would be able to proceed.
    """
    path = os.path.abspath(os.path.expanduser(str(path)))
    out: Dict = {
        "path": path, "platform": "unknown", "layout": None, "confidence": 0.0,
        "files": {}, "issues": [], "notes": [], "contained_samples": [],
        "ingestable": False,
    }
    if not os.path.isdir(path):
        out["issues"].append(
            f"{path} is not a directory"
            + ("" if os.path.exists(path) else " and does not exist"))
        return out

    xen = _detect_xenium(path)
    vis = _detect_visium(path)
    mex = _detect_mex(path)

    best = max((xen, vis, mex), key=lambda d: d["confidence"])
    if best["confidence"] > 0.0:
        out.update(best)
        out["path"] = path
    else:
        # Nothing here. Is it a container of samples?
        for name in _listdir(path):
            child = os.path.join(path, name)
            if not os.path.isdir(child):
                continue
            sub = max((_detect_xenium(child), _detect_visium(child),
                       _detect_mex(child)), key=lambda d: d["confidence"])
            if sub["confidence"] >= 0.5:
                out["contained_samples"].append(
                    {"path": child, "name": name, "platform": sub["platform"],
                     "layout": sub["layout"], "confidence": sub["confidence"]})
        if len(out["contained_samples"]) >= 1:
            out["platform"] = "multi_sample"
            out["layout"] = "container"
            out["confidence"] = 0.9
            out["notes"].append(
                f"{path} holds no sample of its own but "
                f"{len(out['contained_samples'])} of its children are samples; "
                "ingest each child directory separately."
            )
        else:
            out["issues"].append(
                "no recognised layout. Expected one of: a Xenium bundle "
                "(cells.csv.gz + cell_feature_matrix[.tar.gz] + gene_panel.json), "
                "a Space Ranger output (spatial/ with tissue_positions + "
                "scalefactors_json.json, plus filtered_feature_bc_matrix/), or a "
                "MEX triple (matrix.mtx + features.tsv + barcodes.tsv) beside a "
                "coordinates table."
            )
        return out

    # A recognised sample may ALSO contain other samples (data_cache/visium_gbm
    # is the 10x parent section at the root with GBM1/2/3 beside it). Say so.
    own = {out["files"].get("spatial_dir"), out["files"].get("matrix_dir")}
    for name in _listdir(path):
        child = os.path.join(path, name)
        if not os.path.isdir(child) or child in own:
            continue
        sub = max((_detect_xenium(child), _detect_visium(child)),
                  key=lambda d: d["confidence"])
        if sub["confidence"] >= 0.5:
            out["contained_samples"].append(
                {"path": child, "name": name, "platform": sub["platform"],
                 "layout": sub["layout"], "confidence": sub["confidence"]})
    if out["contained_samples"]:
        out["notes"].append(
            f"this directory is BOTH a sample and a container: "
            f"{len(out['contained_samples'])} child directories are separate "
            "samples. Only the root sample is ingested from this path."
        )
    return out


def _detect_xenium(path: str) -> Dict:
    files = resolve_xenium_files(path)
    issues: List[str] = []
    notes: List[str] = []
    if files["cells"] is None:
        return _empty_detection(path)

    # A Space Ranger bundle has no per-cell centroid table, so a cells file is
    # already strong evidence. Confirm by looking for the centroid columns.
    cols = _peek_columns(files["cells"])
    if cols is not None and not all(c in cols for c in _XENIUM_CELL_COLUMNS):
        # Might still be a generic per-cell table -> let the MEX branch have it.
        return _empty_detection(path)

    confidence = 0.6
    layout = None
    if files["matrix_dir"] is not None:
        confidence += 0.25
        layout = "xenium_bundle_extracted"
        missing = _mex_missing(files["matrix_dir"])
        if missing:  # cannot happen via _is_mex_dir, kept as a guard
            issues.append(f"MEX triple incomplete: missing {missing}")
    elif files["matrix_tar"] is not None:
        confidence += 0.2
        layout = "xenium_bundle_tar"
        notes.append(
            f"cell_feature_matrix is still packed in "
            f"{os.path.basename(files['matrix_tar'])}; ingest() extracts it "
            "in place on first load (one-off, ~15 s for a 240 MB tarball)."
        )
    else:
        # A matrix DIRECTORY that exists but is short a file is by far the
        # commonest real failure (an interrupted download, a partial copy), and
        # "no cell-feature matrix" is actively misleading about it - the
        # directory is right there. Name the missing FILE instead.
        partial = None
        for cand in ("cell_feature_matrix", "filtered_feature_bc_matrix",
                     "raw_feature_bc_matrix"):
            cand_path = os.path.join(path, cand)
            if os.path.isdir(cand_path):
                gap = _mex_missing(cand_path)
                if gap:
                    partial = (cand, gap)
                    break
        if partial is not None:
            issues.append(
                f"{partial[0]}/ is INCOMPLETE: missing {', '.join(partial[1])}. "
                "A MEX triple needs all three of matrix.mtx[.gz], "
                "features.tsv[.gz] and barcodes.tsv[.gz]; with any one absent "
                "the counts cannot be aligned to genes and cells. Re-download "
                "or re-extract that directory."
            )
        else:
            issues.append(
                "no cell-feature matrix: neither an extracted "
                "cell_feature_matrix/ MEX triple nor a "
                "*cell_feature_matrix.tar.gz. Expression is REQUIRED - without "
                "it no phenotyping, no ECM field and no capability report are "
                "possible."
            )
        layout = "xenium_cells_only"

    if files["panel"] is not None:
        confidence += 0.15
    else:
        issues.append(
            "no gene_panel.json: the panel identity (and therefore the dataset "
            "provenance) cannot be read, so `source` will say UNREGISTERED. "
            "The gene list still comes from features.tsv.gz, so the capability "
            "report is unaffected."
        )

    if os.path.basename(files["cells"]).startswith("Xenium_GBM_"):
        notes.append(
            "resolved to the repo-canonical Xenium_GBM_* basenames, which is "
            "what data.real.singlecell_spatial.load used before this module "
            "existed - so this section ingests bit-identically."
        )

    return {
        "path": path, "platform": "xenium", "layout": layout,
        "confidence": round(min(confidence, 1.0), 3), "files": files,
        "issues": issues, "notes": notes, "contained_samples": [],
        "ingestable": files["matrix_dir"] is not None or files["matrix_tar"] is not None,
    }


def _detect_visium(path: str) -> Dict:
    dirs = resolve_visium_dirs(path)
    issues: List[str] = []
    notes: List[str] = []
    has_spatial = dirs["spatial_dir"] is not None
    has_matrix = dirs["matrix_dir"] is not None
    if not has_spatial and dirs["spatial_tar"] is None:
        return _empty_detection(path)

    confidence = 0.0
    layout = "spaceranger_outs"
    if has_spatial:
        confidence += 0.5
        sp = dirs["spatial_dir"]
        names = {n.lower() for n in _listdir(sp)}
        if not any(n.startswith("tissue_positions") for n in names):
            issues.append(f"{sp}: no tissue_positions[_list].csv")
        if "scalefactors_json.json" not in names:
            issues.append(
                f"{sp}: no scalefactors_json.json - spot_diameter_fullres is "
                "the only way to convert pixels to micrometres, so without it "
                "the coordinates cannot be put in um at all."
            )
    else:
        layout = "spaceranger_tar"
        confidence += 0.3
        issues.append(
            f"spatial/ is still packed in "
            f"{os.path.basename(dirs['spatial_tar'])}; extract it before "
            "ingesting (this module does not unpack Visium tarballs, because "
            "Space Ranger tars are not consistently rooted)."
        )

    if has_matrix:
        confidence += 0.5
    elif dirs["matrix_tar"] is not None:
        layout = "spaceranger_tar"
        confidence += 0.2
        issues.append(
            f"the feature-barcode matrix is still packed in "
            f"{os.path.basename(dirs['matrix_tar'])}; extract it before "
            "ingesting."
        )
    else:
        issues.append(
            "no feature-barcode matrix: neither filtered_feature_bc_matrix/ "
            "(a MEX triple) nor its tarball."
        )

    return {
        "path": path, "platform": "visium", "layout": layout,
        "confidence": round(min(confidence, 1.0), 3), "files": dirs,
        "issues": issues, "notes": notes, "contained_samples": [],
        "ingestable": bool(has_spatial and has_matrix),
    }


def _detect_mex(path: str) -> Dict:
    """A bare MEX triple plus SOME table that carries coordinates."""
    matrix_dir = None
    if _is_mex_dir(path):
        matrix_dir = path
    else:
        for name in _listdir(path):
            p = os.path.join(path, name)
            if _is_mex_dir(p):
                matrix_dir = p
                break
    if matrix_dir is None:
        return _empty_detection(path)

    coords, coord_cols = _find_coordinates_file(path)
    issues: List[str] = []
    notes: List[str] = []
    confidence = 0.5
    if coords is not None:
        confidence += 0.35
        notes.append(
            f"coordinates read from {os.path.basename(coords)} using columns "
            f"x={coord_cols[0]!r}, y={coord_cols[1]!r}, id={coord_cols[2]!r}."
        )
        if coord_cols[2] is None:
            confidence -= 0.15
            issues.append(
                f"{os.path.basename(coords)} has no barcode/cell-id column, so "
                "its rows are assumed to be in the SAME ORDER as barcodes.tsv. "
                "That assumption is not verifiable from the bytes."
            )
    else:
        issues.append(
            "MEX triple found but no coordinates table. A spatial sample needs "
            "per-cell (x, y): supply a CSV/TSV beside the matrix with an "
            f"x column (one of {list(_X_COLUMNS[:4])}), a y column and ideally "
            "a barcode/cell_id column."
        )

    return {
        "path": path, "platform": "mex", "layout": "mex_plus_coordinates",
        "confidence": round(min(confidence, 1.0), 3),
        "files": {"matrix_dir": matrix_dir, "coordinates": coords,
                  "coordinate_columns": list(coord_cols) if coords else None},
        "issues": issues, "notes": notes, "contained_samples": [],
        "ingestable": coords is not None,
    }


def _empty_detection(path: str) -> Dict:
    return {"path": path, "platform": "unknown", "layout": None,
            "confidence": 0.0, "files": {}, "issues": [], "notes": [],
            "contained_samples": [], "ingestable": False}


# --------------------------------------------------------------------------- #
# Small readers (stdlib; used for cheap peeks and for the generic MEX path)
# --------------------------------------------------------------------------- #
def _open_text(path: str):
    if str(path).endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8",
                                errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def _peek_columns(path: str) -> Optional[List[str]]:
    """Header of a (possibly gzipped) delimited file, or ``None``."""
    try:
        with _open_text(path) as fh:
            first = fh.readline()
    except OSError:
        return None
    if not first:
        return None
    delim = "\t" if first.count("\t") > first.count(",") else ","
    return [c.strip().strip('"') for c in first.rstrip("\r\n").split(delim)]


def _find_coordinates_file(path: str) -> Tuple[Optional[str], Tuple]:
    """First delimited file under ``path`` (depth <= 1) with an x/y column pair."""
    candidates: List[str] = []
    for name in _listdir(path):
        p = os.path.join(path, name)
        if os.path.isfile(p) and _looks_delimited(name):
            candidates.append(p)
    for name in _listdir(path):
        p = os.path.join(path, name)
        if os.path.isdir(p):
            for sub in _listdir(p):
                sp = os.path.join(p, sub)
                if os.path.isfile(sp) and _looks_delimited(sub):
                    candidates.append(sp)
    for cand in candidates:
        if os.path.basename(cand).startswith(("features", "genes", "barcodes",
                                              "matrix")):
            continue
        cols = _peek_columns(cand)
        if not cols:
            continue
        low = [c.lower() for c in cols]
        xcol = next((cols[low.index(c)] for c in _X_COLUMNS if c in low), None)
        ycol = next((cols[low.index(c)] for c in _Y_COLUMNS if c in low), None)
        idcol = next((cols[low.index(c)] for c in _ID_COLUMNS if c in low), None)
        if xcol and ycol and xcol != ycol:
            return cand, (xcol, ycol, idcol)
    return None, (None, None, None)


def _looks_delimited(name: str) -> bool:
    n = name.lower()
    for suffix in (".gz",):
        if n.endswith(suffix):
            n = n[: -len(suffix)]
    return n.endswith((".csv", ".tsv", ".txt"))


def panel_genes(path: str, detection: Optional[Dict] = None) -> List[str]:
    """Gene symbols on this sample's panel, WITHOUT reading the count matrix.

    Reads ``features.tsv.gz`` only. For a Xenium bundle this returns the
    ``Gene Expression`` features (negative-control and codeword rows dropped,
    matching ``singlecell_spatial.load(genes_only=True)``); for Visium it
    returns every feature symbol. Empty list when the matrix directory is not
    extracted - the panel then genuinely is not knowable cheaply.
    """
    det = detection if detection is not None else detect(path)
    matrix_dir = (det.get("files") or {}).get("matrix_dir")
    if not matrix_dir:
        return []
    features = None
    for name in _MEX_FEATURES:
        p = os.path.join(matrix_dir, name)
        if os.path.exists(p):
            features = p
            break
    if features is None:
        return []
    genes: List[str] = []
    genes_only = det.get("platform") == "xenium"
    try:
        with _open_text(features) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                if genes_only and len(parts) >= 3:
                    if parts[2].strip() != "Gene Expression":
                        continue
                genes.append(parts[1])
    except OSError:
        return []
    return genes


# --------------------------------------------------------------------------- #
# The record
# --------------------------------------------------------------------------- #
@dataclass
class SampleRecord:
    """One ingested specimen, in the shape the existing loaders return.

    The four array fields are exactly what
    ``data.real.singlecell_spatial.load`` produces for a Xenium section, so
    downstream reconstruction code needs no changes:

    ``cell_xy_um``   (n, 2) float64 centroids in micrometres.
    ``expression``   (n, n_genes) scipy CSR of raw counts, rows aligned to
                     ``cell_xy_um`` and columns to ``genes``.
    ``genes``        (n_genes,) str ndarray of symbols.
    ``cell_ids``     (n,) str ndarray, aligned to the rows.

    ``native`` is the untouched dict the platform loader returned (so a Visium
    record still carries ``coords_um``, ``umpp``, ``barcodes`` and the rest under
    their original names). ``capability`` is the report from
    :func:`capability_report`; ``validation`` records every check that was run.
    """

    sample_id: str
    platform: str
    path: str
    cell_xy_um: np.ndarray
    expression: object
    genes: np.ndarray
    cell_ids: np.ndarray
    coordinate_units: str
    source: Dict
    capability: Dict = field(default_factory=dict)
    detection: Dict = field(default_factory=dict)
    validation: Dict = field(default_factory=dict)
    native: Dict = field(default_factory=dict)

    @property
    def n_cells(self) -> int:
        return int(self.cell_xy_um.shape[0])

    @property
    def n_genes(self) -> int:
        return int(np.asarray(self.genes).size)

    @property
    def extent_um(self) -> Tuple[float, float]:
        xy = self.cell_xy_um
        return (float(xy[:, 0].max() - xy[:, 0].min()),
                float(xy[:, 1].max() - xy[:, 1].min()))

    def as_loader_dict(self) -> Dict:
        """The record as the Xenium loader's dict, for code that expects it."""
        out = dict(self.native)
        out.update({
            "cell_xy_um": self.cell_xy_um,
            "expression": self.expression,
            "genes": self.genes,
            "cell_ids": self.cell_ids,
            "coordinate_units": self.coordinate_units,
            "source": self.source,
        })
        return out

    def summary(self) -> Dict:
        """JSON-safe summary (no arrays), for an API response."""
        ex, ey = self.extent_um
        return {
            "sample_id": self.sample_id,
            "platform": self.platform,
            "path": self.path,
            "n_cells": self.n_cells,
            "n_genes": self.n_genes,
            "extent_um": [ex, ey],
            "coordinate_units": self.coordinate_units,
            "source": self.source,
            "detection": {k: v for k, v in self.detection.items()
                          if k != "files"},
            "validation": self.validation,
            "capability": self.capability,
        }


# --------------------------------------------------------------------------- #
# ingest()
# --------------------------------------------------------------------------- #
def ingest(
    path: str,
    sample_id: Optional[str] = None,
    detection: Optional[Dict] = None,
    measure_detection_rates: bool = True,
    strict: bool = True,
) -> SampleRecord:
    """Load ANY recognised spatial sample directory into a :class:`SampleRecord`.

    Parameters
    ----------
    path : str
        The sample directory. :func:`detect` decides what it is.
    sample_id : str, optional
        Identifier for the record. Defaults to the directory basename.
    detection : dict, optional
        A :func:`detect` result to reuse (avoids a second directory scan).
    measure_detection_rates : bool
        Measure per-gene detection RATES from the loaded count matrix for the
        capability report. This is what turns "MKI67 is on the panel" into
        "MKI67 fires in 1.42% of cells, below the 5% floor". Cheap once the
        matrix is in memory; set False only to skip it deliberately.
    strict : bool
        Raise :class:`IngestError` on a validation failure (default). With
        ``strict=False`` the failures are recorded in ``record.validation`` and
        the record is returned anyway - for diagnosing a bad specimen, never for
        running one.

    Raises
    ------
    IngestError
        With the specific reason: unrecognised layout, a missing required
        component, or a failed validation check. Never returns a record that
        would produce a broken ``Microenvironment``.
    """
    path = os.path.abspath(os.path.expanduser(str(path)))
    det = detection if detection is not None else detect(path)
    sid = sample_id or os.path.basename(path.rstrip(os.sep)) or "sample"

    if det["platform"] == "multi_sample":
        children = ", ".join(c["name"] for c in det["contained_samples"])
        raise IngestError(
            f"{path} is a CONTAINER, not a sample: it holds no sample of its "
            f"own and {len(det['contained_samples'])} of its children are "
            f"samples ({children}). Ingest one of them."
        )
    if det["platform"] == "unknown" or not det["ingestable"]:
        reason = "; ".join(det["issues"]) or "no recognised layout"
        raise IngestError(f"cannot ingest {path}: {reason}")

    if det["platform"] == "xenium":
        record = _ingest_xenium(path, sid, det)
    elif det["platform"] == "visium":
        record = _ingest_visium(path, sid, det)
    else:
        record = _ingest_mex(path, sid, det)

    record.detection = {k: v for k, v in det.items() if k != "files"}
    record.detection["files"] = {
        k: v for k, v in (det.get("files") or {}).items() if isinstance(v, (str, list))
    }
    record.validation = validate_record(record)
    if strict and not record.validation["ok"]:
        raise IngestError(
            f"{path} loaded but FAILED validation, so it would produce a "
            "broken Microenvironment:\n  - "
            + "\n  - ".join(record.validation["errors"])
        )
    record.capability = capability_report(
        record.genes,
        expression=record.expression if measure_detection_rates else None,
        platform=record.platform,
        sample_id=record.sample_id,
    )
    return record


def _ingest_xenium(path: str, sample_id: str, det: Dict) -> SampleRecord:
    """Reuse ``data.real.singlecell_spatial.load`` verbatim.

    That loader now resolves its filenames through :func:`resolve_xenium_files`,
    so it reads a stock 10x bundle as happily as a repo cache - and for a repo
    cache it resolves to exactly the paths it always used, which is why the
    generic path reproduces the specific one bit-for-bit.
    """
    from .real.singlecell_spatial import load as load_xenium

    data = load_xenium(path=path, genes_only=True, dense=False)
    return SampleRecord(
        sample_id=sample_id,
        platform="xenium",
        path=path,
        cell_xy_um=data["cell_xy_um"],
        expression=data["expression"],
        genes=np.asarray(data["genes"], dtype=str),
        cell_ids=np.asarray(data["cell_ids"], dtype=str),
        coordinate_units=data["coordinate_units"],
        source=dict(data["source"], ingest_path=path, ingest_platform="xenium"),
        native=data,
    )


def _ingest_visium(path: str, sample_id: str, det: Dict) -> SampleRecord:
    """Reuse ``data.real.visium.load_from_dirs`` verbatim.

    ``load_sample(sample_id)`` is that same function plus a registry lookup for
    the two directories, so a registered sample ingests to identical arrays.

    ``coordinate_units`` is ``"micrometer"``: the loader converts pixels using
    ``55 um / spot_diameter_fullres`` from ``scalefactors_json.json`` and shifts
    the origin to (0, 0). ``cell_ids`` are SPOT barcodes - a Visium row is a
    ~55 um spot pooling several cells, not a cell, and nothing downstream may
    read it as one.
    """
    from .real.visium import load_from_dirs

    dirs = det["files"]
    data = load_from_dirs(dirs["spatial_dir"], dirs["matrix_dir"],
                          sample_id=sample_id, marker_report=True)
    source = {
        "dataset": f"Space Ranger output ingested from {path}",
        "accession": "unregistered (ingested by path)",
        "sample": sample_id,
        "assay": "Visium spatial gene expression (spot-level)",
        "panel": f"whole transcriptome ({data['n_genes']} features, 55 um spots)",
        "license": "unknown - verify before redistribution",
        "landing_page": None,
        "has_outcome_label": False,
        "cache_dir": path,
        "n_panel_genes": int(data["n_genes"]),
        "spatial_dir": dirs["spatial_dir"],
        "matrix_dir": dirs["matrix_dir"],
        "umpp": float(data["umpp"]),
        "ingest_path": path,
        "ingest_platform": "visium",
        "resolution_caveat": (
            "SPOT-level: each row is a ~55 um spot pooling several cells, NOT a "
            "cell. Per-cell routes (state_mode='markers', "
            "oxygen_mode='markers', sensitivity_mode='determinant') are "
            "unavailable on this sample by construction, not by panel content."
        ),
    }
    return SampleRecord(
        sample_id=sample_id,
        platform="visium",
        path=path,
        cell_xy_um=data["coords_um"],
        expression=data["expression"],
        genes=np.asarray(data["genes"], dtype=str),
        cell_ids=np.asarray(data["barcodes"], dtype=str),
        coordinate_units="micrometer",
        source=source,
        native=data,
    )


def _ingest_mex(path: str, sample_id: str, det: Dict) -> SampleRecord:
    """Load a bare MEX triple plus a coordinates table.

    The escape hatch for a platform with no dedicated loader. The MEX read is
    ``singlecell_spatial._read_mex`` (the same function the Xenium path uses);
    only the coordinate join is new, and every assumption it makes is recorded
    in ``source``.
    """
    from scipy import sparse

    from .real.singlecell_spatial import _read_mex

    matrix_dir = det["files"]["matrix_dir"]
    coord_path = det["files"]["coordinates"]
    xcol, ycol, idcol = det["files"]["coordinate_columns"]

    mat, features, barcodes = _read_mex(matrix_dir)  # features x cells
    bc = np.asarray(barcodes["cell_id"].to_numpy(), dtype=str)

    rows = _read_coordinate_table(coord_path, xcol, ycol, idcol)
    assumptions: List[str] = []
    if idcol is None:
        if len(rows["x"]) != bc.size:
            raise IngestError(
                f"{coord_path} has {len(rows['x'])} rows but the matrix has "
                f"{bc.size} columns, and the table has no barcode/cell_id "
                "column to join on. Add an id column."
            )
        order = np.arange(bc.size)
        assumptions.append(
            "coordinates joined BY ROW ORDER (the table carries no id column). "
            "If the table is not in barcodes.tsv order, every cell is at the "
            "wrong position and nothing downstream can detect it."
        )
    else:
        index = {cid: i for i, cid in enumerate(rows["id"])}
        missing = [b for b in bc if b not in index]
        if missing:
            raise IngestError(
                f"{len(missing)} of {bc.size} matrix barcodes are absent from "
                f"{os.path.basename(coord_path)} (e.g. {missing[:3]}). Every "
                "cell in the matrix must have a coordinate; a partial join "
                "would silently drop cells."
            )
        order = np.array([index[b] for b in bc], dtype=int)

    xy = np.column_stack([np.asarray(rows["x"], dtype=np.float64)[order],
                          np.asarray(rows["y"], dtype=np.float64)[order]])

    keep = (features["type"] == "Gene Expression").to_numpy() \
        if "type" in features and features["type"].notna().any() \
        else np.ones(len(features), dtype=bool)
    if not keep.any():
        keep = np.ones(len(features), dtype=bool)
    mat = mat[keep, :]
    genes = features.loc[keep, "symbol"].to_numpy().astype(str)
    expr = sparse.csr_matrix(mat.T)

    source = {
        "dataset": f"generic MEX + coordinates ingested from {path}",
        "accession": "unregistered (ingested by path)",
        "sample": sample_id,
        "assay": "unknown (generic MEX triple)",
        "panel": f"{genes.size} features from features.tsv",
        "license": "unknown - verify before redistribution",
        "landing_page": None,
        "has_outcome_label": False,
        "cache_dir": path,
        "n_panel_genes": int(genes.size),
        "matrix_dir": matrix_dir,
        "coordinates_file": coord_path,
        "coordinate_columns": {"x": xcol, "y": ycol, "id": idcol},
        "ingest_path": path,
        "ingest_platform": "mex",
        "assumptions": assumptions or [
            "coordinate UNITS are assumed to be micrometres: a generic table "
            "carries no unit metadata. validate_record() range-checks the "
            "extent, which catches pixels and nanometres but not millimetres "
            "on a small crop."
        ],
    }
    return SampleRecord(
        sample_id=sample_id, platform="mex", path=path,
        cell_xy_um=xy, expression=expr, genes=genes, cell_ids=bc,
        coordinate_units="micrometer (ASSUMED - see source.assumptions)",
        source=source,
        native={"matrix_dir": matrix_dir, "coordinates_file": coord_path,
                "features": features, "barcodes": bc},
    )


def _read_coordinate_table(path: str, xcol: str, ycol: str,
                           idcol: Optional[str]) -> Dict[str, List]:
    out: Dict[str, List] = {"x": [], "y": [], "id": []}
    with _open_text(path) as fh:
        header = fh.readline()
    delim = "\t" if header.count("\t") > header.count(",") else ","
    with _open_text(path) as fh:
        for row in csv.DictReader(fh, delimiter=delim):
            try:
                x = float(row[xcol])
                y = float(row[ycol])
            except (TypeError, ValueError, KeyError):
                continue
            out["x"].append(x)
            out["y"].append(y)
            if idcol is not None:
                out["id"].append(str(row.get(idcol, "")).strip())
    return out


# --------------------------------------------------------------------------- #
# Validation
# --------------------------------------------------------------------------- #
def validate_record(record: SampleRecord) -> Dict:
    """Check a loaded record is usable, and say exactly what is wrong if not.

    Every check is a property of the loaded bytes, and every failure names the
    measured value beside the bound it violated. Returns
    ``{ok, errors[], warnings[], measured{}}``; :func:`ingest` raises on
    ``ok == False`` rather than handing back something that would build a broken
    ``Microenvironment``.
    """
    from scipy import sparse

    errors: List[str] = []
    warnings: List[str] = []
    measured: Dict = {}

    xy = np.asarray(record.cell_xy_um, dtype=np.float64)
    expr = record.expression
    genes = np.asarray(record.genes, dtype=str)
    ids = np.asarray(record.cell_ids, dtype=str)

    # --- shape + alignment ------------------------------------------------- #
    if xy.ndim != 2 or xy.shape[1] != 2:
        errors.append(f"cell_xy_um must be (n, 2); got shape {xy.shape}")
        return {"ok": False, "errors": errors, "warnings": warnings,
                "measured": measured}
    n = int(xy.shape[0])
    measured["n_cells"] = n
    measured["n_genes"] = int(genes.size)

    if expr.shape[0] != n:
        errors.append(
            f"expression has {expr.shape[0]} rows but there are {n} cells; the "
            "matrix is not aligned to the coordinates, so every cell would "
            "carry another cell's transcripts."
        )
    if expr.shape[1] != genes.size:
        errors.append(
            f"expression has {expr.shape[1]} columns but {genes.size} gene "
            "symbols; marker lookup would read the wrong gene."
        )
    if ids.size != n:
        errors.append(f"{ids.size} cell_ids for {n} cells")

    # --- cell count --------------------------------------------------------- #
    if n < MIN_CELLS:
        errors.append(
            f"only {n} cells (floor {MIN_CELLS}). Below this there is no "
            "tissue to reconstruct: phenotyping thresholds, the vascular field "
            "and every band fraction would be dominated by single cells."
        )

    # --- coordinates -------------------------------------------------------- #
    if not np.isfinite(xy).all():
        bad = int((~np.isfinite(xy)).any(axis=1).sum())
        errors.append(f"{bad} cells have non-finite coordinates")
    else:
        ex = float(xy[:, 0].max() - xy[:, 0].min())
        ey = float(xy[:, 1].max() - xy[:, 1].min())
        measured["extent_um"] = [ex, ey]
        if ex < MIN_EXTENT_UM or ey < MIN_EXTENT_UM:
            errors.append(
                f"extent {ex:.1f} x {ey:.1f} um is below the {MIN_EXTENT_UM:g} "
                "um floor - fewer than ~4 voxels across at the production "
                "dx = 12-30 um. Either this is not a tissue section or the "
                "coordinates are not in micrometres."
            )
        elif ex > MAX_EXTENT_UM or ey > MAX_EXTENT_UM:
            errors.append(
                f"extent {ex:.1f} x {ey:.1f} um exceeds {MAX_EXTENT_UM:g} um "
                "(10 cm), which is larger than any histology section. The "
                "coordinates are almost certainly PIXELS or nanometres, not "
                "micrometres. Convert them before ingesting."
            )
        elif ex < _SUSPECT_SMALL_EXTENT_UM or ey < _SUSPECT_SMALL_EXTENT_UM:
            warnings.append(
                f"extent {ex:.1f} x {ey:.1f} um is small: the sections in this "
                "repo span 5,100-24,000 um. If these coordinates are "
                "millimetres the whole reconstruction is 1000x too small and "
                "nothing downstream will notice."
            )
        density = n / max(ex * ey / 1e6, 1e-12)  # cells per mm^2
        measured["cell_density_per_mm2"] = float(density)
        if density > 1e6:
            warnings.append(
                f"{density:,.0f} cells/mm^2 is far above the ~1-10k/mm^2 of "
                "real tissue; check the coordinate units."
            )

    # --- expression --------------------------------------------------------- #
    if sparse.issparse(expr):
        data = expr.data
        total = float(data.sum())
        nnz = int(expr.nnz)
    else:
        arr = np.asarray(expr)
        data = arr[arr != 0]
        total = float(arr.sum())
        nnz = int(data.size)
    measured["total_counts"] = total
    measured["nonzero_entries"] = nnz
    if n and genes.size:
        measured["matrix_density"] = float(nnz) / float(n * genes.size)
    if total <= 0:
        errors.append(
            "the expression matrix is all zeros: no transcript was measured, "
            "so no phenotype, no ECM field and no capability verdict is "
            "possible."
        )
    elif data.size and float(np.min(data)) < 0:
        errors.append(
            "the expression matrix contains negative values; raw COUNTS are "
            "required (this pipeline log1p-transforms and z-scores them)."
        )
    if n:
        measured["median_counts_per_cell"] = float(
            np.median(np.asarray(expr.sum(axis=1)).ravel()))
        if measured["median_counts_per_cell"] < 5 and total > 0:
            warnings.append(
                f"median {measured['median_counts_per_cell']:.1f} counts per "
                "cell is very low; marker z-scores will be dominated by "
                "sampling noise."
            )

    if genes.size and len(set(genes.tolist())) < genes.size:
        dup = genes.size - len(set(genes.tolist()))
        warnings.append(
            f"{dup} duplicate gene symbols; marker lookup takes the FIRST "
            "column with a given symbol."
        )

    return {"ok": not errors, "errors": errors, "warnings": warnings,
            "measured": measured}


# --------------------------------------------------------------------------- #
# Capability reporting - the part that matters
# --------------------------------------------------------------------------- #
#: Canonical endothelial markers, beyond the two panel maps' choices. Reported
#: separately because ``fields.vascular_density`` - the oxygen SOURCE - is built
#: from the ENDOTHELIAL cells, so a panel with no endothelial marker produces an
#: oxygen field with no supply term (``reconstruct/marker_panels.py`` docstring).
ENDOTHELIAL_MARKER_SUPERSET: Tuple[str, ...] = (
    "PECAM1", "VWF", "FLT1", "PLVAP", "CD34", "CLDN5", "CDH5", "KDR", "ESAM",
    "EGFL7", "ENG", "TEK",
)

#: MGMT is the temozolomide sensitivity determinant. It is on NEITHER 10x Xenium
#: panel (``molecular/mgmt.py``: 0/739 union genes, and no panel gene surrogates
#: for it), so TMZ sensitivity there is a no-resistance UPPER BOUND unless a
#: per-patient scalar is supplied. Reported explicitly because its absence is
#: silent otherwise.
MGMT_GENES: Tuple[str, ...] = ("MGMT",)

#: Every detection rate reported here is over ALL cells, because the ingester
#: runs BEFORE phenotyping - cell types come from ``reconstruct.real_xenium``,
#: which needs the reconstruction that this report exists to guard. The run-time
#: figure in ``reconstruct.cell_state`` is over TUMOR cells only
#: (``scoring_population='tumor'``), so the two differ slightly: on the 319-gene
#: flagship panel this report measures 1.34% and the run-time report 1.42%.
#: Both sit an order of magnitude below the 5% floor, so the VERDICT is the same
#: - but the numbers are not interchangeable and are labelled so.
_ALL_CELLS_NOTE = (
    "rate is over ALL cells: phenotyping has not run yet at ingest time. "
    "reconstruct.cell_state re-measures it over TUMOR cells only "
    "(scoring_population='tumor'), which moves it slightly (319-gene panel: "
    "1.34% here vs 1.42% at run time). Treat this as the PRE-RUN estimate that "
    "decides whether the route is worth attempting."
)


def capability_report(
    genes: Sequence[str],
    expression=None,
    platform: str = "xenium",
    sample_id: str = "sample",
    detection_floor: Optional[float] = None,
) -> Dict:
    """What can this specimen support, and what can it NOT?

    Answers, for the panel in ``genes``:

    * which cell-type marker map fits (``reconstruct.marker_panels``), how well,
      and which :class:`CellType` classes have NO marker and therefore vanish;
    * whether ``state_mode='markers'`` is usable - panel presence AND, when
      ``expression`` is given, the measured detection RATE against
      ``cell_state.PER_CELL_DETECTION_FLOOR``;
    * whether hypoxia is measured at all (``cell_state.HYPOXIA_MARKERS`` and
      ``oxygen_proxy.HYPOXIA_MARKER_GENES``);
    * whether MGMT is on the panel;
    * whether the ECM field will be REAL (matrix transcripts) or a stromal proxy;
    * whether the vascular field - the oxygen SOURCE - has any endothelial
      marker behind it;
    * whether PER-CELL SPATIAL statistics are valid at all, which is a property
      of the PLATFORM and not of the panel: on spot data they are refused
      outright (``results/immune_sink_impact.md`` §5.6).

    Parameters
    ----------
    genes : sequence of str
        Panel gene symbols.
    expression : array or sparse matrix, optional
        (n_cells, n_genes) counts aligned to ``genes``. Supplying it upgrades
        every verdict from PANEL PRESENCE to a measured DETECTION RATE. Without
        it, rate-dependent verdicts are reported as ``UNMEASURED`` with the
        reason - never assumed.
    platform : str
        ``'visium'`` disables the per-cell routes by construction (a spot is not
        a cell), independently of panel content.
    detection_floor : float, optional
        Override ``cell_state.PER_CELL_DETECTION_FLOOR``.

    Returns
    -------
    dict
        ``{sample_id, platform, n_genes, panel, cell_types, proliferation,
        hypoxia, mgmt, ecm, endothelial, supports, headline}``. ``supports`` maps
        a capability name to ``{available: bool, reason: str, evidence: str}``;
        ``headline`` is the ordered list of one-line verdicts that
        :func:`format_capability_report` prints.
    """
    from ..physics.oxygen_proxy import HYPOXIA_MARKER_GENES
    from ..reconstruct.cell_state import (
        CELL_CYCLE_SPECIFIC_MARKERS,
        HYPOXIA_MARKERS,
        HYPOXIA_SPECIFIC_MARKERS,
        LOW_SPECIFICITY_MARKERS,
        PER_CELL_DETECTION_FLOOR,
        PROLIFERATION_MARKERS,
    )
    from ..reconstruct.marker_panels import PANEL_REGISTRY, marker_coverage, select_marker_map
    from ..reconstruct.real_pipeline import ECM_MARKERS

    floor = float(PER_CELL_DETECTION_FLOOR if detection_floor is None
                  else detection_floor)
    gene_list = [str(g) for g in np.asarray(genes, dtype=str).tolist()]
    present = set(gene_list)
    is_spot = str(platform).strip().lower() == "visium"

    rates = _detection_rates(gene_list, expression)
    rate_known = rates is not None

    def split(markers: Sequence[str]) -> Tuple[List[str], List[str]]:
        return ([m for m in markers if m in present],
                [m for m in markers if m not in present])

    def union_rate(markers: Sequence[str]) -> Optional[float]:
        if not rate_known:
            return None
        cols = [gene_list.index(m) for m in markers if m in present]
        if not cols:
            return 0.0
        return _union_detection(expression, cols)

    # --- panel / cell-type map --------------------------------------------- #
    marker_map, panel_name, coverage = select_marker_map(gene_list)
    cell_types: Dict[str, Dict] = {}
    empty_classes: List[str] = []
    for ct, markers in marker_map.items():
        got, missing = split(markers)
        cell_types[ct.name] = {
            "markers": list(markers), "present": got, "missing": missing,
            "detection_rate": union_rate(markers),
        }
        if not got:
            empty_classes.append(ct.name)

    best_cov = coverage.get(panel_name, 0.0)
    panel_block = {
        "selected_map": panel_name,
        "coverage": {k: float(v) for k, v in coverage.items()},
        "selected_coverage": float(best_cov),
        "fits": bool(best_cov >= 0.5),
        "n_genes": len(gene_list),
        "registry": [name for name, _ in PANEL_REGISTRY],
    }
    if best_cov >= 0.99:
        panel_block["verdict"] = (
            f"marker map {panel_name!r} fits exactly (coverage 1.00)")
    elif best_cov >= 0.5:
        panel_block["verdict"] = (
            f"marker map {panel_name!r} is the best fit at coverage "
            f"{best_cov:.2f} - {1 - best_cov:.0%} of its markers are NOT on "
            "this panel and are silently dropped"
        )
    else:
        panel_block["verdict"] = (
            f"NO REGISTERED MARKER MAP FITS: the best is {panel_name!r} at "
            f"coverage {best_cov:.2f} (of "
            f"{ {k: round(v, 2) for k, v in coverage.items()} }). "
            "Phenotyping would run on a map written for a different panel; "
            "add a map for this panel to reconstruct/marker_panels.py before "
            "trusting any composition from this section."
        )

    # --- proliferation ------------------------------------------------------ #
    prol_present, prol_missing = split(PROLIFERATION_MARKERS)
    cyc_present, cyc_missing = split(CELL_CYCLE_SPECIFIC_MARKERS)
    call_set = cyc_present or prol_present
    prol_rate = union_rate(call_set)
    proliferation = {
        "markers": list(PROLIFERATION_MARKERS),
        "present": prol_present, "missing": prol_missing,
        "cycle_specific_present": cyc_present,
        "cycle_specific_missing": cyc_missing,
        "low_specificity_present": split(LOW_SPECIFICITY_MARKERS)[0],
        "call_set": call_set,
        "detection_rate": prol_rate,
        "floor": floor,
        "per_gene_rate": ({m: rates[m] for m in prol_present} if rate_known
                          else None),
        "scoring_population": "all cells",
        "scoring_population_note": _ALL_CELLS_NOTE,
    }

    # --- hypoxia ------------------------------------------------------------ #
    hyp_present, hyp_missing = split(HYPOXIA_MARKERS)
    hypoxia = {
        "markers": list(HYPOXIA_MARKERS),
        "present": hyp_present, "missing": hyp_missing,
        "specific_present": split(HYPOXIA_SPECIFIC_MARKERS)[0],
        "oxygen_proxy_present": split(HYPOXIA_MARKER_GENES)[0],
        "oxygen_proxy_missing": split(HYPOXIA_MARKER_GENES)[1],
        "detection_rate": union_rate(HYPOXIA_MARKERS),
        "per_gene_rate": ({m: rates[m] for m in hyp_present} if rate_known
                          else None),
        "scoring_population": "all cells",
        "scoring_population_note": _ALL_CELLS_NOTE,
    }

    # --- MGMT / ECM / endothelial ------------------------------------------ #
    mgmt_present, mgmt_missing = split(MGMT_GENES)
    mgmt = {
        "present": mgmt_present, "missing": mgmt_missing,
        "detection_rate": union_rate(MGMT_GENES),
    }
    ecm_present, ecm_missing = split(ECM_MARKERS)
    ecm = {"markers": list(ECM_MARKERS), "present": ecm_present,
           "missing": ecm_missing, "detection_rate": union_rate(ECM_MARKERS)}
    endo_present, endo_missing = split(ENDOTHELIAL_MARKER_SUPERSET)
    endothelial = {
        "markers": list(ENDOTHELIAL_MARKER_SUPERSET),
        "present": endo_present, "missing": endo_missing,
        "map_markers_present": cell_types.get("ENDOTHELIAL", {}).get("present", []),
        "detection_rate": union_rate(
            cell_types.get("ENDOTHELIAL", {}).get("present", []) or endo_present),
    }

    # --- capability verdicts ------------------------------------------------ #
    supports: Dict[str, Dict] = {}

    def cap(name: str, available: Optional[bool], reason: str, evidence: str):
        supports[name] = {"available": available, "reason": reason,
                          "evidence": evidence}

    # hypoxia measurement
    if not hyp_present:
        cap("hypoxia_measured", False,
            "hypoxia UNMEASURED: this panel carries NO gene of "
            f"{list(HYPOXIA_MARKERS)}. HYPOXIC would be assigned to NO cell, "
            "which is NOT the same statement as 'computed and found to be "
            "zero'. Do not quote a hypoxic fraction, or an absolute kill, from "
            "this sample (results/cell_state_impact.md §5).",
            "panel presence")
    else:
        cap("hypoxia_measured", True,
            f"hypoxia readable from {hyp_present} - hypoxia-RESPONSE "
            "transcripts, which have real non-hypoxic expression. The HYPOXIC "
            "band is an ENRICHMENT label whose SIZE is UNCALIBRATED."
            + ("" if hypoxia["specific_present"] else
               " None of these is CA9, the only near-hypoxia-restricted gene of "
               "the set, so the call rests entirely on non-specific markers."),
            "panel presence" + (
                f" + detection rate {hypoxia['detection_rate']:.2%}"
                if hypoxia["detection_rate"] is not None else ""))

    # state_mode='markers'
    if is_spot:
        cap("state_mode:markers", False,
            "state_mode='markers' needs per-cell expression. This is SPOT-level "
            "data: a row is a ~55 um spot pooling several cells, so the route "
            "raises rather than degrading. Use 'oxygen' or 'vessel_distance'.",
            "platform")
    elif not prol_present and not hyp_present:
        cap("state_mode:markers", False,
            "state_mode='markers' unavailable: the panel carries neither a "
            f"proliferation transcript ({list(PROLIFERATION_MARKERS)}) nor a "
            f"hypoxia-response transcript ({list(HYPOXIA_MARKERS)}). NOTHING "
            "would be measured, so the route has no content here.",
            "panel presence")
    elif not call_set:
        cap("state_mode:markers", False,
            "state_mode='markers' would measure hypoxia only: no proliferation "
            "transcript is on the panel, so PROLIFERATING and QUIESCENT are "
            "UNMEASURED.",
            "panel presence")
    elif prol_rate is None:
        cap("state_mode:markers", None,
            f"state_mode='markers' UNVERIFIED: {call_set} are on the panel, but "
            "whether a PER-CELL call is allowed depends on their union "
            f"detection rate reaching the {floor:.0%} floor, which is a "
            "property of the COUNT MATRIX and was not read. Re-scan with "
            "detection-rate measurement, or read it from "
            "provenance.cell_state.detail.coverage after a run.",
            "panel presence only (matrix not read)")
    elif prol_rate < floor:
        cap("state_mode:markers", False,
            f"state_mode='markers' unavailable: proliferation detection "
            f"{prol_rate:.2%} ({'+'.join(call_set)}) is below the "
            f"{floor:.1%} per-cell floor. The module would fall back to a "
            "NEIGHBOURHOOD enrichment rule whose ORDERING is measured but whose "
            "LEVEL moves 8.79x with the window radius "
            "(results/cell_state_impact.md §1.1, §2.4), so no absolute "
            "proliferating fraction from this route is quotable.",
            f"measured detection rate {prol_rate:.4%}")
    else:
        cap("state_mode:markers", True,
            f"state_mode='markers' supported: proliferation detection "
            f"{prol_rate:.2%} ({'+'.join(call_set)}) clears the {floor:.1%} "
            "per-cell floor."
            + ("" if cyc_present else
               " NOTE: no MKI67 or TOP2A, so the call rests on "
               "low-specificity markers only."),
            f"measured detection rate {prol_rate:.4%}")

    # oxygen_mode='markers'
    if is_spot:
        cap("oxygen_mode:markers", False,
            "oxygen_mode='markers' is Xenium only: spot-level data carries no "
            "per-cell transcripts.", "platform")
    elif not hypoxia["oxygen_proxy_present"]:
        cap("oxygen_mode:markers", False,
            "oxygen_mode='markers' unavailable: the panel carries none of "
            f"{list(HYPOXIA_MARKER_GENES)}. On this sample hypoxia is "
            "UNMEASURED - use 'pde' or 'vessel_distance'.", "panel presence")
    else:
        cap("oxygen_mode:markers", True,
            f"oxygen_mode='markers' available from "
            f"{hypoxia['oxygen_proxy_present']} - a hypoxia-RESPONSE readout, "
            "rank-normalised, NOT an oxygen tension.", "panel presence")

    # cell typing
    if empty_classes:
        cap("cell_typing", False if "ENDOTHELIAL" in empty_classes else None,
            f"cell types {empty_classes} have NO marker on this panel and will "
            "be assigned to ZERO cells."
            + (" ENDOTHELIAL is one of them, and fields.vascular_density - the "
               "oxygen SOURCE - is built from the endothelial cells, so the "
               "oxygen field would have no supply term at all."
               if "ENDOTHELIAL" in empty_classes else ""),
            "panel presence")
    else:
        cap("cell_typing", panel_block["fits"] or None, panel_block["verdict"],
            f"marker-map coverage {best_cov:.2f}")

    # vasculature
    if not endo_present:
        cap("vascular_field", False,
            "NO endothelial marker on this panel "
            f"(none of {list(ENDOTHELIAL_MARKER_SUPERSET)}). "
            "vascular_density would be identically zero, the oxygen PDE would "
            "have no source, and vessel_distance routes would have no vessels. "
            "This sample cannot support ANY oxygen or delivery result.",
            "panel presence")
    else:
        cap("vascular_field", True,
            f"vasculature measurable from {endo_present[:6]}"
            + (f" (map uses {endothelial['map_markers_present']})"
               if endothelial["map_markers_present"] else ""),
            "panel presence" + (
                f" + endothelial detection {endothelial['detection_rate']:.2%}"
                if endothelial["detection_rate"] is not None else ""))

    # ECM
    if ecm_present:
        cap("ecm_field", True,
            f"ECM field will be the REAL expression-weighted KDE of "
            f"{ecm_present}.", "panel presence")
    else:
        cap("ecm_field", None,
            "no ECM transcript on this panel "
            f"(none of {list(ECM_MARKERS)}), so fields.ecm_density falls back "
            "to the STROMAL-cell density PROXY. The barrier term is then a "
            "proxy for a proxy.", "panel presence")

    # MGMT
    if mgmt_present:
        cap("mgmt_from_panel", True,
            "MGMT IS on this panel, so temozolomide sensitivity can be read "
            "spatially rather than supplied as a per-patient scalar. Note the "
            "calibration in molecular/mgmt.py is fitted to BULK log2(CPM+1), "
            "not to targeted per-cell counts, so a per-cell application is an "
            "extrapolation.",
            "panel presence" + (f" + detection {mgmt['detection_rate']:.2%}"
                                if mgmt["detection_rate"] is not None else ""))
    else:
        cap("mgmt_from_panel", False,
            "MGMT is NOT on this panel, so TMZ resistance is unmeasurable here. "
            "Supply a per-patient MGMT value (RunSpec.mgmt, log2 CPM+1) or "
            "every TMZ kill is a no-resistance UPPER BOUND.",
            "panel presence")

    # sensitivity determinants
    if is_spot:
        cap("sensitivity_mode:determinant", False,
            "spot-level data carries no per-cell expression, so sensitivity "
            "stays uniform 1.0 and the kill is a no-resistance upper bound.",
            "platform")
    else:
        det_present = [g for g in ("TOP2A", "MKI67", "MGMT") if g in present]
        cap("sensitivity_mode:determinant", bool(det_present),
            (f"per-drug determinants measurable: {det_present}. NOT comparable "
             "ACROSS drugs - a drug whose determinant is absent keeps 1.0."
             if det_present else
             "no sensitivity determinant (TOP2A / MKI67 / MGMT) on this panel; "
             "'determinant' mode would be identical to 'uniform'."),
            "panel presence")

    # per-cell spatial statistics: the conclusion of
    # results/immune_sink_impact.md §5.6. A spot is not a cell, and the Visium
    # reconstruction does not measure where the cells inside one are - it
    # SYNTHESISES reconstruct.real_visium._CELLS_PER_SPOT (=12) pseudo-cells per
    # spot at the spot centroid with uniform jitter inside a spot-sized disc.
    # Any statistic over the per-cell spatial distribution therefore measures
    # that jitter, not the tissue, so it must be refused rather than degraded.
    if is_spot:
        cap("spatial_spread", False,
            "PER-CELL SPATIAL STATISTICS ARE INVALID ON THIS SAMPLE. A Visium "
            "spot is a ~55 um well pooling several cells, and this pipeline "
            "does not measure where they sit: reconstruct.real_visium places "
            "12 SYNTHESISED pseudo-cells per spot at the spot centroid with "
            "uniform jitter inside a spot-sized disc. The per-cell attenuation "
            "spread (p90/p10) and the immune-neighbourhood structure ratio "
            "would therefore measure that jitter, not the tissue - "
            "results/immune_sink_impact.md §5.6 declares the statistic INVALID "
            "on spot data and the generator refuses to compute it there. THE "
            "DELIVERY MAP IS INVALID HERE TOO, AND FOR THE SAME REASON: the "
            "exposure field is solved on the grid, but its SOURCE is a density "
            "over ENDOTHELIAL CELL POSITIONS, and on a spot section those are "
            "synthesised like every other cell's. The falloff bins by distance "
            "to the nearest endothelial cell, the perivascular reference is "
            "that curve's first bin, and the delivery spread is a percentile "
            "ratio over the same field - so all three measure the synthesis, "
            "not this tissue's vascular architecture. (An earlier version of "
            "this sentence, written 2026-08-21, claimed the map survived "
            "because it is solved on the grid; that confused where the PDE is "
            "discretised with where its source comes from. The same sentence "
            "previously named the cross-drug RANKING as what survived, which "
            "was withdrawn the same day - results/pk_ranking_robustness.md, "
            "results/capacity_ranking_robustness.md.)",
            "platform")
    else:
        cap("spatial_spread", True,
            "cell positions are MEASURED per cell, so per-cell spatial "
            "statistics (attenuation spread, immune-neighbourhood structure) "
            "are computable on this sample. They must still be read against a "
            "spatially-aware null, not against a well-mixed 1.0 "
            "(results/immune_sink_impact.md §5.1).",
            "platform")

    # absolute kill: the conclusion of results/cell_state_impact.md §5
    if not hyp_present:
        cap("absolute_kill", False,
            "ABSOLUTE KILL IS NOT QUOTABLE FROM THIS SAMPLE UNDER ANY "
            "state_mode. With no hypoxia measurement, no state route clears the "
            "bar: 'oxygen' moves the VIABLE band 0.00-10.04% with the "
            "tissue-mask convention and the mesh, 'markers' has no hypoxia "
            "content, and 'vessel_distance' is a geometric proxy calling "
            "91-100% of tumour cells VIABLE. Kill moves up to 0.3252 absolute "
            "(1.77x) across the three routes. The CROSS-DRUG RANKING is "
            "invariant and is the defensible output.",
            "panel presence")
    else:
        cap("absolute_kill", None,
            "hypoxia is measurable here, so an absolute kill is not "
            "disqualified on panel grounds - but the HYPOXIC band is an "
            "uncalibrated enrichment label, so its magnitude is still not a "
            "measured hypoxic fraction. Prefer the cross-drug ranking.",
            "panel presence")

    headline = _headline(supports)
    return {
        "sample_id": sample_id,
        "platform": platform,
        "n_genes": len(gene_list),
        "detection_rates_measured": bool(rate_known),
        "detection_rate_note": (
            "detection rates measured from the loaded count matrix"
            if rate_known else
            "DETECTION RATES NOT MEASURED: only panel PRESENCE was read "
            "(features.tsv). A gene being on the panel does not mean it fires "
            "often enough for a per-cell call."),
        "panel": panel_block,
        "cell_types": cell_types,
        "proliferation": proliferation,
        "hypoxia": hypoxia,
        "mgmt": mgmt,
        "ecm": ecm,
        "endothelial": endothelial,
        "supports": supports,
        "headline": headline,
    }


#: Order the capability lines are reported in: the things that invalidate a
#: result first, the things that merely degrade it after.
_HEADLINE_ORDER = (
    "vascular_field", "cell_typing", "hypoxia_measured", "absolute_kill",
    "spatial_spread", "state_mode:markers", "oxygen_mode:markers",
    "sensitivity_mode:determinant", "mgmt_from_panel", "ecm_field",
)


def _headline(supports: Dict[str, Dict]) -> List[Dict]:
    out = []
    for name in _HEADLINE_ORDER:
        s = supports.get(name)
        if s is None:
            continue
        status = ("YES" if s["available"] is True
                  else "NO" if s["available"] is False else "CAVEAT")
        out.append({"capability": name, "status": status,
                    "reason": s["reason"], "evidence": s["evidence"]})
    return out


def _detection_rates(gene_list: Sequence[str], expression) -> Optional[Dict[str, float]]:
    """Per-gene fraction of cells with a nonzero count, or ``None``."""
    if expression is None:
        return None
    from scipy import sparse

    n = expression.shape[0]
    if not n:
        return None
    if sparse.issparse(expression):
        counts = np.asarray((expression > 0).sum(axis=0)).ravel()
    else:
        counts = (np.asarray(expression) > 0).sum(axis=0)
    rates = {}
    for j, g in enumerate(gene_list):
        if j < counts.size:
            rates.setdefault(g, float(counts[j]) / float(n))
    return rates


def _union_detection(expression, cols: Sequence[int]) -> float:
    """Fraction of cells with a nonzero count in ANY of ``cols``."""
    from scipy import sparse

    n = expression.shape[0]
    if not n or not len(cols):
        return 0.0
    sub = expression[:, list(cols)]
    if sparse.issparse(sub):
        hit = np.asarray((sub > 0).sum(axis=1)).ravel() > 0
    else:
        hit = (np.asarray(sub) > 0).any(axis=1)
    return float(hit.sum()) / float(n)


def format_capability_report(report: Dict, width: int = 78) -> str:
    """Render a :func:`capability_report` as the plain-text CAPABILITY REPORT.

    Deliberately blunt: the point of this block is that a user reads the NO
    lines before running anything, not after quoting a number.
    """
    import textwrap

    lines: List[str] = []
    push = lines.append
    rule = "=" * width

    push(rule)
    push(f"CAPABILITY REPORT  --  {report['sample_id']}")
    push(rule)
    push(f"platform          : {report['platform']}")
    push(f"panel             : {report['n_genes']} genes; "
         f"map {report['panel']['selected_map']} "
         f"(coverage {report['panel']['selected_coverage']:.2f})")
    push(f"detection rates   : "
         f"{'MEASURED from the count matrix' if report['detection_rates_measured'] else 'NOT MEASURED (panel presence only)'}")
    push("")

    push("CELL TYPES")
    for name, block in report["cell_types"].items():
        rate = block["detection_rate"]
        rate_s = f"  detected in {rate:.2%} of cells" if rate is not None else ""
        got = ", ".join(block["present"]) if block["present"] else "NONE ON PANEL"
        push(f"  {name:<12} {got}{rate_s}")
        if block["missing"]:
            push(f"  {'':<12} (missing: {', '.join(block['missing'])})")
    push("")

    push("MARKER AXES")
    for label, block, extra in (
        ("proliferation", report["proliferation"], "call_set"),
        ("hypoxia", report["hypoxia"], None),
        ("MGMT", report["mgmt"], None),
        ("ECM", report["ecm"], None),
        ("endothelial", report["endothelial"], None),
    ):
        got = block["present"]
        rate = block.get("detection_rate")
        rate_s = (f"   union detection {rate:.2%}" if rate is not None else
                  "   (rate UNMEASURED)")
        push(f"  {label:<14} {', '.join(got) if got else 'NONE ON PANEL'}{rate_s}")
    push("")

    push("VERDICTS")
    for row in report["headline"]:
        mark = {"YES": "[ OK ]", "NO": "[ NO ]", "CAVEAT": "[CAVE]"}[row["status"]]
        head = f"  {mark} {row['capability']}"
        push(head)
        for chunk in textwrap.wrap(row["reason"], width - 9):
            push(f"         {chunk}")
        push(f"         evidence: {row['evidence']}")
    if report["detection_rates_measured"]:
        push("")
        for chunk in textwrap.wrap("NOTE ON RATES: " + _ALL_CELLS_NOTE, width):
            push(chunk)
    else:
        push("")
        for chunk in textwrap.wrap(
                "NOTE ON RATES: " + report["detection_rate_note"], width):
            push(chunk)
    push(rule)
    return "\n".join(lines)


# --------------------------------------------------------------------------- #
# scan: detect + capability, without paying for a full load
# --------------------------------------------------------------------------- #
def scan_report(path: str, measure_detection_rates: bool = False) -> Dict:
    """``detect()`` plus a capability report, for a directory nobody has loaded.

    The cheap path (``measure_detection_rates=False``, the default) reads only
    ``features.tsv.gz``, so it answers in milliseconds even for the 851k-cell
    section - at the cost of reporting every RATE-dependent verdict as
    UNMEASURED rather than guessing it. The expensive path performs a full
    :func:`ingest` and reports measured rates.

    Never raises for a bad directory: the failure is the answer, and it is
    returned in ``issues`` / ``error``.
    """
    det = detect(path)
    out: Dict = {
        "detection": {k: v for k, v in det.items() if k != "files"},
        "files": {k: v for k, v in (det.get("files") or {}).items()
                  if isinstance(v, (str, list))},
        "capability": None,
        "capability_text": None,
        "sample": None,
        "error": None,
    }
    if det["platform"] in ("unknown", "multi_sample") or not det["ingestable"]:
        out["error"] = "; ".join(det["issues"]) or (
            "directory is a container of samples, not a sample"
            if det["platform"] == "multi_sample" else "not ingestable")
        return out

    sid = os.path.basename(os.path.abspath(path).rstrip(os.sep)) or "sample"
    if measure_detection_rates:
        try:
            rec = ingest(path, sample_id=sid, detection=det)
        except IngestError as exc:
            # A DELIBERATE refusal. Informative on purpose, rendered verbatim
            # by both the CLI and the Studio so the two can never disagree
            # about a sample's limits. Never redacted.
            out["error"] = str(exc)
            return out
        except Exception as exc:  # noqa: BLE001 - the failure IS the answer
            # UNEXPECTED. The full detail is kept here because `vivum ingest`
            # prints it to an operator sitting at the console; the HTTP layer
            # (app/server/ingest_registry.scan) replaces it for the browser,
            # keyed off this flag. See cache.redact_exception.
            out["error"] = f"{type(exc).__name__}: {exc}"
            out["error_is_internal"] = True
            return out
        out["capability"] = rec.capability
        out["sample"] = rec.summary()
    else:
        genes = panel_genes(path, det)
        if not genes:
            out["error"] = (
                "panel gene list unreadable without extracting the matrix "
                "tarball; re-scan with measure_detection_rates=true (which "
                "extracts it) to get a capability report."
            )
            return out
        out["capability"] = capability_report(
            genes, expression=None, platform=det["platform"], sample_id=sid)
    out["capability_text"] = format_capability_report(out["capability"])
    return out


# --------------------------------------------------------------------------- #
def _extract_tar(tar_path: str, dest: str) -> None:  # pragma: no cover - helper
    with tarfile.open(tar_path, "r:*") as tf:
        tf.extractall(dest)


if __name__ == "__main__":  # pragma: no cover - smoke test on the real caches
    import sys

    targets = sys.argv[1:] or [
        os.path.join("data_cache", "xenium_gbm"),
        os.path.join("data_cache", "xenium_gbm_io"),
        os.path.join("data_cache", "gse324268_xenium"),
        os.path.join("data_cache", "visium_gbm"),
        os.path.join("data_cache", "visium_gbm", "GBM1"),
        os.path.join("data_cache", "visium_gbm", "GBM2"),
        os.path.join("data_cache", "visium_gbm", "GBM3"),
    ]
    for t in targets:
        d = detect(t)
        print(f"\n{t}\n  platform={d['platform']} layout={d['layout']} "
              f"confidence={d['confidence']} ingestable={d['ingestable']}")
        for i in d["issues"]:
            print(f"  ISSUE: {i}")
        rec = ingest(t)
        print(f"  n_cells={rec.n_cells:,} n_genes={rec.n_genes} "
              f"extent={rec.extent_um[0]:.0f} x {rec.extent_um[1]:.0f} um")
        print(format_capability_report(rec.capability))
