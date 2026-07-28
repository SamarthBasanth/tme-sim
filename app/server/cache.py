"""Section registry + the reconstruction cache.

Two responsibilities, both about the EXPENSIVE half of a run:

1. **Which specimens are actually loadable.** :func:`list_sections` enumerates
   the sections this checkout can reconstruct by looking at the bytes in
   ``data_cache/``. A section whose cache is missing is reported with
   ``available: false`` and a ``reason``; it never raises and never disappears
   from the list, so the UI can show *why* something is unusable.

2. **Caching reconstructions by ``(section, dx_um, vascular_mode,
   oxygen_mode)``.** A reconstruction takes 1 s (40k-cell brain section at
   dx = 30 um) to ~90 s (851k-cell GSE324268 section at dx = 12 um) and is
   IDENTICAL across every drug run on the same geometry, so caching it is the
   single biggest responsiveness win available. The cache is in-process (the
   objects are large numpy arrays; a disk round-trip would cost most of the
   saving) and LRU-bounded by entry count.

Everything here WRAPS the production pipeline. ``reconstruct_real`` is the only
thing that builds a microenvironment; the one place this module reaches past it
is ``vascular_mode='kde'``, which ``reconstruct_real`` does not expose as a
parameter — see :func:`build_reconstruction` for exactly what is re-run and why.
"""
from __future__ import annotations

import json
import os
import threading
import time
from collections import OrderedDict
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np

__all__ = [
    "SECTIONS",
    "VASCULAR_MODES",
    "OXYGEN_MODES",
    "SENSITIVITY_MODES",
    "ReconKey",
    "ReconstructionCache",
    "RECON_CACHE",
    "build_reconstruction",
    "get_section",
    "list_sections",
    "section_genes",
    "repo_root",
    "runs_dir",
]

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #
_HERE = os.path.dirname(os.path.abspath(__file__))


def repo_root() -> str:
    """Absolute path to the repository root (the parent of ``app/``)."""
    return os.path.abspath(os.path.join(_HERE, "..", ".."))


def runs_dir() -> str:
    """Absolute path to ``app/runs/`` (created on demand)."""
    path = os.path.join(repo_root(), "app", "runs")
    os.makedirs(path, exist_ok=True)
    return path


def _abs(*parts: str) -> str:
    return os.path.join(repo_root(), *parts)


# --------------------------------------------------------------------------- #
# Enumerations mirrored from the engine (imported lazily where they are used so
# importing this module stays cheap; these literals are asserted against the
# engine's own tuples in tests/test_app_server.py).
# --------------------------------------------------------------------------- #
VASCULAR_MODES: Tuple[str, ...] = ("discrete", "kde")
OXYGEN_MODES: Tuple[str, ...] = ("pde", "vessel_distance", "markers")
SENSITIVITY_MODES: Tuple[str, ...] = ("uniform", "determinant")


# --------------------------------------------------------------------------- #
# Section registry
# --------------------------------------------------------------------------- #
#: Every specimen this backend knows how to reconstruct. ``source`` and
#: ``sample`` are exactly the arguments handed to
#: ``reconstruct.real_pipeline.reconstruct_real``; nothing is re-derived.
SECTIONS: "OrderedDict[str, Dict]" = OrderedDict([
    ("xenium_319", {
        "id": "xenium_319",
        "label": "Xenium GBM · 319-panel",
        "source": "xenium",
        "sample": None,  # loader default -> data_cache/xenium_gbm
        "cache_dir": os.path.join("data_cache", "xenium_gbm"),
        "panel": "hBrain_254g",
        "assay": "Xenium In Situ (single-cell)",
        "accession": "10x Genomics Xenium Human Brain Preview (GBM section)",
        "has_outcome_label": False,
    }),
    ("xenium_io", {
        "id": "xenium_io",
        "label": "Xenium GBM · IO 480-panel",
        "source": "xenium",
        "sample": os.path.join("data_cache", "xenium_gbm_io"),
        "cache_dir": os.path.join("data_cache", "xenium_gbm_io"),
        "panel": "hImmune_v1 + Add-on B",
        "assay": "Xenium In Situ (single-cell)",
        "accession": "10x Genomics FFPE Human Brain Cancer, Immuno-Oncology panel",
        "has_outcome_label": False,
    }),
    ("gse324268", {
        "id": "gse324268",
        "label": "Xenium GBM · GSE324268 tumor 07",
        "source": "xenium",
        "sample": os.path.join("data_cache", "gse324268_xenium"),
        "cache_dir": os.path.join("data_cache", "gse324268_xenium"),
        "panel": None,  # read from the shipped gene_panel.json at enumeration
        "assay": "Xenium In Situ (single-cell)",
        "accession": "GEO GSE324268 (GSM9572843, GBMtumor07)",
        "has_outcome_label": False,
    }),
    ("visium_10x_parent_gbm", {
        "id": "visium_10x_parent_gbm",
        "label": "Visium GBM · 10x parent section",
        "source": "visium",
        "sample": "10x_parent_gbm",
        "cache_dir": os.path.join("data_cache", "visium_gbm"),
        "panel": "whole transcriptome (Visium v1, 55 um spots)",
        "assay": "Visium spatial gene expression (spot-level)",
        "accession": "Parent_Visium_Human_Glioblastoma",
        "has_outcome_label": False,
    }),
    ("visium_gbm1", {
        "id": "visium_gbm1",
        "label": "Visium GBM · GBM1 (primary)",
        "source": "visium",
        "sample": "GBM1",
        "cache_dir": os.path.join("data_cache", "visium_gbm", "GBM1"),
        "panel": "whole transcriptome (Visium v1, 55 um spots)",
        "assay": "Visium spatial gene expression (spot-level)",
        "accession": "GEO GSE194329 / GSM5833533",
        "has_outcome_label": False,
    }),
    ("visium_gbm2", {
        "id": "visium_gbm2",
        "label": "Visium GBM · GBM2 (recurrent)",
        "source": "visium",
        "sample": "GBM2",
        "cache_dir": os.path.join("data_cache", "visium_gbm", "GBM2"),
        "panel": "whole transcriptome (Visium v1, 55 um spots)",
        "assay": "Visium spatial gene expression (spot-level)",
        "accession": "GEO GSE194329 / GSM5833534",
        "has_outcome_label": False,
    }),
    ("visium_gbm3", {
        "id": "visium_gbm3",
        "label": "Visium GBM · GBM3 (primary)",
        "source": "visium",
        "sample": "GBM3",
        "cache_dir": os.path.join("data_cache", "visium_gbm", "GBM3"),
        "panel": "whole transcriptome (Visium v1, 55 um spots)",
        "assay": "Visium spatial gene expression (spot-level)",
        "accession": "GEO GSE194329 / GSM5833535",
        "has_outcome_label": False,
    }),
])

#: Files a Xenium cache directory must contain to be loadable by
#: ``tmesim.data.real.singlecell_spatial.load``.
_XENIUM_REQUIRED = ("Xenium_GBM_cells.csv.gz",)
_XENIUM_MATRIX_FILES = ("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz")

#: On-disk memo of the numbers that need a file read (cell counts, extents,
#: panel gene lists). Keyed by section id + a fingerprint of the source files,
#: so it self-invalidates when the cache directory changes.
_META_FILE = "_section_meta.json"

# Re-entrant so _save_meta can prime the memo through _load_meta first.
_meta_lock = threading.RLock()
_meta_cache: Optional[Dict] = None


def _load_meta() -> Dict:
    global _meta_cache
    with _meta_lock:
        if _meta_cache is None:
            path = os.path.join(runs_dir(), _META_FILE)
            try:
                with open(path) as fh:
                    _meta_cache = json.load(fh)
            except Exception:  # noqa: BLE001 - a corrupt memo is not fatal
                _meta_cache = {}
        return _meta_cache


def _save_meta(section_id: str, entry: Dict) -> None:
    _load_meta()  # prime the memo so a save never drops other sections
    with _meta_lock:
        meta = _meta_cache if _meta_cache is not None else {}
        meta[section_id] = entry
        globals()["_meta_cache"] = meta
        path = os.path.join(runs_dir(), _META_FILE)
        try:
            tmp = path + ".tmp"
            with open(tmp, "w") as fh:
                json.dump(meta, fh, indent=1, sort_keys=True)
            os.replace(tmp, path)
        except Exception:  # noqa: BLE001 - the memo is an optimisation only
            pass


def _fingerprint(paths: List[str]) -> str:
    """(size, mtime) digest of the files a section's metadata was read from."""
    parts = []
    for p in sorted(paths):
        try:
            st = os.stat(p)
            parts.append(f"{os.path.basename(p)}:{st.st_size}:{int(st.st_mtime)}")
        except OSError:
            parts.append(f"{os.path.basename(p)}:missing")
    return "|".join(parts)


# --------------------------------------------------------------------------- #
# Availability + cheap metadata, per source
# --------------------------------------------------------------------------- #
def _xenium_paths(rec: Dict) -> Dict[str, str]:
    root = _abs(rec["cache_dir"])
    return {
        "root": root,
        "cells": os.path.join(root, "Xenium_GBM_cells.csv.gz"),
        "matrix_dir": os.path.join(root, "cell_feature_matrix"),
        "tar": os.path.join(root, "Xenium_GBM_cell_feature_matrix.tar.gz"),
        "panel": os.path.join(root, "Xenium_GBM_gene_panel.json"),
    }


def _xenium_available(rec: Dict) -> Tuple[bool, Optional[str]]:
    p = _xenium_paths(rec)
    if not os.path.isdir(p["root"]):
        return False, (
            f"cache directory {rec['cache_dir']} does not exist; see "
            "docs/real_data/singlecell_spatial.md for the download list"
        )
    if not os.path.exists(p["cells"]):
        return False, (
            f"missing {os.path.relpath(p['cells'], repo_root())} "
            "(the per-cell centroid table)"
        )
    has_mex = all(
        os.path.exists(os.path.join(p["matrix_dir"], f))
        for f in _XENIUM_MATRIX_FILES
    )
    if not has_mex and not os.path.exists(p["tar"]):
        return False, (
            "missing the cell-feature matrix: neither the extracted MEX triple "
            f"in {rec['cache_dir']}/cell_feature_matrix nor the tarball"
        )
    return True, None


def _xenium_panel_name(panel_json_path: str) -> Optional[str]:
    """Panel identity as shipped in the section's ``gene_panel.json``.

    10x records the add-on panel at ``payload.panel.identity.name`` and the base
    panel it extends at ``payload.panel.type.data.base_panel_identity.name``, so
    the honest label is "base + add-on". Returns ``None`` (and the registry's
    static name is kept) if the file is absent or shaped differently.
    """
    if not os.path.exists(panel_json_path):
        return None
    try:
        with open(panel_json_path) as fh:
            doc = json.load(fh)
        panel = ((doc.get("payload") or {}).get("panel") or {})
        addon = (panel.get("identity") or {}).get("name")
        base = ((((panel.get("type") or {}).get("data") or {})
                 .get("base_panel_identity") or {}).get("name"))
        parts = [p for p in (base, addon) if p]
        return " + ".join(parts) or None
    except Exception:  # noqa: BLE001 - the registry name stays
        return None


def _xenium_metadata(rec: Dict) -> Dict:
    """n_cells, measured extent, panel name and gene list, read from disk.

    The extent is the RAW measured centroid bounding-box span, which is
    independent of ``dx_um`` — unlike ``Domain.extent_um``, which is the span
    rounded up to a whole number of grid cells and therefore moves with the
    mesh.
    """
    import gzip

    import pandas as pd

    p = _xenium_paths(rec)
    cells = pd.read_csv(p["cells"], usecols=["x_centroid", "y_centroid"])
    x = cells["x_centroid"].to_numpy(dtype=np.float64)
    y = cells["y_centroid"].to_numpy(dtype=np.float64)

    genes: List[str] = []
    features = os.path.join(p["matrix_dir"], "features.tsv.gz")
    if os.path.exists(features):
        with gzip.open(features, "rt") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 3 and parts[2].strip() == "Gene Expression":
                    genes.append(parts[1])

    panel = _xenium_panel_name(p["panel"]) or rec.get("panel")

    return {
        "n_cells": int(x.size),
        "extent_um": [float(x.max() - x.min()), float(y.max() - y.min())],
        "n_genes": len(genes) or None,
        "genes": genes,
        "panel": panel,
        "n_spots": None,
        "measured_from": "cells.csv.gz centroids + features.tsv.gz",
    }


def _visium_available(rec: Dict) -> Tuple[bool, Optional[str]]:
    try:
        from tmesim.data.real import visium
    except Exception as exc:  # noqa: BLE001
        return False, f"visium loader unimportable: {type(exc).__name__}: {exc}"
    try:
        present = set(visium.list_samples(present_only=True))
    except Exception as exc:  # noqa: BLE001
        return False, f"visium cache scan failed: {type(exc).__name__}: {exc}"
    if rec["sample"] in present:
        return True, None
    return False, (
        f"Visium sample {rec['sample']!r} is not cached under "
        f"{rec['cache_dir']}; see docs/real_data/visium_gbm.md"
    )


def _visium_metadata(rec: Dict) -> Dict:
    """Spot count + measured extent, without loading the count matrix.

    ``n_cells`` is deliberately ``None``: a Visium reconstruction does not have
    measured cells, it PLACES cells from the per-spot deconvolved fractions, so
    the count does not exist until ``reconstruct_from_visium`` has run. It is
    filled in from the reconstruction the first time one is built.
    """
    out = {
        "n_cells": None,
        "n_cells_reason": (
            "Visium is spot-level: cells are placed by deconvolution at "
            "reconstruction time, so no cell count exists in the cached bytes. "
            "Populated after this section is reconstructed once."
        ),
        "extent_um": None,
        "n_genes": None,
        "genes": [],
        "panel": rec.get("panel"),
        "n_spots": None,
        "measured_from": "tissue_positions + scalefactors_json",
    }
    try:
        from tmesim.data.real import visium

        spatial_dir, _matrix_dir = visium._locate(rec["sample"])
        coords_px = visium._read_positions(spatial_dir)
        umpp, _diam = visium._read_umpp(spatial_dir)
        pts = np.asarray(list(coords_px.values()), dtype=np.float64) * float(umpp)
        if pts.size:
            out["n_spots"] = int(pts.shape[0])
            out["extent_um"] = [
                float(pts[:, 0].max() - pts[:, 0].min()),
                float(pts[:, 1].max() - pts[:, 1].min()),
            ]
    except Exception as exc:  # noqa: BLE001 - degrade to nulls, never crash
        out["extent_um_reason"] = (
            f"spot positions unreadable: {type(exc).__name__}: {exc}"
        )
    return out


def _metadata_source_files(rec: Dict) -> List[str]:
    if rec["source"] == "xenium":
        p = _xenium_paths(rec)
        return [p["cells"], os.path.join(p["matrix_dir"], "features.tsv.gz"),
                p["panel"]]
    return [_abs(rec["cache_dir"])]


def _section_metadata(rec: Dict) -> Dict:
    """Cheap on-disk metadata for one section, memoised in ``app/runs``."""
    fp = _fingerprint(_metadata_source_files(rec))
    memo = _load_meta().get(rec["id"])
    if memo and memo.get("_fingerprint") == fp:
        return memo

    if rec["source"] == "xenium":
        meta = _xenium_metadata(rec)
    else:
        meta = _visium_metadata(rec)
    meta["_fingerprint"] = fp
    if memo:
        # Preserve anything only a completed reconstruction can know.
        for k in ("n_cells_observed", "grid_observed"):
            if memo.get(k) is not None:
                meta[k] = memo[k]
    _save_meta(rec["id"], meta)
    return meta


def record_observed_reconstruction(section_id: str, report: Dict) -> None:
    """Remember what a completed reconstruction measured for this section.

    Lets ``/api/sections`` report a real Visium cell count (which cannot be
    known from the cached bytes) instead of ``null`` once one has been built.
    Only ever records numbers that were actually produced by
    ``reconstruction_report``.
    """
    rec = SECTIONS.get(section_id)
    if rec is None:
        return
    meta = dict(_section_metadata(rec))
    meta["n_cells_observed"] = int(report.get("n_cells") or 0) or None
    grid = report.get("grid")
    meta["grid_observed"] = list(grid) if grid else None
    _save_meta(section_id, meta)


def section_genes(section_id: str) -> List[str]:
    """Gene symbols on this section's panel (empty list if not knowable cheaply).

    Used by ``/api/validate`` to reject ``oxygen_mode='markers'`` on a panel
    that carries no hypoxia-response transcript BEFORE a 40 s run discovers it.
    """
    rec = SECTIONS.get(section_id)
    if rec is None:
        return []
    try:
        return list(_section_metadata(rec).get("genes") or [])
    except Exception:  # noqa: BLE001
        return []


def get_section(section_id: str) -> Optional[Dict]:
    """The registry record for ``section_id``, or ``None`` if unknown."""
    return SECTIONS.get(section_id)


def describe_section(section_id: str) -> Dict:
    """One ``Section`` object exactly as ``GET /api/sections`` returns it."""
    rec = SECTIONS[section_id]
    if rec["source"] == "xenium":
        available, reason = _xenium_available(rec)
    else:
        available, reason = _visium_available(rec)

    out = {
        "id": rec["id"],
        "label": rec["label"],
        "n_cells": None,
        "extent_um": None,
        "panel": rec.get("panel"),
        "available": bool(available),
        # --- additive provenance (the contract's four fields are above) ---
        "source": rec["source"],
        "sample": rec["sample"],
        "cache_dir": rec["cache_dir"],
        "assay": rec.get("assay"),
        "accession": rec.get("accession"),
        "has_outcome_label": rec.get("has_outcome_label", False),
        "reason": reason,
        "n_spots": None,
        "n_genes": None,
    }
    if not available:
        out["n_cells_reason"] = reason
        out["extent_um_reason"] = reason
        return out

    try:
        meta = _section_metadata(rec)
    except Exception as exc:  # noqa: BLE001 - a metadata read must not 500
        out["available"] = False
        out["reason"] = (
            f"cache present but unreadable: {type(exc).__name__}: {exc}"
        )
        out["n_cells_reason"] = out["reason"]
        out["extent_um_reason"] = out["reason"]
        return out

    n_cells = meta.get("n_cells")
    if n_cells is None:
        n_cells = meta.get("n_cells_observed")
    out["n_cells"] = n_cells
    if n_cells is None and meta.get("n_cells_reason"):
        out["n_cells_reason"] = meta["n_cells_reason"]
    out["extent_um"] = meta.get("extent_um")
    if out["extent_um"] is None and meta.get("extent_um_reason"):
        out["extent_um_reason"] = meta["extent_um_reason"]
    out["n_spots"] = meta.get("n_spots")
    out["n_genes"] = meta.get("n_genes")
    if meta.get("panel"):
        out["panel"] = meta["panel"]
    return out


def list_sections() -> List[Dict]:
    """Every registered section, available or not. Never raises."""
    out = []
    for sid in SECTIONS:
        try:
            out.append(describe_section(sid))
        except Exception as exc:  # noqa: BLE001
            rec = SECTIONS[sid]
            out.append({
                "id": sid, "label": rec["label"], "n_cells": None,
                "extent_um": None, "panel": rec.get("panel"),
                "available": False,
                "reason": f"enumeration failed: {type(exc).__name__}: {exc}",
                "source": rec["source"], "sample": rec["sample"],
            })
    return out


# --------------------------------------------------------------------------- #
# Reconstruction
# --------------------------------------------------------------------------- #
ReconKey = Tuple[str, float, str, str]


def make_key(section: str, dx_um: float, vascular_mode: str,
             oxygen_mode: str) -> ReconKey:
    """The cache key. ``dx_um`` is rounded so 12.0 and 12.000000001 share an
    entry rather than silently rebuilding a 90 s reconstruction."""
    return (str(section), round(float(dx_um), 6),
            str(vascular_mode).strip().lower(), str(oxygen_mode).strip().lower())


def _max_darcy_speed(micro) -> Optional[float]:
    """Peak |vx| + |vy| (um/s) of the interstitial Darcy velocity field.

    Computed with the production ``physics.pressure.solve_ifp`` — the same call
    the engine makes — so the advection CFL that ``/api/validate`` reports is
    the one ``physics.diffusion.simulate_drug_transport`` will actually apply.
    Returns ``None`` (never a guess) if the solve fails.
    """
    try:
        from tmesim.physics.pressure import solve_ifp

        _ifp, velocity = solve_ifp(micro)
        v = np.asarray(velocity, dtype=np.float64)
        return float(np.max(np.abs(v[0]) + np.abs(v[1]))) if v.size else 0.0
    except Exception:  # noqa: BLE001
        return None


def build_reconstruction(
    section_id: str,
    dx_um: float,
    vascular_mode: str,
    oxygen_mode: str,
    log: Optional[Callable[[str], None]] = None,
) -> Dict:
    """Reconstruct one specimen through the production pipeline.

    ``reconstruct_real`` does the whole job for ``vascular_mode='discrete'``
    (its own default) at any ``oxygen_mode``, so that path is the unmodified
    production call and its numbers are bit-identical to the CLI.

    ``vascular_mode='kde'`` is the one case ``reconstruct_real`` cannot express:
    it hardcodes ``vascular_density_from_cells(domain, cells)``, i.e. the
    discrete segmented-vessel field, with no parameter to select the historical
    KDE. Rather than reimplement the assembler, this wrapper runs the production
    reconstruction and then swaps in
    ``real_fields.vascular_density_from_cells(..., mode='kde')`` — the same
    public function, the same measured endothelial cells — and re-installs
    oxygen through ``real_pipeline._install_oxygen`` because the ``'pde'`` and
    ``'vessel_distance'`` routes are driven BY the vascular field.
    ``oxygen_mode='markers'`` reads hypoxia off transcripts and does not depend
    on vasculature, so it is left alone. The cost of the swap (a discarded
    discrete field, and one discarded oxygen solve) is logged, not hidden.

    Returns ``{micro, report, build_time_s, max_darcy_speed_um_s, notes}``.
    """
    from tmesim.core import OxygenConfig
    from tmesim.reconstruct import real_pipeline
    from tmesim.reconstruct.real_fields import vascular_density_from_cells

    def _log(msg: str) -> None:
        if log is not None:
            log(msg)

    rec = SECTIONS.get(section_id)
    if rec is None:
        raise KeyError(
            f"unknown section {section_id!r}; known: {sorted(SECTIONS)}"
        )
    vmode = str(vascular_mode).strip().lower()
    omode = str(oxygen_mode).strip().lower()
    if vmode not in VASCULAR_MODES:
        raise ValueError(
            f"unknown vascular_mode {vascular_mode!r}; expected one of "
            f"{list(VASCULAR_MODES)}"
        )

    notes: List[str] = []
    t0 = time.time()

    sample = rec["sample"]
    # Xenium sample paths in the registry are repo-relative; make them absolute
    # so the server works regardless of the process working directory.
    if rec["source"] == "xenium" and sample is not None:
        sample = _abs(sample)

    _log(f"reconstruct_real(source={rec['source']!r}, sample={sample!r}, "
         f"dx_um={dx_um}, oxygen_mode={omode!r})")
    micro = real_pipeline.reconstruct_real(
        source=rec["source"], sample=sample, dx_um=float(dx_um),
        oxygen_mode=omode,
    )

    if vmode == "kde":
        if rec["source"] != "xenium":
            notes.append(
                "vascular_mode is not applicable to a Visium section: "
                "reconstruct_from_visium builds vascular_density from the "
                "measured ENDOTHELIAL spot fraction and exposes no mode "
                "switch. The requested 'kde' was NOT applied."
            )
            _log("vascular_mode='kde' ignored (Visium has no mode switch)")
        else:
            _log("re-deriving vascular_density with mode='kde' "
                 "(reconstruct_real hardcodes 'discrete'); the discrete field "
                 "and its oxygen solve are discarded")
            micro.fields.vascular_density = vascular_density_from_cells(
                micro.domain, micro.cells, mode="kde"
            )
            meta = dict(micro.metadata)
            rp = dict(meta.get("real_pipeline", {}))
            rp["vascular_source"] = "KDE of REAL ENDOTHELIAL cell positions"
            if omode in ("pde", "vessel_distance"):
                _log(f"re-solving oxygen ({omode}) on the KDE vascular field")
                rp.update(real_pipeline._install_oxygen(
                    micro, omode, OxygenConfig()
                ))
            else:
                notes.append(
                    "oxygen_mode='markers' reads hypoxia off measured "
                    "transcripts and does not depend on the vascular field, "
                    "so it was not re-solved after the KDE swap."
                )
            meta["real_pipeline"] = rp
            micro.metadata = meta
            notes.append(
                "vascular_density was re-derived as the historical 40 um KDE "
                "of the measured endothelial cells AFTER the production "
                "reconstruction, because reconstruct_real() does not expose a "
                "vascular_mode parameter. The ECM field and the cell "
                "reconstruction are the unmodified production ones."
            )
    else:
        meta = dict(micro.metadata)
        rp = dict(meta.get("real_pipeline", {}))
        rp.setdefault("vascular_source", "DISCRETE segmented vessels")
        meta["real_pipeline"] = rp
        micro.metadata = meta

    report = real_pipeline.reconstruction_report(micro)
    report["vascular_mode"] = vmode
    report["vascular_source"] = micro.metadata.get(
        "real_pipeline", {}).get("vascular_source")
    build_time = time.time() - t0
    _log(f"reconstruction done in {build_time:.1f}s: "
         f"{report['n_cells']} cells, grid {report['grid']}")

    speed = _max_darcy_speed(micro)
    try:
        record_observed_reconstruction(section_id, report)
    except Exception:  # noqa: BLE001 - bookkeeping must not fail a run
        pass

    return {
        "micro": micro,
        "report": report,
        "build_time_s": float(build_time),
        "max_darcy_speed_um_s": speed,
        "notes": notes,
        # Snapshot of the two things a run mutates on the shared cached object,
        # so every run can start from exactly what reconstruct_real produced.
        "pristine_sensitivity": np.array(micro.cells.sensitivity, copy=True),
        "pristine_metadata": dict(micro.metadata),
    }


# --------------------------------------------------------------------------- #
# The cache
# --------------------------------------------------------------------------- #
class ReconstructionCache:
    """Thread-safe LRU of reconstructions keyed by
    ``(section, dx_um, vascular_mode, oxygen_mode)``.

    Two jobs that want the same key do NOT both reconstruct: the second blocks
    on a per-key lock and then finds the first one's entry (single-flight). The
    default capacity is small because each entry holds the full cell population
    and several float64 grids — the 851k-cell section at dx = 12 um is on the
    order of a hundred MB.

    A tiny, disk-backed sidecar records each key's measured peak Darcy speed so
    ``/api/validate`` can still check the advection CFL after a restart, before
    anything has been rebuilt. Only measured values are written; a key that has
    never been built has no entry, and validation says so.
    """

    #: Sidecar file inside ``app/runs/`` (see the class docstring).
    CFL_FILE = "_cfl_speeds.json"

    def __init__(self, maxsize: int = 4):
        self.maxsize = int(maxsize)
        self._entries: "OrderedDict[ReconKey, Dict]" = OrderedDict()
        self._lock = threading.RLock()
        self._key_locks: Dict[ReconKey, threading.Lock] = {}
        self.hits = 0
        self.misses = 0
        self.evictions = 0
        self._speeds: Optional[Dict[str, float]] = None

    # -- speed sidecar ----------------------------------------------------- #
    @staticmethod
    def _speed_key(key: ReconKey) -> str:
        return "|".join(str(k) for k in key)

    def _speed_path(self) -> str:
        return os.path.join(runs_dir(), self.CFL_FILE)

    def _load_speeds(self) -> Dict[str, float]:
        if self._speeds is None:
            try:
                with open(self._speed_path()) as fh:
                    self._speeds = json.load(fh)
            except Exception:  # noqa: BLE001
                self._speeds = {}
        return self._speeds

    def _store_speed(self, key: ReconKey, speed: Optional[float]) -> None:
        if speed is None:
            return
        speeds = self._load_speeds()
        speeds[self._speed_key(key)] = float(speed)
        try:
            path = self._speed_path()
            tmp = path + ".tmp"
            with open(tmp, "w") as fh:
                json.dump(speeds, fh, indent=1, sort_keys=True)
            os.replace(tmp, path)
        except Exception:  # noqa: BLE001
            pass

    def known_speed(self, key: ReconKey) -> Optional[float]:
        """Measured peak Darcy speed for this exact key, or ``None``."""
        with self._lock:
            entry = self._entries.get(key)
            if entry is not None and entry.get("max_darcy_speed_um_s") is not None:
                return float(entry["max_darcy_speed_um_s"])
            return self._load_speeds().get(self._speed_key(key))

    def speeds_for_section(self, section_id: str,
                           vascular_mode: Optional[str] = None
                           ) -> Dict[str, float]:
        """Every measured peak Darcy speed recorded for ``section_id``, keyed by
        the stringified cache key.

        ``vascular_mode`` narrows it to keys whose vascular field was built the
        same way — the interstitial-pressure solve is driven by
        ``fields.vascular_density``, so 'discrete' and 'kde' produce genuinely
        different velocity fields and their speeds should not be mixed. The
        remaining variation across such keys is ``dx`` (and ``oxygen_mode``,
        which ``physics.pressure.solve_ifp`` does not read at all).

        Used to reason about the CFL for a key that has not been built yet;
        every message built from this is labelled with the key it came from and
        never presented as a measurement of the key being validated.
        """
        with self._lock:
            speeds = dict(self._load_speeds())
        prefix = f"{section_id}|"
        out = {k: v for k, v in speeds.items() if k.startswith(prefix)}
        if vascular_mode is None:
            return out
        vm = str(vascular_mode).strip().lower()
        return {k: v for k, v in out.items() if k.split("|")[2:3] == [vm]}

    # -- cache ------------------------------------------------------------- #
    def peek(self, key: ReconKey) -> Optional[Dict]:
        """Return a cached entry without building it (and without counting a
        miss). Refreshes its LRU position on a hit."""
        with self._lock:
            entry = self._entries.get(key)
            if entry is not None:
                self._entries.move_to_end(key)
            return entry

    def _key_lock(self, key: ReconKey) -> threading.Lock:
        with self._lock:
            lock = self._key_locks.get(key)
            if lock is None:
                lock = threading.Lock()
                self._key_locks[key] = lock
            return lock

    def get_or_build(
        self,
        key: ReconKey,
        log: Optional[Callable[[str], None]] = None,
    ) -> Dict:
        """Return the reconstruction for ``key``, building it if necessary.

        The returned dict is the one held in the cache — callers must treat
        ``micro`` as read-only shared state. ``simulate_response`` already works
        on a private copy, and ``apply_drug_sensitivity`` writes only
        ``cells.sensitivity``, which every run sets before simulating.
        """
        entry = self.peek(key)
        if entry is not None:
            with self._lock:
                self.hits += 1
            if log is not None:
                log(f"reconstruction CACHE HIT for {key} "
                    f"(saved ~{entry['build_time_s']:.1f}s)")
            return entry

        with self._key_lock(key):
            entry = self.peek(key)
            if entry is not None:  # another job built it while we waited
                with self._lock:
                    self.hits += 1
                if log is not None:
                    log(f"reconstruction CACHE HIT (after wait) for {key}")
                return entry

            with self._lock:
                self.misses += 1
            if log is not None:
                log(f"reconstruction CACHE MISS for {key}; building")
            section, dx_um, vmode, omode = key
            built = build_reconstruction(section, dx_um, vmode, omode, log=log)
            built["key"] = key
            built["created_at"] = time.time()
            with self._lock:
                self._entries[key] = built
                self._entries.move_to_end(key)
                while len(self._entries) > self.maxsize:
                    old_key, _old = self._entries.popitem(last=False)
                    self.evictions += 1
                    if log is not None:
                        log(f"evicted reconstruction {old_key} (LRU, "
                            f"maxsize={self.maxsize})")
            self._store_speed(key, built.get("max_darcy_speed_um_s"))
            return built

    def stats(self) -> Dict:
        with self._lock:
            return {
                "hits": self.hits,
                "misses": self.misses,
                "evictions": self.evictions,
                "size": len(self._entries),
                "maxsize": self.maxsize,
                "keys": [list(k) for k in self._entries],
            }

    def clear(self) -> None:
        with self._lock:
            self._entries.clear()
            self._key_locks.clear()


#: Process-wide cache used by the job worker and by ``/api/validate``.
RECON_CACHE = ReconstructionCache(
    maxsize=int(os.environ.get("TMESIM_RECON_CACHE_SIZE", "4"))
)
