"""Adapter: Akoya CODEX (PhenoCycler) multiplex-protein OME-TIFF -> an
``ingest.SampleRecord``-shaped object usable by
``tmesim.reconstruct.registration_3d.register_serial_sections``.

Background (full account in ``results/gbm_3d_registration.md``): Greenwald,
Darnell, Hoefflin et al. (Tirosh/Suva labs), *Cell* 2024 (PMID 38653236)
profiled 13 IDH-wildtype GBM tumours by 10x Visium AND, on near-adjacent
sections of the SAME blocks, by CODEX (40-marker spatial proteomics), then
cross-registered the two with STalign. The Visium layer ingests through
``tmesim/data/ingest.py`` UNCHANGED once its flat Zenodo-tarball layout is
restructured into a Space Ranger ``outs/`` shape (a directory-only fix -
symlinks, no new Python; see the results doc). The CODEX layer is genuinely
new format-scope for this project, and this module is the ingester for it.

What CODEX data in this specific deposit actually is - confirmed by direct
inspection with ``tifffile``, not assumed from the README
--------------------------------------------------------------------------
* One ``*.ome.tif`` per tumour/region: a pyramidal, LZW-compressed,
  ``uint8`` multi-channel image (``imagecodecs`` is required to decode it -
  added to ``pyproject.toml``). For the file this adapter was built and
  tested against (``ZH1007_NEC_B.ome.tif``, Zenodo 10.5281/zenodo.12624860):
  43 channels, 7330 x 8205 px, ``PhysicalSizeX``/``PhysicalSizeY`` = 0.50670
  um/px (read from the OME-XML, not assumed), 3 pyramid levels.
* **The channel names embedded in the OME-XML are DYE names, not marker
  names** - confirmed by direct inspection of ``TiffFile.ome_metadata``:
  channel 0 is ``"DAPI"``; the other 42 cycle through ``"ATTO 550"``,
  ``"FITC"``, ``"Cy5"`` (the three non-DAPI CODEX fluorescence channels,
  reused every cycle for a different antibody). This is the standard way a
  CODEX/PhenoCycler acquisition pipeline writes its raw per-cycle image, and
  the marker identity of each (cycle, channel) pair lives in a SEPARATE
  experiment manifest (an Akoya "channel names" / cycle sheet) that this
  specific OME-TIFF does not embed and that this pass did not find published
  as a standalone file in the Zenodo deposit's file listing (the two
  QuPath project ``.zip`` archives, 1.8-2.2 GB each, very likely carry it -
  QuPath project files rename channels on import - but downloading either
  was out of this pass's disk budget; see the results doc). **This module
  therefore does NOT claim to know which of the 43 channels is which of the
  39 antibodies** ``results/gbm_3d_registration.md`` §CODEX-panel lists
  (read directly from the paper's own ``CODEX_antibody_validations.pdf``,
  Zenodo, 1.4 MB) - :data:`CODEX_MARKER_ALIAS_TO_GENE` is the DESIGNED panel,
  used only for a declared-panel capability check, never for a per-cell
  value. Every per-cell array this module returns for non-DAPI channels
  carries a ``channel_XX_<dye>`` label, never a marker name, so nothing
  downstream can mistake a dye channel for the antibody it happened to carry
  that cycle.
* One thing IS embedded and unambiguous: channel 0 is always DAPI (the
  universal nuclear counterstain, imaged every cycle for drift correction) -
  no manifest is needed to know that. This module uses channel 0 alone to
  find approximate nucleus positions, which is enough for the GEOMETRIC
  registration this dataset was fetched for (``register_serial_sections``
  only ever consumes (x, y) point clouds - it does not need marker identity
  at all).
* Also confirmed directly and worth recording plainly: this file's own
  embedded OME-XML ``UUID``/``TiffData`` blocks self-reference the filename
  ``"ZH1017_INF_B.ome.tif"`` throughout - NOT ``"ZH1007_NEC_B.ome.tif"``,
  the name Zenodo actually served it under (patient 1017 vs 1007, region INF
  vs NEC). This looks like a stale metadata artefact from an earlier
  pipeline/renaming step on the depositors' side, not a bug in this reader -
  flagged here exactly as found, the same class of real upstream
  data-serving quirk ``results/registration_3d_breast_poc.md`` already
  documented for the breast Xenium tarballs. It does not affect this
  module's output (channel 0 is still DAPI either way) but is recorded so a
  future reader is not surprised by it.

What "cell position" means here - an honest downgrade from real segmentation
------------------------------------------------------------------------------
This is NOT single-cell segmentation. There is no nucleus boundary, no
cell-body assignment, no doublet/overlap handling. :func:`extract_nucleus_positions`
Gaussian-smooths the DAPI channel and takes local intensity maxima
(``skimage.feature.peak_local_max``) as approximate nucleus centres - the
same class of "cheap, honest stand-in" ``registration_3d.tissue_boundary_points``
already uses for a tissue outline, applied here one level down (per-nucleus
rather than per-tissue-silhouette). It will both miss real nuclei (dim,
touching, or out-of-focus ones) and hallucinate some (noise peaks above
threshold) - no claim of single-cell-accurate counts or positions is made,
and none of this module's outputs should be read as a validated segmentation.
It is a real point cloud with real micron coordinates, dense enough and
distributed correctly enough across the tissue's actual footprint to drive
the SAME density/boundary-based registration method
``registration_3d.register_serial_sections`` already uses - which is exactly
what this dataset was fetched to test.
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

__all__ = [
    "CODEX_MARKER_ALIAS_TO_GENE",
    "CodexSection",
    "read_ome_metadata",
    "extract_nucleus_positions",
    "load_section",
]

#: The 39-antibody CODEX panel this dataset's paper validated, read directly
#: from ``CODEX_antibody_validations.pdf`` (Zenodo 10.5281/zenodo.12624860,
#: "Antibody validations" figure legend) - NOT read from any per-channel
#: metadata in the OME-TIFF (which carries none; see module docstring).
#: Mapped to the canonical HGNC gene symbol this codebase's marker lists
#: (``tmesim.data.ingest.ENDOTHELIAL_MARKER_SUPERSET`` /
#: ``tmesim.reconstruct.cell_state.HYPOXIA_MARKERS`` / ``MGMT_GENES``) use, so
#: :func:`tmesim.data.ingest.capability_report` can be asked "does the
#: DESIGNED panel cover this axis" the same way it already is for Xenium/
#: Visium gene panels. A small number of antibody clone names map to a
#: specific member of a gene family (CD3 -> CD3E, HLA-DR -> HLA-DRA, NeuN ->
#: RBFOX3) where the antibody does not itself specify which paralog/subunit;
#: that ambiguity is real and is not resolved here, only carried forward
#: explicitly in this comment rather than silently picking one.
CODEX_MARKER_ALIAS_TO_GENE: Dict[str, str] = {
    "CD45": "PTPRC",
    "CD3": "CD3E",          # T-cell receptor complex; clone does not specify subunit
    "CD69": "CD69",
    "CD90": "THY1",
    "CD11c": "ITGAX",
    "CD8": "CD8A",
    "CD4": "CD4",
    "HLA-DR": "HLA-DRA",    # HLA-DR is an alpha/beta heterodimer; alpha chain used
    "FN1": "FN1",
    "SOX10": "SOX10",
    "CD31": "PECAM1",
    "VIM": "VIM",
    "PLP1": "PLP1",
    "PDPN": "PDPN",
    "GFAP": "GFAP",
    "EGFR": "EGFR",
    "CD279": "PDCD1",       # PD-1
    "SOX4": "SOX4",
    "CD14": "CD14",
    "DCX": "DCX",
    "AQP4": "AQP4",
    "OLIG2": "OLIG2",
    "NeuN": "RBFOX3",
    "PDGFRA": "PDGFRA",
    "CA9": "CA9",
    "GLUT1": "SLC2A1",
    "BCAN": "BCAN",
    "CD44": "CD44",
    "p53": "TP53",
    "GAP43": "GAP43",
    "CD163": "CD163",
    "NDRG1": "NDRG1",
    "APOD": "APOD",
    "SOX2": "SOX2",
    "S100B": "S100B",
    "APOE": "APOE",
    "CHI3L1": "CHI3L1",
    "CD19": "CD19",
    "MAP2": "MAP2",
}

_SOURCE_LABEL = (
    "Greenwald, Darnell, Hoefflin et al. (Tirosh/Suva labs), Cell 187(10):"
    "2485-2501 (2024), PMID 38653236 -- CODEX (Akoya PhenoCycler-Fusion, "
    "40-marker spatial proteomics) OME-TIFF, Zenodo 10.5281/zenodo.12624860 "
    "-- adapted by tmesim.data.adapters.codex. Cell positions are "
    "APPROXIMATE (DAPI local-maxima, not a validated segmentation); "
    "per-channel marker IDENTITY is UNRESOLVED (see module docstring) - "
    "only the DESIGNED panel (from the paper's antibody-validation PDF) is "
    "carried, never a per-cell marker value."
)


def _decode_ome_channel_names(ome_xml: str) -> List[str]:
    """Every ``Channel Name=\"...\"`` attribute, in channel order.

    Deliberately a plain regex, not a full OME-XML parse (``xml.etree``
    would work too, but the document has no namespaced attributes that
    matter here and this keeps the dependency footprint at ``tifffile``
    alone). Returns dye names for this deposit's files (see module
    docstring) - callers must not treat the result as marker names.
    """
    return re.findall(r'<Channel[^>]*\bName="([^"]*)"', ome_xml)


@dataclass
class CodexSection:
    """One CODEX OME-TIFF, read to the point ``register_serial_sections``
    needs: real micron-scale nucleus positions, plus everything this module
    could and could not determine about channel/marker identity.

    ``nucleus_xy_um`` is the point cloud to feed to
    ``registration_3d.register_serial_sections`` as ``xy_ref``/``xy_mov``.
    ``channel_names`` are DYE names (see module docstring), one per raw
    OME-TIFF channel, index-aligned to ``channel_intensities`` columns.
    ``channel_intensities`` is ``(n_nuclei, n_channels)`` - the RAW
    (unsmoothed) intensity of every raw channel AT each detected nucleus
    position, kept for completeness/future use once a real channel manifest
    is available. This module deliberately does NOT build an
    ``ingest.SampleRecord`` for CODEX sections (that shape's ``expression``
    field means "a named marker's per-cell value", which nothing here can
    honestly provide yet - see the module docstring); ``load_section``
    returns this ``CodexSection`` instead, and
    :func:`declared_panel_capability_report` covers the one thing that
    genuinely can be asked of an unresolved-channel panel: does the
    DESIGNED panel cover a given axis at all.
    """

    path: str
    sample_id: str
    n_channels: int
    image_shape_px: Tuple[int, int]
    um_per_px: float
    channel_names: List[str]
    pyramid_level_used: int
    nucleus_xy_um: np.ndarray
    channel_intensities: np.ndarray
    detection_params: Dict
    filename_metadata_mismatch: Optional[str]

    @property
    def n_nuclei(self) -> int:
        return int(self.nucleus_xy_um.shape[0])

    def summary(self) -> Dict:
        return {
            "path": self.path,
            "sample_id": self.sample_id,
            "n_channels": self.n_channels,
            "image_shape_px": list(self.image_shape_px),
            "um_per_px": self.um_per_px,
            "n_dye_names": len(self.channel_names),
            "pyramid_level_used": self.pyramid_level_used,
            "n_nuclei_detected": self.n_nuclei,
            "detection_params": self.detection_params,
            "marker_identity_resolved": False,
            "filename_metadata_mismatch": self.filename_metadata_mismatch,
        }


def read_ome_metadata(path: str) -> Dict:
    """Read shape/pixel-size/channel-name metadata WITHOUT decoding pixels.

    Raises ``ValueError`` if the file carries no OME-XML (refuses rather
    than guessing a pixel size) or no ``PhysicalSizeX``.
    """
    import tifffile

    with tifffile.TiffFile(path) as tf:
        if not tf.is_ome or tf.ome_metadata is None:
            raise ValueError(
                f"{path}: no OME-XML metadata block found - cannot recover "
                "a pixel-to-micrometre scale, refusing to guess one."
            )
        xml = tf.ome_metadata
        series = tf.series[0]
        shape = tuple(int(v) for v in series.shape)  # (C, Y, X)
        axes = series.axes
        n_pages = len(tf.pages)
        n_levels = len(series.levels)
        level_shapes = [tuple(int(v) for v in lvl.shape) for lvl in series.levels]

    m = re.search(r'PhysicalSizeX="([0-9.eE+-]+)"', xml)
    if m is None:
        raise ValueError(
            f"{path}: OME-XML has no PhysicalSizeX - cannot convert pixel "
            "coordinates to micrometres, refusing to assume one."
        )
    um_per_px = float(m.group(1))
    unit_m = re.search(r'PhysicalSizeXUnit="([^"]*)"', xml)
    unit = unit_m.group(1) if unit_m else None
    if unit not in (None, "µm", "um", "micron", "micrometer"):
        raise ValueError(
            f"{path}: PhysicalSizeXUnit={unit!r} is not micrometres - "
            "refusing to assume a unit conversion this module was not "
            "built to handle."
        )

    channel_names = _decode_ome_channel_names(xml)[: shape[0] if axes[0] == "C" else 0]

    filename_mismatch = None
    fn_matches = set(re.findall(r'UUID FileName="([^"]+)"', xml))
    served_basename = os.path.basename(path)
    if fn_matches and served_basename not in fn_matches:
        filename_mismatch = (
            f"OME-XML TiffData/UUID blocks self-reference {sorted(fn_matches)}, "
            f"which does NOT include the filename this file was actually "
            f"served/saved under ({served_basename!r}). Recorded, not "
            "corrected - does not affect channel 0 = DAPI."
        )

    return {
        "path": os.path.abspath(path),
        "shape_cyx": shape,
        "axes": axes,
        "n_pages": n_pages,
        "n_pyramid_levels": n_levels,
        "level_shapes": level_shapes,
        "um_per_px": um_per_px,
        "channel_names": channel_names,
        "filename_metadata_mismatch": filename_mismatch,
    }


def extract_nucleus_positions(
    path: str,
    dapi_channel: int = 0,
    level: int = 0,
    gaussian_sigma_px: float = 1.5,
    min_distance_px: int = 6,
    threshold_abs: float = 25.0,
    max_nuclei: Optional[int] = 400_000,
) -> Tuple[np.ndarray, np.ndarray, Dict]:
    """Approximate nucleus centroids from the DAPI channel's local maxima.

    Parameters chosen and recorded, not hidden: ``gaussian_sigma_px=1.5`` and
    ``min_distance_px=6`` were picked by inspecting real output on
    ``ZH1007_NEC_B.ome.tif`` (0.507 um/px -> ~3 um smoothing, ~3 um minimum
    peak spacing - short of a real ~10-15 um nucleus diameter, deliberately
    permissive so touching/dim nuclei are not merged away entirely, at the
    cost of some over-counting; see module docstring for what this method
    does and does not claim). ``threshold_abs=25`` was chosen against this
    file's own measured background floor (min=4, mean=17.9 across the whole
    field) - comfortably above background, comfortably below the bright
    end (max=255).

    Returns ``(xy_um, intensities, params)``: ``xy_um`` is ``(n, 2)``,
    ``intensities`` is ``(n, n_channels)`` (every raw channel's RAW,
    unsmoothed pixel value at each detected position - DAPI included, in
    DAPI's own column; the Gaussian smoothing is used only to find peaks,
    not to report intensities), ``params`` records every argument used plus
    ``n_detected`` (before any ``max_nuclei`` subsample) and ``n_returned``.

    Raises ``ImportError`` if ``imagecodecs`` is not installed (the real
    file this module targets is LZW-compressed; tifffile raises a plain
    ``ValueError`` on decode without it, which this function re-raises with
    the fix attached rather than leaving a cryptic message).
    """
    import tifffile
    from scipy import ndimage
    from skimage.feature import peak_local_max

    meta = read_ome_metadata(path)
    n_channels = meta["shape_cyx"][0]
    if not (0 <= dapi_channel < n_channels):
        raise ValueError(
            f"dapi_channel={dapi_channel} out of range for {n_channels} "
            f"channels in {path}"
        )
    scale = 4 ** level  # this deposit's pyramid levels are power-of-4 downsamples
    um_per_px = meta["um_per_px"] * scale

    with tifffile.TiffFile(path) as tf:
        series = tf.series[0]
        try:
            level_arr = series.levels[level].asarray()
        except ValueError as exc:
            if "imagecodecs" in str(exc):
                raise ImportError(
                    f"{path} is LZW/JPEG-compressed and needs the "
                    "'imagecodecs' package to decode (pip install "
                    "imagecodecs; added to pyproject.toml's core "
                    f"dependencies for this reason). Original error: {exc}"
                ) from exc
            raise
    # level_arr is (C, Y, X) for this deposit's OME-TIFFs (axes == "CYX").
    dapi = level_arr[dapi_channel].astype(np.float32)
    smoothed = ndimage.gaussian_filter(dapi, sigma=gaussian_sigma_px)
    coords_yx = peak_local_max(
        smoothed, min_distance=min_distance_px, threshold_abs=threshold_abs,
        exclude_border=False,
    )
    n_before = int(coords_yx.shape[0])
    if max_nuclei is not None and n_before > max_nuclei:
        rng = np.random.default_rng(0)
        idx = rng.choice(n_before, size=max_nuclei, replace=False)
        idx.sort()
        coords_yx = coords_yx[idx]

    xy_um = np.column_stack([
        coords_yx[:, 1].astype(np.float64) * um_per_px,
        coords_yx[:, 0].astype(np.float64) * um_per_px,
    ])

    # Sample every channel's RAW (unsmoothed) intensity at each detected
    # position, for completeness - see CodexSection docstring for why this
    # is not exposed as "expression". Reuses the SAME decoded ``level_arr``
    # from above (one decode of the whole level, not one per channel -
    # decoding is the expensive step for a real LZW-compressed file).
    rows, cols = coords_yx[:, 0], coords_yx[:, 1]
    intensities = level_arr[:, rows, cols].astype(np.float32).T  # (n, n_channels)

    params = {
        "dapi_channel": dapi_channel, "level": level,
        "gaussian_sigma_px": gaussian_sigma_px, "min_distance_px": min_distance_px,
        "threshold_abs": threshold_abs, "um_per_px_at_level": um_per_px,
        "n_detected": n_before, "n_returned": int(coords_yx.shape[0]),
        "max_nuclei": max_nuclei,
    }
    return xy_um, intensities, params


def load_section(
    path: str,
    sample_id: Optional[str] = None,
    level: int = 0,
    **extract_kwargs,
) -> CodexSection:
    """End-to-end: read metadata, extract approximate nucleus positions.

    ``sample_id`` defaults to the file's basename (without extension).
    """
    meta = read_ome_metadata(path)
    sid = sample_id or os.path.splitext(os.path.basename(path))[0]
    xy_um, intensities, params = extract_nucleus_positions(path, level=level, **extract_kwargs)
    return CodexSection(
        path=meta["path"],
        sample_id=sid,
        n_channels=meta["shape_cyx"][0],
        image_shape_px=(meta["shape_cyx"][1], meta["shape_cyx"][2]),
        um_per_px=meta["um_per_px"],
        channel_names=meta["channel_names"],
        pyramid_level_used=level,
        nucleus_xy_um=xy_um,
        channel_intensities=intensities,
        detection_params=params,
        filename_metadata_mismatch=meta["filename_metadata_mismatch"],
    )


def declared_panel_capability_report(sample_id: str = "codex_section") -> Dict:
    """The capability report for the DESIGNED CODEX panel (39 antibodies,
    from the paper's own antibody-validation PDF), via the exact same
    ``ingest.capability_report`` every Xenium/Visium sample in this repo
    goes through - ``expression=None`` throughout, because no per-cell
    marker-identity-resolved value exists to measure a detection RATE from
    (see module docstring). Every rate-dependent verdict will read
    UNMEASURED; PRESENCE verdicts (does the panel cover CD31/PECAM1?
    hypoxia? MGMT?) are real and answer the task's own question honestly.
    """
    from .. import ingest as ing

    genes = sorted(set(CODEX_MARKER_ALIAS_TO_GENE.values()))
    return ing.capability_report(
        genes=genes, expression=None, platform="codex", sample_id=sample_id,
    )
