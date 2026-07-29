"""Tests for tmesim.data.adapters.codex — the CODEX (Akoya PhenoCycler)
OME-TIFF ingester built for the Greenwald/Tirosh 2024 GBM dataset
(results/gbm_3d_registration.md).

REAL-ONLY, like the rest of the suite: these tests read the actual
downloaded ``ZH1007_NEC_B.ome.tif`` (Zenodo 10.5281/zenodo.12624860, ~920 MB,
``data_cache/greenwald_gbm_codex/``) and the matching materialized Visium
directory (``data_cache/greenwald_gbm_visium/GBM_ZH1007nec/outs``), and skip
cleanly if either is not present in a given checkout. Slower, full-resolution
(level=0) extraction is exercised once; the rest use ``level=1`` (4x
downsampled) purely for test speed - both code paths are the same function,
just a different pyramid level.
"""
from __future__ import annotations

import os

import numpy as np
import pytest

from tmesim.data import ingest as ing
from tmesim.data.adapters import codex

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(_HERE, ".."))
CODEX_PATH = os.path.join(
    REPO, "data_cache", "greenwald_gbm_codex", "ZH1007_NEC_B.ome.tif")
VISIUM_OUTS = os.path.join(
    REPO, "data_cache", "greenwald_gbm_visium", "GBM_ZH1007nec", "outs")

# Fast knobs for tests that don't specifically need full resolution.
_FAST_KWARGS = dict(level=1, gaussian_sigma_px=0.8, min_distance_px=2)


def _skip_if_missing_codex():
    if not os.path.isfile(CODEX_PATH):
        pytest.skip(f"{CODEX_PATH} not present in this checkout")


def _skip_if_missing_visium():
    if not os.path.isdir(VISIUM_OUTS):
        pytest.skip(f"{VISIUM_OUTS} not present in this checkout")


# --------------------------------------------------------------------------- #
# read_ome_metadata(): no pixel decode, just header/OME-XML
# --------------------------------------------------------------------------- #
def test_read_ome_metadata_real_file():
    _skip_if_missing_codex()
    meta = codex.read_ome_metadata(CODEX_PATH)

    assert meta["shape_cyx"][0] == 43  # n_channels
    assert meta["axes"] == "CYX"
    assert meta["n_pyramid_levels"] == 3
    # Physical pixel size read directly from the OME-XML, not assumed.
    assert 0.1 < meta["um_per_px"] < 2.0
    assert meta["channel_names"][0] == "DAPI"
    assert len(meta["channel_names"]) == 43
    # Confirmed real quirk in this deposit's file (see module docstring):
    # the embedded OME-XML self-references a DIFFERENT filename than the one
    # this file is actually served/saved under.
    assert meta["filename_metadata_mismatch"] is not None
    assert "ZH1007_NEC_B.ome.tif" not in meta["filename_metadata_mismatch"].split(
        "which does NOT include")[0]


def test_read_ome_metadata_rejects_a_non_ome_tiff(tmp_path):
    import tifffile

    bad = tmp_path / "plain.tif"
    tifffile.imwrite(str(bad), np.zeros((8, 8), dtype=np.uint8))
    with pytest.raises(ValueError, match="OME-XML"):
        codex.read_ome_metadata(str(bad))


# --------------------------------------------------------------------------- #
# extract_nucleus_positions(): the actual pixel work
# --------------------------------------------------------------------------- #
def test_extract_nucleus_positions_fast_level():
    _skip_if_missing_codex()
    xy_um, intensities, params = codex.extract_nucleus_positions(
        CODEX_PATH, **_FAST_KWARGS)

    assert xy_um.ndim == 2 and xy_um.shape[1] == 2
    n = xy_um.shape[0]
    assert n > 5000, "expected a real, non-trivial number of detected nuclei"
    assert intensities.shape == (n, 43)
    assert params["n_detected"] == n
    assert params["level"] == 1

    # Positions must fall within the level's real physical field of view -
    # the whole point of doing this in micrometres, not raw pixels.
    meta = codex.read_ome_metadata(CODEX_PATH)
    max_x_um = meta["shape_cyx"][2] * meta["um_per_px"]
    max_y_um = meta["shape_cyx"][1] * meta["um_per_px"]
    assert xy_um[:, 0].min() >= 0.0
    assert xy_um[:, 1].min() >= 0.0
    assert xy_um[:, 0].max() <= max_x_um
    assert xy_um[:, 1].max() <= max_y_um


def test_extract_nucleus_positions_deterministic():
    _skip_if_missing_codex()
    xy1, _, _ = codex.extract_nucleus_positions(CODEX_PATH, **_FAST_KWARGS)
    xy2, _, _ = codex.extract_nucleus_positions(CODEX_PATH, **_FAST_KWARGS)
    assert xy1.shape == xy2.shape
    assert np.array_equal(xy1, xy2)


def test_extract_nucleus_positions_full_resolution_matches_module_finding():
    """The exact level=0 numbers results/gbm_3d_registration.md reports -
    pinned so a future dependency bump can't silently change the detector's
    real output without a test noticing."""
    _skip_if_missing_codex()
    xy_um, intensities, params = codex.extract_nucleus_positions(CODEX_PATH)
    assert params["level"] == 0
    assert params["n_detected"] == 40340
    assert xy_um.shape == (40340, 2)
    assert intensities.shape == (40340, 43)


def test_max_nuclei_subsamples_deterministically():
    _skip_if_missing_codex()
    xy_full, _, params_full = codex.extract_nucleus_positions(
        CODEX_PATH, max_nuclei=None, **_FAST_KWARGS)
    xy_sub, _, params_sub = codex.extract_nucleus_positions(
        CODEX_PATH, max_nuclei=1000, **_FAST_KWARGS)
    assert params_full["n_detected"] > 1000
    assert xy_sub.shape[0] == 1000
    assert params_sub["n_detected"] == params_full["n_detected"]  # same detector output
    # subsample rows must be an actual subset of the full detection
    full_set = {tuple(row) for row in xy_full}
    assert all(tuple(row) in full_set for row in xy_sub)


# --------------------------------------------------------------------------- #
# load_section(): the CodexSection wrapper
# --------------------------------------------------------------------------- #
def test_load_section_wraps_metadata_and_positions():
    _skip_if_missing_codex()
    sec = codex.load_section(CODEX_PATH, sample_id="test_codex", **_FAST_KWARGS)

    assert sec.sample_id == "test_codex"
    assert sec.n_channels == 43
    assert sec.n_nuclei == sec.nucleus_xy_um.shape[0] > 5000
    assert sec.channel_names[0] == "DAPI"
    assert sec.pyramid_level_used == 1

    summary = sec.summary()
    assert summary["marker_identity_resolved"] is False
    assert summary["n_nuclei_detected"] == sec.n_nuclei


# --------------------------------------------------------------------------- #
# CODEX_MARKER_ALIAS_TO_GENE + declared_panel_capability_report()
# --------------------------------------------------------------------------- #
def test_marker_alias_table_matches_antibody_pdf_panel():
    """39 markers, read directly from CODEX_antibody_validations.pdf
    (Zenodo 10.5281/zenodo.12624860) - see results/gbm_3d_registration.md."""
    assert len(codex.CODEX_MARKER_ALIAS_TO_GENE) == 39
    # A handful of the ones this task specifically asked about.
    assert codex.CODEX_MARKER_ALIAS_TO_GENE["CD31"] == "PECAM1"
    assert codex.CODEX_MARKER_ALIAS_TO_GENE["CA9"] == "CA9"
    assert codex.CODEX_MARKER_ALIAS_TO_GENE["GLUT1"] == "SLC2A1"
    assert codex.CODEX_MARKER_ALIAS_TO_GENE["NDRG1"] == "NDRG1"
    assert "MGMT" not in codex.CODEX_MARKER_ALIAS_TO_GENE.values()


def test_declared_panel_capability_report_confirms_vascular_and_hypoxia():
    cap = codex.declared_panel_capability_report("test_codex_capability")
    assert cap["platform"] == "codex"
    assert cap["supports"]["vascular_field"]["available"] is True
    assert "PECAM1" in cap["endothelial"]["present"]
    assert cap["supports"]["hypoxia_measured"]["available"] is True
    for m in ("CA9", "NDRG1", "SLC2A1"):
        assert m in cap["hypoxia"]["present"]
    # MGMT genuinely is not on this protein panel - honest negative.
    assert cap["supports"]["mgmt_from_panel"]["available"] is False
    assert cap["mgmt"]["present"] == []
    # No expression was supplied - every rate must read UNMEASURED, not a
    # fabricated number.
    assert cap["hypoxia"]["detection_rate"] is None


# --------------------------------------------------------------------------- #
# End-to-end: register the CODEX layer against its near-adjacent Visium spots
# --------------------------------------------------------------------------- #
def test_register_codex_against_visium_real_pair():
    _skip_if_missing_codex()
    _skip_if_missing_visium()
    from tmesim.reconstruct.registration_3d import register_serial_sections

    det = ing.detect(VISIUM_OUTS)
    assert det["platform"] == "visium"
    assert det["confidence"] == 1.0, det["issues"]
    vis_rec = ing.ingest(VISIUM_OUTS, sample_id="test_zh1007nec_visium", detection=det)
    assert vis_rec.n_cells > 500  # real spot count for this section

    sec = codex.load_section(CODEX_PATH, sample_id="test_zh1007nec_codex", **_FAST_KWARGS)

    caveat = (
        "test caveat: z-spacing between Visium and CODEX near-adjacent "
        "sections is unverified (see results/gbm_3d_registration.md)."
    )
    reg = register_serial_sections(
        vis_rec.cell_xy_um, sec.nucleus_xy_um, z_spacing_um=10.0,
        z_spacing_caveat=caveat,
        ref_sample_id="test_zh1007nec_visium", mov_sample_id="test_zh1007nec_codex",
        coarse_n_angles=36, coarse_top_k=3, icp_sample_size=3000,
        icp_max_iterations=25, seed=0,
    )

    assert reg.z_spacing_caveat == caveat  # NOT the breast-specific module default
    assert reg.points_xyz.shape == (vis_rec.n_cells + sec.n_nuclei, 3)
    assert set(np.unique(reg.section_id).tolist()) == {0, 1}
    assert (reg.points_xyz[reg.section_id == 1, 2] == 10.0).all()
    assert (reg.points_xyz[reg.section_id == 0, 2] == 0.0).all()

    # The actual claim: registration measurably improves overlap on this
    # real cross-platform (Visium spots vs CODEX nuclei) pair, by the same
    # honest before/after metrics the breast Xenium pass used. Not
    # hard-pinned to exact numbers (real download content could shift
    # slightly across a re-fetch) - just that the direction is real.
    q = reg.quality
    assert q.chamfer_um_after < q.chamfer_um_before
    assert q.hausdorff95_um_after < q.hausdorff95_um_before


def test_z_spacing_caveat_defaults_to_module_constant_when_not_passed():
    """Regression guard for the additive registration_3d.py change this task
    made: NOT passing z_spacing_caveat must still behave exactly as before
    (the breast-specific default), so the breast PoC's own behaviour is
    unchanged."""
    from tmesim.reconstruct.registration_3d import (
        Z_SPACING_CAVEAT,
        register_serial_sections,
    )

    rng = np.random.default_rng(0)
    xy_ref = rng.normal(size=(50, 2)) * 100.0
    xy_mov = rng.normal(size=(50, 2)) * 100.0
    reg = register_serial_sections(
        xy_ref, xy_mov, z_spacing_um=10.0,
        coarse_n_angles=8, coarse_top_k=1, icp_sample_size=50,
        icp_max_iterations=5,
    )
    assert reg.z_spacing_caveat == Z_SPACING_CAVEAT
