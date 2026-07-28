"""Tests for tmesim.reconstruct.registration_3d — serial-section 3D registration.

REAL-ONLY where real data is available, like the rest of the suite: the
cross-section tests ingest the two actual Janesick et al. 2023 Xenium breast-
cancer replicate sections (via the EXISTING generic ingester,
``tmesim.data.ingest`` - no new ingestion code) and skip cleanly if that cache
is not present on disk (``data_cache/xenium_breast_janesick/{rep1,rep2}``).

The recoverability test uses REAL cell geometry (a subsample of Rep1's actual
centroids) put through a KNOWN, code-applied similarity transform, so the
"ground truth" is a controlled transform on real tissue shape, not fabricated
biology - the same pattern ``tests/conftest.py``'s ``radial_micro_factory``
uses for the physics solvers (controlled geometry, no invented outcome).
"""
from __future__ import annotations

import os

import numpy as np
import pytest

from tmesim.data import ingest as ing
from tmesim.reconstruct.registration_3d import (
    Z_SPACING_UM_ASSUMED,
    SimilarityTransform,
    chamfer_distance,
    density_correlation,
    icp_similarity,
    register_serial_sections,
    tissue_boundary_points,
)

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(_HERE, ".."))
REP1_DIR = os.path.join(REPO, "data_cache", "xenium_breast_janesick", "rep1")
REP2_DIR = os.path.join(REPO, "data_cache", "xenium_breast_janesick", "rep2")

# Reduced search/ICP knobs for test speed. The full-resolution defaults were
# checked by hand (results/registration_3d_breast_poc.md) to converge to the
# same near-identity answer; these settings reproduce that in well under a
# second instead of several seconds.
_FAST_KWARGS = dict(
    coarse_n_angles=36, coarse_top_k=3, icp_sample_size=3000,
    icp_max_iterations=25,
)


def _skip_if_missing():
    if not (os.path.isdir(REP1_DIR) and os.path.isdir(REP2_DIR)):
        pytest.skip(
            f"real breast-cancer serial-section cache not present at "
            f"{REP1_DIR} / {REP2_DIR} - see results/registration_3d_breast_poc.md "
            "for the download commands."
        )


@pytest.fixture(scope="module")
def rep1_xy() -> np.ndarray:
    _skip_if_missing()
    rec = ing.ingest(REP1_DIR, sample_id="xenium_breast_janesick_rep1")
    return rec.cell_xy_um


@pytest.fixture(scope="module")
def rep2_xy() -> np.ndarray:
    _skip_if_missing()
    rec = ing.ingest(REP2_DIR, sample_id="xenium_breast_janesick_rep2")
    return rec.cell_xy_um


# --------------------------------------------------------------------------- #
# Ingestion sanity (this is the "zero new ingestion code" claim, checked)
# --------------------------------------------------------------------------- #
def test_existing_ingester_reads_both_real_sections(rep1_xy, rep2_xy):
    """Both real Janesick Xenium replicates load via the pre-existing ingester.

    No new parsing code was written for this dataset - resolve_xenium_files /
    _ingest_xenium already handled it. This just pins down that the ingested
    arrays look like a real tissue section (matches ``validate_record``'s own
    floors), not degenerate output.
    """
    for xy in (rep1_xy, rep2_xy):
        assert xy.ndim == 2 and xy.shape[1] == 2
        assert xy.shape[0] > 50_000  # both replicates are >100k cells
        assert np.isfinite(xy).all()
        extent = xy.max(axis=0) - xy.min(axis=0)
        assert (extent > 1000.0).all()  # real tissue-scale extent in um


# --------------------------------------------------------------------------- #
# Recoverability: known synthetic transform on REAL geometry
# --------------------------------------------------------------------------- #
def test_icp_recovers_known_transform_given_a_reasonable_seed(rep1_xy):
    """Core correctness check on the ICP/Umeyama machinery: given a KNOWN
    rotation/scale/translation applied to real cell coordinates (so an exact
    point-to-point correspondence genuinely exists) and a seed reasonably
    close to the true answer, :func:`icp_similarity` converges to it almost
    exactly.

    Deliberately calls :func:`icp_similarity` directly with a close seed,
    NOT the full multi-start :func:`register_serial_sections` pipeline. That
    distinction matters and was found the hard way while building this
    module: on a dense, roughly-homogeneous cell point cloud (no
    distinguishing per-point features, fairly uniform nearest-neighbour
    spacing), trimmed point-to-point ICP has a NARROW basin of convergence -
    seeded even a few degrees off the true rotation on real Xenium cell
    density, it can converge to a self-consistent but WRONG local optimum
    rather than the true correspondence (confirmed empirically: seeding 3
    degrees off a true 37-degree synthetic rotation converged to -40 degrees,
    not -37). See results/registration_3d_breast_poc.md for the numbers and
    what this means for the method's scope. The coarse-search
    grid-plus-multi-start in the full pipeline mitigates this for the actual
    Janesick Rep1/Rep2 pair (whose true offset turned out to be close to
    identity), but this test isolates and verifies the piece that DOES have
    an exact, checkable right answer: local refinement from a good seed.
    """
    rng = np.random.default_rng(0)
    n = rep1_xy.shape[0]
    idx = rng.choice(n, size=min(20_000, n), replace=False)
    xy_ref = rep1_xy[idx]

    true_theta = np.radians(37.0)
    true_scale = 1.03
    true_R = np.array([[np.cos(true_theta), -np.sin(true_theta)],
                        [np.sin(true_theta), np.cos(true_theta)]])
    true_t = np.array([450.0, -300.0])
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        # Benign Accelerate/numpy matmul warning on this machine, not a real
        # numerical fault - see the matching comment in registration_3d.py.
        xy_mov = true_scale * (xy_ref @ true_R.T) + true_t

    # Seed near the TRUE inverse (rotation -37 deg, scale 1/1.03) - a stand-in
    # for "the coarse search landed close enough", which is exactly the
    # condition under which icp_similarity is expected to work.
    seed_theta = np.radians(-35.0)  # 2 degrees off the true inverse, not 0
    r0 = np.array([[np.cos(seed_theta), -np.sin(seed_theta)],
                   [np.sin(seed_theta), np.cos(seed_theta)]])
    s0 = 1.0 / true_scale
    t0 = xy_ref.mean(axis=0) - s0 * (r0 @ xy_mov.mean(axis=0))
    seed = SimilarityTransform(rotation=r0, scale=s0, translation=t0)

    transform, history = icp_similarity(
        xy_ref, xy_mov, init=seed, max_iterations=60, sample_size=8000,
        trim_fraction=0.95, seed=0, tol_um=1e-4,
    )
    recovered_xy = transform.apply(xy_mov)
    residual = np.linalg.norm(recovered_xy - xy_ref, axis=1)
    assert np.median(residual) < 1.0  # um - near-exact given a true correspondence
    assert transform.rotation_deg == pytest.approx(-37.0, abs=0.1)
    assert transform.scale == pytest.approx(1.0 / true_scale, abs=0.01)


# --------------------------------------------------------------------------- #
# Determinism
# --------------------------------------------------------------------------- #
def test_deterministic_given_fixed_seed(rep1_xy, rep2_xy):
    r1 = register_serial_sections(
        rep1_xy, rep2_xy, z_spacing_um=Z_SPACING_UM_ASSUMED,
        ref_sample_id="rep1", mov_sample_id="rep2", seed=0, **_FAST_KWARGS,
    )
    r2 = register_serial_sections(
        rep1_xy, rep2_xy, z_spacing_um=Z_SPACING_UM_ASSUMED,
        ref_sample_id="rep1", mov_sample_id="rep2", seed=0, **_FAST_KWARGS,
    )
    assert r1.transform.rotation_deg == pytest.approx(r2.transform.rotation_deg)
    assert r1.transform.scale == pytest.approx(r2.transform.scale)
    np.testing.assert_allclose(r1.transform.translation, r2.transform.translation)
    np.testing.assert_array_equal(r1.points_xyz, r2.points_xyz)
    assert r1.quality.to_dict() == r2.quality.to_dict()


def test_different_seed_still_converges_close(rep1_xy, rep2_xy):
    """The ICP subsample changes with the seed, but the answer should be
    stable (this section pair's true answer is close to identity - a highly
    seed-sensitive result would itself be a red flag worth surfacing)."""
    r_a = register_serial_sections(
        rep1_xy, rep2_xy, z_spacing_um=Z_SPACING_UM_ASSUMED,
        seed=0, **_FAST_KWARGS,
    )
    r_b = register_serial_sections(
        rep1_xy, rep2_xy, z_spacing_um=Z_SPACING_UM_ASSUMED,
        seed=1, **_FAST_KWARGS,
    )
    assert abs(r_a.transform.rotation_deg - r_b.transform.rotation_deg) < 2.0
    assert abs(r_a.transform.scale - r_b.transform.scale) < 0.02


# --------------------------------------------------------------------------- #
# Output shape / z-assignment contract
# --------------------------------------------------------------------------- #
def test_output_shape_and_z_assignment(rep1_xy, rep2_xy):
    z_spacing = 12.5
    result = register_serial_sections(
        rep1_xy, rep2_xy, z_spacing_um=z_spacing,
        ref_sample_id="rep1", mov_sample_id="rep2", seed=0, **_FAST_KWARGS,
    )
    n_ref, n_mov = rep1_xy.shape[0], rep2_xy.shape[0]
    assert result.points_xyz.shape == (n_ref + n_mov, 3)
    assert result.section_id.shape == (n_ref + n_mov,)
    assert result.cell_ids.shape == (n_ref + n_mov,)

    ref_rows = result.section_id == 0
    mov_rows = result.section_id == 1
    assert ref_rows.sum() == n_ref
    assert mov_rows.sum() == n_mov
    np.testing.assert_array_equal(result.points_xyz[ref_rows, 2], 0.0)
    np.testing.assert_array_equal(result.points_xyz[mov_rows, 2], z_spacing)

    # Reference section's (x, y) is untouched by registration.
    np.testing.assert_array_equal(result.points_xyz[ref_rows][:, :2], rep1_xy)

    assert result.z_spacing_um == z_spacing
    assert "UNVERIFIED" in result.z_spacing_caveat


# --------------------------------------------------------------------------- #
# The honesty check: does the quality metric actually improve, and is it
# measured the same way before and after?
# --------------------------------------------------------------------------- #
def test_registration_quality_improves_on_real_sections(rep1_xy, rep2_xy):
    """The core "is this doing something real" check from the task brief.

    Asserts the density-field correlation and BOTH boundary-outline metrics
    (Chamfer mean and Hausdorff-95th-percentile) improve (or at worst do not
    get worse) after registration relative to the raw, unregistered
    coordinate frames - using real Janesick et al. 2023 Xenium Rep1/Rep2
    cells. See results/registration_3d_breast_poc.md for the actual numbers
    and what they do/do not establish (this pair's true alignment turned out
    to be close to identity plus a translation - a small, real, honestly
    modest improvement, not a dramatic one).
    """
    result = register_serial_sections(
        rep1_xy, rep2_xy, z_spacing_um=Z_SPACING_UM_ASSUMED,
        ref_sample_id="rep1", mov_sample_id="rep2", seed=0, **_FAST_KWARGS,
    )
    q = result.quality
    assert q.density_correlation_after >= q.density_correlation_before - 1e-9
    assert q.chamfer_um_after <= q.chamfer_um_before + 1e-6
    assert q.hausdorff95_um_after <= q.hausdorff95_um_before + 1e-6
    # And the winning candidate must actually have been evaluated against
    # plausible alternatives, not accepted blind - the multi-start record
    # must contain more than one seed and the reported quality.after must
    # match the best (lowest-chamfer) candidate that was tried.
    assert len(result.candidates) >= 3
    best_cand_chamfer = min(c["final_chamfer_um"] for c in result.candidates)
    assert q.chamfer_um_after == pytest.approx(best_cand_chamfer)


# --------------------------------------------------------------------------- #
# Small unit-level checks on the building blocks
# --------------------------------------------------------------------------- #
def test_chamfer_distance_zero_for_identical_point_sets():
    rng = np.random.default_rng(0)
    pts = rng.normal(size=(200, 2)) * 100.0
    assert chamfer_distance(pts, pts) == pytest.approx(0.0, abs=1e-9)


def test_density_correlation_perfect_for_identical_clouds():
    rng = np.random.default_rng(0)
    pts = rng.normal(size=(5000, 2)) * 500.0
    assert density_correlation(pts, pts, bin_um=40.0) == pytest.approx(1.0, abs=1e-6)


def test_similarity_transform_roundtrip():
    rng = np.random.default_rng(0)
    pts = rng.normal(size=(50, 2)) * 200.0
    theta = np.radians(15.0)
    r = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    t = SimilarityTransform(rotation=r, scale=1.2, translation=np.array([10.0, -5.0]))
    out = t.apply(pts)
    assert out.shape == pts.shape
    assert t.rotation_deg == pytest.approx(15.0)
