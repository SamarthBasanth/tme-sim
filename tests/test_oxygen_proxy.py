"""Tests for the oxygen SURROGATE routes and for the default staying unchanged.

Two jobs, and the second is the more important one:

1. :mod:`tmesim.physics.oxygen_proxy` behaves — the mapping is the documented
   ``exp(-d/L)``, it is monotone, its calibration lands the hypoxia band edge
   exactly on the published onset distance, the marker route refuses honestly on a
   panel that cannot support it, and misalignment is a hard error.
2. **``reconstruct_real``'s default is BIT-IDENTICAL to the pre-option pipeline.**
   The whole point of adding ``oxygen_mode`` was to make an alternative available
   WITHOUT moving any published number, so the default path is asserted against a
   direct ``physics.solve_oxygen`` call with ``np.array_equal`` (not
   ``allclose``), and the historical ``metadata['real_pipeline']`` keys are
   asserted still present with their historical values.

Everything runs on the real cached 10x Xenium GBM section at a coarse grid, via
the session-scoped ``real_micro`` fixture where possible, so the suite stays fast.
"""
from __future__ import annotations

import numpy as np
import pytest

from tmesim.core import CellState, CellType, Microenvironment, OxygenConfig
from tmesim.physics.oxygen import solve_oxygen
from tmesim.physics.oxygen_proxy import (
    HYPOXIA_MARKER_GENES,
    HYPOXIA_ONSET_UM,
    OXYGEN_MODES,
    diffusion_length_um,
    hypoxia_from_markers,
    oxygen_from_vessel_distance,
    pde_diffusion_length_um,
    vessel_distance_um,
)
from tmesim.reconstruct.real_pipeline import (
    OXYGEN_MODE_DEFAULT,
    reconstruct_real,
    reconstruction_report,
)

from .conftest import REAL_DX_UM, REAL_SOURCE


# --------------------------------------------------------------------------- #
# The calibration is derived, not fitted
# --------------------------------------------------------------------------- #
def test_diffusion_length_puts_the_band_edge_on_the_published_onset():
    """L is defined so exp(-onset / L) == hypoxia_threshold, exactly."""
    t = OxygenConfig().hypoxia_threshold
    L = diffusion_length_um()
    assert L > 0.0
    assert np.exp(-HYPOXIA_ONSET_UM / L) == pytest.approx(t, abs=1e-12)
    # Documented value at the shipping defaults (150 um onset, 0.2 threshold).
    assert L == pytest.approx(93.2, abs=0.1)


def test_diffusion_length_rejects_impossible_anchors():
    for bad in (0.0, 1.0, -0.1, 1.5):
        with pytest.raises(ValueError):
            diffusion_length_um(HYPOXIA_ONSET_UM, bad)
    with pytest.raises(ValueError):
        diffusion_length_um(0.0, 0.2)


# --------------------------------------------------------------------------- #
# The surrogate field
# --------------------------------------------------------------------------- #
def test_surrogate_is_the_documented_exponential_of_distance(real_micro):
    field, info = oxygen_from_vessel_distance(real_micro, detail=True)
    dist = vessel_distance_um(real_micro)

    assert field.shape == real_micro.domain.shape == dist.shape
    assert np.isfinite(field).all()
    assert field.min() >= 0.0 and field.max() <= 1.0
    assert np.array_equal(field, np.exp(-dist / info["diffusion_length_um"]))
    assert info["hypoxic_beyond_um"] == pytest.approx(HYPOXIA_ONSET_UM)
    assert info["kind"].startswith("SURROGATE")


def test_surrogate_is_monotone_decreasing_in_distance(real_micro):
    """Monotonicity is load-bearing: it is why the mapping cannot be tuned to
    improve any rank statistic, which the report leans on."""
    field = oxygen_from_vessel_distance(real_micro)
    dist = vessel_distance_um(real_micro)
    order = np.argsort(dist.ravel(), kind="stable")
    assert np.all(np.diff(field.ravel()[order]) <= 1e-12)


def test_surrogate_is_deterministic(real_micro):
    a = oxygen_from_vessel_distance(real_micro)
    b = oxygen_from_vessel_distance(real_micro)
    assert np.array_equal(a, b)


def test_surrogate_bands_hypoxic_exactly_beyond_the_onset_distance(real_micro):
    """A grid cell is HYPOXIC under the shipping threshold iff it is further than
    the published onset distance from a vessel. That equivalence is the whole
    point of deriving L rather than picking it."""
    cfg = OxygenConfig()
    field = oxygen_from_vessel_distance(real_micro, cfg)
    dist = vessel_distance_um(real_micro)
    below = field < cfg.hypoxia_threshold
    beyond = dist > HYPOXIA_ONSET_UM
    assert np.array_equal(below, beyond)


def test_vessel_network_source_uses_a_subset_of_the_endothelium(real_micro):
    _, plain = oxygen_from_vessel_distance(real_micro, detail=True)
    _, seg = oxygen_from_vessel_distance(
        real_micro, source="vessel_network", detail=True
    )
    assert seg["n_points"] <= plain["n_points"]
    assert seg["n_vessels"] >= 1
    # Discarding points can only move a nearest-neighbour distance up.
    d_plain = vessel_distance_um(real_micro)
    d_seg = vessel_distance_um(real_micro, source="vessel_network")
    assert np.all(d_seg >= d_plain - 1e-9)


def test_unknown_vessel_source_raises(real_micro):
    with pytest.raises(ValueError, match="unknown vessel source"):
        vessel_distance_um(real_micro, source="nope")
    with pytest.raises(ValueError, match="unknown 'at'"):
        vessel_distance_um(real_micro, at="nope")


def test_no_endothelium_gives_a_zero_surrogate_not_a_crash(real_micro):
    """No vessel anywhere => no supply anywhere => maximally hypoxic. That is the
    physically honest answer, and it must not be a silent fallback to something
    else."""
    cells = real_micro.cells.copy()
    cells.cell_type = np.full_like(np.asarray(cells.cell_type), int(CellType.TUMOR))
    stripped = Microenvironment(
        domain=real_micro.domain, fields=real_micro.fields, cells=cells,
        patient_id="no-endothelium",
    )
    dist = vessel_distance_um(stripped)
    assert np.isinf(dist).all()
    field = oxygen_from_vessel_distance(stripped)
    assert np.array_equal(field, np.zeros(real_micro.domain.shape))


def test_pde_length_scale_is_reported_and_finite(real_micro):
    d = pde_diffusion_length_um(real_micro)
    assert d["mean_rho_tum"] > 0.0
    assert np.isfinite(d["L_mean_um"]) and d["L_mean_um"] > 0.0


# --------------------------------------------------------------------------- #
# The marker route refuses honestly where the panel cannot support it
# --------------------------------------------------------------------------- #
def test_marker_route_returns_none_with_a_reason_on_the_brain_panel(real_micro):
    """The 319-gene Xenium brain panel carries no VEGFA / CA9 / ADM / NDRG1 /
    SLC2A1, so hypoxia is UNMEASURED there — and saying so is the requirement."""
    from tmesim.data.real.singlecell_spatial import load

    data = load(dense=False)
    genes = np.asarray(data["genes"], dtype=str)
    assert not (set(HYPOXIA_MARKER_GENES) & set(genes.tolist()))

    assert hypoxia_from_markers(real_micro, data["expression"], genes) is None

    detail = hypoxia_from_markers(
        real_micro, data["expression"], genes, detail=True
    )
    assert detail["field"] is None
    assert detail["genes_used"] == []
    assert "UNMEASURED" in detail["reason"]
    for g in ("VEGFA", "CA9", "ADM", "NDRG1", "SLC2A1"):
        assert g in detail["genes_missing"]


def test_marker_route_builds_a_field_when_the_gene_is_present(real_micro):
    """Drive the success branch with a gene the brain panel DOES carry, so the
    branch is covered without needing the 527 MB IO section in the fast suite."""
    from tmesim.data.real.singlecell_spatial import load

    data = load(dense=False)
    genes = np.asarray(data["genes"], dtype=str)
    present = str(genes[0])

    field = hypoxia_from_markers(
        real_micro, data["expression"], genes, gene_set=(present,)
    )
    assert field is not None
    assert field.shape == real_micro.domain.shape
    assert np.isfinite(field).all()
    assert field.min() >= 0.0 and field.max() <= 1.0

    detail = hypoxia_from_markers(
        real_micro, data["expression"], genes, gene_set=(present,), detail=True
    )
    assert detail["genes_used"] == [present]
    assert detail["reason"] is None
    assert detail["per_cell"].shape == (real_micro.cells.n,)
    assert "not an oxygen tension" in detail["kind"]


def test_marker_route_rejects_misaligned_expression(real_micro):
    from tmesim.data.real.singlecell_spatial import load

    data = load(dense=False)
    genes = np.asarray(data["genes"], dtype=str)
    truncated = data["expression"][:-1]
    with pytest.raises(ValueError, match="not aligned"):
        hypoxia_from_markers(real_micro, truncated, genes, gene_set=(str(genes[0]),))


# --------------------------------------------------------------------------- #
# THE DEFAULT DID NOT CHANGE
# --------------------------------------------------------------------------- #
def test_default_oxygen_mode_is_pde():
    assert OXYGEN_MODE_DEFAULT == "pde"
    assert OXYGEN_MODES[0] == "pde"


def test_default_oxygen_mode_is_bit_identical(real_micro):
    """The load-bearing regression guard: with the default ``oxygen_mode``, the
    field ``reconstruct_real`` installs must be BIT-for-BIT the direct
    ``physics.solve_oxygen`` output — i.e. exactly the pre-option pipeline. Every
    published number in ``results/`` depends on this."""
    reference = solve_oxygen(real_micro, OxygenConfig())
    assert np.array_equal(reference, np.asarray(real_micro.fields.oxygen))

    explicit = reconstruct_real(source=REAL_SOURCE, dx_um=REAL_DX_UM,
                                oxygen_mode="pde")
    assert np.array_equal(np.asarray(explicit.fields.oxygen),
                          np.asarray(real_micro.fields.oxygen))


def test_default_reconstruction_metadata_is_unchanged(real_micro):
    """The historical metadata keys and values survive the additive change."""
    rp = real_micro.metadata["real_pipeline"]
    assert rp["source"] == "xenium"
    assert rp["vascular_source"] == "KDE of REAL ENDOTHELIAL cell positions"
    assert rp["oxygen_source"] == (
        "physics.solve_oxygen driven by the real vascular field"
    )
    assert rp["oxygen_mode"] == "pde"
    assert rp["oxygen_is_surrogate"] is False


def test_report_records_the_oxygen_mode(real_micro):
    rep = reconstruction_report(real_micro)
    assert rep["oxygen_mode"] == "pde"
    assert rep["oxygen"]["mode"] == "pde"
    assert rep["oxygen"]["is_surrogate"] is False
    assert "FALSIFIED AS VALUE-ADD" in rep["oxygen"]["caveat"]
    # The historical keys are all still there.
    for key in ("mean_oxygen", "mean_vascular_density", "mean_ecm_density",
                "celltype_fractions", "sensitivity"):
        assert key in rep


# --------------------------------------------------------------------------- #
# The opt-in modes, through the pipeline
# --------------------------------------------------------------------------- #
def test_vessel_distance_mode_installs_the_surrogate_and_says_so():
    micro = reconstruct_real(source=REAL_SOURCE, dx_um=REAL_DX_UM,
                             oxygen_mode="vessel_distance")
    expected = oxygen_from_vessel_distance(micro, OxygenConfig())
    assert np.array_equal(np.asarray(micro.fields.oxygen), expected)

    rep = reconstruction_report(micro)
    assert rep["oxygen_mode"] == "vessel_distance"
    assert rep["oxygen"]["is_surrogate"] is True
    assert "NOT a solved oxygen tension" in rep["oxygen"]["caveat"]
    assert rep["oxygen"]["detail"]["diffusion_length_um"] == pytest.approx(
        diffusion_length_um()
    )


def test_vessel_distance_mode_forwards_oxygen_kwargs():
    micro = reconstruct_real(
        source=REAL_SOURCE, dx_um=REAL_DX_UM, oxygen_mode="vessel_distance",
        oxygen_kwargs={"source": "vessel_network"},
    )
    rep = reconstruction_report(micro)
    assert rep["oxygen"]["detail"]["source"] == "vessel_network"


def test_markers_mode_raises_with_a_reason_on_the_brain_panel():
    with pytest.raises(ValueError, match="UNMEASURED"):
        reconstruct_real(source=REAL_SOURCE, dx_um=REAL_DX_UM,
                         oxygen_mode="markers")


def test_unknown_oxygen_mode_raises():
    with pytest.raises(ValueError, match="unknown oxygen_mode"):
        reconstruct_real(source=REAL_SOURCE, dx_um=REAL_DX_UM,
                         oxygen_mode="magic")


def test_pde_mode_takes_no_options():
    with pytest.raises(TypeError):
        reconstruct_real(source=REAL_SOURCE, dx_um=REAL_DX_UM,
                         oxygen_mode="pde", oxygen_kwargs={"source": "endothelial"})


# --------------------------------------------------------------------------- #
# The surrogate is usable by the frozen downstream code
# --------------------------------------------------------------------------- #
def test_surrogate_drives_the_unmodified_state_banding(real_micro):
    """``cells.agents.update_cell_states`` must accept the surrogate exactly as it
    accepts the solved field — the surrogate is a drop-in [0,1] field."""
    from tmesim.cells.agents import update_cell_states

    cfg = OxygenConfig()
    field = oxygen_from_vessel_distance(real_micro, cfg)
    work = Microenvironment(
        domain=real_micro.domain, fields=real_micro.fields,
        cells=real_micro.cells.copy(), patient_id="banding-test",
    )
    update_cell_states(work, field, cfg)
    states = np.asarray(work.cells.state)
    assert set(np.unique(states)).issubset(
        {int(CellState.HYPOXIC), int(CellState.QUIESCENT), int(CellState.VIABLE),
         int(CellState.DEAD)}
    )
    # Every band edge respects the field the banding was given.
    col, row = real_micro.domain.to_index(work.cells.x, work.cells.y)
    local = field[row, col]
    assert np.all(local[states == int(CellState.HYPOXIC)] < cfg.hypoxia_threshold)
    assert np.all(local[states == int(CellState.VIABLE)] >= 2 * cfg.hypoxia_threshold)
