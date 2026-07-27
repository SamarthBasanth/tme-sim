"""End-to-end tests on the REAL microenvironment (no synthetic data).

These exercise the two load-bearing halves of the pipeline against a real,
open 10x Xenium GBM section reconstructed by
:func:`tmesim.reconstruct.real_pipeline.reconstruct_real`:

1. **Reconstruction** — ``reconstruct_real`` yields a well-formed
   :class:`~tmesim.core.Microenvironment` with real cells at real positions,
   shaped tissue fields, and realistic cell-type fractions (including a sane,
   single-digit-to-low-teens endothelial fraction).
2. **Simulation** — :func:`tmesim.sim.engine.simulate_response` runs on that
   real microenvironment and returns a valid, populated
   :class:`~tmesim.core.SimulationResult`, deterministically.
"""
from __future__ import annotations

import numpy as np
import pytest

from tmesim.core import (
    CellType,
    Microenvironment,
    PDEConfig,
    SimConfig,
    SimulationResult,
)
from tmesim.reconstruct.real_pipeline import reconstruction_report
from tmesim.sim.engine import simulate_response


def _fast_sim_config(drug, seed: int = 0) -> SimConfig:
    """A SimConfig with a deliberately short PDE window for a fast smoke test."""
    cfg = SimConfig(seed=seed)
    cfg.drug = drug
    cfg.pde = PDEConfig(dt=20.0, total_time=1200.0, convection=True, save_every=10)
    return cfg


# --------------------------------------------------------------------------- #
# Reconstruction: reconstruct_real yields a valid real microenvironment
# --------------------------------------------------------------------------- #
def test_reconstruct_real_is_valid_microenvironment(real_micro):
    assert isinstance(real_micro, Microenvironment)

    # Real cells at real positions.
    n = int(real_micro.cells.n)
    assert n > 0
    assert real_micro.cells.x.shape == (n,)
    assert real_micro.cells.y.shape == (n,)

    # Tumor cells present and a strict subset of all cells.
    ct = np.asarray(real_micro.cells.cell_type)
    n_tumor = int((ct == int(CellType.TUMOR)).sum())
    assert 0 < n_tumor <= n

    # Fields are shaped to the domain and finite.
    shape = real_micro.domain.shape
    f = real_micro.fields
    for name, arr in (("ecm_density", f.ecm_density),
                      ("vascular_density", f.vascular_density),
                      ("oxygen", f.oxygen)):
        assert arr is not None, name
        assert np.asarray(arr).shape == shape, name
        assert np.all(np.isfinite(np.asarray(arr))), name


def test_reconstruct_real_celltype_fractions_realistic(real_micro):
    report = reconstruction_report(real_micro)
    fr = report["celltype_fractions"]

    # Fractions form a proper distribution over the cell-type enum.
    assert all(0.0 <= v <= 1.0 for v in fr.values())
    assert sum(fr.values()) == pytest.approx(1.0, abs=1e-6)

    # This is a GBM tumor section: malignant cells dominate but do not saturate.
    assert 0.3 <= fr["TUMOR"] <= 0.95

    # Endothelium is genuinely reconstructed and sits in a biologically sane
    # single-digit-to-low-teens percentage range (not absent, not implausible).
    assert 0.02 <= fr["ENDOTHELIAL"] <= 0.20, fr

    # Immune + stromal compartments are present (whole-niche, not tumor-only).
    assert fr["IMMUNE"] > 0.0
    assert fr["STROMAL"] > 0.0


# --------------------------------------------------------------------------- #
# simulate_response on the REAL micro: valid, populated result
# --------------------------------------------------------------------------- #
def test_simulate_response_on_real_micro(real_micro, real_drug):
    cfg = _fast_sim_config(real_drug, seed=0)

    result = simulate_response(real_micro, cfg)

    assert isinstance(result, SimulationResult)
    # predicted_response is a probability.
    assert 0.0 <= float(result.predicted_response) <= 1.0
    # viable / kill fractions are complementary and in range.
    assert 0.0 <= float(result.viable_tumor_fraction) <= 1.0
    assert 0.0 <= float(result.kill_fraction) <= 1.0
    assert float(result.kill_fraction) == pytest.approx(
        1.0 - float(result.viable_tumor_fraction), abs=1e-6)
    # history + times populated and aligned.
    assert np.asarray(result.viable_history).ndim == 1
    assert len(result.viable_history) == len(result.times)
    assert len(result.times) >= 1
    # visualization fields populated and shaped to the real domain.
    assert isinstance(result.fields, dict) and len(result.fields) > 0
    for name, arr in result.fields.items():
        assert np.asarray(arr).shape == real_micro.domain.shape, name
    # extras carry the diagnostics the engine documents.
    assert 0.0 <= float(result.extras["bbb_penetrance"]) <= 1.0
    assert float(result.extras["n_tumor"]) > 0.0


# --------------------------------------------------------------------------- #
# simulate_response: determinism (same seed -> identical result)
# --------------------------------------------------------------------------- #
def test_simulate_response_deterministic(real_micro, real_drug):
    cfg = _fast_sim_config(real_drug, seed=7)

    # simulate_response works on a private copy, so the same shared micro is safe.
    r_a = simulate_response(real_micro, cfg)
    r_b = simulate_response(real_micro, cfg)

    assert float(r_a.predicted_response) == float(r_b.predicted_response)
    assert float(r_a.viable_tumor_fraction) == float(r_b.viable_tumor_fraction)
    np.testing.assert_array_equal(
        np.asarray(r_a.viable_history), np.asarray(r_b.viable_history))


# --------------------------------------------------------------------------- #
# Tissue mask: the GEOMETRIC mask is grid-independent (the legacy one is not)
# --------------------------------------------------------------------------- #
def test_geometric_tissue_mask_coverage_is_dx_independent():
    """The geometric tissue mask covers the same PHYSICAL fraction of the
    bounding box at every grid spacing, while the legacy occupancy mask
    ("grid cells containing >= 1 segmented cell") collapses as dx shrinks.

    This is the defect the geometric mask exists to fix: ``real_fields._kde_on_grid``
    multiplies the vascular field by the tissue mask, so a mask whose covered area
    depends on dx makes the total O2 SUPPLY — and therefore mean oxygen and the
    hypoxic fraction — a function of the mesh rather than of the tissue.
    """
    from tmesim.reconstruct.real_pipeline import _load_xenium_retry
    from tmesim.reconstruct.real_xenium import (
        TISSUE_DISC_RADIUS_UM,
        tissue_mask_from_cells,
    )

    data = _load_xenium_retry(None)
    xy = np.asarray(data["cell_xy_um"], dtype=np.float64)
    x = xy[:, 0] - xy[:, 0].min()
    y = xy[:, 1] - xy[:, 1].min()

    dxs = (30.0, 20.0, 12.0)
    geo, occ = [], []
    for dx in dxs:
        nx = int(np.ceil(x.max() / dx)) + 1
        ny = int(np.ceil(y.max() / dx)) + 1
        g = tissue_mask_from_cells(x, y, nx, ny, dx, mode="geometric",
                                   radius_um=TISSUE_DISC_RADIUS_UM)
        o = tissue_mask_from_cells(x, y, nx, ny, dx, mode="occupancy")
        assert g.shape == (ny, nx) and g.dtype == bool
        assert o.shape == (ny, nx) and o.dtype == bool
        # The geometric mask is a superset of the occupancy mask: every grid cell
        # holding a centroid is within one radius of that centroid.
        assert bool((o & ~g).sum() == 0)
        geo.append(float(g.mean()))
        occ.append(float(o.mean()))

    # GEOMETRIC: covered fraction invariant across a 2.5x range of dx.
    assert max(geo) / min(geo) < 1.05, f"geometric mask is dx-dependent: {geo}"

    # LEGACY OCCUPANCY: covered fraction collapses (this is the bug, asserted so
    # the test fails loudly if someone quietly reinstates it as the default).
    assert max(occ) / min(occ) > 2.0, f"occupancy mask unexpectedly stable: {occ}"

    # And the production default is the geometric one.
    from tmesim.reconstruct.real_xenium import TISSUE_MASK_MODE
    assert TISSUE_MASK_MODE == "geometric"


def test_reconstruct_real_uses_geometric_tissue_mask(real_micro):
    report = reconstruction_report(real_micro)
    assert report["tissue_mask_mode"] == "geometric"
    assert report["tissue_disc_radius_um"] == pytest.approx(30.0)
    # A geometric mask covers a substantial, physically meaningful fraction of the
    # section rather than a scattering of occupied tiles.
    assert 0.4 < report["tissue_coverage_fraction"] < 0.95


# --------------------------------------------------------------------------- #
# Sensitivity modes: 'uniform' is fair across drugs; 'determinant' is not
# --------------------------------------------------------------------------- #
def test_sensitivity_uniform_mode_is_identical_across_drugs(real_micro):
    """mode='uniform' (the default) gives every drug the same sensitivity, so it
    cancels out of a cross-drug ranking."""
    from tmesim.reconstruct.real_pipeline import apply_drug_sensitivity

    work = real_micro.copy() if hasattr(real_micro, "copy") else real_micro
    seen = []
    for name in ("doxorubicin", "cisplatin", "temozolomide"):
        info = apply_drug_sensitivity(work, name)          # default mode
        assert info["mode"] == "uniform"
        assert info["is_uniform"] and info["mean"] == 1.0
        seen.append(np.asarray(work.cells.sensitivity, dtype=np.float64).copy())
    for arr in seen[1:]:
        np.testing.assert_array_equal(seen[0], arr)
    np.testing.assert_allclose(seen[0], 1.0)


def test_sensitivity_determinant_mode_flags_unmeasured_drugs(real_micro):
    """mode='determinant' scales drugs that HAVE a panel determinant but leaves
    drugs that do not at 1.0 — so the latter must be flagged UPPER BOUND."""
    from tmesim.reconstruct.real_pipeline import apply_drug_sensitivity

    work = real_micro.copy() if hasattr(real_micro, "copy") else real_micro

    dox = apply_drug_sensitivity(work, "doxorubicin", mode="determinant")
    assert dox["determinant_gene"] == "TOP2A"
    assert dox["upper_bound"] is False
    assert dox["mean"] < 1.0

    tmz = apply_drug_sensitivity(work, "temozolomide", mode="determinant")
    assert tmz["determinant_gene"] is None
    assert tmz["upper_bound"] is True
    assert "UPPER BOUND" in tmz["comparability"]
    assert tmz["mean"] == 1.0

    # This is exactly the unfair comparison: same tissue, same physics, yet TMZ
    # keeps a ~1/0.585 ~ 1.7x sensitivity advantage for a purely technical reason.
    assert tmz["mean"] / dox["mean"] > 1.5

    with pytest.raises(ValueError):
        apply_drug_sensitivity(work, "doxorubicin", mode="nonsense")
