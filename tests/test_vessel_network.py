"""Regression tests for the DISCRETE vessel-network source field.

Guards the two properties the whole A/B study in `results/vessel_network_impact.md`
rests on: (1) the discrete field is EXACTLY zero away from a vessel, and (2) the
KDE default is unchanged by the additive `mode=` keyword. Everything here runs on
a small synthetic population — no data download, no transport solve.
"""
import numpy as np
import pytest

from tmesim.core import CellPopulation, CellType, Domain
from tmesim.reconstruct.real_fields import (
    VASCULAR_MODES,
    build_real_fields,
    vascular_density_from_cells,
)
from tmesim.reconstruct.vessel_network import (
    LINK_DISTANCE_UM,
    VESSEL_RADIUS_UM,
    discrete_vascular_field,
    link_distance_sweep,
    segment_vessels,
    segmentation_report,
)

DX = 12.0


def _population():
    """Three well-separated vessel-like chains of endothelial cells + tumour fill.

    Each chain is 4 cells at 15 um spacing (< the 20 um linking distance), and the
    chains are 400 um apart (>> linking distance), so the correct answer is
    unambiguous: 3 vessels of 4 cells each.
    """
    chains = [(100.0, 100.0), (500.0, 100.0), (100.0, 500.0)]
    ex, ey = [], []
    for x0, y0 in chains:
        for k in range(4):
            ex.append(x0 + 15.0 * k)
            ey.append(y0)
    rng = np.random.default_rng(0)
    tx = rng.uniform(20.0, 700.0, 300)
    ty = rng.uniform(20.0, 700.0, 300)
    x = np.concatenate([ex, tx])
    y = np.concatenate([ey, ty])
    ctype = np.concatenate([
        np.full(len(ex), int(CellType.ENDOTHELIAL), dtype=np.int8),
        np.full(tx.size, int(CellType.TUMOR), dtype=np.int8),
    ])
    cells = CellPopulation(x=x, y=y, cell_type=ctype,
                           state=np.zeros(x.size, dtype=np.int8),
                           sensitivity=np.ones(x.size))
    n = int(np.ceil(760.0 / DX))
    return Domain(nx=n, ny=n, dx=DX), cells


def test_segment_vessels_recovers_known_chains():
    domain, cells = _population()
    v = segment_vessels(cells, domain, link_distance_um=LINK_DISTANCE_UM)
    assert v["n_endothelial_cells"] == 12
    assert v["n_vessels"] == 3, "three separated 4-cell chains must give three vessels"
    assert sorted(v["cells_per_vessel"].tolist()) == [4, 4, 4]
    assert v["singleton_fraction"] == 0.0
    # A 4-cell chain at 15 um spacing spans 45 um and is maximally elongated.
    assert v["median_extent_multi_um"] == pytest.approx(45.0, abs=1e-6)
    assert v["median_aspect_multi"] > 10.0
    # Distances are computed from cell coordinates, so they are dx-invariant.
    coarse = Domain(nx=domain.nx // 2, ny=domain.ny // 2, dx=DX * 2)
    v2 = segment_vessels(cells, coarse, link_distance_um=LINK_DISTANCE_UM)
    assert v2["n_vessels"] == v["n_vessels"]
    assert v2["tumor_to_vessel"]["median"] == pytest.approx(
        v["tumor_to_vessel"]["median"], rel=1e-12)


def test_link_distance_controls_merging():
    domain, cells = _population()
    # Below the 15 um intra-chain spacing every cell is its own object...
    tight = segment_vessels(cells, domain, link_distance_um=10.0)
    assert tight["n_vessels"] == 12
    # ...and far above the 400 um inter-chain spacing everything merges.
    loose = segment_vessels(cells, domain, link_distance_um=500.0)
    assert loose["n_vessels"] == 1
    sweep = link_distance_sweep(cells, domain, link_distances=(10.0, 20.0, 30.0))
    counts = [r["n_vessels"] for r in sweep]
    assert counts == sorted(counts, reverse=True), "merging must be monotone in eps"


def test_discrete_field_is_exactly_zero_away_from_vessels():
    """THE load-bearing property: no Gaussian tail, no floor."""
    domain, cells = _population()
    v = segment_vessels(cells, domain, link_distance_um=LINK_DISTANCE_UM)
    f = discrete_vascular_field(domain, v, radius_um=VESSEL_RADIUS_UM)

    assert f.min() == 0.0
    assert set(np.unique(f)) <= {0.0, 1.0}, "'unit' amplitude must be binary"
    assert (f == 0).mean() > 0.9, "a 7.5 um radius must leave most of the grid empty"
    assert f[~np.asarray(domain.tissue_mask, dtype=bool)].sum() == 0.0

    # Every nonzero pixel is within radius of an endothelial cell, and vice versa.
    ny, nx = domain.shape
    gx = (np.arange(nx) + 0.5) * domain.dx
    gy = (np.arange(ny) + 0.5) * domain.dx
    GX, GY = np.meshgrid(gx, gy)
    xy = v["xy"]
    d = np.min(np.hypot(GX[..., None] - xy[:, 0], GY[..., None] - xy[:, 1]), axis=-1)
    assert np.all(f[d > VESSEL_RADIUS_UM] == 0.0), "source leaks beyond the radius"

    # Larger radius => strictly more source, and the field stays binary.
    big = discrete_vascular_field(domain, v, radius_um=30.0)
    assert (big > 0).sum() > (f > 0).sum()


def test_match_total_to_preserves_the_integral():
    """The dose control: same total drug entering the tissue, different geometry."""
    domain, cells = _population()
    v = segment_vessels(cells, domain, link_distance_um=LINK_DISTANCE_UM)
    kde = vascular_density_from_cells(domain, cells)
    matched = discrete_vascular_field(domain, v, match_total_to=kde)
    tissue = np.asarray(domain.tissue_mask, dtype=bool)
    assert matched[tissue].sum() == pytest.approx(kde[tissue].sum(), rel=1e-12)
    # Still exactly zero away from vessels — rescaling must not fill anything in.
    assert (matched == 0).mean() > 0.9
    assert matched.min() == 0.0


def test_kde_default_is_unchanged_by_the_new_mode_keyword():
    """The additive `mode=` keyword must not perturb any existing caller."""
    domain, cells = _population()
    explicit = vascular_density_from_cells(domain, cells, mode="kde")
    default = vascular_density_from_cells(domain, cells)
    np.testing.assert_array_equal(default, explicit)
    # The KDE is the smooth field the discrete one replaces: nonzero nearly
    # everywhere, which is exactly the defect vessel_network.py addresses.
    tissue = np.asarray(domain.tissue_mask, dtype=bool)
    disc = vascular_density_from_cells(domain, cells, mode="discrete")
    assert (default[tissue] > 0).mean() > (disc[tissue] > 0).mean()
    assert disc.min() == 0.0
    assert "kde" in VASCULAR_MODES and "discrete" in VASCULAR_MODES

    # build_real_fields threads the mode through and defaults to KDE.
    np.testing.assert_array_equal(
        build_real_fields(domain, cells).vascular_density, default)
    np.testing.assert_array_equal(
        build_real_fields(domain, cells, vascular_mode="discrete").vascular_density,
        disc)

    with pytest.raises(ValueError):
        vascular_density_from_cells(domain, cells, mode="nonsense")
    with pytest.raises(TypeError):
        vascular_density_from_cells(domain, cells, radius_um=5.0)  # kde takes none


def test_segmentation_report_scores_against_published_values():
    domain, cells = _population()
    v = segment_vessels(cells, domain, link_distance_um=LINK_DISTANCE_UM)
    rep = segmentation_report(v, domain)
    keys = {r["key"] for r in rep["rows"]}
    assert {"mvd_per_mm2", "vascular_area_fraction_pct"} <= keys
    for r in rep["rows"]:
        assert r["pmid"].isdigit(), "every published comparison must carry a PMID"
        assert r["verdict"] in {"PASS", "HIGH", "LOW", "n/a"}
    # The area fraction is measured off the field, so it must match it.
    f = discrete_vascular_field(domain, v)
    tissue = np.asarray(domain.tissue_mask, dtype=bool)
    assert rep["vessel_area_fraction_pct"] == pytest.approx(
        100.0 * (f > 0)[tissue].mean(), rel=1e-9)


def test_empty_endothelium_is_handled():
    domain, cells = _population()
    tumor_only = CellPopulation(
        x=cells.x, y=cells.y,
        cell_type=np.full(cells.n, int(CellType.TUMOR), dtype=np.int8),
        state=np.zeros(cells.n, dtype=np.int8),
        sensitivity=np.ones(cells.n),
    )
    v = segment_vessels(tumor_only, domain)
    assert v["n_vessels"] == 0 and v["n_endothelial_cells"] == 0
    f = discrete_vascular_field(domain, v)
    assert f.shape == domain.shape and f.sum() == 0.0
