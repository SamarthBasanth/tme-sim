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
from tmesim.reconstruct.vessel_geometry import (
    EC_AXIAL_PITCH_UM,
    MEAN_VESSEL_RADIUS_UM,
    SINGLE_EC_VESSEL_RADIUS_UM,
    subgrid_coverage_field,
    vascular_area_fraction,
    vessel_geometry,
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

#: The project's grid-independence convention (see results/numerical_rigor.md and
#: results/vessel_geometry_fix.md): a quantity that claims to describe the tissue
#: must not move by more than this factor across the production dx sweep.
GRID_INDEPENDENCE_TOL = 1.05
DX_SWEEP = (30.0, 20.0, 12.0, 8.0, 6.0)


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


def test_discrete_is_the_default_and_kde_is_still_reachable():
    """The default vascular construction is the DISCRETE segmented-vessel source
    (`results/vessel_promotion.md`), and `mode='kde'` still returns the historical
    smooth field unchanged."""
    domain, cells = _population()
    disc = vascular_density_from_cells(domain, cells, mode="discrete")
    default = vascular_density_from_cells(domain, cells)
    np.testing.assert_array_equal(default, disc)

    # The KDE remains available and is still the smooth field: nonzero nearly
    # everywhere, which is exactly the defect vessel_network.py addresses.
    kde = vascular_density_from_cells(domain, cells, mode="kde")
    tissue = np.asarray(domain.tissue_mask, dtype=bool)
    assert (kde[tissue] > 0).mean() > (disc[tissue] > 0).mean()
    assert disc.min() == 0.0
    assert "kde" in VASCULAR_MODES and "discrete" in VASCULAR_MODES

    # build_real_fields threads the mode through and follows the same default.
    np.testing.assert_array_equal(
        build_real_fields(domain, cells).vascular_density, disc)
    np.testing.assert_array_equal(
        build_real_fields(domain, cells, vascular_mode="kde").vascular_density,
        kde)

    with pytest.raises(ValueError):
        vascular_density_from_cells(domain, cells, mode="nonsense")
    with pytest.raises(TypeError):
        vascular_density_from_cells(domain, cells, mode="kde", radius_um=5.0)


def test_sparse_source_does_not_degenerate_the_vessel_mask():
    """A discrete source is exactly zero over most of the tissue, so a p95
    threshold on it is 0.0 and `vasc >= 0` would call the WHOLE tissue a vessel.
    `vessel_mask_from_density` must fall back to `field > 0` instead — while
    leaving the (strictly positive) KDE path on its percentile rule."""
    from tmesim.validation.vessel_profile import vessel_mask_from_density

    domain, cells = _population()
    tissue = np.asarray(domain.tissue_mask, dtype=bool)
    disc = vascular_density_from_cells(domain, cells, mode="discrete")
    kde = vascular_density_from_cells(domain, cells, mode="kde")

    vd = vessel_mask_from_density(disc, tissue)
    assert vd["threshold"] == 0.0
    assert "field > 0" in vd["rule"]
    assert vd["area_fraction"] < 0.5          # NOT the whole tissue
    np.testing.assert_array_equal(vd["mask"], (disc > 0) & tissue)

    vk = vessel_mask_from_density(kde, tissue)
    assert vk["threshold"] > 0.0
    assert "threshold on the field" in vk["rule"]
    np.testing.assert_array_equal(vk["mask"], (kde >= vk["threshold"]) & tissue)


def test_segmentation_report_scores_against_published_values():
    domain, cells = _population()
    v = segment_vessels(cells, domain, link_distance_um=LINK_DISTANCE_UM)
    rep = segmentation_report(v, domain)
    keys = {r["key"] for r in rep["rows"]}
    assert {"mvd_per_mm2", "vascular_area_fraction_pct"} <= keys
    for r in rep["rows"]:
        assert r["pmid"].isdigit(), "every published comparison must carry a PMID"
        assert r["verdict"] in {"PASS", "HIGH", "LOW", "n/a"}
    # The LEGACY number is still measured off the legacy binary field...
    f = discrete_vascular_field(domain, v)
    tissue = np.asarray(domain.tissue_mask, dtype=bool)
    assert rep["legacy_binary_area_fraction_pct"] == pytest.approx(
        100.0 * (f > 0)[tissue].mean(), rel=1e-9)
    # ...but the headline area fraction is now the sub-grid, data-derived one
    # (results/vessel_geometry_fix.md), which uses no hardcoded radius.
    assert rep["area_fraction_geometry"] == "subgrid"
    assert rep["vessel_area_fraction_pct"] == pytest.approx(
        vascular_area_fraction(domain, v)["area_fraction_pct"], rel=1e-9)
    assert rep["area_fraction_detail"]["radius_mode"] == "data"


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


# --------------------------------------------------------------------------- #
# The geometry fix: results/vessel_geometry_fix.md
# --------------------------------------------------------------------------- #
#: A physical extent that every dx in DX_SWEEP divides exactly, so that re-meshing
#: changes ONLY the mesh. (760 um does not divide by 30/12/6, which would make the
#: domain — and hence the denominator of the area fraction — a function of dx and
#: confound the very thing under test.)
SWEEP_EXTENT_UM = 720.0


def _population_at(dx):
    """The same physical section as `_population`, on a grid of spacing `dx`."""
    _, cells = _population()
    keep = (cells.x < SWEEP_EXTENT_UM) & (cells.y < SWEEP_EXTENT_UM)
    cells = CellPopulation(
        x=cells.x[keep], y=cells.y[keep], cell_type=cells.cell_type[keep],
        state=cells.state[keep], sensitivity=cells.sensitivity[keep])
    n = int(round(SWEEP_EXTENT_UM / dx))
    assert n * dx == SWEEP_EXTENT_UM
    return Domain(nx=n, ny=n, dx=dx), cells


def _isolated_vessels(dx, spacing=180.0, n_per_axis=4):
    """Single endothelial cells far enough apart to be one vessel each.

    Each is a k = 1 object, i.e. a vessel of radius r1 = 3.1 um and footprint
    ~30 um^2 — two orders of magnitude below a dx = 30 um grid cell (900 um^2).
    """
    g = np.arange(n_per_axis) * spacing + spacing / 2.0
    X, Y = np.meshgrid(g, g)
    x, y = X.ravel(), Y.ravel()
    cells = CellPopulation(
        x=x, y=y,
        cell_type=np.full(x.size, int(CellType.ENDOTHELIAL), dtype=np.int8),
        state=np.zeros(x.size, dtype=np.int8), sensitivity=np.ones(x.size))
    extent = spacing * n_per_axis
    n = int(round(extent / dx))
    assert n * dx == extent
    return Domain(nx=n, ny=n, dx=dx), cells, extent


def test_area_fraction_is_grid_independent():
    """THE acceptance test for the geometry fix.

    The vascular area fraction must be a property of the tissue, so re-meshing the
    SAME section must not move it. This pins the sub-grid, data-derived geometry
    across the production dx sweep; it is what stops a future change from
    re-introducing a mesh-imposed area fraction.
    """
    vals = []
    for dx in DX_SWEEP:
        domain, cells = _population_at(dx)
        v = segment_vessels(cells, domain, link_distance_um=LINK_DISTANCE_UM)
        vals.append(vascular_area_fraction(domain, v)["area_fraction_pct"])
    vals = np.asarray(vals)
    assert np.all(vals > 0)
    spread = float(vals.max() / vals.min())
    assert spread <= GRID_INDEPENDENCE_TOL, (
        f"area fraction moved {spread:.3f}x across dx={DX_SWEEP}: {vals}")
    # In fact the painter is grid-exact, not merely inside the convention: the
    # sub-lattice is anchored to physical coordinates, so the painted area does
    # not depend on dx at all.
    assert spread < 1.005, f"expected grid-exact, got {spread:.4f}x: {vals}"


def test_subgrid_paints_vessels_smaller_than_a_grid_cell():
    """The property the binary centre-in-disc rule could not have.

    A k = 1 vessel has radius 3.1 um and a footprint of ~30 um^2. At dx = 30 um a
    grid cell is 900 um^2, so the binary rule can only ever award it 0 or 900 um^2
    — a 30x error in either direction. The sub-grid painter must return its actual
    area at every dx.
    """
    domain, cells, extent = _isolated_vessels(30.0)
    v = segment_vessels(cells, domain, link_distance_um=LINK_DISTANCE_UM)
    geom = vessel_geometry(v)
    assert geom["n_vessels"] == cells.n
    assert np.allclose(geom["radius_um"], SINGLE_EC_VESSEL_RADIUS_UM)
    cell_area = domain.dx ** 2
    assert geom["footprint_area_um2"].max() < cell_area / 20.0, (
        "this test is only meaningful when every vessel is far sub-cell")

    cov = subgrid_coverage_field(domain, geom)
    assert cov.min() == 0.0, "must still be exactly zero away from a vessel"
    assert 0.0 < cov.max() < 1.0, "every cell must be only PARTLY vessel"
    # Every vessel is accounted for, even though none can fill a cell.
    analytic = float(geom["footprint_area_um2"].sum())
    painted = float((cov * cell_area).sum())
    assert painted == pytest.approx(analytic, rel=0.05)

    # The legacy binary field on the SAME vessels is quantised to whole cells: it
    # can award each vessel only 0 or 900 um^2 against a true 30 um^2, so whatever
    # the lattice phase happens to give, it is wrong by an order of magnitude.
    # (Here it awards 0 — none of the 16 vessels contains a grid-cell centre, so
    # the whole vasculature disappears from the source field.)
    binary = discrete_vascular_field(domain, v, radius_um=VESSEL_RADIUS_UM)
    binary_area = float((binary > 0).sum()) * cell_area
    assert set(np.unique(binary)) <= {0.0, 1.0}
    assert binary_area % cell_area == 0.0, "binary painting is quantised to cells"
    assert abs(binary_area - analytic) / analytic > 0.9
    assert abs(painted - analytic) / analytic < 0.05

    # ...and it is the quantisation, not the choice of radius, that does it: the
    # SAME 7.5 um radius painted sub-grid lands within 10% of its own exact area.
    g75 = vessel_geometry(v, radius_mode="constant", constant_radius_um=7.5)
    sub75 = float((subgrid_coverage_field(domain, g75) * cell_area).sum())
    assert sub75 == pytest.approx(float(g75["footprint_area_um2"].sum()), rel=0.10)


def test_vessel_radius_is_derived_not_assumed():
    """No constant in the radius may be a free choice.

    r1 is the measured human-cortex mean vessel radius divided by the measured
    mean endothelial cells per cross-section; the per-object radius is r1 times
    that object's own cells-per-cross-section. Nothing is tuned to 2.1 %.
    """
    # r_bar follows from VV = L_V . pi r^2 on the two H01 measurements alone.
    from tmesim.reconstruct.vessel_geometry import (
        CORTEX_VASCULAR_VOLUME_FRACTION,
        MEAN_EC_PER_CROSS_SECTION,
        VESSEL_LENGTH_DENSITY_MM_PER_MM3,
    )
    expect_rbar = 1e3 * np.sqrt(
        CORTEX_VASCULAR_VOLUME_FRACTION
        / (np.pi * VESSEL_LENGTH_DENSITY_MM_PER_MM3))
    assert MEAN_VESSEL_RADIUS_UM == pytest.approx(expect_rbar, rel=1e-12)
    assert SINGLE_EC_VESSEL_RADIUS_UM == pytest.approx(
        MEAN_VESSEL_RADIUS_UM / MEAN_EC_PER_CROSS_SECTION, rel=1e-12)
    assert EC_AXIAL_PITCH_UM == pytest.approx(50.0, rel=1e-12)

    # Each 4-cell chain spans 45 um < the 50 um axial pitch, so all four cells are
    # read as one cross-section: k = 4, r = 4 r1 -- straight out of the data.
    domain, cells = _population()
    v = segment_vessels(cells, domain, link_distance_um=LINK_DISTANCE_UM)
    geom = vessel_geometry(v)
    assert np.allclose(geom["ec_per_cross_section"], 4.0)
    assert np.allclose(geom["radius_um"], 4.0 * SINGLE_EC_VESSEL_RADIUS_UM)

    # Spreading the SAME four cells over 400 um (>> the axial pitch) makes them an
    # in-plane run instead, so k falls and the vessel gets thinner. The radius
    # therefore responds to the measured geometry, not to a constant.
    x = np.array([100.0, 233.0, 366.0, 500.0])
    y = np.full(4, 100.0)
    stretched = CellPopulation(
        x=x, y=y,
        cell_type=np.full(4, int(CellType.ENDOTHELIAL), dtype=np.int8),
        state=np.zeros(4, dtype=np.int8), sensitivity=np.ones(4))
    v2 = segment_vessels(stretched, domain, link_distance_um=200.0)
    g2 = vessel_geometry(v2)
    assert g2["n_vessels"] == 1
    assert g2["ec_per_cross_section"][0] < 1.1
    assert g2["radius_um"][0] < geom["radius_um"][0]
    assert g2["axis_length_um"][0] > 300.0, "an in-plane run keeps its length"


def test_footprint_is_one_per_vessel_object_not_one_per_endothelial_cell():
    """The defect that actually moved the number.

    The legacy field paints a disc around every endothelial NUCLEUS, so its area
    counts nuclei — a property of the marker panel and the typing threshold. Four
    nuclei lining one vessel and four nuclei in four separate vessels are the same
    picture to it. They must not be the same picture to the fixed geometry.
    """
    dx = 6.0
    # (a) four nuclei spread over four well-separated vessels.
    apart_dom, apart_cells, _ = _isolated_vessels(dx, spacing=180.0, n_per_axis=2)
    v_apart = segment_vessels(apart_cells, apart_dom,
                              link_distance_um=LINK_DISTANCE_UM)
    # (b) the same four nuclei lining ONE vessel cross-section.
    ang = np.arange(4) * np.pi / 2.0
    x = 180.0 + 6.0 * np.cos(ang)
    y = 180.0 + 6.0 * np.sin(ang)
    together = CellPopulation(
        x=x, y=y,
        cell_type=np.full(4, int(CellType.ENDOTHELIAL), dtype=np.int8),
        state=np.zeros(4, dtype=np.int8), sensitivity=np.ones(4))
    v_together = segment_vessels(together, apart_dom,
                                 link_distance_um=LINK_DISTANCE_UM)

    assert v_apart["n_vessels"] == 4 and v_together["n_vessels"] == 1
    assert (v_apart["n_endothelial_cells"]
            == v_together["n_endothelial_cells"] == 4)

    g_apart = vessel_geometry(v_apart)
    g_together = vessel_geometry(v_together)
    # ONE footprint per vessel object — the count follows the vasculature, not the
    # nucleus count. This is the invariant the per-nucleus painting did not have.
    assert g_apart["footprint_area_um2"].size == 4
    assert g_together["footprint_area_um2"].size == 1
    # Four capillaries are four k=1 vessels; one four-nucleus profile is a single
    # k=4 vessel, whose radius is 4 r1 by circumferential tiling.
    assert np.allclose(g_apart["ec_per_cross_section"], 1.0)
    assert np.allclose(g_together["ec_per_cross_section"], 4.0)
    assert np.allclose(g_together["radius_um"], 4.0 * SINGLE_EC_VESSEL_RADIUS_UM)

    # The legacy binary field cannot tell the two apart per nucleus: it paints the
    # same fixed 7.5 um disc around each of the four either way, so its painted
    # area per NUCLEUS is the same to within the overlap of the clustered case.
    per_nucleus = []
    for dom, vv in ((apart_dom, v_apart), (apart_dom, v_together)):
        f = discrete_vascular_field(dom, vv, radius_um=VESSEL_RADIUS_UM)
        per_nucleus.append(float((f > 0).sum()) * dx ** 2
                           / vv["n_endothelial_cells"])
    assert per_nucleus[1] == pytest.approx(per_nucleus[0], rel=0.5)


def test_subgrid_field_is_opt_in_and_default_is_unchanged():
    """The corrected geometry must not silently move any published result."""
    domain, cells = _population()
    v = segment_vessels(cells, domain, link_distance_um=LINK_DISTANCE_UM)
    default = discrete_vascular_field(domain, v)
    binary = discrete_vascular_field(domain, v, geometry="binary")
    np.testing.assert_array_equal(default, binary)
    assert set(np.unique(default)) <= {0.0, 1.0}

    sub = discrete_vascular_field(domain, v, geometry="subgrid")
    assert sub.min() == 0.0 and sub.max() <= 1.0
    assert np.any((sub > 0.0) & (sub < 1.0)), "sub-grid field must be fractional"
    assert not np.array_equal(sub, binary)
    assert sub[~np.asarray(domain.tissue_mask, dtype=bool)].sum() == 0.0

    with pytest.raises(ValueError):
        discrete_vascular_field(domain, v, geometry="nonsense")
    with pytest.raises(ValueError):
        vessel_geometry(v, radius_mode="nonsense")
    with pytest.raises(ValueError):
        vessel_geometry(v, radius_mode="constant")


def test_mvd_and_intervessel_distance_do_not_use_the_radius():
    """The downstream quantities that were flagged HIGH are radius-free.

    MVD and the inter-vessel distance come from cell COORDINATES and the tissue
    area, so correcting the vessel radius cannot move them and cannot explain
    their HIGH verdicts. Pinning that here stops the fix being credited with a
    change it did not make.
    """
    domain, cells = _population()
    v = segment_vessels(cells, domain, link_distance_um=LINK_DISTANCE_UM)
    base = segmentation_report(v, domain)
    for kw in ({"radius_um": 30.0},
               {"radius_mode": "stereological"},
               {"area_fraction_geometry": "binary"}):
        rep = segmentation_report(v, domain, **kw)
        for key in ("mvd_per_mm2", "intervessel_spacing_um",
                    "median_tumor_to_vessel_um"):
            a = next(r for r in base["rows"] if r["key"] == key)
            b = next(r for r in rep["rows"] if r["key"] == key)
            assert a["reconstructed"] == pytest.approx(b["reconstructed"], rel=1e-12)
            assert a["verdict"] == b["verdict"]


def test_empty_endothelium_is_handled_by_the_fixed_geometry():
    domain, cells = _population()
    tumor_only = CellPopulation(
        x=cells.x, y=cells.y,
        cell_type=np.full(cells.n, int(CellType.TUMOR), dtype=np.int8),
        state=np.zeros(cells.n, dtype=np.int8), sensitivity=np.ones(cells.n))
    v = segment_vessels(tumor_only, domain)
    geom = vessel_geometry(v)
    assert geom["n_vessels"] == 0 and geom["total_footprint_area_um2"] == 0.0
    assert subgrid_coverage_field(domain, geom).sum() == 0.0
    assert vascular_area_fraction(domain, v)["area_fraction_pct"] == 0.0
    assert discrete_vascular_field(domain, v, geometry="subgrid").sum() == 0.0
