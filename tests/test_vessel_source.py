"""Area conservation of the VESSEL SOURCE rasterisation (``geometry='area'``).

``results/fstate_mechanism.md`` attributed ~100 % of the ``F_state`` grid drift
(``cov/var = +1.057``) to the write side: the shipped binary centre-in-disc
painting deposits a WHOLE grid cell of vessel (area ``dx^2``) whenever a grid-cell
centre lands inside a 7.5 um disc, so the total painted vessel area is a function
of the mesh -- 4.73 %, non-monotone, over ``dx = 30...5 um`` on the real GBM
section. ``physics/oxygen.py`` documents and fixes exactly that defect on the
SINK (``_tumor_density``: count -> areal density); the source never received it.

``vessel_network.discrete_vascular_field(geometry='area')`` is that fix. These
tests guard the three things that make it a fix rather than a change:

1. it CONSERVES -- the painted area is identical across a 6x dx range, not merely
   convergent (``test_area_geometry_conserves_painted_area``);
2. it paints the SAME physical object -- the union of the same discs at the same
   radius, so the absolute value is the true disc area and overlaps are not
   double-counted (``test_area_geometry_measures_the_true_disc_area``,
   ``test_area_geometry_takes_the_union_not_the_sum``);
3. it is OPT-IN -- the ``'binary'`` default is bit-identical to the historical
   field and nothing routes through the new branch by accident
   (``test_binary_default_is_untouched``).

Everything runs on a small synthetic population: no data download, no PDE solve.
"""
import numpy as np
import pytest

from tmesim.core import CellPopulation, CellType, Domain
from tmesim.reconstruct.real_fields import vascular_density_from_cells
from tmesim.reconstruct.vessel_network import (
    LINK_DISTANCE_UM,
    SOURCE_SUBPIXEL_PITCH_UM,
    VESSEL_RADIUS_UM,
    _area_coverage_field,
    discrete_vascular_field,
    segment_vessels,
)

#: Every rung is an exact integer multiple of ``SOURCE_SUBPIXEL_PITCH_UM``, which
#: is the condition under which the sub-pixel lattice is literally the same point
#: set at every rung and conservation is EXACT rather than approximate. This
#: mirrors the production ladder of ``results/fstate_mechanism.md``.
DX_LADDER = (30.0, 24.0, 18.0, 16.0, 12.0, 10.0, 8.0, 6.0, 5.0)

#: Physical extent of the synthetic section (um). Divisible by every rung above,
#: so the compared grids cover exactly the same rectangle and the test measures
#: the painting rather than the domain's rounding.
EXTENT_UM = 720.0


def _cells():
    """Three vessel-like chains of endothelial cells, plus tumour fill.

    Chain 1 is isolated. Chains 2 and 3 are placed so that their discs OVERLAP
    (12 um apart, against a 15 um disc diameter), which is what forces the
    painting to take a union rather than a sum.
    """
    ex, ey = [], []
    for k in range(4):                      # isolated chain
        ex.append(100.0 + 15.0 * k)
        ey.append(100.0)
    for k in range(4):                      # overlapping pair of chains
        ex.append(400.0 + 15.0 * k)
        ey.append(400.0)
        ex.append(400.0 + 15.0 * k)
        ey.append(412.0)
    rng = np.random.default_rng(0)
    tx = rng.uniform(20.0, EXTENT_UM - 20.0, 300)
    ty = rng.uniform(20.0, EXTENT_UM - 20.0, 300)
    x = np.concatenate([ex, tx])
    y = np.concatenate([ey, ty])
    ctype = np.concatenate([
        np.full(len(ex), int(CellType.ENDOTHELIAL), dtype=np.int8),
        np.full(tx.size, int(CellType.TUMOR), dtype=np.int8),
    ])
    return CellPopulation(x=x, y=y, cell_type=ctype,
                          state=np.zeros(x.size, dtype=np.int8),
                          sensitivity=np.ones(x.size))


def _domain(dx: float) -> Domain:
    n = int(round(EXTENT_UM / dx))
    return Domain(nx=n, ny=n, dx=float(dx))


def _painted_area(field: np.ndarray, dx: float) -> float:
    """Tissue-integrated painted vessel area (um^2) -- the conserved quantity."""
    return float(np.asarray(field).sum()) * dx * dx


# --------------------------------------------------------------------------- #
# 1. Conservation -- the defect being fixed
# --------------------------------------------------------------------------- #
def test_binary_painting_does_not_conserve_painted_area():
    """The DEFECT, pinned. If this ever stops failing to conserve, the premise of
    ``results/vessel_source_conservation.md`` has changed and the report is stale.
    """
    cells = _cells()
    areas = []
    for dx in DX_LADDER:
        dom = _domain(dx)
        v = segment_vessels(cells, dom, link_distance_um=LINK_DISTANCE_UM)
        areas.append(_painted_area(discrete_vascular_field(dom, v), dx))
    span = (max(areas) - min(areas)) / float(np.mean(areas))
    assert span > 0.02, (
        "the binary centre-in-disc painting is expected NOT to conserve the "
        f"painted vessel area under dx (measured span {span:.4f})")


def test_area_geometry_conserves_painted_area():
    """THE acceptance test: the same physical vessels, painted on a 6x range of
    meshes, must deposit the SAME total vessel area."""
    cells = _cells()
    areas = []
    for dx in DX_LADDER:
        dom = _domain(dx)
        v = segment_vessels(cells, dom, link_distance_um=LINK_DISTANCE_UM)
        areas.append(_painted_area(
            discrete_vascular_field(dom, v, geometry="area"), dx))
    span = (max(areas) - min(areas)) / float(np.mean(areas))
    # Exact to floating point, because every rung is an integer multiple of the
    # sub-pixel pitch and the sub-lattice is therefore the same point set.
    assert span < 1e-9, (
        f"geometry='area' must conserve the painted area under dx; span={span:.3e}, "
        f"areas={areas}")


def test_area_geometry_is_far_better_off_the_pitch_lattice_too():
    """Conservation must not DEPEND on dx dividing the sub-pixel pitch exactly.

    Off-lattice rungs genuinely lose the exactness above -- the sub-lattices are
    different point sets -- so the guarantee weakens to a Riemann-sum error of
    order ``perimeter x pitch`` per disc, which on a 12-endothelial-cell toy is
    a couple of percent because there are too few discs for those errors to
    average out. (On the real section's 4,295 endothelial cells they do: the
    ``pitch = 1.0`` vs ``0.25 um`` totals agree to 0.01 %.) The claim that must
    hold on ANY configuration is the relative one: the fix must be far better
    than the defect on the very same rungs.

    The binary arm is genuinely catastrophic here -- at ``dx = 45 um`` a 7.5 um
    disc is 1/8.6 of a grid cell and this configuration paints EXACTLY ZERO
    vessel -- which is the same quantisation failure, in its limiting form, that
    ``results/fstate_mechanism.md`` measures as 4.73 % on the production ladder.
    """
    cells = _cells()
    rungs = (7.5, 11.25, 22.5, 45.0)
    spans = {}
    for geom in ("binary", "area"):
        areas = []
        for dx in rungs:
            dom = _domain(dx)
            v = segment_vessels(cells, dom, link_distance_um=LINK_DISTANCE_UM)
            areas.append(_painted_area(
                discrete_vascular_field(dom, v, geometry=geom), dx))
        spans[geom] = (max(areas) - min(areas)) / float(np.mean(areas))
    assert spans["area"] < 0.05, f"off-lattice span {spans['area']:.4f} too large"
    assert spans["area"] < spans["binary"] / 10.0, (
        f"the fix must dominate the defect off-lattice too: {spans}")


# --------------------------------------------------------------------------- #
# 2. It paints the SAME object -- absolute value, and overlaps
# --------------------------------------------------------------------------- #
def test_area_geometry_measures_the_true_disc_area():
    """A single isolated endothelial cell must paint pi*r^2, at every dx.

    This is what makes the fix a fix rather than a re-scaling: the conserved
    value is the physically correct one, not an arbitrary constant.
    """
    exact = np.pi * VESSEL_RADIUS_UM ** 2
    for dx in DX_LADDER:
        dom = _domain(dx)
        f = _area_coverage_field(dom, np.array([[301.3, 297.7]]), VESSEL_RADIUS_UM)
        assert _painted_area(f, dx) == pytest.approx(exact, rel=0.01), (
            f"dx={dx}: painted disc area is not pi*r^2")


def test_area_geometry_takes_the_union_not_the_sum():
    """Overlapping discs must contribute their UNION, or the field is not an area."""
    d = 5.0
    r = VESSEL_RADIUS_UM
    lens = 2 * r * r * np.arccos(d / (2 * r)) - 0.5 * d * np.sqrt(4 * r * r - d * d)
    exact_union = 2 * np.pi * r * r - lens
    dom = _domain(12.0)
    f = _area_coverage_field(dom, np.array([[300.0, 300.0], [300.0 + d, 300.0]]), r,
                             0.25)
    got = _painted_area(f, 12.0)
    assert got == pytest.approx(exact_union, rel=0.01)
    assert got < 0.8 * 2 * np.pi * r * r, "overlap was double-counted"


def test_area_field_is_a_fraction_and_stays_inside_the_radius():
    cells = _cells()
    dom = _domain(12.0)
    v = segment_vessels(cells, dom, link_distance_um=LINK_DISTANCE_UM)
    f = discrete_vascular_field(dom, v, geometry="area")
    assert f.min() == 0.0
    assert f.max() <= 1.0
    assert np.any((f > 0.0) & (f < 1.0)), (
        "the whole point is a PARTIAL-area weight; a field that is still 0/1 has "
        "not fixed anything")
    assert f[~np.asarray(dom.tissue_mask, dtype=bool)].sum() == 0.0

    # No source anywhere a disc cannot reach: the nearest endothelial cell to any
    # painted grid cell must be within radius + the grid cell's own half-diagonal.
    ny, nx = dom.shape
    gx = (np.arange(nx) + 0.5) * dom.dx
    gy = (np.arange(ny) + 0.5) * dom.dx
    GX, GY = np.meshgrid(gx, gy)
    xy = v["xy"]
    dist = np.min(np.hypot(GX[..., None] - xy[:, 0], GY[..., None] - xy[:, 1]),
                  axis=-1)
    reach = VESSEL_RADIUS_UM + dom.dx * np.sqrt(2.0) / 2.0
    assert np.all(f[dist > reach] == 0.0), "source leaks beyond the vessel radius"


def test_area_geometry_scales_with_radius():
    cells = _cells()
    dom = _domain(12.0)
    v = segment_vessels(cells, dom, link_distance_um=LINK_DISTANCE_UM)
    small = _painted_area(discrete_vascular_field(dom, v, geometry="area",
                                                  radius_um=5.0), 12.0)
    big = _painted_area(discrete_vascular_field(dom, v, geometry="area",
                                                radius_um=15.0), 12.0)
    assert big > small


def test_subpixel_pitch_refines_towards_the_exact_area():
    """The estimator's own error is a Riemann error and must shrink with pitch."""
    exact = np.pi * VESSEL_RADIUS_UM ** 2
    dom = _domain(12.0)
    xy = np.array([[301.3, 297.7]])
    errs = [abs(_painted_area(_area_coverage_field(dom, xy, VESSEL_RADIUS_UM, p),
                              12.0) - exact)
            for p in (4.0, 1.0, 0.25)]
    assert errs[2] < errs[0], f"finer sub-pixel pitch did not help: {errs}"
    assert errs[2] / exact < 2e-3


# --------------------------------------------------------------------------- #
# 3. It is OPT-IN -- nothing shipped moves
# --------------------------------------------------------------------------- #
def test_binary_default_is_untouched():
    """The shipped painting must remain binary, and must remain the default.

    ``vascular_density_from_cells`` (hence every reconstruction and every
    published result) routes through ``discrete_vascular_field``, so a silent
    flip of this default would silently move all of them.
    """
    cells = _cells()
    dom = _domain(12.0)
    v = segment_vessels(cells, dom, link_distance_um=LINK_DISTANCE_UM)
    default = discrete_vascular_field(dom, v)
    explicit = discrete_vascular_field(dom, v, geometry="binary")
    assert np.array_equal(default, explicit)
    assert set(np.unique(default)) <= {0.0, 1.0}
    # and the production entry point still lands on it
    prod = vascular_density_from_cells(dom, cells)
    assert np.array_equal(prod, default)
    # the new branch must be reachable ONLY by asking for it
    assert not np.array_equal(
        discrete_vascular_field(dom, v, geometry="area"), default)


def test_match_total_to_restores_the_shipped_total_exactly():
    """The distribution-only arm: same total source, different placement.

    This is the arm that separates "the fix changed the total" from "the fix
    changed where the source is", which is the question
    ``results/vessel_source_conservation.md`` turns on.
    """
    cells = _cells()
    dom = _domain(12.0)
    v = segment_vessels(cells, dom, link_distance_um=LINK_DISTANCE_UM)
    binary = discrete_vascular_field(dom, v)
    matched = discrete_vascular_field(dom, v, geometry="area",
                                      match_total_to=binary)
    tissue = np.asarray(dom.tissue_mask, dtype=bool)
    assert float(matched[tissue].sum()) == pytest.approx(
        float(binary[tissue].sum()), rel=1e-12)
    assert not np.array_equal(matched, binary), "the placement must still differ"


def test_area_geometry_rejects_the_caliber_amplitude():
    cells = _cells()
    dom = _domain(12.0)
    v = segment_vessels(cells, dom, link_distance_um=LINK_DISTANCE_UM)
    with pytest.raises(ValueError, match="not defined for geometry='area'"):
        discrete_vascular_field(dom, v, geometry="area",
                                amplitude="per_vessel_cells")
    with pytest.raises(ValueError, match="unknown geometry"):
        discrete_vascular_field(dom, v, geometry="nonsense")


def test_empty_population_paints_nothing():
    dom = _domain(12.0)
    empty = {"xy": np.zeros((0, 2)), "labels": np.zeros(0, dtype=int),
             "cells_per_vessel": np.zeros(0, dtype=int)}
    f = discrete_vascular_field(dom, empty, geometry="area")
    assert f.shape == dom.shape
    assert f.sum() == 0.0
    assert _area_coverage_field(dom, np.zeros((0, 2)), VESSEL_RADIUS_UM).sum() == 0.0


def test_subpixel_pitch_constant_divides_the_production_ladder():
    """The exactness argument in :data:`SOURCE_SUBPIXEL_PITCH_UM` is a claim about
    the pitch and the ladder together; check it rather than trust the docstring."""
    for dx in DX_LADDER:
        q = dx / SOURCE_SUBPIXEL_PITCH_UM
        assert abs(q - round(q)) < 1e-9, (
            f"dx={dx} is not an integer multiple of the sub-pixel pitch, so "
            "conservation on that rung is approximate, not exact")
