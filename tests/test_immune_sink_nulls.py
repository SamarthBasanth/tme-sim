"""Tests for the HONEST nulls in ``validation/immune_sink_impact.py``.

Why this file exists. ``results/immune_sink_impact.md`` used to publish a 2.0x
exposure spread "against a well-mixed model's exactly 1.0", scored against a
single label-permutation null. That null destroys the immune compartment's
clumping, which is the very quantity the structure ratio measures, so it
produces a sink field far LESS variable than the observed one and loses a
variance-driven contest almost regardless of whether any spatial biology
exists. The generator now computes spatially-aware nulls alongside it and leads
with the excess over those; these tests pin the three things that correction
rests on:

1. the clumpiness-matched permutation and the block shuffle really do sit ABOVE
   the label-permutation null (and carry a higher sink-field CV, which is why);
2. the CV diagnostic is emitted for the observed field and for every null, so
   the demotion of the label permutation is checkable rather than asserted;
3. the whole statistic is REFUSED on spot data, where no immune cell was ever
   localised.

SPEED. The real §5 run is ~110 full transport re-solves on a 255x406 grid and
takes over an hour, which is not a test. Everything here runs the IDENTICAL code
paths (``_spatial_context`` / ``null_family`` / ``spatial_diagnostics``) on a
small SYNTHETIC section with a deliberately clumped immune compartment, injected
through the context builder's ``micro`` / ``drug`` / ``pde`` arguments. The
synthetic section is a fixture, not a fit: nothing here is tuned to reproduce
the flagship numbers, only the ORDERING and the emission contract.
"""
from __future__ import annotations

import numpy as np
import pytest

from tmesim.core import (
    CellPopulation,
    CellState,
    CellType,
    Domain,
    DrugParams,
    Microenvironment,
    PDEConfig,
    TissueFields,
)
from tmesim.validation.immune_sink_impact import (
    BLOCK_UM,
    CLUMP_SIZE,
    NULL_FAMILIES,
    SPOT_JITTER_UM,
    SpotDataError,
    _draw_block_shuffle,
    _draw_clumpiness_matched,
    _draw_label_permutation,
    _ensure_tree,
    _spatial_context,
    honest_spatial_summary,
    null_family,
    reconstruction_provenance,
    require_measured_cell_positions,
    spatial_diagnostics,
)


# --------------------------------------------------------------------------- #
# Fixtures — a small synthetic section with a CLUMPED immune compartment
# --------------------------------------------------------------------------- #
@pytest.fixture(scope="module")
def synth_micro() -> Microenvironment:
    """4,000 tumour cells and 800 immune cells in 40 tight clumps.

    Metadata is stamped exactly as ``reconstruct/real_pipeline.py`` stamps a
    Xenium reconstruction, so the single-cell guard sees the same field it sees
    in production.
    """
    rng = np.random.default_rng(20260728)
    domain = Domain(nx=60, ny=54, dx=20.0)
    mask = domain.tissue_mask
    mask[:3, :] = False           # a real tissue edge, so depth is defined
    mask[-3:, :] = False
    mask[:, :3] = False
    mask[:, -3:] = False
    w, h = domain.extent_um

    fields = TissueFields(
        ecm_density=rng.random(domain.shape) * 0.6,
        vascular_density=rng.random(domain.shape) ** 2,
    )

    n_tumor = 4000
    tx = rng.uniform(80.0, w - 80.0, n_tumor)
    ty = rng.uniform(80.0, h - 80.0, n_tumor)

    n_clump, per_clump = 40, 20
    cx = rng.uniform(120.0, w - 120.0, n_clump)
    cy = rng.uniform(120.0, h - 120.0, n_clump)
    ix = (cx[:, None] + rng.normal(0.0, 22.0, (n_clump, per_clump))).ravel()
    iy = (cy[:, None] + rng.normal(0.0, 22.0, (n_clump, per_clump))).ravel()
    ix = np.clip(ix, 80.0, w - 80.0)
    iy = np.clip(iy, 80.0, h - 80.0)

    x = np.concatenate([tx, ix])
    y = np.concatenate([ty, iy])
    ct = np.array([int(CellType.TUMOR)] * tx.size + [int(CellType.IMMUNE)] * ix.size)
    cells = CellPopulation(
        x=x, y=y, cell_type=ct,
        state=np.full(x.size, int(CellState.VIABLE), dtype=np.int8),
        sensitivity=np.ones(x.size),
    )
    micro = Microenvironment(domain=domain, fields=fields, cells=cells)
    micro.metadata = {
        "assay": "synthetic-xenium-like",
        "real_pipeline": {"source": "xenium"},
    }
    return micro


@pytest.fixture(scope="module")
def synth_drug() -> DrugParams:
    return DrugParams(name="probe", diffusivity=200.0, clearance=2e-4,
                      binding=5e-4, ic50=1.0)


@pytest.fixture(scope="module")
def synth_pde() -> PDEConfig:
    # short, non-convective: this is a null-machinery test, not a physics test
    return PDEConfig(dt=10.0, total_time=1200.0, convection=False, save_every=60)


@pytest.fixture(scope="module")
def ctx(synth_micro, synth_drug, synth_pde):
    return _spatial_context(strength=10.0, micro=synth_micro, drug=synth_drug,
                            pde=synth_pde, verbose=False)


@pytest.fixture(scope="module")
def families(ctx):
    """Three draws per family — enough to order the families, cheap enough to run."""
    return {f: null_family(ctx, f, n_draws=3, verbose=False)
            for f in ("label_permutation", "clumpiness_matched", "block_shuffle")}


@pytest.fixture(scope="module")
def spot_micro(synth_micro) -> Microenvironment:
    """The same geometry, stamped as a Visium (spot) reconstruction."""
    micro = Microenvironment(
        domain=synth_micro.domain, fields=synth_micro.fields,
        cells=synth_micro.cells.copy(),
    )
    micro.metadata = {
        "source": "real_visium",
        "sample_id": "AT10-BRA-5-FO-2_1",
        "n_spots": 3212,
        "real_pipeline": {"source": "visium"},
    }
    return micro


# --------------------------------------------------------------------------- #
# 1. The primary nulls sit ABOVE the label permutation
# --------------------------------------------------------------------------- #
def test_primary_nulls_run_and_exceed_the_label_permutation(families):
    """The whole correction: an equally-clumped null is a much harder null.

    The label permutation moves the immune label one cell at a time, so the
    clumping is destroyed; the clumpiness-matched permutation moves it in clumps
    of ``CLUMP_SIZE`` measured neighbours and the block shuffle relocates whole
    blocks' worth of immune count. Both must therefore produce a LARGER
    structure ratio than the label permutation — that gap is exactly the part of
    the published 2.0x that was never immune geometry.
    """
    label = families["label_permutation"]
    for name in ("clumpiness_matched", "block_shuffle"):
        fam = families[name]
        assert fam["n_draws"] == 3
        assert np.isfinite(fam["sr_mean"]) and fam["sr_mean"] > 0.0
        assert fam["is_primary"] is True
        assert fam["sr_mean"] > label["sr_mean"], (
            f"{name} mean SR {fam['sr_mean']:.4f} does not exceed the "
            f"label-permutation null {label['sr_mean']:.4f}; the label "
            f"permutation is supposed to be the WEAKEST null, and if it is not "
            f"then the demotion in results/immune_sink_impact.md §5.3 is wrong"
        )
    assert label["is_primary"] is False


def test_every_null_keeps_immune_cells_on_measured_positions(ctx):
    """A null that invents positions is not a null, it is a different tissue."""
    _ensure_tree(ctx)
    measured = set(zip(ctx["x"].tolist(), ctx["y"].tolist()))
    for draw in (_draw_label_permutation, _draw_clumpiness_matched,
                 _draw_block_shuffle):
        rng = np.random.default_rng(7)
        px, py = draw(ctx, rng)
        assert px.size == ctx["n_immune"], (
            f"{draw.__name__} placed {px.size} cells, not {ctx['n_immune']}"
        )
        got = set(zip(px.tolist(), py.tolist()))
        assert got <= measured, (
            f"{draw.__name__} placed immune cells at positions no cell was "
            f"measured at"
        )


def test_block_shuffle_uses_the_declared_block_size(ctx):
    """The block shuffle must actually MOVE count between 400 um blocks."""
    rng = np.random.default_rng(3)
    px, py = _draw_block_shuffle(ctx, rng)
    key = (np.floor(px / BLOCK_UM).astype(int) * 100003
           + np.floor(py / BLOCK_UM).astype(int))
    obs = ctx["immune"]
    okey = (np.floor(ctx["x"][obs] / BLOCK_UM).astype(int) * 100003
            + np.floor(ctx["y"][obs] / BLOCK_UM).astype(int))
    # the per-block count VECTOR should be a rearrangement, not the identity
    a = np.bincount(np.unique(key, return_inverse=True)[1])
    b = np.bincount(np.unique(okey, return_inverse=True)[1])
    assert not (a.size == b.size and np.array_equal(a, b)), (
        "the block shuffle returned the observed block-count layout unchanged"
    )


# --------------------------------------------------------------------------- #
# 2. The CV diagnostic — why the label permutation was too weak
# --------------------------------------------------------------------------- #
def test_cv_of_the_sink_field_is_emitted_everywhere(ctx, families):
    """CV must be emitted for the observed field, per family and per draw."""
    assert "sink_cv" in ctx["observed"]
    assert np.isfinite(ctx["observed"]["sink_cv"])
    for name, fam in families.items():
        for key in ("cv_mean", "cv_min", "cv_max"):
            assert key in fam and np.isfinite(fam[key]), f"{name} missing {key}"
        for draw in fam["draws"]:
            assert np.isfinite(draw["sink_cv"])
        assert fam["cv_min"] <= fam["cv_mean"] <= fam["cv_max"]


def test_label_permutation_under_clumps_the_field_it_is_compared_against(
        ctx, families):
    """The measured reason the label permutation is only a floor check.

    It is not that it loses; it is that its sink field is markedly less variable
    than the observed one, so it was never being asked a fair question. The
    clumpiness-matched null closes most of that gap.
    """
    obs_cv = ctx["observed"]["sink_cv"]
    assert families["label_permutation"]["cv_mean"] < obs_cv
    assert (families["clumpiness_matched"]["cv_mean"]
            > families["label_permutation"]["cv_mean"])


# --------------------------------------------------------------------------- #
# 3. The floors, the depth confound and the honest excess
# --------------------------------------------------------------------------- #
def test_diagnostics_emit_measured_floors_and_the_edge_confound(ctx):
    diag = spatial_diagnostics(ctx, verbose=False)

    floor = diag["uniform_sink_floor"]
    assert np.isfinite(floor["structure_ratio"])
    assert np.isfinite(floor["a_p90_over_p10"])
    # a genuinely well-mixed sink has essentially no immune-density contrast
    assert floor["structure_ratio"] < ctx["observed"]["structure_ratio"]

    edge = diag["edge_distance"]
    assert np.isfinite(edge["stats"]["structure_ratio"])
    assert -1.0 <= edge["corr_phi_edge_pearson"] <= 1.0
    assert -1.0 <= edge["corr_phi_edge_spearman"] <= 1.0

    # the torus null is diagnosed, never scored: it loses cells off the mask
    assert 0.0 < diag["torus"]["fraction_in_tissue_mean"] <= 1.0

    # the CV sweep is the mechanism, so it must actually span CV
    cvs = [r["sink_cv"] for r in diag["cv_sweep"]]
    assert max(cvs) > min(cvs)


def test_honest_excess_is_the_headline_quantity(ctx, families):
    diag = spatial_diagnostics(ctx, verbose=False)
    hon = honest_spatial_summary({"observed": ctx["observed"]}, families, diag)
    for key in ("observed_structure_ratio", "observed_sink_cv",
                "honest_excess_min", "honest_excess_max",
                "well_mixed_sr_floor", "well_mixed_spread_floor",
                "edge_only_structure_ratio", "p_floor"):
        assert key in hon and np.isfinite(hon[key]), f"missing {key}"
    # the excess is the observed value divided by a null mean, so it must be
    # smaller than the raw ratio the report used to lead with
    assert hon["honest_excess_max"] < hon["observed_structure_ratio"]
    # the p-value floor must follow the number of draws ACTUALLY run, so a
    # short run can never be reported at a long run's floor
    n = families["clumpiness_matched"]["n_draws"]
    assert hon["p_floor"] == pytest.approx(1.0 / (1.0 + n))
    assert hon["p_primary_max"] >= hon["p_floor"]


# --------------------------------------------------------------------------- #
# 4. THE GUARD — the statistic is refused on spot data
# --------------------------------------------------------------------------- #
def test_provenance_is_read_from_metadata_not_from_the_path(synth_micro,
                                                            spot_micro):
    single = reconstruction_provenance(synth_micro)
    assert single["single_cell"] is True and single["spot_based"] is False

    spot = reconstruction_provenance(spot_micro)
    assert spot["spot_based"] is True and spot["single_cell"] is False
    assert spot["n_spots"] == 3212
    assert "_place_cells" in spot["reason"]
    assert f"{SPOT_JITTER_UM:g} um" in spot["reason"]


def test_spot_data_guard_fires(spot_micro, synth_micro):
    with pytest.raises(SpotDataError) as excinfo:
        require_measured_cell_positions(spot_micro)
    msg = str(excinfo.value)
    assert "INVALID" in msg and "5.6" in msg
    # and it does NOT fire on measured single-cell data
    assert require_measured_cell_positions(synth_micro)["single_cell"] is True


def test_the_context_builder_refuses_spot_data(spot_micro, synth_drug,
                                               synth_pde):
    """The refusal is in the shared path, so no null family can bypass it."""
    with pytest.raises(SpotDataError):
        _spatial_context(strength=10.0, micro=spot_micro, drug=synth_drug,
                         pde=synth_pde, verbose=False)

    # allow_invalid exists only to DEMONSTRATE the invalidity, and the result it
    # returns carries the invalidity with it
    forced = _spatial_context(strength=10.0, micro=spot_micro, drug=synth_drug,
                              pde=synth_pde, verbose=False, allow_invalid=True)
    assert forced["provenance"]["spot_based"] is True
    assert forced["provenance"]["single_cell"] is False


def test_unknown_provenance_is_refused_too():
    """Silence is not consent: an unstamped reconstruction is not single-cell."""
    domain = Domain(nx=8, ny=8, dx=20.0)
    cells = CellPopulation(
        x=np.array([10.0]), y=np.array([10.0]),
        cell_type=np.array([int(CellType.TUMOR)]),
        state=np.array([int(CellState.VIABLE)]), sensitivity=np.ones(1),
    )
    micro = Microenvironment(
        domain=domain,
        fields=TissueFields(ecm_density=domain.zeros(),
                            vascular_density=domain.zeros()),
        cells=cells,
    )
    assert reconstruction_provenance(micro)["single_cell"] is False
    with pytest.raises(SpotDataError):
        require_measured_cell_positions(micro)


# --------------------------------------------------------------------------- #
# 5. The families the report enumerates are the families the code can run
# --------------------------------------------------------------------------- #
def test_declared_null_families_are_all_runnable(ctx):
    for fam in NULL_FAMILIES:
        out = null_family(ctx, fam, n_draws=1, verbose=False)
        assert out["label"]
        assert np.isfinite(out["sr_mean"]) and np.isfinite(out["cv_mean"])
    with pytest.raises(ValueError):
        null_family(ctx, "not_a_null", n_draws=1, verbose=False)
    assert CLUMP_SIZE >= 2
