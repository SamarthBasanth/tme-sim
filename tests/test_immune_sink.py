"""Tests for ``tmesim.molecular.immune_sink`` — the TAM drug-sequestration sink.

The load-bearing test in this file is
:func:`test_off_path_is_bit_identical_to_the_shipped_solver`. The whole feature
is opt-in, and a caller who does not supply a sink field must get exactly the
numbers the shipped ``physics.diffusion.simulate_drug_transport`` gives — same
AUC array, same final concentration, same history, same penetration index, bit
for bit.

That test is ALSO the guard against drift. ``transport_with_immune_sink``
mirrors the shipped solver's time loop (it imports the solver's numerical
primitives unchanged but reproduces the stepping so it can carry a per-voxel
loss array). If anyone edits ``physics/diffusion.py``'s numerics, this test
fails loudly instead of the immune-sink report silently comparing against a
stale solver.

Every identity check is paired with a NON-VACUITY check, so it cannot pass by
the term being unreachable.
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
    SimConfig,
    TissueFields,
)
from tmesim.molecular.immune_sink import (
    IMMUNE_CELL_RADIUS_UM,
    IMMUNE_KERNEL_BANDWIDTH_UM,
    MILLER2015_DEPLETION_DROP,
    MILLER2015_TAM_DRUG_FRACTION,
    SEQUESTRATION_SWEEP,
    SUBSTRATE_TABLE,
    diverted_fraction_from_strength,
    field_summary,
    immune_area_fraction,
    immune_sink_detail,
    immune_sink_rate,
    immune_sink_strength_for,
    k_max_from_strength,
    strength_from_diverted_fraction,
    substrate_summary,
    transport_with_immune_sink,
)
from tmesim.physics.diffusion import simulate_drug_transport


# --------------------------------------------------------------------------- #
# Fixtures
# --------------------------------------------------------------------------- #
@pytest.fixture(scope="module")
def probe_micro() -> Microenvironment:
    """A small controlled grid with a clustered immune compartment."""
    rng = np.random.default_rng(20260728)
    domain = Domain(nx=31, ny=27, dx=9.0)
    mask = domain.tissue_mask
    mask[:, 0] = False  # exercise the non-tissue path
    fields = TissueFields(
        ecm_density=rng.random(domain.shape) * 0.6,
        vascular_density=rng.random(domain.shape),
    )
    # Two tight immune clusters plus one lone tumour cell far from both.
    ix = np.concatenate([rng.normal(60.0, 8.0, 40), rng.normal(200.0, 8.0, 40)])
    iy = np.concatenate([rng.normal(60.0, 8.0, 40), rng.normal(180.0, 8.0, 40)])
    tx = np.array([130.0, 60.0])
    ty = np.array([120.0, 60.0])
    x = np.concatenate([ix, tx])
    y = np.concatenate([iy, ty])
    ct = np.array([int(CellType.IMMUNE)] * ix.size + [int(CellType.TUMOR)] * 2)
    cells = CellPopulation(
        x=x, y=y, cell_type=ct,
        state=np.full(x.size, int(CellState.VIABLE), dtype=np.int8),
        sensitivity=np.ones(x.size),
    )
    return Microenvironment(domain=domain, fields=fields, cells=cells)


@pytest.fixture(scope="module")
def probe_drug() -> DrugParams:
    return DrugParams(name="probe", diffusivity=200.0, clearance=2e-4,
                      binding=5e-4, ic50=1.0)


# --------------------------------------------------------------------------- #
# THE load-bearing test: OFF is bit-identical
# --------------------------------------------------------------------------- #
def test_off_path_is_bit_identical_to_the_shipped_solver(probe_micro):
    """``sink_rate=None`` and ``sink_rate=zeros`` reproduce the shipped solver.

    Bit-for-bit, via ``.tobytes()`` and ``np.array_equal`` — never ``allclose``.
    Covers both loss parameterisations, convection on and off, two timesteps and
    the time-resolved-source branch.
    """
    rng = np.random.default_rng(7)
    source = np.asarray(probe_micro.fields.vascular_density) * 0.3
    velocity = (rng.random((2,) + probe_micro.domain.shape) - 0.5) * 0.02
    zeros = probe_micro.domain.zeros()

    n_checked = 0
    for clearance, binding in ((2e-4, 5e-4), (6.3e-5, 3.3e-2), (0.0, 0.0)):
        drug = DrugParams(name="probe", diffusivity=200.0,
                          clearance=clearance, binding=binding)
        for convection in (False, True):
            for dt, total in ((0.5, 60.0), (10.0, 1200.0)):
                pde = PDEConfig(dt=dt, total_time=total,
                                convection=convection, save_every=20)
                vel = velocity if convection else None
                ref = simulate_drug_transport(probe_micro, drug, source, pde,
                                              velocity=vel)
                for sink in (None, zeros):
                    got = transport_with_immune_sink(
                        probe_micro, drug, source, pde, velocity=vel,
                        sink_rate=sink)
                    assert got.auc.tobytes() == ref.auc.tobytes(), (
                        f"AUC drifted at clearance={clearance} binding={binding} "
                        f"convection={convection} dt={dt} sink={sink is None}"
                    )
                    assert (got.final_concentration.tobytes()
                            == ref.final_concentration.tobytes())
                    assert np.array_equal(got.times, ref.times)
                    assert np.array_equal(got.mean_concentration_history,
                                          ref.mean_concentration_history)
                    assert got.penetration_index == ref.penetration_index
                    n_checked += 1

                # time-resolved source branch
                prof = np.abs(rng.random(int(round(total / dt)) + 5))
                ref_t = simulate_drug_transport(probe_micro, drug, source, pde,
                                                velocity=vel, source_profile=prof)
                got_t = transport_with_immune_sink(
                    probe_micro, drug, source, pde, velocity=vel,
                    source_profile=prof, sink_rate=None)
                assert got_t.auc.tobytes() == ref_t.auc.tobytes()
                n_checked += 1

    assert n_checked == 36


def test_a_positive_sink_is_not_vacuous(probe_micro, probe_drug):
    """NON-VACUITY: a positive sink must move the answer, downward, unevenly."""
    pde = PDEConfig(dt=1.0, total_time=600.0, convection=False, save_every=100)
    source = np.asarray(probe_micro.fields.vascular_density) * 0.3
    ref = simulate_drug_transport(probe_micro, probe_drug, source, pde)

    phi = immune_area_fraction(probe_micro)
    sink = immune_sink_rate(probe_micro, k_max_from_strength(probe_drug, 10.0),
                            phi=phi)
    got = transport_with_immune_sink(probe_micro, probe_drug, source, pde,
                                     sink_rate=sink)

    assert not np.array_equal(got.auc, ref.auc), "the sink did nothing"
    assert np.all(got.auc <= ref.auc + 1e-15), "a LOSS term raised exposure"

    # ... and the attenuation is SPATIALLY VARYING, which is the whole point:
    # a well-mixed sink would give one number everywhere.
    tissue = np.asarray(probe_micro.domain.tissue_mask, dtype=bool)
    atten = (got.auc[tissue] + 1e-300) / (ref.auc[tissue] + 1e-300)
    assert float(np.ptp(atten)) > 1e-6, "the sink acted uniformly"


def test_sink_is_monotone_in_strength(probe_micro, probe_drug):
    """More sequestration must never mean more exposure."""
    pde = PDEConfig(dt=1.0, total_time=600.0, convection=False, save_every=100)
    source = np.asarray(probe_micro.fields.vascular_density) * 0.3
    tissue = np.asarray(probe_micro.domain.tissue_mask, dtype=bool)
    phi = immune_area_fraction(probe_micro)

    means = []
    for S in SEQUESTRATION_SWEEP:
        sink = None if S == 0.0 else immune_sink_rate(
            probe_micro, k_max_from_strength(probe_drug, S), phi=phi)
        tr = transport_with_immune_sink(probe_micro, probe_drug, source, pde,
                                        sink_rate=sink)
        means.append(float(tr.auc[tissue].mean()))

    assert all(b <= a + 1e-12 for a, b in zip(means, means[1:])), means
    assert means[-1] < means[0], "the largest swept strength changed nothing"


# --------------------------------------------------------------------------- #
# The field
# --------------------------------------------------------------------------- #
def test_immune_area_fraction_is_dimensionless_and_bounded(probe_micro):
    phi = immune_area_fraction(probe_micro)
    tissue = np.asarray(probe_micro.domain.tissue_mask, dtype=bool)
    assert phi.shape == probe_micro.domain.shape
    assert np.all(phi >= 0.0) and np.all(phi <= 1.0)
    assert np.all(phi[~tissue] == 0.0), "phi must be zero outside tissue"
    assert float(phi[tissue].max()) > 0.0, "no immune signal at all"
    # clustered input -> a right-skewed field (p99 well above the mean)
    s = field_summary(phi, tissue)
    assert s["p99"] > s["mean"]


def test_no_immune_cells_gives_an_exactly_zero_field_and_an_inert_sink(probe_micro):
    """A section with no immune cells must be EXACTLY the unmodified model."""
    empty = Microenvironment(
        domain=probe_micro.domain, fields=probe_micro.fields,
        cells=CellPopulation(
            x=np.array([100.0]), y=np.array([100.0]),
            cell_type=np.array([int(CellType.TUMOR)]),
            state=np.array([int(CellState.VIABLE)], dtype=np.int8),
            sensitivity=np.ones(1),
        ),
        patient_id=probe_micro.patient_id,
    )
    phi = immune_area_fraction(empty)
    assert np.array_equal(phi, empty.domain.zeros())

    drug = DrugParams(name="probe", diffusivity=200.0)
    pde = PDEConfig(dt=1.0, total_time=300.0, convection=False, save_every=100)
    source = np.asarray(empty.fields.vascular_density) * 0.3
    sink = immune_sink_rate(empty, k_max_from_strength(drug, 1e6), phi=phi)
    ref = simulate_drug_transport(empty, drug, source, pde)
    got = transport_with_immune_sink(empty, drug, source, pde, sink_rate=sink)
    assert got.auc.tobytes() == ref.auc.tobytes()


def test_positions_override_preserves_count_for_the_permutation_null(probe_micro):
    """The null's entry point: same number of cells, different places."""
    rng = np.random.default_rng(3)
    x = np.asarray(probe_micro.cells.x)
    y = np.asarray(probe_micro.cells.y)
    n_imm = int((np.asarray(probe_micro.cells.cell_type)
                 == int(CellType.IMMUNE)).sum())
    pick = rng.choice(x.size, size=n_imm, replace=False)

    phi_obs = immune_area_fraction(probe_micro)
    phi_perm = immune_area_fraction(probe_micro, positions=(x[pick], y[pick]))
    tissue = np.asarray(probe_micro.domain.tissue_mask, dtype=bool)

    assert not np.array_equal(phi_obs, phi_perm)
    # The TOTAL amount of immune material is conserved up to kernel truncation
    # at the domain edge; the SHAPE is what the null destroys.
    assert phi_perm.shape == phi_obs.shape
    assert float(phi_perm[tissue].max()) > 0.0


def test_sink_rate_validates_its_input(probe_micro, probe_drug):
    pde = PDEConfig(dt=1.0, total_time=60.0, convection=False, save_every=10)
    source = np.asarray(probe_micro.fields.vascular_density) * 0.3
    with pytest.raises(ValueError, match="non-negative"):
        transport_with_immune_sink(
            probe_micro, probe_drug, source, pde,
            sink_rate=-np.ones(probe_micro.domain.shape) * 1e-6)
    with pytest.raises(ValueError, match="shape"):
        transport_with_immune_sink(probe_micro, probe_drug, source, pde,
                                   sink_rate=np.zeros((3, 3)))
    with pytest.raises(ValueError, match="non-finite"):
        bad = probe_micro.domain.zeros()
        bad[0, 0] = np.nan
        transport_with_immune_sink(probe_micro, probe_drug, source, pde,
                                   sink_rate=bad)


# --------------------------------------------------------------------------- #
# The magnitude: swept, never fitted
# --------------------------------------------------------------------------- #
def test_every_drug_has_strength_exactly_zero():
    """The module's identity: nothing unproven gets a non-zero sink.

    No drug carries a MEASURED free-drug macrophage:tumour uptake ratio, which
    is the finding of the literature review, so the library-level strength must
    be exactly 0.0 for everything — including a name that is not in the table.
    """
    for name in list(SUBSTRATE_TABLE) + ["not_in_the_table", "IMATINIB"]:
        assert immune_sink_strength_for(name) == 0.0

    summary = substrate_summary()
    assert summary["n_with_measured_ratio"] == 0
    assert summary["n_drugs_assessed"] == len(SUBSTRATE_TABLE)


def test_detail_explains_which_reason_applies():
    for name in SUBSTRATE_TABLE:
        det = immune_sink_detail(name)
        assert det["strength"] == 0.0
        assert det["applied"] is False
        assert det["ratio"] is None
        assert det["sign"]
        assert isinstance(det["pmids"], list)
    unknown = immune_sink_detail("not_in_the_table")
    assert unknown["strength"] == 0.0
    assert "SUBSTRATE_TABLE" in unknown["reason"]

    # doxorubicin must carry BOTH the supporting mechanism and the evidence
    # against it — a one-sided entry would be the failure mode this guards.
    dox = immune_sink_detail("doxorubicin")
    assert "25544758" in dox["pmids"]     # lysosomal sequestration of weak bases
    assert "36105861" in dox["against"]   # free DOX uptake below detection
    assert "9067330" in dox["against"]    # TAM uptake does not explain it
    # cisplatin's only quantified TAM result runs in the DEPOT direction
    assert "depot" in immune_sink_detail("cisplatin")["sign"]
    # temozolomide has no claim at all
    assert immune_sink_detail("temozolomide")["pmids"] == []


def test_k_max_is_relative_to_the_drugs_own_sink():
    drug = DrugParams(name="probe", clearance=2e-4, binding=5e-4)
    base = drug.clearance + drug.binding
    assert k_max_from_strength(drug, 0.0) == 0.0
    assert k_max_from_strength(drug, 1.0) == base
    assert k_max_from_strength(drug, 7.5) == pytest.approx(7.5 * base)


def test_anchor_round_trips_and_matches_the_published_partition():
    """The Miller 2015 anchor conversion is exactly invertible."""
    for f in (MILLER2015_DEPLETION_DROP, MILLER2015_TAM_DRUG_FRACTION, 0.1):
        for phi in (0.01, 0.05, 0.2):
            S = strength_from_diverted_fraction(f, phi)
            assert diverted_fraction_from_strength(S, phi) == pytest.approx(f)
    assert strength_from_diverted_fraction(0.0, 0.05) == 0.0
    assert diverted_fraction_from_strength(0.0, 0.05) == 0.0
    # a section with no immune cells can divert nothing at any strength
    assert strength_from_diverted_fraction(0.3, 0.0) == 0.0
    with pytest.raises(ValueError):
        strength_from_diverted_fraction(1.0, 0.05)


def test_declared_constants_are_what_the_report_documents():
    """Pin the two declared constants and the literature-anchored ones."""
    assert IMMUNE_CELL_RADIUS_UM == 7.5
    assert IMMUNE_KERNEL_BANDWIDTH_UM == 25.0
    assert MILLER2015_TAM_DRUG_FRACTION == 0.30
    assert MILLER2015_DEPLETION_DROP == 0.24
    assert SEQUESTRATION_SWEEP[0] == 0.0, "the sweep must contain the OFF point"
    assert all(b > a for a, b in zip(SEQUESTRATION_SWEEP, SEQUESTRATION_SWEEP[1:]))


# --------------------------------------------------------------------------- #
# End-to-end on the REAL reconstruction
# --------------------------------------------------------------------------- #
def test_default_off_reproduces_the_engine_on_the_real_section(real_micro, real_drug):
    """The engine's kill fraction is untouched by this module existing.

    The strongest available statement: run ``simulate_response`` on the real
    reconstruction, then rebuild the same transport through
    ``transport_with_immune_sink`` with the sink OFF and require the AUC field to
    be bit-identical to the engine's own.
    """
    from tmesim.cells.agents import update_cell_states
    from tmesim.physics.oxygen import solve_oxygen
    from tmesim.physics.pressure import solve_ifp
    from tmesim.sim.engine import plasma_profile_safe, simulate_response
    from tmesim.treatment.drugs import delivery_source

    pde = PDEConfig(dt=10.0, total_time=600.0, convection=True, save_every=100)
    cfg = SimConfig(drug=real_drug, pde=pde, seed=0)
    res = simulate_response(real_micro, cfg)

    work = Microenvironment(
        domain=real_micro.domain,
        fields=TissueFields(ecm_density=real_micro.fields.ecm_density,
                            vascular_density=real_micro.fields.vascular_density),
        cells=real_micro.cells.copy(),
        patient_id=real_micro.patient_id, metadata=real_micro.metadata,
    )
    oxygen = solve_oxygen(work, cfg.oxygen)
    work.fields.oxygen = oxygen
    update_cell_states(work, oxygen, cfg.oxygen)
    ifp, velocity = solve_ifp(work)
    work.fields.ifp, work.fields.velocity = ifp, velocity
    times = np.arange(0.0, cfg.pde.total_time, cfg.pde.dt)
    plasma = plasma_profile_safe(cfg, times)
    source = delivery_source(work, real_drug, float(plasma.mean()))

    off = transport_with_immune_sink(work, real_drug, source, cfg.pde,
                                     velocity=velocity, sink_rate=None)
    assert off.auc.tobytes() == res.transport.auc.tobytes()

    # ... and a real sink DOES bite on the real section, in the right direction.
    phi = immune_area_fraction(real_micro)
    tissue = np.asarray(real_micro.domain.tissue_mask, dtype=bool)
    assert float(phi[tissue].mean()) > 0.0, "no immune cells on the real section"
    sink = immune_sink_rate(real_micro, k_max_from_strength(real_drug, 10.0),
                            phi=phi)
    on = transport_with_immune_sink(work, real_drug, source, cfg.pde,
                                    velocity=velocity, sink_rate=sink)
    assert float(on.auc[tissue].mean()) < float(off.auc[tissue].mean())
