"""Regression tests for the two fixes derived from measured human glioma PK
(`results/glioma_regional_pk.md`), plus the permeability ceiling surfaced by the
GIST imatinib run (`results/gist_imatinib.md`).

FIX 1  `molecular.permeability.bbb_penetrance` aggregates the four CNS-MPO
       desirabilities with a GEOMETRIC mean, so a disqualifying descriptor can
       VETO the molecule. An arithmetic mean let two favourable properties buy
       back two disqualifying ones (eribulin: predicted 0.334 vs measured
       0.0051, a 66x over-prediction and the worst miss in the validation).

FIX 2  `reconstruct.barrier` adds the measured-derived additive paracellular
       leak `clip(bbb_penetrance + L * bbb_disruption, 0, 1)`, L = 0.238. It is
       wired through `treatment.drugs.delivery_source` behind a default-off
       argument, and MUST be exactly backward compatible when off.

FIX 3  `molecular.permeability.membrane_permeability` is bounded by the
       unstirred aqueous boundary layer. It used to return 4.23e+07 um/s
       (42 m/s) for imatinib.
"""
import math

import numpy as np
import pytest

from tmesim.core import (
    CellPopulation, CellState, CellType, Domain, DrugParams, Microenvironment,
    TissueFields,
)
from tmesim.molecular import library, permeability
from tmesim.reconstruct import barrier
from tmesim.treatment.drugs import delivery_source


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _micro(nx=23, ny=17, dx=12.0, seed=0):
    domain = Domain(nx=nx, ny=ny, dx=dx)
    rng = np.random.default_rng(seed)
    fields = TissueFields(
        ecm_density=domain.ones() * 0.5,
        vascular_density=rng.random(domain.shape),
    )
    cells = CellPopulation(
        x=np.array([10.0]), y=np.array([10.0]),
        cell_type=np.array([CellType.TUMOR]), state=np.array([CellState.VIABLE]),
        sensitivity=np.array([0.5]),
    )
    return Microenvironment(domain=domain, fields=fields, cells=cells)


def _prefix_source(micro, drug, plasma):
    """The delivery source EXACTLY as it was written before FIX 2."""
    vascular = np.asarray(micro.fields.vascular_density, dtype=np.float64)
    return vascular * drug.bbb_penetrance * float(plasma)


# =========================================================================== #
# FIX 2 (i) -- default / zero disruption reproduces the prior source EXACTLY
# =========================================================================== #
@pytest.mark.parametrize("bbb", [0.0, 0.0051, 0.125, 0.4595, 0.84, 1.0])
@pytest.mark.parametrize("plasma", [0.0, 1.0, 37.3])
def test_default_disruption_reproduces_prior_delivery_source_bit_exactly(bbb, plasma):
    """Backward compatibility is BIT-exact, not approximate.

    Four routes must all return the identical array to the pre-fix expression:
    the default argument, an explicit scalar 0, an explicit zero field, and the
    barrier module's own helper on a microenvironment with no region info.
    """
    micro = _micro()
    drug = DrugParams(name="probe", bbb_penetrance=bbb)
    expected = _prefix_source(micro, drug, plasma)

    routes = {
        "default": delivery_source(micro, drug, plasma),
        "scalar_zero": delivery_source(micro, drug, plasma, bbb_disruption=0.0),
        "zero_field": delivery_source(
            micro, drug, plasma, bbb_disruption=np.zeros(micro.domain.shape)),
        "barrier_helper": barrier.barrier_delivery_source(micro, drug, plasma),
        "explicit_field": barrier.barrier_delivery_source(
            micro, drug, plasma, bbb_disruption=barrier.bbb_disruption_field(micro)),
    }
    for name, got in routes.items():
        assert got.shape == expected.shape, name
        assert got.dtype == expected.dtype, name
        # np.array_equal, not allclose: identical floating-point values.
        assert np.array_equal(got, expected), (
            f"{name} is not bit-identical to the pre-fix source "
            f"(max abs diff {np.abs(got - expected).max():.3e})"
        )


def test_bbb_disruption_field_is_zero_without_region_information():
    """A microenvironment with no region labels must NOT acquire a fabricated
    barrier defect."""
    micro = _micro()
    d = barrier.bbb_disruption_field(micro)
    assert d.shape == micro.domain.shape
    assert np.array_equal(d, np.zeros(micro.domain.shape))


# =========================================================================== #
# FIX 2 -- the leak itself
# =========================================================================== #
def test_leak_constant_is_the_mean_of_its_three_cited_fits():
    """L is documented as the mean of three independently fitted values; the
    docstring and the constant must not drift apart."""
    vals = list(barrier.FITTED_LEAK_VALUES.values())
    assert len(vals) == 3
    assert barrier.L_PARACELLULAR_LEAK == pytest.approx(float(np.mean(vals)), abs=5e-4)
    # drug-INDEPENDENCE is the whole reason the additive form was chosen
    assert max(vals) / min(vals) < 1.4


def test_effective_penetrance_is_additive_and_clipped():
    L = barrier.L_PARACELLULAR_LEAK
    assert float(barrier.effective_penetrance(0.25, 0.0)) == 0.25
    assert float(barrier.effective_penetrance(0.25, 1.0)) == pytest.approx(0.25 + L)
    assert float(barrier.effective_penetrance(0.25, 0.5)) == pytest.approx(0.25 + 0.5 * L)
    # never exceeds the free-drug ceiling
    assert float(barrier.effective_penetrance(0.95, 1.0)) == 1.0
    with pytest.raises(ValueError):
        barrier.effective_penetrance(0.25, 1.5)


def test_additive_leak_helps_an_excluded_drug_far_more_than_a_penetrant_one():
    """The signature of a size-nonselective paracellular opening, and the reason
    a MULTIPLICATIVE factor was rejected (3.2x spread across drugs vs 1.3x).

    Measured CE/non-enhancing: cefazolin (excluded) 4.26x, levetiracetam
    (penetrant) 1.45x — a 2.9x spread from ONE barrier state.
    """
    micro = _micro()
    disrupted = np.ones(micro.domain.shape)
    intact = np.zeros(micro.domain.shape)

    folds = {}
    for name, bbb in (("excluded", 0.0488), ("penetrant", 0.8399)):
        drug = DrugParams(name=name, bbb_penetrance=bbb)
        hi = delivery_source(micro, drug, 1.0, bbb_disruption=disrupted).sum()
        lo = delivery_source(micro, drug, 1.0, bbb_disruption=intact).sum()
        folds[name] = hi / lo

    assert folds["excluded"] > folds["penetrant"] > 1.0
    # right side of the measured anchors: excluded ~4-6x, penetrant ~1.1-1.5x
    assert 3.0 < folds["excluded"] < 8.0
    assert 1.05 < folds["penetrant"] < 1.6


def test_region_disruption_weights_and_field_assembly():
    """Perivascular (microvascular proliferation) is the contrast-enhancing
    anchor at 1.0; margin (leading edge) is the intact reference at 0.0;
    necrosis sits BETWEEN them, which is what was measured."""
    R = barrier.REGION_DISRUPTION
    assert R["perivascular"] == 1.0
    assert R["margin"] == 0.0
    assert 0.0 < R["necrotic"] < R["perivascular"]

    micro = _micro()
    ny, nx = micro.domain.shape
    peri = np.zeros((ny, nx)); peri[:, : nx // 2] = 1.0
    marg = np.zeros((ny, nx)); marg[:, nx // 2:] = 1.0
    d = barrier.bbb_disruption_field(micro, {"perivascular": peri, "margin": marg})
    assert np.allclose(d[:, : nx // 2], 1.0)
    assert np.allclose(d[:, nx // 2:], 0.0)
    assert d.min() >= 0.0 and d.max() <= 1.0

    # a 50/50 mixture interpolates
    mix = barrier.bbb_disruption_field(
        micro, {"perivascular": np.full((ny, nx), 0.5),
                "margin": np.full((ny, nx), 0.5)})
    assert np.allclose(mix, 0.5)

    # an unknown zone is a loud error, not a silent zero
    with pytest.raises(KeyError):
        barrier.bbb_disruption_field(micro, {"not_a_zone": peri})


def test_barrier_leak_only_acts_where_there_are_vessels():
    """A paracellular opening can only leak from a vessel. Avascular tissue must
    stay at zero however broken the barrier is."""
    micro = _micro()
    micro.fields.vascular_density = np.zeros(micro.domain.shape)
    drug = DrugParams(name="probe", bbb_penetrance=0.05)
    src = delivery_source(micro, drug, 10.0,
                          bbb_disruption=np.ones(micro.domain.shape))
    assert np.all(src == 0.0)


# =========================================================================== #
# FIX 1 -- the veto
# =========================================================================== #
def _arithmetic_reference(d):
    """The pre-fix aggregation, for the before/after assertions below."""
    mw, logp = float(d["mw"]), float(d["logp"])
    tpsa, hbd = float(d["tpsa"]), float(d["hbd"])
    ds = (1.0 / (1.0 + math.exp((mw - 450.0) / 60.0)),
          math.exp(-((logp - 2.5) ** 2) / (2.0 * 1.7 ** 2)),
          1.0 / (1.0 + math.exp((tpsa - 90.0) / 12.0)),
          1.0 / (1.0 + math.exp((hbd - 3.0) / 0.8)))
    return 1.0 / (1.0 + math.exp(-8.0 * (sum(ds) / 4.0 - 0.5)))


def test_disqualifying_mw_and_tpsa_veto_despite_favourable_logp_and_hbd():
    """THE eribulin case: MW 730 and TPSA 146 are disqualifying, logP 3.44 and
    2 HBD are favourable. A mean rescues it; a veto must not.

    Measured peritumoral brain:plasma is 0.0051 (Eroglu 2024, PMID 39422741).
    """
    eribulin = dict(mw=729.9, logp=3.44, tpsa=146.4, hbd=2, hba=12,
                    is_antibody=False)
    score = permeability.bbb_penetrance(eribulin)
    prefix = _arithmetic_reference(eribulin)

    assert prefix == pytest.approx(0.334, abs=0.01), "pre-fix reference drifted"
    assert score < 0.10, f"two disqualifying descriptors must veto; got {score:.4f}"
    assert score < prefix / 5.0, "the veto must be a large effect, not a nudge"
    # and it must land the right side of the measurement it was derived from
    assert score / 0.0051 < 10.0, "must be inside an order of magnitude of measured"


@pytest.mark.parametrize("bad", [
    dict(mw=900.0, tpsa=60.0),    # MW alone disqualifies
    dict(mw=300.0, tpsa=220.0),   # TPSA alone disqualifies
])
def test_a_single_disqualifying_descriptor_vetoes_a_favourable_molecule(bad):
    """One near-zero desirability must drive the score down even when every
    other descriptor is in the CNS sweet spot. This is the property an
    arithmetic mean structurally cannot have."""
    favourable = dict(mw=300.0, logp=2.5, tpsa=60.0, hbd=1, hba=4,
                      is_antibody=False)
    good = permeability.bbb_penetrance(favourable)
    assert good > 0.9, f"the control molecule should score high; got {good:.4f}"

    vetoed = permeability.bbb_penetrance({**favourable, **bad})
    assert vetoed < 0.35, f"{bad} should veto; got {vetoed:.4f}"
    assert vetoed < good / 3.0


def test_geometric_and_arithmetic_agree_on_a_uniform_molecule():
    """The change adds NO tunable parameter: for EQUAL desirabilities the
    geometric and arithmetic means coincide, so the veto bites ONLY on
    heterogeneous descriptor profiles.

    MW 450 / TPSA 90 / HBD 3 sit at their sigmoid centres (desirability 0.5
    each); logP 0.498 is the point where the Gaussian equals 0.5 as well
    (2.5 - 1.7*sqrt(2*ln2)).
    """
    logp_half = 2.5 - 1.7 * math.sqrt(2.0 * math.log(2.0))
    balanced = dict(mw=450.0, logp=logp_half, tpsa=90.0, hbd=3, hba=5,
                    is_antibody=False)
    assert permeability.bbb_penetrance(balanced) == pytest.approx(
        _arithmetic_reference(balanced), abs=1e-9)

    # ...and they DIVERGE on a heterogeneous one with the same arithmetic mean.
    hetero = dict(mw=900.0, logp=2.5, tpsa=30.0, hbd=0, hba=3, is_antibody=False)
    assert permeability.bbb_penetrance(hetero) < _arithmetic_reference(hetero) / 2.0


def test_veto_preserves_the_documented_orderings():
    """The fix must not buy eribulin at the cost of the orderings the function
    is actually used for."""
    bev = library.get_molecule("bevacizumab").descriptors
    assert permeability.bbb_penetrance(bev) == pytest.approx(0.0, abs=1e-9)

    scores = {n: permeability.bbb_penetrance(library.get_molecule(n).descriptors)
              for n in ("temozolomide", "lomustine", "carmustine", "methotrexate",
                        "doxorubicin", "etoposide", "vincristine")}
    assert all(0.0 <= v <= 1.0 for v in scores.values())

    # Small lipophilic nitrosoureas beat every large polar drug, by a wide margin.
    small_lipophilic = ("lomustine", "carmustine")
    large_polar = ("methotrexate", "doxorubicin", "etoposide", "vincristine")
    assert min(scores[a] for a in small_lipophilic) > 10.0 * max(
        scores[b] for b in large_polar)

    # Methotrexate is the least penetrant of the six drugs in the measured
    # human-glioma scorecard (`results/glioma_regional_pk.md`) — the ordering the
    # brief pins. Doxorubicin is NOT in that set; in the library it edges below
    # methotrexate (0.0206 vs 0.0225) because it is bigger with one more HBD,
    # which is the veto behaving correctly, not a regression.
    scorecard = {
        "eribulin": dict(mw=729.9, logp=3.44, tpsa=146.4, hbd=2),
        "cefazolin": dict(mw=454.5, logp=-0.64, tpsa=156.1, hbd=2),
        "methotrexate": dict(mw=454.4, logp=0.27, tpsa=210.5, hbd=5),
        "temozolomide": dict(mw=194.2, logp=-2.08, tpsa=108.2, hbd=1),
        "levetiracetam": dict(mw=170.2, logp=-0.13, tpsa=63.4, hbd=1),
        "adavosertib": dict(mw=500.6, logp=2.89, tpsa=104.3, hbd=2),
    }
    sc = {k: permeability.bbb_penetrance({**v, "is_antibody": False})
          for k, v in scorecard.items()}
    assert min(sc, key=sc.get) == "methotrexate"
    assert max(sc, key=sc.get) == "levetiracetam"
    # BBB-penetrant vs BBB-excluded probe separation (Oh 2026 measured 4.66x
    # pooled, 9.3x in intact-barrier tissue); pre-fix the model gave 3.65x.
    assert sc["levetiracetam"] / sc["cefazolin"] > 3.65


# =========================================================================== #
# FIX 3 -- the permeability ceiling
# =========================================================================== #
def test_lipophilic_molecule_does_not_produce_an_absurd_permeability():
    """`membrane_permeability` returned 4.23e+07 um/s (42 m/s) for imatinib,
    because it exponentiated a Crippen logP of 4.59 with no bound
    (`results/gist_imatinib.md` item 7)."""
    imatinib = dict(mw=493.6, logp=4.59, tpsa=86.3, hbd=2, hba=7, is_antibody=False)
    p = permeability.membrane_permeability(imatinib)
    assert math.isfinite(p) and p > 0.0
    assert p < 100.0, f"absurd permeability {p:.4g} um/s"
    # sanity against the physical scale: metres per second is not a permeability
    assert p < 1e6 / 3.0e8 * 1e6  # trivially true, but states the intent: << c


def test_permeability_saturates_rather_than_climbing_with_logp():
    """The unstirred aqueous boundary layer is in SERIES with the membrane, so
    it is a hard ceiling: past some lipophilicity, more logP buys nothing."""
    base = dict(mw=400.0, tpsa=70.0, hbd=1, hba=5, is_antibody=False)
    ps = [permeability.membrane_permeability({**base, "logp": lp})
          for lp in (2.0, 3.0, 4.0, 6.0, 10.0, 20.0)]
    assert all(b >= a - 1e-12 for a, b in zip(ps, ps[1:])), "must stay monotone"
    assert ps[-1] == pytest.approx(ps[2], rel=1e-6), "must saturate"
    assert ps[-1] < 100.0


def test_whole_library_permeability_is_physically_bounded():
    for name in ("temozolomide", "lomustine", "carmustine", "procarbazine",
                 "methotrexate", "doxorubicin", "cisplatin", "etoposide",
                 "vincristine", "bevacizumab"):
        d = library.get_molecule(name).descriptors
        p = permeability.membrane_permeability(d)
        assert math.isfinite(p), name
        assert 0.0 <= p < 100.0, f"{name}: {p:.4g} um/s"
    # the antibody remains a hard zero
    assert permeability.membrane_permeability(
        library.get_molecule("bevacizumab").descriptors) == 0.0
