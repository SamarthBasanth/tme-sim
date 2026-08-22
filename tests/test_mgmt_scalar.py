"""MGMT as a per-patient temozolomide-sensitivity scalar.

The load-bearing test in this file is
:func:`test_omitting_mgmt_is_bit_identical_to_previous_behaviour`: the whole
feature is opt-in, and a caller who does not supply an MGMT value must get
*exactly* the numbers they got before it existed — same sensitivity array, same
info block, same simulated kill fraction, bit for bit.

The rest asserts the mapping's shape (direction, bounds, anchor), that the
frozen calibration constants still match a live refit from the cached primary
bytes, and that a supplied MGMT value is always RECORDED rather than silently
absorbed.
"""
from __future__ import annotations

import numpy as np
import pytest

from tmesim.molecular.mgmt import (
    MGMT_BINARY_SENSITIVITY,
    MGMT_CALIBRATED_DRUGS,
    MGMT_LAMBDA,
    MGMT_SENSITIVITY_FLOOR,
    is_mgmt_calibrated,
    mgmt_sensitivity,
)
from tmesim.reconstruct.real_pipeline import (
    apply_drug_sensitivity,
    reconstruction_report,
)


# --------------------------------------------------------------------------- #
# THE backward-compatibility test
# --------------------------------------------------------------------------- #
def test_omitting_mgmt_is_bit_identical_to_previous_behaviour(real_micro):
    """No MGMT supplied -> byte-for-byte the pre-MGMT behaviour.

    ``mgmt=None`` is the default, and the MGMT branch is only entered when the
    caller passes something. This asserts the observable consequences: the
    sensitivity array, every key of the returned ``info``, and the report block
    are identical whether or not the (defaulted) argument is spelled out.
    """
    work = real_micro.copy() if hasattr(real_micro, "copy") else real_micro

    for mode in ("uniform", "determinant"):
        for drug in ("temozolomide", "doxorubicin", "cisplatin"):
            implicit = apply_drug_sensitivity(work, drug, mode=mode)
            arr_implicit = np.asarray(work.cells.sensitivity,
                                      dtype=np.float64).copy()
            rep_implicit = reconstruction_report(work)["sensitivity"]

            explicit = apply_drug_sensitivity(work, drug, mode=mode, mgmt=None)
            arr_explicit = np.asarray(work.cells.sensitivity,
                                      dtype=np.float64).copy()
            rep_explicit = reconstruction_report(work)["sensitivity"]

            # bit-identical sensitivity
            assert arr_implicit.tobytes() == arr_explicit.tobytes(), (
                f"{drug}/{mode}: passing mgmt=None changed the sensitivity array"
            )
            assert implicit == explicit
            assert rep_implicit == rep_explicit

            # and the *old* contract still holds exactly
            if mode == "uniform":
                assert implicit["mean"] == 1.0
                assert implicit["is_uniform"] is True
                np.testing.assert_array_equal(arr_implicit, 1.0)
            assert implicit["upper_bound"] == (
                implicit["determinant_gene"] is None
            )
            # no MGMT key ever asserts an applied scalar
            assert implicit.get("mgmt_applied", False) is False
            assert implicit.get("mgmt_scalar") is None
            assert rep_implicit["mgmt_applied"] is False
            assert rep_implicit["mgmt_scalar"] is None


def test_omitting_mgmt_reproduces_the_same_simulated_kill(real_micro, real_drug):
    """End-to-end: the engine's kill fraction is unchanged when mgmt is omitted."""
    from tmesim.core import PDEConfig, SimConfig
    from tmesim.sim.engine import simulate_response

    work = real_micro.copy() if hasattr(real_micro, "copy") else real_micro
    pde = PDEConfig(dt=10.0, total_time=600.0, convection=False, save_every=100)

    apply_drug_sensitivity(work, "temozolomide", mode="uniform")
    a = simulate_response(work, SimConfig(drug=real_drug, pde=pde, seed=0))

    apply_drug_sensitivity(work, "temozolomide", mode="uniform", mgmt=None)
    b = simulate_response(work, SimConfig(drug=real_drug, pde=pde, seed=0))

    assert a.kill_fraction == b.kill_fraction

    # ... and a supplied MGMT value DOES move it, in the resistant direction.
    apply_drug_sensitivity(work, "temozolomide", mode="uniform", mgmt=4.0)
    c = simulate_response(work, SimConfig(drug=real_drug, pde=pde, seed=0))
    assert c.kill_fraction <= a.kill_fraction


# --------------------------------------------------------------------------- #
# The mapping
# --------------------------------------------------------------------------- #
def test_mgmt_silent_is_the_anchor():
    """s(MGMT = 0) == 1.0 exactly: the uniform default is the null special case."""
    assert mgmt_sensitivity(0.0) == 1.0
    assert mgmt_sensitivity(0.0, units="cpm") == 1.0
    assert mgmt_sensitivity(0.0, model="linear") == pytest.approx(1.0)


def test_direction_high_mgmt_is_resistant():
    """High MGMT -> repaired O6-meG -> resistant -> LOW sensitivity."""
    xs = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    s = [mgmt_sensitivity(x) for x in xs]
    assert all(s[i] > s[i + 1] for i in range(len(s) - 1)), s
    assert mgmt_sensitivity(status="methylated") > mgmt_sensitivity(
        status="unmethylated"
    )


def test_bounded_in_unit_interval_including_far_extrapolation():
    for x in (0.0, 0.5, 4.47, 10.0, 50.0):
        for model in ("exponential", "linear"):
            v = mgmt_sensitivity(x, model=model)
            assert MGMT_SENSITIVITY_FLOOR <= v <= 1.0, (x, model, v)


def test_units_conversion_is_consistent():
    import math

    assert mgmt_sensitivity(15.33, units="cpm") == pytest.approx(
        mgmt_sensitivity(math.log2(15.33 + 1.0))
    )


def test_binary_status_aliases_and_values():
    assert mgmt_sensitivity(status="methylated") == 1.0
    assert mgmt_sensitivity(status="M") == 1.0
    assert mgmt_sensitivity(status="unmethylated") == pytest.approx(
        MGMT_BINARY_SENSITIVITY["unmethylated"]
    )
    assert mgmt_sensitivity(status="UM") == mgmt_sensitivity(status="unmethylated")
    with pytest.raises(ValueError):
        mgmt_sensitivity(status="probably")


def test_exactly_one_input_required():
    with pytest.raises(ValueError):
        mgmt_sensitivity()
    with pytest.raises(ValueError):
        mgmt_sensitivity(1.0, status="methylated")
    with pytest.raises(ValueError):
        mgmt_sensitivity(-1.0, units="cpm")
    with pytest.raises(ValueError):
        mgmt_sensitivity(1.0, units="fpkm")
    with pytest.raises(ValueError):
        mgmt_sensitivity(1.0, model="quadratic")


def test_detail_block_carries_the_uncertainty_and_extrapolation_flag():
    d = mgmt_sensitivity(4.0, detail=True)
    assert d["n"] == 17
    lo, hi = d["ci95"]
    assert lo <= d["sensitivity"] <= hi
    assert d["extrapolated"] is False
    far = mgmt_sensitivity(4.73, detail=True)
    assert far["extrapolated"] is True
    assert "ABOVE the calibration ceiling" in far["extrapolation_note"]


def test_lambda_matches_a_live_refit_from_the_primary_bytes():
    """The frozen constants are auditable, not asserted."""
    from tmesim.validation.mgmt_scalar import verify_frozen_constants

    chk = verify_frozen_constants()
    assert chk["ok"], chk["mismatches"]
    assert chk["fit"]["n"] == 17
    assert chk["fit"]["lambda"] == pytest.approx(MGMT_LAMBDA)
    # the direction and the significance the whole thing rests on
    assert chk["fit"]["spearman_x_vs_viability"][0] > 0
    assert chk["fit"]["spearman_x_vs_viability"][1] < 0.05
    # one covariate must still beat the constant model out-of-sample
    assert chk["fit"]["loo_mae_exponential"] < chk["fit"]["loo_mae_constant"]


# --------------------------------------------------------------------------- #
# Visium recovery
# --------------------------------------------------------------------------- #
def test_visium_pseudobulk_mgmt_is_recovered_and_labelled_pseudobulk():
    from tmesim.data.real.visium import list_samples
    from tmesim.molecular.mgmt import mgmt_from_visium_all

    # A MISSING DATASET IS NOT A FINDING. data_cache/ is real tissue data, is
    # gitignored, and is absent on every CI runner -- so without this guard the
    # test below fails there for a reason that has nothing to do with MGMT, and
    # pollutes a job (.github/workflows/docker.yml, the linux/amd64 reference
    # job) whose only purpose is float reproducibility.
    #
    # Ask the loader, not a hardcoded path: `list_samples()` defaults to
    # present_only=True, i.e. it already answers "whose bytes are locally
    # cached". So this skips exactly when there is nothing to read and never
    # when there is -- the assertion below is unchanged and still runs on every
    # machine that has the cache.
    if not list_samples():
        pytest.skip("no Visium GBM section is cached in this checkout "
                    "(see data_cache/ and `vivum sections`)")

    table = mgmt_from_visium_all()
    assert table, "no cached Visium GBM sections found"
    for sid, rec in table.items():
        if not rec.get("mgmt_in_reference"):
            continue
        assert 0.0 < rec["pct_spots_mgmt_positive"] < 100.0
        assert rec["pseudobulk_cpm"] > 0.0
        assert MGMT_SENSITIVITY_FLOOR <= rec["sensitivity"] <= 1.0
        # honesty: never sold as a spatial field
        assert "PSEUDOBULK" in rec["resolution"]
        assert "not a spatial field" in rec["caveat"].lower()


# --------------------------------------------------------------------------- #
# Wiring
# --------------------------------------------------------------------------- #
def test_supplying_mgmt_records_it_and_drops_the_upper_bound_flag(real_micro):
    work = real_micro.copy() if hasattr(real_micro, "copy") else real_micro

    plain = apply_drug_sensitivity(work, "temozolomide", mode="determinant")
    assert plain["upper_bound"] is True
    assert "UPPER BOUND" in plain["comparability"]

    informed = apply_drug_sensitivity(work, "temozolomide", mode="determinant",
                                      mgmt=4.029)
    assert informed["mgmt_applied"] is True
    assert informed["determinant_gene"] == "MGMT"
    assert informed["upper_bound"] is False
    assert "UPPER BOUND" not in informed["comparability"]
    assert informed["mgmt_scalar"] == pytest.approx(mgmt_sensitivity(4.029))
    assert informed["mgmt_calibration_n"] == 17
    assert informed["mgmt_scalar_ci95"] is not None
    np.testing.assert_allclose(work.cells.sensitivity, informed["mgmt_scalar"])

    # and the report surfaces it, so a result can never carry it silently
    block = reconstruction_report(work)["sensitivity"]
    assert block["mgmt_applied"] is True
    assert block["mgmt_scalar"] == pytest.approx(informed["mgmt_scalar"])
    assert block["upper_bound"] is False
    assert "SCALAR" in block["mgmt_resolution"]


def test_mgmt_is_a_recorded_noop_for_uncalibrated_drugs(real_micro):
    """A panel loop can pass mgmt for every drug; only TMZ uses it."""
    work = real_micro.copy() if hasattr(real_micro, "copy") else real_micro

    for drug in ("doxorubicin", "lomustine", "procarbazine"):
        assert not is_mgmt_calibrated(drug)
        info = apply_drug_sensitivity(work, drug, mode="uniform", mgmt=4.0)
        assert info["mgmt_supplied"] is True
        assert info["mgmt_applied"] is False
        assert "not in MGMT_CALIBRATED_DRUGS" in info["mgmt_note"]
        np.testing.assert_allclose(work.cells.sensitivity, 1.0)

    assert MGMT_CALIBRATED_DRUGS == ("temozolomide",)
    assert is_mgmt_calibrated("Temozolomide")


def test_binary_status_can_be_passed_through_the_pipeline(real_micro):
    work = real_micro.copy() if hasattr(real_micro, "copy") else real_micro
    m = apply_drug_sensitivity(work, "temozolomide", mode="uniform",
                               mgmt="unmethylated")
    assert m["mgmt_applied"] is True
    assert m["mgmt_scalar"] == pytest.approx(
        MGMT_BINARY_SENSITIVITY["unmethylated"]
    )
    u = apply_drug_sensitivity(work, "temozolomide", mode="uniform",
                               mgmt="methylated")
    assert u["mgmt_scalar"] == 1.0
    np.testing.assert_allclose(work.cells.sensitivity, 1.0)
