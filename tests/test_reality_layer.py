"""Regression tests for the calibration + binding-site-barrier layers."""
import numpy as np
import pytest

from tmesim.molecular import library
from tmesim.molecular.binding_site_barrier import (
    apply_barrier_to_drugparams, binding_site_barrier_factor)
from tmesim.validation.calibration import (
    calibrate_alpha, calibrated_drugparams, penetration_length_um)


def test_alpha_calibration_hits_anchor():
    alpha = calibrate_alpha(target_um=45.0, drug="doxorubicin")
    assert alpha > 0 and np.isfinite(alpha)
    # order of magnitude matches the "~300x too slow" diagnosis
    assert 50 < alpha < 2000
    dp = calibrated_drugparams("doxorubicin", alpha=alpha)
    L = penetration_length_um(dp)
    assert 40 <= L <= 50  # doxorubicin sits at its calibration anchor by construction


def test_calibration_generalizes_direction_to_temozolomide():
    # out-of-sample: weakly-bound TMZ must penetrate much further than doxorubicin
    alpha = calibrate_alpha(target_um=45.0, drug="doxorubicin")
    L_dox = penetration_length_um(calibrated_drugparams("doxorubicin", alpha=alpha))
    L_tmz = penetration_length_um(calibrated_drugparams("temozolomide", alpha=alpha))
    assert L_tmz > 3 * L_dox


def test_binding_site_barrier_traps_antibody_only():
    alpha = calibrate_alpha(target_um=45.0, drug="doxorubicin")
    bev = library.get_molecule("bevacizumab")
    tmz = library.get_molecule("temozolomide")
    dp_bev = calibrated_drugparams("bevacizumab", alpha=alpha)
    dp_tmz = calibrated_drugparams("temozolomide", alpha=alpha)
    # antibody penetration collapses under the barrier; small molecule ~unchanged
    L_bev_before = penetration_length_um(dp_bev)
    L_bev_after = penetration_length_um(apply_barrier_to_drugparams(dp_bev, bev))
    L_tmz_after = penetration_length_um(apply_barrier_to_drugparams(dp_tmz, tmz))
    assert L_bev_after < L_bev_before          # barrier shortens antibody reach
    assert L_bev_after < L_tmz_after / 10      # antibody << small molecule
    # small molecule essentially unaffected
    f_tmz = binding_site_barrier_factor(tmz)
    assert not f_tmz["is_barrier_limited"]
