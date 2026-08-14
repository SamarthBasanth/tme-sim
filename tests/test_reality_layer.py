"""Regression tests for the calibration + binding-site-barrier layers."""
import numpy as np
import pytest

from tmesim.molecular import library
from tmesim.molecular.binding_site_barrier import (
    apply_barrier_to_drugparams, binding_site_barrier_factor)
from tmesim.validation.calibration import (
    DOX_EFOLD_BAND_UM, DOX_TARGET_EFOLD_UM, calibrate_alpha, calibrated_drugparams,
    penetration_length_um)


def test_alpha_calibration_hits_anchor():
    # The published doxorubicin anchor is a 40-50 um HALF-distance; the sim's L is an
    # E-FOLDING length, so the fit target is d_half/ln2 (band ~58-72 um, mid ~64.9 um).
    # calibrate_alpha defaults to exactly that converted target.
    alpha = calibrate_alpha(drug="doxorubicin")
    assert alpha > 0 and np.isfinite(alpha)
    # Order of magnitude matches the "~2 orders too slow" diagnosis.
    #
    # RE-PINNED at the two-route promotion (results/two_route_promotion.md).
    # SUPERSEDED BAND: 50 < alpha < 2000, when alpha was 148.6082. Doxorubicin's
    # compiled D fell 6.616x (184.7486 -> 27.9256 um^2/s) because the shipping
    # compile now applies the two-route cell term, and alpha is fitted in closed
    # form from L = sqrt(D_eff/k_eff) against the SAME measured 64.92 um anchor,
    # so it had to follow. The band is widened at the bottom rather than moved,
    # because what it is guarding is the ORDER OF MAGNITUDE of the diagnosis --
    # that the structure-only binding prior is ~1-2 orders too small -- and that
    # diagnosis is unchanged: alpha is still an amplification of ~22x, i.e. the
    # prior is still more than an order of magnitude short. The anchor assertion
    # below is the one that actually pins correctness, and it is untouched: L
    # must still land in the measured band, which is the whole point of the fit.
    assert 5 < alpha < 2000
    dp = calibrated_drugparams("doxorubicin", alpha=alpha)
    L = penetration_length_um(dp)
    # sits at its calibration anchor by construction: the e-fold-equivalent band
    assert DOX_EFOLD_BAND_UM[0] <= L <= DOX_EFOLD_BAND_UM[1]


def test_calibrate_alpha_target_is_an_efold_length():
    # Explicit guard on the convention: passing the RAW half-distance (45 um) as if it
    # were an e-folding length inflates alpha by exactly 1/ln(2)^2 = 2.081x. This is
    # the historical 309.6-vs-148.6 inconsistency.
    a_efold = calibrate_alpha(drug="doxorubicin")
    a_raw_half = calibrate_alpha(target_efold_um=45.0, drug="doxorubicin")
    # Each calibration lands on its own target by construction, and the two targets
    # differ by exactly ln2 (so the alphas differ by ~1/ln^2(2), up to the tiny
    # un-amplified passive floor).
    L_efold = penetration_length_um(calibrated_drugparams("doxorubicin", alpha=a_efold))
    L_raw = penetration_length_um(calibrated_drugparams("doxorubicin", alpha=a_raw_half))
    assert L_efold == pytest.approx(DOX_TARGET_EFOLD_UM, rel=1e-6)
    assert L_raw == pytest.approx(45.0, rel=1e-6)
    assert L_efold / L_raw == pytest.approx(1.0 / np.log(2.0), rel=1e-6)
    assert a_raw_half / a_efold == pytest.approx(1.0 / np.log(2.0) ** 2, rel=1e-2)
    assert DOX_TARGET_EFOLD_UM == pytest.approx(45.0 / np.log(2.0), rel=1e-9)


def test_calibration_generalizes_direction_to_temozolomide():
    # out-of-sample: weakly-bound TMZ must penetrate much further than doxorubicin
    alpha = calibrate_alpha(drug="doxorubicin")
    L_dox = penetration_length_um(calibrated_drugparams("doxorubicin", alpha=alpha))
    L_tmz = penetration_length_um(calibrated_drugparams("temozolomide", alpha=alpha))
    assert L_tmz > 3 * L_dox


def test_binding_site_barrier_traps_antibody_only():
    alpha = calibrate_alpha(drug="doxorubicin")
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
