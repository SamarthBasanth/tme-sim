"""The half-distance estimator's bin-width bias, and the additions that measure it.

Companion to ``results/spatial_noise_decomposition.md`` and
``experiments/spatial_noise_decomposition.py``.

``results/kyle2004_penetration.md`` §8 found that on the reconstructed section the
compound difference sits inside its own methodological noise, while on the clean
slab it does not — so the instability is in the SPATIAL REDUCTION, not the
transport physics. The decomposition attributes the dominant part of that spread
to a single stage: the DISTANCE BIN WIDTH, acting through
:func:`empirical_half_distance_um`.

The mechanism is closed-form, and these tests pin it as a property of the code
rather than as a claim in a report:

* the shipped estimator overstates a pure exponential's half-distance by exactly
  half a bin width, because it normalises by the FIRST BIN'S MEAN (a value that
  sits at ``bin_um/2``) and then reports the crossing as a distance from ZERO;
* :func:`anchor_corrected_half_distance_um` removes that term exactly, by
  measuring the crossing from the anchor's own location instead of from zero;
* the additive ``distance_mean_um`` / ``n_empty_bins`` keys report where each
  bin's value actually sits and how much of the profile the lattice left empty.

Nothing here asserts that the corrected estimator produces a BETTER number
against any measurement. It asserts that the shipped one has a derivable,
bin-width-proportional offset and that the correction removes it — which is a
statement about arithmetic, not about the tissue.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from tmesim.validation.vessel_profile import (
    anchor_corrected_half_distance_um,
    distance_profile,
    empirical_half_distance_um,
)

LN2 = math.log(2.0)

#: Bin widths swept by the registered arm-B nuisance grid
#: (``results/kyle2004_penetration.md`` §8).
NUISANCE_BINS_UM = (5.0, 10.0, 20.0)


def _exponential_cloud(L: float = 50.0, floor: float = 0.0,
                       pitch: float = 0.05, extent: float = 400.0):
    """A dense 1-D pixel cloud, uniform in distance, carrying
    ``C(d) = (1-floor)*exp(-d/L) + floor``.

    Uniform-in-distance is the geometry under which the closed-form bias below is
    exact, so a failure here is an arithmetic failure and not a geometry
    artefact. ``pitch`` is far finer than any swept bin width, so the lattice
    itself contributes nothing.
    """
    d = np.arange(0.0, extent, pitch) + 0.5 * pitch
    c = (1.0 - floor) * np.exp(-d / L) + floor
    return d, c


def _profile(d, c, bin_um, extent=400.0):
    """Bin the cloud through the shipped ``distance_profile``, with the whole
    cloud inside the tissue mask so nothing is filtered."""
    n = d.size
    return distance_profile(
        field=c.reshape(1, n),
        distance_um=d.reshape(1, n),
        tissue_mask=np.ones((1, n), dtype=bool),
        bin_um=bin_um,
        max_distance_um=extent,
    )


# --------------------------------------------------------------------------- #
# 1. distance_profile's pre-existing contract is untouched
# --------------------------------------------------------------------------- #
def test_distance_profile_keeps_its_original_keys_and_values():
    """The additive diagnostics must not have disturbed anything a caller already
    reads. ``delivery_metrics``, ``vessel_network_impact`` and ``vessel_profile``
    itself all index this dict by key."""
    d, c = _exponential_cloud()
    prof = _profile(d, c, 10.0)
    for key in ("distance_um", "mean", "median", "std", "q25", "q75", "count",
                "bin_um"):
        assert key in prof, f"{key} disappeared from distance_profile"
    # Bin centres are still the geometric centres, unchanged.
    assert prof["distance_um"][0] == pytest.approx(5.0)
    assert prof["distance_um"][1] == pytest.approx(15.0)
    # The first bin's mean is still the plain mean of the values in [0, 10).
    expected = float(c[(d >= 0.0) & (d < 10.0)].mean())
    assert prof["mean"][0] == pytest.approx(expected, rel=1e-12)


def test_distance_mean_um_is_where_the_bin_actually_sits():
    """The new key is the bin's COUNT-WEIGHTED MEAN distance.

    This test is written in two halves ON PURPOSE. The first half — a
    uniform-in-distance cloud, where the mean and the centre coincide — is what
    makes the ``b/2`` offset of the tests below a property of the ANCHOR rather
    than of the labelling. But that half alone would pass against an
    implementation that simply returned the bin centre, so it does not test what
    the key claims. The second half uses a cloud whose pixel density RISES with
    distance — which is the real 2-D situation, and the flagship section's actual
    behaviour (bin centre 10.0 µm against an observed mean of 12.69 µm at
    ``b = 20``) — and there the two must differ.
    """
    d, c = _exponential_cloud()
    prof = _profile(d, c, 10.0)
    assert "distance_mean_um" in prof
    assert prof["distance_mean_um"].shape == prof["distance_um"].shape
    np.testing.assert_allclose(prof["distance_mean_um"], prof["distance_um"],
                               atol=0.05)

    # Density rising linearly with distance: sample d with weight proportional
    # to d, so each bin's mean sits above its centre by a computable amount.
    rng = np.random.default_rng(0)
    u = rng.random(400_000)
    d2 = 100.0 * np.sqrt(u)                       # p(d) ∝ d on [0, 100]
    c2 = np.exp(-d2 / 50.0)
    prof2 = _profile(d2, c2, 20.0, extent=100.0)
    centres = prof2["distance_um"]
    means = prof2["distance_mean_um"]
    assert np.all(means > centres), (centres, means)
    # For p(d) ∝ d the mean of [lo, hi] is (2/3)(hi^3-lo^3)/(hi^2-lo^2); for the
    # first bin [0, 20) that is 13.333 µm against a centre of 10.
    assert means[0] == pytest.approx(40.0 / 3.0, abs=0.1), means[0]
    assert centres[0] == pytest.approx(10.0)


def test_n_empty_bins_reports_lattice_aliasing_that_was_previously_silent():
    """``distance_profile`` skips empty bins with ``if n == 0: continue``, so an
    aliased profile used to be indistinguishable from a well-sampled one.

    A distance map sampled on a lattice of pitch ``dx`` puts pixels on a comb of
    distances. Binning that comb finer than the pitch leaves gaps. This is the
    ``dx30_native`` arm's situation at the shipped ``BIN_UM = 10``.
    """
    # A coarse comb: distances only ever land on multiples of 30 um.
    d = np.repeat(np.arange(0.0, 300.0, 30.0), 50)
    c = np.exp(-d / 50.0)
    prof = _profile(d, c, bin_um=10.0, extent=300.0)
    assert prof["n_empty_bins"] > 0, (
        "a 10 um binning of a 30 um-pitch distance comb must leave empty bins")
    # Binned at or above the pitch, every bin is populated.
    prof30 = _profile(d, c, bin_um=30.0, extent=300.0)
    assert prof30["n_empty_bins"] == 0


# --------------------------------------------------------------------------- #
# 2. The shipped estimator's bias is exactly half a bin width
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("bin_um", NUISANCE_BINS_UM)
def test_shipped_half_distance_overstates_by_exactly_half_a_bin(bin_um):
    """``d_half_shipped = L*ln2 + bin_um/2`` for a pure exponential.

    Derivation (also in ``anchor_corrected_half_distance_um``'s docstring): the
    bin means of ``C0*exp(-d/L)`` are ``C0*(L/b)*(1-exp(-b/L))*exp(-j*b/L)``, so
    ``v_j/v_0 = exp(-j*b/L)`` — the bin-averaging prefactor cancels — the profile
    halves at fractional bin index ``L*ln2/b``, and the estimator reports that at
    the bin CENTRE, i.e. at ``L*ln2 + b/2``.
    """
    L = 50.0
    d, c = _exponential_cloud(L=L)
    prof = _profile(d, c, bin_um)
    got = empirical_half_distance_um(prof["distance_um"], prof["mean"])
    assert got is not None
    predicted = L * LN2 + 0.5 * bin_um
    # The only residual is SECOND order: the estimator interpolates LINEARLY
    # between bin centres across a convex curve, and the chord lies above it. For
    # a curve of radius-of-curvature scale L sampled at pitch b that error is
    # bounded by b^2/(8L) — 0.06 um at b = 5, 1.0 um at b = 20. The tolerance is
    # that bound, derived, not a number chosen to make the test pass.
    tol = bin_um ** 2 / (8.0 * L)
    assert got == pytest.approx(predicted, abs=tol), (
        f"bin {bin_um}: got {got:.4f}, closed form {predicted:.4f}, "
        f"2nd-order bound {tol:.4f}")


def test_the_bias_is_what_makes_the_readout_bin_width_dependent():
    """The whole point: the shipped readout moves across the registered nuisance
    bin widths by ~half the widest bin, on a profile whose true half-distance is
    a single fixed number."""
    L = 50.0
    d, c = _exponential_cloud(L=L)
    vals = [empirical_half_distance_um(*(lambda p: (p["distance_um"], p["mean"]))(
        _profile(d, c, b))) for b in NUISANCE_BINS_UM]
    span = max(vals) - min(vals)
    # First order the span is 20/2 - 5/2 = 7.5 um of pure estimator offset, on a
    # 34.66 um quantity; the second-order chord term (bounded by b^2/(8L)) can
    # only widen it.
    lo, hi = 7.5, 7.5 + 20.0 ** 2 / (8.0 * L)
    assert lo <= span <= hi, (vals, span)
    assert span / (L * LN2) > 0.2, (
        "the bin-width offset should be >20% of the quantity being reported")


# --------------------------------------------------------------------------- #
# 3. The anchor correction removes it, exactly
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("bin_um", NUISANCE_BINS_UM)
def test_anchor_correction_recovers_the_true_half_distance(bin_um):
    L = 50.0
    d, c = _exponential_cloud(L=L)
    prof = _profile(d, c, bin_um)
    got = anchor_corrected_half_distance_um(prof["distance_um"], prof["mean"])
    # The b/2 term is gone exactly; what is left is the same second-order chord
    # bound as above. Nothing between the two estimators is fitted.
    tol = bin_um ** 2 / (8.0 * L)
    assert got == pytest.approx(L * LN2, abs=tol), (
        f"bin {bin_um}: corrected {got:.4f} vs true {L * LN2:.4f}, "
        f"2nd-order bound {tol:.4f}")


def test_anchor_correction_is_an_order_of_magnitude_less_bin_sensitive():
    """Same cloud, same bins, both estimators — compare their spans directly.

    This is the mechanical claim, not a score: the correction removes a term
    proportional to the bin width, so what is left is the second-order
    interpolation residual.
    """
    d, c = _exponential_cloud(L=50.0)
    shipped, corrected = [], []
    for b in NUISANCE_BINS_UM:
        p = _profile(d, c, b)
        shipped.append(empirical_half_distance_um(p["distance_um"], p["mean"]))
        corrected.append(
            anchor_corrected_half_distance_um(p["distance_um"], p["mean"]))
    span_shipped = max(shipped) - min(shipped)
    span_corrected = max(corrected) - min(corrected)
    assert span_corrected < span_shipped / 8.0, (
        f"shipped span {span_shipped:.3f} um, corrected {span_corrected:.3f} um")


def test_anchor_correction_uses_the_supplied_anchor_when_given():
    d, c = _exponential_cloud(L=50.0)
    prof = _profile(d, c, 10.0)
    default = anchor_corrected_half_distance_um(prof["distance_um"], prof["mean"])
    explicit = anchor_corrected_half_distance_um(
        prof["distance_um"], prof["mean"],
        anchor_distance_um=float(prof["distance_mean_um"][0]))
    # For a uniform-in-distance cloud the two anchors coincide.
    assert explicit == pytest.approx(default, abs=0.05)
    shifted = anchor_corrected_half_distance_um(
        prof["distance_um"], prof["mean"], anchor_distance_um=0.0)
    assert shifted == pytest.approx(
        empirical_half_distance_um(prof["distance_um"], prof["mean"]))


# --------------------------------------------------------------------------- #
# 4. The correction inherits the shipped estimator's honest refusal
# --------------------------------------------------------------------------- #
def test_correction_returns_none_when_the_profile_never_halves():
    """``empirical_half_distance_um`` returns ``None`` BY DESIGN for a profile
    that never falls to half its perivascular value (``D_HALF_NEVER``); the
    correction must propagate that rather than manufacture a number."""
    d = np.arange(0.0, 200.0, 1.0) + 0.5
    c = 1.0 - 0.001 * d          # falls, but never to half
    prof = _profile(d, c, 10.0, extent=200.0)
    assert empirical_half_distance_um(prof["distance_um"], prof["mean"]) is None
    assert anchor_corrected_half_distance_um(prof["distance_um"],
                                             prof["mean"]) is None


def test_correction_is_defined_on_a_profile_with_a_floor_and_stays_finite():
    """The real profile is not a pure exponential — the shipped fit reports a
    floor fraction of ~0.12 on the flagship section. The floor makes the
    normalisation prefactor stop cancelling, so a residual bin-width dependence
    survives the correction. That residual is a DESIGN question (what "half"
    means under a background), not an arithmetic error, and this test records
    that it is bounded rather than pretending it is absent.
    """
    d, c = _exponential_cloud(L=50.0, floor=0.12)
    shipped, corrected = [], []
    for b in NUISANCE_BINS_UM:
        p = _profile(d, c, b)
        shipped.append(empirical_half_distance_um(p["distance_um"], p["mean"]))
        corrected.append(
            anchor_corrected_half_distance_um(p["distance_um"], p["mean"]))
    assert all(v is not None for v in corrected)
    assert (max(corrected) - min(corrected)) < (max(shipped) - min(shipped))
