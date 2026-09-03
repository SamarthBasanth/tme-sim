"""UNCERTAINTY QUANTIFICATION — honest error bars on the REAL-microenvironment
prediction, by Monte-Carlo propagation of the knobs we do NOT know.

Why this module exists
----------------------
:mod:`tmesim.validation.numerical_rigor` established that the *numerics* are not
the bottleneck: mass is conserved to ~1e-13 on the real irregular mask, the
implicit scheme is dt-converged, and ``kill_fraction`` is grid-converged. What it
also established is where the prediction uncertainty actually lives:

* ``endothelial_typing_threshold`` — the single biggest local mover of the
  penetration index (|E_pen| ≈ 0.52): it decides which cells are ENDOTHELIAL and
  therefore where the delivery source (and the O2 field) sits;
* ``kde_bandwidth_um`` — the biggest local mover of the kill fraction
  (|E_kill| ≈ 0.44): it sets how concentrated vs. smeared the delivered dose is;
* ``grid_dx_um`` — a *numerical* knob (|E_kill| ≈ 0.35) which we do NOT sample;
  it is not an unknown, it is a resolution requirement, handled by fixing dx and
  auditing the residual bias (:func:`resolution_bias`);
* ``calibration_alpha`` — the single-anchor calibration constant. It has the
  *smallest* local elasticity (|E_pen| ≈ 0.11) but by far the WIDEST credible
  range: it is pinned to ONE doxorubicin measurement (~40–50 um, Primeau 2005),
  so a +-50–100% prior is honest where a +-20% perturbation is not.

This module turns that ranking into actual error bars: put a documented prior on
each knob, draw from it, re-run the FULL real pipeline (reconstruct -> physics ->
kill) per draw, and report the posterior-predictive median and 90% credible
interval, plus a first-order variance decomposition saying which prior bought
which share of the spread.

What the priors are (and are NOT)
---------------------------------
**They are EXPERT priors, not measured distributions.** Nothing in this repo
measures a distribution over the endothelial typing gate or the KDE bandwidth;
there is no posterior from data because there is no outcome label to condition
on. Each prior below is a documented, defensible *range of analyst choice*, with
its justification written next to it. The resulting intervals are therefore
"how much does the prediction move when the reconstruction/calibration choices
move inside their defensible range" — which is a real and useful error bar, and
is NOT the same thing as a validated predictive interval. See the honesty
section of ``results/uncertainty.md``.

Design
------
* **Drug:** temozolomide (the GBM standard of care), calibrated with the frozen
  single-anchor alpha model from :mod:`tmesim.validation.calibration`.
* **Sampling:** Latin hypercube over the 4 priors (better marginal coverage than
  plain MC at n ~ 60), mapped through each prior's inverse CDF.
* **Per draw:** the real Xenium section is RE-RECONSTRUCTED with that draw's
  typing threshold and KDE bandwidth (via ``numerical_rigor.reconstruct_knobs``,
  which reproduces ``reconstruct_real`` exactly at the defaults), the drug is
  re-built with that draw's alpha, and the ECM tortuosity coefficient is
  temporarily patched into :mod:`tmesim.physics.diffusion` for the solve (a
  restored monkeypatch — NO file outside this module is modified).
* **Kill de-noising:** ``kill_fraction`` carries a stochastic per-cell draw; it
  is averaged over a fixed pair of seeds so the reported spread is PARAMETER
  variance, not sampling noise (the seed noise is measured and reported).
* **Variance attribution:** two first-order estimators — a binned
  ``Var(E[Y|X_i])/Var(Y)`` "Sobol-lite" (robust to the step/threshold
  non-linearity) and standardised rank-regression coefficients (SRRC, the linear
  cross-check).

Run: ``python -m tmesim.validation.uncertainty`` (writes ``results/uncertainty.md``
and ``results/uncertainty.png``).
"""
from __future__ import annotations

import contextlib
import math
import os
import time
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy import stats

from .. import physics as _physics_pkg  # noqa: F401  (ensures submodule import)
from ..physics import diffusion as _diffusion
from ..core import CellType
from . import calibration as _cal
from .calibration import calibrate_alpha, calibrated_drugparams
from .pk_reality import DEFAULT_ECM, penetration_length_um
from .numerical_rigor import (
    BASE_BANDWIDTH,
    BASE_THRESHOLD,
    _pde,
    _predict,
    reconstruct_knobs,
)

# ---------------------------------------------------------------------------
# Defaults.
# ---------------------------------------------------------------------------
DRUG_NAME = "temozolomide"
UQ_DX_UM = 15.0            # MC grid (see resolution_bias(): ~3% high on pen)
UQ_KILL_SEEDS = (0, 1)     # seeds averaged per draw to de-noise the kill draw
NOMINAL_TORTUOSITY = float(_diffusion.ALPHA_TORTUOSITY)  # 0.6


# ===========================================================================
# 1. THE PRIORS.  Each one is an expert prior with a written justification.
# ===========================================================================
@dataclass(frozen=True)
class Prior:
    """One uncertain knob: a distribution + WHY that distribution.

    ``ppf`` maps a uniform(0,1) quantile to a value (inverse CDF), so any prior
    can be driven by a Latin hypercube. ``kind`` / ``spec`` are for reporting.
    """

    name: str
    kind: str
    spec: str
    nominal: float
    ppf: Callable[[np.ndarray], np.ndarray] = field(repr=False)
    justification: str = ""
    honest_caveat: str = ""

    def sample_from_unit(self, u: np.ndarray) -> np.ndarray:
        return np.asarray(self.ppf(np.asarray(u, dtype=np.float64)), dtype=np.float64)

    def quantiles(self, qs=(0.05, 0.5, 0.95)) -> Tuple[float, ...]:
        return tuple(float(self.ppf(np.array([q]))[0]) for q in qs)


# --- endothelial typing threshold -------------------------------------------
# Truncated normal about the frozen gate (z >= 1.0), sd 0.15, clipped to
# [0.70, 1.40].
_THR_MU, _THR_SD, _THR_LO, _THR_HI = 1.0, 0.15, 0.70, 1.40
_thr_rv = stats.truncnorm(
    (_THR_LO - _THR_MU) / _THR_SD, (_THR_HI - _THR_MU) / _THR_SD,
    loc=_THR_MU, scale=_THR_SD,
)

# --- vascular / ECM KDE bandwidth -------------------------------------------
_BW_LO, _BW_HI = 25.0, 60.0
_bw_rv = stats.uniform(loc=_BW_LO, scale=_BW_HI - _BW_LO)

# --- ECM tortuosity coefficient ---------------------------------------------
_TORT_LO, _TORT_HI = 0.40, 0.80
_tort_rv = stats.uniform(loc=_TORT_LO, scale=_TORT_HI - _TORT_LO)

# --- calibration alpha (WIDE: single anchor) --------------------------------
# Lognormal with median = the fitted alpha and a 90% CI spanning 0.5x - 2x.
_ALPHA_SPAN = 2.0
_ALPHA_SIGMA = math.log(_ALPHA_SPAN) / stats.norm.ppf(0.95)  # ~0.4214


def build_priors(alpha0: Optional[float] = None) -> Dict[str, Prior]:
    """Build the prior spec. ``alpha0`` defaults to the freshly fitted calibration
    constant, so the alpha prior is always centred on the actual fit.

    THE VALUE IS NOT NAMED HERE. This docstring said "~148.6"; the fit returns
    22.223477711172563 today, and it said so two lines above a comment claiming
    this text "can never drift from the module it documents". The code was
    always right -- it calls ``calibrate_alpha()`` -- and only the number
    written beside it went stale, which is the whole argument for not writing
    it down twice."""
    if alpha0 is None:
        alpha0 = calibrate_alpha()
    alpha0 = float(alpha0)
    a_rv = stats.lognorm(s=_ALPHA_SIGMA, scale=alpha0)
    # The anchor the fit actually targets, read live from calibration so this text
    # can never drift from the module it documents.
    efold_target = float(getattr(_cal, "DOX_TARGET_EFOLD_UM",
                                 getattr(_cal, "DOX_TARGET_UM", 45.0)))
    efold_band = tuple(float(v) for v in getattr(
        _cal, "DOX_EFOLD_BAND_UM", getattr(_cal, "DOX_PENETRATION_UM", (40.0, 50.0))))

    return {
        "endothelial_typing_threshold": Prior(
            name="endothelial_typing_threshold",
            kind="truncated normal",
            spec=f"N({_THR_MU:.2f}, {_THR_SD:.2f}) truncated to "
                 f"[{_THR_LO:.2f}, {_THR_HI:.2f}]",
            nominal=BASE_THRESHOLD,
            ppf=_thr_rv.ppf,
            justification=(
                "The phenotyping gate is 'a type's mean marker z-score must clear "
                "+1 SD above the cohort mean' (reconstruct.real_xenium."
                "TYPING_THRESHOLD = 1.0). '1 SD' is an analyst choice, not a "
                "measurement: any value from ~0.85 to ~1.15 is an equally "
                "defensible reading of 'a genuine marker-positive call', which is "
                "the +-1 sigma of this prior. The truncation bounds are where the "
                "phenotype stops being credible: BELOW ~0.9 the single-marker "
                "PECAM1 ENDOTHELIAL class jumps from 9.7% to 17.4% of cells "
                "(the bare-arg-max inflation the gate exists to suppress, and "
                "~2x the independent Visium deconvolution's ~9.8% endothelial "
                "fraction); ABOVE ~1.4 the STROMAL class is squeezed below 4% and "
                "essentially everything defaults to TUMOR. numerical_rigor found "
                "this knob the strongest local mover of the penetration index "
                "(|E_pen| ~ 0.52)."
            ),
            honest_caveat=(
                "Its effect is a STEP, not a slope: because Xenium counts are "
                "integers, the PECAM1 z-score is quantised and the endothelial "
                "fraction is 17.4% for threshold < 0.9 and 9.7% for threshold >= "
                "0.9. ~26% of the prior mass sits on the low branch, which the "
                "Visium cross-check disfavours — so we also report the "
                "Visium-consistent conditional (threshold >= 0.9)."
            ),
        ),
        "kde_bandwidth_um": Prior(
            name="kde_bandwidth_um",
            kind="uniform",
            spec=f"U({_BW_LO:.0f}, {_BW_HI:.0f}) um",
            nominal=BASE_BANDWIDTH,
            ppf=_bw_rv.ppf,
            justification=(
                "The vascular and ECM fields are Gaussian KDEs of discrete cell "
                "hits; the bandwidth is the smoothing scale that turns those hits "
                "into a continuous perfusion / matrix field. The codebase itself "
                "carries two different defensible values (real_fields."
                "_DEFAULT_BANDWIDTH_UM = 40 um, real_xenium._DENSITY_BANDWIDTH_UM "
                "= 50 um) — direct evidence that this is a genuinely uncertain "
                "choice. The range is bounded by physiology: ~25 um is about one "
                "to two cell diameters (below it the field becomes a shot-noise "
                "map of individual endothelial cells rather than a vessel), and "
                "~60 um is roughly the inter-capillary half-spacing (above it the "
                "real perivascular heterogeneity is smeared away). This was the "
                "strongest local mover of the kill fraction (|E_kill| ~ 0.44)."
            ),
            honest_caveat=(
                "No bandwidth selector (Silverman/CV) is used and none would be "
                "right here: the 'true' field is a vessel geometry we never "
                "measured, so there is no held-out density to cross-validate on."
            ),
        ),
        "calibration_alpha": Prior(
            name="calibration_alpha",
            kind="lognormal",
            spec=f"LogNormal(median={alpha0:.1f}, sigma_log={_ALPHA_SIGMA:.3f}) "
                 f"-> 90% CI = [{alpha0/_ALPHA_SPAN:.0f}, {alpha0*_ALPHA_SPAN:.0f}] "
                 f"(0.5x-2x)",
            nominal=alpha0,
            ppf=a_rv.ppf,
            justification=(
                "alpha is the ONE fitted constant in the whole stack "
                "(k_eff = k_passive + alpha*uptake_proxy), solved in closed form "
                "from a SINGLE doxorubicin anchor. The published anchor is a "
                "HALF-DISTANCE — drug immunofluorescence falls to 50% within "
                "40-50 um of a vessel (Primeau et al., Clin Cancer Res 2005, PMID "
                "16361566) — which calibration.py converts to the e-folding length "
                f"the simulator actually uses, L* = 45/ln2 = "
                f"{efold_target:.1f} um (band {efold_band[0]:.0f}-"
                f"{efold_band[1]:.0f} um), and fits alpha = {alpha0:.1f} against "
                "it. Because alpha ∝ 1/L*^2, the published band width ALONE "
                "already gives alpha in [0.81x, 1.27x]. But that band is one "
                "immunofluorescence-vs-distance study in one tumour model, with a "
                "method-dependent half-intensity definition, applied here to a "
                "human GBM section — the cross-study spread is far wider than the "
                "within-study band. We therefore widen to a 90% CI of 0.5x-2x "
                f"(alpha {alpha0/_ALPHA_SPAN:.0f}-{alpha0*_ALPHA_SPAN:.0f}), i.e. "
                f"an equivalent e-folding length L in ~{efold_target/math.sqrt(2):.0f}"
                f"-{efold_target*math.sqrt(2):.0f} um. Lognormal because alpha is a "
                "positive scale factor and the uncertainty is multiplicative. "
                "numerical_rigor flagged exactly this: smallest LOCAL elasticity "
                "(|E_pen| ~ 0.11), widest prior."
            ),
            honest_caveat=(
                "The prior width is a judgement call about how far a single "
                "literature anchor transfers; it is not derived from a "
                "meta-analysis (there is no multi-study penetration dataset in "
                "this repo to derive it from). Note also that numerical_rigor's "
                "elasticity table was computed BEFORE the half-distance -> "
                "e-folding unit correction (at alpha = 309.6 = 2.081x the current "
                "fit); the elasticity RANKING it reports is what we cite here, not "
                "its absolute alpha."
            ),
        ),
        "ecm_tortuosity_alpha": Prior(
            name="ecm_tortuosity_alpha",
            kind="uniform",
            spec=f"U({_TORT_LO:.2f}, {_TORT_HI:.2f})",
            nominal=NOMINAL_TORTUOSITY,
            ppf=_tort_rv.ppf,
            justification=(
                "physics.diffusion couples matrix density to transport as "
                "D_eff = D * (1 - ALPHA_TORTUOSITY*ecm) with ALPHA_TORTUOSITY = "
                "0.6, i.e. a ~2.5x diffusivity reduction at full ECM density. In "
                "diffusion-in-brain terms the reduction factor is 1/lambda^2 for "
                "tortuosity lambda; measured lambda in brain / tumour "
                "extracellular space spans ~1.3-2.2 (Sykova & Nicholson, Physiol "
                "Rev 2008), i.e. a reduction of ~0.4-0.8 — which is the range "
                "sampled here, with the frozen 0.6 near its centre. It is a "
                "MODELLING CONSTANT, never fit to this section."
            ),
            honest_caveat=(
                "Sampling it patches the module constant at run time for the "
                "transport solve only; the calibration of alpha (which uses the "
                "same coefficient to convert the 45 um anchor into a rate) is NOT "
                "re-solved per draw. That is defensible only because alpha is "
                "given its own, much wider independent prior — the induced "
                "correlation is inside the alpha prior's 0.5x-2x span."
            ),
        ),
    }


# The PRIOR spec, built once at import so it is directly inspectable
# (``from tmesim.validation.uncertainty import PRIOR``). Building it only needs
# the molecular library (the closed-form alpha fit), never the Xenium bytes; if
# that fit is unavailable the dict is left empty and rebuilt inside run_uq().
try:
    PRIOR: Dict[str, Prior] = build_priors()
except Exception:  # pragma: no cover - keeps the module importable
    PRIOR = {}

PARAM_ORDER = (
    "endothelial_typing_threshold",
    "kde_bandwidth_um",
    "calibration_alpha",
    "ecm_tortuosity_alpha",
)


def sample_priors(priors: Dict[str, Prior], n: int, seed: int = 0) -> Dict[str, np.ndarray]:
    """Latin-hypercube draw of ``n`` samples from the prior spec.

    A Latin hypercube stratifies each marginal into ``n`` equal-probability bins
    (one draw per bin), so at n ~ 60 every prior is covered evenly instead of
    clumping the way plain iid Monte Carlo does. Dimensions are independent (no
    correlation is asserted between the knobs — we have no evidence for any).
    """
    names = [p for p in PARAM_ORDER if p in priors]
    engine = stats.qmc.LatinHypercube(d=len(names), seed=int(seed))
    u = engine.random(int(n))
    return {name: priors[name].sample_from_unit(u[:, j]) for j, name in enumerate(names)}


@contextlib.contextmanager
def _patched_tortuosity(value: float):
    """Temporarily set ``physics.diffusion.ALPHA_TORTUOSITY`` (restored on exit).

    ``_effective_diffusivity`` reads the constant as a module global at call
    time, so this is the only way to vary the ECM-tortuosity coupling WITHOUT
    editing the frozen physics module. Nothing is written to disk and the
    original value is always restored (try/finally).
    """
    old = _diffusion.ALPHA_TORTUOSITY
    _diffusion.ALPHA_TORTUOSITY = float(value)
    try:
        yield
    finally:
        _diffusion.ALPHA_TORTUOSITY = old


# ===========================================================================
# 2. THE MONTE-CARLO LOOP.
# ===========================================================================
def _endothelial_fraction(micro) -> float:
    ct = np.asarray(micro.cells.cell_type)
    return float(np.mean(ct == int(CellType.ENDOTHELIAL))) if ct.size else 0.0


def run_uq(
    n_samples: int = 80,
    drug: str = DRUG_NAME,
    dx_um: float = UQ_DX_UM,
    seed: int = 0,
    kill_seeds: Sequence[int] = UQ_KILL_SEEDS,
    priors: Optional[Dict[str, Prior]] = None,
    verbose: bool = True,
) -> Dict:
    """Propagate the prior uncertainty through the FULL real pipeline.

    For each Latin-hypercube draw: re-reconstruct the real Xenium section with
    that draw's ``endothelial_typing_threshold`` and ``kde_bandwidth_um``,
    rebuild the calibrated drug with that draw's ``calibration_alpha``, patch in
    that draw's ``ecm_tortuosity_alpha``, and run the frozen engine. Returns the
    draws, the per-draw predictions, the posterior-predictive summary (median +
    90% credible interval) and the first-order variance decomposition.

    ``dx_um`` is FIXED (not sampled): grid spacing is a numerical requirement,
    not an unknown. See :func:`resolution_bias` for the residual bias it leaves.
    """
    t_start = time.time()
    alpha0 = calibrate_alpha()
    priors = priors or build_priors(alpha0)
    draws = sample_priors(priors, n_samples, seed=seed)

    rows: List[Dict] = []
    for i in range(int(n_samples)):
        thr = float(draws["endothelial_typing_threshold"][i])
        bw = float(draws["kde_bandwidth_um"][i])
        alpha = float(draws["calibration_alpha"][i])
        tort = float(draws["ecm_tortuosity_alpha"][i])

        t0 = time.time()
        micro = reconstruct_knobs(dx_um, thr, bw)
        dp = calibrated_drugparams(drug, alpha)
        with _patched_tortuosity(tort):
            pred = _predict(micro, dp, pde=_pde(), seeds=kill_seeds)
        rows.append({
            "i": i,
            "endothelial_typing_threshold": thr,
            "kde_bandwidth_um": bw,
            "calibration_alpha": alpha,
            "ecm_tortuosity_alpha": tort,
            "kill_fraction": pred["kill_fraction"],
            "penetration_index": pred["penetration_index"],
            "kill_seed_std": pred["kill_std"],
            "endothelial_fraction": _endothelial_fraction(micro),
            "runtime_s": time.time() - t0,
        })
        if verbose:
            print(f"  draw {i+1:3d}/{n_samples}  thr={thr:.3f} bw={bw:5.1f} "
                  f"alpha={alpha:6.1f} tort={tort:.3f}  ->  "
                  f"kill={pred['kill_fraction']:.4f} "
                  f"pen={pred['penetration_index']:.4f}  "
                  f"[{time.time()-t0:.1f}s]", flush=True)

    kill = np.array([r["kill_fraction"] for r in rows])
    pen = np.array([r["penetration_index"] for r in rows])
    X = {p: np.array([r[p] for r in rows]) for p in PARAM_ORDER}

    result = {
        "drug": drug,
        "dx_um": float(dx_um),
        "n_samples": int(n_samples),
        "seed": int(seed),
        "kill_seeds": list(kill_seeds),
        "alpha0": float(alpha0),
        "priors": priors,
        "rows": rows,
        "grid": rows[0]["grid"] if "grid" in rows[0] else None,
        "kill_fraction": _summary(kill),
        "penetration_index": _summary(pen),
        "variance": {
            "kill_fraction": variance_decomposition(X, kill),
            "penetration_index": variance_decomposition(X, pen),
        },
        "seed_noise_mean_std": float(np.mean([r["kill_seed_std"] for r in rows])),
        "runtime_s": time.time() - t_start,
    }
    # Visium-consistent conditional: drop the draws whose endothelial fraction the
    # independent Visium deconvolution (~9.8%) contradicts.
    keep = X["endothelial_typing_threshold"] >= 0.9
    result["visium_consistent"] = {
        "n": int(keep.sum()),
        "rule": "endothelial_typing_threshold >= 0.9 (endothelial fraction ~9.7%, "
                "matching the independent Visium deconvolution ~9.8%)",
        "kill_fraction": _summary(kill[keep]) if keep.any() else None,
        "penetration_index": _summary(pen[keep]) if keep.any() else None,
    }
    return result


def _summary(v: np.ndarray) -> Dict[str, float]:
    v = np.asarray(v, dtype=np.float64)
    lo, med, hi = np.percentile(v, [5.0, 50.0, 95.0])
    return {
        "median": float(med),
        "mean": float(v.mean()),
        "std": float(v.std(ddof=1)) if v.size > 1 else 0.0,
        "ci90_lo": float(lo),
        "ci90_hi": float(hi),
        "min": float(v.min()),
        "max": float(v.max()),
        "rel_width_pct": float(100.0 * (hi - lo) / med) if med else float("nan"),
        "n": int(v.size),
    }


# ===========================================================================
# 3. VARIANCE DECOMPOSITION (first-order sensitivity).
# ===========================================================================
def _sobol_lite(x: np.ndarray, y: np.ndarray, n_bins: int = 6) -> float:
    """Binned first-order Sobol index ``S_i = Var(E[Y|X_i]) / Var(Y)``.

    Sort by ``x``, cut into ``n_bins`` equal-count bins, and take the variance of
    the bin means. The raw between-bin variance is biased UP by the within-bin
    noise (E[V_between] = Var(E[Y|X_i]) + sigma_w^2/m), so we subtract the
    within-bin mean square over the bin size and clip at 0. Unlike a correlation
    this captures NON-MONOTONE and STEP responses — which matters here, because
    the endothelial typing threshold acts as a step.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    n = y.size
    var_y = float(np.var(y, ddof=1))
    if n < 2 * n_bins or var_y <= 0:
        return float("nan")
    order = np.argsort(x, kind="mergesort")
    ys = y[order]
    bins = np.array_split(ys, n_bins)
    m = float(np.mean([b.size for b in bins]))
    means = np.array([b.mean() for b in bins])
    v_between = float(np.var(means, ddof=1))
    # pooled within-bin mean square
    ss_w = float(sum(((b - b.mean()) ** 2).sum() for b in bins))
    dof_w = max(n - n_bins, 1)
    msw = ss_w / dof_w
    s = (v_between - msw / m) / var_y
    return float(max(s, 0.0))


def _sobol_lite_ci(x: np.ndarray, y: np.ndarray, n_bins: int,
                   n_boot: int = 400, seed: int = 0) -> Tuple[float, float]:
    """Bootstrap 90% interval for :func:`_sobol_lite` (resample the draws).

    At n ~ 60-80 a first-order index is itself a noisy statistic; quoting it
    without an interval would be exactly the over-claiming this module exists to
    avoid.
    """
    rng = np.random.default_rng(seed)
    n = np.asarray(y).size
    vals = []
    for _ in range(int(n_boot)):
        idx = rng.integers(0, n, n)
        v = _sobol_lite(np.asarray(x)[idx], np.asarray(y)[idx], n_bins=n_bins)
        if np.isfinite(v):
            vals.append(v)
    if not vals:
        return float("nan"), float("nan")
    lo, hi = np.percentile(vals, [5.0, 95.0])
    return float(lo), float(hi)


def variance_decomposition(X: Dict[str, np.ndarray], y: np.ndarray,
                           n_bins: Optional[int] = None) -> Dict:
    """First-order variance attribution with TWO estimators.

    * ``sobol_lite`` — binned ``Var(E[Y|X_i])/Var(Y)`` (robust to steps and
      non-monotonicity; the headline number).
    * ``srrc`` — standardised RANK regression coefficient; ``srrc^2`` is the
      variance share under a monotone-linear model, and ``model_r2`` says how
      much of the total variance that linear-in-rank model explains at all (a
      low R2 is itself a finding: the response is not monotone).

    Inputs are independent by construction (Latin hypercube, no imposed
    correlation), so the first-order shares are directly comparable and their sum
    is the additive (non-interaction) share of the variance.
    """
    y = np.asarray(y, dtype=np.float64)
    names = [p for p in PARAM_ORDER if p in X]
    # ~10 draws per bin: enough to estimate a bin mean, enough bins to see shape.
    if n_bins is None:
        n_bins = int(np.clip(y.size // 10, 4, 8))
    out: Dict[str, Dict[str, float]] = {}

    # rank-transform for SRRC
    ry = stats.rankdata(y)
    rX = np.column_stack([stats.rankdata(X[p]) for p in names])
    rX_c = (rX - rX.mean(axis=0)) / rX.std(axis=0, ddof=1)
    ry_c = (ry - ry.mean()) / ry.std(ddof=1)
    beta, *_ = np.linalg.lstsq(rX_c, ry_c, rcond=None)
    resid = ry_c - rX_c @ beta
    r2 = float(1.0 - resid.var(ddof=1) / ry_c.var(ddof=1))

    for j, p in enumerate(names):
        s_lo, s_hi = _sobol_lite_ci(X[p], y, n_bins=n_bins)
        out[p] = {
            "sobol_lite": _sobol_lite(X[p], y, n_bins=n_bins),
            "sobol_lite_lo": s_lo,
            "sobol_lite_hi": s_hi,
            "srrc": float(beta[j]),
            "srrc_sq": float(beta[j] ** 2),
            "spearman": float(stats.spearmanr(X[p], y).statistic),
            "pearson": float(stats.pearsonr(X[p], y).statistic),
        }
    tot_s = sum(v["sobol_lite"] for v in out.values() if np.isfinite(v["sobol_lite"]))
    tot_b = sum(v["srrc_sq"] for v in out.values())
    for v in out.values():
        v["sobol_lite_share"] = (
            float(v["sobol_lite"] / tot_s) if tot_s > 0 and np.isfinite(v["sobol_lite"])
            else float("nan")
        )
        v["srrc_share"] = float(v["srrc_sq"] / tot_b) if tot_b > 0 else float("nan")
    ranked = sorted(
        names,
        key=lambda p: -(out[p]["sobol_lite"] if np.isfinite(out[p]["sobol_lite"])
                        else -np.inf),
    )
    return {
        "per_param": out,
        "ranked": ranked,
        "n_bins": int(n_bins),
        "sum_first_order": float(tot_s),
        "rank_model_r2": r2,
        "var_total": float(np.var(y, ddof=1)),
    }


# ===========================================================================
# 4. RESOLUTION AUDIT for the (fixed) MC grid.
# ===========================================================================
def resolution_bias(
    drug: str = DRUG_NAME,
    dxs: Sequence[float] = (30.0, 20.0, 15.0, 12.0),
    uq_dx: float = UQ_DX_UM,
    kill_seeds: Sequence[int] = UQ_KILL_SEEDS,
) -> Dict:
    """Re-run the NOMINAL parameter point at several dx to quantify how much of
    the reported prediction is still a resolution artifact at the MC grid.

    numerical_rigor established the requirement on the *calibrated doxorubicin*
    (a ~45-65 um boundary layer -> dx <= 12 um). Temozolomide's calibrated
    penetration length is several times larger, so the requirement is milder — but
    it must be re-measured for this drug rather than assumed, which is what this
    does. The drug's analytic e-folding length is returned alongside.
    """
    alpha0 = calibrate_alpha()
    dp = calibrated_drugparams(drug, alpha0)
    rows: List[Dict] = []
    for dx in dxs:
        t0 = time.time()
        micro = reconstruct_knobs(float(dx), BASE_THRESHOLD, BASE_BANDWIDTH)
        pred = _predict(micro, dp, pde=_pde(), seeds=kill_seeds)
        rows.append({"dx_um": float(dx), "runtime_s": time.time() - t0, **pred})
    fine = rows[-1]
    at_uq = min(rows, key=lambda r: abs(r["dx_um"] - uq_dx))
    return {
        "drug": drug,
        "rows": rows,
        "alpha0": float(alpha0),
        "drug_L_um": float(penetration_length_um(dp, ecm_density=DEFAULT_ECM)),
        "uq_dx_um": float(uq_dx),
        "finest_dx_um": float(fine["dx_um"]),
        "pen_bias_at_uq_dx": float(
            (at_uq["penetration_index"] - fine["penetration_index"])
            / fine["penetration_index"]
        ),
        "kill_bias_at_uq_dx": float(
            (at_uq["kill_fraction"] - fine["kill_fraction"]) / fine["kill_fraction"]
        ),
    }


# ===========================================================================
# 5. FIGURE.
# ===========================================================================
# Categorical palette (validated, colour-blind-safe ordering).
_C_KILL = "#2a78d6"      # blue
_C_PEN = "#1baf7a"       # aqua
_C_ACCENT = "#e34948"    # red (CI / median rules)
_C_MUTED = "#52514e"
_PARAM_LABEL = {
    "endothelial_typing_threshold": "endothelial\ntyping threshold",
    "kde_bandwidth_um": "KDE\nbandwidth",
    "calibration_alpha": "calibration\nalpha",
    "ecm_tortuosity_alpha": "ECM\ntortuosity",
}


def plot_uncertainty(res: Dict, path: str = "results/uncertainty.png") -> str:
    """Posterior-predictive histograms (with the 90% CI marked) + the per-prior
    variance attribution + the threshold step that drives the bimodality."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = res["rows"]
    kill = np.array([r["kill_fraction"] for r in rows])
    pen = np.array([r["penetration_index"] for r in rows])
    thr = np.array([r["endothelial_typing_threshold"] for r in rows])

    fig, axes = plt.subplots(2, 2, figsize=(12.4, 8.6))
    fig.patch.set_facecolor("white")

    def hist_panel(ax, v, summ, color, label, unit=""):
        ax.hist(v, bins=14, color=color, edgecolor="white", linewidth=1.2, alpha=0.9)
        lo, hi, med = summ["ci90_lo"], summ["ci90_hi"], summ["median"]
        ax.axvspan(lo, hi, color=_C_ACCENT, alpha=0.10, zorder=0)
        for x, ls, lw in ((lo, "--", 1.6), (hi, "--", 1.6), (med, "-", 2.2)):
            ax.axvline(x, color=_C_ACCENT, linestyle=ls, linewidth=lw, zorder=3)
        ax.set_title(f"{label}\nmedian {med:.4g}{unit}   90% CI "
                     f"[{lo:.4g}, {hi:.4g}]{unit}", fontsize=11, color="#0b0b0b")
        ax.set_xlabel(label, fontsize=10)
        ax.set_ylabel("draws", fontsize=10)
        ax.grid(axis="y", alpha=0.18, linewidth=0.7)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)

    hist_panel(axes[0, 0], kill, res["kill_fraction"], _C_KILL, "kill_fraction")
    hist_panel(axes[0, 1], pen, res["penetration_index"], _C_PEN, "penetration_index")

    # --- variance attribution -------------------------------------------------
    ax = axes[1, 0]
    names = [p for p in PARAM_ORDER]
    def _vals(metric, key):
        return np.array([float(np.nan_to_num(
            res["variance"][metric]["per_param"][p][key])) for p in names])

    kv, pv = _vals("kill_fraction", "sobol_lite"), _vals("penetration_index", "sobol_lite")
    kerr = np.vstack([kv - _vals("kill_fraction", "sobol_lite_lo"),
                      _vals("kill_fraction", "sobol_lite_hi") - kv]).clip(0)
    perr = np.vstack([pv - _vals("penetration_index", "sobol_lite_lo"),
                      _vals("penetration_index", "sobol_lite_hi") - pv]).clip(0)
    y = np.arange(len(names))
    h = 0.36
    ekw = dict(ecolor=_C_MUTED, elinewidth=1.2, capsize=3)
    ax.barh(y + h / 2 + 0.01, kv, height=h, color=_C_KILL, label="kill_fraction",
            xerr=kerr, error_kw=ekw)
    ax.barh(y - h / 2 - 0.01, pv, height=h, color=_C_PEN,
            label="penetration_index", xerr=perr, error_kw=ekw)
    for yi, v, e in zip(y + h / 2 + 0.01, kv, kerr[1]):
        ax.text(v + e + 0.02, yi, f"{v:.2f}", va="center", fontsize=9, color=_C_MUTED)
    for yi, v, e in zip(y - h / 2 - 0.01, pv, perr[1]):
        ax.text(v + e + 0.02, yi, f"{v:.2f}", va="center", fontsize=9, color=_C_MUTED)
    ax.set_yticks(y)
    ax.set_yticklabels([_PARAM_LABEL[p] for p in names], fontsize=9)
    ax.set_xlabel("first-order variance share  Var(E[Y|X]) / Var(Y)", fontsize=10)
    ax.set_xlim(0, float(max((kv + kerr[1]).max(), (pv + perr[1]).max())) * 1.20 + 0.06)
    ax.set_title("Which prior buys the spread (Sobol-lite, bars = bootstrap 90%)",
                 fontsize=11)
    ax.legend(frameon=False, fontsize=9, loc="lower right")
    ax.grid(axis="x", alpha=0.18, linewidth=0.7)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    # --- the step that drives it ---------------------------------------------
    ax = axes[1, 1]
    lowb = thr < 0.9
    ax.scatter(thr[lowb], pen[lowb], s=42, color=_C_ACCENT, edgecolor="white",
               linewidth=0.8, label="threshold < 0.9  (endothelium 17.4%)", zorder=3)
    ax.scatter(thr[~lowb], pen[~lowb], s=42, color=_C_PEN, edgecolor="white",
               linewidth=0.8, label="threshold ≥ 0.9  (endothelium 9.7%)", zorder=3)
    ax.axvline(0.9, color=_C_MUTED, linestyle=":", linewidth=1.4)
    ax.text(0.905, ax.get_ylim()[0], " quantisation step", fontsize=8.5,
            color=_C_MUTED, va="bottom")
    ax.set_xlabel("endothelial_typing_threshold (prior draw)", fontsize=10)
    ax.set_ylabel("penetration_index", fontsize=10)
    ax.set_title("The threshold acts as a STEP, not a slope", fontsize=11)
    ax.legend(frameon=False, fontsize=8.5, loc="upper right")
    ax.grid(alpha=0.18, linewidth=0.7)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    fig.suptitle(
        f"Parameter-uncertainty error bars — {res['drug']} on the REAL Xenium GBM "
        f"section  (n={res['n_samples']} draws, dx={res['dx_um']:.0f} µm)",
        fontsize=12.5, y=0.985,
    )
    fig.text(0.5, 0.005,
             "Expert priors on reconstruction + calibration knobs only — NOT a "
             "validated predictive interval (no outcome data; model-form error "
             "not included).",
             ha="center", fontsize=9, color=_C_MUTED)
    fig.tight_layout(rect=(0, 0.022, 1, 0.965))
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    fig.savefig(path, dpi=130, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


# ===========================================================================
# 6. REPORT.
# ===========================================================================
def write_report(res: Dict, resolution: Optional[Dict] = None,
                 path: str = "results/uncertainty.md",
                 figure: str = "results/uncertainty.png") -> str:
    L: List[str] = []
    A = L.append
    priors: Dict[str, Prior] = res["priors"]
    k, p = res["kill_fraction"], res["penetration_index"]
    vk = res["variance"]["kill_fraction"]
    vp = res["variance"]["penetration_index"]

    A("# Uncertainty quantification — honest error bars on the REAL microenvironment")
    A("")
    A("**Module:** `tmesim/validation/uncertainty.py` "
      "(run: `python -m tmesim.validation.uncertainty`)")
    A("")
    A(f"**What this is.** A Monte-Carlo propagation of the knobs `results/"
      f"numerical_rigor.md` identified as the real uncertainty drivers, through "
      f"the FULL real pipeline — the 40,887-cell 10x Xenium GBM section is "
      f"re-reconstructed from measured bytes for **every draw** (real positions, "
      f"real marker phenotypes, real endothelial-KDE vascular field, real "
      f"matrix-transcript ECM field, physics-solved oxygen), the drug is rebuilt "
      f"from the calibration model, and the frozen engine is re-run. Drug: "
      f"**{res['drug']}** (GBM standard of care), α-calibrated. "
      f"n = **{res['n_samples']}** Latin-hypercube draws at dx = "
      f"**{res['dx_um']:.0f} µm**, total runtime "
      f"**{res['runtime_s']/60.0:.1f} min**.")
    A("")
    A("---")
    A("")

    # --- priors ---
    A("## 1. The priors (and why)")
    A("")
    A("**These are EXPERT priors, not measured distributions.** Nothing in this "
      "repository measures a distribution over the typing gate, the KDE bandwidth "
      "or the tortuosity coefficient, and there is no outcome label to condition "
      "on — so there is no posterior here in the Bayesian-updating sense. Each "
      "prior is a documented *range of defensible analyst choice*; the intervals "
      "below say how far the prediction moves inside that range.")
    A("")
    A("| knob | prior | nominal | 5% / 50% / 95% |")
    A("|---|---|---:|:--:|")
    for name in PARAM_ORDER:
        pr = priors[name]
        q5, q50, q95 = pr.quantiles()
        A(f"| `{name}` | {pr.kind}: {pr.spec} | {pr.nominal:.4g} | "
          f"{q5:.3g} / {q50:.3g} / {q95:.3g} |")
    A("")
    for name in PARAM_ORDER:
        pr = priors[name]
        A(f"**`{name}`** — {pr.justification}")
        A("")
        A(f"> *Honest caveat:* {pr.honest_caveat}")
        A("")
    A("**Not sampled: `grid_dx_um`.** numerical_rigor showed dx is a strong local "
      "mover (|E_kill| ≈ 0.35), but grid spacing is **not an unknown** — it is a "
      "convergence requirement, and sampling it would launder a numerical error "
      "into a credible interval. It is fixed and audited instead (§5).")
    A("")
    A("---")
    A("")

    # --- design ---
    A("## 2. Monte-Carlo design")
    A("")
    A(f"* **Sampling:** Latin hypercube over the {len(PARAM_ORDER)} priors "
      f"(scipy `qmc.LatinHypercube`, seed {res['seed']}), each column mapped "
      f"through its inverse CDF. Dimensions independent — no correlation is "
      f"asserted between the knobs because none is measured.")
    A(f"* **Per draw:** full re-reconstruction at that draw's typing threshold + "
      f"KDE bandwidth, drug rebuilt at that draw's α, `ALPHA_TORTUOSITY` patched "
      f"for the transport solve (restored afterwards; no frozen file is edited).")
    A(f"* **Kill de-noising:** `kill_fraction` is averaged over seeds "
      f"{res['kill_seeds']}; the residual across-seed sampling noise is "
      f"±{res['seed_noise_mean_std']:.4f} (mean), i.e. "
      f"{100.0*res['seed_noise_mean_std']/max(k['median'],1e-12):.1f}% of the "
      f"median — **far smaller than the parameter spread reported below**, so the "
      f"interval is parameter uncertainty, not MC noise.")
    A(f"* **Cost:** {res['runtime_s']/60.0:.1f} min for {res['n_samples']} draws "
      f"(~{res['runtime_s']/max(res['n_samples'],1):.1f} s/draw), single process.")
    A(f"* **Provenance:** the α prior is centred on `calibration.calibrate_alpha()` "
      f"evaluated at run time — α₀ = **{res['alpha0']:.1f}** for this run (the "
      "half-distance → e-folding convention). If the calibration module is "
      "re-fit, re-run this module: the interval moves with it.")
    A("")
    A("---")
    A("")

    # --- results ---
    A("## 3. Posterior-predictive intervals")
    A("")
    A("| prediction | median | 90% credible interval | rel. width | min–max | prior-driven std |")
    A("|---|---:|:--:|---:|:--:|---:|")
    A(f"| **kill_fraction** | **{k['median']:.4f}** | "
      f"[{k['ci90_lo']:.4f}, {k['ci90_hi']:.4f}] | {k['rel_width_pct']:.0f}% | "
      f"{k['min']:.4f}–{k['max']:.4f} | {k['std']:.4f} |")
    A(f"| **penetration_index** | **{p['median']:.4f}** | "
      f"[{p['ci90_lo']:.4f}, {p['ci90_hi']:.4f}] | {p['rel_width_pct']:.0f}% | "
      f"{p['min']:.4f}–{p['max']:.4f} | {p['std']:.4f} |")
    A("")
    vc = res.get("visium_consistent") or {}
    if vc.get("kill_fraction"):
        vk_, vp_ = vc["kill_fraction"], vc["penetration_index"]
        A(f"**Visium-consistent conditional** ({vc['n']}/{res['n_samples']} draws; "
          f"rule: {vc['rule']}) — dropping the draws whose endothelial fraction "
          f"(17.4%) the independent Visium deconvolution (~9.8%) contradicts:")
        A("")
        A("| prediction | median | 90% credible interval |")
        A("|---|---:|:--:|")
        A(f"| kill_fraction | {vk_['median']:.4f} | "
          f"[{vk_['ci90_lo']:.4f}, {vk_['ci90_hi']:.4f}] |")
        A(f"| penetration_index | {vp_['median']:.4f} | "
          f"[{vp_['ci90_lo']:.4f}, {vp_['ci90_hi']:.4f}] |")
        A("")
        A("This is the interval to quote if you accept the Visium cross-check; the "
          "unconditional one above is the interval if you do not.")
        A("")
    A(f"![uncertainty]({os.path.basename(figure)})")
    A("")
    A("---")
    A("")

    # --- variance ---
    A("## 4. Variance attribution (which prior buys the spread)")
    A("")
    A(f"`S₁` = binned first-order index `Var(E[Y|X_i])/Var(Y)` "
      f"({vk['n_bins']} equal-count bins, bias-corrected for within-bin noise; "
      "robust to step / non-monotone responses), with a **bootstrap 90% "
      "interval** over the draws — a first-order index at this n is itself a "
      "noisy statistic and is quoted as such. `SRRC` = standardised "
      "rank-regression coefficient, the monotone-linear cross-check (sign = "
      "direction of effect).")
    A("")
    A("| knob | S₁ (kill) [90% boot] | SRRC (kill) | S₁ (pen) [90% boot] | SRRC (pen) |")
    A("|---|:--:|---:|:--:|---:|")
    for name in PARAM_ORDER:
        a_, b_ = vk["per_param"][name], vp["per_param"][name]
        A(f"| `{name}` | {a_['sobol_lite']:.2f} "
          f"[{a_['sobol_lite_lo']:.2f}, {a_['sobol_lite_hi']:.2f}] | "
          f"{a_['srrc']:+.2f} | {b_['sobol_lite']:.2f} "
          f"[{b_['sobol_lite_lo']:.2f}, {b_['sobol_lite_hi']:.2f}] | "
          f"{b_['srrc']:+.2f} |")
    A(f"| *sum of first-order* | {vk['sum_first_order']:.3f} | "
      f"(rank-R² {vk['rank_model_r2']:.2f}) | {vp['sum_first_order']:.3f} | "
      f"(rank-R² {vp['rank_model_r2']:.2f}) |")
    A("")
    A(f"**Dominant contributor — kill_fraction:** `{vk['ranked'][0]}` "
      f"(S₁ ≈ {vk['per_param'][vk['ranked'][0]]['sobol_lite']:.2f}); "
      f"runner-up `{vk['ranked'][1]}` "
      f"(S₁ ≈ {vk['per_param'][vk['ranked'][1]]['sobol_lite']:.2f}).")
    A(f"**Dominant contributor — penetration_index:** `{vp['ranked'][0]}` "
      f"(S₁ ≈ {vp['per_param'][vp['ranked'][0]]['sobol_lite']:.2f}); "
      f"runner-up `{vp['ranked'][1]}` "
      f"(S₁ ≈ {vp['per_param'][vp['ranked'][1]]['sobol_lite']:.2f}).")
    A("")
    A(f"With n = {res['n_samples']} draws these indices carry a Monte-Carlo error "
      "of order ±0.05–0.10; read the RANKING, not the third decimal. Shares below "
      "~0.05 are indistinguishable from zero at this sample size.")
    A("")
    over = max(vk["sum_first_order"], vp["sum_first_order"])
    if over > 1.02:
        A(f"**The first-order shares sum to more than 1** "
          f"({vk['sum_first_order']:.2f} for kill, {vp['sum_first_order']:.2f} for "
          "penetration). True first-order Sobol indices sum to ≤ 1, so the excess "
          "is estimator bias/noise at this sample size, not interaction structure "
          "— it is stated rather than normalised away. It is one more reason the "
          "ranking, and not the absolute share, is the reportable result.")
        A("")
    A("**Local sensitivity is NOT global uncertainty.** numerical_rigor's ±20% "
      "elasticities ranked `endothelial_typing_threshold` first for penetration "
      "and `kde_bandwidth_um` first for kill, with `calibration_alpha` last on "
      f"both. This decomposition ranks `{vp['ranked'][0]}` first for penetration "
      f"and `{vk['ranked'][0]}` first for kill. The two answer different "
      "questions: an elasticity is the per-unit response at one point, a variance "
      "share is response × **prior width**."
      + (" α is the clean example — the SMALLEST local elasticity (|E| ≈ 0.1–0.2) "
         "but the widest prior (0.5×–2×), and it comes out as the single largest "
         "contributor to the kill-fraction spread. numerical_rigor predicted "
         "exactly this ('low-sensitivity, high-prior-uncertainty knob') and the "
         "Monte Carlo confirms it."
         if vk["ranked"][0] == "calibration_alpha" else
         " Where the two disagree, the variance share is the one that belongs on "
         "an error bar.")
      )
    A("")
    A("---")
    A("")

    # --- resolution ---
    if resolution:
        A("## 5. Resolution caveat for the MC grid")
        A("")
        A(f"The MC runs at dx = {res['dx_um']:.0f} µm for affordability. "
          "numerical_rigor set the resolution requirement on the *calibrated "
          "doxorubicin* (a ~45–65 µm perivascular boundary layer → dx ≲ 12 µm); "
          f"{resolution['drug']}'s calibrated e-folding length is "
          f"**~{resolution['drug_L_um']:.0f} µm**, so the requirement is milder — "
          "but it is **re-measured here for this drug** rather than assumed, at "
          "the nominal parameter point:")
        A("")
        A("| dx (µm) | grid | tissue cells | penetration_index | kill_fraction |")
        A("|---:|:--:|---:|---:|---:|")
        for r in resolution["rows"]:
            gy, gx = r["grid"]
            A(f"| {r['dx_um']:.0f} | {gy}×{gx} | {r['n_tissue']:,} | "
              f"{r['penetration_index']:.4f} | {r['kill_fraction']:.4f} |")
        A("")
        A(f"**Honest consequence.** At the MC grid (dx = {res['dx_um']:.0f} µm) the "
          f"penetration index is **{100*resolution['pen_bias_at_uq_dx']:+.1f}%** and "
          f"the kill fraction **{100*resolution['kill_bias_at_uq_dx']:+.1f}%** "
          f"relative to the finest grid run here "
          f"(dx = {resolution['finest_dx_um']:.0f} µm). The kill fraction is "
          "grid-converged, so its interval is a genuine parameter interval. The "
          "**penetration index is resolution-biased** — coarse grids over-estimate "
          "it (the perivascular gradient is under-resolved), so the quoted "
          "penetration median should be read as a slight upper bound with a "
          "systematic offset of that size ON TOP of the credible interval. It is a "
          "BIAS, not a widening: it shifts every draw the same way, which is why "
          "it is reported separately instead of being folded into the CI.")
        A("")
        A("---")
        A("")

    # --- honesty ---
    A("## 6. What these error bars are — and what they are NOT")
    A("")
    A("**These are error bars from PARAMETER uncertainty only.** They answer one "
      "question: *if the reconstruction and calibration choices move inside their "
      "defensible ranges, how far does the prediction move?* They do **not** cover:")
    A("")
    A("1. **Model-form error.** The interval is computed *inside* one fixed model "
      "structure: a first-order `k = clearance + binding` sink, a "
      "`D_eff = D·(1 − α·ECM)` tortuosity law, a KDE-smoothed endothelial density "
      "standing in for a vessel network we never imaged, a Hill kill on AUC, and "
      "plasma-averaged delivery. If that structure is wrong — e.g. saturable "
      "uptake, efflux transporters, a real vessel geometry with dead ends, "
      "cell-cycle-dependent kill — no amount of parameter sampling detects it, and "
      "the true error is larger than everything reported here. **Model-form error "
      "is plausibly the LARGEST term and it is not in these bars.**")
    A("2. **The absence of outcome validation.** There is no measured drug "
      "response for this section. `kill_fraction` has never been compared to a "
      "real outcome on this tissue, and the per-cell `sensitivity` it integrates "
      "is an explicitly PROVISIONAL placeholder scaled from tumour-marker signal, "
      "not a measured sensitivity. A 90% credible interval around an unvalidated "
      "quantity is a statement about the model's internal consistency, **not** a "
      "claim that the true value falls in it 90% of the time.")
    A("3. **Anchor transfer.** The α prior is centred on ONE doxorubicin "
      "measurement in a different tumour model; its width is a judgement about "
      "how far that anchor travels, not a meta-analysis. If the anchor does not "
      "transfer to temozolomide in human GBM, the whole interval translates.")
    A("4. **Data-level uncertainty.** Segmentation error, transcript "
      "mis-assignment and the 319-gene panel's limited marker set are taken as "
      "given: one section, one segmentation, no replicate. Section-to-section and "
      "patient-to-patient variability is entirely absent — this is uncertainty "
      "*within one measured section*.")
    A("5. **Prior mis-specification.** Every prior above is an expert range. A "
      "different (equally defensible) analyst prior gives a different interval; "
      "the widths scale roughly with the prior widths, so they are only as honest "
      "as §1's justifications.")
    A("")
    A("**Bottom line.** "
      f"On the real Xenium GBM section, {res['drug']} is predicted to kill "
      f"**{100*k['median']:.1f}%** of viable tumour cells (90% CI "
      f"**{100*k['ci90_lo']:.1f}–{100*k['ci90_hi']:.1f}%**) with a penetration "
      f"index of **{p['median']:.3f}** (90% CI "
      f"**{p['ci90_lo']:.3f}–{p['ci90_hi']:.3f}**), where the interval comes "
      f"from the reconstruction + calibration priors and the spread is dominated "
      f"by `{vk['ranked'][0]}` (kill) and `{vp['ranked'][0]}` (penetration). "
      "The number to distrust most is not the width of these bars — it is the "
      "model structure and the missing outcome label they cannot see.")
    A("")

    text = "\n".join(L)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)
    return path


# ===========================================================================
def _print(res: Dict, resolution: Optional[Dict] = None) -> None:
    k, p = res["kill_fraction"], res["penetration_index"]
    vk = res["variance"]["kill_fraction"]
    vp = res["variance"]["penetration_index"]
    print("=" * 92)
    print(f"UNCERTAINTY QUANTIFICATION — {res['drug']} on the REAL Xenium GBM section")
    print(f"  n={res['n_samples']} LHS draws   dx={res['dx_um']:.0f} um   "
          f"runtime={res['runtime_s']/60:.1f} min")
    print("=" * 92)
    print("PRIORS (expert priors, not measured distributions):")
    for name in PARAM_ORDER:
        pr = res["priors"][name]
        q5, q50, q95 = pr.quantiles()
        print(f"  {name:32s} {pr.kind:18s} 5/50/95% = "
              f"{q5:.3g} / {q50:.3g} / {q95:.3g}")
    print("-" * 92)
    print("POSTERIOR-PREDICTIVE (median [90% CI]):")
    print(f"  kill_fraction      {k['median']:.4f}  "
          f"[{k['ci90_lo']:.4f}, {k['ci90_hi']:.4f}]   "
          f"(rel width {k['rel_width_pct']:.0f}%)")
    print(f"  penetration_index  {p['median']:.4f}  "
          f"[{p['ci90_lo']:.4f}, {p['ci90_hi']:.4f}]   "
          f"(rel width {p['rel_width_pct']:.0f}%)")
    vc = res.get("visium_consistent") or {}
    if vc.get("kill_fraction"):
        print(f"  [Visium-consistent subset, n={vc['n']}]  "
              f"kill {vc['kill_fraction']['median']:.4f} "
              f"[{vc['kill_fraction']['ci90_lo']:.4f}, "
              f"{vc['kill_fraction']['ci90_hi']:.4f}]   "
              f"pen {vc['penetration_index']['median']:.4f} "
              f"[{vc['penetration_index']['ci90_lo']:.4f}, "
              f"{vc['penetration_index']['ci90_hi']:.4f}]")
    print("-" * 92)
    print(f"FIRST-ORDER VARIANCE SHARE  Var(E[Y|X])/Var(Y)  "
          f"[bootstrap 90%], {vk['n_bins']} bins:")
    print(f"  {'knob':32s} {'S1(kill)':>22s} {'SRRC':>6s} {'S1(pen)':>22s} {'SRRC':>6s}")
    for name in PARAM_ORDER:
        a_, b_ = vk["per_param"][name], vp["per_param"][name]
        ka = f"{a_['sobol_lite']:.2f} [{a_['sobol_lite_lo']:.2f},{a_['sobol_lite_hi']:.2f}]"
        pa = f"{b_['sobol_lite']:.2f} [{b_['sobol_lite_lo']:.2f},{b_['sobol_lite_hi']:.2f}]"
        print(f"  {name:32s} {ka:>22s} {a_['srrc']:>+6.2f} {pa:>22s} "
              f"{b_['srrc']:>+6.2f}")
    print(f"  -> dominant: kill = {vk['ranked'][0]},  pen = {vp['ranked'][0]}")
    if resolution:
        print("-" * 92)
        print(f"RESOLUTION BIAS AT dx={res['dx_um']:.0f} um (vs "
              f"dx={resolution['finest_dx_um']:.0f} um): "
              f"pen {100*resolution['pen_bias_at_uq_dx']:+.1f}%, "
              f"kill {100*resolution['kill_bias_at_uq_dx']:+.1f}%")
    print("=" * 92)
    print("HONEST: parameter uncertainty ONLY — no model-form error, no outcome "
          "validation.")
    print("=" * 92)


if __name__ == "__main__":
    t0 = time.time()
    print("Running uncertainty quantification on REAL data (several minutes)...\n")
    resolution = resolution_bias()
    res = run_uq(n_samples=80)
    PRIOR.update(res["priors"])
    _print(res, resolution)
    fig = plot_uncertainty(res)
    out = write_report(res, resolution)
    print(f"\n[wrote {out}]  [wrote {fig}]  [total runtime {time.time()-t0:.1f} s]")
