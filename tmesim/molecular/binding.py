"""Structure-based estimates of protein binding, tissue retention, and clearance.

This module supplies the three pharmacokinetic scalars the tissue engine needs
that are NOT transport/diffusion:

  * plasma_free_fraction  -> the unbound fraction driving delivery
  * tissue_binding_rate   -> DrugParams.binding   (1/s, first-order retention)
  * clearance_rate        -> DrugParams.clearance (1/s, first-order metabolism)

CONFIDENCE / HONESTY NOTE (read this)
-------------------------------------
Every number here is a **low-confidence, structure-only estimate**. Real plasma
protein binding, tissue partitioning, and metabolic clearance depend on the
specific binding proteins (albumin/AAG), transporters, and CYP isoforms
involved, none of which are derivable from bulk physicochemical descriptors
alone. What we CAN defend is the *direction* and *rough magnitude* of the trends:

  - lipophilic molecules (high logP) bind plasma albumin more -> lower free
    fraction, and stick to tissue lipids/membranes more -> higher tissue binding;
  - hydrogen-bond donors add to non-specific tissue retention;
  - metabolic clearance is mildly higher for lipophilic small molecules (CYP
    substrates) and mildly lower for large molecules.

These are implemented as smooth, monotonic, well-behaved maps that stay inside
physically plausible ranges. Treat them as priors to be replaced by measured
PK values, not as predictions. The provenance strings in the compiled profile
should flag them as `heuristic (low confidence)`.

Units
-----
free_fraction        : dimensionless in [0.01, 1.0]
tissue_binding_rate  : 1/s   (kept in ~1e-5 .. 5e-3, matches DrugParams.binding)
clearance_rate       : 1/s   (kept in ~1e-5 .. 1e-3, matches DrugParams.clearance)
"""
from __future__ import annotations

import math

# ---------------------------------------------------------------------------
# Documented parameter ranges (also the clamp bounds), all defensible priors.
# ---------------------------------------------------------------------------
_FF_MIN, _FF_MAX = 0.01, 1.0        # free fraction, dimensionless
_BIND_MIN, _BIND_MAX = 1e-5, 5e-3   # tissue binding rate, 1/s
_CLR_MIN, _CLR_MAX = 1e-5, 1e-3     # clearance rate, 1/s


def _sigmoid(x: float) -> float:
    """Numerically-stable logistic sigmoid, maps R -> (0, 1)."""
    if x >= 0:
        z = math.exp(-x)
        return 1.0 / (1.0 + z)
    z = math.exp(x)
    return z / (1.0 + z)


def _geom_interp(lo: float, hi: float, s: float) -> float:
    """Geometric (log-linear) interpolation between lo and hi for s in [0,1].

    Returns lo*(hi/lo)**s, which stays in [lo, hi] for s in [0,1]. Using a
    geometric (not linear) blend keeps rates sensible across their 2-3 decade
    span: equal changes in the driving score move the rate by equal *factors*,
    which is the natural scale for first-order kinetic constants.
    """
    return lo * (hi / lo) ** s


def plasma_free_fraction(descriptors: dict) -> float:
    """Unbound fraction of drug in plasma, dimensionless in [0.01, 1.0].

    Heuristic (lipophilicity -> albumin binding)
    --------------------------------------------
    Highly lipophilic drugs bind plasma proteins (mainly albumin) strongly, so
    their *free* (pharmacologically active, tissue-delivering) fraction is low.
    We use a smooth logistic in logP:

        free_fraction = 1 / (1 + 10**((logP - k) / s))

    with midpoint k = 3.0 (free fraction = 0.5 at logP = 3) and slope s = 1.5
    (log10 units per e-fold-ish step). This reproduces the qualitative textbook
    picture: very hydrophilic drugs (logP <~ 0) are almost entirely free
    (~0.97-1.0), moderately lipophilic drugs (logP ~ 3) are about half bound,
    and very lipophilic drugs (logP >~ 5) are >99% bound. The result is clamped
    to [0.01, 1.0] so nothing is ever reported as fully unbound or fully bound.

    CONFIDENCE: low. logP is only a weak proxy for real protein binding (which
    is protein- and site-specific). Direction and rough magnitude only.
    """
    logp = float(descriptors.get("logp", 1.0))
    k, s = 3.0, 1.5
    ff = 1.0 / (1.0 + 10.0 ** ((logp - k) / s))
    return max(_FF_MIN, min(_FF_MAX, ff))


def tissue_binding_rate(descriptors: dict) -> float:
    """First-order tissue retention/binding rate, 1/s (feeds DrugParams.binding).

    Heuristic (lipophilic + H-bond donors stick to tissue)
    ------------------------------------------------------
    Non-specific tissue binding (to membrane lipids, ECM, intracellular protein)
    grows with lipophilicity and with hydrogen-bond donor count. We build a
    dimensionless "stickiness" score and squash it through a logistic, then map
    it geometrically onto the physically plausible band [1e-5, 5e-3] 1/s:

        score = (logP - 2.5) / 2.5  +  (HBD - 2) / 4
        s     = sigmoid(score)
        rate  = 1e-5 * (5e-3 / 1e-5) ** s          # in [1e-5, 5e-3] 1/s

    So a hydrophilic, low-donor molecule retains weakly (~2-3e-5 /s) while a
    lipophilic, multi-donor molecule retains strongly (~1e-3 /s). Monotonically
    increasing in both logP and HBD; always inside the documented range.

    The 1e-5..5e-3 /s band corresponds to retention timescales of roughly
    minutes (fast) to a day (slow), bracketing DrugParams.binding's default
    (5e-4 /s). CONFIDENCE: low; magnitudes are order-of-magnitude priors.
    """
    logp = float(descriptors.get("logp", 1.0))
    hbd = float(descriptors.get("hbd", 1.0))
    score = (logp - 2.5) / 2.5 + (hbd - 2.0) / 4.0
    rate = _geom_interp(_BIND_MIN, _BIND_MAX, _sigmoid(score))
    return max(_BIND_MIN, min(_BIND_MAX, rate))


def clearance_rate(descriptors: dict) -> float:
    """First-order tissue clearance/metabolism rate, 1/s (feeds DrugParams.clearance).

    Heuristic (mild size / lipophilicity dependence)
    ------------------------------------------------
    Local tissue clearance (metabolism + washout) is only weakly predictable
    from structure, so we keep the dependence deliberately MILD: lipophilic
    small molecules are somewhat better CYP substrates (faster metabolism),
    while large molecules are cleared/metabolised more slowly. Score and map:

        score = 0.4 * (logP - 2) / 3  -  0.4 * (MW - 300) / 700
        s     = sigmoid(score)
        rate  = 1e-5 * (1e-3 / 1e-5) ** s          # in [1e-5, 1e-3] 1/s

    The small (0.4) weights make the spread gentle: typical small molecules land
    near the middle of the band (~1e-4 /s, close to DrugParams.clearance's
    default of 2e-4 /s), lipophilic ones a bit higher, large molecules a bit
    lower. CONFIDENCE: low; this is a weak prior, not a metabolism prediction.
    """
    logp = float(descriptors.get("logp", 1.0))
    mw = float(descriptors.get("mw", 300.0))
    score = 0.4 * (logp - 2.0) / 3.0 - 0.4 * (mw - 300.0) / 700.0
    rate = _geom_interp(_CLR_MIN, _CLR_MAX, _sigmoid(score))
    return max(_CLR_MIN, min(_CLR_MAX, rate))


if __name__ == "__main__":
    print("=== binding self-test ===")

    # Contrasting descriptor dicts.
    lipophilic = {"name": "lipophilic_big", "mw": 540.0, "logp": 4.5, "hbd": 4}
    hydrophilic = {"name": "hydrophilic_small", "mw": 150.0, "logp": -1.0, "hbd": 1}
    midrange = {"name": "midrange", "mw": 300.0, "logp": 2.0, "hbd": 2}

    rows = []
    for d in (lipophilic, hydrophilic, midrange):
        ff = plasma_free_fraction(d)
        bind = tissue_binding_rate(d)
        clr = clearance_rate(d)
        rows.append((d["name"], ff, bind, clr))
        print(f"{d['name']:>18}  logP={d['logp']:>5}  HBD={d['hbd']}  MW={d['mw']:>6.0f}  "
              f"free_fraction={ff:.3f}  binding={bind:.2e} /s  clearance={clr:.2e} /s")

    # --- range checks (documented bands) ---
    for name, ff, bind, clr in rows:
        assert _FF_MIN <= ff <= _FF_MAX, f"{name}: free_fraction {ff} out of range"
        assert _BIND_MIN <= bind <= _BIND_MAX, f"{name}: binding {bind} out of range"
        assert _CLR_MIN <= clr <= _CLR_MAX, f"{name}: clearance {clr} out of range"
    print("range checks OK")

    ff_lip, bind_lip, clr_lip = rows[0][1], rows[0][2], rows[0][3]
    ff_hyd, bind_hyd, clr_hyd = rows[1][1], rows[1][2], rows[1][3]

    # --- directional checks (the core honest claims) ---
    assert ff_lip < ff_hyd, (
        f"lipophilic free_fraction ({ff_lip:.3f}) should be LOWER than "
        f"hydrophilic ({ff_hyd:.3f})")
    assert bind_lip > bind_hyd, (
        f"lipophilic tissue_binding ({bind_lip:.2e}) should be HIGHER than "
        f"hydrophilic ({bind_hyd:.2e})")
    print(f"directional checks OK: "
          f"free_fraction {ff_lip:.3f} < {ff_hyd:.3f}, "
          f"binding {bind_lip:.2e} > {bind_hyd:.2e}")

    # --- monotonicity spot-checks ---
    ff_series = [plasma_free_fraction({"logp": lp}) for lp in (-2, 0, 2, 4, 6)]
    assert all(ff_series[i] >= ff_series[i + 1] for i in range(len(ff_series) - 1)), \
        "free_fraction must be non-increasing in logP"
    bind_series = [tissue_binding_rate({"logp": lp, "hbd": 2}) for lp in (-2, 0, 2, 4, 6)]
    assert all(bind_series[i] <= bind_series[i + 1] for i in range(len(bind_series) - 1)), \
        "tissue_binding_rate must be non-decreasing in logP"
    print("monotonicity checks OK")

    print("ALL CHECKS PASSED")
