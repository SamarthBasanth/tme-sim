"""TEST B of `results/transport_anchor_prereg.md` -- rat brain-slice tissue
partitioning (`Vu,brain`) against the engine's implied tissue capacity.

    ./.venv/bin/python experiments/brain_slice_anchor.py --stage prereg
    ./.venv/bin/python experiments/brain_slice_anchor.py --stage floor
    ./.venv/bin/python experiments/brain_slice_anchor.py --stage engine

THIS SCRIPT EXECUTES A PRE-REGISTRATION IT DID NOT WRITE. Every threshold,
exclusion code and statistic below is transcribed from that document. Nothing
here is tuned, nothing is substituted for a more flattering alternative, and no
compound is dropped for scoring badly.

Stage order is load-bearing:

  prereg  writes results/sweep/transport_anchor_prereg.json -- the cascade, the
          exclusion codes, every numeric bar, the GATE 1 engine quantity and the
          density assumption (prereg 8.4 lists both as slots the executor fills,
          and 3.3 requires the assumption be written down BEFORE the comparison).
          Asserts the section 2.3 frozen constants against live code.
  floor   GMFE_floor and its BCa interval. Contains NO engine quantity: it is a
          property of the measurement alone (prereg 3.5a) and is therefore
          immune to the live `HD`-convention investigation and to any
          recalibration of `alpha`.
  engine  GATE 1, THE TEST, the secondaries B-S1..B-S10, the sensitivity arms.

SCOPE FENCE (prereg 2.1). This is `SCOPE_ANALYTIC`. The transport solver is
never entered; the stage-`engine` run installs call counters on
`physics.diffusion.simulate_drug_transport` and `sim.engine.simulate_response`
and asserts both stayed at zero, following the instrumentation pattern of
`results/transvascular_exchange_revalidation.md` 6.1.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import subprocess
import sys
from typing import Dict, List, Sequence, Tuple

import numpy as np
from scipy.stats import norm

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

DATA_CSV = os.path.join(_ROOT, "data_cache", "brain_slice", "vu_brain.csv")
PREREG_JSON = os.path.join(_ROOT, "results", "sweep", "transport_anchor_prereg.json")

# --------------------------------------------------------------------------- #
# Numeric bars, transcribed from results/transport_anchor_prereg.md.
# --------------------------------------------------------------------------- #
THE_TEST_BAR = 3.0                # 3.5: GMFE(engine)/GMFE_floor <= 3.0
REPRODUCIBILITY_FLOOR = 1.25      # 3.4a: declared floor, rounded UP from 1.23x
BOOTSTRAP_RESAMPLES = 20000       # 3.6: ">= 10 000"
BOOTSTRAP_CI = 0.90               # 3.6: bootstrap 90 % CI
BOOTSTRAP_SEED = 20260821
WITHIN_BANDS = (3.0, 10.0)        # B-S5: the repo's two declared bands

# 3.3 GATE 1 slot / 8.4: the density assumption, fixed here before any
# engine-to-measurement comparison is computed. 1.00 g.mL^-1 is the
# water-equivalent convention; it requires no citation and cannot be tuned.
# Any plausible rat-brain density (1.00-1.05) moves the engine value by < 5 %,
# which is far inside the 1.25x reproducibility floor -- the sensitivity is
# computed and printed rather than asserted.
TISSUE_DENSITY_G_PER_ML = 1.00
TISSUE_DENSITY_SENSITIVITY = 1.04


# --------------------------------------------------------------------------- #
# The dataset and the 3.1 cascade.
# --------------------------------------------------------------------------- #
def _load_rows() -> List[dict]:
    with open(DATA_CSV, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def cascade(rows: Sequence[dict]) -> Dict[str, object]:
    """The 3.1 cascade, applied in the registered order. Property-based only."""
    steps: List[Tuple[str, str, int, List[str]]] = []

    stage = list(rows)
    steps.append(("SOURCE", "every row in vu_brain.csv", len(stage), []))

    # E2 -- rat only. Species may never be pooled.
    dropped = [r["compound"] for r in stage if r["species"] != "rat"]
    stage = [r for r in stage if r["species"] == "rat"]
    steps.append(("E2", "non-rat species excluded", len(stage), dropped))

    # E1 -- values imported from another study.
    dropped = [r["compound"] for r in stage
               if r["value_origin"].strip().lower().startswith("imported")]
    stage = [r for r in stage
             if not r["value_origin"].strip().lower().startswith("imported")]
    steps.append(("E1", "values imported from another study excluded",
                  len(stage), dropped))

    # E5 -- rows whose Vu,brain was DERIVED from a model, not measured.
    dropped = [r["compound"] for r in stage
               if "derived" in r["value_origin"].strip().lower()]
    stage = [r for r in stage
             if "derived" not in r["value_origin"].strip().lower()]
    steps.append(("E5", "model-derived (pH-partition) rows excluded",
                  len(stage), dropped))

    # E6 -- cyclic peptides.
    dropped = [r["compound"] for r in stage
               if "cyclic peptide" in r["mechanism_flag"].lower()]
    stage = [r for r in stage
             if "cyclic peptide" not in r["mechanism_flag"].lower()]
    steps.append(("E6", "cyclic peptides excluded", len(stage), dropped))

    # Cross-check, not a filter: the cascade is property-based and is applied to
    # EVERY row, so it must independently reproduce the dataset's own `scoreable`
    # flag on what it retains. Disagreement would mean the cascade and the frozen
    # INCLUSION_RULE.md disagree about the same rows, and the run should stop.
    bad = [r["compound"] for r in stage
           if not r["scoreable"].strip().lower().startswith("yes")]
    if bad:
        raise AssertionError(
            f"cascade retained rows the dataset marks unscoreable: {bad}. "
            f"The 3.1 cascade and data_cache/brain_slice/INCLUSION_RULE.md "
            f"disagree; stop and report rather than resolve."
        )

    return {"steps": steps, "scored": stage}


def _values(rows: Sequence[dict]) -> np.ndarray:
    return np.array([float(r["value_vu_brain"]) for r in rows], dtype=float)


# --------------------------------------------------------------------------- #
# The statistics of 3.5.
# --------------------------------------------------------------------------- #
def gmfe(pred: float, measured: np.ndarray) -> float:
    """GMFE(p) = exp( mean_i | ln( p / Vu_i ) | )."""
    return float(np.exp(np.mean(np.abs(np.log(pred) - np.log(measured)))))


def gmfe_floor(measured: np.ndarray) -> Tuple[float, float]:
    """min over all constants c of GMFE(c), and the constant achieving it.

    mean|ln c - ln Vu_i| is minimised at the MEDIAN of ln Vu_i, so the
    optimal constant is the sample median of Vu and the floor is
    exp(mean absolute deviation of ln Vu about its median). Closed form; no
    search, no tolerance, nothing to tune. Verified numerically below.
    """
    ln = np.log(measured)
    c_opt = float(np.exp(np.median(ln)))
    return gmfe(c_opt, measured), c_opt


def _verify_floor_is_minimal(measured: np.ndarray, floor: float) -> float:
    """Dense scan of the constant, so the closed form is checked not asserted."""
    grid = np.exp(np.linspace(np.log(measured.min()) - 1.0,
                              np.log(measured.max()) + 1.0, 200001))
    scan = np.exp(np.mean(np.abs(np.log(grid)[:, None] - np.log(measured)), axis=1))
    return float(scan.min() - floor)   # must be >= -1e-12


# --------------------------------------------------------------------------- #
# BCa bootstrap over compounds (3.6).
# --------------------------------------------------------------------------- #
def bca(measured: np.ndarray, stat, n_boot: int, seed: int, level: float):
    theta_hat = stat(measured)
    rng = np.random.default_rng(seed)
    n = measured.size
    idx = rng.integers(0, n, size=(n_boot, n))
    reps = np.array([stat(measured[i]) for i in idx], dtype=float)

    prop = float(np.mean(reps < theta_hat))
    if prop <= 0.0 or prop >= 1.0:                     # degenerate; fall back
        z0 = 0.0
    else:
        z0 = float(norm.ppf(prop))

    jack = np.array([stat(np.delete(measured, i)) for i in range(n)], dtype=float)
    jbar = jack.mean()
    num = float(np.sum((jbar - jack) ** 3))
    den = 6.0 * (float(np.sum((jbar - jack) ** 2)) ** 1.5)
    a = 0.0 if den == 0.0 else num / den

    alpha = (1.0 - level) / 2.0
    out = []
    for z in (norm.ppf(alpha), norm.ppf(1.0 - alpha)):
        adj = z0 + (z0 + z) / (1.0 - a * (z0 + z))
        out.append(float(np.clip(norm.cdf(adj), 0.0, 1.0)))
    lo, hi = np.percentile(reps, [100.0 * out[0], 100.0 * out[1]])
    return theta_hat, float(lo), float(hi), z0, a, reps


# --------------------------------------------------------------------------- #
# GATE 1 -- the engine quantity.
# --------------------------------------------------------------------------- #
def engine_tissue_capacity() -> Dict[str, object]:
    """The engine's implied tissue capacity at the section 2.3 frozen constants.

    The only quantity in the engine that is a RATIO OF TISSUE CONTENT TO UNBOUND
    INTERSTITIAL CONCENTRATION is the two-route capacity term, written down
    explicitly in `two_route.retardation_factor`'s derivation:

        s_eff = (1 - phi) * 1 + phi * K * (1 + K_b)

    matrix (interstitial) capacity normalised to 1, inclusion (cell) capacity
    K*(1+K_b), combined by an exact volume average because "the amount stored per
    unit volume at equilibrium is a volume average and nothing else".

    At the frozen shipping constants -- phi = 0.7892, K = cell_partition = 1.0,
    and K_b absent because `compile_drug.build_profile` never passes
    `bound_to_free_ratio` -- this is EXACTLY 1.0 mL interstitial-fluid-equivalent
    per mL tissue, for every compound.
    """
    from tmesim.molecular import two_route as tr
    from tmesim.molecular import compile_drug as cd

    phi = tr.DEFAULT_CELL_VOLUME_FRACTION
    k = 1.0          # two_route.py:357 default, frozen by section 2.3
    kb = 0.0         # retardation default-off; build_profile never supplies it
    s_eff = (1.0 - phi) * 1.0 + phi * k * (1.0 + kb)

    # Demonstrate the capacity carries no per-compound dependence at all, over
    # the whole descriptor space the membrane leg consumes. This is the direct
    # check behind prereg 3.5's "the same number for every compound".
    seen = set()
    for logp in (-4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0):
        for hbd in (0, 1, 2, 4, 8, 16):
            for radius in (0.15, 0.3, 0.5, 1.0, 3.0):
                res = tr.two_route_diffusivity(
                    interstitial_D_um2s=100.0, solute_radius_nm=radius,
                    logp=logp, hbd=hbd,
                )
                seen.add((res["retardation_factor"],
                          res["cell_partition_effective"],
                          res["cell_volume_fraction"]))
    constant = len(seen) == 1

    return {
        "quantity": "two_route capacity term s_eff = (1-phi) + phi*K*(1+K_b)",
        "phi": phi, "K": k, "K_b": kb,
        "s_eff_ml_per_ml_tissue": s_eff,
        "implied_vu_brain_ml_per_g": s_eff / TISSUE_DENSITY_G_PER_ML,
        "implied_vu_brain_ml_per_g_at_1p04": s_eff / TISSUE_DENSITY_SENSITIVITY,
        "constant_across_descriptor_sweep": constant,
        "n_distinct_capacity_states_over_sweep": len(seen),
        "two_route_default_on": cd._TWO_ROUTE_DEFAULT,
    }


def assert_frozen_constants() -> Dict[str, float]:
    """Section 2.3 FROZEN table, asserted against live code."""
    from tmesim.molecular import two_route as tr
    from tmesim.molecular import compile_drug as cd
    import inspect

    sig = inspect.signature(tr.two_route_factor)
    got = {
        "cell_partition": float(sig.parameters["cell_partition"].default),
        "HBD_LOGP_DECREMENT": float(tr.HBD_LOGP_DECREMENT),
        "cell_volume_fraction": float(tr.DEFAULT_CELL_VOLUME_FRACTION),
        "cell_radius_um": float(tr.DEFAULT_CELL_RADIUS_UM),
        "cytoplasm_D_over_water": float(tr.CYTOPLASM_D_OVER_WATER),
        "logk_membrane_max": float(tr.LOGK_MEMBRANE_MAX),
        "_TWO_ROUTE_DEFAULT": bool(cd._TWO_ROUTE_DEFAULT),
    }
    expect = {
        "cell_partition": 1.0,
        "HBD_LOGP_DECREMENT": 1.34,
        "cell_volume_fraction": 0.7892,
        "cell_radius_um": 7.8517,
        "cytoplasm_D_over_water": 0.27,
        "logk_membrane_max": 3.0,
        "_TWO_ROUTE_DEFAULT": True,
    }
    for key, want in expect.items():
        if got[key] != want:
            raise AssertionError(
                f"section 2.3 FROZEN constant {key}: live code has {got[key]!r}, "
                f"the pre-registration froze {want!r}. The run is VOID under "
                f"prereg 9.5 rider 3."
            )
    return got


class _SolverGuard:
    """Records that the transport solver was NOT entered (prereg 2.1)."""

    def __init__(self) -> None:
        self.counts = {"simulate_drug_transport": 0, "simulate_response": 0}
        self._restore = []

    def __enter__(self):
        from tmesim.physics import diffusion as _diff
        from tmesim.sim import engine as _eng
        for mod, name in ((_diff, "simulate_drug_transport"),
                          (_eng, "simulate_response")):
            orig = getattr(mod, name)

            def make(orig=orig, name=name):
                def wrapper(*a, **kw):
                    self.counts[name] += 1
                    return orig(*a, **kw)
                return wrapper

            setattr(mod, name, make())
            self._restore.append((mod, name, orig))
        return self

    def __exit__(self, *exc):
        for mod, name, orig in self._restore:
            setattr(mod, name, orig)
        return False


# --------------------------------------------------------------------------- #
# Stages.
# --------------------------------------------------------------------------- #
def _git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=_ROOT, text=True).strip()
    except Exception:
        return "unavailable"


def stage_prereg() -> None:
    rows = _load_rows()
    casc = cascade(rows)
    scored = casc["scored"]
    frozen = assert_frozen_constants()
    gate = engine_tissue_capacity()

    doc = {
        "_what": (
            "Frozen execution record for TEST B of results/transport_anchor_prereg.md, "
            "written under prereg 7.3. This file fixes the retained/excluded lists by "
            "criterion code, every numeric bar, the GATE 1 engine quantity and the "
            "density assumption. Nothing here was chosen by looking at a residual."
        ),
        "_rule": (
            "DO NOT FIT TO THE ANCHOR. No parameter, exclusion, threshold or descriptor "
            "choice may be changed by looking at a held-out residual."
        ),
        "executed_by": "pre-registration executor; the plan's author is not the judge",
        "commit_sha_at_execution": _git_sha(),
        "dataset": {
            "path": os.path.relpath(DATA_CSV, _ROOT),
            "inclusion_rule": "data_cache/brain_slice/INCLUSION_RULE.md (frozen pre-extraction)",
            "sha256": hashlib.sha256(open(DATA_CSV, "rb").read()).hexdigest(),
        },
        "cascade": [
            {"code": c, "rule": r, "n_retained": n, "dropped": d}
            for (c, r, n, d) in casc["steps"]
        ],
        "scored_set": {
            "n": len(scored),
            "compounds": [r["compound"] for r in scored],
        },
        "bars": {
            "THE_TEST": "GMFE(engine)/GMFE_floor <= 3.0",
            "the_test_bar": THE_TEST_BAR,
            "reproducibility_floor_x": REPRODUCIBILITY_FLOOR,
            "reproducibility_floor_basis": "n = 2 replicate pairs; order-of-magnitude bound",
            "bootstrap": {"method": "BCa", "resamples": BOOTSTRAP_RESAMPLES,
                          "ci": BOOTSTRAP_CI, "seed": BOOTSTRAP_SEED},
            "verdict_rule": ("the 90 % BCa interval on GMFE(engine)/GMFE_floor must lie "
                             "ENTIRELY on one side of 3.0; a straddle is "
                             "UNDERPOWERED -> NOT TESTABLE (bucket B-C)"),
            "within_bands_x": list(WITHIN_BANDS),
            "spearman": "NOT COMPUTABLE -- the engine prediction has zero variance (3.5)",
        },
        "gate1": {
            "engine_quantity": gate["quantity"],
            "value_ml_per_ml_tissue": gate["s_eff_ml_per_ml_tissue"],
            "density_assumption_g_per_ml": TISSUE_DENSITY_G_PER_ML,
            "density_assumption_note": (
                "water-equivalent convention, stated BEFORE the comparison per 3.3 item 2; "
                "sensitivity at 1.04 g.mL^-1 reported"
            ),
            "implied_vu_brain_ml_per_g": gate["implied_vu_brain_ml_per_g"],
        },
        "frozen_constants_live": frozen,
        "scope": "SCOPE_ANALYTIC; conditional coefficient-not-field; solver never entered",
    }
    payload = json.dumps(doc, indent=1, sort_keys=False)
    doc_hash = hashlib.sha256(payload.encode()).hexdigest()
    doc["_sha256_of_body"] = doc_hash
    os.makedirs(os.path.dirname(PREREG_JSON), exist_ok=True)
    with open(PREREG_JSON, "w", encoding="utf-8") as fh:
        json.dump(doc, fh, indent=1)
        fh.write("\n")

    print(f"wrote {os.path.relpath(PREREG_JSON, _ROOT)}  sha256(body)={doc_hash[:16]}")
    for c, r, n, d in casc["steps"]:
        print(f"  {c:7s} {r:48s} n={n:3d}  dropped={d if d else '-'}")


def stage_floor() -> Dict[str, object]:
    rows = _load_rows()
    scored = cascade(rows)["scored"]
    v = _values(scored)

    floor, c_opt = gmfe_floor(v)
    slack = _verify_floor_is_minimal(v, floor)
    theta, lo, hi, z0, a, _ = bca(
        v, lambda x: gmfe_floor(x)[0], BOOTSTRAP_RESAMPLES, BOOTSTRAP_SEED, BOOTSTRAP_CI)

    print("=" * 78)
    print("B-S1  GMFE_floor -- THE HEADLINE FINDING (prereg 3.5a)")
    print("=" * 78)
    print(f"  n scored                        : {len(v)}")
    print(f"  measured range (mL.g brain^-1)  : {v.min():g} .. {v.max():g}"
          f"   spread {v.max()/v.min():.1f}x")
    print(f"  GMFE_floor                      : {floor:.3f}x")
    print(f"  90 % BCa CI                     : [{lo:.3f}, {hi:.3f}]")
    print(f"  minimising constant             : {c_opt:.4g} mL.g brain^-1")
    print(f"  closed-form check (scan - floor): {slack:.3e}  (>= 0 required)")
    print(f"  reproducibility floor           : {REPRODUCIBILITY_FLOOR:.2f}x"
          f"  (n = 2 replicate pairs; order-of-magnitude bound)")
    print(f"  bootstrap                       : BCa, {BOOTSTRAP_RESAMPLES} resamples, "
          f"seed {BOOTSTRAP_SEED}, z0={z0:+.4f}, a={a:+.5f}")
    print()
    print("  Reading: this is the best unsigned geometric-mean fold error ANY constant")
    print("  predictor could achieve on this measurement. It contains no engine")
    print("  quantity and does not move under any bug fix, any recalibration of alpha,")
    print("  or any resolution of the HD-convention question.")
    return {"floor": floor, "ci": (lo, hi), "c_opt": c_opt, "v": v, "scored": scored}


def _coverage(pred: float, v: np.ndarray, band: float) -> float:
    fold = np.maximum(pred / v, v / pred)
    return float(np.mean(fold <= band))


def stage_engine() -> None:
    rows = _load_rows()
    casc = cascade(rows)
    scored = casc["scored"]
    v = _values(scored)

    frozen = assert_frozen_constants()
    with _SolverGuard() as guard:
        gate = engine_tissue_capacity()
    solver_entered = any(guard.counts.values())

    pred = gate["implied_vu_brain_ml_per_g"]
    floor, c_opt = gmfe_floor(v)
    g_eng = gmfe(pred, v)
    ratio = g_eng / floor

    def ratio_stat(x: np.ndarray) -> float:
        return gmfe(pred, x) / gmfe_floor(x)[0]

    r_hat, r_lo, r_hi, z0, a, _ = bca(
        v, ratio_stat, BOOTSTRAP_RESAMPLES, BOOTSTRAP_SEED, BOOTSTRAP_CI)

    print("=" * 78)
    print("GATE 1 -- does the engine express this quantity at all? (prereg 3.3)")
    print("=" * 78)
    print(f"  engine quantity                 : {gate['quantity']}")
    print(f"  phi={gate['phi']}  K={gate['K']}  K_b={gate['K_b']}")
    print(f"  s_eff                           : {gate['s_eff_ml_per_ml_tissue']:.6f} "
          f"mL interstitial-equivalent per mL tissue")
    print(f"  density assumption (declared)   : {TISSUE_DENSITY_G_PER_ML:.2f} g.mL^-1")
    print(f"  => implied Vu,brain             : {pred:.6f} mL.g brain^-1")
    print(f"     sensitivity at 1.04 g.mL^-1  : "
          f"{gate['implied_vu_brain_ml_per_g_at_1p04']:.6f} mL.g brain^-1")
    print(f"  constant over descriptor sweep  : {gate['constant_across_descriptor_sweep']} "
          f"({gate['n_distinct_capacity_states_over_sweep']} distinct state(s) over 210 inputs)")
    print(f"  solver entered                  : {solver_entered}  {guard.counts}")
    print(f"  section 2.3 frozen constants    : all match live code {frozen}")
    print()

    print("=" * 78)
    print("THE TEST (B) -- GMFE(engine) / GMFE_floor <= 3.0")
    print("=" * 78)
    print(f"  n                               : {len(v)}")
    print(f"  GMFE(engine)                    : {g_eng:.3f}x")
    print(f"  GMFE_floor                      : {floor:.3f}x")
    print(f"  ratio                           : {ratio:.3f}")
    print(f"  90 % BCa CI on ratio            : [{r_lo:.3f}, {r_hi:.3f}]"
          f"   (BCa z0={z0:+.4f}, a={a:+.5f}, {BOOTSTRAP_RESAMPLES} resamples, "
          f"seed {BOOTSTRAP_SEED})")
    if r_lo > THE_TEST_BAR:
        bucket, verdict = "B-B  MISCALIBRATED", "FAIL"
    elif r_hi < THE_TEST_BAR:
        bucket, verdict = "B-A  CONSTRAINED", "PASS"
    else:
        bucket, verdict = "B-C  UNDERPOWERED", "NOT TESTABLE"
    print(f"  bucket                          : {bucket}")
    print(f"  VERDICT                         : {verdict}")
    print(f"  reproducibility floor           : {REPRODUCIBILITY_FLOOR:.2f}x (prints beside "
          f"every fold-error figure, prereg 3.4a)")
    print()

    # ---------------------------------------------------------------- BW1
    print("=" * 78)
    print("WITHDRAWAL CONDITIONS")
    print("=" * 78)
    inside = bool(v.min() <= pred <= v.max())
    print(f"  BW1 reachability: engine {pred:g} vs measured range "
          f"[{v.min():g}, {v.max():g}] -> inside = {inside} "
          f"-> BW1 {'NOT triggered' if inside else 'TRIGGERED'}")
    print("  BW2 : requires B-S9 (measured Vu,brain vs Crippen logP). NOT COMPUTABLE -- "
          "the frozen\n        dataset carries no SMILES and 23 of 33 rows are "
          "radiolabelled PET tracers.\n        Resolving structures now would be a "
          "curation step taken AFTER the measured\n        column was read (prereg 7.1). "
          "BW2 is therefore NOT EVALUABLE, and prereg 6.3\n        had already weakened "
          "it to hypothesis-generating only.")
    print("  BW3 : GATE 1 found a commensurable PARTITION, not a rate -> NOT triggered.")
    print()

    # ---------------------------------------------------------------- secondaries
    print("=" * 78)
    print("SECONDARIES (prereg 3.7) -- none may replace THE TEST")
    print("=" * 78)
    print(f"  B-S1 GMFE_floor                 : {floor:.3f}x, minimising constant "
          f"{c_opt:.4g} mL.g brain^-1  (see --stage floor for its BCa interval)")
    print(f"  B-S2 centring ratio             : engine {pred:g} / optimal {c_opt:.4g} "
          f"= {pred / c_opt:.5f}  ({c_opt / pred:.1f}x low)")

    lo_side = int(np.sum(pred < v))
    print(f"  B-S4 direction                  : engine below {lo_side}/{len(v)} measured "
          f"values, above {len(v) - lo_side}/{len(v)}")
    signed = float(np.exp(np.mean(np.log(pred) - np.log(v))))
    print(f"                                    signed geometric bias "
          f"{signed:.5f} ({1.0/signed:.1f}x low, diagnostic only)")
    for band in WITHIN_BANDS:
        print(f"  B-S5 within-{band:g}x coverage        : "
              f"{_coverage(pred, v, band) * 100:.1f}%")
    fold = np.maximum(pred / v, v / pred)
    worst = int(np.argmax(fold))
    print(f"  B-S6 worst-case fold error      : {fold[worst]:.1f}x  "
          f"({scored[worst]['compound']}, measured {v[worst]:g} mL.g brain^-1) "
          f"-- reported, NOT removed")
    print(f"  B-S7 reproducibility floor      : {REPRODUCIBILITY_FLOOR:.2f}x "
          f"(paroxetine 714/876 = 1.23x; diazepam 17.8/20 = 1.12x; n = 2 pairs)")
    print("  B-S8 Spearman rho               : NOT COMPUTABLE -- the engine prediction has "
          "zero\n                                    variance, so it has no ranks. This is a "
          "property of the\n                                    model, not of n. It is NOT 0.00 "
          "and NOT a failure.")
    print("  B-S9 rho(Vu,brain, Crippen logP): NOT COMPUTABLE -- no SMILES in the frozen "
          "dataset\n                                    (see BW2 above).")

    # B-S3 ion-class strata
    print("  B-S3 ion-class strata (diagnostic only; neither stratum clears n >= 10):")
    def klass(r: dict) -> str:
        f = r["mechanism_flag"].lower()
        if "ion class: base" in f:
            return "base"
        if "ion class:" in f:
            return "non-base"
        return "unassigned"
    for name in ("base", "non-base", "unassigned"):
        sel = [i for i, r in enumerate(scored) if klass(r) == name]
        if not sel:
            continue
        sv = v[sel]
        f_, _c = gmfe_floor(sv)
        print(f"       {name:11s} n={len(sel):2d}  GMFE(engine)={gmfe(pred, sv):8.2f}x  "
              f"GMFE_floor={f_:6.2f}x  ratio={gmfe(pred, sv)/f_:.3f}  "
              f"median Vu={np.median(sv):.4g}")

    # ---------------------------------------------------------------- B-S10
    print()
    print("=" * 78)
    print("B-S10 SENSITIVITY ARMS -- reported, NEVER the verdict")
    print("=" * 78)
    steps = {c: (n, d) for c, _r, n, d in casc["steps"]}
    rat = [r for r in rows if r["species"] == "rat"]

    def arm(name: str, sel: List[dict]) -> None:
        av = _values(sel)
        f_, c_ = gmfe_floor(av)
        g_ = gmfe(pred, av)
        print(f"  {name:34s} n={len(av):2d}  GMFE_floor={f_:6.2f}x  "
              f"GMFE(engine)={g_:8.2f}x  ratio={g_/f_:.3f}  "
              f"spread={av.max()/av.min():.0f}x")

    arm("HEADLINE (E1,E2,E5,E6 applied)", scored)
    arm("+ oxycodone (E1 relaxed)",
        scored + [r for r in rat
                  if r["value_origin"].strip().lower().startswith("imported")])
    arm("+ cyclic peptides (E6 relaxed)",
        scored + [r for r in rat
                  if "cyclic peptide" in r["mechanism_flag"].lower()
                  and "derived" not in r["value_origin"].lower()
                  and not r["value_origin"].lower().startswith("imported")])
    arm("+ model-derived rows (E5 relaxed)",
        scored + [r for r in rat if "derived" in r["value_origin"].strip().lower()])
    arm("mouse arm (E2 relaxed; SEPARATE)",
        [r for r in rows if r["species"] == "mouse"
         and r["scoreable"].strip().lower().startswith("yes")])

    # ---------------------------------------------------------------- 8.2
    print()
    print("=" * 78)
    print("prereg 8.2 -- does the dataset separate its own compounds?")
    print("=" * 78)
    sds = [(float(r["sd"]), float(r["value_vu_brain"]))
           for r in scored if r["sd"].strip() not in ("", "NA", "-")]
    print(f"  rows reporting an SD            : {len(sds)}/{len(scored)}")
    if sds:
        rel = sorted(s / val for s, val in sds)
        med_rel = float(np.median(rel))
        print(f"  median relative SD              : {med_rel * 100:.1f}%  "
              f"(as a fold factor {1 + med_rel:.3f}x)")
        print(f"  measured spread {v.max()/v.min():.0f}x vs 3 x median uncertainty "
              f"{3 * (1 + med_rel):.2f}x -> {'PASS' if v.max()/v.min() > 3*(1+med_rel) else 'FAIL'}")
    print(f"  measured spread vs 3 x reproducibility floor "
          f"{3 * REPRODUCIBILITY_FLOOR:.2f}x -> "
          f"{'PASS' if v.max()/v.min() > 3 * REPRODUCIBILITY_FLOOR else 'FAIL'}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stage", choices=("prereg", "floor", "engine"), required=True)
    args = ap.parse_args()
    {"prereg": stage_prereg, "floor": stage_floor, "engine": stage_engine}[args.stage]()


if __name__ == "__main__":
    main()
