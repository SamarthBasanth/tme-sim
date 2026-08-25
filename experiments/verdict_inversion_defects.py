"""Plant each of ``verdict_inversion_census``'s four checks' defects and watch it bite.

A check that has never been watched to fail is a claim, not a check. This file
plants one defect per check, runs the census on a two-compound subset at the
IDENTICAL shipped configuration (nothing is scaled, nothing shipped is changed —
every plant is a monkeypatch on ``experiments/verdict_inversion_census`` itself,
never on ``tmesim/``), and prints the verbatim failure.

    python experiments/verdict_inversion_defects.py                # all four
    python experiments/verdict_inversion_defects.py --defect inert

The subset is ``bevacizumab`` (the zero-source degenerate case, needed by check
D) plus one live-source compound. Defect 2 is the only one that needs the
``reproduce`` stage and pays for its own four solves.

**Every plant is expected to FAIL a check.** A plant that leaves every check
green is reported as ``BLIND`` and is the interesting outcome, not the boring
one.
"""
from __future__ import annotations

import argparse
import os
import sys
import time
import traceback
from typing import Dict, List, Optional

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
for _p in (_ROOT, _HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import verdict_inversion_census as V  # noqa: E402

SUBSET = ("bevacizumab", "doxorubicin")


def _pristine():
    return {"_exchange_arm": V._exchange_arm,
            "_churn_mask": V._churn_mask,
            "_is_degenerate": V._is_degenerate,
            "compile_reproduce": None}


def _restore(saved) -> None:
    V._exchange_arm = saved["_exchange_arm"]
    V._churn_mask = saved["_churn_mask"]
    V._is_degenerate = saved["_is_degenerate"]


# --------------------------------------------------------------------------- #
# The four plants
# --------------------------------------------------------------------------- #
def plant_inert() -> str:
    """DEFECT 1 -> check A. The exchange flag is dropped on its way from
    PDEConfig to the solver, so the ON arm IS the OFF arm. This is the real
    defect ``tmesim/validation/propagation.py`` records happening (auc_sum
    bit-identical between the two arms), and under it this file would report a
    0 % inversion rate with a completely green run."""
    V._exchange_arm = lambda exchange: False
    return "A"


def plant_wrong_arm_churn() -> str:
    """DEFECT 2 -> check C. The per-cell churn is diffed against the arm it came
    from instead of against the other arm. Reports 0 changed cells for every
    compound, which reads as 'the split moves nothing'."""
    V._churn_mask = lambda lab_off, lab_on: lab_off != lab_off
    return "C"


def plant_degenerate_missed() -> str:
    """DEFECT 3 -> check D. The zero-source test is written as ``< 0`` instead
    of ``== 0``, so bevacizumab is never detected as degenerate. Its DELIVERY
    verdict — which is an INPUT — then enters the denominator as a compound with
    a perfectly stable verdict and deflates the inversion rate."""
    V._is_degenerate = lambda dp: float(dp.bbb_penetrance) < 0.0
    return "D"


def plant_wrong_compile_path() -> str:
    """DEFECT 4 -> check B. The reproduction gate compiles §4.2's compounds
    through the PRODUCT's path (per-compound IC50) instead of §4.2's own
    (shared IC50). Plausible — it is the path the rest of this file uses — and
    it makes the gate compare a different measurement to the published one."""
    import isomer_penetration_map as IPM

    _real = IPM.compile_all
    IPM._ORIG_COMPILE_ALL = _real

    def _fake(alpha, two_route=True):
        out = _real(alpha, two_route=two_route)
        for name in out:
            out[name] = dict(out[name])
            out[name]["drug_params"] = V.compile_shipped(name, alpha)
        return out

    IPM.compile_all = _fake
    return "B"


PLANTS = {
    "inert": (plant_inert, "census"),
    "churn": (plant_wrong_arm_churn, "census"),
    "degenerate": (plant_degenerate_missed, "census"),
    "compile": (plant_wrong_compile_path, "reproduce"),
}


# --------------------------------------------------------------------------- #
def run_one(key: str) -> Dict:
    plant, stage = PLANTS[key]
    saved = _pristine()
    print("\n" + "=" * 82)
    print(f"PLANT '{key}'  ({plant.__doc__.splitlines()[0].strip()})")
    print("=" * 82, flush=True)
    t0 = time.time()
    target = plant()
    rec = {"defect": key, "stage": stage, "targets_check": target}
    try:
        if stage == "reproduce":
            repro = V.stage_reproduce(verbose=False)
            V.check_reproduces_published(repro)
        else:
            census = V.stage_census(verbose=False, names=SUBSET)
            V.run_checks(census, repro=None, verbose=False)
        rec["bit"] = False
        rec["verbatim"] = None
        print(f"  *** BLIND *** plant '{key}' produced NO failure. The check "
              f"aimed at it ({target}) does not cover it.")
    except Exception as exc:  # noqa: BLE001
        rec["bit"] = True
        rec["exception"] = type(exc).__name__
        rec["verbatim"] = str(exc)
        print(f"  BIT: {type(exc).__name__}")
        print("  " + "\n  ".join(str(exc).splitlines()))
    finally:
        _restore(saved)
        if key == "compile":
            import isomer_penetration_map as IPM
            IPM.compile_all = IPM._ORIG_COMPILE_ALL
    rec["wall_s"] = time.time() - t0
    print(f"  ({rec['wall_s']:.0f}s)", flush=True)
    return rec


def main(argv: Optional[List[str]] = None) -> List[Dict]:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--defect", action="append", default=None,
                    choices=sorted(PLANTS))
    args = ap.parse_args(argv)
    keys = args.defect or list(PLANTS)

    print(f"config: {V._tissue_cfg()}")
    print(f"subset: {SUBSET}")
    V._assert_frozen("defects/top")

    out = [run_one(k) for k in keys]

    V._assert_frozen("defects/bottom")
    print("\n" + "=" * 82)
    print("SUMMARY")
    print("=" * 82)
    for r in out:
        print(f"  {r['defect']:<12} -> check {r['targets_check']}  "
              f"{'BIT' if r['bit'] else '*** BLIND ***'}  ({r['wall_s']:.0f}s)")
    return out


if __name__ == "__main__":
    main()
