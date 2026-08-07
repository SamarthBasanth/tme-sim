"""Replace the blood-BRAIN-barrier permeability with a MEASURED breast one.

THE DEFECT THIS TESTS
---------------------
``PS = vascular_density * bbb_penetrance``, and ``bbb_penetrance`` is a
structure-derived blood-BRAIN-barrier term. Breast has no such barrier. This
is the second of three defects diagnosed in
``results/transvascular_exchange_revalidation.md``; the first (a one-way
delivery law) and third (a dimensionless partition used as a rate) are fixed.

WHAT THE MEASURED DATA SAYS, AND WHY THE FIX IS DRUG-INDEPENDENT
-----------------------------------------------------------------
In breast tumour the measured 5-FU:methotrexate permeability ratio is
**1.02x** (PMID 9205062 / 9679959) — essentially drug-independent — against
28.1x from ``bbb_penetrance``. The measured ratio rises monotonically with
barrier tightness: ~1.0x in breast tumour, ~2x in skeletal muscle, 6-12x in
normal brain. So for breast the honest replacement is a single TISSUE-level
permeability, not a per-drug one.

Measured human breast tumour PS: **0.06-0.10 min^-1** (DCE-MRI
distributed-parameter fits; 0.22 as an upper bound from n=18 dogs with
spontaneous mammary tumours, down-weighted because its extraction fraction
E = 0.56 breaches the E > 0.2 regime where K_trans underestimates PS).

HOW IT IS APPLIED
-----------------
Spatial structure is preserved and only the LEVEL is replaced: set a
drug-independent ``bbb_penetrance = PS_target / mean(vascular_density)`` over
the tissue mask, so every drug's tissue-mean PS equals the measured value
while vessels stay where the vessels are. No engine change is needed.

PREDICTION, STATED BEFORE RUNNING SO IT IS FALSIFIABLE
-------------------------------------------------------
This is expected to make BOTH anchored drugs WORSE, not better. Reason: to
reach the measured Kp of 0.60 the model needs ``PS = 1.5 * (clearance +
binding)``, which is 8.4e-3 /s for 5-FU and 4.0e-2 /s for methotrexate —
**5x to 40x above the measured breast PS of 1.0-1.7e-3 /s**.

If that holds, the interesting consequence is NOT "the fix failed". It is
that with a MEASURED permeability and a MEASURED anchor the system becomes
OVER-determined and inconsistent, which localises the remaining error to the
third parameter: the alpha-scaled binding sink, fit once on GBM doxorubicin
penetration and reused verbatim on breast. That would be a much sharper
result than "underdetermined", and it is checkable.
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import dataclasses  # noqa: E402

from tmesim.core import PDConfig, PDEConfig, PKConfig, SimConfig  # noqa: E402
from tmesim.molecular.binding_site_barrier import (  # noqa: E402
    apply_barrier_to_drugparams,
)
from tmesim.molecular.breast_library import (  # noqa: E402
    BREAST_DRUGS,
    get_molecule_breast,
)
from tmesim.molecular.compile_drug import molecule_to_drugparams  # noqa: E402
from tmesim.molecular.library import get_molecule  # noqa: E402
from tmesim.reconstruct.real_pipeline import (  # noqa: E402
    apply_drug_sensitivity,
    reconstruct_real,
)
from tmesim.sim.engine import simulate_response  # noqa: E402
from tmesim.validation import delivery_metrics  # noqa: E402
from tmesim.validation.calibration import calibrated_drugparams  # noqa: E402

ALPHA = 148.60815765445545
SECTION = "data_cache/xenium_breast_janesick/rep1"
T = 3600.0 * 6

#: Measured human breast tumour PS, min^-1 -> 1/s. See module docstring.
PS_TARGETS_PER_MIN = {"low (0.06)": 0.06, "mid (0.10)": 0.10, "upper (0.22)": 0.22}

#: drug -> (measured interstitial:plasma AUC ratio, SE, PMID)
ANCHOR = {"5-fluorouracil": (0.61, 0.11, "9205062"),
          "methotrexate": (0.60, 0.20, "9679959")}


def run(micro, name, bbb_override, s0, st0, m0):
    micro.cells.sensitivity = np.array(s0, copy=True)
    micro.cells.state = np.array(st0, copy=True)
    micro.metadata = dict(m0)
    apply_drug_sensitivity(micro, name, mode="uniform")
    mol = get_molecule_breast(name) if name in BREAST_DRUGS else get_molecule(name)
    dp = apply_barrier_to_drugparams(
        calibrated_drugparams(molecule_to_drugparams(mol), ALPHA), mol)
    if bbb_override is not None:
        dp = dataclasses.replace(dp, bbb_penetrance=float(bbb_override))
    cfg = SimConfig(drug=dp, pk=PKConfig(),
                    pde=PDEConfig(dt=10.0, total_time=T, convection=True,
                                  transvascular_exchange=True),
                    pd=PDConfig(), seed=0)
    res = simulate_response(micro, cfg)
    summ = delivery_metrics.tissue_auc_summary(micro, res)
    pl = float(res.extras["plasma_mean"]) * T
    return {"ratio": summ["mean_tumor_auc"] / pl,
            "bbb": float(dp.bbb_penetrance),
            "k_loss": float(dp.clearance + dp.binding)}


def main():
    print("reconstructing real breast section ...", flush=True)
    micro = reconstruct_real(source="xenium", sample=SECTION, dx_um=12.0,
                             oxygen_mode="pde", state_mode="oxygen")
    tis = micro.domain.tissue_mask
    mean_vasc = float(micro.fields.vascular_density[tis].mean())
    print(f"  tissue-mean vascular density = {mean_vasc:.5f}", flush=True)
    s0 = np.array(micro.cells.sensitivity, copy=True)
    st0 = np.array(micro.cells.state, copy=True)
    m0 = dict(micro.metadata)

    out = {"mean_vascular_density": mean_vasc, "arms": {}}

    print("\nbaseline: bbb_penetrance as shipped (blood-BRAIN-barrier)", flush=True)
    for d in ANCHOR:
        r = run(micro, d, None, s0, st0, m0)
        out["arms"].setdefault("baseline", {})[d] = r
        print(f"  {d:16s} bbb={r['bbb']:.4f} ratio={r['ratio']:.4f}", flush=True)

    for label, per_min in PS_TARGETS_PER_MIN.items():
        ps = per_min / 60.0                       # min^-1 -> s^-1
        bbb_new = ps / mean_vasc                  # drug-INDEPENDENT
        print(f"\nmeasured breast PS {label} = {ps:.4e}/s "
              f"-> bbb_penetrance {bbb_new:.4f} (same for every drug)", flush=True)
        for d in ANCHOR:
            r = run(micro, d, bbb_new, s0, st0, m0)
            out["arms"].setdefault(label, {})[d] = r
            print(f"  {d:16s} ratio={r['ratio']:.4f}", flush=True)

    print("\n" + "=" * 74)
    print("PREDICTED vs MEASURED interstitial:plasma AUC ratio, real breast")
    print("=" * 74)
    hdr = f"{'arm':16s}" + "".join(f"{d[:14]:>16s}" for d in ANCHOR)
    print(hdr)
    for arm, rows in out["arms"].items():
        print(f"{arm:16s}" + "".join(
            f"{rows[d]['ratio']:9.4f} ({rows[d]['ratio'] / ANCHOR[d][0]:4.2f}x)"
            for d in ANCHOR))
    print(f"{'MEASURED':16s}" + "".join(
        f"{ANCHOR[d][0]:9.2f} (PMID {ANCHOR[d][2]})" for d in ANCHOR))

    print("\nTHE DIAGNOSTIC -- what PS would each drug NEED to hit its anchor?")
    print("  (steady state Kp = PS/(PS + clearance + binding) -> PS = Kp*k/(1-Kp))")
    lo, hi = min(PS_TARGETS_PER_MIN.values()) / 60.0, max(PS_TARGETS_PER_MIN.values()) / 60.0
    for d in ANCHOR:
        k = out["arms"]["baseline"][d]["k_loss"]
        kp = ANCHOR[d][0]
        need = kp * k / (1.0 - kp)
        print(f"  {d:16s} k_loss={k:.4e}/s  needs PS={need:.4e}/s  "
              f"= {need / hi:5.1f}x the TOP of the measured breast range")
    print(f"\n  measured breast PS range: {lo:.4e} .. {hi:.4e} /s")
    print("  If 'needs PS' >> measured PS for BOTH drugs, permeability is NOT the\n"
          "  remaining free parameter -- the alpha-scaled binding sink is.")

    p = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                      "..", "results",
                                      "breast_tissue_permeability_data.json"))
    with open(p, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\nwrote {p}")


if __name__ == "__main__":
    main()
