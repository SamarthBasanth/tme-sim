"""Breast transvascular-delivery anchor test: predicted vs MEASURED
interstitial:plasma AUC ratio on the real Xenium breast section.

This is the breast analogue of what ``tmesim.validation.calibration`` /
``results/pk_reality.md`` do for GBM, but against a different measured
quantity. GBM's external anchor is a PENETRATION DEPTH (doxorubicin's
half-distance from vessels, PMID 16361566). No breast-specific equivalent of
that exists -- see ``results/breast_transport_anchors.md`` section 3 for the
honest null, and for why the published breast "penetration" numbers are not
mutually convertible.

What DOES exist for breast is better suited to this engine anyway: two
independent human cohorts measured drug concentration directly in primary
breast tumour interstitium by in vivo microdialysis.

    5-fluorouracil   interstitial AUC = 61 +/- 11 % of plasma AUC   PMID 9205062 (n=10)
    methotrexate     interstitial AUC = 60 +/- 20 % of plasma AUC   PMID 9679959 (n=9)

Two chemically unrelated molecules -- a 130 Da hydrophilic antimetabolite and
a 454 Da polar antifolate -- landing on the same number.

THE ENGINE QUANTITY THAT CORRESPONDS TO THAT MEASUREMENT

    ratio = mean_tumor_auc [uM*s] / (plasma_mean [uM] * total_time [s])

i.e. free interstitial AUC sampled at tumour cells, over plasma AUC across the
same window. The transport PDE's state variable is the free interstitial
concentration, which is what a microdialysis probe samples, so these are the
same quantity and not merely analogous.

IC50 is irrelevant here and is deliberately not controlled: AUC is a pure
transport quantity and never enters the Hill kill term.

Nothing under ``tmesim/`` is modified by this script. It only calls the
shipping pipeline. See ``results/breast_kp_validation.md`` for the result and
what it diagnoses.
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

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

#: The alpha fit ONCE on GBM doxorubicin penetration depth and reused verbatim
#: everywhere, breast included (see results/kill_calibration_breast.md
#: limitation 9). Deliberately NOT refit here -- refitting it to this anchor
#: would be exactly the kind of rigging this project forbids.
ALPHA = 148.60815765445545
SECTION_PATH = "data_cache/xenium_breast_janesick/rep1"
DX_UM = 12.0
TOTAL_TIME_S = 3600.0 * 6

#: drug -> (measured interstitial:plasma AUC ratio, SE, n, PMID)
ANCHOR = {
    "5-fluorouracil": (0.61, 0.11, 10, "9205062"),
    "methotrexate": (0.60, 0.20, 9, "9679959"),
}

#: The full breast panel plus the cross-tissue drugs the breast reports
#: already reuse, so the SPREAD is visible and not just the two anchored rows.
DRUGS = sorted(BREAST_DRUGS) + ["methotrexate", "doxorubicin", "cisplatin"]


def compile_drug(name):
    """Same compilation path experiments/kill_calibration_breast.py uses."""
    mol = get_molecule_breast(name) if name in BREAST_DRUGS else get_molecule(name)
    dp = calibrated_drugparams(molecule_to_drugparams(mol), ALPHA)
    return apply_barrier_to_drugparams(dp, mol), mol


def main():
    print("reconstructing real breast section ...", flush=True)
    micro = reconstruct_real(source="xenium", sample=SECTION_PATH, dx_um=DX_UM,
                             oxygen_mode="pde", state_mode="oxygen")
    print(f"  {int(micro.cells.n)} cells", flush=True)
    sens0 = np.array(micro.cells.sensitivity, copy=True)
    state0 = np.array(micro.cells.state, copy=True)
    meta0 = dict(micro.metadata)

    rows = []
    for name in DRUGS:
        # Reset the shared reconstruction so runs cannot contaminate each other.
        micro.cells.sensitivity = np.array(sens0, copy=True)
        micro.cells.state = np.array(state0, copy=True)
        micro.metadata = dict(meta0)
        apply_drug_sensitivity(micro, name, mode="uniform")

        dp, mol = compile_drug(name)
        cfg = SimConfig(
            drug=dp, pk=PKConfig(),
            pde=PDEConfig(dt=10.0, total_time=TOTAL_TIME_S, convection=True),
            pd=PDConfig(), seed=0,
        )
        res = simulate_response(micro, cfg)
        summ = delivery_metrics.tissue_auc_summary(micro, res)
        plasma_auc = float(res.extras["plasma_mean"]) * TOTAL_TIME_S
        rows.append({
            "drug": name, "mw": float(mol.get("mw")),
            "bbb_penetrance": float(dp.bbb_penetrance),
            "diffusivity": float(dp.diffusivity), "binding": float(dp.binding),
            "mean_tumor_auc": summ["mean_tumor_auc"],
            "plasma_auc_uM_s": plasma_auc,
            "ratio": summ["mean_tumor_auc"] / plasma_auc,
        })
        print(f"  {name:18s} bbb={dp.bbb_penetrance:.4f}  "
              f"ratio={rows[-1]['ratio']:.4f}", flush=True)

    rows.sort(key=lambda r: r["ratio"])

    print("\n" + "=" * 74)
    print("PREDICTED interstitial:plasma AUC ratio, real breast section")
    print("MEASURED anchor: 5-FU 0.61+/-0.11 (n=10), MTX 0.60+/-0.20 (n=9)")
    print("=" * 74)
    print(f"{'drug':18s} {'MW':>7s} {'bbb_pen':>9s} {'ratio':>10s} {'vs 0.60':>9s}")
    for r in rows:
        print(f"{r['drug']:18s} {r['mw']:7.1f} {r['bbb_penetrance']:9.4f} "
              f"{r['ratio']:10.4f} {r['ratio'] / 0.60:8.2f}x")

    lo, hi = rows[0]["ratio"], rows[-1]["ratio"]
    print(f"\nmodel spread across panel : {hi / lo:,.0f}x  ({lo:.4f} .. {hi:.4f})")
    print("measured spread, 2 anchors: 1.02x (0.60 vs 0.61)")

    print("\nVERDICT on the anchored drugs")
    by = {r["drug"]: r for r in rows}
    for drug, (m, se, n, pmid) in ANCHOR.items():
        pred = by[drug]["ratio"]
        lo_b, hi_b = m - se, m + se
        verdict = "INSIDE" if lo_b <= pred <= hi_b else "OUTSIDE"
        print(f"  {drug:18s} predicted {pred:8.4f}  vs {m:.2f}+/-{se:.2f} "
              f"[{lo_b:.2f},{hi_b:.2f}] (n={n}, PMID {pmid})  -> {verdict} "
              f"({pred / m:.2f}x)")
    r5, rm = by["5-fluorouracil"]["ratio"], by["methotrexate"]["ratio"]
    print(f"\n  model says 5-FU/MTX differ by {r5 / rm:.1f}x; "
          f"measurement says {0.61 / 0.60:.2f}x")

    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "..", "results", "breast_kp_anchor_data.json")
    with open(os.path.normpath(out), "w") as fh:
        json.dump(rows, fh, indent=2)
    print(f"\nwrote {os.path.normpath(out)}")


if __name__ == "__main__":
    main()
