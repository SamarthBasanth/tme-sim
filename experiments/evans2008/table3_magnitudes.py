"""Convert Evans 2008 Table 3 EF5-binding slopes into mmHg/um via the Table 1 ladder.

This is the ASSUMPTION-BEARING arm of the comparison and is reported, never scored
(pre-registration section 3, step 3). Table 3 gives, per patient and per perivascular
oxygenation class, the intercept (%CRB) and slope (%CRB/um x 1000) of the
binding-vs-distance line. Table 1 gives five (%CRB, pO2) class boundaries. A power law
fitted to the four finite boundaries converts a %CRB slope at a given intercept into
mmHg/um by the chain rule.

All Table 1 and Table 3 values below are transcribed from the PMC XML tables printed by
``dump_tables.py`` -- no figure was digitised.
"""
from __future__ import annotations

import numpy as np

# Table 1: the four finite (%CRB, pO2 mmHg) class boundaries. The 100 %CRB -> 0 mmHg
# row is a limit, not a datum, and is excluded from the fit.
T1_CRB = np.array([1.0, 3.0, 10.0, 30.0])
T1_PO2 = np.array([100.0, 20.0, 4.0, 0.75])

# Table 3, bands 1-5: (patient, class, intercept %CRB, slope %CRB/um x 1000).
T3 = [
    ("A/HPCTA", "oxygenated", 0.3, -0.9),
    ("B/OA", "oxygenated", 1.2, -0.7),
    ("C/AO", "oxygenated", 1.0, -5.8),
    ("D/GB", "oxygenated", 1.7, -2.4),
    ("D/GB", "mild", 5.9, -5.3),
    ("D/GB", "moderate", 13.2, -43.6),
    ("E/GB", "oxygenated", 1.2, 3.9),
    ("E/GB", "mild", 5.0, -52.1),
    ("E/GB", "moderate", 18.0, -153.4),
    ("F/r-GB", "oxygenated", 2.4, 5.9),
    ("F/r-GB", "mild", 5.5, 16.5),
    ("F/r-GB", "moderate", 13.9, 71.4),
    ("F/r-GB", "severe", 33.6, 9.9),
]


def main() -> None:
    k, b = np.polyfit(np.log(T1_CRB), np.log(T1_PO2), 1)
    a = float(np.exp(b))
    pred = a * T1_CRB ** k
    print(f"Table 1 power law: pO2 = {a:.4f} * CRB^{k:.4f}")
    print("  relative residuals:", np.round(pred / T1_PO2 - 1.0, 4).tolist())
    print()
    print(f"{'patient':10s} {'class':12s} {'icpt %CRB':>10s} {'pO2@icpt':>9s} "
          f"{'slope %CRB/um':>14s} {'dpO2/dCRB':>10s} {'mmHg/um':>10s}")
    out = []
    for pat, cls, icpt, slope_k in T3:
        crb_slope = slope_k / 1000.0                  # %CRB per um
        po2 = a * icpt ** k                           # mmHg at the intercept
        dpo2_dcrb = a * k * icpt ** (k - 1.0)         # mmHg per %CRB, at the intercept
        mmhg_per_um = dpo2_dcrb * crb_slope
        out.append((pat, cls, mmhg_per_um))
        print(f"{pat:10s} {cls:12s} {icpt:10.1f} {po2:9.3f} "
              f"{crb_slope:14.5f} {dpo2_dcrb:10.3f} {mmhg_per_um:10.5f}")
    mags = np.array([abs(m) for _, _, m in out])
    print()
    print(f"|mmHg/um| over the 13 Table 3 rows: min {mags.min():.5f}  "
          f"median {np.median(mags):.5f}  max {mags.max():.5f}")
    neg = [(p, c, m) for p, c, m in out if m < 0]
    print(f"rows implying pO2 FALLS away from the vessel (mmHg/um < 0): "
          f"{len(neg)}/{len(out)}")


if __name__ == "__main__":
    main()
