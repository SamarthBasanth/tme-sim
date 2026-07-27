"""Capstone: the FULLY INTEGRATED realistic pipeline.

For each drug: structure -> (Stokes-Einstein + Brownian + Ogston-Renkin mesh)
effective diffusivity -> measured GDSC IC50 -> calibrated tissue-uptake loss
(alpha fit ONCE on doxorubicin's 40-50 um half-distance anchor,
converted to the e-folding length L = d_half/ln2 ~= 65 um) -> binding-site barrier (antibodies)
-> predicted tissue penetration length. Shows the panel now reproduces realistic,
literature-consistent penetration instead of the uncalibrated ~700 um overshoot.
"""
import warnings; warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from tmesim.molecular import library
from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
from tmesim.validation.calibration import calibrate_alpha, calibrated_drugparams, penetration_length_um

PANEL = ["temozolomide", "lomustine", "cisplatin", "methotrexate", "doxorubicin", "bevacizumab"]

alpha = calibrate_alpha(drug="doxorubicin")  # e-fold target 45/ln2 ~= 64.9 um
print(f"calibrated alpha (fit once on doxorubicin, 40-50um half-distance "
      f"-> e-fold L~=64.9um) = {alpha:.1f}\n")

rows = []
for name in PANEL:
    mol = library.get_molecule(name)
    dp_cal = calibrated_drugparams(name, alpha=alpha)          # + calibrated tissue loss
    dp_full = apply_barrier_to_drugparams(dp_cal, mol)          # + binding-site barrier
    L_cal = penetration_length_um(dp_cal)
    L_full = penetration_length_um(dp_full)
    rows.append(dict(name=name, mw=mol.get("mw"), D=dp_full.diffusivity, ic50=dp_full.ic50,
                     bbb=dp_full.bbb_penetrance, L_precal_note="", L_cal=L_cal, L_full=L_full))
    print(f"{name:14s} MW={mol.get('mw'):8.0f}  D_eff={dp_full.diffusivity:7.1f}  "
          f"ic50={dp_full.ic50:7.2f}uM  L(calibrated)={L_cal:7.1f}um  L(+barrier)={L_full:7.1f}um")

# ---- figure: realistic penetration length per drug (log scale) ----
fig, ax = plt.subplots(figsize=(9, 5))
names = [r["name"] for r in rows]
L = [r["L_full"] for r in rows]
colors = ["#d62728" if n == "bevacizumab" else "#1f77b4" for n in names]
bars = ax.bar(names, L, color=colors)
ax.axhspan(40, 50, color="green", alpha=0.15, label="doxorubicin anchor ~40-50 um (Primeau 2005)")
ax.set_yscale("log")
ax.set_ylabel("predicted tissue penetration length (um, log)")
ax.set_title("Realistic pipeline: structure -> mesh -> measured IC50 -> calibrated loss -> barrier")
for b, v in zip(bars, L):
    ax.text(b.get_x()+b.get_width()/2, v*1.1, f"{v:.0f}", ha="center", fontsize=9)
ax.legend(); plt.xticks(rotation=30, ha="right"); plt.tight_layout()
fig.savefig("results/realistic_penetration.png", dpi=130)
print("\nwrote results/realistic_penetration.png")

# ---- honest checks ----
byname = {r["name"]: r for r in rows}
assert 40 <= byname["doxorubicin"]["L_full"] <= 50, "doxorubicin should sit at its calibration anchor"
assert byname["temozolomide"]["L_full"] > 3*byname["doxorubicin"]["L_full"], "TMZ should penetrate much further"
assert byname["bevacizumab"]["L_full"] < byname["temozolomide"]["L_full"]/10, "antibody must be barrier-trapped"
print("checks passed: dox at anchor; TMZ penetrates far; antibody trapped by barrier")
