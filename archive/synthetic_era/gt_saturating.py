"""Ground-truth motif #7: SATURATING dose-response with tissue binding.

Biological assumption
---------------------
Delivery sets the *local* area-under-curve (AUC) of free drug that the tumour
actually sees; the cytotoxic effect then follows a SATURATING (Emax / Hill)
dose-response in that local exposure -- once enough drug arrives, adding more
buys almost nothing. Tissue binding modulates the plateau: strongly-bound
(large, lipophilic) drugs are sequestered in tissue, so they (a) reach their
half-maximal effect at a LOWER delivered exposure (they "plateau early") but
(b) their attainable ceiling (Emax) is capped lower because only the free
fraction hits the target.

This is deliberately NOT the multiplicative delivery*potency form. Kill is a
concave, saturating function of a single delivered-exposure variable, so two
drugs that both clear the EC50 look nearly identical no matter how much more
potent one is -- delivery, not raw potency, decides most head-to-heads.

Construction (all deterministic, seeded noise)
----------------------------------------------
binding b[d]        : min-max normalised drug_mw in [0,1]  (large drug = strongly bound)
local exposure AUC  : dose * BBB penetration * perfusion(vasc), attenuated by an
                      ECM diffusion barrier that hurts large molecules more:
                         AUC[p,d] = dose * bbb[d] * vasc[p] * exp(-k_ecm * ecm[p] * (0.5+b[d]))
EC50 (exposure)     : from intrinsic potency & patient sensitivity, lowered by binding
                         EC50[p,d] = (ic50[d]/sens_scale[p]) * (1 - delta*b[d]) * ref
Emax (ceiling)      : capped by binding sequestration
                         Emax[d]   = Emax0 * (1 - gamma*b[d])
kill (Hill, h=1.3)  : Emax[d] * AUC^h / (EC50^h + AUC^h)  + mild seeded noise
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import selection_lib

c = selection_lib.load_cache()
vasc, ecm, sens = c["vasc"], c["ecm"], c["sens"]
mw, ic50, bbb = c["drug_mw"], c["drug_ic50"], c["drug_bbb"]
dose = float(c["dose"])
P, D = len(vasc), len(mw)

rng = np.random.default_rng(7)

# --- drug binding: strongly-bound == large molecule (min-max normalised MW) ---
b = (mw - mw.min()) / (mw.max() - mw.min())          # [0,1], doxorubicin~1, temo~0

# --- local delivered exposure (AUC of free drug at the tumour) ---
# perfusion delivers, BBB gates entry, dense ECM blocks diffusion (worse for large drugs)
k_ecm = 2.0
AUC = dose * bbb[None, :] * vasc[:, None] * np.exp(-k_ecm * ecm[:, None] * (0.5 + b[None, :]))

# --- saturating dose-response parameters ---
# EC50 in exposure units: less potent drug (high ic50) or less sensitive patient -> higher EC50.
# Normalise ic50 to O(1) so EC50 sits on the same scale as AUC; binding lowers EC50 (plateau early).
ic50n = ic50 / np.median(ic50)                       # dimensionless potency
sens_scale = sens / sens.mean()                      # patient sensitivity (>1 == more sensitive)
ref = 0.6                                             # exposure scale so some drugs saturate, some don't
delta = 0.35                                          # binding pulls EC50 down (earlier half-max)
EC50 = (ic50n[None, :] / sens_scale[:, None]) * (1.0 - delta * b[None, :]) * ref
EC50 = np.clip(EC50, 1e-4, None)

# Emax ceiling capped by tissue sequestration of strongly-bound drugs
gamma = 0.45
Emax = 1.0 * (1.0 - gamma * b)                        # per-drug plateau in [1-gamma, 1]

h = 1.3                                               # Hill coefficient
kill = Emax[None, :] * AUC**h / (EC50**h + AUC**h)

# mild, fixed-seed multiplicative noise, then clip to [0,1]
kill = kill * (1.0 + 0.03 * rng.standard_normal((P, D)))
gt_kill = np.clip(kill, 0.0, 1.0)

result = selection_lib.score_policies(gt_kill, seed=0)

# quick internal diagnostics to stderr (not parsed)
print("mean kill per drug:", np.round(gt_kill.mean(0), 3), file=sys.stderr)
print("mean AUC per drug :", np.round(AUC.mean(0), 3), file=sys.stderr)
print("mean EC50 per drug:", np.round(EC50.mean(0), 3), file=sys.stderr)
print("binding b         :", np.round(b, 3), file=sys.stderr)
print("best_distribution :", result["diag"]["best_distribution"], file=sys.stderr)

selection_lib.print_json(result)
