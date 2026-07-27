"""Ground-truth motif #2: SHARP DELIVERY THRESHOLD.

Biological assumption
---------------------
In GBM, a systemically dosed drug only exerts cytotoxicity if its *perivascular
delivery* -- the local concentration it can build up near the tumor vasculature
after crossing the blood-brain barrier and diffusing through the extracellular
matrix -- exceeds a sharp cutoff. Below that cutoff the drug simply never reaches
a cytotoxic concentration, so kill collapses to ~0 REGARDLESS of intrinsic
potency. Above it, kill is governed by potency (IC50 vs dose, tumor sensitivity).

This is the classic "potent but undeliverable" trap: doxorubicin is by far the
most potent agent in the panel (IC50 0.5) but has the lowest BBB permeability
(0.051) and the highest molecular weight (543 Da), so its perivascular delivery
sits below threshold in most patients and it is effectively useless -- while a
moderately potent, well-delivered agent (e.g. lomustine, BBB 0.97) wins.

Model
-----
delivery[p,d]  built from observable transport determinants:
    + vasc[p]        (vascular density -> more delivery surface)
    + drug_bbb[d]    (BBB permeability -> more drug crosses)
    - ecm[p]         (dense matrix -> slower interstitial diffusion)
    - drug_mw[d]     (large molecules -> slower diffusion / lower flux)

  raw = vasc * bbb * exp(-beta_ecm*(ecm-ecm0)) * exp(-beta_mw*(mw-mw0)/mw0)
  then z-scored across the whole P x D panel so the cutoff is interpretable.

gate[p,d] = sigmoid(k * (delivery_z - cutoff))     # SHARP: large k
  -> ~1 above threshold, ~0 below. This is the motif's defining nonlinearity.

potency[p,d] = sens[p] * dose/(dose + ic50[d])      # Hill-ish, delivery-agnostic

gt_kill[p,d] = gate * potency  (+ mild seeded noise), clipped to [0,1].

NOTE: K_sim is NEVER used to build the ground truth (that would be circular).
"""
import os
import numpy as np

import selection_lib

SEED = 2


def main():
    c = selection_lib.load_cache()
    vasc = c["vasc"].astype(float)          # (P,)
    ecm = c["ecm"].astype(float)            # (P,)
    sens = c["sens"].astype(float)          # (P,)
    mw = c["drug_mw"].astype(float)         # (D,)
    ic50 = c["drug_ic50"].astype(float)     # (D,)
    bbb = c["drug_bbb"].astype(float)       # (D,)
    dose = float(c["dose"])
    P, D = len(vasc), len(mw)

    # ---- delivery term (perivascular concentration proxy) --------------------
    ecm0, mw0 = 0.60, 300.0
    beta_ecm, beta_mw = 6.0, 1.2
    ecm_pen = np.exp(-beta_ecm * (ecm - ecm0))            # (P,)
    mw_pen = np.exp(-beta_mw * (mw - mw0) / mw0)          # (D,)
    raw = (vasc[:, None] * ecm_pen[:, None]) * (bbb[None, :] * mw_pen[None, :])

    # z-score across the whole panel so the threshold cutoff is meaningful
    delivery_z = (raw - raw.mean()) / (raw.std() + 1e-12)

    # ---- SHARP threshold gate ------------------------------------------------
    k = 6.0          # steepness: large -> near step function
    cutoff = 0.15    # slightly above panel mean -> a real fraction fall below
    gate = 1.0 / (1.0 + np.exp(-k * (delivery_z - cutoff)))

    # ---- potency term (delivery-agnostic) ------------------------------------
    potency_eff = sens[:, None] * (dose / (dose + ic50[None, :]))

    # ---- combine + mild seeded noise ----------------------------------------
    rng = np.random.default_rng(SEED)
    noise = rng.normal(0.0, 0.02, size=(P, D))
    gt_kill = np.clip(gate * potency_eff + noise, 0.0, 1.0)

    result = selection_lib.score_policies(gt_kill, seed=0)

    # a few motif sanity numbers on stderr-ish (still stdout, non RESULT line)
    print("delivery_z range: %.2f .. %.2f  frac_below_cutoff=%.2f" % (
        delivery_z.min(), delivery_z.max(), float((delivery_z < cutoff).mean())))
    print("mean gate per drug:", np.round(gate.mean(0), 3).tolist(),
          "(", c["panel"].tolist(), ")")
    print("mean potency per drug:", np.round(potency_eff.mean(0), 3).tolist())
    selection_lib.print_json(result)


if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    main()
