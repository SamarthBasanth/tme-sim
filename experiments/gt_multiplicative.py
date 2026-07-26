"""
Ground-truth motif #0: MULTIPLICATIVE (generic delivery x potency).

Biological assumption
---------------------
For each GBM patient p and drug d, the fraction of tumour cells killed is the
INDEPENDENT PRODUCT of three separable factors -- there is no feedback between
delivery and the pharmacodynamic response:

    kill[p,d] = sens_factor(p) * delivery(p,d) * Hill(dose, ic50[d]) * noise

1. sens_factor(p) = sens[p]
   Intrinsic tumour chemosensitivity (patient-level, drug-agnostic).

2. delivery(p,d) = perfusion(p) * bbb[d] * ecm_penetration(p,d)
   How much drug actually reaches the tumour cells:
     - perfusion(p): more vascularised tumours are better perfused
       (blood-borne drug arrives). perfusion = clip(vasc / 0.45, 0, 1).
     - bbb[d]: fraction of drug crossing the blood-brain barrier
       (drug_bbb, an intrinsic molecular property).
     - ecm_penetration(p,d): interstitial transport through dense ECM.
       Large molecules diffuse slowly (Stokes-Einstein: D ~ 1/radius ~ 1/mw^(1/3)),
       and dense ECM compounds the penalty. Modelled as
       exp(-beta * ecm[p] * mw_norm[d]) with mw_norm = mw/mw_max, beta = 2.2.

3. Hill(dose, ic50[d]) = dose / (dose + ic50[d])   (Hill coefficient h = 1)
   Concentration-response potency at the delivered site. dose is fixed (30).

This is deliberately the canonical multiplicative form: delivery and potency
are orthogonal factors that multiply. The tension it creates is real -- the
most POTENT drug (doxorubicin, ic50=0.5) is the WORST delivered (bbb=0.05,
largest MW), so the globally-most-potent drug is frequently NOT the best drug
per patient. That is exactly the delivery-limited regime a transport-aware
mechanistic sim is supposed to exploit, so this motif is a fair test where the
sim CAN legitimately win over the potency shortcut.

Everything is deterministic given seed=0.
"""
import numpy as np
import selection_lib as sl


def build_gt_kill(cache):
    vasc = cache["vasc"]            # (P,)
    ecm = cache["ecm"]             # (P,)
    sens = cache["sens"]            # (P,)
    mw = cache["drug_mw"]        # (D,)
    ic50 = cache["drug_ic50"]      # (D,)
    bbb = cache["drug_bbb"]       # (D,)
    dose = float(cache["dose"])

    P = vasc.shape[0]
    D = mw.shape[0]

    # ---- delivery term ----
    perfusion = np.clip(vasc / 0.45, 0.0, 1.0)          # (P,)
    mw_norm = mw / mw.max()                              # (D,)
    beta = 2.2
    # ecm_penetration[p,d] = exp(-beta * ecm[p] * mw_norm[d])
    ecm_pen = np.exp(-beta * np.outer(ecm, mw_norm))     # (P,D)
    delivery = perfusion[:, None] * bbb[None, :] * ecm_pen  # (P,D)

    # ---- potency term (Hill, h=1) ----
    hill = dose / (dose + ic50)                          # (D,)

    # ---- multiplicative combination ----
    sens_factor = sens[:, None]                          # (P,1)
    kill = sens_factor * delivery * hill[None, :]        # (P,D)

    # mild multiplicative log-normal noise, fixed seed
    rng = np.random.default_rng(0)
    noise = np.exp(rng.normal(0.0, 0.06, size=(P, D)))
    kill = kill * noise

    return np.clip(kill, 0.0, 1.0)


def main():
    cache = sl.load_cache()
    gt_kill = build_gt_kill(cache)
    result = sl.score_policies(gt_kill, seed=0)
    sl.print_json(result)


if __name__ == "__main__":
    main()
