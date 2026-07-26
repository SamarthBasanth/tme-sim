"""
Ground-truth motif #10: NOISY_REALISTIC (messy additive mixture, no dominant factor).

Biological assumption
---------------------
Real per-patient GBM drug response is not governed by any single clean law. It is
a MIXTURE of comparable-magnitude contributions plus appreciable biological /
measurement noise. No one factor (delivery, potency, or subtype) dominates, so the
best drug for a patient emerges from the *interaction* of several moderate effects.

We combine three normalised terms ADDITIVELY (not multiplicatively) with roughly
balanced weights, so that none can by itself decide the argmax:

    kill[p,d] = clip( w_del*Del~ + w_pot*Pot~ + w_sub*Sub~ + noise , 0, 1 )

1. Delivery  Del(p,d) = perfusion(p) * bbb[d] * exp(-beta*ecm[p]*mw_norm[d])
     - perfusion(p) = clip(vasc/0.45, 0, 1): vascularised tumours are better perfused.
     - bbb[d]: intrinsic blood-brain-barrier permeability of the drug.
     - ecm penetration: large drugs (Stokes-Einstein D~1/mw^(1/3)) penetrate dense
       ECM poorly; modelled exp(-beta*ecm*mw_norm), beta=2.0.
   Del is min-max normalised across the whole matrix -> Del~ in [0,1].

2. Potency   Pot(p,d) = sens[p] * Hill(dose, ic50[d]),  Hill = dose/(dose+ic50),
   the delivery-agnostic pharmacodynamic response scaled by patient chemosensitivity.
   Min-max normalised -> Pot~ in [0,1].

3. Subtype   Sub(p,d) = affinity[d, subtype(p)], a loose drug x GBM-subtype
   preference matrix (proneural / classical / mesenchymal). Encodes that, e.g.,
   alkylators track proneural/classical while platinum leans mesenchymal.
   Min-max normalised -> Sub~ in [0,1].

Weights w_del=0.40, w_pot=0.35, w_sub=0.25 -- deliberately balanced so the motif's
defining property (no single dominant factor) holds. Then MODERATE additive Gaussian
noise (sigma=0.10) is applied with a fixed seed, mimicking assay/biological scatter.

Because delivery and subtype both matter and doxorubicin (most potent, ic50=0.5) is
badly delivered (bbb=0.05, largest MW), the globally-most-potent drug is usually NOT
the per-patient best -- but neither is delivery alone decisive. This is a genuinely
mixed, noisy regime: the mechanistic sim (transport+potency) captures two of the
three signals but is blind to the subtype term and the noise, while the feature-ML
baseline sees subtype one-hots + vasc/ecm/sens but must learn from only 30 patients.
A fair, honest test with no built-in winner.

Everything is deterministic given seed=0.
"""
import numpy as np
import selection_lib as sl


def _minmax(x):
    lo, hi = float(x.min()), float(x.max())
    if hi - lo < 1e-12:
        return np.zeros_like(x)
    return (x - lo) / (hi - lo)


def build_gt_kill(cache):
    vasc = cache["vasc"]        # (P,)
    ecm = cache["ecm"]          # (P,)
    sens = cache["sens"]        # (P,)
    mw = cache["drug_mw"]       # (D,)
    ic50 = cache["drug_ic50"]   # (D,)
    bbb = cache["drug_bbb"]     # (D,)
    dose = float(cache["dose"])
    feat = cache["feat"]        # (P,F); cols 4,5,6 = proneural, classical, mesenchymal
    subtype = feat[:, 4:7].argmax(1)   # (P,) in {0,1,2}

    P = vasc.shape[0]
    D = mw.shape[0]

    # ---- (1) delivery term ----
    perfusion = np.clip(vasc / 0.45, 0.0, 1.0)             # (P,)
    mw_norm = mw / mw.max()                                 # (D,)
    beta = 2.0
    ecm_pen = np.exp(-beta * np.outer(ecm, mw_norm))        # (P,D)
    delivery = perfusion[:, None] * bbb[None, :] * ecm_pen  # (P,D)
    Del = _minmax(delivery)

    # ---- (2) potency term ----
    hill = dose / (dose + ic50)                             # (D,)
    potency = sens[:, None] * hill[None, :]                 # (P,D)
    Pot = _minmax(potency)

    # ---- (3) subtype affinity term ----
    # rows = drugs [tmz, lomustine, cisplatin, methotrexate, doxorubicin]
    # cols = subtype [proneural, classical, mesenchymal]
    affinity = np.array([
        [0.90, 0.70, 0.40],   # temozolomide  -> proneural/classical (alkylator)
        [0.60, 0.85, 0.60],   # lomustine     -> classical
        [0.45, 0.55, 0.90],   # cisplatin     -> mesenchymal (platinum)
        [0.70, 0.50, 0.65],   # methotrexate  -> proneural-ish
        [0.60, 0.70, 0.75],   # doxorubicin   -> mesenchymal-ish
    ])                                                      # (D,3)
    subtype_term = affinity[:, subtype].T                   # (P,D)
    Sub = _minmax(subtype_term)

    # ---- balanced additive blend (no single dominant factor) ----
    w_del, w_pot, w_sub = 0.40, 0.35, 0.25
    score = w_del * Del + w_pot * Pot + w_sub * Sub         # (P,D), ~[0,1]

    # ---- moderate additive noise, fixed seed ----
    rng = np.random.default_rng(0)
    score = score + rng.normal(0.0, 0.10, size=(P, D))

    return np.clip(score, 0.0, 1.0)


def main():
    cache = sl.load_cache()
    gt_kill = build_gt_kill(cache)
    result = sl.score_policies(gt_kill, seed=0)
    sl.print_json(result)


if __name__ == "__main__":
    main()
