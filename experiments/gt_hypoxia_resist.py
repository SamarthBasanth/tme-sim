"""
Ground-truth motif #3: "hypoxia_resist"

Biological assumption
---------------------
In GBM, a drug's realised kill is the product of three gates:

  1. DELIVERY: the fraction of dose that actually reaches viable tumor cells.
     Drug must be perfused in from vasculature (vasc), diffuse through the
     extracellular matrix (ecm barrier), cross the blood-brain barrier
     (drug_bbb), and is size-limited (drug_mw). Crucially, HYPOXIC tumor
     regions are poorly perfused, so hypoxia further throttles delivery.

  2. POTENCY: intrinsic cytotoxicity, a Hill term dose/(dose+ic50) scaled by
     the patient's intrinsic drug sensitivity (sens).

  3. HYPOXIA-RESISTANCE GATE (the motif): hypoxic cells resist MOST drugs
     because they are quiescent (out of cycle) and up-regulate survival/efflux
     programs and lack the O2/ROS many agents need. But this is DRUG-SPECIFIC:
     oxygen-independent alkylators (temozolomide, lomustine/CCNU) keep working
     in hypoxia, whereas O2/proliferation-dependent agents (cisplatin,
     methotrexate, doxorubicin) are strongly blunted. So the gate is
     (1 - h_sens[d] * hypoxic[p]) with a small per-drug h_sens for alkylators
     and a large one for the rest.

Consequence: in strongly hypoxic patients the intrinsically most-potent drug
(doxorubicin, ic50=0.5) is doubly punished -- it barely crosses the BBB
(bbb=0.05) AND it is highly hypoxia-sensitive -- so an oxygen-independent,
BBB-penetrant nitrosourea (lomustine, bbb=0.97) becomes the right pick. This is
a delivery/biology-limited regime the mechanistic sim COULD win, unless the sim
already captures delivery well enough that plain features suffice.

Key parameters
--------------
  delivery = bbb^1 * mw_pen * ecm_barrier * perfusion
    mw_pen      = exp(-(mw-194)/300)          (size-limited diffusion)
    ecm_barrier = exp(-2.0*ecm)               (matrix impedance)
    perfusion   = vasc * (1 - 0.8*hypoxic)    (hypoxia starves perfusion)
  potency  = sens * dose/(dose+ic50)          (Hill)
  hyp_gate = 1 - h_sens[d]*hypoxic            (drug-specific hypoxia resistance)
    h_sens = [TMZ 0.15, lomustine 0.10, cisplatin 0.85, MTX 0.70, doxo 0.90]
  gt = normalize(delivery_norm * potency * hyp_gate) + N(0, 0.02^2), clipped.
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import selection_lib


def build_gt(cache):
    vasc = cache["vasc"]
    ecm = cache["ecm"]
    sens = cache["sens"]
    feat = cache["feat"]
    feat_names = list(cache["feat_names"])
    hypoxic = feat[:, feat_names.index("hypoxic")]

    mw = cache["drug_mw"]
    ic50 = cache["drug_ic50"]
    bbb = cache["drug_bbb"]
    dose = float(cache["dose"])
    panel = list(cache["panel"])

    P = len(vasc)
    D = len(panel)

    # --- per-drug hypoxia sensitivity (the motif) ---
    # oxygen-independent alkylators keep working in hypoxia -> small h_sens.
    h_sens_map = {
        "temozolomide": 0.15,
        "lomustine": 0.10,
        "cisplatin": 0.85,
        "methotrexate": 0.70,
        "doxorubicin": 0.90,
    }
    h_sens = np.array([h_sens_map[name] for name in panel])

    # --- delivery term (P,D) ---
    mw_pen = np.exp(-(mw - mw.min()) / 300.0)          # (D,)
    ecm_barrier = np.exp(-2.0 * ecm)                    # (P,)
    perfusion = vasc * (1.0 - 0.8 * hypoxic)            # (P,) hypoxia starves perfusion
    delivery = (perfusion[:, None]
                * ecm_barrier[:, None]
                * bbb[None, :]
                * mw_pen[None, :])                      # (P,D)
    delivery_norm = delivery / delivery.max()

    # --- potency term (P,D): Hill x patient sensitivity ---
    hill = dose / (dose + ic50)                         # (D,)
    potency_t = sens[:, None] * hill[None, :]           # (P,D)

    # --- hypoxia-resistance gate (P,D) ---
    hyp_gate = 1.0 - h_sens[None, :] * hypoxic[:, None]  # (P,D), in (0,1]

    raw = delivery_norm * potency_t * hyp_gate

    rng = np.random.default_rng(3)  # motif #3, fixed seed
    raw = raw + rng.normal(0.0, 0.02, size=raw.shape)
    raw = np.clip(raw, 0.0, None)
    gt = raw / raw.max()
    gt = np.clip(gt, 0.0, 1.0)
    return gt


def main():
    cache = selection_lib.load_cache()
    gt = build_gt(cache)
    result = selection_lib.score_policies(gt, seed=0)
    selection_lib.print_json(result)


if __name__ == "__main__":
    main()
