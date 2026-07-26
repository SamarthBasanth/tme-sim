"""Ground-truth motif #5: BBB-gated drug delivery in GBM.

Biological assumption
---------------------
The blood-brain barrier (BBB) is the dominant gate on drug delivery in
glioblastoma. A drug only reaches tumor cells if EITHER it is intrinsically
BBB-penetrant (drug_bbb high) OR the tumor vasculature is leaky enough to
disrupt the barrier. Aggressive, highly-vascularized GBMs have a leaky /
neovascular BBB, so even poorly-penetrant drugs can seep in. Well-differentiated
tumors (low vascularization) keep an intact BBB, so low-BBB drugs -- however
potent in a dish -- never reach a lethal concentration and fail.

A secondary, motif-consistent delivery penalty: once past the vessel wall, a
drug must diffuse through dense ECM to reach cells. Large molecules (high
drug_mw) are hindered more strongly in stiff / dense matrix (high ecm), a
Stokes-Einstein-style size sieving. This further penalizes the big, potent,
low-BBB agents (doxorubicin, methotrexate).

Kill = DELIVERY(BBB gate x ECM size-sieve) x POTENCY(Hill dose-response x sens).
This is intentionally a delivery-GATED product: potency is necessary but not
sufficient -- delivery is the limiting factor the motif is about.

Key parameters
--------------
  vasc_mid=0.34, vasc_scale=0.03  -> leakiness sigmoid over the cohort's vasc range
  k_ecm=1.2                        -> ECM x size sieving strength
  dose=30 (from cache)             -> Hill saturating potency dose/(dose+ic50)
  noise_std=0.02, seed=1234
"""
import numpy as np
import selection_lib


def build_gt_kill():
    c = selection_lib.load_cache()
    vasc = c["vasc"]            # (P,)
    ecm = c["ecm"]             # (P,)
    sens = c["sens"]           # (P,)
    bbb = c["drug_bbb"]        # (D,)
    mw = c["drug_mw"]          # (D,)
    ic50 = c["drug_ic50"]      # (D,)
    dose = float(c["dose"])
    P = vasc.shape[0]
    D = bbb.shape[0]

    # --- DELIVERY term 1: BBB gate modulated by vascular leakiness ---
    # leak in [0,1]: well-differentiated (low vasc) -> intact BBB -> ~0
    #                aggressive/neovascular (high vasc) -> leaky -> ~1
    vasc_mid, vasc_scale = 0.34, 0.03
    leak = 1.0 / (1.0 + np.exp(-(vasc - vasc_mid) / vasc_scale))   # (P,)

    # access[p,d]: a drug reaches tumor if it is BBB-penetrant OR barrier is leaky.
    # High-BBB drugs are always ~fully delivered; low-BBB drugs ride on leak[p].
    access = bbb[None, :] + (1.0 - bbb[None, :]) * leak[:, None]   # (P,D) in [0,1]

    # --- DELIVERY term 2: ECM size-sieving (large molecules hindered in dense ECM) ---
    mw_norm = (mw - mw.min()) / (mw.max() - mw.min())             # (D,) in [0,1]
    k_ecm = 1.2
    ecm_sieve = np.exp(-k_ecm * ecm[:, None] * mw_norm[None, :])   # (P,D) in (0,1]

    delivery = access * ecm_sieve                                 # (P,D)

    # --- POTENCY term: Hill saturating dose-response scaled by patient sensitivity ---
    hill = dose / (dose + ic50)                                   # (D,) in (0,1)
    potency = sens[:, None] * hill[None, :]                       # (P,D)

    # --- GATED kill ---
    kill = delivery * potency                                    # (P,D) in [0,1]

    rng = np.random.default_rng(1234)
    kill = kill + rng.normal(0.0, 0.02, size=kill.shape)
    kill = np.clip(kill, 0.0, 1.0)
    return kill


if __name__ == "__main__":
    gt = build_gt_kill()
    result = selection_lib.score_policies(gt, seed=0)
    selection_lib.print_json(result)
