"""Ground-truth motif #9: DELIVERY-DOMINANT drug response in GBM.

Biological assumption
---------------------
In glioblastoma the rate-limiting step for a cytotoxic agent is *getting to the
tumor cell*: crossing the blood-brain barrier and then diffusing through dense
extracellular matrix. Once a drug actually arrives, the delivered concentration
in these regimens sits near the SATURATING top of the cell-kill dose-response
curve, so the intrinsic in-vitro potency (ic50) barely distinguishes one agent
from another at the target -- the classic "the drug that works is the drug that
gets there" picture. We encode this by COMPRESSING the ic50 spread with a strong
concave transform (gamma small), so the potency term is nearly flat across the
panel. The winner for a given patient is therefore set almost entirely by the
DELIVERY term:

  delivery = ACCESS(BBB penetrance OR leaky neovasculature)  x  ECM size-sieve.

  - ACCESS[p,d]: a drug reaches the tumor if it is intrinsically BBB-penetrant
    (drug_bbb high) OR the patient's vasculature is leaky enough (high vasc) to
    have disrupted the barrier. access = bbb + (1-bbb)*leak(vasc).
  - ECM size-sieve[p,d]: after extravasation the drug must percolate through the
    matrix; large molecules (high drug_mw) are hindered ~exponentially in dense
    ECM (high ecm). Stokes-Einstein-style hindered diffusion.

  kill = delivery x potency_compressed x sens

sens is a patient-level scalar (scales all drugs equally, does not change the
per-patient winner). potency_compressed is nearly flat across drugs by design,
so the per-patient argmax is driven by delivery -- exactly the motif. This is a
regime where a mechanistic transport PDE (the sim) SHOULD have an edge over a
delivery-agnostic potency shortcut.

Key parameters
--------------
  vasc_mid=0.34, vasc_scale=0.03  -> leakiness sigmoid over the cohort vasc range
  k_ecm=1.5                        -> ECM x molecular-size sieving strength
  pmin=0.90                        -> potency compressed into a narrow [0.90,1.0]
                                       band (ic50 spread flattened -> does not decide)
  dose=30 (from cache)             -> Hill saturating dose/(dose+ic50) before compression
  noise_std=0.02, seed=909
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

    # ---------- DELIVERY (the dominant term) ----------
    # 1) BBB access gated by vascular leakiness.
    vasc_mid, vasc_scale = 0.34, 0.03
    leak = 1.0 / (1.0 + np.exp(-(vasc - vasc_mid) / vasc_scale))        # (P,)
    access = bbb[None, :] + (1.0 - bbb[None, :]) * leak[:, None]        # (P,D) in [0,1]

    # 2) ECM size-sieving: big molecules hindered in dense matrix.
    mw_norm = (mw - mw.min()) / (mw.max() - mw.min())                   # (D,) in [0,1]
    k_ecm = 1.5
    ecm_sieve = np.exp(-k_ecm * ecm[:, None] * mw_norm[None, :])        # (P,D) in (0,1]

    delivery = access * ecm_sieve                                      # (P,D), wide spread

    # ---------- POTENCY (compressed / flattened so it does NOT decide) ----------
    # Faithful to the motif: flatten the ic50 spread into a narrow band [pmin,1].
    # Raw saturating Hill, min-max normalized, then squeezed so the whole panel
    # differs by at most (1-pmin) in potency -- delivery (which spans ~0.02..1.0)
    # dwarfs it and decides the per-patient winner.
    hill = dose / (dose + ic50)                                        # (D,) in (0,1)
    hill_norm = (hill - hill.min()) / (hill.max() - hill.min())        # (D,) in [0,1]
    pmin = 0.90
    potency_drug = pmin + (1.0 - pmin) * hill_norm                     # (D,) in [0.90,1.0]
    potency = sens[:, None] * potency_drug[None, :]                    # (P,D)

    # ---------- KILL ----------
    kill = delivery * potency                                         # (P,D) in [0,1]

    rng = np.random.default_rng(909)
    kill = kill + rng.normal(0.0, 0.02, size=kill.shape)
    kill = np.clip(kill, 0.0, 1.0)
    return kill


if __name__ == "__main__":
    gt = build_gt_kill()
    # quick self-diagnostic on how flat potency is vs how spread delivery is
    result = selection_lib.score_policies(gt, seed=0)
    selection_lib.print_json(result)
