"""GROUND-TRUTH motif #1: Liebig's law of the minimum.

Biological assumption
---------------------
A cytotoxic drug can only kill GBM cells to the extent permitted by the SINGLE
WORST-limiting factor along the causal chain, not their product. Two independent
bottlenecks:

  1. DELIVERY term  -- what fraction of administered drug actually reaches tumor
     cells. Governed by (a) BBB/blood-tumor-barrier penetrance (drug_bbb),
     (b) vascular perfusion that carries drug to the region (vasc), and
     (c) interstitial diffusion through dense ECM, which penalises large
     molecules more (ecm x molecular weight).

  2. POTENCY term   -- pharmacodynamic efficacy of the drug that DID arrive,
     i.e. a Hill response of the delivered exposure vs the drug's IC50.

Liebig: the achievable effect is capped by whichever term is worse
      kill = sens * min(delivery_term, potency_term)
A super-potent drug that cannot cross the BBB (doxorubicin) is wasted; a highly
deliverable but pharmacologically weak drug (temozolomide) is equally wasted.
Overall tumor chemosensitivity (sens) scales the ceiling.

This is a genuine minimum, NOT the multiplicative sim form -- the two bottlenecks
do not compound smoothly; one of them alone dictates the outcome.

Key parameters: Hill coeff h=1.3, dose=30 (from cache), ECM-diffusion rate
k_ecm=1.6, size reference mw_ref=min(mw), perfusion reference vasc_ref=0.40,
delivery gain G=2.2 (calibrates delivered exposure into the physiological range
where the min-bottleneck genuinely switches drug-to-drug), noise sigma=0.02
(seed=1).
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import selection_lib

H = 1.3
K_ECM = 1.6
VASC_REF = 0.40
GAIN = 2.2
NOISE_SIGMA = 0.02
SEED = 1


def hill(exposure, ic50, h=H):
    exposure = np.maximum(exposure, 0.0)
    return exposure**h / (exposure**h + ic50**h)


def main():
    c = selection_lib.load_cache()
    vasc, ecm, sens = c["vasc"], c["ecm"], c["sens"]
    mw, ic50, bbb = c["drug_mw"], c["drug_ic50"], c["drug_bbb"]
    dose = float(c["dose"])
    P, D = len(vasc), len(mw)

    mw_ref = float(mw.min())
    size = mw / mw_ref                      # >=1, larger molecules diffuse worse

    # --- DELIVERY term (P,D) in [0,1]: fraction of drug reaching tumor cells ---
    # perfusion factor (per patient), saturating around ref perfusion
    perf = np.clip(vasc / VASC_REF, 0.0, 1.0)                       # (P,)
    # ECM interstitial penetration, size-weighted (per patient x drug)
    penet = np.exp(-K_ECM * np.outer(ecm, size))                    # (P,D)
    delivery = GAIN * bbb[None, :] * perf[:, None] * penet           # (P,D)
    delivery = np.clip(delivery, 0.0, 1.0)

    # --- POTENCY term (1,D) broadcast: PD efficacy of what arrives ---
    potency_term = hill(dose, ic50)[None, :]                        # (1,D)

    # --- Liebig minimum, scaled by chemosensitivity ---
    kill = sens[:, None] * np.minimum(delivery, potency_term)       # (P,D)

    rng = np.random.default_rng(SEED)
    kill = kill + rng.normal(0.0, NOISE_SIGMA, size=kill.shape)
    gt_kill = np.clip(kill, 0.0, 1.0)

    # brief provenance so the run is auditable
    panel = c["panel"]
    print("delivery col means:", {panel[d]: round(float(delivery[:, d].mean()), 3) for d in range(D)})
    print("potency term:      ", {panel[d]: round(float(potency_term[0, d]), 3) for d in range(D)})
    binding = np.array(["deliv" if delivery[p, d] < potency_term[0, d] else "poten"
                        for p in range(P) for d in range(D)]).reshape(P, D)
    print("delivery is the bottleneck fraction:", round(float((binding == "deliv").mean()), 3))

    result = selection_lib.score_policies(gt_kill, seed=0)
    selection_lib.print_json(result)


if __name__ == "__main__":
    main()
