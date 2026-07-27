"""Ground-truth motif #4: ECM sieving (size-dependent exclusion).

Biological assumption
---------------------
GBM interstitium is a dense, tortuous ECM mesh. Drug transport into the tumor
core is governed by *steric sieving*: a molecule of hydrodynamic radius r must
diffuse through pores whose effective size collapses as ECM density rises.
Large drugs (high MW) are therefore excluded from dense-ECM tumors almost
independently of how potent they are -- a highly potent large molecule that
cannot physically reach the cells kills nothing. Small drugs pass the sieve
freely regardless of ECM density.

Model
-----
kill[p,d] = potency[p,d] * delivery[p,d] * noise

- size proxy  s_d in [0,1] from drug_mw (min-max); r ~ MW so we use MW directly.
- steric sieving penetration (dominant delivery term, motif core):
      sieve[p,d] = exp( -LAM * s_d * ecm_p^ECM_POW )
  large s_d AND high ecm -> exponential exclusion; small s_d -> sieve ~ 1
  independent of ecm. This is the size-dependent ECM sieve.
- a mild vascular/BBB access factor (does NOT drive per-patient drug ranking
  much; vasc is drug-independent so it cannot change the per-patient argmax;
  bbb kept small so the ECM sieve, not BBB, is the delivery mechanism):
      access[p,d] = (0.5 + 0.5*vasc_p) * (0.7 + 0.3*bbb_d)
- potency: saturating Hill in dose vs ic50, scaled by intrinsic sensitivity:
      potency[p,d] = sens_p * dose/(dose + ic50_d)

Because the sieve term has a large dynamic range keyed to MW, the best drug in a
dense-ECM patient is a smaller, less-potent molecule -- delivery-limited by
design. This is a motif where a transport-aware sim SHOULD do well; we report
the sim_vs_gt rank correlation honestly as a circularity check.
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import selection_lib

LAM = 6.0        # sieve strength
ECM_POW = 2.0    # ECM density sharpens the sieve
SEED = 4

def build_gt():
    c = selection_lib.load_cache()
    vasc, ecm, sens = c["vasc"], c["ecm"], c["sens"]
    mw, ic50, bbb = c["drug_mw"], c["drug_ic50"], c["drug_bbb"]
    dose = float(c["dose"])
    P, D = len(vasc), len(mw)

    # drug size proxy in [0,1]
    s = (mw - mw.min()) / (mw.max() - mw.min())          # (D,)

    # steric ECM sieve: exp(-LAM * size * ecm^pow)  -> (P,D)
    sieve = np.exp(-LAM * s[None, :] * (ecm[:, None] ** ECM_POW))

    # mild vascular / BBB access (drug-ranking-neutral in vasc; small bbb effect)
    access = (0.5 + 0.5 * vasc)[:, None] * (0.7 + 0.3 * bbb)[None, :]

    delivery = sieve * access                             # (P,D)

    # saturating potency (Hill), scaled by patient sensitivity
    potency = sens[:, None] * (dose / (dose + ic50))[None, :]   # (P,D)

    kill = potency * delivery

    # mild multiplicative log-normal noise, fixed seed
    rng = np.random.default_rng(SEED)
    kill = kill * np.exp(0.05 * rng.standard_normal(kill.shape))

    # normalize into [0,1]
    kill = np.clip(kill, 0.0, None)
    kill = kill / kill.max()
    return kill

def main():
    gt = build_gt()
    res = selection_lib.score_policies(gt, seed=0)
    selection_lib.print_json(res)

if __name__ == "__main__":
    main()
