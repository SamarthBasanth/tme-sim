"""Adversarial ground truth: POTENCY-SHORTCUT-WINS regime.

Goal of this probe
------------------
Find a biologically plausible kill matrix in which drug DELIVERY barely matters,
so that the delivery-agnostic POTENCY shortcut is near-oracle and the expensive
mechanistic transport sim adds NOTHING (indeed loses). This maps the boundary
where the sim stops paying for itself.

Biological assumption
---------------------
A well-perfused, thin-stroma tumour cohort. The neovasculature is leaky and the
extracellular matrix is not dense, so EVERY agent in the panel extravasates and
percolates to broadly therapeutic fractions -- the delivery term saturates near 1
for all drugs and acts only as a narrow multiplier. Meanwhile the panel's
intrinsic pharmacodynamic potency spans ~3 orders of magnitude in IC50
(doxorubicin 0.5 -> temozolomide 606 uM), so the Hill dose-response term is
widely separated between drugs. When delivery is a ~[dmin,1] multiplier and Hill
values differ by more than that band, the per-patient drug ranking is fixed by
potency alone. This is the classic "the drugs are so different in potency that
where they go can't reorder them" picture -- the mirror image of the
delivery-limited GBM regime.

    kill[p,d] = sens[p] * Hill(dose, ic50[d]) * delivery[p,d] * noise

  Hill(dose,ic50) = dose/(dose+ic50)          -- DOMINANT, drug-specific term
  delivery[p,d]   in [dmin, 1.0]               -- MINOR multiplier, real transport
                                                  biology (vasc perfusion, BBB,
                                                  ECM size-sieve) compressed into
                                                  a narrow band.
  sens[p]         -- patient-level chemosensitivity, does not reorder drugs.
  hypoxic[p]      -- mild patient-level resistance (drug-agnostic), does not
                     reorder drugs; included for biological realism only.

Because doxorubicin's Hill value (0.98) so dominates, even at the worst delivery
(dmin) it out-kills the next drug at perfect delivery, so the globally-most-potent
drug is the per-patient best for nearly everyone -> potency shortcut ~ oracle.

NOTE: K_sim is NEVER read here. The ground truth is built only from patient
biology (vasc/ecm/sens/hypoxic/subtype) and drug properties (mw/ic50/bbb).

Deterministic given seed.
"""
import json
import numpy as np
import selection_lib as sl


def build_gt_kill(cache, dmin=0.85, k_ecm=0.8, noise_std=0.03, seed=0):
    vasc = cache["vasc"]            # (P,)
    ecm = cache["ecm"]             # (P,)
    sens = cache["sens"]           # (P,)
    feat = cache["feat"]           # (P,7): [...,3] = hypoxic
    mw = cache["drug_mw"]          # (D,)
    ic50 = cache["drug_ic50"]      # (D,)
    bbb = cache["drug_bbb"]        # (D,)
    dose = float(cache["dose"])
    hypoxic = feat[:, 3]           # (P,)

    # ---- DOMINANT potency term: Hill dose-response, h=1 ----
    hill = dose / (dose + ic50)                            # (D,)

    # ---- real transport biology, then COMPRESS to a narrow [dmin,1] band ----
    # leaky-vessel access: BBB-penetrant OR leaky neovasculature lets drug in.
    vmin, vmax = vasc.min(), vasc.max()
    leak = (vasc - vmin) / (vmax - vmin + 1e-12)           # (P,) in [0,1]
    access = bbb[None, :] + (1.0 - bbb[None, :]) * leak[:, None]   # (P,D)
    # mild ECM size-sieve (thin stroma -> small k_ecm)
    mw_norm = (mw - mw.min()) / (mw.max() - mw.min())      # (D,)
    ecm_sieve = np.exp(-k_ecm * ecm[:, None] * mw_norm[None, :])   # (P,D)
    raw = access * ecm_sieve                               # (P,D), full-swing transport
    raw_norm = (raw - raw.min()) / (raw.max() - raw.min() + 1e-12)
    delivery = dmin + (1.0 - dmin) * raw_norm              # (P,D) in [dmin,1]

    # ---- patient-level, drug-agnostic modulation (does not reorder drugs) ----
    resist = 1.0 - 0.3 * hypoxic                           # mild hypoxic resistance
    patient = sens * resist                                # (P,)

    kill = patient[:, None] * hill[None, :] * delivery     # (P,D)

    rng = np.random.default_rng(seed)
    kill = kill * np.exp(rng.normal(0.0, noise_std, size=kill.shape))
    return np.clip(kill, 0.0, 1.0)


def main():
    cache = sl.load_cache()

    # Iterate on the delivery-compression band dmin: higher dmin => delivery is a
    # weaker multiplier => more strongly potency-dominated. We search for the
    # setting that makes the potency shortcut win most cleanly (highest pot acc,
    # sim not beating potency), and report the boundary.
    best = None
    sweep = []
    for dmin in [0.60, 0.70, 0.80, 0.85, 0.90, 0.95]:
        gt = build_gt_kill(cache, dmin=dmin)
        r = sl.score_policies(gt, seed=0)
        row = {
            "dmin": dmin,
            "acc_sim": r["acc"]["sim"],
            "acc_pot": r["acc"]["potency"],
            "acc_ml": r["acc"]["feature_ml"],
            "sim_beats_potency": r["sim_beats_potency"],
            "frac_not_most_potent": r["diag"]["frac_best_not_most_potent"],
        }
        sweep.append(row)
        # pick the run where potency wins hardest over sim while staying plausible
        score = (r["acc"]["potency"] - r["acc"]["sim"])
        if best is None or score > best[0]:
            best = (score, dmin, r)

    _, best_dmin, r = best
    out = {
        "chosen_dmin": best_dmin,
        "sweep": sweep,
        "result": r,
    }
    print("SWEEP_JSON " + json.dumps(sweep, default=float))
    print("RESULT_JSON " + json.dumps(r, default=float))
    print("CHOSEN_DMIN " + str(best_dmin))
    return out


if __name__ == "__main__":
    main()
