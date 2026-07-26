"""Ground-truth motif #11: ANTI-DELIVERY (adversarial) drug response in GBM.

Biological assumption
---------------------
This is a deliberately ADVERSARIAL motif for a delivery-aware PDE sim. In this
tumor biology the drug that works best is the one that is WORST at crossing the
blood-brain barrier -- an efflux / priming quirk.

Mechanism (efflux-saturation + size trapping): the GBM neurovascular unit is
armored with ABC efflux transporters (P-gp / BCRP). The drugs that look most
"deliverable" on paper (small, BBB-penetrant: temozolomide, lomustine,
cisplatin) are precisely the ones the pumps clear efficiently and continuously,
so their intratumoral exposure stays sub-lethal. By contrast, a large,
poorly-penetrant, strong-efflux-SUBSTRATE agent SATURATES the pumps at the dose
used; the fraction that does slip past the barrier is then trapped
intracellularly (high-MW ion trapping / slow efflux) and builds a sustained
lethal AUC. So paradoxically, LOW drug_bbb + HIGH drug_mw predicts the best
effective kill -- the opposite of a delivery-gated model.

There is a molecular-weight "sweet spot", not a monotone bigger-is-better: the
very largest molecules (doxorubicin) aggregate / bind serum protein and lose
some trapped-fraction advantage, while a mid-large agent (methotrexate) sits in
the retention optimum for most tumors. Denser ECM (high ecm) shifts the trapping
optimum toward larger molecules (a stiff matrix holds even the biggest agent),
so which poorly-deliverable drug wins varies patient-to-patient with ecm --
giving a non-degenerate best-drug distribution that an ML model on observable
features can partly learn.

Potency (Hill dose-response x patient sens) still modulates kill, but it is
deliberately COMPRESSED so the anti-delivery retention term drives the argmax --
faithful to the motif that delivery-quirk, not raw potency, decides the winner.

Why the sim should FAIL: a mechanistic delivery-aware sim rewards deliverable
(high-BBB, small) drugs; here those are the worst. The sim's ranking should
anti-correlate with this ground truth. Reported honestly.

Key parameters
--------------
  retention peak  = 0.62 + 0.34*ecm_scaled   (ECM shifts MW trapping optimum)
  retention width = 0.30                       (Gaussian MW-retention kernel)
  efflux pref     = (1 - drug_bbb)             (poor penetrants saturate/starve pumps)
  potency compression: eff = sens*(0.55 + 0.45*hill)
  dose=30 (cache), noise_std=0.02, seed=4711
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

    # normalized molecular weight in [0,1] across the panel
    mw_norm = (mw - mw.min()) / (mw.max() - mw.min())              # (D,)

    # --- ANTI-DELIVERY term 1: efflux-substrate preference ---
    # Poorly BBB-penetrant drugs are strong efflux substrates: at this dose they
    # saturate the pumps and starve competing clearance -> retention advantage.
    efflux_pref = 1.0 - bbb                                        # (D,) high for dox, mtx

    # --- ANTI-DELIVERY term 2: MW size-trapping with a per-patient sweet spot ---
    # Trapping peaks at an intermediate-large MW; denser ECM pushes the optimum
    # toward the largest molecules (stiff matrix holds even big agents).
    ecm_scaled = (ecm - ecm.min()) / (ecm.max() - ecm.min() + 1e-9)  # (P,) in [0,1]
    peak = 0.62 + 0.34 * ecm_scaled                                # (P,) MW-optimum
    width = 0.30
    retention = np.exp(-((mw_norm[None, :] - peak[:, None]) / width) ** 2)  # (P,D)

    # combined paradoxical delivery advantage: high ONLY for low-BBB, large-MW drugs
    anti_deliv = efflux_pref[None, :] * retention                 # (P,D)

    # --- POTENCY term (compressed so it does not dominate the argmax) ---
    hill = dose / (dose + ic50)                                   # (D,) in (0,1)
    eff_potency = sens[:, None] * (0.55 + 0.45 * hill[None, :])   # (P,D)

    kill = anti_deliv * eff_potency                              # (P,D)

    # normalize to a sensible [0,1] dynamic range before noise
    kill = kill / (kill.max() + 1e-9)

    rng = np.random.default_rng(4711)
    kill = kill + rng.normal(0.0, 0.02, size=kill.shape)
    kill = np.clip(kill, 0.0, 1.0)
    return kill


if __name__ == "__main__":
    gt = build_gt_kill()
    result = selection_lib.score_policies(gt, seed=0)
    selection_lib.print_json(result)
