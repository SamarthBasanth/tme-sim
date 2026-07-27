"""Precompute the EXPENSIVE, ground-truth-INDEPENDENT quantities for the
drug-selection study, so hundreds of ground-truth models can be swept cheaply.

Saves, for a fixed cohort and drug panel:
  - K_sim[p,d]        : mechanistic PDE-sim predicted kill fraction (the model)
  - potency[p,d]      : delivery-agnostic potency score sens*Hill(DOSE, IC50)
  - feat[p,:]         : OBSERVABLE patient features (vasc, ecm, sens, hypoxic, subtype)
  - vasc[p], ecm[p]   : tumor-region perfusion / ECM (drive ground-truth delivery)
  - drug meta (mw, ic50, bbb)
None of these depend on any ground-truth model, so the sweep never re-runs the sim.
"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np

from tmesim.core import SimConfig, PKConfig, PDEConfig
from tmesim.data.synthetic import generate_synthetic_cohort
from tmesim.molecular import library
from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
from tmesim.validation.calibration import calibrate_alpha, calibrated_drugparams
from tmesim.sim.engine import simulate_response

PANEL = ["temozolomide", "lomustine", "cisplatin", "methotrexate", "doxorubicin"]
N = 30
DOSE = 30.0
TOTAL_TIME = 600.0
SUBTYPES = ["proneural", "classical", "mesenchymal"]


def hill(e, ic50, h=1.3):
    e = np.maximum(e, 0.0)
    return e**h / (e**h + ic50**h)


def main():
    alpha = calibrate_alpha(target_um=45.0, drug="doxorubicin")
    dmeta = {}
    dparams = {}
    for name in PANEL:
        mol = library.get_molecule(name)
        dp = apply_barrier_to_drugparams(calibrated_drugparams(name, alpha=alpha), mol)
        dparams[name] = dp
        dmeta[name] = dict(mw=float(mol.get("mw")), ic50=float(dp.ic50), bbb=float(dp.bbb_penetrance))

    cohort = generate_synthetic_cohort(n=N, seed=0)
    P, D = N, len(PANEL)
    K_sim = np.zeros((P, D)); potency = np.zeros((P, D))
    vasc = np.zeros(P); ecm = np.zeros(P); sens = np.zeros(P)
    hypoxic = np.zeros(P); subtype_idx = np.zeros(P, dtype=int)

    print(f"precompute: {P} patients x {D} drugs = {P*D} sims", flush=True)
    for p, rec in enumerate(cohort):
        m = rec.micro
        tum = m.cells.tumor_mask()
        col, row = m.domain.to_index(m.cells.x[tum], m.cells.y[tum])
        vasc[p] = float(np.mean(m.fields.vascular_density[row, col]))
        ecm[p] = float(np.mean(m.fields.ecm_density[row, col]))
        sens[p] = float(np.mean(m.cells.sensitivity[tum]))
        hypoxic[p] = float(rec.features.get("hypoxic_fraction", 0.0))
        st = str(rec.micro.metadata.get("subtype", rec.features.get("subtype", "proneural")))
        subtype_idx[p] = SUBTYPES.index(st) if st in SUBTYPES else 0
        for d, name in enumerate(PANEL):
            cfg = SimConfig(drug=dparams[name], pk=PKConfig(plasma_peak=1.0),
                            pde=PDEConfig(total_time=TOTAL_TIME), seed=0)
            K_sim[p, d] = float(simulate_response(m, cfg).kill_fraction)
            potency[p, d] = float(sens[p] * hill(DOSE, dmeta[name]["ic50"]))
        if (p + 1) % 5 == 0:
            print(f"  {p+1}/{P} patients", flush=True)

    # observable feature matrix (what a cheap ML baseline is allowed to use)
    onehot = np.eye(len(SUBTYPES))[subtype_idx]
    feat = np.column_stack([vasc, ecm, sens, hypoxic, onehot])
    feat_names = ["vasc", "ecm", "sens", "hypoxic"] + [f"subtype_{s}" for s in SUBTYPES]

    drug_mw = np.array([dmeta[n]["mw"] for n in PANEL])
    drug_ic50 = np.array([dmeta[n]["ic50"] for n in PANEL])
    drug_bbb = np.array([dmeta[n]["bbb"] for n in PANEL])

    np.savez("results/selection_cache.npz",
             panel=np.array(PANEL), K_sim=K_sim, potency=potency, feat=feat,
             feat_names=np.array(feat_names), vasc=vasc, ecm=ecm, sens=sens,
             drug_mw=drug_mw, drug_ic50=drug_ic50, drug_bbb=drug_bbb, dose=DOSE, alpha=alpha)
    print("saved results/selection_cache.npz")
    print(f"vasc spread [{vasc.min():.3f},{vasc.max():.3f}] ecm [{ecm.min():.3f},{ecm.max():.3f}] "
          f"sens [{sens.min():.3f},{sens.max():.3f}]")
    print("sim kill matrix column means (per drug):",
          {PANEL[d]: round(float(K_sim[:, d].mean()), 3) for d in range(D)})


if __name__ == "__main__":
    main()
