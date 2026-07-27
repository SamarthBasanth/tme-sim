"""The decisive Stage 0 test, run in the TRANSPORT-LIMITED regime (E ~ IC50),
where poorly-perfused tumor core survives and spatial drug delivery — the thing
the mechanistic sim models and the intrinsic-sensitivity baseline cannot see —
should become predictive.

Same fixed cohort as the saturating sweep, so the only change is the dose.
"""
import warnings; warnings.filterwarnings("ignore")
import json
import numpy as np
from tmesim.core import SimConfig, PKConfig, PDEConfig
from tmesim.data.synthetic import generate_synthetic_cohort
from tmesim.validation.experiment import run_stage0_experiment

N = 40
TOTAL_TIME = 600.0
DOSE_SWEEP = [0.1, 0.3, 0.6, 1.0, 2.0]   # plasma_peak uM; IC50=25 -> E/ic50 ~ 0.05..3

cohort = generate_synthetic_cohort(n=N, seed=0)
y = np.array([r.observed_response for r in cohort])
print(f"cohort n={N}: observed mean={y.mean():.3f} std={y.std():.3f}\n")

hdr = f"{'plasma_pk':>9} {'mean_kill':>9} {'sim_AUC':>8} {'sim_rho':>8} {'intrinsic_AUC':>13} {'logreg_AUC':>10} {'delta':>7} {'beat':>5}"
print(hdr); print("-"*len(hdr))
rows = []
for pp in DOSE_SWEEP:
    cfg = SimConfig(pk=PKConfig(plasma_peak=pp), pde=PDEConfig(total_time=TOTAL_TIME), seed=0)
    exp = run_stage0_experiment(cohort, cfg, out_dir=f"results/tl_pp{pp}".replace(".","p"))
    sim = exp.model_scores
    bs = exp.baseline_scores
    best_name = max(bs, key=lambda k: bs[k].get("roc_auc", 0))
    delta = sim["roc_auc"] - bs[best_name]["roc_auc"]
    mean_kill = float(np.mean([p.get("sim_pred", np.nan) for p in exp.per_patient]))
    rows.append(dict(pp=pp, mean_kill=mean_kill, sim_auc=sim["roc_auc"], sim_rho=sim["spearman"],
                     intrinsic_auc=bs.get("intrinsic_sensitivity",{}).get("roc_auc"),
                     logreg_auc=bs.get("expression_logreg",{}).get("roc_auc"),
                     best=best_name, delta=delta, beat=bool(exp.beat_baseline), verdict=exp.verdict))
    print(f"{pp:>9.2f} {mean_kill:>9.3f} {sim['roc_auc']:>8.3f} {sim['spearman']:>8.3f} "
          f"{bs.get('intrinsic_sensitivity',{}).get('roc_auc',float('nan')):>13.3f} "
          f"{bs.get('expression_logreg',{}).get('roc_auc',float('nan')):>10.3f} {delta:>+7.3f} "
          f"{'BEAT' if exp.beat_baseline else 'no':>5}")

best = max(rows, key=lambda r: r["delta"])
print(f"\nBest regime for transport: pp={best['pp']} sim_AUC={best['sim_auc']:.3f} "
      f"vs intrinsic {best['intrinsic_auc']:.3f} (delta {best['delta']:+.3f}) beat={best['beat']}")
print("VERDICT:", best["verdict"])
with open("results/transport_limited.json","w") as f:
    json.dump(rows, f, indent=2, default=float)
print("saved results/transport_limited.json ; BEST_PP", best["pp"])
