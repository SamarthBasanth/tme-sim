"""Stage 0 contrast experiment: does spatial drug-transport modeling add value
OVER a non-spatial intrinsic-sensitivity baseline, and does that depend on the
dose regime?

Hypothesis: when the dose saturates kill (exposure >> IC50 everywhere), spatial
delivery differences are erased and the mechanistic sim can only recover what the
intrinsic-sensitivity baseline already has. When the dose is transport-LIMITED
(exposure near IC50, so poorly-perfused tumor core survives), the sim's transport
physics should become predictive and beat the baseline.

We sweep plasma_peak (dose) on ONE fixed synthetic cohort and tabulate, per
regime: mean kill (saturation indicator), sim ROC-AUC/Spearman, the best baseline,
and the sim-minus-best-baseline delta.
"""
import json
import numpy as np

from tmesim.core import SimConfig, PKConfig, PDEConfig
from tmesim.data.synthetic import generate_synthetic_cohort
from tmesim.validation.experiment import run_stage0_experiment

N = 40
TOTAL_TIME = 1200.0
DOSE_SWEEP = [4.0, 10.0, 25.0, 60.0, 140.0]   # plasma_peak uM; IC50(temozolomide-like)=25

print(f"Building fixed synthetic cohort n={N} (seed=0) ...", flush=True)
cohort = generate_synthetic_cohort(n=N, seed=0)
y_true = np.array([r.observed_response for r in cohort])
print(f"  observed_response: mean={y_true.mean():.3f} std={y_true.std():.3f} "
      f"responders(>median)={int((y_true> np.median(y_true)).sum())}/{N}", flush=True)

rows = []
for pp in DOSE_SWEEP:
    cfg = SimConfig(pk=PKConfig(plasma_peak=pp), pde=PDEConfig(total_time=TOTAL_TIME), seed=0)
    exp = run_stage0_experiment(cohort, cfg, out_dir=f"results/regime_pp{int(pp)}")
    sim = exp.model_scores
    # best baseline by roc_auc
    best_name = max(exp.baseline_scores, key=lambda k: exp.baseline_scores[k].get("roc_auc", 0))
    best = exp.baseline_scores[best_name]
    intr = exp.baseline_scores.get("intrinsic_sensitivity", {})
    mean_kill = float(np.mean([p.get("sim_pred", np.nan) for p in exp.per_patient]))
    delta = sim.get("roc_auc", 0) - best.get("roc_auc", 0)
    rows.append(dict(plasma_peak=pp, mean_kill=mean_kill,
                     sim_auc=sim.get("roc_auc"), sim_spear=sim.get("spearman"),
                     intrinsic_auc=intr.get("roc_auc"),
                     best_baseline=best_name, best_auc=best.get("roc_auc"),
                     delta_auc=delta, beat=exp.beat_baseline, verdict=exp.verdict))
    print(f"  pp={pp:6.1f}  mean_kill={mean_kill:.3f}  sim_auc={sim.get('roc_auc'):.3f}  "
          f"best_baseline={best_name}({best.get('roc_auc'):.3f})  delta={delta:+.3f}  "
          f"{'BEAT' if exp.beat_baseline else 'no'}", flush=True)

print("\n================ DOSE-REGIME CONTRAST ================")
hdr = f"{'plasma_peak':>11} {'mean_kill':>9} {'sim_AUC':>8} {'sim_rho':>8} {'intrinsic_AUC':>13} {'delta_AUC':>9} {'beat':>5}"
print(hdr); print("-"*len(hdr))
for r in rows:
    print(f"{r['plasma_peak']:>11.1f} {r['mean_kill']:>9.3f} {r['sim_auc']:>8.3f} "
          f"{r['sim_spear']:>8.3f} {r['intrinsic_auc']:>13.3f} {r['delta_auc']:>+9.3f} "
          f"{str(r['beat']):>5}")

best_regime = max(rows, key=lambda r: r["delta_auc"])
print(f"\nMost favorable regime for transport modeling: plasma_peak={best_regime['plasma_peak']} "
      f"-> sim_AUC={best_regime['sim_auc']:.3f} vs intrinsic={best_regime['intrinsic_auc']:.3f} "
      f"(delta {best_regime['delta_auc']:+.3f})")
with open("results/dose_regime_contrast.json", "w") as f:
    json.dump(rows, f, indent=2, default=float)
print("Saved results/dose_regime_contrast.json")
print("BEST_REGIME_PP", best_regime["plasma_peak"])
