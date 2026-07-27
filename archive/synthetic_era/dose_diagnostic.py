"""Why is kill invariant to a 35x dose change? Trace dose -> tissue AUC ->
effective exposure E -> kill on ONE patient across a wide dose range, and check
whether we can reach a genuinely transport-limited regime (E near IC50, so the
poorly-perfused core survives) by lowering dose/time/penetrance.
"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np
from tmesim.core import SimConfig, PKConfig, PDEConfig, DrugParams, DRUG_LIBRARY
from tmesim.data.synthetic import generate_synthetic_patient
from tmesim.sim.engine import simulate_response

rec = generate_synthetic_patient("diag", seed=1)
drug = DRUG_LIBRARY["temozolomide-like"]
print(f"drug ic50={drug.ic50} uM, bbb={drug.bbb_penetrance}, clearance={drug.clearance}/s")
print(f"{'plasma_pk':>9} {'window_s':>8} {'mean_AUC':>10} {'meanE=AUC/T':>11} {'E/ic50':>7} {'kill':>6} {'pen_idx':>7}")
for pp, tt in [(140,1200),(25,1200),(4,1200),(1,1200),(0.2,600),(0.05,300),(0.01,300)]:
    cfg = SimConfig(drug=drug, pk=PKConfig(plasma_peak=pp), pde=PDEConfig(total_time=tt), seed=0)
    res = simulate_response(rec.micro, cfg)
    auc = res.transport.auc
    meanauc = float(auc[auc>0].mean()) if (auc>0).any() else 0.0
    meanE = meanauc/tt
    print(f"{pp:>9.2f} {tt:>8.0f} {meanauc:>10.1f} {meanE:>11.3f} {meanE/drug.ic50:>7.2f} "
          f"{res.kill_fraction:>6.3f} {res.transport.penetration_index:>7.3f}")

# a low-penetrance drug (doxorubicin-like, bbb=0.05) at low dose -> transport-limited?
print("\n-- doxorubicin-like (bbb=0.05, ic50=0.5), low dose, short window --")
dox = DRUG_LIBRARY["doxorubicin-like"]
for pp, tt in [(5,600),(1,600),(0.2,300)]:
    cfg = SimConfig(drug=dox, pk=PKConfig(plasma_peak=pp), pde=PDEConfig(total_time=tt), seed=0)
    res = simulate_response(rec.micro, cfg)
    auc = res.transport.auc; meanauc = float(auc[auc>0].mean()) if (auc>0).any() else 0.0
    print(f"  pp={pp:>5.2f} tt={tt:>4.0f}  meanE={meanauc/tt:>8.4f}  E/ic50={meanauc/tt/dox.ic50:>6.2f}  "
          f"kill={res.kill_fraction:>6.3f}  pen={res.transport.penetration_index:.3f}")
