"""Generate the canonical Stage 0 report + figures at a representative
(non-saturated) dose, so the drug-penetration figure is actually informative.
"""
import warnings; warnings.filterwarnings("ignore")
from tmesim.core import SimConfig, PKConfig, PDEConfig
from tmesim.data.synthetic import generate_synthetic_cohort
from tmesim.sim.engine import simulate_response
from tmesim.validation.experiment import run_stage0_experiment
from tmesim.viz.report import build_report

cfg = SimConfig(pk=PKConfig(plasma_peak=1.0), pde=PDEConfig(total_time=600.0), seed=0)
cohort = generate_synthetic_cohort(n=40, seed=0)
print("running experiment for report ...", flush=True)
exp = run_stage0_experiment(cohort, cfg, out_dir="results")
example = simulate_response(cohort[0].micro, cfg)
path = build_report(exp, example, cohort[0].micro, out_dir="results")
print("report:", path)
print("verdict:", exp.verdict)
