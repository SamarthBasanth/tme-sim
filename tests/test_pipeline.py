"""End-to-end smoke + determinism tests for the simulation engine and the
Stage 0 validation experiment.

These exercise the *integrator* wave (``tmesim.sim.engine.simulate_response``)
and the *validation harness* (``tmesim.validation.experiment.run_stage0_experiment``)
against real small synthetic microenvironments.

NOTE: at the time of writing, ``tmesim/sim/engine.py`` and
``tmesim/validation/experiment.py`` are not yet implemented (only stub package
``__init__.py`` files exist). Each test therefore ``importorskip``s its target so
the suite stays green until those wave-2 modules land, at which point these tests
activate automatically and assert the real contract.
"""
from __future__ import annotations

import numpy as np
import pytest

from tmesim.core import (
    ExperimentResult,
    PDEConfig,
    SimConfig,
    SimulationResult,
)
from tmesim.data.synthetic import generate_synthetic_patient


def _fast_sim_config(seed: int = 0) -> SimConfig:
    """A SimConfig with a deliberately short PDE window for a fast smoke test."""
    cfg = SimConfig(seed=seed)
    cfg.pde = PDEConfig(dt=20.0, total_time=1200.0, convection=True, save_every=10)
    return cfg


# --------------------------------------------------------------------------- #
# simulate_response: single-patient smoke test
# --------------------------------------------------------------------------- #
def test_simulate_response_smoke():
    engine = pytest.importorskip(
        "tmesim.sim.engine",
        reason="tmesim.sim.engine.simulate_response not implemented yet (wave 2)",
    )
    rec = generate_synthetic_patient("SYN-SMOKE", seed=3, grid=(40, 40), dx=15.0)
    cfg = _fast_sim_config(seed=0)

    result = engine.simulate_response(rec.micro, cfg)

    assert isinstance(result, SimulationResult)
    # predicted_response is a probability.
    assert 0.0 <= float(result.predicted_response) <= 1.0
    # viable / kill fractions are complementary and in range.
    assert 0.0 <= float(result.viable_tumor_fraction) <= 1.0
    assert 0.0 <= float(result.kill_fraction) <= 1.0
    assert float(result.kill_fraction) == pytest.approx(
        1.0 - float(result.viable_tumor_fraction), abs=1e-6)
    # history + times populated and aligned.
    assert np.asarray(result.viable_history).ndim == 1
    assert len(result.viable_history) == len(result.times)
    assert len(result.times) >= 1
    # visualization fields populated (contract: oxygen/ifp/auc/vascular).
    assert isinstance(result.fields, dict) and len(result.fields) > 0
    for name, arr in result.fields.items():
        assert np.asarray(arr).shape == rec.micro.domain.shape, name


# --------------------------------------------------------------------------- #
# simulate_response: determinism
# --------------------------------------------------------------------------- #
def test_simulate_response_deterministic():
    engine = pytest.importorskip(
        "tmesim.sim.engine",
        reason="tmesim.sim.engine.simulate_response not implemented yet (wave 2)",
    )
    cfg = _fast_sim_config(seed=7)

    # Regenerate identical inputs each run (simulate_response mutates cells).
    rec_a = generate_synthetic_patient("SYN-DET", seed=5, grid=(40, 40), dx=15.0)
    rec_b = generate_synthetic_patient("SYN-DET", seed=5, grid=(40, 40), dx=15.0)

    r_a = engine.simulate_response(rec_a.micro, cfg)
    r_b = engine.simulate_response(rec_b.micro, cfg)

    assert float(r_a.predicted_response) == float(r_b.predicted_response)
    assert float(r_a.viable_tumor_fraction) == float(r_b.viable_tumor_fraction)


# --------------------------------------------------------------------------- #
# run_stage0_experiment: cohort-level smoke test
# --------------------------------------------------------------------------- #
def test_run_stage0_experiment_smoke(tmp_path, cohort_factory):
    experiment = pytest.importorskip(
        "tmesim.validation.experiment",
        reason="tmesim.validation.experiment.run_stage0_experiment not implemented "
               "yet (wave 2)",
    )
    records = cohort_factory(n=6, seed=0, grid=(36, 36))
    cfg = _fast_sim_config(seed=0)

    exp = experiment.run_stage0_experiment(records, cfg, out_dir=str(tmp_path))

    assert isinstance(exp, ExperimentResult)
    # Model scores present and finite.
    assert isinstance(exp.model_scores, dict) and exp.model_scores
    for metric, val in exp.model_scores.items():
        assert np.isfinite(val), f"model score {metric} is not finite: {val}"
    # Baselines scored.
    assert isinstance(exp.baseline_scores, dict) and exp.baseline_scores
    for bname, scores in exp.baseline_scores.items():
        assert isinstance(scores, dict) and scores
        for metric, val in scores.items():
            assert np.isfinite(val), f"baseline {bname}/{metric} not finite: {val}"
    # Primary metric is scored for the model.
    assert exp.primary_metric in exp.model_scores
    # Verdict is a non-empty human-readable string; beat_baseline is a bool.
    assert isinstance(exp.verdict, str) and exp.verdict.strip()
    assert isinstance(exp.beat_baseline, bool)
    # Per-patient rows, if present, align with the cohort.
    if exp.per_patient:
        assert len(exp.per_patient) == len(records)
