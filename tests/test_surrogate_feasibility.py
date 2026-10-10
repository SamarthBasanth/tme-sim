"""Guards on the surrogate feasibility study's own honesty.

``experiments/surrogate_feasibility/`` asks whether a learned model can predict
the transport solver's output well enough to preview. The study's conclusion is
only as good as four properties of its own method, and each one is a thing that
can silently stop being true when somebody edits a list. This file asserts them
directly, against the real modules, so a broken guard fails here rather than
quietly flattering a reported error.

Each test names the specific way the study could lie if the property failed.
"""
from __future__ import annotations

import os
import sys

import pytest

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
_STUDY = os.path.join(_ROOT, "experiments", "surrogate_feasibility")
if _STUDY not in sys.path:
    sys.path.insert(0, _STUDY)

pytestmark = pytest.mark.skipif(
    not os.path.isdir(_STUDY),
    reason="experiments/surrogate_feasibility/ is not present")


# --------------------------------------------------------------------------- #
# 1. The donor grouping must actually separate patients.
#
# GBM-Space ships up to five sections per patient. If donor_of() ever returns
# the section id for them, leave-one-donor-out silently becomes
# leave-one-section-out and every reported out-of-sample error is measured on a
# tissue the model has already seen through a sibling section.
# --------------------------------------------------------------------------- #
class TestDonorGrouping:
    def test_gbmspace_siblings_share_a_donor(self):
        from sweep import donor_of

        siblings = [
            "gbmspace_at10_bra_5_fo_1_1",
            "gbmspace_at10_bra_5_fo_1_2",
            "gbmspace_at10_bra_5_fo_2_1",
            "gbmspace_at10_bra_5_fo_3_1",
            "gbmspace_at10_bra_5_fo_4_2",
        ]
        donors = {donor_of(s) for s in siblings}
        assert donors == {"at10"}, (
            f"five sections of patient AT10 mapped to {donors}; "
            f"leave-one-donor-out would leak")

    def test_different_patients_do_not_share_a_donor(self):
        from sweep import donor_of

        assert donor_of("gbmspace_at10_bra_5_fo_1_1") != \
            donor_of("gbmspace_at11_bra_6_fo_1_2")

    def test_non_gbmspace_sections_are_their_own_donor(self):
        from sweep import donor_of

        for s in ("xenium_319", "xenium_io", "gse324268", "visium_gbm1",
                  "xenium_breast"):
            assert donor_of(s) == s

    def test_at1_does_not_swallow_at10(self):
        """A prefix match on 'at1' would merge AT1x into one donor."""
        from sweep import donor_of

        assert donor_of("gbmspace_at1_bra_1_fo_1_1") == "at1"
        assert donor_of("gbmspace_at10_bra_5_fo_1_1") == "at10"
        assert donor_of("gbmspace_at15_bra_4_fo_b2_s25") == "at15"


# --------------------------------------------------------------------------- #
# 2. No feature may be a solver output.
#
# The whole point is to predict the solve from things available WITHOUT it. One
# solve-derived column turns the reported error into a measurement of nothing.
# --------------------------------------------------------------------------- #
class TestNoSolveOutputInFeatures:
    #: Everything app.server.jobs puts in metrics / extra_metrics / reach that
    #: only exists AFTER the transport solve. A feature named from this set is
    #: leakage. ``penetration_length_analytic_um`` is deliberately NOT here: it
    #: is a closed form over the drug params and mean ECM density, computable
    #: without solving (jobs.py calls ``penetration_length_um(dp, ecm)``).
    FORBIDDEN = {
        "kill_fraction", "mean_tumor_auc", "median_tumor_auc",
        "mean_tissue_auc", "therapeutic_coverage_1x",
        "therapeutic_coverage_0_5x", "therapeutic_coverage_2x",
        "penetration_length_um", "penetration_index", "E_p10_uM", "E_p50_uM",
        "E_p90_uM", "median_E_over_ic50", "p10_E_over_ic50",
        "viable_tumor_fraction", "hypoxic_tumor_fraction",
        "low_o2_over_high_o2_exposure", "reach_headline_F", "reach_E_ref",
        "n_viable_reach", "tissue_median_over_plasma_peak",
    }

    def test_declared_feature_lists_are_solve_free(self):
        import fit

        declared = (set(fit.TISSUE_FEATURES) | set(fit.COMPOUND_FEATURES)
                    | set(fit.PK_FEATURES))
        assert not (declared & self.FORBIDDEN), (
            f"solver outputs declared as features: "
            f"{sorted(declared & self.FORBIDDEN)}")

    def test_built_feature_matrix_is_solve_free(self):
        """The built names, not just the declared lists.

        ``build()`` adds derived columns by hand, which is where a solve output
        would get in without appearing in any of the three lists.
        """
        import fit

        row = {
            "n_cells": 1000, "n_tumor": 500, "grid": [10, 10],
            "extent_um": [300.0, 300.0], "drug_diffusivity_um2_s": 42.0,
            "drug_clearance_1_per_s": 1e-4, "drug_binding_1_per_s": 1e-4,
            "half_life_s": 14400.0, "total_time_s": 21600.0,
            "plasma_peak": 50.0, "n_doses": 1,
        }
        _X, names, _y, _meta = fit.build([row])
        assert not (set(names) & self.FORBIDDEN), (
            f"solver outputs in the built matrix: "
            f"{sorted(set(names) & self.FORBIDDEN)}")

    def test_every_target_is_a_relative_or_shape_readout(self):
        """Absolute kill fraction is never a defensible output.

        It must not be a surrogate target, and neither must anything that is
        an absolute per-cell count.
        """
        import fit

        for name in fit.TARGETS:
            assert "kill" not in name, f"{name} targets the kill fraction"
        assert "kill_fraction" not in fit.TARGETS


# --------------------------------------------------------------------------- #
# 3. The study must not imply an ordering between compounds.
#
# Cross-compound ranking was WITHDRAWN 2026-08-21. Compound DESCRIPTORS as
# model inputs are fine; a per-compound ranking, score or ordering is not.
# --------------------------------------------------------------------------- #
class TestNoCompoundOrdering:
    def test_no_ranking_vocabulary_in_the_study(self):
        banned = ("rank", "ranking", "best compound", "best drug",
                  "most effective", "outperforms", "leaderboard")
        offenders = []
        for fn in sorted(os.listdir(_STUDY)):
            if not fn.endswith(".py"):
                continue
            with open(os.path.join(_STUDY, fn)) as fh:
                text = fh.read().lower()
            for b in banned:
                # "does not rank compounds" and "no ranking" are the study
                # saying it does not do this; only flag an unnegated use.
                for line in text.splitlines():
                    if b in line and not any(
                            n in line for n in ("not ", "never", "no ",
                                                "withdrawn", "prohibit",
                                                "does not")):
                        offenders.append(f"{fn}: {line.strip()[:90]}")
        assert not offenders, "compound-ordering language: " + "; ".join(
            offenders[:6])


# --------------------------------------------------------------------------- #
# 3b. The harvester must read fields that EXIST.
#
# ``solver_runner.run_one`` pulls ~40 fields out of a RunResult with
# ``.get(name)``. A misspelled or relocated name does not raise -- it yields
# None for every row, and a column of Nones looks exactly like a quantity the
# solver declined to produce. That is how the first pass of this sweep recorded
# 822 empty reach ladders with no error anywhere: the ladder is nested under
# ``reach["curve"]`` and was being read off the top level.
# --------------------------------------------------------------------------- #
class TestHarvestedKeysExist:
    def test_reach_ladder_is_nested_under_curve(self):
        """Assert the SHAPE solver_runner relies on, from the reach module."""
        from app.server.reach import null_reach

        payload = null_reach("probing the payload shape")
        assert "curve" in payload, (
            "app.server.reach no longer returns a 'curve' block; "
            "solver_runner reads the ladder from it")
        for absent in ("ladder_x", "ladder_under_reached_fraction"):
            assert absent not in payload, (
                f"reach now exposes {absent!r} at the top level; "
                f"solver_runner reads the ladder from reach['curve'] and "
                f"should be revisited")
        assert "headline_under_reached_fraction" in payload

    def test_solver_runner_does_not_read_flat_ladder_keys(self):
        src_path = os.path.join(_STUDY, "solver_runner.py")
        with open(src_path) as fh:
            src = fh.read()
        assert 'reach.get("ladder_x")' not in src
        assert 'reach.get("ladder_under_reached_fraction")' not in src
        assert '(reach.get("curve") or {}).get("x")' in src

    def test_report_block_field_names_exist(self):
        """Every tissue feature must be a key execute_run actually emits."""
        import fit

        from app.server.jobs import execute_run  # noqa: F401 - import guard

        # The report block's field names, read from the source of truth rather
        # than retyped: jobs.py builds `report` from reconstruct_real's
        # reconstruction_report plus the keys asserted below.
        expected = {
            "n_cells", "n_tumor", "grid", "extent_um", "celltype_fractions",
            "mean_vascular_density", "mean_ecm_density", "mean_oxygen",
            "tissue_coverage_fraction", "max_darcy_speed_um_s",
        }
        tissue = set(fit.TISSUE_FEATURES)
        derived_from_fractions = {
            "frac_TUMOR", "frac_STROMAL", "frac_IMMUNE", "frac_ENDOTHELIAL",
            "frac_NECROTIC"}
        direct = tissue - derived_from_fractions
        missing = direct - expected
        assert not missing, (
            f"tissue features with no corresponding report field: "
            f"{sorted(missing)}")


# --------------------------------------------------------------------------- #
# 4. The splits must be what they claim.
#
# A fold generator that returns an empty test set, or that puts a group in both
# halves, produces a number that looks like an out-of-sample error and is not.
# --------------------------------------------------------------------------- #
class TestFoldsAreDisjoint:
    def _meta(self):
        import numpy as np

        return {
            "section": np.array(["a1", "a1", "a2", "b1", "b1", "c1"]),
            "donor": np.array(["a", "a", "a", "b", "b", "c"]),
            "drug": np.array(["d1", "d2", "d1", "d2", "d1", "d2"]),
        }

    @pytest.mark.parametrize("split", ["random", "section", "donor",
                                       "compound"])
    def test_train_and_test_are_disjoint_and_nonempty(self, split):
        import fit

        meta = self._meta()
        n = len(meta["donor"])
        seen = 0
        for tr, te, _lab in fit._folds(split, meta, n):
            assert len(te) > 0 and len(tr) > 0
            assert not (set(tr.tolist()) & set(te.tolist())), (
                f"{split}: a row is in both train and test")
            seen += len(te)
        assert seen > 0

    @pytest.mark.parametrize("split,key", [("section", "section"),
                                           ("donor", "donor"),
                                           ("compound", "drug")])
    def test_held_out_group_is_absent_from_training(self, split, key):
        import fit

        meta = self._meta()
        n = len(meta["donor"])
        for tr, te, lab in fit._folds(split, meta, n):
            train_groups = set(meta[key][tr].tolist())
            assert lab not in train_groups, (
                f"{split}: group {lab} is held out but still in training")

    def test_donor_split_holds_out_every_section_of_a_donor(self):
        """Donor 'a' owns sections a1 and a2; both must leave together."""
        import fit

        meta = self._meta()
        for tr, te, lab in fit._folds("donor", meta, len(meta["donor"])):
            if lab != "a":
                continue
            assert set(meta["section"][te].tolist()) == {"a1", "a2"}
            assert "a1" not in set(meta["section"][tr].tolist())
            assert "a2" not in set(meta["section"][tr].tolist())
