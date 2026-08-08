"""The guard against a config flag that never reaches the solver.

``PDEConfig.transvascular_exchange`` only means something if something reads it
and hands ``exchange_k`` to :func:`tmesim.physics.diffusion.simulate_drug_transport`.
Exactly one function does: :func:`tmesim.sim.engine.simulate_response`. On every
path that calls the solver directly the flag is a **silent no-op** -- the run
completes, the numbers are plausible, and they are BIT-IDENTICAL to the flag
being off, so "this physics change has no effect" is both false and
unfalsifiable from the output.

Five of eleven re-validations in one session would have silently reported
"nothing changed" for this reason. Each brief happened to carry a hand-written
warning. That does not scale, so this file replaces the hand-written warning
with a mechanical check.

Layout
------
1. Proof of correctness -- the checker classifies a real example of EACH
   category, including the two that both look like "nothing changed"
   (``INERT``, the bug) and the two that both look like "the flag is fine"
   (``PROPAGATED``, ``REACHED_UNCHANGED``).
2. The declared table -- every known path lands in the category it was
   declared in, and ``NO_PDE`` entries carry a written reason.
3. The sweep -- no file anywhere in the repo calls the solver bare without
   being declared, so a NEW module cannot join the failure mode unnoticed.

Cost
----
Everything here runs on the small synthetic builders from ``conftest.py``
(``flat_micro_factory``, ``radial_micro_factory``, ``inert_drug``) at a 1200 s
window on a <=40x40 grid, or on closed-form analytic functions. No real-tissue
reconstruction is performed, and sections 2 and 3 parse source without executing
it. The whole file is well under a second.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

from tmesim.core import DrugParams
from tmesim.core.config import PDEConfig, SimConfig
from tmesim.sim.engine import simulate_response
from tmesim.validation.propagation import (
    DECLARATIONS,
    EXEMPTIONS,
    WIRING_POINTS,
    Category,
    Declaration,
    check_flag_propagation,
    declaration_for,
    digest,
    instrument,
    scan_call_sites,
    sweep_repo,
    undeclared_direct_callers,
    wiring_points,
)

REPO = Path(__file__).resolve().parents[1]

# ``experiments/`` is a script directory, not a package; one of the modules
# under proof lives there.
if str(REPO / "experiments") not in sys.path:
    sys.path.insert(0, str(REPO / "experiments"))

#: Short window, coarse grid, no convection: enough steps for the exchange sink
#: to bite, few enough to run in milliseconds.
PDE_KW = dict(dt=20.0, total_time=1200.0, convection=False, save_every=1000)

#: A perfusable probe drug. ``bbb_penetrance`` must be non-zero or ``PS`` is
#: zero everywhere and even the engine path would show no effect -- which is a
#: real category (see the avascular test below), just not the one being proved.
PROBE_DRUG = DrugParams(
    name="propagation-probe",
    diffusivity=300.0,
    clearance=1e-4,
    binding=1e-4,
    bbb_penetrance=0.5,
)


@pytest.fixture
def vascular_micro(radial_micro_factory):
    """Small perfused tissue: ``PS`` is non-zero, so exchange can matter."""
    micro, _, _ = radial_micro_factory(ny=40, nx=40, dx=15.0, n_tumor=300, seed=0)
    return micro


# =========================================================================== #
# 1. Proof of correctness: one real example of each category
# =========================================================================== #
def test_category_A_engine_path_propagates(vascular_micro):
    """(A) Through ``simulate_response`` the flag reaches the solver and the
    numbers move. This is the only path in the repo that wires it up."""

    def run(flag: bool) -> float:
        cfg = SimConfig(
            drug=PROBE_DRUG,
            pde=PDEConfig(transvascular_exchange=flag, **PDE_KW),
            seed=0,
        )
        return float(np.asarray(simulate_response(vascular_micro, cfg).transport.auc).sum())

    v = check_flag_propagation(run, label="sim.engine.simulate_response")

    assert v.category is Category.PROPAGATED, v
    assert v.exchange_k_calls_on == v.solver_calls_on == 1
    assert v.results_differ and v.trustworthy
    # The flag being off must NOT hand the solver an exchange field.
    assert v.exchange_coefficient_calls_on == 1


def test_category_B_direct_solver_call_is_inert(vascular_micro):
    """(B) ``experiments/dose_schedule.py`` -- the module the failure was
    MEASURED on. ``unit_transport`` takes the ``PDEConfig`` and calls the solver
    directly, so the flag never crosses the boundary.

    The checker must report INERT, and must do so from the STRUCTURE (exchange_k
    was never passed) rather than from the identical numbers, because identical
    numbers are a symptom that other things could also produce.
    """
    import dose_schedule

    def run(flag: bool) -> float:
        pde = PDEConfig(transvascular_exchange=flag, **PDE_KW)
        tr = dose_schedule.unit_transport(vascular_micro, PROBE_DRUG, pde, None)
        return float(np.asarray(tr.auc).sum())

    v = check_flag_propagation(run, label="experiments/dose_schedule.unit_transport")

    assert v.category is Category.INERT, v
    assert v.solver_calls_on == 1, "the solver DID run -- this is not NO_PDE"
    assert v.exchange_k_calls_on == 0, "and it was never handed exchange_k"
    assert not v.results_differ
    assert not v.trustworthy, "an on-vs-off comparison here compares a run to itself"

    # The failure in its raw form: the two runs are bit-identical.
    assert run(False) == run(True)
    assert digest(run(False)) == digest(run(True))


@pytest.mark.parametrize(
    "module_name,call",
    [
        # Slab delivery efficiency + the profile CV it is defined against.
        (
            "tmesim/validation/gist_imatinib.py",
            lambda m: (m.delivery_efficiency(50.0, 60.0), m.slab_profile_cv(50.0, 60.0, n=201)),
        ),
        # Descriptor-derived BBB penetrance: arithmetic over a molecule record.
        (
            "tmesim/validation/glioma_regional_pk.py",
            lambda m: m.molecular_profile("temozolomide"),
        ),
    ],
    ids=["gist_imatinib", "glioma_regional_pk"],
)
def test_category_C_analytic_modules_run_no_pde(module_name, call):
    """(C) Analytic closed forms. The solver is never called, so the flag is
    legitimately unreachable -- and that is a claim about this module's SCOPE,
    not evidence that the exchange term is negligible.

    NO_PDE is only an acceptable outcome because these two modules are declared
    as such in ``DECLARATIONS`` with a written reason (asserted below).
    """
    import importlib

    dotted = module_name[: -len(".py")].replace("/", ".")
    mod = importlib.import_module(dotted)

    def run(flag: bool):
        # A re-validation would set the flag exactly like this and then be
        # baffled that nothing moved. Nothing moves because nothing solves.
        _ = PDEConfig(transvascular_exchange=flag, **PDE_KW)
        return call(mod)

    v = check_flag_propagation(run, label=dotted)

    assert v.category is Category.NO_PDE, v
    assert v.solver_calls_off == v.solver_calls_on == 0
    assert not v.opaque, "digest could not fully see the result; verdict is weakened"

    d = declaration_for(module_name)
    assert d.category is Category.NO_PDE
    assert len(d.reason.strip()) >= 25


def test_inert_and_no_pde_are_distinguished(vascular_micro):
    """The whole point. Both produce "nothing changed"; only one is honest.

    A checker that cannot separate these is not solving the problem, so this
    asserts the separation directly rather than leaving it implicit in two
    tests that happen to pass.
    """
    import dose_schedule
    from tmesim.validation import glioma_regional_pk

    def inert(flag: bool) -> float:
        pde = PDEConfig(transvascular_exchange=flag, **PDE_KW)
        return float(
            np.asarray(dose_schedule.unit_transport(vascular_micro, PROBE_DRUG, pde, None).auc).sum()
        )

    def no_pde(flag: bool):
        _ = PDEConfig(transvascular_exchange=flag, **PDE_KW)
        return glioma_regional_pk.molecular_profile("temozolomide")

    v_inert = check_flag_propagation(inert, label="inert")
    v_nopde = check_flag_propagation(no_pde, label="no_pde")

    # Indistinguishable on the evidence a human would look at ...
    assert not v_inert.results_differ and not v_nopde.results_differ
    # ... and cleanly separated by the checker.
    assert v_inert.category is not v_nopde.category
    assert v_inert.solver_calls_on > 0 and v_nopde.solver_calls_on == 0


def test_reached_unchanged_is_not_no_pde(flat_micro_factory):
    """"Robust" and "unreachable" are different claims, and the checker must not
    collapse them either.

    An avascular tissue has ``vascular_density == 0``, so ``PS`` is zero
    everywhere: the engine builds a valid all-zero ``exchange_k``, hands it to
    the solver, and the answer is bit-identical. The flag was genuinely applied
    and genuinely did nothing. That is REACHED_UNCHANGED -- real evidence about
    the physics -- and it must not be reported as NO_PDE, which is no evidence
    at all.
    """
    micro = flat_micro_factory(ny=30, nx=30, dx=15.0)  # vascular_density = 0

    def run(flag: bool) -> float:
        cfg = SimConfig(
            drug=PROBE_DRUG,
            pde=PDEConfig(transvascular_exchange=flag, **PDE_KW),
            seed=0,
        )
        return float(np.asarray(simulate_response(micro, cfg).transport.auc).sum())

    v = check_flag_propagation(run, label="engine on avascular tissue")

    assert v.category is Category.REACHED_UNCHANGED, v
    assert v.exchange_k_calls_on == 1, "the flag DID reach the solver"
    assert not v.results_differ
    assert v.trustworthy, "'no effect' is an honest finding here, unlike INERT"
    assert v.category is not Category.NO_PDE


def test_reimplemented_solver_is_inert_not_no_pde(vascular_micro, inert_drug):
    """The worst blind spot, closed.

    ``molecular/immune_sink.py::transport_with_immune_sink`` integrates the same
    PDE with its own time loop (it imports ``physics.diffusion``'s primitives
    but never calls ``simulate_drug_transport``) and takes no ``exchange_k``
    argument at all. Watching only the shipped solver would report this path as
    NO_PDE -- "the flag is legitimately unreachable" -- when in fact a PDE WAS
    integrated and the flag could not reach it. That is INERT, and calling it
    NO_PDE would be exactly the conflation this module exists to prevent.
    """
    from tmesim.molecular import immune_sink
    from tmesim.treatment.drugs import delivery_source

    src = delivery_source(vascular_micro, inert_drug, 1.0)

    def run(flag: bool) -> float:
        pde = PDEConfig(transvascular_exchange=flag, **PDE_KW)
        # Attribute access, as a module-level ``from ... import`` would compile
        # to. A name bound into a LOCAL before the instrumented block starts is
        # not patchable and would escape the probe -- see the limitation noted
        # in ``instrument``.
        tr = immune_sink.transport_with_immune_sink(
            vascular_micro, inert_drug, src, pde
        )
        return float(np.asarray(tr.auc).sum())

    v = check_flag_propagation(run, label="transport_with_immune_sink")

    assert v.solver_calls_on == 0, "the shipped solver genuinely was not called"
    assert v.category is Category.INERT, v
    assert v.category is not Category.NO_PDE
    assert not v.trustworthy
    assert "RE-IMPLEMENTED" in v.detail


def test_nondeterministic_path_is_refused_not_guessed():
    """A stochastic path would show an off-vs-on difference that has nothing to
    do with the flag. Reporting that as PROPAGATED would be the same class of
    false conclusion this module exists to stop, so it is refused instead."""

    def run(flag: bool) -> float:
        return float(np.random.default_rng().random())

    v = check_flag_propagation(run, label="unseeded")
    assert v.category is Category.NONDETERMINISTIC, v
    assert not v.trustworthy


# =========================================================================== #
# 1b. The instrumentation itself
# =========================================================================== #
def test_instrument_sees_direct_callers_and_restores(vascular_micro):
    """A direct caller binds the solver by ``from ... import``, so patching only
    the defining module would miss precisely the calls that matter."""
    import dose_schedule
    from tmesim.physics import diffusion

    before = dose_schedule.simulate_drug_transport
    assert before is diffusion.simulate_drug_transport

    pde = PDEConfig(**PDE_KW)
    with instrument() as probe:
        dose_schedule.unit_transport(vascular_micro, PROBE_DRUG, pde, None)
        assert dose_schedule.simulate_drug_transport is not before, "not patched"
    assert probe.solver_calls == 1
    assert probe.exchange_k_calls == 0
    assert probe.of("delivery_source"), "delivery_source was not observed"

    # Fully restored, in both the defining module and the direct caller.
    assert dose_schedule.simulate_drug_transport is before
    assert diffusion.simulate_drug_transport is before


def test_instrument_restores_on_exception():
    from tmesim.physics import diffusion

    before = diffusion.simulate_drug_transport
    with pytest.raises(RuntimeError):
        with instrument():
            raise RuntimeError("boom")
    assert diffusion.simulate_drug_transport is before


def test_digest_is_bit_level_not_approximate():
    """The failure mode produces bit-identical output, so the comparison must be
    bit-level; an ``allclose`` digest would be blind to the distinction it is
    being asked about."""
    a = np.array([1.0, 2.0, 3.0])
    assert digest(a) == digest(a.copy())
    assert digest(a) != digest(a + 1e-16 * 4)  # one ulp at 3.0
    assert digest(0.0) != digest(-0.0)
    assert digest({"x": 1}) != digest({"x": 1.0})
    assert digest([1, 2]) != digest((1, 2))


# =========================================================================== #
# 2. The declared table
# =========================================================================== #
def test_every_declaration_has_a_written_reason():
    """``NO_PDE`` in particular can never be an inferred default: it exists only
    where a person wrote down why. ``Declaration`` enforces this at
    construction, so this test also proves the enforcement is live."""
    assert DECLARATIONS, "the allow-list must not be empty"
    for d in DECLARATIONS:
        assert len(d.reason.strip()) >= 25, d.name
    for e in EXEMPTIONS:
        assert len(e.reason.strip()) >= 25, e.path

    with pytest.raises(ValueError):
        Declaration(path="x.py", category=Category.NO_PDE, reason="")
    with pytest.raises(ValueError):
        Declaration(path="x.py", category=Category.NO_PDE, reason="unreachable")


def test_no_pde_entries_are_explicitly_declared():
    """There must be at least one NO_PDE entry and it must name a real file that
    genuinely contains no solver call site."""
    no_pde = [d for d in DECLARATIONS if d.category is Category.NO_PDE]
    assert no_pde, "category C must be an explicit allow-list, not an absence"
    for d in no_pde:
        assert (REPO / d.path).is_file(), d.path
        cs = scan_call_sites(REPO / d.path, d.symbol)
        assert not cs.direct_calls, f"{d.name} is declared NO_PDE but calls the solver"
        assert not cs.engine_calls, f"{d.name} is declared NO_PDE but calls the engine"


@pytest.mark.parametrize("decl", DECLARATIONS, ids=[d.name for d in DECLARATIONS])
def test_declared_path_lands_in_its_declared_category(decl: Declaration):
    """Every known path, held to its word. Source-level, so it covers all of the
    named modules -- including the expensive ones -- without a reconstruction.

    ``vessel_profile.py`` is declared per-function precisely because a
    whole-file answer would average an engine-path function together with a
    direct-solver one and produce something false about both.
    """
    cs = scan_call_sites(REPO / decl.path, decl.symbol)
    assert cs.predicted is decl.category, (
        f"{decl.name}: declared {decl.category.value} but the source says "
        f"{cs.predicted.value} (direct={cs.direct_calls}, "
        f"bare={cs.bare_direct_calls}, engine={cs.engine_calls}).\n"
        f"Declared reason: {decl.reason}"
    )


def test_vessel_profile_is_split_not_averaged():
    """The mixed file, called out on its own: one function is safe and one is
    not, and the guard must say so rather than issuing a file-level verdict."""
    path = REPO / "tmesim/validation/vessel_profile.py"
    assert scan_call_sites(path, "analyze_drug").predicted is Category.PROPAGATED
    assert scan_call_sites(path, "isolated_vessel_profile").predicted is Category.INERT
    # A whole-file glance would be dominated by the unsafe half, which is the
    # safe direction to fail in but is not a description of either function.
    assert scan_call_sites(path).predicted is Category.INERT


# =========================================================================== #
# 3. The sweep: a new module cannot join the failure mode unnoticed
# =========================================================================== #
def test_no_undeclared_direct_solver_callers():
    """The part that does not depend on anyone remembering.

    Any file that calls ``simulate_drug_transport`` without ``exchange_k`` and
    is neither declared nor exempted is a code path on which a ``PDEConfig``
    physics flag is silently inert and nobody has written it down.
    """
    undeclared = undeclared_direct_callers(REPO)
    assert not undeclared, (
        "these files call simulate_drug_transport without exchange_k and have "
        "no propagation declaration, so PDEConfig physics flags are silent "
        "no-ops there:\n  " + "\n  ".join(undeclared) + "\n\n"
        "Add a Declaration (with a reason) to "
        "tmesim.validation.propagation.DECLARATIONS, or pass "
        "exchange_coefficient(micro, drug) as exchange_k=. Widening EXEMPTIONS "
        "is not the fix."
    )


def _wired_paths():
    return {
        path
        for path, cs in sweep_repo(REPO).items()
        if cs.direct_calls_with_exchange_k
    }


def test_engine_is_the_only_place_that_wires_the_flag_up():
    """The original invariant, unchanged and still strict.

    If a SECOND PRODUCTION wiring point ever appears, the mental model behind
    this whole guard has changed and every verdict above needs revisiting --
    two places deciding whether the delivery law is active is how a flag ends
    up on in one path and off in another.
    """
    production = wiring_points(production_only=True)
    assert production == {"tmesim/sim/engine.py"}, production
    assert _wired_paths() >= production, "the production wiring point vanished"


def test_every_wiring_point_is_declared():
    """Validation experiments that measure an INERT path have to inject
    ``exchange_k`` by hand -- the flag alone would be a silent no-op there, so
    doing the obvious thing would yield a confident, meaningless "no change".
    They are legitimate, but each one is written down individually rather than
    exempted by directory, because a blanket rule would let an undeclared
    no-op back in exactly where nobody would look for it.
    """
    undeclared = _wired_paths() - wiring_points()
    assert not undeclared, (
        "these files pass exchange_k to the solver but are not in "
        "tmesim.validation.propagation.WIRING_POINTS:\n  "
        + "\n  ".join(sorted(undeclared))
        + "\n\nAdd a WiringPoint (production=False for an experiment) with a "
        "written reason. Do NOT exempt experiments/ wholesale."
    )


def test_declared_wiring_points_still_exist():
    """A stale entry is a hole: it would silently permit a future file at that
    path to wire the flag up with nobody re-reading the reason."""
    missing = {w.path for w in WIRING_POINTS} - _wired_paths()
    assert not missing, f"declared wiring points that no longer wire anything: {missing}"


def test_no_unaccounted_reimplemented_solver():
    """A hand-rolled time loop calls nothing the checker watches, so without
    this sweep a new one would silently register as NO_PDE -- "the flag is
    unreachable" -- when a PDE was in fact solved without it.

    Every file that borrows ``physics.diffusion``'s private numerics must be
    either watched (its solver is in ``ALT_SOLVERS``) or already declared for
    calling the shipped solver bare.
    """
    from tmesim.validation.propagation import ALT_SOLVERS, reimplementation_candidates

    # Where the watched alternate solvers live.
    alt_homes = {
        "tmesim/molecular/immune_sink.py",       # transport_with_immune_sink
        "tmesim/validation/moncayo_explant.py",  # bath_diffusion
    }
    declared = {d.path for d in DECLARATIONS}

    candidates = reimplementation_candidates(REPO)
    assert candidates, "the primitive sweep found nothing; it has gone blind"

    unaccounted = set(candidates) - alt_homes - declared - {e.path for e in EXEMPTIONS}
    assert not unaccounted, (
        "these files build their own transport loop from physics.diffusion's "
        "private numerics, so a PDEConfig flag cannot reach them and the "
        "checker would call them NO_PDE:\n  "
        + "\n  ".join(sorted(unaccounted))
        + "\n\nAdd the solver to ALT_SOLVERS/WATCHED, or declare the file."
    )

    # The watched alternate solvers still exist and still take no exchange_k --
    # if one grows the parameter it is no longer a blind spot and should move.
    import inspect

    from tmesim.molecular.immune_sink import transport_with_immune_sink
    from tmesim.validation.moncayo_explant import bath_diffusion

    for fn in (transport_with_immune_sink, bath_diffusion):
        assert fn.__name__ in ALT_SOLVERS
        assert "exchange_k" not in inspect.signature(fn).parameters, (
            f"{fn.__name__} now takes exchange_k; it is no longer an "
            "unreachable solver and ALT_SOLVERS is out of date"
        )


def test_sweep_covers_every_known_affected_module():
    """The seven paths that were actually caught out, still visible to the
    sweep. Guards against the scanner silently going blind."""
    known_affected = {
        "tmesim/validation/pk_reality.py",
        "tmesim/validation/kill_calibration.py",
        "tmesim/validation/immune_sink_impact.py",
        "tmesim/validation/vessel_profile.py",
        "tmesim/validation/numerical_rigor.py",
        "experiments/regimen.py",
        "experiments/dose_schedule.py",
    }
    swept = sweep_repo(REPO)
    missing = known_affected - set(swept)
    assert not missing, missing

    #: ``vessel_profile.py`` is declared per-function, so a file-level lookup
    #: MUST fail for it -- the two halves are in different categories and there
    #: is no honest whole-file answer. Every other affected file is declared
    #: whole and INERT.
    split = {"tmesim/validation/vessel_profile.py": "isolated_vessel_profile"}

    for path in known_affected:
        assert swept[path].bare_direct_calls, path
        symbol = split.get(path)
        if symbol is None:
            assert declaration_for(path).category is Category.INERT, path
        else:
            with pytest.raises(LookupError):
                declaration_for(path)
            assert declaration_for(path, symbol).category is Category.INERT
