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
4. The SECOND flag -- ``compile_drug.build_profile(two_route=...)``, which
   landed with the same silent-no-op shape and, until the ``Flag``
   parameterisation, had no table and no wiring invariant of its own. Sections
   2 and 3 are parameterised over every flag rather than duplicated, which is
   also the proof that the generalisation is real and not a rename.

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
    EXCHANGE,
    EXEMPTIONS,
    FLAGS,
    TWO_ROUTE,
    Category,
    Declaration,
    Flag,
    check_flag_propagation,
    declaration_for,
    digest,
    instrument,
    scan_call_sites,
    sweep_repo,
    undeclared_direct_callers,
    undeclared_paths,
    undeclared_wiring,
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
#
# Parameterised over FLAGS. Every one of these was written for
# ``transvascular_exchange``; that they hold ``two_route`` to the same word
# without a second copy is the point of the ``Flag`` parameterisation.
# =========================================================================== #
_ALL_DECLARATIONS = [(f, d) for f in FLAGS for d in f.declarations]
_DECL_IDS = [f"{f.name}:{d.name}" for f, d in _ALL_DECLARATIONS]


@pytest.mark.parametrize("flag", FLAGS, ids=[f.name for f in FLAGS])
def test_every_declaration_has_a_written_reason(flag: Flag):
    """``NO_PDE`` in particular can never be an inferred default: it exists only
    where a person wrote down why. ``Declaration`` enforces this at
    construction, so this test also proves the enforcement is live."""
    assert flag.declarations, f"{flag.name}: the allow-list must not be empty"
    for d in flag.declarations:
        assert len(d.reason.strip()) >= 25, d.name
    for e in flag.exemptions:
        assert len(e.reason.strip()) >= 25, e.path
    # The flag itself carries the same obligation as everything it declares.
    assert len(flag.reason.strip()) >= 25, flag.name

    with pytest.raises(ValueError):
        Declaration(path="x.py", category=Category.NO_PDE, reason="")
    with pytest.raises(ValueError):
        Declaration(path="x.py", category=Category.NO_PDE, reason="unreachable")


def test_a_flag_cannot_be_added_without_a_reason_or_a_way_to_fail():
    """The two ways a new ``Flag`` would be born blind, both refused.

    A flag with no reason is an undocumented switch. A CARRIER-LESS flag with
    no precursors is worse: it can never report anything, because a call is an
    arrival and there are no bare calls to find, so the guard would be green by
    construction. That is precisely the state ``two_route`` was in.
    """
    with pytest.raises(ValueError):
        Flag(name="x", destination="d", carrier="c", reason="short")
    with pytest.raises(ValueError):
        Flag(name="x", destination="d", carrier=None, precursors=(),
             reason="a perfectly adequate written reason, but no precursors")


def test_unreachable_is_the_same_verdict_as_no_pde():
    """The general name and the PDE spelling must be ONE member, not two.

    Two names for one claim is how "unreachable" and "no transport solve" drift
    into meaning subtly different things -- the exact conflation this module
    exists to prevent, one level up.
    """
    assert Category.UNREACHABLE is Category.NO_PDE
    assert Category.UNREACHABLE.value == "NO_PDE"


@pytest.mark.parametrize("flag", FLAGS, ids=[f.name for f in FLAGS])
def test_no_pde_entries_are_explicitly_declared(flag: Flag):
    """There must be at least one unreachable entry per flag and it must name a
    real file that genuinely never reaches that flag's destination."""
    no_pde = [d for d in flag.declarations if d.category is Category.NO_PDE]
    assert no_pde, "category C must be an explicit allow-list, not an absence"
    for d in no_pde:
        assert (REPO / d.path).is_file(), d.path
        cs = scan_call_sites(REPO / d.path, d.symbol, flag=flag)
        assert not cs.direct_calls, f"{d.name} is declared NO_PDE but calls the destination"
        assert not cs.engine_calls, f"{d.name} is declared NO_PDE but calls an entry point"
        assert not cs.precursor_calls, (
            f"{d.name} is declared NO_PDE but computes the governed quantity "
            "out of the precursors -- that is INERT, which is a different claim")


@pytest.mark.parametrize("flag,decl", _ALL_DECLARATIONS, ids=_DECL_IDS)
def test_declared_path_lands_in_its_declared_category(flag: Flag, decl: Declaration):
    """Every known path, held to its word. Source-level, so it covers all of the
    named modules -- including the expensive ones -- without a reconstruction.

    ``vessel_profile.py`` is declared per-function precisely because a
    whole-file answer would average an engine-path function together with a
    direct-solver one and produce something false about both.
    ``molecular/ecm_mesh.py`` is the same shape for the two-route flag.
    """
    cs = scan_call_sites(REPO / decl.path, decl.symbol, flag=flag)
    assert cs.predicted is decl.category, (
        f"{flag.name}: {decl.name}: declared {decl.category.value} but the "
        f"source says {cs.predicted.value} (destination={cs.direct_calls}, "
        f"bare={cs.bare_direct_calls}, entry={cs.engine_calls}, "
        f"precursor={cs.precursor_calls}).\n"
        f"Declared reason: {decl.reason}"
    )


def test_ecm_mesh_is_split_not_averaged():
    """The two-route twin of ``test_vessel_profile_is_split_not_averaged``.

    ``tissue_hindrance`` composes the IMPERMEABLE-cell law and cannot see the
    flag; ``tissue_hindrance_two_route`` reaches it. One file, two categories,
    and a file-level verdict would be false about one of them.
    """
    path = REPO / "tmesim/molecular/ecm_mesh.py"
    assert scan_call_sites(
        path, "tissue_hindrance_two_route", flag=TWO_ROUTE).predicted is Category.PROPAGATED
    assert scan_call_sites(
        path, "tissue_hindrance", flag=TWO_ROUTE).predicted is Category.INERT


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
@pytest.mark.parametrize("flag", FLAGS, ids=[f.name for f in FLAGS])
def test_no_undeclared_paths(flag: Flag):
    """The part that does not depend on anyone remembering.

    A file that reaches the destination without the carrier -- or that computes
    the governed quantity out of the precursors and never reaches the
    destination at all -- is a code path on which the flag is silently inert
    and nobody has written it down.
    """
    undeclared = undeclared_paths(REPO, flag=flag)
    assert not undeclared, (
        f"the {flag.name} flag is a silent no-op in these files and they have "
        "no declaration:\n  " + "\n  ".join(undeclared) + "\n\n"
        f"Add a Declaration (with a reason) to the {flag.name} table in "
        "tmesim.validation.propagation, or route the computation through "
        + (" / ".join(flag.entry_points) or "the destination")
        + ". Widening the exemptions is not the fix."
    )


def test_the_exchange_alias_still_answers_the_original_question():
    """``undeclared_direct_callers`` is the name every existing failure message
    and brief uses. It must keep meaning exactly what it meant."""
    assert undeclared_direct_callers(REPO) == undeclared_paths(REPO, flag=EXCHANGE)


def _wired_paths(flag: Flag = EXCHANGE):
    return {
        path
        for path, cs in sweep_repo(REPO, flag=flag).items()
        if cs.arrivals
    }


#: One production wiring point per flag, named. A SECOND one means two places
#: decide whether the physics is active, which is how a flag ends up silently
#: on in one path and off in another.
_PRODUCTION = {
    "transvascular_exchange": {"tmesim/sim/engine.py"},
    "two_route": {"tmesim/molecular/compile_drug.py"},
}


@pytest.mark.parametrize("flag", FLAGS, ids=[f.name for f in FLAGS])
def test_there_is_exactly_one_production_wiring_point(flag: Flag):
    """The original invariant, unchanged and still strict -- now for both flags.

    If a SECOND PRODUCTION wiring point ever appears, the mental model behind
    this whole guard has changed and every verdict above needs revisiting.
    """
    production = wiring_points(production_only=True, flag=flag)
    assert production == _PRODUCTION[flag.name], production
    assert len(production) == 1
    assert _wired_paths(flag) >= production, "the production wiring point vanished"


def test_engine_is_the_only_place_that_wires_the_flag_up():
    """Kept under its original name: this is THE invariant of this module."""
    assert wiring_points(production_only=True) == {"tmesim/sim/engine.py"}


@pytest.mark.parametrize("flag", FLAGS, ids=[f.name for f in FLAGS])
def test_every_wiring_point_in_scope_is_declared(flag: Flag):
    """Validation experiments that measure an INERT path have to inject
    ``exchange_k`` by hand -- the flag alone would be a silent no-op there, so
    doing the obvious thing would yield a confident, meaningless "no change".
    They are legitimate, but each one is written down individually rather than
    exempted by directory, because a blanket rule would let an undeclared
    no-op back in exactly where nobody would look for it.

    ``two_route`` narrows the SCOPE of that rule to ``tmesim/`` and the
    argument is at ``Flag.wiring_scope``: it has no carrier to inject, so a
    script calling the law directly is not faking a production path.
    """
    undeclared = undeclared_wiring(REPO, flag=flag)
    assert not undeclared, (
        f"these files reach the {flag.name} destination inside "
        f"{flag.wiring_scope} but are not declared wiring points:\n  "
        + "\n  ".join(undeclared)
        + "\n\nAdd a WiringPoint (production=False for anything that is not "
        "the shipped path) with a written reason. Do NOT exempt a directory "
        "wholesale."
    )


@pytest.mark.parametrize("flag", FLAGS, ids=[f.name for f in FLAGS])
def test_declared_wiring_points_still_exist(flag: Flag):
    """A stale entry is a hole: it would silently permit a future file at that
    path to wire the flag up with nobody re-reading the reason."""
    missing = {w.path for w in flag.wiring_points} - _wired_paths(flag)
    assert not missing, f"{flag.name}: wiring points that no longer wire anything: {missing}"


@pytest.mark.parametrize("flag", FLAGS, ids=[f.name for f in FLAGS])
def test_no_declaration_or_exemption_is_stale(flag: Flag):
    """A table nobody can delete an entry from stops meaning anything.

    Every declared file must still exist, and every exemption must still be
    excusing something real -- an exemption that matches nothing is a standing
    permission for whatever is written at that path next.
    """
    for d in flag.declarations:
        assert (REPO / d.path).is_file(), f"{flag.name}: declared file is gone: {d.name}"
    swept = sweep_repo(REPO, flag=flag)
    for e in flag.exemptions:
        cs = swept.get(e.path)
        # An exemption excuses one of two things: an undeclarable call site
        # (bare or stranded), or reaching the destination inside the wiring
        # scope without being a WiringPoint. Either is enough to keep it alive;
        # neither means the file has moved on and the entry is now a standing
        # permission for whatever gets written there next.
        excusing = cs is not None and (
            cs.undeclarable or (cs.arrivals and flag.in_wiring_scope(e.path))
        )
        assert excusing, (
            f"{flag.name}: {e.path} no longer has anything to excuse -- delete "
            "the exemption rather than leaving a standing permission")


# =========================================================================== #
# 4. The second flag, dynamically
#
# The static sweep says a file's call sites are in order. These two run the
# thing and watch the boundary, which is the stronger instrument and the one
# that separates "the law was applied and did nothing" from "the law was never
# applied". Cost: four molecular compiles at the module's own reduced
# particle/step defaults, ~1.5 s total. No tissue, no PDE.
# =========================================================================== #
def test_two_route_propagates_through_build_profile():
    """``build_profile(two_route=...)`` really does reach the law, and the
    number really does move. Without this the static table could be perfectly
    self-consistent about a parameter that does nothing."""
    from tmesim.molecular import compile_drug
    from tmesim.molecular.library import get_molecule

    mol = get_molecule("doxorubicin")

    def run(flag: bool) -> float:
        return float(
            compile_drug.build_profile(mol, two_route=flag).effective_diffusivity
        )

    v = check_flag_propagation(run, label="compile_drug.build_profile", flag=TWO_ROUTE)

    assert v.category is Category.PROPAGATED, v
    assert v.solver_calls_on == 1, "the law ran exactly once with the flag ON"
    assert v.solver_calls_off == 0, "and not at all with it OFF"
    assert v.exchange_k_calls_on == 1, "a call to a carrier-less destination IS an arrival"
    assert v.results_differ and v.trustworthy


def test_two_route_stranded_path_is_inert_not_unreachable():
    """The failure this flag can actually have, and the distinction that matters.

    Computing an effective diffusivity the pre-promotion way -- Brownian
    micro-tortuosity times the Ogston-Renkin mesh gate -- produces a number that
    is bit-identical whatever ``two_route`` is set to, because the law is never
    reached. That is INERT. Reporting it as UNREACHABLE would say "there is no
    diffusivity here for the flag to change", which is false: there is one, and
    it is the falsified one.
    """
    from tmesim.molecular import brownian, ecm_mesh, stokes_einstein

    def run(flag: bool) -> float:
        # Exactly steps 2-3b of build_profile, and then a stop.
        free_D = stokes_einstein.free_diffusivity(0.539)
        sim = brownian.simulate_effective_diffusivity(
            free_D, 0.4, n_particles=60, steps=400, seed=0, particle_radius_nm=0.539
        )
        mesh = ecm_mesh.ogston_renkin_hindrance(0.539)
        return float(sim["effective_D"]) * float(mesh["effective_hindrance"])

    v = check_flag_propagation(run, label="pre-promotion D_eff chain", flag=TWO_ROUTE)

    assert v.category is Category.INERT, v
    assert v.solver_calls_on == 0, "the law genuinely was not called"
    assert v.category is not Category.UNREACHABLE
    assert not v.trustworthy
    assert "RE-IMPLEMENTED" in v.detail
    # And the raw symptom a human would see: the two runs are the same run.
    assert run(False) == run(True)


def test_a_carrier_less_flag_has_no_bare_calls_by_construction():
    """The polarity, asserted rather than left implicit.

    For ``exchange_k`` the defect is a call that does not carry the flag. For
    ``two_route`` no such call can exist, so a guard copied across unchanged
    would have been green forever. ``stranded_calls`` is the dual, and this
    pins that the two halves do not overlap.
    """
    for path, cs in sweep_repo(REPO, flag=TWO_ROUTE).items():
        assert cs.bare_direct_calls == (), path
    stranded = {p: cs.stranded_calls
                for p, cs in sweep_repo(REPO, flag=TWO_ROUTE).items()
                if cs.stranded_calls}
    assert stranded, "the stranded sweep found nothing; the flag is blind again"

    # ...and the exchange flag is the mirror image: bare calls exist, and it
    # declares no precursors, so nothing is ever stranded.
    assert EXCHANGE.precursors == ()
    bare = {p for p, cs in sweep_repo(REPO).items() if cs.bare_direct_calls}
    assert bare, "the bare sweep found nothing; the original flag is blind"
    assert not any(cs.stranded_calls for cs in sweep_repo(REPO).values())


def test_the_two_route_sweep_actually_reads_the_modules_it_claims_to():
    """A sweep over an empty file list passes vacuously. Pin the floor."""
    swept = sweep_repo(REPO, flag=TWO_ROUTE)
    for expected in ("tmesim/molecular/compile_drug.py",
                     "tmesim/molecular/ecm_mesh.py",
                     "tmesim/molecular/two_route.py",
                     "tmesim/molecular/brownian.py",
                     "tmesim/validation/calibration.py"):
        assert expected in swept, f"{expected} vanished from the two-route sweep"
    assert swept["tmesim/molecular/compile_drug.py"].arrivals, (
        "the production wiring point stopped reaching the law")


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
        "app/server/reach.py",                   # steady_state_concentration
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
