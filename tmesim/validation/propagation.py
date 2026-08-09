"""Does a config flag actually REACH the solver on this code path?

Why this module exists
======================
``PDEConfig`` carries physics flags -- ``transvascular_exchange``,
``time_resolved`` -- that only have meaning once something reads them and wires
them into the solver. In this codebase exactly ONE thing does that:
:func:`tmesim.sim.engine.simulate_response`. It turns
``cfg.pde.transvascular_exchange`` into an ``exchange_k`` field and passes it to
:func:`tmesim.physics.diffusion.simulate_drug_transport`.

Many validation modules and experiment scripts call ``simulate_drug_transport``
DIRECTLY, bypassing the engine. On those paths ``PDEConfig(transvascular_exchange=True)``
is a **silent no-op**: the solver never receives ``exchange_k``, the run completes,
the numbers look plausible, and the result is BIT-IDENTICAL to the flag being off.
The solver is not at fault -- it was simply never told. It was measured, from
``experiments/dose_schedule.py``::

    exchange OFF                       auc_sum 1.554447793e+09
    flag set, exchange_k NOT passed     auc_sum 1.554447793e+09   <-- bit-identical
    honest ON (exchange_k passed)       auc_sum 6.86608671e+08

Anyone comparing "flag on" against "flag off" on such a path concludes "this
physics change has no effect". That conclusion is false and, critically, it is
**unfalsifiable from the output alone** -- the two runs are indistinguishable
because they are the same run.

What this module provides
=========================
:func:`check_flag_propagation` runs a callable twice (flag off, flag on) with
:func:`instrument` watching the solver boundary, and returns a
:class:`Verdict` in one of these categories:

``PROPAGATED``
    The solver ran and received ``exchange_k``, and the results differ. The
    comparison is real.
``INERT``
    The solver ran but was NEVER handed ``exchange_k`` even with the flag on.
    **This is the bug.** Any "flag on vs flag off" claim on this path is a
    comparison of a run to itself.
``REACHED_UNCHANGED``
    The solver ran AND received ``exchange_k``, and the results are still
    bit-identical. The flag genuinely reached the physics and genuinely did
    nothing (e.g. an avascular tissue where ``PS`` is zero everywhere). This is
    a claim about ROBUSTNESS.
``NO_PDE``
    The solver never ran at all. The flag is legitimately unreachable, because
    there is no transport solve on this path (analytic closed forms). This is a
    claim about UNREACHABILITY.
``NONDETERMINISTIC``
    Two flag-OFF runs of the same callable disagreed, so no off-vs-on comparison
    can be trusted. Undecidable, reported as such rather than guessed.

``REACHED_UNCHANGED`` and ``NO_PDE`` are deliberately different verdicts.
"Robust" and "unreachable" are completely different claims and collapsing them
is precisely the error this module exists to prevent. So is INERT vs NO_PDE:
both produce "nothing changed", and only one of them is honest.

Note that INERT is diagnosed STRUCTURALLY (the solver ran, and ``exchange_k``
was ``None`` on every call) rather than from bit-identical output. That matters:
identical output is a symptom that a coincidence could also produce, whereas
"the flag never crossed the function boundary" is the cause itself.

The declared allow-list
=======================
:data:`DECLARATIONS` is an EXPLICIT table: every known code path, its expected
category, and a written reason. ``NO_PDE`` can only ever be reached by writing
it down here with a justification -- it is never an inferred default, and there
is no path by which an unlisted module gets a free pass. ``tests/test_propagation.py``
holds the table to its word.

Static companion
================
:func:`scan_call_sites` answers the same question by reading source instead of
running it, which is cheap enough to sweep every module in the repo. It is a
weaker instrument -- see its docstring -- and where the two disagree the dynamic
:class:`Verdict` wins.

This module is purely observational. It imports the code it watches and
monkeypatches it only for the duration of a ``with`` block; it changes the
behaviour of nothing.
"""
from __future__ import annotations

import ast
import enum
import hashlib
import inspect
import sys
from contextlib import contextmanager
from dataclasses import dataclass, field, is_dataclass
from pathlib import Path
from typing import Any, Callable, Dict, Iterator, List, Optional, Sequence, Tuple

import numpy as np

__all__ = [
    "Category",
    "CallRecord",
    "Probe",
    "Verdict",
    "Declaration",
    "Exemption",
    "ALT_SOLVERS",
    "DECLARATIONS",
    "WIRING_POINTS",
    "WiringPoint",
    "wiring_points",
    "EXEMPTIONS",
    "sweep_repo",
    "undeclared_direct_callers",
    "declaration_for",
    "instrument",
    "check_flag_propagation",
    "digest",
    "CallSites",
    "scan_call_sites",
]


# --------------------------------------------------------------------------- #
# The solver boundary being watched
# --------------------------------------------------------------------------- #
#: ``(module, attribute)`` pairs instrumented by :func:`instrument`. The first
#: is THE solver -- the single function a physics flag must reach to mean
#: anything. The other two build the fields the flag is supposed to produce;
#: they are counted so that a diagnosis can distinguish "the exchange field was
#: never even built" from "it was built and then dropped on the floor".
WATCHED: Tuple[Tuple[str, str], ...] = (
    ("tmesim.physics.diffusion", "simulate_drug_transport"),
    ("tmesim.treatment.drugs", "delivery_source"),
    ("tmesim.treatment.drugs", "exchange_coefficient"),
    # Alternate solvers -- see ALT_SOLVERS.
    ("tmesim.molecular.immune_sink", "transport_with_immune_sink"),
    ("tmesim.validation.moncayo_explant", "bath_diffusion"),
)

SOLVER = "simulate_drug_transport"

#: Functions that integrate the transport PDE themselves instead of calling
#: :data:`SOLVER`. They import ``physics.diffusion``'s numerical primitives
#: (``_effective_diffusivity``, ``_build_diffusion_operator``,
#: ``_advective_divergence``) and run their own time loop, and NONE of them
#: takes an ``exchange_k`` argument at all.
#:
#: These matter because without them the checker would be actively misleading:
#: a path through one of these solves a PDE but never touches
#: ``simulate_drug_transport``, so it would be reported as NO_PDE -- "the flag
#: is legitimately unreachable" -- when the truth is "a PDE was solved and the
#: flag could not reach it", which is INERT. That is the exact conflation this
#: module exists to prevent, so they are watched and the NO_PDE verdict is
#: withheld whenever one of them ran.
ALT_SOLVERS: Tuple[str, ...] = (
    "transport_with_immune_sink",
    "bath_diffusion",
)


class Category(str, enum.Enum):
    """Verdicts. See the module docstring for what each one licenses you to say."""

    PROPAGATED = "PROPAGATED"
    INERT = "INERT"
    REACHED_UNCHANGED = "REACHED_UNCHANGED"
    NO_PDE = "NO_PDE"
    NONDETERMINISTIC = "NONDETERMINISTIC"


# --------------------------------------------------------------------------- #
# Digest -- a structural hash that is sensitive to every float in a result
# --------------------------------------------------------------------------- #
_MAX_DEPTH = 12


class _Digester:
    """Recursive, order-sensitive hash over arbitrary result objects.

    Every float that a run produced must influence the hash, or a real
    difference could be missed and an INERT path mistaken for a PROPAGATED one.
    Arrays are hashed from their raw bytes, so the comparison is BIT-level and
    not ``allclose`` -- the failure mode being hunted produces bit-identical
    output, so anything looser would be unable to see the distinction it is
    being asked about.

    Objects the walker does not understand are recorded by type name only and
    set :attr:`opaque`, which is surfaced on the :class:`Verdict` rather than
    hidden: an opaque payload means the digest saw less than the whole result.
    """

    def __init__(self) -> None:
        self._h = hashlib.blake2b(digest_size=16)
        self.opaque = False

    def _feed(self, tag: str, payload: bytes = b"") -> None:
        self._h.update(tag.encode("utf-8"))
        self._h.update(b"\x1f")
        self._h.update(payload)
        self._h.update(b"\x1e")

    def add(self, obj: Any, depth: int = 0) -> "_Digester":
        if depth > _MAX_DEPTH:
            self._feed("depth-cut")
            self.opaque = True
            return self

        if obj is None:
            self._feed("none")
        elif isinstance(obj, (bool, np.bool_)):
            self._feed("bool", b"\x01" if obj else b"\x00")
        elif isinstance(obj, (str,)):
            self._feed("str", obj.encode("utf-8"))
        elif isinstance(obj, (bytes, bytearray)):
            self._feed("bytes", bytes(obj))
        elif isinstance(obj, (int, np.integer)):
            self._feed("int", repr(int(obj)).encode())
        elif isinstance(obj, (float, np.floating)):
            # float64 bytes: distinguishes -0.0 from 0.0 and preserves NaN payloads.
            self._feed("float", np.float64(obj).tobytes())
        elif isinstance(obj, (complex, np.complexfloating)):
            self._feed("complex", np.complex128(obj).tobytes())
        elif isinstance(obj, np.ndarray):
            arr = np.ascontiguousarray(obj)
            self._feed("ndarray", f"{arr.dtype.str}|{arr.shape}".encode())
            if arr.dtype == object:
                self.opaque = True
                for item in arr.ravel().tolist():
                    self.add(item, depth + 1)
            else:
                self._feed("bytes", arr.tobytes())
        elif isinstance(obj, dict):
            self._feed("dict", repr(len(obj)).encode())
            for key in obj:  # insertion order is part of the result
                self.add(key, depth + 1)
                self.add(obj[key], depth + 1)
        elif isinstance(obj, (list, tuple)):
            self._feed(type(obj).__name__, repr(len(obj)).encode())
            for item in obj:
                self.add(item, depth + 1)
        elif isinstance(obj, (set, frozenset)):
            self._feed("set", repr(len(obj)).encode())
            for item in sorted(obj, key=repr):
                self.add(item, depth + 1)
        elif is_dataclass(obj) and not isinstance(obj, type):
            self._feed("dataclass", type(obj).__name__.encode())
            for f in obj.__dataclass_fields__:  # type: ignore[attr-defined]
                self.add(f, depth + 1)
                self.add(getattr(obj, f, None), depth + 1)
        elif hasattr(obj, "__dict__") and vars(obj):
            self._feed("object", type(obj).__name__.encode())
            for key, val in sorted(vars(obj).items()):
                self.add(key, depth + 1)
                self.add(val, depth + 1)
        else:
            # Unknown and unwalkable. Record only the type -- never ``repr``,
            # which for most objects embeds a memory address and would make two
            # identical runs look different.
            self._feed("opaque", type(obj).__name__.encode())
            self.opaque = True
        return self

    def hexdigest(self) -> str:
        return self._h.hexdigest()


def digest(obj: Any) -> str:
    """Bit-level structural hash of a result object. See :class:`_Digester`."""
    return _Digester().add(obj).hexdigest()


# --------------------------------------------------------------------------- #
# Instrumentation
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class CallRecord:
    """One observed call across the instrumented boundary."""

    target: str
    #: Only meaningful for the solver: was ``exchange_k`` supplied and not None?
    exchange_k: bool = False
    #: Only meaningful for the solver: was a time profile supplied (the
    #: ``time_resolved`` flag's equivalent hand-off)?
    source_profile: bool = False
    #: Digest of what the call returned.
    result: str = ""


@dataclass
class Probe:
    """Everything observed at the solver boundary during one run."""

    calls: List[CallRecord] = field(default_factory=list)
    opaque: bool = False

    def of(self, target: str) -> List[CallRecord]:
        return [c for c in self.calls if c.target == target]

    @property
    def solver_calls(self) -> int:
        """How many times the transport PDE was actually solved."""
        return len(self.of(SOLVER))

    @property
    def exchange_k_calls(self) -> int:
        """Solver calls that were handed a transvascular-exchange field."""
        return sum(1 for c in self.of(SOLVER) if c.exchange_k)

    @property
    def source_profile_calls(self) -> int:
        """Solver calls that were handed a time-resolved plasma profile."""
        return sum(1 for c in self.of(SOLVER) if c.source_profile)

    @property
    def exchange_coefficient_calls(self) -> int:
        """Times the ``PS`` field was BUILT (which does not mean it was used)."""
        return len(self.of("exchange_coefficient"))

    @property
    def alt_solver_calls(self) -> int:
        """Transport solves that went through a re-implemented solver instead.

        Non-zero means a PDE WAS integrated even if :data:`SOLVER` was never
        called -- so "no PDE ran" would be a false statement.
        """
        return sum(len(self.of(name)) for name in ALT_SOLVERS)

    def solver_digest(self) -> str:
        """Hash of every solver return value, in call order."""
        return digest([c.result for c in self.of(SOLVER)])


def _resolve(mod_name: str, attr: str) -> Any:
    import importlib

    return getattr(importlib.import_module(mod_name), attr)


def _make_wrapper(original: Callable[..., Any], name: str, probe: Probe) -> Callable[..., Any]:
    import functools

    try:
        sig: Optional[inspect.Signature] = inspect.signature(original)
    except (TypeError, ValueError):  # pragma: no cover - defensive
        sig = None

    @functools.wraps(original)
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        exchange_k = False
        source_profile = False
        if name == SOLVER and sig is not None:
            try:
                bound = sig.bind(*args, **kwargs)
                exchange_k = bound.arguments.get("exchange_k", None) is not None
                source_profile = bound.arguments.get("source_profile", None) is not None
            except TypeError:
                # Let the real call raise the real error; record what we can.
                exchange_k = kwargs.get("exchange_k", None) is not None
                source_profile = kwargs.get("source_profile", None) is not None
        out = original(*args, **kwargs)
        d = _Digester().add(out)
        probe.opaque = probe.opaque or d.opaque
        probe.calls.append(
            CallRecord(
                target=name,
                exchange_k=exchange_k,
                source_profile=source_profile,
                result=d.hexdigest(),
            )
        )
        return out

    wrapper.__wrapped__ = original  # type: ignore[attr-defined]
    return wrapper


@contextmanager
def instrument(watched: Sequence[Tuple[str, str]] = WATCHED) -> Iterator[Probe]:
    """Count and fingerprint every call across the solver boundary.

    Patches the defining module AND every already-imported module that did
    ``from ... import simulate_drug_transport`` -- a from-import binds a new
    name, so patching only the source module would miss exactly the direct
    callers this exists to catch. Modules imported *during* the block pick the
    wrapper up from the (patched) source module; on exit ``sys.modules`` is
    rescanned so those late bindings are restored too. Lazy, function-scoped
    imports are fine for the same reason -- they run inside the block.

    LIMITATION. A reference bound into a LOCAL variable or a closure BEFORE the
    block opens cannot be reached (there is no module attribute to rebind), so a
    call through it is invisible and the path would look quieter than it is.
    Every solver call site in this repo resolves through a module attribute, so
    this is currently theoretical -- but it is the shape of the hole.

    Restores everything it touched, including on exception.
    """
    probe = Probe()
    installed: List[Tuple[Any, str, Any, Any]] = []  # (module, attr, original, wrapper)

    try:
        for mod_name, attr in watched:
            try:
                original = _resolve(mod_name, attr)
            except Exception:  # pragma: no cover - a watched target vanished
                continue
            wrapper = _make_wrapper(original, attr, probe)
            for module in list(sys.modules.values()):
                if module is None:
                    continue
                try:
                    if getattr(module, attr, None) is original:
                        setattr(module, attr, wrapper)
                        installed.append((module, attr, original, wrapper))
                except Exception:  # pragma: no cover - exotic module objects
                    continue
        yield probe
    finally:
        for module, attr, original, wrapper in installed:
            try:
                if getattr(module, attr, None) is wrapper:
                    setattr(module, attr, original)
            except Exception:  # pragma: no cover
                pass
        # Anything imported mid-block that bound the wrapper.
        seen = {(id(m), a) for m, a, _, _ in installed}
        for _, attr, original, wrapper in list(installed):
            for module in list(sys.modules.values()):
                if module is None or (id(module), attr) in seen:
                    continue
                try:
                    if getattr(module, attr, None) is wrapper:
                        setattr(module, attr, original)
                except Exception:  # pragma: no cover
                    pass


# --------------------------------------------------------------------------- #
# The verdict
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Verdict:
    """The mechanical answer to "am I comparing a run to itself?"."""

    category: Category
    label: str = ""
    solver_calls_off: int = 0
    solver_calls_on: int = 0
    exchange_k_calls_on: int = 0
    exchange_coefficient_calls_on: int = 0
    digest_off: str = ""
    digest_on: str = ""
    #: True when the digest could not fully see some part of a result.
    opaque: bool = False
    detail: str = ""

    @property
    def results_differ(self) -> bool:
        return self.digest_off != self.digest_on

    @property
    def trustworthy(self) -> bool:
        """Is an "on vs off" comparison on this path meaningful at all?

        False for INERT (the comparison is a run against itself) and for
        NONDETERMINISTIC (the comparison is noise). True for NO_PDE and
        REACHED_UNCHANGED -- in both of those "nothing changed" is an honest
        statement, though of two different kinds.
        """
        return self.category not in (Category.INERT, Category.NONDETERMINISTIC)

    def __str__(self) -> str:  # pragma: no cover - human output
        head = f"{self.category.value:<18} {self.label}".rstrip()
        return f"{head}\n    {self.detail}"


def _run_digest(probe: Probe, returned: Any) -> str:
    """Fingerprint a run from BOTH ends.

    The solver's own outputs are the primary evidence -- they sit exactly at the
    boundary the flag has to cross. The callable's return value is folded in as
    well, because a path can post-process the solver output and a difference
    might only be visible downstream (the ``auc_sum`` in the measured example
    above is a returned summary, not a raw solver field).
    """
    return digest((probe.solver_digest(), digest(returned)))


def check_flag_propagation(
    run: Callable[[bool], Any],
    *,
    label: str = "",
    determinism_check: bool = True,
) -> Verdict:
    """Run ``run`` with the flag off and on; report whether the flag got through.

    Parameters
    ----------
    run : callable taking one bool
        ``run(False)`` must do the work with the physics flag OFF and
        ``run(True)`` the same work with it ON -- typically by building a
        ``PDEConfig(transvascular_exchange=flag)`` and handing it onward exactly
        as the real code path does. Whatever it returns is folded into the
        comparison, so returning a summary the module actually reports (an AUC
        sum, a scorecard dict) makes the check strictly sharper. Returning
        ``None`` is fine: the solver's own outputs still carry the comparison.
    label : str
        Name used in the verdict, for test output.
    determinism_check : bool, default True
        Run the OFF case twice first and require the two to agree. Without this
        a stochastic path would show a difference between off and on that has
        nothing to do with the flag, and be reported as PROPAGATED -- a false
        pass, which is the same class of error this module exists to stop.

    Returns
    -------
    Verdict
    """
    with instrument() as p_off:
        r_off = run(False)
    d_off = _run_digest(p_off, r_off)

    if determinism_check:
        with instrument() as p_ctl:
            r_ctl = run(False)
        d_ctl = _run_digest(p_ctl, r_ctl)
        if d_ctl != d_off:
            return Verdict(
                category=Category.NONDETERMINISTIC,
                label=label,
                solver_calls_off=p_off.solver_calls,
                digest_off=d_off,
                digest_on=d_ctl,
                opaque=p_off.opaque or p_ctl.opaque,
                detail=(
                    "two flag-OFF runs of this callable disagreed, so an "
                    "off-vs-on comparison cannot separate the flag from the "
                    "noise. Seed the path (or hold the RNG fixed) and re-check; "
                    "no propagation claim is decidable until then."
                ),
            )

    with instrument() as p_on:
        r_on = run(True)
    d_on = _run_digest(p_on, r_on)

    common = dict(
        label=label,
        solver_calls_off=p_off.solver_calls,
        solver_calls_on=p_on.solver_calls,
        exchange_k_calls_on=p_on.exchange_k_calls,
        exchange_coefficient_calls_on=p_on.exchange_coefficient_calls,
        digest_off=d_off,
        digest_on=d_on,
        opaque=p_off.opaque or p_on.opaque,
    )

    # (C) No transport solve on this path at all. Legitimately unreachable --
    #     but only ever an honest verdict when it has been DECLARED as such,
    #     and only when no RE-IMPLEMENTED solver ran either.
    if p_off.solver_calls == 0 and p_on.solver_calls == 0:
        if p_on.alt_solver_calls:
            return Verdict(
                category=Category.INERT,
                detail=(
                    f"simulate_drug_transport never ran, but {p_on.alt_solver_calls} "
                    "transport solve(s) went through a RE-IMPLEMENTED solver "
                    f"({', '.join(ALT_SOLVERS)}) that takes no exchange_k "
                    "argument at all. A PDE was integrated and the flag could "
                    "not reach it. This is INERT, not NO_PDE: reporting it as "
                    "'no PDE ran' would be false."
                ),
                **common,
            )
        return Verdict(
            category=Category.NO_PDE,
            detail=(
                "simulate_drug_transport was never called, with the flag off or "
                "on. There is no transport solve on this path for a PDE flag to "
                "reach (analytic closed form, or the work happens elsewhere). "
                "'Nothing changed' here means UNREACHABLE, not robust -- it is "
                "no evidence at all about the physics."
            ),
            **common,
        )

    # (B) The solver ran and the flag never crossed the boundary. THE BUG.
    if p_on.exchange_k_calls == 0:
        if d_off != d_on:
            return Verdict(
                category=Category.PROPAGATED,
                detail=(
                    f"results differ across {p_on.solver_calls} solver call(s) "
                    "even though exchange_k was never passed, so the flag is "
                    "reaching the physics by some OTHER route. Real, but not "
                    "the route this checker models -- find it before trusting it."
                ),
                **common,
            )
        built = p_on.exchange_coefficient_calls
        return Verdict(
            category=Category.INERT,
            detail=(
                f"simulate_drug_transport ran {p_on.solver_calls} time(s) with "
                "the flag ON and was handed exchange_k=None every time"
                + (
                    f" (exchange_coefficient WAS built {built} time(s) and then "
                    "dropped)"
                    if built
                    else " (exchange_coefficient was never even built)"
                )
                + ". The flag is a silent no-op on this path: the two runs are "
                "bit-identical because they are the same run. Any 'on vs off' "
                "conclusion drawn here is a comparison of a run to itself. Fix "
                "by routing through sim.engine.simulate_response, or by passing "
                "exchange_coefficient(micro, drug) as exchange_k= explicitly."
            ),
            **common,
        )

    # (A) The flag crossed the boundary and moved the numbers.
    if d_off != d_on:
        return Verdict(
            category=Category.PROPAGATED,
            detail=(
                f"exchange_k reached the solver on {p_on.exchange_k_calls} of "
                f"{p_on.solver_calls} call(s) and the results changed. An "
                "on-vs-off comparison on this path is real."
            ),
            **common,
        )

    # The flag reached the physics and genuinely did nothing. ROBUST, not unreachable.
    return Verdict(
        category=Category.REACHED_UNCHANGED,
        detail=(
            f"exchange_k reached the solver on {p_on.exchange_k_calls} of "
            f"{p_on.solver_calls} call(s), and the results are still "
            "bit-identical. The flag was genuinely applied and genuinely had no "
            "effect (e.g. PS is zero everywhere on this geometry). This is a "
            "ROBUSTNESS result -- distinct from NO_PDE, where the flag was never "
            "applied at all."
        ),
        **common,
    )


# --------------------------------------------------------------------------- #
# Static companion: where are the call sites?
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class CallSites:
    """What a source file says about how it reaches the solver.

    Read straight off the AST of one file (optionally one function inside it).
    Cheap enough to sweep the whole repo, and it needs no fixtures, but it is
    strictly weaker than a :class:`Verdict`: it sees call SITES, not call PATHS.
    It cannot tell whether a site is reached at runtime, it does not follow calls
    into helper modules, and a dynamic ``exchange_k=`` built behind a conditional
    still counts as passed. Use it to triage; use
    :func:`check_flag_propagation` to conclude.
    """

    target: str
    #: Line numbers of direct ``simulate_drug_transport(...)`` calls.
    direct_calls: Tuple[int, ...] = ()
    #: Subset of ``direct_calls`` that pass an ``exchange_k=`` keyword.
    direct_calls_with_exchange_k: Tuple[int, ...] = ()
    #: Line numbers of ``simulate_response(...)`` calls (the engine path).
    engine_calls: Tuple[int, ...] = ()

    @property
    def bare_direct_calls(self) -> Tuple[int, ...]:
        """Direct solver calls that do NOT pass ``exchange_k``. Each one is a
        line at which a ``PDEConfig`` physics flag stops existing."""
        withk = set(self.direct_calls_with_exchange_k)
        return tuple(ln for ln in self.direct_calls if ln not in withk)

    @property
    def predicted(self) -> Category:
        """Best static guess at the dynamic category.

        A bare direct call outranks an engine call: if a path can reach the
        solver without the flag, the honest prediction is INERT even when some
        other line in the same file does it correctly. Under-claiming here is
        the safe direction -- the cost of a wrong INERT prediction is one
        investigation, the cost of a wrong PROPAGATED prediction is a published
        false conclusion.
        """
        if self.bare_direct_calls:
            return Category.INERT
        if self.direct_calls_with_exchange_k or self.engine_calls:
            return Category.PROPAGATED
        return Category.NO_PDE


def _called_name(node: ast.Call) -> str:
    func = node.func
    if isinstance(func, ast.Name):
        return func.id
    if isinstance(func, ast.Attribute):
        return func.attr
    return ""


def scan_call_sites(path: str | Path, symbol: Optional[str] = None) -> CallSites:
    """AST-scan one file (or one function in it) for solver call sites.

    ``symbol`` restricts the scan to a single top-level ``def``/``async def``,
    which is what makes a file like ``validation/vessel_profile.py`` -- one
    function on the engine path and another calling the solver directly --
    describable without averaging the two into a meaningless whole-file answer.
    """
    path = Path(path)
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))

    root: ast.AST = tree
    if symbol is not None:
        found = None
        for node in ast.walk(tree):
            if (
                isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
                and node.name == symbol
            ):
                found = node
                break
        if found is None:
            raise LookupError(f"{path}: no function named {symbol!r}")
        root = found

    direct: List[int] = []
    with_k: List[int] = []
    engine: List[int] = []
    for node in ast.walk(root):
        if not isinstance(node, ast.Call):
            continue
        name = _called_name(node)
        if name == SOLVER:
            direct.append(node.lineno)
            if any(kw.arg == "exchange_k" for kw in node.keywords):
                with_k.append(node.lineno)
        elif name == "simulate_response":
            engine.append(node.lineno)

    target = str(path) if symbol is None else f"{path}::{symbol}"
    return CallSites(
        target=target,
        direct_calls=tuple(sorted(direct)),
        direct_calls_with_exchange_k=tuple(sorted(with_k)),
        engine_calls=tuple(sorted(engine)),
    )


# --------------------------------------------------------------------------- #
# The declared allow-list
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Declaration:
    """One code path, its expected category, and WHY.

    A ``reason`` is mandatory and enforced at construction. For ``NO_PDE`` in
    particular this is the whole point: "the flag cannot reach here" must be an
    asserted, argued claim that someone signed, never a silence that a checker
    filled in by default. A path nobody has thought about produces no entry, and
    a path with no entry is a test failure -- not a pass.
    """

    #: Repo-relative path to the file.
    path: str
    #: Top-level function to scan, or None for the whole file.
    symbol: Optional[str] = None
    #: Expected category.
    category: Category = Category.INERT
    #: Why. Required.
    reason: str = ""

    def __post_init__(self) -> None:
        if len(self.reason.strip()) < 25:
            raise ValueError(
                f"{self.name}: a Declaration needs a written reason "
                "(>=25 chars). 'Unreachable' and 'robust' are different claims "
                "and neither may be asserted by default."
            )
        if not isinstance(self.category, Category):
            raise TypeError(f"{self.name}: category must be a Category")

    @property
    def name(self) -> str:
        return self.path if self.symbol is None else f"{self.path}::{self.symbol}"


#: Every known code path on which ``PDEConfig.transvascular_exchange`` might be
#: set, with the category it is asserted to be in. Adding a module to the repo
#: that touches the solver and NOT adding it here is a test failure.
DECLARATIONS: Tuple[Declaration, ...] = (
    # --- (A) engine path: the flag propagates, results MUST change ------------
    Declaration(
        path="tmesim/sim/engine.py",
        symbol="simulate_response",
        category=Category.PROPAGATED,
        reason=(
            "The one and only place that reads cfg.pde.transvascular_exchange, "
            "turns it into exchange_coefficient(work, drug) and passes it to the "
            "solver as exchange_k. Every honest 'flag on' number in this repo "
            "comes through here."
        ),
    ),
    Declaration(
        path="tmesim/validation/vessel_profile.py",
        symbol="analyze_drug",
        category=Category.PROPAGATED,
        reason=(
            "Routes through sim.engine.simulate_response, so the flag reaches "
            "the solver. This is the half of vessel_profile.py that is safe; the "
            "other half is not, which is why this file is declared per-function."
        ),
    ),
    # --- (B) direct-solver path: the flag is INERT unless exchange_k is passed -
    Declaration(
        path="tmesim/validation/pk_reality.py",
        category=Category.INERT,
        reason=(
            "Calls simulate_drug_transport directly (pk_reality.py:218) with no "
            "exchange_k. Setting the PDEConfig flag here changes nothing; a "
            "re-validation must pass exchange_coefficient(micro, dp) explicitly."
        ),
    ),
    Declaration(
        path="tmesim/validation/kill_calibration.py",
        category=Category.INERT,
        reason=(
            "Calls simulate_drug_transport directly (kill_calibration.py:681) "
            "with no exchange_k. The calibration constant it fits is therefore "
            "an exchange-OFF constant whatever the config says."
        ),
    ),
    Declaration(
        path="tmesim/validation/immune_sink_impact.py",
        category=Category.INERT,
        reason=(
            "Calls simulate_drug_transport directly at three sites (252, 272, "
            "286) with no exchange_k. Its own per-voxel immune sink is passed by "
            "hand, which is exactly the pattern exchange_k would need and does "
            "not get."
        ),
    ),
    Declaration(
        path="tmesim/validation/vessel_profile.py",
        symbol="isolated_vessel_profile",
        category=Category.INERT,
        reason=(
            "Builds its own PDEConfig and calls simulate_drug_transport directly "
            "(vessel_profile.py:646) with no exchange_k, unlike analyze_drug in "
            "the same file. A whole-file verdict would average these two into "
            "something false."
        ),
    ),
    Declaration(
        path="tmesim/validation/numerical_rigor.py",
        category=Category.INERT,
        reason=(
            "Calls simulate_drug_transport directly (314, 456) with no "
            "exchange_k, so its convergence and CFL results characterise the "
            "exchange-OFF operator only. It also calls simulate_response, which "
            "is precisely why a file-level 'it uses the engine' glance is unsafe."
        ),
    ),
    Declaration(
        path="experiments/regimen.py",
        category=Category.INERT,
        reason=(
            "Calls simulate_drug_transport directly at three sites (415, 826, "
            "853) with no exchange_k, alongside engine calls elsewhere in the "
            "file. The direct sites win: a bare call is a line at which the flag "
            "ceases to exist."
        ),
    ),
    Declaration(
        path="experiments/dose_schedule.py",
        category=Category.INERT,
        reason=(
            "The module the failure was MEASURED on. unit_transport (line 307) "
            "calls the solver directly with no exchange_k: auc_sum was "
            "1.554447793e+09 with the flag off and bit-identical with it on, "
            "against 6.86608671e+08 when exchange_k was actually passed."
        ),
    ),
    # --- (C) no PDE at all: DECLARED unreachable, with a reason ---------------
    Declaration(
        path="tmesim/validation/gist_imatinib.py",
        category=Category.NO_PDE,
        reason=(
            "Analytic closed forms only -- slab delivery efficiency and a "
            "tissue-partition Kpu. No transport solve exists on this path, so a "
            "PDE flag has nothing to reach. 'Nothing changed' here is a statement "
            "about this module's scope and is NOT evidence that exchange is "
            "negligible."
        ),
    ),
    Declaration(
        path="tmesim/validation/glioma_regional_pk.py",
        category=Category.NO_PDE,
        reason=(
            "Analytic/statistical only: descriptor-derived BBB penetrance and "
            "permutation tests over measured regional ratios. Never solves a PDE, "
            "so this module can neither confirm nor refute the exchange term."
        ),
    ),
)


@dataclass(frozen=True)
class Exemption:
    """A file with bare solver calls that is deliberately outside the sweep."""

    path: str
    reason: str

    def __post_init__(self) -> None:
        if len(self.reason.strip()) < 25:
            raise ValueError(f"{self.path}: an Exemption needs a written reason")


#: Files whose bare ``simulate_drug_transport`` calls are NOT propagation
#: failures. Kept deliberately tiny -- every entry is a hole in the sweep.
EXEMPTIONS: Tuple[Exemption, ...] = (
    Exemption(
        path="tmesim/physics/diffusion.py",
        reason=(
            "The module that DEFINES the solver, plus its own __main__ "
            "self-test. The solver is the destination a flag has to reach, so "
            "'does the flag reach it from here' is not a meaningful question; "
            "these calls construct exchange_k by hand where they mean to."
        ),
    ),
)


@dataclass(frozen=True)
class WiringPoint:
    """A file that deliberately PASSES ``exchange_k`` into the solver."""

    path: str
    #: True only for the shipped runtime path. See :data:`WIRING_POINTS`.
    production: bool
    reason: str

    def __post_init__(self) -> None:
        if len(self.reason.strip()) < 25:
            raise ValueError(f"{self.path}: a WiringPoint needs a written reason")


#: Every file that hands ``exchange_k`` to the solver, and why.
#:
#: The invariant this table protects is unchanged and still strict: there is
#: EXACTLY ONE production wiring point. A second one would mean two places
#: decide whether the gradient-driven delivery law is active, which is how a
#: flag ends up silently on in one path and off in another -- the whole failure
#: mode this module exists to catch.
#:
#: What the table adds is a second, non-production category. A validation
#: experiment that measures an INERT path has no choice but to inject
#: ``exchange_k`` by hand: setting ``PDEConfig(transvascular_exchange=True)``
#: on such a path is a silent no-op, so an experiment that did the obvious
#: thing would produce a clean, confident, meaningless "nothing changed"
#: result. Those scripts are doing exactly the right thing and must not be
#: driven to do the wrong thing by a guard.
#:
#: They are enumerated individually rather than exempted by directory. A rule
#: like "anything under experiments/ may wire it up" would let an undeclared
#: no-op creep back in precisely where it is least likely to be noticed, and
#: the point of this module is that every such path is written down.
WIRING_POINTS: Tuple[WiringPoint, ...] = (
    WiringPoint(
        path="tmesim/sim/engine.py",
        production=True,
        reason=(
            "The shipped runtime path, and the only place that translates "
            "PDEConfig.transvascular_exchange into an exchange_k argument. "
            "Everything a user runs reaches the solver through here."
        ),
    ),
    WiringPoint(
        path="experiments/vessel_profile_arms.py",
        production=False,
        reason=(
            "Re-measures results/vessel_profile.md under both delivery laws. "
            "Its subject is declared INERT, so the PDEConfig flag alone would "
            "be a silent no-op and the experiment would report a false "
            "'no change'. It injects exchange_k deliberately, and asserts the "
            "two arms' fields are not bit-identical before drawing any "
            "conclusion."
        ),
    ),
    WiringPoint(
        path="experiments/dose_schedule_arms.py",
        production=False,
        reason=(
            "Re-measures results/dose_schedule.md under both delivery laws. "
            "Its subject, experiments/dose_schedule.py, is the module the "
            "silent no-op was first MEASURED on and is declared INERT above, "
            "so the PDEConfig flag alone would change nothing there. This "
            "script injects exchange_coefficient(work, drug) as exchange_k at "
            "a single call site, asserts per solve that exchange_k crossed the "
            "boundary on every call in the ON arm and none in the OFF arm, "
            "requires the two arms' AUC fields to differ bit-wise, and "
            "additionally reconciles the injection against the PRODUCTION "
            "wiring by re-running sim.engine.simulate_response with "
            "transvascular_exchange=True."
        ),
    ),
    WiringPoint(
        path="experiments/boundary_layer_resolution.py",
        production=False,
        reason=(
            "Measures the discretisation error of the perivascular boundary "
            "layer delta = sqrt(D_eff/(PS + clearance + binding)) against dx. "
            "PS is the exchange term, so delta only EXISTS under the corrected "
            "law -- the experiment cannot pose its own question without "
            "exchange_k, and on its controlled slab there is no PDEConfig "
            "path to route through. It injects at a single call site and its "
            "stage-0 trap check asserts the ON and OFF fields differ before "
            "any error law is fitted."
        ),
    ),
)


def wiring_points(production_only: bool = False) -> frozenset:
    """Paths permitted to pass ``exchange_k``; the production one if asked."""
    return frozenset(
        w.path for w in WIRING_POINTS if w.production or not production_only
    )


def declaration_for(path: str, symbol: Optional[str] = None) -> Declaration:
    """Look up a declaration, or raise. There is no default."""
    name = path if symbol is None else f"{path}::{symbol}"
    for d in DECLARATIONS:
        if d.name == name:
            return d
    raise LookupError(
        f"{name} has no propagation declaration. Add one to "
        "tmesim.validation.propagation.DECLARATIONS with a written reason -- "
        "an undeclared path is a failure, not a pass."
    )


#: Directories swept by :func:`sweep_repo`.
SWEPT_ROOTS: Tuple[str, ...] = ("tmesim", "experiments", "scripts", "app")


def sweep_repo(root: str | Path = ".") -> Dict[str, CallSites]:
    """Static-scan the repo for every file containing a solver call site.

    Returns ``{repo-relative path: CallSites}``. Cheap (AST only) and the reason
    the guard does not depend on someone remembering to register a new module:
    a file that starts calling the solver shows up here whether or not anyone
    thought about it.
    """
    root = Path(root)
    found: Dict[str, CallSites] = {}
    for base in SWEPT_ROOTS:
        base_dir = root / base
        if not base_dir.is_dir():
            continue
        for py in sorted(base_dir.rglob("*.py")):
            try:
                cs = scan_call_sites(py)
            except (SyntaxError, UnicodeDecodeError):  # pragma: no cover
                continue
            if cs.direct_calls or cs.engine_calls:
                found[str(py.relative_to(root))] = cs
    return found


#: The private numerics a re-implemented solver has to borrow. Importing one of
#: these is the signature of a hand-rolled time loop -- the thing that produces
#: a PDE solve the flag cannot reach.
DIFFUSION_PRIMITIVES: Tuple[str, ...] = (
    "_build_diffusion_operator",
    "_advective_divergence",
)


def reimplementation_candidates(root: str | Path = ".") -> Dict[str, Tuple[str, ...]]:
    """Files that import ``physics.diffusion``'s private numerics.

    A file here is building its own transport loop rather than calling
    :data:`SOLVER`, which means a physics flag cannot reach it no matter how
    the config is set. Returns ``{path: (primitives imported,)}``.

    This is how :data:`ALT_SOLVERS` is kept from going stale: a new hand-rolled
    solver shows up in this sweep even though it calls nothing the checker
    watches.
    """
    root = Path(root)
    out: Dict[str, Tuple[str, ...]] = {}
    for base in SWEPT_ROOTS:
        base_dir = root / base
        if not base_dir.is_dir():
            continue
        for py in sorted(base_dir.rglob("*.py")):
            try:
                tree = ast.parse(py.read_text(encoding="utf-8"), filename=str(py))
            except (SyntaxError, UnicodeDecodeError):  # pragma: no cover
                continue
            hits = {
                alias.name
                for node in ast.walk(tree)
                if isinstance(node, ast.ImportFrom)
                for alias in node.names
                if alias.name in DIFFUSION_PRIMITIVES
            }
            if hits:
                out[str(py.relative_to(root))] = tuple(sorted(hits))
    return out


def undeclared_direct_callers(root: str | Path = ".") -> List[str]:
    """Files that call the solver WITHOUT ``exchange_k`` and are neither
    declared in :data:`DECLARATIONS` nor listed in :data:`EXEMPTIONS`.

    Any name this returns is a code path on which a ``PDEConfig`` physics flag
    is currently a silent no-op that nobody has written down. The correct
    response is to add a declaration (or fix the call), never to widen the
    exemption list.
    """
    declared = {d.path for d in DECLARATIONS}
    exempt = {e.path for e in EXEMPTIONS}
    return [
        path
        for path, cs in sweep_repo(root).items()
        if cs.bare_direct_calls and path not in declared and path not in exempt
    ]
