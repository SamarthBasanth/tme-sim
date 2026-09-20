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

More than one flag
==================
Everything above was written for ONE flag reaching ONE solver. It is now
parameterised by :class:`Flag`, because a second switch of the identical shape
landed with ``fd2a9e7``: ``compile_drug.build_profile(two_route: bool = True)``
decides whether the two-route cell term is applied to a compiled drug's
diffusivity, and until this parameterisation it had no declaration table and no
"exactly one production wiring point" invariant. It was not red; it was blind.

The two flags are DUALS, and the generalisation turns on noticing that:

``transvascular_exchange``
    The destination (``simulate_drug_transport``) is called from all over the
    repo, and the question is whether the flag CAME ALONG -- it travels as the
    ``exchange_k`` keyword. A call without the carrier is the defect.
``two_route``
    The destination (``two_route.two_route_diffusivity``) has no carrier: the
    flag's whole content is whether the destination is called AT ALL. So a call
    is an arrival and there is no such thing as a bare one. The defect is the
    other polarity -- a file that computes an effective diffusivity out of
    ``brownian``/``ecm_mesh`` primitives and never reaches the cell term. Those
    are :attr:`Flag.precursors`, and a file holding them with no arrival is
    STRANDED -- it produced a pre-promotion number that the flag could not have
    changed. Stranded lands in ``INERT``, which is the same verdict for the
    same reason: a PDE was solved / a diffusivity was computed, and the flag
    could not reach it. It is emphatically not ``NO_PDE``.

:attr:`Flag.carrier` is what expresses the difference. ``carrier=None`` means
"arrival is the call", and every polarity above falls out of that one field
rather than out of a second copy of the machinery.

One thing did NOT generalise, and it is data rather than logic:
:attr:`Flag.wiring_scope`. See that field for the argument.

What this module does NOT guard
===============================
Code, only. It decides reachability and identity questions about SOURCE. The
prose that quotes this engine's results is guarded separately, by
``tests/test_claim_guard.py``, and that guard is narrower still -- see its
"What this cannot do" section. Neither reaches a published NUMBER and neither
reaches an INTERPRETATION. A guard is not coverage.

This module is purely observational. It imports the code it watches and
monkeypatches it only for the duration of a ``with`` block; it changes the
behaviour of nothing.
"""
from __future__ import annotations

import ast
import enum
import hashlib
import inspect
import re
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
    "undeclared_paths",
    "declaration_for",
    "instrument",
    "check_flag_propagation",
    "digest",
    "CallSites",
    "scan_call_sites",
    "Flag",
    "FLAGS",
    "EXCHANGE",
    "TWO_ROUTE",
    "flag_named",
    "WatchedTargetMissing",
    "unresolvable_watched",
    "unwatched_alternates",
    "cited_lines",
    "recognised_lines",
    "stale_citations",
    "undeclared_sites",
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
    # Lives in app/, which is an OPTIONAL extra, so this is a soft reference:
    # WATCHED holds strings and ``instrument`` skips a target it cannot resolve,
    # so a checkout without '.[app]' simply does not watch it. That is safe here
    # and only here -- such a checkout cannot run the reach map either, so there
    # is no path that solves transport unwatched. It is NOT a licence to point
    # WATCHED at other optional modules.
    ("app.server.reach", "steady_state_concentration"),
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
    # Not a time loop at all: assembles (K - L + Adv) c = S from the same
    # primitives and solves it once. It is listed here for exactly the reason
    # above -- it SOLVES TRANSPORT, so letting it read as NO_PDE ("no transport
    # solve happens here") would be false in the one direction this module
    # exists to prevent. It takes exchange_k explicitly and defaults it to None.
    "steady_state_concentration",
)


class Category(str, enum.Enum):
    """Verdicts. See the module docstring for what each one licenses you to say."""

    PROPAGATED = "PROPAGATED"
    INERT = "INERT"
    REACHED_UNCHANGED = "REACHED_UNCHANGED"
    NO_PDE = "NO_PDE"
    NONDETERMINISTIC = "NONDETERMINISTIC"

    #: ``NO_PDE`` under its general name. The claim is "the flag is legitimately
    #: unreachable on this path"; for a PDE flag that is spelled "no transport
    #: solve happens here", and for :data:`TWO_ROUTE` it is spelled "no
    #: diffusivity is computed here". Declared as an Enum ALIAS deliberately:
    #: ``Category.UNREACHABLE is Category.NO_PDE`` is True, so the two spellings
    #: cannot drift into two subtly different verdicts, which is the exact
    #: failure this module exists to prevent one level down.
    UNREACHABLE = "NO_PDE"


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
    #: Only meaningful for the destination: did the flag ARRIVE on this call?
    #: For a carrier-bearing flag that means the carrier was supplied and not
    #: None; for a carrier-less flag the call itself is the arrival, so this is
    #: True on every destination call. Named for the first flag it served.
    exchange_k: bool = False
    #: Only meaningful for the solver: was a time profile supplied (the
    #: ``time_resolved`` flag's equivalent hand-off)?
    source_profile: bool = False
    #: Digest of what the call returned.
    result: str = ""


@dataclass
class Probe:
    """Everything observed at one flag's destination boundary during one run."""

    calls: List[CallRecord] = field(default_factory=list)
    opaque: bool = False
    #: Which flag's boundary this is. ``None`` means :data:`EXCHANGE`, resolved
    #: lazily so the class can be defined before the flag tables exist.
    flag: Optional["Flag"] = None

    @property
    def _flag(self) -> "Flag":
        return self.flag or EXCHANGE

    def of(self, target: str) -> List[CallRecord]:
        return [c for c in self.calls if c.target == target]

    @property
    def solver_calls(self) -> int:
        """How many times the destination actually ran."""
        return len(self.of(self._flag.destination))

    @property
    def exchange_k_calls(self) -> int:
        """Destination calls at which the flag ARRIVED."""
        return sum(1 for c in self.of(self._flag.destination) if c.exchange_k)

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
        """Runs that went through a re-implementation of the destination.

        Non-zero means the work the flag governs WAS done even if the
        destination was never called -- so "it never happened here" would be a
        false statement, and ``UNREACHABLE`` must be withheld.
        """
        return sum(len(self.of(name)) for name in self._flag.alt_destinations)

    @property
    def precursor_calls(self) -> int:
        """Primitives of the governed quantity that ran on this path.

        Withholds ``UNREACHABLE`` for the same reason
        :attr:`alt_solver_calls` does, one step lower down: a run that computed
        an effective diffusivity out of ``brownian`` and ``ecm_mesh`` and never
        reached the two-route law did NOT fail to produce the governed
        quantity. It produced the pre-promotion one. Saying "there was nothing
        here for the flag to change" would be false.
        """
        return sum(len(self.of(name)) for name in self._flag.precursors)

    def solver_digest(self) -> str:
        """Hash of every destination return value, in call order."""
        return digest([c.result for c in self.of(self._flag.destination)])


def _resolve(mod_name: str, attr: str) -> Any:
    import importlib

    return getattr(importlib.import_module(mod_name), attr)


class WatchedTargetMissing(RuntimeError):
    """A :data:`WATCHED` target stopped resolving, and it is not an absent extra.

    :data:`WATCHED` holds STRINGS, so a watched function that is renamed, moved
    or deleted does not break an import -- it just stops being patched. That
    failure is silent and its consequence is the exact inversion this module
    exists to prevent: nothing observes the alternate solver, so
    :attr:`Probe.alt_solver_calls` is zero, so ``UNREACHABLE`` is no longer
    withheld, and the guard publishes "no transport solve happens here" about a
    path that solves transport.

    Exactly one degradation is tolerated, and it is narrow: the target's own
    PACKAGE is not importable at all, which is the ``app/`` optional-extra case
    argued at :data:`WATCHED`. A module that imports fine and no longer has the
    attribute is a RENAME, not an absent extra, and raises this.
    """


def _optional_package_absent(mod_name: str, exc: BaseException) -> bool:
    """Is ``exc`` the module itself not being installed, rather than a rename?

    ``ModuleNotFoundError.name`` is the module that could not be found. It
    counts as "this optional extra is not installed" only when it IS the
    watched module or a package containing it. A ModuleNotFoundError raised by
    a broken dependency INSIDE an installed module is a real failure and must
    stay loud.
    """
    if not isinstance(exc, ModuleNotFoundError):
        return False
    missing = getattr(exc, "name", None)
    if not missing:
        return False
    return mod_name == missing or mod_name.startswith(missing + ".")


def unresolvable_watched(
    flag: Optional["Flag"] = None,
) -> Dict[Tuple[str, str], str]:
    """Watched targets that no longer resolve for a reason that is NOT an
    absent optional package. Empty is the only healthy answer.

    This is :func:`instrument`'s tolerance turned into an assertion. The
    ``except Exception: continue`` in the patch loop is correct for a checkout
    without ``.[app]`` and catastrophic for a rename, and nothing distinguished
    the two until this existed.
    """
    flag = flag or EXCHANGE
    out: Dict[Tuple[str, str], str] = {}
    for mod_name, attr in flag.watched:
        try:
            _resolve(mod_name, attr)
        except Exception as exc:  # noqa: BLE001 - the point is to classify it
            if _optional_package_absent(mod_name, exc):
                continue
            out[(mod_name, attr)] = f"{type(exc).__name__}: {exc}"
    return out


def unwatched_alternates(flag: Optional["Flag"] = None) -> Tuple[str, ...]:
    """Names registered as alternates/precursors that nothing INSTRUMENTS.

    :attr:`Flag.alt_destinations` and :attr:`Flag.precursors` only do their job
    -- withholding ``UNREACHABLE`` -- if :func:`instrument` actually wraps
    them, and it wraps :attr:`Flag.watched`. A name in one table and not the
    other reads as covered and counts nothing: :meth:`Probe.of` returns an
    empty list for a function that was never patched, so the verdict falls
    through to ``UNREACHABLE`` exactly as if the alternate did not exist.
    """
    flag = flag or EXCHANGE
    watched_attrs = {attr for _, attr in flag.watched}
    registered = set(flag.alt_destinations) | set(flag.precursors)
    return tuple(sorted(registered - watched_attrs))


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
        flag = probe._flag
        if name == flag.destination and flag.carrier is None:
            # No carrier: reaching the destination IS the flag being applied.
            exchange_k = True
        elif name == flag.destination and sig is not None:
            try:
                bound = sig.bind(*args, **kwargs)
                exchange_k = bound.arguments.get(flag.carrier, None) is not None
                source_profile = bound.arguments.get("source_profile", None) is not None
            except TypeError:
                # Let the real call raise the real error; record what we can.
                exchange_k = kwargs.get(flag.carrier, None) is not None
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
def instrument(
    watched: Optional[Sequence[Tuple[str, str]]] = None,
    flag: Optional["Flag"] = None,
) -> Iterator[Probe]:
    """Count and fingerprint every call across one flag's destination boundary.

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
    flag = flag or EXCHANGE
    watched = flag.watched if watched is None else watched
    probe = Probe(flag=flag)
    installed: List[Tuple[Any, str, Any, Any]] = []  # (module, attr, original, wrapper)

    try:
        for mod_name, attr in watched:
            try:
                original = _resolve(mod_name, attr)
            except Exception as exc:  # noqa: BLE001 - classified, not swallowed
                # An absent optional package is the ONE tolerated degradation
                # (see WATCHED). Anything else -- above all a module that
                # imports fine and no longer has the attribute -- is a rename,
                # and a renamed watch degrades silently into a false
                # UNREACHABLE. Fail loudly instead.
                if _optional_package_absent(mod_name, exc):
                    continue
                raise WatchedTargetMissing(
                    f"watched target {mod_name}.{attr} did not resolve "
                    f"({type(exc).__name__}: {exc}). WATCHED holds strings, so "
                    "a rename here does not break an import -- it silently "
                    "stops the target being instrumented, and an alternate "
                    "solver that is not instrumented cannot withhold the "
                    "UNREACHABLE verdict. Update WATCHED (and ALT_SOLVERS if "
                    "the name changed); do NOT let this degrade to 'not "
                    "watched'."
                ) from exc
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
    flag: Optional["Flag"] = None,
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

    flag : Flag, optional
        Which switch is being checked. Defaults to :data:`EXCHANGE`.

    Returns
    -------
    Verdict
    """
    flag = flag or EXCHANGE
    with instrument(flag=flag) as p_off:
        r_off = run(False)
    d_off = _run_digest(p_off, r_off)

    if determinism_check:
        with instrument(flag=flag) as p_ctl:
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

    with instrument(flag=flag) as p_on:
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
    dest = flag.destination
    carried = f"{flag.carrier}" if flag.carrier else "the flag"
    route = " or ".join(flag.entry_points) or "a sanctioned entry point"
    if p_off.solver_calls == 0 and p_on.solver_calls == 0:
        if p_on.alt_solver_calls or p_on.precursor_calls:
            if p_on.alt_solver_calls:
                how = (
                    f"{p_on.alt_solver_calls} call(s) went through a "
                    "RE-IMPLEMENTED destination "
                    f"({', '.join(flag.alt_destinations)}) that the flag "
                    "cannot reach at all"
                )
            else:
                how = (
                    f"{p_on.precursor_calls} call(s) went to the RE-IMPLEMENTED "
                    "primitives the governed quantity used to be built out of "
                    f"({', '.join(flag.precursors)}), and the chain stopped "
                    "there"
                )
            return Verdict(
                category=Category.INERT,
                detail=(
                    f"{dest} never ran, but {how}. The work the flag governs "
                    "WAS done and the flag could not reach it. This is INERT, "
                    "not NO_PDE: reporting it as 'nothing happened here' would "
                    "be false -- a number was produced, and it is the flag-OFF "
                    "number whatever the caller set."
                ),
                **common,
            )
        return Verdict(
            category=Category.NO_PDE,
            detail=(
                f"{dest} was never called, with the flag off or on. There is "
                "nothing on this path for the flag to reach (analytic closed "
                "form, or the work happens elsewhere). 'Nothing changed' here "
                "means UNREACHABLE, not robust -- it is no evidence at all "
                "about the physics."
            ),
            **common,
        )

    # (B) The destination ran and the flag never crossed the boundary. THE BUG.
    if p_on.exchange_k_calls == 0:
        if d_off != d_on:
            return Verdict(
                category=Category.PROPAGATED,
                detail=(
                    f"results differ across {p_on.solver_calls} {dest} call(s) "
                    f"even though {carried} was never passed, so the flag is "
                    "reaching the physics by some OTHER route. Real, but not "
                    "the route this checker models -- find it before trusting it."
                ),
                **common,
            )
        built = p_on.exchange_coefficient_calls
        return Verdict(
            category=Category.INERT,
            detail=(
                f"{dest} ran {p_on.solver_calls} time(s) with "
                f"the flag ON and was handed {carried}=None every time"
                + (
                    f" (exchange_coefficient WAS built {built} time(s) and then "
                    "dropped)"
                    if built
                    else " (exchange_coefficient was never even built)"
                    if flag is EXCHANGE
                    else ""
                )
                + ". The flag is a silent no-op on this path: the two runs are "
                "bit-identical because they are the same run. Any 'on vs off' "
                "conclusion drawn here is a comparison of a run to itself. Fix "
                f"by routing through {route}, or by passing {carried} explicitly."
            ),
            **common,
        )

    # (A) The flag crossed the boundary and moved the numbers.
    if d_off != d_on:
        return Verdict(
            category=Category.PROPAGATED,
            detail=(
                f"{carried} reached {dest} on {p_on.exchange_k_calls} of "
                f"{p_on.solver_calls} call(s) and the results changed. An "
                "on-vs-off comparison on this path is real."
            ),
            **common,
        )

    # The flag reached the physics and genuinely did nothing. ROBUST, not unreachable.
    return Verdict(
        category=Category.REACHED_UNCHANGED,
        detail=(
            f"{carried} reached {dest} on {p_on.exchange_k_calls} of "
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
    """What a source file says about how it reaches one flag's destination.

    Read straight off the AST of one file (optionally one function inside it).
    Cheap enough to sweep the whole repo, and it needs no fixtures, but it is
    strictly weaker than a :class:`Verdict`: it sees call SITES, not call PATHS.
    It cannot tell whether a site is reached at runtime, it does not follow calls
    into helper modules, and a dynamic ``exchange_k=`` built behind a conditional
    still counts as passed. Use it to triage; use
    :func:`check_flag_propagation` to conclude.

    The field names are the ``transvascular_exchange`` spelling, kept because
    every existing caller uses them; the generic aliases beside each one are
    what the second flag reads.
    """

    target: str
    #: Line numbers of direct calls to the flag's destination.
    direct_calls: Tuple[int, ...] = ()
    #: Subset of ``direct_calls`` that CARRY the flag. When the flag has no
    #: carrier (:attr:`Flag.carrier` is None) a call IS an arrival, so this
    #: equals ``direct_calls``.
    direct_calls_with_exchange_k: Tuple[int, ...] = ()
    #: Line numbers of calls to a sanctioned entry point (the engine path).
    engine_calls: Tuple[int, ...] = ()
    #: Line numbers of calls to a :attr:`Flag.precursor` -- a primitive that
    #: computes the quantity the flag governs, WITHOUT the flag. Empty for a
    #: carrier-bearing flag, where being a direct caller is already the
    #: statement that this file is on the path.
    precursor_calls: Tuple[int, ...] = ()
    #: Line numbers of calls to a :attr:`Flag.alt_destinations` -- a
    #: RE-IMPLEMENTED destination that answers the same question and that the
    #: flag cannot reach. See :attr:`unreached_alt_calls` for why these are
    #: scanned at all.
    alt_calls: Tuple[int, ...] = ()
    #: Which flag this was scanned for, for the failure message.
    flag: str = "transvascular_exchange"

    # -- generic aliases ---------------------------------------------------- #
    @property
    def destination_calls(self) -> Tuple[int, ...]:
        return self.direct_calls

    @property
    def arrivals(self) -> Tuple[int, ...]:
        """Destination calls at which the flag actually arrived."""
        return self.direct_calls_with_exchange_k

    @property
    def entry_calls(self) -> Tuple[int, ...]:
        return self.engine_calls

    @property
    def bare_direct_calls(self) -> Tuple[int, ...]:
        """Destination calls that do NOT carry the flag. Each one is a line at
        which the flag stops existing.

        Always empty for a carrier-less flag: there, calling the destination is
        the whole of applying the flag, so a call cannot fail to carry it. The
        equivalent defect for those flags is :attr:`stranded_calls`.
        """
        withk = set(self.direct_calls_with_exchange_k)
        return tuple(ln for ln in self.direct_calls if ln not in withk)

    @property
    def stranded_calls(self) -> Tuple[int, ...]:
        """Precursor calls in a file that never reaches the destination.

        The dual of :attr:`bare_direct_calls`. The quantity the flag governs is
        computed here out of the primitives, and the flag cannot touch it -- so
        a number taken from this file is the flag-OFF number no matter what the
        caller set.
        """
        if self.direct_calls or self.engine_calls:
            return ()
        return self.precursor_calls

    @property
    def unreached_alt_calls(self) -> Tuple[int, ...]:
        """Alternate-destination calls in a file that reaches nothing else.

        The third polarity, and the one the sweep was blind to. ``ALT_SOLVERS``
        was only ever consulted DYNAMICALLY -- :attr:`Probe.alt_solver_calls`
        withholds ``UNREACHABLE`` when a re-implemented solver ran. But
        :func:`sweep_repo` looked for the destination, the entry points and the
        precursors and nothing else, so a file whose ONLY route to the physics
        is an alternate solver appeared in no sweep at all. It solved transport,
        the flag could not reach it, and the guard never asked the question --
        which is the same silence as an undeclared ``NO_PDE``, arrived at from
        the other side.

        Suppressed once the file also reaches the destination or an entry
        point, exactly as :attr:`stranded_calls` is: such a file is already on
        the path and already gets a verdict.
        """
        if self.direct_calls or self.engine_calls:
            return ()
        return self.alt_calls

    @property
    def undeclarable(self) -> Tuple[int, ...]:
        """Every line at which this file needs a declaration to stay honest."""
        return tuple(sorted(
            set(self.bare_direct_calls)
            | set(self.stranded_calls)
            | set(self.unreached_alt_calls)
        ))

    @property
    def predicted(self) -> Category:
        """Best static guess at the dynamic category.

        A bare direct call outranks an entry call: if a path can reach the
        destination without the flag, the honest prediction is INERT even when
        some other line in the same file does it correctly. Under-claiming here
        is the safe direction -- the cost of a wrong INERT prediction is one
        investigation, the cost of a wrong PROPAGATED prediction is a published
        false conclusion. A stranded precursor call is ranked the same way and
        for the same reason.
        """
        if self.bare_direct_calls:
            return Category.INERT
        if self.direct_calls_with_exchange_k or self.engine_calls:
            return Category.PROPAGATED
        if self.precursor_calls:
            return Category.INERT
        return Category.UNREACHABLE


def _called_name(node: ast.Call) -> str:
    func = node.func
    if isinstance(func, ast.Name):
        return func.id
    if isinstance(func, ast.Attribute):
        return func.attr
    return ""


def scan_call_sites(
    path: str | Path,
    symbol: Optional[str] = None,
    flag: Optional["Flag"] = None,
) -> CallSites:
    """AST-scan one file (or one function in it) for one flag's call sites.

    ``symbol`` restricts the scan to a single top-level ``def``/``async def``,
    which is what makes a file like ``validation/vessel_profile.py`` -- one
    function on the engine path and another calling the solver directly --
    describable without averaging the two into a meaningless whole-file answer.
    ``molecular/ecm_mesh.py`` is the same shape for :data:`TWO_ROUTE`.

    ``flag`` defaults to :data:`EXCHANGE`, so existing callers are unchanged.
    """
    flag = flag or EXCHANGE
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
    precursor: List[int] = []
    alt: List[int] = []
    for node in ast.walk(root):
        if not isinstance(node, ast.Call):
            continue
        name = _called_name(node)
        if name == flag.destination:
            direct.append(node.lineno)
            # No carrier => arriving IS being called. With a carrier, only a
            # call that actually passes it counts as an arrival.
            if flag.carrier is None or any(
                kw.arg == flag.carrier for kw in node.keywords
            ):
                with_k.append(node.lineno)
        elif name in flag.entry_points:
            engine.append(node.lineno)
        else:
            # Not elif-chained with the precursor arm: for TWO_ROUTE the two
            # tables OVERLAP on purpose (cell_exclusion is both), and a name
            # that fell into only one of them would make the other table's
            # sweep silently narrower than it reads.
            if name in flag.precursors:
                precursor.append(node.lineno)
            if name in flag.alt_destinations:
                alt.append(node.lineno)

    target = str(path) if symbol is None else f"{path}::{symbol}"
    return CallSites(
        target=target,
        direct_calls=tuple(sorted(direct)),
        direct_calls_with_exchange_k=tuple(sorted(with_k)),
        engine_calls=tuple(sorted(engine)),
        precursor_calls=tuple(sorted(precursor)),
        alt_calls=tuple(sorted(alt)),
        flag=flag.name,
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
            "unit_transport calls simulate_drug_transport directly "
            "(kill_calibration.py:844) with no exchange_k. The calibration "
            "constant it fits is therefore an exchange-OFF constant whatever "
            "the config says. (This citation had drifted 163 lines and was "
            "corrected on 2026-08-30 by stale_citations; the claim was still "
            "true, but nothing had re-checked it since it was written.)"
        ),
    ),
    Declaration(
        path="tmesim/validation/immune_sink_impact.py",
        category=Category.INERT,
        reason=(
            "Calls simulate_drug_transport directly at three sites (258, 278, "
            "292) with no exchange_k. Its own per-voxel immune sink is passed by "
            "hand, which is exactly the pattern exchange_k would need and does "
            "not get. (These citations had drifted by 6 and were corrected on "
            "2026-08-30 by stale_citations.)"
        ),
    ),
    Declaration(
        path="tmesim/validation/vessel_profile.py",
        symbol="isolated_vessel_profile",
        category=Category.INERT,
        reason=(
            "Builds its own PDEConfig and calls simulate_drug_transport directly "
            "(vessel_profile.py:729) with no exchange_k, unlike analyze_drug in "
            "the same file. A whole-file verdict would average these two into "
            "something false. (This citation had drifted 83 lines and was "
            "corrected on 2026-08-30 by stale_citations. The stale value had "
            "landed inside the OTHER half of the split -- the half declared "
            "PROPAGATED -- so the evidence for 'this function is INERT' was "
            "pointing at the function that is not.)"
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
    # app/server/reach.py is NOT declared here: it is the HOME of the
    # steady_state_concentration alt solver, so it is registered in ALT_SOLVERS
    # and WATCHED alongside immune_sink.py and moncayo_explant.py. Declaring it
    # instead was tried and was wrong -- the source scan correctly answers
    # NO_PDE for a file that calls no solver because it IS one, and accepting
    # that verdict would have published "no transport solve happens here" about
    # the module that solves the transport steady state.
    # --- the 2026-09-20 verification experiments ------------------------- #
    # Both call the solver BARE, and for both that is the point rather than an
    # oversight: their subject is the solver itself, not the product pipeline
    # the flag travels through. Routing them through simulate_response would
    # destroy what they measure by burying the solve under a reconstruction.
    Declaration(
        path="experiments/mms_transport_verification.py",
        category=Category.INERT,
        reason=(
            "Method of manufactured solutions. It calls the diffusion solver "
            "directly with a source term DERIVED from a chosen exact solution, "
            "so the exchange flag is deliberately absent: an exchange term "
            "would add a sink the manufactured source does not balance, and "
            "the measured order of convergence would then be an order of "
            "convergence to the wrong equation. The verdict it publishes is "
            "about discretisation -- second order for diffusion, first for "
            "upwind advection, no order loss at either boundary -- and it "
            "says in its own 'what this does not show' that it verifies the "
            "SOLVER and not the model. Nothing it reports transfers to a "
            "claim about transvascular exchange."
        ),
    ),
    Declaration(
        path="experiments/conservation_audit.py",
        category=Category.INERT,
        reason=(
            "Physical-invariant audit: mass balance, non-negativity, "
            "monotonicity in dose and window, determinism, scale invariance. "
            "It builds its own transport loops from the private numerics on "
            "purpose, because an invariant must be checked on the operator "
            "rather than on whatever the product pipeline happened to "
            "configure -- a budget that closes only under the shipped default "
            "is not a conservation result. Exchange is a SOURCE term, so "
            "switching it on would change the budget the audit is balancing "
            "and every invariant would have to be restated around it. Its "
            "findings are about the operator as discretised and carry no "
            "exchange claim."
        ),
    ),
    Declaration(
        path="experiments/steady_state_reach.py",
        category=Category.INERT,
        reason=(
            "march() calls simulate_drug_transport BARE "
            "(steady_state_reach.py:233) with a hand-built PDEConfig and no "
            "exchange_k. That is the load-bearing fact and not an incidental "
            "one: the finite-window arm is the object the steady state is "
            "measured AGAINST, so both arms of that comparison are "
            "exchange-OFF and the comparison is internally consistent -- but "
            "neither arm licenses a claim about exchange. The steady solve it "
            "drives takes exchange_k explicitly and is passed None."
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
    Exemption(
        path="tmesim/validation/moncayo_explant.py",
        reason=(
            "The module that DEFINES the bath_diffusion alternate solver "
            "(line 528) and then exercises it. Its own calls are the solver "
            "being run by its home, which is the same situation "
            "physics/diffusion.py is in one level up; the question this sweep "
            "asks -- 'did a caller reach the physics without the flag' -- is "
            "not meaningful at the definition site. Callers ELSEWHERE are not "
            "exempt and are still swept."
        ),
    ),
    Exemption(
        path="app/server/reach.py",
        reason=(
            "The module that DEFINES the steady_state_concentration alternate "
            "solver (line 559) and then drives it from its own request handler "
            "(line 952). Same standing as diffusion.py and moncayo_explant.py: "
            "the definition site is where a flag ARRIVES, not a place that can "
            "fail to carry it, and this solver takes exchange_k explicitly. It "
            "is exempted rather than declared for the reason written against "
            "DECLARATIONS: a Declaration here would have to name a category, "
            "and NO_PDE -- 'no transport solve happens here' -- is exactly the "
            "false statement to publish about the module that solves the "
            "transport steady state."
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


#: Directories swept by :func:`sweep_repo`.
SWEPT_ROOTS: Tuple[str, ...] = ("tmesim", "experiments", "scripts", "app")


# --------------------------------------------------------------------------- #
# The flag itself
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Flag:
    """One physics switch, the destination it must reach, and its own tables.

    Everything above this line is written against a :class:`Flag` rather than
    against ``transvascular_exchange``, so a second switch of the same shape
    gets the same guard instead of a second copy of it.
    """

    #: The switch as a person would name it in a config or a signature.
    name: str
    #: The function the flag has to reach to mean anything.
    destination: str
    #: The keyword argument that CARRIES the flag across that boundary, or
    #: ``None`` when there is no carrier and the call itself is the arrival.
    #: This one field is the whole difference between the two flags -- see the
    #: module docstring.
    carrier: Optional[str]
    #: Sanctioned routes: calling one of these is how the flag is SUPPOSED to
    #: reach the destination, and a file that does so is on the honest path.
    entry_points: Tuple[str, ...] = ()
    #: Primitives that compute the quantity the flag governs WITHOUT it. Only
    #: meaningful for a carrier-less flag; see :attr:`CallSites.stranded_calls`.
    precursors: Tuple[str, ...] = ()
    #: ``(module, attribute)`` pairs :func:`instrument` patches for this flag.
    watched: Tuple[Tuple[str, str], ...] = ()
    #: Re-implementations that answer the same question the destination does,
    #: and that the flag cannot reach. Their presence WITHHOLDS the
    #: ``UNREACHABLE`` verdict -- see :data:`ALT_SOLVERS`.
    alt_destinations: Tuple[str, ...] = ()
    declarations: Tuple[Declaration, ...] = ()
    exemptions: Tuple[Exemption, ...] = ()
    wiring_points: Tuple[WiringPoint, ...] = ()
    #: Roots swept for undeclared paths.
    swept_roots: Tuple[str, ...] = SWEPT_ROOTS
    #: Roots in which a file that reaches the destination must be a declared
    #: :class:`WiringPoint` (or an :class:`Exemption`).
    #:
    #: THIS IS THE ONE THING THAT DID NOT GENERALISE, and it is data, not a
    #: special case in the logic.
    #:
    #: For ``transvascular_exchange`` the scope is the whole sweep, and the
    #: argument for that is in :data:`WIRING_POINTS`: passing ``exchange_k`` by
    #: hand is how an experiment FAKES the production path, so an undeclared
    #: injection is an experiment that may be measuring something other than
    #: what it says. That danger is real and it lives mostly in
    #: ``experiments/``.
    #:
    #: For ``two_route`` that danger does not exist. There is no carrier to
    #: inject, so calling ``two_route_diffusivity`` in a script is not faking a
    #: path -- it is calling a pure function with explicit arguments, which is
    #: how every one of this law's probe experiments is written and how the
    #: next one will be. Requiring a WiringPoint per probe would charge a
    #: written justification for the ordinary way of asking the law a question,
    #: and a rule that fires on ordinary work is a rule people re-stamp without
    #: reading. So the scope is narrowed to the SHIPPED PACKAGE, where the
    #: invariant that actually matters lives: exactly one place inside
    #: ``tmesim/`` decides whether a compiled drug gets the cell term.
    #:
    #: Narrowing the scope costs nothing that was being guarded, because the
    #: exhaustive sweep for STRANDED paths still covers every root. An
    #: ``experiments/`` script that computes a diffusivity without the cell
    #: term is still a failure; only "called the law directly" is untaxed.
    wiring_scope: Tuple[str, ...] = SWEPT_ROOTS
    #: Why this flag is watched at all. Required, like every other reason here.
    reason: str = ""

    def __post_init__(self) -> None:
        if len(self.reason.strip()) < 25:
            raise ValueError(
                f"{self.name}: a Flag needs a written reason (>=25 chars) "
                "saying what goes silently wrong when it does not propagate."
            )
        if self.carrier is None and not self.precursors:
            raise ValueError(
                f"{self.name}: a carrier-less flag has no bare calls to find, "
                "so without precursors this Flag can never report the defect "
                "it exists to report. Declare the primitives that compute the "
                "quantity it governs."
            )

    def in_wiring_scope(self, path: str) -> bool:
        norm = str(path).replace("\\", "/")
        return any(norm == r or norm.startswith(r + "/") for r in self.wiring_scope)


def wiring_points(
    production_only: bool = False, flag: Optional[Flag] = None
) -> frozenset:
    """Paths permitted to carry the flag; the production one if asked."""
    flag = flag or EXCHANGE
    return frozenset(
        w.path for w in flag.wiring_points if w.production or not production_only
    )


def declaration_for(
    path: str, symbol: Optional[str] = None, flag: Optional[Flag] = None
) -> Declaration:
    """Look up a declaration, or raise. There is no default."""
    flag = flag or EXCHANGE
    name = path if symbol is None else f"{path}::{symbol}"
    for d in flag.declarations:
        if d.name == name:
            return d
    raise LookupError(
        f"{name} has no {flag.name} propagation declaration. Add one to the "
        f"declarations of tmesim.validation.propagation.{flag.name.upper()} "
        "with a written reason -- an undeclared path is a failure, not a pass."
    )


def sweep_repo(root: str | Path = ".", flag: Optional[Flag] = None) -> Dict[str, CallSites]:
    """Static-scan the repo for every file on one flag's path.

    Returns ``{repo-relative path: CallSites}``. Cheap (AST only) and the reason
    the guard does not depend on someone remembering to register a new module:
    a file that starts calling the destination -- or that starts computing the
    governed quantity out of the primitives -- shows up here whether or not
    anyone thought about it.
    """
    flag = flag or EXCHANGE
    root = Path(root)
    found: Dict[str, CallSites] = {}
    for base in flag.swept_roots:
        base_dir = root / base
        if not base_dir.is_dir():
            continue
        for py in sorted(base_dir.rglob("*.py")):
            try:
                cs = scan_call_sites(py, flag=flag)
            except (SyntaxError, UnicodeDecodeError):  # pragma: no cover
                continue
            if (cs.direct_calls or cs.engine_calls or cs.precursor_calls
                    or cs.alt_calls):
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


def undeclared_sites(
    root: str | Path = ".", flag: Optional[Flag] = None
) -> Dict[str, Tuple[int, ...]]:
    """Undeclared silent-no-op LINES, per file. See :func:`undeclared_paths`.

    Line-aware on purpose. The path-keyed version of this check treated ANY
    declaration for a file as covering the whole file, which meant a
    per-function declaration -- the mechanism that exists precisely because
    ``vessel_profile.py`` has one safe half and one unsafe half -- also silently
    licensed every OTHER function in that file, including ones written later.
    A new bare call added beside a declared one produced no failure at all.

    A ``symbol`` declaration now covers only the lines inside that symbol; the
    residue is reported.
    """
    flag = flag or EXCHANGE
    root = Path(root)
    exempt = {e.path for e in flag.exemptions}
    whole_file = {d.path for d in flag.declarations if d.symbol is None}
    by_symbol: Dict[str, List[str]] = {}
    for d in flag.declarations:
        if d.symbol is not None:
            by_symbol.setdefault(d.path, []).append(d.symbol)

    out: Dict[str, Tuple[int, ...]] = {}
    for path, cs in sweep_repo(root, flag=flag).items():
        if path in exempt or path in whole_file:
            continue
        residue = set(cs.undeclarable)
        if not residue:
            continue
        for symbol in by_symbol.get(path, ()):
            # A declared symbol that has since been renamed away is caught by
            # test_declared_path_lands_in_its_declared_category, which scans it
            # directly and raises. Here it simply covers nothing, which is the
            # safe direction.
            try:
                sub = scan_call_sites(root / path, symbol, flag=flag)
            except LookupError:
                continue
            residue -= set(sub.undeclarable)
        if residue:
            out[path] = tuple(sorted(residue))
    return out


def undeclared_paths(root: str | Path = ".", flag: Optional[Flag] = None) -> List[str]:
    """Files on which the flag is a silent no-op and nobody has written it down.

    Three shapes, one per polarity (see :class:`CallSites`):

    * BARE -- the destination is called without the carrier.
    * STRANDED -- the governed quantity is computed out of the precursors and
      the destination is never reached.
    * UNREACHED ALTERNATE -- the only route to the physics is a RE-IMPLEMENTED
      destination the flag cannot reach.

    Either way the run completes, the numbers are plausible, and they are the
    flag-OFF numbers whatever the caller set. The correct response is to add a
    declaration (or fix the call), never to widen the exemption list.
    """
    return sorted(undeclared_sites(root, flag=flag))


# --------------------------------------------------------------------------- #
# Line citations inside a written reason
# --------------------------------------------------------------------------- #
#: The three ways a reason in this module cites a line. Deliberately narrow --
#: a looser pattern would harvest digits out of ``1.554447793e+09`` and
#: ``R^2 = 0.022`` and the check would fail on prose instead of on drift.
_CITE_PATH_LINE = re.compile(r"\b([\w./-]+\.py):(\d+)\b")
_CITE_LINE_WORD = re.compile(r"\blines?\s+(\d+(?:\s*,\s*\d+)*)")
_CITE_TUPLE = re.compile(r"\((\d+(?:\s*,\s*\d+)+)\)")


def cited_lines(reason: str) -> Tuple[Tuple[Optional[str], int], ...]:
    """Every ``(file-or-None, line)`` a written reason points at.

    Recognises ``pk_reality.py:218``, ``(line 307)`` and ``(252, 272, 286)``,
    which are the three spellings the tables in this module actually use.
    """
    out: List[Tuple[Optional[str], int]] = []
    for m in _CITE_PATH_LINE.finditer(reason):
        out.append((m.group(1), int(m.group(2))))
    for rx in (_CITE_LINE_WORD, _CITE_TUPLE):
        for m in rx.finditer(reason):
            for chunk in m.group(1).split(","):
                out.append((None, int(chunk.strip())))
    return tuple(sorted(set(out), key=lambda t: (t[0] or "", t[1])))


def recognised_lines(path: str | Path, flag: Optional[Flag] = None) -> frozenset:
    """Lines in ``path`` this flag's scanner can honestly point at.

    Call sites for the flag (destination, entry points, precursors, alternate
    destinations) plus the ``def`` lines of those same names, so that a table
    entry may cite either "the call that strands the flag" or "the function
    that defines the alternate solver".
    """
    flag = flag or EXCHANGE
    path = Path(path)
    cs = scan_call_sites(path, flag=flag)
    lines = (set(cs.direct_calls) | set(cs.engine_calls)
             | set(cs.precursor_calls) | set(cs.alt_calls))
    names = {flag.destination, *flag.entry_points, *flag.precursors,
             *flag.alt_destinations}
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    for node in ast.walk(tree):
        if (isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
                and node.name in names):
            lines.add(node.lineno)
    return frozenset(lines)


def stale_citations(
    root: str | Path = ".", flag: Optional[Flag] = None
) -> Dict[str, Tuple[str, ...]]:
    """Table entries whose written reason cites a line that has MOVED.

    The reason is the whole product of this module: ``NO_PDE`` and ``INERT``
    are claims someone signed, and the signature is the sentence saying WHERE
    the flag stops existing. A reason that names a line number is making a
    checkable statement, and until this function existed nothing checked it --
    so the line drifted with every edit above it and the entry went on reading
    as though it had been verified. Three of them had already drifted when this
    was written, and one had drifted onto the OTHER half of a split file: the
    reason for "this function is INERT" pointed into the function declared
    PROPAGATED.

    Returns ``{table entry name: (problem, ...)}``. Empty is the healthy answer.
    """
    flag = flag or EXCHANGE
    root = Path(root)
    entries: List[Tuple[str, str, str]] = [
        (d.name, d.path, d.reason) for d in flag.declarations
    ]
    entries += [(f"exemption:{e.path}", e.path, e.reason) for e in flag.exemptions]
    entries += [(f"wiring:{w.path}", w.path, w.reason) for w in flag.wiring_points]

    cache: Dict[str, frozenset] = {}
    out: Dict[str, Tuple[str, ...]] = {}
    for name, path, reason in entries:
        problems: List[str] = []
        for named_file, line in cited_lines(reason):
            if named_file is not None and Path(named_file).name != Path(path).name:
                problems.append(
                    f"cites {named_file}:{line}, which is not the file this "
                    f"entry is about ({path}). Cite lines in the declared file."
                )
                continue
            if path not in cache:
                f = root / path
                cache[path] = (
                    recognised_lines(f, flag=flag) if f.is_file() else frozenset()
                )
            if line not in cache[path]:
                known = sorted(cache[path])
                problems.append(
                    f"cites line {line}, which is not a {flag.name} call site "
                    f"or definition in {path}. The lines the scanner "
                    f"recognises there are {known or 'none at all'}."
                )
        if problems:
            out[name] = out.get(name, ()) + tuple(problems)
    return out


def undeclared_direct_callers(root: str | Path = ".") -> List[str]:
    """:func:`undeclared_paths` for :data:`EXCHANGE`. Kept as the name every
    existing caller and failure message uses."""
    return undeclared_paths(root, flag=EXCHANGE)


def undeclared_wiring(root: str | Path = ".", flag: Optional[Flag] = None) -> List[str]:
    """Files inside :attr:`Flag.wiring_scope` that reach the destination WITH
    the flag and are not declared as a :class:`WiringPoint` (nor exempted).

    Every such file is a second place that decides whether the physics is
    active, which is how a flag ends up silently on in one path and off in
    another.
    """
    flag = flag or EXCHANGE
    permitted = wiring_points(flag=flag) | {e.path for e in flag.exemptions}
    return sorted(
        path
        for path, cs in sweep_repo(root, flag=flag).items()
        if cs.arrivals and flag.in_wiring_scope(path) and path not in permitted
    )


# --------------------------------------------------------------------------- #
# Flag 2: molecular.compile_drug.build_profile(two_route=...)
#
# The second switch of the identical shape, added by fd2a9e7 with no table of
# its own. `build_profile(two_route: bool = _TWO_ROUTE_DEFAULT)` decides whether
# a compiled drug's interstitial diffusivity is multiplied by the two-route cell
# term. Everything below is the same discipline the exchange flag already had,
# in the dual polarity described in the module docstring.
# --------------------------------------------------------------------------- #

#: What :func:`instrument` patches for :data:`TWO_ROUTE`. The first is THE
#: destination. The rest are the precursors, counted so a diagnosis can tell
#: "the interstitial leg was never even computed" from "it was computed and the
#: cell term was never applied to it".
TWO_ROUTE_WATCHED: Tuple[Tuple[str, str], ...] = (
    ("tmesim.molecular.two_route", "two_route_diffusivity"),
    ("tmesim.molecular.brownian", "simulate_effective_diffusivity"),
    ("tmesim.molecular.ecm_mesh", "ogston_renkin_hindrance"),
    # Alternate destinations -- see TWO_ROUTE_ALT. BOTH of them: tissue_hindrance
    # was listed in TWO_ROUTE_ALT and TWO_ROUTE_PRECURSORS but was never watched,
    # so Probe.of("tissue_hindrance") returned [] however many times it ran, and
    # neither alt_solver_calls nor precursor_calls could count it. It only failed
    # to produce a false UNREACHABLE because it happens to call cell_exclusion
    # internally; inlining that would have flipped the verdict. Found by
    # unwatched_alternates() on 2026-08-30.
    ("tmesim.molecular.ecm_mesh", "cell_exclusion"),
    ("tmesim.molecular.ecm_mesh", "tissue_hindrance"),
)

#: The IMPERMEABLE-cell law: it answers the same question the two-route law
#: answers (how much do cells hinder a solute?) and it cannot see the flag,
#: because it takes no ``logp`` and no ``hbd`` -- it is identical for every drug.
#: ``two_route`` reduces to it exactly at zero membrane permeability.
#:
#: It is registered here for the same reason ``ALT_SOLVERS`` is: a path that
#: computes cell hindrance through ``cell_exclusion`` and never reaches
#: ``two_route_diffusivity`` would otherwise be reported UNREACHABLE -- "the
#: flag has nothing to reach here" -- when the truth is that the governed
#: quantity WAS computed and the flag could not reach it. That is INERT.
TWO_ROUTE_ALT: Tuple[str, ...] = (
    "cell_exclusion",
    "tissue_hindrance",
)

#: The primitives the pre-promotion diffusivity chain was built out of. A file
#: that calls these and never reaches ``two_route_diffusivity`` (directly or
#: through ``build_profile``) produces the MOLECULAR-WEIGHT-ONLY number that
#: ``results/diffusivity_falsification.md`` scored at R^2 = 0.022 -- and it
#: produces it whatever anyone sets ``two_route`` to.
TWO_ROUTE_PRECURSORS: Tuple[str, ...] = (
    "simulate_effective_diffusivity",
    "ogston_renkin_hindrance",
    "cell_exclusion",
    "tissue_hindrance",
)

TWO_ROUTE_DECLARATIONS: Tuple[Declaration, ...] = (
    # --- (A) the compile path: the flag propagates ---------------------------
    Declaration(
        path="tmesim/molecular/compile_drug.py",
        symbol="build_profile",
        category=Category.PROPAGATED,
        reason=(
            "The one and only place that reads the two_route parameter, calls "
            "two_route.two_route_diffusivity with the logp and hbd descriptors, "
            "and multiplies the interstitial D by the returned factor. Every "
            "shipped DrugParams.diffusivity in this repo comes through here, "
            "and alpha moved 148.60815765445545 -> 22.223477711172563 because "
            "of it."
        ),
    ),
    Declaration(
        path="tmesim/molecular/ecm_mesh.py",
        symbol="tissue_hindrance_two_route",
        category=Category.PROPAGATED,
        reason=(
            "Reaches two_route.two_route_diffusivity, so the cell term is "
            "genuinely applied. This is the half of ecm_mesh.py that is on the "
            "law; the other half is not, which is why this file is declared "
            "per-function exactly as vessel_profile.py is for the exchange flag."
        ),
    ),
    # --- (B) stranded: the governed quantity is computed without the flag ----
    Declaration(
        path="tmesim/molecular/ecm_mesh.py",
        symbol="tissue_hindrance",
        category=Category.INERT,
        reason=(
            "Composes ogston_renkin_hindrance with cell_exclusion, the "
            "IMPERMEABLE-cell law, and never reaches two_route_diffusivity. Its "
            "number is identical for two drugs of the same size, which is the "
            "decisive failure results/diffusivity_falsification.md measured. It "
            "is kept and untouched on purpose as the pre-promotion comparison, "
            "and it must never be mistaken for the shipping law."
        ),
    ),
    Declaration(
        path="tmesim/viz/explainer.py",
        category=Category.INERT,
        reason=(
            "The first thing this sweep caught, and it is a real split rather "
            "than a false alarm. explainer.py calls ecm_mesh."
            "ogston_renkin_hindrance directly (explainer.py:278) to publish "
            "mesh_partition and mesh_hindrance, which are the ECM-ONLY leg and "
            "carry no cell term. The tissue_D column beside them comes from "
            "calibration.calibrated_drugparams, which does route through "
            "compile_drug.molecule_to_drugparams and so IS the two-route "
            "number. The two columns are therefore products of different laws: "
            "free_D x mesh_hindrance does NOT reproduce tissue_D, and nothing "
            "in the figure says so. Labelled correctly today; declared here so "
            "that if the mesh leg is ever presented as the tissue answer, the "
            "sweep is already holding the distinction."
        ),
    ),
    # --- (C) unreachable: no diffusivity of the governed kind exists here ----
    Declaration(
        path="tmesim/molecular/stokes_einstein.py",
        category=Category.UNREACHABLE,
        reason=(
            "Stokes-Einstein in a DILUTE fluid: there are no cells and no ECM "
            "in the model at all, so a cell-permeability term has nothing to "
            "reach. free_diffusivity is bit-identical under either setting, and "
            "that is a statement about this module's scope -- NOT evidence that "
            "the cell term is negligible."
        ),
    ),
)

TWO_ROUTE_EXEMPTIONS: Tuple[Exemption, ...] = (
    Exemption(
        path="tmesim/molecular/two_route.py",
        reason=(
            "The module that DEFINES the law, plus its own __main__ self-test. "
            "The destination is where a flag has to arrive, so 'does the flag "
            "reach it from here' is not a meaningful question; these calls are "
            "the law being exercised directly with explicit arguments."
        ),
    ),
    Exemption(
        path="tmesim/molecular/brownian.py",
        reason=(
            "The module that DEFINES simulate_effective_diffusivity, plus its "
            "own __main__ self-test. It computes the ECM micro-tortuosity leg "
            "and is SUPPOSED to stop there -- the cell term is applied one "
            "level up by compile_drug. A stranded verdict here would be the "
            "checker complaining that a component is a component."
        ),
    ),
)

#: Every file inside the SHIPPED PACKAGE that reaches the two-route law, and why.
#:
#: The invariant is the same one WIRING_POINTS protects and it is just as
#: strict: there is EXACTLY ONE production wiring point. A second place inside
#: tmesim/ deciding whether a compiled drug gets the cell term is how the same
#: drug ends up with two different diffusivities depending on which module
#: compiled it.
#:
#: The scope is tmesim/ rather than the whole sweep. The argument is written out
#: at Flag.wiring_scope; in one line, this flag has no carrier to inject, so an
#: experiments/ script calling the law directly is not faking a production path
#: and charging it a written justification would tax ordinary work.
TWO_ROUTE_WIRING_POINTS: Tuple[WiringPoint, ...] = (
    WiringPoint(
        path="tmesim/molecular/compile_drug.py",
        production=True,
        reason=(
            "The shipped compile path, and the only place inside tmesim/ that "
            "turns the two_route parameter into a call to "
            "two_route.two_route_diffusivity. Everything a user runs -- every "
            "DrugParams the engine is handed -- reaches the law through here."
        ),
    ),
    WiringPoint(
        path="tmesim/molecular/ecm_mesh.py",
        production=False,
        reason=(
            "tissue_hindrance_two_route, a convenience wrapper for analysis. "
            "Its own docstring says it is NOT the shipping path: it recomputes "
            "the interstitial leg analytically from Ogston-Renkin instead of "
            "from the Brownian sim compile_drug uses, so its result is NOT "
            "bit-identical to a compile and must never be quoted as one. It is "
            "declared rather than exempted precisely because the two look alike."
        ),
    ),
)

TWO_ROUTE = Flag(
    name="two_route",
    destination="two_route_diffusivity",
    # No carrier. The flag's entire content is whether the law is invoked, so a
    # call IS an arrival and there is no such thing as a bare one.
    carrier=None,
    entry_points=("build_profile", "molecule_to_drugparams"),
    precursors=TWO_ROUTE_PRECURSORS,
    watched=TWO_ROUTE_WATCHED,
    alt_destinations=TWO_ROUTE_ALT,
    declarations=TWO_ROUTE_DECLARATIONS,
    exemptions=TWO_ROUTE_EXEMPTIONS,
    wiring_points=TWO_ROUTE_WIRING_POINTS,
    wiring_scope=("tmesim",),
    reason=(
        "compile_drug.build_profile(two_route=...) decides whether a compiled "
        "drug's diffusivity carries the two-route cell term. A path that "
        "computes an effective diffusivity out of brownian/ecm_mesh primitives "
        "without reaching two_route_diffusivity gets the pre-promotion "
        "MOLECULAR-WEIGHT-ONLY number -- the model "
        "results/diffusivity_falsification.md scored at R^2 = 0.022 against 71 "
        "compounds -- and gets it whatever anyone set the flag to. The run "
        "completes, the number is plausible, and it is flat by construction at "
        "matched molecular weight."
    ),
)

EXCHANGE = Flag(
    name="transvascular_exchange",
    destination=SOLVER,
    carrier="exchange_k",
    entry_points=("simulate_response",),
    # None. For a carrier-bearing flag, being a direct caller of the solver is
    # already the statement that this file is on the path; there is no separate
    # "computed it out of primitives" shape to look for at a CALL SITE.
    # (reimplementation_candidates covers that shape by IMPORT instead.)
    precursors=(),
    watched=WATCHED,
    alt_destinations=ALT_SOLVERS,
    declarations=DECLARATIONS,
    exemptions=EXEMPTIONS,
    wiring_points=WIRING_POINTS,
    wiring_scope=SWEPT_ROOTS,
    reason=(
        "PDEConfig.transvascular_exchange only means something once something "
        "turns it into an exchange_k field and hands it to the solver. On every "
        "path that calls the solver directly it is a silent no-op and the two "
        "runs are BIT-IDENTICAL, measured at auc_sum 1.554447793e+09 either way "
        "against 6.86608671e+08 when exchange_k was actually passed."
    ),
)

#: Every flag under guard. A new physics switch of this shape belongs here, not
#: in a second copy of this module.
FLAGS: Tuple[Flag, ...] = (EXCHANGE, TWO_ROUTE)


def flag_named(name: str) -> Flag:
    for f in FLAGS:
        if f.name == name:
            return f
    raise LookupError(
        f"no Flag named {name!r}; known flags are "
        + ", ".join(f.name for f in FLAGS)
    )
