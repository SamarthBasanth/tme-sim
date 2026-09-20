"""The guard against a document that says something the code stopped doing.

``results/CLAIM_INVENTORY.md`` says it outright: ``results/*.md`` **is** this
product's claim surface. It is the thing a buyer reads, and it is the only part
of the repo with no CI.

The recurrence
--------------
When the two-route promotion (``fd2a9e7``, 2026-08-14) landed it was disciplined
about staleness in CODE -- it re-pinned the reference JSON, preserved every
superseded value under a dated key, wrote a ``_repin_history``. That discipline
was never extended to the DOCUMENTS. The engine has tests, so its stale state
fails loudly; the documents have prose, so theirs fails silently, and nine of
them had to be corrected by hand. A sample of what had gone false:

    results/logp_provenance.md   "logP is not an input to transport at any point"
    results/two_route_transport.md  "Nothing is wired into compile_drug.py ...
                                     calibrate_alpha() still returns 148.60815765445545"
    docs/QUICKSTART.md           "default-OFF and wired into nothing"
    results/logp_provenance.md   resolve_logp "has no callers"

Every one of those is a **negative reachability or identity assertion that
happens to be typed in English**. They are decidable. This file decides them.

How a claim gets guarded
------------------------
A document DECLARES its claim in :data:`CLAIMS`. Nothing is inferred by
scanning prose, and that is deliberate rather than lazy:

* This repo's house convention is to CORRECT a stale document with a dated
  banner and leave the superseded sentence standing underneath it. So
  ``results/two_route_transport.md`` contains, on purpose, the literal words
  ``calibrate_alpha() still returns 148.60815765445545`` -- inside a quotation
  of what is now false. A text scan would fire on the correction.
* Most files quoting ``alpha = 148.6`` are historical run records, and they are
  RIGHT: that was alpha when the run happened. Past tense and present tense are
  not statically separable.

So a claim is opt-in, and the price of opting in is an :attr:`Claim.anchor` --
a short load-bearing substring that IS the claim. The anchor is checked to still
be present, so a declaration cannot outlive the sentence it speaks for, and a
document that drops a claim must drop its declaration in the same edit.

A failure names BOTH sides:

    results/logp_provenance.md asserts `logp` is not read on the path to
    two_route_diffusivity; tmesim/molecular/compile_drug.py:193 reads it.

Five verbs
----------
``assert_call`` / ``assert_no_call``
    Does module A reach module B at all? (:class:`CallEdge`)
``assert_const``
    Does this function still return this literal? (:class:`Const`) The highest
    yield per line here: the stale ``alpha`` literal appeared in at least four
    documents and was wrong in all of them the moment ``fd2a9e7`` landed.
``assert_default``
    Is this switch still shipped in this position, in the constant, at runtime,
    AND in the signature that is supposed to read it? (:class:`Default`)
``assert_no_callers``
    Is this function still inert? (:class:`NoCallers`) This one would have
    fired on ``resolve_logp`` BEFORE someone wired it up trusting a document
    that said it was inert.
``assert_descriptor`` / ``assert_no_descriptor``
    Does this molecular descriptor key reach this physics law? (:class:`Descriptor`)

Both polarities exist for three of the five because the corrections were
polarity flips: the document that used to say "nothing is wired in" now says
"this is wired in", and the guard has to be able to hold it to either.

Why the reason is mandatory
---------------------------
The known failure mode of a guard like this is **reflexive re-stamping** --
someone silences it without reading. The mitigation is the one this repo already
uses in ``tmesim/validation/propagation.py``: a written reason, validated for
length at construction. The friction is the mechanism. Editing a number in
:data:`CLAIMS` is one keystroke; editing it while a 25-character justification
sits beside it is a decision.

WHAT THIS CANNOT DO -- read this before trusting it
---------------------------------------------------
This guard reaches REACHABILITY and IDENTITY. It does not reach numbers, and it
does not reach meaning.

1. **Numbers are out of scope.** ``R² = 0.718`` cannot be checked this way; the
   only way to know it is still 0.718 is to re-run the scoring, which is a CI
   budget problem and not a static one. :class:`Const` is the narrow exception,
   and only because ``calibrate_alpha()`` is a closed-form solve that costs
   0.4 s. Nothing here re-runs a Brownian sweep, a reconstruction, or a PDE.

2. **Interpretations rot while every number stays correct, and NO static check
   reaches that.** ``results/external_anchor_survey.md`` said the engine
   "hinders diffusion 7× too weakly". The 7× was an accurate number. The
   DIAGNOSIS attached to it was dissolved by a later commit, and the sentence
   went false without a single digit in it changing. This guard would have been
   green throughout. That failure mode is the more common one and the more
   dangerous one, and the only thing that catches it is a person re-reading the
   document.

3. It is opt-in, so it covers exactly the claims someone declared. An undeclared
   document is not guarded, and green here is not a statement about it.

A guard is not coverage. Four decidable claims held to their word is what this
file is; it is not a certificate that ``results/`` is true.

Cost
----
Sections 1-3 parse source without executing it. :class:`Const` calls
``calibrate_alpha()`` once (0.4 s, closed form, memoised across the module) and
:class:`Default` imports two modules. No reconstruction, no PDE, no sweep. The
whole file is about a second.
"""
from __future__ import annotations

import ast
import importlib
import inspect
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

import pytest

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))


# =========================================================================== #
# A source index: enough import resolution to answer "who calls what"
# =========================================================================== #
@dataclass(frozen=True)
class Module:
    """One parsed module, with its import bindings resolved to dotted names.

    Resolving imports is what separates this from grep. ``_two_route.
    two_route_diffusivity(...)`` and ``two_route_diffusivity(...)`` and
    ``two_route.two_route_diffusivity(...)`` are the same edge written three
    ways, and a claim about reachability has to see all three. A mention inside
    a docstring or a string literal is NOT an edge, and the AST never confuses
    the two -- which matters here, because these modules discuss their own call
    graph in prose at length.
    """

    dotted: str
    path: Path
    tree: ast.Module
    #: local name -> dotted module  (``from . import binding`` -> binding)
    module_aliases: Dict[str, str]
    #: local name -> (dotted module, symbol)  (``from .x import f`` -> f)
    symbol_aliases: Dict[str, Tuple[str, str]]
    #: name -> def node, at any nesting depth; first definition wins
    functions: Dict[str, ast.AST]

    @property
    def rel(self) -> str:
        return str(self.path.relative_to(REPO)).replace("\\", "/")


def _dotted_to_path(dotted: str) -> Optional[Path]:
    base = REPO / Path(*dotted.split("."))
    for candidate in (base.with_suffix(".py"), base / "__init__.py"):
        if candidate.is_file():
            return candidate
    return None


def _resolve_relative(current: str, level: int, module: Optional[str]) -> str:
    """``from ..core import X`` inside ``tmesim.molecular.compile_drug``."""
    parts = current.split(".")
    # A module's own name is level 1's starting point (its package).
    base = parts[: len(parts) - level]
    if module:
        base = base + module.split(".")
    return ".".join(base)


_MODULE_CACHE: Dict[str, Optional[Module]] = {}


def load_module(dotted: str) -> Optional[Module]:
    if dotted in _MODULE_CACHE:
        return _MODULE_CACHE[dotted]
    path = _dotted_to_path(dotted)
    if path is None:
        _MODULE_CACHE[dotted] = None
        return None
    mod = _parse(dotted, path)
    _MODULE_CACHE[dotted] = mod
    return mod


def _parse(dotted: str, path: Path) -> Module:
    return parse_source(dotted, path.read_text(encoding="utf-8"), path)


def parse_source(dotted: str, source: str, path: Path) -> Module:
    """Index one module from SOURCE TEXT.

    Split out from :func:`_parse` so the proof tests can exercise the import
    resolver on synthetic modules without writing a file into ``tmesim/``. A
    file dropped there, even for a millisecond, is visible to
    ``propagation.sweep_repo`` and would make a concurrent run of
    ``test_propagation.py`` fail for a reason that has nothing to do with it.
    """
    tree = ast.parse(source, filename=str(path))
    module_aliases: Dict[str, str] = {}
    symbol_aliases: Dict[str, Tuple[str, str]] = {}
    functions: Dict[str, ast.AST] = {}

    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                local = alias.asname or alias.name.split(".")[0]
                target = alias.name if alias.asname else alias.name.split(".")[0]
                module_aliases[local] = target
        elif isinstance(node, ast.ImportFrom):
            base = (
                _resolve_relative(dotted, node.level, node.module)
                if node.level
                else (node.module or "")
            )
            for alias in node.names:
                local = alias.asname or alias.name
                child = f"{base}.{alias.name}" if base else alias.name
                if _dotted_to_path(child) is not None:
                    # ``from . import two_route`` -- a MODULE, not a symbol.
                    module_aliases[local] = child
                else:
                    symbol_aliases[local] = (base, alias.name)
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            functions.setdefault(node.name, node)

    return Module(
        dotted=dotted,
        path=path,
        tree=tree,
        module_aliases=module_aliases,
        symbol_aliases=symbol_aliases,
        functions=functions,
    )


def call_target(node: ast.Call, mod: Module) -> Optional[Tuple[str, str]]:
    """``(dotted module, symbol)`` this call resolves to, or None.

    None means "not resolvable from this module's imports" -- a method on a
    local object, a dynamic dispatch, a builtin. Unresolved calls are COUNTED
    by the call-graph walker rather than ignored, because a walk that silently
    loses edges makes a negative claim look true.
    """
    func = node.func
    if isinstance(func, ast.Name):
        if func.id in mod.symbol_aliases:
            return mod.symbol_aliases[func.id]
        if func.id in mod.functions:
            return (mod.dotted, func.id)
        return None
    if isinstance(func, ast.Attribute) and isinstance(func.value, ast.Name):
        owner = func.value.id
        if owner in mod.module_aliases:
            return (mod.module_aliases[owner], func.attr)
        return None
    return None


def calls_in(node: ast.AST) -> List[ast.Call]:
    return [n for n in ast.walk(node) if isinstance(n, ast.Call)]


#: Trees swept when a claim asks "does anything call this?".
SOURCE_ROOTS: Tuple[str, ...] = ("tmesim", "app", "experiments", "scripts")


def python_files(roots: Sequence[str]) -> List[Path]:
    out: List[Path] = []
    for root in roots:
        base = REPO / root
        if base.is_dir():
            out.extend(sorted(base.rglob("*.py")))
    return out


def module_for_path(path: Path) -> Optional[Module]:
    try:
        rel = path.relative_to(REPO)
    except ValueError:  # pragma: no cover
        return None
    dotted = ".".join(rel.with_suffix("").parts)
    if dotted.endswith(".__init__"):
        dotted = dotted[: -len(".__init__")]
    try:
        return _parse(dotted, path)
    except (SyntaxError, UnicodeDecodeError):  # pragma: no cover
        return None


# =========================================================================== #
# Descriptor-key reads
# =========================================================================== #
@dataclass(frozen=True)
class Read:
    """One place a string-literal key is read out of a mapping."""

    key: str
    path: str
    lineno: int
    how: str  # "[]" or ".get()"

    def where(self) -> str:
        return f"{self.path}:{self.lineno}"


def key_reads(node: ast.AST, rel: str) -> List[Read]:
    """Every ``x["k"]`` and ``x.get("k", ...)`` under ``node``.

    Deliberately blind to WHICH mapping is being read. Proving that the object
    is the descriptor dict would need type inference, and for the keys anyone
    asks about here -- ``logp``, ``hbd``, ``tpsa``, ``mw`` -- there is no other
    mapping in this package that uses the same word as a key. The claim being
    decided is "does this descriptor reach this law", and a key-name match
    answers it without pretending to more precision than it has.
    """
    out: List[Read] = []
    for n in ast.walk(node):
        if isinstance(n, ast.Subscript):
            sl = n.slice
            # Python 3.9 wraps a simple subscript in ast.Index on some builds.
            sl = getattr(sl, "value", sl) if isinstance(sl, ast.Index) else sl
            if isinstance(sl, ast.Constant) and isinstance(sl.value, str):
                out.append(Read(sl.value, rel, n.lineno, "[]"))
        elif (
            isinstance(n, ast.Call)
            and isinstance(n.func, ast.Attribute)
            and n.func.attr == "get"
            and n.args
            and isinstance(n.args[0], ast.Constant)
            and isinstance(n.args[0].value, str)
        ):
            out.append(Read(n.args[0].value, rel, n.lineno, ".get()"))
    return out


@dataclass
class Walk:
    """The result of a bounded call-graph walk."""

    reads: List[Read]
    visited: Set[Tuple[str, str]]
    unresolved: Set[str]


def walk_call_graph(entry: str, package: str, max_nodes: int = 400) -> Walk:
    """Descriptor-key reads reachable from ``entry``, staying inside ``package``.

    Bounded on purpose. The chain this exists for --
    ``two_route_diffusivity`` and the handful of ``molecular`` helpers under it
    -- is short and lives in one package, and a walk that wandered into numpy
    would be answering a different question badly.

    Unresolved callees are RECORDED, not swallowed: a negative claim proved by
    a walk that quietly lost half its edges is worth nothing, and
    ``test_the_walk_reaches_the_chain_it_claims_to`` pins the floor.
    """
    mod_name, _, func_name = entry.rpartition(".")
    reads: List[Read] = []
    visited: Set[Tuple[str, str]] = set()
    unresolved: Set[str] = set()
    queue: List[Tuple[str, str]] = [(mod_name, func_name)]

    while queue and len(visited) < max_nodes:
        key = queue.pop()
        if key in visited:
            continue
        visited.add(key)
        mod = load_module(key[0])
        if mod is None:
            continue
        node = mod.functions.get(key[1])
        if node is None:
            continue
        reads.extend(key_reads(node, mod.rel))
        for call in calls_in(node):
            target = call_target(call, mod)
            if target is None:
                name = getattr(call.func, "attr", None) or getattr(call.func, "id", "?")
                unresolved.add(str(name))
                continue
            if target[0] == package or target[0].startswith(package + "."):
                queue.append(target)
    return Walk(reads=reads, visited=visited, unresolved=unresolved)


# =========================================================================== #
# The claim types
# =========================================================================== #
@dataclass(frozen=True)
class Claim:
    """One statement a document makes about this codebase, made decidable.

    ``anchor`` is a short substring that IS the claim as the document phrases
    it. It is checked to still be present, so a declaration cannot outlive its
    sentence -- delete the sentence and the declaration fails as stale, which
    is the only way an opt-in table stays honest.

    ``reason`` is mandatory and validated at construction, exactly as
    ``propagation.Declaration`` and ``propagation.WiringPoint`` do. It is not
    decoration: the known failure mode of this file is somebody re-stamping a
    red assertion without reading it, and a justification that has to be
    rewritten alongside the value is the only friction that reliably stops it.
    """

    document: str
    anchor: str
    reason: str

    def __post_init__(self) -> None:
        if len(self.reason.strip()) < 25:
            raise ValueError(
                f"{self.document}: a Claim needs a written reason (>=25 chars) "
                "saying what goes false in this document if the code moves. "
                "Silencing an assertion without writing that down is the "
                "failure mode this table exists to make expensive."
            )
        if len(self.anchor.strip()) < 8:
            raise ValueError(
                f"{self.document}: an anchor must be long enough to be the "
                "claim, not merely to occur in the file."
            )

    # -- subclass hooks ----------------------------------------------------- #
    def subject(self) -> str:
        """The claim, in one clause, for the failure message."""
        raise NotImplementedError

    def evidence(self) -> Optional[str]:
        """``None`` if the claim holds; otherwise what the code actually does."""
        raise NotImplementedError

    # -- the check ---------------------------------------------------------- #
    def anchor_present(self) -> bool:
        path = REPO / self.document
        if not path.is_file():
            return False
        return self.anchor in path.read_text(encoding="utf-8")

    def check(self) -> Optional[str]:
        got = self.evidence()
        if got is None:
            return None
        return f"{self.document} asserts {self.subject()}; {got}."


@dataclass(frozen=True)
class CallEdge(Claim):
    """"Module A does (not) reach module B at all."

    The shape of *"Nothing is wired into ``compile_drug.py``"*.
    """

    caller: str = ""   #: dotted module making the calls
    callee: str = ""   #: dotted module being called
    present: bool = True

    def subject(self) -> str:
        verb = "reaches" if self.present else "does not reach"
        return f"`{self.caller}` {verb} `{self.callee}`"

    def _edges(self) -> List[Tuple[str, int]]:
        mod = load_module(self.caller)
        if mod is None:
            return []
        found = []
        for call in calls_in(mod.tree):
            target = call_target(call, mod)
            if target and target[0] == self.callee:
                found.append((f"{mod.rel}:{call.lineno}", call.lineno))
        return found

    def evidence(self) -> Optional[str]:
        edges = self._edges()
        if self.present and not edges:
            return (
                f"`{self.caller}` no longer calls anything in `{self.callee}` "
                "-- the wiring this document describes is gone"
            )
        if not self.present and edges:
            shown = ", ".join(w for w, _ in edges[:4])
            return f"it is called at {shown}"
        return None


@dataclass(frozen=True)
class Const(Claim):
    """"This function still returns this literal."

    The highest-yield verb in the file. ``calibrate_alpha()`` is quoted to
    seventeen figures in at least six documents, and every one of them went
    false in the same instant when ``fd2a9e7`` moved it.

    Compared with ``==`` and no tolerance, because that is what the documents
    say: they quote all seventeen digits and ``tests/test_two_route.py`` already
    pins the value bit-exactly. A tolerance here would be this file quietly
    disagreeing with the claim it is guarding.
    """

    dotted: str = ""
    expected: Any = None
    call: bool = True

    def subject(self) -> str:
        suffix = "()" if self.call else ""
        return f"`{self.dotted}{suffix} == {self.expected!r}`"

    def evidence(self) -> Optional[str]:
        mod_name, _, attr = self.dotted.rpartition(".")
        obj = getattr(importlib.import_module(mod_name), attr)
        got = _memo_call(self.dotted, obj) if self.call else obj
        if got == self.expected:
            return None
        return (
            f"it returns {got!r} (a difference of "
            f"{got - self.expected!r})"
            if isinstance(got, (int, float)) and isinstance(self.expected, (int, float))
            else f"it is {got!r}"
        )


_CALL_MEMO: Dict[str, Any] = {}


def _memo_call(key: str, fn: Any) -> Any:
    """One evaluation per constant per session. Six documents quote ``alpha``;
    solving it six times would be six identical closed-form solves."""
    if key not in _CALL_MEMO:
        _CALL_MEMO[key] = fn()
    return _CALL_MEMO[key]


@dataclass(frozen=True)
class Default(Claim):
    """"This switch still ships in this position."

    The shape of *"default-OFF and wired into nothing"*. Checked three ways,
    because they fail separately:

    * the LITERAL in the source, so an edit to the file is seen even if the
      module cannot be imported in this environment;
    * the RUNTIME value, so a constant that is computed rather than written is
      seen too;
    * the SIGNATURE default of the function that is supposed to read it, which
      is the one a source-literal check misses entirely -- a constant can stay
      ``True`` while ``build_profile`` stops defaulting to it.
    """

    dotted: str = ""
    expected: Any = None
    #: ``(dotted function, parameter)`` that must default to this constant.
    bound_to: Optional[Tuple[str, str]] = None

    def subject(self) -> str:
        return f"`{self.dotted} is {self.expected!r}`"

    def _literal(self) -> Tuple[bool, Any]:
        mod_name, _, attr = self.dotted.rpartition(".")
        mod = load_module(mod_name)
        if mod is None:
            return False, None
        for node in mod.tree.body:
            targets = (
                node.targets
                if isinstance(node, ast.Assign)
                else ([node.target] if isinstance(node, ast.AnnAssign) else [])
            )
            for t in targets:
                if isinstance(t, ast.Name) and t.id == attr and node.value is not None:
                    try:
                        return True, ast.literal_eval(node.value)
                    except ValueError:
                        return False, None
        return False, None

    def evidence(self) -> Optional[str]:
        mod_name, _, attr = self.dotted.rpartition(".")
        found, literal = self._literal()
        if not found:
            return (
                f"there is no module-level literal named `{attr}` in "
                f"`{mod_name}` any more"
            )
        if literal != self.expected:
            return f"the source says `{attr} = {literal!r}`"

        runtime = getattr(importlib.import_module(mod_name), attr)
        if runtime != self.expected:
            return f"at runtime it is {runtime!r}"

        if self.bound_to is not None:
            fn_dotted, param = self.bound_to
            fn_mod, _, fn_name = fn_dotted.rpartition(".")
            fn = getattr(importlib.import_module(fn_mod), fn_name)
            params = inspect.signature(fn).parameters
            if param not in params:
                return f"`{fn_dotted}` no longer takes a `{param}` parameter"
            if params[param].default != self.expected:
                return (
                    f"`{fn_dotted}({param}=...)` now defaults to "
                    f"{params[param].default!r}, so the constant is no longer "
                    "what ships"
                )
        return None


@dataclass(frozen=True)
class NoCallers(Claim):
    """"This function is still inert / is now wired in."

    This is the verb that would have fired on ``resolve_logp`` BEFORE somebody
    wired it up trusting a document that said it was transport-inert. The scope
    is the SHIPPED trees; a call from ``experiments/`` is a probe, not wiring,
    and treating one as the other would make this fire on ordinary analysis
    work -- which is how a guard gets routed around.
    """

    dotted: str = ""
    scope: Tuple[str, ...] = ("tmesim", "app")
    present: bool = False  #: False asserts there are NO callers

    def subject(self) -> str:
        state = "is called" if self.present else "has no caller"
        return f"`{self.dotted}` {state} in {'/'.join(self.scope)}"

    def _callers(self) -> List[str]:
        home, _, symbol = self.dotted.rpartition(".")
        out: List[str] = []
        for path in python_files(self.scope):
            mod = module_for_path(path)
            if mod is None or mod.dotted == home:
                continue  # the defining module's own self-test is not wiring
            for call in calls_in(mod.tree):
                target = call_target(call, mod)
                if target == (home, symbol):
                    out.append(f"{mod.rel}:{call.lineno}")
        return sorted(out)

    def evidence(self) -> Optional[str]:
        callers = self._callers()
        if self.present and not callers:
            return f"nothing in {'/'.join(self.scope)} calls it any more"
        if not self.present and callers:
            return (
                "it is called at " + ", ".join(callers[:4])
                + " -- this document's advice about it is now advice about "
                "live code"
            )
        return None


@dataclass(frozen=True)
class Descriptor(Claim):
    """"This molecular descriptor key does (not) reach this physics law."

    The hard one, and the one that catches the defect the corrections were
    about: *"logP is not an input to transport at any point."*

    HOW THE QUESTION IS POSED, because the obvious posing is wrong. Walking the
    call graph from ``build_profile`` and asking whether ``logp`` is read
    anywhere on it answers YES both before and after ``fd2a9e7`` -- ``build_
    profile`` has always called ``permeability.membrane_permeability``, which
    has always read ``descriptors["logp"]``. That reading is TRUE and it is not
    what the document meant; a guard built on it would have been red against a
    correct document, which is the worst thing a guard can be.

    What the document meant is narrower and is exactly decidable: does the key
    reach the TRANSPORT LAW? So the entry point is the law itself, and the keys
    that reach it are

      (a) those read inside the ARGUMENT EXPRESSIONS of each call to it, which
          is where ``fd2a9e7`` put ``logp`` and ``hbd``; plus
      (b) those read anywhere in its own call graph, so routing a descriptor in
          through a helper does not launder it.

    Before ``fd2a9e7`` there were no call sites at all, so the set was empty and
    the negative claim was true -- which is the correct answer, and the reason
    this posing is the right one.
    """

    entry: str = ""                     #: dotted function -- the law
    key: str = ""                       #: descriptor key
    callers: Tuple[str, ...] = ()       #: dotted modules that call the law
    package: str = "tmesim.molecular"   #: walk boundary
    present: bool = True

    def subject(self) -> str:
        verb = "is read on the path to" if self.present else "is not read on the path to"
        law = self.entry.rpartition(".")[2]
        return f"`{self.key}` {verb} `{law}`"

    def reads(self) -> List[Read]:
        """(a) argument expressions at each call site, plus (b) the law's own walk."""
        found: List[Read] = []
        symbol = self.entry.rpartition(".")[2]
        home = self.entry.rpartition(".")[0]
        for caller in self.callers:
            mod = load_module(caller)
            if mod is None:
                continue
            for call in calls_in(mod.tree):
                if call_target(call, mod) != (home, symbol):
                    continue
                for arg in list(call.args) + [kw.value for kw in call.keywords]:
                    found.extend(key_reads(arg, mod.rel))
        found.extend(walk_call_graph(self.entry, self.package).reads)
        return found

    def evidence(self) -> Optional[str]:
        hits = [r for r in self.reads() if r.key == self.key]
        if self.present and not hits:
            return (
                f"nothing on that path reads `{self.key}` any more -- the "
                "input this document describes is gone"
            )
        if not self.present and hits:
            first = hits[0]
            more = f" (and {len(hits) - 1} more)" if len(hits) > 1 else ""
            return f"{first.where()} reads it{more}"
        return None


# =========================================================================== #
# 1. The declared claims
#
# Each entry: the document, the sentence in it, and the reason that sentence is
# load-bearing. A red assertion here means one of exactly two things -- the code
# moved and a document is now lying, or the document was rewritten and this
# entry is stale. Both are decisions. Neither is a re-stamp.
# =========================================================================== #
ALPHA = 22.223477711172563
ALPHA_SUPERSEDED = 148.60815765445545

_ALPHA_REASON = (
    "This document quotes `alpha` as the CURRENT value, not as a historical "
    "run record. It was 148.60815765445545 until fd2a9e7 re-solved it against "
    "the same Primeau anchor on a two-route diffusivity, and the stale literal "
    "then stood in at least four documents at once. calibrate_alpha() is a "
    "closed-form solve of one scalar equation, so it moves whenever any input "
    "to a compiled doxorubicin moves -- which is to say, whenever the "
    "molecular pipeline is touched at all. Every absolute exposure number in "
    "this repo rides on it."
)


CLAIMS: Tuple[Claim, ...] = (
    # ---------------- results/two_route_transport.md ----------------------- #
    CallEdge(
        document="results/two_route_transport.md",
        anchor="`compile_drug.build_profile` now calls",
        caller="tmesim.molecular.compile_drug",
        callee="tmesim.molecular.two_route",
        present=True,
        reason=(
            "This file's banner exists because the header underneath it says "
            "'Nothing is wired into compile_drug.py'. That sentence was true "
            "on 2026-08-13 and false on 2026-08-14, and nothing noticed for "
            "months. The banner asserts the opposite edge; if the edge is ever "
            "removed the banner becomes the wrong correction and the header "
            "becomes right again, which is a state no reader could untangle."
        ),
    ),
    Const(
        document="results/two_route_transport.md",
        anchor="**22.223477711172563**",
        dotted="tmesim.validation.calibration.calibrate_alpha",
        expected=ALPHA,
        reason=_ALPHA_REASON + " Here it is the banner's own correction of the "
        "superseded 148.60815765445545 quoted directly beneath it, so a stale "
        "value would leave two contradictory numbers on one page.",
    ),
    # ---------------- results/logp_provenance.md --------------------------- #
    Descriptor(
        document="results/logp_provenance.md",
        anchor="`logp` is an input to `effective_diffusivity`",
        entry="tmesim.molecular.two_route.two_route_diffusivity",
        key="logp",
        callers=("tmesim.molecular.compile_drug",),
        present=True,
        reason=(
            "The exact claim that went false. This file's finding 2 said 'logP "
            "is not an input to transport at any point' and enumerated logP's "
            "only three routes to an output; fd2a9e7 added a fourth, the "
            "transport coefficient itself. The banner asserts the corrected "
            "direction. If the cell term ever stops consuming logp the banner "
            "is wrong and the original finding is right again."
        ),
    ),
    NoCallers(
        document="results/logp_provenance.md",
        anchor="still has no caller",
        dotted="tmesim.molecular.lipophilicity.resolve_logp",
        scope=("tmesim", "app"),
        present=False,
        reason=(
            "The assertion this whole file was written to protect. "
            "resolve_logp holds measured logP values that supersede Crippen "
            "for eleven drugs; the banner warns that wiring it up is now a "
            "change to the TRANSPORT COEFFICIENT of every compound it "
            "supersedes and needs a revalidation run, not a provenance review. "
            "The moment someone wires it in, that warning has to be acted on "
            "-- and the person doing the wiring is exactly the person who will "
            "have read this document's older finding 2 and concluded it was "
            "inert."
        ),
    ),
    Const(
        document="results/logp_provenance.md",
        anchor="**22.223477711172563**",
        dotted="tmesim.validation.calibration.calibrate_alpha",
        expected=ALPHA,
        reason=_ALPHA_REASON + " Here it is load-bearing twice over: the "
        "banner uses it to explain that section 5.3's 148.6082 -> 68.6368 "
        "measurement no longer even starts from the right number.",
    ),
    # ---------------- results/descriptor_falsification.md ------------------ #
    Default(
        document="results/descriptor_falsification.md",
        anchor="`_TWO_ROUTE_DEFAULT = True`",
        dotted="tmesim.molecular.compile_drug._TWO_ROUTE_DEFAULT",
        expected=True,
        bound_to=("tmesim.molecular.compile_drug.build_profile", "two_route"),
        reason=(
            "This file's entire subject is that the SHIPPING descriptor layer "
            "has never been scored, and that argument only bites because the "
            "two-route term is applied on every compile. If the default flips "
            "back to False the file is describing a hypothetical rather than "
            "the shipping engine, and its central complaint quietly stops "
            "being about anything a user runs."
        ),
    ),
    NoCallers(
        document="results/descriptor_falsification.md",
        anchor="has no callers",
        dotted="tmesim.molecular.lipophilicity.resolve_logp",
        scope=("tmesim", "app"),
        present=False,
        reason=(
            "This file's section 2 turns on resolve_logp being available and "
            "unapplied -- it argues the measured-vs-Crippen question is "
            "ORTHOGONAL to the descriptor-convention question it is actually "
            "about. If resolve_logp gets wired in, that orthogonality argument "
            "is about a configuration that no longer exists, and the file's "
            "recommendation ordering goes with it."
        ),
    ),
    # ---------------- docs/QUICKSTART.md ----------------------------------- #
    Default(
        document="docs/QUICKSTART.md",
        anchor="comes from a two-route transport model",
        dotted="tmesim.molecular.compile_drug._TWO_ROUTE_DEFAULT",
        expected=True,
        bound_to=("tmesim.molecular.compile_drug.build_profile", "two_route"),
        reason=(
            "The first thing a new user reads about D_eff. It said the law was "
            "default-OFF and wired into nothing long after it had become the "
            "shipping default, so the document that introduces the product "
            "described a different product. Every caveat in this callout -- "
            "R^2 = 0.718, 24 of 31 pairs, 'read a ranking as a hypothesis' -- "
            "is a caveat about a model the reader is only running if this "
            "default is True."
        ),
    ),
    Descriptor(
        document="docs/QUICKSTART.md",
        anchor="hydrogen-bond donor count now reach transport",
        entry="tmesim.molecular.two_route.two_route_diffusivity",
        key="hbd",
        callers=("tmesim.molecular.compile_drug",),
        present=True,
        reason=(
            "The sentence that tells a user WHY two similar-mass compounds can "
            "now rank differently. hbd is checked rather than logp because it "
            "is the weaker of the two claims: the donor count reaches "
            "transport only through the membrane-permeability leg, and this "
            "same callout goes on to disclose that the shipping donor count "
            "uses the other of RDKit's two conventions. If hbd stops reaching "
            "transport, that disclosure is describing a risk the engine no "
            "longer takes and the ranking caveat above it is over-strict."
        ),
    ),
    # ---------------- results/CLAIM_INVENTORY.md --------------------------- #
    Const(
        document="results/CLAIM_INVENTORY.md",
        anchor="**`22.223477711172563`**",
        dotted="tmesim.validation.calibration.calibrate_alpha",
        expected=ALPHA,
        reason=_ALPHA_REASON + " This document is the INDEX of the claim "
        "surface -- the file that says results/*.md is what this product "
        "claims. A stale number here is stale at the top of the tree of "
        "trust, and every file it points at inherits it.",
    ),
    # ---------------- results/area_conservation_promotion.md --------------- #
    Const(
        document="results/area_conservation_promotion.md",
        anchor="`calibrate_alpha()` returns **22.223477711172563** before and after the flip",
        dotted="tmesim.validation.calibration.calibrate_alpha",
        expected=ALPHA,
        reason=_ALPHA_REASON + " This file's promotion argument rests on alpha "
        "NOT moving across the area-conservation flip -- it is quoted as an "
        "invariant that was checked, twice, on either side of the change. An "
        "invariant quoted at a value that is no longer the value is worse than "
        "no invariant.",
    ),
    # ---------------- results/rebrand_backend.md --------------------------- #
    Const(
        document="results/rebrand_backend.md",
        anchor="returns **22.223477711172563** — exact",
        dotted="tmesim.validation.calibration.calibrate_alpha",
        expected=ALPHA,
        reason=_ALPHA_REASON + " Here it is the evidence that a pure rename "
        "changed no physics: the whole point of quoting it is that it is "
        "bit-identical across the rebrand. If it moves for an unrelated "
        "reason, this file is asserting a no-op that was not one.",
    ),
    # ---------------- results/transport_anchor_prereg.md ------------------- #
    Const(
        document="results/transport_anchor_prereg.md",
        anchor="| **22.223477711172563** |",
        dotted="tmesim.validation.calibration.calibrate_alpha",
        expected=ALPHA,
        reason=_ALPHA_REASON + " This is a PRE-REGISTRATION: it fixes the "
        "values a future measurement will be scored against, before the "
        "measurement is made. A pre-registered constant that silently drifts "
        "away from the engine destroys the only thing pre-registration buys.",
    ),
)


# =========================================================================== #
# 2. Proof of correctness -- the checker decides what it claims to decide
#
# Every one of these runs the real checker against the real tree, in BOTH
# directions, so a verb that had quietly stopped being able to fail would be
# caught here rather than by being green forever.
# =========================================================================== #
def test_the_verb_that_would_have_caught_the_headline_defect():
    """The exact sentence from the brief, produced by the real checker.

    ``results/logp_provenance.md`` said *"logP is not an input to transport at
    any point"*. Posed as a claim, that is ``present=False`` -- and against
    today's tree it must fail, naming the document and the line.
    """
    stale = Descriptor(
        document="results/logp_provenance.md",
        anchor="`logp` is an input to `effective_diffusivity`",
        entry="tmesim.molecular.two_route.two_route_diffusivity",
        key="logp",
        callers=("tmesim.molecular.compile_drug",),
        present=False,
        reason=(
            "The pre-fd2a9e7 claim, kept here as the proof that this verb can "
            "fail. It is NOT in CLAIMS -- the document no longer says it."
        ),
    )
    msg = stale.check()
    assert msg is not None, "the verb has gone blind to the defect it exists for"
    assert msg.startswith(
        "results/logp_provenance.md asserts `logp` is not read on the path to "
        "`two_route_diffusivity`;"
    ), msg
    assert "tmesim/molecular/compile_drug.py:" in msg, msg
    assert re.search(r"compile_drug\.py:\d+ reads it", msg), (
        "the verb must cite the LINE it found the read on, not merely the "
        f"file: {msg}")


def test_the_descriptor_verb_is_not_fooled_by_the_permeability_leg():
    """The posing that would have made this verb useless, refused.

    ``build_profile`` has ALWAYS read ``descriptors["logp"]`` -- it hands the
    whole dict to ``permeability.membrane_permeability``. A verb rooted at
    ``build_profile`` would therefore have been red against the ORIGINAL,
    CORRECT document. Rooting it at the transport law is what makes it decide
    the question the document actually asked.
    """
    from_law = walk_call_graph(
        "tmesim.molecular.two_route.two_route_diffusivity", "tmesim.molecular")
    from_compile = walk_call_graph(
        "tmesim.molecular.compile_drug.build_profile", "tmesim.molecular")

    assert "logp" in {r.key for r in from_compile.reads}, (
        "build_profile does read logp -- via the permeability leg -- which is "
        "precisely why it is the wrong root for this question")
    assert "logp" not in {r.key for r in from_law.reads}, (
        "the law takes logp as a plain float argument; it reads no descriptors")
    # And the claim as declared still finds it, because (a) covers the call site.
    declared = [c for c in CLAIMS
                if isinstance(c, Descriptor) and c.key == "logp"][0]
    assert "logp" in {r.key for r in declared.reads()}


def test_the_descriptor_verb_sees_a_key_routed_through_a_helper():
    """Moving the read one call away must not launder it."""
    walk = walk_call_graph(
        "tmesim.molecular.compile_drug.build_profile", "tmesim.molecular")
    # membrane_permeability lives in a different module and is reached only
    # through the call graph, so its keys prove the walk crosses modules.
    assert ("tmesim.molecular.permeability", "membrane_permeability") in walk.visited
    assert any(r.path.endswith("permeability.py") and r.key == "logp"
               for r in walk.reads)


def test_the_walk_reaches_the_chain_it_claims_to():
    """A walk that silently loses edges makes a negative claim look true.

    Pin the floor: the real chain, and a bound on how much went unresolved.
    """
    walk = walk_call_graph(
        "tmesim.molecular.compile_drug.build_profile", "tmesim.molecular")
    for expected in (
        ("tmesim.molecular.compile_drug", "build_profile"),
        ("tmesim.molecular.two_route", "two_route_diffusivity"),
        ("tmesim.molecular.permeability", "membrane_permeability"),
        ("tmesim.molecular.binding", "plasma_free_fraction"),
        ("tmesim.molecular.stokes_einstein", "free_diffusivity"),
        ("tmesim.molecular.ecm_mesh", "ogston_renkin_hindrance"),
        ("tmesim.molecular.brownian", "simulate_effective_diffusivity"),
    ):
        assert expected in walk.visited, f"the walk no longer reaches {expected}"
    assert len(walk.reads) > 10, "the walk found almost no descriptor reads"


def test_a_mention_in_prose_is_not_a_call():
    """These modules discuss their own call graph at length in docstrings.

    ``ecm_mesh.py``'s module docstring contains the words
    ``compile_drug.build_profile calls two_route.two_route_diffusivity``. A
    grep-based guard would read that as an edge. The AST does not.
    """
    src = (REPO / "tmesim/molecular/ecm_mesh.py").read_text(encoding="utf-8")
    assert "compile_drug.build_profile` calls" in src
    edge = CallEdge(
        document="results/CLAIM_INVENTORY.md",
        anchor="claim surface",
        caller="tmesim.molecular.ecm_mesh",
        callee="tmesim.molecular.compile_drug",
        present=True,
        reason="proof fixture: ecm_mesh only MENTIONS compile_drug, in prose",
    )
    assert edge.evidence() is not None, (
        "a prose mention was counted as a call edge")


def test_the_call_edge_verb_sees_every_spelling_of_the_same_edge():
    """The same edge written four ways, and one lookalike that is not it.

    ``_tr.f()``, ``two_route.f()``, a bare ``f()`` from a ``from ... import``,
    and an absolute ``from tmesim.molecular.two_route import f`` all resolve to
    one edge. A local variable that merely happens to be called ``two_route``
    resolves to nothing -- which is what stops this from being grep.

    Parsed from source rather than written to disk: a probe file dropped into
    ``tmesim/molecular/`` is visible to ``propagation.sweep_repo`` and would
    make a concurrent test run fail for an unrelated reason.
    """
    import textwrap

    mod = parse_source("tmesim.molecular._probe", textwrap.dedent('''
        from . import two_route as _tr
        from . import two_route
        from .two_route import two_route_diffusivity
        from tmesim.molecular.two_route import two_route_diffusivity as tr4

        def a(): return _tr.two_route_diffusivity(1)
        def b(): return two_route.two_route_diffusivity(2)
        def c(): return two_route_diffusivity(3)
        def d(): return tr4(4)
        def e(two_route): return two_route.two_route_diffusivity(5)
    '''), REPO / "tmesim/molecular/_probe.py")

    targets = [call_target(c, mod) for c in calls_in(mod.tree)]
    edge = ("tmesim.molecular.two_route", "two_route_diffusivity")
    # a, b, c, d are the four real spellings. e is the LIMITATION, asserted
    # here so it is recorded rather than discovered: its parameter shadows the
    # module name, the resolver keys off the module-level import table, and so
    # it reports an edge that is not one. That direction is the safe one for a
    # NEGATIVE claim (it can only over-report a caller, never miss one), and
    # nothing in tmesim/ shadows a module name this way.
    assert targets.count(edge) == 5, targets
    assert len(targets) == 5


def test_the_no_callers_verb_can_fail_in_both_directions():
    """It is inert today; asserting the opposite must go red, and a function
    that IS wired in must satisfy the positive form."""
    inert = NoCallers(
        document="results/logp_provenance.md",
        anchor="still has no caller",
        dotted="tmesim.molecular.lipophilicity.resolve_logp",
        present=True,  # deliberately wrong
        reason="proof fixture: resolve_logp is inert, so present=True must fail",
    )
    assert inert.evidence() is not None
    wired = NoCallers(
        document="results/logp_provenance.md",
        anchor="still has no caller",
        dotted="tmesim.molecular.compile_drug.molecule_to_drugparams",
        present=True,
        reason="proof fixture: molecule_to_drugparams is genuinely wired in",
    )
    assert wired.evidence() is None


def test_the_defining_modules_own_self_test_is_not_a_caller():
    """``lipophilicity.py``'s ``__main__`` block calls ``resolve_logp`` five
    times. Counting those would make the claim permanently, uselessly red."""
    src = (REPO / "tmesim/molecular/lipophilicity.py").read_text(encoding="utf-8")
    assert src.count("resolve_logp(") >= 4, "the self-test calls vanished"
    claim = [c for c in CLAIMS if isinstance(c, NoCallers)][0]
    assert claim._callers() == []


def test_the_const_verb_reports_the_superseded_value_when_it_is_wrong():
    stale = Const(
        document="results/two_route_transport.md",
        anchor="**22.223477711172563**",
        dotted="tmesim.validation.calibration.calibrate_alpha",
        expected=ALPHA_SUPERSEDED,   # the pre-fd2a9e7 literal
        reason="proof fixture: the literal that stood in four documents at once",
    )
    msg = stale.check()
    assert msg is not None
    assert repr(ALPHA) in msg
    assert msg.startswith("results/two_route_transport.md asserts")


def test_the_default_verb_checks_the_signature_and_not_only_the_constant():
    """A constant can stay ``True`` while the function stops defaulting to it.
    That is the failure a source-literal check cannot see."""
    from tmesim.molecular import compile_drug

    good = [c for c in CLAIMS if isinstance(c, Default)][0]
    assert good.evidence() is None
    assert (inspect.signature(compile_drug.build_profile)
            .parameters["two_route"].default is compile_drug._TWO_ROUTE_DEFAULT)

    wrong_param = Default(
        document="docs/QUICKSTART.md",
        anchor="comes from a two-route transport model",
        dotted="tmesim.molecular.compile_drug._TWO_ROUTE_DEFAULT",
        expected=True,
        bound_to=("tmesim.molecular.compile_drug.build_profile", "no_such_param"),
        reason="proof fixture: the signature half of the check must be live",
    )
    assert "no longer takes" in (wrong_param.evidence() or "")


def test_a_claim_cannot_be_declared_without_a_reason():
    """The friction IS the mechanism. Construction refuses a bare silencing."""
    with pytest.raises(ValueError):
        Const(document="d.md", anchor="a long enough anchor", reason="",
              dotted="x.y", expected=1)
    with pytest.raises(ValueError):
        Const(document="d.md", anchor="a long enough anchor", reason="updated",
              dotted="x.y", expected=1)
    with pytest.raises(ValueError):
        Const(document="d.md", anchor="short",
              reason="a perfectly good and sufficiently long written reason",
              dotted="x.y", expected=1)


# =========================================================================== #
# 3. The claims themselves
# =========================================================================== #
@pytest.mark.parametrize(
    "claim", CLAIMS,
    ids=[f"{c.document.rsplit('/', 1)[-1]}:{type(c).__name__}:{i}"
         for i, c in enumerate(CLAIMS)])
def test_every_declared_claim_still_holds(claim: Claim):
    failure = claim.check()
    assert failure is None, (
        failure + "\n\n"
        "A document on this product's claim surface now says something the "
        "code does not do. Two honest resolutions, and re-stamping this "
        "assertion is neither:\n"
        "  (a) the code moved -- correct the document (this repo's convention "
        "is a dated banner that leaves the superseded sentence standing), then "
        "update the value here;\n"
        "  (b) the code should not have moved -- fix the code.\n\n"
        f"Why this claim is load-bearing: {claim.reason}"
    )


@pytest.mark.parametrize(
    "claim", CLAIMS,
    ids=[f"{c.document.rsplit('/', 1)[-1]}:{i}" for i, c in enumerate(CLAIMS)])
def test_no_claim_is_stale(claim: Claim):
    """The document must still make the claim, in the words declared here.

    An opt-in table has exactly one way to rot: a document is rewritten, the
    sentence goes, and an entry survives asserting something nobody claims any
    more. The anchor closes it -- and a reword that keeps the claim keeps the
    anchor, because the anchor is chosen to BE the claim rather than to be
    near it.
    """
    assert (REPO / claim.document).is_file(), f"{claim.document} is gone"
    assert claim.anchor_present(), (
        f"{claim.document} no longer contains {claim.anchor!r}.\n"
        "Either the claim was withdrawn -- delete this entry -- or it was "
        "reworded, in which case re-read it and move the anchor. Do not widen "
        "the anchor until it matches; an anchor that is not the claim guards "
        "nothing.\n"
        f"Declared reason: {claim.reason}"
    )


def test_every_claim_carries_a_re_checkable_reason():
    """25 characters is what construction enforces; this is the standard.

    A reason has to survive being read by someone who did not write it, six
    months later, at the moment the assertion is red and they are in a hurry.
    """
    for claim in CLAIMS:
        assert len(claim.reason.strip()) >= 60, (
            f"{claim.document}: {type(claim).__name__} is declared with a "
            "one-liner. Write down what goes false in this document if the "
            "code moves.")


def test_the_table_covers_every_verb_and_more_than_one_document():
    """A table that has quietly collapsed onto one verb or one file is not the
    guard this was built to be. Pin the floor rather than the exact contents,
    so ordinary additions do not fire this."""
    kinds = {type(c).__name__ for c in CLAIMS}
    assert kinds == {"CallEdge", "Const", "Default", "NoCallers", "Descriptor"}
    documents = {c.document for c in CLAIMS}
    assert len(documents) >= 6, documents
    # The alpha literal is the recurrence this file was built for: it stood
    # stale in at least four documents at once. Guarding it in one would be
    # guarding the wrong thing.
    alpha_docs = {c.document for c in CLAIMS
                  if isinstance(c, Const) and c.expected == ALPHA}
    assert len(alpha_docs) >= 4, alpha_docs


def test_no_claim_asserts_the_superseded_alpha():
    """A guard that pins the WRONG number is worse than no guard: it would hold
    the repo to the value the corrections were made to remove."""
    for claim in CLAIMS:
        if isinstance(claim, Const):
            assert claim.expected != ALPHA_SUPERSEDED, claim.document
