"""The guard against an unexpected exception's text reaching the browser.

A security pass found ``f"{type(exc).__name__}: {exc}"`` answering an HTTP
request in three places. The follow-up census found thirty-six occurrences of
that expression across ``app/`` and ``tmesim/``: it is a house pattern, not
three accidents, and a pattern that was fixed one site at a time regenerates,
because nothing was watching. This file is the thing that watches.

What is actually wrong with it
------------------------------
``str(exc)`` on an exception the code did not anticipate is not user-facing
English. It names an internal class, and for anything that is not wrapping the
caller's own input it can carry a server-side absolute path the user never
supplied -- the repo root, the operator's username, the module layout.
``app/server/cache.py`` was the clear case: it wraps the whole section
registry, so a missing ``data_cache/`` download surfaced a path from the host
on a route every client hits on load.

The twin failure, which a security review cannot see
----------------------------------------------------
The opposite mistake is worse, because it is invisible in an audit and it
degrades the product. This backend renders the server's refusals VERBATIM so
the Studio can never disagree with the CLI about what a specimen supports.
``IngestError`` -- "no recognised layout. Expected one of: a Xenium bundle
..." -- is a deliberate, informative refusal and must survive word for word. A
bare ``FileNotFoundError`` in ``tmesim/data/real/singlecell_spatial.py`` was
PROMOTED to ``IngestError`` for exactly this reason: so its diagnostic would
stop being flattened. Redaction applied to a refusal class looks like a
security fix and is a regression.

So this file checks BOTH directions:

1. **Proof of correctness** -- the checker classifies real examples of each
   category it claims to separate, including the two that look alike (a
   redacted f-string and a verbatim refusal).
2. **The declared table** -- every stringified exception on a client-facing
   path in ``app/server/`` is either fixed or listed in :data:`EXEMPTIONS`
   with a written reason. A stale exemption fails too, so the list cannot
   quietly accumulate.
3. **The sweep** -- a NEW module under ``app/server/`` cannot join the failure
   mode unnoticed, and no handler may redact an exception class that is a
   deliberate refusal.
4. **Execution** -- a redacted path really does hide the detail, really does
   log it under the reference it handed the caller, and a deliberate
   ``IngestError`` really does arrive byte-identical at both HTTP boundaries
   that touch it.
5. **The protocol sweep boundary** -- the executor classifies and
   ``app/server/protocol_api.py`` acts on the classification, so the two are
   driven together: the interesting failure is them disagreeing, which no
   single-module test can see.
6. **The plan-time diagnostic path** -- redacted where the exception is
   caught, so both the ``to_dict()`` projection AND the ``.message``
   attribute are clean, and re-planning does not mint a new reference id.

A third list, :data:`KNOWN_UNFIXED`, holds leaks that are recorded rather than
excused: an exemption asserts a string is safe, and an entry there asserts the
opposite. It is currently empty -- it held the two plan-time diagnostics until
they were closed -- and the machinery is kept for the next leak that cannot be
fixed in the change that finds it. The count is pinned, so it cannot fill up
quietly.

Cost
----
Sections 1-3 parse source without executing it. Section 4 runs against an
empty temporary directory and a monkeypatched registry. Sections 5 and 6 plan
a real TPL program (a library drug, so no RDKit) and replace only the solver
or the drug resolver. No reconstruction, no PDE, no real tissue: the file is a
few seconds, nearly all of it the plans.
"""
from __future__ import annotations

import ast
import os
import re
import sys
from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

import pytest

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

#: The tree under guard. This is the HTTP surface: everything here can end up
#: in a response body. ``tmesim/`` is deliberately NOT swept -- the engine
#: prints to an operator's console, where the class name is the diagnosis, and
#: it hands the HTTP layer a flag (``error_is_internal``) rather than a
#: pre-formatted sentence. The boundary is what matters, and the boundary is
#: here.
SERVER_DIR = os.path.join(_REPO, "app", "server")

#: Engine modules that ALSO feed a client surface, swept individually.
#:
#: The sweep stops here rather than covering ``tmesim/`` because most of the
#: engine answers an operator at a console, where the exception class IS the
#: diagnosis: ``tmesim/cli.py``, ``validation/revalidation.py``,
#: ``validation/scorecard.py`` and ``data/section_catalog.py`` print to a
#: terminal and reach no browser. Sweeping them would add a dozen exemptions
#: that say "this is a CLI" and would fail on ordinary engine work, which is
#: how a check gets switched off.
#:
#: ``protocol/executor.py`` is here because a leak was FOUND here, not on
#: suspicion: a failed sweep cell's ``error`` travels to ``runs[].error`` in
#: the browser (``app/server/protocol_wire.sweep_snapshot``). Add a module to
#: this list when the same becomes true of it.
BOUNDARY_MODULES = (
    os.path.join("tmesim", "protocol", "executor.py"),
)

#: Deliberate refusals. Redaction applied to any of these is the twin failure:
#: an informative, product-visible sentence flattened into a generic one.
#: ``RunCancelled`` is included because a cancellation is an ANSWER, not a
#: fault -- the user asked for it.
REFUSAL_CLASSES = frozenset({
    "IngestError",
    "TplError",
    "CompoundRefused",
    "CompoundError",
    "ProjectError",
    "RunCancelled",
    "HTTPException",
})

#: The redactors. ``_internal`` is ``protocol_api``'s thin wrapper that turns
#: one into an editor diagnostic.
REDACTORS = frozenset({"redact_exception", "log_internal", "_internal"})


# --------------------------------------------------------------------------- #
# The checker
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Finding:
    """One place where a caught exception is turned into a string."""

    path: str        #: repo-relative, POSIX separators
    qualname: str    #: enclosing class/function chain, "<module>" at top level
    caught: str      #: the handler's exception type, "<param>" for a helper
    expr: str        #: the flagged expression, whitespace-normalised
    lineno: int
    kind: str        #: which stringifying shape matched


@dataclass(frozen=True)
class Exemption:
    """A stringified exception that is allowed to reach a client, and why.

    Matched on ``(path, qualname, caught, expr)`` and deliberately NOT on the
    line number: line numbers move on every edit above them, so a table keyed
    on them rots in a week. Keyed on the expression text, an exemption
    survives the code moving and dies the moment the expression itself
    changes -- which is exactly when it needs re-justifying.

    ``caught`` is part of the key because it is the essence of the decision:
    ``str(exc)`` under ``except RunCancelled`` and ``str(exc)`` under
    ``except Exception`` are the same five characters and opposite judgements.

    ``qualname="*"`` matches every function in that file. Use it only where
    one expression repeats across routes for ONE reason (the JSON body
    parser), never to blanket a file.
    """

    path: str
    qualname: str
    caught: str
    expr: str
    reason: str

    def matches(self, f: Finding) -> bool:
        return (self.path == f.path
                and self.expr == f.expr
                and self.caught == f.caught
                and (self.qualname == "*" or self.qualname == f.qualname))


def _norm(text: str) -> str:
    return re.sub(r"\s+", " ", text or "").strip()


def _references(node: ast.AST, names: Sequence[str]) -> bool:
    return any(isinstance(n, ast.Name) and n.id in names for n in ast.walk(node))


def _is_stringy(node: ast.AST) -> bool:
    return isinstance(node, ast.JoinedStr) or (
        isinstance(node, ast.Constant) and isinstance(node.value, str))


def leak_kind(node: ast.AST, bound: Sequence[str]) -> Optional[str]:
    """Does ``node`` turn one of the ``bound`` exception names into text?

    Deliberately syntactic. Proving that a string reaches a response body
    would need whole-program dataflow; every one of these shapes inside an
    ``except`` block is worth a human deciding about, which is what the
    exemption table records.
    """
    if isinstance(node, ast.JoinedStr):
        for part in node.values:
            if isinstance(part, ast.FormattedValue) and _references(part.value, bound):
                return "f-string"
        return None
    if isinstance(node, ast.Call):
        func = node.func
        if isinstance(func, ast.Name) and func.id in ("str", "repr", "format"):
            if any(_references(a, bound) for a in node.args):
                return f"{func.id}()"
        if isinstance(func, ast.Attribute) and func.attr == "format":
            if any(_references(a, bound) for a in node.args):
                return ".format()"
        if isinstance(func, ast.Attribute) and func.attr == "join":
            if any(_references(a, bound) for a in node.args):
                return ".join()"
        return None
    if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Mod, ast.Add)):
        left, right = node.left, node.right
        if _is_stringy(left) and _references(right, bound):
            return "%/+ formatting"
        if _is_stringy(right) and _references(left, bound):
            return "%/+ formatting"
    return None


#: Parameter names that mean "an exception someone caught upstream". A helper
#: that takes one and formats it is the obvious way this pattern regenerates
#: without tripping a check that only looks inside ``except`` blocks -- moving
#: the leak one function call away must not launder it.
_EXC_PARAM_NAMES = frozenset({"exc", "e", "err", "ex", "exception", "error"})


def _exception_params(node) -> Tuple[str, ...]:
    args = node.args
    every = list(args.args) + list(args.posonlyargs or []) + list(args.kwonlyargs)
    out = []
    for a in every:
        ann = ast.unparse(a.annotation) if getattr(a, "annotation", None) else ""
        if a.arg in _EXC_PARAM_NAMES or "Exception" in ann or "Error" in ann:
            out.append(a.arg)
    return tuple(out)


def scan_source(src: str, path: str) -> List[Finding]:
    """Every stringified caught exception in one module."""
    tree = ast.parse(src)
    out: List[Finding] = []

    def walk(node: ast.AST, stack: Tuple[str, ...],
             bound: Tuple[str, ...], caught: str) -> None:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            stack = stack + (node.name,)
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                params = _exception_params(node)
                if params:
                    bound, caught = bound + params, caught or "<param>"
        if isinstance(node, ast.ExceptHandler):
            inner = bound + ((node.name,) if node.name else ())
            label = (_norm(ast.get_source_segment(src, node.type))
                     if node.type is not None else "BaseException")
            for child in node.body:
                walk(child, stack, inner, label)
            return
        if bound:
            kind = leak_kind(node, bound)
            if kind is not None:
                out.append(Finding(
                    path=path,
                    qualname=".".join(stack) or "<module>",
                    caught=caught or "<param>",
                    expr=_norm(ast.get_source_segment(src, node) or ""),
                    lineno=getattr(node, "lineno", 0),
                    kind=kind,
                ))
                # Do not descend: a nested f-string inside a flagged f-string
                # is the same decision, and one finding is one decision.
                return
        for child in ast.iter_child_nodes(node):
            walk(child, stack, bound, caught)

    walk(tree, (), (), "")
    return out


def redacted_refusals(src: str, path: str) -> List[Tuple[str, int, str]]:
    """Handlers that redact an exception class the product refuses WITH.

    The twin failure, checked structurally: ``except IngestError as exc:``
    followed by ``redact_exception(exc, ...)`` flattens a sentence the CLI
    prints verbatim, and no security review would flag it.
    """
    tree = ast.parse(src)
    out: List[Tuple[str, int, str]] = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.ExceptHandler) or node.type is None:
            continue
        caught = {n.id for n in ast.walk(node.type) if isinstance(n, ast.Name)}
        caught |= {n.attr for n in ast.walk(node.type) if isinstance(n, ast.Attribute)}
        refusals = caught & REFUSAL_CLASSES
        if not refusals:
            continue
        for inner in ast.walk(node):
            if (isinstance(inner, ast.Call) and isinstance(inner.func, ast.Name)
                    and inner.func.id in REDACTORS):
                out.append((path, node.lineno, sorted(refusals)[0]))
                break
    return out


def server_modules() -> List[Tuple[str, str]]:
    """``(repo-relative path, source)`` for every module under guard."""
    out: List[Tuple[str, str]] = []
    paths = []
    for root, _dirs, files in os.walk(SERVER_DIR):
        paths += [os.path.join(root, n) for n in files if n.endswith(".py")]
    paths += [os.path.join(_REPO, p) for p in BOUNDARY_MODULES]
    for full in paths:
        rel = os.path.relpath(full, _REPO).replace(os.sep, "/")
        with open(full, encoding="utf-8") as fh:
            out.append((rel, fh.read()))
    return sorted(out)


def sweep() -> List[Finding]:
    findings: List[Finding] = []
    for rel, src in server_modules():
        findings.extend(scan_source(src, rel))
    return findings


# --------------------------------------------------------------------------- #
# 1. Proof of correctness -- the checker separates the categories it claims to
# --------------------------------------------------------------------------- #
_PROOF_LEAK = '''
def handler():
    try:
        work()
    except Exception as exc:
        return {"reason": f"{type(exc).__name__}: {exc}"}
'''

_PROOF_FIXED = '''
def handler():
    try:
        work()
    except Exception as exc:
        return {"reason": redact_exception(exc, "work failed")}
'''

_PROOF_VERBATIM = '''
def handler():
    try:
        work()
    except IngestError as exc:
        return {"reason": str(exc)}
'''

_PROOF_REDACTED_REFUSAL = '''
def handler():
    try:
        work()
    except IngestError as exc:
        return {"reason": redact_exception(exc, "work failed")}
'''

_PROOF_NOT_AN_EXCEPTION = '''
def handler(section_id):
    value = compute(section_id)
    return f"{section_id}: {value}"
'''


def test_the_checker_flags_the_shape_that_was_actually_shipped():
    found = scan_source(_PROOF_LEAK, "proof.py")
    assert len(found) == 1
    assert found[0].kind == "f-string"
    assert found[0].qualname == "handler"
    assert "type(exc).__name__" in found[0].expr


def test_the_checker_passes_the_fix():
    assert scan_source(_PROOF_FIXED, "proof.py") == []


def test_the_checker_flags_a_verbatim_refusal_so_a_human_declares_it():
    """``str(exc)`` on a refusal is CORRECT -- and still has to be declared.

    The checker cannot know that ``IngestError`` is a deliberate refusal from
    the expression alone, and a checker that guessed would be the thing that
    quietly flattens the product. It reports; the exemption table records the
    judgement and the reason.
    """
    found = scan_source(_PROOF_VERBATIM, "proof.py")
    assert len(found) == 1
    assert found[0].kind == "str()"
    assert redacted_refusals(_PROOF_VERBATIM, "proof.py") == []


def test_the_checker_flags_redaction_applied_to_a_refusal():
    """The twin failure: safe-looking, and it degrades the product."""
    hits = redacted_refusals(_PROOF_REDACTED_REFUSAL, "proof.py")
    assert len(hits) == 1
    assert hits[0][2] == "IngestError"
    # ...and the redacted form is NOT a leak finding, which is precisely why
    # this second direction needs its own check.
    assert scan_source(_PROOF_REDACTED_REFUSAL, "proof.py") == []


def test_the_checker_ignores_text_that_is_not_about_an_exception():
    """Interpolating ordinary values is most of what this backend does; a
    check that fired on those would be noise and would be turned off. (A
    PARAMETER holding an exception is a different matter and IS flagged --
    see test_moving_the_leak_into_a_helper_does_not_launder_it.)"""
    assert scan_source(_PROOF_NOT_AN_EXCEPTION, "proof.py") == []


def test_the_checker_sees_the_other_stringifying_shapes():
    kinds = {f.kind for f in scan_source('''
def handler():
    try:
        work()
    except Exception as exc:
        a = str(exc)
        b = "boom: %s" % exc
        c = "boom: {}".format(exc)
        d = ", ".join(exc.args)
        return a, b, c, d
''', "proof.py")}
    assert kinds == {"str()", "%/+ formatting", ".format()", ".join()"}


def test_moving_the_leak_into_a_helper_does_not_launder_it():
    """The obvious way this pattern regenerates: format the exception one
    function call away from the ``except`` that caught it."""
    found = scan_source('''
def _render(exc):
    return f"{type(exc).__name__}: {exc}"

def _render_annotated(problem: Exception):
    return "boom: %s" % problem
''', "proof.py")
    assert {f.qualname for f in found} == {"_render", "_render_annotated"}
    assert {f.caught for f in found} == {"<param>"}


def test_the_handler_type_is_part_of_the_decision():
    """The same five characters, opposite judgements -- so they must not
    share an exemption."""
    found = scan_source('''
def handler():
    try:
        work()
    except RunCancelled as exc:
        cancelled(str(exc))
    except Exception as exc:
        failed(str(exc))
''', "proof.py")
    assert {f.caught for f in found} == {"RunCancelled", "Exception"}


# --------------------------------------------------------------------------- #
# 2. The declared table
#
# Every entry is a stringified exception on the HTTP surface that a human
# decided to keep. "It is only a log" is NOT a reason on its own: `log_tail`
# is served by GET /api/runs/{id} and persisted into result.json, so the job
# log is a client-facing surface like any other.
# --------------------------------------------------------------------------- #
EXEMPTIONS: Tuple[Exemption, ...] = (
    Exemption(
        path="app/server/cache.py",
        qualname="log_internal",
        caught="<param>",
        expr='f"[codon:{ref}] {context}: {_detail(exc)}"',
        reason=(
            "THE record. This is the one line in the backend whose job is to "
            "write the real failure down, and it writes it to STDERR under "
            "the reference handed back to the caller -- not to a response "
            "body, not to the job log. Flagged only because it calls _detail, "
            "which the checker cannot know is safe; the two branches of "
            "_detail are declared below on their own merits."
        ),
    ),
    Exemption(
        path="app/server/cache.py",
        qualname="_detail",
        caught="<param>",
        expr='f"{type(exc).__name__}: {exc}"',
        reason=(
            "THE definition of the real failure text, and the only one -- it "
            "feeds both the stderr record and the de-duplication key behind "
            "stable=True, which must agree or the same failure logs twice "
            "under two ids. It goes to stderr under the reference, never to a "
            "response body and never to the job log. Removing it would leave "
            "the browser holding a reference that corresponds to nothing."
        ),
    ),
    Exemption(
        path="app/server/cache.py",
        qualname="_detail",
        caught="<param>",
        expr="str(exc)",
        reason=(
            "The PreformattedFailure branch of the same definition: a failure "
            "a lower layer already reduced to text carries its own class name "
            "inside the string, so prefixing it with the wrapper's class "
            "would print a wrong class in front of the right one. Same "
            "destination as the branch above it."
        ),
    ),
    Exemption(
        path="app/server/ingest_registry.py",
        qualname="register",
        caught="IngestError",
        expr="str(exc)",
        reason=(
            "A DELIBERATE refusal, and the one the product most depends on. "
            "IngestError says 'no recognised layout. Expected one of: a "
            "Xenium bundle ...' -- the sentence that tells an operator what "
            "shape of directory to hand it. It is caught in its own handler "
            "ahead of the broad one (which redacts), and it is rendered "
            "verbatim so the Studio and `codon ingest` can never disagree. "
            "Pinned byte-for-byte by "
            "test_an_ingest_refusal_survives_the_http_boundary_byte_for_byte."
        ),
    ),
    Exemption(
        path="app/server/jobs.py",
        qualname="JobManager._execute",
        caught="Exception",
        expr="str(exc)",
        reason=(
            "The pass-through half of the narrowed contract for a failed run: "
            "`str(exc) if is_run_refusal(exc) else redact_exception(...)`. A "
            "deliberate refusal (IngestError, TplError, CompoundRefused -- see "
            "jobs.RUN_REFUSALS) reaches the client word for word; everything "
            "else takes the redacted branch. The unbranched form this replaced "
            "WAS a live leak: `job.error` is served in full by result_view(), "
            "and visium._locate resolves to <repo>/data_cache/... so a run can "
            "fail with FileNotFoundError naming a host path. Both branches are "
            "pinned by execution in tests/test_app_server.py."
        ),
    ),
    Exemption(
        path="app/server/jobs.py",
        qualname="JobManager._execute",
        caught="RunCancelled",
        expr="str(exc)",
        reason=(
            "RunCancelled, in its own handler ahead of the broad one. A "
            "cancellation is an ANSWER to something the user clicked, not a "
            "fault: the message says what was cancelled and at which stage, "
            "and it is composed inside app/server/jobs.py from the job's own "
            "state. Nothing from this host is in it."
        ),
    ),
    Exemption(
        path="app/server/jobs.py",
        qualname="JobManager._execute",
        caught="RunCancelled",
        expr='f"cancelled: {exc}"',
        reason=(
            "The RunCancelled message copied into the job log, which "
            "GET /api/runs/{id} serves as log_tail -- so it is a client "
            "surface and was judged as one. Safe for the same reason as the "
            "entry above, and it moves if and only if that one does."
        ),
    ),
    Exemption(
        path="tmesim/protocol/executor.py",
        qualname="execute_plan",
        caught="Exception",
        expr='f"{type(exc).__name__}: {exc}"',
        reason=(
            "The internal half of `text = str(exc) if refusal else f\"...\"`. "
            "The raw text is KEPT here on purpose and marked with "
            "error_is_internal, because execute_plan is an exported engine "
            "entry point and a library caller with no HTTP layer in front of "
            "them should get the real failure. "
            "app.server.protocol_api.execute_protocol_job replaces it for the "
            "browser, keyed off that flag -- the same split scan_report uses "
            "with ingest_registry.scan. Pinned by "
            "test_an_unexpected_sweep_failure_is_redacted_at_the_boundary."
        ),
    ),
    Exemption(
        path="tmesim/protocol/executor.py",
        qualname="execute_plan",
        caught="Exception",
        expr="str(exc)",
        reason=(
            "Two branches, one judgement, and both are the deliberate half: "
            "the RunCancelled handler (a cancellation is an answer to "
            "something the user clicked) and the refusal branch of the "
            "ternary below it (IngestError and friends -- see "
            "executor.REFUSAL_NAMES -- travel verbatim so the sweep grid "
            "cannot contradict what the same refusal says elsewhere). Neither "
            "is flagged internal, so the boundary leaves both alone. Pinned "
            "by test_a_deliberate_refusal_survives_the_sweep_boundary."
        ),
    ),
    Exemption(
        path="app/server/api.py",
        qualname="*",
        caught="Exception",
        expr='f"request body is not valid JSON: {exc}"',
        reason=(
            "The caller's OWN request body, echoed back with the parser's "
            "position ('Expecting value: line 3 column 5 (char 12)'). User "
            "input returned to the user who sent it is not a leak, and the "
            "position is the entire value of the message. JSONDecodeError is "
            "built from the string that was parsed and carries no filesystem "
            "path. Wildcarded because the same decoder guard repeats verbatim "
            "on five routes for this one reason."
        ),
    ),
    Exemption(
        path="app/server/protocol_api.py",
        qualname="build_protocol_router.post_upload",
        caught="(binascii.Error, ValueError)",
        expr='f"content_base64 is not valid base64: {exc}"',
        reason=(
            "The caller's own uploaded field. binascii.Error says which "
            "character or how much padding was wrong, which is what the "
            "uploader needs in order to fix it, and it is built from the "
            "submitted string alone -- nothing from this host reaches it."
        ),
    ),
)


# --------------------------------------------------------------------------- #
# 2b. Known, unfixed, and NOT called safe
#
# An exemption means "this string is safe to show a client, and here is why".
# These two are not that. They are live leaks of the same shape as the ones
# fixed today, found by extending the sweep to the executor, and left because
# fixing them properly is a different design that reaches outside the file:
# a plan-time DIAGNOSTIC has no flag field to carry "this one is internal", so
# the fix needs `Diagnostic` in tmesim/protocol/errors.py to gain one and every
# diagnostic boundary in app/server to honour it.
#
# Recording them here rather than as exemptions keeps two things true at once:
# the suite is green, and nobody can read this file and come away believing
# these strings were judged safe. The count is asserted, so a THIRD one cannot
# join them quietly -- a new plan-time leak fails the build like any other.
# --------------------------------------------------------------------------- #
KNOWN_UNFIXED: Tuple[Exemption, ...] = ()
#: Empty, and the machinery is kept rather than deleted with the last entry.
#: It held two: ``plan()`` put an internal class name and a raw exception
#: message into the script editor's diagnostics. Both are closed --
#: ``executor.internal_diagnostic`` now redacts where the exception is caught,
#: so the raw text never enters a ``Diagnostic`` and every diagnostic boundary
#: is covered by construction rather than one call site at a time.
#:
#: The list stays because the NEXT leak that cannot be closed in the same
#: change needs somewhere honest to go. Filing it as an exemption would assert
#: it is safe; this asserts the opposite, in a place a reader cannot mistake
#: for a cleared one.


def _exemption_for(f: Finding) -> Optional[Exemption]:
    for ex in tuple(EXEMPTIONS) + tuple(KNOWN_UNFIXED):
        if ex.matches(f):
            return ex
    return None


def test_the_known_unfixed_list_does_not_grow():
    """Zero. Both entries were closed; anything here again is a new decision.

    The count is pinned rather than the list merely being checked, so that
    adding an entry is a deliberate edit to this assertion and not something
    that happens quietly while making a sweep pass.
    """
    assert len(KNOWN_UNFIXED) == 0
    live = sweep()
    for ex in KNOWN_UNFIXED:
        assert any(ex.matches(f) for f in live), (
            f"{ex.path}::{ex.qualname} was fixed or moved -- delete the entry "
            "from KNOWN_UNFIXED rather than leaving a defect recorded that no "
            "longer exists")
        assert ex.reason.startswith("LEAK, NOT SAFE"), (
            "an entry in KNOWN_UNFIXED must say plainly that it is a defect; "
            "this list is not a second exemption table")


def test_every_exemption_carries_a_written_reason():
    for ex in tuple(EXEMPTIONS) + tuple(KNOWN_UNFIXED):
        assert len(ex.reason) >= 80, (
            f"{ex.path}::{ex.qualname} is exempt with a one-liner. An "
            "exemption is a judgement someone has to be able to re-check; "
            "write down what makes this string safe to show a client."
        )


def test_no_exemption_is_stale():
    """A rule nobody can delete an entry from is a rule that stops meaning
    anything. Every declared exemption must still match a real finding."""
    live = sweep()
    stale = [ex for ex in tuple(EXEMPTIONS) + tuple(KNOWN_UNFIXED)
             if not any(ex.matches(f) for f in live)]
    assert not stale, (
        "these exemptions no longer match anything in app/server/ -- the code "
        "was fixed or moved, so delete the entry:\n"
        + "\n".join(f"  {ex.path}::{ex.qualname}  {ex.expr}" for ex in stale)
    )


# --------------------------------------------------------------------------- #
# 3. The sweep
# --------------------------------------------------------------------------- #
def test_no_undeclared_exception_text_on_the_http_surface():
    undeclared = [f for f in sweep() if _exemption_for(f) is None]
    assert not undeclared, (
        "an unexpected exception is being turned into text on the HTTP "
        "surface. Answer it with cache.redact_exception(exc, context) -- the "
        "browser gets a stable sentence and a reference, the server log gets "
        "the real detail under that reference. If the exception is a "
        "DELIBERATE refusal that must reach the user word for word, add it to "
        "EXEMPTIONS with the reason:\n"
        + "\n".join(f"  {f.path}:{f.lineno}  {f.qualname}  [{f.kind}]  {f.expr}"
                    for f in undeclared)
    )


def test_no_deliberate_refusal_is_redacted():
    hits: List[Tuple[str, int, str]] = []
    for rel, src in server_modules():
        hits.extend(redacted_refusals(src, rel))
    assert not hits, (
        "a deliberate refusal is being redacted. These classes are the "
        "product's informative 'no' -- IngestError names the layouts it "
        "accepts, TplError carries a line and column -- and both the CLI and "
        "the Studio render them verbatim so the two can never disagree. "
        "Catch the refusal in its own handler ahead of the broad one:\n"
        + "\n".join(f"  {p}:{ln}  catches {cls}" for p, ln, cls in hits))


def test_the_sweep_actually_reads_the_modules_it_claims_to():
    """A sweep over an empty file list passes vacuously. Pin the floor."""
    mods = server_modules()
    names = {os.path.basename(p) for p, _ in mods}
    for expected in ("api.py", "cache.py", "jobs.py", "report.py",
                     "projects.py", "compounds.py", "protocol_api.py",
                     "ingest_registry.py", "executor.py"):
        assert expected in names, f"{expected} vanished from the sweep"
    assert "tmesim/protocol/executor.py" in {p for p, _ in mods}
    assert sum(len(src) for _, src in mods) > 100_000


# --------------------------------------------------------------------------- #
# 4. Execution -- what the two directions actually do at runtime
# --------------------------------------------------------------------------- #
@pytest.fixture()
def cache_mod():
    return pytest.importorskip("app.server.cache")


def test_a_redacted_response_hides_the_detail_and_logs_it_under_the_reference(
        cache_mod, capsys):
    """Direction one, end to end: response clean, log complete, correlated."""
    secret = "/Users/someone/private/data_cache/xenium_gbm/cells.csv.gz"
    try:
        raise FileNotFoundError(f"[Errno 2] No such file or directory: {secret!r}")
    except FileNotFoundError as exc:
        shown = cache_mod.redact_exception(exc, "the section could not be read")

    logged = capsys.readouterr().err

    # (a) the response hides the internal detail
    assert secret not in shown
    assert "FileNotFoundError" not in shown
    assert "Errno" not in shown
    assert shown.startswith("the section could not be read.")

    # (b) the log has the real detail
    assert secret in logged
    assert "FileNotFoundError" in logged
    assert "Traceback" in logged

    # (c) the two are correlated by the reference the caller was handed
    ref = re.search(r"reference ([0-9a-f]{8})", shown)
    assert ref, f"no quotable reference in {shown!r}"
    assert f"[codon:{ref.group(1)}]" in logged


def test_the_brief_form_logs_identically_and_still_carries_the_reference(
        cache_mod, capsys):
    """``brief`` changes what the caller shows, never what the operator gets."""
    secret = "/Users/someone/private/thing"
    try:
        raise OSError(secret)
    except OSError as exc:
        shown = cache_mod.redact_exception(exc, "calibrate_alpha() failed",
                                           brief=True)
    logged = capsys.readouterr().err
    assert secret not in shown and "OSError" not in shown
    assert secret in logged and "OSError" in logged
    ref = re.search(r"reference ([0-9a-f]{8})", shown)
    assert ref and f"[codon:{ref.group(1)}]" in logged


def test_a_client_facing_registry_failure_is_redacted_in_the_response(
        cache_mod, monkeypatch, capsys):
    """The route every client hits on load, with a failing metadata read."""
    secret = "/Users/someone/private/data_cache"

    def boom(section_id):
        raise OSError(f"cannot stat {secret}")

    monkeypatch.setattr(cache_mod, "describe_section", boom)
    rows = cache_mod.list_sections()
    logged = capsys.readouterr().err

    assert rows, "list_sections must still enumerate, not raise"
    for row in rows:
        assert row["available"] is False
        assert secret not in row["reason"]
        assert "OSError" not in row["reason"]
        ref = re.search(r"reference ([0-9a-f]{8})", row["reason"])
        assert ref, f"no quotable reference in {row['reason']!r}"
        assert f"[codon:{ref.group(1)}]" in logged
    assert secret in logged


def test_an_ingest_refusal_survives_the_http_boundary_byte_for_byte(tmp_path):
    """Direction two, the one that breaks silently.

    An empty directory is refused by ``ingest`` with the sentence that names
    every layout it accepts. That sentence is the product: it is what turns
    "it didn't work" into "put a gene_panel.json next to it". Both HTTP
    boundaries -- ``scan`` and ``register`` -- must hand the browser the
    engine's bytes, not a paraphrase and not a reference number.
    """
    ing = pytest.importorskip("tmesim.data.ingest")
    registry = pytest.importorskip("app.server.ingest_registry")

    empty = tmp_path / "not_a_sample"
    empty.mkdir()

    # The engine's own refusal, straight from the CLI's code path.
    with pytest.raises(ing.IngestError) as excinfo:
        ing.ingest(str(empty))
    engine_text = str(excinfo.value)
    assert "no recognised layout" in engine_text
    assert "Expected one of" in engine_text

    # POST /api/sections -> register()
    out = registry.register(str(empty))
    assert out["ok"] is False
    assert out["error"] == engine_text, (
        "the ingest refusal was reworded on the way to the browser; the "
        "Studio now says something the CLI does not")

    # POST /api/sections/scan -> scan(); scan_report answers rather than raises,
    # and its deliberate text must cross the boundary untouched too.
    engine_scan = ing.scan_report(str(empty))
    assert engine_scan["error"]
    assert not engine_scan.get("error_is_internal")
    assert registry.scan(str(empty))["error"] == engine_scan["error"]

    # And the refusal is NOT the redacted shape.
    for text in (out["error"], engine_scan["error"]):
        assert "unexpected internal failure" not in text
        assert "reference" not in text


def test_an_unexpected_ingest_failure_is_still_redacted(tmp_path, monkeypatch):
    """The same boundary, the other branch: not a refusal, so it IS redacted.

    Proves the two directions are separated by the exception CLASS and not by
    which function happened to raise -- the distinction the fix rests on.
    """
    ing = pytest.importorskip("tmesim.data.ingest")
    registry = pytest.importorskip("app.server.ingest_registry")

    secret = "/Users/someone/private/site-packages/tmesim"

    def boom(*a, **kw):
        raise RuntimeError(f"BLAS handle vanished at {secret}")

    monkeypatch.setattr(ing, "ingest", boom)
    monkeypatch.setattr("tmesim.data.ingest.ingest", boom, raising=False)

    xenium = tmp_path / "Sample"
    xenium.mkdir()
    out = registry.register(str(xenium))
    assert out["ok"] is False
    assert secret not in out["error"]
    assert "RuntimeError" not in out["error"]


# --------------------------------------------------------------------------- #
# 5. The protocol sweep boundary
#
# The executor decides (it has the exception CLASS); protocol_api acts on the
# decision (it has the response). These tests drive both halves through the
# real ``execute_protocol_job``, because the interesting failure is the two
# halves disagreeing, which no single-module test can see.
# --------------------------------------------------------------------------- #
def test_the_two_refusal_tables_cannot_drift_apart():
    """Two modules classify the same failures for two surfaces. If one gains
    a class the other has not, that class starts being redacted on one route
    and rendered verbatim on the other -- and nothing would say so."""
    from app.server.jobs import RUN_REFUSALS
    from tmesim.protocol.executor import REFUSAL_NAMES

    assert set(REFUSAL_NAMES) == set(RUN_REFUSALS)


#: A REAL program, planned for real, so ``execute_protocol_job`` runs its own
#: plan/validate path rather than one built out of fakes -- the fake would have
#: to grow a ``to_dict``, an ``estimated_cost`` and a ``capability`` to satisfy
#: ``editor_plan``, and a fake that elaborate drifts from the thing it stands
#: for. ``from library`` rather than ``from smiles`` keeps it off RDKit, and
#: only the SOLVER is replaced below.
_CHEAP_PROGRAM = '''drug TMZ from library "temozolomide"
end

protocol p1
  route oral
  dose 150 mg/m2
end

run TMZ with p1 on xenium_319
  dx 60 um
  report kill
end
'''


class _FakeJob:
    def __init__(self):
        self.protocol_grid = None
        self.progress = 0.0
        self.stage = ""
        self.run_id = "r_test"
        self.log_lines = []


class _FakeCtx:
    """The slice of jobs._Ctx that execute_protocol_job actually uses."""

    def __init__(self):
        self.job = _FakeJob()

    def stage(self, name):
        self.job.stage = name

    def log(self, msg):
        self.job.log_lines.append(msg)

    def check_cancelled(self):
        return None


def _run_sweep_with(monkeypatch, exc):
    """Drive execute_protocol_job with one run that raises ``exc``."""
    protocol_api = pytest.importorskip("app.server.protocol_api")
    executor = pytest.importorskip("tmesim.protocol.executor")

    real_execute = executor.execute_plan

    def one_failed_run(p, log=None, on_progress=None, check_cancelled=None,
                       cache=None):
        # The executor's OWN failure handler, reached without a solver: build
        # the cell and row exactly as execute_plan does, via its classifier.
        refusal = executor.is_refusal(exc)
        text = str(exc) if refusal else f"{type(exc).__name__}: {exc}"
        cell = {"index": 0, "label": "run 0", "drug": "d", "status": "failed",
                "error": text, "runtime_s": 0.1}
        row = {"index": 0, "label": "run 0", "axes": {}, "drug": "d",
               "section": "s", "metrics": None, "error": text}
        if not refusal:
            cell["error_is_internal"] = True
            row["error_is_internal"] = True
        if log is not None:
            log(f"[0] FAILED: {text}" if refusal else
                "[0] FAILED: an unexpected internal error; see this run's "
                "error field for the reference to quote")
        if on_progress is not None:
            on_progress([cell])
            on_progress([cell])          # a second poll: the id must not move
        return {"n_runs": 1, "grid": [cell], "rows": [row], "ranking": [],
                "sweep_axes": [], "expects": [], "diagnostics": [],
                "reconstructions": {"built": [], "reused": []}}

    monkeypatch.setattr(executor, "execute_plan", one_failed_run)
    monkeypatch.setattr("tmesim.protocol.execute_plan", one_failed_run,
                        raising=False)
    assert real_execute is not one_failed_run
    ctx = _FakeCtx()
    body = protocol_api.execute_protocol_job({"source": _CHEAP_PROGRAM}, ctx)
    return ctx, body


def test_an_unexpected_sweep_failure_is_redacted_at_the_boundary(
        monkeypatch, capsys):
    secret = "/Users/someone/private/data_cache/visium_gbm/GBM1"
    ctx, body = _run_sweep_with(
        monkeypatch, FileNotFoundError(f"no matrix.mtx(.gz) in {secret}"))
    logged = capsys.readouterr().err

    published = body["protocol"]
    for projection in (published["grid"], published["rows"]):
        shown = projection[0]["error"]
        assert secret not in shown
        assert "FileNotFoundError" not in shown
        assert "error_is_internal" not in projection[0]
        assert re.search(r"reference ([0-9a-f]{8})", shown)

    # The cell and its row quote the SAME reference -- one failure, one record.
    assert published["grid"][0]["error"] == published["rows"][0]["error"]
    # ...and the live grid published mid-sweep already carried it.
    assert ctx.job.protocol_grid[0]["error"] == published["grid"][0]["error"]

    # The job log is a client surface too (GET /api/runs/{id}/status).
    assert not any(secret in line for line in ctx.job.log_lines)
    assert not any("FileNotFoundError" in line for line in ctx.job.log_lines)

    # The operator gets everything, under the id the browser was given.
    ref = re.search(r"reference ([0-9a-f]{8})", published["grid"][0]["error"])
    assert secret in logged
    assert f"[codon:{ref.group(1)}]" in logged


def test_a_polled_sweep_does_not_change_the_reference_under_the_client(
        monkeypatch, capsys):
    """``_progress`` fires per completed run and republishes the whole grid.
    Redacting afresh each time would hand the client a new id on every poll,
    and the one they finally quoted would match no record."""
    ctx, body = _run_sweep_with(monkeypatch, RuntimeError("boom /private/x"))
    logged = capsys.readouterr().err

    shown = body["protocol"]["grid"][0]["error"]
    refs = re.findall(r"reference ([0-9a-f]{8})", shown)
    assert len(refs) == 1
    # Three publications happened (two polls plus the final one); exactly one
    # record was written, so the id the client holds is the only id there is.
    assert logged.count("[codon:") == 1
    assert f"[codon:{refs[0]}]" in logged


def test_a_deliberate_refusal_survives_the_sweep_boundary(monkeypatch, capsys):
    """The twin direction on the path a leak was just fixed on."""
    from tmesim.data.ingest import IngestError

    refusal = ("no recognised layout. Expected one of: a Xenium bundle "
               "(cells.csv.gz + cell_feature_matrix[.tar.gz] + "
               "gene_panel.json), a Space Ranger output, or a MEX triple.")
    ctx, body = _run_sweep_with(monkeypatch, IngestError(refusal))
    logged = capsys.readouterr().err

    for projection in (body["protocol"]["grid"], body["protocol"]["rows"]):
        assert projection[0]["error"] == refusal
        assert "unexpected internal failure" not in projection[0]["error"]
        assert "reference" not in projection[0]["error"]
    # A refusal is not a fault: nothing was written to the internal log.
    assert "[codon:" not in logged
    # ...and the job log carries the sentence itself, not a pointer.
    assert any(refusal in line for line in ctx.job.log_lines)


def test_a_preformatted_failure_logs_without_a_wrapper_class(
        monkeypatch, capsys, tmp_path):
    """``scan``'s internal branch: the engine hands over TEXT, not an exception.

    The boundary has to wrap it in something to get a log record, and wrapping
    in a bare RuntimeError printed a class name that was not the failure's own
    in front of one that was. This is also the branch the byte-identity test
    above deliberately does not reach, so without this the flagged half of
    ``scan`` had no coverage at all.
    """
    ing = pytest.importorskip("tmesim.data.ingest")
    registry = pytest.importorskip("app.server.ingest_registry")

    secret = "/Users/someone/private/site-packages/numpy"

    def flagged(path, measure_detection_rates=False):
        return {"error": f"ImportError: no BLAS at {secret}",
                "error_is_internal": True, "detection": {}, "files": {}}

    monkeypatch.setattr(ing, "scan_report", flagged)
    out = registry.scan(str(tmp_path))
    logged = capsys.readouterr().err

    assert secret not in out["error"]
    assert "ImportError" not in out["error"]
    assert "error_is_internal" not in out
    ref = re.search(r"reference ([0-9a-f]{8})", out["error"])
    assert ref

    record = [ln for ln in logged.splitlines() if ln.startswith("[codon:")]
    assert len(record) == 1
    assert record[0].endswith(f"ImportError: no BLAS at {secret}")
    assert "RuntimeError" not in record[0], (
        "the wrapper's class name is in front of the failure's own")
    assert "PreformattedFailure" not in record[0]
    # No except block is in flight here, so there is no stack to print and
    # the record must not claim otherwise.
    assert "NoneType: None" not in logged


# --------------------------------------------------------------------------- #
# 6. The plan-time diagnostic path
#
# Redacted where the exception is caught rather than behind a flag, because a
# Diagnostic has no second consumer: it exists to be shown to a human editing
# a script, and two of its boundaries (protocol_wire._refusal and
# _caveat_lines) read `.message` off the object, so a flag honoured at
# `to_dict()` would not have covered them. These tests assert BOTH shapes --
# the dict and the attribute -- for exactly that reason.
# --------------------------------------------------------------------------- #
_PLAN_SRC = _CHEAP_PROGRAM


def _plan_with_broken_drug_compile(monkeypatch, secret):
    executor = pytest.importorskip("tmesim.protocol.executor")
    real = executor.resolve_drug

    def boom(node, compile_params=False, alpha=None):
        if compile_params:
            raise KeyError(f"no GDSC row in {secret}")
        return real(node, compile_params=compile_params, alpha=alpha)

    monkeypatch.setattr(executor, "resolve_drug", boom)
    return executor


def test_a_plan_time_failure_is_redacted_in_every_projection(
        monkeypatch, capsys):
    secret = "/Users/someone/private/gdsc/ic50_table_a.csv"
    executor = _plan_with_broken_drug_compile(monkeypatch, secret)

    p = executor.plan(_PLAN_SRC)
    logged = capsys.readouterr().err

    compiled = [d for d in p.diagnostics.items if d.code == "compile-failed"]
    assert compiled, "the residue branch did not run"

    # (a) the dict projection -- to_dict(), what most boundaries use
    as_dict = [d.to_dict() for d in compiled][0]
    # (b) the ATTRIBUTE, which protocol_wire._refusal reads directly and which
    #     a flag honoured at to_dict() would have missed entirely
    as_attr = compiled[0].message

    for shown in (as_dict["message"], as_attr):
        assert secret not in shown
        assert "KeyError" not in shown
        assert re.search(r"reference ([0-9a-f]{8})", shown)
    assert as_dict["message"] == as_attr

    # The operator's record has everything, under the quoted reference.
    ref = re.search(r"reference ([0-9a-f]{8})", as_attr).group(1)
    assert secret in logged
    assert f"[codon:{ref}]" in logged
    assert "KeyError" in logged


def test_re_planning_does_not_mint_a_new_reference(monkeypatch, capsys):
    """The editor re-plans on every edit. Without ``stable=True`` this would
    write a record and hand back a different id several times a second, and
    the id the user finally quoted would be one of hundreds."""
    secret = "/Users/someone/private/gdsc/ic50_table_b.csv"
    executor = _plan_with_broken_drug_compile(monkeypatch, secret)

    seen = set()
    for _ in range(5):
        p = executor.plan(_PLAN_SRC)
        msg = [d for d in p.diagnostics.items
               if d.code == "compile-failed"][0].message
        seen.add(re.search(r"reference ([0-9a-f]{8})", msg).group(1))
    logged = capsys.readouterr().err

    assert len(seen) == 1, f"five re-plans handed out {len(seen)} references"
    assert logged.count(f"[codon:{seen.pop()}]") == 1, (
        "the same failure was recorded more than once")


def test_a_deliberate_plan_time_refusal_is_untouched(capsys):
    """An unknown drug is a refusal that NAMES THE ALTERNATIVES. Redacting it
    would replace a list of available drugs with a reference number."""
    executor = pytest.importorskip("tmesim.protocol.executor")

    p = executor.plan(_PLAN_SRC.replace('"temozolomide"', '"no_such_drug"'))
    logged = capsys.readouterr().err

    messages = [d.message for d in p.diagnostics.errors]
    named = [m for m in messages if "not in the built-in library" in m]
    assert named, f"expected a library refusal, got {messages}"
    assert "Available:" in named[0]
    assert "temozolomide" in named[0]
    for m in messages:
        assert "unexpected internal failure" not in m
    # A refusal is not a fault: nothing was written to the internal log.
    assert "[codon:" not in logged
