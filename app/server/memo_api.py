"""``POST /api/memo`` -- the HTTP surface of :mod:`app.server.memo`.

WHY THIS IS A SECOND FILE AND NOT THE BOTTOM OF ``memo.py``
-----------------------------------------------------------
``app/server/jobs.py`` imports ``memo`` to stamp every run with its key, and
``jobs`` runs wherever the engine runs -- the CLI, ``experiments/``, a fresh
checkout. ``fastapi`` is in the OPTIONAL ``[app]`` extra of ``pyproject.toml``
("fastapi>=0.115", declared beside uvicorn and httpx, absent from
``constraints.txt``), so an install without that extra can run the engine and
cannot import fastapi. A module-level ``from fastapi import ...`` inside
``memo.py`` would therefore make the memo stamp -- and with it every run --
depend on the web stack.

Keeping the import inside ``build_memo_router`` is the other obvious fix and it
does not work: ``memo.py`` uses ``from __future__ import annotations``, so a
route handler's ``request: Request`` annotation is a STRING that FastAPI
resolves against the handler's MODULE globals. A ``Request`` imported into a
function's locals is not there, and the route comes back 422 asking for a
query parameter named ``request``. Measured, not assumed: that is exactly what
the first version of this route did.

So the split is along the dependency, which is also where this repository
already splits it -- ``browse.py`` and ``workspace.py`` both import fastapi at
module level and are both imported lazily from inside ``create_app``.

THE MOUNT, which is one line inside ``create_app`` in ``app/server/api.py``::

    app.include_router(build_memo_router(jobs=jobs))

PLACEMENT IS LOAD-BEARING, for the reason that file's own mounting block gives:
``include_router`` appends and ``_mount_web`` mounts StaticFiles at ``/``,
which matches every path, so the line must sit BEFORE ``_mount_web(app)``.
Next to ``app.include_router(build_workspace_router(...))`` is the right place.
``tests/test_memo.py::TestRoute`` asserts both halves -- that the route answers
when it is registered ahead of the web mount, and that it is never reached when
it is appended after.
"""
from __future__ import annotations

from typing import Dict

from fastapi import APIRouter, Request

from .jobs import jsonable
from .memo import _DISCLOSURE, Unmemoisable, answer_spec

__all__ = ["build_memo_router"]


def build_memo_router(*, jobs) -> APIRouter:
    """Mount ``POST /api/memo``. One route, no state, nothing mutated.

    Gated by the same middleware as every other ``/api/`` route: the gate in
    ``create_app`` tests whether the route path starts with ``/api/``, so a
    router mounted there is authenticated by construction and cannot be
    forgotten.

    This route deliberately does NOT hook ``POST /api/runs``. A memo on the
    submit path would make a second click return the first click's ``run_id``,
    and changing what a shipped surface records is not a decision a cache gets
    to make on its own. Wiring it in is a few lines in ``JobManager.submit``
    whenever that is wanted, and it is a product decision rather than a
    correctness one.
    """
    router = APIRouter()

    @router.post("/api/memo")
    async def post_memo(request: Request) -> Dict:
        """Answer a RunSpec without solving it, or say that nothing can.

        Body is either a bare RunSpec or
        ``{"spec": <RunSpec>, "allow_scaling": <bool>}``. A RunSpec has no
        ``spec`` key, so the two forms are unambiguous.

        ``allow_scaling`` defaults to FALSE. A caller that asked whether the
        answer is already known, and was handed an algebraic rearrangement of a
        different run instead, would have been told something it did not ask
        about -- and that distinction is the one thing this surface must never
        blur.

        Always 200. A miss is an ANSWER ("nothing here answers this; submit
        it"), not an error, and every shape of non-answer carries the same
        ``mode``/``reason``/``disclosure`` fields so a client has one branch to
        write. The spec is put through ``/api/validate``'s own validator
        first, so this surface can never answer a spec the run surface would
        reject.
        """
        from .api import _json_body, validate_spec

        body = await _json_body(request)
        if "spec" in body:
            spec = body.get("spec")
            allow = bool(body.get("allow_scaling", False))
        else:
            spec, allow = body, False
        if not isinstance(spec, dict):
            return jsonable({
                "mode": "miss", "is_solve": False,
                "reason": ("the request body must be a RunSpec object, or an "
                           "object with a 'spec' key holding one."),
                "disclosure": _DISCLOSURE["miss"],
            })
        verdict = validate_spec(spec)
        if not verdict.get("ok"):
            return jsonable({
                "mode": "miss", "is_solve": False, "spec": spec,
                "reason": ("this RunSpec is rejected by /api/validate, so "
                           "there is nothing to recall and nothing to scale. "
                           "Its errors are below."),
                "validate": verdict,
                "disclosure": _DISCLOSURE["miss"],
            })
        try:
            answer = answer_spec(jobs, spec, allow_scaling=allow)
        except Unmemoisable as refusal:
            # The sentence is one ``memo.py`` authored, for this caller, about
            # the caller's own spec. It names no path, no exception class and
            # no internal state -- see ``Unmemoisable``'s docstring. Nothing
            # from the exception is stringified here.
            return jsonable({
                "mode": "miss", "is_solve": False, "spec": spec,
                "reason": (refusal.args[0] if refusal.args else
                           "this RunSpec cannot be given a comparable key."),
                "disclosure": _DISCLOSURE["miss"],
            })
        answer["validate_warnings"] = verdict.get("warnings") or []
        return jsonable(answer)

    return router
