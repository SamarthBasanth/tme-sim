"""``/api/protocol/*`` — the TPL routes (see API_CONTRACT.md).

===========  ================================  =================================
Method       Path                              Body -> Returns
===========  ================================  =================================
POST         /api/protocol/validate            {source} -> {ok, diagnostics[]}
POST         /api/protocol/plan                {source} -> the expanded run list
POST         /api/protocol/run                 {source} -> {run_id}
GET          /api/protocol/runs/{id}           results + the sweep grid
POST         /api/protocol/upload              a drug file -> {path, preview}
===========  ================================  =================================

Two things this module deliberately does NOT do:

* **It does not own a queue.** ``POST /api/protocol/run`` submits to the SAME
  :class:`app.server.jobs.JobManager` the ``/api/runs`` path uses, so one worker
  thread serialises every simulation in the process and the reconstruction cache
  is shared. A sweep over dose on one specimen reconstructs once — and once,
  process-wide, even if a plain ``/api/runs`` job wanted the same key.
* **It does not re-implement validation.** Every check lives in
  :mod:`tmesim.protocol`; these routes are transport.

A failure a user caused comes back as ``{ok: false, diagnostics: [...]}`` with a
line and a column on every entry — never a 500 with a traceback.
"""
from __future__ import annotations

import base64
import binascii
import json
import os
import re
import time
from typing import Dict, List, Optional

from fastapi import APIRouter, HTTPException, Request

from .jobs import JobManager, jsonable

__all__ = ["build_protocol_router", "execute_protocol_job", "PROTOCOL_JOB_KIND"]

#: ``spec['kind']`` that marks a job as a TPL sweep rather than a RunSpec.
#: ``app.server.jobs.execute_run`` dispatches on it.
PROTOCOL_JOB_KIND = "protocol"

#: Extensions ``POST /api/protocol/upload`` accepts.
_ALLOWED_UPLOAD_EXT = (".smi", ".smiles", ".mol", ".sdf", ".json", ".yaml", ".yml")

#: Largest upload accepted, decoded. A structure file is kilobytes; anything
#: larger is not a drug definition.
_MAX_UPLOAD_BYTES = 2 * 1024 * 1024

_SAFE_NAME = re.compile(r"[^A-Za-z0-9._-]+")


# --------------------------------------------------------------------------- #
# The job body (called from app.server.jobs.execute_run)
# --------------------------------------------------------------------------- #
def execute_protocol_job(spec: Dict, ctx, cache=None) -> Dict:
    """Run one TPL program as a single queued job.

    A whole sweep is ONE job, not N: the runs share a reconstruction and must be
    serialised anyway, and the UI wants one progress grid rather than N poll
    targets. Per-run status is published live on ``ctx.job.protocol_grid`` so
    ``GET /api/protocol/runs/{id}`` can show the grid while it fills.
    """
    from tmesim.protocol import execute_plan, plan
    from tmesim.protocol.errors import TplError

    source = spec.get("source") or ""
    max_runs = int(spec.get("max_runs") or 0) or None
    max_window_s = spec.get("max_window_s")

    ctx.stage("plan")
    ctx.job.protocol_grid = []
    kwargs: Dict = {}
    if max_runs:
        kwargs["max_runs"] = max_runs
    if max_window_s:
        kwargs["max_window_s"] = float(max_window_s)
    p = plan(source, **kwargs)
    ctx.log(f"planned {p.n_runs} run(s); "
            f"{p.estimated_cost.get('n_reconstructions')} reconstruction(s)")
    if not p.ok:
        raise TplError(
            "the protocol was rejected at plan time: "
            + "; ".join(d.message for d in p.diagnostics.errors),
            diagnostics=list(p.diagnostics.items))

    n = max(p.n_runs, 1)
    done = {"i": 0}

    def _progress(grid: List[Dict]) -> None:
        ctx.job.protocol_grid = jsonable(grid)
        finished = sum(1 for c in grid if c["status"] in
                       ("done", "failed", "cancelled"))
        if finished != done["i"]:
            done["i"] = finished
            ctx.job.progress = 0.05 + 0.95 * (finished / n)
            ctx.job.stage = f"run {finished}/{n}"

    ctx.stage("execute")
    result = execute_plan(p, log=ctx.log, on_progress=_progress,
                          check_cancelled=ctx.check_cancelled, cache=cache)
    ctx.job.protocol_grid = jsonable(result["grid"])
    ctx.stage("done")
    return {
        "protocol": jsonable({
            "source": source,
            "plan": p.to_dict(),
            **result,
        }),
        # Keys the generic RunResult view expects, so a protocol job is not a
        # broken-looking entry in GET /api/runs.
        "report": None,
        "metrics": None,
        "provenance": {"kind": PROTOCOL_JOB_KIND,
                       "n_runs": result["n_runs"],
                       "reconstructions": result["reconstructions"]},
        "_payloads": {},
    }


# --------------------------------------------------------------------------- #
# Routes
# --------------------------------------------------------------------------- #
def build_protocol_router(jobs: JobManager) -> APIRouter:
    """Build the ``/api/protocol`` router bound to an existing job manager."""
    router = APIRouter(prefix="/api/protocol", tags=["protocol"])

    @router.post("/validate")
    async def post_validate(request: Request) -> Dict:
        """Parse + static checks only. No compilation, no compute."""
        from tmesim.protocol import validate_source
        from tmesim.protocol.errors import TplError

        source = await _source(request)
        try:
            return jsonable(validate_source(source))
        except TplError as exc:
            return {"ok": False, "diagnostics": exc.to_dicts()}
        except Exception as exc:  # noqa: BLE001 - never a bare traceback
            return {"ok": False, "diagnostics": [_internal(exc)]}

    @router.post("/plan")
    async def post_plan(request: Request) -> Dict:
        """The expanded run list, run count and estimated cost. Executes nothing."""
        from tmesim.protocol import plan as make_plan
        from tmesim.protocol.errors import TplError

        body = await _body(request)
        source = _source_from(body)
        kwargs: Dict = {}
        if body.get("max_runs"):
            kwargs["max_runs"] = int(body["max_runs"])
        if body.get("max_window_s"):
            kwargs["max_window_s"] = float(body["max_window_s"])
        t0 = time.time()
        try:
            p = make_plan(source, compile_drugs=True, **kwargs)
        except TplError as exc:
            return {"ok": False, "diagnostics": exc.to_dicts(), "runs": [],
                    "n_runs": 0}
        except Exception as exc:  # noqa: BLE001
            return {"ok": False, "diagnostics": [_internal(exc)], "runs": [],
                    "n_runs": 0}
        out = p.to_dict()
        out["plan_time_s"] = round(time.time() - t0, 3)
        return jsonable(out)

    @router.post("/run", status_code=202)
    async def post_run(request: Request) -> Dict:
        """Plan, then enqueue on the SHARED job queue. Returns immediately."""
        from tmesim.protocol import plan as make_plan
        from tmesim.protocol.errors import TplError

        body = await _body(request)
        source = _source_from(body)
        kwargs: Dict = {}
        if body.get("max_runs"):
            kwargs["max_runs"] = int(body["max_runs"])
        if body.get("max_window_s"):
            kwargs["max_window_s"] = float(body["max_window_s"])
        try:
            # Cheap pass first so a bad script is a 400 with positions rather
            # than a job that fails 200 ms after being accepted.
            p = make_plan(source, compile_drugs=False, **kwargs)
        except TplError as exc:
            raise HTTPException(status_code=400, detail={
                "message": "protocol rejected", "ok": False,
                "diagnostics": exc.to_dicts()}) from None
        except Exception as exc:  # noqa: BLE001
            raise HTTPException(status_code=400, detail={
                "message": "protocol rejected", "ok": False,
                "diagnostics": [_internal(exc)]}) from None
        if not p.ok:
            raise HTTPException(status_code=400, detail={
                "message": "protocol rejected", "ok": False,
                "diagnostics": p.diagnostics.to_dicts()})

        spec = {"kind": PROTOCOL_JOB_KIND, "source": source}
        for key in ("max_runs", "max_window_s"):
            if body.get(key):
                spec[key] = body[key]
        job = jobs.submit(spec)
        return jsonable({
            "run_id": job.run_id,
            "status": job.status,
            "n_runs": p.n_runs,
            "max_runs": p.max_runs,
            "allow_runs": p.allow_runs,
            "warnings": [d.to_dict() for d in p.diagnostics.warnings],
            "estimated_cost": p.estimated_cost,
        })

    @router.get("/runs/{run_id}")
    def get_protocol_run(run_id: str) -> Dict:
        """Results + the sweep grid. Safe to poll while the job is running."""
        job = jobs.get(run_id)
        if job is None:
            raise HTTPException(status_code=404,
                                detail=f"no such run {run_id!r}")
        if job.spec.get("kind") != PROTOCOL_JOB_KIND:
            raise HTTPException(status_code=409, detail={
                "message": f"run {run_id!r} is not a protocol run",
                "kind": job.spec.get("kind"),
                "reason": "use GET /api/runs/{id} for a RunSpec run"})
        body = (job.body or {}).get("protocol") or {}
        grid = body.get("grid") or getattr(job, "protocol_grid", None) or []
        return jsonable({
            "run_id": run_id,
            "status": job.status,
            "stage": job.stage,
            "progress": round(float(job.progress), 4),
            "elapsed_s": round(job.elapsed_s, 3),
            "error": job.error,
            "source": job.spec.get("source"),
            "n_runs": body.get("n_runs", len(grid)),
            "grid": grid,
            "rows": body.get("rows", []),
            "ranking": body.get("ranking"),
            "sweep_axes": body.get("sweep_axes", []),
            "expects": body.get("expects", []),
            "diagnostics": body.get("diagnostics", []),
            "reconstructions": body.get("reconstructions"),
            "plan": body.get("plan"),
            "log_tail": job.log_tail(40),
        })

    @router.post("/upload")
    async def post_upload(request: Request) -> Dict:
        """Store a drug file and return its path plus a parsed preview.

        Body (JSON, so no multipart dependency is added to this checkout)::

            {"filename": "compound_47.sdf",
             "content": "<text>"}            # or "content_base64": "<b64>"

        The returned ``path`` is exactly what a ``drug X from file "<path>"``
        statement should use.
        """
        from tmesim.protocol.errors import TplError
        from tmesim.protocol.resolve import load_drug_file, uploads_dir

        body = await _body(request)
        raw_name = str(body.get("filename") or "").strip()
        if not raw_name:
            raise HTTPException(status_code=400, detail={
                "message": "filename is required",
                "reason": "the extension selects the parser"})
        name = _SAFE_NAME.sub("_", os.path.basename(raw_name))
        ext = os.path.splitext(name)[1].lower()
        if ext not in _ALLOWED_UPLOAD_EXT:
            raise HTTPException(status_code=400, detail={
                "message": f"unsupported file type {ext or '(none)'}",
                "allowed": list(_ALLOWED_UPLOAD_EXT),
                "reason": ("structures are .smi/.mol/.sdf; parameter files are "
                           ".json/.yaml")})

        if body.get("content_base64") is not None:
            try:
                data = base64.b64decode(str(body["content_base64"]),
                                        validate=True)
            except (binascii.Error, ValueError) as exc:
                raise HTTPException(status_code=400, detail={
                    "message": f"content_base64 is not valid base64: {exc}"
                }) from None
        elif body.get("content") is not None:
            data = str(body["content"]).encode("utf-8")
        else:
            raise HTTPException(status_code=400, detail={
                "message": "either `content` (text) or `content_base64` is "
                           "required"})
        if len(data) > _MAX_UPLOAD_BYTES:
            raise HTTPException(status_code=413, detail={
                "message": f"file is {len(data)} B, above the "
                           f"{_MAX_UPLOAD_BYTES} B limit",
                "reason": "a drug definition is kilobytes; this is not one"})

        target = os.path.join(uploads_dir(), name)
        tmp = target + ".tmp"
        with open(tmp, "wb") as fh:
            fh.write(data)
        os.replace(tmp, target)

        rel = os.path.join("uploads", name)
        try:
            preview = load_drug_file(rel)
        except TplError as exc:
            os.remove(target)
            raise HTTPException(status_code=400, detail={
                "message": "the uploaded file could not be parsed and was not "
                           "kept",
                "diagnostics": exc.to_dicts()}) from None
        return jsonable({
            "ok": True,
            "path": rel,
            "bytes": len(data),
            "use_in_script": f'drug MYCPD from file "{rel}"',
            "preview": preview,
        })

    return router


# --------------------------------------------------------------------------- #
# Body helpers
# --------------------------------------------------------------------------- #
async def _body(request: Request) -> Dict:
    """Parse the request body as JSON, or wrap raw text as ``{"source": ...}``."""
    raw = await request.body()
    if not raw:
        return {}
    text = raw.decode("utf-8", errors="replace")
    try:
        doc = json.loads(text)
    except Exception:  # noqa: BLE001 - a plain-text script body is allowed
        return {"source": text}
    if isinstance(doc, str):
        return {"source": doc}
    if not isinstance(doc, dict):
        return {"source": text}
    return doc


def _source_from(body: Dict) -> str:
    return str(body.get("source") or body.get("script") or body.get("tpl") or "")


async def _source(request: Request) -> str:
    return _source_from(await _body(request))


def _internal(exc: Exception) -> Dict:
    """An unexpected failure, rendered as a diagnostic instead of a 500."""
    return {
        "severity": "error", "line": 0, "col": 0, "end_col": 0,
        "code": "internal",
        "message": (f"{type(exc).__name__}: {exc}. This is a backend fault, "
                    "not a mistake in the script."),
        "hint": None,
    }
