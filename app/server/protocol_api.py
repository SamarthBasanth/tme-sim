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
from typing import Dict, List

from fastapi import APIRouter, HTTPException, Request

from .cache import (
    PRIVATE_DIR_MODE,
    PRIVATE_FILE_MODE,
    PreformattedFailure,
    own_only,
    redact_exception,
)
from .jobs import JobManager, jsonable
from .protocol_wire import (
    NO_RUN_ID_REASON,
    editor_diagnostics,
    editor_plan,
    editor_validate,
    sweep_snapshot,
)
from .protocol_wire import _refusal as editor_refusal

__all__ = ["build_protocol_router", "execute_protocol_job", "PROTOCOL_JOB_KIND"]

#: ``spec['kind']`` that marks a job as a TPL sweep rather than a RunSpec.
#: ``app.server.jobs.execute_run`` dispatches on it.
PROTOCOL_JOB_KIND = "protocol"

#: Extensions ``POST /api/protocol/upload`` accepts.
_ALLOWED_UPLOAD_EXT = (".smi", ".smiles", ".mol", ".sdf", ".json", ".yaml", ".yml")

#: Largest upload accepted, decoded. A structure file is kilobytes; anything
#: larger is not a drug definition.
_MAX_UPLOAD_BYTES = 2 * 1024 * 1024

#: SECURITY: hard ceiling on the RAW request body for every protocol route,
#: enforced by :func:`_capped_body` before a byte is buffered. Larger than
#: ``_MAX_UPLOAD_BYTES`` because a multipart envelope carries headers and a
#: boundary around the 2 MB payload; small enough that the parser's transient
#: copies stay bounded. There is no reverse proxy in deployment mode 1a, so
#: this is the only body limit in the product.
_MAX_BODY_BYTES = 4 * 1024 * 1024

_SAFE_NAME = re.compile(r"[^A-Za-z0-9._-]+")

#: RFC 2046 §5.1.1: 1-70 characters from this set, not ending in a space.
_BOUNDARY_RE = re.compile(r"^[0-9A-Za-z'()+_,\-./:=?]([0-9A-Za-z'()+_,\-./:=? ]{0,68}"
                          r"[0-9A-Za-z'()+_,\-./:=?])?$")


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
    #: run index -> the sentence a failed cell shows. Memoised because
    #: ``_progress`` fires on every completed run and the grid is re-published
    #: each time: without this, a client polling the sweep would watch the
    #: reference id of an already-failed run change under it, and the number
    #: the operator finally quoted would match no record on the server.
    _shown: Dict[int, str] = {}

    def _published(cells: List[Dict]) -> List[Dict]:
        """The client's copy of the grid (or rows): JSON-safe, and redacted.

        ``jsonable`` deep-copies, so the executor's own structures keep the
        real failure text for a library caller; only what leaves this process
        is rewritten. The executor sets ``error_is_internal`` on exactly the
        failures it could not classify as a deliberate refusal -- it is the
        only frame that still has the exception CLASS to judge by, which is
        why the flag comes from there rather than being guessed at here.
        """
        out = jsonable(cells)
        for cell in out:
            if not cell.pop("error_is_internal", False):
                continue
            index = cell.get("index")
            if index not in _shown:
                _shown[index] = redact_exception(
                    PreformattedFailure(cell.get("error") or "unknown"),
                    f"run {index} of this sweep failed")
            cell["error"] = _shown[index]
        return out

    def _progress(grid: List[Dict]) -> None:
        ctx.job.protocol_grid = _published(grid)
        finished = sum(1 for c in grid if c["status"] in
                       ("done", "failed", "cancelled"))
        if finished != done["i"]:
            done["i"] = finished
            ctx.job.progress = 0.05 + 0.95 * (finished / n)
            ctx.job.stage = f"run {finished}/{n}"

    ctx.stage("execute")
    result = execute_plan(p, log=ctx.log, on_progress=_progress,
                          check_cancelled=ctx.check_cancelled, cache=cache)
    # Both projections of the same failures, through the same memo, so a cell
    # and its row quote the same reference.
    published = dict(result,
                     grid=_published(result["grid"]),
                     rows=_published(result["rows"]))
    ctx.job.protocol_grid = published["grid"]
    ctx.stage("done")
    return {
        "protocol": jsonable({
            "source": source,
            # The editor projection is stored, not just the contract one, so
            # GET /sweeps/{id} can serve axes / coords / caveats / assumptions
            # from a finished job without re-planning it.
            "plan": editor_plan(p.to_dict(), p),
            **published,
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
        from tmesim.protocol import plan as make_plan
        from tmesim.protocol.errors import TplError

        source = await _source(request)
        try:
            p = make_plan(source, compile_drugs=False)
            body = {"ok": p.ok, "diagnostics": p.diagnostics.to_dicts(),
                    "n_runs": p.n_runs, "max_runs": p.max_runs}
            return jsonable(editor_validate(body, p))
        except TplError as exc:
            return _diag_only(exc.to_dicts())
        except Exception as exc:  # noqa: BLE001 - never a bare traceback
            return _diag_only([_internal(exc)])

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
            return _diag_only(exc.to_dicts(), plan_shape=True)
        except Exception as exc:  # noqa: BLE001
            return _diag_only([_internal(exc)], plan_shape=True)
        out = p.to_dict()
        out["plan_time_s"] = round(time.time() - t0, 3)
        return jsonable(editor_plan(out, p))

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
            errors, warnings = editor_diagnostics(p.diagnostics.to_dicts())
            raise HTTPException(status_code=400, detail={
                "message": "protocol rejected", "ok": False,
                "diagnostics": p.diagnostics.to_dicts(),
                "errors": errors, "warnings": warnings,
                "refusal": editor_refusal(p)})

        # Plan-before-run: the operator confirmed a specific number, so a
        # program that now expands to a different one must NOT quietly launch.
        confirm = body.get("confirm_run_count")
        if confirm is not None and int(confirm) != p.n_runs:
            raise HTTPException(status_code=409, detail={
                "message": (f"plan/confirm mismatch: you confirmed "
                            f"{int(confirm)} runs, this program now expands to "
                            f"{p.n_runs}. Re-plan before running."),
                "confirmed": int(confirm), "run_count": p.n_runs})

        spec = {"kind": PROTOCOL_JOB_KIND, "source": source}
        for key in ("max_runs", "max_window_s"):
            if body.get(key):
                spec[key] = body[key]
        job = jobs.submit(spec)
        return jsonable({
            "run_id": job.run_id,
            "sweep_id": job.run_id,          # the editor's name for the same id
            "status": job.status,
            "n_runs": p.n_runs,
            "run_count": p.n_runs,
            "max_runs": p.max_runs,
            "allow_runs": p.allow_runs,
            "warnings": [d.to_dict() for d in p.diagnostics.warnings],
            "estimated_cost": p.estimated_cost,
            "runs": [{"index": r.index, "run_id": None, "status": "queued",
                      "run_id_reason": NO_RUN_ID_REASON} for r in p.runs],
        })

    def _protocol_job(run_id: str):
        job = jobs.get(run_id)
        if job is None:
            raise HTTPException(status_code=404,
                                detail=f"no such run {run_id!r}")
        if job.spec.get("kind") != PROTOCOL_JOB_KIND:
            raise HTTPException(status_code=409, detail={
                "message": f"run {run_id!r} is not a protocol run",
                "kind": job.spec.get("kind"),
                "reason": "use GET /api/runs/{id} for a RunSpec run"})
        return job, ((job.body or {}).get("protocol") or {})

    @router.get("/runs/{run_id}")
    def get_protocol_run(run_id: str) -> Dict:
        """Results + the sweep grid. Safe to poll while the job is running."""
        job, body = _protocol_job(run_id)
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
            # No `ranking` key: `rank by` and cross-compound ordering were
            # withdrawn on 2026-08-21 (PROTOCOL_LANGUAGE.md §5). Each row
            # carries its own metrics; nothing orders the rows for the caller.
            "sweep_axes": body.get("sweep_axes", []),
            "expects": body.get("expects", []),
            "diagnostics": body.get("diagnostics", []),
            "reconstructions": body.get("reconstructions"),
            "plan": body.get("plan"),
            "log_tail": job.log_tail(40),
            # The editor's projection of the same job, under its own key names
            # (app/web/js/protocol/api.js). See app/server/protocol_wire.py for
            # why both live in one payload.
            **sweep_snapshot(job, body),
        })

    @router.get("/sweeps/{sweep_id}")
    def get_sweep(sweep_id: str) -> Dict:
        """The editor's alias of ``GET /runs/{id}``. A TPL sweep IS one job."""
        job, body = _protocol_job(sweep_id)
        return jsonable(sweep_snapshot(job, body))

    @router.delete("/sweeps/{sweep_id}")
    def cancel_sweep(sweep_id: str) -> Dict:
        """Cancel a sweep.

        Cooperative, exactly like ``DELETE /api/runs/{id}``: a queued sweep dies
        at once, a running one stops at the next run boundary, and a point
        already inside ``simulate_response`` finishes that solve first — the
        engine is not modified to accept an interrupt.
        """
        job, body = _protocol_job(sweep_id)
        jobs.cancel(sweep_id)
        out = sweep_snapshot(job, (job.body or {}).get("protocol") or body)
        out["cancelled"] = job.status not in ("done", "failed")
        out["reason"] = (
            "cancel flag set. The sweep stops at the next run boundary; a run "
            "already inside simulate_response finishes that solve first.")
        return jsonable(out)

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
        if not raw_name and body.get("_multipart_error"):
            raise HTTPException(status_code=400, detail={
                "message": body["_multipart_error"]})
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

        if body.get("content_bytes") is not None:
            # multipart: already raw bytes, never base64 round-tripped
            data = bytes(body["content_bytes"])
        elif body.get("content_base64") is not None:
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

        updir = uploads_dir()
        own_only(updir, PRIVATE_DIR_MODE)
        target = os.path.join(updir, name)
        tmp = target + ".tmp"
        with open(tmp, "wb") as fh:
            fh.write(data)
        own_only(tmp, PRIVATE_FILE_MODE)
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
            # --- the editor's projection (app/web/js/protocol/api.js) ------ #
            "file": {"name": name, "size": len(data),
                     "format": ext.lstrip(".")},
            "drug": _upload_drug_block(rel, name, preview),
            "warnings": [],
            "caveats": [
                "These are the fields the pipeline can derive from what you "
                "uploaded. A field marked flagged-default is a placeholder, not "
                "a property of your compound.",
                "Potency is NOT derivable from structure. Unless you supply an "
                "ic50 or the structure matches a drug in the GDSC glioma "
                "screen, the IC50 here is a flagged low-confidence default.",
            ],
        })

    return router


# --------------------------------------------------------------------------- #
# Body helpers
# --------------------------------------------------------------------------- #
async def _capped_body(request: Request) -> bytes:
    """Read the request body, refusing anything over :data:`_MAX_BODY_BYTES`
    **before** it is buffered.

    This used to be a bare ``await request.body()``, with the
    ``_MAX_UPLOAD_BYTES`` check only afterwards -- by which point the body had
    already been buffered, split, base64-encoded and decoded, roughly 5.7x the
    request size resident. A ~500 MB POST cost ~2.8 GB and a handful of them
    OOM'd the process, unauthenticated, with no rate limit. A limit that fires
    after the allocation is not a limit.

    Content-Length is checked first because it lets us refuse without reading
    a byte, but it is client-supplied and may be absent or a lie, so the
    stream is counted as it arrives and aborted the moment it crosses.
    """
    declared = request.headers.get("content-length")
    if declared is not None and declared.isdigit() and int(declared) > _MAX_BODY_BYTES:
        raise HTTPException(status_code=413, detail={
            "message": f"request body is {int(declared)} B, above the "
                       f"{_MAX_BODY_BYTES} B limit",
            "reason": "a protocol script or a drug definition is kilobytes"})
    chunks: List[bytes] = []
    total = 0
    async for chunk in request.stream():
        total += len(chunk)
        if total > _MAX_BODY_BYTES:
            raise HTTPException(status_code=413, detail={
                "message": f"request body exceeds the {_MAX_BODY_BYTES} B limit",
                "reason": "a protocol script or a drug definition is kilobytes"})
        chunks.append(chunk)
    raw = b"".join(chunks)
    # Populate starlette's own body cache so anything downstream that calls
    # request.body()/stream() again sees the same bytes rather than an empty
    # or already-consumed stream -- request.body() would have set this.
    request._body = raw
    return raw


async def _body(request: Request) -> Dict:
    """Parse the request body as JSON, multipart, or raw text.

    Multipart is parsed HERE rather than with ``request.form()`` because that
    needs ``python-multipart``, which this checkout does not have and which the
    contract says must not be added ("it must run from a clean checkout with
    only the existing venv plus fastapi"). One file part is all these routes
    ever receive, so the ~20 lines below are cheaper than a dependency.

    SECURITY: the body is capped BEFORE it is buffered -- see
    :func:`_capped_body`.
    """
    raw = await _capped_body(request)
    if not raw:
        return {}
    ctype = request.headers.get("content-type") or ""
    if ctype.startswith("multipart/form-data"):
        return _parse_multipart(raw, ctype)
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


def _parse_multipart(raw: bytes, content_type: str) -> Dict:
    """Extract the first file part of a multipart body -> ``{filename, content_bytes}``.

    Returns the payload as raw ``bytes``. It used to be base64-encoded here
    and immediately decoded again by the caller, which added 2.66x to the peak
    allocation for no benefit.
    """
    marker = "boundary="
    if marker not in content_type:
        return {"_multipart_error": "multipart body with no boundary in the "
                                    "Content-Type header"}
    boundary = content_type.split(marker, 1)[1].split(";")[0].strip().strip('"')
    # SECURITY: RFC 2046 caps a boundary at 70 characters from a fixed set.
    # Unvalidated, `boundary=` (empty) made `sep` == b"--", which split the
    # body on every "--" in it -- 4.2 million fragment objects out of an 8 MB
    # body, ~35 MB of peak allocation, from one request.
    if not _BOUNDARY_RE.match(boundary):
        return {"_multipart_error": "malformed multipart boundary (RFC 2046: "
                                    "1-70 characters from the boundary charset)"}
    sep = b"--" + boundary.encode()
    for part in raw.split(sep):
        if b"\r\n\r\n" not in part:
            continue
        head, _, payload = part.partition(b"\r\n\r\n")
        headers = head.decode("utf-8", errors="replace")
        if "filename=" not in headers:
            continue
        filename = headers.split("filename=", 1)[1].split("\r\n")[0]
        filename = filename.strip().strip('"').strip("'")
        # Drop the ONE CRLF the encoder puts before the next boundary. Exactly
        # one: a loop would eat a trailing blank line that is part of the file.
        if payload.endswith(b"\r\n"):
            payload = payload[:-2]
        return {"filename": filename, "content_bytes": payload}
    return {"_multipart_error": "multipart body carried no file part "
                                "(expected a part with a filename)"}


def _upload_drug_block(path: str, name: str, preview: Dict) -> Dict:
    """Per-field provenance for an uploaded file, as the editor's `drug` block.

    Resolved through the REAL resolver, so what the upload panel shows is what a
    `drug X from file "<path>"` statement would actually compile to — including
    every flagged default. A file that cannot be compiled on its own (a bare
    parameter set missing diffusivity, say) reports the reason instead of a row
    of invented numbers.
    """
    from tmesim.protocol.ast import DrugDef
    from tmesim.protocol.errors import TplError
    from tmesim.protocol.resolve import resolve_drug

    node = DrugDef(line=1, col=1, end_col=2, name=os.path.splitext(name)[0],
                   source_kind="file", source_value=path)
    try:
        resolved = resolve_drug(node, compile_params=True)
    except TplError as exc:
        return {"name": node.name, "source": path, "fields": {},
                "unresolved": [d.to_dict() for d in exc.diagnostics]}
    except Exception as exc:  # noqa: BLE001
        return {"name": node.name, "source": path, "fields": {},
                "unresolved": [_internal(
                    exc, "the backend failed while resolving this drug file")]}
    fields = {
        k: {"value": v["value"], "unit": v["unit"],
            "provenance": v["provenance"], "source": v["source"],
            "note": v["source"]}
        for k, v in resolved.provenance.items() if not k.startswith("_")
    }
    return {"name": resolved.alias, "source": path, "fields": fields,
            "identity": resolved.identity,
            "identity_route": resolved.identity_route,
            "smiles": resolved.smiles, "notes": resolved.notes}


def _diag_only(diags: List[Dict], plan_shape: bool = False) -> Dict:
    """A failure body that still has the SHAPE both clients expect."""
    errors, warnings = editor_diagnostics(diags)
    out: Dict = {"ok": False, "diagnostics": diags, "errors": errors,
                 "warnings": warnings, "n_runs": 0, "run_count": 0,
                 "max_runs": 0,
                 "symbols": {"drugs": [], "protocols": [], "sweeps": []}}
    if plan_shape:
        out.update({"runs": [], "axes": [], "metrics": [], "assumptions": [],
                    "caveats": [], "estimated_cost_s": 0, "cost_note": None,
                    "refusal": None, "expects": []})
    return out


def _internal(exc: Exception,
              context: str = "the backend failed while processing this "
                             "program") -> Dict:
    """An unexpected failure, rendered as a diagnostic instead of a 500.

    Every caller catches ``TplError`` FIRST: a rejected script is a deliberate,
    positioned refusal and reaches the editor verbatim. This is only the
    residue -- an exception the protocol layer did not anticipate -- so the
    class name and message go to stderr under a reference and the editor gets
    the reference. ``redact_exception`` already says "not a problem with what
    you asked for", which is the sentence this diagnostic was carrying by hand.
    """
    return {
        "severity": "error", "line": 0, "col": 0, "end_col": 0,
        "code": "internal",
        "message": redact_exception(exc, context),
        "hint": None,
    }
