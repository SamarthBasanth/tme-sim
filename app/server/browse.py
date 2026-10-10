"""``GET /api/browse`` — navigate the SERVER's filesystem to find tissue.

WHY THIS EXISTS
---------------
Registering a tissue section means typing an absolute filesystem path into a
text box. That is the whole interaction: ``/ingest`` has one ``<input
type=text>`` and the operator is expected to know, from memory, where their
Xenium export lives on the machine running this server.

It cannot be replaced by a file picker, and that is not a limitation of this
module — it is what the product does. A Xenium or Space Ranger bundle is
multiple gigabytes of matrices and **the server reads it in place**; a section
is registered BY PATH and persisted BY PATH (see
``ingest_registry``'s guarantee 4), so the bytes never move and nothing is ever
transferred. A browser ``<input type="file">`` cannot express "that directory
over there on the host" at all: it would either be a lie about what happens or
a multi-gigabyte transfer of data that is already on the right machine.

So the fix is the other half of the same capability: let the person NAVIGATE
the server's filesystem instead of reciting it. Same input, same confinement,
same registration — no typing, and no false metaphor.

WHAT IT RETURNS, AND WHAT IT DELIBERATELY DOES NOT
--------------------------------------------------
A listing is **names and a type**. Not sizes, not modification times, not
owners, not inode counts — nothing that is a function of a file's CONTENT or
of anything but its presence under a directory the caller already named. A
size column would make this route a read primitive over directories it is only
supposed to let you walk.

Entries carry **no paths**. The caller is handed ``path`` (the one directory it
asked about, resolved) and a list of child NAMES, and composes the next request
itself. One absolute path per response, and it is the caller's own.

CONFINEMENT — this is the sharp edge
------------------------------------
A route that lists directories from a client-supplied path is a traversal and a
disclosure primitive if it is built carelessly, and this codebase has already
had exactly that: ``projects.py::_dir`` built a path from a request-supplied id
and ``%2e%2e`` turned it into a read-and-delete primitive.

**The confinement is ``TMESIM_STUDIO_DATA_ROOTS``, through
:func:`cache.confine_scan_path`, and there is no second mechanism here.** Every
directory this module opens has been through that one function. Read
``cache.py``'s block comment above :data:`cache.DATA_ROOTS_VAR` for why unset
means unconfined (bring-your-own-tissue points outside this repo, and
``__main__.py`` auto-confines to ``<repo>/data_cache`` the moment the bind is
not loopback, so a reachable instance IS confined) and why the check is on
``realpath`` rather than on the text.

Four decisions specific to this module, each one a thing that could reasonably
have gone the other way:

1. **The allowed ROOTS are disclosed; the 403 refusal still does not disclose
   them.** ``confine_scan_path``'s refusal deliberately never names the
   allowlist. A browser cannot start without it: you cannot navigate to a place
   you are not told exists, and a browser that opens on an empty box has not
   replaced the text field at all. So ``GET /api/browse`` with no ``path``
   returns the roots as navigation anchors — to a caller that already holds the
   API key, on a route whose only purpose is navigation. Everything BELOW a
   root is disclosure the key already had: ``POST /api/sections/scan`` returns
   ``detection.contained_samples`` and the resolved path for any allowed
   directory. The root paths themselves are the one new bit, they are printed
   on the server's own startup banner, and the refusal path is untouched —
   pinned by a test.
2. **Symlinks are never followed to compute an entry's type.** ``is_dir()``
   follows by default, which would ``stat`` the target — and for a link
   pointing out of the allowlist that is a stat of something outside the roots,
   i.e. precisely the disclosure the roots exist to prevent (does it exist? is
   it a directory?). A symlink is reported as ``kind: "link"`` and nothing
   else. Clicking it issues an ordinary request, which goes through
   ``confine_scan_path`` like everything else and is refused if it leaves.
3. **The diagnosis is carried, not re-derived.** ``ingest_registry.diagnose``
   runs the registry's own ``detect`` and strips the absolute paths out of its
   sentences. A browser that said "not ingestable" while the engine had
   "spatial: no scalefactors_json.json" would be hiding the answer.
4. **Nothing here registers anything.** Registration stays ``POST
   /api/sections/scan`` with ``register: true`` — the route that already owns
   the ingest, the refusals and the confinement. This module adds a way to
   FIND a directory and not one new thing you can do to it.

MOUNT
-----
One line in ``app/server/api.py``, anywhere among the other ``include_router``
calls (before ``_mount_web(app)``, like all of them)::

    app.include_router(build_browse_router())
"""
from __future__ import annotations

import os
import sys
from typing import Dict, List, Optional

from fastapi import APIRouter, HTTPException, Request
from fastapi.responses import FileResponse

from .cache import (DATA_ROOTS_VAR, PathNotAllowed, confine_scan_path,
                    data_roots, repo_root)

__all__ = ["build_browse_router", "MAX_ENTRIES", "BROWSE_PAGE"]

#: Hard cap on entries in one listing. A directory of 200k files must not turn
#: a navigation click into a 20 MB response, and a route anyone with the key
#: can hit in a loop must have a bounded answer. Over the cap the response says
#: ``truncated: true`` and the client narrows with ``q=``; it does not pretend
#: the directory is small.
MAX_ENTRIES = 500

#: The page this router serves at ``/tissue``. Relative to ``app/web/``.
BROWSE_PAGE = "ingest.html"

#: One warning per process. ``cache._warn_unconfined`` already fires for the
#: scan route and is worded for it; this names THIS capability, because an
#: operator reading "scan is unconfined" has not been told that a directory
#: browser is unconfined too.
_UNCONFINED_WARNED = False


def _warn_unconfined_browse() -> None:
    """Say once that browsing is unconfined, and name the control.

    Deliberately separate from ``cache._warn_unconfined``: that one is a true
    statement about ``POST /api/sections/scan``, fires at most once per
    process, and would therefore be the ONLY thing an operator ever saw if
    this route happened to be the first one used. The risk is the same risk;
    the capability named has to be the one that was actually exercised.
    """
    global _UNCONFINED_WARNED
    if _UNCONFINED_WARNED:
        return
    _UNCONFINED_WARNED = True
    print(
        f"[vivum] GET /api/browse is UNCONFINED: it will list the NAMES in any "
        f"directory this process can read, and {DATA_ROOTS_VAR} is unset. On "
        f"loopback that is the documented default (bring-your-own-tissue "
        f"points at directories outside this repo). On anything reachable it "
        f"is read access to this host as this uid: restart with "
        f"{DATA_ROOTS_VAR} set to a {os.pathsep!r}-separated list of absolute "
        f"directories. This warning prints once per process.",
        file=sys.stderr,
    )
    sys.stderr.flush()


def _place(path: str, why: str) -> Dict:
    """A navigation anchor: where a browse session can start."""
    return {"name": os.path.basename(path.rstrip(os.sep)) or path,
            "path": path, "why": why}


def _places(roots: Optional[List[str]]) -> List[Dict]:
    """Where the browser may start, and why each one is offered.

    Confined: the allowlist, and nothing else — there is nowhere else to go.
    Unconfined: four directories a bring-your-own-tissue operator plausibly
    keeps an export in, which is a convenience and not a boundary (unconfined
    means every path is reachable by typing it; see this module's docstring).
    Only directories that exist are offered, de-duplicated on ``realpath`` so
    a repo checked out inside the home directory does not list twice.
    """
    if roots is not None:
        return [_place(r, "allowed by " + DATA_ROOTS_VAR) for r in roots
                if os.path.isdir(r)]
    repo = repo_root()
    candidates = [
        (os.path.join(repo, "data_cache"), "this checkout's data cache"),
        (repo, "this checkout"),
        (os.path.expanduser("~"), "home directory"),
        (os.getcwd(), "the server's working directory"),
    ]
    out: List[Dict] = []
    seen = set()
    for path, why in candidates:
        try:
            real = os.path.realpath(path)
        except (OSError, ValueError):
            continue
        if real in seen or not os.path.isdir(real):
            continue
        seen.add(real)
        out.append(_place(real, why))
    return out


def _crumbs(resolved: str, roots: Optional[List[str]]) -> List[Dict]:
    """The ancestors of ``resolved`` that the caller may navigate to.

    CLIPPED AT THE ROOT when confinement is in force. Offering ``/`` as a
    breadcrumb on a confined instance would be a row of links that all 403 —
    an interface advertising a capability the server will refuse. Every crumb
    returned is a directory that a subsequent request will actually accept,
    which is also why this walks UP from ``resolved`` rather than splitting a
    string: the containing root is found by the same ``_within`` test that
    allowed the path.
    """
    chain: List[str] = []
    cur = resolved
    while True:
        chain.append(cur)
        parent = os.path.dirname(cur)
        if parent == cur:
            break
        cur = parent
    chain.reverse()
    if roots is not None:
        # Keep the longest allowed prefix and everything below it. A root that
        # IS an ancestor of `resolved` is in `chain` by construction, because
        # the path was allowed on the strength of being under one.
        keep = [p for p in chain if any(_contains(root, p) for root in roots)]
        chain = keep or [resolved]
    return [{"name": os.path.basename(p.rstrip(os.sep)) or p, "path": p}
            for p in chain]


def _contains(root: str, target: str) -> bool:
    """Is ``target`` the root itself or something under it?

    THE SAME TEST, from the same place: ``cache._within`` is the implementation
    the confinement check itself uses, and it compares PATH COMPONENTS via
    ``commonpath`` so that ``/data`` does not contain ``/data_other`` — which a
    ``startswith`` gets wrong. Imported rather than re-written because a second
    copy of a boundary test is how the two of them come to disagree. This is
    used ONLY for presentation (which breadcrumbs to offer); every directory
    actually opened went through :func:`cache.confine_scan_path`.
    """
    from .cache import _within

    return _within(target, root)


def _entries(resolved: str, q: str, show_hidden: bool) -> Dict:
    """Names and a type, sorted directories-first, capped at :data:`MAX_ENTRIES`.

    ``kind`` is one of ``dir``, ``file`` or ``link``, and a ``link`` is NOT
    resolved — see decision 2 in this module's docstring. ``os.scandir`` is
    used for the type because it reads it from the directory entry the kernel
    already returned, so a listing of 500 names costs one syscall rather than
    501.
    """
    needle = q.strip().lower()
    dirs: List[Dict] = []
    files: List[Dict] = []
    total = 0
    try:
        with os.scandir(resolved) as it:
            for entry in it:
                name = entry.name
                if not show_hidden and name.startswith("."):
                    continue
                if needle and needle not in name.lower():
                    continue
                total += 1
                if len(dirs) + len(files) >= MAX_ENTRIES:
                    continue
                try:
                    if entry.is_symlink():
                        kind = "link"
                    elif entry.is_dir(follow_symlinks=False):
                        kind = "dir"
                    else:
                        kind = "file"
                except OSError:
                    # The type could not be read. Say so rather than guessing:
                    # a race (the entry went away) and a permission problem
                    # look the same from here, and neither is a file.
                    kind = "unknown"
                (dirs if kind in ("dir", "link") else files).append(
                    {"name": name, "kind": kind})
    except PermissionError:
        raise HTTPException(status_code=403, detail={
            "code": "directory_unreadable",
            "message": "this server cannot read that directory: the user "
                       "account running it has no permission on it. Nothing "
                       "about its contents is known here."})
    except NotADirectoryError:
        raise HTTPException(status_code=404, detail={
            "code": "not_a_directory",
            "message": "that is a file, not a directory. Browsing lists "
                       "directories; point this at the directory holding the "
                       "export."})
    except OSError:
        raise HTTPException(status_code=404, detail={
            "code": "directory_unreadable",
            "message": "that directory could not be listed on this server."})
    dirs.sort(key=lambda e: e["name"].lower())
    files.sort(key=lambda e: e["name"].lower())
    shown = dirs + files
    return {"entries": shown, "n_entries": total,
            "truncated": total > len(shown)}


def build_browse_router() -> APIRouter:
    """The browse router. Mounts with one line; see this module's docstring."""
    router = APIRouter()

    @router.get("/api/browse")
    async def get_browse(request: Request) -> Dict:
        """List the NAMES in one directory on the server, confined to the
        allowlisted roots.

        Query parameters, all optional:

        ``path``
            The directory to list. Absent (or blank) returns the navigation
            anchors and no listing — the allowlisted roots when
            ``TMESIM_STUDIO_DATA_ROOTS`` is set, otherwise a handful of
            plausible starting points. Resolved with ``realpath`` and then
            re-checked against the allowlist, so ``..``, a ``~`` and every
            symlink are followed BEFORE the decision, and the directory that
            is opened is the directory that was allowed.
        ``q``
            Case-insensitive substring filter on the names, applied BEFORE the
            cap, so a big directory stays navigable instead of merely
            truncated.
        ``show_hidden``
            Include dot-entries. Off by default: ``.DS_Store`` and ``.git``
            between a person and their sample is noise, not information.

        Returns ``{confined, places, path, crumbs, entries,
        n_entries, truncated, diagnosis}``. An entry is ``{name, kind}`` and
        nothing else — no size, no timestamp, no path (see this module's
        docstring). ``diagnosis`` is the registry's own ``detect`` verdict for
        the listed directory, with the absolute paths stripped out of its
        sentences, plus the recognised bundles among the children BY NAME.

        A path outside the allowlist is a 403 in the same ``{code, message,
        path}`` shape ``POST /api/sections/scan`` already refuses with, naming
        the caller's own text, the environment variable and the way out — and
        never the roots.
        """
        params = request.query_params
        raw = params.get("path") or ""
        show_hidden = str(params.get("show_hidden", "")).lower() in (
            "1", "true", "yes", "on")
        q = params.get("q") or ""
        roots = data_roots()
        out: Dict = {
            # `confined` plus an EMPTY `places` is the set-but-unusable state,
            # which refuses everything by design: an operator who set the
            # variable asked for confinement, and a typo or an `-e VAR=` must
            # not hand back the unconfined behaviour they were switching off.
            # Published as two fields rather than one so a client can say
            # "this server is configured to allow nothing" instead of showing
            # an empty disk for a misconfiguration.
            "confined": roots is not None,
            "places": _places(roots),
            "path": None, "crumbs": [], "entries": [], "n_entries": 0,
            "truncated": False, "diagnosis": None,
        }
        if not raw.strip():
            if roots is None:
                _warn_unconfined_browse()
            return out

        # NUL FIRST, and explicitly, because the obvious assumption is wrong:
        # `os.path.realpath` does NOT raise for an embedded null byte on
        # CPython 3.10+ (posixpath swallows the ValueError its own lstat
        # raises) and returns the string with the NUL still in it. So
        # `?path=%00` sailed through the confinement check, failed `isdir`, and
        # came back as a 404 carrying a NUL in a JSON string -- a byte that
        # truncates a C string in whatever reads the log next. It is refused
        # here, and the message echoes no part of the input: a byte string that
        # cannot be a path is not something to hand back.
        if "\x00" in raw:
            raise HTTPException(status_code=400, detail={
                "code": "path_not_a_path",
                "message": "that is not a usable filesystem path on this "
                           "server (it contains a byte a path cannot hold).",
            })

        try:
            confined = confine_scan_path(raw)
        except PathNotAllowed as exc:
            # `from None`: the chain would carry the resolved path into a
            # 500-level log line on some ASGI servers, and the resolved path
            # is deliberately not the client's business.
            raise HTTPException(status_code=403,
                                detail=exc.as_detail()) from None
        except ValueError:
            # Belt and braces for the NUL check above: CPython's posixpath
            # swallows the ValueError today, but `os.path.realpath` is
            # documented to raise it for an embedded null byte and other
            # implementations do. If that ever becomes true here, it is still
            # an answer and not a 500 — with the same wording, so the two
            # paths cannot tell a caller different things about one input.
            raise HTTPException(status_code=400, detail={
                "code": "path_not_a_path",
                "message": "that is not a usable filesystem path on this "
                           "server (it contains a byte a path cannot hold).",
            }) from None
        except OSError:
            raise HTTPException(status_code=400, detail={
                "code": "path_not_a_path",
                "message": "that path could not be resolved on this server "
                           "(it may be too long for this filesystem).",
            }) from None

        if confined is not None:
            # Confinement in force: list the path that was CHECKED, so the
            # directory allowed and the directory opened cannot be two
            # different directories.
            resolved = confined
        else:
            _warn_unconfined_browse()
            # Unconfined, so there is no allowlist to be consistent with --
            # but a browser still needs a canonical path for breadcrumbs and
            # for composing children, and resolving it here is what makes
            # `..` collapse the way the person clicking expects.
            try:
                resolved = os.path.realpath(os.path.expanduser(raw))
            except (OSError, ValueError):
                raise HTTPException(status_code=400, detail={
                    "code": "path_not_a_path",
                    "message": "that is not a usable filesystem path on this "
                               "server.",
                }) from None

        if not os.path.isdir(resolved):
            # The caller's own text, never the resolved path: a caller that
            # sent `~/x` must not be told the operator's home directory, which
            # is the rule `cache.PathNotAllowed` already follows.
            raise HTTPException(status_code=404, detail={
                "code": "no_such_directory",
                "message": "no such directory on the machine running this "
                           "server.",
                "path": raw})

        out["path"] = resolved
        out["crumbs"] = _crumbs(resolved, roots)
        out.update(_entries(resolved, q, show_hidden))
        if not out["truncated"]:
            # The registry's `detect` walks every child to find the bundles
            # among them. On a directory too big to list it is also too big to
            # walk, and a browser that stalls for ten seconds on a click into
            # /usr/lib is worse than one that says "narrow this down".
            from . import ingest_registry

            out["diagnosis"] = ingest_registry.diagnose(resolved)
        return out

    @router.get("/tissue", include_in_schema=False)
    def tissue_page():
        """The bring-your-own-tissue page, at a name that says what it is.

        ``/ingest`` is the ORIGINAL page and is untouched by this router: it is
        registered directly on the app before any router is included, and
        starlette matches in registration order, so it keeps winning that path
        and every link and bookmark to it keeps working. This is the new page
        at a new URL. ``app/web/ingest.html`` is also served as a static asset
        at ``/ingest.html`` by the catch-all mount; this route exists so the
        URL a person is given has no file extension in it.
        """
        page = os.path.join(repo_root(), "app", "web", BROWSE_PAGE)
        if not os.path.exists(page):
            # The front end is not built. A 404 is the honest answer; the API
            # above is unaffected and `/ingest` still serves the original
            # self-contained form.
            raise HTTPException(status_code=404, detail={
                "code": "page_not_built",
                "message": f"app/web/{BROWSE_PAGE} is not present in this "
                           "deployment. The API this page drives is at "
                           "GET /api/browse, and the original form is at "
                           "/ingest."})
        return FileResponse(page, media_type="text/html",
                            headers={"Cache-Control": "no-cache"})

    return router
