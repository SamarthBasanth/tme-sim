"""Section registry + the reconstruction cache.

Two responsibilities, both about the EXPENSIVE half of a run:

1. **Which specimens are actually loadable.** :func:`list_sections` enumerates
   the sections this checkout can reconstruct by looking at the bytes in
   ``data_cache/``. A section whose cache is missing is reported with
   ``available: false`` and a ``reason``; it never raises and never disappears
   from the list, so the UI can show *why* something is unusable.

2. **Caching reconstructions by ``(section, dx_um, vascular_mode, oxygen_mode,
   state_mode)``.** A reconstruction takes 1 s (40k-cell brain section at
   dx = 30 um) to ~90 s (851k-cell GSE324268 section at dx = 12 um) and is
   IDENTICAL across every drug run on the same geometry, so caching it is the
   single biggest responsiveness win available. The cache is in-process (the
   objects are large numpy arrays; a disk round-trip would cost most of the
   saving) and LRU-bounded by entry count.

Everything here WRAPS the production pipeline. ``reconstruct_real`` is the only
thing that builds a microenvironment; the one place this module reaches past it
is ``vascular_mode='kde'``, which ``reconstruct_real`` does not expose as a
parameter — see :func:`build_reconstruction` for exactly what is re-run and why.
"""
from __future__ import annotations

import json
import os
import sys
import threading
import time
import traceback
import uuid
from collections import OrderedDict
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np

# Imported unconditionally on purpose: section_catalog is numpy-free and
# app-free precisely so this module can depend on it without dragging the
# engine into the server's import path.
from tmesim.data.section_catalog import remedy_for

__all__ = [
    "SECTIONS",
    "all_sections",
    "lookup_section",
    "VASCULAR_MODES",
    "OXYGEN_MODES",
    "SENSITIVITY_MODES",
    "STATE_MODES",
    "ReconKey",
    "ReconstructionCache",
    "RECON_CACHE",
    "RECON_CACHE_SIZE_DEFAULT",
    "recon_cache_size",
    "build_reconstruction",
    "get_section",
    "list_sections",
    "section_genes",
    "state_route_support",
    "repo_root",
    "log_internal",
    "PreformattedFailure",
    "redact_exception",
    "runs_dir",
    "DATA_ROOTS_VAR",
    "PathNotAllowed",
    "data_roots",
    "confine_scan_path",
]

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #
_HERE = os.path.dirname(os.path.abspath(__file__))


def repo_root() -> str:
    """Absolute path to the repository root (the parent of ``app/``)."""
    return os.path.abspath(os.path.join(_HERE, "..", ".."))


# --------------------------------------------------------------------------- #
# SECURITY: owner-only modes for everything this app writes
#
# The runtime trees hold unpublished SMILES (app/compounds/), patient tissue
# reconstructions and specimen identifiers (app/runs/, app/uploads/), and the
# operator's project state. Written with the default umask they land 0644 in
# 0755 directories -- readable by every local account on the machine. These
# helpers make them owner-only, and repair the mode of trees that predate this.
# --------------------------------------------------------------------------- #
PRIVATE_DIR_MODE = 0o700
PRIVATE_FILE_MODE = 0o600


def own_only(path: str, mode: int) -> None:
    """Best-effort ``chmod``. A store on a filesystem with no POSIX modes (a
    Windows share, some FUSE mounts) must not break storage outright."""
    try:
        os.chmod(path, mode)
    except OSError:
        pass


class PreformattedFailure(RuntimeError):
    """A failure a LOWER layer already reduced to text, carried up to be logged.

    Two boundaries redact something that was never an exception object by the
    time they saw it: ``ingest_registry.scan`` and
    ``protocol_api.execute_protocol_job`` both receive a string plus an
    ``error_is_internal`` flag, because the engine keeps the real text for a
    console or library caller and marks it rather than pre-formatting a
    sentence. They still need a log record and a reference, so they wrap.

    Wrapping in a bare ``RuntimeError`` made the record read
    ``RuntimeError: FileNotFoundError: no matrix.mtx(.gz) in /Users/...`` --
    a class name that is not the failure's class, in front of one that is.
    This type exists to say "the text is already the whole diagnosis", so
    :func:`log_internal` prints it unprefixed.
    """


def _detail(exc: BaseException) -> str:
    """The real failure as one line: the class name, then the message.

    One definition, used both by the stderr record and by the de-duplication
    key behind ``redact_exception(stable=True)``. If those two ever computed
    the detail differently, the key would stop matching the record it is
    supposed to stand for and the same failure would log twice under two ids.
    """
    if isinstance(exc, PreformattedFailure):
        return str(exc)
    return f"{type(exc).__name__}: {exc}"


def log_internal(exc: BaseException, context: str) -> str:
    """Write the real failure to stderr under a short reference; return the
    reference.

    Split out of :func:`redact_exception` so the two are one implementation:
    every redacted sentence in this backend is correlated with exactly one
    stderr record, in exactly one format. A caller that needs to compose its
    own wording still goes through here, so the log side can never drift from
    what the browser was told to quote.
    """
    ref = uuid.uuid4().hex[:8]
    print(f"[vivum:{ref}] {context}: {_detail(exc)}", file=sys.stderr)
    # Only when there IS one. Called from a boundary that was handed a string
    # rather than from an ``except`` block, ``print_exc`` writes the literal
    # line "NoneType: None", which reads like a second, failed diagnosis.
    if sys.exc_info()[0] is not None:
        traceback.print_exc(file=sys.stderr)
    sys.stderr.flush()
    return ref


#: (context, detail) -> reference, for ``redact_exception(stable=True)``.
#: Bounded and FIFO-evicted: this is a de-duplicator, not a record. The record
#: is the stderr line, which is written once per distinct failure.
_STABLE_REFS: "OrderedDict[Tuple[str, str], str]" = OrderedDict()
_STABLE_REFS_MAX = 64
_STABLE_LOCK = threading.RLock()


def redact_exception(exc: BaseException, context: str, *,
                     brief: bool = False, stable: bool = False) -> str:
    """Log the real failure server-side; return a stable, path-free sentence.

    SECURITY / product: three places in this backend answered an *unexpected*
    exception with ``f"{type(exc).__name__}: {exc}"`` straight into an HTTP
    response body. That is not user-facing English, it names an internal
    exception class, and — for anything that is not wrapping the caller's own
    input — the message can carry server-side absolute paths the user never
    supplied. ``list_sections`` is the clear case: it wraps the whole section
    registry, so a missing ``data_cache/`` download or a filesystem error
    surfaces a path from this machine on a route every client hits on load.

    This is deliberately NOT a way to make failures quieter. The full detail,
    including the traceback, goes to stderr under a short reference, and the
    same reference is handed to the caller, so a field tester who hits an
    unexpected error can quote it and the operator can find it exactly.

    It applies ONLY to unexpected exceptions. The product's *deliberate*
    refusals -- ``IngestError``, the capability report, every validation
    sentence -- are informative on purpose, are rendered verbatim so the UI
    can never disagree with the CLI, and must pass through untouched.

    ``brief=True`` returns a one-clause form for the handful of callers that
    embed the redaction INSIDE a larger sentence they own (a warning on a
    compiled compound, a per-row ``reason`` that already names the field that
    went missing). The full three-sentence form nested inside another sentence
    reads as boilerplate and buries the reference; the short form keeps the
    reference, which is the part a field tester has to quote. Both forms log
    identically -- ``brief`` changes what the caller shows, never what the
    operator gets.

    ``stable=True`` is for a surface that RE-RUNS the same failing work rather
    than republishing one result. The script editor re-plans on every edit, so
    a backend fault that survives a keystroke would otherwise write a fresh
    record and hand back a fresh id several times a second: the operator would
    get a log full of near-identical entries and the user would quote whichever
    id happened to be on screen. Keyed on ``(context, detail)``, an identical
    recurrence reuses the first reference and writes nothing further -- the id
    the user quotes still resolves, and it resolves to one record. Recurrence
    is deliberately not re-logged; the entry already says what happened.
    """
    if stable:
        key = (context, _detail(exc))
        with _STABLE_LOCK:
            ref = _STABLE_REFS.get(key)
            if ref is None:
                ref = log_internal(exc, context)
                _STABLE_REFS[key] = ref
                while len(_STABLE_REFS) > _STABLE_REFS_MAX:
                    _STABLE_REFS.popitem(last=False)
        return _sentence(context, ref, brief)
    return _sentence(context, log_internal(exc, context), brief)


def _sentence(context: str, ref: str, brief: bool) -> str:
    """The two wordings, in one place, so ``stable`` cannot grow a third."""
    if brief:
        return (f"{context} — an unexpected internal failure, logged under "
                f"reference {ref}")
    return (f"{context}. This is an unexpected internal failure, not a "
            f"problem with what you asked for. The server logged the details "
            f"under reference {ref} — quote it to whoever runs this server.")


def private_makedirs(path: str) -> str:
    """``os.makedirs`` owner-only. ``mode=`` applies only to directories that
    are actually created, so chmod as well to tighten a pre-existing tree."""
    os.makedirs(path, mode=PRIVATE_DIR_MODE, exist_ok=True)
    own_only(path, PRIVATE_DIR_MODE)
    return path


def state_dir() -> Optional[str]:
    """The operator-chosen persistent root, or ``None``.

    Read per call rather than cached at import, so a test can set it
    and so a value exported after the module loaded still takes.
    """
    raw = os.environ.get("TMESIM_STATE_DIR", "").strip()
    return os.path.abspath(raw) if raw else None


def state_dir_problem() -> Optional[str]:
    """Why ``TMESIM_STATE_DIR`` cannot be used, as one operator sentence, or
    ``None`` if it is unset or usable.

    MEASURED 2026-09-20, in a container, on the configuration this variable
    exists to enable. A platform disk is mounted at a path the image does not
    contain, so the mount point is created ROOT-OWNED; the image runs as uid
    10001; and the first thing the process does with the value is
    ``private_makedirs(<root>/runs)`` while printing its startup banner. The
    result was ``PermissionError: [Errno 13] ... '/var/lib/vivum/runs'`` and a
    dead container -- the whole service down, at boot, because of a setting
    whose only job is to make projects survive a redeploy.

    A traceback is the wrong answer to that. It names ``os.makedirs`` and a
    frozen importlib frame, and the one fact the operator needs -- the running
    uid is not the mount's owner -- appears nowhere in it. This reduces it to
    the sentence, and the caller decides whether to die.

    NOT a permission to continue silently. Falling back to the ephemeral
    default here would leave the operator believing work is being kept while
    every redeploy destroys it, which is strictly worse than not starting: the
    product asks people to organise work into projects, and a workspace that
    forgets betrays that invitation quietly. Refuse loudly, name the fix.
    """
    root = state_dir()
    if root is None:
        return None

    # The check is a real write, not os.access(). os.access answers for the
    # REAL uid against the permission bits and is wrong under exactly the
    # conditions that matter here -- a root-squashing network mount, a
    # read-only filesystem, an ACL the mode bits do not describe. The only
    # trustworthy test of "can this process create its stores here" is to
    # create something and remove it.
    probe = os.path.join(root, ".vivum-write-probe")
    try:
        os.makedirs(root, mode=PRIVATE_DIR_MODE, exist_ok=True)
        with open(probe, "w", encoding="utf-8") as fh:
            fh.write("")
        os.unlink(probe)
    except OSError as exc:
        try:
            st = os.stat(root)
            owner = f"uid {st.st_uid}:gid {st.st_gid}, mode {oct(st.st_mode & 0o777)}"
        except OSError:
            owner = "could not stat it"
        return (
            f"TMESIM_STATE_DIR is set to {root}, and this process cannot "
            f"write there: {type(exc).__name__}: {exc}. "
            f"It runs as uid {os.getuid()}:gid {os.getgid()}; that directory "
            f"is owned by {owner}. "
            "A platform disk mounted at a path the image does not already "
            "contain is created root-owned, which is this failure exactly. "
            "Fix it by mounting the disk where the image pre-creates an "
            "owned directory (this image: /var/lib/vivum), or by chowning the "
            "mount to the running uid. "
            "Refusing to start rather than falling back, because the fallback "
            "is destroyed on every redeploy and you would not be told."
        )
    return None


def runs_dir() -> str:
    """Absolute path to the run store (created on demand).

    A HOSTED INSTANCE KEEPS NOTHING WITHOUT THIS.

    Measured 2026-09-20 in a real container: this path and
    ``projects.projects_dir()`` were both hardcoded off ``repo_root()``, which
    puts them at ``/app/app/runs`` and ``/app/app/projects`` -- siblings under
    ``/app/app``, whose parent also holds ``server/`` and ``web/``. Render
    allows ONE persistent disk per service, and mounting it at the only
    directory that contains both would hide the application. So every run and
    every project a hosted user created was destroyed on the next deploy.

    That is survivable for runs, which are re-derivable from a RunSpec. It is
    not survivable for PROJECTS, which are the thing a user is asked to
    organise their work into -- a workspace that forgets is worse than no
    workspace, because it invites the investment and then loses it.

    ``TMESIM_STATE_DIR`` is one directory an operator can point a single disk
    at; ``runs/`` and ``projects/`` become children of it, so one mount covers
    both and neither shadows the code. Unset, the behaviour is byte-identical
    to before -- the local single-user workflow this was built for is
    untouched, and no existing installation moves.
    """
    root = state_dir()
    if root:
        return private_makedirs(os.path.join(root, "runs"))
    return private_makedirs(os.path.join(repo_root(), "app", "runs"))


def _abs(*parts: str) -> str:
    return os.path.join(repo_root(), *parts)


# --------------------------------------------------------------------------- #
# SECURITY: allowlisted roots for POST /api/sections/scan
#
# That route takes a filesystem path in the request body and inspects whatever
# directory it names. The API key AUTHENTICATES it and does not CONFINE it, so
# with no allowlist the key is equivalent to read access as this process's uid
# -- which docs/DEPLOYMENT.md §2 has stated honestly ("treat the key as shell
# read access") since it was deferred as a MEDIUM. That was an acceptable trade
# for one operator on their own laptop; it is not acceptable once the key is
# handed to an external field tester.
#
# TMESIM_STUDIO_DATA_ROOTS is the opt-in confinement. Two decisions worth
# stating, because both could reasonably have gone the other way:
#
# 1. UNSET MEANS UNCONFINED, not "confined to data_cache/". Bring-your-own
#    tissue is the product's main input and the documented workflow points the
#    scanner at wherever the operator's export actually lives -- §0 tells a
#    field tester to copy their sample to a study-code directory, which is
#    almost never inside this repo. A secure-by-default root would break every
#    existing local workflow silently, and this repo's precedent for exactly
#    that trade is the unauthenticated loopback default in ``__main__.py``:
#    keep the single-user path working, and make the risk LOUD instead. So the
#    first unconfined scan in a process prints a warning naming the variable
#    (see :func:`_warn_unconfined`).
#
# 2. RESOLVE, THEN CHECK. The check is on ``realpath`` of the request's path,
#    never on the text. A prefix test over the raw string is defeated by ``..``
#    and by a symlink; so, more subtly, is a test over ``abspath`` alone, which
#    collapses ``..`` TEXTUALLY: for ``/allowed/link/../x`` with
#    ``link -> /elsewhere/sub``, ``abspath`` says ``/allowed/x`` -- inside the
#    root -- while the directory that path really designates is
#    ``/elsewhere/x``. (This was not hypothetical: the first version of this
#    module resolved ``realpath(abspath(...))`` and had exactly that hole;
#    tests/test_scan_confinement.py caught it.) That is also why
#    :func:`confine_scan_path` returns the resolved path for the caller to use
#    -- when confinement is in force, the path that was checked is the path
#    that is scanned.
# --------------------------------------------------------------------------- #
#: Environment variable holding ``os.pathsep``-separated allowlisted roots.
DATA_ROOTS_VAR = "TMESIM_STUDIO_DATA_ROOTS"

#: Sentinel for "read the environment", distinct from an explicit ``None``
#: (which means the variable is unset, i.e. unconfined).
_ENV = object()

#: One warning per process, not per scan. Reset in tests.
_UNCONFINED_WARNED = False


class PathNotAllowed(Exception):
    """A caller-supplied path resolved OUTSIDE every allowlisted root.

    Carries a client-safe diagnostic (``message``) and the caller's own path
    verbatim (``path``). The RESOLVED path is deliberately not part of it: a
    caller who sent ``~/x`` would otherwise be told the operator's home
    directory, and one who sent a symlink would be told where it points. The
    resolved path goes to the server's stderr instead, which is this backend's
    standing rule for detail that belongs to the operator and not to the
    client (see :func:`redact_exception`).
    """

    def __init__(self, message: str, *, path: str, resolved: str) -> None:
        super().__init__(message)
        self.message = message
        self.path = path
        self.resolved = resolved

    def as_detail(self) -> Dict:
        """The ``detail`` body for a 403, in this app's ``{code, message}`` shape."""
        return {"code": "path_outside_data_roots",
                "message": self.message,
                "path": self.path}


def data_roots(raw=_ENV) -> Optional[List[str]]:
    """The allowlisted scan roots, fully resolved.

    Returns ``None`` when :data:`DATA_ROOTS_VAR` is unset — the documented
    default, meaning no confinement. Returns a list otherwise, which may be
    EMPTY if the value names nothing usable; an empty list confines scanning to
    nothing, so every scan is refused. That is deliberate fail-closed
    behaviour: an operator who set this variable asked for confinement, and a
    typo or an ``-e VAR=`` with no value must not silently hand back the
    unconfined behaviour they were trying to switch off.

    Entries are separated by ``os.pathsep`` (``:`` on POSIX), NOT by commas
    like ``TMESIM_STUDIO_ALLOWED_HOSTS`` and ``_ALLOWED_ORIGINS`` — a directory
    name may legitimately contain a comma; a hostname or an origin may not.
    """
    if raw is _ENV:
        raw = os.environ.get(DATA_ROOTS_VAR)
    if raw is None:
        return None
    roots: List[str] = []
    for part in str(raw).split(os.pathsep):
        part = part.strip()
        if not part:
            continue
        root = os.path.realpath(os.path.expanduser(part))
        if root not in roots:
            roots.append(root)
    return roots


def _within(target: str, root: str) -> bool:
    """Is ``target`` the root itself or something under it?

    ``os.path.commonpath`` compares PATH COMPONENTS, so ``/data`` does not
    contain ``/data_other`` — which a ``startswith`` test would get wrong.
    """
    try:
        return os.path.commonpath([target, root]) == root
    except ValueError:
        # Different drives, or a mix of absolute and relative. Both arguments
        # are absolute here, so this cannot happen on POSIX; refuse anyway.
        return False


def _warn_unconfined(resolved: str) -> None:
    """Say once per process that scanning is unconfined, and what it means.

    The unconfined default is a deliberate choice (see the block comment
    above) but a silent one is not defensible for a capability that reads any
    directory this process can read. This fires at the moment the capability
    is actually exercised and names the path, which a startup banner line
    cannot do.
    """
    global _UNCONFINED_WARNED
    if _UNCONFINED_WARNED:
        return
    _UNCONFINED_WARNED = True
    print(
        f"[vivum] POST /api/sections/scan is UNCONFINED: it just resolved and "
        f"read {resolved!r}, and it will read ANY directory this process can "
        f"read. {DATA_ROOTS_VAR} is unset, so anyone holding the API key has "
        f"read access to this host as this uid (docs/DEPLOYMENT.md §2). To "
        f"confine it, restart with {DATA_ROOTS_VAR} set to a "
        f"{os.pathsep!r}-separated list of absolute directories, e.g. "
        f"{DATA_ROOTS_VAR}={os.path.join(repo_root(), 'data_cache')!r}. "
        f"This warning prints once per process.",
        file=sys.stderr,
    )
    sys.stderr.flush()


def confine_scan_path(path: str, *, raw=_ENV) -> Optional[str]:
    """Canonicalise ``path`` and confine it to :data:`DATA_ROOTS_VAR`.

    Returns the RESOLVED path when confinement is in force, and ``None`` when
    it is not (the variable is unset) — in which case the caller must pass the
    operator's original path through unchanged, so the default workflow is
    byte-for-byte what it was before this guard existed.

    Raises :class:`PathNotAllowed` when confinement is in force and the
    resolved path is outside every root.

    Two things this does NOT do, stated so nobody reads more into it:

    * It confines the DIRECTORY that is scanned, not the files inside it. A
      symlink *within* an allowed sample directory still points wherever it
      points, and the ingester follows it.
    * It is a check at one instant. A symlink re-pointed between this call and
      the read is not caught (TOCTOU). Both residuals are far smaller than the
      "any path at all" they replace, and neither is closable here.
    """
    # NOT realpath(abspath(...)): abspath collapses ``..`` TEXTUALLY first, so
    # ``<root>/link/../x`` would arrive here already rewritten to ``<root>/x``
    # and be allowed on the strength of a component that link does not have.
    # realpath takes the string itself, resolves each component in order
    # against the filesystem, and makes it absolute on the way -- so a relative
    # path (docs/DEPLOYMENT.md §1d's own example) still resolves against the
    # server's working directory, exactly as the ingester would.
    resolved = os.path.realpath(os.path.expanduser(str(path)))
    roots = data_roots(raw)
    if roots is None:
        _warn_unconfined(resolved)
        return None
    if any(_within(resolved, root) for root in roots):
        return resolved

    # The operator gets the real detail on the console; the client gets a
    # sentence that names its own input, the control, and the way out.
    print(
        f"[vivum] REFUSED scan of {str(path)!r} (resolves to {resolved!r}): "
        f"outside every root in {DATA_ROOTS_VAR}.",
        file=sys.stderr,
    )
    sys.stderr.flush()
    if not roots:
        raise PathNotAllowed(
            f"this server refuses every scan: {DATA_ROOTS_VAR} is set but "
            f"names no usable directory, so the allowlist of directories "
            f"/api/sections/scan may read is empty. Ask whoever runs this "
            f"server to set {DATA_ROOTS_VAR} to a {os.pathsep!r}-separated "
            f"list of absolute directories and restart it, or to unset it "
            f"entirely to go back to unconfined scanning.",
            path=str(path), resolved=resolved)
    raise PathNotAllowed(
        f"{str(path)!r} is outside the directories this server allows "
        f"/api/sections/scan to read. The path was fully resolved first "
        f"(a leading ~, any .. components and every symlink), so it is the "
        f"directory it really points at that is outside the allowlist, not "
        f"the text you sent. The allowlist is the {DATA_ROOTS_VAR} "
        f"environment variable on the server, and its contents are not "
        f"disclosed here. To scan this sample: move or copy it under a "
        f"directory that is already allowed, or ask whoever runs the server "
        f"to add its root to {DATA_ROOTS_VAR} (a {os.pathsep!r}-separated "
        f"list of absolute directories) and restart.",
        path=str(path), resolved=resolved)


# --------------------------------------------------------------------------- #
# Enumerations mirrored from the engine (imported lazily where they are used so
# importing this module stays cheap; these literals are asserted against the
# engine's own tuples in tests/test_app_server.py).
# --------------------------------------------------------------------------- #
VASCULAR_MODES: Tuple[str, ...] = ("discrete", "kde")
OXYGEN_MODES: Tuple[str, ...] = ("pde", "vessel_distance", "markers")
SENSITIVITY_MODES: Tuple[str, ...] = ("uniform", "determinant")
#: How ``cells.state`` is labelled. Mirrors
#: ``reconstruct.cell_state.STATE_MODES``; asserted equal in
#: ``tests/test_app_server.py::test_engine_mode_tuples_match_the_registry``.
STATE_MODES: Tuple[str, ...] = ("oxygen", "vessel_distance", "markers")


# --------------------------------------------------------------------------- #
# Section registry
# --------------------------------------------------------------------------- #
#: Every specimen this backend knows how to reconstruct. ``source`` and
#: ``sample`` are exactly the arguments handed to
#: ``reconstruct.real_pipeline.reconstruct_real``; nothing is re-derived.
SECTIONS: "OrderedDict[str, Dict]" = OrderedDict([
    ("xenium_319", {
        "id": "xenium_319",
        "label": "Xenium GBM · 319-panel",
        "source": "xenium",
        "sample": None,  # loader default -> data_cache/xenium_gbm
        "cache_dir": os.path.join("data_cache", "xenium_gbm"),
        "panel": "hBrain_254g",
        "assay": "Xenium In Situ (single-cell)",
        "accession": "10x Genomics Xenium Human Brain Preview (GBM section)",
        "has_outcome_label": False,
    }),
    ("xenium_io", {
        "id": "xenium_io",
        "label": "Xenium GBM · IO 480-panel",
        "source": "xenium",
        "sample": os.path.join("data_cache", "xenium_gbm_io"),
        "cache_dir": os.path.join("data_cache", "xenium_gbm_io"),
        "panel": "hImmune_v1 + Add-on B",
        "assay": "Xenium In Situ (single-cell)",
        "accession": "10x Genomics FFPE Human Brain Cancer, Immuno-Oncology panel",
        "has_outcome_label": False,
    }),
    ("gse324268", {
        "id": "gse324268",
        "label": "Xenium GBM · GSE324268 tumor 07",
        "source": "xenium",
        "sample": os.path.join("data_cache", "gse324268_xenium"),
        "cache_dir": os.path.join("data_cache", "gse324268_xenium"),
        "panel": None,  # read from the shipped gene_panel.json at enumeration
        "assay": "Xenium In Situ (single-cell)",
        "accession": "GEO GSE324268 (GSM9572843, GBMtumor07)",
        "has_outcome_label": False,
    }),
    ("visium_10x_parent_gbm", {
        "id": "visium_10x_parent_gbm",
        "label": "Visium GBM · 10x parent section",
        "source": "visium",
        "sample": "10x_parent_gbm",
        "cache_dir": os.path.join("data_cache", "visium_gbm"),
        "panel": "whole transcriptome (Visium v1, 55 um spots)",
        "assay": "Visium spatial gene expression (spot-level)",
        "accession": "Parent_Visium_Human_Glioblastoma",
        "has_outcome_label": False,
    }),
    ("visium_gbm1", {
        "id": "visium_gbm1",
        "label": "Visium GBM · GBM1 (primary)",
        "source": "visium",
        "sample": "GBM1",
        "cache_dir": os.path.join("data_cache", "visium_gbm", "GBM1"),
        "panel": "whole transcriptome (Visium v1, 55 um spots)",
        "assay": "Visium spatial gene expression (spot-level)",
        "accession": "GEO GSE194329 / GSM5833533",
        "has_outcome_label": False,
    }),
    ("visium_gbm2", {
        "id": "visium_gbm2",
        "label": "Visium GBM · GBM2 (recurrent)",
        "source": "visium",
        "sample": "GBM2",
        "cache_dir": os.path.join("data_cache", "visium_gbm", "GBM2"),
        "panel": "whole transcriptome (Visium v1, 55 um spots)",
        "assay": "Visium spatial gene expression (spot-level)",
        "accession": "GEO GSE194329 / GSM5833534",
        "has_outcome_label": False,
    }),
    ("visium_gbm3", {
        "id": "visium_gbm3",
        "label": "Visium GBM · GBM3 (primary)",
        "source": "visium",
        "sample": "GBM3",
        "cache_dir": os.path.join("data_cache", "visium_gbm", "GBM3"),
        "panel": "whole transcriptome (Visium v1, 55 um spots)",
        "assay": "Visium spatial gene expression (spot-level)",
        "accession": "GEO GSE194329 / GSM5833535",
        "has_outcome_label": False,
    }),
])

#: Files a Xenium cache directory must contain to be loadable by
#: ``tmesim.data.real.singlecell_spatial.load``.
_XENIUM_REQUIRED = ("Xenium_GBM_cells.csv.gz",)
_XENIUM_MATRIX_FILES = ("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz")


# --------------------------------------------------------------------------- #
# Built-in + ingested sections
#
# :data:`SECTIONS` stays exactly what it was: the seven specimens this repo
# ships evidence for, in their original order. Sections ingested from a
# directory (``app/server/ingest_registry``, backed by ``tmesim.data.ingest``)
# are appended AFTER them and can never replace one — ``register`` refuses a
# built-in id. So every lookup below is "built-in first, ingested second", and
# the built-in behaviour is unchanged whether or not anything has been ingested.
# --------------------------------------------------------------------------- #
def all_sections() -> "OrderedDict[str, Dict]":
    """Built-in sections, then ingested ones. Never raises."""
    out: "OrderedDict[str, Dict]" = OrderedDict(SECTIONS)
    try:
        from .ingest_registry import ingested_sections

        for sid, rec in ingested_sections().items():
            out.setdefault(sid, rec)
    except Exception:  # noqa: BLE001 - a broken store must not hide the built-ins
        pass
    return out


def lookup_section(section_id: str) -> Optional[Dict]:
    """The registry record for ``section_id``, built-in or ingested."""
    rec = SECTIONS.get(section_id)
    if rec is not None:
        return rec
    try:
        from .ingest_registry import ingested_sections

        return ingested_sections().get(section_id)
    except Exception:  # noqa: BLE001
        return None

#: On-disk memo of the numbers that need a file read (cell counts, extents,
#: panel gene lists). Keyed by section id + a fingerprint of the source files,
#: so it self-invalidates when the cache directory changes.
_META_FILE = "_section_meta.json"

# Re-entrant so _save_meta can prime the memo through _load_meta first.
_meta_lock = threading.RLock()
_meta_cache: Optional[Dict] = None


def _load_meta() -> Dict:
    global _meta_cache
    with _meta_lock:
        if _meta_cache is None:
            path = os.path.join(runs_dir(), _META_FILE)
            try:
                with open(path) as fh:
                    _meta_cache = json.load(fh)
            except Exception:  # noqa: BLE001 - a corrupt memo is not fatal
                _meta_cache = {}
        return _meta_cache


def _save_meta(section_id: str, entry: Dict) -> None:
    _load_meta()  # prime the memo so a save never drops other sections
    with _meta_lock:
        meta = _meta_cache if _meta_cache is not None else {}
        meta[section_id] = entry
        globals()["_meta_cache"] = meta
        path = os.path.join(runs_dir(), _META_FILE)
        try:
            tmp = path + ".tmp"
            with open(tmp, "w") as fh:
                json.dump(meta, fh, indent=1, sort_keys=True)
            os.replace(tmp, path)
        except Exception:  # noqa: BLE001 - the memo is an optimisation only
            pass


def _fingerprint(paths: List[str]) -> str:
    """(size, mtime) digest of the files a section's metadata was read from."""
    parts = []
    for p in sorted(paths):
        try:
            st = os.stat(p)
            parts.append(f"{os.path.basename(p)}:{st.st_size}:{int(st.st_mtime)}")
        except OSError:
            parts.append(f"{os.path.basename(p)}:missing")
    return "|".join(parts)


# --------------------------------------------------------------------------- #
# Availability + cheap metadata, per source
# --------------------------------------------------------------------------- #
def _xenium_paths(rec: Dict) -> Dict[str, str]:
    """Resolve a Xenium cache's components.

    Delegates to ``tmesim.data.ingest.resolve_xenium_files``, which matches the
    repo-canonical ``Xenium_GBM_*`` basenames FIRST — so all three built-in
    caches resolve to exactly the paths this function used to hardcode, and an
    ingested stock 10x bundle (``cells.csv.gz`` / ``cell_feature_matrix/`` /
    ``gene_panel.json``) resolves too. The hardcoded names are kept as the
    fallback so a MISSING file still produces the same "missing
    Xenium_GBM_cells.csv.gz" message a user has seen before.
    """
    root = _abs(rec["cache_dir"])
    found: Dict[str, Optional[str]] = {}
    try:
        from tmesim.data.ingest import resolve_xenium_files

        found = resolve_xenium_files(root)
    except Exception:  # noqa: BLE001 - fall back to the historical names
        found = {}
    return {
        "root": root,
        "cells": found.get("cells")
        or os.path.join(root, "Xenium_GBM_cells.csv.gz"),
        "matrix_dir": found.get("matrix_dir")
        or os.path.join(root, "cell_feature_matrix"),
        "tar": found.get("matrix_tar")
        or os.path.join(root, "Xenium_GBM_cell_feature_matrix.tar.gz"),
        "panel": found.get("panel")
        or os.path.join(root, "Xenium_GBM_gene_panel.json"),
    }


def _xenium_available(rec: Dict) -> Tuple[bool, Optional[str]]:
    p = _xenium_paths(rec)
    if not os.path.isdir(p["root"]):
        return False, (
            f"cache directory {rec['cache_dir']} does not exist; see "
            "docs/real_data/singlecell_spatial.md for the download list"
        )
    if not os.path.exists(p["cells"]):
        return False, (
            f"missing {os.path.relpath(p['cells'], repo_root())} "
            "(the per-cell centroid table)"
        )
    has_mex = all(
        any(os.path.exists(os.path.join(p["matrix_dir"], f))
            for f in (name, name[:-3] if name.endswith(".gz") else name + ".gz"))
        for name in _XENIUM_MATRIX_FILES
    )
    if not has_mex and not os.path.exists(p["tar"]):
        return False, (
            "missing the cell-feature matrix: neither the extracted MEX triple "
            f"in {rec['cache_dir']}/cell_feature_matrix nor the tarball"
        )
    return True, None


def _xenium_panel_name(panel_json_path: str) -> Optional[str]:
    """Panel identity as shipped in the section's ``gene_panel.json``.

    10x records the add-on panel at ``payload.panel.identity.name`` and the base
    panel it extends at ``payload.panel.type.data.base_panel_identity.name``, so
    the honest label is "base + add-on". Returns ``None`` (and the registry's
    static name is kept) if the file is absent or shaped differently.
    """
    if not os.path.exists(panel_json_path):
        return None
    try:
        with open(panel_json_path) as fh:
            doc = json.load(fh)
        panel = ((doc.get("payload") or {}).get("panel") or {})
        addon = (panel.get("identity") or {}).get("name")
        base = ((((panel.get("type") or {}).get("data") or {})
                 .get("base_panel_identity") or {}).get("name"))
        parts = [p for p in (base, addon) if p]
        return " + ".join(parts) or None
    except Exception:  # noqa: BLE001 - the registry name stays
        return None


def _xenium_metadata(rec: Dict) -> Dict:
    """n_cells, measured extent, panel name and gene list, read from disk.

    The extent is the RAW measured centroid bounding-box span, which is
    independent of ``dx_um`` — unlike ``Domain.extent_um``, which is the span
    rounded up to a whole number of grid cells and therefore moves with the
    mesh.
    """
    import gzip

    import pandas as pd

    p = _xenium_paths(rec)
    cells = pd.read_csv(p["cells"], usecols=["x_centroid", "y_centroid"])
    x = cells["x_centroid"].to_numpy(dtype=np.float64)
    y = cells["y_centroid"].to_numpy(dtype=np.float64)

    genes: List[str] = []
    features = next(
        (f for f in (os.path.join(p["matrix_dir"], n)
                     for n in ("features.tsv.gz", "features.tsv"))
         if os.path.exists(f)), None)
    if features is not None:
        opener = gzip.open if features.endswith(".gz") else open
        with opener(features, "rt") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 3 and parts[2].strip() == "Gene Expression":
                    genes.append(parts[1])

    panel = _xenium_panel_name(p["panel"]) or rec.get("panel")

    return {
        "n_cells": int(x.size),
        "extent_um": [float(x.max() - x.min()), float(y.max() - y.min())],
        "n_genes": len(genes) or None,
        "genes": genes,
        "panel": panel,
        "n_spots": None,
        "measured_from": "cells.csv.gz centroids + features.tsv.gz",
    }


def _visium_available(rec: Dict) -> Tuple[bool, Optional[str]]:
    try:
        from tmesim.data.real import visium
    except Exception as exc:  # noqa: BLE001
        return False, redact_exception(
            exc, "the Visium loader could not be imported")
    try:
        present = set(visium.list_samples(present_only=True))
    except Exception as exc:  # noqa: BLE001
        return False, redact_exception(
            exc, "the Visium cache could not be scanned")
    if rec["sample"] in present:
        return True, None
    return False, (
        f"Visium sample {rec['sample']!r} is not cached under "
        f"{rec['cache_dir']}; see docs/real_data/visium_gbm.md"
    )


def _visium_metadata(rec: Dict) -> Dict:
    """Spot count + measured extent, without loading the count matrix.

    ``n_cells`` is deliberately ``None``: a Visium reconstruction does not have
    measured cells, it PLACES cells from the per-spot deconvolved fractions, so
    the count does not exist until ``reconstruct_from_visium`` has run. It is
    filled in from the reconstruction the first time one is built.
    """
    out = {
        "n_cells": None,
        "n_cells_reason": (
            "Visium is spot-level: cells are placed by deconvolution at "
            "reconstruction time, so no cell count exists in the cached bytes. "
            "Populated after this section is reconstructed once."
        ),
        "extent_um": None,
        "n_genes": None,
        "genes": [],
        "panel": rec.get("panel"),
        "n_spots": None,
        "measured_from": "tissue_positions + scalefactors_json",
    }
    try:
        from tmesim.data.real import visium

        spatial_dir, _matrix_dir = visium._locate(rec["sample"])
        coords_px = visium._read_positions(spatial_dir)
        umpp, _diam = visium._read_umpp(spatial_dir)
        pts = np.asarray(list(coords_px.values()), dtype=np.float64) * float(umpp)
        if pts.size:
            out["n_spots"] = int(pts.shape[0])
            out["extent_um"] = [
                float(pts[:, 0].max() - pts[:, 0].min()),
                float(pts[:, 1].max() - pts[:, 1].min()),
            ]
    except Exception as exc:  # noqa: BLE001 - degrade to nulls, never crash
        out["extent_um_reason"] = redact_exception(
            exc, "the spot positions for this section could not be read")
    return out


def _ingested_available(rec: Dict) -> Tuple[bool, Optional[str]]:
    """Is an ingested section's directory still there and still ingestable?

    Re-checked on every enumeration rather than trusted from the store: the
    store records a PATH, and a path can go away. Cheap — ``detect`` only
    stats filenames.
    """
    path = rec.get("cache_dir")
    if not path or not os.path.isdir(path):
        return False, (
            f"ingested directory {path!r} no longer exists. Re-ingest it with "
            f"POST /api/sections/scan, or DELETE /api/sections/{rec['id']}."
        )
    if not rec.get("runnable", True):
        return False, (
            f"platform {rec.get('platform')!r} is ingestable and inspectable "
            "but not reconstructable: reconstruct_real accepts 'xenium' or "
            "'visium' only. Its capability report is still valid."
        )
    try:
        from tmesim.data.ingest import detect

        det = detect(path)
    except Exception as exc:  # noqa: BLE001
        return False, redact_exception(
            exc, "this ingested directory could not be re-detected")
    if not det.get("ingestable"):
        return False, (
            f"{path} is no longer ingestable: "
            + ("; ".join(det.get("issues") or []) or "layout not recognised")
        )
    return True, None


def _metadata_source_files(rec: Dict) -> List[str]:
    if rec["source"] == "xenium":
        p = _xenium_paths(rec)
        return [p["cells"], os.path.join(p["matrix_dir"], "features.tsv.gz"),
                p["panel"]]
    return [_abs(rec["cache_dir"])]


def _ingested_metadata(rec: Dict) -> Dict:
    """Metadata for an INGESTED section, taken from the ingest itself.

    An ingested section was fully loaded and validated at registration, so its
    cell count, extent and gene list are already MEASURED — re-deriving them
    from disk would repeat the work and, for ``platform='mex'``, would need a
    reader the built-in paths do not have. Nothing here is estimated.
    """
    ex = rec.get("extent_um_ingested")
    genes = list(rec.get("genes") or [])
    return {
        "n_cells": rec.get("n_cells_ingested"),
        "extent_um": list(ex) if ex else None,
        "n_genes": rec.get("n_genes_ingested") or (len(genes) or None),
        "genes": genes,
        "panel": rec.get("panel"),
        "n_spots": (rec.get("n_cells_ingested")
                    if rec.get("platform") == "visium" else None),
        "measured_from": (
            f"tmesim.data.ingest.ingest({rec.get('cache_dir')!r}) at "
            f"registration time"
        ),
    }


def _section_metadata(rec: Dict) -> Dict:
    """Cheap on-disk metadata for one section, memoised in ``app/runs``."""
    if rec.get("ingested"):
        return _ingested_metadata(rec)
    fp = _fingerprint(_metadata_source_files(rec))
    memo = _load_meta().get(rec["id"])
    if memo and memo.get("_fingerprint") == fp:
        return memo

    if rec["source"] == "xenium":
        meta = _xenium_metadata(rec)
    else:
        meta = _visium_metadata(rec)
    meta["_fingerprint"] = fp
    if memo:
        # Preserve anything only a completed reconstruction can know.
        for k in ("n_cells_observed", "grid_observed"):
            if memo.get(k) is not None:
                meta[k] = memo[k]
    _save_meta(rec["id"], meta)
    return meta


def record_observed_reconstruction(section_id: str, report: Dict) -> None:
    """Remember what a completed reconstruction measured for this section.

    Lets ``/api/sections`` report a real Visium cell count (which cannot be
    known from the cached bytes) instead of ``null`` once one has been built.
    Only ever records numbers that were actually produced by
    ``reconstruction_report``.
    """
    rec = lookup_section(section_id)
    if rec is None:
        return
    meta = dict(_section_metadata(rec))
    meta["n_cells_observed"] = int(report.get("n_cells") or 0) or None
    grid = report.get("grid")
    meta["grid_observed"] = list(grid) if grid else None
    _save_meta(section_id, meta)


def section_genes(section_id: str) -> List[str]:
    """Gene symbols on this section's panel (empty list if not knowable cheaply).

    Used by ``/api/validate`` to reject ``oxygen_mode='markers'`` on a panel
    that carries no hypoxia-response transcript BEFORE a 40 s run discovers it.
    """
    rec = lookup_section(section_id)
    if rec is None:
        return []
    try:
        return list(_section_metadata(rec).get("genes") or [])
    except Exception:  # noqa: BLE001
        return []


def state_route_support(section_id: str) -> Dict:
    """Which cell-state axes this section's PANEL can carry, from its gene list.

    Cheap (it reads the memoised ``features.tsv.gz`` gene symbols, never the count
    matrix), so ``POST /api/validate`` can answer before a run is paid for. It
    reports **presence on the panel** only — a gene being on the panel does not
    mean it fires often enough for a per-cell call, which is a detection RATE and
    can only be measured by reading the matrix. That distinction is exactly what
    disqualifies the 319-gene flagship panel: MKI67 and TOP2A ARE on it, and they
    fire in 1.42% of tumour cells, 13x below the published 19-40% Ki-67 range
    (``results/cell_state_impact.md`` §1.1). ``genes_known`` says whether any
    verdict here is evidence at all.

    ``endothelial_present`` is the equivalent gate for ``state_mode='vessel_
    distance'``: without an endothelial marker no cell can be typed ENDOTHELIAL,
    the distance transform has no point to measure from, and every cell would be
    labelled HYPOXIC by construction. It is checked across EVERY registered panel
    map, because which map is used is chosen from the panel that was actually
    loaded (``marker_panels.select_marker_map``).

    Returns ``{genes_known, source, proliferation_present, proliferation_missing,
    hypoxia_present, hypoxia_missing, cycle_specific_present,
    endothelial_present, endothelial_missing}``. Never raises — a malformed
    ``section_id`` (the code editor can POST a list or an object) degrades to
    "nothing is known", so ``POST /api/validate`` cannot 500 on it.
    """
    from tmesim.core import CellType
    from tmesim.reconstruct.cell_state import (
        CELL_CYCLE_SPECIFIC_MARKERS,
        HYPOXIA_MARKERS,
        PROLIFERATION_MARKERS,
    )
    from tmesim.reconstruct.marker_panels import PANEL_REGISTRY

    if not isinstance(section_id, str):
        return {
            "genes_known": False, "source": None,
            "proliferation_present": [], "proliferation_missing": [],
            "cycle_specific_present": [],
            "hypoxia_present": [], "hypoxia_missing": [],
            "endothelial_present": [], "endothelial_missing": [],
            "reason": (
                f"section must be a string id, got {type(section_id).__name__}"),
        }

    rec = lookup_section(section_id) or {}
    genes = set(section_genes(section_id))
    endo_markers = sorted({
        g for _name, mmap in PANEL_REGISTRY
        for g in (mmap.get(CellType.ENDOTHELIAL) or ())
    })
    present = lambda ms: [g for g in ms if g in genes]      # noqa: E731
    missing = lambda ms: [g for g in ms if g not in genes]  # noqa: E731
    return {
        "genes_known": bool(genes),
        "source": rec.get("source"),
        "proliferation_present": present(PROLIFERATION_MARKERS),
        "proliferation_missing": missing(PROLIFERATION_MARKERS),
        "cycle_specific_present": present(CELL_CYCLE_SPECIFIC_MARKERS),
        "hypoxia_present": present(HYPOXIA_MARKERS),
        "hypoxia_missing": missing(HYPOXIA_MARKERS),
        "endothelial_present": present(endo_markers),
        "endothelial_missing": missing(endo_markers),
    }


def get_section(section_id: str) -> Optional[Dict]:
    """The registry record for ``section_id``, or ``None`` if unknown."""
    return lookup_section(section_id)


def describe_section(section_id: str) -> Dict:
    """One ``Section`` object exactly as ``GET /api/sections`` returns it."""
    rec = lookup_section(section_id)
    if rec is None:
        raise KeyError(section_id)
    if rec.get("ingested"):
        available, reason = _ingested_available(rec)
    elif rec["source"] == "xenium":
        available, reason = _xenium_available(rec)
    else:
        available, reason = _visium_available(rec)

    out = {
        "id": rec["id"],
        "label": rec["label"],
        "n_cells": None,
        "extent_um": None,
        "panel": rec.get("panel"),
        "available": bool(available),
        # --- additive provenance (the contract's four fields are above) ---
        "source": rec["source"],
        "sample": rec["sample"],
        "cache_dir": rec["cache_dir"],
        "assay": rec.get("assay"),
        "accession": rec.get("accession"),
        "has_outcome_label": rec.get("has_outcome_label", False),
        "reason": reason,
        # A reason without a remedy is a dead end, and this is the exact screen
        # where a new user hits one: a fresh clone has no ``data_cache/``, so
        # the picker is seven sections that all say "unavailable" and name a
        # doc rather than an action. ``remedy_for`` returns the literal command
        # that fixes THIS section (or the shell block, where the data is not a
        # plain GET), ``[]`` when it does not know, and never raises -- so a
        # section the catalog has not heard of degrades to today's behaviour
        # instead of breaking enumeration. The catalog is asserted against
        # SECTIONS in tests/test_cli.py, so drift fails the suite rather than
        # silently producing a remedy that points at the wrong place.
        "remedy": remedy_for(section_id) if not available else [],
        "n_spots": None,
        "n_genes": None,
    }
    if rec.get("ingested"):
        # Everything an ingested section knows that a built-in one does not.
        # ``capability`` is the whole point: the UI can show what the specimen
        # cannot support BEFORE a run, not after a meaningless number.
        out.update({
            "ingested": True,
            "platform": rec.get("platform"),
            "layout": rec.get("layout"),
            "detect_confidence": rec.get("detect_confidence"),
            "detect_issues": rec.get("detect_issues") or [],
            "runnable": rec.get("runnable", True),
            "capability": rec.get("capability"),
            "validation": rec.get("validation"),
        })
    if not available:
        out["n_cells_reason"] = reason
        out["extent_um_reason"] = reason
        return out

    try:
        meta = _section_metadata(rec)
    except Exception as exc:  # noqa: BLE001 - a metadata read must not 500
        out["available"] = False
        out["reason"] = redact_exception(
            exc, "this section's cache is present but could not be read")
        out["n_cells_reason"] = out["reason"]
        out["extent_um_reason"] = out["reason"]
        return out

    n_cells = meta.get("n_cells")
    if n_cells is None:
        n_cells = meta.get("n_cells_observed")
    out["n_cells"] = n_cells
    if n_cells is None and meta.get("n_cells_reason"):
        out["n_cells_reason"] = meta["n_cells_reason"]
    out["extent_um"] = meta.get("extent_um")
    if out["extent_um"] is None and meta.get("extent_um_reason"):
        out["extent_um_reason"] = meta["extent_um_reason"]
    out["n_spots"] = meta.get("n_spots")
    out["n_genes"] = meta.get("n_genes")
    if meta.get("panel"):
        out["panel"] = meta["panel"]
    return out


def list_sections() -> List[Dict]:
    """Every registered section, available or not. Never raises."""
    out = []
    registry = all_sections()
    for sid in registry:
        try:
            out.append(describe_section(sid))
        except Exception as exc:  # noqa: BLE001
            rec = registry[sid]
            out.append({
                "id": sid, "label": rec["label"], "n_cells": None,
                "extent_um": None, "panel": rec.get("panel"),
                "available": False,
                "reason": redact_exception(
                    exc, f"section {sid!r} could not be enumerated"),
                "source": rec["source"], "sample": rec["sample"],
            })
    return out


# --------------------------------------------------------------------------- #
# Reconstruction
# --------------------------------------------------------------------------- #
ReconKey = Tuple[str, float, str, str, str]

#: The state route a key defaults to when a caller omits it. Kept as a literal
#: rather than imported at module scope so importing this module stays cheap; it
#: is asserted against ``real_pipeline.STATE_MODE_DEFAULT`` in the tests.
_STATE_MODE_DEFAULT = "oxygen"


def make_key(section: str, dx_um: float, vascular_mode: str,
             oxygen_mode: str,
             state_mode: str = _STATE_MODE_DEFAULT) -> ReconKey:
    """The cache key. ``dx_um`` is rounded so 12.0 and 12.000000001 share an
    entry rather than silently rebuilding a 90 s reconstruction.

    ``state_mode`` is part of the key because it CHANGES ``cells.state`` on the
    cached ``Microenvironment`` (``'vessel_distance'`` and ``'markers'`` both
    overwrite it; ``'oxygen'`` leaves it as the reconstructor produced it). Two
    specs that differ only in ``state_mode`` are therefore two different
    reconstructions, and serving one from the other's entry would report the
    wrong band occupancies — which is the largest single term in every absolute
    kill number this project produces (``results/cell_state_impact.md`` §3.2:
    up to 0.3252 absolute). A cached reconstruction under a different
    ``state_mode`` is WRONG to reuse, not merely stale.
    """
    return (str(section), round(float(dx_um), 6),
            str(vascular_mode).strip().lower(), str(oxygen_mode).strip().lower(),
            str(state_mode).strip().lower())


def _max_darcy_speed(micro) -> Optional[float]:
    """Peak |vx| + |vy| (um/s) of the interstitial Darcy velocity field.

    Computed with the production ``physics.pressure.solve_ifp`` — the same call
    the engine makes — so the advection CFL that ``/api/validate`` reports is
    the one ``physics.diffusion.simulate_drug_transport`` will actually apply.
    Returns ``None`` (never a guess) if the solve fails.
    """
    try:
        from tmesim.physics.pressure import solve_ifp

        _ifp, velocity = solve_ifp(micro)
        v = np.asarray(velocity, dtype=np.float64)
        return float(np.max(np.abs(v[0]) + np.abs(v[1]))) if v.size else 0.0
    except Exception:  # noqa: BLE001
        return None


def build_reconstruction(
    section_id: str,
    dx_um: float,
    vascular_mode: str,
    oxygen_mode: str,
    state_mode: str = _STATE_MODE_DEFAULT,
    log: Optional[Callable[[str], None]] = None,
) -> Dict:
    """Reconstruct one specimen through the production pipeline.

    ``reconstruct_real`` does the whole job for ``vascular_mode='discrete'``
    (its own default) at any ``oxygen_mode``, so that path is the unmodified
    production call and its numbers are bit-identical to the CLI.

    ``vascular_mode='kde'`` is the one case ``reconstruct_real`` cannot express:
    it hardcodes ``vascular_density_from_cells(domain, cells)``, i.e. the
    discrete segmented-vessel field, with no parameter to select the historical
    KDE. Rather than reimplement the assembler, this wrapper runs the production
    reconstruction and then swaps in
    ``real_fields.vascular_density_from_cells(..., mode='kde')`` — the same
    public function, the same measured endothelial cells — and re-installs
    oxygen through ``real_pipeline._install_oxygen`` because the ``'pde'`` and
    ``'vessel_distance'`` routes are driven BY the vascular field.
    ``oxygen_mode='markers'`` reads hypoxia off transcripts and does not depend
    on vasculature, so it is left alone. The cost of the swap (a discarded
    discrete field, and one discarded oxygen solve) is logged, not hidden.

    Returns ``{micro, report, build_time_s, max_darcy_speed_um_s, notes}``.
    """
    from tmesim.core import OxygenConfig
    from tmesim.reconstruct import real_pipeline
    from tmesim.reconstruct.real_fields import vascular_density_from_cells

    def _log(msg: str) -> None:
        if log is not None:
            log(msg)

    rec = lookup_section(section_id)
    if rec is None:
        raise KeyError(
            f"unknown section {section_id!r}; known: {sorted(all_sections())}"
        )
    vmode = str(vascular_mode).strip().lower()
    omode = str(oxygen_mode).strip().lower()
    smode = str(state_mode).strip().lower()
    if vmode not in VASCULAR_MODES:
        raise ValueError(
            f"unknown vascular_mode {vascular_mode!r}; expected one of "
            f"{list(VASCULAR_MODES)}"
        )
    if smode not in STATE_MODES:
        raise ValueError(
            f"unknown state_mode {state_mode!r}; expected one of "
            f"{list(STATE_MODES)}"
        )

    notes: List[str] = []
    t0 = time.time()

    sample = rec["sample"]
    # Xenium sample paths in the registry are repo-relative; make them absolute
    # so the server works regardless of the process working directory.
    if rec["source"] == "xenium" and sample is not None:
        sample = _abs(sample)

    _log(f"reconstruct_real(source={rec['source']!r}, sample={sample!r}, "
         f"dx_um={dx_um}, oxygen_mode={omode!r}, state_mode={smode!r})")
    micro = real_pipeline.reconstruct_real(
        source=rec["source"], sample=sample, dx_um=float(dx_um),
        oxygen_mode=omode, state_mode=smode,
    )

    if vmode == "kde":
        if rec["source"] != "xenium":
            notes.append(
                "vascular_mode is not applicable to a Visium section: "
                "reconstruct_from_visium builds vascular_density from the "
                "measured ENDOTHELIAL spot fraction and exposes no mode "
                "switch. The requested 'kde' was NOT applied."
            )
            _log("vascular_mode='kde' ignored (Visium has no mode switch)")
        else:
            _log("re-deriving vascular_density with mode='kde' "
                 "(reconstruct_real hardcodes 'discrete'); the discrete field "
                 "and its oxygen solve are discarded")
            micro.fields.vascular_density = vascular_density_from_cells(
                micro.domain, micro.cells, mode="kde"
            )
            meta = dict(micro.metadata)
            rp = dict(meta.get("real_pipeline", {}))
            rp["vascular_source"] = "KDE of REAL ENDOTHELIAL cell positions"
            if omode in ("pde", "vessel_distance"):
                _log(f"re-solving oxygen ({omode}) on the KDE vascular field")
                rp.update(real_pipeline._install_oxygen(
                    micro, omode, OxygenConfig()
                ))
            else:
                notes.append(
                    "oxygen_mode='markers' reads hypoxia off measured "
                    "transcripts and does not depend on the vascular field, "
                    "so it was not re-solved after the KDE swap."
                )
            meta["real_pipeline"] = rp
            micro.metadata = meta
            notes.append(
                "vascular_density was re-derived as the historical 40 um KDE "
                "of the measured endothelial cells AFTER the production "
                "reconstruction, because reconstruct_real() does not expose a "
                "vascular_mode parameter. The ECM field and the cell "
                "reconstruction are the unmodified production ones."
            )
            if smode != "oxygen":
                notes.append(
                    f"cells.state was labelled by state_mode={smode!r} inside the "
                    "production reconstruction and was NOT re-derived after the "
                    "KDE swap. That is correct for both non-default routes and "
                    "not an omission: 'markers' reads transcripts, and "
                    "'vessel_distance' measures distance to the ENDOTHELIAL CELL "
                    "POSITIONS (physics.oxygen_proxy.oxygen_from_vessel_distance), "
                    "not to fields.vascular_density — so neither depends on which "
                    "vascular field is in force."
                )
    else:
        meta = dict(micro.metadata)
        rp = dict(meta.get("real_pipeline", {}))
        rp.setdefault("vascular_source", "DISCRETE segmented vessels")
        meta["real_pipeline"] = rp
        micro.metadata = meta

    report = real_pipeline.reconstruction_report(micro)
    report["vascular_mode"] = vmode
    report["vascular_source"] = micro.metadata.get(
        "real_pipeline", {}).get("vascular_source")
    report["state_mode"] = micro.metadata.get(
        "real_pipeline", {}).get("state_mode", smode)
    build_time = time.time() - t0
    _log(f"reconstruction done in {build_time:.1f}s: "
         f"{report['n_cells']} cells, grid {report['grid']}")

    speed = _max_darcy_speed(micro)
    try:
        record_observed_reconstruction(section_id, report)
    except Exception:  # noqa: BLE001 - bookkeeping must not fail a run
        pass

    return {
        "micro": micro,
        "report": report,
        "build_time_s": float(build_time),
        "max_darcy_speed_um_s": speed,
        "notes": notes,
        # Snapshot of what a run could mutate on the shared cached object, so
        # every run starts from exactly what reconstruct_real produced.
        # ``state`` is snapshotted too now that ``state_mode`` can install a
        # non-default labelling: it is the largest single term in an absolute
        # kill number, so it must not be able to drift across runs sharing a key.
        "pristine_sensitivity": np.array(micro.cells.sensitivity, copy=True),
        "pristine_state": np.array(micro.cells.state, copy=True),
        "pristine_metadata": dict(micro.metadata),
    }


# --------------------------------------------------------------------------- #
# The cache
# --------------------------------------------------------------------------- #
class ReconstructionCache:
    """Thread-safe LRU of reconstructions keyed by
    ``(section, dx_um, vascular_mode, oxygen_mode, state_mode)``.

    Two jobs that want the same key do NOT both reconstruct: the second blocks
    on a per-key lock and then finds the first one's entry (single-flight). The
    default capacity is small because each entry holds the full cell population
    and several float64 grids — the 851k-cell section at dx = 12 um is on the
    order of a hundred MB.

    A tiny, disk-backed sidecar records each key's measured peak Darcy speed so
    ``/api/validate`` can still check the advection CFL after a restart, before
    anything has been rebuilt. Only measured values are written; a key that has
    never been built has no entry, and validation says so.
    """

    #: Sidecar file inside ``app/runs/`` (see the class docstring).
    CFL_FILE = "_cfl_speeds.json"

    def __init__(self, maxsize: int = 4):
        self.maxsize = int(maxsize)
        self._entries: "OrderedDict[ReconKey, Dict]" = OrderedDict()
        self._lock = threading.RLock()
        self._key_locks: Dict[ReconKey, threading.Lock] = {}
        self.hits = 0
        self.misses = 0
        self.evictions = 0
        self._speeds: Optional[Dict[str, float]] = None

    # -- speed sidecar ----------------------------------------------------- #
    @staticmethod
    def _speed_key(key: ReconKey) -> str:
        return "|".join(str(k) for k in key)

    def _speed_path(self) -> str:
        return os.path.join(runs_dir(), self.CFL_FILE)

    def _load_speeds(self) -> Dict[str, float]:
        if self._speeds is None:
            try:
                with open(self._speed_path()) as fh:
                    self._speeds = json.load(fh)
            except Exception:  # noqa: BLE001
                self._speeds = {}
        return self._speeds

    def _store_speed(self, key: ReconKey, speed: Optional[float]) -> None:
        if speed is None:
            return
        speeds = self._load_speeds()
        speeds[self._speed_key(key)] = float(speed)
        try:
            path = self._speed_path()
            tmp = path + ".tmp"
            with open(tmp, "w") as fh:
                json.dump(speeds, fh, indent=1, sort_keys=True)
            os.replace(tmp, path)
        except Exception:  # noqa: BLE001
            pass

    def known_speed(self, key: ReconKey) -> Optional[float]:
        """Measured peak Darcy speed for this exact key, or ``None``.

        ``state_mode`` is ignored when matching. ``physics.pressure.solve_ifp``
        reads ``fields.vascular_density`` and the domain and never reads
        ``cells.state``, so a speed recorded under one state route is an EXACT
        measurement of the same quantity under another at the same
        ``(section, dx, vascular_mode, oxygen_mode)``. Requiring the state route
        to match too would throw that exact evidence away and fall back to the
        cross-mesh tier, which takes the MAXIMUM speed over every recorded dx and
        would then reject a dx the exact measurement says is fine.

        The same is not true of the other four fields, which all change the
        pressure solve or the mesh it runs on, so they are matched exactly.
        """
        with self._lock:
            entry = self._entries.get(key)
            if entry is not None and entry.get("max_darcy_speed_um_s") is not None:
                return float(entry["max_darcy_speed_um_s"])
            speeds = self._load_speeds()
            exact = speeds.get(self._speed_key(key))
            if exact is not None:
                return float(exact)
            # Any state route at the same section/dx/vascular/oxygen key.
            prefix = "|".join(str(k) for k in tuple(key)[:4]) + "|"
            for k, v in speeds.items():
                if k.startswith(prefix):
                    return float(v)
            # Sidecar entries written before state_mode joined the key have only
            # four parts; they are the same measurement and are still valid.
            legacy = "|".join(str(k) for k in tuple(key)[:4])
            if legacy in speeds:
                return float(speeds[legacy])
            return None

    def speeds_for_section(self, section_id: str,
                           vascular_mode: Optional[str] = None
                           ) -> Dict[str, float]:
        """Every measured peak Darcy speed recorded for ``section_id``, keyed by
        the stringified cache key.

        ``vascular_mode`` narrows it to keys whose vascular field was built the
        same way — the interstitial-pressure solve is driven by
        ``fields.vascular_density``, so 'discrete' and 'kde' produce genuinely
        different velocity fields and their speeds should not be mixed. The
        remaining variation across such keys is ``dx`` (and ``oxygen_mode``,
        which ``physics.pressure.solve_ifp`` does not read at all).

        Used to reason about the CFL for a key that has not been built yet;
        every message built from this is labelled with the key it came from and
        never presented as a measurement of the key being validated.
        """
        with self._lock:
            speeds = dict(self._load_speeds())
        prefix = f"{section_id}|"
        out = {k: v for k, v in speeds.items() if k.startswith(prefix)}
        if vascular_mode is None:
            return out
        vm = str(vascular_mode).strip().lower()
        return {k: v for k, v in out.items() if k.split("|")[2:3] == [vm]}

    # -- cache ------------------------------------------------------------- #
    def peek(self, key: ReconKey) -> Optional[Dict]:
        """Return a cached entry without building it (and without counting a
        miss). Refreshes its LRU position on a hit."""
        with self._lock:
            entry = self._entries.get(key)
            if entry is not None:
                self._entries.move_to_end(key)
            return entry

    def _key_lock(self, key: ReconKey) -> threading.Lock:
        with self._lock:
            lock = self._key_locks.get(key)
            if lock is None:
                lock = threading.Lock()
                self._key_locks[key] = lock
            return lock

    def get_or_build(
        self,
        key: ReconKey,
        log: Optional[Callable[[str], None]] = None,
    ) -> Dict:
        """Return the reconstruction for ``key``, building it if necessary.

        The returned dict is the one held in the cache — callers must treat
        ``micro`` as read-only shared state. ``simulate_response`` already works
        on a private copy, and ``apply_drug_sensitivity`` writes only
        ``cells.sensitivity``, which every run sets before simulating.
        """
        entry = self.peek(key)
        if entry is not None:
            with self._lock:
                self.hits += 1
            if log is not None:
                log(f"reconstruction CACHE HIT for {key} "
                    f"(saved ~{entry['build_time_s']:.1f}s)")
            return entry

        with self._key_lock(key):
            entry = self.peek(key)
            if entry is not None:  # another job built it while we waited
                with self._lock:
                    self.hits += 1
                if log is not None:
                    log(f"reconstruction CACHE HIT (after wait) for {key}")
                return entry

            with self._lock:
                self.misses += 1
            if log is not None:
                log(f"reconstruction CACHE MISS for {key}; building")
            section, dx_um, vmode, omode, smode = key
            built = build_reconstruction(section, dx_um, vmode, omode, smode,
                                         log=log)
            built["key"] = key
            built["created_at"] = time.time()
            with self._lock:
                self._entries[key] = built
                self._entries.move_to_end(key)
                while len(self._entries) > self.maxsize:
                    old_key, _old = self._entries.popitem(last=False)
                    self.evictions += 1
                    if log is not None:
                        log(f"evicted reconstruction {old_key} (LRU, "
                            f"maxsize={self.maxsize})")
            self._store_speed(key, built.get("max_darcy_speed_um_s"))
            return built

    def stats(self) -> Dict:
        with self._lock:
            return {
                "hits": self.hits,
                "misses": self.misses,
                "evictions": self.evictions,
                "size": len(self._entries),
                "maxsize": self.maxsize,
                "keys": [list(k) for k in self._entries],
            }

    def clear(self) -> None:
        with self._lock:
            self._entries.clear()
            self._key_locks.clear()


#: LRU capacity used when ``TMESIM_RECON_CACHE_SIZE`` is unset or unusable.
RECON_CACHE_SIZE_DEFAULT = 4


def recon_cache_size(raw: Optional[str] = None,
                     default: int = RECON_CACHE_SIZE_DEFAULT) -> int:
    """``TMESIM_RECON_CACHE_SIZE`` as a usable LRU capacity.

    This used to be ``int(os.environ.get("TMESIM_RECON_CACHE_SIZE", "4"))``
    inline, which is three separate defects for a variable an operator sets by
    hand:

    * ``TMESIM_RECON_CACHE_SIZE=`` (an ``-e VAR=`` with no value, or a shell
      variable that expanded to nothing) and any non-numeric string killed the
      PROCESS AT IMPORT with a bare ``ValueError: invalid literal for int()
      with base 10: ''`` and a traceback naming this module -- no mention of
      the variable, the value or what would have been acceptable. In Docker
      that is a restart loop whose only symptom is a stack trace from a cache
      module, for a dial that is documented as "not a security control".
    * ``0`` and ``-1`` were accepted silently and are not capacities at all:
      ``get_or_build``'s eviction loop (``while len(self._entries) >
      self.maxsize``) drops the entry it has just inserted, so every run
      rebuilds a reconstruction that costs 1 s to ~90 s. A cache that is
      configured off while reporting ``maxsize=0`` is worse than one that
      refuses the setting.

    So: unset is the documented default and says nothing, a usable value is
    used verbatim, and anything else comes back as a DIAGNOSTIC on stderr
    naming the variable, the offending value and the valid range -- this
    project's "a failure a user caused comes back as a diagnostic, never a
    silent default" applied to a variable that previously did the opposite in
    both directions at once.

    It diagnoses and continues on the documented default rather than exiting,
    deliberately:

    * this module is imported by the CLI, the test suite and every library
      caller, not only by the server, so a hard exit turns a cosmetic
      misconfiguration into an unimportable package;
    * the failure is contained. ``maxsize`` decides how much memory the
      process spends re-using work. It changes no result, no protocol number
      and no access decision (``docs/DEPLOYMENT.md`` §4: "Not a security
      control"), so continuing cannot make anything wrong -- only slower;
    * the operator is told, unmissably and before the server's own
      ``recon cache maxsize=`` banner line, which then shows the value that
      really took effect.

    A variable whose bad value COULD change a result or open access should
    exit instead; this one cannot.
    """
    if raw is None:
        raw = os.environ.get("TMESIM_RECON_CACHE_SIZE")
    if raw is None:
        return default  # unset -> the documented default. Not a user error.

    text = str(raw).strip()
    if text == "":
        problem = "it is empty"
    else:
        try:
            value = int(text)
        except ValueError:
            problem = f"{text!r} is not an integer"
        else:
            if value >= 1:
                return value
            problem = (
                f"{value} is below the minimum of 1; a cache that cannot hold "
                "one entry evicts every reconstruction the moment it is built"
            )

    print(
        f"[vivum] TMESIM_RECON_CACHE_SIZE={raw!r} is not usable: {problem}. "
        f"Valid values are integers >= 1 (each entry holds a full cell "
        f"population and several float64 grids; 1 is the minimum, {default} "
        f"is the default). Continuing with {default}. Set --cache-size or "
        f"TMESIM_RECON_CACHE_SIZE to an integer >= 1 to choose it yourself.",
        file=sys.stderr,
    )
    sys.stderr.flush()
    return default


#: Process-wide cache used by the job worker and by ``/api/validate``.
RECON_CACHE = ReconstructionCache(maxsize=recon_cache_size())
