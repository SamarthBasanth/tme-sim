"""Sections ingested from a directory, alongside the seven built-in ones.

``cache.SECTIONS`` is a hardcoded table: seven specimens whose paths are typed
into the source. That is fine for the specimens this repo ships evidence for and
useless for a new patient sample, which is the whole point of
:mod:`tmesim.data.ingest`. This module is the bridge — it turns an ingestable
directory into a section record of exactly the shape ``cache.SECTIONS`` holds,
so :func:`cache.list_sections`, ``/api/validate`` and ``build_reconstruction``
treat it identically to a built-in one.

Guarantees, in order of importance:

1. **A built-in section can never be shadowed.** :func:`register` refuses an id
   that is already in ``cache.SECTIONS``, and refuses an id already ingested at
   a different path. Nothing an ingested section does can change what a
   published number was measured from.
2. **Registration is REFUSED, not degraded, when the sample is not usable.**
   The directory is fully ingested first (which runs every check in
   ``ingest.validate_record``); a failure returns the reason and registers
   nothing.
3. **The capability report is stored with the section**, so the UI can show what
   the specimen cannot support before a run rather than after.
4. **Persistence is by PATH, not by content.** The store in
   ``app/runs/_ingested_sections.json`` records the directory and the metadata
   measured from it; on reload the path is re-checked and a section whose bytes
   have gone is reported ``available: false`` with a reason, exactly like a
   built-in one with a missing cache.
"""
from __future__ import annotations

import json
import os
import re
import threading
from collections import OrderedDict
from typing import Dict, List, Optional

from .cache import PreformattedFailure
from .cache import redact_exception as _internal

__all__ = [
    "STORE_FILE",
    "ingested_sections",
    "register",
    "unregister",
    "scan",
    "load_store",
    "diagnose",
    "registered_at",
]

#: Sidecar in ``app/runs/`` holding the ingested-section table.
STORE_FILE = "_ingested_sections.json"

_lock = threading.RLock()
_store: "Optional[OrderedDict[str, Dict]]" = None

_ID_RE = re.compile(r"[^a-z0-9_]+")


def _store_path() -> str:
    from .cache import runs_dir

    return os.path.join(runs_dir(), STORE_FILE)


def load_store() -> "OrderedDict[str, Dict]":
    """The persisted ingested-section table, read once per process.

    "Once per process" is load-bearing, not just an optimisation: if a
    SEPARATE process registers a new section (a data-ingestion script, or a
    second `python -m app.server` briefly started to run `/api/sections/scan`)
    while a long-running Studio server is already up, that server's `_store`
    is stale until it restarts -- confirmed live: a Studio instance running
    since before this session's GBM-Space ingestion work only exposed the 7
    built-in sections via `/api/sections`, not the 30 sections a fresh
    process could see immediately from the same on-disk file. There is no
    live-reload here; restart the server to pick up sections registered by
    another process.
    """
    global _store
    with _lock:
        if _store is None:
            try:
                with open(_store_path()) as fh:
                    raw = json.load(fh)
                _store = OrderedDict(
                    (k, v) for k, v in raw.items() if isinstance(v, dict))
            except Exception:  # noqa: BLE001 - a corrupt store is not fatal
                _store = OrderedDict()
            for sid, rec in list(_store.items()):
                _wire_visium(sid, rec)
        return _store


def _save() -> None:
    with _lock:
        path = _store_path()
        try:
            tmp = path + ".tmp"
            with open(tmp, "w") as fh:
                json.dump(_store or {}, fh, indent=1, sort_keys=True)
            os.replace(tmp, path)
        except Exception:  # noqa: BLE001 - persistence is best effort
            pass


def _wire_visium(section_id: str, rec: Dict) -> None:
    """Make an ingested Visium section reachable by ``reconstruct_from_visium``.

    That reconstructor takes a ``sample_id`` and resolves it through
    ``data.real.visium``'s registry, so an ingested Visium directory has to be
    registered there under the id this section uses. Xenium needs nothing: its
    reconstructor already takes a path.
    """
    if rec.get("source") != "visium":
        return
    try:
        from tmesim.data.real import visium

        visium.register_sample(
            rec["sample"], rec["cache_dir"],
            source=f"ingested Space Ranger output at {rec['cache_dir']}",
            accession=rec.get("accession") or "unregistered",
            diagnosis="unknown (ingested by path)",
        )
    except KeyError:
        # Already registered elsewhere; the section stays but will report
        # unavailable through the normal availability check.
        pass
    except Exception:  # noqa: BLE001
        pass


def _make_id(path: str, sample_id: Optional[str]) -> str:
    base = sample_id or os.path.basename(os.path.abspath(path).rstrip(os.sep))
    slug = _ID_RE.sub("_", str(base).strip().lower()).strip("_")
    return slug or "ingested_sample"


def ingested_sections() -> "OrderedDict[str, Dict]":
    """Ingested sections, in ``cache.SECTIONS`` record shape."""
    return OrderedDict(load_store())


def scan(path: str, measure_detection_rates: bool = False) -> Dict:
    """``ingest.scan_report`` — detection plus a capability report, no side effects.

    This is the HTTP boundary, so an UNEXPECTED exception's text is replaced
    here with a stable sentence and a log reference. ``scan_report`` keeps the
    raw text for ``vivum ingest``, where an operator is at the console, and
    flags it as internal; every deliberate refusal (the detection issues, the
    capability report, ``IngestError``) is unflagged and passes through
    verbatim so this surface can never disagree with the CLI.
    """
    from tmesim.data.ingest import scan_report

    out = scan_report(path, measure_detection_rates=measure_detection_rates)
    if out.pop("error_is_internal", False):
        # PreformattedFailure, not RuntimeError: scan_report already reduced
        # the exception to text, so wrapping in a class that is not the
        # failure's own would put a wrong class name in front of the right one
        # in the log. The sentence the browser sees is unchanged.
        out["error"] = _internal(
            PreformattedFailure(out.get("error") or "unknown"),
            "this sample could not be read")
    return out


def register(
    path: str,
    section_id: Optional[str] = None,
    label: Optional[str] = None,
) -> Dict:
    """Ingest ``path`` and add it to the section registry.

    Returns ``{ok, section_id, section, capability, capability_text, error}``.
    ``ok`` is False with a populated ``error`` whenever the directory is not
    ingestable, fails validation, or the id collides — nothing is registered in
    any of those cases.
    """
    from tmesim.data.ingest import IngestError, detect, format_capability_report, ingest

    from .cache import SECTIONS

    abs_path = os.path.abspath(os.path.expanduser(str(path)))
    sid = _make_id(abs_path, section_id)

    if sid in SECTIONS:
        return {
            "ok": False, "section_id": sid, "section": None,
            "capability": None, "capability_text": None,
            "error": (
                f"section id {sid!r} is a BUILT-IN section "
                f"({SECTIONS[sid]['cache_dir']}). Refusing to shadow it — an "
                "ingested section must never change what a published number "
                "was measured from. Pass a different section_id."
            ),
        }
    store = load_store()
    existing = store.get(sid)
    if existing is not None and existing["cache_dir"] != abs_path:
        return {
            "ok": False, "section_id": sid, "section": None,
            "capability": None, "capability_text": None,
            "error": (
                f"section id {sid!r} is already ingested from "
                f"{existing['cache_dir']}. Pass a different section_id, or "
                f"DELETE /api/sections/{sid} first."
            ),
        }

    det = detect(abs_path)
    try:
        record = ingest(abs_path, sample_id=sid, detection=det)
    except IngestError as exc:
        return {"ok": False, "section_id": sid, "section": None,
                "capability": None, "capability_text": None, "error": str(exc)}
    except Exception as exc:  # noqa: BLE001 - report, never 500
        # UNEXPECTED, as opposed to the IngestError branch above, which is a
        # deliberate refusal and is returned verbatim.
        return {"ok": False, "section_id": sid, "section": None,
                "capability": None, "capability_text": None,
                "error": _internal(exc, "this sample could not be ingested")}

    ex, ey = record.extent_um
    rec: Dict = {
        "id": sid,
        "label": label or f"{record.platform.title()} · {sid} (ingested)",
        "source": record.platform if record.platform in ("xenium", "visium")
        else "xenium",
        # Xenium reconstructs from a PATH; Visium from a registered sample id.
        "sample": abs_path if record.platform == "xenium" else sid,
        "cache_dir": abs_path,
        "panel": record.source.get("panel"),
        "assay": record.source.get("assay"),
        "accession": record.source.get("accession"),
        "has_outcome_label": bool(record.source.get("has_outcome_label", False)),
        # --- ingested-only fields ------------------------------------------ #
        "ingested": True,
        "platform": record.platform,
        "layout": det.get("layout"),
        "detect_confidence": det.get("confidence"),
        "detect_issues": list(det.get("issues") or []),
        "n_cells_ingested": record.n_cells,
        "n_genes_ingested": record.n_genes,
        "extent_um_ingested": [ex, ey],
        "genes": [str(g) for g in record.genes.tolist()],
        "capability": record.capability,
        "validation": record.validation,
    }
    if record.platform == "mex":
        rec["detect_issues"] = list(rec["detect_issues"]) + [
            "platform 'mex' is a generic MEX + coordinates bundle. It is "
            "ingestable and its capability report is real, but the "
            "reconstruction pipeline has no route for it: reconstruct_real "
            "accepts 'xenium' or 'visium' only. This section is registered as "
            "INSPECTABLE, not runnable."
        ]
        rec["runnable"] = False
    else:
        rec["runnable"] = True

    with _lock:
        store[sid] = rec
        _save()
    _wire_visium(sid, rec)

    return {
        "ok": True, "section_id": sid, "section": rec,
        "capability": record.capability,
        "capability_text": format_capability_report(record.capability),
        "error": None,
    }


def unregister(section_id: str) -> Dict:
    """Remove an ingested section. Built-in sections are never removable."""
    from .cache import SECTIONS

    if section_id in SECTIONS:
        return {"ok": False, "removed": False,
                "error": f"{section_id!r} is a BUILT-IN section and cannot be "
                         "removed."}
    store = load_store()
    with _lock:
        rec = store.pop(section_id, None)
        if rec is None:
            return {"ok": False, "removed": False,
                    "error": f"no ingested section {section_id!r}"}
        _save()
    if rec.get("source") == "visium":
        try:
            from tmesim.data.real import visium

            visium.unregister_sample(rec["sample"])
        except Exception:  # noqa: BLE001
            pass
    return {"ok": True, "removed": True, "section_id": section_id, "error": None}


def ingested_genes(section_id: str) -> List[str]:
    """Panel gene symbols recorded at ingest time, or ``[]``."""
    rec = load_store().get(section_id)
    return list(rec.get("genes") or []) if rec else []


# --------------------------------------------------------------------------- #
# The diagnosis, for a surface that is BROWSING rather than registering
#
# `scan()` above is the full answer: it detects the layout AND builds the
# capability report, which costs a read of features.tsv.gz (and, with
# `measure_detection_rates`, a whole ingest). That is the right cost once a
# person has decided WHICH directory they mean.
#
# A directory browser asks a cheaper question about every directory it shows:
# is this a recognised bundle, and if not, why not? `ingest.detect` already
# answers exactly that from filenames alone -- the registry has always
# produced this diagnosis, and nothing but `vivum ingest` and the scan route
# ever surfaced it, so a browser would show a bare "not ingestable" where the
# engine had a sentence naming the missing file.
#
# Two things this adds on top of `detect`:
#
# 1. **Every absolute path is removed from the text.** `_detect_visium` writes
#    issues like f"{sp}: no scalefactors_json.json" where `sp` is an absolute
#    path; `_detect_xenium` names a tarball by basename only. Browsing is a
#    read of NAMES (see app/server/browse.py) and a diagnosis that smuggled a
#    resolved absolute path back out would defeat that -- and would be swept
#    up by tests/test_response_redaction.py for good reason. The stripped form
#    is also simply better copy: "spatial: no scalefactors_json.json" is the
#    sentence a person reading a file list wants.
# 2. **It says whether this directory is ALREADY a section.** The store is
#    keyed by path, so this is the one place that can answer it, and it is the
#    difference between a browser that offers to register the same sample
#    twice and one that tells you it is already done.
# --------------------------------------------------------------------------- #

#: Keys of ``detect()`` that are safe to publish verbatim -- scalars and
#: enumerations, no paths. ``files`` and ``path`` are deliberately absent.
_DIAGNOSIS_KEYS = ("platform", "layout", "confidence")


def _strip_paths(text: str, *prefixes: str) -> str:
    """Remove the directory's own absolute path from a diagnosis sentence.

    Longest prefix first, so a realpath and an abspath of the same directory
    (which differ only when a symlink is involved) cannot leave half of one
    behind. ``<prefix>/x`` becomes ``x`` and a bare ``<prefix>`` becomes
    ``this directory``, which is what the sentence meant in both cases.
    """
    out = str(text)
    for prefix in sorted({p for p in prefixes if p}, key=len, reverse=True):
        out = out.replace(prefix + os.sep, "")
        out = out.replace(prefix, "this directory")
    return out


def registered_at(path: str) -> Optional[str]:
    """The id of the ingested section already registered at ``path``, or None.

    Compared on ``realpath`` at BOTH ends: the store records the path the
    operator registered and a browser arrives with the path it resolved, and
    those are the same directory reached two ways whenever a symlink is in
    play. A string comparison would offer to re-register a sample that is
    already in the list.
    """
    try:
        target = os.path.realpath(os.path.expanduser(str(path)))
    except (OSError, ValueError):
        return None
    for sid, rec in load_store().items():
        stored = rec.get("cache_dir")
        if not stored:
            continue
        try:
            same = os.path.realpath(stored) == target
        except (OSError, ValueError):
            continue
        if same:
            return sid
    return None


def diagnose(path: str) -> Dict:
    """Why ``path`` is, or is not, a recognised tissue bundle. Names only.

    Returns ``{platform, layout, confidence, ingestable, issues, notes,
    samples, registered_as, error}``. ``samples`` maps a CHILD DIRECTORY NAME
    to its own detection summary for every child that is itself a recognised
    bundle -- which is what ``detect`` already walks the directory to find, so
    a browser gets it for the same cost rather than probing each child again.

    Never raises: an unreadable directory is an answer. An UNEXPECTED failure
    is redacted (``cache.redact_exception``) because this is an HTTP surface;
    the detection issues are deliberate refusals and pass through verbatim,
    minus the absolute paths described above.
    """
    abs_path = os.path.abspath(os.path.expanduser(str(path)))
    out: Dict = {
        "platform": "unknown", "layout": None, "confidence": 0.0,
        "ingestable": False, "issues": [], "notes": [], "samples": {},
        "registered_as": registered_at(abs_path), "error": None,
    }
    try:
        real_path = os.path.realpath(abs_path)
    except (OSError, ValueError):
        real_path = abs_path
    try:
        from tmesim.data.ingest import detect

        det = detect(abs_path)
    except Exception as exc:  # noqa: BLE001 - report, never 500
        out["error"] = _internal(exc, "this directory could not be inspected")
        return out

    for key in _DIAGNOSIS_KEYS:
        if key in det:
            out[key] = det[key]
    out["ingestable"] = bool(det.get("ingestable"))
    out["issues"] = [_strip_paths(s, real_path, abs_path)
                     for s in (det.get("issues") or [])]
    out["notes"] = [_strip_paths(s, real_path, abs_path)
                    for s in (det.get("notes") or [])]
    for child in det.get("contained_samples") or []:
        name = child.get("name")
        if not name:
            continue
        out["samples"][str(name)] = {
            "platform": child.get("platform"),
            "layout": child.get("layout"),
            "confidence": child.get("confidence"),
            "registered_as": registered_at(os.path.join(abs_path, str(name))),
        }
    return out
