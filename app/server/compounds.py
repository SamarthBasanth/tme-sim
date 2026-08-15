"""BRING YOUR OWN COMPOUND -- the Studio's side of ``molecular.byo``.

A researcher's own compound series is the point of the product: they want to rank
THEIR molecules, not the twelve in the shipped library. This module owns the two
things the Studio needs beyond the compile itself:

* **A compile endpoint's payload.** :func:`compile_preview` runs
  ``molecular.byo.compile_user_compound`` with the project's calibrated ``alpha``
  and returns the per-field provenance table plus the can/cannot-conclude lists,
  so the UI never has to decide for itself whether a field is a prediction.

* **A store.** A compiled compound must survive the page reload and the server
  restart, because a run record that names a compound whose structure has been
  forgotten is a run nobody can reproduce. The store keeps the SMILES and the
  supplied IC50 -- the inputs -- and recompiles from them on demand, so a change
  to the physics can never leave a stale DrugParams on disk pretending to be
  current.

Files land in ``app/compounds/<id>/compound.json``, written atomically, exactly
like ``app/server/projects.py`` writes a project. Mutations are serialised by one
lock.

Nothing here decides what a compiled compound MEANS. That vocabulary lives in
``molecular.byo`` (``CAN_CONCLUDE`` / ``CANNOT_CONCLUDE`` / the per-field
``method`` strings) and is passed through verbatim, so the wording on screen and
the wording in the module can never drift.
"""
from __future__ import annotations

import json
import os
import re
import threading
import time
from typing import Dict, List, Optional

from tmesim.molecular.byo import (  # noqa: F401 - CompoundRefused is re-exported
    CompiledCompound,
    CompoundRefused,
    compile_user_compound,
)

from .cache import PRIVATE_FILE_MODE, own_only, private_makedirs, repo_root

__all__ = [
    "CompoundError",
    "CompoundStore",
    "COMPOUND_STORE",
    "compile_preview",
    "compounds_dir",
    "calibrated_alpha",
    "user_compound_names",
    "user_drugparams",
    "user_compound_caveats",
]

#: A compound id is derived from its name and must be a safe path component.
COMPOUND_ID_RE = re.compile(r"^[a-z0-9][a-z0-9_.+-]{0,79}$")

#: Hard ceiling on how many compounds one deployment will hold. A user compound
#: costs a directory and roughly a second of compile time; a picker with
#: thousands of entries is a different product. Refused loudly rather than
#: silently truncated.
MAX_COMPOUNDS = 200


class CompoundError(Exception):
    """A bad request against the compound store (mapped to 400/404/409)."""

    def __init__(self, message: str, status: int = 400,
                 detail: Optional[Dict] = None, code: str = "compound_error"):
        super().__init__(message)
        self.message = message
        self.status = int(status)
        self.detail = detail or {}
        self.code = code


def compounds_dir() -> str:
    """Absolute path to ``app/compounds/`` (created on demand).

    SECURITY: owner-only. This tree holds unpublished SMILES -- with the
    default umask it would land 0644 in a 0755 directory, readable by every
    local account on the machine.
    """
    return private_makedirs(os.path.join(repo_root(), "app", "compounds"))


def compound_id(name: str) -> str:
    """Stable id for a compound name: lowercased, non-safe characters folded."""
    slug = re.sub(r"[^a-z0-9_.+-]+", "-", str(name).strip().lower()).strip("-")
    return slug[:80] or "compound"


# --------------------------------------------------------------------------- #
# alpha
# --------------------------------------------------------------------------- #
_alpha_lock = threading.Lock()
_alpha_cache: Optional[float] = None


def calibrated_alpha(force: bool = False) -> float:
    """The project's single fitted loss constant, memoised.

    ``calibrate_alpha()`` compiles doxorubicin through five Brownian sims, so it
    is not free and it is not per-request. It is the SAME number every library
    drug is compiled with (``app/server/jobs.py`` calls ``calibrate_alpha()``
    once per run), which is the whole reason a user compound must be compiled
    with it too: two compounds calibrated differently cannot be ranked against
    each other.
    """
    global _alpha_cache
    with _alpha_lock:
        if _alpha_cache is not None and not force:
            return _alpha_cache
        from tmesim.validation.calibration import calibrate_alpha

        _alpha_cache = float(calibrate_alpha())
        return _alpha_cache


# --------------------------------------------------------------------------- #
# Reserved names
# --------------------------------------------------------------------------- #
def _curated_names() -> List[str]:
    """Every name a user compound must not shadow: both curated libraries, plus
    mitoxantrone (which ``calibration._compile`` special-cases by name and which
    anchors the out-of-sample revalidation).

    Deliberately derived from the two libraries rather than listed, so a drug
    added to either is reserved the moment it lands.
    """
    from tmesim.molecular.breast_library import BREAST_DRUGS
    from tmesim.molecular.library import GBM_DRUGS

    names = set(GBM_DRUGS) | set(BREAST_DRUGS)
    names.add("mitoxantrone")
    return sorted(names)


# --------------------------------------------------------------------------- #
# Compile
# --------------------------------------------------------------------------- #
def compile_preview(name: str, smiles: str, ic50_uM=None) -> Dict:
    """Compile without storing. Raises :class:`CompoundRefused` on refusal.

    The returned document is exactly what the Studio renders: the compiled
    ``drug_params`` a run would simulate, one provenance row per field, the
    caveats, and the can/cannot-conclude lists.
    """
    try:
        alpha = calibrated_alpha()
        alpha_error = None
    except Exception as exc:  # noqa: BLE001
        alpha, alpha_error = None, f"{type(exc).__name__}: {exc}"

    compiled: CompiledCompound = compile_user_compound(
        name, smiles, ic50_uM=ic50_uM, alpha=alpha,
        reserved_names=tuple(_curated_names()),
    )
    doc = compiled.as_dict()
    doc["id"] = compound_id(compiled.name)
    doc["is_user_compound"] = True
    if alpha_error:
        doc["warnings"] = list(doc.get("warnings", [])) + [
            f"The alpha calibration is unavailable ({alpha_error}), so the loss "
            "rates shown are the raw structure heuristics and are NOT the ones "
            "a run would use. Do not compare these to a library drug."
        ]
    return doc


# --------------------------------------------------------------------------- #
# Store
# --------------------------------------------------------------------------- #
class CompoundStore:
    """User-compiled compounds on disk.

    Only the INPUTS are persisted (name, SMILES, supplied IC50). Every compiled
    number is recomputed on read, so a physics change propagates instead of
    leaving a stale-but-authoritative-looking parameter set in a JSON file.
    """

    def __init__(self, root: Optional[str] = None):
        self.root = root or compounds_dir()
        private_makedirs(self.root)
        self._lock = threading.RLock()
        #: (id, name, smiles, ic50, alpha) -> compiled document. In-process only;
        #: see compile_stored for why this must never reach disk.
        self._compiled: Dict[tuple, Dict] = {}

    # -- paths -------------------------------------------------------------- #
    def _dir(self, cid: str) -> str:
        return os.path.join(self.root, cid)

    def _path(self, cid: str) -> str:
        return os.path.join(self._dir(cid), "compound.json")

    # -- io ----------------------------------------------------------------- #
    def _read(self, cid: str) -> Optional[Dict]:
        try:
            with open(self._path(cid)) as fh:
                doc = json.load(fh)
        except FileNotFoundError:
            return None
        except Exception as exc:  # noqa: BLE001 - a corrupt file is reported
            raise CompoundError(
                f"compound {cid!r} is on disk but unreadable: "
                f"{type(exc).__name__}: {exc}", status=500)
        if not isinstance(doc, dict) or "smiles" not in doc:
            raise CompoundError(
                f"compound {cid!r} is not a compound document", status=500)
        doc.setdefault("id", cid)
        return doc

    def _write(self, doc: Dict) -> Dict:
        d = self._dir(doc["id"])
        private_makedirs(d)
        tmp = os.path.join(d, "compound.json.tmp")
        with open(tmp, "w") as fh:
            json.dump(doc, fh, indent=1, sort_keys=False)
        # chmod the temp file BEFORE the rename, so the SMILES are never
        # world-readable even for the width of the os.replace.
        own_only(tmp, PRIVATE_FILE_MODE)
        os.replace(tmp, os.path.join(d, "compound.json"))
        return doc

    def _ids(self) -> List[str]:
        try:
            names = os.listdir(self.root)
        except OSError:
            return []
        return sorted(n for n in names if COMPOUND_ID_RE.match(n)
                      and os.path.exists(self._path(n)))

    # -- api ---------------------------------------------------------------- #
    def list(self) -> List[Dict]:
        """Stored input records, newest first. Cheap: no compile."""
        with self._lock:
            out = []
            for cid in self._ids():
                doc = self._read(cid)
                if doc:
                    out.append(doc)
        out.sort(key=lambda d: d.get("created_at", 0.0), reverse=True)
        return out

    def names(self) -> List[str]:
        return sorted({str(d["name"]).strip().lower() for d in self.list()})

    def get(self, name_or_id: str) -> Optional[Dict]:
        cid = compound_id(name_or_id)
        with self._lock:
            return self._read(cid)

    def create(self, name: str, smiles: str, ic50_uM=None,
               note: Optional[str] = None) -> Dict:
        """Compile (validating everything) and store the INPUTS.

        Returns the full compiled document, not just the stored record: the
        caller has just paid for the compile and the UI needs it.
        """
        with self._lock:
            existing = self._ids()
            if len(existing) >= MAX_COMPOUNDS:
                raise CompoundError(
                    f"this deployment already holds {len(existing)} user "
                    f"compounds (limit {MAX_COMPOUNDS}). Delete one first.",
                    status=409, code="too_many_compounds")
            existing_names = []
            for cid in existing:
                rec = self._read(cid)
                if rec:                       # a directory can vanish under us
                    existing_names.append(str(rec.get("name", "")))
            reserved = tuple(_curated_names()) + tuple(existing_names)
            alpha = None
            alpha_error = None
            try:
                alpha = calibrated_alpha()
            except Exception as exc:  # noqa: BLE001
                alpha_error = f"{type(exc).__name__}: {exc}"

            compiled = compile_user_compound(
                name, smiles, ic50_uM=ic50_uM, alpha=alpha,
                reserved_names=reserved)

            cid = compound_id(compiled.name)
            record = {
                "id": cid,
                "name": compiled.name,
                "smiles": compiled.smiles,
                "ic50_uM": float(ic50_uM) if compiled.ic50_supplied else None,
                "note": (str(note).strip()[:500] if note else None),
                "created_at": time.time(),
                "origin": "user_smiles",
            }
            self._write(record)

        doc = compiled.as_dict()
        doc.update({"id": cid, "is_user_compound": True,
                    "created_at": record["created_at"],
                    "note": record["note"]})
        if alpha_error:
            doc["warnings"] = list(doc.get("warnings", [])) + [
                f"The alpha calibration was unavailable at compile time "
                f"({alpha_error}); the loss rates shown are raw structure "
                "heuristics, not the ones a run would use."]
        return doc

    def delete(self, name_or_id: str) -> bool:
        cid = compound_id(name_or_id)
        with self._lock:
            path = self._path(cid)
            if not os.path.exists(path):
                return False
            os.remove(path)
            try:
                os.rmdir(self._dir(cid))
            except OSError:
                pass
            for key in [k for k in self._compiled if k[0] == cid]:
                del self._compiled[key]
            return True

    def compile_stored(self, name_or_id: str) -> Dict:
        """Recompile a stored compound from its inputs. 404 -> CompoundError.

        Memoised IN PROCESS ONLY, keyed on the stored inputs. A full compile runs
        five Brownian sims (~0.5 s), and ``GET /api/drugs`` compiles every stored
        compound on every call, so without this a picker with fifty compounds
        would take half a minute to open. The memo is deliberately NOT written to
        disk: a compiled number on disk goes stale the moment the physics moves
        and then looks authoritative and is wrong, which is the exact failure
        this whole feature exists to avoid. A process restart -- which is what a
        code change means -- empties it.
        """
        rec = self.get(name_or_id)
        if rec is None:
            raise CompoundError(f"no user compound named {name_or_id!r}",
                                status=404, code="compound_not_found")
        try:
            alpha = calibrated_alpha()
        except Exception:  # noqa: BLE001
            alpha = None

        key = (rec["id"], rec["name"], rec["smiles"], rec.get("ic50_uM"), alpha)
        with self._lock:
            hit = self._compiled.get(key)
        if hit is not None:
            return dict(hit)

        # Reserved-name check is skipped on recompile: the name is already ours.
        compiled = compile_user_compound(
            rec["name"], rec["smiles"], ic50_uM=rec.get("ic50_uM"), alpha=alpha)
        doc = compiled.as_dict()
        doc.update({"id": rec["id"], "is_user_compound": True,
                    "created_at": rec.get("created_at"),
                    "note": rec.get("note")})
        with self._lock:
            if len(self._compiled) > MAX_COMPOUNDS * 2:
                self._compiled.clear()
            self._compiled[key] = doc
        return dict(doc)


#: Process-wide store, mirroring PROJECT_STORE in app/server/projects.py.
COMPOUND_STORE = CompoundStore()


# --------------------------------------------------------------------------- #
# Integration points used by api.py / jobs.py
# --------------------------------------------------------------------------- #
def user_compound_names(store: Optional[CompoundStore] = None) -> List[str]:
    """Lower-cased names a RunSpec may legally put in ``drug``."""
    try:
        return (store or COMPOUND_STORE).names()
    except Exception:  # noqa: BLE001 - a broken store must not break /validate
        return []


def user_drugparams(name: str, alpha: Optional[float] = None,
                    store: Optional[CompoundStore] = None):
    """``(DrugParams, Molecule, document)`` for a stored user compound, or None.

    Called by the run pipeline in place of the library lookup. Returns ``None``
    (rather than raising) when the name is not a user compound, so the caller can
    fall through to the curated libraries.
    """
    st = store or COMPOUND_STORE
    rec = st.get(name)
    if rec is None or str(rec.get("name", "")).strip().lower() != str(name).strip().lower():
        return None
    compiled = compile_user_compound(
        rec["name"], rec["smiles"], ic50_uM=rec.get("ic50_uM"),
        alpha=alpha if alpha is not None else calibrated_alpha())
    doc = compiled.as_dict()
    doc.update({"id": rec["id"], "is_user_compound": True,
                "created_at": rec.get("created_at"), "note": rec.get("note")})
    return compiled.drug_params, compiled.molecule, doc


def user_compound_caveats(doc: Dict) -> List[str]:
    """The caveats a RUN carrying this compound must state, verbatim."""
    return [str(c) for c in doc.get("caveats", [])]
