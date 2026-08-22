"""Resolve a TPL ``drug`` block into a ``core.DrugParams`` with PER-FIELD
provenance, and check a specimen's capability before anything is executed.

Three drug sources, in the spec's order of decreasing evidence:

``from smiles "..."``
    compiled through the molecular layer
    (:func:`tmesim.molecular.compile_drug.molecule_to_drugparams`).
``from library "name"``
    :func:`tmesim.molecular.library.get_molecule`.
``from file "path"``
    ``.smi`` / ``.mol`` / ``.sdf`` structures (parsed with RDKit, then compiled
    exactly like a SMILES) or ``.json`` / ``.yaml`` parameter files (every field
    taken verbatim and labelled ``user-supplied``).

Provenance tags, one per DrugParams field:

``measured``
    a wet-lab number — currently only ``ic50`` via ``ic50 from gdsc``.
``computed``
    derived from structure by the molecular layer (diffusivity, permeability,
    binding, clearance, BBB) or by the calibration/barrier stages.
``user-supplied``
    written in the script or in an uploaded parameter file.
``flagged-default``
    a placeholder the pipeline had no evidence for — the Hill coefficient, and
    the 25 uM IC50 fallback.

The identity problem, and why it is solved rather than papered over
-------------------------------------------------------------------
``ic50 from gdsc`` needs a NAME to look up; ``drug TMZ from smiles "CN1C(=O)..."``
supplies a structure and a local alias. :func:`resolve_identity` recovers the
generic name by canonicalising the SMILES with RDKit and matching it against
:data:`tmesim.molecular.library.GBM_DRUGS` — a STRUCTURAL match, not a name
match. When a structure is supplied and matches nothing curated, the identity is
unknown and the GDSC synonym table is deliberately NOT consulted: a TPL
identifier is a local alias, and letting ``drug TMZ from smiles "<anything>"``
inherit temozolomide's measured potency would be exactly the silent substitution
this package exists to prevent. The route that established the identity is
always recorded, and ``ic50 from gdsc`` refuses when there is none.
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from .errors import ResolutionError, collect

__all__ = [
    "PROVENANCE_TAGS", "ResolvedDrug", "resolve_identity", "resolve_drug",
    "load_drug_file", "read_structure_file", "read_parameter_file",
    "check_capability", "resolve_upload_path", "uploads_dir",
]

PROVENANCE_TAGS = ("measured", "computed", "user-supplied", "flagged-default")

#: File extensions ``from file`` accepts.
STRUCTURE_EXTENSIONS = (".smi", ".smiles", ".mol", ".sdf")
PARAMETER_EXTENSIONS = (".json", ".yaml", ".yml")


# --------------------------------------------------------------------------- #
# Result container
# --------------------------------------------------------------------------- #
@dataclass
class ResolvedDrug:
    """One compiled drug plus everything needed to defend each of its numbers."""

    alias: str = ""                     # the TPL identifier
    identity: Optional[str] = None      # generic name, when recoverable
    identity_route: str = ""            # how the identity was established
    source_kind: str = ""               # smiles | library | file
    source_value: str = ""
    smiles: Optional[str] = None
    drug_params: object = None          # core.DrugParams
    molecule: object = None             # molecular.types.Molecule (or None)
    provenance: Dict[str, Dict] = field(default_factory=dict)
    notes: List[str] = field(default_factory=list)
    ic50_resolution: Optional[Dict] = None
    alpha: Optional[float] = None

    def field_table(self) -> List[Dict]:
        """One row per DrugParams field, for the UI's provenance panel."""
        rows = []
        for name, block in self.provenance.items():
            rows.append({"field": name, **block})
        return rows

    def to_dict(self) -> Dict:
        dp = self.drug_params
        return {
            "alias": self.alias,
            "identity": self.identity,
            "identity_route": self.identity_route,
            "source_kind": self.source_kind,
            "source_value": self.source_value,
            "smiles": self.smiles,
            "alpha": self.alpha,
            "drug_params": None if dp is None else {
                "name": dp.name,
                "ic50_uM": float(dp.ic50),
                "hill": float(dp.hill),
                "diffusivity_um2_s": float(dp.diffusivity),
                "binding_1_per_s": float(dp.binding),
                "clearance_1_per_s": float(dp.clearance),
                "bbb_penetrance": float(dp.bbb_penetrance),
                "molecular_weight": float(dp.molecular_weight),
            },
            "provenance": dict(self.provenance),
            "ic50_resolution": self.ic50_resolution,
            "notes": list(self.notes),
        }


def _prov(value, unit: str, tag: str, source: str) -> Dict:
    assert tag in PROVENANCE_TAGS, tag
    return {"value": value, "unit": unit, "provenance": tag, "source": source}


# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #
def _repo_root() -> str:
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def uploads_dir() -> str:
    """``<repo>/app/uploads`` — where ``POST /api/protocol/upload`` writes."""
    path = os.path.join(_repo_root(), "app", "uploads")
    os.makedirs(path, exist_ok=True)
    return path


def resolve_upload_path(path: str, line: int = 0, col: int = 0) -> str:
    """Resolve a ``from file`` path, refusing anything outside the repository.

    Tried in order: as given (relative to the process cwd), ``<repo>/<path>``,
    ``<repo>/app/<path>``. The last one makes the spec's example
    ``"uploads/compound_47.sdf"`` resolve to what
    :func:`uploads_dir` produced.
    """
    root = _repo_root()
    candidates = [path, os.path.join(root, path), os.path.join(root, "app", path)]
    for cand in candidates:
        full = os.path.abspath(cand)
        if not full.startswith(root + os.sep) and full != root:
            continue
        if os.path.isfile(full):
            return full
    tried = ", ".join(os.path.relpath(os.path.abspath(c), root)
                      for c in candidates)
    raise ResolutionError(
        f"no such file {path!r}. Looked for: {tried} (all relative to the "
        "repository root). Upload one with POST /api/protocol/upload, which "
        "returns the path to use here.",
        line=line, col=col, code="file-not-found")


# --------------------------------------------------------------------------- #
# File readers
# --------------------------------------------------------------------------- #
def read_structure_file(full_path: str, line: int = 0, col: int = 0
                        ) -> Tuple[str, Optional[str], str]:
    """Read a ``.smi`` / ``.mol`` / ``.sdf`` file. Returns ``(smiles, name, how)``.

    RDKit is the parser for ``.mol`` / ``.sdf``; ``.smi`` is a text format
    (``SMILES [name]`` per line) and is read directly, then canonicalised with
    RDKit so the identity match downstream is structural rather than textual.
    """
    ext = os.path.splitext(full_path)[1].lower()
    try:
        from rdkit import Chem
        from rdkit import RDLogger

        RDLogger.DisableLog("rdApp.*")
    except Exception as exc:  # pragma: no cover - RDKit is installed here
        raise ResolutionError(
            f"RDKit is required to read {ext} structure files and is "
            f"unavailable ({exc}).", line=line, col=col, code="rdkit-missing")

    if ext in (".smi", ".smiles"):
        with open(full_path) as fh:
            for raw in fh:
                text = raw.split("#", 1)[0].strip()
                if not text:
                    continue
                parts = text.split()
                smi = parts[0]
                name = parts[1] if len(parts) > 1 else None
                mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    raise ResolutionError(
                        f"{os.path.basename(full_path)}: RDKit could not parse "
                        f"the SMILES {smi!r} on the first record.",
                        line=line, col=col, code="bad-structure")
                return Chem.MolToSmiles(mol), name, f"{ext} first record"
        raise ResolutionError(
            f"{os.path.basename(full_path)}: the file contains no SMILES "
            "record.", line=line, col=col, code="empty-structure")

    if ext == ".mol":
        with open(full_path) as fh:
            block = fh.read()
        mol = Chem.MolFromMolBlock(block)
        if mol is None:
            raise ResolutionError(
                f"{os.path.basename(full_path)}: RDKit could not parse this "
                "MDL molfile.", line=line, col=col, code="bad-structure")
        name = (mol.GetProp("_Name").strip()
                if mol.HasProp("_Name") and mol.GetProp("_Name").strip() else None)
        return Chem.MolToSmiles(mol), name, ".mol (MDL molfile)"

    if ext == ".sdf":
        with open(full_path) as fh:
            text = fh.read()
        records = [r for r in text.split("$$$$") if r.strip()]
        if not records:
            raise ResolutionError(
                f"{os.path.basename(full_path)}: the SDF contains no records.",
                line=line, col=col, code="empty-structure")
        mol = Chem.MolFromMolBlock(records[0])
        if mol is None:
            raise ResolutionError(
                f"{os.path.basename(full_path)}: RDKit could not parse the "
                "first SDF record.", line=line, col=col, code="bad-structure")
        name = (mol.GetProp("_Name").strip()
                if mol.HasProp("_Name") and mol.GetProp("_Name").strip() else None)
        how = f".sdf (first of {len(records)} record(s))"
        if len(records) > 1:
            how += " — the remaining records are NOT read; one TPL drug is one molecule"
        return Chem.MolToSmiles(mol), name, how

    raise ResolutionError(
        f"unsupported structure extension {ext!r}. Supported: "
        f"{', '.join(STRUCTURE_EXTENSIONS + PARAMETER_EXTENSIONS)}.",
        line=line, col=col, code="unsupported-file")


def _parse_flat_yaml(text: str, where: str) -> Dict:
    """A deliberately minimal ``key: value`` reader.

    PyYAML is not a dependency of this checkout, and a drug parameter file is a
    flat mapping of numbers and strings. Anything more structured (a list, a
    nested block, an anchor) is REFUSED with a message saying so, rather than
    half-parsed.
    """
    out: Dict = {}
    for lineno, raw in enumerate(text.replace("\r\n", "\n").split("\n"), 1):
        stripped = raw.split("#", 1)[0].rstrip()
        if not stripped.strip():
            continue
        if stripped[:1] in (" ", "\t", "-"):
            raise ResolutionError(
                f"{where}, line {lineno}: this reader accepts only a FLAT "
                "`key: value` mapping (PyYAML is not installed in this "
                "checkout). Indented blocks, lists and anchors are refused "
                "rather than half-parsed — use a .json file for anything "
                "structured.", code="unsupported-yaml")
        if ":" not in stripped:
            raise ResolutionError(
                f"{where}, line {lineno}: expected `key: value`, got "
                f"{stripped.strip()!r}.", code="bad-yaml")
        key, _, value = stripped.partition(":")
        value = value.strip().strip('"').strip("'")
        key = key.strip()
        try:
            out[key] = float(value) if value.lower() not in ("true", "false") \
                else (value.lower() == "true")
        except ValueError:
            out[key] = value
    return out


def read_parameter_file(full_path: str, line: int = 0, col: int = 0) -> Dict:
    """Read a ``.json`` / ``.yaml`` drug parameter file into a flat dict."""
    ext = os.path.splitext(full_path)[1].lower()
    with open(full_path) as fh:
        text = fh.read()
    if ext == ".json":
        try:
            doc = json.loads(text)
        except Exception as exc:
            raise ResolutionError(
                f"{os.path.basename(full_path)}: invalid JSON ({exc}).",
                line=line, col=col, code="bad-json") from None
        if not isinstance(doc, dict):
            raise ResolutionError(
                f"{os.path.basename(full_path)}: expected a JSON object "
                "mapping parameter names to values.",
                line=line, col=col, code="bad-json")
        return doc
    return _parse_flat_yaml(text, os.path.basename(full_path))


def load_drug_file(path: str, line: int = 0, col: int = 0) -> Dict:
    """Read an uploaded drug file and return a preview dict.

    ``{"path", "kind", "smiles", "name", "how", "parameters", "descriptors"}``
    — ``kind`` is ``"structure"`` or ``"parameters"``.
    """
    full = resolve_upload_path(path, line, col)
    ext = os.path.splitext(full)[1].lower()
    rel = os.path.relpath(full, _repo_root())
    if ext in STRUCTURE_EXTENSIONS:
        smiles, name, how = read_structure_file(full, line, col)
        from ..molecular import descriptors as _desc

        mol = _desc.from_smiles(name or os.path.splitext(
            os.path.basename(full))[0], smiles)
        return {"path": rel, "kind": "structure", "smiles": smiles,
                "name": name, "how": how,
                "descriptors": dict(mol.descriptors), "parameters": None}
    if ext in PARAMETER_EXTENSIONS:
        params = read_parameter_file(full, line, col)
        unknown = [k for k in params
                   if k not in _PARAM_FILE_KEYS and k != "name"]
        return {"path": rel, "kind": "parameters", "smiles": params.get("smiles"),
                "name": params.get("name"), "how": f"{ext} parameter file",
                "descriptors": None, "parameters": params,
                "unrecognised_keys": unknown}
    raise ResolutionError(
        f"unsupported file type {ext!r} for {os.path.basename(full)}. "
        f"Structures: {', '.join(STRUCTURE_EXTENSIONS)}. Parameters: "
        f"{', '.join(PARAMETER_EXTENSIONS)}.",
        line=line, col=col, code="unsupported-file")


#: Keys a ``.json`` / ``.yaml`` drug file may set. ``smiles`` switches the file
#: from a pure parameter set to a structure the molecular layer can compile.
_PARAM_FILE_KEYS = ("smiles", "ic50", "hill", "diffusivity", "binding",
                    "clearance", "bbb_penetrance", "molecular_weight")


# --------------------------------------------------------------------------- #
# Identity
# --------------------------------------------------------------------------- #
def _canonical(smiles: str) -> Optional[str]:
    try:
        from rdkit import Chem, RDLogger

        RDLogger.DisableLog("rdApp.*")
        mol = Chem.MolFromSmiles(smiles)
        return Chem.MolToSmiles(mol) if mol is not None else None
    except Exception:  # pragma: no cover
        return None


def _alias_reverse_map() -> Dict[str, str]:
    from ..molecular.potency import _GDSC_ALIASES

    out: Dict[str, str] = {}
    for canonical, aliases in _GDSC_ALIASES.items():
        out[canonical] = canonical
        for a in aliases:
            out.setdefault(a.lower(), canonical)
    return out


def resolve_identity(alias: str, smiles: Optional[str],
                     library_name: Optional[str] = None
                     ) -> Tuple[Optional[str], str]:
    """Recover a generic drug name. Returns ``(identity_or_None, how)``.

    Order: an explicit library name; a canonical-SMILES match against the
    curated GBM library; and — **only when no structure was supplied** — the
    GDSC alias table keyed on the TPL identifier.

    The last restriction is the point of this function. A TPL identifier is a
    LOCAL ALIAS: ``drug TMZ from smiles "<something else>"`` is a compound the
    author chose to call TMZ, and matching that label against the GDSC synonym
    list would attach a measured temozolomide potency to a molecule that is not
    temozolomide. When a structure is given, the STRUCTURE is the claim; if it
    matches nothing curated, the identity is unknown and says so.
    """
    if library_name:
        return library_name.strip().lower(), "declared `from library`"

    from ..molecular.library import GBM_DRUGS

    if smiles:
        canon = _canonical(smiles)
        if canon:
            for name, lib_smiles in GBM_DRUGS.items():
                lib_canon = _canonical(lib_smiles)
                if lib_canon and lib_canon == canon:
                    return name, ("canonical-SMILES match against "
                                  "molecular.library.GBM_DRUGS (structural "
                                  "identity, not a name match)")
        return None, (
            f"no identity established: the supplied structure matches no drug "
            f"in molecular.library.GBM_DRUGS. The identifier {alias!r} is a "
            "LOCAL ALIAS and is deliberately NOT matched against the drug "
            "synonym table when a structure is given — a name cannot certify "
            "the identity of a molecule that differs from it.")

    key = (alias or "").strip().lower()
    if key in GBM_DRUGS:
        return key, "the TPL identifier is a curated library drug name"
    rev = _alias_reverse_map()
    if key in rev:
        return rev[key], (f"GDSC alias table: {alias!r} is a recognised "
                          f"synonym for {rev[key]!r} (no structure was "
                          "supplied to check it against)")
    return None, ("no identity established — no structure was supplied and the "
                  "name is not a recognised drug synonym")


# --------------------------------------------------------------------------- #
# The resolver
# --------------------------------------------------------------------------- #
def resolve_drug(node, compile_params: bool = True,
                 alpha: Optional[float] = None,
                 cache_dir: str = "data_cache") -> ResolvedDrug:
    """Resolve one :class:`tmesim.protocol.ast.DrugDef` into a ResolvedDrug.

    ``compile_params=False`` performs every CHEAP check (the library name
    exists, the SMILES parses, the file is readable, ``ic50 from gdsc``
    resolves to a measurement) without running the molecular pipeline, so
    ``POST /api/protocol/validate`` can answer in milliseconds.
    """
    out = ResolvedDrug(alias=node.name, source_kind=node.source_kind,
                       source_value=node.source_value)
    line, col = node.source_line or node.line, node.source_col or node.col

    # ---- 1. the structure / parameter source ------------------------------ #
    file_params: Dict = {}
    library_name: Optional[str] = None
    mol = None
    from ..molecular import descriptors as _desc

    if node.source_kind == "library":
        from ..molecular.library import GBM_DRUGS, get_molecule

        library_name = node.source_value.strip().lower()
        if library_name not in GBM_DRUGS:
            raise ResolutionError(
                f"drug {node.name}: {node.source_value!r} is not in the "
                f"built-in library. Available: {', '.join(sorted(GBM_DRUGS))}.",
                line=line, col=col, code="unknown-drug")
        mol = get_molecule(library_name)
        out.smiles = mol.smiles
    elif node.source_kind == "smiles":
        smiles = node.source_value
        if _canonical(smiles) is None:
            raise ResolutionError(
                f"drug {node.name}: RDKit could not parse the SMILES "
                f"{smiles!r}. A structure that does not parse cannot be "
                "compiled, and no fallback descriptor set is substituted for a "
                "drug the script named by structure.",
                line=line, col=col, code="bad-smiles")
        out.smiles = smiles
        mol = _desc.from_smiles(node.name, smiles)
    elif node.source_kind == "file":
        preview = load_drug_file(node.source_value, line, col)
        out.notes.append(f"loaded {preview['path']} ({preview['how']})")
        if preview["kind"] == "structure":
            out.smiles = preview["smiles"]
            mol = _desc.from_smiles(preview["name"] or node.name,
                                    preview["smiles"])
        else:
            file_params = dict(preview["parameters"] or {})
            if preview.get("unrecognised_keys"):
                out.notes.append(
                    "ignored unrecognised keys in the parameter file: "
                    f"{preview['unrecognised_keys']}")
            smiles = file_params.pop("smiles", None)
            if smiles:
                if _canonical(smiles) is None:
                    raise ResolutionError(
                        f"drug {node.name}: the `smiles` field in "
                        f"{preview['path']} does not parse ({smiles!r}).",
                        line=line, col=col, code="bad-smiles")
                out.smiles = smiles
                mol = _desc.from_smiles(
                    str(file_params.get("name") or node.name), smiles)
            file_params.pop("name", None)
    else:  # pragma: no cover - the parser restricts this
        raise ResolutionError(
            f"drug {node.name}: unknown source {node.source_kind!r}.",
            line=line, col=col, code="unknown-drug-source")

    out.molecule = mol
    out.identity, out.identity_route = resolve_identity(
        node.name, out.smiles, library_name)

    # ---- 2. explicit field overrides from the script ---------------------- #
    from .ast import OVERRIDABLE_FIELDS

    overrides: Dict[str, Tuple[float, str]] = {}
    wants_gdsc: Optional[object] = None
    for fld in node.fields:
        if fld.source is not None:
            wants_gdsc = fld.source
            continue
        q = fld.value
        value, record = q.canonical()
        spec = OVERRIDABLE_FIELDS[fld.name]
        overrides[fld.name] = (float(value), record)
        if fld.name == "bbb_penetrance" and not (0.0 <= value <= 1.0):
            raise ResolutionError(
                f"drug {node.name}: bbb_penetrance must be in [0, 1], got "
                f"{value:g}.", line=q.line, col=q.col, end_col=q.end_col,
                code="out-of-range")
        if value <= 0 and fld.name in ("ic50", "diffusivity", "hill",
                                       "molecular_weight"):
            raise ResolutionError(
                f"drug {node.name}: {fld.name} must be > 0, got {value:g} "
                f"{spec['engine_unit']}.", line=q.line, col=q.col,
                end_col=q.end_col, code="out-of-range")

    # ---- 3. potency ------------------------------------------------------- #
    from ..molecular.potency import resolve_ic50

    ic50_value: Optional[float] = None
    if "ic50" in overrides:
        ic50_value, record = overrides["ic50"]
        out.ic50_resolution = {
            "ic50_uM": ic50_value, "confidence": "user-supplied",
            "source": f"written in the script: {record}", "note": record}
        if wants_gdsc is not None:
            raise ResolutionError(
                f"drug {node.name}: both `ic50 from gdsc` and an explicit "
                "`ic50` value are given. Pick one — a measurement and an "
                "override cannot both be the provenance of one number.",
                line=wants_gdsc.line, col=wants_gdsc.col,
                end_col=wants_gdsc.end_col, code="conflicting-ic50")
    elif wants_gdsc is not None:
        if out.identity is None:
            raise ResolutionError(
                f"drug {node.name}: `ic50 from gdsc` needs a drug identity and "
                f"none could be established ({out.identity_route}). The GDSC "
                "glioma screen is keyed by generic drug name; a novel structure "
                "is not in it. Supply an explicit `ic50 <N> uM` instead — it "
                "will be recorded as user-supplied.",
                line=wants_gdsc.line, col=wants_gdsc.col,
                end_col=wants_gdsc.end_col, code="gdsc-no-identity")
        res = resolve_ic50(out.identity, cache_dir=cache_dir)
        if res.get("confidence") != "measured":
            raise ResolutionError(
                f"drug {node.name}: `ic50 from gdsc` requires a MEASURED value "
                f"and {out.identity!r} is absent from the cached GDSC glioma "
                f"screen (resolve_ic50 fell back to '{res.get('source')}' / "
                f"'{res.get('confidence')}'). Refusing to pass a "
                f"{res.get('confidence')} number off as a measurement. Either "
                "drop `from gdsc` (the run then uses the best available tier "
                "and says so) or supply `ic50 <N> uM` explicitly.",
                line=wants_gdsc.line, col=wants_gdsc.col,
                end_col=wants_gdsc.end_col, code="gdsc-absent")
        ic50_value = float(res["ic50_uM"])
        out.ic50_resolution = dict(res)
        out.ic50_resolution["confidence"] = "measured"
    elif "ic50" in file_params:
        ic50_value = float(file_params["ic50"])
        out.ic50_resolution = {
            "ic50_uM": ic50_value, "confidence": "user-supplied",
            "source": f"uploaded parameter file {node.source_value!r}",
            "note": "taken verbatim from the file"}
    elif out.identity is not None:
        res = resolve_ic50(out.identity, cache_dir=cache_dir)
        ic50_value = float(res["ic50_uM"])
        out.ic50_resolution = dict(res)
    else:
        out.ic50_resolution = {
            "ic50_uM": None, "confidence": "default_LOW",
            "source": "default",
            "note": ("no identity and no explicit ic50: the molecular layer's "
                     "flagged 25 uM default is used. Potency is NOT derivable "
                     "from structure.")}

    # A parameter file with no structure has nothing for the molecular layer to
    # compute from, so the fields it cannot derive must be present. Checked
    # BEFORE the compile gate so `POST /api/protocol/validate` catches it too.
    if mol is None:
        missing = [k for k in ("diffusivity", "molecular_weight")
                   if k not in file_params and k not in overrides]
        if missing:
            raise ResolutionError(
                f"drug {node.name}: this source supplies no structure, so "
                f"{', '.join(missing)} cannot be computed and must be given "
                "explicitly (in the file or in the drug block). No default is "
                "invented for a compound this pipeline has never seen.",
                line=line, col=col, code="incomplete-parameters")

    if not compile_params:
        return out

    # ---- 4. compile ------------------------------------------------------- #
    import dataclasses

    from ..core import DrugParams
    from ..molecular.binding_site_barrier import apply_barrier_to_drugparams
    from ..molecular.compile_drug import molecule_to_drugparams
    from ..validation.calibration import calibrate_alpha, calibrated_drugparams

    if alpha is None:
        alpha = float(calibrate_alpha())
    out.alpha = float(alpha)

    if mol is not None:
        dp_raw = molecule_to_drugparams(mol, ic50=ic50_value)
        computed_src = (
            "molecular.compile_drug.molecule_to_drugparams: Stokes-Einstein "
            "free diffusivity -> Brownian-dynamics ECM hindrance -> "
            "Ogston-Renkin fiber-matrix size gate; permeability/BBB from "
            "CNS-MPO-style physicochemical rules; binding and clearance from "
            "logP/HBD heuristics")
        dp_cal = calibrated_drugparams(dp_raw, alpha)
        dp = apply_barrier_to_drugparams(dp_cal, mol)
        base_tag = "computed"
    else:
        # A pure parameter file: nothing was computed, everything is the user's.
        dp = DrugParams(
            name=node.name,
            diffusivity=float(file_params.get("diffusivity",
                                              DrugParams(name="x").diffusivity)),
            molecular_weight=float(file_params.get("molecular_weight", 0.0)),
            ic50=float(ic50_value if ic50_value is not None else 25.0),
            hill=float(file_params.get("hill", DrugParams(name="x").hill)),
            binding=float(file_params.get("binding", DrugParams(name="x").binding)),
            clearance=float(file_params.get("clearance",
                                            DrugParams(name="x").clearance)),
            bbb_penetrance=float(file_params.get(
                "bbb_penetrance", DrugParams(name="x").bbb_penetrance)),
        )
        computed_src = "engine dataclass default"
        base_tag = "user-supplied"
        out.notes.append(
            "this drug has NO structure, so nothing was computed by the "
            "molecular layer and no perivascular-sink calibration or "
            "binding-site barrier was applied. Every value below is either the "
            "file's or the engine's dataclass default.")
        dp_cal = dp

    # ---- 5. the script's overrides win, applied LAST ---------------------- #
    replaced = {}
    for name, (value, _record) in overrides.items():
        if name == "ic50":
            continue  # already folded in above
        replaced[name] = value
    if replaced:
        dp = dataclasses.replace(dp, **replaced)
    if ic50_value is not None:
        dp = dataclasses.replace(dp, ic50=float(ic50_value))
    dp = dataclasses.replace(dp, name=node.name)
    out.drug_params = dp

    # ---- 6. per-field provenance ------------------------------------------ #
    def _tag(fname: str, computed_note: str, unit: str) -> Dict:
        if fname in overrides:
            return _prov(getattr(dp, fname), unit, "user-supplied",
                         f"written in the script: {overrides[fname][1]}")
        if fname in file_params:
            return _prov(getattr(dp, fname), unit, "user-supplied",
                         f"uploaded parameter file {node.source_value!r}")
        return _prov(getattr(dp, fname), unit, base_tag, computed_note)

    ic50_conf = (out.ic50_resolution or {}).get("confidence")
    if "ic50" in overrides or ic50_conf == "user-supplied":
        ic50_block = _prov(float(dp.ic50), "uM", "user-supplied",
                           (out.ic50_resolution or {}).get("source", "script"))
    elif ic50_conf == "measured":
        n = (out.ic50_resolution or {}).get("n_cell_lines")
        ic50_block = _prov(
            float(dp.ic50), "uM", "measured",
            f"GDSC glioma screen, median across {n} cell lines "
            f"({(out.ic50_resolution or {}).get('note')})")
    elif ic50_conf == "literature":
        ic50_block = _prov(float(dp.ic50), "uM", "flagged-default",
                           "curated LITERATURE anchor, not a measurement: "
                           f"{(out.ic50_resolution or {}).get('note')}")
    else:
        ic50_block = _prov(float(dp.ic50), "uM", "flagged-default",
                           (out.ic50_resolution or {}).get("note")
                           or "flagged low-confidence default")

    out.provenance = {
        "ic50": ic50_block,
        "hill": (_prov(float(dp.hill), "", "user-supplied",
                       "written in the script")
                 if "hill" in overrides else
                 _prov(float(dp.hill), "", "flagged-default",
                       "compile_drug._DEFAULT_HILL = 1.2 — a placeholder; the "
                       "Hill slope is not derived from structure and was not "
                       "measured for this compound")),
        "diffusivity": _tag("diffusivity", computed_src, "um2/s"),
        "molecular_weight": _tag(
            "molecular_weight",
            "RDKit Descriptors.MolWt on the parsed structure", "g/mol"),
        "bbb_penetrance": _tag(
            "bbb_penetrance",
            "molecular.permeability.bbb_penetrance (CNS-MPO-style "
            "physicochemical rules), then molecular.binding_site_barrier",
            ""),
        "binding": _tag(
            "binding",
            (f"validation.calibration.calibrated_drugparams: alpha x "
             f"uptake_proxy with alpha = {out.alpha:.4f} (fit ONCE on the "
             "doxorubicin 40-50 um penetration half-distance, Primeau 2005, "
             "PMID 16361566), plus any binding-site-barrier retention"
             if mol is not None else "engine dataclass default"),
            "1/s"),
        "clearance": _tag(
            "clearance",
            ("validation.calibration.k_passive: the compiled passive "
             "metabolic-loss term" if mol is not None
             else "engine dataclass default"),
            "1/s"),
    }
    if mol is not None:
        out.provenance["_pipeline"] = {
            "value": None, "unit": "", "provenance": "computed",
            "source": ("molecule -> molecule_to_drugparams -> "
                       "calibrated_drugparams(alpha) -> "
                       "apply_barrier_to_drugparams; the identical chain "
                       "app/server/jobs.py runs for a library drug"),
        }
    return out


# --------------------------------------------------------------------------- #
# Capability
# --------------------------------------------------------------------------- #
def check_capability(section_id: str, state_mode: str,
                     sensitivity_mode: str = "uniform",
                     oxygen_mode: str = "pde",
                     line: int = 0, col: int = 0) -> Tuple[Dict, collect]:
    """Validate the requested modes against the specimen's PANEL, before compute.

    Returns ``(capability_report, diagnostics)``. Errors here are raised at plan
    time, so a request the panel cannot support never reaches the solver.

    The rule for ``state_mode='markers'`` is STRICTER than ``POST /api/validate``
    (which warns), and deliberately so: ``PROTOCOL_LANGUAGE.md`` §4 says an
    unsupported mode "is an **error at validation time**, before anything runs".
    A panel carrying no hypoxia-response transcript cannot support the marker
    route — ``results/cell_state_impact.md`` §1.3, §5 — and on the 319-gene
    flagship panel NO route clears the bar at all, which §5 calls "the most
    important sentence in this report".
    """
    diags = collect()
    from app.server.cache import (
        SECTIONS,
        STATE_MODES,
        describe_section,
        state_route_support,
    )

    rec = SECTIONS.get(section_id)
    if rec is None:
        diags.error(
            f"unknown specimen {section_id!r}. Known sections: "
            f"{', '.join(SECTIONS)}. An ingested sample id must be registered "
            "before a protocol can name it.",
            code="unknown-section", line=line, col=col)
        return {"section": section_id, "known": False}, diags

    row = describe_section(section_id)
    support = state_route_support(section_id)
    report = {
        "section": section_id,
        "known": True,
        "label": row.get("label"),
        "panel": row.get("panel"),
        "source": rec.get("source"),
        "available": bool(row.get("available")),
        "n_cells": row.get("n_cells"),
        "state_mode": state_mode,
        "oxygen_mode": oxygen_mode,
        "sensitivity_mode": sensitivity_mode,
        "state_route_support": support,
        "has_outcome_label": rec.get("has_outcome_label", False),
    }

    if not row.get("available"):
        diags.error(
            f"specimen {section_id!r} is not loadable in this checkout: "
            f"{row.get('reason')}",
            code="section-unavailable", line=line, col=col)
        return report, diags

    if state_mode not in STATE_MODES:
        diags.error(f"state_mode {state_mode!r} is not one of "
                    f"{list(STATE_MODES)}.",
                    code="unknown-state-mode", line=line, col=col)
        return report, diags

    if state_mode == "markers":
        if rec.get("source") == "visium":
            diags.error(
                f"state_mode='markers' is unavailable on {section_id!r}: "
                "Visium is spot-level and carries no per-cell transcripts, so "
                "the marker route has nothing to read. Use 'oxygen' or "
                "'vessel_distance'.",
                code="capability", line=line, col=col)
        elif not support.get("genes_known"):
            diags.error(
                f"state_mode='markers' cannot be verified for {section_id!r}: "
                "this checkout cannot read the panel's gene list, so whether "
                "the route is supported is unknown. Refusing rather than "
                "running and finding out.",
                code="capability", line=line, col=col)
        else:
            hypox = support.get("hypoxia_present") or []
            prolif = support.get("proliferation_present") or []
            if not hypox:
                diags.error(
                    f"state_mode='markers' is NOT supported by the panel on "
                    f"{section_id!r} ({row.get('panel')}): it carries none of "
                    f"the hypoxia-response transcripts "
                    f"{support.get('hypoxia_missing')}, so HYPOXIC is "
                    "UNMEASURED and is assigned to NO cell. On this panel NO "
                    "state route clears the bar at all "
                    "(results/cell_state_impact.md §1.3, §5): kill moves by up "
                    "to 0.3252 absolute purely with the choice of route. "
                    "Refusing. Use 'oxygen' (the shipping default) and read the "
                    "coverage MAP of ONE compound. Within a single compound on "
                    "a single tissue the RELATIVE exposure between regions is "
                    "the defensible output; the absolute level is not, and "
                    "comparison BETWEEN compounds is not — compound ordering "
                    "was withdrawn on 2026-08-21 "
                    "(results/pk_ranking_robustness.md, "
                    "results/capacity_ranking_robustness.md).",
                    code="capability", line=line, col=col)
            if not prolif:
                diags.error(
                    f"state_mode='markers' on {section_id!r}: the panel carries "
                    f"none of the proliferation transcripts "
                    f"{support.get('proliferation_missing')}, so the route has "
                    "no proliferation content either.",
                    code="capability", line=line, col=col)
            if hypox and prolif and not support.get("cycle_specific_present"):
                diags.warn(
                    "state_mode='markers': this panel carries no MKI67 or "
                    "TOP2A, so the proliferation call rests on low-specificity "
                    "markers only.",
                    code="capability", line=line, col=col)
    elif state_mode == "vessel_distance":
        diags.warn(
            "state_mode='vessel_distance' bands a geometric SURROGATE "
            "(distance to the nearest measured endothelial cell) — the vessel "
            "positions are measured, the state is not.",
            code="capability", line=line, col=col)

    if state_mode != "oxygen":
        diags.warn(
            f"state_mode={state_mode!r} relabels cells.state on the "
            "RECONSTRUCTION. sim.engine.simulate_response re-bands every living "
            "cell from its own oxygen solve before applying the kill, so this "
            "moves the reported band occupancies and their provenance, NOT this "
            "run's kill_fraction. It is also part of the reconstruction cache "
            "key, so selecting it forces a rebuild.",
            code="capability", line=line, col=col)

    if support.get("genes_known") and not support.get("hypoxia_present"):
        diags.warn(
            f"ABSOLUTE KILL MUST NOT BE QUOTED FROM {section_id!r} UNDER ANY "
            "STATE ROUTE. Its panel carries no hypoxia-response transcript, so "
            "no route clears the bar (results/cell_state_impact.md §5). What "
            "IS defensible is narrower and does not cross compounds: WITHIN "
            "ONE COMPOUND ON THIS ONE TISSUE, the RELATIVE exposure between "
            "regions of the map — this region against that one, in the same "
            "run — is the reading this section supports. The absolute level is "
            "not, and the comparison BETWEEN compounds is not: compound "
            "ordering was withdrawn on 2026-08-21 "
            "(results/pk_ranking_robustness.md, "
            "results/capacity_ranking_robustness.md). "
            "TWO MEASURED INVARIANCES hold that ratio up, and they hold it to "
            "DIFFERENT strengths, which is why they are stated apart. (1) A "
            "per-compound plasma factor is an exact linear rescale of the "
            "whole AUC field — worst relative deviation 1.41e-14 over plasma "
            "factors 0.005 to 26, a 5,200x span, at the shipping "
            "configuration (results/pk_ranking_robustness.md §1.2) — so dose "
            "scale divides out of a within-map ratio exactly, finite window "
            "and all. (2) The tissue-capacity term this engine omits enters as "
            "a pure time rescale: it leaves the CONVERGED gradient shape, the "
            "penetration length and every steady-state number exactly unmoved "
            "(1e-9, results/capacity_ranking_robustness.md §2b) — but NOT a "
            "finite-window quantity, and this run's window is finite, so an "
            "unconverged map carries a residual capacity dependence the 1e-9 "
            "does not cover. "
            "NEITHER INVARIANCE PROTECTS THE TRANSPORT LAW ITSELF. Both are "
            "invariances to SCALE. A wrong effective diffusivity changes the "
            "SHAPE of the gradient and not merely its level, and a ratio "
            "between two regions is precisely what shape moves; the descriptor "
            "path feeding that law is UNSCORED — the shipping engine feeds the "
            "two-route law RDKit descriptors that have never been scored "
            "against measurement (results/descriptor_falsification.md). "
            "The perturbations ACTUALLY SWEPT are: the three cell-state routes "
            "(results/cell_state_impact.md), a uniform and a per-compound "
            "tissue capacity across the measured 3.2-4184 range "
            "(results/capacity_ranking_robustness.md), and the per-compound "
            "plasma factor (results/pk_ranking_robustness.md). That is a LIST, "
            "not a general claim of invariance: an axis absent from it was not "
            "swept, and silence about it is not coverage.",
            code="caveat", line=line, col=col)

    if sensitivity_mode == "uniform":
        diags.warn(
            "sensitivity uniform: every cell is an equal "
            "NO-INTRINSIC-RESISTANCE UPPER BOUND, so the absolute kill is not "
            "a prediction of response. What the flat sensitivity buys is "
            "narrower than the wording it replaces claimed. It removes a "
            "per-cell confound from the comparison of REGIONS WITHIN ONE "
            "COMPOUND'S MAP on this one tissue — the relative exposure between "
            "regions, which is the defensible output. It buys nothing for the "
            "comparison BETWEEN compounds, which was withdrawn on 2026-08-21: "
            "results/pk_ranking_robustness.md measures a 113.1x per-compound "
            "plasma exposure spread against a published GBM order that first "
            "breaks at a spread of 1.53x, and "
            "results/capacity_ranking_robustness.md moves the breast top-1 "
            "pick under this repo's own structural capacity rule. A within-map "
            "ratio is invariant to the dose/plasma SCALE exactly (1.41e-14 "
            "over a 5,200x span of plasma factors, "
            "results/pk_ranking_robustness.md §1.2) and to the omitted tissue "
            "capacity where the field has CONVERGED (1e-9, "
            "results/capacity_ranking_robustness.md §2b) — capacity is a time "
            "rescale, so over a finite unconverged window that invariance is "
            "not exact. It is NOT invariant to error in the transport law: a "
            "wrong effective diffusivity moves the SHAPE of the gradient, and "
            "the descriptor inputs to that law are unscored "
            "(results/descriptor_falsification.md). Those axes, plus the three "
            "cell-state routes (results/cell_state_impact.md), are the "
            "perturbations that were swept; the list is a list, and an axis "
            "not on it was not tested.",
            code="caveat", line=line, col=col)
    elif sensitivity_mode == "determinant" and rec.get("source") != "xenium":
        diags.error(
            f"sensitivity determinant needs per-cell expression and "
            f"{section_id!r} is a {rec.get('source')} section, which has none. "
            "Sensitivity would silently stay at uniform 1.0.",
            code="capability", line=line, col=col)

    if not rec.get("has_outcome_label", False):
        diags.warn(
            f"{section_id!r} is an untreated resection with NO outcome label. "
            "Nothing produced from it has been scored against a measured "
            "treatment response. NOT FOR CLINICAL USE.",
            code="caveat", line=line, col=col)
    return report, diags
