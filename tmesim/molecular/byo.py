"""BRING YOUR OWN COMPOUND: compile a user-supplied SMILES into the same
``core.DrugParams`` a shipped library drug is compiled into -- and return, beside
every field, the honest statement of where that field came from.

Why this module exists
----------------------
``descriptors.from_smiles -> compile_drug.build_profile -> molecule_to_drugparams``
already turns a structure into engine parameters. Two things stop that chain from
being safe to hand to a user directly:

1. **It never refuses.** ``from_smiles`` is deliberately crash-proof: an
   unparseable SMILES with an unknown name yields a Molecule with ``mw=0.0`` and
   every descriptor zero, ``source='fallback_table'``. Downstream that becomes a
   ``DrugParams`` with a nonsense diffusivity rather than an error. That
   behaviour is right for a pipeline that must not die mid-panel; it is exactly
   wrong for a form field. :func:`compile_user_compound` therefore calls
   ``compute_descriptors`` (which DOES raise) and adds the range checks the
   library never needed, because the library only ever contained compounds
   somebody had already vetted.

2. **It does not say which numbers are measurements.** For a user compound the
   answer is stark and worth stating plainly:

   =====================  =====================================================
   field                  where it comes from
   =====================  =====================================================
   mw                     COMPUTED exactly from the supplied graph. The one
                          field here that is not a model output.
   logp                   PREDICTED -- RDKit Crippen ``MolLogP``, an atom-type
                          contribution model, not a measurement. EVERY logP in
                          this repo is a Crippen prediction
                          (``results/logp_provenance.md`` section 1); the
                          Hansch/experimental route is empty for all nine drugs
                          that audit examined, so a user compound is on exactly
                          the same footing as a library one here. Crippen is
                          STEREOCHEMISTRY-BLIND: doxorubicin and epirubicin,
                          which differ only at one centre, produce bit-identical
                          descriptor dicts. If your series varies by
                          stereochemistry, this field cannot see the difference.
                          logP is also the widest-reaching estimate on the list:
                          it feeds free fraction, tissue binding, clearance,
                          membrane permeability AND BBB penetrance.
   diffusivity            PREDICTED, and MEASURED AGAINST DATA. The chain is
                          ``r = 0.066*MW^(1/3)`` -> Stokes-Einstein -> Brownian
                          hindrance -> a two-route cell term that consumes logP
                          and hydrogen-bond donor count through a membrane
                          permeability. The MW-ONLY version of this chain was
                          scored against 71 compounds measured in HT29
                          multicellular layers (Pruijn et al. 2008) and
                          REFUTED: R^2 = 0.022, flat by construction at matched
                          mass where the measurement moves 3.7-3.8x
                          (``results/diffusivity_falsification.md``). The
                          shipping chain scores R^2 = 0.718 on the same 71
                          compounds, blind and unfitted, and orders 24 of 31
                          exact-mass pairs correctly against the old model's
                          0 of 31 (``results/two_route_promotion.md``).
                          WHAT IS STILL NOT TRUE: it reproduces 8.3x of a
                          measured 29.7x spread -- about 24 % -- so a similar-
                          mass series is separated WEAKLY rather than not at
                          all, the whole test is ONE cell line, and no absolute
                          number here is defensible.
   binding / clearance    ASSUMED. Structure heuristics, then OVERWRITTEN by the
                          single global constant ``alpha`` fitted on
                          doxorubicin's penetration anchor. The fit is one drug.
   bbb_penetrance         PREDICTED -- CNS-MPO-style physicochemical rules,
                          known to under-predict temozolomide.
   ic50                   YOURS or ABSENT. Structure cannot give potency. If you
                          do not supply one, a flagged 25 uM placeholder is used
                          and every potency-dependent output is meaningless in
                          absolute terms.
   =====================  =====================================================

What a user may conclude
------------------------
See :data:`CAN_CONCLUDE` / :data:`CANNOT_CONCLUDE`. The short version, which is
the same claim the rest of this project makes for its own library: a compiled
compound may be RANKED against other compounds on the same tissue section under
the same settings. Its absolute AUC, kill fraction and penetration depth may not
be read as predictions.

This module does not modify ``core``, the tissue engine, ``descriptors``,
``compile_drug`` or ``validation.calibration``. It composes them.

Units: MW g/mol, diffusivity um^2/s, IC50 uM, rates 1/s.
"""
from __future__ import annotations

import math
import re
from typing import Dict, List, Optional, Tuple

from ..core import DrugParams
from . import descriptors as _descriptors
from .binding_site_barrier import apply_barrier_to_drugparams
from .compile_drug import build_profile, molecule_to_drugparams
from .types import Molecule, MolecularProfile

__all__ = [
    "CompoundRefused",
    "CompiledCompound",
    "compile_user_compound",
    "check_smiles",
    "CAN_CONCLUDE",
    "CANNOT_CONCLUDE",
    "MW_MIN",
    "MW_MAX",
    "MW_CALIBRATED_MIN",
    "MW_CALIBRATED_MAX",
    "BIOLOGIC_MW_THRESHOLD",
    "MAX_NAME_LEN",
    "MAX_SMILES_LEN",
]


# --------------------------------------------------------------------------- #
# Range of validity -- every bound below is a repo fact, not a round number
# picked for comfort.
# --------------------------------------------------------------------------- #

#: Above this MW ``descriptors.hydrodynamic_radius_nm`` silently switches from
#: the small-molecule density-sphere law to the Erickson globular-PROTEIN fit,
#: and the rest of this chain (Crippen logP, CNS-MPO BBB rules, the
#: solubility-diffusion permeability) has no meaning for a folded protein. This
#: is ``descriptors._LARGE_MW_THRESHOLD``, read from that module so the two can
#: never drift apart.
BIOLOGIC_MW_THRESHOLD = float(_descriptors._LARGE_MW_THRESHOLD)  # 10 000 g/mol

#: The MW window ``compile_drug``'s own bridge self-test asserts every compiled
#: small molecule must land inside (``50 <= mw <= 2000``). Outside it, nothing in
#: this repo has ever checked that the compiled DrugParams is even finite and
#: sane, so a compile is refused rather than reported.
MW_MIN = 50.0
MW_MAX = 2000.0

#: The MW span actually spanned by the two curated libraries this project's
#: transport claims were built and validated on: 5-fluorouracil (130.078,
#: ``molecular.breast_library``) to paclitaxel (853.918, same). Inside this band
#: the MW -> radius law is interpolating between compounds the repo has looked
#: at; outside it (but still inside ``[MW_MIN, MW_MAX]``) it is extrapolating,
#: which is a warning, not a refusal.
MW_CALIBRATED_MIN = 130.078
MW_CALIBRATED_MAX = 853.918

#: Input hygiene. A SMILES longer than this is not a small molecule; a name
#: longer than this is not a name. Both are refused before RDKit is asked, so a
#: paste of a whole FASTA file fails fast with a sentence rather than after a
#: multi-second parse.
MAX_SMILES_LEN = 600
MAX_NAME_LEN = 64

#: A compound name must be usable as a RunSpec ``drug`` value and as a path
#: component in the compound store. Deliberately narrow.
_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9 ._+\-]*$")

#: Sequence-like pastes that are NOT SMILES but often parse into something.
#: A one-letter amino-acid sequence such as "EVQLVESGGGLVQPGG..." is a valid
#: SMILES string character-by-character (C, N, O, S, P, B, F, I are all elements)
#: and RDKit will happily build a nonsense hydrocarbon out of it. Catching that
#: before it becomes a "compound" is the single most likely way a biologic gets
#: silently compiled as a small molecule.
_PROTEIN_LETTERS = set("ACDEFGHIKLMNPQRSTVWY")


# --------------------------------------------------------------------------- #
# What a user may and may not conclude. Stated once, here, and rendered
# verbatim by the API and the Studio -- so the wording cannot drift between the
# docs and the screen.
# --------------------------------------------------------------------------- #
CAN_CONCLUDE: Tuple[str, ...] = (
    "Read this compound against ITSELF on one tissue section: the RELATIVE "
    "exposure between regions of its own map. A per-compound plasma factor is "
    "an exact linear rescale of the whole field, so it cancels from a ratio "
    "between two points of the same map -- which is why a within-map reading "
    "survives the uncertainties an absolute level does not. "
    "DO NOT read an ORDERING across compounds. This entry said to do exactly "
    "that until 2026-08-22, and the cross-compound ranking was measured and "
    "WITHDRAWN on 2026-08-21 (results/pk_ranking_robustness.md, "
    "results/capacity_ranking_robustness.md): per-compound plasma exposure "
    "spans 113.1x across the retained compounds, and there is no per-compound "
    "tissue-capacity term at all. sensitivity_mode='uniform' does make the "
    "per-cell sensitivity term identical for every compound, so it genuinely "
    "cancels -- but cancelling the sensitivity term does not rescue an "
    "ordering that the plasma and capacity terms already broke.",
    "See WHICH property separates it from another compound -- its diffusivity, "
    "its BBB penetrance, its loss rate, its IC50 -- because each is reported "
    "separately and each is a stated function of what you supplied. A difference "
    "driven by an IC50 you measured is on far firmer ground than one driven by "
    "transport.",
    "Detect that two of your compounds are INDISTINGUISHABLE to this model. If "
    "their compiled fields match, the model cannot rank them, and it says so "
    "rather than inventing a difference. This is the most reliable thing the "
    "compiler does, and for a congeneric series it is often the honest answer.",
)

CANNOT_CONCLUDE: Tuple[str, ...] = (
    "Its ABSOLUTE kill fraction. Kill is Hill-driven off an IC50 that, for a "
    "self-compiled compound, is either a number you supplied for your own assay "
    "or a flagged 25 uM placeholder -- and it is further scaled by a per-cell "
    "sensitivity that is 1.0 for every cell. Absolute kill is not a defensible "
    "output of this project for ANY compound, and supplying the molecule "
    "yourself does not make it one.",
    "Its ABSOLUTE AUC or exposure. The transvascular delivery law and the single "
    "alpha calibration set the level; both are project-wide choices fitted on one "
    "drug (doxorubicin), not on yours.",
    "Its ABSOLUTE penetration depth in um. The DIMENSIONLESS hindrance factor "
    "now lands inside its measured band (10.3-fold shipping, against a measured "
    "12-fold in a V79-171b multicellular layer, PMID 9806526) -- but the "
    "ABSOLUTE diffusivity does not: tirapazamine compiles to 49.3 um^2/s against "
    "a measured 74-130 um^2/s, i.e. 1.5x too LOW, having previously been 2.2x "
    "too high. The two anchors disagree about where the engine sits and that is "
    "unresolved (results/two_route_promotion.md).",
    "A TRANSPORT RANKING OF ANALOGUES OF SIMILAR MOLECULAR WEIGHT, AS AN "
    "ABSOLUTE OR A FINE-GRAINED ONE. The MW-only model this engine used to ship "
    "could not resolve such a series AT ALL -- flat by construction at matched "
    "mass, R^2 = 0.022 against 71 compounds in HT29 multicellular layers "
    "(Pruijn et al. 2008; results/diffusivity_falsification.md). The shipping "
    "two-route model does resolve it, at R^2 = 0.718 blind on those same 71 "
    "compounds with 24 of 31 exact-mass pairs ordered correctly "
    "(results/two_route_promotion.md). But it reproduces 8.3x of a measured "
    "29.7x spread -- roughly a QUARTER -- so it under-separates a congeneric "
    "series systematically, and it has been tested against exactly ONE cell "
    "line. Read the ORDER over a mass-matched series as weak evidence, not as a "
    "measurement, and read no absolute diffusivity from it at all.",
    "Anything about a stereoisomer, a salt form, a prodrug or a formulation. "
    "Crippen logP is stereochemistry-blind and there is no metabolism step, so an "
    "enantiomer pair compiles to bit-identical parameters and a prodrug compiles "
    "as its parent.",
    "Anything clinical. Nothing in this project has been scored against a patient "
    "outcome.",
)


# --------------------------------------------------------------------------- #
# Refusal
# --------------------------------------------------------------------------- #
class CompoundRefused(ValueError):
    """A user compound this pipeline will not compile, with a reason.

    ``code`` is a stable machine token (the UI keys help text off it);
    ``message`` is a complete sentence written for the person at the form;
    ``field`` names the input at fault so the form can mark it invalid.
    """

    def __init__(self, code: str, message: str, field: str = "smiles",
                 detail: Optional[Dict] = None):
        super().__init__(message)
        self.code = code
        self.message = message
        self.field = field
        self.detail = detail or {}

    def as_dict(self) -> Dict:
        return {"code": self.code, "message": self.message,
                "field": self.field, "detail": self.detail}


# --------------------------------------------------------------------------- #
# Input gates
# --------------------------------------------------------------------------- #
def _check_name(name: str, reserved: Tuple[str, ...] = ()) -> str:
    if not isinstance(name, str) or not name.strip():
        raise CompoundRefused(
            "name_empty",
            "A compound needs a name. It labels the compound everywhere it "
            "appears -- in the drug picker, in the run record and in any report "
            "you export.",
            field="name")
    n = name.strip()
    if len(n) > MAX_NAME_LEN:
        raise CompoundRefused(
            "name_too_long",
            f"Name is {len(n)} characters; the limit is {MAX_NAME_LEN}.",
            field="name")
    if not _NAME_RE.match(n):
        raise CompoundRefused(
            "name_charset",
            "Name must start with a letter or digit and contain only letters, "
            "digits, spaces, and . _ + - characters.",
            field="name")
    low = n.strip().lower()
    if low in {r.strip().lower() for r in reserved}:
        raise CompoundRefused(
            "name_collision",
            f"{n!r} is already the name of a curated compound in this build. "
            "Pick a different name: if a user compound could shadow a library "
            "one, a run record would no longer say unambiguously which structure "
            "it simulated.",
            field="name")
    return n


def _check_logp_tier(mol: Molecule) -> None:
    """Backstop against a NAME-keyed measured-logP lookup, if one ever lands.

    Today this is a no-op: ``descriptors.from_smiles`` derives logP purely from
    the structure (RDKit Crippen), so a name cannot influence it and there is
    nothing to guard. The guard exists anyway because a measured-logP tier keyed
    on drug NAME has been proposed for this repo more than once, and it is safe
    for a curated library (somebody checked that the name and the structure are
    the same drug) and NOT safe for a form field: a user could type any
    structure, call it "aspirin", and silently receive aspirin's published logP
    on a molecule that is not aspirin -- which then propagates into five
    estimators (free fraction, tissue binding, clearance, membrane permeability,
    BBB penetrance) before anything else in the pipeline runs.

    So if a ``logp_tier`` of ``'experimental'`` ever appears on a user compound,
    it is refused rather than compiled. The remedy is a rename, which costs
    nothing; the alternative is a measurement attributed to the wrong graph,
    which cannot be undone downstream. Costs one dict lookup per compile.
    """
    if str(mol.get("logp_tier") or "") == "experimental":
        raise CompoundRefused(
            "name_matches_measured_logp",
            f"{mol.name!r} matches an entry in this build's MEASURED logP table, "
            "which is looked up by NAME. This compiler cannot verify that the "
            "name you typed describes the structure you typed, so it will not "
            "attach somebody else's published measurement to your molecular "
            "graph. Rename the compound and it will compile with the RDKit "
            "Crippen estimate, which is what a self-compiled compound is "
            "entitled to.",
            field="name",
            detail={"logp_source": mol.get("logp_source")})


def _check_ic50(ic50_uM) -> Optional[float]:
    if ic50_uM is None or ic50_uM == "":
        return None
    try:
        v = float(ic50_uM)
    except (TypeError, ValueError):
        raise CompoundRefused(
            "ic50_not_a_number",
            f"IC50 must be a number in uM, got {ic50_uM!r}. Leave it blank if "
            "you do not have one.",
            field="ic50_uM")
    if not math.isfinite(v):
        raise CompoundRefused(
            "ic50_not_finite",
            f"IC50 must be a finite number in uM, got {v!r}.",
            field="ic50_uM")
    if v <= 0.0:
        raise CompoundRefused(
            "ic50_not_positive",
            f"IC50 must be > 0 uM, got {v}. The Hill kill term divides by it.",
            field="ic50_uM")
    if v > 1.0e6:
        raise CompoundRefused(
            "ic50_absurd",
            f"IC50 of {v:g} uM is 1 M. Check the units -- this field is uM, not "
            "nM and not M.",
            field="ic50_uM")
    return v


def _looks_like_a_protein_sequence(s: str) -> bool:
    """True for a one-letter amino-acid paste masquerading as a SMILES.

    Only ever consulted AFTER RDKit has already failed to parse the string, so
    it can never reject something that is a valid structure -- its whole job is
    to turn "unparseable" into a message that names the actual mistake.

    A peptide/protein sequence is long, all letters from the 20-letter alphabet,
    carries none of the punctuation any real SMILES of that length would have
    (rings, branches, bonds, brackets), and -- the discriminating test -- uses
    MANY DISTINCT letters. A long alkane (``CCCCCC...``) satisfies every other
    condition with a single distinct letter, so the distinct-letter floor is what
    keeps this from firing on real chemistry.
    """
    if len(s) < 25:
        return False
    if any(c in s for c in "()[]=#@+-.%/\\0123456789"):
        return False
    if not (s.isalpha() and all(c in _PROTEIN_LETTERS for c in s.upper())):
        return False
    return len(set(s.upper())) >= 8


def check_smiles(smiles: str) -> "object":
    """Parse ``smiles`` or raise :class:`CompoundRefused`. Returns the RDKit mol.

    Every refusal here is a case where the downstream chain would otherwise
    produce a number that looks like an answer.
    """
    if not isinstance(smiles, str) or not smiles.strip():
        raise CompoundRefused(
            "smiles_empty",
            "Enter a SMILES string. Without a structure there is nothing to "
            "compile -- this pipeline derives every transport parameter from the "
            "molecular graph.")
    s = smiles.strip()
    if len(s) > MAX_SMILES_LEN:
        raise CompoundRefused(
            "smiles_too_long",
            f"That SMILES is {len(s)} characters ({MAX_SMILES_LEN} is the "
            "limit). Anything that long is a biologic or a polymer, and the "
            "descriptor chain here is small-molecule only.")
    if not getattr(_descriptors, "_HAVE_RDKIT", False):
        raise CompoundRefused(
            "rdkit_unavailable",
            "RDKit is not available in this deployment, so no structure can be "
            "parsed. This is a server-side problem, not a problem with your "
            "SMILES.")

    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdMolDescriptors

    RDLogger.DisableLog("rdApp.*")     # RDKit's parse errors go in our message
    try:
        mol = Chem.MolFromSmiles(s)
    except Exception as exc:  # noqa: BLE001 - RDKit can throw on exotic input
        raise CompoundRefused(
            "smiles_unparseable",
            f"RDKit could not parse that SMILES ({type(exc).__name__}: {exc}).")
    finally:
        RDLogger.EnableLog("rdApp.*")

    if mol is None:
        # Name the actual mistake where we can recognise it. A pasted one-letter
        # amino-acid sequence is by far the most likely way somebody tries to
        # put a biologic through a small-molecule compiler, and "unparseable"
        # would send them looking for a typo that is not there.
        if _looks_like_a_protein_sequence(s):
            raise CompoundRefused(
                "smiles_is_a_sequence",
                "That is a one-letter amino-acid sequence, not a SMILES. "
                "Peptides, proteins and antibodies are out of range: every "
                "descriptor in this compiler (Crippen logP, the CNS-MPO BBB "
                "rules, the MW-to-radius law) is a small-molecule model, and "
                "none of them means anything for a folded chain.")
        raise CompoundRefused(
            "smiles_unparseable",
            "RDKit could not parse that SMILES. Check for an unclosed ring "
            "index, an unbalanced bracket, or a copy/paste that picked up a "
            "name or an InChI instead of a SMILES.")
    if mol.GetNumHeavyAtoms() == 0:
        raise CompoundRefused(
            "smiles_no_atoms",
            "That SMILES parses but contains no heavy atoms, so there is no "
            "molecule to compute descriptors from.")

    frags = Chem.GetMolFrags(mol)
    if len(frags) > 1:
        raise CompoundRefused(
            "smiles_multi_component",
            f"That SMILES has {len(frags)} disconnected components (a salt, a "
            "co-crystal or a mixture). Every transport parameter here is derived "
            "from a single molecular weight, so a counter-ion would be carried "
            "into the diffusivity as if it travelled bonded to the drug. Supply "
            "the free acid/base -- the largest component -- on its own.")

    # Crippen logP is an atom-type contribution model. Atoms outside its
    # parameterisation contribute exactly (0.0, 0.0), i.e. the model has nothing
    # to say about them and silently says zero. That is the one failure mode of
    # this descriptor that produces a plausible-looking number from no
    # information at all, so it is refused rather than flagged.
    contribs = rdMolDescriptors._CalcCrippenContribs(mol)
    unparameterised = sorted({
        a.GetSymbol() for a, (logp, mr) in zip(mol.GetAtoms(), contribs)
        if logp == 0.0 and mr == 0.0
    })
    if unparameterised:
        raise CompoundRefused(
            "crippen_unparameterised_atom",
            "The Crippen logP model has no atom type for "
            + ", ".join(unparameterised)
            + ". It does not fail on such an atom -- it assigns it a "
            "contribution of exactly zero, so the compound would receive a "
            "confident-looking logP computed as if those atoms were not there. "
            "Refused rather than reported.",
            detail={"elements": unparameterised})
    return mol


def _check_mw(mw: float, is_antibody: bool) -> List[str]:
    """Refuse out-of-range MW; return warnings for in-range-but-extrapolated."""
    warnings: List[str] = []
    if is_antibody or mw >= BIOLOGIC_MW_THRESHOLD:
        raise CompoundRefused(
            "biologic_out_of_range",
            f"MW {mw:,.0f} g/mol is at or above {BIOLOGIC_MW_THRESHOLD:,.0f}, "
            "where this pipeline stops treating a molecule as a small molecule: "
            "the hydrodynamic radius switches to a globular-protein fit, and "
            "Crippen logP, the CNS-MPO BBB rules and the solubility-diffusion "
            "permeability have no meaning for a folded protein. The one antibody "
            "in the shipped library (bevacizumab) is not compiled from a "
            "structure at all -- it is a hand-entered row of known protein "
            "values, which is exactly what a SMILES cannot reproduce. Antibodies "
            "and biologics are out of range.")
    if mw > MW_MAX:
        raise CompoundRefused(
            "mw_above_range",
            f"MW {mw:,.1f} g/mol is above the {MW_MAX:,.0f} g/mol ceiling this "
            "bridge asserts for a compiled small molecule. Nothing in this repo "
            "has ever checked that the MW -> radius -> diffusivity chain stays "
            "sane up there.")
    if mw < MW_MIN:
        raise CompoundRefused(
            "mw_below_range",
            f"MW {mw:,.1f} g/mol is below the {MW_MIN:,.0f} g/mol floor this "
            "bridge asserts. The MW -> radius prefactor (0.066) was set to "
            "reproduce ~0.4 nm at ~200 g/mol; a fragment this small is far off "
            "the end of that calibration.")
    if mw < MW_CALIBRATED_MIN or mw > MW_CALIBRATED_MAX:
        warnings.append(
            f"MW {mw:,.1f} g/mol is OUTSIDE the {MW_CALIBRATED_MIN:g}-"
            f"{MW_CALIBRATED_MAX:g} g/mol span of the two curated libraries this "
            "project's transport claims were built on (5-fluorouracil to "
            "paclitaxel). The MW -> hydrodynamic-radius law is extrapolating "
            "here, not interpolating.")
    return warnings


# --------------------------------------------------------------------------- #
# The compile
# --------------------------------------------------------------------------- #
class CompiledCompound:
    """A user compound compiled to engine parameters, plus per-field provenance.

    Attributes
    ----------
    name, smiles : str
    molecule : Molecule
    profile : MolecularProfile
    drug_params : DrugParams
        The EXACT object a run would simulate, alpha-calibrated and passed
        through the binding-site barrier -- not the raw bridge output. Reporting
        the raw output would show the user numbers the engine never sees.
    fields : list of dict
        One row per compiled field: ``{key, label, value, unit, kind, method}``
        where ``kind`` is the Studio's provenance vocabulary
        (``computed`` / ``assumed``; nothing here is ``measured``).
    warnings, caveats : list of str
    """

    __slots__ = ("name", "smiles", "molecule", "profile", "drug_params",
                 "fields", "warnings", "caveats", "alpha", "ic50_supplied")

    def __init__(self, **kw):
        for k in self.__slots__:
            setattr(self, k, kw.get(k))

    def as_dict(self) -> Dict:
        dp = self.drug_params
        return {
            "name": self.name,
            "smiles": self.smiles,
            "descriptors": {k: (bool(v) if k == "is_antibody" else v)
                            for k, v in self.molecule.descriptors.items()},
            "drug_params": {
                "name": dp.name,
                "diffusivity_um2_s": float(dp.diffusivity),
                "molecular_weight": float(dp.molecular_weight),
                "ic50_uM": float(dp.ic50),
                "hill": float(dp.hill),
                "binding_per_s": float(dp.binding),
                "clearance_per_s": float(dp.clearance),
                "bbb_penetrance": float(dp.bbb_penetrance),
            },
            "fields": self.fields,
            "warnings": list(self.warnings),
            "caveats": list(self.caveats),
            "can_conclude": list(CAN_CONCLUDE),
            "cannot_conclude": list(CANNOT_CONCLUDE),
            "ic50_supplied": bool(self.ic50_supplied),
            "alpha": self.alpha,
            # How many of the compiled fields are model predictions rather than
            # exact functions of the graph the user typed. The UI leads with
            # this number: "8 of 11 fields are predictions" is the honest
            # one-line summary of a self-compiled compound.
            "predicted_field_count": sum(
                1 for f in self.fields if f.get("is_prediction")),
        }


def _fields(mol: Molecule, prof: MolecularProfile, dp: DrugParams,
            ic50_supplied: bool, alpha: Optional[float]) -> List[Dict]:
    """The per-field provenance table. Every row states its method in words the
    person at the form can check, and carries the Studio provenance ``kind``.

    Nothing here is ``measured``: no instrument touched this compound in this
    build. ``computed`` is reserved for the two quantities that are exact
    functions of the graph the user typed (MW, and the descriptor counts);
    everything downstream of an empirical fit is ``assumed``.
    """
    d = mol.descriptors
    rows: List[Dict] = [
        {
            "key": "mw", "label": "Molecular weight", "unit": "g/mol",
            "value": _f(d.get("mw")), "kind": "computed",
            "is_prediction": False,
            "method": "Exact sum of atomic masses over the graph you supplied "
                      "(RDKit Descriptors.MolWt). This is the only field on this "
                      "list that is not a model output -- and it is also the "
                      "field the whole transport chain hangs on.",
        },
        {
            "key": "logp", "label": "logP (octanol/water)", "unit": "",
            "value": _f(d.get("logp")), "kind": "assumed",
            "is_prediction": True,
            "method": "PREDICTION, not a measurement: RDKit Crippen MolLogP, an "
                      "atom-type contribution model (Wildman & Crippen 1999). "
                      "Crippen is STEREOCHEMISTRY-BLIND -- two enantiomers give "
                      "bit-identical descriptors, so a series that varies by "
                      "stereochemistry is invisible to this field. It is also "
                      "the WEAKER of the two tiers this build carries: several "
                      "curated library drugs use a published octanol/water "
                      "MEASUREMENT instead, and the gap is large (methotrexate: "
                      "Crippen +0.268 vs measured -1.85). A self-compiled "
                      "compound cannot reach that tier, because the measurement "
                      "is looked up by NAME and nothing here can confirm your "
                      "name describes your structure. logP feeds free fraction, "
                      "tissue binding, clearance, membrane permeability AND BBB "
                      "penetrance, so this is the single widest-reaching "
                      "estimate on this list.",
        },
        {
            "key": "tpsa", "label": "TPSA", "unit": "A^2",
            "value": _f(d.get("tpsa")), "kind": "computed",
            "is_prediction": False,
            "method": "Topological polar surface area: a deterministic count of "
                      "polar-atom contributions over the graph (Ertl fragment "
                      "table). Exact given the structure.",
        },
        {
            "key": "radius_hydrodynamic_nm", "label": "Hydrodynamic radius",
            "unit": "nm", "value": _f(d.get("radius_hydrodynamic_nm")),
            "kind": "assumed", "is_prediction": True,
            "method": "PREDICTED FROM MW ALONE: r = 0.066 * MW^(1/3), a "
                      "density-sphere estimate whose prefactor was chosen to "
                      "give ~0.4 nm at ~200 g/mol. Shape, logP, charge and "
                      "H-bonding do not enter it. Two isomers of equal mass get "
                      "an identical radius by construction.",
        },
        {
            "key": "free_diffusivity", "label": "Free diffusivity",
            "unit": "um^2/s", "value": _f(prof.free_diffusivity),
            "kind": "assumed", "is_prediction": True,
            "method": "Stokes-Einstein D = kB*T/(6*pi*eta*r) at 310 K with "
                      "interstitial eta = 1.2e-3 Pa.s. The physics is exact; the "
                      "radius it consumes is the MW-only prediction above, so "
                      "this inherits that prediction.",
        },
        {
            "key": "diffusivity", "label": "Tissue diffusivity (engine input)",
            "unit": "um^2/s", "value": _f(dp.diffusivity),
            "kind": "assumed", "is_prediction": True,
            "method": "MEASURED AGAINST DATA, PARTLY SUPPORTED -- read this "
                      "before you rank on it. Free diffusivity x 2D "
                      "Brownian-dynamics micro-tortuosity x an Ogston-Renkin "
                      "ECM fibre-mesh size gate x a two-route cell term "
                      "(paracellular between cells + transcellular through "
                      "membranes, Hasselman & Johnson 1987), which is what "
                      "brings your compound's logP and hydrogen-bond donor "
                      "count into its diffusivity. Scored BLIND and UNFITTED "
                      "against 71 compounds measured in HT29 multicellular "
                      "layers (Pruijn et al. 2008): R^2 = 0.718, Spearman "
                      "+0.821, 99% within 3x of measurement, geometric bias "
                      "1.03x, and 24 of 31 exact-mass pairs ordered correctly. "
                      "The MW-only model this replaced scored R^2 = 0.022 and "
                      "ordered 0 of 31 (results/diffusivity_falsification.md, "
                      "results/two_route_promotion.md). WHAT IS STILL NOT TRUE, "
                      "and it is the part that matters for a congeneric series: "
                      "the model reproduces 8.3x of a measured 29.7x spread, "
                      "i.e. about a QUARTER, so it UNDER-SEPARATES similar-mass "
                      "analogues systematically rather than separating them "
                      "correctly; the entire test is ONE cell line; and the "
                      "absolute value here is not defensible -- the same drug's "
                      "compiled D_eff misses the measured multicellular-layer "
                      "band (74-130 um^2/s for tirapazamine) by 1.5x on the low "
                      "side. The dimensionless hindrance factor, by contrast, "
                      "now lands INSIDE its measured band (10.3-fold against a "
                      "measured 12-fold, results/cellular_exclusion.md).",
        },
        {
            "key": "bbb_penetrance", "label": "BBB penetrance", "unit": "frac",
            "value": _f(dp.bbb_penetrance), "kind": "assumed",
            "is_prediction": True,
            "method": "PREDICTION from CNS-MPO-style physicochemical rules on "
                      "MW/logP/TPSA/HBD. Structure-only: it cannot see active "
                      "transport or efflux, and it is known to under-predict "
                      "temozolomide. It enters the solve as a LINEAR multiplier "
                      "on the vascular source, so an error here scales exposure "
                      "one-for-one.",
        },
        {
            "key": "binding", "label": "Tissue loss: binding", "unit": "1/s",
            "value": _f(dp.binding), "kind": "assumed", "is_prediction": True,
            "method": (
                "ASSUMED. The structure heuristic (logP + H-bond donors) supplies "
                "only an ORDERING; its absolute scale is then overwritten by "
                "alpha * uptake_proxy, where alpha"
                + (f" = {alpha:.4g}" if alpha is not None else "")
                + " is the project's SINGLE free constant, fitted on ONE drug "
                "(doxorubicin's 40-50 um half-distance anchor). Your compound "
                "inherits that one fit unchanged."
            ),
        },
        {
            "key": "clearance", "label": "Tissue loss: clearance", "unit": "1/s",
            "value": _f(dp.clearance), "kind": "assumed", "is_prediction": True,
            "method": "ASSUMED. A mild logP/MW heuristic standing in for "
                      "metabolism. It is a weak prior, not a metabolic "
                      "prediction, and no metabolism or prodrug step exists "
                      "anywhere in this pipeline.",
        },
        {
            "key": "ic50", "label": "IC50", "unit": "uM", "value": _f(dp.ic50),
            "kind": "literature" if ic50_supplied else "assumed",
            "is_prediction": not ic50_supplied,
            "method": (
                "YOURS. Supplied with this compound and used exactly as given; "
                "this pipeline did not measure it and cannot check it. It is "
                "potency in YOUR assay, on YOUR cells -- not on this section."
                if ic50_supplied else
                "ABSENT. You supplied no IC50 and structure alone cannot give "
                "potency without docking or an assay, so a flagged 25 uM "
                "PLACEHOLDER is in use. It is not a property of your compound. "
                "Any kill fraction or IC50-coverage number computed from it is "
                "meaningless in absolute terms, and a ranking against other "
                "compounds is only a DELIVERY ranking, not a potency one."
            ),
        },
        {
            "key": "hill", "label": "Hill coefficient", "unit": "",
            "value": _f(dp.hill), "kind": "assumed", "is_prediction": True,
            "method": "ASSUMED constant 1.2 for every compiled compound. It is "
                      "not derived from your structure and was not fitted to "
                      "your compound.",
        },
    ]
    return rows


def _f(v):
    if v is None:
        return None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None


def compile_user_compound(
    name: str,
    smiles: str,
    ic50_uM=None,
    alpha: Optional[float] = None,
    ecm_reference: float = 0.4,
    reserved_names: Tuple[str, ...] = (),
) -> CompiledCompound:
    """Compile a user SMILES into the DrugParams a run would actually simulate.

    Refuses -- with a :class:`CompoundRefused` carrying a specific sentence --
    on: an empty/over-long/unparseable SMILES, a pasted amino-acid sequence, a
    multi-component (salt/mixture) SMILES, an element Crippen has no atom type
    for, an antibody or anything at/above 10 kDa, an MW outside [50, 2000]
    g/mol, a bad name, a name that collides with a curated compound, a name that
    reaches the MEASURED logP table (see :func:`_check_logp_tier` -- that lookup
    is keyed by name and would otherwise attach somebody else's published
    measurement to this graph), and a non-positive or non-finite IC50.

    Parameters
    ----------
    alpha : float, optional
        The project's fitted loss-rate constant. If given, the returned
        ``drug_params`` is calibrated with it and passed through the
        binding-site barrier, i.e. it is byte-for-byte the object the run
        pipeline builds for a library drug. If ``None`` the raw bridge output is
        returned and a caveat says so.
    ecm_reference : float
        ECM volume fraction the effective diffusivity is computed at. Defaults to
        the same 0.4 the bridge uses.
    """
    nm = _check_name(name, reserved=reserved_names)
    ic50 = _check_ic50(ic50_uM)
    smi = (smiles or "").strip()

    check_smiles(smi)                     # raises CompoundRefused

    try:
        desc = _descriptors.compute_descriptors(smi)
    except ValueError as exc:             # belt-and-braces; check_smiles caught it
        raise CompoundRefused("smiles_unparseable", str(exc))

    warnings = _check_mw(float(desc.get("mw", 0.0) or 0.0),
                         bool(desc.get("is_antibody", False)))

    # Build the Molecule the same way the library does, so a user compound and a
    # library compound differ ONLY in where the SMILES came from.
    mol = _descriptors.from_smiles(nm, smi)
    if mol.source != "rdkit":
        raise CompoundRefused(
            "descriptor_fallback",
            "The descriptor build fell back to the curated table instead of "
            "reading your structure, which would mean the compiled parameters "
            "are not your compound's. Refused.")
    _check_logp_tier(mol)

    try:
        prof = build_profile(mol, ecm_reference=ecm_reference, ic50=ic50)
        dp = molecule_to_drugparams(prof, ic50=ic50 if ic50 is not None else None)
    except Exception as exc:  # noqa: BLE001
        raise CompoundRefused(
            "compile_failed",
            f"The molecular pipeline failed on this structure "
            f"({type(exc).__name__}: {exc}). It parsed, but a downstream physics "
            "step could not produce a finite parameter for it.")

    caveats: List[str] = []
    if alpha is not None:
        # The SAME two steps the run pipeline applies to a library drug:
        # calibration.calibrated_drugparams (documented to pass a DrugParams
        # through _compile unchanged and then rewrite the loss terms) followed by
        # the binding-site barrier. Composed here rather than reimplemented, so
        # a user compound can never drift from a library one.
        from ..validation.calibration import calibrated_drugparams

        dp = apply_barrier_to_drugparams(calibrated_drugparams(dp, alpha), mol)
    else:
        caveats.append(
            "This preview was compiled WITHOUT the project's alpha calibration, "
            "so its binding/clearance are the raw structure heuristics rather "
            "than the loss rates a run would use.")

    # Final sanity: never hand back a DrugParams the engine would choke on.
    for fld in ("diffusivity", "molecular_weight", "ic50", "hill", "binding",
                "clearance", "bbb_penetrance"):
        v = getattr(dp, fld)
        if not (isinstance(v, (int, float)) and math.isfinite(float(v))):
            raise CompoundRefused(
                "compile_nonfinite",
                f"The compiled parameter {fld!r} is not a finite number "
                f"({v!r}). This structure parses but does not produce usable "
                "physics; it is refused rather than simulated.")
    if dp.diffusivity <= 0.0:
        raise CompoundRefused(
            "compile_nonfinite",
            f"The compiled tissue diffusivity is {dp.diffusivity!r} um^2/s. A "
            "non-positive diffusivity is not a transport parameter; refused.")

    # The wording is deliberate: the Studio's caveat renderer promotes a caveat
    # to a HARD STOP (tinted field, stop glyph) on the phrases "NOT DEFENSIBLE"
    # and "MUST NOT". A self-compiled compound's absolute outputs are exactly
    # that, so the caveat is written to trip that rule rather than to read
    # gently.
    caveats.append(
        "TRANSPORT RANKING OVER SIMILAR-MASS ANALOGUES IS WEAK EVIDENCE, NOT A "
        "MEASUREMENT. Diffusivity here consumes molecular weight, logP and "
        "hydrogen-bond donor count through a two-route cell term. Scored blind "
        "against 71 compounds measured in HT29 multicellular layers (Pruijn et "
        "al. 2008, results/two_route_promotion.md) it reaches R^2 = 0.718 and "
        "orders 24 of 31 exact-mass pairs correctly -- against R^2 = 0.022 and "
        "0 of 31 for the molecular-weight-only model this replaced "
        "(results/diffusivity_falsification.md). It nonetheless reproduces only "
        "8.3x of the measured 29.7x spread, so it under-separates such a series "
        "systematically, and it has been tested in ONE cell line. Treat the "
        "order as a hypothesis to test, and do not quote the absolute "
        "diffusivity at all."
    )
    caveats.append(
        "logP here is an RDKit Crippen PREDICTION, not a measurement -- as it is "
        "for every drug in this repo (results/logp_provenance.md). Crippen is "
        "stereochemistry-blind, so an enantiomer pair compiles identically, and "
        "logP drives free fraction, tissue binding, clearance, membrane "
        "permeability and BBB penetrance."
    )
    caveats.append(
        "SELF-COMPILED COMPOUND: every parameter except molecular weight and "
        "TPSA is a PREDICTION from the structure you supplied, not a "
        "measurement. Its absolute kill fraction, AUC and penetration depth are "
        "NOT DEFENSIBLE and MUST NOT be quoted as results. Use it to RANK "
        "compounds against each other on one tissue section; that ranking is the "
        "claim this project defends."
    )
    if ic50 is None:
        caveats.append(
            "NO IC50 WAS SUPPLIED for this compound, so a flagged 25 uM "
            "placeholder stands in. It is not a property of your molecule: any "
            "kill fraction or IC50-coverage number derived from it MUST NOT be "
            "read as potency, and a ranking that includes this compound is a "
            "DELIVERY ranking only.")
    else:
        caveats.append(
            f"The IC50 of {ic50:g} uM is YOURS, taken as given. It was measured "
            "on your cells in your assay, not on this section, and this "
            "pipeline neither verified it nor knows the assay conditions.")
    caveats.extend(warnings)

    return CompiledCompound(
        name=nm, smiles=smi, molecule=mol, profile=prof, drug_params=dp,
        fields=_fields(mol, prof, dp, ic50 is not None, alpha),
        warnings=warnings, caveats=caveats, alpha=alpha,
        ic50_supplied=ic50 is not None,
    )
