"""Physicochemical descriptor computation for the molecular layer.

Turns a drug SMILES into the cheminformatics descriptors the downstream
physics estimators need (Stokes-Einstein radius, membrane permeability, BBB
rules, protein binding). RDKit is the primary engine; every entry point
degrades gracefully to a curated fallback table so the pipeline never crashes
on a bad/absent SMILES (e.g. an antibody marker with no small-molecule graph).

Units: mw g/mol, tpsa A^2, radius_hydrodynamic_nm in nm.
"""
from __future__ import annotations

from typing import Dict

from .types import Molecule

# --- RDKit (primary path) -------------------------------------------------
try:
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors

    _HAVE_RDKIT = True
except Exception:  # pragma: no cover - RDKit is installed in this env
    _HAVE_RDKIT = False


# --- MW -> hydrodynamic radius --------------------------------------------
# Two documented, distinct scalings:
#   small molecules: r_nm = 0.066 * mw**(1/3)
#       density-based sphere estimate (r ~ (3M/4*pi*rho*N_A)^(1/3)); the 0.066
#       prefactor reproduces ~0.4 nm for a ~200 g/mol drug and ~0.5 nm for
#       doxorubicin (543 g/mol), matching measured small-molecule radii.
#   proteins/antibodies: r_nm = 0.0515 * mw**0.392
#       Erickson 2009 empirical relation for compact globular proteins;
#       gives ~5.5 nm for a 149 kDa IgG, matching the measured IgG radius.
_LARGE_MW_THRESHOLD = 10_000.0


def hydrodynamic_radius_nm(mw: float, is_antibody: bool = False) -> float:
    """Empirical MW -> hydrodynamic radius (nm). See module docstring for the
    two scaling laws (small-molecule density sphere vs. Erickson protein fit)."""
    if is_antibody or mw >= _LARGE_MW_THRESHOLD:
        return 0.0515 * float(mw) ** 0.392
    return 0.066 * float(mw) ** (1.0 / 3.0)


# --- curated fallback descriptor table ------------------------------------
# Used only when RDKit cannot parse a SMILES (or none exists, e.g. antibody).
# Values pre-computed from RDKit for the library drugs so a fallback still
# yields sane physics. Keyed by lowercased molecule name.
_FALLBACK_TABLE: Dict[str, Dict[str, float]] = {
    "temozolomide": {"mw": 194.15, "logp": -2.08, "tpsa": 108.17, "hbd": 1, "hba": 7, "rotatable_bonds": 1, "aromatic_rings": 2, "formal_charge": 0, "is_antibody": False},
    "lomustine": {"mw": 233.70, "logp": 2.25, "tpsa": 61.77, "hbd": 1, "hba": 3, "rotatable_bonds": 4, "aromatic_rings": 0, "formal_charge": 0, "is_antibody": False},
    "carmustine": {"mw": 214.05, "logp": 1.16, "tpsa": 61.77, "hbd": 1, "hba": 3, "rotatable_bonds": 5, "aromatic_rings": 0, "formal_charge": 0, "is_antibody": False},
    "procarbazine": {"mw": 221.30, "logp": 1.05, "tpsa": 53.16, "hbd": 3, "hba": 3, "rotatable_bonds": 5, "aromatic_rings": 1, "formal_charge": 0, "is_antibody": False},
    "methotrexate": {"mw": 454.45, "logp": 0.27, "tpsa": 210.54, "hbd": 5, "hba": 10, "rotatable_bonds": 9, "aromatic_rings": 3, "formal_charge": 0, "is_antibody": False},
    "doxorubicin": {"mw": 543.53, "logp": 0.00, "tpsa": 206.07, "hbd": 6, "hba": 12, "rotatable_bonds": 5, "aromatic_rings": 2, "formal_charge": 0, "is_antibody": False},
    "cisplatin": {"mw": 298.03, "logp": 0.20, "tpsa": 52.04, "hbd": 2, "hba": 2, "rotatable_bonds": 0, "aromatic_rings": 0, "formal_charge": 0, "is_antibody": False},
    "etoposide": {"mw": 588.56, "logp": 1.34, "tpsa": 160.83, "hbd": 3, "hba": 13, "rotatable_bonds": 5, "aromatic_rings": 2, "formal_charge": 0, "is_antibody": False},
    "vincristine": {"mw": 824.97, "logp": 3.52, "tpsa": 171.17, "hbd": 3, "hba": 12, "rotatable_bonds": 8, "aromatic_rings": 3, "formal_charge": 0, "is_antibody": False},
    "irinotecan": {"mw": 586.689, "logp": 4.091, "tpsa": 114.2, "hbd": 1, "hba": 9, "rotatable_bonds": 4, "aromatic_rings": 3, "formal_charge": 0, "is_antibody": False},
    # Macromolecule marker: therapeutic IgG monoclonal antibody (~149 kDa).
    # No small-molecule SMILES; descriptors are population/known values.
    "bevacizumab": {"mw": 149000.0, "logp": -1.0, "tpsa": 0.0, "hbd": 0, "hba": 0, "rotatable_bonds": 0, "aromatic_rings": 0, "formal_charge": 0, "is_antibody": True},
}


def compute_descriptors(smiles: str) -> Dict[str, float]:
    """RDKit -> physicochemical descriptor dict.

    Returns keys: mw, logp, tpsa, hbd, hba, rotatable_bonds, aromatic_rings,
    formal_charge, is_antibody (always False here; small molecules only).

    Raises ValueError on unparseable SMILES.
    """
    if not _HAVE_RDKIT:
        raise ValueError("RDKit unavailable; cannot compute descriptors from SMILES")
    if not smiles or not isinstance(smiles, str):
        raise ValueError(f"empty/invalid SMILES: {smiles!r}")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"unparseable SMILES: {smiles!r}")
    return {
        "mw": round(Descriptors.MolWt(mol), 3),
        "logp": round(Crippen.MolLogP(mol), 3),
        "tpsa": round(rdMolDescriptors.CalcTPSA(mol), 3),
        "hbd": int(Lipinski.NumHDonors(mol)),
        "hba": int(Lipinski.NumHAcceptors(mol)),
        "rotatable_bonds": int(rdMolDescriptors.CalcNumRotatableBonds(mol)),
        "aromatic_rings": int(rdMolDescriptors.CalcNumAromaticRings(mol)),
        "formal_charge": int(Chem.GetFormalCharge(mol)),
        "is_antibody": False,
    }


def from_smiles(name: str, smiles: str) -> Molecule:
    """Build a Molecule from name + SMILES.

    Computes RDKit descriptors and additionally derives radius_hydrodynamic_nm
    from mw via `hydrodynamic_radius_nm`. If RDKit parsing fails (or no valid
    SMILES exists, e.g. an antibody marker), falls back to the curated
    descriptor table keyed by name and sets source='fallback_table'.

    Never raises: an unknown-name parse failure yields a minimal, honest
    Molecule rather than crashing the pipeline.
    """
    source = "rdkit"
    descriptors: Dict[str, float]
    try:
        descriptors = compute_descriptors(smiles)
    except ValueError:
        key = (name or "").strip().lower()
        if key in _FALLBACK_TABLE:
            descriptors = dict(_FALLBACK_TABLE[key])
            source = "fallback_table"
        else:
            # Nothing parseable and nothing curated: emit an honest, minimal
            # record instead of throwing. Radius/physics downstream will be
            # degenerate but the pipeline stays alive.
            descriptors = {
                "mw": 0.0, "logp": 0.0, "tpsa": 0.0, "hbd": 0, "hba": 0,
                "rotatable_bonds": 0, "aromatic_rings": 0, "formal_charge": 0,
                "is_antibody": False,
            }
            source = "fallback_table"

    mw = float(descriptors.get("mw", 0.0) or 0.0)
    is_ab = bool(descriptors.get("is_antibody", False))
    descriptors["radius_hydrodynamic_nm"] = round(hydrodynamic_radius_nm(mw, is_ab), 4)

    return Molecule(name=name, smiles=smiles, descriptors=descriptors, source=source)


if __name__ == "__main__":
    # Quick smoke test of the descriptor path.
    for nm, sm in [("temozolomide", "CN1C(=O)N2C=NC(=C2N=N1)C(=O)N"),
                   ("doxorubicin", "COc1cccc2c1C(=O)c1c(O)c3c(c(O)c1C2=O)CC(O)(C(=O)CO)CC3OC1CC(N)C(O)C(C)O1")]:
        m = from_smiles(nm, sm)
        print(nm, "->", {k: m.get(k) for k in ("mw", "logp", "tpsa", "radius_hydrodynamic_nm")}, "src", m.source)
    # Fallback path (bad SMILES, name in table):
    m = from_smiles("temozolomide", "not-a-smiles")
    assert m.source == "fallback_table", m.source
    print("fallback ok:", m.get("mw"), m.source)
    try:
        compute_descriptors("Xy!!zzz")
        raise SystemExit("expected ValueError")
    except ValueError:
        print("compute_descriptors raises on bad SMILES: ok")
