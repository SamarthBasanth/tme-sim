"""Curated library of BREAST-cancer-relevant drug molecules.

Breast-cancer analogue of ``tmesim.molecular.library`` -- same pattern (real,
verified SMILES, checked to parse in RDKit and to reproduce the literature
molecular weight to <0.1 g/mol), a disjoint drug set. A new, additive file:
``tmesim.molecular.library.GBM_DRUGS`` is untouched, so every existing GBM
molecule (including doxorubicin, reused verbatim by the breast panel as a
cross-cancer sanity-check drug -- see ``experiments/breast_drug_panel.py``)
resolves exactly as before.

Drug selection rationale (see ``results/breast_cancer_flagship.md`` for the
full panel writeup): three drugs with real GDSC2 breast (BRCA) cell-line IC50
data (``tmesim.molecular.potency_breast``), spanning standard breast-cancer
chemotherapy classes --

  paclitaxel        taxane (microtubule stabiliser), standard in TNBC/HER2+
                     and hormone-receptor-positive regimens.
  cyclophosphamide   alkylating pro-drug, backbone of AC/CMF regimens.
  tamoxifen          selective ER modulator (SERM), standard endocrine
                      therapy for ER+ disease -- and this specific tumour
                      (Janesick et al. 2023 Sample #1) is ER+/HER2+/PR-, so
                      tamoxifen is not an arbitrary panel filler here.

Trastuzumab (anti-HER2 antibody) was evaluated and DELIBERATELY EXCLUDED: this
tumour is HER2+ and ERBB2 is a measured transcript on the panel (real target
expression exists), but trastuzumab is a biologic with no cytotoxic IC50 in
GDSC (GDSC screens small-molecule dose-response, not antibody dose-response),
and no curated single-agent cytotoxic IC50 for it exists in the literature in
the same sense the small molecules here have one (its mechanism is
receptor blockade / ADCC, not a Hill-type cytotoxic dose-response curve this
pipeline's IC50-Hill kill model is built for). Compiling it anyway would mean
inventing a potency number, which the task's own standard forbids. See
``results/breast_cancer_flagship.md`` for the full honest writeup of this
boundary.
"""
from __future__ import annotations

from typing import Dict, List

from .descriptors import from_smiles
from .types import Molecule

# name -> SMILES (verified via RDKit; MW matches literature to <0.1 g/mol,
# checked directly against PubChem/DrugBank reference values):
#   paclitaxel        853.92 vs lit 853.9  (C47H51NO14)
#   cyclophosphamide  261.09 vs lit 261.09 (C7H15Cl2N2O2P)
#   tamoxifen         371.52 vs lit 371.51 (C26H29NO)
BREAST_DRUGS: Dict[str, str] = {
    "paclitaxel": (
        "CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)"
        "NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C"
    ),
    "cyclophosphamide": "ClCCN(CCCl)P1(=O)NCCCO1",
    "tamoxifen": "CCC(=C(c1ccccc1)c1ccc(OCCN(C)C)cc1)c1ccccc1",
}


def get_molecule_breast(name: str) -> Molecule:
    """Return a Molecule for a breast-library drug (via descriptors.from_smiles).

    Raises KeyError (listing available names) for unknown drugs. Does NOT fall
    through to ``tmesim.molecular.library.GBM_DRUGS`` -- the two libraries are
    kept disjoint on purpose so a typo can't silently resolve a breast-panel
    call to a GBM molecule (or vice versa).
    """
    key = name.strip().lower()
    if key not in BREAST_DRUGS:
        available = ", ".join(sorted(BREAST_DRUGS))
        raise KeyError(f"unknown breast-panel drug {name!r}; available: {available}")
    return from_smiles(key, BREAST_DRUGS[key])


def all_molecules_breast() -> List[Molecule]:
    """Build and return every molecule in the breast library."""
    return [get_molecule_breast(n) for n in BREAST_DRUGS]


if __name__ == "__main__":
    mols = all_molecules_breast()
    cols = ("mw", "logp", "tpsa", "hbd", "hba", "radius_hydrodynamic_nm", "is_antibody")
    header = f"{'name':18s} {'src':14s} " + " ".join(f"{c:>10s}" for c in cols)
    print(header)
    print("-" * len(header))
    for m in mols:
        row = f"{m.name:18s} {m.source:14s} "
        cells = []
        for c in cols:
            v = m.get(c)
            if isinstance(v, bool):
                cells.append(f"{str(v):>10s}")
            elif isinstance(v, (int, float)):
                cells.append(f"{v:>10.2f}")
            else:
                cells.append(f"{str(v):>10s}")
        print(row + " ".join(cells))

    by = {m.name: m for m in mols}
    ptx = by["paclitaxel"]
    assert 850.0 <= ptx.get("mw") <= 858.0, f"paclitaxel MW off: {ptx.get('mw')}"
    cyc = by["cyclophosphamide"]
    assert 259.0 <= cyc.get("mw") <= 263.0, f"cyclophosphamide MW off: {cyc.get('mw')}"
    tam = by["tamoxifen"]
    assert 369.0 <= tam.get("mw") <= 374.0, f"tamoxifen MW off: {tam.get('mw')}"
    for m in mols:
        assert m.source == "rdkit", f"{m.name} should parse via RDKit, got {m.source}"
        assert m.get("is_antibody") is False

    try:
        get_molecule_breast("trastuzumab")
        raise SystemExit("expected KeyError for trastuzumab (not in breast library)")
    except KeyError:
        pass

    print("\nAll self-tests passed.")
