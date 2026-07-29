"""Curated library of GBM-relevant drug molecules.

Real, verified SMILES for the small-molecule agents used against glioblastoma,
plus one macromolecule marker (bevacizumab, a ~149 kDa anti-VEGF IgG) to expose
the size extreme that dominates transport. Every SMILES here has been checked to
parse in RDKit and to reproduce the literature molecular weight.

The antibody has NO small-molecule SMILES: it is handled specially so its
descriptors come from known protein values (is_antibody=True, mw ~149 kDa) with
a protein hydrodynamic radius, never a (doomed) SMILES parse.
"""
from __future__ import annotations

from typing import Dict, List

from .descriptors import from_smiles
from .types import Molecule

# name -> SMILES (verified via RDKit; MW matches literature to <1 g/mol for
# small molecules). See module docstring for the antibody exception.
GBM_DRUGS: Dict[str, str] = {
    # Alkylating / methylating agents (first-line GBM chemotherapy).
    "temozolomide": "CN1C(=O)N2C=NC(=C2N=N1)C(=O)N",              # MW ~194
    "lomustine":    "O=C(NC1CCCCC1)N(CCCl)N=O",                    # CCNU, MW ~234
    "carmustine":   "O=C(NCCCl)N(CCCl)N=O",                        # BCNU, MW ~214
    "procarbazine": "CNNCc1ccc(cc1)C(=O)NC(C)C",                   # MW ~221
    # Antimetabolite.
    "methotrexate": "CN(Cc1cnc2nc(N)nc(N)c2n1)c1ccc(cc1)C(=O)NC(CCC(=O)O)C(=O)O",  # MW ~454
    # Anthracycline (intercalator).
    "doxorubicin":  "COc1cccc2c1C(=O)c1c(O)c3c(c(O)c1C2=O)CC(O)(C(=O)CO)CC3OC1CC(N)C(O)C(C)O1",  # MW ~544
    # Platinum crosslinker (contains Pt; charge-neutral square-planar complex).
    "cisplatin":    "N[Pt](N)(Cl)Cl",                             # MW ~298
    # Topoisomerase-II poison.
    "etoposide":    "COc1cc(cc(OC)c1O)C1c2cc3OCOc3cc2C(OC2OC4COC(C)OC4C(O)C2O)C2COC(=O)C12",  # MW ~589
    # Vinca alkaloid (microtubule inhibitor); PubChem canonical SMILES.
    "vincristine":  "CCC1(CC2CC(C3=C(CCN(C2)C1)C4=CC=CC=C4N3)(C5=C(C=C6C(=C5)C78CCN9C7C(C=CC9)(C(C(C8N6C=O)(C(=O)OC)O)OC(=O)C)CC)OC)C(=O)OC)O",  # MW ~825
    # Topoisomerase-I inhibitor prodrug (camptothecin derivative); PubChem CID
    # 60838 canonical SMILES. MW ~586.7 (verified in RDKit, matches the
    # literature free-base MW of 586.68). NOTE (honesty): irinotecan itself is
    # a prodrug -- it is cleaved by tissue/serum carboxylesterases to the ~1000x
    # more potent topoisomerase-I poison SN-38. This pipeline has no prodrug
    # metabolism step, so the compiled DrugParams describe the PARENT molecule
    # (irinotecan) only; both its GDSC2 potency (cell assays dose irinotecan
    # directly, so intracellular conversion is baked into the measured IC50) and
    # its structure-derived BBB penetrance are for irinotecan, not SN-38. See
    # ``results/real_drug_panel.md`` for the honesty note this produces.
    "irinotecan":   "CCC1=C2CN3C(=CC4=C(C3=O)COC(=O)C4(CC)O)C2=NC5=C1C=C(C=C5)OC(=O)N6CCC(CC6)N7CCCCC7",  # MW ~587
    # Macromolecule marker: anti-VEGF monoclonal antibody (~149 kDa).
    # Sentinel string, NOT a parseable SMILES -> forces the antibody path.
    "bevacizumab":  "ANTIBODY:bevacizumab",
}

# Names whose "SMILES" are markers; resolved via the descriptors fallback table.
_ANTIBODIES = {"bevacizumab"}


def get_molecule(name: str) -> Molecule:
    """Return a Molecule for a library drug (built via descriptors.from_smiles).

    Antibody markers resolve through from_smiles' fallback table (the sentinel
    string does not parse), yielding is_antibody=True and a protein radius.

    Raises KeyError (listing available names) for unknown drugs.
    """
    key = name.strip().lower()
    if key not in GBM_DRUGS:
        available = ", ".join(sorted(GBM_DRUGS))
        raise KeyError(f"unknown drug {name!r}; available: {available}")
    return from_smiles(key, GBM_DRUGS[key])


def all_molecules() -> List[Molecule]:
    """Build and return every molecule in the library."""
    return [get_molecule(n) for n in GBM_DRUGS]


if __name__ == "__main__":
    mols = all_molecules()

    cols = ("mw", "logp", "tpsa", "hbd", "hba", "radius_hydrodynamic_nm", "is_antibody")
    header = f"{'name':14s} {'src':14s} " + " ".join(f"{c:>10s}" for c in cols)
    print(header)
    print("-" * len(header))
    for m in mols:
        row = f"{m.name:14s} {m.source:14s} "
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

    # --- assertions -------------------------------------------------------
    tmz = by["temozolomide"]
    assert 193.0 <= tmz.get("mw") <= 195.0, f"TMZ MW off: {tmz.get('mw')}"

    dox = by["doxorubicin"]
    assert dox.get("mw") > 500.0, f"doxorubicin MW should exceed 500: {dox.get('mw')}"

    cis = by["cisplatin"]
    assert "Pt" in cis.smiles, "cisplatin SMILES must contain Pt"
    assert cis.get("mw") < 320.0, f"cisplatin should be small: {cis.get('mw')}"

    bev = by["bevacizumab"]
    assert bev.get("is_antibody") is True, "bevacizumab must be flagged is_antibody"
    assert bev.get("mw") > 100_000.0, f"antibody MW should be ~149 kDa: {bev.get('mw')}"
    assert bev.source == "fallback_table", f"antibody should use fallback: {bev.source}"
    # Antibody radius (~5.5 nm) must dwarf a small molecule (~0.4 nm).
    assert bev.get("radius_hydrodynamic_nm") > 4.0, bev.get("radius_hydrodynamic_nm")
    assert tmz.get("radius_hydrodynamic_nm") < 1.0, tmz.get("radius_hydrodynamic_nm")

    # KeyError-safe unknown lookup.
    try:
        get_molecule("aspirin_not_here")
        raise SystemExit("expected KeyError for unknown drug")
    except KeyError:
        pass

    print("\nAll self-tests passed.")
