"""Curated library of BREAST-cancer-relevant drug molecules.

Breast-cancer analogue of ``tmesim.molecular.library`` -- same pattern (real,
verified SMILES, checked to parse in RDKit and to reproduce the literature
molecular weight to <0.1 g/mol), a disjoint drug set. A new, additive file:
``tmesim.molecular.library.GBM_DRUGS`` is untouched, so every existing GBM
molecule (including doxorubicin, reused verbatim by the breast panel as a
cross-cancer sanity-check drug -- see ``experiments/breast_drug_panel.py``)
resolves exactly as before.

Drug selection rationale (see ``results/breast_cancer_flagship.md`` for the
full panel writeup): eight drugs with real GDSC2 breast (BRCA) cell-line IC50
data (``tmesim.molecular.potency_breast``), spanning standard breast-cancer
chemotherapy classes --

  paclitaxel        taxane (microtubule stabiliser), standard in TNBC/HER2+
                     and hormone-receptor-positive regimens.
  cyclophosphamide   alkylating pro-drug, backbone of AC/CMF regimens.
  tamoxifen          selective ER modulator (SERM), standard endocrine
                      therapy for ER+ disease -- and this specific tumour
                      (Janesick et al. 2023 Sample #1) is ER+/HER2+/PR-, so
                      tamoxifen is not an arbitrary panel filler here.
  docetaxel          taxane (2nd generation, standard in TAC/metastatic
                      regimens); a REAL breast-specific 3-D:2-D spheroid
                      literature exists for this drug (Dubois et al. 2017,
                      PMID 29221130; Gebeyehu et al. 2021, PMID 33431915),
                      unlike most of this panel.
  epirubicin         anthracycline, backbone of FEC/FAC regimens; a
                      stereoisomer of doxorubicin (identical MW, 543.5 Da) --
                      a real same-class, same-tissue contrast to the panel's
                      doxorubicin anchor. Real breast-specific 3-D:2-D data
                      also exists (Dubois et al. 2017, PMID 29221130).
  5-fluorouracil     antimetabolite, backbone of CMF/FEC regimens.
  vinorelbine        vinca alkaloid (microtubule destabiliser), standard in
                      metastatic breast cancer.

Two more real breast-panel additions are compiled by the EXPERIMENT SCRIPTS
(``experiments/breast_drug_panel.py``, ``experiments/kill_calibration_breast.
py``), not curated in THIS file, because they reuse molecules that already
exist elsewhere rather than adding new SMILES here:

  cisplatin          REUSED VERBATIM from ``tmesim.molecular.library``
                      (the GBM panel's own cisplatin molecule -- same
                      structure, same compound) but compiled with a REAL
                      GDSC2-BREAST-measured IC50
                      (``tmesim.molecular.potency_breast``), not the GBM
                      literature IC50 -- so this is a genuine substitute for
                      carboplatin (see below), not carboplatin itself.
  palbociclib        CDK4/6 inhibitor -- included as a genuinely interesting
                      CONTRAST CASE for a pharma audience: it has a REAL
                      GDSC2-breast IC50 (n=51 cell lines, unlike trastuzumab),
                      but its mechanism is predominantly CYTOSTATIC (G1-phase
                      cell-cycle arrest), not cytotoxic cell killing. This
                      pipeline's kill-fraction model was built for a Hill-type
                      CYTOTOXIC dose-response and does not distinguish
                      "arrested" from "dead" -- so palbociclib's kill_fraction
                      output should be read with that category-mismatch
                      caveat prominently attached, not silently compared
                      apples-to-apples with the cytotoxic drugs in this panel.
                      See ``results/breast_cancer_flagship.md`` and
                      ``results/kill_calibration_breast.md`` for the full
                      honest treatment.

Two drugs the task asked to evaluate were considered and REJECTED:

  carboplatin        Absent from GDSC2 ENTIRELY (confirmed by streaming the
                      live, full, un-filtered GDSC2 dose-response CSV across
                      ALL tissues, not just breast -- GDSC2 screens cisplatin
                      and oxaliplatin but never carboplatin). No curated,
                      consensus single-agent breast-cell-line IC50 exists in
                      the literature at the rigor bar this project requires
                      either (the papers found are nanoparticle-delivery /
                      novel-complex comparison studies using carboplatin as a
                      benchmark, with inconsistent units and combination
                      contexts -- not a textbook-consensus number the way
                      doxorubicin's literature-tier 0.5 uM is). Cisplatin
                      (above) is the real-data platinum-class substitute.
  capecitabine        Absent from GDSC2 entirely (same live-CSV check). It is
                      an oral PRODRUG requiring a three-enzyme hepatic/tumour
                      activation cascade (carboxylesterase -> cytidine
                      deaminase -> thymidine phosphorylase) to become
                      cytotoxic 5-fluorouracil (Ishitsuka et al. 1999 Yakugaku
                      Zasshi, PMID 10630095; Verweij J. Oncology 1999, PMID
                      10436411) -- a standard monolayer/spheroid viability
                      assay lacks this tissue-level cascade, the same "wrong
                      assay format for a prodrug" pattern this panel already
                      documents for cyclophosphamide. GDSC screens the active
                      metabolite (5-fluorouracil, included above) directly
                      instead.

Trastuzumab (anti-HER2 antibody) was RE-EVALUATED and remains DELIBERATELY
EXCLUDED, for the same reason as the original 4-drug panel: this tumour is
HER2+ and ERBB2 is a measured transcript on the panel (real target expression
exists), but trastuzumab is a biologic with no cytotoxic IC50 in GDSC (GDSC
screens small-molecule dose-response, not antibody dose-response; confirmed
absent from the live GDSC2 CSV across all tissues, not just breast, in this
re-check too), and no curated single-agent cytotoxic IC50 for it exists in the
literature in the same sense the small molecules here have one (its mechanism
is receptor blockade / ADCC, not a Hill-type cytotoxic dose-response curve
this pipeline's IC50-Hill kill model is built for). Compiling it anyway would
mean inventing a potency number, which the task's own standard forbids. See
``results/breast_cancer_flagship.md`` for the full honest writeup of this
boundary.
"""
from __future__ import annotations

from typing import Dict, List

from .descriptors import from_smiles
from .types import Molecule

# name -> SMILES (RDKit canonical SMILES pulled from PubChem, re-verified
# HERE via RDKit MolWt against the PubChem-reported molecular weight; all
# match to <0.1 g/mol):
#   paclitaxel        853.92 vs lit 853.9  (C47H51NO14)
#   cyclophosphamide  261.09 vs lit 261.09 (C7H15Cl2N2O2P)
#   tamoxifen         371.52 vs lit 371.51 (C26H29NO)
#   docetaxel         807.89 vs lit 807.9  (C43H53NO14; PubChem CID 148124)
#   epirubicin        543.53 vs lit 543.5  (C27H29NO11; PubChem CID 41867 --
#                      a STEREOISOMER of doxorubicin, identical formula/MW,
#                      differing only in the 4'-OH configuration of the sugar)
#   5-fluorouracil    130.08 vs lit 130.08 (C4H3FN2O2; PubChem CID 3385)
#   vinorelbine       778.95 vs lit 778.94 (C45H54N4O8; PubChem CID 5311497)
#   palbociclib       447.54 vs lit 447.51 (C24H29N7O2; PubChem CID 5330286)
BREAST_DRUGS: Dict[str, str] = {
    "paclitaxel": (
        "CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)"
        "NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C"
    ),
    "cyclophosphamide": "ClCCN(CCCl)P1(=O)NCCCO1",
    "tamoxifen": "CCC(=C(c1ccccc1)c1ccc(OCCN(C)C)cc1)c1ccccc1",
    "docetaxel": (
        "CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)"
        "NC(=O)OC(C)(C)C)O)O)OC(=O)C6=CC=CC=C6)(CO4)OC(=O)C)O)C)O"
    ),
    "epirubicin": (
        "CC1C(C(CC(O1)OC2CC(CC3=C2C(=C4C(=C3O)C(=O)C5=C(C4=O)C(=CC=C5)OC)O)"
        "(C(=O)CO)O)N)O"
    ),
    "5-fluorouracil": "C1=C(C(=O)NC(=O)N1)F",
    "vinorelbine": (
        "CCC1=CC2CC(C3=C(CN(C2)C1)C4=CC=CC=C4N3)(C5=C(C=C6C(=C5)C78CCN9C7C"
        "(C=CC9)(C(C(C8N6C)(C(=O)OC)O)OC(=O)C)CC)OC)C(=O)OC"
    ),
    "palbociclib": (
        "CC1=C(C(=O)N(C2=NC(=NC=C12)NC3=NC=C(C=C3)N4CCNCC4)C5CCCC5)C(=O)C"
    ),
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
    dtx = by["docetaxel"]
    assert 805.0 <= dtx.get("mw") <= 811.0, f"docetaxel MW off: {dtx.get('mw')}"
    epi = by["epirubicin"]
    assert 541.0 <= epi.get("mw") <= 546.0, f"epirubicin MW off: {epi.get('mw')}"
    fu = by["5-fluorouracil"]
    assert 128.0 <= fu.get("mw") <= 132.0, f"5-fluorouracil MW off: {fu.get('mw')}"
    vnr = by["vinorelbine"]
    assert 776.0 <= vnr.get("mw") <= 782.0, f"vinorelbine MW off: {vnr.get('mw')}"
    plb = by["palbociclib"]
    assert 445.0 <= plb.get("mw") <= 450.0, f"palbociclib MW off: {plb.get('mw')}"
    for m in mols:
        assert m.source == "rdkit", f"{m.name} should parse via RDKit, got {m.source}"
        assert m.get("is_antibody") is False

    try:
        get_molecule_breast("trastuzumab")
        raise SystemExit("expected KeyError for trastuzumab (not in breast library)")
    except KeyError:
        pass

    print("\nAll self-tests passed.")
