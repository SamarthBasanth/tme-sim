"""Panel-aware marker -> :class:`CellType` maps for the real Xenium sections.

Why this module exists
----------------------
``reconstruct.real_xenium`` phenotypes every measured cell from a hardcoded
marker map that was written for ONE section: the 10x Xenium **Human Brain**
panel (hBrain_254g + 65-gene add-on = 319 targets). That map is built on
neural-lineage and vessel markers — ``SOX2``, ``OLIG2``, ``AQP4``, ``GJA1``,
``MOG``, ``PECAM1``.

A second real GBM section is now cached (``data_cache/xenium_gbm_io``): 10x
**Xenium Human Brain Cancer FFPE, Immuno-Oncology panel** (hImmune_v1 380 genes
+ Add-on B 100 genes = 480 targets, 816,769 cells). It carries a far deeper
immune / vascular / hypoxia readout — crucially **VEGFA** and **HIF1A**, which
the brain panel does not measure at all — but it carries NONE of
``PECAM1 / SOX2 / OLIG2 / AQP4 / GJA1 / MOG / TOP2A``.

Applied blindly to that section the brain map degenerates: only ``EGFR`` and
``MKI67`` survive for TUMOR, the whole ENDOTHELIAL and STROMAL classes vanish
(zero markers present), and the reconstruction collapses to **TUMOR 93.1% /
IMMUNE 6.9%** with no vasculature at all. Since ``vascular_density`` (the
oxygen SOURCE) is a KDE of the ENDOTHELIAL cells, that failure silently
produces an oxygen field with no supply term.

So the map must follow the PANEL, not the code path. :func:`select_marker_map`
picks the map whose markers the loaded panel actually measures.

Backward compatibility (load-bearing)
-------------------------------------
:data:`BRAIN_PANEL_MARKER_MAP` is the byte-for-byte content, and the same key
ORDER, of the map that used to live in ``real_xenium.MARKER_MAP`` (key order
matters: it sets the arg-max tie-break in ``phenotype_cells``).
``real_xenium.MARKER_MAP`` is now an alias for it, and :func:`select_marker_map`
returns *that same dict object* for the 319-gene brain panel. The 319-panel
reconstruction is therefore unchanged: TUMOR 73.80 / STROMAL 7.92 / IMMUNE 8.54
/ ENDOTHELIAL 9.74 percent, verified against the pre-change numbers.

The IO-panel map, and why these genes
-------------------------------------
Each entry is the closest measurable analogue of the brain-panel class:

  TUMOR       <- EGFR, PDGFRA, IL13RA2, MET
      Glioma driver receptors. ``EGFR`` is shared with the brain map; ``PDGFRA``
      is the proneural-subtype driver; ``IL13RA2`` and ``MET`` are
      glioma-restricted (mesenchymal) markers. Replaces the unmeasurable
      ``SOX2``/``OLIG2``. TUMOR is also the ungated default class, so this map
      determines only which cells are *positively* called tumour.
  IMMUNE      <- PTPRC, CD68, AIF1, CD163
      IDENTICAL to the brain map — all four are on both panels (this is the IO
      panel's strong suit).
  ENDOTHELIAL <- VWF, FLT1, PLVAP
      ``PECAM1`` is absent. ``VWF`` (von Willebrand factor, Weibel-Palade
      bodies), ``FLT1`` (VEGFR1) and ``PLVAP`` (plasmalemma vesicle-associated
      protein, the fenestral diaphragm protein strongly induced in GBM's leaky
      neovasculature) are all canonical endothelial genes. Three markers rather
      than the brain map's single ``PECAM1``, which also removes the
      single-vs-multi-marker imbalance discussed in ``real_xenium``.
  STROMAL     <- SLC1A3, S100B, PTN
      The direct analogue of the brain map's astro-lineage ``AQP4``/``GJA1``
      (+ oligodendrocyte ``MOG``): astrocytic glutamate transporter, the
      astro/glial calcium-binding protein, and pleiotrophin (glial/perivascular
      niche matrix).

NECROTIC is still never assigned: neither panel carries a necrosis marker.

MURAL / PERICYTE — measured, evaluated, and deliberately NOT used
----------------------------------------------------------------
The IO panel does measure the pericyte/smooth-muscle programme (``RGS5`` 12.5%,
``PDGFRB`` 9.5%, ``NOTCH3`` 24.9%, ``ACTA2`` 38.0% of cells). ``CellType`` has
no MURAL member, so the only options were to fold those genes into STROMAL or
drop them. Folding them in was measured and rejected: ``phenotype_cells`` scores
a type by the MEAN z-score of its markers, so averaging two disjoint programmes
(astro + mural) dilutes BOTH, and a cell strongly positive for one is dragged
below the gate by the other. Measured cell-type distributions on the IO section:

    STROMAL markers                            TUMOR  STROM  IMMUNE  ENDO
    SLC1A3,S100B,PTN            (chosen)       79.09   9.98    6.74  4.18
    RGS5,PDGFRB,NOTCH3,SLC1A3,S100B,PTN        86.57   2.68    6.85  3.90
    RGS5,PDGFRB,NOTCH3                         86.18   4.02    6.76  3.04

i.e. the union collapses STROMAL to 2.7% — worse than either programme alone.
The astro-lineage programme was chosen because it is the like-for-like analogue
of the brain-panel STROMAL class, keeping the two sections comparable. HONEST
CONSEQUENCE: pericytes on the IO section are NOT separately labelled; they fall
into whichever class wins, usually the ungated TUMOR default. Fixing this
properly needs a MURAL member on ``CellType``, which is a core-schema change and
is out of scope here.

CAVEAT on the astro-lineage programme in GBM: ``SLC1A3``/``S100B``/``PTN`` are
expressed by *malignant* astrocyte-like glioma cells as well as by reactive
astrocytes, so this STROMAL class is not a clean non-malignant compartment. The
brain panel's ``AQP4``/``GJA1`` have exactly the same ambiguity, so the two
sections are at least consistently wrong here rather than differently wrong.
"""
from __future__ import annotations

from typing import Dict, List, Sequence, Tuple

import numpy as np

from ..core.microenvironment import CellType

# ---------------------------------------------------------------------------
# Panel maps.
#
# KEY ORDER IS LOAD-BEARING: ``phenotype_cells`` builds its score matrix in this
# order and resolves an exact arg-max tie by taking the first column, so
# reordering these keys can silently change cell-type assignments. Do not
# reorder BRAIN_PANEL_MARKER_MAP.
# ---------------------------------------------------------------------------

#: The ORIGINAL 319-gene Xenium Human Brain panel map. Byte-for-byte (and in the
#: same key order as) the map that used to be defined inline in
#: ``real_xenium.MARKER_MAP``, which is now an alias for this object.
BRAIN_PANEL_MARKER_MAP: Dict[CellType, List[str]] = {
    CellType.TUMOR: ["EGFR", "SOX2", "OLIG2", "MKI67"],
    CellType.IMMUNE: ["PTPRC", "CD68", "AIF1", "CD163"],
    CellType.ENDOTHELIAL: ["PECAM1"],
    CellType.STROMAL: ["MOG", "AQP4", "GJA1"],
}

#: The 480-gene Xenium Human Immuno-Oncology Profiling panel map (see the module
#: docstring for the per-gene rationale and for the mural/pericyte decision).
IO_PANEL_MARKER_MAP: Dict[CellType, List[str]] = {
    CellType.TUMOR: ["EGFR", "PDGFRA", "IL13RA2", "MET"],
    CellType.IMMUNE: ["PTPRC", "CD68", "AIF1", "CD163"],
    CellType.ENDOTHELIAL: ["VWF", "FLT1", "PLVAP"],
    CellType.STROMAL: ["SLC1A3", "S100B", "PTN"],
}

#: Measured hypoxia-RESPONSE genes on each panel. This is what makes the IO
#: section usable as an independent check on the SIMULATED oxygen field: the
#: brain panel measures neither gene, so the check is impossible there.
#:
#: NOT an oxygen tension. VEGFA and HIF1A are downstream transcriptional
#: responses to hypoxia (HIF1A is itself regulated mainly post-translationally,
#: so its mRNA is the weaker of the two readouts), and both are induced by
#: non-hypoxic stimuli as well. See ``tmesim/validation/hypoxia_validation.py``.
HYPOXIA_RESPONSE_GENES: Dict[str, List[str]] = {
    "brain_319": [],
    "io_480": ["VEGFA", "HIF1A"],
}

#: Registry, in PREFERENCE ORDER. :func:`select_marker_map` walks this list and
#: keeps the first map with the highest marker coverage, so the brain panel wins
#: any tie and remains the default for an unrecognised panel.
PANEL_REGISTRY: Tuple[Tuple[str, Dict[CellType, List[str]]], ...] = (
    ("brain_319", BRAIN_PANEL_MARKER_MAP),
    ("io_480", IO_PANEL_MARKER_MAP),
)


def marker_coverage(
    marker_map: Dict[CellType, List[str]], genes: Sequence[str]
) -> float:
    """Fraction of ``marker_map``'s markers that the panel ``genes`` measures.

    A map is only usable to the extent its markers were actually assayed. 1.0
    means every marker is on the panel; 0.5 means half of them are silently
    dropped (and any class whose markers are ALL missing simply disappears).
    """
    present = set(np.asarray(genes, dtype=str).tolist())
    markers = [m for ms in marker_map.values() for m in ms]
    if not markers:
        return 0.0
    return sum(m in present for m in markers) / len(markers)


def select_marker_map(
    genes: Sequence[str],
) -> Tuple[Dict[CellType, List[str]], str, Dict[str, float]]:
    """Pick the marker map best matched to the panel that was actually loaded.

    Selection is by :func:`marker_coverage` — the fraction of a map's markers
    the panel measures — with ties broken by :data:`PANEL_REGISTRY` order, so
    the 319-gene brain map wins any tie and is the fallback for an unrecognised
    panel (preserving the historical behaviour exactly).

    Measured coverages (the selection is unambiguous, not a close call):

        panel loaded        brain map   IO map    selected
        brain 319-gene        1.00        0.50    brain_319
        IO 480-gene           0.50        1.00    io_480

    Parameters
    ----------
    genes : sequence of str
        The panel gene symbols, i.e. ``load(...)["genes"]``.

    Returns
    -------
    marker_map : dict CellType -> list of str
        The selected map. For the brain panel this is the *same object* as
        :data:`BRAIN_PANEL_MARKER_MAP` (== ``real_xenium.MARKER_MAP``).
    panel_name : str
        Registry key of the selected panel (``"brain_319"`` / ``"io_480"``).
    coverage : dict str -> float
        Coverage of EVERY registered map, so a call site can report how close
        the decision was and flag a panel that fits nothing well.
    """
    coverage = {name: marker_coverage(mm, genes) for name, mm in PANEL_REGISTRY}
    best_name, best_map = PANEL_REGISTRY[0]
    best_cov = coverage[best_name]
    for name, mm in PANEL_REGISTRY[1:]:
        if coverage[name] > best_cov:
            best_name, best_map, best_cov = name, mm, coverage[name]
    return best_map, best_name, coverage


def hypoxia_genes_for(genes: Sequence[str]) -> List[str]:
    """Hypoxia-response genes measured by the panel that ``genes`` came from.

    Returns ``[]`` when the panel measures none (the 319-gene brain panel), which
    is the honest answer: on that section there is no measured hypoxia readout to
    validate the simulated oxygen field against.
    """
    _, panel_name, _ = select_marker_map(genes)
    present = set(np.asarray(genes, dtype=str).tolist())
    return [g for g in HYPOXIA_RESPONSE_GENES.get(panel_name, []) if g in present]


if __name__ == "__main__":  # pragma: no cover - smoke test on the real caches
    from ..data.real.singlecell_spatial import load

    for path in ("data_cache/xenium_gbm", "data_cache/xenium_gbm_io"):
        d = load(path=path)
        g = np.asarray(d["genes"], dtype=str)
        mm, name, cov = select_marker_map(g)
        print(f"{path}  ({len(g)} genes)")
        print(f"  selected panel : {name}   coverage {cov}")
        print(f"  is brain map   : {mm is BRAIN_PANEL_MARKER_MAP}")
        for ct, markers in mm.items():
            present = [m for m in markers if m in set(g.tolist())]
            print(f"    {ct.name:<12}: {present}")
        print(f"  hypoxia genes  : {hypoxia_genes_for(g)}")
