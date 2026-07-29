"""Real potency (IC50) resolution for the BREAST drug panel.

Breast-cancer analogue of ``tmesim.molecular.potency`` -- SAME three-tier
trust order (measured GDSC > curated literature > flagged default), SAME
median-in-log-space aggregation, only pointed at the breast-filtered GDSC2
cache (``data_cache/gdsc_breast_dose_response.csv``, built by
``tmesim.data.acquire_breast``) instead of the glioma one. A new, additive
file: nothing in ``tmesim.molecular.potency`` is modified, and no existing GBM
IC50 resolution changes.

Why a separate module rather than parameterising the existing one: the glioma
module's ``_GDSC_ALIASES`` / ``LITERATURE_IC50`` tables are drug-name keyed and
were curated specifically against the GBM library (temozolomide, lomustine,
...); breast needs its own alias/literature tables for a different drug set
(paclitaxel, cyclophosphamide, tamoxifen, ...), so a parallel module keeps both
curation efforts independently auditable rather than interleaving two tissue
contexts in one table.

GDSC coverage on the live breast (BRCA) screen, verified 2026-07-28 by
streaming the real GDSC2 CSV and counting: paclitaxel, cyclophosphamide,
tamoxifen, docetaxel, epirubicin, 5-fluorouracil, vinorelbine, palbociclib and
cisplatin are all present with 50-51 breast cell lines each (re-verified
2026-07-29 for the panel-expansion pass -- vinorelbine 50, palbociclib 51,
cisplatin 50). Doxorubicin is NOT in GDSC2 at all (neither the glioma nor the
breast screen carries it -- it is simply not one of the ~200+ drugs GDSC2
tested), which is why the breast drug panel
(``experiments/breast_drug_panel.py``) reuses the GBM panel's own
literature-tier doxorubicin IC50 (0.5 uM, ``tmesim.molecular.potency.
LITERATURE_IC50``) verbatim rather than inventing a second, breast-specific
number for the same drug and the same literature evidence.

Carboplatin and capecitabine were checked against the FULL, un-filtered live
GDSC2 CSV (every tissue, not just BRCA) and confirmed ABSENT entirely -- see
``tmesim/molecular/breast_library.py``'s module docstring for the full
rejection rationale for both.

Units: IC50 in micromolar (uM) throughout.
"""
from __future__ import annotations

import csv
import math
import os
import statistics
from typing import Dict, List, Optional

# ---------------------------------------------------------------------------
# GDSC name matching (case-insensitive), breast panel drugs.
# ---------------------------------------------------------------------------
_GDSC_ALIASES: Dict[str, List[str]] = {
    "paclitaxel": ["paclitaxel", "taxol"],
    "cyclophosphamide": ["cyclophosphamide", "endoxan", "cytoxan", "ctx"],
    "tamoxifen": ["tamoxifen", "nolvadex"],
    "docetaxel": ["docetaxel", "taxotere"],
    "epirubicin": ["epirubicin"],
    "5-fluorouracil": ["5-fluorouracil", "5-fu", "fluorouracil"],
    "vinorelbine": ["vinorelbine", "navelbine"],
    "palbociclib": ["palbociclib", "ibrance", "pd-0332991", "pd0332991"],
    # cisplatin: the molecule itself is reused VERBATIM from
    # tmesim.molecular.library (the GBM panel's own compiled cisplatin), but
    # its BREAST potency is resolved through THIS module (measured GDSC2
    # BRCA, n=50 lines) rather than the GBM literature IC50 -- see
    # experiments/breast_drug_panel.py's / kill_calibration_breast.py's
    # compile_breast_drug() special case.
    "cisplatin": ["cisplatin", "cddp", "ddp", "cis-platinum", "cisplatinum"],
}

# ---------------------------------------------------------------------------
# Tier 2: curated literature fallback (breast-relevant IC50s, uM), used ONLY
# if the GDSC cache is missing/unreachable. Each entry cites a source; used to
# degrade gracefully, not as the primary route (GDSC hits are used in the
# actual panel run since the live cache is present).
# ---------------------------------------------------------------------------
LITERATURE_IC50: Dict[str, Dict] = {
    "paclitaxel": {
        "ic50_uM": 0.01,
        "note": (
            "Taxane, microtubule stabiliser; low-nM to ~0.01 uM in sensitive "
            "breast cancer cell lines (e.g. MCF-7). Backup only -- GDSC-measured "
            "normally wins."
        ),
    },
    "cyclophosphamide": {
        "ic50_uM": 200.0,
        "note": (
            "Alkylating pro-drug requiring hepatic CYP2B6 activation to its "
            "active metabolite (4-hydroxycyclophosphamide); the parent compound "
            "is weakly cytotoxic in a standard in-vitro viability assay, so "
            "cell-line IC50s are hundreds of uM -- the same 'weak pro-drug in "
            "vitro' pattern the GBM panel already documents for procarbazine. "
            "Backup only -- GDSC-measured normally wins."
        ),
    },
    "tamoxifen": {
        "ic50_uM": 10.0,
        "note": (
            "Selective ER modulator; cytotoxic (not just cytostatic) IC50 in "
            "breast cancer cell-line panels is typically single-digit to tens "
            "of uM. Backup only -- GDSC-measured normally wins."
        ),
    },
    "docetaxel": {
        "ic50_uM": 0.05,
        "note": (
            "2nd-generation taxane; low-nM to low-uM in breast cancer cell "
            "lines (Dubois et al. 2017, PMID 29221130, report 2-D IC50 "
            "0.2-1.8 uM in MDA-MB-231 monolayer). Backup only -- GDSC-measured "
            "normally wins."
        ),
    },
    "epirubicin": {
        "ic50_uM": 0.3,
        "note": (
            "Anthracycline, a stereoisomer of doxorubicin (same MW/formula); "
            "breast cell-line 2-D IC50 in the same sub-uM range as doxorubicin "
            "(Dubois et al. 2017, PMID 29221130, report 0.1-0.6 uM in "
            "MDA-MB-231/SUM1315 monolayer). Backup only -- GDSC-measured "
            "normally wins."
        ),
    },
    "5-fluorouracil": {
        "ic50_uM": 20.0,
        "note": (
            "Antimetabolite (thymidylate synthase inhibitor); typically "
            "single-digit to tens of uM in monolayer cell-line assays (Koch "
            "et al. 2021, PMID 33395433, report 7.5-75 uM across four "
            "colorectal lines). Backup only -- GDSC-measured normally wins."
        ),
    },
    "vinorelbine": {
        "ic50_uM": 0.05,
        "note": (
            "Vinca alkaloid (microtubule destabiliser); low-nM to low-uM in "
            "sensitive cell lines, similar order to other tubulin-targeting "
            "agents in this panel (paclitaxel, docetaxel). Backup only -- "
            "GDSC-measured normally wins."
        ),
    },
    "palbociclib": {
        "ic50_uM": 1.0,
        "note": (
            "CDK4/6 inhibitor; predominantly CYTOSTATIC (G1-arrest), not "
            "cytotoxic -- see this module's and breast_library.py's own "
            "category-mismatch caveat. Cell-viability-assay IC50 in "
            "RB-proficient breast lines typically sub-uM to low-uM. Backup "
            "only -- GDSC-measured normally wins."
        ),
    },
    "cisplatin": {
        "ic50_uM": 10.0,
        "note": (
            "Platinum crosslinker; molecule REUSED VERBATIM from "
            "tmesim.molecular.library (identical to the GBM panel's "
            "cisplatin). Backup value matches the GBM module's own cisplatin "
            "literature backup (tmesim.molecular.potency.LITERATURE_IC50) for "
            "consistency; breast cell-line IC50 typically single-digit to "
            "tens of uM. Backup only -- GDSC-measured (BRCA, n=50) normally "
            "wins."
        ),
    },
}

_CACHE_FILENAME = "gdsc_breast_dose_response.csv"
_DEFAULT_IC50_UM = 25.0
_DEFAULT_NOTE = (
    f"flagged default {_DEFAULT_IC50_UM} uM; confidence LOW -- no measured (GDSC) "
    "or curated-literature IC50 for this drug on the breast panel. Potency is "
    "NOT derivable from structure without docking/assay. Replace with a "
    "measured IC50."
)


def _gdsc_path(cache_dir: str) -> str:
    return os.path.join(cache_dir, _CACHE_FILENAME)


def _aliases_for(drug_name: str) -> List[str]:
    key = drug_name.strip().lower()
    aliases = set(_GDSC_ALIASES.get(key, []))
    aliases.add(key)
    return [a.lower() for a in aliases]


def measured_ic50(drug_name: str, cache_dir: str = "data_cache") -> Optional[Dict]:
    """Look up a drug's MEASURED IC50 from the cached GDSC BREAST (BRCA) screen.

    Identical aggregation rule to ``tmesim.molecular.potency.measured_ic50``:
    per-cell-line median LN_IC50 (collapsing repeat screens), then the
    across-cell-line median, converted back to uM with ``exp``.
    """
    path = _gdsc_path(cache_dir)
    if not os.path.exists(path):
        return None

    aliases = set(_aliases_for(drug_name))

    per_line: Dict[str, List[float]] = {}
    matched_names: set = set()
    try:
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                gdsc_name = (row.get("DRUG_NAME") or "").strip()
                if gdsc_name.lower() not in aliases:
                    continue
                line = row.get("CELL_LINE_NAME") or ""
                if not line:
                    continue
                try:
                    ln = float(row["LN_IC50"])
                except (KeyError, ValueError, TypeError):
                    continue
                if not math.isfinite(ln):
                    continue
                per_line.setdefault(line, []).append(ln)
                matched_names.add(gdsc_name)
    except Exception:  # noqa: BLE001 - degrade gracefully; caller falls back
        return None

    if not per_line:
        return None

    per_cell = [statistics.median(vs) for vs in per_line.values()]
    median_ln = float(statistics.median(per_cell))
    ic50_uM = float(math.exp(median_ln))

    return {
        "ic50_uM": ic50_uM,
        "n_cell_lines": len(per_line),
        "source": "GDSC",
        "confidence": "measured",
        "raw_median_ln_ic50": median_ln,
        "matched_gdsc_name": "|".join(sorted(matched_names)),
    }


def resolve_ic50(drug_name: str, cache_dir: str = "data_cache") -> Dict:
    """Resolve an IC50 (uM) for a BREAST-panel drug: measured GDSC (BRCA) ->
    curated literature -> flagged default. Same contract as
    ``tmesim.molecular.potency.resolve_ic50``."""
    m = measured_ic50(drug_name, cache_dir=cache_dir)
    if m is not None:
        return {
            "ic50_uM": m["ic50_uM"],
            "source": m["source"],
            "confidence": m["confidence"],
            "note": (
                f"median across {m['n_cell_lines']} breast (BRCA) cell lines "
                f"(GDSC2 {m['matched_gdsc_name']}); raw median LN_IC50="
                f"{m['raw_median_ln_ic50']:.4f}; ic50_uM=exp(median LN_IC50)."
            ),
            "n_cell_lines": m["n_cell_lines"],
            "raw_median_ln_ic50": m["raw_median_ln_ic50"],
        }

    key = drug_name.strip().lower()
    lit = LITERATURE_IC50.get(key)
    if lit is not None:
        return {
            "ic50_uM": float(lit["ic50_uM"]),
            "source": "literature",
            "confidence": "literature",
            "note": lit["note"],
        }

    return {
        "ic50_uM": _DEFAULT_IC50_UM,
        "source": "default",
        "confidence": "default_LOW",
        "note": _DEFAULT_NOTE,
    }


def ic50_for(drug_name: str, cache_dir: str = "data_cache") -> float:
    """Thin helper: just the resolved IC50 (uM) as a float."""
    return float(resolve_ic50(drug_name, cache_dir=cache_dir)["ic50_uM"])


if __name__ == "__main__":
    print("=" * 74)
    print("potency_breast.resolve_ic50 self-test: real GDSC breast (BRCA) potency")
    print("=" * 74)
    _ALL_MEASURED = ("paclitaxel", "cyclophosphamide", "tamoxifen", "docetaxel",
                     "epirubicin", "5-fluorouracil", "vinorelbine",
                     "palbociclib", "cisplatin")
    for name in _ALL_MEASURED + ("totally_unknown_drug",):
        res = resolve_ic50(name)
        n = res.get("n_cell_lines", "-")
        print(f"{name:20s} ic50={res['ic50_uM']:>10.4f} uM  {res['source']:>10s} "
              f"{res['confidence']:>12s}  n_lines={n}")
        assert res["ic50_uM"] > 0.0 and math.isfinite(res["ic50_uM"])
    for name in _ALL_MEASURED:
        r = resolve_ic50(name)
        assert r["source"] == "GDSC" and r["confidence"] == "measured", name
        assert r["n_cell_lines"] > 0
    unk = resolve_ic50("totally_unknown_drug")
    assert unk["confidence"] == "default_LOW"
    print("\nAll self-tests passed.")
