"""Real potency (IC50) resolution for the molecular layer.

Closing the "potency gap": in the molecular contract, IC50 is the single weakest
number — it is NOT derivable from structure without docking/assay, so
``compile_drug.build_profile`` uses a flagged low-confidence DEFAULT (25 uM) when
no IC50 is supplied. This module makes potency REAL wherever measured data exists,
without touching ``core``, ``types``, the engine, or any existing molecular module.

Three tiers, in strict order of trust:

  1. MEASURED  (``measured_ic50``)  -- median IC50 across glioma cell lines from
     the cached GDSC2 dose-response screen (``data_cache/gdsc_glioma_dose_response.csv``,
     downloaded by ``tmesim.data.acquire``). This is genuine wet-lab ground truth.
     GDSC reports ``LN_IC50`` = natural log of the fitted IC50 in micromolar; we
     aggregate per drug across glioma (GBM/LGG) cell lines by median in log-space
     and convert back with ``exp``. confidence = 'measured'.

  2. LITERATURE (``LITERATURE_IC50``) -- a small curated table of GBM-relevant
     IC50 ranges for library drugs that GDSC does NOT screen (e.g. the
     nitrosoureas, doxorubicin, etoposide, vincristine). Each entry cites a
     rough well-known range. confidence = 'literature'.

  3. DEFAULT   (flagged) -- the same 25 uM low-confidence placeholder the bridge
     already uses, returned only when neither measured nor literature data exists.
     confidence = 'default_LOW'.

Units: IC50 in micromolar (uM) throughout.

GDSC name matching (case-insensitive, documented)
--------------------------------------------------
Our library uses lowercase generic names (temozolomide, cisplatin, ...). GDSC
``DRUG_NAME`` uses Title-Case generic names plus, occasionally, code names or
salt forms. We match by:
  * exact case-insensitive equality of the library name to a GDSC DRUG_NAME, then
  * a curated synonym/alias list per library drug (brand/abbreviation/salt forms,
    e.g. TMZ, CDDP/DDP for cisplatin, MTX, VP-16 for etoposide, VCR for
    vincristine, CCNU/BCNU for the nitrosoureas).
For the current cached glioma screen the matches that actually resolve are:
    temozolomide -> 'Temozolomide', cisplatin -> 'Cisplatin',
    methotrexate -> 'Methotrexate'.
Everything else in the library (doxorubicin, etoposide, vincristine, lomustine,
carmustine, procarbazine, bevacizumab) is absent from this glioma screen and
falls through to the literature table or the flagged default.
"""
from __future__ import annotations

import csv
import math
import os
import statistics
from typing import Dict, List, Optional

# ---------------------------------------------------------------------------
# GDSC name matching: library canonical name -> accepted GDSC DRUG_NAME aliases.
# All comparisons are case-insensitive. The canonical name itself is always an
# implicit alias; the extra strings cover abbreviations / brand / salt forms so
# the matcher stays robust if GDSC's naming shifts between releases.
# ---------------------------------------------------------------------------
_GDSC_ALIASES: Dict[str, List[str]] = {
    "temozolomide": ["temozolomide", "tmz", "temodar", "temodal"],
    "cisplatin": ["cisplatin", "cddp", "ddp", "cis-platinum", "cisplatinum"],
    "methotrexate": ["methotrexate", "mtx", "amethopterin"],
    "doxorubicin": ["doxorubicin", "dox", "adriamycin", "hydroxydaunorubicin"],
    "etoposide": ["etoposide", "vp-16", "vp16"],
    "vincristine": ["vincristine", "vcr", "oncovin", "leurocristine"],
    "lomustine": ["lomustine", "ccnu"],
    "carmustine": ["carmustine", "bcnu"],
    "procarbazine": ["procarbazine", "matulane"],
    # bevacizumab (antibody) has no cytotoxic IC50 in a cell-line viability screen.
}


# ---------------------------------------------------------------------------
# Tier 2: curated literature fallback (GBM-relevant IC50s, uM). Used ONLY when
# GDSC has no measured value. These are deliberately rough order-of-magnitude
# anchors from well-known glioma / cancer cell-line ranges, each with a source
# note; confidence is 'literature', never 'measured'.
# ---------------------------------------------------------------------------
LITERATURE_IC50: Dict[str, Dict] = {
    "doxorubicin": {
        "ic50_uM": 0.5,
        "note": (
            "Anthracycline; sub-micromolar in glioma cell lines (~0.1-1 uM). "
            "Poor BBB penetration limits clinical GBM use. Range from cancer "
            "cell-line viability screens (e.g. NCI-60 / GDSC anthracyclines)."
        ),
    },
    "etoposide": {
        "ic50_uM": 5.0,
        "note": (
            "Topoisomerase-II poison (VP-16); glioma cell-line IC50 ~1-10 uM. "
            "Literature range from GBM in-vitro cytotoxicity studies."
        ),
    },
    "vincristine": {
        "ic50_uM": 0.05,
        "note": (
            "Vinca alkaloid (microtubule inhibitor); highly potent, low-nM to "
            "~0.1 uM in cell lines. Historically part of the PCV regimen for "
            "glioma. Range from vinca-alkaloid cytotoxicity literature."
        ),
    },
    "carmustine": {
        "ic50_uM": 30.0,
        "note": (
            "Nitrosourea (BCNU), lipophilic alkylator used in GBM (Gliadel "
            "wafers); glioma cell-line IC50 tens of uM. Literature order-of-"
            "magnitude anchor."
        ),
    },
    "lomustine": {
        "ic50_uM": 50.0,
        "note": (
            "Nitrosourea (CCNU), oral GBM alkylator (PCV / lomustine arms); "
            "glioma cell-line IC50 tens of uM. Literature order-of-magnitude "
            "anchor."
        ),
    },
    "procarbazine": {
        "ic50_uM": 500.0,
        "note": (
            "Methylating pro-drug requiring hepatic activation; weak direct "
            "cytotoxicity in vitro (hundreds of uM). Literature order-of-"
            "magnitude anchor."
        ),
    },
    # temozolomide / cisplatin / methotrexate intentionally omitted here: they
    # resolve to MEASURED GDSC values. Backups are provided below only so the
    # literature tier still degrades gracefully if the cache is missing.
    "temozolomide": {
        "ic50_uM": 100.0,
        "note": (
            "Backup only (GDSC-measured normally wins). TMZ is a weak in-vitro "
            "methylator; glioma cell-line IC50 is large (tens to hundreds of uM), "
            "consistent with the measured GDSC glioma median."
        ),
    },
    "cisplatin": {
        "ic50_uM": 10.0,
        "note": (
            "Backup only (GDSC-measured normally wins). Platinum crosslinker; "
            "glioma cell-line IC50 ~single-digit to tens of uM."
        ),
    },
    "methotrexate": {
        "ic50_uM": 10.0,
        "note": (
            "Backup only (GDSC-measured normally wins). Antifolate; cell-line "
            "IC50 spans widely with folate context, ~single-digit to tens of uM."
        ),
    },
}


# ---------------------------------------------------------------------------
# The flagged default (mirrors compile_drug._DEFAULT_IC50_UM). Kept in sync by
# value; imported lazily so this module has no hard dependency on the bridge.
# ---------------------------------------------------------------------------
_DEFAULT_IC50_UM = 25.0
_DEFAULT_NOTE = (
    f"flagged default {_DEFAULT_IC50_UM} uM; confidence LOW -- no measured (GDSC) "
    "or curated-literature IC50 for this drug. Potency is NOT derivable from "
    "structure without docking/assay. Replace with a measured IC50."
)

_CACHE_FILENAME = "gdsc_glioma_dose_response.csv"


def _gdsc_path(cache_dir: str) -> str:
    return os.path.join(cache_dir, _CACHE_FILENAME)


def _aliases_for(drug_name: str) -> List[str]:
    """Return the lowercase accepted GDSC aliases for a library drug.

    Always includes the drug's own lowercased name, so drugs not in the curated
    alias map still match an identically-named GDSC entry.
    """
    key = drug_name.strip().lower()
    aliases = set(_GDSC_ALIASES.get(key, []))
    aliases.add(key)
    return [a.lower() for a in aliases]


def measured_ic50(drug_name: str, cache_dir: str = "data_cache") -> Optional[Dict]:
    """Look up a drug's MEASURED IC50 from the cached GDSC glioma screen.

    GDSC reports ``LN_IC50`` = natural log of the fitted IC50 in micromolar, one
    row per (cell line, drug). We keep only rows whose ``DRUG_NAME`` matches this
    library drug (case-insensitive; via the curated alias map), collapse any
    repeat screens of the same cell line by median, then take the median
    ``LN_IC50`` across the distinct glioma cell lines and convert to uM::

        ic50_uM = exp(median_over_cell_lines(LN_IC50))

    (Median in log-space == geometric median IC50, the standard robust summary
    for dose-response data.)

    Parameters
    ----------
    drug_name : str
        Library drug name (e.g. "temozolomide"); case-insensitive.
    cache_dir : str
        Directory holding ``gdsc_glioma_dose_response.csv``.

    Returns
    -------
    dict or None
        ``{ic50_uM, n_cell_lines, source:'GDSC', confidence:'measured',
        raw_median_ln_ic50, matched_gdsc_name}`` if the drug is found with at
        least one usable glioma cell line, else ``None`` (cache missing, drug not
        screened, or no parseable LN_IC50).
    """
    path = _gdsc_path(cache_dir)
    if not os.path.exists(path):
        return None

    aliases = set(_aliases_for(drug_name))

    # cell_line -> list of LN_IC50 for the matched drug; also remember the exact
    # GDSC DRUG_NAME string(s) we matched (for provenance / reporting).
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

    # One representative LN_IC50 per distinct cell line (median over any repeats),
    # then the across-cell-line median.
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
    """Resolve an IC50 (uM) for a drug using the best available evidence tier.

    Order of trust: measured GDSC -> curated literature -> flagged default. The
    return value ALWAYS has the keys ``ic50_uM``, ``source``, ``confidence``,
    ``note``. Measured results additionally carry ``n_cell_lines`` and
    ``raw_median_ln_ic50`` so callers can report the provenance.

    Parameters
    ----------
    drug_name : str
        Library drug name; case-insensitive.
    cache_dir : str
        Directory holding the GDSC cache (for the measured tier).

    Returns
    -------
    dict
        Always ``{ic50_uM, source, confidence, note, ...}``.
    """
    # Tier 1: measured GDSC.
    m = measured_ic50(drug_name, cache_dir=cache_dir)
    if m is not None:
        return {
            "ic50_uM": m["ic50_uM"],
            "source": m["source"],  # 'GDSC'
            "confidence": m["confidence"],  # 'measured'
            "note": (
                f"median across {m['n_cell_lines']} glioma cell lines "
                f"(GDSC2 {m['matched_gdsc_name']}); raw median LN_IC50="
                f"{m['raw_median_ln_ic50']:.4f}; ic50_uM=exp(median LN_IC50)."
            ),
            "n_cell_lines": m["n_cell_lines"],
            "raw_median_ln_ic50": m["raw_median_ln_ic50"],
        }

    # Tier 2: curated literature.
    key = drug_name.strip().lower()
    lit = LITERATURE_IC50.get(key)
    if lit is not None:
        return {
            "ic50_uM": float(lit["ic50_uM"]),
            "source": "literature",
            "confidence": "literature",
            "note": lit["note"],
        }

    # Tier 3: flagged default.
    return {
        "ic50_uM": _DEFAULT_IC50_UM,
        "source": "default",
        "confidence": "default_LOW",
        "note": _DEFAULT_NOTE,
    }


def ic50_for(drug_name: str, cache_dir: str = "data_cache") -> float:
    """Thin helper: just the resolved IC50 (uM) as a float.

    Lets the bridge obtain a real potency WITHOUT modifying ``compile_drug``:
    call ``ic50_for("cisplatin")`` and pass the result as the ``ic50`` argument
    to ``build_profile`` / ``molecule_to_drugparams``.
    """
    return float(resolve_ic50(drug_name, cache_dir=cache_dir)["ic50_uM"])


if __name__ == "__main__":
    print("=" * 74)
    print("potency.resolve_ic50 self-test: closing the IC50 gap with real data")
    print("=" * 74)

    demo = ["temozolomide", "cisplatin", "doxorubicin", "methotrexate",
            "totally_unknown_drug"]

    measured_hits: List[str] = []
    fallbacks: List[str] = []

    hdr = f"{'drug':22s} {'ic50_uM':>12s} {'source':>12s} {'confidence':>13s} {'n_lines':>8s}"
    print(hdr)
    print("-" * len(hdr))
    for name in demo:
        res = resolve_ic50(name)
        n = res.get("n_cell_lines", "-")
        print(f"{name:22s} {res['ic50_uM']:>12.4f} {res['source']:>12s} "
              f"{res['confidence']:>13s} {str(n):>8s}")
        # sanity: value always positive & finite
        assert res["ic50_uM"] > 0.0 and math.isfinite(res["ic50_uM"]), name
        if res["confidence"] == "measured":
            measured_hits.append(name)
        else:
            fallbacks.append(f"{name} ({res['confidence']})")

    # --- assertions -------------------------------------------------------
    # The three GDSC-screened glioma drugs must resolve to MEASURED values.
    for name in ("temozolomide", "cisplatin", "methotrexate"):
        r = resolve_ic50(name)
        assert r["source"] == "GDSC" and r["confidence"] == "measured", (
            f"{name} should be measured, got {r['source']}/{r['confidence']}")
        assert r["ic50_uM"] > 0.0 and math.isfinite(r["ic50_uM"])
        assert r["n_cell_lines"] > 0

    # Doxorubicin is NOT in this glioma screen -> literature fallback.
    dox = resolve_ic50("doxorubicin")
    assert dox["confidence"] == "literature", dox
    assert measured_ic50("doxorubicin") is None

    # Unknown drug -> flagged default.
    unk = resolve_ic50("totally_unknown_drug")
    assert unk["confidence"] == "default_LOW" and unk["ic50_uM"] == _DEFAULT_IC50_UM

    # ic50_for returns a bare float and agrees with resolve_ic50.
    assert isinstance(ic50_for("cisplatin"), float)
    assert abs(ic50_for("cisplatin") - resolve_ic50("cisplatin")["ic50_uM"]) < 1e-9

    print("\nREAL measured (GDSC) IC50s:")
    for name in measured_hits:
        r = resolve_ic50(name)
        print(f"  {name:16s} ic50 = {r['ic50_uM']:9.3f} uM  "
              f"(n={r['n_cell_lines']} glioma cell lines, "
              f"raw median LN_IC50={r['raw_median_ln_ic50']:.4f})")
    print("\nFallback (no GDSC glioma measurement):")
    for f in fallbacks:
        print(f"  {f}")

    print("\nAll potency self-tests passed.")
