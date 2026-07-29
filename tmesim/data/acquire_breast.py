"""Best-effort LIVE download of the GDSC2 BREAST-cancer (BRCA) dose-response
screen into ``./data_cache/``.

This is the breast-cancer analogue of ``tmesim.data.acquire._acquire_gdsc_glioma``
-- SAME source (Sanger GDSC2 fitted dose-response CSV), SAME streaming-filter
pattern, only the kept ``TCGA_DESC`` value changes (``BRCA`` instead of
``GBM``/``LGG``). A new, additive file: nothing in ``tmesim.data.acquire`` is
modified, and the existing ``data_cache/gdsc_glioma_dose_response.csv`` cache
(and every GBM result built on it) is untouched.

Why a separate cache file rather than re-filtering the glioma one: the glioma
acquirer only ever writes GBM/LGG rows to disk, so the BRCA rows are not
recoverable from that file -- the full ~40 MB GDSC2 CSV has to be streamed and
re-filtered once, independently.

Verified live 2026-07-28 (same GDSC release 8.5 as the glioma acquirer):
the BRCA-filtered stream yields 13,106 rows across 50-51 breast cell lines per
drug for paclitaxel / cyclophosphamide / tamoxifen / docetaxel / epirubicin /
5-fluorouracil (doxorubicin is NOT screened in GDSC2 at all, for any tissue --
same gap the glioma module already documents; the breast panel handles this by
reusing the GBM-panel's literature-tier doxorubicin IC50 verbatim, see
``experiments/breast_drug_panel.py``).

Run standalone::

    python -m tmesim.data.acquire_breast
"""
from __future__ import annotations

import os
from typing import Dict

from .acquire import GDSC2_DOSE_RESPONSE_CSV, _get, _log

#: Which TCGA_DESC rows to keep. GDSC's breast-cancer code is "BRCA" (The
#: Cancer Genome Atlas project abbreviation for breast invasive carcinoma) --
#: confirmed by streaming the live CSV and counting TCGA_DESC values directly
#: (13,106 BRCA rows, the 4th-most-represented tissue in the whole screen).
BREAST_TCGA_DESC = ("BRCA",)

_CACHE_FILENAME = "gdsc_breast_dose_response.csv"


def acquire_gdsc_breast(cache_dir: str = "data_cache", timeout: int = 60) -> Dict:
    """Stream the GDSC2 CSV and keep only BRCA (breast cancer) rows.

    Mirrors ``tmesim.data.acquire._acquire_gdsc_glioma`` exactly (same source
    URL, same streaming-filter-write pattern, same never-raise contract), with
    the kept ``TCGA_DESC`` set swapped from ``{GBM, LGG}`` to ``{BRCA}``.

    Returns a status dict: ``{status, path, note, n_records}``. Never raises --
    a failed/partial download returns ``status='failed'`` so a caller can
    degrade gracefully (fall back to the curated literature table).
    """
    os.makedirs(cache_dir, exist_ok=True)
    out_path = os.path.join(cache_dir, _CACHE_FILENAME)

    r = _get(GDSC2_DOSE_RESPONSE_CSV, timeout=max(timeout, 120), stream=True)
    if r is None:
        return {"status": "failed", "path": None,
                "note": "GDSC CSV download failed", "n_records": 0}

    n_rows = 0
    tcga_idx = None
    try:
        with open(out_path, "w") as fh:
            for i, raw in enumerate(r.iter_lines(decode_unicode=True)):
                if raw is None:
                    continue
                if i == 0:
                    header = raw.split(",")
                    try:
                        tcga_idx = header.index("TCGA_DESC")
                    except ValueError:
                        tcga_idx = 6  # known position, fallback
                    fh.write(raw + "\n")
                    continue
                if not raw:
                    continue
                fields = raw.split(",")
                if len(fields) > tcga_idx and fields[tcga_idx] in BREAST_TCGA_DESC:
                    fh.write(raw + "\n")
                    n_rows += 1
    except Exception as exc:  # noqa: BLE001
        r.close()
        if n_rows == 0:
            return {"status": "failed", "path": None,
                    "note": f"GDSC stream failed: {exc}", "n_records": 0}
    finally:
        r.close()

    if n_rows == 0:
        return {"status": "failed", "path": None,
                "note": "no BRCA rows in GDSC stream", "n_records": 0}
    _log(f"gdsc_breast: {n_rows} breast (BRCA) cell-line dose-response rows -> {out_path}")
    return {
        "status": "ok",
        "path": out_path,
        "note": f"{n_rows} breast (BRCA) cell-line dose-response rows",
        "n_records": n_rows,
    }


if __name__ == "__main__":
    result = acquire_gdsc_breast()
    print(result)
