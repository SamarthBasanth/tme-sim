"""What every built-in tissue section needs on disk, and how to get it.

``app/server/cache.py`` owns the section *registry* — the seven specimens the
Studio can reconstruct — and answers "is this section available?" with an
honest reason when it is not. That reason names a document. **This module
names an action**: for every section id, the exact command that turns
"unavailable" into "runnable", and whether that command is something
``vivum acquire`` can run for you or something you must run by hand.

Two classes of section:

``AUTO``
    Every file is a plain, open, unauthenticated HTTP GET with a stable URL and
    a known byte count. ``vivum acquire --sections <id>`` fetches it.

``MANUAL``
    Open, but not a copy-paste-and-go: the deposit needs extraction, a format
    conversion, or a companion dataset. The remedy here is the literal shell
    block to run, not a pointer to one.

Nothing in this module invents a URL. Every ``AUTO`` URL below was checked live
(HTTP status + ``content-length``) on **2026-08-15**, and every byte count is
that response's ``content-length``, which also matches the counts recorded
independently in ``docs/REPRODUCE.md`` §2b. ``xenium_319`` was additionally
downloaded end-to-end and reconstructed from a clean clone on that date; the
others are verified as *fetchable* and not as *reconstructable-in-this-pass*,
and each recipe says which of the two it is in ``verified``.

This module is deliberately importable with only the core install (no ``app``
extra, no FastAPI): ``vivum sections`` has to work in the ten minutes before a
user has decided whether to install anything else.
"""
from __future__ import annotations

import os
import shutil
from collections import OrderedDict
from typing import Callable, Dict, List, Optional, Tuple

#: 10x Genomics public sample CDN. Open, no login, no DUA.
_CF = "https://cf.10xgenomics.com/samples"

_XENIUM_319_BASE = (
    f"{_CF}/xenium/1.3.0/Xenium_V1_FFPE_Human_Brain_Glioblastoma_With_Addon/"
    "Xenium_V1_FFPE_Human_Brain_Glioblastoma_With_Addon"
)
_XENIUM_IO_BASE = (
    f"{_CF}/xenium/2.0.0/Xenium_V1_Human_Brain_GBM_FFPE/"
    "Xenium_V1_Human_Brain_GBM_FFPE"
)
_VISIUM_PARENT_BASE = (
    f"{_CF}/spatial-exp/1.2.0/Parent_Visium_Human_Glioblastoma/"
    "Parent_Visium_Human_Glioblastoma"
)

AUTO = "auto"
MANUAL = "manual"

#: The section every "first result" path in this repo goes through.
DEFAULT_SECTION = "xenium_319"

#: The Visium deconvolution reference. Without it the Visium sections load but
#: cannot be reconstructed, so no Visium remedy is complete without it.
_SCREF_BLOCK = """mkdir -p data_cache/scref_gbm
curl -L -o data_cache/scref_gbm/GSE84465_GBM_All_data.csv.gz \\
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84465/suppl/GSE84465_GBM_All_data.csv.gz
curl -L -o data_cache/scref_gbm/GSE84465_series_matrix.txt.gz \\
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84465/matrix/GSE84465_series_matrix.txt.gz"""


def _geo_visium_block(gsm: str, sample: str) -> str:
    return (
        "mkdir -p data_cache/visium_gbm && cd data_cache/visium_gbm\n"
        f"curl -L -O https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5833nnn/{gsm}/"
        f"suppl/{gsm}_{sample}_spaceranger_out.tar.gz\n"
        f"mkdir -p {sample} && tar xzf {gsm}_{sample}_spaceranger_out.tar.gz "
        f"-C {sample} --strip-components=1\n"
        "cd ../.."
    )


#: id -> recipe. Ids are exactly ``app/server/cache.py``'s ``SECTIONS`` keys;
#: :func:`section_status` asserts nothing about that, it simply looks them up,
#: so a section this module does not know about degrades to "no recipe" rather
#: than disappearing.
SECTION_RECIPES: "OrderedDict[str, Dict]" = OrderedDict([
    ("xenium_319", {
        "id": "xenium_319",
        "label": "Xenium GBM · 319-panel",
        "platform": "xenium",
        "cache_dir": os.path.join("data_cache", "xenium_gbm"),
        "mode": AUTO,
        "bytes": 16_598_845,
        "files": [
            (f"{_XENIUM_319_BASE}_cells.csv.gz",
             "Xenium_GBM_cells.csv.gz", 2_167_253),
            (f"{_XENIUM_319_BASE}_cell_feature_matrix.tar.gz",
             "Xenium_GBM_cell_feature_matrix.tar.gz", 14_279_458),
            (f"{_XENIUM_319_BASE}_gene_panel.json",
             "Xenium_GBM_gene_panel.json", 152_134),
        ],
        "note": "the smallest runnable section, and the one every example in "
                "README.md and docs/QUICKSTART.md uses",
        "verified": "downloaded end-to-end from a clean clone on 2026-08-15 "
                    "and reconstructed + simulated from the result",
    }),
    ("xenium_io", {
        "id": "xenium_io",
        "label": "Xenium GBM · IO 480-panel",
        "platform": "xenium",
        "cache_dir": os.path.join("data_cache", "xenium_gbm_io"),
        "mode": AUTO,
        "bytes": 287_623_295,
        "files": [
            (f"{_XENIUM_IO_BASE}_cells.csv.gz",
             "Xenium_GBM_cells.csv.gz", 44_907_408),
            (f"{_XENIUM_IO_BASE}_cell_feature_matrix.tar.gz",
             "Xenium_GBM_cell_feature_matrix.tar.gz", 242_459_483),
            (f"{_XENIUM_IO_BASE}_gene_panel.json",
             "Xenium_GBM_gene_panel.json", 256_404),
        ],
        "note": "20x more cells than xenium_319 and correspondingly slower to "
                "reconstruct; not needed for a first run",
        "verified": "downloaded end-to-end from a clean clone on 2026-08-15 "
                    "(3 files, 287,623,295 B, 4 m 14 s) and reconstructed from "
                    "the result",
    }),
    ("gse324268", {
        "id": "gse324268",
        "label": "Xenium GBM · GSE324268 tumor 07",
        "platform": "xenium",
        "cache_dir": os.path.join("data_cache", "gse324268_xenium"),
        "mode": MANUAL,
        "bytes": None,
        "manual": (
            "# GEO deposit (open, no DUA) — file list:\n"
            "#   https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9572nnn/GSM9572843/suppl/\n"
            "# GEO ships this section as .parquet.gz + a .h5 matrix. The loader\n"
            "# reads a CSV cell table + an MEX/tarball matrix, so this deposit\n"
            "# needs a format conversion that this repo does not ship a script\n"
            "# for. There is no verified one-command recipe; do not expect one.\n"
            "# Route that does work on an arbitrary export:\n"
            "vivum ingest /path/to/your/converted/export --register"
        ),
        "note": "no automated recipe exists in this repo — see 'manual' above",
        "verified": "GEO file listing read live 2026-08-15; the conversion step "
                    "is NOT scripted here and was NOT re-run in that pass",
    }),
    ("visium_10x_parent_gbm", {
        "id": "visium_10x_parent_gbm",
        "label": "Visium GBM · 10x parent section",
        "platform": "visium",
        "cache_dir": os.path.join("data_cache", "visium_gbm"),
        "mode": MANUAL,
        "bytes": 63_269_841,
        "manual": (
            "mkdir -p data_cache/visium_gbm && cd data_cache/visium_gbm\n"
            f"curl -L -O {_VISIUM_PARENT_BASE}_filtered_feature_bc_matrix.tar.gz\n"
            f"curl -L -O {_VISIUM_PARENT_BASE}_spatial.tar.gz\n"
            "tar xzf Parent_Visium_Human_Glioblastoma_filtered_feature_bc_matrix.tar.gz\n"
            "tar xzf Parent_Visium_Human_Glioblastoma_spatial.tar.gz\n"
            "cd ../..\n"
            "# every Visium section ALSO needs the deconvolution reference:\n"
            + _SCREF_BLOCK
        ),
        "needs": "data_cache/scref_gbm (Darmanis GSE84465) — the Visium path "
                 "cannot run at all without it",
        "note": "spot-level, not single-cell: the spatial-attribution statistic "
                "is refused in code on this section",
        "verified": "URLs and byte counts checked live 2026-08-15 (HTTP 200, "
                    "content-length); extraction NOT re-run in that pass",
    }),
    ("visium_gbm1", {
        "id": "visium_gbm1",
        "label": "Visium GBM · GBM1 (primary)",
        "platform": "visium",
        "cache_dir": os.path.join("data_cache", "visium_gbm", "GBM1"),
        "mode": MANUAL,
        "bytes": 258_814_702,
        "manual": _geo_visium_block("GSM5833533", "GBM1")
                  + "\n# and the deconvolution reference:\n" + _SCREF_BLOCK,
        "needs": "data_cache/scref_gbm (Darmanis GSE84465)",
        "note": "GEO GSE194329; needed for the cross-tissue stage",
        "verified": "URL and byte count checked live 2026-08-15 (HTTP 200, "
                    "content-length); extraction NOT re-run in that pass",
    }),
    ("visium_gbm2", {
        "id": "visium_gbm2",
        "label": "Visium GBM · GBM2 (recurrent)",
        "platform": "visium",
        "cache_dir": os.path.join("data_cache", "visium_gbm", "GBM2"),
        "mode": MANUAL,
        "bytes": None,
        "manual": _geo_visium_block("GSM5833534", "GBM2")
                  + "\n# and the deconvolution reference:\n" + _SCREF_BLOCK,
        "needs": "data_cache/scref_gbm (Darmanis GSE84465)",
        "note": "GEO GSE194329; needed for the cross-tissue stage",
        "verified": "URL pattern from docs/REPRODUCE.md §2b; the GBM1 sibling "
                    "on the same host was checked live 2026-08-15, this "
                    "accession was NOT",
    }),
    ("visium_gbm3", {
        "id": "visium_gbm3",
        "label": "Visium GBM · GBM3 (primary)",
        "platform": "visium",
        "cache_dir": os.path.join("data_cache", "visium_gbm", "GBM3"),
        "mode": MANUAL,
        "bytes": None,
        "manual": _geo_visium_block("GSM5833535", "GBM3")
                  + "\n# and the deconvolution reference:\n" + _SCREF_BLOCK,
        "needs": "data_cache/scref_gbm (Darmanis GSE84465)",
        "note": "GEO GSE194329; needed for the cross-tissue stage",
        "verified": "URL pattern from docs/REPRODUCE.md §2b; the GBM1 sibling "
                    "on the same host was checked live 2026-08-15, this "
                    "accession was NOT",
    }),
])


def auto_section_ids() -> List[str]:
    """Section ids ``vivum acquire --sections`` can fetch, in catalog order."""
    return [k for k, v in SECTION_RECIPES.items() if v["mode"] == AUTO]


def manual_section_ids() -> List[str]:
    """Section ids that need a documented human step."""
    return [k for k, v in SECTION_RECIPES.items() if v["mode"] == MANUAL]


def total_auto_bytes(ids: Optional[List[str]] = None) -> int:
    """Sum of the declared byte counts of the given auto sections."""
    ids = list(SECTION_RECIPES) if ids is None else ids
    return sum(int(SECTION_RECIPES[i]["bytes"] or 0)
               for i in ids
               if i in SECTION_RECIPES and SECTION_RECIPES[i]["mode"] == AUTO)


def human_bytes(n: Optional[int]) -> str:
    """``16598845`` -> ``'16.6 MB'``. Never rounds an unknown into a number."""
    if not n:
        return "unknown size"
    for unit, scale in (("GB", 1e9), ("MB", 1e6), ("kB", 1e3)):
        if n >= scale:
            return f"{n / scale:.1f} {unit}"
    return f"{n} B"


# --------------------------------------------------------------------------- #
# The remedy: what to DO about an unavailable section
# --------------------------------------------------------------------------- #
def remedy_for(section_id: str) -> List[str]:
    """Printable lines telling a user how to make ``section_id`` runnable.

    Returns ``[]`` for a section with no recipe (an ingested one, say — those
    already exist on disk by construction), so a caller can fall through to
    whatever it printed before. Never raises.
    """
    rec = SECTION_RECIPES.get(section_id)
    if rec is None:
        return []
    out: List[str] = []
    if rec["mode"] == AUTO:
        out.append(f"fetch it:  vivum acquire --sections {section_id}")
        out.append(f"           {human_bytes(rec.get('bytes'))}, "
                   f"{len(rec['files'])} files, into {rec['cache_dir']}/")
    else:
        out.append("this section is NOT auto-fetchable. Run, from the repo root:")
        out.extend("    " + ln for ln in str(rec["manual"]).splitlines())
    if rec.get("needs"):
        out.append(f"also needs: {rec['needs']}")
    return out


def remedy_text(section_id: str, indent: str = "  ") -> str:
    """:func:`remedy_for` as one indented block. Empty string when unknown."""
    lines = remedy_for(section_id)
    return "\n".join(indent + ln for ln in lines)


# --------------------------------------------------------------------------- #
# Availability, without needing the app extra
# --------------------------------------------------------------------------- #
def _xenium_available(cache_dir: str) -> Tuple[bool, Optional[str]]:
    """Mirror of ``app.server.cache._xenium_available`` on the same primitives.

    Both call ``tmesim.data.ingest.resolve_xenium_files``, so the two agree by
    construction rather than by a copied file list.
    """
    if not os.path.isdir(cache_dir):
        return False, f"cache directory {cache_dir} does not exist"
    try:
        from tmesim.data.ingest import resolve_xenium_files

        found = resolve_xenium_files(cache_dir)
    except Exception as exc:  # noqa: BLE001
        return False, f"could not scan {cache_dir}: {type(exc).__name__}: {exc}"
    if not found.get("cells"):
        return False, f"{cache_dir} has no per-cell centroid table"
    if not (found.get("matrix_dir") or found.get("matrix_tar")):
        return False, (f"{cache_dir} has neither an extracted cell-feature "
                       "matrix nor its tarball")
    return True, None


def _visium_available(cache_dir: str, sample: Optional[str]) -> Tuple[bool, Optional[str]]:
    try:
        from tmesim.data.real import visium
    except Exception as exc:  # noqa: BLE001
        return False, f"visium loader unimportable: {type(exc).__name__}: {exc}"
    try:
        present = set(visium.list_samples(present_only=True))
    except Exception as exc:  # noqa: BLE001
        return False, f"visium cache scan failed: {type(exc).__name__}: {exc}"
    if sample in present:
        return True, None
    return False, f"Visium sample {sample!r} is not cached under {cache_dir}"


#: ``id -> sample`` for the Visium recipes, matching ``app/server/cache.py``.
_VISIUM_SAMPLE = {
    "visium_10x_parent_gbm": "10x_parent_gbm",
    "visium_gbm1": "GBM1",
    "visium_gbm2": "GBM2",
    "visium_gbm3": "GBM3",
}

#: ``id -> --sample`` argument, matching ``app/server/cache.py``'s ``SECTIONS``.
#: ``None`` means "the loader default", which is what ``xenium_319`` is.
_SAMPLE_ARG = {
    "xenium_319": None,
    "xenium_io": os.path.join("data_cache", "xenium_gbm_io"),
    "gse324268": os.path.join("data_cache", "gse324268_xenium"),
    **_VISIUM_SAMPLE,
}


def run_command(section_id: str, verb: str = "simulate",
                extra: str = "--drug temozolomide") -> Optional[str]:
    """The literal ``vivum`` command that runs ``verb`` on this section.

    Uses ``--dx 30`` because this string is shown to someone who has not run
    anything yet and a first result should arrive in seconds; every number in
    ``results/`` is at ``--dx 12`` and the docs say so where they quote one.
    """
    rec = SECTION_RECIPES.get(section_id)
    if rec is None:
        return None
    sample = _SAMPLE_ARG.get(section_id)
    cmd = f"vivum {verb} --source {rec['platform']} --dx 30"
    if sample:
        cmd += f" --sample {sample}"
    return f"{cmd} {extra}".strip()


def section_status(section_id: str) -> Dict:
    """``{id, label, available, reason, remedy, mode}`` for one catalog id.

    Uses the same availability primitives the Studio uses, so a section that
    reads ``available`` here reads ``available`` in ``GET /api/sections`` too.
    Never raises.
    """
    rec = SECTION_RECIPES.get(section_id)
    if rec is None:
        return {"id": section_id, "label": section_id, "available": False,
                "reason": "not a built-in section", "remedy": [], "mode": None}
    if rec["platform"] == "visium":
        ok, why = _visium_available(rec["cache_dir"],
                                    _VISIUM_SAMPLE.get(section_id))
    else:
        ok, why = _xenium_available(rec["cache_dir"])
    return {
        "id": section_id,
        "label": rec["label"],
        "available": ok,
        "reason": why,
        "remedy": [] if ok else remedy_for(section_id),
        "mode": rec["mode"],
        "cache_dir": rec["cache_dir"],
    }


def all_section_status() -> List[Dict]:
    """:func:`section_status` for every catalog section, in catalog order."""
    return [section_status(sid) for sid in SECTION_RECIPES]


def n_available() -> Tuple[int, int]:
    """``(available, total)`` over the catalog. Never raises."""
    st = all_section_status()
    return sum(1 for s in st if s["available"]), len(st)


# --------------------------------------------------------------------------- #
# Fetching an AUTO section
# --------------------------------------------------------------------------- #
def _download(url: str, dest: str, expect_bytes: int, timeout: int,
              log: Callable[[str], None]) -> Tuple[bool, str]:
    """Stream ``url`` to ``dest``, verifying the byte count. Never raises.

    A file already on disk at exactly ``expect_bytes`` is left alone — re-running
    ``acquire`` after a partial failure resumes at file granularity rather than
    re-downloading what already landed. A short or long file is DELETED and
    reported failed: a truncated ``.csv.gz`` that stays on disk would surface
    later as an unreadable-cache error a long way from its cause.
    """
    from .acquire import _get  # local import: keeps the import graph one-way

    name = os.path.basename(dest)
    if os.path.exists(dest) and os.path.getsize(dest) == expect_bytes:
        return True, f"{name}: already present ({expect_bytes:,} B)"

    part = dest + ".part"
    r = _get(url, timeout=max(timeout, 120), stream=True)
    if r is None:
        return False, f"{name}: download failed (see the [acquire] line above)"
    try:
        with open(part, "wb") as fh:
            for chunk in r.iter_content(chunk_size=1 << 16):
                if chunk:
                    fh.write(chunk)
    except Exception as exc:  # noqa: BLE001
        _rm(part)
        return False, f"{name}: write failed: {exc}"
    finally:
        r.close()

    got = os.path.getsize(part)
    if got != expect_bytes:
        _rm(part)
        return False, (f"{name}: size mismatch — got {got:,} B, expected "
                       f"{expect_bytes:,} B. Nothing was kept; re-run to retry.")
    try:
        shutil.move(part, dest)
    except Exception as exc:  # noqa: BLE001
        _rm(part)
        return False, f"{name}: could not finalise: {exc}"
    log(f"  {name}: {got:,} B")
    return True, f"{name}: {got:,} B"


def _rm(path: str) -> None:
    try:
        os.remove(path)
    except OSError:
        pass


def fetch_section(section_id: str, cache_root: str = "data_cache",
                  timeout: int = 60,
                  log: Optional[Callable[[str], None]] = None) -> Dict:
    """Download one AUTO section into ``cache_root``. Never raises.

    Returns the same status shape the other acquire sources use::

        {"status": "ok"|"failed"|"skipped", "path": ..., "note": ...,
         "n_records": <files landed>}

    ``skipped`` is what a MANUAL section returns, and its ``note`` is the
    instruction — a caller that prints notes prints the remedy for free.
    """
    log = log or (lambda _m: None)
    rec = SECTION_RECIPES.get(section_id)
    if rec is None:
        return {"status": "failed", "path": None, "n_records": 0,
                "note": f"unknown section {section_id!r}; known: "
                        f"{', '.join(SECTION_RECIPES)}"}
    if rec["mode"] != AUTO:
        return {"status": "skipped", "path": None, "n_records": 0,
                "note": "not auto-fetchable — run `vivum sections "
                        f"{section_id}` for the exact commands"}

    # The recipe's cache_dir is repo-relative by construction; honour a caller
    # that redirected the cache root (``acquire --out``) by swapping the head.
    tail = os.path.relpath(rec["cache_dir"], "data_cache")
    dest_dir = os.path.join(cache_root, tail)
    os.makedirs(dest_dir, exist_ok=True)

    notes: List[str] = []
    n_ok = 0
    for url, basename, size in rec["files"]:
        ok, msg = _download(url, os.path.join(dest_dir, basename), size,
                            timeout, log)
        notes.append(msg)
        n_ok += int(ok)

    if n_ok != len(rec["files"]):
        return {"status": "failed", "path": dest_dir, "n_records": n_ok,
                "note": f"{n_ok}/{len(rec['files'])} files — "
                        + "; ".join(notes)}
    return {"status": "ok", "path": dest_dir, "n_records": n_ok,
            "note": f"{n_ok} files, {human_bytes(rec['bytes'])}"}


def resolve_section_arg(arg: Optional[str]) -> Tuple[List[str], List[str]]:
    """Parse ``--sections``. Returns ``(ids, unknown)``.

    ``None`` -> the default section. ``"none"`` -> nothing. ``"all"`` -> every
    AUTO section. Otherwise a comma-separated list, which MAY name a MANUAL
    section: :func:`fetch_section` then reports it as skipped with its remedy,
    which is more useful than an argparse error.
    """
    if arg is None:
        return [DEFAULT_SECTION], []
    raw = [s.strip() for s in arg.split(",") if s.strip()]
    if len(raw) == 1 and raw[0].lower() == "none":
        return [], []
    if len(raw) == 1 and raw[0].lower() == "all":
        return auto_section_ids(), []
    known = [s for s in raw if s in SECTION_RECIPES]
    unknown = [s for s in raw if s not in SECTION_RECIPES]
    return known, unknown
