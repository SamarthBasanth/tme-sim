"""Breast-cancer analogue of ``tmesim.validation.kill_calibration`` (the GBM
3-D:2-D absolute-kill calibration module).

This is a DRIVER script only. It imports and re-uses the frozen, unmodified
functions in ``tmesim.validation.kill_calibration`` (the Hill mapping, the
delivery/state/sensitivity decomposition, the ratio test, the anchoring
machinery) and ``experiments.breast_drug_panel`` (the breast drug-compilation
pipeline) exactly as they exist. Nothing under ``tmesim/`` is modified.
Nothing in this file re-implements any physics, PK/PD, or calibration code --
every number here comes from calling an existing function.

Tissue: the real Xenium breast section already used by the breast flagship
(``data_cache/xenium_breast_janesick/rep1``; Janesick et al. 2023, Nat Commun
14:8353, PMID 38114474), reconstructed via the identical
``tmesim.reconstruct.real_pipeline.reconstruct_real`` call GBM uses, at the
same grid-converged ``dx = 12 um``.

Drugs: the breast flagship's panel -- doxorubicin (IC50 reused verbatim from
the GBM library, literature tier), paclitaxel / cyclophosphamide / tamoxifen
(IC50 measured live from GDSC2's breast (BRCA) screen, 50 cell lines each).

Run:    ./.venv/bin/python experiments/kill_calibration_breast.py
Writes: results/kill_calibration_breast.md
"""
from __future__ import annotations

import csv
import math
import os
import sys
import time
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, ".."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from tmesim.core import PDConfig, PDEConfig, PKConfig, SimConfig  # noqa: E402
from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams  # noqa: E402
from tmesim.molecular.breast_library import get_molecule_breast  # noqa: E402
from tmesim.molecular.compile_drug import molecule_to_drugparams  # noqa: E402
from tmesim.molecular.library import get_molecule as get_molecule_gbm  # noqa: E402
from tmesim.molecular.potency_breast import resolve_ic50 as resolve_ic50_breast  # noqa: E402
from tmesim.reconstruct.real_pipeline import (  # noqa: E402
    apply_drug_sensitivity,
    reconstruct_real,
    reconstruction_report,
)
from tmesim.sim.engine import simulate_response  # noqa: E402
from tmesim.validation.calibration import calibrate_alpha, calibrated_drugparams  # noqa: E402

# Every one of these is imported, not re-implemented: the Hill kill mapping,
# the delivery/state/sensitivity decomposition, the literature bridge, the
# ratio test, the anchoring/absolute-kill machinery, the cell-state-label
# sensitivity sweep, and the model-prologue plumbing.
from tmesim.validation.kill_calibration import (  # noqa: E402
    BATH_MULTIPLES,
    BATH_REFERENCE_PCT,
    LIT_2D3D,
    absolute_kills,
    anchor_intrinsic,
    decompose,
    expected_kill,
    ic50_fold_to_kill_ratio,
    kill_ratio_to_ic50_fold,
    kill_with,
    literature_bridge,
    p_hill,
    prepare_tissue,
    ratio_test,
    resistance_vector,
    state_breakdown,
    state_dependency,
    unit_transport,
    window_plasma_mean,
)

RESULTS_DIR = os.path.join(_ROOT, "results")
OUT_MD = os.path.join(RESULTS_DIR, "kill_calibration_breast.md")
#: On-disk cache of the (expensive) unit transport solves, gitignored -- same
#: pattern as ``tmesim.validation.kill_calibration.RUN_CACHE``, kept in a
#: separate file so the two never collide.
RUN_CACHE = os.path.join(_ROOT, "data_cache", "kill_calibration_breast_runs.json")

BREAST_SAMPLE_PATH = os.path.join(_ROOT, "data_cache", "xenium_breast_janesick", "rep1")
GDSC_BREAST_CSV = os.path.join(_ROOT, "data_cache", "gdsc_breast_dose_response.csv")

# Matches tmesim.validation.kill_calibration exactly, so the two reports are
# numerically comparable (same grid, same PDE settings, same nominal dose).
DX_UM = 12.0
DT_S = 10.0
WINDOW_S = 3600.0 * 6
SAVE_EVERY = 10
SEED = 0
NOMINAL_PEAK_UM = 50.0
NOMINAL_HALF_LIFE_S = 3600.0 * 4

DRUGS: Tuple[str, ...] = ("doxorubicin", "paclitaxel", "cyclophosphamide", "tamoxifen")
IC50_SOURCE_TAG = {
    "doxorubicin": "literature (reused verbatim from GBM library; not in GDSC2 at all)",
    "paclitaxel": "GDSC2 breast (BRCA), n=50 cell lines, measured",
    "cyclophosphamide": "GDSC2 breast (BRCA), n=50 cell lines, measured",
    "tamoxifen": "GDSC2 breast (BRCA), n=50 cell lines, measured",
}

GDSC_ALIASES = {
    "paclitaxel": ["paclitaxel", "taxol"],
    "cyclophosphamide": ["cyclophosphamide", "endoxan", "cytoxan", "ctx"],
    "tamoxifen": ["tamoxifen", "nolvadex"],
}


# --------------------------------------------------------------------------- #
# Drug compilation -- identical pipeline to experiments/breast_drug_panel.py.
# Reproduced here (rather than imported) only because that module's function
# is a private-by-convention helper inside a script guarded by
# ``if __name__ == "__main__"``; every call it makes is to the SAME frozen
# tmesim functions, none of which are modified or re-implemented.
# --------------------------------------------------------------------------- #
def compile_breast_drug(name: str, alpha: float):
    key = name.strip().lower()
    if key == "doxorubicin":
        dp = calibrated_drugparams("doxorubicin", alpha)
        mol = get_molecule_gbm("doxorubicin")
        ic50_info = {
            "ic50_uM": float(dp.ic50), "source": "literature",
            "confidence": "literature", "n_cell_lines": None,
        }
    else:
        mol = get_molecule_breast(key)
        ic50_info = resolve_ic50_breast(key)
        dp_uncal = molecule_to_drugparams(mol, ic50=ic50_info["ic50_uM"])
        dp = calibrated_drugparams(dp_uncal, alpha)
    dp = apply_barrier_to_drugparams(dp, mol)
    return dp, mol, ic50_info


def gdsc_breast_ic50_band(drug_name: str) -> Optional[Dict]:
    """Between-cell-line IC50 spread for a drug in GDSC2 BREAST (BRCA) -- the
    breast analogue of ``kill_calibration.gdsc_glioma_ic50_band``. GDSC is
    itself a monolayer assay, so its fitted IC50 IS the 2-D-appropriate
    potency; the between-line IQR is the honest error bar on it. Reads the
    same cache file ``tmesim.molecular.potency_breast.measured_ic50`` reads,
    but keeps the per-line LN_IC50 values (that function only returns the
    aggregated median) so the IQR can be computed."""
    if not os.path.exists(GDSC_BREAST_CSV):
        return None
    aliases = set(GDSC_ALIASES.get(drug_name.strip().lower(), [drug_name.strip().lower()]))
    aliases.add(drug_name.strip().lower())
    per_line: Dict[str, List[float]] = {}
    with open(GDSC_BREAST_CSV, newline="") as fh:
        for row in csv.DictReader(fh):
            if (row.get("DRUG_NAME") or "").strip().lower() not in aliases:
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
    if not per_line:
        return None
    import statistics
    per_cell = [statistics.median(v) for v in per_line.values()]
    v = np.exp(np.asarray(per_cell, dtype=np.float64))
    return {"n_lines": int(v.size), "median_uM": float(np.median(v)),
            "q1_uM": float(np.percentile(v, 25)), "q3_uM": float(np.percentile(v, 75)),
            "min_uM": float(v.min()), "max_uM": float(v.max())}


def build_anchors_breast(res: Dict) -> Dict:
    """Per-drug 2-D potency anchor + uncertainty band, breast panel.

    No breast-specific patient-derived ex-vivo viability dataset exists in
    this repo's cache (unlike GBM's GSE232173 / Ntafoulis et al.), so -- unlike
    temozolomide in the GBM report -- NONE of these four drugs get a
    re-anchored potency. All four use an ALREADY-2-D anchor and report its
    between-line (or quoted-literature) spread as the band, exactly the
    treatment the GBM report gives cisplatin and doxorubicin.
    """
    anchors: Dict[str, Dict] = {}
    for name, d in res["per_drug"].items():
        if name == "doxorubicin":
            anchors[name] = {
                "kind": "literature 2-D range (~0.1-1 uM, not breast- or "
                        "GBM-specific; doxorubicin is not in GDSC2 at all)",
                "ic50_point": d["ic50_uM"], "ic50_lo": 0.1, "ic50_hi": 1.0,
                "c_ref_uM": d["ic50_uM"],
                "band_source": "quoted literature range 0.1-1 uM "
                               "(tmesim.molecular.potency provenance, reused verbatim)",
            }
        else:
            g = gdsc_breast_ic50_band(name)
            anchors[name] = {
                "kind": "GDSC2 breast (BRCA) monolayer median IC50 (already a 2-D anchor)",
                "ic50_point": d["ic50_uM"],
                "ic50_lo": (g["q1_uM"] if g else d["ic50_uM"] * 0.5),
                "ic50_hi": (g["q3_uM"] if g else d["ic50_uM"] * 2.0),
                "c_ref_uM": d["ic50_uM"],
                "band_source": (f"between-line IQR, n={g['n_lines']} GDSC2 breast (BRCA) lines"
                                if g else "assumed 2x"),
                "gdsc": g,
            }
    return anchors


# --------------------------------------------------------------------------- #
# Orchestration -- mirrors kill_calibration.run(), on the breast tissue.
# --------------------------------------------------------------------------- #
def run(dx_um: float = DX_UM, drugs: Sequence[str] = DRUGS, verbose: bool = True,
        verify_engine: bool = True, verify_drug: str = "doxorubicin",
        use_cache: bool = True) -> Dict:
    import json
    t0 = time.time()

    def log(msg: str) -> None:
        if verbose:
            print(f"[kill_calibration_breast] {msg}", flush=True)

    out: Dict = {"dx_um": float(dx_um), "dt_s": DT_S, "window_s": WINDOW_S,
                 "seed": SEED, "drugs": list(drugs)}

    log(f"reconstructing REAL Xenium BREAST section (Rep1) at dx = {dx_um} um ...")
    micro = reconstruct_real(source="xenium", sample=BREAST_SAMPLE_PATH, dx_um=float(dx_um))
    out["reconstruction"] = {k: v for k, v in reconstruction_report(micro).items()
                             if isinstance(v, (int, float, str, bool))}

    pde = PDEConfig(dt=DT_S, total_time=WINDOW_S, convection=True, save_every=SAVE_EVERY)
    pk = PKConfig(plasma_peak=NOMINAL_PEAK_UM, half_life=NOMINAL_HALF_LIFE_S,
                  n_doses=1, dose_interval=3600.0 * 24)
    pd_cfg = PDConfig()  # fresh import -- see the report's own note on what was observed
    base_cfg = SimConfig(pde=pde, pk=pk, pd=pd_cfg, seed=SEED)
    out["pd_config"] = {"hypoxia_resistance": pd_cfg.hypoxia_resistance,
                        "quiescent_resistance": pd_cfg.quiescent_resistance,
                        "max_kill_rate": pd_cfg.max_kill_rate}
    out["nominal_plasma_mean_uM"] = window_plasma_mean(pk, WINDOW_S, DT_S)

    log("engine prologue (oxygen -> cell states -> IFP), once ...")
    work, labeled, velocity = prepare_tissue(micro, base_cfg)
    tumor = labeled.tumor_mask()
    state_t = np.asarray(labeled.state)[tumor]
    out["n_cells"] = int(labeled.n)
    out["n_tumor"] = int(tumor.sum())
    out["state_breakdown"] = state_breakdown(state_t, pd_cfg)

    alpha = calibrate_alpha()
    out["alpha"] = float(alpha)
    log(f"alpha = {alpha:.4f} (fit once on GBM doxorubicin, reused verbatim)")

    cache: Dict = {}
    if use_cache and os.path.exists(RUN_CACHE):
        try:
            cache = json.load(open(RUN_CACHE))
        except Exception:
            cache = {}

    per_drug: Dict[str, Dict] = {}
    for name in drugs:
        log(f"--- {name} ---")
        drug, mol, ic50_info = compile_breast_drug(name, alpha)
        d = {"ic50_uM": float(drug.ic50), "hill": float(drug.hill),
             "bbb_penetrance": float(drug.bbb_penetrance),
             "diffusivity": float(drug.diffusivity),
             "binding": float(drug.binding), "clearance": float(drug.clearance),
             "ic50_source": ic50_info["source"], "ic50_confidence": ic50_info["confidence"],
             "ic50_n_cell_lines": ic50_info.get("n_cell_lines")}

        m2 = apply_drug_sensitivity_on_copy(work, labeled, name)
        sens_uniform = np.asarray(m2.cells.sensitivity, dtype=np.float64)[tumor]
        d["sensitivity_uniform_mean"] = float(sens_uniform.mean())

        ckey = f"{name}|{dx_um}|{DT_S}|{WINDOW_S}|{SAVE_EVERY}|ic50={drug.ic50:.6g}"
        tr = None
        if ckey in cache:
            log("  unit transport: cache hit")
            e_unit = np.asarray(cache[ckey]["e_unit"], dtype=np.float64)
            if "kill_nominal_engine_path" in cache[ckey]:
                d["kill_nominal_engine_path"] = float(cache[ckey]["kill_nominal_engine_path"])
            if "kill_full_engine" in cache[ckey]:
                d["kill_full_engine"] = float(cache[ckey]["kill_full_engine"])
        else:
            log("  unit transport solve ...")
            ts = time.time()
            tr = unit_transport(work, drug, pde, velocity)
            log(f"  ... {time.time()-ts:.1f}s")
            col, row = work.domain.to_index(labeled.x, labeled.y)
            auc = np.asarray(tr.auc, dtype=np.float64)[row, col]
            duration = float(tr.times[-1]) if tr.times.size and tr.times[-1] > 0 else 1.0
            e_unit = (auc / duration)[tumor]
            cache[ckey] = {"e_unit": e_unit.tolist()}

        d["e_unit_p99"] = float(np.percentile(e_unit, 99.0))
        d["e_unit_max"] = float(e_unit.max())
        d["e_unit_median"] = float(np.median(e_unit))
        d["e_unit_mean"] = float(e_unit.mean())

        baths = [m * drug.ic50 for m in BATH_MULTIPLES]
        labels = [f"{m:g}xIC50" for m in BATH_MULTIPLES]
        d["decomposition"] = []
        for lab, cb in zip(labels, baths):
            dec = decompose(e_unit, sens_uniform, resistance_vector(state_t, pd_cfg),
                            cb, drug.ic50, drug.hill, pd_cfg.max_kill_rate)
            dec["label"] = lab
            dec["implied_ic50_fold"] = kill_ratio_to_ic50_fold(
                dec["ratio_exact"], cb, drug.ic50, drug.hill)
            d["decomposition"].append(dec)
        d["headline"] = d["decomposition"][BATH_MULTIPLES.index(1.0)]

        pm = out["nominal_plasma_mean_uM"]
        if tr is not None:
            d["kill_nominal_engine_path"] = kill_with(labeled, work, tr, drug, pd_cfg, pm)
            cache[ckey]["kill_nominal_engine_path"] = d["kill_nominal_engine_path"]

        if verify_engine and name == verify_drug and "kill_full_engine" not in d:
            log("  full-engine verification run (simulate_response) ...")
            from tmesim.core import Microenvironment
            m_v = Microenvironment(domain=micro.domain, fields=micro.fields,
                                   cells=micro.cells.copy(), patient_id=micro.patient_id,
                                   metadata=micro.metadata)
            apply_drug_sensitivity(m_v, name, mode="uniform")
            full = simulate_response(m_v, SimConfig(drug=drug, pde=pde, pk=pk, pd=pd_cfg, seed=SEED))
            d["kill_full_engine"] = float(full.kill_fraction)
            cache[ckey]["kill_full_engine"] = d["kill_full_engine"]
            log(f"  full engine kill = {d['kill_full_engine']:.10f}")

        d["_e_unit"] = e_unit
        d["_sens_uniform"] = sens_uniform
        per_drug[name] = d

    if use_cache:
        try:
            os.makedirs(os.path.dirname(RUN_CACHE), exist_ok=True)
            json.dump(cache, open(RUN_CACHE, "w"))
        except Exception as exc:  # pragma: no cover
            log(f"cache write failed: {exc}")

    out["per_drug"] = per_drug
    out["_work"] = work
    out["_labeled"] = labeled
    out["_tumor"] = tumor
    out["_pd_cfg"] = pd_cfg
    out["elapsed_s"] = time.time() - t0
    log(f"done in {out['elapsed_s']/60:.1f} min")
    return out


def apply_drug_sensitivity_on_copy(work, labeled, name):
    from tmesim.core import Microenvironment
    m2 = Microenvironment(domain=work.domain, fields=work.fields, cells=labeled.copy(),
                          patient_id=work.patient_id, metadata=work.metadata)
    apply_drug_sensitivity(m2, name, mode="uniform")
    return m2


if __name__ == "__main__":
    res = run()
    print(f"\nreconstruction: {res['reconstruction']['n_cells']:,} cells, "
          f"{res['reconstruction']['n_tumor']:,} tumor")
    print(f"alpha = {res['alpha']:.4f}")
    print(f"state breakdown: {res['state_breakdown']}")
    for name, d in res["per_drug"].items():
        print(f"{name:16s} ic50={d['ic50_uM']:10.4g} uM ({d['ic50_source']:>10s})  "
              f"kill_nominal={d['kill_nominal_engine_path']:.4f}  "
              f"headline_3D2D_ratio={d['headline']['ratio_exact']:.4f}")
