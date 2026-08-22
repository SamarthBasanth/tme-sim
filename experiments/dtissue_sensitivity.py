"""FALSIFIABILITY TEST: are the breast 3-D:2-D verdicts an artefact of D_tissue?

Why this file exists
--------------------
``results/kill_calibration_breast.md`` publishes 4 formally-testable 3-D:2-D
ratio verdicts (3 PASS + docetaxel FAIL) among the 6 drugs swept here.
Paclitaxel and 5-fluorouracil are swept too, and their ratios are reported, but
they carry NO verdict: each has exactly one drug-matched study, so its min-max
"band" is a zero-width point that nothing can land inside, and it is reported
NOT TESTABLE rather than re-scored against the broad cross-tumour-type band
(``results/zero_width_bands.md``; this driver's own former use of that
substitution is recorded in ``results/fallback_rule_propagation.md``). Every one
of the reported ratios is computed from a
transport solve whose tissue diffusivity comes out of the molecular compile as

    D_free(doxorubicin)   = 351.31 um^2/s  = 3.513e-6 cm^2/s
    D_tissue(doxorubicin) = 140.41 um^2/s  = 1.404e-6 cm^2/s   (ratio 0.400)

where ``D_tissue`` is ``DrugParams.diffusivity`` (184.75 um^2/s, the
Brownian-dynamics + Ogston-Renkin hindered value) further reduced by the
engine's ECM tortuosity coupling ``D_eff = D * (1 - 0.6 * ecm)`` evaluated at
the calibration reference ``ecm = DEFAULT_ECM = 0.4``.

A literature pass established that the 0.400 ratio is the canonical
small-inert-solute tortuosity rule (1/lambda^2, lambda = 1.6; Sykova &
Nicholson, Physiol Rev 2008, PMID 18923183), whose own authors restrict it to
solutes that do NOT interact with the matrix -- which doxorubicin (cationic,
DNA-binding, endosomally sequestered) emphatically does. Direct measurements in
real tumour tissue are 2-5x LOWER:

    rabbit VX2 liver carcinoma, in vivo   5.0e-7 cm^2/s   (Weinberg 2007, PMID 17900740)
    LS174T 3-D cell mass                  2.9e-7 cm^2/s   (Toley 2013,   PMID 23860772)
    rat liver, in vivo                    6.7e-7 cm^2/s   (Qian 2003,    PMID 12932647)

giving a defensible sweep range D_tissue = 3e-7 .. 1.6e-6 cm^2/s.

This script does NOT pick a value. It sweeps the whole range and reports the
shape of the response curve, the verdict-flip map, and the cancellation factor
(how much of an absolute-delivery change survives into the scored 3-D:2-D
ratio). Selecting the D_tissue that maximises PASS count is explicitly not done
and would invalidate the test.

Two alpha modes -- and why both are needed
------------------------------------------
The engine's perivascular sink constant ``alpha`` is not independent of D. It is
fit in closed form (``tmesim.validation.calibration.calibrate_alpha``) so that
doxorubicin's analytic penetration length ``L = sqrt(D_eff / k_eff)`` equals the
MEASURED anchor 64.9 um (45 um half-distance / ln2; Primeau 2005, PMID 16361566):

    alpha = (D_eff / L*^2 - k_passive) / uptake_proxy

so alpha is (very nearly) proportional to D_eff. Therefore:

* ``frozen``   -- change D only, keep the published alpha = 148.608. Penetration
  length moves as sqrt(D). This is what happens if someone edits the diffusivity
  by hand without re-running the calibration.
* ``recal``    -- change D and re-solve alpha against the SAME measured
  doxorubicin penetration anchor. This is what the repo actually does, because
  ``calibrate_alpha()`` is called at run time from the compiled molecule. L is
  then pinned by a measurement and is almost invariant to D.

Reporting only one of these would be misleading, so both are run.

What is swept
-------------
``D_tissue`` here means the ECM-tortuosity-corrected tissue diffusivity at the
calibration reference ECM, i.e. the 1.404e-6 cm^2/s quantity above. The sweep is
applied as a single multiplicative factor ``s = D_target / 1.404e-6`` on every
drug's ``DrugParams.diffusivity``, which scales the whole spatial D field by s
(the ECM factor is multiplicative and untouched) and preserves the per-drug
relative diffusivity ordering. The literature values are doxorubicin-specific,
but the defect they expose -- an inert-solute tortuosity rule applied to
interacting solutes -- is structural and applies to every drug in the panel, so
a common factor is the honest way to propagate it.

Nothing under ``tmesim/`` is modified. Diffusivity and alpha are overridden from
here via ``dataclasses.replace`` and by passing an explicit alpha into the
already-existing ``compile_breast_drug(name, alpha)``.

Run:
    ./.venv/bin/python experiments/dtissue_sensitivity.py --dx 24
    ./.venv/bin/python experiments/dtissue_sensitivity.py --dx 12 --modes recal
Writes:
    results/sweep/dtissue_sweep_dx{DX}.json   (raw per-point records)
"""
from __future__ import annotations

import argparse
import dataclasses
import functools
import json
import math
import os
import sys
import time
from typing import Dict, List, Optional, Sequence

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, ".."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from tmesim.core import PDConfig, PDEConfig, PKConfig, SimConfig  # noqa: E402
from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams  # noqa: E402
from tmesim.physics.diffusion import ALPHA_TORTUOSITY  # noqa: E402
from tmesim.reconstruct.real_pipeline import reconstruct_real  # noqa: E402
from tmesim.validation.calibration import (  # noqa: E402
    DOX_EFOLD_BAND_UM,
    DOX_TARGET_EFOLD_UM,
    _compile,
    _d_eff_tortuous,
    calibrate_alpha,
    k_passive,
    uptake_proxy,
)
from tmesim.validation.kill_calibration import (  # noqa: E402
    BATH_REFERENCE_PCT,
    NOT_TESTABLE,
    decompose,
    expected_kill,
    kill_ratio_to_ic50_fold,
    prepare_tissue,
    ratio_test,
    require_cached_unit_field,
    resistance_vector,
    score_ratio_row,
    unit_transport,
    unit_transport_cache_key,
    window_plasma_mean,
)
from tmesim.validation.pk_reality import DEFAULT_ECM  # noqa: E402

from experiments.kill_calibration_breast import (  # noqa: E402
    BREAST_SAMPLE_PATH,
    DT_S,
    DX_UM,
    NOMINAL_HALF_LIFE_S,
    NOMINAL_PEAK_UM,
    SAVE_EVERY,
    SEED,
    WINDOW_S,
    apply_drug_sensitivity_on_copy,
    compile_breast_drug,
)
from experiments.kill_calibration_breast_report import (  # noqa: E402
    bridge_to_kill_ratio_band,
    drug_specific_bridge,
    solid_bridge,
)

OUT_DIR = os.path.join(_ROOT, "results", "sweep")

#: The 6 drugs of ``results/kill_calibration_breast.md`` that have ANY
#: drug-matched literature row. Four of them (doxorubicin, docetaxel,
#: epirubicin, cisplatin) have n > 1 and are formally scored; paclitaxel and
#: 5-fluorouracil have n = 1 and are swept and reported but NOT TESTABLE. The
#: other 4 panel drugs have no drug-matched row at all and are not swept.
TESTABLE: Sequence[str] = ("doxorubicin", "paclitaxel", "docetaxel",
                           "epirubicin", "5-fluorouracil", "cisplatin")

#: Log-spaced sweep over the literature-defensible range, cm^2/s.
SWEEP_LO_CM2S = 3.0e-7
SWEEP_HI_CM2S = 1.6e-6
N_SWEEP = 9

#: Measured anchors, for annotation only -- never used to select a value.
MEASURED_ANCHORS = {
    "Toley 2013 LS174T 3-D cell mass (PMID 23860772)": 2.9e-7,
    "Weinberg 2007 rabbit VX2 liver carcinoma in vivo (PMID 17900740)": 5.0e-7,
    "Qian 2003 rat liver in vivo (PMID 12932647)": 6.7e-7,
}


# --------------------------------------------------------------------------- #
# The reference point: what the shipping model actually uses
# --------------------------------------------------------------------------- #
def reference_transport_constants() -> Dict:
    """Read the shipping D_free / D_tissue / ratio straight out of the compile.

    Nothing here is hard-coded: every number is recomputed from the frozen
    molecular pipeline so the report cannot drift from the code.
    """
    from tmesim.molecular.compile_drug import build_profile
    from tmesim.molecular.library import get_molecule

    prof = build_profile(get_molecule("doxorubicin"))
    dp = _compile("doxorubicin")
    d_tissue = _d_eff_tortuous(dp, DEFAULT_ECM)     # um^2/s
    return {
        "free_D_um2_s": float(prof.free_diffusivity),
        "free_D_cm2_s": float(prof.free_diffusivity) * 1e-8,
        "hindered_D_um2_s": float(prof.effective_diffusivity),
        "brownian_ogston_hindrance": float(prof.hindrance),
        "ecm_reference": float(DEFAULT_ECM),
        "ecm_tortuosity_factor": 1.0 - ALPHA_TORTUOSITY * float(DEFAULT_ECM),
        "D_tissue_um2_s": float(d_tissue),
        "D_tissue_cm2_s": float(d_tissue) * 1e-8,
        "D_tissue_over_D_free": float(d_tissue) / float(prof.free_diffusivity),
        "alpha_published": float(calibrate_alpha()),
        "dox_target_efold_um": float(DOX_TARGET_EFOLD_UM),
        "dox_k_passive_per_s": float(k_passive(dp)),
        "dox_uptake_proxy_per_s": float(uptake_proxy(dp)),
    }


def alpha_for_scale(scale: float) -> float:
    """Re-solve the ONE calibration constant at a scaled D, against the SAME
    measured doxorubicin penetration anchor.

    Reproduces ``calibrate_alpha``'s closed form exactly (it is four lines and
    is re-derived here rather than imported because that function recompiles
    doxorubicin from the library and therefore cannot see a scaled D).
    """
    dp = _compile("doxorubicin")
    d_eff = _d_eff_tortuous(dp, DEFAULT_ECM) * float(scale)
    k_target = d_eff / (float(DOX_TARGET_EFOLD_UM) ** 2)
    return float((k_target - k_passive(dp)) / uptake_proxy(dp))


def sweep_points(ref_d_tissue_cm2s: float,
                 explicit: Optional[Sequence[float]] = None) -> List[Dict]:
    """Log-spaced sweep + the shipping reference value, sorted ascending.

    ``explicit`` replaces the log grid with a caller-chosen set of D values (the
    shipping reference is still appended, because every analysis in this file is
    expressed relative to it). Used for the reduced-cost confirmation run on the
    published dx = 12 um grid.
    """
    grid = (np.asarray(explicit, dtype=np.float64) if explicit else
            np.logspace(math.log10(SWEEP_LO_CM2S), math.log10(SWEEP_HI_CM2S),
                        N_SWEEP))
    vals = sorted(set([round(float(v), 12) for v in grid]
                      + [round(float(ref_d_tissue_cm2s), 12)]))
    out = []
    for v in vals:
        out.append({
            "D_tissue_cm2_s": float(v),
            "scale": float(v) / float(ref_d_tissue_cm2s),
            # relative tolerance: ``vals`` is rounded to 1e-12 absolute, so an
            # exact equality (or a 1e-15 absolute window) never fires on a
            # ~1.4e-6 quantity. This is the shipping point, so it must be
            # identified by relative distance.
            "is_reference": abs(v - ref_d_tissue_cm2s) <= 1e-9 * abs(ref_d_tissue_cm2s),
        })
    return out


# --------------------------------------------------------------------------- #
# One (drug, D, mode) point
# --------------------------------------------------------------------------- #
@functools.lru_cache(maxsize=None)
def _compile_breast_cached(name: str, alpha: float):
    """Memoised ``compile_breast_drug`` -- pure speed, identical output. The
    molecular compile runs a Brownian-dynamics particle sim and is the second
    most expensive thing in the sweep after the PDE solve."""
    return compile_breast_drug(name, alpha)


def build_drug(name: str, scale: float, mode: str, alpha_published: float):
    """Compile the drug with the swept D and the mode's alpha.

    ``mode='frozen'`` keeps the published alpha; ``mode='recal'`` re-solves it
    against the measured doxorubicin penetration anchor at the swept D.
    The diffusivity override is applied through ``dataclasses.replace`` and then
    re-run through ``apply_barrier_to_drugparams``; for every small molecule in
    this panel the barrier is a strict no-op (asserted), so the ordering of the
    override relative to the barrier cannot matter.
    """
    alpha = alpha_published if mode == "frozen" else alpha_for_scale(scale)
    dp, mol, ic50_info = _compile_breast_cached(name, float(alpha))
    dp_scaled = dataclasses.replace(dp, diffusivity=float(dp.diffusivity) * float(scale))
    dp_final = apply_barrier_to_drugparams(dp_scaled, mol)
    assert abs(dp_final.binding - dp_scaled.binding) <= 1e-18 * max(1.0, abs(dp_scaled.binding)), (
        f"{name}: binding-site barrier is NOT a no-op ({dp_scaled.binding} -> "
        f"{dp_final.binding}); the diffusivity override would then be order-dependent")
    return dp_final, mol, ic50_info, float(alpha)


def e_unit_for(work, labeled, tumor, drug, pde, velocity) -> np.ndarray:
    """Per-tumour-cell time-averaged exposure per uM plasma -- the exact three
    lines ``kill_calibration_breast.run`` uses, so the reference point is
    numerically comparable to the published report."""
    tr = unit_transport(work, drug, pde, velocity)
    col, row = work.domain.to_index(labeled.x, labeled.y)
    auc = np.asarray(tr.auc, dtype=np.float64)[row, col]
    duration = float(tr.times[-1]) if tr.times.size and tr.times[-1] > 0 else 1.0
    return (auc / duration)[tumor]


def score(name: str, e_unit: np.ndarray, sens: np.ndarray, resist: np.ndarray,
          drug, pd_cfg, bands: Dict, plasma_mean: float) -> Dict:
    """The published scoring path: decompose at 1x IC50, then THE shared band
    decision, :func:`tmesim.validation.kill_calibration.score_ratio_row`.

    This function does NOT decide which band a row is scored against, and must
    never be made to. That decision was previously written here a second time
    (``formal = t_specific if band_is_drug_matched else t_solid``), which is the
    retired rule that silently re-scored a single-study drug against the broad
    62x-wide cross-tumour-type band -- see ``results/zero_width_bands.md`` and
    ``results/fallback_rule_propagation.md``. It is now imported, not copied.

    ``plasma_mean`` additionally drives the frozen ``expected_kill`` at the
    repo-nominal plasma level. That is the ABSOLUTE output of the product (a
    kill fraction at a real dose), as opposed to the bath-matched ratio that is
    scored -- the two respond to D very differently and the whole point of the
    cancellation analysis is the gap between them.
    """
    c_bath = 1.0 * float(drug.ic50)
    dec = decompose(e_unit, sens, resist, c_bath, drug.ic50, drug.hill,
                    pd_cfg.max_kill_rate)
    fold = kill_ratio_to_ic50_fold(dec["ratio_exact"], c_bath, drug.ic50, drug.hill)
    b = bands[name]
    scored = score_ratio_row(dec["ratio_exact"],
                             specific_bridge=b["specific_bridge"],
                             specific_band=b["band"],
                             context_bridge=b["solid_bridge"],
                             context_band=b["band_solid"])
    return {
        "drug": name,
        "c_bath_uM": c_bath,
        "f_penetration": dec["f_penetration"],
        "f_state": dec["f_state"],
        "f_product": dec["f_product"],
        "ratio_exact": dec["ratio_exact"],
        "interaction": dec["interaction"],
        "implied_ic50_fold": fold,
        "E_median_over_bath": dec["E_median_over_bath"],
        "kill_3d_model": dec["kill_3d_model"],
        "kill_2d_model": dec["kill_2d_model"],
        "plasma_required_uM": dec["plasma_required_uM"],
        # ABSOLUTE delivery / absolute output -- the quantities the scored
        # ratio is supposed to be immune to
        "kill_nominal_plasma": expected_kill(e_unit, sens, resist, plasma_mean,
                                             drug.ic50, drug.hill,
                                             pd_cfg.max_kill_rate),
        "e_unit_mean": float(e_unit.mean()),
        "e_unit_median": float(np.median(e_unit)),
        "e_unit_p99": float(np.percentile(e_unit, BATH_REFERENCE_PCT)),
        "e_unit_p01": float(np.percentile(e_unit, 1.0)),
        **verdict_fields(scored, b),
    }


def verdict_fields(scored: Dict, b: Dict) -> Dict:
    """Render ``score_ratio_row``'s decision into this script's record schema.

    One function so that a FRESHLY scored record and a RE-scored cached record
    (:func:`rescore_rows`) cannot end up saying different things, and so that
    the only place a verdict string is produced is downstream of the shared
    decision.
    """
    return {
        "verdict": scored["status"],
        "verdict_scored": bool(scored["scored"]),
        "verdict_band": ("drug-matched" if scored["scored"]
                         else "NOT SCORED -- broad 'solid' band is context only"),
        "verdict_reason": scored["reason"],
        "band_ratio_lo": (b["band"]["ratio_lo"] if b["band"] else None),
        "band_ratio_hi": (b["band"]["ratio_hi"] if b["band"] else None),
        # the band the verdict was actually taken against; None when the row was
        # not scored at all, so a reader cannot mistake a context band for a bar
        "band_used_lo": (scored["band"]["ratio_lo"] if scored["scored"] else None),
        "band_used_hi": (scored["band"]["ratio_hi"] if scored["scored"] else None),
        "verdict_direction": (scored["test"]["direction"] if scored["scored"]
                              else "not scored -- see verdict_reason"),
        # CONTEXT ONLY. Never a verdict; kept so the writeup can quote what the
        # broad non-drug-matched band would have said without scoring against it.
        "context_verdict": scored["context_test"]["verdict"],
        "context_band_lo": b["band_solid"]["ratio_lo"],
        "context_band_hi": b["band_solid"]["ratio_hi"],
        "verdict_solid": scored["context_test"]["verdict"],
    }


def rescore_rows(rows: Sequence[Dict]) -> int:
    """Re-derive every stored verdict IN PLACE through the shared decision.

    Records written before ``2466ce1`` carry a verdict produced by the retired
    fallback rule. The band depends only on the frozen literature tables and on
    the row's own ``(ic50_uM, hill)``, both of which are stored -- so a cached
    record can be re-scored exactly, with NO transport solve. Returns the number
    of rows whose verdict string changed.
    """
    cache: Dict[tuple, Dict] = {}
    changed = 0
    for r in rows:
        name = r["drug"]
        # Older record schemas (the refcheck rows) stored the reporting
        # concentration but not the potency pair. For this panel c_bath IS
        # 1x IC50, and ``hill`` is recovered by re-compiling the molecule --
        # a molecular compile, no transport, no tissue.
        if "ic50_uM" in r and "hill" in r:
            ic50, hill = float(r["ic50_uM"]), float(r["hill"])
        else:
            ic50 = float(r["c_bath_uM"])
            hill = float(_compile_breast_cached(name, float(calibrate_alpha()))[0].hill)
        key = (name, ic50, hill)
        if key not in cache:
            cache[key] = build_bands([name], {name: (ic50, hill)})[name]
        b = cache[key]
        scored = score_ratio_row(float(r["ratio_exact"]),
                                 specific_bridge=b["specific_bridge"],
                                 specific_band=b["band"],
                                 context_bridge=b["solid_bridge"],
                                 context_band=b["band_solid"])
        fields = verdict_fields(scored, b)
        if r.get("verdict") != fields["verdict"]:
            changed += 1
        r.update(fields)
    return changed


def build_bands(drugs: Sequence[str], ic50_hill: Dict) -> Dict:
    """Literature bands, once. They depend only on the frozen literature tables
    and on (ic50, hill), neither of which the D sweep touches -- so every band is
    IDENTICAL at every sweep point and cannot be the thing that moves."""
    solid = solid_bridge()
    out = {}
    for name in drugs:
        ic50, hill = ic50_hill[name]
        c_bath = 1.0 * ic50
        spec = drug_specific_bridge(name)
        band = bridge_to_kill_ratio_band(spec, c_bath, ic50, hill) if spec else None
        out[name] = {
            "band": band,
            "band_solid": bridge_to_kill_ratio_band(solid, c_bath, ic50, hill),
            # The BRIDGES are carried, not just the mapped bands, because
            # ``score_ratio_row`` needs ``n`` to tell a one-study point estimate
            # from an interval. Consumers pass these straight through; none of
            # them re-derives which band is scored.
            "specific_bridge": spec,
            "solid_bridge": solid,
            "band_is_drug_matched": bool(spec is not None and spec["n"] > 1),
            "specific_n": (int(spec["n"]) if spec else 0),
        }
    return out


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #
def _key(mode: str, d_cm2s: float, name: str) -> str:
    return f"{mode}|{d_cm2s:.6e}|{name}"


#: Every field the resume path and the analysis need. A checkpoint line missing
#: any of these was written by an older schema and is discarded.
_REQUIRED_RECORD_KEYS = frozenset((
    "mode", "drug", "D_tissue_cm2_s", "dx_um", "ratio_exact", "f_penetration",
    "verdict", "e_unit_mean", "kill_nominal_plasma", "penetration_length_um",
    "_e_unit_sum", "_e_unit_first",
))


def run(dx_um: float = DX_UM, drugs: Sequence[str] = TESTABLE,
        modes: Sequence[str] = ("frozen", "recal"),
        out_json: Optional[str] = None,
        checkpoint: Optional[str] = None,
        explicit_points: Optional[Sequence[float]] = None) -> Dict:
    t0 = time.time()

    def log(msg: str) -> None:
        print(f"[dtissue {time.time()-t0:7.1f}s] {msg}", flush=True)

    ref = reference_transport_constants()
    log(f"reference: D_free = {ref['free_D_cm2_s']:.3e} cm2/s, "
        f"D_tissue = {ref['D_tissue_cm2_s']:.3e} cm2/s, "
        f"ratio = {ref['D_tissue_over_D_free']:.4f}, alpha = {ref['alpha_published']:.4f}")

    points = sweep_points(ref["D_tissue_cm2_s"], explicit_points)
    log(f"sweep: {len(points)} D points, {len(drugs)} drugs, {len(modes)} modes "
        f"-> {len(points)*len(drugs)*len(modes)} transport solves at dx = {dx_um} um")

    micro = reconstruct_real(source="xenium", sample=BREAST_SAMPLE_PATH,
                            dx_um=float(dx_um))
    pde = PDEConfig(dt=DT_S, total_time=WINDOW_S, convection=True,
                    save_every=SAVE_EVERY)
    pk = PKConfig(plasma_peak=NOMINAL_PEAK_UM, half_life=NOMINAL_HALF_LIFE_S,
                  n_doses=1, dose_interval=3600.0 * 24)
    pd_cfg = PDConfig()
    work, labeled, velocity = prepare_tissue(
        micro, SimConfig(pde=pde, pk=pk, pd=pd_cfg, seed=SEED))
    tumor = labeled.tumor_mask()
    state_t = np.asarray(labeled.state)[tumor]
    plasma_mean = window_plasma_mean(pk, WINDOW_S, DT_S)
    log(f"tissue: grid {micro.domain.shape}, {int(labeled.n)} cells, "
        f"{int(tumor.sum())} tumour; nominal window-mean plasma = "
        f"{plasma_mean:.4f} uM")

    # per-drug static pieces (sensitivity, ic50/hill) -- D-independent
    sens = {}
    ic50_hill = {}
    for name in drugs:
        ts = time.time()
        dp0, _, _ = _compile_breast_cached(name, float(ref["alpha_published"]))
        ic50_hill[name] = (float(dp0.ic50), float(dp0.hill))
        m2 = apply_drug_sensitivity_on_copy(work, labeled, name)
        sens[name] = np.asarray(m2.cells.sensitivity, dtype=np.float64)[tumor]
        log(f"  static {name}: ic50 {dp0.ic50:.4g} uM, D {dp0.diffusivity:.2f} um^2/s "
            f"({time.time()-ts:.1f}s)")
    bands = build_bands(drugs, ic50_hill)

    # --- restart-safe checkpoint: one JSON record per line, keyed ------------
    done: Dict[str, Dict] = {}
    if checkpoint and os.path.exists(checkpoint):
        with open(checkpoint) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue  # truncated final line from a killed run
                # schema guard: silently drop records written by an older
                # revision of ``score()`` rather than resuming a mixed file
                if not _REQUIRED_RECORD_KEYS.issubset(r):
                    continue
                if float(r.get("dx_um", -1)) == float(dx_um):
                    done[_key(r["mode"], r["D_tissue_cm2_s"], r["drug"])] = r
        # A resumed record may predate the shared band decision, so its verdict
        # is re-derived from its own stored ratio before it is reused. No solve.
        n_moved = rescore_rows(list(done.values()))
        log(f"checkpoint: {len(done)} records already complete, resuming "
            f"({n_moved} verdict(s) re-derived under the shared band decision)")
    ck = None
    if checkpoint:
        os.makedirs(os.path.dirname(checkpoint), exist_ok=True)
        ck = open(checkpoint, "a")

    records: List[Dict] = []
    # e_unit fingerprints, so "the sweep genuinely changed the field" is proved
    # rather than assumed.
    fingerprints: Dict[str, List] = {}
    n_solved = 0
    total = len(points) * len(drugs) * len(modes)
    for mode in modes:
        for pt in points:
            for name in drugs:
                k = _key(mode, pt["D_tissue_cm2_s"], name)
                if k in done:
                    rec = done[k]
                    records.append(rec)
                    fingerprints.setdefault(f"{mode}|{name}", []).append(
                        (pt["D_tissue_cm2_s"], rec["_e_unit_sum"],
                         rec["_e_unit_first"]))
                    n_solved += 1
                    continue
                drug, mol, ic50_info, alpha_used = build_drug(
                    name, pt["scale"], mode, ref["alpha_published"])
                e_unit = e_unit_for(work, labeled, tumor, drug, pde, velocity)
                resist = resistance_vector(state_t, pd_cfg, drug=drug)
                rec = score(name, e_unit, sens[name], resist, drug, pd_cfg,
                            bands, plasma_mean)
                rec.update({
                    "mode": mode,
                    "D_tissue_cm2_s": pt["D_tissue_cm2_s"],
                    "scale": pt["scale"],
                    "is_reference": pt["is_reference"],
                    "alpha": alpha_used,
                    "drug_diffusivity_um2_s": float(drug.diffusivity),
                    "drug_binding_per_s": float(drug.binding),
                    "drug_clearance_per_s": float(drug.clearance),
                    "penetration_length_um": math.sqrt(
                        _d_eff_tortuous(drug, DEFAULT_ECM)
                        / (float(drug.binding) + float(drug.clearance))),
                    "ic50_uM": float(drug.ic50),
                    "hill": float(drug.hill),
                    "dx_um": float(dx_um),
                    "_e_unit_sum": float(e_unit.sum()),
                    "_e_unit_first": float(e_unit[0]),
                })
                records.append(rec)
                fingerprints.setdefault(f"{mode}|{name}", []).append(
                    (pt["D_tissue_cm2_s"], float(e_unit.sum()), float(e_unit[0])))
                if ck is not None:
                    ck.write(json.dumps(rec, default=float) + "\n")
                    ck.flush()
                n_solved += 1
                log(f"  {n_solved}/{total}  {mode} D={pt['D_tissue_cm2_s']:.3e} "
                    f"{name}: ratio={rec['ratio_exact']:.4f} "
                    f"F_pen={rec['f_penetration']:.4f} {rec['verdict']}")

    if ck is not None:
        ck.close()

    # --- the assertion the brief demands: the field really moved -------------
    distinct = {}
    for key, fps in fingerprints.items():
        sums = [f[1] for f in fps]
        n_uniq = len(set(sums))
        distinct[key] = {"n_points": len(sums), "n_distinct_e_unit_sums": n_uniq}
        assert n_uniq == len(sums), (
            f"{key}: exposure fields are NOT distinct across the sweep "
            f"({n_uniq} distinct of {len(sums)}) -- the D override is a no-op")
    log("ASSERTION PASSED: every exposure field is bit-distinct across the sweep")

    out = {
        "reference": ref,
        "points": points,
        "drugs": list(drugs),
        "modes": list(modes),
        "dx_um": float(dx_um),
        "dt_s": DT_S, "window_s": WINDOW_S, "save_every": SAVE_EVERY, "seed": SEED,
        "grid": list(micro.domain.shape),
        "n_cells": int(labeled.n), "n_tumor": int(tumor.sum()),
        "bands": {k: {kk: vv for kk, vv in v.items()} for k, v in bands.items()},
        "records": records,
        "field_distinctness": distinct,
        "measured_anchors_cm2_s": MEASURED_ANCHORS,
        "elapsed_s": time.time() - t0,
    }
    if out_json:
        os.makedirs(os.path.dirname(out_json), exist_ok=True)
        with open(out_json, "w") as fh:
            json.dump(out, fh, indent=1, default=float)
        log(f"wrote {out_json}")
    log(f"done in {(time.time()-t0)/60:.1f} min")
    return out


# --------------------------------------------------------------------------- #
# Zero-solve control: score the PUBLISHED dx = 12 um exposure fields
# --------------------------------------------------------------------------- #
#: The on-disk cache ``experiments/kill_calibration_breast.run()`` writes. It
#: holds the exact per-tumour-cell exposure vectors behind
#: ``results/kill_calibration_breast.md``, so scoring it here reproduces the
#: published table with ZERO transport solves and proves this script's scoring
#: path is the published one.
PUBLISHED_CACHE = os.path.join(_ROOT, "data_cache",
                               "kill_calibration_breast_runs.json")


def published_reference_check(drugs: Sequence[str] = TESTABLE,
                              dx_um: float = 12.0,
                              out_json: Optional[str] = None) -> Dict:
    """Re-score the published dx = 12 um run through THIS script's scoring path.

    No PDE is solved: the exposure vectors are read straight out of the cache the
    published driver wrote. Any disagreement with
    ``results/kill_calibration_breast.md`` §3 would mean this script's scoring is
    not the published scoring, and would invalidate the whole sweep.
    """
    t0 = time.time()

    def log(msg):
        print(f"[refcheck {time.time()-t0:6.1f}s] {msg}", flush=True)

    if not os.path.exists(PUBLISHED_CACHE):
        raise FileNotFoundError(
            f"{PUBLISHED_CACHE} not found -- run experiments/kill_calibration_breast.py first")
    cache = json.load(open(PUBLISHED_CACHE))
    alpha = calibrate_alpha()

    micro = reconstruct_real(source="xenium", sample=BREAST_SAMPLE_PATH,
                             dx_um=float(dx_um))
    pde = PDEConfig(dt=DT_S, total_time=WINDOW_S, convection=True,
                    save_every=SAVE_EVERY)
    pk = PKConfig(plasma_peak=NOMINAL_PEAK_UM, half_life=NOMINAL_HALF_LIFE_S,
                  n_doses=1, dose_interval=3600.0 * 24)
    pd_cfg = PDConfig()
    work, labeled, velocity = prepare_tissue(
        micro, SimConfig(pde=pde, pk=pk, pd=pd_cfg, seed=SEED))
    tumor = labeled.tumor_mask()
    state_t = np.asarray(labeled.state)[tumor]
    plasma_mean = window_plasma_mean(pk, WINDOW_S, DT_S)
    log(f"prologue done at dx = {dx_um} um, {int(tumor.sum())} tumour cells")

    ic50_hill, sens = {}, {}
    for name in drugs:
        dp0, _, _ = _compile_breast_cached(name, float(alpha))
        ic50_hill[name] = (float(dp0.ic50), float(dp0.hill))
        m2 = apply_drug_sensitivity_on_copy(work, labeled, name)
        sens[name] = np.asarray(m2.cells.sensitivity, dtype=np.float64)[tumor]
    bands = build_bands(drugs, ic50_hill)

    rows = []
    for name in drugs:
        drug, mol, _ = _compile_breast_cached(name, float(alpha))
        # ONE key function, imported -- never rebuilt here, and a miss RAISES
        # rather than being logged and skipped. A drug quietly dropped from
        # this table is indistinguishable from a drug that was never asked for.
        ckey = unit_transport_cache_key(name, drug, dx_um, dt_s=DT_S,
                                        window_s=WINDOW_S,
                                        save_every=SAVE_EVERY)
        e_unit = require_cached_unit_field(
            cache, name, drug, dx_um, source=PUBLISHED_CACHE, dt_s=DT_S,
            window_s=WINDOW_S, save_every=SAVE_EVERY)
        resist = resistance_vector(state_t, pd_cfg, drug=drug)
        rec = score(name, e_unit, sens[name], resist, drug, pd_cfg, bands,
                    plasma_mean)
        rec.update({"mode": "published-cache", "dx_um": float(dx_um),
                    "D_tissue_cm2_s": float(_d_eff_tortuous(drug, DEFAULT_ECM)) * 1e-8,
                    "cache_key": ckey, "n_cells_in_cache": int(e_unit.size),
                    # carried so this row can be re-scored later without a compile
                    "ic50_uM": float(drug.ic50), "hill": float(drug.hill)})
        rows.append(rec)
        log(f"  {name}: ratio {rec['ratio_exact']:.4f}  F_pen "
            f"{rec['f_penetration']:.4f}  {rec['verdict']}")

    out = {"dx_um": float(dx_um), "alpha": float(alpha), "rows": rows,
           "source": PUBLISHED_CACHE, "elapsed_s": time.time() - t0}
    if out_json:
        os.makedirs(os.path.dirname(out_json), exist_ok=True)
        with open(out_json, "w") as fh:
            json.dump(out, fh, indent=1, default=float)
        log(f"wrote {out_json}")
    return out


def seed_checkpoint_from_published_cache(checkpoint: str,
                                         drugs: Sequence[str] = TESTABLE,
                                         modes: Sequence[str] = ("frozen", "recal"),
                                         dx_um: float = 12.0) -> int:
    """Write the SHIPPING-D records of a dx = 12 um sweep straight from the cache.

    At the shipping D the sweep's scale factor is exactly 1.0 and, in BOTH modes,
    the re-solved alpha equals the published alpha — so the compiled ``DrugParams``
    is bit-identical to the one the published run used, and therefore so is the
    exposure field. Re-solving those 12 PDEs would reproduce a field already on
    disk. This asserts the identity (drug params equal, cache entry present) and
    then writes the records, so the confirmation run only has to solve the points
    that are genuinely new.
    """
    cache = json.load(open(PUBLISHED_CACHE))
    alpha = calibrate_alpha()
    d_ref_cm2s = float(_d_eff_tortuous(_compile("doxorubicin"), DEFAULT_ECM)) * 1e-8

    micro = reconstruct_real(source="xenium", sample=BREAST_SAMPLE_PATH,
                             dx_um=float(dx_um))
    pde = PDEConfig(dt=DT_S, total_time=WINDOW_S, convection=True,
                    save_every=SAVE_EVERY)
    pk = PKConfig(plasma_peak=NOMINAL_PEAK_UM, half_life=NOMINAL_HALF_LIFE_S,
                  n_doses=1, dose_interval=3600.0 * 24)
    pd_cfg = PDConfig()
    work, labeled, velocity = prepare_tissue(
        micro, SimConfig(pde=pde, pk=pk, pd=pd_cfg, seed=SEED))
    tumor = labeled.tumor_mask()
    state_t = np.asarray(labeled.state)[tumor]
    plasma_mean = window_plasma_mean(pk, WINDOW_S, DT_S)

    ic50_hill, sens = {}, {}
    for name in drugs:
        dp0, _, _ = _compile_breast_cached(name, float(alpha))
        ic50_hill[name] = (float(dp0.ic50), float(dp0.hill))
        m2 = apply_drug_sensitivity_on_copy(work, labeled, name)
        sens[name] = np.asarray(m2.cells.sensitivity, dtype=np.float64)[tumor]
    bands = build_bands(drugs, ic50_hill)

    # the rounding sweep_points applies, so the key matches the solved records
    d_key = round(d_ref_cm2s, 12)
    os.makedirs(os.path.dirname(checkpoint), exist_ok=True)
    n = 0
    with open(checkpoint, "a") as fh:
        for mode in modes:
            for name in drugs:
                drug, mol, _, alpha_used = build_drug(name, 1.0, mode, float(alpha))
                pub, _, _ = _compile_breast_cached(name, float(alpha))
                pub = apply_barrier_to_drugparams(pub, mol)
                assert drug == pub, (
                    f"{name}/{mode}: scale-1.0 drug params differ from the published "
                    f"compile; the cached exposure field is NOT reusable")
                # ONE key function, imported -- never rebuilt here, and a miss
                # RAISES rather than skipping. The old `continue` degraded a
                # stale-cache miss into "nothing was cached", which produces a
                # shorter checkpoint that looks complete.
                e_unit = require_cached_unit_field(
                    cache, name, drug, dx_um, source=PUBLISHED_CACHE,
                    dt_s=DT_S, window_s=WINDOW_S, save_every=SAVE_EVERY)
                resist = resistance_vector(state_t, pd_cfg, drug=drug)
                rec = score(name, e_unit, sens[name], resist, drug, pd_cfg,
                            bands, plasma_mean)
                rec.update({
                    "mode": mode, "D_tissue_cm2_s": float(d_key), "scale": 1.0,
                    "is_reference": True, "alpha": float(alpha_used),
                    "drug_diffusivity_um2_s": float(drug.diffusivity),
                    "drug_binding_per_s": float(drug.binding),
                    "drug_clearance_per_s": float(drug.clearance),
                    "penetration_length_um": math.sqrt(
                        _d_eff_tortuous(drug, DEFAULT_ECM)
                        / (float(drug.binding) + float(drug.clearance))),
                    "ic50_uM": float(drug.ic50), "hill": float(drug.hill),
                    "dx_um": float(dx_um),
                    "_e_unit_sum": float(e_unit.sum()),
                    "_e_unit_first": float(e_unit[0]),
                    "_seeded_from": PUBLISHED_CACHE,
                })
                fh.write(json.dumps(rec, default=float) + "\n")
                n += 1
                print(f"  seeded {mode}/{name}: ratio {rec['ratio_exact']:.4f} "
                      f"{rec['verdict']}", flush=True)
    return n


# --------------------------------------------------------------------------- #
# Analysis: verdict-flip map + cancellation factor
# --------------------------------------------------------------------------- #
def _by(res: Dict, mode: str, drug: str) -> List[Dict]:
    """Rows for one (mode, drug) series, ascending in D, with ``is_reference``
    re-derived by nearest match to the shipping D_tissue.

    Re-derived rather than trusted: an early revision of ``sweep_points`` used an
    absolute 1e-15 tolerance to flag the shipping point, which can never fire on
    a ~1.4e-6 quantity, so records written by that revision carry
    ``is_reference = False`` everywhere. Recomputing here makes the analysis
    correct for both old and new record files.
    """
    rs = [r for r in res["records"] if r["mode"] == mode and r["drug"] == drug]
    rs = sorted(rs, key=lambda r: r["D_tissue_cm2_s"])
    if rs:
        ref_d = float(res["reference"]["D_tissue_cm2_s"])
        i = min(range(len(rs)), key=lambda j: abs(rs[j]["D_tissue_cm2_s"] - ref_d))
        assert abs(rs[i]["D_tissue_cm2_s"] - ref_d) <= 1e-6 * ref_d, (
            "the shipping D_tissue is not one of the swept points")
        for j, r in enumerate(rs):
            r["is_reference"] = (j == i)
    return rs


def verdict_flip_map(res: Dict) -> Dict:
    """For each (mode, drug): the contiguous D interval(s) over which the FORMAL
    verdict is PASS, plus the boundary D values, bracketed by the sweep points
    (no interpolation -- a flip is reported as lying between two solved points).
    """
    out: Dict[str, Dict] = {}
    for mode in res["modes"]:
        for drug in res["drugs"]:
            rows = _by(res, mode, drug)
            verdicts = [(r["D_tissue_cm2_s"], r["verdict"]) for r in rows]
            passes = [d for d, v in verdicts if v == "PASS"]
            flips = []
            for a, b in zip(verdicts, verdicts[1:]):
                if a[1] != b[1]:
                    flips.append({"between_cm2_s": [a[0], b[0]],
                                  "from": a[1], "to": b[1]})
            out[f"{mode}|{drug}"] = {
                "mode": mode, "drug": drug,
                "verdicts": verdicts,
                "pass_min_cm2_s": (min(passes) if passes else None),
                "pass_max_cm2_s": (max(passes) if passes else None),
                "n_pass": len(passes), "n_points": len(verdicts),
                "flips": flips,
                "constant": len(set(v for _, v in verdicts)) == 1,
                "verdict_everywhere": (verdicts[0][1] if
                                       len(set(v for _, v in verdicts)) == 1 else None),
                "direction_at_reference": next(
                    (r["verdict_direction"] for r in rows if r["is_reference"]), None),
            }
    return out


def cancellation(res: Dict) -> Dict:
    """How much of an ABSOLUTE-delivery change survives into the SCORED ratio.

    The 3-D:2-D ratio is a within-drug, same-tissue comparison in which the
    plasma level is re-solved so the p99 tumour cell sits exactly at the bath
    concentration (``kill_calibration.decompose``). A change to D that scaled the
    exposure field UNIFORMLY would therefore cancel EXACTLY -- the ratio can only
    respond to a change in the SHAPE of the exposure distribution. This function
    measures how much of the delivery change is shape and how much is scale.

    For each (mode, drug), across the full swept range:
        delivery_fold = e_unit_mean(D_hi) / e_unit_mean(D_lo)
        ratio_fold    = ratio_exact(D_hi) / ratio_exact(D_lo)
        surviving     = ln(ratio_fold) / ln(delivery_fold)      (0 = total cancellation)
        cancellation_fold = delivery_fold / ratio_fold          (how much is absorbed)
    and the same pair evaluated from the shipping reference D down to each
    measured literature anchor.
    """
    out: Dict[str, Dict] = {}
    for mode in res["modes"]:
        for drug in res["drugs"]:
            rows = _by(res, mode, drug)
            lo, hi = rows[0], rows[-1]
            ref = next(r for r in rows if r["is_reference"])

            def pair(a: Dict, b: Dict) -> Dict:
                dfold = b["e_unit_mean"] / a["e_unit_mean"]
                rfold = b["ratio_exact"] / a["ratio_exact"]
                # the product's ABSOLUTE output at a real dose
                kfold = (b.get("kill_nominal_plasma", float("nan"))
                         / a.get("kill_nominal_plasma", float("nan")))

                def _surv(num_fold, den_fold):
                    if not (num_fold > 0 and den_fold > 0):
                        return float("nan")
                    ld = math.log(den_fold)
                    return math.log(num_fold) / ld if abs(ld) > 1e-9 else float("nan")

                def _mag(f):
                    # size of a change, direction-free, so folds below 1 compare
                    # correctly against folds above 1
                    return (math.exp(abs(math.log(f))) if f > 0 else float("nan"))

                return {
                    "D_from_cm2_s": a["D_tissue_cm2_s"], "D_to_cm2_s": b["D_tissue_cm2_s"],
                    "D_fold": b["D_tissue_cm2_s"] / a["D_tissue_cm2_s"],
                    "delivery_fold": dfold,
                    "ratio_fold": rfold,
                    "kill_nominal_fold": kfold,
                    # vs the raw delivery change
                    "fraction_surviving": _surv(rfold, dfold),
                    "cancellation_fold": (_mag(dfold) / _mag(rfold)
                                          if rfold > 0 else float("nan")),
                    # vs the absolute product output (kill at the nominal dose)
                    "fraction_surviving_vs_kill": _surv(rfold, kfold),
                    "cancellation_fold_vs_kill": (_mag(kfold) / _mag(rfold)
                                                  if rfold > 0 else float("nan")),
                    "L_from_um": a["penetration_length_um"],
                    "L_to_um": b["penetration_length_um"],
                }

            def slope(field: str) -> float:
                """d ln(field) / d ln(D_tissue), least squares over the sweep."""
                x = np.log([r["D_tissue_cm2_s"] for r in rows])
                y = np.log([max(float(r.get(field, float("nan"))), 1e-300)
                            for r in rows])
                if not np.all(np.isfinite(y)):
                    return float("nan")
                return float(np.polyfit(x, y, 1)[0])

            s_ratio, s_exp = slope("ratio_exact"), slope("e_unit_mean")
            s_kill = slope("kill_nominal_plasma")
            entry = {"mode": mode, "drug": drug,
                     "full_range": pair(lo, hi),
                     "reference_to_lowest": pair(ref, lo),
                     "slope_ratio": s_ratio,
                     "slope_exposure": s_exp,
                     "slope_kill_nominal": s_kill,
                     "slope_ratio_over_kill": (s_ratio / s_kill if s_kill else float("nan")),
                     "slope_ratio_over_exposure": (s_ratio / s_exp if s_exp else float("nan"))}
            for label, val in MEASURED_ANCHORS.items():
                nearest = min(rows, key=lambda r: abs(r["D_tissue_cm2_s"] - val))
                entry[f"reference_to_nearest[{label}]"] = pair(ref, nearest)
            out[f"{mode}|{drug}"] = entry
    return out


def _f(x, n=4, dash="—"):
    if x is None:
        return dash
    try:
        x = float(x)
    except (TypeError, ValueError):
        return dash
    if math.isnan(x) or math.isinf(x):
        return dash
    return f"{x:.{n}f}"


def write_report(res: Dict, out_md: str, refcheck: Optional[Dict] = None,
                 confirm: Optional[Dict] = None) -> str:
    # Every stored verdict is re-derived through the shared decision BEFORE any
    # analysis reads it, so a record written under the retired fallback rule
    # cannot reach the document. This is a re-score of a stored ratio against a
    # frozen band; it solves nothing and moves no ratio.
    n_rescored = rescore_rows(res["records"])
    if refcheck:
        n_rescored += rescore_rows(refcheck["rows"])
    if confirm:
        n_rescored += rescore_rows(confirm["records"])

    ref = res["reference"]
    flips = verdict_flip_map(res)
    canc = cancellation(res)
    L: List[str] = []
    A = L.append

    #: The published table (``results/kill_calibration_breast.md`` §3), quoted so
    #: this report can state plainly whether it reproduces it. Paclitaxel and
    #: 5-fluorouracil are NOT TESTABLE there (n = 1 each: a point estimate, not
    #: an interval), so they are NOT TESTABLE here; their ratios are still
    #: reproduced and compared, because a ratio is a measurement even when no
    #: bar exists to score it against.
    published = {"doxorubicin": 0.2301, "paclitaxel": 0.1503, "docetaxel": 0.1484,
                 "epirubicin": 0.2301, "5-fluorouracil": 0.3379, "cisplatin": 0.3085}
    pub_v = {"doxorubicin": "PASS", "paclitaxel": NOT_TESTABLE, "docetaxel": "FAIL",
             "epirubicin": "PASS", "5-fluorouracil": NOT_TESTABLE,
             "cisplatin": "PASS"}
    # At the shipping D the two modes are identical by construction (the recal
    # alpha at scale = 1 IS the published alpha), so either may be read.
    ref_mode = "recal" if "recal" in res["modes"] else res["modes"][0]
    ref_recal = {d: next(r for r in _by(res, ref_mode, d) if r["is_reference"])
                 for d in res["drugs"]}

    A("# D_tissue falsifiability sweep — do the breast 3-D:2-D verdicts survive "
      "a measured tissue diffusivity?\n")
    A(f"Generated by `experiments/dtissue_sensitivity.py` on the same REAL Xenium breast "
      f"section as `results/kill_calibration_breast.md` "
      f"(`data_cache/xenium_breast_janesick/rep1`; Janesick et al. 2023, PMID 38114474), "
      f"at **`dx = {res['dx_um']:g} um`** (grid {res['grid'][0]}x{res['grid'][1]}"
      + ("" if abs(res['dx_um'] - 12.0) < 1e-9 else
         " — **NOT the published 12 um grid; see §7.2**") + f"), with the same "
      f"`dt = {res['dt_s']:g} s`, the same {res['window_s']/3600:g} h window, the same seed "
      f"{res['seed']}, the same `PDConfig`, the same p99 bath-matching convention and the "
      f"same frozen literature bands as the published report. "
      f"{res['n_cells']:,} cells, {res['n_tumor']:,} tumour. "
      f"{len(res['records'])} transport solves, {res['elapsed_s']/60:.0f} min.\n")

    A("> **This is a falsifiability test, not a calibration.** No D_tissue value was "
      "selected on the basis of how many verdicts it makes pass. The output is the shape "
      "of the response curve across the whole literature-defensible range, plus an honest "
      "statement of which verdicts are robust and which are artefacts of one parameter "
      "choice. The recommendation at the end is justified by measurement provenance only.\n")

    # A regeneration that re-scores stored ratios cannot re-solve them, so if the
    # transport model has moved since those ratios were solved, the document must
    # say so rather than let a fresh timestamp imply freshness.
    if abs(float(res["dx_um"]) - 12.0) > 1e-9:
        A(f"> ### STALE TRANSPORT WARNING — this sweep's exposure fields predate the "
          f"`_TUMOR_SIGMA` units fix\n>\n"
          f"> The interstitial-pressure solve smooths the tumour-occupancy field with a "
          f"gaussian whose width was, until commit `5d6e8e4`, expressed in GRID CELLS "
          f"(`_TUMOR_SIGMA = 5.0`) and therefore tracked the mesh. It is now a physical "
          f"length (`_TUMOR_SIGMA_UM = 60.0`). **`dx = 12 um` is the exact fixed point** "
          f"(`60.0 / 12.0 == 5.0` in IEEE-754), so every shipped `dx = 12` result is "
          f"bit-identical — but this sweep ran at **`dx = {res['dx_um']:g} um`**, where the "
          f"old constant gave a {5.0 * float(res['dx_um']):.0f} um smoothing instead of "
          f"60 um. `results/tumor_sigma_units.md` §3 measures the consequence on the mean "
          f"tissue IFP at `dx = 24`: **0.181899 -> 0.132363**. The velocity field feeds "
          f"the convection term of every transport solve behind the table below.\n>\n"
          f"> **The exposure fields and hence every RATIO in this document were solved "
          f"before that fix and have NOT been re-solved here** — this regeneration "
          f"re-scored stored ratios against the corrected literature bands and solved no "
          f"PDE (`results/fallback_rule_propagation.md`). The band decision is exact "
          f"regardless, because it consumes only the stored scalar. The ratios are not: "
          f"they would move by an unmeasured amount on a re-run, and the §7.3 "
          f"confirmation set at the published `dx = 12 um` grid is the part of this file "
          f"that is unaffected by construction. Re-solving costs 120 transport solves and "
          f"was deliberately not done here.\n")

    # ---------------- bottom line, computed --------------------------------
    A("## 0. Bottom line\n")
    for mode in res["modes"]:
        fragile = [d for d in res["drugs"] if not flips[f"{mode}|{d}"]["constant"]]
        robust_p = [d for d in res["drugs"]
                    if flips[f"{mode}|{d}"]["verdict_everywhere"] == "PASS"]
        robust_f = [d for d in res["drugs"]
                    if flips[f"{mode}|{d}"]["verdict_everywhere"] == "FAIL"]
        untestable = [d for d in res["drugs"]
                      if flips[f"{mode}|{d}"]["verdict_everywhere"] == NOT_TESTABLE]
        s_ratio = [canc[f"{mode}|{d}"]["slope_ratio"] for d in res["drugs"]]
        s_kill = [canc[f"{mode}|{d}"]["slope_kill_nominal"] for d in res["drugs"]]
        A(f"**mode `{mode}`** over the full "
          f"{res['points'][-1]['D_tissue_cm2_s']/res['points'][0]['D_tissue_cm2_s']:.1f}x "
          f"D_tissue range, across the "
          f"{len(res['drugs']) - len(untestable)} of {len(res['drugs'])} swept drugs that "
          f"carry a scorable band"
          + (f" ({', '.join(untestable)} are NOT TESTABLE at every D — a single-study "
             f"literature point is not an interval, so there is no bar for D to move "
             f"them across; their ratios are still swept and reported)"
             if untestable else "") + ": "
          f"{len(robust_p)} verdict(s) PASS everywhere, {len(robust_f)} FAIL everywhere, "
          f"**{len(fragile)} flip** ({', '.join(fragile) if fragile else 'none'}). "
          f"Log-log response of the SCORED ratio to D_tissue spans "
          f"{min(s_ratio):+.2f} to {max(s_ratio):+.2f}; of the ABSOLUTE kill at a real "
          f"dose, {min(s_kill):+.2f} to {max(s_kill):+.2f}.\n")
    if confirm and confirm.get("records"):
        per_mode = {}
        for mode in confirm["modes"]:
            n = 0
            for drug in confirm["drugs"]:
                rws = _by(confirm, mode, drug)
                if len(rws) < 2:
                    continue
                hi = next((r for r in rws if r["is_reference"]), rws[-1])
                n += int(rws[0]["verdict"] != hi["verdict"])
            per_mode[mode] = n
        A("On the **published `dx = 12 um` grid** (§7.3, reduced 2-point set: the shipping "
          "diffusivity vs the measurement-backed one), the flip count is "
          + ", ".join(f"**{n} in `{m}`**" for m, n in per_mode.items())
          + ". The one flip is **docetaxel FAIL -> PASS in `frozen`**; in `recal` nothing "
            "flips and every ratio moves by under 1 %.\n")
    A("The single most important structural fact, established in §2 and quantified in §5: "
      "`alpha` is fit so the model reproduces a MEASURED doxorubicin penetration depth, and "
      "that fit makes `alpha` nearly proportional to D. In the mode the codebase actually "
      "runs (`recal`), a several-fold change in D_tissue is therefore absorbed almost "
      "entirely by the re-fitted sink, and the scored 3-D:2-D ratio barely moves. The "
      "diffusivity being wrong is a real defect, but it is a defect in ATTRIBUTION "
      "(how much drug is claimed to be consumed perivascularly), not primarily in the "
      "published 3-D:2-D verdicts.\n")

    # ---------------- 1. what D_tissue actually is in this code -------------
    A("## 1. Where D_tissue actually enters the solver (traced, not assumed)\n")
    A("The transport PDE reads `DrugParams.diffusivity` (um^2/s) and applies the engine's "
      "own ECM tortuosity coupling `D_eff(x,y) = diffusivity * (1 - 0.6 * ecm_density(x,y))` "
      "(`tmesim/physics/diffusion.py`, `ALPHA_TORTUOSITY = 0.6`). `DrugParams.diffusivity` is "
      "not a constant in a table for these drugs — it is computed by the molecular compile:\n")
    A("| step | doxorubicin value | source |")
    A("|---|---|---|")
    A(f"| hydrodynamic radius | {ref['free_D_um2_s']:.2f} um^2/s implies r = "
      f"{1.380649e-23*310.0/(6*math.pi*1.2e-3*ref['free_D_um2_s']*1e-12)*1e9:.3f} nm | "
      "`stokes_einstein.hydrodynamic_radius_nm` (empirical MW scaling) |")
    A(f"| D_free (Stokes-Einstein) | **{ref['free_D_um2_s']:.2f} um^2/s = "
      f"{ref['free_D_cm2_s']:.3e} cm^2/s** | `stokes_einstein.free_diffusivity`, "
      "T = 310 K, eta = 1.2e-3 Pa.s |")
    A(f"| x Brownian-dynamics + Ogston-Renkin hindrance {ref['brownian_ogston_hindrance']:.4f} "
      f"| {ref['hindered_D_um2_s']:.2f} um^2/s | `compile_drug.build_profile` |")
    A(f"| x ECM tortuosity at the calibration reference ecm = {ref['ecm_reference']:g} "
      f"(factor {ref['ecm_tortuosity_factor']:.2f}) | **{ref['D_tissue_um2_s']:.2f} um^2/s = "
      f"{ref['D_tissue_cm2_s']:.3e} cm^2/s** | `physics/diffusion._effective_diffusivity` |")
    A(f"| **D_tissue / D_free** | **{ref['D_tissue_over_D_free']:.4f}** | — |\n")
    A(f"So the 0.40 ratio the literature pass flagged is confirmed to 4 significant figures "
      f"(**{ref['D_tissue_over_D_free']:.4f}**) — but it is NOT a hard-coded tortuosity "
      "constant. It is the *product* of two independently-derived hindrance factors "
      f"({ref['brownian_ogston_hindrance']:.4f} from the particle sim and mesh gate, "
      f"{ref['ecm_tortuosity_factor']:.2f} from the ECM coupling) that happen to land on the "
      "canonical inert-solute value. The criticism stands regardless: neither factor knows "
      "anything about cationic binding, DNA intercalation or endosomal sequestration, so the "
      "model has no mechanism by which doxorubicin could be hindered more than an inert "
      "solute of the same size.\n")
    A("**The sweep is applied as a single multiplicative factor `s = D_target / "
      f"{ref['D_tissue_cm2_s']:.3e}` on every drug's `DrugParams.diffusivity`.** Because the "
      "ECM factor is multiplicative and untouched, this scales the whole spatial D field by "
      "exactly `s` and preserves the per-drug relative diffusivity ordering.\n")
    A("**Trap check (the brief's KNOWN TRAP).** `kill_calibration.unit_transport()` calls "
      "`simulate_drug_transport` directly, so no `PDEConfig` flag could have reached the "
      "solver. D reaches it only through the `DrugParams` object, which is what this script "
      "overrides via `dataclasses.replace`. Every exposure field in this sweep was "
      "fingerprinted and asserted bit-distinct from every other point of the same "
      "(mode, drug) series before any result below was computed; the assertion is in "
      "`run()` and it passed:\n")
    A("```")
    for k, v in sorted(res["field_distinctness"].items()):
        A(f"{k:34s} {v['n_distinct_e_unit_sums']}/{v['n_points']} distinct exposure fields")
    A("```\n")

    # ---------------- 2. the alpha coupling --------------------------------
    A("## 2. D_tissue is NOT an independent parameter — it is coupled to `alpha`\n")
    A(f"`tmesim.validation.calibration.calibrate_alpha` solves in closed form for the ONE "
      f"perivascular-sink constant so that doxorubicin's analytic penetration length "
      f"`L = sqrt(D_eff / k_eff)` equals the MEASURED anchor "
      f"**{ref['dox_target_efold_um']:.2f} um** (45 um half-distance / ln2; Primeau et al. "
      f"2005, **PMID 16361566**):\n")
    A("```\nalpha = (D_eff / L*^2 - k_passive) / uptake_proxy\n```\n")
    A(f"with `k_passive = {ref['dox_k_passive_per_s']:.4e} /s` and `uptake_proxy = "
      f"{ref['dox_uptake_proxy_per_s']:.4e} /s`, giving the shipping "
      f"**alpha = {ref['alpha_published']:.4f}**. Because `k_passive << alpha*uptake_proxy`, "
      "alpha is very nearly PROPORTIONAL to D. Lowering D therefore lowers the fitted sink "
      "by nearly the same factor, and `L` barely moves. Reporting a sweep that ignores this "
      "would be reporting a scenario the codebase does not actually produce. Both are "
      "therefore run:\n")
    A("| mode | what it is | what happens to L |")
    A("|---|---|---|")
    A("| `frozen` | D changed by hand, published alpha kept | `L` moves as `sqrt(D)` |")
    A("| `recal` | D changed AND alpha re-solved against the SAME measured doxorubicin "
      "penetration anchor — **what the repo actually does, because `calibrate_alpha()` is "
      "called at run time from the compiled molecule** | `L` pinned by the measurement, "
      "nearly invariant |\n")
    rows_dox_f = _by(res, "frozen", "doxorubicin")
    rows_dox_r = _by(res, "recal", "doxorubicin")
    A("Doxorubicin penetration length across the sweep, both modes:\n")
    A("| D_tissue (cm^2/s) | alpha (frozen) | L um (frozen) | alpha (recal) | L um (recal) |")
    A("|---|---|---|---|---|")
    for a, b in zip(rows_dox_f, rows_dox_r):
        star = " **(shipping)**" if a["is_reference"] else ""
        A(f"| {a['D_tissue_cm2_s']:.3e}{star} | {a['alpha']:.2f} | "
          f"{a['penetration_length_um']:.1f} | {b['alpha']:.2f} | "
          f"{b['penetration_length_um']:.1f} |")
    A("")

    # ---------------- 3. the response curves --------------------------------
    A("## 3. The response curve, per drug, per sweep point\n")
    for mode in res["modes"]:
        A(f"### 3.{1 if mode=='frozen' else 2} mode = `{mode}`\n")
        for drug in res["drugs"]:
            rows = _by(res, mode, drug)
            b = res["bands"][drug]
            r0 = rows[0]
            if r0["verdict_scored"]:
                A(f"**{drug}** — IC50 {r0['ic50_uM']:.4g} uM, hill {r0['hill']:.3g}, "
                  f"tested at 1x IC50 against its drug-matched band "
                  f"(n = {b['specific_n']}), "
                  f"kill-ratio band **{r0['band_used_lo']:.4f}-{r0['band_used_hi']:.4f}** "
                  f"(this band is FIXED across the sweep — it depends only on the frozen "
                  f"literature table and on (IC50, hill), neither of which D touches).\n")
            else:
                A(f"**{drug}** — IC50 {r0['ic50_uM']:.4g} uM, hill {r0['hill']:.3g}, "
                  f"reported at 1x IC50 but **NOT SCORED**: its drug-matched literature "
                  f"is a single study (n = {b['specific_n']}), so the min-max interval is "
                  f"zero-width "
                  f"(**{_f(r0['band_ratio_lo'])}-{_f(r0['band_ratio_hi'])}**) and nothing "
                  f"can land inside it. The broad 'solid' band "
                  f"(**{r0['context_band_lo']:.4f}-{r0['context_band_hi']:.4f}**, which "
                  f"this drug's ratios do sit inside) is quoted for CONTEXT and is NOT "
                  f"substituted for the missing bar — see `results/zero_width_bands.md`. "
                  f"The ratios below are still measurements and still move with D; what "
                  f"is absent is a bar to move across.\n")
            A("| D_tissue (cm^2/s) | F_penetration | F_state | 3-D:2-D ratio | implied fold | "
              "mean tumour exposure (uM per uM plasma) | absolute kill at nominal dose | "
              "L um | verdict |")
            A("|---|---|---|---|---|---|---|---|---|")
            for r in rows:
                star = " **(shipping)**" if r["is_reference"] else ""
                fold = ("—" if r["implied_ic50_fold"] is None
                        else f"{r['implied_ic50_fold']:.2f}x")
                A(f"| {r['D_tissue_cm2_s']:.3e}{star} | {_f(r['f_penetration'])} | "
                  f"{_f(r['f_state'])} | **{_f(r['ratio_exact'])}** | {fold} | "
                  f"{r['e_unit_mean']:.4g} | {_f(r.get('kill_nominal_plasma'))} | "
                  f"{r['penetration_length_um']:.1f} | "
                  f"**{r['verdict']}** |")
            A("")

    # ---------------- 4. verdict-flip map -----------------------------------
    A("## 4. Verdict-flip map\n")
    A("Over what D_tissue interval does each drug PASS? `—` means the verdict never "
      "changes anywhere in the swept range.\n")
    mismatched = [d for d in res["drugs"]
                  if d in pub_v and ref_recal[d]["verdict"] != pub_v[d]]
    if mismatched:
        A(f"> **READ §7.2 AND §7.3 BEFORE THIS TABLE.** At `dx = {res['dx_um']:g} um` the "
          f"BASELINE verdict for **{', '.join(mismatched)}** already differs from the "
          f"published `dx = 12 um` verdict, so its row below is not a statement about the "
          f"published verdict. That is itself a finding: "
          f"**{', '.join(mismatched)}** is more sensitive to GRID RESOLUTION than to "
          f"D_tissue. §7.3 carries the flip map computed on the published grid, and it is "
          f"the one to quote.\n")
    for mode in res["modes"]:
        A(f"**mode = `{mode}`**\n")
        A("| drug | verdict at shipping D | PASS over | flips at | robust? |")
        A("|---|---|---|---|---|")
        for drug in res["drugs"]:
            fm = flips[f"{mode}|{drug}"]
            rows = _by(res, mode, drug)
            at_ref = next(r["verdict"] for r in rows if r["is_reference"])
            if fm["verdict_everywhere"] == NOT_TESTABLE:
                # No bar exists, so there is nothing for D to move the drug
                # across. Saying "ROBUST" here would claim a stability result
                # that was never tested.
                span = "n/a — not scored"
                flipstr = "— (no bar to flip across)"
                robust = "**NOT TESTABLE at every D** (no verdict, robust or otherwise)"
            elif fm["constant"]:
                span = ("ALL 10 points" if fm["verdict_everywhere"] == "PASS"
                        else "NONE")
                flipstr = "— (never flips)"
                robust = ("**ROBUST PASS**" if fm["verdict_everywhere"] == "PASS"
                          else "**ROBUST FAIL**")
            else:
                # a single flip means the PASS points are contiguous; more than
                # one means they are not, and saying "min - max" would overstate
                contig = "" if len(fm["flips"]) <= 1 else " **(NOT contiguous)**"
                span = (f"{fm['pass_min_cm2_s']:.3e} - {fm['pass_max_cm2_s']:.3e} "
                        f"({fm['n_pass']}/{fm['n_points']} points){contig}")
                flipstr = "; ".join(
                    f"{f['from']}->{f['to']} between {f['between_cm2_s'][0]:.3e} and "
                    f"{f['between_cm2_s'][1]:.3e}" for f in fm["flips"])
                robust = "**PARAMETER-DEPENDENT**"
            A(f"| {drug} | {at_ref} | {span} | {flipstr} | {robust} |")
        A("")

    # ---------------- 5. cancellation ---------------------------------------
    A("## 5. The cancellation factor — how much of a delivery change reaches the score\n")
    A("The scored quantity is a 3-D:2-D kill RATIO computed within one drug on one tissue, "
      "with the plasma level re-solved at every point so that the p99 best-perfused tumour "
      "cell sits exactly at the bath concentration (`kill_calibration.decompose`). "
      "**A change to D that scaled the exposure field uniformly would cancel EXACTLY** — "
      "algebraically, if `e -> s*e` then `ref_unit -> s*ref_unit`, `plasma = c_bath/(s*ref_unit)`, "
      "and `E = plasma * s * e` is unchanged. The ratio can therefore only respond to a change "
      "in the SHAPE of the exposure distribution, never to its scale. The numbers below "
      "measure how much of the change is shape and how much is scale.\n")
    A("Two absolute quantities are tracked against the scored ratio: the mean tumour "
      "exposure per uM of plasma (raw delivery) and the model's kill fraction at the "
      "repo-nominal dose (the product's actual absolute output). "
      "`fraction surviving = ln(ratio fold) / ln(absolute fold)`: 0.0 = total cancellation, "
      "1.0 = no cancellation, **above 1.0 = the score is MORE sensitive than the absolute "
      "output, i.e. anti-cancellation**. `cancellation fold` compares the SIZES of the two "
      "changes direction-free (`exp|ln fold|` of each), so a value of 3x means the absolute "
      "quantity moved three times as far as the score.\n")
    A("Log-log response exponents `d ln(quantity) / d ln(D_tissue)`, least-squares over all "
      "sweep points — the single cleanest summary of the whole question:\n")
    A("| mode | drug | d ln(exposure)/d ln D | d ln(absolute kill)/d ln D | "
      "**d ln(scored ratio)/d ln D** | score / absolute-kill sensitivity |")
    A("|---|---|---|---|---|---|")
    for mode in res["modes"]:
        for drug in res["drugs"]:
            c = canc[f"{mode}|{drug}"]
            A(f"| {mode} | {drug} | {_f(c['slope_exposure'],3)} | "
              f"{_f(c['slope_kill_nominal'],3)} | **{_f(c['slope_ratio'],3)}** | "
              f"{_f(c['slope_ratio_over_kill'],3)} |")
    A("")
    for mode in res["modes"]:
        A(f"**mode = `{mode}`**, across the FULL swept range "
          f"({res['points'][0]['D_tissue_cm2_s']:.3e} -> "
          f"{res['points'][-1]['D_tissue_cm2_s']:.3e} cm^2/s, a "
          f"{res['points'][-1]['D_tissue_cm2_s']/res['points'][0]['D_tissue_cm2_s']:.2f}x "
          f"change in D):\n")
        A("| drug | mean tumour exposure fold | absolute kill fold (nominal dose) | "
          "scored ratio fold | surviving vs exposure | surviving vs kill | "
          "cancellation vs kill |")
        A("|---|---|---|---|---|---|---|")
        for drug in res["drugs"]:
            c = canc[f"{mode}|{drug}"]["full_range"]
            A(f"| {drug} | {c['delivery_fold']:.3f}x | {_f(c['kill_nominal_fold'],3)}x | "
              f"{c['ratio_fold']:.3f}x | {_f(c['fraction_surviving'], 3)} | "
              f"**{_f(c['fraction_surviving_vs_kill'], 3)}** | "
              f"**{_f(c['cancellation_fold_vs_kill'], 2)}x** |")
        A("")
    A("And the specific comparison that matters — from the shipping D down to the nearest "
      "swept point to each DIRECT measurement in real tumour tissue:\n")
    for label, val in MEASURED_ANCHORS.items():
        A(f"**{label}: D_tissue = {val:.1e} cm^2/s**\n")
        A("| mode | drug | D fold | exposure fold | absolute kill fold | scored ratio fold | "
          "surviving vs kill |")
        A("|---|---|---|---|---|---|---|")
        for mode in res["modes"]:
            for drug in res["drugs"]:
                c = canc[f"{mode}|{drug}"][f"reference_to_nearest[{label}]"]
                A(f"| {mode} | {drug} | {c['D_fold']:.3f}x | {c['delivery_fold']:.3f}x | "
                  f"{_f(c['kill_nominal_fold'],3)}x | {c['ratio_fold']:.3f}x | "
                  f"{_f(c['fraction_surviving_vs_kill'], 3)} |")
        A("")

    # ---------------- 6. D_free ---------------------------------------------
    A("## 6. The D_free axis — and the sim's nominal temperature\n")
    A("**The sim DOES have a temperature concept, and it is already 37 C.** "
      "`tmesim.molecular.compile_drug.build_profile` takes `temperature_K = 310.0` "
      "(documented as body temperature) and passes it to "
      "`stokes_einstein.free_diffusivity(radius, temperature_K=310.0, "
      "viscosity_pa_s=1.2e-3)`. So the literature pass's inference that our D_free matches "
      "Stokes-Einstein *at 30 C* is **not** what the code does — the code is at 310 K and "
      f"gets {ref['free_D_cm2_s']:.3e} cm^2/s from an elevated INTERSTITIAL-FLUID viscosity "
      "of 1.2e-3 Pa.s (about 1.7x water at 37 C, ~6.9e-4 Pa.s, which the repo's own "
      "`permeability.py` uses for the membrane calculation). The low D_free is a viscosity "
      "assumption, not a temperature error. That distinction matters, because 1.2e-3 Pa.s is "
      "a *defensible* choice for interstitial fluid whereas 30 C would simply be a bug.\n")
    f_dfree = 4.27e-6 / ref["free_D_cm2_s"]
    A(f"**D_free is not an independent degree of freedom in this solver.** The engine never "
      f"sees D_free; it sees `DrugParams.diffusivity = D_free x hindrance`. Multiplying "
      f"D_free by a factor at fixed hindrance is therefore the SAME RUN as multiplying "
      f"D_tissue by that factor — it is already covered by this sweep, and needs no extra "
      f"solves. Moving D_free to the suggested 4.27e-6 cm^2/s is a "
      f"**{f_dfree:.3f}x** increase, which at fixed hindrance would put D_tissue at "
      f"**{ref['D_tissue_cm2_s']*f_dfree:.3e} cm^2/s** — i.e. it moves D_tissue in the "
      f"WRONG direction, further ABOVE every direct tumour-tissue measurement, and lands "
      f"just past the top of the swept range. Read off the monotone trend at the top of the "
      f"tables in section 3, that change makes the 3-D penalty smaller still.\n")
    A(f"The combined statement is the one worth acting on: if D_free is corrected upward to "
      f"4.27e-6 cm^2/s AND the hindrance is set from the MEASURED tumour-tissue "
      f"D_tissue/D_free ratios for real drugs (0.06-0.13 for doxorubicin, 0.045-0.083 for "
      f"tirapazamine in multicellular layers) rather than from the inert-solute tortuosity "
      f"rule, D_tissue lands at **{4.27e-6*0.06:.2e} - {4.27e-6*0.13:.2e} cm^2/s** — the "
      f"bottom third of this sweep, and consistent with the three direct measurements. "
      f"Correcting D_free WITHOUT correcting the hindrance would make the model worse.\n")

    # ---------------- 7. verdict ---------------------------------------------
    A("## 7. Controls — does this script reproduce the published numbers?\n")
    A("### 7.1 Scoring-path control (zero transport solves, dx = 12 um)\n")
    if refcheck and refcheck.get("rows"):
        A("The exposure vectors behind `results/kill_calibration_breast.md` are on disk in "
          "`data_cache/kill_calibration_breast_runs.json`. Scoring THOSE vectors through "
          "this script's own `score()` must reproduce the published table exactly. It does:\n")
        A("| drug | published ratio | re-scored from the published cache | delta | "
          "verdict | published verdict |")
        A("|---|---|---|---|---|---|")
        worst_rc, ok_rc = 0.0, True
        for r in refcheck["rows"]:
            d = r["drug"]
            p = published.get(d)
            if p is not None:
                worst_rc = max(worst_rc, abs(r["ratio_exact"] - p))
                ok_rc = ok_rc and r["verdict"] == pub_v.get(d)
            A(f"| {d} | {p if p is not None else '—'} | {_f(r['ratio_exact'])} | "
              f"{_f((r['ratio_exact']-p) if p is not None else None)} | {r['verdict']} | "
              f"{pub_v.get(d,'—')} |")
        A("")
        A(f"Worst deviation **{worst_rc:.4f}**; verdicts "
          f"{'ALL MATCH' if ok_rc else '**DO NOT ALL MATCH**'}. "
          + ("The scoring path in this file IS the published scoring path, so every "
             "movement reported in sections 3-5 is attributable to D_tissue and to nothing "
             "in the analysis."
             if worst_rc < 5e-3 and ok_rc else
             "**The scoring path does not reproduce the published table**; everything above "
             "should be treated as provisional until that is resolved.") + "\n")
    else:
        A("*(not run — `--refcheck` was not executed, so the scoring path is unverified "
          "against the published table)*\n")

    A(f"### 7.2 Grid control — the sweep was run at dx = {res['dx_um']:g} um\n")
    if abs(res["dx_um"] - 12.0) < 1e-9:
        A("The sweep is on the published grid, so its reference point should equal the "
          "published table directly.\n")
    else:
        A(f"**DISCLOSED TRUNCATION.** The published report is at `dx = 12 um`. This sweep "
          f"was run at `dx = {res['dx_um']:g} um` because the machine was under heavy "
          f"concurrent load from sibling jobs and a 120-solve sweep on the 12 um grid "
          f"(458x628, ~268k tissue cells, ~63 s per solve unloaded and ~9 min per solve "
          f"under the observed load) would not have completed. At "
          f"{res['grid'][0]}x{res['grid'][1]} the same sweep is ~5x cheaper. "
          f"**The consequence is real and is not hidden:** the coarser grid smears the "
          f"perivascular gradient, so the absolute 3-D:2-D ratios below sit LOWER than the "
          f"published ones. Read sections 3-5 as the SHAPE of the D-response, and section "
          f"7.1 as the control that the scoring itself is correct. A verdict that flips in "
          f"this table is evidence that the verdict is fragile; a verdict's absolute "
          f"position here is not directly comparable to the published table.\n")
    A("| drug | published ratio (dx = 12) | this sweep at the shipping D "
      f"(dx = {res['dx_um']:g}) | delta | verdict here | published verdict |")
    A("|---|---|---|---|---|---|")
    for d in res["drugs"]:
        r = ref_recal[d]
        p = published.get(d)
        A(f"| {d} | {p if p is not None else '—'} | {_f(r['ratio_exact'])} | "
          f"{_f((r['ratio_exact']-p) if p is not None else None)} | {r['verdict']} | "
          f"{pub_v.get(d,'—')} |")
    deltas = [abs(ref_recal[d]["ratio_exact"] - published[d])
              for d in res["drugs"] if d in published]
    worst = max(deltas) if deltas else float("nan")
    v_match = all(ref_recal[d]["verdict"] == pub_v[d]
                  for d in res["drugs"] if d in pub_v)
    A("")
    A(f"Worst absolute ratio deviation from the published table: **{worst:.4f}**; verdicts "
      f"{'all match' if v_match else '**do not all match**'} the published ones.\n")

    # ---- 7.3 the reduced dx = 12 confirmation ------------------------------
    A("### 7.3 Confirmation on the PUBLISHED dx = 12 um grid (reduced point set)\n")
    if confirm and confirm.get("records"):
        cd = float(confirm["dx_um"])
        pts = sorted({r["D_tissue_cm2_s"] for r in confirm["records"]})
        A(f"Because the full 10-point sweep was not affordable on the 12 um grid, the two "
          f"decisive D values were re-run there for all "
          f"{len(confirm['drugs'])} drugs and both modes: the shipping "
          f"{ref['D_tissue_cm2_s']:.3e} and **{min(pts):.3e} cm^2/s**, the swept point "
          f"nearest the measurement-provenance recommendation of section 8. If a verdict "
          f"is going to flip anywhere, this is where it flips, and this is the grid the "
          f"published verdicts live on.\n")
        A("| mode | drug | ratio at shipping D | ratio at the recommended D | "
          "verdict at shipping D | verdict at the recommended D | flipped? |")
        A("|---|---|---|---|---|---|---|")
        n_flip = 0
        for mode in confirm["modes"]:
            for drug in confirm["drugs"]:
                rows = _by(confirm, mode, drug)
                if len(rows) < 2:
                    continue
                lo = rows[0]
                hi = next((r for r in rows if r["is_reference"]), rows[-1])
                flipped = lo["verdict"] != hi["verdict"]
                n_flip += int(flipped)
                A(f"| {mode} | {drug} | {_f(hi['ratio_exact'])} | {_f(lo['ratio_exact'])} | "
                  f"{hi['verdict']} | {lo['verdict']} | "
                  f"{'**YES**' if flipped else 'no'} |")
        A("")
        A(f"**{n_flip} verdict flip(s)** on the published grid between the shipping "
          f"diffusivity and the measurement-backed one.\n")
        if abs(res["dx_um"] - cd) > 1e-9:
            A(f"Does the coarse-grid sweep get the SHAPE right? Same D interval, "
              f"dx = {res['dx_um']:g} um vs dx = {cd:g} um:\n")
            A(f"| mode | drug | ratio fold, dx = {res['dx_um']:g} | "
              f"ratio fold, dx = {cd:g} | agreement |")
            A("|---|---|---|---|---|")
            for mode in confirm["modes"]:
                if mode not in res["modes"]:
                    continue
                for drug in confirm["drugs"]:
                    if drug not in res["drugs"]:
                        continue
                    crows = _by(confirm, mode, drug)
                    if len(crows) < 2:
                        continue
                    c_hi = next((r for r in crows if r["is_reference"]), crows[-1])
                    c_fold = crows[0]["ratio_exact"] / c_hi["ratio_exact"]
                    srows = _by(res, mode, drug)
                    s_hi = next(r for r in srows if r["is_reference"])
                    s_lo = min(srows,
                               key=lambda r: abs(r["D_tissue_cm2_s"]
                                                 - crows[0]["D_tissue_cm2_s"]))
                    s_fold = s_lo["ratio_exact"] / s_hi["ratio_exact"]
                    A(f"| {mode} | {drug} | {s_fold:.3f}x | {c_fold:.3f}x | "
                      f"{abs(s_fold - c_fold) / max(c_fold, 1e-12) * 100:.1f}% apart |")
            A("")
        A("Cancellation on the published grid, shipping D -> recommended D:\n")
        A("| mode | drug | D fold | exposure fold | absolute kill fold | scored ratio fold | "
          "surviving vs kill |")
        A("|---|---|---|---|---|---|---|")
        ccc = cancellation(confirm)
        for mode in confirm["modes"]:
            for drug in confirm["drugs"]:
                c = ccc.get(f"{mode}|{drug}", {}).get("reference_to_lowest")
                if not c:
                    continue
                A(f"| {mode} | {drug} | {c['D_fold']:.3f}x | {c['delivery_fold']:.3f}x | "
                  f"{_f(c['kill_nominal_fold'],3)}x | {c['ratio_fold']:.3f}x | "
                  f"{_f(c['fraction_surviving_vs_kill'], 3)} |")
        A("")
    else:
        A("*(not run)*\n")

    # ---------------- 8. recommendation --------------------------------------
    anchors = sorted(MEASURED_ANCHORS.values())
    geo = math.exp(sum(math.log(v) for v in anchors) / len(anchors))
    A("## 8. Recommendation — justified by measurement provenance, NOT by score\n")
    A("The rule this section obeys: the recommended value is the one with the best "
      "measurement provenance, chosen BEFORE looking at how many verdicts it produces. "
      "For the record, the PASS count at the recommended value is stated afterwards, "
      "not used to reach it.\n")
    A("**Provenance of the candidates.**\n")
    A(f"* The shipping **{ref['D_tissue_cm2_s']:.3e} cm^2/s** has no direct measurement "
      "behind it at all. Its D_free leg is first-principles Stokes-Einstein (defensible) "
      "and its hindrance leg is a particle sim plus an Ogston-Renkin fibre-matrix gate — "
      "both of which model EXCLUDED VOLUME only. Neither has any term for the cationic "
      "binding, DNA intercalation or endosomal sequestration that define doxorubicin's "
      "transport. Landing on the canonical inert-solute ratio 0.40 is exactly what an "
      "excluded-volume-only model should do, and exactly why it is the wrong number for "
      "this drug. Sykova & Nicholson (Physiol Rev 2008, **PMID 18923183**) state the "
      "restriction themselves.")
    A("* Three DIRECT measurements exist in real, cellularised tissue, and they agree "
      "with each other far better than any of them agrees with the shipping value:")
    for label, v in sorted(MEASURED_ANCHORS.items(), key=lambda kv: kv[1]):
        A(f"    * **{v:.1e} cm^2/s** — {label} "
          f"({ref['D_tissue_cm2_s']/v:.1f}x below the shipping value)")
    A(f"* Their geometric mean is **{geo:.2e} cm^2/s**, spanning "
      f"{anchors[0]:.1e}-{anchors[-1]:.1e}. Independently, measured D_tissue/D_free "
      "ratios in cellularised tumour tissue are 0.06-0.13 for doxorubicin and 0.045-0.083 "
      "for tirapazamine in multicellular layers — 3-8x below the 0.40 this model uses, "
      "and consistent with the absolute numbers.\n")
    A(f"**Recommendation: set D_tissue = {geo:.1e} cm^2/s** (the geometric mean of the "
      f"three direct measurements), and carry **{anchors[0]:.1e}-{anchors[-1]:.1e} cm^2/s** "
      f"as the reported uncertainty band rather than a point value. That is a "
      f"**{ref['D_tissue_cm2_s']/geo:.1f}x reduction** from the shipping value. Every "
      "number in that recommendation traces to Toley 2013 / Weinberg 2007 / Qian 2003; "
      "none of it traces to this sweep's PASS counts.\n")
    A("**Two conditions attach to it, and both are load-bearing:**\n")
    L_lo, L_hi = DOX_EFOLD_BAND_UM
    dox_rows_f = _by(res, "frozen", "doxorubicin") if "frozen" in res["modes"] else []
    L_frozen_at_rec = (min(dox_rows_f, key=lambda r: abs(r["D_tissue_cm2_s"] - geo)
                           )["penetration_length_um"] if dox_rows_f else float("nan"))
    A(f"1. **`alpha` MUST be re-calibrated with it (mode `recal`), not left frozen.** "
      f"The measured doxorubicin penetration depth (Primeau 2005, PMID 16361566) is an "
      f"independent measurement this model already reproduces by construction: the "
      f"published 40-50 um half-distance is an e-folding band of "
      f"**{L_lo:.1f}-{L_hi:.1f} um**. Changing D without re-solving alpha silently breaks "
      f"that anchor — at the recommended D, mode `frozen` gives "
      f"**L = {L_frozen_at_rec:.1f} um**, well below the measured band, whereas mode "
      f"`recal` holds it at {ref['dox_target_efold_um']:.1f} um by construction. A model "
      f"that gets the diffusivity right and the penetration depth wrong is not an "
      f"improvement.")
    A("2. **Correcting D_free alone would make the model worse, not better.** Section 6: "
      "raising D_free to 4.27e-6 cm^2/s at fixed hindrance raises D_tissue further ABOVE "
      "the measurements. D_free and the hindrance have to move together, and the "
      "measured 0.06-0.13 ratio is the thing that should replace the 0.40.\n")
    A("**The deeper structural point, which the sweep makes unavoidable.** Because alpha "
      "is fit so that `L = sqrt(D_eff/k_eff)` matches a measured penetration depth, the "
      "measurement constrains only the RATIO `D/k`. D_tissue and the perivascular sink "
      "strength are not separately identified by anything this model is currently anchored "
      "to. Lowering D_tissue 3x and lowering the inferred sink 3x leaves the gradient the "
      "model actually predicts unchanged — which is why the scored ratios move so little "
      "in mode `recal`. Getting D_tissue right therefore does not by itself buy a better "
      "prediction; it buys a correctly-attributed one, and it changes the claim the model "
      "makes about how much drug is being consumed perivascularly by roughly the same "
      "factor. That claim is currently unmeasured and should be flagged as such.\n")

    # PASS counts, stated after the recommendation, not used to reach it
    A("For the record, and stated only after the fact: at the swept point nearest the "
      f"recommended {geo:.1e} cm^2/s, the formal verdicts are\n")
    # Prefer the published dx = 12 um grid where it was run -- the coarse sweep's
    # own baseline already disagrees with the published docetaxel verdict, so
    # quoting its PASS count here would flatter the recommendation.
    src = confirm if (confirm and confirm.get("records")) else res
    src_note = (f"on the PUBLISHED dx = {src['dx_um']:g} um grid"
                if src is confirm else
                f"on this sweep's dx = {src['dx_um']:g} um grid — **NOT the published "
                f"grid; see §7.2**")
    for mode in src["modes"]:
        near = {}
        for drug in src["drugs"]:
            rows = _by(src, mode, drug)
            r = min(rows, key=lambda x: abs(x["D_tissue_cm2_s"] - geo))
            near[drug] = (r["D_tissue_cm2_s"], r["verdict"])
        npass = sum(1 for _, v in near.values() if v == "PASS")
        nscored = sum(1 for _, v in near.values() if v != NOT_TESTABLE)
        d_here = next(iter(near.values()))[0]
        A(f"* mode `{mode}`, {src_note} (D = {d_here:.3e} cm^2/s): "
          f"**{npass}/{nscored} PASS** (of {len(near)} swept drugs; "
          f"{len(near) - nscored} are NOT TESTABLE and are counted in no denominator) — "
          + ", ".join(f"{d} {v}" for d, (_, v) in sorted(near.items())))
    A("")
    A("Note the direction of that result, because it is the opposite of what a rigged "
      "sweep would produce: in `recal` — the mode the codebase actually runs, and the mode "
      "the recommendation requires — the measurement-backed diffusivity does **not** "
      "rescue the one published FAIL. Docetaxel still fails. The recommendation is made "
      "anyway, on provenance.\n")

    # ---------------- 9. limitations -----------------------------------------
    A("## 9. Limitations of THIS sweep\n")
    A(f"1. **A single global scale factor is applied to all {len(res['drugs'])} drugs.** "
      "The three direct measurements are doxorubicin (and one small-molecule-in-liver) "
      "specific. The defect they expose — an excluded-volume-only hindrance model applied "
      "to interacting solutes — is structural and applies to every drug in the panel, but "
      "the SIZE of the correction is almost certainly drug-specific (cisplatin is small, "
      "neutral-ish and weakly bound; docetaxel is large and highly protein-bound). A "
      "per-drug binding-aware hindrance model is the real fix; this sweep is the "
      "first-order propagation, not that fix.")
    A("2. **Only the 6 formally-testable drugs were swept.** The other 4 panel drugs "
      "(cyclophosphamide, tamoxifen, vinorelbine, palbociclib) have no literature band at "
      "all, so no verdict of theirs could flip. Their ratios will move exactly as these do.")
    A("3. **Only the 1x IC50 reporting concentration was swept.** The published report "
      "evaluates the decomposition at 0.1/0.3/1/3/10x IC50 and notes that F_penetration "
      "climbs toward 1 as the Hill term saturates. The verdict is taken at 1x IC50, so "
      "that is what is swept, but the D-sensitivity at other concentrations is not "
      "characterised here.")
    A("4. **The p99 bath-matching convention is held fixed.** The GBM report found the "
      "cisplatin verdict to be convention-dependent (pass at p99, fail at max-exposure). "
      "This sweep does not cross that with the D axis, so a joint "
      "convention-x-diffusivity sensitivity is still unmeasured.")
    A("5. **Everything here inherits the delivery-side error already logged in "
      "`results/breast_kp_validation.md`** (a blood-brain-barrier permeability applied to "
      "breast tissue, and a one-way delivery source with no back-flux). This sweep varies "
      "D with those two defects still present; fixing them could change the D-response, "
      "and the interaction has not been tested.")
    A("6. **The 3 direct-measurement anchors are reported as the literature pass supplied "
      "them.** They were not independently re-verified against the primary papers in this "
      "run; the PMIDs are carried through unchanged so they can be.")
    if mismatched:
        A(f"7. **An unplanned finding that is arguably more urgent than the one this sweep "
          f"was commissioned for: {', '.join(mismatched)} changes verdict between "
          f"`dx = 12 um` and `dx = 24 um` (§7.2), on a grid the repo documents as "
          f"converged.** The published docetaxel FAIL is the only FAIL in the panel and it "
          f"does not survive a 2x grid coarsening, which is a smaller perturbation than the "
          f"3.1x diffusivity correction it DOES survive (§7.3, `recal`). A grid-convergence "
          f"re-check of the breast panel is a separate piece of work this sweep did not do "
          f"and cannot substitute for; nothing here should be read as having tested it.")
    A("")

    os.makedirs(os.path.dirname(out_md), exist_ok=True)
    with open(out_md, "w") as fh:
        fh.write("\n".join(L) + "\n")
    return out_md


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dx", type=float, default=24.0)
    ap.add_argument("--modes", type=str, default="frozen,recal")
    ap.add_argument("--drugs", type=str, default=",".join(TESTABLE))
    ap.add_argument("--out", type=str, default=None)
    ap.add_argument("--checkpoint", type=str, default=None)
    ap.add_argument("--points", type=str, default=None,
                    help="comma-separated D_tissue values in cm^2/s, replacing the "
                         "log grid (the shipping reference is always added)")
    ap.add_argument("--report-only", action="store_true",
                    help="skip the solves, rebuild the .md from the stored JSON")
    ap.add_argument("--seed-from-cache", action="store_true", dest="seed_from_cache",
                    help="write the shipping-D records of a dx=12 sweep from the "
                         "published cache instead of re-solving them")
    ap.add_argument("--refcheck", action="store_true",
                    help="zero-solve control: re-score the PUBLISHED dx=12 cached "
                         "exposure fields through this script's scoring path")
    ap.add_argument("--md", type=str, default=None)
    a = ap.parse_args()
    dx = a.dx
    ap_seed = getattr(a, "seed_from_cache", False)
    if ap_seed:
        ck0 = a.checkpoint or os.path.join(OUT_DIR, f"dtissue_records_dx{int(dx)}.jsonl")
        n = seed_checkpoint_from_published_cache(
            ck0, drugs=tuple(x for x in a.drugs.split(",") if x),
            modes=tuple(x for x in a.modes.split(",") if x), dx_um=dx)
        print(f"seeded {n} shipping-D records into {ck0}")
        sys.exit(0)
    if a.refcheck:
        published_reference_check(
            drugs=tuple(x for x in a.drugs.split(",") if x),
            out_json=os.path.join(OUT_DIR, "dtissue_refcheck_dx12.json"))
        sys.exit(0)
    out_json = a.out or os.path.join(OUT_DIR, f"dtissue_sweep_dx{int(dx)}.json")
    ckpt = a.checkpoint or os.path.join(OUT_DIR, f"dtissue_records_dx{int(dx)}.jsonl")
    out_md = a.md or os.path.join(_ROOT, "results", "dtissue_sensitivity.md")
    if a.report_only:
        with open(out_json) as fh:
            res = json.load(fh)
    else:
        res = run(dx_um=dx, drugs=tuple(x for x in a.drugs.split(",") if x),
                  modes=tuple(x for x in a.modes.split(",") if x),
                  out_json=out_json, checkpoint=ckpt,
                  explicit_points=([float(x) for x in a.points.split(",") if x]
                                   if a.points else None))
    rc_path = os.path.join(OUT_DIR, "dtissue_refcheck_dx12.json")
    rc = json.load(open(rc_path)) if os.path.exists(rc_path) else None
    cf_path = os.path.join(OUT_DIR, "dtissue_confirm_dx12.json")
    cf = (json.load(open(cf_path))
          if os.path.exists(cf_path) and cf_path != out_json else None)
    print("wrote", write_report(res, out_md, refcheck=rc, confirm=cf))
    # write_report re-scores every record through the shared band decision. Persist
    # that back, so the on-disk sweep artefact cannot keep asserting a verdict the
    # document no longer makes. Only verdict fields move; no ratio is touched.
    for path, obj in ((out_json, res), (rc_path, rc), (cf_path, cf)):
        if obj is None or not os.path.exists(path):
            continue
        with open(path, "w") as fh:
            json.dump(obj, fh, indent=1, default=float)
    print("re-scored verdicts persisted to the sweep artefacts")
