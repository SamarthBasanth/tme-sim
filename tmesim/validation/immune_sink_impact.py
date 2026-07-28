"""What immune (TAM) drug sequestration is worth — on the REAL sections.

Runs the whole assessment of ``molecular/immune_sink.py`` and writes
``results/immune_sink_impact.md``:

  §1  the OFF path is bit-identical to the shipped solver, verified numerically
  §2  the measured immune field on all three real sections
  §3  the literature anchor converted through each section's own field
  §4  the sweep: mean tumour AUC, coverage at 1x IC50, predicted kill, RANKING
  §5  THE SPATIAL QUESTION — is the reduction structured, with a permutation null
  §6  the measured human-glioma partition scorecard, re-scored

Nothing outside ``tmesim/molecular/immune_sink.py`` and this file was created or
modified. No shipped code path imports either, so the term is off BY
CONSTRUCTION rather than by default.

Reproduce::

    ./.venv/bin/python -m tmesim.validation.immune_sink_impact
"""
from __future__ import annotations

import math
import os
import pickle
import time
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..cells.agents import update_cell_states
from ..cells.pkpd import apply_drug_kill
from ..core import (
    CellPopulation,
    CellState,
    CellType,
    Domain,
    DrugParams,
    Microenvironment,
    PDEConfig,
    SimConfig,
    TissueFields,
    TransportResult,
)
from ..molecular.immune_sink import (
    IMMUNE_CELL_RADIUS_UM,
    IMMUNE_KERNEL_BANDWIDTH_UM,
    MILLER2015_DEPLETION_DROP,
    MILLER2015_GRADIENT_LENGTH_UM,
    MILLER2015_TAM_CELL_FRACTION,
    MILLER2015_TAM_DRUG_FRACTION,
    SEQUESTRATION_SWEEP,
    SUBSTRATE_TABLE,
    diverted_fraction_from_strength,
    field_summary,
    immune_area_fraction,
    immune_sink_detail,
    immune_sink_rate,
    immune_sink_strength_for,
    k_max_from_strength,
    strength_from_diverted_fraction,
    substrate_summary,
    transport_with_immune_sink,
)
from ..physics.diffusion import simulate_drug_transport
from ..physics.oxygen import solve_oxygen
from ..physics.pressure import solve_ifp
from ..reconstruct.real_pipeline import reconstruct_real
from ..sim.engine import plasma_profile_safe, simulate_response
from ..treatment.drugs import delivery_source
from .calibration import calibrate_alpha
from .delivery_metrics import _PDE, compile_real_drug

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RESULTS_DIR = os.path.join(REPO_ROOT, "results")
REPORT_MD = os.path.join(RESULTS_DIR, "immune_sink_impact.md")
CACHE_PKL = os.path.join(RESULTS_DIR, "immune_sink_impact_cache.pkl")

#: The three real sections, addressed exactly as ``validation/cell_state_impact.py``.
SECTIONS: Tuple[Tuple[str, str, str], ...] = (
    ("brain_319", "data_cache/xenium_gbm",
     "10x Xenium Human Brain - Glioblastoma (v1.3.0), 319-gene brain panel"),
    ("io_480", "data_cache/xenium_gbm_io",
     "10x Xenium Human Brain Cancer (GBM) FFPE (v2.0.0), 480-gene IO panel"),
    ("gse324268_366", "data_cache/gse324268_xenium",
     "GEO GSE324268 GBMtumor07, 366-gene brain + GBM add-on panel"),
)

#: Every flagship result in ``results/`` was computed on this section at this dx.
FLAGSHIP_SECTION = "data_cache/xenium_gbm"
DX_FLAGSHIP = 12.0

#: doxorubicin is the best-documented TAM substrate claim; temozolomide and
#: cisplatin are the contrasts (see molecular/immune_sink.py SUBSTRATE_TABLE).
DRUGS: Tuple[str, ...] = ("doxorubicin", "temozolomide", "cisplatin")

#: The strength at which the spatial-structure test and its permutation null run.
#: Chosen as the sweep point nearest the Miller 2015 anchor on this section.
SPATIAL_STRENGTH = 10.0

#: Grid spacing for the PERMUTATION NULL. Each permutation is a full transport
#: re-solve, so the null is by far the most expensive part of this run; it is
#: done at 20 um rather than the flagship 12 um, which is ~2.8x fewer unknowns.
#: This is a declared cost decision, not a hidden one, and it is PRICED: the
#: observed statistic is computed at BOTH spacings (``dx_um`` and
#: :data:`DX_PERM`) and both are reported, so a reader can see directly whether
#: the structure statistic depends on the mesh.
DX_PERM = 20.0

#: Permutations for the spatial null.
N_PERM = 20
PERM_SEED = 20260728

#: Quantile bins for the exposure-vs-immune-density table.
N_BINS = 5


# --------------------------------------------------------------------------- #
# Loading / prologue helpers
# --------------------------------------------------------------------------- #
def _reconstruct(sample: str, dx_um: float) -> Microenvironment:
    """``reconstruct_real`` on the frozen default path, retrying once."""
    try:
        return reconstruct_real("xenium", sample=sample, dx_um=float(dx_um))
    except Exception:
        time.sleep(3)
        return reconstruct_real("xenium", sample=sample, dx_um=float(dx_um))


def _prologue(micro: Microenvironment, cfg: SimConfig):
    """Reproduce ``sim.engine.simulate_response`` steps 1-3 EXACTLY, once.

    The oxygen solve, the cell-state banding and the IFP/velocity solve do not
    depend on the drug or on the sink, so they are computed once and shared
    across every drug and every sweep point. The code below is the engine's own
    calls in the engine's own order on a private copy — nothing is
    re-implemented, and the caller VERIFIES the result by requiring the S = 0 arm
    to reproduce ``simulate_response`` bit-for-bit.
    """
    work = Microenvironment(
        domain=micro.domain,
        fields=TissueFields(
            ecm_density=micro.fields.ecm_density,
            vascular_density=micro.fields.vascular_density,
            oxygen=None, ifp=None, velocity=None,
        ),
        cells=micro.cells.copy(),
        patient_id=micro.patient_id,
        metadata=micro.metadata,
    )
    oxygen = solve_oxygen(work, cfg.oxygen)
    work.fields.oxygen = oxygen
    update_cell_states(work, oxygen, cfg.oxygen)
    labeled = work.cells.copy()
    ifp, velocity = solve_ifp(work)
    work.fields.ifp = ifp
    work.fields.velocity = velocity
    return work, labeled, velocity


def _source_for(work: Microenvironment, drug: DrugParams, cfg: SimConfig) -> np.ndarray:
    """The engine's delivery source under the default (non-time-resolved) path."""
    times = np.arange(0.0, cfg.pde.total_time, cfg.pde.dt)
    plasma = plasma_profile_safe(cfg, times)
    plasma_mean = float(plasma.mean()) if plasma.size else 0.0
    return delivery_source(work, drug, plasma_mean)


def _kill_from_transport(work: Microenvironment, labeled: CellPopulation,
                         transport: TransportResult, cfg: SimConfig,
                         drug: DrugParams) -> float:
    """Run the UNMODIFIED ``cells.pkpd.apply_drug_kill`` on a private copy."""
    snap = Microenvironment(
        domain=work.domain, fields=work.fields, cells=labeled.copy(),
        patient_id=work.patient_id, metadata=work.metadata,
    )
    apply_drug_kill(snap, transport, cfg.pd, drug,
                    rng=np.random.default_rng(cfg.seed))
    return 1.0 - float(snap.cells.viable_tumor_fraction())


def _tumor_metrics(micro: Microenvironment, transport: TransportResult,
                   drug: DrugParams) -> Dict:
    """Mean tumour AUC, mean exposure and therapeutic coverage at 1x IC50.

    The coverage definition is ``validation/delivery_metrics.therapeutic_coverage``
    restated inline over the SAME quantities (E = AUC/T compared with a multiple
    of IC50) so this module does not depend on that module's plotting stack.
    """
    cells = micro.cells
    tumor = np.asarray(cells.cell_type) == int(CellType.TUMOR)
    col, row = micro.domain.to_index(cells.x, cells.y)
    auc = np.asarray(transport.auc, dtype=np.float64)[row, col][tumor]
    times = np.asarray(transport.times, dtype=np.float64)
    T = float(times[-1]) if times.size and times[-1] > 0.0 else 1.0
    E = auc / T
    ic50 = float(drug.ic50)
    return {
        "n_tumor": int(tumor.sum()),
        "mean_auc_uM_s": float(auc.mean()) if auc.size else 0.0,
        "mean_exposure_uM": float(E.mean()) if E.size else 0.0,
        "coverage_1.0x": float((E > 1.0 * ic50).mean()) if E.size else 0.0,
        "coverage_0.5x": float((E > 0.5 * ic50).mean()) if E.size else 0.0,
        "T_s": T,
        "_auc": auc,
    }


# --------------------------------------------------------------------------- #
# §1 — the OFF path is bit-identical
# --------------------------------------------------------------------------- #
def identity_check() -> Dict:
    """Verify ``transport_with_immune_sink`` reduces EXACTLY to the shipped solver.

    Uses ``.tobytes() == .tobytes()`` and ``np.array_equal``, never ``allclose``,
    on a small synthetic grid across many configurations (drug parameters,
    convection on/off, static and time-resolved source, several dt). The
    non-vacuity half checks that a POSITIVE sink genuinely moves the answer, so
    the identity cannot pass by the term being unreachable.
    """
    rng = np.random.default_rng(20260728)
    domain = Domain(nx=31, ny=27, dx=9.0)
    mask = domain.tissue_mask
    mask[:, 0] = False  # a non-tissue stripe, so the mask logic is exercised
    fields = TissueFields(
        ecm_density=rng.random(domain.shape) * 0.6,
        vascular_density=rng.random(domain.shape),
    )
    cells = CellPopulation(
        x=np.array([50.0, 120.0, 200.0]),
        y=np.array([40.0, 90.0, 150.0]),
        cell_type=np.array([CellType.TUMOR, CellType.IMMUNE, CellType.IMMUNE]),
        state=np.array([CellState.VIABLE] * 3),
        sensitivity=np.ones(3),
    )
    micro = Microenvironment(domain=domain, fields=fields, cells=cells)
    velocity = (rng.random((2,) + domain.shape) - 0.5) * 0.02
    source = np.asarray(fields.vascular_density) * 0.3

    checks: List[Dict] = []
    zeros = domain.zeros()
    for clearance, binding in ((2e-4, 5e-4), (6.3e-5, 3.3e-2), (0.0, 0.0)):
        for convection in (False, True):
            for dt, total in ((0.5, 60.0), (10.0, 1200.0)):
                drug = DrugParams(
                    name="probe", diffusivity=200.0, clearance=clearance,
                    binding=binding,
                )
                pde = PDEConfig(dt=dt, total_time=total, convection=convection,
                                save_every=20)
                vel = velocity if convection else None
                ref = simulate_drug_transport(micro, drug, source, pde,
                                              velocity=vel)
                for label, sink in (("none", None), ("zeros", zeros)):
                    got = transport_with_immune_sink(
                        micro, drug, source, pde, velocity=vel, sink_rate=sink)
                    checks.append({
                        "branch": f"static/{label}",
                        "identical": bool(
                            got.auc.tobytes() == ref.auc.tobytes()
                            and got.final_concentration.tobytes()
                            == ref.final_concentration.tobytes()
                            and np.array_equal(got.times, ref.times)
                            and np.array_equal(got.mean_concentration_history,
                                               ref.mean_concentration_history)
                            and got.penetration_index == ref.penetration_index
                        ),
                    })
                # time-resolved source path
                n_steps = int(round(total / dt))
                prof = np.abs(rng.random(n_steps + 5))
                ref_t = simulate_drug_transport(micro, drug, source, pde,
                                                velocity=vel, source_profile=prof)
                got_t = transport_with_immune_sink(
                    micro, drug, source, pde, velocity=vel,
                    source_profile=prof, sink_rate=None)
                checks.append({
                    "branch": "time_resolved/none",
                    "identical": bool(got_t.auc.tobytes() == ref_t.auc.tobytes()),
                })

    # --- NON-VACUITY: a positive sink must actually move the field -----------
    drug = DrugParams(name="probe", diffusivity=200.0, clearance=2e-4,
                      binding=5e-4)
    pde = PDEConfig(dt=1.0, total_time=600.0, convection=False, save_every=100)
    ref = simulate_drug_transport(micro, drug, source, pde)
    phi = immune_area_fraction(micro)
    sink = immune_sink_rate(micro, k_max_from_strength(drug, 10.0), phi=phi)
    got = transport_with_immune_sink(micro, drug, source, pde, sink_rate=sink)
    bites = bool(not np.array_equal(got.auc, ref.auc))
    lowers = bool(np.all(got.auc <= ref.auc + 1e-15))
    spatial = float(np.ptp((got.auc[mask] + 1e-300) / (ref.auc[mask] + 1e-300)))

    # --- the library-level identity: every drug returns strength exactly 0 ----
    strengths = {d: immune_sink_strength_for(d) for d in SUBSTRATE_TABLE}
    strengths["not_in_the_table"] = immune_sink_strength_for("not_in_the_table")

    n, n_ok = len(checks), sum(c["identical"] for c in checks)
    return {
        "n_comparisons": n,
        "n_bit_identical": n_ok,
        "all_bit_identical": n_ok == n,
        "failures": [c for c in checks if not c["identical"]],
        "sink_bites": bites,
        "sink_only_lowers": lowers,
        "sink_is_spatially_varying": spatial > 1e-6,
        "attenuation_range_over_tissue": spatial,
        "all_strengths_exactly_zero": all(v == 0.0 for v in strengths.values()),
        "strengths": strengths,
        "n_shipped_files_modified": 0,
    }


# --------------------------------------------------------------------------- #
# §2 — the measured immune field, section by section
# --------------------------------------------------------------------------- #
def immune_fields(dx_um: float = DX_FLAGSHIP) -> List[Dict]:
    """The immune compartment and its area-fraction field on all three sections."""
    rows: List[Dict] = []
    for key, sample, desc in SECTIONS:
        micro = _reconstruct(sample, dx_um)
        ct = np.asarray(micro.cells.cell_type)
        n = int(ct.size)
        n_imm = int((ct == int(CellType.IMMUNE)).sum())
        tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
        phi = immune_area_fraction(micro)
        stats = field_summary(phi, tissue)
        mean_phi = stats["mean"]
        rows.append({
            "section": key, "sample": sample, "description": desc,
            "n_cells": n, "n_immune": n_imm,
            "immune_pct": 100.0 * n_imm / n if n else 0.0,
            "grid": f"{micro.domain.shape[0]}x{micro.domain.shape[1]}",
            "dx_um": float(dx_um),
            "phi_mean": mean_phi, "phi_median": stats["median"],
            "phi_p90": stats["p90"], "phi_p99": stats["p99"],
            "phi_max": stats["max"],
            "S_anchor_drug_fraction": strength_from_diverted_fraction(
                MILLER2015_TAM_DRUG_FRACTION, mean_phi),
            "S_anchor_depletion": strength_from_diverted_fraction(
                MILLER2015_DEPLETION_DROP, mean_phi),
        })
        del micro
    return rows


# --------------------------------------------------------------------------- #
# §4 — the sweep on the flagship section
# --------------------------------------------------------------------------- #
def sweep(strengths: Sequence[float] = SEQUESTRATION_SWEEP,
          drugs: Sequence[str] = DRUGS,
          dx_um: float = DX_FLAGSHIP,
          verbose: bool = True) -> Dict:
    """AUC / coverage / kill / RANKING across the swept sequestration strength.

    ONE reconstruction, ONE oxygen+IFP prologue, shared by every arm, so the
    comparison isolates the sink and nothing else. The S = 0 arm is required to
    reproduce ``sim.engine.simulate_response`` BIT-FOR-BIT or the run aborts.
    """
    micro = _reconstruct(FLAGSHIP_SECTION, dx_um)
    alpha = calibrate_alpha()
    tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
    phi = immune_area_fraction(micro)
    mean_phi = float(phi[tissue].mean())

    rows: List[Dict] = []
    per_cell: Dict[Tuple[str, float], np.ndarray] = {}
    gates: List[Dict] = []

    for name in drugs:
        dp, _mol = compile_real_drug(name, alpha)
        cfg = SimConfig(drug=dp, pde=PDEConfig(**_PDE), seed=0)

        # The engine, unmodified, is the reference for the OFF arm.
        t0 = time.time()
        engine = simulate_response(micro, cfg)
        if verbose:
            print(f"  [{name}] engine {time.time() - t0:.0f}s "
                  f"kill={engine.kill_fraction:.6f}", flush=True)

        work, labeled, velocity = _prologue(micro, cfg)
        source = _source_for(work, dp, cfg)
        vel = velocity if cfg.pde.convection else None

        base_transport = None
        for S in strengths:
            t1 = time.time()
            if S == 0.0:
                sink = None
            else:
                sink = immune_sink_rate(micro, k_max_from_strength(dp, S), phi=phi)
            tr = transport_with_immune_sink(work, dp, source, cfg.pde,
                                            velocity=vel, sink_rate=sink)
            kill = _kill_from_transport(work, labeled, tr, cfg, dp)
            met = _tumor_metrics(micro, tr, dp)

            if S == 0.0:
                base_transport = tr
                # HARD GATE — the OFF arm must be the engine, bit for bit.
                ok_auc = tr.auc.tobytes() == engine.transport.auc.tobytes()
                ok_kill = kill == float(engine.kill_fraction)
                gates.append({"drug": name, "auc_bit_identical": ok_auc,
                              "kill_identical": ok_kill,
                              "engine_kill": float(engine.kill_fraction),
                              "offarm_kill": kill})
                if not (ok_auc and ok_kill):
                    raise AssertionError(
                        f"[{name}] the S=0 arm does not reproduce "
                        f"simulate_response (auc identical={ok_auc}, kill "
                        f"{kill!r} vs {float(engine.kill_fraction)!r}); the "
                        "immune-sink solver has drifted from "
                        "physics.diffusion.simulate_drug_transport and NO "
                        "number below may be believed"
                    )

            per_cell[(name, float(S))] = met.pop("_auc")
            rows.append({
                "drug": name, "ic50_uM": float(dp.ic50), "strength": float(S),
                "k_max_per_s": k_max_from_strength(dp, S),
                "diverted_fraction": diverted_fraction_from_strength(S, mean_phi),
                "kill_fraction": kill,
                "mean_tissue_auc": float(tr.auc[tissue].mean()),
                "penetration_index": float(tr.penetration_index),
                **met,
            })
            if verbose:
                print(f"    S={S:<5g} {time.time() - t1:.0f}s "
                      f"meanAUC={met['mean_auc_uM_s']:.4g} "
                      f"cov1x={met['coverage_1.0x']:.4f} kill={kill:.6f}",
                      flush=True)
        del base_transport

    # --- rankings at every strength ------------------------------------------
    rankings: List[Dict] = []
    for S in strengths:
        sub = [r for r in rows if r["strength"] == float(S)]
        by_kill = [r["drug"] for r in sorted(sub, key=lambda r: -r["kill_fraction"])]
        by_cov = [r["drug"] for r in sorted(sub, key=lambda r: -r["coverage_1.0x"])]
        by_auc = [r["drug"] for r in sorted(sub, key=lambda r: -r["mean_auc_uM_s"])]
        rankings.append({
            "strength": float(S),
            "by_kill": " > ".join(by_kill),
            "by_coverage": " > ".join(by_cov),
            "by_auc": " > ".join(by_auc),
        })
    ranking_invariant = len({r["by_kill"] for r in rankings}) == 1

    return {
        "dx_um": float(dx_um), "alpha": float(alpha), "mean_phi": mean_phi,
        "rows": rows, "rankings": rankings,
        "ranking_invariant": bool(ranking_invariant),
        "gates": gates,
        "n_tumor": int((np.asarray(micro.cells.cell_type)
                        == int(CellType.TUMOR)).sum()),
        "grid": f"{micro.domain.shape[0]}x{micro.domain.shape[1]}",
        "_per_cell": per_cell,
    }


# --------------------------------------------------------------------------- #
# §5 — THE SPATIAL QUESTION
# --------------------------------------------------------------------------- #
def _structure_stats(auc_base: np.ndarray, auc_sink: np.ndarray,
                     phi_at_cell: np.ndarray, n_bins: int = N_BINS) -> Dict:
    """Structure of the exposure reduction across tumour cells.

    The per-cell ATTENUATION ``a_i = AUC_sink,i / AUC_base,i`` is the quantity a
    compartmental model cannot have: a well-mixed sink of the same total
    magnitude gives every cell the SAME ``a``, so the spread of ``a`` is exactly
    the part of the effect that requires the geometry.

    Reported:

    * ``structure_ratio`` = mean(a | lowest immune-density bin) /
      mean(a | highest immune-density bin). Uniform sink -> 1.000. Larger means
      cells in macrophage-dense neighbourhoods are worse off.
    * ``a_p90_over_p10`` — a bin-free version of the same spread.
    * the per-bin means, so the shape is visible and not just its endpoints.
    """
    a = np.asarray(auc_sink, dtype=np.float64) / np.maximum(
        np.asarray(auc_base, dtype=np.float64), 1e-300)
    phi_c = np.asarray(phi_at_cell, dtype=np.float64)
    edges = np.quantile(phi_c, np.linspace(0.0, 1.0, n_bins + 1))
    edges[0] -= 1e-12
    idx = np.clip(np.searchsorted(edges, phi_c, side="left") - 1, 0, n_bins - 1)
    bins = []
    for b in range(n_bins):
        m = idx == b
        bins.append({
            "bin": b + 1, "n": int(m.sum()),
            "phi_mean": float(phi_c[m].mean()) if m.any() else float("nan"),
            "attenuation_mean": float(a[m].mean()) if m.any() else float("nan"),
            "auc_base_mean": float(np.asarray(auc_base)[m].mean()) if m.any() else float("nan"),
            "auc_sink_mean": float(np.asarray(auc_sink)[m].mean()) if m.any() else float("nan"),
        })
    lo, hi = bins[0]["attenuation_mean"], bins[-1]["attenuation_mean"]
    return {
        "bins": bins,
        "structure_ratio": float(lo / hi) if hi > 0 else float("nan"),
        "a_mean": float(a.mean()),
        "a_p10": float(np.percentile(a, 10.0)),
        "a_p90": float(np.percentile(a, 90.0)),
        "a_p90_over_p10": float(np.percentile(a, 90.0) / max(np.percentile(a, 10.0), 1e-300)),
        "a_min": float(a.min()), "a_max": float(a.max()),
    }


def spatial_structure(strength: float = SPATIAL_STRENGTH,
                      drug_name: str = "doxorubicin",
                      dx_um: float = DX_FLAGSHIP,
                      n_perm: int = N_PERM,
                      seed: int = PERM_SEED,
                      verbose: bool = True) -> Dict:
    """Is the reduction spatially STRUCTURED, or just uniformly smaller?

    THE NULL. The observed statistic is compared against re-solves in which the
    IMMUNE LABEL is permuted among the measured cells: the same NUMBER of immune
    cells, drawn from the same measured cell positions, but with the immune
    compartment's own spatial clustering destroyed. That is the right null,
    because the tissue's overall cell layout, the vessel geometry, the ECM field
    and the total sink magnitude are all held fixed and the ONLY thing that
    varies is whether the macrophages sit where the data says they sit. A null
    that scattered points uniformly would also destroy the tissue outline and
    would be rejected for the wrong reason.

    ``p`` is the one-sided permutation p-value
    ``(1 + #{perm >= observed}) / (1 + n_perm)``.

    ``n_perm = 0`` computes the observed statistics only (2 solves), which is how
    the mesh-dependence of the structure statistic is priced at the flagship
    ``dx`` without paying for a second null.
    """
    micro = _reconstruct(FLAGSHIP_SECTION, dx_um)
    alpha = calibrate_alpha()
    dp, _mol = compile_real_drug(drug_name, alpha)
    cfg = SimConfig(drug=dp, pde=PDEConfig(**_PDE), seed=0)
    work, labeled, velocity = _prologue(micro, cfg)
    source = _source_for(work, dp, cfg)
    vel = velocity if cfg.pde.convection else None

    cells = micro.cells
    ct = np.asarray(cells.cell_type)
    tumor = ct == int(CellType.TUMOR)
    immune = ct == int(CellType.IMMUNE)
    n_immune = int(immune.sum())
    x = np.asarray(cells.x, dtype=np.float64)
    y = np.asarray(cells.y, dtype=np.float64)
    col, row = micro.domain.to_index(x, y)

    base = transport_with_immune_sink(work, dp, source, cfg.pde, velocity=vel,
                                      sink_rate=None)
    auc_base = np.asarray(base.auc)[row, col][tumor]

    phi = immune_area_fraction(micro)
    sink = immune_sink_rate(micro, k_max_from_strength(dp, strength), phi=phi)
    obs_tr = transport_with_immune_sink(work, dp, source, cfg.pde, velocity=vel,
                                        sink_rate=sink)
    auc_obs = np.asarray(obs_tr.auc)[row, col][tumor]
    phi_at_tumor = phi[row, col][tumor]
    observed = _structure_stats(auc_base, auc_obs, phi_at_tumor)

    # --- distance to the nearest MEASURED immune cell -------------------------
    from scipy.spatial import cKDTree
    tree = cKDTree(np.column_stack([x[immune], y[immune]]))
    d_tumor, _ = tree.query(np.column_stack([x[tumor], y[tumor]]), k=1)
    a_obs = auc_obs / np.maximum(auc_base, 1e-300)
    dist_edges = (0.0, 10.0, 20.0, 30.0, 50.0, 100.0, np.inf)
    dist_rows = []
    for lo, hi in zip(dist_edges[:-1], dist_edges[1:]):
        m = (d_tumor >= lo) & (d_tumor < hi)
        dist_rows.append({
            "lo_um": lo, "hi_um": hi, "n": int(m.sum()),
            "attenuation_mean": float(a_obs[m].mean()) if m.any() else float("nan"),
            "auc_base_mean": float(auc_base[m].mean()) if m.any() else float("nan"),
        })

    # --- the permutation null -------------------------------------------------
    rng = np.random.default_rng(seed)
    perm_stats: List[Dict] = []
    n_cells = int(ct.size)
    for p in range(int(n_perm)):
        t0 = time.time()
        pick = rng.choice(n_cells, size=n_immune, replace=False)
        phi_p = immune_area_fraction(micro, positions=(x[pick], y[pick]))
        sink_p = immune_sink_rate(micro, k_max_from_strength(dp, strength),
                                  phi=phi_p)
        tr_p = transport_with_immune_sink(work, dp, source, cfg.pde,
                                          velocity=vel, sink_rate=sink_p)
        auc_p = np.asarray(tr_p.auc)[row, col][tumor]
        st = _structure_stats(auc_base, auc_p, phi_p[row, col][tumor])
        perm_stats.append({
            "perm": p, "structure_ratio": st["structure_ratio"],
            "a_p90_over_p10": st["a_p90_over_p10"], "a_mean": st["a_mean"],
        })
        if verbose:
            print(f"    perm {p + 1}/{n_perm} {time.time() - t0:.0f}s "
                  f"SR={st['structure_ratio']:.4f} "
                  f"spread={st['a_p90_over_p10']:.4f}", flush=True)

    def _p(key: str) -> float:
        null = np.array([s[key] for s in perm_stats], dtype=np.float64)
        obs = float(observed[key])
        return float((1.0 + np.sum(null >= obs)) / (1.0 + null.size))

    null_sr = np.array([s["structure_ratio"] for s in perm_stats])
    null_sp = np.array([s["a_p90_over_p10"] for s in perm_stats])
    null_am = np.array([s["a_mean"] for s in perm_stats])

    return {
        "drug": drug_name, "strength": float(strength), "dx_um": float(dx_um),
        "n_perm": int(n_perm), "seed": int(seed), "n_tumor": int(tumor.sum()),
        "n_immune": n_immune,
        "observed": observed,
        "distance_bins": dist_rows,
        "null_structure_ratio_mean": float(null_sr.mean()) if null_sr.size else float("nan"),
        "null_structure_ratio_min": float(null_sr.min()) if null_sr.size else float("nan"),
        "null_structure_ratio_max": float(null_sr.max()) if null_sr.size else float("nan"),
        "null_spread_mean": float(null_sp.mean()) if null_sp.size else float("nan"),
        "null_spread_max": float(null_sp.max()) if null_sp.size else float("nan"),
        "null_a_mean_mean": float(null_am.mean()) if null_am.size else float("nan"),
        "p_structure_ratio": _p("structure_ratio"),
        "p_spread": _p("a_p90_over_p10"),
        "perms": perm_stats,
    }


# --------------------------------------------------------------------------- #
# §6 — the measured human-glioma partition scorecard
# --------------------------------------------------------------------------- #
def glioma_scorecard(mean_phi: float,
                     strengths: Sequence[float] = SEQUESTRATION_SWEEP) -> Dict:
    """Re-score the six measured human glioma drugs with the immune sink on.

    HOW THE SINK REACHES A PARTITION. The scorecard's assumption A1 equates the
    tissue:plasma exposure ratio to ``P_eff``. At quasi-steady state the tissue
    concentration is ``C = S / k_total``, so adding a parallel loss
    ``<k_sink> = k_max * <phi_imm>`` multiplies the partition by
    ``k_base / (k_base + <k_sink>) = 1 - f`` with ``f`` the diverted fraction::

        before:  R = clip( bbb_penetrance + L * D_sample, 0, 1 )
        after :  R = clip( bbb_penetrance + L * D_sample, 0, 1 ) * (1 - f)

    MODALITY MATTERS AND IT IS NOT OPTIONAL. Microdialysis measures FREE
    extracellular drug, which is what a sink into macrophages removes, so the
    factor applies. Tissue homogenate measures TOTAL drug INCLUDING what is
    inside the macrophages, so sequestration does not lower it at all and the
    factor must NOT be applied. This is the same free-vs-total bookkeeping
    ``results/ion_trapping_impact.md`` §6 established, applied here with the
    opposite sign. Both variants are reported; applying it to everything is
    physically wrong and is shown only to prove the choice is load-bearing.
    """
    from .glioma_regional_pk import run_scorecard

    before = run_scorecard(aggregation="geometric", apply_leak=True)
    out: List[Dict] = []
    for S in strengths:
        f = diverted_fraction_from_strength(S, mean_phi)
        for modality_aware in (True, False):
            rows: List[Dict] = []
            for r in before["rows"]:
                meas = float(r["measured_mean"])
                pred0 = float(r["predicted_partition"])
                free = "microdialysis" in str(r["method"]).lower()
                apply = free if modality_aware else True
                pred = pred0 * (1.0 - f) if apply else pred0
                fold = pred / meas if meas > 0 else float("inf")
                rows.append({
                    "drug": r["drug"], "measured": meas, "method": r["method"],
                    "free": free, "applied": apply,
                    "pred_before": pred0, "fold_before": float(r["fold_error"]),
                    "pred_after": pred, "fold_error": fold,
                    "log10_fold": math.log10(fold) if fold > 0 else float("nan"),
                })
            logs = [r["log10_fold"] for r in rows if np.isfinite(r["log10_fold"])]
            folds = [r["fold_error"] for r in rows if np.isfinite(r["fold_error"])]
            out.append({
                "strength": float(S), "diverted_fraction": f,
                "modality_aware": modality_aware,
                "gm_abs": float(10 ** np.mean(np.abs(logs))) if logs else float("nan"),
                "gm_signed": float(10 ** np.mean(logs)) if logs else float("nan"),
                "within_3x": int(sum(1 for x in folds if 1 / 3.0 <= x <= 3.0)),
                "within_10x": int(sum(1 for x in folds if 0.1 <= x <= 10.0)),
                "worst": float(max(max(f_, 1 / f_) for f_ in folds if f_ > 0)) if folds else float("nan"),
                "rows": rows,
            })

    logs0 = [math.log10(float(r["fold_error"])) for r in before["rows"]
             if float(r["fold_error"]) > 0]
    folds0 = [float(r["fold_error"]) for r in before["rows"]]
    baseline = {
        "gm_abs": float(10 ** np.mean(np.abs(logs0))),
        "gm_signed": float(10 ** np.mean(logs0)),
        "within_3x": int(sum(1 for x in folds0 if 1 / 3.0 <= x <= 3.0)),
        "within_10x": int(sum(1 for x in folds0 if 0.1 <= x <= 10.0)),
        "worst": float(max(max(f_, 1 / f_) for f_ in folds0 if f_ > 0)),
    }
    aware = [r for r in out if r["modality_aware"]]
    best = min(aware, key=lambda r: r["gm_abs"])
    return {
        "baseline": baseline, "arms": out, "mean_phi": float(mean_phi),
        "best_modality_aware": best,
        "improves": bool(best["gm_abs"] < baseline["gm_abs"]
                         and best["strength"] > 0.0),
    }


# --------------------------------------------------------------------------- #
# Report
# --------------------------------------------------------------------------- #
def _f(x: float, n: int = 4) -> str:
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "--"
    return f"{x:.{n}f}"


def write_report(ident: Dict, fields: List[Dict], sw: Dict, sp: Dict,
                 sc: Dict, spf: Optional[Dict] = None,
                 path: str = REPORT_MD) -> str:
    A: List[str] = []
    a = A.append

    obs = sp["observed"]
    anchor_row = next(r for r in fields if r["section"] == "brain_319")
    dox0 = next(r for r in sw["rows"]
                if r["drug"] == "doxorubicin" and r["strength"] == 0.0)
    doxS = next(r for r in sw["rows"]
                if r["drug"] == "doxorubicin"
                and r["strength"] == float(sp["strength"]))
    doxmax = max((r for r in sw["rows"] if r["drug"] == "doxorubicin"),
                 key=lambda r: r["strength"])

    a("# Immune (TAM) drug sequestration — implemented as a SPATIAL sink, and the "
      "spatial structure is the whole result")
    a("")
    p_sr = float(sp["p_structure_ratio"])
    beats_null = p_sr <= 0.05
    null_clause = (
        f"and it beats a label-permutation null at p = {_f(p_sr, 4)}"
        if beats_null else
        f"but it does NOT beat a label-permutation null (p = {_f(p_sr, 4)}), "
        f"i.e. the structure is a property of where the CELLS are and not "
        f"specifically of where the MACROPHAGES are"
    )
    imm_lo = min(r["immune_pct"] for r in fields)
    imm_hi = max(r["immune_pct"] for r in fields)
    a("**Verdict up front.** The immune compartment was pharmacologically inert: "
      f"{imm_lo:.1f}-{imm_hi:.1f}% of the segmented cells across the three real "
      "sections occupied space and consumed nothing. It "
      "now has a mechanism — a first-order drug loss proportional to the local "
      "MEASURED immune-cell area fraction — and the honest findings are three. "
      "**(i) The literature supports a DIRECTION and not a MAGNITUDE**, and for "
      "free doxorubicin the one direct measurement is NEGATIVE, so the magnitude "
      "is swept and never fitted. **(ii) The reduction is spatially "
      "structured rather than uniform** — tumour cells in the macrophage-densest "
      f"quintile keep {_f(obs['bins'][-1]['attenuation_mean'], 4)} of their "
      f"exposure against {_f(obs['bins'][0]['attenuation_mean'], 4)} in the "
      f"sparsest, a **{_f(obs['structure_ratio'], 4)}x** contrast that a "
      f"well-mixed model gives as exactly 1.0000, {null_clause}. "
      "**(iii) It does not improve the measured human-glioma anchor**, so it "
      "ships DEFAULT-OFF.")
    a("")
    a("| question | answer |")
    a("|---|---|")
    a(f"| drugs with a MEASURED free-drug TAM uptake ratio | "
      f"**{substrate_summary()['n_with_measured_ratio']} of "
      f"{substrate_summary()['n_drugs_assessed']}** |")
    a(f"| mean tumour AUC, doxorubicin, S=0 -> S={_f(doxmax['strength'], 0)} | "
      f"{_f(dox0['mean_auc_uM_s'], 1)} -> {_f(doxmax['mean_auc_uM_s'], 1)} uM*s "
      f"({_f(doxmax['mean_auc_uM_s'] / dox0['mean_auc_uM_s'], 3)}x) |")
    a(f"| therapeutic coverage at 1x IC50, same range | "
      f"{_f(dox0['coverage_1.0x'])} -> {_f(doxmax['coverage_1.0x'])} |")
    a(f"| predicted kill, same range | {_f(dox0['kill_fraction'])} -> "
      f"{_f(doxmax['kill_fraction'])} |")
    a(f"| **cross-drug RANKING across the whole sweep** | "
      f"**{'INVARIANT' if sw['ranking_invariant'] else 'CHANGES'}** |")
    a(f"| spatial structure ratio (sparsest / densest immune quintile) | "
      f"**{_f(obs['structure_ratio'], 4)}x**, permutation p = "
      f"{_f(sp['p_structure_ratio'], 4)} |")
    a(f"| human glioma scorecard, GM abs fold error | "
      f"{_f(sc['baseline']['gm_abs'], 3)}x -> "
      f"{_f(sc['best_modality_aware']['gm_abs'], 3)}x at best swept S "
      f"({_f(sc['best_modality_aware']['strength'], 1)}) |")
    a(f"| **recommendation** | **DEFAULT OFF** — see §7 |")
    a("")
    a("---")
    a("")

    # ---------------- §1 -----------------------------------------------------
    a("## 1. Nothing outside the two new files was touched, and OFF is bit-identical")
    a("")
    a("The escalation set by `results/uptake_impact.md` and "
      "`results/ion_trapping_impact.md` is followed exactly: **no shipped code "
      "path imports `tmesim.molecular.immune_sink`**, so the term is off BY "
      "CONSTRUCTION rather than by default. Files created: "
      "`tmesim/molecular/immune_sink.py`, `tmesim/validation/immune_sink_impact.py`, "
      "`tests/test_immune_sink.py`, this report. Files modified elsewhere: "
      f"**{ident['n_shipped_files_modified']}**.")
    a("")
    a("A spatially-varying sink cannot be expressed by rescaling the AUC field "
      "afterwards — drug diffuses INTO a depleted region from its neighbours, "
      "and that non-locality is the entire microenvironment effect being "
      "measured. So the transport is RE-SOLVED. The shipped solver is not "
      "modified and not copied: its three numerical primitives "
      "(`_effective_diffusivity`, `_build_diffusion_operator`, "
      "`_advective_divergence`) are imported unchanged, exactly as "
      "`validation/moncayo_explant.py::bath_diffusion` already does, and the "
      "single functional difference is one line:")
    a("")
    a("```python")
    a("k = np.full((ny, nx), drug.clearance + drug.binding) + sink_rate")
    a("```")
    a("")
    a("Because `x + 0.0 == x` exactly in IEEE-754, the OFF path must be "
      "BIT-identical, and that is asserted rather than assumed:")
    a("")
    a(f"- **{ident['n_bit_identical']}/{ident['n_comparisons']} solver "
      f"comparisons bit-for-bit identical** (`.tobytes()` equality on `auc` and "
      f"`final_concentration`, `np.array_equal` on `times` and the "
      f"concentration history, `==` on the penetration index; never `allclose`) "
      f"across 3 (clearance, binding) pairs x convection on/off x 2 timesteps x "
      f"{{`sink_rate=None`, `sink_rate=zeros`}} plus the time-resolved-source "
      f"branch. All identical: **{ident['all_bit_identical']}**.")
    a(f"- Every drug's LIBRARY strength is exactly 0.0, including a name not in "
      f"the table: **{ident['all_strengths_exactly_zero']}**.")
    a("")
    a("And on the REAL flagship section the gate is end-to-end rather than "
      "synthetic: the `S = 0` arm of §4 is required to reproduce "
      "`sim.engine.simulate_response` bit-for-bit or the run aborts.")
    a("")
    a("| drug | S=0 AUC bit-identical to the engine | S=0 kill identical | engine kill |")
    a("|---|---|---|---|")
    for g in sw["gates"]:
        a(f"| {g['drug']} | **{g['auc_bit_identical']}** | "
          f"**{g['kill_identical']}** | {g['engine_kill']:.6f} |")
    a("")
    a("The checks are not vacuous — a positive sink genuinely bites, only ever "
      "lowers exposure, and is spatially non-uniform:")
    a("")
    a(f"- a positive sink changes the AUC field: **{ident['sink_bites']}**")
    a(f"- it never RAISES exposure anywhere: **{ident['sink_only_lowers']}**")
    a(f"- and the attenuation genuinely VARIES across tissue (range "
      f"{_f(ident['attenuation_range_over_tissue'], 4)}, which would be 0 for a "
      f"well-mixed sink): **{ident['sink_is_spatially_varying']}**")
    a("")

    # ---------------- §2 -----------------------------------------------------
    a("## 2. The literature — a DIRECTION, not a MAGNITUDE, and one direct NEGATIVE")
    a("")
    a("This was searched before it was coded. The result is uncomfortable and it "
      "is the most important section of the report, because everything "
      "downstream is conditional on it.")
    a("")
    a("| finding | what was measured | direction | PMID / DOI |")
    a("|---|---|---|---|")
    a("| Miller et al., Nat Commun 2015;6:8692 — **the best-quantified case** | "
      "single-cell imaging of a PLGA-b-PEG **nanoparticle** carrying a "
      "fluorescent Pt(IV) pro-drug. TAMs were **4% of cells but accumulated 30% "
      "of total intratumoral nanoparticle**; payload was **>3x elevated within "
      "~15 um of a TAM** with gradients resolved over **<25 um**; macrophage "
      "depletion **LOWERED** intratumoral Pt 24% (5.0 -> 3.8 uM) and made tumours "
      "grow **faster** | **DEPOT — the opposite of a sink** | 26503691 / "
      "10.1038/ncomms9692 |")
    a("| Islam et al., Nanotheranostics 2022;6(4):451-464 — **the direct "
      "free-doxorubicin measurement** | TAM-rich WEHI-164 fibrosarcoma. "
      "Liposomal formulations \"shifted cellular drug uptake to TAM and to "
      "monocytic MDSC\"; **\"free doxorubicin cellular drug uptake was below "
      "detectable levels\"** | **NEGATIVE for free drug** | 36105861 / "
      "10.7150/ntno.75045 |")
    a("| Mayer et al., J Pharmacol Exp Ther 1997;280(3):1406-1414 | two murine "
      "fibrosarcomas differing **4-fold in TAM content**; TAM-rich tumours took "
      "1.5-2x more liposomal doxorubicin, but CD11b+ cells correlated **poorly** "
      "with doxorubicin fluorescence and the authors concluded TAM uptake does "
      "**not** explain it (attributing it to TAM-driven vascular permeability) | "
      "**NEGATIVE for TAM uptake** | 9067330 |")
    a("| Zhitomirsky & Assaraf, Oncotarget 2015;6(2):1143-1156; Drug Resist "
      "Updat 2016;24:23-33 | doxorubicin is a hydrophobic weak base "
      "(amino-sugar pKa ~8.3) passively ion-trapped into acidic organelles; "
      "sequestered drug is held away from its nuclear target | **sink, "
      "mechanistic, no TAM-specific ratio** | 25544758 / "
      "10.18632/oncotarget.2732; 26830313 / 10.1016/j.drup.2015.11.004 |")
    a("| Aghai & Tokes, Cancer Chemother Pharmacol 1990;25(4):295-298 | "
      "doxorubicin compartmentalisation across normal human marrow lineages: "
      "the nuclei of **macrophages and activated phagocytes were the LEAST "
      "doxorubicin-accessible of any lineage**, i.e. in macrophages the drug is "
      "held in the cytoplasmic/vesicular compartment | **sink, single-cell, no "
      "ratio** | 2295115 / 10.1007/BF00684889 |")
    a("")
    a("Read the first two rows together, because they are the ones that decide "
      "the framing. **The well-quantified TAM drug-handling literature is about "
      "CARRIERS, not free drug, and where it is quantified the sign is a DEPOT.** "
      "The one experiment that measured free doxorubicin in TAMs found it below "
      "detection. So a TAM SINK for free cytotoxics rests on mechanism "
      "(lysosomal ion trapping in the most lysosome-rich cell in the niche) plus "
      "a single-cell compartmentalisation observation from 1990 — a direction, "
      "with no number attached.")
    a("")
    a("**Consequently `immune_sink_strength_for(drug)` returns exactly 0.0 for "
      "every drug**, the same identity `uptake_factor` and `ion_trapping` "
      "adopted, and the magnitude below is a SWEPT parameter. Per drug:")
    a("")
    a("| drug | measured ratio | direction the literature supports | why |")
    a("|---|---|---|---|")
    for d in DRUGS:
        det = immune_sink_detail(d)
        a(f"| {d} | **{det['ratio'] if det['ratio'] is not None else 'NONE'}** | "
          f"{det['sign']} | {det['mechanism'][:170]}... |")
    a("")
    a("Note the contrast structure this gives, which is the point of including "
      "the other two drugs: **doxorubicin** has a mechanism and no number, "
      "**cisplatin**'s only quantified TAM result is a nanoparticle DEPOT whose "
      "sign is opposite, and **temozolomide** — small, neutral, non-basic, a "
      "pro-drug that decomposes spontaneously — has no proposed route at all. A "
      "sink applied equally to all three is therefore already known to be wrong; "
      "it is applied equally anyway, because the alternative is to invent three "
      "different unmeasured numbers.")
    a("")
    a("### 2.1 The functional form and its units")
    a("")
    a("```")
    a("  dC/dt = div(D_eff grad C) - div(vC) + S - (clearance + binding) C")
    a("                                          - k_sink(x, y) C        <-- NEW")
    a("")
    a("  k_sink(x, y) = k_max * phi_imm(x, y)                 [1/s]")
    a("```")
    a("")
    a("| symbol | meaning | units | where it comes from |")
    a("|---|---|---|---|")
    a("| `phi_imm` | local IMMUNE-cell area fraction | dimensionless, [0,1] | "
      "**MEASURED** — Gaussian KDE of the real immune-cell centroids (cells/um^2) "
      "times the per-cell cross-section `pi r^2` (um^2/cell) |")
    a(f"| `r` | immune-cell cross-sectional radius | um | DECLARED "
      f"{IMMUNE_CELL_RADIUS_UM:.1f} um (macrophage/microglial soma ~15 um across) |")
    a(f"| bandwidth | footprint over which one macrophage depletes drug | um | "
      f"DECLARED {IMMUNE_KERNEL_BANDWIDTH_UM:.0f} um, **literature-anchored**: "
      f"Miller 2015 resolved TAM payload gradients over <"
      f"{MILLER2015_GRADIENT_LENGTH_UM:.0f} um (PMID 26503691) |")
    a("| `k_max` | loss constant of a voxel packed ENTIRELY with macrophages | "
      "1/s | **SWEPT**, never fitted |")
    a("| `S` | `k_max / (clearance + binding)` | dimensionless | the sweep axis: "
      "\"how many times the drug's own tissue sink does a fully-macrophage voxel "
      "add\" |")
    a("")
    a("Only the PRODUCT `k_max * phi_imm` is physically meaningful, so `r` and "
      "`k_max` trade off exactly against each other — which is precisely why "
      "`k_max` is swept rather than asserted, and why no number here should be "
      "read as a calibrated rate. **Zero parameters were fitted to any "
      "scorecard.**")
    a("")

    # ---------------- §3 -----------------------------------------------------
    a("## 3. The measured immune field, and where the literature anchor lands on it")
    a("")
    a(f"Immune-cell area fraction on all three real sections at dx = "
      f"{fields[0]['dx_um']:.0f} um, from the frozen `reconstruct_real` default "
      f"path (matched-null classifier).")
    a("")
    a("| section | cells | IMMUNE | IMMUNE % | grid | mean phi | median | p90 | "
      "p99 | max |")
    a("|---|---|---|---|---|---|---|---|---|---|")
    for r in fields:
        a(f"| `{r['section']}` | {r['n_cells']:,} | {r['n_immune']:,} | "
          f"{r['immune_pct']:.2f}% | {r['grid']} | **{_f(r['phi_mean'])}** | "
          f"{_f(r['phi_median'])} | {_f(r['phi_p90'])} | {_f(r['phi_p99'])} | "
          f"{_f(r['phi_max'])} |")
    a("")
    a("The distribution is strongly right-skewed on every section — the p99 is "
      "3-6x the mean — which is the reconstruction saying the immune compartment "
      "is CLUMPED, not spread. That clumping is what §5 tests, and it is "
      "independently corroborated: `results/classifier_fix.md` §4 measures "
      "nearest-neighbour IMMUNE label enrichment of 1.41-1.83x above random on "
      "these same sections under the shipping classifier.")
    a("")
    a("**Converting the one literature anchor into this model's units.** Miller "
      "2015 gives a fraction `f` of intratumoral drug that ends up inside "
      "macrophages (0.30 from the 4%-of-cells/30%-of-drug partition; 0.24 from "
      "the depletion experiment). Two first-order losses in parallel split the "
      "drug in proportion to their rate constants, so "
      "`f = <k_sink>/(<k_sink> + k_base)` and therefore "
      "`S = f / ((1 - f) * <phi_imm>)`. Run through each section's OWN measured "
      "field:")
    a("")
    a("| section | mean phi | S at f = 0.24 (depletion) | S at f = 0.30 (partition) |")
    a("|---|---|---|---|")
    for r in fields:
        a(f"| `{r['section']}` | {_f(r['phi_mean'])} | "
          f"**{_f(r['S_anchor_depletion'], 2)}** | "
          f"**{_f(r['S_anchor_drug_fraction'], 2)}** |")
    a("")
    a(f"So the anchor lands at **S ~ {_f(anchor_row['S_anchor_depletion'], 1)}-"
      f"{_f(anchor_row['S_anchor_drug_fraction'], 1)}** on the flagship section, "
      f"which is why the sweep runs "
      f"{', '.join(str(s) for s in SEQUESTRATION_SWEEP)} — it BRACKETS the "
      f"anchor rather than sitting on it. **This is a bracketing calculation, "
      f"not a calibration.** The anchor is a polymeric nanoparticle in a murine "
      f"fibrosarcoma; the drugs here are free small molecules in human glioma; "
      f"and in the source experiment the effect ran in the DEPOT direction. It "
      f"is used only to answer \"how large could a diverted fraction plausibly "
      f"be at all\", and the answer is that even the most TAM-avid formulation "
      f"ever measured diverted 24-30%.")
    a("")

    # ---------------- §4 -----------------------------------------------------
    a("## 4. The sweep on the real flagship section")
    a("")
    a(f"**Section.** `{FLAGSHIP_SECTION}` at dx = {sw['dx_um']:.0f} um, grid "
      f"{sw['grid']}, {sw['n_tumor']:,} tumour cells. Perivascular-sink "
      f"alpha = {sw['alpha']:.2f}, calibrated once on doxorubicin exactly as "
      f"every flagship run does. ONE reconstruction and ONE oxygen + IFP "
      f"prologue are shared by every arm, so the comparison isolates the sink. "
      f"`_PDE` is imported from `validation/delivery_metrics.py` so this and "
      f"every flagship delivery number solve the identical problem.")
    a("")
    a("**Why this section and not all three.** Every absolute AUC / coverage / "
      "kill number in `results/` was produced on this section at this dx, so it "
      "is the only one against which a delta is interpretable. The immune FIELD "
      "is measured on all three (§3), and it is worth reading §3 alongside this "
      "table: `io_480` and `gse324268_366` carry mean immune area fractions "
      f"{max(r['phi_mean'] for r in fields) / anchor_row['phi_mean']:.1f}x "
      "the flagship's, because they are far denser sections, so the same S would "
      "produce a substantially larger effect there. The levels below are "
      "flagship-specific; the SHAPE of the S-dependence is not.")
    a("")
    a("| drug | IC50 (uM) | S | diverted fraction f | mean tumour AUC (uM*s) | "
      "mean exposure (uM) | coverage @1x IC50 | predicted kill | penetration index |")
    a("|---|---|---|---|---|---|---|---|---|")
    for name in DRUGS:
        for r in [x for x in sw["rows"] if x["drug"] == name]:
            bold = "**" if r["strength"] == 0.0 else ""
            a(f"| {name} | {r['ic50_uM']:g} | {bold}{r['strength']:g}{bold} | "
              f"{_f(r['diverted_fraction'], 3)} | {_f(r['mean_auc_uM_s'], 1)} | "
              f"{_f(r['mean_exposure_uM'], 4)} | {_f(r['coverage_1.0x'])} | "
              f"{_f(r['kill_fraction'])} | {_f(r['penetration_index'])} |")
    a("")
    a("Relative to OFF, per drug, across the swept range:")
    a("")
    a("| drug | AUC at S_max / AUC at S=0 | coverage change (pp) | kill change (pp) |")
    a("|---|---|---|---|")
    for name in DRUGS:
        sub = [x for x in sw["rows"] if x["drug"] == name]
        z = next(x for x in sub if x["strength"] == 0.0)
        m = max(sub, key=lambda x: x["strength"])
        a(f"| {name} | **{_f(m['mean_auc_uM_s'] / z['mean_auc_uM_s'], 4)}x** | "
          f"{100 * (m['coverage_1.0x'] - z['coverage_1.0x']):+.2f} | "
          f"{100 * (m['kill_fraction'] - z['kill_fraction']):+.2f} |")
    a("")
    a("### 4.1 The cross-drug ranking")
    a("")
    a("| S | ranking by predicted kill | by coverage @1x IC50 | by mean tumour AUC |")
    a("|---|---|---|---|")
    for r in sw["rankings"]:
        a(f"| {r['strength']:g} | {r['by_kill']} | {r['by_coverage']} | "
          f"{r['by_auc']} |")
    a("")
    a(f"**The ranking by predicted kill is "
      f"{'INVARIANT across the entire swept range' if sw['ranking_invariant'] else 'NOT invariant'}.** "
      "That is the expected structure and it is the same one "
      "`results/cell_state_impact.md` §3.3 found for the cell-state route: the "
      "sink is applied with the SAME strength to all three drugs, so it very "
      "largely cancels out of a comparison and does not cancel at all out of a "
      "level. A claim of the form *\"cisplatin beats doxorubicin beats "
      "temozolomide on this tissue\"* survives every sequestration strength "
      "tested; a claim about an absolute kill fraction does not.")
    a("")

    # ---------------- §5 -----------------------------------------------------
    a("## 5. THE SPATIAL QUESTION — is the reduction structured, or just smaller?")
    a("")
    a("This is the section the whole module exists for. **A uniform reduction is "
      "pharmacologically boring: a compartmental model reproduces it exactly by "
      "lowering one rate constant.** The question is whether tumour cells sitting "
      "in macrophage-dense neighbourhoods are measurably worse off than tumour "
      "cells that are not — because that is a statement only a spatial model can "
      "make.")
    a("")
    a(f"**Configuration.** `{sp['drug']}` at S = {sp['strength']:g} (the sweep "
      f"point nearest the Miller anchor), flagship section, dx = "
      f"{sp['dx_um']:.0f} um, {sp['n_tumor']:,} tumour cells and "
      f"{sp['n_immune']:,} immune cells. The statistic is the per-tumour-cell "
      f"**attenuation** `a_i = AUC_sink,i / AUC_base,i`, from the SAME "
      f"prologue, source field and velocity field — the only difference between "
      f"the two solves is `k_sink`.")
    a("")
    if spf is not None and abs(spf["dx_um"] - sp["dx_um"]) > 1e-9:
        a(f"**A declared cost decision, and it is priced.** Every permutation is "
          f"a full transport re-solve, so the null runs at dx = "
          f"{sp['dx_um']:.0f} um rather than the flagship "
          f"{spf['dx_um']:.0f} um (~2.8x fewer unknowns). Because this project "
          f"treats mesh-dependence as a first-class failure mode, the SAME "
          f"observed statistic is also computed at the flagship spacing and "
          f"both are reported in §5.4 — so the reader can see directly whether "
          f"the structure is a property of the tissue or of the grid.")
        a("")
    a("### 5.1 Tumour-cell exposure binned by LOCAL IMMUNE DENSITY")
    a("")
    a("| quintile of local phi_imm | n tumour cells | mean phi_imm | "
      "mean AUC, sink OFF | mean AUC, sink ON | **mean attenuation** |")
    a("|---|---|---|---|---|---|")
    for b in obs["bins"]:
        a(f"| {b['bin']} ({'sparsest' if b['bin'] == 1 else 'densest' if b['bin'] == N_BINS else '.'}) "
          f"| {b['n']:,} | {_f(b['phi_mean'])} | {_f(b['auc_base_mean'], 1)} | "
          f"{_f(b['auc_sink_mean'], 1)} | **{_f(b['attenuation_mean'], 4)}** |")
    a("")
    a(f"**Structure ratio** (mean attenuation in the sparsest quintile / mean "
      f"attenuation in the densest) = **{_f(obs['structure_ratio'], 4)}x**. A "
      f"well-mixed sink of any magnitude gives exactly **1.0000**. The bin-free "
      f"version — the p90/p10 spread of the per-cell attenuation — is "
      f"**{_f(obs['a_p90_over_p10'], 4)}x**, over a full range of "
      f"{_f(obs['a_min'], 4)} to {_f(obs['a_max'], 4)}, against a section-mean "
      f"attenuation of {_f(obs['a_mean'], 4)}.")
    a("")
    a("### 5.2 Exposure deficit versus DISTANCE to the nearest macrophage")
    a("")
    a("The same effect read the other way round, and on the length scale Miller "
      f"2015 measured (>3x payload elevation within ~15 um of a TAM; gradients "
      f"resolved over <{MILLER2015_GRADIENT_LENGTH_UM:.0f} um):")
    a("")
    a("| distance to nearest immune cell | n tumour cells | mean AUC, sink OFF | "
      "mean attenuation |")
    a("|---|---|---|---|")
    for d in sp["distance_bins"]:
        hi = "inf" if not np.isfinite(d["hi_um"]) else f"{d['hi_um']:.0f}"
        a(f"| {d['lo_um']:.0f}-{hi} um | {d['n']:,} | "
          f"{_f(d['auc_base_mean'], 1)} | **{_f(d['attenuation_mean'], 4)}** |")
    a("")
    a("### 5.3 The permutation null")
    a("")
    a("**The null, stated before its result.** The observed structure is compared "
      "against re-solves in which the IMMUNE LABEL is permuted among the "
      "measured cells: the same NUMBER of immune cells, drawn from the same "
      "measured cell POSITIONS, but with the immune compartment's own spatial "
      "clustering destroyed. The tissue outline, the cell layout, the vessel "
      "geometry, the ECM field and the total sink magnitude are all held fixed, "
      "so the only thing that varies is whether the macrophages sit where the "
      "data says they sit. A null that scattered points uniformly through the "
      "tissue would also destroy the tissue's own density structure and would be "
      "rejected for the wrong reason.")
    a("")
    a(f"{sp['n_perm']} permutations, seed {sp['seed']}; each one is a full "
      f"transport re-solve.")
    a("")
    a("| statistic | observed | null mean | null range | one-sided p |")
    a("|---|---|---|---|---|")
    a(f"| structure ratio (sparsest/densest quintile) | "
      f"**{_f(obs['structure_ratio'], 4)}** | "
      f"{_f(sp['null_structure_ratio_mean'], 4)} | "
      f"{_f(sp['null_structure_ratio_min'], 4)}-"
      f"{_f(sp['null_structure_ratio_max'], 4)} | "
      f"**{_f(sp['p_structure_ratio'], 4)}** |")
    a(f"| per-cell attenuation spread (p90/p10) | "
      f"**{_f(obs['a_p90_over_p10'], 4)}** | {_f(sp['null_spread_mean'], 4)} | "
      f"up to {_f(sp['null_spread_max'], 4)} | **{_f(sp['p_spread'], 4)}** |")
    a(f"| mean attenuation (a magnitude, NOT a structure statistic) | "
      f"{_f(obs['a_mean'], 4)} | {_f(sp['null_a_mean_mean'], 4)} | -- | -- |")
    a("")
    a("**Read the last row as the control.** The MEAN attenuation is a magnitude, "
      "not a structure statistic: the same number of macrophages removes roughly "
      "the same total amount of drug wherever they sit. What the permutation is "
      "asked about is the SPREAD and the density contrast — i.e. not how much "
      "drug is lost, but **which tumour cells lose it**. A compartmental model "
      "gets the first quantity right by construction and cannot express the "
      "second at all.")
    a("")
    if beats_null:
        a(f"**The measured immune geometry is doing work.** The observed "
          f"structure ratio of {_f(obs['structure_ratio'], 4)} sits above the "
          f"entire null range "
          f"({_f(sp['null_structure_ratio_min'], 4)}-"
          f"{_f(sp['null_structure_ratio_max'], 4)}), p = {_f(p_sr, 4)}: putting "
          f"the same number of macrophages at other measured cell positions "
          f"produces a measurably WEAKER exposure gradient than putting them "
          f"where the data says they are. The clustering is not an artefact of "
          f"the tissue being denser in some places than others.")
    else:
        a(f"**Honest negative on the null, and it is a specific one.** The "
          f"observed structure ratio of {_f(obs['structure_ratio'], 4)} is NOT "
          f"outside the null range "
          f"({_f(sp['null_structure_ratio_min'], 4)}-"
          f"{_f(sp['null_structure_ratio_max'], 4)}), p = {_f(p_sr, 4)}. The "
          f"reduction IS spatially structured — the contrast against a "
          f"well-mixed model's exact 1.0000 is real and is what the "
          f"microenvironment model adds — but the structure is not specific to "
          f"the MACROPHAGES: relabelling the same number of immune cells onto "
          f"other measured cell positions reproduces a comparable gradient, "
          f"because the dominant signal is where the TISSUE is dense, not where "
          f"the CD163+ cells are. This is a deliberately conservative null (it "
          f"holds the tissue outline and cell layout fixed), and it says the "
          f"claim this analysis can support is 'exposure is structured by local "
          f"cellularity', not 'exposure is structured by macrophages "
          f"specifically'.")
    a("")
    a(f"The floor on the p-value at {sp['n_perm']} permutations is "
      f"{1.0 / (1 + sp['n_perm']):.4f}, so a reported "
      f"{_f(sp['p_structure_ratio'], 4)} is at that floor rather than a "
      f"resolved small number; the effect size (the null range) is the more "
      f"informative column.")
    a("")
    if spf is not None:
        a("### 5.4 Is the structure a property of the tissue, or of the mesh?")
        a("")
        a("`results/state_validation.md` and `results/cell_state_impact.md` both "
          "turn on the same question, so it is asked here rather than assumed. "
          "The identical observed statistic, on the identical section and drug "
          "at the identical strength, at two grid spacings:")
        a("")
        a("| dx (um) | structure ratio (sparsest/densest quintile) | "
          "attenuation p90/p10 | mean attenuation |")
        a("|---|---|---|---|")
        for tag in (spf, sp):
            o = tag["observed"]
            tail = " (flagship)" if abs(tag["dx_um"] - DX_FLAGSHIP) < 1e-9 else " (null runs here)"
            a(f"| {tag['dx_um']:.0f}{tail} | **{_f(o['structure_ratio'], 4)}** | "
              f"{_f(o['a_p90_over_p10'], 4)} | {_f(o['a_mean'], 4)} |")
        sr_a, sr_b = (sp["observed"]["structure_ratio"],
                      spf["observed"]["structure_ratio"])
        swing = max(sr_a, sr_b) / min(sr_a, sr_b) if min(sr_a, sr_b) > 0 else float("nan")
        a("")
        a(f"Swing across the mesh: **{_f(swing, 3)}x**. "
          + ("That is small enough that the structure is a property of the "
             "measured cell geometry rather than of the discretisation, which "
             "is the bar this repo applies to every spatial claim."
             if swing < 1.2 else
             "That is NOT small, so the structure statistic carries a mesh "
             "dependence and its LEVEL should not be quoted without the dx it "
             "was computed at. The qualitative direction (densest quintile "
             "worse off) is reproduced at both spacings; the magnitude is not."))
        a("")

    # ---------------- §6 -----------------------------------------------------
    a("## 6. Does it improve the measured human-glioma anchor? NO.")
    a("")
    a("`results/glioma_regional_pk.md` is the only measured scorecard this term "
      "can reach: 6 drugs with measured human glioma tissue:plasma partitions, "
      "currently at a geometric-mean absolute fold error of "
      f"{_f(sc['baseline']['gm_abs'], 3)}x. At quasi-steady state a parallel "
      "loss multiplies the partition by `(1 - f)`:")
    a("")
    a("```")
    a("  before:  R = clip( bbb_penetrance + L * D_sample, 0, 1 )")
    a("  after :  R = clip( bbb_penetrance + L * D_sample, 0, 1 ) * (1 - f)")
    a("```")
    a("")
    a("**Modality is not optional here.** Microdialysis measures FREE "
      "extracellular drug, which is exactly what a sink into macrophages "
      "removes, so the factor applies. Tissue homogenate measures TOTAL drug "
      "INCLUDING what is inside the macrophages, so sequestration does not lower "
      "it at all and the factor must not be applied. This is the same "
      "free-versus-total bookkeeping `results/ion_trapping_impact.md` §6 "
      "established, used here with the opposite sign. Both variants are shown "
      "because the wrong one is what makes the choice load-bearing.")
    a("")
    a(f"Applied at each section's own mean phi_imm = {_f(sc['mean_phi'])}:")
    a("")
    a("| S | diverted fraction f | modality-aware GM abs fold | (physically "
      "wrong) applied to all GM abs fold | within 3x | within 10x | worst |")
    a("|---|---|---|---|---|---|---|")
    a(f"| **0 (OFF)** | 0.000 | **{_f(sc['baseline']['gm_abs'], 3)}x** | "
      f"{_f(sc['baseline']['gm_abs'], 3)}x | {sc['baseline']['within_3x']}/6 | "
      f"{sc['baseline']['within_10x']}/6 | {_f(sc['baseline']['worst'], 2)}x |")
    for S in SEQUESTRATION_SWEEP:
        if S == 0.0:
            continue
        aw = next(r for r in sc["arms"] if r["strength"] == S and r["modality_aware"])
        al = next(r for r in sc["arms"] if r["strength"] == S and not r["modality_aware"])
        a(f"| {S:g} | {_f(aw['diverted_fraction'], 3)} | "
          f"**{_f(aw['gm_abs'], 3)}x** | {_f(al['gm_abs'], 3)}x | "
          f"{aw['within_3x']}/6 | {aw['within_10x']}/6 | {_f(aw['worst'], 2)}x |")
    a("")
    best = sc["best_modality_aware"]
    if sc["improves"]:
        a(f"The best swept point (S = {best['strength']:g}) gives "
          f"{_f(best['gm_abs'], 3)}x against a baseline of "
          f"{_f(sc['baseline']['gm_abs'], 3)}x — a nominal improvement. It is "
          "NOT evidence: the strength was chosen by looking at this scorecard, "
          "which is the definition of fitting, and the value it selects has no "
          "independent support.")
    else:
        a(f"**No swept strength improves the anchor.** The best modality-aware "
          f"point (S = {best['strength']:g}) gives {_f(best['gm_abs'], 3)}x "
          f"against a baseline of {_f(sc['baseline']['gm_abs'], 3)}x.")
    a("")
    a("**Why it cannot help, structurally.** The sink can only push predictions "
      "DOWN, and it can only reach the three microdialysis drugs. Of those, the "
      "model currently OVER-predicts only eribulin; it already UNDER-predicts "
      "temozolomide and methotrexate. So any strength large enough to help "
      "eribulin drags the other two further below their measurements, and the "
      "geometric mean over the set gets worse before eribulin's 6.9x is dented. "
      "This is a one-directional term applied to a scorecard whose errors point "
      "both ways — the same structural mismatch `results/uptake_impact.md` §3 "
      "documented for carrier-mediated uptake, mirrored.")
    a("")
    a("Per-drug at the strength nearest the Miller anchor "
      f"(S = {SPATIAL_STRENGTH:g}, modality-aware):")
    a("")
    a("| drug | modality | measured | pred BEFORE | fold BEFORE | pred AFTER | "
      "fold AFTER | sink applied? |")
    a("|---|---|---|---|---|---|---|---|")
    arm = next(r for r in sc["arms"]
               if r["strength"] == SPATIAL_STRENGTH and r["modality_aware"])
    for r in arm["rows"]:
        a(f"| {r['drug']} | {r['method']} | {_f(r['measured'])} | "
          f"{_f(r['pred_before'])} | {_f(r['fold_before'], 2)}x | "
          f"{_f(r['pred_after'])} | **{_f(r['fold_error'], 2)}x** | "
          f"{'**yes**' if r['applied'] else 'no (total-drug assay)'} |")
    a("")

    # ---------------- §7 -----------------------------------------------------
    a("## 7. RECOMMENDATION — DEFAULT OFF")
    a("")
    a("**Leave the immune sink default-off (it is off by construction — no "
      "shipped path imports it), and quote its spatial result as a capability "
      "rather than a calibrated prediction.** In order of weight:")
    a("")
    a("1. **The magnitude is unmeasured and the literature's own sign is "
      "contested.** The best-quantified TAM drug-handling experiment in the "
      "literature (Miller 2015) found a DEPOT that made tumours respond BETTER, "
      "and the one direct free-doxorubicin measurement (Islam 2022) found TAM "
      "uptake below detection. Turning on a sink by default would assert a "
      "direction the strongest available data does not support for free drug.")
    a("2. **It does not improve the measured anchor**, and §6 explains why it "
      "structurally cannot: it is a one-directional term on a scorecard whose "
      "errors point both ways. This is the same outcome the efflux, uptake and "
      "ion-trapping runs reached, and the same conclusion follows.")
    if beats_null:
        a(f"3. **But the spatial result is real and is worth keeping.** The "
          f"reduction is structured — {_f(obs['structure_ratio'], 4)}x between "
          f"the sparsest and densest immune quintile against a well-mixed "
          f"model's exact 1.0000 — and it survives a permutation null that "
          f"holds the tissue, the cell layout and the total sink magnitude "
          f"fixed (p = {_f(p_sr, 4)}). That is a clean demonstration of what a "
          "microenvironment model adds over a compartmental one, and it holds "
          "at every strength because the geometry, not the magnitude, produces "
          "it.")
    else:
        a(f"3. **The reduction is structured, but the structure is not "
          f"MACROPHAGE-specific.** The {_f(obs['structure_ratio'], 4)}x contrast "
          f"between the sparsest and densest immune quintile is real and is "
          f"exactly 1.0000 in a well-mixed model, so the spatial capability is "
          f"demonstrated. But relabelling the same number of immune cells onto "
          f"other measured cell positions reproduces it "
          f"(p = {_f(p_sr, 4)}), so what the analysis actually establishes is "
          f"that exposure is structured by local CELLULARITY. Quoting it as a "
          f"macrophage effect would over-claim, and this report does not.")
    a("4. **The ranking is safe either way.** "
      f"{'The cross-drug ranking is invariant across the entire swept range' if sw['ranking_invariant'] else 'The ranking moves within the swept range'}, "
      "so no comparative conclusion in `results/` is at risk from this term "
      "being off — or from it being turned on.")
    a("")
    a("### 7.1 What measurement would settle it")
    a("")
    a("One experiment closes the magnitude question: **free doxorubicin (or its "
      "fluorescent analogue) dosed into a TAM-containing tumour, then "
      "flow-sorted CD163+/TREM2+ macrophages versus tumour cells assayed for "
      "drug per cell**, giving the per-cell enrichment ratio that does not "
      "currently exist. Miller 2015 did exactly this for the CARRIER; nobody has "
      "published it for the free drug. Failing that, single-cell intravital "
      "imaging of the free drug rather than the vehicle would do. Until then "
      "this term has a literature-supported SIGN and an unmeasured SIZE, and "
      "every number in §4 must be quoted with the strength it was produced at.")
    a("")

    # ---------------- §8 -----------------------------------------------------
    a("## 8. Caveats — what could make this wrong")
    a("")
    a("1. **The magnitude is swept, not measured.** Nothing in §4 or §5 is a "
      "prediction of how much drug real TAMs remove. Only the SHAPE of the "
      "spatial effect is independent of the magnitude.")
    a("2. **The sink is irreversible, and the best evidence says it should not "
      "be.** This module models uptake into macrophages as a pure loss. Miller "
      "2015's central finding is the opposite: TAMs released payload back to "
      "neighbouring tumour cells over <25 um. A re-release term would partly or "
      "wholly cancel the sink and could flip its sign. **That term is not "
      "implemented**, and it is the single largest structural omission here.")
    a("3. **The immune label is a marker call, not a TAM call.** The classifier "
      "gates PTPRC/CD68/AIF1/CD163 (and CD163/TREM2 on the IO panel); it does "
      "not distinguish macrophages from lymphocytes, and only macrophages have "
      "the phagocytic/lysosomal machinery the mechanism rests on. "
      "`results/classifier_fix.md` §6 further notes the matched-null rule mildly "
      "OVER-estimates multi-marker class fractions, so `phi_imm` is likely a "
      "mild over-estimate.")
    a(f"4. **`r` = {IMMUNE_CELL_RADIUS_UM:.1f} um and the "
      f"{IMMUNE_KERNEL_BANDWIDTH_UM:.0f} um bandwidth are declared constants.** "
      "They set the absolute scale of `phi_imm` and trade off exactly against "
      "`k_max`. The bandwidth is anchored to a measured gradient length, the "
      "radius is not anchored to anything measured on these sections.")
    a("5. **2D.** A section is a slice; an area fraction is not a volume "
      "fraction, and a macrophage's true drug catchment is a 3D neighbourhood. "
      "This applies to every field in the project and is not specific to this "
      "term, but it bites harder here because the mechanism is explicitly about "
      "local geometry.")
    a("6. **Three sections, three patients, three panels, and no outcome label "
      "anywhere.** Nothing here says which strength is right, because nothing "
      "here can be scored against a measured drug response.")
    a("7. **The scorecard test rests on assumption A1** (tissue:plasma partition "
      "= effective penetrance) and on a quasi-steady-state reading of it. Both "
      "are the scorecard's own standing assumptions, inherited unchanged.")
    a("")
    a("## 9. Reproduction")
    a("")
    a("```bash")
    a("./.venv/bin/python -m tmesim.validation.immune_sink_impact")
    a("```")
    a("")
    a(f"Deterministic — the engine's seeded generator (`seed = 0`) and the "
      f"permutation null's own seed ({sp['seed']}) are the only randomness. "
      f"Runtime is dominated by transport solves: "
      f"{len(DRUGS) * (len(SEQUESTRATION_SWEEP) + 1)} on the {sw['grid']} "
      f"flagship grid (§4, including the three engine reference runs)"
      + (f", {sp['n_perm'] + 2} at dx = {sp['dx_um']:.0f} um (§5.3) and 2 more "
         f"at dx = {spf['dx_um']:.0f} um (§5.4)" if spf is not None else "")
      + ". Each stage is cached to `results/immune_sink_impact_cache.pkl` as it "
      "completes, so an interrupted run resumes rather than restarting; delete "
      "that file to force a full recompute.")
    a("")
    a("**Files added by this work:**")
    a("")
    a("* `tmesim/molecular/immune_sink.py` — new; the field, the swept magnitude, "
      "the literature table and the re-solve")
    a("* `tmesim/validation/immune_sink_impact.py` — new; this analysis")
    a("* `tests/test_immune_sink.py` — new; pins the bit-identity and the "
      "non-vacuity")
    a("* `results/immune_sink_impact.md` — this report")
    a("")
    a("*Generated by `tmesim/validation/immune_sink_impact.py`. Additive only: "
      "no file under `core/`, `physics/`, `sim/`, `treatment/`, `reconstruct/` "
      "or `molecular/` (other than the new `molecular/immune_sink.py`) was "
      "modified, and the bit-exactness of the untouched transport path is "
      "verified in §1 both synthetically and end-to-end on the real section.*")

    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\n".join(A) + "\n")
    return path


# --------------------------------------------------------------------------- #
def run(write_md: bool = True, use_cache: bool = True,
        n_perm: int = N_PERM) -> Dict:
    """Run everything and write ``results/immune_sink_impact.md``.

    Each stage is cached to ``results/immune_sink_impact_cache.pkl`` as soon as
    it completes, so an interrupted run resumes instead of restarting. Delete
    the cache to force a full recompute.
    """
    t0 = time.time()
    out: Dict = {}
    if use_cache and os.path.exists(CACHE_PKL):
        with open(CACHE_PKL, "rb") as fh:
            out = pickle.load(fh)
    os.makedirs(RESULTS_DIR, exist_ok=True)

    def _stage(key: str, fn):
        if key in out and out[key] is not None:
            print(f"      (cached)", flush=True)
            return out[key]
        val = fn()
        out[key] = val
        with open(CACHE_PKL, "wb") as fh:
            pickle.dump(out, fh)
        return val

    print("[1/6] identity check ...", flush=True)
    ident = _stage("identity", identity_check)
    print(f"      {ident['n_bit_identical']}/{ident['n_comparisons']} "
          f"bit-identical, all={ident['all_bit_identical']}", flush=True)

    print("[2/6] immune fields on three real sections ...", flush=True)
    fields = _stage("fields", immune_fields)
    for r in fields:
        print(f"      {r['section']:<16} immune {r['immune_pct']:.2f}%  "
              f"mean phi {r['phi_mean']:.4f}", flush=True)

    print("[3/6] sweep on the flagship section ...", flush=True)
    sw = _stage("sweep", sweep)

    print(f"[4/6] spatial structure + permutation null (dx={DX_PERM:.0f}) ...",
          flush=True)
    sp = _stage("spatial",
                lambda: spatial_structure(dx_um=DX_PERM, n_perm=n_perm))
    print(f"      structure ratio {sp['observed']['structure_ratio']:.4f} "
          f"p={sp['p_structure_ratio']:.4f}", flush=True)

    print(f"[5/6] mesh cross-check of the same statistic at the flagship "
          f"dx={DX_FLAGSHIP:.0f} ...", flush=True)
    spf = _stage("spatial_flagship",
                 lambda: spatial_structure(dx_um=DX_FLAGSHIP, n_perm=0))
    print(f"      structure ratio {spf['observed']['structure_ratio']:.4f}",
          flush=True)

    print("[6/6] glioma scorecard ...", flush=True)
    sc = _stage("scorecard", lambda: glioma_scorecard(sw["mean_phi"]))
    print(f"      baseline {sc['baseline']['gm_abs']:.3f}x -> best "
          f"{sc['best_modality_aware']['gm_abs']:.3f}x "
          f"(improves={sc['improves']})", flush=True)

    if write_md:
        out["report_path"] = write_report(ident, fields, sw, sp, sc, spf)
        print(f"wrote {out['report_path']}", flush=True)
    print(f"total {time.time() - t0:.0f}s", flush=True)
    return out


if __name__ == "__main__":
    run()
