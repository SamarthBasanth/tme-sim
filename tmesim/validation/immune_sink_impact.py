"""What immune (TAM) drug sequestration is worth — on the REAL sections.

Runs the whole assessment of ``molecular/immune_sink.py`` and writes
``results/immune_sink_impact.md``:

  §1  the OFF path is bit-identical to the shipped solver, verified numerically
  §2  the measured immune field on all three real sections
  §3  the literature anchor converted through each section's own field
  §4  the sweep: mean tumour AUC, coverage at 1x IC50, predicted kill, RANKING
  §5  THE SPATIAL QUESTION — is the reduction structured, and by HOW MUCH once
      the nulls are honest: the measured well-mixed floors, the depth-into-tissue
      confound, the CV of the sink field, a clumpiness-matched permutation and a
      block shuffle as the primary comparators (the original label permutation
      is kept but DEMOTED to a floor check), and a hard refusal to compute any
      of it on spot data, where no immune cell was ever localised
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
from ..reconstruct.real_visium import _CELLS_PER_SPOT as _VISIUM_CELLS_PER_SPOT
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
# §5 (corrected) — THE HONEST NULLS
#
# The label permutation above is kept and still computed, but it is DEMOTED to a
# floor check: it destroys the immune compartment's clumping entirely, and
# clumping is the very quantity the structure ratio measures, so it produces a
# sink field far less variable than the observed one and loses a variance-driven
# contest almost regardless of whether any spatial biology exists. The CV of the
# sink field is emitted for the observed data and for every null so that this is
# visible rather than asserted.
#
# The two PRIMARY comparators keep every immune cell on a real MEASURED cell
# position, and differ from the label permutation only in that they are allowed
# to be as clumped as the data:
#
#   * clumpiness-matched permutation — the immune label moves in clumps of
#     CLUMP_SIZE mutually-nearest measured cells rather than one cell at a time;
#   * BLOCK_UM block shuffle — the per-block immune COUNTS are permuted across
#     blocks, so the observed block-scale count heterogeneity is inherited and
#     only its placement is randomised.
#
# The jitter reshuffle is a NEGATIVE CONTROL, not a null: it is the operation
# ``reconstruct/real_visium.py::_place_cells`` performs by construction, applied
# to the Xenium section. It is the control that entirely reproduces the statistic
# on spot data, so its COLLAPSE here is what separates the two platforms.
# --------------------------------------------------------------------------- #

#: Cells per clump in the clumpiness-matched permutation.
CLUMP_SIZE = 5

#: Block edge for the block shuffle (um).
BLOCK_UM = 400.0

#: The Visium synthesis, read off ``reconstruct/real_visium.py``: 12 pseudo-cells
#: per ~100 um spot, all at the spot centroid with uniform-in-disc jitter of at
#: most ``min(27.5, 0.45 * dx)`` um. The jitter reshuffle control applies exactly
#: that operation to the measured Xenium positions.
SPOT_PITCH_UM = 100.0
SPOT_JITTER_UM = 27.5

#: Clump sizes for the CV-vs-structure-ratio sweep (one draw each).
CLUMP_SWEEP_SIZES: Tuple[int, ...] = (1, 5, 20, 50)

#: Random shifts used for the torus-null FEASIBILITY diagnostic. The torus null
#: is never scored — the diagnostic exists to show why.
N_TORUS_PROBE = 20

#: The null families, in report order. ``jitter_reshuffle`` is a control.
NULL_FAMILIES: Tuple[str, ...] = (
    "clumpiness_matched", "block_shuffle", "label_permutation",
    "jitter_reshuffle", "spot_emulation",
)

_PRIMARY_NULLS: Tuple[str, ...] = ("clumpiness_matched", "block_shuffle")

NULL_SEEDS: Dict[str, int] = {
    "label_permutation": PERM_SEED,
    "clumpiness_matched": PERM_SEED + 1,
    "block_shuffle": PERM_SEED + 2,
    "jitter_reshuffle": PERM_SEED + 3,
    "spot_emulation": PERM_SEED + 4,
}

NULL_LABELS: Dict[str, str] = {
    "clumpiness_matched": "clumpiness-matched permutation",
    "block_shuffle": f"{BLOCK_UM:.0f} um block shuffle",
    "label_permutation": "label permutation (as originally published)",
    "jitter_reshuffle": "jitter reshuffle (re-draw the Visium jitter only)",
    "spot_emulation": "spot emulation (the whole Visium synthesis)",
}


class SpotDataError(RuntimeError):
    """Raised when a single-cell-only statistic is asked for on spot data."""


#: Reconstruction sources that DO measure individual cell positions.
_SINGLE_CELL_SOURCES = frozenset({"xenium", "real_xenium", "cosmx", "merfish"})
#: Reconstruction sources that synthesise pseudo-cells from spot deconvolution.
_SPOT_SOURCES = frozenset({"visium", "real_visium", "spot"})


def reconstruction_provenance(micro: Microenvironment) -> Dict:
    """Where the cell POSITIONS in ``micro`` came from, from its own metadata.

    This reads ``metadata['real_pipeline']['source']`` (set by
    ``reconstruct/real_pipeline.py`` for every reconstruction it returns) and the
    spot-level fingerprints the Visium reconstructor leaves behind
    (``metadata['source'] == 'real_visium'``, ``metadata['n_spots']``). It never
    looks at the sample PATH, so a section moved or renamed on disk cannot fool
    it.
    """
    md = dict(micro.metadata or {})
    rp = md.get("real_pipeline") or {}
    source = str(rp.get("source") or md.get("source") or "unknown").strip().lower()
    n_spots = md.get("n_spots")
    spot_based = bool(source in _SPOT_SOURCES or n_spots is not None)
    single_cell = bool(source in _SINGLE_CELL_SOURCES and not spot_based)
    if spot_based:
        reason = (
            f"reconstruction source {source!r}"
            + (f" with {int(n_spots)} spots" if n_spots is not None else "")
            + " — cell positions are SYNTHESISED: "
            f"reconstruct/real_visium.py::_place_cells puts a fixed number of "
            f"pseudo-cells at each spot centroid with uniform-in-disc jitter of "
            f"at most {SPOT_JITTER_UM:g} um, apportioned by deconvolved "
            f"fractions. No individual immune cell was ever localised."
        )
    elif single_cell:
        reason = f"reconstruction source {source!r} — measured single-cell centroids"
    else:
        reason = (
            f"reconstruction source {source!r} is not a recognised single-cell "
            f"platform; positional provenance cannot be established"
        )
    return {
        "source": source, "n_spots": (int(n_spots) if n_spots is not None else None),
        "spot_based": spot_based, "single_cell": single_cell, "reason": reason,
    }


def require_measured_cell_positions(
        micro: Microenvironment,
        what: str = "the immune-density structure ratio") -> Dict:
    """Refuse to compute ``what`` unless the immune POSITIONS were measured.

    The structure ratio asks whether tumour cells near MEASURED macrophages lose
    more drug. On a spot platform there are no measured macrophage positions to
    be near: everything below the ~100 um spot pitch is a random draw, and the
    sink field is a 25 um-bandwidth KDE of those draws, so every plausible null
    is comparing a field to itself. The statistic is not weak there, it is
    undefined, so it is refused rather than degraded.
    """
    prov = reconstruction_provenance(micro)
    if not prov["single_cell"]:
        raise SpotDataError(
            f"{what} is INVALID on this reconstruction and will not be "
            f"computed: {prov['reason']}. See results/immune_sink_impact.md "
            f"§5.6. Pass allow_invalid=True only to DEMONSTRATE the invalidity; "
            f"the result is then stamped INVALID and must never be quoted."
        )
    return prov


def _field_cv(field: np.ndarray, tissue: np.ndarray) -> float:
    """Coefficient of variation of a field over the tissue mask."""
    v = np.asarray(field, dtype=np.float64)[np.asarray(tissue, dtype=bool)]
    m = float(v.mean())
    return float(v.std() / m) if m > 0 else float("nan")


def _spatial_context(strength: float = SPATIAL_STRENGTH,
                     drug_name: str = "doxorubicin",
                     dx_um: float = DX_PERM,
                     sample: str = FLAGSHIP_SECTION,
                     allow_invalid: bool = False,
                     verbose: bool = True,
                     micro: Optional[Microenvironment] = None,
                     drug: Optional[DrugParams] = None,
                     pde: Optional[PDEConfig] = None) -> Dict:
    """Everything every null shares: the section, the prologue, and 2 solves.

    Built ONCE and reused by every null family and every diagnostic, so the only
    thing that differs between the observed value and any null is `k_sink`.

    ``micro`` / ``drug`` / ``pde`` are injection points for the tests, which run
    the identical code on a small synthetic section; leaving them at ``None``
    gives the flagship path exactly as §5 runs it.
    """
    if micro is None:
        micro = _reconstruct(sample, dx_um)
    prov = reconstruction_provenance(micro)
    if not allow_invalid:
        require_measured_cell_positions(micro)
    if drug is None:
        alpha = calibrate_alpha()
        drug, _mol = compile_real_drug(drug_name, alpha)
    dp = drug
    cfg = SimConfig(drug=dp, pde=(pde if pde is not None
                                  else PDEConfig(**_PDE)), seed=0)
    work, labeled, velocity = _prologue(micro, cfg)
    source = _source_for(work, dp, cfg)
    vel = velocity if cfg.pde.convection else None

    cells = micro.cells
    ct = np.asarray(cells.cell_type)
    tumor = ct == int(CellType.TUMOR)
    immune = ct == int(CellType.IMMUNE)
    x = np.asarray(cells.x, dtype=np.float64)
    y = np.asarray(cells.y, dtype=np.float64)
    col, row = micro.domain.to_index(x, y)
    tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)

    t0 = time.time()
    base = transport_with_immune_sink(work, dp, source, cfg.pde, velocity=vel,
                                      sink_rate=None)
    auc_base = np.asarray(base.auc)[row, col][tumor]
    phi = immune_area_fraction(micro)
    k_max = k_max_from_strength(dp, strength)
    sink = immune_sink_rate(micro, k_max, phi=phi)
    obs_tr = transport_with_immune_sink(work, dp, source, cfg.pde, velocity=vel,
                                        sink_rate=sink)
    auc_obs = np.asarray(obs_tr.auc)[row, col][tumor]
    if verbose:
        print(f"    context: 2 solves in {time.time() - t0:.0f}s", flush=True)

    # depth into tissue — the confound §5.1a measures, computed once
    from scipy.ndimage import distance_transform_edt
    edge_um = distance_transform_edt(tissue) * float(micro.domain.dx)

    observed = _structure_stats(auc_base, auc_obs, phi[row, col][tumor])
    observed["sink_cv"] = _field_cv(sink, tissue)
    observed["phi_mean_tissue"] = float(phi[tissue].mean())

    return {
        "micro": micro, "work": work, "cfg": cfg, "drug": dp, "source": source,
        "vel": vel, "tissue": tissue, "row": row, "col": col, "x": x, "y": y,
        "tumor": tumor, "immune": immune, "n_cells": int(ct.size),
        "n_immune": int(immune.sum()), "n_tumor": int(tumor.sum()),
        "phi": phi, "sink": sink, "k_max": float(k_max),
        "auc_base": auc_base, "auc_obs": auc_obs,
        "edge_um": edge_um, "edge_at_tumor": edge_um[row, col][tumor],
        "phi_at_tumor": phi[row, col][tumor],
        "observed": observed, "provenance": prov,
        "strength": float(strength), "drug_name": drug_name,
        "dx_um": float(dx_um), "sample": sample,
    }


def _stats_for_sink(ctx: Dict, sink_field: np.ndarray,
                    bin_by: Optional[np.ndarray] = None) -> Dict:
    """One transport re-solve, and the structure statistics that come out of it.

    ``bin_by`` is the per-tumour-cell variable the quintiles are formed on and
    defaults to the sink field's OWN local immune density — i.e. every null is
    scored on exactly the statistic the observed value was scored on.
    """
    tr = transport_with_immune_sink(ctx["work"], ctx["drug"], ctx["source"],
                                    ctx["cfg"].pde, velocity=ctx["vel"],
                                    sink_rate=sink_field)
    auc = np.asarray(tr.auc)[ctx["row"], ctx["col"]][ctx["tumor"]]
    if bin_by is None:
        bin_by = sink_field[ctx["row"], ctx["col"]][ctx["tumor"]]
    st = _structure_stats(ctx["auc_base"], auc, np.asarray(bin_by))
    st["sink_cv"] = _field_cv(sink_field, ctx["tissue"])
    st["sink_mean_tissue"] = float(sink_field[ctx["tissue"]].mean())
    return st


def _stats_for_positions(ctx: Dict, px: np.ndarray, py: np.ndarray) -> Dict:
    """The null's own immune positions -> its own phi -> its own sink -> stats."""
    phi_p = immune_area_fraction(ctx["micro"], positions=(px, py))
    sink_p = immune_sink_rate(ctx["micro"], ctx["k_max"], phi=phi_p)
    st = _stats_for_sink(ctx, sink_p, bin_by=phi_p[ctx["row"], ctx["col"]][ctx["tumor"]])
    st["n_placed"] = int(np.asarray(px).size)
    st["phi_mean_tissue"] = float(phi_p[ctx["tissue"]].mean())
    st["corr_phi_edge"] = float(np.corrcoef(
        phi_p[ctx["row"], ctx["col"]][ctx["tumor"]], ctx["edge_at_tumor"])[0, 1])
    return st


# --- the position generators ------------------------------------------------ #
def _draw_label_permutation(ctx: Dict, rng) -> Tuple[np.ndarray, np.ndarray]:
    """The ORIGINAL null: the immune label moves one cell at a time."""
    pick = rng.choice(ctx["n_cells"], size=ctx["n_immune"], replace=False)
    return ctx["x"][pick], ctx["y"][pick]


def _draw_clumpiness_matched(ctx: Dict, rng,
                             clump_size: int = CLUMP_SIZE
                             ) -> Tuple[np.ndarray, np.ndarray]:
    """The immune label moves in CLUMPS of mutually-nearest measured cells.

    A seed cell is drawn at random and its ``clump_size`` nearest MEASURED
    neighbours (of any type) take the immune label together, until the observed
    immune count is reached. Every immune cell therefore still sits on a real
    measured cell position and the total count is exactly preserved — the only
    thing that changes relative to the label permutation is that the null is
    allowed to be as clumped as the data.
    """
    tree = ctx["_tree_all"]
    pts = ctx["_pts_all"]
    n_cells, n_immune = ctx["n_cells"], ctx["n_immune"]
    taken = np.zeros(n_cells, dtype=bool)
    picked: List[int] = []
    k = max(int(clump_size), 1)
    while len(picked) < n_immune:
        need = n_immune - len(picked)
        n_seed = max(1, int(np.ceil(need / k)))
        seeds = rng.choice(n_cells, size=min(n_seed, n_cells), replace=False)
        _d, nbr = tree.query(pts[seeds], k=k)
        nbr = np.atleast_2d(nbr)
        for j in nbr.reshape(-1):
            if not taken[j]:
                taken[j] = True
                picked.append(int(j))
                if len(picked) >= n_immune:
                    break
    idx = np.asarray(picked[:n_immune], dtype=int)
    return ctx["x"][idx], ctx["y"][idx]


def _draw_block_shuffle(ctx: Dict, rng, block_um: float = BLOCK_UM
                        ) -> Tuple[np.ndarray, np.ndarray]:
    """Permute the per-block immune COUNTS across ``block_um`` blocks.

    The observed section is tiled into ``block_um`` squares. Each block's
    observed immune count is reassigned to a randomly chosen OTHER block, which
    then draws that many of its own measured cells to carry the immune label.
    The block-scale count heterogeneity of the real section is therefore
    inherited exactly (only relocated), every immune cell stays on a measured
    cell position, and the total count is preserved up to block capacity — any
    shortfall is topped up among blocks with spare capacity, and the realised
    count is reported.
    """
    x, y = ctx["x"], ctx["y"]
    bx = np.floor(x / float(block_um)).astype(np.int64)
    by = np.floor(y / float(block_um)).astype(np.int64)
    key = bx * 1000003 + by
    uniq, inv = np.unique(key, return_inverse=True)
    n_blocks = int(uniq.size)
    members = [np.flatnonzero(inv == b) for b in range(n_blocks)]
    imm = ctx["immune"]
    counts = np.array([int(imm[m].sum()) for m in members], dtype=int)
    cap = np.array([int(m.size) for m in members], dtype=int)

    order = rng.permutation(n_blocks)
    want = counts[order]                       # block b receives block order[b]'s count
    take = np.minimum(want, cap)
    picked: List[np.ndarray] = []
    for b in range(n_blocks):
        if take[b] <= 0:
            continue
        picked.append(rng.choice(members[b], size=int(take[b]), replace=False))
    chosen = np.concatenate(picked) if picked else np.zeros(0, dtype=int)

    deficit = int(ctx["n_immune"] - chosen.size)
    if deficit > 0:                            # top up where there is room
        free = np.ones(ctx["n_cells"], dtype=bool)
        free[chosen] = False
        spare = np.flatnonzero(free)
        extra = rng.choice(spare, size=min(deficit, spare.size), replace=False)
        chosen = np.concatenate([chosen, extra])
    return x[chosen], y[chosen]


def _draw_jitter_reshuffle(ctx: Dict, rng, jitter_um: float = SPOT_JITTER_UM
                           ) -> Tuple[np.ndarray, np.ndarray]:
    """Re-draw the synthetic JITTER, and nothing else.

    Every measured immune cell is displaced by a uniform-in-disc offset of at
    most ``jitter_um`` — the exact amplitude ``real_visium.py::_place_cells``
    applies (``jitter_um = min(27.5, 0.45 dx)``). Nothing about which
    neighbourhood a macrophage is in changes; only positional detail AT THE
    VISIUM JITTER SCALE is re-randomised.

    On a spot reconstruction this control is a no-op by construction — the
    positions WERE a centroid plus this jitter — so it reproduces the observed
    statistic there and leaves nothing. That is why its behaviour here is
    decisive: if the Xenium number were the same artefact, this would reproduce
    it here too.
    """
    x = ctx["x"][ctx["immune"]]
    y = ctx["y"][ctx["immune"]]
    ang = rng.uniform(0.0, 2.0 * np.pi, size=x.size)
    rad = float(jitter_um) * np.sqrt(rng.uniform(0.0, 1.0, size=x.size))
    w, h = ctx["micro"].domain.extent_um
    return (np.clip(x + rad * np.cos(ang), 0.0, w),
            np.clip(y + rad * np.sin(ang), 0.0, h))


def _draw_spot_emulation(ctx: Dict, rng, pitch_um: float = SPOT_PITCH_UM,
                         jitter_um: float = SPOT_JITTER_UM
                         ) -> Tuple[np.ndarray, np.ndarray]:
    """The whole Visium synthesis, applied to the measured Xenium immune cells.

    Every measured immune cell is snapped to the centroid of its ``pitch_um``
    lattice cell and re-placed with uniform-in-disc jitter of at most
    ``jitter_um`` — literally what ``real_visium.py::_place_cells`` does. The
    per-spot immune COUNT is preserved exactly; everything finer is destroyed and
    redrawn. This answers "what would this section's statistic be if it had been
    measured at spot resolution", and it is reported as a CONTROL rather than a
    null because snapping to a lattice CONCENTRATES the point pattern and
    therefore biases the sink-field CV upward, in the same mechanical way the
    torus shift does.
    """
    x = ctx["x"][ctx["immune"]]
    y = ctx["y"][ctx["immune"]]
    cx = (np.floor(x / float(pitch_um)) + 0.5) * float(pitch_um)
    cy = (np.floor(y / float(pitch_um)) + 0.5) * float(pitch_um)
    ang = rng.uniform(0.0, 2.0 * np.pi, size=x.size)
    rad = float(jitter_um) * np.sqrt(rng.uniform(0.0, 1.0, size=x.size))
    w, h = ctx["micro"].domain.extent_um
    return (np.clip(cx + rad * np.cos(ang), 0.0, w),
            np.clip(cy + rad * np.sin(ang), 0.0, h))


def _draw_uniform_in_tissue(ctx: Dict, rng) -> Tuple[np.ndarray, np.ndarray]:
    """Scatter the same number of immune cells uniformly INSIDE the tissue mask.

    Contains no measured immune structure whatsoever, which is what makes it the
    right control for the monotonicity question (§5.4).
    """
    tissue = ctx["tissue"]
    dx = float(ctx["micro"].domain.dx)
    rr, cc = np.nonzero(tissue)
    pick = rng.choice(rr.size, size=ctx["n_immune"], replace=True)
    px = (cc[pick] + rng.uniform(0.0, 1.0, size=pick.size)) * dx
    py = (rr[pick] + rng.uniform(0.0, 1.0, size=pick.size)) * dx
    return px, py


_DRAWERS = {
    "label_permutation": _draw_label_permutation,
    "clumpiness_matched": _draw_clumpiness_matched,
    "block_shuffle": _draw_block_shuffle,
    "jitter_reshuffle": _draw_jitter_reshuffle,
    "spot_emulation": _draw_spot_emulation,
}


def _ensure_tree(ctx: Dict) -> None:
    if "_tree_all" not in ctx:
        from scipy.spatial import cKDTree
        pts = np.column_stack([ctx["x"], ctx["y"]])
        ctx["_pts_all"] = pts
        ctx["_tree_all"] = cKDTree(pts)


def null_family(ctx: Dict, family: str, n_draws: int = N_PERM,
                seed: Optional[int] = None, verbose: bool = True) -> Dict:
    """Run one null family: ``n_draws`` full transport re-solves.

    Emits, per draw and in the summary, the CV of the sink field — the diagnostic
    that shows WHICH nulls are being asked a fair question and which are not.
    """
    if family not in _DRAWERS:
        raise ValueError(f"unknown null family {family!r}; "
                         f"expected one of {sorted(_DRAWERS)}")
    _ensure_tree(ctx)
    rng = np.random.default_rng(NULL_SEEDS[family] if seed is None else seed)
    draw = _DRAWERS[family]
    obs = ctx["observed"]
    rows: List[Dict] = []
    for i in range(int(n_draws)):
        t0 = time.time()
        px, py = draw(ctx, rng)
        st = _stats_for_positions(ctx, px, py)
        st["draw"] = i
        rows.append(st)
        if verbose:
            print(f"    [{family}] {i + 1}/{n_draws} {time.time() - t0:.0f}s "
                  f"SR={st['structure_ratio']:.4f} CV={st['sink_cv']:.3f} "
                  f"n={st['n_placed']}", flush=True)

    def _arr(k: str) -> np.ndarray:
        return np.array([r[k] for r in rows], dtype=np.float64)

    sr, sp_, cv = _arr("structure_ratio"), _arr("a_p90_over_p10"), _arr("sink_cv")
    am = _arr("a_mean")
    quint = np.array([[b["attenuation_mean"] for b in r["bins"]] for r in rows],
                     dtype=np.float64)
    qmean = quint.mean(axis=0) if quint.size else np.full(N_BINS, np.nan)
    return {
        "family": family, "label": NULL_LABELS[family], "n_draws": int(n_draws),
        "seed": int(NULL_SEEDS[family] if seed is None else seed),
        "is_primary": family in _PRIMARY_NULLS,
        "sr_mean": float(sr.mean()), "sr_min": float(sr.min()),
        "sr_max": float(sr.max()), "sr_std": float(sr.std()),
        "spread_mean": float(sp_.mean()), "spread_min": float(sp_.min()),
        "spread_max": float(sp_.max()),
        "cv_mean": float(cv.mean()), "cv_min": float(cv.min()),
        "cv_max": float(cv.max()),
        "a_mean_mean": float(am.mean()),
        "n_placed_mean": float(np.mean([r["n_placed"] for r in rows])),
        "corr_phi_edge_mean": float(np.mean([r["corr_phi_edge"] for r in rows])),
        "quintile_attenuation_mean": [float(v) for v in qmean],
        "monotone": bool(np.all(np.diff(qmean) < 0.0)),
        # one-sided permutation p-values; at n draws the FLOOR is 1/(1+n)
        "p_structure_ratio": float((1.0 + np.sum(sr >= obs["structure_ratio"]))
                                   / (1.0 + sr.size)),
        "p_spread": float((1.0 + np.sum(sp_ >= obs["a_p90_over_p10"]))
                          / (1.0 + sp_.size)),
        "p_floor": float(1.0 / (1.0 + sr.size)),
        "excess_over_null": float(obs["structure_ratio"] / sr.mean())
        if sr.mean() > 0 else float("nan"),
        "draws": rows,
    }


def spatial_diagnostics(ctx: Dict, verbose: bool = True) -> Dict:
    """The floors, the depth confound, the CV sweep and the torus feasibility.

    Costs 2 + ``len(CLUMP_SWEEP_SIZES)`` transport re-solves. Nothing here is a
    scored null; every entry exists to stop a statistic being quoted against an
    ASSUMED floor of 1.0.
    """
    _ensure_tree(ctx)
    tissue, phi = ctx["tissue"], ctx["phi"]
    obs = ctx["observed"]

    # --- floor 1: a genuinely WELL-MIXED sink of the same tissue-mean magnitude
    t0 = time.time()
    uniform = np.zeros_like(phi)
    uniform[tissue] = float(ctx["sink"][tissue].mean())
    uni = _stats_for_sink(ctx, uniform, bin_by=ctx["phi_at_tumor"])
    if verbose:
        print(f"    [floor] uniform sink {time.time() - t0:.0f}s "
              f"SR={uni['structure_ratio']:.4f} "
              f"p90/p10={uni['a_p90_over_p10']:.4f}", flush=True)

    # --- the depth-into-tissue confound, with NO immune input at all ----------
    from scipy.stats import spearmanr
    edge = ctx["edge_at_tumor"]
    edge_stats = _structure_stats(ctx["auc_base"], ctx["auc_obs"], edge)
    corr_pearson = float(np.corrcoef(ctx["phi_at_tumor"], edge)[0, 1])
    corr_spearman = float(spearmanr(ctx["phi_at_tumor"], edge).statistic)

    # within-edge-stratum immune structure ratio
    q = np.quantile(edge, np.linspace(0.0, 1.0, N_BINS + 1))
    q[0] -= 1e-12
    strat_idx = np.clip(np.searchsorted(q, edge, side="left") - 1, 0, N_BINS - 1)
    strata = []
    for s in range(N_BINS):
        m = strat_idx == s
        if m.sum() < 5 * N_BINS:
            strata.append(float("nan"))
            continue
        st = _structure_stats(ctx["auc_base"][m], ctx["auc_obs"][m],
                              ctx["phi_at_tumor"][m])
        strata.append(float(st["structure_ratio"]))

    # --- the monotonicity control: uniform-in-tissue scatter ------------------
    t0 = time.time()
    rng = np.random.default_rng(PERM_SEED + 11)
    ux, uy = _draw_uniform_in_tissue(ctx, rng)
    scatter = _stats_for_positions(ctx, ux, uy)
    if verbose:
        print(f"    [control] uniform-in-tissue scatter {time.time() - t0:.0f}s "
              f"SR={scatter['structure_ratio']:.4f} "
              f"CV={scatter['sink_cv']:.3f}", flush=True)

    # --- CV sweep: the structure ratio very largely TRACKS sink-field CV ------
    sweep_rows: List[Dict] = []
    for k in CLUMP_SWEEP_SIZES:
        t0 = time.time()
        r = np.random.default_rng(PERM_SEED + 100 + int(k))
        px, py = _draw_clumpiness_matched(ctx, r, clump_size=int(k))
        st = _stats_for_positions(ctx, px, py)
        sweep_rows.append({"clump_size": int(k), "sink_cv": st["sink_cv"],
                           "structure_ratio": st["structure_ratio"]})
        if verbose:
            print(f"    [cv-sweep] clump={k} {time.time() - t0:.0f}s "
                  f"CV={st['sink_cv']:.3f} SR={st['structure_ratio']:.4f}",
                  flush=True)

    # --- torus FEASIBILITY (never scored) ------------------------------------
    r = np.random.default_rng(PERM_SEED + 7)
    w, h = ctx["micro"].domain.extent_um
    xi = ctx["x"][ctx["immune"]]
    yi = ctx["y"][ctx["immune"]]
    inside, phi_ratio, cv_ratio = [], [], []
    for _ in range(N_TORUS_PROBE):
        sx, sy = r.uniform(0.0, w), r.uniform(0.0, h)
        tx, ty = np.mod(xi + sx, w), np.mod(yi + sy, h)
        c, rw = ctx["micro"].domain.to_index(tx, ty)
        keep = tissue[rw, c]
        inside.append(float(keep.mean()))
        phi_t = immune_area_fraction(ctx["micro"],
                                     positions=(tx[keep], ty[keep]))
        phi_ratio.append(float(phi_t[tissue].mean() / phi[tissue].mean()))
        cv_ratio.append(_field_cv(phi_t, tissue) / obs["sink_cv"])

    # --- the tissue's own density envelope (why depth carries a gradient) -----
    rim_um = max(IMMUNE_KERNEL_BANDWIDTH_UM, float(ctx["micro"].domain.dx))
    rim = ctx["edge_um"][tissue] <= rim_um
    from scipy.ndimage import gaussian_filter
    counts = ctx["micro"].domain.zeros()
    np.add.at(counts, (ctx["row"], ctx["col"]), 1.0)
    dens = gaussian_filter(counts, sigma=max(IMMUNE_KERNEL_BANDWIDTH_UM
                                             / float(ctx["micro"].domain.dx), 1e-6),
                           mode="constant")
    rim_mean = float(dens[tissue][rim].mean()) if rim.any() else float("nan")
    core_mean = float(dens[tissue][~rim].mean()) if (~rim).any() else float("nan")

    return {
        "uniform_sink_floor": uni,
        "edge_distance": {
            "stats": edge_stats,
            "corr_phi_edge_pearson": corr_pearson,
            "corr_phi_edge_spearman": corr_spearman,
            "within_stratum_sr": strata,
            "density_rim_mean": rim_mean, "density_core_mean": core_mean,
            "density_envelope_ratio": (core_mean / rim_mean
                                       if rim_mean > 0 else float("nan")),
            "rim_um": float(rim_um),
        },
        "uniform_in_tissue_scatter": scatter,
        "cv_sweep": sweep_rows,
        "torus": {
            "n_probe": N_TORUS_PROBE,
            "fraction_in_tissue_mean": float(np.mean(inside)),
            "fraction_in_tissue_min": float(np.min(inside)),
            "fraction_in_tissue_max": float(np.max(inside)),
            "phi_mean_ratio": float(np.mean(phi_ratio)),
            "cv_ratio": float(np.mean(cv_ratio)),
        },
    }


def honest_spatial_summary(ctx_meta: Dict, families: Dict[str, Dict],
                           diag: Dict) -> Dict:
    """Fold the families and the diagnostics into the numbers §5 leads with."""
    obs = ctx_meta["observed"]
    primary = [families[f] for f in _PRIMARY_NULLS if f in families]
    excess = [p["excess_over_null"] for p in primary]
    return {
        "observed_structure_ratio": float(obs["structure_ratio"]),
        "observed_spread": float(obs["a_p90_over_p10"]),
        "observed_sink_cv": float(obs["sink_cv"]),
        "honest_excess_min": float(min(excess)) if excess else float("nan"),
        "honest_excess_max": float(max(excess)) if excess else float("nan"),
        "primary_families": list(_PRIMARY_NULLS),
        "well_mixed_sr_floor": float(diag["uniform_sink_floor"]["structure_ratio"]),
        "well_mixed_spread_floor": float(diag["uniform_sink_floor"]["a_p90_over_p10"]),
        "spread_excess_over_floor": float(
            obs["a_p90_over_p10"] / diag["uniform_sink_floor"]["a_p90_over_p10"]),
        "edge_only_structure_ratio": float(
            diag["edge_distance"]["stats"]["structure_ratio"]),
        # the floor is 1/(1+draws) at the number of draws ACTUALLY run
        "p_floor": float(max(p["p_floor"] for p in primary))
        if primary else float("nan"),
        "p_primary_max": float(max(p["p_structure_ratio"] for p in primary))
        if primary else float("nan"),
        "above_all_primary_null_draws": bool(all(
            obs["structure_ratio"] > p["sr_max"] for p in primary)),
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
                 hc: Optional[Dict] = None,
                 families: Optional[Dict[str, Dict]] = None,
                 diag: Optional[Dict] = None,
                 hon: Optional[Dict] = None,
                 path: str = REPORT_MD) -> str:
    if not (hc and families and diag and hon):
        raise ValueError(
            "write_report requires the honest-null results (hc, families, "
            "diag, hon). Refusing to regenerate results/immune_sink_impact.md "
            "from the label-permutation null alone: that is exactly the claim "
            "the §5 correction removed, and a report that silently reverts to "
            "it is worse than no report. Run tmesim.validation."
            "immune_sink_impact.run(), which computes them."
        )
    fams: Dict[str, Dict] = families
    dg: Dict = diag

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
    beats_null = bool(hon["above_all_primary_null_draws"])
    _lab, _clump, _block = (fams["label_permutation"],
                            fams["clumpiness_matched"], fams["block_shuffle"])
    _floor = dg["uniform_sink_floor"]
    _edge = dg["edge_distance"]
    ex_lo, ex_hi = hon["honest_excess_min"], hon["honest_excess_max"]
    p_floor = hon["p_floor"]
    null_clause = (
        f"the observed value still sits above their entire ranges "
        f"(p = {_f(p_floor, 4)}, the {_lab['n_draws']}-draw floor)"
        if beats_null else
        f"the observed value does NOT clear their ranges, so at "
        f"{_lab['n_draws']} draws the measured macrophage positions are not "
        f"distinguishable from an equally-clumped rearrangement of them"
    )

    # ---- the CORRECTION NOTICE, emitted by the generator that computes it ----
    a("> ## CORRECTION NOTICE — the headline spatial number is smaller than this "
      "report originally claimed")
    a(">")
    a(f"> An audit of §5 found that the "
      f"{_f(obs['structure_ratio'], 4)}x structure ratio was being quoted "
      f"against a null that is too weak to be informative, and that a large part "
      f"of the spread is a transport gradient with depth into tissue rather than "
      f"immune biology. **The honest excess attributable to where the "
      f"macrophages actually sit is ~{_f(ex_lo, 2)}-{_f(ex_hi, 2)}x, not "
      f"{_f(obs['structure_ratio'], 1)}x.** "
      + ("The effect is still there and it still clears an honest null — but "
         "the excess over that null is "
         f"{_f(ex_lo - 1.0, 2)}-{_f(ex_hi - 1.0, 2)}, against the "
         f"{_f(obs['structure_ratio'] - 1.0, 1)} excess that "
         f"\"{_f(obs['structure_ratio'], 1)}x versus a well-mixed 1.0\" "
         f"advertised."
         if beats_null else
         "The effect does not clear an honest null at all."))
    a(">")
    a("> The nulls, the floors, the depth confound and the CV diagnostic below "
      "are all computed by "
      "`tmesim/validation/immune_sink_impact.py`, so re-running the command in "
      "§9 reproduces this correction rather than reverting it. What changed, in "
      "one place:")
    a(">")
    a("> | claim as originally published | corrected |")
    a("> |---|---|")
    a(f"> | structure ratio {_f(obs['structure_ratio'], 4)}x \"against a "
      f"well-mixed model's exactly 1.0\" | observed "
      f"{_f(obs['structure_ratio'], 4)}x is real and reproduces to the digit. "
      f"But the defensible comparator is a **spatially-aware** null, not a "
      f"well-mixed one: clumpiness-matched permutation "
      f"**{_f(_clump['sr_mean'], 4)}**, {BLOCK_UM:.0f} um block shuffle "
      f"**{_f(_block['sr_mean'], 4)}**. **Excess = "
      f"{_f(_clump['excess_over_null'], 3)}x and "
      f"{_f(_block['excess_over_null'], 3)}x.** |")
    a(f"> | the label-permutation null (mean {_f(_lab['sr_mean'], 4)}) is the "
      f"evidence | that null under-clumps the sink field it is being compared "
      f"against — CV **{_f(_lab['cv_mean'], 2)}** vs the observed "
      f"**{_f(hon['observed_sink_cv'], 3)}** — so it loses almost "
      f"automatically. It is **demoted** to a floor check (§5.3). |")
    a(f"> | the quintile gradient is an immune effect | binning the same tumour "
      f"cells by **distance to tissue edge alone**, with no reference to immune "
      f"density at all, already gives SR "
      f"**{_f(_edge['stats']['structure_ratio'], 4)}** (§5.1a). |")
    a("> | monotonicity across all five quintiles is corroboration | **every** "
      "control is monotone, including uniform-in-tissue scatter. Monotonicity "
      "is not evidence (§5.4). |")
    a(f"> | p90/p10 spread {_f(obs['a_p90_over_p10'], 4)}x against a well-mixed "
      f"1.0 | the well-mixed floor for that statistic is "
      f"**{_f(_floor['a_p90_over_p10'], 4)}**, not 1.0 (§5.1). |")
    a("> | (not previously stated) | **the statistic is INVALID on Visium/spot "
      "data**; the generator now REFUSES to compute it there (§5.6). |")
    a(">")
    a("> Nothing in §1-§4 or §6 changes. Those sections were not implicated by "
      "the audit.")
    a("")
    imm_lo = min(r["immune_pct"] for r in fields)
    imm_hi = max(r["immune_pct"] for r in fields)
    a("**Verdict up front.** The immune compartment was pharmacologically inert: "
      f"{imm_lo:.1f}-{imm_hi:.1f}% of the segmented cells across the three real "
      "sections occupied space and consumed nothing. It "
      "now has a mechanism — a first-order drug loss proportional to the local "
      "MEASURED immune-cell area fraction — and the honest findings are three. "
      "**(i) The literature supports a DIRECTION and not a MAGNITUDE**, and for "
      "free doxorubicin the one direct measurement is NEGATIVE, so the magnitude "
      "is swept and never fitted. **(ii) The reduction is spatially structured "
      "rather than uniform, but the structure attributable to immune geometry is "
      f"~{_f(ex_lo, 2)}-{_f(ex_hi, 2)}x, not the "
      f"{_f(obs['structure_ratio'], 1)}x this report previously led with.** The "
      f"raw quintile contrast is {_f(obs['structure_ratio'], 4)}x (densest "
      f"quintile keeps {_f(obs['bins'][-1]['attenuation_mean'], 4)} of its "
      f"exposure against {_f(obs['bins'][0]['attenuation_mean'], 4)} in the "
      f"sparsest), but most of that is reproduced by nulls that keep the "
      f"macrophages' clumpiness ({_f(_clump['sr_mean'], 4)} and "
      f"{_f(_block['sr_mean'], 4)}) and by a depth-into-tissue gradient with no "
      f"immune content at all ({_f(_edge['stats']['structure_ratio'], 4)}); "
      f"against honest spatially-aware nulls {null_clause}, which is what "
      "remains defensible. **(iii) It does not improve the measured "
      "human-glioma anchor**, so it ships DEFAULT-OFF.")
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
    a(f"| spatial structure ratio, raw (sparsest / densest immune quintile) | "
      f"{_f(obs['structure_ratio'], 4)}x |")
    _p_prim = hon["p_primary_max"]
    a(f"| **same ratio, excess over an honest spatially-aware null** | "
      f"**{_f(ex_lo, 2)}-{_f(ex_hi, 2)}x**, p = {_f(_p_prim, 4)}"
      + (f" (= the {_lab['n_draws']}-permutation floor)"
         if abs(_p_prim - p_floor) < 1e-12 else "") + " |")
    a(f"| same ratio from tissue-EDGE DISTANCE alone, no immune input | "
      f"{_f(_edge['stats']['structure_ratio'], 4)}x — see §5.1a |")
    a(f"| validity on Visium / spot transcriptomics | **INVALID — refused by "
      f"the generator (§5.6)** |")
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
    a("**Read the three drugs against each other — the asymmetry is the useful "
      "part, and it is not about the sink.** The AUC column is nearly identical "
      "across drugs, because the sink is applied with the same `S` to all three "
      "and the transport is close to linear in it. What differs is how that same "
      "fractional exposure loss TRANSLATES, and the translator is each drug's "
      "position on its own Hill curve:")
    a("")
    tbl = []
    for name in DRUGS:
        sub = [x for x in sw["rows"] if x["drug"] == name]
        z = next(x for x in sub if x["strength"] == 0.0)
        m = max(sub, key=lambda x: x["strength"])
        auc_rel = m["mean_auc_uM_s"] / z["mean_auc_uM_s"] if z["mean_auc_uM_s"] else float("nan")
        kill_rel = m["kill_fraction"] / z["kill_fraction"] if z["kill_fraction"] else float("nan")
        tbl.append((name, z["mean_exposure_uM"] / z["ic50_uM"] if z["ic50_uM"] else float("nan"),
                    auc_rel, kill_rel, z["coverage_1.0x"], m["coverage_1.0x"]))
    a("| drug | mean exposure / IC50 at S=0 | AUC retained at S_max | "
      "**kill retained** at S_max | coverage @1x IC50, S=0 -> S_max |")
    a("|---|---|---|---|---|")
    for name, ratio, auc_rel, kill_rel, c0, c1 in tbl:
        a(f"| {name} | {_f(ratio, 4)} | {_f(auc_rel, 4)}x | "
          f"**{_f(kill_rel, 4)}x** | {_f(c0)} -> {_f(c1)} |")
    a("")
    a("A drug sitting FAR BELOW its IC50 is on the steep, near-proportional part "
      "of the Hill curve, so a given fractional exposure loss costs it nearly the "
      "same fraction of its kill. A drug sitting AT or ABOVE its IC50 for much of "
      "the tissue is partly saturated, so the same fractional exposure loss costs "
      "it proportionally less kill — but costs it a great deal of COVERAGE, "
      "because coverage is a threshold count and the threshold is exactly where "
      "its cells are. Temozolomide is the extreme of the first case (its "
      "coverage at 1x IC50 is 0.0000 at every strength — 606 uM is simply "
      "unreachable in this tissue, so the metric cannot move) and doxorubicin of "
      "the second (coverage collapses while kill is comparatively defended). "
      "**The practical consequence: which metric is sensitive to a delivery "
      "insult depends on the drug, and reporting only one of AUC / coverage / "
      "kill would have hidden two-thirds of this effect.**")
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
    lab = fams["label_permutation"]
    clump = fams["clumpiness_matched"]
    block = fams["block_shuffle"]
    jit = fams.get("jitter_reshuffle")
    spot_em = fams.get("spot_emulation")
    floor = dg["uniform_sink_floor"]
    edge = dg["edge_distance"]
    scat = dg["uniform_in_tissue_scatter"]
    ex_lo, ex_hi = hon["honest_excess_min"], hon["honest_excess_max"]
    p_floor = hon["p_floor"]

    a("## 5. THE SPATIAL QUESTION — is the reduction structured, or just smaller?")
    a("")
    a("This is the section the whole module exists for, and it is the section "
      "the audit corrected. **A uniform reduction is pharmacologically boring: "
      "a compartmental model reproduces it exactly by lowering one rate "
      "constant.** The question is whether tumour cells sitting in "
      "macrophage-dense neighbourhoods are measurably worse off than tumour "
      "cells that are not — because that is a statement only a spatial model "
      "can make.")
    a("")
    a(f"**Read the answer first, then the working.** The raw quintile contrast "
      f"is {_f(obs['structure_ratio'], 4)}x and reproduces to the digit. But "
      f"three things have to be subtracted from it before it means what this "
      f"report originally claimed it meant: a depth-into-tissue transport "
      f"gradient that carries {_f(edge['stats']['structure_ratio'], 4)}x on its "
      f"own with no immune input (§5.1a), a null that under-clumps the sink "
      f"field and therefore loses without contest (§5.3), and a monotonicity "
      f"argument that every control also satisfies (§5.4). What survives all "
      f"three is an excess of **{_f(ex_lo, 2)}-{_f(ex_hi, 2)}x** over honest "
      f"spatially-aware nulls"
      + (f", above their entire ranges at p = {_f(p_floor, 4)} (the "
         f"{lab['n_draws']}-draw floor)."
         if hon["above_all_primary_null_draws"] else
         f" — and it does NOT clear their entire ranges, so the effect is not "
         f"separable from them at {lab['n_draws']} draws."))
    a("")
    a(f"**Configuration.** `{sp['drug']}` at S = {sp['strength']:g} (the sweep "
      f"point nearest the Miller anchor), flagship section `{FLAGSHIP_SECTION}`, "
      f"dx = {sp['dx_um']:.0f} um, transvascular exchange off, "
      f"{sp['n_tumor']:,} tumour cells and {sp['n_immune']:,} immune cells. The "
      f"statistic is the per-tumour-cell **attenuation** "
      f"`a_i = AUC_sink,i / AUC_base,i`, from the SAME prologue, source field "
      f"and velocity field — the only difference between the two solves is "
      f"`k_sink`.")
    a("")
    a(f"The observed statistic is computed twice by two independent code paths "
      f"in this file (`spatial_structure` for §5.5's mesh check and "
      f"`_spatial_context` for the nulls below); they agree to "
      f"{_f(abs(hc['observed']['structure_ratio'] - obs['structure_ratio']), 6)} "
      f"absolute"
      + (", which is the internal consistency check on everything that follows."
         if abs(hc['observed']['structure_ratio'] - obs['structure_ratio']) < 1e-6
         else " — which is NOT agreement to round-off, and until that is "
              "explained no §5 number should be quoted."))
    a("")
    if spf is not None and abs(spf["dx_um"] - sp["dx_um"]) > 1e-9:
        a(f"**A declared cost decision, and it is priced.** Every null draw is a "
          f"full transport re-solve, so the nulls run at dx = "
          f"{sp['dx_um']:.0f} um rather than the flagship "
          f"{spf['dx_um']:.0f} um (~2.8x fewer unknowns). Because this project "
          f"treats mesh-dependence as a first-class failure mode, the SAME "
          f"observed statistic is also computed at the flagship spacing and "
          f"both are reported in §5.5 — so the reader can see directly whether "
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
      f"attenuation in the densest) = **{_f(obs['structure_ratio'], 4)}x**. This "
      f"number is correct and reproduces to the digit; what it is compared "
      f"against is where this report previously went wrong. Read §5.1a, §5.3 "
      f"and §5.4 before quoting it.")
    a("")
    a("**The two floors, MEASURED rather than asserted.** The claim \"a "
      "well-mixed model gives exactly 1.0\" is not assumed here: a genuinely "
      "uniform sink field of the same tissue-mean magnitude is built, re-solved, "
      "and binned identically.")
    a("")
    a("| statistic | observed | well-mixed floor, MEASURED | what the excess "
      "over the floor actually is |")
    a("|---|---|---|---|")
    a(f"| structure ratio (sparsest/densest quintile) | "
      f"**{_f(obs['structure_ratio'], 4)}** | "
      f"**{_f(floor['structure_ratio'], 4)}** | "
      + ("the analytic 1.0000 is a FAIR comparator for this binned statistic — "
         "the measured floor sits slightly below it, so quoting 1.0 is if "
         "anything mildly conservative"
         if floor["structure_ratio"] <= 1.0 else
         f"the measured floor is ABOVE 1.0, so the analytic 1.0 OVERSTATES the "
         f"excess; read {_f(obs['structure_ratio'], 4)} against "
         f"{_f(floor['structure_ratio'], 4)}")
      + " |")
    a(f"| per-cell attenuation spread (p90/p10) | "
      f"**{_f(obs['a_p90_over_p10'], 4)}** | "
      f"**{_f(floor['a_p90_over_p10'], 4)}** | "
      f"**the analytic 1.0 is NOT the right floor here.** "
      f"{_f(obs['a_p90_over_p10'], 4)} must be read against "
      f"{_f(floor['a_p90_over_p10'], 4)}, i.e. an excess of "
      f"{_f(hon['spread_excess_over_floor'], 3)}x |")
    a("")
    a(f"The p90/p10 figure of {_f(obs['a_p90_over_p10'], 4)}x (over a full range "
      f"of {_f(obs['a_min'], 4)} to {_f(obs['a_max'], 4)}, against a "
      f"section-mean attenuation of {_f(obs['a_mean'], 4)}) was previously "
      f"quoted as though its floor were 1.0. **It is "
      f"{_f(floor['a_p90_over_p10'], 4)}.** Any external use of the p90/p10 "
      f"number must carry that floor with it.")
    a("")
    a("The well-mixed floor is the easy comparator, and it is not the "
      "informative one. A uniform sink is not a serious rival hypothesis — the "
      "serious rivals are \"the macrophages are clumped, but not clumped *where "
      "the data says*\" (§5.3) and \"this is depth into tissue, not immune "
      "density\" (§5.1a).")
    a("")
    a("#### 5.1a A large part of the spread is DEPTH INTO TISSUE, not immune biology")
    a("")
    a("This is the largest single correction in this report and it was not "
      "tested in the original version.")
    a("")
    a(f"Take the same {sp['n_tumor']:,} tumour cells, the same two solves, and "
      f"bin them by **distance to the tissue edge alone** — no reference to "
      f"immune density anywhere in the binning:")
    a("")
    a("| binning variable | sparsest/shallowest bin | densest/deepest bin | "
      "**structure ratio** |")
    a("|---|---|---|---|")
    a(f"| local immune density `phi_imm` (table above) | "
      f"{_f(obs['bins'][0]['attenuation_mean'], 4)} | "
      f"{_f(obs['bins'][-1]['attenuation_mean'], 4)} | "
      f"**{_f(obs['structure_ratio'], 4)}** |")
    a(f"| **distance to tissue edge, no immune input at all** | "
      f"**{_f(edge['stats']['bins'][0]['attenuation_mean'], 4)}** | "
      f"**{_f(edge['stats']['bins'][-1]['attenuation_mean'], 4)}** | "
      f"**{_f(edge['stats']['structure_ratio'], 4)}** |")
    a("")
    a(f"So a variable containing **zero immune information** already reaches "
      f"{_f(edge['stats']['structure_ratio'], 4)} of the headline "
      f"{_f(obs['structure_ratio'], 4)} — "
      f"{100.0 * edge['stats']['structure_ratio'] / obs['structure_ratio']:.0f}% "
      f"of the ratio. The two binning variables are correlated: "
      f"**corr(phi_imm, edge-distance) = "
      f"{_f(edge['corr_phi_edge_pearson'], 3)}** (Spearman "
      f"{_f(edge['corr_phi_edge_spearman'], 3)}).")
    a("")
    a(f"**The depth gradient comes from the measured tissue, not from the "
      f"kernel.** The segmented cell-density envelope rises "
      f"**{_f(edge['density_envelope_ratio'], 1)}x** from the outer "
      f"{edge['rim_um']:.0f} um band into the interior. Cells deep in the "
      f"section are simply in denser tissue, which carries both more macrophages "
      f"and different transport. The uniform-in-tissue scatter control (§5.4) "
      f"passes the identical KDE and reaches only "
      f"{_f(scat['structure_ratio'], 4)}, so the KDE is not manufacturing the "
      f"gradient by itself.")
    a("")
    a(f"**And the published label-permutation null cannot remove this confound, "
      f"because it PRESERVES it.** Permuting the immune label among measured "
      f"cell positions keeps every immune cell on a real cell, so it inherits "
      f"the same density envelope: the permuted fields carry "
      f"corr(phi_imm, edge-distance) = "
      f"**{_f(lab['corr_phi_edge_mean'], 3)}** against the observed "
      f"{_f(edge['corr_phi_edge_pearson'], 3)}. The null reproduces the depth "
      f"gradient it should be controlling for.")
    a("")
    a("**How much survives after stratifying it out.** Splitting the tumour "
      "cells into five edge-distance strata and computing the immune-density "
      "structure ratio *within* each:")
    a("")
    a("| edge-distance stratum | " + " | ".join(
        f"{i + 1}{' (shallowest)' if i == 0 else ' (deepest)' if i == N_BINS - 1 else ''}"
        for i in range(N_BINS)) + " |")
    a("|---|" + "---|" * N_BINS)
    a("| within-stratum structure ratio | " + " | ".join(
        f"**{_f(v, 2)}**" for v in edge["within_stratum_sr"]) + " |")
    a("")
    _strat = [v for v in edge["within_stratum_sr"] if np.isfinite(v)]
    a(f"An immune-density contrast persists at every depth, so the effect is not "
      f"purely depth. But it is {_f(min(_strat), 2)}-{_f(max(_strat), 2)}x "
      f"within strata rather than {_f(obs['structure_ratio'], 4)}x overall, and "
      f"it is not monotone in depth — which is what one expects when two "
      f"partly-collinear variables are being separated on a single section.")
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
    a("**This table is NOT independent corroboration of Miller's length scale, "
      "and must not be read as such.** The kernel bandwidth was CHOSEN as "
      f"{IMMUNE_KERNEL_BANDWIDTH_UM:.0f} um from that same measurement (§2.1), "
      "so the model was told what length scale to have and the fact that it has "
      "it is arithmetic, not evidence.")
    a("")
    a("**Nor is the monotonicity of this column evidence of anything** — §5.4 "
      "shows every control reproduces monotonicity, including controls with no "
      "measured immune structure in them at all. The original version of this "
      "report leaned on that monotonicity; it should not have.")
    a("")
    a("What the table does show that is not circular is that the deficit does "
      "not saturate until well beyond the bandwidth (drug also diffuses in from "
      "further away, which is exactly the non-local effect an AUC rescale could "
      "not have produced), and that the worst-off cells are a small minority — "
      f"{sp['distance_bins'][0]['n']:,} of {sp['n_tumor']:,} tumour cells lie "
      "within 10 um of a macrophage. Note also that the sink-OFF AUC column "
      "falls with distance too: cells near macrophages are in denser, "
      "better-perfused tissue and started with MORE drug, so the attenuation "
      "column is the right comparison and the raw exposure column would point "
      "the other way. That same tissue-density envelope is the confound "
      "quantified in §5.1a.")
    a("")
    a("### 5.3 The nulls — the published one was too weak, and here are the honest ones")
    a("")
    a("**What was wrong.** The original version of this report tested the "
      "observed structure against ONE null: permuting the IMMUNE LABEL among "
      "the measured cell positions. That null does hold the tissue outline, the "
      "cell layout, the vessel geometry, the ECM field and the total sink "
      "magnitude fixed — but it destroys the immune compartment's clumping "
      "entirely, and clumping is the very quantity the structure ratio "
      "measures.")
    a("")
    a(f"**Measured, not argued.** The sink field the label-permutation null "
      f"produces has coefficient of variation **{_f(lab['cv_mean'], 3)}** "
      f"(range {_f(lab['cv_min'], 3)}-{_f(lab['cv_max'], 3)}). The observed sink "
      f"field has CV **{_f(hon['observed_sink_cv'], 3)}**. A null whose field is "
      f"{100.0 * (1.0 - lab['cv_mean'] / hon['observed_sink_cv']):.0f}% less "
      f"variable than the observed field will lose a variance-driven contest "
      f"almost regardless of whether any real spatial biology exists — it is not "
      f"being asked a fair question.")
    a("")
    a("And the dependence is direct rather than incidental. Sweeping the clump "
      "size of the permutation, holding the immune COUNT and therefore the total "
      "sink magnitude fixed:")
    a("")
    a("| clump size (cells) | " + " | ".join(
        f"{r['clump_size']}" for r in dg["cv_sweep"]) + " |")
    a("|---|" + "---|" * len(dg["cv_sweep"]))
    a("| sink-field CV | " + " | ".join(
        f"{_f(r['sink_cv'], 2)}" for r in dg["cv_sweep"]) + " |")
    a("| resulting structure ratio | " + " | ".join(
        f"**{_f(r['structure_ratio'], 2)}**" for r in dg["cv_sweep"]) + " |")
    a("")
    a(f"**The structure ratio very largely tracks sink-field heterogeneity.** "
      f"Any null that does not match the observed CV of "
      f"{_f(hon['observed_sink_cv'], 3)} is measuring CV, not geometry.")
    a("")
    a("#### The primary comparators")
    a("")
    a("Two nulls fix this. Both keep **every immune cell on a real measured cell "
      "position** — nothing is scattered into empty space, and nothing leaves "
      "the tissue mask:")
    a("")
    a(f"- **clumpiness-matched permutation** — the immune label moves in clumps "
      f"of {CLUMP_SIZE} mutually-nearest MEASURED cells instead of one cell at a "
      f"time, so the CV asymmetry above is largely removed "
      f"({_f(clump['cv_mean'], 3)} against the observed "
      f"{_f(hon['observed_sink_cv'], 3)}) and only the *location* of the clumps "
      f"varies.")
    a(f"- **{BLOCK_UM:.0f} um block shuffle** — the section is tiled into "
      f"{BLOCK_UM:.0f} um blocks and each block's observed immune COUNT is "
      f"reassigned to another block, which draws that many of its own measured "
      f"cells. The block-scale count heterogeneity of the real section is "
      f"inherited exactly and only relocated.")
    a("")
    a(f"{lab['n_draws']} draws each, each one a full transport re-solve.")
    a("")
    a("| null | rank | null mean SR | null range | sink-field CV | observed | "
      "one-sided p | **observed / null mean** |")
    a("|---|---|---|---|---|---|---|---|")
    _order = [("clumpiness_matched", "**PRIMARY**"),
              ("block_shuffle", "**PRIMARY**"),
              ("label_permutation", "DEMOTED — floor check only"),
              ("jitter_reshuffle", "control, see below"),
              ("spot_emulation", "control, see below")]
    for key, rank in _order:
        f_ = fams.get(key)
        if f_ is None:
            continue
        star = "**" if f_["is_primary"] else ""
        a(f"| {star}{f_['label']}{star} | {rank} | {star}{_f(f_['sr_mean'], 4)}{star} | "
          f"{_f(f_['sr_min'], 4)}-{_f(f_['sr_max'], 4)} | "
          f"{_f(f_['cv_mean'], 3)} | {_f(obs['structure_ratio'], 4)} | "
          f"{_f(f_['p_structure_ratio'], 4)} | "
          f"{star}{_f(f_['excess_over_null'], 3)}x{star} |")
    a(f"| torus shift | **UNUSABLE on this mask** — see below | -- | -- | "
      f"{_f(dg['torus']['cv_ratio'] * hon['observed_sink_cv'], 3)} | "
      f"{_f(obs['structure_ratio'], 4)} | -- | -- |")
    a("")
    a(f"**The honest number is the last column: "
      f"{_f(obs['structure_ratio'], 4)} / {_f(clump['sr_mean'], 4)} = "
      f"{_f(clump['excess_over_null'], 3)}x and "
      f"{_f(obs['structure_ratio'], 4)} / {_f(block['sr_mean'], 4)} = "
      f"{_f(block['excess_over_null'], 3)}x.** That — an excess of roughly "
      f"{_f(ex_lo, 2)} to {_f(ex_hi, 2)}x — is what is attributable to where the "
      f"macrophages actually are, once a null is allowed to be as clumped as the "
      f"data. It is not {_f(obs['structure_ratio'], 1)}x. This report previously "
      f"said {_f(obs['structure_ratio'], 1)}x, and that was quoting a raw ratio "
      f"against a floor of 1.0 as though the floor were the alternative "
      f"hypothesis.")
    a("")
    if hon["above_all_primary_null_draws"]:
        a(f"**What still holds, and it is not nothing.** The observed "
          f"{_f(obs['structure_ratio'], 4)} sits **above the entire range of "
          f"both primary nulls** ({_f(clump['sr_min'], 4)}-"
          f"{_f(clump['sr_max'], 4)} and {_f(block['sr_min'], 4)}-"
          f"{_f(block['sr_max'], 4)}), at p = {_f(p_floor, 4)}. A "
          f"clumpiness-matched null and a block-shuffle null both fail to reach "
          f"the observed value in {lab['n_draws']} out of {lab['n_draws']} "
          f"draws. Putting the same number of equally-clumped macrophages "
          f"*somewhere else* on real measured cell positions produces a "
          f"measurably weaker exposure gradient than putting them where the data "
          f"says they are. A {_f(ex_lo, 2)}-{_f(ex_hi, 2)}x effect that survives "
          f"an honest null is worth more than a "
          f"{_f(obs['structure_ratio'], 1)}x effect resting on a weak one.")
    else:
        a(f"**And it does not clear the honest nulls.** The observed "
          f"{_f(obs['structure_ratio'], 4)} does NOT sit above the entire range "
          f"of both primary nulls ({_f(clump['sr_min'], 4)}-"
          f"{_f(clump['sr_max'], 4)} and {_f(block['sr_min'], 4)}-"
          f"{_f(block['sr_max'], 4)}). Once a null is allowed to be as clumped "
          f"as the data, the measured immune POSITIONS are not distinguishable "
          f"from an equally-clumped rearrangement of them at "
          f"{lab['n_draws']} draws. Nothing here should be quoted as a "
          f"macrophage-geometry effect.")
    a("")
    a("#### The platform controls — what happens if the positional detail goes away")
    a("")
    def _control_verdict(f_: Dict, what: str, scale: str) -> str:
        """Read a control by how close it lands to the observed value."""
        r = f_["sr_mean"] / obs["structure_ratio"]
        if r >= 0.95:
            return (f"it **reproduces the observed value** "
                    f"({_f(f_['sr_mean'], 4)} against "
                    f"{_f(obs['structure_ratio'], 4)}, "
                    f"{_f(r, 3)}x). {what[0].upper() + what[1:]} therefore "
                    f"carries essentially none of "
                    f"the statistic: the structure ratio is built from immune "
                    f"density at scales COARSER than {scale}, not from it. This "
                    f"control does not discriminate, and nothing is rested on "
                    f"it.")
        if r >= 0.8:
            return (f"it lands close to but below the observed value "
                    f"({_f(f_['sr_mean'], 4)} against "
                    f"{_f(obs['structure_ratio'], 4)}, {_f(r, 3)}x), so {what} "
                    f"carries a minority of the statistic.")
        return (f"it **collapses** ({_f(f_['sr_mean'], 4)} against "
                f"{_f(obs['structure_ratio'], 4)}, {_f(r, 3)}x), so the "
                f"statistic genuinely depends on {what}.")

    if jit is not None:
        a(f"**Jitter reshuffle — re-draw the Visium jitter and nothing else.** "
          f"Every measured immune cell is displaced by a uniform-in-disc offset "
          f"of at most {SPOT_JITTER_UM:g} um, the exact amplitude "
          f"`reconstruct/real_visium.py::_place_cells` applies "
          f"(`jitter_um = min(27.5, 0.45 * dx_um)`). Nothing about which "
          f"neighbourhood a macrophage is in changes. On a SPOT reconstruction "
          f"this control is a no-op by construction — the positions were "
          f"generated that way. Here, {lab['n_draws']} draws give "
          f"**{_f(jit['sr_mean'], 4)}** (range {_f(jit['sr_min'], 4)}-"
          f"{_f(jit['sr_max'], 4)}, CV {_f(jit['cv_mean'], 3)} against the "
          f"observed {_f(hon['observed_sink_cv'], 3)}): "
          + _control_verdict(jit, "positional detail below the jitter amplitude",
                             f"~{SPOT_JITTER_UM:g} um"))
        a("")
    if spot_em is not None:
        a(f"**Spot emulation — the whole Visium synthesis.** Snapping every "
          f"measured immune cell to its {SPOT_PITCH_UM:.0f} um lattice centroid "
          f"and re-jittering it (literally what `_place_cells` does, keeping the "
          f"per-spot immune COUNT and destroying everything finer) gives "
          f"**{_f(spot_em['sr_mean'], 4)}** at CV {_f(spot_em['cv_mean'], 3)}: "
          + _control_verdict(spot_em,
                             f"positional detail below the "
                             f"{SPOT_PITCH_UM:.0f} um spot pitch",
                             f"the {SPOT_PITCH_UM:.0f} um spot pitch"))
        a("")
    if jit is not None and jit["sr_mean"] >= 0.95 * obs["structure_ratio"]:
        a("**This does NOT reproduce a claim the §5 audit made, and the audit's "
          "version should be dropped rather than repeated.** The audit reported "
          "a jitter reshuffle at 1.4188 on this section and read its collapse as "
          "the decisive evidence separating Xenium from Visium. Re-implemented "
          "here — displace every measured immune cell by at most "
          f"{SPOT_JITTER_UM:g} um, change nothing else — the control returns "
          f"{_f(jit['sr_mean'], 4)}, i.e. the observed value, and the full spot "
          f"synthesis returns {_f(spot_em['sr_mean'], 4)} if run as well. Either "
          f"the audit's control was a different operation or its number does not "
          f"reproduce; either way **this report does not rest the Xenium result "
          f"on a jitter control**, and no downstream document should. What the "
          f"Xenium result rests on is the two primary nulls in the table above, "
          f"which do not involve the platform question at all. The scope limit "
          f"in §5.6 rests on PROVENANCE — spot platforms never localise an "
          f"immune cell — and not on this control.")
        a("")
    a("#### Two methodological notes that must travel with these numbers")
    a("")
    a(f"1. **Torus / Moran-style nulls are unusable on this tissue mask and "
      f"nothing rests on them.** Circularly shifting the immune point pattern "
      f"lands only **{100.0 * dg['torus']['fraction_in_tissue_mean']:.1f}%** of "
      f"shifted cells back inside tissue (over {dg['torus']['n_probe']} random "
      f"shifts; range {100.0 * dg['torus']['fraction_in_tissue_min']:.1f}-"
      f"{100.0 * dg['torus']['fraction_in_tissue_max']:.1f}%). The surviving "
      f"field has mean `phi_imm` at "
      f"**{_f(dg['torus']['phi_mean_ratio'], 2)}x** the observed and CV inflated "
      f"to **{_f(dg['torus']['cv_ratio'], 2)}x** the observed — i.e. the shift "
      f"simultaneously weakens and sharpens the sink, which biases the structure "
      f"ratio upward for entirely mechanical reasons. **This is measured here as "
      f"a FEASIBILITY diagnostic and the torus null is never scored.** The "
      f"clumpiness-matched permutation and the block shuffle are the defensible "
      f"nulls, precisely because both keep every immune cell on a real measured "
      f"position inside the mask.")
    a(f"2. **p = {_f(p_floor, 4)} is the {lab['n_draws']}-draw floor** "
      f"(1/{lab['n_draws'] + 1}). It appears in every row above for that reason "
      f"and **must never be quoted as a resolved p-value**. It means only \"not "
      f"exceeded in {lab['n_draws']} draws\". The informative columns are the "
      f"null range and the observed/null ratio.")
    a("")
    a(f"**One control that is not about structure at all.** Mean attenuation is "
      f"a magnitude, not a structure statistic: observed "
      f"{_f(obs['a_mean'], 4)} against a label-permutation mean of "
      f"{_f(lab['a_mean_mean'], 4)} and a clumpiness-matched mean of "
      f"{_f(clump['a_mean_mean'], 4)}. The same number of macrophages removes "
      f"roughly the same total amount of drug wherever they sit. A compartmental "
      f"model gets that quantity right by construction; what it cannot express "
      f"is **which** tumour cells lose the drug. That distinction is still the "
      f"point of the module — it is the *size* of the spatial term that has been "
      f"corrected, not its existence.")
    a("")
    a("### 5.4 Monotonicity across the quintiles is NOT evidence, and was "
      "previously treated as though it were")
    a("")
    a(f"The original report pointed to the attenuation falling smoothly across "
      f"all five immune quintiles "
      f"({_f(obs['bins'][0]['attenuation_mean'], 4)} -> "
      f"{_f(obs['bins'][-1]['attenuation_mean'], 4)}) as corroboration. It is "
      f"not. Every control tested is also monotone across all five quintiles:")
    a("")
    a("| field | quintile 1 | quintile 5 | monotone across all 5? | ratio |")
    a("|---|---|---|---|---|")
    _obs_mono = all(obs["bins"][i]["attenuation_mean"]
                    > obs["bins"][i + 1]["attenuation_mean"]
                    for i in range(N_BINS - 1))
    a(f"| **observed** | {_f(obs['bins'][0]['attenuation_mean'], 4)} | "
      f"{_f(obs['bins'][-1]['attenuation_mean'], 4)} | "
      f"{'yes' if _obs_mono else 'no'} | {_f(obs['structure_ratio'], 4)} |")
    for key in ("label_permutation", "clumpiness_matched", "block_shuffle",
                "jitter_reshuffle", "spot_emulation"):
        f_ = fams.get(key)
        if f_ is None:
            continue
        q = f_["quintile_attenuation_mean"]
        a(f"| {f_['label']} | {_f(q[0], 4)} | {_f(q[-1], 4)} | "
          f"**{'yes' if f_['monotone'] else 'no'}** | {_f(f_['sr_mean'], 2)} |")
    _scat_mono = all(scat["bins"][i]["attenuation_mean"]
                     > scat["bins"][i + 1]["attenuation_mean"]
                     for i in range(N_BINS - 1))
    a(f"| uniform-in-tissue scatter | "
      f"{_f(scat['bins'][0]['attenuation_mean'], 4)} | "
      f"{_f(scat['bins'][-1]['attenuation_mean'], 4)} | "
      f"**{'yes' if _scat_mono else 'no'}** | "
      f"{_f(scat['structure_ratio'], 2)} |")
    a("")
    a("Uniform-in-tissue scatter contains **no measured immune structure "
      "whatsoever** and is still monotone. Monotonicity is a generic consequence "
      "of binning a smooth field by its own local value — it would be surprising "
      "only in its absence. **Only the RATIO carries information, and only "
      "against a null that matches the observed field's clumpiness (§5.3).** Any "
      "downstream document quoting \"monotone across all five quintiles\" as "
      "support should drop that phrasing.")
    a("")
    if spf is not None:
        a("### 5.5 Is the structure a property of the tissue, or of the mesh?")
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
            tail = " (flagship)" if abs(tag["dx_um"] - DX_FLAGSHIP) < 1e-9 else " (nulls run here)"
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
             "was computed at.")
          + " Note that mesh-stability is a necessary check, not a sufficient "
            "one: a discretisation-stable number can still be a depth artefact "
            "(§5.1a) or be quoted against the wrong null (§5.3), and here it "
            "was both.")
        a("")
    a("### 5.6 SCOPE LIMIT — this statistic is INVALID on Visium / spot "
      "transcriptomics, and the generator now REFUSES to compute it there")
    a("")
    a("This is a hard limit, not a caution. **The structure ratio computed on "
      "Visium or any other spot-based platform is meaningless, because the "
      "platform does not measure immune-cell positions at all.**")
    a("")
    a(f"**The mechanism, read off the reconstruction code.** "
      f"`tmesim/reconstruct/real_visium.py::_place_cells` builds the cell "
      f"population from a deconvolved spot table by synthesising a fixed "
      f"**{_VISIUM_CELLS_PER_SPOT} pseudo-cells per ~{SPOT_PITCH_UM:.0f} um "
      f"spot** (`_CELLS_PER_SPOT = {_VISIUM_CELLS_PER_SPOT}`), apportioned "
      f"across types by the deconvolved fractions, and placing all of them **at "
      f"the spot centroid with uniform-in-disc jitter of at most "
      f"{SPOT_JITTER_UM:g} um** (`jitter_um = min(27.5, 0.45 * dx_um)`). There "
      f"is no measured position for any individual immune cell anywhere in that "
      f"pipeline. The only spatial information present is the per-spot immune "
      f"FRACTION on a {SPOT_PITCH_UM:.0f} um lattice; everything finer is drawn "
      f"from a random number generator.")
    a("")
    a(f"**Why every plausible null is therefore a no-op.** The sink field is a "
      f"{IMMUNE_KERNEL_BANDWIDTH_UM:.0f} um-bandwidth KDE. Sub-spot jitter of "
      f"<={SPOT_JITTER_UM:g} um is at or below that bandwidth, so re-drawing it "
      f"— or reshuffling labels within a spot — changes the KDE almost not at "
      f"all. Any control that preserves the per-spot immune count is, by "
      f"construction, comparing a field to itself.")
    a("")
    a("**The refusal is in the code, not in this prose.** "
      "`validation/immune_sink_impact.py::require_measured_cell_positions` "
      "raises `SpotDataError` for any reconstruction whose own metadata says "
      "the positions were synthesised from spots "
      "(`metadata['real_pipeline']['source']`, `metadata['n_spots']`); the "
      "shared context builder every null family runs through calls it first, so "
      "no null can bypass it, and an unstamped reconstruction of unknown "
      "provenance is refused too. It reads the reconstruction's metadata rather "
      "than the sample PATH, so a renamed section cannot get past it. "
      "`tests/test_immune_sink_nulls.py` pins the refusal.")
    a("")
    if spot_em is not None:
        _r = spot_em["sr_mean"] / obs["structure_ratio"]
        a(f"**And the cost of ignoring the limit is measured here rather than "
          f"asserted.** Passing the measured single-cell immune positions "
          f"through the identical spot synthesis (§5.3, spot emulation) returns "
          f"**{_f(spot_em['sr_mean'], 4)}** against the observed "
          f"{_f(obs['structure_ratio'], 4)} — {_f(_r, 3)}x. "
          + ("**A spot-resolution field reproduces essentially the entire "
             "claimed effect.** Read that as a warning and not as a "
             "reassurance: it means the structure ratio is carried by immune "
             "density at the spot scale and coarser, so on a platform where "
             "the per-spot immune content is a DECONVOLUTION ESTIMATE rather "
             "than a count of localised cells, the statistic will report "
             "whatever the deconvolution wrote — at full apparent strength, "
             "and with no null available on that platform able to tell the two "
             "apart, because every control that preserves per-spot counts is a "
             "no-op there."
             if _r >= 0.95 else
             "The spot pipeline does not reproduce the measured section's "
             "value, which is consistent with the statistic needing positional "
             "detail a spot platform does not have.")
          + " Any use of a Visium-derived structure ratio in a deck, a landing "
            "page or a customer-facing claim is unsupported.")
        a("")
    a("**This does not retract the single-cell result, and it does not rest on "
      "a platform control either.** What that result rests on is the two "
      "primary nulls of §5.3, which never leave the measured section: an "
      "equally-clumped rearrangement of the same macrophages onto other "
      "measured positions "
      + ("does not reach the observed value in any draw"
         if hon["above_all_primary_null_draws"] else
         "reaches the observed value in at least one draw")
      + ". The scope limit here is about PROVENANCE — a spot platform never "
        "localises an immune cell, so there is nothing on it for such a null "
        "to rearrange — and the generator enforces it by refusing to run.")
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
        a(f"3. **The spatial result is real but is SMALLER than this report "
          f"originally claimed, and must be quoted at its corrected size.** The "
          f"raw quintile contrast is {_f(obs['structure_ratio'], 4)}x, but the "
          f"defensible statement is: **an excess of "
          f"{_f(ex_lo, 2)}-{_f(ex_hi, 2)}x over spatially-aware nulls that "
          f"match the observed clumpiness (§5.3), above the entire range of "
          f"both in {_lab['n_draws']} draws.** Most of the raw "
          f"{_f(obs['structure_ratio'], 4)} is accounted for by things that are "
          f"not immune geometry: depth into tissue alone reaches "
          f"{_f(_edge['stats']['structure_ratio'], 4)}x with no immune input "
          f"whatsoever (§5.1a), and simply matching the observed clumpiness "
          f"reaches {_f(_clump['sr_mean'], 4)}x (§5.3). The correct external "
          f"phrasing is *\"tumour cells in macrophage-dense neighbourhoods lose "
          f"measurably more drug than an equally-clumped random arrangement of "
          f"the same macrophages would predict, by "
          f"~{_f(ex_lo, 2)}-{_f(ex_hi, 2)}x\"* — not *\"a "
          f"{_f(obs['structure_ratio'], 1)}x spread that a compartmental model "
          f"gives as 1.0\"*. It remains a genuine demonstration that geometry "
          f"carries information a compartmental model cannot express; it is not "
          f"a {_f(obs['structure_ratio'], 1)}x demonstration. **And it is "
          f"single-cell-only: on spot data the statistic is invalid and the "
          f"generator refuses to compute it (§5.6).**")
    else:
        a(f"3. **The reduction is structured, but the structure is not "
          f"MACROPHAGE-specific once the null is allowed to be as clumped as "
          f"the data.** The raw {_f(obs['structure_ratio'], 4)}x contrast is "
          f"real, but the clumpiness-matched permutation reaches "
          f"{_f(_clump['sr_mean'], 4)} and the block shuffle "
          f"{_f(_block['sr_mean'], 4)} without using the measured macrophage "
          f"positions at all, and depth into tissue alone reaches "
          f"{_f(_edge['stats']['structure_ratio'], 4)}. What the analysis "
          f"establishes is that exposure is structured by local CELLULARITY and "
          f"by depth. Quoting it as a macrophage effect would over-claim, and "
          f"this report does not.")
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
    a(f"7. **The spatial structure ratio is partly confounded with depth into "
      f"tissue, and one section cannot fully separate them.** §5.1a measures "
      f"corr(phi_imm, edge-distance) = "
      f"{_f(_edge['corr_phi_edge_pearson'], 3)} on a cell-density envelope that "
      f"rises {_f(_edge['density_envelope_ratio'], 1)}x from the tissue rim "
      f"inward. Stratifying gives "
      f"{_f(min(v for v in _edge['within_stratum_sr'] if np.isfinite(v)), 2)}-"
      f"{_f(max(v for v in _edge['within_stratum_sr'] if np.isfinite(v)), 2)}x "
      f"within strata, but with {sp['n_tumor']:,} tumour cells on ONE section "
      f"there is no way to decide how much of the residual is immune geometry "
      f"and how much is residual collinearity. A second single-cell section "
      f"with a different immune-versus-depth arrangement is the cheapest test "
      f"that would resolve it, and it has not been run.")
    a(f"8. **Everything in §5 is one section, one drug, one strength.** The "
      f"corrected {_f(ex_lo, 2)}-{_f(ex_hi, 2)}x is a single measurement with a "
      f"p-value pinned at its {_lab['n_draws']}-draw floor "
      f"({_f(p_floor, 4)}). It is not an interval and no confidence should be "
      f"attached to the second decimal place.")
    a("9. **The statistic does not transfer to spot platforms** (§5.6), so the "
      "module's flagship spatial claim is available only where single-cell "
      "centroids are measured. That restricts it to Xenium/CosMx-class data and "
      "excludes the Visium sections in `data_cache/`; the generator raises "
      "rather than reporting a number there.")
    a("10. **The scorecard test rests on assumption A1** (tissue:plasma "
      "partition = effective penetrance) and on a quasi-steady-state reading of "
      "it. Both are the scorecard's own standing assumptions, inherited "
      "unchanged.")
    a("")
    a("## 9. Reproduction")
    a("")
    a("```bash")
    a("./.venv/bin/python -m tmesim.validation.immune_sink_impact")
    a("```")
    a("")
    _n_null_solves = sum(f_["n_draws"] for f_ in fams.values())
    a(f"Deterministic — the engine's seeded generator (`seed = 0`) and the "
      f"nulls' own seeds ({', '.join(str(v) for v in sorted(set(NULL_SEEDS.values())))}) "
      f"are the only randomness. Runtime is dominated by transport solves: "
      f"{len(DRUGS) * (len(SEQUESTRATION_SWEEP) + 1)} on the {sw['grid']} "
      f"flagship grid (§4, including the three engine reference runs)"
      + (f", {sp['n_perm'] + 2} for the demoted label-permutation null and "
         f"{_n_null_solves + len(CLUMP_SWEEP_SIZES) + 4} more at dx = "
         f"{sp['dx_um']:.0f} um for the honest nulls, the floors and the CV "
         f"sweep (§5.1-§5.4), and 2 at dx = {spf['dx_um']:.0f} um (§5.5)"
         if spf is not None else "")
      + ". Each stage is cached to `results/immune_sink_impact_cache.pkl` as it "
      "completes, so an interrupted run resumes rather than restarting; delete "
      "that file to force a full recompute. **Every §5 number above — the "
      "nulls, the CV diagnostic, the measured floors, the edge-distance "
      "confound, the monotonicity controls and the platform controls — is "
      "produced by this command.** `write_report` refuses to write this file at "
      "all if the honest nulls are missing, so the corrected claim cannot be "
      "silently reverted by a regeneration.")
    a("")
    a("**Files added by this work:**")
    a("")
    a("* `tmesim/molecular/immune_sink.py` — new; the field, the swept magnitude, "
      "the literature table and the re-solve")
    a("* `tmesim/validation/immune_sink_impact.py` — new; this analysis, "
      "including the honest nulls of §5.3 and the spot-data refusal of §5.6")
    a("* `tests/test_immune_sink.py` — new; pins the bit-identity and the "
      "non-vacuity")
    a("* `tests/test_immune_sink_nulls.py` — new; pins that the primary nulls "
      "exceed the label permutation, that the CV diagnostic is emitted, and "
      "that the spot-data guard fires")
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

    print("[1/8] identity check ...", flush=True)
    ident = _stage("identity", identity_check)
    print(f"      {ident['n_bit_identical']}/{ident['n_comparisons']} "
          f"bit-identical, all={ident['all_bit_identical']}", flush=True)

    print("[2/8] immune fields on three real sections ...", flush=True)
    fields = _stage("fields", immune_fields)
    for r in fields:
        print(f"      {r['section']:<16} immune {r['immune_pct']:.2f}%  "
              f"mean phi {r['phi_mean']:.4f}", flush=True)

    print("[3/8] sweep on the flagship section ...", flush=True)
    sw = _stage("sweep", sweep)

    print(f"[4/8] spatial structure + permutation null (dx={DX_PERM:.0f}) ...",
          flush=True)
    sp = _stage("spatial",
                lambda: spatial_structure(dx_um=DX_PERM, n_perm=n_perm))
    print(f"      structure ratio {sp['observed']['structure_ratio']:.4f} "
          f"p={sp['p_structure_ratio']:.4f}", flush=True)

    print(f"[5/8] mesh cross-check of the same statistic at the flagship "
          f"dx={DX_FLAGSHIP:.0f} ...", flush=True)
    spf = _stage("spatial_flagship",
                 lambda: spatial_structure(dx_um=DX_FLAGSHIP, n_perm=0))
    print(f"      structure ratio {spf['observed']['structure_ratio']:.4f}",
          flush=True)

    # --- the honest nulls: built on ONE shared context, cached per family -----
    _ctx: List[Dict] = []

    def ctx() -> Dict:
        if not _ctx:
            print("      building the shared spatial context "
                  f"(dx={DX_PERM:.0f}) ...", flush=True)
            _ctx.append(_spatial_context(dx_um=DX_PERM))
        return _ctx[0]

    print(f"[6/8] honest-null diagnostics: floors, depth confound, CV sweep ...",
          flush=True)
    diag = _stage("honest_diag", lambda: spatial_diagnostics(ctx()))
    print(f"      well-mixed SR floor "
          f"{diag['uniform_sink_floor']['structure_ratio']:.4f}, "
          f"edge-only SR "
          f"{diag['edge_distance']['stats']['structure_ratio']:.4f}", flush=True)

    print(f"[7/8] the null families ({N_PERM} re-solves each) ...", flush=True)
    families: Dict[str, Dict] = {}
    for fam in NULL_FAMILIES:
        families[fam] = _stage(f"null_{fam}",
                               lambda f=fam: null_family(ctx(), f, n_draws=n_perm))
        print(f"      {fam:<20} mean SR {families[fam]['sr_mean']:.4f} "
              f"[{families[fam]['sr_min']:.4f}-{families[fam]['sr_max']:.4f}] "
              f"CV {families[fam]['cv_mean']:.3f} "
              f"excess {families[fam]['excess_over_null']:.3f}x", flush=True)

    hc = _stage("honest_ctx", lambda: {
        "observed": ctx()["observed"],
        "provenance": ctx()["provenance"],
        "n_tumor": ctx()["n_tumor"], "n_immune": ctx()["n_immune"],
        "n_cells": ctx()["n_cells"], "dx_um": ctx()["dx_um"],
        "strength": ctx()["strength"], "drug_name": ctx()["drug_name"],
        "sample": ctx()["sample"],
        "corr_phi_edge_observed": float(np.corrcoef(
            ctx()["phi_at_tumor"], ctx()["edge_at_tumor"])[0, 1]),
    })
    hon = honest_spatial_summary(hc, families, diag)
    out["honest_summary"] = hon
    print(f"      honest excess {hon['honest_excess_min']:.3f}-"
          f"{hon['honest_excess_max']:.3f}x over the primary nulls", flush=True)
    _ctx.clear()

    print("[8/8] glioma scorecard ...", flush=True)
    sc = _stage("scorecard", lambda: glioma_scorecard(sw["mean_phi"]))
    print(f"      baseline {sc['baseline']['gm_abs']:.3f}x -> best "
          f"{sc['best_modality_aware']['gm_abs']:.3f}x "
          f"(improves={sc['improves']})", flush=True)

    if write_md:
        out["report_path"] = write_report(ident, fields, sw, sp, sc, spf,
                                          hc, families, diag, hon)
        print(f"wrote {out['report_path']}", flush=True)
    print(f"total {time.time() - t0:.0f}s", flush=True)
    return out


if __name__ == "__main__":
    run()
