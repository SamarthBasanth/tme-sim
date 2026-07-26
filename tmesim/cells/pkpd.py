"""Pharmacodynamics: map local drug exposure to per-cell kill.

This module consumes a :class:`~tmesim.core.results.TransportResult` (produced by
``tmesim.physics.diffusion``) and decides, cell by cell, whether each cell is
killed by the drug. It is the point where *spatial drug delivery* meets
*intrinsic tumor biology* — the crux of the Stage 0 kill-test.

Exposure -> kill mapping
------------------------
The transport solve returns ``transport.auc`` (ny, nx): the time-integrated drug
concentration ``integral C dt`` at each grid point, in units of ``uM * s``. Each
cell samples the AUC at its own grid index.

We first reduce the AUC (uM*s) to a *time-averaged effective concentration*::

    E = AUC / T          # uM, where T is the exposure window duration (s)

with ``T = transport.times[-1]`` (the length of the simulated window). Working in
a concentration allows us to reuse the drug's IC50 directly as the half-effect
point, which is how IC50 is defined pharmacologically.

The effective concentration is turned into a kill PROBABILITY with a Hill
function (a sigmoidal dose-response curve)::

    p_hill = max_kill_rate * E^h / (E50^h + E^h)

where

* ``h  = drug.hill``                 — the Hill coefficient (curve steepness),
* ``E50 = drug.ic50``                — **E50 derivation**: because ``E`` is a
  time-averaged concentration in uM, the drug's IC50 (the concentration giving
  half-maximal effect, also uM) IS the exposure that yields half of
  ``max_kill_rate``. No unit conversion is needed; ``E50 = drug.ic50`` directly.
* ``max_kill_rate = pd.max_kill_rate`` — the asymptotic kill ceiling (E -> inf).

The Hill probability is then modulated multiplicatively::

    p_kill = p_hill * sensitivity * R_hypoxia * R_quiescence

* ``sensitivity`` — per-cell intrinsic sensitivity in [0,1] (MGMT status,
  resistance clones, ...). A fully resistant cell (0) is never killed.
* ``R_hypoxia   = pd.hypoxia_resistance``   if the cell is HYPOXIC else 1.0.
* ``R_quiescence = pd.quiescent_resistance`` if the cell is QUIESCENT else 1.0.

Both resistance factors are < 1, so hypoxic / quiescent cells survive more often
than viable cells at identical exposure (the classic drug-resistance of poorly
perfused, non-cycling tumor regions).

``p_kill`` is clipped to [0, 1] and a Bernoulli draw with the supplied RNG
decides which non-dead cells die. Given the same RNG the result is fully
deterministic.
"""
from __future__ import annotations

from typing import Optional

import numpy as np

from ..core import (
    CellState,
    DrugParams,
    Microenvironment,
    PDConfig,
    TransportResult,
)


def apply_drug_kill(
    micro: Microenvironment,
    transport: TransportResult,
    pd: PDConfig,
    drug: DrugParams,
    rng: Optional[np.random.Generator] = None,
) -> np.ndarray:
    """Apply pharmacodynamic drug kill to the cells of ``micro`` in place.

    For each cell, the local drug AUC is sampled from ``transport.auc`` at the
    cell's grid index (via ``micro.domain.to_index``), converted to a kill
    probability through a Hill dose-response on the time-averaged effective
    concentration, and modulated by per-cell sensitivity plus hypoxia /
    quiescence resistance (see the module docstring for the exact mapping and the
    E50 derivation). Killed non-dead cells are set to ``CellState.DEAD`` and each
    cell's accumulated AUC exposure is written into ``micro.cells.exposure``.

    Parameters
    ----------
    micro : Microenvironment
        The microenvironment whose ``cells`` are killed in place.
    transport : TransportResult
        Drug-transport solution; ``transport.auc`` (ny, nx) supplies per-cell
        exposure and ``transport.times`` supplies the exposure-window duration.
    pd : PDConfig
        Pharmacodynamic knobs: ``max_kill_rate``, ``hypoxia_resistance``,
        ``quiescent_resistance``.
    drug : DrugParams
        Drug parameters; ``drug.ic50`` sets E50 and ``drug.hill`` sets the Hill
        coefficient.
    rng : np.random.Generator, optional
        RNG for the Bernoulli kill draws. If ``None`` a fresh
        ``np.random.default_rng()`` is used. Pass a seeded generator for
        deterministic results.

    Returns
    -------
    np.ndarray
        Per-cell kill-probability array (length ``micro.cells.n``), for
        diagnostics. Already-dead cells are reported with probability 0.0.
    """
    if rng is None:
        rng = np.random.default_rng()

    cells = micro.cells
    n = cells.n

    # --- sample local AUC exposure at each cell's grid index -----------------
    col, row = micro.domain.to_index(cells.x, cells.y)
    auc_per_cell = np.asarray(transport.auc, dtype=np.float64)[row, col]  # (n,) uM*s

    # store accumulated exposure (AUC) back onto the cells
    cells.exposure = auc_per_cell.copy()

    # --- reduce AUC (uM*s) to a time-averaged effective concentration (uM) ---
    times = np.asarray(transport.times, dtype=np.float64)
    duration = float(times[-1]) if times.size and times[-1] > 0.0 else 1.0
    E = auc_per_cell / duration  # (n,) uM

    # --- Hill dose-response -> kill probability ------------------------------
    h = float(drug.hill)
    E50 = float(drug.ic50)  # E50 == IC50 because E is a concentration (see docstring)
    E_h = np.power(np.maximum(E, 0.0), h)
    E50_h = E50 ** h
    p_hill = pd.max_kill_rate * E_h / (E50_h + E_h + 1e-300)

    # --- modulation: sensitivity + hypoxia / quiescence resistance -----------
    resistance = np.ones(n, dtype=np.float64)
    resistance[cells.state == CellState.HYPOXIC] = pd.hypoxia_resistance
    resistance[cells.state == CellState.QUIESCENT] = pd.quiescent_resistance

    p_kill = p_hill * cells.sensitivity * resistance
    p_kill = np.clip(p_kill, 0.0, 1.0)

    # dead cells cannot be killed again; report probability 0 for them
    alive = cells.state != CellState.DEAD
    p_kill = np.where(alive, p_kill, 0.0)

    # --- Bernoulli kill draws (deterministic given rng) ----------------------
    draws = rng.random(n)
    killed = alive & (draws < p_kill)
    cells.state[killed] = np.int8(CellState.DEAD)

    return p_kill


if __name__ == "__main__":
    from ..core import (
        CellPopulation,
        CellType,
        Domain,
        TissueFields,
    )

    rng = np.random.default_rng(12345)

    # ---- build a domain with a spatially varying AUC field ------------------
    # AUC increases smoothly along +x, so a cell's column controls its exposure.
    nx, ny, dx = 40, 40, 10.0
    domain = Domain(nx=nx, ny=ny, dx=dx)
    xs = np.arange(nx)
    # AUC ranges ~0 (left) to a high value (right), in uM*s.
    duration = 1000.0  # s
    ic50 = 5.0
    # peak time-averaged conc reaches ~4x IC50 on the right edge
    auc_row = (xs / (nx - 1)) * (4.0 * ic50) * duration
    auc = np.tile(auc_row, (ny, 1))  # (ny, nx)
    transport = TransportResult(
        auc=auc,
        final_concentration=np.zeros((ny, nx)),
        times=np.array([0.0, duration]),
        mean_concentration_history=np.array([0.0, 0.0]),
    )

    drug = DrugParams(name="test-drug", ic50=ic50, hill=1.5)
    pd_cfg = PDConfig(max_kill_rate=1.0, hypoxia_resistance=0.35, quiescent_resistance=0.5)

    def make_cells(n_per, cols, sensitivity, state):
        """n_per cells at each given column, spread over rows, given sens/state."""
        x_list, y_list = [], []
        for c in cols:
            xc = (c + 0.5) * dx
            for r in range(n_per):
                x_list.append(xc)
                y_list.append((r % ny + 0.5) * dx)
        n = len(x_list)
        return CellPopulation(
            x=np.array(x_list),
            y=np.array(y_list),
            cell_type=np.full(n, CellType.TUMOR),
            state=np.full(n, state),
            sensitivity=np.full(n, sensitivity),
        )

    def kill_fraction(cells_before, drug, pd_cfg, seed=7):
        micro = Microenvironment(
            domain=domain,
            fields=TissueFields(ecm_density=domain.zeros(), vascular_density=domain.zeros()),
            cells=cells_before.copy(),
        )
        apply_drug_kill(micro, transport, pd_cfg, drug, rng=np.random.default_rng(seed))
        dead = int((micro.cells.state == CellState.DEAD).sum())
        return dead / micro.cells.n

    n_per = 400

    # (1) higher AUC => higher kill fraction (same sensitivity, viable) --------
    low_auc = make_cells(n_per, cols=[2, 3, 4], sensitivity=0.8, state=CellState.VIABLE)
    high_auc = make_cells(n_per, cols=[35, 36, 37], sensitivity=0.8, state=CellState.VIABLE)
    kf_low = kill_fraction(low_auc, drug, pd_cfg)
    kf_high = kill_fraction(high_auc, drug, pd_cfg)
    print(f"kill_fraction  low-AUC={kf_low:.3f}   high-AUC={kf_high:.3f}")
    assert kf_high > kf_low, "higher AUC must give higher kill fraction"

    # (2) higher sensitivity => higher kill fraction (same AUC, viable) --------
    lo_sens = make_cells(n_per, cols=[20, 21, 22], sensitivity=0.2, state=CellState.VIABLE)
    hi_sens = make_cells(n_per, cols=[20, 21, 22], sensitivity=0.9, state=CellState.VIABLE)
    kf_los = kill_fraction(lo_sens, drug, pd_cfg)
    kf_his = kill_fraction(hi_sens, drug, pd_cfg)
    print(f"kill_fraction  low-sens={kf_los:.3f}  high-sens={kf_his:.3f}")
    assert kf_his > kf_los, "higher sensitivity must give higher kill fraction"

    # (3) hypoxic cells survive more than viable at EQUAL exposure -------------
    viable = make_cells(n_per, cols=[25, 26, 27], sensitivity=0.8, state=CellState.VIABLE)
    hypoxic = make_cells(n_per, cols=[25, 26, 27], sensitivity=0.8, state=CellState.HYPOXIC)
    kf_viable = kill_fraction(viable, drug, pd_cfg)
    kf_hypoxic = kill_fraction(hypoxic, drug, pd_cfg)
    print(f"kill_fraction  viable={kf_viable:.3f}   hypoxic={kf_hypoxic:.3f}")
    assert kf_hypoxic < kf_viable, "hypoxic cells must survive more than viable"

    # (4) exposure is recorded on the cells ------------------------------------
    micro = Microenvironment(
        domain=domain,
        fields=TissueFields(ecm_density=domain.zeros(), vascular_density=domain.zeros()),
        cells=high_auc.copy(),
    )
    p = apply_drug_kill(micro, transport, pd_cfg, drug, rng=np.random.default_rng(1))
    assert micro.cells.exposure.shape == (micro.cells.n,)
    assert np.all(micro.cells.exposure > 0), "high-AUC cells must record positive exposure"
    assert p.shape == (micro.cells.n,) and np.all((p >= 0) & (p <= 1))

    # (5) determinism: same rng seed => identical outcome ----------------------
    kf_a = kill_fraction(high_auc, drug, pd_cfg, seed=99)
    kf_b = kill_fraction(high_auc, drug, pd_cfg, seed=99)
    assert kf_a == kf_b, "must be deterministic given the rng seed"

    print("All pkpd self-tests passed.")
