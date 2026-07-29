"""Sensitivity-band sweep for ``OxygenConfig.consumption`` x ``sink_form``.

Regenerates the grid reported in ``results/oxygen_sensitivity_band.md``. Builds
the flagship real reconstruction ONCE
(``reconstruct_real('xenium', dx_um=12.0)``, the grid-converged production
spacing used throughout ``results/``) and re-solves
:func:`tmesim.physics.oxygen.solve_oxygen` many times against it, for:

1. The literature-defensible ``consumption`` grid derived in
   ``tmesim/physics/oxygen.py``'s module docstring, bracketed on one side by
   Watabe et al. 2017 (rat C6 glioma 15O-PET TMRO2, PMID 28210996 -- implies
   ~12-25x HIGHER than the shipping default) and on the other by the
   Diepart-derived first-order rate constant cited in Milotti et al. 2017
   (DOI 10.1038/srep39762 -- implies ~4-100x LOWER), crossed with both sink
   forms ("linear", the shipping default, and "michaelis_menten" at a fixed,
   literature-derived half-saturation constant).
2. A secondary sweep of the Michaelis-Menten half-saturation constant ``K`` at
   the current default ``consumption``, to show that axis's own sensitivity
   and cross-check this production Newton-Raphson solver against the
   standalone Picard prototype in
   ``results/sweep/oxygen_functional_form_prototype_2026-07-29T05-54-00Z.md``.

Run with the project's own interpreter: ``./.venv/bin/python
scripts/oxygen_sensitivity_band.py``. Takes a few minutes (dominated by the
Michaelis-Menten Newton solves at ~286k unknowns/grid point); everything here
is deterministic (no seeds).
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from tmesim.core import CellType, OxygenConfig
from tmesim.physics.oxygen import solve_oxygen
from tmesim.reconstruct.real_pipeline import reconstruct_real

DEFAULT_CONSUMPTION = 0.02
HYPOXIA_THRESHOLD = OxygenConfig().hypoxia_threshold

# Literature-defensible band (see tmesim/physics/oxygen.py module docstring):
#   Diepart/Milotti route: 4x-100x LOWER than default -> [default/100, default/4]
#   Watabe route:          12x-25x HIGHER than default -> [default*12, default*25]
CONSUMPTION_GRID = [
    0.0002,  # = default / 100  (Diepart/Milotti lower bound)
    0.0005,
    0.001,
    0.002,
    0.005,   # = default / 4    (Diepart/Milotti upper bound)
    0.01,
    0.02,    # = current shipping default (UNCHANGED)
    0.05,
    0.1,
    0.24,    # = default * 12   (Watabe lower bound)
    0.5,     # = default * 25   (Watabe upper bound)
]

# Grimes et al. 2014 (PMID 26064525) base-case Km = 1 mmHg, expressed against
# this project's own O=1 anchor (oxygen.REFERENCE_PO2_MMHG = 35 mmHg,
# Erecinska & Silver 2001, PMID 11718758): K = 1 / 35.
MM_K_DEFAULT = 1.0 / 35.0

K_GRID = [0.01, 0.025, 0.05, 0.10, 0.25, 0.50]  # secondary K-sensitivity sweep


def hypoxic_fraction(oxygen: np.ndarray, row: np.ndarray, col: np.ndarray) -> float:
    local = oxygen[row, col]
    return float((local < HYPOXIA_THRESHOLD).mean())


def main() -> None:
    print("Reconstructing flagship section (xenium, dx_um=12.0)...")
    t0 = time.time()
    micro = reconstruct_real("xenium", dx_um=12.0)
    print(f"  done in {time.time()-t0:.2f}s, grid={micro.domain.shape}, "
          f"n_cells={micro.cells.n}")

    tum_mask = micro.cells.cell_type == int(CellType.TUMOR)
    n_tumor = int(tum_mask.sum())
    print(f"  n_tumor_cells={n_tumor}")
    row_col = micro.domain.to_index(micro.cells.x[tum_mask], micro.cells.y[tum_mask])
    col, row = row_col

    results = []
    print("\n=== LINEAR sink (current default functional form) ===")
    for c in CONSUMPTION_GRID:
        cfg = OxygenConfig(consumption=c, sink_form="linear")
        t0 = time.time()
        oxygen = solve_oxygen(micro, cfg)
        dt = time.time() - t0
        hf = hypoxic_fraction(oxygen, row, col)
        mo_t = float(oxygen[row, col].mean())
        mo_tissue = float(oxygen[micro.domain.tissue_mask].mean())
        print(f"  consumption={c:<8.4g} hypoxic_frac={hf:.4f} "
              f"mean_O2_tumor={mo_t:.4f} mean_O2_tissue={mo_tissue:.4f} time={dt:.2f}s")
        results.append(dict(sink_form="linear", consumption=c, K=None,
                             hypoxic_fraction=hf, mean_o2_tumor=mo_t,
                             mean_o2_tissue=mo_tissue, solve_time_s=dt))

    print(f"\n=== MICHAELIS-MENTEN sink (K={MM_K_DEFAULT:.4f}) ===")
    for c in CONSUMPTION_GRID:
        cfg = OxygenConfig(consumption=c, sink_form="michaelis_menten",
                            mm_half_saturation=MM_K_DEFAULT)
        t0 = time.time()
        oxygen = solve_oxygen(micro, cfg)
        dt = time.time() - t0
        hf = hypoxic_fraction(oxygen, row, col)
        mo_t = float(oxygen[row, col].mean())
        mo_tissue = float(oxygen[micro.domain.tissue_mask].mean())
        print(f"  consumption={c:<8.4g} hypoxic_frac={hf:.4f} "
              f"mean_O2_tumor={mo_t:.4f} mean_O2_tissue={mo_tissue:.4f} time={dt:.2f}s")
        results.append(dict(sink_form="michaelis_menten", consumption=c, K=MM_K_DEFAULT,
                             hypoxic_fraction=hf, mean_o2_tumor=mo_t,
                             mean_o2_tissue=mo_tissue, solve_time_s=dt))

    print("\n=== K-sensitivity at consumption=0.02 (current default) ===")
    k_results = []
    for K in K_GRID:
        cfg = OxygenConfig(consumption=DEFAULT_CONSUMPTION, sink_form="michaelis_menten",
                            mm_half_saturation=K)
        t0 = time.time()
        oxygen = solve_oxygen(micro, cfg)
        dt = time.time() - t0
        hf = hypoxic_fraction(oxygen, row, col)
        mo_t = float(oxygen[row, col].mean())
        mo_tissue = float(oxygen[micro.domain.tissue_mask].mean())
        print(f"  K={K:<6.3f} hypoxic_frac={hf:.4f} mean_O2_tumor={mo_t:.4f} "
              f"mean_O2_tissue={mo_tissue:.4f} time={dt:.2f}s")
        k_results.append(dict(K=K, hypoxic_fraction=hf, mean_o2_tumor=mo_t,
                               mean_o2_tissue=mo_tissue, solve_time_s=dt))

    out_path = Path(__file__).resolve().parent.parent / "results" / "sweep" / \
        "oxygen_sensitivity_band_raw.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(dict(consumption_sweep=results, k_sweep=k_results,
                        n_tumor=n_tumor, n_cells=micro.cells.n,
                        grid=list(micro.domain.shape)), f, indent=2)
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
