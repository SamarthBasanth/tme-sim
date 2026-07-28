"""Export a REAL reconstructed microenvironment to compact JSON for the viewer.

Subsamples cells for browser performance but keeps the real spatial structure and
the real field grids. Everything here is measured data or a field computed from it.
"""
import warnings; warnings.filterwarnings("ignore")
import json, sys
import numpy as np

from tmesim.reconstruct.real_pipeline import reconstruct_real, reconstruction_report
from tmesim.core import CellType

DX = 20.0            # coarse enough to export quickly, fine enough to look right
MAX_CELLS = 18000    # browser-friendly

print("reconstructing real Xenium microenvironment ...", flush=True)
micro = reconstruct_real(source="xenium", dx_um=DX)
rep = reconstruction_report(micro)

cells = micro.cells
n = cells.n
rng = np.random.default_rng(0)
idx = np.arange(n) if n <= MAX_CELLS else rng.choice(n, MAX_CELLS, replace=False)
idx = np.sort(idx)

# round positions to 1 um to shrink the payload
x = np.round(cells.x[idx], 1).tolist()
y = np.round(cells.y[idx], 1).tolist()
t = cells.cell_type[idx].astype(int).tolist()

def grid(field, decimals=3):
    if field is None:
        return None
    a = np.asarray(field, dtype=float)
    return [[round(float(v), decimals) for v in row] for row in a]

f = micro.fields
payload = {
    "meta": {
        "patient_id": micro.patient_id,
        "dx_um": DX,
        "grid": list(micro.domain.shape),
        "extent_um": [round(v, 1) for v in micro.domain.extent_um],
        "n_cells_total": int(n),
        "n_cells_shown": int(len(idx)),
        "source": "10x Xenium human GBM, 319-gene panel (CC BY 4.0)",
    },
    "report": {k: (round(v, 4) if isinstance(v, float) else v)
               for k, v in rep.items() if isinstance(v, (int, float, str))},
    "celltypes": {int(c.value): c.name for c in CellType},
    "cells": {"x": x, "y": y, "t": t},
    "fields": {
        "vascular_density": grid(f.vascular_density),
        "ecm_density": grid(f.ecm_density),
        "oxygen": grid(f.oxygen),
    },
}

counts = {CellType(v).name: int((cells.cell_type == v).sum()) for v in set(cells.cell_type.tolist())}
payload["counts"] = counts

out = "results/viewer_data.json"
with open(out, "w") as fh:
    json.dump(payload, fh, separators=(",", ":"))
import os
print(f"wrote {out}  ({os.path.getsize(out)/1e6:.1f} MB)")
print("grid:", micro.domain.shape, "extent um:", payload["meta"]["extent_um"])
print("cells shown:", len(idx), "of", n)
print("counts:", counts)
