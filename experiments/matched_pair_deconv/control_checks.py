"""Controls run BEFORE the stage-A numbers are interpreted.

1. Does the repo's own ``deconvolve_sample`` still reproduce the Visium vector
   published in ``results/revalidation.md`` (82.66 / 1.04 / 6.53 / 9.77)?
2. Does MY call path reproduce the same thing on the same section with the same
   Darmanis reference? If yes, the pipeline is not the explanation for the
   FFPE-section result and the degeneracy is a property of the data.
3. How degenerate are the per-spot NNLS solutions on each FFPE section?
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

from tmesim.reconstruct.deconvolve import (  # noqa: E402
    build_signatures, deconvolve_sample, deconvolve_spot, iter_spot_vectors,
)
from tmesim.data.real.scref_gbm import reference_gene_symbols  # noqa: E402
import io_loaders as L  # noqa: E402
from run_test import PAIRS, dedup_first  # noqa: E402

print("=" * 72)
print("CONTROL 1 -- repo's own deconvolve_sample on the published section")
coords, fr, names = deconvolve_sample("10x_parent_gbm")
print(f"  spots={fr.shape[0]} types={names}")
m = fr.mean(axis=0)
for n, v in zip(names, m):
    print(f"  {n:12s} {100*v:6.2f}%")
print("  published in results/revalidation.md: TUMOR 82.66 STROMAL 1.04 "
      "IMMUNE 6.53 ENDOTHELIAL 9.77")

print("=" * 72)
print("CONTROL 3 -- per-spot degeneracy on the FFPE matched-pair sections")
ref_syms = set(reference_gene_symbols())
for pid, vstem, _xg, _xs, _d, _a in PAIRS:
    vdir_m = os.path.join("/tmp/vmp/ex_visium", f"{vstem}_filtered_feature_bc_matrix")
    vdir_s = os.path.join("/tmp/vmp/ex_visium", f"{vstem}_spatial")
    vx, vgenes, vbcs = L.read_10x(vdir_m)
    coords_map, _um = L.read_visium_spatial(vdir_s)
    keep = np.asarray([b in coords_map for b in vbcs], dtype=bool)
    vx = vx[keep].tocsr()
    seen, shared, cols = set(), [], []
    for j, g in enumerate(vgenes):
        if g in ref_syms and g not in seen:
            seen.add(g)
            shared.append(g)
            cols.append(j)
    sig, nm = build_signatures(shared)
    sub = vx[:, np.asarray(cols, dtype=int)].tocsr()
    frs = np.zeros((sub.shape[0], sig.shape[1]))
    for i, spot in enumerate(iter_spot_vectors(sub, len(shared))):
        frs[i] = deconvolve_spot(spot, sig)
    argmax_share = np.bincount(frs.argmax(axis=1), minlength=len(nm)) / frs.shape[0]
    print(f"  {pid} {vstem}")
    print(f"    mean total counts/spot = {np.asarray(vx.sum(axis=1)).mean():,.0f}")
    print(f"    spots where max fraction > 0.99 : "
          f"{100*np.mean(frs.max(axis=1) > 0.99):5.1f}%")
    print(f"    spots where max fraction > 0.90 : "
          f"{100*np.mean(frs.max(axis=1) > 0.90):5.1f}%")
    print(f"    mean entropy of the 4-vector    : "
          f"{float(np.mean(-(frs*np.log(np.clip(frs,1e-12,1))).sum(axis=1))):.3f}"
          f"  (uniform = {np.log(len(nm)):.3f})")
    print("    per-spot argmax share: "
          + "  ".join(f"{n}={100*s:5.1f}%" for n, s in zip(nm, argmax_share)))
    print("    mean fraction       : "
          + "  ".join(f"{n}={100*v:5.2f}%" for n, v in zip(nm, frs.mean(axis=0))))
