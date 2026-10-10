"""Is the per-spot NNLS degeneracy a property of the FFPE sections, or of the
estimator itself? Measured on the PUBLISHED fresh-frozen sections for a baseline,
and on the signature matrices' own conditioning.
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

from tmesim.reconstruct.deconvolve import build_signatures, deconvolve_sample  # noqa: E402
from tmesim.data.real.scref_gbm import reference_gene_symbols  # noqa: E402
import io_loaders as L  # noqa: E402
from run_test import dedup_first, signatures_sparse  # noqa: E402


def report(name, fr, names):
    ent = float(np.mean(-(fr * np.log(np.clip(fr, 1e-12, 1))).sum(axis=1)))
    print(f"  {name:>18s}  >0.99 {100*np.mean(fr.max(1) > .99):5.1f}%   "
          f">0.90 {100*np.mean(fr.max(1) > .90):5.1f}%   entropy {ent:.3f}   "
          + " ".join(f"{n}={100*v:5.2f}" for n, v in zip(names, fr.mean(0))))


print("=" * 78)
print("BASELINE -- per-spot degeneracy on the PUBLISHED fresh-frozen sections")
print("  (uniform entropy = 1.386)")
for sid in ("10x_parent_gbm", "GBM1", "GBM2", "GBM3"):
    try:
        _c, fr, nm = deconvolve_sample(sid)
        report(sid, fr, nm)
    except Exception as exc:  # noqa: BLE001
        print(f"  {sid}: {type(exc).__name__}: {exc}")

print("=" * 78)
print("SIGNATURE CONDITIONING -- how separable are the four columns?")
sc, sc_genes, _ = L.read_10x("/tmp/vmp/ex_scrna")
labels = list(np.load("/tmp/vmp/out/scrna_labels.npy", allow_pickle=True))
vx, vgenes, _ = L.read_10x(
    "/tmp/vmp/ex_visium/GSM9497624_TWBK-ST073021-2_filtered_feature_bc_matrix")

sc_dedup, _ = dedup_first(sc_genes)
ref_sets = {"own_b189": set(sc_dedup), "darmanis": set(reference_gene_symbols())}
for refname, refsyms in ref_sets.items():
    seen, shared = set(), []
    for g in vgenes:
        if g in refsyms and g not in seen:
            seen.add(g)
            shared.append(g)
    if refname == "own_b189":
        sig, nm = signatures_sparse(sc, sc_genes, labels, shared)
    else:
        sig, nm = build_signatures(shared)
    norm = sig / np.maximum(np.linalg.norm(sig, axis=0, keepdims=True), 1e-300)
    cos = norm.T @ norm
    print(f"\n  {refname}  ({len(shared)} genes)  cond = "
          f"{np.linalg.cond(sig):.3e}")
    print("        " + "".join(f"{n[:5]:>8s}" for n in nm))
    for i, n in enumerate(nm):
        print(f"  {n[:6]:>6s}" + "".join(f"{cos[i, j]:8.4f}" for j in range(len(nm))))
    print("  column sums: " + " ".join(f"{n}={s:.4f}" for n, s in zip(nm, sig.sum(0))))
