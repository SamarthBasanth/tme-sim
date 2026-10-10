"""Tabulate stage_a.json + stage_b.json against the pre-registered criteria."""
from __future__ import annotations

import json
import os
import sys

OUT = sys.argv[1] if len(sys.argv) > 1 else "/tmp/vmp/out"
A = json.load(open(os.path.join(OUT, "stage_a.json")))
B = json.load(open(os.path.join(OUT, "stage_b.json")))
TAU_R, TAU_C = 0.1778, 0.4884
ORDER = ["TUMOR", "STROMAL", "IMMUNE", "ENDOTHELIAL"]

print("### Footprint reproduction vs results/data_scouting_threads_closed §1.4")
print("| pair | spots | Visium mm2 | Xenium cells | Xenium mm2 | area x | cells/mm2 |")
print("|---|---:|---:|---:|---:|---:|---:|")
for pid, r in sorted(B["pairs"].items()):
    s = A["sections"][pid]
    print(f"| {pid} {r['xenium_stem']} | {s['n_spots_in_tissue_used']:,} | "
          f"{s['visium_area_mm2']:.2f} | {r['n_cells']:,} | {r['xenium_area_mm2']:.2f} | "
          f"**{r['area_ratio_xen_over_vis']:.2f}** | "
          f"{r['xenium_cell_density_per_mm2']:,.0f} |")

print("\n### Compositions, percent")
print("| pair | side | TUMOR | STROMAL | IMMUNE | ENDOTHELIAL |")
print("|---|---|---:|---:|---:|---:|")
for pid, r in sorted(B["pairs"].items()):
    rows = [("Xenium whole field (truth)", r["xenium_whole_field_pct"])]
    for k, lbl in (("own_b189", "Visium NNLS, B189 own ref"),
                   ("darmanis", "Visium NNLS, Darmanis ref")):
        if k in r:
            rows.append((lbl, r[k]["visium_composition_pct"]))
    for lbl, d in rows:
        print(f"| {pid} | {lbl} | " + " | ".join(f"{d[n]:.2f}" for n in ORDER) + " |")

print("\n### Region-mismatch band and the criteria")
print("| pair | area x | windows | L1 band 50 | L1 band 95 | ref | L1_full | L1_best "
      "| mean|d| best (pp) | C1 fail? | C2 met? | C3 |")
print("|---|---:|---:|---:|---:|---|---:|---:|---:|---|---|---|")
for pid, r in sorted(B["pairs"].items()):
    for k in ("own_b189", "darmanis"):
        if k not in r:
            continue
        x = r[k]
        print(f"| {pid} | {r['area_ratio_xen_over_vis']:.2f} | {r['windows_accepted']} "
              f"| {r['l1_wf_50']:.4f} | {r['l1_wf_95']:.4f} | {k} "
              f"| {x['l1_full']:.4f} | {x['l1_best']:.4f} "
              f"| {x['mean_abs_diff_pp_best']:.2f} "
              f"| {'YES' if x['C1_fail_worse_than_current'] else 'no'} "
              f"| {'yes' if x['C2_meets_retracted_era'] else 'NO'} "
              f"| {x['l1_full_vs_band95']} |")

print("\n### C4 non-tumour 3-vector")
print("| pair | band 50 | band 95 | ref | L1nt_full | L1nt_best |")
print("|---|---:|---:|---|---:|---:|")
for pid, r in sorted(B["pairs"].items()):
    for k in ("own_b189", "darmanis"):
        if k not in r:
            continue
        print(f"| {pid} | {r['l1nt_wf_50']:.4f} | {r['l1nt_wf_95']:.4f} | {k} "
              f"| {r[k]['l1nt_full']:.4f} | {r[k]['l1nt_best']:.4f} |")

print("\n### C5 Spearman rho (reported, not a gate)")
for pid, r in sorted(B["pairs"].items()):
    bits = [f"{k}={r[k]['spearman_rho_full']:+.2f}" for k in ("own_b189", "darmanis")
            if k in r]
    print(f"- {pid}: " + "  ".join(bits))

print("\n### C6 own-reference vs Darmanis, same pair, same windows")
print("| pair | L1_best own | L1_best Darmanis | own better? | predicted |")
print("|---|---:|---:|---|---|")
for pid, r in sorted(B["pairs"].items()):
    if "own_b189" in r and "darmanis" in r:
        o, d = r["own_b189"]["l1_best"], r["darmanis"]["l1_best"]
        print(f"| {pid} | {o:.4f} | {d:.4f} | "
              f"{'YES' if o < d else 'NO'} | own < Darmanis |")

print("\n### C7 unmixing-only arm (Xenium pseudo-spots, same B189 signature)")
for pid, r in sorted(B["pairs"].items()):
    c7 = r.get("C7_pseudospots")
    if not c7:
        continue
    print(f"- **{pid}**  n={c7['n_pseudospots']:,} pseudo-spots, {c7['n_genes']} genes; "
          f"mean per-spot L1 {c7['mean_l1_per_spot']:.4f} "
          f"(median {c7['median_l1_per_spot']:.4f}); "
          f"pooled L1 {c7['pooled_l1']:.4f} = {c7['pooled_mean_abs_diff_pp']:.2f} pp")
    print(f"    - truth  " + "  ".join(f"{n}={c7['pooled_truth_pct'][n]:5.2f}%"
                                       for n in ORDER))
    print(f"    - pred   " + "  ".join(f"{n}={c7['pooled_pred_pct'][n]:5.2f}%"
                                       for n in ORDER))

print("\n### Provenance")
print(f"- prereg sha256 recorded in stage_b: {B.get('prereg_sha256')}")
print(f"- scRNA: {A['scrna_n_cells']:,} cells x {A['scrna_n_genes']:,} genes; "
      f"classes {A['scrna_class_counts']}")
print(f"- marker map {A['marker_map_panel']!r}, coverage {A['marker_map_coverage']}, "
      f"missing {A['markers_missing_from_panel']}")
print(f"- sparse-vs-repo signature max|diff| = "
      f"{A['signature_equivalence_max_abs_diff']:.3e}")

print("\n### Within-donor reproducibility, B189 P1 vs P2 (no Xenium truth needed "
      "for the Visium rows)")
print("| quantity | L1 | mean abs diff (pp) |")
print("|---|---:|---:|")
import numpy as _np
p1, p2 = B["pairs"]["P1"], B["pairs"]["P2"]
rows = [("Xenium truth, P1 vs P2",
         [p1["xenium_whole_field_pct"][n] for n in ORDER],
         [p2["xenium_whole_field_pct"][n] for n in ORDER])]
for k, lbl in (("own_b189", "Visium NNLS B189 own ref, P1 vs P2"),
               ("darmanis", "Visium NNLS Darmanis ref, P1 vs P2")):
    rows.append((lbl,
                 [p1[k]["visium_composition_pct"][n] for n in ORDER],
                 [p2[k]["visium_composition_pct"][n] for n in ORDER]))
for lbl, a, b in rows:
    d = float(_np.abs(_np.array(a) - _np.array(b)).sum())
    print(f"| {lbl} | {d/100:.4f} | {d/4:.2f} |")
