"""Print every number the write-up is allowed to quote, straight from the data.

The report is written by hand, but no figure in it is typed from memory: this
prints the exact values out of ``sweep_data.json``, ``fit_results.json``,
``analyse_results.json`` and ``coarse_baseline.json``, labelled with the name
the report uses for each one. A number that does not appear here does not go in
the report.
"""
from __future__ import annotations

import json
import os
from typing import Dict, Optional

_HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name: str) -> Optional[Dict]:
    p = os.path.join(_HERE, name)
    if not os.path.exists(p):
        print(f"  (missing: {name})")
        return None
    with open(p) as fh:
        return json.load(fh)


def sec(title: str) -> None:
    print(f"\n{'=' * 78}\n{title}\n{'=' * 78}")


def main() -> int:
    sw = _load("sweep_data.json")
    an = _load("analyse_results.json")
    ft = _load("fit_results.json")
    cb = _load("coarse_baseline.json")

    # ------------------------------------------------------------ the dataset
    if sw:
        sec("1. THE DATASET (sweep_data.json)")
        rows = sw["rows"]
        ok = [r for r in rows if r.get("status") == "ok"]
        bad = [r for r in rows if r.get("status") != "ok"]
        secs = sorted({r["section"] for r in ok})
        dons = sorted({r["donor"] for r in ok})
        valid_secs = sorted({r["section"] for r in ok
                             if r["spatial_spread_valid"]})
        valid_dons = sorted({r["donor"] for r in ok
                             if r["spatial_spread_valid"]})
        print(f"alpha                       {sw.get('alpha')!r}")
        print(f"alpha_end                   {sw.get('alpha_end')!r}")
        print(f"dx_um                       {sw.get('dx_um')}")
        print(f"rows total                  {len(rows)}")
        print(f"rows ok                     {len(ok)}")
        print(f"rows failed                 {len(bad)}")
        print(f"sections attempted          {len(sw.get('sections') or {})}")
        print(f"sections with ok rows       {len(secs)}")
        print(f"donor groups                {len(dons)}")
        print(f"sections spatial_spread OK  {len(valid_secs)}  {valid_secs}")
        print(f"donors spatial_spread OK    {len(valid_dons)}  {valid_dons}")
        print(f"compounds                   {sw.get('compounds')}")
        print(f"pk_shapes                   {sw.get('pk_shapes')}")
        print(f"sweep elapsed_s             {sw.get('elapsed_s')}")
        if sw.get("elapsed_s"):
            print(f"sweep elapsed_min           "
                  f"{float(sw['elapsed_s']) / 60:.1f}")
        for r in bad[:8]:
            print(f"  FAILED {r['section']}/{r.get('drug')}: "
                  f"{r.get('error_type')}: {str(r.get('error'))[:110]}")
        fails = {}
        for r in bad:
            fails.setdefault(r["section"], 0)
            fails[r["section"]] += 1
        if fails:
            print(f"failures by section         {fails}")

    # -------------------------------------------------------- axis + spread
    if an:
        sec("2. AXIS SENSITIVITY — which inputs move a target at all")
        for axis, rep in an["axis_sensitivity"].items():
            print(f"\n  axis {axis}  "
                  f"(n comparison groups {rep['n_comparison_groups']})")
            for t, v in rep["max_fold_change_per_target"].items():
                print(f"    {t:26s} max fold change {v['max']!r}  "
                      f"median {v['median']!r}  n={v['n']}")

        sec("3. REACH LADDER — is it even a target?")
        for k, v in an["reach_ladder"].items():
            print(f"  {k:46s} {v!r}")

        sec("4. ACROSS-DONOR SPREAD of each target (shipped config, TMZ)")
        for t, v in an["across_donor_spread"].items():
            print(f"\n  {t}")
            print(f"    config          {v['config']}")
            print(f"    n_donors        {v['n_donors']}")
            print(f"    min / max       {v['min']!r} / {v['max']!r}")
            print(f"    fold_spread     {v['fold_spread']!r}")
            print(f"    abs_range       {v['abs_range']!r}")
            print(f"    valid-only med  {v['valid_only_per_donor_median']}")
            print(f"    valid-only fold {v['valid_only_fold_spread']!r}")

        sec("5. SOLVER COST at dx=30 in this sweep")
        for k, v in an["solver_cost_dx30"].items():
            print(f"  {k:28s} {v!r}")

    # ------------------------------------------------------------- the models
    if ft:
        sec("6. HELD-OUT ERROR, by split (fit_results.json)")
        print(f"n_rows_ok {ft['n_rows_ok']}  donors {ft['n_donors']}  "
              f"sections {ft['n_sections']}  compounds {ft['n_compounds']}")
        print(f"n_features {len(ft['feature_names'])}")
        for e in ft["evals"]:
            print(f"\n  {e['target']}  split={e['split']}  "
                  f"subset={e['subset']}  n={e['n']}  folds={e['n_folds']}")
            print(f"    unit: {e['unit']}")
            print(f"    y range            {e['y_min']!r} .. {e['y_max']!r}  "
                  f"median {e['y_median']!r}")
            print(f"    across-group med   {e['across_group_median_min']!r} .."
                  f" {e['across_group_median_max']!r}")
            print(f"    across-group fold  "
                  f"{e['across_group_fold_spread']!r}   abs range "
                  f"{e['across_group_abs_range']!r}")
            for m in ("mean", "analytic", "ridge", "hgb"):
                r = e["models"].get(m)
                if not r or not r.get("n_scored"):
                    continue
                print(f"    {m:9s} MAE={r['mae']!r} medAE={r['median_ae']!r} "
                      f"p90AE={r['p90_ae']!r} maxAE={r['max_ae']!r}")
                print(f"    {'':9s} foldErr med={r['median_fold_error']!r} "
                      f"p90={r['p90_fold_error']!r} R2={r['r2']!r} "
                      f"skill_vs_mean={r.get('skill_vs_mean')!r} "
                      f"skill_vs_analytic={r.get('skill_vs_analytic')!r}")

        sec("7. FEATURE-GROUP ABLATION (donor split)")
        for t, rec in (ft.get("ablation_donor_split") or {}).items():
            print(f"\n  {t}")
            for cname, r in rec.items():
                if "error" in r:
                    print(f"    {cname:14s} ERROR {r['error'][:80]}")
                    continue
                print(f"    {cname:14s} nfeat={r['n_features']:3d} "
                      f"ridge MAE={(r.get('ridge') or {}).get('mae')!r} "
                      f"skill={(r.get('ridge') or {}).get('skill_vs_mean')!r}")
                print(f"    {'':14s} {'':9s} "
                      f"hgb   MAE={(r.get('hgb') or {}).get('mae')!r} "
                      f"skill={(r.get('hgb') or {}).get('skill_vs_mean')!r}")

        sec("8. COST OF A FINAL MODEL")
        for t, rec in (ft.get("final_model_cost") or {}).items():
            print(f"\n  {t}  (n_train {rec['n_train']})")
            for m in ("ridge", "hgb"):
                c = rec.get(m) or {}
                print(f"    {m:6s} train={c.get('train_s')!r}s  "
                      f"pickle={c.get('pickle_bytes')!r} B  "
                      f"predict1 median={c.get('single_row_predict_ms_median')!r} ms  "
                      f"p90={c.get('single_row_predict_ms_p90')!r} ms")

    # --------------------------------------------------------- coarse vs fine
    if cb:
        sec("9. COARSE SOLVE vs FINE SOLVE (the non-ML baseline)")
        print(f"preview dx {cb['preview_dx_um']}  shipped dx "
              f"{cb['shipped_dx_um']}")
        for t, s in (cb.get("coarse_vs_fine") or {}).items():
            print(f"\n  {t}")
            for k, v in s.items():
                print(f"    {k:18s} {v!r}")
        print(f"\n  timing {json.dumps(cb.get('timing') or {}, indent=2)}")
        print("\n  per-pair detail:")
        for pr in cb["pairs"]:
            c, f = pr.get("coarse") or {}, pr.get("fine") or {}
            if c.get("status") != "ok" or f.get("status") != "ok":
                print(f"    {pr['section'][:26]:26s} {pr['drug'][:13]:13s} "
                      f"coarse={c.get('status')} fine={f.get('status')}")
                continue
            print(f"    {pr['section'][:26]:26s} {pr['drug'][:13]:13s} "
                  f"E_p90/p10 {c.get('E_p90_over_p10')!r:22.22s} -> "
                  f"{f.get('E_p90_over_p10')!r:22.22s}  "
                  f"L {c.get('penetration_length_um')!r:20.20s} -> "
                  f"{f.get('penetration_length_um')!r:20.20s}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
