"""Scratch probe: can a section run here, how long, and what varies with dose?

Not a result-producing script. It answers two things before the sweep is
designed:

1.  Does an ingested GBM-Space section reconstruct and solve from this worktree,
    and what does it cost?
2.  Is the candidate target ``E_p90/E_p10`` invariant to ``plasma_peak``? If
    the transport chain is linear in the plasma amplitude then the dose axis
    carries NO information for a ratio target, and a dose slider needs an exact
    scaling law rather than a learned model. That is a structural fact about
    the solver and it has to be measured, not assumed.
"""
from __future__ import annotations

import sys

from solver_runner import assert_alpha, make_spec, run_one, section_record

SECTION = sys.argv[1] if len(sys.argv) > 1 else "gbmspace_at10_bra_5_fo_1_1"


def main() -> int:
    from app.server.cache import ReconstructionCache

    assert_alpha("probe start")
    rec = section_record(SECTION)
    print(f"section={SECTION} source={rec.get('source')} "
          f"label={rec.get('label')}")
    for c in (rec.get("caveats") or []):
        print(f"  caveat {c.get('capability')}: {c.get('status')}")

    cache = ReconstructionCache(maxsize=1)
    for peak in (50.0, 100.0):
        spec = make_spec(SECTION, "temozolomide", 30.0, 21600.0,
                         peak, 14400.0, 1, 86400.0)
        out = run_one(SECTION, spec, cache, run_id=f"peak{peak:g}")
        if out["status"] != "ok":
            print(f"  FAILED {out['error_type']}: {out['error']}")
            print(out.get("traceback", "")[-1500:])
            return 1
        print(f"  peak={peak:g}  wall={out['wall_clock_s']:.2f}s "
              f"recon={out['recon_build_s']} hit={out['recon_cache_hit']}")
        print(f"    n_cells={out['n_cells']} grid={out['grid']} "
              f"E_p10={out['E_p10_uM']!r}")
        print(f"    E_p90/p10={out['E_p90_uM'] / out['E_p10_uM']!r}")
        print(f"    L={out['penetration_length_um']!r} "
              f"pen_idx={out['penetration_index']!r}")
        print(f"    median_auc/peak="
              f"{out['median_tumor_auc'] / out['plasma_peak_uM']!r}")
        print(f"    reach_F={out['reach_headline_F']!r} "
              f"null={out['reach_null_reason']!r}")
    assert_alpha("probe end")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
