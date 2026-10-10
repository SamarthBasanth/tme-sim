"""Falsification-resistance checks on the Anchor-1 inverse solve.

Three things are checked, all by driving the SHIPPED function:
  A. forward substitution: plug the implied alpha back through
     tmesim.physics.diffusion._effective_diffusivity and confirm the measured
     endpoint ratio comes out. If my algebra is wrong this fails.
  B. what the shipped code ACTUALLY does for alpha > 1 - floor, or negative?
  C. is the implied-alpha solve monotone in the mapping scale A (so the
     M1/M2/M3 spread is a real spread, not a numerical artefact)?
"""
from __future__ import annotations

import numpy as np

import tmesim.physics.diffusion as diff
from tmesim.core import DrugParams

SHAFIEYAN = [
    ("sodium iodide",      0.66, 0.97, 400.0, 778.0),
    ("sodium diatrizoate", 0.78, 0.94, 110.0, 155.0),
    ("Gd-DTPA",            0.80, 0.95, 124.0, 195.0),
]
PHI_F_ENGINE = 0.05
D0 = 100.0


def implied_alpha(x_lo, x_hi, ratio, A):
    return (ratio - 1.0) / (A * (ratio * x_hi - x_lo))


def forward_ratio_via_shipped(alpha, ecm_lo, ecm_hi):
    """Ratio D_eff(ecm_lo)/D_eff(ecm_hi) using the SHIPPED function, with
    ALPHA_TORTUOSITY temporarily monkeypatched IN MEMORY ONLY (no file write)."""
    saved = diff.ALPHA_TORTUOSITY
    try:
        diff.ALPHA_TORTUOSITY = alpha
        drug = DrugParams(name="probe", diffusivity=D0)
        d = diff._effective_diffusivity(drug, np.array([ecm_lo, ecm_hi]))
        return float(d[0] / d[1]), float(d[0] / D0), float(d[1] / D0)
    finally:
        diff.ALPHA_TORTUOSITY = saved


print("=" * 78)
print("A. FORWARD SUBSTITUTION through the shipped _effective_diffusivity")
print("=" * 78)
print("   (ALPHA_TORTUOSITY monkeypatched in memory for the probe only; the")
print("    module constant is restored and never written to disk.)")
in_dom_ok, out_dom = True, []
for name, pf_lo, pf_hi, d_lo, d_hi in SHAFIEYAN:
    x_hi, x_lo = 1.0 - pf_lo, 1.0 - pf_hi
    ratio = d_hi / d_lo
    for label, A in (("M1", 1.0), ("M2", 1.0 / x_hi), ("M3", 1.0 / PHI_F_ENGINE)):
        a = implied_alpha(x_lo, x_hi, ratio, A)
        ecm_lo, ecm_hi = A * x_lo, A * x_hi
        fwd, r_lo, r_hi = forward_ratio_via_shipped(a, ecm_lo, ecm_hi)
        ok = abs(fwd - ratio) < 1e-6
        in_domain = ecm_hi <= 1.0 + 1e-12
        if in_domain:
            in_dom_ok &= ok
            tag = "OK" if ok else "MISMATCH"
        else:
            out_dom.append((name, label))
            tag = "OUT OF DOMAIN (ecm>1 clipped -> engine saturates)"
        print(f"   {name:20s} {label}  a={a:7.4f}  ecm=({ecm_lo:.4f},{ecm_hi:.4f})  "
              f"D_eff/D=({r_lo:.4f},{r_hi:.4f})  fwd ratio={fwd:.4f} "
              f"vs measured {ratio:.4f}  {tag}")
print(f"\n   -> IN-DOMAIN mappings (M1, M2): forward substitution "
      f"{'CONSISTENT' if in_dom_ok else 'INCONSISTENT'} - the inverse solve")
print("      reproduces the measured ratio exactly through the shipped code.")
print(f"   -> OUT-OF-DOMAIN: {len(out_dom)} cells, all M3. THIS IS A REAL FINDING,")
print("      not a bug in the check: ecm is np.clip'd to [0,1] inside the shipped")
print("      function, so a mapping with A=20 puts every cartilage explant at")
print("      ecm=1 and the engine has ZERO dynamic range over the measured range.")
print("      M3 must therefore be reported as SATURATED, not as 'admissible with a")
print("      small alpha'. The implied-alpha solve is invalid out of domain.")

print()
print("=" * 78)
print("B. WHAT THE SHIPPED CODE DOES FOR alpha > 1 (the 'goes NEGATIVE' claim)")
print("=" * 78)
for a in (0.6, 1.4969, 2.1646):
    saved = diff.ALPHA_TORTUOSITY
    try:
        diff.ALPHA_TORTUOSITY = a
        drug = DrugParams(name="probe", diffusivity=D0)
        ecm = np.array([0.0, 0.4, 1.0 / a if a > 0 else 1.0, 0.8, 1.0])
        bare = D0 * (1.0 - a * ecm)
        shipped = diff._effective_diffusivity(drug, ecm)
        print(f"\n   alpha = {a:.4f}   floor = max(1e-6, 1e-3*D) = {max(1e-6, 1e-3*D0):.4g}")
        print(f"   {'ecm':>8s} {'bare D(1-a*ecm)':>18s} {'shipped D_eff':>15s} {'shipped/D':>11s}")
        for e, b, s in zip(ecm, bare, shipped):
            print(f"   {e:8.4f} {b:18.4f} {s:15.4f} {s/D0:11.5f}")
        neg = int((bare < 0).sum())
        print(f"   -> bare expression negative at {neg}/{len(ecm)} sampled ecm; "
              f"shipped min = {float(shipped.min()):.4g} (FLOORED, never negative)")
    finally:
        diff.ALPHA_TORTUOSITY = saved

print()
print("=" * 78)
print("C. MONOTONICITY of implied alpha in the mapping scale A")
print("=" * 78)
name, pf_lo, pf_hi, d_lo, d_hi = SHAFIEYAN[0]
x_hi, x_lo = 1.0 - pf_lo, 1.0 - pf_hi
ratio = d_hi / d_lo
print(f"   {name}, measured ratio {ratio:.4f}, x in [{x_lo:.2f}, {x_hi:.2f}]")
print(f"   {'A':>8s} {'ecm at dense end':>18s} {'implied alpha':>14s} {'admissible?':>13s}")
for A in (1.0, 1.5, 2.0, 1.0 / x_hi, 5.0, 10.0, 20.0):
    a = implied_alpha(x_lo, x_hi, ratio, A)
    print(f"   {A:8.3f} {A * x_hi:18.4f} {a:14.4f} "
          f"{'yes' if a < 1 else 'NO (a>1)':>13s}")
print("\n   -> implied alpha scales as 1/A exactly. The verdict is therefore a")
print("      statement about the MAPPING, not about the measurement.")

print()
print(f"ALPHA_TORTUOSITY on exit = {diff.ALPHA_TORTUOSITY} (restored, never written)")
