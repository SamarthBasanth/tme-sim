"""Score the shipped ECM/diffusivity coupling against two external anchors.

SCORING ONLY. Nothing is fitted, nothing under tmesim/ is written. The implied
coefficients printed here are DIAGNOSTICS, never values to install.

Pre-registration: experiments/ecm_coupling_score/PREREGISTRATION.md (frozen
before this script was run).

Anchor 1 - Shafieyan 2014, Biophys J 107(2):485-492, PMID 25028890, PMC4104049.
  Verbatim from the PMC full text (Results, Fig. 6 paragraph):
  "Diffusivity increased from 400 +/- 123 to 778 +/- 111, 110 +/- 19 to
   155 +/- 14, and 124 +/- 28 to 195 +/- 26 um2/s when fluid volume fraction
   increased from 0.66 to 0.97, 0.78 to 0.94, and 0.80 to 0.95 for sodium
   iodide, sodium diatrizoate, and Gd-DTPA, respectively (Fig. 6)."

Anchor 2 - Shi, Zeng, Sun & Fu 2014, J Biomech Eng 136(3):031005, PMID 24193698.
  Verbatim from the PubMed abstract:
  "The corresponding mean ratio of Deff to the free diffusion coefficient
   Dfree, Deff/Dfree, were 0.46, 0.19, 0.12, 0.12, 0.11, and 0.11 for these
   solutes."
"""
from __future__ import annotations

import numpy as np

from tmesim.physics.diffusion import ALPHA_TORTUOSITY, _effective_diffusivity
from tmesim.molecular import compile_drug, library
from tmesim.core import DrugParams

SEP = "=" * 78


# --------------------------------------------------------------------------
# 0. Re-verify the three structural facts the pre-registration rests on
# --------------------------------------------------------------------------
def verify_structure() -> dict:
    print(SEP)
    print("0. STRUCTURAL RE-VERIFICATION (by running the shipped functions)")
    print(SEP)
    print(f"ALPHA_TORTUOSITY (imported from tmesim.physics.diffusion) = {ALPHA_TORTUOSITY}")

    # Drive the real _effective_diffusivity over the full ecm range.
    drug = DrugParams(name="probe", diffusivity=100.0)
    ecm = np.linspace(0.0, 1.0, 11)
    d = _effective_diffusivity(drug, ecm)
    ratio = d / drug.diffusivity
    print("\n ecm   D_eff/D   (from the shipped _effective_diffusivity)")
    for e, r in zip(ecm, ratio):
        print(f" {e:4.2f}   {r:7.4f}")
    lo, hi = float(ratio.min()), float(ratio.max())
    print(f"\n  -> reachable band D_eff/D = [{lo:.4f}, {hi:.4f}]")

    # Out-of-range ecm is clipped, so the band cannot be escaped by a wild field.
    wild = _effective_diffusivity(drug, np.array([-5.0, 0.0, 1.0, 5.0])) / drug.diffusivity
    print(f"  -> ecm = [-5, 0, 1, 5] clips to D_eff/D = {np.round(wild, 4).tolist()}")
    print("     (ecm is np.clip'd to [0,1] inside the function, so the band is hard)")
    return {"band": (lo, hi)}


# --------------------------------------------------------------------------
# 1. Anchor 1 - Shafieyan 2014, FORM level, under three declared mappings
# --------------------------------------------------------------------------
# (solute, phi_fluid_low, phi_fluid_high, D_at_low_phi, D_at_high_phi)
SHAFIEYAN = [
    ("sodium iodide",      0.66, 0.97, 400.0, 778.0),
    ("sodium diatrizoate", 0.78, 0.94, 110.0, 155.0),
    ("Gd-DTPA",            0.80, 0.95, 124.0, 195.0),
]

# Engine's own hard-wired matrix volume fraction (ecm_mesh.ogston_renkin_hindrance
# default, called with radius only at compile_drug.py:174).
PHI_F_ENGINE = 0.05


def implied_alpha(x_lo: float, x_hi: float, ratio: float, A: float) -> float:
    """Solve ratio = (1 - a*A*x_lo) / (1 - a*A*x_hi) for a.

    x = solid volume fraction; x_hi is the DENSE end (low D), x_lo the dilute
    end (high D). ratio = D(x_lo)/D(x_hi) > 1.
    Rearranged:  a = (ratio - 1) / (A * (ratio*x_hi - x_lo))
    """
    denom = A * (ratio * x_hi - x_lo)
    return (ratio - 1.0) / denom


def anchor1() -> dict:
    print()
    print(SEP)
    print("1. ANCHOR 1 - Shafieyan 2014 (cartilage, D vs fluid volume fraction)")
    print("   CRITERION C1a (confirmatory) + C1b (pre-registered, new)")
    print(SEP)

    rows = []
    for name, pf_lo, pf_hi, d_lo, d_hi in SHAFIEYAN:
        # x = solid fraction. Low fluid fraction = HIGH solid = dense end.
        x_hi = 1.0 - pf_lo       # dense end  (e.g. 0.34)
        x_lo = 1.0 - pf_hi       # dilute end (e.g. 0.03)
        ratio = d_hi / d_lo      # measured D ratio, dilute/dense

        A1 = 1.0                  # M1: ecm == solid volume fraction
        A2 = 1.0 / x_hi           # M2: max-normalised, as real_fields.py:114 does
        A3 = 1.0 / PHI_F_ENGINE   # M3: ecm = 1 at the engine's own phi_f = 0.05

        a1 = implied_alpha(x_lo, x_hi, ratio, A1)
        a2 = implied_alpha(x_lo, x_hi, ratio, A2)
        a3 = implied_alpha(x_lo, x_hi, ratio, A3)

        # What the shipped alpha = 0.6 can actually deliver over this interval,
        # under each mapping.
        def engine_ratio(A):
            num = 1.0 - ALPHA_TORTUOSITY * A * x_lo
            den = 1.0 - ALPHA_TORTUOSITY * A * x_hi
            return num / den if den > 0 else float("nan")

        rows.append(dict(
            name=name, pf_lo=pf_lo, pf_hi=pf_hi, x_lo=x_lo, x_hi=x_hi,
            ratio=ratio, a_M1=a1, a_M2=a2, a_M3=a3,
            eng_M1=engine_ratio(A1), eng_M2=engine_ratio(A2), eng_M3=engine_ratio(A3),
            A2=A2,
        ))

    print("\nMeasured endpoints and the solid-fraction interval they span:")
    print(f"{'solute':20s} {'phi_fluid':>12s} {'x=phi_solid':>14s} {'D ratio':>9s}")
    for r in rows:
        print(f"{r['name']:20s} {r['pf_lo']:.2f}->{r['pf_hi']:.2f}   "
              f"{r['x_hi']:.2f}->{r['x_lo']:.2f}      {r['ratio']:7.4f}")

    print("\nC1a CONFIRMATORY - implied coefficient under M1 (ecm == phi_solid):")
    print("   (the scouting record published 1.50 / 1.64 / 2.16)")
    pub = [1.50, 1.64, 2.16]
    c1a_ok = True
    for r, p in zip(rows, pub):
        ok = abs(r["a_M1"] - p) <= 0.02
        c1a_ok &= ok
        print(f"   {r['name']:20s} a_M1 = {r['a_M1']:7.4f}   published {p:.2f}   "
              f"{'MATCH' if ok else 'MISMATCH'}")
    print(f"   -> C1a: {'PASS (arithmetic confirmed)' if c1a_ok else 'FAIL'}")

    print("\nC1b PRE-REGISTERED - implied coefficient under all three mappings:")
    print(f"{'solute':20s} {'M1 (A=1)':>10s} {'M2 (A=1/xmax)':>15s} {'M3 (A=20)':>11s}")
    for r in rows:
        print(f"{r['name']:20s} {r['a_M1']:10.4f} {r['a_M2']:15.4f} {r['a_M3']:11.4f}")

    # DOMAIN CHECK (added after the forward-substitution check in
    # verify_inverse.py caught M3 out of domain): the shipped function clips ecm
    # to [0,1]. A mapping that sends a measured endpoint above ecm=1 is NOT
    # admissible-with-a-small-alpha; it is OUT OF DOMAIN, and the engine then has
    # zero dynamic range over the measured interval.
    print("\n   DOMAIN CHECK - does the mapping keep both endpoints inside ecm<=1?")
    for r in rows:
        for label, A in (("M1", 1.0), ("M2", r["A2"]), ("M3", 1.0 / PHI_F_ENGINE)):
            e_lo, e_hi = A * r["x_lo"], A * r["x_hi"]
            ok = e_hi <= 1.0 + 1e-12
            r[f"domain_{label}"] = ok
            print(f"   {r['name']:20s} {label}: ecm=({e_lo:.4f}, {e_hi:.4f})  "
                  f"{'in domain' if ok else 'OUT OF DOMAIN (clipped to 1.0)'}")

    verdicts = {}
    for key, label in (("a_M1", "M1"), ("a_M2", "M2"), ("a_M3", "M3")):
        vals = [r[key] for r in rows]
        in_domain = all(r[f"domain_{label}"] for r in rows)
        if not in_domain:
            verdicts[label] = "OUT OF DOMAIN (ecm clips to 1; engine saturates)"
        else:
            all_gt1 = all(v > 1.0 for v in vals)
            all_lt1 = all(v < 1.0 for v in vals)
            verdicts[label] = ("INADMISSIBLE (a>1)" if all_gt1
                               else "ADMISSIBLE (a<1)" if all_lt1
                               else "SPLIT within mapping")
        print(f"   {label}: a in [{min(vals):.4f}, {max(vals):.4f}]  -> {verdicts[label]}")

    distinct = set(verdicts.values())
    if len(distinct) > 1:
        c1b = "MAPPING-DEPENDENT -> form claim is a NON-RESULT"
    elif "INADMISSIBLE (a>1)" in distinct:
        c1b = "LINEAR FORM REFUTED under every declared mapping"
    else:
        c1b = "slope ADMISSIBLE under every declared mapping"
    print(f"\n   -> C1b VERDICT: {c1b}")

    # Where does the verdict flip? a_M1 * x_ref = 1.
    print("\n   Flip point (the solid fraction that ecm=1 must denote for a=1):")
    for r in rows:
        print(f"   {r['name']:20s} x_ref = 1/a_M1 = {1.0 / r['a_M1']:.4f}")
    print("   -> above that solid fraction the slope is inadmissible; below it, admissible.")

    # Amplitude statement under the code-faithful mapping M2.
    print("\n   Under M2 (the mapping faithful to density /= dmax), the shipped 0.6 is:")
    for r in rows:
        print(f"   {r['name']:20s} 0.6 / {r['a_M2']:.4f} = {ALPHA_TORTUOSITY / r['a_M2']:.2f}x "
              f"too steep")
    print("   Under M3 (engine's own phi_f = 0.05):")
    for r in rows:
        print(f"   {r['name']:20s} 0.6 / {r['a_M3']:.4f} = {ALPHA_TORTUOSITY / r['a_M3']:.2f}x "
              f"too steep")

    # ---------------------------------------------------------------
    # HEADLINE statistic: no inverse solve at all. Just run the SHIPPED
    # alpha = 0.6 forward through the shipped function and compare the
    # predicted endpoint ratio against the measured one, under each mapping.
    # Nothing is fitted; this is pure scoring.
    # ---------------------------------------------------------------
    print("\n   HEADLINE (no inverse solve, no fitting): the SHIPPED alpha = 0.6")
    print("   run FORWARD through _effective_diffusivity, predicted vs measured ratio.")
    print(f"\n{'solute':20s} {'mapping':8s} {'ecm_lo':>8s} {'ecm_hi':>8s} "
          f"{'pred ratio':>11s} {'measured':>9s} {'error':>9s} {'sign':>7s}")
    drug = DrugParams(name="probe", diffusivity=100.0)
    headline = []
    for r in rows:
        for label, A in (("M1", 1.0), ("M2", r["A2"]), ("M3", 1.0 / PHI_F_ENGINE)):
            e = np.array([A * r["x_lo"], A * r["x_hi"]])
            d = _effective_diffusivity(drug, e)
            pred = float(d[0] / d[1])
            meas = r["ratio"]
            if pred >= meas:
                err, sign = pred / meas, "STEEP"
            else:
                err, sign = meas / pred, "FLAT"
            e_cl = np.clip(e, 0.0, 1.0)
            print(f"{r['name']:20s} {label:8s} {e_cl[0]:8.4f} {e_cl[1]:8.4f} "
                  f"{pred:11.4f} {meas:9.4f} {err:8.3f}x {sign:>7s}")
            headline.append(dict(solute=r["name"], mapping=label, pred=pred,
                                 meas=meas, err=err, sign=sign))

    errs = [h["err"] for h in headline]
    signs = sorted({h["sign"] for h in headline})
    print(f"\n   -> discrepancy MAGNITUDE spans [{min(errs):.3f}x, {max(errs):.3f}x] "
          f"across all 9 solute x mapping cells.")
    print(f"   -> discrepancy SIGN takes values {signs}: it is set by the MAPPING,")
    print("      not by the measurement. Same data, opposite conclusion.")
    for label in ("M1", "M2", "M3"):
        sub = [h for h in headline if h["mapping"] == label]
        s = sorted({h['sign'] for h in sub})
        print(f"      {label}: {s}, magnitude "
              f"[{min(h['err'] for h in sub):.3f}x, {max(h['err'] for h in sub):.3f}x]")

    return {"rows": rows, "c1a": c1a_ok, "c1b": c1b, "verdicts": verdicts,
            "headline": headline, "err_span": (min(errs), max(errs)), "signs": signs}


# --------------------------------------------------------------------------
# 2. Anchor 2 - Shi 2014, BAND level, with the DENOMINATOR fixed
# --------------------------------------------------------------------------
SHI_RATIOS = [
    ("sodium fluorescein (376 Da)", 0.46),
    ("FITC-dextran 4k",             0.19),
    ("FITC-dextran 20k",            0.12),
    ("FITC-dextran 40k",            0.12),
    ("FITC-dextran 70k",            0.11),
    ("Alexa-488 IgG (~160 kDa)",    0.11),
]

FLAGSHIP_MEAN_ECM = 0.076   # results/heldout_patient.md, section 1 (11-marker panel)


def anchor2(band: tuple) -> dict:
    print()
    print(SEP)
    print("2. ANCHOR 2 - Shi 2014 (rat brain, Deff/Dfree around single microvessels)")
    print("   CRITERION C2a (confirmatory) + C2b (pre-registered, new)")
    print(SEP)

    lo, hi = band
    print(f"\nC2a CONFIRMATORY - spatial factor alone spans [{lo:.2f}, {hi:.2f}]:")
    n_out = 0
    for name, r in SHI_RATIOS:
        inside = lo <= r <= hi
        n_out += (not inside)
        print(f"   {name:30s} Deff/Dfree = {r:.2f}   {'inside' if inside else 'OUTSIDE'}")
    c2a_ok = (n_out == 4)
    print(f"   -> {n_out} of 6 outside; claim was 4. C2a: {'PASS' if c2a_ok else 'FAIL'}")

    print("\nC2b PRE-REGISTERED - the engine's FULL hindrance chain.")
    print("   Shi's denominator is FREE-WATER Dfree. The engine's band denominator")
    print("   is the already-hindered COMPILED diffusivity. Comparable quantity:")
    print("     (Deff/Dfree)_engine = (compiled_D / free_D) * (1 - 0.6*ecm)")

    chain = {}
    for drug_name in ("doxorubicin", "temozolomide"):
        mol = library.get_molecule(drug_name)
        prof = compile_drug.build_profile(mol)
        free_D = float(prof.free_diffusivity)
        comp_D = float(prof.effective_diffusivity)
        h = comp_D / free_D
        print(f"\n   {drug_name}:")
        print(f"     free_D (Stokes-Einstein, T=310K) = {free_D:9.2f} um^2/s")
        print(f"     compiled effective_D             = {comp_D:9.2f} um^2/s")
        print(f"     compiled_D / free_D              = {h:9.5f}  "
              f"({1.0 / h:.2f}x total pre-spatial hindrance)")
        full_lo, full_hi = h * lo, h * hi
        at_flag = h * (1.0 - ALPHA_TORTUOSITY * FLAGSHIP_MEAN_ECM)
        print(f"     FULL-CHAIN Deff/Dfree band       = [{full_lo:.5f}, {full_hi:.5f}]")
        print(f"     FULL-CHAIN at flagship ecm={FLAGSHIP_MEAN_ECM:.3f}  = {at_flag:.5f}")
        chain[drug_name] = dict(free_D=free_D, comp_D=comp_D, h=h,
                                full_lo=full_lo, full_hi=full_hi, at_flag=at_flag)

    ref = chain["doxorubicin"]
    print("\n   C2b-i: does any Shi ratio fall inside the full-chain band "
          f"[{ref['full_lo']:.5f}, {ref['full_hi']:.5f}] (doxorubicin)?")
    inside_any = [n for n, r in SHI_RATIOS if ref["full_lo"] <= r <= ref["full_hi"]]
    print(f"     -> {inside_any if inside_any else 'NONE. All six measured ratios are ABOVE the band.'}")

    smallest = min(r for _, r in SHI_RATIOS)
    print(f"     -> smallest measured ratio {smallest:.2f} vs engine max "
          f"{ref['full_hi']:.5f}: engine is {smallest / ref['full_hi']:.2f}x SLOWER")

    shi_small = 0.46
    factor = shi_small / ref["at_flag"]
    print(f"\n   C2b-ii: Shi small-solute 0.46 vs engine full chain at flagship ecm "
          f"({ref['at_flag']:.5f}):")
    print(f"     -> factor = {factor:.2f}x. Engine is "
          f"{'SLOWER (too hindered)' if factor > 1 else 'FASTER (too free)'} "
          f"than the measurement.")
    print(f"     Scouting record claimed '~2.1x too FAST'. Sign here is "
          f"{'OPPOSITE' if factor > 1.2 else 'consistent'}.")

    # Consequence for the 6x Deff/k partition.
    DEFF_OVER_K = 6.0
    print(f"\n   Consequence for the repo's own Deff/k ~{DEFF_OVER_K:.0f}x-too-large "
          "(results/kyle2004_penetration.md):")
    if factor > 1.0:
        k_factor = DEFF_OVER_K * factor
        print(f"     If Deff is {factor:.2f}x too SMALL, then to keep Deff/k {DEFF_OVER_K:.0f}x")
        print(f"     too LARGE, k must be too small by {k_factor:.1f}x "
              f"- NOT the ~2.9x the scouting record derived.")
    else:
        print(f"     Deff {1/factor:.2f}x too large leaves "
              f"{DEFF_OVER_K * factor:.2f}x in k.")

    print("\n   C2c: Shi's tissue is NORMAL RAT BRAIN, not tumour; his small solute is")
    print("        sodium fluorescein (376 Da), not a repo compound. Cross-species,")
    print("        cross-tissue, cross-compound. Stated, not hidden.")

    return {"c2a": c2a_ok, "chain": chain, "factor": factor,
            "inside_any": inside_any}


def main():
    s = verify_structure()
    a1 = anchor1()
    a2 = anchor2(s["band"])

    print()
    print(SEP)
    print("SUMMARY AGAINST THE PRE-REGISTERED CRITERIA")
    print(SEP)
    print(f"C1a (confirmatory, M1 arithmetic)  : {'PASS' if a1['c1a'] else 'FAIL'}")
    print(f"C1b (pre-registered, mapping test) : {a1['c1b']}")
    print(f"C2a (confirmatory, band claim)     : {'PASS' if a2['c2a'] else 'FAIL'}")
    print(f"C2b (pre-registered, full chain)   : engine {a2['factor']:.2f}x "
          f"{'SLOWER' if a2['factor'] > 1 else 'FASTER'} than Shi's small solute")
    print(f"      Shi ratios inside full-chain band: "
          f"{a2['inside_any'] if a2['inside_any'] else 'NONE'}")
    print()
    print("NO PARAMETER WAS TUNED. ALPHA_TORTUOSITY is still "
          f"{ALPHA_TORTUOSITY} and was not written to.")


if __name__ == "__main__":
    main()
