"""The BINDING-SITE BARRIER: why high-affinity antibodies penetrate tumors WORSE
than small molecules, even though the engine would otherwise only slow them down
via their (already low) diffusivity.

Motivation
----------
The rest of the molecular layer hands the tissue engine an antibody with a small
*effective diffusivity* (a ~5.5 nm IgG is hindered far more by the ECM than a
~0.4 nm drug). Low diffusivity alone, though, only *slows* transport -- given
enough time an antibody would still fill the tissue. That is NOT what is seen
experimentally. Real anti-tumor antibodies pile up in a thin perivascular cuff
and leave an antibody-free tumor core (Fujimori, Covell, Fletcher & Weinstein,
*J Nucl Med* 1990; Weinstein & van Osdol, *Cancer Res* 1992; Juweid et al. 1992
"direct experimental evidence for a binding site barrier"). The mechanism is
SATURABLE, HIGH-AFFINITY ANTIGEN BINDING: the first antibody to arrive binds the
nearest antigen essentially irreversibly, so it cannot advance until the local
antigen is saturated. Binding is a *sink*, not merely a *drag*.

Reaction-diffusion / Thiele-modulus framing
--------------------------------------------
Treat antibody concentration ``C`` near a vessel as steady reaction-diffusion
with a first-order antigen-binding loss ``k_ret`` (1/s):

    D d^2C/dx^2  -  k_ret C  =  0        =>        C(x) = C0 * exp(-x / L)

with the characteristic **penetration length**

    L  =  sqrt( D / k_ret )                          [um]

This is exactly the Thiele-modulus / Damkohler picture used for antibody tumor
penetration (Thurber, Schmidt & Wittrup, *Adv Drug Deliv Rev* 2008; *J Nucl Med*
2007 "Theoretic Criteria for Antibody Penetration"): define a tissue Thiele
modulus ``phi = R / L = R * sqrt(k_ret / D)`` over a length scale ``R`` (e.g. the
inter-vessel half-distance). ``phi >> 1`` means binding consumes the antibody
long before it diffuses across the tissue -> a **binding-site-barrier-limited**,
perivascular-only distribution. ``phi <~ 1`` means the drug saturates the tissue.

The barrier retention rate (pseudo-first-order antigen binding)
---------------------------------------------------------------
For an antibody meeting antigen at local concentration ``[Ag]`` with association
rate ``kon`` and dissociation rate ``koff`` (Kd = koff/kon):

    k_forward       =  kon * [Ag]                    (pseudo-first-order, 1/s)
    irreversibility  =  k_forward / (k_forward + koff)   in [0, 1]
    k_ret            =  k_forward * irreversibility   +   baseline_loss

* ``k_forward`` is the rate at which free antibody is captured perivascularly.
  It rises with antigen density -> denser targets build a stronger barrier.
* ``irreversibility`` encodes AFFINITY: a high-affinity binder (small ``koff``,
  low Kd) stays bound once captured, so binding is a near-perfect sink
  (irreversibility -> 1). A low-affinity binder releases and "hops" forward
  (irreversibility -> 0), which is why raising ``koff`` / lowering ``kon`` is the
  textbook way to *relieve* the binding-site barrier and homogenise uptake
  (Thurber 2008). This is the mechanistic gate that makes *high-affinity*
  binders penetrate worse.

Literature-anchored defaults (Thurber, *J Nucl Med* 2007; Schmidt & Wittrup):
  antigen [Ag] ~ 0.3 uM, kon ~ 1e5 1/(M*s), Kd ~ 1 nM (=> koff ~ 1e-4 1/s) for a
  typical IgG against a solid-tumor antigen. These give k_forward ~ 0.03 1/s and
  L ~ 5-15 um for an antibody with tissue D ~ 1-10 um^2/s -- i.e. a perivascular
  cuff a few cell diameters thick, matching the observed barrier.

Small molecules
---------------
A small molecule with NO specific, high-capacity, high-affinity target does not
experience this sink: its only tissue loss is the generic, low-capacity binding
/ clearance already in ``DrugParams`` (~1e-4 1/s). ``k_ret`` adds essentially
nothing, its penetration length stays diffusion/clearance-limited (hundreds of
um to mm), and ``is_barrier_limited`` is False. That is the whole point: the
barrier hurts the antibody, not the small molecule.

Scope / honesty
---------------
This is a lumped, steady-state *effective-medium* description (one perivascular
sink term), not a spatially resolved moving-front / shrinking-core solve. It
captures the correct mechanism and scaling (L ~ sqrt(D/k_ret), affinity- and
antigen-density-dependent) and the correct RANK (antibody << small molecule),
which is what the tissue engine needs to stop treating an antibody as merely a
slow small molecule. It does NOT modify the frozen engine, core, or any existing
molecular module -- it only produces a new ``core.DrugParams`` with an increased
first-order ``binding`` sink, which the engine already integrates as
``k = clearance + binding`` in its reaction term.

Units: length um, time s, diffusivity um^2/s, concentration mol/L (M).
"""
from __future__ import annotations

import math
from dataclasses import replace
from typing import Dict, Optional, Union

from ..core import DrugParams
from .types import Molecule

# ---------------------------------------------------------------------------
# Literature-anchored defaults (Thurber & Wittrup lineage). See module docstring.
# ---------------------------------------------------------------------------
_DEFAULT_ANTIGEN_M = 3.0e-7      # local perivascular antigen conc, mol/L (~0.3 uM)
_DEFAULT_KON = 1.0e5             # antibody-antigen association rate, 1/(M*s)
_DEFAULT_KOFF = 1.0e-4           # dissociation rate, 1/s (Kd = koff/kon = 1 nM)
_DEFAULT_BASELINE_LOSS = 5.0e-4  # generic clearance+binding when no barrier, 1/s

# Antibody effective tissue diffusivity fallback (um^2/s) if none supplied. IgG
# in tumor interstitium is ~1-10 um^2/s (heavily hindered); small molecules
# ~100-500 um^2/s. Only used when we cannot read a diffusivity from the caller.
_ANTIBODY_D_FALLBACK = 5.0
_SMALLMOL_D_FALLBACK = 300.0

# Length scale for the tissue Thiele modulus: an inter-vessel half-distance in
# tumor tissue is ~100-200 um. Used only to report phi / is_barrier_limited.
_TISSUE_LENGTH_UM = 150.0

# A high-affinity binder is one whose Kd is at/below this (mol/L). Antibodies are
# always treated as barrier candidates; a small molecule is only treated as one
# if the caller explicitly flags it (high_affinity=True) or passes a target.
_HIGH_AFFINITY_KD_M = 1.0e-8     # 10 nM


def _descriptors_of(molecule_or_descriptors: Union[Molecule, Dict]) -> Dict:
    """Accept either a Molecule or a raw descriptors dict; return descriptors."""
    if isinstance(molecule_or_descriptors, Molecule):
        return molecule_or_descriptors.descriptors
    if isinstance(molecule_or_descriptors, dict):
        return molecule_or_descriptors
    raise TypeError(
        "expected a Molecule or a descriptors dict, got "
        f"{type(molecule_or_descriptors).__name__}"
    )


def binding_site_barrier_factor(
    molecule_or_descriptors: Union[Molecule, Dict],
    antigen_density_M: float = _DEFAULT_ANTIGEN_M,
    kon_per_M_s: float = _DEFAULT_KON,
    koff_per_s: float = _DEFAULT_KOFF,
    diffusivity_um2_s: Optional[float] = None,
    baseline_loss_per_s: float = _DEFAULT_BASELINE_LOSS,
    high_affinity: Optional[bool] = None,
    tissue_length_um: float = _TISSUE_LENGTH_UM,
) -> Dict:
    """Perivascular retention from saturable high-affinity antigen binding.

    Implements the binding-site-barrier / Thiele-modulus model (see module
    docstring). For antibodies and explicitly high-affinity binders it returns a
    strong effective tissue-retention rate ``k_ret`` and the resulting
    penetration length ``L = sqrt(D / k_ret)``. Small molecules with no specific
    high-capacity target get ``k_ret ~ 0`` (negligible extra retention) and keep
    a diffusion/clearance-limited penetration length.

    Parameters
    ----------
    molecule_or_descriptors : Molecule | dict
        Carries ``is_antibody`` (and optionally ``radius_hydrodynamic_nm``, ``mw``).
    antigen_density_M : float
        Local perivascular antigen concentration [Ag], mol/L. Denser target ->
        stronger barrier. Default ~0.3 uM (Thurber, J Nucl Med 2007).
    kon_per_M_s, koff_per_s : float
        Antibody-antigen association / dissociation rates. Kd = koff/kon. Default
        kon 1e5, koff 1e-4 -> Kd 1 nM (a typical high-affinity IgG).
    diffusivity_um2_s : float or None
        Effective tissue diffusivity. If None, a size-appropriate fallback is used
        (antibody ~5, small molecule ~300 um^2/s).
    baseline_loss_per_s : float
        Generic first-order tissue loss (clearance+binding) present even without a
        specific target; sets the small-molecule (barrier-free) penetration length.
    high_affinity : bool or None
        Force whether to treat this molecule as a barrier-forming high-affinity
        binder. If None, defaults to ``is_antibody`` (antibodies -> True). Set True
        to model a small high-affinity targeted agent; leave None/False for a plain
        cytotoxic small molecule.
    tissue_length_um : float
        Inter-vessel length scale ``R`` for the Thiele modulus phi = R / L.

    Returns
    -------
    dict with keys
        retention_rate_1_per_s : float
            Effective first-order perivascular retention/loss rate ``k_ret``
            attributable to saturable antigen binding (0 for non-barrier drugs).
            This is the quantity ``apply_barrier_to_drugparams`` adds to binding.
        penetration_length_um : float
            L = sqrt(D / (k_ret + baseline)) -- how far the drug reaches from a
            vessel before the sink consumes it.
        is_barrier_limited : bool
            True when antigen binding dominates the loss AND phi = R/L > 1, i.e.
            the drug is trapped in a perivascular cuff much thinner than the tissue.
        note : str
            Human-readable provenance / explanation.
        (plus diagnostics: kd_M, k_forward_1_per_s, irreversibility, thiele_modulus,
         baseline_penetration_length_um, diffusivity_um2_s, is_antibody, is_high_affinity)
    """
    d = _descriptors_of(molecule_or_descriptors)
    is_antibody = bool(d.get("is_antibody", False))

    # Who experiences the barrier? Antibodies always qualify; a small molecule
    # only if the caller explicitly flags a high-affinity specific target.
    if high_affinity is None:
        is_high_affinity = is_antibody
    else:
        is_high_affinity = bool(high_affinity) or is_antibody

    # Effective tissue diffusivity (um^2/s).
    if diffusivity_um2_s is not None:
        D = float(diffusivity_um2_s)
    else:
        D = _ANTIBODY_D_FALLBACK if is_antibody else _SMALLMOL_D_FALLBACK
    D = max(D, 1e-9)

    baseline = max(float(baseline_loss_per_s), 0.0)
    kd_M = (float(koff_per_s) / float(kon_per_M_s)) if kon_per_M_s > 0 else float("inf")

    if not is_high_affinity:
        # No specific high-capacity target: binding adds essentially nothing.
        k_forward = 0.0
        irreversibility = 0.0
        k_ret = 0.0
        note = (
            "No specific high-capacity/high-affinity target: binding-site barrier "
            "does not apply. Penetration stays diffusion/clearance-limited "
            f"(L = sqrt(D/baseline), D={D:.1f} um^2/s, baseline={baseline:.1e}/s). "
            "Small molecules are essentially unaffected by the barrier."
        )
    else:
        # Pseudo-first-order antigen capture rate.
        k_forward = float(kon_per_M_s) * float(antigen_density_M)  # 1/s
        # Affinity gate: fraction that stays bound instead of releasing forward.
        denom = k_forward + float(koff_per_s)
        irreversibility = (k_forward / denom) if denom > 0 else 0.0
        k_ret = k_forward * irreversibility
        affinity_word = "high" if kd_M <= _HIGH_AFFINITY_KD_M else "moderate/low"
        note = (
            f"Binding-site barrier ACTIVE ({'antibody' if is_antibody else 'flagged high-affinity binder'}). "
            f"Saturable antigen binding: [Ag]={antigen_density_M:.2e} M, kon={kon_per_M_s:.1e}/(M*s), "
            f"koff={koff_per_s:.1e}/s -> Kd={kd_M:.2e} M ({affinity_word} affinity). "
            f"k_forward={k_forward:.3e}/s, irreversibility={irreversibility:.3f} -> "
            f"retention k_ret={k_ret:.3e}/s. Penetration L=sqrt(D/k_ret) collapses "
            "to a thin perivascular cuff (Fujimori/Weinstein binding-site barrier; "
            "Thurber/Wittrup Thiele-modulus criterion)."
        )

    total_loss = k_ret + baseline
    pen_len = math.sqrt(D / total_loss) if total_loss > 0 else float("inf")
    baseline_len = math.sqrt(D / baseline) if baseline > 0 else float("inf")

    # Tissue Thiele modulus phi = R / L; >1 => binding consumes drug before it
    # crosses the inter-vessel distance (perivascular-limited).
    thiele = (float(tissue_length_um) / pen_len) if math.isfinite(pen_len) and pen_len > 0 else 0.0
    # Barrier-limited iff the specific sink dominates baseline AND phi > 1.
    is_barrier_limited = bool(k_ret > baseline and thiele > 1.0)

    return {
        "retention_rate_1_per_s": k_ret,
        "penetration_length_um": pen_len,
        "is_barrier_limited": is_barrier_limited,
        "note": note,
        # diagnostics
        "kd_M": kd_M,
        "k_forward_1_per_s": k_forward,
        "irreversibility": irreversibility,
        "thiele_modulus": thiele,
        "baseline_penetration_length_um": baseline_len,
        "diffusivity_um2_s": D,
        "baseline_loss_1_per_s": baseline,
        "is_antibody": is_antibody,
        "is_high_affinity": is_high_affinity,
    }


def apply_barrier_to_drugparams(
    drug_params: DrugParams,
    molecule: Union[Molecule, Dict],
    **barrier_kw,
) -> DrugParams:
    """Return a NEW ``core.DrugParams`` with the binding-site-barrier retention
    folded into the first-order ``binding`` sink.

    The engine integrates ``dC/dt = div(D grad C) + S - (clearance + binding) C``,
    so adding the perivascular antigen-binding retention ``k_ret`` to ``binding``
    makes the reaction term dominate transport for an antibody -- its penetration
    length collapses from ``sqrt(D/baseline)`` (hundreds of um) to
    ``sqrt(D/(baseline + k_ret))`` (a few um). Small molecules get ``k_ret ~ 0``
    and are returned essentially unchanged.

    The molecule's own ``drug_params.diffusivity`` is used as the effective
    tissue D (so the reported penetration length is consistent with the engine),
    and ``drug_params.binding + drug_params.clearance`` is used as the barrier-free
    baseline loss. Any keyword in ``barrier_kw`` (antigen_density_M, kon_per_M_s,
    koff_per_s, high_affinity, ...) is forwarded to ``binding_site_barrier_factor``.

    Does NOT mutate the input; uses ``dataclasses.replace``.
    """
    baseline = float(drug_params.binding) + float(drug_params.clearance)
    factor = binding_site_barrier_factor(
        molecule,
        diffusivity_um2_s=drug_params.diffusivity,
        baseline_loss_per_s=baseline,
        **barrier_kw,
    )
    k_ret = float(factor["retention_rate_1_per_s"])
    if k_ret <= 0.0:
        # Small molecule / no barrier: return an equivalent copy (unchanged binding).
        return replace(drug_params)
    return replace(drug_params, binding=float(drug_params.binding) + k_ret)


if __name__ == "__main__":
    from . import library
    from .compile_drug import molecule_to_drugparams

    print("=" * 74)
    print("BINDING-SITE BARRIER self-test: antibody trapped perivascularly vs")
    print("small molecule (diffusion/clearance-limited). L = sqrt(D / k_loss).")
    print("=" * 74)

    bev_mol = library.get_molecule("bevacizumab")   # ~149 kDa anti-VEGF IgG
    tmz_mol = library.get_molecule("temozolomide")  # ~194 Da small molecule

    # Compile each molecule through the real bridge to get engine DrugParams
    # (this is what the tissue engine would otherwise consume, barrier-free).
    bev_dp = molecule_to_drugparams(bev_mol)
    tmz_dp = molecule_to_drugparams(tmz_mol)

    print(f"\nCompiled diffusivities (hindered, barrier-free):")
    print(f"  bevacizumab  D_eff = {bev_dp.diffusivity:8.2f} um^2/s  "
          f"(binding={bev_dp.binding:.2e}, clearance={bev_dp.clearance:.2e})")
    print(f"  temozolomide D_eff = {tmz_dp.diffusivity:8.2f} um^2/s  "
          f"(binding={tmz_dp.binding:.2e}, clearance={tmz_dp.clearance:.2e})")

    # ------------------------------------------------------------------ #
    # WITHOUT the barrier: penetration limited ONLY by diffusivity + the
    # generic clearance/binding baseline already in DrugParams.
    # ------------------------------------------------------------------ #
    bev_base_loss = bev_dp.binding + bev_dp.clearance
    tmz_base_loss = tmz_dp.binding + tmz_dp.clearance
    bev_no_barrier = math.sqrt(bev_dp.diffusivity / bev_base_loss)
    tmz_no_barrier = math.sqrt(tmz_dp.diffusivity / tmz_base_loss)

    # ------------------------------------------------------------------ #
    # WITH the barrier.
    # ------------------------------------------------------------------ #
    bev_bar = binding_site_barrier_factor(
        bev_mol, diffusivity_um2_s=bev_dp.diffusivity, baseline_loss_per_s=bev_base_loss
    )
    tmz_bar = binding_site_barrier_factor(
        tmz_mol, diffusivity_um2_s=tmz_dp.diffusivity, baseline_loss_per_s=tmz_base_loss
    )

    bev_with = bev_bar["penetration_length_um"]
    tmz_with = tmz_bar["penetration_length_um"]

    print("\n--- penetration length L (um) : WITHOUT barrier -> WITH barrier ---")
    print(f"  bevacizumab (antibody)   : {bev_no_barrier:10.1f}  ->  {bev_with:10.2f}   "
          f"barrier_limited={bev_bar['is_barrier_limited']}  Thiele={bev_bar['thiele_modulus']:.1f}")
    print(f"  temozolomide (small mol) : {tmz_no_barrier:10.1f}  ->  {tmz_with:10.2f}   "
          f"barrier_limited={tmz_bar['is_barrier_limited']}  Thiele={tmz_bar['thiele_modulus']:.2f}")
    print(f"\n  bevacizumab barrier note : {bev_bar['note']}")
    print(f"  temozolomide barrier note: {tmz_bar['note']}")

    # DrugParams effect (what the engine actually sees).
    bev_dp2 = apply_barrier_to_drugparams(bev_dp, bev_mol)
    tmz_dp2 = apply_barrier_to_drugparams(tmz_dp, tmz_mol)
    print("\n--- DrugParams.binding : before -> after barrier ---")
    print(f"  bevacizumab  : {bev_dp.binding:.3e} -> {bev_dp2.binding:.3e} /s  "
          f"(added k_ret = {bev_bar['retention_rate_1_per_s']:.3e} /s)")
    print(f"  temozolomide : {tmz_dp.binding:.3e} -> {tmz_dp2.binding:.3e} /s  "
          f"(added k_ret = {tmz_bar['retention_rate_1_per_s']:.3e} /s)")

    # ------------------------------------------------------------------ #
    # ASSERTIONS
    # ------------------------------------------------------------------ #
    # 1. Without the barrier, the antibody is limited ONLY by (low) diffusivity:
    #    given its small baseline loss it would still reach a long way from a vessel.
    assert bev_no_barrier > 50.0, (
        f"barrier-free antibody should reach far (>50 um), got {bev_no_barrier:.1f} um")

    # 2. The barrier collapses the antibody's penetration length dramatically.
    assert bev_with < bev_no_barrier, "barrier must SHORTEN antibody penetration"
    assert bev_with < 0.25 * bev_no_barrier, (
        f"barrier should collapse antibody penetration by >4x: "
        f"{bev_no_barrier:.1f} -> {bev_with:.1f} um")
    assert bev_bar["is_barrier_limited"], "antibody must be flagged barrier-limited"

    # 3. Antibody WITH barrier penetrates MUCH LESS than a small molecule.
    assert bev_with < tmz_with, (
        f"antibody L ({bev_with:.2f}) must be << small-molecule L ({tmz_with:.1f})")
    assert bev_with < 0.05 * tmz_with, (
        f"antibody L should be <<5% of small-molecule L: "
        f"{bev_with:.2f} vs {tmz_with:.1f} um")

    # 4. Small molecule is essentially UNAFFECTED by the barrier: negligible extra
    #    retention, and DrugParams.binding unchanged.
    assert tmz_bar["retention_rate_1_per_s"] == 0.0, (
        f"small-molecule barrier retention must be ~0, got "
        f"{tmz_bar['retention_rate_1_per_s']:.3e}")
    assert not tmz_bar["is_barrier_limited"], "small molecule must NOT be barrier-limited"
    assert abs(tmz_with - tmz_no_barrier) < 1e-6, (
        "small-molecule penetration must be unchanged by the barrier")
    assert tmz_dp2.binding == tmz_dp.binding, "small-molecule DrugParams.binding must be unchanged"

    # 5. Antibody DrugParams.binding is strongly increased (drives the engine collapse).
    assert bev_dp2.binding > bev_dp.binding, "antibody binding must increase under barrier"
    assert bev_dp2.binding > 100.0 * bev_dp.binding, (
        "antibody binding sink should dominate its tiny baseline")

    print("\nAll binding-site-barrier assertions passed.")
