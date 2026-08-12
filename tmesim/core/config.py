"""Configuration dataclasses and the drug library schema.

All physical parameters carry explicit units in comments. Defaults are chosen to
be biologically plausible for a small-molecule chemotherapeutic in brain tumor
tissue (order-of-magnitude, not patient-calibrated — calibration is exactly the
Stage-1+ work this scaffold is built to enable).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional


@dataclass
class DrugParams:
    """Transport + pharmacodynamic parameters for one agent.

    diffusivity : effective free diffusion coefficient in tissue, um^2/s.
        (Small molecules ~ 100-500 um^2/s; antibodies ~ 1-10 um^2/s.)
    molecular_weight : Da (informational / for scaling heuristics).
    ic50 : concentration giving half-maximal kill effect, uM.
    hill : Hill coefficient for the dose-response curve.
    binding : first-order tissue binding/uptake rate, 1/s (retention).
    clearance : first-order clearance/metabolism rate in tissue, 1/s.
    bbb_penetrance : dimensionless brain:plasma PARTITION coefficient in [0, 1]
        — the quasi-equilibrium tissue:plasma ratio implied by the drug's
        structure (blood-brain-barrier factor — central to glioblastoma).
        ``molecular/permeability.py`` computes it from MW / logP / TPSA / HBD and
        calibrates it against measured tissue:plasma RATIOS, so it is **NOT a
        permeability** and **NOT a rate**: it carries no units of its own, and
        neither does ``vascular_density``, which is a volume fraction. The
        first-order transvascular rate the transport PDE needs is built from
        both as ``PS = PS_SCALE_S_INV * vascular_density * bbb_penetrance``
        (1/s), where ``treatment.drugs.PS_SCALE_S_INV`` carries the whole
        ``1/s`` — see that constant's docstring, ``reconstruct/barrier.py``, and
        ``results/transvascular_exchange_revalidation.md`` section 6.2 for why
        the two readings must not be conflated.
    hypoxia_resistance : optional PER-DRUG override of
        ``PDConfig.hypoxia_resistance`` (see that field's docstring for the
        general unvalidated-prior caveat). ``None`` (the default for every drug
        except temozolomide, see :data:`DRUG_LIBRARY`) means "use the global
        ``PDConfig`` default" — i.e. this field changes nothing for a drug that
        does not set it. It exists because the literature basis found in
        ``results/sweep/state_multiplier_literature_check_*.md`` is
        drug-specific (real GBM-cell-line hypoxia-vs-normoxia TMZ data), not a
        generic chemotherapy constant, and ``PDConfig`` has no drug axis of its
        own. ``tmesim.cells.pkpd.apply_drug_kill`` (and every module that
        recomputes the same resistance vector directly, e.g.
        ``validation.kill_calibration.resistance_vector``) resolves the
        effective multiplier as ``drug.hypoxia_resistance if drug.hypoxia_resistance
        is not None else pd.hypoxia_resistance``.
    """

    name: str
    diffusivity: float = 250.0
    molecular_weight: float = 300.0
    ic50: float = 5.0
    hill: float = 1.5
    binding: float = 5e-4
    clearance: float = 2e-4
    bbb_penetrance: float = 0.3
    hypoxia_resistance: Optional[float] = None


@dataclass
class PKConfig:
    """Systemic pharmacokinetics setting the vascular (blood) concentration that
    drives tissue delivery. A simple 1-compartment exponential by default.

    plasma_peak : peak plasma concentration after a dose, uM.
    half_life : plasma elimination half-life, seconds.
    n_doses : number of doses in the treatment course.
    dose_interval : seconds between doses.
    """

    plasma_peak: float = 50.0
    half_life: float = 3600.0 * 4
    n_doses: int = 1
    dose_interval: float = 3600.0 * 24


@dataclass
class PDEConfig:
    """Numerical settings for the drug-transport PDE solve."""

    dt: float = 0.5              # time step, seconds
    total_time: float = 3600.0 * 12  # simulated exposure window, seconds
    convection: bool = True      # include interstitial-flow advection
    save_every: int = 200        # store history every N steps
    #: Drive the transport source with the full time-resolved plasma curve
    #: ``plasma(t)`` instead of its window time-average. ``False`` (default)
    #: preserves the historical Stage-0 plasma-averaging behaviour BIT-FOR-BIT;
    #: ``True`` lets tissue concentration genuinely rise and fall, which is what
    #: makes dosing SCHEDULE (fractionation, infusion duration) resolvable at
    #: all. See :func:`tmesim.sim.engine.simulate_response` and the
    #: ``source_profile`` argument of
    #: :func:`tmesim.physics.diffusion.simulate_drug_transport`.
    time_resolved: bool = False
    #: Treat vascular delivery as a genuine gradient-driven EXCHANGE,
    #: ``S = PS*(Cp - C)``, instead of the historical one-way leak ``S = PS*Cp``.
    #:
    #: The one-way form has no path back from interstitium to blood, so its
    #: steady state is ``PS*Cp/(clearance + binding)`` -- unbounded above, and
    #: measurably wrong: against human breast microdialysis it puts
    #: 5-fluorouracil at 7.80x plasma (measured 0.61x) and methotrexate at
    #: 0.057x (measured 0.60x), i.e. wrong in OPPOSITE directions, so no single
    #: constant can absorb it (``results/breast_kp_validation.md``).
    #:
    #: With this on, the same permeability field is subtracted as a first-order
    #: loss, giving ``C_ss = PS*Cp/(PS + clearance + binding) <= Cp``: bounded
    #: by the plasma level by construction. No new parameter is introduced --
    #: it is the SAME ``PS`` the source already uses, applied to both halves of
    #: the gradient instead of only one.
    #:
    #: ``False`` (default) preserves the historical behaviour BIT-FOR-BIT so
    #: every stored reference result is unaffected. This defaults off because
    #: turning it on changes every absolute delivery number this project has
    #: ever published, which is a re-validation exercise, not a flag flip.
    #:
    #: NOTE this fixes the FORM of the delivery law, not the VALUE of the
    #: permeability. ``PS`` is still driven by ``bbb_penetrance``, a
    #: blood-BRAIN-barrier term that is simply wrong for non-brain tissue.
    #: Both defects were diagnosed together and are independent.
    #:
    #: NUMERICS. Verified against the closed-form steady state
    #: ``C_ss = PS*Cp/(PS + clearance + binding)`` to 6.3e-15 worst case across
    #: three orders of magnitude in PS, with the discrete mass budget closing
    #: to round-off. PS enters ``K`` inside the IMPLICIT operator
    #: ``A = I - dt*L + dt*K``, so it adds a non-negative diagonal to an
    #: already-M-matrix and imposes NO new timestep limit: solutions are
    #: invariant to 1e-4 even at ``dt*PS = 1245``, the advective CFL trip point
    #: is bit-identical, and the solver's zero-clamp was observed never to fire.
    #:
    #: RESOLUTION. The term introduces a new perivascular length scale
    #: ``sqrt(D_eff/(PS + clearance + binding))`` inside perfused cells, which
    #: collapses to ~15-17 um for the most permeable drugs (``bbb_penetrance``
    #: >~ 0.6: lomustine, carmustine, procarbazine, ribociclib) -- barely one
    #: grid cell per e-fold at dx = 12. An earlier version of this note
    #: predicted from that geometry alone that dx = 12 would be inadequate and
    #: that finer grids would push values DOWN. **Measured, and both halves
    #: were wrong**: sweeping dx = 12/10/8 um on real GBM tissue moves those
    #: drugs' interstitial:plasma ratios UP by only 1.3-2.6 % (lomustine +1.9,
    #: carmustine +1.4, procarbazine +1.8, ribociclib +2.6, cisplatin +1.3),
    #: with a low-permeability doxorubicin control flat at -0.2 %. **dx = 12 um
    #: is adequate with this flag on.** Kills move +1.0 to +2.0 % and the
    #: cross-drug reordering this change causes persists unchanged at dx = 8.
    #: See ``results/transvascular_exchange_revalidation.md`` section 5.
    transvascular_exchange: bool = False


@dataclass
class OxygenConfig:
    """Reaction-diffusion parameters for the steady-state oxygen field.

    ``consumption`` is UNCALIBRATED against any measured pO2 (see
    ``tmesim/physics/oxygen.py`` module docstring and
    ``results/oxygen_sensitivity_band.md`` for the quantified, literature-
    bracketed uncertainty range this default sits inside, rather than at the
    center of). ``sink_form`` selects between two functional forms for how the
    sink depends on local oxygen -- see :func:`tmesim.physics.oxygen.solve_oxygen`.
    Both ``consumption`` and ``sink_form`` default to the ORIGINAL behavior so
    nothing existing changes silently; the alternative sink form and the
    consumption sensitivity band are opt-in / reported-not-adopted.
    """

    diffusivity: float = 2000.0  # um^2/s
    consumption: float = 0.02    # per-tumor-cell consumption rate (relative)
    supply: float = 0.05         # vascular supply rate (relative)
    hypoxia_threshold: float = 0.2  # oxygen below this => HYPOXIC state

    #: {"linear", "michaelis_menten"}. "linear" (default) is the ORIGINAL sink,
    #: ``consumption * rho_tum * O`` -- proportional to O everywhere, so it
    #: self-throttles as O falls and never dominates the (bounded) supply term.
    #: "michaelis_menten" is ``consumption * rho_tum * O / (O + mm_half_saturation)``
    #: (Grimes, Fletcher & Partridge 2014, PMID 26064525) -- saturates toward a
    #: CONSTANT rate as O falls, so it does not self-throttle. See
    #: ``tmesim.physics.oxygen.solve_oxygen`` and ``results/oxygen_sensitivity_band.md``
    #: for why this axis is NOT a minor correction on this project's real flagship
    #: reconstruction (hypoxic tumor fraction swings 0.0 -> 0.69-1.0 across the
    #: literature-plausible half-saturation range at the current ``consumption``).
    sink_form: str = "linear"

    #: Dimensionless Michaelis-Menten half-saturation constant K (only used when
    #: ``sink_form == "michaelis_menten"``). Default = Grimes' own literature
    #: base case, ``Km = 1 mmHg`` (their sensitivity ceiling is 10 mmHg), expressed
    #: against this project's own O=1 anchor (``oxygen.REFERENCE_PO2_MMHG = 35``
    #: mmHg, Erecinska & Silver 2001, PMID 11718758): ``K = 1 / 35 = 0.0286``.
    #: This is A defensible choice, not THE resolved value -- see
    #: ``results/oxygen_sensitivity_band.md``.
    mm_half_saturation: float = 1.0 / 35.0

    #: Newton-Raphson iteration cap for the Michaelis-Menten sink (unused by
    #: "linear", which is a single direct sparse solve). Typical convergence is
    #: a handful of iterations; see ``tmesim/physics/oxygen.py``.
    mm_max_iter: int = 50

    #: Newton-Raphson convergence tolerance (max absolute change in O between
    #: iterations) for the Michaelis-Menten sink.
    mm_tol: float = 1e-9


@dataclass
class PDConfig:
    """Pharmacodynamics: how accumulated exposure maps to cell kill.

    ``hypoxia_resistance`` (this GLOBAL default) and ``quiescent_resistance`` are
    still UNVALIDATED PRIORS with no citation (flagged in-repo at
    ``validation/kill_calibration.py:1954,1997``). A literature search
    (``results/sweep/state_multiplier_literature_check_*.md``) found real,
    GBM-cell-line, hypoxia-vs-normoxia temozolomide data (Musah-Eroje & Watson
    2019, PMID 30694423; Wang et al. 2017, PMID 28801626; Ahmed et al. 2018, PMID
    29492900) that converges on a much SMALLER hypoxia effect for TMZ
    specifically than 0.35 implies -- IC50 shifts of only ~1.1-1.4x, and one clean
    null result in the matched-2D-oxygen comparison closest to what this
    multiplier represents.

    **TMZ-specific recalibration (done).** That IC50-fold evidence has been
    converted, Hill-curve-consistently, into an equivalent kill-probability
    multiplier at this project's own documented operating point (temozolomide,
    ``results/dose_schedule.md``: median tumour-cell exposure/IC50 = 0.056 at 1x
    nominal dose, drug Hill coefficient h=1.2 -- the "not saturated" sub-clinical
    regime this project's real runs actually occupy). For a hypoxia-only IC50
    shift by factor f at fixed absolute exposure E (ratio r = E/IC50):

        R_equivalent(r, f, h) = (1 + r^h) / (f^h + r^h)

    which as r -> 0 (this project's regime) tends to the closed form
    ``R_equivalent ~ 1 / f^h``. Using a literature-representative f = 1.2 (within
    the reported 1.1-1.4x TMZ range; Wang et al. 2017) gives R_equivalent ~
    0.81 at r=0.056, essentially unchanged (0.806-0.811) across this project's
    whole documented sub-saturating dose band (r = 0.034-0.078, the 0.6x-1.4x
    "approximate clinical band" in ``dose_schedule.md``). The chosen value,
    **0.80**, is stored as a PER-DRUG override on temozolomide's ``DrugParams``
    (see ``DrugParams.hypoxia_resistance`` and ``DRUG_LIBRARY``), NOT here --
    this global default stays at the original, still-uncited 0.35 for cisplatin
    and doxorubicin, for which no comparably quantitative literature basis was
    found in that search (cisplatin: real, GBM-specific, statistically
    significant hypoxia effect confirmed in direction, magnitude unextracted;
    doxorubicin: a real 4-10x fold-change exists but is prostate tissue, not
    glioma, and a mismatched exposure protocol).

    For ``quiescent_resistance``, no quantitative fold-change was found anywhere
    in the literature for any of this project's three drugs -- only mechanistic
    support for the direction (<1). It has NOT been changed: inventing a
    precise number with no measurement behind it would be worse than leaving
    the honestly-labeled unvalidated prior in place.
    """

    hypoxia_resistance: float = 0.35  # multiplier on kill for hypoxic cells (<1)
    quiescent_resistance: float = 0.5
    max_kill_rate: float = 1.0        # scaling on the Hill kill term


@dataclass
class SimConfig:
    """Top-level knobs bundling the per-stage configs."""

    drug: DrugParams = field(default_factory=lambda: DrugParams(name="temozolomide-like"))
    pk: PKConfig = field(default_factory=PKConfig)
    pde: PDEConfig = field(default_factory=PDEConfig)
    oxygen: OxygenConfig = field(default_factory=OxygenConfig)
    pd: PDConfig = field(default_factory=PDConfig)
    seed: int = 0


#: TMZ's literature-informed hypoxia-resistance override -- see
#: ``PDConfig.hypoxia_resistance`` docstring for the full operating-point
#: derivation (R_equivalent = (1+r^h)/(f^h+r^h), f~1.2 from Wang et al. 2017
#: PMID 28801626, r=0.056 from ``results/dose_schedule.md``). Applied both to
#: the synthetic "temozolomide-like" library entry below AND to the real,
#: GDSC-calibrated molecular-compile pathway
#: (``validation.calibration._compile``), which is what every real-data report
#: (``kill_calibration.md``, ``real_drug_panel.md``, ``dose_schedule.md``,
#: ``cell_state_impact.md``) actually uses for temozolomide.
TMZ_HYPOXIA_RESISTANCE: float = 0.80

# A default glioblastoma-relevant drug library. Agents may extend this.
#
# PROVENANCE NOTE (read before adding/consuming an entry). The first three
# entries below (temozolomide-like, doxorubicin-like, cisplatin-like) are the
# ORIGINAL Stage-0 entries: biologically-plausible, HAND-SET numbers, not
# computed from any molecular pipeline or GDSC data, and their ``binding`` /
# ``clearance`` predate the perivascular-sink calibration in
# ``validation.calibration`` (they carry the Stage-0 ~1e-4 1/s loss rate that
# ``validation.pk_reality`` showed over-predicts penetration by ~15x for
# doxorubicin -- see that module's docstring). They are kept exactly as they
# were; nothing about them changes here.
#
# The entries added below them (lomustine-like ... irinotecan-like) are
# DIFFERENT IN KIND: each is a FROZEN SNAPSHOT of
# ``validation.calibration.calibrated_drugparams(name, alpha)`` --
# i.e. real RDKit/Stokes-Einstein/Brownian-dynamics structure -> transport
# (``molecular.compile_drug``), real GDSC2-glioma-measured or curated-
# literature potency (``molecular.potency.resolve_ic50`` -- see that module for
# exact sourcing/citations per drug), and the SAME single doxorubicin-fit
# perivascular-sink constant every real-data report in this repo uses
# (``alpha = 148.608...``, frozen at compile time below; re-derive with
# ``validation.calibration.calibrate_alpha()`` if the doxorubicin anchor is
# ever revised -- these constants would then need refreezing to match).
# ``hill=1.2`` for every one of them is not a per-drug fit: it is
# ``molecular.compile_drug._DEFAULT_HILL``, the pipeline's constant, carried
# through unchanged -- an honest reflection of what the real pipeline actually
# outputs, not an invented per-drug shape parameter.
#
# CONSEQUENCE OF MIXING THE TWO STYLES IN ONE DICT: the new entries' loss
# rates (``binding``/``clearance``) are ~30-300x FASTER than the original
# three's Stage-0 defaults, because they carry the calibration fix the
# original three never received. A cross-drug comparison built directly on
# this dict (e.g. iterating ``DRUG_LIBRARY.values()``) is therefore NOT an
# apples-to-apples potency/penetration comparison across ALL ten entries --
# only within each style group. For a fair, validated cross-drug comparison,
# do what ``experiments/real_drug_panel.py`` already does for every drug
# (including the original three, by name): call
# ``validation.calibration.calibrated_drugparams`` directly rather than
# reading this dict. ``DRUG_LIBRARY`` remains useful as a quick, named,
# frozen-in-time reference (e.g. for ``treatment.drugs.get_drug`` / tests /
# simple demos) where reproducibility of a fixed snapshot matters more than
# picking up upstream recalibration automatically.
#
# ``bevacizumab`` (anti-VEGF IgG, ~149 kDa) was evaluated for this expansion
# and DELIBERATELY EXCLUDED, for a reason distinct from "no data found": it is
# not directly cytotoxic (no cell-kill IC50 exists for it in a kill-fraction
# sense -- GDSC/any cytotoxicity screen has no entry for it, and none should
# be expected, since it acts by VEGF sequestration/anti-angiogenesis, not by
# killing tumour cells at a target concentration). Forcing a number into
# ``ic50`` would misrepresent its mechanism inside a Hill-curve kill model
# built for cytotoxics. It remains available and useful elsewhere in this
# repo as the delivery-limited/antibody size-extreme demonstration
# (``molecular.library.GBM_DRUGS["bevacizumab"]``, ``results/molecular_findings.md``:
# effective_D ~1 um^2/s, bbb_penetrance ~0.000, "a large macromolecule
# essentially cannot reach the tumor by passive transport") -- that is the
# correct home for it, not this kill-fraction drug library.
DRUG_LIBRARY: Dict[str, DrugParams] = {
    "temozolomide-like": DrugParams(
        name="temozolomide-like", diffusivity=300.0, molecular_weight=194.0,
        ic50=25.0, hill=1.2, bbb_penetrance=0.4, clearance=3e-4,
        hypoxia_resistance=TMZ_HYPOXIA_RESISTANCE),
    "doxorubicin-like": DrugParams(
        name="doxorubicin-like", diffusivity=150.0, molecular_weight=543.0,
        ic50=0.5, hill=1.5, bbb_penetrance=0.05, binding=2e-3, clearance=1e-4),
    "cisplatin-like": DrugParams(
        name="cisplatin-like", diffusivity=350.0, molecular_weight=300.0,
        ic50=8.0, hill=1.3, bbb_penetrance=0.15, clearance=2e-4),

    # --- Real molecular-compile + GDSC/literature-calibrated entries -------
    # Frozen snapshots of calibrated_drugparams(name, alpha=148.60815765445545)
    # (alpha fit once on doxorubicin's e-fold penetration length; see
    # validation.calibration). Potency provenance (GDSC2 glioma-measured vs.
    # curated literature, with citations) is in molecular.potency; run
    # ``python -m tmesim.molecular.potency`` or see molecular/potency.py's
    # LITERATURE_IC50 table for the exact sources.

    # Nitrosourea (CCNU), oral GBM alkylator (PCV / lomustine regimens).
    # IC50: literature order-of-magnitude anchor (50 uM); real comparative
    # potency vs TMZ/ACNU across GBM lines in Yamamuro S et al., Cancer Sci
    # 2021;112(11):4736-47, PMID 34536314 (numeric values not independently
    # extracted from the paywalled figure -- see molecular/potency.py note).
    "lomustine-like": DrugParams(
        name="lomustine-like", diffusivity=259.2, molecular_weight=233.7,
        ic50=50.0, hill=1.2, bbb_penetrance=0.973, binding=1.941e-2,
        clearance=1.086e-4),

    # Nitrosourea (BCNU), lipophilic alkylator (Gliadel wafers). IC50:
    # literature order-of-magnitude anchor (30 uM); real BCNU/ACNU
    # cytotoxicity (apoptosis + DNA double-strand breaks) in GBM lines U87MG/
    # U138MG confirmed in Batista LFZ et al., Cancer Res 2007;67(24):11886-95,
    # PMID 18089819 (no numeric IC50 in that abstract; see molecular/potency.py
    # note on the wide 30-615 uM spread found across sources for this drug).
    "carmustine-like": DrugParams(
        name="carmustine-like", diffusivity=268.3, molecular_weight=214.1,
        ic50=30.0, hill=1.2, bbb_penetrance=0.955, binding=1.038e-2,
        clearance=9.298e-5),

    # Methylating pro-drug requiring metabolic activation (PCV regimen).
    # IC50: 500 uM, the WEAKEST-SOURCED literature entry in this library --
    # no citable glioma-cell-line IC50 with a verifiable PMID was found; this
    # is a mechanism-based order-of-magnitude guess (see molecular/potency.py
    # note, PMID 1596582 for the activation mechanism in a non-glioma system).
    "procarbazine-like": DrugParams(
        name="procarbazine-like", diffusivity=264.7, molecular_weight=221.3,
        ic50=500.0, hill=1.2, bbb_penetrance=0.885, binding=1.998e-2,
        clearance=9.101e-5),

    # Vinca alkaloid (microtubule inhibitor); historically part of PCV.
    # IC50: 0.05 uM (50 nM), a GENERAL vinca-alkaloid potency anchor, not a
    # glioma-specific measurement -- real GBM activity (U87/U251) confirmed in
    # Filippi-Chiela EC et al., Toxicol In Vitro 2022;85:105472, PMID 36116745,
    # but that paper reports no numeric IC50 (see molecular/potency.py note).
    # Structure-derived bbb_penetrance is correctly very low (0.023): vincristine
    # is large (MW 825) and does not cross the BBB well.
    "vincristine-like": DrugParams(
        name="vincristine-like", diffusivity=155.0, molecular_weight=825.0,
        ic50=0.05, hill=1.2, bbb_penetrance=0.023, binding=8.906e-2,
        clearance=8.938e-5),

    # Antifolate antimetabolite. IC50 12.0 uM is REAL, MEASURED GDSC2 glioma
    # data (median across n=49 glioma cell lines; molecular.potency). Note the
    # structure-derived bbb_penetrance is correctly very low (0.0225): methotrexate
    # is large, highly polar (TPSA ~211, 5 H-bond donors) and is well known
    # clinically to cross the BBB poorly (it is dosed intrathecally for CNS
    # disease) -- the model reproduces this qualitative fact from structure alone.
    "methotrexate-like": DrugParams(
        name="methotrexate-like", diffusivity=198.9, molecular_weight=454.4,
        ic50=12.0, hill=1.2, bbb_penetrance=0.0225, binding=2.663e-2,
        clearance=6.946e-5),

    # Topoisomerase-II poison (VP-16). IC50 5.0 uM is a representative point
    # inside a REAL measured range: 48h IC50 spans 0.6-29 uM across GBM lines
    # M059K/M059J/T98G depending on the tumour's own TOP2A expression (Sevim H
    # et al., J Cancer Res Clin Oncol 2011;137(11):1705-12, PMID 21904904) --
    # this drug's real glioma potency is a ~48x range, not a point; see
    # molecular/potency.py for the full note.
    "etoposide-like": DrugParams(
        name="etoposide-like", diffusivity=178.6, molecular_weight=588.6,
        ic50=5.0, hill=1.2, bbb_penetrance=0.039, binding=2.385e-2,
        clearance=7.484e-5),

    # Topoisomerase-I inhibitor PRO-DRUG (camptothecin analogue; clinically
    # combined with bevacizumab for recurrent GBM). IC50 11.58 uM is REAL,
    # MEASURED GDSC2 glioma data (median across n=50 glioma cell lines).
    # HONESTY: irinotecan is cleaved by carboxylesterases to the ~1000x more
    # potent SN-38; this pipeline has no prodrug-metabolism step, so diffusivity
    # / bbb_penetrance describe the PARENT molecule (irinotecan), not SN-38 --
    # the GDSC assay doses/measures irinotecan directly, so its cytotoxicity
    # tier is unaffected. Structure-derived bbb_penetrance (0.151) lands close
    # to the REAL measured irinotecan CSF:plasma AUC ratio (14 +/- 3%, nonhuman
    # primates; Blaney S et al., Cancer Chemother Pharmacol 1998;41(6):464-8,
    # PMID 9554590) -- a genuine, if modest, external cross-check on this
    # project's structure-only BBB rule. See molecular/potency.py for the full note.
    "irinotecan-like": DrugParams(
        name="irinotecan-like", diffusivity=178.9, molecular_weight=586.7,
        ic50=11.58, hill=1.2, bbb_penetrance=0.151, binding=6.012e-2,
        clearance=1.141e-4),

    # CDK4/6 inhibitors. Added following Neiswender/Shibue/Vazquez et al.,
    # "A dependency map enhanced with next-generation 3D cancer models",
    # Nature 2026, DOI 10.1038/s41586-026-10843-7 -- their patient-derived GBM
    # organoid CRISPR + drug panel found CDK6 dependency, and CDK4/6-inhibitor
    # sensitivity, ENRICHED IN CDKN2A-DEFICIENT glioblastoma specifically (Fig.
    # 3g-i, Extended Data Fig. 6h-k; GR_AOC across 7 patient-derived models:
    # 5 CDKN2A-loss [BT145, BT179, BT286, BT320, BT444], 2 CDKN2A-neutral
    # [BT187, BT224]; abemaciclib P=7.1e-4, palbociclib P=0.021, ribociclib
    # P=0.051 (NS)).
    #
    # Verified against that paper's OWN deposited Source Data (Fig. 3h sheet,
    # Springer MOESM10), not just its prose: the per-model GR_AOC values
    # recompute to the published group means (loss vs neutral 0.358/0.058
    # abemaciclib, 0.252/-0.033 palbociclib, 0.154/-0.070 ribociclib). Three
    # things that recomputation surfaced which the prose does not say:
    #   (1) RIBOCICLIB's loss group is n=4, not 5 -- BT444 has no ribociclib
    #       data at all. Its P=0.051 is a 4-vs-2 comparison.
    #   (2) Several GR50s are RIGHT-CENSORED at the 10 uM top dose (BT320
    #       palbo; BT187 and BT224 palbo+ribo), so "10" there is a bound, not
    #       a fitted value.
    #   (3) Culture format is NOT balanced against CDKN2A status (4 of 5 loss
    #       models are 3-D spheroid/ULA; the 2 neutral models split 1 spheroid
    #       / 1 adherent-Matrigel) -- i.e. a geometry/transport variable is
    #       entangled with the genotype contrast. Checked directly rather than
    #       assumed: the effect holds WITHIN both format strata and in the
    #       same direction for every drug, and the entanglement runs AGAINST
    #       the finding (the mostly-spheroid loss group is the MORE sensitive
    #       one, whereas diffusion limitation would make 3-D look LESS
    #       sensitive). So the confound does not manufacture the result; it
    #       makes it conservative. Recorded because it is exactly the kind of
    #       transport-vs-genotype entanglement this simulator exists to
    #       separate.
    #
    # HONESTY on what the single ic50 below can and can't represent: it is the
    # REAL, MEASURED GDSC2 glioma IC50 (median across n=51 GBM/LGG cell lines
    # per drug; tmesim.molecular.potency), an UNSELECTED population average
    # that does NOT stratify by CDKN2A status. The Nature 2026 finding above is
    # that potency is genotype-dependent -- CDKN2A-loss tumours are more
    # sensitive than this number implies, CDKN2A-neutral tumours less so. This
    # library, and the Hill kill model it feeds, has no per-tumour genotype
    # covariate, so a single ic50 necessarily collapses that real
    # heterogeneity; do not read a run's kill_fraction for these two drugs as
    # CDKN2A-status-aware.
    #
    # palbociclib's STRUCTURE is reused from ``molecular.breast_library.
    # BREAST_DRUGS["palbociclib"]`` (same real drug -- the two libraries are
    # kept disjoint on purpose, see that file and ``molecular.library``'s
    # comment at the CDK4/6 entry), compiled here with THIS library's own
    # glioma-context GDSC2 IC50 (52.19 uM) instead of the breast library's
    # GDSC2-breast IC50 (59.98 uM, see ``results/kill_calibration_breast.md``)
    # -- same molecule, genuinely different measured potency per tissue's own
    # cell-line panel, the same cross-library reuse pattern already
    # established for doxorubicin/cisplatin.
    "palbociclib-like": DrugParams(
        name="palbociclib-like", diffusivity=200.11076033801277,
        molecular_weight=447.543, ic50=52.186247378767696, hill=1.2,
        bbb_penetrance=0.578756414591833, binding=0.044354404560997165,
        clearance=0.00010525462510012649),
    "ribociclib-like": DrugParams(
        name="ribociclib-like", diffusivity=202.4826355721452,
        molecular_weight=434.548, ic50=50.042017373304844, hill=1.2,
        bbb_penetrance=0.7995183083447064, binding=0.03993226707215974,
        clearance=0.0001034244781996349),
    # abemaciclib is DELIBERATELY EXCLUDED despite being the strongest hit in
    # the Nature 2026 CDKN2A analysis (P=7.1e-4): no glioma cell-viability
    # IC50 (GDSC-measured or citable literature) was found for it -- only its
    # CDK4/CDK6 biochemical kinase-inhibition constants (Ki ~0.002/0.010 uM),
    # a different unit of measurement (enzyme affinity, not cell-kill potency)
    # that would misrepresent this drug's Hill-curve potency if forced in.
}
