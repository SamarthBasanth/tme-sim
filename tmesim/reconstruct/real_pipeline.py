"""Whole-niche REAL microenvironment assembler (the flagship reconstruction).

This is the single entry point that turns a REAL spatial-transcriptomics sample
into an assembled, physics-ready :class:`~tmesim.core.Microenvironment` with NO
synthetic data anywhere in the path (see ``docs/REAL_BUILD_PLAN.md``). It wires
together the Wave-2 real modules:

  * cells + phenotypes  — :func:`reconstruct.real_xenium.reconstruct_from_xenium`
    (single-cell Xenium, threshold-gated marker typing) or
    :func:`reconstruct.real_visium.reconstruct_from_visium` (spot deconvolution);
  * REAL vascular field — KDE of measured ENDOTHELIAL cells
    (:func:`reconstruct.real_fields.vascular_density_from_cells`);
  * REAL ECM field      — for Xenium, an EXPRESSION-WEIGHTED KDE of the panel's
    real collagen / matrix transcripts (COL1A2, POSTN, VCAN, DCN, LAMA2, THBS1,
    ...) via :func:`reconstruct.real_fields.ecm_density_from_positions`; for
    Visium (no per-cell transcripts here) the stromal-fraction proxy that
    ``reconstruct_from_visium`` already builds;
  * REAL-perfusion oxygen — the unmodified physics solver driven by the real
    vascular field (:func:`reconstruct.real_fields.add_oxygen`).

Every field derives from a measured byte on disk. The per-cell ``sensitivity``
carried through from the Xenium reconstruction is **uniform 1.0** — the geometry
is drug-agnostic and asserts no intrinsic per-cell resistance, so a kill fraction
computed on it is pure delivery + potency physics. Drug-specific sensitivity is
opt-in via :func:`apply_drug_sensitivity`, which uses a measured panel determinant
only where one genuinely exists (TOP2A for TOP2 poisons, MKI67 for cell-cycle-
dependent cytotoxics) and leaves sensitivity at 1.0 otherwise (temozolomide's
determinant, MGMT, is not on this panel). ABSOLUTE kill fractions are CONDITIONAL
on whichever assumption is in force; this section is an untreated resection with
no outcome label. Nothing here modifies core / engine / physics / molecular.

MGMT (temozolomide) — an OPT-IN per-patient escape from that uniform assumption.
:func:`apply_drug_sensitivity` accepts an optional ``mgmt=`` value (measured
expression, or a binary promoter status). MGMT is not a spatial field we can
resolve — it is on neither Xenium panel and has no panel surrogate — but it IS a
per-patient scalar, and the mapping from MGMT level to TMZ sensitivity is
calibrated on measured ex-vivo response (:mod:`tmesim.molecular.mgmt`;
``results/mgmt_scalar.md``). **Omitting ``mgmt`` reproduces the previous
behaviour exactly**; supplying it flips ``upper_bound`` to False and records the
value, so a TMZ result can no longer silently carry an UPPER BOUND label.
"""
from __future__ import annotations

import time
from typing import Dict, List, Optional

import numpy as np

from ..core import (
    CellPopulation,
    CellType,
    Microenvironment,
    OxygenConfig,
    TissueFields,
)
from ..data.real.singlecell_spatial import load as load_singlecell_spatial
from ..molecular.mgmt import (
    MGMT_CALIBRATED_DRUGS,
    is_mgmt_calibrated,
    mgmt_sensitivity,
)
from ..physics.oxygen_proxy import (
    OXYGEN_MODES,
    hypoxia_from_markers,
    oxygen_from_vessel_distance,
)
from .real_fields import (
    add_oxygen,
    ecm_density_from_cells,
    ecm_density_from_positions,
    vascular_density_from_cells,
)
from .real_visium import reconstruct_from_visium
from .real_xenium import (
    TISSUE_DISC_RADIUS_UM,
    TISSUE_MASK_MODE,
    per_cell_sensitivity,
    reconstruct_from_xenium,
)

# REAL extracellular-matrix / matricellular transcripts carried on the Xenium
# Human Brain panel. The ECM field is an expression-weighted KDE of these (only
# the ones actually present in a given panel are used at runtime).
ECM_MARKERS: List[str] = [
    "COL1A2", "COL12A1", "COL25A1", "POSTN", "VCAN",
    "DCN", "FBLN1", "LAMA2", "THBS1", "SPON1", "TGFBI",
]

_DEFAULT_VISIUM_SAMPLE = "10x_parent_gbm"

#: How the oxygen field is produced. **The default is and stays ``'pde'``** so no
#: existing result moves; ``tests/test_oxygen_proxy.py`` asserts the default path
#: is BIT-IDENTICAL to the pre-option pipeline.
#:
#: ``'pde'`` (DEFAULT, unchanged)
#:     ``physics.solve_oxygen`` — the steady-state reaction-diffusion solve driven
#:     by the real vascular KDE. **Falsified as value-add**
#:     (``results/state_validation.md`` §1.2, ``results/oxygen_demotion.md``): it
#:     is out-predicted by a distance transform built from its own source cells,
#:     and by its own undiffused source. Kept as the default only because flipping
#:     it silently would invalidate every published number in ``results/``.
#: ``'vessel_distance'``
#:     :func:`physics.oxygen_proxy.oxygen_from_vessel_distance` — a monotone
#:     SURROGATE from distance to the nearest endothelial cell. Better against
#:     measured markers, ~10^2-10^3x cheaper, and honestly labelled as not being a
#:     solved oxygen tension.
#: ``'markers'``
#:     :func:`physics.oxygen_proxy.hypoxia_from_markers` — hypoxia read directly
#:     off MEASURED hypoxia-response transcripts. Only available where the panel
#:     carries them; on the 319-gene brain panel it raises with an explicit reason
#:     rather than silently degrading, because on that section hypoxia is
#:     UNMEASURED. Xenium only (Visium has no per-cell transcripts here).
OXYGEN_MODE_DEFAULT = "pde"


def _load_xenium_retry(sample: Optional[str]):
    """Load the real Xenium data, retrying once on any transient/IO error."""
    kwargs = {"dense": False}
    if sample is not None:
        kwargs["path"] = sample
    try:
        return load_singlecell_spatial(**kwargs)
    except Exception:
        time.sleep(3)
        return load_singlecell_spatial(**kwargs)


def _ecm_weights_from_expression(expression, genes: np.ndarray):
    """Per-cell REAL ECM abundance = summed log1p counts of the panel's matrix
    genes. Returns ``(weights (n,), used_markers)``; weights are all-zero (and the
    caller should fall back to the stromal proxy) if no ECM marker is on the panel.
    """
    genes = np.asarray(genes, dtype=str)
    gene_index = {g: j for j, g in enumerate(genes)}
    used = [g for g in ECM_MARKERS if g in gene_index]
    if not used:
        return np.zeros(expression.shape[0], dtype=np.float64), used
    cols = [gene_index[g] for g in used]
    from scipy import sparse

    if sparse.issparse(expression):
        sub = expression[:, cols].toarray()
    else:
        sub = np.asarray(expression)[:, cols]
    weights = np.log1p(sub.astype(np.float64)).sum(axis=1)
    return weights, used


def _install_oxygen(
    micro: Microenvironment,
    oxygen_mode: str,
    cfg: OxygenConfig,
    expression=None,
    genes=None,
    oxygen_kwargs: Optional[Dict] = None,
) -> Dict:
    """Fill ``micro.fields.oxygen`` by the requested route and describe what it is.

    ADDITIVE and DEFAULT-PRESERVING: with ``oxygen_mode='pde'`` (the default) the
    only call made is the historical ``add_oxygen(micro, cfg)``, so the field is
    bit-identical to the pre-option pipeline. Returns a provenance block that
    :func:`reconstruction_report` surfaces, so an alternative oxygen route can
    never be carried silently.

    Raises ``ValueError`` for an unknown mode, and for ``'markers'`` on a panel
    with no hypoxia-response gene — the honest failure, since on such a section
    hypoxia is UNMEASURED (see
    :func:`physics.oxygen_proxy.hypoxia_from_markers`).
    """
    mode = str(oxygen_mode).strip().lower()
    kw = dict(oxygen_kwargs or {})
    if mode not in OXYGEN_MODES:
        raise ValueError(
            f"unknown oxygen_mode {oxygen_mode!r}; expected one of {OXYGEN_MODES}"
        )

    if mode == "pde":
        if kw:
            raise TypeError(
                f"unexpected oxygen_kwargs {sorted(kw)} for oxygen_mode='pde'")
        add_oxygen(micro, cfg)
        return {
            "oxygen_mode": "pde",
            "oxygen_source": "physics.solve_oxygen driven by the real vascular field",
            "oxygen_is_surrogate": False,
            "oxygen_caveat": (
                "FALSIFIED AS VALUE-ADD (results/state_validation.md §1.2, "
                "results/oxygen_demotion.md): out-predicted against measured "
                "hypoxia markers by a distance transform on its own source cells "
                "(rho -0.1503 vs -0.1132) and by its own undiffused source "
                "(-0.1270). Retained as the default so no published number moves "
                "silently, NOT because it is the best available route."
            ),
            "oxygen_detail": None,
        }

    if mode == "vessel_distance":
        field, info = oxygen_from_vessel_distance(micro, cfg, detail=True, **kw)
        micro.fields.oxygen = field
        return {
            "oxygen_mode": "vessel_distance",
            "oxygen_source": (
                "SURROGATE: exp(-d/L) in distance to the nearest "
                f"{info['source']} vessel point (L = "
                f"{info['diffusion_length_um']:.2f} um)"
            ),
            "oxygen_is_surrogate": True,
            "oxygen_caveat": (
                "NOT a solved oxygen tension. A monotone geometric proxy, "
                "calibrated to the Thomlinson & Gray 1955 oxygen-diffusion "
                "distance (PMID 13304213) so its hypoxia band edge sits at 150 um "
                "from a vessel. Its ABSOLUTE level is uncalibrated — as is the "
                "PDE's."
            ),
            "oxygen_detail": info,
        }

    # mode == 'markers'
    if expression is None or genes is None:
        raise ValueError(
            "oxygen_mode='markers' needs per-cell expression; it is only "
            "available for source='xenium' (Visium carries no per-cell transcripts)"
        )
    detail = hypoxia_from_markers(micro, expression, genes, detail=True, **kw)
    if detail["field"] is None:
        raise ValueError(f"oxygen_mode='markers' is unavailable: {detail['reason']}")
    micro.fields.oxygen = detail["field"]
    slim = {k: v for k, v in detail.items() if k not in ("field", "per_cell")}
    return {
        "oxygen_mode": "markers",
        "oxygen_source": (
            "MEASURED hypoxia-response transcripts "
            f"({', '.join(detail['genes_used'])}), kernel-regressed onto the grid "
            "and rank-inverted to an oxygen-like [0, 1] field"
        ),
        "oxygen_is_surrogate": True,
        "oxygen_caveat": (
            "MEASURED, but a hypoxia-RESPONSE readout, not an oxygen tension. "
            "Rank-normalised, so its absolute level and its band fractions are "
            "artefacts of the rank map, not measurements. Validating a state model "
            "built on this against these same markers would be CIRCULAR."
        ),
        "oxygen_detail": slim,
    }


def _reconstruct_xenium(
    sample: Optional[str],
    dx_um: float,
    tissue_mask_mode: str = TISSUE_MASK_MODE,
    tissue_disc_radius_um: float = TISSUE_DISC_RADIUS_UM,
    oxygen_mode: str = OXYGEN_MODE_DEFAULT,
    oxygen_kwargs: Optional[Dict] = None,
) -> Microenvironment:
    """Assemble the real single-cell Xenium microenvironment with REAL fields.

    Cells + threshold-gated phenotypes + domain come from
    :func:`reconstruct_from_xenium` (the canonical, fixed-typing cell
    reconstructor). We then rebuild the physical fields as REAL fields: vascular
    from the measured endothelial cells, ECM from the real matrix-transcript KDE,
    and oxygen from the physics solver driven by that real vascular field.
    """
    micro = reconstruct_from_xenium(
        sample_path=sample, dx_um=dx_um,
        tissue_mask_mode=tissue_mask_mode,
        tissue_disc_radius_um=tissue_disc_radius_um,
    )
    domain, cells = micro.domain, micro.cells

    # REAL vascular density: KDE of the measured ENDOTHELIAL cell positions.
    vascular = vascular_density_from_cells(domain, cells)

    # REAL ECM: expression-weighted KDE of the panel's collagen / matrix genes.
    # Load once more (aligned row-for-row to the deterministic reconstruction) to
    # read the per-cell matrix-transcript abundance the CellPopulation drops.
    data = _load_xenium_retry(sample)
    weights, ecm_markers = _ecm_weights_from_expression(
        data["expression"], data["genes"]
    )
    if weights.size == cells.n and float(weights.sum()) > 0.0:
        ecm = ecm_density_from_positions(domain, cells.x, cells.y, weights=weights)
        ecm_source = (
            "expression-weighted KDE of REAL matrix transcripts "
            f"({', '.join(ecm_markers)})"
        )
    else:
        # Panel carries no ECM transcript (or row mismatch): stromal-cell proxy.
        ecm = ecm_density_from_cells(domain, cells)
        ecm_source = "STROMAL-cell density proxy (no ECM transcript available)"

    micro.fields = TissueFields(
        ecm_density=ecm, vascular_density=vascular, oxygen=None
    )
    oxy_info = _install_oxygen(
        micro, oxygen_mode, OxygenConfig(),
        expression=data["expression"], genes=data["genes"],
        oxygen_kwargs=oxygen_kwargs,
    )

    micro.metadata = dict(micro.metadata)
    micro.metadata["real_pipeline"] = {
        "source": "xenium",
        "vascular_source": "KDE of REAL ENDOTHELIAL cell positions",
        "ecm_source": ecm_source,
        "ecm_markers_used": ecm_markers,
        "tissue_mask_mode": str(tissue_mask_mode),
        "tissue_disc_radius_um": (
            float(tissue_disc_radius_um)
            if str(tissue_mask_mode).lower() == "geometric" else None
        ),
        **oxy_info,
    }
    return micro


def _reconstruct_visium(
    sample: Optional[str],
    dx_um: float,
    oxygen_mode: str = OXYGEN_MODE_DEFAULT,
    oxygen_kwargs: Optional[Dict] = None,
) -> Microenvironment:
    """Assemble the real spot-level Visium microenvironment and add oxygen.

    ``reconstruct_from_visium`` already builds vascular density (from the real
    ENDOTHELIAL spot fraction) and an ECM proxy (from the real STROMAL fraction);
    per-cell collagen transcripts are not available at spot resolution, so the
    stromal proxy is the honest ECM here. We only add the real-perfusion oxygen.
    """
    sample_id = sample or _DEFAULT_VISIUM_SAMPLE
    # Visium spot pitch is ~100 um; a fine (xenium-style) dx would over-resolve
    # the spots, so fall back to the loader default when dx is left small.
    visium_dx = float(dx_um) if dx_um and dx_um >= 50.0 else 100.0
    micro = reconstruct_from_visium(sample_id, dx_um=visium_dx)
    mode = str(oxygen_mode).strip().lower()
    if mode == "pde":
        # UNCHANGED path: only solve if the Visium reconstructor left oxygen unset.
        oxy_info = {
            "oxygen_mode": "pde",
            "oxygen_source": "physics.solve_oxygen driven by the real vascular field",
            "oxygen_is_surrogate": False,
            "oxygen_caveat": None,
            "oxygen_detail": None,
        }
        if micro.fields.oxygen is None:
            oxy_info = _install_oxygen(micro, "pde", OxygenConfig())
    else:
        # A non-default mode always (re)installs, so the caller gets what they asked
        # for rather than whatever the Visium reconstructor happened to leave behind.
        oxy_info = _install_oxygen(
            micro, mode, OxygenConfig(), oxygen_kwargs=oxygen_kwargs
        )
    micro.metadata = dict(micro.metadata)
    micro.metadata["real_pipeline"] = {
        "source": "visium",
        "vascular_source": "KDE of REAL ENDOTHELIAL spot fraction",
        "ecm_source": "STROMAL-fraction proxy (spot resolution; no per-cell ECM transcript)",
        "ecm_markers_used": [],
        **oxy_info,
    }
    return micro


def reconstruct_real(
    source: str = "xenium",
    sample: Optional[str] = None,
    dx_um: float = 30.0,
    tissue_mask_mode: str = TISSUE_MASK_MODE,
    tissue_disc_radius_um: float = TISSUE_DISC_RADIUS_UM,
    oxygen_mode: str = OXYGEN_MODE_DEFAULT,
    oxygen_kwargs: Optional[Dict] = None,
) -> Microenvironment:
    """Reconstruct a whole-niche REAL microenvironment from open spatial data.

    Parameters
    ----------
    source : {'xenium', 'visium'}
        'xenium' -> single-cell-resolution 10x Xenium GBM section (threshold-gated
        marker typing, expression-weighted ECM). 'visium' -> spot-level 10x Visium
        GBM (reference-based deconvolution, stromal-proxy ECM).
    sample : str or None
        Xenium: cache directory (None -> loader default). Visium: sample id (None
        -> ``'10x_parent_gbm'``).
    dx_um : float
        Continuum grid spacing (um). Default 30 um suits the single-cell Xenium
        section; for Visium a value < 50 falls back to the 100 um spot pitch.
    tissue_mask_mode : {'geometric', 'occupancy'}
        Xenium only. ``'geometric'`` (default) builds ``Domain.tissue_mask`` as
        the union of ``tissue_disc_radius_um`` discs around the measured cell
        centroids, so the simulated tissue domain — and therefore the total O2
        supply the vascular KDE is masked to — does NOT depend on ``dx_um``.
        ``'occupancy'`` reproduces the legacy grid-dependent rule for comparison
        (see :func:`reconstruct.real_xenium.tissue_mask_from_cells`).
    tissue_disc_radius_um : float
        Disc radius for the geometric mask (default 30 um; justified against the
        measured 16.7 um median nearest-neighbour cell spacing).
    oxygen_mode : {'pde', 'vessel_distance', 'markers'}
        **ADDITIVE, OPT-IN, DEFAULT UNCHANGED.** How ``fields.oxygen`` is built —
        see :data:`OXYGEN_MODE_DEFAULT` for what each route is and what it is
        worth. The default ``'pde'`` runs exactly the historical call and is
        BIT-IDENTICAL to the pre-option pipeline (asserted by
        ``tests/test_oxygen_proxy.py::test_default_oxygen_mode_is_bit_identical``),
        so no existing result in ``results/`` moves.

        ``'pde'`` is retained as the default despite being falsified as value-add
        (``results/oxygen_demotion.md``) because flipping it would silently
        invalidate every number that has already been published against it. The
        recommendation and the re-run list are in that report.
    oxygen_kwargs : dict, optional
        Extra keyword arguments forwarded to the chosen route
        (:func:`physics.oxygen_proxy.oxygen_from_vessel_distance` or
        :func:`~physics.oxygen_proxy.hypoxia_from_markers`) — e.g.
        ``{'source': 'vessel_network'}`` to measure distance to a SEGMENTED vessel
        rather than to any endothelial cell. Must be empty (or ``None``) for
        ``oxygen_mode='pde'``, which takes no options.

    Returns
    -------
    Microenvironment
        Real cells at real positions, REAL vascular / ECM / oxygen fields, ready
        for :func:`tmesim.sim.engine.simulate_response`.

        NOTE for callers using a non-default ``oxygen_mode``:
        :func:`tmesim.sim.engine.simulate_response` RE-SOLVES the PDE internally
        and overwrites ``fields.oxygen``, so passing a surrogate reconstruction
        straight into the engine does NOT change its kill. Oxygen reaches the kill
        only through ``cells.agents.update_cell_states`` -> the ``PDConfig``
        resistance multipliers, so the way to propagate a surrogate is to re-band
        the states and re-run ``cells.pkpd.apply_drug_kill`` on the engine's own
        transport result. ``validation.oxygen_benchmark`` does exactly that and
        verifies it reproduces the engine bit-for-bit in the ``'pde'`` arm.
    """
    src = str(source).strip().lower()
    if src == "xenium":
        return _reconstruct_xenium(
            sample, dx_um, tissue_mask_mode, tissue_disc_radius_um,
            oxygen_mode=oxygen_mode, oxygen_kwargs=oxygen_kwargs,
        )
    if src == "visium":
        return _reconstruct_visium(
            sample, dx_um, oxygen_mode=oxygen_mode, oxygen_kwargs=oxygen_kwargs
        )
    raise ValueError(f"unknown source {source!r}; expected 'xenium' or 'visium'")


#: The two per-cell-sensitivity modes, and what each one is FOR.
#:
#: ``'uniform'`` — every drug gets sensitivity ``1.0`` for every cell. **This is
#:   the default, and it is the only FAIR mode for a CROSS-DRUG comparison**: it
#:   isolates DELIVERY x POTENCY with no sensitivity confound, because the
#:   sensitivity factor is identical for every drug and therefore cancels out of
#:   any ranking. Absolute kills are then a no-intrinsic-resistance upper bound
#:   for every drug equally.
#:
#: ``'determinant'`` — each drug uses its own measured panel determinant where
#:   one exists (:data:`reconstruct.real_xenium.SENSITIVITY_DETERMINANTS`:
#:   TOP2A for TOP2 poisons, MKI67 for cell-cycle-dependent cytotoxics) and
#:   uniform 1.0 where none does. Legitimate for a WITHIN-drug question ("where
#:   on this section does doxorubicin work?"), but **NOT comparable ACROSS
#:   drugs**: a drug whose determinant is on the panel is scaled down (mean
#:   sensitivity ~0.585 on this section) while a drug whose determinant is
#:   absent keeps 1.0, so the unmeasured drug is flattered by a factor ~1.7 for
#:   a purely technical reason. Temozolomide is exactly that case — its
#:   determinant MGMT is not on the 319-gene panel — so in this mode TMZ is an
#:   UPPER BOUND while its competitors are not, and every such drug is flagged
#:   ``upper_bound=True`` in the returned ``info``.
SENSITIVITY_MODES = ("uniform", "determinant")

_UPPER_BOUND_FLAG = (
    "UPPER BOUND — no intrinsic resistance is modelled for this drug (its "
    "determinant is not measurable on this panel), while other drugs in the same "
    "comparison ARE scaled down by theirs. Do NOT rank this drug against them in "
    "'determinant' mode; use mode='uniform' for cross-drug ranking."
)


def _apply_mgmt_scalar(
    micro: Microenvironment,
    info: Dict,
    mgmt,
    mgmt_units: str,
    mgmt_model: str,
) -> Dict:
    """Scale the per-cell sensitivity by a MEASURED per-patient MGMT scalar.

    Only called when the caller actually supplies ``mgmt``; with ``mgmt=None``
    this function is never reached and behaviour is bit-identical to before.

    The scalar is a **per-patient**, not a per-cell, quantity: MGMT is on neither
    Xenium panel and no panel gene surrogates for it
    (``results/gbm_tmz_exvivo.md`` §B2), while Visium recovers it only as
    pseudobulk. So it multiplies the whole sensitivity array uniformly, and the
    returned ``info`` says so.
    """
    name = info["drug"]
    if not is_mgmt_calibrated(name):
        info["mgmt_supplied"] = True
        info["mgmt_applied"] = False
        info["mgmt_note"] = (
            f"an MGMT value was supplied but {name!r} is not in "
            f"MGMT_CALIBRATED_DRUGS {MGMT_CALIBRATED_DRUGS}: the calibration is a "
            "temozolomide dose-response (GSE232173, n = 17) and is NOT "
            "extrapolated to other alkylators. Sensitivity was left unchanged."
        )
        return info

    if isinstance(mgmt, str):
        det = mgmt_sensitivity(status=mgmt, detail=True)
    else:
        det = mgmt_sensitivity(float(mgmt), units=mgmt_units, model=mgmt_model,
                               detail=True)
    scalar = float(det["sensitivity"])

    base = (np.ones(micro.cells.n, dtype=np.float64)
            if info.get("is_uniform", False)
            else np.asarray(micro.cells.sensitivity, dtype=np.float64))
    sens = np.clip(base * scalar, 0.0, 1.0)
    micro.cells.sensitivity = sens

    lam = det.get("lambda")
    mapping = (
        f"s = exp(-{lam:.4f} * log2(CPM+1))" if lam is not None
        else "measured median-kill ratio between the MGMT-expressing and "
             "MGMT-silent culture groups (Mann-Whitney p = 0.015)"
    )
    info["determinant_gene"] = "MGMT"
    info["kind"] = ("measured per-patient MGMT scalar (calibrated on measured "
                    "ex-vivo TMZ response, GSE232173 n = 17)")
    info["rationale"] = (
        "MGMT reverses the O6-methylguanine lesion temozolomide depends on "
        "(Hegi PMID 15758010, Esteller PMID 11070098). It is absent from both "
        "Xenium panels and has no panel surrogate, but it is a PER-PATIENT "
        "SCALAR, not a spatial field — so a measured MGMT value (clinical assay, "
        f"or Visium pseudobulk) scales TMZ sensitivity directly. Mapping: {mapping}, "
        "anchored at s = 1.0 for MGMT-silent, calibrated on measured ex-vivo "
        "% viability at 100 uM TMZ in 17 paired patient-derived GSC cultures."
    )
    info["mgmt_supplied"] = True
    info["mgmt_applied"] = True
    info["mgmt_input"] = det["input"]
    info["mgmt_model"] = det["model"]
    info["mgmt_scalar"] = scalar
    info["mgmt_scalar_ci95"] = det.get("ci95")
    info["mgmt_calibration_n"] = det.get("n")
    info["mgmt_extrapolated"] = bool(det.get("extrapolated", False))
    info["mgmt_resolution"] = (
        "PER-PATIENT SCALAR (uniform over cells). MGMT is not resolvable "
        "spatially on either Xenium panel or at Visium spot depth; no per-cell "
        "MGMT heterogeneity is asserted."
    )
    info["mean"] = float(sens.mean())
    info["std"] = float(sens.std())
    info["min"] = float(sens.min())
    info["max"] = float(sens.max())
    info["is_uniform"] = bool(np.allclose(sens, sens.flat[0])) if sens.size else True
    info["upper_bound"] = False
    info["comparability"] = (
        f"MGMT-INFORMED (scalar {scalar:.4f}). This drug's resistance determinant "
        "is now a MEASURED input, so its kill is no longer a "
        "no-intrinsic-resistance upper bound. It remains NOT directly comparable "
        "with drugs still carrying uniform 1.0 (those are still upper bounds) — "
        "compare like with like, and see results/mgmt_scalar.md for the "
        "calibration, its n = 17, and its wide CI."
    )
    return info


def apply_drug_sensitivity(
    micro: Microenvironment,
    drug_name: str,
    sample: Optional[str] = None,
    mode: str = "uniform",
    mgmt=None,
    mgmt_units: str = "log2cpm",
    mgmt_model: str = "exponential",
) -> Dict:
    """Set ``micro.cells.sensitivity`` for one drug under an EXPLICIT mode.

    In place, deterministic, no RNG.

    Parameters
    ----------
    micro : Microenvironment
        A reconstruction from :func:`reconstruct_real`. Modified in place.
    drug_name : str
        Case-insensitive drug name.
    sample : str, optional
        Xenium cache path (only used in ``'determinant'`` mode).
    mode : {'uniform', 'determinant'}
        See :data:`SENSITIVITY_MODES`.

        * ``'uniform'`` (**default**) — sensitivity is 1.0 everywhere, for every
          drug. This is the mode that makes a CROSS-DRUG ranking fair: the
          sensitivity factor is identical across drugs, so the comparison
          isolates delivery x potency. Every drug is then equally an upper bound.
        * ``'determinant'`` — each drug uses its own measured panel determinant
          where one exists and uniform 1.0 where none does
          (:func:`reconstruct.real_xenium.per_cell_sensitivity`). Any drug that
          falls through to uniform 1.0 is flagged ``upper_bound=True`` with a
          human-readable ``comparability`` string, because in this mode it is
          NOT comparable against drugs that were scaled down by a determinant.
    mgmt : float or str, optional
        **OPTIONAL, OPT-IN, BACKWARD-COMPATIBLE.** A MEASURED per-patient MGMT
        value for a drug whose resistance determinant is MGMT (temozolomide —
        see :data:`tmesim.molecular.mgmt.MGMT_CALIBRATED_DRUGS`). Either a
        continuous expression level (``float``, in ``log2(CPM + 1)`` by default;
        pass ``mgmt_units='cpm'`` for plain CPM) or a binary promoter status
        (``'methylated'`` / ``'unmethylated'``). Continuous expression is
        PREFERRED: in the source cohort MGMT *expression* was prognostic
        (CoxPH HR 1.85, p = 6.1e-05, n = 56) while promoter *methylation* was
        not (p = 0.16).

        **Default ``None`` leaves every code path exactly as it was** — no MGMT
        term is applied and the result is bit-identical to the pre-MGMT
        behaviour (asserted by ``tests/test_mgmt_scalar.py``).

        When supplied and applicable, the sensitivity array is multiplied by
        :func:`tmesim.molecular.mgmt.mgmt_sensitivity`, the drug is no longer an
        upper bound (``upper_bound=False``), and ``info`` records
        ``mgmt_applied``, ``mgmt_scalar``, ``mgmt_input``, ``mgmt_model``,
        ``mgmt_scalar_ci95`` and ``mgmt_calibration_n`` so the assumption can
        never be carried silently. Supplying ``mgmt`` for a drug that is not
        MGMT-calibrated is a NO-OP recorded as ``mgmt_applied=False`` with a
        reason (so a panel loop can pass it for every drug safely).

        MGMT is applied as ONE PER-PATIENT SCALAR over all cells: it is absent
        from both Xenium panels, has no panel surrogate (0/396 pass FDR), and is
        only pseudobulk-recoverable from Visium. No spatial MGMT field is
        asserted.
    mgmt_units : {'log2cpm', 'cpm'}
        Unit of a numeric ``mgmt``. Default is the calibration unit.
    mgmt_model : {'exponential', 'linear'}
        Which calibrated mapping to use; ``'exponential'`` is the default and the
        primary. See :mod:`tmesim.molecular.mgmt`.

    Notes
    -----
    This lives here rather than in the engine on purpose: ``core``, ``sim.engine``
    and ``physics`` are untouched — the engine keeps consuming a plain
    ``CellPopulation.sensitivity`` array and knows nothing about drugs or genes.

    ``'determinant'`` is only meaningful for the Xenium (single-cell,
    per-cell-expression) source; a Visium reconstruction has no per-cell
    transcripts, so sensitivity stays uniform 1.0 and that is reported (and
    flagged as an upper bound).

    Returns
    -------
    dict
        The ``info`` block (mode, determinant gene, kind, rationale, marker
        detection rate, resulting distribution, ``upper_bound``,
        ``comparability``), also stashed in
        ``micro.metadata['sensitivity_model_applied']`` so
        :func:`reconstruction_report` picks it up.
    """
    name = str(drug_name).strip().lower()
    m = str(mode).strip().lower()
    if m not in SENSITIVITY_MODES:
        raise ValueError(
            f"unknown sensitivity mode {mode!r}; expected one of {SENSITIVITY_MODES}"
        )

    if m == "uniform":
        micro.cells.sensitivity = np.ones(micro.cells.n, dtype=np.float64)
        info = {
            "drug": name,
            "mode": "uniform",
            "determinant_gene": None,
            "kind": "uniform 1.0 (identical for every drug)",
            "rationale": (
                "FAIR CROSS-DRUG MODE: sensitivity is uniform 1.0 for every drug, "
                "so it cancels out of any ranking and the comparison isolates "
                "DELIVERY x POTENCY with no sensitivity confound."
            ),
            "marker_detection_rate": None,
            "mean": 1.0, "std": 0.0, "min": 1.0, "max": 1.0,
            "is_uniform": True,
            "upper_bound": True,
            "comparability": (
                "COMPARABLE across drugs. Every drug is an equal "
                "no-intrinsic-resistance upper bound, so absolute kills are upper "
                "bounds but the RANKING is not confounded by sensitivity."
            ),
        }
        if mgmt is not None:
            info = _apply_mgmt_scalar(micro, info, mgmt, mgmt_units, mgmt_model)
        micro.metadata = dict(micro.metadata)
        micro.metadata["sensitivity_model_applied"] = info
        return info

    # --- mode == 'determinant' ------------------------------------------------
    source = micro.metadata.get("real_pipeline", {}).get("source")
    if source != "xenium":
        info = {
            "drug": name,
            "mode": "determinant",
            "determinant_gene": None,
            "kind": "none (uniform 1.0)",
            "rationale": f"source {source!r} carries no per-cell expression; "
                         "sensitivity left uniform 1.0",
            "marker_detection_rate": None,
            "mean": 1.0, "std": 0.0, "min": 1.0, "max": 1.0,
            "is_uniform": True,
        }
    else:
        data = _load_xenium_retry(sample)
        sens, info = per_cell_sensitivity(
            data["expression"], data["genes"], drug_name
        )
        if sens.size != micro.cells.n:
            raise ValueError(
                f"expression rows ({sens.size}) != cells ({micro.cells.n}); "
                "the reconstruction and the expression matrix are not aligned"
            )
        micro.cells.sensitivity = sens
        info = dict(info)

    info["mode"] = "determinant"
    info["upper_bound"] = bool(info.get("determinant_gene") is None)
    info["comparability"] = (
        _UPPER_BOUND_FLAG if info["upper_bound"] else
        "Scaled by its own measured determinant "
        f"({info['determinant_gene']}, mean {info['mean']:.3f}). NOT comparable "
        "against a drug whose determinant is unmeasured and therefore left at "
        "1.0 — use mode='uniform' for cross-drug ranking."
    )

    if mgmt is not None:
        info = _apply_mgmt_scalar(micro, info, mgmt, mgmt_units, mgmt_model)

    micro.metadata = dict(micro.metadata)
    micro.metadata["sensitivity_model_applied"] = info
    return info


def reconstruction_report(micro: Microenvironment) -> Dict:
    """Summarise a reconstructed real microenvironment: cell-type fractions, mean
    vascular / ECM / oxygen (over the tissue mask), the grid, and the per-cell
    SENSITIVITY MODEL in force (which every absolute kill fraction is conditional
    on -- see ``sensitivity`` in the returned dict)."""
    ct = np.asarray(micro.cells.cell_type)
    n = int(micro.cells.n)
    counts = {c.name: int((ct == int(c)).sum()) for c in CellType}
    fractions = {
        c.name: (counts[c.name] / n if n else 0.0) for c in CellType
    }

    tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)

    def _mean(field: Optional[np.ndarray]):
        if field is None:
            return None
        field = np.asarray(field, dtype=np.float64)
        return float(field[tissue].mean()) if tissue.any() else float(field.mean())

    # --- the per-cell sensitivity assumption (reported, never silent) --------
    sens = np.asarray(micro.cells.sensitivity, dtype=np.float64)
    tumor = ct == int(CellType.TUMOR)
    applied = micro.metadata.get("sensitivity_model_applied")
    sensitivity_block = {
        "mode": (applied.get("mode") if applied else None),
        "upper_bound": (applied.get("upper_bound") if applied else True),
        "comparability": (
            applied.get("comparability") if applied else
            "COMPARABLE across drugs — the reconstruction default is uniform 1.0 "
            "for every drug, so sensitivity cancels out of any ranking (and every "
            "absolute kill is an equal no-resistance upper bound)."
        ),
        "model": (
            f"mode={applied.get('mode')} ({applied['drug']} -> "
            f"{applied['determinant_gene'] or 'no determinant on panel'})"
            if applied else "uniform 1.0 (drug-agnostic reconstruction default)"
        ),
        "determinant_gene": applied["determinant_gene"] if applied else None,
        "kind": applied["kind"] if applied else "none (uniform 1.0)",
        "rationale": applied["rationale"] if applied else (
            "no per-cell intrinsic-resistance information is asserted; kill is pure "
            "delivery + potency physics"
        ),
        "marker_detection_rate": (
            applied["marker_detection_rate"] if applied else None
        ),
        # MGMT per-patient scalar (temozolomide). Present ONLY when the caller
        # supplied a measured MGMT value; None means the old uniform assumption
        # is in force and the TMZ kill is still an upper bound.
        "mgmt_applied": bool(applied.get("mgmt_applied", False)) if applied else False,
        "mgmt_scalar": (applied.get("mgmt_scalar") if applied else None),
        "mgmt_input": (applied.get("mgmt_input") if applied else None),
        "mgmt_model": (applied.get("mgmt_model") if applied else None),
        "mgmt_scalar_ci95": (applied.get("mgmt_scalar_ci95") if applied else None),
        "mgmt_calibration_n": (applied.get("mgmt_calibration_n") if applied else None),
        "mgmt_resolution": (applied.get("mgmt_resolution") if applied else None),
        "mgmt_note": (applied.get("mgmt_note") if applied else None),
        "is_uniform": bool(np.allclose(sens, sens[0])) if sens.size else True,
        "mean_all_cells": float(sens.mean()) if sens.size else None,
        "std_all_cells": float(sens.std()) if sens.size else None,
        "mean_tumor": float(sens[tumor].mean()) if tumor.any() else None,
        "std_tumor": float(sens[tumor].std()) if tumor.any() else None,
        "min": float(sens.min()) if sens.size else None,
        "max": float(sens.max()) if sens.size else None,
        "caveat": (
            "ABSOLUTE kill_fraction is CONDITIONAL on this sensitivity assumption. "
            "It is deterministic and documented, but it is NOT a measured drug "
            "response -- this section has no outcome label."
        ),
    }

    # --- which oxygen route produced fields.oxygen (never silent) ------------
    rp = micro.metadata.get("real_pipeline", {})
    oxygen_block = {
        "mode": rp.get("oxygen_mode", OXYGEN_MODE_DEFAULT),
        "source": rp.get("oxygen_source"),
        "is_surrogate": bool(rp.get("oxygen_is_surrogate", False)),
        "caveat": rp.get("oxygen_caveat"),
        "detail": rp.get("oxygen_detail"),
    }

    f = micro.fields
    return {
        "patient_id": micro.patient_id,
        "source": rp.get("source"),
        "n_cells": n,
        "n_tumor": int((ct == int(CellType.TUMOR)).sum()),
        "grid": tuple(micro.domain.shape),
        "dx_um": float(micro.domain.dx),
        "extent_um": tuple(float(v) for v in micro.domain.extent_um),
        "celltype_counts": counts,
        "celltype_fractions": fractions,
        "mean_vascular_density": _mean(f.vascular_density),
        "mean_ecm_density": _mean(f.ecm_density),
        "mean_oxygen": _mean(f.oxygen),
        "oxygen": oxygen_block,
        "oxygen_mode": oxygen_block["mode"],
        "ecm_markers_used": rp.get("ecm_markers_used", []),
        "tissue_mask_mode": micro.metadata.get("tissue_mask_mode"),
        "tissue_disc_radius_um": micro.metadata.get("tissue_disc_radius_um"),
        "tissue_coverage_fraction": float(tissue.mean()),
        "sensitivity": sensitivity_block,
    }


if __name__ == "__main__":
    micro = reconstruct_real("xenium")
    report = reconstruction_report(micro)
    print("=== REAL whole-niche microenvironment (xenium) ===")
    for k, v in report.items():
        print(f"  {k}: {v}")
