"""Regional blood-brain-barrier INTEGRITY — the additive paracellular leak term.

Why this module exists
----------------------
Before this module the model had **no barrier-integrity term anywhere**. The
delivery source was ``vascular_density * bbb_penetrance * plasma``, and
``bbb_penetrance`` is a pure function of the drug's 2D descriptors
(:func:`tmesim.molecular.permeability.bbb_penetrance`). Vessel *density* is not
barrier *leakiness*: two regions can carry identical vessel density and a
completely different permeability-surface-area product. So the model's predicted
contrast-enhancing-vs-non-enhancing fold change was **exactly 1.000 for every
drug**, against measured within-patient fold changes of **1.45x** (levetiracetam)
and **4.26x** (cefazolin), both at within-case permutation p < 0.0001
(``results/glioma_regional_pk.md`` section 4). That is the model failing to
express the answer, not mis-tuning it.

The functional FORM was derived from the data, not assumed
----------------------------------------------------------
Two candidate one-parameter region terms were fitted to the CE-vs-non-enhancing
contrast across **three independent drug/assay combinations**. A term that is
genuinely a property of the BARRIER rather than of the DRUG must come out
drug-INDEPENDENT, so the spread across drugs is the diagnostic:

===============  ==========================  ==============  ==============
drug             source                      multiplicative  additive leak
                                             f = CE/nonCE    L = CE - nonCE
===============  ==========================  ==============  ==============
levetiracetam    Oh 2026 (homogenate)        1.45            0.272
cefazolin        Oh 2026 (homogenate)        4.26            0.211
methotrexate     Blakeley 2008 (microdial.)  4.65            0.230
===============  ==========================  ==============  ==============

* multiplicative ``f`` spread across drugs: **3.2x**  -> drug-DEPENDENT, rejected
* additive ``L``  spread across drugs: **1.3x**  -> drug-independent, accepted

A multiplicative barrier factor **cannot** work: it would have to be ~4.7x for
the barrier-excluded drugs and only ~1.5x for the already-penetrant one, because
a drug that is already near-freely permeable has no headroom left to gain. The
additive form is what a **size-nonselective paracellular opening** should look
like — a parallel leak conductance that adds roughly a constant increment to the
partition regardless of the molecule's own physicochemistry.

Hence :data:`L_PARACELLULAR_LEAK` = 0.238 (the mean of 0.211 / 0.230 / 0.272) and

    effective_penetrance = clip(bbb_penetrance + L * bbb_disruption, 0, 1)

Sources
-------
* **Oh et al. 2026** *Neuro-Oncology* ``noaf204``, PMC12962645, Table S3 —
  125 surgical samples / 42 high-grade-glioma cases, levetiracetam (BBB-penetrant
  probe) and cefazolin (BBB-excluded probe) assayed in the SAME MRI-region-targeted
  resection specimens.
* **Blakeley et al. 2008** *J Neurooncol*, PMID 18787762 — methotrexate
  intracerebral microdialysis in contrast-enhancing vs non-enhancing tissue.

Both regional effects survive a mannitol-naive sensitivity analysis (mannitol
osmotically opens the barrier and is the obvious confounder): CE/non-enhancing
1.45x for levetiracetam and 4.38x for cefazolin in mannitol-naive samples.

Scope / honesty
---------------
* This is the **partition / delivery** axis only. Nothing here says anything
  about response, cell kill, or survival.
* ``L`` is fitted to *tissue:plasma ratio* units, which under the engine's
  Stage-0 quasi-equilibrium approximation is the same dimensionless quantity as
  ``bbb_penetrance``. It is NOT a permeability with physical units. The
  first-order rate the transport PDE needs is
  ``PS = PS_SCALE_S_INV * vascular_density * effective_penetrance`` (1/s), with
  :data:`tmesim.treatment.drugs.PS_SCALE_S_INV` carrying the whole ``1/s`` —
  that constant is where a reader should land for the units, and
  ``results/transvascular_exchange_revalidation.md`` section 6.2 for why.
* Two of the three fitted values come from tissue homogenate (total drug), one
  from microdialysis (free drug). They agree to 1.3x anyway, which is the reason
  to believe the additive form; but the mixing of modalities is a real caveat.
* The disruption FIELD below is a documented **proxy**, not a measurement — see
  :data:`REGION_DISRUPTION`.

Architecture: this module reads ``core`` and ``molecular`` and writes nothing
back. ``core``, ``sim.engine`` and ``physics`` are NOT modified.
"""
from __future__ import annotations

from typing import Dict, Mapping, Optional, Sequence, Union

import numpy as np

from ..core import DrugParams, Microenvironment
from ..core.domain import Domain

__all__ = [
    "L_PARACELLULAR_LEAK",
    "FITTED_LEAK_VALUES",
    "REGION_DISRUPTION",
    "effective_penetrance",
    "bbb_disruption_field",
    "disruption_from_region_labels",
    "barrier_delivery_source",
]


# ---------------------------------------------------------------------------
# The fitted constant
# ---------------------------------------------------------------------------

#: The three independently fitted additive leaks, in tissue:plasma ratio units.
#: ``L = mean(contrast-enhancing) - mean(non-enhancing)`` per drug/assay.
FITTED_LEAK_VALUES: Dict[str, float] = {
    # drug | study | assay
    "cefazolin | Oh 2026 Neuro-Oncol PMC12962645 | tissue homogenate": 0.2111,
    "methotrexate | Blakeley 2008 J Neurooncol PMID 18787762 | microdialysis": 0.2300,
    "levetiracetam | Oh 2026 Neuro-Oncol PMC12962645 | tissue homogenate": 0.2724,
}

#: Additive paracellular leak conductance, dimensionless (tissue:plasma ratio
#: units), for FULLY disrupted barrier (``bbb_disruption == 1``).
#:
#: Value = mean of :data:`FITTED_LEAK_VALUES` = mean(0.211, 0.230, 0.272)
#: = **0.238**, spread across the three drugs 1.3x (0.211 -> 0.272).
#:
#: This is a FITTED constant, not a first-principles one. It is fitted to three
#: measurements and it is the only free parameter this module introduces.
L_PARACELLULAR_LEAK: float = 0.238


# ---------------------------------------------------------------------------
# The disruption field — an HONEST PROXY, not a measurement
# ---------------------------------------------------------------------------

#: Per-Ivy-GAP-zone barrier disruption in [0, 1] (0 = intact, 1 = fully
#: disrupted), used to rasterise a ``bbb_disruption`` field from the region
#: labels produced by :mod:`tmesim.reconstruct.region_label`.
#:
#: **THIS IS A PROXY, AND IT IS NOT THE THING THAT WAS MEASURED.** The measured
#: quantity is an *MRI* label (T1-Gd contrast enhancement) on a resection
#: specimen. What this dictionary maps is a *transcriptional* zone label scored
#: against Ivy GAP laser-microdissected structure signatures (Puchalski et al.
#: 2018 *Science*). The bridge between them is anatomic, not measured in the same
#: samples: contrast enhancement on T1-Gd is produced by gadolinium leaking
#: through disrupted vessels, and the Ivy GAP structure whose defining feature is
#: exactly that vasculature is **microvascular proliferation (CTmvp)** — this
#: module's ``perivascular`` zone. The non-enhancing tissue the contrast is
#: measured against is the **leading edge / infiltrating tumour**, where vessels
#: are co-opted and largely intact — this module's ``margin`` zone. No study
#: assayed Ivy GAP zone labels and drug partition in the same tissue, so the
#: correspondence is an ANATOMIC ARGUMENT and should be reported as one.
#:
#: Values, and where each comes from:
#:
#: * ``perivascular`` = **1.00** — the CE anchor BY CONSTRUCTION. ``L`` was fitted
#:   as (CE - non-enhancing), so disruption must be 1 in the enhancing substrate
#:   and 0 in non-enhancing for the fit to reproduce itself.
#: * ``margin`` = **0.00** — the non-enhancing anchor, by the same construction.
#: * ``necrotic`` = **0.48** — MEASURED, not assumed. Necrosis-vs-non-enhancing
#:   additive leak is 0.0717 (levetiracetam) and 0.1553 (cefazolin), mean 0.1135;
#:   divided by L = 0.238 that is 0.477. Note this contradicts the naive
#:   expectation that necrosis is the LEAST perfused zone: measured, necrosis sits
#:   BETWEEN contrast-enhancing and non-enhancing, because necrotic cores sit
#:   inside the enhancing rim and share its disrupted barrier even though their
#:   perfusion is poor. Barrier integrity is not tissue viability.
#: * ``core`` = **0.50** — **UNMEASURED**. Cellular tumour spans the enhancing
#:   boundary, and no drug/region measurement isolates it. This is an explicit
#:   midpoint interpolation between the two anchors and is the weakest number
#:   here. Override it via the ``region_disruption`` argument if you have better
#:   information for a given tumour.
#:
#: An alternative proxy the validation also considered is a VEGFA/HIF1A
#: expression score off the IO panel. It was not chosen because VEGFA/HIF1A track
#: HYPOXIA, and the measurement above shows barrier disruption and hypoxia
#: dissociate — necrosis is the most hypoxic zone but only ~half as leaky as the
#: enhancing rim. The region-label route at least keeps the perivascular
#: (vessel-defined) zone as the disruption anchor, which is what enhancement
#: actually reports.
REGION_DISRUPTION: Dict[str, float] = {
    "perivascular": 1.00,
    "necrotic": 0.48,
    "core": 0.50,
    "margin": 0.00,
}


# ---------------------------------------------------------------------------
# The term itself
# ---------------------------------------------------------------------------
def effective_penetrance(
    bbb_penetrance: Union[float, np.ndarray],
    bbb_disruption: Union[float, np.ndarray],
    leak: float = L_PARACELLULAR_LEAK,
) -> np.ndarray:
    """Barrier-integrity-corrected penetrance, in ``[0, 1]``.

    .. math::

        P_\\mathrm{eff}(x, y) =
            \\mathrm{clip}\\bigl(P_\\mathrm{bbb} + L \\cdot D(x, y),\\; 0,\\; 1\\bigr)

    where :math:`P_\\mathrm{bbb}` is the drug's structure-only penetrance,
    :math:`D` is the per-location ``bbb_disruption`` in [0, 1], and :math:`L` is
    :data:`L_PARACELLULAR_LEAK`.

    The two routes are in PARALLEL: the drug's intrinsic transcellular crossing
    plus a drug-independent paracellular leak that opens where the barrier is
    broken. The clip at 1 enforces the free-drug ceiling — a partition above 1
    requires tissue binding or active uptake, neither of which this term models.

    ``bbb_disruption = 0`` returns ``bbb_penetrance`` unchanged (exactly:
    ``x + L*0.0 == x`` in IEEE-754, and the clip is a no-op because
    ``bbb_penetrance`` is already in [0, 1]).

    Parameters
    ----------
    bbb_penetrance : float or ndarray
        Structure-only penetrance in [0, 1] (broadcasts against ``bbb_disruption``).
    bbb_disruption : float or ndarray
        Barrier disruption in [0, 1]; 0 = intact, 1 = fully disrupted.
    leak : float
        The additive leak constant. Defaults to the fitted
        :data:`L_PARACELLULAR_LEAK`; exposed so a sensitivity sweep over the
        fitted range (0.211-0.272) can be run without monkeypatching.

    Returns
    -------
    np.ndarray
        Effective penetrance in [0, 1], broadcast shape.
    """
    p = np.asarray(bbb_penetrance, dtype=np.float64)
    d = np.asarray(bbb_disruption, dtype=np.float64)
    if np.any(d < -1e-12) or np.any(d > 1.0 + 1e-12):
        raise ValueError(
            f"bbb_disruption must lie in [0, 1]; got range "
            f"[{float(d.min()):.4g}, {float(d.max()):.4g}]"
        )
    return np.clip(p + float(leak) * d, 0.0, 1.0)


# ---------------------------------------------------------------------------
# Building the field
# ---------------------------------------------------------------------------
def _weighted_disruption(
    region_fields: Mapping[str, np.ndarray],
    region_disruption: Mapping[str, float],
    shape: Sequence[int],
) -> np.ndarray:
    """Fraction-weighted disruption from per-region fraction fields.

    ``region_field`` returns fields that sum to ~1 where any labeled location has
    support and to 0 elsewhere. Renormalising by that support total (rather than
    assuming it is 1) keeps unsupported grid points at disruption 0 instead of
    letting an arbitrary partial sum leak through.
    """
    num = np.zeros(tuple(shape), dtype=np.float64)
    den = np.zeros(tuple(shape), dtype=np.float64)
    unknown = [r for r in region_fields if r not in region_disruption]
    if unknown:
        raise KeyError(
            f"No disruption weight for region(s) {sorted(unknown)}. Known regions: "
            f"{sorted(region_disruption)}. Pass `region_disruption=` to extend."
        )
    for name, frac in region_fields.items():
        f = np.asarray(frac, dtype=np.float64)
        if f.shape != tuple(shape):
            raise ValueError(
                f"region field {name!r} has shape {f.shape}, expected {tuple(shape)}"
            )
        num += f * float(region_disruption[name])
        den += f
    return np.clip(np.where(den > 1e-12, num / np.maximum(den, 1e-12), 0.0), 0.0, 1.0)


def bbb_disruption_field(
    micro: Microenvironment,
    region_fields: Optional[Mapping[str, np.ndarray]] = None,
    region_disruption: Optional[Mapping[str, float]] = None,
) -> np.ndarray:
    """Per-location barrier disruption ``D(x, y)`` in ``[0, 1]`` on ``micro``'s grid.

    0 = intact barrier (non-enhancing / leading-edge-like tissue), 1 = fully
    disrupted (microvascular-proliferation / contrast-enhancing-like tissue).

    Sources of the region fields, in order of precedence:

    1. the explicit ``region_fields`` argument;
    2. ``micro.metadata["region_fields"]``, where the reconstruction pipeline
       may have stashed the output of
       :func:`tmesim.reconstruct.region_label.region_field`;
    3. neither -> an all-zero field (intact barrier everywhere), which makes
       every downstream call reduce EXACTLY to the pre-fix behaviour. This is
       deliberate: a microenvironment with no region information must not
       silently acquire a fabricated barrier defect.

    Parameters
    ----------
    micro : Microenvironment
        Supplies the target grid (``micro.domain.shape``) and, optionally,
        ``micro.metadata['region_fields']``.
    region_fields : mapping str -> (ny, nx) ndarray, optional
        Per-region spatial FRACTION fields in [0, 1], as returned by
        :func:`tmesim.reconstruct.region_label.region_field`.
    region_disruption : mapping str -> float, optional
        Override for :data:`REGION_DISRUPTION`. Merged over the defaults, so you
        may pass only the zones you want to change.

    Returns
    -------
    np.ndarray
        ``(ny, nx)`` disruption field in [0, 1].
    """
    shape = tuple(micro.domain.shape)
    if region_fields is None:
        region_fields = (micro.metadata or {}).get("region_fields")
    if not region_fields:
        return np.zeros(shape, dtype=np.float64)

    weights = dict(REGION_DISRUPTION)
    if region_disruption:
        weights.update({str(k): float(v) for k, v in region_disruption.items()})
    return _weighted_disruption(region_fields, weights, shape)


def disruption_from_region_labels(
    domain: Domain,
    coords_um: np.ndarray,
    region_labels: Sequence[str],
    bandwidth_um: float = 100.0,
    region_disruption: Optional[Mapping[str, float]] = None,
) -> np.ndarray:
    """Convenience: per-location Ivy GAP zone labels -> a ``bbb_disruption`` field.

    Rasterises ``region_labels`` with
    :func:`tmesim.reconstruct.region_label.region_field` and applies
    :data:`REGION_DISRUPTION`. Labels not in :data:`REGION_DISRUPTION` (e.g.
    ``"unlabeled"``) contribute no support and leave those grid points at 0.

    Parameters
    ----------
    domain : Domain
        Target grid.
    coords_um : ndarray (n_loc, 2)
        Location coordinates in micrometres, aligned to ``region_labels``.
    region_labels : sequence of str
        Zone per location, from :func:`tmesim.reconstruct.region_label.label_regions`.
    bandwidth_um : float
        Gaussian rasterisation bandwidth, passed through.
    region_disruption : mapping str -> float, optional
        Override for :data:`REGION_DISRUPTION`.

    Returns
    -------
    np.ndarray
        ``(ny, nx)`` disruption field in [0, 1].
    """
    from .region_label import region_field

    fields = region_field(domain, coords_um, region_labels, bandwidth_um=bandwidth_um)
    weights = dict(REGION_DISRUPTION)
    if region_disruption:
        weights.update({str(k): float(v) for k, v in region_disruption.items()})
    return _weighted_disruption(fields, weights, tuple(domain.shape))


# ---------------------------------------------------------------------------
# The delivery-source helper
# ---------------------------------------------------------------------------
def barrier_delivery_source(
    micro: Microenvironment,
    drug: DrugParams,
    plasma_conc: float,
    bbb_disruption: Optional[Union[float, np.ndarray]] = None,
    leak: float = L_PARACELLULAR_LEAK,
) -> np.ndarray:
    """Delivery source with the regional barrier term active.

    ::

        S(x, y) = vascular_density(x, y)
                  * clip(bbb_penetrance + L * bbb_disruption(x, y), 0, 1)
                  * plasma_conc

    which is :func:`tmesim.treatment.drugs.delivery_source` with the barrier term
    switched on. ``bbb_disruption=None`` resolves via
    :func:`bbb_disruption_field` (which yields all-zero when the microenvironment
    carries no region information), so the default reproduces the pre-fix source
    exactly.

    Note the leak multiplies ``vascular_density`` too: a paracellular opening can
    only leak where there is a vessel to leak from. Avascular tissue gets nothing
    regardless of how broken the barrier is.

    Parameters
    ----------
    micro : Microenvironment
        Supplies ``micro.fields.vascular_density`` and the grid.
    drug : DrugParams
        Supplies ``bbb_penetrance``.
    plasma_conc : float
        Blood concentration in uM at this instant.
    bbb_disruption : float or (ny, nx) ndarray, optional
        Barrier disruption in [0, 1]. ``None`` -> :func:`bbb_disruption_field`.
    leak : float
        Additive leak constant; defaults to :data:`L_PARACELLULAR_LEAK`.

    Returns
    -------
    np.ndarray
        ``(ny, nx)`` non-negative delivery-rate source, same units as
        :func:`tmesim.treatment.drugs.delivery_source`.
    """
    if bbb_disruption is None:
        bbb_disruption = bbb_disruption_field(micro)
    vascular = np.asarray(micro.fields.vascular_density, dtype=np.float64)
    p_eff = effective_penetrance(drug.bbb_penetrance, bbb_disruption, leak=leak)
    return vascular * p_eff * float(plasma_conc)


# ---------------------------------------------------------------------------
# self-test
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import numpy as np

    from ..core import (
        CellPopulation, CellState, CellType, Domain as _D, Microenvironment as _M,
        TissueFields,
    )

    # The fitted constant really is the mean of the three cited values.
    vals = list(FITTED_LEAK_VALUES.values())
    assert abs(float(np.mean(vals)) - L_PARACELLULAR_LEAK) < 5e-4, np.mean(vals)
    assert max(vals) / min(vals) < 1.4, "additive leak should be drug-independent to ~1.3x"

    # effective_penetrance: additive, clipped, and a no-op at disruption 0.
    assert effective_penetrance(0.25, 0.0) == 0.25
    assert abs(float(effective_penetrance(0.25, 1.0)) - 0.488) < 1e-12
    assert float(effective_penetrance(0.95, 1.0)) == 1.0  # clipped at the free ceiling

    domain = _D(nx=24, ny=20, dx=15.0)
    vasc = np.full(domain.shape, 0.4)
    fields = TissueFields(ecm_density=domain.ones() * 0.5, vascular_density=vasc)
    cells = CellPopulation(
        x=np.array([10.0]), y=np.array([10.0]),
        cell_type=np.array([CellType.TUMOR]), state=np.array([CellState.VIABLE]),
        sensitivity=np.array([0.5]),
    )
    micro = _M(domain=domain, fields=fields, cells=cells)
    drug = DrugParams(name="probe", bbb_penetrance=0.05)

    # No region info -> all-zero disruption -> pre-fix behaviour, bit-exact.
    from ..treatment.drugs import delivery_source
    assert np.array_equal(bbb_disruption_field(micro), np.zeros(domain.shape))
    assert np.array_equal(
        barrier_delivery_source(micro, drug, 10.0),
        delivery_source(micro, drug, 10.0),
    )

    # Region fields -> a graded field; perivascular is the maximum.
    half = domain.zeros(); half[:, :12] = 1.0
    other = domain.zeros(); other[:, 12:] = 1.0
    d_peri = bbb_disruption_field(micro, {"perivascular": half, "margin": other})
    assert np.allclose(d_peri[:, :12], 1.0) and np.allclose(d_peri[:, 12:], 0.0)

    src = barrier_delivery_source(micro, drug, 10.0, bbb_disruption=d_peri)
    ratio = src[:, :12].mean() / src[:, 12:].mean()
    # 0.05 -> 0.05 + 0.238 = 0.288, i.e. 5.76x for a barrier-EXCLUDED drug.
    assert abs(ratio - (0.05 + 0.238) / 0.05) < 1e-9, ratio

    # A penetrant drug gains far less in fold terms (the point of an ADDITIVE term).
    drug_pen = DrugParams(name="penetrant", bbb_penetrance=0.84)
    src_p = barrier_delivery_source(micro, drug_pen, 10.0, bbb_disruption=d_peri)
    ratio_p = src_p[:, :12].mean() / src_p[:, 12:].mean()
    assert 1.0 < ratio_p < 1.5, ratio_p
    assert ratio_p < ratio

    print("barrier.py self-test PASSED")
    print(f"  L_PARACELLULAR_LEAK = {L_PARACELLULAR_LEAK} "
          f"(mean of {sorted(round(v, 4) for v in vals)}, spread "
          f"{max(vals) / min(vals):.2f}x)")
    print(f"  CE/non-enhancing fold, excluded drug  (bbb=0.05): {ratio:.2f}x "
          f"(measured cefazolin 4.26x)")
    print(f"  CE/non-enhancing fold, penetrant drug (bbb=0.84): {ratio_p:.2f}x "
          f"(measured levetiracetam 1.45x)")
