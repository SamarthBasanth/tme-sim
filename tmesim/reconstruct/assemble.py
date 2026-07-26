"""Assemble validated pieces into a `Microenvironment` (`reconstruct` stage).

This is the final step of reconstruction: take a `Domain`, a `CellPopulation`,
and the continuum `TissueFields` (from `reconstruct.spatial` / `.vasculature`)
and bundle them into the central `Microenvironment` object the `sim` package
consumes. It performs the consistency checks that every downstream solver
implicitly assumes, so failures surface here with a clear message rather than as
a confusing shape error deep inside a PDE loop.
"""
from __future__ import annotations

from typing import Dict, Optional

import numpy as np

from ..core.domain import Domain
from ..core.microenvironment import (
    CellPopulation,
    Microenvironment,
    TissueFields,
)


def _check_field(name: str, arr: Optional[np.ndarray], shape: tuple[int, int]) -> None:
    """Validate an optional (ny, nx) field against the expected grid shape."""
    if arr is None:
        return
    arr = np.asarray(arr)
    if arr.shape != shape:
        raise ValueError(
            f"fields.{name} shape {arr.shape} != domain.shape {shape}"
        )
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"fields.{name} contains non-finite values (NaN/Inf)")


def assemble_microenvironment(
    domain: Domain,
    cells: CellPopulation,
    fields: TissueFields,
    patient_id: str = "reconstructed",
    metadata: Optional[Dict] = None,
) -> Microenvironment:
    """Bundle a domain, cells, and tissue fields into a `Microenvironment`.

    Validates shape/consistency before constructing the object:

    - every continuum field in ``fields`` (ECM, vascular, and any of
      oxygen/ifp/velocity that are set) matches ``domain.shape``;
    - the velocity field, if present, has shape ``(2, ny, nx)``;
    - all field values are finite;
    - every cell coordinate lies within the physical domain extent.

    Parameters
    ----------
    domain : Domain
        The grid the fields live on.
    cells : CellPopulation
        Discrete cells; their (x, y) must lie within ``domain.extent_um``.
    fields : TissueFields
        Continuum fields (from `reconstruct.spatial` / `.vasculature`).
    patient_id : str
        Identifier stored on the microenvironment.
    metadata : dict, optional
        Non-spatial patient features (subtype, mgmt_methylated, expression, ...)
        used by baselines and the sensitivity model. Defaults to ``{}``.

    Returns
    -------
    Microenvironment

    Raises
    ------
    ValueError
        On any shape mismatch, non-finite field, or out-of-bounds cell.
    """
    shape = domain.shape

    # --- required continuum fields ------------------------------------------
    _check_field("ecm_density", fields.ecm_density, shape)
    _check_field("vascular_density", fields.vascular_density, shape)

    # --- optional physics fields (may be None until physics runs) ------------
    _check_field("oxygen", fields.oxygen, shape)
    _check_field("ifp", fields.ifp, shape)
    if fields.velocity is not None:
        vel = np.asarray(fields.velocity)
        expected = (2, shape[0], shape[1])
        if vel.shape != expected:
            raise ValueError(
                f"fields.velocity shape {vel.shape} != expected {expected}"
            )
        if not np.all(np.isfinite(vel)):
            raise ValueError("fields.velocity contains non-finite values (NaN/Inf)")

    # --- cell coordinates within the physical extent -------------------------
    if cells.n > 0:
        width, height = domain.extent_um
        if not (np.all(np.isfinite(cells.x)) and np.all(np.isfinite(cells.y))):
            raise ValueError("cell coordinates contain non-finite values")
        if cells.x.min() < 0.0 or cells.x.max() > width:
            raise ValueError(
                f"cell x-coords [{cells.x.min():.1f}, {cells.x.max():.1f}] "
                f"outside domain extent [0, {width:.1f}] um"
            )
        if cells.y.min() < 0.0 or cells.y.max() > height:
            raise ValueError(
                f"cell y-coords [{cells.y.min():.1f}, {cells.y.max():.1f}] "
                f"outside domain extent [0, {height:.1f}] um"
            )

    return Microenvironment(
        domain=domain,
        fields=fields,
        cells=cells,
        patient_id=patient_id,
        metadata=metadata if metadata is not None else {},
    )


if __name__ == "__main__":
    # Self-test: Domain + random CellPopulation -> fields -> Microenvironment.
    from .spatial import build_tissue_fields

    rng = np.random.default_rng(0)
    dom = Domain(nx=120, ny=120, dx=15.0)
    n = 5000
    xmax, ymax = dom.extent_um
    pop = CellPopulation(
        x=rng.uniform(0, xmax, size=n),
        y=rng.uniform(0, ymax, size=n),
        cell_type=rng.integers(0, 5, size=n).astype(np.int8),
        state=np.zeros(n, dtype=np.int8),
        sensitivity=rng.uniform(0, 1, size=n),
    )
    fields = build_tissue_fields(dom, pop)
    micro = assemble_microenvironment(
        dom, pop, fields,
        patient_id="selftest-001",
        metadata={"subtype": "mesenchymal", "mgmt_methylated": True},
    )
    print("summary:")
    for k, v in micro.summary().items():
        print(f"  {k}: {v}")

    # Negative check: an out-of-bounds cell must raise.
    bad = pop.copy()
    bad.x[0] = xmax + 500.0
    try:
        assemble_microenvironment(dom, bad, fields)
    except ValueError as e:
        print("out-of-bounds correctly rejected:", str(e).split(" outside")[0], "...")

    print("assemble.py self-test OK")
