"""Synthetic glioblastoma (GBM) cohort generator.

This module fabricates biologically-plausible GBM microenvironments with a
*known* generative ground truth, so the whole tmesim pipeline
(reconstruct -> physics -> cells -> sim -> validation) can be exercised and
scored end-to-end before any real patient data exists.

Design goals (a FAIR test, not a giveaway)
------------------------------------------
* A **dense tumor core** with an **infiltrative margin** (glioma cells spread
  beyond the core along an exponential tail), embedded in reactive **stroma**,
  with **immune** cells (TAM/T-cell-like) enriched at the margin, and
  **endothelial** cells tracing a tortuous **vascular network**.
* Vasculature is **sparser inside the core** -> seeds hypoxia and poor drug
  delivery (the crux of GBM chemotherapy).
* **ECM density** varies spatially and is **higher in the core** (raises
  tortuosity, lowers effective diffusivity).
* Per-cell ``sensitivity`` in [0,1] is correlated with metadata
  (``mgmt_methylated``, ``subtype``) plus a spatially-localized **resistant
  sub-clone**.
* ``observed_response`` is a *noisy* function of BOTH transport-limitation
  (core cells receive less drug) AND intrinsic sensitivity, so it is **not**
  trivially recoverable from any single feature. See
  :func:`_compute_observed_response` for the exact generative formula.

Everything is deterministic given a seed (``np.random.default_rng``) and
vectorized.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.ndimage import gaussian_filter

from tmesim.core import (
    CellPopulation,
    CellState,
    CellType,
    Domain,
    Microenvironment,
    PatientRecord,
    TissueFields,
)

# ---------------------------------------------------------------------------
# Biological constants for the generative model
# ---------------------------------------------------------------------------
SUBTYPES: Tuple[str, ...] = ("proneural", "classical", "mesenchymal")

# Baseline intrinsic drug sensitivity by transcriptional subtype.
# Mesenchymal GBM is the most treatment-resistant; proneural the most sensitive.
_SUBTYPE_BASE_SENS: Dict[str, float] = {
    "proneural": 0.65,
    "classical": 0.50,
    "mesenchymal": 0.38,
}

# MGMT promoter methylation silences the repair enzyme -> alkylator (TMZ) works
# much better. Methylated => more sensitive.
_MGMT_METHYLATED_BONUS: float = 0.15
_MGMT_UNMETHYLATED_PENALTY: float = -0.05

# Response-model (transport x sensitivity) parameters.
_ECM_DELIVERY_ATTEN: float = 0.5   # how strongly ECM chokes local delivery
_DELIVERY_KD: float = 0.10         # half-saturation of the exposure->kill Hill
_DELIVERY_HILL: float = 1.4        # Hill coefficient of the exposure term
_RESPONSE_NOISE_SIGMA: float = 0.06  # additive per-patient measurement noise


# ---------------------------------------------------------------------------
# Small field helpers
# ---------------------------------------------------------------------------
def _radius_field(domain: Domain, cx_um: float, cy_um: float) -> np.ndarray:
    """Radial distance (um) of every grid cell from a center, shape (ny, nx)."""
    X, Y = domain.coords_um()
    return np.hypot(X - cx_um, Y - cy_um)


def _smooth_noise(rng: np.random.Generator, shape: Tuple[int, int],
                  sigma: float) -> np.ndarray:
    """A smooth random field in [0,1] via blurred white noise (min-max scaled)."""
    raw = rng.standard_normal(shape)
    blurred = gaussian_filter(raw, sigma=sigma, mode="reflect")
    lo, hi = blurred.min(), blurred.max()
    if hi - lo < 1e-12:
        return np.zeros(shape)
    return (blurred - lo) / (hi - lo)


def _build_vascular_field(rng: np.random.Generator, domain: Domain,
                          core_r_um: float, cx_um: float, cy_um: float,
                          base_level: float) -> np.ndarray:
    """A tortuous vascular-density field in [0,1], SPARSER in the tumor core.

    Vessels are rasterized random walks (tortuous, angiogenic-looking), blurred
    into a density, then suppressed inside the core radius to model the
    hypoxic, poorly-perfused GBM center.
    """
    ny, nx = domain.shape
    field = np.zeros((ny, nx))

    n_vessels = int(rng.integers(14, 24))
    steps = 70
    for _ in range(n_vessels):
        r = rng.integers(0, ny)
        c = rng.integers(0, nx)
        angle = rng.uniform(0, 2 * np.pi)
        for _ in range(steps):
            angle += rng.normal(0, 0.5)  # tortuosity
            r = int(np.clip(r + np.sin(angle) * 1.4, 0, ny - 1))
            c = int(np.clip(c + np.cos(angle) * 1.4, 0, nx - 1))
            field[r, c] += 1.0

    field = gaussian_filter(field, sigma=2.2, mode="reflect")
    # Add a diffuse capillary background so bulk tissue is perfused (a real
    # capillary bed, not just the large vessels).
    field = field + (0.4 + 0.5 * base_level) * _smooth_noise(rng, (ny, nx),
                                                             sigma=6.0)
    # Robust normalization by a high percentile (not the max), so a typical
    # well-perfused region reads near 1 and the field is realistically perfused.
    p98 = np.percentile(field, 98)
    if p98 > 0:
        field = field / p98
    field = np.clip(field, 0.0, 1.0)

    # Core suppression: vasculature drops to ~25% inside the core, ramping up
    # across the margin (smoothstep in radius).
    r_field = _radius_field(domain, cx_um, cy_um)
    ramp = np.clip((r_field - 0.6 * core_r_um) / (0.8 * core_r_um), 0.0, 1.0)
    core_factor = 0.25 + 0.75 * (ramp * ramp * (3 - 2 * ramp))  # smoothstep
    field = field * core_factor
    return np.clip(field, 0.0, 1.0)


def _build_ecm_field(rng: np.random.Generator, domain: Domain,
                     core_r_um: float, cx_um: float, cy_um: float) -> np.ndarray:
    """ECM density in [0,1], HIGHER in the dense core (raises tortuosity)."""
    r_field = _radius_field(domain, cx_um, cy_um)
    # High plateau in the core, decaying outward.
    core_ecm = 0.35 + 0.5 * np.exp(-(r_field / (1.1 * core_r_um)) ** 2)
    noise = 0.15 * _smooth_noise(rng, domain.shape, sigma=5.0)
    return np.clip(core_ecm + noise, 0.0, 1.0)


def _sample_cells_from_field(rng: np.random.Generator, domain: Domain,
                             density: np.ndarray, n: int
                             ) -> Tuple[np.ndarray, np.ndarray]:
    """Draw ``n`` cell positions (um) with probability proportional to a density
    field. Vectorized: weighted grid-cell choice + sub-cell jitter."""
    ny, nx = domain.shape
    weights = np.clip(density.astype(np.float64), 0.0, None).ravel()
    total = weights.sum()
    if total <= 0 or n <= 0:
        return np.empty(0), np.empty(0)
    probs = weights / total
    idx = rng.choice(ny * nx, size=n, replace=True, p=probs)
    row, col = np.divmod(idx, nx)
    x = (col + 0.5) * domain.dx + rng.uniform(-0.5, 0.5, size=n) * domain.dx
    y = (row + 0.5) * domain.dx + rng.uniform(-0.5, 0.5, size=n) * domain.dx
    x = np.clip(x, 0, domain.extent_um[0] - 1e-6)
    y = np.clip(y, 0, domain.extent_um[1] - 1e-6)
    return x, y


def _sample_field_at_cells(domain: Domain, field: np.ndarray,
                           x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Nearest-grid-cell lookup of a (ny,nx) field at cell positions."""
    col, row = domain.to_index(x, y)
    return field[row, col]


# ---------------------------------------------------------------------------
# Ground-truth response model
# ---------------------------------------------------------------------------
def _compute_observed_response(
    rng: np.random.Generator,
    sensitivity: np.ndarray,
    vascular_at_cell: np.ndarray,
    ecm_at_cell: np.ndarray,
    bbb: float,
) -> Tuple[float, float, float]:
    """Generate the noisy ground-truth drug response for the tumor cells.

    Generative formula (documented so the validation harness is honest)
    ------------------------------------------------------------------
    For each tumor cell *i*:

        delivery_i  = vascular_i * (1 - ALPHA_ecm * ecm_i)   # proxy local drug AUC
        efficacy    = 0.5 + bbb                              # BBB penetrance gate
        kill_i      = sensitivity_i * efficacy
                          * delivery_i**h / (delivery_i**h + KD**h)

    where ``vascular_i`` and ``ecm_i`` are the local field values at the cell
    (so core cells, with sparse vasculature + dense ECM, get little drug),
    ``bbb`` is the patient's blood-brain-barrier penetrance, and
    ``sensitivity_i`` is the intrinsic per-cell sensitivity (MGMT/subtype/clone).

    The population truth is the mean kill probability, and the *observed* label
    adds a small additive Gaussian measurement noise:

        true_response     = mean_i(kill_i)
        observed_response = clip(true_response + N(0, SIGMA), 0, 1)

    Because the response multiplies a spatial **transport** term by an intrinsic
    **sensitivity** term and then adds noise, it cannot be recovered from either
    the mean sensitivity alone (the intrinsic-sensitivity baseline) or the mean
    vascular fraction alone (a pure-transport reading). Beating those baselines
    requires actually resolving the spatial drug transport — which is the Stage 0
    hypothesis under test.

    Returns
    -------
    (observed_response, true_response, mean_delivery) all floats.
    """
    delivery = vascular_at_cell * (1.0 - _ECM_DELIVERY_ATTEN * ecm_at_cell)
    delivery = np.clip(delivery, 0.0, None)
    efficacy = 0.5 + bbb  # BBB penetrance gates overall drug efficacy
    hill = delivery ** _DELIVERY_HILL / (
        delivery ** _DELIVERY_HILL + _DELIVERY_KD ** _DELIVERY_HILL
    )
    kill = sensitivity * efficacy * hill
    true_response = float(kill.mean()) if kill.size else 0.0
    observed = float(np.clip(
        true_response + rng.normal(0.0, _RESPONSE_NOISE_SIGMA), 0.0, 1.0))
    return observed, true_response, float(delivery.mean())


# ---------------------------------------------------------------------------
# Per-cell sensitivity
# ---------------------------------------------------------------------------
def _assign_sensitivity(
    rng: np.random.Generator,
    n_tumor: int,
    subtype: str,
    mgmt_methylated: bool,
    in_resistant_clone: np.ndarray,
) -> np.ndarray:
    """Per-tumor-cell intrinsic sensitivity in [0,1].

    base(subtype) + MGMT effect, overridden by a low-sensitivity resistant
    sub-clone, plus per-cell biological noise.
    """
    base = _SUBTYPE_BASE_SENS[subtype]
    mgmt = _MGMT_METHYLATED_BONUS if mgmt_methylated else _MGMT_UNMETHYLATED_PENALTY
    sens = np.full(n_tumor, base + mgmt, dtype=np.float64)
    sens += rng.normal(0.0, 0.08, size=n_tumor)
    # Resistant sub-clone: collapse to near-resistant regardless of subtype.
    sens[in_resistant_clone] = 0.10 + rng.normal(0.0, 0.05,
                                                 size=int(in_resistant_clone.sum()))
    return np.clip(sens, 0.0, 1.0)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
def generate_synthetic_patient(patient_id: str, seed: int = 0,
                               grid: Tuple[int, int] = (120, 120),
                               dx: float = 15.0,
                               subtype: Optional[str] = None) -> PatientRecord:
    """Build a biologically-plausible synthetic GBM ``PatientRecord``.

    Creates a tumor mass + infiltrative margin, reactive stroma, immune and
    endothelial cells, a tortuous vascular network (sparser in the core ->
    hypoxia + poor delivery), spatially-varying ECM (denser in the core), and
    per-cell ``sensitivity`` correlated with metadata (``mgmt_methylated``,
    ``subtype``, a resistant sub-clone). ``observed_response`` is a noisy
    function of transport-limitation AND intrinsic sensitivity (see
    :func:`_compute_observed_response`). ``features`` holds non-spatial scalars
    for the baselines.

    Parameters
    ----------
    patient_id : identifier stored on the record and microenvironment.
    seed : reproducibility seed (``np.random.default_rng``).
    grid : (nx, ny) grid size.
    dx : grid spacing in micrometers.
    subtype : force a transcriptional subtype; if ``None`` one is drawn.

    Returns
    -------
    PatientRecord with a populated ``Microenvironment`` and ground-truth
    ``observed_response``.
    """
    rng = np.random.default_rng(seed)
    nx, ny = int(grid[0]), int(grid[1])
    domain = Domain(nx=nx, ny=ny, dx=dx)
    width_um, height_um = domain.extent_um

    # --- patient-level metadata -------------------------------------------
    if subtype is None:
        subtype = str(rng.choice(SUBTYPES))
    elif subtype not in _SUBTYPE_BASE_SENS:
        raise ValueError(f"unknown subtype {subtype!r}; choose from {SUBTYPES}")
    mgmt_methylated = bool(rng.random() < 0.45)
    bbb = float(np.clip(rng.normal(0.35, 0.08), 0.1, 0.7))  # BBB penetrance

    # --- tumor geometry ---------------------------------------------------
    cx = width_um * (0.5 + rng.uniform(-0.08, 0.08))
    cy = height_um * (0.5 + rng.uniform(-0.08, 0.08))
    core_r_um = min(width_um, height_um) * rng.uniform(0.22, 0.30)
    margin_um = core_r_um * rng.uniform(0.7, 1.1)  # infiltrative tail scale

    r_field = _radius_field(domain, cx, cy)
    tissue_noise = 0.6 + 0.4 * _smooth_noise(rng, domain.shape, sigma=4.0)

    # Tumor density: flat, dense core; exponential infiltrative margin beyond.
    tumor_density = np.where(
        r_field <= core_r_um, 1.0, np.exp(-(r_field - core_r_um) / margin_um))
    tumor_density = tumor_density * tissue_noise

    # --- continuum fields --------------------------------------------------
    vascular_base = float(np.clip(rng.normal(0.5, 0.15), 0.15, 0.9))
    vascular_density = _build_vascular_field(rng, domain, core_r_um, cx, cy,
                                             vascular_base)
    ecm_density = _build_ecm_field(rng, domain, core_r_um, cx, cy)

    # Cheap oxygen proxy (physics.oxygen will solve it for real) — used only to
    # report a hypoxic_fraction feature.
    oxygen_proxy = gaussian_filter(vascular_density, sigma=3.0, mode="reflect")
    if oxygen_proxy.max() > 0:
        oxygen_proxy = oxygen_proxy / oxygen_proxy.max()

    # --- cell populations --------------------------------------------------
    # Counts scale with grid area so density is grid-independent.
    area_factor = (nx * ny) / (120 * 120)
    n_tumor = int(rng.integers(2200, 2800) * area_factor)
    n_stroma = int(rng.integers(900, 1300) * area_factor)
    n_immune = int(rng.integers(500, 900) * area_factor)
    n_endo = int(rng.integers(400, 700) * area_factor)

    tx, ty = _sample_cells_from_field(rng, domain, tumor_density, n_tumor)

    # Reactive stroma: baseline everywhere + enhanced ring around the tumor.
    stroma_ring = np.exp(-((r_field - 1.15 * core_r_um) / (0.5 * core_r_um)) ** 2)
    stroma_density = (0.3 + 0.7 * stroma_ring) * tissue_noise
    sx, sy = _sample_cells_from_field(rng, domain, stroma_density, n_stroma)

    # Immune (TAM/T-cell): enriched at the infiltrative margin.
    margin_band = np.exp(-((r_field - core_r_um) / (0.6 * margin_um)) ** 2)
    immune_density = (0.25 + 0.75 * margin_band) * tissue_noise
    ix, iy = _sample_cells_from_field(rng, domain, immune_density, n_immune)

    # Endothelial: trace the vasculature (hence sparse in the core).
    ex, ey = _sample_cells_from_field(rng, domain, vascular_density, n_endo)

    # --- per-cell sensitivity (tumor) -------------------------------------
    # Resistant sub-clone: a blob offset inside the tumor.
    clone_cx = cx + rng.uniform(-0.6, 0.6) * core_r_um
    clone_cy = cy + rng.uniform(-0.6, 0.6) * core_r_um
    clone_r = core_r_um * rng.uniform(0.25, 0.45)
    clone_frac = rng.uniform(0.12, 0.28)  # target fraction (spatial gate below)
    dist_to_clone = np.hypot(tx - clone_cx, ty - clone_cy)
    in_clone = dist_to_clone < clone_r
    # Thin the clone so it is not a hard 100% blob (partial penetrance).
    keep = rng.random(n_tumor) < (clone_frac / max(in_clone.mean(), 1e-6))
    in_clone = in_clone & keep

    tumor_sens = _assign_sensitivity(rng, n_tumor, subtype, mgmt_methylated,
                                     in_clone)

    # Bystander cells: modest, irrelevant to tumor kill scoring.
    stroma_sens = np.clip(rng.normal(0.35, 0.1, n_stroma), 0, 1)
    immune_sens = np.clip(rng.normal(0.4, 0.1, n_immune), 0, 1)
    endo_sens = np.clip(rng.normal(0.3, 0.1, n_endo), 0, 1)

    # --- assemble CellPopulation ------------------------------------------
    x = np.concatenate([tx, sx, ix, ex])
    y = np.concatenate([ty, sy, iy, ey])
    cell_type = np.concatenate([
        np.full(n_tumor, CellType.TUMOR, dtype=np.int8),
        np.full(n_stroma, CellType.STROMAL, dtype=np.int8),
        np.full(n_immune, CellType.IMMUNE, dtype=np.int8),
        np.full(n_endo, CellType.ENDOTHELIAL, dtype=np.int8),
    ])
    sensitivity = np.concatenate([tumor_sens, stroma_sens, immune_sens, endo_sens])
    state = np.full(x.size, CellState.VIABLE, dtype=np.int8)

    cells = CellPopulation(x=x, y=y, cell_type=cell_type, state=state,
                           sensitivity=sensitivity)
    fields = TissueFields(ecm_density=ecm_density,
                          vascular_density=vascular_density)

    # --- ground-truth response --------------------------------------------
    vasc_at_tumor = _sample_field_at_cells(domain, vascular_density, tx, ty)
    ecm_at_tumor = _sample_field_at_cells(domain, ecm_density, tx, ty)
    oxy_at_tumor = _sample_field_at_cells(domain, oxygen_proxy, tx, ty)
    hypoxic_fraction = float((oxy_at_tumor < 0.2).mean()) if n_tumor else 0.0

    observed_response, true_response, mean_delivery = _compute_observed_response(
        rng, tumor_sens, vasc_at_tumor, ecm_at_tumor, bbb)

    # --- non-spatial features for baselines -------------------------------
    vascular_fraction = float(vascular_density.mean())
    mean_sensitivity = float(tumor_sens.mean())
    tumor_burden = float(n_tumor / (nx * ny))

    # Noisy "expression" scalars weakly tied to biology (partial signal for a
    # logistic-regression baseline, deliberately not fully informative).
    expr_mgmt = float(np.clip(0.7 - 0.5 * mgmt_methylated + rng.normal(0, 0.15),
                              0, 1.5))
    expr_vegf = float(np.clip(0.2 + 0.9 * hypoxic_fraction + rng.normal(0, 0.15),
                              0, 1.5))
    expr_egfr = float(np.clip(0.3 + 0.5 * (subtype == "classical")
                              + rng.normal(0, 0.15), 0, 1.5))
    expr_stem = float(np.clip(0.3 + 0.5 * (subtype == "proneural")
                              + rng.normal(0, 0.15), 0, 1.5))

    features: Dict[str, float] = {
        "mean_sensitivity": mean_sensitivity,
        "mgmt": float(mgmt_methylated),
        "subtype_proneural": float(subtype == "proneural"),
        "subtype_classical": float(subtype == "classical"),
        "subtype_mesenchymal": float(subtype == "mesenchymal"),
        "vascular_fraction": vascular_fraction,
        "hypoxic_fraction": hypoxic_fraction,
        "tumor_burden": tumor_burden,
        "bbb_penetrance": bbb,
        "expr_mgmt": expr_mgmt,
        "expr_vegf": expr_vegf,
        "expr_egfr": expr_egfr,
        "expr_stem": expr_stem,
    }

    metadata: Dict = {
        "subtype": subtype,
        "mgmt_methylated": mgmt_methylated,
        "bbb_penetrance": bbb,
        "tumor_center_um": (cx, cy),
        "core_radius_um": core_r_um,
        "resistant_clone_fraction": float(in_clone.mean()),
        "hypoxic_fraction": hypoxic_fraction,
        "true_response": true_response,   # noiseless label (for diagnostics)
        "mean_delivery": mean_delivery,
        "expression": {
            "MGMT": expr_mgmt, "VEGFA": expr_vegf,
            "EGFR": expr_egfr, "stemness": expr_stem,
        },
    }

    micro = Microenvironment(domain=domain, fields=fields, cells=cells,
                             patient_id=patient_id, metadata=metadata)

    return PatientRecord(patient_id=patient_id, micro=micro,
                         observed_response=observed_response, features=features)


def generate_synthetic_cohort(n: int = 40, seed: int = 0,
                              drug_names: Optional[List[str]] = None
                              ) -> List[PatientRecord]:
    """Generate a cohort spanning subtypes / vascularization / MGMT status.

    Subtypes are cycled so the cohort is balanced; per-patient seeds are drawn
    from a master ``default_rng(seed)`` so the whole cohort is deterministic.

    Parameters
    ----------
    n : number of patients.
    seed : master seed.
    drug_names : optional list of drug names (from ``DRUG_LIBRARY``) recorded in
        each patient's metadata; the first is used as the treating agent. The
        ground-truth response here is drug-agnostic (a TMZ-like alkylator), but
        the field is carried so downstream modules can key on it.

    Returns
    -------
    list[PatientRecord]
    """
    master = np.random.default_rng(seed)
    drug = (drug_names[0] if drug_names else "temozolomide-like")

    records: List[PatientRecord] = []
    for i in range(n):
        patient_seed = int(master.integers(0, 2**31 - 1))
        subtype = SUBTYPES[i % len(SUBTYPES)]  # balance subtypes across cohort
        rec = generate_synthetic_patient(
            patient_id=f"SYN-{i:03d}", seed=patient_seed, subtype=subtype)
        rec.micro.metadata["drug"] = drug
        records.append(rec)
    return records


def render_image(micro: Microenvironment) -> np.ndarray:
    """Render a multi-channel (H, W, C) image of the microenvironment.

    Produces a fluorescence/IHC-like stack the segmentation module can be
    exercised on. Cells are splatted as Gaussian blobs; each channel carries a
    marker whose intensity identifies a phenotype:

        channel 0 : nuclei (DAPI)      — all cells
        channel 1 : tumor marker       — TUMOR (e.g. SOX2)
        channel 2 : stroma marker      — STROMAL (e.g. GFAP)
        channel 3 : immune marker      — IMMUNE (e.g. IBA1/CD45)
        channel 4 : endothelial marker — ENDOTHELIAL (e.g. CD31)

    Image resolution matches the domain grid (pixel size == ``domain.dx`` um)
    so segmentation results land in the same coordinate frame. Returns a float32
    array in [0, 1].
    """
    domain = micro.domain
    cells = micro.cells
    ny, nx = domain.shape
    n_channels = 5
    img = np.zeros((ny, nx, n_channels), dtype=np.float64)

    col, row = domain.to_index(cells.x, cells.y)

    # Nuclei channel: every cell contributes.
    np.add.at(img[:, :, 0], (row, col), 1.0)

    marker_of_type = {
        int(CellType.TUMOR): 1,
        int(CellType.STROMAL): 2,
        int(CellType.IMMUNE): 3,
        int(CellType.ENDOTHELIAL): 4,
    }
    for ctype, chan in marker_of_type.items():
        m = cells.cell_type == ctype
        if m.any():
            np.add.at(img[:, :, chan], (row[m], col[m]), 1.0)

    # Blur into blobs (nuclei tighter than diffuse cytoplasmic markers).
    img[:, :, 0] = gaussian_filter(img[:, :, 0], sigma=0.9, mode="reflect")
    for chan in range(1, n_channels):
        img[:, :, chan] = gaussian_filter(img[:, :, chan], sigma=1.2,
                                          mode="reflect")

    # Normalize each channel to [0,1] and add mild imaging noise.
    rng = np.random.default_rng(abs(hash(micro.patient_id)) % (2**31))
    for chan in range(n_channels):
        c = img[:, :, chan]
        if c.max() > 0:
            c = c / c.max()
        c = np.clip(c + rng.normal(0, 0.01, size=c.shape), 0.0, 1.0)
        img[:, :, chan] = c

    return img.astype(np.float32)


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import numpy as _np

    print("=== single synthetic patient ===")
    rec = generate_synthetic_patient("SYN-DEMO", seed=1)
    print("summary:", rec.micro.summary())
    print("observed_response:", round(rec.observed_response, 4),
          "| true_response:", round(rec.micro.metadata["true_response"], 4))
    print("features:")
    for k, v in rec.features.items():
        print(f"    {k:22s} {v:.4f}")

    img = render_image(rec.micro)
    print("render_image shape:", img.shape, "dtype:", img.dtype,
          "range:", (float(img.min()), round(float(img.max()), 3)))
    assert img.ndim == 3 and img.shape[2] == 5, "image must be (H,W,C=5)"
    assert img.shape[:2] == rec.micro.domain.shape

    print("\n=== cohort of 12 ===")
    cohort = generate_synthetic_cohort(n=12, seed=0)
    responses = _np.array([r.observed_response for r in cohort])
    subtypes = [r.micro.metadata["subtype"] for r in cohort]
    mgmt = [r.micro.metadata["mgmt_methylated"] for r in cohort]
    print("n patients:", len(cohort))
    print("subtypes:", subtypes)
    print("mgmt methylated:", mgmt)
    print("responses:", _np.round(responses, 3).tolist())
    print(f"response  mean={responses.mean():.3f}  std={responses.std():.3f}  "
          f"min={responses.min():.3f}  max={responses.max():.3f}")

    # Sensitivity spread within a tumor (subtype + clone + noise).
    ts = cohort[0].micro.cells
    tsens = ts.sensitivity[ts.tumor_mask()]
    print(f"tumor sensitivity  mean={tsens.mean():.3f} "
          f"min={tsens.min():.3f} max={tsens.max():.3f}")

    # --- assertions: the test must be non-trivial -------------------------
    assert responses.std() > 0.03, "response must vary across the cohort"
    assert responses.min() < responses.max(), "responses collapsed"
    assert _np.all((responses >= 0) & (responses <= 1)), "response out of [0,1]"

    # Vasculature must be sparser in the core than the periphery.
    m = cohort[0].micro
    dom = m.domain
    cx, cy = m.metadata["tumor_center_um"]
    r = _radius_field(dom, cx, cy)
    core_mask = r < m.metadata["core_radius_um"]
    v = m.fields.vascular_density
    print(f"vascular core={v[core_mask].mean():.3f} "
          f"periphery={v[~core_mask].mean():.3f}")
    assert v[core_mask].mean() < v[~core_mask].mean(), \
        "vasculature must be sparser in the core"

    # ECM must be denser in the core.
    e = m.fields.ecm_density
    assert e[core_mask].mean() > e[~core_mask].mean(), \
        "ECM must be denser in the core"

    # Response should not be a pure function of mean sensitivity (needs transport).
    ms = _np.array([r.features["mean_sensitivity"] for r in cohort])
    corr = _np.corrcoef(ms, responses)[0, 1]
    print(f"corr(mean_sensitivity, response) = {corr:.3f} "
          "(should be <1: transport matters)")
    assert corr < 0.98, "response trivially recoverable from sensitivity alone"

    print("\nALL SELF-TESTS PASSED")
