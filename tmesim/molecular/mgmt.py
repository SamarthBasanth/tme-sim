"""MGMT as a **per-patient scalar** for temozolomide sensitivity.

Why this module exists
----------------------
``results/gbm_tmz_exvivo.md`` measured what our spatial model forfeits by leaving
temozolomide (TMZ) sensitivity at a uniform ``1.0``:

* MGMT expression predicts MEASURED ex-vivo TMZ response at Spearman
  **rho = +0.516 (p = 0.034, n = 17)** in patient-derived GBM stem-cell cultures
  (GEO **GSE232173** / Ntafoulis *et al.*, Br J Cancer 2023, **PMID 37620410**),
  independently reproducing the paper's own rho = +0.71 from raw counts.
* That is **~27 % of the between-model rank variance** (rho^2), an **81
  percentage-point** viability spread, and a leave-one-out MAE of 17.0 pp
  (constant model) vs 13.2 pp (one MGMT covariate).
* MGMT is on **neither** Xenium panel (739-gene union), **no** panel gene
  surrogates for it (0/396 pass both FDR bars; panel MGMT-correlation base rate
  14.1 % vs a 12.3 % genome-wide base rate, i.e. **no enrichment**), and GDSC is
  structurally blind to it (rho = -0.068, p = 0.71 -- GDSC dosed TMZ over
  0.03-30 uM against a 25-300 uM effective range).

The actionable insight: **MGMT does not need to be spatial.** It is a
per-patient scalar. Visium recovers it as pseudobulk (5-29 % of spots positive
across our four real sections) and clinical MGMT testing is routine. So the model
can accept a measured per-patient MGMT value and scale TMZ sensitivity by it,
turning an unmodelled determinant into a **measured input**.

What is calibrated here
-----------------------
:func:`mgmt_sensitivity` maps a per-patient MGMT level onto the multiplicative
``sensitivity`` factor the engine already consumes
(``cells/pkpd.py``: ``p_kill = p_hill * sensitivity * R_hypoxia * R_quiescence``).

**Primary (continuous) mapping -- exponential in kill space.** Because
``sensitivity`` enters the PD model *multiplicatively*, the natural functional
form is multiplicative too. Let ``x = log2(MGMT_CPM + 1)`` and let
``kill = 1 - viability/100`` be the measured ex-vivo TMZ kill at 100 uM. Fitting
``ln(kill) = c - lambda * x`` by OLS over the n = 17 paired cultures gives

    lambda = 0.3429 per log2 unit   (95 % CI 0.0059-0.6798, p = 0.0466, n = 17)
    exp(c) = 0.2438                 (the fitted kill of an MGMT-silent culture)

and therefore

    s(x) = exp(-lambda * x)          in (0, 1]

anchored so that **an MGMT-silent tumour (CPM = 0 -> x = 0) has s = 1.0**. That
anchor is deliberate: it makes the current uniform-1.0 default *exactly* the
MGMT-null special case of this mapping, so supplying no MGMT value changes
nothing and supplying a measured one can only scale kill **down**.

DIRECTION (the biology, Esteller PMID 11070098 / Hegi PMID 15758010):
**high MGMT -> O6-meG lesions repaired -> RESISTANT -> LOW sensitivity.**

**Secondary mapping -- linear in viability.** The report's own OLS,
``viability% = 65.72 + 6.754 * x`` (slope 95 % CI 1.47-12.04, p = 0.0156,
r^2 = 0.331), inverted the same way: ``s = (1 - v_hat/100) / 0.34285``. It is the
better-fitting model in viability space (LOO MAE 13.2 vs 14.5 pp) but it crosses
**zero kill at x = 5.08** and goes negative above it, i.e. it is not usable as a
bounded multiplier at the top of the observed Visium range. It is offered as
``model='linear'`` for sensitivity analysis, clipped at the floor, and it is NOT
the default.

**Binary mapping.** For a clinically-reported methylated / unmethylated status,
the measured group contrast is used directly rather than the regression: median
kill 0.2944 (MGMT-silent, n = 9) vs 0.1155 (MGMT-expressing, n = 8),
Mann-Whitney **p = 0.0152**, Cohen d = 1.30 -> a ratio of **0.392**
(bootstrap 95 % CI 0.097-0.629). Reassuringly this is an *independent* check on
the continuous fit: averaging ``s(x)`` within the same two groups gives 0.942
(silent) and 0.309 (expressing), a ratio of 0.328 -- the two routes agree.

PREFER CONTINUOUS EXPRESSION OVER BINARY METHYLATION. In this same cohort
(MOESM3 ``CoxPH_tissue_005``, n = 56 parental tumours) MGMT **expression** was
prognostic for overall survival (**HR 1.85, p = 6.1e-05**, q = 0.0079, externally
replicated in GLIOTRAIN n = 70 at HR 1.21, p = 0.025) while MGMT **promoter
methylation** -- the routine clinical biomarker -- was **not** (p = 0.16, same
cohort). Methylation is a noisy, thresholded proxy for the thing that actually
matters, which is how much MGMT protein the tumour can make. Use
``mgmt_sensitivity(level)`` where an expression number exists and fall back to
``mgmt_sensitivity(status='unmethylated')`` only when it does not.

Scope
-----
The calibration is **temozolomide-specific** (see :data:`MGMT_CALIBRATED_DRUGS`).
The nitrosoureas (lomustine, carmustine) and procarbazine are also O6-alkylating
agents whose lesions MGMT reverses, but no measured dose-response was fitted for
them here, so this module refuses to extrapolate to them.

Honest limits (repeated in ``results/mgmt_scalar.md``)
------------------------------------------------------
* **n = 17.** The 95 % CI on ``lambda`` spans 0.006-0.680, i.e. a factor of >100
  in slope. The *direction* is solid; the *magnitude* is approximate.
* The response variable is the operator-**reproducibility** table (%CV up to
  43 %, median 8.3 %), which attenuates the fit.
* The cultures are **dissociated 2-D GSC monolayers** -- no architecture, no
  vasculature, no diffusion gradient. Nothing spatial is calibrated here.
* Visium MGMT is **pseudobulk**, not a resolvable spatial field
  (:func:`mgmt_from_visium`), and it is a different assay chemistry from the bulk
  RNA-seq the mapping was fitted on.
* A **clinically measured** MGMT (the routine assay, ideally quantitative
  expression) is the proper input to this function; the Visium numbers are a
  demonstration that a real per-patient value is recoverable from our own data.

Nothing here modifies ``core`` / ``physics`` / ``sim`` / ``cells``. The constants
below are re-derivable from the cached primary bytes by
:func:`tmesim.validation.mgmt_scalar.calibrate_from_exvivo`, and
``tests/test_mgmt_scalar.py`` asserts they still agree.
"""
from __future__ import annotations

import math
from typing import Dict, Optional, Sequence, Union

# --------------------------------------------------------------------------- #
# Provenance of every constant below
# --------------------------------------------------------------------------- #
#: Primary source for the calibration.
CALIBRATION_SOURCE = {
    "citation": ("Ntafoulis I, Kleijn A, Ju J, et al. Br J Cancer "
                 "2023;129(8):1327-1338"),
    "pmid": "37620410",
    "geo": "GSE232173",
    "response_metric": "mean % cell viability at 100 uM TMZ vs DMSO, 6-day "
                       "ATP (CellTiter-Glo) assay -- LOWER = MORE SENSITIVE",
    "viability_table": "BJC2023_MOESM2.docx, Supplemental Table 2 "
                       "(operator-reproducibility table, n = 40 cultures)",
    "expression": "GSE232173_count_data_24_cells.txt.gz, raw counts, "
                  "34,461 genes x 24 GSC cultures, simple CPM",
    "n_paired": 17,
    "report": "results/gbm_tmz_exvivo.md",
}

#: Expression unit the mapping is calibrated in.
MGMT_UNITS = "log2(CPM + 1)"

#: n of the calibration join (cultures with BOTH a published viability number
#: and paired RNA-seq).
MGMT_CALIBRATION_N = 17

#: Observed MGMT range in the calibration set, in CPM and in log2(CPM+1).
MGMT_OBSERVED_CPM_RANGE = (0.0, 21.119620297386547)
MGMT_OBSERVED_LOG2CPM_RANGE = (0.0, 4.467254715551343)

# --- PRIMARY MAPPING: exponential decay of kill with MGMT ------------------- #
#: Decay constant of measured TMZ kill per log2 unit of MGMT CPM.
#: OLS of ln(kill) on log2(CPM+1), n = 17.
MGMT_LAMBDA = 0.34285104563373975
#: 95 % CI on ``MGMT_LAMBDA`` (t-based, n-2 = 15 df). WIDE -- see honest limits.
MGMT_LAMBDA_CI95 = (0.0058882519035368674, 0.6798138393639426)
#: Bootstrap (10,000 resamples, seed 20260727) 95 % CI on ``MGMT_LAMBDA``.
MGMT_LAMBDA_BOOT_CI95 = (0.012536403419776635, 0.6448370887808658)
#: Two-sided p for lambda != 0.
MGMT_LAMBDA_P = 0.04659205084528332
#: r^2 of the log-kill fit.
MGMT_LAMBDA_R2 = 0.23870390039059564
#: Fitted kill of an MGMT-silent culture, ``exp(c)``. Reported for provenance;
#: the sensitivity mapping is a RATIO so this cancels out.
MGMT_KILL_AT_ZERO = 0.24375679877143605
#: Kill floor used when taking ln(kill): 2/17 cultures measured >= 99 % viability
#: (one at 104.2 %, i.e. more cells than the DMSO control), which is assay noise,
#: not negative kill. Those two are pinned at 1 % kill rather than dropped.
MGMT_KILL_FLOOR = 0.01

# --- SECONDARY MAPPING: linear in viability (the report's own OLS) ---------- #
#: ``viability% = MGMT_VIAB_INTERCEPT + MGMT_VIAB_SLOPE * log2(CPM+1)``.
MGMT_VIAB_SLOPE = 6.754163030238055
MGMT_VIAB_INTERCEPT = 65.7154902291493
MGMT_VIAB_SLOPE_CI95 = (1.4714628688087128, 12.036863191667397)
MGMT_VIAB_SLOPE_P = 0.015648635171249592
MGMT_VIAB_R2 = 0.331146930922554
#: log2(CPM+1) at which the linear model predicts 100 % viability, i.e. ZERO
#: kill. Above this the linear mapping is meaningless -- which is exactly why it
#: is not the default (three of our four Visium sections sit at x >= 4.0).
MGMT_VIAB_ZERO_KILL_X = (100.0 - MGMT_VIAB_INTERCEPT) / MGMT_VIAB_SLOPE

# --- BINARY MAPPING: measured group contrast -------------------------------- #
#: CPM at/above which a culture was called MGMT-EXPRESSING in the calibration.
MGMT_EXPRESSED_CPM_THRESHOLD = 1.0
#: Median measured kill in each group (n = 9 silent, n = 8 expressing).
MGMT_BINARY_KILL = {"silent": 0.2944, "expressing": 0.11550000000000005}
#: Sensitivity multiplier for a binary status. 'methylated' == promoter
#: methylated == MGMT silenced == the SENSITIVE group -> 1.0 (the same anchor the
#: continuous mapping uses). 'unmethylated' == MGMT expressed == RESISTANT.
MGMT_BINARY_SENSITIVITY = {
    "methylated": 1.0,
    "unmethylated": 0.39232336956521796,
}
#: Bootstrap 95 % CI (10,000 resamples, seed 20260727) on the unmethylated value.
MGMT_BINARY_CI95 = (0.0968114602587799, 0.629230769230769)
#: Mann-Whitney p for the two-group viability contrast, and Cohen's d.
MGMT_BINARY_P = 0.015220074043603454
MGMT_BINARY_COHENS_D = 1.30

#: Cross-check: averaging the CONTINUOUS s(x) within the same two groups.
#: Independent of the binary numbers above; they agree (0.309/0.942 = 0.328 vs
#: the measured 0.392).
MGMT_BINARY_FROM_CONTINUOUS = {"methylated": 0.942210523045451,
                               "unmethylated": 0.3090353780799962}

#: Floor on the returned multiplier. A measured-high MGMT tumour is resistant,
#: not immortal: MGMT reverses the O6-meG lesion but TMZ also produces N7-meG and
#: N3-meA adducts that MGMT does not touch, and the assay itself never observed a
#: zero-kill culture. A hard 0 would let one expression number assert perfect
#: resistance, so the mapping is floored. Declared constant, not fitted.
MGMT_SENSITIVITY_FLOOR = 0.05

#: Drugs this calibration may be applied to. TMZ only: the fit is a TMZ
#: dose-response. The nitrosoureas and procarbazine are also O6-alkylators whose
#: lesions MGMT reverses, but no measured curve was fitted for them here and this
#: module will not extrapolate.
MGMT_CALIBRATED_DRUGS = ("temozolomide",)

#: Aliases accepted by the ``status`` argument, mapped onto the two groups.
_STATUS_ALIASES = {
    "methylated": "methylated", "m": "methylated",
    "mgmt-methylated": "methylated", "mgmt_methylated": "methylated",
    "silenced": "methylated", "silent": "methylated", "negative": "methylated",
    "unmethylated": "unmethylated", "u": "unmethylated", "um": "unmethylated",
    "mgmt-unmethylated": "unmethylated", "mgmt_unmethylated": "unmethylated",
    "expressed": "unmethylated", "expressing": "unmethylated",
    "positive": "unmethylated",
}

MGMT_MODELS = ("exponential", "linear")


# --------------------------------------------------------------------------- #
# The mapping
# --------------------------------------------------------------------------- #
def _to_log2cpm(level: float, units: str) -> float:
    u = str(units).strip().lower()
    if u in ("log2cpm", "log2(cpm+1)", "log2_cpm", "log2cpm1"):
        x = float(level)
    elif u in ("cpm", "tpm", "linear_cpm"):
        if float(level) < 0.0:
            raise ValueError(f"MGMT CPM must be >= 0, got {level!r}")
        x = math.log2(float(level) + 1.0)
    else:
        raise ValueError(
            f"unknown units {units!r}; expected 'log2cpm' (the calibration unit, "
            f"{MGMT_UNITS}) or 'cpm'"
        )
    if x < 0.0:
        raise ValueError(
            f"MGMT level maps to log2(CPM+1) = {x:.4g} < 0, which is impossible; "
            "check the units"
        )
    return x


def mgmt_sensitivity(
    mgmt_level: Optional[float] = None,
    *,
    status: Optional[str] = None,
    units: str = "log2cpm",
    model: str = "exponential",
    floor: float = MGMT_SENSITIVITY_FLOOR,
    detail: bool = False,
) -> Union[float, Dict]:
    """Per-patient MGMT level -> temozolomide sensitivity multiplier in ``[0, 1]``.

    Exactly one of ``mgmt_level`` (preferred) or ``status`` must be given.

    Parameters
    ----------
    mgmt_level : float, optional
        Continuous MGMT expression. By default in the calibration unit
        ``log2(CPM + 1)``; pass ``units='cpm'`` to supply plain CPM/TPM. This is
        the PREFERRED input -- MGMT *expression* was prognostic in the source
        cohort (CoxPH HR 1.85, p = 6.1e-05, n = 56) while promoter *methylation*
        was not (p = 0.16).
    status : {'methylated', 'unmethylated'}, optional
        A clinically reported binary MGMT promoter status (aliases in
        :data:`_STATUS_ALIASES`). ``'methylated'`` = silenced = sensitive = 1.0;
        ``'unmethylated'`` = expressed = resistant = 0.392 (measured median-kill
        ratio, bootstrap 95 % CI 0.097-0.629). Use only when no expression
        number exists.
    units : {'log2cpm', 'cpm'}
        Unit of ``mgmt_level``. Default is the calibration unit.
    model : {'exponential', 'linear'}
        ``'exponential'`` (default, PRIMARY): ``s = exp(-lambda * x)``, bounded,
        monotone, multiplicative -- matching how ``sensitivity`` enters the PD
        model. ``'linear'`` (SECONDARY, sensitivity analysis only): invert the
        report's viability OLS. The linear model reaches zero kill at
        ``x = 5.08`` and is negative above it, so it is clipped at ``floor``.
    floor : float
        Lower clip. Default :data:`MGMT_SENSITIVITY_FLOOR` = 0.05 -- a
        measured-high-MGMT tumour is resistant, not immortal (TMZ's N7-meG /
        N3-meA adducts are not MGMT substrates, and no zero-kill culture was
        observed in n = 40).
    detail : bool
        If True return a dict with the value, the model, the inputs, the fit
        constants, and the extrapolation flag, instead of a bare float.

    Returns
    -------
    float in ``[floor, 1.0]`` (or dict if ``detail``).

    Notes
    -----
    ``mgmt_sensitivity(0.0) == 1.0`` exactly, by construction: an MGMT-silent
    tumour is the anchor, so the model's existing uniform-1.0 TMZ default is the
    MGMT-null special case of this mapping and nothing changes when no MGMT value
    is supplied.

    Examples
    --------
    >>> round(mgmt_sensitivity(0.0), 6)          # MGMT-silent -> the anchor
    1.0
    >>> round(mgmt_sensitivity(4.0290), 4)       # Visium GBM1 pseudobulk
    0.2512
    >>> round(mgmt_sensitivity(status='unmethylated'), 4)
    0.3923
    """
    if (mgmt_level is None) == (status is None):
        raise ValueError(
            "supply exactly one of mgmt_level (preferred -- expression is the "
            "prognostic variable) or status (binary methylation)"
        )
    if not (0.0 <= float(floor) <= 1.0):
        raise ValueError(f"floor must be in [0, 1], got {floor!r}")

    if status is not None:
        key = _STATUS_ALIASES.get(str(status).strip().lower())
        if key is None:
            raise ValueError(
                f"unknown MGMT status {status!r}; expected 'methylated' "
                f"(silenced/sensitive) or 'unmethylated' (expressed/resistant)"
            )
        value = float(MGMT_BINARY_SENSITIVITY[key])
        value = min(1.0, max(float(floor), value))
        if not detail:
            return value
        return {
            "sensitivity": value,
            "input": {"status": key},
            "model": "binary (measured median-kill ratio)",
            "n": MGMT_CALIBRATION_N,
            "ci95": (None if key == "methylated" else MGMT_BINARY_CI95),
            "p": MGMT_BINARY_P,
            "cohens_d": MGMT_BINARY_COHENS_D,
            "extrapolated": False,
            "direction": "high MGMT -> resistant -> LOW sensitivity",
            "caveat": (
                "BINARY input. Promoter methylation was NOT prognostic in the "
                "source cohort (p = 0.16, n = 56) while MGMT EXPRESSION was "
                "(HR 1.85, p = 6.1e-05). Prefer mgmt_level."
            ),
            "source": CALIBRATION_SOURCE,
        }

    m = str(model).strip().lower()
    if m not in MGMT_MODELS:
        raise ValueError(f"unknown model {model!r}; expected one of {MGMT_MODELS}")
    x = _to_log2cpm(float(mgmt_level), units)

    if m == "exponential":
        raw = math.exp(-MGMT_LAMBDA * x)
        lo = math.exp(-MGMT_LAMBDA_CI95[1] * x)
        hi = math.exp(-MGMT_LAMBDA_CI95[0] * x)
    else:  # linear-in-viability, inverted
        v_hat = MGMT_VIAB_INTERCEPT + MGMT_VIAB_SLOPE * x
        raw = (1.0 - v_hat / 100.0) / (1.0 - MGMT_VIAB_INTERCEPT / 100.0)
        lo = (1.0 - (MGMT_VIAB_INTERCEPT + MGMT_VIAB_SLOPE_CI95[1] * x) / 100.0) / (
            1.0 - MGMT_VIAB_INTERCEPT / 100.0)
        hi = (1.0 - (MGMT_VIAB_INTERCEPT + MGMT_VIAB_SLOPE_CI95[0] * x) / 100.0) / (
            1.0 - MGMT_VIAB_INTERCEPT / 100.0)

    value = min(1.0, max(float(floor), float(raw)))
    if not detail:
        return value

    extrapolated = x > MGMT_OBSERVED_LOG2CPM_RANGE[1]
    return {
        "sensitivity": value,
        "sensitivity_unclipped": float(raw),
        "input": {"mgmt_level": float(mgmt_level), "units": units,
                  "log2cpm": x},
        "model": m,
        "n": MGMT_CALIBRATION_N,
        "lambda": MGMT_LAMBDA if m == "exponential" else None,
        "ci95": (min(1.0, max(float(floor), lo)),
                 min(1.0, max(float(floor), hi))),
        "ci95_basis": ("t-based 95 % CI on the fitted slope, n-2 = 15 df -- it is "
                       "WIDE (lambda 0.006-0.680); direction is solid, magnitude "
                       "is approximate"),
        "p": MGMT_LAMBDA_P if m == "exponential" else MGMT_VIAB_SLOPE_P,
        "r2": MGMT_LAMBDA_R2 if m == "exponential" else MGMT_VIAB_R2,
        "clipped": value != float(raw),
        "floor": float(floor),
        "extrapolated": bool(extrapolated),
        "extrapolation_note": (
            f"log2(CPM+1) = {x:.3f} is ABOVE the calibration ceiling "
            f"{MGMT_OBSERVED_LOG2CPM_RANGE[1]:.3f}; the value is an extrapolation"
            if extrapolated else None
        ),
        "direction": "high MGMT -> resistant -> LOW sensitivity",
        "anchor": "s(MGMT = 0 CPM) = 1.0, so the uniform-1.0 default is the "
                  "MGMT-null special case",
        "source": CALIBRATION_SOURCE,
    }


def is_mgmt_calibrated(drug_name: str) -> bool:
    """Is the MGMT scalar calibrated for this drug? (Temozolomide only.)"""
    return str(drug_name).strip().lower() in MGMT_CALIBRATED_DRUGS


# --------------------------------------------------------------------------- #
# Recovering a per-patient MGMT level from our OWN real data
# --------------------------------------------------------------------------- #
def mgmt_from_visium(sample_id: str) -> Dict:
    """Per-sample **pseudobulk** MGMT level from a REAL cached Visium GBM section.

    Reads the section with :func:`tmesim.data.real.visium.load_sample` and sums
    the MGMT row over all in-tissue spots::

        pseudobulk_cpm = 1e6 * sum(MGMT counts) / sum(total UMI)
        log2cpm        = log2(pseudobulk_cpm + 1)

    THIS IS PSEUDOBULK, NOT A RESOLVABLE SPATIAL FIELD. MGMT is detected in only
    5-29 % of spots across our four sections and each positive spot typically
    carries 1-2 UMI, so per-spot MGMT is at the detection floor -- any "MGMT map"
    built from it would be counting-noise. The honest product is one number per
    section, which is precisely the per-patient scalar this module needs.

    Parameters
    ----------
    sample_id : str
        A cached Visium GBM sample id (``'10x_parent_gbm'``, ``'GBM1'``,
        ``'GBM2'``, ``'GBM3'``); see ``tmesim.data.real.visium.list_samples()``.

    Returns
    -------
    dict with ``sample_id``, ``mgmt_in_reference``, ``n_spots``,
    ``n_spots_mgmt_positive``, ``pct_spots_mgmt_positive``, ``total_mgmt_counts``,
    ``total_umi``, ``median_umi_per_spot``, ``pseudobulk_cpm``,
    ``log2cpm``, ``sensitivity`` (the calibrated TMZ multiplier),
    ``sensitivity_ci95``, ``resolution``, ``caveat``, ``accession``.

    Raises
    ------
    FileNotFoundError / KeyError if the sample's bytes are not cached.
    """
    import numpy as np
    from ..data.real.visium import load_sample, sample_info

    data = load_sample(sample_id, marker_report=False)
    genes = list(data["genes"])
    info = sample_info(sample_id)
    if "MGMT" not in genes:
        return {
            "sample_id": sample_id,
            "mgmt_in_reference": False,
            "n_spots": int(data["n_spots"]),
            "reason": "MGMT is not in this section's feature reference",
            "accession": info.get("accession"),
        }

    col = genes.index("MGMT")
    expr = data["expression"]
    row = np.asarray(expr[:, col].todense()).ravel().astype(float)
    totals = np.asarray(data["total_counts"], dtype=float)
    total_umi = float(totals.sum())
    mgmt_counts = float(row.sum())
    cpm = 1e6 * mgmt_counts / total_umi if total_umi > 0 else 0.0
    x = math.log2(cpm + 1.0)
    det = mgmt_sensitivity(x, detail=True)

    return {
        "sample_id": sample_id,
        "mgmt_in_reference": True,
        "n_spots": int(data["n_spots"]),
        "n_spots_mgmt_positive": int((row > 0).sum()),
        "pct_spots_mgmt_positive": float(100.0 * (row > 0).mean()),
        "total_mgmt_counts": mgmt_counts,
        "total_umi": total_umi,
        "median_umi_per_spot": float(np.median(totals)),
        "mean_mgmt_counts_per_spot": float(row.mean()),
        "pseudobulk_cpm": float(cpm),
        "log2cpm": float(x),
        "sensitivity": float(det["sensitivity"]),
        "sensitivity_ci95": det["ci95"],
        "extrapolated": bool(det["extrapolated"]),
        "resolution": "PSEUDOBULK (one scalar per section)",
        "caveat": (
            "Pseudobulk, NOT a spatial field: MGMT is detected in "
            f"{100.0 * (row > 0).mean():.1f} % of spots at ~1-2 UMI each, i.e. at "
            "the detection floor, so a per-spot MGMT map would be counting noise. "
            "This is also a different assay (3' spot-based Visium) from the bulk "
            "RNA-seq the mapping was calibrated on; CPM is nominally comparable "
            "but not identically calibrated. A clinically measured MGMT would be "
            "the proper input."
        ),
        "source": info.get("source"),
        "accession": info.get("accession"),
    }


def mgmt_from_visium_all(sample_ids: Optional[Sequence[str]] = None) -> Dict[str, Dict]:
    """:func:`mgmt_from_visium` over every cached Visium GBM section."""
    from ..data.real.visium import list_samples

    ids = list(sample_ids) if sample_ids is not None else list_samples()
    return {sid: mgmt_from_visium(sid) for sid in ids}


if __name__ == "__main__":  # pragma: no cover - manual inspection
    print("MGMT -> temozolomide sensitivity  "
          f"(exp(-{MGMT_LAMBDA:.4f} * log2(CPM+1)), n = {MGMT_CALIBRATION_N})")
    print(f"  95 % CI on lambda: {MGMT_LAMBDA_CI95[0]:.4f}-{MGMT_LAMBDA_CI95[1]:.4f}"
          f"  p = {MGMT_LAMBDA_P:.4f}")
    for cpm in (0.0, 1.0, 3.0, 7.0, 15.0, 25.0):
        s = mgmt_sensitivity(cpm, units="cpm")
        print(f"  MGMT {cpm:5.1f} CPM (log2cpm+1 = {math.log2(cpm + 1):.3f})"
              f" -> sensitivity {s:.4f}")
    print(f"  binary methylated   -> {MGMT_BINARY_SENSITIVITY['methylated']:.4f}")
    print(f"  binary unmethylated -> {MGMT_BINARY_SENSITIVITY['unmethylated']:.4f}"
          f"  (95 % CI {MGMT_BINARY_CI95[0]:.3f}-{MGMT_BINARY_CI95[1]:.3f})")
    print("\nREAL Visium GBM sections (pseudobulk, NOT a spatial field):")
    for sid, rec in mgmt_from_visium_all().items():
        if not rec.get("mgmt_in_reference"):
            print(f"  {sid:16s} MGMT absent from reference")
            continue
        print(f"  {sid:16s} {rec['n_spots']:5d} spots  "
              f"{rec['pct_spots_mgmt_positive']:5.1f}% MGMT+  "
              f"pseudobulk {rec['pseudobulk_cpm']:6.2f} CPM  "
              f"-> TMZ sensitivity {rec['sensitivity']:.4f}")
