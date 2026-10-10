"""Answer a RunSpec WITHOUT solving it, in the two cases where that is exact.

There are exactly two of them, and neither is a model.

1. **An exact memo.** There is no RunSpec-keyed run cache in this backend --
   only the reconstruction cache (:data:`app.server.cache.RECON_CACHE`) -- so
   submitting a spec that has already been solved re-solves it from scratch. If
   the spec, the server-side solver constants and the engine are the same, the
   answer is the same, and it can be served from the stored record at zero
   error.

2. **Algebraic scaling along dose and half-life.** In the shipped
   (``time_resolved=False``) path every PK number a user can set reaches the
   transport solve through ONE float. ``tmesim/sim/engine.py`` computes
   ``plasma_mean = plasma_profile(pk, times).mean()`` and then
   ``source = delivery_source(work, drug, plasma_mean)``; nothing else about
   the schedule is passed down. ``plasma_profile`` is exactly linear in
   ``pk.plasma_peak``, and the transport PDE, its first-order sink and the
   steady-state reach solve are all linear in the source with ``C(0) = 0``. So
   the whole concentration field is homogeneous of degree 1 in that one
   scalar: every RATIO of exposures is invariant, and every exposure MAGNITUDE
   is proportional. That is an identity, not an approximation.

   What is NOT linear, and is refused: anything compared against a THRESHOLD
   (``therapeutic_coverage_*``), anything through the Hill kill model
   (``kill_fraction``, ``viable_tumor_fraction``), and the exponential-plus-floor
   FIT (``penetration_length_um``), whose round-off ``app/server/reach.py``'s
   ``SCALE_INVARIANCE_NOTE`` already measures at up to 2.9e-05 relative at
   dx = 30 under a pure global multiplier on the source.

THE CORRECTNESS PROBLEM THIS MODULE IS ABOUT
--------------------------------------------
A cache that serves a stale number after the engine changed is far worse than
no cache. This repository asserts bit-exact stored references and quotes
figures to ten significant figures; a silently stale hit would corrupt that
and leave no trace. So the key is the whole design, and it is composed of FIVE
parts, each answering a different way the stored answer could stop being the
answer:

``memo_format``
    This module's own key format. Bumping it invalidates every stored key, so
    a change to the canonicaliser can never be mistaken for a match.

``spec``
    The RunSpec, canonicalised. Protects against two specs being read as one.
    Canonicalisation is TOTAL (every key participates, nulls included) and
    CONSERVATIVE: it normalises only what provably cannot change the run --
    key order, and the int/float spelling of one number, because
    ``execute_run`` coerces with ``float()``. It does NOT resolve omitted
    fields to their defaults, so ``{}`` and ``{"timelapse_frames": 24}`` are
    two keys for one run. That costs a MISS, which is free. The opposite
    error -- two runs reading as one -- is the one that puts a wrong number on
    a screen, so every doubtful case is resolved towards the miss.

``solver_constants``
    ``PDE_DT_S``, ``PDE_SAVE_EVERY``, ``PDE_CONVECTION``, ``SIM_SEED``,
    ``DEFAULT_TIMELAPSE_FRAMES`` and ``TIMELAPSE_BUDGET_BYTES``. A frozen
    RunSpec carries no solver settings -- ``app/server/jobs.py`` says so where
    it fixes them -- so the spec alone does not pin the run. The last two are
    in the key because an omitted ``timelapse_frames`` resolves against the
    default and the budget can reduce the realised frame count, so both change
    what a run produces.

``engine``
    The hard part: something that changes when the engine's NUMBERS change.
    See :func:`engine_fingerprint` for what the choice does and does not
    protect against, stated plainly.

``tissue``
    A (size, mtime) digest of the files the section's metadata is read from,
    via ``app/server/cache.py``'s own ``_fingerprint``. Tissue registers by
    filesystem path in this product, so the path's contents are an input to
    the answer like any other.

``claim_set_version``
    ``app/server/jobs.py`` stamps ``CLAIM_SET_VERSION`` into every run so a
    STORED record banners itself when a claim is later withdrawn. A memo hit
    must not launder that. The version is in the KEY, so a withdrawal makes
    every memo a miss and a repeat re-solves into compliant caveats; and
    :func:`memo_lookup` ALSO re-runs ``claim_set_status`` on the candidate and
    refuses anything SUPERSEDED or CONTRADICTED. Two guards, because they
    catch different things: the key catches the era, and the live check
    catches a record that is stamped current and still carries a retired
    sentence -- which is a real specimen in this store, not a hypothetical.

WHAT THIS MODULE DOES NOT DO
----------------------------
It does not serve, scale or compare an absolute kill fraction: that is never a
defensible output here. It answers ONE spec per call and produces no ordering
between compounds; the cross-compound ranking was withdrawn on 2026-08-21 and
nothing here reinstates it. It makes no prediction of any kind -- there is no
model, nothing is fitted, and both paths return either the engine's own stored
numbers or an exact algebraic rearrangement of them.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import sys
from typing import Dict, List, Optional, Tuple

from .cache import _fingerprint, _metadata_source_files, lookup_section

__all__ = [
    "MEMO_FORMAT_VERSION",
    "Unmemoisable",
    "canonical",
    "spec_digest",
    "solver_constants",
    "engine_fingerprint",
    "tissue_fingerprint",
    "run_keys",
    "SCALED_AXES",
    "READOUT_RULES",
    "plasma_scalars",
    "payload_status",
    "memo_lookup",
    "scale_lookup",
    "answer_spec",
]

#: Bump on ANY change to the canonicalisation or to the key's composition.
#: Every stored key then stops matching, which is the safe direction: a stale
#: miss costs one solve, a stale hit costs the truth.
MEMO_FORMAT_VERSION = 1

#: Replaces the value of each scalable axis in the SHAPE key, so two specs
#: differing only along those axes land on one key and can find each other.
_AXIS_SENTINEL = "<memo:scaled-axis>"

#: The two PK fields whose change is answerable by arithmetic. Deliberately
#: NOT ``n_doses`` or ``dose_interval_s``: the derivation above says they
#: should scale too, because they also reach the solver only through
#: ``plasma_mean`` -- but every row of the measurement this module is built on
#: carries ``n_doses = 1``, so that is DERIVED and UNMEASURED, and an
#: unmeasured identity is not an exact answer. And deliberately NOT
#: ``total_time_s``, which is not a prefactor at all: it changes the
#: integration window as well as the plasma average, and is measured to move
#: ``median_tumor_auc`` by 8.93x and ``E_p10`` by 3.97x under a 6x change.
SCALED_AXES = ("plasma_peak", "half_life_s")

#: Largest int that survives ``float()`` without losing an integer step. A
#: spec carrying something bigger is REFUSED rather than canonicalised through
#: a float, which would merge two distinct integers into one key.
_EXACT_INT = 2 ** 53


class Unmemoisable(Exception):
    """This spec cannot be given a key that is safe to compare.

    Carries an authored sentence for the caller. Never raised from another
    exception's text, so nothing internal reaches the wire through it.
    """


# --------------------------------------------------------------------------- #
# 1. Canonicalisation
# --------------------------------------------------------------------------- #
def canonical(value):
    """A JSON-serialisable, order-free, type-tagged image of ``value``.

    Every tag is explicit, so ``True`` and ``1`` cannot collide and neither can
    ``None`` and ``"null"``. Numbers are tagged by their float VALUE because
    ``execute_run`` coerces every numeric spec field with ``float()``, which
    makes ``12`` and ``12.0`` the same run -- the one normalisation this
    function performs, and the only one it can justify.
    """
    if value is None:
        return ["null"]
    if isinstance(value, bool):
        # Before int: bool is an int subclass, and ``resistance_off: true`` and
        # ``resistance_off: 1`` reach ``bool()`` identically -- but keeping them
        # apart only costs a miss, so they are kept apart.
        return ["bool", value]
    if isinstance(value, int):
        if abs(value) > _EXACT_INT:
            raise Unmemoisable(
                "a RunSpec field holds an integer too large to compare "
                "exactly as a floating-point number, so this spec cannot be "
                "given a key that is safe to match. Submit it as a run."
            )
        return ["num", repr(float(value))]
    if isinstance(value, float):
        if not math.isfinite(value):
            raise Unmemoisable(
                "a RunSpec field holds a non-finite number (NaN or an "
                "infinity), which has no canonical form and no well-defined "
                "equality, so this spec cannot be given a key. Submit it as a "
                "run; /api/validate will say what is wrong with it."
            )
        return ["num", repr(value)]
    if isinstance(value, str):
        return ["str", value]
    if isinstance(value, dict):
        names = [str(k) for k in value]
        if len(set(names)) != len(names):
            raise Unmemoisable(
                "a RunSpec object holds two keys that are distinct objects "
                "with the same text, so canonicalising it would merge them. "
                "Submit it as a run."
            )
        return ["obj", [[str(k), canonical(v)]
                        for k, v in sorted(value.items(), key=lambda kv: str(kv[0]))]]
    if isinstance(value, (list, tuple)):
        return ["arr", [canonical(v) for v in value]]
    raise Unmemoisable(
        "a RunSpec field holds a value of a type this backend has no "
        "canonical form for, so it cannot be matched against a stored run. "
        "Submit it as a run."
    )


def _digest(obj) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, separators=(",", ":"),
                   allow_nan=False).encode("utf-8")
    ).hexdigest()


def spec_digest(spec: Dict) -> str:
    """The canonical digest of one RunSpec."""
    return _digest(canonical(spec))


def _shape_spec(spec: Dict) -> Dict:
    """``spec`` with both scalable axes blanked, so two specs differing ONLY
    along them share a key.

    The sentinel is INSERTED when the field is absent as well as when it is
    present, because an absent ``plasma_peak`` means the engine's own default
    and that is a value like any other. The effective values are recovered
    separately by :func:`plasma_scalars` and then CHECKED against the stored
    run's own ``plasma_mean_uM`` bit-for-bit, so a mistake in recovering them
    cannot silently become a wrong scale factor.
    """
    out = dict(spec)
    pk = dict(out.get("pk") or {})
    for axis in SCALED_AXES:
        pk[axis] = _AXIS_SENTINEL
    out["pk"] = pk
    return out


# --------------------------------------------------------------------------- #
# 2. What the answer depends on besides the spec
# --------------------------------------------------------------------------- #
def solver_constants() -> Dict:
    """The server-side settings a frozen RunSpec does not carry.

    Read live off ``app/server/jobs.py`` rather than copied, so there is one
    definition and this module cannot drift from the thing it is describing.
    """
    from . import jobs

    return {
        "pde_dt_s": float(jobs.PDE_DT_S),
        "pde_save_every": int(jobs.PDE_SAVE_EVERY),
        "pde_convection": bool(jobs.PDE_CONVECTION),
        "sim_seed": int(jobs.SIM_SEED),
        "default_timelapse_frames": int(jobs.DEFAULT_TIMELAPSE_FRAMES),
        "timelapse_budget_bytes": int(jobs.TIMELAPSE_BUDGET_BYTES),
    }


#: Which source trees are hashed into the engine fingerprint. ``tmesim/`` whole,
#: because the run path imports lazily from most of it and an import-time
#: snapshot would miss a module not yet loaded; plus the server modules that
#: compute or shape a run's numbers.
_ENGINE_TREES = ("tmesim",)
_ENGINE_FILES = (
    "app/server/jobs.py",
    "app/server/reach.py",
    "app/server/cache.py",
    "app/server/compounds.py",
    "app/server/delivery_law.py",
    "app/server/payloads.py",
    "app/server/memo.py",
    "app/server/memo_api.py",
)

_FINGERPRINT: Optional[Dict] = None


def _repo_root() -> str:
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def engine_fingerprint() -> Dict:
    """A digest that changes when the engine's numbers can change.

    WHAT IT IS. The sha256 of (a) the bytes of every ``.py`` file under
    ``tmesim/`` and of the server modules in :data:`_ENGINE_FILES`, by
    repo-relative path; (b) the live value of ``calibrate_alpha()``, to its
    full repr; and (c) the numpy, scipy and CPython versions in this process.

    WHY BOTH HALVES. ``calibrate_alpha()`` alone is the obvious candidate and
    it is not enough: it is one constant fit on one penetration anchor, and the
    PDE stepper, the kill model, the delivery law, the metric reductions and
    the reach convention can all move a run's numbers without moving it. A
    source hash alone is not enough either: it cannot see a changed
    calibration INPUT, and the calibration is data-driven. Together they cover
    the two ways the engine's output moves -- its code and its one fitted
    constant.

    WHAT IT DOES NOT PROTECT AGAINST, stated plainly because a guard whose
    gaps are unlisted reads as a guarantee:

    * **Source on disk is not bytecode in memory.** This is computed ONCE per
      process and cached, which is deliberate: the fingerprint should describe
      the code that is RUNNING, and that is whatever was imported. Edit a
      module without restarting and this digest keeps describing the modules in
      memory, which is correct; restart and it changes, which is also correct.
      What it cannot do is notice a reloader that re-imported one module in
      place.
    * **Input data.** No expression matrix, compound library file or ingest
      registry entry is hashed here. The section's metadata files are covered
      separately by :func:`tissue_fingerprint`, and that is a (size, mtime)
      check over METADATA files only -- it would not see an edit to a matrix
      inside a Visium directory. A memo could be made immune to this by
      rebuilding the reconstruction (cheap, and cached anyway) and requiring
      the stored ``report`` block to match bit-for-bit before the metrics are
      served. That is the closure for this gap and it is NOT built here.
    * **The numerical environment.** A different BLAS, a different thread
      count or a different CPU can change a reduction order. Nothing in this
      digest sees that.
    WHAT IT COSTS. ``memo.py`` and ``memo_api.py`` are in the hashed set, so
    editing a COMMENT in either invalidates every stored memo. That is a
    deliberate trade of cache longevity for correctness: a change to the
    canonicaliser or the key composition MUST invalidate, and
    :data:`MEMO_FORMAT_VERSION` only does that when a human remembers to bump
    it. Hashing the file removes the human from the loop. The cost is paid in
    re-solves, which is the safe direction; the alternative is paid in wrong
    numbers. Observed in development: an edit to this docstring turned every
    hit in a test store into a miss, which is the mechanism working.

    * **Anything outside the hashed set.** A module added under ``app/server/``
      that reaches a run's numbers and is not in :data:`_ENGINE_FILES` is
      invisible to this. ``tests/test_memo.py`` pins the list against the
      ``engine_entrypoints`` a run records, so a new entry point in the chain
      fails the suite rather than quietly weakening the key.

    So: this makes a stale hit require either an unhashed input, an unhashed
    module or a changed numerical environment. It does not make one
    impossible.
    """
    global _FINGERPRINT
    if _FINGERPRINT is not None:
        return _FINGERPRINT

    from tmesim.validation.calibration import calibrate_alpha

    root = _repo_root()
    parts: List[Tuple[str, str]] = []
    for tree in _ENGINE_TREES:
        base = os.path.join(root, tree)
        for dirpath, dirnames, filenames in os.walk(base):
            dirnames[:] = sorted(d for d in dirnames
                                 if d not in ("__pycache__",))
            for name in sorted(filenames):
                if not name.endswith(".py"):
                    continue
                full = os.path.join(dirpath, name)
                rel = os.path.relpath(full, root).replace(os.sep, "/")
                with open(full, "rb") as fh:
                    parts.append((rel, hashlib.sha256(fh.read()).hexdigest()))
    for rel in _ENGINE_FILES:
        full = os.path.join(root, rel)
        if not os.path.isfile(full):
            # A declared file that is gone is itself a change, recorded as one.
            parts.append((rel, "absent"))
            continue
        with open(full, "rb") as fh:
            parts.append((rel, hashlib.sha256(fh.read()).hexdigest()))

    import numpy
    try:
        import scipy
        scipy_version = scipy.__version__
    except ImportError:          # pragma: no cover - scipy is a hard dep
        scipy_version = "absent"

    alpha = float(calibrate_alpha())
    body = {
        "sources": sorted(set(parts)),
        "alpha": repr(alpha),
        "numpy": numpy.__version__,
        "scipy": scipy_version,
        "python": sys.version.split()[0],
    }
    _FINGERPRINT = {
        "digest": _digest(body),
        "alpha": alpha,
        "n_source_files": len(body["sources"]),
        "numpy": body["numpy"],
        "scipy": body["scipy"],
        "python": body["python"],
    }
    return _FINGERPRINT


def tissue_fingerprint(section) -> str:
    """A (size, mtime) digest of the files this section's metadata is read from.

    Tissue registers BY FILESYSTEM PATH here, so what is at the path is an
    input to the answer. Reuses ``app/server/cache.py``'s own ``_fingerprint``
    so there is one definition of what "this specimen's files changed" means.
    An unknown section digests to a distinct constant rather than to the same
    value as a known one, so a renamed or de-registered specimen cannot inherit
    another's memo. The raw digest is never returned on the wire -- it is built
    from absolute host paths; only its sha256 travels.
    """
    rec = lookup_section(str(section)) if section is not None else None
    if rec is None:
        return _digest(["unregistered-section", str(section)])
    return _digest(["section", rec.get("id"), _fingerprint(
        _metadata_source_files(rec))])


def run_keys(spec: Dict) -> Dict:
    """``{"run_key", "shape_key", "components"}`` for one RunSpec.

    ``run_key`` matches a stored run that is the SAME run. ``shape_key``
    matches a stored run that differs at most along :data:`SCALED_AXES`.
    ``components`` is what the keys are made of, with every host path already
    reduced to a digest so it is safe to return to a client.
    """
    from .jobs import CLAIM_SET_VERSION

    engine = engine_fingerprint()
    components = {
        "memo_format": MEMO_FORMAT_VERSION,
        "spec_digest": spec_digest(spec),
        "solver_constants": solver_constants(),
        "engine": engine,
        "tissue_digest": tissue_fingerprint(spec.get("section")),
        "claim_set_version": int(CLAIM_SET_VERSION),
    }
    common = {k: v for k, v in components.items() if k != "spec_digest"}
    return {
        "run_key": _digest([canonical(common), components["spec_digest"]]),
        "shape_key": _digest([canonical(common),
                              spec_digest(_shape_spec(spec)), "shape"]),
        "components": components,
    }


# --------------------------------------------------------------------------- #
# 3. The algebra: which readouts are invariant, which scale, and by what
# --------------------------------------------------------------------------- #
#: ``{block: {leaf: (a, b, max_measured_ulp, note)}}`` -- the served value is
#: ``stored * mean_ratio**a * peak_ratio**b``, where ``mean_ratio`` is the
#: engine's own ``plasma_mean`` ratio and ``peak_ratio`` the plasma-peak ratio.
#:
#: ``(0, 0)`` is an INVARIANT and is copied through. ``max_measured_ulp`` is
#: the worst deviation this rule showed against a real re-solve in
#: ``experiments/exact_memo/probe_scaling.py`` -- 5 (section, compound) cases x
#: 15 dose arms and 10 half-life arms, at dose ratios 4 (a power of two) and
#: 1.5 and 3 (not), and half-life ratios 4 and 1.37. 0 means BIT-IDENTICAL in
#: every arm of every case.
#:
#: EVERY leaf a run reports and this table does not name is REFUSED. That is
#: the fail-closed half and it is the point: a metric added to ``execute_run``
#: is unmeasured until somebody measures it, and unmeasured is not invariant.
READOUT_RULES: Dict[str, Dict[str, Tuple[int, int, float, str]]] = {
    "metrics": {
        "mean_tumor_auc": (1, 0, 2.0,
                           "mean tumour-cell exposure: linear in the source"),
    },
    "extra_metrics": {
        # --- exact invariants (bit-identical in every measured arm) -------- #
        "exposure_window_s": (0, 0, 0.0, "the window, which neither axis moves"),
        "hypoxic_tumor_fraction": (0, 0, 0.0,
                                   "banded from the oxygen field, before any drug"),
        "n_tumor": (0, 0, 0.0, "a cell count from the reconstruction"),
        "penetration_length_analytic_um": (
            0, 0, 0.0, "closed form in the compound and the ECM; no PK input"),
        # --- ratios: invariant because numerator and denominator scale alike #
        "low_o2_over_high_o2_exposure": (
            0, 0, 2.0, "a ratio of two exposures over the same field"),
        "penetration_index": (
            0, 0, 34.0, "a ratio of exposures; round-off only"),
        "tissue_median_over_plasma_mean": (
            0, 0, 35.0, "E_p50 and plasma_mean scale together"),
        # --- magnitudes: proportional to the engine's own PK scalar -------- #
        "mean_tissue_auc": (1, 0, 2.0, "mean tissue exposure"),
        "median_tumor_auc": (1, 0, 53.0, "a percentile of the exposure field"),
        "E_p10_uM": (1, 0, 27.0, "a percentile of the exposure field"),
        "E_p50_uM": (1, 0, 54.0, "a percentile of the exposure field"),
        "E_p90_uM": (1, 0, 21.0, "a percentile of the exposure field"),
        "median_E_over_ic50": (1, 0, 53.0, "E_p50 over a fixed IC50"),
        "p10_E_over_ic50": (1, 0, 27.0, "E_p10 over a fixed IC50"),
        "plasma_mean_uM": (1, 0, 0.0, "the scalar itself"),
        # --- driven by the PEAK rather than by the window average ---------- #
        "plasma_peak_uM": (0, 1, 0.0, "the prescribed peak"),
        "plasma_peak_over_ic50": (0, 1, 0.0, "the peak over a fixed IC50"),
        "tissue_median_over_plasma_peak": (
            1, -1, 17.0, "E_p50 scales with the window average, the "
                         "denominator with the peak"),
    },
}

#: Why each refused readout is refused, in the words a caller needs. Measured
#: fold changes are from the same probe. A leaf absent from BOTH tables is
#: refused with a sentence saying it was never measured.
REFUSED_READOUTS: Dict[str, str] = {
    "kill_fraction": (
        "REFUSED, for two independent reasons. An absolute kill fraction is "
        "never a defensible output of this product, so it is not a scaling "
        "feature. And it is not scalable even in principle: it is a Hill "
        "function of exposure, measured to move 1.0021x-1.9959x under dose "
        "changes of 1.5x-4x, i.e. neither invariant nor proportional. Re-solve."
    ),
    "therapeutic_coverage_1x": (
        "REFUSED: coverage counts cells whose exposure crosses a FIXED IC50 "
        "threshold, so it is not homogeneous in the source. Measured to move "
        "up to 8.1109x under a 4x dose. Re-solve."
    ),
    "therapeutic_coverage_0_5x": (
        "REFUSED: a fixed-threshold count, not homogeneous in the source. "
        "Measured up to 1.3185x under a 4x dose. Re-solve."
    ),
    "therapeutic_coverage_2x": (
        "REFUSED: a fixed-threshold count, not homogeneous in the source. "
        "Measured up to 3.4403x under a 4x dose. Re-solve."
    ),
    "viable_tumor_fraction": (
        "REFUSED: it is a post-kill cell count, so it inherits the kill "
        "model's non-linearity. Measured to move to 0.2111x under a 4x dose. "
        "Re-solve."
    ),
    "penetration_length_um": (
        "REFUSED. The underlying length IS invariant in exact arithmetic, but "
        "it is read off an exponential-plus-floor LEAST-SQUARES FIT, which "
        "amplifies round-off: app/server/reach.py's SCALE_INVARIANCE_NOTE "
        "measures this same fitted length moving 2.9e-05 relative at dx = 30 "
        "under a pure global multiplier on the source. A readout whose "
        "agreement with a re-solve depends on the grid is not an exact "
        "answer. Re-solve."
    ),
}

#: Non-numeric entries that ARE reproduced, because they have no PK input at
#: all by construction rather than by measurement.
#:
#: ``verdict_stability`` is ``verdict_stability(drug_name)`` -- a lookup in the
#: verdict-inversion census keyed on the compound name and nothing else. It
#: says whether THIS compound keeps its delivery/dose verdict when the
#: transvascular exchange term the transport law lacks is switched on, and it
#: has to travel with any delivery number read off this answer: a reader handed
#: a scaled exposure without it would be reading one of the seven compounds in
#: twelve whose verdict inverts and have nothing saying so.
COPIED_VERBATIM: Tuple[str, ...] = ("verdict_stability",)

#: Whole blocks a scaled answer does not reproduce, and why.
REFUSED_BLOCKS: Dict[str, str] = {
    "mechanism": (
        "REFUSED: mechanism attribution is a partition of the KILL outcome, "
        "which does not scale. It also has a binary payload, and a scaled "
        "answer has no payload of its own."
    ),
    "timelapse": (
        "REFUSED: a time-lapse is a stored binary payload. A scaled answer "
        "produces no file, and pointing at the base run's would hand over "
        "frames computed at the other dose."
    ),
}


def plasma_scalars(pk: Dict, total_time_s: float) -> Dict:
    """``plasma_peak`` and the window-average ``plasma_mean`` for one PK block.

    Computed by calling the engine's OWN ``plasma_profile`` on the engine's own
    time grid, so the factor below is the engine's number and not this module's
    reading of a formula. ``total_time`` and ``dt`` are taken from the spec and
    from ``jobs.PDE_DT_S``, exactly as ``execute_run`` assembles them.
    """
    import numpy as np

    from tmesim.core import PKConfig
    from tmesim.treatment.drugs import plasma_profile

    from .jobs import PDE_DT_S

    default = PKConfig()
    cfg = PKConfig(
        plasma_peak=float((pk or {}).get("plasma_peak", default.plasma_peak)),
        half_life=float((pk or {}).get("half_life_s", default.half_life)),
        n_doses=int((pk or {}).get("n_doses", default.n_doses)),
        dose_interval=float((pk or {}).get("dose_interval_s",
                                           default.dose_interval)),
    )
    times = np.arange(0.0, float(total_time_s), float(PDE_DT_S))
    profile = plasma_profile(cfg, times)
    return {
        "plasma_peak": float(cfg.plasma_peak),
        "plasma_mean": float(profile.mean()) if profile.size else 0.0,
        "n_doses": int(cfg.n_doses),
        "dose_interval_s": float(cfg.dose_interval),
        "half_life_s": float(cfg.half_life),
    }


def _is_power_of_two(x: float) -> bool:
    """Is multiplying a double by ``x`` an exact operation?

    The one case in which a scaled magnitude is BIT-IDENTICAL to a re-solve
    rather than merely equal to within round-off. The study this module is
    built on swept the dose at a ratio of exactly 4.0 and reported "exactly
    4.0"; that is true, and it is true BECAUSE 4 is a power of two. At 1.5x
    the same readouts land 1-54 ulp away.
    """
    if not (math.isfinite(x) and x > 0.0):
        return False
    m, e = math.frexp(x)
    del e
    return m == 0.5


# --------------------------------------------------------------------------- #
# 4. Is a stored run safe to re-serve?
# --------------------------------------------------------------------------- #
#: Every binary payload a completed run can hold, by the name the API serves it
#: under. Derived from ``payloads.FIELD_NAMES`` so a new field joins the check
#: automatically.
def _payload_names() -> Tuple[str, ...]:
    from .payloads import FIELD_NAMES

    return ("geometry",) + tuple(f"field_{n}" for n in FIELD_NAMES) + (
        "mechanism", "timelapse")


def payload_status(mgr, run_id: str, body: Dict) -> Dict:
    """Which of this run's payloads are ACTUALLY on disk, and which are not.

    Three states, and the distinction is the whole point:

    ``present``
        the run produced it and the file is there.
    ``pruned``
        the run produced it and ``scripts/prune_run_payloads.py`` deleted it,
        with a date and a byte count in the record.
    ``missing``
        the run produced it and the file is gone with no prune recorded --
        a bare ``rm``, a partial copy, a failed restore.
    ``not_produced``
        the run never had one.

    ``missing`` exists because a prune RECORD is not the same as the file's
    absence, and the prune script's own docstring is about exactly that gap:
    ``rm`` reclaims the disk and leaves the server explaining confidently that
    the payload was never captured. So this stats the filesystem instead of
    trusting the record, and a memo hit reports what it found.
    """
    pruned = body.get("pruned_payloads") or {}
    pruned = pruned if isinstance(pruned, dict) else {}
    drug = bool((body.get("spec") or {}).get("drug"))
    out: Dict[str, Dict] = {}
    for name in _payload_names():
        path = mgr.payload_path(run_id, name)
        on_disk = os.path.isfile(path)
        if on_disk:
            state = "present"
        elif name in pruned:
            state = "pruned"
        else:
            # Did this run ever have one? ``geometry`` and the reconstruction
            # fields always; ``field_drug_auc`` only with a drug; mechanism and
            # timelapse only when the run says they are available.
            if name == "mechanism":
                produced = bool(body.get("mechanism"))
            elif name == "timelapse":
                produced = bool((body.get("timelapse") or {}).get("available"))
            elif name == "field_drug_auc":
                produced = drug
            else:
                produced = True
            state = "missing" if produced else "not_produced"
        rec: Dict = {"state": state,
                     "available": state == "present",
                     "endpoint": f"/api/runs/{run_id}/{name}"}
        if state == "pruned":
            entry = pruned.get(name) or {}
            rec["reason"] = (
                "this run DID produce it; the file was deleted on "
                f"{entry.get('at', 'an unrecorded date')} to reclaim disk "
                f"({int(entry.get('bytes', 0)):,} bytes). Re-run the RunSpec "
                "to reproduce it."
            )
        elif state == "missing":
            rec["reason"] = (
                "this run produced it and the file is NOT on disk, with no "
                "prune recorded. Treat this memo hit as incomplete and "
                "re-run the RunSpec if you need the payload."
            )
        elif state == "not_produced":
            rec["reason"] = "this run never produced it."
        out[name] = rec
    complete = all(r["state"] in ("present", "not_produced")
                   for r in out.values())
    return {"complete": bool(complete), "payloads": out}


# --------------------------------------------------------------------------- #
# 5. The two lookups
# --------------------------------------------------------------------------- #
def _done_runs(mgr) -> List:
    """Completed runs, NEWEST FIRST -- which is the order ``mgr.list()`` gives.

    When two stored runs carry the same key they are the same execution twice
    and either answers, so the choice is free and the newest is taken: it is
    the least likely to have had its payloads pruned, and the oldest copy of a
    duplicated spec is the one a reader is least likely to be looking at.
    """
    return [j for j in mgr.list()
            if j.status == "done" and isinstance(j.body, dict)]


def _stored_key(job, field: str) -> Optional[str]:
    prov = (job.body or {}).get("provenance") or {}
    memo = prov.get("memo") or {}
    if not isinstance(memo, dict):
        return None
    if int(memo.get("format") or 0) != MEMO_FORMAT_VERSION:
        # A run stamped by an older key format is never a candidate. Every run
        # written before this module existed has no stamp at all and is
        # likewise never a candidate -- which is the only honest answer: its
        # engine state was not recorded, so it cannot be shown to match.
        return None
    value = memo.get(field)
    return str(value) if value else None


def _claim_set_gate(job) -> Optional[str]:
    """``None`` when this record may be re-served, else why it may not."""
    from .jobs import claim_set_status

    prov = (job.body or {}).get("provenance") or {}
    status = claim_set_status(prov.get("claim_set"), prov.get("caveats"))
    if status.get("contradicted"):
        return (
            "the stored run matching this spec CONTRADICTS the claim set it is "
            "stamped with: its caveats still contain a claim that claim set "
            "had already withdrawn. It is served by /api/runs under a banner "
            "saying so, and it is not served from here at all. Re-run the "
            "spec to obtain a compliant record."
        )
    if status.get("superseded"):
        return (
            "the stored run matching this spec was written under a SUPERSEDED "
            "claim set, so re-serving it here would put a withdrawn claim in "
            "front of you without the banner /api/runs attaches. Re-run the "
            "spec."
        )
    return None


def memo_lookup(mgr, spec: Dict) -> Dict:
    """The exact memo. Either a stored run that IS this spec's answer, or a miss.

    Never approximate and never partial: a hit is the whole stored record,
    served through ``Job.result_view()`` so the claim-set banner and the
    pruned-payload reconciliation are computed on THIS read exactly as they are
    for ``GET /api/runs/{id}``.
    """
    keys = run_keys(spec)
    # The key match is a sha256 comparison. For the SPEC half that is not
    # necessary -- the stored run carries its own spec -- so it is checked for
    # EQUALITY as well, and the hash function is taken out of the trust chain
    # for the one component where it can be. It cannot be for the others: an
    # engine fingerprint is a digest by necessity.
    want = canonical(spec)
    candidates = [j for j in _done_runs(mgr)
                  if _stored_key(j, "run_key") == keys["run_key"]
                  and canonical(j.spec or {}) == want]
    refusal = None
    best = None
    for job in candidates:
        gate = _claim_set_gate(job)
        if gate:
            refusal = refusal or gate
            continue
        payloads = payload_status(mgr, job.run_id, job.body or {})
        # Two stored runs on one key are the same execution twice, so the
        # choice is free on the numbers and is made on the PAYLOADS: a copy
        # whose binaries are still on disk spares the caller more than one
        # whose have been pruned. ``_done_runs`` is newest-first, so among
        # equally complete copies the newest wins.
        if best is None or (payloads["complete"] and not best[1]["complete"]):
            best = (job, payloads)
        if best[1]["complete"]:
            break
    if best is not None:
        job, payloads = best
        return {
            "mode": "exact_memo",
            "keys": keys,
            "source_run_id": job.run_id,
            "result": job.result_view(),
            "payloads": payloads,
        }
    if refusal:
        return {"mode": "miss", "reason": refusal, "keys": keys}
    return {
        "mode": "miss",
        "reason": (
            "no stored run carries this spec's key under the current solver "
            "constants, engine fingerprint and claim set. Submit it to "
            "POST /api/runs."
        ),
        "keys": keys,
    }


def _axis_delta(base_spec: Dict, spec: Dict) -> Dict:
    """Which of the two scalable axes actually differ, by effective value."""
    base = plasma_scalars(base_spec.get("pk") or {},
                          float(base_spec["total_time_s"]))
    new = plasma_scalars(spec.get("pk") or {}, float(spec["total_time_s"]))
    changed = [a for a, k in (("plasma_peak", "plasma_peak"),
                              ("half_life_s", "half_life_s"))
               if base[k] != new[k]]
    return {"base": base, "new": new, "changed": changed}


def scale_lookup(mgr, spec: Dict) -> Dict:
    """The algebraic answer. NOT A SOLVE, and the payload says so.

    Finds a stored run differing from ``spec`` at most along
    :data:`SCALED_AXES`, then rearranges its readouts by the identity the
    module docstring derives. Refuses everything the identity does not cover.
    """
    keys = run_keys(spec)
    if spec.get("time_resolved"):
        return {"mode": "miss", "keys": keys, "reason": (
            "this spec asks for time_resolved transport, where the plasma "
            "curve drives the solver step by step instead of through its "
            "window average. A half-life change is then not a prefactor at "
            "all, and the dose identity on that path has not been measured "
            "here. Submit it to POST /api/runs."
        )}
    if spec.get("kind") == "protocol":
        return {"mode": "miss", "keys": keys, "reason": (
            "this is a TPL protocol sweep, not a RunSpec; it has no single PK "
            "block to scale. Submit it to POST /api/runs."
        )}
    # Rank the candidates by how FEW axes they differ along, newest first
    # within a tie. Every candidate is equally exact -- the identity does not
    # care how far the dose moved -- so the ranking is about the ANSWER a
    # reader gets: a one-axis derivation is shorter to check, and a dose-only
    # one can be bit-exact where a two-axis one never is. Making it a ranking
    # rather than "the first match" also makes the answer deterministic as the
    # store fills up.
    refusal = None
    ranked = []
    for index, job in enumerate(_done_runs(mgr)):
        if _stored_key(job, "shape_key") != keys["shape_key"]:
            continue
        base_spec = job.spec or {}
        # Same hardening as in memo_lookup: the shape key is a digest, and the
        # two specs are both in hand, so the thing the key ASSERTS is checked
        # for equality. Everything outside the two scaled axes must be
        # identical, or the candidate is not a dose change of this spec.
        if canonical(_shape_spec(base_spec)) != canonical(_shape_spec(spec)):
            continue
        if canonical(base_spec) == canonical(spec):
            continue       # the same spec: that is the exact memo's business
        gate = _claim_set_gate(job)
        if gate:
            refusal = refusal or gate
            continue
        try:
            n_changed = len(_axis_delta(base_spec, spec)["changed"])
        except (KeyError, TypeError, ValueError):
            continue
        if n_changed == 0:
            # Same shape key and both axes resolving to the same effective
            # values: this is the SAME run spelled differently (an omitted
            # field against its explicit default). It is not a scaling
            # question, and it must not be ranked first and then refused --
            # that would hide a genuine one-axis candidate behind it. Skipped,
            # and the miss below says what to do.
            continue
        ranked.append((n_changed, index, job, base_spec))
    if ranked:
        _, _, job, base_spec = min(ranked, key=lambda r: (r[0], r[1]))
        return _scale_from(job, base_spec, spec, keys)
    if refusal:
        return {"mode": "miss", "reason": refusal, "keys": keys}
    return {
        "mode": "miss",
        "keys": keys,
        "reason": (
            "no stored run differs from this spec ONLY in pk.plasma_peak "
            "and/or pk.half_life_s under the current solver constants, engine "
            "fingerprint and claim set. Any other difference -- including "
            "n_doses, dose_interval_s and total_time_s -- is not scalable and "
            "is deliberately not approximated. Submit it to POST /api/runs."
        ),
    }


def _scale_from(job, base_spec: Dict, spec: Dict, keys: Dict) -> Dict:
    body = job.body or {}
    delta = _axis_delta(base_spec, spec)
    stored_mean = ((body.get("extra_metrics") or {}).get("plasma_mean_uM"))

    # THE SELF-CHECK. The scale factor is a ratio of two numbers this module
    # computed; the base one must be the number the stored run actually used,
    # BIT FOR BIT. If it is not, then the PK model, the time grid or the
    # default resolution has moved under the record, and the factor would be
    # wrong in a way nothing downstream could detect. Refuse instead.
    if stored_mean is None or float(stored_mean) != delta["base"]["plasma_mean"]:
        return {"mode": "miss", "keys": keys, "reason": (
            "the window-average plasma level recomputed for the stored run's "
            "own PK block does not match the value that run recorded, so the "
            "scale factor cannot be trusted. Submit this spec to "
            "POST /api/runs."
        )}
    if not delta["changed"]:
        return {"mode": "miss", "keys": keys, "reason": (
            "this spec differs from the stored run somewhere other than "
            "pk.plasma_peak or pk.half_life_s. Submit it to POST /api/runs."
        )}
    if delta["base"]["plasma_mean"] <= 0.0 or delta["base"]["plasma_peak"] <= 0.0:
        return {"mode": "miss", "keys": keys, "reason": (
            "the stored run's plasma level is zero, so there is no ratio to "
            "scale by. Submit this spec to POST /api/runs."
        )}

    mean_ratio = delta["new"]["plasma_mean"] / delta["base"]["plasma_mean"]
    peak_ratio = delta["new"]["plasma_peak"] / delta["base"]["plasma_peak"]

    readouts: Dict[str, Dict] = {}
    refused: Dict[str, str] = {}
    copied: Dict[str, object] = {}
    worst_ulp = 0.0
    for block, rules in READOUT_RULES.items():
        stored = body.get(block) or {}
        if not isinstance(stored, dict):
            continue
        out: Dict[str, Dict] = {}
        for leaf, value in stored.items():
            rule = rules.get(leaf)
            if rule is None:
                if leaf in COPIED_VERBATIM:
                    copied[f"{block}.{leaf}"] = value
                    continue
                refused[f"{block}.{leaf}"] = REFUSED_READOUTS.get(leaf) or (
                    "REFUSED: this readout's response to a dose or half-life "
                    "change has not been measured, so it is not scaled. "
                    "Re-solve."
                    if isinstance(value, (int, float)) and
                    not isinstance(value, bool) else
                    "REFUSED: this is not a single number, and the algebra "
                    "below rearranges numbers. It is not reproduced here "
                    "even if it would be unchanged. Re-solve, or read it "
                    "off the run named in source_run_id."
                )
                continue
            a, b, ulp, note = rule
            if value is None or not isinstance(value, (int, float)) or \
                    isinstance(value, bool) or not math.isfinite(float(value)):
                # The stored run declined to produce it; a scaled answer must
                # not invent one. Its own reason travels in metric_reasons.
                out[leaf] = {"value": None, "basis": "stored_null",
                             "note": "the stored run reported no value here; "
                                     "see its metric_reasons."}
                continue
            factor = (mean_ratio ** a) * (peak_ratio ** b)
            out[leaf] = {
                "value": float(value) * factor if factor != 1.0
                else float(value),
                "factor": factor,
                "exponents": {"plasma_mean": a, "plasma_peak": b},
                "basis": "invariant" if (a, b) == (0, 0) else "proportional",
                "max_measured_ulp": ulp,
                "note": note,
            }
            worst_ulp = max(worst_ulp, ulp)
        readouts[block] = out
    for block, reason in REFUSED_BLOCKS.items():
        if body.get(block) is not None:
            refused[block] = reason

    bit_exact = (not delta["changed"] or
                 (delta["changed"] == ["plasma_peak"]
                  and _is_power_of_two(peak_ratio)
                  and _is_power_of_two(mean_ratio)))
    return {
        "mode": "algebraic_scale",
        "keys": keys,
        "source_run_id": job.run_id,
        "derivation": {
            "axes_changed": delta["changed"],
            "plasma_peak_from": delta["base"]["plasma_peak"],
            "plasma_peak_to": delta["new"]["plasma_peak"],
            "half_life_s_from": delta["base"]["half_life_s"],
            "half_life_s_to": delta["new"]["half_life_s"],
            "plasma_mean_from": delta["base"]["plasma_mean"],
            "plasma_mean_to": delta["new"]["plasma_mean"],
            "plasma_mean_ratio": mean_ratio,
            "plasma_peak_ratio": peak_ratio,
            "bit_exact": bool(bit_exact),
            "identity": (
                "In the shipped path the entire PK block reaches the solver as "
                "one float, plasma_mean = plasma_profile(pk, t).mean(), and "
                "the transport chain is linear in it with C(0) = 0. So every "
                "exposure ratio is invariant and every exposure magnitude is "
                "proportional. plasma_profile is exactly linear in "
                "plasma_peak, which is why a dose-only change factors as the "
                "peak ratio."
            ),
            "exactness": (
                "BIT-IDENTICAL to a re-solve: the change is dose-only and the "
                "ratio is a power of two, so the floating-point scaling is "
                "exact."
                if bit_exact else
                "EXACT IN EXACT ARITHMETIC, and in floating point agreeing "
                f"with a re-solve to at most {worst_ulp:g} units in the last "
                "place on the readouts above (measured, per readout, in "
                "max_measured_ulp). It is bit-identical only when the change "
                "is dose-only and the ratio is a power of two -- which is why "
                "the 4x dose sweep in "
                "results/surrogate_feasibility_2026-10-09.md reads exactly "
                "4.0 and a 1.5x change does not."
            ),
            "measured_by": "experiments/exact_memo/probe_scaling.py",
        },
        "readouts": readouts,
        "copied_verbatim": copied,
        "refused": refused,
        "report": body.get("report"),
        "drug_params": body.get("drug_params"),
        "reach": body.get("reach"),
        "reach_basis": (
            "copied unchanged. app/server/reach.py solves the steady-state "
            "field from a UNIT source (delivery_source(micro, drug, 1.0)), so "
            "the whole reach block -- E_ref, the headline and all seven rungs "
            "-- carries no PK input at all. That is structural, not a measured "
            "invariance, and it is why the ladder is bit-identical under a "
            "simultaneous dose, half-life and window change."
        ),
        "metric_reasons": body.get("metric_reasons"),
        "provenance": (job.result_view().get("provenance")),
    }


def answer_spec(mgr, spec: Dict, allow_scaling: bool = False) -> Dict:
    """Answer ``spec`` from the store if that is exact; otherwise say so.

    Always reports ``is_solve: False``. Scaling is OPT-IN: a caller that asked
    for an exact memo and silently received an algebraic rearrangement would
    have been told something it did not ask about, and the distinction is the
    one thing this surface must never blur.
    """
    hit = memo_lookup(mgr, spec)
    if hit["mode"] == "exact_memo":
        return _envelope(spec, hit, allow_scaling)
    if not allow_scaling:
        reason = hit.get("reason", "")
        return _envelope(spec, {
            "mode": "miss", "keys": hit["keys"],
            "reason": reason + (
                " Algebraic scaling along dose and half-life was NOT "
                "attempted: pass allow_scaling=true to ask for it, and read "
                "what comes back as a rearrangement of another run rather "
                "than as a solve."
            ),
        }, allow_scaling)
    scaled = scale_lookup(mgr, spec)
    if scaled["mode"] == "algebraic_scale":
        return _envelope(spec, scaled, allow_scaling)
    return _envelope(spec, {
        "mode": "miss", "keys": hit["keys"],
        "reason": hit.get("reason", "") + " " + scaled.get("reason", ""),
    }, allow_scaling)


#: One sentence per mode, for a consumer with one place to put it.
_DISCLOSURE = {
    "exact_memo": (
        "NOT A SOLVE: this is the stored record of run {run_id}, recalled "
        "because its RunSpec, the server-side solver constants, the engine "
        "fingerprint and the claim-set version all match. The numbers are the "
        "engine's own and are bit-identical to what a re-solve would print. "
        "The run_id, the timestamps and the payload endpoints belong to that "
        "ORIGINAL execution, not to a new one."
    ),
    "algebraic_scale": (
        "NOT A SOLVE AND NOT A STORED RUN: these numbers are run {run_id}'s, "
        "rearranged by an exact algebraic identity for a change of dose "
        "and/or half-life. Readouts the identity does not cover are listed "
        "under 'refused' and were NOT estimated. There is no new run record "
        "and no payload."
    ),
    "miss": (
        "NOT ANSWERED: nothing in the store answers this spec exactly, and "
        "nothing was approximated. Submit it to POST /api/runs."
    ),
}


def _envelope(spec: Dict, inner: Dict, allow_scaling: bool) -> Dict:
    mode = inner["mode"]
    out = {
        "mode": mode,
        # WHY THIS FIELD EXISTS. A memo hit and a fresh solve agree on every
        # NUMBER -- that is the claim, and tests/test_memo.py asserts it
        # bit-for-bit. They do not agree on provenance: a hit carries the
        # original execution's run_id and timestamps, and its payload
        # endpoints point at that run's files. A client that recorded a hit as
        # a new run would mint a second citation to one execution, which is
        # exactly the kind of drift this repository keeps finding in itself.
        # So the difference is disclosed rather than hidden, and it is
        # disclosed as "this is not a solve" rather than as "cached", because
        # "cached" invites the reading that the number might be stale -- and
        # the whole key exists to make that impossible.
        "is_solve": False,
        "scaling_requested": bool(allow_scaling),
        "spec": spec,
        "key": inner["keys"]["run_key"],
        "shape_key": inner["keys"]["shape_key"],
        "key_components": inner["keys"]["components"],
        "disclosure": _DISCLOSURE[mode].format(
            run_id=inner.get("source_run_id") or "?"),
    }
    out.update({k: v for k, v in inner.items() if k not in ("mode", "keys")})
    return out
