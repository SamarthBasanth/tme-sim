"""account.py — the two routes behind Studio's Account panel.

WHY THIS IS A SEPARATE MODULE AND NOT SIX MORE LINES IN ``api.py``
------------------------------------------------------------------
``api.py`` is 2,500 lines and every route in it answers a question about
TISSUE. These two answer questions about the DEPLOYMENT: who is calling, what
this instance will refuse, what it has done, and who pays for it. Keeping them
apart means the account surface can be read in one screen, which matters
because it is the surface most likely to grow a claim nobody checked.

``build_account_router`` takes its limits as an argument rather than importing
them from ``api.py``: ``api.py`` imports THIS module to mount the router, so
importing back would be a cycle. Passing them in also means the numbers on
``/api/limits`` are physically the same objects the refusal path uses — they
cannot drift into a second, stale copy, which is this project's most frequent
defect (the document that was true when it was written).

BOTH ROUTES ARE GATED. They live under ``/api/`` so the middleware in
``api.py`` requires ``X-API-Key`` like every other data route. That is
deliberate: a fingerprint, a usage count and a deployment's refusal thresholds
are operational detail about a machine. The "how do I get a key" copy that an
UNAUTHENTICATED visitor needs lives in the client (``app/web/js/ui/apikey.js``),
because by definition that visitor cannot call this.
"""

from __future__ import annotations

import hashlib
import hmac
import os
import time
from typing import Dict, Optional

from fastapi import APIRouter, Request

#: Domain-separation label for the key fingerprint. NOT a secret and not a
#: salt in the security sense -- it exists so this digest can never collide
#: with a digest of the same key computed for some other purpose elsewhere.
_FP_LABEL = b"vivum-studio-api-key-fingerprint-v1"

#: Fingerprint length in hex characters. 12 hex = 48 bits, which is far more
#: than "tell two keys apart in a UI" needs and far less than anything that
#: would help an attacker who already has to guess a 128-bit key.
_FP_HEX = 12


def key_fingerprint(key: str) -> Optional[str]:
    """A short, stable, non-reversing identifier for an API key.

    Returns ``None`` for an empty key so a caller cannot accidentally render a
    fingerprint of "".

    HMAC rather than a bare ``sha256(key)``: a bare digest of a low-entropy
    secret is a dictionary attack, and while the keys this server generates are
    128-bit, an operator-supplied ``TMESIM_STUDIO_API_KEY`` can be anything at
    all -- including a word. The label makes a precomputed table useless.

    This is an IDENTIFIER, not a verifier. Nothing authenticates against it;
    the gate in ``api.py`` still compares the full key with
    ``secrets.compare_digest``.
    """
    key = (key or "").strip()
    if not key:
        return None
    digest = hmac.new(_FP_LABEL, key.encode("utf-8"), hashlib.sha256).hexdigest()
    return digest[:_FP_HEX]


def _key_source() -> Optional[str]:
    """Where the configured key came from, or ``None`` if there isn't one.

    ``__main__.py`` sets ``TMESIM_STUDIO_API_KEY`` itself when a non-loopback
    ``--host`` is requested without an explicit ``--api-key``, and records
    which of its three sources won in ``TMESIM_STUDIO_API_KEY_SOURCE``. The
    distinction that matters is between "the operator chose this" and "the
    process invented this at boot and printed it exactly once" -- an operator
    who has LOST the key needs to know which of those they are in, and the
    startup banner that says so is scrollback they no longer have.

    An unset source variable with a key present is reported as ``environment``
    and that is the literal truth: the variable was already in the process
    environment when the app was built (Docker ENV, systemd, a shell export)
    and ``__main__.py`` never ran. An UNRECOGNISED value is reported as
    ``unknown`` rather than mapped to a plausible guess.
    """
    if not os.environ.get("TMESIM_STUDIO_API_KEY", "").strip():
        return None
    raw = os.environ.get("TMESIM_STUDIO_API_KEY_SOURCE", "").strip()
    if not raw:
        return "environment"
    return {
        "generated": "generated-at-startup",
        "flag": "--api-key",
        "environment": "environment",
    }.get(raw, "unknown")


def build_account_router(
    *,
    jobs,
    store,
    reports,
    limits: Dict,
    cache=None,
    started_at: Optional[float] = None,
) -> APIRouter:
    """Mount ``GET /api/limits`` and ``GET /api/account``.

    ``limits`` is the live limits dict from ``api.py``; see the module
    docstring for why it is passed rather than imported.
    """
    router = APIRouter()
    boot = float(started_at) if started_at is not None else time.time()

    @router.get("/api/limits")
    def get_limits() -> Dict:
        """What this instance will refuse, BEFORE you spend a run finding out.

        The Studio spec editor reads this so a spec that the server would
        reject is visibly rejected while the user is still typing. Every value
        here is the same object the refusal path in ``validate_spec`` reads --
        see ``build_account_router``'s docstring.

        ``grid_cells_formula`` is stated in words because the client has to
        reimplement it: the client knows the section extent and ``dx_um`` and
        cannot call ``validate_spec`` on every keystroke. Two implementations
        of one formula CAN drift; this field is where the client's copy is
        checked against the server's intent, and
        ``tests/test_app_server.py::test_limits_formula_matches_validate``
        checks that the numbers agree on a real section.
        """
        return dict(limits)

    @router.get("/api/account")
    def get_account(request: Request) -> Dict:
        """Who is calling, what this deployment is, what it has done, who pays.

        Reached only with a valid key (the middleware in ``api.py`` gates every
        ``/api/`` route), so ``auth.matches`` is always true when a caller sees
        this. It is reported anyway rather than assumed, because the day this
        route stops being gated is the day the field starts mattering, and a
        field that is always ``true`` is at least a field that is WRONG
        visibly rather than absent silently.
        """
        configured = os.environ.get("TMESIM_STUDIO_API_KEY", "").strip()
        presented = request.headers.get("x-api-key", "")

        by_status: Dict[str, int] = {}
        for j in jobs.list():
            by_status[j.status] = by_status.get(j.status, 0) + 1

        try:
            n_reports = len(reports.list())
        except Exception:  # noqa: BLE001 - a report store fault is not an auth fault
            n_reports = None

        return {
            "contract": "v1",
            "auth": {
                "required": bool(configured),
                "header": "X-API-Key",
                "scheme": "shared-secret",
                "source": _key_source(),
                # Fingerprint of the SERVER's key and of the CALLER's key,
                # separately. They are equal here by construction; showing both
                # is what lets a user with two saved keys tell which one this
                # browser actually sent, which is the question the panel exists
                # to answer.
                "key_fingerprint": key_fingerprint(configured),
                "caller_fingerprint": key_fingerprint(presented),
                "matches": bool(configured) and presented == configured,
                # There is no sign-up, no user table, and no way to issue a
                # second key. Saying so in the payload keeps the client from
                # inventing an onboarding flow that the server cannot honour.
                "self_service_signup": False,
                "rotatable_via_api": False,
                "how_to_rotate": (
                    "Set TMESIM_STUDIO_API_KEY on the host and restart the "
                    "process. There is no rotation endpoint: issuing keys "
                    "requires a key store this deployment does not have."
                ),
            },
            "deployment": {
                "single_shared_key": bool(configured),
                "data_roots_confined": bool(
                    os.environ.get("TMESIM_STUDIO_DATA_ROOTS", "").strip()),
                "limits": dict(limits),
                "started_at_unix": boot,
                "uptime_s": max(0.0, time.time() - boot),
            },
            "usage": {
                # Every number here is a COUNT OF THINGS THAT EXIST ON DISK
                # RIGHT NOW, not a meter. Deleting a run decrements it. It is
                # not a billing record and must never be presented as one.
                "kind": "current inventory, not a meter",
                "runs_by_status": by_status,
                "runs_total": sum(by_status.values()),
                "projects": len(store.list()),
                "reports": n_reports,
                "reconstruction_cache": (
                    cache.stats() if cache is not None else None),
                # The honest limitation, stated in the payload rather than
                # left for the UI to remember: one key, no identities.
                "attribution": (
                    "none — this deployment authenticates with a single "
                    "shared key, so no run can be attributed to a person"
                ),
            },
            "billing": _billing(),
        }

    return router


def _billing() -> Dict:
    """The billing block. Scaffolding, and labelled as scaffolding.

    Everything a billing integration needs a home for has a home here, and
    every one of those homes is empty and says it is empty. The alternative --
    omitting the block until a provider exists -- means the client has no shape
    to render and the shape gets invented at the last minute, in the UI, where
    it is one step from becoming a claim.

    ``enabled`` is driven by ``TMESIM_BILLING_PROVIDER`` so wiring a provider
    later is a config change plus a real implementation, and CANNOT be
    mistaken for already working: with a provider named but no implementation
    mounted, ``status`` reads ``configured-but-not-implemented``, which is a
    state an operator will notice.
    """
    provider = os.environ.get("TMESIM_BILLING_PROVIDER", "").strip() or None
    return {
        "enabled": False,
        "provider": provider,
        "status": ("configured-but-not-implemented" if provider
                   else "not-configured"),
        "plan": None,
        "seats": None,
        "period": None,
        "metered_quantities": [],
        "invoices": [],
        # A metering PRECONDITION, not a roadmap item. Usage cannot be billed
        # to anyone until requests carry an identity, and today they carry a
        # single shared secret. Anyone reading this block needs that fact
        # before they read anything else in it.
        "blocked_on": [
            "per-caller identity: the API authenticates with one shared key, "
            "so there is no principal to bill",
            "a durable usage record: usage counts are an inventory of files "
            "on disk and go down when a run is deleted",
        ],
        "note": (
            "Billing is not enabled on this deployment. No charge is made, no "
            "payment method is stored, and no usage is metered."
        ),
    }


# --------------------------------------------------------------------------- #
# Host memory, and whether this server's ceiling actually protects it
# --------------------------------------------------------------------------- #
#
# THE MEASUREMENT THIS IS BUILT ON, AND ITS LIMITS.
#
# app/server/api.py records two points, both taken on a 512 MB instance
# (Render's free tier), flagship Xenium GBM, in a real container with no bind
# mount:
#
#     46,070 voxels (dx = 30 um)  -> peaks at 220 MiB, comfortable
#    286,200 voxels (dx = 12 um)  -> the container is OOM-KILLED
#
# TWO POINTS AT ONE INSTANCE SIZE IS NOT A LAW. It cannot be turned into a
# bytes-per-voxel constant and it must not be: 220 MiB across 46,070 voxels is
# ~5 KiB a voxel, which is nothing like the 8 bytes a float64 field costs, so
# the peak is dominated by the cell population and the interpreter and NOT by
# the grid. Fitting a line through two points of a quantity that is mostly not
# the independent variable would produce a confident number with nothing under
# it — the exact failure this repository keeps finding.
#
# So nothing here derives a ceiling. This code detects what the host actually
# grants, and reports one comparison that is a plain fact rather than a model:
# does this server's configured ceiling still permit the spec that was MEASURED
# to kill an instance of this size? If it does, the operator is one ordinary
# request away from a dead container, and that is worth saying out loud at
# startup and on the Account panel instead of leaving it to be discovered.

#: cgroup v2, then v1. Absent on macOS and on a bare host, which is not an
#: error — it means "no container limit", and the answer is None.
_CGROUP_V2 = "/sys/fs/cgroup/memory.max"
_CGROUP_V1 = "/sys/fs/cgroup/memory/memory.limit_in_bytes"

#: cgroup v1 writes a sentinel rather than "max" when there is no limit. The
#: exact value varies with page size, so anything absurdly large is treated as
#: "unlimited" rather than compared against a magic constant.
_UNLIMITED_ABOVE = 1 << 50  # 1 PiB

#: The two recorded points. Values, not a model. See the block above.
MEASURED_INSTANCE_MB = 512
MEASURED_OK_VOXELS = 46_070
MEASURED_OK_PEAK_MIB = 220
MEASURED_OOM_VOXELS = 286_200


def _read_int(path: str) -> Optional[int]:
    try:
        with open(path) as fh:
            raw = fh.read().strip()
    except OSError:
        return None
    if raw == "max":
        return None
    try:
        return int(raw)
    except ValueError:
        return None


def memory_limit() -> Dict:
    """What this host actually grants the process, and how we know.

    ``TMESIM_HOST_MEMORY_MB`` overrides everything. It exists because a
    platform can enforce a limit it does not expose through cgroups, and
    because a test needs to drive this without a container -- an operator who
    knows the number can state it rather than have the code shrug.

    ``bytes`` is ``None`` when nothing declares a limit. That is reported as
    "not declared", never as "unlimited": a host can be capped by something
    this process cannot see, and the difference matters to anyone reading the
    warning below.
    """
    override = os.environ.get("TMESIM_HOST_MEMORY_MB", "").strip()
    if override:
        try:
            mb = float(override)
            if mb > 0:
                return {"bytes": int(mb * 1024 * 1024),
                        "source": "TMESIM_HOST_MEMORY_MB"}
        except ValueError:
            pass  # a malformed override falls through to detection
    for path, source in ((_CGROUP_V2, "cgroup v2 memory.max"),
                         (_CGROUP_V1, "cgroup v1 memory.limit_in_bytes")):
        v = _read_int(path)
        if v is not None and 0 < v < _UNLIMITED_ABOVE:
            return {"bytes": v, "source": source}
    return {"bytes": None, "source": "not declared"}


def oom_exposure(max_grid_cells: int, ceiling_is_default: bool) -> Dict:
    """Is this server's ceiling still above the spec measured to kill it?

    Returns the recorded measurements verbatim plus two booleans, both of which
    are comparisons rather than predictions:

    ``ceiling_permits_measured_oom``
        the configured ceiling is at or above the voxel count that was measured
        to OOM-kill a 512 MB instance. True says the refusal path will NOT stop
        that spec. It says nothing about whether THIS host would die of it.

    ``host_at_or_below_measured_instance``
        this host's declared memory is at or below the 512 MB the measurement
        was taken on. Unknown memory reads as None, not False -- an undetected
        limit is not the same as a generous one.

    When both are true the operator is one ordinary request from a dead
    container, and that is stated. When memory is undetected the warning is
    softened to name the uncertainty rather than dropped: the ceiling is still
    above a spec known to have killed something.
    """
    mem = memory_limit()
    mb = (mem["bytes"] / 1024 / 1024) if mem["bytes"] else None
    at_or_below = None if mb is None else mb <= MEASURED_INSTANCE_MB
    permits = int(max_grid_cells) >= MEASURED_OOM_VOXELS
    return {
        "host_memory_bytes": mem["bytes"],
        "host_memory_source": mem["source"],
        "ceiling_is_default": bool(ceiling_is_default),
        "ceiling_permits_measured_oom": permits,
        "host_at_or_below_measured_instance": at_or_below,
        "measured": {
            "instance_mb": MEASURED_INSTANCE_MB,
            "ok_voxels": MEASURED_OK_VOXELS,
            "ok_peak_mib": MEASURED_OK_PEAK_MIB,
            "oom_voxels": MEASURED_OOM_VOXELS,
            "note": ("Two points at one instance size, recorded in "
                     "app/server/api.py. NOT a bytes-per-voxel law: the "
                     "measured peak is dominated by the cell population, not "
                     "the grid, so no ceiling is derived from it."),
        },
    }


def oom_warning(max_grid_cells: int, ceiling_is_default: bool) -> Optional[str]:
    """The startup line, or ``None`` when there is nothing to warn about.

    Returns a STRING rather than logging, so ``__main__.py`` prints it in the
    banner where an operator is already looking and the tests can assert on it
    without capturing a logger.
    """
    ex = oom_exposure(max_grid_cells, ceiling_is_default)
    if not ex["ceiling_permits_measured_oom"]:
        return None  # the ceiling already refuses the spec that killed a box
    mb = ex["host_memory_bytes"]
    if ex["host_at_or_below_measured_instance"]:
        return (
            f"  !!  this host declares {mb / 1024 / 1024:,.0f} MB "
            f"({ex['host_memory_source']}), at or below the "
            f"{MEASURED_INSTANCE_MB} MB instance where {MEASURED_OOM_VOXELS:,} "
            f"voxels OOM-KILLED the container -- and the grid ceiling is "
            f"{max_grid_cells:,}, which still permits it.\n"
            "  !!  One ordinary request can take the whole service down. Set "
            "TMESIM_MAX_GRID_CELLS to what this machine can afford (100000 "
            "permits dx >= 20.25 um on the flagship section and refuses both "
            "dx = 12 and dx = 20)."
        )
    if mb is None and ceiling_is_default:
        return (
            f"  !   no memory limit is declared to this process, and the grid "
            f"ceiling is the {max_grid_cells:,} default -- which permits the "
            f"{MEASURED_OOM_VOXELS:,}-voxel spec that OOM-killed a "
            f"{MEASURED_INSTANCE_MB} MB container. If this host is small, set "
            "TMESIM_MAX_GRID_CELLS."
        )
    return None
