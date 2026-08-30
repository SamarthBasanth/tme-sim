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
