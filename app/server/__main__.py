"""``./.venv/bin/python -m app.server`` — start the tmesim Studio server.

Binds to 127.0.0.1 by default: this is a local research workstation application,
not a network service, and the data it exposes is a reconstruction of a
de-identified tissue section. It performs no diagnosis and is not built for
identifiable patient data. Pass ``--host`` explicitly if you really mean to
expose it -- doing so with no ``--api-key`` and no ``TMESIM_STUDIO_API_KEY`` set
auto-generates one and prints it once, rather than serving unauthenticated on a
non-loopback address. See docs/DEPLOYMENT.md (§0 for the de-identified-tissue
scope).

That auto-generation is keyed on the BIND ADDRESS, which is not the same thing
as who can reach the process: a loopback bind fronted by an ngrok/Cloudflare
tunnel or a reverse proxy is world-reachable while this process still sees
``127.0.0.1`` and generates nothing. The default is deliberately left alone --
unauthenticated loopback is the documented single-user workflow -- but every
start with no key now prints ``_no_auth_warning()`` so the operator cannot
reach that situation without having been told.
"""
from __future__ import annotations

import argparse
import os
import secrets
import sys

_LOOPBACK_HOSTS = {"127.0.0.1", "localhost", "::1"}


def _ensure_repo_on_path() -> str:
    """Put the repository root on ``sys.path`` so ``import tmesim`` works from
    any working directory."""
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    if root not in sys.path:
        sys.path.insert(0, root)
    return root


def _no_auth_warning(host: str) -> str:
    """The startup banner block printed whenever no API key is in force.

    SECURITY: the auto-generate rule above reasons about the bind address, and
    the bind address is not reachability. Tunnelling a loopback service
    (``ngrok http 8000``, ``cloudflared``, an SSH ``-R`` forward, a reverse
    proxy with no auth of its own) is the ordinary way people demo software,
    and it is exactly the case where the two differ: ``is_loopback`` is true,
    no key is generated, and the server is on the internet unauthenticated.
    Nothing in the process can detect that, so the only honest fix is to say
    it out loud on every unauthenticated start rather than to guess.

    Kept as a function so the wording is asserted once, in
    tests/test_app_server.py, instead of being duplicated per call site.
    """
    return "\n".join((
        "",
        "  !!  WARNING: NO API KEY -- every /api/ route is UNAUTHENTICATED.",
        "  !!",
        f"  !!  Binding {host} is NOT protection. This process knows only what",
        "  !!  it bound to, not who can reach it. Put an ngrok/Cloudflare",
        "  !!  tunnel, an SSH -R forward or a reverse proxy in front of a",
        "  !!  loopback bind and it is reachable from the internet -- while",
        "  !!  this process still sees 127.0.0.1 and auto-generates no key,",
        "  !!  because that rule reads the bind address, not reachability.",
        "  !!  If anything at all forwards traffic here, you are exposed.",
        "  !!",
        "  !!  Unauthenticated, POST /api/sections/scan takes an arbitrary",
        "  !!  filesystem path and reports on any directory this process can",
        "  !!  read -- that is filesystem reconnaissance, not just data access.",
        "  !!",
        "  !!  Fix: restart with TMESIM_STUDIO_API_KEY=<a-real-secret> set, or",
        "  !!  pass --api-key <a-real-secret>. See docs/DEPLOYMENT.md §1e.",
        "  !!  A key authenticates /api/sections/scan; it does not confine the",
        "  !!  path it will read. Only give it to people you'd give a shell to.",
        "",
    ))


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="python -m app.server",
        description="Codon Studio backend (API_CONTRACT.md v1).",
    )
    p.add_argument("--host", default="127.0.0.1",
                   help="bind address (default 127.0.0.1, local only)")
    p.add_argument("--port", type=int, default=8000, help="bind port")
    p.add_argument("--reload", action="store_true",
                   help="uvicorn autoreload (development)")
    p.add_argument("--log-level", default="info",
                   choices=("critical", "error", "warning", "info", "debug"))
    p.add_argument("--cache-size", type=int, default=None,
                   help="max cached reconstructions (default 4; each holds a "
                        "full cell population and several float64 grids)")
    p.add_argument("--api-key", default=None,
                   help="require this value in every request's X-API-Key "
                        "header on /api/* routes. If omitted and --host is "
                        "non-loopback, one is generated and printed once. "
                        "Pass --api-key '' to explicitly disable (not "
                        "recommended off loopback).")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    root = _ensure_repo_on_path()

    if args.cache_size is not None:
        os.environ["TMESIM_RECON_CACHE_SIZE"] = str(args.cache_size)

    is_loopback = args.host in _LOOPBACK_HOSTS
    api_key = args.api_key
    generated = False
    if api_key is None:
        # --api-key wasn't passed. Respect an already-configured environment
        # (Docker ENV, systemd, a shell export) rather than silently clobber
        # it -- only decide a value ourselves when the environment has none.
        api_key = os.environ.get("TMESIM_STUDIO_API_KEY", "")
        if not api_key and not is_loopback:
            api_key = secrets.token_urlsafe(24)
            generated = True
    os.environ["TMESIM_STUDIO_API_KEY"] = api_key

    # SECURITY: on the default loopback bind the only legitimate Host values
    # are the loopback names, so pin them -- that blocks DNS rebinding, where
    # a hostile domain resolves to 127.0.0.1 and thereby becomes same-origin
    # with this server (CORS cannot help there). Left unset for a non-loopback
    # bind, where we do not know the operator's hostname and refusing it would
    # break the deployment; set TMESIM_STUDIO_ALLOWED_HOSTS yourself there.
    #
    # IPv6 IS EXCLUDED, and not as an oversight. Starlette's TrustedHostMiddleware
    # matches on ``headers["host"].split(":")[0]`` to strip the port. That is
    # correct for a name or an IPv4 literal and wrong for an IPv6 one: a browser
    # sends ``Host: [::1]:8000``, which splits to ``[``. No allowlist entry can
    # ever equal ``[``, so pinning here would make ``--host ::1`` reject 100 % of
    # its own traffic with "400 Invalid host header" -- a supported flag broken
    # by its own protection. Enabling the middleware on a bind it cannot express
    # buys no safety and costs the whole server, so IPv6 loopback keeps the
    # non-loopback treatment: unset, and the operator may set it explicitly.
    _ipv6_bind = ":" in args.host
    if is_loopback and not _ipv6_bind and "TMESIM_STUDIO_ALLOWED_HOSTS" not in os.environ:
        os.environ["TMESIM_STUDIO_ALLOWED_HOSTS"] = "localhost,127.0.0.1"

    import uvicorn

    from .cache import RECON_CACHE, runs_dir

    print(f"Codon Studio backend")
    print(f"  repo        {root}")
    print(f"  run history {runs_dir()}")
    print(f"  recon cache maxsize={RECON_CACHE.maxsize}")
    print(f"  serving     http://{args.host}:{args.port}")
    # Bring-your-own-tissue is the whole value proposition and it was
    # reachable only by hand-POSTing JSON (results/first_user_readiness.md).
    # Printing the page here is how a first user finds out it exists.
    print(f"  your tissue http://{args.host}:{args.port}/ingest")
    if api_key:
        # Name the source that actually won. There are three, and this line is
        # exactly where an operator looks to work out which: it used to print
        # "set via --api-key" for a key that came from the environment, which
        # is the one situation where the answer is not obvious (an exported
        # TMESIM_STUDIO_API_KEY and an --api-key flag can disagree, and the
        # flag wins) and the one place the banner could have said so.
        if generated:
            label = "generated (not persisted -- save it now)"
        elif args.api_key is not None:
            label = "set via --api-key"
        else:
            label = "set via TMESIM_STUDIO_API_KEY"
        print(f"  api key     {api_key}  [{label}]")
        print(f"              every request needs header: X-API-Key: {api_key}")
    else:
        # No key, for ANY reason -- loopback default, `--api-key ''`, or an
        # empty TMESIM_STUDIO_API_KEY in the environment. The loopback case
        # used to print nothing at all, which is the defect: the risk was
        # invisible precisely where it is easiest to hit (see
        # _no_auth_warning). The default is unchanged -- still no key, still
        # serves -- only now it says so.
        print("  api key     NONE -- no authentication on any /api/ route")
        if not is_loopback:
            print(f"              bound to {args.host}, which is not a "
                  "loopback address")
        print(_no_auth_warning(args.host), end="")

    if args.reload:
        uvicorn.run("app.server.api:app", host=args.host, port=args.port,
                    reload=True, log_level=args.log_level)
        return 0

    from .api import create_app

    uvicorn.run(create_app(), host=args.host, port=args.port,
                log_level=args.log_level)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
