"""``./.venv/bin/python -m app.server`` — start the tmesim Studio server.

Binds to 127.0.0.1 by default: this is a local diagnostic workstation
application, not a network service, and the data it exposes is a patient tissue
reconstruction. Pass ``--host`` explicitly if you really mean to expose it --
doing so with no ``--api-key`` and no ``TMESIM_STUDIO_API_KEY`` set auto-
generates one and prints it once, rather than serving unauthenticated on a
non-loopback address. See docs/DEPLOYMENT.md.
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


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="python -m app.server",
        description="tmesim Studio backend (API_CONTRACT.md v1).",
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

    import uvicorn

    from .cache import RECON_CACHE, runs_dir

    print(f"tmesim Studio backend")
    print(f"  repo        {root}")
    print(f"  run history {runs_dir()}")
    print(f"  recon cache maxsize={RECON_CACHE.maxsize}")
    print(f"  serving     http://{args.host}:{args.port}")
    if api_key:
        label = "generated (not persisted -- save it now)" if generated else "set via --api-key"
        print(f"  api key     {api_key}  [{label}]")
        print(f"              every request needs header: X-API-Key: {api_key}")
    elif not is_loopback:
        print("  api key     NONE -- serving on a non-loopback address with "
              "no authentication. Only do this on a trusted network.")

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
