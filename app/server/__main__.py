"""``./.venv/bin/python -m app.server`` — start the tmesim Studio server.

Binds to 127.0.0.1 by default: this is a local diagnostic workstation
application, not a network service, and the data it exposes is a patient tissue
reconstruction. Pass ``--host`` explicitly if you really mean to expose it.
"""
from __future__ import annotations

import argparse
import os
import sys


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
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    root = _ensure_repo_on_path()

    if args.cache_size is not None:
        os.environ["TMESIM_RECON_CACHE_SIZE"] = str(args.cache_size)

    import uvicorn

    from .cache import RECON_CACHE, runs_dir

    print(f"tmesim Studio backend")
    print(f"  repo        {root}")
    print(f"  run history {runs_dir()}")
    print(f"  recon cache maxsize={RECON_CACHE.maxsize}")
    print(f"  serving     http://{args.host}:{args.port}")

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
