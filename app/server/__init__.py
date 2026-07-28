"""tmesim Studio backend — the local application server.

Serves the frozen ``API_CONTRACT.md`` v1 interface over FastAPI and drives the
REAL simulation engine underneath it. Run it with::

    ./.venv/bin/python -m app.server            # http://127.0.0.1:8000

Layout:

``payloads.py``
    TMSG / TMSF binary wire formats (little-endian, typed-array-ready).
``cache.py``
    Section registry (what is actually loadable from ``data_cache/``) and the
    LRU cache of reconstructions keyed by
    ``(section, dx_um, vascular_mode, oxygen_mode)``.
``jobs.py``
    Run execution against the production engine entry points, plus the queue,
    the worker thread, and on-disk run history in ``app/runs/``.
``api.py``
    The routes, the drug/defaults registries, and ``POST /api/validate``.

This package WRAPS ``tmesim``. It contains no science: every number it reports
comes from ``reconstruct.real_pipeline``, ``sim.engine``,
``validation.calibration``, ``validation.delivery_metrics``,
``molecular.potency`` or ``molecular.binding_site_barrier``, unchanged.
"""
from __future__ import annotations

__all__ = ["create_app", "main"]

__version__ = "1.0.0"


def create_app(*args, **kwargs):
    """Re-export of :func:`app.server.api.create_app` (imported lazily)."""
    from .api import create_app as _create_app

    return _create_app(*args, **kwargs)


def main(argv=None) -> int:
    """Re-export of :func:`app.server.__main__.main` (imported lazily)."""
    from .__main__ import main as _main

    return _main(argv)
