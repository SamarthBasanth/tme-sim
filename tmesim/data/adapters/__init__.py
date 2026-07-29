"""Format adapters: turn a foreign on-disk deposit into a layout
``tmesim/data/ingest.py`` already recognises, rather than teaching ``ingest.py``
a new detector for every new source.

See :mod:`tmesim.data.adapters.gbmspace` for the first (and so far only)
adapter here, and ``results/gbmspace_ingestion.md`` / ``results/gbmspace_adapter.md``
for why an adapter was the right shape for that source specifically.
"""
from __future__ import annotations
