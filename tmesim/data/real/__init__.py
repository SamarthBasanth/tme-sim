"""Loaders for REAL, open glioblastoma datasets cached under ``data_cache/``.

Each module in this package exposes a small, dependency-light loader that reads
a genuinely-downloaded open dataset (no synthetic stand-ins) and returns plain
numpy/pandas objects. See ``docs/real_data/`` for the provenance of each cache.
"""
from __future__ import annotations

from .scref_gbm import load_reference
from .singlecell_spatial import load as load_singlecell_spatial
from .visium import list_samples, load_sample, sample_info
from .ivygap import load_structures, structure_signatures
from .drug_response import (
    SOURCES as DRUG_RESPONSE_SOURCES,
    load_gdsc,
    load_hgcc,
    load_hgcc_clinical,
    summarize as summarize_drug_response,
)

__all__ = [
    "load_reference",
    "load_singlecell_spatial",
    "list_samples",
    "load_sample",
    "sample_info",
    "load_structures",
    "structure_signatures",
    "load_hgcc",
    "load_hgcc_clinical",
    "load_gdsc",
    "DRUG_RESPONSE_SOURCES",
    "summarize_drug_response",
]
