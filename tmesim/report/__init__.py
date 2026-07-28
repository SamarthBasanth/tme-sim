"""Reporting primitives shared by the CLI experiments and the Studio backend.

Nothing in here simulates anything. It is the *explanation* layer: given a set
of drug rows that were already produced by the production engine, decompose each
drug's standing into exactly additive log10 terms and check that decomposition
against the simulator's own ranking.

It lives in ``tmesim/`` rather than in ``experiments/`` or ``app/`` so that
``experiments/patient_report.py`` (the CLI artifact) and ``app/server/report.py``
(the in-product artifact) run the SAME arithmetic. Two implementations of an
attribution would be two sciences.
"""
from .attribution import (
    TERM_LABEL,
    TERMS,
    attribute,
    gap_reason,
    pct_rank,
    spearman,
)

__all__ = [
    "TERMS",
    "TERM_LABEL",
    "attribute",
    "gap_reason",
    "pct_rank",
    "spearman",
]
