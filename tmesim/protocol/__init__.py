"""TPL — the Tumour Protocol Language: lexer, parser, planner and executor.

``PROTOCOL_LANGUAGE.md`` is the frozen spec; this package implements it, and
``app/web/js/ui/protocol.js`` implements the editor against the same document.

Typical use::

    from tmesim.protocol import parse, plan, validate_source, execute_plan

    verdict = validate_source(src)          # {ok, diagnostics[] with line/col}
    p = plan(src)                           # expanded run list, cost, refusals
    result = execute_plan(p)                # metrics + provenance per run

Layout
------
``lexer``      tokens with line/col
``parser``     recursive descent -> ``ast``, with error recovery
``ast``        typed nodes, every one positioned
``units``      the frozen unit table, dimensions, explicit conversions
``dosing``     dose -> ``core.PKConfig``, with every assumption written out
``resolve``    drug sources (smiles / library / file) + per-field provenance,
               and the specimen capability check
``executor``   sweep expansion, the run budget, and the wrapped engine call
``errors``     positioned diagnostics; every refusal has a line number
"""
from __future__ import annotations

from .ast import METRICS, OVERRIDABLE_FIELDS, Program, SWEEPABLE
from .dosing import BSA_M2, DoseConversion, convert_dose, schedule_times
from .errors import (
    BudgetError,
    CapabilityError,
    Diagnostic,
    LexError,
    ParseError,
    ResolutionError,
    TplError,
    UnitError,
)
from .executor import DEFAULT_MAX_RUNS, Plan, RunPlan, execute_plan, plan, validate_source
from .lexer import Token, tokenize
from .parser import parse
from .resolve import ResolvedDrug, check_capability, load_drug_file, resolve_drug
from .units import UNITS, Quantity, convert

__all__ = [
    # language
    "tokenize", "Token", "parse", "Program", "METRICS", "SWEEPABLE",
    "OVERRIDABLE_FIELDS",
    # units
    "UNITS", "Quantity", "convert",
    # dosing
    "BSA_M2", "DoseConversion", "convert_dose", "schedule_times",
    # drugs + capability
    "ResolvedDrug", "resolve_drug", "load_drug_file", "check_capability",
    # planning / execution
    "plan", "Plan", "RunPlan", "execute_plan", "validate_source",
    "DEFAULT_MAX_RUNS",
    # errors
    "Diagnostic", "TplError", "LexError", "ParseError", "UnitError",
    "ResolutionError", "CapabilityError", "BudgetError",
]
