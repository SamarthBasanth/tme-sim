"""The typed TPL syntax tree. **Every node carries a line and a column.**

That is not decoration: the editor renders diagnostics inline against these
positions, and the executor attributes a refusal ("this sweep is 96 runs",
"this panel cannot support state_mode=markers") to the exact token that caused
it. A node without a position cannot be reported, so :class:`Node` gives every
subclass ``line`` / ``col`` / ``end_col`` and :meth:`Node.to_dict` always emits
them.

Positions are 1-based; ``end_col`` is exclusive, matching
:class:`tmesim.protocol.errors.Diagnostic`.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

from .units import Quantity

__all__ = [
    "Node", "Program", "DrugDef", "DrugField", "Ic50FromGdsc", "ProtocolDef",
    "Route", "Dose", "Every", "Repeat", "Rest", "Cycles", "Sweep", "RunStmt",
    "CompareStmt", "ExpectStmt", "AllowRuns", "RunOptions",
    "DRUG_SOURCES", "ROUTES", "OVERRIDABLE_FIELDS", "SWEEPABLE", "METRICS",
    "METRIC_ALIASES", "STATE_MODES", "SENSITIVITY_MODES",
]

#: ``drug X from <source> "<value>"`` — the three sources, in the spec's order
#: of decreasing evidence.
DRUG_SOURCES = ("smiles", "library", "file")

#: ``route`` values the dosing model can honour.
ROUTES = ("oral", "iv", "ip")

#: ``DrugParams`` fields a drug block may override, and whether the value is a
#: PHYSICAL quantity (unit mandatory) or a dimensionless number (unit refused).
#: The unit column is deliberately restricted to the frozen unit table: the
#: engine's ``binding`` / ``clearance`` are 1/s and ``molecular_weight`` is
#: g/mol, none of which the spec's unit list contains, so those three take a
#: bare number in the engine's own unit and the executor records that.
OVERRIDABLE_FIELDS: Dict[str, Dict] = {
    "ic50": {"dimension": "concentration", "engine_unit": "uM",
             "note": "half-maximal effect concentration"},
    "hill": {"dimension": None, "engine_unit": "",
             "note": "Hill coefficient (dimensionless)"},
    "diffusivity": {"dimension": "diffusivity", "engine_unit": "um2/s",
                    "note": "effective tissue diffusion coefficient"},
    "binding": {"dimension": None, "engine_unit": "1/s",
                "note": "first-order tissue binding/uptake rate, in 1/s "
                        "(1/s is not in the TPL unit table, so this field "
                        "takes a bare number in the engine's own unit)"},
    "clearance": {"dimension": None, "engine_unit": "1/s",
                  "note": "first-order tissue clearance rate, in 1/s (bare "
                          "number; see `binding`)"},
    "bbb_penetrance": {"dimension": "fraction", "engine_unit": "",
                       "note": "blood-brain-barrier penetrance in [0, 1]; a "
                               "bare number or a percentage"},
    "molecular_weight": {"dimension": None, "engine_unit": "g/mol",
                         "note": "molecular weight in g/mol (bare number; "
                                 "g/mol is not in the TPL unit table)"},
}

#: Parameters a ``sweep`` may vary, with the dimension each requires.
SWEEPABLE: Dict[str, Dict] = {
    "dose": {"dimension": "dose", "target": "protocol.dose"},
    "over": {"dimension": "time", "target": "protocol.over"},
    "every": {"dimension": "time", "target": "protocol.every"},
    "repeat": {"dimension": None, "target": "protocol.repeat", "integer": True},
    "cycles": {"dimension": None, "target": "protocol.cycles", "integer": True},
    "dx": {"dimension": "length", "target": "run.dx_um"},
}

#: ``report`` / ``expect`` / ``rank by`` metric names -> the RunResult key.
METRICS: Dict[str, str] = {
    "kill": "kill_fraction",
    "coverage": "therapeutic_coverage_1x",
    "mean_auc": "mean_tumor_auc",
    "penetration_length": "penetration_length_um",
}
#: Accepted spellings that map onto :data:`METRICS`.
METRIC_ALIASES: Dict[str, str] = {
    "kill_fraction": "kill",
    "therapeutic_coverage_1x": "coverage",
    "coverage_1x": "coverage",
    "mean_tumor_auc": "mean_auc",
    "auc": "mean_auc",
    "penetration_length_um": "penetration_length",
    "penetration": "penetration_length",
}

#: Mirrors ``app.server.cache.STATE_MODES`` / ``SENSITIVITY_MODES``; asserted
#: equal in ``tests/test_protocol.py`` so a divergence is caught, not shipped.
STATE_MODES = ("oxygen", "vessel_distance", "markers")
SENSITIVITY_MODES = ("uniform", "determinant")


# --------------------------------------------------------------------------- #
# Base
# --------------------------------------------------------------------------- #
@dataclass
class Node:
    """Every AST node is positioned. No exceptions."""

    line: int = 0
    col: int = 0
    end_col: int = 0

    def pos(self) -> Dict[str, int]:
        return {"line": int(self.line), "col": int(self.col),
                "end_col": int(self.end_col or self.col + 1)}

    def to_dict(self) -> Dict:
        return {"kind": type(self).__name__, **self.pos()}


def _q(q: Optional[Quantity]) -> Optional[Dict]:
    if q is None:
        return None
    return {"value": q.value, "unit": q.unit, "line": q.line, "col": q.col,
            "end_col": q.end_col}


# --------------------------------------------------------------------------- #
# Drug definition
# --------------------------------------------------------------------------- #
@dataclass
class Ic50FromGdsc(Node):
    """``ic50 from gdsc`` — a MEASURED potency, or a refusal.

    The executor calls :func:`tmesim.molecular.potency.resolve_ic50` and ERRORS
    if the resolved confidence is not ``'measured'``. It never silently accepts
    the literature or flagged-default tiers that function falls back to: asking
    for a measurement and being handed a placeholder is exactly the failure the
    spec forbids.
    """

    def to_dict(self) -> Dict:
        return {"kind": "Ic50FromGdsc", **self.pos()}


@dataclass
class DrugField(Node):
    """One ``<field> <value> [<unit>]`` override inside a ``drug`` block."""

    name: str = ""
    value: Optional[Quantity] = None          # None when `source` is set
    source: Optional[Ic50FromGdsc] = None     # `ic50 from gdsc`

    def to_dict(self) -> Dict:
        return {"kind": "DrugField", "name": self.name,
                "value": _q(self.value),
                "source": self.source.to_dict() if self.source else None,
                **self.pos()}


@dataclass
class DrugDef(Node):
    """``drug IDENT from <source> "<value>" ... end``."""

    name: str = ""
    source_kind: str = ""     # one of DRUG_SOURCES
    source_value: str = ""
    source_line: int = 0
    source_col: int = 0
    fields: List[DrugField] = field(default_factory=list)

    def to_dict(self) -> Dict:
        return {"kind": "DrugDef", "name": self.name,
                "source_kind": self.source_kind,
                "source_value": self.source_value,
                "fields": [f.to_dict() for f in self.fields], **self.pos()}


# --------------------------------------------------------------------------- #
# Protocol
# --------------------------------------------------------------------------- #
@dataclass
class Route(Node):
    value: str = "iv"

    def to_dict(self) -> Dict:
        return {"kind": "Route", "value": self.value, **self.pos()}


@dataclass
class Dose(Node):
    """``dose N <dose-unit> [over T <time-unit>]``.

    ``over`` is an INFUSION DURATION, not a bolus. The executor honours the
    distinction — see :mod:`tmesim.protocol.dosing`.
    """

    amount: Optional[Quantity] = None
    over: Optional[Quantity] = None

    def to_dict(self) -> Dict:
        return {"kind": "Dose", "amount": _q(self.amount),
                "over": _q(self.over), **self.pos()}


@dataclass
class Every(Node):
    interval: Optional[Quantity] = None

    def to_dict(self) -> Dict:
        return {"kind": "Every", "interval": _q(self.interval), **self.pos()}


@dataclass
class Repeat(Node):
    count: int = 1

    def to_dict(self) -> Dict:
        return {"kind": "Repeat", "count": self.count, **self.pos()}


@dataclass
class Rest(Node):
    duration: Optional[Quantity] = None

    def to_dict(self) -> Dict:
        return {"kind": "Rest", "duration": _q(self.duration), **self.pos()}


@dataclass
class Cycles(Node):
    count: int = 1

    def to_dict(self) -> Dict:
        return {"kind": "Cycles", "count": self.count, **self.pos()}


@dataclass
class ProtocolDef(Node):
    """``protocol IDENT ... end``."""

    name: str = ""
    route: Optional[Route] = None
    dose: Optional[Dose] = None
    every: Optional[Every] = None
    repeat: Optional[Repeat] = None
    rest: Optional[Rest] = None
    cycles: Optional[Cycles] = None

    def statements(self) -> List[Node]:
        return [s for s in (self.route, self.dose, self.every, self.repeat,
                            self.rest, self.cycles) if s is not None]

    def to_dict(self) -> Dict:
        return {"kind": "ProtocolDef", "name": self.name,
                "route": self.route.to_dict() if self.route else None,
                "dose": self.dose.to_dict() if self.dose else None,
                "every": self.every.to_dict() if self.every else None,
                "repeat": self.repeat.to_dict() if self.repeat else None,
                "rest": self.rest.to_dict() if self.rest else None,
                "cycles": self.cycles.to_dict() if self.cycles else None,
                **self.pos()}


# --------------------------------------------------------------------------- #
# Sweep
# --------------------------------------------------------------------------- #
@dataclass
class Sweep(Node):
    """``sweep <param> from A to B <unit> (step S | in N steps)`` or
    ``sweep <param> V <unit>`` (a fixed override, one point)."""

    param: str = ""
    kind: str = "range"          # "range" | "fixed"
    start: Optional[Quantity] = None
    stop: Optional[Quantity] = None
    step: Optional[float] = None
    n_steps: Optional[int] = None
    fixed: Optional[Quantity] = None
    param_col: int = 0

    def to_dict(self) -> Dict:
        return {"kind": "Sweep", "param": self.param, "sweep_kind": self.kind,
                "start": _q(self.start), "stop": _q(self.stop),
                "step": self.step, "n_steps": self.n_steps,
                "fixed": _q(self.fixed), **self.pos()}


# --------------------------------------------------------------------------- #
# Run / compare / expect / allow
# --------------------------------------------------------------------------- #
@dataclass
class RunOptions(Node):
    """The ``run`` / ``compare`` body: ``dx``, ``state_mode``, ``sensitivity``,
    ``report``, and (``compare`` only) ``rank by``."""

    dx: Optional[Quantity] = None
    state_mode: Optional[str] = None
    state_mode_line: int = 0
    state_mode_col: int = 0
    sensitivity: Optional[str] = None
    sensitivity_line: int = 0
    sensitivity_col: int = 0
    report: List[str] = field(default_factory=list)
    report_line: int = 0
    rank_by: Optional[str] = None
    rank_by_line: int = 0
    rank_by_col: int = 0

    def to_dict(self) -> Dict:
        return {"kind": "RunOptions", "dx": _q(self.dx),
                "state_mode": self.state_mode,
                "sensitivity": self.sensitivity, "report": list(self.report),
                "rank_by": self.rank_by, **self.pos()}


@dataclass
class RunStmt(Node):
    """``run DRUG with PROTOCOL on SECTION ... end``."""

    drug: str = ""
    drug_col: int = 0
    protocol: str = ""
    protocol_col: int = 0
    section: str = ""
    section_col: int = 0
    options: RunOptions = field(default_factory=RunOptions)

    @property
    def drugs(self) -> List[str]:
        return [self.drug]

    def to_dict(self) -> Dict:
        return {"kind": "RunStmt", "drug": self.drug,
                "protocol": self.protocol, "section": self.section,
                "options": self.options.to_dict(), **self.pos()}


@dataclass
class CompareStmt(Node):
    """``compare A, B with PROTOCOL on SECTION ... end`` — an N-drug run."""

    drugs: List[str] = field(default_factory=list)
    drug_cols: List[int] = field(default_factory=list)
    protocol: str = ""
    protocol_col: int = 0
    section: str = ""
    section_col: int = 0
    options: RunOptions = field(default_factory=RunOptions)

    def to_dict(self) -> Dict:
        return {"kind": "CompareStmt", "drugs": list(self.drugs),
                "protocol": self.protocol, "section": self.section,
                "options": self.options.to_dict(), **self.pos()}


@dataclass
class ExpectStmt(Node):
    """``expect METRIC <cmp> N [UNIT]``.

    Reported as pass/fail, never enforced: a failed expectation marks the run,
    it does not error. The point is to state a hypothesis before running it.
    """

    metric: str = ""
    metric_col: int = 0
    comparator: str = ">"
    value: Optional[Quantity] = None

    def to_dict(self) -> Dict:
        return {"kind": "ExpectStmt", "metric": self.metric,
                "comparator": self.comparator, "value": _q(self.value),
                **self.pos()}


@dataclass
class AllowRuns(Node):
    """``allow runs N`` — the explicit consent that lifts the sweep budget."""

    count: int = 0

    def to_dict(self) -> Dict:
        return {"kind": "AllowRuns", "count": self.count, **self.pos()}


@dataclass
class Program(Node):
    """A whole TPL source file."""

    drugs: List[DrugDef] = field(default_factory=list)
    protocols: List[ProtocolDef] = field(default_factory=list)
    sweeps: List[Sweep] = field(default_factory=list)
    runs: List[Node] = field(default_factory=list)      # RunStmt | CompareStmt
    expects: List[ExpectStmt] = field(default_factory=list)
    allow_runs: Optional[AllowRuns] = None
    source: str = ""

    def to_dict(self) -> Dict:
        return {
            "kind": "Program",
            "drugs": [d.to_dict() for d in self.drugs],
            "protocols": [p.to_dict() for p in self.protocols],
            "sweeps": [s.to_dict() for s in self.sweeps],
            "runs": [r.to_dict() for r in self.runs],
            "expects": [e.to_dict() for e in self.expects],
            "allow_runs": self.allow_runs.to_dict() if self.allow_runs else None,
            **self.pos(),
        }
