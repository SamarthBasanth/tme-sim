"""Diagnostics for TPL — every refusal carries a position and a reason.

The editor (``app/web/js/ui/protocol.js``) renders these inline, so the shape is
part of the contract:

.. code-block:: jsonc

    { "severity": "error",      // "error" | "warning"
      "line": 7,                 // 1-based
      "col": 8,                  // 1-based, of the offending token
      "end_col": 11,             // 1-based, exclusive
      "code": "missing-unit",
      "message": "dose: a unit is required ...",
      "hint": "write e.g. `dose 150 mg/m2`" }

Two rules this module exists to enforce:

* **No bare stack traces.** Every failure a user can cause is raised as a
  :class:`TplError` carrying at least one :class:`Diagnostic`. The routes turn
  those into a 200 ``{ok: false, diagnostics: [...]}`` or a 400 with the same
  body — never a 500 with a traceback.
* **Refusal, not default.** The spec's non-negotiable list (unknown drug,
  unsupported mode for the specimen, missing unit, an over-budget sweep, an
  ``ic50 from gdsc`` for a drug absent from the screen) all land here with a
  ``code`` the UI can key on.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

__all__ = [
    "Diagnostic",
    "TplError",
    "LexError",
    "ParseError",
    "UnitError",
    "ResolutionError",
    "CapabilityError",
    "BudgetError",
    "collect",
]

#: Severity levels the editor knows how to render.
SEVERITIES = ("error", "warning")


@dataclass
class Diagnostic:
    """One positioned, actionable message."""

    message: str
    line: int = 0
    col: int = 0
    end_col: int = 0
    severity: str = "error"
    code: str = "error"
    hint: Optional[str] = None

    def __post_init__(self) -> None:
        self.line = int(self.line or 0)
        self.col = int(self.col or 0)
        self.end_col = int(self.end_col or 0) or (self.col + 1 if self.col else 0)
        if self.severity not in SEVERITIES:
            self.severity = "error"

    # -- views -------------------------------------------------------------- #
    def to_dict(self) -> Dict:
        return {
            "severity": self.severity,
            "line": self.line,
            "col": self.col,
            "end_col": self.end_col,
            "code": self.code,
            "message": self.message,
            "hint": self.hint,
        }

    def __str__(self) -> str:  # pragma: no cover - convenience only
        where = f"line {self.line}, col {self.col}" if self.line else "input"
        out = f"{self.severity} [{self.code}] at {where}: {self.message}"
        return out + (f"  ({self.hint})" if self.hint else "")


class TplError(Exception):
    """Base class for every TPL failure. Always carries >= 1 diagnostic."""

    #: Default code used when a subclass does not supply one.
    code = "error"

    def __init__(self, message: str, line: int = 0, col: int = 0,
                 end_col: int = 0, hint: Optional[str] = None,
                 code: Optional[str] = None,
                 diagnostics: Optional[List[Diagnostic]] = None):
        super().__init__(message)
        if diagnostics:
            self.diagnostics: List[Diagnostic] = list(diagnostics)
        else:
            self.diagnostics = [Diagnostic(
                message=message, line=line, col=col, end_col=end_col,
                severity="error", code=code or self.code, hint=hint,
            )]

    def to_dicts(self) -> List[Dict]:
        return [d.to_dict() for d in self.diagnostics]

    @property
    def line(self) -> int:
        return self.diagnostics[0].line

    @property
    def col(self) -> int:
        return self.diagnostics[0].col


class LexError(TplError):
    """A character the tokenizer cannot begin a token with."""

    code = "lex"


class ParseError(TplError):
    """A token sequence the grammar does not accept."""

    code = "parse"


class UnitError(TplError):
    """A missing, unknown, or dimensionally-wrong unit."""

    code = "unit"


class ResolutionError(TplError):
    """A name that does not resolve: unknown drug, protocol, section, metric,
    or an ``ic50 from gdsc`` for a drug absent from the screen."""

    code = "resolution"


class CapabilityError(TplError):
    """A mode the chosen specimen's panel cannot support."""

    code = "capability"


class BudgetError(TplError):
    """A sweep larger than ``max_runs`` without an explicit ``allow runs N``."""

    code = "budget"


@dataclass
class collect:
    """A diagnostic accumulator, so a whole file's errors surface at once.

    The editor wants every problem in one round trip, not the first one. Passes
    that can continue after a failure append here; passes that cannot raise.
    """

    items: List[Diagnostic] = field(default_factory=list)

    def error(self, message: str, node=None, code: str = "error",
              hint: Optional[str] = None, line: int = 0, col: int = 0,
              end_col: int = 0) -> None:
        self._add("error", message, node, code, hint, line, col, end_col)

    def warn(self, message: str, node=None, code: str = "warning",
             hint: Optional[str] = None, line: int = 0, col: int = 0,
             end_col: int = 0) -> None:
        self._add("warning", message, node, code, hint, line, col, end_col)

    def extend(self, other: "collect") -> None:
        self.items.extend(other.items)

    def add(self, diag: Diagnostic) -> None:
        self.items.append(diag)

    def _add(self, severity: str, message: str, node, code: str,
             hint: Optional[str], line: int, col: int, end_col: int) -> None:
        if node is not None:
            line = line or int(getattr(node, "line", 0) or 0)
            col = col or int(getattr(node, "col", 0) or 0)
            end_col = end_col or int(getattr(node, "end_col", 0) or 0)
        self.items.append(Diagnostic(
            message=message, line=line, col=col, end_col=end_col,
            severity=severity, code=code, hint=hint,
        ))

    # -- views -------------------------------------------------------------- #
    @property
    def errors(self) -> List[Diagnostic]:
        return [d for d in self.items if d.severity == "error"]

    @property
    def warnings(self) -> List[Diagnostic]:
        return [d for d in self.items if d.severity == "warning"]

    @property
    def ok(self) -> bool:
        return not self.errors

    def to_dicts(self) -> List[Dict]:
        return [d.to_dict() for d in self.items]

    def raise_if_errors(self) -> None:
        if self.errors:
            first = self.errors[0]
            raise TplError(first.message, diagnostics=list(self.items))
