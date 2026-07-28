"""Units as a subsystem, not string matching.

Every physical quantity in TPL carries a unit; a missing one is a parse ERROR
(``PROTOCOL_LANGUAGE.md`` §"Units are mandatory"). This module is the single
place that knows what a unit *is*:

* a **dimension** (what kind of quantity it measures), and
* a **factor to that dimension's canonical unit** (what the engine consumes).

Nothing downstream compares unit strings. A conversion is performed by
:meth:`Quantity.to` / :meth:`Quantity.canonical`, which return the number **and
a human-readable record of the arithmetic**, so every conversion can be shown in
provenance rather than trusted.

Recognised units are exactly the frozen list in the spec::

    mg/m2  mg/kg  mg  uM  nM  mM  um  mm  s  min  h  d  um2/s  %

Canonical units (what the engine's ``core`` dataclasses use):

===============  ==========  =====================================
dimension        canonical   used by
===============  ==========  =====================================
``dose_bsa``     ``mg/m2``   protocol dose, body-surface-area based
``dose_weight``  ``mg/kg``   protocol dose, weight based
``mass``         ``mg``      protocol dose, flat
``concentration````uM``      ``DrugParams.ic50``, ``PKConfig.plasma_peak``
``length``       ``um``      ``dx``
``time``         ``s``       ``PKConfig.dose_interval``, ``PDEConfig.total_time``
``diffusivity``  ``um2/s``   ``DrugParams.diffusivity``
``fraction``     (unitless)  ``%`` -> a fraction in [0, 1]
===============  ==========  =====================================

The three dose dimensions are deliberately kept distinct: ``mg/m2`` and
``mg/kg`` are not interconvertible without a patient, and turning either into a
milligram amount is a modelling assumption (body-surface area / body weight)
that belongs in :mod:`tmesim.protocol.dosing` where it is recorded, not in a
silent unit conversion here.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from .errors import UnitError

__all__ = [
    "UNITS",
    "DIMENSIONS",
    "CANONICAL",
    "UnitDef",
    "Quantity",
    "unit_def",
    "is_unit",
    "dimension_of",
    "units_for",
    "parse_unit",
    "convert",
]


@dataclass(frozen=True)
class UnitDef:
    """One recognised unit."""

    symbol: str
    dimension: str
    #: multiply a value in this unit by ``factor`` to get the canonical unit
    factor: float
    description: str


#: The frozen unit table. Adding a symbol here is a change to the LANGUAGE and
#: must be mirrored in PROTOCOL_LANGUAGE.md and the editor.
UNITS: Dict[str, UnitDef] = {
    # --- dose ------------------------------------------------------------- #
    "mg/m2": UnitDef("mg/m2", "dose_bsa", 1.0,
                     "milligrams per square metre of body-surface area"),
    "mg/kg": UnitDef("mg/kg", "dose_weight", 1.0,
                     "milligrams per kilogram of body weight"),
    "mg": UnitDef("mg", "mass", 1.0, "milligrams (flat amount)"),
    # --- concentration ---------------------------------------------------- #
    "uM": UnitDef("uM", "concentration", 1.0, "micromolar"),
    "nM": UnitDef("nM", "concentration", 1e-3, "nanomolar"),
    "mM": UnitDef("mM", "concentration", 1e3, "millimolar"),
    # --- length ----------------------------------------------------------- #
    "um": UnitDef("um", "length", 1.0, "micrometres"),
    "mm": UnitDef("mm", "length", 1e3, "millimetres"),
    # --- time ------------------------------------------------------------- #
    "s": UnitDef("s", "time", 1.0, "seconds"),
    "min": UnitDef("min", "time", 60.0, "minutes"),
    "h": UnitDef("h", "time", 3600.0, "hours"),
    "d": UnitDef("d", "time", 86400.0, "days"),
    # --- diffusivity ------------------------------------------------------ #
    "um2/s": UnitDef("um2/s", "diffusivity", 1.0,
                     "square micrometres per second"),
    # --- fraction --------------------------------------------------------- #
    "%": UnitDef("%", "fraction", 0.01, "percent (converted to a fraction)"),
}

#: dimension -> canonical unit symbol.
CANONICAL: Dict[str, str] = {
    "dose_bsa": "mg/m2",
    "dose_weight": "mg/kg",
    "mass": "mg",
    "concentration": "uM",
    "length": "um",
    "time": "s",
    "diffusivity": "um2/s",
    "fraction": "",  # dimensionless in [0, 1]
}

#: dimension -> the symbols that carry it (stable order, for error messages).
DIMENSIONS: Dict[str, List[str]] = {}
for _sym, _u in UNITS.items():
    DIMENSIONS.setdefault(_u.dimension, []).append(_sym)

#: The three dose dimensions, which a `dose` statement may use.
DOSE_DIMENSIONS = ("dose_bsa", "dose_weight", "mass")


def is_unit(symbol: str) -> bool:
    """Is ``symbol`` a recognised unit?  Case-sensitive on purpose: ``mM`` and
    ``mm`` are a thousand-fold apart and mean different dimensions."""
    return symbol in UNITS


def unit_def(symbol: str) -> UnitDef:
    """The :class:`UnitDef` for ``symbol``; raises :class:`UnitError`."""
    try:
        return UNITS[symbol]
    except KeyError:
        raise UnitError(
            f"unknown unit {symbol!r}. Recognised units: "
            f"{', '.join(UNITS)}.",
            code="unknown-unit",
            hint=("units are case-sensitive — 'mM' is millimolar, 'mm' is "
                  "millimetres"),
        ) from None


def dimension_of(symbol: str) -> str:
    """The dimension ``symbol`` measures."""
    return unit_def(symbol).dimension


def units_for(dimension: str) -> List[str]:
    """Every recognised symbol carrying ``dimension``."""
    return list(DIMENSIONS.get(dimension, ()))


def parse_unit(symbol: str, dimension: Optional[str] = None,
               line: int = 0, col: int = 0, end_col: int = 0,
               what: str = "value") -> UnitDef:
    """Validate ``symbol``, optionally against a required ``dimension``.

    Raises a positioned :class:`UnitError` naming the acceptable symbols — the
    message is the one the editor shows inline, so it has to be actionable.
    """
    if symbol not in UNITS:
        raise UnitError(
            f"{what}: unknown unit {symbol!r}. Recognised units are "
            f"{', '.join(UNITS)}.",
            line=line, col=col, end_col=end_col, code="unknown-unit",
            hint="units are case-sensitive ('mM' millimolar vs 'mm' millimetres)",
        )
    u = UNITS[symbol]
    if dimension is not None and u.dimension != dimension:
        raise UnitError(
            f"{what}: {symbol!r} measures {u.dimension} but a "
            f"{dimension} unit is required here. Use one of "
            f"{', '.join(units_for(dimension))}.",
            line=line, col=col, end_col=end_col, code="wrong-dimension",
        )
    return u


@dataclass(frozen=True)
class Quantity:
    """A number with a unit and a source position.

    ``unit`` is ``None`` only for a genuinely dimensionless value (a Hill
    coefficient, a repeat count). A *physical* quantity with ``unit=None`` never
    reaches this class — the parser refuses it first.
    """

    value: float
    unit: Optional[str] = None
    line: int = 0
    col: int = 0
    end_col: int = 0

    # -- properties --------------------------------------------------------- #
    @property
    def dimension(self) -> Optional[str]:
        return UNITS[self.unit].dimension if self.unit else None

    @property
    def is_dimensionless(self) -> bool:
        return self.unit is None

    def __str__(self) -> str:
        return f"{_fmt(self.value)} {self.unit}" if self.unit else _fmt(self.value)

    # -- conversion --------------------------------------------------------- #
    def to(self, target: str) -> Tuple[float, str]:
        """Convert to ``target``. Returns ``(value, record)``.

        ``record`` is the arithmetic in words, e.g.
        ``"24 h -> 86400 s (x 3600)"``. It is what goes into provenance; no
        conversion in this package is performed without one.
        """
        if self.unit is None:
            raise UnitError(
                f"cannot convert a dimensionless value {_fmt(self.value)} to "
                f"{target!r}: it has no unit.",
                line=self.line, col=self.col, end_col=self.end_col,
                code="missing-unit",
            )
        src = unit_def(self.unit)
        dst = unit_def(target)
        if src.dimension != dst.dimension:
            raise UnitError(
                f"cannot convert {self.unit!r} ({src.dimension}) to "
                f"{target!r} ({dst.dimension}): different dimensions.",
                line=self.line, col=self.col, end_col=self.end_col,
                code="wrong-dimension",
            )
        ratio = src.factor / dst.factor
        out = float(self.value) * ratio
        if self.unit == target:
            record = f"{_fmt(self.value)} {self.unit} (no conversion)"
        else:
            record = (f"{_fmt(self.value)} {self.unit} -> {_fmt(out)} {target} "
                      f"(x {_fmt(ratio)})")
        return out, record

    def canonical(self) -> Tuple[float, str]:
        """Convert to this dimension's canonical unit. Returns ``(value, record)``.

        ``%`` is the one dimension whose canonical form has no symbol: it
        becomes a plain fraction in [0, 1].
        """
        if self.unit is None:
            return float(self.value), f"{_fmt(self.value)} (dimensionless)"
        dim = self.dimension
        target = CANONICAL[dim]
        if target == "":  # fraction
            out = float(self.value) * UNITS[self.unit].factor
            return out, (f"{_fmt(self.value)}% -> {_fmt(out)} "
                         "(fraction, x 0.01)")
        return self.to(target)

    def value_in(self, target: str) -> float:
        """Just the converted number (the record is discarded)."""
        return self.to(target)[0]


def convert(value: float, source_unit: str, target_unit: str) -> Tuple[float, str]:
    """Free-function form of :meth:`Quantity.to`. Returns ``(value, record)``."""
    return Quantity(float(value), source_unit).to(target_unit)


def _fmt(x: float) -> str:
    """Compact number formatting for provenance strings."""
    v = float(x)
    if v == int(v) and abs(v) < 1e15:
        return str(int(v))
    return f"{v:.6g}"
