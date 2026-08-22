"""TPL tokenizer.

Line-oriented, because TPL is: statements end at a newline and blocks close with
``end``. Indentation is not significant. Comments are ``#`` to end of line.

Every token carries ``line`` (1-based), ``col`` (1-based) and ``end_col``
(1-based, exclusive) so the parser can attribute any refusal to the exact
characters that caused it.

Token kinds
-----------
``WORD``
    ``[A-Za-z_][A-Za-z0-9_.-]*`` optionally followed by ONE ``/`` and a second
    such run. That single rule covers identifiers (``stupp_adjuvant``),
    keywords (``dose``), section ids (``xenium_319``) and the compound units
    (``mg/m2``, ``um2/s``) — the parser, not the lexer, decides which a WORD is,
    because ``d`` is both "days" and a legal identifier.
``NUMBER``
    integer / decimal / scientific. A leading ``-`` is part of the number.
``STRING``
    double- or single-quoted, with ``\\`` escapes. An unterminated string is a
    positioned :class:`~tmesim.protocol.errors.LexError`, not a crash.
``PERCENT`` ``COMMA`` ``COMPARATOR`` ``NEWLINE`` ``EOF``
    one-character / operator tokens. ``COMPARATOR`` is ``>= <= == != > <``.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from typing import List, Optional

from .errors import LexError

__all__ = ["Token", "tokenize", "KEYWORDS", "COMPARATORS"]

#: Words the grammar gives meaning to. A WORD that is not one of these is an
#: identifier (or a unit). Listed for the editor's syntax highlighter.
KEYWORDS = frozenset({
    "drug", "from", "smiles", "library", "file", "end",
    "protocol", "route", "dose", "over", "every", "repeat", "then", "rest",
    "cycles", "sweep", "to", "step", "in", "steps",
    "run", "with", "on", "compare", "expect", "allow", "runs",
    "ic50", "gdsc", "report", "dx", "state_mode", "sensitivity",
})
# `rank` and `by` were in this set until 2026-08-21. `rank by` was withdrawn
# from TPL (PROTOCOL_LANGUAGE.md §5) and `by` had no other use, so neither word
# carries grammar meaning any more and neither must highlight as live syntax.
# `tokenize` does not read this set — the parser still recognises the WORD
# `rank` in a run body so it can refuse it by name rather than as a typo.

COMPARATORS = (">=", "<=", "==", "!=", ">", "<")

_NUMBER = re.compile(r"-?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][+-]?\d+)?")
_WORD = re.compile(r"[A-Za-z_][A-Za-z0-9_.\-]*(?:/[A-Za-z0-9_][A-Za-z0-9_.\-]*)?")


@dataclass(frozen=True)
class Token:
    kind: str
    text: str
    line: int
    col: int
    end_col: int

    @property
    def value(self) -> str:
        return self.text

    def is_word(self, *words: str) -> bool:
        return self.kind == "WORD" and self.text in words

    def __str__(self) -> str:  # pragma: no cover - diagnostics only
        return f"{self.kind}({self.text!r}) @ {self.line}:{self.col}"


def tokenize(source: str) -> List[Token]:
    """Tokenize ``source``. Raises :class:`LexError` with a position.

    A single trailing ``NEWLINE`` is guaranteed before ``EOF`` so the parser can
    treat "end of statement" uniformly.
    """
    if source is None:
        source = ""
    tokens: List[Token] = []
    for lineno, raw in enumerate(source.replace("\r\n", "\n").split("\n"), start=1):
        i = 0
        n = len(raw)
        emitted = False
        while i < n:
            ch = raw[i]
            if ch in " \t":
                i += 1
                continue
            if ch == "#":
                break  # comment to end of line
            col = i + 1

            # --- string --------------------------------------------------- #
            if ch in "\"'":
                text, consumed = _scan_string(raw, i, lineno)
                tokens.append(Token("STRING", text, lineno, col, col + consumed))
                i += consumed
                emitted = True
                continue

            # --- number (a '-' only starts one if a value can start here) -- #
            if ch.isdigit() or (ch == "." and i + 1 < n and raw[i + 1].isdigit()) \
                    or (ch == "-" and i + 1 < n and (raw[i + 1].isdigit() or raw[i + 1] == ".")
                        and _minus_starts_a_number(tokens)):
                m = _NUMBER.match(raw, i)
                if m:
                    tokens.append(Token("NUMBER", m.group(0), lineno, col,
                                        col + len(m.group(0))))
                    i = m.end()
                    emitted = True
                    continue

            # --- comparator ------------------------------------------------ #
            hit = _match_comparator(raw, i)
            if hit:
                tokens.append(Token("COMPARATOR", hit, lineno, col,
                                    col + len(hit)))
                i += len(hit)
                emitted = True
                continue

            # --- word / unit ----------------------------------------------- #
            m = _WORD.match(raw, i)
            if m:
                tokens.append(Token("WORD", m.group(0), lineno, col,
                                    col + len(m.group(0))))
                i = m.end()
                emitted = True
                continue

            # --- punctuation ------------------------------------------------ #
            if ch == "%":
                tokens.append(Token("PERCENT", "%", lineno, col, col + 1))
                i += 1
                emitted = True
                continue
            if ch == ",":
                tokens.append(Token("COMMA", ",", lineno, col, col + 1))
                i += 1
                emitted = True
                continue

            raise LexError(
                f"unexpected character {ch!r}. TPL has no operators beyond "
                "comparators in `expect`; a value is a number with a unit, a "
                "quoted string, or a name.",
                line=lineno, col=col, end_col=col + 1, code="unexpected-char",
            )
        if emitted:
            tokens.append(Token("NEWLINE", "\\n", lineno, len(raw) + 1,
                                len(raw) + 2))
    last_line = tokens[-1].line if tokens else 1
    if not tokens or tokens[-1].kind != "NEWLINE":
        tokens.append(Token("NEWLINE", "\\n", last_line, 1, 2))
    tokens.append(Token("EOF", "", last_line + 1, 1, 2))
    return tokens


def _minus_starts_a_number(tokens: List[Token]) -> bool:
    """``-`` begins a negative literal unless it directly follows a value.

    TPL has no arithmetic, so this only ever has to distinguish ``step -50``
    from a stray dash; but keeping the rule explicit means ``expect kill > -1``
    lexes and ``foo-bar`` stays one identifier (handled by ``_WORD``).
    """
    if not tokens:
        return True
    return tokens[-1].kind not in ("NUMBER", "STRING")


def _match_comparator(raw: str, i: int) -> Optional[str]:
    for op in COMPARATORS:
        if raw.startswith(op, i):
            return op
    return None


def _scan_string(raw: str, i: int, lineno: int) -> "tuple[str, int]":
    """Scan a quoted string starting at ``raw[i]``. Returns (value, n_consumed)."""
    quote = raw[i]
    out: List[str] = []
    j = i + 1
    while j < len(raw):
        c = raw[j]
        if c == "\\" and j + 1 < len(raw):
            nxt = raw[j + 1]
            out.append({"n": "\n", "t": "\t", "\\": "\\",
                        '"': '"', "'": "'"}.get(nxt, "\\" + nxt))
            j += 2
            continue
        if c == quote:
            return "".join(out), (j - i + 1)
        out.append(c)
        j += 1
    raise LexError(
        f"unterminated string: no closing {quote} on this line.",
        line=lineno, col=i + 1, end_col=len(raw) + 1, code="unterminated-string",
        hint="TPL strings do not span lines",
    )
