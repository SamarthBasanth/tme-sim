"""Recursive-descent parser for TPL -> :mod:`tmesim.protocol.ast`.

Implements the grammar in ``PROTOCOL_LANGUAGE.md`` exactly::

    program     := (drug_def | protocol_def | sweep | run_stmt | compare_stmt
                    | expect_stmt | allow_stmt)*
    drug_def    := "drug" IDENT "from" source NEWLINE field* "end"
    source      := "smiles" STRING | "library" STRING | "file" STRING
    protocol_def:= "protocol" IDENT NEWLINE proto_stmt* "end"
    proto_stmt  := "route" ("oral"|"iv"|"ip") | "dose" NUM UNIT ["over" NUM UNIT]
                 | "every" NUM UNIT | "repeat" INT | "then" "rest" NUM UNIT
                 | "cycles" INT
    sweep       := "sweep" IDENT ("from" NUM "to" NUM UNIT
                                  ("step" NUM | "in" INT "steps") | NUM UNIT)
    run_stmt    := "run" IDENT "with" IDENT "on" IDENT NEWLINE run_opt* "end"
    compare_stmt:= "compare" IDENT ("," IDENT)* "with" IDENT "on" IDENT
                   NEWLINE run_opt* "end"
    expect_stmt := "expect" IDENT COMPARATOR NUM [UNIT]
    allow_stmt  := "allow" "runs" INT

Two properties the editor depends on:

1. **Positions on everything.** Each node records the token that introduced it,
   so a diagnostic can point at ``mg/m3`` rather than "somewhere on line 7".
2. **Error recovery.** A bad statement is reported and the parser resynchronises
   at the next top-level keyword, so one typo does not hide the other four
   problems in the file. :func:`parse` therefore returns *both* a program and a
   diagnostic list, and only raises when there is nothing usable at all.

Units are mandatory on physical quantities and REFUSED on dimensionless ones
(``hill``, ``repeat``). Both directions are parse errors, because both are a
user saying something the engine cannot honour.
"""
from __future__ import annotations

from typing import List, Optional, Tuple

from . import ast as A
from .errors import ParseError, TplError, UnitError, collect
from .lexer import Token, tokenize
from .units import DOSE_DIMENSIONS, UNITS, Quantity, is_unit, unit_def, units_for

__all__ = ["parse", "Parser"]

#: Keywords that may begin a top-level statement (also the resync set).
TOP_KEYWORDS = ("drug", "protocol", "sweep", "run", "compare", "expect", "allow")


def parse(source: str) -> Tuple[A.Program, collect]:
    """Parse ``source``. Returns ``(program, diagnostics)``.

    ``diagnostics.ok`` is ``False`` when the program must not be executed. A
    program is still returned in that case (partially built) so the editor can
    keep highlighting.
    """
    diags = collect()
    try:
        tokens = tokenize(source)
    except TplError as exc:
        for d in exc.diagnostics:
            diags.add(d)
        prog = A.Program(line=1, col=1, end_col=2)
        prog.source = source or ""
        return prog, diags
    return Parser(tokens, source or "", diags).parse_program(), diags


class Parser:
    """One parse of one token stream. Not reusable."""

    def __init__(self, tokens: List[Token], source: str, diags: collect):
        self.toks = tokens
        self.i = 0
        self.source = source
        self.diags = diags

    # ------------------------------------------------------------------ #
    # Token helpers
    # ------------------------------------------------------------------ #
    def peek(self, offset: int = 0) -> Token:
        j = min(self.i + offset, len(self.toks) - 1)
        return self.toks[j]

    def next(self) -> Token:
        tok = self.toks[self.i]
        if tok.kind != "EOF":
            self.i += 1
        return tok

    def at(self, kind: str) -> bool:
        return self.peek().kind == kind

    def at_word(self, *words: str) -> bool:
        return self.peek().is_word(*words)

    def skip_newlines(self) -> None:
        while self.at("NEWLINE"):
            self.next()

    def _fail(self, tok: Token, message: str, code: str = "parse",
              hint: Optional[str] = None) -> "ParseError":
        return ParseError(message, line=tok.line, col=tok.col,
                          end_col=tok.end_col, code=code, hint=hint)

    def expect_newline(self, what: str) -> None:
        tok = self.peek()
        if tok.kind == "NEWLINE":
            self.next()
            return
        if tok.kind == "EOF":
            return
        raise self._fail(
            tok, f"{what}: unexpected {_describe(tok)} at the end of the "
                 "statement. TPL statements end at the newline.",
            code="trailing-tokens")

    def expect_ident(self, what: str) -> Token:
        tok = self.peek()
        if tok.kind != "WORD":
            raise self._fail(tok, f"{what}: expected a name, found "
                                  f"{_describe(tok)}.", code="expected-name")
        self.next()
        return tok

    def expect_keyword(self, word: str, what: str) -> Token:
        tok = self.peek()
        if not tok.is_word(word):
            raise self._fail(tok, f"{what}: expected {word!r}, found "
                                  f"{_describe(tok)}.", code="expected-keyword")
        self.next()
        return tok

    def expect_string(self, what: str) -> Token:
        tok = self.peek()
        if tok.kind != "STRING":
            raise self._fail(
                tok, f"{what}: expected a quoted string, found "
                     f"{_describe(tok)}.", code="expected-string",
                hint='write e.g. from library "cisplatin"')
        self.next()
        return tok

    def expect_number(self, what: str) -> Tuple[float, Token]:
        tok = self.peek()
        if tok.kind != "NUMBER":
            raise self._fail(tok, f"{what}: expected a number, found "
                                  f"{_describe(tok)}.", code="expected-number")
        self.next()
        return float(tok.text), tok

    def expect_integer(self, what: str, minimum: int = 1) -> Tuple[int, Token]:
        value, tok = self.expect_number(what)
        if value != int(value):
            raise self._fail(tok, f"{what}: expected a whole number, got "
                                  f"{tok.text}.", code="expected-integer")
        if int(value) < minimum:
            raise self._fail(tok, f"{what}: must be >= {minimum}, got "
                                  f"{int(value)}.", code="out-of-range")
        return int(value), tok

    # ------------------------------------------------------------------ #
    # Quantities and units
    # ------------------------------------------------------------------ #
    def _peek_unit(self) -> Optional[Token]:
        """The next token if it is a unit (a recognised WORD, or ``%``)."""
        tok = self.peek()
        if tok.kind == "PERCENT":
            return tok
        if tok.kind == "WORD" and is_unit(tok.text):
            return tok
        return None

    def quantity(self, what: str, dimension: Optional[str],
                 unit_required: bool = True) -> Quantity:
        """A number plus its unit.

        ``dimension=None`` means the value is DIMENSIONLESS and a unit is a
        parse error. ``dimension='dose'`` accepts any of ``mg/m2``, ``mg/kg``,
        ``mg``. Otherwise the unit must carry exactly that dimension.
        """
        value, num_tok = self.expect_number(what)
        unit_tok = self._peek_unit()

        if dimension is None:
            if unit_tok is not None:
                raise UnitError(
                    f"{what} is dimensionless and must not carry a unit; "
                    f"found {unit_tok.text!r}.",
                    line=unit_tok.line, col=unit_tok.col,
                    end_col=unit_tok.end_col, code="unexpected-unit",
                    hint=f"write `{what.split(':')[0]} {num_tok.text}`",
                )
            return Quantity(value, None, num_tok.line, num_tok.col,
                            num_tok.end_col)

        if unit_tok is None:
            nxt = self.peek()
            allowed = (", ".join(_dose_units()) if dimension == "dose"
                       else ", ".join(units_for(dimension)))
            found = ("nothing" if nxt.kind in ("NEWLINE", "EOF")
                     else f"{_describe(nxt)}")
            raise UnitError(
                f"{what}: a unit is REQUIRED and none was given (found "
                f"{found}). Every physical quantity in TPL carries a unit — "
                f"`{what.split(':')[0]} {num_tok.text}` is an error. Use one "
                f"of: {allowed}.",
                line=num_tok.line, col=num_tok.col, end_col=num_tok.end_col,
                code="missing-unit",
                hint=f"write `{num_tok.text} {allowed.split(',')[0].strip()}`",
            )

        self.next()
        symbol = "%" if unit_tok.kind == "PERCENT" else unit_tok.text
        if symbol not in UNITS:
            raise UnitError(
                f"{what}: unknown unit {symbol!r}. Recognised units are "
                f"{', '.join(UNITS)}.",
                line=unit_tok.line, col=unit_tok.col, end_col=unit_tok.end_col,
                code="unknown-unit",
                hint="units are case-sensitive ('mM' millimolar vs 'mm' "
                     "millimetres)",
            )
        dim = unit_def(symbol).dimension
        ok = (dim in DOSE_DIMENSIONS) if dimension == "dose" else (dim == dimension)
        if not ok:
            allowed = (", ".join(_dose_units()) if dimension == "dose"
                       else ", ".join(units_for(dimension)))
            raise UnitError(
                f"{what}: {symbol!r} measures {dim}, but this position needs a "
                f"{dimension} unit. Use one of: {allowed}.",
                line=unit_tok.line, col=unit_tok.col, end_col=unit_tok.end_col,
                code="wrong-dimension",
            )
        return Quantity(value, symbol, num_tok.line, num_tok.col,
                        unit_tok.end_col)

    # ------------------------------------------------------------------ #
    # Program
    # ------------------------------------------------------------------ #
    def parse_program(self) -> A.Program:
        prog = A.Program(line=1, col=1, end_col=2)
        prog.source = self.source
        while True:
            self.skip_newlines()
            tok = self.peek()
            if tok.kind == "EOF":
                break
            try:
                self._statement(prog, tok)
            except TplError as exc:
                for d in exc.diagnostics:
                    self.diags.add(d)
                self._resync()
        return prog

    def _statement(self, prog: A.Program, tok: Token) -> None:
        if tok.kind != "WORD":
            raise self._fail(
                tok, f"expected a statement, found {_describe(tok)}. A TPL "
                     f"statement starts with one of: {', '.join(TOP_KEYWORDS)}.",
                code="expected-statement")
        word = tok.text
        if word == "drug":
            prog.drugs.append(self.drug_def())
        elif word == "protocol":
            prog.protocols.append(self.protocol_def())
        elif word == "sweep":
            prog.sweeps.append(self.sweep())
        elif word == "run":
            prog.runs.append(self.run_stmt())
        elif word == "compare":
            prog.runs.append(self.compare_stmt())
        elif word == "expect":
            prog.expects.append(self.expect_stmt())
        elif word == "allow":
            allow = self.allow_stmt()
            if prog.allow_runs is not None:
                raise self._fail(
                    tok, "`allow runs` is given more than once; the sweep "
                         "budget must be stated exactly once.",
                    code="duplicate-allow")
            prog.allow_runs = allow
        else:
            raise self._fail(
                tok, f"unknown statement {word!r}. A TPL statement starts with "
                     f"one of: {', '.join(TOP_KEYWORDS)}. TPL is a protocol "
                     "language — it has no assignments, loops or functions.",
                code="unknown-statement")

    def _resync(self) -> None:
        """Skip to the start of the next plausible top-level statement."""
        depth_end_seen = False
        while True:
            tok = self.peek()
            if tok.kind == "EOF":
                return
            if tok.kind == "NEWLINE":
                self.next()
                nxt = self.peek()
                if nxt.kind == "EOF":
                    return
                if nxt.kind == "WORD" and nxt.text in TOP_KEYWORDS:
                    return
                if depth_end_seen:
                    return
                continue
            if tok.is_word("end"):
                depth_end_seen = True
            self.next()

    # ------------------------------------------------------------------ #
    # drug
    # ------------------------------------------------------------------ #
    def drug_def(self) -> A.DrugDef:
        kw = self.next()  # 'drug'
        name_tok = self.expect_ident("drug: the name after `drug`")
        _reject_keyword_as_name(self, name_tok, "drug")
        self.expect_keyword("from", f"drug {name_tok.text}")
        src_tok = self.peek()
        if not (src_tok.kind == "WORD" and src_tok.text in A.DRUG_SOURCES):
            raise self._fail(
                src_tok,
                f"drug {name_tok.text}: source must be one of "
                f"{', '.join(A.DRUG_SOURCES)}, found {_describe(src_tok)}.",
                code="unknown-drug-source",
                hint='e.g. `from smiles "CN1C(=O)N2C=NC(=C2N=N1)C(=O)N"`')
        self.next()
        val_tok = self.expect_string(f"drug {name_tok.text}: `from "
                                     f"{src_tok.text}`")
        node = A.DrugDef(line=kw.line, col=kw.col, end_col=kw.end_col,
                         name=name_tok.text, source_kind=src_tok.text,
                         source_value=val_tok.text,
                         source_line=src_tok.line, source_col=src_tok.col)
        self.expect_newline(f"drug {name_tok.text}")

        seen: dict = {}
        while True:
            self.skip_newlines()
            tok = self.peek()
            if tok.kind == "EOF":
                raise self._fail(
                    tok, f"drug {name_tok.text}: block is never closed — "
                         "add `end`.", code="unclosed-block")
            if tok.is_word("end"):
                self.next()
                self.expect_newline("end")
                node.end_col = tok.end_col
                return node
            if tok.kind != "WORD":
                raise self._fail(
                    tok, f"drug {name_tok.text}: expected a field name or "
                         f"`end`, found {_describe(tok)}.",
                    code="expected-field")
            fld = self.drug_field(name_tok.text)
            if fld.name in seen:
                raise self._fail(
                    tok, f"drug {name_tok.text}: field {fld.name!r} is set "
                         f"twice (first on line {seen[fld.name]}).",
                    code="duplicate-field")
            seen[fld.name] = tok.line
            node.fields.append(fld)

    def drug_field(self, drug_name: str) -> A.DrugField:
        tok = self.next()
        fname = tok.text
        if fname not in A.OVERRIDABLE_FIELDS:
            raise self._fail(
                tok, f"drug {drug_name}: {fname!r} is not an overridable "
                     f"field. Overridable: "
                     f"{', '.join(A.OVERRIDABLE_FIELDS)}.",
                code="unknown-field")
        spec = A.OVERRIDABLE_FIELDS[fname]

        # `ic50 from gdsc`
        if fname == "ic50" and self.at_word("from"):
            self.next()
            g = self.peek()
            if not g.is_word("gdsc"):
                raise self._fail(
                    g, "ic50: the only recognised measured source is `gdsc` "
                       f"(the cached GDSC glioma screen), found "
                       f"{_describe(g)}.", code="unknown-ic50-source")
            self.next()
            self.expect_newline("ic50 from gdsc")
            return A.DrugField(
                line=tok.line, col=tok.col, end_col=g.end_col, name="ic50",
                source=A.Ic50FromGdsc(line=tok.line, col=tok.col,
                                      end_col=g.end_col))

        q = self.quantity(f"{fname}: {drug_name}.{fname}", spec["dimension"])
        self.expect_newline(fname)
        return A.DrugField(line=tok.line, col=tok.col, end_col=q.end_col,
                           name=fname, value=q)

    # ------------------------------------------------------------------ #
    # protocol
    # ------------------------------------------------------------------ #
    def protocol_def(self) -> A.ProtocolDef:
        kw = self.next()
        name_tok = self.expect_ident("protocol: the name after `protocol`")
        _reject_keyword_as_name(self, name_tok, "protocol")
        node = A.ProtocolDef(line=kw.line, col=kw.col, end_col=kw.end_col,
                             name=name_tok.text)
        self.expect_newline(f"protocol {name_tok.text}")
        label = f"protocol {name_tok.text}"

        while True:
            self.skip_newlines()
            tok = self.peek()
            if tok.kind == "EOF":
                raise self._fail(
                    tok, f"{label}: block is never closed — add `end`.",
                    code="unclosed-block")
            if tok.is_word("end"):
                self.next()
                self.expect_newline("end")
                node.end_col = tok.end_col
                return node
            if tok.kind != "WORD":
                raise self._fail(
                    tok, f"{label}: expected a protocol statement or `end`, "
                         f"found {_describe(tok)}.", code="expected-statement")
            self._proto_stmt(node, label)

    def _proto_stmt(self, node: A.ProtocolDef, label: str) -> None:
        tok = self.next()
        word = tok.text

        def _once(attr: str) -> None:
            existing = getattr(node, attr)
            if existing is not None:
                raise self._fail(
                    tok, f"{label}: `{word}` is given more than once "
                         f"(first on line {existing.line}).",
                    code="duplicate-statement")

        if word == "route":
            _once("route")
            r = self.peek()
            if not (r.kind == "WORD" and r.text in A.ROUTES):
                raise self._fail(
                    r, f"{label}: route must be one of "
                       f"{', '.join(A.ROUTES)}, found {_describe(r)}.",
                    code="unknown-route")
            self.next()
            self.expect_newline("route")
            node.route = A.Route(line=tok.line, col=tok.col,
                                 end_col=r.end_col, value=r.text)
        elif word == "dose":
            _once("dose")
            amount = self.quantity(f"{label}: dose", "dose")
            over = None
            if self.at_word("over"):
                over_tok = self.next()
                over = self.quantity(f"{label}: dose ... over", "time")
                if over.value <= 0:
                    raise self._fail(
                        over_tok, f"{label}: infusion duration must be > 0.",
                        code="out-of-range")
            self.expect_newline("dose")
            if amount.value <= 0:
                raise self._fail(tok, f"{label}: dose must be > 0.",
                                 code="out-of-range")
            node.dose = A.Dose(line=tok.line, col=tok.col,
                               end_col=(over or amount).end_col,
                               amount=amount, over=over)
        elif word == "every":
            _once("every")
            q = self.quantity(f"{label}: every", "time")
            self.expect_newline("every")
            if q.value <= 0:
                raise self._fail(tok, f"{label}: `every` must be > 0.",
                                 code="out-of-range")
            node.every = A.Every(line=tok.line, col=tok.col, end_col=q.end_col,
                                 interval=q)
        elif word == "repeat":
            _once("repeat")
            n, ntok = self.expect_integer(f"{label}: repeat", minimum=1)
            self.expect_newline("repeat")
            node.repeat = A.Repeat(line=tok.line, col=tok.col,
                                   end_col=ntok.end_col, count=n)
        elif word == "cycles":
            _once("cycles")
            n, ntok = self.expect_integer(f"{label}: cycles", minimum=1)
            self.expect_newline("cycles")
            node.cycles = A.Cycles(line=tok.line, col=tok.col,
                                   end_col=ntok.end_col, count=n)
        elif word == "then":
            _once("rest")
            self.expect_keyword("rest", f"{label}: `then`")
            q = self.quantity(f"{label}: then rest", "time")
            self.expect_newline("then rest")
            if q.value <= 0:
                raise self._fail(tok, f"{label}: rest must be > 0.",
                                 code="out-of-range")
            node.rest = A.Rest(line=tok.line, col=tok.col, end_col=q.end_col,
                               duration=q)
        elif word == "rest":
            raise self._fail(
                tok, f"{label}: a rest period is written `then rest {'{N}'} d`.",
                code="expected-keyword", hint="add the `then`")
        else:
            raise self._fail(
                tok, f"{label}: unknown protocol statement {word!r}. Valid: "
                     "route, dose, every, repeat, then rest, cycles.",
                code="unknown-statement")

    # ------------------------------------------------------------------ #
    # sweep
    # ------------------------------------------------------------------ #
    def sweep(self) -> A.Sweep:
        kw = self.next()
        p = self.expect_ident("sweep: the parameter after `sweep`")
        if p.text not in A.SWEEPABLE:
            raise self._fail(
                p, f"sweep: {p.text!r} is not a sweepable parameter. "
                   f"Sweepable: {', '.join(A.SWEEPABLE)}.",
                code="unknown-sweep-param")
        spec = A.SWEEPABLE[p.text]
        dim = spec["dimension"]
        node = A.Sweep(line=kw.line, col=kw.col, end_col=p.end_col,
                       param=p.text, param_col=p.col)

        if self.at_word("from"):
            self.next()
            start_val, start_tok = self.expect_number(
                f"sweep {p.text}: the `from` value")
            self.expect_keyword("to", f"sweep {p.text}")
            stop = self.quantity(f"sweep {p.text}: the `to` value", dim)
            node.kind = "range"
            node.start = Quantity(start_val, stop.unit, start_tok.line,
                                  start_tok.col, start_tok.end_col)
            node.stop = stop
            nxt = self.peek()
            if nxt.is_word("step"):
                self.next()
                step, step_tok = self.expect_number(f"sweep {p.text}: `step`")
                if step <= 0:
                    raise self._fail(step_tok,
                                     f"sweep {p.text}: `step` must be > 0.",
                                     code="out-of-range")
                node.step = step
                node.end_col = step_tok.end_col
            elif nxt.is_word("in"):
                self.next()
                n, ntok = self.expect_integer(f"sweep {p.text}: `in N steps`",
                                              minimum=1)
                self.expect_keyword("steps", f"sweep {p.text}: `in {n}`")
                node.n_steps = n
                node.end_col = ntok.end_col
            else:
                raise self._fail(
                    nxt, f"sweep {p.text}: a range needs either `step S` or "
                         f"`in N steps`, found {_describe(nxt)}.",
                    code="incomplete-sweep",
                    hint="`in N steps` is linear; `step S` walks by S")
            if node.stop.value == node.start.value:
                raise self._fail(
                    stop.line and start_tok or start_tok,
                    f"sweep {p.text}: `from` and `to` are the same value "
                    f"({start_val}); use `sweep {p.text} {start_val}"
                    f"{' ' + (stop.unit or '') if stop.unit else ''}` for a "
                    "fixed override.", code="empty-sweep")
        else:
            fixed = self.quantity(f"sweep {p.text}", dim)
            node.kind = "fixed"
            node.fixed = fixed
            node.end_col = fixed.end_col
        self.expect_newline(f"sweep {p.text}")
        if spec.get("integer"):
            for q in (node.start, node.stop, node.fixed):
                if q is not None and q.value != int(q.value):
                    raise self._fail(
                        Token("NUMBER", str(q.value), q.line, q.col, q.end_col),
                        f"sweep {p.text}: {p.text} is a whole count, got "
                        f"{q.value}.", code="expected-integer")
        return node

    # ------------------------------------------------------------------ #
    # run / compare
    # ------------------------------------------------------------------ #
    def run_stmt(self) -> A.RunStmt:
        kw = self.next()
        d = self.expect_ident("run: the drug after `run`")
        self.expect_keyword("with", f"run {d.text}")
        p = self.expect_ident(f"run {d.text}: the protocol after `with`")
        self.expect_keyword("on", f"run {d.text} with {p.text}")
        s = self.expect_ident(f"run {d.text}: the specimen after `on`")
        node = A.RunStmt(line=kw.line, col=kw.col, end_col=s.end_col,
                         drug=d.text, drug_col=d.col, protocol=p.text,
                         protocol_col=p.col, section=s.text, section_col=s.col)
        self.expect_newline("run")
        node.options = self._run_body(f"run {d.text}", allow_rank=False)
        return node

    def compare_stmt(self) -> A.CompareStmt:
        kw = self.next()
        drugs: List[str] = []
        cols: List[int] = []
        while True:
            d = self.expect_ident("compare: a drug name")
            drugs.append(d.text)
            cols.append(d.col)
            if self.at("COMMA"):
                self.next()
                continue
            break
        if len(drugs) < 2:
            raise self._fail(
                self.peek(),
                "compare: needs at least two drugs (`compare A, B with ...`); "
                "use `run` for a single drug.", code="incomplete-compare")
        self.expect_keyword("with", "compare")
        p = self.expect_ident("compare: the protocol after `with`")
        self.expect_keyword("on", "compare")
        s = self.expect_ident("compare: the specimen after `on`")
        node = A.CompareStmt(line=kw.line, col=kw.col, end_col=s.end_col,
                             drugs=drugs, drug_cols=cols, protocol=p.text,
                             protocol_col=p.col, section=s.text,
                             section_col=s.col)
        self.expect_newline("compare")
        node.options = self._run_body("compare", allow_rank=True)
        return node

    def _run_body(self, label: str, allow_rank: bool) -> A.RunOptions:
        opts = A.RunOptions()
        seen: dict = {}
        while True:
            self.skip_newlines()
            tok = self.peek()
            if tok.kind == "EOF":
                raise self._fail(
                    tok, f"{label}: block is never closed — add `end`.",
                    code="unclosed-block")
            if tok.is_word("end"):
                self.next()
                self.expect_newline("end")
                return opts
            if tok.kind != "WORD":
                raise self._fail(
                    tok, f"{label}: expected a run option or `end`, found "
                         f"{_describe(tok)}.", code="expected-option")
            self.next()
            word = tok.text
            valid = ["dx", "state_mode", "sensitivity", "report"]
            if allow_rank:
                valid.append("rank by")
            if word in seen:
                raise self._fail(
                    tok, f"{label}: `{word}` is given more than once (first on "
                         f"line {seen[word]}).", code="duplicate-option")
            seen[word] = tok.line

            if word == "dx":
                opts.dx = self.quantity(f"{label}: dx", "length")
                self.expect_newline("dx")
            elif word == "state_mode":
                v = self.expect_ident(f"{label}: state_mode")
                if v.text not in A.STATE_MODES:
                    raise self._fail(
                        v, f"{label}: state_mode {v.text!r} is not one of "
                           f"{', '.join(A.STATE_MODES)}.",
                        code="unknown-state-mode")
                opts.state_mode = v.text
                opts.state_mode_line, opts.state_mode_col = v.line, v.col
                self.expect_newline("state_mode")
            elif word == "sensitivity":
                v = self.expect_ident(f"{label}: sensitivity")
                if v.text not in A.SENSITIVITY_MODES:
                    raise self._fail(
                        v, f"{label}: sensitivity {v.text!r} is not one of "
                           f"{', '.join(A.SENSITIVITY_MODES)}.",
                        code="unknown-sensitivity-mode")
                opts.sensitivity = v.text
                opts.sensitivity_line, opts.sensitivity_col = v.line, v.col
                self.expect_newline("sensitivity")
            elif word == "report":
                opts.report_line = tok.line
                while True:
                    m = self.expect_ident(f"{label}: report")
                    opts.report.append(_metric(self, m))
                    if self.at("COMMA"):
                        self.next()
                        continue
                    break
                self.expect_newline("report")
            elif word == "rank" and allow_rank:
                self.expect_keyword("by", f"{label}: `rank`")
                m = self.expect_ident(f"{label}: rank by")
                opts.rank_by = _metric(self, m)
                opts.rank_by_line, opts.rank_by_col = m.line, m.col
                self.expect_newline("rank by")
            else:
                raise self._fail(
                    tok, f"{label}: unknown option {word!r}. Valid options: "
                         f"{', '.join(valid)}.", code="unknown-option")

    # ------------------------------------------------------------------ #
    # expect / allow
    # ------------------------------------------------------------------ #
    def expect_stmt(self) -> A.ExpectStmt:
        kw = self.next()
        m = self.expect_ident("expect: the metric after `expect`")
        metric = _metric(self, m)
        c = self.peek()
        if c.kind != "COMPARATOR":
            raise self._fail(
                c, f"expect {metric}: expected a comparator "
                   f"(>, <, >=, <=, ==, !=), found {_describe(c)}.",
                code="expected-comparator")
        self.next()
        value, vtok = self.expect_number(f"expect {metric}")
        unit_tok = self._peek_unit()
        unit = None
        if unit_tok is not None:
            self.next()
            unit = "%" if unit_tok.kind == "PERCENT" else unit_tok.text
        self.expect_newline("expect")
        end = (unit_tok or vtok).end_col
        return A.ExpectStmt(line=kw.line, col=kw.col, end_col=end,
                            metric=metric, metric_col=m.col,
                            comparator=c.text,
                            value=Quantity(value, unit, vtok.line, vtok.col,
                                           end))

    def allow_stmt(self) -> A.AllowRuns:
        kw = self.next()
        self.expect_keyword("runs", "allow")
        n, ntok = self.expect_integer("allow runs", minimum=1)
        self.expect_newline("allow runs")
        return A.AllowRuns(line=kw.line, col=kw.col, end_col=ntok.end_col,
                           count=n)


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _metric(parser: Parser, tok: Token) -> str:
    name = tok.text
    if name in A.METRICS:
        return name
    if name in A.METRIC_ALIASES:
        return A.METRIC_ALIASES[name]
    raise parser._fail(
        tok, f"unknown metric {name!r}. Reportable metrics: "
             f"{', '.join(A.METRICS)}.", code="unknown-metric")


def _reject_keyword_as_name(parser: Parser, tok: Token, what: str) -> None:
    if tok.text in TOP_KEYWORDS or tok.text == "end":
        raise parser._fail(
            tok, f"{what}: {tok.text!r} is a reserved word and cannot be used "
                 "as a name.", code="reserved-name")


def _dose_units() -> List[str]:
    out: List[str] = []
    for dim in DOSE_DIMENSIONS:
        out.extend(units_for(dim))
    return out


def _describe(tok: Token) -> str:
    if tok.kind == "EOF":
        return "end of file"
    if tok.kind == "NEWLINE":
        return "end of line"
    if tok.kind == "STRING":
        return f"the string {tok.text!r}"
    if tok.kind == "NUMBER":
        return f"the number {tok.text}"
    return f"{tok.text!r}"
