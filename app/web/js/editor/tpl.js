// tpl.js — the TPL lexer and block model, implemented against PROTOCOL_LANGUAGE.md v1.
//
// This module knows the SHAPE of the language (what a token is, where a block
// opens and closes, how a line indents). It deliberately does NOT know whether
// a program is *valid*: that is the executor's job and it is asked over
// POST /api/protocol/validate. Duplicating semantic rules on the client is how
// an editor ends up disagreeing with the engine, and when those two disagree
// the operator believes the one that is wrong.
//
// It does track block STRUCTURE, because the editor needs it for three things
// that are not validation at all: matching an `end` to its opener for the
// block-highlight, deciding what a new line indents to, and snapping a typed
// `end` back to its opener's column.
//
// Units are lexed as their OWN token class, not folded into the number, because
// in TPL the unit is semantic: `dose 150` is an error and `dose 150 mg/m2` is
// not. The editor colours them separately for exactly that reason.

export const INDENT = '  ';

/* ------------------------------------------------------------ vocabulary -- */

/** Structural keywords: statement heads, block delimiters, connectors. */
export const KEYWORDS = new Set([
  'drug', 'protocol', 'sweep', 'run', 'compare', 'expect', 'allow', 'end',
  // No `by`: its only use was `rank by`, withdrawn 2026-08-21 with the rest of
  // compound ordering (PROTOCOL_LANGUAGE.md §5). Mirrors tmesim/protocol/lexer.py.
  'from', 'with', 'on', 'over', 'to', 'step', 'in', 'steps', 'then',
  'smiles', 'library', 'file',
]);

/** Field / option names — the things that take a value. */
export const FIELDS = new Set([
  // drug block
  'ic50', 'hill', 'diffusivity', 'binding', 'clearance', 'bbb_penetrance',
  'molecular_weight',
  // protocol block
  'route', 'dose', 'every', 'repeat', 'rest', 'cycles',
  // run / compare block. No `rank`: `rank by` was withdrawn from TPL on
  // 2026-08-21 (PROTOCOL_LANGUAGE.md §5), so it must not highlight as a live
  // option — the parser refuses it and the highlighter has to agree.
  'dx', 'state_mode', 'sensitivity', 'report', 'runs',
]);

/** Enumerated values and reportable metric names. */
export const ENUMS = new Set([
  'oral', 'iv', 'ip',
  'uniform', 'determinant',
  'oxygen', 'vessel_distance', 'markers', 'pde',
  'discrete', 'kde',
  'gdsc',
  'coverage', 'kill', 'mean_auc', 'penetration_length',
]);

/**
 * Recognised units, longest first so that `mg/m2` is never lexed as `mg`
 * followed by garbage. This list is the spec's list, verbatim.
 */
export const UNITS = [
  'mg/m2', 'mg/kg', 'um2/s', 'min', 'uM', 'nM', 'mM', 'mg', 'um', 'mm',
  '%', 's', 'h', 'd',
];

/** Fields whose value is a bare count, so a missing unit is NOT an error. */
export const UNITLESS_FIELDS = new Set([
  'repeat', 'cycles', 'runs', 'hill', 'binding', 'bbb_penetrance',
]);

export const BLOCK_OPENERS = new Set(['drug', 'protocol', 'run', 'compare']);

/* ---------------------------------------------------------------- lexer -- */

const ID_CH = /[A-Za-z0-9_]/;
const ID_START = /[A-Za-z_]/;
const DIGIT = /[0-9]/;

function matchUnit(line, i) {
  for (const u of UNITS) {
    if (!line.startsWith(u, i)) continue;
    const nxt = line[i + u.length];
    if (nxt !== undefined && ID_CH.test(nxt)) continue;   // `4 steps` is not `4 s`
    return u;
  }
  return null;
}

/**
 * Lex ONE line into tokens. Tokens tile the line exactly, so a renderer can
 * walk them and reproduce the source character for character.
 *
 * @returns {{t:string,s:number,e:number,v:string}[]}
 *   t ∈ ws | comment | string | number | unit | kw | fld | enum | ident | op | bad
 */
export function tokenizeLine(line) {
  const out = [];
  let i = 0;
  // True while we are still in the "a unit may legally follow" window, i.e.
  // directly after a number, across intervening spaces only.
  let afterNumber = false;

  while (i < line.length) {
    const c = line[i];

    if (c === '#') {
      out.push({ t: 'comment', s: i, e: line.length, v: line.slice(i) });
      break;
    }

    if (c === ' ' || c === '\t') {
      let j = i;
      while (j < line.length && (line[j] === ' ' || line[j] === '\t')) j++;
      out.push({ t: 'ws', s: i, e: j, v: line.slice(i, j) });
      i = j;
      continue;                                   // whitespace does not close the unit window
    }

    if (c === '"') {
      let j = i + 1;
      while (j < line.length && line[j] !== '"') j++;
      const closed = line[j] === '"';
      const e = closed ? j + 1 : line.length;
      out.push({ t: closed ? 'string' : 'bad', s: i, e, v: line.slice(i, e) });
      i = e;
      afterNumber = false;
      continue;
    }

    if (DIGIT.test(c) || ((c === '-' || c === '+' || c === '.') && DIGIT.test(line[i + 1] || ''))) {
      const m = /^[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?/.exec(line.slice(i));
      const e = i + (m ? m[0].length : 1);
      out.push({ t: 'number', s: i, e, v: line.slice(i, e) });
      i = e;
      afterNumber = true;
      continue;
    }

    if (afterNumber) {
      const u = matchUnit(line, i);
      if (u) {
        out.push({ t: 'unit', s: i, e: i + u.length, v: u });
        i += u.length;
        afterNumber = false;
        continue;
      }
    }

    if (ID_START.test(c)) {
      let j = i;
      while (j < line.length && ID_CH.test(line[j])) j++;
      const w = line.slice(i, j);
      const t = KEYWORDS.has(w) ? 'kw'
        : FIELDS.has(w) ? 'fld'
          : ENUMS.has(w) ? 'enum'
            : 'ident';
      out.push({ t, s: i, e: j, v: w });
      i = j;
      afterNumber = false;
      continue;
    }

    if ('<>=!,:/*+-'.includes(c)) {
      let j = i + 1;
      if ((c === '<' || c === '>' || c === '=' || c === '!') && line[j] === '=') j++;
      out.push({ t: 'op', s: i, e: j, v: line.slice(i, j) });
      i = j;
      afterNumber = false;
      continue;
    }

    out.push({ t: 'bad', s: i, e: i + 1, v: c });
    i += 1;
    afterNumber = false;
  }
  return out;
}

/** Lex a whole buffer. Returns one token array per line. */
export function tokenize(text) {
  return text.split('\n').map(tokenizeLine);
}

/** First token that is not whitespace, or null for a blank / comment-only line. */
export function headToken(tokens) {
  for (const t of tokens) {
    if (t.t === 'ws') continue;
    if (t.t === 'comment') return null;
    return t;
  }
  return null;
}

/* ------------------------------------------------------------- structure -- */

/**
 * Pair up block openers with their `end`.
 *
 * @returns {{
 *   partner: Map<number, number>,   // 1-based line -> its matching line
 *   depth: number[],                // 0-based index -> nesting depth of that line
 *   errors: {line:number, col:number, end_col:number, message:string, severity:string, source:string}[]
 * }}
 */
export function blockModel(text) {
  const lines = text.split('\n');
  const stack = [];
  const partner = new Map();
  const depth = new Array(lines.length).fill(0);
  const errors = [];

  lines.forEach((line, idx) => {
    const toks = tokenizeLine(line);
    const head = headToken(toks);
    depth[idx] = stack.length;
    if (!head || head.t !== 'kw') return;

    if (head.v === 'end') {
      const open = stack.pop();
      depth[idx] = stack.length;
      if (open === undefined) {
        errors.push({
          line: idx + 1, col: head.s, end_col: head.e, severity: 'error', source: 'local',
          message: '`end` closes nothing — there is no open drug / protocol / run / compare block above it.',
        });
      } else {
        partner.set(open.line, idx + 1);
        partner.set(idx + 1, open.line);
      }
      return;
    }

    if (BLOCK_OPENERS.has(head.v)) stack.push({ line: idx + 1, kw: head.v, col: head.s });
  });

  for (const open of stack) {
    errors.push({
      line: open.line, col: open.col, end_col: open.col + open.kw.length,
      severity: 'error', source: 'local',
      message: `\`${open.kw}\` block is never closed — every block ends with \`end\` on its own line.`,
    });
  }
  return { partner, depth, errors };
}

/* ------------------------------------------------------------ indentation -- */

export const indentOf = (line) => (/^[ \t]*/.exec(line) || [''])[0].replace(/\t/g, INDENT);

/**
 * The indent a NEW line should get, given the line the caret just left.
 * Opening a block indents one level; nothing else changes the level, because
 * TPL indentation is cosmetic and the language says so.
 */
export function indentAfter(line) {
  const head = headToken(tokenizeLine(line));
  const base = indentOf(line);
  if (head && head.t === 'kw' && BLOCK_OPENERS.has(head.v)) return base + INDENT;
  return base;
}

/** The indent the line containing `end` should snap back to. */
export function indentForEnd(text, lineNo) {
  const lines = text.split('\n');
  let level = 0;
  for (let i = lineNo - 2; i >= 0; i--) {
    const head = headToken(tokenizeLine(lines[i]));
    if (!head || head.t !== 'kw') continue;
    if (head.v === 'end') { level += 1; continue; }
    if (BLOCK_OPENERS.has(head.v)) {
      if (level === 0) return indentOf(lines[i]);
      level -= 1;
    }
  }
  return '';
}

/* --------------------------------------------------------------- symbols -- */

/**
 * Names the program declares, so the editor can offer them and the panel can
 * say what a plan will touch. Purely lexical — no semantic checking.
 */
export function symbols(text) {
  const drugs = [];
  const protocols = [];
  const sweeps = [];
  text.split('\n').forEach((line, idx) => {
    const toks = tokenizeLine(line).filter((t) => t.t !== 'ws' && t.t !== 'comment');
    const head = toks[0];
    if (!head || head.t !== 'kw') return;
    if (head.v === 'drug' && toks[1]) drugs.push({ name: toks[1].v, line: idx + 1 });
    if (head.v === 'protocol' && toks[1]) protocols.push({ name: toks[1].v, line: idx + 1 });
    if (head.v === 'sweep' && toks[1]) sweeps.push({ name: toks[1].v, line: idx + 1 });
  });
  return { drugs, protocols, sweeps };
}
