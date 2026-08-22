// protocol/stub.js — a LOCAL stand-in for the TPL parser/executor.
//
// It exists so the protocol editor is exercisable before `tmesim/protocol/` and
// the `/api/protocol/*` routes land, and so the workstation degrades honestly
// rather than dying when they are absent. Everything it produces is labelled
// STUB in the UI, in words, on every surface.
//
// WHAT IT DOES HONESTLY
//   * lexes and parses TPL exactly as PROTOCOL_LANGUAGE.md v1 describes;
//   * enforces mandatory units, block closure, declared-before-used, and the
//     max_runs refusal, all with line and column;
//   * expands a sweep cross product and states the run count and a cost estimate;
//   * converts dose -> PKConfig and WRITES OUT the conversion assumptions,
//     because the spec says an unstated conversion is not acceptable;
//   * executes by expanding each planned run into a v1 RunSpec and submitting it
//     through the real POST /api/runs, so a sweep cell opens a real run in the
//     viewer. It does not simulate anything itself.
//
// WHAT IT REFUSES TO DO
//   * compile a structure (`from smiles` / `from file`) — that is
//     `tmesim.molecular`, it is not reachable from a browser, and guessing a
//     surrogate compound would be inventing a number. It errors and says so.
//   * invent a DrugParams field it cannot derive from an uploaded file. Such a
//     field is reported with provenance `flagged-default` and a null value.

import { getSections, getDrugs, postRun, getStatus, getRun, cancelRun } from '../api.js';
import { tokenizeLine } from '../editor/tpl.js';

/* ------------------------------------------------------------- tables --- */

const TIME_S = { s: 1, min: 60, h: 3600, d: 86400 };
const LEN_UM = { um: 1, mm: 1000 };
const CONC_UM = { uM: 1, nM: 1e-3, mM: 1e3 };

const NEED_UNIT = {
  ic50: ['uM', 'nM', 'mM'],
  diffusivity: ['um2/s'],
  dose: ['mg/m2', 'mg/kg', 'mg'],
  over: ['s', 'min', 'h', 'd'],
  every: ['s', 'min', 'h', 'd'],
  rest: ['s', 'min', 'h', 'd'],
  dx: ['um', 'mm'],
};
const BARE_OK = new Set(['hill', 'binding', 'clearance', 'bbb_penetrance',
  'molecular_weight', 'repeat', 'cycles', 'runs']);

const DRUG_FIELDS = new Set(['ic50', 'hill', 'diffusivity', 'binding', 'clearance',
  'bbb_penetrance', 'molecular_weight']);
const PROTO_FIELDS = new Set(['route', 'dose', 'every', 'repeat', 'then', 'rest', 'cycles']);
const RUN_FIELDS = new Set(['dx', 'state_mode', 'sensitivity', 'report']);
// `rank` is NOT a run field any more. It is still recognised in runOpt() so the
// refusal names the withdrawal instead of reading as a typo. Kept identical to
// tmesim/protocol/parser.py:WITHDRAWN_RANK_BY — an editor that accepts what the
// engine refuses (or the reverse) is worse than either rule on its own.
const WITHDRAWN_RANK_BY = '`rank by` has been WITHDRAWN from TPL (2026-08-21). Ordering compounds by predicted effect is no longer an output of this repo, because the ordering was measured and does not hold: results/pk_ranking_robustness.md — per-compound plasma exposure spans 113.1x across the compounds the panel ranked, and the published GBM order first breaks at a spread of only 1.53x; results/capacity_ranking_robustness.md — the breast top-1 pick changes under this repo\'s own structural capacity rule. Report the metric with `report <metric>` and read the map for ONE compound; use `compare` for two or more RUNS of one compound (a dose ladder), not for two compounds.';
const WITHDRAWN_CROSS_COMPOUND = 'ordering two different compounds by predicted effect has been WITHDRAWN (2026-08-21). It was measured and it does not hold: results/pk_ranking_robustness.md — per-compound plasma exposure spans 113.1x across the compounds the panel ranked, while the published GBM order first breaks at a spread of only 1.53x; results/capacity_ranking_robustness.md — the breast top-1 pick changes under this repo\'s own structural capacity rule. `compare` now takes two or more RUNS OF ONE COMPOUND on one section — a dose ladder, or any other `sweep`. Run each compound in its own program if you need both maps; the map for one compound is the product, the ordering between them is not.';
const ROUTES = ['oral', 'iv', 'ip'];
const STATE_MODES = ['oxygen', 'vessel_distance', 'markers'];
const SENSITIVITIES = ['uniform', 'determinant'];
const SWEEPABLE = new Set(['dose', 'every', 'repeat', 'cycles', 'dx', 'over', 'rest']);

export const METRICS = {
  coverage: { spec: 'therapeutic_coverage_1x', label: 'Coverage @ 1x IC50', unit: '%', fmt: 'pct', better: 'higher' },
  kill: { spec: 'kill_fraction', label: 'Kill fraction', unit: '%', fmt: 'pct', better: 'higher' },
  mean_auc: { spec: 'mean_tumor_auc', label: 'Mean tumour AUC', unit: 'uM.s', fmt: 'sci', better: 'higher' },
  penetration_length: { spec: 'penetration_length_um', label: 'Penetration length', unit: 'um', fmt: 'num1', better: 'higher' },
};

const DEFAULT_MAX_RUNS = 64;

/* --------------------------------------------------------- catalogue ---- */

let catalogue = null;
async function loadCatalogue() {
  if (catalogue) return catalogue;
  const [sections, drugs] = await Promise.all([getSections(), getDrugs()]);
  catalogue = { sections: sections || [], drugs: drugs || [] };
  return catalogue;
}

/** What a specimen's panel can back, per API_CONTRACT + results/cell_state_impact. */
function capability(sec) {
  const panel = (sec && sec.panel) || '';
  const deconvolved = /deconv|anatomic|ISH/i.test(panel);
  return {
    markers: !deconvolved,
    markers_complete: false,      // no panel in this deployment carries a hypoxia transcript
    note: deconvolved
      ? `${sec.id}: the panel is "${panel}" — there is no per-cell transcript to band a cell state from.`
      : `${sec.id}: marker banding is available but the panel carries no hypoxia-response transcript, so HYPOXIC is UNMEASURED.`,
  };
}

/* ------------------------------------------------------------- parser --- */

function sigOf(line) {
  return tokenizeLine(line).filter((t) => t.t !== 'ws' && t.t !== 'comment');
}

function parse(source) {
  const lines = source.split('\n');
  const errors = [];
  const warnings = [];
  const prog = {
    drugs: new Map(), protocols: new Map(), sweeps: [],
    statements: [], expects: [], allow: null,
  };

  const D = (bag, line, tok, message, extra) => bag.push({
    line,
    col: tok ? tok.s : 0,
    end_col: tok ? tok.e : undefined,
    severity: bag === errors ? 'error' : 'warning',
    source: 'stub',
    message,
    ...(extra || {}),
  });
  const E = (line, tok, message, extra) => D(errors, line, tok, message, extra);
  const W = (line, tok, message, extra) => D(warnings, line, tok, message, extra);

  /** Read NUM [UNIT] at sig[i]; enforces the field's unit requirement. */
  function qty(sig, i, field, line) {
    const n = sig[i];
    if (!n || n.t !== 'number') {
      E(line, n || sig[i - 1], `\`${field}\` expects a number here.`);
      return null;
    }
    const value = parseFloat(n.v);
    const u = sig[i + 1];
    const need = NEED_UNIT[field];
    if (need) {
      if (!u || u.t !== 'unit') {
        E(line, n, `\`${field} ${n.v}\` has no unit. Units are MANDATORY in TPL — write \`${field} ${n.v} ${need[0]}\`. Recognised here: ${need.join(', ')}.`,
          {
            code: 'missing_unit',
            remedy: `${field} ${n.v} ${need[0]}`,
            remedy_text: `${field} ${n.v} ${need[0]}`,
            remedy_kind: 'replace_line',
            remedy_line: line,
          });
        return { value, unit: null, next: i + 1, tok: n };
      }
      if (!need.includes(u.v)) {
        E(line, u, `\`${u.v}\` is not a valid unit for \`${field}\`. Expected one of: ${need.join(', ')}.`,
          { code: 'bad_unit' });
        return { value, unit: null, next: i + 2, tok: n };
      }
      return { value, unit: u.v, next: i + 2, tok: n };
    }
    if (u && u.t === 'unit') {
      E(line, u, `\`${field}\` is a bare quantity and takes no unit; \`${u.v}\` is not meaningful here.`);
      return { value, unit: null, next: i + 2, tok: n };
    }
    return { value, unit: null, next: i + 1, tok: n };
  }

  let block = null;

  lines.forEach((raw, idx) => {
    const no = idx + 1;
    const sig = sigOf(raw);
    if (!sig.length) return;
    const head = sig[0];

    /* ---- block close ------------------------------------------------- */
    if (head.t === 'kw' && head.v === 'end') {
      if (!block) {
        E(no, head, '`end` closes nothing — there is no open block above it.');
      } else {
        if (block.kind === 'protocol' && !block.obj.dose) {
          E(block.line, null, `protocol \`${block.obj.name}\` declares no \`dose\`. A protocol without a dose cannot be converted to a PK profile.`);
        }
        if (block.kind === 'drug' && !block.obj.ic50 && !block.obj.ic50FromGdsc) {
          W(block.line, null, `drug \`${block.obj.name}\` gives no \`ic50\`. It will be DERIVED by the molecular layer and labelled \`computed\`, not measured.`);
        }
        block = null;
      }
      if (sig.length > 1) E(no, sig[1], '`end` takes nothing after it.');
      return;
    }

    /* ---- inside a block ---------------------------------------------- */
    if (block) {
      parseBlockLine(block, sig, no);
      return;
    }

    /* ---- top level ---------------------------------------------------- */
    switch (head.t === 'kw' ? head.v : '') {
      case 'drug': return openDrug(sig, no);
      case 'protocol': return openProtocol(sig, no);
      case 'run': return openRun(sig, no, 'run');
      case 'compare': return openRun(sig, no, 'compare');
      case 'sweep': return parseSweep(sig, no);
      case 'expect': return parseExpect(sig, no);
      case 'allow': return parseAllow(sig, no);
      default:
        E(no, head, `\`${head.v}\` is not a TPL statement. A program is made of drug / protocol / sweep / run / compare / expect / allow statements.`);
    }
  });

  if (block) {
    E(block.line, null, `\`${block.kind}\` block is never closed — add \`end\` on its own line.`,
      { code: 'unclosed_block' });
  }

  /* ------------------------------------------------------ statement fns -- */

  function openDrug(sig, no) {
    const name = sig[1];
    if (!name || name.t !== 'ident') { E(no, name || sig[0], '`drug` needs a name: `drug TMZ from library "temozolomide"`.'); return; }
    if (!sig[2] || sig[2].v !== 'from') { E(no, sig[2] || name, 'expected `from` after the drug name.'); return; }
    const kind = sig[3];
    if (!kind || !['smiles', 'library', 'file'].includes(kind.v)) {
      E(no, kind || sig[2], 'a drug comes `from smiles "..."`, `from library "..."`, or `from file "..."`.');
      return;
    }
    const str = sig[4];
    if (!str || str.t !== 'string') { E(no, str || kind, `\`from ${kind.v}\` needs a quoted string.`); return; }
    if (prog.drugs.has(name.v)) E(no, name, `drug \`${name.v}\` is already defined.`);
    const obj = {
      name: name.v, line: no, sourceKind: kind.v, sourceValue: str.v.slice(1, -1),
      fields: {}, ic50: null, ic50FromGdsc: false, srcTok: kind,
    };
    if (kind.v !== 'library') {
      E(no, kind, `STUB LIMIT: compiling a drug \`from ${kind.v}\` runs through tmesim.molecular, which this browser stub cannot reach. This is a limitation of STUB MODE, not of TPL — the real executor does it. Define the drug \`from library "..."\` to work offline.`,
        { code: 'stub_unsupported' });
    }
    prog.drugs.set(name.v, obj);
    block = { kind: 'drug', obj, line: no };
  }

  function openProtocol(sig, no) {
    const name = sig[1];
    if (!name || name.t !== 'ident') { E(no, name || sig[0], '`protocol` needs a name: `protocol stupp_adjuvant`.'); return; }
    if (sig.length > 2) E(no, sig[2], '`protocol NAME` takes nothing else on its line.');
    if (prog.protocols.has(name.v)) E(no, name, `protocol \`${name.v}\` is already defined.`);
    const obj = { name: name.v, line: no, route: null, dose: null, over: null, every: null, repeat: null, rest: null, cycles: null };
    prog.protocols.set(name.v, obj);
    block = { kind: 'protocol', obj, line: no };
  }

  function openRun(sig, no, kind) {
    const drugs = [];
    let i = 1;
    while (sig[i] && sig[i].t === 'ident') {
      drugs.push({ name: sig[i].v, tok: sig[i] });
      i += 1;
      if (sig[i] && sig[i].v === ',') { i += 1; continue; }
      break;
    }
    if (!drugs.length) { E(no, sig[1] || sig[0], `\`${kind}\` needs at least one drug name.`); return; }
    // Arity is checked in semantics(), not here: whether a `compare` plans two
    // or more runs depends on the sweeps, which are a whole-program property.
    if (kind === 'run' && drugs.length > 1) E(no, drugs[1].tok, '`run` takes ONE compound. ' + WITHDRAWN_CROSS_COMPOUND);
    if (!sig[i] || sig[i].v !== 'with') { E(no, sig[i] || sig[0], 'expected `with <protocol>`.'); return; }
    const proto = sig[i + 1];
    if (!proto || proto.t !== 'ident') { E(no, proto || sig[i], '`with` needs a protocol name.'); return; }
    if (!sig[i + 2] || sig[i + 2].v !== 'on') { E(no, sig[i + 2] || proto, 'expected `on <section>`.'); return; }
    const sec = sig[i + 3];
    if (!sec || sec.t !== 'ident') { E(no, sec || sig[i + 2], '`on` needs a section id.'); return; }
    const obj = {
      kind, line: no, drugs, protocol: proto.v, protoTok: proto,
      section: sec.v, sectionTok: sec,
      opts: { dx: null, state_mode: null, sensitivity: null, report: [] },
    };
    prog.statements.push(obj);
    block = { kind, obj, line: no };
  }

  function parseSweep(sig, no) {
    const p = sig[1];
    if (!p || !(p.t === 'fld' || p.t === 'ident')) { E(no, p || sig[0], '`sweep` needs a parameter name.'); return; }
    if (!SWEEPABLE.has(p.v)) {
      E(no, p, `\`${p.v}\` is not a sweepable parameter. Sweepable: ${[...SWEEPABLE].join(', ')}.`);
      return;
    }
    // fixed override: `sweep dx 12 um`
    if (sig[2] && sig[2].t === 'number') {
      const q = qty(sig, 2, p.v, no);
      if (q) prog.sweeps.push({ kind: 'fixed', param: p.v, line: no, value: q.value, unit: q.unit });
      return;
    }
    if (!sig[2] || sig[2].v !== 'from') { E(no, sig[2] || p, 'a sweep is `sweep P from A to B UNIT step S`, `... in N steps`, or a fixed `sweep P V UNIT`.'); return; }
    const a = sig[3];
    if (!a || a.t !== 'number') { E(no, a || sig[2], '`from` needs a number.'); return; }
    if (!sig[4] || sig[4].v !== 'to') { E(no, sig[4] || a, 'expected `to`.'); return; }
    const b = sig[5];
    if (!b || b.t !== 'number') { E(no, b || sig[4], '`to` needs a number.'); return; }
    const need = NEED_UNIT[p.v];
    let unit = null;
    let i = 6;
    if (need) {
      const u = sig[6];
      if (!u || u.t !== 'unit') {
        E(no, b, `sweep of \`${p.v}\` has no unit. Write \`sweep ${p.v} from ${a.v} to ${b.v} ${need[0]} ...\`.`,
          { code: 'missing_unit' });
      } else if (!need.includes(u.v)) {
        E(no, u, `\`${u.v}\` is not a valid unit for \`${p.v}\`. Expected: ${need.join(', ')}.`);
      } else { unit = u.v; i = 7; }
    }
    const mode = sig[i];
    if (mode && mode.v === 'step') {
      const s = sig[i + 1];
      if (!s || s.t !== 'number') { E(no, s || mode, '`step` needs a number.'); return; }
      const step = parseFloat(s.v);
      if (!(step > 0)) { E(no, s, '`step` must be > 0.'); return; }
      prog.sweeps.push({ kind: 'range', param: p.v, line: no, from: parseFloat(a.v), to: parseFloat(b.v), unit, step, tok: p });
      return;
    }
    if (mode && mode.v === 'in') {
      const s = sig[i + 1];
      if (!s || s.t !== 'number') { E(no, s || mode, '`in` needs a step count.'); return; }
      const n = parseInt(s.v, 10);
      if (!(n >= 2)) { E(no, s, '`in N steps` needs N >= 2.'); return; }
      if (!sig[i + 2] || sig[i + 2].v !== 'steps') E(no, sig[i + 2] || s, 'expected `steps` after the count.');
      prog.sweeps.push({ kind: 'range', param: p.v, line: no, from: parseFloat(a.v), to: parseFloat(b.v), unit, count: n, tok: p });
      return;
    }
    E(no, mode || b, 'a range sweep ends with `step S` or `in N steps`.');
  }

  function parseExpect(sig, no) {
    const m = sig[1];
    if (!m || !METRICS[m.v]) {
      E(no, m || sig[0], `\`expect\` needs a reported metric. Known: ${Object.keys(METRICS).join(', ')}.`);
      return;
    }
    const op = sig[2];
    if (!op || op.t !== 'op' || !['>', '<', '>=', '<=', '==', '!='].includes(op.v)) {
      E(no, op || m, '`expect` needs a comparator: > < >= <= == !=.');
      return;
    }
    const n = sig[3];
    if (!n || n.t !== 'number') { E(no, n || op, '`expect` needs a number to compare against.'); return; }
    prog.expects.push({ metric: m.v, op: op.v, value: parseFloat(n.v), unit: sig[4] && sig[4].t === 'unit' ? sig[4].v : null, line: no });
  }

  function parseAllow(sig, no) {
    if (!sig[1] || sig[1].v !== 'runs') { E(no, sig[1] || sig[0], 'the only `allow` form is `allow runs N`.'); return; }
    const n = sig[2];
    if (!n || n.t !== 'number') { E(no, n || sig[1], '`allow runs` needs an integer.'); return; }
    prog.allow = { n: parseInt(n.v, 10), line: no };
  }

  /* ------------------------------------------------------- block bodies -- */

  function parseBlockLine(blk, sig, no) {
    const head = sig[0];
    if (blk.kind === 'drug') return drugField(blk.obj, sig, head, no);
    if (blk.kind === 'protocol') return protoField(blk.obj, sig, head, no);
    return runOpt(blk.obj, sig, head, no);
  }

  function drugField(obj, sig, head, no) {
    if (!DRUG_FIELDS.has(head.v)) {
      E(no, head, `\`${head.v}\` is not an overridable DrugParams field. Overridable: ${[...DRUG_FIELDS].join(', ')}.`);
      return;
    }
    if (head.v === 'ic50' && sig[1] && sig[1].v === 'from') {
      const src = sig[2];
      if (!src || src.v !== 'gdsc') { E(no, src || sig[1], 'the only measured IC50 source is `gdsc`.'); return; }
      obj.ic50FromGdsc = true;
      obj.ic50Line = no;
      obj.ic50Tok = src;
      return;
    }
    const q = qty(sig, 1, head.v, no);
    if (!q) return;
    obj.fields[head.v] = { value: q.value, unit: q.unit, line: no };
    if (head.v === 'ic50') obj.ic50 = obj.fields.ic50;
  }

  function protoField(obj, sig, head, no) {
    if (head.v === 'then') {
      if (!sig[1] || sig[1].v !== 'rest') { E(no, sig[1] || head, 'the only `then` form is `then rest N d`.'); return; }
      const q = qty(sig, 2, 'rest', no);
      if (q) obj.rest = { value: q.value, unit: q.unit, line: no };
      return;
    }
    if (!PROTO_FIELDS.has(head.v)) {
      E(no, head, `\`${head.v}\` is not a protocol statement. Allowed: route, dose, every, repeat, then rest, cycles.`);
      return;
    }
    if (head.v === 'route') {
      const r = sig[1];
      if (!r || !ROUTES.includes(r.v)) { E(no, r || head, `route must be one of ${ROUTES.join(' | ')}.`); return; }
      obj.route = r.v;
      return;
    }
    if (head.v === 'dose') {
      const q = qty(sig, 1, 'dose', no);
      if (!q) return;
      obj.dose = { value: q.value, unit: q.unit, line: no };
      const nx = sig[q.next];
      if (nx && nx.v === 'over') {
        const o = qty(sig, q.next + 1, 'over', no);
        if (o) obj.over = { value: o.value, unit: o.unit, line: no };
      } else if (nx) {
        E(no, nx, `unexpected \`${nx.v}\` after the dose. The only continuation is \`over N h\` (an infusion duration).`);
      }
      return;
    }
    const q = qty(sig, 1, head.v, no);
    if (!q) return;
    if (head.v === 'every') obj.every = { value: q.value, unit: q.unit, line: no };
    if (head.v === 'repeat') obj.repeat = { value: Math.round(q.value), line: no };
    if (head.v === 'cycles') obj.cycles = { value: Math.round(q.value), line: no };
    if (head.v === 'rest') obj.rest = { value: q.value, unit: q.unit, line: no };
  }

  function runOpt(obj, sig, head, no) {
    if (head.v === 'rank') { E(no, head, WITHDRAWN_RANK_BY, { code: 'withdrawn_ranking' }); return; }
    if (!RUN_FIELDS.has(head.v)) {
      E(no, head, `\`${head.v}\` is not a run option. Allowed: dx, state_mode, sensitivity, report.`);
      return;
    }
    if (head.v === 'dx') {
      const q = qty(sig, 1, 'dx', no);
      if (q) obj.opts.dx = { value: q.value, unit: q.unit, line: no };
      return;
    }
    if (head.v === 'state_mode') {
      const v = sig[1];
      if (!v || !STATE_MODES.includes(v.v)) { E(no, v || head, `state_mode must be one of ${STATE_MODES.join(' | ')}.`); return; }
      obj.opts.state_mode = { value: v.v, line: no, tok: v };
      return;
    }
    if (head.v === 'sensitivity') {
      const v = sig[1];
      if (!v || !SENSITIVITIES.includes(v.v)) { E(no, v || head, `sensitivity must be one of ${SENSITIVITIES.join(' | ')}.`); return; }
      obj.opts.sensitivity = { value: v.v, line: no, tok: v };
      return;
    }
    // report a, b, c
    let i = 1;
    while (sig[i]) {
      const m = sig[i];
      if (m.v === ',') { i += 1; continue; }
      if (!METRICS[m.v]) { E(no, m, `\`${m.v}\` is not a reportable metric. Known: ${Object.keys(METRICS).join(', ')}.`); return; }
      obj.opts.report.push(m.v);
      i += 1;
    }
    if (!obj.opts.report.length) E(no, head, '`report` needs at least one metric.');
  }

  return { prog, errors, warnings };
}

/* ------------------------------------------------------ semantic pass --- */

function semantics(prog, cat, errors, warnings) {
  const E = (line, tok, message, extra) => errors.push({
    line, col: tok ? tok.s : 0, end_col: tok ? tok.e : undefined,
    severity: 'error', source: 'stub', message, ...(extra || {}),
  });
  const W = (line, tok, message, extra) => warnings.push({
    line, col: tok ? tok.s : 0, end_col: tok ? tok.e : undefined,
    severity: 'warning', source: 'stub', message, ...(extra || {}),
  });

  for (const d of prog.drugs.values()) {
    if (d.sourceKind !== 'library') continue;
    const known = cat.drugs.find((x) => x.name === d.sourceValue);
    if (!known) {
      E(d.line, d.srcTok, `\`${d.sourceValue}\` is not in the built-in library. Known: ${cat.drugs.map((x) => x.name).join(', ')}.`,
        { code: 'unknown_drug' });
      continue;
    }
    d.resolved = known;
    if (d.ic50FromGdsc && !/measured/i.test(known.ic50_source || '')) {
      E(d.ic50Line, d.ic50Tok, `\`ic50 from gdsc\` refuses: ${d.sourceValue} has no measured GDSC entry (source is "${known.ic50_source}"). Give an explicit \`ic50 <n> uM\` instead — it will be recorded as user-supplied.`,
        { code: 'ic50_not_measured' });
    }
  }

  for (const s of prog.statements) {
    for (const d of s.drugs) {
      if (!prog.drugs.has(d.name)) {
        E(s.line, d.tok, `unknown drug \`${d.name}\` — declare it with a \`drug ${d.name} from ...\` block above.`,
          { code: 'undeclared_drug' });
      }
    }
    if (!prog.protocols.has(s.protocol)) {
      E(s.line, s.protoTok, `unknown protocol \`${s.protocol}\` — declare it with a \`protocol ${s.protocol}\` block above.`,
        { code: 'undeclared_protocol' });
    }
    const sec = cat.sections.find((x) => x.id === s.section);
    if (!sec) {
      E(s.line, s.sectionTok, `unknown section \`${s.section}\`. Available: ${cat.sections.filter((x) => x.available).map((x) => x.id).join(', ')}.`,
        { code: 'unknown_section' });
    } else if (!sec.available) {
      E(s.line, s.sectionTok, `section \`${s.section}\` is not available in this deployment (controlled access).`,
        { code: 'section_unavailable' });
    } else {
      const cap = capability(sec);
      const sm = s.opts.state_mode;
      if (sm && sm.value === 'markers') {
        if (!cap.markers) {
          E(sm.line, sm.tok, `\`state_mode markers\` is unsupported on \`${s.section}\`: ${cap.note}`, { code: 'unsupported_mode' });
        } else {
          W(sm.line, sm.tok, `\`state_mode markers\` on \`${s.section}\`: ${cap.note} An unmeasured state is assigned to NO cell — that is not "computed and found to be zero".`);
        }
      }
      if (s.opts.sensitivity && s.opts.sensitivity.value === 'uniform') {
        W(s.opts.sensitivity.line, s.opts.sensitivity.tok,
          '`sensitivity uniform` sets every cell to 1.0. Kill and coverage from this run are a NO-RESISTANCE UPPER BOUND, not an expected response.');
      }
      if (!s.opts.sensitivity) {
        W(s.line, null, 'no `sensitivity` given; the executor defaults to `uniform`, which makes kill an upper bound. State it explicitly.');
      }
    }
    if (!s.opts.report.length) {
      W(s.line, null, 'no `report` list; the executor will report its default metric set.');
    }
    if (s.kind === 'compare') {
      // Mirrors tmesim/protocol/executor.py:_check_compare, in the same order.
      const distinct = [...new Set(s.drugs.map((d) => d.name))];
      if (distinct.length > 1) {
        E(s.line, s.drugs[1].tok, `compare ${distinct.join(', ')}: ${WITHDRAWN_CROSS_COMPOUND}`,
          { code: 'withdrawn_ranking' });
      } else if (s.drugs.length > 1) {
        E(s.line, s.drugs[1].tok, `compare ${distinct[0]}, ${distinct[0]}: ${distinct[0]} is named ${s.drugs.length} times, so every run of this compare is the same run. A compare over one compound varies a swept PARAMETER, not the compound — write \`compare ${distinct[0]} with ...\` and add a \`sweep\`.`,
          { code: 'incomplete_compare' });
      } else if (plannedRuns(prog) < 2) {
        E(s.line, null, `compare ${distinct[0]}: this program plans only one run, and one run is not a comparison. Add a \`sweep\` so there are two or more runs of ${distinct[0]} to put side by side (e.g. \`sweep dose from 50 to 400 mg/m2 step 50\`), or use \`run\` for a single run.`,
          { code: 'incomplete_compare' });
      }
    }
  }

  for (const sw of prog.sweeps) {
    if (sw.kind !== 'range') continue;
    if (sw.to <= sw.from) E(sw.line, sw.tok, `sweep range is empty: to (${sw.to}) must exceed from (${sw.from}).`);
    const proto = [...prog.protocols.values()][0];
    if (sw.param === 'dose' && !proto) W(sw.line, sw.tok, 'sweeping `dose` with no protocol declared has nothing to override.');
  }
  if (!prog.statements.length) {
    warnings.push({ line: 1, severity: 'warning', source: 'stub', message: 'the program declares no `run` or `compare`, so it plans zero runs.' });
  }
}

/* --------------------------------------------------------- expansion ---- */

/**
 * Runs ONE statement expands to: the cross-product of the range sweeps. A fixed
 * sweep is one point, so it does not turn a single run into a comparison.
 * Used by semantics() to tell "compare over a dose ladder" from "compare that
 * plans exactly one run".
 */
const plannedRuns = (prog) => prog.sweeps.reduce(
  (n, sw) => n * (sw.kind === 'range' ? Math.max(1, axisValues(sw).length) : 1), 1);

const axisValues = (sw) => {
  const out = [];
  if (sw.count) {
    for (let i = 0; i < sw.count; i++) out.push(sw.from + (sw.to - sw.from) * (i / (sw.count - 1)));
    return out;
  }
  for (let v = sw.from; v <= sw.to + 1e-9; v += sw.step) out.push(Math.round(v * 1e9) / 1e9);
  return out;
};

const BSA_M2 = 1.73;
const MASS_KG = 70;
const VD_L = 30;
const HALF_LIFE_S = 14400;
const TIME_BUDGET_S = 2.5e6;

function toRunSpec(params, stmt, drugDef, sec) {
  const assumptions = [];
  const mw = drugDef && drugDef.resolved ? drugDef.resolved.mw : null;
  const doseMg = params.dose
    ? (params.dose.unit === 'mg' ? params.dose.value
      : params.dose.unit === 'mg/kg' ? params.dose.value * MASS_KG
        : params.dose.value * BSA_M2)
    : null;

  if (params.dose && params.dose.unit === 'mg/m2') {
    assumptions.push(`body-surface area = ${BSA_M2} m2 (adult mean) — ${params.dose.value} mg/m2 -> ${doseMg.toFixed(1)} mg`);
  }
  if (params.dose && params.dose.unit === 'mg/kg') {
    assumptions.push(`body mass = ${MASS_KG} kg — ${params.dose.value} mg/kg -> ${doseMg.toFixed(1)} mg`);
  }
  assumptions.push(`volume of distribution = ${VD_L} L, bioavailability F = 1.0 for route \`${params.route}\` — FLAGGED DEFAULT, not a measured PK parameter for this compound`);
  assumptions.push(`molecular weight = ${mw === null ? 'UNKNOWN' : mw + ' g/mol'} (${mw === null ? 'unresolved' : 'catalogue, measured'}) for the mg -> uM conversion`);
  assumptions.push(`plasma half-life = ${HALF_LIFE_S} s — FLAGGED DEFAULT; the catalogue carries no per-drug half-life`);

  const peak = doseMg !== null && mw ? (doseMg / VD_L) / mw * 1000 : 1.0;
  if (doseMg === null || !mw) assumptions.push('plasma_peak fell back to 1.0 uM because the dose or the molecular weight could not be resolved — FLAGGED DEFAULT');

  const everyS = params.every ? params.every.value * TIME_S[params.every.unit] : 86400;
  const repeat = params.repeat || 1;
  const cycles = params.cycles || 1;
  const nDoses = Math.max(1, Math.round(repeat * cycles));
  let halfLife = HALF_LIFE_S;
  if (params.over) {
    halfLife = Math.max(HALF_LIFE_S, params.over.value * TIME_S[params.over.unit]);
    assumptions.push(`\`over ${params.over.value} ${params.over.unit}\` is an INFUSION. PROTOCOL_LANGUAGE requires PDEConfig.time_resolved = true for it; RunSpec v1 carries no such field, so the STUB approximates the sustained profile by widening half_life_s to ${halfLife} s. The real executor sets time_resolved and says so.`);
  }
  const courseS = nDoses * everyS;
  const totalS = Math.min(TIME_BUDGET_S, Math.max(600, courseS));
  if (courseS > TIME_BUDGET_S) {
    assumptions.push(`the schedule spans ${(courseS / 86400).toFixed(1)} d but the engine's integration budget is ${(TIME_BUDGET_S / 86400).toFixed(0)} d — total_time_s is CAPPED, so the reported AUC covers the first ${(totalS / 86400).toFixed(0)} d only.`);
  }
  if (params.rest && cycles > 1) {
    assumptions.push(`PKConfig expresses a UNIFORM dose interval. The ${repeat}-on / ${params.rest.value} ${params.rest.unit}-rest x ${cycles}-cycle structure is FLATTENED to n_doses=${nDoses} at ${everyS} s. The real executor honours the rest periods; this stub does not.`);
  }

  const dxUm = params.dx ? params.dx.value * LEN_UM[params.dx.unit || 'um'] : 12;
  const spec = {
    section: stmt.section,
    dx_um: dxUm,
    drug: drugDef && drugDef.resolved ? drugDef.resolved.name : null,
    sensitivity_mode: (stmt.opts.sensitivity && stmt.opts.sensitivity.value) || 'uniform',
    vascular_mode: 'discrete',
    oxygen_mode: 'pde',
    state_mode: (stmt.opts.state_mode && stmt.opts.state_mode.value) || 'oxygen',
    mgmt: null,
    pk: {
      plasma_peak: Math.round(peak * 1e6) / 1e6,
      half_life_s: halfLife,
      n_doses: nDoses,
      dose_interval_s: everyS,
    },
    total_time_s: totalS,
  };
  const estS = sec
    ? 6 + (sec.n_cells / 40000) * Math.pow(12 / Math.max(2, dxUm), 2) * 11 * (nDoses > 1 ? 1.35 : 1)
    : 20;
  return { spec, assumptions, est_s: Math.round(estS * 10) / 10 };
}

function fmtQty(q) {
  if (!q) return '—';
  const v = Math.abs(q.value) >= 100 ? q.value.toFixed(0) : q.value.toFixed(Math.abs(q.value) < 1 ? 3 : 1);
  return q.unit ? `${v} ${q.unit}` : v;
}

function expand(prog, cat) {
  const axes = prog.sweeps.filter((s) => s.kind === 'range').map((s) => ({
    name: s.param, unit: s.unit, line: s.line, values: axisValues(s),
  }));
  const fixed = prog.sweeps.filter((s) => s.kind === 'fixed');
  const combos = axes.reduce((acc, ax) => {
    const out = [];
    for (const base of acc) for (let i = 0; i < ax.values.length; i++) out.push(base.concat(i));
    return out;
  }, [[]]);

  const runs = [];
  const assumptionSet = new Set();
  for (const stmt of prog.statements) {
    const proto = prog.protocols.get(stmt.protocol);
    if (!proto) continue;
    const sec = cat.sections.find((x) => x.id === stmt.section);
    for (const d of stmt.drugs) {
      const drugDef = prog.drugs.get(d.name);
      for (const coords of combos) {
        const params = {
          route: proto.route || 'oral',
          dose: proto.dose ? { ...proto.dose } : null,
          over: proto.over ? { ...proto.over } : null,
          every: proto.every ? { ...proto.every } : null,
          repeat: proto.repeat ? proto.repeat.value : 1,
          cycles: proto.cycles ? proto.cycles.value : 1,
          rest: proto.rest ? { ...proto.rest } : null,
          dx: stmt.opts.dx ? { ...stmt.opts.dx } : null,
        };
        for (const f of fixed) applyParam(params, f.param, f.value, f.unit);
        coords.forEach((vi, ai) => applyParam(params, axes[ai].name, axes[ai].values[vi], axes[ai].unit));
        if (!params.dx) params.dx = { value: 12, unit: 'um', flagged: true };

        const { spec, assumptions, est_s } = toRunSpec(params, stmt, drugDef, sec);
        assumptions.forEach((a) => assumptionSet.add(a));
        const label = [
          stmt.drugs.length > 1 || prog.statements.length > 1 ? d.name : null,
          ...coords.map((vi, ai) => `${axes[ai].name} ${fmtQty({ value: axes[ai].values[vi], unit: axes[ai].unit })}`),
        ].filter(Boolean).join(' · ') || `${d.name} · ${stmt.protocol}`;

        runs.push({
          index: runs.length,
          coords,
          drug: d.name,
          statement_line: stmt.line,
          label,
          params: {
            dose: params.dose, every: params.every, dx: params.dx, over: params.over,
            repeat: params.repeat, cycles: params.cycles, route: params.route,
          },
          spec,
          est_s,
          report: stmt.opts.report.length ? stmt.opts.report : Object.keys(METRICS),
        });
      }
    }
  }
  return { axes, runs, assumptions: [...assumptionSet] };
}

function applyParam(params, name, value, unit) {
  switch (name) {
    case 'dose': params.dose = { value, unit: unit || (params.dose && params.dose.unit) || 'mg/m2' }; break;
    case 'every': params.every = { value, unit: unit || 'h' }; break;
    case 'over': params.over = { value, unit: unit || 'h' }; break;
    case 'rest': params.rest = { value, unit: unit || 'd' }; break;
    case 'dx': params.dx = { value, unit: unit || 'um' }; break;
    case 'repeat': params.repeat = Math.round(value); break;
    case 'cycles': params.cycles = Math.round(value); break;
    default: break;
  }
}

/* ---------------------------------------------------------- caveats ----- */

function standingCaveats(prog, runs) {
  const out = [
    'STUB EXECUTOR. tmesim/protocol/ is not reachable from this browser. The program below was parsed and expanded locally and each run was submitted as a v1 RunSpec through POST /api/runs. The dose -> PK arithmetic is the stub\'s, not the executor\'s.',
  ];
  const anyUniform = runs.some((r) => r.spec.sensitivity_mode === 'uniform');
  if (anyUniform) {
    out.push('sensitivity uniform sets every cell to 1.0: kill and coverage are a NO-RESISTANCE UPPER BOUND, not an expected response. Differences BETWEEN cells of this sweep are differences between upper bounds.');
  }
  if (runs.some((r) => r.spec.section === 'xenium_319')) {
    out.push('ABSOLUTE KILL IS NOT DEFENSIBLE ON xenium_319 UNDER ANY STATE ROUTE: the 319-gene panel carries no hypoxia-response transcript, so no cell-state route clears the bar and kill moves by up to 0.3252 absolute (1.77x) across routes. What this grid supports is the RELATIVE reading within ONE compound on this one section — how the map changes from one point of the swept axis to the next, and how one region of a map stands against another. The absolute level does not, and neither does any comparison BETWEEN compounds: that was withdrawn on 2026-08-21 (results/pk_ranking_robustness.md, results/capacity_ranking_robustness.md).');
  }
  out.push('No section in this deployment carries a measured treatment-outcome label. Nothing on this grid has been scored against a real drug response.');
  return out;
}

/* ------------------------------------------------------- public: validate */

export async function validate(source) {
  const cat = await loadCatalogue();
  const { prog, errors, warnings } = parse(source);
  semantics(prog, cat, errors, warnings);
  return {
    ok: errors.length === 0,
    errors,
    warnings,
    stub: true,
    symbols: {
      drugs: [...prog.drugs.keys()],
      protocols: [...prog.protocols.keys()],
      sweeps: prog.sweeps.map((s) => s.param),
    },
  };
}

/* ----------------------------------------------------------- public: plan */

export async function plan(source) {
  const cat = await loadCatalogue();
  const { prog, errors, warnings } = parse(source);
  semantics(prog, cat, errors, warnings);
  const maxRuns = prog.allow ? prog.allow.n : DEFAULT_MAX_RUNS;
  if (errors.length) {
    // A plan that cannot be made still answers with the shape of a plan, so the
    // panel renders "0 runs, of 64 max_runs" rather than a row of NaNs.
    return {
      ok: false, stub: true, errors, warnings,
      runs: [], axes: [], metrics: [], assumptions: [],
      run_count: 0, max_runs: maxRuns, estimated_cost_s: 0,
      cost_note: null, caveats: [], refusal: null,
    };
  }
  const { axes, runs, assumptions } = expand(prog, cat);
  const cost = runs.reduce((a, r) => a + r.est_s, 0);

  let refusal = null;
  if (runs.length > maxRuns) {
    const lastSweep = prog.sweeps.length ? prog.sweeps[prog.sweeps.length - 1].line : 1;
    refusal = {
      code: 'max_runs_exceeded',
      line: prog.allow ? prog.allow.line : lastSweep,
      message: `This program expands to ${runs.length} runs, over the ${maxRuns}-run ceiling. A two-axis sweep is trivially a hundred PDE solves; the executor refuses to launch that without you saying so in the program.`,
      remedy: `allow runs ${runs.length}`,
      remedy_text: `allow runs ${runs.length}`,
      remedy_line: prog.allow ? prog.allow.line : lastSweep,
      run_count: runs.length,
      max_runs: maxRuns,
    };
  }

  return {
    ok: !refusal,
    stub: true,
    run_count: runs.length,
    max_runs: maxRuns,
    estimated_cost_s: Math.round(cost),
    // NOT measured, on this machine or anywhere else: `cost` is the hardcoded
    // polynomial in this file. The previous wording ("on this machine's
    // measured reconstruction timings") claimed a measurement that never
    // happened, inside the string whose job is to disclose that it is a guess.
    cost_note: `GUESSED at ${runs.length ? Math.round(cost / runs.length) : 0} s per run by a hardcoded formula in js/protocol/stub.js over section size and grid pitch. Nothing was timed: this machine's reconstruction speed has not been measured. It is not an estimate from data, it is an ASSUMPTION.`,
    axes,
    runs: runs.map((r) => ({
      index: r.index, coords: r.coords, label: r.label, drug: r.drug,
      params: r.params, spec: r.spec, est_s: r.est_s, report: r.report,
    })),
    metrics: (runs[0] ? runs[0].report : Object.keys(METRICS)).map((k) => ({ key: k, ...METRICS[k] })),
    expects: prog.expects,
    assumptions,
    caveats: standingCaveats(prog, runs),
    refusal,
    errors: [],
    warnings,
  };
}

/* ------------------------------------------------------ public: execute -- */

const TERMINAL = new Set(['done', 'failed', 'cancelled']);
const CONCURRENCY = 4;
const sweeps = new Map();
let sweepSeq = 0;

export async function start(source, confirmRunCount) {
  const p = await plan(source);
  if (!p.ok) {
    const err = new Error(p.refusal ? p.refusal.message : 'the program does not validate');
    err.status = 422;
    err.payload = p;
    throw err;
  }
  if (confirmRunCount !== undefined && confirmRunCount !== p.run_count) {
    const err = new Error(`plan/confirm mismatch: you confirmed ${confirmRunCount} runs, the program now expands to ${p.run_count}. Re-plan before running.`);
    err.status = 409;
    throw err;
  }
  sweepSeq += 1;
  const id = 's_' + String(sweepSeq).padStart(4, '0');
  const sw = {
    sweep_id: id,
    status: 'running',
    started_at: Date.now() / 1000,
    axes: p.axes,
    metrics: p.metrics,
    caveats: p.caveats,
    assumptions: p.assumptions,
    expects: p.expects,
    max_runs: p.max_runs,
    estimated_cost_s: p.estimated_cost_s,
    cancelled: false,
    runs: p.runs.map((r) => ({
      index: r.index, coords: r.coords, label: r.label, drug: r.drug,
      params: r.params, spec: r.spec, est_s: r.est_s, report: r.report,
      run_id: null, status: 'queued', progress: 0, metrics: null, result: null,
      error: null, runtime_s: null, submitting: false,
    })),
  };
  sweeps.set(id, sw);
  pump(sw);
  return { sweep_id: id, run_count: sw.runs.length, stub: true, runs: sw.runs.map((r) => ({ index: r.index, run_id: null, status: 'queued' })) };
}

function mapMetrics(m) {
  if (!m) return null;
  const out = {};
  for (const [k, meta] of Object.entries(METRICS)) out[k] = m[meta.spec] ?? null;
  return out;
}

async function pump(sw) {
  if (sw.cancelled) return;
  // "Active" includes a run whose POST is still in flight, otherwise the next
  // tick would over-subscribe the queue before the first submission returns.
  const active = sw.runs.filter((r) => r.submitting || (r.run_id && !TERMINAL.has(r.status)));
  const inflight = active.filter((r) => r.run_id && !TERMINAL.has(r.status));
  let slots = CONCURRENCY - active.length;

  for (const r of sw.runs) {
    if (slots <= 0) break;
    if (r.run_id || r.submitting || TERMINAL.has(r.status)) continue;
    r.submitting = true;
    slots -= 1;
    postRun(r.spec)
      .then((res) => { r.run_id = (res && (res.run_id || res.id)) || null; r.submitting = false; if (!r.run_id) { r.status = 'failed'; r.error = 'POST /api/runs returned no run_id'; } })
      .catch((err) => { r.submitting = false; r.status = 'failed'; r.error = String((err && err.message) || err); });
  }

  await Promise.all(inflight.map(async (r) => {
    try {
      const s = await getStatus(r.run_id);
      r.status = s.status || r.status;
      r.progress = s.progress || 0;
      r.stage = s.stage || null;
      if (r.status === 'done' && !r.metrics) {
        const full = await getRun(r.run_id);
        r.result = full;
        r.metrics = mapMetrics(full.metrics);
        r.runtime_s = full.runtime_s;
      } else if (r.status === 'failed' && !r.error) {
        const full = await getRun(r.run_id).catch(() => null);
        r.error = (full && full.error) || 'the run failed and returned no message';
      }
    } catch (err) {
      // A transient poll failure is not a run failure. Leave the status alone.
      r.poll_error = String((err && err.message) || err);
    }
  }));

  const done = sw.runs.every((r) => TERMINAL.has(r.status));
  if (done) {
    sw.status = sw.runs.some((r) => r.status === 'failed') ? 'partial' : 'done';
    sw.finished_at = Date.now() / 1000;
    return;
  }
  setTimeout(() => pump(sw), 900);
}

export async function sweep(id) {
  const sw = sweeps.get(id);
  if (!sw) { const e = new Error(`no such sweep ${id}`); e.status = 404; throw e; }
  return snapshot(sw);
}

function snapshot(sw) {
  const counts = { queued: 0, running: 0, done: 0, failed: 0, cancelled: 0 };
  for (const r of sw.runs) counts[r.status] = (counts[r.status] || 0) + 1;
  return {
    sweep_id: sw.sweep_id,
    status: sw.status,
    stub: true,
    run_count: sw.runs.length,
    counts,
    axes: sw.axes,
    metrics: sw.metrics,
    caveats: sw.caveats,
    assumptions: sw.assumptions,
    expects: sw.expects,
    estimated_cost_s: sw.estimated_cost_s,
    elapsed_s: (sw.finished_at || Date.now() / 1000) - sw.started_at,
    runs: sw.runs.map((r) => ({
      index: r.index, coords: r.coords, label: r.label, drug: r.drug,
      params: r.params, spec: r.spec, run_id: r.run_id, status: r.status,
      progress: r.progress, stage: r.stage, metrics: r.metrics, result: r.result,
      error: r.error, runtime_s: r.runtime_s, est_s: r.est_s, report: r.report,
    })),
  };
}

export async function cancel(id) {
  const sw = sweeps.get(id);
  if (!sw) { const e = new Error(`no such sweep ${id}`); e.status = 404; throw e; }
  sw.cancelled = true;
  sw.status = 'cancelled';
  await Promise.all(sw.runs.map(async (r) => {
    if (TERMINAL.has(r.status)) return;
    if (r.run_id) { await cancelRun(r.run_id).catch(() => {}); }
    r.status = 'cancelled';
  }));
  return snapshot(sw);
}

/* ------------------------------------------------------- public: upload -- */

const ATOMIC = {
  H: 1.008, B: 10.81, C: 12.011, N: 14.007, O: 15.999, F: 18.998, Na: 22.99,
  Mg: 24.305, Si: 28.085, P: 30.974, S: 32.06, Cl: 35.45, K: 39.098, Ca: 40.078,
  Fe: 55.845, Br: 79.904, I: 126.904, Pt: 195.084,
};

/** Heavy-atom composition of a SMILES string. Implicit hydrogens are NOT counted. */
function smilesAtoms(smi) {
  const counts = {};
  const re = /(Cl|Br|Si|Na|Mg|Ca|Fe|Pt|[BCNOPSFIKH])|\[([A-Z][a-z]?)[^\]]*\]/g;
  let m;
  let bracketed = 0;
  while ((m = re.exec(smi))) {
    const el = m[2] || m[1];
    if (m[2]) bracketed += 1;
    if (!el || !ATOMIC[el]) continue;
    counts[el] = (counts[el] || 0) + 1;
  }
  // Lower-case aromatic atoms: c n o s p
  for (const ch of smi.replace(/\[[^\]]*\]/g, '')) {
    const up = { c: 'C', n: 'N', o: 'O', s: 'S', p: 'P' }[ch];
    if (up) counts[up] = (counts[up] || 0) + 1;
  }
  return { counts, bracketed };
}

/** Atom block of an MDL molfile / SDF — element symbols are explicit there. */
function molAtoms(text) {
  const lines = text.split(/\r?\n/);
  if (lines.length < 4) return null;
  const counts = {};
  const n = parseInt(lines[3].slice(0, 3), 10);
  if (!Number.isFinite(n) || n <= 0) return null;
  for (let i = 4; i < 4 + n && i < lines.length; i++) {
    const el = lines[i].slice(31, 34).trim();
    if (!el) continue;
    counts[el] = (counts[el] || 0) + 1;
  }
  return { counts, n_atoms: n };
}

function mass(counts) {
  let m = 0;
  let unknown = [];
  for (const [el, n] of Object.entries(counts)) {
    if (!ATOMIC[el]) { unknown.push(el); continue; }
    m += ATOMIC[el] * n;
  }
  return { m, unknown };
}

const FIELD_ORDER = [
  ['molecular_weight', 'g/mol'],
  ['ic50', 'uM'],
  ['hill', ''],
  ['diffusivity', 'um2/s'],
  ['binding', ''],
  ['clearance', '1/s'],
  ['bbb_penetrance', ''],
];

function flatYaml(text) {
  const out = {};
  for (const raw of text.split(/\r?\n/)) {
    const line = raw.replace(/#.*$/, '').trim();
    const m = /^([A-Za-z_][\w.]*)\s*:\s*(.+)$/.exec(line);
    if (!m) continue;
    const v = m[2].trim().replace(/^["']|["']$/g, '');
    out[m[1]] = Number.isFinite(parseFloat(v)) && /^[-+0-9.eE]+$/.test(v) ? parseFloat(v) : v;
  }
  return out;
}

export async function upload(file) {
  const text = await file.text();
  const ext = (file.name.split('.').pop() || '').toLowerCase();
  const warnings = [];
  const caveats = [
    'STUB PARSER. The molecular layer (tmesim.molecular) is not reachable from a browser, so every field it would DERIVE is reported here as `flagged-default` with a null value rather than guessed. The real /api/protocol/upload resolves them.',
  ];

  let name = file.name.replace(/\.[^.]+$/, '');
  let smiles = null;
  let counts = null;
  let explicit = {};
  let format = ext;

  if (ext === 'smi' || ext === 'smiles') {
    const first = text.split(/\r?\n/).find((l) => l.trim() && !l.trim().startsWith('#'));
    if (!first) { const e = new Error('the .smi file has no SMILES line'); e.status = 422; throw e; }
    const parts = first.trim().split(/\s+/);
    smiles = parts[0];
    if (parts[1]) name = parts[1];
    const a = smilesAtoms(smiles);
    counts = a.counts;
    warnings.push({ message: 'Molecular weight is a HEAVY-ATOM sum: implicit hydrogens are not counted by the stub parser, so it is a lower bound on the true mass.' });
  } else if (ext === 'mol' || ext === 'sdf') {
    const first = text.split(/\r?\n/)[0];
    if (first && first.trim()) name = first.trim();
    const a = molAtoms(text);
    if (!a) { const e = new Error('not a readable MDL molfile (no counts line)'); e.status = 422; throw e; }
    counts = a.counts;
  } else if (ext === 'json') {
    let o;
    try { o = JSON.parse(text); } catch (err) { const e = new Error('JSON does not parse: ' + err.message); e.status = 422; throw e; }
    explicit = o && typeof o === 'object' ? o : {};
    if (explicit.name) name = explicit.name;
    if (explicit.smiles) { smiles = explicit.smiles; counts = smilesAtoms(smiles).counts; }
    format = 'json';
  } else if (ext === 'yaml' || ext === 'yml') {
    explicit = flatYaml(text);
    if (explicit.name) name = explicit.name;
    if (explicit.smiles) { smiles = explicit.smiles; counts = smilesAtoms(smiles).counts; }
    format = 'yaml';
  } else {
    const e = new Error(`unsupported file type ".${ext}" — accepted: .smi .mol .sdf .json .yaml`);
    e.status = 415;
    throw e;
  }

  const fields = {};
  for (const [key, unit] of FIELD_ORDER) {
    if (Object.prototype.hasOwnProperty.call(explicit, key) && Number.isFinite(Number(explicit[key]))) {
      fields[key] = {
        value: Number(explicit[key]), unit,
        provenance: 'user-supplied',
        source: `given explicitly in ${file.name}`,
        note: 'A value you supplied. It overrides anything the molecular layer would derive, and it is not checked against a measurement.',
      };
      continue;
    }
    if (key === 'molecular_weight' && counts) {
      const { m, unknown } = mass(counts);
      if (unknown.length) warnings.push({ message: `unrecognised element symbol(s) ${unknown.join(', ')} were skipped in the mass sum.` });
      fields[key] = {
        value: Math.round(m * 1000) / 1000, unit,
        provenance: 'computed',
        source: smiles ? 'heavy-atom sum over the parsed SMILES' : 'atom block of the MDL molfile',
        note: smiles
          ? 'Heavy-atom sum only — implicit hydrogens are NOT counted by the stub parser. This is a lower bound; the real molecular layer returns the exact mass.'
          : 'Sum over the explicit atom block. Molfiles usually omit hydrogens, in which case this is a lower bound.',
      };
      continue;
    }
    fields[key] = {
      value: null, unit,
      provenance: 'flagged-default',
      source: 'not resolved',
      note: 'NOT RESOLVED. Nothing in the uploaded file gives this field and the browser stub cannot derive it. The executor either computes it through the molecular layer or substitutes a flagged default — either way it must be visible here, before the run, not discovered in a result.',
    };
  }

  if (fields.ic50.provenance === 'flagged-default') {
    warnings.push({ message: 'No IC50 in the file. Add `ic50 <n> uM` to the drug block (recorded as user-supplied) or define the drug `from library "..."` and use `ic50 from gdsc` (recorded as measured with its n).' });
  }

  return {
    stub: true,
    file: { name: file.name, size: file.size, format },
    drug: {
      name: name.replace(/[^A-Za-z0-9_]/g, '_').slice(0, 24) || 'UPLOADED',
      source: smiles ? `smiles "${smiles}"` : `file "${file.name}"`,
      smiles,
      composition: counts,
      fields,
    },
    warnings,
    caveats,
  };
}
