// editor.js — a hand-written code editor for TPL. No CodeMirror, no Monaco,
// no npm, no CDN (a strict CSP blocks external hosts, and a focused editor for
// one small language is smaller than the loader for a general one).
//
// ARCHITECTURE
// ------------
// Three layers share one coordinate system:
//
//   .ce-gutter   line numbers + diagnostic markers + block-match marks
//   pre.ce-hl    the syntax-highlighted rendering (pointer-transparent)
//   textarea.ce-ta   a real textarea, transparent text, visible caret
//
// The textarea stays a textarea. That keeps native selection, IME, spellcheck
// suppression, screen-reader behaviour, and — critically — the native undo
// stack, which a contenteditable re-implementation always gets wrong. Edits we
// make on the user's behalf (Tab, auto-indent, comment toggle, quick fixes) go
// through `execCommand('insertText')` where it is available precisely so they
// land on that same undo stack.
//
// The highlight layer is translated by the textarea's scroll offset rather than
// being scrolled itself, so there is exactly one scroller and no drift.
//
// Diagnostics arrive from POST /api/protocol/validate carrying a 1-based line
// (and optionally a column range). They land in three places at once: a gutter
// marker on that exact line, a tinted band across the line, and an underline
// under the offending columns. A refusal is a diagnostic, not a toast.

import { h, fill, debounce } from '../util/dom.js';
import { tokenizeLine, blockModel, indentAfter, indentForEnd, INDENT } from './tpl.js';

const MOD = navigator.platform && /mac/i.test(navigator.platform) ? 'Cmd' : 'Ctrl';

const SHORTCUTS = [
  [`${MOD}+Enter`, 'plan + run'],
  [`${MOD}+S`, 'validate now'],
  ['Tab / Shift+Tab', 'indent / outdent'],
  [`${MOD}+/`, 'toggle comment'],
  ['F8 / Shift+F8', 'next / previous diagnostic'],
  [`${MOD}+G`, 'go to line'],
  ['Esc', 'leave the editor'],
];

/**
 * @param {object} opts
 * @param {string} [opts.value]
 * @param {string} [opts.title]
 * @param {(text:string)=>void} [opts.onChange]  debounced text change
 * @param {(text:string)=>void} [opts.onInput]   every keystroke
 * @param {()=>void} [opts.onRun]
 * @param {()=>void} [opts.onValidateNow]
 */
export function createEditor(opts = {}) {
  const onChange = opts.onChange || (() => {});
  const onInput = opts.onInput || (() => {});

  let diagnostics = [];
  let charW = 7;      // measured below; monospace, so one number is enough
  let lineH = 16;

  /* ------------------------------------------------------------- DOM ---- */

  const ta = h('textarea.ce-ta', {
    spellcheck: 'false',
    autocomplete: 'off',
    autocapitalize: 'off',
    autocorrect: 'off',
    wrap: 'off',
    'aria-label': opts.ariaLabel || 'TPL protocol source',
    'aria-describedby': 'ce-legend',
  });
  ta.value = opts.value || '';

  const hl = h('pre.ce-hl', { 'aria-hidden': 'true' });
  const ruler = h('span.ce-ruler', { 'aria-hidden': 'true' }, '0'.repeat(40));
  const gutter = h('div.ce-gutter', { 'aria-hidden': 'true' });
  const code = h('div.ce-code', hl, ruler, ta);
  const scroll = h('div.ce-scroll', gutter, code);

  const posEl = h('span.mono.dim', 'ln 1, col 1');
  const statEl = h('span.ce-stat', 'not yet validated');
  const head = h('div.ce-head',
    h('span.ce-title', opts.title || 'protocol.tpl'),
    h('span.grow'),
    posEl,
    h('span.vrule', { 'aria-hidden': 'true' }),
    statEl);

  const legend = h('div#ce-legend.ce-legend', { role: 'note', 'aria-label': 'Editor keyboard shortcuts' },
    SHORTCUTS.map(([k, what]) => h('span.ce-key', h('kbd', k), h('span', what))));

  const el = h('div.ce', head, scroll, legend);

  /* --------------------------------------------------------- geometry --- */

  function measure() {
    const r = ruler.getBoundingClientRect();
    if (r.width > 0) charW = r.width / 40;
    const lr = hl.firstChild && hl.firstChild.getBoundingClientRect();
    if (lr && lr.height > 0) lineH = lr.height;
  }

  /* -------------------------------------------------------- rendering --- */

  const lineDiags = () => {
    const map = new Map();
    for (const d of diagnostics) {
      if (!d || !d.line) continue;
      const cur = map.get(d.line);
      const sev = d.severity === 'error' ? 2 : d.severity === 'warning' ? 1 : 0;
      if (!cur || sev > cur.rank) map.set(d.line, { rank: sev, list: (cur ? cur.list : []).concat(d) });
      else cur.list.push(d);
    }
    return map;
  };

  function caretLine() {
    return ta.value.slice(0, ta.selectionStart).split('\n').length;
  }

  function paint() {
    const text = ta.value;
    const lines = text.split('\n');
    const dmap = lineDiags();
    const { partner } = blockModel(text);
    const cur = caretLine();
    const mate = partner.get(cur) || 0;

    /* --- highlight layer -------------------------------------------- */
    const frag = document.createDocumentFragment();
    lines.forEach((line, i) => {
      const n = i + 1;
      const d = dmap.get(n);
      let cls = 'ce-l';
      if (n === cur) cls += ' cur';
      if (mate && (n === cur || n === mate)) cls += ' blk';
      if (d) cls += d.rank === 2 ? ' has-err' : d.rank === 1 ? ' has-warn' : ' has-info';
      const span = h('span.' + cls.split(' ').join('.'));

      if (line.length === 0) {
        span.appendChild(document.createTextNode('​'));
      } else {
        for (const t of tokenizeLine(line)) {
          if (t.t === 'ws') { span.appendChild(document.createTextNode(line.slice(t.s, t.e))); continue; }
          span.appendChild(h('span.tk-' + t.t, line.slice(t.s, t.e)));
        }
      }

      // Column-precise underline for any diagnostic that gave us a range.
      if (d) {
        for (const one of d.list) {
          if (!Number.isFinite(one.col)) continue;
          const c0 = Math.max(0, one.col);
          const c1 = Math.max(c0 + 1, Number.isFinite(one.end_col) ? one.end_col : c0 + 1);
          const u = h('i.ce-under.' + (one.severity === 'warning' ? 'warn' : 'err'));
          u.style.left = (c0 * charW) + 'px';
          u.style.width = ((c1 - c0) * charW) + 'px';
          span.appendChild(u);
        }
      }
      frag.appendChild(span);
    });
    fill(hl, frag);

    /* --- gutter ------------------------------------------------------ */
    const gfrag = document.createDocumentFragment();
    lines.forEach((_, i) => {
      const n = i + 1;
      const d = dmap.get(n);
      const isBlk = mate && (n === cur || n === mate);
      const g = h('div.ce-g' + (n === cur ? '.cur' : '') + (isBlk ? '.blk' : '')
        + (d ? (d.rank === 2 ? '.err' : d.rank === 1 ? '.warn' : '.info') : ''),
        { dataset: { line: String(n) } },
        h('span.ce-gm', d ? (d.rank === 2 ? '×' : d.rank === 1 ? '!' : 'i') : ''),
        h('span.ce-gn', String(n)));
      if (d) {
        g.title = d.list.map((x) => `ln ${x.line}: ${x.message}`).join('\n\n');
        g.setAttribute('role', 'button');
        g.tabIndex = 0;
      }
      gfrag.appendChild(g);
    });
    fill(gutter, gfrag);

    syncScroll();
    posEl.textContent = `ln ${cur}, col ${caretCol()}  ·  ${lines.length} ln`;
  }

  /** Cheap repaint: only the current-line / block-match decoration moved. */
  function paintCaret() {
    const cur = caretLine();
    const { partner } = blockModel(ta.value);
    const mate = partner.get(cur) || 0;
    const ls = hl.children;
    const gs = gutter.children;
    for (let i = 0; i < ls.length; i++) {
      const n = i + 1;
      ls[i].classList.toggle('cur', n === cur);
      ls[i].classList.toggle('blk', !!mate && (n === cur || n === mate));
      if (gs[i]) {
        gs[i].classList.toggle('cur', n === cur);
        gs[i].classList.toggle('blk', !!mate && (n === cur || n === mate));
      }
    }
    posEl.textContent = `ln ${cur}, col ${caretCol()}  ·  ${ls.length} ln`;
  }

  function caretCol() {
    const s = ta.value.lastIndexOf('\n', ta.selectionStart - 1);
    return ta.selectionStart - s;
  }

  function syncScroll() {
    hl.style.transform = `translate(${-ta.scrollLeft}px, ${-ta.scrollTop}px)`;
    gutter.scrollTop = ta.scrollTop;
  }

  /* ------------------------------------------------------------ edits --- */

  /**
   * Replace [from,to) with `text` and leave the caret where `sel` says.
   * Routed through execCommand when available so the browser's own undo stack
   * records it — an editor that loses undo on Tab is worse than no editor.
   */
  function replaceRange(from, to, text, selStart, selEnd) {
    ta.focus();
    ta.setSelectionRange(from, to);
    let ok = false;
    try { ok = document.execCommand('insertText', false, text); } catch { ok = false; }
    if (!ok) {
      ta.value = ta.value.slice(0, from) + text + ta.value.slice(to);
    }
    const a = selStart === undefined ? from + text.length : selStart;
    const b = selEnd === undefined ? a : selEnd;
    ta.setSelectionRange(a, b);
    afterEdit();
  }

  function lineBounds(pos) {
    const s = ta.value.lastIndexOf('\n', pos - 1) + 1;
    let e = ta.value.indexOf('\n', pos);
    if (e < 0) e = ta.value.length;
    return [s, e];
  }

  function selectedLineRange() {
    const [s] = lineBounds(ta.selectionStart);
    const [, e] = lineBounds(ta.selectionEnd);
    return [s, e];
  }

  function indentSelection(out) {
    const [s, e] = selectedLineRange();
    const lines = ta.value.slice(s, e).split('\n');
    let first = 0;
    let total = 0;
    const next = lines.map((ln, i) => {
      if (out) {
        const m = /^(\t| {1,2})/.exec(ln);
        const cut = m ? m[0].length : 0;
        if (i === 0) first = cut;
        total += cut;
        return ln.slice(cut);
      }
      if (!ln.length && lines.length > 1) return ln;   // never indent a blank line
      if (i === 0) first = INDENT.length;
      total += INDENT.length;
      return INDENT + ln;
    }).join('\n');
    const sign = out ? -1 : 1;
    const a = Math.max(s, ta.selectionStart + sign * first);
    const b = Math.max(a, ta.selectionEnd + sign * total);
    replaceRange(s, e, next, a, b);
  }

  function toggleComment() {
    const [s, e] = selectedLineRange();
    const lines = ta.value.slice(s, e).split('\n');
    const meaningful = lines.filter((l) => l.trim().length);
    const allCommented = meaningful.length > 0 && meaningful.every((l) => /^\s*#/.test(l));
    const next = lines.map((l) => {
      if (!l.trim().length) return l;
      if (allCommented) return l.replace(/^(\s*)#\s?/, '$1');
      const ind = /^\s*/.exec(l)[0];
      return `${ind}# ${l.slice(ind.length)}`;
    }).join('\n');
    replaceRange(s, e, next, s, s + next.length);
  }

  function newline() {
    const [s] = lineBounds(ta.selectionStart);
    const before = ta.value.slice(s, ta.selectionStart);
    replaceRange(ta.selectionStart, ta.selectionEnd, '\n' + indentAfter(before));
  }

  /** After typing `end`, snap the line back to its opener's indent. */
  function maybeSnapEnd() {
    const [s, e] = lineBounds(ta.selectionStart);
    const line = ta.value.slice(s, e);
    if (!/^\s*end\s*$/.test(line)) return;
    const n = ta.value.slice(0, s).split('\n').length;
    const want = indentForEnd(ta.value, n) + 'end';
    if (line === want) return;
    const caretAtEnd = ta.selectionStart === e;
    replaceRange(s, e, want, caretAtEnd ? s + want.length : ta.selectionStart);
  }

  const emitChange = debounce(() => onChange(ta.value), opts.debounceMs || 350);

  function afterEdit() {
    paint();
    onInput(ta.value);
    emitChange();
  }

  /* --------------------------------------------------------- listeners -- */

  let snapping = false;
  ta.addEventListener('input', () => {
    if (!snapping) {
      snapping = true;
      try { maybeSnapEnd(); } finally { snapping = false; }
    }
    paint();
    onInput(ta.value);
    emitChange();
  });
  ta.addEventListener('scroll', syncScroll);
  ta.addEventListener('click', paintCaret);
  ta.addEventListener('keyup', (e) => {
    if (e.key.startsWith('Arrow') || e.key === 'Home' || e.key === 'End'
      || e.key === 'PageUp' || e.key === 'PageDown') paintCaret();
  });
  ta.addEventListener('select', paintCaret);

  // Tabs would break the column arithmetic the underlines and the gutter rely
  // on, and TPL indentation is two spaces by convention. Normalise on paste.
  ta.addEventListener('paste', (e) => {
    const txt = e.clipboardData && e.clipboardData.getData('text/plain');
    if (!txt || !txt.includes('\t')) return;
    e.preventDefault();
    replaceRange(ta.selectionStart, ta.selectionEnd, txt.replace(/\t/g, INDENT));
  });

  ta.addEventListener('keydown', (e) => {
    const mod = e.ctrlKey || e.metaKey;

    if (e.key === 'Tab') {
      e.preventDefault();
      if (e.shiftKey) { indentSelection(true); return; }
      if (ta.selectionStart !== ta.selectionEnd) { indentSelection(false); return; }
      replaceRange(ta.selectionStart, ta.selectionEnd, INDENT);
      return;
    }
    if (e.key === 'Enter' && !mod) {
      e.preventDefault();
      newline();
      return;
    }
    if (mod && e.key === 'Enter') {
      e.preventDefault(); e.stopPropagation();
      if (opts.onRun) opts.onRun();
      return;
    }
    if (mod && (e.key === 's' || e.key === 'S')) {
      e.preventDefault(); e.stopPropagation();
      emitChange.flush();
      if (opts.onValidateNow) opts.onValidateNow();
      return;
    }
    if (mod && e.key === '/') {
      e.preventDefault(); e.stopPropagation();
      toggleComment();
      return;
    }
    if (mod && (e.key === 'g' || e.key === 'G')) {
      e.preventDefault(); e.stopPropagation();
      const n = window.prompt('Go to line');
      const v = parseInt(n, 10);
      if (Number.isFinite(v)) gotoLine(v);
      return;
    }
    if (e.key === 'F8') {
      e.preventDefault();
      stepDiagnostic(e.shiftKey ? -1 : 1);
      return;
    }
    if (e.key === 'Escape') {
      // Let the workspace handle it; just make sure the viewer's Esc (unpin)
      // does not also fire from inside the editor.
      e.stopPropagation();
      ta.blur();
    }
  });

  gutter.addEventListener('click', (e) => {
    const g = e.target.closest('.ce-g');
    if (g && g.dataset.line) gotoLine(Number(g.dataset.line));
  });
  gutter.addEventListener('keydown', (e) => {
    const g = e.target.closest('.ce-g');
    if (g && g.dataset.line && (e.key === 'Enter' || e.key === ' ')) {
      e.preventDefault();
      gotoLine(Number(g.dataset.line));
    }
  });

  /* ------------------------------------------------------ public API ---- */

  function offsetOfLine(n) {
    const lines = ta.value.split('\n');
    let pos = 0;
    for (let i = 0; i < n - 1 && i < lines.length; i++) pos += lines[i].length + 1;
    return pos;
  }

  function gotoLine(n, col, endCol) {
    const lines = ta.value.split('\n');
    const line = Math.max(1, Math.min(lines.length, n || 1));
    const pos = offsetOfLine(line);
    const c0 = Number.isFinite(col) ? Math.max(0, col) : 0;
    const c1 = Number.isFinite(endCol) ? endCol : (lines[line - 1] || '').length;
    ta.focus();
    ta.setSelectionRange(pos + Math.min(c0, (lines[line - 1] || '').length),
      pos + Math.min(c1, (lines[line - 1] || '').length));
    const view = scroll.clientHeight || 240;
    const target = (line - 1) * lineH - view / 2 + lineH;
    ta.scrollTop = Math.max(0, target);
    syncScroll();
    paintCaret();
  }

  let diagCursor = -1;
  function stepDiagnostic(dir) {
    const ordered = diagnostics.filter((d) => d && d.line)
      .slice().sort((a, b) => a.line - b.line);
    if (!ordered.length) return;
    diagCursor = (diagCursor + dir + ordered.length) % ordered.length;
    const d = ordered[diagCursor];
    gotoLine(d.line, d.col, d.end_col);
  }

  function setDiagnostics(list) {
    diagnostics = Array.isArray(list) ? list : [];
    diagCursor = -1;
    paint();
  }

  function setStatus(text, kind) {
    statEl.textContent = text;
    statEl.className = 'ce-stat' + (kind ? ' ' + kind : '');
  }

  function setValue(text) {
    if (ta.value === text) return;
    ta.value = String(text == null ? '' : text).replace(/\t/g, INDENT);
    paint();
    onInput(ta.value);
    emitChange();
  }

  /** Insert a line of TPL after `afterLine` (0 = at the very top). */
  function insertLine(afterLine, text) {
    const lines = ta.value.split('\n');
    const at = Math.max(0, Math.min(lines.length, afterLine));
    const pos = offsetOfLine(at + 1);
    replaceRange(pos, pos, text + '\n', pos + text.length);
    gotoLine(at + 1);
  }

  // First paint has to wait for layout so the ruler can be measured.
  requestAnimationFrame(() => { measure(); paint(); });
  paint();

  // charW/lineH are font metrics, not content: they cannot change between
  // keystrokes, only when the editor is actually resized (a window resize, a
  // splitter/pane drag, or the workspace sheet going from hidden to visible,
  // which ResizeObserver reports as a 0->nonzero box change same as a real
  // resize). measure() used to run inside paint() — a getBoundingClientRect()
  // pair forced on every single keystroke for a number that was, in the
  // overwhelming common case, identical to the one already held.
  if (typeof ResizeObserver === 'function') new ResizeObserver(measure).observe(el);

  return {
    el,
    ta,
    get value() { return ta.value; },
    set value(v) { setValue(v); },
    setValue,
    setDiagnostics,
    setStatus,
    gotoLine,
    insertLine,
    focus: () => ta.focus(),
    flush: () => emitChange.flush(),
    repaint: paint,
  };
}

export { SHORTCUTS };
