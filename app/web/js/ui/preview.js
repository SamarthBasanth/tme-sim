// preview.js — the Setup rail's coarse-solve preview.
//
// WHERE THIS LIVES, AND WHY THAT IS THE DESIGN
//
// In the COMMIT group of the SETUP rail, directly above `Active job`. Not in
// the Read rail. That placement is the main honesty mechanism in this file and
// it is worth more than any of the wording below.
//
// The Read rail is where this product reports results: a run id, a provenance
// chip per number, a caveat list, an export, a copy button that writes a
// quotable sentence into the clipboard. Putting a preview number there would
// put an UNQUOTABLE number inside the surface whose whole job is to produce
// quotable ones, one panel away from a copy button. The Setup rail is where
// you decide what to spend 59 seconds on, and "what will this roughly say"
// is a question about the decision, not about the result.
//
// So:
//   * a preview number NEVER appears in the Read rail, in a report, on a
//     project card, or in `/api/runs` — the server mints no run_id for it and
//     stores nothing (app/server/preview.py);
//   * the preview is REPLACED, not kept beside. `supersede()` drops the number
//     the instant a real run is submitted and leaves the run id in its place.
//     Keeping the two numbers side by side was the alternative and it is the
//     worse one: two values for one readout, 0.04–1.76 % apart, invites the
//     reader to treat the gap as a result and to quote whichever is handier.
//     The preview's job ends when the real solve starts.
//
// THE BROWSER DOES NOT DECIDE WHAT IS PREVIEWABLE.
//
// There is no copy of the Visium rule, the readout rule or the measured-pair
// rule in this file. The button always asks, and whatever the server answers —
// a number or a refusal — is rendered in full. A second copy of a refusal in
// the client is exactly the drift this repository keeps finding, and a refusal
// that the client can get wrong is a refusal that can be relaxed by accident.
//
// EVERY RENDERING OF THE NUMBER CARRIES TWO THINGS, OR IT IS NOT RENDERED.
//
//   1. the word PREVIEW, at no smaller than the number's own size;
//   2. the MEASURED disagreement with the real solve.
//
// Both ride on `.caveat.sev-2`, whose `.cv-key` css already resolves to
// `--fs-lg` — the same token `.metric .m-val` uses — precisely so that a hard
// stop is "never smaller than the number it qualifies"
// (docs/WORKSPACE_CONTRACT.md §6; app/web/css/components.css, the nine measured
// violations table). That is why this file adds no stylesheet: the rule is
// already encoded in a class, and spelling a font-size here would be a second
// copy of it.
//
// `assertCaveatNotSmaller()` checks it in the live browser rather than trusting
// the class, because the class is what a future edit overrides.

import { h, fill, clear } from '../util/dom.js';
import { num, int, dur } from '../util/fmt.js';
import { sectTitle } from './rail.js';
import { api, ApiError } from '../api.js';

/** The number's own size token. The caveat beside it may not go below this. */
const VALUE_FS = 'var(--fs-lg)';

/* --------------------------------------------------------------- fetch --- */

/**
 * POST /api/preview, with the same two credentials every other call carries.
 *
 * Written here rather than added to js/api.js because that module is shared and
 * this is the only caller. It mirrors `rawFetch`'s credential handling exactly:
 * the session cookie rides via `credentials: 'same-origin'` and the stored key
 * via `x-api-key`, and the server decides which it honours
 * (docs/AUTH_CONTRACT.md §0).
 */
async function postPreview(spec, { signal } = {}) {
  const headers = { 'content-type': 'application/json', accept: 'application/json' };
  const key = api.sendKey ? api.apiKey : '';
  if (key) headers['x-api-key'] = key;
  let res;
  try {
    res = await fetch(api.base + '/api/preview', {
      method: 'POST',
      headers,
      body: JSON.stringify(spec || {}),
      signal,
      cache: 'no-store',
      credentials: 'same-origin',
    });
  } catch (e) {
    if (e && e.name === 'AbortError') throw e;
    throw new ApiError((e && e.message) || 'the network request did not complete',
      { status: 0, path: '/api/preview', cause: e });
  }
  let payload = null;
  try { payload = await res.json(); } catch { /* the status is still the story */ }
  if (!res.ok) {
    const d = payload && payload.detail;
    const msg = (d && typeof d === 'object' && d.error) || (typeof d === 'string' && d)
      || `${res.status} ${res.statusText}`;
    throw new ApiError(msg, { status: res.status, path: '/api/preview', payload });
  }
  return payload;
}

/* ---------------------------------------------------------------- parts -- */

/**
 * A hard stop, in the house markup, so it inherits the house sizes.
 *
 * `.caveat.sev-2` → `.cv-key` at --fs-lg, `.cv-body` at --fs-sm. Nothing here
 * sets a font-size; that is the point.
 */
function hardStop(subject, body) {
  return h('div.caveat.sev-2',
    h('div.cv-key', h('span.cv-subj', subject)),
    h('div.cv-body', body));
}

/**
 * A prose note, at --fs-sm.
 *
 * NOT --fs-disclosure (10.5px). --fs-disclosure is the FLOOR for a sentence,
 * and on top of the floor sits the second rule: never smaller than the number
 * it qualifies. Every sentence in this panel qualifies the preview value, and
 * --fs-sm (11.5px) is the size app/web/css/components.css settled on for
 * exactly that reason after measuring nine violations of it. Shaped like
 * js/ui/reach.js::figNote, which carries its own styles inline for the same
 * reason this does: there is no class to drift from.
 */
function note(...kids) {
  return h('div.figmode', {
    style: {
      fontFamily: 'var(--font-ui)', fontSize: 'var(--fs-sm)',
      color: 'var(--fg-1)', lineHeight: 'var(--lh-body)',
      borderTop: '1px solid var(--line-1)',
      padding: 'var(--sp-2) var(--sp-5) 0', margin: 'var(--sp-2) 0 var(--sp-4)',
      overflowWrap: 'anywhere',
    },
  }, ...kids);
}

/**
 * The agreement, in the units a reader can check, with its n and its source.
 *
 * NOT "±1.8 %". The three order statistics and the absolute worst case, because
 * a single number here would be a choice about which one to show and the
 * median is the flattering one.
 */
function agreementBlock(ag) {
  if (!ag) return null;
  const rows = [
    ['median disagreement', `${num(ag.median_rel_err_pct, 4)} %`],
    ['p90', `${num(ag.p90_rel_err_pct, 4)} %`],
    ['worst measured', `${num(ag.max_rel_err_pct, 4)} %`],
    ['worst, absolute', `${num(ag.max_abs_err, 4)}`],
    ['over', `${int(ag.n_pairs)} pairs · ${int(ag.n_sections)} sections`
      + ` × ${int(ag.n_compounds)} compounds`],
  ];
  return h('div',
    // THE ERROR BAR IS ITSELF A HARD STOP, so its subject line resolves to
    // --fs-lg — the same size as the number it qualifies. A measured
    // disagreement set smaller than the figure it is about is the defect
    // app/web/css/components.css measured nine times and fixed.
    hardStop(
      `MAY DIFFER FROM THE REAL SOLVE BY UP TO ${num(ag.max_rel_err_pct, 4)} %`,
      `Measured, pair by pair, not estimated: median `
      + `${num(ag.median_rel_err_pct, 4)} %, p90 ${num(ag.p90_rel_err_pct, 4)} %, `
      + `worst ${num(ag.max_rel_err_pct, 4)} % — `
      + `${num(ag.max_abs_err, 4)} in absolute terms — over `
      + `${int(ag.n_pairs)} genuine dx ${num(ag.preview_dx_um, 0)} um / `
      + `dx ${num(ag.target_dx_um, 0)} um pairs on ${int(ag.n_sections)} `
      + `sections and ${int(ag.n_compounds)} compounds. This figure covers `
      + `${(ag.sections || []).join(', ')} with `
      + `${(ag.compounds || []).join(', ')} and nothing else; outside that set `
      + 'the server refuses a preview rather than stretching the figure over '
      + 'it.'),
    h('div.tblwrap', { style: { padding: '0 var(--sp-5)' } },
      h('table.tbl.zebra',
        h('thead', h('tr',
          h('th', `dx ${num(ag.preview_dx_um, 0)} um vs dx ${num(ag.target_dx_um, 0)} um`),
          h('th.right', 'value'))),
        h('tbody', rows.map(([k, v]) => h('tr',
          h('td.mono', k), h('td.right.num.mono', v)))))),
    // The two comparators that say whether that error bar is small, each named
    // with the record it comes from. Both are within-compound or across-donor;
    // neither is an ordering between compounds.
    note(`The worst case, ${num(ag.max_abs_err, 4)}, is `
      + `${num((ag.max_abs_err / ag.across_donor_range) * 100, 0)} % of this `
      + `readout's whole across-donor range (${num(ag.across_donor_range, 4)}) `
      + `and ${num((ag.max_abs_err / ag.nuisance_band_width) * 100, 0)} % of the `
      + `width of its own published nuisance band (${num(ag.nuisance_band_width, 3)}). `
      + 'Neither comparator is an ordering between compounds: the first is a '
      + 'range across donors, the second one compound\'s own span across every '
      + 'nuisance setting tested.'),
    note(`Cost: the coarse solve's median is ${dur(ag.coarse_median_s)} against `
      + `${dur(ag.target_median_s)} for the real run — `
      + `${num(ag.cost_ratio, 2)}x cheaper, measured at loadavg `
      + `${num(ag.loadavg_1m_min, 1)}–${num(ag.loadavg_1m_max, 1)} on `
      + `${int((ag.n_cpu || [])[0])} CPUs. The RATIO is the quotable part: this `
      + "repository's performance record measures five different wall-clocks "
      + 'for one unchanged RunSpec and concludes that the shape is stable and '
      + 'the absolutes are not.'),
    note(`Source: ${ag.source} · ${ag.record}. Recomputed from the stored pairs `
      + "on every request, never copied into the server's source — so if the "
      + 'measurement is gone the feature is gone, rather than the error bar '
      + `being gone. ${int(ag.pairs_discarded_same_grid)} stored pairs were `
      + 'discarded because both of their runs produced the same grid: a run '
      + 'compared with itself has an error of 0.000 % by construction, and '
      + 'mixing those in is how the study behind this feature first reported a '
      + 'perfect result.'));
}

/** What a coarse solve will NOT preview, so the absence is readable. */
function refusedBlock(refused) {
  if (!refused || !refused.length) return null;
  return h('details.sect',
    h('summary', sectTitle('Not previewable'),
      h('span.hint', `${int(refused.length)} readouts`)),
    h('div.sect-body',
      refused.map((r) => hardStop(r.readout, r.why))));
}

/* ---------------------------------------------------------------- panel -- */

/**
 * @param {object}   o
 * @param {() => object} o.getSpec  the RunSpec the form currently holds
 */
export function createPreviewPanel({ getSpec }) {
  const out = h('div');
  let busy = false;
  let inflight = null;

  const btn = h('button.btn', {
    type: 'button',
    title: 'POST /api/preview — the same solver on a 2.5x coarser grid. '
      + 'Nothing is stored and no run_id is minted.',
  }, 'Preview (coarse solve)');

  const head = h('div', {
    style: { display: 'flex', alignItems: 'center', gap: 'var(--sp-4)',
             padding: 'var(--sp-4) var(--sp-5)' },
  }, btn, h('span.grow'));

  /** The idle copy. Says what the button buys and what it does not. */
  function idle() {
    fill(out, note(
      'Runs the SAME solver on a 2.5x coarser grid and shows one readout — the '
      + 'delivery spread — with the measured difference from the real run '
      + 'beside it. It is not a model, nothing is fitted, and nothing is '
      + 'stored: a preview has no run_id and cannot be opened, exported or '
      + 'reported. The real solve is still the thing you end up with.'));
  }

  function renderRefusal(err) {
    const d = (err && err.payload && err.payload.detail) || {};
    fill(out,
      hardStop('NO PREVIEW HERE', err.message || String(err)),
      d.remedy ? note(
        `What to do instead: ${d.remedy}.`) : null,
      note(
        'The real solve is unaffected — Submit run below is the route to a '
        + 'number for this spec.'));
  }

  /**
   * The number, its label and its error bar, in that order in the DOM.
   *
   * ORDER IS LOAD-BEARING: the label is a SIBLING immediately after the value,
   * not a footnote further down a scroller, so there is no scroll position at
   * which the number is on screen and the word PREVIEW is not.
   */
  function renderPreview(p) {
    const ag = p.agreement || {};
    const value = p.value === null || p.value === undefined
      ? null
      : `${num(p.value, 2)}${p.unit || ''}`;

    fill(out,
      value === null
        ? hardStop('NO NUMBER FROM THE COARSE SOLVE', p.null_reason
            || 'the engine produced no delivery spread on this section.')
        : h('div', { style: { display: 'flex', alignItems: 'baseline',
                              gap: 'var(--sp-4)', padding: '0 var(--sp-5)' } },
          h('span.mono', {
            id: 'preview-value',
            // The ONE place this file names a size, and it names it so the
            // assertion below has something to compare the caveat against.
            style: { fontSize: VALUE_FS, color: 'var(--fg-0)',
                     fontVariantNumeric: 'tabular-nums' },
          }, value),
          h('span', { style: { fontSize: 'var(--fs-disclosure)', color: 'var(--fg-1)',
                               lineHeight: 'var(--lh-body)' } },
            h('b', 'DELIVERY SPREAD, PREVIEWED'), ' — the best-perfused tenth '
            + 'of tumour cells receives this many times the exposure of the '
            + 'worst-perfused tenth, on the coarse grid.')),

      // (1) the label, immediately after the value, at the value's own size.
      hardStop('PREVIEW — NOT A RESULT', p.preview_label),

      // (2) the measured agreement, which is the other half of the number.
      agreementBlock(ag),

      // Both terms' provenance and the denominator, per the house rule that a
      // fold-change travels with what it is a fraction of.
      note(
        `${p.derivation} Over `
        + (Number.isFinite(p.n_tumor) ? `${int(p.n_tumor)} tumour cells` : 'an unrecorded number of tumour cells')
        + `, grid ${(p.preview_grid || []).join('x') || 'unrecorded'}, solved in `
        + `${dur(p.preview_wall_clock_s)} here.`),

      refusedBlock(p.refused),

      (p.caveats || []).length
        ? h('details.sect',
          h('summary', sectTitle('The engine\'s own limits'),
            h('span.hint', `${int(p.caveats.length)} carried verbatim`)),
          h('div.sect-body', p.caveats.map((c) => hardStop('LIMIT', c))))
        : null);

    assertCaveatNotSmaller();
  }

  /**
   * docs/WORKSPACE_CONTRACT.md §6, checked in the browser.
   *
   * The rule — a caveat is never set smaller than the number it qualifies — is
   * carried by `.caveat.sev-2 .cv-key`'s --fs-lg, which is a CSS fact and
   * therefore something a later edit can change without touching this file. So
   * it is measured rather than assumed, and a violation is logged loudly and
   * marks the DOM, which is what the browser test reads.
   */
  function assertCaveatNotSmaller() {
    const val = out.querySelector('#preview-value');
    const key = out.querySelector('.caveat.sev-2 .cv-key');
    if (!val || !key || typeof getComputedStyle !== 'function') return;
    const vs = parseFloat(getComputedStyle(val).fontSize) || 0;
    const ks = parseFloat(getComputedStyle(key).fontSize) || 0;
    const ok = ks >= vs;
    out.dataset.caveatFloor = ok ? 'ok' : 'violated';
    if (!ok) {
      // eslint-disable-next-line no-console
      console.error(
        `preview: the caveat subject is ${ks}px against a ${vs}px number. `
        + 'docs/WORKSPACE_CONTRACT.md §6 — a caveat is never set smaller than '
        + 'the number it qualifies.');
    }
  }

  async function run() {
    if (busy) return;
    busy = true;
    btn.disabled = true;
    const prev = btn.textContent;
    btn.textContent = 'solving…';
    fill(out, note(
      'Solving on the coarse grid. The real run is the one you submit below; '
      + 'this is only ahead of it.'));
    const ctrl = new AbortController();
    inflight = ctrl;
    try {
      // The spec is sent AS THE FORM HOLDS IT. The server substitutes its own
      // dx_um and refuses anything it has not measured; the browser does not
      // pre-edit the spec, so what is previewed is the spec you are about to
      // submit.
      const p = await postPreview(getSpec(), { signal: ctrl.signal });
      renderPreview(p);
    } catch (err) {
      if (err && err.name === 'AbortError') { idle(); } else { renderRefusal(err); }
    } finally {
      busy = false;
      btn.disabled = false;
      btn.textContent = prev;
      inflight = null;
    }
  }

  btn.addEventListener('click', run);

  const el = h('details.sect', { open: true, id: 'sect-preview' },
    h('summary', sectTitle('Preview'), h('span.hint', 'coarse solve · not a result')),
    h('div.sect-body', { style: { padding: '0 0 var(--sp-4)' } }, head, out));

  idle();

  return {
    el,
    /**
     * The real solve has started: DROP THE PREVIEW NUMBER.
     *
     * Not "annotate it", not "keep it beside the real one". A readout with two
     * live values 0.04–1.76 % apart is an invitation to read the gap as a
     * finding and to quote whichever is nearer to hand. The preview answered a
     * question about a decision that has now been taken.
     */
    supersede(runId) {
      if (inflight) { try { inflight.abort(); } catch { /* already done */ } }
      fill(out, note(
        'Superseded'
        + (runId ? ` by ${runId}` : '')
        + ' — the real solve is running and its number is the one this product '
        + 'stands behind. The preview is dropped rather than shown beside it: '
        + 'one readout with two live values is an invitation to quote the '
        + 'convenient one.'));
      delete out.dataset.caveatFloor;
    },
    /** Back to the idle copy, e.g. after the spec changes. */
    reset: idle,
    clear() { clear(out); idle(); },
  };
}
