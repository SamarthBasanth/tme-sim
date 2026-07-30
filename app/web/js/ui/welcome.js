// welcome.js — a short, multi-step orientation tour for a first-time visitor.
//
// Distinct from apikey.js's gate (which BLOCKS everything until resolved,
// because nothing works without a key): this is advisory only, dismissible
// at any step (its own button, the backdrop, or Escape), and never blocks a
// returning or impatient user from just using the tool.
//
// Shown automatically once per browser (localStorage flag), and replayable
// any time afterward from the "tour" button next to "keys" in the top bar
// (see main.js's wiring and index.html's #btn-tour).

import { h } from '../util/dom.js';

const SEEN_KEY = 'tmesim_welcome_seen_v1';

export function shouldShowWelcome() {
  try { return !localStorage.getItem(SEEN_KEY); }
  catch { return false; } // no localStorage (private mode) -- don't nag every load
}

function markSeen() {
  try { localStorage.setItem(SEEN_KEY, '1'); } catch { /* fine, just re-shows next time */ }
}

/** Each step: a title, body content (array of nodes/strings for h()), and an
 * optional CSS selector to visually point the viewer at while this step is
 * showing (adds `.tour-highlight` to the matched element; removed on the
 * next step or on close). Selectors are best-effort -- if an element isn't
 * found (different viewport, a panel collapsed), the step still shows, just
 * without the highlight. */
const STEPS = [
  {
    title: 'Welcome to Vivum',
    selector: null,
    body: () => [
      h('p',
        'This dashboard runs a real physics simulation on a real, measured ' +
        'tumor tissue section — the cell viewer, the mechanism breakdown, ' +
        'and every metric are computed, not illustrative. This tour is ' +
        `${STEPS.length} short steps; skip it any time.`),
    ],
  },
  {
    title: '1 · Pick a tissue section',
    selector: '#f-section',
    body: () => [
      h('p',
        'The ', h('b', 'SECTION'), ' dropdown (top-left of the setup panel) ' +
        'is not one fixed demo tumor — it currently lists ',
        h('b', '38 real sections across 13 patients'),
        ', spanning glioblastoma (multiple assays, 12 patients) and breast ' +
        'cancer. Switching sections loads genuinely different measured ' +
        'tissue, not a re-skin of the same one.'),
    ],
  },
  {
    title: '2 · Configure and run',
    selector: '#btn-submit-run',
    body: () => [
      h('p',
        'Pick a drug, adjust dose/schedule if you want, then ',
        h('b', 'Submit run'),
        '. The viewer, mechanism breakdown, and metrics populate once it ' +
        'finishes — usually well under a minute at a coarse grid spacing.'),
    ],
  },
  {
    title: '3 · Read the CAVEATS panel',
    selector: '.caveats',
    body: () => [
      h('p',
        'The ', h('b', 'CAVEATS'), ' panel (right rail) states plainly what ' +
        'each number can and can’t be read as for the run you have open ' +
        '— read it before quoting anything. This project treats an honest ' +
        '"we don’t know" as a valid result, not a gap to hide.'),
    ],
  },
  {
    title: '4 · Script it — MATLAB-style',
    selector: '#btn-protocol',
    body: () => [
      h('p',
        'The ', h('b', '"protocol"'), ' tab (top bar) is a real code editor ' +
        'for defining and sweeping drug/dose experiments — 6 working ' +
        'templates to start from, syntax highlighting, inline diagnostics.'),
    ],
  },
];

function clearHighlight() {
  document.querySelectorAll('.tour-highlight').forEach((el) => el.classList.remove('tour-highlight'));
}

function highlight(selector) {
  clearHighlight();
  if (!selector) return;
  const el = document.querySelector(selector);
  if (el) el.classList.add('tour-highlight');
}

/**
 * @param {() => void} onGoToProtocol - called if the user reaches the final
 *   step and clicks its shortcut action, in addition to closing the tour.
 */
export function createWelcomeOverlay({ onGoToProtocol } = {}) {
  let i = 0;

  function close() {
    markSeen();
    clearHighlight();
    overlay.remove();
    document.removeEventListener('keydown', onKey);
  }
  function onKey(ev) {
    if (ev.key === 'Escape') close();
    else if (ev.key === 'ArrowRight') next();
    else if (ev.key === 'ArrowLeft') back();
  }
  function render() {
    const step = STEPS[i];
    titleEl.textContent = step.title;
    bodyEl.innerHTML = '';
    for (const node of step.body()) bodyEl.appendChild(node);
    dotsEl.innerHTML = '';
    STEPS.forEach((_, idx) => {
      dotsEl.appendChild(h(`span.tour-dot${idx === i ? '.active' : ''}`));
    });
    backBtn.disabled = i === 0;
    nextBtn.textContent = i === STEPS.length - 1 ? 'Open the protocol editor' : 'Next';
    skipBtn.style.display = i === STEPS.length - 1 ? 'none' : '';
    highlight(step.selector);
  }
  function next() {
    if (i === STEPS.length - 1) {
      close();
      if (onGoToProtocol) onGoToProtocol();
      return;
    }
    i = Math.min(i + 1, STEPS.length - 1);
    render();
  }
  function back() { i = Math.max(i - 1, 0); render(); }

  const titleEl = h('h3');
  const bodyEl = h('div.tour-body');
  const dotsEl = h('div.tour-dots');
  const backBtn = h('button.btn.ghost', { type: 'button', onclick: back }, 'Back');
  const nextBtn = h('button.btn.primary', { type: 'button', onclick: next }, 'Next');
  const skipBtn = h('button.btn.ghost', { type: 'button', onclick: close }, 'Skip tour');

  const card = h('div.gate-card.notice.tour-card', { onclick: (ev) => ev.stopPropagation() },
    titleEl,
    bodyEl,
    dotsEl,
    h('div.btn-row', { style: { marginTop: 'var(--sp-5)', justifyContent: 'space-between' } },
      h('div.btn-row', backBtn, nextBtn),
      skipBtn,
    ),
  );

  const overlay = h('div.gate-overlay.welcome-overlay', { onclick: close }, card);
  document.addEventListener('keydown', onKey);
  render();

  return { node: overlay, close };
}
