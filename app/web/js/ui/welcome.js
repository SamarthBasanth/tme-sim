// welcome.js — a short, multi-step orientation tour for a first-time visitor.
//
// A floating popover anchored to the real UI element each step is about
// (the SECTION dropdown, the Submit run button, ...), NOT a full-screen
// modal -- the rest of the dashboard stays fully visible underneath so
// "look at the dropdown top-left" actually has a real dropdown to look at
// in context, rather than a blurred backdrop covering it. Distinct from
// apikey.js's gate, which legitimately DOES need to block everything
// (nothing works without a key); this is advisory only.
//
// Shown automatically once per browser (localStorage flag), and replayable
// any time afterward from the "tour" button next to "keys" in the top bar
// (see main.js's wiring and index.html's #btn-tour).

import { h } from '../util/dom.js';

const SEEN_KEY = 'tmesim_welcome_seen_v1';
const MARGIN = 14;

export function shouldShowWelcome() {
  try { return !localStorage.getItem(SEEN_KEY); }
  catch { return false; } // no localStorage (private mode) -- don't nag every load
}

function markSeen() {
  try { localStorage.setItem(SEEN_KEY, '1'); } catch { /* fine, just re-shows next time */ }
}

/** Each step: a title, body content, and an optional CSS selector for the
 * real element this step is anchored to (`.tour-highlight` outline + the
 * popover positions itself next to it). No selector -> centered, for the
 * intro step, which isn't about any one control. Selectors are best-effort:
 * if not found (different viewport, a panel collapsed), the step still
 * shows, centered, rather than anchoring to nothing. */
const STEPS = [
  {
    title: 'Welcome to Vivum',
    selector: null,
    body: () => [
      h('p',
        'This dashboard runs a real physics simulation on a real, measured ' +
        'tumor tissue section — the cell viewer, the mechanism breakdown, ' +
        'and every metric are computed, not illustrative. This tour is ' +
        `${STEPS.length} short steps and points at the real controls as it goes; skip it any time.`),
    ],
  },
  {
    title: '1 · Pick a tissue section',
    selector: '#f-section',
    side: 'bottom',
    body: () => [
      h('p',
        'This ', h('b', 'SECTION'), ' dropdown is not one fixed demo tumor — ' +
        'it currently lists ',
        h('b', '38 real sections across 13 patients'),
        ', spanning glioblastoma (multiple assays, 12 patients) and breast ' +
        'cancer. Switching sections loads genuinely different measured ' +
        'tissue, not a re-skin of the same one.'),
    ],
  },
  {
    title: '2 · Configure and run',
    selector: '#btn-submit-run',
    side: 'top',
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
    side: 'left',
    body: () => [
      h('p',
        h('b', 'This panel'), ' states plainly what each number can and ' +
        'can’t be read as for the run you have open — read it before ' +
        'quoting anything. This project treats an honest "we don’t know" ' +
        'as a valid result, not a gap to hide.'),
    ],
  },
  {
    title: '4 · Script it — MATLAB-style',
    selector: '#btn-protocol',
    side: 'bottom',
    body: () => [
      h('p',
        'This ', h('b', '"protocol"'), ' tab is a real code editor for ' +
        'defining and sweeping drug/dose experiments — 6 working templates ' +
        'to start from, syntax highlighting, inline diagnostics.'),
    ],
  },
];

function clearHighlight() {
  document.querySelectorAll('.tour-highlight').forEach((el) => el.classList.remove('tour-highlight'));
}

function highlight(selector) {
  clearHighlight();
  if (!selector) return null;
  const el = document.querySelector(selector);
  if (el) el.classList.add('tour-highlight');
  return el;
}

/** Position `card` (already in the DOM, so its real measured size is
 * available) next to `targetEl`'s current bounding rect, preferring `side`
 * but flipping or clamping to whatever actually fits in the viewport --
 * a fixed side would run the card off-screen on a small window or a target
 * near an edge. Centers the card (no anchor) when there is no target. */
function positionCard(card, targetEl, side) {
  const cw = card.offsetWidth;
  const ch = card.offsetHeight;
  const vw = window.innerWidth;
  const vh = window.innerHeight;

  if (!targetEl) {
    card.style.top = `${Math.max(MARGIN, (vh - ch) / 2)}px`;
    card.style.left = `${Math.max(MARGIN, (vw - cw) / 2)}px`;
    return;
  }

  const r = targetEl.getBoundingClientRect();
  const fits = {
    bottom: r.bottom + MARGIN + ch <= vh,
    top: r.top - MARGIN - ch >= 0,
    left: r.left - MARGIN - cw >= 0,
    right: r.right + MARGIN + cw <= vw,
  };
  // Preferred side first, then any other side that actually fits, then fall
  // back to the preferred side clamped into the viewport rather than hidden.
  const order = [side, 'bottom', 'top', 'right', 'left'].filter((s, i, a) => a.indexOf(s) === i);
  const chosen = order.find((s) => fits[s]) || side || 'bottom';

  let top, left;
  if (chosen === 'bottom') { top = r.bottom + MARGIN; left = r.left; }
  else if (chosen === 'top') { top = r.top - MARGIN - ch; left = r.left; }
  else if (chosen === 'left') { top = r.top; left = r.left - MARGIN - cw; }
  else { top = r.top; left = r.right + MARGIN; }

  top = Math.min(Math.max(top, MARGIN), vh - ch - MARGIN);
  left = Math.min(Math.max(left, MARGIN), vw - cw - MARGIN);
  card.style.top = `${top}px`;
  card.style.left = `${left}px`;
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
    window.removeEventListener('resize', reposition);
    card.remove();
    document.removeEventListener('keydown', onKey);
  }
  function onKey(ev) {
    if (ev.key === 'Escape') close();
    else if (ev.key === 'ArrowRight') next();
    else if (ev.key === 'ArrowLeft') back();
  }
  function reposition() {
    const step = STEPS[i];
    const el = step.selector ? document.querySelector(step.selector) : null;
    positionCard(card, el, step.side);
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
    const target = highlight(step.selector);
    // Reflow first so offsetWidth/Height below reflect this step's content
    // (body text length differs per step), then anchor to the real element.
    positionCard(card, target, step.side);
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
  const closeX = h('button.btn.ghost.icon.sm', { type: 'button', 'aria-label': 'Close tour', onclick: close }, '×');
  const backBtn = h('button.btn.ghost', { type: 'button', onclick: back }, 'Back');
  const nextBtn = h('button.btn.primary', { type: 'button', onclick: next }, 'Next');
  const skipBtn = h('button.btn.ghost', { type: 'button', onclick: close }, 'Skip tour');

  const card = h('div.tour-popover.notice',
    h('div.tour-head', titleEl, closeX),
    bodyEl,
    dotsEl,
    h('div.btn-row', { style: { marginTop: 'var(--sp-5)', justifyContent: 'space-between' } },
      h('div.btn-row', backBtn, nextBtn),
      skipBtn,
    ),
  );

  document.body.appendChild(card);
  document.addEventListener('keydown', onKey);
  window.addEventListener('resize', reposition);
  render();

  return { node: card, close };
}
