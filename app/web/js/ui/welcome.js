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
    title: 'Welcome to Codon',
    selector: null,
    body: () => [
      h('p',
        'This dashboard runs a real physics simulation on a real, measured ' +
        'tumor tissue section — the cell viewer, the mechanism breakdown, ' +
        'and every metric are computed, not illustrative. This tour is ' +
        `${STEPS.length} short steps and points at the real controls as it goes; skip it any time.`),
      // The tour was the one first-run surface in the Studio that stated no
      // limit at all — it introduced a computed tumour response to someone who
      // had never seen the product and handed them straight to Submit run. The
      // claim below is the project's own standing limit (see
      // ui/provenance.js::STANDING_LIMITS) and it belongs in the first thing a
      // new user reads, not only in the panel they reach afterwards.
      h('p',
        h('b', 'Nothing here is a clinical or treatment-planning tool. '),
        'No number in this Studio has been validated against a patient ' +
        'outcome: these sections carry no treatment-response label, so every ' +
        'result is a prediction of the physical model that has never been ' +
        'scored against one. The ranking across compounds is the defensible ' +
        'output; the absolute numbers are not. NOT FOR CLINICAL USE.'),
    ],
  },
  {
    title: '1 · Pick a tissue section',
    selector: '#f-section',
    side: 'bottom',
    body: () => [
      /*
       * DO NOT PUT A SECTION COUNT BACK IN THIS SENTENCE.
       *
       * It said "38 real sections across 13 patients ... spanning
       * glioblastoma and breast cancer". That was true on the machine it was
       * written on and false on a fresh install, which is the only machine
       * whose reader is being addressed here. Two independent reasons:
       *
       *   (1) `app/runs/` is gitignored, so the ingest registry does not ship.
       *       A clean checkout has SEVEN built-in sections -- all glioblastoma
       *       (Xenium 319-panel, Xenium IO 480-panel, GSE324268, and four
       *       Visium). There is no built-in breast section. The other 31 and
       *       the breast tissue were locally ingested by the developer.
       *   (2) The count is now user-mutable: the "Your tissue" panel changes
       *       it, so any literal is stale the moment a field tester ingests.
       *
       * A first-run tour is the one surface a new user cannot check against
       * anything else, so a wrong number here is caught immediately and costs
       * the trust every other disclosure depends on. Describe the KIND of
       * thing in the list and let the list speak for its own length.
       */
      h('p',
        'This ', h('b', 'SECTION'), ' dropdown is not one fixed demo tumor. ' +
        'It lists the built-in specimens this project ships measured evidence ' +
        'for — real glioblastoma tissue across two assay families, Xenium and ' +
        'Visium — alongside any section you ingest yourself. Switching ' +
        'sections loads genuinely different measured tissue, not a re-skin of ' +
        'the same one.'),
    ],
  },
  {
    title: '2 · Or bring your own tissue',
    selector: '#sect-tissue',
    side: 'right',
    body: () => [
      h('p',
        'The list is not the limit. ', h('b', 'Your tissue'),
        ' takes a 10x Xenium, Space Ranger (Visium) or MEX export from a ' +
        'directory on the machine running this server — nothing is uploaded ' +
        'from the browser — and tells you what that specimen can and ' +
        'cannot support before anything is run.'),
      h('p',
        'Two things it will say to you, and both are load-bearing: a section ' +
        'you ingest is ', h('b', 'not validated'), ' — the built-in specimens ' +
        'are the ones this project ships measured evidence for — and Codon ' +
        'may be pointed at ', h('b', 'de-identified tissue only'), '.'),
    ],
  },
  {
    title: '3 · Configure and run',
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
    title: '4 · Read the CAVEATS panel',
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
    title: '5 · Script it — MATLAB-style',
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
    if (ev.key === 'Escape') { close(); return; }
    // Arrow keys advance the tour ONLY when the operator is not typing. This
    // listener is on `document`, and the tour is deliberately non-modal (see
    // the module comment) — so the rest of the Studio stays live underneath
    // it, including the RunSpec editor and every form field. Without this
    // guard, moving the caret with ArrowLeft inside the JSON editor stepped
    // the tour backwards, which reads as the tour fighting the user for the
    // keyboard. Escape is exempt: closing the tour is what Escape should do
    // from anywhere, and the fields that care already stopPropagation on it.
    const t = ev.target;
    const typing = t && (t.isContentEditable
      || /^(INPUT|TEXTAREA|SELECT)$/.test(t.tagName || ''));
    if (typing || ev.metaKey || ev.ctrlKey || ev.altKey) return;
    if (ev.key === 'ArrowRight') next();
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

  // `role="dialog"` with aria-modal FALSE, which is the honest description:
  // this popover deliberately does not block the app behind it (that is the
  // whole point of anchoring it to real controls), so claiming modality would
  // tell a screen reader the rest of the page is inert when it is not.
  // `aria-labelledby` points at the step title, so each step announces itself.
  titleEl.id = 'tour-title';
  const card = h('div.tour-popover.notice', {
    role: 'dialog', 'aria-modal': 'false', 'aria-labelledby': 'tour-title',
  },
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
  // Land the keyboard on the tour once, at creation — not on every render,
  // which would yank focus back on each step and make Back/Next unusable from
  // the keyboard. Without this a keyboard user arriving at a first run had to
  // tab past the entire application (445 focusable elements were counted on
  // this workstation) before reaching `Next` or `Skip tour`.
  queueMicrotask(() => nextBtn.focus());

  return { node: card, close };
}
