// welcome.js — a one-time, dismissible orientation overlay for a first-time
// visitor. Unlike apikey.js's gate (which BLOCKS everything until resolved,
// because nothing works without a key), this is advisory only: dismissible
// by its own button, the backdrop, or Escape, and never blocks a returning
// or impatient user from just using the tool.
//
// Shown once per browser (localStorage flag), not once per session, so it
// doesn't reappear on every reload during someone's first real work session.

import { h } from '../util/dom.js';

const SEEN_KEY = 'tmesim_welcome_seen_v1';

export function shouldShowWelcome() {
  try { return !localStorage.getItem(SEEN_KEY); }
  catch { return false; } // no localStorage (private mode) -- don't nag every load
}

function markSeen() {
  try { localStorage.setItem(SEEN_KEY, '1'); } catch { /* fine, just re-shows next time */ }
}

/**
 * @param {() => void} onGoToProtocol - called if the user clicks the
 *   "Open the protocol editor" shortcut, in addition to dismissing.
 */
export function createWelcomeOverlay({ onGoToProtocol } = {}) {
  function dismiss() {
    markSeen();
    overlay.remove();
    document.removeEventListener('keydown', onKey);
  }
  function onKey(ev) { if (ev.key === 'Escape') dismiss(); }

  const card = h('div.gate-card.notice', { onclick: (ev) => ev.stopPropagation() },
    h('h3', 'Welcome to Vivum'),
    h('p',
      'This dashboard is running a real physics simulation on a real, ' +
      'measured tumor tissue section — the cell viewer, the mechanism ' +
      'breakdown, and every metric are computed, not illustrative. The ',
      h('b', 'CAVEATS'),
      ' panel on the right states plainly what each number can and can’t ' +
      'be read as; it’s worth reading before quoting anything from a run.'),
    h('p',
      'To script your own drug / dose / sweep experiments — MATLAB-style, ' +
      'not just this form — open the ',
      h('b', '"protocol"'),
      ' tab in the top bar. It starts with 6 working templates.'),
    h('div.btn-row', { style: { marginTop: 'var(--sp-5)' } },
      h('button.btn.primary', {
        onclick: () => { dismiss(); if (onGoToProtocol) onGoToProtocol(); },
      }, 'Open the protocol editor'),
      h('button.btn.ghost', { onclick: dismiss }, 'Explore on my own'),
    ),
  );

  const overlay = h('div.gate-overlay.welcome-overlay', { onclick: dismiss }, card);
  document.addEventListener('keydown', onKey);

  return { node: overlay, dismiss };
}
