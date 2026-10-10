// main.js — boots the bring-your-own-tissue page.
//
// Two panes and one rule between them: the LEFT pane (browser.js) owns which
// directory is in view, and the RIGHT pane (report.js) renders whatever the
// left pane last reported. The right pane never fetches a listing and the left
// pane never judges a sample, so "which directory" has exactly one owner and
// the two can never be looking at different places.
//
// The access-key gate is here rather than in either pane because either one
// can be the first to meet a 401, and a page with two key prompts on it is a
// page that teaches people to ignore both.

import { h, fill, $ } from '../util/dom.js';
import { createBrowser } from './browser.js';
import { createReport } from './report.js';
import { net, setApiKey } from './net.js';

function createKeyGate({ onKey }) {
  const input = h('input#kg-key.inp', {
    type: 'password', spellcheck: 'false', autocomplete: 'off',
    'aria-label': 'Access key',
  });
  const persist = h('input#kg-persist', { type: 'checkbox', checked: true });
  const go = h('button#kg-go.btn.primary.sm', { type: 'button' }, 'Use this key');
  const msg = h('p#kg-msg.kg-msg');

  const el = h('div#kg.kg', { hidden: true },
    h('h2.kg-title', 'This server needs an access key'),
    msg,
    h('div.kg-row', input, go),
    h('label.kg-check', persist, 'stay signed in on this browser'),
    h('p.kg-hint',
      'The key is whatever the server printed when it started, or whatever '
      + 'was set in its environment. It is the same key the Studio asks for, '
      + 'and entering it here signs you in there too.'));

  function submit() {
    const v = input.value.trim();
    if (!v) return;
    setApiKey(v, persist.checked);
    input.value = '';
    el.hidden = true;
    if (onKey) onKey();
  }
  go.addEventListener('click', submit);
  input.addEventListener('keydown', (ev) => {
    if (ev.key === 'Enter') { ev.preventDefault(); submit(); }
  });

  return {
    el,
    show(message) {
      fill(msg, message || '');
      el.hidden = false;
      input.focus();
    },
  };
}

function boot() {
  const browseHost = $('#tissue-browse');
  const reportHost = $('#tissue-report');
  const gateHost = $('#tissue-gate');
  if (!browseHost || !reportHost || !gateHost) return;

  const report = createReport({
    onRegistered: () => browser.reload(),
  });

  const gate = createKeyGate({ onKey: () => browser.start() });

  const browser = createBrowser({
    onPlace: (state) => report.show(state),
    onRefusal: (err) => {
      // Only the key is the page's problem. A 403 from the allowlist is an
      // answer about the path and belongs where the path is, which is the
      // browser pane's own status line — it already renders the server's
      // sentence there.
      if (err.status === 401) gate.show(err.message);
    },
  });

  fill(gateHost, gate.el);
  fill(browseHost, browser.el);
  fill(reportHost, report.el);

  // A key already in storage (from the Studio, or from a previous visit) is
  // used without asking, and with no key at all the first call still goes
  // out: the loopback default has no key, so prompting before the server has
  // said it wants one is a gate in front of an open door. `net` is imported
  // for that contract and read by the gate's own retry path.
  void net;
  browser.start();
}

if (document.readyState === 'loading') {
  document.addEventListener('DOMContentLoaded', boot, { once: true });
} else {
  boot();
}
