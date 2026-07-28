// toast.js — transient notices. Errors are NEVER auto-dismissed: a failure
// stays on screen until the operator dismisses it.

import { h, clear } from '../util/dom.js';

let host = null;

function ensure() {
  if (!host) {
    host = h('div.toasts', { role: 'log', 'aria-live': 'polite', 'aria-relevant': 'additions' });
    document.body.appendChild(host);
  }
  return host;
}

/**
 * @param {'info'|'ok'|'warn'|'err'} kind
 * @param {string} title
 * @param {string} [detail]
 * @param {number} [ttl] ms; 0 = sticky. Errors default to sticky.
 */
export function toast(kind, title, detail, ttl) {
  const el = h('div.toast' + (kind === 'info' ? '' : '.' + kind),
    h('div.bar'),
    h('div.body',
      h('b', title),
      detail ? h('code', detail) : null),
    h('button.btn.ghost.icon.sm.x', {
      title: 'Dismiss', 'aria-label': 'Dismiss notification',
      onclick: () => el.remove(),
    }, '×'),
  );
  ensure().appendChild(el);
  const life = ttl === undefined ? (kind === 'err' ? 0 : kind === 'warn' ? 9000 : 4500) : ttl;
  if (life > 0) setTimeout(() => el.remove(), life);
  return el;
}

export const notifyError = (title, err) => {
  const msg = err && (err.message || String(err));
  const st = err && err.status ? ` [HTTP ${err.status}]` : '';
  // eslint-disable-next-line no-console
  console.error(title, err);
  return toast('err', title, (msg || 'unknown error') + st);
};

export function clearToasts() { if (host) clear(host); }
