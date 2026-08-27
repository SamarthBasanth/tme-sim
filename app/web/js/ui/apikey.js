// apikey.js — the access-key gate shown when the backend requires
// TMESIM_STUDIO_API_KEY and this browser doesn't have a working one stored.
//
// Distinct from the WebGL-unavailable notice in main.js: that one deliberately
// does NOT cover the rest of the shell (the run control / history / provenance
// panels still work without a viewport). This one does, because nothing in
// Studio is meaningfully usable — every panel is backed by the same gated API
// — until a valid key is entered. It is the only overlay in this app that
// covers the whole viewport rather than sitting inside the canvas.

import { h } from '../util/dom.js';
import { setApiKey } from '../api.js';

/**
 * @param {(key: string, remember: boolean) => Promise<boolean>} onSubmit - try the key against
 *   the real backend; resolve true if it worked (the gate then removes
 *   itself), false if the server still rejects it (shows an inline error).
 */
export function createApiKeyGate({ onSubmit }) {
  const input = h('input.inp', {
    type: 'password',
    autocomplete: 'off',
    spellcheck: 'false',
    placeholder: 'X-API-Key value',
    style: { width: '100%', marginTop: 'var(--sp-4)' },
  });
  // Default ON. The key was already persisted unconditionally before this
  // existed, so the checkbox does not change the default behaviour -- it makes
  // the behaviour VISIBLE and gives a shared or public machine a way out.
  const remember = h('input', { type: 'checkbox', checked: true });
  const error = h('div.gate-error', { style: { display: 'none' } });
  const submitBtn = h('button.btn.primary', { type: 'submit' }, 'Connect');

  async function submit(ev) {
    ev.preventDefault();
    const key = input.value.trim();
    if (!key) { input.focus(); return; }
    submitBtn.disabled = true;
    submitBtn.textContent = 'Connecting…';
    error.style.display = 'none';
    try {
      const ok = await onSubmit(key, remember.checked);
      if (!ok) {
        error.textContent = 'That key was rejected by the server. Check it and try again.';
        error.style.display = '';
        input.select();
      }
      // On success the caller removes this gate; nothing more to do here.
    } finally {
      submitBtn.disabled = false;
      submitBtn.textContent = 'Connect';
    }
  }

  const form = h('form.gate-card.notice', { onsubmit: submit },
    h('h3', 'Access key required'),
    h('p', 'This is a private deployment. Every /api/ request needs a key, '
      + 'so that whoever runs this instance is not paying for anyone who '
      + 'finds the URL.'),
    h('p', { style: { color: 'var(--fg-2)' } },
      'If someone shared this Studio with you, ask them for the key — there '
      + 'is no self-service sign-up. If you are running it yourself, it is '
      + 'the TMESIM_STUDIO_API_KEY you set on the host; with none set, the '
      + 'server generates one and prints it once at startup next to '
      + '"api key".'),
    input,
    h('label.gate-remember', { style: { display: 'flex', gap: 'var(--sp-2)',
                                        alignItems: 'center',
                                        marginTop: 'var(--sp-3)',
                                        fontSize: 'var(--fs-xs)',
                                        color: 'var(--fg-2)' } },
      remember, 'Stay signed in on this browser'),
    error,
    h('div.btn-row', { style: { marginTop: 'var(--sp-5)' } }, submitBtn),
  );

  const overlay = h('div.gate-overlay', form);

  // Submitting via Enter in the input works through the form's own submit
  // event above; no separate keydown handler needed.
  queueMicrotask(() => input.focus());

  return {
    node: overlay,
    remove: () => overlay.remove(),
  };
}

/**
 * Wires the gate to api.js's setApiKey + a real probe against the backend.
 * `probeFn` is api.probe (passed in rather than imported directly so this
 * module has no circular import on api.js's own use of it).
 */
export async function ensureAuthenticated(probeFn) {
  const first = await probeFn();
  if (!first.authRequired) return; // no key needed, or genuinely unreachable (stub mode handles that)

  await new Promise((resolve) => {
    const gate = createApiKeyGate({
      onSubmit: async (key, remember) => {
        setApiKey(key, remember !== false);
        const retry = await probeFn();
        if (retry.authRequired) {
          // Confirmed wrong by the server. Clear it from BOTH stores so a
          // reload does not silently retry a key already known to be rejected.
          setApiKey('');
          return false;
        }
        gate.remove();
        resolve();
        return true;
      },
    });
    document.body.appendChild(gate.node);
  });
}
