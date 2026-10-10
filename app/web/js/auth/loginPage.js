// loginPage.js — /login. The page a person meets instead of a key prompt.
//
// WHAT THIS REPLACES
// ------------------
// A password field labelled `X-API-Key value`, which is a prompt for a shared
// secret. It was the first thing a visitor saw, it could not be answered by
// anyone who had not been handed the operator's secret out of band, and it had
// no concept of a person at all. This page has two modes and neither of them
// asks a human for a machine credential:
//
//   signin   e-mail + password -> POST /api/auth/login   (§1)
//   claim    the first account -> POST /api/auth/signup  (§4, bootstrap token)
//
// §0 IS NOT UNDONE BY ANY OF IT. `X-API-Key` keeps working for every
// non-browser caller -- workflows, render.yaml, the documented curl flows --
// and this page says so out loud rather than leaving a reader to wonder
// whether their scripts just broke. Two credentials, one gate.
//
// THE SENTENCE THIS PAGE EXISTS TO SAY
// ------------------------------------
// §5: a session does not partition anything. The store is shared, every
// account sees the same projects and runs, and nothing is attributed. Every
// other login in the world has taught this reader the opposite, so the page
// says it in the first screenful, above the form and not under it, out of
// js/auth/shared.js so that it cannot drift from what the in-app gate says.
// There is no avatar on this page, no "your projects", no owner anywhere, and
// no "private workspace" -- there is no server behind any of them.
//
// WHAT IT DOES NOT DO
// -------------------
// No reset mail (§7 forbids it and this page states the absence), no OAuth, no
// roles, no invitations. The one route in §1 with no front end here is
// `POST /api/auth/password`: changing a password is a real need but it was not
// in this pass, and a form built against a body shape the contract does not
// pin is how two agents disagree in production.

import { h, fill, clear } from '../util/dom.js';
import * as api from '../api.js';
import {
  MIN_PASSWORD, wayInNote, keyStillWorksNote, noResetNote, sessionLifetimeLine,
} from './shared.js';
import { signOut, signOutWording } from './gate.js';

/* ============================================================ helpers === */

/**
 * Where to go after a successful sign-in.
 *
 * VALIDATED, not trusted. `?next=` arrives in a URL anyone can send, and the
 * naive check -- "does it start with a slash" -- accepts `//evil.example`,
 * which is a protocol-relative URL that the browser follows clean off this
 * origin. A sign-in page that forwards to an attacker's copy of itself is the
 * classic open redirect, and it is worth four lines to not have one.
 */
export function safeNext(search = window.location.search) {
  let raw = '';
  try { raw = new URLSearchParams(search).get('next') || ''; } catch { raw = ''; }
  if (!raw || raw[0] !== '/') return '/';
  if (raw.startsWith('//') || raw.includes('\\')) return '/';
  if (/(^|\/)login(\.html)?\/?($|\?|#)/.test(raw)) return '/';   // no loop
  return raw;
}

/** Which mode the URL asks for. `?mode=claim`, or `#claim` for a link that
 *  does not want to re-request the page. Absent means sign-in. */
function requestedMode(loc = window.location) {
  let m = '';
  try { m = new URLSearchParams(loc.search).get('mode') || ''; } catch { m = ''; }
  if (m === 'claim' || loc.hash === '#claim') return 'claim';
  return 'signin';
}

/** The server's own reason, turned into something a reader can act on.
 *  Status by status, because "it did not work" sends everyone to the same
 *  wrong place. */
function failureText(err, mode) {
  const status = err && typeof err.status === 'number' ? err.status : 0;
  const said = (err && err.message) ? String(err.message) : '';
  if (status === 401 || status === 403) {
    return mode === 'claim'
      ? 'The server refused that token. It is the ' +
        'TMESIM_STUDIO_API_KEY this instance was started with — on a hosted '
        + 'instance the dashboard holds it, and locally the server prints it '
        + 'once at startup next to "api key". ' + said
      : 'That e-mail and password were refused. ' + (said || 'Nothing else is wrong: '
        + 'the server is answering, it simply does not accept that pair.');
  }
  if (status === 409) {
    return 'That e-mail already has an account on this server. Sign in with '
      + 'it instead. ' + said;
  }
  if (status === 429) {
    return 'Too many attempts. This server throttles failed sign-ins per '
      + 'e-mail address, so the next attempt has to wait — the delay is the '
      + 'throttle doing its job, not a fault. ' + said;
  }
  if (status === 400 || status === 422) {
    return said || 'The server rejected the details as malformed.';
  }
  if (status === 404 || status === 405) {
    return 'This server does not serve ' + (mode === 'claim'
      ? '/api/auth/signup' : '/api/auth/login')
      + '. It answered, so it is running — it is running a build without the '
      + 'sign-in routes, and no password will work against it until that is '
      + 'deployed.';
  }
  if (status === 0) {
    return 'The request never reached a server that answered. That is DNS, '
      + 'the network in between, or an instance that is not running. If this '
      + 'instance sleeps when idle, the first request pays the wake and a '
      + 'second attempt often lands immediately.';
  }
  if (status >= 500) {
    return 'The server answered ' + status + ' and the sign-in did not '
      + 'complete. That is a fault on the server, not in what you typed. ' + said;
  }
  return said || 'The sign-in did not complete.';
}

/* =============================================================== page === */

/**
 * Build the page into `root`.
 *
 * Renders SYNCHRONOUSLY first and asks the server afterwards. The form is on
 * screen in the first paint, before `GET /api/auth/me` has answered, because
 * this instance can take 20-35 s to wake and a sign-in page that is blank for
 * half a minute is a sign-in page nobody reaches. The answer, when it comes,
 * can only ADD: it switches to the already-signed-in panel, or preselects the
 * claim mode, or reports that the routes are missing.
 *
 * @param {Element} root
 */
export function mountLoginPage(root) {
  if (!root) throw new Error('mountLoginPage: no root element');

  const next = safeNext();
  let mode = requestedMode();
  let busy = false;
  /** The last thing the server said, as nodes. Never a stale success. */
  let status = null;
  /** `api.resolveSession()`'s answer, or null until it has answered. */
  let me = null;
  /** Why /api/auth/me did not answer, if it did not. */
  let meError = null;

  const host = h('div.auth-host');
  fill(root, host);

  /* ----------------------------------------------------------- fields --- */

  const email = h('input#auth-email.inp', {
    type: 'email', name: 'email', autocomplete: 'username',
    spellcheck: 'false', autocapitalize: 'none', required: true,
    placeholder: 'you@example.org',
  });
  const password = h('input#auth-password.inp', {
    type: 'password', name: 'password', autocomplete: 'current-password',
    required: true, minlength: String(MIN_PASSWORD),
  });
  // The bootstrap token (§4). `type=password` so it is not left on screen in a
  // room with other people in it, and `autocomplete=off` so no browser offers
  // to remember the operator's deployment secret as if it were a password.
  const token = h('input#auth-token.inp', {
    type: 'password', name: 'bootstrap_token', autocomplete: 'off',
    spellcheck: 'false', required: true,
  });

  const submitBtn = h('button#auth-submit.btn.primary', { type: 'submit' }, 'Sign in');
  const statusHost = h('div#auth-status.auth-status', { role: 'alert', 'aria-live': 'assertive' });

  function field(labelText, input, hintNodes) {
    return h('div.auth-field',
      h('label.auth-label', { for: input.id }, labelText),
      input,
      hintNodes ? h('div.auth-hint', hintNodes) : null);
  }

  /* ----------------------------------------------------------- submit --- */

  async function submit(ev) {
    if (ev) ev.preventDefault();
    if (busy) return;

    const e = email.value.trim();
    const p = password.value;
    const t = token.value.trim();

    // Client-side checks are a courtesy to save a round-trip, never a policy:
    // the server decides, and anything this skips is simply sent and refused.
    if (!e) { email.focus(); return say('warn', 'An e-mail address is needed.'); }
    if (!p) { password.focus(); return say('warn', 'A password is needed.'); }
    if (mode === 'claim' && !t) {
      token.focus();
      return say('warn', 'The bootstrap token is needed to create the first account.');
    }
    if (mode === 'claim' && p.length < MIN_PASSWORD) {
      password.focus();
      return say('warn', `The password floor on this server is ${MIN_PASSWORD} `
        + 'characters, so this one would be refused. There are no character '
        + 'rules — length is the whole requirement.');
    }

    busy = true;
    submitBtn.disabled = true;
    const label = submitBtn.textContent;
    submitBtn.textContent = mode === 'claim' ? 'Creating…' : 'Signing in…';
    say('info', mode === 'claim'
      ? 'Creating the account…'
      : 'Checking those details with the server…');

    try {
      const session = mode === 'claim'
        ? await api.signup({ email: e, password: p, token: t })
        : await api.login(e, p);
      // The cookie is set. Leave for `next` rather than rendering a signed-in
      // page the reader did not ask for and would have to click through.
      say('ok', session && session.email
        ? `Signed in as ${session.email}. Opening the workspace…`
        : 'Signed in. Opening the workspace…');
      password.value = '';
      token.value = '';
      window.location.assign(next);
      return undefined;
    } catch (err) {
      say('crit', failureText(err, mode));
      if (mode === 'claim') token.select();
      else password.select();
      return undefined;
    } finally {
      busy = false;
      submitBtn.disabled = false;
      submitBtn.textContent = label;
    }
  }

  function say(kind, text) {
    status = { kind, text };
    paintStatus();
    return undefined;
  }

  function paintStatus() {
    clear(statusHost);
    if (!status) return;
    statusHost.appendChild(h('div.auth-msg.' + status.kind, status.text));
  }

  /* ----------------------------------------------------------- render --- */

  function modeSwitch() {
    const toClaim = mode !== 'claim';
    return h('p.auth-switch',
      toClaim
        ? 'Nobody has claimed this instance yet? '
        : 'This instance already has an account? ',
      h('button.btn.ghost.sm', {
        type: 'button',
        id: 'auth-mode',
        onclick: () => {
          mode = toClaim ? 'claim' : 'signin';
          status = null;
          render();
          (mode === 'claim' ? token : password).focus();
        },
      }, toClaim ? 'Set up the first account' : 'Back to signing in'));
  }

  function signedInPanel() {
    const who = me && me.email ? me.email : null;
    const via = me && me.via ? String(me.via) : null;
    return h('div.auth-card',
      h('h1.auth-h1', 'You are already signed in'),
      who
        ? h('p.auth-who', 'This browser is signed in as ', h('b', who),
            via ? h('span.auth-via', ` (${via})`) : null, '.')
        : h('p.auth-who', 'This browser is signed in. The server did not say '
            + 'which account, so this page will not guess one.'),
      wayInNote(),
      h('p.auth-aside', sessionLifetimeLine()),
      h('div.btn-row',
        h('a.btn.primary', { href: next }, next === '/' ? 'Open the workspace' : 'Continue'),
        h('a.btn.ghost', { href: '/studio' }, 'Open the Studio'),
        h('button#auth-signout.btn.ghost', {
          type: 'button',
          onclick: async () => {
            say('info', 'Signing out…');
            const out = await signOut();
            me = null;
            status = { kind: out.ok ? 'ok' : 'crit', text: signOutWording(out) };
            render();
          },
        }, 'Sign out')),
      statusHost);
  }

  function formPanel() {
    const claim = mode === 'claim';
    const heading = claim ? 'Set up the first account' : 'Sign in';

    submitBtn.textContent = claim ? 'Create the account and sign in' : 'Sign in';
    password.setAttribute('autocomplete', claim ? 'new-password' : 'current-password');

    const form = h('form#auth-form.auth-form', { onsubmit: submit, novalidate: true },
      field('E-mail', email,
        'Stored lower-cased and matched lower-cased, so capitals do not '
        + 'matter. It is an identifier here and nothing else: this server '
        + 'sends no mail to it, ever.'),
      field('Password', password,
        claim
          ? `At least ${MIN_PASSWORD} characters. No required mixture of `
            + 'symbols and digits — those push everyone towards the same '
            + 'handful of passwords, and the length does more.'
          : null),
      claim
        ? field('Bootstrap token', token,
            h('span', 'The ', h('span.mono', 'TMESIM_STUDIO_API_KEY'),
              ' this server runs with. Sign-up is not open: an instance on '
              + 'the public internet with open registration belongs to '
              + 'whoever finds it first, so creating an account requires '
              + 'proof that you run this one. The token is sent to this '
              + 'server and nowhere else, and it is not stored in this '
              + 'browser.'))
        : null,
      h('div.btn-row', submitBtn),
      statusHost);

    return h('div.auth-card',
      h('h1.auth-h1', heading),
      h('p.auth-lede',
        'This deployment asks a person to sign in instead of handing round a '
        + 'shared secret. Nobody has to paste a key into a browser any more.'),
      wayInNote(),
      form,
      modeSwitch(),
      h('div.auth-notes',
        h('p.auth-aside', sessionLifetimeLine()),
        noResetNote(),
        keyStillWorksNote(),
        meNote()),
    );
  }

  /** What the `GET /api/auth/me` probe found, when it is worth saying. It is
   *  never the reason the form is missing: the form is always there. */
  function meNote() {
    if (meError) {
      const s = meError && typeof meError.status === 'number' ? meError.status : 0;
      if (s === 404 || s === 405) {
        return h('p.auth-aside.warn',
          'This server answered, but it does not serve ',
          h('span.mono', '/api/auth/me'), '. It is running a build without the '
          + 'sign-in routes, so the form below will be refused until that '
          + 'build is deployed. Nothing you type is at fault.');
      }
      return h('p.auth-aside.warn',
        'This page could not ask the server who is signed in (',
        h('span.mono', String((meError && meError.message) || 'no reason given')),
        '). The form still works — it is one request and it either succeeds '
        + 'or says why.');
    }
    if (me && me.firstRun === true && mode !== 'claim') {
      return h('p.auth-aside',
        'This server reports that it holds no accounts yet, so signing in '
        + 'cannot succeed until the first one is created with the bootstrap '
        + 'token.');
    }
    if (me && me.firstRun === false && mode === 'claim') {
      return h('p.auth-aside',
        'This server already holds at least one account. The token still '
        + 'works for creating another — sign-up stays token-only, and it does '
        + 'not open up once the first account exists.');
    }
    return null;
  }

  function render() {
    clear(host);
    host.appendChild(me && me.authenticated ? signedInPanel() : formPanel());
    paintStatus();
  }

  render();
  queueMicrotask(() => { try { email.focus(); } catch { /* not focusable yet */ } });

  /* ------------------------------------------------------------- probe --- */

  // §1 guarantees this never 401s, which is exactly what makes it safe to call
  // with no credential at all and before anything else is known.
  (async () => {
    try {
      me = await api.resolveSession();
      meError = null;
      if (me && me.firstRun === true && requestedMode() === 'signin' && !me.authenticated) {
        // An unclaimed instance: the only thing that CAN be done here is the
        // claim, so open on it rather than on a form that cannot succeed.
        mode = 'claim';
      }
    } catch (err) {
      me = null;
      meError = err;
    }
    render();
  })();

  return { node: host, render, get mode() { return mode; } };
}
