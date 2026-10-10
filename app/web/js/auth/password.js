// password.js — the front end for `POST /api/auth/password` (contract §1).
//
// WHY THIS DID NOT SHIP WITH THE PORTAL
// -------------------------------------
// `js/auth/loginPage.js` says it, in its own header: the route was the one in
// §1 with no front end, because the contract pins the ROUTE and not the shape
// of its request body, and "a form built against a body shape the contract
// does not pin is how two agents disagree in production". That was the right
// call. The route is now written, so the shape is no longer a guess — it is
// read off `app/server/auth.py::auth_password`, which does exactly this:
//
//     old_password = _text_field(body, "old_password")
//     new_password = _text_field(body, "new_password")
//
// Two fields, both required, both plain strings, and `_text_field` 400s on a
// missing one or on anything that is not a non-empty string. That is the whole
// contract of the request, and `REQUEST_FIELDS` below is the single place this
// module names them, so the pairing is asserted against the server's own
// source rather than against a copy of it (tests/test_account_recovery.py).
//
// SESSION ONLY, AND NOT THE KEY (§1, and the route's own docstring)
// ----------------------------------------------------------------
// This is the one call in the product that must NOT carry `X-API-Key`. The key
// identifies a program, a program has no "own password", and the route refuses
// it on purpose — letting the key set any account's password would make it a
// master credential over identities, which is more than §4 gave it. So the
// fetch below sends the cookie (`credentials: 'same-origin'`) and no key
// header, which is also why it does not go through `api.js::call`: that helper
// attaches the stored key to everything, by design and correctly, for every
// other route in the app.
//
// THE TRAP THIS FILE IS MOST LIKELY TO FALL INTO
// ----------------------------------------------
// A `<form>` whose inputs have `name` attributes and whose submit is not
// prevented does a GET navigation to the current URL with every field value in
// the query string. That would put a live password in the address bar, in the
// browser's history, in the server's access log and in any proxy in between —
// from code that looks completely ordinary. Three independent guards, because
// one is a single edit away from being gone:
//
//   1. `event.preventDefault()` is the FIRST statement of the submit handler,
//      outside every try block, so nothing that can throw runs before it;
//   2. no input carries a `name` attribute, so there is nothing for a form
//      submission to serialise even if one happened;
//   3. the request is built by hand as a POST with a JSON body.
//
// §5 — WHAT THIS SURFACE MAY NOT SAY
// ----------------------------------
// An account is a way in, not a private space. There is no "your password" in
// here (it is "the password for this account"), no avatar, no owner, and
// nothing that implies changing a password changes who can see the workspace.
// It cannot: the store is shared and a password is only the door.

import { h, clear } from '../util/dom.js';
import * as api from '../api.js';
import { MIN_PASSWORD, SESSION_DAYS } from './shared.js';

/** The route (§1). */
export const PASSWORD_PATH = '/api/auth/password';

/** The request body's field names, read off the route (see the header). The
 *  ONE place in the front end that names them. */
export const REQUEST_FIELDS = Object.freeze(['old_password', 'new_password']);

/**
 * `POST /api/auth/password`.
 *
 * Resolves with the server's body on success, which carries the two fields
 * this surface is obliged to report rather than paraphrase:
 * `other_sessions_invalidated` (always `false`) and `other_sessions_note`.
 *
 * Rejects with an `Error` carrying `status`, `message` and `reason` lifted out
 * of the server's `{detail: {message, reason, retry_after_s}}` envelope.
 *
 * WHY IT READS `detail` ITSELF: `api.js::detailOf` handles a `detail` that is
 * a string or an array, which is FastAPI's own two shapes, and returns '' for
 * an object. Every refusal from `auth.py` is an object
 * (`AuthError.payload()`), so routed through that helper the server's sentence
 * is dropped and the reader gets "401 Unauthorized" instead of "your current
 * password is incorrect". That helper is not this agent's file to change, and
 * reading the envelope here is correct in its own right: this module knows
 * exactly which server wrote it.
 */
export async function changePassword(oldPassword, newPassword) {
  const res = await fetch(api.api.base + PASSWORD_PATH, {
    method: 'POST',
    headers: {
      'content-type': 'application/json',
      accept: 'application/json',
    },
    // The session cookie is HttpOnly (§2), so no script here can read it or
    // even tell that it exists -- the browser attaches it and the server
    // decides. `same-origin` and not `include`, matching api.js: the app and
    // the API are one origin, so `include` would buy nothing and would start
    // sending this cookie to a cross-origin base URL.
    credentials: 'same-origin',
    cache: 'no-store',
    // THE SECRETS GO HERE, IN THE BODY, AND NOWHERE ELSE. Never a query
    // string: see the header.
    body: JSON.stringify({
      old_password: String(oldPassword == null ? '' : oldPassword),
      new_password: String(newPassword == null ? '' : newPassword),
    }),
  });

  let payload = null;
  try { payload = await res.json(); } catch { payload = null; }

  if (!res.ok) {
    const detail = (payload && payload.detail) || null;
    const said = detail && typeof detail.message === 'string'
      ? detail.message : '';
    const err = new Error(said || `the server answered ${res.status}`);
    err.status = res.status;
    err.reason = detail && typeof detail.reason === 'string'
      ? detail.reason : '';
    err.retryAfterS = detail && typeof detail.retry_after_s === 'number'
      ? detail.retry_after_s : null;
    throw err;
  }
  return payload || {};
}

/** The server's own reason, as something a reader can act on. Status by
 *  status, the same approach `loginPage.js::failureText` takes, because "it
 *  did not work" sends every cause to the same wrong remedy. */
export function failureText(err) {
  const status = err && typeof err.status === 'number' ? err.status : 0;
  const said = (err && err.message) ? String(err.message) : '';
  const why = (err && err.reason) ? ' ' + String(err.reason) : '';
  if (status === 401) {
    // Two genuinely different 401s on this route: no session at all, and a
    // wrong current password. The server distinguishes them in words and may
    // do so here -- the caller already proved they hold this account's
    // session, so naming the field discloses nothing they did not have. That
    // is the route's own reasoning, and it is why this surface is allowed to
    // be specific where the sign-in page is not.
    return said || 'The server refused that.';
  }
  if (status === 429) {
    return said + ' This server throttles failed attempts per address, so the '
      + 'next one has to wait. The delay is the throttle working, not a fault.';
  }
  if (status === 400 || status === 422) {
    return said + why;
  }
  if (status === 404 || status === 405) {
    return 'This server does not serve ' + PASSWORD_PATH + '. It answered, so '
      + 'it is running — it is running a build without the password route.';
  }
  if (status === 0) {
    return 'The request never reached a server that answered, so nothing was '
      + 'changed. That is the network or an instance that is not running.';
  }
  if (status >= 500) {
    return 'The server answered ' + status + ' and the change did not '
      + 'complete, so the old password still works. ' + said;
  }
  return said || 'The change did not complete.';
}

/** One labelled password field. The id is on the input and the label points at
 *  it, so the control is reachable by its label; `h()` requires the id BEFORE
 *  the classes in a selector and throws on the other order. */
function field(id, label, autocomplete, hint) {
  const input = h('input#' + id + '.inp', {
    type: 'password',
    autocomplete,
    spellcheck: 'false',
    // NO `name` ATTRIBUTE, deliberately -- guard 2 of the three in the
    // header. A name is what a form submission serialises into a URL.
    style: { width: '100%' },
  });
  const wrap = h('div.acct-pwd-field',
    h('label.hint', { for: id }, label),
    input,
    hint ? h('div.hint', hint) : null);
  return { wrap, input };
}

/**
 * The change-password section for the Account panel.
 *
 * Returns one element to append. It owns its own status line and never calls
 * the panel's `render()`: a re-render would rebuild these three inputs, and on
 * success that would wipe the one sentence the reader has to read (what the
 * change did NOT end) at the moment it appears.
 *
 * @returns {Element}
 */
export function passwordChangeSection() {
  const current = field('pwd-current', 'Current password', 'current-password');
  const next = field('pwd-new', 'New password', 'new-password',
    `At least ${MIN_PASSWORD} characters. There are no composition rules — `
    + 'no required digit, no required symbol — because those push people '
    + 'toward Password1! and length does more.');
  const again = field('pwd-again', 'New password again', 'new-password');

  const statusHost = h('div.acct-pwd-status', {
    role: 'status', 'aria-live': 'polite',
  });
  const submit = h('button.btn.primary.sm', { type: 'submit' },
    'Change password');

  function say(kind, text) {
    clear(statusHost);
    statusHost.appendChild(h('div.caveat.' + kind, h('span', text)));
  }

  async function run() {
    const oldPassword = current.input.value;
    const newPassword = next.input.value;
    const repeated = again.input.value;

    // Checked here only to save a round trip on the three mistakes that need
    // no server to diagnose. THE SERVER IS THE AUTHORITY on all of them --
    // `password_problem` in auth.py is the one set of rules, and this never
    // decides that a password is GOOD, only that it is not worth sending.
    if (!oldPassword || !newPassword || !repeated) {
      say('crit', 'All three fields are needed. Nothing was sent.');
      return;
    }
    if (newPassword !== repeated) {
      say('crit', 'The two new entries do not match. Nothing was sent — '
        + 'and nothing has changed, so the current password still works.');
      return;
    }
    if (newPassword.length < MIN_PASSWORD) {
      say('crit', `A password needs at least ${MIN_PASSWORD} characters. `
        + 'Nothing was sent.');
      return;
    }
    if (newPassword === oldPassword) {
      say('crit', 'The new password is the same as the current one. Nothing '
        + 'was sent.');
      return;
    }

    submit.disabled = true;
    say('info', 'Changing it…');
    try {
      const out = await changePassword(oldPassword, newPassword);
      current.input.value = '';
      next.input.value = '';
      again.input.value = '';
      clear(statusHost);
      statusHost.appendChild(h('div.caveat.info',
        h('span',
          h('b', 'PASSWORD CHANGED, AND NOTHING WAS SIGNED OUT. '),
          // THE SERVER'S OWN SENTENCE, RENDERED VERBATIM. It is on the
          // response as `other_sessions_note` precisely so this surface does
          // not re-type it and drift; the fallback is only for a server that
          // predates the field. Reporting this as "signed out everywhere"
          // would be the comfortable lie the route refuses to tell.
          (typeof out.other_sessions_note === 'string'
            && out.other_sessions_note)
            ? out.other_sessions_note
            : 'Sessions here are signed tokens with no server-side record, so '
              + 'changing a password cannot end the ones already issued. They '
              + `expire on their own within ${SESSION_DAYS} days.`)));
    } catch (err) {
      say('crit', failureText(err));
    } finally {
      submit.disabled = false;
    }
  }

  const form = h('form.acct-form', {
    // No `action` and no `method`: this never navigates. See guard 1 in the
    // header -- preventDefault() is the first statement and is outside the
    // try, so a throw anywhere in `run()` cannot let the browser fall back to
    // submitting these three fields as a GET query string.
    onsubmit: (event) => {
      event.preventDefault();
      run();
      return false;
    },
  },
    current.wrap, next.wrap, again.wrap,
    h('div.btn-row', submit),
    statusHost);

  const wrap = h('div.acct-pwd');
  wrap.appendChild(h('div.hr-lbl', 'Change the password for this account'));
  wrap.appendChild(h('div.hint',
    'This changes how you sign in. It changes nothing about what you can '
    + 'see: the projects, runs and registered specimens are one shared set '
    + 'and every account reaches the same one.'));
  wrap.appendChild(form);
  return wrap;
}
