// gate.js — what the app does when the server says 401.
//
// It used to take a password field. `ui/apikey.js` put up a full-viewport card
// whose input was labelled `X-API-Key value` and asked a visitor to paste a
// shared secret, and that is what the product owner was objecting to: it is a
// credential prompt for a machine, shown to a person. This module replaces it
// with a route to /login (docs/AUTH_CONTRACT.md §1). It contains no password
// field, no key field, and no way to enter either -- a reader who meets this
// gate is sent to the sign-in page and never asked for a secret in place.
//
// THE KEY PATH IS NOT TOUCHED. §0 and §7 both forbid weakening it, and nothing
// here weakens it: a browser that holds a working key never produces the 401
// that brings this module into play, so an operator who only ever uses the key
// does not see a login page at all. That is §0's compatibility bar, and it is
// met by arithmetic rather than by care -- there is no code path from "the key
// worked" to anything in this file.
//
// TWO SHAPES, AND WHY THE PAGES DO NOT AGREE
// ------------------------------------------
//   ensureAuthenticated()  the Studio. Nothing in the Studio is usable without
//                          a credential -- every panel reads the same gated
//                          API -- so the reader is NAVIGATED to /login. The
//                          overlay it paints first is not a dialog to dismiss;
//                          it is what is on screen during the navigation, and
//                          it is also the only thing a reader sees if the
//                          navigation is blocked.
//   signInRoute()          the workspace home page. It is the landing page: it
//                          has a working header, a footer, and honest copy
//                          about what this deployment is, and bouncing the
//                          root URL straight to a login form throws all of
//                          that away before anyone has seen the product. So
//                          the home page keeps its chrome and shows a sign-in
//                          route inside its own status region.
// The asymmetry is the point rather than an inconsistency: one page has
// something to show a signed-out reader and the other does not.

import { h } from '../util/dom.js';
import * as api from '../api.js';
import { WAY_IN_LINE, sessionLifetimeLine } from './shared.js';

/** Where to come back to after signing in. */
function hereHref() {
  return window.location.pathname + window.location.search;
}

/** The primary action, every time this module offers one. An `<a href>` and
 *  deliberately not a button with a click handler: it is a navigation, so it
 *  should be middle-clickable, keyboard-activatable and visible in the status
 *  bar before it is followed, and it must still work if this module's own JS
 *  has gone wrong somewhere after this node was built. */
export function signInLink(label = 'Sign in', extraClass = 'primary') {
  return h('a.btn.' + extraClass, {
    href: api.signInHref(hereHref()),
    rel: 'nofollow',
  }, label);
}

/** The sentence that explains a refusal, given what credential was on the
 *  wire. Three genuinely different situations, and conflating them is how a
 *  reader with a dead key spends an afternoon re-typing a good password. */
export function refusalLines({ keyRejected = false } = {}) {
  const lines = [
    h('p', 'The server is there and answering. It refused this request '
      + 'because nothing it accepts was attached to it — so nothing here can '
      + 'be read until somebody signs in.'),
  ];
  if (keyRejected) {
    lines.push(h('p',
      h('b', 'This browser also sent a stored '), h('span.mono', 'X-API-Key'),
      h('b', ' and the server refused that too. '),
      'So the key in this browser is wrong, revoked, or for a different '
      + 'deployment. It has not been deleted — it is still stored here and '
      + 'will be sent again — but it is not what gets you in. Signing in '
      + 'does.'));
  }
  lines.push(h('p.auth-aside.muted', WAY_IN_LINE
    + ' This server keeps one shared set of projects and runs, and every '
    + 'account sees the same one.'));
  return lines;
}

/**
 * The home page's replacement for the access-key banner: nodes to drop into
 * whatever container that page already uses for its status region, plus the
 * action row. No field of any kind.
 *
 * @returns {{body: Element[], actions: Element[]}}
 */
export function signInRoute({ keyRejected = false } = {}) {
  return {
    body: refusalLines({ keyRejected }),
    actions: [signInLink('Sign in'), h('a.btn.ghost', { href: '/api/health' }, 'Check server health')],
  };
}

/**
 * The Studio's gate. Same name and same signature as the `ensureAuthenticated`
 * it replaces in `ui/apikey.js`, so the call site in `js/main.js` did not have
 * to learn anything new.
 *
 * Resolves when the reader is already allowed through. When they are not it
 * paints the overlay, starts the navigation, and returns a promise that NEVER
 * settles -- which is not a leak and not a hang: boot must not continue
 * painting a Studio over a page that is in the middle of leaving, and the one
 * thing that outlives this promise is the page itself. `ui/apikey.js` did the
 * same thing with an awaited promise that only a correct key resolved.
 *
 * @param {() => Promise<{authRequired?: boolean}>} probeFn - api.probe
 */
export async function ensureAuthenticated(probeFn) {
  const first = await probeFn();
  // Not refused: either no credential is needed (the ordinary local case) or
  // the backend is genuinely unreachable, which is stub mode's problem and
  // not this module's. Unchanged from the behaviour this replaces.
  if (!first || !first.authRequired) return;

  const overlay = signInOverlay({ keyRejected: !!(api.api && api.api.keyRejected) });
  document.body.appendChild(overlay);

  // If the navigation is refused (a test harness, an embedded view), the
  // overlay above is the whole experience and it carries a working link, so
  // there is still exactly one way forward and it is visible.
  api.requireSignIn({ next: hereHref() });
  return new Promise(() => {});
}

/** The full-viewport card. Reuses `.gate-overlay` / `.gate-card` / `.notice`
 *  from css/components.css -- the same classes the key gate used, because the
 *  SHAPE of this screen (one thing, covering everything, keyboard-reachable)
 *  did not change even though its content did. Only the new `.auth-*` classes
 *  inside it come from css/auth.css. */
export function signInOverlay({ keyRejected = false } = {}) {
  return h('div.gate-overlay',
    h('div.gate-card.notice',
      h('h3', 'Sign in to open the Studio'),
      ...refusalLines({ keyRejected }),
      h('p.auth-aside.muted', sessionLifetimeLine()),
      h('div.btn-row', { style: { marginTop: 'var(--sp-5)' } },
        signInLink('Go to sign-in'),
        h('a.btn.ghost', { href: '/' }, 'Back to the workspace'))));
}

/**
 * End the session, from anywhere.
 *
 * Returns what actually happened rather than a boolean, because the caller has
 * to be able to SAY it: on a browser that also held a key, signing out of the
 * session alone would leave the reader able to read everything while being
 * told they are signed out, so `api.logout` forgets the key too by default
 * (its own comment argues the trade). A UI that does not report the second
 * half is back to implying something false, which is the failure this whole
 * contract is written against.
 */
export async function signOut() {
  try {
    const out = await api.logout();
    return { ok: true, ...out };
  } catch (err) {
    return { ok: false, error: err, sessionCleared: false, keyForgotten: false };
  }
}

/** What to say after signOut(), as one sentence, matching what it did. */
export function signOutWording(out) {
  if (!out || !out.ok) {
    return 'The server did not confirm the sign-out, so treat this browser as '
      + 'still signed in until a reload says otherwise.';
  }
  return out.keyForgotten
    ? 'The session is closed and the stored X-API-Key was removed from this '
      + 'browser too — both credentials are gone, so nothing here can read '
      + 'this server now.'
    : 'The session is closed. Nothing is stored in this browser any more.';
}
