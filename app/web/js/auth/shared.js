// shared.js — the words every sign-in surface has to say, in one place.
//
// WHY THESE SENTENCES ARE A MODULE AND NOT TYPED TWICE
// ---------------------------------------------------
// docs/AUTH_CONTRACT.md §5 is the clause this feature is most likely to break
// by accident, and it breaks by OMISSION rather than by error: nobody writes
// "your private workspace", they simply fail to mention that the workspace is
// shared, and the reader supplies the wrong half themselves. Every product
// with a login has taught them that signing in is how you reach your own
// things. Here it is not. Projects, runs and specimens are one store that
// every account sees, and `GET /api/workspace` returns the same workspace to
// everyone.
//
// §5 says the half-built version is the dangerous one, and a UI that implies
// privacy over a shared store is worse than no login at all: it does not just
// fail to inform, it actively teaches a false belief about who can see what.
// So the correction is not a disclaimer at the bottom of one page. It is the
// text these exports hold, rendered on every surface that asks for a
// credential, and asserted on each of them by tests/test_login_portal.py
// against this module rather than against a copy.
//
// WHAT MAY NOT APPEAR, anywhere downstream of here (§5): "your projects", an
// owner column, an avatar, a "private workspace". None of them has a server
// behind it. There is no per-user partitioning to display.

import { h } from '../util/dom.js';

/** The password floor (§3). Twelve characters, and no composition rules --
 *  §3's own reasoning is that those push people toward `Password1!` and the
 *  length floor does more. Mirrored here only so the form can say the rule
 *  before a round-trip; the server is the authority and this never decides
 *  that a password is GOOD, only that it is too short to bother sending. */
export const MIN_PASSWORD = 12;

/** How long a session lasts (§2): 14 days, with no refresh-on-use in v1. */
export const SESSION_DAYS = 14;

/** One sentence, the load-bearing one. Short enough to put in a banner. */
export const WAY_IN_LINE = 'An account here is a way in, not a private space.';

/**
 * The §5 paragraph, as nodes. Used by the sign-in page and by the in-app
 * sign-in route, so the two cannot drift apart.
 *
 * It states the shared store as a fact about this server rather than as a
 * warning, because it is not a defect and nothing is going wrong: one store,
 * deliberately, and the thing worth knowing is simply that it is one store.
 */
export function wayInNote() {
  return h('div.auth-shared',
    h('h2.auth-shared-h', WAY_IN_LINE),
    h('p', 'Signing in does not give you an area of your own. This server '
      + 'keeps ', h('b', 'one'), ' set of projects, runs and registered '
      + 'specimens, and every account sees the same one. Anything you create '
      + 'is visible to everyone who can sign in here, and anything they '
      + 'create is visible to you.'),
    h('p', 'Nothing is attributed, either. The store records no author, so '
      + 'after the fact there is no way to tell from a project or a run which '
      + 'account made it. If that matters for what you are about to do, it '
      + 'matters before you do it rather than after.'),
    h('p.auth-shared-why', 'The sign-in exists so that whoever runs this '
      + 'instance is not paying for everyone who finds the URL. That is the '
      + 'whole of what it is for.'));
}

/**
 * §0, for the reader who already uses the key and is now looking at a
 * password field wondering what happened to it. Nothing happened to it.
 */
export function keyStillWorksNote() {
  return h('p.auth-aside',
    'Scripts and CI do not sign in: they keep sending ',
    h('span.mono', 'X-API-Key'),
    ', exactly as before, and nothing about that has changed or is going to '
    + 'be asked to change. This page is for a person at a keyboard, which is '
    + 'the one caller a shared secret was always a poor fit for.');
}

/** §7, said plainly, because the absence is the surprising part. */
export function noResetNote() {
  return h('p.auth-aside',
    'There is no reset e-mail. This server sends no mail at all — it has no '
    + 'mail configuration and no address for you — so a forgotten password is '
    + 'reset by whoever runs the instance. Half-built e-mail recovery is '
    + 'worse than none, so there is none.');
}

/** The one-line truth about how long the cookie lasts (§2). */
export function sessionLifetimeLine() {
  return `You stay signed in on this browser for ${SESSION_DAYS} days. `
    + 'There is no sliding renewal: the clock starts when you sign in and '
    + 'does not move, so on the fifteenth day you sign in again.';
}
