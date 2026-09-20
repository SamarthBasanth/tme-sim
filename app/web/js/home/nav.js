// nav.js — the top bar shared by the workspace home (`/`) and the Studio
// (`/studio`).
//
// docs/WORKSPACE_CONTRACT.md §5 puts this file in js/home/, but nothing in it
// is home-specific: it is the CONTEXT zone of docs/UI_ARCHITECTURE.md §2.1,
// rendered by ONE module instead of by whichever page happened to need it
// first. That rule is already written down and was already broken once —
// "three different modules injected buttons into it independently and one of
// them did so from a dead code path" (§2.1) — so this module is deliberately
// ADOPTIVE rather than assertive:
//
//   * on the Studio, whose header js/main.js already owns, it adds only the
//     one control the contract adds (§8: "a link back to /") plus the project
//     context, and touches nothing else;
//   * it never builds a second theme button on a page that already has
//     `#btn-theme`;
//   * it never builds a second account button on a page that already has
//     `#btn-account` (js/ui/account.js injects its own);
//   * a second `createNav()` on the same page returns the first one.
//
// THE HOME PAGE MUST NOT RENDER ITS OWN `.topbar`. On `/` this module builds
// the bar and mounts it; a hand-written header in home.html would be a second
// owner of the one zone the architecture says has exactly one.
//
// WHAT THIS FILE IS NOT ALLOWED TO BECOME
// ---------------------------------------
// A sign-in bar. WORKSPACE_CONTRACT.md §9 and `GET /api/account` agree: one
// shared key, no identities, no attribution. There is therefore no avatar, no
// name, no "your workspace" and no owner column anywhere in here. The access
// control that DOES exist is a key this browser sends, and the entry point for
// it is js/ui/account.js — reused, not reimplemented (§9 forbids a new key
// mechanism).

import { h, $ } from '../util/dom.js';
import * as api from '../api.js';
import { createAccountPanel } from '../ui/account.js';

/* ============================================================== theme ==== */

/**
 * The storage key the Studio already writes.
 *
 * GREPPED, not guessed — js/main.js:29 is `const THEME_KEY = 'tmesim.theme';`
 * and js/main.js:42 is `localStorage.setItem(THEME_KEY, t)`. Contract §7 ("the
 * two pages share localStorage, so a user who picked light in the Studio lands
 * on a light home page") is satisfied by this one string being identical and by
 * nothing else, so it is exported: anyone who needs it imports it from here
 * rather than retyping the literal a third time.
 */
export const THEME_KEY = 'tmesim.theme';

const themeSubs = new Set();

/** 'light' | 'dark' — what the document is actually wearing right now. */
export function currentTheme() {
  return document.documentElement.getAttribute('data-theme') === 'light'
    ? 'light' : 'dark';
}

/**
 * Set the theme, persist it, and repaint every toggle that names it.
 *
 * Deliberately mirrors js/main.js `applyTheme` down to the button label, so the
 * two pages' toggles are the same control rather than two controls that happen
 * to do the same thing (UI_ARCHITECTURE.md §5). If both modules are ever live
 * on one page — they are not today, but §2.1's history says assume they will
 * be — the `#btn-theme` lookup below finds the SAME element main.js updates, so
 * neither can leave the other's label stale.
 */
export function applyTheme(t) {
  const theme = t === 'light' ? 'light' : 'dark';
  document.documentElement.setAttribute('data-theme', theme);
  try { localStorage.setItem(THEME_KEY, theme); } catch { /* private mode */ }
  for (const b of document.querySelectorAll('#btn-theme, [data-theme-toggle]')) {
    b.textContent = 'theme: ' + theme;
    b.setAttribute('aria-pressed', String(theme === 'light'));
  }
  for (const fn of themeSubs) {
    // One bad subscriber must not stop the others repainting: a half-themed
    // page is worse than a logged error.
    try { fn(theme); } catch (e) { console.error('theme subscriber failed', e); }
  }
  return theme;
}

/** Flip it. */
export function toggleTheme() {
  return applyTheme(currentTheme() === 'light' ? 'dark' : 'light');
}

/**
 * Adopt the stored theme, falling back to the OS preference.
 *
 * Same precedence as js/main.js `initTheme`, so landing on `/` and landing on
 * `/studio` cannot disagree about what "no stored preference" means.
 */
export function initTheme() {
  let t = null;
  try { t = localStorage.getItem(THEME_KEY); } catch { /* private mode */ }
  if (t !== 'light' && t !== 'dark') {
    t = (window.matchMedia && window.matchMedia('(prefers-color-scheme: light)').matches)
      ? 'light' : 'dark';
  }
  return applyTheme(t);
}

/**
 * Run `fn(theme)` after every flip. Returns an unsubscribe.
 *
 * For anything that CACHES a colour rather than inheriting it — a project
 * card's composition bar drawn into a canvas is that case on the home page,
 * exactly as the viewer and the mechanism legend are in the Studio. Those two
 * drifted apart once (js/main.js:51-55: the click path repainted the legend and
 * the key path did not), which is why this is a subscription and not a second
 * copy of somebody's repaint list.
 */
export function onThemeChange(fn) {
  themeSubs.add(fn);
  return () => themeSubs.delete(fn);
}

/* ================================================================ nav ==== */

/** `.hidden` from css/base.css:140 is `display: none !important`, which is the
 *  only form that survives `.topbar > * { display: flex }` in css/layout.css.
 *  The bare `hidden` ATTRIBUTE does not: the UA's `[hidden] {display:none}` is
 *  a single-attribute selector and loses to that class rule, so a group hidden
 *  by the attribute alone stays on screen. Both are set — the class to hide it,
 *  the attribute so it leaves the accessibility tree too. */
function setHidden(el, on) {
  el.classList.toggle('hidden', !!on);
  if (on) el.setAttribute('hidden', '');
  else el.removeAttribute('hidden');
}

/**
 * Build (or adopt) the top bar. Idempotent per page.
 *
 * @param {object} [o]
 * @param {'home'|'studio'} [o.page]  which page this is. Decides the product
 *   mark's second word, which link carries `aria-current="page"`, and whether
 *   an existing `.topbar` is adopted rather than replaced.
 * @param {Element} [o.host]          js/home/main.js's name for the mount
 *   point; it calls `createNav({page, host, onThemeToggle})` and then appends
 *   the returned `.el` into that same host, so both halves of the handshake
 *   are honoured — this mounts itself into `host` AND returns `{el}`.
 * @param {Element} [o.mount]         same thing under this module's own name.
 *   Default: first child of `<body>`.
 * @param {Function} [o.onThemeToggle] the host page's own theme flip. Used in
 *   place of this module's when supplied, so the `data-theme` stamp keeps a
 *   single owner; the bar then re-asserts the result to relabel every toggle
 *   and notify subscribers, which is idempotent.
 * @param {string} [o.version]        the version mark. Default 'v1' — what
 *   index.html carries and what the contract calls the server contract.
 * @param {{id: string, name?: string}} [o.project]  project context, or null.
 * @param {boolean} [o.account]       wire js/ui/account.js's entry point.
 *   Default true.
 * @param {boolean} [o.watchBackend]  keep the LED in sync with api.js's
 *   stub/live mode. Default true, and it refuses to steal the hook if another
 *   module already owns it.
 * @returns {{el: Element, setProject: Function, setBackend: Function,
 *            themeButton: Element|null, accountPanel: object|null,
 *            destroy: Function}}
 */
export function createNav(o = {}) {
  const page = o.page === 'studio' ? 'studio' : 'home';
  const version = o.version || 'v1';
  const homeHref = o.homeHref || '/';
  const studioHref = o.studioHref || '/studio';

  // `host` is js/home/main.js's name for it, `mount` is this module's. A host
  // INSIDE an existing banner (home.html mounts at `#home-nav-mount`, which
  // sits inside `<header class="topbar home-topbar">`) changes what gets
  // built: a second `.topbar` nested in the first would paint a second
  // background and a second bottom border, and a second `role="banner"` would
  // give the page two banner landmarks. So in that case the bar is a plain
  // flex row that carries only its contents, and the host keeps being the bar.
  const host = o.host || o.mount || null;
  const nested = !!(host && typeof host.closest === 'function'
    && host.closest('.topbar, header, [role="banner"]'));

  let project = o.project && o.project.id ? o.project : null;

  /* --- already here? ---------------------------------------------------- */

  const mounted = $('.nav-bar');
  if (mounted && mounted.__navController) return mounted.__navController;

  // The Studio's own header. It already carries the mark, the backend LED, the
  // theme toggle and the launcher group.
  if (page === 'studio') {
    const bar = $('.topbar');
    if (bar) return adoptExistingBar(bar, { homeHref, project });
  }

  /* --- brand ------------------------------------------------------------ */

  const brandEl = h('div.brand.nav-brand',
    h('span.name', h('b', 'Vivum'), ' ',
      h('span', page === 'studio' ? 'Studio' : 'Workspace')),
    h('span.ver', version));

  /* --- the link between the two pages ----------------------------------- */
  // Both links are always rendered, and the current page's is still a real
  // link rather than a disabled button: `aria-current="page"` is what says "you
  // are here". A control present on one page and absent on the other is the
  // thing that makes two pages feel like two products.

  const homeLink = h('a.btn.ghost.sm.nav-link', {
    href: homeHref,
    title: 'Workspace: every project on this server, and the runs across them',
  }, 'workspace');
  const studioLink = h('a.btn.ghost.sm.nav-link', {
    href: studioHref,
    title: 'Studio: reconstruct a specimen, run a therapy on it, read the mechanism',
  }, 'studio');
  const here = page === 'home' ? homeLink : studioLink;
  here.setAttribute('aria-current', 'page');
  here.classList.add('nav-link-on');

  const linksEl = h('nav.topbar-group.keep.nav-links',
    { 'aria-label': 'Pages' }, homeLink, studioLink);

  /* --- project context (contract §8) ------------------------------------ */

  const projIdEl = h('span.mono.strong.nav-project-id');
  const projNameEl = h('span.dim.nav-project-name');
  const projEl = h('div.topbar-group.nav-project',
    h('span.lbl', 'Project'), projIdEl, projNameEl);

  /* --- backend state ----------------------------------------------------- */
  // Same two elements, same two words, same colours as the Studio's — see
  // js/main.js `paintBackendState`. The ids are claimed only if free, so this
  // can never shadow the Studio's own.

  const ledEl = h('span.led.nav-led', { 'aria-hidden': 'true' });
  const backendEl = h('span.mono.nav-backend',
    { title: 'Backend connection state' }, 'probing…');
  if (!document.getElementById('hdr-led')) ledEl.id = 'hdr-led';
  if (!document.getElementById('hdr-backend')) backendEl.id = 'hdr-backend';
  const statusEl = h('div.topbar-group.keep.nav-status', ledEl, backendEl);

  /* --- actions ----------------------------------------------------------- */

  let themeButton = null;
  const actionsEl = h('div.topbar-group.last.nav-actions');
  // If nothing has set the theme yet, set it here rather than let the bar paint
  // a toggle for a state the page is not in.
  if (!document.documentElement.hasAttribute('data-theme')) initTheme();
  if (!document.getElementById('btn-theme')) {
    // Built from the LIVE theme, never from a hardcoded default. Written as a
    // constant first and caught by a browser test: `initTheme()` legitimately
    // runs before `createNav()` — main.js does exactly that — so a button born
    // saying "theme: dark" on a light page is not an edge case, it is the
    // ordinary boot order, and it is the same stale-label drift §2.1 records.
    const t0 = currentTheme();
    themeButton = h('button.btn.ghost.sm.nav-theme', {
      id: 'btn-theme',
      type: 'button',
      'data-theme-toggle': '',
      'aria-pressed': String(t0 === 'light'),
      // Labelled `theme: <state>`, which is what the Studio's button says.
      // UI_ARCHITECTURE.md §5 asks a toggle to state the state it produces,
      // and asks two controls that do the same thing to carry the same label.
      // Those pull opposite ways here and the second wins: one relabelled
      // toggle across two pages is two controls again.
      title: 'Toggle light / dark theme. The choice is shared with the Studio.',
      onclick: () => {
        if (typeof o.onThemeToggle === 'function') {
          o.onThemeToggle();
          // Re-assert rather than flip: the host has already changed
          // `data-theme`, and this relabels every toggle on the page and wakes
          // the subscribers. Writing the same value twice is free; flipping
          // twice would cancel the user's click.
          applyTheme(currentTheme());
        } else {
          toggleTheme();
        }
      },
    }, 'theme: ' + t0);
    actionsEl.appendChild(themeButton);
  }

  const groups = [brandEl, linksEl, projEl, h('div.topbar-spacer'), statusEl, actionsEl];
  const el = nested
    ? h('div.nav-bar.nav-bar-nested', groups)
    : h('header.topbar.nav-bar', { role: 'banner' }, groups);

  if (nested) {
    // `.topbar > *` in css/layout.css is what makes each group a flex row, and
    // it is a CHILD combinator: nested one level down, none of these groups
    // would get it. Set here rather than asking home.css for a rule, because
    // this module has to render correctly on any page that mounts it.
    el.style.display = 'flex';
    el.style.alignItems = 'center';
    el.style.flexWrap = 'wrap';
    el.style.width = '100%';
    el.style.minWidth = '0';
    for (const g of groups) { g.style.display = 'flex'; g.style.alignItems = 'center'; }
  } else {
    // `.topbar` is `grid-row: 1` of the Studio's shell grid; on a page that is
    // not that grid it collapses to content height, so the one height the
    // design system names for this bar is stated here rather than left to
    // whichever page mounts it. --h-bar is the token, never a literal.
    el.style.minHeight = 'var(--h-bar)';
  }

  // Nested, the host IS the bar and this fills a slot inside it, so append.
  // Otherwise this bar IS the top bar and must be the mount's FIRST child:
  // css/home.css styles `.home-shell > .topbar` as `position: sticky; top: 0`,
  // and a sticky bar appended after `<main>` sticks to the bottom of the
  // scroll rather than the top of the page.
  const mount = host || document.body;
  if (nested || !mount.firstChild) mount.appendChild(el);
  else mount.insertBefore(el, mount.firstChild);

  /* The account entry point, mounted AFTER the bar is in the document:
     js/ui/account.js finds its slot with `$('.topbar-group.last')` and silently
     does nothing if that element is not there yet. */
  let accountPanel = null;
  if (o.account !== false && !document.getElementById('btn-account')) {
    accountPanel = createAccountPanel();
  }

  /* --- state -------------------------------------------------------- */

  function setProject(p) {
    project = p && p.id ? p : null;
    setHidden(projEl, !project);
    projIdEl.textContent = project ? project.id : '';
    projNameEl.textContent = project ? (project.name || '') : '';
    projEl.title = project ? `${project.id} · ${project.name || ''}` : '';
    studioLink.href = project
      ? `${studioHref}?project=${encodeURIComponent(project.id)}`
      : studioHref;
    studioLink.title = project
      ? `Studio, scoped to ${project.id}. Its runs open attached to this project.`
      : 'Studio: reconstruct a specimen, run a therapy on it, read the mechanism';
  }

  /**
   * @param {object} [a]  api.js's `api` object: `{stub, stubReason, reachable}`.
   *
   * STUB is labelled, loudly, in the same words the Studio uses (contract §6).
   * This bar STATES it; it does not set `data-mode` or latch the provenance
   * chips — those are session-wide effects owned by the page's own boot, and
   * two modules writing them is exactly the §2.1 failure this file avoids.
   */
  function setBackend(a) {
    const s = a || api.api;
    if (s && s.stub) {
      ledEl.className = 'led nav-led warn';
      backendEl.textContent = 'STUB DATA · ' + (s.stubReason || 'backend unreachable');
      backendEl.style.color = 'var(--synthetic)';
      backendEl.title =
        'The FastAPI backend did not answer (' + (s.stubReason || 'backend unreachable')
        + '). Every response is being synthesised locally by js/mock.js from a '
        + 'seeded PRNG. There is no tissue, no solve and no measurement behind '
        + 'any of it. Projects are not synthesised at all, so creating, '
        + 'renaming, duplicating and deleting a project are refused until the '
        + 'backend answers.';
    } else if (s && s.reachable) {
      ledEl.className = 'led nav-led ok';
      backendEl.textContent = 'live · ' + location.host;
      backendEl.style.color = 'var(--fg-1)';
      backendEl.title = 'Connected to the FastAPI backend.';
    } else {
      ledEl.className = 'led nav-led';
      backendEl.textContent = 'probing…';
      backendEl.style.color = 'var(--fg-1)';
      backendEl.title = 'Backend connection state';
    }
  }

  let releaseModeHook = null;
  if (o.watchBackend !== false && !api.api.onModeChange) {
    api.api.onModeChange = setBackend;
    releaseModeHook = () => {
      if (api.api.onModeChange === setBackend) api.api.onModeChange = null;
    };
  }

  setProject(project);
  setBackend(api.api);

  const controller = {
    el,
    setProject,
    setBackend,
    themeButton,
    accountPanel,
    destroy() {
      if (releaseModeHook) releaseModeHook();
      delete el.__navController;
      el.remove();
    },
  };
  el.__navController = controller;
  return controller;
}

/**
 * The Studio path: its header already exists and js/main.js owns it. Add the
 * link back to `/` (contract §8 step 2) and the project context beside it, and
 * touch nothing else.
 */
function adoptExistingBar(bar, { homeHref, project }) {
  let group = bar.querySelector('.nav-links');
  let link = group && group.querySelector('.nav-link-home');
  if (!link) {
    link = h('a.btn.ghost.sm.nav-link.nav-link-home', {
      href: homeHref,
      title: 'Back to the workspace: every project on this server',
    }, 'workspace');
    group = h('nav.topbar-group.keep.nav-links', { 'aria-label': 'Pages' }, link);
    const brand = bar.querySelector('.brand');
    // insertBefore(node, null) appends, which is the right fallback when the
    // brand is the last child — and `bar` is the parent in both branches, so
    // the reference node is always a child of it.
    bar.insertBefore(group, brand ? brand.nextSibling : bar.firstChild);
  }

  let projEl = bar.querySelector('.nav-project');
  let projIdEl, projNameEl;
  if (!projEl) {
    projIdEl = h('span.mono.strong.nav-project-id');
    projNameEl = h('span.dim.nav-project-name');
    projEl = h('div.topbar-group.nav-project',
      h('span.lbl', 'Project'), projIdEl, projNameEl);
    bar.insertBefore(projEl, group.nextSibling);
  } else {
    projIdEl = projEl.querySelector('.nav-project-id');
    projNameEl = projEl.querySelector('.nav-project-name');
  }

  function setProject(p) {
    const cur = p && p.id ? p : null;
    setHidden(projEl, !cur);
    if (projIdEl) projIdEl.textContent = cur ? cur.id : '';
    if (projNameEl) projNameEl.textContent = cur ? (cur.name || '') : '';
    projEl.title = cur ? `${cur.id} · ${cur.name || ''}` : '';
  }
  setProject(project);

  return {
    el: bar,
    setProject,
    // The Studio paints its own LED from js/main.js `paintBackendState`.
    // Overwriting it from here would give one strip two owners.
    setBackend: () => {},
    themeButton: document.getElementById('btn-theme'),
    accountPanel: null,
    destroy() {
      projEl.remove();
      if (group) group.remove();
    },
  };
}
