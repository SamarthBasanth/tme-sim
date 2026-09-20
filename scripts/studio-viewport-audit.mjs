#!/usr/bin/env node
/**
 * Is every control in the Studio shell REACHABLE, at every width we serve?
 *
 * WHY THIS EXISTS AS A SCRIPT AND NOT AS A PARAGRAPH IN A CSS COMMENT.
 * `app/web/css/layout.css` is ~1,100 lines of shell rules and nearly every one
 * of them is justified by a number someone measured by hand once: "the grid's
 * floor was 912px", "the bar's content is 924px wide", "MEASURED at 720x450".
 * Those numbers were right when they were taken and there was no way to take
 * them again, so the same class of defect has now shipped three times in two
 * days:
 *
 *   1. a 912px three-rail grid floor inside `html, body { overflow: hidden }`,
 *      which DESTROYED up to 552px of the READ rail rather than scrolling it.
 *      Invisible to `scrollWidth`, because a clip is not an overflow.
 *   2. a first-run tour card with no height bound, whose `Next` / `Skip tour`
 *      controls left the bottom of a 450px-tall window.
 *   3. `.topbar-group.last` — six buttons, 335.1px, `position: sticky;
 *      right: 0` with an OPAQUE fill and a higher stacking order — painting
 *      over the `← workspace` link, which is the only route back to `/`.
 *      Below 336px it covered the whole bar at every scroll position, so a
 *      320px phone that opened a project could not leave it.
 *
 * Every one of those is the same bug: something is on screen in the DOM, has a
 * non-zero rectangle, answers `getBoundingClientRect()` honestly, and CANNOT BE
 * CLICKED. None of them move `scrollWidth`. None of them throw. Reading the CSS
 * does not find them; only asking the browser "if a finger landed here, what
 * would it hit?" finds them.
 *
 * WHAT IT MEASURES
 *
 *   REACHABILITY  For every focusable control whose box is on screen, sample a
 *                 grid of points inside it and ask `document.elementFromPoint`.
 *                 That is literally what a click does, so an element that fails
 *                 here is an element a user cannot operate. Controls named in
 *                 `SOLE_ROUTES` must be 100% reachable — they are the only way
 *                 to do the thing they do, and half a button is not a route.
 *                 Everything else must be reachable at its own centre.
 *   CLIPPING      Every control, against each `overflow: hidden` ancestor. This
 *                 is defect (1) above, and it is the one a scroll-width check
 *                 cannot see: the pixels are not pushed out of the box, they
 *                 are erased inside it.
 *   OVERFLOW      `documentElement` and `.shell` must never scroll
 *                 horizontally. The shell's own contract is that it IS the
 *                 viewport.
 *
 * AND IT MEASURES THEM WHILE SCROLLING. Below 900px the shell is one column and
 * `.body` is the scroller, with a sticky region head, a sticky rail navigator
 * and a sticky commit bar riding on it. A sticky element that fits at
 * `scrollTop: 0` can still clamp over the form it commits further down — that
 * is exactly how the commit bar came to paint over its own fields at 720x450.
 * So the body scroller is sampled across its whole range and the WORST case is
 * reported.
 *
 * AND IN BOTH FIRST-RUN STATES. The onboarding tour is `position: fixed;
 * z-index: 1001` and is shown once per browser, so a first-time visitor and a
 * returning one are two different stacking situations over the same shell.
 * Measured before the exit-lane fix: with the tour up, the `← workspace` link
 * was 0% reachable at 415px and below in states where it was still 8%
 * reachable without it. Auditing only the returning visitor would have called
 * the worse of the two states clean.
 *
 * USAGE
 *   node scripts/studio-viewport-audit.mjs
 *   VPS=320x568,900x800 node scripts/studio-viewport-audit.mjs
 *   SAMPLES=12 node scripts/studio-viewport-audit.mjs   # body-scroll samples
 *   PORT=3187 node scripts/studio-viewport-audit.mjs
 *
 * Exits 1 if any control is covered, clipped or off-screen-unreachable, or if
 * the shell scrolls horizontally. Exits 2 if it CANNOT RUN — not being able to
 * run the check is not the same as the check passing.
 *
 * Playwright is deliberately not a dependency of this repo (there is no
 * package.json at the root and adding one for a ~100MB browser download would
 * be a poor trade). It is resolved from wherever it already exists, exactly as
 * landing-3d/scripts/viewport-audit.mjs does.
 */

import { createServer } from "node:http";
import { readFile } from "node:fs/promises";
import { existsSync } from "node:fs";
import { createRequire } from "node:module";
import { dirname, extname, join, normalize, resolve } from "node:path";
import { fileURLToPath, pathToFileURL } from "node:url";

const REPO = resolve(dirname(fileURLToPath(import.meta.url)), "..");
const WEB = join(REPO, "app", "web");
const require = createRequire(import.meta.url);

const PORT = Number(process.env.PORT ?? 3187);
const SAMPLES = Number(process.env.SAMPLES ?? 8);

/**
 * The standing set. Each is a real class of device or a boundary this shell has
 * been wrong at before. The heights are as load-bearing as the widths — three
 * of the shell's four region sizes are height constants, and the commit-bar
 * pin and the tour card both failed on the height axis while every width was
 * fine.
 *
 *   320x568   iPhone SE 1. The narrowest viewport served, and the one where a
 *             335.1px sticky launcher did not fit the screen it was pinned to:
 *             the exit was covered at EVERY scroll position.
 *   360x740   common Android portrait
 *   375x667   iPhone SE 2/3, 8
 *   390x844   iPhone 12-15
 *   412x915   Pixel
 *   464x800   the width at which the exit link was ENTIRELY under the launcher
 *   505x800   the width at which the exit link's own MIDPOINT was covered —
 *             i.e. where clicking the middle of it stopped working, which is
 *             the width the defect was first reported at
 *   538x800   FIRST CONTACT. The launcher reaches the exit at
 *             203.6 + 335.1 = 538.7px, so 539 is the last clean width and 538
 *             is the first broken one — by 0.7px, which is why sole routes are
 *             sampled at 1px and not on a grid.
 *   720x450   a 1440x900 laptop at 200% zoom. The viewport the dock's height
 *             budget, the commit bar's pin and the tour card all failed at,
 *             and the reason this set is not width-only.
 *   768x1024  iPad portrait
 *   899x800   one pixel below the stacked breakpoint
 *   900x800   the breakpoint itself: the narrowest three-rail workstation,
 *             where the elastic tracks are under maximum pressure
 *   1280x720  the 16:9 laptop
 *   1440x900  the viewport the shell was authored at
 */
const DEFAULT_VIEWPORTS = [
  "320x568", "360x740", "375x667", "390x844", "412x915",
  "464x800", "505x800", "538x800", "720x450", "768x1024",
  "899x800", "900x800", "1280x720", "1440x900",
];

const VIEWPORTS = (process.env.VPS ?? DEFAULT_VIEWPORTS.join(",")).split(",").map((label) => {
  const [width, height] = label.trim().split("x").map(Number);
  if (!Number.isFinite(width) || !Number.isFinite(height)) {
    throw new Error(`bad viewport "${label}" — expected WIDTHxHEIGHT`);
  }
  return { label: label.trim(), width, height };
});

// ---------------------------------------------------------------------------
// A browser, without taking a dependency on one
// ---------------------------------------------------------------------------

async function loadChromium() {
  const candidates = [
    process.env.PLAYWRIGHT_MODULE,
    join(REPO, "landing-3d", "node_modules", "playwright", "index.mjs"),
    ...(process.env.HOME
      ? [join(process.env.HOME, ".npm/_npx/e41f203b7505f1fb/node_modules/playwright/index.mjs")]
      : []),
  ].filter(Boolean);
  for (const c of candidates) {
    if (!existsSync(c)) continue;
    try { return (await import(pathToFileURL(c).href)).chromium; } catch { /* next */ }
  }
  try { return require("playwright").chromium; } catch { /* fall through */ }
  console.error(
    "studio viewport audit CANNOT RUN: no playwright found.\n" +
    "  Tried:\n" + candidates.map((c) => `    ${c}`).join("\n") +
    "\n  Set PLAYWRIGHT_MODULE=/path/to/playwright/index.mjs, or `npx playwright install chromium`.\n" +
    "  Not being able to run the check is not the same as the check passing.",
  );
  process.exit(2);
}

// ---------------------------------------------------------------------------
// Serving app/web/ with NO backend
// ---------------------------------------------------------------------------
//
// Deliberately not `uvicorn app.server.api:app`. This audit is about the shell
// — boxes, stacking and hit testing — and a backend adds a python environment,
// a data cache and whichever projects happen to exist on the machine to the
// list of things that can make it red for reasons that are not the subject.
//
// Every `/api/*` request is answered 503, which is a state the front end is
// built for: js/api.js flips `api.stub = true` and the Studio renders from
// js/mock.js with a STUB DATA badge. The shell is the same shell.
//
// The one place that matters is the project context, because the exit link is
// mounted by BOTH of js/main.js's branches — `showProjectContext` on success
// and `showProjectUnavailable` on failure — and in both it is the FIRST child
// of `#hdr-project-group`, so its box is in the same place either way. The
// audit asserts the link exists at all (`NO EXIT LINK` is a failure, not a
// skip), so a refactor that drops it from one branch cannot pass quietly.

const MEDIA = {
  ".html": "text/html; charset=utf-8", ".js": "text/javascript; charset=utf-8",
  ".css": "text/css; charset=utf-8", ".woff2": "font/woff2",
  ".json": "application/json", ".svg": "image/svg+xml", ".txt": "text/plain",
};

async function startServer() {
  if (!existsSync(join(WEB, "index.html"))) {
    console.error(`studio viewport audit CANNOT RUN: ${join(WEB, "index.html")} does not exist.`);
    process.exit(2);
  }
  const server = createServer(async (req, res) => {
    let path = decodeURIComponent(new URL(req.url, "http://x").pathname);
    if (path === "/studio") path = "/index.html";          // app/server/api.py::_mount_web
    if (path === "/") path = "/home.html";
    if (path.startsWith("/api/")) { res.writeHead(503, MEDIA[".json"]); res.end("{}"); return; }
    const file = join(WEB, normalize(path).replace(/^(\.\.[/\\])+/, ""));
    if (!file.startsWith(WEB)) { res.writeHead(403); res.end("no"); return; }
    try {
      const body = await readFile(file);
      res.writeHead(200, { "content-type": MEDIA[extname(file)] ?? "application/octet-stream" });
      res.end(body);
    } catch { res.writeHead(404); res.end("not found"); }
  });
  await new Promise((r) => server.listen(PORT, "127.0.0.1", r));
  return { base: `http://127.0.0.1:${PORT}`, stop: () => server.close() };
}

// ---------------------------------------------------------------------------
// The in-page measurement
// ---------------------------------------------------------------------------

/**
 * The controls that are the ONLY route to the thing they do. These must be
 * whole: a button that is 60% clickable is a button that fails for whoever
 * aims at the covered 40%, and a route out of a screen that works if you aim
 * left is not a route out.
 *
 * Keyed by CSS selector so the list reads as the claim it is. `href="/"` rather
 * than a class because js/main.js's `workspaceLink()` gives the element no
 * class of its own.
 */
const SOLE_ROUTES = {
  '#hdr-project-group > a[href="/"]': "the route back to the workspace home",
  "#btn-protocol": "the protocol workspace (also P, which needs a keyboard)",
  "#btn-account": "the account panel / API key",
  "#btn-projects": "the projects panel",
  "#btn-help": "the keyboard map",
  "#btn-theme": "the theme toggle (also T, which needs a keyboard)",
  "#btn-tour": "replaying the orientation tour",
};

/**
 * THE PERSISTENT CHROME. Everything here is furniture that does not scroll
 * away: it is on screen the whole time you use the Studio, so if something is
 * painted over it, it is painted over it permanently. That is the shape of all
 * three shipped defects, and it is why the general "is this clickable" rule is
 * scoped to these rather than sprayed over all ~450 focusable elements — a
 * control halfway down a scrolling rail that is briefly under a sticky header
 * is a different (and much weaker) claim, and mixing the two produces a check
 * nobody believes.
 *
 * `.rail-left > :last-child` is the commit bar: `position: sticky` and, per
 * UI_ARCHITECTURE.md 2.2, "the one control whose reachability must not depend
 * on scroll position". It has already shipped clamped over the form it submits.
 */
const CHROME = [".topbar", ".stage-toolbar", ".statusbar", ".rail-foot", ".rail-left > :last-child"];

/**
 * Overlays that are ON TOP BY DESIGN and dismissible. Coverage BY one of these
 * is not a defect — that is what an overlay is. Two consequences, and the
 * second is the one that matters:
 *   * they are subtracted from the hit test, so the shell underneath is judged
 *     on its own geometry rather than on whatever notification happened to be
 *     up;
 *   * their own controls are held to the STRICTEST rule in this file — on
 *     screen and clickable, every time. Nothing legitimately covers an
 *     overlay, and an overlay whose dismiss button is off the bottom of the
 *     window (which is exactly how the first-run tour shipped) traps the
 *     reader behind it.
 */
const OVERLAYS = [".tour-popover", ".toast", "dialog[open]", ".gate-overlay"];

const MEASURE = (cfg) => {
  const de = document.documentElement;
  const VW = de.clientWidth;
  const VH = de.clientHeight;
  const out = { controls: [], overflow: [], hOverflow: de.scrollWidth - de.clientWidth };
  const shell = document.querySelector(".shell");
  out.shellHOverflow = shell ? shell.scrollWidth - shell.clientWidth : 0;

  const name = (el) => {
    if (!el) return "(nothing)";
    if (el.id) return `#${el.id}`;
    const cls = typeof el.className === "string" && el.className
      ? "." + el.className.trim().split(/\s+/).slice(0, 2).join(".") : "";
    const txt = (el.textContent || "").trim().replace(/\s+/g, " ").slice(0, 18);
    return `${el.tagName.toLowerCase()}${cls}${txt ? ` "${txt}"` : ""}`;
  };

  const inOverlay = (el) => !!el && OVERLAYS_SEL.some((s) => el.closest(s));
  const OVERLAYS_SEL = cfg.overlays;

  /** Does anything between `el` and `stop` scroll on this axis? */
  const scrollsBetween = (el, stop, axis) => {
    for (let p = el.parentElement; p && p !== stop; p = p.parentElement) {
      if (axis === "x" && p.scrollWidth - p.clientWidth > 1) return true;
      if (axis === "y" && p.scrollHeight - p.clientHeight > 1) return true;
    }
    return false;
  };

  /**
   * What is actually LEFT of an element's box after every ancestor that clips
   * has had its say, and — the part that matters — WHY the rest went.
   *
   * These are three different states that a single rectangle cannot tell
   * apart, and conflating them is how an audit ends up either blind or so
   * noisy nobody reads it:
   *
   *   SCROLLED   an ancestor clipped it AND that ancestor scrolls on that
   *              axis. The pixels are one gesture away. `.stage-toolbar` is
   *              833px of controls in a 768px stage and is MEANT to be —
   *              it scrolls, and js/util/dom.js paints an edge cue on it.
   *              Not a defect, and an audit that calls it one will be
   *              switched off within a day.
   *   CLIPPED    an ancestor clipped it and that ancestor does NOT scroll.
   *              The pixels are gone. Nothing in the DOM says so, and
   *              `scrollWidth` does not move. This is the 912px grid floor
   *              that erased 552px of the READ rail.
   *   COVERED    the box is on screen and unclipped, and something else is
   *              painted over it. This is the sticky launcher over the exit.
   */
  const visibleRect = (el) => {
    const r = el.getBoundingClientRect();
    let box = { l: r.left, r: r.right, t: r.top, b: r.bottom };
    let clip = null, scrolled = false;
    const eat = (px, by, isScroll) => {
      if (px <= 1) return;
      if (isScroll) scrolled = true;
      else if (!clip || px > clip.px) clip = { px: +px.toFixed(1), by };
    };
    for (let p = el.parentElement; p; p = p.parentElement) {
      const cs = getComputedStyle(p);
      const pr = p.getBoundingClientRect();
      const clipsX = cs.overflowX !== "visible";
      const clipsY = cs.overflowY !== "visible";
      if (!clipsX && !clipsY) continue;
      const sx = p.scrollWidth - p.clientWidth > 1 || scrollsBetween(el, p, "x");
      const sy = p.scrollHeight - p.clientHeight > 1 || scrollsBetween(el, p, "y");
      if (clipsX) {
        eat(pr.left - box.l, name(p), sx); eat(box.r - pr.right, name(p), sx);
        box = { ...box, l: Math.max(box.l, pr.left), r: Math.min(box.r, pr.right) };
      }
      if (clipsY) {
        eat(pr.top - box.t, name(p), sy); eat(box.b - pr.bottom, name(p), sy);
        box = { ...box, t: Math.max(box.t, pr.top), b: Math.min(box.b, pr.bottom) };
      }
    }
    // The viewport is the last clipper, and it never scrolls: `.shell` is
    // `position: fixed; inset: 0` and `html, body { overflow: hidden }`.
    box = { l: Math.max(box.l, 0), r: Math.min(box.r, VW), t: Math.max(box.t, 0), b: Math.min(box.b, VH) };
    const whole = box.l <= r.left + 1 && box.r >= r.right - 1 &&
                  box.t <= r.top + 1 && box.b >= r.bottom - 1;
    return { rect: r, box, clip, scrolled, whole, empty: box.r - box.l < 1 || box.b - box.t < 1 };
  };

  /**
   * What a finger landing at (x, y) would actually operate.
   * `elementsFromPoint` gives the whole paint stack, so this answers the real
   * question — "is `el` reached before anything that would swallow the event"
   * — instead of the approximate one a single `elementFromPoint` answers.
   * `el`'s own ancestors and descendants are not blockers (a <td> inside the
   * <tr> being tested IS the <tr>), and neither is an overlay unless `el` is
   * itself inside one.
   */
  const blockerAt = (el, x, y, ignoreOverlays) => {
    if (x < 0 || y < 0 || x >= VW || y >= VH) return { off: true };
    for (const n of document.elementsFromPoint(x, y)) {
      if (n === el || el.contains(n) || n.contains(el)) return null;   // reached it
      if (ignoreOverlays && inOverlay(n)) continue;
      return { by: name(n) };
    }
    return { by: "(nothing painted here)" };
  };

  const measure = (el, kind, sole) => {
    const v = visibleRect(el);
    // Overlays are never blockers — not for the shell underneath them, and not
    // for each other. A dismissible card on top of a dismissible toast is two
    // things to dismiss, not a trap; the shell underneath is judged on its own
    // geometry either way. What an overlay's controls ARE held to is being on
    // screen and whole, which is checked below and is the rule the first-run
    // tour broke when `Next` and `Skip tour` left the bottom of the window.
    const ignoreOverlays = true;
    let hits = 0, total = 0, coverer = null, centreOk = true, centreBy = null;
    if (!v.empty) {
      const x0 = v.box.l + 1, x1 = v.box.r - 1, y0 = v.box.t + 1, y1 = v.box.b - 1;
      // A SOLE ROUTE IS SAMPLED AT 1px, EVERYTHING ELSE ON A COARSE GRID, and
      // the difference is the check's RESOLUTION — which is a number this file
      // has to be honest about rather than leave to be discovered.
      // A 9-point grid across a 76px button steps 8.4px, so it cannot see a
      // sliver thinner than that: the launcher first touches the exit at 538px
      // and the overlap there is 0.7px. A coarse grid calls 538 clean, which
      // is true of what it sampled and false of the shell. At 1px the first
      // covered pixel is the first covered pixel.
      // It costs ~230 hit tests per sole route per sample; there are seven of
      // them, which is affordable, and spraying the same density over ~450
      // focusable elements would not be.
      const nx = sole
        ? Math.max(2, Math.ceil(x1 - x0) + 1)
        : Math.max(2, Math.min(9, Math.floor((x1 - x0) / 4) + 2));
      const ny = Math.max(2, Math.min(4, Math.floor((y1 - y0) / 4) + 2));
      for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) {
        const x = x1 <= x0 ? x0 : x0 + ((x1 - x0) * i) / (nx - 1);
        const y = y1 <= y0 ? y0 : y0 + ((y1 - y0) * j) / (ny - 1);
        total++;
        const bl = blockerAt(el, x, y, ignoreOverlays);
        if (!bl) hits++; else if (!coverer && bl.by) coverer = bl.by;
      }
      const cb = blockerAt(el, (v.box.l + v.box.r) / 2, (v.box.t + v.box.b) / 2, ignoreOverlays);
      centreOk = !cb; centreBy = cb && cb.by ? cb.by : null;
    }
    return {
      n: name(el), kind, sole,
      onScreen: !v.empty, whole: v.whole, scrolled: v.scrolled, clip: v.clip,
      frac: total ? +(hits / total).toFixed(3) : 0,
      centreOk, centreBy, coverer,
      rect: { l: +v.rect.left.toFixed(1), r: +v.rect.right.toFixed(1),
              t: +v.rect.top.toFixed(1), b: +v.rect.bottom.toFixed(1) },
    };
  };

  const FOCUSABLE = 'a[href], button, input, select, textarea, summary, [tabindex]:not([tabindex="-1"])';
  const visible = (el) => {
    const cs = getComputedStyle(el);
    if (cs.display === "none" || cs.visibility === "hidden" || Number(cs.opacity) === 0) return false;
    if (el.closest("[hidden], .hidden")) return false;
    if (el.closest("dialog:not([open])")) return false;
    // The skip link is parked outside the viewport until it takes focus —
    // the standard pattern, and base.css's `.skip:focus` brings it back. It
    // is off screen on purpose, which is the one thing this audit is not
    // looking for.
    if (el.classList.contains("skip")) return false;
    const r = el.getBoundingClientRect();
    return r.width >= 1 && r.height >= 1;
  };

  const done = new Set();
  // (a) the sole routes, by name, so a MISSING one is a failure and not a skip
  for (const sel of cfg.sole) {
    const el = document.querySelector(sel);
    if (!el || !visible(el)) { out.controls.push({ n: sel, sole: sel, kind: "sole", missing: true }); continue; }
    done.add(el);
    out.controls.push(measure(el, "sole", sel));
  }
  // (b) the persistent chrome, and (c) whatever is inside an open overlay
  for (const [kind, roots] of [["chrome", cfg.chrome], ["overlay", cfg.overlays]]) {
    for (const root of document.querySelectorAll(roots.join(","))) {
      for (const el of root.querySelectorAll(FOCUSABLE)) {
        if (done.has(el) || !visible(el)) continue;
        done.add(el);
        out.controls.push(measure(el, kind, null));
      }
    }
  }
  // (d) clipping, over EVERY focusable element in the shell — this is the axis
  //     a scroll-width check is blind to, so it is not scoped to the chrome.
  for (const el of document.querySelectorAll(FOCUSABLE)) {
    if (done.has(el) || !visible(el)) continue;
    const v = visibleRect(el);
    if (!v.clip) continue;
    done.add(el);
    out.controls.push({
      n: name(el), kind: "body", sole: null, onScreen: !v.empty, whole: v.whole,
      scrolled: v.scrolled, clip: v.clip, frac: 1, centreOk: true, centreBy: null, coverer: null,
      rect: { l: +v.rect.left.toFixed(1), r: +v.rect.right.toFixed(1),
              t: +v.rect.top.toFixed(1), b: +v.rect.bottom.toFixed(1) },
    });
  }

  // (e) anything sticking out of the viewport with nothing to scroll it back
  for (const el of document.querySelectorAll(".shell *")) {
    const r = el.getBoundingClientRect();
    if (r.width === 0 && r.height === 0) continue;
    if (r.right <= VW + 0.5 && r.left >= -0.5) continue;
    if (getComputedStyle(el).position === "fixed") continue;
    if (scrollsBetween(el, null, "x")) continue;      // reachable with a scrollbar
    out.overflow.push({ n: name(el), l: +r.left.toFixed(1), r: +r.right.toFixed(1) });
  }
  out.overflow = out.overflow.slice(0, 6);
  return out;
};

// ---------------------------------------------------------------------------

const worse = (a, b) => {
  // Keep the worst observation of a control across every scroll position.
  if (!a) return b;
  if (b.missing) return b;
  const score = (r) => (r.missing ? -2 : 0) + (r.onScreen ? r.frac : -1) +
    (r.centreOk ? 0.001 : 0) - (r.clip ? 1 : 0);
  return score(b) < score(a) ? b : a;
};

async function auditViewport(browser, base, vp, firstRun) {
  const ctx = await browser.newContext({
    viewport: { width: vp.width, height: vp.height },
    deviceScaleFactor: 1,
  });
  // A returning visitor has already been offered the tour; a first-time one has
  // a `position: fixed; z-index: 1001` card over the shell. Both are real.
  if (!firstRun) {
    await ctx.addInitScript(() => {
      try { localStorage.setItem("tmesim_welcome_seen_v1", "1"); } catch { /* private mode */ }
    });
  }
  const page = await ctx.newPage();
  try {
    await page.goto(`${base}/studio?project=p_audit`, { waitUntil: "networkidle" });
    // The self-hosted faces and the rails land a frame or two late; measuring
    // before they do measures the fallback metrics, which is a different shell.
    await page.waitForTimeout(1800);

    const cfg = { sole: Object.keys(SOLE_ROUTES), chrome: CHROME, overlays: OVERLAYS };
    const worstBy = new Map();
    const merge = (m) => {
      for (const c of m.controls) {
        const k = `${c.kind}|${c.sole ?? c.n}`;
        worstBy.set(k, worse(worstBy.get(k), c));
      }
    };

    // Sample the one scroller the shell has in stacked mode. In workstation
    // mode `.body` does not scroll and the loop collapses to a single sample.
    const range = await page.evaluate(() => {
      const b = document.getElementById("body");
      return b ? Math.max(0, b.scrollHeight - b.clientHeight) : 0;
    });
    let overflow = [], hOverflow = 0, shellH = 0;
    for (let i = 0; i <= (range > 0 ? SAMPLES : 0); i++) {
      const y = range > 0 ? Math.round((range * i) / SAMPLES) : 0;
      await page.evaluate((yy) => {
        const b = document.getElementById("body");
        if (b) b.scrollTop = yy;
      }, y);
      await page.waitForTimeout(30);
      const m = await page.evaluate(MEASURE, cfg);
      merge(m);
      if (m.overflow.length > overflow.length) overflow = m.overflow;
      hOverflow = Math.max(hOverflow, m.hOverflow);
      shellH = Math.max(shellH, m.shellHOverflow ?? 0);
    }
    return { vp: vp.label, firstRun, controls: [...worstBy.values()], overflow, hOverflow, shellH, range };
  } finally {
    await ctx.close();
  }
}

async function main() {
  const chromium = await loadChromium();
  const { base, stop } = await startServer();
  const browser = await chromium.launch();
  const failures = [];

  try {
    for (const vp of VIEWPORTS) {
      for (const firstRun of [false, true]) {
        const r = await auditViewport(browser, base, vp, firstRun);
        const tag = `${r.vp} ${firstRun ? "first-run" : "returning"}`.padEnd(24);

        const bad = [];
        for (const c of r.controls) {
          if (c.missing) { bad.push({ c, why: `MISSING from the DOM — ${SOLE_ROUTES[c.sole]}` }); continue; }
          // CLIPPED is a failure for anything, anywhere. Pixels erased by an
          // `overflow: hidden` box that does not scroll are pixels no gesture
          // brings back, and nothing else in this repo can see them.
          if (c.clip) {
            bad.push({ c, why: `CLIPPED ${c.clip.px}px by ${c.clip.by} — overflow hidden with nothing to scroll it back, so these pixels are ERASED` });
          }
          if (c.sole) {
            // The strictest rule in the file, and the one the exit failed.
            // A route that is the only route must be WHOLE and CLICKABLE: a
            // 60%-reachable button fails for whoever aims at the other 40%,
            // and "scroll the bar sideways and it comes out from under the
            // launcher" is not a route anyone finds.
            if (!c.onScreen) bad.push({ c, why: `NOT ON SCREEN — ${SOLE_ROUTES[c.sole]}` });
            else if (c.frac < 1) bad.push({ c, why: `COVERED ${Math.round((1 - c.frac) * 100)}% by ${c.coverer} — ${SOLE_ROUTES[c.sole]}` });
            // AND IT MUST BE WHOLE, with no "it is only one scroll away"
            // escape. That escape is legitimate for the stage toolbar; it is
            // not legitimate here, for the reason layout.css already gives
            // about the launcher: "a control you have to discover a horizontal
            // scroll to reach is a control most operators never find". It is
            // also the rule that catches the OTHER half of the shipped defect
            // — `.topbar-group.last` is 335.1px pinned to the right edge of a
            // 320px window, so its own `theme` button hung 15.1px off the
            // screen, and no amount of scrolling brought it back because the
            // group is pinned and the scroll takes the pin with it.
            else if (!c.whole) bad.push({ c, why: `only partly on screen (clipped by ${c.clip ? c.clip.by : "the viewport"}) — ${SOLE_ROUTES[c.sole]}` });
          } else if (c.kind === "overlay") {
            // An overlay is on top by definition, so the only way one of its
            // controls fails is by leaving the window — which is exactly how
            // the first-run tour shipped with `Next` and `Skip tour` below the
            // bottom of a 450px-tall viewport.
            if (!c.onScreen) bad.push({ c, why: "an OVERLAY control that is not on screen — there is nothing to scroll it back" });
            else if (!c.whole) bad.push({ c, why: `an OVERLAY control cut off by ${c.clip ? c.clip.by : "the viewport"}` });
            else if (!c.centreOk) bad.push({ c, why: `an OVERLAY control covered by ${c.centreBy}, which is not itself an overlay` });
          } else if (c.kind === "chrome") {
            // Persistent furniture. It may legitimately be one scroll away
            // (the stage toolbar is 833px of controls in a 768px stage, with
            // an edge cue), but it may not be painted over where it sits.
            if (c.onScreen && !c.centreOk) bad.push({ c, why: `centre COVERED by ${c.centreBy}` });
          }
        }
        const ok = bad.length === 0 && r.hOverflow <= 0 && r.shellH <= 0 && r.overflow.length === 0;

        console.log(
          `${tag} ${String(r.controls.filter((c) => c.sole).length).padStart(2)} sole routes` +
          `  hOverflow ${r.hOverflow}px  bodyScroll ${r.range}px  ` +
          (ok ? "all controls reachable"
              : `${bad.length + r.overflow.length + (r.hOverflow > 0 ? 1 : 0) + (r.shellH > 0 ? 1 : 0)} problem(s)  <-- FAIL`),
        );
        for (const { c, why } of bad) {
          console.log(`      ${(c.n || c.sole).padEnd(34)} ${c.rect ? `x ${c.rect.l}..${c.rect.r} y ${c.rect.t}..${c.rect.b}` : ""}`);
          console.log(`      ${"".padEnd(34)} ${why}`);
        }
        for (const o of r.overflow) console.log(`      OVERFLOWS VIEWPORT  ${o.n}  ${o.l}..${o.r}`);
        if (r.hOverflow > 0) console.log(`      DOCUMENT SCROLLS HORIZONTALLY by ${r.hOverflow}px`);
        if (r.shellH > 0) console.log(`      .shell SCROLLS HORIZONTALLY by ${r.shellH}px`);
        if (!ok) failures.push(tag.trim());
      }
    }
  } finally {
    await browser.close();
    stop();
  }

  console.log("");
  if (failures.length) {
    console.log(`studio viewport audit FAILED at ${failures.length} viewport/state(s):`);
    for (const f of failures) console.log(`  ${f}`);
    process.exit(1);
  }
  console.log(`studio viewport audit passed — ${VIEWPORTS.length} viewports x 2 first-run states`);
}

main().catch((e) => { console.error(e); process.exit(1); });
