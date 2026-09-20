#!/usr/bin/env node
/**
 * Does every beat fit inside the pinned viewport, at every scroll position?
 *
 * WHY THIS EXISTS AS A SCRIPT AND NOT AS A PARAGRAPH IN A CSS COMMENT.
 * `app/globals.css` carries roughly five hundred lines of responsive rules and
 * every one of them is justified by a number someone measured by hand:
 * "beat 02 top 76.0 bottom 846.6, clipped 132.6px off the bottom". Those
 * numbers were right when they were taken and there was no way to take them
 * again. So the fitting problem was re-discovered once per pass, at a
 * different viewport each time, usually by reading rather than measuring — and
 * reading is how 320x568 shipped with five of six beats clipped and 320x881
 * with three, both invisible to every check the project had.
 *
 * WHAT IT MEASURES, AND WHY THE ANSWER IS NOT "DOES IT FIT AT REST".
 * `.scroll-pin` is `height: 100svh; overflow: hidden`, so a beat taller than
 * the viewport is not scrolled — it is CUT, and what gets cut off the bottom of
 * beats 02 and 05 is the last rows of their ledgers, i.e. the registered
 * limits, which are the one thing on this page that must never be the part that
 * goes missing. On top of that `sectionVisual` (lib/timeline.ts) translates a
 * beat by `(1 - opacity) * 26px` — down while it arrives, up while it leaves —
 * so a block that fits at peak opacity is still cut by up to 26px on the way in
 * and out. This therefore samples the WHOLE timeline and reports the worst case
 * over every position at which a beat is visible at all (opacity > 0.02, the
 * same threshold ScrollExperience uses to drop a beat out of the tab order).
 *
 * It also measures the two floating instrument panels against the copy beside
 * them, as a true rectangle intersection rather than a vertical overlap: at
 * desktop widths the copy and the panel are side by side and share a band of
 * height on purpose, so a vertical test reports a collision on every viewport
 * and is worse than no test.
 *
 * AND IT REPORTS WHICH LAYOUT RENDERED. Below `MIN_PINNED_VIEWPORT_HEIGHT`
 * (hooks/useShortViewport.ts) the page falls through to normal document flow,
 * where nothing can be clipped; there the check is horizontal overflow and
 * whether the hero copy fits the stage that clips it.
 *
 * USAGE
 *   node scripts/viewport-audit.mjs                 # the standing set
 *   VPS=320x568,360x740 node scripts/viewport-audit.mjs
 *   SAMPLES=240 node scripts/viewport-audit.mjs     # denser timeline sampling
 *   PORT=3199 node scripts/viewport-audit.mjs
 *
 * Exits non-zero if any beat is clipped, any panel collides with copy, or the
 * page scrolls horizontally.
 *
 * Playwright is deliberately NOT a dependency of this package, for the reason
 * `scripts/network-audit.mjs` gives at length: it is a ~100 MB browser download
 * in service of a check that runs by hand and in CI. It is resolved from
 * wherever it already exists.
 */

import { spawn } from "node:child_process";
import { createRequire } from "node:module";
import { dirname, join } from "node:path";
import { fileURLToPath, pathToFileURL } from "node:url";
import { existsSync } from "node:fs";

const APP_ROOT = join(dirname(fileURLToPath(import.meta.url)), "..");
const require = createRequire(import.meta.url);

const PORT = Number(process.env.PORT ?? 3142);
const SAMPLES = Number(process.env.SAMPLES ?? 120);

/**
 * The standing set. Each one is here because it is a real class of device or
 * because it is a boundary this page has been wrong at before:
 *
 *   320x568   iPhone SE 1 — the shortest viewport worth serving, and the one
 *             that had five of six beats clipped until 2026-09-20
 *   375x667   iPhone SE 2/3, 8 — below the gate, so it must report FLOW
 *   360x700   the gate's own boundary: the shortest PINNED viewport allowed
 *   360x740   a common Android portrait
 *   390x844   iPhone 12–15
 *   412x915   Pixel — taller than the 880px tier ceiling, so it exercises the
 *             desktop composition at a phone width
 *   320x881   narrow AND tall: one pixel past that ceiling, which is where
 *             three beats were clipped until the tier gained its second
 *             condition
 *   768x1024  iPad portrait
 *   1280x720  the 16:9 laptop the "short and wide" block exists for
 *   1440x900  the tightest desktop on the page (beat 02 lands ~13px clear)
 *   1512x950  the viewport the composition was authored at
 */
const DEFAULT_VIEWPORTS = [
  "320x568",
  "375x667",
  "360x700",
  "360x740",
  "390x844",
  "412x915",
  "320x881",
  "768x1024",
  "1280x720",
  "1440x900",
  "1512x950",
];

const VIEWPORTS = (process.env.VPS ?? DEFAULT_VIEWPORTS.join(",")).split(",").map((label) => {
  const [width, height] = label.trim().split("x").map(Number);
  if (!Number.isFinite(width) || !Number.isFinite(height)) {
    throw new Error(`bad viewport "${label}" — expected WIDTHxHEIGHT`);
  }
  return { label: label.trim(), width, height };
});

// ---------------------------------------------------------------------------
// Resolving a browser without taking a dependency on one
// ---------------------------------------------------------------------------

async function loadChromium() {
  const candidates = [
    process.env.PLAYWRIGHT_MODULE,
    join(APP_ROOT, "node_modules", "playwright", "index.mjs"),
    ...(process.env.HOME
      ? [join(process.env.HOME, ".npm/_npx/e41f203b7505f1fb/node_modules/playwright/index.mjs")]
      : []),
  ].filter(Boolean);
  for (const c of candidates) {
    if (!existsSync(c)) continue;
    try {
      return (await import(pathToFileURL(c).href)).chromium;
    } catch {
      /* try the next one */
    }
  }
  try {
    return require("playwright").chromium;
  } catch {
    /* fall through */
  }
  console.error(
    "viewport audit CANNOT RUN: no playwright found.\n" +
      "  Tried:\n" +
      candidates.map((c) => `    ${c}`).join("\n") +
      "\n  Set PLAYWRIGHT_MODULE=/path/to/playwright/index.mjs, or `npx playwright install chromium`.\n" +
      "  Not being able to run the check is not the same as the check passing.",
  );
  process.exit(2);
}

// ---------------------------------------------------------------------------
// Serving the BUILT page. A dev build lays out differently enough to lie.
// ---------------------------------------------------------------------------

async function startServer() {
  if (!existsSync(join(APP_ROOT, ".next", "BUILD_ID"))) {
    console.error(
      "viewport audit CANNOT RUN: no .next/BUILD_ID.\n" +
        "  NEXT_PUBLIC_STUDIO_URL=https://example.org npm run build   # then re-run\n" +
        "  A stale .next is worse than none: check the build's EXIT CODE, not its log.",
    );
    process.exit(2);
  }
  const child = spawn("npx", ["next", "start", "--port", String(PORT)], {
    cwd: APP_ROOT,
    env: {
      ...process.env,
      NEXT_PUBLIC_STUDIO_URL: process.env.NEXT_PUBLIC_STUDIO_URL ?? "https://example.org",
    },
    stdio: ["ignore", "ignore", "pipe"],
  });
  const base = `http://127.0.0.1:${PORT}`;
  const deadline = Date.now() + 60_000;
  for (;;) {
    if (Date.now() > deadline) {
      child.kill("SIGTERM");
      throw new Error("`next start` did not come up within 60s");
    }
    try {
      const r = await fetch(base);
      if (r.ok) break;
    } catch {
      /* not up yet */
    }
    await new Promise((r) => setTimeout(r, 400));
  }
  return { base, stop: () => child.kill("SIGTERM") };
}

// ---------------------------------------------------------------------------
// The two in-page measurements
// ---------------------------------------------------------------------------

/** Rects of every beat and both panels, relative to the pin that clips them. */
const MEASURE_PINNED = () => {
  const pin = document.querySelector(".scroll-pin");
  if (!pin) return null;
  const p = pin.getBoundingClientRect();
  const rel = (el) => {
    const r = el.getBoundingClientRect();
    return {
      opacity: Number(getComputedStyle(el).opacity),
      top: r.top - p.top,
      bottom: r.bottom - p.top,
      left: r.left,
      right: r.right,
    };
  };
  const out = { pin: p.height, sections: [], panels: {} };
  for (const [key, sel] of [
    ["attribution", ".scroll-pin > .attribution"],
    ["metric", ".scroll-pin > .metric-overlay"],
  ]) {
    const el = document.querySelector(sel);
    if (el) out.panels[key] = rel(el);
  }
  document.querySelectorAll(".narrative-layer > .narrative-section").forEach((el) => {
    // The BODY, not the section: the section is the full-height positioning box
    // and its padding is the fitting instrument, not the thing being fitted.
    out.sections.push(rel(el.querySelector(".narrative-section__body") ?? el));
  });
  return out;
};

/** The flow layout cannot clip vertically, so the question is width. */
const MEASURE_FLOW = () => {
  const de = document.documentElement;
  const overflow = [];
  document.querySelectorAll("main *").forEach((el) => {
    const r = el.getBoundingClientRect();
    if (r.width === 0 && r.height === 0) return;
    if (r.right > de.clientWidth + 0.5 || r.left < -0.5) {
      const name =
        typeof el.className === "string" && el.className ? el.className.split(" ")[0] : el.tagName;
      overflow.push({ name, left: Math.round(r.left), right: Math.round(r.right) });
    }
  });

  // The hero stage is the one thing in this layout that DOES clip.
  const stage = document.querySelector(".reduced-motion-experience__stage");
  const hero = document.querySelector(".reduced-motion-experience__hero-copy");
  let heroFit = null;
  if (stage && hero) {
    const sr = stage.getBoundingClientRect();
    let top = Infinity;
    let bottom = -Infinity;
    for (const n of hero.children) {
      const r = n.getBoundingClientRect();
      if (r.height === 0) continue;
      top = Math.min(top, r.top - sr.top);
      bottom = Math.max(bottom, r.bottom - sr.top);
    }
    heroFit = { stage: sr.height, top, bottom };
  }
  return {
    hOverflow: de.scrollWidth - de.clientWidth,
    overflow: overflow.slice(0, 12),
    heroFit,
  };
};

// ---------------------------------------------------------------------------

const r1 = (n) => Math.round(n * 10) / 10;

async function auditViewport(browser, base, vp, reduced) {
  const ctx = await browser.newContext({
    viewport: { width: vp.width, height: vp.height },
    reducedMotion: reduced ? "reduce" : "no-preference",
    deviceScaleFactor: 1,
  });
  const page = await ctx.newPage();
  try {
    await page.goto(base, { waitUntil: "networkidle" });
    // Fonts and the GLB settle a frame or two late; measuring before they do
    // measures the fallback metrics, which is a different page.
    await page.waitForTimeout(1500);

    const pinned = await page.evaluate(() => Boolean(document.querySelector(".scroll-pin")));
    if (!pinned) {
      return { vp: vp.label, reduced, mode: "flow", flow: await page.evaluate(MEASURE_FLOW) };
    }

    const scrollable = await page.evaluate(
      () => document.documentElement.scrollHeight - window.innerHeight,
    );
    const beats = new Map();
    const panels = new Map();

    for (let i = 0; i <= SAMPLES; i++) {
      await page.evaluate((y) => window.scrollTo(0, y), Math.round((scrollable * i) / SAMPLES));
      await page.waitForTimeout(12);
      const m = await page.evaluate(MEASURE_PINNED);
      if (!m) continue;

      m.sections.forEach((s, idx) => {
        if (s.opacity <= 0.02) return;
        const cur = beats.get(idx) ?? { idx, top: Infinity, bottom: -Infinity, pin: m.pin };
        cur.top = Math.min(cur.top, s.top);
        cur.bottom = Math.max(cur.bottom, s.bottom);
        cur.pin = m.pin;
        beats.set(idx, cur);
      });

      // The panel sits under the beat whose evidence it is: attribution under
      // 03, metric under 04 (lib/timeline.ts drives both windows).
      for (const [key, idx] of [["attribution", 2], ["metric", 3]]) {
        const pan = m.panels[key];
        const sec = m.sections[idx];
        if (!pan || !sec || pan.opacity <= 0.02 || sec.opacity <= 0.02) continue;
        const dy = Math.min(sec.bottom, pan.bottom) - Math.max(sec.top, pan.top);
        const dx = Math.min(sec.right, pan.right) - Math.max(sec.left, pan.left);
        // Positive only when the rectangles genuinely intersect. Negative is
        // the clearance, which is the number worth printing.
        const overlap = dy > 0 && dx > 0 ? Math.min(dy, dx) : -Math.max(-dy, -dx);
        const cur = panels.get(key) ?? { key, overlap: -Infinity };
        if (overlap > cur.overlap) Object.assign(cur, { overlap, top: pan.top, bottom: pan.bottom });
        panels.set(key, cur);
      }
    }

    const hOverflow = await page.evaluate(
      () => document.documentElement.scrollWidth - document.documentElement.clientWidth,
    );

    return {
      vp: vp.label,
      reduced,
      mode: "pinned",
      pin: vp.height,
      hOverflow,
      beats: [...beats.values()].map((b) => ({
        beat: String(b.idx + 1).padStart(2, "0"),
        top: r1(b.top),
        bottom: r1(b.bottom),
        clipTop: b.top < 0 ? r1(-b.top) : 0,
        clipBottom: b.bottom > b.pin ? r1(b.bottom - b.pin) : 0,
        spare: r1(b.pin - b.bottom),
      })),
      panels: [...panels.values()].map((p) => ({
        key: p.key,
        top: r1(p.top),
        bottom: r1(p.bottom),
        overlap: r1(p.overlap),
      })),
    };
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
      for (const reduced of [false, true]) {
        const r = await auditViewport(browser, base, vp, reduced);
        const tag = `${r.vp} ${reduced ? "reduced-motion" : "default"}`.padEnd(28);

        if (r.mode === "flow") {
          const h = r.flow.heroFit;
          const heroClipped = h && (h.bottom > h.stage + 0.5 || h.top < -0.5);
          const ok = r.flow.hOverflow <= 0 && r.flow.overflow.length === 0 && !heroClipped;
          console.log(
            `${tag} FLOW    hOverflow ${r.flow.hOverflow}px` +
              (h ? `  hero copy ${r1(h.top)}..${r1(h.bottom)} of stage ${r1(h.stage)}` : "") +
              (ok ? "" : "   <-- FAIL"),
          );
          for (const o of r.flow.overflow) {
            console.log(`      OVERFLOWS  ${o.name}  ${o.left}..${o.right}`);
          }
          if (heroClipped) console.log("      HERO COPY CLIPPED by the stage");
          if (!ok) failures.push(`${r.vp} ${reduced ? "reduced" : "default"} (flow)`);
          continue;
        }

        const clipped = r.beats.filter((b) => b.clipTop || b.clipBottom);
        const collided = r.panels.filter((p) => p.overlap > 0);
        const ok = clipped.length === 0 && collided.length === 0 && r.hOverflow <= 0;
        console.log(
          `${tag} PINNED  pin ${r.pin}px  hOverflow ${r.hOverflow}px  ` +
            (ok ? "all beats fit" : `${clipped.length} clipped, ${collided.length} collided  <-- FAIL`),
        );
        for (const b of r.beats) {
          const flag = b.clipTop || b.clipBottom
            ? `  <-- CLIPPED top ${b.clipTop} bottom ${b.clipBottom}`
            : `  (${b.spare}px spare)`;
          console.log(`      beat ${b.beat}  ${b.top} .. ${b.bottom}${flag}`);
        }
        for (const p of r.panels) {
          console.log(
            `      panel ${p.key.padEnd(11)} ${p.top} .. ${p.bottom}  ` +
              (p.overlap > 0
                ? `COLLIDES with copy by ${p.overlap}px  <-- FAIL`
                : `clears copy by ${-p.overlap}px`),
          );
        }
        if (!ok) failures.push(`${r.vp} ${reduced ? "reduced" : "default"} (pinned)`);
      }
    }
  } finally {
    await browser.close();
    stop();
  }

  console.log("");
  if (failures.length) {
    console.log(`viewport audit FAILED at ${failures.length} viewport(s):`);
    for (const f of failures) console.log(`  ${f}`);
    process.exit(1);
  }
  console.log(`viewport audit passed — ${VIEWPORTS.length} viewports x 2 motion settings`);
}

main().catch((e) => {
  console.error(e);
  process.exit(1);
});
