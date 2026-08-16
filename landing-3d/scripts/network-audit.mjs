#!/usr/bin/env node
/**
 * The execution-level proof behind the claim that this page contacts nothing
 * but its own server.
 *
 * `scripts/audit-origins.mjs` reads the build output as text. Text cannot tell
 * you whether a URL is fetched or merely printed in an error message, and it
 * cannot see a host assembled at runtime from fragments. So it settles a
 * narrower question: is the set of external hosts NAMED in the output the set
 * that has been reviewed.
 *
 * This file settles the actual question, by loading the built page in a real
 * browser and watching every request the browser makes.
 *
 * WHY IT SCROLLS THE WHOLE PAGE AND WAITS FOR THE MODEL
 * ====================================================
 * The bug this exists for was CONDITIONAL. drei's `useGLTF` attaches a
 * `DRACOLoader` pointed at `https://www.gstatic.com/draco/...`, and that
 * decoder is requested only when a glTF actually carries
 * `KHR_draco_mesh_compression`. An audit that loaded the page and looked at the
 * request log without the 3D scene ever having loaded would have reported a
 * clean result no matter what -- it would have been measuring the absence of a
 * request it never gave the page a chance to make.
 *
 * So this audit asserts, before it will report anything, that the conditional
 * work actually happened: the WebGL context is live, `/models/tumor.glb` came
 * back 200, the loading overlay dismissed, and the full scroll timeline ran.
 * A zero that was never at risk of being non-zero is not evidence, and
 * `assertNotVacuous` below refuses to let this print a pass without it.
 *
 * WHAT IT DOES NOT PROVE
 * ======================
 * A Draco-COMPRESSED asset is never decoded here, because producing one needs
 * a Draco encoder and none is vendored (the decoder is; see
 * public/draco/README.md). What is checked instead is that the pinned decoder
 * path is genuinely served by this origin -- all three files, 200, same origin
 * -- so a compressed re-export resolves locally. The end-to-end decode of a
 * compressed mesh is untested.
 *
 * USAGE
 *   node scripts/network-audit.mjs                 # builds nothing; starts `next start`
 *   node scripts/network-audit.mjs --url http://localhost:3000
 *   PLAYWRIGHT_MODULE=/path/to/playwright/index.mjs node scripts/network-audit.mjs
 *
 * Playwright is intentionally NOT a dependency of this package. It is a
 * ~100 MB browser download in service of a check that runs by hand and in CI,
 * and adding it would put a browser in the install path of a marketing page.
 * It is resolved from wherever it already exists.
 */

import { spawn } from "node:child_process";
import { createRequire } from "node:module";
import { dirname, join } from "node:path";
import { fileURLToPath, pathToFileURL } from "node:url";
import { existsSync } from "node:fs";

const APP_ROOT = join(dirname(fileURLToPath(import.meta.url)), "..");
const require = createRequire(import.meta.url);

const argv = process.argv.slice(2);
const urlArg = argv.find((a) => a.startsWith("--url="))?.slice(6);
const PORT = Number(argv.find((a) => a.startsWith("--port="))?.slice(7) ?? 3141);

/**
 * Both widths and both experiences. `prefers-reduced-motion` mounts
 * `ReducedMotionExperience`, a different component tree that still loads the
 * same GLB through the same loader -- so it is a second, independent chance
 * for a decoder request to fire, and auditing only the scroll experience would
 * leave it unmeasured.
 */
const VIEWPORTS = [
  { label: "cinematic 1512x950", width: 1512, height: 950, reduced: false },
  { label: "cinematic 390x844", width: 390, height: 844, reduced: false },
  { label: "reduced-motion 1512x950", width: 1512, height: 950, reduced: true },
  { label: "reduced-motion 390x844", width: 390, height: 844, reduced: true },
];

// ---------------------------------------------------------------------------
// Resolving a browser without taking a dependency on one
// ---------------------------------------------------------------------------

async function loadChromium() {
  const candidates = [
    process.env.PLAYWRIGHT_MODULE,
    join(APP_ROOT, "node_modules", "playwright", "index.mjs"),
    // The npx cache, where this repo's other browser checks have run from.
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
    "network audit CANNOT RUN: no playwright found.\n" +
      "  Tried:\n" +
      candidates.map((c) => `    ${c}`).join("\n") +
      "\n  Set PLAYWRIGHT_MODULE=/path/to/playwright/index.mjs, or `npx playwright install chromium`.\n" +
      "  Not being able to run the proof is not the same as the proof passing.",
  );
  process.exit(2);
}

// ---------------------------------------------------------------------------
// Serving the BUILT page (not the dev server)
// ---------------------------------------------------------------------------

async function startServer() {
  if (urlArg) return { base: urlArg, stop: async () => {} };
  if (!existsSync(join(APP_ROOT, ".next", "BUILD_ID"))) {
    console.error("network audit CANNOT RUN: no production build. Run `npm run build` first.");
    console.error("  The dev server is a different bundle; auditing it would prove nothing about what ships.");
    process.exit(2);
  }
  const base = `http://127.0.0.1:${PORT}`;
  const child = spawn("npx", ["next", "start", "-p", String(PORT)], {
    cwd: APP_ROOT,
    stdio: ["ignore", "pipe", "pipe"],
    env: process.env,
  });
  const deadline = Date.now() + 60_000;
  for (;;) {
    if (Date.now() > deadline) {
      child.kill("SIGKILL");
      throw new Error(`next start did not come up on ${base} within 60s`);
    }
    try {
      const r = await fetch(base, { method: "HEAD" });
      if (r.status < 500) break;
    } catch {
      /* not up yet */
    }
    await new Promise((r) => setTimeout(r, 300));
  }
  return {
    base,
    stop: async () => {
      child.kill("SIGTERM");
      await new Promise((r) => setTimeout(r, 400));
      if (!child.killed) child.kill("SIGKILL");
    },
  };
}

// ---------------------------------------------------------------------------
// The audit
// ---------------------------------------------------------------------------

/**
 * Schemes that are not network requests to a host at all. `blob:` matters
 * specifically: `DRACOLoader` and three's loaders run their workers from Blob
 * URLs, so a filter that only understood http(s) would either miss them or
 * miscount them as foreign.
 */
const NON_NETWORK = new Set(["data:", "blob:", "about:", "chrome-extension:", "file:"]);

async function auditViewport(chromium, base, viewport) {
  const origin = new URL(base).origin;
  const browser = await chromium.launch();
  const context = await browser.newContext({
    viewport: { width: viewport.width, height: viewport.height },
    deviceScaleFactor: 1,
    reducedMotion: viewport.reduced ? "reduce" : "no-preference",
  });

  /** Every request the browser attempts, from any frame, worker or service worker. */
  const requests = [];
  /**
   * Flipped once this audit starts probing `/draco/` itself. Those three
   * requests are the audit's, not the page's, and counting them as page
   * traffic would overstate what the page does.
   */
  let probing = false;
  const record = (req) => {
    let u;
    try {
      u = new URL(req.url());
    } catch {
      return;
    }
    requests.push({
      url: req.url(),
      origin: NON_NETWORK.has(u.protocol) ? u.protocol : u.origin,
      type: req.resourceType(),
      probe: probing,
      from: req.frame()?.url?.() ?? "(worker)",
    });
  };
  // Both, deduped by identity below: page-level catches subresources, the
  // context-level listener also catches workers and service workers.
  context.on("request", record);

  const page = await context.newPage();
  const consoleErrors = [];
  page.on("console", (m) => m.type() === "error" && consoleErrors.push(m.text().slice(0, 200)));
  page.on("pageerror", (e) => consoleErrors.push("pageerror: " + String(e).slice(0, 200)));

  const responses = new Map();
  page.on("response", (r) => responses.set(new URL(r.url()).pathname, r.status()));

  await page.goto(base, { waitUntil: "networkidle", timeout: 120_000 });
  await page.evaluate(() => document.fonts.ready);

  // Give the availability HEAD check, the Suspense swap and the GLB parse time.
  await page.waitForTimeout(2500);

  // Drive the whole scroll timeline: every beat, every lazily-mounted overlay,
  // the region separation, the drug fronts. This is the "interacted with" part
  // -- the request under test only fires under a condition, and the condition
  // is reached by using the page.
  await page.evaluate(async () => {
    const max = document.body.scrollHeight - innerHeight;
    for (let i = 0; i <= 40; i++) {
      window.scrollTo(0, (max * i) / 40);
      await new Promise((r) => requestAnimationFrame(() => requestAnimationFrame(r)));
      await new Promise((r) => setTimeout(r, 40));
    }
    window.scrollTo(0, 0);
    await new Promise((r) => setTimeout(r, 200));
  });
  await page.waitForTimeout(1500);

  // Evidence that the conditional work actually happened.
  const scene = await page.evaluate(() => {
    const canvas = document.querySelector("canvas");
    const gl = canvas && (canvas.getContext("webgl2") || canvas.getContext("webgl"));
    return {
      canvas: !!canvas,
      canvasPixels: canvas ? canvas.width * canvas.height : 0,
      webgl: !!gl,
      renderer: gl
        ? (() => {
            const ext = gl.getExtension("WEBGL_debug_renderer_info");
            return ext ? gl.getParameter(ext.UNMASKED_RENDERER_WEBGL) : gl.getParameter(gl.RENDERER);
          })()
        : null,
      loaderGone: !document.querySelector(".scene-loader:not(.scene-loader--done)"),
      scrollHeight: document.body.scrollHeight,
    };
  });

  // The pinned decoder path must actually be served by this origin. Fetched
  // from inside the page so it goes through the same origin the page has.
  probing = true;
  const decoder = await page.evaluate(async () => {
    const out = {};
    for (const f of ["draco_wasm_wrapper.js", "draco_decoder.wasm", "draco_decoder.js"]) {
      try {
        const r = await fetch(`/draco/${f}`, { method: "GET" });
        out[f] = { status: r.status, sameOrigin: new URL(r.url).origin === location.origin };
      } catch (e) {
        out[f] = { status: "threw", error: String(e) };
      }
    }
    return out;
  });

  await context.close();
  await browser.close();
  return { origin, requests, consoleErrors, responses, scene, decoder };
}

/**
 * Refuse to report a pass on a run that could not have failed.
 *
 * If the scene never loaded, "zero external requests" is a statement about a
 * page that never did the thing under test.
 */
function assertNotVacuous({ scene, responses, requests }) {
  const problems = [];
  if (!scene.canvas) problems.push("no <canvas> was rendered -- the 3D scene never mounted");
  if (!scene.webgl) problems.push("the canvas has no live WebGL context -- nothing was drawn");
  if (!scene.canvasPixels) problems.push("the canvas has zero pixels");
  if (!scene.loaderGone) problems.push("the loading overlay never dismissed -- assets did not finish");

  const glb = responses.get("/models/tumor.glb");
  if (glb !== 200) {
    problems.push(
      `/models/tumor.glb did not load (status ${glb ?? "never requested"}). That GLB is the ` +
        `asset whose properties decide whether a decoder is fetched; without it this audit ` +
        `is measuring the absence of a request the page never had cause to make.`,
    );
  }
  const fonts = requests.filter((r) => /\.woff2?($|\?)/.test(r.url));
  if (!fonts.length) {
    problems.push(
      "no webfont was requested. next/font/google is supposed to self-host and serve woff2 " +
        "from this origin; zero font requests means the self-hosting path was not exercised.",
    );
  }
  return problems;
}

// ---------------------------------------------------------------------------

const chromium = await loadChromium();
const { base, stop } = await startServer();
let exitCode = 0;

try {
  console.log(`network audit -- built page at ${base}`);
  console.log(`asserting: zero requests to any host but the page's own origin\n`);

  for (const viewport of VIEWPORTS) {
    const result = await auditViewport(chromium, base, viewport);
    const { origin, requests, consoleErrors, scene, decoder } = result;

    const byOrigin = new Map();
    for (const r of requests) {
      if (!byOrigin.has(r.origin)) byOrigin.set(r.origin, []);
      byOrigin.get(r.origin).push(r);
    }
    const foreign = [...byOrigin.entries()].filter(
      ([o]) => o !== origin && !NON_NETWORK.has(o),
    );

    const pageRequests = requests.filter((r) => !r.probe);
    console.log(`### ${viewport.label}`);
    console.log(`    renderer          ${scene.renderer}`);
    console.log(`    page height       ${scene.scrollHeight}px, scrolled end to end`);
    console.log(
      `    requests observed ${requests.length} ` +
        `(${pageRequests.length} by the page, ${requests.length - pageRequests.length} this audit's /draco/ probe)`,
    );
    for (const [o, rs] of [...byOrigin].sort()) {
      const tag = o === origin ? "OWN ORIGIN" : NON_NETWORK.has(o) ? "not a network request" : "FOREIGN";
      console.log(`      ${String(rs.length).padStart(4)}  ${o.padEnd(24)} ${tag}`);
    }
    // Print the log rather than only its size. A count is a summary; the list
    // is the evidence, and it is short enough to read.
    console.log(`    every request the page made:`);
    for (const r of pageRequests) {
      console.log(`      ${r.type.padEnd(10)} ${r.url.replace(origin, "")}`);
    }

    const vacuous = assertNotVacuous(result);
    if (vacuous.length) {
      exitCode = 1;
      console.log(`    VACUOUS -- this run proves nothing:`);
      for (const p of vacuous) console.log(`      - ${p}`);
    } else {
      console.log(
        `    not vacuous       canvas ${scene.canvasPixels}px, WebGL live, ` +
          `/models/tumor.glb 200, loader dismissed`,
      );
    }

    const decoderBad = Object.entries(decoder).filter(
      ([, v]) => v.status !== 200 || !v.sameOrigin,
    );
    if (decoderBad.length) {
      exitCode = 1;
      console.log(`    FAIL  vendored decoder not served from this origin at /draco/:`);
      for (const [f, v] of decoderBad) console.log(`      - ${f}: ${JSON.stringify(v)}`);
    } else {
      console.log(`    decoder           /draco/ x3 served 200 from own origin`);
    }

    if (consoleErrors.length) {
      exitCode = 1;
      console.log(`    FAIL  ${consoleErrors.length} console/page error(s):`);
      for (const e of consoleErrors.slice(0, 6)) console.log(`      - ${e}`);
    }

    if (foreign.length) {
      exitCode = 1;
      console.log(`    FAIL  ${foreign.length} foreign origin(s):`);
      for (const [o, rs] of foreign) {
        console.log(`      ${o}`);
        for (const r of rs.slice(0, 6)) console.log(`        ${r.type.padEnd(10)} ${r.url}`);
      }
    } else {
      console.log(`    PASS  zero requests to any origin but ${origin}`);
    }
    console.log();
  }
} finally {
  await stop();
}

console.log(exitCode === 0 ? "network audit passed" : "network audit FAILED");
process.exit(exitCode);
