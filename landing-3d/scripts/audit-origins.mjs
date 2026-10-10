#!/usr/bin/env node
/**
 * The guard against a published privacy claim that quietly stops being true.
 *
 * WHAT IS BEING PROTECTED
 * =======================
 * `lib/copy.ts` states, as a product claim, that the solve "sends no tissue or
 * compound data anywhere", and the page itself is meant to contact nothing but
 * its own origin. That was true when someone audited it, and it was true for a
 * reason nobody was checking: `@react-three/drei`'s `useGLTF` points a
 * `DRACOLoader` at `https://www.gstatic.com/draco/versioned/decoders/1.5.5/`,
 * and that decoder is fetched ONLY when a glTF actually carries
 * `KHR_draco_mesh_compression`. `public/models/tumor.glb` carries none, so the
 * request never fired.
 *
 * Re-export a 660 KB mesh with Draco compression -- an ordinary size
 * optimisation -- and the page keeps rendering, nothing errors, and the claim
 * is false. That is the shape of failure this repository keeps hitting: the
 * engine stays sound while a claim surface drifts out from under it. The fix
 * (vendoring the decoder into `public/draco/`, loading through
 * `components/canvas/tumorGLTF.ts` instead of drei) removes today's instance.
 * This file is what stops the NEXT one -- an analytics snippet, a CDN font, a
 * `next/image` remote loader, a re-introduced drei default -- from landing
 * unnoticed.
 *
 * WHAT THIS FILE CAN AND CANNOT PROVE
 * ===================================
 * It reads build output as text. Text cannot tell you whether
 * `"https://example.com"` is passed to `fetch()` or printed in an error
 * message, and it cannot see a host assembled at runtime from fragments. So
 * the honest statement of what this guard does is:
 *
 *   IT IS A CHANGE DETECTOR, NOT A PROOF.
 *
 * It proves that the set of external hosts named anywhere in the shipped
 * output is exactly the set that has been looked at and written down. Every
 * one of those hosts carries a written reason for being tolerable. A new one
 * -- from any source -- fails the build and has to be justified by a human.
 *
 * The PROOF that none of them is actually contacted is a different artefact:
 * `scripts/network-audit.mjs` loads the built page in a real browser, drives
 * the scroll so the 3D scene fully loads, and asserts zero requests to any
 * host but the page's own origin. Static "this string is not here" and runtime
 * "this request did not happen" are different claims about different things.
 * Collapsing them would be the same error this repo's
 * `tests/test_propagation.py` exists to prevent (REACHED_UNCHANGED is not
 * NO_PDE), so they are kept apart and both are run.
 *
 * LAYOUT (after tests/test_propagation.py)
 * ========================================
 *   1. Proof of correctness -- the classifier is run against a fixture of
 *      EVERY category, including the pair that look alike from the outside
 *      (a forbidden CDN host vs a doc link in an error string), before it is
 *      allowed to judge the real build. A guard nobody has watched fail is
 *      not a guard, so it is made to fail here on every run.
 *   2. The declared table -- every host in the real output lands in a
 *      declared entry, and every entry carries a written reason.
 *   3. The sweep -- nothing undeclared anywhere in the shipped output.
 *   4. The asset checks -- the vendored decoder is present and unmodified,
 *      the loader is still pinned at it, and no model declares an extension
 *      whose decoder we do not serve ourselves.
 *
 * Runs in well under a second, reads only files, executes nothing it reads.
 */

import { createHash } from "node:crypto";
import { readFileSync, readdirSync, existsSync, statSync } from "node:fs";
import { join, extname, relative, dirname, sep } from "node:path";
import { fileURLToPath } from "node:url";

const APP_ROOT = join(dirname(fileURLToPath(import.meta.url)), "..");

// --------------------------------------------------------------------------
// What counts as "shipped output"
// --------------------------------------------------------------------------

/**
 * Directories whose contents a browser receives. `.next/static` and
 * `.next/server` are the built client chunks and the prerendered HTML;
 * `public` is served verbatim at the origin root.
 *
 * `.next/dev` is deliberately absent. It is the dev server's own unminified
 * build, it is never deployed, and it lags behind source until the dev server
 * rebuilds -- auditing it would report stale findings about code that no
 * longer exists.
 */
const SHIPPED_DIRS = [".next/static", ".next/server", "public"];

/**
 * Extensions a browser PARSES AND ACTS ON as part of the page. A URL inside
 * one of these can become a request; a URL in anything else cannot without a
 * human first navigating to that file directly.
 */
const PARSED_EXTENSIONS = new Set([".js", ".mjs", ".cjs", ".css", ".html", ".json", ".rsc"]);

/**
 * Declared exclusions. Not silent skips -- each says what it is and why a URL
 * inside it cannot become a page request.
 */
const EXCLUDED = [
  {
    match: (p) => extname(p) === ".map",
    reason:
      "Source maps. Never referenced by the page's own execution; a devtools " +
      "user fetches them explicitly. They also embed the entire original " +
      "source, so every URL in a comment would surface here as a finding.",
  },
  {
    match: (p) => extname(p) === ".md" || p.endsWith(`${sep}LICENSE`),
    reason:
      "Documentation shipped alongside a vendored asset (public/draco/). " +
      "Served if requested, but the browser never parses it as part of the " +
      "page. public/draco/README.md necessarily NAMES the gstatic URL in " +
      "order to explain why it is no longer reached; auditing prose would " +
      "make an explanation of the fix indistinguishable from the bug.",
  },
];

// --------------------------------------------------------------------------
// Extraction
// --------------------------------------------------------------------------

/**
 * Absolute `http(s)://host` literals.
 *
 * Deliberately loose on what follows the host: the point is to enumerate
 * HOSTS, not to validate URLs.
 */
const ABSOLUTE = /https?:\/\/([A-Za-z0-9._-]+)/g;

/**
 * Protocol-relative `//host/...`, which resolves against the page's scheme and
 * is a completely ordinary way to reference a CDN -- `<script src="//cdn...">`.
 * A pattern anchored on `https?:` alone would miss it entirely.
 *
 * Anchored on a preceding string/attribute delimiter so that `a//b` in
 * minified code, and the `//` of a URL path, do not register. Requires a dot
 * in the host so bare `//foo` comment starts do not either.
 */
const PROTOCOL_RELATIVE = /["'`(=,\s]\/\/([A-Za-z0-9-]+(?:\.[A-Za-z0-9-]+)+)/g;

/**
 * Every host named in `text`, with a snippet of surrounding context for each
 * so a failure report can be acted on without a second search.
 */
function hostsIn(text, label) {
  const found = new Map();
  for (const re of [ABSOLUTE, PROTOCOL_RELATIVE]) {
    re.lastIndex = 0;
    for (const m of text.matchAll(re)) {
      const host = m[1].toLowerCase().replace(/\.$/, "");
      if (!found.has(host)) found.set(host, []);
      const sites = found.get(host);
      if (sites.length < 2) {
        const at = m.index ?? 0;
        sites.push({
          file: label,
          context: text
            .slice(Math.max(0, at - 60), at + 100)
            .replace(/\s+/g, " ")
            .trim(),
        });
      }
    }
  }
  return found;
}

function mergeHosts(into, from) {
  for (const [host, sites] of from) {
    if (!into.has(host)) into.set(host, []);
    const acc = into.get(host);
    for (const s of sites) if (acc.length < 4) acc.push(s);
  }
}

function walkShipped(root) {
  const out = [];
  const stack = [root];
  while (stack.length) {
    const dir = stack.pop();
    for (const entry of readdirSync(dir, { withFileTypes: true })) {
      const p = join(dir, entry.name);
      if (entry.isDirectory()) {
        stack.push(p);
      } else if (entry.isFile()) {
        if (EXCLUDED.some((e) => e.match(p))) continue;
        if (!PARSED_EXTENSIONS.has(extname(p))) continue;
        out.push(p);
      }
    }
  }
  return out;
}

// --------------------------------------------------------------------------
// The declared table
// --------------------------------------------------------------------------

/**
 * Hosts whose PRESENCE IS THE FAILURE. These are the CDN defaults that
 * three.js-adjacent libraries and web fonts reach for on their own, plus the
 * usual analytics endpoints. No entry here can be exempted: if one of these
 * strings is in the output, something in the dependency graph is prepared to
 * fetch from it, and the claim is at best conditionally true again.
 *
 * `www.gstatic.com` heads the list because it is the one that was actually
 * here. It is gone because `components/canvas/tumorGLTF.ts` loads through
 * `useLoader` rather than drei's `useGLTF`, so drei's `Gltf` module -- which
 * initialises its `decoderPath` variable to that literal -- is never imported.
 * Repointing the path with `useGLTF.setDecoderPath` would have fixed the
 * runtime behaviour while leaving the string in the bundle, and then this
 * entry would have had to become an exemption instead of a prohibition. An
 * exemption is a hole a real regression can hide in.
 *
 * THE TWO FONT ENTRIES' REASONS WERE REWRITTEN ON 2026-10-10, AND THE RULE WAS
 * NOT TOUCHED. They used to read "next/font/google self-hosts at build time; a
 * runtime hit here means it stopped", which described a page that no longer
 * exists: `app/layout.tsx` moved to `next/font/local` against four .woff2 files
 * committed in `app/fonts/`, so Google is not contacted at BUILD time either
 * and "a runtime hit" is no longer the only thing these strings could mean.
 * Leaving the old reason in place would have told whoever next sees this fail
 * to look for a runtime fetch, when the likelier cause is now a
 * `next/font/google` import coming back — which would also silently un-pin the
 * outlines, on a page whose tightest beat has 8.9px of vertical headroom
 * (`npm run audit:viewports`, 360x700 beat 06), and ship the binaries without
 * the OFL licence text that app/fonts/ carries. Same verdict, same
 * no-exemption rule, accurate instruction.
 */
const FORBIDDEN = [
  ["www.gstatic.com", "drei's default DRACO decoder path, and Google's asset CDN generally."],
  [
    "fonts.gstatic.com",
    "Google Fonts binaries. All four faces are vendored in app/fonts/ and loaded " +
      "through next/font/local, so nothing reaches Google at build time either -- " +
      "this string means a next/font/google import or a hand-written @font-face src " +
      "came back, not merely that a runtime fetch slipped through.",
  ],
  [
    "fonts.googleapis.com",
    "Google Fonts CSS. Same provenance, same conclusion: the faces are local, so " +
      "this appearing means a stylesheet @import or a next/font/google call was " +
      "reintroduced.",
  ],
  ["cdn.jsdelivr.net", "drei's default KTX2/Basis transcoder path, and a general-purpose CDN."],
  ["unpkg.com", "General-purpose npm CDN; a common default for WASM and decoder assets."],
  ["cdnjs.cloudflare.com", "General-purpose script CDN."],
  ["www.google-analytics.com", "Analytics beacon."],
  ["www.googletagmanager.com", "Tag manager; loads arbitrary third-party script."],
  ["vitals.vercel-insights.com", "Vercel analytics beacon."],
  ["vercel.live", "Vercel toolbar/feedback widget; injected by some deployment settings."],
];

/**
 * Kinds a tolerable host can be. These are kept distinct rather than merged
 * into one "allowed" bucket because they are different claims: "this URL is
 * inside an English sentence" and "a human may click this URL" are not the
 * same statement, and reviewing them as if they were is how a real request
 * gets waved through.
 */
const KIND = {
  /** Named only inside human-readable text: an error message, a code comment,
   *  a licence field. Nothing dereferences it. */
  DIAGNOSTIC: "DIAGNOSTIC",
  /** An XML/SVG namespace URI. Identifiers by specification; no conforming
   *  browser fetches them. */
  NAMESPACE: "NAMESPACE",
  /** A synthetic base handed to `new URL(relative, base)`, or a literal in a
   *  spec-conformance feature test. Parsed, never opened. */
  URL_PARSE_BASE: "URL_PARSE_BASE",
  /** The target of an anchor a person clicks. Not an automatic request: the
   *  page issues nothing, and a navigation is a decision the reader makes. */
  NAVIGATION: "NAVIGATION",
};

class Declaration {
  constructor(host, kind, reason) {
    if (!Object.values(KIND).includes(kind)) {
      throw new Error(`${host}: unknown kind ${kind}`);
    }
    if (reason.trim().length < 40) {
      throw new Error(`${host}: a declaration needs a written reason`);
    }
    this.host = host;
    this.kind = kind;
    this.reason = reason;
  }
}

const DECLARED = [
  new Declaration(
    "www.w3.org",
    KIND.NAMESPACE,
    "XHTML and SVG namespace URIs, from GSAP's SVG helpers and React's DOM " +
      "creation path. Namespace URIs are identifiers, not locations; nothing " +
      "resolves them.",
  ),
  new Declaration(
    "nextjs.org",
    KIND.DIAGNOSTIC,
    "Documentation links inside Next.js runtime error messages " +
      "(nextjs.org/docs/messages/...). Text in a thrown Error.",
  ),
  new Declaration(
    "react.dev",
    KIND.DIAGNOSTIC,
    "React's minified-error decoder link, concatenated into the message of an " +
      "Error that React throws. Text, not a target.",
  ),
  new Declaration(
    "docs.pmnd.rs",
    KIND.DIAGNOSTIC,
    "react-three-fiber's 'did you forget to extend?' developer warning points " +
      "at its own docs. Text in a console warning.",
  ),
  new Declaration(
    "gsap.com",
    KIND.DIAGNOSTIC,
    "GSAP's 'target not found' console warning includes its homepage. Text in " +
      "a console warning.",
  ),
  new Declaration(
    "github.com",
    KIND.DIAGNOSTIC,
    "core-js stamps its own licence and source URLs into a metadata object it " +
      "exposes for tooling. Strings in an object literal; core-js never " +
      "requests them.",
  ),
  new Declaration(
    "jcgt.org",
    KIND.DIAGNOSTIC,
    "A paper citation in a GLSL comment inside three.js's PBR shader chunks " +
      "(the GGX VNDF sampling reference). Inside shader source text.",
  ),
  new Declaration(
    "datatracker.ietf.org",
    KIND.DIAGNOSTIC,
    "An RFC 3986 citation in a comment in Turbopack's runtime, explaining why " +
      "it trims a query string before a fragment.",
  ),
  new Declaration(
    "a",
    KIND.URL_PARSE_BASE,
    "core-js's URL/URLSearchParams conformance tests construct new URL(...) " +
      "against the placeholder origin 'https://a' and compare the parse " +
      "result. Never opened; also not a resolvable host.",
  ),
  new Declaration(
    "x",
    KIND.URL_PARSE_BASE,
    "Same core-js conformance suite, checking that new URL('https://x').host " +
      "parses to 'x'. A parser assertion.",
  ),
  new Declaration(
    "n",
    KIND.URL_PARSE_BASE,
    "Next.js's client router resolves relative hrefs against a throwaway base " +
      "('http://n') when there is no window, i.e. during server rendering, so " +
      "that URL parsing works. Discarded after parsing.",
  ),
];

/**
 * The Studio link is the one host on the page a request can legitimately reach,
 * and only because a reader clicks it. It is deployment-configured, so its
 * host is read from the same env var the build inlines rather than hardcoded
 * -- a hardcoded value would either fail every non-default deployment or, if
 * loosened to "any host", stop guarding anything.
 *
 * It is a NAVIGATION, and the distinction is the whole point: following a link
 * carries no tissue or compound data anywhere, because the page never had any.
 * Sections register by filesystem path; there is nothing here to send.
 */
function studioDeclaration() {
  const raw = process.env.NEXT_PUBLIC_STUDIO_URL ?? "http://localhost:8000";
  let host;
  try {
    host = new URL(raw).hostname.toLowerCase();
  } catch {
    throw new Error(
      `NEXT_PUBLIC_STUDIO_URL is not a URL: ${JSON.stringify(raw)}. ` +
        `It is inlined into the bundle as an href, so it has to parse.`,
    );
  }
  return new Declaration(
    host,
    KIND.NAVIGATION,
    `The "Open the Studio" call to action, from NEXT_PUBLIC_STUDIO_URL ` +
      `(${raw}). An <a href>: reached only when a reader clicks it, and it is ` +
      `the product this page is about. Read from the environment so this ` +
      `audit tracks whatever the deployment actually shipped.`,
  );
}

// --------------------------------------------------------------------------
// Classification
// --------------------------------------------------------------------------

const VERDICT = { FORBIDDEN: "FORBIDDEN", UNDECLARED: "UNDECLARED", DECLARED: "DECLARED" };

function classify(host, declarations) {
  const forbidden = FORBIDDEN.find(([h]) => h === host);
  if (forbidden) return { verdict: VERDICT.FORBIDDEN, note: forbidden[1] };
  const declared = declarations.find((d) => d.host === host);
  if (declared) return { verdict: VERDICT.DECLARED, note: declared.kind };
  return { verdict: VERDICT.UNDECLARED, note: "" };
}

// --------------------------------------------------------------------------
// 1. Proof of correctness
// --------------------------------------------------------------------------

/**
 * Run the extractor and the classifier over a fixture of every category before
 * letting them judge the real build.
 *
 * The two that matter most are the two that look identical from the outside: a
 * gstatic decoder path and a github.com link in a licence string are both
 * "an external URL in a JS chunk". One is the bug and one is fine, and a check
 * that cannot tell them apart is either useless or a nuisance. Both are here.
 *
 * This also means the failure path executes on every single run. A guard whose
 * failure branch has never been observed is a guard nobody has tested.
 */
function proveClassifier() {
  const declarations = [...DECLARED, studioDeclaration()];
  const cases = [
    {
      name: "drei's DRACO default -- the actual bug this file exists for",
      text: `let decoderPath='https://www.gstatic.com/draco/versioned/decoders/1.5.5/';`,
      host: "www.gstatic.com",
      expect: VERDICT.FORBIDDEN,
    },
    {
      name: "a CDN font added to the stylesheet",
      text: `@import url("https://fonts.googleapis.com/css2?family=Inter");`,
      host: "fonts.googleapis.com",
      expect: VERDICT.FORBIDDEN,
    },
    {
      name: "an analytics snippet",
      text: `<script src="https://www.googletagmanager.com/gtag/js?id=G-X"></script>`,
      host: "www.googletagmanager.com",
      expect: VERDICT.FORBIDDEN,
    },
    {
      name: "a protocol-relative CDN reference, which an https?: pattern misses",
      text: `<script src="//cdn.jsdelivr.net/npm/thing/dist/thing.min.js"></script>`,
      host: "cdn.jsdelivr.net",
      expect: VERDICT.FORBIDDEN,
    },
    {
      name: "an external host nobody has looked at yet",
      text: `fetch("https://telemetry.example.net/collect",{method:"POST"})`,
      host: "telemetry.example.net",
      expect: VERDICT.UNDECLARED,
    },
    {
      name: "a doc link in an error string -- looks the same, is not the same",
      text: `{license:"https://github.com/zloirock/core-js/blob/v3/LICENSE"}`,
      host: "github.com",
      expect: VERDICT.DECLARED,
    },
    {
      name: "an SVG namespace URI",
      text: `createElementNS("http://www.w3.org/2000/svg","svg")`,
      host: "www.w3.org",
      expect: VERDICT.DECLARED,
    },
    {
      name: "a synthetic base for new URL()",
      text: `new URL(f, typeof window>"u"?"http://n":window.location.href)`,
      host: "n",
      expect: VERDICT.DECLARED,
    },
  ];

  const failures = [];
  for (const c of cases) {
    const hosts = hostsIn(c.text, "<fixture>");
    if (!hosts.has(c.host)) {
      failures.push(`extractor missed ${c.host} in fixture: ${c.name}`);
      continue;
    }
    const got = classify(c.host, declarations).verdict;
    if (got !== c.expect) {
      failures.push(`${c.name}: expected ${c.expect} for ${c.host}, got ${got}`);
    }
  }

  // The extractor must NOT invent hosts out of same-origin references, or
  // every build would fail and the guard would be turned off within a week.
  const quiet = [
    `loader.setDecoderPath("/draco/")`,
    `<img src="/models/preview.png">`,
    `background:url(data:image/png;base64,iVBORw0KGgo=)`,
    `const a = 10 //2; const b = c//d;`,
  ];
  for (const text of quiet) {
    const hosts = hostsIn(text, "<fixture>");
    if (hosts.size > 0) {
      failures.push(
        `extractor invented ${[...hosts.keys()].join(", ")} from a same-origin ` +
          `or inert reference: ${text}`,
      );
    }
  }

  // A declaration without a written reason must be impossible to add.
  let rejected = false;
  try {
    new Declaration("example.com", KIND.DIAGNOSTIC, "fine");
  } catch {
    rejected = true;
  }
  if (!rejected) failures.push("Declaration accepted an unwritten reason");

  return failures;
}

// --------------------------------------------------------------------------
// 4. Asset checks
// --------------------------------------------------------------------------

/**
 * The vendored decoder, by digest. Copied from
 * `node_modules/three/examples/jsm/libs/draco/gltf/` as shipped in
 * three@0.185.1; Apache-2.0. See `public/draco/README.md`.
 *
 * Checked by content rather than by existence: "a file is at that path" and
 * "the decoder we reviewed and licensed is at that path" are different
 * statements, and only the second one is what the origin promises to serve.
 */
const VENDORED_DRACO = {
  "draco_wasm_wrapper.js": "8bb2952d2ba7d67e1414f8df819410cb0434a666be53f671fff75f68843d76f6",
  "draco_decoder.wasm": "a680d927bed9cb864ddbd63521868891af2bfbe755092761b4837487618df8ac",
  "draco_decoder.js": "8625489da79a805f4f2a7d511c3e52d8b4085608a9d2a4d5f4f9de5db0aea04f",
};

function checkVendoredDecoder() {
  const failures = [];
  const dir = join(APP_ROOT, "public", "draco");
  for (const [name, expected] of Object.entries(VENDORED_DRACO)) {
    const p = join(dir, name);
    if (!existsSync(p)) {
      failures.push(
        `public/draco/${name} is missing. The DRACO decoder is served from ` +
          `this origin so that a Draco-compressed model does not send anyone ` +
          `to a CDN. Re-copy it from ` +
          `node_modules/three/examples/jsm/libs/draco/gltf/.`,
      );
      continue;
    }
    const actual = createHash("sha256").update(readFileSync(p)).digest("hex");
    if (actual !== expected) {
      failures.push(
        `public/draco/${name} does not match its recorded digest.\n` +
          `      expected ${expected}\n` +
          `      actual   ${actual}\n` +
          `      If this was a deliberate three.js upgrade, update the digest ` +
          `here and in public/draco/README.md. Otherwise this origin is ` +
          `serving a decoder nobody reviewed.`,
      );
    }
  }
  return failures;
}

/**
 * The decoder is only reached from this origin if the loader is still pointed
 * at it. Vendoring the files and then loading through something that resolves
 * elsewhere would leave every file check passing and the claim false, so the
 * pin is asserted in source too.
 */
function checkDecoderPin() {
  const failures = [];
  const p = join(APP_ROOT, "components", "canvas", "tumorGLTF.ts");
  if (!existsSync(p)) {
    return [
      `components/canvas/tumorGLTF.ts is gone. That module is what keeps ` +
        `glTF loading off drei's CDN defaults; without it useGLTF's gstatic ` +
        `path comes straight back.`,
    ];
  }
  const src = readFileSync(p, "utf8");
  const m = src.match(/DRACO_DECODER_PATH\s*=\s*["'`]([^"'`]+)["'`]/);
  if (!m) {
    failures.push(`tumorGLTF.ts no longer defines DRACO_DECODER_PATH.`);
  } else if (!m[1].startsWith("/") || m[1].startsWith("//")) {
    failures.push(
      `DRACO_DECODER_PATH is ${JSON.stringify(m[1])}, which is not a ` +
        `same-origin absolute path. It must resolve against this origin.`,
    );
  } else if (m[1] !== "/draco/") {
    failures.push(
      `DRACO_DECODER_PATH is ${JSON.stringify(m[1])} but the decoder is ` +
        `vendored at public/draco/. The loader would 404.`,
    );
  }

  // A component reaching for drei's useGLTF re-imports the gstatic default.
  // The host sweep would catch it, but this says so in the language of the
  // mistake rather than as a mysterious string finding.
  for (const dir of ["components", "app", "hooks", "lib"]) {
    const root = join(APP_ROOT, dir);
    if (!existsSync(root)) continue;
    const stack = [root];
    while (stack.length) {
      const d = stack.pop();
      for (const e of readdirSync(d, { withFileTypes: true })) {
        const fp = join(d, e.name);
        if (e.isDirectory()) {
          stack.push(fp);
          continue;
        }
        if (![".ts", ".tsx"].includes(extname(fp))) continue;
        const s = readFileSync(fp, "utf8");
        if (/\buseGLTF\b|\bGltf\b\s*[,}]/.test(s) && /@react-three\/drei/.test(s)) {
          failures.push(
            `${relative(APP_ROOT, fp)} imports useGLTF from @react-three/drei. ` +
              `That attaches a DRACOLoader pointed at ` +
              `https://www.gstatic.com/draco/... Use useTumorGLTF from ` +
              `components/canvas/tumorGLTF.ts, which is pinned at /draco/.`,
          );
        }
      }
    }
  }
  return failures;
}

/**
 * glTF extensions that need a DECODER, and where ours comes from. Extensions
 * outside this table (materials, texture transforms, lights) are handled by
 * three.js itself and fetch nothing, so they are not this check's business.
 *
 * This is the check that closes the original trap from the asset side: a model
 * declaring an extension whose only decoder is a CDN default fails here, at
 * build time, with the reason -- instead of shipping and quietly issuing the
 * request.
 */
const DECODER_EXTENSIONS = {
  KHR_draco_mesh_compression: {
    ok: true,
    how: "Decoder vendored at public/draco/ and pinned in components/canvas/tumorGLTF.ts.",
  },
  EXT_meshopt_compression: {
    ok: true,
    how: "three-stdlib's MeshoptDecoder carries its WASM inline as base64 and fetches nothing; wired up in components/canvas/tumorGLTF.ts.",
  },
  KHR_texture_basisu: {
    ok: false,
    how:
      "Needs a KTX2/Basis transcoder. three's KTX2Loader has no transcoder " +
      "until one is set, and drei's default points at cdn.jsdelivr.net. " +
      "Nothing on this page is textured, so rather than vendor a transcoder " +
      "for a case that does not exist, this fails: either export without " +
      "Basis textures, or vendor the transcoder into public/ the way " +
      "public/draco/ was done and wire it up in tumorGLTF.ts.",
  },
};

function glbExtensions(buf) {
  if (buf.length < 20 || buf.readUInt32LE(0) !== 0x46546c67) return null; // 'glTF'
  const chunkLength = buf.readUInt32LE(12);
  const chunkType = buf.readUInt32LE(16);
  if (chunkType !== 0x4e4f534a) return null; // 'JSON'
  const json = JSON.parse(buf.subarray(20, 20 + chunkLength).toString("utf8"));
  return [
    ...new Set([...(json.extensionsUsed ?? []), ...(json.extensionsRequired ?? [])]),
  ];
}

function checkModelExtensions() {
  const failures = [];
  const dir = join(APP_ROOT, "public", "models");
  if (!existsSync(dir)) return failures;
  for (const name of readdirSync(dir)) {
    if (![".glb", ".gltf"].includes(extname(name))) continue;
    const p = join(dir, name);
    let exts;
    try {
      exts =
        extname(name) === ".glb"
          ? glbExtensions(readFileSync(p))
          : [
              ...new Set([
                ...(JSON.parse(readFileSync(p, "utf8")).extensionsUsed ?? []),
                ...(JSON.parse(readFileSync(p, "utf8")).extensionsRequired ?? []),
              ]),
            ];
    } catch (err) {
      failures.push(`public/models/${name}: could not read glTF header (${err.message})`);
      continue;
    }
    if (exts === null) {
      failures.push(`public/models/${name}: not a readable GLB container.`);
      continue;
    }
    for (const ext of exts) {
      const known = DECODER_EXTENSIONS[ext];
      if (known && !known.ok) {
        failures.push(`public/models/${name} declares ${ext}. ${known.how}`);
      }
    }
  }
  return failures;
}

// --------------------------------------------------------------------------
// Run
// --------------------------------------------------------------------------

function main() {
  const out = [];
  const say = (s = "") => out.push(s);
  let failed = false;
  const fail = (heading, lines) => {
    failed = true;
    say(`FAIL  ${heading}`);
    for (const l of lines) say(`      ${l}`);
    say();
  };

  say("origin audit -- every host named in the shipped output must be declared");
  say();

  // 1. Proof of correctness.
  const selfTest = proveClassifier();
  if (selfTest.length) {
    fail("the classifier does not classify its own fixtures correctly", selfTest);
    say("The guard is broken, so its verdict on the real build means nothing.");
    process.stdout.write(out.join("\n") + "\n");
    process.exit(1);
  }
  say("  ok    classifier proved on 8 fixtures (forbidden CDN, protocol-relative");
  say("        CDN, undeclared host, doc link, namespace, parse base) + 4 quiet");
  say("        cases + reason-required");

  // 2 & 3. The table and the sweep.
  const declarations = [...DECLARED, studioDeclaration()];
  const hosts = new Map();
  let fileCount = 0;
  const missingDirs = [];
  for (const rel of SHIPPED_DIRS) {
    const root = join(APP_ROOT, rel);
    if (!existsSync(root) || !statSync(root).isDirectory()) {
      missingDirs.push(rel);
      continue;
    }
    for (const f of walkShipped(root)) {
      fileCount += 1;
      mergeHosts(hosts, hostsIn(readFileSync(f, "utf8"), relative(APP_ROOT, f)));
    }
  }

  if (missingDirs.includes(".next/static")) {
    say();
    say("FAIL  no build output at .next/static -- run `next build` first.");
    say("      An audit of nothing is not a passing audit.");
    process.stdout.write(out.join("\n") + "\n");
    process.exit(1);
  }
  say(`  ok    scanned ${fileCount} shipped files in ${SHIPPED_DIRS.join(", ")}`);
  say();

  const forbidden = [];
  const undeclared = [];
  const declaredSeen = [];
  for (const [host, sites] of [...hosts].sort()) {
    const { verdict } = classify(host, declarations);
    if (verdict === VERDICT.FORBIDDEN) forbidden.push([host, sites]);
    else if (verdict === VERDICT.UNDECLARED) undeclared.push([host, sites]);
    else declaredSeen.push(host);
  }

  if (forbidden.length) {
    fail(
      `${forbidden.length} forbidden host(s) in the shipped output`,
      forbidden.flatMap(([h, sites]) => [
        `${h} -- ${FORBIDDEN.find(([f]) => f === h)[1]}`,
        ...sites.map((s) => `  ${s.file}: ${s.context}`),
        "",
      ]),
    );
    say("      These are not exemptable. The page claims it contacts nothing but");
    say("      its own origin; a CDN default in the bundle makes that claim");
    say("      conditional on an asset property, which is how it broke before.");
    say();
  }

  if (undeclared.length) {
    fail(
      `${undeclared.length} undeclared external host(s) in the shipped output`,
      undeclared.flatMap(([h, sites]) => [
        h,
        ...sites.map((s) => `  ${s.file}: ${s.context}`),
        "",
      ]),
    );
    say("      Either remove it, or add a Declaration in scripts/audit-origins.mjs");
    say("      saying what it is and why it is never fetched. Writing the reason");
    say("      is the review -- that is the point of the table.");
    say();
  }

  if (!forbidden.length && !undeclared.length) {
    say(`  ok    ${declaredSeen.length} external host(s), all declared:`);
    for (const h of declaredSeen) {
      const d = declarations.find((x) => x.host === h);
      say(`          ${h.padEnd(24)} ${d.kind}`);
    }
    say();
  }

  // 4. Asset checks.
  const decoder = checkVendoredDecoder();
  if (decoder.length) fail("vendored DRACO decoder", decoder);
  else say(`  ok    vendored DRACO decoder present and unmodified (3 files, 745 KiB)`);

  const pin = checkDecoderPin();
  if (pin.length) fail("decoder path pin", pin);
  else say("  ok    glTF loading pinned at /draco/, no drei useGLTF import");

  const models = checkModelExtensions();
  if (models.length) fail("model extensions", models);
  else say("  ok    no model declares an extension needing a decoder we do not serve");

  say();
  if (failed) {
    say("origin audit FAILED");
    say();
    say("Note what this check does and does not settle. It proves the set of");
    say("external hosts NAMED in the output is the set that has been reviewed.");
    say("It cannot prove none is contacted -- a host can be assembled at runtime");
    say("from fragments it will never see. `node scripts/network-audit.mjs` is");
    say("the execution-level proof; run that too.");
  } else {
    say("origin audit passed");
    say();
    say("This says every external host named in the output is one that has been");
    say("looked at and written down -- not that none is contacted. For that,");
    say("`node scripts/network-audit.mjs` loads the built page in a browser and");
    say("asserts zero cross-origin requests with the scene fully loaded.");
  }
  process.stdout.write(out.join("\n") + "\n");
  process.exit(failed ? 1 : 0);
}

main();
