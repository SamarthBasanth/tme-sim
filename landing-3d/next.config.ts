import type { NextConfig } from "next";

/**
 * RESPONSE HEADERS.
 *
 * The page's stated property is that it contacts nothing but its own origin --
 * `scripts/audit-origins.mjs` proves no foreign host is NAMED in the build,
 * and `scripts/network-audit.mjs` drives the built page in a real browser and
 * proves none is CONTACTED. Both are audits of what we shipped. Neither stops
 * a browser from loading something else if a string ever gets injected into
 * the page. The Content-Security-Policy below is the enforcement half: the
 * same statement, expressed to the browser, applying to whatever the page
 * turns out to contain rather than to what we inspected.
 *
 * Every directive here was verified by loading the started production build in
 * headless Chromium and asserting zero CSP violations with the WebGL scene
 * live and the page scrolled end to end -- not by reasoning about it. The
 * non-obvious allowances are annotated below; each one exists because removing
 * it broke something observable.
 */

/**
 * WHY 'unsafe-inline' IN script-src, AND WHY IT IS NOT REPLACED BY A NONCE.
 *
 * The App Router streams its RSC payload as inline <script> blocks
 * (`self.__next_f.push(...)`). The production HTML for `/` carries five of
 * them. Without `'unsafe-inline'` the browser blocks all five, React never
 * hydrates, and the page renders as static server HTML with no scroll
 * experience and no canvas.
 *
 * Next's documented alternative is a per-request nonce issued from `proxy.ts`.
 * Its own docs state the constraint plainly: "you must use dynamic rendering
 * to add nonces". This route is currently fully prerendered (`○ (Static)` in
 * the build output, `x-nextjs-prerender: 1` on the response), which is what
 * lets a CDN serve it with no function invocation. Adopting nonces would trade
 * that for a server render on every view. That is a real cost, and the thing
 * being bought is small here: this page has no user input, no query
 * parameters, no forms, no third-party script, and no `dangerouslySetInnerHTML`
 * -- there is no injection vector for a nonce to close. `object-src 'none'`,
 * `base-uri 'self'` and the absence of any remote host in script-src remove
 * the escalations that usually make `'unsafe-inline'` dangerous.
 *
 * If this page ever grows a form, a search param it echoes, or any user- or
 * CMS-supplied string, that calculus flips and the nonce becomes worth the
 * dynamic render. This comment is the marker for that moment.
 *
 * `'unsafe-eval'` is deliberately NOT present. React uses `eval` only in
 * development; the production bundle does not, and the audit run below
 * confirms nothing reaches for it.
 */
const SCRIPT_SRC = "'self' 'unsafe-inline' 'wasm-unsafe-eval'";

/**
 * `'wasm-unsafe-eval'` is for the vendored DRACO decoder in `public/draco/`.
 * `WebAssembly.instantiate` counts as script evaluation under CSP, so
 * `draco_decoder.wasm` cannot start without it.
 *
 * Honesty about what has and has not been executed: `public/models/tumor.glb`
 * declares no `KHR_draco_mesh_compression`, so `DRACOLoader` never actually
 * fetches or instantiates the decoder today, and this allowance is therefore
 * NOT exercised by the verification run. It is here so that the day someone
 * re-exports the mesh with Draco compression -- an ordinary size optimisation,
 * and the exact scenario `components/canvas/tumorGLTF.ts` and
 * `scripts/audit-origins.mjs` were written against -- the scene does not go
 * blank because of a header nobody thought about. `worker-src blob:` is the
 * other half of that: DRACOLoader decodes in a Worker built from a Blob URL.
 */

/**
 * `'unsafe-inline'` in style-src covers the `style` ATTRIBUTES React writes for
 * the values that are genuinely computed rather than authored -- scroll-driven
 * heights, progress-bar widths, the swatch colours in the specimen caption.
 * The production HTML carries 22 of them before hydration and more after.
 * There is no `style-src-attr`-only formulation that keeps those and drops
 * inline <style> blocks without also needing hashes Next does not emit.
 */
const STYLE_SRC = "'self' 'unsafe-inline'";

/**
 * `data:` in img-src is for the SVG feTurbulence grain in `app/globals.css`,
 * which is a data: URI rather than a fetched asset precisely so that it is not
 * a request. `blob:` is for canvas readback (`toDataURL`/`toBlob`) paths in
 * three.js. Neither reaches the network.
 */
const CSP = [
  "default-src 'self'",
  `script-src ${SCRIPT_SRC}`,
  `style-src ${STYLE_SRC}`,
  "img-src 'self' data: blob:",
  "font-src 'self'",
  // The page makes no cross-origin request of any kind. `'self'` here is the
  // enforced form of the claim the two audit scripts merely observe.
  "connect-src 'self'",
  "media-src 'self'",
  "worker-src 'self' blob:",
  "manifest-src 'self'",
  "object-src 'none'",
  "base-uri 'self'",
  "form-action 'self'",
  "frame-ancestors 'none'",
  // Deliberately NOT set: `upgrade-insecure-requests`. The one outbound link
  // on the page is the Studio call to action, whose URL is deployment-supplied
  // (NEXT_PUBLIC_STUDIO_URL) and is `http://localhost:8000` in local and
  // default builds. Nothing else on the page is loaded over a scheme this
  // could fix -- every subresource is a same-origin relative path -- so the
  // directive would buy nothing and would silently break the local link.
].join("; ");

const SECURITY_HEADERS = [
  { key: "Content-Security-Policy", value: CSP },
  // Redundant with `frame-ancestors 'none'` for any browser released this
  // decade, and kept anyway: it is what old scanners and old browsers read.
  { key: "X-Frame-Options", value: "DENY" },
  { key: "X-Content-Type-Options", value: "nosniff" },
  { key: "Referrer-Policy", value: "strict-origin-when-cross-origin" },
  // No feature on this page needs any of these, and a page about tissue data
  // should be visibly unable to ask for them.
  {
    key: "Permissions-Policy",
    value: "camera=(), microphone=(), geolocation=(), browsing-topics=(), interest-cohort=()",
  },
  // Two years, WITHOUT `includeSubDomains` and WITHOUT `preload`. The Studio
  // is expected to live on a sibling subdomain of whatever domain this ends up
  // on; `includeSubDomains` from the marketing page would impose HTTPS on
  // hosts this deployment does not control and cannot verify, and `preload` is
  // effectively irreversible. Add them from the apex once every subdomain is
  // known to be HTTPS-only -- not before.
  { key: "Strict-Transport-Security", value: "max-age=63072000" },
];

/**
 * STATIC EXPORT, AND WHO SETS THE HEADERS THEN.
 *
 * `VIVUM_STATIC_EXPORT=1` builds this page to plain files in `out/` so it can
 * be served by the Python Studio at `/landing` -- a deploy path that needs a
 * git push and nothing else, which is the only path available when the hosting
 * dashboard is not reachable. Vercel builds are unaffected: without the flag
 * this file is byte-for-byte the configuration that was verified.
 *
 * `headers()` is NOT honoured by `output: "export"` -- there is no server left
 * to send them. Next errors on the combination rather than dropping them
 * silently, which is why the two are mutually exclusive below. The headers do
 * not disappear: `app/server/landing.py` sends the same list, and a test
 * asserts the two agree rather than trusting that they do.
 */
const STATIC_EXPORT = process.env.VIVUM_STATIC_EXPORT === "1";

const nextConfig: NextConfig = {
  /**
   * `X-Powered-By: Next.js` on every response is free reconnaissance and buys
   * nothing.
   */
  poweredByHeader: false,

  ...(STATIC_EXPORT
    ? {
        output: "export" as const,
        // The Studio owns `/`. Served under a prefix, every emitted asset URL
        // has to carry it or the page loads and the scene never does; Next
        // rewrites `/_next/*` from `basePath` alone, so this one line is the
        // whole adaptation. It applies ONLY to the static build -- a Vercel
        // deploy still owns its own root.
        basePath: "/landing",
      }
    : {
        async headers() {
          return [
            {
              // Everything, including `/models/*.glb`, `/draco/*` and the
              // static chunks. Headers are matched before the filesystem, so
              // `public/` assets are covered too.
              source: "/:path*",
              headers: SECURITY_HEADERS,
            },
          ];
        },
      }),
};

export { SECURITY_HEADERS };

export default nextConfig;
