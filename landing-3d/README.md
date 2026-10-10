# landing-3d — the Vivum marketing site

The public page. Next.js, deployed to Vercel. It is a claim surface, so the
rules below are not style preferences.

## Running it

```bash
npm run dev          # http://localhost:3000
```

A **production** build refuses without `NEXT_PUBLIC_STUDIO_URL`:

```bash
NEXT_PUBLIC_STUDIO_URL=https://vivum.onrender.com npm run build
```

That variable is the product's origin, and it is inlined into the browser
bundle at build time — changing it later needs a redeploy, not a settings
change. Unset, the page's calls to action would point at `localhost:8000`, i.e.
at the visitor's own machine, and every check would still pass; `lib/envGuard.ts`
exists to make that a failed build instead. See `DEPLOY_VERCEL.md`.

**Do not judge a build by grepping its log for a success string.** Check the
exit code. A previous session grepped for "Compiled successfully" while the
build exited 1 throughout.

## The four checks

```bash
npm run build                   # also runs the origin audit
node scripts/network-audit.mjs  # loads the built page in a browser
node scripts/viewport-audit.mjs # does every beat fit the viewport that clips it
grep -r localhost .next/static/
```

* **`scripts/audit-origins.mjs`** reads the build output as text and asserts
  that every external host NAMED in it has been reviewed and written down.
* **`scripts/network-audit.mjs`** loads the built page in Chromium at four
  viewport/motion combinations, scrolls the whole timeline, and asserts zero
  cross-origin requests with the 3D scene fully loaded. It refuses to report a
  pass unless the conditional work it exists to catch actually happened. It
  needs Playwright, which is deliberately not a dependency here — see the file.
* **`scripts/viewport-audit.mjs`** serves the built page at eleven viewports in
  both motion settings, samples the whole scroll timeline, and asserts that no
  beat is cut off by `.scroll-pin`'s `overflow: hidden`, that neither floating
  panel intersects the copy beside it, that no beat travels under the fixed
  wordmark, and that nothing scrolls horizontally.
  It exists because every responsive rule in `app/globals.css` is justified by
  a number someone measured once by hand and could not measure again, and
  because 320x568 shipped with five of six beats clipped — losing, in two of
  them, the last rows of a ledger. Also needs Playwright.
* The `localhost` grep must return only core-js's URL polyfill, which is
  third-party and unavoidable. Anything else is the defect `lib/env.ts`
  documents at length.

## Where the claims live

`lib/copy.ts` holds every rendered sentence and, above them, the guardrails
that say what may and may not be claimed — what was retracted, what measured
it, and what must never drift back. **Read it before editing any copy.** Every
number on the page has to trace to a file in `../results/`; a figure the
product computed once is not a citation.

Two standing prohibitions, both absolute: no cross-compound ordering in any
form, and no absolute kill fraction. Tissue registers by filesystem path — the
word "upload" does not appear.

## Layout

The page has two layouts and they carry identical copy:

* **pinned** — `components/scroll/ScrollExperience.tsx`, one ScrollTrigger
  driving six beats over a 3D scene.
* **flow** — `components/scroll/ReducedMotionExperience.tsx`, normal document
  order. Reached by `prefers-reduced-motion: reduce` OR by a viewport shorter
  than `MIN_PINNED_VIEWPORT_HEIGHT` (`hooks/useShortViewport.ts`), because the
  pinned layout clips at `overflow: hidden` and a clipped ledger row is a
  disclosure that went missing.

`app/globals.css` carries responsive tiers keyed on BOTH width and height, each
with the measurement that produced it in a comment above it. **If you change
copy length, run `scripts/viewport-audit.mjs`** — the tiers were tuned to real
numbers, several beats land within 10px of the edge that cuts them, and a
clipped ledger row is the one failure this page cannot have.

## Fonts

Vendored, not fetched: the four `.woff2` files are committed in `app/fonts/`
with their OFL licences and loaded through `next/font/local`. Newsreader
(display), Atkinson Hyperlegible Next (body, roman + italic) and Atkinson
Hyperlegible Mono (instruments). Chosen in `app/layout.tsx`, which also records
why Geist and Archivo were rejected and why the fallback metrics are supplied
by hand in `app/globals.css`.

This is what lets the network audit assert zero third-party requests at
runtime, and it is also why the build needs no egress and why the outlines are
byte-pinned — the tightest beat on the page has 8.9px of vertical headroom, so a
family re-cut upstream could change a line count with nothing in the diff.
