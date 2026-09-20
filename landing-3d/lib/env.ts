/**
 * Client-readable environment configuration.
 *
 * `NEXT_PUBLIC_*` vars are inlined into the browser bundle AT BUILD TIME, so
 * the value baked here is permanent for that deployment — setting it after a
 * deploy does nothing until the next build. That property is what makes the
 * old silent fallback dangerous rather than convenient.
 *
 * REFUSES IN A PRODUCTION BUILD, 2026-08-22. `STUDIO_URL` used to be
 *
 *     process.env.NEXT_PUBLIC_STUDIO_URL ?? "http://localhost:8000"
 *
 * with no check. A deploy that forgot the variable was GREEN — the build
 * passed, the page rendered, both the origin audit and the network audit
 * passed — and shipped a public marketing site whose "Open the workspace"
 * calls to action pointed at `localhost:8000`, i.e. at the visitor's own
 * machine.
 * Nothing anywhere would have reported it, and the fix would have required a
 * rebuild rather than a settings change.
 *
 * So: dev still works with no setup, and a production build REFUSES rather
 * than defaulting. This is the same rule the rest of this project applies to
 * any configuration fact that cannot be inferred and changes what a user sees
 * — an empty `TMESIM_STUDIO_DATA_ROOTS` fails closed rather than being
 * helpfully replaced, and the reach map refuses when told nothing about
 * convection instead of guessing. A default that is right on a laptop and
 * wrong in production is not a default, it is a trap with a plausible excuse.
 *
 * THE REFUSALS THEMSELVES ARE NOT IN THIS FILE, AND THAT SEPARATION IS
 * DELIBERATE. This module is imported by client components, so anything in it
 * is compiled into the browser bundle — and the loopback check needs a regex
 * containing the literal word `localhost`, which then sits permanently in
 * `.next/static/` and makes the grep that watches for exactly that string
 * return a hit on every clean build. They live in `lib/envGuard.ts`, called
 * from the server component in `app/page.tsx`; that file carries the full
 * argument.
 */

/**
 * The product's ORIGIN. Localhost only ever in development.
 *
 * Still named `STUDIO_URL` after the page it used to point at. As of
 * docs/WORKSPACE_CONTRACT.md §1 the product serves the workspace home at `/`
 * and the Studio at `/studio`, so this value is the origin and the page's call
 * to action lands on the workspace. The variable keeps its name because
 * `NEXT_PUBLIC_STUDIO_URL` is set in a hosting dashboard, outside this repo:
 * renaming it would silently unset it on every existing deployment, and an
 * unset value is the exact failure the refusals below exist to catch.
 *
 * READS `process.env` DIRECTLY RATHER THAN VIA A LOCAL, AND THAT IS
 * LOAD-BEARING. Next inlines `process.env.NEXT_PUBLIC_STUDIO_URL` as a string
 * literal at build time, so in a production build this whole expression folds
 * to the real URL and the localhost fallback is dead-code-eliminated — it never
 * reaches the bundle. Verified on the shipped output: `grep -rc "localhost:8000"
 * .next/static/` returns nothing.
 *
 * Writing `RAW ?? "..."` against a `const RAW` broke that, because `RAW` was
 * also referenced inside the guards' throw messages and the minifier could no
 * longer prove the fallback unreachable. The literal survived into
 * `.next/static/chunks/*.js`, and `scripts/audit-origins.mjs` correctly failed
 * the build with "1 undeclared external host: localhost".
 *
 * Which is the joke worth leaving in the file: the guard written to keep
 * localhost out of production was the thing that put localhost into the
 * production bundle. The origin audit caught it; a grep for
 * "Compiled successfully" did not. Moving the guards out of this module (see
 * the note above) removes the coupling that made that possible at all.
 */
export const STUDIO_URL =
  process.env.NEXT_PUBLIC_STUDIO_URL ?? "http://localhost:8000";
