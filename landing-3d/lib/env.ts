/**
 * Client-readable environment configuration.
 *
 * `NEXT_PUBLIC_*` vars are inlined into the browser bundle AT BUILD TIME, so
 * the value baked here is permanent for that deployment — setting it after a
 * deploy does nothing until the next build. That property is what makes the
 * old silent fallback dangerous rather than convenient.
 *
 * REFUSES IN A PRODUCTION BUILD, 2026-08-22. This used to be
 *
 *     process.env.NEXT_PUBLIC_STUDIO_URL ?? "http://localhost:8000"
 *
 * unconditionally. A deploy that forgot the variable was GREEN — the build
 * passed, the page rendered, both the origin audit and the network audit
 * passed — and shipped a public marketing site whose "Open the Studio" calls
 * to action pointed at `localhost:8000`, i.e. at the visitor's own machine.
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
 */
const RAW = process.env.NEXT_PUBLIC_STUDIO_URL;

/** Hosts that only ever mean "this machine". */
const LOOPBACK = /^(https?:\/\/)?(localhost|127\.0\.0\.1|0\.0\.0\.0|\[::1\])(:|\/|$)/i;

if (RAW && LOOPBACK.test(RAW.trim()) && process.env.NODE_ENV === "production") {
  // A SET-BUT-LOOPBACK VALUE IS WORSE THAN AN UNSET ONE, and the unset check
  // below cannot see it. `.env.local.example` ships
  // `NEXT_PUBLIC_STUDIO_URL=http://localhost:8000` because that is correct for
  // local development -- so importing that file into a hosting dashboard, which
  // is the obvious thing to do with it, would satisfy the unset check and ship
  // the precise defect the unset check exists to prevent: a public page whose
  // primary call to action points at the visitor's own machine.
  throw new Error(
    `NEXT_PUBLIC_STUDIO_URL is set to a loopback address (${RAW}) in a ` +
      "production build.\n\n" +
      "That address means 'this computer' to whoever loads the page, so the " +
      "'Open the Studio' links would point at the VISITOR'S machine, not at " +
      "the Studio. It is inlined into the browser bundle at build time.\n\n" +
      "This usually means .env.local.example was imported into a hosting " +
      "dashboard. That file is correct for `npm run dev` and wrong for a " +
      "deployment. Set the variable to the Studio's real public URL and " +
      "rebuild.\n\n" +
      "See landing-3d/DEPLOY_VERCEL.md.",
  );
}

if (!RAW && process.env.NODE_ENV === "production") {
  throw new Error(
    "NEXT_PUBLIC_STUDIO_URL is not set, and this is a production build.\n\n" +
      "It is inlined into the browser bundle at build time, so leaving it " +
      "unset would ship a page whose 'Open the Studio' links point at " +
      "http://localhost:8000 — the VISITOR'S machine, not the Studio. The " +
      "build would pass and both audits would pass; nothing would report it.\n\n" +
      "Set it to wherever the Studio is actually hosted, then rebuild:\n" +
      "  Vercel  Project Settings -> Environment Variables -> " +
      "NEXT_PUBLIC_STUDIO_URL\n" +
      "  Local   NEXT_PUBLIC_STUDIO_URL=https://studio.example.org npm run build\n\n" +
      "Changing it later requires a REDEPLOY, not just a settings change.\n" +
      "See landing-3d/DEPLOY_VERCEL.md.",
  );
}

/**
 * Where the Studio is hosted. Localhost only ever in development.
 *
 * READS `process.env` DIRECTLY RATHER THAN `RAW`, AND THAT IS LOad-BEARING.
 * Next inlines `process.env.NEXT_PUBLIC_STUDIO_URL` as a string literal at
 * build time, so in a production build this whole expression folds to the real
 * URL and the localhost fallback is dead-code-eliminated -- it never reaches
 * the bundle.
 *
 * Writing `RAW ?? "..."` instead broke that, because `RAW` is also referenced
 * inside the throw messages above and the minifier could no longer prove the
 * fallback unreachable. The literal survived into
 * `.next/static/chunks/*.js`, and `scripts/audit-origins.mjs` correctly failed
 * the build with "1 undeclared external host: localhost".
 *
 * Which is the joke worth leaving in the file: the guard written to keep
 * localhost out of production was the thing that put localhost into the
 * production bundle. The origin audit caught it; a grep for
 * "Compiled successfully" did not.
 */
export const STUDIO_URL =
  process.env.NEXT_PUBLIC_STUDIO_URL ?? "http://localhost:8000";
