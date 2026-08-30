/**
 * The build-time refusals for `NEXT_PUBLIC_STUDIO_URL`.
 *
 * WHY THIS IS A SEPARATE FILE FROM `lib/env.ts`, AND IT IS NOT TIDINESS.
 *
 * `lib/env.ts` is imported by client components, so everything in it is
 * compiled into the browser bundle. Both checks below used to live there, and
 * the loopback one carries a regular expression with the literal string
 * `localhost` in it. Measured on the shipped output of a passing build:
 *
 *   .next/static/chunks/1jndf5v81c0gx.js
 *     let bu="https://studio.example.org";
 *     if(bu&&/^(https?:\/\/)?(localhost|127\.0\.0\.1|0\.0\.0\.0|\[::1\])(:|\/|$)/i
 *        .test(bu.trim())&&1)throw ...
 *
 * Two things are wrong with that and neither is cosmetic.
 *
 * (1) IT SHIPS AND RE-RUNS. `process.env.NODE_ENV === "production"` folds to the
 *     `&&1` you can see, so this build-time assertion was also being evaluated
 *     in every visitor's browser, on a value that is a compile-time constant and
 *     therefore cannot have changed. It can never fire and it is never useful.
 *
 * (2) IT BREAKS THE CHECK THAT WATCHES IT. The invariant this project actually
 *     runs against the output is "no `localhost` in `.next/static/`" — the
 *     header comment on `STUDIO_URL` in lib/env.ts records the time a fallback
 *     defeated dead-code elimination and put `http://localhost:8000` into the
 *     bundle. A guard that permanently plants the word `localhost` in the same
 *     directory makes that grep return a hit on every clean build, and a check
 *     that is always red is a check people stop reading. This repository has
 *     the failure written down twice already (see DEPLOY_VERCEL.md §2 on the
 *     origin audit's setup failure being indistinguishable from a finding).
 *
 * So the ASSERTION moves to where it was always meant to run — once, at build
 * time, on the server — and the VALUE stays in lib/env.ts where the client can
 * read it. `app/page.tsx` is a server component and calls this in its body, so
 * it executes while the page is prerendered: exactly when it executed before,
 * with the same message and the same non-zero exit, and nothing of it reaches
 * `.next/static/`.
 *
 * DO NOT import this from a client component, and do not make it a bare
 * side-effect import. A call in a server component's body is what makes it
 * unremovable; an unused `import "./envGuard"` is something a bundler is
 * entitled to drop.
 */

/** Hosts that only ever mean "this machine". */
const LOOPBACK = /^(https?:\/\/)?(localhost|127\.0\.0\.1|0\.0\.0\.0|\[::1\])(:|\/|$)/i;

/**
 * Refuses a production build whose Studio URL is missing or points at the
 * machine loading the page. A no-op in development, where localhost is right.
 */
export function assertStudioUrl(): void {
  const raw = process.env.NEXT_PUBLIC_STUDIO_URL;

  if (process.env.NODE_ENV !== "production") return;

  if (raw && LOOPBACK.test(raw.trim())) {
    // A SET-BUT-LOOPBACK VALUE IS WORSE THAN AN UNSET ONE, and the unset check
    // below cannot see it. `.env.local.example` ships
    // `NEXT_PUBLIC_STUDIO_URL=http://localhost:8000` because that is correct for
    // local development -- so importing that file into a hosting dashboard, which
    // is the obvious thing to do with it, would satisfy the unset check and ship
    // the precise defect the unset check exists to prevent: a public page whose
    // primary call to action points at the visitor's own machine.
    throw new Error(
      `NEXT_PUBLIC_STUDIO_URL is set to a loopback address (${raw}) in a ` +
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

  if (!raw) {
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
}
