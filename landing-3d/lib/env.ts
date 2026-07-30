/**
 * Client-readable environment configuration.
 *
 * `NEXT_PUBLIC_*` vars are inlined into the browser bundle at build time
 * (see Next.js docs on bundling env vars for the browser), so each
 * deployment of this landing page must set `NEXT_PUBLIC_STUDIO_URL` to
 * wherever the actual Studio app is hosted. Falls back to the local Studio
 * dev server so `npm run dev` works out of the box without any setup.
 */
export const STUDIO_URL = process.env.NEXT_PUBLIC_STUDIO_URL ?? "http://localhost:8000";
