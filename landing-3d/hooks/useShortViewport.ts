"use client";

import { useSyncExternalStore } from "react";

/**
 * THE SHORTEST VIEWPORT THE PINNED EXPERIENCE IS ALLOWED TO RUN AT.
 *
 * `.scroll-pin` is `height: 100svh; overflow: hidden`, so a beat that is taller
 * than the viewport is not scrolled — it is CUT, and what gets cut off the
 * bottom of beats 02 and 05 is the last rows of their ledgers, i.e. the
 * registered limits, which are the one thing on this page that must never be
 * the part that goes missing.
 *
 * Every fit rule in globals.css buys height back by spending leading, gaps and
 * type size. That works down to a floor and then it stops: at some point the
 * copy is simply taller than the window and no amount of tightening is honest,
 * because the only thing left to spend is a disclosure.
 *
 * WHERE THE FLOOR IS. Measured against a production `next start` build on the
 * commit immediately BEFORE this gate existed — which is the last build in
 * which these viewports could be measured at all, since they now fall through —
 * 240 scroll positions per viewport, worst case over every opacity, content
 * edges against the pin:
 *
 *     360x700   beat 05 bottom 676.0 of 700   fits (24px spare)
 *     360x680   beat 05 bottom 676.1 of 680   fits (3.9px spare)
 *     360x660   beat 05 bottom 676.0 of 660   CLIPPED 16.0px, beat 02 by 1.7px
 *     360x640   beats 02, 05, 06 clipped by 21.8, 36.1 and 9.1px
 *     320x568   FIVE of six beats clipped; beat 02 by 180.4px
 *
 * So the pinned layout bottomed out at ~677px, and the threshold is 700: the
 * smallest round height measured to fit, with more headroom than a ledger
 * detail line (13.83px), so the page can gain a row without silently
 * re-breaking its own floor. 680 also fitted, on 3.9px, which is less than a
 * line of anything and would have made the next copy edit a regression nobody
 * would see.
 *
 * WHAT THE BOUNDARY MEASURES TODAY, re-measured after this pass's copy changes
 * with `scripts/viewport-audit.mjs`, which is now the way to take this number
 * again: at 360x700 the tightest beats are 05 at 677.8 (22.2px spare) and 06 at
 * 686.8 (13.2px spare). 06 is inside one line of the edge. THE NEXT ADDITION TO
 * BEAT 06 HAS TO BE MEASURED, not estimated.
 *
 * BELOW IT THE PAGE FALLS THROUGH TO THE FLOW LAYOUT, WHICH CANNOT CLIP. That
 * is the whole argument for doing this structurally rather than with a fourth
 * tier of CSS: `ReducedMotionExperience` puts every beat in normal document
 * order, so its height is unbounded and no ledger row can be cut off no matter
 * how short the window is. It carries the same copy, the same ledgers and the
 * same two instrument panels. What a short-viewport visitor loses is the
 * choreography; what they stop losing is the disclosures.
 *
 * A MEDIA QUERY, NOT `window.innerHeight`, AND THE TRADE IS RECORDED BECAUSE
 * IT IS NOT FREE. `innerHeight` is the quantity that actually matches what the
 * pin renders at (it tracks the small viewport, toolbars included), but on a
 * phone it GROWS as the browser chrome collapses during scroll — so gating on
 * it would tear the experience down and rebuild it mid-scroll. A media query
 * is stable, re-evaluates correctly when a desktop window is resized, and is
 * the same mechanism every other responsive decision on this page uses.
 *
 * WHAT THAT LEAVES OPEN, STATED PLAINLY: `@media (height)` sees the LARGE
 * viewport, so a phone reporting 740px renders this pin at roughly 640px and
 * still passes a 700px gate. That gap is the one globals.css already documents
 * where it takes every ≤820px-wide viewport up to 880px tall through its
 * tightest tier; closing it properly needs a measurement on a real device,
 * which has not been made. This gate fixes the case that IS measurable and
 * does not pretend to fix the one that is not.
 */
export const MIN_PINNED_VIEWPORT_HEIGHT = 700;

const QUERY = `(min-height: ${MIN_PINNED_VIEWPORT_HEIGHT}px)`;

function subscribe(callback: () => void): () => void {
  const mql = window.matchMedia(QUERY);
  mql.addEventListener("change", callback);
  return () => mql.removeEventListener("change", callback);
}

function getSnapshot(): boolean {
  return !window.matchMedia(QUERY).matches;
}

/**
 * Server snapshot: `false`, i.e. "tall enough".
 *
 * The same choice `usePrefersReducedMotion` makes, for the same reason — the
 * server has no viewport, and React requires a snapshot that is stable across
 * the hydration render. The client corrects it on the first commit, which is
 * the identical swap the reduced-motion path already performs.
 */
function getServerSnapshot(): boolean {
  return false;
}

/**
 * True when the viewport is too short for the pinned experience to hold a beat
 * without clipping it. See `MIN_PINNED_VIEWPORT_HEIGHT` for the measurement.
 */
export function useShortViewport(): boolean {
  return useSyncExternalStore(subscribe, getSnapshot, getServerSnapshot);
}
