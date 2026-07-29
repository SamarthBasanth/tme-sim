"use client";

import { useSyncExternalStore } from "react";

function getServerSnapshot(): boolean {
  return false;
}

/**
 * Coarse device-tier signal used to scale particle counts, geometry detail,
 * and DPR down for phones/tablets per the spec's responsive requirements.
 * Deliberately simple (a single width breakpoint) rather than full UA
 * sniffing — the goal is a performance budget, not device detection.
 *
 * Built on `useSyncExternalStore` (rather than `useState` + `useEffect`)
 * so reading `window.matchMedia` never trips the "setState synchronously
 * in an effect" hooks-lint rule this Next.js version enforces.
 */
export function useIsMobile(breakpoint = 820): boolean {
  const query = `(max-width: ${breakpoint}px)`;

  return useSyncExternalStore(
    (callback) => {
      const mql = window.matchMedia(query);
      mql.addEventListener("change", callback);
      return () => mql.removeEventListener("change", callback);
    },
    () => window.matchMedia(query).matches,
    getServerSnapshot,
  );
}
