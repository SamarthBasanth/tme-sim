"use client";

import { useSyncExternalStore } from "react";

const QUERY = "(prefers-reduced-motion: reduce)";

function subscribe(callback: () => void): () => void {
  const mql = window.matchMedia(QUERY);
  mql.addEventListener("change", callback);
  return () => mql.removeEventListener("change", callback);
}

function getSnapshot(): boolean {
  return window.matchMedia(QUERY).matches;
}

function getServerSnapshot(): boolean {
  return false;
}

/**
 * Tracks `(prefers-reduced-motion: reduce)` via `useSyncExternalStore` —
 * the React-blessed way to subscribe to a browser API that can change
 * outside of React, without the "setState synchronously in an effect"
 * anti-pattern flagged by the stricter hooks lint rules that ship with
 * this Next.js version. Consumers use this to swap the pinned
 * scroll-scrubbed cinematic experience for a simplified, non-parallax
 * static layout (see ScrollExperience).
 */
export function usePrefersReducedMotion(): boolean {
  return useSyncExternalStore(subscribe, getSnapshot, getServerSnapshot);
}
