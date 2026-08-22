"use client";

import { useSyncExternalStore } from "react";

/**
 * Can this browser give us a WebGL context at all?
 *
 * WHY THIS IS A PROBE AND NOT AN ERROR BOUNDARY. `SceneErrorBoundary` catches
 * a failed model load, because `useLoader` re-throws during render. It cannot
 * catch this one: R3F builds the renderer inside an `async function run()`
 * called from a layout effect (see CanvasImpl), so a `WebGLRenderer` that
 * cannot get a context rejects a promise rather than throwing through React.
 * Measured with `getContext("webgl*")` stubbed to return null, that surfaced as
 * an uncaught `pageerror`, the Canvas's own error state was never set, and the
 * boundary above it was never invoked. An error boundary here would have been a
 * check blind to its own subject.
 *
 * What made this worth detecting rather than tolerating: the loading overlay is
 * `position: fixed; inset: 0` on an opaque background and waits for a model.
 * With no context nothing is ever requested, so it waited forever -- a black
 * page over six beats of prose, both CTAs and every published limit, all of
 * which are plain DOM and were rendering correctly underneath it.
 *
 * The detection canvas is discarded immediately via `WEBGL_lose_context`, so
 * this does not hold a second context open on machines with a small budget for
 * them. The result is cached at module scope: it cannot change during a page
 * load, and `useSyncExternalStore` requires a stable snapshot.
 */
let cached: boolean | null = null;

function detect(): boolean {
  if (cached !== null) return cached;
  try {
    const canvas = document.createElement("canvas");
    const gl =
      (canvas.getContext("webgl2") as WebGLRenderingContext | null) ??
      (canvas.getContext("webgl") as WebGLRenderingContext | null);
    cached = gl !== null;
    (gl?.getExtension("WEBGL_lose_context") as { loseContext(): void } | null)?.loseContext();
  } catch {
    cached = false;
  }
  return cached;
}

/** No-op: WebGL support cannot change within a document. */
function subscribe(): () => void {
  return () => {};
}

/**
 * `true` on the server and during hydration, so the markup React renders on
 * the client matches what it sent — the same shape `usePrefersReducedMotion`
 * uses. The real answer is read on the commit after hydration; on the
 * overwhelmingly common path (WebGL present) it agrees and nothing re-renders.
 */
function getServerSnapshot(): boolean {
  return true;
}

export function useWebGLSupport(): boolean {
  return useSyncExternalStore(subscribe, detect, getServerSnapshot);
}
