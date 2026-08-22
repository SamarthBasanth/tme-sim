"use client";

import { Component, type ErrorInfo, type ReactNode } from "react";

interface SceneErrorBoundaryProps {
  children: ReactNode;
  /** Rendered in place of `children` once they have thrown. */
  fallback?: ReactNode;
  /** Prefix for the single console warning emitted on first failure. */
  label: string;
  /**
   * Called once, on the first failure. Exists because something OUTSIDE the
   * failed subtree has to know: when the whole `<Canvas>` fails to get a WebGL
   * context, no model will ever be requested, and the loading overlay -- which
   * is a sibling of the canvas, not a child -- would otherwise wait for a
   * download that is never coming and sit opaque over the entire page forever.
   * Measured: blocking `getContext("webgl*")` produced exactly that, a black
   * screen with a spinner over a page whose copy was intact underneath.
   */
  onFailure?: () => void;
}

interface SceneErrorBoundaryState {
  failed: boolean;
}

/**
 * A React error boundary for the 3D scene, and the reason it exists is a
 * MEASURED whole-page failure rather than defensive habit.
 *
 * `useTumorGLTF` (see tumorGLTF.ts) loads `/models/tumor.glb` through
 * `useLoader`, which is Suspense-based: a rejected load is re-thrown during
 * render. R3F's `<Canvas>` wraps its own children in an internal error
 * boundary, but that boundary does not absorb the error -- it lifts it into
 * the host component and re-throws it OUTSIDE the canvas
 * (`if (error) throw error` in CanvasImpl). With no boundary above that point
 * the error reached the Next.js root and unmounted the entire page.
 *
 * Measured against the built page (127.0.0.1, `next start`), intercepting
 * `/models/tumor.glb`:
 *
 *   GET returns a truncated body -> "This page couldn't load", <main> empty,
 *                                   0 headings, no copy, no CTA.
 *   GET returns HTTP 500          -> identical.
 *
 * That is the worst possible failure for this particular page: a stale CDN
 * object or a partial deploy of one 660 KB binary took away every sentence,
 * every published limit and both calls to action, none of which need WebGL to
 * be read. The tissue render is the least load-bearing thing on the page and
 * it was the single point of failure for all of it.
 *
 * WHY THE FALLBACK FOR A FAILED MODEL IS NOTHING, NOT THE PROCEDURAL TUMOUR.
 * `ProceduralTumor` is a synthetic noise-sampled blob. Substituting it for a
 * failed load would put a generated shape under a hero that states the section
 * is 40,887 measured cells, next to a caption naming the panel, the extent and
 * the cell-type breakdown of a specific real specimen. Silently drawing
 * something else there is precisely the claim-surface drift this repository
 * keeps having, so a broken measured asset renders NO tissue: the scene keeps
 * its lighting and particles, the page keeps its argument, and nothing on
 * screen asserts to be measured data that is not.
 */
export default class SceneErrorBoundary extends Component<
  SceneErrorBoundaryProps,
  SceneErrorBoundaryState
> {
  state: SceneErrorBoundaryState = { failed: false };

  static getDerivedStateFromError(): SceneErrorBoundaryState {
    return { failed: true };
  }

  componentDidCatch(error: Error, info: ErrorInfo) {
    // One warning, not a throw: the page is expected to keep working without
    // this subtree, and a silent swallow would make a broken asset invisible
    // to whoever deployed it.
    console.warn(`[codon] ${this.props.label} failed; continuing without it.`, error, info);
    this.props.onFailure?.();
  }

  render() {
    if (this.state.failed) return this.props.fallback ?? null;
    return this.props.children;
  }
}
