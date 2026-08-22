"use client";

import { useEffect, useId, useState } from "react";
import { useProgress } from "@react-three/drei";
import type { ModelAvailability } from "@/hooks/useModelAvailability";

/**
 * Visible loading state for the 3D scene.
 *
 * WHAT WAS WRONG WITH THE PREVIOUS VERSION, TRACED RATHER THAN REASONED.
 *
 * It decided it was finished with `total > 0 ? !active && progress >= 100 :
 * true` -- i.e. "nothing has been queued with three's LoadingManager, so there
 * is nothing to wait for". At mount nothing HAS been queued yet, because the
 * availability probe and the glTF load both start after the first commit. So
 * the overlay was born in its dismissing state and stayed there. Sampled every
 * 20 ms against the built page (1512x950, `next start`, headless Chromium),
 * with the model body deliberately delayed 4 s:
 *
 *     48 ms  class="scene-loader scene-loader--done"  "Loading model — 0%"
 *   1740 ms  HEAD /models/tumor.glb
 *   1825 ms  GET  /models/tumor.glb
 *   1981 ms  overlay removed          <- while the model still had ~3.8 s to go
 *
 * It read 0% for the whole of its life, was already fading before the request
 * existed, and left. On a fast local load it happened to linger long enough to
 * look right, because the ~1.6 s startup task delayed its own dismissal timer;
 * that coincidence was the entire illusion. A visitor on a slow connection got
 * an empty black frame with no indication anything was coming.
 *
 * HOW IT DECIDES NOW. The scene knows whether a model is expected -- the same
 * availability probe that picks the mesh (see useModelAvailability) -- so it is
 * passed in rather than inferred from the LoadingManager's idle state:
 *
 *   expectModel === null   probe unresolved. Hold. We do not yet know whether
 *                          anything is coming.
 *   expectModel === false  no GLB on the server; the procedural placeholder
 *                          needs no download. Dismiss.
 *   expectModel === true   hold until the LoadingManager has actually seen the
 *                          item (`total > 0`, which drei's store only ever
 *                          moves upward) AND finished it. The gap between
 *                          "probe answered" and "request registered" can no
 *                          longer be mistaken for completion, because in that
 *                          gap `total` is still 0 and this reads false rather
 *                          than the true it used to short-circuit to.
 *
 * FAILURE DISMISSES TOO, from either direction. `useProgress().errors` is
 * non-empty once a load has failed, and `sceneFailed` covers the case where
 * there is no canvas to load into at all; in neither case will anything ever
 * reach 100%. SceneErrorBoundary has
 * already dropped the tissue subtree at that point and the page behind this
 * overlay is intact and readable, so holding a spinner over it forever would be
 * strictly worse than getting out of the way.
 *
 * `dismissed` is only ever set from inside the timeout callback (never
 * synchronously in the effect body) and only ever transitions false -> true.
 */
export default function Loader({
  expectModel,
  sceneFailed = false,
}: {
  expectModel: ModelAvailability;
  /**
   * The canvas itself failed to start (no WebGL context). Nothing will ever be
   * requested, so waiting for a download is waiting forever -- and this overlay
   * is opaque and covers the viewport, so "forever" means a black page over
   * copy that is perfectly readable. Measured: without this, blocking
   * `getContext("webgl*")` left the spinner up permanently.
   */
  sceneFailed?: boolean;
}) {
  const { progress, active, total, errors } = useProgress();
  const [dismissed, setDismissed] = useState(false);
  const labelId = useId();

  /**
   * The one genuine hang this can have: the availability probe never answers,
   * so `expectModel` stays `null` forever and the overlay covers a page that
   * is otherwise perfectly readable. A same-origin HEAD to the server that
   * just delivered this HTML does not take twelve seconds unless something is
   * badly wrong, and at that point getting out of the way beats a spinner over
   * six beats of prose.
   *
   * Scoped to the unresolved-probe case ONLY. Once `expectModel` is `true`
   * there is real progress information and the overlay should keep showing it
   * however long the download takes -- a ceiling there would recreate exactly
   * the bug this file documents, just later.
   */
  const [probeTimedOut, setProbeTimedOut] = useState(false);
  useEffect(() => {
    if (expectModel !== null) return;
    const timeout = setTimeout(() => setProbeTimedOut(true), 12_000);
    return () => clearTimeout(timeout);
  }, [expectModel]);

  const failed = errors.length > 0 || sceneFailed;
  const isDone =
    failed ||
    probeTimedOut ||
    expectModel === false ||
    (expectModel === true && total > 0 && !active && progress >= 100);

  useEffect(() => {
    if (!isDone) return;
    const timeout = setTimeout(() => setDismissed(true), 420);
    return () => clearTimeout(timeout);
  }, [isDone]);

  if (dismissed) return null;

  const pct = Math.min(100, progress);

  return (
    <div className={`scene-loader ${isDone ? "scene-loader--done" : ""}`}>
      <div className="scene-loader__mark" aria-hidden="true" />
      {/*
        The bar carries the ARIA, not the wrapper: a `progressbar` needs a
        name and a value, and both already exist on screen. `aria-labelledby`
        points at the visible label so the announced text and the rendered
        text cannot drift apart -- there is no second copy of the string to
        keep in sync.
      */}
      <div
        className="scene-loader__bar"
        role="progressbar"
        aria-labelledby={labelId}
        aria-valuemin={0}
        aria-valuemax={100}
        aria-valuenow={Math.round(pct)}
      >
        <div className="scene-loader__bar-fill" style={{ width: `${pct}%` }} />
      </div>
      <span className="scene-loader__label" id={labelId}>
        Loading model — {Math.round(pct)}%
      </span>
    </div>
  );
}
