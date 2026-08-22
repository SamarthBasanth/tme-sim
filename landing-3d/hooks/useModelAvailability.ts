"use client";

import { useEffect, useState } from "react";

export const TUMOR_MODEL_URL = "/models/tumor.glb";

/**
 * `null` while the probe is still in flight -- not `false`.
 *
 * The distinction is load-bearing, and it was measured. This used to return a
 * plain boolean seeded `false`, so the very first render of every visit
 * asserted "there is no model" and `TumorModel` mounted the full procedural
 * placeholder: ~12,400 noise-sampled points across three clouds, five
 * `TubeGeometry` vessels, and four materials that were compiled and uploaded
 * to the GPU -- all of it thrown away a moment later, on every load, because
 * `public/models/tumor.glb` has existed since this page shipped.
 *
 * That cost was not just wasted: it was BLOCKING. Traced on the built page
 * (1512x950, `next start`, headless Chromium), the sequence was
 *
 *     48 ms   first paint, loading overlay up
 *   1740 ms   HEAD /models/tumor.glb   <- the probe, held behind a ~1.6 s task
 *   1825 ms   GET  /models/tumor.glb
 *
 * The effect that issues this probe cannot run until React commits, and the
 * commit it was queued behind was busy sampling and compiling a placeholder
 * for a file that was already there. So the 660 KB measured asset -- the one
 * thing on the page that actually has to arrive -- did not begin downloading
 * for nearly two seconds.
 *
 * `null` lets `TumorModel` render nothing until the answer is known. The
 * loading overlay is up over that window anyway, so nothing is lost visually,
 * and the placeholder is still mounted the moment the probe comes back
 * `false` -- which is the case it was written for.
 */
export type ModelAvailability = boolean | null;

/**
 * Checks (once, via a lightweight HEAD request) whether a real `tumor.glb`
 * has been dropped into `public/models/`. `SceneCanvas` uses this to decide,
 * with zero code changes required, whether to render the shipped procedural
 * placeholder or load the real asset. This is what makes the placeholder
 * "easy to replace": drop the file in, reload, done. It is also what tells
 * `Loader` whether a model is coming, so the overlay can wait for one instead
 * of guessing (see Loader.tsx).
 */
/**
 * One probe per URL per page load, shared across mounts.
 *
 * `SceneCanvas` is mounted more than once on a normal visit under
 * prefers-reduced-motion: `usePrefersReducedMotion` is a `useSyncExternalStore`
 * whose server snapshot is `false`, so the cinematic experience renders once
 * during hydration and is then replaced by the fallback. Traced on the built
 * page, that produced two identical HEAD requests 0 ms apart. Caching the
 * promise makes the second mount read the first mount's answer instead of
 * asking again.
 */
const probes = new Map<string, Promise<boolean>>();

function probe(url: string): Promise<boolean> {
  let inFlight = probes.get(url);
  if (!inFlight) {
    inFlight = fetch(url, { method: "HEAD" }).then(
      (res) => res.ok,
      () => false,
    );
    probes.set(url, inFlight);
  }
  return inFlight;
}

export function useModelAvailability(url: string = TUMOR_MODEL_URL): ModelAvailability {
  const [available, setAvailable] = useState<ModelAvailability>(null);

  useEffect(() => {
    let cancelled = false;
    probe(url).then((ok) => {
      if (!cancelled) setAvailable(ok);
    });
    return () => {
      cancelled = true;
    };
  }, [url]);

  return available;
}
