"use client";

import { useEffect, useState } from "react";

export const TUMOR_MODEL_URL = "/models/tumor.glb";

/**
 * Checks (once, via a lightweight HEAD request) whether a real `tumor.glb`
 * has been dropped into `public/models/`. TumorModel uses this to decide,
 * with zero code changes required, whether to render the shipped procedural
 * placeholder or load the real asset via `useGLTF`. This is what makes the
 * placeholder "easy to replace": drop the file in, reload, done.
 */
export function useModelAvailability(url: string = TUMOR_MODEL_URL): boolean {
  const [available, setAvailable] = useState(false);

  useEffect(() => {
    let cancelled = false;
    fetch(url, { method: "HEAD" })
      .then((res) => {
        if (!cancelled) setAvailable(res.ok);
      })
      .catch(() => {
        if (!cancelled) setAvailable(false);
      });
    return () => {
      cancelled = true;
    };
  }, [url]);

  return available;
}
