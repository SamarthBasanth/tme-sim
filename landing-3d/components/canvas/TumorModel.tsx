"use client";

import { Suspense, type ReactNode } from "react";
import ProceduralTumor from "./ProceduralTumor";
import GLTFTumor from "./GLTFTumor";
import SceneErrorBoundary from "./SceneErrorBoundary";
import { TUMOR_MODEL_URL } from "@/hooks/useModelAvailability";
import type { ModelAvailability } from "@/hooks/useModelAvailability";

interface TumorModelProps {
  /**
   * Result of the availability probe, owned by SceneCanvas so that the Loader
   * outside the canvas can read the same answer. `null` = not yet known.
   */
  available: ModelAvailability;
  lowDetail?: boolean;
  reducedMotion?: boolean;
  /**
   * Forwarded to ProceduralTumor's per-region group children (see that
   * component). Not passed to GLTFTumor: a real scan-derived model has no
   * equivalent core/edema/margin group split, so these are specific to the
   * placeholder path and will need their own approach once a real GLB
   * exists.
   */
  coreChildren?: ReactNode;
  edemaChildren?: ReactNode;
  marginChildren?: ReactNode;
}

/**
 * Entry point for the tumor mesh.
 *
 * `public/models/tumor.glb` exists, so the path that actually ships is
 * `GLTFTumor` and the real Xenium section. `ProceduralTumor` remains the
 * documented drop-out: remove the GLB and the page still renders a tumour.
 *
 * THREE STATES, NOT TWO, AND THAT IS THE POINT.
 *
 *   available === null   the probe has not answered. Render NOTHING. The
 *                        loading overlay covers this window, and mounting a
 *                        placeholder here used to delay the real asset's
 *                        download by ~1.75 s (see useModelAvailability).
 *   available === false  no GLB on the server. The placeholder is the whole
 *                        point; render it.
 *   available === true   load the measured asset.
 *
 * The Suspense fallback is `null` rather than the placeholder for the same
 * reason it is null above, plus an honesty one: the window it covers is the
 * window the loading overlay owns, and flashing a synthetic blob underneath a
 * hero that says the tissue is measured is not an improvement on an empty
 * frame.
 *
 * `SceneErrorBoundary` is what stops a bad 660 KB binary from taking the page
 * down with it; see that file for the measured failure it was written against.
 */
export default function TumorModel({
  available,
  lowDetail = false,
  reducedMotion = false,
  coreChildren,
  edemaChildren,
  marginChildren,
}: TumorModelProps) {
  if (available === null) return null;

  if (!available) {
    return (
      <ProceduralTumor
        lowDetail={lowDetail}
        reducedMotion={reducedMotion}
        coreChildren={coreChildren}
        edemaChildren={edemaChildren}
        marginChildren={marginChildren}
      />
    );
  }

  return (
    <SceneErrorBoundary label={`tissue model ${TUMOR_MODEL_URL}`} fallback={null}>
      <Suspense fallback={null}>
        <GLTFTumor url={TUMOR_MODEL_URL} reducedMotion={reducedMotion} />
      </Suspense>
    </SceneErrorBoundary>
  );
}
