"use client";

import { Suspense, useEffect, type ReactNode } from "react";
import { useGLTF } from "@react-three/drei";
import ProceduralTumor from "./ProceduralTumor";
import GLTFTumor from "./GLTFTumor";
import { TUMOR_MODEL_URL, useModelAvailability } from "@/hooks/useModelAvailability";

interface TumorModelProps {
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
 * Entry point for the tumor mesh. No `/models/tumor.glb` has been supplied
 * yet, so by default this renders `ProceduralTumor` — a distorted-sphere /
 * metaball-style placeholder built to be visually convincing on its own.
 *
 * The moment a real GLB is placed at `public/models/tumor.glb`, this
 * component swaps to `GLTFTumor` automatically (via a HEAD-request
 * availability check) with no code changes required. `useGLTF.preload` is
 * called eagerly so the swap has no loading pop-in.
 */
export default function TumorModel({
  lowDetail = false,
  reducedMotion = false,
  coreChildren,
  edemaChildren,
  marginChildren,
}: TumorModelProps) {
  const modelAvailable = useModelAvailability();

  useEffect(() => {
    if (modelAvailable) useGLTF.preload(TUMOR_MODEL_URL);
  }, [modelAvailable]);

  const procedural = (
    <ProceduralTumor
      lowDetail={lowDetail}
      reducedMotion={reducedMotion}
      coreChildren={coreChildren}
      edemaChildren={edemaChildren}
      marginChildren={marginChildren}
    />
  );

  if (!modelAvailable) {
    return procedural;
  }

  return (
    <Suspense fallback={procedural}>
      <GLTFTumor url={TUMOR_MODEL_URL} reducedMotion={reducedMotion} />
    </Suspense>
  );
}
