"use client";

import { Suspense, useEffect } from "react";
import { useGLTF } from "@react-three/drei";
import ProceduralTumor from "./ProceduralTumor";
import GLTFTumor from "./GLTFTumor";
import { TUMOR_MODEL_URL, useModelAvailability } from "@/hooks/useModelAvailability";

interface TumorModelProps {
  lowDetail?: boolean;
  reducedMotion?: boolean;
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
export default function TumorModel({ lowDetail = false, reducedMotion = false }: TumorModelProps) {
  const modelAvailable = useModelAvailability();

  useEffect(() => {
    if (modelAvailable) useGLTF.preload(TUMOR_MODEL_URL);
  }, [modelAvailable]);

  if (!modelAvailable) {
    return <ProceduralTumor lowDetail={lowDetail} reducedMotion={reducedMotion} />;
  }

  return (
    <Suspense fallback={<ProceduralTumor lowDetail={lowDetail} reducedMotion={reducedMotion} />}>
      <GLTFTumor url={TUMOR_MODEL_URL} reducedMotion={reducedMotion} />
    </Suspense>
  );
}
