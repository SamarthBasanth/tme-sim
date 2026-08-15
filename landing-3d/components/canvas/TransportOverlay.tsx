"use client";

import { useRef, type ReactNode } from "react";
import { useFrame } from "@react-three/fiber";
import * as THREE from "three";
import { useScrollProgressRef } from "@/components/scroll/ScrollProgressContext";
import { sectionOffsetAmount } from "@/lib/timeline";
import { damp } from "@/lib/scrollMath";

/**
 * Keeps everything drawn in the section's own coordinates travelling with the
 * section.
 *
 * The transport furniture -- pressure contours, the concentration field,
 * vessel stand-ins, the fronts leaving them, the concentration source -- is
 * all positioned in the measured section's world coordinates, which is the
 * whole reason it registers with the tissue. But the tissue itself is moved
 * aside in the final scene to clear room for the CTA (see TumorModel's
 * children: `position.x = sectionOffsetAmount(p) * 1.15`), and anything that
 * does not move with it detaches and floats in empty space -- which is
 * exactly what the vessels used to do.
 *
 * These children cannot simply be parented to the tumour: the real asset is
 * loaded by GLTFTumor, which takes no scene children, and hanging them off it
 * would also inherit its continuous sway. So the offset is re-derived here
 * from the same pure timeline function with the same damping constant.
 * `sectionOffsetAmount` is a function of scroll progress only, so the two stay
 * in lockstep frame for frame rather than merely close.
 *
 * Only the translation is mirrored, deliberately. GLTFTumor's shallow sway is
 * an accumulated per-frame angle, not a function of progress, and would drift
 * out of phase the moment the model suspended on its own GLB load; a few
 * degrees of parallax between illustrative overlays and the tissue is a much
 * smaller lie than a metre of lateral drift.
 */
export default function TransportOverlay({ children }: { children?: ReactNode }) {
  const progressRef = useScrollProgressRef();
  const groupRef = useRef<THREE.Group>(null);
  const offsetSmoothed = useRef(0);

  useFrame((state, delta) => {
    // Same aspect the tissue and the camera read, so the three stay in
    // lockstep at every viewport shape rather than only at the authored one.
    const aspect = state.size.width / Math.max(state.size.height, 1);
    offsetSmoothed.current = damp(
      offsetSmoothed.current,
      sectionOffsetAmount(progressRef.current, aspect),
      3.2,
      delta,
    );
    if (groupRef.current) groupRef.current.position.x = offsetSmoothed.current * 1.15;
  });

  return <group ref={groupRef}>{children}</group>;
}
