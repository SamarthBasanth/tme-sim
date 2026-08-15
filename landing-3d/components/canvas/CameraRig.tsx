"use client";

import { useRef } from "react";
import { useFrame } from "@react-three/fiber";
import * as THREE from "three";
import { getCameraPose } from "@/lib/timeline";
import { damp3 } from "@/lib/scrollMath";
import { useScrollProgressRef } from "@/components/scroll/ScrollProgressContext";

/**
 * Drives the camera purely from scroll progress. Position, look-at target,
 * and field of view are looked up from `getCameraPose(progress)` — an
 * explicit function of scroll, not a triggered tween — and then eased
 * toward with per-frame exponential damping so the camera never jitters or
 * snaps even if `progress` updates in an uneven burst. Because the target
 * pose is a pure function of `progress`, scrolling backward simply drives
 * the camera back along the same path with no special-casing.
 *
 * The camera is read from `useFrame`'s own `state.camera` (R3F's per-frame
 * loop, which runs entirely outside React's render/commit cycle) rather
 * than destructured from `useThree()` at the top of the component — R3F's
 * whole rendering model is imperative mutation of scene objects inside
 * this loop, so the object is intentionally mutated here every frame.
 */
export default function CameraRig() {
  const progressRef = useScrollProgressRef();

  const smoothedPos = useRef<[number, number, number]>([0, 0, 4.4]);
  const smoothedLookAt = useRef<[number, number, number]>([0, 0, 0]);
  const smoothedFov = useRef(32);
  const lookAtTarget = useRef(new THREE.Vector3());

  useFrame((state, delta) => {
    const { camera } = state;
    /*
     * The pose is aspect-corrected, not raw. `fov` in three.js is vertical, so
     * a portrait viewport sees a fraction of the horizontal extent the
     * keyframes were composed against — see the aspect block in lib/timeline.
     * Read from `state.size` each frame rather than captured once, so a
     * rotation or a resize recomposes rather than staying wrong until reload.
     */
    const aspect = state.size.width / Math.max(state.size.height, 1);
    const pose = getCameraPose(progressRef.current, aspect);
    const lambda = 4.2;

    smoothedPos.current = damp3(smoothedPos.current, pose.position, lambda, delta);
    smoothedLookAt.current = damp3(smoothedLookAt.current, pose.lookAt, lambda, delta);
    smoothedFov.current = THREE.MathUtils.damp(smoothedFov.current, pose.fov, lambda, delta);

    camera.position.set(...smoothedPos.current);
    lookAtTarget.current.set(...smoothedLookAt.current);
    camera.lookAt(lookAtTarget.current);

    if (camera instanceof THREE.PerspectiveCamera) {
      if (Math.abs(camera.fov - smoothedFov.current) > 0.001) {
        camera.fov = smoothedFov.current;
        camera.updateProjectionMatrix();
      }
    }
  });

  return null;
}
