"use client";

import { useMemo, useRef } from "react";
import { useFrame } from "@react-three/fiber";
import * as THREE from "three";
import { useScrollProgressRef } from "@/components/scroll/ScrollProgressContext";
import { particleOpacity } from "@/lib/timeline";
import { damp } from "@/lib/scrollMath";
import { mulberry32 } from "@/lib/prng";

/**
 * Ambient dust/particulate field for depth in the hero and early scenes.
 * Count is halved on mobile per the spec's responsive performance
 * requirements. Pure decoration — opacity fades down (never fully to zero,
 * for continued depth) as the narrative moves past the hero.
 *
 * `reducedMotion` skips the per-frame work entirely rather than slowing it.
 * Both animated quantities are constant in that mode anyway — the rotation is
 * a function of the wall clock, which is the one thing the fallback exists to
 * stop, and the opacity is a function of scroll progress, which is pinned at 0
 * there because the fallback has no scroll timeline. Not subscribing at all is
 * what lets SceneCanvas put the whole canvas on `frameloop="demand"`: a single
 * `useFrame` anywhere in the tree is enough to keep a "static" scene
 * recomputing forever.
 */
export default function Particles({
  count = 900,
  reducedMotion = false,
}: {
  count?: number;
  reducedMotion?: boolean;
}) {
  const progressRef = useScrollProgressRef();
  const pointsRef = useRef<THREE.Points>(null);
  const materialRef = useRef<THREE.PointsMaterial>(null);
  const opacitySmoothed = useRef(1);

  const positions = useMemo(() => {
    const random = mulberry32(1337);
    const arr = new Float32Array(count * 3);
    for (let i = 0; i < count; i++) {
      const radius = 2.4 + random() * 3.2;
      const theta = random() * Math.PI * 2;
      const phi = Math.acos(2 * random() - 1);
      arr[i * 3] = radius * Math.sin(phi) * Math.cos(theta);
      arr[i * 3 + 1] = radius * Math.sin(phi) * Math.sin(theta) * 0.6;
      arr[i * 3 + 2] = radius * Math.cos(phi);
    }
    return arr;
  }, [count]);

  useFrame((state, delta) => {
    if (reducedMotion) return;
    const target = particleOpacity(progressRef.current);
    opacitySmoothed.current = damp(opacitySmoothed.current, target, 2, delta);
    if (materialRef.current) materialRef.current.opacity = opacitySmoothed.current * 0.45;
    if (pointsRef.current) {
      pointsRef.current.rotation.y = state.clock.elapsedTime * 0.015;
    }
  });

  return (
    <points ref={pointsRef}>
      <bufferGeometry>
        <bufferAttribute attach="attributes-position" args={[positions, 3]} />
      </bufferGeometry>
      <pointsMaterial
        ref={materialRef}
        color="#cfe8f0"
        size={0.014}
        sizeAttenuation
        transparent
        opacity={0.45}
        depthWrite={false}
        blending={THREE.AdditiveBlending}
      />
    </points>
  );
}
