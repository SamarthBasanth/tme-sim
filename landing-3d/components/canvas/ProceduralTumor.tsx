"use client";

import { useMemo, useRef } from "react";
import { useFrame } from "@react-three/fiber";
import { MeshDistortMaterial } from "@react-three/drei";
import * as THREE from "three";
import { useScrollProgressRef } from "@/components/scroll/ScrollProgressContext";
import { finalOffsetAmount, separationAmount } from "@/lib/timeline";
import { damp } from "@/lib/scrollMath";

interface ProceduralTumorProps {
  lowDetail?: boolean;
  /** Slows/removes scroll-coupled and ambient rotation for prefers-reduced-motion. */
  reducedMotion?: boolean;
}

const CORE_COLOR = "#b23a3a";
const EDEMA_COLOR = "#8fd6e8";
const MARGIN_COLOR = "#e7ecef";

/**
 * Attractive stand-in for a real `/models/tumor.glb`: a cluster of distorted,
 * organically-noised spheres (drei's `MeshDistortMaterial`) rather than a
 * single perfect sphere, so the silhouette reads as an irregular mass. It is
 * built from three concentric groups — core / edema / margin — because the
 * "Segmentation" scene needs to pull those apart into distinct, separately
 * materialed regions; a single fused mesh could not do that. Swap this
 * component out (see TumorModel.tsx) once a real scan-derived GLB exists.
 */
export default function ProceduralTumor({
  lowDetail = false,
  reducedMotion = false,
}: ProceduralTumorProps) {
  const progressRef = useScrollProgressRef();

  const coreGroup = useRef<THREE.Group>(null);
  const edemaGroup = useRef<THREE.Group>(null);
  const marginGroup = useRef<THREE.Group>(null);
  const rootGroup = useRef<THREE.Group>(null);

  const sepSmoothed = useRef(0);
  const offsetSmoothed = useRef(0);
  const autoRotation = useRef(0);
  const scrollRotation = useRef(0);

  const detail = lowDetail ? 3 : 5;

  // Small off-axis "lobes" bolted onto the core so the mass silhouette reads
  // as irregular/organic rather than a perfect sphere — a cheap approximation
  // of metaball blending without the cost of an actual iso-surface march.
  const lobes = useMemo(
    () => [
      { position: [0.42, 0.18, 0.1] as const, scale: 0.46 },
      { position: [-0.32, -0.28, 0.22] as const, scale: 0.4 },
      { position: [0.08, 0.4, -0.3] as const, scale: 0.36 },
    ],
    [],
  );

  useFrame((state, delta) => {
    const p = progressRef.current;
    sepSmoothed.current = damp(sepSmoothed.current, separationAmount(p), 3.2, delta);
    offsetSmoothed.current = damp(offsetSmoothed.current, finalOffsetAmount(p), 3.2, delta);

    const sep = sepSmoothed.current;
    if (coreGroup.current) coreGroup.current.position.set(0, -0.55 * sep, 0);
    if (edemaGroup.current) edemaGroup.current.position.set(0.05 * sep, 0.62 * sep, 0.1 * sep);
    if (marginGroup.current) marginGroup.current.position.set(-0.08 * sep, 0.05 * sep, -0.6 * sep);

    if (rootGroup.current) {
      // Gentle continuous autonomous rotation, plus a slow scroll-coupled
      // turn so the model always feels connected to scroll input even
      // within a single section (not just at section boundaries). Reduced
      // to a slow idle spin with no scroll coupling under reduced-motion.
      autoRotation.current += delta * (reducedMotion ? 0.035 : 0.12);
      scrollRotation.current = reducedMotion
        ? 0
        : damp(scrollRotation.current, p * 0.6, 2.4, delta);
      rootGroup.current.rotation.y = autoRotation.current + scrollRotation.current;
      rootGroup.current.position.x = reducedMotion ? 0 : offsetSmoothed.current * 1.15;
    }

    void state;
  });

  return (
    <group ref={rootGroup}>
      {/* Core — dense, saturated, the primary mass */}
      <group ref={coreGroup}>
        <mesh castShadow={false} receiveShadow={false}>
          <icosahedronGeometry args={[0.62, detail]} />
          <MeshDistortMaterial
            color={CORE_COLOR}
            distort={0.32}
            speed={1.1}
            roughness={0.35}
            metalness={0.1}
            emissive={CORE_COLOR}
            emissiveIntensity={0.12}
          />
        </mesh>
        {lobes.map((lobe, i) => (
          <mesh key={i} position={lobe.position} scale={lobe.scale}>
            <icosahedronGeometry args={[0.62, Math.max(2, detail - 1)]} />
            <MeshDistortMaterial
              color={CORE_COLOR}
              distort={0.4}
              speed={0.9}
              roughness={0.4}
              metalness={0.08}
              emissive={CORE_COLOR}
              emissiveIntensity={0.1}
            />
          </mesh>
        ))}
      </group>

      {/* Edema — translucent halo surrounding the core */}
      <group ref={edemaGroup}>
        <mesh>
          <icosahedronGeometry args={[0.98, detail]} />
          <MeshDistortMaterial
            color={EDEMA_COLOR}
            distort={0.22}
            speed={0.6}
            roughness={0.2}
            metalness={0}
            transparent
            opacity={0.22}
            depthWrite={false}
          />
        </mesh>
      </group>

      {/* Margin — faint outer boundary of surrounding tissue */}
      <group ref={marginGroup}>
        <mesh>
          <icosahedronGeometry args={[1.28, Math.max(2, detail - 1)]} />
          <meshStandardMaterial
            color={MARGIN_COLOR}
            wireframe
            transparent
            opacity={0.08}
            depthWrite={false}
          />
        </mesh>
      </group>
    </group>
  );
}
