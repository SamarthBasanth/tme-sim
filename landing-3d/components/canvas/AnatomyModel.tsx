"use client";

import { useRef } from "react";
import { useFrame } from "@react-three/fiber";
import * as THREE from "three";
import { useScrollProgressRef } from "@/components/scroll/ScrollProgressContext";
import { anatomyOpacity } from "@/lib/timeline";
import { damp } from "@/lib/scrollMath";

/**
 * Simplified translucent anatomical shell surrounding the tumor — stands in
 * for a brain/organ boundary. Opacity is driven entirely by scroll (fades
 * in through "Anatomical Context", recedes again toward the final scene),
 * never toggled — so it reads as the camera revealing context, not a UI
 * element appearing.
 */
export default function AnatomyModel({ lowDetail = false }: { lowDetail?: boolean }) {
  const progressRef = useScrollProgressRef();
  const materialRef = useRef<THREE.MeshPhysicalMaterial>(null);
  const groupRef = useRef<THREE.Group>(null);
  const opacitySmoothed = useRef(0);

  useFrame((_, delta) => {
    const target = anatomyOpacity(progressRef.current);
    opacitySmoothed.current = damp(opacitySmoothed.current, target, 3, delta);
    if (materialRef.current) {
      materialRef.current.opacity = opacitySmoothed.current;
    }
    if (groupRef.current) {
      groupRef.current.visible = opacitySmoothed.current > 0.005;
      groupRef.current.rotation.y += delta * 0.02;
    }
  });

  return (
    <group ref={groupRef} scale={[1.9, 2.3, 1.7]} position={[0.1, 0.1, -0.1]}>
      <mesh>
        <icosahedronGeometry args={[1, lowDetail ? 2 : 4]} />
        <meshPhysicalMaterial
          ref={materialRef}
          color="#dfe6ea"
          transparent
          opacity={0}
          transmission={0.55}
          roughness={0.35}
          thickness={0.6}
          ior={1.2}
          side={THREE.DoubleSide}
          depthWrite={false}
        />
      </mesh>
    </group>
  );
}
