"use client";

import { useMemo, useRef } from "react";
import { useFrame } from "@react-three/fiber";
import * as THREE from "three";
import { useScrollProgressRef } from "@/components/scroll/ScrollProgressContext";
import { criticalStructureOpacity, highlightOpacity } from "@/lib/timeline";
import { damp } from "@/lib/scrollMath";

/** Small emissive marker that "ignites" as the AI-analysis scan plane sweeps past it. */
export function HighlightMarker({ position }: { position: readonly [number, number, number] }) {
  const progressRef = useScrollProgressRef();
  const materialRef = useRef<THREE.MeshStandardMaterial>(null);
  const ringRef = useRef<THREE.Mesh>(null);
  const opacitySmoothed = useRef(0);

  useFrame((state, delta) => {
    const target = highlightOpacity(progressRef.current);
    opacitySmoothed.current = damp(opacitySmoothed.current, target, 4, delta);
    const o = opacitySmoothed.current;
    if (materialRef.current) {
      materialRef.current.opacity = o;
      materialRef.current.emissiveIntensity = 0.6 + o * 1.4;
    }
    if (ringRef.current) {
      ringRef.current.scale.setScalar(1 + Math.sin(state.clock.elapsedTime * 2) * 0.06 * o);
      const ringMat = ringRef.current.material as THREE.MeshBasicMaterial;
      ringMat.opacity = o * 0.5;
    }
  });

  return (
    <group position={position}>
      <mesh>
        <sphereGeometry args={[0.05, 16, 16]} />
        <meshStandardMaterial
          ref={materialRef}
          color="#f2b23a"
          emissive="#f2b23a"
          transparent
          opacity={0}
          depthWrite={false}
        />
      </mesh>
      <mesh ref={ringRef} rotation={[Math.PI / 2, 0, 0]}>
        <ringGeometry args={[0.09, 0.11, 32]} />
        <meshBasicMaterial color="#f2b23a" transparent opacity={0} side={THREE.DoubleSide} depthWrite={false} />
      </mesh>
    </group>
  );
}

/** Faint translucent structures (stand-ins for nearby vessels) revealed during treatment planning. */
export function CriticalStructures() {
  const progressRef = useScrollProgressRef();
  const groupRef = useRef<THREE.Group>(null);
  const opacitySmoothed = useRef(0);

  const curves = useMemo(
    () => [
      new THREE.CatmullRomCurve3([
        new THREE.Vector3(1.5, 0.9, -0.3),
        new THREE.Vector3(0.9, 0.5, -0.5),
        new THREE.Vector3(0.3, -0.1, -0.4),
        new THREE.Vector3(-0.4, -0.6, -0.3),
      ]),
      new THREE.CatmullRomCurve3([
        new THREE.Vector3(-1.3, -0.7, 0.6),
        new THREE.Vector3(-0.7, -0.3, 0.5),
        new THREE.Vector3(-0.1, 0.2, 0.4),
      ]),
    ],
    [],
  );

  const geometries = useMemo(
    () => curves.map((curve) => new THREE.TubeGeometry(curve, 32, 0.014, 6, false)),
    [curves],
  );

  useFrame((_, delta) => {
    const target = criticalStructureOpacity(progressRef.current);
    opacitySmoothed.current = damp(opacitySmoothed.current, target, 3, delta);
    if (groupRef.current) {
      groupRef.current.children.forEach((child) => {
        const mesh = child as THREE.Mesh;
        const material = mesh.material as THREE.MeshStandardMaterial;
        material.opacity = opacitySmoothed.current * 0.5;
      });
    }
  });

  return (
    <group ref={groupRef}>
      {geometries.map((geometry, i) => (
        <mesh key={i} geometry={geometry}>
          <meshStandardMaterial
            color="#9fc6d8"
            transparent
            opacity={0}
            depthWrite={false}
            roughness={0.4}
          />
        </mesh>
      ))}
    </group>
  );
}
