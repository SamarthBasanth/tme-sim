"use client";

import { useEffect, useMemo, useRef } from "react";
import { useFrame } from "@react-three/fiber";
import { Billboard } from "@react-three/drei";
import * as THREE from "three";
import { useScrollProgressRef } from "@/components/scroll/ScrollProgressContext";
import { criticalStructureOpacity, highlightOpacity } from "@/lib/timeline";
import { damp } from "@/lib/scrollMath";

/**
 * Small emissive marker that "ignites" as the AI-analysis scan plane sweeps
 * past it — a flagged-region indicator during AI Analysis / Treatment
 * Planning, not decoration.
 *
 * Previously a filled sphere with a second ring orbiting it in a fixed
 * horizontal plane (`rotation={[Math.PI / 2, 0, 0]}`): from most camera
 * angles a plane-locked ring around a small sphere reads exactly like a
 * miniature Saturn rather than a clinical marker — flagged (and left
 * unresolved) in the prior visual pass, then independently re-reported.
 * Rebuilt as a camera-facing (`Billboard`, so it never tilts into a
 * Saturn-ring silhouette regardless of viewing angle) crosshair/reticle: a
 * small pulsing dot, a thin ring tight against it, and four short radial
 * tick marks — the standard "flagged point" convention in medical and
 * technical imaging overlays.
 */
export function HighlightMarker({ position }: { position: readonly [number, number, number] }) {
  const progressRef = useScrollProgressRef();
  const dotMaterialRef = useRef<THREE.MeshBasicMaterial>(null);
  const ringMaterialRef = useRef<THREE.MeshBasicMaterial>(null);
  const tickMaterialRef = useRef<THREE.LineBasicMaterial>(null);
  const billboardRef = useRef<THREE.Group>(null);
  const opacitySmoothed = useRef(0);

  // Four short radial ticks just outside the ring, at 45/135/225/315deg —
  // a crosshair/reticle silhouette, built once as a single LineSegments
  // geometry (disconnected pairs) rather than four separate Line elements.
  const tickGeometry = useMemo(() => {
    const inner = 0.052;
    const outer = 0.08;
    const positions: number[] = [];
    for (let i = 0; i < 4; i++) {
      const angle = (i / 4) * Math.PI * 2 + Math.PI / 4;
      const cx = Math.cos(angle);
      const cy = Math.sin(angle);
      positions.push(cx * inner, cy * inner, 0, cx * outer, cy * outer, 0);
    }
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute("position", new THREE.Float32BufferAttribute(positions, 3));
    return geometry;
  }, []);

  useEffect(() => () => tickGeometry.dispose(), [tickGeometry]);

  useFrame((state, delta) => {
    const target = highlightOpacity(progressRef.current);
    opacitySmoothed.current = damp(opacitySmoothed.current, target, 4, delta);
    const o = opacitySmoothed.current;
    const pulse = 1 + Math.sin(state.clock.elapsedTime * 2.4) * 0.08 * o;

    if (dotMaterialRef.current) dotMaterialRef.current.opacity = o;
    if (ringMaterialRef.current) ringMaterialRef.current.opacity = o * 0.85;
    if (tickMaterialRef.current) tickMaterialRef.current.opacity = o * 0.7;
    if (billboardRef.current) billboardRef.current.scale.setScalar(pulse);
  });

  return (
    <group position={position}>
      <Billboard ref={billboardRef}>
        <mesh>
          <circleGeometry args={[0.03, 20]} />
          <meshBasicMaterial ref={dotMaterialRef} color="#f2b23a" transparent opacity={0} depthWrite={false} />
        </mesh>
        <mesh>
          <ringGeometry args={[0.044, 0.05, 32]} />
          <meshBasicMaterial
            ref={ringMaterialRef}
            color="#f2b23a"
            transparent
            opacity={0}
            depthWrite={false}
            side={THREE.DoubleSide}
          />
        </mesh>
        <lineSegments geometry={tickGeometry}>
          <lineBasicMaterial ref={tickMaterialRef} color="#f2b23a" transparent opacity={0} depthWrite={false} />
        </lineSegments>
      </Billboard>
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
