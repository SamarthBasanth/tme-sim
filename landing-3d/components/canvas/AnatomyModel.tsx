"use client";

import { useMemo, useRef } from "react";
import { useFrame } from "@react-three/fiber";
import { Line } from "@react-three/drei";
import type { Line2, LineSegments2 } from "three-stdlib";
import * as THREE from "three";
import { useScrollProgressRef } from "@/components/scroll/ScrollProgressContext";
import { anatomyOpacity } from "@/lib/timeline";
import { damp } from "@/lib/scrollMath";
import { fbm3 } from "@/lib/noise";
import { mulberry32 } from "@/lib/prng";

type Plane = "axial" | "sagittal" | "coronal";

interface RingSpec {
  plane: Plane;
  level: number;
  radiusA: number;
  radiusB: number;
  seed: number;
}

// Three axial (horizontal) slice contours at different heights, plus one
// sagittal and one coronal reference plane — the standard three-plane
// convention every radiology/volume viewer uses to give anatomical context,
// rather than one enclosing globe.
const RING_SPECS: readonly RingSpec[] = [
  { plane: "axial", level: -0.52, radiusA: 0.78, radiusB: 0.7, seed: 11 },
  { plane: "axial", level: 0.02, radiusA: 0.94, radiusB: 0.88, seed: 47 },
  { plane: "axial", level: 0.58, radiusA: 0.6, radiusB: 0.52, seed: 83 },
  { plane: "sagittal", level: 0.04, radiusA: 0.9, radiusB: 0.86, seed: 129 },
  { plane: "coronal", level: -0.05, radiusA: 0.88, radiusB: 0.9, seed: 171 },
];

const SEGMENTS = 72;

/**
 * One irregular slice-contour ring — an approximation of a segmentation
 * outline traced on a single axial/sagittal/coronal plane, not a perfect
 * circle: radius is jittered by a noise walk around the loop so it reads as
 * traced from real anatomy rather than a geometric primitive.
 */
function buildRing(spec: RingSpec): Array<[number, number, number]> {
  const rand = mulberry32(spec.seed);
  const phase = rand() * 30;
  const pts: Array<[number, number, number]> = [];
  for (let i = 0; i <= SEGMENTS; i++) {
    const t = (i / SEGMENTS) * Math.PI * 2;
    const n = fbm3(Math.cos(t) * 1.6 + phase, Math.sin(t) * 1.6 + phase, spec.level * 2 + phase, 3);
    const rA = spec.radiusA * (1 + n * 0.07);
    const rB = spec.radiusB * (1 + n * 0.07);
    const a = Math.cos(t) * rA;
    const b = Math.sin(t) * rB;
    if (spec.plane === "axial") pts.push([a, spec.level, b]);
    else if (spec.plane === "sagittal") pts.push([spec.level, a, b]);
    else pts.push([a, b, spec.level]);
  }
  return pts;
}

/**
 * Anatomical-context indicator, reworked away from a translucent faceted
 * icosahedron shell (read as a generic sci-fi "network globe") toward
 * something grounded in how real medical volume viewers actually show
 * anatomical context: a sparse set of axial/sagittal/coronal slice-contour
 * outlines (see tumorGeometry.ts's BraTS/radiology-viewer citations), each
 * an irregular traced loop rather than a geometric primitive. Opacity is
 * still driven entirely by scroll (fades in through "Anatomical Context",
 * recedes again toward the final scene), never toggled, exactly as the
 * shell it replaces.
 */
export default function AnatomyModel({ lowDetail = false }: { lowDetail?: boolean }) {
  const progressRef = useScrollProgressRef();
  const groupRef = useRef<THREE.Group>(null);
  const lineRefs = useRef<Array<Line2 | LineSegments2 | null>>([]);
  const opacitySmoothed = useRef(0);

  const rings = useMemo(() => {
    const specs = lowDetail ? RING_SPECS.slice(0, 3) : RING_SPECS;
    return specs.map(buildRing);
  }, [lowDetail]);

  useFrame((_, delta) => {
    const target = anatomyOpacity(progressRef.current);
    opacitySmoothed.current = damp(opacitySmoothed.current, target, 3, delta);
    const o = opacitySmoothed.current;
    lineRefs.current.forEach((line) => {
      if (!line) return;
      const material = line.material as THREE.Material & { opacity: number; transparent: boolean };
      material.opacity = o * 0.62;
      material.transparent = true;
    });
    if (groupRef.current) {
      groupRef.current.visible = o > 0.005;
      groupRef.current.rotation.y += delta * 0.02;
    }
  });

  return (
    <group ref={groupRef} scale={[1.9, 2.3, 1.7]} position={[0.1, 0.1, -0.1]}>
      {rings.map((points, i) => (
        <Line
          key={i}
          ref={(el) => {
            lineRefs.current[i] = el;
          }}
          points={points}
          color="#8fd8e8"
          lineWidth={1.1}
          transparent
          opacity={0}
        />
      ))}
    </group>
  );
}
