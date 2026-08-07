"use client";

import { useMemo, useRef, type ReactNode } from "react";
import { useFrame } from "@react-three/fiber";
import { Line } from "@react-three/drei";
import type { Line2, LineSegments2 } from "three-stdlib";
import * as THREE from "three";
import { useScrollProgressRef } from "@/components/scroll/ScrollProgressContext";
import { fieldContourOpacity } from "@/lib/timeline";
import { damp } from "@/lib/scrollMath";
import { fbm3 } from "@/lib/noise";
import { mulberry32 } from "@/lib/prng";
import { SECTION_HALF_EXTENT } from "./shaders/viridis";

interface ContourSpec {
  /** Fraction of the section's extent this iso-level encloses. */
  level: number;
  seed: number;
}

/**
 * Four nested iso-level contours of the interstitial pressure field.
 *
 * This replaces a set of axial / sagittal / coronal slice outlines -- the
 * three-plane convention a radiology volume viewer uses -- which was the right
 * furniture for a product that read scans and the wrong furniture for one that
 * solves transport. Three orthogonal cutting planes around a flat measured
 * section were also geometrically incoherent: the asset has no third axis to
 * cut through, so the old rings enclosed it like a wireframe globe.
 *
 * Nested closed contours in the plane of the section are what a scalar field
 * over a 2D section actually looks like, and interstitial pressure is a field
 * this product solves (high in the interior where drainage fails, falling
 * toward the rim -- which is why the outer levels are the ones near the
 * section edge). Nothing here is a solved field; it is contour furniture, and
 * the annotation that points at it says so.
 */
const CONTOUR_SPECS: readonly ContourSpec[] = [
  { level: 1.0, seed: 47 },
  { level: 0.74, seed: 11 },
  { level: 0.5, seed: 83 },
  { level: 0.29, seed: 129 },
];

const SEGMENTS = 84;

/**
 * The pressure peak sits inside the tissue rather than at the geometric centre
 * of the bounding box -- a field maximum is a property of the tissue, not of
 * the crop. Offset in the section's own world units.
 */
const PEAK = { x: -0.18, y: -0.06 };

/** Drop the innermost level on mobile, per the responsive detail budget. */
const LOW_DETAIL_CONTOURS = 3;

/**
 * One closed iso-contour: an ellipse around the field peak, with its radius
 * walked by noise around the loop so it reads as a traced level set of a
 * lumpy field rather than a geometric primitive. Inner levels are jittered
 * more, the way a contour near a peak is rougher than a smooth outer one.
 */
function buildContour(spec: ContourSpec): Array<[number, number, number]> {
  const rand = mulberry32(spec.seed);
  const phase = rand() * 30;
  const jitter = 0.06 + (1 - spec.level) * 0.07;
  const pts: Array<[number, number, number]> = [];
  for (let i = 0; i <= SEGMENTS; i++) {
    const t = (i / SEGMENTS) * Math.PI * 2;
    const n = fbm3(Math.cos(t) * 1.6 + phase, Math.sin(t) * 1.6 + phase, spec.level * 2 + phase, 3);
    const rx = SECTION_HALF_EXTENT.x * 0.88 * spec.level * (1 + n * jitter);
    const ry = SECTION_HALF_EXTENT.y * 0.88 * spec.level * (1 + n * jitter);
    pts.push([PEAK.x + Math.cos(t) * rx, PEAK.y + Math.sin(t) * ry, 0]);
  }
  return pts;
}

interface FieldContoursProps {
  lowDetail?: boolean;
  /**
   * Rendered INSIDE the same swaying group as the contours, not as a sibling
   * in the parent scene. This group sways every frame, so any annotation
   * meant to sit on a specific contour point must be a child and let
   * Three.js's normal parent-child transform composition keep it correct,
   * instead of a hand-computed world coordinate that is only ever right at
   * the one instant it was computed for. See SceneCanvas.tsx for the anchor
   * point itself (in this group's LOCAL space).
   */
  children?: ReactNode;
}

/**
 * Just in front of the section plane. Behind it, the 40k-point cloud swallowed
 * the lines almost entirely and the annotation pointed at nothing; a field
 * overlay is drawn ON a section anyway, which is where the Studio draws it.
 */
const CONTOUR_Z = 0.04;

export default function FieldContours({ lowDetail = false, children }: FieldContoursProps) {
  const progressRef = useScrollProgressRef();
  const groupRef = useRef<THREE.Group>(null);
  const lineRefs = useRef<Array<Line2 | LineSegments2 | null>>([]);
  const opacitySmoothed = useRef(0);

  const contours = useMemo(() => {
    const specs = lowDetail ? CONTOUR_SPECS.slice(0, LOW_DETAIL_CONTOURS) : CONTOUR_SPECS;
    return specs.map(buildContour);
  }, [lowDetail]);

  useFrame((state, delta) => {
    const target = fieldContourOpacity(progressRef.current);
    opacitySmoothed.current = damp(opacitySmoothed.current, target, 3, delta);
    const o = opacitySmoothed.current;
    lineRefs.current.forEach((line) => {
      if (!line) return;
      const material = line.material as THREE.Material & { opacity: number; transparent: boolean };
      material.opacity = o * 0.85;
      material.transparent = true;
    });
    if (groupRef.current) {
      groupRef.current.visible = o > 0.005;
      // The same shallow sway law GLTFTumor gives a planar asset (a flat thing
      // must never turn edge-on or it vanishes), so the contours stay roughly
      // registered with the section they annotate instead of spinning free of
      // it. Driven from the clock rather than shared state: these are
      // illustrative contours, and a few degrees of drift costs nothing.
      groupRef.current.rotation.y = Math.sin(state.clock.elapsedTime * 0.22) * 0.15;
    }
  });

  return (
    <group ref={groupRef} position={[0, 0, CONTOUR_Z]}>
      {contours.map((points, i) => (
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
      {children}
    </group>
  );
}
