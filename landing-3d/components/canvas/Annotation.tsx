"use client";

import { useRef } from "react";
import { useFrame } from "@react-three/fiber";
import { Html, Line } from "@react-three/drei";
import type { Line2, LineSegments2 } from "three-stdlib";
import * as THREE from "three";
import { damp } from "@/lib/scrollMath";
import { useScrollProgressRef } from "@/components/scroll/ScrollProgressContext";

interface AnnotationProps {
  anchor: readonly [number, number, number];
  labelOffset: readonly [number, number, number];
  label: string;
  detail?: string;
  /** Pure function of global scroll progress -> 0..1 visibility. */
  getOpacity: (progress: number) => number;
  align?: "left" | "right";
  /**
   * Render the reticle ONLY — no leader line, no HTML label.
   *
   * Exists for narrow viewports. The label is a projection of a world-space
   * point, so its screen position sweeps as the camera moves and cannot be
   * kept out of a copy column that runs full-bleed; measured at 390px it
   * overlapped beat 02's body text at every scroll position where both were
   * visible. Rather than hand-tune a 3D offset against a moving projection,
   * the label is dropped where there is no room for it and the beat's own
   * in-flow ledger carries the words — which it does far better, naming the
   * quantity AND its three registered verdicts instead of two lines of chip.
   *
   * The reticle survives because it is the part that still works without
   * space: it marks the anchor on the tissue. The leader line does not,
   * because a line routed to a label that is not there points at nothing.
   */
  reticleOnly?: boolean;
}

/**
 * A clinical-instrument-style callout: a small reticle marking the exact
 * anchor point on the tumour, and an elbowed (dogleg) leader line into the
 * label rather than one straight diagonal spoke. This is deliberately not
 * "a line radiating from a hub to a node" — the convention it replaces,
 * which combined with the old wireframe shell read as a generic network
 * globe — but the routed-callout convention used for point annotations in
 * technical and medical figures. Reads scroll progress itself and derives
 * its own visibility via `getOpacity` (one of the authored curves in
 * lib/timeline.ts), eased here per-frame so callouts never pop.
 */
export default function Annotation({
  anchor,
  labelOffset,
  label,
  detail,
  getOpacity,
  align = "right",
  reticleOnly = false,
}: AnnotationProps) {
  const progressRef = useScrollProgressRef();
  const groupOpacity = useRef(0);
  const lineRef = useRef<Line2 | LineSegments2>(null);
  const reticleRef = useRef<THREE.Mesh>(null);
  const htmlWrapRef = useRef<HTMLDivElement>(null);

  const labelPoint: [number, number, number] = [
    anchor[0] + labelOffset[0],
    anchor[1] + labelOffset[1],
    anchor[2] + labelOffset[2],
  ];
  // Dogleg: travel horizontally from the anchor first, then into the label
  // — a leader line that reads as deliberately routed rather than a spoke
  // pointing straight back at a central hub.
  const elbow: [number, number, number] = [labelPoint[0], anchor[1], anchor[2]];

  useFrame((_, delta) => {
    groupOpacity.current = damp(groupOpacity.current, getOpacity(progressRef.current), 5, delta);
    const o = groupOpacity.current;
    if (lineRef.current) {
      const material = lineRef.current.material as THREE.Material & { opacity: number; transparent: boolean };
      material.opacity = o * 0.85;
      material.transparent = true;
    }
    if (reticleRef.current) {
      const material = reticleRef.current.material as THREE.Material & { opacity: number };
      material.opacity = o * 0.9;
      reticleRef.current.scale.setScalar(0.85 + o * 0.15);
    }
    if (htmlWrapRef.current) {
      htmlWrapRef.current.style.opacity = String(o);
      htmlWrapRef.current.style.transform = `translateY(${(1 - o) * 6}px)`;
    }
  });

  return (
    <>
      <mesh ref={reticleRef} position={anchor} rotation={[Math.PI / 2, 0, 0]}>
        <ringGeometry args={[0.026, 0.036, 24]} />
        <meshBasicMaterial
          color="#bdeaf7"
          transparent
          opacity={0}
          depthWrite={false}
          side={THREE.DoubleSide}
        />
      </mesh>
      {!reticleOnly && (
        <>
          <Line
            ref={lineRef}
            points={[anchor, elbow, labelPoint]}
            color="#bdeaf7"
            lineWidth={1.25}
            transparent
            opacity={0}
          />
          <Html position={labelPoint} center={false} zIndexRange={[10, 0]} sprite={false}>
            <div
              ref={htmlWrapRef}
              className={`annotation-label ${align === "left" ? "annotation-label--left" : ""}`}
              style={{ opacity: 0 }}
            >
              <span className="annotation-label__title">{label}</span>
              {detail ? <span className="annotation-label__detail">{detail}</span> : null}
            </div>
          </Html>
        </>
      )}
    </>
  );
}
