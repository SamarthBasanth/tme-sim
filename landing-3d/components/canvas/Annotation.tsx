"use client";

import { useRef } from "react";
import { useFrame } from "@react-three/fiber";
import { Html, Line } from "@react-three/drei";
import type { Line2, LineMaterial, LineSegments2 } from "three-stdlib";
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
}

/**
 * A clinical-style leader line + label anchored to a 3D point (tumor core,
 * edema margin, ...). Reads scroll progress itself and derives its own
 * visibility via `getOpacity` (one of the authored curves in
 * lib/timeline.ts), eased here per-frame so labels never pop.
 */
export default function Annotation({
  anchor,
  labelOffset,
  label,
  detail,
  getOpacity,
  align = "right",
}: AnnotationProps) {
  const progressRef = useScrollProgressRef();
  const groupOpacity = useRef(0);
  const lineRef = useRef<Line2 | LineSegments2>(null);
  const htmlWrapRef = useRef<HTMLDivElement>(null);

  const labelPoint: [number, number, number] = [
    anchor[0] + labelOffset[0],
    anchor[1] + labelOffset[1],
    anchor[2] + labelOffset[2],
  ];

  useFrame((_, delta) => {
    groupOpacity.current = damp(groupOpacity.current, getOpacity(progressRef.current), 5, delta);
    const o = groupOpacity.current;
    if (lineRef.current) {
      const material = lineRef.current.material as LineMaterial;
      material.opacity = o * 0.9;
      material.transparent = true;
    }
    if (htmlWrapRef.current) {
      htmlWrapRef.current.style.opacity = String(o);
      htmlWrapRef.current.style.transform = `translateY(${(1 - o) * 6}px)`;
    }
  });

  return (
    <>
      <Line ref={lineRef} points={[anchor, labelPoint]} color="#bdeaf7" lineWidth={1.75} />
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
  );
}
