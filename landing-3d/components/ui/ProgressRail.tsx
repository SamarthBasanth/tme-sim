"use client";

import type { RefObject } from "react";

interface ProgressRailProps {
  fillRef: RefObject<HTMLDivElement | null>;
  labels: readonly string[];
}

/**
 * Minimal wayfinding rail along the viewport edge — a track with a fill bar
 * whose `scaleY` is set imperatively (by ScrollExperience's ScrollTrigger
 * `onUpdate`) directly to scroll progress. Decorative/orientation only, so
 * it's marked `aria-hidden`.
 */
export default function ProgressRail({ fillRef, labels }: ProgressRailProps) {
  return (
    <nav className="progress-rail" aria-hidden="true">
      <div className="progress-rail__track">
        <div ref={fillRef} className="progress-rail__fill" />
      </div>
      <ol className="progress-rail__ticks">
        {labels.map((label) => (
          <li key={label} className="progress-rail__tick" />
        ))}
      </ol>
    </nav>
  );
}
