"use client";

import type { RefObject } from "react";

interface ProgressRailProps {
  fillRef: RefObject<HTMLDivElement | null>;
  /**
   * The tick list. ScrollExperience writes `data-active` on it every scroll
   * frame (see its ScrollTrigger `onUpdate`) with the index of the beat the
   * viewer is currently in.
   */
  ticksRef: RefObject<HTMLOListElement | null>;
  /**
   * How many beats. This used to be `labels: readonly string[]`, and the
   * labels were never rendered — they were passed in from SECTION_COPY and
   * used only as React keys, so the rail advertised an API it did not have.
   * They are deliberately still not rendered: measured at 1512x950 the rail
   * sits at 98% of frame width and beats 03 and 05 set their copy to 94%, so
   * any label growing leftward out of the rail lands in running text.
   */
  count: number;
}

/**
 * Minimal wayfinding rail along the viewport edge — a track with a fill bar
 * whose `scaleY` is set imperatively (by ScrollExperience's ScrollTrigger
 * `onUpdate`) directly to scroll progress. Decorative/orientation only, so
 * it's marked `aria-hidden`.
 *
 * The ticks used to be six identical hollow dots that never changed, which
 * made the rail a bar chart of one number: it told you how far down the page
 * you were and nothing about WHERE, even though the page is explicitly six
 * discrete beats and every other element on screen knows which one is
 * current. The active tick now fills — one attribute written per frame, no
 * React state in the scroll loop, consistent with how everything else here is
 * driven.
 */
export default function ProgressRail({ fillRef, ticksRef, count }: ProgressRailProps) {
  return (
    <nav className="progress-rail" aria-hidden="true">
      <div className="progress-rail__track">
        <div ref={fillRef} className="progress-rail__fill" />
      </div>
      <ol ref={ticksRef} className="progress-rail__ticks" data-active="0">
        {Array.from({ length: count }, (_, i) => (
          <li key={i} className="progress-rail__tick" />
        ))}
      </ol>
    </nav>
  );
}
