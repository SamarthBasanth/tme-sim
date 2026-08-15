"use client";

import { forwardRef, type ReactNode } from "react";

interface NarrativeSectionProps {
  kicker: string;
  align?: "left" | "right" | "center";
  children: ReactNode;
  className?: string;
  /**
   * First-paint opacity, before ScrollTrigger's first `onUpdate` writes a
   * scroll-derived value over it. Everything but the hero (index 0) starts
   * at 0 (off-screen/blurred-in); the hero should read correctly even in
   * the one frame before JS has run a scroll calculation.
   */
  initialOpacity?: number;
}

/**
 * One full-viewport narrative beat. Purely presentational — its
 * opacity/translate/blur are set imperatively by ScrollExperience via the
 * forwarded ref (see the comment on ScrollExperience's single ScrollTrigger
 * for why), so this component owns layout and typography only.
 */
const NarrativeSection = forwardRef<HTMLDivElement, NarrativeSectionProps>(
  function NarrativeSection(
    { kicker, align = "left", children, className = "", initialOpacity = 0 },
    ref,
  ) {
    return (
      <div
        ref={ref}
        className={`narrative-section narrative-section--${align} ${className}`}
        /*
         * `visibility` is set here, not only in the scroll loop, because the
         * scroll loop does not run until the page is scrolled.
         *
         * ScrollExperience writes `visibility: hidden` on any beat below 0.02
         * opacity, and its comment explains why: `pointer-events: none` hides a
         * section from the mouse but NOT from the tab order, so a keyboard user
         * tabs into a fully transparent beat and the focus ring appears around
         * nothing. That fix lives in ScrollTrigger's `onUpdate`, which fires on
         * a progress CHANGE — so at progress 0, on first load, it has never
         * fired, and every beat is still `visibility: visible` at opacity 0.
         *
         * Measured on the loaded page at scroll 0: Tab reached the hero's "Open
         * the Studio" (correct), then beat 06's "Open the Studio" at y = 572 and
         * its email link at y = 637 — both invisible, both on screen, both
         * ringed. The bug the loop fixes was therefore live for exactly as long
         * as a visitor had not yet scrolled, which is the entire window in which
         * a keyboard user first presses Tab.
         *
         * Deriving it from `initialOpacity` keeps the two in agreement by
         * construction, and the loop overwrites it from the first frame of
         * scroll onward.
         */
        style={{
          opacity: initialOpacity,
          visibility: initialOpacity > 0.02 ? "visible" : "hidden",
        }}
      >
        <span className="narrative-section__kicker">{kicker}</span>
        <div className="narrative-section__body">{children}</div>
      </div>
    );
  },
);

export default NarrativeSection;
