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
        style={{ opacity: initialOpacity }}
      >
        <span className="narrative-section__kicker">{kicker}</span>
        <div className="narrative-section__body">{children}</div>
      </div>
    );
  },
);

export default NarrativeSection;
