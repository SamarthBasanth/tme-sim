"use client";

import { useEffect, useRef } from "react";
import SceneCanvas from "@/components/canvas/SceneCanvas";
import MetricOverlay from "@/components/ui/MetricOverlay";
import { SECTION_COPY } from "@/lib/copy";
import { STUDIO_URL } from "@/lib/env";

/**
 * `prefers-reduced-motion: reduce` fallback. No Lenis, no ScrollTrigger, no
 * pinned/parallax camera — the model gets one small static render in the
 * hero (slow idle rotation only, see TumorModel's `reducedMotion` prop) and
 * every section simply flows down the page in normal document order. Each
 * section fades in (opacity only — no slide/blur/parallax) the first time
 * it enters the viewport, via IntersectionObserver rather than scroll math.
 */
export default function ReducedMotionExperience({ isMobile }: { isMobile: boolean }) {
  const sectionRefs = useRef<Array<HTMLElement | null>>([]);

  useEffect(() => {
    const observer = new IntersectionObserver(
      (entries) => {
        entries.forEach((entry) => {
          if (entry.isIntersecting) {
            entry.target.classList.add("is-visible");
          }
        });
      },
      { threshold: 0.2 },
    );

    sectionRefs.current.forEach((el) => el && observer.observe(el));
    return () => observer.disconnect();
  }, []);

  return (
    <div className="reduced-motion-experience">
      <div className="reduced-motion-experience__stage">
        <SceneCanvas isMobile={isMobile} mode="static" />
        <div className="reduced-motion-experience__hero-copy">
          <span className="narrative-section__kicker">
            {SECTION_COPY[0].kicker} — {SECTION_COPY[0].eyebrow}
          </span>
          <h1>{SECTION_COPY[0].headline}</h1>
          <p>{SECTION_COPY[0].body}</p>
          <a className="cta-button" href={STUDIO_URL}>
            Explore the model
          </a>
        </div>
      </div>

      {SECTION_COPY.slice(1).map((section, i) => (
        <section
          key={section.headline}
          id={`section-${i + 2}`}
          ref={(el) => {
            sectionRefs.current[i] = el;
          }}
          className={`reduced-motion-section reduced-motion-section--${section.align} fade-in-on-view`}
        >
          <span className="narrative-section__kicker">
            {section.kicker} — {section.eyebrow}
          </span>
          <h2>{section.headline}</h2>
          <p>{section.body}</p>
          {i === 2 && (
            <div className="reduced-motion-section__metric">
              <MetricOverlay initiallyHidden={false} />
            </div>
          )}
          {i === SECTION_COPY.length - 2 && (
            <a
              className="cta-button cta-button--primary"
              href="mailto:hello@example.com?subject=Demo%20request"
            >
              Request a Demo
            </a>
          )}
        </section>
      ))}
    </div>
  );
}
