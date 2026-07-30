"use client";

import { useEffect, useRef } from "react";
import gsap from "gsap";
import { ScrollTrigger } from "gsap/ScrollTrigger";
import type Lenis from "lenis";
import SceneCanvas from "@/components/canvas/SceneCanvas";
import NarrativeSection from "@/components/ui/NarrativeSection";
import MetricOverlay from "@/components/ui/MetricOverlay";
import ProgressRail from "@/components/ui/ProgressRail";
import ReducedMotionExperience from "./ReducedMotionExperience";
import { ScrollProgressProvider, useScrollProgressRef } from "./ScrollProgressContext";
import { useLenis } from "@/hooks/useLenis";
import { useIsMobile } from "@/hooks/useIsMobile";
import { usePrefersReducedMotion } from "@/hooks/usePrefersReducedMotion";
import { SECTION_COUNT, metricPanelOpacity, sectionVisual } from "@/lib/timeline";
import { SECTION_COPY } from "@/lib/copy";
import { STUDIO_URL } from "@/lib/env";

if (typeof window !== "undefined") {
  gsap.registerPlugin(ScrollTrigger);
}

export default function ScrollExperience() {
  return (
    <ScrollProgressProvider>
      <Experience />
    </ScrollProgressProvider>
  );
}

function Experience() {
  const progressRef = useScrollProgressRef();
  const reducedMotion = usePrefersReducedMotion();
  const isMobile = useIsMobile();

  const wrapperRef = useRef<HTMLDivElement>(null);
  const sectionElRefs = useRef<Array<HTMLDivElement | null>>([]);
  const metricElRef = useRef<HTMLDivElement>(null);
  const railFillRef = useRef<HTMLDivElement>(null);
  const lenisRef = useRef<Lenis | null>(null);

  useLenis(reducedMotion, (lenis) => {
    lenisRef.current = lenis;
  });

  // ---------------------------------------------------------------------
  // Why one ScrollTrigger + imperative DOM writes, instead of Drei's
  // ScrollControls or a swarm of per-section GSAP tweens:
  //
  // ScrollControls renders its own virtual scroll container inside the
  // R3F canvas tree, which fights with mixing real DOM narrative sections
  // (with their own text layout/accessibility needs) into the same scroll
  // space, and it doesn't have an obvious hook into Lenis. A single
  // ScrollTrigger spanning the whole pinned experience gives one
  // authoritative `progress` number per frame, tied 1:1 to scroll offset
  // (`scrub: true`, no discrete "enters viewport" triggers per section —
  // required by the brief). That number drives *everything*: the 3D scene
  // reads it from a ref inside useFrame (see ScrollProgressContext), and
  // narrative text/metrics are updated here via direct `el.style` writes
  // rather than `gsap.to()` tweens or React state — both would re-run on
  // every scroll frame (60+ times/sec) for no benefit, since every visual
  // value is already a pure, precomputed function of `progress`
  // (lib/timeline.ts) and R3F applies its own per-frame easing/damping on
  // top. Because nothing here is stateful (no "has this fired yet"), the
  // whole timeline is trivially correct scrolling in either direction.
  // ---------------------------------------------------------------------
  useEffect(() => {
    if (reducedMotion || !wrapperRef.current) return;

    const ctx = gsap.context(() => {
      const trigger = ScrollTrigger.create({
        trigger: wrapperRef.current,
        start: "top top",
        end: "bottom bottom",
        scrub: true,
        onUpdate: (self) => {
          const p = self.progress;
          progressRef.current = p;

          sectionElRefs.current.forEach((el, i) => {
            if (!el) return;
            const { opacity, translateY, blurPx } = sectionVisual(p, i);
            el.style.opacity = String(opacity);
            el.style.transform = `translate3d(0, ${translateY}px, 0)`;
            el.style.filter = blurPx > 0.05 ? `blur(${blurPx}px)` : "";
            el.style.pointerEvents = opacity > 0.4 ? "auto" : "none";
          });

          if (metricElRef.current) {
            const o = metricPanelOpacity(p);
            metricElRef.current.style.opacity = String(o);
            metricElRef.current.style.transform = `translate3d(0, ${(1 - o) * 14}px, 0)`;
            metricElRef.current.style.pointerEvents = o > 0.4 ? "auto" : "none";
          }

          if (railFillRef.current) {
            railFillRef.current.style.transform = `scaleY(${p})`;
          }
        },
      });

      // Fonts / layout can settle a frame late; one refresh keeps the
      // trigger's measured start/end honest without polling.
      requestAnimationFrame(() => ScrollTrigger.refresh());

      return () => trigger.kill();
    }, wrapperRef);

    // Cleans up the ScrollTrigger instance (and anything else created
    // inside the context) on unmount or when reducedMotion flips true.
    return () => ctx.revert();
  }, [reducedMotion, progressRef]);

  if (reducedMotion) {
    return <ReducedMotionExperience isMobile={isMobile} />;
  }

  return (
    <div
      ref={wrapperRef}
      className="scroll-space"
      style={{ height: `${SECTION_COUNT * 100}vh` }}
    >
      <div className="scroll-pin">
        <SceneCanvas isMobile={isMobile} mode="scroll" />

        <div className="narrative-layer">
          <NarrativeSection
            ref={(el) => {
              sectionElRefs.current[0] = el;
            }}
            kicker={`${SECTION_COPY[0].kicker} — ${SECTION_COPY[0].eyebrow}`}
            align={SECTION_COPY[0].align}
            className="narrative-section--hero"
            initialOpacity={1}
          >
            <h1>{SECTION_COPY[0].headline}</h1>
            <p>{SECTION_COPY[0].body}</p>
            <a className="cta-button" href={STUDIO_URL}>
              Explore the model
            </a>
          </NarrativeSection>

          <NarrativeSection
            ref={(el) => {
              sectionElRefs.current[1] = el;
            }}
            kicker={`${SECTION_COPY[1].kicker} — ${SECTION_COPY[1].eyebrow}`}
            align={SECTION_COPY[1].align}
          >
            <h2>{SECTION_COPY[1].headline}</h2>
            <p>{SECTION_COPY[1].body}</p>
          </NarrativeSection>

          <NarrativeSection
            ref={(el) => {
              sectionElRefs.current[2] = el;
            }}
            kicker={`${SECTION_COPY[2].kicker} — ${SECTION_COPY[2].eyebrow}`}
            align={SECTION_COPY[2].align}
          >
            <h2>{SECTION_COPY[2].headline}</h2>
            <p>{SECTION_COPY[2].body}</p>
          </NarrativeSection>

          <NarrativeSection
            ref={(el) => {
              sectionElRefs.current[3] = el;
            }}
            kicker={`${SECTION_COPY[3].kicker} — ${SECTION_COPY[3].eyebrow}`}
            align={SECTION_COPY[3].align}
          >
            <h2>{SECTION_COPY[3].headline}</h2>
            <p>{SECTION_COPY[3].body}</p>
          </NarrativeSection>

          <NarrativeSection
            ref={(el) => {
              sectionElRefs.current[4] = el;
            }}
            kicker={`${SECTION_COPY[4].kicker} — ${SECTION_COPY[4].eyebrow}`}
            align={SECTION_COPY[4].align}
          >
            <h2>{SECTION_COPY[4].headline}</h2>
            <p>{SECTION_COPY[4].body}</p>
          </NarrativeSection>

          <NarrativeSection
            ref={(el) => {
              sectionElRefs.current[5] = el;
            }}
            kicker={`${SECTION_COPY[5].kicker} — ${SECTION_COPY[5].eyebrow}`}
            align={SECTION_COPY[5].align}
            className="narrative-section--cta"
          >
            <h2>{SECTION_COPY[5].headline}</h2>
            <p>{SECTION_COPY[5].body}</p>
            <a className="cta-button cta-button--primary" href="mailto:hello@example.com?subject=Demo%20request">
              Request a Demo
            </a>
          </NarrativeSection>
        </div>

        <MetricOverlay ref={metricElRef} />
        <ProgressRail fillRef={railFillRef} labels={SECTION_COPY.map((s) => s.eyebrow)} />
      </div>
    </div>
  );
}
