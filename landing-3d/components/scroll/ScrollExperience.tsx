"use client";

import { useEffect, useRef } from "react";
import gsap from "gsap";
import { ScrollTrigger } from "gsap/ScrollTrigger";
import type Lenis from "lenis";
import SceneCanvas from "@/components/canvas/SceneCanvas";
import NarrativeSection from "@/components/ui/NarrativeSection";
import MetricOverlay from "@/components/ui/MetricOverlay";
import AttributionOverlay from "@/components/ui/AttributionOverlay";
import ProgressRail from "@/components/ui/ProgressRail";
import ReducedMotionExperience from "./ReducedMotionExperience";
import { ScrollProgressProvider, useScrollProgressRef } from "./ScrollProgressContext";
import { useLenis } from "@/hooks/useLenis";
import { useIsMobile } from "@/hooks/useIsMobile";
import { usePrefersReducedMotion } from "@/hooks/usePrefersReducedMotion";
import {
  SECTION_COUNT,
  activeSectionIndex,
  attributionPanelOpacity,
  heroCopyOpacity,
  metricPanelOpacity,
  sectionVisual,
} from "@/lib/timeline";
import { SECTION_COPY } from "@/lib/copy";
import { STUDIO_URL } from "@/lib/env";
import SpecimenCaption from "@/components/ui/SpecimenCaption";
import SectionLedger from "@/components/ui/SectionLedger";
import ScrollCue from "@/components/ui/ScrollCue";

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
  const attributionElRef = useRef<HTMLDivElement>(null);
  const railFillRef = useRef<HTMLDivElement>(null);
  const railTicksRef = useRef<HTMLOListElement>(null);
  const scrollCueRef = useRef<HTMLDivElement>(null);
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
            /*
             * `pointer-events: none` hides a section from the mouse but NOT
             * from the tab order, so tabbing through this page walked into
             * "Request a Demo" while beat 06 was fully transparent — the focus
             * ring added below then appeared around nothing. `visibility`
             * removes it from the tab order and the accessibility tree.
             *
             * The threshold is 0.02, not the 0.4 used for pointer-events: at
             * 0.4 a section would blink out while it was still plainly on
             * screen. At 0.02 the section is already indistinguishable from
             * absent, so this only ever fires where the eye agrees.
             */
            el.style.visibility = opacity > 0.02 ? "visible" : "hidden";
          });

          // Both instrument panels are driven the same way and with the same
          // easing, so they read as one recurring device appearing twice
          // rather than two unrelated cards.
          for (const [el, opacity] of [
            [metricElRef.current, metricPanelOpacity(p)],
            [attributionElRef.current, attributionPanelOpacity(p)],
          ] as const) {
            if (!el) continue;
            el.style.opacity = String(opacity);
            el.style.transform = `translate3d(0, ${(1 - opacity) * 14}px, 0)`;
            el.style.pointerEvents = opacity > 0.4 ? "auto" : "none";
          }

          // Dies with the hero copy it belongs to, on that copy's own curve.
          if (scrollCueRef.current) {
            const cue = heroCopyOpacity(p);
            scrollCueRef.current.style.opacity = String(cue);
            scrollCueRef.current.style.visibility = cue > 0.02 ? "visible" : "hidden";
          }

          if (railFillRef.current) {
            railFillRef.current.style.transform = `scaleY(${p})`;
          }

          // Which beat the rail should mark as current. Written as an
          // attribute rather than a class toggle over children so the DOM
          // write is a single string compare per frame; the CSS below selects
          // the matching tick.
          if (railTicksRef.current) {
            const active = String(activeSectionIndex(p));
            if (railTicksRef.current.dataset.active !== active) {
              railTicksRef.current.dataset.active = active;
            }
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

        {/*
          ONE MAPPING OVER SECTION_COPY, NOT SIX HAND-WRITTEN BLOCKS.

          This was six near-identical copies differing only in a literal index,
          which meant the beat->element wiring was retyped six times and the
          index appeared four times per beat (kicker, align, headline, body,
          ledger) plus once more in the ref callback. The reduced-motion
          fallback has already been bitten once by exactly this class of
          mistake — see the "KEYED ON THE BEAT, NOT ON THE SLICE INDEX" note in
          ReducedMotionExperience.tsx, where an off-by-one put the wrong
          instrument under two beats. Deriving every index from one loop
          variable removes the opportunity.

          The two beats that are genuinely different stay explicitly different:
          the hero (index 0) carries the h1, the specimen caption and
          `initialOpacity={1}`, and the closing beat carries the action pair.
          Both are keyed off the position rather than off a class name, so the
          rendered output is what it was before, element for element.
        */}
        <div className="narrative-layer">
          {SECTION_COPY.map((section, i) => {
            const isHero = i === 0;
            const isClosing = i === SECTION_COPY.length - 1;
            const Heading = isHero ? "h1" : "h2";
            return (
              <NarrativeSection
                key={section.kicker}
                ref={(el) => {
                  sectionElRefs.current[i] = el;
                }}
                kicker={`${section.kicker} — ${section.eyebrow}`}
                align={section.align}
                className={
                  isHero
                    ? "narrative-section--hero"
                    : isClosing
                      ? "narrative-section--cta"
                      : ""
                }
                initialOpacity={isHero ? 1 : 0}
              >
                <Heading>{section.headline}</Heading>
                <p>{section.body}</p>
                {section.ledger ? <SectionLedger {...section.ledger} /> : null}
                {isHero ? (
                  <>
                    <a className="cta-button" href={STUDIO_URL}>
                      Open the Studio
                    </a>
                    {/* Turns the tissue behind this text from decoration into
                        data: without a legend the measured cell types read as
                        confetti and "the section is real" lands as an
                        assertion. */}
                    <SpecimenCaption />
                  </>
                ) : null}
                {isClosing ? (
                  /*
                   * The page ends by handing over the software, not by asking
                   * for an introduction. This was `mailto:` — the closing beat
                   * of a page whose whole argument is "run it and check the
                   * numbers yourself" terminated in a contact form, which puts
                   * a human gate in front of the one thing a visitor can
                   * verify unaided. Same STUDIO_URL the hero opens, so the
                   * first and last frames offer the same door.
                   *
                   * Email survives as the secondary route, deliberately
                   * quieter: someone who wants to talk rather than run software
                   * still has a way through, and it is a text link rather than
                   * a second pill so it cannot compete with the primary action.
                   */
                  <div className="cta-row">
                    <a className="cta-button cta-button--primary" href={STUDIO_URL}>
                      Open the Studio
                    </a>
                    <a
                      className="cta-link"
                      href="mailto:samarthbasanth@gmail.com?subject=Codon%20enquiry"
                    >
                      or email us about your tissue
                    </a>
                  </div>
                ) : null}
              </NarrativeSection>
            );
          })}
        </div>

        {/* Section 03's copy is set on the right, so its readout sits left —
            mirroring the metric panel one beat later, which sits right
            opposite left-set copy. */}
        <ScrollCue ref={scrollCueRef} />
        <AttributionOverlay ref={attributionElRef} />
        <MetricOverlay ref={metricElRef} />
        <ProgressRail fillRef={railFillRef} ticksRef={railTicksRef} count={SECTION_COUNT} />
      </div>
    </div>
  );
}
