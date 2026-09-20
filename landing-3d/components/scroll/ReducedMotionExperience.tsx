"use client";

import { useEffect, useRef } from "react";
import SceneCanvas from "@/components/canvas/SceneCanvas";
import MetricOverlay from "@/components/ui/MetricOverlay";
import AttributionOverlay from "@/components/ui/AttributionOverlay";
import { SECTION_COPY } from "@/lib/copy";
import SectionLedger from "@/components/ui/SectionLedger";
import SpecimenCaption from "@/components/ui/SpecimenCaption";
import { STUDIO_URL } from "@/lib/env";

/**
 * The FLOW layout. No Lenis, no ScrollTrigger, no pinned/parallax camera — the
 * model gets one small static render in the hero (slow idle rotation only, see
 * TumorModel's `reducedMotion` prop) and every section simply flows down the
 * page in normal document order. Each section fades in (opacity only — no
 * slide/blur/parallax) the first time it enters the viewport, via
 * IntersectionObserver rather than scroll math.
 *
 * REACHED BY TWO INDEPENDENT ROUTES, AND THE SECOND ONE IS NEWER THAN THE
 * NAME OF THIS FILE.
 *
 *   1. `prefers-reduced-motion: reduce` — what the visitor asked for.
 *   2. A viewport shorter than `MIN_PINNED_VIEWPORT_HEIGHT` — what the pinned
 *      layout can physically hold. See hooks/useShortViewport.ts.
 *
 * The component is unchanged by that and so are its class names: renaming
 * `.reduced-motion-*` across a 3,000-line stylesheet would be a large diff
 * that changes nothing a visitor can see, and the CSS comments that explain
 * each rule would all have to be re-pointed. What DID have to change is any
 * rule that keyed on the PREFERENCE to style this layout — see the
 * `body:has(.reduced-motion-experience) .site-header` block in globals.css,
 * which would otherwise have left prose running under the wordmark on exactly
 * the short screens route 2 exists to rescue.
 *
 * The property that makes it the right fallback for route 2 is the one it was
 * never chosen for: it is in normal document flow, so its height is unbounded
 * and nothing it contains can be clipped — which is precisely what the pinned
 * layout cannot promise at 568px.
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
          {/* Same door, same label as the cinematic hero: `/` on the
              product's origin is the workspace home, not the Studio. See the
              note where the cinematic hero's button is rendered. */}
          <a className="cta-button" href={STUDIO_URL}>
            Open the workspace
          </a>
        </div>
      </div>

      {/*
        The specimen caption belongs to the tissue above it and the cinematic
        experience has always carried it; this fallback did not, so a
        reduced-motion visitor was shown 40,887 coloured points with no legend,
        no extent and no source — which is exactly the "decorative confetti"
        reading SpecimenCaption.tsx exists to prevent, and it left the hero's
        claim that the section is REAL as an unbacked assertion for that
        visitor alone.

        It sits BELOW the stage rather than inside the hero copy on purpose:
        `.reduced-motion-experience__stage` clips at overflow:hidden, and at
        390x844 the hero copy plus this block exceeds 88svh. A caption that is
        silently cropped is worse than one placed a beat lower.
      */}
      <div className="reduced-motion-experience__specimen">
        <SpecimenCaption />
      </div>

      {SECTION_COPY.slice(1).map((section, i) => {
        /*
          KEYED ON THE BEAT, NOT ON THE SLICE INDEX, AND THAT IS A BUG FIX.
          These read `i === 2` and `i === 3`, but this list is
          SECTION_COPY.slice(1), so index i is beat i+2 — every panel landed
          one beat LATE. Measured on the rendered fallback: the mechanism
          readout ("Why they survived", 0.172 / 0.825 / 0.003) sat under beat
          04, whose subject is the macrophage spread, and the exposure
          readout (~0.43 / ~0.86 / 2.0x) sat under beat 05, the stated-limits
          spec sheet. So beat 03 asked "Which one killed the compound?" with
          no decomposition under it, and two beats were evidenced by numbers
          answering a different question.

          In the cinematic experience the coupling is explicit and correct:
          attributionPanelOpacity peaks over 2/6..3/6 (beat 03) and
          metricPanelOpacity over 3/6..4/6 (beat 04) — see lib/timeline.ts.
          The two experiences are required not to drift, and this was drift
          in the one direction that matters, which number backs which claim.

          Matching on the authored kicker removes the arithmetic that caused
          it: the panel is attached to a BEAT, not to a position in a
          derived array.
        */
        const evidence = section.ledger ? (
          <SectionLedger {...section.ledger} />
        ) : section.kicker === "03" ? (
          <div className="reduced-motion-section__metric">
            <AttributionOverlay initiallyHidden={false} />
          </div>
        ) : section.kicker === "04" ? (
          <div className="reduced-motion-section__metric">
            <MetricOverlay initiallyHidden={false} />
          </div>
        ) : null;

        return (
          <section
            key={section.headline}
            id={`section-${i + 2}`}
            ref={(el) => {
              sectionRefs.current[i] = el;
            }}
            /*
              `--paired` is set from the CONTENT, not from a width or a beat
              number: a beat gets the two-column desktop layout exactly when it
              has something to put in the second column. Beat 06 has no
              instrument, so it stays the single centred block it was composed
              as. See the `--paired` block in globals.css for why the second
              column exists at all.
            */
            className={[
              "reduced-motion-section",
              `reduced-motion-section--${section.align}`,
              evidence ? "reduced-motion-section--paired" : "",
              "fade-in-on-view",
            ]
              .filter(Boolean)
              .join(" ")}
          >
            {/*
              The argument is now a wrapper rather than three loose children,
              because at desktop widths it becomes one grid cell and the
              evidence becomes the other. Below that breakpoint this changes
              nothing: it is still a flex column with the same 1rem rhythm.
            */}
            <div className="reduced-motion-section__argument">
              <span className="narrative-section__kicker">
                {section.kicker} — {section.eyebrow}
              </span>
              <h2>{section.headline}</h2>
              <p>{section.body}</p>
              {/* Mirrors the cinematic experience's closing beat exactly: the
                  primary action is the software, email is the quiet secondary
                  route that says where the hosted instance's access key comes
                  from. The two experiences must offer the same door, and must
                  state the same condition on it. */}
              {i === SECTION_COPY.length - 2 && (
                <div className="cta-row">
                  <a className="cta-button cta-button--primary" href={STUDIO_URL}>
                    Open the workspace
                  </a>
                  <a
                    className="cta-link"
                    href="mailto:samarthbasanth@gmail.com?subject=Vivum%20access%20key"
                  >
                    or email us for an access key
                  </a>
                </div>
              )}
            </div>

            {evidence ? (
              <div className="reduced-motion-section__evidence">{evidence}</div>
            ) : null}
          </section>
        );
      })}
    </div>
  );
}
