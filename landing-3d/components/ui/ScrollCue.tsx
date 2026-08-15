"use client";

import { forwardRef } from "react";

/**
 * The one thing the hero was not doing without interaction: saying that it is
 * a scroll story.
 *
 * WHY THIS EXISTS. The entire page is a single pinned viewport scrubbed by
 * scroll — six beats, a camera path, a solved field that advances. Before this,
 * a visitor's first frame was a headline, a paragraph, a button and a specimen
 * caption over a still-ish tissue, with nothing anywhere on screen indicating
 * that any of the rest existed. The progress rail is the closest thing, and it
 * is `aria-hidden`, hidden entirely below 900px, and reads as decoration. A
 * visitor who clicks "Open the Studio" from the first frame never sees the
 * argument the page was built to make.
 *
 * WHY IT IS THIS SMALL. It is wayfinding, not a call to action — the hero
 * already has one of those, and a second competing affordance would be worse
 * than none. So: the same mono/tracked register as the kickers and the
 * specimen caption, one hairline that draws down and resets, no chevron, no
 * bounce, no circle. It is set in `--ink-faint`, below the caption in the
 * hierarchy, because it should be findable rather than loud.
 *
 * Its opacity is written by ScrollExperience from `heroCopyOpacity`, the same
 * curve the hero's own copy uses, so it cannot outlive the beat it belongs to
 * — the failure mode this page has already had twice with floating elements.
 * `pointer-events: none` throughout: it is not a control, and a visitor who
 * aims at it should hit the scene, not a dead target.
 */
const ScrollCue = forwardRef<HTMLDivElement>(function ScrollCue(_props, ref) {
  return (
    <div ref={ref} className="scroll-cue" aria-hidden="true">
      <span className="scroll-cue__label">Scroll</span>
      <span className="scroll-cue__line" />
    </div>
  );
});

export default ScrollCue;
