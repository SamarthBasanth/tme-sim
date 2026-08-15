/**
 * Authored scroll timeline.
 *
 * Every value the 3D scene and the narrative text need is expressed here as
 * a pure function of a single global scroll progress `p` in [0, 1]. Nothing
 * in this file touches the DOM, GSAP, or React — it is just data + math, so
 * it is trivial to reason about, unit-test, and reuse identically whether
 * the user is scrolling down, scrolling up, or (for prefers-reduced-motion)
 * jumping straight to a section.
 */

import { clamp, lerp, mapRange, sampleCurve, smoothstep } from "./scrollMath";

export const SECTION_COUNT = 6;

/**
 * Section ids/labels track lib/copy.ts. They used to track an earlier
 * product -- "AI Analysis", "Treatment Planning", "Segmentation" -- long after
 * the copy above them stopped saying any of that, which is how the scene ended
 * up animating an imaging viewer under drug-transport prose.
 */
export const SECTIONS = [
  { id: "hero", label: "Hero", kicker: "01" },
  { id: "transport", label: "Transport Physics", kicker: "02" },
  { id: "mechanism", label: "Delivery or Potency", kicker: "03" },
  { id: "attribution", label: "Spatial Attribution", kicker: "04" },
  { id: "limits", label: "Stated Limits", kicker: "05" },
  { id: "cta", label: "Final", kicker: "06" },
] as const;

export type SectionId = (typeof SECTIONS)[number]["id"];

/** [start, end) progress range owned by section `index`. */
export const sectionRange = (index: number): [number, number] => {
  const seg = 1 / SECTION_COUNT;
  return [index * seg, (index + 1) * seg];
};

/** Linear (unsmoothed) local progress of `global` inside section `index`. */
export const sectionLocalProgress = (global: number, index: number): number => {
  const [start, end] = sectionRange(index);
  return mapRange(global, start, end, 0, 1, false);
};

const SECTION_INDEX: Record<SectionId, number> = {
  hero: 0,
  transport: 1,
  mechanism: 2,
  attribution: 3,
  limits: 4,
  cta: 5,
};

export const activeSectionIndex = (p: number): number =>
  clamp(Math.floor(clamp(p) * SECTION_COUNT), 0, SECTION_COUNT - 1);

// ---------------------------------------------------------------------------
// Camera path — 7 authored keyframes (one per section boundary, including a
// final resting pose) interpolated with smoothstep easing. The camera is
// therefore always an explicit, continuous function of scroll progress: no
// tweening "toward the last viewport-entered section", no discrete jumps.
// ---------------------------------------------------------------------------

export interface CameraKeyframe {
  position: readonly [number, number, number];
  lookAt: readonly [number, number, number];
  fov: number;
}

export const CAMERA_KEYFRAMES: readonly CameraKeyframe[] = [
  { position: [0, 0.05, 4.4], lookAt: [0, 0, 0], fov: 32 }, // hero
  /*
   * Beat 02 USED TO BE THE FARTHEST AND WIDEST POSE ON THE PATH — [1.7, 0.7,
   * 6.3] at fov 36 — on the one beat whose headline is "The shape is solved."
   * Measured at 1512x950 the section then spanned about 40% of the frame
   * width, and the left 56% of that frame is under the copy's directional
   * scrim, so the beat that argues the field has a solved SHAPE was showing
   * the smallest, dimmest view of that shape on the whole page. It is now the
   * closest pose of the first half instead: the section fills roughly 59% of
   * the frame, and the iso-contours (whose lobe sits in the section's
   * lower-left, where this tissue's vasculature is densest) are legible while
   * their own copy is still on screen.
   *
   * The lateral offset comes down with it. At z=6.3 an x of 1.7 was a modest
   * angle; at z=5.0 the same number would swing the section hard across frame,
   * and the pose is a truck-and-rotate (lookAt is near origin), not a pan.
   */
  { position: [0.95, 0.5, 5.0], lookAt: [0, -0.05, 0], fov: 31 }, // transport physics establishes
  { position: [-2.4, 0.3, 4.2], lookAt: [0, 0, 0], fov: 33 }, // delivery or potency
  { position: [0.35, 1.2, 3.0], lookAt: [0, 0.1, 0], fov: 29 }, // spatial attribution
  { position: [2.6, 1.0, 3.3], lookAt: [0.15, 0, 0], fov: 31 }, // stated limits
  { position: [-0.7, 0.4, 4.6], lookAt: [0.35, 0, 0], fov: 29 }, // cta enters
  { position: [0, 0.15, 5.3], lookAt: [1.05, 0, 0], fov: 27 }, // cta rests
] as const;

// ---------------------------------------------------------------------------
// Aspect correction.
//
// THE KEYFRAMES ABOVE ARE AUTHORED FOR ONE FRAME SHAPE and were catastrophic
// outside it. A three.js perspective camera's `fov` is VERTICAL, so the
// horizontal extent it sees is fov * aspect: at 1512x950 (aspect 1.59) the
// hero pose sees 4.01 world units across, and at 390x844 (aspect 0.46) the
// same numbers see 1.17 — less than half the width of the 2.6-unit section.
// Measured at 390px, four of the six beats showed no tissue at all: beat 02
// was an empty black frame with a few cells clipped to the right margin, and
// the hero — whose entire claim is "the section behind this text is a real
// one" — had the section pushed off the right edge by `heroOffsetAmount`,
// which is a desktop-sized shift applied to a phone-sized frustum.
//
// So the composition is re-derived rather than re-authored. One factor `k`,
// damped with an exponent so portrait does not overcorrect into a fisheye,
// (a) widens the vertical fov until the horizontal coverage is comparable and
// (b) divides every LATERAL offset — camera x and the section's own x slide —
// by the same amount, because a lateral offset only means anything relative
// to how wide the frame is. Vertical offsets are left alone: portrait has
// vertical room to spare, and damping them would flatten the one genuinely
// expressive move on the path (the elevated look-down at beat 04).
//
// At the reference aspect this is EXACTLY the identity (the clamp floor is 1),
// so every value verified at 1512x950 is untouched.
// ---------------------------------------------------------------------------

/** The frame shape the keyframes above were composed against (1512x950). */
export const REFERENCE_ASPECT = 1.59;

/**
 * How much narrower than the authored frame this viewport is, damped.
 * 1.0 at or above the reference aspect; ~2.16 at 390x844.
 */
export const framingFactor = (aspect: number): number =>
  Math.pow(clamp(REFERENCE_ASPECT / Math.max(aspect, 0.25), 1, 4), 0.62);

const RAD = Math.PI / 180;

/** Vertical fov, in degrees, that preserves `fov`'s coverage widened by `k`. */
const widenFov = (fov: number, k: number): number =>
  k === 1 ? fov : (2 * Math.atan(Math.tan((fov * RAD) / 2) * k)) / RAD;

export const getCameraPose = (
  progress: number,
  aspect: number = REFERENCE_ASPECT,
): CameraKeyframe => {
  const scaled = clamp(progress) * SECTION_COUNT;
  const i = Math.min(Math.floor(scaled), SECTION_COUNT - 1);
  const t = smoothstep(scaled - i);
  const a = CAMERA_KEYFRAMES[i];
  const b = CAMERA_KEYFRAMES[i + 1];
  const k = framingFactor(aspect);
  return {
    position: [
      lerp(a.position[0], b.position[0], t) / k,
      lerp(a.position[1], b.position[1], t),
      lerp(a.position[2], b.position[2], t),
    ],
    lookAt: [
      lerp(a.lookAt[0], b.lookAt[0], t) / k,
      lerp(a.lookAt[1], b.lookAt[1], t),
      lerp(a.lookAt[2], b.lookAt[2], t),
    ],
    fov: widenFov(lerp(a.fov, b.fov, t), k),
  };
};

// ---------------------------------------------------------------------------
// Scene 1 — Hero
// ---------------------------------------------------------------------------

export const particleOpacity = (p: number): number =>
  sampleCurve(p, [
    [0, 1],
    [1 / 6, 0.85],
    [2 / 6, 0.3],
    [1, 0.16],
  ]);

export const heroCopyOpacity = (p: number): number =>
  sampleCurve(p, [
    [0, 1],
    [1 / 6 - 0.03, 1],
    [1 / 6, 0],
  ]);

// ---------------------------------------------------------------------------
// Scene 2 — Transport physics
// ---------------------------------------------------------------------------

/**
 * Iso-level contour overlay (see FieldContours.tsx).
 *
 * THE CONTOURS USED TO PEAK A BEAT LATE, which is the same editorial error
 * documented (and fixed) for `fieldLabelOpacity` below — it was simply missed
 * on the overlay those labels point at. The curve read
 * `[[1/6, 0], [2/6, 0.34], [4/6, 0.26], ...]`: maximum at p = 1/3, which is
 * the BOUNDARY between beats 02 and 03, and beat 02's own copy is already
 * fully gone by then (see sectionVisual: fade = 0.32 of the segment). So the
 * solved pressure field reached full strength exactly as the paragraph
 * describing it disappeared, and then held near-maximum across the whole of
 * beat 03 and beat 04, where it belongs to nothing. Screenshotted at p=0.25
 * the contours were invisible; at p=0.42, under beat 03's headline about
 * compound structure, they were the brightest thing on screen.
 *
 * Now it rises with beat 02's copy, peaks at p ≈ 0.27 (the beat's own centre,
 * where the headline and ledger are at full opacity), and is down to a
 * residual before beat 03's copy rises at p ≈ 0.34. The residual is
 * deliberate and not zero: the field is genuinely still solved underneath the
 * later beats, and ConcentrationField — which paints the SAME viridis ramp
 * for a different quantity — takes over at 3/6, so the contours must be well
 * clear of it by then or the two read as one overlay.
 */
export const fieldContourOpacity = (p: number): number =>
  sampleCurve(p, [
    [1 / 6, 0],
    [1 / 6 + (1 / 6) * 0.38, 0.4],
    [1 / 6 + (1 / 6) * 0.72, 0.4],
    [2 / 6, 0.14],
    [4 / 6, 0.1],
    [5 / 6, 0.06],
    [1, 0.03],
  ]);

/*
 * THE CALLOUT MUST DIE WITH ITS OWN BEAT. This window used to be
 * `trapezoid(p, 1/6 + 0.02, 2/6 - 0.05, 2/6 + 0.05, 3/6 - 0.02)`, which held
 * the "Interstitial pressure" label at FULL opacity from p = 0.383 to 0.480
 * — i.e. through the middle of beat 03, half a beat after the beat it
 * describes. Beat 02's own copy is gone by p = 0.333 (see sectionVisual:
 * segment [1/6, 2/6], fade = 0.32 of the segment), so the label reached full
 * strength exactly as the text it belonged to disappeared, and was still lit
 * over beat 03's headline.
 *
 * That was not only wrong editorially, it was the cause of a measured text
 * collision at BOTH viewport sizes: the label is a camera projection of a
 * world-space anchor, and over this window it sweeps into the copy column
 * (desktop: right edge reaches x = 1051 against beat 03's body at x = 877;
 * mobile: the copy is full-bleed, so any mid-viewport label overlaps it).
 * Earlier passes missed it because it is a TRANSIENT — at each beat's settled
 * peak-opacity position the two are apart, and the overlap only exists while
 * one is rising and the other falling.
 *
 * The window is now inside beat 02: it rises with the copy and is fully out
 * by p = 1/3, before beat 03's copy begins to rise at p ≈ 0.34. Verified by
 * sweeping 61 scroll positions at 1512x950 and 390x844.
 */
export const fieldLabelOpacity = (p: number): number =>
  trapezoid(p, 1 / 6 + 0.015, 1 / 6 + 0.055, 2 / 6 - 0.045, 2 / 6 - 0.004);

function trapezoid(t: number, r0: number, r1: number, f0: number, f1: number) {
  return sampleCurve(t, [
    [r0, 0],
    [r1, 1],
    [f0, 1],
    [f1, 0],
  ]);
}

// ---------------------------------------------------------------------------
// Scene 3 — Delivery or potency
// ---------------------------------------------------------------------------

/**
 * 0 = fused mass, 1 = fully separated core / edema / margin.
 *
 * The fall-back used to run all the way to p=1, i.e. the exact end of
 * scroll. But `ctaCopyOpacity` below reaches full opacity much earlier
 * (by 5/6 + (1/6)*0.35), so for a wide stretch of the final section the
 * "Bring this into your own workflow" headline and "Request a Demo"
 * button were fully legible while the core point cloud was still sitting
 * in its separated position — offset -0.55 in Y (see ProceduralTumor.tsx's
 * `coreGroup.position.set(0, -0.55 * sep, 0)`) — low enough to visibly
 * drop across and obscure the CTA button. Finishing the fall-back well
 * before the copy is fully opaque (with margin for the per-frame damping
 * lag in ProceduralTumor/GLTFTumor) lets the mass fully recompose before
 * the final CTA becomes the visual focus.
 */
export const separationAmount = (p: number): number =>
  sampleCurve(p, [
    [2 / 6, 0],
    [3 / 6, 1],
    [5 / 6, 1],
    [5 / 6 + (1 / 6) * 0.2, 0],
  ]);

export const regionLabelOpacity = (p: number): number =>
  trapezoid(p, 2 / 6 + 0.03, 3 / 6 - 0.06, 3 / 6 + 0.1, 3 / 6 + 0.2);

// ---------------------------------------------------------------------------
// Scene 4 — Spatial attribution
// ---------------------------------------------------------------------------

/** Raw (unsmoothed) local progress through the attribution section, 0..1. */
export const frontLocalProgress = (p: number): number =>
  clamp(sectionLocalProgress(p, SECTION_INDEX.attribution));

/**
 * 0..1 — how far the drug has penetrated out of the vasculature (see
 * ConcentrationField.tsx, which turns this into a decay length).
 *
 * This used to be `scanPlaneY`, the world-space height of a plane sweeping
 * top-to-bottom through the model for the "AI Analysis" beat. The sweep is
 * kept; what advances is now concentration rather than a scanner. It stays a
 * pure function of scroll, so scrolling back up retracts the front.
 */
export const frontAdvance = (p: number): number => smoothstep(frontLocalProgress(p));

/** Drives the concentration marker and the point-cloud category highlight. */
export const exposureHighlight = (p: number): number =>
  sampleCurve(p, [
    [3 / 6 + (1 / 6) * 0.4, 0],
    [3 / 6 + (1 / 6) * 0.62, 1],
    [5 / 6, 1],
    [1, 0.12],
  ]);

export const metricPanelOpacity = (p: number): number =>
  sampleCurve(p, [
    [3 / 6, 0],
    [3 / 6 + (1 / 6) * 0.22, 1],
    [4 / 6 - (1 / 6) * 0.18, 1],
    [4 / 6, 0],
  ]);

/**
 * Section 03's mechanism-attribution readout. Same shape as the metric panel
 * one beat later, shifted to the 2/6..3/6 window, so the two instrument
 * panels enter and leave on an identical rhythm rather than each having its
 * own hand-tuned timing.
 */
export const attributionPanelOpacity = (p: number): number =>
  sampleCurve(p, [
    [2 / 6, 0],
    [2 / 6 + (1 / 6) * 0.22, 1],
    [3 / 6 - (1 / 6) * 0.18, 1],
    [3 / 6, 0],
  ]);

export const concentrationFieldOpacity = (p: number): number =>
  trapezoid(p, 3 / 6, 3 / 6 + 0.02, 4 / 6 - 0.04, 4 / 6);

// ---------------------------------------------------------------------------
// Scene 5 — Stated limits
// ---------------------------------------------------------------------------

/** Draw-on progress of the fronts leaving the vasculature (DrugFront.tsx). */
export const drugFrontReveal = (p: number): number =>
  smoothstep(clamp(sectionLocalProgress(p, SECTION_INDEX.limits)));

/**
 * Independent of the draw-on reveal above: once fully traced, the fronts
 * recede again as the final CTA settles, so the closing composition reads as
 * "the section, alone" rather than staying cluttered (the spec calls for a
 * "smooth, deliberate composition" at the very end, not an abrupt stop).
 */
export const drugFrontOpacity = (p: number): number =>
  sampleCurve(p, [
    [5 / 6, 1],
    [5 / 6 + (1 / 6) * 0.5, 0.15],
    [1, 0.15],
  ]);

/** Vessel stand-ins (Markers.tsx), revealed late in the scroll. */
export const vesselOpacity = (p: number): number =>
  sampleCurve(p, [
    [4 / 6, 0],
    [4 / 6 + (1 / 6) * 0.3, 1],
    [5 / 6, 1],
    [1, 0.2],
  ]);

// ---------------------------------------------------------------------------
// Scene 6 — Final CTA
// ---------------------------------------------------------------------------

/** 0..1 — reassembled tumor drifting to one side of frame. */
export const finalOffsetAmount = (p: number): number =>
  sampleCurve(p, [
    [5 / 6, 0],
    [1, 1],
  ]);

/**
 * HERO composition offset — the section slides RIGHT so the opening copy has
 * its own ground instead of being set on top of the densest part of the
 * specimen. Decays to 0 by the end of scene 1, after which the section is
 * centred and the scroll story proceeds exactly as before.
 *
 * Negative sign convention matches `finalOffsetAmount`'s consumer, which
 * multiplies by a positive constant: this returns a POSITIVE value for a
 * rightward shift.
 */
export const heroOffsetAmount = (p: number): number =>
  sampleCurve(p, [
    // 0.78, not 1.0: the section is 2.6 units wide and the frustum is ~2.5
    // half-widths at z=0, so a full unit of shift pushes its right edge past
    // the viewport and the crop reads as a bug rather than as a specimen
    // continuing off-frame. This leaves a margin on the right.
    [0, 0.78],
    [1 / 6, 0],
  ]);

/**
 * BEAT 02 composition offset — the section slides right so the solved
 * iso-contours are not sitting under the copy that describes them.
 *
 * The contour set is in the section's LOWER-LEFT, because that is where this
 * tissue's vasculature is densest (see the anchor comment in SceneCanvas), and
 * beat 02's copy is left-set. Measured at 1512x950 with the beat's own camera
 * pose, the contour lobe projects to about 50% of frame width — inside the
 * directional scrim, which is still at ~0.28 there. So the one beat that
 * argues the field has a solved SHAPE was drawing that shape behind its own
 * text scrim. Half a unit to the right clears it (the scrim reaches 0 at 56%)
 * without pushing the section's right edge out of frame.
 *
 * Returns to 0 by 2/6: beat 03 composes the other way round — copy right,
 * tissue left — and a rightward shift held into it would push the tissue
 * toward the text rather than away from it.
 */
export const transportOffsetAmount = (p: number): number =>
  sampleCurve(p, [
    [1 / 6, 0],
    [1 / 6 + (1 / 6) * 0.35, 0.44],
    [2 / 6 - (1 / 6) * 0.18, 0.44],
    [2 / 6, 0],
  ]);

/**
 * The section's total x offset, in the units its consumers multiply by 1.15.
 *
 * Exists so the tissue and everything positioned in ITS coordinates
 * (TransportOverlay) read ONE function rather than each summing the terms
 * themselves — that is exactly how the overlays previously drifted out of
 * lockstep with the model. Add any future offset term here, not at a call
 * site.
 *
 * `aspect` divides the total by the same `framingFactor` the camera path uses
 * (see the aspect-correction block above). Both terms here are lateral
 * COMPOSITION offsets — "slide the section out from under the copy", "slide it
 * aside for the CTA" — expressed in world units against a wide frame. Fed
 * unchanged into a portrait frustum a third that width, 0.78 stopped meaning
 * "clear of the text" and started meaning "off the right edge of the screen",
 * which is exactly what it measured at 390px. Identity at the reference
 * aspect.
 */
export const sectionOffsetAmount = (
  p: number,
  aspect: number = REFERENCE_ASPECT,
): number =>
  (finalOffsetAmount(p) + heroOffsetAmount(p) + transportOffsetAmount(p)) /
  framingFactor(aspect);

export const ctaCopyOpacity = (p: number): number =>
  sampleCurve(p, [
    [5 / 6, 0],
    [5 / 6 + (1 / 6) * 0.35, 1],
    [1, 1],
  ]);

// ---------------------------------------------------------------------------
// Narrative text choreography — shared by every NarrativeSection instance.
// ---------------------------------------------------------------------------

export interface SectionVisual {
  opacity: number;
  translateY: number;
  blurPx: number;
}

/**
 * Fade / slide / blur curve for narrative text, index-aware so the hero
 * (index 0) starts fully visible instead of fading in from nothing, and the
 * final CTA (last index) stays fully visible once reached instead of fading
 * back out at the very bottom of the page.
 */
export const sectionVisual = (p: number, index: number): SectionVisual => {
  const [start, end] = sectionRange(index);
  const segLen = end - start;
  const fade = segLen * 0.32;

  let opacity: number;
  if (index === 0) {
    opacity = sampleCurve(p, [
      [start, 1],
      [end - fade, 1],
      [end, 0],
    ]);
  } else if (index === SECTION_COUNT - 1) {
    opacity = sampleCurve(p, [
      [start, 0],
      [start + fade, 1],
      [1, 1],
    ]);
  } else {
    opacity = sampleCurve(p, [
      [start, 0],
      [start + fade, 1],
      [end - fade, 1],
      [end, 0],
    ]);
  }

  const center = (start + end) / 2;
  const direction = p < center ? 1 : -1;
  return {
    opacity,
    translateY: (1 - opacity) * 26 * direction,
    blurPx: (1 - opacity) * 6,
  };
};

export { SECTION_INDEX };
