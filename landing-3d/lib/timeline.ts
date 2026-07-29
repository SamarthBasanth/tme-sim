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

export const SECTIONS = [
  { id: "hero", label: "Hero", kicker: "01" },
  { id: "anatomy", label: "Anatomical Context", kicker: "02" },
  { id: "segmentation", label: "Segmentation", kicker: "03" },
  { id: "analysis", label: "AI Analysis", kicker: "04" },
  { id: "planning", label: "Treatment Planning", kicker: "05" },
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
  anatomy: 1,
  segmentation: 2,
  analysis: 3,
  planning: 4,
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
  { position: [1.7, 0.7, 6.3], lookAt: [0, -0.05, 0], fov: 36 }, // anatomy establishes
  { position: [-2.4, 0.3, 4.2], lookAt: [0, 0, 0], fov: 33 }, // segmentation
  { position: [0.35, 1.2, 3.0], lookAt: [0, 0.1, 0], fov: 29 }, // analysis
  { position: [2.6, 1.0, 3.3], lookAt: [0.15, 0, 0], fov: 31 }, // planning
  { position: [-0.7, 0.4, 4.6], lookAt: [0.35, 0, 0], fov: 29 }, // cta enters
  { position: [0, 0.15, 5.3], lookAt: [1.05, 0, 0], fov: 27 }, // cta rests
] as const;

export const getCameraPose = (progress: number): CameraKeyframe => {
  const scaled = clamp(progress) * SECTION_COUNT;
  const i = Math.min(Math.floor(scaled), SECTION_COUNT - 1);
  const t = smoothstep(scaled - i);
  const a = CAMERA_KEYFRAMES[i];
  const b = CAMERA_KEYFRAMES[i + 1];
  return {
    position: [
      lerp(a.position[0], b.position[0], t),
      lerp(a.position[1], b.position[1], t),
      lerp(a.position[2], b.position[2], t),
    ],
    lookAt: [
      lerp(a.lookAt[0], b.lookAt[0], t),
      lerp(a.lookAt[1], b.lookAt[1], t),
      lerp(a.lookAt[2], b.lookAt[2], t),
    ],
    fov: lerp(a.fov, b.fov, t),
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
// Scene 2 — Anatomical context
// ---------------------------------------------------------------------------

export const anatomyOpacity = (p: number): number =>
  sampleCurve(p, [
    [1 / 6, 0],
    [2 / 6, 0.34],
    [4 / 6, 0.26],
    [5 / 6, 0.1],
    [1, 0.04],
  ]);

export const anatomyLabelOpacity = (p: number): number =>
  trapezoid(p, 1 / 6 + 0.02, 2 / 6 - 0.05, 2 / 6 + 0.05, 3 / 6 - 0.02);

function trapezoid(t: number, r0: number, r1: number, f0: number, f1: number) {
  return sampleCurve(t, [
    [r0, 0],
    [r1, 1],
    [f0, 1],
    [f1, 0],
  ]);
}

// ---------------------------------------------------------------------------
// Scene 3 — Segmentation
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

export const segmentationLabelOpacity = (p: number): number =>
  trapezoid(p, 2 / 6 + 0.03, 3 / 6 - 0.06, 3 / 6 + 0.1, 3 / 6 + 0.2);

// ---------------------------------------------------------------------------
// Scene 4 — AI analysis
// ---------------------------------------------------------------------------

/** Raw (unsmoothed) local progress through the analysis section, 0..1. */
export const scanLocalProgress = (p: number): number =>
  clamp(sectionLocalProgress(p, SECTION_INDEX.analysis));

/** World-space Y the scan plane sweeps across, top to bottom. */
export const scanPlaneY = (p: number): number =>
  lerp(1.35, -1.35, smoothstep(scanLocalProgress(p)));

export const highlightOpacity = (p: number): number =>
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

export const scanPlaneOpacity = (p: number): number =>
  trapezoid(p, 3 / 6, 3 / 6 + 0.02, 4 / 6 - 0.04, 4 / 6);

// ---------------------------------------------------------------------------
// Scene 5 — Treatment planning
// ---------------------------------------------------------------------------

export const trajectoryReveal = (p: number): number =>
  smoothstep(clamp(sectionLocalProgress(p, SECTION_INDEX.planning)));

/**
 * Independent of the draw-on reveal above: once fully traced, the path
 * recedes again as the final CTA settles, so the closing composition reads
 * as "tumor, reassembled" rather than staying cluttered with planning
 * artifacts (the spec calls for a "smooth, deliberate composition" at the
 * very end, not an abrupt stop).
 */
export const trajectoryOpacity = (p: number): number =>
  sampleCurve(p, [
    [5 / 6, 1],
    [5 / 6 + (1 / 6) * 0.5, 0.15],
    [1, 0.15],
  ]);

export const criticalStructureOpacity = (p: number): number =>
  sampleCurve(p, [
    [4 / 6, 0],
    [4 / 6 + (1 / 6) * 0.3, 1],
    [5 / 6, 1],
    [1, 0.2],
  ]);

export const planningCopyOpacity = (p: number): number =>
  trapezoid(p, 4 / 6 + 0.02, 4 / 6 + 0.12, 5 / 6 - 0.08, 5 / 6 - 0.01);

// ---------------------------------------------------------------------------
// Scene 6 — Final CTA
// ---------------------------------------------------------------------------

/** 0..1 — reassembled tumor drifting to one side of frame. */
export const finalOffsetAmount = (p: number): number =>
  sampleCurve(p, [
    [5 / 6, 0],
    [1, 1],
  ]);

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
