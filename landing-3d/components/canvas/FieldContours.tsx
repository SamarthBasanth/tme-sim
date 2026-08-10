"use client";

import { useMemo, useRef, type ReactNode } from "react";
import { useFrame } from "@react-three/fiber";
import { Line } from "@react-three/drei";
import type { Line2, LineSegments2 } from "three-stdlib";
import * as THREE from "three";
import { useScrollProgressRef } from "@/components/scroll/ScrollProgressContext";
import { fieldContourOpacity } from "@/lib/timeline";
import { damp } from "@/lib/scrollMath";
import { IFP_LEVELS } from "@/lib/ifpContours";

/**
 * Iso-contours of the SOLVED interstitial fluid pressure field over the
 * measured section.
 *
 * This used to be four procedurally-generated rings — an ellipse per level
 * with its radius walked by fBm noise — carrying the annotation "Illustrative,
 * not solved". They were the right SHAPE (nested closed contours in the plane
 * of the section are what a scalar field over a 2D section looks like) but
 * they were furniture, and they were the most prominent labelled element in
 * the one beat of the page that claims the simulator actually solves
 * something. The disclaimer was honest; deleting it would have been the
 * dishonest fix. So the picture was made real instead.
 *
 * Every vertex now comes from `tmesim.physics.pressure.solve_ifp` run on the
 * same reconstruction the hero point cloud is exported from — see
 * `lib/ifpContours.ts` for what is solved, what the one display choice is
 * (a stated 250 um local mean, forced by the fact that the raw field's
 * per-vessel structure is smaller than a screen pixel at this zoom), and why
 * the contours register with the cells rather than merely sitting near them.
 *
 * The peak is no longer placed by hand. The old version offset its rings by
 * `PEAK = { x: -0.18, y: -0.06 }` with a comment explaining that a field
 * maximum is a property of the tissue and not of the crop — which was true,
 * and is now actually the case: the lobe sits in the lower-left of the section
 * because that is where this tissue's vasculature is densest, on top of the
 * endothelial-density maximum `shaders/viridis.ts` derived independently from
 * the GLB as `ENTRY_POINTS[0]`.
 */

/**
 * Just in front of the section plane. Behind it, the 40k-point cloud swallowed
 * the lines almost entirely and the annotation pointed at nothing; a field
 * overlay is drawn ON a section anyway, which is where the Studio draws it.
 */
const CONTOUR_Z = 0.04;

/**
 * Peak of the authored scroll envelope in `lib/timeline.ts`, sampled rather
 * than hard-coded so the two cannot drift apart if that curve is retimed.
 *
 * The envelope is renormalised by this below. Its absolute peak (0.34) was
 * authored for a single light cyan line; these lines are coloured by the
 * viridis ramp, whose low end is dark by construction — at 0.34 the two lower
 * iso-levels are simply not on screen, which would hide half the field. The
 * curve's SHAPE — when it rises, holds and falls across the six sections — is
 * what the timeline is actually authoring, and renormalising preserves it
 * exactly while letting the ramp's dark end survive on a `#07080a` background.
 */
const ENVELOPE_PEAK = (() => {
  let peak = 0;
  for (let i = 0; i <= 240; i++) peak = Math.max(peak, fieldContourOpacity(i / 240));
  return peak > 0 ? peak : 1;
})();

/** Opacity of a contour at the envelope's peak. */
const MAX_OPACITY = 1;

/**
 * Wider than the 1.1 the procedural rings used.
 *
 * Two reasons, both consequences of the lines being real. There are far fewer
 * of them (four solved levels over one lobe, not four rings spanning the whole
 * section), so there is room. And the timeline's envelope only reaches its
 * peak at the very end of section 02 — at that section's settled centre the
 * contours are at roughly half of it, which the ramp's dark end does not
 * survive at hairline width. Width is the one legibility lever here that costs
 * no honesty: it does not move a vertex or change a level's colour.
 */
const LINE_WIDTH = 2.4;

/**
 * Number of iso-levels dropped on mobile, from the LOW end.
 *
 * The old version dropped the innermost level, because its contours were
 * 84-segment noise walks and the count was a vertex budget. These are ~30
 * vertices per level — the whole solved set is about a third of one old ring —
 * so the budget no longer binds and this is purely about legibility on a
 * narrow viewport. It drops the OUTERMOST level instead: that is the lowest
 * pressure, therefore the darkest point on the ramp and the least legible, and
 * it is also the widest, so it is the one that runs into the section's edge
 * first at 390px.
 */
const LOW_DETAIL_DROPPED_LEVELS = 1;

/** Matches GLTFTumor's planar sway, so real contours track the real cells. */
const SWAY_RATE = 0.22;
const SWAY_AMPLITUDE = 0.16;
const SWAY_TILT = 0.05;
const SCROLL_ROTATION = 0.35;

interface FieldContoursProps {
  lowDetail?: boolean;
  /**
   * Rendered INSIDE the same swaying group as the contours, not as a sibling
   * in the parent scene. This group sways every frame, so any annotation
   * meant to sit on a specific contour point must be a child and let
   * Three.js's normal parent-child transform composition keep it correct,
   * instead of a hand-computed world coordinate that is only ever right at
   * the one instant it was computed for. See SceneCanvas.tsx for the anchor
   * point itself (in this group's LOCAL space).
   */
  children?: ReactNode;
}

export default function FieldContours({ lowDetail = false, children }: FieldContoursProps) {
  const progressRef = useScrollProgressRef();
  const groupRef = useRef<THREE.Group>(null);
  const lineRefs = useRef<Array<Line2 | LineSegments2 | null>>([]);
  const opacitySmoothed = useRef(0);
  const scrollRotation = useRef(0);

  /**
   * Flattened to one entry per PATH, not per level: an iso-level of a real
   * field is not guaranteed to be a single loop (the innermost level of this
   * section resolves into two separate closed loops around two local maxima,
   * which the old one-ellipse-per-level model could not express at all).
   */
  const lines = useMemo(() => {
    const levels = lowDetail ? IFP_LEVELS.slice(LOW_DETAIL_DROPPED_LEVELS) : IFP_LEVELS;
    return levels.flatMap((level) =>
      level.paths.map((points) => ({ points, color: level.colorHex })),
    );
  }, [lowDetail]);

  useFrame((state, delta) => {
    const progress = progressRef.current;
    const target = fieldContourOpacity(progress);
    opacitySmoothed.current = damp(opacitySmoothed.current, target, 3, delta);
    const o = Math.min(1, (opacitySmoothed.current / ENVELOPE_PEAK) * MAX_OPACITY);
    lineRefs.current.forEach((line) => {
      if (!line) return;
      const material = line.material as THREE.Material & { opacity: number; transparent: boolean };
      material.opacity = o;
      material.transparent = true;
    });
    if (groupRef.current) {
      groupRef.current.visible = opacitySmoothed.current > 0.005;
      // The same shallow sway law GLTFTumor gives a planar asset (a flat thing
      // must never turn edge-on or it vanishes), so the contours stay
      // registered with the section they annotate instead of spinning free of
      // it. This used to be a loose approximation with the note that "a few
      // degrees of drift costs nothing" because the contours were
      // illustrative. They are a solved field now and registration is the
      // claim, so the law is mirrored properly: same rate, same amplitude,
      // same tilt, and — the term that was missing and is worth the most,
      // being a pure function of scroll rather than of a clock — the same
      // scroll-driven yaw, damped with the same constant.
      //
      // What still cannot be shared is the PHASE of the idle sway: GLTFTumor
      // accumulates its own angle from the moment its GLB resolves, so it has
      // no common origin with this clock (the same reason TransportOverlay
      // mirrors only the translation). That leaves a few degrees of parallax,
      // which is far smaller than the feature the contours enclose.
      scrollRotation.current = damp(scrollRotation.current, progress * SCROLL_ROTATION, 2.4, delta);
      const sway = Math.sin(state.clock.elapsedTime * SWAY_RATE) * SWAY_AMPLITUDE;
      groupRef.current.rotation.y = sway + scrollRotation.current;
      groupRef.current.rotation.x = Math.cos(state.clock.elapsedTime * SWAY_RATE * 0.6) * SWAY_TILT;
    }
  });

  return (
    <group ref={groupRef} position={[0, 0, CONTOUR_Z]}>
      {lines.map((line, i) => (
        <Line
          key={i}
          ref={(el) => {
            lineRefs.current[i] = el;
          }}
          points={line.points}
          /* Colour encodes the iso-level: viridis(level / field max) under the
             product's own ramp table, baked by the exporter. IFP is not drug
             concentration, which ConcentrationField paints with the same ramp
             — the annotation in SceneCanvas.tsx names the quantity and its
             units so the two are never read as one scale. */
          color={line.color}
          lineWidth={LINE_WIDTH}
          transparent
          opacity={0}
        />
      ))}
      {children}
    </group>
  );
}
