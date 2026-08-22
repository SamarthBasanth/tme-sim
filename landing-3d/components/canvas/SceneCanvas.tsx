"use client";

import { Suspense, useCallback, useEffect, useState } from "react";
import { Canvas, useThree } from "@react-three/fiber";
import { Preload } from "@react-three/drei";
import * as THREE from "three";
import CameraRig from "./CameraRig";
import SceneErrorBoundary from "./SceneErrorBoundary";
import { useModelAvailability } from "@/hooks/useModelAvailability";
import { useWebGLSupport } from "@/hooks/useWebGLSupport";
import SceneLighting from "./SceneLighting";
import TumorModel from "./TumorModel";
import FieldContours from "./FieldContours";
import Particles from "./Particles";
import ConcentrationField from "./ConcentrationField";
import DrugFront from "./DrugFront";
import TransportOverlay from "./TransportOverlay";
import Annotation from "./Annotation";
import { ConcentrationSource, VesselNetwork } from "./Markers";
import Loader from "./Loader";
import { ENTRY_POINTS } from "./shaders/viridis";
import {
  fieldLabelOpacity,
  regionLabelOpacity,
} from "@/lib/timeline";

interface SceneCanvasProps {
  isMobile: boolean;
  /**
   * "scroll" mounts the full scroll-scrubbed scene (camera path,
   * concentration field, region separation, drug fronts, annotations).
   * "static" — used under prefers-reduced-motion — mounts only the tissue +
   * lighting with a fixed camera and a slow idle rotation, per the spec's
   * simplified fallback. Nothing that animates a physical process is mounted
   * there: a transport front that cannot advance would be worse than absent.
   */
  mode: "scroll" | "static";
}

/**
 * Under `frameloop="demand"` R3F only draws when something asks it to. The
 * static scene has no animation left to ask (see the `mode === "static"`
 * branches in GLTFTumor and Particles), so without this the very first frame
 * — and every frame after a Suspense boundary resolves — could be missed.
 *
 * No dependency array on purpose: it runs after EVERY commit inside the
 * canvas, which is exactly the set of moments at which the scene's content
 * can have changed. Resize is already handled by R3F itself.
 */
function DrawOnCommit() {
  const invalidate = useThree((s) => s.invalidate);
  useEffect(() => {
    invalidate();
  });
  return null;
}

/**
 * Owns the R3F `<Canvas>` element and every render-loop performance knob:
 * DPR is capped (lower on mobile), shadows are off entirely (see
 * SceneLighting), geometry detail is reduced on mobile, and under
 * prefers-reduced-motion the render loop is switched off outright.
 *
 * It also owns the model-availability probe, because two different consumers
 * need the same answer: `TumorModel` (which mesh to mount) and `Loader`
 * (whether a download is coming at all). The overlay lives outside the Canvas
 * as plain DOM, so the probe cannot live inside it.
 */
export default function SceneCanvas({ isMobile, mode }: SceneCanvasProps) {
  const modelAvailable = useModelAvailability();
  const isStatic = mode === "static";

  /*
   * Set when something inside the canvas tree threw through React. The
   * narrative, the ledgers and both CTAs are plain DOM outside this component
   * and stay fully readable, so the right response is to render no scene and
   * get everything else out of the way -- including the loading overlay, which
   * is `position: fixed; inset: 0` on an opaque background and would otherwise
   * wait forever for a model that is not coming, i.e. a black page.
   */
  const [boundaryFailed, setBoundaryFailed] = useState(false);
  const onSceneFailure = useCallback(() => setBoundaryFailed(true), []);
  /*
   * No WebGL at all is the other way this scene does not exist, and it does NOT
   * arrive through the boundary — R3F builds its renderer asynchronously, so
   * the failure is an unhandled rejection rather than a React error. See
   * useWebGLSupport for the trace. Detected up front instead, which also means
   * the dead canvas is never mounted.
   */
  const webglSupported = useWebGLSupport();

  /*
   * DPR CEILING, AND WHY DESKTOP IS NOT 2.
   *
   * Nothing drawn in this canvas has a hard edge that a second device pixel
   * resolves: the tissue is an additively-blended point cloud, the field is a
   * smooth viridis ramp, the contour lines are 2.4 px wide and the reticles are
   * rings. Every piece of type on the page is DOM, at full device resolution,
   * and is unaffected by this number.
   *
   * What it does cost is fill rate, quadratically. At 1512x950 on a 2x display
   * the drawing buffer was 3024x1900 = 5.74 Mpx of mostly transparent,
   * heavily-overdrawn geometry, redrawn every frame for the entire visit. 1.75
   * makes it 4.40 Mpx -- 23% fewer fragments for a difference that has to be
   * hunted for at 200% zoom on a stationary frame, on a scene that is never
   * stationary.
   *
   * Kept as a RANGE rather than a fixed number: R3F clamps the real
   * devicePixelRatio into it, so a 1x display still renders at 1x rather than
   * being upsampled to the ceiling.
   */
  const dpr: [number, number] = isMobile ? [1, 1.5] : [1, 1.75];

  /*
   * No context, no canvas. Returning before the Canvas is mounted at all is
   * what keeps three.js from constructing a renderer it cannot use — measured,
   * that attempt logged two uncaught errors and left a dead <canvas> in the
   * DOM. `sceneFailed` tells the overlay to stand down; everything the page
   * actually argues is plain DOM in the layers above this component and is
   * unaffected.
   */
  if (!webglSupported) {
    return (
      <div className="scene-canvas">
        <Loader expectModel={modelAvailable} sceneFailed />
      </div>
    );
  }

  return (
    <div className="scene-canvas">
      {/*
        aria-hidden on the canvas CONTAINER, not on the <canvas> element, and
        that is deliberate: drei's <Html> portals the 3D annotation labels into
        `gl.domElement.parentNode`, so those labels are siblings of the canvas
        inside this container. They are decorative duplicates of information the
        beats already carry in flow — the mobile breakpoint drops them entirely
        for exactly that reason (see Annotation's `reticleOnly`) — and unhidden
        they were read out to assistive technology at every scroll position,
        including while at 0% opacity. One attribute here removes the canvas and
        those labels together. The loading overlay is OUTSIDE this container and
        keeps its progressbar semantics.
      */}
      <SceneErrorBoundary label="WebGL scene" fallback={null} onFailure={onSceneFailure}>
        <Canvas
          aria-hidden="true"
          dpr={dpr}
          /*
           * The whole point of the reduced-motion fallback is that nothing
           * moves. A canvas that keeps compositing a full-viewport WebGL frame
           * sixty times a second to display an unchanging image is not a
           * reduced-motion experience, it is the same power draw with the
           * motion hidden. On demand, R3F draws once and then stops.
           */
          frameloop={isStatic ? "demand" : "always"}
          gl={{ antialias: true, powerPreference: "high-performance", alpha: false }}
          camera={{ position: [0, 0.05, 4.4], fov: 32, near: 0.1, far: 30 }}
          shadows={false}
          onCreated={({ scene, gl }) => {
            scene.background = new THREE.Color("#07080a");
            scene.fog = new THREE.Fog("#07080a", 6, 14);
            gl.setClearColor("#07080a", 1);
          }}
        >
          {isStatic && <DrawOnCommit />}
          <Suspense fallback={null}>
            <SceneLighting />
            <TumorModel
              available={modelAvailable}
              lowDetail={isMobile}
                reducedMotion={isStatic}
              coreChildren={
              mode === "scroll" && (
                /*
                  Nested INSIDE coreGroup (see ProceduralTumor.tsx), not a
                  sibling with a hand-computed world coordinate: coreGroup
                  sits under rootGroup, which spins continuously forever
                  (autoRotation += delta * 0.12) and also translates on
                  scroll-driven separation, so no static world anchor can
                  stay correct — it was previously computed for one instant
                  and drifted from there, ending up pointing at empty space
                  within moments (the actual bug being fixed here). As a
                  child, Three.js's own parent-child transform composition
                  keeps this correct automatically. Anchor is in coreGroup's
                  own local space, well inside its ~0.62 mass radius.
                */
                <Annotation
                  anchor={[0.16, 0.21, 0.42]}
                  labelOffset={[-0.85, -0.2, 0.35]}
                  label="Tumor core"
                  detail="High interstitial pressure"
                  align="left"
                  getOpacity={regionLabelOpacity}
                />
              )
            }
            edemaChildren={
              mode === "scroll" && (
                // Nested inside edemaGroup; local anchor at the shell's
                // actual populated radius (edemaRadius 0.92 * shell [0.55,1]
                // -> points fill [0.51, 0.92]), not the shell's empty
                // center, which the previous version pointed at.
                <Annotation
                  anchor={[0.54, 0.45, 0.27]}
                  labelOffset={[0.55, 0.2, 0]}
                  label="Diffusion front"
                  detail="Concentration gradient"
                  getOpacity={regionLabelOpacity}
                />
              )
            }
            marginChildren={
              mode === "scroll" && (
                // Nested inside marginGroup; local anchor at its populated
                // shell radius (marginRadius 1.34 * shell [0.82,1] ->
                // [1.10, 1.34]).
                <Annotation
                  anchor={[-0.39, 0.2, 1.12]}
                  labelOffset={[-0.85, -0.1, 0.55]}
                  label="Perivascular rim"
                  detail="Where drug enters"
                  align="left"
                  getOpacity={regionLabelOpacity}
                />
              )
            }
          />

          {mode === "scroll" && (
            <>
              <CameraRig />
              <Particles count={isMobile ? 320 : 900} />
              {/* Everything below is positioned in the measured section's own
                  coordinates, so it all has to travel with the section when
                  the final scene slides it aside. See TransportOverlay. */}
              <TransportOverlay>
              <FieldContours lowDetail={isMobile}>
                {/*
                  Nested inside FieldContours' own <group> (which sways every
                  frame) rather than a sibling with a pre-multiplied world
                  coordinate — same bug/fix as the core/edema/margin
                  annotations above. Anchor is a real vertex of the SECOND
                  drawn iso-level (0.36 mmHg) in that group's own local space:
                  its upper-right shoulder, which is the outermost point of the
                  contour set that is still legible. Anchoring on the true
                  outermost level (0.22 mmHg) was tried first and read as the
                  old bug — that level sits at viridis 0.30, dark by
                  construction, so at the opacity the timeline gives section
                  02's centre the reticle appeared to float in empty space
                  next to a line nobody could see. The contour set sits in the
                  LOWER-LEFT of the section because that is where this tissue's
                  vasculature is densest, so the anchor moved there with it —
                  the previous one was on the right, where the procedural rings
                  used to be.

                  The label leads RIGHT, across the section, deliberately:
                  section 02's copy is left-aligned, and the solved lobe is on
                  the same side, so a leftward label would sit under the
                  headline.

                  detail= names the quantity, the ONE display transformation
                  applied to it (a 250 um local mean of the solved field —
                  see lib/ifpContours.ts) and the mmHg range the viridis ramp
                  spans, so the contours can never be read as the drug
                  concentration ConcentrationField paints with the same ramp.
                  It replaces "Illustrative, not solved", which was true of
                  the noise-walk ellipses it labelled and is not true of these.

                  The label is `white-space: nowrap` and about 250px wide at
                  its full length, which is two thirds of a 390px viewport and
                  ran off the right edge there. Below the 820px breakpoint it
                  therefore drops the level range and is re-placed. What it
                  keeps is the load-bearing half: that the
                  field is SOLVED and that what is drawn is a 250 um mean of
                  it. The numeric range is the part a narrow screen can afford
                  to lose — the ramp is still monotone without it, and the
                  quantity is still named.

                  The mobile x offset is not a smaller version of the desktop
                  one, it is LARGER relative to the label's travel: over the
                  window this callout is visible the projected label sweeps
                  about 140px leftward at 390px wide (the section rotates and
                  slides under it), which took it off the left edge while it
                  was still at 38% opacity. The offset is set from the left
                  extreme of that sweep, not from its resting position.
                */}
                {/*
                  THE NUMERIC RANGE WAS REMOVED, and it must not come back.

                  The label used to end with the iso-contour levels
                  (`IFP_LEVEL_RANGE_LABEL`, 0.22-0.65 mmHg). Two independent
                  reasons that was wrong, both measured:

                  1. It quoted the wrong quantity. Those are contour LEVELS of
                     a 250 um coarse-grained field; the raw field spans
                     0-12.93 mmHg with a median in tissue of 0.0036. This
                     module's own header (ifpContours.ts, "the page must never
                     quote one as the other") warned against exactly this, and
                     the annotation did it anyway.
                  2. The magnitude is not defensible. results/ifp_validation.md
                     scored the solved field against wick-in-needle
                     measurements and it comes out 41x LOW on breast against
                     15 +/- 9 mmHg (PMID 1423283) -- and refining the grid
                     makes it worse, not better. The GBM comparison is
                     NOT TESTABLE rather than a pass: Boucher's intracranial
                     2.0 +/- 2.5 mmHg (SD, PMID 9062403) is statistically
                     indistinguishable from normal brain, so it cannot
                     separate this engine from one with no tumour in it.

                  What survives is the SHAPE, which is what this beat is for:
                  the field is genuinely solved on measured geometry, and
                  showing where pressure rises is a real claim. The absolute
                  scale is not, so the label now says "relative, uncalibrated"
                  instead of a number the project cannot stand behind.
                */}
                <Annotation
                  anchor={[-0.831, -0.33, 0]}
                  labelOffset={isMobile ? [0.48, 0.5, 0] : [0.78, 0.58, 0]}
                  label="Interstitial pressure"
                  detail={
                    isMobile
                      ? "Solved · relative"
                      : "Solved · relative, uncalibrated · 250 µm mean"
                  }
                  getOpacity={fieldLabelOpacity}
                  /* At 390px the copy column is full-bleed and beat 02 now
                     carries a four-row ledger, so there is no band left for a
                     floating callout — measured, the label overlapped the body
                     text at every scroll position where both were visible. The
                     ledger names the quantity and its three registered
                     verdicts, so the chip is redundant there rather than lost.
                     See Annotation's `reticleOnly` prop. */
                  reticleOnly={isMobile}
                />
              </FieldContours>
                <ConcentrationField />
                {/* Vessels, then the fronts that leave them, in that order
                    down the scroll — the two halves of one idea, and the
                    reason the old TrajectoryPath's draw-on effect was worth
                    reviving under a different name (see DrugFront.tsx)
                    rather than deleting. */}
                <VesselNetwork />
                <DrugFront />
                {/* On a measured endothelial density maximum (the first
                    ENTRY_POINT), lifted just off the section's plane so it is
                    not co-planar with the vessel stand-in running through the
                    same point. */}
                <ConcentrationSource
                  position={[ENTRY_POINTS[0][0], ENTRY_POINTS[0][1], 0.05]}
                />
              </TransportOverlay>
            </>
          )}

            {isStatic && <Particles count={isMobile ? 120 : 260} reducedMotion />}
          </Suspense>
          <Preload all />
        </Canvas>
      </SceneErrorBoundary>
      <Loader expectModel={modelAvailable} sceneFailed={boundaryFailed} />
    </div>
  );
}
