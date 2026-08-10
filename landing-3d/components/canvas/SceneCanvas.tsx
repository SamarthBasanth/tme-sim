"use client";

import { Suspense } from "react";
import { Canvas } from "@react-three/fiber";
import { Preload } from "@react-three/drei";
import * as THREE from "three";
import CameraRig from "./CameraRig";
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
import { IFP_LEVEL_RANGE_LABEL } from "@/lib/ifpContours";
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
 * Owns the R3F `<Canvas>` element and every render-loop performance knob:
 * DPR is capped (lower on mobile), shadows are off entirely (see
 * SceneLighting), and geometry detail is reduced on mobile. The loading
 * overlay lives outside the Canvas (plain DOM) but reads the same
 * three.js LoadingManager via drei's `useProgress`.
 */
export default function SceneCanvas({ isMobile, mode }: SceneCanvasProps) {
  const dpr: [number, number] = isMobile ? [1, 1.5] : [1, 2];

  return (
    <div className="scene-canvas">
      <Canvas
        dpr={dpr}
        gl={{ antialias: true, powerPreference: "high-performance", alpha: false }}
        camera={{ position: [0, 0.05, 4.4], fov: 32, near: 0.1, far: 30 }}
        shadows={false}
        onCreated={({ scene, gl }) => {
          scene.background = new THREE.Color("#07080a");
          scene.fog = new THREE.Fog("#07080a", 6, 14);
          gl.setClearColor("#07080a", 1);
        }}
      >
        <Suspense fallback={null}>
          <SceneLighting />
          <TumorModel
            lowDetail={isMobile}
            reducedMotion={mode === "static"}
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
                <Annotation
                  anchor={[-0.831, -0.33, 0]}
                  labelOffset={isMobile ? [0.48, 0.5, 0] : [0.78, 0.58, 0]}
                  label="Interstitial pressure"
                  detail={
                    isMobile
                      ? "Solved · 250 µm mean"
                      : `Solved · 250 µm mean · ${IFP_LEVEL_RANGE_LABEL}`
                  }
                  getOpacity={fieldLabelOpacity}
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

          {mode === "static" && <Particles count={isMobile ? 120 : 260} />}
        </Suspense>
        <Preload all />
      </Canvas>
      <Loader />
    </div>
  );
}
