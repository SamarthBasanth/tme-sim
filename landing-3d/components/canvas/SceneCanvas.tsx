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
                  annotations above. Anchor is a point on the outermost
                  contour itself (level 1.0, so radius 0.88 of the section's
                  measured half-extent — see CONTOUR_SPECS), in that group's
                  own local space, at about 35 degrees round the loop — high
                  enough that the label clears the right-aligned headline in
                  section 03, which it collided with when anchored at the
                  contour's widest point.
                */}
                <Annotation
                  anchor={[0.757, 0.351, 0]}
                  labelOffset={[0.42, 0.52, 0]}
                  label="Pressure contours"
                  detail="Illustrative, not solved"
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
