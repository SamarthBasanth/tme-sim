"use client";

import { Suspense } from "react";
import { Canvas } from "@react-three/fiber";
import { Preload } from "@react-three/drei";
import * as THREE from "three";
import CameraRig from "./CameraRig";
import SceneLighting from "./SceneLighting";
import TumorModel from "./TumorModel";
import AnatomyModel from "./AnatomyModel";
import Particles from "./Particles";
import ScanPlane from "./ScanPlane";
import TrajectoryPath from "./TrajectoryPath";
import Annotation from "./Annotation";
import { HighlightMarker, CriticalStructures } from "./Markers";
import Loader from "./Loader";
import {
  anatomyLabelOpacity,
  segmentationLabelOpacity,
} from "@/lib/timeline";

interface SceneCanvasProps {
  isMobile: boolean;
  /**
   * "scroll" mounts the full scroll-scrubbed scene (camera path, scan
   * plane, segmentation, trajectory, annotations). "static" — used under
   * prefers-reduced-motion — mounts only the tumor + lighting with a fixed
   * camera and a slow idle rotation, per the spec's simplified fallback.
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
                  detail="Primary mass"
                  align="left"
                  getOpacity={segmentationLabelOpacity}
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
                  label="Peritumoral edema"
                  detail="Surrounding swelling"
                  getOpacity={segmentationLabelOpacity}
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
                  label="Surrounding tissue"
                  detail="Margin boundary"
                  align="left"
                  getOpacity={segmentationLabelOpacity}
                />
              )
            }
          />

          {mode === "scroll" && (
            <>
              <CameraRig />
              <AnatomyModel lowDetail={isMobile}>
                {/*
                  Nested inside AnatomyModel's own <group> (scale
                  [1.9,2.3,1.7], continuous rotation.y += delta*0.02) rather
                  than a sibling with a pre-multiplied world coordinate —
                  same bug/fix as the core/edema/margin annotations above.
                  Anchor is a point on the coronal ring itself (level -0.05,
                  radiusA 0.88, radiusB 0.9 — see RING_SPECS), in the ring's
                  own local space, before that group's scale/position is
                  applied (which now happens automatically since this is a
                  child of the same group the ring <Line> is).
                */}
                <Annotation
                  anchor={[0.875, 0.09, -0.05]}
                  labelOffset={[0.32, 0.25, 0]}
                  label="Left temporal lobe"
                  detail="Simplified anatomical context"
                  getOpacity={anatomyLabelOpacity}
                />
              </AnatomyModel>
              <Particles count={isMobile ? 320 : 900} />
              <ScanPlane />
              <TrajectoryPath />
              <CriticalStructures />
              <HighlightMarker position={[0.3, 0.25, 0.15]} />
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
