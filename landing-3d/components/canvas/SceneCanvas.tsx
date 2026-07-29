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
          <TumorModel lowDetail={isMobile} reducedMotion={mode === "static"} />

          {mode === "scroll" && (
            <>
              <CameraRig />
              <AnatomyModel lowDetail={isMobile} />
              <Particles count={isMobile ? 320 : 900} />
              <ScanPlane />
              <TrajectoryPath />
              <CriticalStructures />
              <HighlightMarker position={[0.3, 0.25, 0.15]} />

              <Annotation
                anchor={[1.3, 0.9, 0.3]}
                labelOffset={[0.55, 0.28, 0]}
                label="Left temporal lobe"
                detail="Simplified anatomical context"
                getOpacity={anatomyLabelOpacity}
              />
              <Annotation
                anchor={[0, -0.55, 0]}
                labelOffset={[-0.85, -0.2, 0.35]}
                label="Tumor core"
                detail="Primary mass"
                align="left"
                getOpacity={segmentationLabelOpacity}
              />
              <Annotation
                anchor={[0.05, 0.62, 0.1]}
                labelOffset={[0.8, 0.22, 0]}
                label="Peritumoral edema"
                detail="Surrounding swelling"
                getOpacity={segmentationLabelOpacity}
              />
              <Annotation
                anchor={[-0.08, 0.05, -0.6]}
                labelOffset={[-0.9, -0.05, -0.35]}
                label="Surrounding tissue"
                detail="Margin boundary"
                align="left"
                getOpacity={segmentationLabelOpacity}
              />
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
