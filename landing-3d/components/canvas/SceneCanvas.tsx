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

              {/*
                Anchor points below are chosen to sit ON populated geometry,
                not just at a convenient-looking coordinate:

                - "Left temporal lobe" used to sit at a bare [1.3,0.9,0.3],
                  a leftover from before AnatomyModel was reworked from a
                  translucent icosahedron shell to explicit axial/sagittal/
                  coronal slice-contour rings (see AnatomyModel.tsx) — that
                  point falls between ring planes, in empty space. It now
                  sits on the coronal ring (level -0.05, radiusA 0.88,
                  radiusB 0.9) at its scroll-quarter angle, transformed
                  through that group's own scale/position.
                - "Peritumoral edema" / "Surrounding tissue" used to sit
                  exactly at the *unrotated local origin* of the edema/
                  margin point-cloud groups — i.e. the translation target
                  those groups animate to during separation (see
                  ProceduralTumor.tsx's coreGroup/edemaGroup/marginGroup
                  `.position.set(...)`, which uses these same numbers). Both
                  clouds are hollow shells (sampleHaloPointCloud fills only
                  a `shell` fraction of the local radius outward), so that
                  origin is the one place guaranteed to have zero points —
                  a leader line that always finds the exact center of the
                  doughnut hole. Anchors now extend further out, to each
                  shell's populated radius (edema ~0.71, margin ~1.22),
                  along the same direction — landing on actual points
                  regardless of the model's continuous idle rotation, since
                  both shells are isotropic (uniform random direction).
              */}
              <Annotation
                anchor={[1.55, 1.14, -0.19]}
                labelOffset={[0.5, 0.25, 0]}
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
                anchor={[0.65, 0.62, 0.3]}
                labelOffset={[0.55, 0.2, 0]}
                label="Peritumoral edema"
                detail="Surrounding swelling"
                getOpacity={segmentationLabelOpacity}
              />
              <Annotation
                anchor={[-0.24, 0.15, -1.8]}
                labelOffset={[-0.85, -0.1, 0.55]}
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
