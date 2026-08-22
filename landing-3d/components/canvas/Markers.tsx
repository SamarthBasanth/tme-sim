"use client";

import { useEffect, useMemo, useRef } from "react";
import { useFrame } from "@react-three/fiber";
import { Billboard } from "@react-three/drei";
import * as THREE from "three";
import { useScrollProgressRef } from "@/components/scroll/ScrollProgressContext";
import { exposureHighlight, vesselOpacity } from "@/lib/timeline";
import { damp } from "@/lib/scrollMath";
import { ENTRY_POINTS } from "./shaders/viridis";

/**
 * Concentration ramp stops used by the WebGL materials in this file, sampled
 * from the product's own viridis table (app/web/js/gl/ramps.js) at roughly
 * t = 0.95 / 0.55 / 0.2. Mirrors `--conc-high` in globals.css; three.js cannot
 * read a CSS custom property, so the anchor value is repeated, not forked.
 */
const CONC_HIGH = "#d8e219";
const CONC_MID = "#28ae80";
const CONC_LOW = "#424086";

/**
 * A point source of drug, sitting on the vasculature it enters from.
 *
 * This was an amber crosshair/reticle -- dot, ring and four radial ticks --
 * the standard "flagged point" convention in imaging overlays, built for an
 * "AI flagged region" beat this product does not have. Amber carried no
 * meaning beyond "look here", and the page reserved it site-wide for exactly
 * that.
 *
 * Rebuilt as three concentric bands stepping DOWN the concentration ramp
 * outward -- bright at the source, mid, then low -- so the marker is a
 * miniature of the quantity the whole scene is about instead of a flag. The
 * ticks are gone with the reticle they belonged to. It stays a camera-facing
 * `Billboard`, which is what keeps a ringed marker from reading as a tilted
 * Saturn from off-axis (the reason the original stopped being a plain sphere
 * plus an orbiting ring).
 */
export function ConcentrationSource({
  position,
}: {
  position: readonly [number, number, number];
}) {
  const progressRef = useScrollProgressRef();
  const coreMaterialRef = useRef<THREE.MeshBasicMaterial>(null);
  const midMaterialRef = useRef<THREE.MeshBasicMaterial>(null);
  const outerMaterialRef = useRef<THREE.MeshBasicMaterial>(null);
  const billboardRef = useRef<THREE.Group>(null);
  const opacitySmoothed = useRef(0);

  useFrame((state, delta) => {
    const target = exposureHighlight(progressRef.current);
    opacitySmoothed.current = damp(opacitySmoothed.current, target, 4, delta);
    const o = opacitySmoothed.current;
    const pulse = 1 + Math.sin(state.clock.elapsedTime * 2.4) * 0.08 * o;

    if (coreMaterialRef.current) coreMaterialRef.current.opacity = o;
    if (midMaterialRef.current) midMaterialRef.current.opacity = o * 0.7;
    if (outerMaterialRef.current) outerMaterialRef.current.opacity = o * 0.45;
    if (billboardRef.current) {
      billboardRef.current.scale.setScalar(pulse);
      // A `Billboard` re-aims itself at the camera every frame; three circles
      // at opacity 0 were being oriented, sorted and blended for the whole
      // scroll to draw nothing.
      billboardRef.current.visible = o > 0.004;
    }
  });

  return (
    <group position={position}>
      <Billboard ref={billboardRef}>
        <mesh>
          <circleGeometry args={[0.026, 24]} />
          <meshBasicMaterial
            ref={coreMaterialRef}
            color={CONC_HIGH}
            transparent
            opacity={0}
            depthWrite={false}
          />
        </mesh>
        <mesh>
          <ringGeometry args={[0.036, 0.048, 40]} />
          <meshBasicMaterial
            ref={midMaterialRef}
            color={CONC_MID}
            transparent
            opacity={0}
            depthWrite={false}
            side={THREE.DoubleSide}
          />
        </mesh>
        <mesh>
          <ringGeometry args={[0.058, 0.074, 40]} />
          <meshBasicMaterial
            ref={outerMaterialRef}
            color={CONC_LOW}
            transparent
            opacity={0}
            depthWrite={false}
            side={THREE.DoubleSide}
          />
        </mesh>
      </Billboard>
    </group>
  );
}

/**
 * Faint translucent vessel stand-ins, revealed late in the scroll.
 *
 * Kept, and moved onto the tissue. Vasculature is the one piece of scene
 * furniture that is entirely on-message for a transport product -- it is where
 * the drug comes from -- but these were laid out for a volumetric placeholder
 * mass: long curves swinging through +/-0.6 of depth around a section that has
 * no depth, which rendered as grey scaffolding poles crossing the frame in
 * front of a flat slide.
 *
 * Every control point below is now a measured endothelial density maximum from
 * the section itself (ENTRY_POINTS), in the section's own plane. The positions
 * are measured; the connectivity between them -- which point joins which -- is
 * drawn, not reconstructed, which is why these stay thin and faint and are
 * called stand-ins.
 */
const VESSEL_PATHS: ReadonlyArray<ReadonlyArray<readonly [number, number]>> = [
  [ENTRY_POINTS[0], ENTRY_POINTS[1], ENTRY_POINTS[3]],
  [ENTRY_POINTS[2], ENTRY_POINTS[4], ENTRY_POINTS[5]],
];

/** Just in front of the section plane (the asset itself is flat at z = 0). */
const VESSEL_Z = 0.01;

export function VesselNetwork() {
  const progressRef = useScrollProgressRef();
  const groupRef = useRef<THREE.Group>(null);
  const opacitySmoothed = useRef(0);

  const geometries = useMemo(
    () =>
      VESSEL_PATHS.map((path) => {
        const curve = new THREE.CatmullRomCurve3(
          path.map(([x, y]) => new THREE.Vector3(x, y, VESSEL_Z)),
        );
        return new THREE.TubeGeometry(curve, 48, 0.007, 6, false);
      }),
    [],
  );

  // Built outside R3F's declarative tree, so they need an explicit dispose.
  useEffect(() => () => geometries.forEach((g) => g.dispose()), [geometries]);

  useFrame((_, delta) => {
    const target = vesselOpacity(progressRef.current);
    opacitySmoothed.current = damp(opacitySmoothed.current, target, 3, delta);
    if (groupRef.current) {
      groupRef.current.visible = opacitySmoothed.current > 0.004;
      groupRef.current.children.forEach((child) => {
        const mesh = child as THREE.Mesh;
        const material = mesh.material as THREE.MeshStandardMaterial;
        material.opacity = opacitySmoothed.current * 0.55;
      });
    }
  });

  return (
    <group ref={groupRef} visible={false}>
      {geometries.map((geometry, i) => (
        <mesh key={i} geometry={geometry}>
          <meshStandardMaterial
            color="#9fc6d8"
            transparent
            opacity={0}
            depthWrite={false}
            roughness={0.4}
          />
        </mesh>
      ))}
    </group>
  );
}
