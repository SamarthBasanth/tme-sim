"use client";

import { useEffect, useMemo, useRef } from "react";
import { useFrame } from "@react-three/fiber";
import * as THREE from "three";
import { useScrollProgressRef } from "@/components/scroll/ScrollProgressContext";
import { finalOffsetAmount, separationAmount } from "@/lib/timeline";
import { damp } from "@/lib/scrollMath";
import {
  createEdemaGeometry,
  createLobulatedMass,
  createMarginGeometry,
  createVesselCurves,
} from "./tumorGeometry";
import { EDEMA_FRAGMENT, EDEMA_VERTEX } from "./shaders/edemaShader";

interface ProceduralTumorProps {
  lowDetail?: boolean;
  /** Slows/removes scroll-coupled and ambient rotation for prefers-reduced-motion. */
  reducedMotion?: boolean;
}

const EDEMA_COLOR = "#7fc2d6";
const MARGIN_COLOR = "#cfd6d9";
const VESSEL_COLOR = "#6b1f2c";
const RIM_GLOW_COLOR = "#ff8a72";

/**
 * Stand-in for a real `/models/tumor.glb`, redesigned to actually read as
 * tumor tissue rather than a glossy sphere. Grounded in how glioblastoma
 * looks on contrast MRI and in surgical-planning 3D reconstructions: an
 * irregular, lobulated mass (domain-warped multi-octave noise displacement,
 * not a smooth metaball), stylized "ring enhancement" — a darker necrotic
 * cavity ringed by a brighter enhancing band, baked into vertex colors — a
 * matte, low-gloss surface (high roughness, no clearcoat, a cheap
 * backface-shell fresnel rim standing in for subsurface scattering), thin
 * procedural vessels, and a genuinely diffuse/soft-edged edema volume
 * instead of a second hard sphere. See `tumorGeometry.ts` for the geometry
 * math and its citations. Built from three groups — core / edema / margin —
 * because the "Segmentation" scene pulls those apart on scroll; swap this
 * component out (see TumorModel.tsx) once a real scan-derived GLB exists.
 */
export default function ProceduralTumor({
  lowDetail = false,
  reducedMotion = false,
}: ProceduralTumorProps) {
  const progressRef = useScrollProgressRef();

  const coreGroup = useRef<THREE.Group>(null);
  const edemaGroup = useRef<THREE.Group>(null);
  const marginGroup = useRef<THREE.Group>(null);
  const rootGroup = useRef<THREE.Group>(null);
  const edemaMaterialRef = useRef<THREE.ShaderMaterial>(null);

  const sepSmoothed = useRef(0);
  const offsetSmoothed = useRef(0);
  const autoRotation = useRef(0);
  const scrollRotation = useRef(0);

  const coreRadius = 0.62;
  const edemaRadius = 0.92;
  const marginRadius = 1.34;

  const massDetail = lowDetail ? 6 : 9;
  const edemaDetail = lowDetail ? 4 : 6;
  const marginDetail = lowDetail ? 3 : 5;

  const massGeometry = useMemo(
    () => createLobulatedMass(coreRadius, massDetail),
    [massDetail],
  );
  const edemaGeometry = useMemo(
    () => createEdemaGeometry(edemaRadius, edemaDetail),
    [edemaDetail],
  );
  const marginGeometry = useMemo(
    () => createMarginGeometry(marginRadius, marginDetail),
    [marginDetail],
  );

  const vesselGeometries = useMemo(() => {
    const curves = createVesselCurves({ count: lowDetail ? 3 : 5, radius: coreRadius });
    return curves.map((curve) => new THREE.TubeGeometry(curve, 28, 0.011, 5, false));
  }, [lowDetail]);

  // Dispose everything created above on unmount / whenever a memo above
  // recomputes (e.g. the mobile lowDetail toggle) — R3F does not do this
  // automatically for objects constructed outside its own JSX tree.
  useEffect(() => {
    return () => {
      massGeometry.dispose();
    };
  }, [massGeometry]);
  useEffect(() => {
    return () => {
      edemaGeometry.dispose();
    };
  }, [edemaGeometry]);
  useEffect(() => {
    return () => {
      marginGeometry.dispose();
    };
  }, [marginGeometry]);
  useEffect(() => {
    return () => {
      vesselGeometries.forEach((g) => g.dispose());
    };
  }, [vesselGeometries]);

  const edemaUniforms = useMemo(
    () => ({
      uColor: { value: new THREE.Color(EDEMA_COLOR) },
      uOpacity: { value: 0.8 },
      uTime: { value: 0 },
    }),
    [],
  );

  useFrame((state, delta) => {
    const p = progressRef.current;
    sepSmoothed.current = damp(sepSmoothed.current, separationAmount(p), 3.2, delta);
    offsetSmoothed.current = damp(offsetSmoothed.current, finalOffsetAmount(p), 3.2, delta);

    const sep = sepSmoothed.current;
    if (coreGroup.current) coreGroup.current.position.set(0, -0.55 * sep, 0);
    if (edemaGroup.current) edemaGroup.current.position.set(0.05 * sep, 0.62 * sep, 0.1 * sep);
    if (marginGroup.current) marginGroup.current.position.set(-0.08 * sep, 0.05 * sep, -0.6 * sep);

    if (edemaMaterialRef.current) {
      edemaMaterialRef.current.uniforms.uTime.value = state.clock.elapsedTime;
    }

    if (rootGroup.current) {
      // Gentle continuous autonomous rotation, plus a slow scroll-coupled
      // turn so the model always feels connected to scroll input even
      // within a single section (not just at section boundaries). Reduced
      // to a slow idle spin with no scroll coupling under reduced-motion.
      autoRotation.current += delta * (reducedMotion ? 0.035 : 0.12);
      scrollRotation.current = reducedMotion
        ? 0
        : damp(scrollRotation.current, p * 0.6, 2.4, delta);
      rootGroup.current.rotation.y = autoRotation.current + scrollRotation.current;
      rootGroup.current.position.x = reducedMotion ? 0 : offsetSmoothed.current * 1.15;
    }
  });

  return (
    <group ref={rootGroup}>
      {/* Core — the tumor mass itself: lobulated, vertex-colored necrotic
          core / enhancing rim heterogeneity, matte surface, a thin
          fresnel-rim "glow" shell standing in for subsurface scattering,
          and a few thin vessel curves hugging the surface. */}
      <group ref={coreGroup}>
        <mesh geometry={massGeometry} castShadow={false} receiveShadow={false}>
          <meshStandardMaterial
            vertexColors
            color="#ffffff"
            roughness={0.88}
            metalness={0.03}
            emissive="#200a0a"
            emissiveIntensity={0.06}
          />
        </mesh>

        {/* Backface-shell fresnel rim: a slightly larger, inverted copy of
            the mass rendered additively so silhouette edges pick up a warm
            glow — a cheap, robust stand-in for subsurface scattering that
            doesn't require touching Three's internal shader chunks. */}
        <mesh geometry={massGeometry} scale={1.02}>
          <meshBasicMaterial
            color={RIM_GLOW_COLOR}
            transparent
            opacity={0.12}
            side={THREE.BackSide}
            depthWrite={false}
            blending={THREE.AdditiveBlending}
          />
        </mesh>

        {vesselGeometries.map((geometry, i) => (
          <mesh key={i} geometry={geometry}>
            <meshStandardMaterial
              color={VESSEL_COLOR}
              roughness={0.55}
              metalness={0.05}
              emissive="#3a0f18"
              emissiveIntensity={0.2}
            />
          </mesh>
        ))}
      </group>

      {/* Edema — diffuse, soft-edged halo. Opacity is driven almost
          entirely by a fresnel + noise shader (see edemaShader.ts) so the
          volume actually dissolves at its boundary instead of presenting a
          crisp transparent-sphere silhouette. */}
      <group ref={edemaGroup}>
        <mesh geometry={edemaGeometry}>
          <shaderMaterial
            ref={edemaMaterialRef}
            vertexShader={EDEMA_VERTEX}
            fragmentShader={EDEMA_FRAGMENT}
            uniforms={edemaUniforms}
            transparent
            depthWrite={false}
            side={THREE.DoubleSide}
          />
        </mesh>
      </group>

      {/* Margin — faint irregular boundary shell of surrounding tissue. */}
      <group ref={marginGroup}>
        <mesh geometry={marginGeometry}>
          <meshStandardMaterial
            color={MARGIN_COLOR}
            roughness={1}
            metalness={0}
            transparent
            opacity={0.06}
            side={THREE.DoubleSide}
            depthWrite={false}
          />
        </mesh>
      </group>
    </group>
  );
}
