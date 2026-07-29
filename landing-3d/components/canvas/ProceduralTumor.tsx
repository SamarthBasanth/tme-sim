"use client";

import { useEffect, useMemo, useRef } from "react";
import { useFrame, useThree } from "@react-three/fiber";
import * as THREE from "three";
import { useScrollProgressRef } from "@/components/scroll/ScrollProgressContext";
import { finalOffsetAmount, highlightOpacity, separationAmount } from "@/lib/timeline";
import { damp } from "@/lib/scrollMath";
import {
  createVesselCurves,
  sampleEdemaPointCloud,
  sampleMarginPointCloud,
  sampleTumorMassPointCloud,
} from "./tumorGeometry";
import { POINT_CLOUD_FRAGMENT, POINT_CLOUD_VERTEX } from "./shaders/pointCloudShader";

interface ProceduralTumorProps {
  lowDetail?: boolean;
  /** Slows/removes scroll-coupled and ambient rotation for prefers-reduced-motion. */
  reducedMotion?: boolean;
}

const VESSEL_COLOR = "#6b1f2c";

/**
 * Stand-in for a real `/models/tumor.glb`, rendered the way the real
 * product renders a specimen — as a point cloud of individual cell-like
 * sprites (see `app/web/js/gl/shaders.js` POINTS_VS/POINTS_FS, and
 * `app/web/js/gl/viewer.js`'s gl.POINTS draw call) — rather than a solid
 * displaced mesh. Each of the three regions (core, edema, margin) is its
 * own point cloud, sampled once in `useMemo` from domain-warped noise (see
 * tumorGeometry.ts) so the cloud itself carries the lobulated, irregular
 * tumour silhouette instead of relying on an invisible mesh boundary. A
 * per-point category attribute (core / necrotic / enhancing-rim within the
 * mass cloud) drives colour through `pointCloudShader.ts`, the R3F analogue
 * of the real viewer's category-attribute channel; `uHighlight` nudges the
 * enhancing category brighter while the "AI Analysis" scan plane is active,
 * bridging the point-cloud category concept to that section the same way
 * the real product's mechanism-attribution colouring works. Built from
 * three groups — core / edema / margin — because the "Segmentation" scene
 * pulls those apart on scroll; swap this component out (see TumorModel.tsx)
 * once a real scan-derived GLB exists.
 */
export default function ProceduralTumor({
  lowDetail = false,
  reducedMotion = false,
}: ProceduralTumorProps) {
  const progressRef = useScrollProgressRef();
  const size = useThree((s) => s.size);
  const dpr = useThree((s) => s.viewport.dpr);

  const coreGroup = useRef<THREE.Group>(null);
  const edemaGroup = useRef<THREE.Group>(null);
  const marginGroup = useRef<THREE.Group>(null);
  const rootGroup = useRef<THREE.Group>(null);

  const massMaterialRef = useRef<THREE.ShaderMaterial>(null);
  const edemaMaterialRef = useRef<THREE.ShaderMaterial>(null);
  const marginMaterialRef = useRef<THREE.ShaderMaterial>(null);

  const sepSmoothed = useRef(0);
  const offsetSmoothed = useRef(0);
  const highlightSmoothed = useRef(0);
  const autoRotation = useRef(0);
  const scrollRotation = useRef(0);

  const coreRadius = 0.62;
  const edemaRadius = 0.92;
  const marginRadius = 1.34;

  const massCount = lowDetail ? 4500 : 9000;
  const edemaCount = lowDetail ? 1100 : 2200;
  const marginCount = lowDetail ? 600 : 1200;

  const massCloud = useMemo(() => sampleTumorMassPointCloud(coreRadius, { count: massCount }), [massCount]);
  const edemaCloud = useMemo(() => sampleEdemaPointCloud(edemaRadius, edemaCount), [edemaCount]);
  const marginCloud = useMemo(() => sampleMarginPointCloud(marginRadius, marginCount), [marginCount]);

  const vesselGeometries = useMemo(() => {
    const curves = createVesselCurves({ count: lowDetail ? 3 : 5, radius: coreRadius });
    return curves.map((curve) => new THREE.TubeGeometry(curve, 28, 0.011, 5, false));
  }, [lowDetail]);

  // Vessel tube geometries are built outside R3F's declarative tree, so
  // (unlike the point clouds' <bufferGeometry> JSX, which R3F disposes
  // automatically) they need an explicit dispose on unmount / recompute.
  useEffect(() => {
    return () => {
      vesselGeometries.forEach((g) => g.dispose());
    };
  }, [vesselGeometries]);

  const massUniforms = useMemo(
    () => ({
      uColor0: { value: new THREE.Color("#8a3230") }, // core tissue
      uColor1: { value: new THREE.Color("#3a1c18") }, // necrotic
      uColor2: { value: new THREE.Color("#e2685c") }, // enhancing rim
      uSize: { value: 0.046 },
      uSizeScale: { value: 1 },
      uOpacity: { value: 1 },
      uHighlight: { value: 0 },
      uHighlightCategory: { value: 2 },
    }),
    [],
  );
  const edemaUniforms = useMemo(
    () => ({
      uColor0: { value: new THREE.Color("#7fc2d6") },
      uColor1: { value: new THREE.Color("#7fc2d6") },
      uColor2: { value: new THREE.Color("#a7dcec") },
      uSize: { value: 0.05 },
      uSizeScale: { value: 1 },
      uOpacity: { value: 0.55 },
      uHighlight: { value: 0 },
      uHighlightCategory: { value: -1 },
    }),
    [],
  );
  const marginUniforms = useMemo(
    () => ({
      uColor0: { value: new THREE.Color("#cfd6d9") },
      uColor1: { value: new THREE.Color("#cfd6d9") },
      uColor2: { value: new THREE.Color("#cfd6d9") },
      uSize: { value: 0.032 },
      uSizeScale: { value: 1 },
      uOpacity: { value: 0.3 },
      uHighlight: { value: 0 },
      uHighlightCategory: { value: -1 },
    }),
    [],
  );

  // Point size follows three.js's own PointsMaterial formula
  // (gl_PointSize = size * scale / -viewZ, scale = canvasHeightPx * 0.5) —
  // kept in sync with the canvas here rather than recomputed per-frame,
  // since it only needs to change on resize. Mutated through the material
  // refs, not the `massUniforms` etc. literals passed as JSX props: R3F
  // clones each uniform entry onto the material the first time it applies
  // the `uniforms` prop (see its ShaderMaterial special-case), so the
  // material's live `uniforms.uSizeScale` is a *different* object from
  // `massUniforms.uSizeScale` after that first commit — mutating the
  // original literal silently does nothing.
  useEffect(() => {
    const scale = size.height * dpr * 0.5;
    [massMaterialRef, edemaMaterialRef, marginMaterialRef].forEach((ref) => {
      if (ref.current) ref.current.uniforms.uSizeScale.value = scale;
    });
  }, [size, dpr]);

  useFrame((_, delta) => {
    const p = progressRef.current;
    sepSmoothed.current = damp(sepSmoothed.current, separationAmount(p), 3.2, delta);
    offsetSmoothed.current = damp(offsetSmoothed.current, finalOffsetAmount(p), 3.2, delta);
    highlightSmoothed.current = damp(highlightSmoothed.current, highlightOpacity(p), 4, delta);

    const sep = sepSmoothed.current;
    if (coreGroup.current) coreGroup.current.position.set(0, -0.55 * sep, 0);
    if (edemaGroup.current) edemaGroup.current.position.set(0.05 * sep, 0.62 * sep, 0.1 * sep);
    if (marginGroup.current) marginGroup.current.position.set(-0.08 * sep, 0.05 * sep, -0.6 * sep);

    // A gentle glow on the enhancing/margin category while the AI-analysis
    // scan is active — the same curve that ignites Markers.tsx's amber
    // marker — so "AI Analysis" reads as flagging a cell category rather
    // than just sweeping a plane past an inert mass.
    if (massMaterialRef.current) {
      massMaterialRef.current.uniforms.uHighlight.value = highlightSmoothed.current;
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
      {/* Core — the tumor mass itself: a dense point cloud carrying its own
          lobulated silhouette, with a per-point category (core / necrotic /
          enhancing rim) standing in for BraTS-style ring enhancement, plus
          a few thin vessel curves hugging the surface. */}
      <group ref={coreGroup}>
        <points frustumCulled={false}>
          <bufferGeometry>
            <bufferAttribute attach="attributes-position" args={[massCloud.positions, 3]} />
            <bufferAttribute attach="attributes-aCategory" args={[massCloud.categories, 1]} />
            <bufferAttribute attach="attributes-aJitter" args={[massCloud.jitter, 1]} />
          </bufferGeometry>
          <shaderMaterial
            ref={massMaterialRef}
            vertexShader={POINT_CLOUD_VERTEX}
            fragmentShader={POINT_CLOUD_FRAGMENT}
            uniforms={massUniforms}
            transparent
            depthWrite
          />
        </points>

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

      {/* Edema — a sparser, differently-coloured point cloud standing in
          for the diffuse peritumoral swelling: additive blending and a
          lower per-point opacity so it reads as a soft halo rather than a
          second solid cell population. */}
      <group ref={edemaGroup}>
        <points frustumCulled={false}>
          <bufferGeometry>
            <bufferAttribute attach="attributes-position" args={[edemaCloud.positions, 3]} />
            <bufferAttribute attach="attributes-aCategory" args={[edemaCloud.categories, 1]} />
            <bufferAttribute attach="attributes-aJitter" args={[edemaCloud.jitter, 1]} />
          </bufferGeometry>
          <shaderMaterial
            ref={edemaMaterialRef}
            vertexShader={POINT_CLOUD_VERTEX}
            fragmentShader={POINT_CLOUD_FRAGMENT}
            uniforms={edemaUniforms}
            transparent
            depthWrite={false}
            blending={THREE.AdditiveBlending}
          />
        </points>
      </group>

      {/* Margin — very sparse, faint point cloud standing in for
          surrounding tissue at the outer boundary. */}
      <group ref={marginGroup}>
        <points frustumCulled={false}>
          <bufferGeometry>
            <bufferAttribute attach="attributes-position" args={[marginCloud.positions, 3]} />
            <bufferAttribute attach="attributes-aCategory" args={[marginCloud.categories, 1]} />
            <bufferAttribute attach="attributes-aJitter" args={[marginCloud.jitter, 1]} />
          </bufferGeometry>
          <shaderMaterial
            ref={marginMaterialRef}
            vertexShader={POINT_CLOUD_VERTEX}
            fragmentShader={POINT_CLOUD_FRAGMENT}
            uniforms={marginUniforms}
            transparent
            depthWrite={false}
          />
        </points>
      </group>
    </group>
  );
}
