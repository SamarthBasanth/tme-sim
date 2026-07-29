"use client";

import { useMemo, useRef } from "react";
import { useFrame } from "@react-three/fiber";
import * as THREE from "three";
import { useScrollProgressRef } from "@/components/scroll/ScrollProgressContext";
import { scanPlaneOpacity, scanPlaneY } from "@/lib/timeline";
import { damp } from "@/lib/scrollMath";

const VERTEX = /* glsl */ `
  varying vec2 vUv;
  void main() {
    vUv = uv;
    gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
  }
`;

const FRAGMENT = /* glsl */ `
  varying vec2 vUv;
  uniform float uOpacity;
  uniform vec3 uColor;

  void main() {
    vec2 grid = abs(fract(vUv * 14.0 - 0.5) - 0.5) / fwidth(vUv * 14.0);
    float line = 1.0 - min(min(grid.x, grid.y), 1.0);

    // Radial vignette instead of a square edge fade: a plane this size is
    // much larger than the tumor it sweeps through, so a fade that only
    // approaches (but never reaches) zero at the quad's literal edges left
    // a faint, uniformly-lit square visible against the dark background —
    // read as an unstyled debug quad floating on the model. Fading fully to
    // zero, and doing it radially rather than per-axis, removes both the
    // residual brightness and the square silhouette; what remains reads as
    // a circular scan disc, which also fits a round tumor better.
    vec2 centered = vUv - 0.5;
    float radial = length(centered) * 2.0;
    float radialFade = smoothstep(1.0, 0.5, radial);

    float alpha = line * 0.7 * uOpacity * radialFade;
    gl_FragColor = vec4(uColor, alpha);
  }
`;

/**
 * Thin scanning plane that sweeps through the model during "AI Analysis",
 * its Y position an explicit function of scroll progress (`scanPlaneY`) so
 * scrolling back up pulls the scan back with it, not forward again.
 */
export default function ScanPlane() {
  const progressRef = useScrollProgressRef();
  const meshRef = useRef<THREE.Mesh>(null);
  const materialRef = useRef<THREE.ShaderMaterial>(null);
  const ySmoothed = useRef(1.35);
  const opacitySmoothed = useRef(0);

  const uniforms = useMemo(
    () => ({
      uOpacity: { value: 0 },
      uColor: { value: new THREE.Color("#8fe0ff") },
    }),
    [],
  );

  useFrame((_, delta) => {
    const p = progressRef.current;
    ySmoothed.current = damp(ySmoothed.current, scanPlaneY(p), 5, delta);
    opacitySmoothed.current = damp(opacitySmoothed.current, scanPlaneOpacity(p), 4, delta);

    if (meshRef.current) meshRef.current.position.y = ySmoothed.current;
    if (materialRef.current) materialRef.current.uniforms.uOpacity.value = opacitySmoothed.current;
  });

  return (
    <mesh ref={meshRef} rotation={[-Math.PI / 2, 0, 0]}>
      <planeGeometry args={[2.6, 2.6, 1, 1]} />
      <shaderMaterial
        ref={materialRef}
        vertexShader={VERTEX}
        fragmentShader={FRAGMENT}
        uniforms={uniforms}
        transparent
        depthWrite={false}
        side={THREE.DoubleSide}
      />
    </mesh>
  );
}
