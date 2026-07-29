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
    float edge = smoothstep(0.0, 0.12, vUv.x) * smoothstep(1.0, 0.88, vUv.x)
      * smoothstep(0.0, 0.12, vUv.y) * smoothstep(1.0, 0.88, vUv.y);
    float alpha = (line * 0.7 + 0.06) * uOpacity * mix(0.35, 1.0, edge);
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
