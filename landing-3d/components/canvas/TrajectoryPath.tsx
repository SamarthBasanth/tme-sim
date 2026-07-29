"use client";

import { useMemo, useRef } from "react";
import { useFrame } from "@react-three/fiber";
import * as THREE from "three";
import { useScrollProgressRef } from "@/components/scroll/ScrollProgressContext";
import { trajectoryOpacity, trajectoryReveal } from "@/lib/timeline";
import { damp } from "@/lib/scrollMath";

const VERTEX = /* glsl */ `
  varying float vAlong;
  void main() {
    vAlong = uv.x;
    gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
  }
`;

const FRAGMENT = /* glsl */ `
  varying float vAlong;
  uniform float uProgress;
  uniform float uFade;
  uniform vec3 uColor;

  void main() {
    if (vAlong > uProgress) discard;
    float head = smoothstep(uProgress - 0.05, uProgress, vAlong);
    vec3 col = uColor + head * vec3(0.5, 0.55, 0.6);
    float alpha = (0.85 - head * 0.15) * uFade;
    gl_FragColor = vec4(col, alpha);
  }
`;

/**
 * A proposed surgical/treatment trajectory drawn on as the "Treatment
 * Planning" section scrolls by — a tube whose material discards fragments
 * beyond `uProgress` along its length, so it reads as being traced in real
 * time rather than fading in as a whole. `uProgress` is a direct function of
 * scroll (via `trajectoryReveal`), so scrolling backward erases the path
 * from the same end it was drawn from.
 */
export default function TrajectoryPath() {
  const progressRef = useScrollProgressRef();
  const materialRef = useRef<THREE.ShaderMaterial>(null);
  const revealSmoothed = useRef(0);
  const fadeSmoothed = useRef(1);

  const geometry = useMemo(() => {
    const curve = new THREE.CatmullRomCurve3([
      new THREE.Vector3(1.9, 1.1, 0.9),
      new THREE.Vector3(1.1, 0.6, 0.55),
      new THREE.Vector3(0.55, 0.25, 0.3),
      new THREE.Vector3(0.1, 0.02, 0.05),
      new THREE.Vector3(-0.15, -0.05, -0.05),
    ]);
    return new THREE.TubeGeometry(curve, 64, 0.018, 8, false);
  }, []);

  const uniforms = useMemo(
    () => ({
      uProgress: { value: 0 },
      uFade: { value: 1 },
      uColor: { value: new THREE.Color("#7fb8c9") },
    }),
    [],
  );

  useFrame((_, delta) => {
    const p = progressRef.current;
    revealSmoothed.current = damp(revealSmoothed.current, trajectoryReveal(p), 3.4, delta);
    fadeSmoothed.current = damp(fadeSmoothed.current, trajectoryOpacity(p), 3, delta);
    if (materialRef.current) {
      materialRef.current.uniforms.uProgress.value = revealSmoothed.current;
      materialRef.current.uniforms.uFade.value = fadeSmoothed.current;
    }
  });

  return (
    <mesh geometry={geometry}>
      <shaderMaterial
        ref={materialRef}
        vertexShader={VERTEX}
        fragmentShader={FRAGMENT}
        uniforms={uniforms}
        transparent
      />
    </mesh>
  );
}
