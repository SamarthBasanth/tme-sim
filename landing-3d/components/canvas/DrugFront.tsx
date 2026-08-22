"use client";

import { useEffect, useMemo, useRef } from "react";
import { useFrame } from "@react-three/fiber";
import * as THREE from "three";
import { useScrollProgressRef } from "@/components/scroll/ScrollProgressContext";
import { drugFrontOpacity, drugFrontReveal } from "@/lib/timeline";
import { damp } from "@/lib/scrollMath";
import { ENTRY_POINTS, VIRIDIS_GLSL } from "./shaders/viridis";

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

  ${VIRIDIS_GLSL}

  void main() {
    if (vAlong > uProgress) discard;

    // Concentration falls with distance from the vessel it left, so the same
    // ramp that paints ConcentrationField paints this: bright at the wall,
    // dropping toward the front. A single flat-coloured line would have said
    // "a path was drawn"; this says "less arrives further out", which is the
    // whole reason the product is spatial.
    float c = 1.0 - 0.78 * vAlong;
    float head = smoothstep(uProgress - 0.06, uProgress, vAlong);
    vec3 col = viridis(c) + head * 0.2;
    float alpha = (0.82 - 0.32 * vAlong) * uFade;
    gl_FragColor = vec4(col, alpha);
  }
`;

/**
 * Drug leaving the vasculature and advancing into tissue.
 *
 * This is the old TrajectoryPath, repurposed rather than deleted. Its draw-on
 * material -- a tube whose fragments are discarded beyond `uProgress` along
 * its own length, so it reads as being traced in real time and un-traced when
 * you scroll back up -- is exactly the right effect for a front advancing, and
 * was wrong only for what it used to trace (a planned surgical/needle path
 * into the mass, which this product has no business depicting). The curves,
 * the colour and the name are new; the mechanism is the one that was already
 * here.
 *
 * Each stub starts at a measured endothelial density maximum (ENTRY_POINTS --
 * the same six the concentration field is sourced from) and reaches a short
 * way outward, all within the section's measured bounds. It sits just in front
 * of the section's own plane so it reads as movement through the tissue rather
 * than in the space around it. The extent it reaches is illustrative; the
 * points it starts from are not.
 */

/** [start (a measured entry point), bend, tip] per stub, in section world XY. */
const FRONTS: ReadonlyArray<ReadonlyArray<readonly [number, number]>> = [
  [ENTRY_POINTS[0], [-1.04, -0.56], [-1.19, -0.66]],
  [ENTRY_POINTS[1], [-0.52, 0.26], [-0.6, 0.5]],
  [ENTRY_POINTS[4], [0.78, -0.46], [1.0, -0.56]],
  [ENTRY_POINTS[5], [0.73, 0.27], [0.92, 0.46]],
];

/**
 * Filament-thin, and thinner than the vessel stand-ins they leave. At the
 * radius this inherited from the surgical path it drew before, four stubs
 * rendered as flat-shaded slabs -- opaque bars laid over the tissue, reading
 * as annotation rather than as anything moving through it.
 */
const FRONT_RADIUS = 0.005;

/** Just in front of the section plane (the asset itself is flat at z = 0). */
const FRONT_Z = 0.02;

export default function DrugFront() {
  const progressRef = useScrollProgressRef();
  const revealSmoothed = useRef(0);
  const fadeSmoothed = useRef(1);

  const geometries = useMemo(
    () =>
      FRONTS.map((stub) => {
        const curve = new THREE.CatmullRomCurve3(
          stub.map(([x, y]) => new THREE.Vector3(x, y, FRONT_Z)),
        );
        return new THREE.TubeGeometry(curve, 48, FRONT_RADIUS, 6, false);
      }),
    [],
  );

  // Built outside R3F's declarative tree, so these need an explicit dispose.
  useEffect(() => () => geometries.forEach((g) => g.dispose()), [geometries]);

  const uniforms = useMemo(
    () => ({ uProgress: { value: 0 }, uFade: { value: 1 } }),
    [],
  );

  /**
   * One ref per stub, and every one of them written each frame. R3F clones a
   * ShaderMaterial's uniform objects onto the material the first time it
   * applies the `uniforms` prop, so the four materials below hold four
   * independent uniform sets: writing through the literal above, or through a
   * single ref, would animate at most one stub and freeze the rest.
   */
  const materialRefs = useRef<Array<THREE.ShaderMaterial | null>>([]);
  const groupRef = useRef<THREE.Group>(null);

  useFrame((_, delta) => {
    const p = progressRef.current;
    revealSmoothed.current = damp(revealSmoothed.current, drugFrontReveal(p), 3.4, delta);
    fadeSmoothed.current = damp(fadeSmoothed.current, drugFrontOpacity(p), 3, delta);
    materialRefs.current.forEach((material) => {
      if (!material) return;
      material.uniforms.uProgress.value = revealSmoothed.current;
      material.uniforms.uFade.value = fadeSmoothed.current;
    });

    /*
     * Four transparent tubes, 48 segments x 6 radial each, drawn on every
     * frame of the entire scroll regardless of `uFade`. The fragments were
     * cheap individually but they are `transparent` and `depthWrite: false`,
     * so they were sorted and blended for six beats to contribute nothing for
     * five of them. `uFade` starts at 1 and damps DOWN toward the timeline, so
     * the gate reads the same smoothed value the shader does rather than
     * assuming anything about which beat is current.
     */
    if (groupRef.current) {
      groupRef.current.visible =
        fadeSmoothed.current > 0.004 && revealSmoothed.current > 0.002;
    }
  });

  return (
    <group ref={groupRef} visible={false}>
      {geometries.map((geometry, i) => (
        <mesh key={i} geometry={geometry}>
          <shaderMaterial
            ref={(el) => {
              materialRefs.current[i] = el;
            }}
            vertexShader={VERTEX}
            fragmentShader={FRAGMENT}
            uniforms={uniforms}
            transparent
            depthWrite={false}
          />
        </mesh>
      ))}
    </group>
  );
}
