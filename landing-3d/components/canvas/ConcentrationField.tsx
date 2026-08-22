"use client";

import { useMemo, useRef } from "react";
import { useFrame } from "@react-three/fiber";
import * as THREE from "three";
import { useScrollProgressRef } from "@/components/scroll/ScrollProgressContext";
import { concentrationFieldOpacity, frontAdvance } from "@/lib/timeline";
import { damp } from "@/lib/scrollMath";
import { ENTRY_POINTS, SECTION_HALF_EXTENT, VIRIDIS_GLSL } from "./shaders/viridis";

const VERTEX = /* glsl */ `
  varying vec2 vUv;
  void main() {
    vUv = uv;
    gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
  }
`;

const entry = (i: number) => `vec2(${ENTRY_POINTS[i][0].toFixed(3)}, ${ENTRY_POINTS[i][1].toFixed(3)})`;

const FRAGMENT = /* glsl */ `
  varying vec2 vUv;
  uniform float uOpacity;
  uniform float uAdvance;
  uniform vec2 uHalfExtent;

  ${VIRIDIS_GLSL}

  // Contribution of one vessel entry point at distance d, decaying over the
  // current penetration length. Not a solved field -- an exponential falloff
  // is the shape a diffusion-limited source has, which is the point being
  // made, but the number that matters (how far drug actually gets) is what
  // the solver returns, not what this draws.
  float fromVessel(vec2 p, vec2 s, float lambda) {
    return exp(-length(p - s) / lambda);
  }

  void main() {
    // vUv -> the section's own world coordinates, so the entry points below
    // are literal measured positions rather than shader-space guesses.
    vec2 p = (vUv * 2.0 - 1.0) * uHalfExtent;

    // The front advances by the drug reaching further, not by a plane moving
    // across the tissue: penetration length grows with scroll.
    float lambda = mix(0.045, 0.34, uAdvance);

    float sum = 0.0;
    sum += fromVessel(p, ${entry(0)}, lambda);
    sum += fromVessel(p, ${entry(1)}, lambda);
    sum += fromVessel(p, ${entry(2)}, lambda);
    sum += fromVessel(p, ${entry(3)}, lambda);
    sum += fromVessel(p, ${entry(4)}, lambda);
    sum += fromVessel(p, ${entry(5)}, lambda);
    float c = 1.0 - exp(-sum);

    // Iso-concentration contours, drawn the way a field is drawn and not the
    // way a scanner is: evenly spaced levels of the scalar itself, made
    // screen-space thin with fwidth so they stay one line wide at any zoom.
    float levels = c * 6.0;
    float d = abs(fract(levels - 0.5) - 0.5) / max(fwidth(levels), 1e-5);
    float line = (1.0 - min(d, 1.0)) * smoothstep(0.03, 0.14, c);

    // Fade out at the quad border so the field never ends on a hard rectangle
    // edge (the section's own outline is ragged; a crisp box reads as a debug
    // quad, which is exactly what the grid this replaced used to read as).
    vec2 e = min(vUv, 1.0 - vUv) * 2.0;
    float edge = smoothstep(0.0, 0.22, min(e.x, e.y));

    float alpha = (c * 0.42 + line * 0.75) * uOpacity * edge;
    gl_FragColor = vec4(viridis(c), alpha);
  }
`;

/**
 * The drug reaching the tissue, painted on the section it is reaching.
 *
 * This replaces the sweeping "AI Analysis" scan plane. The sweep is kept --
 * it was good motion -- but what sweeps is now the thing the product actually
 * computes: concentration spreading out of the vasculature. Six exponential
 * sources sit on measured endothelial density maxima (see ENTRY_POINTS), the
 * penetration length grows with scroll progress (`frontAdvance`), and the
 * result is drawn as viridis-ramped iso-concentration contours in the
 * product's own ramp.
 *
 * The plane is sized to the section's own measured extent and laid in the same
 * XY plane GLTFTumor stands the planar asset up into, so the field reads as
 * painted ON the tissue rather than as a slab passing through it. That matters
 * here: the asset is one measured section with nothing along the third axis,
 * so a plane travelling through depth would be animating an axis the data does
 * not have.
 *
 * As always in this scene, the visual state is an explicit function of scroll
 * progress, so scrolling back up retracts the front instead of replaying it.
 */
/**
 * Below this smoothed opacity the quad contributes nothing a display can show
 * (alpha peaks at `c * 0.42 + line * 0.75` times this, on a #07080a ground) so
 * it is skipped entirely rather than drawn transparent. See the note on
 * `meshRef` below for why that is worth doing here in particular.
 */
const VISIBLE_EPSILON = 0.004;

export default function ConcentrationField() {
  const progressRef = useScrollProgressRef();
  const meshRef = useRef<THREE.Mesh>(null);
  const materialRef = useRef<THREE.ShaderMaterial>(null);
  const advanceSmoothed = useRef(0);
  const opacitySmoothed = useRef(0);

  const uniforms = useMemo(
    () => ({
      uOpacity: { value: 0 },
      uAdvance: { value: 0 },
      uHalfExtent: {
        value: new THREE.Vector2(SECTION_HALF_EXTENT.x, SECTION_HALF_EXTENT.y),
      },
    }),
    [],
  );

  useFrame((_, delta) => {
    const p = progressRef.current;
    advanceSmoothed.current = damp(advanceSmoothed.current, frontAdvance(p), 5, delta);
    opacitySmoothed.current = damp(opacitySmoothed.current, concentrationFieldOpacity(p), 4, delta);

    if (materialRef.current) {
      materialRef.current.uniforms.uAdvance.value = advanceSmoothed.current;
      materialRef.current.uniforms.uOpacity.value = opacitySmoothed.current;
    }

    /*
     * THE MOST EXPENSIVE SINGLE SURFACE IN THE SCENE, AND IT WAS DRAWN FOR THE
     * WHOLE SCROLL.
     *
     * This quad spans the section's full measured extent and its fragment
     * shader is not cheap: six `exp()` source terms, a viridis lookup and an
     * `fwidth()` derivative per pixel, all `transparent` and therefore
     * unculled by the depth buffer. `uOpacity` is authored to be zero across
     * most of the timeline — the field belongs to one beat — but a zero
     * uniform does not save a single fragment. Every one of those pixels was
     * still shaded, blended and discarded on every frame of every beat.
     *
     * `visible` is the flag that actually removes the draw call. FieldContours
     * already did this for its lines; this surface is far larger and had no
     * equivalent.
     */
    if (meshRef.current) {
      meshRef.current.visible = opacitySmoothed.current > VISIBLE_EPSILON;
    }
  });

  return (
    <mesh ref={meshRef} visible={false} position={[0, 0, -0.03]} renderOrder={-1}>
      <planeGeometry
        args={[SECTION_HALF_EXTENT.x * 2, SECTION_HALF_EXTENT.y * 2, 1, 1]}
      />
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
