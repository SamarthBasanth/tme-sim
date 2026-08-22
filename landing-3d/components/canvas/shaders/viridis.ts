/**
 * The drug-concentration colour ramp, shared by every scene element that
 * paints a concentration.
 *
 * These are NOT invented brand colours. They are matplotlib's viridis control
 * points, copied from the real product's own ramp table
 * (`app/web/js/gl/ramps.js`, `STOPS.viridis`), which is what the Studio uses
 * to paint scalar overlays — so a reader who has seen a Vivum field overlay
 * sees the same ramp here rather than a third, decorative palette invented for
 * the landing page. Nine of the seventeen stops are enough at the sizes this
 * is drawn at; the ramp stays perceptually monotone either way, which is the
 * only property the scene depends on.
 *
 * Written as an unrolled `mix` chain rather than a `const vec3[]` lookup on
 * purpose: three.js compiles `ShaderMaterial` as GLSL ES 1.00 unless
 * `glslVersion` is set, and that version only permits constant-index
 * expressions into arrays.
 */
export const VIRIDIS_GLSL = /* glsl */ `
  vec3 viridis(float t) {
    t = clamp(t, 0.0, 1.0);
    // ramps.js STOPS.viridis, indices 0,2,4,...,16, /255.
    vec3 v0 = vec3(0.2667, 0.0039, 0.3294);
    vec3 v1 = vec3(0.2784, 0.1765, 0.4824);
    vec3 v2 = vec3(0.2314, 0.3216, 0.5451);
    vec3 v3 = vec3(0.1725, 0.4471, 0.5569);
    vec3 v4 = vec3(0.1294, 0.5686, 0.5490);
    vec3 v5 = vec3(0.1569, 0.6824, 0.5020);
    vec3 v6 = vec3(0.3686, 0.7882, 0.3843);
    vec3 v7 = vec3(0.6784, 0.8627, 0.1882);
    vec3 v8 = vec3(0.9922, 0.9059, 0.1451);

    float s = t * 8.0;
    vec3 c = v0;
    c = mix(c, v1, clamp(s - 0.0, 0.0, 1.0));
    c = mix(c, v2, clamp(s - 1.0, 0.0, 1.0));
    c = mix(c, v3, clamp(s - 2.0, 0.0, 1.0));
    c = mix(c, v4, clamp(s - 3.0, 0.0, 1.0));
    c = mix(c, v5, clamp(s - 4.0, 0.0, 1.0));
    c = mix(c, v6, clamp(s - 5.0, 0.0, 1.0));
    c = mix(c, v7, clamp(s - 6.0, 0.0, 1.0));
    c = mix(c, v8, clamp(s - 7.0, 0.0, 1.0));
    return c;
  }
`;

/**
 * Where drug enters the section.
 *
 * Six points, one per third-by-half block of the section, each the centroid of
 * the densest bin of a 26x16 histogram over the 4,295 measured ENDOTHELIAL
 * cell positions in `public/models/tumor.glb` (mesh "ENDOTHELIAL", POSITION
 * accessor). Coordinates are in the same world units as the rendered section
 * and are already expressed in the standing-up frame GLTFTumor puts the
 * planar asset in (the exporter lays the section in XZ with y == 0; a -90 deg
 * turn about X sends its z to y).
 *
 * So these are measured positions, not art direction: drug enters the tissue
 * where the vasculature actually is. What is illustrative is everything the
 * scene then does with them — the decay length, the front timing, and the
 * curves drawn between them.
 */
export const ENTRY_POINTS: ReadonlyArray<readonly [number, number]> = [
  [-0.865, -0.465],
  [-0.428, 0.061],
  [-0.249, -0.776],
  [-0.262, 0.468],
  [0.564, -0.333],
  [0.552, 0.072],
] as const;

/** Measured extent of the section in world units (GLB POSITION accessor bounds). */
export const SECTION_HALF_EXTENT = { x: 1.3, y: 0.8142 } as const;
