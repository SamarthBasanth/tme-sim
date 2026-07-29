/**
 * Cell-sprite point-cloud shader — the R3F analogue of the real product's
 * gl.POINTS cell viewer (see `app/web/js/gl/shaders.js`, `POINTS_VS` /
 * `POINTS_FS`): per-point category drives colour, and a circular sprite is
 * cut out of the square point quad in the fragment stage via a hard
 * `discard` at the unit circle — the same technique the real viewer uses,
 * and necessary here (not just an optimisation) because several of these
 * clouds render fully opaque with blending disabled, where alpha has no
 * effect on the framebuffer at all: a soft alpha-only edge would silently
 * paint the whole square quad, which is exactly what a hard discard avoids.
 * A softer alpha falloff right at the edge additionally gives clouds that
 * DO blend (edema, margin) a less jagged edge. Categories are picked with a
 * couple of `mix()`s rather than an indexed uniform array, since there is
 * no per-cell channel-visibility toggle to drive here.
 *
 * Point size follows exactly the formula three.js's own PointsMaterial uses
 * internally (`ShaderLib/points.glsl.js`): `size * (scale / -mvPosition.z)`
 * with `scale = canvasHeightPx * 0.5`, fed in as `uSizeScale` and kept in
 * sync with the canvas from the owning component (see ProceduralTumor.tsx).
 */
export const POINT_CLOUD_VERTEX = /* glsl */ `
  attribute float aCategory;
  attribute float aJitter;

  uniform vec3 uColor0;
  uniform vec3 uColor1;
  uniform vec3 uColor2;
  uniform float uSize;
  uniform float uSizeScale;
  uniform float uHighlight;
  uniform float uHighlightCategory;

  varying vec3 vColor;
  varying float vAlpha;

  void main() {
    vec3 base = mix(uColor0, uColor1, step(0.5, aCategory));
    base = mix(base, uColor2, step(1.5, aCategory));

    float emphasis = uHighlight * (1.0 - smoothstep(0.0, 0.6, abs(aCategory - uHighlightCategory)));
    vColor = base + emphasis * vec3(0.6, 0.52, 0.3);
    vAlpha = 0.72 + aJitter * 0.28;

    vec4 mvPosition = modelViewMatrix * vec4(position, 1.0);
    gl_Position = projectionMatrix * mvPosition;

    float size = uSize * (0.65 + aJitter * 0.65);
    gl_PointSize = size * (uSizeScale / -mvPosition.z);
  }
`;

export const POINT_CLOUD_FRAGMENT = /* glsl */ `
  uniform float uOpacity;
  varying vec3 vColor;
  varying float vAlpha;

  void main() {
    vec2 cxy = 2.0 * gl_PointCoord - 1.0;
    float r2 = dot(cxy, cxy);
    // A hard \`discard\` at the unit circle, not just a soft alpha falloff:
    // several of these point clouds render opaque (no blending, for the
    // dense core cloud) where alpha has no effect on the framebuffer at
    // all, so the circular shape has to come from \`discard\` itself, the
    // same approach the real product's viewer uses (see POINTS_FS in
    // app/web/js/gl/shaders.js). The softer edge alpha below only matters
    // for the clouds that DO blend (edema, margin).
    if (r2 > 1.0) discard;

    // A faint inner lift so each sprite reads as a small rounded cell body
    // rather than a flat, evenly-lit disc.
    float core = 1.0 - smoothstep(0.0, 0.6, r2);
    vec3 col = vColor + core * 0.12;

    float edgeAlpha = 1.0 - smoothstep(0.78, 1.0, r2);
    gl_FragColor = vec4(col, edgeAlpha * vAlpha * uOpacity);
  }
`;
