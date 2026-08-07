// shaders.js — GLSL ES 3.00 sources for the cell viewer.
//
// Positions arrive from the TMSG payload as f32 micrometres. They are
// re-quantised on upload to normalised uint16 over the section extent, so the
// vertex shader receives a_pos in [0,1] and multiplies by u_extent. At the
// 8.1 mm extent of a Xenium section that is a 0.124 um step — roughly two
// orders of magnitude below the ~10 um cell pitch, i.e. invisible — and it
// halves the per-vertex fetch, which is the only thing that matters at
// 850k points.

export const POINTS_VS = `#version 300 es
precision highp float;

layout(location = 0) in vec2  a_pos;   // normalised [0,1] over the extent
layout(location = 1) in float a_type;  // CellType 0..4
// Colour-by-attribute channel. Currently carries the mechanism-attribution
// label (tmesim.sim.mechanism, 0..7); it is a generic per-cell u8 slot, so a
// future per-cell category rides it without another buffer. Disabled attributes
// read the constant set by gl.vertexAttrib1f, so this costs nothing when no
// attribute is loaded.
layout(location = 2) in float a_label;

uniform vec2  u_extent;      // um
uniform vec2  u_scale;       // world -> clip
uniform vec2  u_offset;
uniform float u_pointSize;   // device px
uniform vec3  u_typeColor[5];
uniform float u_typeVis[5];
uniform int   u_colorMode;   // 0 = cell type, 1 = field value at cell, 2 = label
uniform int   u_hasField;
uniform int   u_hasLabel;
uniform vec3  u_labelColor[8];
uniform float u_labelVis[8];
uniform sampler2D u_field;
uniform sampler2D u_ramp;
uniform vec2  u_fieldRange;  // min, max in field units
uniform vec4  u_force;       // rgb + flag(w>0.5) : override colour

out vec4 v_color;

void main() {
  int t = int(a_type + 0.5);
  int l = int(a_label + 0.5);
  vec2 world = a_pos * u_extent;

  bool labelMode = (u_colorMode == 2 && u_hasLabel == 1);
  if (u_typeVis[t] < 0.5 || (labelMode && u_labelVis[l] < 0.5)) {
    // Cull invisible channels in the vertex stage: zero-area point pushed
    // outside the frustum costs no rasterisation at all. A mechanism category
    // toggled off in the legend culls by exactly the same route.
    gl_Position  = vec4(2.0, 2.0, 2.0, 1.0);
    gl_PointSize = 0.0;
    v_color = vec4(0.0);
    return;
  }

  gl_Position  = vec4(world * u_scale + u_offset, 0.0, 1.0);
  gl_PointSize = u_pointSize;

  vec3 c;
  if (u_force.w > 0.5) {
    c = u_force.rgb;
  } else if (labelMode) {
    c = u_labelColor[l];
  } else if (u_colorMode == 1 && u_hasField == 1) {
    float f  = texture(u_field, a_pos).r;
    float t01 = clamp((f - u_fieldRange.x) / max(1e-12, u_fieldRange.y - u_fieldRange.x), 0.0, 1.0);
    c = texture(u_ramp, vec2(t01, 0.5)).rgb;
  } else {
    c = u_typeColor[t];
  }
  v_color = vec4(c, 1.0);
}`;

export const POINTS_FS = `#version 300 es
precision mediump float;

in  vec4 v_color;
out vec4 o_color;

uniform float u_round;  // >0.5 -> circular sprite (costs a discard)
uniform float u_alpha;
uniform float u_ring;   // >0.5 -> hollow ring (selection marker)

void main() {
  if (u_round > 0.5 || u_ring > 0.5) {
    vec2 d = gl_PointCoord - 0.5;
    float r2 = dot(d, d);
    if (r2 > 0.25) discard;
    if (u_ring > 0.5 && r2 < 0.135) discard;
  }
  o_color = vec4(v_color.rgb, u_alpha);
}`;

/* --------------------------------------------------------------- field --- */

export const FIELD_VS = `#version 300 es
precision highp float;

layout(location = 0) in vec2 a_uv;   // unit quad

uniform vec2 u_extent;
uniform vec2 u_scale;
uniform vec2 u_offset;

out vec2 v_uv;

void main() {
  v_uv = a_uv;
  vec2 world = a_uv * u_extent;
  gl_Position = vec4(world * u_scale + u_offset, 0.0, 1.0);
}`;

export const FIELD_FS = `#version 300 es
precision highp float;

in  vec2 v_uv;
out vec4 o_color;

uniform sampler2D u_field;
uniform sampler2D u_ramp;
uniform vec2  u_range;    // min, max in field units
uniform float u_opacity;
uniform float u_gamma;
uniform float u_alphaRamp; // >0.5 -> alpha tracks the value; 0 -> flat alpha

void main() {
  float f = texture(u_field, v_uv).r;
  float t = clamp((f - u_range.x) / max(1e-12, u_range.y - u_range.x), 0.0, 1.0);
  t = pow(t, u_gamma);
  vec3 c = texture(u_ramp, vec2(t, 0.5)).rgb;

  // ALPHA TRACKS THE VALUE, LINEARLY.
  //
  // With flat alpha this quad painted ramp(0) — for viridis a #440154 dark
  // purple — over the ENTIRE section extent, including every pixel outside the
  // tissue and every pixel where the concentration is exactly zero. On a real
  // Xenium section that is ~100% of the viewport at t=0 (the first time-lapse
  // frame is all zeros by construction: the drug has not entered yet), so the
  // cell layer was being read through a near-opaque wash and the tissue was
  // invisible. It also stated something false — a flat, textured colour where
  // the data is a constant zero reads as "a quantity was measured here".
  //
  // alpha = opacity * t is the honest encoding for a concentration: no drug,
  // no ink. It is LINEAR in the same t the hue uses, so it adds no contrast
  // stretch and invents no structure — a pixel at half the window is at half
  // the opacity, and the legend swatch is drawn with the same alpha ramp so
  // what is on screen is what the legend says.
  //
  // MISLEADING-CASE, and it is why the legend prints "0 = transparent": if the
  // operator narrows the window so its floor is ABOVE zero, everything below
  // the floor goes fully transparent rather than being clamped to the ramp's
  // dark end. Sub-window data then reads as absent. The overlay legend already
  // declares a narrowed window ("windowed · data lo → hi"); u_alphaRamp is the
  // escape hatch if a flat wash is ever wanted back.
  float a = u_opacity * mix(1.0, t, step(0.5, u_alphaRamp));
  o_color = vec4(c * a, a);   // premultiplied
}`;

/* -------------------------------------------------------------- outline -- */

export const LINE_VS = `#version 300 es
precision highp float;
layout(location = 0) in vec2 a_uv;
uniform vec2 u_extent;
uniform vec2 u_scale;
uniform vec2 u_offset;
void main() {
  vec2 world = a_uv * u_extent;
  gl_Position = vec4(world * u_scale + u_offset, 0.0, 1.0);
}`;

export const LINE_FS = `#version 300 es
precision mediump float;
uniform vec4 u_color;
out vec4 o_color;
void main() { o_color = u_color; }`;
