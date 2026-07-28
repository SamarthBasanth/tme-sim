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

uniform vec2  u_extent;      // um
uniform vec2  u_scale;       // world -> clip
uniform vec2  u_offset;
uniform float u_pointSize;   // device px
uniform vec3  u_typeColor[5];
uniform float u_typeVis[5];
uniform int   u_colorMode;   // 0 = cell type, 1 = field value at cell
uniform int   u_hasField;
uniform sampler2D u_field;
uniform sampler2D u_ramp;
uniform vec2  u_fieldRange;  // min, max in field units
uniform vec4  u_force;       // rgb + flag(w>0.5) : override colour

out vec4 v_color;

void main() {
  int t = int(a_type + 0.5);
  vec2 world = a_pos * u_extent;

  if (u_typeVis[t] < 0.5) {
    // Cull invisible channels in the vertex stage: zero-area point pushed
    // outside the frustum costs no rasterisation at all.
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

void main() {
  float f = texture(u_field, v_uv).r;
  float t = clamp((f - u_range.x) / max(1e-12, u_range.y - u_range.x), 0.0, 1.0);
  t = pow(t, u_gamma);
  vec3 c = texture(u_ramp, vec2(t, 0.5)).rgb;
  o_color = vec4(c * u_opacity, u_opacity);   // premultiplied
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
