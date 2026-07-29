/**
 * Fresnel + noise "diffuse cloud" shader for the peritumoral edema halo.
 *
 * Real peritumoral edema (T2/FLAIR hyperintensity) has no hard boundary —
 * it's a diffuse, irregular swelling that fades into surrounding tissue,
 * not a second solid, crisply-edged sphere. A flat-opacity transparent mesh
 * (the previous implementation) still reads as a hard shell because its
 * *silhouette* is always crisp even when its fill is faint. This shader
 * fixes that: opacity is driven mostly by a fresnel term (dense right at
 * the grazing-angle rim, nearly see-through where the surface faces the
 * camera) modulated by slow low-frequency noise "wisps", so the volume
 * genuinely dissolves at its edges instead of just fading its fill color.
 */
export const EDEMA_VERTEX = /* glsl */ `
  varying vec3 vNormal;
  varying vec3 vViewPosition;
  varying vec3 vWorldPosition;

  void main() {
    vNormal = normalize(normalMatrix * normal);
    vec4 mvPosition = modelViewMatrix * vec4(position, 1.0);
    vViewPosition = -mvPosition.xyz;
    vWorldPosition = (modelMatrix * vec4(position, 1.0)).xyz;
    gl_Position = projectionMatrix * mvPosition;
  }
`;

// Ashima Arts / Stefan Gustavson public-domain-style webgl simplex noise
// (the same widely-reused `snoise` implementation found throughout
// shader/creative-coding references) — evaluated per-fragment here since
// this is the one place in the scene that genuinely needs a soft,
// ever-so-slightly-animated volumetric edge.
const SIMPLEX_3D = /* glsl */ `
  vec3 mod289(vec3 x) { return x - floor(x * (1.0 / 289.0)) * 289.0; }
  vec4 mod289(vec4 x) { return x - floor(x * (1.0 / 289.0)) * 289.0; }
  vec4 permute(vec4 x) { return mod289(((x * 34.0) + 1.0) * x); }
  vec4 taylorInvSqrt(vec4 r) { return 1.79284291400159 - 0.85373472095314 * r; }

  float snoise(vec3 v) {
    const vec2 C = vec2(1.0 / 6.0, 1.0 / 3.0);
    const vec4 D = vec4(0.0, 0.5, 1.0, 2.0);

    vec3 i  = floor(v + dot(v, C.yyy));
    vec3 x0 = v - i + dot(i, C.xxx);

    vec3 g = step(x0.yzx, x0.xyz);
    vec3 l = 1.0 - g;
    vec3 i1 = min(g.xyz, l.zxy);
    vec3 i2 = max(g.xyz, l.zxy);

    vec3 x1 = x0 - i1 + C.xxx;
    vec3 x2 = x0 - i2 + C.yyy;
    vec3 x3 = x0 - D.yyy;

    i = mod289(i);
    vec4 p = permute(permute(permute(
      i.z + vec4(0.0, i1.z, i2.z, 1.0))
      + i.y + vec4(0.0, i1.y, i2.y, 1.0))
      + i.x + vec4(0.0, i1.x, i2.x, 1.0));

    float n_ = 0.142857142857;
    vec3 ns = n_ * D.wyz - D.xzx;

    vec4 j = p - 49.0 * floor(p * ns.z * ns.z);

    vec4 x_ = floor(j * ns.z);
    vec4 y_ = floor(j - 7.0 * x_);

    vec4 x = x_ * ns.x + ns.yyyy;
    vec4 y = y_ * ns.x + ns.yyyy;
    vec4 h = 1.0 - abs(x) - abs(y);

    vec4 b0 = vec4(x.xy, y.xy);
    vec4 b1 = vec4(x.zw, y.zw);

    vec4 s0 = floor(b0) * 2.0 + 1.0;
    vec4 s1 = floor(b1) * 2.0 + 1.0;
    vec4 sh = -step(h, vec4(0.0));

    vec4 a0 = b0.xzyw + s0.xzyw * sh.xxyy;
    vec4 a1 = b1.xzyw + s1.xzyw * sh.zzww;

    vec3 p0 = vec3(a0.xy, h.x);
    vec3 p1 = vec3(a0.zw, h.y);
    vec3 p2 = vec3(a1.xy, h.z);
    vec3 p3 = vec3(a1.zw, h.w);

    vec4 norm = taylorInvSqrt(vec4(dot(p0, p0), dot(p1, p1), dot(p2, p2), dot(p3, p3)));
    p0 *= norm.x;
    p1 *= norm.y;
    p2 *= norm.z;
    p3 *= norm.w;

    vec4 m = max(0.6 - vec4(dot(x0, x0), dot(x1, x1), dot(x2, x2), dot(x3, x3)), 0.0);
    m = m * m;
    return 42.0 * dot(m * m, vec4(dot(p0, x0), dot(p1, x1), dot(p2, x2), dot(p3, x3)));
  }
`;

export const EDEMA_FRAGMENT = /* glsl */ `
  varying vec3 vNormal;
  varying vec3 vViewPosition;
  varying vec3 vWorldPosition;

  uniform vec3 uColor;
  uniform float uOpacity;
  uniform float uTime;

  ${SIMPLEX_3D}

  float fbm(vec3 p) {
    float f = 0.0;
    f += 0.5 * snoise(p); p *= 2.02;
    f += 0.25 * snoise(p); p *= 2.03;
    f += 0.125 * snoise(p);
    return f;
  }

  void main() {
    vec3 viewDir = normalize(vViewPosition);
    // \`abs()\`, not \`max(...,0.0)\`: this material renders both faces
    // (DoubleSide, so the far side of the volume shows through the near
    // side, for depth). Clamping negative dot-products to 0 made every
    // back-facing triangle (the entire far hemisphere) evaluate to the
    // *maximum* fresnel term regardless of its actual grazing angle, which
    // filled in as a near-opaque ball instead of a thin halo. \`abs()\` gives
    // both faces the same "thin at direct-facing, brighter at grazing"
    // falloff.
    float ndotv = dot(normalize(vNormal), viewDir);
    float fresnel = pow(1.0 - abs(ndotv), 2.4);

    float wisp = fbm(vWorldPosition * 1.5 + vec3(0.0, uTime * 0.025, 0.0));
    float density = smoothstep(-0.35, 0.55, wisp);

    float alpha = uOpacity * (fresnel * 0.9 + 0.02) * mix(0.3, 1.0, density);
    alpha = clamp(alpha, 0.0, 0.34);

    gl_FragColor = vec4(uColor, alpha);
  }
`;
