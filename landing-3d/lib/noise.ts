/**
 * Deterministic 3D Perlin noise (Ken Perlin's "improved noise" reference
 * algorithm) plus fbm (fractal Brownian motion) and domain-warped fbm built
 * on top of it. CPU-side only — used once, inside `useMemo`, to displace and
 * color the procedural tumor geometry (see `tumorGeometry.ts`), never
 * per-frame, so raw performance is not a concern here.
 *
 * The permutation table is shuffled once with the project's existing seeded
 * PRNG (`mulberry32`, see prng.ts) rather than `Math.random()`, for the same
 * "pure function of its inputs" reason particles.tsx already cares about:
 * this module's exports are deterministic across reloads/renders.
 */
import { mulberry32 } from "./prng";

function buildPermutation(seed: number): Uint8Array {
  const p = new Uint8Array(256);
  for (let i = 0; i < 256; i++) p[i] = i;
  const rand = mulberry32(seed);
  for (let i = 255; i > 0; i--) {
    const j = Math.floor(rand() * (i + 1));
    const tmp = p[i];
    p[i] = p[j];
    p[j] = tmp;
  }
  const perm = new Uint8Array(512);
  for (let i = 0; i < 512; i++) perm[i] = p[i & 255];
  return perm;
}

const PERM = buildPermutation(7919);

function fade(t: number): number {
  return t * t * t * (t * (t * 6 - 15) + 10);
}

function lerp(a: number, b: number, t: number): number {
  return a + t * (b - a);
}

function grad(hash: number, x: number, y: number, z: number): number {
  const h = hash & 15;
  const u = h < 8 ? x : y;
  const v = h < 4 ? y : h === 12 || h === 14 ? x : z;
  return ((h & 1) === 0 ? u : -u) + ((h & 2) === 0 ? v : -v);
}

/** Classic Perlin noise, roughly in [-1, 1]. */
export function perlin3(x: number, y: number, z: number): number {
  const X = Math.floor(x) & 255;
  const Y = Math.floor(y) & 255;
  const Z = Math.floor(z) & 255;
  const xf = x - Math.floor(x);
  const yf = y - Math.floor(y);
  const zf = z - Math.floor(z);
  const u = fade(xf);
  const v = fade(yf);
  const w = fade(zf);

  const A = PERM[X] + Y;
  const AA = PERM[A] + Z;
  const AB = PERM[A + 1] + Z;
  const B = PERM[X + 1] + Y;
  const BA = PERM[B] + Z;
  const BB = PERM[B + 1] + Z;

  return lerp(
    lerp(
      lerp(grad(PERM[AA], xf, yf, zf), grad(PERM[BA], xf - 1, yf, zf), u),
      lerp(grad(PERM[AB], xf, yf - 1, zf), grad(PERM[BB], xf - 1, yf - 1, zf), u),
      v,
    ),
    lerp(
      lerp(grad(PERM[AA + 1], xf, yf, zf - 1), grad(PERM[BA + 1], xf - 1, yf, zf - 1), u),
      lerp(grad(PERM[AB + 1], xf, yf - 1, zf - 1), grad(PERM[BB + 1], xf - 1, yf - 1, zf - 1), u),
      v,
    ),
    w,
  );
}

/** Fractal Brownian motion: several octaves of Perlin noise summed and normalized to ~[-1, 1]. */
export function fbm3(
  x: number,
  y: number,
  z: number,
  octaves = 4,
  lacunarity = 2,
  gain = 0.5,
): number {
  let amp = 0.5;
  let freq = 1;
  let sum = 0;
  let norm = 0;
  for (let i = 0; i < octaves; i++) {
    sum += amp * perlin3(x * freq, y * freq, z * freq);
    norm += amp;
    amp *= gain;
    freq *= lacunarity;
  }
  return norm === 0 ? 0 : sum / norm;
}

/**
 * Domain-warped fbm: the sample point is pushed around by two more fbm
 * fields before the final fbm lookup. This is what turns "bumpy sphere"
 * into genuinely irregular, lobulated, non-radially-symmetric blobs —
 * plain fbm displacement along a normal still reads as a lumpy ball because
 * every octave is centered on the same origin; warping breaks that symmetry.
 */
export function warpedFbm3(
  x: number,
  y: number,
  z: number,
  warpAmount = 0.6,
  octaves = 4,
): number {
  const wx = fbm3(x + 5.2, y + 1.3, z + 7.1, 3);
  const wy = fbm3(x + 0.5, y + 3.9, z + 2.3, 3);
  const wz = fbm3(x + 8.1, y + 6.2, z + 0.7, 3);
  return fbm3(x + wx * warpAmount, y + wy * warpAmount, z + wz * warpAmount, octaves);
}
