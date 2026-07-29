/**
 * Small, dependency-free math helpers used to derive every animated value
 * (camera pose, material opacity, separation offsets, label visibility...)
 * as a *pure function* of a single normalized scroll progress number.
 *
 * Keeping these pure (no hidden state) is what makes the whole experience
 * work identically forwards and backwards: scrolling up simply feeds a
 * smaller `progress` back in, and every curve below re-evaluates cleanly.
 */

export const clamp = (value: number, min = 0, max = 1): number =>
  Math.min(max, Math.max(min, value));

export const lerp = (a: number, b: number, t: number): number =>
  a + (b - a) * t;

/** Classic smoothstep easing, used to soften every linear segment below. */
export const smoothstep = (t: number): number => {
  const x = clamp(t, 0, 1);
  return x * x * (3 - 2 * x);
};

export const mapRange = (
  value: number,
  inMin: number,
  inMax: number,
  outMin: number,
  outMax: number,
  ease: boolean = true,
): number => {
  const span = inMax - inMin;
  const t = clamp((value - inMin) / (span === 0 ? 1e-6 : span), 0, 1);
  return lerp(outMin, outMax, ease ? smoothstep(t) : t);
};

/** A [progress, value] control point for {@link sampleCurve}. */
export type CurvePoint = readonly [t: number, value: number];

/**
 * Piecewise-linear (smoothstep-eased per segment) curve sampler.
 * This is the single primitive every authored animation curve in
 * `lib/timeline.ts` is built from: give it a handful of (progress, value)
 * keyframes and it interpolates smoothly between them, clamping outside.
 */
export const sampleCurve = (t: number, points: readonly CurvePoint[]): number => {
  if (points.length === 0) return 0;
  if (t <= points[0][0]) return points[0][1];
  const last = points[points.length - 1];
  if (t >= last[0]) return last[1];
  for (let i = 0; i < points.length - 1; i++) {
    const [t0, v0] = points[i];
    const [t1, v1] = points[i + 1];
    if (t >= t0 && t <= t1) {
      return mapRange(t, t0, t1, v0, v1);
    }
  }
  return last[1];
};

/** A curve that rises 0->1, plateaus at 1, then falls back to 0. */
export const trapezoid = (
  t: number,
  riseStart: number,
  riseEnd: number,
  fallStart: number,
  fallEnd: number,
): number =>
  sampleCurve(t, [
    [riseStart, 0],
    [riseEnd, 1],
    [fallStart, 1],
    [fallEnd, 0],
  ]);

/**
 * Frame-rate independent exponential damping (scalar). Used every render
 * frame to ease 3D objects toward their scroll-derived target so nothing
 * jitters or snaps, even if the scroll signal itself updates in bursts.
 */
export const damp = (
  current: number,
  target: number,
  lambda: number,
  dt: number,
): number => lerp(current, target, 1 - Math.exp(-lambda * dt));

export const damp3 = (
  current: readonly [number, number, number],
  target: readonly [number, number, number],
  lambda: number,
  dt: number,
): [number, number, number] => [
  damp(current[0], target[0], lambda, dt),
  damp(current[1], target[1], lambda, dt),
  damp(current[2], target[2], lambda, dt),
];
