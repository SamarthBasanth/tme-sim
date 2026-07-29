/**
 * Deterministic PRNG (mulberry32) used anywhere geometry needs "random"
 * scatter — e.g. the ambient particle field. A seeded generator is a pure
 * function of its seed, which matters for React's render purity rules
 * (`Math.random()` inside `useMemo` is flagged as impure since the same
 * inputs would not reliably reproduce the same output); this keeps the
 * particle field both lint-clean and stable across re-renders.
 */
export function mulberry32(seed: number): () => number {
  let a = seed;
  return function next() {
    a |= 0;
    a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
