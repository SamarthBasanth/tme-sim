"use client";

/**
 * Soft, restrained "medical imaging suite" lighting: a cool key light, a
 * warm-neutral fill, and a dim ambient floor. Real-time shadow maps are
 * deliberately skipped everywhere in this scene (per the spec's "minimal
 * real-time shadows" guidance) — depth instead comes from material
 * transmission/emissive treatment and the rim light below, which is far
 * cheaper on mobile GPUs.
 */
export default function SceneLighting() {
  return (
    <>
      <ambientLight intensity={0.35} />
      <directionalLight position={[3.5, 4, 2.5]} intensity={1.1} color="#f4f7fa" />
      <directionalLight position={[-3, -1.5, -2]} intensity={0.35} color="#6fa8c9" />
      <pointLight position={[0, 1.4, 2.2]} intensity={0.5} color="#bfe3ee" distance={8} decay={2} />
    </>
  );
}
