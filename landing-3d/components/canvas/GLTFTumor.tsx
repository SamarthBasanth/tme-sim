"use client";

import { useEffect, useMemo, useRef } from "react";
import { useFrame } from "@react-three/fiber";
import { useGLTF } from "@react-three/drei";
import * as THREE from "three";
import { useScrollProgressRef } from "@/components/scroll/ScrollProgressContext";
import { finalOffsetAmount, separationAmount } from "@/lib/timeline";
import { damp } from "@/lib/scrollMath";

/**
 * Real-asset path for a scan-derived `/models/tumor.glb`. Loaded with
 * `useGLTF` (suspends the nearest `<Suspense>`, i.e. the Loader in
 * SceneCanvas). Preloaded from TumorModel.tsx once the file is known to
 * exist so there is no pop-in when this component mounts.
 *
 * Segmentation is opt-in: if the authored GLB contains child objects named
 * "core", "edema", and "margin"/"tissue" (case-insensitive, as separate
 * meshes or groups), those are driven by the same separation timeline as
 * the procedural placeholder. Otherwise the whole model is treated as one
 * rigid body — segmentation simply becomes a no-op — so a first-pass export
 * with no naming convention still renders correctly, just without the
 * "regions pull apart" beat in Scene 3.
 */
export default function GLTFTumor({
  url,
  reducedMotion = false,
}: {
  url: string;
  reducedMotion?: boolean;
}) {
  const { scene } = useGLTF(url);
  const progressRef = useScrollProgressRef();

  const cloned = useMemo(() => scene.clone(true), [scene]);

  const rootGroup = useRef<THREE.Group>(null);
  const autoRotation = useRef(0);
  const scrollRotation = useRef(0);
  const sepSmoothed = useRef(0);
  const offsetSmoothed = useRef(0);

  const regions = useMemo(() => {
    const find = (name: string) =>
      cloned.getObjectByName(name) ??
      cloned.children.find((c) => c.name.toLowerCase().includes(name));
    return {
      core: find("core"),
      edema: find("edema"),
      margin: find("margin") ?? find("tissue"),
    };
  }, [cloned]);

  useEffect(() => {
    return () => {
      cloned.traverse((obj) => {
        if (obj instanceof THREE.Mesh) {
          obj.geometry.dispose();
          const mat = obj.material;
          if (Array.isArray(mat)) mat.forEach((m) => m.dispose());
          else mat?.dispose();
        }
      });
    };
  }, [cloned]);

  useFrame((_, delta) => {
    const p = progressRef.current;
    sepSmoothed.current = damp(sepSmoothed.current, separationAmount(p), 3.2, delta);
    offsetSmoothed.current = damp(offsetSmoothed.current, finalOffsetAmount(p), 3.2, delta);
    const sep = sepSmoothed.current;

    if (regions.core) regions.core.position.set(0, -0.4 * sep, 0);
    if (regions.edema) regions.edema.position.set(0.04 * sep, 0.45 * sep, 0.08 * sep);
    if (regions.margin) regions.margin.position.set(-0.06 * sep, 0.04 * sep, -0.45 * sep);

    if (rootGroup.current) {
      autoRotation.current += delta * (reducedMotion ? 0.035 : 0.12);
      scrollRotation.current = reducedMotion
        ? 0
        : damp(scrollRotation.current, p * 0.6, 2.4, delta);
      rootGroup.current.rotation.y = autoRotation.current + scrollRotation.current;
      rootGroup.current.position.x = reducedMotion ? 0 : offsetSmoothed.current * 1.15;
    }
  });

  return (
    <group ref={rootGroup}>
      <primitive object={cloned} />
    </group>
  );
}
