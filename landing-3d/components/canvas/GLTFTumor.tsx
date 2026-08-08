"use client";

import { useEffect, useMemo, useRef } from "react";
import { useFrame } from "@react-three/fiber";
import { useGLTF } from "@react-three/drei";
import * as THREE from "three";
import { useScrollProgressRef } from "@/components/scroll/ScrollProgressContext";
import { sectionOffsetAmount, separationAmount } from "@/lib/timeline";
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
/** World-unit point size for point-cloud assets; see the fixup effect below. */
const POINT_SIZE_UNITS = 0.009;

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

  /**
   * Point-cloud assets need their material fixed up on load. The glTF spec has
   * no point-size property, so GLTFLoader always builds `PointsMaterial` with
   * `size = 1, sizeAttenuation = false` -- i.e. 1 screen pixel per point,
   * which renders 40k real cells as a barely-visible dust on a dark ground.
   *
   * POINT_SIZE_UNITS is in WORLD units with attenuation on, so it is a real
   * physical size rather than a screen artefact: the exporter's scale is
   * ~3111 um per scene unit, and the reconstruction grid is dx = 12 um, so one
   * grid cell is ~0.0039 units. The value below is a few grid cells across,
   * which closes the gaps between neighbouring cells without smearing the
   * tissue structure (the section outline, the tears, and the vessel tracks
   * traced by endothelial cells all have to stay legible -- that visible
   * structure is the entire argument for using a real section here).
   */
  /**
   * Is this a PLANAR asset? The real Xenium export is one measured tissue
   * section, so every point has y == 0 exactly -- there is no measured depth
   * and none was invented. Laid in the XZ plane (the exporter's convention) it
   * is edge-on to the camera and renders as a flat sliver, which reads as a
   * broken model rather than as what it is.
   *
   * A tissue section IS a flat slice, so the honest presentation is to stand
   * it up facing the viewer, like a slide, rather than to fake volume. Detect
   * planarity from the data instead of hardcoding it, so a genuinely
   * volumetric asset dropped in later is left alone.
   */
  const isPlanar = useMemo(() => {
    const box = new THREE.Box3().setFromObject(cloned);
    const size = box.getSize(new THREE.Vector3());
    const largest = Math.max(size.x, size.y, size.z);
    return largest > 0 && Math.min(size.x, size.y, size.z) / largest < 0.02;
  }, [cloned]);

  useEffect(() => {
    cloned.traverse((obj) => {
      if (obj instanceof THREE.Points) {
        const mat = obj.material as THREE.PointsMaterial;
        mat.size = POINT_SIZE_UNITS;
        mat.sizeAttenuation = true;
        mat.transparent = true;
        /*
         * ADDITIVE, and this is the difference between the section reading as
         * tissue and reading as noise.
         *
         * Under normal blending each point simply paints over whatever is
         * behind it, so 40,887 overlapping cells average out to one flat
         * wash -- and since 58.2% of them are tumour red, that wash was a
         * muddy maroon in which none of the structure was visible. Additive
         * accumulates instead: where cells are dense the luminance sums and
         * the region brightens, where they are sparse it stays near the
         * ground. Density becomes something you can see.
         *
         * That is an honest encoding, not a cosmetic glow. Local density is
         * measured data -- it is the same quantity the reconstruction meshes
         * onto its grid -- so brightness here is carrying real information
         * rather than asserting any. It is also how the fluorescence imaging
         * this section came from actually looks, which is the register the
         * page is claiming to be in.
         *
         * The HUES are untouched. They remain exactly the CELLTYPE_COLORS
         * that `tmesim.viz.explainer` uses for every figure in results/, so
         * the page and the papers still agree and the SpecimenCaption swatches
         * remain the key. Only accumulation changed.
         */
        mat.blending = THREE.AdditiveBlending;
        /*
         * Per-point alpha has to drop a long way under additive or the dense
         * interior clips to white and destroys exactly the structure this was
         * meant to reveal -- the vessel tracks and the tears go first,
         * because they read as EDGES against dense tissue. At 0.55 the
         * densest region sits just below clipping, so the full range from
         * section outline to core stays inside the ramp.
         */
        mat.opacity = 0.55;
        // Real cells overlap in projection; without this the far side of the
        // section punches holes through the near side as points sort badly.
        // Additive is order-independent, so this is now doubly correct.
        mat.depthWrite = false;
        mat.needsUpdate = true;
      }
    });
  }, [cloned]);

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
        // Points as well as Mesh: the real Xenium asset is a point cloud, and
        // matching only Mesh here leaked its 40k-vertex geometry on unmount.
        if (obj instanceof THREE.Mesh || obj instanceof THREE.Points) {
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
    offsetSmoothed.current = damp(offsetSmoothed.current, sectionOffsetAmount(p), 3.2, delta);
    const sep = sepSmoothed.current;

    if (regions.core) regions.core.position.set(0, -0.4 * sep, 0);
    if (regions.edema) regions.edema.position.set(0.04 * sep, 0.45 * sep, 0.08 * sep);
    if (regions.margin) regions.margin.position.set(-0.06 * sep, 0.04 * sep, -0.45 * sep);

    if (rootGroup.current) {
      if (isPlanar) {
        // A flat section must never turn edge-on, or it vanishes. Instead of
        // a full spin, OSCILLATE within a shallow angle so the slice keeps
        // facing the viewer and the parallax still reads as depth.
        autoRotation.current += delta * (reducedMotion ? 0.06 : 0.22);
        const sway = Math.sin(autoRotation.current) * (reducedMotion ? 0.05 : 0.16);
        scrollRotation.current = reducedMotion
          ? 0
          : damp(scrollRotation.current, p * 0.35, 2.4, delta);
        rootGroup.current.rotation.y = sway + scrollRotation.current;
        rootGroup.current.rotation.x = Math.cos(autoRotation.current * 0.6) * 0.05;
      } else {
        autoRotation.current += delta * (reducedMotion ? 0.035 : 0.12);
        scrollRotation.current = reducedMotion
          ? 0
          : damp(scrollRotation.current, p * 0.6, 2.4, delta);
        rootGroup.current.rotation.y = autoRotation.current + scrollRotation.current;
      }
      rootGroup.current.position.x = reducedMotion ? 0 : offsetSmoothed.current * 1.15;
    }
  });

  return (
    <group ref={rootGroup}>
      {/* Stand a planar section up to face the camera. The exporter lays the
          measured section in the XZ plane (y == 0); -90 deg about X turns that
          "floor" into a facing slide without touching a single coordinate, so
          every distance on screen is still a true measured distance. */}
      <group rotation={isPlanar ? [-Math.PI / 2, 0, 0] : [0, 0, 0]}>
        <primitive object={cloned} />
      </group>
    </group>
  );
}
