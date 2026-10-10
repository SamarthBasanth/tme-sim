"use client";

import { ENTRY_POINTS } from "./shaders/viridis";

/**
 * The rig for a specimen on a stage, and an honest account of how little of
 * this scene it reaches.
 *
 * WHAT IS ACTUALLY LIT, COUNTED RATHER THAN ASSUMED. On the path that ships —
 * `public/models/tumor.glb` present, so `GLTFTumor` mounts — the scene
 * contains exactly ONE material that responds to a light: `VesselNetwork`'s
 * `meshStandardMaterial` in Markers.tsx, two `TubeGeometry` tubes of radius
 * 0.007 world units at up to 0.55 opacity, revealed late in the scroll.
 * Everything else is unlit by construction and cannot be changed by anything
 * in this file:
 *
 *     the tissue        THREE.Points / PointsMaterial, additively blended
 *     the particles     PointsMaterial
 *     the field, fronts ShaderMaterial
 *     the contours      Line2 / LineMaterial
 *     the source marker meshBasicMaterial
 *     the annotations   meshBasicMaterial + drei <Html>
 *
 * So the previous rig here — a 0.35 ambient, a 1.1 key, a 0.35 fill and a 0.5
 * point — was four lights shading a pair of 5-pixel translucent noodles. That
 * is worth writing down because it sets the ceiling on what tuning lights can
 * do for this page, and it is the reason the colour work that actually changed
 * the image is in SceneCanvas's `onCreated` (tone mapping and output colour
 * space) rather than here. Do not expect a lighting change to show up in a
 * screenshot of the hero; it cannot.
 *
 * THE SECOND CONSUMER, WHICH IS WHY THIS IS NOT DEAD CODE. `ProceduralTumor`
 * is a real mesh with a `meshStandardMaterial`, and it is the documented
 * drop-out: remove the GLB and the page still renders a tumour. On that path
 * every light below does visible work, and a flat ambient was the thing
 * flattening it.
 *
 * WHAT CHANGED, AND WHY EACH ONE IS MOTIVATED RATHER THAN CONVENTIONAL.
 *
 *   AMBIENT -> HEMISPHERE. `ambientLight` adds the same value to every surface
 *   at every orientation; it is the one light that can only ever remove form.
 *   A hemisphere costs the same (one uniform, no extra pass) and carries a
 *   direction: cool from above, near-ground below. The ground colour is the
 *   scene background rather than black, so an unlit underside falls toward the
 *   page's own ink instead of toward a hole.
 *
 *   KEY, WARM, AND OFF-AXIS. The old key sat at #f4f7fa — a cool white — and
 *   the fill at #6fa8c9, also cool. Two cool sources is not a key and a fill,
 *   it is one light arriving twice, and it is the specific thing that reads as
 *   a default rig. This key is a low-saturation WARM white against a clearly
 *   cool fill, so there is a real temperature split across the form. It is
 *   warm-neutral and not amber on purpose: amber was removed from this scene
 *   site-wide when the reticle went (see Markers.tsx), and no new hue is
 *   being reintroduced here.
 *
 *   FILL, COOLER AND WEAKER, AND FROM BEHIND. Moved behind the section plane
 *   as well as to the opposite side, so it doubles as the rim that separates a
 *   translucent tube from a near-black ground. That separation is the cue a
 *   shadow would normally provide, and this scene correctly has no shadow maps
 *   (SceneCanvas explains why: there is nothing here that can cast one). Key
 *   to fill is 4.3 : 1 before the hemisphere, which is a lit subject rather
 *   than an evenly-exposed one.
 *
 *   PRACTICAL, WITH FALLOFF, PLACED ON THE MEASUREMENT. The old point light
 *   floated at [0, 1.4, 2.2] — in front of and above the section, pointing at
 *   nothing in particular, with `distance: 8` on a subject 2.6 units across,
 *   i.e. falloff so slow it acted as a second ambient. This one sits just off
 *   the plane at `ENTRY_POINTS[0]`, the measured endothelial density maximum
 *   that `ConcentrationSource` also marks and that the vessel stand-ins are
 *   drawn through, with `distance` cut to the scale of the thing it is lighting.
 *   It is the only light in the scene whose position is a measurement rather
 *   than a composition, which is the register this page is in.
 *
 * COST. Four lights before, four after; no shadow maps before or after; no
 * environment map, no light probe, no second pass. The hemisphere replaces the
 * ambient one-for-one in the shader's uniform block. This is a re-aim, not an
 * addition, and it was measured not to move frame time — see the note with
 * this change.
 */

/** The measured endothelial density maximum the vessel stand-ins run through. */
const SOURCE = ENTRY_POINTS[0];

export default function SceneLighting() {
  return (
    <>
      {/* Sky cool, ground the page's own background rather than black. */}
      <hemisphereLight args={["#b9d6e6", "#0b1116", 0.5]} />

      {/* Key: warm-neutral, high and front-right of the slide. */}
      <directionalLight position={[2.9, 3.6, 2.4]} intensity={1.3} color="#fff2e4" />

      {/* Fill/rim: distinctly cooler, 4.3x weaker, and BEHIND the plane. */}
      <directionalLight position={[-3.2, -1.1, -2.6]} intensity={0.3} color="#5f9ed0" />

      {/*
        Practical, on the measured source, with falloff at the scale of the
        section (half-extent 1.3 x 0.81) instead of the old 8.
      */}
      <pointLight
        position={[SOURCE[0], SOURCE[1], 0.5]}
        intensity={0.85}
        color="#cfeaf4"
        distance={2.4}
        decay={2}
      />
    </>
  );
}
