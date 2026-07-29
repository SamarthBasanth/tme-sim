/**
 * Procedural, GBM-referenced tumor geometry + vertex-color builders.
 *
 * Grounded in how glioblastoma actually reads on contrast-enhanced MRI and
 * in gross/3D-printed surgical-planning specimens (BraTS challenge
 * documentation; see PMC11945730 "Advancing Precision: A Comprehensive
 * Review of MRI Segmentation Datasets from BraTS Challenges", and the
 * multi-material 3D-print surgical planning literature, e.g. PubMed
 * 34656248): a GBM mass is irregular and lobulated (not a smooth sphere),
 * and is famous for *ring enhancement* — a necrotic, non-enhancing core
 * (hypointense on T1Gd) surrounded by a brighter, irregular enhancing rim
 * (BraTS label convention: NCR=necrotic core, ET=enhancing tumor,
 * ED=peritumoral edema). Since we can't literally render the interior of an
 * opaque mesh, that anatomy is stylized here as one or two shallow,
 * noise-jittered necrotic "cavities" sunk into the lobulated surface, each
 * ringed by a brighter enhancing band right at its edge — the single most
 * recognizable GBM signature, made visible on the surface itself instead of
 * requiring a cutaway.
 *
 * All noise is evaluated once per vertex inside `useMemo` (never per
 * frame) using the deterministic CPU noise in `lib/noise.ts`.
 */
import * as THREE from "three";
import { mergeVertices } from "three/examples/jsm/utils/BufferGeometryUtils.js";
import { fbm3, warpedFbm3 } from "@/lib/noise";
import { mulberry32 } from "@/lib/prng";

function smoothstep(edge0: number, edge1: number, x: number): number {
  const t = THREE.MathUtils.clamp((x - edge0) / (edge1 - edge0), 0, 1);
  return t * t * (3 - 2 * t);
}

export interface NecroticPit {
  center: THREE.Vector3;
  angularRadius: number;
  seed: number;
}

/** Two off-center, differently sized cavities so the result reads as organic, not symmetric. */
export const DEFAULT_PITS: readonly NecroticPit[] = [
  { center: new THREE.Vector3(0.55, 0.4, 0.62).normalize(), angularRadius: 0.4, seed: 11.3 },
  { center: new THREE.Vector3(-0.48, -0.32, 0.36).normalize(), angularRadius: 0.24, seed: 47.1 },
];

function pitInfluence(dir: THREE.Vector3, pits: readonly NecroticPit[]): number {
  let best = 0;
  for (const pit of pits) {
    const d = THREE.MathUtils.clamp(dir.dot(pit.center), -1, 1);
    const angle = Math.acos(d);
    const jitter = fbm3(dir.x * 2.4 + pit.seed, dir.y * 2.4, dir.z * 2.4 + pit.seed, 3) * 0.16;
    const edge = pit.angularRadius + jitter;
    const infl = 1 - smoothstep(edge * 0.35, edge, angle);
    if (infl > best) best = infl;
  }
  return best;
}

const NECROTIC = new THREE.Color("#42221e");
const BASE_TISSUE = new THREE.Color("#8a2e2e");
const ENHANCING = new THREE.Color("#e2685c");

export interface MassOptions {
  lobeAmp?: number;
  fineAmp?: number;
  pitDepth?: number;
  pits?: readonly NecroticPit[];
}

/**
 * The tumor mass proper (necrotic core + enhancing rim, conceptually one
 * contiguous solid — real ring enhancement is the outer boundary of the
 * core, not a separate object). Returns a merged/indexed, smooth-normaled
 * geometry with a `color` attribute for `vertexColors`.
 */
export function createLobulatedMass(
  radius: number,
  detail: number,
  opts: MassOptions = {},
): THREE.BufferGeometry {
  const { lobeAmp = 0.14, fineAmp = 0.03, pitDepth = 0.09, pits = DEFAULT_PITS } = opts;

  let geometry: THREE.BufferGeometry = new THREE.IcosahedronGeometry(radius, detail);
  geometry = mergeVertices(geometry);

  const posAttr = geometry.attributes.position;
  const colorArr = new Float32Array(posAttr.count * 3);
  const dir = new THREE.Vector3();
  const tmpColor = new THREE.Color();

  for (let i = 0; i < posAttr.count; i++) {
    dir.set(posAttr.getX(i), posAttr.getY(i), posAttr.getZ(i)).normalize();

    const coarse = warpedFbm3(dir.x * 1.1, dir.y * 1.1, dir.z * 1.1, 0.7, 4);
    const fine = fbm3(dir.x * 5.5, dir.y * 5.5, dir.z * 5.5, 3);
    const pit = pitInfluence(dir, pits);

    const displacement = radius * (coarse * lobeAmp + fine * fineAmp - pit * pitDepth);
    const r = radius + displacement;
    posAttr.setXYZ(i, dir.x * r, dir.y * r, dir.z * r);

    // Base tissue, mottled (heterogeneous, not one flat color).
    const mottle = fbm3(dir.x * 3.2 + 91, dir.y * 3.2, dir.z * 3.2 + 91, 3);
    tmpColor.copy(BASE_TISSUE).lerp(ENHANCING, 0.16 + mottle * 0.12);

    // Enhancing rim: a bright band peaking right at the pit boundary.
    const rimBand = smoothstep(0.1, 0.45, pit) * smoothstep(0.95, 0.45, pit);
    tmpColor.lerp(ENHANCING, THREE.MathUtils.clamp(rimBand * 1.8, 0, 1));

    // Necrotic fill inside the pit — desaturated/darker, not black: this is
    // a shallow dimple read as tissue, not a crater or a rendering glitch.
    const necroticAmount = smoothstep(0.45, 0.85, pit);
    tmpColor.lerp(NECROTIC, necroticAmount * 0.72);

    // Cheap curvature-based AO: crevices (negative fine noise) read very
    // slightly darker — kept subtle so it reads as surface texture, not
    // self-shadowed cavities going to black.
    const ao = THREE.MathUtils.mapLinear(THREE.MathUtils.clamp(fine, -1, 1), -1, 1, 0.88, 1.04);
    tmpColor.multiplyScalar(ao);

    colorArr[i * 3] = tmpColor.r;
    colorArr[i * 3 + 1] = tmpColor.g;
    colorArr[i * 3 + 2] = tmpColor.b;
  }

  posAttr.needsUpdate = true;
  geometry.setAttribute("color", new THREE.BufferAttribute(colorArr, 3));
  geometry.computeVertexNormals();
  return geometry;
}

/**
 * Soft, low-frequency-only irregular blob for the edema halo — deliberately
 * lower amplitude / lower frequency than the mass so it reads as a diffuse
 * swelling, not a second hard lobulated sphere. Rendered by a separate
 * fresnel-driven shader (see ProceduralTumor) that fades it to nothing well
 * before its literal mesh boundary, so this geometry only needs to avoid a
 * perfectly spherical *silhouette* — the translucency does the rest.
 */
export function createEdemaGeometry(radius: number, detail: number): THREE.BufferGeometry {
  let geometry: THREE.BufferGeometry = new THREE.IcosahedronGeometry(radius, detail);
  geometry = mergeVertices(geometry);
  const posAttr = geometry.attributes.position;
  const dir = new THREE.Vector3();
  for (let i = 0; i < posAttr.count; i++) {
    dir.set(posAttr.getX(i), posAttr.getY(i), posAttr.getZ(i)).normalize();
    const n = warpedFbm3(dir.x * 0.75 + 40, dir.y * 0.75, dir.z * 0.75 + 40, 0.8, 3);
    const r = radius * (1 + n * 0.16);
    posAttr.setXYZ(i, dir.x * r, dir.y * r, dir.z * r);
  }
  posAttr.needsUpdate = true;
  geometry.computeVertexNormals();
  return geometry;
}

/**
 * Faint outer margin shell — very subtle irregular boundary indicator
 * (replaces a crystalline wireframe cage, which read as artificial/plastic).
 */
export function createMarginGeometry(radius: number, detail: number): THREE.BufferGeometry {
  let geometry: THREE.BufferGeometry = new THREE.IcosahedronGeometry(radius, detail);
  geometry = mergeVertices(geometry);
  const posAttr = geometry.attributes.position;
  const dir = new THREE.Vector3();
  for (let i = 0; i < posAttr.count; i++) {
    dir.set(posAttr.getX(i), posAttr.getY(i), posAttr.getZ(i)).normalize();
    const n = warpedFbm3(dir.x * 0.6 + 80, dir.y * 0.6, dir.z * 0.6 + 80, 0.7, 3);
    const r = radius * (1 + n * 0.08);
    posAttr.setXYZ(i, dir.x * r, dir.y * r, dir.z * r);
  }
  posAttr.needsUpdate = true;
  geometry.computeVertexNormals();
  return geometry;
}

export interface VesselCurveOptions {
  count?: number;
  seed?: number;
  radius?: number;
}

/**
 * A handful of thin, branching vessel-like curves hugging just outside the
 * mass surface — a cheap stand-in for surface vascularity (real GBM is
 * famously hypervascular). Built from a seeded PRNG so it is deterministic
 * across renders.
 */
export function createVesselCurves(opts: VesselCurveOptions = {}): THREE.CatmullRomCurve3[] {
  const { count = 5, seed = 2024, radius = 0.66 } = opts;
  const rand = mulberry32(seed);
  const curves: THREE.CatmullRomCurve3[] = [];

  for (let i = 0; i < count; i++) {
    const theta0 = rand() * Math.PI * 2;
    const phi0 = Math.acos(2 * rand() - 1);
    const start = new THREE.Vector3(
      Math.sin(phi0) * Math.cos(theta0),
      Math.sin(phi0) * Math.sin(theta0),
      Math.cos(phi0),
    );

    const points: THREE.Vector3[] = [];
    const segments = 4 + Math.floor(rand() * 3);
    let dir = start.clone();
    for (let s = 0; s <= segments; s++) {
      const wobble = new THREE.Vector3(
        (rand() - 0.5) * 0.22,
        (rand() - 0.5) * 0.22,
        (rand() - 0.5) * 0.22,
      );
      dir = dir.clone().add(wobble).normalize();
      // Deliberately well inside the nominal mass radius: the mass surface
      // itself is heavily noise-displaced (see createLobulatedMass) and can
      // dip well below `radius` locally, so a curve pinned just outside
      // `radius` would poke through empty space wherever the surface dips
      // — read as a stray floating line rather than a vessel under the
      // tissue. Sitting further inside keeps it plausibly "just under the
      // surface" everywhere instead.
      const r = radius * (0.62 + rand() * 0.12);
      points.push(dir.clone().multiplyScalar(r));
    }
    curves.push(new THREE.CatmullRomCurve3(points));
  }

  return curves;
}
